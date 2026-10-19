"""How much of a reported LINE-1 span is supported by read alignment?

Brandler, 2026-10-09: "we expect some line-1 insertions to be full length but this
should be based on read alignment and not imputed."

Background. In `_sample_insertion_span_and_orientation` (mei_support.py) there are
branches that set the far end of the consensus footprint to the subfamily's consensus
length when one side has a strong anchor and the other has a strong poly-A run:

    elif l_strong and r_poly_strong:
        tlen = max(l_target_len, r_target_len)
        raw_end = max(tlen, max(l_start, l_end))

A span produced that way ends exactly at the consensus length without a read ever
reaching it. A genuinely full-length insertion ends there too, so "ends at the
consensus length" is necessary but NOT sufficient evidence of imputation. This script
measures how large the ambiguous set is; it cannot resolve it from VCF fields alone.

Consensus lengths come from the panels, not from the data. This script used to take
each subfamily's consensus length to be the largest MEI_3P observed for it, because
`retrotransposon_db/dfam/dfam_human_mei_l1_alu_sva.fasta` was not on disk. That proxy
made the headline statistic true by construction: no call can end past an axis defined
as the furthest call, so "ends at the consensus length" could not distinguish a
genuine full-length insertion from an imputed one, and an overrun could not be seen
at all. The panels are now built and hash-pinned
(`scripts/fetch_dfam_mei_panel.py`, see `retrotransposon_db/dfam/*.manifest.json`),
and lengths are read from them.

Two axes exist and they are not interchangeable. The caller's MEI_5P/MEI_3P
coordinates are on the full-length axis (`L1HS_5end` calls reach MEI_3P 6018, the
`L1HS_full` length, not their own 2136 bp fragment), so a label is resolved to its own
`{subfamily}_full` record when the panel has one, else to its family canonical
(`L1HS_full`, the same fallback the annotate projection uses), else -- only for a
label with no full record -- to its Dfam fragment. `consensus_source` records which,
per call. Labels the panels do not know keep the old data-derived estimate and are
counted and named in the output.

Usage:
    python line1_span_support.py <dir-of-vcfs> [--out-dir DIR] [--min-full-span 5500]
                                 [--public-data-dir DIR] [--mei-panel FASTA]...
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import sys
from pathlib import Path

import pandas as pd

_PACKAGE_ROOTS = tuple(
    candidate
    for parent in Path(__file__).resolve().parents
    for candidate in (parent / "src", parent / "retrotransposon-miner" / "src")
)
for _root in _PACKAGE_ROOTS:
    if (_root / "retro_miner" / "__init__.py").exists():
        if str(_root) not in sys.path:
            sys.path.insert(0, str(_root))
        break

from retro_miner.mei_panel_index import (  # noqa: E402
    find_panel_paths,
    resolve_consensus_lengths,
)

SUPPORT_KEYS = ("SR_L", "SR_R", "DPE_L", "DPE_R", "MEI_MAPPED", "polyA_MAPPED")
NUMERIC_INFO = ("MEI_SPAN", "MEI_5P", "MEI_3P", "POLYA_MIN_BP")

PROXY_SOURCE = "proxy_max_MEI_3P"


def load_cohort(vcf_dir: str | Path) -> pd.DataFrame:
    """One row per call across every *.vcf in the directory, with SUPPORT parsed out."""
    from retro_miner.gene_annotation import read_vcf

    paths = sorted(glob.glob(os.path.join(str(vcf_dir), "*.vcf")))
    if not paths:
        raise FileNotFoundError(f"no *.vcf in {vcf_dir}")
    frames = []
    for path in paths:
        _, records = read_vcf(path)
        frame = pd.DataFrame([{"chrom": r.chrom, "pos": r.pos, **r.info} for r in records])
        frame["sample"] = os.path.basename(path).split("_")[0]
        frames.append(frame)
    calls = pd.concat(frames, ignore_index=True)
    for col in NUMERIC_INFO:
        calls[col] = pd.to_numeric(calls.get(col), errors="coerce")
    support = calls.get("SUPPORT", pd.Series("", index=calls.index)).fillna("")
    for key in SUPPORT_KEYS:
        calls[key] = pd.to_numeric(
            support.str.extract(rf"(?:^|\|){re.escape(key)}:(\d+)")[0], errors="coerce").fillna(0)
    calls["polyA_side"] = support.str.extract(r"polyA_side:(\w)")[0]
    calls["subfamily"] = calls["MEISUBFAMILY"].fillna("").str.split("#").str[0]
    calls["sentinel"] = (calls["MEI_SPAN"] <= 0) | (calls["MEI_5P"] < 0)
    return calls


def estimate_consensus_lengths(calls: pd.DataFrame) -> pd.Series:
    """Largest MEI_3P seen per subfamily. A PROXY, used only where no panel knows."""
    usable = calls[~calls.sentinel]
    return usable.groupby("subfamily")["MEI_3P"].max().rename("consensus_len")


def panel_consensus_lengths(
    subfamilies: "pd.Series | list[str]",
    *,
    panel_paths: "list[tuple[str, Path]]" = (),
    public_data_dir: str | Path | None = None,
    explicit_panels: "list[Path]" = (),
) -> tuple[pd.DataFrame, list[tuple[str, Path]]]:
    """Resolve every subfamily to a panel length; report how each label resolved."""
    panels = panel_paths or find_panel_paths(explicit_panels, public_data_dir=public_data_dir)
    labels = sorted({str(s) for s in subfamilies if str(s) and str(s) != "nan"})
    resolved = resolve_consensus_lengths(labels, panels)
    rows = []
    for label in labels:
        hit = resolved.get(label)
        rows.append({
            "subfamily": label,
            "consensus_len": hit.bp if hit else None,
            "consensus_target": hit.target if hit else "",
            "consensus_axis": hit.axis if hit else "",
            "consensus_basis": hit.basis if hit else "",
            "consensus_source": f"panel:{hit.basis}" if hit else PROXY_SOURCE,
        })
    return pd.DataFrame(rows).set_index("subfamily"), list(panels)


def annotate(calls: pd.DataFrame, lengths: pd.DataFrame | None = None) -> pd.DataFrame:
    """Join each call to a consensus length: panel record first, data proxy as fallback."""
    usable = calls[~calls.sentinel].copy()
    if lengths is None or lengths.empty:
        usable["consensus_len"] = usable["subfamily"].map(estimate_consensus_lengths(calls))
        usable["consensus_source"] = PROXY_SOURCE
    else:
        joined = usable.join(lengths, on="subfamily")
        proxy = estimate_consensus_lengths(calls)
        missing = joined["consensus_len"].isna()
        joined.loc[missing, "consensus_len"] = joined.loc[missing, "subfamily"].map(proxy)
        joined["consensus_source"] = joined["consensus_source"].fillna(PROXY_SOURCE)
        usable = joined
    usable["ends_at_consensus"] = usable["MEI_3P"] == usable["consensus_len"]
    usable["span_fraction"] = usable["MEI_SPAN"] / usable["consensus_len"]
    usable["overruns_axis"] = usable["MEI_3P"] > usable["consensus_len"]
    return usable


def summarise(calls: pd.DataFrame, usable: pd.DataFrame, min_full_span: int) -> dict:
    line1 = usable[usable.MEIFAMILY == "LINE1"]
    full = line1[line1.MEI_SPAN >= min_full_span]
    return {
        "genomes": int(calls["sample"].nunique()),
        "calls_total": len(calls),
        "calls_by_family": calls.MEIFAMILY.value_counts().to_dict(),
        "sentinels_by_family": calls[calls.sentinel].MEIFAMILY.value_counts().to_dict(),
        "line1_usable": len(line1),
        "line1_apparently_full": len(full),
        "line1_full_ending_at_consensus": int(full.ends_at_consensus.sum()),
        "line1_full_per_genome": full.groupby("sample").size().to_dict(),
        "axis_sources": line1.consensus_source.value_counts().to_dict(),
        "line1_overrunning_axis": int(line1.overruns_axis.sum()),
        "line1_overrun_bp_max": int((line1.MEI_3P - line1.consensus_len).max()) if len(line1) else 0,
        "line1_overrun_bp_median": float((line1.MEI_3P - line1.consensus_len)
                                         [line1.overruns_axis].median()) if line1.overruns_axis.any() else 0.0,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("vcf_dir")
    ap.add_argument("--out-dir", default="line1_span_out")
    ap.add_argument("--min-full-span", type=int, default=5500,
                    help="span at or above which a LINE-1 call is treated as apparently full length")
    ap.add_argument("--public-data-dir", default=None,
                    help="root holding retrotransposon_db/{dfam,full_consensus}; defaults to "
                         "$RTM_PUBLIC_DATA_DIR or the working directory")
    ap.add_argument("--mei-panel", action="append", default=[], type=Path,
                    help="explicit Dfam MEI panel FASTA (repeatable); overrides discovery")
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    calls = load_cohort(args.vcf_dir)
    lengths, panels = panel_consensus_lengths(
        calls["subfamily"],
        public_data_dir=args.public_data_dir,
        explicit_panels=list(args.mei_panel),
    )
    usable = annotate(calls, lengths)
    stats = summarise(calls, usable, args.min_full_span)
    lengths.reset_index().to_csv(out / "consensus_lengths.tsv", sep="\t", index=False)

    line1 = usable[usable.MEIFAMILY == "LINE1"]
    bands = [(20, 99), (100, 999), (1000, args.min_full_span - 1), (args.min_full_span, 10 ** 9)]
    rows = []
    for lo, hi in bands:
        sel = line1[line1.MEI_SPAN.between(lo, hi)]
        if len(sel):
            rows.append({"band": f"{lo}-{hi}" if hi < 10 ** 9 else f">={lo}",
                         "n": len(sel),
                         "ends_at_consensus": int(sel.ends_at_consensus.sum()),
                         "pct": round(100 * sel.ends_at_consensus.mean(), 1)})
    pd.DataFrame(rows).to_csv(out / "line1_span_bands.tsv", sep="\t", index=False)
    line1.to_csv(out / "line1_calls.tsv", sep="\t", index=False)

    print(f"genomes {stats['genomes']} | calls {stats['calls_total']} | {stats['calls_by_family']}")
    print(f"sentinels {stats['sentinels_by_family']}")
    print(f"\nconsensus axis: {len(panels)} panel(s) "
          f"{[str(path) for _kind, path in panels] or '(none found)'}")
    if not lengths.empty:
        unknown = lengths.index[lengths.consensus_len.isna()].tolist()
        print(f"  subfamilies resolved: {int(lengths.consensus_len.notna().sum())}"
              f"/{len(lengths)} | axes {sorted(set(lengths.consensus_axis) - {''})}"
              f" | bases {lengths.consensus_basis.value_counts().to_dict()}")
        if unknown:
            print(f"  not in any panel, using the data-derived proxy: {len(unknown)} "
                  f"{unknown[:10]}{'...' if len(unknown) > 10 else ''}")
    else:
        print("  no panel found: every length falls back to the data-derived proxy")
    print("\nLINE-1 with a usable footprint: "
          f"{stats['line1_usable']}")
    print(f"apparently full length (span >= {args.min_full_span}): {stats['line1_apparently_full']}")
    print(f"  of which the 3' end sits exactly at the consensus length: "
          f"{stats['line1_full_ending_at_consensus']}")
    if stats["line1_usable"]:
        print(f"LINE-1 calls whose 3' end is past their own consensus: "
              f"{stats['line1_overrunning_axis']}"
              + (f" (overrun median {stats['line1_overrun_bp_median']:.0f} bp, "
                 f"max {stats['line1_overrun_bp_max']} bp)" if stats["line1_overrunning_axis"] else ""))
        print(f"  by axis source: {stats['axis_sources']}")
        print("  (a small overrun is projection slack on the fragment-to-full map, not a "
              "span longer than the element; with the proxy axis this count was zero by "
              "construction)")
    print("\nby span band:")
    for row in rows:
        print(f"  {row['band']:>10s}  n={row['n']:5d}  ends at consensus {row['ends_at_consensus']:5d} "
              f"({row['pct']}%)")
    print("\nper genome, apparently full-length LINE-1:")
    for sample, n in sorted(stats["line1_full_per_genome"].items()):
        print(f"  {sample}  {n}")
    print(f"\nwrote {out}/line1_span_bands.tsv, {out}/line1_calls.tsv "
          f"and {out}/consensus_lengths.tsv")


if __name__ == "__main__":
    main()
