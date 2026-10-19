"""Can the cohort's own data estimate how often an SVA is 5'-transduced?

The marker map node ``b_trans_sva`` carries Damert 2009's figures: 220/2,398
TSD-carrying SVAs (9.2%) carry a transduced segment, those segments run 14-2,161 bp
(mean 223), and one source group (the MAST2 group) feeds >=32% of SVA-F.  This
script asks what our own insertion data can say about them, instead of repeating
the snippet.

5' transduction leaves one footprint in a DNA catalog: an external promoter reads
*upstream* of the source element, so the new copy carries sequence from the source
locus's upstream flank at its own 5' end.  The observable is therefore "does this
insertion contain bases before its element's 5' end".  Two of the three published
numbers are of that kind; the third (MAST2 >=32% of SVA-F) is a claim about *which
source element* fed the copy, which needs a donor, not a length.

Three arms, in decreasing order of what they can prove:

  A  floor test (call level, cohort-wide).  MEI_5P is a coordinate on the subfamily's
     panel consensus and MEI_SPAN the consensus-alignment span.  Both are bounded by
     the panel target, whose 5' end is coordinate 1, so if the axis is clamped no call
     can carry a 5'-transduced segment and the rate is not estimable from these fields
     -- a representational fact, not a biological one.  This arm tests the clamp with
     counts and compares the implied expectation from the rate the repo already stores.
  B  junction inventory (call level).  How many calls carry 5'-junction evidence at all
     (a split read on the element's 5' side, an alignment that reaches the element's own
     5' end).  These are the calls in which a 5' transduction could hide; the fields
     cannot separate a donor-derived upstream segment from filler sequence the caller
     never resolved, so this arm sizes the ambiguity rather than resolving it.  It also
     counts the fields that could name a donor, because the MAST2 figure needs one.
  C  read-level architecture (fixture scale).  At the SVA loci the per-read detail
     carries a consensus match (mei_start/mei_end) and a soft clip whose tip sits at the
     junction.  For a read whose match starts at the element's 5' end, the clip outside
     that match is the only place a transduced segment can be, so the clip length is an
     upper bound on it.  This arm shows whether the channel exists in our data, at how
     many reads, and at what resolution.

Definitions, stated because they change the answer:

  element 5' end    coordinate 1 of the subfamily's consensus axis.  A call's MEI_5P is
                    its alignment's start on that axis, so MEI_5P == 1 means the read
                    evidence reaches the element's own 5' end -- exactly where a
                    transduced upstream segment would abut and stay unreported.
  5' extension      bases upstream of coordinate 1, i.e. MEI_SPAN > axis length with the
                    alignment anchored at 1.  Not expressible in a coordinate system
                    that starts at 1; the span test is the part that is not circular.
  X5                the read-level candidate.  The L1 arm of this test (see
                    l1_readthrough_probe.py) computes X = clip_len - consensus bases in
                    the clip - tail run, so a 3'-junction clip has two measurable
                    subtrahends.  At a 5'-junction clip both are zero by construction:
                    the caller's match begins *outside* the clip rather than inside it,
                    and the only run the fields record is the 3' poly(A) run.  So
                    X5 == clip_len and every 5' clip is "unexplained"; X5 is an upper
                    bound on a transduced segment and can never be positively attributed.

Caveats that are load-bearing:

  - Arm A tests the fields, not the biology.  Zero span overruns mean the catalog cannot
    express the signal, NOT that 5' transduction is absent.
  - Arm C runs on the SVA loci among the 22 fixtures, in one sample (HG00100).  Three of
    them are SVA loci and one of the three has reads at the element's own 5' end, so arm C
    rests on a single SVA 5' junction.  Nothing in arm C is a rate; it is a demonstration
    of what the read level can resolve, with n stated.
  - The published 9.2% is over 2,398 *TSD-carrying* SVAs, so the nearest denominator here
    is the SVA calls that carry a TSD, reported alongside the full call count.
  - These calls come from one caller.  Every count below is a statement about calls.

Usage:
    python sva_5prime_transduction_probe.py [--vcf-dir nested_analysis/callsets]
                                            [--calls sva_5prime_out/sva_calls.tsv]
                                            [--loci rtm-longread-cohort/tests/fixtures/loci]
                                            [--out-dir sva_5prime_out]
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import re
import sys
from pathlib import Path

import pandas as pd

# The sibling package holds the call loader, the panel index and the marker module.  The
# L1 script of this pair imports from it the same way, so both run with plain python3.
sys.path.insert(0, str(Path(__file__).resolve().parent / "retrotransposon-miner" / "src"))

# The repo's own literature constants; imported rather than retyped so the two cannot
# drift.  They live with the marker module, which is not on every branch (the L1 script
# made the same fallback): the stored values are quoted here for the case where the
# package is not importable from this working directory.
try:  # pragma: no cover - import convenience
    from retro_miner.mei_markers import (  # type: ignore[attr-defined]
        SVA_MAST2_GROUP_FRACTION_SVA_F,
        SVA_TRANSDUCTION_RATE,
    )
except ImportError:  # pragma: no cover
    SVA_TRANSDUCTION_RATE = 0.08
    SVA_MAST2_GROUP_FRACTION_SVA_F = 0.32

#: Damert 2009 as recorded in the marker catalog: 220 of 2,398 TSD-carrying SVAs.
DAMERT_TRANSDUCED, DAMERT_TSD_DENOMINATOR = 220, 2398

#: Segment lengths in Damert 2009 (Fig. 4), for the read-level channel comparison.
DAMERT_SEGMENT_RANGE = (14, 2161)
DAMERT_SEGMENT_MEAN = 223

NUMERIC_CALL_COLS = ("MEI_SPAN", "MEI_5P", "MEI_3P", "POLYA_MIN_BP")
SUPPORT_KEYS = ("SR_L", "SR_R", "DPE_L", "DPE_R", "MEI_MAPPED", "polyA_MAPPED", "VNTR_MAPPED")
#: Fields that name an element.  All of them point outside this call set (an external
#: catalog's record id); none names the source element that fed the copy.
EXTERNAL_ID_FIELDS = ("G1K", "LR", "KNOWN_ID", "KNOWN_SRC")
#: Fields that would state how long the *inserted copy* is.  A transduced copy is longer
#: than its element, so without one of these the excess has nowhere to be written.
LENGTH_FIELDS = ("SVLEN", "INSLEN", "SVLEN_BP", "INSERTION_LEN")
DETAIL_NUMERIC = ("mei_start", "mei_end", "soft_clip_len", "clip_len", "clip_poly_at_run",
                  "anchor_poly_at_run", "soft_clip_pos")

#: Coordinate of the element's own 5' end on a panel consensus axis.
FIVE_PRIME_END = 1


def support_value(support: pd.Series, key: str) -> pd.Series:
    """A ``KEY:<n>`` field out of the caller's pipe-separated SUPPORT string."""
    return pd.to_numeric(
        support.str.extract(rf"(?:^|\|){re.escape(key)}:(\d+)")[0], errors="coerce")


def load_calls(path: str | Path) -> pd.DataFrame:
    calls = pd.read_csv(path, sep="\t", low_memory=False)
    for col in NUMERIC_CALL_COLS:
        if col in calls.columns:
            calls[col] = pd.to_numeric(calls[col], errors="coerce")
    if "subfamily" not in calls.columns and "MEISUBFAMILY" in calls.columns:
        calls["subfamily"] = calls["MEISUBFAMILY"].fillna("").str.split("#").str[0]
    support = calls.get("SUPPORT", pd.Series("", index=calls.index)).fillna("")
    for key in SUPPORT_KEYS:
        calls[key] = support_value(support, key).fillna(0)
    calls["polyA_side"] = support.str.extract(r"polyA_side:(\w)")[0]
    calls["orient"] = calls.get("ORIENT", pd.Series("", index=calls.index)).fillna("")
    calls["at_5p_end"] = calls["MEI_5P"] == FIVE_PRIME_END
    calls["near_5p_end"] = calls["MEI_5P"].between(FIVE_PRIME_END, 10)
    calls["unresolved_5p"] = calls["MEI_5P"] < 0
    calls["has_tsd"] = calls.get("TSD", pd.Series("", index=calls.index)).fillna("").astype(str).str.len() > 0
    # The element's own 5' junction is the genomic left end for a '+' call and the right
    # end for a '-' call: a '-' SVA's upstream flank is the reference downstream of it.
    calls["support_5p"] = [
        row[f"SR_{'L' if row['orient'] == '+' else 'R'}"] for _, row in calls.iterrows()]
    calls["dpe_5p"] = [
        row[f"DPE_{'L' if row['orient'] == '+' else 'R'}"] for _, row in calls.iterrows()]
    return calls


def build_sva_calls(vcf_dir: str | Path) -> tuple[pd.DataFrame, dict]:
    """Every SVA call in the cohort's VCFs, with the caller's fields as columns.

    Built here rather than shipped, so the probe's inputs are reproducible from the
    cohort call sets.  ``retro_miner.gene_annotation.read_vcf`` is the loader the
    sibling analyses use; a call file whose schema has no SUPPORT cannot answer arm B
    and is counted instead of silently mixed in.
    """
    from retro_miner.gene_annotation import read_vcf

    counts: dict = {"families_all": {}, "files": {}, "dropped_no_support": 0}
    frames: list[pd.DataFrame] = []
    for path in sorted(glob.glob(os.path.join(str(vcf_dir), "*.vcf"))):
        name = os.path.basename(path)
        _, records = read_vcf(path)
        rows = []
        for record in records:
            info = dict(record.info)
            counts["families_all"][info.get("MEIFAMILY", "")] = (
                counts["families_all"].get(info.get("MEIFAMILY", ""), 0) + 1)
            if info.get("MEIFAMILY") != "SVA":
                continue
            if "SUPPORT" not in info:
                counts["dropped_no_support"] += 1
                continue
            rows.append({"chrom": record.chrom, "pos": record.pos,
                         "sample": name.split("_")[0], "vcf": name, **info})
        if rows:
            frames.append(pd.DataFrame(rows))
            counts["files"][name] = len(rows)
        else:
            counts["files"][name] = 0
    if not frames:
        raise SystemExit(f"no SVA call with a SUPPORT field in {vcf_dir}")
    return pd.concat(frames, ignore_index=True), counts


def panel_axis(labels: list[str]) -> tuple[dict[str, tuple[int, str]], dict]:
    """Consensus length per subfamily from the Dfam panel, or a data-derived fallback.

    The number travels with its basis, because the caller's coordinates are on the
    *full* axis (a fragment record projects onto its family canonical) -- a length read
    off the wrong axis is how a span appears to exceed its own consensus.
    """
    axis: dict[str, tuple[int, str]] = {}
    notes: dict = {"panel": "none", "basis": {}}
    try:  # pragma: no cover - depends on the working directory
        from retro_miner.mei_panel_index import (
            find_panel_paths,
            resolve_consensus_lengths,
        )
    except ImportError:  # pragma: no cover
        notes["panel"] = "unavailable (retro_miner not importable)"
        return axis, notes

    # An empty panel list is not an error here: resolve_consensus_lengths answers None for
    # every label, and every caller of this function treats None as "no panel knows".
    panels = find_panel_paths()
    resolved = resolve_consensus_lengths([f"{label}_full" for label in labels], panels)
    notes["panel"] = ",".join(f"{kind}:{path.name}" for kind, path in panels) or "none"
    for label in labels:
        hit = resolved.get(f"{label}_full")
        if hit is not None:
            axis[label] = (hit.bp, f"panel_{hit.axis}:{hit.basis}:{hit.target}")
            notes["basis"][label] = axis[label][1]
    return axis, notes


def axis_table(calls: pd.DataFrame, axis: dict[str, tuple[int, str]]) -> pd.DataFrame:
    """Arm A. Per-subfamily axis length, how many calls sit at its 5' end, any overrun.

    Without a panel the axis is estimated the way span_work/line1_span_support.py did it
    (the largest MEI_3P observed for the subfamily), which makes "at the ceiling" partly
    circular; the span test is the part that survives that.
    """
    calls = calls.copy()
    fallback = calls.groupby("subfamily")["MEI_3P"].max().to_dict()
    calls["axis_len"] = [axis.get(s, (int(fallback.get(s, 0)), "data_derived_max_mei_3p"))[0]
                         for s in calls["subfamily"]]
    calls["axis_basis"] = [axis.get(s, (0, "data_derived_max_mei_3p"))[1]
                           for s in calls["subfamily"]]
    calls["span_over_axis"] = calls["MEI_SPAN"] > calls["axis_len"]
    calls["anchored_at_5p"] = calls["at_5p_end"] & calls["span_over_axis"]

    per_subfamily = (calls.groupby(["subfamily", "axis_basis"])
                     .agg(n=("MEI_5P", "size"),
                          axis_len=("axis_len", "max"),
                          n_at_5p_end=("at_5p_end", "sum"),
                          n_within_10bp=("near_5p_end", "sum"),
                          n_span_over_axis=("span_over_axis", "sum"),
                          max_mei_span=("MEI_SPAN", "max"),
                          median_5p=("MEI_5P", "median"),
                          median_3p=("MEI_3P", "median"),
                          n_with_tsd=("has_tsd", "sum"))
                     .reset_index()
                     .sort_values("n", ascending=False))
    per_subfamily["pct_at_5p_end"] = (
        100 * per_subfamily.n_at_5p_end / per_subfamily.n).round(1)
    return calls, per_subfamily


def junction_table(calls: pd.DataFrame) -> pd.DataFrame:
    """Arm B. Which calls hold 5'-junction evidence at all, and which name a donor."""
    external = pd.Series(0, index=calls.index)
    for field in EXTERNAL_ID_FIELDS:
        if field in calls.columns:
            external = external | calls[field].fillna("").astype(str).str.len().gt(0)
    tail = calls.get("polyA_MAPPED", pd.Series(0, index=calls.index)).fillna(0)
    rows = [
        {"class": "calls", "n": len(calls), "pct": 100.0},
        {"class": "5' end resolved by the caller (MEI_5P >= 1)",
         "n": int((~calls.unresolved_5p).sum()),
         "pct": round(100 * (~calls.unresolved_5p).mean(), 1)},
        {"class": "alignment reaches the element's 5' end (MEI_5P == 1)",
         "n": int(calls.at_5p_end.sum()),
         "pct": round(100 * calls.at_5p_end.mean(), 1)},
        {"class": "alignment within 10 bp of the element's 5' end",
         "n": int(calls.near_5p_end.sum()),
         "pct": round(100 * calls.near_5p_end.mean(), 1)},
        {"class": "split read on the element's 5' side (SR_L for '+', SR_R for '-')",
         "n": int((calls.support_5p >= 1).sum()),
         "pct": round(100 * (calls.support_5p >= 1).mean(), 1)},
        {"class": "discordant pair on the element's 5' side",
         "n": int((calls.dpe_5p >= 1).sum()),
         "pct": round(100 * (calls.dpe_5p >= 1).mean(), 1)},
        {"class": "5' end reached AND 5'-side split read",
         "n": int((calls.at_5p_end & (calls.support_5p >= 1)).sum()),
         "pct": round(100 * (calls.at_5p_end & (calls.support_5p >= 1)).mean(), 1)},
        {"class": "3'-side tail read mapped (the measured run sits there)",
         "n": int((tail >= 1).sum()),
         "pct": round(100 * (tail >= 1).mean(), 1)},
        {"class": "call carries a TSD sequence (Damert's denominator)",
         "n": int(calls.has_tsd.sum()),
         "pct": round(100 * calls.has_tsd.mean(), 1)},
        {"class": "call matches an external catalog record (G1K/LR/KNOWN_ID)",
         "n": int(external.sum()),
         "pct": round(100 * external.mean(), 1)},
        {"class": "call names the *source* element that fed it",
         "n": 0, "pct": 0.0},
    ]
    return pd.DataFrame(rows)


def read_level_5p(loci: str | Path, *, five_prime_window: int = 25,
                  junction_window: int = 50) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Arm C. 5'-anchored SVA reads and the clip that could carry a transduced segment.

    A row counts when the caller's per-read consensus match starts within
    ``five_prime_window`` bp of the element's 5' end (coordinate 1) and the read's soft
    clip tip sits within ``junction_window`` bp of the locus's expected breakpoint.
    Identical rows are collapsed first: the detail files can list one read twice, and a
    row is not a read, and the per-locus funnel comes back with the rows because how many
    SVA loci have any read at the element's 5' end at all is part of the answer.

    The detail files carry the same read under both halves of the design (``sample`` is
    ``disease`` on one row and ``control`` on the other, everything else identical), so
    duplicates are collapsed on every column except that one.
    """
    rows: list[dict] = []
    funnel: list[dict] = []
    for path in sorted(glob.glob(os.path.join(str(loci), "*/supporting_reads_detail.mei.tsv"))):
        locus_dir = Path(path).parent
        if not (locus_dir / "gold_locus.tsv").exists():
            continue
        detail = pd.read_csv(path, sep="\t", low_memory=False)
        detail = detail.drop_duplicates(
            subset=[c for c in detail.columns if c != "sample"]).reset_index(drop=True)
        if "mei_hit" not in detail.columns:
            continue
        for col in DETAIL_NUMERIC:
            if col in detail.columns:
                detail[col] = pd.to_numeric(detail[col], errors="coerce")
        hit = detail[detail["mei_hit"].astype(str).str.lower().isin({"true", "1"})]
        sva = hit[hit.get("mei_target", pd.Series("", index=hit.index)).fillna("")
                  .str.startswith("SVA")].copy()  # panel labels are SVA_<subfamily>_full
        if sva.empty:
            continue
        manifest_path = locus_dir / "manifest.json"
        manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
        breakpoint = manifest.get("expected_breakpoint")
        anchored = sva[sva["mei_start"].between(FIVE_PRIME_END, five_prime_window)].copy()
        if breakpoint is None:
            anchored["at_junction"] = False
        else:
            anchored["at_junction"] = anchored["soft_clip_pos"].sub(
                int(breakpoint)).abs() <= junction_window
        junction = anchored[anchored["at_junction"] & (anchored["soft_clip_len"] > 0)]
        funnel.append({
            "locus": locus_dir.name,
            "sva_target_rows": len(sva),
            "five_prime_anchored_rows": len(anchored),
            "junction_clip_rows": len(junction),
            "junction_reads": int(junction.read_name.nunique()) if "read_name" in junction else 0,
            "clip_min": int(junction.soft_clip_len.min()) if len(junction) else 0,
            "clip_median": float(junction.soft_clip_len.median()) if len(junction) else 0.0,
            "clip_max": int(junction.soft_clip_len.max()) if len(junction) else 0,
        })
        for _, row in junction.iterrows():
            rows.append({
                "locus": locus_dir.name,
                "read_name": row.get("read_name", ""),
                "evidence_type": row.get("evidence_type", ""),
                "mei_target": row.get("mei_target", ""),
                "mei_start": int(row["mei_start"]),
                "mei_end": int(row["mei_end"]),
                "clip_bp": int(row["soft_clip_len"]),
                "clip_side": row.get("soft_clip_side", ""),
                "clip_tip": (int(row["soft_clip_pos"]) if pd.notna(row["soft_clip_pos"]) else None),
                # The two subtrahends the L1 arm had, and whether this row has them:
                # element bases inside the clip (the match starts at the clip, not in it)
                # and a landmark run (clip_poly_at_run records the 3' poly(A) run only).
                "element_bp_in_clip": 0,
                "clip_landmark_run_bp": (int(row["clip_poly_at_run"])
                                         if pd.notna(row.get("clip_poly_at_run")) else 0),
                "X5_bp": int(row["soft_clip_len"]),
                "anchor_run_bp": (int(row["anchor_poly_at_run"])
                                  if pd.notna(row.get("anchor_poly_at_run")) else 0),
            })
    return pd.DataFrame(rows), pd.DataFrame(funnel)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf-dir", default="nested_analysis/callsets")
    ap.add_argument("--calls", default="sva_5prime_out/sva_calls.tsv",
                    help="SVA call table; built from --vcf-dir when it does not exist")
    ap.add_argument("--loci", default="rtm-longread-cohort/tests/fixtures/loci")
    ap.add_argument("--out-dir", default="sva_5prime_out")
    ap.add_argument("--five-prime-window", type=int, default=25,
                    help="a read is 5'-anchored when its consensus match starts within "
                         "this many bp of the element's 5' end")
    ap.add_argument("--junction-window", type=int, default=50,
                    help="bp around the locus's expected breakpoint in which a clip tip "
                         "counts as a junction clip")
    ap.add_argument("--rebuild-calls", action="store_true",
                    help="rebuild the call table from the VCFs even if it exists")
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    counts: dict = {}
    if args.rebuild_calls or not Path(args.calls).exists():
        calls, counts = build_sva_calls(args.vcf_dir)
        calls.to_csv(args.calls, sep="\t", index=False)
    calls = load_calls(args.calls)

    labels = sorted({s for s in calls["subfamily"] if s})
    axis, axis_notes = panel_axis(labels)
    calls, per_subfamily = axis_table(calls, axis)
    junction = junction_table(calls)
    per_subfamily.to_csv(out / "five_prime_axis_by_subfamily.tsv", sep="\t", index=False)
    junction.to_csv(out / "junction_support.tsv", sep="\t", index=False)

    n = len(calls)
    at_5p = int(calls.at_5p_end.sum())
    overruns = int(calls.span_over_axis.sum())
    with_tsd = int(calls.has_tsd.sum())
    rate, mast2 = SVA_TRANSDUCTION_RATE, SVA_MAST2_GROUP_FRACTION_SVA_F
    damert = DAMERT_TRANSDUCED / DAMERT_TSD_DENOMINATOR
    sva_f = calls[calls["subfamily"].str.startswith("SVA_F")]
    expected = {
        "stored constant on all calls": rate * n,
        "stored constant on TSD-carrying calls": rate * with_tsd,
        "Damert ratio on TSD-carrying calls": damert * with_tsd,
    }

    print(f"calls {n} | subfamilies {calls.subfamily.nunique()} "
          f"| samples {calls.get('sample', pd.Series(dtype=str)).nunique()}")
    if counts:
        print(f"built from {args.vcf_dir}: files with SVA calls "
              f"{sum(1 for v in counts['files'].values() if v)}")
        print(f"  families in that directory: {counts['families_all']}; SVA calls kept "
              f"{n}, dropped for a missing SUPPORT {counts['dropped_no_support']}")
    print(f"panel used for the axis: {axis_notes['panel']}")
    if any(b.startswith("data_derived") for b in axis_notes.get("basis", {}).values()) or not axis:
        print("  NOTE: no panel on disk for some labels -- those axes are the largest "
              "MEI_3P observed here, which makes 'at the ceiling' partly circular.")

    print("\nARM A . can any call carry sequence before its element's 5' end?")
    print(f"  calls whose consensus span exceeds the panel axis length "
          f"(MEI_SPAN > axis): {overruns}")
    print(f"  calls with both the 5' anchor and a span past the axis (the exact "
          f"5'-extension shape): {int(calls.anchored_at_5p.sum())}")
    print(f"  calls whose alignment reaches the element's 5' end (MEI_5P == 1): "
          f"{at_5p} ({100 * at_5p / n:.1f}%)")
    print(f"  calls with no resolved 5' coordinate (MEI_5P < 0): "
          f"{int(calls.unresolved_5p.sum())}")
    print(f"  fields in this call table that state the inserted copy's length: "
          f"{sorted(f for f in LENGTH_FIELDS if f in calls.columns) or 'none'} "
          f"(MEI_SPAN is the consensus-matched span, bounded by the axis)")
    print("  (a 5'-transduced segment sits upstream of coordinate 1, so it cannot appear "
          "in MEI_5P, in an axis-bounded span, or in the call's length field)")
    print(f"  if SVA 5' transduction ran at the stored rate {rate:.2f}, this call set "
          f"would hold {expected['stored constant on all calls']:.0f} transduced copies "
          f"({expected['stored constant on TSD-carrying calls']:.0f} of the {with_tsd} "
          f"TSD-carrying calls); at Damert's own ratio {damert:.3f}, "
          f"{expected['Damert ratio on TSD-carrying calls']:.0f}")
    print(f"  SVA-F calls: {len(sva_f)}; at the stored MAST2 figure {mast2:.2f} of SVA-F, "
          f"{mast2 * len(sva_f):.0f} of them would belong to one source group")
    print("\n  per-subfamily 5' floor (top rows):")
    print(per_subfamily.head(10).to_string(index=False))

    print("\nARM B . 5'-junction evidence in the call fields")
    print(junction.to_string(index=False))

    reads, funnel = read_level_5p(args.loci, five_prime_window=args.five_prime_window,
                                  junction_window=args.junction_window)
    reads.to_csv(out / "read_level_5p_clips.tsv", sep="\t", index=False)
    print(f"\nARM C . read-level architecture at the SVA fixture loci "
          f"({args.five_prime_window} bp 5' window, {args.junction_window} bp junction "
          f"window)")
    if not funnel.empty:
        print("  per-locus funnel (rows, then distinct reads):")
        print(funnel.to_string(index=False))
    if reads.empty:
        print("  no read-level row carries both a 5'-anchored consensus match and a "
              "junction clip")
    else:
        print(f"  reads at an SVA 5' junction: {len(reads)} rows at {reads.locus.nunique()} "
              f"loci ({reads.read_name.nunique()} distinct read names)")
        print(f"  clip length range: {int(reads.clip_bp.min())}-{int(reads.clip_bp.max())} bp "
              f"(median {reads.clip_bp.median():.0f})")
        print(f"  rows with a landmark run inside the clip: "
              f"{int((reads.clip_landmark_run_bp > 0).sum())}; rows whose *other* "
              f"junction carries the recorded run: {int((reads.anchor_run_bp > 0).sum())}")
        print(f"  X5 == clip_bp for every row (nothing to subtract), so the candidate "
              f"upper bound is {int(reads.clip_bp.max())} bp in the largest clip and "
              f"0 bp is attributable to a donor")
        seg_lo, seg_hi = DAMERT_SEGMENT_RANGE
        print(f"  Damert's segments run {seg_lo}-{seg_hi} bp (mean {DAMERT_SEGMENT_MEAN}): "
              f"{int((reads.clip_bp >= seg_lo).sum())} of {len(reads)} clips are long "
              f"enough to host the short end of that range, "
              f"{int((reads.clip_bp >= DAMERT_SEGMENT_MEAN).sum())} to host its mean")
    print(f"\nwrote {out}/five_prime_axis_by_subfamily.tsv, {out}/junction_support.tsv, "
          f"{out}/read_level_5p_clips.tsv"
          + (f", {args.calls}" if counts else ""))
    print("NOTE: arm C covers one SVA 5' junction in one sample.  Nothing here is a rate.")


if __name__ == "__main__":
    main()
