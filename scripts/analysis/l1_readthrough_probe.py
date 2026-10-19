"""Can the cohort's own data estimate how often L1 transcription reads through?

The map node b_l1term carries a snippet saying 21% of L1 transcription terminates
downstream of the element.  That is a statement about where RNA polymerase and the
polyadenylation machinery stopped.  This script asks what our own data can say about
it, instead of repeating the snippet.

Read-through leaves exactly one footprint in a DNA catalog: when the extended
transcript is the one that retrotransposes, the new insertion carries a *transduced*
segment -- sequence from the donor locus's downstream region -- between the element's
own 3' end and the poly(A) tail.  So the observable is "does this insertion contain
bases past its element's 3' end".

Three arms, in decreasing order of what they can prove:

  A  ceiling test (call level, cohort-wide).  MEI_3P is a coordinate on the
     subfamily's full-length panel reference and MEI_SPAN the consensus-alignment
     span.  Both are bounded by the panel target, so if the axis is clamped no call can
     carry a transduced tail and the rate is not estimable from these fields -- a
     representational fact, not a biological one.  This arm tests the clamp with
     counts and compares the implied expectation from the literature band the repo
     already stores (L1 3' transduction 8.8-15%).
  B  junction inventory (call level).  How many calls carry 3'-junction evidence at
     all (a mapped poly-A/T read, a reported tail, a 3' end at the axis ceiling).
     These are the calls in which read-through could hide; the fields cannot separate
     a tail placed at the element's own signal from one placed after a transduced
     segment, so this arm sizes the ambiguity rather than resolving it.
  C  read-level architecture (fixture scale).  Where a read's junction clip carries
     both a consensus match and a poly-A/T run, the caller's own fields give the
     element's contribution (mei_end - mei_start + 1) and the tail (clip_poly_at_run),
     so the unexplained remainder of the clip is measurable:
         X = clip_len - (element bases in the clip) - (tail run)
     X is candidate transduced sequence sitting between the element's 3' end and the
     tail.  This arm shows whether the channel exists in our data, on how much of the
     cohort, and at what resolution.

Definitions, stated because they change the answer:

  axis length L(s)   largest MEI_3P observed for subfamily s across the cohort.  A
                     proxy for the Dfam panel length, data-derived (same convention as
                     span_work/line1_span_support.py), so "at the ceiling" is partly
                     circular and is reported as such.
  span overrun       MEI_SPAN > L(s).  MEI_3P > L(s) is impossible by construction of
                     L(s), so only the span test carries information.
  X                  clip bases left after removing the consensus-matched part and the
                     tail run.  Alignment slack, mismatches and clip trimming all land
                     in X, so X is an upper bound on transduced bases and a positive X
                     is a candidate, never a call.

A measurement this script does NOT make, and why.  A naive alternative is to measure,
inside a junction soft clip, the distance from the alignment boundary to the first
poly-A/T run.  That cannot discriminate: the transcript's tail abuts the flanking
reference in BOTH the terminated and the read-through case (order is always
[element][(donor)][tail][flank]), so the tail's position relative to the junction is
the same either way.  Only the identity of the sequence between the element's 3' end
and the tail separates the two, which is why arm C needs the consensus match.

Caveats that are load-bearing:

  - Arm A tests the field, not the biology.  Zero span overruns mean the catalog cannot
    express the signal, NOT that read-through is absent.
  - Arm C runs on the 22 locus fixtures under tests/fixtures/loci/, which is one LINE-1
    locus.  Nothing in arm C is a rate; it is a demonstration of what the read level can
    resolve, with n stated.
  - The published rate is for *detected 3' transduction*, a filtered lower bound on
    transcription read-through: read-through must also reach a usable downstream
    poly(A) site, be the transcript that retrotransposes, and be long enough to detect.
    A 21% read-through fraction and an 8.8-15% transduction rate are therefore not in
    tension, and neither confirms the other.
  - These calls come from one caller.  Every count below is a statement about calls.

Usage:
    python l1_readthrough_probe.py [--calls span_work/out/line1_calls.tsv]
                                   [--loci rtm-longread-cohort/tests/fixtures/loci]
                                   [--out-dir l1_readthrough_out]
"""

from __future__ import annotations

import argparse
import glob
import os
from pathlib import Path

import pandas as pd

# The repo's own literature constants; imported rather than retyped so the two cannot
# drift.  Falls back to the stored numbers if the package is not importable here.
try:  # pragma: no cover - import convenience
    from retro_miner.mei_markers import L1_TRANSDUCTION_RATE_RANGE
except Exception:  # pragma: no cover
    L1_TRANSDUCTION_RATE_RANGE = (0.088, 0.15)

NUMERIC_CALL_COLS = ("MEI_SPAN", "MEI_5P", "MEI_3P", "POLYA_MIN_BP", "MEI_MAPPED",
                     "polyA_MAPPED")
DETAIL_NUMERIC = ("mei_start", "mei_end", "clip_len", "clip_poly_at_run")


def family_from_target(target: str) -> str:
    t = (target or "").strip()
    if t.startswith("L1"):
        return "LINE1"
    if t.startswith("Alu"):
        return "ALU"
    if t.startswith("SVA"):
        return "SVA"
    return ""


def load_calls(path: str | Path) -> pd.DataFrame:
    calls = pd.read_csv(path, sep="\t", low_memory=False)
    for col in NUMERIC_CALL_COLS:
        if col in calls.columns:
            calls[col] = pd.to_numeric(calls[col], errors="coerce")
    if "subfamily" not in calls.columns and "MEISUBFAMILY" in calls.columns:
        calls["subfamily"] = calls["MEISUBFAMILY"].fillna("").str.split("#").str[0]
    return calls


def ceiling_table(calls: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Arm A. Per-subfamily axis length, how many calls sit on it, and any overrun."""
    calls = calls.copy()
    axis = calls.groupby("subfamily")["MEI_3P"].max().rename("axis_len")
    calls = calls.join(axis, on="subfamily")
    calls["at_3p_ceiling"] = calls["MEI_3P"] == calls["axis_len"]
    calls["reaches_3p_end"] = calls["MEI_3P"] >= calls["axis_len"] - 50
    calls["span_overrun"] = calls["MEI_SPAN"] > calls["axis_len"]

    per_subfamily = (calls.groupby("subfamily")
                     .agg(n=("MEI_3P", "size"),
                          axis_len=("axis_len", "max"),
                          n_at_ceiling=("at_3p_ceiling", "sum"),
                          n_reach_3p_end=("reaches_3p_end", "sum"),
                          n_span_overrun=("span_overrun", "sum"),
                          max_mei_span=("MEI_SPAN", "max"),
                          n_span_at_or_over_axis=("MEI_SPAN", "size"),
                          median_5p=("MEI_5P", "median"),
                          median_polya_bp=("POLYA_MIN_BP", "median"))
                     .reset_index()
                     .sort_values("n", ascending=False))
    per_subfamily["pct_at_ceiling"] = (
        100 * per_subfamily.n_at_ceiling / per_subfamily.n).round(1)
    per_subfamily = per_subfamily.drop(columns=["n_span_at_or_over_axis"])
    return calls, per_subfamily


def junction_table(calls: pd.DataFrame) -> pd.DataFrame:
    """Arm B. Which calls hold 3'-junction evidence at all, by architecture class."""
    mapped = calls.get("polyA_MAPPED", pd.Series(0, index=calls.index)).fillna(0)
    tail_min = calls.get("POLYA_MIN_BP", pd.Series(0, index=calls.index)).fillna(0)
    side = calls.get("polyA_side", pd.Series("", index=calls.index)).fillna("")
    rows = [
        {"class": "calls", "n": len(calls), "pct": 100.0},
        {"class": "any poly-A/T read mapped",
         "n": int((mapped >= 1).sum()),
         "pct": round(100 * (mapped >= 1).mean(), 1)},
        {"class": "reported tail >= 20 bp",
         "n": int((tail_min >= 20).sum()),
         "pct": round(100 * (tail_min >= 20).mean(), 1)},
        {"class": "poly-A/T run on the right of the window (polyA_side=R)",
         "n": int((side == "R").sum()),
         "pct": round(100 * (side == "R").mean(), 1)},
        {"class": "3' end within 50 bp of the axis ceiling",
         "n": int(calls["reaches_3p_end"].sum()),
         "pct": round(100 * calls["reaches_3p_end"].mean(), 1)},
        {"class": "3' end at the axis ceiling AND a mapped poly-A/T read",
         "n": int((calls["at_3p_ceiling"] & (mapped >= 1)).sum()),
         "pct": round(100 * (calls["at_3p_ceiling"] & (mapped >= 1)).mean(), 1)},
    ]
    return pd.DataFrame(rows)


def read_level_x(detail: pd.DataFrame, locus: str) -> list[dict]:
    """Arm C. X = clip bases that are neither consensus-matched nor tail.

    Only reads whose consensus match reaches the furthest coordinate observed for its
    own target in this locus file are kept: a match that stops short of the target's end
    cannot speak about the junction, and a 5'-end subfamily target never does.
    """
    if detail.empty or "mei_hit" not in detail.columns:
        return []
    work = detail.copy()
    for col in DETAIL_NUMERIC:
        if col in work.columns:
            work[col] = pd.to_numeric(work[col], errors="coerce")
    work["mei_target"] = work.get("mei_target", pd.Series("", index=work.index)).fillna("")
    hit = work[work.get("mei_hit", False).astype(str).str.lower().isin({"true", "1"})]
    if hit.empty:
        return []
    axis_end = hit.groupby("mei_target")["mei_end"].max()
    rows: list[dict] = []
    for _, r in hit.iterrows():
        if pd.isna(r.get("mei_end")) or pd.isna(r.get("clip_len")) or r["clip_len"] <= 0:
            continue
        if pd.isna(axis_end.get(r["mei_target"])) or r["mei_end"] < axis_end[r["mei_target"]]:
            continue
        if pd.isna(r.get("clip_poly_at_run")) or r["clip_poly_at_run"] <= 0:
            continue
        element_bp = r["mei_end"] - r["mei_start"] + 1
        tail_bp = r["clip_poly_at_run"]
        non_tail_clip = max(int(r["clip_len"]) - int(tail_bp), 0)
        rows.append({
            "locus": locus,
            "family": family_from_target(r["mei_target"]),
            "mei_target": r["mei_target"],
            "read_name": r.get("read_name", ""),
            "clip_len": int(r["clip_len"]),
            "element_bp_in_clip": int(element_bp),
            "tail_run_bp": int(tail_bp),
            "X_bp": int(r["clip_len"] - element_bp - tail_bp),
            # X is only about the junction when the caller's own match already accounts
            # for the clip's non-tail bases.  Where the match is a short seed, X is
            # unmatched element sequence and must not be read as a transduced segment.
            "match_explains_clip": bool(element_bp >= 0.8 * non_tail_clip),
        })
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--calls", default="span_work/out/line1_calls.tsv")
    ap.add_argument("--loci", default="rtm-longread-cohort/tests/fixtures/loci")
    ap.add_argument("--out-dir", default="l1_readthrough_out")
    ap.add_argument("--candidate-x-bp", type=int, default=20,
                    help="X at or above which a read is listed as a read-through candidate")
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    calls = load_calls(args.calls)
    calls, per_subfamily = ceiling_table(calls)
    junction = junction_table(calls)
    per_subfamily.to_csv(out / "axis_ceiling_by_subfamily.tsv", sep="\t", index=False)
    junction.to_csv(out / "junction_support.tsv", sep="\t", index=False)

    n = len(calls)
    at_ceiling = int(calls.at_3p_ceiling.sum())
    span_overruns = int(calls.span_overrun.sum())
    lo, hi = L1_TRANSDUCTION_RATE_RANGE
    expected = (n * lo, n * hi)

    rows: list[dict] = []
    loci_covered = 0
    detail_files = sorted(glob.glob(os.path.join(args.loci, "*/supporting_reads_detail.mei.tsv")))
    for path in detail_files:
        locus = Path(path).parent.name
        rows.extend(read_level_x(pd.read_csv(path, sep="\t", low_memory=False), locus))
        loci_covered += 1
    reads = pd.DataFrame(rows)
    reads.to_csv(out / "read_level_junction_x.tsv", sep="\t", index=False)
    resolved = reads[reads.match_explains_clip] if not reads.empty else reads

    print(f"calls {n} | subfamilies {calls.subfamily.nunique()} "
          f"| genomes {calls.get('sample', pd.Series(dtype=str)).nunique()}")
    print("\nARM A . does any call carry sequence past its element's 3' end?")
    print(f"  calls whose 3' end sits exactly on the subfamily axis ceiling: "
          f"{at_ceiling} ({100 * at_ceiling / n:.1f}%)")
    print(f"  calls whose consensus span overruns the axis (MEI_SPAN > axis length): "
          f"{span_overruns}")
    print("  (MEI_3P > axis length is impossible by construction of the axis, so only")
    print("   the span test carries information.)")
    print(f"  if 3' transduction ran at the stored literature band "
          f"{100 * lo:.1f}-{100 * hi:.1f}%, this call set would hold "
          f"{expected[0]:.0f}-{expected[1]:.0f} transduced insertions")
    print("\n  per-subfamily ceiling (top rows):")
    print(per_subfamily.head(10).to_string(index=False))
    print("\nARM B . junction evidence in the call fields")
    print(junction.to_string(index=False))
    print(f"\nARM C . read-level junction architecture ({loci_covered} loci with a "
          f"supporting_reads detail file)")
    if reads.empty:
        print("  no read carried both a full consensus match and a tail in one clip")
    else:
        print(f"  reads with a full consensus match plus a tail in one clip: {len(reads)}")
        print(f"  of which the caller's match explains the clip's non-tail bases: "
              f"{len(resolved)}")
        print("\n  resolved reads only (X is about the junction only here):")
        by_locus = (resolved.groupby(["locus", "family"])
                    .agg(reads=("X_bp", "size"), median_X=("X_bp", "median"),
                         min_X=("X_bp", "min"), max_X=("X_bp", "max"))
                    .reset_index())
        by_locus["n_X_ge_candidate"] = [
            int((resolved[(resolved.locus == l) & (resolved.family == f)].X_bp
                 >= args.candidate_x_bp).sum())
            for l, f in zip(by_locus.locus, by_locus.family)]
        print(by_locus.sort_values(["n_X_ge_candidate", "reads"],
                                   ascending=False).to_string(index=False))
        pooled = (resolved.groupby("family")
                  .agg(reads=("X_bp", "size"), median_X=("X_bp", "median"),
                       min_X=("X_bp", "min"), max_X=("X_bp", "max")).reset_index())
        pooled["n_X_ge_candidate"] = [
            int((resolved[resolved.family == f].X_bp >= args.candidate_x_bp).sum())
            for f in pooled.family]
        print("\n  pooled by family over resolved reads (read counts, not a rate):")
        print(pooled.to_string(index=False))
        print(f"\n  X >= {args.candidate_x_bp} bp is the candidate cut; negative X means the")
        print("  match plus tail already exceed the clip, i.e. alignment slack, which sets")
        print("  the resolution floor for this measurement.")
        dropped = reads[~reads.match_explains_clip]
        if not dropped.empty:
            print(f"  {len(dropped)} reads were dropped as unresolved: their consensus")
            print("  match covers too little of the clip for X to mean anything.")
    print(f"\nwrote {out}/axis_ceiling_by_subfamily.tsv, {out}/junction_support.tsv, "
          f"{out}/read_level_junction_x.tsv")
    print("NOTE: arm C covers one LINE-1 locus.  Nothing here is a rate.")


if __name__ == "__main__":
    main()
