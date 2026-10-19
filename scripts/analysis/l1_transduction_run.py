"""Run the per-call transduction detector on real data, and check the reference reader.

Two things happen here, and the first is a prerequisite for trusting the second:

1. **Reference check.**  ``retro_miner.transduction_detect.TwoBitReference`` reads
   hg38 from ``ref/hg38.2bit`` because no 2bit library is installed in this
   environment.  A bit-order or offset mistake in that reader would silently
   change every donor flank the detector searches, so the reader is tested
   against evidence it did not produce: the aligned bases of real BAM reads.  If
   the reader is right the mismatch rate is sequencing-plus-mapping error
   (roughly 0.2-3%); if it is wrong the comparison is noise (roughly 65-75%).

   Strand convention is measured, not assumed.  SAM stores SEQ already oriented
   to the reference (a reverse-strand read is stored reverse-complemented
   relative to the original read), so an aligned base is compared as stored.
   The complemented tally is kept as a deliberate control: had it come out near
   zero instead, that would mean these files store original-read orientation.
   A reader error sits near 0.65-0.75 in *both* columns, so neither can hide
   one.

2. **Detection.**  A source-element index is built from the RepeatMasker track
   (full-length L1 copies that carry their 3' end, full-length SVAs that carry
   their 5' end) with donor flanks taken from the reference, then the detector
   classifies every junction clip recovered from the fixture BAMs: partnered
   (a donor is named), orphan (transduction-like, absent from the reference),
   overhang (the segment is the reference next to the insertion) or
   not_transduction (poly(A) extension, low complexity, too short).

The index is genome-wide by default and should stay that way.  An earlier version
of this script restricted it to three assemblies to keep the run small, and that
restriction *manufactured a donor*: a clip matched one chr18 L1HS flank over a
28 bp span, and with the competing copies absent from the index the detector
named it.  Those 28 bp sit in 46 genome-wide 3'-competent L1 flanks, so the
singular donor was an artifact of scope, not a link.  The real cost of doing it
properly is measured, not assumed: 9,360 sources, 15.4 Mbp of flank, about 140 s
and ~2 GB peak RSS.  ``--chroms`` still exists for a deliberate subset, and any
run that uses it sees a donor set truncated by assembly order.

Scope, stated because it bounds the result: the fixture loci are one sample at 22
loci, four of them LINE-1, against 2,186 cohort calls.  Nothing here is a rate.

Usage:
    python l1_transduction_run.py [--loci rtm-longread-cohort/tests/fixtures/loci]
                                  [--chroms chr18,chr19,chr22]   # optional subset
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import resource
import sys
import time
from collections.abc import Iterator
from pathlib import Path

import pandas as pd
import pysam

sys.path.insert(0, str(Path(__file__).resolve().parent / "retrotransposon-miner" / "src"))

from retro_miner.transduction_detect import (  # noqa: E402
    DEFAULT_MIN_MATCHED_BP,
    DEFAULT_MIN_SEGMENT_BP,
    TwoBitReference,
    annotate_transduction,
    build_source_index,
    candidate_source_elements,
    rmsk_rows,
)

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def _peak_rss_gb() -> float:
    """Peak resident set size in GB.

    ``ru_maxrss`` is bytes on macOS and kilobytes on Linux, so a single divisor
    labels one platform's number wrong by three orders of magnitude (a 38 MB
    peak printed as "38.81 GB", which is what prompted this helper).
    """
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return rss / (1024 ** 3) if sys.platform == "darwin" else rss / (1024 ** 2)


# --- 1. reference reader check -------------------------------------------------


#: Forward reads are compared one way; a reverse read's stored bases are
#: compared as stored (the reference-oriented convention) and again
#: complemented, as a control on the strand convention itself.
CHECK_TALLIES = (
    "fwd_match", "fwd_mismatch",
    "rev_stored_match", "rev_stored_mismatch",
    "rev_complemented_match", "rev_complemented_mismatch",
)


def _rate(counts: dict[str, int], prefix: str) -> float | None:
    matched = counts.get(f"{prefix}_match", 0)
    total = matched + counts.get(f"{prefix}_mismatch", 0)
    return round(1 - matched / total, 5) if total else None


def check_reference(path: str, loci: str, max_reads: int) -> pd.DataFrame:
    """Compare 2bit-fetched reference bases against real BAM alignments.

    Returns per-locus mismatch rates for forward reads, for reverse reads as
    stored (the convention SAM defines, and the operative number), and for
    reverse reads complemented (the control).  The control is the point: it is
    the number a *wrong* strand convention would drive to ~0, so seeing it sit
    at ~1 while the as-stored column sits at ~0.01 is positive evidence that the
    files follow SAM's convention and that the reader decoded them correctly.
    """
    ref = TwoBitReference(path)
    rows = []
    for locus_dir in sorted(glob.glob(os.path.join(loci, "*"))):
        bams = sorted(glob.glob(os.path.join(locus_dir, "*.bam")))
        if not bams:
            continue
        counts = {key: 0 for key in CHECK_TALLIES}
        with pysam.AlignmentFile(bams[0], "rb") as bam:
            for read in bam.fetch(until_eof=True):
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                if read.mapping_quality < 20 or read.cigartuples is None:
                    continue
                seq = read.query_sequence or ""
                for query_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
                    if query_pos is None or ref_pos is None or query_pos >= len(seq):
                        continue
                    base = seq[query_pos].upper()
                    fetch_start = ref_pos + 1
                    reference = ref.fetch(read.reference_name, fetch_start, fetch_start)
                    if not reference or base not in "ACGT" or reference == "N":
                        continue
                    if read.is_reverse:
                        counts["rev_stored_match" if base == reference
                               else "rev_stored_mismatch"] += 1
                        counts["rev_complemented_match" if _revcomp(base) == reference
                               else "rev_complemented_mismatch"] += 1
                    else:
                        counts["fwd_match" if base == reference else "fwd_mismatch"] += 1
                if sum(counts.values()) >= max_reads * 150:
                    break
        rows.append({
            "locus": os.path.basename(locus_dir),
            "ref_bases_compared": sum(counts.values()),
            "forward_reads": counts["fwd_match"] + counts["fwd_mismatch"],
            "reverse_reads": (counts["rev_stored_match"] + counts["rev_stored_mismatch"]),
            "forward_mismatch_rate": _rate(counts, "fwd"),
            "reverse_mismatch_rate_as_stored": _rate(counts, "rev_stored"),
            "reverse_mismatch_rate_complemented_control": _rate(counts, "rev_complemented"),
        })
    return pd.DataFrame(rows)


# --- 2. junction clips from the fixture BAMs -----------------------------------


def junction_clips(loci: str, ref: TwoBitReference, *, min_clip: int = 20,
                   breakpoint_window: int = 50, local_bp: int = 300) -> pd.DataFrame:
    rows = []
    for locus_dir in sorted(glob.glob(os.path.join(loci, "*"))):
        bams = sorted(glob.glob(os.path.join(locus_dir, "*.bam")))
        gold_path = os.path.join(locus_dir, "gold_locus.tsv")
        manifest_path = os.path.join(locus_dir, "manifest.json")
        if not bams or not os.path.exists(manifest_path):
            continue
        manifest = json.loads(Path(manifest_path).read_text())
        breakpoint = manifest.get("expected_breakpoint")
        family = ""
        subfamily = ""
        if os.path.exists(gold_path):
            gold = pd.read_csv(gold_path, sep="\t", low_memory=False)
            family = str(gold.get("consensus_mei_family", pd.Series([""]))[0]).strip()
            subfamily = str(gold.get("consensus_mei_subfamily", pd.Series([""]))[0]).strip()
        locus_id = os.path.basename(locus_dir)
        seen: set[tuple[str, str]] = set()
        with pysam.AlignmentFile(bams[0], "rb") as bam:
            for read in bam.fetch(until_eof=True):
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                if read.is_duplicate or read.is_qcfail or read.mapping_quality < 20:
                    continue
                if read.cigartuples is None or not read.cigartuples:
                    continue
                seq = (read.query_sequence or "").upper()
                clip_len, side, tip = 0, "", None
                first_op, first_len = read.cigartuples[0]
                last_op, last_len = read.cigartuples[-1]
                if first_op == 4 and first_len >= min_clip:
                    clip_len, side, tip = int(first_len), "L", int(read.reference_start) + 1
                if last_op == 4 and last_len > clip_len:
                    clip_len, side = int(last_len), "R"
                    tip = int(read.reference_end) if read.reference_end is not None else None
                if not side or tip is None:
                    continue
                if breakpoint is not None and abs(tip - int(breakpoint)) > breakpoint_window:
                    continue
                key = (read.query_name, side)
                if key in seen:
                    continue
                seen.add(key)
                clip = seq[:clip_len] if side == "L" else seq[-clip_len:]
                if side == "L":
                    local = _revcomp(ref.fetch(read.reference_name,
                                              max(1, tip - local_bp + 1), tip))
                else:
                    local = ref.fetch(read.reference_name, tip, tip + local_bp - 1)
                rows.append({
                    "call_id": f"{locus_id}:{read.query_name}:{side}",
                    "locus": locus_id,
                    "family": family,
                    "subfamily": subfamily,
                    "clip_side": side,
                    "clip_bp": clip_len,
                    "clip_tip": tip,
                    "junction_seq": clip,
                    "local_reference_seq": local,
                })
    return pd.DataFrame(rows)


def _rows_for(rmsk: str, assemblies: set[str]) -> Iterator[dict[str, str]]:
    """Annotation rows for the assemblies a run indexes (all of them if empty).

    Genome-wide is the default because a subset index is a different claim: a
    donor named from a filtered index has been chosen from an arbitrary slice of
    the copies that carry the segment, and section 4 exists to measure that.
    """
    if not assemblies:
        yield from rmsk_rows(rmsk)
        return
    for row in rmsk_rows(rmsk):
        if str(row.get("genoName", "")) in assemblies:
            yield row


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rmsk", default="nested_analysis/data/rmsk.txt.gz")
    ap.add_argument("--reference", default="ref/hg38.2bit")
    ap.add_argument("--loci", default="rtm-longread-cohort/tests/fixtures/loci")
    ap.add_argument("--chroms", default="",
                    help="optional comma-separated subset of assemblies to index; "
                         "home of the 28 bp artifact, so empty (all) is the default")
    ap.add_argument("--max-elements", type=int, default=None,
                    help="cap on indexed candidate sources; a cap truncates the "
                         "donor set in annotation order, so leave it unset for a "
                         "result anyone could call a donor search")
    ap.add_argument("--min-matched-bp", type=int, default=DEFAULT_MIN_MATCHED_BP,
                    help="shortest seed-supported span a donor may be named from")
    ap.add_argument("--reuse-reference-check", action="store_true",
                    help="reuse a cached reference_reader_check.tsv instead of "
                         "re-reading every aligned base; the 2bit reader is not "
                         "changed by this run, so the cached verdict still applies")
    ap.add_argument("--out-dir", default="l1_transduction_out")
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    ref = TwoBitReference(args.reference)

    check_path = out / "reference_reader_check.tsv"
    print("1. REFERENCE CHECK . 2bit reader vs real BAM alignments")
    if args.reuse_reference_check and check_path.exists():
        check = pd.read_csv(check_path, sep="\t")
        print(f"  reused {check_path} | loci {len(check)} | bases compared "
              f"{int(check.ref_bases_compared.sum()):,}")
    else:
        check = check_reference(args.reference, args.loci, max_reads=400)
        check.to_csv(check_path, sep="\t", index=False)
        pooled_bases = int(check.ref_bases_compared.sum())
        print(f"  loci checked {len(check)} | reference bases compared {pooled_bases:,}"
              f" (forward {int(check.forward_reads.sum()):,},"
              f" reverse {int(check.reverse_reads.sum()):,})")
        for label, column in (
            ("forward reads", "forward_mismatch_rate"),
            ("reverse reads, as stored (SAM convention)", "reverse_mismatch_rate_as_stored"),
            ("reverse reads, complemented (control)",
             "reverse_mismatch_rate_complemented_control"),
        ):
            rates = check[column].dropna()
            if rates.empty:
                print(f"  {label}: no comparable bases")
                continue
            print(f"  {label}: mismatch rate median {rates.median():.4f} "
                  f"min {rates.min():.4f} max {rates.max():.4f}")
        print("  (a decoding error sits near 0.65-0.75 in every column; ~0 in the")
        print("   control column instead would mean a non-SAM strand convention)")

    print("\n2. SOURCE INDEX . full-length L1/SVA copies from the RepeatMasker track")
    chroms = {c.strip() for c in args.chroms.split(",") if c.strip()}
    print(f"  assemblies indexed: {','.join(sorted(chroms)) if chroms else 'all'}")

    t0 = time.time()
    candidates = candidate_source_elements(_rows_for(args.rmsk, chroms), ref,
                                          families=("LINE1", "SVA"),
                                          max_elements=args.max_elements)
    index = build_source_index(candidates)
    elements = candidates
    by_family: dict[str, int] = {}
    by_dir: dict[str, int] = {}
    for element in elements:
        by_family[element.family] = by_family.get(element.family, 0) + 1
        by_dir[element.flank_dir] = by_dir.get(element.flank_dir, 0) + 1
    print(f"  candidate sources indexed: {len(elements)} {by_family} | flank directions {by_dir}")
    print(f"  extraction + index build: {time.time() - t0:.0f}s | peak RSS "
          f"{_peak_rss_gb():.2f} GB (measured, because it decides whether the "
          f"index may be genome-wide)")

    print("\n3. DETECTION . junction clips from the fixture BAMs")
    clips = junction_clips(args.loci, ref)
    clips.to_csv(out / "junction_clips.tsv", sep="\t", index=False)
    print(f"  junction clips recovered: {len(clips)} at {clips.locus.nunique()} loci")
    if not clips.empty:
        print("  LINE1 clips: "
              f"{int((clips.family == 'LINE1').sum())} at "
              f"{clips[clips.family == 'LINE1'].locus.nunique()} loci")

    results = annotate_transduction(clips.to_dict("records"), index,
                                    min_segment_bp=DEFAULT_MIN_SEGMENT_BP,
                                    min_matched_bp=args.min_matched_bp)
    frame = pd.DataFrame([r.to_dict() for r in results])
    if not frame.empty:
        # `family` is already a column here; assign rather than concat, because a
        # shared name would leave two columns called `family` and break selection.
        frame.insert(0, "locus", clips["locus"].to_numpy())
        frame.insert(1, "clip_side", clips["clip_side"].to_numpy())
        frame.insert(2, "clip_bp", clips["clip_bp"].to_numpy())
    frame.to_csv(out / "transduction_calls.tsv", sep="\t", index=False)
    print("\n  outcome counts:")
    if frame.empty:
        print("    (no clips to classify)")
    else:
        for outcome, count in frame["outcome"].value_counts().items():
            print(f"    {outcome}: {count}")
        print("  reason counts (why a clip was not called a partnered transduction):")
        for reason, count in frame["reason"].value_counts().items():
            print(f"    {reason}: {count}")
        for fam in sorted(frame["family"].unique()):
            sub = frame[frame["family"] == fam]
            counts = ", ".join(f"{outcome} {count}" for outcome, count
                              in sub["outcome"].value_counts().items())
            print(f"  {fam or '(unlabelled family)'}: {len(sub)} clips -> {counts}")
        print("  named donors:")
        named = frame[frame.donor_source_id.astype(str).str.len() > 0]
        if named.empty:
            print("    none -- no junction clip in this fixture matched an indexed donor flank")
        else:
            for _, row in named.iterrows():
                print(f"    {row['call_id']} -> {row['donor_source_id']} "
                      f"({row['donor_matched_start0']}-{row['donor_matched_end0']}, "
                      f"{row['donor_orientation']}, {row['transduced_bp']} bp)")
        refused = frame[frame["reason"] == "matched_span_below_specificity_floor"]
        if not refused.empty:
            print("  matched a flank, but too short a span to name a donor:")
            for _, row in refused.iterrows():
                print(f"    {row['call_id']} -> span {row['transduced_bp']} bp,"
                      f" {row['n_donor_candidates']} carrier(s)")
    print("\n4. INDEX SCOPE CONTROL . the same clips, a truncated index and the full one")
    # The first version of this script indexed three assemblies and capped the
    # source list, and that is the difference between naming a donor and showing
    # that forty loci carry the same short piece.  Measured, not argued: the same
    # clips, the same detector, only the index scope and the span floor move.
    scope = {"chr18", "chr19", "chr22"}
    scoped_elements = candidate_source_elements(_rows_for(args.rmsk, scope), ref,
                                                families=("LINE1", "SVA"), max_elements=4000)
    scoped_index = build_source_index(scoped_elements)
    for label, scoped in (
        (f"3 assemblies, cap 4000 ({len(scoped_elements)} sources)", scoped_index),
        (f"all assemblies ({len(elements)} sources)", index),
    ):
        for floor in (0, DEFAULT_MIN_MATCHED_BP):
            alt = annotate_transduction(clips.to_dict("records"), scoped,
                                        min_segment_bp=DEFAULT_MIN_SEGMENT_BP,
                                        min_matched_bp=floor)
            alt_frame = pd.DataFrame([call.to_dict() for call in alt])
            named = alt_frame[alt_frame.donor_source_id.astype(str).str.len() > 0]
            spans = sorted({int(bp) for bp in
                            pd.to_numeric(named["transduced_bp"], errors="coerce").dropna()})
            donors = sorted(set(named.donor_source_id.astype(str)))
            reasons = ", ".join(f"{reason} {count}" for reason, count
                                in alt_frame["reason"].value_counts().items()
                                if reason in ("donor_match", "donor_tie",
                                              "matched_span_below_specificity_floor"))
            print(f"  {label}, min_matched_bp={floor}:")
            print(f"    named {donors or ['none']} | spans {spans}"
                  f" | {reasons or 'no clip matched a flank'}")

    print(f"\nwrote {out}/reference_reader_check.tsv, {out}/junction_clips.tsv, "
          f"{out}/transduction_calls.tsv")


if __name__ == "__main__":
    main()
