#!/usr/bin/env python3
"""Generate synthetic paired-end BAM files for retrotransposon-miner testing.

Produces a deterministic (seeded), coordinate-sorted BAM whose reads mix:

  - normal paired-end reads,
  - LINE-1 insertion reads (soft-clipped split anchors + discordant pairs),
  - Alu retrotransposon insertion reads, and
  - alignment edge cases
      * secondary alignments (flag 256),
      * supplementary alignments (flag 2048),
      * unmapped mates (flag 8), and
      * zero mapping quality (MAPQ 0).

The companion ``*.bai`` index is written automatically via ``pysam.index()``.

Example:
    python3 scripts/generate_synthetic_bam.py \
      --output /tmp/syn_chr22.bam \
      --num-reads 1000 \
      --seed 42 \
      --chromosomes chr1,chr22 \
      --n-l1-insertions 2 \
      --n-alu-insertions 3
"""

from __future__ import annotations

import argparse
import random
from pathlib import Path

import pysam

READ_LENGTH = 100

HG38_CHROM_LENGTHS = {
    "chr1": 248956422,
    "chr2": 242193529,
    "chr3": 198295559,
    "chr4": 190214555,
    "chr5": 181538259,
    "chr6": 170805979,
    "chr7": 159345973,
    "chr8": 145138636,
    "chr9": 138394717,
    "chr10": 133797422,
    "chr11": 135086622,
    "chr12": 133275309,
    "chr13": 114364328,
    "chr14": 107043718,
    "chr15": 101991189,
    "chr16": 90338345,
    "chr17": 83257441,
    "chr18": 80373285,
    "chr19": 58617616,
    "chr20": 64444167,
    "chr21": 46709983,
    "chr22": 50818468,
    "chrX": 156040895,
    "chrY": 57227415,
    "chrM": 16569,
}


def default_chromosomes() -> list[str]:
    """Return the full canonical hg38 chromosome set (chr1-22, X, Y, M)."""
    return ["chr%d" % i for i in range(1, 23)] + ["chrX", "chrY", "chrM"]


def validate_chromosomes(chromosomes: list[str]) -> list[str]:
    """Validate hg38-style chromosome names and reject duplicates/unknown contigs."""
    seen: set[str] = set()
    for chrom in chromosomes:
        if chrom not in HG38_CHROM_LENGTHS:
            raise ValueError("not an hg38 contig")
        if chrom in seen:
            raise ValueError("duplicate chromosome in target list: %s" % chrom)
        seen.add(chrom)
    return list(chromosomes)


def _normalize_chromosomes(
    chromosomes: list[str] | str | None,
) -> list[str]:
    """Accept None, a comma-separated string, or a list/tuple."""
    if chromosomes is None:
        return default_chromosomes()
    if isinstance(chromosomes, str):
        names = [c.strip() for c in chromosomes.split(",") if c.strip()]
        if not names:
            raise ValueError("chromosomes string is empty")
        return validate_chromosomes(names)
    if isinstance(chromosomes, (list, tuple)):
        return validate_chromosomes([str(c) for c in chromosomes])
    raise ValueError(
        "chromosomes must be None, a comma-separated string, or a list/tuple"
    )


def _random_sequence(rng: random.Random, length: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(length))


def _make_segment(
    name: str,
    flag: int,
    tid: int,
    pos: int,
    mapq: int,
    cigar: list[tuple[int, int]],
    sequence: str,
    mate_tid: int | None = None,
    mate_pos: int | None = None,
) -> pysam.AlignedSegment:
    """Build a single AlignedSegment; deterministic given inputs."""
    seg = pysam.AlignedSegment()
    seg.query_name = name
    seg.query_sequence = sequence
    seg.flag = flag
    seg.reference_id = tid
    seg.reference_start = pos
    seg.mapping_quality = mapq
    seg.cigar = cigar
    seg.mrnm = mate_tid if mate_tid is not None else -1
    seg.mpos = mate_pos if mate_pos is not None else -1
    return seg


def _pair_base(
    rng: random.Random,
    chrom: str,
    tid: int,
    chrom_len: int,
    start: int,
    insert_size: int,
    tag: str,
    index: int,
    mapq: int = 60,
) -> tuple[pysam.AlignedSegment, pysam.AlignedSegment]:
    """Build a concordant 99/147 pair with the given insert size."""
    seq1 = _random_sequence(rng, READ_LENGTH)
    seq2 = _random_sequence(rng, READ_LENGTH)
    r1 = _make_segment(
        "%s_%06d_1" % (tag, index),
        99,
        tid,
        start,
        mapq,
        [(0, READ_LENGTH)],
        seq1,
        mate_tid=tid,
        mate_pos=start + insert_size,
    )
    r2 = _make_segment(
        "%s_%06d_2" % (tag, index),
        147,
        tid,
        start + insert_size,
        mapq,
        [(0, READ_LENGTH)],
        seq2,
        mate_tid=tid,
        mate_pos=start,
    )
    return r1, r2


def _insertion_span(element_len: int, chrom_len: int) -> int:
    """Clamp insertion element span so it fits inside short contigs (e.g. chrM)."""
    return max(READ_LENGTH, min(element_len, chrom_len // 2))


def _breakpoint(
    rng: random.Random, chrom_len: int, span: int
) -> tuple[int, int]:
    """Pick a valid breakpoint/mate window inside the chromosome."""
    lo = 3 * READ_LENGTH
    hi = chrom_len - span - 2 * READ_LENGTH - 500
    if hi < lo:
        lo = max(0, (lo + hi) // 2)
    pos = rng.randint(lo, hi) if hi > lo else lo
    return pos, lo


def _hard_clip_split_read(
    rng: random.Random,
    chrom: str,
    tid: int,
    chrom_len: int,
    breakpoint: int,
    tag: str,
    index: int,
) -> tuple[pysam.AlignedSegment, pysam.AlignedSegment]:
    """Build a hard-clipped split-read pair (CIGAR 'H').

    Hard clips (H) consume query bases that are NOT present in the query
    sequence field.  The query sequence therefore only covers the aligned
    (M) portion: 30H + 40M + 30H  ->  query length = 40.
    """
    seq1 = _random_sequence(rng, 40)
    seq2 = _random_sequence(rng, 60)
    r1 = _make_segment(
        "%s_hardclip_%06d_1" % (tag, index),
        99,
        tid,
        breakpoint - 40,
        60,
        [(5, 30), (0, 40), (5, 30)],  # 30H 40M 30H
        seq1,
        mate_tid=tid,
        mate_pos=breakpoint + 100,
    )
    r2 = _make_segment(
        "%s_hardclip_%06d_2" % (tag, index),
        147,
        tid,
        breakpoint + 100,
        60,
        [(5, 20), (0, 60), (5, 20)],  # 20H 60M 20H
        seq2,
        mate_tid=tid,
        mate_pos=breakpoint - 40,
    )
    return r1, r2


def build_read_plan(
    rng: random.Random,
    num_reads: int,
    chromosomes: list[str],
    n_l1_insertions: int,
    n_alu_insertions: int,
    ref_id: dict[str, int],
    ref_len: dict[str, int],
    sv_chrom: str,
    n_hard_clip_reads: int = 0,
    n_interchrom_pairs: int = 0,
) -> list[pysam.AlignedSegment]:
    """Assemble exactly ``num_reads`` alignment records."""
    sv_tid = ref_id[sv_chrom]
    sv_len = ref_len[sv_chrom]

    fixed_edge = 7
    fixed_insertions = 4 * (n_l1_insertions + n_alu_insertions)
    fixed_hard = 2 * n_hard_clip_reads
    fixed_interchrom = 2 * n_interchrom_pairs
    fixed = fixed_edge + fixed_insertions + fixed_hard + fixed_interchrom
    if num_reads < fixed:
        raise ValueError(
            "num_reads=%d is below minimum of %d required to include all edge cases "
            "(%d L1 + %d Alu + %d hard-clip + %d interchrom pairs)"
            % (
                num_reads,
                fixed,
                n_l1_insertions,
                n_alu_insertions,
                n_hard_clip_reads,
                n_interchrom_pairs,
            )
        )

    reads: list[pysam.AlignedSegment] = []

    # --- Secondary alignment (flag 256) ---------------------------------
    pos = rng.randint(3 * READ_LENGTH, sv_len - READ_LENGTH - 100)
    reads.append(
        _make_segment(
            "secondary_000001",
            256,
            sv_tid,
            pos,
            10,
            [(0, READ_LENGTH)],
            _random_sequence(rng, READ_LENGTH),
        )
    )

    # --- Supplementary alignment (flag 2048), chimeric read pair ---------
    seq = _random_sequence(rng, 2 * READ_LENGTH)
    pos = rng.randint(3 * READ_LENGTH, sv_len - 2 * 2 * READ_LENGTH - 100)
    reads.append(
        _make_segment(
            "suppl_000001",
            0,
            sv_tid,
            pos,
            60,
            [(0, 40), (4, 60)],
            seq[:READ_LENGTH],
        )
    )
    reads.append(
        _make_segment(
            "suppl_000001",
            2048,
            sv_tid,
            pos + 2 * READ_LENGTH,
            60,
            [(4, 60), (0, 40)],
            seq[READ_LENGTH:],
        )
    )

    # --- Unmapped mates (flag 8): both reads mapped, mate reported unmapped
    pos = rng.randint(3 * READ_LENGTH, sv_len - 2 * READ_LENGTH - 100)
    reads.append(
        _make_segment(
            "mate_unmapped_000001",
            73,
            sv_tid,
            pos,
            60,
            [(0, READ_LENGTH)],
            _random_sequence(rng, READ_LENGTH),
        )
    )
    reads.append(
        _make_segment(
            "mate_unmapped_000001",
            133,
            sv_tid,
            pos + READ_LENGTH,
            60,
            [(0, READ_LENGTH)],
            _random_sequence(rng, READ_LENGTH),
        )
    )

    # --- Zero mapping quality (MAPQ 0) ----------------------------------
    pos = rng.randint(3 * READ_LENGTH, sv_len - 3 * READ_LENGTH - 100)
    r1, r2 = _pair_base(
        rng, sv_chrom, sv_tid, sv_len, pos, READ_LENGTH, "mapq0", 1, mapq=0
    )
    reads.extend([r1, r2])

    # --- MEI insertions with SAM tags ------------------------------------
    count = 1
    for _element_len, _tag, _count in (
        (6000, "l1", n_l1_insertions),
        (300, "alu", n_alu_insertions),
    ):
        for _ in range(_count):
            span = _insertion_span(_element_len, sv_len)
            bp, _lo = _breakpoint(rng, sv_len, span)

            clipped_seq = _random_sequence(rng, READ_LENGTH)
            sr1 = _make_segment(
                "syn_%s_insertion_%06d" % (_tag, count),
                99,
                sv_tid,
                bp - 40,
                60,
                [(0, 50), (4, 50)],
                clipped_seq,
                mate_tid=sv_tid,
                mate_pos=bp + span,
            )
            sr1.set_tag("ME", "L1HS")
            sr1.set_tag("SA", "chr1,123456,+,60M,60,0;")
            reads.append(sr1)

            sr2 = _make_segment(
                "syn_%s_insertion_%06d" % (_tag, count),
                147,
                sv_tid,
                bp + span,
                60,
                [(0, READ_LENGTH)],
                _random_sequence(rng, READ_LENGTH),
                mate_tid=sv_tid,
                mate_pos=bp - 40,
            )
            reads.append(sr2)

            d1 = _make_segment(
                "syn_%s_insertion_%06d" % (_tag, count),
                97,
                sv_tid,
                bp - READ_LENGTH,
                60,
                [(0, READ_LENGTH)],
                _random_sequence(rng, READ_LENGTH),
                mate_tid=sv_tid,
                mate_pos=bp + span + READ_LENGTH,
            )
            reads.append(d1)

            d2 = _make_segment(
                "syn_%s_insertion_%06d" % (_tag, count),
                145,
                sv_tid,
                bp + span + READ_LENGTH,
                60,
                [(0, READ_LENGTH)],
                _random_sequence(rng, READ_LENGTH),
                mate_tid=sv_tid,
                mate_pos=bp - READ_LENGTH,
            )
            reads.append(d2)
            count += 1

    # --- Hard-clipped split reads (CIGAR 'H') ----------------------------
    for i in range(1, n_hard_clip_reads + 1):
        bp, _lo = _breakpoint(rng, sv_len, 200)
        hr1, hr2 = _hard_clip_split_read(rng, sv_chrom, sv_tid, sv_len, bp, "hardclip", i)
        reads.extend([hr1, hr2])

    # --- Inter-chromosomal mate pairs -------------------------------------
    if n_interchrom_pairs > 0 and len(chromosomes) >= 2:
        chrom_a = chromosomes[0]
        chrom_b = chromosomes[1]
        tid_a = ref_id[chrom_a]
        tid_b = ref_id[chrom_b]
        len_a = ref_len[chrom_a]
        len_b = ref_len[chrom_b]
        for i in range(1, n_interchrom_pairs + 1):
            pos_a = rng.randint(3 * READ_LENGTH, len_a - 2 * READ_LENGTH - 100)
            pos_b = rng.randint(3 * READ_LENGTH, len_b - 2 * READ_LENGTH - 100)
            ir1 = _make_segment(
                "interchrom_%06d_1" % i,
                99,
                tid_a,
                pos_a,
                60,
                [(0, READ_LENGTH)],
                _random_sequence(rng, READ_LENGTH),
                mate_tid=tid_b,
                mate_pos=pos_b,
            )
            ir2 = _make_segment(
                "interchrom_%06d_2" % i,
                147,
                tid_b,
                pos_b,
                60,
                [(0, READ_LENGTH)],
                _random_sequence(rng, READ_LENGTH),
                mate_tid=tid_a,
                mate_pos=pos_a,
            )
            reads.extend([ir1, ir2])

    # --- Normal paired-end reads fill the rest of the budget --------------
    n_normal = num_reads - fixed
    n_pairs = n_normal // 2
    leftover = n_normal % 2
    for i in range(1, n_pairs + 1):
        chrom = rng.choice(chromosomes)
        tid = ref_id[chrom]
        chrom_len = ref_len[chrom]
        insert_size = rng.randint(200, 600)
        start = rng.randint(3 * READ_LENGTH, chrom_len - insert_size - READ_LENGTH - 100)
        r1, r2 = _pair_base(rng, chrom, tid, chrom_len, start, insert_size, "normal", i)
        reads.extend([r1, r2])

    if leftover:
        chrom = rng.choice(chromosomes)
        tid = ref_id[chrom]
        chrom_len = ref_len[chrom]
        pos = rng.randint(3 * READ_LENGTH, chrom_len - READ_LENGTH - 100)
        reads.append(
            _make_segment(
                "single_end_000001",
                0,
                tid,
                pos,
                60,
                [(0, READ_LENGTH)],
                _random_sequence(rng, READ_LENGTH),
            )
        )

    return reads


def generate_synthetic_bam(
    output: Path | str,
    num_reads: int = 1000,
    seed: int = 42,
    chromosomes: list[str] | None = None,
    n_l1_insertions: int = 1,
    n_alu_insertions: int = 1,
    n_hard_clip_reads: int = 0,
    n_interchrom_pairs: int = 0,
) -> tuple[Path, Path]:
    """Write a synthetic, indexed BAM and return ``(bam_path, bai_path)``."""
    if num_reads < 5:
        raise ValueError("below minimum required (5)")

    chroms = _normalize_chromosomes(chromosomes)
    out = Path(output)
    out.parent.mkdir(parents=True, exist_ok=True)

    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": c, "LN": HG38_CHROM_LENGTHS[c]} for c in chroms],
    }
    ref_id = {c: i for i, c in enumerate(chroms)}
    ref_len = {c: HG38_CHROM_LENGTHS[c] for c in chroms}
    sv_chrom = max(chroms, key=lambda c: ref_len[c])

    rng = random.Random(seed)
    reads = build_read_plan(
        rng=rng,
        num_reads=num_reads,
        chromosomes=chroms,
        n_l1_insertions=n_l1_insertions,
        n_alu_insertions=n_alu_insertions,
        ref_id=ref_id,
        ref_len=ref_len,
        sv_chrom=sv_chrom,
        n_hard_clip_reads=n_hard_clip_reads,
        n_interchrom_pairs=n_interchrom_pairs,
    )
    reads.sort(key=lambda r: (r.reference_id, r.reference_start))

    with pysam.AlignmentFile(str(out), "wb", header=header) as outf:
        for seg in reads:
            outf.write(seg)

    pysam.index(str(out))
    bai = Path(str(out) + ".bai")
    return out, bai


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Generate a synthetic, indexed BAM with L1/Alu insertion reads "
        "and alignment edge cases for retrotransposon-miner testing."
    )
    p.add_argument("--output", type=Path, required=True, help="Output BAM path (e.g. /tmp/syn.bam)")
    p.add_argument("--num-reads", type=int, default=1000, help="Total number of reads to write")
    p.add_argument("--seed", type=int, default=42, help="Random seed for reproducible output")
    p.add_argument(
        "--chromosomes",
        default=None,
        help="Comma-separated hg38 contigs (default: all chr1-22,chrX,chrY,chrM)",
    )
    p.add_argument("--n-l1-insertions", type=int, default=1, help="Number of LINE-1 insertion loci to model")
    p.add_argument("--n-alu-insertions", type=int, default=1, help="Number of Alu insertion loci to model")
    p.add_argument("--n-hard-clip-reads", type=int, default=0, help="Number of hard-clipped split-read pairs")
    p.add_argument("--n-interchrom-pairs", type=int, default=0, help="Number of inter-chromosomal mate pairs")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.num_reads < 5:
        raise SystemExit("--num-reads must be at least 5")
    if args.n_l1_insertions < 0 or args.n_alu_insertions < 0:
        raise SystemExit("--n-l1-insertions and --n-alu-insertions must be non-negative")
    if args.n_hard_clip_reads < 0 or args.n_interchrom_pairs < 0:
        raise SystemExit("--n-hard-clip-reads and --n-interchrom-pairs must be non-negative")
    chroms = (
        [c.strip() for c in args.chromosomes.split(",") if c.strip()]
        if args.chromosomes
        else default_chromosomes()
    )
    bam, bai = generate_synthetic_bam(
        output=args.output,
        num_reads=args.num_reads,
        seed=args.seed,
        chromosomes=chroms,
        n_l1_insertions=args.n_l1_insertions,
        n_alu_insertions=args.n_alu_insertions,
        n_hard_clip_reads=args.n_hard_clip_reads,
        n_interchrom_pairs=args.n_interchrom_pairs,
    )
    print(
        "wrote %d reads -> %s\nindex           -> %s" % (args.num_reads, bam, bai)
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
