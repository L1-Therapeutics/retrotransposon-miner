"""One-pass extract must match the historical two-pass tables."""

from __future__ import annotations

import shutil
from pathlib import Path

import pandas as pd
import pytest

from retro_miner.evidence_extract import (
    extract_discordant_evidence,
    extract_split_and_discordant_evidence,
    extract_split_evidence,
)


def _have_pysam_and_samtools() -> bool:
    if shutil.which("samtools") is None:
        return False
    try:
        import pysam  # noqa: F401
    except ImportError:
        return False
    return True


def _write_pair(
    bam,
    *,
    qname: str,
    chrom: str,
    pos: int,
    mate_chrom: str,
    mate_pos: int,
    seq: str,
    cigar: str,
    is_read1: bool,
    is_reverse: bool,
    mate_is_reverse: bool,
    is_proper_pair: bool,
    tlen: int,
    mapq: int = 60,
    nm: int = 0,
) -> None:
    import pysam

    a = pysam.AlignedSegment()
    a.query_name = qname
    a.query_sequence = seq
    flag = 1 | (64 if is_read1 else 128)
    if is_proper_pair:
        flag |= 2
    if is_reverse:
        flag |= 16
    if mate_is_reverse:
        flag |= 32
    a.flag = flag
    a.reference_id = bam.get_tid(chrom)
    a.reference_start = pos
    a.next_reference_id = bam.get_tid(mate_chrom)
    a.next_reference_start = mate_pos
    a.template_length = tlen
    a.mapping_quality = mapq
    a.cigarstring = cigar
    a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
    a.set_tag("NM", nm)
    bam.write(a)


def _build_extract_bam(tmp_path: Path) -> Path:
    import pysam

    raw = tmp_path / "raw.bam"
    sorted_bam = tmp_path / "extract.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "SQ": [{"SN": "chr22", "LN": 50_000_000}, {"SN": "chr1", "LN": 50_000_000}],
    }
    clip = "A" * 25 + "G" * 25
    normal = "ACGT" * 15  # 60bp

    with pysam.AlignmentFile(str(raw), "wb", header=header) as bam:
        # Proper pairs spanning a range of insert sizes (quantile + floor).
        for i, tlen in enumerate((300, 350, 400, 450, 500, 550, 600, 800, 1200, 2500)):
            qn = f"proper_{i}"
            _write_pair(
                bam,
                qname=qn,
                chrom="chr22",
                pos=1000 + i * 200,
                mate_chrom="chr22",
                mate_pos=1000 + i * 200 + abs(tlen) - 60,
                seq=normal,
                cigar="60M",
                is_read1=True,
                is_reverse=False,
                mate_is_reverse=True,
                is_proper_pair=True,
                tlen=tlen,
            )
            _write_pair(
                bam,
                qname=qn,
                chrom="chr22",
                pos=1000 + i * 200 + abs(tlen) - 60,
                mate_chrom="chr22",
                mate_pos=1000 + i * 200,
                seq=normal,
                cigar="60M",
                is_read1=False,
                is_reverse=True,
                mate_is_reverse=False,
                is_proper_pair=True,
                tlen=-tlen,
            )

        # Soft-clipped split read (left 25S).
        _write_pair(
            bam,
            qname="split_left",
            chrom="chr22",
            pos=20_000,
            mate_chrom="chr22",
            mate_pos=20_300,
            seq=clip,
            cigar="25S25M",
            is_read1=True,
            is_reverse=False,
            mate_is_reverse=True,
            is_proper_pair=True,
            tlen=300,
        )
        _write_pair(
            bam,
            qname="split_left",
            chrom="chr22",
            pos=20_300,
            mate_chrom="chr22",
            mate_pos=20_000,
            seq=normal,
            cigar="60M",
            is_read1=False,
            is_reverse=True,
            mate_is_reverse=False,
            is_proper_pair=True,
            tlen=-300,
        )

        # Interchrom discordant pair.
        _write_pair(
            bam,
            qname="disc_inter",
            chrom="chr22",
            pos=30_000,
            mate_chrom="chr1",
            mate_pos=5_000,
            seq=normal,
            cigar="60M",
            is_read1=True,
            is_reverse=False,
            mate_is_reverse=False,
            is_proper_pair=False,
            tlen=0,
        )
        _write_pair(
            bam,
            qname="disc_inter",
            chrom="chr1",
            pos=5_000,
            mate_chrom="chr22",
            mate_pos=30_000,
            seq="T" * 60,
            cigar="60M",
            is_read1=False,
            is_reverse=False,
            mate_is_reverse=False,
            is_proper_pair=False,
            tlen=0,
        )

        # Same-strand only (weak; should be counted then dropped).
        _write_pair(
            bam,
            qname="weak_ss",
            chrom="chr22",
            pos=40_000,
            mate_chrom="chr22",
            mate_pos=40_200,
            seq=normal,
            cigar="60M",
            is_read1=True,
            is_reverse=False,
            mate_is_reverse=False,
            is_proper_pair=False,
            tlen=200,
        )

    pysam.sort("-o", str(sorted_bam), str(raw))
    pysam.index(str(sorted_bam))
    return sorted_bam


def _read_tables(outdir: Path, sample: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    split = pd.read_csv(outdir / f"split_evidence.{sample}.tsv", sep="\t")
    disc = pd.read_csv(outdir / f"discordant_evidence.{sample}.tsv", sep="\t")
    return split, disc


def _assert_tables_equal(left: pd.DataFrame, right: pd.DataFrame, *, label: str) -> None:
    left_s = left.sort_values(list(left.columns), kind="mergesort").reset_index(drop=True)
    right_s = right.sort_values(list(right.columns), kind="mergesort").reset_index(drop=True)
    pd.testing.assert_frame_equal(left_s, right_s, check_dtype=False, obj=label)


@pytest.mark.skipif(not _have_pysam_and_samtools(), reason="pysam and samtools required")
@pytest.mark.parametrize("fetch_mate_seq", [False, True])
def test_one_pass_matches_two_pass_tables(tmp_path: Path, fetch_mate_seq: bool) -> None:
    bam = _build_extract_bam(tmp_path)
    two = tmp_path / "two"
    one = tmp_path / "one"
    two.mkdir()
    one.mkdir()
    two_split = extract_split_evidence(
        bam_path=bam,
        sample_name="disease",
        outdir=two,
        regions="chr22",
        min_mapq=20,
        min_clip_len=20,
    )
    two_disc = extract_discordant_evidence(
        bam_path=bam,
        sample_name="disease",
        outdir=two,
        regions="chr22",
        min_mapq=20,
        insert_quantile=0.995,
        min_abs_tlen=1000,
        fetch_mate_seq=fetch_mate_seq,
        mate_bam_path=bam,
    )
    one_split, one_disc = extract_split_and_discordant_evidence(
        bam_path=bam,
        sample_name="disease",
        outdir=one,
        regions="chr22",
        min_mapq=20,
        min_mapq_discordant=20,
        min_clip_len=20,
        insert_quantile=0.995,
        min_abs_tlen=1000,
        fetch_mate_seq=fetch_mate_seq,
        mate_bam_path=bam,
    )

    two_split_df, two_disc_df = _read_tables(two, "disease")
    one_split_df, one_disc_df = _read_tables(one, "disease")
    _assert_tables_equal(two_split_df, one_split_df, label="split_evidence")
    _assert_tables_equal(two_disc_df, one_disc_df, label="discordant_evidence")

    assert one_split.total_reads_scanned == two_split.total_reads_scanned
    assert one_split.passing_reads == two_split.passing_reads
    assert one_split.split_evidence_rows == two_split.split_evidence_rows
    assert one_disc.passing_reads == two_disc.passing_reads
    assert one_disc.discordant_evidence_rows == two_disc.discordant_evidence_rows
    assert one_disc.insert_size_threshold == two_disc.insert_size_threshold
    assert one_disc.weak_only_discordant_filtered_rows == two_disc.weak_only_discordant_filtered_rows
    assert one_disc.mate_seq_fetched_rows == two_disc.mate_seq_fetched_rows
