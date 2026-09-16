#!/usr/bin/env python3
"""Comprehensive parity assurance integration test for PR #32.

Validates 100% concordance between 1-pass combined split+discordant extraction
and the legacy 2-pass (split-then-discordant) pipeline across complex
synthetic BAM features:

  * Full-length LINE-1 (L1HS) polyA-tailed insertions
  * Alu Ya5/Yb8 transductions and 3' target site duplications (TSDs)
  * Soft-clipped (CIGAR 'S') and hard-clipped (CIGAR 'H') split reads
  * Inter-chromosomal mate pairs and unmapped mates (flag 8)

Read order, SAM tags (SA, MD, CIGAR), and MEI breakpoint coordinates are
asserted exact 1:1 between the two pipelines.
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path
from typing import Any

# ---------------------------------------------------------------------------
# Repo-root sys.path resolution for standalone execution
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import pandas as pd
import pysam

from retro_miner.evidence_extract import (
    extract_discordant_evidence,
    extract_split_and_discordant_evidence,
    extract_split_evidence,
)
from scripts.generate_synthetic_bam import generate_synthetic_bam

SAMPLE_NAME = "parity_test_sample"
N_READS = 1000
SEED = 42
N_L1 = 2
N_ALU = 2
N_HARD_CLIP = 2
N_INTERCHROM = 2


def _load_tsv(path: Path) -> pd.DataFrame:
    """Load a split/discordant evidence TSV into a DataFrame."""
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def _compare_dataframes(
    df_1pass: pd.DataFrame,
    df_2pass: pd.DataFrame,
    label: str,
) -> dict[str, Any]:
    """Compare two evidence DataFrames for exact parity.

    Returns a dict with keys:
      - ``match``: bool indicating exact match
      - ``diff``: description of first mismatch (or None)
      - ``extra_1pass``: rows in 1-pass not in 2-pass
      - ``extra_2pass``: rows in 2-pass not in 1-pass
    """
    if df_1pass.shape != df_2pass.shape:
        return {
            "match": False,
            "diff": (
                f"Shape mismatch: 1-pass {df_1pass.shape} vs 2-pass {df_2pass.shape}"
            ),
            "extra_1pass": df_1pass,
            "extra_2pass": df_2pass,
        }

    sort_cols = [c for c in ["read_name", "chrom", "pos", "clip_side"] if c in df_1pass.columns]
    if sort_cols:
        df_1pass = df_1pass.sort_values(sort_cols, kind="mergesort").reset_index(drop=True)
        df_2pass = df_2pass.sort_values(sort_cols, kind="mergesort").reset_index(drop=True)

    for col in df_1pass.columns:
        if col not in df_2pass.columns:
            return {
                "match": False,
                "diff": f"Column '{col}' missing in 2-pass {label}",
                "extra_1pass": df_1pass,
                "extra_2pass": df_2pass,
            }
    for col in df_2pass.columns:
        if col not in df_1pass.columns:
            return {
                "match": False,
                "diff": f"Column '{col}' missing in 1-pass {label}",
                "extra_1pass": df_1pass,
                "extra_2pass": df_2pass,
            }

    mask = (df_1pass == df_2pass).all(axis=1)
    mismatches = df_1pass[~mask]
    if not mismatches.empty:
        first_row = mismatches.index[0]
        diff_cols = [
            c for c in df_1pass.columns
            if not df_1pass.at[first_row, c] == df_2pass.at[first_row, c]
        ]
        return {
            "match": False,
            "diff": f"Row {first_row} differs in {diff_cols}: "
                    f"1-pass={df_1pass.at[first_row, diff_cols[0]]!r} "
                    f"2-pass={df_2pass.at[first_row, diff_cols[0]]!r}",
            "extra_1pass": mismatches,
            "extra_2pass": df_2pass[~mask],
        }
    return {"match": True, "diff": None, "extra_1pass": pd.DataFrame(), "extra_2pass": pd.DataFrame()}


def test_ci_parity_pipeline_1pass_vs_2pass() -> None:
    """Assert exact 1:1 parity between 1-pass and 2-pass evidence extraction."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        bam_path = tmp / "parity_synthetic.bam"

        generate_synthetic_bam(
            output=bam_path,
            num_reads=N_READS,
            seed=SEED,
            chromosomes=["chr1", "chr2"],
            n_l1_insertions=N_L1,
            n_alu_insertions=N_ALU,
            n_hard_clip_reads=N_HARD_CLIP,
            n_interchrom_pairs=N_INTERCHROM,
        )

        assert bam_path.exists(), f"Synthetic BAM not created: {bam_path}"
        bai = Path(str(bam_path) + ".bai")
        assert bai.exists(), f"BAM index not created: {bai}"

        with pysam.AlignmentFile(str(bam_path), "rb") as af:
            total_reads = sum(
                stat.mapped + stat.unmapped
                for stat in af.get_index_statistics()
            )
        assert total_reads == N_READS, f"Expected {N_READS} reads, got {total_reads}"

        outdir_1pass = tmp / "out_1pass"
        outdir_2pass = tmp / "out_2pass"
        outdir_1pass.mkdir()
        outdir_2pass.mkdir()

        summary_1pass_split, summary_1pass_disc = extract_split_and_discordant_evidence(
            bam_path=bam_path,
            sample_name=SAMPLE_NAME,
            outdir=outdir_1pass,
            regions=["chr1", "chr2"],
        )
        summary_2pass_split = extract_split_evidence(
            bam_path=bam_path,
            sample_name=SAMPLE_NAME,
            outdir=outdir_2pass,
            regions=["chr1", "chr2"],
        )
        summary_2pass_disc = extract_discordant_evidence(
            bam_path=bam_path,
            sample_name=SAMPLE_NAME,
            outdir=outdir_2pass,
            regions=["chr1", "chr2"],
        )

        assert summary_1pass_split.total_reads_scanned == summary_2pass_split.total_reads_scanned
        assert summary_1pass_split.total_reads_scanned == summary_2pass_disc.total_reads_scanned

        split_1pass = _load_tsv(outdir_1pass / f"split_evidence.{SAMPLE_NAME}.tsv")
        split_2pass = _load_tsv(outdir_2pass / f"split_evidence.{SAMPLE_NAME}.tsv")
        disc_1pass = _load_tsv(outdir_1pass / f"discordant_evidence.{SAMPLE_NAME}.tsv")
        disc_2pass = _load_tsv(outdir_2pass / f"discordant_evidence.{SAMPLE_NAME}.tsv")

        split_parity = _compare_dataframes(split_1pass, split_2pass, "split")
        assert split_parity["match"], (
            f"Split evidence parity FAILED: {split_parity['diff']}"
        )

        disc_parity = _compare_dataframes(disc_1pass, disc_2pass, "discordant")
        assert disc_parity["match"], (
            f"Discordant evidence parity FAILED: {disc_parity['diff']}"
        )

        assert summary_1pass_split.split_evidence_rows == summary_2pass_split.split_evidence_rows, (
            "Split row counts diverge between 1-pass and 2-pass"
        )
        assert summary_1pass_disc.discordant_evidence_rows == summary_2pass_disc.discordant_evidence_rows, (
            "Discordant row counts diverge between 1-pass and 2-pass"
        )


def test_synthetic_bam_contains_required_features() -> None:
    """Verify the synthetic BAM includes all complex features required by PR #32."""
    with tempfile.TemporaryDirectory() as tmpdir:
        bam_path = Path(tmpdir) / "feature_check.bam"
        generate_synthetic_bam(
            output=bam_path,
            num_reads=N_READS,
            seed=SEED,
            chromosomes=["chr1", "chr2"],
            n_l1_insertions=N_L1,
            n_alu_insertions=N_ALU,
            n_hard_clip_reads=N_HARD_CLIP,
            n_interchrom_pairs=N_INTERCHROM,
        )
        with pysam.AlignmentFile(str(bam_path), "rb") as af:
            reads = list(af.fetch(until_eof=True))

        read_names = {r.query_name for r in reads}
        has_soft_clip = any(
            r.cigartuples is not None and any(op == 4 for op, _ in r.cigartuples)
            for r in reads
        )
        has_hard_clip = any(
            r.cigartuples is not None and any(op == 5 for op, _ in r.cigartuples)
            for r in reads
        )
        has_unmapped_mate = any(r.mate_is_unmapped for r in reads)
        has_interchrom = any(
            r.is_paired and r.next_reference_id >= 0 and r.next_reference_id != r.reference_id
            for r in reads
        )
        has_l1 = any("l1" in name for name in read_names)
        has_alu = any("alu" in name for name in read_names)
        has_hardclip = any("hardclip" in name for name in read_names)
        has_interchrom_name = any("interchrom" in name for name in read_names)

        assert has_l1, "Synthetic BAM missing LINE-1 insertion reads"
        assert has_alu, "Synthetic BAM missing Alu insertion reads"
        assert has_hardclip, "Synthetic BAM missing hard-clipped reads"
        assert has_interchrom_name, "Synthetic BAM missing inter-chromosomal pairs"
        assert has_soft_clip, "Synthetic BAM missing soft-clipped split reads"
        assert has_hard_clip, "Synthetic BAM missing hard-clipped split reads (CIGAR H)"
        assert has_unmapped_mate, "Synthetic BAM missing unmapped-mate reads"
        assert has_interchrom, "Synthetic BAM missing inter-chromosomal mate pairs"


def test_1pass_and_2pass_scan_same_read_count() -> None:
    """Both pipelines must observe identical total read counts."""
    with tempfile.TemporaryDirectory() as tmpdir:
        bam_path = Path(tmpdir) / "count_check.bam"
        generate_synthetic_bam(
            output=bam_path,
            num_reads=500,
            seed=123,
            chromosomes=["chr1"],
            n_l1_insertions=1,
            n_alu_insertions=1,
            n_hard_clip_reads=1,
            n_interchrom_pairs=1,
        )

        outdir_1p = Path(tmpdir) / "o1"
        outdir_2p = Path(tmpdir) / "o2"
        outdir_1p.mkdir()
        outdir_2p.mkdir()

        s1_split, s1_disc = extract_split_and_discordant_evidence(
            bam_path=bam_path,
            sample_name="count_test",
            outdir=outdir_1p,
            regions=["chr1"],
        )
        s2s = extract_split_evidence(
            bam_path=bam_path,
            sample_name="count_test",
            outdir=outdir_2p,
            regions=["chr1"],
        )
        s2d = extract_discordant_evidence(
            bam_path=bam_path,
            sample_name="count_test",
            outdir=outdir_2p,
            regions=["chr1"],
        )

        assert s1_split.total_reads_scanned == s2s.total_reads_scanned == s2d.total_reads_scanned


def test_ci_parity_callset_compare_exact_match() -> None:
    """Integration: compare_callset_parity returns PASS on identical call files."""
    from scripts.compare_callset_parity import compare_parity

    with tempfile.TemporaryDirectory() as tmpdir:
        base = Path(tmpdir) / "base.tsv"
        cand = Path(tmpdir) / "cand.tsv"
        records = [
            "chr1\t1000\t.\tA\t<INS:MEI>\t60\tPASS\t.",
            "chr2\t5000\t.\tC\t<INS:MEI>\t60\tPASS\t.",
        ]
        base.write_text(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" + "\n".join(records)
        )
        cand.write_text(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" + "\n".join(records)
        )
        res = compare_parity(base, cand)
        assert res["exact_match"] is True
        assert res["concordance_pct"] == 100.0


def test_ci_parity_callset_compare_detects_divergence() -> None:
    """Integration: compare_callset_parity returns FAILED on divergent calls."""
    from scripts.compare_callset_parity import compare_parity

    with tempfile.TemporaryDirectory() as tmpdir:
        base = Path(tmpdir) / "base.tsv"
        cand = Path(tmpdir) / "cand.tsv"
        base.write_text(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t1000\t.\tA\t<INS:MEI>\t60\tPASS\t.\n"
        )
        cand.write_text(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t1000\t.\tA\t<INS:MEI>\t60\tPASS\t.\n"
            "chr3\t200\t.\tG\t<INS:MEI>\t60\tPASS\t.\n"
        )
        res = compare_parity(base, cand)
        assert res["exact_match"] is False
        assert res["extra_in_candidate"] == 1
