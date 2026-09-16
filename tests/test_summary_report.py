"""End-to-end integration test for scientific summary report generation.

Generates synthetic BAM data, extracts evidence, exports annotated VCF,
and verifies the scientific summary report contains non-zero metrics.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import pandas as pd
import pytest

from retro_miner.candidate_loci import build_candidate_loci
from retro_miner.cleavage_motif import score_cleavage_motif
from retro_miner.evidence_extract import extract_split_and_discordant_evidence
from retro_miner.genotyper import GenotypeCall
from retro_miner.subfamily_voter import classify_mei_subfamily
from retro_miner.summary_report import generate_scientific_summary_report
from retro_miner.transduction_detector import detect_3prime_transduction_from_contig
from retro_miner.tsd_refiner import TSDResult
from retro_miner.vcf_export import write_mei_vcf
from scripts.generate_synthetic_bam import generate_synthetic_bam


def test_end_to_end_summary_report() -> None:
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        bam_path = tmp / "summary_test.bam"
        generate_synthetic_bam(
            output=bam_path,
            num_reads=500,
            seed=7,
            chromosomes=["chr1"],
            n_l1_insertions=2,
            n_alu_insertions=1,
            n_hard_clip_reads=1,
            n_interchrom_pairs=1,
        )

        evidence_dir = tmp / "evidence"
        evidence_dir.mkdir()
        extract_split_and_discordant_evidence(
            bam_path=bam_path,
            sample_name="disease",
            outdir=evidence_dir,
            regions=["chr1"],
        )

        summary_path = evidence_dir / "split_evidence.summary.tsv"
        summary_path.write_text(
            "sample\ttotal_reads_scanned\tpassing_reads\tsplit_evidence_rows\tdiscordant_evidence_rows\tinsert_size_threshold\tweak_only_discordant_filtered_rows\n"
            "disease\t500\t500\t10\t5\t200\t0\n"
            "control\t500\t500\t2\t1\t200\t0\n",
            encoding="utf-8",
        )

        (evidence_dir / "split_evidence.control.parquet").write_bytes(
            (evidence_dir / "split_evidence.disease.parquet").read_bytes()
        )
        (evidence_dir / "split_evidence.control.tsv").write_text(
            (evidence_dir / "split_evidence.disease.tsv").read_text(encoding="utf-8"),
            encoding="utf-8",
        )
        (evidence_dir / "discordant_evidence.control.parquet").write_bytes(
            (evidence_dir / "discordant_evidence.disease.parquet").read_bytes()
        )
        (evidence_dir / "discordant_evidence.control.tsv").write_text(
            (evidence_dir / "discordant_evidence.disease.tsv").read_text(encoding="utf-8"),
            encoding="utf-8",
        )

        candidate_tsv = build_candidate_loci(
            evidence_dir=evidence_dir,
            outdir=evidence_dir,
            window_size=200,
        )
        candidates = pd.read_csv(candidate_tsv, sep="\t")
        assert not candidates.empty, "Expected non-empty candidate loci table"

        vcf_records = []
        for _idx, row in enumerate(candidates.head(5).itertuples(index=False)):
            chrom = str(getattr(row, "chrom", "chr1"))
            pos = int(getattr(row, "window_start", 10000)) + 50

            sub_call = classify_mei_subfamily(["ACAGAG", "GACAGAG"], family_hint="L1")
            contig = "TTTTAAAA" + "A" * 12 + "GATTACACATGCAGCTAGCTAGCTAGCTA"
            tr = detect_3prime_transduction_from_contig(contig, mei_alignment_end=8, min_transduction_len=20)
            tsd = TSDResult(
                chrom=chrom,
                pos=pos,
                tsd_seq="ATTGCAG",
                tsd_length=7,
                confidence_score=0.85,
                poly_a_detected=True,
                entropy=1.8,
                locus_id="test",
            )
            vcf_records.append(
                {
                    "chrom": chrom,
                    "pos": pos,
                    "ref_base": "A",
                    "family": "L1",
                    "genotype_call": GenotypeCall(
                        genotype="0/1",
                        vaf=0.50,
                        log_likelihood_ratio=2.3134,
                        genotype_quality=72.0,
                        posterior_probs={"0/0": 0.0, "0/1": 1.0, "1/1": 0.0},
                    ),
                    "tsd_result": tsd,
                    "subfamily_call": sub_call,
                    "k_ref": 12,
                    "k_alt": 13,
                    "transduction_type": tr.transduction_type if tr.has_transduction else "NONE",
                    "transduction_seq": tr.transduction_seq,
                    "transduction_length": int(tr.transduction_length),
                    "tprt_motif_score": float(score_cleavage_motif(contig).pwm_score),
                }
            )

        vcf_path = tmp / "summary_test.vcf"
        write_mei_vcf(vcf_records, vcf_path, sample_name="summary_sample")
        assert vcf_path.exists()

        md_path = tmp / "summary_report.md"
        summary = generate_scientific_summary_report(vcf_path, md_path)

        assert summary.total_candidates == len(vcf_records)
        assert summary.mean_vaf > 0.0
        assert summary.tsds_detected_count > 0
        assert summary.poly_a_tails_detected_count > 0
        assert any(v >= 0.7 for v in [r.get("tprt_motif_score", 0.0) for r in vcf_records])

        md_text = md_path.read_text(encoding="utf-8")
        assert "# MEI Callset Scientific Summary Report" in md_text
        assert "## Genotype Distribution" in md_text
        assert "## Subfamily Breakdown" in md_text
        assert "## Mechanistic Annotations" in md_text
        assert "TPRT" in md_text
        assert "Transductions Detected" in md_text

        json_path = md_path.with_suffix(".json")
        assert json_path.exists()
        payload = __import__("json").loads(json_path.read_text(encoding="utf-8"))
        assert payload["total_candidates"] == len(vcf_records)
        assert payload["mean_vaf"] > 0.0


def test_empty_vcf_produces_zero_summary(tmp_path: Path) -> None:
    vcf = tmp_path / "empty.vcf"
    vcf.write_text("##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n", encoding="utf-8")
    md = tmp_path / "empty.md"
    summary = generate_scientific_summary_report(vcf, md)

    assert summary.total_candidates == 0
    assert summary.mean_vaf == 0.0
    assert summary.tsds_detected_count == 0
    assert summary.poly_a_tails_detected_count == 0
    assert summary.transductions_detected_count == 0
    assert summary.canonical_tprt_count == 0

    json_path = md.with_suffix(".json")
    assert json_path.exists()


def test_vcf_without_mechanistic_annotations_counts_only_subfamilies(tmp_path: Path) -> None:
    vcf = tmp_path / "minimal.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.3\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "chr1\t1000\t.\tA\t<INS:MEI:L1HS>\t60\tPASS\tSUBFAM=L1HS;MEI_LLR=3.5\tGT:GQ:VAF:AD\t0/1:60:0.50:10,12\n"
        "chr2\t2000\t.\tC\t<INS:MEI:ALU>\t55\tPASS\tSUBFAM=AluYa5;MEI_LLR=2.1\tGT:GQ:VAF:AD\t1/1:55:0.95:1,45\n",
        encoding="utf-8",
    )
    md = tmp_path / "minimal.md"
    summary = generate_scientific_summary_report(vcf, md)

    assert summary.total_candidates == 2
    assert summary.subfamily_breakdown == {"L1HS": 1, "AluYa5": 1}
    assert summary.tsds_detected_count == 0
    assert summary.poly_a_tails_detected_count == 0
    assert summary.transductions_detected_count == 0
    assert summary.canonical_tprt_count == 0


def test_vcf_with_malformed_vaf_skips_bad_records(tmp_path: Path) -> None:
    vcf = tmp_path / "bad_vaf.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.3\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "chr1\t1000\t.\tA\t<INS:MEI:L1HS>\t60\tPASS\tSUBFAM=L1HS;TPRT_MOTIF_SCORE=0.9\tGT:GQ:VAF:AD\t0/1:60:not_a_number:10,12\n"
        "chr2\t2000\t.\tC\t<INS:MEI:L1HS>\t55\tPASS\tSUBFAM=L1HS;TPRT_MOTIF_SCORE=0.8\tGT:GQ:VAF:AD\t1/1:55:0.95:1,45\n",
        encoding="utf-8",
    )
    md = tmp_path / "bad_vaf.md"
    summary = generate_scientific_summary_report(vcf, md)

    assert summary.total_candidates == 2
    assert summary.mean_vaf == pytest.approx(0.95)


def test_summary_report_json_roundtrip(tmp_path: Path) -> None:
    from retro_miner.summary_report import load_summary_json

    vcf = tmp_path / "roundtrip.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.3\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "chr1\t1000\t.\tA\t<INS:MEI:L1HS>\t60\tPASS\tSUBFAM=L1HS;TSD=ATTGCAG;TSDLEN=7;POLYA=1;TRANSDUCTION_TYPE=3_PRIME_PARTNERED;TRANSDUCTION_SEQ=GATTACA;TRANSDUCTION_LENGTH=7;TPRT_MOTIF_SCORE=0.85\tGT:GQ:VAF:AD\t0/1:60:0.50:10,12\n",
        encoding="utf-8",
    )
    md = tmp_path / "roundtrip.md"
    summary = generate_scientific_summary_report(vcf, md)
    loaded = load_summary_json(md.with_suffix(".json"))

    assert loaded == summary
