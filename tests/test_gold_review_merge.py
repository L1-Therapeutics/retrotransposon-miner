"""Genome-wide gold-review merge uses the single-chromosome review sort."""

from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd
import pytest

from retro_miner.gold_review_merge import (
    incomplete_chromosomes,
    main,
    merge_gold_review_tables,
    write_merged_gold_review,
)

# Real HG03086 gold-review rows (ALU / LINE1 / SVA) for VCF + reference checks.
_VCF_FIXTURE_DIR = Path(__file__).resolve().parent / "fixtures" / "vcf"
_GOLD_SLICE = _VCF_FIXTURE_DIR / "hg03086_gold_slice.tsv"


def _row(
    chrom: str,
    mei_mapped: int,
    *,
    tier: str,
) -> dict[str, object]:
    support = f"SR_L=1,SR_R=1,DPE_L=0,DPE_R=0,MEI_MAPPED={mei_mapped},polyA_MAPPED=0,VNTR_MAPPED=0"
    return {
        "chrom": chrom,
        "window_start": "1",
        "window_end": "2",
        "consensus_insertion_breakpoint_pos": "10",
        "consensus_mei_family": "Alu",
        "known_mei_polymorphism": "False",
        "analysis_stage_tier": tier,
        "gold_stage_pass": "True" if tier == "gold" else "False",
        "silver_stage_pass": "True" if tier in {"gold", "silver"} else "False",
        "disease_supporting_reads": support,
        "control_supporting_reads": support,
        "coherence_score": "0.8",
        "insertion_model_score": "0.5",
    }


def test_merge_keeps_gold_and_ranks_by_support() -> None:
    chr1 = pd.DataFrame(
        [
            _row("chr1", mei_mapped=2, tier="gold"),
            _row("chr1", mei_mapped=99, tier="silver"),
            _row("chr1", mei_mapped=50, tier="bronze"),
        ]
    )
    chr22 = pd.DataFrame([_row("chr22", mei_mapped=40, tier="gold")])
    ranked = merge_gold_review_tables([chr1, chr22])
    assert ranked["chrom"].tolist() == ["chr22", "chr1"]
    assert ranked["analysis_stage_tier"].tolist() == ["gold", "gold"]
    assert ranked["gold_stage_pass"].tolist() == [1, 1]


def _finish(base: Path, chrom: str, *, done: bool) -> None:
    chrom_dir = base / chrom
    chrom_dir.mkdir(parents=True)
    (chrom_dir / "candidate_loci.mei.gold_review.tsv").write_text("chrom\twindow_start\n", encoding="utf-8")
    log_dir = base / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    line = f"stage=annotate-mei-support done region={chrom}\n" if done else "stage=extract\n"
    (log_dir / f"{chrom}.log").write_text(line, encoding="utf-8")


def test_incomplete_chromosome_blocks_aggregation(tmp_path: Path) -> None:
    _finish(tmp_path, "chr1", done=True)
    _finish(tmp_path, "chr2", done=False)
    assert incomplete_chromosomes(tmp_path, ["chr1", "chr2", "chr3"]) == ["chr2", "chr3"]
    out = tmp_path / "genome.tsv"
    with pytest.raises(SystemExit, match="chr2, chr3"):
        main(["--output", str(out), "--base-outdir", str(tmp_path), "chr1", "chr2", "chr3"])
    assert not out.exists()


def test_write_merged_gold_review_reads_per_chrom_tables(tmp_path: Path) -> None:
    chr1 = tmp_path / "chr1.tsv"
    chr22 = tmp_path / "chr22.tsv"
    pd.DataFrame(
        [
            _row("chr1", mei_mapped=2, tier="gold"),
            _row("chr1", mei_mapped=99, tier="silver"),
        ]
    ).to_csv(chr1, sep="\t", index=False)
    pd.DataFrame([_row("chr22", mei_mapped=40, tier="gold")]).to_csv(chr22, sep="\t", index=False)
    out = tmp_path / "candidate_loci.mei.gold_review.tsv"
    n_rows = write_merged_gold_review([chr1, chr22], out)
    assert n_rows == 2
    written = pd.read_csv(out, sep="\t")
    assert written["chrom"].tolist() == ["chr22", "chr1"]


def test_aggregated_gold_table_writes_vcf_from_pipeline_params(tmp_path: Path) -> None:
    """Aggregate a few real HG03086 gold rows and emit a ##reference VCF."""
    rows = list(csv.DictReader(_GOLD_SLICE.open(), delimiter="\t"))
    # Fixture order: Alu + L1 + SVA (high score) then a low-score L1; keep the first three.
    selected = rows[:3]
    assert [r["consensus_mei_family"] for r in selected] == ["ALU", "LINE1", "SVA"]
    fieldnames = list(rows[0].keys())
    inputs: list[Path] = []
    for row in selected:
        path = tmp_path / f"{row['chrom']}_{row['consensus_insertion_breakpoint_pos']}.tsv"
        with path.open("w", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerow(row)
        inputs.append(path)
    (tmp_path / "pipeline_params.env").write_text(
        "reference_build=hg38\nreference_fasta=/data/reference/hg38/Homo_sapiens_assembly38.fasta\n",
        encoding="utf-8",
    )
    out = tmp_path / "candidate_loci.mei.gold_review.tsv"
    assert main(["--output", str(out), *(str(p) for p in inputs)]) == 0
    text = out.with_suffix(".vcf").read_text(encoding="utf-8")
    assert "##reference=hg38\n" in text
    assert "##assembly=GRCh38\n" in text
    assert ",assembly=" not in text
    assert "L1TX-chr4-66240893-ALU" in text
    assert "L1TX-chr4-102422592-LINE1" in text
    assert "L1TX-chr10-3041566-SVA" in text
