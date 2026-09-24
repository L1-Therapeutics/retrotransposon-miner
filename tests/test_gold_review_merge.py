"""Genome-wide gold-review merge uses the single-chromosome review sort."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from retro_miner.gold_review_merge import merge_gold_review_tables, write_merged_gold_review


def _row(
    chrom: str,
    mei_mapped: int,
    *,
    gold: bool,
) -> dict[str, object]:
    support = f"SR_L=1,SR_R=1,DPE_L=0,DPE_R=0,MEI_MAPPED={mei_mapped},polyA_MAPPED=0,VNTR_MAPPED=0"
    return {
        "chrom": chrom,
        "window_start": "1",
        "window_end": "2",
        "consensus_insertion_breakpoint_pos": "10",
        "consensus_mei_family": "Alu",
        "known_mei_polymorphism": "False",
        "analysis_stage_tier": "gold" if gold else "silver",
        "gold_stage_pass": "True" if gold else "False",
        "silver_stage_pass": "True",
        "disease_supporting_reads": support,
        "control_supporting_reads": support,
        "coherence_score": "0.8",
        "insertion_model_score": "0.5",
    }


def test_merge_ranks_by_stage_then_support_across_chromosomes() -> None:
    chr1 = pd.DataFrame(
        [
            _row("chr1", mei_mapped=2, gold=True),
            _row("chr1", mei_mapped=99, gold=False),
        ]
    )
    chr22 = pd.DataFrame([_row("chr22", mei_mapped=40, gold=True)])
    ranked = merge_gold_review_tables([chr1, chr22])
    assert ranked["chrom"].tolist() == ["chr22", "chr1", "chr1"]
    assert ranked["analysis_stage_tier"].tolist() == ["gold", "gold", "silver"]
    assert ranked["gold_stage_pass"].tolist() == [1, 1, 0]


def test_write_merged_gold_review_reads_per_chrom_tables(tmp_path: Path) -> None:
    chr1 = tmp_path / "chr1.tsv"
    chr22 = tmp_path / "chr22.tsv"
    pd.DataFrame([_row("chr1", mei_mapped=2, gold=True)]).to_csv(chr1, sep="\t", index=False)
    pd.DataFrame([_row("chr22", mei_mapped=40, gold=True)]).to_csv(chr22, sep="\t", index=False)
    out = tmp_path / "candidate_loci.mei.gold_review.tsv"
    n_rows = write_merged_gold_review([chr1, chr22], out)
    assert n_rows == 2
    written = pd.read_csv(out, sep="\t")
    assert written["chrom"].tolist() == ["chr22", "chr1"]
