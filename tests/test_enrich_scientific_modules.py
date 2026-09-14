"""Unit tests for _enrich_candidates_with_scientific_modules in mei_support.py."""

from __future__ import annotations

import pandas as pd
import pytest

from retro_miner.mei_support import _enrich_candidates_with_scientific_modules


def _candidate_df(**rows: object) -> pd.DataFrame:
    if not rows:
        return pd.DataFrame(columns=["chrom", "window_start", "window_end"])
    return pd.DataFrame(rows)


def test_empty_candidate_returns_empty() -> None:
    candidate = pd.DataFrame(columns=["chrom", "window_start", "window_end"])
    result = _enrich_candidates_with_scientific_modules(
        candidate,
        split_disease=pd.DataFrame(),
        split_control=pd.DataFrame(),
        discordant_disease=pd.DataFrame(),
        discordant_control=pd.DataFrame(),
        asm_df=pd.DataFrame(),
    )
    assert result.empty


def test_empty_evidence_returns_candidate_with_defaults() -> None:
    candidate = _candidate_df(
        chrom=["chr1"], window_start=[1000], window_end=[1200]
    )
    result = _enrich_candidates_with_scientific_modules(
        candidate,
        split_disease=pd.DataFrame(),
        split_control=pd.DataFrame(),
        discordant_disease=pd.DataFrame(),
        discordant_control=pd.DataFrame(),
        asm_df=pd.DataFrame(),
    )
    assert len(result) == 1
    assert result.iloc[0]["subfamily"] == "UNKNOWN"
    assert result.iloc[0]["mei_llr"] == 0.0
    assert result.iloc[0]["transduction_type"] == "NONE"
    assert result.iloc[0]["transduction_length"] == 0
    assert result.iloc[0]["tprt_motif_score"] == 0.0


def test_locus_with_no_matching_sequences_gets_defaults() -> None:
    candidate = _candidate_df(
        chrom=["chr1"], window_start=[1000], window_end=[1200]
    )
    split_disease = pd.DataFrame(
        {
            "chrom": ["chr1"],
            "window_start": [1000],
            "window_end": [1200],
            "soft_clip_seq": ["AC"],
        }
    )
    result = _enrich_candidates_with_scientific_modules(
        candidate,
        split_disease=split_disease,
        split_control=pd.DataFrame(),
        discordant_disease=pd.DataFrame(),
        discordant_control=pd.DataFrame(),
        asm_df=pd.DataFrame(),
    )
    assert len(result) == 1
    assert result.iloc[0]["subfamily"] == "L1HS"
    assert result.iloc[0]["transduction_type"] == "NONE"


def test_locus_with_short_sequences_skips_assembly() -> None:
    candidate = _candidate_df(
        chrom=["chr1"], window_start=[1000], window_end=[1200]
    )
    split_disease = pd.DataFrame(
        {
            "chrom": ["chr1"],
            "window_start": [1000],
            "window_end": [1200],
            "soft_clip_seq": ["AC"],
        }
    )
    result = _enrich_candidates_with_scientific_modules(
        candidate,
        split_disease=split_disease,
        split_control=pd.DataFrame(),
        discordant_disease=pd.DataFrame(),
        discordant_control=pd.DataFrame(),
        asm_df=pd.DataFrame(),
    )
    assert len(result) == 1
    assert result.iloc[0]["transduction_type"] == "NONE"
    assert result.iloc[0]["tprt_motif_score"] == 0.0


def test_locus_with_sequences_classifies_subfamily() -> None:
    candidate = _candidate_df(
        chrom=["chr1"], window_start=[1000], window_end=[1200]
    )
    split_disease = pd.DataFrame(
        {
            "chrom": ["chr1"],
            "window_start": [1000],
            "window_end": [1200],
            "soft_clip_seq": ["ACAGAG"],
        }
    )
    result = _enrich_candidates_with_scientific_modules(
        candidate,
        split_disease=split_disease,
        split_control=pd.DataFrame(),
        discordant_disease=pd.DataFrame(),
        discordant_control=pd.DataFrame(),
        asm_df=pd.DataFrame(),
    )
    assert len(result) == 1
    assert result.iloc[0]["subfamily"] == "L1HS"
    assert float(result.iloc[0]["mei_llr"]) > 0.0
