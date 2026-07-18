"""Pure polyA/polyT must never be reported as consensus_tsd_seq."""

from __future__ import annotations

import pandas as pd

from retro_miner.mei_support import (
    _clear_poly_at_artifact_tsd_fields,
    _poly_at_artifact_tsd_mask,
)


def test_poly_at_mask_catches_pure_homopolymers():
    seqs = pd.Series(
        [
            "AAAAAAAAAAAA",
            "TTTTTTTTTTTT",
            "AAGAAAACTCCT",
            "CCGCCTCGGCTTCCCAAAGTGCTGGGATTA",
            "",
            "ATA",  # pure AT but <4 bp — below TSD minimum
            # Long polyA/T with a short GC tip must still be cleared.
            "AAAAAAAAAAAAAAAAAAAAAAAGTTTTAA",
            "TTTTTTTTTTTTTTTTTTGAGATGGAGTCT",
        ]
    )
    mask = _poly_at_artifact_tsd_mask(seqs)
    assert bool(mask.iloc[0])
    assert bool(mask.iloc[1])
    assert not bool(mask.iloc[2])
    assert not bool(mask.iloc[3])
    assert not bool(mask.iloc[4])
    assert not bool(mask.iloc[5])
    assert bool(mask.iloc[6])
    assert bool(mask.iloc[7])


def test_clear_poly_at_applies_to_primary_tsd_disease_source():
    """Regression: filter used to require *_rescue in tsd_evidence_source."""
    out = pd.DataFrame(
        {
            "tsd_seq": ["TTTTTTTTTTTTTTTTTTTTTTTTTTT", "AAGAAAACTCCT", "AAAAAAAAAAAAA"],
            "tsd_len_estimate": [27, 12, 13],
            "tsd_left_breakpoint": [100, 200, 300],
            "tsd_right_breakpoint": [126, 211, 312],
            "tsd_detected": [True, True, True],
            "tsd_evidence_source": ["tsd_disease", "tsd_disease", "tsd_control"],
        }
    )
    mask = _clear_poly_at_artifact_tsd_fields(out)
    assert list(mask.astype(bool)) == [True, False, True]
    assert out.loc[0, "tsd_seq"] == ""
    assert int(out.loc[0, "tsd_len_estimate"]) == 0
    assert out.loc[1, "tsd_seq"] == "AAGAAAACTCCT"
    assert out.loc[2, "tsd_seq"] == ""
    assert "filtered_poly_at_only" in str(out.loc[0, "tsd_evidence_source"])
