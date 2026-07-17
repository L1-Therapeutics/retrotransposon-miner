"""DPE MEI remap must use soft-clipped bases, not ref-matched read bodies."""

from __future__ import annotations

import pandas as pd

from retro_miner.evidence_extract import _soft_clip_query_seq
from retro_miner.mei_support import (
    _DPE_MEI_REMAP_MIN_CLIP_BP,
    _discordant_anchor_mei_query_seq,
    _discordant_mate_mei_query_seq,
)


class TestSoftClipQuerySeq:
    def test_left_and_right_clips(self):
        # 100bp "ref" + 50bp clip on right; left clip on a different read
        seq = ("R" * 100) + ("I" * 50)
        assert _soft_clip_query_seq(seq, "R", 50) == "I" * 50
        assert _soft_clip_query_seq(seq, "L", 50) == "R" * 50
        assert _soft_clip_query_seq(seq, "R", 20) == "I" * 20

    def test_too_short_or_missing(self):
        assert _soft_clip_query_seq("ACGT", "R", 0) == ""
        assert _soft_clip_query_seq("ACGT", "R", 10) == ""


class TestDiscordantMeiQuerySelection:
    def test_min_clip_threshold_is_20(self):
        assert _DPE_MEI_REMAP_MIN_CLIP_BP == 20

    def test_anchor_uses_only_soft_clip_not_full_read(self):
        """User example: 100bp ref + 50bp clip → remap the 50bp clip only."""
        row = pd.Series(
            {
                "read_seq": ("R" * 100) + ("I" * 50),
                "soft_clip_side": "R",
                "soft_clip_len": 50,
                "soft_clip_seq": "I" * 50,
            }
        )
        q = _discordant_anchor_mei_query_seq(row)
        assert q == "I" * 50
        assert "R" not in q

    def test_anchor_without_clip_is_empty(self):
        row = pd.Series(
            {
                "read_seq": "A" * 150,
                "soft_clip_side": "",
                "soft_clip_len": 0,
                "soft_clip_seq": "",
            }
        )
        assert _discordant_anchor_mei_query_seq(row) == ""

    def test_anchor_derives_clip_from_side_len_when_seq_missing(self):
        row = pd.Series(
            {
                "read_seq": ("R" * 100) + ("I" * 50),
                "soft_clip_side": "R",
                "soft_clip_len": 50,
                "soft_clip_seq": "",
            }
        )
        assert _discordant_anchor_mei_query_seq(row) == "I" * 50

    def test_anchor_rejects_short_clip(self):
        row = pd.Series(
            {
                "read_seq": ("R" * 140) + ("I" * 10),
                "soft_clip_side": "R",
                "soft_clip_len": 10,
                "soft_clip_seq": "I" * 10,
            }
        )
        assert _discordant_anchor_mei_query_seq(row) == ""

    def test_mate_clip_only_when_soft_clipped(self):
        """Opposite-junction mate: 50bp clip + 100bp ref → remap clip only."""
        row = pd.Series(
            {
                "mate_seq": ("I" * 50) + ("R" * 100),
                "mate_soft_clip_side": "L",
                "mate_soft_clip_len": 50,
                "mate_soft_clip_seq": "I" * 50,
            }
        )
        q = _discordant_mate_mei_query_seq(row)
        assert q == "I" * 50
        assert q.count("R") == 0

    def test_unclipped_mate_may_use_full_sequence(self):
        row = pd.Series(
            {
                "mate_seq": "M" * 120,
                "mate_soft_clip_side": "",
                "mate_soft_clip_len": 0,
                "mate_soft_clip_seq": "",
            }
        )
        assert _discordant_mate_mei_query_seq(row) == "M" * 120

    def test_short_mate_clip_does_not_fall_back_to_full_body(self):
        row = pd.Series(
            {
                "mate_seq": ("I" * 10) + ("R" * 140),
                "mate_soft_clip_side": "L",
                "mate_soft_clip_len": 10,
                "mate_soft_clip_seq": "I" * 10,
            }
        )
        assert _discordant_mate_mei_query_seq(row) == ""


class TestMergeFetchedMateSequences:
    def test_map_apply_is_vectorized_and_preserves_existing_seq(self):
        from retro_miner.mei_support import _merge_fetched_mate_sequences

        out = pd.DataFrame(
            {
                "read_name": ["a", "b", "c"],
                "mate_seq": ["KEEP", "", ""],
                "mate_ref_start": [0, 0, 0],
                "mate_ref_end": [0, 0, 0],
                "mate_soft_clip_side": ["", "", ""],
                "mate_soft_clip_len": [0, 0, 0],
                "mate_soft_clip_seq": ["", "", ""],
            }
        )
        fetched = {
            "a": ("AAAA", 10, 20, "R", 5, "AAAAA"),
            "b": ("BBBB", 30, 40, "L", 4, "BBBB"),
        }
        merged = _merge_fetched_mate_sequences(out, fetched)
        assert merged.loc[0, "mate_seq"] == "KEEP"  # existing kept
        assert merged.loc[0, "mate_soft_clip_side"] == "R"  # soft-clip refreshed
        assert merged.loc[1, "mate_seq"] == "BBBB"
        assert int(merged.loc[1, "mate_ref_start"]) == 30
        assert merged.loc[2, "mate_seq"] == ""  # not in fetched
