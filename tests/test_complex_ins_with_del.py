"""Same-chrom del-bridge cluster + breakpoint depth drop → COMPLEX_INS_WITH_DEL."""

from __future__ import annotations

import pandas as pd

from retro_miner.mei_support import (
    _aggregate_same_chrom_deletion_dpe_metrics,
    _apply_complex_ins_with_del,
    _deletion_depth_supports_del,
    _deletion_flank_intervals,
    _drop_deletion_cluster_reads,
    _refresh_polya_rescue_excluding_del_cluster,
    _same_chrom_deletion_cluster_member_reads,
)


def _dpe_row(**kwargs):
    base = {
        "chrom": "chr22",
        "window_start": 10000,
        "window_end": 10100,
        "pos": 10040,
        "mate_chrom": "chr22",
        "mate_pos": 15000,
        "read_name": "r",
        "is_reverse": False,
        "mate_is_reverse": True,
        "discordant_reasons": "large_insert",
    }
    base.update(kwargs)
    return base


class TestDeletionFlankIntervals:
    def test_mates_to_the_right_use_left_intact_flank(self):
        intact, interior = _deletion_flank_intervals(10000, 16000, flank_bp=500)
        assert intact == (9500, 9999)
        assert interior == (10001, 10500)

    def test_mates_to_the_left_use_right_intact_flank(self):
        intact, interior = _deletion_flank_intervals(16000, 10000, flank_bp=500)
        assert intact == (16001, 16500)
        assert interior == (15500, 15999)


class TestDeletionDepthDrop:
    def test_half_depth_is_a_drop(self):
        assert _deletion_depth_supports_del(30.0, 15.0)
        assert _deletion_depth_supports_del(30.0, 0.0)
        assert _deletion_depth_supports_del(30.0, 19.5)  # 0.65

    def test_no_drop_when_interior_stays_high(self):
        assert not _deletion_depth_supports_del(30.0, 20.0)  # 0.667

    def test_uninterpretable_low_intact_depth(self):
        assert not _deletion_depth_supports_del(4.0, 0.0)


class TestDeletionClusterFraction:
    def test_minority_cluster_still_reports_fraction(self):
        rows = []
        for i in range(3):
            rows.append(_dpe_row(read_name=f"del{i}", mate_pos=15000 + i * 10, pos=10040 + i))
        for i in range(7):
            rows.append(
                _dpe_row(
                    read_name=f"other{i}",
                    mate_chrom="chr1",
                    mate_pos=2000 + i,
                    discordant_reasons="interchrom",
                )
            )
        out = _aggregate_same_chrom_deletion_dpe_metrics(pd.DataFrame(rows), "disease")
        assert len(out) == 1
        assert int(out.iloc[0]["disease_deletion_cluster_reads"]) == 3
        assert abs(float(out.iloc[0]["disease_deletion_cluster_fraction"]) - 0.3) < 1e-9


class TestDeletionClusterMeiMappedDrop:
    def test_cluster_members_are_excluded_from_mei_support(self):
        rows = []
        for i in range(3):
            rows.append(_dpe_row(read_name=f"del{i}", mate_pos=15000 + i * 10, pos=10040 + i))
        for i in range(7):
            rows.append(
                _dpe_row(
                    read_name=f"mei{i}",
                    mate_chrom="chr1",
                    mate_pos=2000 + i,
                    discordant_reasons="interchrom",
                )
            )
        dpe = pd.DataFrame(rows)
        members = _same_chrom_deletion_cluster_member_reads(dpe)
        assert set(members["read_name"]) == {"del0", "del1", "del2"}
        kept = _drop_deletion_cluster_reads(dpe, members)
        assert set(kept["read_name"]) == {f"mei{i}" for i in range(7)}

    def test_sub_threshold_cluster_keeps_mei_reads(self):
        rows = [_dpe_row(read_name="del0", mate_pos=15000, pos=10040)]
        for i in range(9):
            rows.append(
                _dpe_row(
                    read_name=f"other{i}",
                    mate_chrom="chr1",
                    mate_pos=2000 + i,
                    discordant_reasons="interchrom",
                )
            )
        dpe = pd.DataFrame(rows)
        members = _same_chrom_deletion_cluster_member_reads(dpe)
        assert members.empty
        kept = _drop_deletion_cluster_reads(dpe, members)
        assert len(kept) == 10

    def test_del_cluster_mei_hits_do_not_unlock_polya_rescue(self):
        rows = []
        for i in range(4):
            rows.append(
                _dpe_row(
                    read_name=f"del{i}",
                    mate_pos=15000 + i * 10,
                    pos=10040 + i,
                    mei_hit=True,
                    mate_mei_hit=True,
                    family="ALU",
                    target="AluY#SINE/Alu",
                    polya_rescue=False,
                )
            )
        rows.append(
            _dpe_row(
                read_name="poly",
                mate_pos=20000,
                pos=10050,
                discordant_reasons="large_insert",
                mei_hit=False,
                mate_mei_hit=False,
                mate_seq="A" * 40,
                polya_rescue=True,
                family="ALU",
                target="ALU_polyA_rescue#SINE/Alu",
                mei_hit_source="polya_rescue",
            )
        )
        dpe = pd.DataFrame(rows)
        members = _same_chrom_deletion_cluster_member_reads(dpe)
        out = _refresh_polya_rescue_excluding_del_cluster(dpe, members)
        poly = out.loc[out["read_name"].eq("poly")].iloc[0]
        assert not bool(poly["polya_rescue"])


class TestComplexInsWithDelLabel:
    def test_labels_only_when_depth_drops(self):
        df = pd.DataFrame(
            {
                "insertion_event_class": ["SIMPLE_MEI", "SIMPLE_MEI"],
                "disease_deletion_depth_drop": [True, False],
                "control_deletion_depth_drop": [False, False],
            }
        )
        out = _apply_complex_ins_with_del(df)
        assert bool(out.loc[0, "complex_ins_with_del"])
        assert out.loc[0, "insertion_event_class"] == "COMPLEX_INS_WITH_DEL"
        assert not bool(out.loc[1, "complex_ins_with_del"])
        assert out.loc[1, "insertion_event_class"] == "SIMPLE_MEI"
