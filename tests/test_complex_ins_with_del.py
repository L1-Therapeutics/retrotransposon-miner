"""Same-chrom del-bridge cluster + breakpoint depth drop → COMPLEX_INS_WITH_DEL."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pysam

from retro_miner.mei_support import (
    _aggregate_same_chrom_deletion_dpe_metrics,
    _apply_complex_ins_with_del,
    _best_ungapped_identity,
    _deletion_depth_supports_del,
    _deletion_flank_intervals,
    _drop_deletion_cluster_reads,
    _gapped_local_match,
    _refresh_polya_rescue_excluding_del_cluster,
    _revcomp,
    _same_chrom_deletion_cluster_member_reads,
    _same_chrom_deletion_split_member_reads,
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


def _split_row(**kwargs):
    base = {
        "chrom": "chr22",
        "window_start": 10000,
        "window_end": 10100,
        "pos": 10040,
        "read_name": "split0",
        "clip_seq": "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
        "clip_side": "L",
        "clip_len": 40,
        "read_seq": "",
        "mei_hit": True,
    }
    base.update(kwargs)
    return base


def _indexed_fasta(tmp_path: Path, seq: str, chrom: str = "chr22") -> Path:
    fa = tmp_path / "ref.fa"
    fa.write_text(f">{chrom}\n{seq}\n", encoding="utf-8")
    pysam.faidx(str(fa))
    return fa


class TestDeletionSplitClipMeiDrop:
    def test_clip_matching_far_end_is_excluded(self, tmp_path: Path):
        clip = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        # 0-based 14999 == mate_pos 15000 (1-based cluster median).
        seq = ("N" * 14999) + clip + ("N" * 200)
        fa = _indexed_fasta(tmp_path, seq)
        dpe = pd.DataFrame(
            [
                _dpe_row(read_name="del0", mate_pos=15000, pos=10040),
                _dpe_row(read_name="del1", mate_pos=15010, pos=10041),
            ]
        )
        splits = pd.DataFrame(
            [
                _split_row(
                    read_name="hit",
                    clip_seq=clip,
                    # A conflicting DPE-schema value proves the production
                    # split-evidence column is the one being examined.
                    soft_clip_seq="T" * 40,
                ),
                _split_row(read_name="other", clip_seq="T" * 40),
            ]
        )
        members = _same_chrom_deletion_split_member_reads(splits, dpe, fa)
        assert set(members["read_name"]) == {"hit"}
        kept = _drop_deletion_cluster_reads(splits, members)
        assert set(kept["read_name"]) == {"other"}

    def test_matching_split_read_is_also_excluded_from_dpe_representation(self, tmp_path: Path):
        clip = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        seq = ("N" * 14999) + clip + ("N" * 200)
        fa = _indexed_fasta(tmp_path, seq)
        deletion_dpe = pd.DataFrame(
            [
                _dpe_row(read_name="del0", mate_pos=15000, pos=10040),
                _dpe_row(read_name="del1", mate_pos=15010, pos=10041),
            ]
        )
        splits = pd.DataFrame([_split_row(read_name="duplicated", clip_seq=clip)])
        members = _same_chrom_deletion_split_member_reads(splits, deletion_dpe, fa)
        dpe_mei_hits = pd.DataFrame(
            [
                _dpe_row(read_name="duplicated", mei_hit=True),
                _dpe_row(read_name="independent", mei_hit=True),
            ]
        )

        kept = _drop_deletion_cluster_reads(dpe_mei_hits, members)

        assert set(kept["read_name"]) == {"independent"}

    def test_reverse_complement_clip_is_excluded(self, tmp_path: Path):
        clip = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        seq = ("N" * 14999) + _revcomp(clip) + ("N" * 200)
        fa = _indexed_fasta(tmp_path, seq)
        dpe = pd.DataFrame(
            [
                _dpe_row(read_name="del0", mate_pos=15000, pos=10040),
                _dpe_row(read_name="del1", mate_pos=15010, pos=10041),
            ]
        )
        splits = pd.DataFrame([_split_row(read_name="hit", clip_seq=clip)])
        members = _same_chrom_deletion_split_member_reads(splits, dpe, fa)
        assert set(members["read_name"]) == {"hit"}

    def test_unrelated_clip_is_kept(self, tmp_path: Path):
        clip = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        seq = ("N" * 15200)
        fa = _indexed_fasta(tmp_path, seq)
        dpe = pd.DataFrame(
            [
                _dpe_row(read_name="del0", mate_pos=15000, pos=10040),
                _dpe_row(read_name="del1", mate_pos=15010, pos=10041),
            ]
        )
        splits = pd.DataFrame([_split_row(read_name="keep", clip_seq=clip)])
        members = _same_chrom_deletion_split_member_reads(splits, dpe, fa)
        assert members.empty

    def test_clip_elsewhere_in_deleted_interval_is_not_excluded(self, tmp_path: Path):
        clip = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        # The inferred opposite breakpoint is at 15,000. A copy at 12,000 is
        # inside the broad candidate-to-mate span but outside mate ±500 bp.
        seq = ("N" * 11999) + clip + ("N" * 3200)
        fa = _indexed_fasta(tmp_path, seq)
        dpe = pd.DataFrame(
            [
                _dpe_row(read_name="del0", mate_pos=15000, pos=10040),
                _dpe_row(read_name="del1", mate_pos=15010, pos=10041),
            ]
        )
        splits = pd.DataFrame([_split_row(read_name="interval_only", clip_seq=clip)])
        members = _same_chrom_deletion_split_member_reads(splits, dpe, fa)
        assert members.empty

    def test_missing_split_evidence_clip_column_is_not_treated_as_a_match(self, tmp_path: Path):
        clip = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        seq = ("N" * 14999) + clip + ("N" * 200)
        fa = _indexed_fasta(tmp_path, seq)
        dpe = pd.DataFrame(
            [
                _dpe_row(read_name="del0", mate_pos=15000, pos=10040),
                _dpe_row(read_name="del1", mate_pos=15010, pos=10041),
            ]
        )
        legacy_split = _split_row(read_name="legacy", clip_seq=clip)
        legacy_split.pop("clip_seq")
        legacy_split["soft_clip_seq"] = clip
        splits = pd.DataFrame([legacy_split])
        members = _same_chrom_deletion_split_member_reads(splits, dpe, fa)
        assert members.empty

    def test_no_fasta_does_not_drop(self, tmp_path: Path):
        dpe = pd.DataFrame(
            [
                _dpe_row(read_name="del0", mate_pos=15000, pos=10040),
                _dpe_row(read_name="del1", mate_pos=15010, pos=10041),
            ]
        )
        splits = pd.DataFrame([_split_row()])
        members = _same_chrom_deletion_split_member_reads(splits, dpe, tmp_path / "missing.fa")
        assert members.empty

    def test_identity_helper_exact_and_mismatch(self):
        q = "ACGT" * 5
        assert _best_ungapped_identity(q, "N" * 5 + q + "N" * 5) == 1.0
        assert _best_ungapped_identity(q, "G" * 40) < 0.3
        almost = q[:-1] + "C"
        ident = _best_ungapped_identity(q, almost)
        assert ident == 0.95

    def test_gapped_matcher_allows_small_query_insertion(self):
        reference_clip = "ACGTCAGTGCATGACCTAGCGTACCATGCTAGTCAGTGCA"
        query = reference_clip[:20] + "A" + reference_clip[20:]
        subject = "N" * 25 + reference_clip + "N" * 25
        assert _gapped_local_match(query, subject)

    def test_gapped_matcher_allows_small_query_deletion(self):
        reference_clip = "ACGTCAGTGCATGACCTAGCGTACCATGCTAGTCAGTGCA"
        query = reference_clip[:20] + reference_clip[21:]
        subject = "N" * 25 + reference_clip + "N" * 25
        assert _gapped_local_match(query, subject)

    def test_gapped_matcher_enforces_edit_similarity(self):
        query = "ACGTTGCA" * 5
        three_mismatches = list(query)
        for idx in (3, 17, 31):
            three_mismatches[idx] = "A" if three_mismatches[idx] != "A" else "C"
        subject = "N" * 25 + "".join(three_mismatches) + "N" * 25
        assert not _gapped_local_match(query, subject)
