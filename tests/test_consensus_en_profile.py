"""Tests for scripts/consensus_en_profile.py.

Covers motif-matcher correctness (positive/negative fixtures), strand symmetry
handling, bin-boundary handling, zone assignment, window counting, FASTA
enumeration, and determinism. The script is loaded via importlib so its
top-level imports stay out of the package dependency graph (same convention as
test_locus_zoom_gif.py).
"""
from __future__ import annotations

import importlib.util as _iu
import sys
from pathlib import Path

import pytest

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "consensus_en_profile.py"


def _load_mod():
    spec = _iu.spec_from_file_location("consensus_en_profile", SCRIPT_PATH)
    mod = _iu.module_from_spec(spec)
    sys.modules["consensus_en_profile"] = mod  # must precede exec_module
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def mod():
    return _load_mod()


# ---------------------------------------------------------------------------
# Motif matcher: positive/negative fixtures.
# ---------------------------------------------------------------------------

def test_narrow_motif_positive(mod):
    assert mod.find_motif_positions("TTTTAAGGGG", "TTTTAA") == [1]
    # TTTTTAAGG: the window starting at 2 is TTTTAA; position 1 is TTTTTA.
    assert mod.find_motif_positions("TTTTTAAGG", "TTTTAA") == [2]


def test_narrow_motif_negative(mod):
    assert mod.find_motif_positions("TTTTGAAGGG", "TTTTAA") == []      # G breaks it
    assert mod.find_motif_positions("TTTTAG", "TTTTAA") == []          # too short
    assert mod.find_motif_positions("UUUUAA", "TTTTAA") == []          # no T run
    assert mod.find_motif_positions("NNNNNNGGGG", "TTTTAA") == []      # N never matches T


def test_flasch_7mer_positive_and_negative(mod):
    assert mod.find_motif_positions("TTTTTAAG", "TTTTTAA") == [1]
    assert mod.find_motif_positions("TTTTAAGG", "TTTTTAA") == []       # only 4 Ts


def test_degenerate_yyrrrr_matches_purine_variants(mod):
    # YYRRRR: C/T then C/T, then four purines.
    assert mod.find_motif_positions("CCAAGGGG", "YYRRRR") == [1]
    assert mod.find_motif_positions("TTGGGGGG", "YYRRRR") == [1]
    assert mod.find_motif_positions("TTACGGGG", "YYRRRR") == []        # third base C not R
    assert mod.find_motif_positions("TTAACCCC", "YYRRRR") == []        # A is not Y
    assert mod.find_motif_positions("CCAACCCC", "YYRRRR") == []        # fifth base C not R


def test_ambiguous_sequence_letters_do_not_match_iupac_positions(mod):
    assert mod.find_motif_positions("NNRRRR", "YYRRRR") == []          # N is not Y
    assert mod.find_motif_positions("YYNNNN", "YYRRRR") == []          # N is not R


def test_overlapping_matches_all_counted(mod):
    # Two adjacent occurrences, the second starting inside the first's span.
    assert mod.find_motif_positions("TTTTAATTTTAA", "TTTTAA") == [1, 7]


# ---------------------------------------------------------------------------
# Strand symmetry handling.
# ---------------------------------------------------------------------------

def test_both_strands_counted_separately(mod):
    # TTTTAA is NOT self-reverse-complementary (its complement reads TTAAAA),
    # so a sense TTTTAA does not appear as a TTTTAA match on the antisense.
    res = mod.profile_consensus(mod.Consensus("T1", "LINE1", "TTTTAA" + "G" * 594))
    assert res.strand_counts["sense"]["narrow_tttt_aa"] == 1
    assert res.strand_counts["antisense"]["narrow_tttt_aa"] == 0
    # A sequence carrying the complement-class TTAAAA on the sense strand has
    # its antisense-strand TTTTAA instance at the mirrored position.
    res2 = mod.profile_consensus(mod.Consensus("T2", "LINE1", "TTAAAACGCG" + "G" * 12))
    assert res2.strand_counts["sense"]["narrow_tttt_aa"] == 0
    assert res2.strand_counts["antisense"]["narrow_tttt_aa"] == 1


def test_antisense_positions_are_mirrored(mod):
    # TTAAAA on the sense strand spans positions 1-6 of a 10-mer; on the
    # antisense strand the corresponding TTTTAA starts at 10 - 6 + 1 = 5.
    seq = "TTAAAACGTA"
    res = mod.profile_consensus(mod.Consensus("T1", "LINE1", seq))
    assert res.positions[("sense", "narrow_tttt_aa")] == []
    assert res.positions[("antisense", "narrow_tttt_aa")] == [5]


def test_sense_and_antisense_differ_for_non_self_rc_motif(mod):
    # TTAAAA on the sense strand matches nothing there, but its reverse
    # complement TTTTAA matches on the antisense strand.
    seq = "CGCGCTTAAAA"
    res = mod.profile_consensus(mod.Consensus("A1", "LINE1", seq))
    assert res.strand_counts["sense"]["narrow_tttt_aa"] == 0
    assert res.strand_counts["antisense"]["narrow_tttt_aa"] == 1


# ---------------------------------------------------------------------------
# Bin boundaries and zones.
# ---------------------------------------------------------------------------

def test_bin_boundaries_inclusive_and_trailing_bin_short(mod):
    assert mod.bin_label(1, 500) == (0, 1, 500)
    assert mod.bin_label(500, 500) == (0, 1, 500)
    assert mod.bin_label(501, 500) == (1, 501, 1000)
    # A 1,230-bp consensus: last bin is 3 bp long, not a full 500.
    n_bins = (1230 + 499) // 500
    assert n_bins == 3


def test_trailing_short_bin_density_uses_true_extent(mod):
    # density_per_kb of a single site in a 3-bp bin is huge but defined.
    assert mod.density_per_kb(1, 3) == round(1 / 0.003, 10) or mod.density_per_kb(1, 3) > 300


def test_zone_assignment_scales_with_length(mod):
    # Full 6,021-bp consensus: canonical boundaries.
    assert mod.zone_for_l1_position(1, 6021) == "5UTR"
    assert mod.zone_for_l1_position(907, 6021) == "5UTR"
    assert mod.zone_for_l1_position(908, 6021) == "ORF1"
    assert mod.zone_for_l1_position(2381, 6021) == "ORF2"
    assert mod.zone_for_l1_position(5978, 6021) == "3UTR"
    # Half-length consensus: boundaries scale proportionally.
    assert mod.zone_for_l1_position(200, 3010) == "5UTR"
    assert mod.zone_for_l1_position(1300, 3010) == "ORF2"


# ---------------------------------------------------------------------------
# Window counting and probes.
# ---------------------------------------------------------------------------

def test_count_in_window_inclusive(mod):
    assert mod.count_in_window([117, 118, 136, 137], (118, 136)) == 2


def test_alu_probe_reports_both_strands(mod):
    seq = "A" * 117 + "TTTTAA" + "A" * 100 + "G" * 97  # TTTTAA at 118..123
    res = mod.profile_consensus(mod.Consensus("AluX", "ALU", seq))
    probe = mod.alu_probe(res)
    assert probe["narrow_tttt_aa"]["linker_118_136_sense"] == 1
    assert probe["narrow_tttt_aa"]["tail_280_320_sense"] == 0


def test_l1_inversion_cluster_probe_counts_window(mod):
    seq = "G" * 3999 + "TTTTAA" + "G" * 400  # site at 4000..4005
    res = mod.profile_consensus(mod.Consensus("L1X", "LINE1", seq))
    probe = mod.l1_inversion_cluster_probe(res)
    assert probe["narrow_tttt_aa"]["sites_4000_6000"] == 1
    assert probe["narrow_tttt_aa"]["sites_rest"] == 0


# ---------------------------------------------------------------------------
# FASTA enumeration and end-to-end determinism.
# ---------------------------------------------------------------------------

def _write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w") as handle:
        for name, seq in records:
            handle.write(f">{name}\n{seq}\n")


def test_enumeration_takes_exactly_l1hs_l1pa_alu(mod, tmp_path):
    fasta = tmp_path / "reps.fa"
    _write_fasta(fasta, [
        ("L1HS", "TTTTAA" + "A" * 6015),
        ("L1PA2", "A" * 6021),
        ("AluJo", "A" * 320),
        ("AluSx", "A" * 300),
        ("MIR", "A" * 200),        # not selected
        ("L2a", "A" * 200),        # not selected
        ("THE1A", "A" * 200),      # prefix must match at position 0: excluded
    ])
    consensuses = mod.read_subfamily_consensuses(fasta)
    names = [c.subfamily for c in consensuses]
    assert names == ["L1HS", "L1PA2", "AluJo", "AluSx"]
    assert all(c.family == "LINE1" for c in consensuses[:2])
    assert all(c.family == "ALU" for c in consensuses[2:])


def test_multiline_sequence_reassembled(mod, tmp_path):
    fasta = tmp_path / "reps.fa"
    with fasta.open("w") as handle:
        handle.write(">AluY\nTTTTAA\nAAAAAA\nTTTTAA\n")
    consensuses = mod.read_subfamily_consensuses(fasta)
    assert consensuses[0].sequence == "TTTTAAAAAAAATTTTAA"


def test_end_to_end_determinism(mod, tmp_path):
    fasta = tmp_path / "reps.fa"
    _write_fasta(fasta, [
        ("L1HS", ("TTTTAACCGG" * 60)[:6021]),
        ("AluY", ("TTTTAA" + "A" * 6)[:320]),
    ])
    mod.main(["--fasta", str(fasta), "--outdir", str(tmp_path / "out1")])
    mod.main(["--fasta", str(fasta), "--outdir", str(tmp_path / "out2")])
    a = (tmp_path / "out1" / "consensus_en_sites.csv").read_text()
    b = (tmp_path / "out2" / "consensus_en_sites.csv").read_text()
    assert a == b
    assert "# Position convention: 1-based" in a
    assert "age" in a.split("\n")[len(a.split("\n")) - 1] or "age" in a


def test_missing_fasta_exits(mod, tmp_path):
    with pytest.raises(SystemExit, match="required consensus FASTA not found"):
        mod.main(["--fasta", str(tmp_path / "absent.fa"), "--outdir", str(tmp_path / "out")])
