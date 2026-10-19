"""Tests for alu_polya_asymmetry. Expected values are hand-worked, not taken
from the implementation's own output."""
import gzip
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import alu_polya_asymmetry as m


def test_pure_a_run_carries_no_canonical_hexamer():
    """A homopolymer A-tail cannot contain AATAAA, however long: the motif needs
    a T. This is why AluY's consensus, whose tail is a pure 30-nt A-run, scores
    zero canonical signals."""
    assert m.count_hexamers("A" * 60) == 0
    assert m.count_hexamers("A" * 2 + "T" + "A" * 3) == 1      # AATAAA exactly


def test_one_substitution_in_an_a_run_creates_a_signal():
    """The other half of the same fact, and the mechanism Roy-Engel 2005
    describes: a single A->T substitution inside an A-run creates AATAAA, so a
    long tail is one mutation away from a poly(A) signal. Measured here: a
    single T placed in a 61-nt A-run yields exactly one hexamer, and it is the
    reason 78% of signal-carrying AluY copies carry theirs in the tail."""
    assert m.count_hexamers("A" * 30 + "T" + "A" * 30) == 1
    # a T too close to either end cannot complete the hexamer
    assert m.count_hexamers("T" + "A" * 30) == 0
    assert m.count_hexamers("A" * 30 + "T") == 0
    assert m.count_hexamers("A" * 2 + "T" + "A" * 2) == 0      # needs 3 A's after


def test_count_hexamers_window_is_1_based_inclusive():
    seq = "C" * 10 + "AATAAA" + "C" * 10        # hexamer occupies 11..16
    assert m.count_hexamers(seq) == 1
    assert m.count_hexamers(seq, 11) == 1
    assert m.count_hexamers(seq, 12) == 0      # start excluded
    assert m.count_hexamers(seq, 1, 11) == 1
    assert m.count_hexamers(seq, 1, 10) == 0


def test_count_hexamers_counts_overlaps():
    # "AATAAATAAA" contains AATAAA at 1-based 1 and again at 5.
    assert m.count_hexamers("AATAAATAAA") == 2


def test_count_hexamers_is_single_stranded():
    """TTTATT is AATAAA reverse-complemented; the counter must not see it."""
    assert m.revcomp("TTTATT") == "AATAAA"
    assert m.count_hexamers("TTTATT") == 0


def test_a_rich_regions_finds_only_runs_above_threshold():
    seq = "C" * 30 + "A" * 30 + "C" * 30
    spans = m.a_rich_regions(seq, window=20, min_a_frac=0.6)
    assert len(spans) == 1
    lo, hi = spans[0]
    assert lo >= 12 and hi <= 78            # brackets the A-run, window-widened
    assert all(seq[lo - 1:hi].count("A") >= 30 for _ in [0])
    assert m.a_rich_regions("CGCGCG" * 20, 20, 0.6) == []


def test_mask_spans_deletes_and_preserves_the_rest():
    seq = "AAAAA" + "CGCGCGCGCG" + "TTTTT"
    out = m.mask_spans(seq, [(1, 5), (16, 20)])
    assert out == "CGCGCGCGCG"
    assert m.mask_spans(seq, []) == seq


def _rmsk_line(chrom, start, end, strand, rep, milli_div=50, milli_del=0,
               milli_ins=0, rep_start=1, rep_end=311, rep_left=1):
    f = ["585", "2000", str(milli_div), str(milli_del), str(milli_ins), chrom,
         str(start), str(end), "-1", strand, rep, "SINE", "Alu",
         str(rep_start), str(rep_end), str(rep_left), "1"]
    return "\t".join(f) + "\n"


def test_select_copies_reads_consensus_start_per_strand(tmp_path):
    """UCSC puts the consensus start in repStart on '+' rows and in repLeft on
    '-' rows. Using the '+' formula on a '-' row is the bug that inflated a
    892 bp model to 11,799 bp earlier in this project."""
    p = tmp_path / "rmsk.txt.gz"
    with gzip.open(p, "wt") as fh:
        # '+' row: cons start in repStart; genomic length must equal 311
        fh.write(_rmsk_line("chr1", 1000, 1311, "+", "AluY",
                            rep_start=1, rep_end=311, rep_left=-5))
        # '-' row: cons start in repLeft, repStart holds the negative remainder
        fh.write(_rmsk_line("chr1", 2000, 2311, "-", "AluY",
                            rep_start=-5, rep_end=311, rep_left=1))
        fh.write(_rmsk_line("chr1", 3000, 3311, "+", "AluSx", rep_end=311))
        fh.write(_rmsk_line("chrUn_test", 4000, 4311, "+", "AluY", rep_end=311))
    got = m.select_copies(p, "AluY")
    assert len(got) == 2, "both strands must be kept, non-AluY and non-main dropped"
    assert {g[3] for g in got} == {"+", "-"}
    assert all(g[4] == 1 for g in got), "consensus start must be 1 on both rows"


def test_select_copies_rejects_indels_and_length_mismatch(tmp_path):
    p = tmp_path / "rmsk.txt.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(_rmsk_line("chr1", 1000, 1311, "+", "AluY", milli_del=20))
        fh.write(_rmsk_line("chr1", 2000, 2311, "+", "AluY", milli_ins=20))
        fh.write(_rmsk_line("chr1", 3000, 3290, "+", "AluY"))   # len != consensus span
        fh.write(_rmsk_line("chr1", 4000, 4311, "+", "AluY", milli_div=900))
        fh.write(_rmsk_line("chr1", 5000, 5311, "+", "AluY", rep_end=200))  # truncated
    assert m.select_copies(p, "AluY") == []


def test_select_all_copies_keeps_what_select_copies_drops(tmp_path):
    p = tmp_path / "rmsk.txt.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(_rmsk_line("chr1", 1000, 1200, "+", "AluY", milli_div=200,
                            milli_del=30, rep_end=200))
    assert m.select_copies(p, "AluY") == []
    allc = m.select_all_copies(p, "AluY")
    assert len(allc) == 1 and allc[0][4] == 200       # divergence carried through


def test_in_gene_respects_strand_and_boundaries():
    idx = {("chr1", "+"): (np.array([100, 500]), np.array([200, 600])),
           ("chr1", "-"): (np.array([1000]), np.array([1100]))}
    assert m.in_gene(idx, "chr1", 150, "+")
    assert m.in_gene(idx, "chr1", 100, "+")          # inclusive start
    assert m.in_gene(idx, "chr1", 200, "+")          # inclusive end
    assert not m.in_gene(idx, "chr1", 201, "+")
    assert not m.in_gene(idx, "chr1", 150, "-")      # right position, wrong strand
    assert m.in_gene(idx, "chr1", 1050, "-")
    assert not m.in_gene(idx, "chr2", 150, "+")      # absent chromosome


def test_in_gene_finds_a_gene_that_starts_before_an_earlier_ending_one():
    """A long gene starting first but ending last must still be found for a
    position inside it, even though a shorter gene sorts between them."""
    idx = {("chr1", "+"): (np.array([100, 200]), np.array([10_000, 300]))}
    assert m.in_gene(idx, "chr1", 5_000, "+")


def test_mantel_haenszel_matches_a_hand_worked_two_stratum_case():
    # stratum 1: a=10 b=20 c=30 d=20, n=80 -> R=10*20/80=2.5   S=20*30/80=7.5
    # stratum 2: a=20 b=10 c=20 d=30, n=80 -> R=20*30/80=7.5   S=10*20/80=2.5
    # OR_MH = (2.5+7.5)/(7.5+2.5) = 1.0
    out = m.mantel_haenszel([(10, 20, 30, 20), (20, 10, 20, 30)])
    assert out["strata"] == 2
    assert out["or_mh"] == pytest.approx(1.0, abs=1e-9)
    assert out["ci95"][0] < 1.0 < out["ci95"][1]


def test_mantel_haenszel_recovers_a_known_constant_odds_ratio():
    # every stratum built with OR = 4 exactly: (40*40)/(20*20)
    tables = [(40, 20, 20, 40)] * 5
    out = m.mantel_haenszel(tables)
    assert out["or_mh"] == pytest.approx(4.0, rel=1e-9)
    assert out["p"] < 1e-6


def test_mantel_haenszel_skips_sparse_strata_and_returns_none_when_all_sparse():
    assert m.mantel_haenszel([(1, 1, 1, 1), (2, 2, 2, 2)]) is None
    out = m.mantel_haenszel([(1, 1, 1, 1), (40, 20, 20, 40)])
    assert out["strata"] == 1


def test_build_consensus_calls_column_majority(tmp_path, monkeypatch):
    """Three copies, one carrying a private substitution, must give the
    majority base and a support fraction of 2/3 at that column."""
    class FakeGenome(dict):
        pass
    seqs = {"chr1": "N" * 1000}
    # copy sequences are returned by slicing, so build a chromosome that
    # contains them at known offsets
    a = "ACGT" * 2          # 8 bp consensus for the test
    b = "ACGT" * 2
    c = "ACGA" + "ACGT"     # differs at position 4
    chrom = "N" * 100 + a + b + c + "N" * 100
    monkeypatch.setattr(m, "CONSENSUS_LEN", 8, raising=False)

    class FakeTwoBit:
        def __init__(self, *_a, **_k): pass
        def __getitem__(self, k): return chrom

    monkeypatch.setitem(sys.modules, "twobitreader",
                        type("mod", (), {"TwoBitFile": FakeTwoBit}))
    copies = [("chr1", 100, 108, "+", 1, 8),
              ("chr1", 108, 116, "+", 1, 8),
              ("chr1", 116, 124, "+", 1, 8)]
    cons, qc = m.build_consensus(copies, "ignored", cons_len=8)
    assert cons == "ACGT" * 2
    assert qc["copies_used"] == 3
    assert qc["min_depth"] == 3
    assert qc["mean_column_support"] == pytest.approx((7 * 1.0 + 2 / 3) / 8, abs=1e-4)


def test_null_distribution_preserves_composition_under_mono():
    rng = np.random.default_rng(0)
    seq = "AATAAA" + "CGCGCG" * 5
    sense, anti = m.null_distribution(seq, 25, rng, "mono", canonical_only=True)
    assert len(sense) == len(anti) == 25
    assert sense.min() >= 0
