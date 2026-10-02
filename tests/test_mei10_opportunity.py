"""Tests for the ten-genome reference-opportunity and consensus helpers."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "mei10_opportunity.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("mei10_opportunity", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


m = _load_module()


# --------------------------------------------------------------------------
# interval algebra
# --------------------------------------------------------------------------
def test_merge_sorts_unions_and_drops_empties():
    assert m.merge([(10, 20), (0, 5), (4, 12)]) == [(0, 20)]
    assert m.merge([(0, 10), (10, 20)]) == [(0, 20)]  # touching merges
    assert m.merge([(5, 5), (7, 6)]) == []


def test_subtract_removes_only_the_covered_bases():
    assert m.subtract([(0, 100)], [(20, 30)]) == [(0, 20), (30, 100)]
    assert m.subtract([(0, 100)], [(0, 100)]) == []
    assert m.subtract([(0, 100)], [(90, 200)]) == [(0, 90)]
    assert m.subtract([(0, 100)], [(-50, 10)]) == [(10, 100)]
    # a span straddling both edges of several masks
    assert m.subtract([(0, 100)], [(10, 20), (40, 50), (80, 90)]) == [
        (0, 10),
        (20, 40),
        (50, 80),
        (90, 100),
    ]


def test_subtract_is_conservative_never_invents_bases():
    span = [(0, 1000)]
    masks = [(i * 37, i * 37 + 11) for i in range(27)]
    kept = m.subtract(span, masks)
    total = sum(e - s for s, e in kept)
    assert 0 < total < 1000
    assert total == 1000 - sum(min(1000, e) - max(0, s) for s, e in m.merge(masks))


def test_overlap_counts_shared_bases():
    spans = [(0, 100), (200, 300)]
    assert m.overlap(spans, 0, 100) == 100
    assert m.overlap(spans, 50, 250) == 50 + 50
    assert m.overlap(spans, 400, 500) == 0


def test_bin_opportunity_absolute_bins():
    # edges 0,10,20,30 -> three 10 bp bins; one 25 bp span covers all three
    out = m.bin_opportunity([(0, 25)], [0, 10, 20, 30])
    assert out.tolist() == [10, 10, 5]
    # nothing outside the edge range counts
    assert m.bin_opportunity([(100, 200)], [0, 10, 20, 30]).tolist() == [0, 0, 0]


def test_bin_opportunity_relative_bins_scale_by_span():
    # every span contributes exactly len(span)/n_bins to each bin, so a 500 bp
    # and a 1000 bp host contribute 250 and 500 per bin
    edges = [0, 50, 100]
    out = m.bin_opportunity([(0, 500), (0, 1000)], edges, relative=True)
    assert out.tolist() == [250 + 500, 250 + 500]


def test_consensus_bin_edges_covers_full_length():
    assert m.consensus_bin_edges(1200, 500) == [0, 500, 1000, 1200]
    assert m.consensus_bin_edges(1000, 500) == [0, 500, 1000]


def test_normalize_family_matches_project_rules():
    assert m.normalize_family("L1MC5a LINE L1") == "LINE1"
    assert m.normalize_family("AluY SINE Alu") == "ALU"
    assert m.normalize_family("SVA_F Retroposon SVA") == "SVA"
    assert m.normalize_family("TAR1 Satellite telo") == ""


# --------------------------------------------------------------------------
# FASTA reading and gap scanning
# --------------------------------------------------------------------------
def _write_fasta(path: Path, seq: str, line_bases: int = 10) -> None:
    lines = []
    for i in range(0, len(seq), line_bases):
        lines.append(seq[i : i + line_bases])
    path.write_text(">chr1\n" + "\n".join(lines) + "\n", encoding="utf-8")


def test_fasta_chunk_reads_wrapped_sequence(tmp_path):
    seq = "ACGT" * 25 + "NN" + "TGCA" * 5  # 120 bases
    fasta = tmp_path / "ref.fa"
    _write_fasta(fasta, seq, line_bases=10)
    # write a matching .fai by hand: length, offset, line_bases, line_width
    (tmp_path / "ref.fa.fai").write_text(
        f"chr1\t{len(seq)}\t6\t10\t11\n", encoding="utf-8"
    )
    meta = m.read_fai(fasta)
    with fasta.open("rb") as handle:
        assert m.fasta_chunk(handle, meta["chr1"], 0, 120) == seq.encode()
        assert m.fasta_chunk(handle, meta["chr1"], 8, 30) == seq[8:30].encode()
        # a slice that straddles a line break
        assert m.fasta_chunk(handle, meta["chr1"], 5, 25) == seq[5:25].encode()


def test_fasta_chunk_rejects_out_of_range(tmp_path):
    fasta = tmp_path / "ref.fa"
    _write_fasta(fasta, "ACGT" * 10, line_bases=10)
    (tmp_path / "ref.fa.fai").write_text("chr1\t40\t6\t10\t11\n", encoding="utf-8")
    meta = m.read_fai(fasta)
    with fasta.open("rb") as handle:
        with pytest.raises(ValueError):
            m.fasta_chunk(handle, meta["chr1"], 30, 10)
        with pytest.raises(ValueError):
            m.fasta_chunk(handle, meta["chr1"], 0, 41)


def test_scan_gaps_finds_every_n_run(tmp_path):
    seq = "ACGT" * 10 + "N" * 7 + "ACGT" * 10 + "NNNN" + "ACGT" * 10
    fasta = tmp_path / "ref.fa"
    _write_fasta(fasta, seq, line_bases=10)
    (tmp_path / "ref.fa.fai").write_text(
        f"chr1\t{len(seq)}\t6\t10\t11\n", encoding="utf-8"
    )
    meta = m.read_fai(fasta)
    gaps = m.scan_gaps(fasta, meta, ["chr1"], min_run=1, chunk_bp=17)
    assert gaps["chr1"].tolist() == [[40, 47], [87, 91]]
    # min_run filters short runs
    assert m.scan_gaps(fasta, meta, ["chr1"], min_run=8, chunk_bp=17)[
        "chr1"
    ].tolist() == []


def test_host_segments_removes_gaps_from_a_host(tmp_path):
    host = m.Host("chr1", 0, 100, "+", "AluY", "ALU")
    gaps = np.asarray([[40, 50]], dtype=np.int64)
    assert m.host_segments(host, gaps) == [(0, 40), (50, 100)]
    # an extra mask removes a further slice
    assert m.host_segments(host, gaps, masks=[(0, 10)]) == [(10, 40), (50, 100)]


# --------------------------------------------------------------------------
# consensus projection
# --------------------------------------------------------------------------
def test_project_to_consensus_maps_offset_forward():
    consensus = "ACGT" * 25  # 100 bp
    host = consensus[:60]
    got = m.project_to_consensus(host, consensus, host_offset=30)
    assert got is not None
    offset, orientation = got
    assert orientation == "forward"
    assert offset == 30


def test_project_to_consensus_detects_reverse_orientation():
    # a consensus that is NOT its own reverse complement, so orientation is decidable
    consensus = "AAAACCCCGGGGTTTT" * 10  # 160 bp
    rc = consensus[::-1].translate(str.maketrans("ACGT", "TGCA"))
    assert rc != consensus
    got = m.project_to_consensus(rc, consensus, host_offset=10)
    assert got is not None
    offset, orientation = got
    assert orientation == "reverse"
    assert offset == len(consensus) - 1 - 10


def test_project_to_consensus_returns_none_for_unrelated_sequence():
    consensus = "ACGT" * 50
    junk = ("C" * 5 + "G" * 3 + "T" * 7) * 6
    assert m.project_to_consensus(junk, consensus, host_offset=10) is None


def test_assign_subfamily_prefers_the_closest_consensus():
    a = "ACGTACGTAA" * 6
    b = "TTTTGGGGCC" * 6
    kmers = {"L1PA2": m.kmer_set(a), "L1PB1": m.kmer_set(b)}
    assert m.assign_subfamily(a, kmers)[0] == "L1PA2"
    assert m.assign_subfamily(b, kmers)[0] == "L1PB1"
    assert m.assign_subfamily("", kmers) == (None, 0.0)
