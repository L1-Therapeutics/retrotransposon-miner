"""Tests for transduction source tagging (transduction_tags.py).

The fixture encodes the Damert 2009 scenario: two SVA insertions share a
unique upstream segment from one source element's flank; a third SVA shares
a *different* equally-strong segment from another source, which must produce
an unambiguous assignment only when one source truly dominates.
"""

from __future__ import annotations

import random

import pytest

from retro_miner.transduction_tags import (
    DEFAULT_MIN_SEED_SITES,
    FLANK_3P,
    FLANK_5P,
    TransducedSegment,
    annotate_callset,
    build_source_index,
    build_transduction_groups,
    extract_source_flank,
    match_call_to_source,
    uniqueness_score,
)

Complement = str.maketrans("ACGT", "TGCA")


def _revcomp(seq: str) -> str:
    return seq.translate(Complement)[::-1]

# --- deterministic random sequence -------------------------------------------

_RNG = random.Random(20260930)


def _random_seq(n: int, rng: random.Random = _RNG) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


# --- uniqueness screen ---------------------------------------------------------

def test_uniqueness_score_bounds():
    assert uniqueness_score("") == 0.0
    # A 60-mer of pure ACGT: all 25-mer windows unique.
    assert uniqueness_score(_random_seq(60)) == pytest.approx(1.0)
    # Perfect 25-mer repeat: zero unique windows.
    unit = _random_seq(25)
    assert uniqueness_score(unit * 3) == 0.0


# --- flank extraction: strand correctness ----------------------------------------

class _DictGenome:
    """Tiny dict-backed SequenceFetcher; coordinates 1-based inclusive."""

    def __init__(self, seqs: dict[str, str]) -> None:
        self.seqs = {k: v.upper() for k, v in seqs.items()}

    def __call__(self, chrom: str, start: int, end: int) -> str:
        # 1-based inclusive -> 0-based half-open
        return self.seqs[chrom][start - 1 : end]


def test_extract_l1_plus_strand_uses_downstream_flank():
    upstream = _random_seq(100)
    body = _random_seq(200)
    downstream = _random_seq(300)
    genome = _DictGenome({"chr1": upstream + body + downstream})
    el = extract_source_flank(
        "L1_source", "LINE1", "chr1", "+", 100, 300, genome
    )
    assert el is not None
    assert el.flank_dir == FLANK_3P
    # Flank must be downstream: starts at element_end0=300 (0-based).
    assert el.flank_start0 == 300
    assert el.flank_seq == downstream[: el.flank_len]


def test_extract_l1_minus_strand_revcomp_upstream():
    upstream = _random_seq(300)
    body = _random_seq(200)
    downstream = _random_seq(100)
    genome = _DictGenome({"chr1": upstream + body + downstream})
    el = extract_source_flank(
        "L1_source", "LINE1", "chr1", "-", 300, 500, genome
    )
    assert el is not None
    assert el.flank_dir == FLANK_3P
    # '-' L1 reads through its 3' end into reference *upstream* of the body,
    # and the flank is the reverse complement of that forward-strand sequence.
    assert el.flank_start0 == 0
    assert el.flank_seq == _revcomp(upstream[-el.flank_len :])


def test_extract_sva_plus_strand_uses_upstream_flank():
    upstream = _random_seq(300)
    body = _random_seq(200)
    downstream = _random_seq(100)
    genome = _DictGenome({"chr1": upstream + body + downstream})
    el = extract_source_flank("SVA_source", "SVA", "chr1", "+", 300, 500, genome)
    assert el is not None
    assert el.flank_dir == FLANK_5P
    assert el.flank_end0 == 300
    assert el.flank_seq == upstream[-el.flank_len :]


def test_extract_sva_minus_strand_revcomp_downstream():
    upstream = _random_seq(100)
    body = _random_seq(200)
    downstream = _random_seq(300)
    genome = _DictGenome({"chr1": upstream + body + downstream})
    el = extract_source_flank("SVA_source", "SVA", "chr1", "-", 100, 300, genome)
    assert el is not None
    assert el.flank_dir == FLANK_5P
    assert el.flank_start0 == 300
    assert el.flank_seq == _revcomp(downstream[: el.flank_len])


def test_extract_rejects_repetitive_flank():
    # ATATAT... flank: uniqueness below the gate -> no source indexed.
    genome = _DictGenome({"chr1": "AT" * 400})
    el = extract_source_flank("rep_source", "SVA", "chr1", "+", 300, 500, genome)
    assert el is None


def test_extract_rejects_unknown_family():
    genome = _DictGenome({"chr1": _random_seq(600)})
    assert (
        extract_source_flank("alu_src", "ALU", "chr1", "+", 300, 500, genome)
        is None
    )


# --- end-to-end Damert scenario ----------------------------------------------------


def _make_source_and_call(family: str, rng: random.Random) -> tuple[str, str, str, str]:
    """Return (genome_seq, flank_seed_seq, call_extra_seq, body_seq).

    The call's extra sequence embeds 120 bp of the source flank; the rest of
    the call sequence is unrelated random DNA.
    """
    upstream = _random_seq(600, rng)
    body = _random_seq(500, rng)
    downstream = _random_seq(300, rng)
    genome = upstream + body + downstream
    flank_seq = upstream[-600:]  # SVA '+' flank
    shared = flank_seq[200:320]  # 120 bp unique piece inside the flank
    call_extra = _random_seq(80, rng) + shared + _random_seq(60, rng)
    return genome, flank_seq, call_extra, body


def test_two_calls_share_source_group_link():
    rng = random.Random(77)
    genome, _flank, call1_extra, body = _make_source_and_call("SVA", rng)
    genome = genome  # same source locus for both calls by construction
    _flank2, call2_extra, _body2 = _flank, _random_seq(90, rng) + _flank[200:320] + _random_seq(40, rng), body
    fetch = _DictGenome({"chrA": genome})
    src = extract_source_flank("SRC1", "SVA", "chrA", "+", 600, 1100, fetch)
    assert src is not None and src.flank_len > 0
    index = build_source_index([src])
    assert "SRC1" in index and len(index) == 1

    hit1 = match_call_to_source("call1", "SVA", call1_extra, index)
    hit2 = match_call_to_source("call2", "SVA", call2_extra, index)
    assert hit1 is not None and hit2 is not None
    assert hit1.source_id == hit2.source_id == "SRC1"
    assert hit1.unique and hit2.unique
    assert hit1.n_seed_sites >= DEFAULT_MIN_SEED_SITES
    # Both calls match the same flank window (criterion 2: overlapping
    # transduced sequence).
    assert hit1.flank_span0 < hit2.flank_span1
    assert hit2.flank_span0 < hit1.flank_span1

    groups = build_transduction_groups(
        [("call1", hit1.source_id), ("call2", hit2.source_id)]
    )
    assert len(groups) == 1
    assert groups[0].size == 2
    assert groups[0].source_ids == {"SRC1"}


def test_ambiguous_tie_refuses_to_name_source():
    rng = random.Random(99)
    # Two distinct sources; the call embeds equal-length unique segments from
    # each flank — a genuine tie.
    flank_a = _random_seq(400, rng)
    flank_b = _random_seq(400, rng)
    body_a = _random_seq(300, rng)
    body_b = _random_seq(300, rng)
    fetch = _DictGenome(
        {
            "chrA": _random_seq(50, rng) + body_a + flank_a,
            "chrB": _random_seq(50, rng) + body_b + flank_b,
        }
    )
    src_a = extract_source_flank("SA", "SVA", "chrA", "+", 50, 350, fetch)
    src_b = extract_source_flank("SB", "SVA", "chrB", "+", 50, 350, fetch)
    assert src_a and src_b
    index = build_source_index([src_a, src_b])
    call_extra = flank_a[100:200] + _random_seq(30, rng) + flank_b[100:200]
    hit = match_call_to_source("callx", "SVA", call_extra, index)
    assert hit is None  # tie at the top -> no source named


def test_different_family_never_matches():
    rng = random.Random(5)
    genome, _flank, call_extra, _body = _make_source_and_call("SVA", rng)
    fetch = _DictGenome({"chrA": genome})
    src = extract_source_flank("SRC_SVA", "SVA", "chrA", "+", 600, 1100, fetch)
    assert src is not None
    index = build_source_index([src])
    assert match_call_to_source("c1", "LINE1", call_extra, index) is None


def test_below_seed_floor_returns_none():
    rng = random.Random(11)
    genome, flank, _call_extra, _body = _make_source_and_call("SVA", rng)
    fetch = _DictGenome({"chrA": genome})
    src = extract_source_flank("SRC1", "SVA", "chrA", "+", 600, 1100, fetch)
    assert src is not None
    index = build_source_index([src])
    # A short shared piece (18 bp < one k-mer window) cannot produce seeds.
    weak_call = _random_seq(40, rng) + flank[200:218] + _random_seq(30, rng)
    assert match_call_to_source("weak", "SVA", weak_call, index) is None


def test_annotate_callset_preserves_order_and_none():
    rng = random.Random(23)
    genome, _flank, call_extra, _body = _make_source_and_call("SVA", rng)
    fetch = _DictGenome({"chrA": genome})
    src = extract_source_flank("SRC1", "SVA", "chrA", "+", 600, 1100, fetch)
    assert src is not None
    index = build_source_index([src])
    calls = [
        {"call_id": "good", "family": "SVA", "extra_seq": call_extra},
        {"call_id": "none_seq", "family": "SVA", "extra_seq": ""},
        {"call_id": "no_match", "family": "SVA", "extra_seq": _random_seq(150, rng)},
    ]
    out = annotate_callset(calls, index)
    assert len(out) == 3
    assert out[0] is not None and out[0].source_id == "SRC1"
    assert out[1] is None and out[2] is None


def test_transduced_segment_spans_sane_offsets():
    rng = random.Random(31)
    genome, _flank, call_extra, _body = _make_source_and_call("SVA", rng)
    fetch = _DictGenome({"chrA": genome})
    src = extract_source_flank("SRC1", "SVA", "chrA", "+", 600, 1100, fetch)
    assert src is not None
    index = build_source_index([src])
    hit = match_call_to_source("c", "SVA", call_extra, index)
    assert hit is not None
    assert 0 <= hit.call_span0 < hit.call_span1 <= len(call_extra)
    assert 0 <= hit.flank_span0 < hit.flank_span1 <= src.flank_len
    assert isinstance(hit, TransducedSegment)


# --- L1 3' transduction end-to-end -------------------------------------------------


def test_l1_three_prime_transduction_link():
    rng = random.Random(47)
    upstream = _random_seq(300, rng)
    body = _random_seq(1000, rng)
    downstream = _random_seq(800, rng)
    genome = upstream + body + downstream
    fetch = _DictGenome({"chrL": genome})
    src = extract_source_flank("L1SRC", "LINE1", "chrL", "+", 300, 1300, fetch)
    assert src is not None and src.flank_dir == FLANK_3P
    index = build_source_index([src])
    shared = src.flank_seq[100:260]  # 160 bp of the source's 3' flank
    call_extra = _random_seq(50, rng) + shared + _random_seq(70, rng)
    hit = match_call_to_source("lc1", "LINE1", call_extra, index)
    assert hit is not None and hit.source_id == "L1SRC"


# --- group edge cases ----------------------------------------------------------------


def test_groups_isolated_pairs_and_singletons():
    groups = build_transduction_groups(
        [("c1", "S1"), ("c2", "S1"), ("c3", "S2")]
    )
    by_id = {g.group_id: g for g in groups}
    assert len(groups) == 2  # {c1,c2,S1} and {c3,S2}
    sizes = sorted(g.size for g in groups)
    assert sizes == [1, 2]
    # Singleton group still records its source.
    singleton = [g for g in groups if g.size == 1][0]
    assert singleton.source_ids == {"S2"}
    assert by_id  # group ids unique
    assert len({g.group_id for g in groups}) == 2


def test_empty_assignments():
    assert build_transduction_groups([]) == []
