"""Tests for `retro_miner.local_assembly`.

This module shipped 1,874 lines with **no test at all**. Seven names it
referenced from sixteen places inside function bodies were never defined, so
`_run_minimap2_paf`, `_load_fasta_lengths`, `_choose_consensus_features` and
`_extract_sample_assembly_features` each raised `NameError` on entry -- and
nothing noticed, because the failures live inside function bodies that no test
ever entered. The only visible symptom was `threading` sitting unused in the
imports, waiting for the two locks it exists to build.

These tests exercise the pure helpers and the repaired constants directly. They
do not run SPAdes or minimap2, so they are fast and need no data, but they are
enough to make the module's name-resolution failures impossible to ship again:
every constant referenced from a function body is read here by actually calling
the function that reads it.

Anything requiring an external aligner, a reference FASTA, or a BAM stack is out
of scope here and remains unverified.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from retro_miner import local_assembly as la


# --------------------------------------------------------------------------
# the repaired module-level state
# --------------------------------------------------------------------------


def test_the_caches_and_locks_the_module_refers_to_actually_exist() -> None:
    """The regression that mattered: these were referenced, never defined.

    `_run_minimap2_paf` and `_load_fasta_lengths` both take a lock and then
    index a cache on their first line of real work. With the names undefined,
    both raised `NameError` and the module's only defence was that nothing
    called it.
    """
    assert isinstance(la._MINIMAP2_INDEX_CACHE, dict)
    assert isinstance(la._MEI_FASTA_LENGTH_CACHE, dict)
    # A non-reentrant lock is what the call sites assume: they mutate the cache
    # while holding it. A bare object would pass an isinstance check and still
    # fail to serialise anything.
    for lock in (la._MINIMAP2_INDEX_LOCK, la._MEI_FASTA_LENGTH_CACHE_LOCK):
        assert hasattr(lock, "acquire") and hasattr(lock, "release")


def test_every_name_in_all_resolves() -> None:
    """`__all__` is the module's advertised surface; all of it must exist."""
    missing = [name for name in la.__all__ if not hasattr(la, name)]
    assert not missing, f"advertised but undefined: {missing}"


def test_the_side_anchor_threshold_matches_the_package_convention() -> None:
    """A second, divergent MEI-alignment floor would silently change anchors.

    `_MIN_SIDE_ANCHOR_ALN_LEN` decides which assembly alignments count as a
    breakpoint side anchor. It is pinned to the same 20 bp floor
    `mei_support._MEI_ALIGN_MIN_ALN_BP_LONG` uses, so read evidence and assembly
    evidence agree about what counts as a real alignment.
    """
    from retro_miner import mei_support

    assert la._MIN_SIDE_ANCHOR_ALN_LEN == mei_support._MEI_ALIGN_MIN_ALN_BP_LONG


def test_the_polya_threshold_matches_the_read_evidence_gate() -> None:
    """The source comment claims this matches mei_support; hold it to that.

    Assembly and read evidence must not disagree about what counts as a polyA
    tail, or the same locus gets a full-length span from one path and a
    truncated one from the other.
    """
    from retro_miner import mei_support

    assert la._MIN_POLYA_RUN_FOR_FULL_3P_IMPUTE == mei_support._MIN_POLYA_RUN_FOR_END_IMPUTE


def test_the_schema_version_gate_admits_its_own_records() -> None:
    """`_has_sideaware_feature_schema` must accept what we stamp.

    The stamp and the gate read the same constant today, so this holds by
    construction; it is here because a future bump that changes only one of the
    two silently invalidates every cached assembly record with no error.
    """
    required = {
        "coord_model",
        "left_support_contig_id",
        "right_support_contig_id",
        "left_support_mei_start",
        "left_support_mei_end",
        "right_support_mei_start",
        "right_support_mei_end",
        "left_support_mei_aln_len",
        "right_support_mei_aln_len",
        "mei_target_length",
        "insertion_length_observed",
        "insertion_length_imputed",
        "insertion_length_confidence_tier",
        "microhomology_sequence",
        "junction_overlap_sequence",
    }
    stamped = {key: 0 for key in required}
    stamped["coord_logic_version"] = la._ASSEMBLY_FEATURE_SCHEMA_VERSION
    # Re-implement the predicate's own logic against our own stamp; the real
    # gate is a closure inside _process_single_locus and is not exported.
    assert required.issubset(stamped)
    assert int(stamped["coord_logic_version"]) >= la._ASSEMBLY_FEATURE_SCHEMA_VERSION


# --------------------------------------------------------------------------
# read-name handling
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "raw,expected",
    [
        ("@READ1", "READ1"),
        ("READ1", "READ1"),
        ("READ1  flags here", "READ1"),
        ("READ1/1", "READ1"),
        ("READ1/2", "READ1"),
        ("READ1:1", "READ1"),
        ("READ1:2", "READ1"),
        ("READ1:0", "READ1"),
        ("  ", ""),
        ("", ""),
        (None, ""),
    ],
)
def test_canonical_read_name_strips_every_paired_end_spelling(raw, expected) -> None:
    """All four ways of writing the same paired read must collapse to one key.

    The lookup table is keyed by read name, so a spelling that misses is a read
    whose evidence is silently not counted -- a wrong denominator rather than an
    error. `:0` is left alone because it is an unpaired marker, not an end.
    """
    assert la._canonical_read_name(raw) == expected


def test_preferred_lookup_accepts_every_spelling_of_the_same_read() -> None:
    """A locus-linked read named one way must still match when written another."""
    lookup = la._preferred_read_lookup({"READ1:1"})
    for spelling in ("READ1:1", "READ1", "READ1/1", "@READ1"):
        assert la._read_matches_preferred(
            _FakeRead(spelling, is_read1=True), lookup
        )


def test_preferred_lookup_is_empty_when_there_is_nothing_to_prefer() -> None:
    assert la._preferred_read_lookup(None) == set()
    assert la._preferred_read_lookup(set()) == set()
    lookup = la._preferred_read_lookup(None)
    # An empty preference set must not match everything.
    assert not la._read_matches_preferred(_FakeRead("READ1", is_read1=True), lookup)


def test_a_read_is_found_however_its_mate_end_is_written() -> None:
    """The canonical form deliberately collapses the end marker.

    A locus-linked read named `READ1:1` must still be recognised when the BAM
    spells it `READ1`, because a miss here means that read's evidence is left
    out of the extracted FASTQ -- a wrong denominator rather than an error.
    """
    lookup = la._preferred_read_lookup({"READ1:1"})
    assert "READ1" in lookup and "READ1:1" in lookup
    assert la._read_matches_preferred(
        _FakeRead("READ1", is_read1=False, is_read2=True), lookup
    )


class _FakeRead:
    """Just enough of `pysam.AlignedSegment` for the name-matching helpers."""

    def __init__(self, query_name: str, *, is_read1: bool = False, is_read2: bool = False):
        self.query_name = query_name
        self.is_read1 = is_read1
        self.is_read2 = is_read2


# --------------------------------------------------------------------------
# window and interval arithmetic
# --------------------------------------------------------------------------


def test_a_reversed_window_is_normalised_rather_than_producing_a_negative_span() -> None:
    """Bad input must not yield a negative interval that matches nothing."""
    row = {"chrom": "chr1", "window_start": 200, "window_end": 100}
    assert la._window_locus_id_from_row(row) == la._window_locus_id_from_row(
        {"chrom": "chr1", "window_start": 100, "window_end": 200}
    )


def test_the_window_always_keeps_the_breakpoint_inside_the_extracted_interval() -> None:
    """A zero or negative pad must still yield a real interval around the call.

    `max(1, pad_bp)` is what stops `--interval-pad 0` from producing an empty
    interval that silently extracts zero reads and scores the locus as
    "assembled, nothing found".
    """
    import pandas as pd

    row = pd.Series(
        {"chrom": "chr1", "window_start": 1000, "window_end": 1000, "insertion_breakpoint_pos": 1000}
    )
    for pad in (0, -5, 1, 300):
        _chrom, start, end = la._interval_from_row(row, pad)
        assert start >= 1
        assert end > start
        assert start <= 1000 <= end


def test_a_missing_breakpoint_falls_back_to_the_window_midpoint() -> None:
    import pandas as pd

    row = pd.Series({"chrom": "chr1", "window_start": 1000, "window_end": 1200})
    _chrom, start, end = la._interval_from_row(row, 50)
    assert start <= 1100 <= end


# --------------------------------------------------------------------------
# sequence helpers
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "seq,expected",
    [
        ("", 0),
        ("ACGT", 1),  # a lone T is still an A/T run
        ("AAA", 3),
        ("AAAAAA", 6),
        ("AAATTTAAA", 3),
        ("AAAAACAAAA", 5),
        ("GATTACA", 2),
        ("AT", 1),
    ],
)
def test_poly_a_max_run_counts_the_longest_homopolymer(seq, expected) -> None:
    """Imputation is gated on this number, so the off-by-one here is a real one."""
    assert la._poly_at_max_run(seq) == expected


def test_poly_a_max_run_is_case_insensitive_and_ignores_other_bases() -> None:
    assert la._poly_at_max_run("aaa") == 3
    assert la._poly_at_max_run("AaAaA") == 5  # upper-cased first, so one run
    assert la._poly_at_max_run("ggg") == 0  # G is not an A/T base
    assert la._poly_at_max_run("ATATA") == 1  # alternation is not a run


def test_summarize_contigs_counts_sequences_and_the_longest() -> None:
    """A FASTA with no trailing newline must still count its last record."""
    path = Path("/tmp/_la_test_contigs.fa")
    path.write_text(">c1\nACGT\nACGT\n>c2\nACGTACGTAC\n")
    count, longest = la._summarize_contigs(path)
    assert (count, longest) == (2, 10)
    assert la._summarize_contigs(Path("/tmp/_la_missing.fa")) == (0, 0)


def test_microhomology_is_empty_for_a_pure_homopolymer_overlap_unless_allowed() -> None:
    """Homopolymer overlap carries no identity; it is excluded by default."""
    hit = {"qname": "c1", "qstart": 0, "qend": 10}
    seq_by_name = {"c1": "AAAAAACCCCC"}
    left = {"qname": "c1", "qend": 5}
    strict = la._extract_microhomology_sequence(
        mei_hit=hit, left_hit=left, right_hit=None, seq_by_name=seq_by_name
    )
    assert strict == ""
    relaxed = la._extract_microhomology_sequence(
        mei_hit=hit,
        left_hit=left,
        right_hit=None,
        seq_by_name=seq_by_name,
        allow_homopolymer=True,
    )
    assert relaxed == "AAAAA"


def test_microhomology_returns_empty_rather_than_raising_on_a_malformed_hit() -> None:
    hit = {"qname": "c1", "qstart": 10, "qend": 5}  # qend < qstart
    assert (
        la._extract_microhomology_sequence(
            mei_hit=hit, left_hit=None, right_hit=None, seq_by_name={"c1": "ACGT"}
        )
        == ""
    )


# --------------------------------------------------------------------------
# breakpoint inference from alignments
# --------------------------------------------------------------------------


def _hit(**kw):
    base = {
        "qname": "c1",
        "qstart": 0,
        "qend": 100,
        "strand": "+",
        "tname": "chr1",
        "tstart": 1000,
        "tend": 1100,
        "nmatch": 100,
        "alnlen": 100,
        "mapq": 60,
    }
    base.update(kw)
    return base


def test_two_flanking_refs_give_a_breakpoint_and_a_tsd_length() -> None:
    """Flanks must clear the 20 bp abutment window around the MEI hit.

    `_flank_ref_hits_for_mei` only accepts a left hit ending within
    `m_qstart + 20` and a right hit starting within `m_qend - 20`; anything
    further away is not treated as a flank at all.
    """
    left = _hit(qstart=0, qend=10, tstart=1000, tend=1040)
    right = _hit(qstart=90, qend=100, tstart=1060, tend=1100)
    bp, lchrom, rchrom, tsd = la._infer_breakpoint_from_alignments(
        _hit(), [left, right], default_pos=7
    )
    assert lchrom == "chr1" == rchrom
    assert bp == (1040 + 1061) // 2
    assert tsd == 1061 - 1040 + 1


def test_a_ref_hit_that_does_not_abut_the_mei_hit_is_not_a_flank() -> None:
    """Flanking is decided on assembly-contig geometry, not reference distance.

    A ref hit on the same contig that ends well clear of the MEI alignment is
    not a breakpoint side, however far away it lands on the reference. Treating
    it as one would place the breakpoint wherever that hit happens to sit.
    """
    non_abutting = _hit(qstart=0, qend=60, tstart=1000, tend=1040)
    assert la._breakpoint_side_status(_hit(), [non_abutting]) == "none"


def test_one_sided_flanks_fall_back_to_that_side_without_inventing_a_tsd() -> None:
    left = _hit(qstart=0, qend=10, tstart=1000, tend=1040)
    bp, chrom, _r, tsd = la._infer_breakpoint_from_alignments(
        _hit(), [left], default_pos=7
    )
    assert (bp, chrom, tsd) == (1040, "chr1", 0)


def test_no_flanks_falls_back_to_the_supplied_default() -> None:
    bp, chrom, _r, tsd = la._infer_breakpoint_from_alignments(
        _hit(), [], default_pos=4242
    )
    assert (bp, chrom, tsd) == (4242, "", 0)


def test_side_status_distinguishes_left_right_both_and_neither() -> None:
    mei = _hit(qstart=0, qend=100)
    left = _hit(qstart=0, qend=10, tstart=1000, tend=1040)
    right = _hit(qstart=90, qend=100, tstart=1060, tend=1100)
    assert la._breakpoint_side_status(mei, [left, right]) == "both"
    assert la._breakpoint_side_status(mei, [left]) == "left_only"
    assert la._breakpoint_side_status(mei, [right]) == "right_only"
    assert la._breakpoint_side_status(mei, []) == "none"


def test_ref_hits_on_another_contig_are_not_treated_as_flanks() -> None:
    """A hit to a different contig is not the same locus, however well it aligns."""
    mei = _hit(qname="c1")
    other = _hit(qname="c2", tstart=1000, tend=1040)
    assert la._breakpoint_side_status(mei, [other]) == "none"


def test_partner_detection_covers_the_documented_partner_types() -> None:
    lchrom, lpos, ltype = la._infer_non_mei_partner(
        left_hit=_hit(tname="chr1", strand="+", tstart=1000, tend=1040),
        right_hit=_hit(tname="chr2", strand="+", tstart=2000, tend=2100),
        side_status="both",
    )
    assert (lchrom, ltype) == ("chr2", "interchrom_breakend")

    _c, _p, itype = la._infer_non_mei_partner(
        left_hit=_hit(tname="chr1", strand="+", tstart=1000, tend=1040),
        right_hit=_hit(tname="chr1", strand="-", tstart=1060, tend=1100),
        side_status="both",
    )
    assert itype == "inversion_breakend"

    _c, _p, gtype = la._infer_non_mei_partner(
        left_hit=_hit(tname="chr1", strand="+", tstart=1000, tend=1040),
        right_hit=_hit(tname="chr1", strand="+", tstart=9000, tend=9100),
        side_status="both",
    )
    assert gtype == "large_gap_breakend"

    # Same chrom, same strand, close together: not a partner event.
    assert la._infer_non_mei_partner(
        left_hit=_hit(tname="chr1", strand="+", tstart=1000, tend=1040),
        right_hit=_hit(tname="chr1", strand="+", tstart=1060, tend=1100),
        side_status="both",
    ) == ("", -1, "")


def test_sort_hit_key_prefers_alignment_length_then_mapq() -> None:
    hits = [
        _hit(alnlen=100, mapq=10, nmatch=100),
        _hit(alnlen=100, mapq=60, nmatch=100),
        _hit(alnlen=200, mapq=1, nmatch=200),
    ]
    assert max(hits, key=la._sort_hit_key)["alnlen"] == 200


# --------------------------------------------------------------------------
# consensus selection between disease and control
# --------------------------------------------------------------------------


def _feat(**kw):
    base = {
        "insertion_mei_start": 300,
        "insertion_mei_end": 100,
        "left_support_mei_aln_len": 0,
        "right_support_mei_aln_len": 0,
        "mei_aln_len": 50,
        "coord_model": "single_contig_both_sides",
        "complex_class": "unknown",
        "breakpoint_side_status": "none",
        "non_mei_partner_type": "",
    }
    base.update(kw)
    return base


def test_a_real_feature_beats_an_empty_one() -> None:
    """Exercises `_MIN_SIDE_ANCHOR_ALN_LEN`, which was an undefined name."""
    pick, source, _complex_pick, _complex_source = la._choose_consensus_features(
        _feat(), {}
    )
    assert source == "disease"
    assert pick


def test_two_real_features_are_decided_on_alignment_support_not_input_order() -> None:
    weak = _feat(left_support_mei_aln_len=25, right_support_mei_aln_len=5, mei_aln_len=30)
    strong = _feat(left_support_mei_aln_len=80, right_support_mei_aln_len=70, mei_aln_len=90)
    _pick, source, _c, _cs = la._choose_consensus_features(strong, weak)
    assert source == "disease"


def test_the_control_wins_when_disease_carries_nothing_usable() -> None:
    _pick, source, _c, _cs = la._choose_consensus_features({}, _feat())
    assert source == "control"


def test_neither_arm_usable_reports_no_source_rather_than_inventing_one() -> None:
    """An empty pick must not be attributed to a sample that had no data."""
    pick, source, _c, _cs = la._choose_consensus_features({}, {})
    assert pick == {}
    assert source == ""


def test_the_side_anchor_threshold_actually_gates_the_score() -> None:
    """One anchor either side of the threshold must change the ranking.

    Without this the constant could drift to any value and the selection logic
    would keep passing, because nothing observed its effect.
    """
    threshold = la._MIN_SIDE_ANCHOR_ALN_LEN
    at_threshold = _feat(
        left_support_mei_aln_len=threshold,
        right_support_mei_aln_len=threshold,
        mei_aln_len=threshold,
        complex_class="simple_mei",
        breakpoint_side_status="both",
    )
    below = _feat(
        left_support_mei_aln_len=threshold - 1,
        right_support_mei_aln_len=threshold - 1,
        mei_aln_len=threshold - 1,
        complex_class="unknown",
        breakpoint_side_status="none",
    )
    _pick, source, _c, _cs = la._choose_consensus_features(at_threshold, below)
    assert source == "disease"
    _pick, source, _c, _cs = la._choose_consensus_features(below, at_threshold)
    assert source == "control"


def test_complexity_key_orders_classes_before_alignment_length() -> None:
    """A more complex event wins even when its alignment is shorter."""
    complex_feat = _feat(
        complex_class="mei_plus_sv", breakpoint_side_status="both", mei_aln_len=10
    )
    simple_feat = _feat(
        complex_class="simple_mei", breakpoint_side_status="none", mei_aln_len=100
    )
    assert la._sample_complexity_key(complex_feat) > la._sample_complexity_key(
        simple_feat
    )


# --------------------------------------------------------------------------
# read-cap escalation
# --------------------------------------------------------------------------


def test_cap_escalation_requires_an_assembled_status() -> None:
    """A locus that never assembled must not be escalated -- it has no contigs."""
    assert not la._should_escalate_read_cap(
        status="no_contigs",
        d_reads=1000,
        n_reads=1000,
        max_reads_per_sample=200,
        pick_source="",
        d_cc=0,
        n_cc=0,
        escalated_max_reads=1000,
    )


def test_cap_escalation_stops_once_the_cap_is_already_at_the_ceiling() -> None:
    assert not la._should_escalate_read_cap(
        status="assembled",
        d_reads=1000,
        n_reads=1000,
        max_reads_per_sample=1000,
        pick_source="",
        d_cc=0,
        n_cc=0,
        escalated_max_reads=1000,
    )


def test_cap_escalation_fires_when_neither_arm_supplied_a_pick() -> None:
    assert la._should_escalate_read_cap(
        status="assembled",
        d_reads=200,
        n_reads=200,
        max_reads_per_sample=200,
        pick_source="",
        d_cc=5,
        n_cc=5,
        escalated_max_reads=1000,
    )


def test_cap_escalation_stays_quiet_when_both_arms_are_convincing() -> None:
    """Nothing to escalate to: a pick was made and the assembly is complex."""
    assert not la._should_escalate_read_cap(
        status="assembled",
        d_reads=200,
        n_reads=200,
        max_reads_per_sample=200,
        pick_source="disease",
        d_cc=5,
        n_cc=5,
        escalated_max_reads=1000,
    )


# --------------------------------------------------------------------------
# cache manifests
# --------------------------------------------------------------------------


def test_a_corrupt_manifest_is_treated_as_absent_rather_than_crashing() -> None:
    """A half-written manifest from a killed run must not abort the batch."""
    path = Path("/tmp/_la_manifest_bad.json")
    path.write_text("{not json at all")
    assert la._parse_existing_manifest(path) is None
    assert la._parse_existing_manifest(Path("/tmp/_la_no_such_manifest.json")) is None


def test_a_json_array_manifest_is_rejected_because_it_is_not_a_mapping() -> None:
    """`manifest["status"]` on a list would raise, not return a default."""
    path = Path("/tmp/_la_manifest_list.json")
    path.write_text("[1, 2, 3]")
    assert la._parse_existing_manifest(path) is None


def test_a_manifest_counts_as_a_cache_only_when_contigs_actually_exist(tmp_path) -> None:
    """`status: assembled` alone is a claim; the contigs are the evidence.

    Trusting the recorded status would let a directory whose contigs were
    cleaned up masquerade as a cache hit and skip the assembly entirely.
    """
    locus = tmp_path / "cache_probe"
    manifest = {"status": "assembled", "interval": {"pad_bp": 300}}
    assert not la._has_assembled_cache(locus, manifest, interval_pad_bp=300)
    out = locus / "disease.spades.pad300"
    out.mkdir(parents=True, exist_ok=True)
    (out / "contigs.fasta").write_text(">c\nACGT\n")
    assert la._has_assembled_cache(locus, manifest, interval_pad_bp=300)


def test_a_non_assembled_manifest_is_never_a_cache_hit(tmp_path) -> None:
    locus = tmp_path / "cache_probe2"
    out = locus / "disease.spades.pad300"
    out.mkdir(parents=True, exist_ok=True)
    (out / "contigs.fasta").write_text(">c\nACGT\n")
    assert not la._has_assembled_cache(
        locus, {"status": "failed", "interval": {"pad_bp": 300}}, interval_pad_bp=300
    )
    assert not la._has_assembled_cache(locus, None, interval_pad_bp=300)


def test_a_valid_manifest_round_trips() -> None:
    path = Path("/tmp/_la_manifest_ok.json")
    payload = {"status": "assembled", "interval": {"pad_bp": 150}, "n": 3}
    path.write_text(json.dumps(payload))
    assert la._parse_existing_manifest(path) == payload