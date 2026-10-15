"""Tests for the per-call transduction detector (transduction_detect.py).

The fixture encodes the distinction the module exists to make: a call whose extra
sequence matches another locus's transducible flank (partnered), one whose extra
sequence is nowhere in the reference (orphan), one whose extra sequence is the
reference next to its own insertion (overhang), and one whose extra sequence is
just a poly(A) tail (terminated at its own signal, not a transduction).
"""

from __future__ import annotations

import gzip
import random
import struct
from pathlib import Path

import pytest

from retro_miner.transduction_detect import (
    DEFAULT_MIN_SEGMENT_BP,
    NOT_TRANSDUCTION,
    ORPHAN_CANDIDATE,
    PARTNERED,
    SOURCE_AMBIGUOUS,
    UNEVALUABLE,
    TwoBitReference,
    annotate_transduction,
    build_annotation_index,
    candidate_source_elements,
    consensus_coords,
    detect_junction_transduction,
    family_of_repclass,
    fasta_fetcher,
    find_donors,
    main,
    rmsk_rows,
    strip_poly_at_ends,
)
from retro_miner.transduction_tags import build_source_index, extract_source_flank

_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def _revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


_RNG = random.Random(20261014)


def _random_seq(n: int, rng: random.Random = _RNG) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


class _DictGenome:
    """Tiny dict-backed SequenceFetcher; coordinates 1-based inclusive."""

    def __init__(self, seqs: dict[str, str]) -> None:
        self.seqs = {k: v.upper() for k, v in seqs.items()}

    def __call__(self, chrom: str, start: int, end: int) -> str:
        return self.seqs[chrom][start - 1 : end]


# --- a minimal 2bit writer, so the reader is tested against the real format ----


def _pack_2bit(seq: str) -> bytes:
    codes = {"T": 0, "C": 1, "A": 2, "G": 3}
    out = bytearray()
    for i in range(0, len(seq), 4):
        byte = 0
        for j in range(4):
            base = seq[i + j] if i + j < len(seq) else "T"  # padding is not visible
            byte = (byte << 2) | codes[base]
        out.append(byte)
    return bytes(out)


def _write_2bit(path: Path, seqs: dict[str, str],
                n_blocks: dict[str, list[tuple[int, int]]] | None = None) -> None:
    """Write a UCSC .2bit file.  ``n_blocks`` are 0-based (start, size) pairs."""
    n_blocks = n_blocks or {}
    header_size = 16 + sum(1 + len(name) + 4 for name in seqs)
    offsets: dict[str, int] = {}
    cursor = header_size
    records: dict[str, bytes] = {}
    for name, seq in seqs.items():
        blocks = n_blocks.get(name, [])
        # Packed bases: N positions carry a T, exactly as UCSC does; the reader
        # is expected to override them from the N-block table.
        packed_seq = "".join("T" if any(s <= i < s + n for s, n in blocks) else ch
                             for i, ch in enumerate(seq))
        record = (
            struct.pack("<I", len(seq))
            + struct.pack("<I", len(blocks))
            + b"".join(struct.pack("<I", s) for s, _n in blocks)
            + b"".join(struct.pack("<I", n) for _s, n in blocks)
            + struct.pack("<I", 0)   # mask block count
            + struct.pack("<I", 0)   # reserved
            + _pack_2bit(packed_seq)
        )
        offsets[name] = cursor
        records[name] = record
        cursor += len(record)
    with path.open("wb") as handle:
        handle.write(struct.pack("<IIII", 0x1A412743, 0, len(seqs), 0))
        for name in seqs:
            handle.write(bytes([len(name)]))
            handle.write(name.encode("ascii"))
            handle.write(struct.pack("<I", offsets[name]))
        for name in seqs:
            handle.write(records[name])


# --- reference access ---------------------------------------------------------


def test_twobit_round_trip_and_bounds(tmp_path: Path):
    seq = _random_seq(137)
    path = tmp_path / "tiny.2bit"
    _write_2bit(path, {"chr1": seq})
    ref = TwoBitReference(path)
    assert ref.chroms() == ("chr1",)
    assert ref.fetch("chr1", 1, len(seq)) == seq
    assert ref.fetch("chr1", 5, 20) == seq[4:20]
    # A request past the end clamps to the sequence, and a backwards request is empty.
    assert ref.fetch("chr1", len(seq) - 4, len(seq) + 500) == seq[-5:]
    assert ref.fetch("chr1", 50, 20) == ""
    assert ref.fetch("nope", 1, 10) == ""


def test_twobit_accepts_bare_and_chr_prefixed_names(tmp_path: Path):
    seq = _random_seq(60)
    path = tmp_path / "tiny.2bit"
    _write_2bit(path, {"chr7": seq})
    ref = TwoBitReference(path)
    assert ref.fetch("7", 1, 10) == seq[:10]
    assert ref.fetch("chr7", 1, 10) == seq[:10]
    # And the reverse direction, when the file itself uses the bare name.
    path2 = tmp_path / "bare.2bit"
    _write_2bit(path2, {"7": seq})
    assert TwoBitReference(path2).fetch("chr7", 1, 10) == seq[:10]


def test_twobit_n_blocks_override_packed_bases(tmp_path: Path):
    seq = _random_seq(80)
    with_n = seq[:20] + "N" * 12 + seq[32:]
    path = tmp_path / "nblocks.2bit"
    _write_2bit(path, {"chr1": with_n}, n_blocks={"chr1": [(20, 12)]})
    ref = TwoBitReference(path)
    assert ref.fetch("chr1", 1, len(with_n)) == with_n
    assert ref.fetch("chr1", 10, 40) == with_n[9:40]
    # A window straddling the block keeps its non-N bases.
    assert ref.fetch("chr1", 30, 34) == with_n[29:34]


def test_twobit_rejects_a_file_that_is_not_2bit(tmp_path: Path):
    path = tmp_path / "not.2bit"
    path.write_bytes(b"this is not a 2bit file at all")
    with pytest.raises(ValueError):
        TwoBitReference(path).fetch("chr1", 1, 10)


def test_fasta_fetcher_adapter(tmp_path: Path):
    seq = _random_seq(70)
    path = tmp_path / "tiny.fa"
    path.write_text(f">chr1\n{seq}\n")
    fetch = fasta_fetcher(path)
    assert fetch("chr1", 1, 10) == seq[:10]          # type: ignore[operator]
    assert fetch("chr1", 60, 70) == seq[59:70]       # type: ignore[operator]


# --- annotation parsing -------------------------------------------------------

_HEADER = ("bin", "swScore", "milliDiv", "milliDel", "milliIns", "genoName", "genoStart",
           "genoEnd", "genoLeft", "strand", "repName", "repClass", "repFamily", "repStart",
           "repEnd", "repLeft", "id")

_L1MC5A = "585\t484\t251\t132\t0\tchr1\t11504\t11675\t-248944747\t-\tL1MC5a\tLINE\tL1\t-2382\t395\t199\t3"


def test_rmsk_rows_handles_both_column_arities(tmp_path: Path):
    with_bin = tmp_path / "rmsk.tsv"
    with_bin.write_text("\t".join(_HEADER) + "\n" + _L1MC5A + "\n")
    rows = list(rmsk_rows(with_bin))
    assert len(rows) == 1 and rows[0]["repName"] == "L1MC5a" and rows[0]["genoStart"] == "11504"

    without_bin = tmp_path / "no_bin.tsv.gz"
    with gzip.open(without_bin, "wt") as handle:
        handle.write("\t".join(_HEADER[1:]) + "\n" + "\t".join(_L1MC5A.split("\t")[1:]) + "\n")
    rows = list(rmsk_rows(without_bin))
    assert len(rows) == 1 and rows[0]["repFamily"] == "L1" and "bin" not in rows[0]


def test_family_of_repclass_maps_only_the_two_transducible_families():
    assert family_of_repclass("LINE", "L1") == "LINE1"
    assert family_of_repclass("Retroposon", "SVA") == "SVA"
    assert family_of_repclass("SINE", "Alu") == ""
    assert family_of_repclass("LINE", "L2") == ""


def test_consensus_coords_plus_row():
    row = dict(zip(_HEADER, ("585", "0", "0", "0", "0", "chr1", "100", "6100", "-5", "+",
                             "L1HS", "LINE", "L1", "1", "6000", "0", "1")))
    copy_start, copy_end, model_len = consensus_coords(row)
    assert (copy_start, copy_end, model_len) == (1, 6000, 6000)


def test_consensus_coords_minus_row_swaps_start_and_remainder():
    # The real L1MC5a '-' row: repStart carries the remainder (-2382) and repLeft
    # the consensus start (199), so reading it with the '+' formula would inflate
    # the model to 395 + 199.
    row = dict(zip(_HEADER, _L1MC5A.split("\t")))
    copy_start, copy_end, model_len = consensus_coords(row)
    assert (copy_start, copy_end, model_len) == (199, 395, 2777)


def _rmsk_row(**overrides: str) -> dict[str, str]:
    row = dict(zip(_HEADER, ("585", "0", "0", "0", "0", "chr1", "200", "6200", "-5", "+",
                             "L1HS", "LINE", "L1", "1", "6000", "0", "1")))
    row.update(overrides)
    return row


# --- candidate source elements -------------------------------------------------


def test_candidate_source_accepts_a_full_length_l1():
    genome = _DictGenome({"chr1": _random_seq(200) + _random_seq(6000) + _random_seq(900)})
    elements = candidate_source_elements([_rmsk_row()], genome)
    assert len(elements) == 1
    assert elements[0].family == "LINE1" and elements[0].flank_dir == "3p"
    assert elements[0].flank_start0 == 6200


def test_candidate_source_rejects_a_fragmentary_copy():
    # 4 kb of a 6 kb model: not a competent source element.
    genome = _DictGenome({"chr1": _random_seq(200) + _random_seq(4000) + _random_seq(900)})
    row = _rmsk_row(genoEnd="4200", repEnd="4000", repLeft="2000")
    assert candidate_source_elements([row], genome) == []


def test_candidate_source_rejects_a_copy_missing_the_transducible_end():
    # Covers 96.7% of the model but stops 200 bp short of its 3' end.
    genome = _DictGenome({"chr1": _random_seq(200) + _random_seq(5800) + _random_seq(900)})
    row = _rmsk_row(genoEnd="6000", repEnd="5800", repLeft="200")
    assert candidate_source_elements([row], genome) == []


def test_candidate_source_sva_needs_its_five_prime_end():
    genome = _DictGenome({"chr1": _random_seq(200) + _random_seq(1600) + _random_seq(600)})
    full = _rmsk_row(repName="SVA_F", repClass="Retroposon", repFamily="SVA",
                     repStart="1", repEnd="1600", repLeft="0", genoEnd="1800")
    assert candidate_source_elements([full], genome) != []
    truncated = _rmsk_row(repName="SVA_F", repClass="Retroposon", repFamily="SVA",
                          repStart="300", repEnd="1600", repLeft="0", genoEnd="1800")
    assert candidate_source_elements([truncated], genome) == []


def test_candidate_source_skips_other_families():
    genome = _DictGenome({"chr1": _random_seq(900)})
    alu = _rmsk_row(repName="AluY", repClass="SINE", repFamily="Alu")
    assert candidate_source_elements([alu], genome) == []


def test_candidate_source_drops_a_repetitive_flank():
    genome = _DictGenome({"chr1": "AT" * 3400})
    assert candidate_source_elements([_rmsk_row(genoEnd="6200")], genome) == []


def test_build_annotation_index_reports_what_went_in(tmp_path: Path):
    genome = _DictGenome({"chr1": _random_seq(200) + _random_seq(6000) + _random_seq(900)})
    path = tmp_path / "rmsk.tsv"
    path.write_text("\t".join(_HEADER) + "\n"
                    + "\t".join(_rmsk_row().values()) + "\n")
    index, elements = build_annotation_index(path, genome)
    assert len(elements) == 1 and len(index) == 1


# --- tail stripping -------------------------------------------------------------


def test_strip_poly_at_ends_cases():
    assert strip_poly_at_ends("") == ("", 0, "")
    unchanged = "ACGT" * 10
    assert strip_poly_at_ends(unchanged) == (unchanged, 0, "")
    core, run, sides = strip_poly_at_ends("A" * 10 + "G" * 30)
    assert (core, run, sides) == ("G" * 30, 10, "5p")
    core, run, sides = strip_poly_at_ends("T" * 12 + "G" * 20 + "A" * 12)
    assert (core, run, sides) == ("G" * 20, 12, "5p,3p")
    # A run shorter than the floor is left alone: no evidence to move the boundary.
    short = "AA" + "G" * 30
    assert strip_poly_at_ends(short) == (short, 0, "")
    # min_run is respected.
    assert strip_poly_at_ends("A" * 10 + "G" * 30, min_run=12)[1] == 0


# --- donor search ----------------------------------------------------------------


def _clean_piece(flank: str, length: int = 200) -> str:
    """A flank window starting and ending on G/C, so no incidental A/T run blurs
    where the tail begins -- otherwise a trailing run of the piece merges with the
    poly(A) tail being tested and the assertion measures the fixture, not the code."""
    for start in range(0, len(flank) - length):
        window = flank[start : start + length]
        if window[0] in "CG" and window[-1] in "CG":
            return window
    raise AssertionError("fixture flank has no clean window")


def _l1_fixture(rng: random.Random) -> tuple[object, object, str]:
    """Return (genome, source element, a unique 200 bp piece of its 3' flank)."""
    body = _random_seq(1000, rng)
    downstream = _random_seq(1200, rng)
    genome = _DictGenome({"chrL": _random_seq(200, rng) + body + downstream})
    source = extract_source_flank("L1SRC", "LINE1", "chrL", "+", 200, 1200, genome)
    assert source is not None
    return genome, source, _clean_piece(source.flank_seq)


def test_find_donors_names_the_source_and_locates_the_match():
    genome, source, piece = _l1_fixture(random.Random(1))
    index = build_source_index([source])
    ranked = find_donors(piece, "LINE1", index)
    assert len(ranked) == 1
    best = ranked[0]
    assert best.source_id == "L1SRC" and best.orientation == "sense"
    assert best.matched_bp == len(piece)
    # The matched region is located inside the donor's flank, not merely named.
    offset = source.flank_seq.find(piece)
    assert offset >= 0
    assert best.donor_start0 == source.flank_start0 + offset
    assert best.donor_end0 == source.flank_start0 + offset + len(piece)


def test_find_donors_reports_the_reverse_orientation():
    genome, source, piece = _l1_fixture(random.Random(2))
    index = build_source_index([source])
    ranked = find_donors(_revcomp(piece), "LINE1", index)
    assert len(ranked) == 1 and ranked[0].orientation == "antisense"


# --- the detector ----------------------------------------------------------------


def test_partnered_transduction_with_tail_is_named():
    genome, source, piece = _l1_fixture(random.Random(3))
    index = build_source_index([source])
    call = detect_junction_transduction("c1", "LINE1", piece + "A" * 25, index,
                                        observed_excess_bp=400)
    assert call.outcome == PARTNERED and call.is_transduction
    assert call.donor is not None and call.donor.source_id == "L1SRC"
    assert call.donor.chrom == "chrL" and call.donor.strand == "+"
    assert call.donor.flank_dir == "3p"
    assert call.tail_run_bp == 25 and call.tail_side == "3p"
    assert call.transduced_bp == len(piece)
    assert call.unexplained_excess_bp == 400 - len(piece)
    assert call.donor.matched_bp == len(piece)


def test_partnered_transduction_detected_when_stored_in_reverse_orientation():
    genome, source, piece = _l1_fixture(random.Random(4))
    index = build_source_index([source])
    call = detect_junction_transduction("c2", "LINE1", "T" * 20 + _revcomp(piece), index)
    assert call.outcome == PARTNERED
    assert call.donor is not None and call.donor.orientation == "antisense"
    assert call.tail_side == "5p" and call.tail_run_bp == 20


def test_minus_strand_donor_is_named_with_its_strand():
    rng = random.Random(5)
    upstream = _random_seq(600, rng)
    genome = _DictGenome({"chrM": upstream + _random_seq(1000, rng) + _random_seq(300, rng)})
    source = extract_source_flank("L1MINUS", "LINE1", "chrM", "-", 600, 1600, genome)
    assert source is not None
    index = build_source_index([source])
    call = detect_junction_transduction("c3", "LINE1", source.flank_seq[50:250], index)
    assert call.outcome == PARTNERED and call.donor is not None
    assert call.donor.strand == "-" and call.donor.orientation == "sense"


def test_pure_poly_a_extension_is_not_a_transduction():
    genome, source, _piece = _l1_fixture(random.Random(6))
    index = build_source_index([source])
    for seq in ("A" * 40, "T" * 40):
        call = detect_junction_transduction("c", "LINE1", seq, index)
        assert call.outcome == NOT_TRANSDUCTION
        assert call.reason == "polyA_extension_only"
        assert call.transduced_bp == 0 and call.tail_run_bp == 40


def test_low_complexity_segment_is_rejected():
    genome, source, _piece = _l1_fixture(random.Random(7))
    index = build_source_index([source])
    call = detect_junction_transduction("c", "LINE1", "G" * 60, index)
    assert call.outcome == NOT_TRANSDUCTION and call.reason == "low_complexity_segment"


def test_repeat_like_segment_is_rejected():
    genome, source, _piece = _l1_fixture(random.Random(8))
    index = build_source_index([source])
    call = detect_junction_transduction("c", "LINE1", "AC" * 60, index)
    assert call.outcome == NOT_TRANSDUCTION and call.reason == "repeat_like_segment"


def test_segment_below_the_length_floor_is_unevaluable():
    genome, source, _piece = _l1_fixture(random.Random(9))
    index = build_source_index([source])
    short = _random_seq(DEFAULT_MIN_SEGMENT_BP - 5)
    call = detect_junction_transduction("c", "LINE1", short, index)
    assert call.outcome == UNEVALUABLE and call.reason == "segment_below_floor"
    with_tail = detect_junction_transduction("c", "LINE1", short + "A" * 20, index)
    assert with_tail.outcome == UNEVALUABLE
    assert with_tail.reason == "segment_below_floor_after_tail_strip"
    assert with_tail.segment_bp == len(short)


def test_missing_inputs_are_unevaluable_never_negative():
    genome, source, piece = _l1_fixture(random.Random(10))
    index = build_source_index([source])
    assert detect_junction_transduction("c", "LINE1", "", index).reason == "no_junction_sequence"
    assert detect_junction_transduction("c", "LINE1", piece, None).reason == "no_source_index"
    empty_index = build_source_index([])
    assert detect_junction_transduction("c", "LINE1", piece, empty_index).reason == "no_source_index"


def test_segment_matching_the_local_reference_is_an_overhang():
    genome, source, _piece = _l1_fixture(random.Random(11))
    index = build_source_index([source])
    local = _random_seq(400, random.Random(12))
    segment = local[100:260]      # the reference continuation, not a transduction
    call = detect_junction_transduction("c", "LINE1", segment, index,
                                        local_reference_seq=local)
    assert call.outcome == NOT_TRANSDUCTION
    assert call.reason == "overhang_matches_local_reference"
    assert call.transduced_bp == 0


def test_segment_absent_from_the_reference_is_an_orphan_candidate():
    genome, source, _piece = _l1_fixture(random.Random(13))
    index = build_source_index([source])
    call = detect_junction_transduction("c", "LINE1", _random_seq(150, random.Random(14)),
                                        index, local_reference_seq=_random_seq(400))
    assert call.outcome == ORPHAN_CANDIDATE and call.reason == "no_reference_match"
    # A candidate is not a named source, and is not reported as a transduction.
    assert not call.is_transduction and call.donor is None


def test_equal_donors_tie_and_no_source_is_named():
    rng = random.Random(15)
    flank_a = _random_seq(500, rng)
    flank_b = _random_seq(500, rng)
    genome = _DictGenome({"chrA": _random_seq(50, rng) + _random_seq(400, rng) + flank_a,
                          "chrB": _random_seq(50, rng) + _random_seq(400, rng) + flank_b})
    src_a = extract_source_flank("SA", "LINE1", "chrA", "+", 50, 450, genome)
    src_b = extract_source_flank("SB", "LINE1", "chrB", "+", 50, 450, genome)
    assert src_a is not None and src_b is not None
    index = build_source_index([src_a, src_b])
    # A poly-C spacer contributes seeds to neither flank, so the two donors really
    # do carry equal evidence -- a random spacer can straddle into one flank's
    # continuation and silently break the tie (it did, before this comment).
    segment = flank_a[100:200] + "C" * 30 + flank_b[100:200]
    call = detect_junction_transduction("c", "LINE1", segment, index)
    assert call.outcome == SOURCE_AMBIGUOUS and call.reason == "donor_tie"
    assert call.donor is None and call.n_donor_candidates == 2
    # Ambiguity is reported as a transduction with no named source, not as absence.
    assert call.is_transduction


def test_a_donor_of_another_family_is_never_named():
    rng = random.Random(16)
    upstream = _random_seq(900, rng)
    genome = _DictGenome({"chrS": upstream + _random_seq(1600, rng) + _random_seq(200, rng)})
    sva = extract_source_flank("SVA_SRC", "SVA", "chrS", "+", 900, 2500, genome)
    assert sva is not None
    index = build_source_index([sva])
    call = detect_junction_transduction("c", "LINE1", sva.flank_seq[50:250], index)
    assert call.outcome == ORPHAN_CANDIDATE and call.donor is None
    # The same segment tagged as SVA does name the SVA source.
    sva_call = detect_junction_transduction("c2", "SVA", sva.flank_seq[50:250], index)
    assert sva_call.outcome == PARTNERED and sva_call.donor is not None
    assert sva_call.donor.family == "SVA" and sva_call.donor.flank_dir == "5p"


def test_observed_excess_is_recorded_even_without_a_donor():
    genome, source, _piece = _l1_fixture(random.Random(17))
    index = build_source_index([source])
    call = detect_junction_transduction("c", "LINE1", _random_seq(120, random.Random(18)),
                                        index, observed_excess_bp=300)
    assert call.observed_excess_bp == 300
    # An orphan still measures: the segment itself is the only candidate, so the
    # excess beyond it stays visible rather than being silently attributed.
    assert call.transduced_bp == 120
    assert call.unexplained_excess_bp == 180


# --- result shape ----------------------------------------------------------------


def test_result_is_frozen_and_serialises_blank_donor_fields():
    genome, source, piece = _l1_fixture(random.Random(19))
    index = build_source_index([source])
    call = detect_junction_transduction("c", "LINE1", piece, index)
    with pytest.raises(Exception):
        call.outcome = "changed"  # type: ignore[misc]
    row = call.to_dict()
    for key in ("call_id", "family", "outcome", "reason", "junction_bp", "segment_bp",
                "tail_run_bp", "tail_side", "transduced_bp", "observed_excess_bp",
                "unexplained_excess_bp", "n_donor_candidates", "donor_source_id",
                "donor_chrom", "donor_strand", "donor_flank_dir",
                "donor_matched_start0", "donor_matched_end0", "donor_orientation",
                "donor_seed_sites", "donor_seed_hits"):
        assert key in row
    assert row["donor_source_id"] == "L1SRC"
    orphan = detect_junction_transduction("c", "LINE1", _random_seq(120, random.Random(20)),
                                          index)
    blank = orphan.to_dict()
    assert blank["donor_source_id"] == "" and blank["donor_matched_start0"] == ""


def test_annotate_transduction_preserves_order_and_tolerates_missing_fields():
    genome, source, piece = _l1_fixture(random.Random(21))
    index = build_source_index([source])
    rows = [
        {"call_id": "a", "family": "LINE1", "junction_seq": piece + "A" * 20,
         "observed_excess_bp": "512"},
        {"call_id": "b", "family": "LINE1", "junction_seq": ""},
        {"call_id": "c", "family": "LINE1", "junction_seq": piece,
         "observed_excess_bp": 180.0, "local_reference_seq": "ACGT" * 50},
        {"call_id": "d", "family": "LINE1", "junction_seq": piece, "observed_excess_bp": ""},
    ]
    out = annotate_transduction(rows, index)
    assert [call.call_id for call in out] == ["a", "b", "c", "d"]
    assert out[0].observed_excess_bp == 512
    assert out[1].outcome == UNEVALUABLE
    assert out[2].observed_excess_bp == 180
    assert out[3].observed_excess_bp is None
    assert all(call.outcome == PARTNERED for call in (out[0], out[2], out[3]))


# --- CLI -----------------------------------------------------------------------


def test_cli_builds_the_index_and_annotates_a_calls_table(tmp_path: Path):
    rng = random.Random(22)
    body = _random_seq(6000, rng)
    downstream = _random_seq(900, rng)
    genome = _random_seq(200, rng) + body + downstream
    ref_path = tmp_path / "tiny.2bit"
    _write_2bit(ref_path, {"chr1": genome})

    rmsk_path = tmp_path / "rmsk.tsv"
    rmsk_path.write_text("\t".join(_HEADER) + "\n" + "\t".join(_rmsk_row().values()) + "\n")

    junction = _clean_piece(downstream) + "A" * 20
    calls_path = tmp_path / "calls.tsv"
    with calls_path.open("w") as handle:
        handle.write("call_id\tfamily\tjunction_seq\tobserved_excess_bp\n")
        handle.write(f"callA\tLINE1\t{junction}\t500\n")
        handle.write(f"callB\tLINE1\t{_random_seq(150, random.Random(23))}\t\n")

    out_path = tmp_path / "out.tsv"
    assert main(["--rmsk", str(rmsk_path), "--reference", str(ref_path),
                 "--calls", str(calls_path), "--out", str(out_path)]) == 0
    import pandas as pd

    result = pd.read_csv(out_path, sep="\t")
    assert list(result["call_id"]) == ["callA", "callB"]
    assert list(result["outcome"]) == [PARTNERED, ORPHAN_CANDIDATE]
    assert result.loc[0, "donor_chrom"] == "chr1"
    assert int(result.loc[0, "transduced_bp"]) == 200
    assert int(result.loc[0, "unexplained_excess_bp"]) == 300
