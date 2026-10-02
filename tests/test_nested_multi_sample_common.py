"""Tests for the shared multi-sample nested-insertion module.

The recurring theme is that this module's job is to refuse to smooth over
awkwardness: a four-value enum that must not become three, a contig filter that
must not widen, a FORMAT column that must never be read, and a multi-carrier
field whose separator, if guessed wrong, quietly reduces a five-genome cohort to
its private sites and looks like a result.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import nested_multi_sample_common as common  # noqa: E402


# --------------------------------------------------------------------------
# the four-value NESTED enum
# --------------------------------------------------------------------------


def test_all_four_canonical_values_survive_a_round_trip():
    for value in ("unnested", "nested_sense", "nested_antisense", "nested_unknown"):
        assert common.parse_nested_state(value) == value


def test_nested_unknown_is_never_folded_into_antisense():
    assert common.parse_nested_state("nested_unknown") != "nested_antisense"
    assert common.parse_nested_state("nested_unknown") == "nested_unknown"


def test_the_legacy_coarse_value_gets_its_own_label_not_an_orientation():
    assert common.parse_nested_state("nested") == common.LEGACY_NESTED_UNLABELED
    assert common.parse_nested_state("nested") not in ("nested_sense", "nested_antisense")


def test_an_unrecognised_nested_value_raises_rather_than_coercing():
    with pytest.raises(common.NestedEnumError):
        common.parse_nested_state("nested_maybe")
    with pytest.raises(common.NestedEnumError):
        common.parse_nested_state("")


def test_the_census_reports_every_value_including_the_empty_ones():
    census = common.nested_census(
        ["unnested", "nested_sense", "nested_antisense", "nested_unknown", "nested"]
    )
    assert census == {
        "unnested": 1,
        "nested_sense": 1,
        "nested_antisense": 1,
        "nested_unknown": 1,
        common.LEGACY_NESTED_UNLABELED: 1,
        "unrecognised_rejected": 0,
    }


def test_the_census_normalises_a_raw_legacy_value_rather_than_dropping_it():
    census = common.nested_census(["nested", "nested", "unnested"])
    assert census[common.LEGACY_NESTED_UNLABELED] == 2
    assert census["unrecognised_rejected"] == 0


def test_the_census_counts_a_genuinely_unrecognised_value_separately():
    census = common.nested_census(["nested_maybe", "nested"])
    assert census["unrecognised_rejected"] == 1
    assert census[common.LEGACY_NESTED_UNLABELED] == 1


# --------------------------------------------------------------------------
# contig scope
# --------------------------------------------------------------------------


@pytest.mark.parametrize("chrom", ["chr1", "chr9", "chr22", "chrX"])
def test_in_scope_contigs_are_accepted(chrom):
    assert common.in_scope_contig(chrom) is True


@pytest.mark.parametrize("chrom", ["chrY", "chrM", "chr23", "1", "", "chrX_random"])
def test_out_of_scope_contigs_are_rejected(chrom):
    assert common.in_scope_contig(chrom) is False


# --------------------------------------------------------------------------
# subfamily canonicalisation
# --------------------------------------------------------------------------


def test_the_class_suffix_is_stripped_so_subfamilies_can_be_compared():
    assert common.canonical_subfamily("AluYb9#SINE/Alu") == "AluYb9"
    assert common.canonical_subfamily("AluYb8#SINE/Alu") == "AluYb8"


def test_subfamilies_differing_only_in_class_path_compare_equal():
    assert common.canonical_subfamily("AluYb9#SINE/Alu") == common.canonical_subfamily(
        "AluYb9#SINE/Other"
    )


def test_a_missing_subfamily_is_none_not_an_empty_string():
    assert common.canonical_subfamily("") is None
    assert common.canonical_subfamily(None) is None


# --------------------------------------------------------------------------
# GT / GQ must never enter
# --------------------------------------------------------------------------


VCF_WITH_GENOTYPES = """##fileformat=VCFv4.2
##INFO=<ID=NESTED,Number=1,Type=String,Description="nested state">
##INFO=<ID=ORIENT,Number=1,Type=String,Description="orientation">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
chr1\t100\tid1\tN\t<INS>\t.\t.\tNESTED=nested;ORIENT=+\tGT:GQ\t1/1:60\t0/0:30
"""


def test_parsing_a_callset_never_reads_the_format_column(tmp_path):
    path = tmp_path / "s.vcf"
    path.write_text(VCF_WITH_GENOTYPES)
    parsed = common.read_callset(path)
    call = parsed["chr1"][0]
    assert call["nested_state"] == common.LEGACY_NESTED_UNLABELED
    assert call["orientation"] == "+"
    # No genotype-shaped key may exist on a parsed call at all.
    assert "GT" not in call
    assert "GQ" not in call
    assert not any(k in ("GT", "GQ", "gt", "gq") for k in call)
    assert "1/1" not in repr(call)
    assert "60" not in repr(call)


def test_the_read_callset_signature_only_pulls_the_columns_it_needs():
    """A structural guard: adding FORMAT parsing here would change the signature."""
    import inspect

    source = inspect.getsource(common.read_callset)
    assert "parts[9]" not in source
    assert "FORMAT" not in source.replace("FORMAT column", "")


# --------------------------------------------------------------------------
# loading the dedup output
# --------------------------------------------------------------------------


UNIQUE_SITES_HEADER = [
    "site_id", "chrom", "representative_pos", "representative_sample", "family",
    "insertion_orientation", "host_name", "host_start0", "host_end0", "host_strand",
    "host_len", "host_offset_5p_0based", "host_selection_rule", "n_carriers",
    "samples", "private", "private_to_sample", "window_bp",
    "tsd_overlap_within_site", "site_allele_count", "site_allele_freq",
    "pooled_site_matched",
]


def _write_unique_sites(path: Path, rows) -> None:
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=UNIQUE_SITES_HEADER)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _site_row(**overrides) -> dict:
    base = {
        "site_id": "US00001", "chrom": "chr1", "representative_pos": "1000",
        "representative_sample": "S1", "family": "ALU", "insertion_orientation": "+",
        "host_name": "AluSx1", "host_start0": "900", "host_end0": "1200",
        "host_strand": "+", "host_len": "300", "host_offset_5p_0based": "100",
        "host_selection_rule": "longest_span", "n_carriers": "1", "samples": "S1",
        "private": "True", "private_to_sample": "S1", "window_bp": "10",
        "tsd_overlap_within_site": "0", "site_allele_count": "1",
        "site_allele_freq": "0.2", "pooled_site_matched": "False",
    }
    base.update(overrides)
    return base


def test_multi_carrier_sites_are_parsed_from_a_pipe_separated_field(tmp_path):
    """The separator is the whole ballgame.

    Splitting on a comma or semicolon yields one bogus token per multi-carrier
    site, which silently drops every shared site from the per-call join and
    leaves only private sites -- a plausible-looking result rather than an error.
    """
    path = tmp_path / "us.csv"
    _write_unique_sites(
        path,
        [
            _site_row(site_id="US00001", samples="S1|S2|S3", n_carriers="3", private="False"),
            _site_row(site_id="US00002", samples="S1", n_carriers="1", private="True"),
        ],
    )
    rows, _ = common.load_unique_sites(path)
    assert rows[0]["carriers"] == ["S1", "S2", "S3"]
    assert rows[0]["n_carriers"] == 3
    assert rows[0]["private"] is False
    assert rows[1]["carriers"] == ["S1"]


def test_out_of_scope_contigs_are_dropped_and_counted(tmp_path):
    path = tmp_path / "us.csv"
    _write_unique_sites(
        path,
        [_site_row(site_id="a", chrom="chr1"), _site_row(site_id="b", chrom="chrY")],
    )
    rows, report = common.load_unique_sites(path)
    assert len(rows) == 1
    assert report["dropped_out_of_scope_contig"] == 1
    assert report["rows_in_scope"] == 1


def test_rows_missing_geometry_are_dropped_and_counted(tmp_path):
    path = tmp_path / "us.csv"
    _write_unique_sites(
        path,
        [
            _site_row(site_id="a"),
            _site_row(site_id="b", host_len=""),
            _site_row(site_id="c", representative_pos=""),
        ],
    )
    rows, report = common.load_unique_sites(path)
    assert len(rows) == 1
    assert report["dropped_missing_geometry"] == 2


def test_typed_fields_are_converted(tmp_path):
    path = tmp_path / "us.csv"
    _write_unique_sites(path, [_site_row()])
    rows, _ = common.load_unique_sites(path)
    assert rows[0]["pos"] == 1000
    assert rows[0]["host_len"] == 300
    assert rows[0]["host_offset"] == 100


# --------------------------------------------------------------------------
# schema drift is an error, not an empty cohort
# --------------------------------------------------------------------------


def test_a_table_missing_the_host_columns_is_rejected(tmp_path):
    """The wrong `unique_sites.csv` must not load.

    Two tables in this workspace share the filename. The multi-sample site table
    has `same_family_host_name` and no host interval at all. Fed to this loader
    it used to fail the geometry gate on every row, so the loader returned 0
    sites from 4,998 rows and both Phase 4 analyses carried on to report an
    empty cohort as a null result. The header below is that file's real header.
    """
    path = tmp_path / "us.csv"
    other_header = [
        "site_id", "chrom", "representative_pos", "family",
        "insertion_orientation", "same_family_nested_state", "raw_NESTED_value",
        "same_family_host_id", "same_family_host_name", "same_family_host_strand",
        "same_family_host_length", "any_MEI_host_id_for_matching",
        "host_offset_5p_0based", "n_carriers", "carrier_count_class", "samples",
        "private", "private_to_sample", "source_record_count", "source_call_ids",
        "match_window_bp",
    ]
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=other_header)
        writer.writeheader()
        writer.writerow({name: "x" for name in other_header})

    with pytest.raises(common.SchemaError) as excinfo:
        common.load_unique_sites(path)
    message = str(excinfo.value)
    for column in ("host_name", "host_start0", "host_end0", "host_len"):
        assert column in message


def test_a_file_where_every_row_fails_the_gate_loads_empty_without_raising(tmp_path):
    """The loader stays permissive; the analysis scripts are what refuse.

    A table of hostless sites legitimately has empty geometry on every row, so
    "all rows dropped" is not by itself an error at parse time. The empty-cohort
    guard lives in `joint_enrichment.py` and `recurrence_test.py`, which is
    where "no effect" would actually be printed.
    """
    path = tmp_path / "us.csv"
    _write_unique_sites(
        path,
        [_site_row(site_id="a", chrom="chrY"), _site_row(site_id="b", chrom="chrY")],
    )
    rows, report = common.load_unique_sites(path)
    assert rows == []
    assert report["rows_read"] == 2
    assert report["dropped_out_of_scope_contig"] == 2


def test_a_header_only_file_is_empty_rather_than_an_error(tmp_path):
    """No rows at all is a legitimately empty input; all rows failing is not."""
    path = tmp_path / "us.csv"
    _write_unique_sites(path, [])
    rows, report = common.load_unique_sites(path)
    assert rows == []
    assert report["rows_read"] == 0


def test_every_required_column_is_one_the_loader_actually_reads():
    """Keep the guard list honest about the loader's own fields.

    The constant is what decides whether a schema change is caught, so a required
    column the fixture header lacks -- or a load-bearing column dropped from the
    list -- would reintroduce the silent failure the guard exists to prevent.
    """
    required = set(common.REQUIRED_UNIQUE_SITE_COLUMNS)
    assert required <= set(UNIQUE_SITES_HEADER)
    for column in ("chrom", "representative_pos", "host_len",
                   "host_offset_5p_0based", "n_carriers", "samples"):
        assert column in required


# --------------------------------------------------------------------------
# the join
# --------------------------------------------------------------------------


def _call(pos: int, orientation: str = "+", subfamily: str = "AluYb9",
          tsd: str = "ACGTACGT") -> dict:
    return {
        "pos": pos,
        "id": "x",
        "family": "ALU",
        "subfamily": subfamily,
        "orientation": orientation,
        "nested_raw": "nested",
        "nested_state": common.LEGACY_NESTED_UNLABELED,
        "tsd": tsd,
        "tsd_len": len(tsd),
        "tsd_reported": bool(tsd),
        "mei_span": 282,
        "known_mei": "",
        "call_tier": "",
    }


def _attach_one(site, calls):
    return common.attach_call_details([typed(site)], {"S1": {"chr1": calls}})[0]


def typed(row, tmp_path=None) -> dict:
    """A site as `load_unique_sites` produces it, from a CSV-shaped row.

    The join operates on loaded sites, not on raw CSV rows, so the tests build
    them through the real loader rather than hand-rolling the dict -- otherwise
    the tests would pass against a shape the pipeline never produces.
    """
    import tempfile

    holder = tmp_path
    if holder is None:
        holder = Path(tempfile.mkdtemp())
    path = holder / "us.csv"
    _write_unique_sites(path, [row])
    rows, _ = common.load_unique_sites(path)
    return rows[0]


def test_nearest_call_within_the_window_is_used(tmp_path):
    site = _site_row(representative_pos="1000", n_carriers="1", samples="S1")
    attached = _attach_one(site, [_call(990), _call(1004)])
    assert len(attached["calls"]) == 1
    assert attached["calls"][0]["pos"] == 1004


def test_an_exact_tie_keeps_the_first_call_deterministically(tmp_path):
    """Ties must resolve deterministically, not by iteration luck."""
    site = _site_row(representative_pos="1000", n_carriers="1", samples="S1")
    first = _attach_one(site, [_call(995), _call(1005)])
    second = _attach_one(site, [_call(995), _call(1005)])
    assert first["calls"][0]["pos"] == second["calls"][0]["pos"]


def test_a_call_outside_the_window_is_not_joined(tmp_path):
    site = _site_row(representative_pos="1000", n_carriers="1", samples="S1")
    attached = _attach_one(site, [_call(1020)])
    assert attached["calls"] == []
    assert attached["unmatched_carrier_samples"] == ["S1"]


def test_a_consistent_join_verifies(tmp_path):
    site = _site_row(representative_pos="1000", n_carriers="1", samples="S1",
                     insertion_orientation="+")
    attached = _attach_one(site, [_call(1000, orientation="+")])
    report = common.verify_join([attached])
    assert report["verdict"] == "join_consistent_with_dedup_output"
    assert report["calls_recovered"] == 1


def test_a_carrier_count_disagreement_fails_the_gate(tmp_path):
    """The join is a gate: a partial join must stop the analysis, not warn."""
    site = _site_row(representative_pos="1000", n_carriers="2", samples="S1")
    attached = _attach_one(site, [_call(1000)])
    report = common.verify_join([attached])
    assert report["verdict"] == "join_disagrees_with_dedup_output"
    assert report["sites_where_carrier_count_disagrees"] == 1


def test_an_orientation_disagreement_fails_the_gate(tmp_path):
    site = _site_row(representative_pos="1000", n_carriers="1", samples="S1",
                     insertion_orientation="+")
    attached = _attach_one(site, [_call(1000, orientation="-")])
    assert common.verify_join([attached])["verdict"] == "join_disagrees_with_dedup_output"


def test_a_call_not_labelled_nested_is_reported_but_is_not_a_join_failure(tmp_path):
    """The legacy `NESTED` field is a different definition from site nesting.

    Gating on it fires on every antisense site on this cohort -- 100% of them
    carry the legacy `unnested` label while the producer, using the host-strand
    rule, calls them `nested_antisense`. That is a definitional difference, not
    a broken join, so it is counted and surfaced rather than made fatal.
    """
    site = _site_row(representative_pos="1000", n_carriers="1", samples="S1")
    call = _call(1000)
    call["nested_state"] = "unnested"
    attached = common.attach_call_details(
        [typed(site)], {"S1": {"chr1": [call]}}
    )[0]
    report = common.verify_join([attached])
    assert report["calls_where_source_nesting_label_differs_from_site"] == 1
    assert report["verdict"] == "join_consistent_with_dedup_output"


def test_a_call_of_the_wrong_family_fails_the_gate(tmp_path):
    """Joining a LINE-1 call onto an Alu site is a broken join, not a label difference.

    The family check is the one that actually proves the right record was
    recovered. `family` is carried on both sides, so this is checkable, and the
    `+/-10 bp` window is wide enough that neighbouring-family calls at the same
    locus are common enough to be worth pinning.
    """
    site = _site_row(representative_pos="1000", n_carriers="1", samples="S1",
                     family="ALU")
    call = _call(1000)
    call["family"] = "LINE1"
    attached = common.attach_call_details(
        [typed(site)], {"S1": {"chr1": [call]}}
    )[0]
    report = common.verify_join([attached])
    assert report["calls_where_family_disagrees"] == 1
    assert report["verdict"] == "join_disagrees_with_dedup_output"


# --------------------------------------------------------------------------
# host identity
# --------------------------------------------------------------------------


def test_the_host_key_is_the_element_interval_not_the_inserted_element():
    site = _site_row(host_start0="900", host_end0="1200", host_name="AluSx1")
    other = _site_row(site_id="US00002", family="SVA", representative_pos="1100")
    assert common.host_key(site) == common.host_key(other)
    assert common.host_key(site) != (
        "chr1", "AluSx1", 901, 1201, "+"
    )


def test_host_family_is_read_from_the_host_name():
    assert common.host_family(_site_row(host_name="AluSx1")) == "Alu"
    assert common.host_family(_site_row(host_name="L1PA2")) == "L1"
    assert common.host_family(_site_row(host_name="SVA_A")) == "SVA"
    assert common.host_family(_site_row(host_name="MIRb")) == "other"


def test_host_family_differs_from_the_inserted_element_family():
    site = _site_row(host_name="AluSx1", family="LINE1")
    assert common.host_family(site) == "Alu"
    assert site["family"] == "LINE1"


# --------------------------------------------------------------------------
# opportunity
# --------------------------------------------------------------------------


def test_the_host_annotation_does_not_mask_its_own_host():
    annotations = [(900, 1200, "AluSx1")]
    assert common.residual_mask_fraction((900, 1200, "AluSx1"), annotations) == 0.0


def test_a_neighbouring_repeat_does_reduce_opportunity():
    annotations = [(900, 1200, "AluSx1"), (1000, 1100, "AluY")]
    fraction = common.residual_mask_fraction((900, 1200, "AluSx1"), annotations)
    assert fraction == pytest.approx(100 / 300)


def test_annotations_outside_the_host_do_not_count():
    annotations = [(0, 100, "AluY"), (2000, 3000, "AluZ")]
    assert common.residual_mask_fraction((900, 1200, "AluSx1"), annotations) == 0.0


def test_a_degenerate_host_is_fully_masked():
    assert common.residual_mask_fraction((100, 100, "x"), []) == 1.0


def test_coverage_is_the_union_not_the_widest_single_annotation():
    """Two neighbours each covering 40% leave 20% of the host, not 60%.

    Taking the maximum overlap instead of the union reported a host with several
    intruding repeats as less masked than it is, which flattered the mappability
    diagnostic in exactly the direction that would have been preferred.
    """
    annotations = [
        (900, 1200, "AluSx1"),
        (900, 1020, "AluY"),
        (1080, 1200, "AluZ"),
    ]
    # 900-1020 and 1080-1200 is 240 of the host's 300.
    assert common.residual_mask_fraction((900, 1200, "AluSx1"), annotations) == pytest.approx(
        240 / 300
    )


def test_overlapping_neighbours_are_not_double_counted():
    annotations = [
        (900, 1200, "AluSx1"),
        (950, 1100, "AluY"),
        (1000, 1150, "AluZ"),
    ]
    # Union of 950-1100 and 1000-1150 is 950-1150, i.e. 200 of 300.
    assert common.residual_mask_fraction((900, 1200, "AluSx1"), annotations) == pytest.approx(
        200 / 300
    )


def test_annotations_are_clipped_to_the_host_before_merging():
    """An annotation hanging off the host contributes only its overlap."""
    annotations = [
        (900, 1200, "AluSx1"),
        (1150, 1400, "AluY"),  # only 1150-1200 lies inside the host
    ]
    assert common.residual_mask_fraction((900, 1200, "AluSx1"), annotations) == pytest.approx(
        50 / 300
    )


# --------------------------------------------------------------------------
# the RepeatMasker stream
# --------------------------------------------------------------------------


def _rmsk_line(contig: str, start: int, end: int, name: str = "AluSx") -> str:
    return "\t".join(
        [
            "17",  # sw ID
            str(start),  # alignment start
            str(end - start),  # alignment size, deliberately unused
            "40",  # score
            "ALU",  # strand
            contig,
            str(start),  # 0-based start
            str(end),  # 0-based end
            "(0)",  # per-rep
            name,
            "(97)",  # class code
        ]
    )


def test_the_rmsk_stream_drops_annotations_away_from_every_wanted_host(tmp_path):
    """The contig filter alone kept the whole genome's annotations in memory.

    `read_rmsk_intervals` is documented to keep only intervals overlapping a
    wanted host. It used to filter on contig only, so the caller received every
    annotation on all 23 in-scope contigs and had to re-filter per host.
    """
    path = tmp_path / "rmsk.txt"
    path.write_text(
        "\n".join(
            [
                _rmsk_line("chr1", 1000, 1100, "AluSx"),  # wanted
                _rmsk_line("chr1", 5000, 5100, "AluY"),  # same contig, elsewhere
                _rmsk_line("chr1", 20000, 20100, "AluZ"),  # same contig, elsewhere
                _rmsk_line("chr2", 1000, 1100, "AluS"),  # other contig, same coords
            ]
        )
        + "\n"
    )
    wanted = {"chr1": [(1000, 1200, "AluSx")]}
    keep = common.read_rmsk_intervals(path, wanted)
    assert keep["chr1"] == [(1000, 1100, "AluSx")]
    assert "chr2" not in keep


def test_the_rmsk_stream_keeps_a_host_touching_a_wanted_host_edge(tmp_path):
    """Half-open overlap: an annotation abutting the host start does not count."""
    path = tmp_path / "rmsk.txt"
    path.write_text(
        "\n".join(
            [
                _rmsk_line("chr1", 800, 1000, "A"),  # ends exactly at host start
                _rmsk_line("chr1", 1000, 1100, "B"),  # starts exactly at host start
            ]
        )
        + "\n"
    )
    keep = common.read_rmsk_intervals(path, {"chr1": [(1000, 1200, "H")]})
    assert keep["chr1"] == [(1000, 1100, "B")]
