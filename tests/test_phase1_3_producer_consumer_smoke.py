"""End-to-end smoke test: the Phase 0 cohort producer feeding Phases 1 and 2b.

Same premise as `test_phase4_producer_consumer_smoke.py`, one phase earlier.
`build_longread_nested_cohort.py` writes `per_call_longread_nested.csv`;
`test_position_enrichment.py` (Phase 1) and `score_tprt_hallmarks.py` (Phase 2b)
both read that one file and both apply the same four-part analysis-set gate:

    nested_in_alu_host == "1"  and  280 <= host_len <= 320
    and consensus_mapping == "unambiguous"  and  consensus_offset is an int

The three read the gate independently. `nested_analysis_set` is documented as
"reused verbatim so a Phase 2b number is never computed on a different
denominator than the Phase 1 number printed beside it" -- but it is a *re-implementation*
of `load_nested_rows`, not a call to it, so the two can and do drift. This
module pins them to each other against a table built from the producer's own
`OUTPUT_COLUMNS`, which is what makes it a producer/consumer test rather than
three fixtures agreeing with each other by construction.

Why the fixture is serialised through `OUTPUT_COLUMNS` rather than written
column-by-column. `build_longread_nested_cohort.main()` cannot be driven here:
it needs a site BCF, a genotype BCF, a 155 MB RepeatMasker, a FASTA, an Alu
consensus and a real aligner. What *defines* the contract is the 39-element
`OUTPUT_COLUMNS` list and the `_csv_cell` escaper, and both are importable. The
write loop itself is four lines reproduced below from `main()`; the risk that
lives in a column *name* is therefore covered by the real constant, and the risk
that lives in those four lines is covered by
`test_the_header_on_disk_is_the_producers_own_column_list`.

Runtime is a few seconds and it needs no workspace data.
"""

from __future__ import annotations

import csv
import importlib.util
import json
import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import score_tprt_hallmarks as tprt  # noqa: E402
import test_position_enrichment as tpe  # noqa: E402


def _load_producer():
    """Import the Phase 0 producer by path, as its own tests do.

    Loaded by path because `scripts/` is not a package and the module's argparse
    defaults point at the real workspace.
    """
    script = SCRIPTS / "build_longread_nested_cohort.py"
    spec = importlib.util.spec_from_file_location("build_longread_nested_cohort", script)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


producer = _load_producer()

#: Columns the two consumers read, with the line each requirement comes from.
#: Asserted against the producer's own `OUTPUT_COLUMNS` so a rename on either
#: side fails here rather than as an empty cohort three phases downstream.
CONSUMER_REQUIRED_COLUMNS = {
    # the four-part analysis-set gate, in both consumers
    "nested_in_alu_host",  # tpe.load_nested_rows, tprt.nested_analysis_set
    "host_len",  # both gates
    "consensus_mapping",  # both gates
    "consensus_offset",  # both gates, and the coordinate the 133 window is in
    "chrom",  # tpe.load_nested_rows
    "pos",  # tpe.load_nested_rows
    "host_start0",  # tpe.load_nested_rows
    "host_end0",  # tpe.load_nested_rows
    "host_strand",  # tpe.load_nested_rows
    "host_name",  # tpe.load_nested_rows
    "consensus_span_bp",  # tpe.load_nested_rows
    "site_allele_count",  # tpe.load_nested_rows
    "insert_strand",  # tpe.load_nested_rows
    "consensus_match_name",  # tpe.load_nested_rows
    "offset_drift_bp",  # tpe.load_nested_rows
    "perc_resolved",  # tpe.load_nested_rows
    "not_canonical",  # tpe.load_nested_rows
    # Phase 2b scores the TPRT triad off these
    "tsd_len",
    "tsd_seq",
    "tsd_len_is_sentinel",
    "polya_len",
    "polya_seq",
    "host_offset_5p_0based",  # tpe --coordinate alternative
    "family",
}

#: The four clauses that decide cohort membership. Distinct from
#: `tpe.REQUIRED_COHORT_COLUMNS`, which is that set plus the four identity
#: fields every downstream cell is keyed on. Both are required columns; only
#: these four can reject a row.
GATE_CLAUSES = (
    "nested_in_alu_host",
    "host_len",
    "consensus_mapping",
    "consensus_offset",
)

#: Every column `load_nested_rows` must find under the default coordinate,
#: tracked against the script's own declaration rather than a local copy, so
#: adding one there fails here. The alternate coordinate column is conditional
#: and is covered by its own test below.
REQUIRED_FOR_PHASE1 = (*tpe.REQUIRED_COHORT_COLUMNS, "consensus_offset")
ALTERNATE_COORDINATE = "host_offset_5p_0based"

assert set(GATE_CLAUSES) <= set(REQUIRED_FOR_PHASE1)
SOFT_COLUMNS = sorted(
    set(CONSUMER_REQUIRED_COLUMNS) - set(REQUIRED_FOR_PHASE1)
)

HOST_LEN = 300  # inside NEAR_FULL_HOST / NEAR_FULL_MIN..MAX
N_HOSTS = 8


def _row(
    *,
    chrom: str,
    pos: int,
    host_start0: int,
    host_len: int = HOST_LEN,
    host_strand: str = "+",
    insert_strand: str = "+",
    consensus_offset=130,
    nested: str = "1",
    consensus_mapping: str = "unambiguous",
    tsd_len: str = "16",
    tsd_seq: str = "ACGTACGTACGTACGT",
    allele_count: str = "1",
) -> dict:
    """One row in the producer's own vocabulary, defaults all in-gate."""
    return {
        "chrom": chrom,
        "pos": str(pos),
        "family": "ALU",
        "insert_strand": insert_strand,
        "conformation": "internal",
        "perc_resolved": "100",
        "not_canonical": "0",
        "rt_len": "300",
        "nested_in_alu_host": nested,
        "host_start0": str(host_start0),
        "host_end0": str(host_start0 + host_len),
        "host_strand": host_strand,
        "host_name": "AluSx1",
        "host_len": str(host_len),
        "host_offset_5p_0based": str(consensus_offset) if consensus_offset != "" else "",
        "host_offset_bin_20bp": "100-120",
        "host_selection_rule": "longest_span",
        "consensus_match_name": "AluSx1",
        "consensus_offset": str(consensus_offset) if consensus_offset != "" else "",
        "consensus_mapping": consensus_mapping,
        "offset_drift_bp": "0",
        "consensus_span_bp": "300",
        "consensus_offset_min": "0",
        "consensus_offset_max": "300",
        "alignment_identity": "99.0",
        "alignment_score": "900",
        "n_unaligned_host_bases": "0",
        "max_indel_bp": "1",
        "tsd_len": tsd_len,
        "tsd_seq": tsd_seq,
        # Sentinel values must read as "unevaluable", never as a real length.
        "tsd_len_is_sentinel": "0",
        "polya_len": "34",
        "polya_seq": "A" * 20,
        "site_allele_count": allele_count,
        "site_allele_freq": "0.5",
        "genotype_gt": "0/1",
        "genotype_state": "het",
        "genotype_state_reason": "",
        "genotype_callability_source": "fixture",
    }


def _fixture_rows() -> list[dict]:
    """In-gate rows on both peak and off-peak, plus one row per gate rejection.

    The excluded rows are the point. A fixture of only in-gate rows would pass
    with every gate inverted, because a gate that rejects nothing is
    indistinguishable from a gate that is not running.
    """
    rows: list[dict] = []
    for index in range(N_HOSTS):
        start0 = 1_000_000 + index * 10_000
        host_strand = "+" if index % 2 == 0 else "-"
        insert_strand = "+" if index % 3 else "-"
        # Three inside the pre-registered 133 +/- 5 peak window, three in the
        # mid-host control arm, so both comparison arms are non-empty.
        for offset in (130, 133, 135):
            rows.append(
                _row(
                    chrom="chr1", pos=start0 + offset + 1, host_start0=start0,
                    host_strand=host_strand, insert_strand=insert_strand,
                    consensus_offset=offset,
                )
            )
        for offset in (20, 60, 250):
            rows.append(
                _row(
                    chrom="chr1", pos=start0 + offset + 1, host_start0=start0,
                    host_strand=host_strand, insert_strand=insert_strand,
                    consensus_offset=offset,
                )
            )
    # One row per gate rejection, all on in-range chromosomes.
    rows.append(_row(chrom="chr1", pos=9_000_001, host_start0=9_000_000, nested="0"))
    rows.append(_row(chrom="chr1", pos=9_100_001, host_start0=9_100_000, host_len=200))
    rows.append(_row(chrom="chr1", pos=9_200_001, host_start0=9_200_000, host_len=400))
    rows.append(
        _row(chrom="chr1", pos=9_300_001, host_start0=9_300_000, consensus_mapping="ambiguous")
    )
    rows.append(_row(chrom="chr1", pos=9_400_001, host_start0=9_400_000, consensus_offset=""))
    return rows


def _write_cohort(rows: list[dict], target: Path) -> Path:
    """Serialize exactly as `build_longread_nested_cohort.main()` does.

    Reproduced from the write block at the end of `main()`: the header is the
    producer's `OUTPUT_COLUMNS`, and every cell goes through its own
    `_csv_cell`. Any column the producer stops emitting, or renames, changes
    this file's header -- which is the drift the tests below look for.
    """
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w") as fh:
        fh.write(",".join(producer.OUTPUT_COLUMNS) + "\n")
        for row in rows:
            fh.write(
                ",".join(producer._csv_cell(row.get(c, "")) for c in producer.OUTPUT_COLUMNS)
                + "\n"
            )
    return target


@pytest.fixture(scope="module")
def cohort(tmp_path_factory):
    root = tmp_path_factory.mktemp("smoke_phase13")
    rows = _fixture_rows()
    return {
        "root": root,
        "rows": rows,
        "cohort_csv": _write_cohort(rows, root / "per_call_longread_nested.csv"),
    }


def _read_back(cohort) -> list[dict[str, str]]:
    with cohort["cohort_csv"].open(newline="") as handle:
        return list(csv.DictReader(handle))


# --------------------------------------------------------------------------
# the schema contract, asserted against the producer's own constant
# --------------------------------------------------------------------------


def test_the_producer_emits_every_column_the_consumers_read():
    missing = sorted(CONSUMER_REQUIRED_COLUMNS - set(producer.OUTPUT_COLUMNS))
    assert not missing, (
        "build_longread_nested_cohort.OUTPUT_COLUMNS does not carry column(s) "
        f"{missing} that the Phase 1/2b consumers read. A consumer that finds "
        "its gate column missing gets an empty analysis set, which reads as a "
        "clean null result rather than a schema break."
    )


def test_the_header_on_disk_is_the_producers_own_column_list(cohort):
    with cohort["cohort_csv"].open() as handle:
        header = handle.readline().rstrip("\n").split(",")
    assert header == list(producer.OUTPUT_COLUMNS), (
        "the fixture header drifted from the producer's column list; the write "
        "loop in this test and the one in build_longread_nested_cohort.main() "
        "must stay identical"
    )


def test_the_fixture_survives_a_round_trip_through_the_producers_escaper(cohort):
    """Every value must land in one cell.

    `_csv_cell` is the only quoting in the pipeline. A value containing a comma
    that escaped it would shift every later column by one and silently shift
    `host_len` onto `consensus_mapping`, which the gates read.
    """
    rows = _read_back(cohort)
    assert len(rows) == len(cohort["rows"])
    for original, round_tripped in zip(cohort["rows"], rows):
        for column in producer.OUTPUT_COLUMNS:
            assert round_tripped[column] == producer._csv_cell(original.get(column, "")), (
                f"column {column} did not survive the producer's CSV escaper"
            )


# --------------------------------------------------------------------------
# the shared denominator
# --------------------------------------------------------------------------


def test_both_consumers_admit_the_same_rows_to_their_analysis_set(cohort):
    """The invariant `nested_analysis_set`'s docstring claims.

    Phase 1 prints a site count and Phase 2b prints a count beside it; the
    comment says the gate is reused verbatim so the two can never disagree on
    the denominator. It is a re-implementation, so that is an intention and not
    a guarantee. This asserts it against a table neither side wrote by hand.
    """
    rows = _read_back(cohort)
    phase1 = tpe.load_nested_rows(rows, coordinate="consensus_offset")
    phase2b = tprt.nested_analysis_set(rows)
    assert len(phase1) == len(phase2b), (
        f"Phase 1 admits {len(phase1)} rows and Phase 2b admits {len(phase2b)} "
        "from the same file; their side-by-side site counts are not comparable"
    )
    # And the gate is actually admitting the in-gate rows, not an empty set.
    assert len(phase1) == N_HOSTS * 6


def test_every_gate_clause_has_a_row_it_alone_rejects(cohort):
    """Each clause of the four-part gate must have something to reject.

    Without this, inverting any clause leaves the counts above unchanged and the
    gate is untested rather than correct. The test is per-row rather than
    per-file: it looks for rows that satisfy every clause *except* one, so a
    clause cannot look exercised merely because some other clause happened to
    drop the same row.
    """
    rows = _read_back(cohort)

    def fails(row, clause: str) -> bool:
        """Whether this row fails one clause of the gate."""
        if clause == "nested_in_alu_host":
            return row.get("nested_in_alu_host") != "1"
        if clause == "host_len":
            try:
                host_len = int(row["host_len"])
            except (KeyError, TypeError, ValueError):
                return True
            return not (tpe.NEAR_FULL_MIN <= host_len <= tpe.NEAR_FULL_MAX)
        if clause == "consensus_mapping":
            return row.get("consensus_mapping") != "unambiguous"
        if clause == "consensus_offset":
            try:
                int(row.get("consensus_offset", ""))
            except (TypeError, ValueError):
                return True
            return False
        raise AssertionError(f"unhandled gate clause {clause!r}")

    # The fixture rows must match this reader exactly, or the exclusion counts
    # below describe the reader rather than the real gate.
    assert len([r for r in rows if not any(fails(r, c) for c in GATE_CLAUSES)]) == (
        N_HOSTS * 6
    )

    for clause in GATE_CLAUSES:
        rejected_only_by = [
            row
            for row in rows
            if fails(row, clause)
            and not any(fails(row, other) for other in GATE_CLAUSES if other != clause)
        ]
        assert rejected_only_by, (
            f"no fixture row is rejected by the {clause!r} clause alone, so that "
            "clause is not exercised by this fixture"
        )


# --------------------------------------------------------------------------
# both consumers run to completion on the producer-shaped table
# --------------------------------------------------------------------------


def test_phase1_position_enrichment_runs_on_a_producer_shaped_table(cohort, tmp_path):
    out = tmp_path / "phase1"
    rc = tpe.main(
        [
            "--cohort", str(cohort["cohort_csv"]),
            "--outdir", str(out),
            # Absent is a supported state, not a crash: Null B records why.
            "--mapability-bed", str(cohort["root"] / "absent.bed"),
            "--replicates", "200",
        ]
    )
    assert rc == 0
    report = json.loads((out / "position_enrichment.json").read_text())
    primary = report["primary_test"]
    assert primary["n_sites"] == N_HOSTS * 6, (
        "Phase 1 scored a different denominator from "
        "`tpe.load_nested_rows` on the same file"
    )
    assert primary["observed"] > 0
    assert primary["expected_a"] is not None
    assert report["null_b"]["ran"] is False
    assert (out / "position_enrichment.md").exists()


def test_phase2b_tprt_hallmarks_runs_on_a_producer_shaped_table(cohort, tmp_path):
    out = tmp_path / "phase2b"
    rc = tprt.main(
        [
            "--cohort", str(cohort["cohort_csv"]),
            "--outdir", str(out),
            "--permutations", "200",
            # No reference: TSD-junction and EN-motif become `unevaluable`
            # rather than absent, which is the honest state for a fixture.
            "--no-reference",
        ]
    )
    assert rc == 0
    report = json.loads((out / "tprt_hallmarks.json").read_text())
    assert report["phase"] == "2b"
    assert report["set_sizes"]["analysis_set"] == N_HOSTS * 6
    assert report["set_sizes"]["peak"] > 0
    assert report["set_sizes"]["mid_host_control"] > 0
    assert report["set_sizes"]["unnested"] > 0, "the unnested comparison arm is empty"
    assert set(report["comparisons"]) == {"peak_vs_rest_of_host", "nested_vs_unnested"}
    for name, comparison in report["comparisons"].items():
        assert comparison["exposed_input"] > 0, f"{name} exposed arm is empty"
        assert comparison["reference_input"] > 0, f"{name} reference arm is empty"
    assert (out / "tprt_hallmarks.md").exists()


def test_both_phases_report_the_same_site_count(cohort, tmp_path):
    """The number a reader compares across the two reports must be one number."""
    phase1_out = tmp_path / "phase1"
    tpe.main(
        [
            "--cohort", str(cohort["cohort_csv"]),
            "--outdir", str(phase1_out),
            "--mapability-bed", str(cohort["root"] / "absent.bed"),
            "--replicates", "50",
        ]
    )
    phase2b_out = tmp_path / "phase2b"
    tprt.main(
        [
            "--cohort", str(cohort["cohort_csv"]),
            "--outdir", str(phase2b_out),
            "--permutations", "50",
            "--no-reference",
        ]
    )
    phase1 = json.loads((phase1_out / "position_enrichment.json").read_text())
    phase2b = json.loads((phase2b_out / "tprt_hallmarks.json").read_text())
    assert phase1["primary_test"]["n_sites"] == phase2b["set_sizes"]["analysis_set"], (
        "Phase 1 and Phase 2b disagree on the analysis-set size for the same "
        "input file, so their headline numbers are not comparable"
    )


# --------------------------------------------------------------------------
# drift: a schema break must be refused, not reported as a clean null
# --------------------------------------------------------------------------


def _drop_column(source: Path, target: Path, column: str) -> None:
    with source.open() as handle:
        header = handle.readline().rstrip("\n").split(",")
        body = handle.readlines()
    index = header.index(column)
    target.write_text(
        ",".join(name for i, name in enumerate(header) if i != index) + "\n"
        + "".join(
            ",".join(cell for i, cell in enumerate(line.rstrip("\n").split(",")) if i != index)
            + "\n"
            for line in body
        )
    )


@pytest.mark.parametrize("column", GATE_CLAUSES)
def test_removing_a_gate_column_stops_both_phases(cohort, tmp_path, column):
    """Neither phase may turn a missing column into a zero denominator.

    The two consumers now stop for different reasons, and that is the honest
    shape of the problem. Phase 1 validates the header and says so
    (`CohortSchemaError`). Phase 2b's `nested_analysis_set` is still a
    re-implementation that reads with `.get(...)`, so it empties silently and is
    caught one layer up by the empty-cohort refusal in `score_tprt_hallmarks.main`.
    Both are refusals; neither is a null result. If Phase 2b is ever given the
    same header check, its expected error changes and this test says so.
    """
    broken = tmp_path / f"drop_{column}.csv"
    _drop_column(cohort["cohort_csv"], broken, column)
    rows = _read_back({"cohort_csv": broken})
    assert tprt.nested_analysis_set(rows) == [], (
        f"dropping {column} should empty the analysis set; if it does not, the "
        "gate no longer reads this column and this test proves nothing"
    )
    with pytest.raises(tpe.CohortSchemaError):
        tpe.main(
            [
                "--cohort", str(broken),
                "--outdir", str(tmp_path / "phase1"),
                "--mapability-bed", str(cohort["root"] / "absent.bed"),
                "--replicates", "10",
            ]
        )
    with pytest.raises(SystemExit):
        tprt.main(
            [
                "--cohort", str(broken),
                "--outdir", str(tmp_path / "phase2b"),
                "--permutations", "10", "--no-reference",
            ]
        )


@pytest.mark.parametrize("column", REQUIRED_FOR_PHASE1)
def test_dropping_any_required_column_is_a_schema_error(cohort, tmp_path, column):
    """One access discipline: absence raises, it never empties the cohort.

    This replaces a pair of tests that asserted the three-way split the function
    used to have -- `row[...]` raised KeyError, `row.get(...)` on a gate column
    emptied the cohort silently, and `row.get(...)` on an evidence column
    degraded. The split was the defect: a missing column could quietly produce a
    zero denominator. Now every required column is validated against the header up
    front and reports the same way.
    """
    broken = tmp_path / f"drop_{column}.csv"
    _drop_column(cohort["cohort_csv"], broken, column)
    rows = _read_back({"cohort_csv": broken})
    with pytest.raises(tpe.CohortSchemaError):
        tpe.load_nested_rows(rows)
    # `main` propagates it rather than converting it into a null result.
    with pytest.raises(tpe.CohortSchemaError):
        tpe.main(
            [
                "--cohort", str(broken),
                "--outdir", str(tmp_path / "phase1"),
                "--mapability-bed", str(cohort["root"] / "absent.bed"),
                "--replicates", "10",
            ]
        )


def test_a_cohort_where_every_row_fails_the_gate_is_refused_not_nulled(
    cohort, tmp_path
):
    """The remaining empty-cohort path, and it is a data outcome.

    Columns present, values in range for the schema, but nothing qualifies. That
    is a legitimate thing for a table to be, so it is refused by `main` rather
    than raised as a schema error -- the two must stay distinguishable, or a real
    "nothing passes the gate" result would be reported as a broken pipeline.
    """
    all_unnested = _write_cohort(
        [{**row, "nested_in_alu_host": "0"} for row in cohort["rows"]],
        tmp_path / "all_unnested.csv",
    )
    rows = list(csv.DictReader(all_unnested.open(newline="")))
    assert tpe.load_nested_rows(rows) == []
    with pytest.raises(SystemExit):
        tpe.main(
            [
                "--cohort", str(all_unnested),
                "--outdir", str(tmp_path / "phase1"),
                "--mapability-bed", str(cohort["root"] / "absent.bed"),
                "--replicates", "10",
            ]
        )


def test_the_alternate_coordinate_column_is_required_only_when_selected(
    cohort, tmp_path
):
    """`host_offset_5p_0based` is conditional, and the check follows the choice.

    Requiring it unconditionally would reject a table that is perfectly readable
    under the default frame; not requiring it at all would let a run selected on
    that coordinate empty its cohort. So the required set is the base columns plus
    whichever coordinate was asked for.
    """
    broken = tmp_path / f"drop_{ALTERNATE_COORDINATE}.csv"
    _drop_column(cohort["cohort_csv"], broken, ALTERNATE_COORDINATE)
    rows = _read_back({"cohort_csv": broken})

    # Not selected: the run is unaffected.
    assert tpe.load_nested_rows(rows, coordinate="consensus_offset")

    # Selected: absence is a schema error, not an empty analysis set.
    with pytest.raises(tpe.CohortSchemaError):
        tpe.load_nested_rows(rows, coordinate=ALTERNATE_COORDINATE)
    with pytest.raises(tpe.CohortSchemaError):
        tpe.main(
            [
                "--cohort", str(broken),
                "--outdir", str(tmp_path / "phase1"),
                "--mapability-bed", str(cohort["root"] / "absent.bed"),
                "--coordinate", ALTERNATE_COORDINATE,
                "--replicates", "10",
            ]
        )


def test_dropping_a_soft_column_degrades_rather_than_empties(cohort, tmp_path):
    """The columns read with `.get()` must not be load-bearing for the cohort.

    `tsd_seq` and friends are evidence, not membership. If removing one emptied
    the analysis set it would be a gate column by accident, and the four-clause
    gate would be a lie.
    """
    for column in SOFT_COLUMNS:
        broken = tmp_path / f"drop_{column}.csv"
        _drop_column(cohort["cohort_csv"], broken, column)
        rows = _read_back({"cohort_csv": broken})
        assert tprt.nested_analysis_set(rows), (
            f"dropping {column} emptied the analysis set, so it is a required "
            "column and belongs in tpe.REQUIRED_COHORT_COLUMNS"
        )
        assert tpe.load_nested_rows(rows), f"dropping {column} emptied Phase 1"