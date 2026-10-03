"""End-to-end smoke test: the dedup producer feeding both Phase 4 consumers.

The unit tests on either side of this boundary were both green while the pipeline
was broken. `dedup_samples.py` wrote host geometry under the `same_family_host_*`
names; `joint_enrichment.py` and `recurrence_test.py` load the same file through
`load_unique_sites`, which requires `host_name`/`host_start0`/`host_end0`/
`host_strand`/`host_len`. On the real cohort both consumers died instantly with
`SchemaError`, and nothing failed for six days, because:

  * the producer's own tests asserted on `site_row()` in isolation, so the table
    it produced was never given to a loader; and
  * each consumer's tests hand-wrote their own `unique_sites.csv` with the
    schema the consumer *wanted*, so the fixture agreed with the reader by
    construction and could not disagree with the writer.

That second point is the important one. A hand-written fixture cannot detect
drift between a producer and a consumer, because it is the consumer's own
belief about the schema rather than the producer's output. Every assertion below
therefore starts from a table produced by the real
`dedup_samples.build_sites_table()` and written to disk exactly as `run()` writes
it, so the CSV on disk is genuinely the producer's, not a transcription of it.

Scope. `dedup_samples.main()` cannot be driven here, and that is deliberate: its
`parse_callset` enforces hard-coded per-sample `BASELINES` counts, a chrY count
of exactly 1 for HG03172/NA18498 and 0 elsewhere, and a genome-wide check, so it
refuses any cohort but the real one. This module exercises the part of the
producer that defines the contract -- `unique_sites()` -> `build_sites_table()` ->
`to_csv` -- and then hands the result to both real consumers through their real
`main()`. That is the whole boundary; the QC gates above it are not the contract
under test and are covered by their own tests.

Runtime is a few seconds and it needs no workspace data, which is the point: the
previous failure mode was a schema break that only surfaced after a full
five-genome cohort run over a 155 MB RepeatMasker table.
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

import joint_enrichment as je  # noqa: E402
import nested_multi_sample_common as common  # noqa: E402
import recurrence_test as rt  # noqa: E402


def _load_producer():
    """Import `dedup_samples.py` by path, as `test_dedup_samples.py` does.

    Loaded by path rather than imported because `scripts/` is not a package and
    the module is a script with module-level argparse defaults pointing at the
    real workspace.
    """
    script = SCRIPTS / "dedup_samples.py"
    spec = importlib.util.spec_from_file_location("dedup_samples", script)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


producer = _load_producer()

WINDOW = producer.PRIMARY_WINDOW

#: Two samples, matching the shape the consumers' join gate expects: a shared
#: site must name both carriers and both must resolve from their own callset.
#:
#: These are two of the five real cohort names because `dedup_samples.SAMPLES` is
#: a module constant that `collapse_within_sample` and `anchored_groups` both
#: index into (`SAMPLES.index(call.sample)`). A sample name outside that tuple
#: does not raise -- it is silently skipped, so every call is dropped and the
#: producer returns zero sites, which reads as "no nested insertions found"
#: rather than "this fixture used the wrong sample name".
SAMPLE_A = "HG03086"
SAMPLE_B = "HG01474"

#: 300 bp hosts so the two pre-specified position bins are both reachable:
#: `P1_primary` is 120-140 and `S1_tail_sense` is 280-300, the latter clipping
#: to offsets 280-299 of a 300 bp host.
HOST_LEN = 300
N_HOSTS = 20

#: (offset within host, insertion orientation, producer-derived nesting state).
#: The nesting state is set from the same rule the producer uses -- compare the
#: insertion strand to the host strand -- rather than hardcoded, so the fixture
#: cannot drift away from `Call.nested_state`.
EVENTS = (
    (130, "+"),  # inside P1_primary's bin, sense
    (285, "+"),  # inside S1_tail_sense's bin, sense
    (10, "-"),   # antisense, so the interaction cells are not vacuous
)


VCF_PREAMBLE = "\n".join(
    [
        "##fileformat=VCFv4.2",
        '##INFO=<ID=NESTED,Number=1,Type=String,Description="nested state">',
        '##INFO=<ID=ORIENT,Number=1,Type=String,Description="orientation">',
        '##INFO=<ID=MEIFAMILY,Number=1,Type=String,Description="family">',
        '##INFO=<ID=MEISUBFAMILY,Number=1,Type=String,Description="subfamily">',
        '##INFO=<ID=TSD,Number=1,Type=String,Description="tsd">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    ]
)

#: The legacy binary label these callsets actually carry. Using it here is
#: deliberate: it is what makes the consumers' NESTED census report
#: `nested_orientation_unlabeled` rather than fabricating an orientation, and it
#: keeps the join gate on the family/geometry checks where it belongs.
RAW_NESTED = "nested"
SUBFAMILY = "AluYb9#SINE/Alu"
TSD = "ACGTACGTACGTACGT"


def _nested_state(orientation: str, host) -> str:
    """The producer's own rule, so the fixture cannot disagree with it."""
    if orientation == host.strand:
        return "nested_sense"
    if orientation in {"+", "-"}:
        return "nested_antisense"
    return "nested_unknown"


def _build_calls() -> list:
    """Synthetic calls carrying real `Host` objects and real nesting states.

    Even-numbered hosts are shared by both samples so the consumers' join gate
    has multi-carrier sites to recover; odd-numbered hosts are private to one
    sample so the carrier-count check has both denominators to disagree about if
    it is going to.
    """
    calls = []
    for index in range(N_HOSTS):
        start0 = 1_000_000 + index * 10_000
        host = producer.Host("chr1", start0, start0 + HOST_LEN, "+", f"AluSx{index}", "ALU")
        carriers = (SAMPLE_A, SAMPLE_B) if index % 2 == 0 else (SAMPLE_A,)
        for offset, orientation in EVENTS:
            # `Call.pos` is 1-based VCF POS and `pos0` is the zero-based start.
            pos = host.start0 + offset + 1
            for sample in carriers:
                calls.append(
                    producer.Call(
                        sample=sample,
                        source_index=len(calls),
                        chrom="chr1",
                        pos=pos,
                        family="ALU",
                        orientation=orientation,
                        raw_nested=RAW_NESTED,
                        info={},
                        sample_column=sample,
                        same_family_host=host,
                        match_host=host,
                        nested_state=_nested_state(orientation, host),
                        source_count=1,
                        source_ids=[f"{sample}-{host.name}-{offset}"],
                    )
                )
    return calls


def _write_callsets(calls: list, callset_dir: Path) -> Path:
    """One VCF per sample, carrying exactly the calls the producer saw.

    `attach_call_details` joins by position within the window rather than by
    record ID, so these IDs need not match `site_id`; they are written
    independently anyway so that the fixture cannot accidentally assert the ID
    relationship the consumers do not rely on.
    """
    callset_dir.mkdir(parents=True, exist_ok=True)
    for sample in (SAMPLE_A, SAMPLE_B):
        lines = [VCF_PREAMBLE, f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}"]
        for call in calls:
            if call.sample != sample:
                continue
            call_id = call.source_ids[0]
            info = (
                f"NESTED={RAW_NESTED};ORIENT={call.orientation};MEIFAMILY={call.family};"
                f"MEISUBFAMILY={SUBFAMILY};TSD={TSD}"
            )
            lines.append(
                f"{call.chrom}\t{call.pos}\t{call_id}\tN\t<INS>\t.\t.\t{info}\tGT\t0/1"
            )
        (callset_dir / f"{sample}.vcf").write_text("\n".join(lines) + "\n")
    return callset_dir


@pytest.fixture(scope="module")
def cohort(tmp_path_factory):
    """A producer-written `unique_sites.csv` plus the callsets behind it.

    Module-scoped: the table and the two consumer runs are immutable inputs, so
    rebuilding them per test would triple the cost of a suite whose whole purpose
    is to be fast. The mutation tests below copy the CSV before editing it.
    """
    root = tmp_path_factory.mktemp("smoke_cohort")
    calls = _build_calls()
    sites, _collapsed = producer.unique_sites(calls, WINDOW)
    assert sites, "fixture produced no sites"

    # Exactly the producer's own write, from `dedup_samples.run`.
    frame = producer.build_sites_table(sites, WINDOW)
    unique_sites = root / "unique_sites.csv"
    frame.to_csv(unique_sites, index=False, quoting=csv.QUOTE_MINIMAL)
    callset_dir = _write_callsets(calls, root / "callsets")
    return {
        "root": root,
        "unique_sites": unique_sites,
        "callset_dir": callset_dir,
        "calls": calls,
        "frame": frame,
    }


# --------------------------------------------------------------------------
# the boundary: producer output is loadable by the consumers' reader
# --------------------------------------------------------------------------


def test_the_producer_table_carries_every_column_the_loader_requires(cohort):
    """The literal schema contract, checked on the producer's real output.

    `REQUIRED_UNIQUE_SITE_COLUMNS` is the loader's own declaration of what it
    reads. Asserting the producer satisfies it, on a table the producer built,
    is the direct statement of the contract that was broken. It is written
    against the header on disk rather than the frame in memory so that it tests
    what a consumer would actually open.
    """
    with cohort["unique_sites"].open(newline="") as handle:
        header = next(csv.reader(handle))
    missing = [c for c in common.REQUIRED_UNIQUE_SITE_COLUMNS if c not in header]
    assert not missing, (
        "dedup_samples.build_sites_table did not emit column(s) "
        f"{missing} that nested_multi_sample_common.load_unique_sites requires.\n"
        f"  producer header: {header}"
    )


def test_the_producer_table_loads_through_the_consumers_reader(cohort):
    sites, report = common.load_unique_sites(cohort["unique_sites"])
    assert sites, "the producer's own table loaded to an empty cohort"
    assert report["rows_read"] == len(cohort["frame"])
    assert report["rows_in_scope"] == len(sites)
    assert report["dropped_not_nested"] == 0, (
        "every fixture site is nested, so the cohort gate must not drop any"
    )


def test_the_load_report_reports_what_it_excluded(cohort):
    """A loader that drops rows must say which rule dropped them.

    The whole incident was a table that loaded as empty, so the report is what
    distinguishes "nothing was in the cohort" from "everything was excluded".
    """
    _sites, report = common.load_unique_sites(cohort["unique_sites"])
    assert report["cohort_rule"] == common.COHORT_RULE
    assert report["cohort_nesting_states"], "the cohort census is empty"
    assert report["dropped_out_of_scope_contig"] == 0


# --------------------------------------------------------------------------
# the join gate, which is what makes the per-call detail trustworthy
# --------------------------------------------------------------------------


def test_the_consumers_join_reproduces_the_producer_s_carrier_counts(cohort):
    sites, _report = common.load_unique_sites(cohort["unique_sites"])
    carriers = sorted({s for site in sites for s in site["carriers"]})
    callsets = common.load_callsets(cohort["callset_dir"], carriers)
    derived = common.verify_join(common.attach_call_details(sites, callsets))

    assert derived["verdict"] == "join_consistent_with_dedup_output", derived
    assert derived["sites_with_no_joined_call"] == 0
    assert derived["calls_recovered"] == derived["calls_declared_by_dedup"]
    assert derived["calls_recovered"] == sum(site["n_carriers"] for site in sites)
    assert derived["sites_where_carrier_count_disagrees"] == 0
    assert derived["calls_where_family_disagrees"] == 0
    assert derived["calls_where_orientation_disagrees"] == 0


# --------------------------------------------------------------------------
# both consumers run to completion on the producer's output
# --------------------------------------------------------------------------


def test_joint_enrichment_runs_on_a_producer_written_table(cohort, tmp_path):
    out = tmp_path / "joint"
    rc = je.main(
        [
            "--unique-sites", str(cohort["unique_sites"]),
            "--callset-dir", str(cohort["callset_dir"]),
            # No RepeatMasker is needed to prove the schema contract; the
            # opportunity diagnostic records its own absence.
            "--rmsk", str(cohort["root"] / "absent.rmsk.gz"),
            "--outdir", str(out),
            "--replicates", "200",
            "--bootstrap", "50",
        ]
    )
    assert rc == 0
    report = json.loads((out / "joint_enrichment.json").read_text())
    assert report["analysis"] == "joint_nested_signature_enrichment"
    assert [cell["cell_id"] for cell in report["cells"]] == [c.cell_id for c in je.CELLS]
    assert report["load_report"]["rows_in_scope"] > 0
    assert report["join_verification"]["verdict"] == "join_consistent_with_dedup_output"
    assert (out / "joint_enrichment.md").exists()


def test_recurrence_runs_on_a_producer_written_table(cohort, tmp_path):
    out = tmp_path / "recurrence"
    rc = rt.main(
        [
            "--unique-sites", str(cohort["unique_sites"]),
            "--callset-dir", str(cohort["callset_dir"]),
            "--outdir", str(out),
            "--chance-draws", "200",
        ]
    )
    assert rc == 0
    report = json.loads((out / "recurrence.json").read_text())
    assert report["analysis"] == "same_host_recurrence_ibd_vs_independent"
    assert report["load_report"]["rows_in_scope"] > 0
    assert report["join_verification"]["verdict"] == "join_consistent_with_dedup_output"
    assert (out / "recurrence.md").exists()
    assert (out / "recurrence_pairs.csv").exists()


# --------------------------------------------------------------------------
# the drift cases. Without these the module only proves the happy path works,
# which is not the claim; the claim is that a schema break is *caught*.
# --------------------------------------------------------------------------


def _rewrite_header(source: Path, target: Path, transform, drop=()) -> None:
    """Copy the producer's CSV with one column renamed or removed.

    Renaming has to touch the row keys as well as the header, or `DictWriter`
    rejects the rows instead of producing the broken table.
    """
    with source.open(newline="") as handle:
        reader = csv.DictReader(handle)
        header = [transform(name) for name in reader.fieldnames if name not in drop]
        rows = [
            {transform(k): v for k, v in row.items() if k not in drop}
            for row in reader
        ]
    with target.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header)
        writer.writeheader()
        writer.writerows(rows)


def test_the_host_length_regression_is_refused_not_silently_nulled(cohort, tmp_path):
    """The original bug, reproduced and pinned.

    Emitting `same_family_host_length` instead of `host_len` is exactly what
    `dedup_samples.site_row` did. The consumers must raise rather than proceed,
    because proceeding on an empty cohort produces a well-formed report full of
    zeros that reads as a decisive negative.
    """
    broken = tmp_path / "no_canonical_host_columns.csv"
    _rewrite_header(
        cohort["unique_sites"],
        broken,
        lambda name: "same_family_host_length" if name == "host_len" else name,
    )
    with pytest.raises(common.SchemaError):
        common.load_unique_sites(broken)

    with pytest.raises(common.SchemaError):
        je.main(
            [
                "--unique-sites", str(broken),
                "--callset-dir", str(cohort["callset_dir"]),
                "--rmsk", str(cohort["root"] / "absent.rmsk.gz"),
                "--outdir", str(tmp_path / "joint"),
                "--replicates", "10", "--bootstrap", "5",
            ]
        )

    with pytest.raises(common.SchemaError):
        rt.main(
            [
                "--unique-sites", str(broken),
                "--callset-dir", str(cohort["callset_dir"]),
                "--outdir", str(tmp_path / "recurrence"),
                "--chance-draws", "10",
            ]
        )


def test_a_table_with_no_cohort_column_is_refused(cohort, tmp_path):
    """Cohort membership must not be able to fail open.

    `same_family_nested_state` decides which sites are in the cohort. If the
    column is absent and the loader tolerates it, the cohort predicate has
    nothing to test, everything reads as unnested, and both consumers report an
    empty cohort -- which is the same failure mode as the host-geometry break,
    one column over.
    """
    broken = tmp_path / "no_cohort_column.csv"
    _rewrite_header(
        cohort["unique_sites"],
        broken,
        lambda name: name,
        drop={"same_family_nested_state"},
    )

    with pytest.raises(common.SchemaError):
        common.load_unique_sites(broken)
    with pytest.raises(common.SchemaError):
        je.main(
            [
                "--unique-sites", str(broken),
                "--callset-dir", str(cohort["callset_dir"]),
                "--rmsk", str(cohort["root"] / "absent.rmsk.gz"),
                "--outdir", str(tmp_path / "joint"),
                "--replicates", "10", "--bootstrap", "5",
            ]
        )


# --------------------------------------------------------------------------
# guards on the fixture itself
# --------------------------------------------------------------------------


def test_the_fixture_is_not_vacuously_nested(cohort):
    """If this drifts, the run-until-completion tests above stop meaning much.

    A smoke test over a cohort that is entirely one orientation class, or
    entirely single-carrier, would still pass while testing almost nothing. The
    four properties below are the ones the pre-registered cells and the join
    gate actually consume.
    """
    sites, _report = common.load_unique_sites(cohort["unique_sites"])
    states = {site["site_nested_state"] for site in sites}
    assert states == {"nested_sense", "nested_antisense"}, states
    assert {site["orientation"] for site in sites} == {"+", "-"}

    carrier_counts = {site["n_carriers"] for site in sites}
    assert carrier_counts == {1, 2}, carrier_counts
    assert any(site["private"] for site in sites)
    assert any(not site["private"] for site in sites)


def test_the_fixture_exercises_both_pre_specified_position_bins(cohort):
    """`P1_primary` is 120-140 and `S1_tail_sense` is 280-300.

    If the fixture only populated one of them, the other cell would be
    structurally zero and the run would pass without ever computing the
    quantity the analysis exists to measure.
    """
    sites, _report = common.load_unique_sites(cohort["unique_sites"])
    offsets = [site["host_offset"] for site in sites]
    assert any(cell.bin_lo <= off <= cell.bin_hi for off in offsets for cell in je.CELLS)
    for cell in (je.CELLS[0], je.CELLS[1]):
        hits = sum(cell.bin_lo <= off <= cell.bin_hi for off in offsets)
        assert hits > 0, f"cell {cell.cell_id} has no fixture events in its bin"