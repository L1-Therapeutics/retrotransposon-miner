"""Frozen-rule conformance tests for the ten-genome production pipeline.

`scripts/analyze_ten_genome_mei.py` produced every number published in
`nested_analysis/results_multi_sample/`, and nothing tested it. The 13 tests in
`test_multi10_analysis.py` pin a re-implementation of the dedup core that no
production path imports, so a rule could be broken in the runner with the suite
still green. These tests drive the runner itself.

The rules below are the ones the analysis specification freezes. Each test names
the rule it defends so a future failure says which guarantee broke.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

REPO = Path(__file__).resolve().parents[1]
SCRIPTS = REPO / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import analyze_ten_genome_mei as ten  # noqa: E402
import mei_reference_opportunity as opportunity  # noqa: E402


@pytest.fixture(scope="module")
def engine():
    """The frozen dedup engine the runner loads, with the 10-sample manifest."""
    return ten.load_engine(str(SCRIPTS))


# ---------------------------------------------------------------------------
# fixtures
# ---------------------------------------------------------------------------
def vcf_line(chrom, pos, info, record_id="rec", extra=()):
    cols = [chrom, str(pos), record_id, "N", "<INS>", ".", ".", info, *extra]
    return "\t".join(cols)


def genome_wide_lines(family="ALU", nested="unnested", start=1000):
    """One call per primary chromosome, which the ingest gate requires."""
    return [
        vcf_line(chrom, start + index, f"MEIFAMILY={family};ORIENT=+;NESTED={nested}",
                 record_id=f"{chrom}-{index}")
        for index, chrom in enumerate(_primary_chroms())
    ]


def _primary_chroms():
    return [f"chr{i}" for i in range(1, 23)] + ["chrX"]


def baseline_lines(sample, hg03086_nested=261):
    """A callset that clears every ingest gate for `sample`.

    `read_calls` enforces genome-wide coverage, the binary NESTED alphabet, and
    an exact raw-NESTED count of 261 for HG03086, in that order.
    """
    lines = genome_wide_lines()
    if sample == "HG03086":
        lines += [vcf_line("chr1", 80_000 + i, "MEIFAMILY=ALU;ORIENT=+;NESTED=nested",
                           record_id=f"n{i}") for i in range(hg03086_nested)]
    return lines


def write_cohort(directory, mutate=None, **kwargs):
    for sample in ten.SAMPLES:
        lines = baseline_lines(sample, **kwargs)
        if mutate is not None:
            lines = mutate(sample, lines)
        write_callset(directory, sample, lines)


def write_callset(directory, sample, lines, genotype="0/1:50"):
    path = Path(directory) / f"{sample}.vcf"
    header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample
    path.write_text("\n".join(["##reference=hg38", header, *lines]) + "\n")
    return path


@pytest.fixture
def cohort(tmp_path, engine, monkeypatch):
    """A 10-sample cohort directory with the count baselines disabled."""
    monkeypatch.setattr(engine, "BASELINES", {})
    return tmp_path


def make_call(engine, sample, pos, *, family="ALU", orientation="+", chrom="chr1",
              raw_nested="unnested", source_id=None, source_index=0):
    return engine.Call(sample, source_index, chrom, pos, family, orientation,
                       raw_nested, {}, sample,
                       source_ids=[source_id or f"{sample}:{pos}"])


def host(engine, start0, end0, *, chrom="chr1", strand="+", name="AluSq", family="ALU"):
    return engine.Host(chrom, start0, end0, strand, name, family)


# ---------------------------------------------------------------------------
# cohort definition and the frozen opportunity caveat
# ---------------------------------------------------------------------------
def test_tier1_caveat_is_the_specified_wording():
    """Tier 1 must carry the caveat verbatim, not a paraphrase of it."""
    assert ten.CAVEAT == (
        "conditional on reference-genome opportunity; "
        "reference host spans are assumed present in the cohort"
    )
    assert ten.CAVEAT == opportunity.CAVEAT


def test_tier2_label_is_the_specified_wording():
    assert opportunity.TIER2 == "excluding known-MEI-overlapping hosts"


def test_cohort_is_ten_genomes_with_the_original_five_intact():
    assert len(ten.SAMPLES) == 10
    assert len(set(ten.SAMPLES)) == 10
    assert ten.ORIGINAL == ten.SAMPLES[:5]
    assert set(ten.ORIGINAL) <= set(ten.SAMPLES)


def test_gene_annotation_fields_are_never_parsed():
    """Gene/snpEff fields are never matching keys, so never read at all."""
    assert "GENE" not in ten.INFO_FIELDS
    assert "GENEID" not in ten.INFO_FIELDS
    assert "CSQ" not in ten.INFO_FIELDS
    assert "MEIFAMILY" in ten.INFO_FIELDS and "ORIENT" in ten.INFO_FIELDS


# ---------------------------------------------------------------------------
# ingest: GT/GQ and chrY
# ---------------------------------------------------------------------------
def test_genotype_columns_cannot_change_the_parse(cohort, engine):
    """GT/GQ are never parsed, so they cannot change a single parsed field."""
    results = []
    for index, genotype in enumerate(("0/1:50", "1/1:99", "./.:0", "0|1:3", "1/1:999999")):
        directory = cohort / f"gt{index}"
        directory.mkdir()

        def add_genotype(sample, lines, genotype=genotype):
            if sample != ten.SAMPLES[1]:
                return lines
            return [line + f"\tGT:GQ\t{genotype}" for line in lines]

        write_cohort(directory, add_genotype)
        calls, _ = ten.read_calls(engine, directory)
        results.append([(c.pos, c.family, c.orientation, c.raw_nested, sorted(c.info))
                        for c in calls[ten.SAMPLES[1]]])
    assert len(results[0]) == 23
    assert all(r == results[0] for r in results)


def test_gene_annotation_fields_never_reach_the_parsed_call(cohort, engine):
    """Gene/snpEff INFO keys are filtered out, not merely unused downstream."""
    directory = cohort / "genes"
    directory.mkdir()
    write_cohort(directory, lambda sample, lines: lines + [
        vcf_line("chr1", 95_000, "MEIFAMILY=ALU;ORIENT=+;NESTED=unnested;GENE=TTN;"
                                "GENEID=23327;CSQ=intron|TTN", record_id="g")])
    calls, _ = ten.read_calls(engine, directory)
    extra = [c for c in calls["HG01474"] if c.pos == 95_000]
    assert len(extra) == 1
    assert set(extra[0].info) <= ten.INFO_FIELDS
    assert "GENE" not in extra[0].info and "CSQ" not in extra[0].info


def test_chr_y_is_qc_counted_but_excluded_from_analysis(cohort, engine):
    """chrY is counted for QC and never reaches matching or analysis."""
    write_cohort(cohort, lambda sample, lines: lines + [
        vcf_line("chrY", 500, "MEIFAMILY=ALU;ORIENT=+;NESTED=unnested", record_id="y1"),
        vcf_line("chrY", 700, "MEIFAMILY=ALU;ORIENT=+;NESTED=unnested", record_id="y2")])
    calls, qc = ten.read_calls(engine, cohort)
    assert all(entry["chrY"] == 2 for entry in qc)
    assert all(c.chrom != "chrY" for sample in calls for c in calls[sample])
    expected = 9 * 23 + (23 + 261)  # nine plain callsets, plus HG03086's nested records
    assert sum(len(v) for v in calls.values()) == expected


def test_ingest_rejects_a_callset_missing_a_primary_chromosome(cohort, engine):
    """The genome-wide gate must actually fire."""
    for sample in ten.SAMPLES:
        kept = [line for line in baseline_lines(sample) if not line.startswith("chrX\t")]
        write_callset(cohort, sample, kept)
    with pytest.raises(engine.InputGateError, match="non-genome-wide"):
        ten.read_calls(engine, cohort)


def test_ingest_rejects_a_non_binary_nested_value(cohort, engine):
    """The VCF NESTED field is binary; anything else is a schema change."""
    write_cohort(cohort, lambda sample, lines: lines + [
        vcf_line("chr1", 90_000, "MEIFAMILY=ALU;ORIENT=+;NESTED=sense", record_id="s")])
    with pytest.raises(engine.InputGateError, match="invalid binary NESTED"):
        ten.read_calls(engine, cohort)


def test_hg03086_raw_nested_gate_is_exactly_261(cohort, engine):
    """The one hard-coded raw-NESTED count is enforced, not assumed."""
    ok = cohort / "ok"
    ok.mkdir()
    write_cohort(ok)
    _, qc = ten.read_calls(engine, ok)
    assert next(e for e in qc if e["sample"] == "HG03086")["raw_nested"] == 261

    short = cohort / "short"
    short.mkdir()
    write_cohort(short, hg03086_nested=260)
    with pytest.raises(engine.InputGateError, match="raw nested gate !=261"):
        ten.read_calls(engine, short)


def test_ingest_rejects_an_unsupported_family(cohort, engine):
    write_cohort(cohort, lambda sample, lines: lines + [
        vcf_line("chr1", 91_000, "MEIFAMILY=ERVL;ORIENT=+;NESTED=unnested")])
    with pytest.raises(engine.InputGateError, match="unsupported family"):
        ten.read_calls(engine, cohort)


# ---------------------------------------------------------------------------
# dedup: one carrier per sample, and no transitive chaining
# ---------------------------------------------------------------------------
def test_repeat_records_from_one_sample_count_as_one_carrier(engine):
    """Exactly one call per sample per site, so a carrier contributes once."""
    calls = [make_call(engine, "HG03086", 1000, source_id="A"),
             make_call(engine, "HG03086", 1005, source_id="B"),
             make_call(engine, "HG03086", 1009, source_id="C")]
    sites, collapsed = ten.dedup(engine, calls)
    assert collapsed == 2
    assert len(sites) == 1
    assert sites[0].samples == ("HG03086",)
    assert sites[0].anchor.source_count == 3


def test_groups_are_anchored_and_never_chain_transitively(engine):
    """A-B match and B-C match must not merge A, B and C into one site."""
    calls = [make_call(engine, "HG03086", 1000),
             make_call(engine, "HG01474", 1008),
             make_call(engine, "HG01566", 1016)]
    sites, _ = ten.dedup(engine, calls, window=10)
    groups = sorted(tuple(sorted(s.samples)) for s in sites)
    assert groups == [("HG01474", "HG03086"), ("HG01566",)]


def test_a_far_sample_may_not_drag_a_near_one_into_a_site(engine):
    """One out-of-window member disqualifies the whole candidate pair."""
    calls = [make_call(engine, "HG03086", 1000),
             make_call(engine, "HG01474", 1005),
             make_call(engine, "HG01566", 1100)]
    sites, _ = ten.dedup(engine, calls, window=10)
    assert sorted(len(s.samples) for s in sites) == [1, 2]


def test_repeated_runs_do_not_mutate_the_caller_s_provenance(engine):
    """`collapse_within_sample` must not write provenance back in place.

    The 5-genome sweep called `unique_sites` five times over the same Call
    objects. In-place aggregation made each pass re-sum its own inflated
    predecessor, so source_count walked 3 -> 4 -> 5 and source_call_ids grew
    "A|B|C" -> "A|B|B|B|B|C".
    """
    calls = [make_call(engine, "HG03086", 1000, source_id="A"),
             make_call(engine, "HG03086", 1003, source_id="B"),
             make_call(engine, "HG03086", 1006, source_id="C")]
    first = ten.dedup(engine, calls)
    second = ten.dedup(engine, calls)
    assert [s.anchor.source_ids for s in first[0]] == [s.anchor.source_ids for s in second[0]]
    assert [s.anchor.source_count for s in first[0]] == [s.anchor.source_count for s in second[0]]
    assert first[0][0].anchor.source_ids == ["A", "B", "C"]
    assert first[0][0].anchor.source_count == 3


# ---------------------------------------------------------------------------
# the four-state nesting enum
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("orientation, expected", [
    ("+", "nested_sense"),
    ("-", "nested_antisense"),
    ("", "nested_unknown"),
    (".", "nested_unknown"),
])
def test_nesting_state_is_recomputed_and_unknown_is_never_antisense(engine, orientation, expected):
    """Nesting is derived from host overlap plus orientation, not from NESTED."""
    hosts = {("chr1", "ALU"): [host(engine, 1000, 2000, strand="+")]}
    call = make_call(engine, "HG03086", 1500, orientation=orientation, raw_nested="nested")
    engine.assign_hosts([call], hosts)
    assert call.nested_state == expected


def test_unknown_orientation_is_a_state_of_its_own(engine):
    """Folding `unknown` into `antisense` would invent polarity."""
    hosts = {("chr1", "ALU"): [host(engine, 1000, 2000, strand="+")]}
    unknown = make_call(engine, "HG03086", 1500, orientation="")
    antisense = make_call(engine, "HG03086", 1500, orientation="-")
    engine.assign_hosts([unknown, antisense], hosts)
    assert unknown.nested_state != antisense.nested_state
    assert unknown.nested_state == "nested_unknown"


def test_uncontained_calls_are_unnested_not_nested_unknown(engine):
    """No containing host means `unnested`, not an unknown polarity."""
    hosts = {("chr1", "ALU"): [host(engine, 1000, 2000)]}
    call = make_call(engine, "HG03086", 9_000, raw_nested="nested")
    engine.assign_hosts([call], hosts)
    assert call.nested_state == "unnested"
    assert call.same_family_host is None


# ---------------------------------------------------------------------------
# host assignment and the matching key
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("pos, contained", [
    (1_000, False),  # start0 < POS fails at the left edge
    (1_001, True),
    (2_000, True),   # POS <= end0 holds at the right edge
    (2_001, False),
])
def test_host_containment_is_start_lt_pos_le_end(engine, pos, contained):
    """VCF POS is 1-based; an off-by-one here silently reassigns every host."""
    hosts = {("chr1", "ALU"): [host(engine, 1000, 2000)]}
    call = make_call(engine, "HG03086", pos)
    engine.assign_hosts([call], hosts)
    assert (call.same_family_host is not None) is contained


def test_the_longest_containing_host_wins(engine):
    short = host(engine, 1000, 1200, name="AluSq")
    long = host(engine, 1100, 1400, name="AluY")
    hosts = {("chr1", "ALU"): [short, long]}
    call = make_call(engine, "HG03086", 1150)
    engine.assign_hosts([call], hosts)
    assert call.same_family_host is long


def test_equal_length_hosts_tie_break_to_the_leftmost(engine):
    """Deterministic tie-break, or site identity depends on file order."""
    left = host(engine, 1000, 1200, name="AluY")
    right = host(engine, 1100, 1300, name="AluSq")
    hosts = {("chr1", "ALU"): [left, right]}
    call = make_call(engine, "HG03086", 1150)
    engine.assign_hosts([call], hosts)
    assert call.same_family_host is left


def test_the_matching_key_separates_family_orientation_and_host(engine):
    """The key is chrom + family + orientation + host element."""
    h = host(engine, 1000, 2000)
    base = make_call(engine, "HG03086", 1500)
    base.match_host = h
    twin = make_call(engine, "HG01474", 1502)
    twin.match_host = h
    assert base.match_key == twin.match_key

    for other in (
        make_call(engine, "HG03086", 1500, family="LINE1"),
        make_call(engine, "HG03086", 1500, orientation="-"),
        make_call(engine, "HG03086", 1500, chrom="chr2"),
    ):
        other.match_host = h
        assert other.match_key != base.match_key

    elsewhere = make_call(engine, "HG03086", 1500)
    elsewhere.match_host = host(engine, 1000, 2000, name="AluSx")
    assert elsewhere.match_key != base.match_key

    unnested = make_call(engine, "HG03086", 1500)
    assert unnested.match_host is None
    assert unnested.match_key != base.match_key


def test_anchored_groups_refuses_a_second_call_from_one_sample(engine):
    """Without the collapse step, one sample must still contribute one carrier."""
    h = host(engine, 1000, 2000)
    calls = [make_call(engine, "HG03086", 1500), make_call(engine, "HG03086", 1508),
             make_call(engine, "HG01474", 1504)]
    for call in calls:
        call.match_host = h
    sites = engine.anchored_groups(calls, 10, allow_same_sample=False)
    shared = [s for s in sites if len(s.samples) >= 2]
    assert len(shared) == 1
    assert sum(1 for member in shared[0].members if member.sample == "HG03086") == 1


def test_a_collapsed_record_is_kept_as_provenance_not_deleted(engine):
    """Dropping the extra record would understate the underlying evidence."""
    h = host(engine, 1000, 2000)
    calls = [make_call(engine, "HG03086", 1500, source_id="a"),
             make_call(engine, "HG03086", 1505, source_id="b"),
             make_call(engine, "HG01474", 1503, source_id="c")]
    for call in calls:
        call.match_host = h
    sites, collapsed = ten.dedup(engine, calls)
    assert collapsed == 1
    assert len(sites) == 1
    assert sites[0].samples == ("HG03086", "HG01474")
    assert sum(m.source_count for m in sites[0].members) == 3
    assert sorted(sites[0].anchor.source_ids) == ["a", "b"]


def test_known_mei_calls_reach_the_tier_two_exclusion_set(engine):
    """Tier 2 drops hosts holding a known-MEI breakpoint, whatever the nesting."""
    h = host(engine, 1000, 2000)
    hosts = {("chr1", "ALU"): [h]}
    known = make_call(engine, "HG03086", 1500, raw_nested="unnested")
    known.info["KNOWNMEI"] = "True"
    other = make_call(engine, "HG01474", 9_000)
    excluded = opportunity.known_exclusions(hosts, [known, other])
    assert excluded == {h.key}
    # exclusion is a separate overlay, not a mutation of the host universe
    assert hosts[("chr1", "ALU")] == [h]


# ---------------------------------------------------------------------------
# published table schema
# ---------------------------------------------------------------------------
def test_site_rows_carry_the_geometry_downstream_loaders_require(engine):
    """`load_unique_sites` rejects a table without these five columns.

    Emitting only the `same_family_host_*` names once made both Phase 4
    consumers die with a SchemaError on real data while unit tests passed.
    """
    import nested_multi_sample_common as common

    hosts = {("chr1", "ALU"): [host(engine, 1000, 2000)]}
    calls = [make_call(engine, "HG03086", 1500), make_call(engine, "HG01474", 1502)]
    engine.assign_hosts(calls, hosts)
    sites, _ = ten.dedup(engine, calls)
    rows = ten.rows_for_sites(engine, sites)
    required = {"host_name", "host_start0", "host_end0", "host_strand", "host_len"}
    assert required <= set(rows[0])
    assert required <= set(common.REQUIRED_UNIQUE_SITE_COLUMNS)
    assert rows[0]["host_name"] == "AluSq"
    assert rows[0]["host_len"] == 1000


def test_site_rows_carry_the_cohort_caveat_and_a_ten_bit_bitmap(engine):
    hosts = {("chr1", "ALU"): [host(engine, 1000, 2000)]}
    calls = [make_call(engine, "HG03086", 1500), make_call(engine, "HG01474", 1502)]
    engine.assign_hosts(calls, hosts)
    sites, _ = ten.dedup(engine, calls)
    row = ten.rows_for_sites(engine, sites)[0]
    assert row["opportunity_caveat"] == ten.CAVEAT
    assert len(row["presence_bitmap"]) == 10
    assert row["presence_bitmap"].count("1") == 2
    per_sample = {k: v for k, v in row.items() if k.startswith("present_")}
    assert sorted(per_sample) == sorted(f"present_{s}" for s in ten.SAMPLES)
    assert sum(per_sample.values()) == row["n_carriers"] == 2


# ---------------------------------------------------------------------------
# layers, Alu profile inputs, delta
# ---------------------------------------------------------------------------
def test_layers_partition_the_dedup_union_without_double_counting(engine):
    calls = [make_call(engine, "HG03086", 1500), make_call(engine, "HG01474", 1502),
             make_call(engine, "HG01566", 9000)]
    hosts = {("chr1", "ALU"): [host(engine, 1000, 2000)]}
    engine.assign_hosts(calls, hosts)
    sites, _ = ten.dedup(engine, calls)
    rows = {row["family"]: row for row in ten.layer_counts(sites)}
    alu = rows["ALU"]
    assert alu["union"] == 2
    assert alu["private"] == 1
    assert alu["shared"] == 1
    assert alu["private"] + alu["shared"] == alu["union"]
    assert set(alu["carrier_histogram"]) == {str(i) for i in range(1, 11)}


def test_alu_profile_accepts_only_hosts_in_the_declared_length_range(engine):
    """The linker/tail windows are defined on 280-320 bp Alu hosts."""
    lengths = (279, 280, 300, 320, 321)
    for length in lengths:
        h = host(engine, 1000, 1000 + length)
        call = make_call(engine, "HG03086", 1100)
        call.same_family_host = h
        call.match_host = h
        assert ten.alu_offsets([call]) == ([99] if 280 <= length <= 320 else [])


def test_alu_profile_honours_the_known_mei_exclusion_set(engine):
    h = host(engine, 1000, 1300)
    call = make_call(engine, "HG03086", 1100)
    call.same_family_host = h
    call.match_host = h
    assert ten.alu_offsets([call]) == [99]
    assert ten.alu_offsets([call], excluded={h.key}) == []


def test_unmasked_calls_drops_only_overlapping_breakpoints(engine):
    calls = [make_call(engine, "HG03086", 1000), make_call(engine, "HG01474", 2000)]
    masks = {"chr1": opportunity.SpanMask([(999, 1001)])}
    assert [c.pos for c in ten.unmasked_calls(calls, masks)] == [2000]
    assert ten.unmasked_calls(calls, None) == calls


# ---------------------------------------------------------------------------
# input/output provenance
#
# Two cohorts wrote `unique_sites.csv` at the same path, so the generation-1
# phase4 reports now name an input that no longer exists and cannot be
# reconstructed. Recording hashes of the inputs and of what a run overwrites
# turns that silent hazard into a recorded fact.
# ---------------------------------------------------------------------------
def fake_args(tmp_path, *, with_mappability=True):
    callsets = tmp_path / "callsets"
    callsets.mkdir()
    for sample in ten.SAMPLES:
        (callsets / f"{sample}.vcf").write_text(f"# {sample}\n")
    rmsk = tmp_path / "rmsk.txt.gz"
    rmsk.write_bytes(b"rmsk")
    fasta = tmp_path / "genome.fa"
    fasta.write_bytes(b">chr1\nACGT")
    consensus = tmp_path / "consensus.fa"
    consensus.write_bytes(b">L1HS\nACGT")
    low = tmp_path / "low.bed"
    if with_mappability:
        low.write_text("chr1\t0\t1\n")
    return SimpleNamespace(callset_dir=callsets, rmsk=rmsk, fasta=fasta,
                           consensus=consensus,
                           low_mappability=low if with_mappability else tmp_path / "absent.bed",
                           outdir=tmp_path / "out", replicates=1000)


def test_input_hashes_cover_every_input_the_run_reads(tmp_path, engine):
    hashes = ten.input_hashes(engine, fake_args(tmp_path))
    assert set(hashes) == {f"callsets/{s}.vcf" for s in ten.SAMPLES} | {
        "rmsk.txt.gz", "hg38.fasta", "consensus.fa", "low_mappability.bed"}
    assert all(len(h) == 64 for h in hashes.values())
    assert len(set(hashes.values())) == len(hashes)


def test_input_hashes_refuse_to_record_a_missing_input(tmp_path, engine):
    """A silently absent input would be recorded as if it had been used."""
    args = fake_args(tmp_path)
    args.fasta = tmp_path / "absent.fa"
    with pytest.raises(FileNotFoundError, match="hg38.fasta"):
        ten.input_hashes(engine, args)


def test_absent_optional_mask_is_omitted_rather_than_faked(tmp_path, engine):
    hashes = ten.input_hashes(engine, fake_args(tmp_path, with_mappability=False))
    assert "low_mappability.bed" not in hashes


def test_superseded_table_fingerprints_the_table_about_to_be_overwritten(tmp_path, engine):
    outdir = tmp_path / "out"
    outdir.mkdir()
    assert ten.superseded_table(engine, outdir) is None

    table = outdir / "unique_sites.csv"
    table.write_text("site_id,chrom\nUS00001,chr1\nUS00002,chr1\n")
    prior = ten.superseded_table(engine, outdir)
    assert prior["data_rows"] == 2
    assert prior["sha256"] == engine.sha256(table)

    # After the rewrite the fingerprint must describe the new bytes, not the old.
    table.write_text("site_id,chrom\nUS00001,chr1\n")
    assert ten.superseded_table(engine, outdir)["data_rows"] == 1
    assert ten.superseded_table(engine, outdir)["sha256"] != prior["sha256"]


def test_superseded_table_counts_a_header_only_table_as_zero_rows(tmp_path, engine):
    outdir = tmp_path / "out"
    outdir.mkdir()
    (outdir / "unique_sites.csv").write_text("site_id,chrom\n")
    assert ten.superseded_table(engine, outdir)["data_rows"] == 0


def test_output_hashes_cover_artifacts_and_never_the_manifest_itself(tmp_path, engine):
    """A file cannot contain its own hash; including it would be a lie."""
    outdir = tmp_path / "out"
    outdir.mkdir()
    (outdir / "unique_sites.csv").write_text("a\n")
    (outdir / "analysis_summary.md").write_text("b\n")
    (outdir / "outputs_sha256.json").write_text("{}")

    hashes = ten.output_hashes(engine, outdir)
    assert set(hashes) == {"unique_sites.csv", "analysis_summary.md"}
    assert hashes["unique_sites.csv"] == engine.sha256(outdir / "unique_sites.csv")


def test_provenance_records_inputs_and_supersession_not_only_scripts(tmp_path, engine):
    args = fake_args(tmp_path)
    outdir = tmp_path / "out"
    outdir.mkdir()
    (outdir / "unique_sites.csv").write_text("site_id\nUS00001\n")
    prior = ten.superseded_table(engine, outdir)

    prov = ten.provenance(Path(SCRIPTS), engine, args, prior)
    assert prov["inputs_sha256"]["hg38.fasta"] == engine.sha256(args.fasta)
    assert set(prov["scripts_sha256"]) == {
        "analyze_ten_genome_mei.py", "mei_reference_opportunity.py", "dedup_samples.py",
        "nested_multi_sample_common.py", "joint_enrichment.py", "recurrence_test.py"}
    assert prov["supersedes_unique_sites"]["data_rows"] == 1
    assert prov["opportunity_caveat"] == ten.CAVEAT
    assert prov["commit"] and len(prov["commit"]) == 40


def test_provenance_records_no_supersession_on_a_first_run(tmp_path, engine):
    prov = ten.provenance(Path(SCRIPTS), engine, fake_args(tmp_path), None)
    assert prov["supersedes_unique_sites"] is None


DELTA_KEYS = ("enrichment_ALU", "enrichment_LINE1", "enrichment_SVA", "linker", "tail",
              "replication_fraction", "union", "private", "shared", "nested_L1", "l1_projected")


def _headline(**overrides):
    base = {key: 1.0 for key in DELTA_KEYS}
    base.update(overrides)
    return {"headline": base}


def test_delta_flags_only_movers_beyond_twenty_percent():
    five = _headline(union=100, private=50, enrichment_ALU=2.0)
    ten_ = _headline(union=121, private=60, enrichment_ALU=2.5)
    rows = {row["headline"]: row for row in ten.delta_rows(five, ten_)}
    assert rows["Dedup union"]["flag_gt20pct"] is True
    assert round(rows["Dedup union"]["percent_change"], 6) == 21.0
    assert rows["Private"]["flag_gt20pct"] is False
    assert rows["Alu same-family enrichment"]["percent_change"] == 25.0
    assert rows["Alu same-family enrichment"]["flag_gt20pct"] is True


def test_delta_explains_a_flagged_mover_and_stays_silent_otherwise():
    rows = {row["headline"]: row
            for row in ten.delta_rows(_headline(union=100), _headline(union=200))}
    assert rows["Dedup union"]["reason"]
    quiet = {row["headline"]: row
             for row in ten.delta_rows(_headline(union=100), _headline(union=105))}
    assert quiet["Dedup union"]["reason"] == ""


def test_delta_reports_a_zero_baseline_as_no_change_not_a_division_error():
    rows = {row["headline"]: row
            for row in ten.delta_rows(_headline(union=0), _headline(union=10))}
    assert rows["Dedup union"]["percent_change"] is None
    assert rows["Dedup union"]["flag_gt20pct"] is False


# ---------------------------------------------------------------------------
# dependent_gate: the failure message must name the check that actually fired
# ---------------------------------------------------------------------------


def _join(**overrides):
    """A `verify_join` result, consistent unless a test says otherwise."""
    result = {
        "sites_total": 10,
        "sites_with_at_least_one_joined_call": 10,
        "sites_with_no_joined_call": 0,
        "calls_recovered": 20,
        "calls_declared_by_dedup": 20,
        "sites_where_carrier_count_disagrees": 0,
        "calls_where_orientation_disagrees": 0,
        "calls_where_family_disagrees": 0,
        "calls_where_source_nesting_label_differs_from_site": 592,
        "verdict": "join_consistent_with_dedup_output",
    }
    result.update(overrides)
    return result


def test_the_failure_reason_names_the_counter_that_fired():
    """The message used to assert a cause the gate no longer tests.

    It claimed the committed consumers require `NESTED=nested` on every
    recovered call. `verify_join` no longer gates on that field, so the message
    would have told a reader to look at the nesting definition and leave the
    consumers alone while the real cause sat in a different counter.
    """
    assert "family disagreement" in ten.join_failure_reason(
        _join(calls_where_family_disagrees=3)
    )
    assert "carrier-count disagreement" in ten.join_failure_reason(
        _join(sites_where_carrier_count_disagrees=1)
    )
    assert "orientation disagreement" in ten.join_failure_reason(
        _join(calls_where_orientation_disagrees=7)
    )


def test_the_failure_reason_reports_every_counter_that_fired():
    reason = ten.join_failure_reason(
        _join(
            sites_where_carrier_count_disagrees=2,
            calls_where_family_disagrees=1,
        )
    )
    assert "carrier-count disagreement" in reason
    assert "family disagreement" in reason


def test_the_failure_reason_never_blames_the_nesting_label():
    """The legacy NESTED comparison is a diagnostic, not a gate.

    592 calls differ on that field on the real cohort and the verdict is still
    `join_consistent_with_dedup_output`. A message that named the label as the
    cause would be wrong for every run that actually fails.
    """
    reason = ten.join_failure_reason(_join(calls_where_family_disagrees=1))
    assert "NESTED" in reason  # it is mentioned, to rule it out
    assert "not part of this verdict" in reason
    assert "not zero recurrence" not in reason


def test_the_failure_reason_refuses_to_guess_when_no_counter_fires():
    """An unexplained failure must not acquire a confident invented cause."""
    reason = ten.join_failure_reason(
        _join(verdict="join_disagrees_with_dedup_output")
    )
    assert "not authoritative" in reason


def test_the_failure_reason_cannot_claim_the_old_nesting_cause():
    """Directly pins the regression: the stale text must not come back."""
    reason = ten.join_failure_reason(_join(calls_where_orientation_disagrees=4))
    for stale in ("require every recovered source call",
                  "nesting-definition/input-contract mismatch",
                  "Consumers were not patched"):
        assert stale not in reason


FAKE_COMMON = '''
"""Stand-in for `nested_multi_sample_common.py`, written to a real path.

`dependent_gate` loads the committed consumer by file path, so the honest way to
test it is to give it a real file to load rather than to patch `importlib`.
"""
SITES = [{"site_id": "US00001", "carriers": ["HG03086"], "chrom": "chr1",
          "pos": 100, "n_carriers": 1}]

def load_unique_sites(path):
    return SITES, {"rows_read": 1}

def load_callsets(callsets, samples):
    return {}

def attach_call_details(sites, callsets):
    return sites

def verify_join(sites):
    return VERDICT
'''


def _write_fake_scripts(tmp_path, verdict):
    """A `scripts/` directory holding a `nested_multi_sample_common.py` stub."""
    scripts = tmp_path / "scripts"
    scripts.mkdir(exist_ok=True)
    (scripts / "nested_multi_sample_common.py").write_text(
        f"VERDICT = {verdict!r}\n" + FAKE_COMMON
    )
    return scripts


def test_a_leftover_scratch_directory_does_not_break_a_rerun(tmp_path, monkeypatch):
    """`mkdir(exist_ok=False)` made the gate single-use after any abnormal exit.

    A kill between creating the scratch directory and moving its contents out
    left it behind, and every later run then failed with a FileExistsError
    naming a temporary directory rather than the cause.
    """
    scripts = _write_fake_scripts(tmp_path, _join())
    outdir = tmp_path / "out"
    outdir.mkdir()
    (outdir / "unique_sites.csv").write_text("site_id\n")
    stale = outdir / "phase4_run"
    stale.mkdir()
    (stale / "half-written.json").write_text("{}")

    monkeypatch.setattr(
        ten.subprocess, "run",
        lambda *a, **k: SimpleNamespace(returncode=0, stdout="", stderr=""),
    )
    result = ten.dependent_gate(scripts, outdir, tmp_path / "callsets")

    assert result["status"] == "ready"
    assert not stale.exists(), "the leftover scratch directory was not cleaned up"
    assert not (outdir / "phase4_run").exists(), "the scratch directory outlived the run"


def test_the_scratch_directory_is_removed_even_when_a_consumer_fails(tmp_path, monkeypatch):
    """A failing consumer must not leave the next run unable to start."""
    scripts = _write_fake_scripts(tmp_path, _join())
    outdir = tmp_path / "out"
    outdir.mkdir()
    (outdir / "unique_sites.csv").write_text("site_id\n")

    monkeypatch.setattr(
        ten.subprocess, "run",
        lambda *a, **k: SimpleNamespace(returncode=1, stdout="", stderr="boom"),
    )
    result = ten.dependent_gate(scripts, outdir, tmp_path / "callsets")

    assert result["status"] == "failed"
    assert not (outdir / "phase4_run").exists()


def test_an_incompatible_join_records_the_fired_counter_and_no_reason_to_blame(
    tmp_path,
):
    """End to end through the real gate, with a real failing verdict."""
    failing = _join(calls_where_family_disagrees=4,
                    verdict="join_disagrees_with_dedup_output")
    scripts = _write_fake_scripts(tmp_path, failing)
    outdir = tmp_path / "out"
    outdir.mkdir()
    (outdir / "unique_sites.csv").write_text("site_id\n")

    result = ten.dependent_gate(scripts, outdir, tmp_path / "callsets")
    assert result["status"] == "incompatible"
    assert "family disagreement" in result["reason"]
    assert "require every recovered source call" not in result["reason"]
    assert not (outdir / "phase4_run").exists(), "the consumers must not run"
    recorded = json.loads((outdir / "phase4_status.json").read_text())
    assert recorded["status"] == "incompatible"
    assert recorded["reason"] == result["reason"], (
        "the provenance artifact must record the same reason the run reported"
    )
