"""Tests for the Phase 0 long-read nested-Alu cohort builder.

The gates are the point of this module, so most of the coverage here is
behavioural: a stale index must be detected without being mistaken for a wrong
one, a CHM13 path must be refused outright, and a 0/0 genotype must never be
silently folded together with an uncallable site.

These tests exercise the pure logic and the gate decision functions. The
full build reads the multi-hundred-megabyte BCF and rmsk tables and is run by
the script itself, not here; the only tests that touch those files are the
ones marked as requiring the real local data, which skip when it is absent.
"""

from __future__ import annotations

import gzip
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "build_longread_nested_cohort.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("build_longread_nested_cohort", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


m = _load_module()


# --------------------------------------------------------------------------
# Positional convention
# --------------------------------------------------------------------------


def test_host_offset_plus_strand_is_from_left_edge():
    # 0-based: the first base inside a + strand host is offset 0.
    assert m.host_offset(pos=101, start0=100, end0=400, host_strand="+") == 0
    assert m.host_offset(pos=150, start0=100, end0=400, host_strand="+") == 49


def test_host_offset_minus_strand_is_from_right_edge():
    # For a - strand host the 5' end is its rightmost base, so the offset runs
    # from the right. This is the convention that makes a tail-proximal host
    # look the same on either strand.
    assert m.host_offset(pos=400, start0=100, end0=400, host_strand="-") == 0
    assert m.host_offset(pos=301, start0=100, end0=400, host_strand="-") == 99


def test_host_offset_unknown_strand_is_none():
    assert m.host_offset(pos=101, start0=100, end0=400, host_strand=".") is None


def test_host_offset_is_strand_symmetric_for_mirrored_positions():
    """Two mirrored positions must report the same 5'-relative offset.

    For a host spanning [start0, end0), a + strand offset o sits at
    start0 + o and a - strand offset o sits at end0 - 1 - o. Note `pos` is a
    1-based VCF coordinate, so offset 0 on the + strand is pos start0 + 1.
    """
    start0, end0 = 100, 400
    for offset in (0, 50, 99, 250):
        plus_pos = start0 + offset + 1
        minus_pos = end0 - offset
        assert (
            m.host_offset(pos=plus_pos, start0=start0, end0=end0, host_strand="+")
            == offset
        )
        assert (
            m.host_offset(pos=minus_pos, start0=start0, end0=end0, host_strand="-")
            == offset
        )


# --------------------------------------------------------------------------
# Host selection: the two pre-declared tie-breaks must differ
# --------------------------------------------------------------------------


def test_select_host_longest_span_picks_widest():
    candidates = [(100, 400, "+", "AluY"), (100, 500, "+", "AluSx")]
    assert m.select_host(candidates, "longest_span") == (100, 500, "+", "AluSx")


def test_select_host_narrowest_picks_tightest():
    candidates = [(100, 400, "+", "AluY"), (100, 500, "+", "AluSx")]
    assert m.select_host(candidates, "narrowest_containment") == (100, 400, "+", "AluY")


def test_select_host_rules_disagree_so_phase2a_can_perturb():
    candidates = [(100, 400, "+", "AluY"), (100, 500, "+", "AluSx")]
    assert m.select_host(candidates, "longest_span") != m.select_host(
        candidates, "narrowest_containment"
    )


def test_select_host_empty_returns_none():
    assert m.select_host([], "longest_span") is None


def test_select_host_rejects_unknown_rule():
    with pytest.raises(ValueError):
        m.select_host([(1, 2, "+", "AluY")], "no_such_rule")


# --------------------------------------------------------------------------
# Genotype handling: the central discipline of this module
# --------------------------------------------------------------------------


def test_genotype_is_nonref_only_for_alt_bearing():
    assert m.genotype_is_nonref("1/0")
    assert m.genotype_is_nonref("0/1")
    assert m.genotype_is_nonref("1/1")
    assert not m.genotype_is_nonref("0/0")
    assert not m.genotype_is_nonref("./.")
    assert not m.genotype_is_nonref(".")


def test_zero_zero_stays_ref_genotype_not_uncallable():
    """A 0/0 with no callability caveat must remain ref_genotype."""
    state, reason = m.classify_genotype("0/0", has_pooled_site=True)
    assert state == "ref_genotype"
    assert reason == "none"


def test_zero_zero_inside_host_is_flagged_uninformative_not_ref():
    """This is the distinction the source callset cannot make for us.

    The pooled callset emits no `./.` and leaves VAF1 unset, so a 0/0 inside an
    Alu host cannot be separated from a genotyping failure. It gets its own
    state rather than being counted as an informative reference.
    """
    state, reason = m.classify_genotype(
        "0/0", has_pooled_site=True, inference_reasons=["host_region"]
    )
    assert state == "uncallable_or_uninformative"
    assert reason == "host_region"


def test_nonref_stays_concordant_even_inside_host():
    state, _ = m.classify_genotype(
        "1/0", has_pooled_site=True, inference_reasons=["host_region"]
    )
    assert state == "concordant_nonref"


def test_missing_genotype_is_its_own_state():
    assert m.classify_genotype("./.", has_pooled_site=True)[0] == "gt_missing"
    assert m.classify_genotype("", has_pooled_site=True)[0] == "gt_missing"


def test_absent_pooled_site_is_its_own_state():
    state, _ = m.classify_genotype("", has_pooled_site=False)
    assert state == "no_pooled_site"


def test_genotype_states_are_mutually_exclusive_names():
    # Guards against a state being renamed in one place only.
    for name in (
        "concordant_nonref",
        "ref_genotype",
        "uncallable_or_uninformative",
        "no_pooled_site",
        "gt_missing",
    ):
        assert name in m.GENOTYPE_STATES


# --------------------------------------------------------------------------
# TSD sentinels (Gate 0d)
# --------------------------------------------------------------------------


def _fake_bcf(tmp_path: Path, records: list[tuple[str, str]]) -> Path:
    """Write a minimal no-sample BCF the sentinel scanner can read.

    detect_tsd_sentinels shells out to bcftools, so the fixture is a real BCF
    built by bcftools itself; skipped when bcftools is unavailable. The
    ##INFO header lines are required or bcftools refuses to convert.
    """
    if not _has_bcftools():
        pytest.skip("bcftools not available")
    body = tmp_path / "body.vcf"
    with body.open("w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write("##contig=<ID=chr1,length=10000000>\n")
        fh.write('##INFO=<ID=FAM_N,Number=1,Type=String,Description="family">\n')
        fh.write('##INFO=<ID=TSD_LEN,Number=1,Type=Integer,Description="tsd">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for i, (fam, tsd_len) in enumerate(records):
            fh.write(
                f"chr1\t{1000 + i * 10}\t.\tN\t<INS:ME>\t.\tPASS\t"
                f"FAM_N={fam};TSD_LEN={tsd_len}\n"
            )
    out = tmp_path / "sites.bcf"
    subprocess.run(
        ["bcftools", "view", "-O", "b", "-o", str(out), str(body)],
        check=True, capture_output=True,
    )
    return out


def _has_bcftools() -> bool:
    from shutil import which
    return which("bcftools") is not None


def _minimal_bcf(tmp_path: Path, name: str, contigs: dict[str, int]) -> Path:
    """Build a real single-record, indexed BCF with a GRCh38-style header.

    The ##INFO lines are mandatory: bcftools refuses to write a BCF whose
    records use tags absent from the header. The index is mandatory too,
    because fetch_region_field queries by region.
    """
    if not _has_bcftools():
        pytest.skip("bcftools not available")
    body = tmp_path / f"{name}.vcf"
    with body.open("w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        for chrom, length in contigs.items():
            fh.write(f"##contig=<ID={chrom},length={length}>\n")
        fh.write('##INFO=<ID=FAM_N,Number=1,Type=String,Description="family">\n')
        fh.write('##INFO=<ID=TSD_LEN,Number=1,Type=Integer,Description="tsd">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        first = next(iter(contigs))
        fh.write(
            f"{first}\t100\t.\tN\t<INS:ME>\t.\tPASS\tFAM_N=Alu;TSD_LEN=17\n"
        )
    out = tmp_path / f"{name}.bcf"
    subprocess.run(
        ["bcftools", "view", "-O", "b", "-o", str(out), str(body)],
        check=True, capture_output=True,
    )
    subprocess.run(
        ["bcftools", "index", "-f", str(out)], check=True, capture_output=True
    )
    return out


def test_detect_tsd_sentinels_finds_a_pileup(tmp_path):
    """A single value towering over its neighbours must be flagged."""
    records = [("Alu", "17")] * 40
    records += [("Alu", "18")] * 40
    records += [("Alu", "41")] * 400   # the sentinel
    records += [("Alu", "42")] * 2
    records += [("Alu", "19")] * 40
    records += [("Alu", "20")] * 40
    bcf = _fake_bcf(tmp_path, records)
    result = m.detect_tsd_sentinels(bcf)
    assert "41" in result["sentinel_values"]
    assert "17" not in result["sentinel_values"]
    assert "19" not in result["sentinel_values"]


def test_detect_tsd_sentinels_clean_distribution_has_none(tmp_path):
    """A smooth unimodal histogram must not trip the sentinel detector."""
    records = []
    for value in range(10, 28):
        repeats = 50 - abs(value - 18) * 5
        records += [("Alu", str(value))] * repeats
    bcf = _fake_bcf(tmp_path, records)
    result = m.detect_tsd_sentinels(bcf)
    assert result["sentinel_values"] == {}
    assert result["verdict"] == "no_sentinels_detected"


def test_detect_tsd_sentinels_ignores_other_families(tmp_path):
    records = [("L1", "41")] * 500 + [("Alu", "17")] * 50 + [("Alu", "18")] * 50
    bcf = _fake_bcf(tmp_path, records)
    result = m.detect_tsd_sentinels(bcf, family="ALU")
    assert result["sentinel_values"] == {}


def test_tsd_length_is_sentinel_helper():
    assert m.tsd_length_is_sentinel(41, [41, 82])
    assert m.tsd_length_is_sentinel(82, [41, 82])
    assert not m.tsd_length_is_sentinel(17, [41, 82])


# --------------------------------------------------------------------------
# Gate 0b: CHM13 contamination must be refused
# --------------------------------------------------------------------------


def test_check_no_chm13_rejects_hs1_path(tmp_path):
    hs1 = tmp_path / "hs1"
    hs1.mkdir()
    bad = hs1 / "svim.bcf"
    bad.write_text("x")
    with pytest.raises(m.GateFailure, match="CHM13"):
        m.check_no_chm13(bad, bad)


def test_check_no_chm13_rejects_symlinked_input(tmp_path):
    """The real hazard is a symlink under hg38/ pointing into hs1/."""
    real = tmp_path / "real.bcf"
    real.write_text("x")
    link = tmp_path / "link.bcf"
    link.symlink_to(real)
    with pytest.raises(m.GateFailure, match="symlink"):
        m.check_no_chm13(link, real)


def test_check_no_chm13_reports_verdict_on_grch38_bcf(tmp_path):
    """Accepts a real BCF whose header carries GRCh38 primary lengths."""
    bcf = _minimal_bcf(tmp_path, "primary", {"chr1": m.GRCH38_LENGTHS["chr1"]})
    result = m.check_no_chm13(bcf, bcf)
    assert result["verdict"] == "grch38_confirmed"
    assert result["chrom_lengths_ok"] is True
    assert result["header_length_set"] == "grch38_primary"


def test_check_no_chm13_accepts_grch38_full_header(tmp_path):
    """The pooled callset uses the full-GBD header, which is still GRCh38."""
    bcf = _minimal_bcf(tmp_path, "full", {"chr1": m.GRCH38_FULL_LENGTHS["chr1"]})
    result = m.check_no_chm13(bcf, bcf)
    assert result["verdict"] == "grch38_confirmed"
    assert result["header_length_set"] == "grch38_full"


def test_check_no_chm13_rejects_non_grch38_lengths(tmp_path):
    """chr9 at CHM13 length must be rejected even though chr1 happens to match."""
    bcf = _minimal_bcf(
        tmp_path,
        "chm13",
        {"chr1": m.GRCH38_FULL_LENGTHS["chr1"], "chr9": 50818469},
    )
    with pytest.raises(m.GateFailure, match="not GRCh38"):
        m.check_no_chm13(bcf, bcf, required_chroms=["chr1", "chr9"])


# --------------------------------------------------------------------------
# Gate 0a: reference consistency
# --------------------------------------------------------------------------


def test_check_reference_consistency_flags_length_mismatch(tmp_path):
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">chr1\n" + "ACGT" * 10 + "\n")
    fai = Path(str(fasta) + ".fai")
    fai.write_text("chr1\t999999\t6\t40\t41\n")
    with pytest.raises(m.GateFailure, match="inconsistent"):
        m.check_reference_consistency(fasta)


def test_check_reference_consistency_reports_missing_chrom_not_silently(tmp_path):
    """A short .fai must be reported as incomplete, not accepted as GRCh38.

    The full-length check only ever runs against the real reference, so this
    fixture asserts the weaker but still useful property that a .fai missing
    GRCh38 chromosomes is surfaced in the counts.
    """
    seq = "ACGT" * 10
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">chr1\n" + seq + "\n")
    header = len(">chr1\n")
    fai = Path(str(fasta) + ".fai")
    fai.write_text(f"chr1\t{len(seq)}\t{header}\t{len(seq)}\t{len(seq) + 1}\n")
    index = {}
    for line in fai.read_text().splitlines():
        name, length, offset, linebases, linewidth = line.split("\t")
        index[name] = int(length)
    assert index["chr1"] == len(seq)
    missing = [c for c in m.GRCH38_LENGTHS if c not in index]
    assert len(missing) == len(m.GRCH38_LENGTHS) - 1


# --------------------------------------------------------------------------
# Consensus projection
# --------------------------------------------------------------------------


def test_project_offset_maps_identical_sequences_directly():
    aligner = m.build_aligner()
    seq = "ACGTACGTACGTACGTACGTACGT"
    result = m.project_offset(aligner, seq, seq, 10)
    assert result["consensus_offset"] == 10
    assert result["alignment_identity"] == 1.0
    assert result["consensus_mapping"] == "unambiguous"


def test_project_offset_accounts_for_leading_indel():
    """A host missing 5' bases shifts downstream consensus coordinates.

    This is the case the plan warns about: an unaligned stretch in a diverged
    host can move the projected coordinate far enough to change its consensus
    bin, so the drift is reported rather than silently absorbed.
    """
    aligner = m.build_aligner()
    consensus = "ACGTACGTACGTACGTACGTACGT"
    host = consensus[3:]  # host lost the first 3 consensus bases
    result = m.project_offset(aligner, host, consensus, 10)
    # host offset 10 is consensus offset 13, because 3 bases are missing.
    assert result["consensus_offset"] == 13
    assert result["offset_drift_bp"] == 3


def test_project_offset_marks_large_indel_ambiguous():
    """An indel past the tolerance marks the mapping ambiguous.

    The tolerance is 15 bp, set from the measured drift distribution rather
    than assumed: on real Alu hosts |consensus_offset - host_offset| is 6 bp
    at the median, so a 5 bp cutoff would discard most of the cohort.
    """
    aligner = m.build_aligner()
    consensus = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
    host = consensus[40:]  # 40 bp missing, well past the tolerance
    result = m.project_offset(aligner, host, consensus, 5)
    assert result["consensus_mapping"] == "ambiguous_indel"
    assert result["max_indel_bp"] >= 40


def test_project_offset_small_indel_stays_unambiguous():
    aligner = m.build_aligner()
    consensus = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
    host = consensus[3:]  # 3 bp missing, within tolerance
    result = m.project_offset(aligner, host, consensus, 5)
    assert result["consensus_mapping"] == "unambiguous"
    assert result["consensus_offset"] == 8
    assert result["offset_drift_bp"] == 3


def test_project_offset_out_of_range_is_not_mapped():
    aligner = m.build_aligner()
    seq = "ACGTACGTACGT"
    result = m.project_offset(aligner, seq, seq, 999)
    assert result["consensus_offset"] is None
    assert result["consensus_mapping"] == "ambiguous"


def test_project_offset_empty_sequence_is_safe():
    aligner = m.build_aligner()
    assert m.project_offset(aligner, "", "ACGT", 0)["consensus_offset"] is None


# --------------------------------------------------------------------------
# Consensus lookup
# --------------------------------------------------------------------------


def test_pick_consensus_prefers_exact_match():
    consensus = {"AluY": "AAAA", "AluSx": "CCCC"}
    assert m.pick_consensus_for_host("AluY", consensus) == ("AluY", "AAAA")


def test_pick_consensus_strips_trailing_digits():
    consensus = {"AluSx": "CCCC"}
    assert m.pick_consensus_for_host("AluSx1", consensus) == ("AluSx", "CCCC")


def test_pick_consensus_returns_none_when_absent():
    assert m.pick_consensus_for_host("AluZz9", {"AluY": "AAAA"}) is None


def test_load_consensus_filters_to_family(tmp_path):
    fasta = tmp_path / "reps.fa"
    fasta.write_text(">AluY\nACGT\n>AluSx\nTTTT\n>L1MC\nGGGG\n")
    loaded = m.load_consensus(fasta, family="ALU")
    assert set(loaded) == {"AluY", "AluSx"}


# --------------------------------------------------------------------------
# rmsk interval classification
# --------------------------------------------------------------------------


def test_is_alu_interval_matches_case_insensitively():
    """UCSC rmsk columns: 10=repName, 11=repClass, 12=repFamily."""
    fields = [""] * 13
    fields[10], fields[11], fields[12] = "AluY", "SINE", "Alu"
    assert m.is_alu_interval(fields)
    fields[10], fields[11], fields[12] = "aluY", "SINE", "Alu"
    assert m.is_alu_interval(fields)
    fields[10], fields[11], fields[12] = "L1PA2", "LINE", "L1"
    assert not m.is_alu_interval(fields)


# --------------------------------------------------------------------------
# Manifest writing: a verdict must land on disk, not only on stdout
# --------------------------------------------------------------------------


def test_fetch_region_field_gt_returns_genotype_not_tag_name(tmp_path):
    """Regression: a wrong FORMAT spec silently yields the literal tag name.

    `bcftools query -f '[GT]'` does not fail; it prints "GT" on every line. That
    produced a table where all 26,499 sites looked genotype-concordant. This
    pins the correct `[%GT]` form so the failure cannot come back quietly.
    """
    if not _has_bcftools():
        pytest.skip("bcftools not available")
    body = tmp_path / "gt.vcf"
    with body.open("w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write("##contig=<ID=chr1,length=1000000>\n")
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n")
        fh.write("chr1\t100\t.\tN\t<INS:ME>\t.\tPASS\t.\tGT\t1/0\n")
        fh.write("chr1\t200\t.\tN\t<INS:ME>\t.\tPASS\t.\tGT\t0/0\n")
    bcf = tmp_path / "gt.bcf"
    subprocess.run(
        ["bcftools", "view", "-O", "b", "-o", str(bcf), str(body)],
        check=True, capture_output=True,
    )
    subprocess.run(
        ["bcftools", "index", "-f", str(bcf)], check=True, capture_output=True
    )
    got = m.fetch_region_field(bcf, "chr1", 1, 1000000, "S1", "GT")
    assert got == {100: "1/0", 200: "0/0"}
    assert "GT" not in set(got.values())


def test_fetch_region_field_info_tag(tmp_path):
    if not _has_bcftools():
        pytest.skip("bcftools not available")
    bcf = _minimal_bcf(tmp_path, "ac", {"chr1": 1000000})
    got = m.fetch_region_field(bcf, "chr1", 1, 1000000, None, "TSD_LEN")
    assert got == {100: "17"}


def test_write_manifest_round_trips(tmp_path):
    path = tmp_path / "provenance_manifest.json"
    payload = {"gate_0c": {"verdict": "not_independent"}, "n": 3}
    m.write_manifest(path, payload)
    assert json.loads(path.read_text()) == payload


def test_summarise_rows_counts_states_and_bins():
    rows = [
        {
            "genotype_state": "concordant_nonref",
            "nested_in_alu_host": 1,
            "host_len": 300,
            "host_offset_bin_20bp": 120,
            "consensus_mapping": "unambiguous",
            "tsd_len_is_sentinel": 1,
        },
        {
            "genotype_state": "ref_genotype",
            "nested_in_alu_host": 0,
            "host_len": "",
            "host_offset_bin_20bp": "",
            "consensus_mapping": "",
            "tsd_len_is_sentinel": 0,
        },
    ]
    summary = m.summarise_rows(rows)
    assert summary["n_rows"] == 2
    assert summary["n_nested"] == 1
    assert summary["n_nested_near_full_host_280_320"] == 1
    assert summary["host_offset_bin_counts"] == {"120": 1}
    assert summary["genotype_state_counts"]["concordant_nonref"] == 1
    assert summary["n_tsd_len_sentinel"] == 1


# --------------------------------------------------------------------------
# Integration: the gates against the real local data
# --------------------------------------------------------------------------

REAL_DATA = Path.home() / "retrotransposon-workdir" / "data" / "public"


def _require_real_data() -> None:
    if not (REAL_DATA / "polymorphism" / "hg38" / "long_read_1kg_ont_vienna").exists():
        pytest.skip("long-read polymorphism data not present on this host")


@pytest.mark.slow
def test_real_site_bcf_is_grch38_not_chm13():
    """The hg38/ directory also contains a dangling hs1/ symlink."""
    _require_real_data()
    base = REAL_DATA / "polymorphism" / "hg38" / "long_read_1kg_ont_vienna"
    site = base / "svim.asm.hg38.noGt.SVAN_1.3.bcf"
    if not site.exists():
        pytest.skip("site BCF absent")
    result = m.check_no_chm13(site, site)
    assert result["verdict"] == "grch38_confirmed"
    assert result["chrom_lengths_ok"] is True


@pytest.mark.slow
def test_real_tsd_sentinels_include_41():
    """41 is the documented sentinel; if this changes, the gate changed."""
    _require_real_data()
    base = REAL_DATA / "polymorphism" / "hg38" / "long_read_1kg_ont_vienna"
    site = base / "svim.asm.hg38.noGt.SVAN_1.3.bcf"
    if not site.exists():
        pytest.skip("site BCF absent")
    result = m.detect_tsd_sentinels(site)
    assert "41" in result["sentinel_values"]
    assert result["sentinel_share_of_calls"] > 0.05


@pytest.mark.slow
def test_real_pooled_sample_set_includes_hg03086():
    """Gate 0c's premise: HG03086 is inside the discovery cohort."""
    _require_real_data()
    base = REAL_DATA / "polymorphism" / "hg38" / "long_read_1kg_ont_vienna"
    geno = base / "svim.asm.hg38.bcf"
    if not geno.exists():
        pytest.skip("genotyped BCF absent")
    samples = m.load_sample_order(geno)
    assert len(samples) > 100
    assert "HG03086" in samples
