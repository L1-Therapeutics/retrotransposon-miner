"""Unit tests for retro_miner.gene_annotation (offline; VEP transport is faked).

Set RTM_LIVE_VEP=1 to additionally run the live rest.ensembl.org check.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from retro_miner import gene_annotation as ga

FIXTURE = Path(__file__).parent / "data" / "chr22_mei.vcf"


# ----------------------------------------------------------------- helpers


def _rec(chrom="chr22", pos=19223382, alt="<INS:ME:LINE1>", info="SVTYPE=INS;SVLEN=6018;END=19223382;MEINFO=L1HS,1,6018,+"):
    return ga.VcfRecord(chrom, pos, ".", "N", alt, ".", ".", ga.parse_info(info), ["GT", "./."])


def _vep_json(input_str, txs=None, most_severe=None):
    rec = {"input": input_str, "transcript_consequences": txs or []}
    if most_severe:
        rec["most_severe_consequence"] = most_severe
    return rec


class FakeVep:
    """Records POST bodies; answers every variant as intronic in a fixed gene unless overridden."""

    def __init__(self, responder=None, statuses=None):
        self.bodies: list[list[str]] = []
        self.responder = responder
        self.statuses = list(statuses or [])
        self.sleeps: list[float] = []

    def post(self, url, body, timeout):
        variants = json.loads(body)["variants"]
        self.bodies.append(variants)
        status = self.statuses.pop(0) if self.statuses else 200
        if status != 200:
            return status, b"slow down", {"retry-after": "2"}
        if self.responder:
            payload = [self.responder(v) for v in variants]
        else:
            payload = [
                _vep_json(v, [{"gene_symbol": "GENEX", "gene_id": "ENSG0", "consequence_terms": ["intron_variant"]}], "intron_variant")
                for v in variants
            ]
        return 200, json.dumps(payload).encode(), {}

    def sleep(self, s):
        self.sleeps.append(s)


# ----------------------------------------------------------------- VCF I/O


def test_fixture_reads_30_records_and_headers():
    headers, records = ga.read_vcf(FIXTURE)
    assert headers[-1].startswith("#CHROM")
    assert len(records) == 30
    assert all(r.alt.startswith("<INS:ME:") for r in records)
    assert {r.info["SVTYPE"] for r in records} == {"INS"}


def test_parse_info_handles_flags_and_values():
    d = ga.parse_info("SVTYPE=INS;IMPRECISE;SVLEN=10")
    assert d == {"SVTYPE": "INS", "IMPRECISE": None, "SVLEN": "10"}
    assert ga.parse_info(".") == {}


def test_record_roundtrip_preserves_columns(tmp_path):
    headers, records = ga.read_vcf(FIXTURE)
    out = tmp_path / "rt.vcf"
    ga.write_vcf(out, headers, records)
    assert [ln.rstrip("\n") for ln in open(FIXTURE) if ln.strip()] == [ln.rstrip("\n") for ln in open(out) if ln.strip()]


def test_add_info_headers_inserts_before_chrom_and_is_idempotent():
    headers = ["##fileformat=VCFv4.4", "#CHROM\tPOS"]
    once = ga.add_info_headers(headers)
    assert once[-1].startswith("#CHROM")
    assert sum(h.startswith("##INFO=<ID=GENE,") for h in once) == 1
    twice = ga.add_info_headers(once)
    assert twice == once


# ----------------------------------------------------------------- SVLEN regression


def test_vep_region_string_drops_svlen_and_pins_end_to_pos():
    """Regression for the VEP span-widening bug: SVLEN must not reach VEP by default."""
    s = ga.vep_region_string(_rec())
    assert s == "chr22 19223382 . N <INS:ME:LINE1> . . SVTYPE=INS;END=19223382"
    assert "SVLEN" not in s


def test_vep_region_string_keep_svlen_opt_in():
    s = ga.vep_region_string(_rec(), keep_svlen=True)
    assert "SVLEN=6018" in s


def test_vep_region_string_overrides_wrong_end_for_symbolic_ins():
    r = _rec(info="SVTYPE=INS;END=19229400")  # END != POS is invalid for a symbolic INS
    assert "END=19223382" in ga.vep_region_string(r)


def test_output_vcf_keeps_svlen_even_though_request_drops_it(tmp_path):
    fake = FakeVep()
    out = tmp_path / "o.vcf"
    ga.annotate_vcf(FIXTURE, out, post=fake.post, sleep=fake.sleep)
    _, recs = ga.read_vcf(out)
    assert all("SVLEN" in r.info for r in recs)
    assert all("SVLEN" not in v for batch in fake.bodies for v in batch)


# ----------------------------------------------------------------- VEP client


def test_query_vep_batches_at_200():
    fake = FakeVep()
    variants = [f"chr22 {1000 + i} . N <INS:ME:ALU> . . SVTYPE=INS" for i in range(201)]
    out = ga.query_vep(variants, post=fake.post, sleep=fake.sleep)
    assert [len(b) for b in fake.bodies] == [200, 1]
    assert len(out) == 201


def test_query_vep_rejects_bad_batch_size():
    with pytest.raises(ValueError):
        ga.query_vep(["x"], batch_size=201, post=FakeVep().post)


def test_query_vep_honours_429_retry_after():
    fake = FakeVep(statuses=[429, 429, 200])
    out = ga.query_vep(["chr22 5 . N <INS:ME:ALU> . . SVTYPE=INS"], post=fake.post, sleep=fake.sleep)
    assert fake.sleeps == [2.0, 2.0]
    assert len(out) == 1


def test_query_vep_raises_on_persistent_error():
    fake = FakeVep(statuses=[500])
    with pytest.raises(RuntimeError, match="HTTP 500"):
        ga.query_vep(["chr22 5 . N <INS:ME:ALU> . . SVTYPE=INS"], post=fake.post, sleep=fake.sleep)


# ----------------------------------------------------------------- consequence summary


def test_most_severe_term_follows_ensembl_ranking():
    assert ga.most_severe_term(["intron_variant", "coding_sequence_variant", "downstream_gene_variant"]) == "coding_sequence_variant"
    assert ga.most_severe_term(["intron_variant", "splice_acceptor_variant"]) == "splice_acceptor_variant"
    assert ga.most_severe_term([]) == "intergenic_variant"
    assert ga.most_severe_term(["made_up_term", "intron_variant"]) == "intron_variant"


def test_severity_table_has_no_duplicates_and_expected_endpoints():
    assert len(set(ga.CONSEQUENCE_SEVERITY)) == len(ga.CONSEQUENCE_SEVERITY)
    assert ga.CONSEQUENCE_SEVERITY[0] == "transcript_ablation"
    assert ga.CONSEQUENCE_SEVERITY[-1] == "sequence_variant"


def test_summarize_dedups_genes_and_orders_terms_by_severity():
    rec = _vep_json(
        "chr22 1 . N <INS:ME:ALU> . . SVTYPE=INS",
        [
            {"gene_symbol": "CLTCL1", "gene_id": "ENSG00000070371", "consequence_terms": ["intron_variant"]},
            {"gene_symbol": "CLTCL1", "gene_id": "ENSG00000070371", "consequence_terms": ["intron_variant", "NMD_transcript_variant"]},
            {"gene_symbol": "OTHER", "gene_id": "ENSG1", "consequence_terms": ["upstream_gene_variant"]},
        ],
        "intron_variant",
    )
    ann = ga.summarize_vep_record(rec)
    assert ann.gene_symbols == ("CLTCL1", "OTHER")
    assert ann.gene_ids == ("ENSG00000070371", "ENSG1")
    assert ann.most_severe == "intron_variant"
    assert ann.terms == ("intron_variant", "NMD_transcript_variant", "upstream_gene_variant")
    assert ann.n_transcripts == 3
    assert ann.info_fields() == {
        "GENE": "CLTCL1,OTHER",
        "GENEID": "ENSG00000070371,ENSG1",
        "CSQ": "intron_variant",
        "CSQ_TERMS": "intron_variant,NMD_transcript_variant,upstream_gene_variant",
        "CSQ_NTX": "3",
    }


def test_summarize_intergenic_has_no_gene_fields():
    rec = {"input": "chr22 1 . N <INS:ME:ALU> . . SVTYPE=INS", "intergenic_consequences": [{"consequence_terms": ["intergenic_variant"]}], "most_severe_consequence": "intergenic_variant"}
    ann = ga.summarize_vep_record(rec)
    assert ann.gene_symbols == () and ann.n_transcripts == 0
    assert "GENE" not in ann.info_fields()
    assert ann.info_fields()["CSQ"] == "intergenic_variant"


# ----------------------------------------------------------------- end-to-end (offline)


def test_annotate_vcf_end_to_end_offline(tmp_path):
    def responder(v):
        pos = int(v.split()[1])
        if pos == 19223382:
            return _vep_json(v, [{"gene_symbol": "CLTCL1", "gene_id": "ENSG00000070371", "consequence_terms": ["intron_variant"]}], "intron_variant")
        if pos == 49029650:
            return {"input": v, "intergenic_consequences": [{"consequence_terms": ["intergenic_variant"]}], "most_severe_consequence": "intergenic_variant"}
        return _vep_json(v, [{"gene_symbol": "G", "gene_id": "E", "consequence_terms": ["intron_variant"]}], "intron_variant")

    fake = FakeVep(responder)
    out, tsv = tmp_path / "ann.vcf", tmp_path / "ann.tsv"
    stats = ga.annotate_vcf(FIXTURE, out, tsv_path=tsv, post=fake.post, sleep=fake.sleep)

    assert stats.n_records == 30 and stats.n_annotated == 30 and stats.n_unmatched == 0
    assert stats.n_with_gene == 29

    headers, recs = ga.read_vcf(out)
    assert len(recs) == 30
    assert any(h.startswith("##INFO=<ID=CSQ,") for h in headers)
    by_pos = {r.pos: r for r in recs}
    assert by_pos[19223382].info["GENE"] == "CLTCL1"
    assert by_pos[19223382].info["CSQ"] == "intron_variant"
    assert by_pos[49029650].info["CSQ"] == "intergenic_variant"
    assert "GENE" not in by_pos[49029650].info
    # original columns untouched
    _, orig = ga.read_vcf(FIXTURE)
    assert [(r.chrom, r.pos, r.id, r.alt, r.rest) for r in orig] == [(r.chrom, r.pos, r.id, r.alt, r.rest) for r in recs]

    rows = [ln.rstrip("\n").split("\t") for ln in open(tsv)]
    assert rows[0] == list(ga.TSV_COLUMNS)
    assert len(rows) == 31
    line1 = next(r for r in rows[1:] if r[1] == "19223382")
    assert line1[rows[0].index("GENE")] == "CLTCL1"


def test_annotate_records_counts_unmatched_instead_of_dropping():
    fake = FakeVep(responder=lambda v: _vep_json("chrZ 1 . N <INS:ME:ALU> . . X"))  # never matches
    _, recs = ga.read_vcf(FIXTURE)
    out, stats = ga.annotate_records(recs, post=fake.post, sleep=fake.sleep)
    assert len(out) == 30
    assert stats.n_unmatched == 30 and stats.n_annotated == 0
    assert all("CSQ" not in r.info for r in out)


# ----------------------------------------------------------------- live (opt-in)


@pytest.mark.skipif(os.environ.get("RTM_LIVE_VEP") != "1", reason="set RTM_LIVE_VEP=1 to hit rest.ensembl.org")
def test_live_vep_chr22_svlen_regression(tmp_path):
    """The measured 2026-09-18 behaviour: with SVLEN dropped, chr22:19223382 is intronic in CLTCL1."""
    _, recs = ga.read_vcf(FIXTURE)
    out, stats = ga.annotate_records(recs)
    assert stats.n_unmatched == 0
    by_pos = {r.pos: r for r in out}
    assert "CLTCL1" in by_pos[19223382].info["GENE"]
    assert by_pos[19223382].info["CSQ"] == "intron_variant"
    assert "coding_sequence_variant" not in by_pos[19223382].info.get("CSQ_TERMS", "")
    assert by_pos[49029650].info["CSQ"] == "intergenic_variant"


# ----------------------------------------------------------------- snpEff backend

SNPEFF_FIXTURE = Path(__file__).parent / "data" / "chr22_mei.snpeff.vcf"


def test_snpeff_fixture_is_real_output_for_same_30_loci():
    _, vep_in = ga.read_vcf(FIXTURE)
    _, snp_out = ga.read_vcf(SNPEFF_FIXTURE)
    assert len(snp_out) == 30
    assert [(r.chrom, r.pos, r.alt) for r in snp_out] == [(r.chrom, r.pos, r.alt) for r in vep_in]
    assert all("ANN" in r.info for r in snp_out)
    assert snp_out[0].chrom == "chr22"  # snpEff echoed the chr prefix back


def test_parse_snpeff_ann_ranks_by_severity_not_snpeff_order():
    """snpEff lists upstream_gene_variant first for ADA2 even though it also reports intron_variant."""
    _, snp_out = ga.read_vcf(SNPEFF_FIXTURE)
    ada2 = next(r for r in snp_out if r.pos == 17224410)
    assert ada2.info["ANN"].split(",")[0].split("|")[1] == "upstream_gene_variant"
    ann = ga.parse_snpeff_ann(ada2.info["ANN"])
    assert ann.most_severe == "intron_variant"
    assert ann.gene_symbols == ("ADA2",)
    assert ann.terms[0] == "intron_variant"


def test_parse_snpeff_ann_intergenic_does_not_leak_flanking_genes():
    ann = ga.parse_snpeff_ann("N|intergenic_region|MODIFIER|CECR3-CECR9|ENSG1-ENSG2|intergenic_region|ENSG1-ENSG2|||||||||")
    assert ann.most_severe == "intergenic_variant"
    assert ann.gene_symbols == () and ann.gene_ids == () and ann.n_transcripts == 0
    assert ann.terms == ("intergenic_variant",)


def test_parse_snpeff_ann_splits_ampersand_terms_and_counts_transcripts():
    ann = ga.parse_snpeff_ann(
        "N|intron_variant&non_coding_transcript_variant|MODIFIER|G|ENSG|transcript|ENST1|lncRNA|||||||||,"
        "N|upstream_gene_variant|MODIFIER|G|ENSG|transcript|ENST2|protein_coding||||||||1234|"
    )
    assert ann.n_transcripts == 2
    assert ann.terms == ("intron_variant", "non_coding_transcript_variant", "upstream_gene_variant")
    assert ann.gene_ids == ("ENSG",)


def test_snpeff_svlen_locus_is_intronic_without_workaround():
    _, snp_out = ga.read_vcf(SNPEFF_FIXTURE)
    l1 = next(r for r in snp_out if r.pos == 19223382)
    ann = ga.parse_snpeff_ann(l1.info["ANN"])
    assert ann.gene_symbols == ("CLTCL1",)
    assert ann.most_severe == "intron_variant"
    assert "coding_sequence_variant" not in ann.terms


def test_snpeff_command_shape():
    cmd = ga.snpeff_command("in.vcf", "GRCh38.99", config="/x/snpEff.config", xmx="4g")
    assert cmd == ["snpEff", "-Xmx4g", "-noStats", "-c", "/x/snpEff.config", "GRCh38.99", "in.vcf"]
    assert "-c" not in ga.snpeff_command("in.vcf", "GRCh38.99")


def test_annotate_vcf_snpeff_offline_with_canned_output(tmp_path):
    import subprocess

    canned = SNPEFF_FIXTURE.read_text()
    calls = []

    def fake_run(cmd, **kw):
        calls.append(cmd)
        return subprocess.CompletedProcess(cmd, 0, stdout=canned, stderr="")

    out, tsv = tmp_path / "s.vcf", tmp_path / "s.tsv"
    stats = ga.annotate_vcf_snpeff(FIXTURE, out, "GRCh38.99", tsv_path=tsv, run=fake_run)
    assert calls and calls[0][0] == "snpEff" and "GRCh38.99" in calls[0]
    assert stats.n_records == 30 and stats.n_annotated == 30 and stats.n_unmatched == 0
    assert stats.n_with_gene == 23  # 7 intergenic loci in GRCh38.99 (VEP/release 116 finds genes at 1 of them)

    headers, recs = ga.read_vcf(out)
    assert any(h.startswith("##INFO=<ID=CSQ,") for h in headers)
    by_pos = {r.pos: r for r in recs}
    assert by_pos[19223382].info["GENE"] == "CLTCL1" and by_pos[19223382].info["CSQ"] == "intron_variant"
    assert by_pos[49029650].info["CSQ"] == "intergenic_variant" and "GENE" not in by_pos[49029650].info
    assert all("ANN" not in r.info for r in recs)         # raw ANN not carried into the shared contract
    assert all("SVLEN" in r.info for r in recs)           # original INFO preserved
    assert tsv.read_text().count("\n") == 31


def test_annotate_vcf_snpeff_reports_nonzero_exit(tmp_path):
    import subprocess

    def fake_run(cmd, **kw):
        return subprocess.CompletedProcess(cmd, 1, stdout="", stderr="java.lang.OutOfMemoryError: Java heap space")

    with pytest.raises(RuntimeError, match="OutOfMemoryError"):
        ga.annotate_vcf_snpeff(FIXTURE, tmp_path / "o.vcf", "GRCh38.99", run=fake_run)


def test_vep_and_snpeff_backends_agree_after_severity_ranking(tmp_path):
    """Cross-backend check on the fixture: same INFO keys, and concordant most-severe term
    wherever snpEff's GRCh38.99 database contains the gene VEP (release 116) reported."""
    import subprocess

    fake_run = lambda cmd, **kw: subprocess.CompletedProcess(cmd, 0, stdout=SNPEFF_FIXTURE.read_text(), stderr="")  # noqa: E731
    ga.annotate_vcf_snpeff(FIXTURE, tmp_path / "s.vcf", "GRCh38.99", run=fake_run)
    _, snp = ga.read_vcf(tmp_path / "s.vcf")
    for r in snp:
        assert set(r.info) >= {"CSQ", "CSQ_NTX"}
    # loci where both databases contain the gene: identical top term expected
    expect_intron = {19223382: "CLTCL1", 42705164: "A4GALT", 42818644: "ARFGAP3", 45595784: "FBLN1", 50351083: "PPP6R2", 17224410: "ADA2"}
    by_pos = {r.pos: r for r in snp}
    for pos, gene in expect_intron.items():
        assert by_pos[pos].info["CSQ"] == "intron_variant", pos
        assert gene in by_pos[pos].info["GENE"], pos


@pytest.mark.skipif(os.environ.get("RTM_LIVE_SNPEFF") != "1", reason="set RTM_LIVE_SNPEFF=1 (and have snpEff on PATH) to run snpEff for real")
def test_live_snpeff_run(tmp_path):
    genome = os.environ.get("RTM_SNPEFF_GENOME", "GRCh38.99")
    stats = ga.annotate_vcf_snpeff(FIXTURE, tmp_path / "o.vcf", genome, config=os.environ.get("RTM_SNPEFF_CONFIG"))
    assert stats.n_annotated == 30
