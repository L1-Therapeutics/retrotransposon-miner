"""Unit tests for retro_miner.gene_annotation.

Set RTM_LIVE_SNPEFF=1 to additionally run snpEff against the chr22 fixture.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from retro_miner import gene_annotation as ga

FIXTURE = Path(__file__).parent / "data" / "chr22_mei.vcf"


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


# ----------------------------------------------------------------- snpEff

SNPEFF_FIXTURE = Path(__file__).parent / "data" / "chr22_mei.snpeff.vcf"

# Recorded snpEff 5.4c GRCh38.115 run of tests/data/chr22_mei.vcf with -noHgvs
# and every splice-site size set to 0. Gene symbols are comma-joined in the
# order parse_snpeff_ann emits them; empty means intergenic.
CHR22_SNPEFF: dict[int, tuple[str, str]] = {
    17224410: ("ADA2", "intron_variant"),
    17289460: ("ENSG00000308779", "intron_variant"),
    17567662: ("SLC25A18,ENSG00000286195", "intron_variant"),
    19223382: ("CLTCL1", "intron_variant"),
    19919244: ("TXNRD2", "intron_variant"),
    20075438: ("DGCR8,ENSG00000236540,Metazoa_SRP", "upstream_gene_variant"),
    20521112: ("MED15,ENSG00000236003", "intron_variant"),
    20595738: ("CCDC74BP1,ENSG00000307604", "upstream_gene_variant"),
    23938127: ("ENSG00000206090,ENSG00000290199", "intron_variant"),
    29236892: ("EMID1", "intron_variant"),
    29239707: ("EMID1", "intron_variant"),
    31355872: ("", "intergenic_variant"),
    31380162: ("", "intergenic_variant"),
    33124372: ("LINC01640", "intron_variant"),
    33132520: ("ENSG00000233632,LINC01640", "intron_variant"),
    34034616: ("LINC01643", "intron_variant"),
    36746494: ("CACNG2-DT", "intron_variant"),
    36752165: ("ENSG00000287773,CACNG2-DT", "intron_variant"),
    37529127: ("ENSG00000305681", "intron_variant"),
    40007330: ("FAM83F", "intron_variant"),
    41050312: ("Y_RNA", "downstream_gene_variant"),
    41051286: ("", "intergenic_variant"),
    41835230: ("SREBF2-AS1,CCDC134,SREBF2", "intron_variant"),
    42705164: ("A4GALT,CYB5R3", "intron_variant"),
    42818644: ("ARFGAP3", "intron_variant"),
    45595784: ("FBLN1", "intron_variant"),
    49029650: ("", "intergenic_variant"),
    49760480: ("RPL5P35", "upstream_gene_variant"),
    50351083: ("PPP6R2", "intron_variant"),
    50495066: ("MIOX", "downstream_gene_variant"),
}


def _chr22_calls(records: list[ga.VcfRecord]) -> dict[int, tuple[str, str]]:
    out: dict[int, tuple[str, str]] = {}
    for rec in records:
        ann = ga.parse_snpeff_ann(rec.info["ANN"]) if "ANN" in rec.info else None
        if ann is None:
            gene = rec.info.get("GENE") or ""
            csq = rec.info.get("CSQ") or ""
        else:
            gene = ",".join(ann.gene_symbols)
            csq = ann.most_severe
        out[rec.pos] = (gene, csq)
    return out


def test_snpeff_fixture_is_real_output_for_same_30_loci():
    _, source = ga.read_vcf(FIXTURE)
    _, snp_out = ga.read_vcf(SNPEFF_FIXTURE)
    assert len(snp_out) == 30
    assert [(r.chrom, r.pos, r.alt) for r in snp_out] == [(r.chrom, r.pos, r.alt) for r in source]
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


def test_chr22_snpeff_annotation_matches_recorded_run():
    """All 30 chr22 calls, from the recorded GRCh38.115 snpEff run."""
    _, snp_out = ga.read_vcf(SNPEFF_FIXTURE)
    assert _chr22_calls(snp_out) == CHR22_SNPEFF
    sva = next(r for r in snp_out if r.pos == 49029650)
    ann = ga.parse_snpeff_ann(sva.info["ANN"])
    assert "ENSG00000299733-ENSG00000307761" in sva.info["ANN"]
    assert ann.gene_symbols == () and ann.gene_ids == () and ann.n_transcripts == 0


def test_snpeff_svlen_locus_is_intronic_without_workaround():
    _, snp_out = ga.read_vcf(SNPEFF_FIXTURE)
    l1 = next(r for r in snp_out if r.pos == 19223382)
    ann = ga.parse_snpeff_ann(l1.info["ANN"])
    assert ann.gene_symbols == ("CLTCL1",)
    assert ann.most_severe == "intron_variant"
    assert "coding_sequence_variant" not in ann.terms


def test_snpeff_command_shape():
    cmd = ga.snpeff_command("in.vcf", "GRCh38.99", config="/x/snpEff.config", xmx="4g")
    assert cmd == [
        "snpEff", "-Xmx4g", "-noStats", "-noHgvs",
        "-spliceSiteSize", "0", "-spliceRegionExonSize", "0",
        "-spliceRegionIntronMin", "0", "-spliceRegionIntronMax", "0",
        "-c", "/x/snpEff.config", "GRCh38.99", "in.vcf",
    ]
    assert "-c" not in ga.snpeff_command("in.vcf", "GRCh38.99")


def test_annotate_vcf_snpeff_offline_with_canned_output(tmp_path):
    import subprocess

    canned = SNPEFF_FIXTURE.read_text()
    calls = []

    def fake_run(cmd, **kw):
        calls.append(cmd)
        return subprocess.CompletedProcess(cmd, 0, stdout=canned, stderr="")

    out, tsv = tmp_path / "s.vcf", tmp_path / "s.tsv"
    stats = ga.annotate_vcf(FIXTURE, out, "GRCh38.99", tsv_path=tsv, run=fake_run)
    assert calls and calls[0][0] == "snpEff" and "GRCh38.99" in calls[0]
    assert stats.n_records == 30 and stats.n_annotated == 30 and stats.n_unmatched == 0
    assert stats.n_with_gene == 26

    headers, recs = ga.read_vcf(out)
    assert any(h.startswith("##INFO=<ID=CSQ,") for h in headers)
    assert _chr22_calls(recs) == CHR22_SNPEFF
    assert all("ANN" not in r.info for r in recs)         # raw ANN not carried into the shared contract
    assert all("SVLEN" in r.info for r in recs)           # original INFO preserved
    assert tsv.read_text().count("\n") == 31


def test_annotate_vcf_snpeff_reports_nonzero_exit(tmp_path):
    import subprocess

    def fake_run(cmd, **kw):
        return subprocess.CompletedProcess(cmd, 1, stdout="", stderr="java.lang.OutOfMemoryError: Java heap space")

    with pytest.raises(RuntimeError, match="OutOfMemoryError"):
        ga.annotate_vcf(FIXTURE, tmp_path / "o.vcf", "GRCh38.99", run=fake_run)


@pytest.mark.skipif(os.environ.get("RTM_LIVE_SNPEFF") != "1", reason="set RTM_LIVE_SNPEFF=1 (and have snpEff on PATH) to run snpEff for real")
def test_live_snpeff_chr22_matches_recorded_run(tmp_path):
    genome = os.environ.get("RTM_SNPEFF_GENOME", "GRCh38.115")
    out = tmp_path / "o.vcf"
    stats = ga.annotate_vcf(
        FIXTURE, out, genome,
        config=os.environ.get("RTM_SNPEFF_CONFIG"),
        xmx=os.environ.get("RTM_SNPEFF_XMX", "3g"),
    )
    assert stats.n_annotated == 30 and stats.n_unmatched == 0
    _, recs = ga.read_vcf(out)
    assert _chr22_calls(recs) == CHR22_SNPEFF
