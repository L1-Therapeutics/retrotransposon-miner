"""Unit tests for retro_miner.gene_annotation.

Set RTM_LIVE_SNPEFF=1 to additionally run snpEff against the human chr22 and
mouse chr19 fixtures (needs the GRCh38.115 and GRCm39.115 databases).
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


def test_load_gene_strands_from_gtf_strips_version(tmp_path):
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(
        '1\thavana\tgene\t1\t10\t.\t-\t.\tgene_id "ENSG00000093072"; gene_version "20";\n'
        '1\thavana\texon\t1\t10\t.\t+\t.\tgene_id "ENSG00000093072";\n'
    )
    assert ga.load_gene_strands(gtf) == {"ENSG00000093072": "-"}
    ann = ga.GeneAnnotation(("ADA2",), ("ENSG00000093072.20",), "intron_variant", ("intron_variant",), 1)
    assert ann.info_fields(ga.load_gene_strands(gtf))["GENE_STRAND"] == "-"


def test_bundled_gene_strands_match_chr22_genes():
    strands = ga.bundled_gene_strands()
    assert strands["ENSG00000093072"] == "-"  # ADA2
    assert strands["ENSG00000182902"] == "+"  # SLC25A18
    assert strands["ENSG00000286195"] == "-"
    # Mouse Ensembl GRCm39.115 gene IDs are also bundled so mm39 VCFs resolve GENE_STRAND out of the box.
    assert strands["ENSMUSG00000104478"] == "+"
    assert strands["ENSMUSG00000101231"] == "-"


# ----------------------------------------------------------------- mouse (mm39 / GRCm39.115)

MOUSE_FIXTURE = Path(__file__).parent / "data" / "chr19_mouse_mei.vcf"
MOUSE_SNPEFF_FIXTURE = Path(__file__).parent / "data" / "chr19_mouse_mei.snpeff.vcf"

# 30 gold chr19 calls from s3://l1tx-data/mouse/results/hazen2016_b2_vs_thymus_chr19/
# candidate_loci.mei.gold_review.tsv (rtm export-vcf --reference-build mm39), and the
# recorded snpEff 5.4c GRCm39.115 run of them with the snpeff_command() flags.
# Values: (GENE, CSQ, GENE_STRAND) as annotate_vcf writes them.
CHR19_MOUSE_SNPEFF: dict[int, tuple[str, str, str]] = {
    3977507: ("Gm59305,Gm50327,Aldh3b1", "intron_variant", "-,-,-"),
    4775325: ("Gm59018,Sptbn2", "intron_variant", "-,+"),
    6107110: ("Fau,Znhit2,Mrpl49,Syvn1,Gm52998,Gm56670", "non_coding_transcript_exon_variant", "+,+,-,+,+,-"),
    6245682: ("Majin", "intron_variant", "+"),
    6401894: ("Map4k2,Gm14966", "non_coding_transcript_exon_variant", "+,-"),
    7262557: ("Mark2", "intron_variant", "-"),
    7838242: ("Slc22a27", "downstream_gene_variant", "-"),
    8407890: ("", "intergenic_variant", ""),
    8520633: ("Gm8189", "upstream_gene_variant", "-"),
    9125748: ("Gm9750", "intron_variant", "+"),
    9835247: ("Gm54995", "downstream_gene_variant", "+"),
    13473506: ("Or5b120", "upstream_gene_variant", "+"),
    17821268: ("", "intergenic_variant", ""),
    18843274: ("Trpm6", "intron_variant", "+"),
    20468640: ("Aldh1a1,C730002L08Rik", "intron_variant", "+,-"),
    22986035: ("C330002G04Rik", "non_coding_transcript_exon_variant", "-"),
    24072174: ("Tjp2,Gm41824", "3_prime_UTR_variant", "-,+"),
    28208847: ("", "intergenic_variant", ""),
    38776704: ("Noc3l,Plce1", "3_prime_UTR_variant", "-,+"),
    38810626: ("Noc3l", "upstream_gene_variant", "-"),
    38962677: ("Cyp2c52-ps,Hells", "upstream_gene_variant", "+,+"),
    40740574: ("Cc2d2b", "intron_variant", "+"),
    41908570: ("Pgam1,Exosc1", "downstream_gene_variant", "+,-"),
    43934594: ("", "intergenic_variant", ""),
    45169757: ("", "intergenic_variant", ""),
    45509186: ("Btrc", "non_coding_transcript_exon_variant", "+"),
    46745048: ("Cnnm2,Gm59200,Gm59203", "upstream_gene_variant", "+,+,-"),
    47550374: ("Gm59192", "upstream_gene_variant", "+"),
    47558300: ("Gm59192", "downstream_gene_variant", "+"),
    60509819: ("Cacul1", "downstream_gene_variant", "-"),
}


def _mouse_calls(records: list[ga.VcfRecord]) -> dict[int, tuple[str, str, str]]:
    return {r.pos: (r.info.get("GENE") or "", r.info["CSQ"], r.info.get("GENE_STRAND") or "") for r in records}


def test_mouse_snpeff_fixture_is_real_output_for_same_30_loci():
    headers, source = ga.read_vcf(MOUSE_FIXTURE)
    _, snp_out = ga.read_vcf(MOUSE_SNPEFF_FIXTURE)
    assert "##reference=mm39" in headers and "##assembly=GRCm39" in headers
    assert len(source) == 30
    assert [(r.chrom, r.pos, r.alt) for r in snp_out] == [(r.chrom, r.pos, r.alt) for r in source]
    assert {r.alt for r in source} == {"<INS:ME:B1>", "<INS:ME:B2>", "<INS:ME:LINE1>", "<INS:ME>"}
    assert all("ANN" in r.info for r in snp_out)
    assert all("ENSMUSG" in r.info["ANN"] and "ENSG0" not in r.info["ANN"] for r in snp_out)


def test_resolve_snpeff_genome_auto_selects_and_guards_species():
    from retro_miner.vcf_export import vcf_header_lines

    mouse_headers, _ = ga.read_vcf(MOUSE_FIXTURE)
    assert ga.vcf_header_species(mouse_headers) == "mouse"
    assert [h for h in vcf_header_lines("mm39") if h.startswith(("##reference=", "##assembly="))] == [
        h for h in mouse_headers if h.startswith(("##reference=", "##assembly="))
    ]
    human_headers = vcf_header_lines("hg38")
    assert ga.vcf_header_species(human_headers) == "human"
    assert ga.resolve_snpeff_genome(human_headers, None) == "GRCh38.99"
    with pytest.raises(ValueError, match="human reference"):
        ga.resolve_snpeff_genome(human_headers, "GRCm39.115")
    # The chr22 fixture predates ##reference headers; it keeps the human default.
    chr22_headers, _ = ga.read_vcf(FIXTURE)
    assert ga.vcf_header_species(chr22_headers) is None
    assert ga.resolve_snpeff_genome(chr22_headers, None) == "GRCh38.99"
    assert ga.resolve_snpeff_genome(mouse_headers, None) == "GRCm39.115"
    assert ga.resolve_snpeff_genome(mouse_headers, "GRCm39.115") == "GRCm39.115"
    with pytest.raises(ValueError, match="mouse reference"):
        ga.resolve_snpeff_genome(mouse_headers, "GRCh38.99")
    with pytest.raises(ValueError, match="mouse reference"):
        ga.resolve_snpeff_genome(mouse_headers, "GRCh38.115")


def test_annotate_vcf_mouse_replays_recorded_grcm39_run(tmp_path):
    import subprocess

    recorded = MOUSE_SNPEFF_FIXTURE.read_text()
    calls = []

    def replay_run(cmd, **kw):
        calls.append(cmd)
        return subprocess.CompletedProcess(cmd, 0, stdout=recorded, stderr="")

    out = tmp_path / "m.vcf"
    stats = ga.annotate_vcf(MOUSE_FIXTURE, out, run=replay_run)
    assert calls and "GRCm39.115" in calls[0]
    assert stats.n_records == 30 and stats.n_annotated == 30 and stats.n_unmatched == 0
    assert stats.n_with_gene == 25
    _, recs = ga.read_vcf(out)
    assert _mouse_calls(recs) == CHR19_MOUSE_SNPEFF
    sptbn2 = next(r for r in recs if r.pos == 4775325)
    assert sptbn2.info["GENEID"] == "ENSMUSG00000126466,ENSMUSG00000067889"


@pytest.mark.skipif(os.environ.get("RTM_LIVE_SNPEFF") != "1", reason="set RTM_LIVE_SNPEFF=1 (and have snpEff with GRCm39.115 on PATH) to run snpEff for real")
def test_live_snpeff_mouse_chr19_matches_recorded_run(tmp_path):
    out = tmp_path / "m.vcf"
    stats = ga.annotate_vcf(
        MOUSE_FIXTURE, out,
        config=os.environ.get("RTM_SNPEFF_CONFIG"),
        xmx=os.environ.get("RTM_SNPEFF_XMX", "4g"),
    )
    assert stats.n_annotated == 30 and stats.n_unmatched == 0
    _, recs = ga.read_vcf(out)
    assert _mouse_calls(recs) == CHR19_MOUSE_SNPEFF


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
    assert any(h.startswith("##INFO=<ID=GENE_STRAND,") for h in headers)
    assert _chr22_calls(recs) == CHR22_SNPEFF
    ada2 = next(r for r in recs if r.pos == 17224410)
    assert ada2.info["GENE_STRAND"] == "-"
    both = next(r for r in recs if r.pos == 17567662)
    assert both.info["GENE"] == "SLC25A18,ENSG00000286195"
    assert both.info["GENE_STRAND"] == "+,-"
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
