import gzip

import numpy as np
import pytest

from exon_distance_annotate import (
    INFO_HEADERS,
    annotate_info,
    annotate_vcf,
    exon_distance,
    gene_strand_at,
    insert_headers,
    load_gencode,
    merge_intervals,
)

EX = np.array([[1000, 1100], [5000, 5200]], dtype=np.int64)

GTF = """\
chr1\tH\tgene\t1\t9000\t.\t+\t.\tgene_id "G1";
chr1\tH\texon\t1001\t1100\t.\t+\t.\tgene_id "G1";
chr1\tH\texon\t5001\t5200\t.\t+\t.\tgene_id "G1";
chr2\tH\tgene\t1\t9000\t.\t-\t.\tgene_id "G2";
chr2\tH\texon\t1001\t1100\t.\t-\t.\tgene_id "G2";
"""

VCF = """\
##fileformat=VCFv4.3
##INFO=<ID=MEIFAMILY,Number=1,Type=String,Description="f">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t900\t.\tN\t<INS>\t.\tPASS\tMEIFAMILY=ALU
chr1\t1050\t.\tN\t<INS>\t.\tPASS\tMEIFAMILY=ALU
chr1\t1200\t.\tN\t<INS>\t.\tPASS\tMEIFAMILY=ALU
chr2\t900\t.\tN\t<INS>\t.\tPASS\tMEIFAMILY=ALU
chr3\t900\t.\tN\t<INS>\t.\tPASS\tMEIFAMILY=ALU
"""


@pytest.fixture
def paths(tmp_path):
    g = tmp_path / "t.gtf"
    g.write_text(GTF)
    v = tmp_path / "in.vcf"
    v.write_text(VCF)
    return g, v, tmp_path / "out.vcf"


def test_merge_intervals_joins_touching_exons():
    a = np.array([[10, 20], [20, 30], [100, 110]], dtype=np.int64)
    assert merge_intervals(a).tolist() == [[10, 30], [100, 110]]


def test_gene_strand_is_none_when_two_genes_disagree():
    assert gene_strand_at([(0, 100, "+"), (50, 150, "-")], 75) is None


def test_gene_strand_resolves_when_both_genes_agree():
    assert gene_strand_at([(0, 100, "+"), (50, 150, "+")], 75) == "+"


def test_inside_an_exon_gives_distance_zero_and_no_side():
    assert exon_distance(EX, 1050, "+") == (0, ".")


def test_upstream_on_a_plus_strand_gene():
    assert exon_distance(EX, 900, "+") == (100, "upstream")


def test_the_same_position_is_downstream_on_a_minus_strand_gene():
    assert exon_distance(EX, 900, "-") == (100, "downstream")


def test_downstream_on_a_plus_strand_gene():
    assert exon_distance(EX, 1200, "+") == (101, "downstream")


def test_a_gene_with_no_exons_returns_the_no_gene_sentinel():
    assert exon_distance(np.empty((0, 2), np.int64), 900, "+") == (-1, ".")


def test_annotate_info_appends_both_keys():
    assert annotate_info("A=1;B=2", 50, "upstream") == "A=1;B=2;EXON_DIST=50;EXON_SIDE=upstream"


def test_annotate_info_replaces_rather_than_duplicates_on_a_second_pass():
    once = annotate_info("A=1", 50, "upstream")
    twice = annotate_info(once, 7, "downstream")
    assert twice.count("EXON_DIST") == 1 and twice.count("EXON_SIDE") == 1
    assert twice == "A=1;EXON_DIST=7;EXON_SIDE=downstream"


def test_annotate_info_leaves_an_unrelated_key_with_a_matching_prefix_alone():
    out = annotate_info("EXON_DISTANCE_OLD=9", 1, "upstream")
    assert "EXON_DISTANCE_OLD=9" in out


def test_insert_headers_adds_both_lines_to_a_fresh_file():
    assert len(insert_headers(["##fileformat=VCFv4.3"])) == 2


def test_insert_headers_is_idempotent():
    assert insert_headers(["##fileformat=VCFv4.3", *INFO_HEADERS]) == []


def test_load_gencode_reads_genes_and_merges_exons(paths):
    g, _, _ = paths
    genes, exons = load_gencode(str(g))
    assert genes["chr1"] == [(0, 9000, "+")]
    assert exons["chr1"].tolist() == [[1000, 1100], [5000, 5200]]


def test_load_gencode_reads_a_gzipped_gtf(tmp_path):
    p = tmp_path / "t.gtf.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(GTF)
    genes, _ = load_gencode(str(p))
    assert genes["chr2"] == [(0, 9000, "-")]


def test_end_to_end_assigns_the_expected_value_to_every_record(paths):
    g, v, o = paths
    stats = annotate_vcf(str(v), str(g), str(o))
    lines = [x for x in o.read_text().splitlines() if not x.startswith("#")]
    info = [x.split("\t")[7] for x in lines]
    assert "EXON_DIST=100;EXON_SIDE=upstream" in info[0]      # chr1:900, + gene
    assert "EXON_DIST=0;EXON_SIDE=." in info[1]               # inside exon
    assert "EXON_DIST=101;EXON_SIDE=downstream" in info[2]    # chr1:1200, + gene
    assert "EXON_DIST=100;EXON_SIDE=downstream" in info[3]    # chr2:900, - gene
    assert "EXON_DIST=-1;EXON_SIDE=." in info[4]              # chr3, no gene
    assert stats == {"records": 5, "in_gene": 4, "in_exon": 1, "no_gene": 1}


def test_headers_are_written_before_the_chrom_line(paths):
    g, v, o = paths
    annotate_vcf(str(v), str(g), str(o))
    out = o.read_text().splitlines()
    chrom_i = next(i for i, x in enumerate(out) if x.startswith("#CHROM"))
    assert all(any(h in x for x in out[:chrom_i]) for h in INFO_HEADERS)


def test_running_twice_does_not_duplicate_headers_or_keys(paths):
    g, v, o = paths
    annotate_vcf(str(v), str(g), str(o))
    second = o.parent / "out2.vcf"
    annotate_vcf(str(o), str(g), str(second))
    text = second.read_text()
    assert text.count("ID=EXON_DIST") == 1
    assert all(x.split("\t")[7].count("EXON_DIST=") == 1
               for x in text.splitlines() if not x.startswith("#"))


def test_columns_other_than_info_are_untouched(paths):
    g, v, o = paths
    annotate_vcf(str(v), str(g), str(o))
    src = [x.split("\t")[:7] for x in VCF.splitlines() if not x.startswith("#")]
    got = [x.split("\t")[:7] for x in o.read_text().splitlines() if not x.startswith("#")]
    assert src == got
