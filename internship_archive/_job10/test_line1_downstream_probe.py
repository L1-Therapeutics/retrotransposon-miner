import gzip

import numpy as np
import pandas as pd
import pytest

from line1_downstream_probe import (
    build,
    containing_gene,
    exon_rank,
    load_genes_and_exons,
    nearest_exon_index,
    test as window_test,
)

EX = np.array([[1000, 1100], [3000, 3100], [5000, 5200]], dtype=np.int64)

GTF = """\
chr1\tH\tgene\t1\t9000\t.\t+\t.\tgene_id "G1";
chr1\tH\texon\t1001\t1100\t.\t+\t.\tgene_id "G1";
chr1\tH\texon\t3001\t3100\t.\t+\t.\tgene_id "G1";
chr1\tH\texon\t5001\t5200\t.\t+\t.\tgene_id "G1";
chr2\tH\tgene\t1\t9000\t.\t-\t.\tgene_id "G2";
chr2\tH\texon\t1001\t1100\t.\t-\t.\tgene_id "G2";
chr2\tH\texon\t5001\t5200\t.\t-\t.\tgene_id "G2";
"""


@pytest.fixture
def ann(tmp_path):
    p = tmp_path / "t.gtf.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(GTF)
    return load_genes_and_exons(str(p))


def test_exon_rank_on_a_plus_strand_gene():
    assert exon_rank(EX, 0, "+") == "first"
    assert exon_rank(EX, 1, "+") == "internal"
    assert exon_rank(EX, 2, "+") == "last"


def test_exon_rank_reverses_on_a_minus_strand_gene():
    assert exon_rank(EX, 0, "-") == "last"
    assert exon_rank(EX, 2, "-") == "first"


def test_a_single_exon_gene_is_labelled_only():
    assert exon_rank(EX[:1], 0, "+") == "only"


def test_nearest_exon_index_picks_the_closest():
    assert nearest_exon_index(EX, 1150, 1250) == 0
    assert nearest_exon_index(EX, 2900, 2950) == 1
    assert nearest_exon_index(EX, 4000, 4100) in (1, 2)


def test_containing_gene_returns_none_when_strands_conflict():
    gl = [(0, 1000, "+", "A"), (500, 1500, "-", "B")]
    assert containing_gene(gl, 700) is None


def test_containing_gene_picks_the_longest_when_one_strand_overlaps():
    gl = [(0, 1000, "+", "A"), (0, 5000, "+", "B")]
    assert containing_gene(gl, 500)[3] == "B"


def _rep(chrom, start, end, strand):
    return pd.DataFrame({"chrom": [chrom], "start": [start], "end": [end],
                         "strand": [strand], "milli_div": [50]})


def test_build_labels_a_copy_downstream_of_the_first_exon(ann):
    genes, gex = ann
    out = build(_rep("chr1", 1150, 1250, "+"), genes, gex)
    r = out.iloc[0]
    assert r.side == "downstream" and r.exon_rank == "first" and r.sense


def test_build_flips_side_and_rank_on_a_minus_strand_gene(ann):
    genes, gex = ann
    out = build(_rep("chr2", 1150, 1250, "+"), genes, gex)
    r = out.iloc[0]
    assert r.side == "upstream" and r.exon_rank == "last" and not r.sense


def test_gene_fraction_is_measured_from_the_tss(ann):
    genes, gex = ann
    plus = build(_rep("chr1", 1150, 1250, "+"), genes, gex).iloc[0]
    minus = build(_rep("chr2", 1150, 1250, "+"), genes, gex).iloc[0]
    assert plus.gene_frac == pytest.approx(0.1333, abs=1e-3)
    assert minus.gene_frac == pytest.approx(0.8667, abs=1e-3)


def test_distance_to_tes_is_strand_aware(ann):
    genes, gex = ann
    plus = build(_rep("chr1", 1150, 1250, "+"), genes, gex).iloc[0]
    minus = build(_rep("chr2", 1150, 1250, "+"), genes, gex).iloc[0]
    assert plus.dist_to_tes == 7800
    assert minus.dist_to_tes == 1200


def test_a_copy_overlapping_an_exon_is_dropped(ann):
    genes, gex = ann
    assert len(build(_rep("chr1", 1050, 1080, "+"), genes, gex)) == 0


def test_a_copy_outside_any_gene_is_dropped(ann):
    genes, gex = ann
    assert len(build(_rep("chr3", 1150, 1250, "+"), genes, gex)) == 0


def _frame(near_anti, near_sense, far_anti, far_sense):
    rows = []
    for n, (d, snс) in enumerate([(50, False)] * near_anti + [(50, True)] * near_sense
                                 + [(9000, False)] * far_anti + [(9000, True)] * far_sense):
        rows.append({"side": "downstream", "dist_to_exon": d, "sense": snс})
    return pd.DataFrame(rows)


def test_window_test_recovers_a_planted_enrichment():
    r = window_test(_frame(90, 10, 500, 500), "downstream")
    assert r["anti_near"] == pytest.approx(0.90)
    assert r["OR"] > 1 and r["p"] < 1e-6


def test_window_test_is_null_when_fractions_match():
    r = window_test(_frame(50, 50, 500, 500), "downstream")
    assert r["OR"] == pytest.approx(1.0) and r["p"] > 0.9
