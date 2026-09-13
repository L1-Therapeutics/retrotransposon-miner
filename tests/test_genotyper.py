import pytest
from retro_miner.genotyper import calculate_mei_genotype, GenotypeCall

def test_homozygous_reference_genotype():
    # 0 alt reads, 30 spanning reference reads
    res = calculate_mei_genotype(k_alt=0, k_ref=30)
    assert isinstance(res, GenotypeCall)
    assert res.genotype == "0/0"
    assert res.vaf == 0.0
    assert res.genotype_quality >= 60.0

def test_heterozygous_mei_genotype():
    # 15 alt reads, 15 spanning reference reads
    res = calculate_mei_genotype(k_alt=15, k_ref=15)
    assert res.genotype == "0/1"
    assert res.vaf == 0.50
    assert res.genotype_quality >= 50.0

def test_homozygous_non_reference_mei_genotype():
    # 29 alt reads, 1 spanning reference read
    res = calculate_mei_genotype(k_alt=29, k_ref=1)
    assert res.genotype == "1/1"
    assert res.vaf >= 0.95
    assert res.genotype_quality >= 50.0

def test_no_coverage_returns_missing():
    res = calculate_mei_genotype(k_alt=0, k_ref=0)
    assert res.genotype == "./."
    assert res.vaf == 0.0
    assert res.genotype_quality == 0.0
