import pytest

from retro_miner.genotyper import calculate_mei_genotype, GenotypeCall


class TestCalculateMeiGenotype:
    """Bayesian MEI genotyper biological validation."""

    def test_homozygous_reference(self):
        call = calculate_mei_genotype(k_alt=0, k_ref=30)
        assert isinstance(call, GenotypeCall)
        assert call.genotype == "0/0"
        assert call.vaf == 0.0
        assert call.genotype_quality >= 60.0

    def test_heterozygous_integration(self):
        call = calculate_mei_genotype(k_alt=15, k_ref=15)
        assert call.genotype == "0/1"
        assert call.vaf == pytest.approx(0.50)
        assert call.genotype_quality >= 60.0

    def test_homozygous_non_reference(self):
        call = calculate_mei_genotype(k_alt=30, k_ref=1)
        assert call.genotype == "1/1"
        assert call.vaf == pytest.approx(30 / 31, rel=1e-3)
        assert call.genotype_quality >= 60.0

    def test_low_coverage_returns_missing(self):
        call = calculate_mei_genotype(k_alt=1, k_ref=1)
        assert call.genotype == "./."
        assert call.genotype_quality < 10.0

    def test_zero_depth_returns_missing(self):
        call = calculate_mei_genotype(k_alt=0, k_ref=0)
        assert call.genotype == "./."
        assert call.vaf == 0.0
        assert call.genotype_quality == 0.0
