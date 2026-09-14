"""Tests for the EM subclonal somatic MEI mixture model and pipeline tagging.

Validates that the 3-component Binomial mixture separates low-coverage
subclonal somatic retrotranspositions (VAF in [0.02, 0.25]) from germline
diploid insertions (VAF ~ 0.50) and chimeric PCR artifacts (VAF < 0.02),
per the Rodriguez-Martin / Evrony somatic mosaicism rationale.
"""

import tempfile
from pathlib import Path

import pytest

from retro_miner.pipeline_engine import run_scientific_mei_pipeline
from retro_miner.somatic_em import (
    ARTIFACT,
    GERMLINE,
    SOMATIC_SUBCLONAL,
    SomaticCall,
    SomaticEMClassifier,
)


def _germline_loci(count: int) -> list[tuple[int, int]]:
    return [(30, 60)] * count


def _somatic_loci(count: int) -> list[tuple[int, int]]:
    return [(4, 80)] * count


def _artifact_loci(count: int) -> list[tuple[int, int]]:
    return [(1, 100)] * count


class TestClassifyLocus:
    def test_somatic_subclone_classification(self):
        clf = SomaticEMClassifier()
        res = clf.classify_locus(k_alt=4, n=60)

        assert isinstance(res, SomaticCall)
        assert res.classification == SOMATIC_SUBCLONAL
        assert res.is_somatic is True
        assert res.somatic_posterior > 0.50

    def test_germline_heterozygous_classification(self):
        res = SomaticEMClassifier().classify_locus(k_alt=30, n=60)

        assert res.classification == GERMLINE
        assert res.is_somatic is False

    def test_artifact_classification(self):
        res = SomaticEMClassifier().classify_locus(k_alt=1, n=200)

        assert res.classification == ARTIFACT
        assert res.is_somatic is False


class TestSomaticEMClassifier:
    def test_somatic_subclonal_high_posterior(self):
        batch = _germline_loci(8) + _somatic_loci(8) + _artifact_loci(8)
        k_alt_list = [k for k, _ in batch]
        n_list = [n for _, n in batch]

        calls = SomaticEMClassifier().fit_predict(k_alt_list, n_list)

        assert len(calls) == len(batch)
        som_index = len(_germline_loci(8))
        call = calls[som_index]

        assert isinstance(call, SomaticCall)
        assert call.classification == SOMATIC_SUBCLONAL
        assert call.is_somatic is True
        assert call.somatic_posterior > 0.90
        assert call.subclone_vaf == pytest.approx(0.05, abs=0.02)

    def test_artifact_noise_filtered_at_vaf_001(self):
        batch = _germline_loci(5) + _artifact_loci(5)
        k_alt_list = [k for k, _ in batch]
        n_list = [n for _, n in batch]

        calls = SomaticEMClassifier().fit_predict(k_alt_list, n_list)

        art_index = len(_germline_loci(5))
        call = calls[art_index]

        assert call.classification == ARTIFACT
        assert call.is_somatic is False
        assert call.somatic_posterior < 0.05

    def test_germline_heterozygous_separation(self):
        batch = _germline_loci(5) + _somatic_loci(5) + _artifact_loci(5)
        k_alt_list = [k for k, _ in batch]
        n_list = [n for _, n in batch]

        calls = SomaticEMClassifier().fit_predict(k_alt_list, n_list)

        call = calls[0]

        assert call.classification == GERMLINE
        assert call.is_somatic is False
        assert call.somatic_posterior < 0.05

    def test_single_locus_spec_examples(self):
        somatic = SomaticEMClassifier().fit_predict([4], [80])[0]
        assert somatic.classification == SOMATIC_SUBCLONAL
        assert somatic.somatic_posterior > 0.90

        artifact = SomaticEMClassifier().fit_predict([1], [100])[0]
        assert artifact.classification == ARTIFACT
        assert artifact.somatic_posterior < 0.05

        germline = SomaticEMClassifier().fit_predict([30], [60])[0]
        assert germline.classification == GERMLINE
        assert germline.is_somatic is False

    def test_em_converges_and_learns_subclone_vaf(self):
        batch = _somatic_loci(10) + _germline_loci(10)
        k_alt_list = [k for k, _ in batch]
        n_list = [n for _, n in batch]

        classifier = SomaticEMClassifier()
        calls = classifier.fit_predict(k_alt_list, n_list)

        assert classifier.converged is True
        assert classifier.is_fitted is True
        assert classifier.n_iter > 0
        assert classifier.theta_som == pytest.approx(0.05, abs=0.02)
        assert any(call.is_somatic for call in calls)

    def test_predict_requires_fitted_model(self):
        classifier = SomaticEMClassifier()
        with pytest.raises(RuntimeError):
            classifier.predict(k_alt=4, n=80)

    def test_rejects_misaligned_or_empty_inputs(self):
        with pytest.raises(ValueError):
            SomaticEMClassifier().fit_predict([1, 2], [10])
        with pytest.raises(ValueError):
            SomaticEMClassifier().fit_predict([], [])


class TestPipelineTagging:
    def _raw_locus(self, chrom: str, pos: int, k_alt: int, k_ref: int) -> dict:
        return {
            "chrom": chrom,
            "pos": pos,
            "family": "L1",
            "left_clips": ["ATTGCAGCTAGCTAG"],
            "right_clips": ["ATTGCAGCTAGCTAGAAAAA"],
            "k_alt": k_alt,
            "k_ref": k_ref,
        }

    def test_pipeline_tags_somatic_and_germline_info_fields(self):
        raw_loci = [
            self._raw_locus("chr1", 100000, k_alt=4, k_ref=76),
            self._raw_locus("chr1", 200000, k_alt=30, k_ref=30),
            self._raw_locus("chr1", 300000, k_alt=1, k_ref=99),
        ]

        with tempfile.TemporaryDirectory() as tmpdir:
            vcf_out = run_scientific_mei_pipeline(raw_loci, Path(tmpdir) / "somatic.vcf")
            content = vcf_out.read_text()

        assert "SOMATIC=1" in content
        assert "SOMATIC_POST=0." in content
        assert "SUBCLONE_VAF=" in content

        somatic_line = next(line for line in content.splitlines() if "100000" in line)
        germline_line = next(line for line in content.splitlines() if "200000" in line)

        assert "SOMATIC=1" in somatic_line
        assert "SOMATIC_POST=0.9" in somatic_line
        assert "SOMATIC=0" in germline_line
        assert "##INFO=<ID=SOMATIC" in content
        assert "##INFO=<ID=SOMATIC_POST" in content
        assert "##INFO=<ID=SUBCLONE_VAF" in content