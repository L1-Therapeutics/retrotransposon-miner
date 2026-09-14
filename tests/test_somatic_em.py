import pytest
from retro_miner.somatic_em import SomaticEMClassifier, SomaticCall

def test_somatic_subclone_classification():
    clf = SomaticEMClassifier()
    # 4 alt reads out of 60 total reads -> VAF = 0.0667 (Somatic subclone)
    res = clf.classify_locus(k_alt=4, n=60)
    assert isinstance(res, SomaticCall)
    assert res.classification == "SOMATIC_SUBCLONAL"
    assert res.is_somatic is True
    assert res.somatic_posterior > 0.50

def test_germline_heterozygous_classification():
    clf = SomaticEMClassifier()
    # 30 alt reads out of 60 total reads -> VAF = 0.50 (Germline)
    res = clf.classify_locus(k_alt=30, n=60)
    assert res.classification == "GERMLINE"
    assert res.is_somatic is False

def test_artifact_classification():
    clf = SomaticEMClassifier()
    # 1 alt read out of 200 total reads -> VAF = 0.005 (Artifact)
    res = clf.classify_locus(k_alt=1, n=200)
    assert res.classification == "ARTIFACT"
    assert res.is_somatic is False
