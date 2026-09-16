from retro_miner.en_cleavage_scorer import ElemCleavageScore, score_en_cleavage_site


def test_canonical_en_motif_scoring():
    res = score_en_cleavage_site("CGATTTTTAAGCTA")
    assert isinstance(res, ElemCleavageScore)
    assert res.is_canonical_en_site is True
    assert res.match_score >= 0.80

def test_non_canonical_motif_scoring():
    res = score_en_cleavage_site("CGACGCGCGCGC")
    assert res.is_canonical_en_site is False
    assert res.match_score < 0.50
