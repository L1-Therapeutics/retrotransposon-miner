from retro_miner.subfamily_voter import SubfamilyCall, classify_mei_subfamily


def test_l1hs_active_motif_classification():
    # Reads containing L1HS 3' UTR diagnostic ACA/GAG sequence
    reads = [
        "ACTGACTGAAAGACAGAGTTTTTTT",
        "GACAGAGACTGATTGCAGCTAGCT",
    ]
    call = classify_mei_subfamily(reads, family_hint="L1")
    assert isinstance(call, SubfamilyCall)
    assert call.top_subfamily == "L1HS"
    assert call.matching_kmers >= 2
    assert call.posterior_prob > 0.50

def test_aluya5_classification():
    reads = [
        "GCTAGCTTCCATCTACTGATCG",
        "CCATCTACTGCAGCTAGCTA",
    ]
    call = classify_mei_subfamily(reads, family_hint="Alu")
    assert call.top_subfamily == "AluYa5"
    assert call.matching_kmers >= 2

def test_empty_reads_returns_unknown():
    call = classify_mei_subfamily([], family_hint="L1")
    assert call.top_subfamily == "UNKNOWN"
    assert call.posterior_prob < 0.50
