import numpy as np
import pytest
from retro_miner.candidate_loci import _assign_discordant_to_seeds

def test_assign_discordant_chunking_empty():
    starts = np.array([], dtype=int)
    pos = np.array([], dtype=int)
    seeds, matched_pos = _assign_discordant_to_seeds(starts, pos, max_dist=200)
    assert len(seeds) == 0
    assert len(matched_pos) == 0

def test_assign_discordant_chunking_large_scale():
    starts = np.sort(np.random.randint(1000, 500000, size=5000))
    pos = np.sort(np.random.randint(1000, 500000, size=10000))
    seeds, matched_pos = _assign_discordant_to_seeds(starts, pos, max_dist=200)
    assert len(seeds) == len(matched_pos)