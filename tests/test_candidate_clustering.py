import numpy as np
from retro_miner.candidate_loci import _assign_discordant_to_seeds

def test_assign_discordant_chunking_empty():
    assigned = _assign_discordant_to_seeds([], [], [], 200)
    assert assigned.size == 0

def test_assign_discordant_contained():
    # position inside a seed interval -> distance 0, assigned to that seed
    assigned = _assign_discordant_to_seeds([100, 500], [200, 600], [150], 200)
    assert assigned.tolist() == [0]

def test_assign_discordant_nearest():
    # position 190: seed0=[100,200] dist 0 wins over seed1=[500,600]
    assigned = _assign_discordant_to_seeds([100, 500], [200, 600], [190], 200)
    assert assigned.tolist() == [0]

def test_assign_discordant_outside_radius():
    # position 900 is > 200 bp from every seed -> unassigned (-1)
    assigned = _assign_discordant_to_seeds([100, 500], [200, 600], [900], 200)
    assert assigned.tolist() == [-1]

def test_assign_discordant_chunking_large_scale():
    rng = np.random.default_rng(0)
    starts = np.sort(rng.integers(1000, 500000, size=5000))
    ends = starts + 200
    pos = np.sort(rng.integers(1000, 500000, size=10000))
    assigned = _assign_discordant_to_seeds(starts.tolist(), ends.tolist(), pos.tolist(), 200)
    assert assigned.size == pos.size
    assert assigned.max() < 5000
    # every valid assignment must be within radius of its nearest interval
    s = starts
    e = ends
    for i, p in enumerate(pos.tolist()):
        ai = int(assigned[i])
        if ai < 0:
            continue
        d = 0 if s[ai] <= p <= e[ai] else min(abs(p - s[ai]), abs(p - e[ai]))
        assert d <= 200
