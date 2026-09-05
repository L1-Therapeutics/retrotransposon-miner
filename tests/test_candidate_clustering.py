"""Unit tests for ``_assign_discordant_to_seeds``.

The production implementation in ``candidate_loci`` assigns each discordant
position to the *nearest* seed interval rather than materializing the dense
``n_seeds x n_positions`` distance matrix (which OOM'd on chr22).  It uses
monotonic left/right sweeps plus an active-seed heap and returns seed indices
(``-1`` for unassigned) parallel to the input ``positions``.

These tests validate correctness against a brute-force reference, plus the
edge cases and scale/memory invariants that motivated the memory-safe rewrite.
"""
from __future__ import annotations

import time

import numpy as np
import pytest

from retro_miner.candidate_loci import _assign_discordant_to_seeds


# ─────────────────────────────────────────────────────────────────────────────
# Brute-force reference implementation
# ─────────────────────────────────────────────────────────────────────────────


def _brute_assign(starts, ends, positions, radius):
    """Reference: O(n*m) scan matching the documented tie-break rule.

    Chosen seed minimizes the key ``(dist, start, end, index)`` where ``dist``
    is 0 when the position is contained in the closed interval ``[s, e]``,
    otherwise the distance to the nearest interval boundary.  Positions whose
    minimal distance exceeds ``radius`` are unassigned (``-1``).
    """
    starts = [int(s) for s in starts]
    ends = [int(e) for e in ends]
    radius = max(0, int(radius))
    out: list[int] = []
    for p in positions:
        p = int(p)
        best_key = None
        best_seed = -1
        for i, (s, e) in enumerate(zip(starts, ends)):
            if s <= p <= e:
                d = 0
            elif p < s:
                d = s - p
            else:
                d = p - e
            if d > radius:
                continue
            key = (d, s, e, i)
            if best_key is None or key < best_key:
                best_key = key
                best_seed = i
        out.append(best_seed)
    return np.asarray(out, dtype=np.int64)


# ─────────────────────────────────────────────────────────────────────────────
# Empty inputs
# ─────────────────────────────────────────────────────────────────────────────


class TestEmptyInputs:
    def test_no_seeds_no_positions(self):
        result = _assign_discordant_to_seeds([], [], [], 200)
        assert isinstance(result, np.ndarray)
        assert result.size == 0

    def test_no_seeds_with_positions_all_unassigned(self):
        result = _assign_discordant_to_seeds([], [], [100, 200, 300], 200)
        assert result.tolist() == [-1, -1, -1]

    def test_seeds_no_positions(self):
        result = _assign_discordant_to_seeds([100], [200], [], 200)
        assert result.size == 0

    def test_negative_radius_clamped_to_zero(self):
        # radius < 0 is clamped to 0: only perfectly contained positions assign
        result = _assign_discordant_to_seeds([100, 500], [200, 600], [150, 520], -5)
        assert result.tolist() == [0, 1]


# ─────────────────────────────────────────────────────────────────────────────
# Basic containment / nearest-seed behavior
# ─────────────────────────────────────────────────────────────────────────────


class TestBasicAssignment:
    def test_inside_interval_distance_zero(self):
        # position 150 is inside seed0=[100,200] -> dist 0
        result = _assign_discordant_to_seeds([100, 500], [200, 600], [150], 200)
        assert result.tolist() == [0]

    def test_nearest_interval_wins(self):
        # position 190: interval [100,200] contains it (dist 0)
        # interval [500,600] is far away -> seed0
        result = _assign_discordant_to_seeds([100, 500], [200, 600], [190], 200)
        assert result.tolist() == [0]

    def test_outside_radius_unassigned(self):
        # position 900 is > 200 bp from both seeds -> -1
        result = _assign_discordant_to_seeds([100, 500], [200, 600], [900], 200)
        assert result.tolist() == [-1]

    def test_exact_partial_overlap(self):
        # position 250 is closer to [100,200] (dist 50) than [500,600] (dist 250)
        result = _assign_discordant_to_seeds([100, 500], [200, 600], [250], 200)
        assert result.tolist() == [0]

    def test_boundary_left_of_first_seed(self):
        # position 50 at dist 50 from seed0 start=100, within radius
        result = _assign_discordant_to_seeds([100, 500], [200, 600], [50], 100)
        assert result.tolist() == [0]

    def test_boundary_right_of_last_seed(self):
        # position 700 at dist 100 from seed1 end=600, within radius
        result = _assign_discordant_to_seeds([100, 500], [200, 600], [700], 100)
        assert result.tolist() == [1]

    def test_single_seed_multiple_positions(self):
        result = _assign_discordant_to_seeds([100], [200], [150, 300, 50], 100)
        # 150 contained -> 0; 300 -> dist 100 (within radius) -> 0; 50 -> dist 50 -> 0
        assert result.tolist() == [0, 0, 0]


# ─────────────────────────────────────────────────────────────────────────────
# Distance boundary conditions
# ─────────────────────────────────────────────────────────────────────────────


class TestDistanceBoundaries:
    def test_exactly_max_dist_assigns(self):
        # position 300 is exactly 100 bp past seed0 end=200 (|s-p| == radius)
        result = _assign_discordant_to_seeds([100], [200], [300], 100)
        assert result.tolist() == [0]

    def test_just_beyond_max_dist_unassigned(self):
        # position 301 is 101 bp past seed0 end=200 (|s-p| == radius + 1) -> -1
        result = _assign_discordant_to_seeds([100], [200], [301], 100)
        assert result.tolist() == [-1]

    def test_exactly_max_dist_left(self):
        # position 0 is exactly 100 bp before seed start=100 -> assigns
        result = _assign_discordant_to_seeds([100], [200], [0], 100)
        assert result.tolist() == [0]

    def test_just_beyond_max_dist_left(self):
        result = _assign_discordant_to_seeds([100], [200], [-1], 100)
        assert result.tolist() == [-1]


# ─────────────────────────────────────────────────────────────────────────────
# Unsorted position inputs
# ─────────────────────────────────────────────────────────────────────────────


class TestUnsortedPositions:
    def test_unsorted_positions_keep_input_order(self):
        # Output must be parallel to the *input* positions, not sorted order.
        starts = [100, 500]
        ends = [200, 600]
        positions = [700, 150, 300, 900]  # deliberately unsorted
        result = _assign_discordant_to_seeds(starts, ends, positions, 100)
        # 700 -> seed1 (dist 100), 150 -> seed0 (contained), 300 -> within 100 of seed0,
        # 900 -> unassigned
        assert result.tolist() == [1, 0, 0, -1]

    def test_unsorted_positions_matches_brute_force(self):
        rng = np.random.default_rng(7)
        starts = np.sort(rng.integers(0, 100000, size=40))
        ends = starts + rng.integers(10, 500, size=40)
        positions = rng.integers(0, 120000, size=200)
        radius = 500
        result = _assign_discordant_to_seeds(starts.tolist(), ends.tolist(), positions.tolist(), radius)
        expected = _brute_assign(starts, ends, positions, radius)
        assert result.tolist() == expected.tolist()


# ─────────────────────────────────────────────────────────────────────────────
# Tie-breaking
# ─────────────────────────────────────────────────────────────────────────────


class TestTieBreaking:
    def test_contained_seed_chosen_over_boundary(self):
        # position 150 contained by seed0, also within radius of seed1 left boundary
        result = _assign_discordant_to_seeds([100, 140], [200, 400], [150], 60)
        # seed0 dist 0; seed1 dist 0 (140<=150<=400) -> tie on dist, smaller start wins
        assert result.tolist() == [0]

    def test_tie_on_equal_distance_smaller_start_wins(self):
        # position 250: dist to seed0 end=200 is 50; dist to seed1 start=300 is 50
        result = _assign_discordant_to_seeds([100, 300], [200, 400], [250], 100)
        # equal dist=50; seed0 has smaller start (100 < 300) -> seed0
        assert result.tolist() == [0]


# ─────────────────────────────────────────────────────────────────────────────
# Brute-force parity over randomized inputs
# ─────────────────────────────────────────────────────────────────────────────


class TestBruteForceParity:
    @pytest.mark.parametrize("seed", [0, 1, 2, 3, 4, 42, 1234])
    def test_randomized_parity(self, seed):
        rng = np.random.default_rng(seed)
        starts = np.sort(rng.integers(0, 500000, size=300))
        ends = starts + rng.integers(1, 1000, size=300)
        positions = rng.integers(0, 600000, size=2000)
        radius = int(rng.integers(0, 2000))
        result = _assign_discordant_to_seeds(starts.tolist(), ends.tolist(), positions.tolist(), radius)
        expected = _brute_assign(starts, ends, positions, radius)
        assert result.tolist() == expected.tolist()

    def test_dense_overlapping_seeds(self):
        # Many heavily overlapping seeds stress the active-seed heap expiry logic.
        rng = np.random.default_rng(9)
        starts = np.sort(rng.integers(0, 5000, size=500))  # tightly clustered
        ends = starts + rng.integers(1, 200, size=500)
        positions = rng.integers(0, 6000, size=3000)
        radius = 500
        result = _assign_discordant_to_seeds(starts.tolist(), ends.tolist(), positions.tolist(), radius)
        expected = _brute_assign(starts, ends, positions, radius)
        assert result.tolist() == expected.tolist()

    def test_contained_only_within_radius(self):
        result = _assign_discordant_to_seeds([100], [200], [150], 0)
        assert result.tolist() == [0]
        result2 = _assign_discordant_to_seeds([100], [200], [250], 0)
        assert result2.tolist() == [-1]


# ─────────────────────────────────────────────────────────────────────────────
# Scale / performance / memory invariants
# ─────────────────────────────────────────────────────────────────────────────


class TestScalePerformance:
    def test_large_scale_fast_and_memory_safe(self):
        # 50k seeds, 200k positions — must run < 1.0s and stay under 50 MB of
        # peak *additional* memory allocation (the old dense O(n*m) matrix was
        # ~50k * 200k * 8 bytes = 80 GB and OOM'd low-RAM hosts).
        rng = np.random.default_rng(2024)
        starts = np.sort(rng.integers(0, 50_000_000, size=50_000))
        ends = starts + rng.integers(100, 1000, size=50_000)
        positions = rng.integers(0, 55_000_000, size=200_000)
        radius = 2000

        baseline = _peak_rss_kb()
        t0 = time.perf_counter()
        result = _assign_discordant_to_seeds(starts.tolist(), ends.tolist(), positions.tolist(), radius)
        elapsed = time.perf_counter() - t0
        peak_kb = max(baseline, _peak_rss_kb()) - baseline

        assert elapsed < 1.0, f"assignment took {elapsed:.3f}s, expected < 1.0s"
        assert peak_kb < 50 * 1024, f"peak extra RSS {peak_kb} KiB exceeds 50 MiB"
        assert result.size == positions.size
        # Every non-negative assignment must be within radius of its seed.
        s, e = starts, ends
        for i in np.flatnonzero(result >= 0):
            ai = int(result[i])
            p = int(positions[i])
            d = 0 if s[ai] <= p <= e[ai] else min(abs(p - int(s[ai])), abs(p - int(e[ai])))
            assert d <= radius

    def test_large_scale_matches_brute_force(self):
        # A smaller large-ish random case cross-checked against the reference.
        rng = np.random.default_rng(55)
        starts = np.sort(rng.integers(0, 1_000_000, size=2000))
        ends = starts + rng.integers(10, 500, size=2000)
        positions = rng.integers(0, 1_200_000, size=5000)
        radius = 800
        result = _assign_discordant_to_seeds(starts.tolist(), ends.tolist(), positions.tolist(), radius)
        expected = _brute_assign(starts, ends, positions, radius)
        assert result.tolist() == expected.tolist()


def _peak_rss_kb() -> int:
    """Return current resident set size in KiB via /proc/self/statm (Linux)."""
    with open("/proc/self/statm") as fh:
        fields = fh.read().split()
    page_size = 4096
    return int(fields[1]) * page_size // 1024
