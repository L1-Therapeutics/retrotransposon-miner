"""Unit tests for retro_miner.local_assembly._poly_at_max_run."""
from __future__ import annotations

from retro_miner.local_assembly import _poly_at_max_run


class TestPolyAtMaxRun:
    def test_empty(self):
        assert _poly_at_max_run("") == 0

    def test_no_poly_base(self):
        assert _poly_at_max_run("CGCG") == 0

    def test_single_homopolymer(self):
        assert _poly_at_max_run("CCCCAAAAAGGGG") == 5

    def test_longest_of_multiple_runs(self):
        assert _poly_at_max_run("AAACCCAAATT") == 3

    def test_interspersed_a_t(self):
        assert _poly_at_max_run("ATAT") == 1

    def test_case_insensitive(self):
        assert _poly_at_max_run("ggggTTTTaaaa") == 4

    def test_long_poly_a(self):
        assert _poly_at_max_run("A" * 30) == 30
