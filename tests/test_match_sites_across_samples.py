"""Tests for cross-sample site matching in scripts/match_sites_across_samples.py.

Covers the pre-registered rule: exact match, +/-10 window edge, opposite-orientation
and different-host non-matches, caller-offset calibration, greedy tie-breaking
determinism, and empty-directory behavior. The script is loaded via importlib so its
top-level imports stay out of the package dependency graph (same convention as
test_locus_zoom_gif.py).
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "match_sites_across_samples.py"


def _load_match_module():
    import importlib.util as _iu

    spec = _iu.spec_from_file_location("match_sites_across_samples", SCRIPT_PATH)
    mod = _iu.module_from_spec(spec)
    sys.modules["match_sites_across_samples"] = mod  # must precede exec_module (dataclasses)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def mod():
    return _load_match_module()


# ---------------------------------------------------------------------------
# Fixture builders.
# ---------------------------------------------------------------------------

def write_rmsk(path: Path, records: list[tuple[str, int, int, str, str, str, str]]) -> None:
    """(chrom, start0, end0, strand, repName, repClass, repFamily) in rmsk.txt layout."""
    with path.open("w") as handle:
        for chrom, start0, end0, strand, name, cls, fam in records:
            handle.write(
                f"0\t1000\t0\t0\t0\t{chrom}\t{start0}\t{end0}\t0\t{strand}\t{name}\t{cls}\t{fam}\t0\t0\t0\n"
            )


def write_vcf(path: Path, records: list[dict]) -> None:
    """Minimal HG03086-schema callset: chrom, pos, MEIFAMILY, ORIENT, NESTED, TSD,
    KNOWNMEI, CALLTIER via INFO."""
    with path.open("w") as handle:
        handle.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for rec in records:
            info = ";".join(
                f"{key}={value}"
                for key, value in (
                    ("MEIFAMILY", rec.get("family", "ALU")),
                    ("ORIENT", rec.get("orient", "+")),
                    ("NESTED", "nested" if rec.get("nested", True) else "unnested"),
                    ("TSD", rec.get("tsd", "ACGT")),
                    ("KNOWNMEI", "True" if rec.get("known", False) else "False"),
                    ("CALLTIER", rec.get("tier", "none")),
                )
            )
            handle.write(f"{rec['chrom']}\t{rec['pos']}\t.\tN\t<INS>\t.\tPASS\t{info}\n")


def make_call(mod, sample="S1", chrom="chr1", pos=1000, family="ALU", orient="+",
              nested=True, tsd="ACGTACGT", known=False, tier="none"):
    return mod.Call(
        sample=sample, chrom=chrom, pos=pos, family=family, orientation=orient,
        nested=nested, tsd_seq=tsd, known_mei=known, call_tier=tier,
    )


def calls_to_sites(mod, calls, window=10):
    sites = mod.match_calls(calls, window_bp=window)
    return sorted(
        tuple(sorted((c.sample, c.pos) for c in members)) for members in sites
    )


# ---------------------------------------------------------------------------
# Matching rule.
# ---------------------------------------------------------------------------

class TestMatchingRule:
    def test_exact_match_same_site_two_samples(self, mod):
        calls = [make_call(mod, "S1", pos=1000), make_call(mod, "S2", pos=1000)]
        assert calls_to_sites(mod, calls) == [(("S1", 1000), ("S2", 1000))]

    def test_window_edge_ten_bp_matches(self, mod):
        calls = [make_call(mod, "S1", pos=1000), make_call(mod, "S2", pos=1010)]
        sites = calls_to_sites(mod, calls)
        assert len(sites) == 1 and len(sites[0]) == 2

    def test_window_edge_eleven_bp_does_not_match(self, mod):
        calls = [make_call(mod, "S1", pos=1000), make_call(mod, "S2", pos=1011)]
        sites = calls_to_sites(mod, calls)
        assert sites == [(("S1", 1000),), (("S2", 1011),)]

    def test_opposite_orientation_does_not_match(self, mod):
        calls = [make_call(mod, "S1", pos=1000, orient="+"), make_call(mod, "S2", pos=1000, orient="-")]
        sites = calls_to_sites(mod, calls)
        assert sites == [(("S1", 1000),), (("S2", 1000),)]

    def test_different_family_does_not_match(self, mod):
        calls = [make_call(mod, "S1", pos=1000, family="ALU"), make_call(mod, "S2", pos=1000, family="LINE1")]
        assert len(calls_to_sites(mod, calls)) == 2

    def test_different_host_does_not_match_within_window(self, mod, tmp_path):
        """Two calls 10 bp apart land in different same-family hosts (overlapping
        intervals; longest span wins) and therefore do not match."""
        rmsk = tmp_path / "rmsk.txt.gz"
        with_rmsk = rmsk.with_suffix("")
        write_rmsk(
            with_rmsk,
            [
                ("chr1", 100, 200, "+", "AluA", "SINE", "Alu"),   # len 100
                ("chr1", 190, 400, "+", "AluB", "SINE", "Alu"),   # len 210, wins at pos 195
            ],
        )
        import gzip as _gzip
        import shutil as _shutil

        with open(with_rmsk, "rb") as src, _gzip.open(rmsk, "wb") as dst:
            _shutil.copyfileobj(src, dst)

        calls = [
            make_call(mod, "S1", pos=186),   # 0-based 185: only AluA contains it
            make_call(mod, "S2", pos=196),   # 0-based 195: both contain; AluB longest
        ]
        assigned = mod.assign_hosts(calls, rmsk)
        for call in calls:
            call.host = assigned.get(call.chrom_pos_key)
        assert {c.host.name for c in calls} == {"AluA", "AluB"}
        assert calls_to_sites(mod, calls) == [(("S1", 186),), (("S2", 196),)]

    def test_unnested_calls_are_out_of_matching_scope(self, mod):
        """Matching is defined for nested calls; unnested rows never form sites."""
        calls = [make_call(mod, "S1", pos=1000, nested=False), make_call(mod, "S2", pos=1000, nested=False)]
        assert calls_to_sites(mod, calls) == []

    def test_no_transitive_chaining_sites_never_merge(self, mod):
        """Two already-formed sites stay separate even though cross pairs (S2-S3 at
        11 bp, S2-S4 at 15 bp... S2@1005 and S3@1016 are 11 bp apart) would allow
        merging: a call joins a site only via its representative, and two assigned
        calls are never fused into one site.
        """
        calls = [
            make_call(mod, "S1", pos=1000),
            make_call(mod, "S2", pos=1005),
            make_call(mod, "S3", pos=1016),
            make_call(mod, "S4", pos=1020),
        ]
        sites = calls_to_sites(mod, calls)
        assert sorted(sites) == [
            (("S1", 1000), ("S2", 1005)),
            (("S3", 1016), ("S4", 1020)),
        ]

    def test_sensitivity_windows_change_results_reportably(self, mod):
        calls = [make_call(mod, "S1", pos=1000), make_call(mod, "S2", pos=1012)]
        assert len(calls_to_sites(mod, calls, window=5)) == 2
        assert len(calls_to_sites(mod, calls, window=10)) == 2
        assert len(calls_to_sites(mod, calls, window=20)) == 1


# ---------------------------------------------------------------------------
# Greedy tie-breaking determinism.
# ---------------------------------------------------------------------------

class TestGreedyDeterminism:
    def test_tsd_overlap_decides_which_site_absorbs_an_ambiguous_call(self, mod):
        """S5@1008 is equidistant (8 bp) from both potential representatives; the
        greedy queue ranks the equal-distance pairs by TSD agreement, so S5's own
        TSD decides which call it seeds a site with. Identical inputs must give
        identical partitions."""
        base = [
            make_call(mod, "S1", pos=1000, tsd="AAAAAAAA"),
            make_call(mod, "S3", pos=1016, tsd="CCCCCCCC"),
        ]
        sites_x = mod.match_calls(base + [make_call(mod, "S5", pos=1008, tsd="AAAAAAAA")])
        sites_y = mod.match_calls(base + [make_call(mod, "S5", pos=1008, tsd="CCCCCCCC")])

        def members_of(sites):
            return next(sorted(c.sample for c in m) for m in sites if "S5" in {c.sample for c in m})

        assert members_of(sites_x) == ["S1", "S5"]          # S5's TSD favors S1's site
        # S5's TSD now favors S3: the greedy queue pairs S5 with S3 first, so S5
        # seeds a site that S1 then joins (S1 is 8 bp from the new representative) --
        # the partition itself changes, which is the tie-break taking effect.
        assert members_of(sites_y) == ["S1", "S3", "S5"]
        again_x = mod.match_calls(base + [make_call(mod, "S5", pos=1008, tsd="AAAAAAAA")])
        assert members_of(again_x) == members_of(sites_x)

    def test_repeated_runs_are_identical(self, mod):
        import copy

        calls = [
            make_call(mod, f"S{i}", pos=1000 + (i % 3) * 4, tsd="ACGT" * (1 + i % 4))
            for i in range(12)
        ]
        snapshot = [sorted((c.sample, c.pos) for c in m) for m in mod.match_calls(calls)]
        for _ in range(3):
            fresh = mod.match_calls(copy.deepcopy(calls))
            assert [sorted((c.sample, c.pos) for c in m) for m in fresh] == snapshot

    def test_tsd_overlap_counts_matching_prefix_only(self, mod):
        assert mod.tsd_overlap("ACGTACGT", "ACGTACGT") == 8
        assert mod.tsd_overlap("ACGTACGT", "ACGTTTTT") == 4
        assert mod.tsd_overlap("ACGT", "") == 0
        assert mod.tsd_overlap(".", "ACGT") == 0


# ---------------------------------------------------------------------------
# Caller-offset calibration.
# ---------------------------------------------------------------------------

class TestOffsetCalibration:
    def test_systematic_offset_estimated_and_applied(self, mod):
        """Other sample is shifted -7 bp relative to reference on 12 shared KNOWNMEI
        sites; median offset +7 must be estimated and align the nested site."""
        ref = [make_call(mod, "REF", pos=p, known=True) for p in range(10000, 10000 + 12 * 200, 200)]
        ref.append(make_call(mod, "REF", pos=40000))
        other = [make_call(mod, "S2", pos=p - 7, known=True) for p in range(10000, 10000 + 12 * 200, 200)]
        other.append(make_call(mod, "S2", pos=40000 - 7))
        estimate = mod.estimate_offset(ref, other, "S2")
        assert estimate.applied is True
        assert estimate.median_offset_bp == 7
        assert estimate.n_shared == 12

        all_calls = ref + other
        mod.apply_calibration("REF", all_calls, {"S2": estimate})
        assert all(c.calibr_pos == c.pos for c in all_calls if c.sample == "REF")
        assert all(c.calibr_pos == c.pos + 7 for c in all_calls if c.sample == "S2")

        sites = mod.match_calls(all_calls)
        assert len(sites) == 13
        assert all(len(members) == 2 for members in sites)

    def test_uncorrected_offset_prevents_match_and_calibration_restores_it(self, mod):
        """With a -14 bp systematic shift the raw positions sit outside +/-10, so the
        calls stay separate until calibration aligns them -- and the correction must
        be estimated from the KNOWNMEI sites, not assumed."""
        ref_known = [make_call(mod, "REF", pos=p, known=True) for p in range(10000, 10000 + 12 * 200, 200)]
        other_known = [make_call(mod, "S2", pos=p - 14, known=True) for p in range(10000, 10000 + 12 * 200, 200)]
        estimate = mod.estimate_offset(ref_known, other_known, "S2")
        assert estimate.applied is True and estimate.median_offset_bp == 14

        ref = [make_call(mod, "REF", pos=20000)]
        other = [make_call(mod, "S2", pos=20000 - 14)]
        for call in ref + other:
            call.calibr_pos = call.pos
        assert len(mod.match_calls(ref + other)) == 2  # 14 bp apart: no match

        mod.apply_calibration("REF", ref + other, {"S2": estimate})
        assert len(mod.match_calls(ref + other)) == 1  # aligned: one site

    def test_insufficient_shared_sites_not_applied(self, mod):
        ref = [make_call(mod, "REF", pos=p, known=True) for p in range(10000, 10000 + 9 * 200, 200)]
        other = [make_call(mod, "S2", pos=p - 7, known=True) for p in range(10000, 10000 + 9 * 200, 200)]
        estimate = mod.estimate_offset(ref, other, "S2")
        assert estimate.applied is False
        assert estimate.reason == "insufficient_knownmei_overlap"

    def test_unstable_offset_not_applied(self, mod):
        """Offsets alternate +/-40 bp: median ~0 but spread far beyond 3 bp MAD cap."""
        ref = [make_call(mod, "REF", pos=p, known=True) for p in range(10000, 10000 + 12 * 200, 200)]
        other = [
            make_call(mod, "S2", pos=p - (40 if i % 2 else -40), known=True)
            for i, p in enumerate(range(10000, 10000 + 12 * 200, 200))
        ]
        estimate = mod.estimate_offset(ref, other, "S2")
        assert estimate.applied is False
        assert estimate.reason == "offset_unstable_not_applied"


# ---------------------------------------------------------------------------
# Host assignment convention (audit_nested_tsd.py verbatim).
# ---------------------------------------------------------------------------

class TestHostAssignment:
    def test_containing_interval_selected_by_longest_span(self, mod, tmp_path):
        rmsk = tmp_path / "rmsk.txt"
        write_rmsk(
            rmsk,
            [
                ("chr1", 100, 200, "+", "AluA", "SINE", "Alu"),
                ("chr1", 150, 320, "+", "AluB", "SINE", "Alu"),
                ("chr1", 180, 190, "+", "AluC", "SINE", "Alu"),
            ],
        )
        calls = [make_call(mod, "S1", pos=185)]
        hosts = mod.assign_hosts(calls, rmsk)
        host = hosts[("chr1", 185)]
        assert (host.name, host.start0, host.end0, host.strand) == ("AluB", 150, 320, "+")
        assert host.offset_5p_0based == 34  # 0-based 184 - 150

    def test_minus_strand_offset_from_3p_end(self, mod, tmp_path):
        rmsk = tmp_path / "rmsk.txt"
        write_rmsk(rmsk, [("chr1", 100, 300, "-", "AluM", "SINE", "Alu")])
        hosts = mod.assign_hosts([make_call(mod, "S1", pos=150)], rmsk)
        assert hosts[("chr1", 150)].offset_5p_0based == 150  # pos0=149; 300-1-149 = 150 bases from the 5' (right) end

    def test_no_same_family_interval_leaves_host_unset(self, mod, tmp_path):
        rmsk = tmp_path / "rmsk.txt"
        write_rmsk(rmsk, [("chr1", 100, 200, "+", "AluA", "SINE", "Alu")])
        calls = [make_call(mod, "S1", pos=5000, family="LINE1")]
        assert mod.assign_hosts(calls, rmsk) == {}

    def test_family_for_rmsk_matches_audit_convention(self, mod):
        assert mod.family_for_rmsk("Alu", "SINE/Alu") == "ALU"
        assert mod.family_for_rmsk("L1", "LINE/L1") == "LINE1"
        assert mod.family_for_rmsk("L1MC", "LINE/L1") == "LINE1"
        assert mod.family_for_rmsk("SVA_A", "SVA") == "SVA"
        assert mod.family_for_rmsk("Simple_repeat", "Low_complexity") == ""


# ---------------------------------------------------------------------------
# End-to-end CLI behavior.
# ---------------------------------------------------------------------------

class TestEndToEnd:
    def _setup_two_samples(self, tmp_path):
        callset_dir = tmp_path / "callsets"
        callset_dir.mkdir()
        rmsk = tmp_path / "rmsk.txt.gz"
        import gzip as _gzip
        import shutil as _shutil

        plain = tmp_path / "rmsk.txt"
        write_rmsk(plain, [("chr1", 100, 2000, "+", "AluA", "SINE", "Alu")])
        with open(plain, "rb") as src, _gzip.open(rmsk, "wb") as dst:
            _shutil.copyfileobj(src, dst)

        write_vcf(
            callset_dir / "SAMPLE_A.vcf",
            [
                {"chrom": "chr1", "pos": 500, "known": True},
                {"chrom": "chr1", "pos": 1200, "tsd": "ACGTACGTAA"},
            ],
        )
        write_vcf(
            callset_dir / "SAMPLE_B.vcf",
            [
                {"chrom": "chr1", "pos": 493, "known": True},    # -7 bp systematic shift
                {"chrom": "chr1", "pos": 1193, "tsd": "ACGTACGTAA"},
            ],
        )
        return callset_dir, rmsk

    def test_full_run_outputs(self, mod, tmp_path):
        callset_dir, rmsk = self._setup_two_samples(tmp_path)
        outdir = tmp_path / "out"
        rc = mod.main([
            "--callset-dir", str(callset_dir),
            "--rmsk", str(rmsk),
            "--outdir", str(outdir),
            "--pooled-bcf", str(tmp_path / "absent.bcf"),
        ])
        assert rc == 0
        unique = pd.read_parquet(outdir / "unique_sites.parquet")
        assert len(unique) == 2
        # Both fixture sites are carried by both samples: A/B differ by the -7 bp
        # shift, which is inside the +/-10 window even without calibration here
        # (only one KNOWNMEI site is shared, below the calibration support floor).
        merged = unique[unique["n_carriers"] == 2]
        assert len(merged) == 2
        assert set(merged["samples"]) == {"SAMPLE_A|SAMPLE_B"}
        # Representative = earliest match_pos member: B@493 anchors the first site.
        assert str(merged.iloc[0]["representative_sample"]) == "SAMPLE_B"
        assert int(merged.iloc[0]["representative_pos"]) == 493
        assert all(bool(v) for v in merged["host_selection_concordance"])
        assert unique[unique["n_carriers"] == 1].empty

        curve = pd.read_csv(outdir / "accumulation_curve.csv")
        assert list(curve["genomes_added"]) == [1, 2]
        assert int(curve.iloc[1]["unique_all_families_sites"]) == 2

        audit = (outdir / "matching_audit.md").read_text()
        assert "SAMPLE_B" in audit and "+/-10" in audit and "offset" in audit
        manifest = (outdir / "matching_manifest.json").read_text()
        assert "SAMPLE_A" in manifest and "SAMPLE_B" in manifest
        assert "primary_window_bp" in manifest

    def test_empty_callset_directory_exits(self, mod, tmp_path):
        empty = tmp_path / "empty"
        empty.mkdir()
        with pytest.raises(SystemExit, match="no per-sample callsets found"):
            mod.main([
                "--callset-dir", str(empty),
                "--rmsk", str(tmp_path / "rmsk.txt"),
                "--outdir", str(tmp_path / "out"),
            ])

    def test_new_samples_are_picked_up_deterministically(self, mod, tmp_path):
        """A third callset dropped into the directory changes the result without any
        configuration change, and sample order follows the ID sort."""
        callset_dir, rmsk = self._setup_two_samples(tmp_path)
        write_vcf(
            callset_dir / "SAMPLE_C.vcf",
            [
                {"chrom": "chr1", "pos": 1200, "tsd": "ACGTACGTAA"},
                {"chrom": "chr1", "pos": 5000},   # private to SAMPLE_C
            ],
        )
        outdir = tmp_path / "out"
        mod.main([
            "--callset-dir", str(callset_dir),
            "--rmsk", str(rmsk),
            "--outdir", str(outdir),
            "--pooled-bcf", str(tmp_path / "absent.bcf"),
        ])
        unique = pd.read_parquet(outdir / "unique_sites.parquet")
        tri = unique[unique["samples"] == "SAMPLE_A|SAMPLE_B|SAMPLE_C"]
        assert len(tri) == 1 and int(tri.iloc[0]["n_carriers"]) == 3

        private = unique[unique["private"] == True]  # noqa: E712 -- pandas boolean column
        assert len(private) == 1
        assert str(private.iloc[0]["private_to_sample"]) == "SAMPLE_C"
        assert int(private.iloc[0]["n_carriers"]) == 1

        curve = pd.read_csv(outdir / "accumulation_curve.csv")
        assert list(curve["sample_added"]) == ["SAMPLE_A", "SAMPLE_B", "SAMPLE_C"]
        # A starts with 2 private calls; B matches both; C adds one new private site.
        assert int(curve.iloc[0]["unique_all_families_sites"]) == 2
        assert int(curve.iloc[1]["unique_all_families_sites"]) == 2
        assert int(curve.iloc[2]["unique_all_families_sites"]) == 3

    def test_vcf_gz_callsets_load(self, mod, tmp_path):
        callset_dir, rmsk = self._setup_two_samples(tmp_path)
        import gzip as _gzip

        plain = callset_dir / "SAMPLE_A.vcf"
        gz = callset_dir / "SAMPLE_A.vcf.gz"
        with open(plain, "rb") as src, _gzip.open(gz, "wb") as dst:
            import shutil as _shutil

            _shutil.copyfileobj(src, dst)
        plain.unlink()
        calls = mod.load_callset(gz, "SAMPLE_A")
        assert len(calls) == 2 and calls[0].family == "ALU" and calls[1].nested
