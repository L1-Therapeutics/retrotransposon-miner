"""Tests for scripts/insertion_size_histograms.py.

The four loader tests exercise ``load_calls`` (parsing, sample attribution,
per-file metadata, missing optional fields). The remaining tests are the real
replacements for the four skipped placeholders: grouping/merge-bp, family
separation, single-sample collapse, filter-vs-group order, and sentinel
handling, plus plotting and the CLI score-column error contract.

All VCFs are written under pytest's ``tmp_path`` by the module-level
``write_vcf`` helper; no ``*.vcf`` is ever written into the repository (the
repo's .gitignore ignores them outside ``docs/examples/``).
"""

from __future__ import annotations

import csv
import importlib.util
import sys
from pathlib import Path

import pytest

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "insertion_size_histograms.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("insertion_size_histograms", SCRIPT_PATH)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["insertion_size_histograms"] = mod  # must precede exec_module
    spec.loader.exec_module(mod)
    return mod


# ---------------------------------------------------------------------------
# Minimal VCF builder (all files go under tmp_path)
# ---------------------------------------------------------------------------
_INFO_HEADER_LINES = [
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
    '##INFO=<ID=MEIFAMILY,Number=1,Type=String,Description="Consensus MEI family">',
    '##INFO=<ID=MEI_5P,Number=1,Type=Integer,Description="5-prime coordinate">',
    '##INFO=<ID=MEI_3P,Number=1,Type=Integer,Description="3-prime coordinate">',
    '##INFO=<ID=MEI_SPAN,Number=1,Type=Integer,Description="Full-length span">',
    '##INFO=<ID=L1TXGOLDSCORE,Number=1,Type=Float,Description="Classifier probability">',
    '##INFO=<ID=CALLTIER,Number=1,Type=String,Description="Call tier">',
]


def write_vcf(path, rows, samples=("S1",)):
    """Write a minimal VCF for the loader under test.

    ``rows`` is an iterable of dicts; recognised keys:
    ``chrom`` (default ``chr22``), ``pos`` (required), ``family`` (default
    ``ALU``), ``span``, ``five_p``, ``three_p``, ``score``, ``tier``,
    ``ref``/``alt``/``id``/``filter``, and ``extra_info`` (dict merged last).

    ``MEI_3P`` defaults to ``five_p + span - 1`` when not given. ``samples``
    controls the ``#CHROM`` sample columns; pass ``()`` to emit a sites-only VCF
    (the loader then falls back to the file stem as the sample name).
    """
    if isinstance(samples, str):
        samples = (samples,)

    lines = ["##fileformat=VCFv4.3", "##contig=<ID=chr22>", *_INFO_HEADER_LINES]
    header = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
    if samples:
        header += ["FORMAT", *samples]
    lines.append("\t".join(header))

    for row in rows:
        chrom = row.get("chrom", "chr22")
        pos = row["pos"]
        family = row.get("family", "ALU")
        # Default to a confident, tiered call; pass score=None / tier=None to
        # omit those INFO keys (used by the score-absence contract tests).
        score = row.get("score", 0.999)
        tier = row.get("tier", "none")
        info = {"SVTYPE": "INS", "MEIFAMILY": family}
        if row.get("five_p") is not None:
            info["MEI_5P"] = str(row["five_p"])
        if row.get("span") is not None:
            info["MEI_SPAN"] = str(row["span"])
        if row.get("three_p") is not None:
            info["MEI_3P"] = str(row["three_p"])
        elif row.get("span") is not None and row.get("five_p") is not None:
            info["MEI_3P"] = str(row["five_p"] + row["span"] - 1)
        if score is not None:
            info["L1TXGOLDSCORE"] = f"{float(score):.4f}"
        if tier is not None:
            info["CALLTIER"] = str(tier)
        info.update({str(k): str(v) for k, v in row.get("extra_info", {}).items()})

        cols = [
            chrom,
            str(pos),
            row.get("id", f"L1TX-{chrom}-{pos}-{family}"),
            row.get("ref", "N"),
            row.get("alt", f"<INS:ME:{family}>"),
            ".",
            row.get("filter", "PASS"),
            ";".join(f"{k}={v}" for k, v in info.items()),
        ]
        if samples:
            cols += ["GT:GQ", *(["./.:."] * len(samples))]
        lines.append("\t".join(cols))

    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


# ---------------------------------------------------------------------------
# Small test-side helpers
# ---------------------------------------------------------------------------
def _groups_for(mod, vcf_dir, merge_bp=50, min_score=0.997, group_first=False):
    """Mirror main()'s two orders without going through argparse."""
    calls, _ = mod.load_calls(vcf_dir)
    kept, sentinels = mod.split_sentinels(calls)
    if group_first:
        groups = mod.keep_unique_across_samples(kept, merge_bp)
        groups, _ = mod._filter_groups_by_confidence(groups, min_score, None)
    else:
        confident, _ = mod.select_confident_calls(kept, min_score=min_score)
        groups = mod.keep_unique_across_samples(confident, merge_bp)
    return groups, sentinels


def _read_tsv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _summary_totals(out_dir):
    """(union_total, private_total, {(family, site_set): row}) from size_summary.tsv."""
    rows = _read_tsv(Path(out_dir) / "size_summary.tsv")
    by_key = {(r["family"], r["site_set"]): r for r in rows}
    union = sum(int(r["n"]) for r in rows if r["site_set"] == "union")
    private = sum(int(r["n"]) for r in rows if r["site_set"] == "private")
    return union, private, by_key


def _run_main(mod, vcf_dir, out_dir, extra=()):
    argv = ["--vcf-dir", str(vcf_dir), "--out-dir", str(out_dir), *extra]
    return mod.main(argv)


@pytest.fixture
def fast_plots(monkeypatch):
    """Module under test with PNG rendering stubbed out.

    ``main`` still computes and writes every TSV, but the (slow) matplotlib
    savefig calls are replaced by a histogram-only stub. Tests that assert on
    the summary/dropped tables use this; the real rendering path stays covered
    by ``test_plot_size_histogram_writes_png_and_handles_empty_family``.
    """
    mod = _load_module()

    def _histogram_only(sizes, family, site_set, bin_width, size_field="MEI_SPAN",
                        edges=None, out_path=None, settings_text=""):
        if edges is None:
            edges = mod._bin_edges(sizes, bin_width)
        counts, edges = mod.np.histogram(sizes, bins=edges)
        return edges, counts

    monkeypatch.setattr(mod, "plot_size_histogram", _histogram_only)
    monkeypatch.setattr(mod, "plot_combined", lambda *args, **kwargs: None)
    return mod


# ---------------------------------------------------------------------------
# Loader tests (kept)
# ---------------------------------------------------------------------------
def test_load_calls_parses_fields_and_attributes_every_sample(tmp_path):
    mod = _load_module()
    write_vcf(
        tmp_path / "two_samples.vcf",
        [dict(pos=1000, family="ALU", span=150, five_p=5, three_p=154, score=0.9981, tier="high_conf_two_sided")],
        samples=("S1", "S2"),
    )
    calls, metas = mod.load_calls(tmp_path)

    assert len(calls) == 2
    assert {c.sample for c in calls} == {"S1", "S2"}
    for call in calls:
        assert call.chrom == "chr22"
        assert call.pos == 1000
        assert call.family == "ALU"
        assert call.size == 150
        assert call.mei_5p == 5
        assert call.tier == "high_conf_two_sided"
        assert call.score == pytest.approx(0.9981)
    assert len(metas) == 1
    assert metas[0].samples == ["S1", "S2"]


def test_load_calls_records_field_presence_per_file(tmp_path):
    mod = _load_module()
    write_vcf(tmp_path / "with_score.vcf", [dict(pos=1000, family="ALU", span=150, five_p=1, score=0.999, tier="none")])
    # Header declares L1TXGOLDSCORE but no record carries it (the gold-review
    # VCF shape): record-level detection must report it as absent.
    write_vcf(tmp_path / "header_only_score.vcf", [dict(pos=1000, family="ALU", span=150, five_p=1, score=None, tier=None)])

    _, metas = mod.load_calls(tmp_path)
    by_name = {Path(m.path).name: m for m in metas}
    assert by_name["with_score.vcf"].score_declared is True
    assert by_name["with_score.vcf"].tier_declared is True
    assert by_name["header_only_score.vcf"].score_declared is False
    assert by_name["header_only_score.vcf"].tier_declared is False
    assert by_name["header_only_score.vcf"].size_declared is True


def test_load_calls_uses_file_stem_when_no_sample_columns(tmp_path):
    mod = _load_module()
    write_vcf(tmp_path / "sites_only.vcf", [dict(pos=1000, family="ALU", span=150, five_p=1, score=0.999)], samples=())
    calls, metas = mod.load_calls(tmp_path)

    assert len(calls) == 1
    assert calls[0].sample == "sites_only"
    assert metas[0].samples == ["sites_only"]


def test_load_calls_missing_optional_fields_are_none(tmp_path):
    mod = _load_module()
    write_vcf(
        tmp_path / "sparse.vcf",
        [dict(pos=1000, family="ALU", span=None, five_p=None, score=None, tier=None)],
    )
    calls, _ = mod.load_calls(tmp_path)
    call = calls[0]

    assert call.size is None
    assert call.mei_5p is None
    assert call.score is None
    assert call.tier is None


# ---------------------------------------------------------------------------
# Grouping / merge-bp  (worked example A)
# ---------------------------------------------------------------------------
def test_grouping_uses_merge_bp_and_counts_private_sites(tmp_path, fast_plots):
    mod = fast_plots
    # S1: ALU chr22:1000.  S2: ALU chr22:1030 and ALU chr22:5000.
    write_vcf(tmp_path / "s1.vcf", [dict(pos=1000, family="ALU", span=150, five_p=1)], samples="S1")
    write_vcf(
        tmp_path / "s2.vcf",
        [dict(pos=1030, family="ALU", span=200, five_p=1), dict(pos=5000, family="ALU", span=300, five_p=1)],
        samples="S2",
    )

    # --merge-bp 50: 1000 and 1030 chain together, 5000 is separate.
    groups50, _ = _groups_for(mod, tmp_path, merge_bp=50)
    assert len(groups50) == 2
    assert mod.size_series(groups50)["ALU"] == [150, 300]  # rep at 1000 has size 150
    merged = next(g for g in groups50 if g.min_pos == 1000)
    assert merged.n_samples == 2
    assert merged.samples == ["S1", "S2"]
    assert len(mod.select_site_set(groups50, "private")) == 1

    # --merge-bp 20: nothing chains; all three calls are their own site.
    groups20, _ = _groups_for(mod, tmp_path, merge_bp=20)
    assert len(groups20) == 3
    assert len(mod.select_site_set(groups20, "private")) == 3

    # And the CLI agrees on the site counts.
    union50, private50, _ = _summary_totals(_main_out(mod, tmp_path, "bp50", ["--merge-bp", "50"]))
    union20, private20, _ = _summary_totals(_main_out(mod, tmp_path, "bp20", ["--merge-bp", "20"]))
    assert (union50, private50) == (2, 1)
    assert (union20, private20) == (3, 3)


def _main_out(mod, vcf_dir, name, extra):
    out_dir = Path(vcf_dir) / f"out_{name}"
    _run_main(mod, vcf_dir, out_dir, extra)
    return out_dir


# ---------------------------------------------------------------------------
# Family separation  (worked example B)
# ---------------------------------------------------------------------------
def test_family_separation_keeps_same_position_sites_apart(tmp_path, fast_plots):
    mod = fast_plots
    write_vcf(tmp_path / "s1.vcf", [dict(pos=1000, family="LINE1", span=150, five_p=1)], samples="S1")
    write_vcf(tmp_path / "s2.vcf", [dict(pos=1000, family="ALU", span=200, five_p=1)], samples="S2")

    groups, _ = _groups_for(mod, tmp_path, merge_bp=50)
    assert len(groups) == 2
    assert {g.family for g in groups} == {"ALU", "LINE1"}
    assert all(g.n_samples == 1 for g in groups)
    assert len(mod.select_site_set(groups, "private")) == 2

    union, private, by_key = _summary_totals(_main_out(mod, tmp_path, "family", []))
    assert (union, private) == (2, 2)
    assert by_key[("ALU", "union")]["n"] == "1"
    assert by_key[("LINE1", "union")]["n"] == "1"


# ---------------------------------------------------------------------------
# Same sample twice  (worked example C)
# ---------------------------------------------------------------------------
def test_calls_from_one_sample_collapse_to_one_private_site(tmp_path):
    mod = _load_module()
    write_vcf(
        tmp_path / "s1.vcf",
        [dict(pos=1000, family="ALU", span=150, five_p=1, score=0.999), dict(pos=1020, family="ALU", span=180, five_p=1, score=0.998)],
        samples="S1",
    )

    groups, _ = _groups_for(mod, tmp_path, merge_bp=50)
    assert len(groups) == 1
    site = groups[0]
    assert site.n_samples == 1
    assert site.samples == ["S1"]
    assert site.representative.pos == 1000  # highest score wins
    assert mod.size_series(groups)["ALU"] == [150]
    assert len(mod.select_site_set(groups, "private")) == 1


# ---------------------------------------------------------------------------
# Order matters  (worked example D)
# ---------------------------------------------------------------------------
def test_group_before_filter_changes_private_sites(tmp_path, fast_plots):
    mod = fast_plots
    # S1 low confidence (0.50, span 280); S2 high confidence (0.999, span 270).
    write_vcf(tmp_path / "s1.vcf", [dict(pos=1000, family="ALU", span=280, five_p=1, score=0.50)], samples="S1")
    write_vcf(tmp_path / "s2.vcf", [dict(pos=1010, family="ALU", span=270, five_p=1, score=0.999)], samples="S2")

    # filter-then-group (default): the weak S1 call is gone before grouping.
    filtered, _ = _groups_for(mod, tmp_path, merge_bp=50, min_score=0.9)
    assert len(filtered) == 1
    assert filtered[0].n_samples == 1
    assert filtered[0].representative.sample == "S2"
    assert mod.size_series(filtered)["ALU"] == [270]
    assert len(mod.select_site_set(filtered, "private")) == 1

    # group-then-filter: both calls group, so n_samples == 2 and it is not private.
    grouped, _ = _groups_for(mod, tmp_path, merge_bp=50, min_score=0.9, group_first=True)
    assert len(grouped) == 1
    assert grouped[0].n_samples == 2
    assert grouped[0].representative.sample == "S2"
    assert mod.size_series(grouped)["ALU"] == [270]
    assert len(mod.select_site_set(grouped, "private")) == 0

    # CLI: same numbers through main() with --min-score 0.9.
    out_ff = _main_out(mod, tmp_path, "ff", ["--min-score", "0.9"])
    out_gf = _main_out(mod, tmp_path, "gf", ["--min-score", "0.9", "--group-first"])
    assert _summary_totals(out_ff)[:2] == (1, 1)
    assert _summary_totals(out_gf)[:2] == (1, 0)
    assert _summary_totals(out_ff)[2][("ALU", "union")]["max"] == "270"
    assert _summary_totals(out_gf)[2][("ALU", "union")]["max"] == "270"


# ---------------------------------------------------------------------------
# Sentinel  (worked example E)
# ---------------------------------------------------------------------------
def test_sentinel_excluded_from_sizes_and_counted(tmp_path, fast_plots):
    mod = fast_plots
    write_vcf(
        tmp_path / "sentinel.vcf",
        [
            # sentinel: 5P/3P == -1, span 0
            dict(pos=1000, family="ALU", five_p=-1, three_p=-1, span=0, score=0.999),
            # usable: span 150
            dict(pos=2000, family="ALU", five_p=1, span=150, score=0.999),
        ],
        samples="S1",
    )

    calls, _ = mod.load_calls(tmp_path)
    kept, sentinels = mod.split_sentinels(calls)

    assert len(kept) == 1
    assert len(sentinels) == 1
    assert sentinels[0].call.pos == 1000
    assert sentinels[0].call.mei_5p == -1

    # Counts per family and per sample.
    assert sum(1 for d in sentinels if d.call.family == "ALU") == 1
    assert sum(1 for d in sentinels if d.call.sample == "S1") == 1

    groups, _ = _groups_for(mod, tmp_path, merge_bp=50)
    sizes = mod.size_series(groups)["ALU"]
    assert sizes == [150]
    assert all(size > 0 for size in sizes)  # never a bar at 0

    # CLI wiring: dropped_sentinel reported, and the dropped record is listed.
    out_dir = _main_out(mod, tmp_path, "sentinel", [])
    _, _, by_key = _summary_totals(out_dir)
    assert by_key[("ALU", "union")]["n"] == "1"
    assert by_key[("ALU", "union")]["dropped_sentinel"] == "1"
    dropped = _read_tsv(Path(out_dir) / "dropped_records.tsv")
    assert len(dropped) == 1
    assert dropped[0]["sample"] == "S1"
    assert dropped[0]["family"] == "ALU"
    assert dropped[0]["stage"] == "sentinel"
    assert "MEI_5P == -1" in dropped[0]["reason"]


# ---------------------------------------------------------------------------
# Representative determinism
# ---------------------------------------------------------------------------
def test_representative_ties_use_sample_then_position(tmp_path):
    mod = _load_module()

    # Same score across samples -> smallest sample name wins.
    case1 = tmp_path / "case1"
    case1.mkdir()
    write_vcf(case1 / "z.vcf", [dict(pos=1000, family="ALU", span=100, five_p=1, score=0.999)], samples="zzz")
    write_vcf(case1 / "a.vcf", [dict(pos=1030, family="ALU", span=200, five_p=1, score=0.999)], samples="aaa")
    group, _ = _groups_for(mod, case1, merge_bp=50)
    assert len(group) == 1
    assert group[0].n_samples == 2
    assert group[0].representative.sample == "aaa"

    # Same score and same sample -> smallest POS wins within the merged group.
    case2 = tmp_path / "case2"
    case2.mkdir()
    write_vcf(
        case2 / "s.vcf",
        [dict(pos=1020, family="ALU", span=180, five_p=1, score=0.999), dict(pos=1000, family="ALU", span=150, five_p=1, score=0.999)],
        samples="same",
    )
    group2, _ = _groups_for(mod, case2, merge_bp=50)
    assert len(group2) == 1
    assert group2[0].n_samples == 1
    assert group2[0].representative.pos == 1000


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def test_plot_size_histogram_writes_png_and_handles_empty_family(tmp_path):
    mod = _load_module()

    out = tmp_path / "alu_union.png"
    edges, counts = mod.plot_size_histogram(
        [50, 60, 70], family="ALU", site_set="union", bin_width=10,
        size_field="MEI_SPAN", out_path=out,
    )
    assert out.is_file() and out.stat().st_size > 0
    assert int(counts.sum()) == 3
    assert int(counts[0]) == 0  # nothing plotted in the 0-10 bin
    assert list(edges[:2]) == [0.0, 10.0]

    empty = tmp_path / "sva_private.png"
    edges_empty, counts_empty = mod.plot_size_histogram(
        [], family="SVA", site_set="private", bin_width=100,
        size_field="MEI_SPAN", out_path=empty,
    )
    assert empty.is_file() and empty.stat().st_size > 0
    assert int(counts_empty.sum()) == 0
    assert list(edges_empty) == [0.0, 100.0]

    combined = tmp_path / "combined.png"
    panels = {
        ("union", "ALU"): ([50, 60], edges),
        ("union", "SVA"): ([], edges_empty),
        ("private", "ALU"): ([], edges),
        ("private", "SVA"): ([], edges_empty),
    }
    mod.plot_combined(panels, ("ALU", "SVA"), ("union", "private"), "MEI_SPAN", combined, "settings")
    assert combined.is_file() and combined.stat().st_size > 0


def test_main_writes_all_documented_outputs_end_to_end(tmp_path):
    mod = _load_module()
    write_vcf(tmp_path / "s1.vcf", [dict(pos=1000, family="ALU", span=50, five_p=1)], samples="S1")
    out_dir = tmp_path / "out"

    assert _run_main(mod, tmp_path, out_dir) == 0
    for name in (
        "size_summary.tsv",
        "bin_counts.tsv",
        "dropped_records.tsv",
        "size_histograms_combined.png",
        "size_hist_ALU_union.png",
        "size_hist_ALU_private.png",
        "size_hist_LINE1_union.png",
        "size_hist_LINE1_private.png",
        "size_hist_SVA_union.png",
        "size_hist_SVA_private.png",
    ):
        artifact = out_dir / name
        assert artifact.is_file(), name
        assert artifact.stat().st_size > 0, name


# ---------------------------------------------------------------------------
# CLI confidence-column contract
# ---------------------------------------------------------------------------
def test_main_errors_when_score_column_absent(tmp_path):
    mod = _load_module()
    write_vcf(tmp_path / "no_score.vcf", [dict(pos=1000, family="ALU", span=150, five_p=1, score=None, tier=None)])

    with pytest.raises(SystemExit) as excinfo:
        _run_main(mod, tmp_path, tmp_path / "out")
    message = str(excinfo.value)
    assert "L1TXGOLDSCORE is absent" in message
    assert "--no-score-filter" in message
    assert "--tier" in message


def test_main_no_score_filter_and_tier_avoid_the_error(tmp_path, fast_plots):
    mod = fast_plots
    write_vcf(tmp_path / "no_score.vcf", [dict(pos=1000, family="ALU", span=150, five_p=1, score=None, tier="none")])

    assert _run_main(mod, tmp_path, tmp_path / "out_noscore", ["--no-score-filter"]) == 0
    assert _run_main(mod, tmp_path, tmp_path / "out_tier", ["--tier", "none"]) == 0
