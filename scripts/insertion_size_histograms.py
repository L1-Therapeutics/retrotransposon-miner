#!/usr/bin/env python3
"""Insertion-size histograms from one or more caller VCFs.

This script reads every ``*.vcf`` in ``--vcf-dir``, picks one length-like INFO
field (``--size-field``, default ``MEI_SPAN``), collapses the same insertion
seen in several samples into one site, and writes per-family histograms plus
summary tables.

What ``MEI_SPAN`` is (see ``audit_notes.md``)
---------------------------------------------
``MEI_SPAN`` is **not** a read-resolved inserted length. It is the width of the
covered footprint on the full-length family consensus:
``MEI_SPAN == MEI_3P - MEI_5P + 1``, where 5P/3P are the extreme consensus
positions placed by supporting reads/read-pairs. Sentinel records
(``MEI_5P == -1`` or ``size <= 0``) carry no usable consensus footprint and are
dropped from the size distribution (they are counted and listed instead). The
x-axis label therefore says "covered footprint on the consensus", never
"insertion length". Swapping ``--size-field`` is possible but the audit found no
other exported VCF field that is a defensible inserted length (``SVLEN`` is
deliberately unset).

Site definition and single linkage
----------------------------------
Two calls are the "same insertion across genomes" when they share ``chrom`` and
``MEIFAMILY`` and are close on the reference. Calls are sorted by ``POS`` and a
new group starts when the gap to the *previous* call exceeds ``--merge-bp``.
Because the comparison is to the immediately preceding call, groups can **chain**
(single linkage): a chain of calls each within ``--merge-bp`` of the next forms a
single group even if its total width is much larger than ``--merge-bp``. A
group's ``n_samples`` is the number of *distinct* sample names in it; the
group's representative is the call with the highest ``L1TXGOLDSCORE``, ties
broken by smallest sample name, then smallest ``POS`` (fully deterministic).

Why filter/group order changes which sites survive
--------------------------------------------------
* ``filter-then-group`` (default): low-confidence calls are removed *before*
  grouping. They cannot join a group, so they neither raise ``n_samples`` nor
  bridge two nearby high-confidence calls. A site seen by one confident sample
  plus one weak sample stays **private**.
* ``group-then-filter`` (``--group-first``): grouping sees every non-sentinel
  call. The weak call can push a site's ``n_samples`` to 2 (removing it from the
  private set), and single-linkage chaining through a weak mediator can merge two
  otherwise-separate union sites into one. A group then survives only if its
  representative (highest score) passes the confidence gate; groups with no
  confident call are dropped.

Sample identity
---------------
A call is attributed to every sample column named in its file's ``#CHROM`` line
(these pipeline VCFs are per-sample, with a ``GT`` column that is intentionally
unset). If a file has no sample columns the file stem is used. Point
``--vcf-dir`` at a directory of per-sample VCFs to exploit cross-genome grouping.
"""

from __future__ import annotations

import argparse
import math
import statistics
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

# ---------------------------------------------------------------------------
# Defaults (every one is a CLI flag; print them on figures and the summary)
# ---------------------------------------------------------------------------
DEFAULT_SIZE_FIELD = "MEI_SPAN"
DEFAULT_MIN_SCORE = 0.997
DEFAULT_MERGE_BP = 50
DEFAULT_FAMILIES = ("ALU", "LINE1", "SVA")
DEFAULT_BIN_WIDTHS = {"ALU": 10, "LINE1": 250, "SVA": 100}
OBSERVED_TIERS = (
    "none",
    "provisional_one_sided",
    "high_conf_two_sided",
    "mei_with_complex",
)
SITE_SETS = ("union", "private")

_FAMILY_COLORS = {"ALU": "#4c72b0", "LINE1": "#dd8452", "SVA": "#55a868"}


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Call:
    """One VCF record attributed to one sample."""

    source_file: str
    sample: str
    chrom: str
    pos: int
    call_id: str
    family: str
    size: int | None
    mei_5p: int | None
    score: float | None
    tier: str | None


@dataclass
class FileMeta:
    """Per-file header facts used to validate the requested confidence gate."""

    path: str
    samples: list[str] = field(default_factory=list)
    score_declared: bool = False
    tier_declared: bool = False
    size_declared: bool = False
    five_declared: bool = False


@dataclass
class Dropped:
    """A record removed from the analysis, with the stage and the reason."""

    call: Call
    stage: str
    reason: str


@dataclass
class SiteGroup:
    """A cluster of calls at one reference locus for one family."""

    chrom: str
    family: str
    calls: list[Call]

    def __post_init__(self) -> None:
        self.calls = sorted(self.calls, key=lambda c: c.pos)
        self.samples = sorted({c.sample for c in self.calls})
        self.n_samples = len(self.samples)
        self.min_pos = self.calls[0].pos
        self.max_pos = self.calls[-1].pos
        self.representative = _representative(self.calls)


def _representative(calls: list[Call]) -> Call:
    """Highest score; ties -> smallest sample, then smallest POS."""

    def key(c: Call):
        score = c.score if c.score is not None else float("-inf")
        return (-score, c.sample, c.pos)

    return min(calls, key=key)


# ---------------------------------------------------------------------------
# Loading helpers
# ---------------------------------------------------------------------------
def _parse_info(raw: str) -> dict[str, str]:
    info: dict[str, str] = {}
    for item in raw.split(";"):
        if not item:
            continue
        key, sep, value = item.partition("=")
        info[key] = value if sep else ""
    return info


def _parse_int(value: str | None) -> int | None:
    if value is None or value == "":
        return None
    try:
        return int(value)
    except ValueError:
        return None


def _parse_float(value: str | None) -> float | None:
    if value is None or value == "":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def load_calls(
    vcf_dir: str | Path,
    size_field: str = DEFAULT_SIZE_FIELD,
    score_field: str = "L1TXGOLDSCORE",
    five_field: str = "MEI_5P",
    family_field: str = "MEIFAMILY",
    tier_field: str = "CALLTIER",
) -> tuple[list[Call], list[FileMeta]]:
    """Load every ``*.vcf`` in ``vcf_dir``.

    Returns ``(calls, file_metadata)``. Each record becomes one ``Call`` per
    sample column in its ``#CHROM`` header (file stem if there are none).
    """

    paths = sorted(Path(vcf_dir).glob("*.vcf"))
    if not paths:
        raise SystemExit(f"ERROR: no *.vcf files found in {vcf_dir!r}")

    calls: list[Call] = []
    metas: list[FileMeta] = []

    for path in paths:
        samples: list[str] = []
        header_fields: set[str] = set()
        record_fields: set[str] = set()

        for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
            if line.startswith("##"):
                if line.startswith("##INFO=<ID="):
                    header_fields.add(line.split("ID=", 1)[1].split(",", 1)[0].split(">", 1)[0])
                continue
            if line.startswith("#CHROM"):
                cols = line.rstrip("\n").split("\t")
                samples = [c for c in cols[9:] if c] or [path.stem]
                continue
            if not line.strip() or line.startswith("#"):
                continue

            cols = line.split("\t")
            if len(cols) < 8:
                continue
            info = _parse_info(cols[7])
            record_fields.update(info)
            family = info.get(family_field, "")
            size = _parse_int(info.get(size_field))
            mei_5p = _parse_int(info.get(five_field))
            score = _parse_float(info.get(score_field))
            tier = info.get(tier_field)
            pos = _parse_int(cols[1])
            if pos is None:
                continue
            for sample in samples:
                calls.append(
                    Call(
                        source_file=str(path),
                        sample=sample,
                        chrom=cols[0],
                        pos=pos,
                        call_id=cols[2],
                        family=family,
                        size=size,
                        mei_5p=mei_5p,
                        score=score,
                        tier=tier,
                    )
                )

        fields = header_fields | record_fields
        # The confidence gates read per-record values, so a header-only
        # declaration (e.g. the gold-review VCF, whose records never carry
        # L1TXGOLDSCORE) must NOT count as "present"; otherwise every record
        # would be silently dropped instead of raising the requested error.
        metas.append(
            FileMeta(
                path=str(path),
                samples=samples,
                score_declared=score_field in record_fields,
                tier_declared=tier_field in record_fields,
                size_declared=size_field in fields,
                five_declared=five_field in fields,
            )
        )

    return calls, metas


# ---------------------------------------------------------------------------
# Sentinels
# ---------------------------------------------------------------------------
def split_sentinels(
    calls: list[Call], size_field: str = DEFAULT_SIZE_FIELD
) -> tuple[list[Call], list[Dropped]]:
    """Separate sentinel records (no usable consensus footprint)."""

    kept: list[Call] = []
    dropped: list[Dropped] = []
    for call in calls:
        no_5p = call.mei_5p is None or call.mei_5p == -1
        if call.size is None:
            reason = f"{size_field} missing"
        elif call.size <= 0:
            reason = f"{size_field}={call.size} <= 0"
        else:
            reason = ""
        if no_5p:
            reason = "MEI_5P == -1" + (f" and {reason}" if reason else "")
        if reason:
            dropped.append(Dropped(call, "sentinel", reason))
        else:
            kept.append(call)
    return kept, dropped


# ---------------------------------------------------------------------------
# Confidence filter
# ---------------------------------------------------------------------------
def _confidence_reason(
    call: Call,
    min_score: float | None,
    allowed_tiers: set[str] | None,
    score_field: str,
    tier_field: str,
) -> str | None:
    """Return why ``call`` fails the confidence gate, or ``None`` if it passes."""

    if allowed_tiers is not None:
        if call.tier in allowed_tiers:
            return None
        return f"{tier_field}={call.tier!r} not in {sorted(allowed_tiers)}"
    if min_score is not None:
        if call.score is not None and call.score >= min_score:
            return None
        got = "missing" if call.score is None else f"{call.score:.4f}"
        return f"{score_field}={got} < {min_score}"
    return None


def select_confident_calls(
    calls: list[Call],
    min_score: float | None = DEFAULT_MIN_SCORE,
    tier: list[str] | None = None,
    score_field: str = "L1TXGOLDSCORE",
    tier_field: str = "CALLTIER",
) -> tuple[list[Call], list[Dropped]]:
    """Apply the confidence gate to raw calls.

    If ``tier`` is provided it **replaces** the score gate (this is what makes
    ``--tier`` useful on files that do not carry ``L1TXGOLDSCORE``). Otherwise
    calls must have ``score >= min_score``. ``min_score=None`` with no ``tier``
    disables the gate entirely.
    """

    allowed = set(tier) if tier else None
    kept: list[Call] = []
    dropped: list[Dropped] = []
    for call in calls:
        reason = _confidence_reason(call, min_score, allowed, score_field, tier_field)
        if reason is None:
            kept.append(call)
        else:
            dropped.append(Dropped(call, "confidence", reason))
    return kept, dropped


def _filter_groups_by_confidence(
    groups: list[SiteGroup],
    min_score: float | None,
    tier: list[str] | None,
    score_field: str = "L1TXGOLDSCORE",
    tier_field: str = "CALLTIER",
) -> tuple[list[SiteGroup], list[Dropped]]:
    """Group-first gate: a group survives iff its representative passes."""

    allowed = set(tier) if tier else None
    kept: list[SiteGroup] = []
    dropped: list[Dropped] = []
    for group in groups:
        rep = group.representative
        reason = _confidence_reason(rep, min_score, allowed, score_field, tier_field)
        if reason is None:
            kept.append(group)
            continue
        for call in group.calls:
            dropped.append(Dropped(call, "confidence(group-first)", reason))
    return kept, dropped


# ---------------------------------------------------------------------------
# Grouping
# ---------------------------------------------------------------------------
def keep_unique_across_samples(
    calls: list[Call],
    merge_bp: int = DEFAULT_MERGE_BP,
    families: tuple[str, ...] = DEFAULT_FAMILIES,
) -> list[SiteGroup]:
    """Collapse calls into one group per distinct insertion site.

    Groups share ``chrom`` and ``MEIFAMILY``; calls are sorted by ``POS`` and a
    new group starts when the gap to the previous call exceeds ``merge_bp``
    (single linkage, so chains are allowed).
    """

    wanted = set(families)
    buckets: dict[tuple[str, str], list[Call]] = defaultdict(list)
    for call in calls:
        if call.family in wanted:
            buckets[(call.chrom, call.family)].append(call)

    groups: list[SiteGroup] = []
    for (chrom, family), members in buckets.items():
        members.sort(key=lambda c: c.pos)
        current: list[Call] = [members[0]]
        for previous, call in zip(members, members[1:]):
            if call.pos - previous.pos > merge_bp:
                groups.append(SiteGroup(chrom, family, current))
                current = [call]
            else:
                current.append(call)
        groups.append(SiteGroup(chrom, family, current))

    groups.sort(key=lambda g: (g.chrom, g.min_pos, g.family))
    return groups


def select_site_set(groups: list[SiteGroup], site_set: str) -> list[SiteGroup]:
    """Return the requested site set: ``union`` (all) or ``private`` (n==1)."""

    if site_set == "union":
        return list(groups)
    if site_set == "private":
        return [g for g in groups if g.n_samples == 1]
    raise ValueError(f"unknown site_set {site_set!r}")


def size_series(
    groups: list[SiteGroup],
    families: tuple[str, ...] = DEFAULT_FAMILIES,
) -> dict[str, list[int]]:
    """Per-family list of representative sizes for the supplied site groups."""

    series: dict[str, list[int]] = {f: [] for f in families}
    for group in groups:
        if group.family not in series:
            continue
        size = group.representative.size
        if size is None:
            continue
        series[group.family].append(int(size))
    return series


# ---------------------------------------------------------------------------
# Statistics, bins, plotting
# ---------------------------------------------------------------------------
def _bin_edges(sizes: list[int], bin_width: int) -> np.ndarray:
    """Edges from 0 to the family maximum rounded up to a whole bin."""

    top = 0
    if sizes:
        top = math.ceil(max(sizes) / bin_width) * bin_width
    if top <= 0:
        top = bin_width
    return np.arange(0, top + bin_width, bin_width, dtype=float)


def _stats(sizes: list[int]) -> dict[str, float | int | None]:
    if not sizes:
        return {"n": 0, "min": None, "median": None, "max": None}
    return {
        "n": len(sizes),
        "min": int(min(sizes)),
        "median": float(statistics.median(sizes)),
        "max": int(max(sizes)),
    }


def _draw_panel(
    ax,
    sizes: list[int],
    edges: np.ndarray,
    family: str,
    site_set: str,
    size_field: str,
) -> None:
    color = _FAMILY_COLORS.get(family, "#8172b3")
    ax.hist(sizes, bins=edges, color=color, edgecolor="white", linewidth=0.4)
    ax.set_xscale("linear")
    ax.set_xlim(edges[0], edges[-1])
    ax.set_title(f"{family} — {site_set}")
    ax.set_xlabel(f"{size_field} — covered footprint on the consensus (bp)")
    ax.set_ylabel("number of sites")
    if sizes:
        label = f"n={len(sizes)}\nmedian={statistics.median(sizes):g} bp"
    else:
        label = "n=0"
    ax.text(
        0.97,
        0.95,
        label,
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=9,
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.8, "edgecolor": "#cccccc"},
    )


def plot_size_histogram(
    sizes: list[int],
    family: str,
    site_set: str,
    bin_width: int,
    size_field: str = DEFAULT_SIZE_FIELD,
    edges: np.ndarray | None = None,
    out_path: str | Path | None = None,
    settings_text: str = "",
) -> tuple[np.ndarray, np.ndarray]:
    """Draw and save a single family/site-set histogram.

    Always creates the figure with ``fig, ax = plt.subplots()`` and saves with
    ``fig.savefig`` (never module-level ``plt.plot`` / ``plt.savefig``). Returns
    ``(edges, counts)`` so the caller can reuse the exact bins.
    """

    if edges is None:
        edges = _bin_edges(sizes, bin_width)
    counts, edges = np.histogram(sizes, bins=edges)

    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    _draw_panel(ax, sizes, edges, family, site_set, size_field)
    if settings_text:
        fig.text(0.005, 0.005, settings_text, fontsize=6, va="bottom", ha="left", color="#444444")
    fig.tight_layout(rect=(0.0, 0.07, 1.0, 1.0))
    if out_path is not None:
        fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return edges, counts


def plot_combined(
    panels: dict[tuple[str, str], tuple[list[int], np.ndarray]],
    families: tuple[str, ...],
    site_sets: tuple[str, ...],
    size_field: str,
    out_path: str | Path,
    settings_text: str,
) -> None:
    """Grid of all family x site-set panels in one figure, with shared footer."""

    nrows, ncols = len(site_sets), len(families)
    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(5.0 * ncols, 3.8 * nrows),
        squeeze=False,
    )
    for row, site_set in enumerate(site_sets):
        for col, family in enumerate(families):
            sizes, edges = panels[(site_set, family)]
            _draw_panel(axes[row][col], sizes, edges, family, site_set, size_field)
    fig.suptitle("Insertion-size (consensus footprint) histograms", fontsize=13)
    fig.text(0.005, 0.005, settings_text, fontsize=6, va="bottom", ha="left", color="#444444")
    fig.tight_layout(rect=(0.0, 0.045, 1.0, 0.96))
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# CLI + orchestration
# ---------------------------------------------------------------------------
def _parse_bin_widths(items: list[str] | None, families: tuple[str, ...]) -> dict[str, int]:
    widths = dict(DEFAULT_BIN_WIDTHS)
    for item in items or []:
        if "=" not in item:
            raise SystemExit(f"ERROR: --bin-width expects FAMILY=BP (got {item!r})")
        family, _, value = item.partition("=")
        family = family.strip()
        try:
            widths[family] = int(value)
        except ValueError:
            raise SystemExit(f"ERROR: --bin-width BP must be an integer (got {item!r})")
    for family in families:
        if family not in widths:
            raise SystemExit(
                f"ERROR: no bin width for family {family!r}; pass --bin-width {family}=<bp>"
            )
    return widths


def _confidence_label(args) -> str:
    if args.tier:
        return f"CALLTIER in {list(args.tier)}"
    if args.no_score_filter:
        return "none (--no-score-filter)"
    return f"L1TXGOLDSCORE >= {args.min_score:g}"


def _settings_lines(args, families: tuple[str, ...], bin_widths: dict[str, int]) -> list[str]:
    order = "group-then-filter" if args.group_first else "filter-then-group"
    return [
        f"vcf_dir={args.vcf_dir}",
        f"size_field={args.size_field}",
        f"confidence={_confidence_label(args)}",
        f"merge_bp={args.merge_bp}",
        f"order={order}",
        f"families={','.join(families)}",
        "bin_width=" + ",".join(f"{f}:{bin_widths[f]}" for f in families),
        f"out_dir={args.out_dir}",
    ]


def _write_rows(path: str | Path, header: list[str], rows: list[list[object]]) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\t".join(header) + "\n")
        for row in rows:
            handle.write("\t".join("" if v is None else str(v) for v in row) + "\n")


def _fmt(value: object) -> str:
    if value is None:
        return "NA"
    if isinstance(value, float):
        return f"{value:g}"
    return str(value)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Insertion-size histograms from caller VCFs (see module docstring).",
    )
    parser.add_argument("--vcf-dir", required=True, help="Directory scanned for *.vcf")
    parser.add_argument(
        "--out-dir",
        default="insertion_size_histograms_out",
        help="Where PNGs and TSVs are written (default: %(default)s)",
    )
    parser.add_argument(
        "--size-field",
        default=DEFAULT_SIZE_FIELD,
        help="INFO field used as the size (default: %(default)s)",
    )
    parser.add_argument(
        "--min-score",
        type=float,
        default=DEFAULT_MIN_SCORE,
        help="Minimum L1TXGOLDSCORE for the confidence gate (default: %(default)s)",
    )
    parser.add_argument(
        "--no-score-filter",
        action="store_true",
        help="Disable the L1TXGOLDSCORE gate (use when the field is absent)",
    )
    parser.add_argument(
        "--tier",
        action="append",
        metavar="CALLTIER",
        help=(
            "Filter by CALLTIER membership instead of L1TXGOLDSCORE (repeatable). "
            "Observed values: " + ", ".join(OBSERVED_TIERS)
        ),
    )
    parser.add_argument(
        "--merge-bp",
        type=int,
        default=DEFAULT_MERGE_BP,
        help="Max gap to the previous call for single-linkage grouping (default: %(default)s)",
    )
    parser.add_argument(
        "--group-first",
        action="store_true",
        help="Group before applying the confidence gate (default: filter, then group)",
    )
    parser.add_argument(
        "--families",
        default=",".join(DEFAULT_FAMILIES),
        help="Comma-separated families to plot (default: %(default)s)",
    )
    parser.add_argument(
        "--bin-width",
        action="append",
        metavar="FAMILY=BP",
        help=(
            "Per-family bin width in bp (repeatable). Defaults: "
            + ", ".join(f"{f}={w}" for f, w in DEFAULT_BIN_WIDTHS.items())
        ),
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    families = tuple(f.strip() for f in args.families.split(",") if f.strip())
    bin_widths = _parse_bin_widths(args.bin_width, families)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    calls, metas = load_calls(args.vcf_dir, size_field=args.size_field)

    # --- column-presence validation (before any filtering) ------------------
    if args.tier:
        missing_tier = [m.path for m in metas if not m.tier_declared]
        if missing_tier:
            raise SystemExit(
                "ERROR: CALLTIER is absent from: "
                + ", ".join(missing_tier)
                + "\n  Drop --tier, or filter with --no-score-filter / --min-score instead."
            )
    if not args.no_score_filter and not args.tier:
        missing_score = [m.path for m in metas if not m.score_declared]
        if missing_score:
            raise SystemExit(
                "ERROR: L1TXGOLDSCORE is absent from: "
                + ", ".join(missing_score)
                + "\n  Re-run with --no-score-filter to disable the score gate, or with "
                "--tier <value> to filter by CALLTIER instead."
            )

    # --- sentinels (always dropped from the size distribution) --------------
    kept, sentinel_dropped = split_sentinels(calls, size_field=args.size_field)

    # --- confidence + grouping, in the chosen order -------------------------
    if args.group_first:
        all_groups = keep_unique_across_samples(kept, args.merge_bp, families)
        groups, confidence_dropped = _filter_groups_by_confidence(
            all_groups,
            None if (args.no_score_filter or args.tier) else args.min_score,
            args.tier,
        )
    else:
        min_score = None if (args.no_score_filter or args.tier) else args.min_score
        confident, confidence_dropped = select_confident_calls(kept, min_score=min_score, tier=args.tier)
        groups = keep_unique_across_samples(confident, args.merge_bp, families)

    # --- records whose family is outside the plotted set --------------------
    out_of_scope = [c for c in kept if c.family not in set(families)]
    scope_dropped = [
        Dropped(c, "out-of-scope", f"MEIFAMILY={c.family!r} not in {list(families)}") for c in out_of_scope
    ]

    # --- build site sets, sizes, bins, per-panel plots ----------------------
    settings_lines = _settings_lines(args, families, bin_widths)
    settings_text = "  |  ".join(settings_lines)

    panels: dict[tuple[str, str], tuple[list[int], np.ndarray]] = {}
    bin_count_rows: list[list[object]] = []
    summary_rows: list[list[object]] = []

    sentinel_by_family = Counter(d.call.family for d in sentinel_dropped)
    sentinel_by_sample = Counter(d.call.sample for d in sentinel_dropped)

    for site_set in SITE_SETS:
        site_groups = select_site_set(groups, site_set)
        series = size_series(site_groups, families)
        for family in families:
            sizes = series.get(family, [])
            edges = _bin_edges(sizes, bin_widths[family])
            edges, counts = plot_size_histogram(
                sizes,
                family=family,
                site_set=site_set,
                bin_width=bin_widths[family],
                size_field=args.size_field,
                edges=edges,
                out_path=out_dir / f"size_hist_{family}_{site_set}.png",
                settings_text=settings_text,
            )
            panels[(site_set, family)] = (sizes, edges)
            for idx, count in enumerate(counts):
                bin_count_rows.append(
                    [family, site_set, int(edges[idx]), int(edges[idx + 1]), int(count)]
                )
            stats = _stats(sizes)
            summary_rows.append(
                [
                    family,
                    site_set,
                    stats["n"],
                    int(sentinel_by_family.get(family, 0)),
                    stats["min"],
                    stats["median"],
                    stats["max"],
                ]
            )

    plot_combined(
        panels,
        families,
        SITE_SETS,
        args.size_field,
        out_dir / "size_histograms_combined.png",
        settings_text,
    )

    # --- dropped records ----------------------------------------------------
    all_dropped = sentinel_dropped + confidence_dropped + scope_dropped
    dropped_rows = [
        [
            d.call.source_file,
            d.call.sample,
            d.call.chrom,
            d.call.pos,
            d.call.call_id,
            d.call.family,
            d.call.size,
            d.call.mei_5p,
            d.call.score,
            d.call.tier,
            d.stage,
            d.reason,
        ]
        for d in all_dropped
    ]

    _write_rows(
        out_dir / "dropped_records.tsv",
        [
            "file",
            "sample",
            "chrom",
            "pos",
            "id",
            "family",
            "size",
            "mei_5p",
            "score",
            "tier",
            "stage",
            "reason",
        ],
        dropped_rows,
    )
    _write_rows(
        out_dir / "size_summary.tsv",
        ["family", "site_set", "n", "dropped_sentinel", "min", "median", "max"],
        summary_rows,
    )
    _write_rows(
        out_dir / "bin_counts.tsv",
        ["family", "site_set", "bin_start", "bin_end", "count"],
        bin_count_rows,
    )

    # --- printed summary (active settings + counts) -------------------------
    print("=== retrotransposon-miner insertion-size histograms ===")
    print("Active settings:")
    for line in settings_lines:
        print(f"  {line}")
    print()
    print(f"VCF files        : {len(metas)}")
    print(f"VCF records      : {len(calls)} calls across {len({c.sample for c in calls})} sample(s)")
    print(f"Kept for grouping: {len(kept)} non-sentinel calls")
    print(f"Dropped sentinel : {len(sentinel_dropped)}  (from the size distribution)")
    print(f"Dropped confidence: {len(confidence_dropped)}")
    print(f"Out-of-scope family: {len(scope_dropped)}")
    print()
    print("Sentinel counts per family:")
    for family in sorted(set(sentinel_by_family) | set(families)):
        print(f"  {family:<8} {sentinel_by_family.get(family, 0)}")
    print("Sentinel counts per sample:")
    if sentinel_by_sample:
        for sample in sorted(sentinel_by_sample):
            print(f"  {sample:<12} {sentinel_by_sample[sample]}")
    else:
        print("  (none)")
    print()
    print(f"Sites grouped (union): {len(groups)}  |  private (n_samples==1): "
          f"{len(select_site_set(groups, 'private'))}")
    print()
    header = f"{'family':<8} {'site_set':<8} {'n':>5} {'sent':>5} {'min':>7} {'median':>8} {'max':>7}"
    print("Summary table")
    print(header)
    print("-" * len(header))
    for row in summary_rows:
        family, site_set, n, sent, lo, med, hi = row
        print(f"{family:<8} {site_set:<8} {n:>5} {sent:>5} {_fmt(lo):>7} {_fmt(med):>8} {_fmt(hi):>7}")
    print()
    print(f"Wrote PNGs and TSVs to {out_dir}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
