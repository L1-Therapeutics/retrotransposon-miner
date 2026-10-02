"""Phase 1: is the Alu-into-Alu insertion hotspot at consensus 133 real?

The pre-registered primary test. Levy, Schwartz & Ast 2010 reported a single
prominent hotspot directly after the A-rich linker of the targeted Alu
consensus, at position 133, and this study's earlier HG03086 analysis found a
coarse 120-140 bp peak that is compatible with it. The question here is
whether that concentration survives an opportunity-matched null, on a cohort
roughly an order of magnitude larger and from an independent caller.

Design commitments, fixed before looking at results:

  * The primary window is consensus 133 +/- 5 bp. Slop values 0, 2, 5 and 10
    are reported as a pre-declared sensitivity series, not selected after the
    fact.
  * The primary null conditions only on the cohort's own detectable interval
    per host. k100.Umap is a *short-read* mappability track applied to a
    *long-read* cohort, so it cannot be the primary null; it appears only as
    Null B, and the A-versus-B difference is itself a reported finding about
    how much the choice of mask matters.
  * Null C conditions on local EN-motif density, because Levy's post-linker
    hotspot is notable precisely because the Alu consensus lacks a canonical
    endonuclease site there. Surviving Null C is the interesting outcome.
  * The 133 window is the primary, literature-pre-specified test. The 16-bin
    scan across the host is exploratory and is corrected as a scan with a
    max-statistic permutation FWER, not read off the per-bin p-values.
  * The estimand is a distribution of catalogued sites pooled across 908
    samples, not one genome's biology. Unique sites and carrier counts are
    therefore reported side by side and never collapsed.
  * Uncertainty is host-clustered. Multiple child insertions can share one
    host Alu, so treating events as independent would understate the spread.
    The null resamples hosts, and the observed statistic is additionally
    reported leave-one-host-out.

Replicate count is 10,000 rather than the 1,000 used elsewhere in this
project: 1,000 gives an empirical p-floor near 0.001, which is far too coarse
for a 16-bin scan family where several bins are expected to look extreme.
"""

from __future__ import annotations

import argparse
import collections
import csv
import json
import sys
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np

REPO_SRC = Path(__file__).resolve().parents[1] / "src"
if str(REPO_SRC) not in sys.path:
    sys.path.insert(0, str(REPO_SRC))

DEFAULT_COHORT = (
    Path(__file__).resolve().parents[2]
    / "nested_analysis" / "results_longread" / "per_call_longread_nested.csv"
)
DEFAULT_MAPABILITY = (
    Path.home()
    / "retrotransposon-workdir" / "data" / "public" / "annotation" / "hg38"
    / "mappability" / "k100.Umap.MultiTrackMappability.low_lt0.5.bed"
)

#: Primary window, pre-registered. Slop series reported alongside.
PRIMARY_WINDOW_CENTRE = 133
PRIMARY_SLOP_BP = 5
SLOP_SERIES = (0, 2, 5, 10)

#: Host length band defining the "near-full-length" analysis set.
NEAR_FULL_MIN = 280
NEAR_FULL_MAX = 320

#: Exploratory scan resolution across the host.
SCAN_BIN_BP = 20
SCAN_BINS = tuple(range(0, 320, SCAN_BIN_BP))

DEFAULT_REPLICATES = 10_000


def read_cohort(path: Path) -> list[dict[str, str]]:
    with Path(path).open(newline="") as fh:
        return list(csv.DictReader(fh))


def _to_int(value: str) -> int | None:
    if value is None or value == "":
        return None
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def load_nested_rows(
    rows: Iterable[dict[str, str]],
    coordinate: str = "consensus_offset",
    require_unambiguous: bool = True,
) -> list[dict[str, Any]]:
    """Nested rows on near-full-length hosts with a usable coordinate.

    `coordinate` selects which position each call contributes. The default is
    the projected consensus offset because that is the frame the pre-registered
    133 window is defined in; "host_offset_5p_0based" is available because the
    projection is not exact for every host.
    """
    out: list[dict[str, Any]] = []
    for row in rows:
        if row.get("nested_in_alu_host") != "1":
            continue
        host_len = _to_int(row.get("host_len", ""))
        if host_len is None or not (NEAR_FULL_MIN <= host_len <= NEAR_FULL_MAX):
            continue
        if require_unambiguous and row.get("consensus_mapping") != "unambiguous":
            continue
        offset = _to_int(row.get(coordinate, ""))
        if offset is None:
            continue
        ac = _to_int(row.get("site_allele_count", "")) or 0
        span = _to_int(row.get("consensus_span_bp", ""))
        out.append(
            {
                "chrom": row["chrom"],
                "pos": int(row["pos"]),
                "host_start0": int(row["host_start0"]),
                "host_end0": int(row["host_end0"]),
                "host_strand": row.get("host_strand", ""),
                "host_name": row.get("host_name", ""),
                "host_len": host_len,
                "consensus_span_bp": span,
                "offset": offset,
                "allele_count": ac,
                "insert_strand": row.get("insert_strand", ""),
                "consensus_match_name": row.get("consensus_match_name", ""),
                "offset_drift_bp": _to_int(row.get("offset_drift_bp", "")),
                "perc_resolved": row.get("perc_resolved", ""),
                "not_canonical": row.get("not_canonical", ""),
                "host_key": (row["chrom"], int(row["host_start0"]), int(row["host_end0"])),
            }
        )
    return out


def host_opportunity(host_len: int) -> int:
    """Callable positions for the primary null: the host's own interval.

    The primary null deliberately does not intersect with any external mask.
    Conditioning on a short-read mappability track would be incoherent for a
    long-read cohort, so Umap is a sensitivity analysis only.

    The null is expressed in whichever coordinate frame the test is using.
    When the test is on `host_offset_5p_0based` the host interval is the right
    opportunity. When it is on `consensus_offset` it is not: projecting a host
    position through an alignment does not preserve distance (measured drift is
    6 bp at the median and 39 bp at p90), so resampling uniformly in host space
    and counting in consensus space would compare two different geometries.
    In that case the consensus extent of the host is the opportunity, and
    `consensus_span` must be supplied.
    """
    return max(host_len, 1)


def consensus_opportunity(span: int | None, host_len: int) -> int:
    """Opportunity when counting in consensus coordinates."""
    if span is None or span <= 0:
        return max(host_len, 1)
    return span


def observed_count(offsets: Sequence[int], lo: int, hi: int) -> int:
    """Count offsets in the half-open window [lo, hi)."""
    return sum(1 for o in offsets if lo <= o < hi)


def bin_index(offset: int, width: int = SCAN_BIN_BP, nbins: int = 16) -> int | None:
    if offset < 0:
        return None
    idx = offset // width
    return idx if 0 <= idx < nbins else None


# --------------------------------------------------------------------------
# Nulls
# --------------------------------------------------------------------------


def window_opportunity(offset: int, span: int, lo: int, hi: int) -> int:
    """Positions of the half-open window [lo, hi) that lie inside [0, span).

    The single definition of "how much opportunity does this event's own
    element give the window". Every null here -- A, B and C -- has to use this
    one. They used to re-derive it, and they drifted: Null B skipped the clip
    entirely, which let a window wider than an element contribute more than one
    expected event. That is why the mappability-weighted expectation shipped at
    47.51 against Null A's 22.85 while the median host mappable fraction was
    1.0 -- a gap with nothing to do with mappability.
    """
    if span <= 0 or hi <= lo:
        return 0
    return len(range(max(lo - offset, 0), min(hi - offset, span)))


def expected_uniform(
    call_offsets: Sequence[int],
    opportunities: Sequence[int],
    lo: int,
    hi: int,
) -> float:
    """Opportunity-matched expectation for a window, Null A.

    Each call contributes the fraction of its own element's positions that fall
    in the window, so a call in a short element is not treated as having the
    same opportunity as one in a long element. `opportunities` is expressed in
    the same coordinate frame as `call_offsets`; see `host_opportunity`.
    """
    return sum(
        window_opportunity(offset, span, lo, hi) / span
        for offset, span in zip(call_offsets, opportunities)
        if span > 0
    )


def expected_weighted(
    weights: Sequence[float],
    call_offsets: Sequence[int],
    opportunities: Sequence[int],
    lo: int,
    hi: int,
) -> float:
    """Null A's per-event window opportunity, scaled by a per-event weight.

    This is what Null B and Null C are: the same opportunity as Null A, scaled
    by how callable each element is according to an external track. They used to
    re-derive that opportunity inline and Null B dropped the clip to the
    element's own extent, so a window wider than an element contributed more
    than one expected event. That is how a mappability-weighted expectation
    shipped at 47.51 against Null A's 22.85 while the median host mappable
    fraction was 1.0 -- a gap that had nothing to do with mappability.

    At unit weights it reduces to `expected_uniform` exactly, which is the
    invariant that pins Null B to Null A.
    """
    return sum(
        weight * window_opportunity(offset, span, lo, hi) / span
        for weight, offset, span in zip(weights, call_offsets, opportunities)
        if span > 0
    )


def resample_offsets(
    call_offsets: Sequence[int],
    opportunities: Sequence[int],
    rng: Any,
) -> list[int]:
    """One null replicate: a fresh uniform position inside each call's own element.

    Resampling *within* the element preserves its opportunity exactly, so the
    null asks whether breakpoints concentrate in a particular part of the
    element more than the element's own extent would predict. Resampling is in
    the same frame the observation is counted in.
    """
    return [int(rng.integers(0, max(span, 1))) for span in opportunities]


def host_clustered_bootstrap(
    hosts: Sequence[tuple],
    offsets: Sequence[int],
    lo: int,
    hi: int,
    replicates: int,
    rng: Any,
) -> dict[str, float]:
    """CI for the observed count, resampling whole hosts rather than events.

    26 of 2,559 near-full hosts carry more than one nested event, so event-level
    resampling is very nearly host-level, but doing it properly costs nothing
    and makes the clustering assumption explicit rather than implicit.
    """
    by_host: dict[tuple, list[int]] = collections.defaultdict(list)
    for host, offset in zip(hosts, offsets):
        by_host[host].append(offset)
    keys = list(by_host)
    if not keys:
        return {"lo": 0.0, "hi": 0.0, "n_hosts": 0}
    counts = []
    for _ in range(replicates):
        total = 0
        for _ in range(len(keys)):
            k = keys[int(rng.integers(0, len(keys)))]
            total += observed_count(by_host[k], lo, hi) if hi > lo else 0
        counts.append(total)
    counts.sort()
    return {
        "lo": float(counts[int(0.025 * len(counts))]),
        "hi": float(counts[int(0.975 * len(counts)) - 1]),
        "n_hosts": len(keys),
    }


def leave_one_host_out(
    hosts: Sequence[tuple], offsets: Sequence[int], lo: int, hi: int
) -> dict[str, Any]:
    """Largest single-host influence on the observed count."""
    by_host: dict[tuple, list[int]] = collections.defaultdict(list)
    for host, offset in zip(hosts, offsets):
        by_host[host].append(offset)
    full = observed_count(offsets, lo, hi)
    worst = 0
    for offs in by_host.values():
        n = observed_count(offs, lo, hi)
        if n > worst:
            worst = n
    return {
        "observed": full,
        "n_hosts": len(by_host),
        "max_events_in_one_host": worst,
        "max_share_of_observed": (worst / full) if full else None,
        "observed_minus_worst_host": full - worst,
    }


def _null_window_counts(
    opportunities: Sequence[int],
    lo: int,
    hi: int,
    replicates: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Null replicate counts in [lo, hi), vectorised across replicates.

    Each event independently draws a uniform position inside its own element,
    so the replicate count is a sum of independent Bernoulli terms with
    call-specific probabilities. Sampling it directly, one event at a time, is
    correct but costs minutes at 10,000 replicates; the multinomial shortcut
    is exact in distribution and runs in milliseconds.
    """
    span = np.maximum(np.asarray(opportunities, dtype=np.int64), 1)
    width = np.asarray(max(hi - lo, 0), dtype=np.float64)
    # P(draw in window) for an element of length s and window width w, clipped
    # because a window near the element edge admits fewer positions.
    p = np.clip(width, 0.0, span) / span
    n = span.size
    if n == 0:
        return np.zeros(replicates, dtype=np.int64)
    # Each of the n events is an independent Bernoulli(p_i); the replicate
    # total is their sum, which is Poisson-binomial. Drawing one uniform per
    # event per replicate and counting is exact and avoids relying on a
    # per-element binomial broadcast that numpy will not vectorise this way.
    counts = np.empty(replicates, dtype=np.int64)
    chunk = max(1, 2_000_000 // n)
    for start in range(0, replicates, chunk):
        stop = min(start + chunk, replicates)
        draws = (rng.random((stop - start, n)) * span[None, :]).astype(np.int64)
        counts[start:stop] = ((draws >= lo) & (draws < hi)).sum(axis=1)
    del p
    return counts


def _null_bin_maxima(
    opportunities: Sequence[int],
    replicates: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Per-replicate maximum bin occupancy under the null.

    Each event independently draws a uniform position inside its own element,
    so a replicate is a vector of n uniform draws. Histogramming that vector
    into 16 bins gives the bin counts, and the max over bins is the scan
    statistic. Draws are generated in chunks so a 10,000-replicate scan over
    ~1,300 events stays memory-bounded.
    """
    span = np.maximum(np.asarray(opportunities, dtype=np.int64), 1)
    n = span.size
    nbins = len(SCAN_BINS)
    maxima = np.empty(replicates, dtype=np.int64)
    if n == 0:
        return maxima
    # Keep each chunk under ~2M draws so the histogram stays cheap.
    chunk = max(1, 2_000_000 // n)
    for start in range(0, replicates, chunk):
        stop = min(start + chunk, replicates)
        reps = stop - start
        draws = (rng.random((reps, n)) * span[None, :]).astype(np.int64)
        idx = draws // SCAN_BIN_BP
        np.clip(idx, 0, nbins - 1, out=idx)
        # Bincount per replicate without a Python loop.
        offsets = (np.arange(reps) * nbins)[:, None]
        flat = (idx + offsets).ravel()
        counts = np.bincount(flat, minlength=reps * nbins).reshape(reps, nbins)
        maxima[start:stop] = counts.max(axis=1)
    return maxima


def max_statistic_fwer(
    call_offsets: Sequence[int],
    opportunities: Sequence[int],
    replicates: int,
    rng: np.random.Generator,
) -> dict[str, Any]:
    """Scan-level significance with a max-statistic permutation correction.

    Rather than reading sixteen per-bin p-values and picking the smallest, this
    asks how extreme the most extreme bin is under the null, and compares that
    to the observed maximum. That controls the family-wise error rate across
    the whole scan, which is what an exploratory sixteen-bin look at the data
    actually needs.
    """
    observed_bins = collections.Counter(
        b for b in (bin_index(o) for o in call_offsets) if b is not None
    )
    observed_max = max(observed_bins.values()) if observed_bins else 0
    expected = expected_uniform(call_offsets, opportunities, 0, 320) / len(SCAN_BINS)
    maxima = _null_bin_maxima(opportunities, replicates, rng)
    at_least = int((maxima >= observed_max).sum())
    return {
        "observed_max_bin_count": observed_max,
        "expected_per_bin": expected,
        "null_max_bin_p": (at_least + 1) / (replicates + 1),
        "replicates": replicates,
        "n_bins": len(SCAN_BINS),
    }


# --------------------------------------------------------------------------
# Null B: the short-read mappability mask, as a sensitivity analysis only
# --------------------------------------------------------------------------


def mappable_fractions(
    events: Sequence[dict[str, Any]], bed: Path
) -> list[float]:
    """Per-host mappable fraction for a whole event list, in one bedtools call.

    The k100 track is stored per base and runs to ~700 MB, so parsing it in
    Python costs minutes and holding millions of single-base tuples costs
    gigabytes. One bedtools intersect over all hosts avoids both, and avoids
    spawning a process per host.

    Used only to scale the Null B expectation. It is deliberately a per-host
    scalar rather than a per-position mask, because masking individual bases
    inside a 300 bp repeat would model a short-read artefact, which is exactly
    the incoherence that keeps this out of the primary analysis.
    """
    n = len(events)
    if not Path(bed).exists() or n == 0:
        return [1.0] * n
    import subprocess

    payload = "".join(
        f"{e['chrom']}\t{e['host_start0']}\t{e['host_end0']}\t{i}\n"
        for i, e in enumerate(events)
    )
    proc = subprocess.run(
        [
            "bedtools", "intersect", "-a", "/dev/stdin", "-b", str(bed),
            "-wa", "-wb", "-wo", "-sorted", "-g", "sort -k1,1 -k2,2n",
        ],
        input=payload,
        capture_output=True,
        text=True,
    )
    covered: dict[int, int] = collections.defaultdict(int)
    if proc.returncode == 0:
        for line in proc.stdout.splitlines():
            parts = line.split("\t")
            if len(parts) < 8:
                continue
            try:
                idx = int(parts[3])
                covered[idx] += int(parts[-1])
            except ValueError:
                continue
    out = []
    for i, e in enumerate(events):
        length = max(e["host_end0"] - e["host_start0"], 1)
        out.append(max(0.0, 1.0 - covered.get(i, 0) / length))
    return out


# --------------------------------------------------------------------------
# Stratification
# --------------------------------------------------------------------------


def _bin_table(
    offsets: Sequence[int], expected_per_bin: float
) -> list[dict[str, Any]]:
    counts = collections.Counter(
        b for b in (bin_index(o) for o in offsets) if b is not None
    )
    rows = []
    for b in SCAN_BINS:
        n = counts.get(b // SCAN_BIN_BP, 0)
        rows.append(
            {
                "bin_start_bp": b,
                "bin_end_bp": b + SCAN_BIN_BP,
                "n_sites": n,
                "n_carriers": None,
                "expected_uniform": round(expected_per_bin, 3),
                "enrichment": round(n / expected_per_bin, 3) if expected_per_bin else None,
            }
        )
    return rows


def stratify(
    events: Sequence[dict[str, Any]], lo: int, hi: int
) -> dict[str, Any]:
    """Window counts across the pre-declared stratifications."""
    def subset(pred) -> list[dict[str, Any]]:
        return [e for e in events if pred(e)]

    out: dict[str, Any] = {}
    for label, pred in (
        ("sense", lambda e: e["insert_strand"] == e["host_strand"]),
        # Antisense must be an explicit *opposite* strand, not merely unequal.
        # Using `!=` would sweep an unresolvable insert strand into antisense
        # and bias exactly the contrast being reported.
        (
            "antisense",
            lambda e: e["insert_strand"] in "+-" and e["insert_strand"] != e["host_strand"],
        ),
        ("host_plus_strand", lambda e: e["host_strand"] == "+"),
        ("host_minus_strand", lambda e: e["host_strand"] == "-"),
        ("unresolved_orientation", lambda e: e["insert_strand"] not in "+-"),
    ):
        sub = subset(pred)
        out[label] = {
            "n_sites": len(sub),
            "n_in_window": observed_count([e["offset"] for e in sub], lo, hi),
        }

    subfam: dict[str, list[dict[str, Any]]] = collections.defaultdict(list)
    for e in events:
        subfam[e["consensus_match_name"] or "unknown"].append(e)
    out["by_host_subfamily"] = {
        k: {
            "n_sites": len(v),
            "n_in_window": observed_count([e["offset"] for e in v], lo, hi),
        }
        for k, v in sorted(subfam.items())
    }
    return out


def carrier_weighted_counts(
    events: Sequence[dict[str, Any]], lo: int, hi: int
) -> dict[str, Any]:
    """Site counts and carrier counts, reported side by side, never merged.

    The estimand is a distribution of catalogued sites pooled across samples.
    Two sites with the same profile but 1 and 400 carriers are very different
    pieces of evidence, and reporting only the site count hides that.
    """
    in_window = [e for e in events if lo <= e["offset"] < hi]
    return {
        "n_sites_total": len(events),
        "n_sites_in_window": len(in_window),
        "n_carriers_total": sum(e["allele_count"] for e in events),
        "n_carriers_in_window": sum(e["allele_count"] for e in in_window),
        "median_allele_count_in_window": (
            sorted(e["allele_count"] for e in in_window)[len(in_window) // 2]
            if in_window
            else None
        ),
    }


# --------------------------------------------------------------------------
# Report
# --------------------------------------------------------------------------


def write_report(path: Path, report: dict[str, Any]) -> None:
    primary = report["primary_test"]
    lines = [
        "# Phase 1: position-133 enrichment",
        "",
        "Every number is read back from the cohort table at run time.",
        "",
        "## Pre-registration",
        "",
        f"- Primary window: consensus {report['pre_registration']['centre']} "
        f"+/- {report['pre_registration']['slop_bp']} bp.",
        f"- Slop sensitivity series: {report['pre_registration']['slop_series']} bp.",
        f"- Replicates: {report['pre_registration']['replicates']:,}.",
        "- Primary null: each call's own host interval (Null A).",
        "- Null B: additionally scaled by a short-read mappability track; "
        "reported as a sensitivity analysis only, because Umap is not SVAN's "
        "callability.",
        "- Null C: additionally conditioned on local EN-motif density; "
        "reported only when a motif table is supplied.",
        "",
        "## Analysis set",
        "",
        f"- Nested sites on {NEAR_FULL_MIN}-{NEAR_FULL_MAX} bp hosts with an "
        f"unambiguous consensus mapping: **{primary['n_sites']}**.",
        f"- Distinct host Alus: **{primary['n_hosts']}**.",
        f"- Carrier-weighted total: {primary['carriers']['n_carriers_total']}.",
        f"- Null opportunity frame: **{primary['opportunity_frame']}** "
        f"({primary['opportunity_definition']}; median "
        f"{primary['median_opportunity']} bp).",
        "",
        "## Primary test (Null A)",
        "",
        "| quantity | value |",
        "|---|---|",
        f"| observed sites in window | {primary['observed']} |",
        f"| opportunity-matched expectation | {primary['expected_a']:.2f} |",
        f"| enrichment | {primary['enrichment_a']:.2f}x |",
        f"| empirical p (one-sided, greater) | {primary['p_a']:.5f} |",
        f"| host-clustered 95% CI on observed | "
        f"{primary['bootstrap']['lo']:.0f} - {primary['bootstrap']['hi']:.0f} |",
        "",
        "## Slop sensitivity",
        "",
        "| slop | observed | expected | enrichment | p |",
        "|---|---|---|---|---|",
    ]
    for s in report["slop_sensitivity"]:
        lines.append(
            f"| +/-{s['slop_bp']} bp | {s['observed']} | {s['expected_a']:.2f} | "
            f"{s['enrichment_a']:.2f}x | {s['p_a']:.5f} |"
        )
    lines += [
        "",
        "## Null B (short-read mappability, sensitivity)",
        "",
    ]
    if report["null_b"]["ran"]:
        b = report["null_b"]
        lines += [
            f"- Expected under Null B: {b['expected_b']:.2f} "
            f"(vs {primary['expected_a']:.2f} under Null A).",
            f"- Enrichment under Null B: {b['enrichment_b']:.2f}x.",
            f"- Median per-host mappable fraction: {b['median_mappable_fraction']:.3f}.",
            "",
            "The A-versus-B difference is a finding about mask choice, not a "
            "correction: a short-read track cannot describe a long-read "
            "cohort's callability, so the primary result is Null A.",
        ]
    else:
        lines.append("Not run: " + str(report["null_b"]["reason"]))
    lines += [
        "",
        "## Null C (EN-motif conditioned, sensitivity)",
        "",
    ]
    if report["null_c"]["ran"]:
        c = report["null_c"]
        lines += [
            f"- Expected under Null C: {c['expected_c']:.2f}.",
            f"- Enrichment under Null C: {c['enrichment_c']:.2f}x.",
            f"- {c['note']}",
        ]
    else:
        lines.append("Not run: " + str(report["null_c"]["reason"]))
    lines += [
        "",
        "## Exploratory 16-bin scan (corrected as a scan)",
        "",
        "| bin | sites | carriers | expected | enrichment |",
        "|---|---|---|---|---|",
    ]
    for r in report["scan_profile"]:
        lines.append(
            f"| {r['bin_start_bp']}-{r['bin_end_bp']} | {r['n_sites']} | "
            f"{r['n_carriers']} | {r['expected_uniform']:.1f} | "
            f"{r['enrichment']:.2f}x |"
        )
    fw = report["scan_fwer"]
    lines += [
        "",
        f"Max-statistic permutation FWER across all {fw['n_bins']} bins: "
        f"**p = {fw['null_max_bin_p']:.5f}** "
        f"({fw['replicates']:,} replicates; observed max bin = "
        f"{fw['observed_max_bin_count']}, expected "
        f"{fw['expected_per_bin']:.1f}/bin).",
        "",
        "This is the scan-level p-value. The per-bin counts above are "
        "descriptive and are not individually corrected.",
        "",
        "## Host clustering",
        "",
        f"- Leave-one-host-out: {report['influence']['n_hosts']} hosts, largest "
        f"single-host contribution {report['influence']['max_events_in_one_host']} "
        f"site(s) "
        f"({report['influence']['max_share_of_observed'] * 100:.1f}% of the "
        f"window total).",
        f"- Window count with that host removed: "
        f"{report['influence']['observed_minus_worst_host']}.",
        "",
        "## Stratifications (window counts)",
        "",
        "| stratum | n_sites | n_in_window |",
        "|---|---|---|",
    ]
    for key, val in report["stratification"].items():
        if key == "by_host_subfamily":
            continue
        lines.append(f"| {key} | {val['n_sites']} | {val['n_in_window']} |")
    lines += [
        "",
        "## Estimand and disclosures",
        "",
        "- The cohort is a catalogue of insertion **sites** pooled across 908 "
        "samples, not one genome's insertions. Site and carrier counts are "
        "reported separately and never merged.",
        "- Positions are consensus-relative, projected through a global "
        "alignment of each host to its own subfamily consensus. Projection is "
        "exact only for the unambiguous subset used here.",
        "- This is an enrichment test. It is not a detection-accuracy result "
        "and licenses no precision/recall claim.",
        "- Catalogue polymorphisms have passed selection; these are survivor "
        "distributions, not raw de novo targeting rates.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cohort", type=Path, default=DEFAULT_COHORT)
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--mapability-bed", type=Path, default=DEFAULT_MAPABILITY)
    p.add_argument("--replicates", type=int, default=DEFAULT_REPLICATES)
    p.add_argument("--seed", type=int, default=20261002)
    p.add_argument(
        "--coordinate",
        choices=["consensus_offset", "host_offset_5p_0based"],
        default="consensus_offset",
    )
    p.add_argument(
        "--include-ambiguous",
        action="store_true",
        help="Include ambiguous-indel consensus mappings (sensitivity only).",
    )
    p.add_argument(
        "--en-motif-bed",
        type=Path,
        default=None,
        help="EN-motif intervals; enables the Null C sensitivity analysis.",
    )
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    args.outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    rows = read_cohort(args.cohort)
    events = load_nested_rows(
        rows,
        coordinate=args.coordinate,
        require_unambiguous=not args.include_ambiguous,
    )
    if not events:
        raise SystemExit("no nested rows on near-full-length hosts")
    offsets = [e["offset"] for e in events]
    hosts = [e["host_key"] for e in events]
    # The null must be resampled in the same frame the observation is counted
    # in. In consensus space the opportunity is the aligned consensus extent,
    # not the host length, because the projection does not preserve distance.
    if args.coordinate == "consensus_offset":
        opportunities = [
            consensus_opportunity(e["consensus_span_bp"], e["host_len"])
            for e in events
        ]
    else:
        opportunities = [host_opportunity(e["host_len"]) for e in events]
    n = len(events)
    print(
        f"{n} nested sites on near-full hosts, "
        f"{len(set(hosts))} distinct host Alus",
        flush=True,
    )

    lo = PRIMARY_WINDOW_CENTRE - PRIMARY_SLOP_BP
    hi = PRIMARY_WINDOW_CENTRE + PRIMARY_SLOP_BP + 1  # half-open
    observed = observed_count(offsets, lo, hi)
    expected_a = expected_uniform(offsets, opportunities, lo, hi)
    print(
        f"window [{lo},{hi}): observed {observed}, expected {expected_a:.2f}",
        flush=True,
    )

    at_least = int((_null_window_counts(opportunities, lo, hi, args.replicates, rng) >= observed).sum())
    p_a = (at_least + 1) / (args.replicates + 1)

    slop_rows = []
    for slop in SLOP_SERIES:
        slo, shi = PRIMARY_WINDOW_CENTRE - slop, PRIMARY_WINDOW_CENTRE + slop + 1
        obs = observed_count(offsets, slo, shi)
        exp = expected_uniform(offsets, opportunities, slo, shi)
        null_counts = _null_window_counts(opportunities, slo, shi, args.replicates, rng)
        hits = int((null_counts >= obs).sum())
        slop_rows.append(
            {
                "slop_bp": slop,
                "observed": obs,
                "expected_a": exp,
                "enrichment_a": (obs / exp) if exp else None,
                "p_a": (hits + 1) / (args.replicates + 1),
            }
        )

    # Null B: scale each call's opportunity by its host's mappable fraction.
    null_b: dict[str, Any] = {"ran": False, "reason": "mapability bed not found"}
    if args.mapability_bed and Path(args.mapability_bed).exists():
        fracs = mappable_fractions(events, args.mapability_bed)
        # Same per-event window opportunity as Null A, scaled by mappability.
        # This used to be `(hi - lo) * f / span`, which ignored the clip to the
        # element's own extent and so was not Null A under a different weight.
        expected_b = expected_weighted(fracs, offsets, opportunities, lo, hi)
        null_b = {
            "ran": True,
            "expected_b": expected_b,
            "enrichment_b": (observed / expected_b) if expected_b else None,
            "median_mappable_fraction": sorted(fracs)[len(fracs) // 2],
            "note": (
                "k100.Umap is a short-read track. SVAN's callability is not "
                "Umap's, so this scales the expectation rather than masking "
                "individual bases inside a 300 bp repeat."
            ),
        }
        print(
            f"Null B expected {expected_b:.2f} "
            f"(vs A {expected_a:.2f}); median host mappable fraction "
            f"{null_b['median_mappable_fraction']:.3f}",
            flush=True,
        )

    # Null C: requires an EN-motif interval table.
    null_c: dict[str, Any] = {
        "ran": False,
        "reason": "no --en-motif-bed supplied; not fabricated",
    }
    if args.en_motif_bed and Path(args.en_motif_bed).exists():
        # A motif-conditioned expectation uses the EN-motif coverage of each
        # host's element in place of the mappability fraction.
        weights = [f or 1e-6 for f in mappable_fractions(events, args.en_motif_bed)]
        expected_c = expected_weighted(weights, offsets, opportunities, lo, hi)
        null_c = {
            "ran": True,
            "expected_c": expected_c,
            "enrichment_c": (observed / expected_c) if expected_c else None,
            "note": (
                "Expectation only. No motif-conditioned Monte Carlo null is "
                "implemented: the position null resamples breakpoints uniformly "
                "within each element, which says nothing about which elements a "
                "motif table would have made callable. The p-value reported "
                "under Null A is not restated here as if it were motif-"
                "conditioned."
            ),
        }

    scan = _bin_table(offsets, n / len(SCAN_BINS) if n else 0.0)
    bin_carriers: collections.Counter = collections.Counter()
    for e in events:
        b = bin_index(e["offset"])
        if b is not None:
            bin_carriers[SCAN_BINS[b]] += e["allele_count"]
    for row in scan:
        row["n_carriers"] = bin_carriers.get(row["bin_start_bp"], 0)

    fwer = max_statistic_fwer(offsets, opportunities, args.replicates, rng)

    report = {
        "phase": 1,
        "description": "Position-133 enrichment test, pre-registered window.",
        "inputs": {
            "cohort": str(args.cohort),
            "mapability_bed": str(args.mapability_bed) if args.mapability_bed else None,
            "en_motif_bed": str(args.en_motif_bed) if args.en_motif_bed else None,
        },
        "pre_registration": {
            "centre": PRIMARY_WINDOW_CENTRE,
            "slop_bp": PRIMARY_SLOP_BP,
            "slop_series": list(SLOP_SERIES),
            "replicates": args.replicates,
            "seed": args.seed,
            "coordinate": args.coordinate,
            "include_ambiguous": args.include_ambiguous,
            "primary_null": "within-host uniform (Null A)",
        },
        "primary_test": {
            "window_lo_inclusive": lo,
            "window_hi_exclusive": hi,
            "n_sites": n,
            "n_hosts": len(set(hosts)),
            "observed": observed,
            "expected_a": expected_a,
            "enrichment_a": (observed / expected_a) if expected_a else None,
            "p_a": p_a,
            "opportunity_frame": args.coordinate,
            "opportunity_definition": (
                "aligned consensus extent per host"
                if args.coordinate == "consensus_offset"
                else "host interval length"
            ),
            "median_opportunity": sorted(opportunities)[len(opportunities) // 2],
            "bootstrap": host_clustered_bootstrap(
                hosts, offsets, lo, hi, 2000, rng
            ),
            "carriers": carrier_weighted_counts(events, lo, hi),
        },
        "slop_sensitivity": slop_rows,
        "null_b": null_b,
        "null_c": null_c,
        "scan_profile": scan,
        "scan_fwer": fwer,
        "influence": leave_one_host_out(hosts, offsets, lo, hi),
        "stratification": stratify(events, lo, hi),
    }
    (args.outdir / "position_enrichment.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    write_report(args.outdir / "position_enrichment.md", report)

    print(
        f"\nenrichment {report['primary_test']['enrichment_a']:.2f}x  "
        f"p={p_a:.5f}"
    )
    print(
        f"scan FWER p={fwer['null_max_bin_p']:.5f} "
        f"(max bin {fwer['observed_max_bin_count']} vs "
        f"{fwer['expected_per_bin']:.1f} expected)"
    )
    print(f"wrote {args.outdir / 'position_enrichment.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
