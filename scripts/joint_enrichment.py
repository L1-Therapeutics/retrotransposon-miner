"""Joint nested-signature enrichment on the multi-sample dedup output.

Phase 4 of the nested-insertion study, extending the frozen plan
(`docs/NESTED_ALU_ANALYSIS_PLAN.md`), which remains authoritative. This script
tests *pre-specified joint cells* only. The full family x position x orientation
table is roughly 1.2 million cells for this cohort, and scanning a table that
size to find whichever cells look good is how a pre-registration stops being one.
Every cell tested here is named in `CELLS` below before any data is read.

The cells
---------
Host family Alu, position measured as `host_offset_5p_0based` (the strand-aware
host-relative breakpoint offset defined in the Phase 0 docstring):

  PRIMARY    {Alu, 120-140, sense}          the joint cell around consensus 133
  SECONDARY  {Alu, 280-300, sense}          the tail region of Phase 2a
  SECONDARY  position x orientation interaction at 120-140
  SECONDARY  position x orientation interaction at 280-300

The three secondary cells are Holm-adjusted within the Alu family. The primary
cell is not adjusted against them: it is the one pre-specified cell, and
adjusting it against cells that were named after it is not a correction.

The null
--------
Host-stratified conditional randomization, per the pre-registration. No
log-linear model is fitted.

  - The analysis unit is the host element -- one physical copy of a repeat
    element, keyed by contig, interval and strand.
  - Within each host the observed number of nested events is held fixed, and
    their breakpoints are redrawn uniformly from that host's opportunity. This
    makes the null conditional on the host, so a host that simply hosts more
    insertions cannot masquerade as positional preference.
  - For the two interaction cells, sense/antisense labels are additionally
    permuted within each host while preserving that host's orientation totals, so
    the interaction is tested against a null in which position and orientation
    carry no information about each other.
  - A position bin with zero opportunity inside a host -- which happens whenever
    the host is shorter than the bin's start -- is EXCLUDED from the null for
    that host. It is never pseudocounted into a probability.

Opportunity
-----------
The frozen plan ranks the host-interval opportunity null as primary and a
mappability-scaled null as sensitivity only, on the grounds that "a short-read
track is not this cohort's callability". No mappability or gap track exists in
this workspace. The RepeatMasker-derived substitute was measured rather than
assumed: a strict repeat mask is degenerate here, because the host element is
itself a RepeatMasker annotation, so masking repeats would delete the entire
opportunity; and once the host's own annotation is excluded, the residual
overlap by *other* repeats is 0.000 of host length for every Alu host measured.
The mappability-scaled variant is therefore identical to the host-interval null
on this data, which is reported as a measured deviation rather than quietly
assumed either way.

What is deliberately not reported
---------------------------------
No multiplied "1.47x * 3.0x = 4.4x" figure is produced anywhere. The two
per-factor estimates appear only as the benchmark the joint estimate is compared
against, and that comparison -- stated as a verdict, not as a number to quote --
is the multiplication test.

Private sites
-------------
Per the plan's estimand statement, private (single-carrier) sites are a separate
exploratory analysis and no primary-test p-value is computed over them. Every
cell is therefore reported twice: once on shared sites, which is the primary
analysis, and once on private sites, which is exploratory and carries no p-value.
Site counts and carrier counts are reported side by side and never summed.
"""

from __future__ import annotations

import argparse
import collections
import csv
import json
import math
import sys
from pathlib import Path
from typing import Any, Sequence

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import nested_multi_sample_common as common  # noqa: E402

DEFAULT_REPLICATES = 10_000
DEFAULT_BOOTSTRAP = 2_000
DEFAULT_SEED = 20261001
DEFAULT_RMSK = common.ANALYSIS_DIR / "data" / "rmsk.txt.gz"

SENSE = common.SENSE
ANTISENSE = common.ANTISENSE


class Cell:
    """One pre-specified cell. Cells are declared, never discovered."""

    def __init__(
        self,
        cell_id: str,
        role: str,
        bin_lo: int,
        bin_hi: int,
        interaction: bool = False,
    ) -> None:
        self.cell_id = cell_id
        self.role = role
        self.bin_lo = bin_lo
        self.bin_hi = bin_hi
        self.interaction = interaction

    @property
    def label(self) -> str:
        if self.interaction:
            return (
                f"Alu x position {self.bin_lo}-{self.bin_hi} x orientation "
                "(sense minus antisense)"
            )
        return f"Alu x position {self.bin_lo}-{self.bin_hi} x sense"

    def as_dict(self) -> dict[str, Any]:
        return {
            "cell_id": self.cell_id,
            "role": self.role,
            "label": self.label,
            "host_family": "Alu",
            "position_bin": [self.bin_lo, self.bin_hi],
            "position_measure": "host_offset_5p_0based",
            "statistic": (
                "n_sense_in_bin - n_antisense_in_bin"
                if self.interaction
                else "n_sense_in_bin"
            ),
            "orientation_permuted_in_null": self.interaction,
        }


#: The complete, pre-specified cell list. Order is the reporting order.
CELLS: tuple[Cell, ...] = (
    Cell("P1_primary", "primary", 120, 140),
    Cell("S1_tail_sense", "secondary", 280, 300),
    Cell("S2_interaction_120_140", "secondary", 120, 140, interaction=True),
    Cell("S3_interaction_280_300", "secondary", 280, 300, interaction=True),
)

#: Structural-zero rule: a host shorter than this contributes no opportunity to a
#: bin and is excluded from that bin's null rather than pseudocounted.
def opportunity_fraction(host_len: int, bin_lo: int, bin_hi: int) -> float:
    """Fraction of a host of `host_len` bp that the bin occupies.

    Zero when the bin lies entirely beyond the host's own end, which is the
    structural zero the pre-registration says to exclude rather than
    pseudocount. The bin is measured in host-relative offset space and the host
    spans offsets [0, host_len).

    Both bin edges are INCLUSIVE, because that is how the null draws and tests
    offsets in `simulate`. So a host of 300 bp against the [120, 140] bin has
    21 qualifying offsets -- 120 through 140 -- and the fraction is 21/300, not
    20/300. Clamping the upper edge to `host_len` and then taking `hi - lo`
    treats the bin edge as exclusive, which understates every host longer than
    `bin_hi` by one offset and puts the closed form below the simulated mean.
    `test_opportunity_fraction_counts_the_same_offsets_the_null_draws` pins the
    two together.
    """
    if host_len <= bin_lo:
        return 0.0
    lo = max(0, bin_lo)
    # Largest offset that both the bin and the host actually contain.
    hi = min(host_len - 1, bin_hi)
    if hi < lo:
        return 0.0
    return (hi - lo + 1) / host_len


# --------------------------------------------------------------------------
# host stratification
# --------------------------------------------------------------------------


def build_host_table(
    sites: Sequence[dict[str, Any]], host_family: str = "Alu"
) -> dict[tuple, dict[str, Any]]:
    """Group events by host copy, restricted to one host family.

    Only nested analysis events are used. The caller has already enforced the
    plan's private/shared split before this point.
    """
    table: dict[tuple, dict[str, Any]] = {}
    for site in sites:
        if common.host_family(site) != host_family:
            continue
        key = common.host_key(site)
        entry = table.get(key)
        if entry is None:
            entry = {
                "key": key,
                "label": common.host_label(key),
                "host_len": int(site["host_len"]),
                "events": [],
                "sites": [],
            }
            table[key] = entry
        entry["events"].append((int(site["host_offset"]), site["orientation"], site))
        entry["sites"].append(site)
    return table


# --------------------------------------------------------------------------
# statistics
# --------------------------------------------------------------------------


def observed_statistic(
    host: dict[str, Any], cell: Cell
) -> int:
    lo, hi = cell.bin_lo, cell.bin_hi
    in_bin_sense = 0
    in_bin_anti = 0
    for offset, orientation, _site in host["events"]:
        if lo <= offset <= hi:
            if orientation == SENSE:
                in_bin_sense += 1
            elif orientation == ANTISENSE:
                in_bin_anti += 1
    if cell.interaction:
        return in_bin_sense - in_bin_anti
    return in_bin_sense


def analytic_expectation(hosts: Sequence[dict[str, Any]], cell: Cell) -> float:
    """Closed-form expectation under the same conditional null.

    Reported alongside the simulated mean so the two can be checked against each
    other. For a non-interaction cell the expectation is exactly
    sum_h n_sense_h * frac_h. For an interaction cell both orientations are
    redrawn uniformly, so the expectation is sum_h frac_h * (n_sense_h -
    n_anti_h), which is zero whenever a host is orientation-balanced.
    """
    total = 0.0
    for host in hosts:
        fraction = opportunity_fraction(host["host_len"], cell.bin_lo, cell.bin_hi)
        if fraction <= 0.0:
            continue
        n_sense = sum(1 for _o, orient, _s in host["events"] if orient == SENSE)
        n_anti = sum(1 for _o, orient, _s in host["events"] if orient == ANTISENSE)
        if cell.interaction:
            total += (n_sense - n_anti) * fraction
        else:
            total += n_sense * fraction
    return total


def simulate(
    hosts: Sequence[dict[str, Any]],
    cell: Cell,
    replicates: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Draw the null distribution of the statistic.

    Vectorised across hosts, which matters: a per-host Python loop would be
    O(replicates x hosts) and would dominate the runtime at 10,000 replicates.

    Each host keeps its observed event count. Breakpoints are redrawn uniformly
    from that host's own interval. For interaction cells the sense/antisense
    labels are additionally permuted *within* each host, preserving that host's
    orientation totals, via a per-host random ordering of the event indices.

    Hosts whose length leaves no opportunity in the bin need no special case:
    every offset they can draw lies below `bin_lo`, so they can never land in
    it. The structural-zero exclusion is therefore enforced by the arithmetic
    itself, and is separately counted for reporting. Note that this holds only
    for hosts shorter than the bin's lower edge -- a host *longer* than the bin
    does have offsets that land in it, which is correct and is why
    `opportunity_fraction` decides eligibility rather than host length alone.
    """
    lengths = np.array([h["host_len"] for h in hosts], dtype=np.int64)
    counts = np.array([len(h["events"]) for h in hosts], dtype=np.int64)
    n_sense_per_host = np.array(
        [sum(1 for _o, o, _s in h["events"] if o == SENSE) for h in hosts],
        dtype=np.int64,
    )
    total_events = int(counts.sum())
    out = np.zeros(replicates, dtype=np.int64)
    if total_events == 0:
        return out

    host_of_event = np.repeat(np.arange(len(hosts), dtype=np.int64), counts)
    length_of_event = np.repeat(lengths, counts)
    # Fixed sense membership for the non-interaction null: positions move, the
    # labels do not.
    is_sense_event = np.zeros(total_events, dtype=bool)
    start = 0
    for i, host in enumerate(hosts):
        n = int(counts[i])
        k = int(n_sense_per_host[i])
        if k:
            # label the first k sense-ordered events of this host as sense
            order = sorted(
                range(n), key=lambda j: 0 if host["events"][j][1] == SENSE else 1
            )
            for j in order[:k]:
                is_sense_event[start + j] = True
        start += n

    lo, hi = cell.bin_lo, cell.bin_hi
    for r in range(replicates):
        offsets = (rng.random(total_events) * length_of_event).astype(np.int64)
        in_bin = (offsets >= lo) & (offsets <= hi)
        if not cell.interaction:
            out[r] = int(np.count_nonzero(in_bin & is_sense_event))
        else:
            # within-host permutation of the orientation labels
            order = np.lexsort((rng.random(total_events), host_of_event))
            sense = np.zeros(total_events, dtype=bool)
            start = 0
            for i in range(len(hosts)):
                n = int(counts[i])
                k = int(n_sense_per_host[i])
                if k:
                    sense[order[start : start + k]] = True
                start += n
            out[r] = int(np.count_nonzero(in_bin & sense)) - int(
                np.count_nonzero(in_bin & ~sense)
            )
    return out


def monte_carlo_p(observed: int, null: np.ndarray) -> float:
    """P(simulated >= observed), with the +1 correction required by the plan."""
    b = null.size
    if b == 0:
        return float("nan")
    return (1.0 + float(np.count_nonzero(null >= observed))) / (b + 1)


def host_clustered_bootstrap(
    hosts: Sequence[dict[str, Any]],
    cell: Cell,
    draws: int,
    rng: np.random.Generator,
) -> tuple[float | None, float | None]:
    """Percentile CI for observed/expected, resampling whole hosts.

    Resampling hosts rather than events keeps the CI honest about the fact that
    events in one host are not independent. On a resample that contains no host
    contributing to the cell the ratio is undefined; those draws are dropped and
    the count of dropped draws is reported by the caller.
    """
    if not hosts:
        return None, None
    ratios: list[float] = []
    n = len(hosts)
    for _ in range(draws):
        pick = rng.integers(0, n, size=n)
        sample = [hosts[i] for i in pick]
        obs = sum(observed_statistic(h, cell) for h in sample)
        exp = analytic_expectation(sample, cell)
        if exp > 0:
            ratios.append(obs / exp)
    if not ratios:
        return None, None
    arr = np.sort(np.array(ratios))
    lo = float(np.quantile(arr, 0.025))
    hi = float(np.quantile(arr, 0.975))
    return lo, hi


def holm(pvalues: Sequence[float | None]) -> list[float | None]:
    """Holm step-down. None (unestimable) entries stay None and are not counted."""
    idx = [i for i, p in enumerate(pvalues) if p is not None]
    out: list[float | None] = [None] * len(pvalues)
    if not idx:
        return out
    order = sorted(idx, key=lambda i: pvalues[i])  # type: ignore[index]
    m = len(order)
    running = 0.0
    for rank, i in enumerate(order):
        value = min(1.0, (m - rank) * pvalues[i])  # type: ignore[operator]
        running = max(running, value)
        out[i] = running
    return out


# --------------------------------------------------------------------------
# per-cell execution
# --------------------------------------------------------------------------


def run_cell(
    cell: Cell,
    shared_sites: Sequence[dict[str, Any]],
    private_sites: Sequence[dict[str, Any]],
    replicates: int,
    bootstrap: int,
    seed: int,
) -> dict[str, Any]:
    hosts_all = build_host_table(shared_sites)
    eligible = [
        h for h in hosts_all.values() if opportunity_fraction(h["host_len"], cell.bin_lo, cell.bin_hi) > 0
    ]
    structural_zero = len(hosts_all) - len(eligible)
    observed = sum(observed_statistic(h, cell) for h in eligible)
    analytic = analytic_expectation(eligible, cell)

    null = np.zeros(replicates, dtype=np.int64)
    p_value: float | None = None
    expected_sim: float | None = None
    rng = np.random.default_rng(seed)
    if eligible and replicates > 0:
        null = simulate(eligible, cell, replicates, rng)
        p_value = monte_carlo_p(observed, null)
        expected_sim = float(null.mean())

    expected = expected_sim if expected_sim is not None else analytic
    effect = (observed / expected) if expected and expected > 0 else None
    # A ratio is meaningless when the denominator approaches zero, and the
    # bootstrap interval shows it: an expected of -0.4 yields a lower bound of
    # -10 and an upper bound of +357. Where that happens the difference is the
    # interpretable quantity and the ratio is flagged rather than printed as if
    # it meant something.
    effect_unstable = bool(
        expected is not None and abs(expected) < 1.0
    )
    ci_lo, ci_hi = host_clustered_bootstrap(eligible, cell, bootstrap, rng)

    n_events = sum(len(h["events"]) for h in eligible)
    n_sites = sum(len(h["sites"]) for h in eligible)
    n_carriers = sum(
        int(s["n_carriers"]) for h in eligible for s in h["sites"]
    )
    # private exploratory pass: no p-value is computed for it, by plan
    priv_hosts = build_host_table(private_sites)
    priv_eligible = [
        h for h in priv_hosts.values() if opportunity_fraction(h["host_len"], cell.bin_lo, cell.bin_hi) > 0
    ]
    priv_observed = sum(observed_statistic(h, cell) for h in priv_eligible)
    priv_expected = analytic_expectation(priv_eligible, cell)

    return {
        **cell.as_dict(),
        "hosts_total": len(hosts_all),
        "hosts_with_opportunity": len(eligible),
        "hosts_structural_zero_excluded": structural_zero,
        "sites_in_scope": n_sites,
        "carrier_observations_in_scope": n_carriers,
        "events_in_scope": n_events,
        "observed": observed,
        "expected_analytic": analytic,
        "expected_null_mean": expected_sim,
        "observed_minus_expected": (
            observed - expected if expected is not None else None
        ),
        "effect_size_observed_over_expected": effect,
        "effect_size_unstable_near_zero_expected": effect_unstable,
        "bootstrap_ci_low": ci_lo,
        "bootstrap_ci_high": ci_hi,
        "replicates": replicates,
        "monte_carlo_p": p_value,
        "private_exploratory_observed": priv_observed,
        "private_exploratory_expected": priv_expected,
        "private_exploratory_sites": sum(len(h["sites"]) for h in priv_eligible),
    }


def per_factor_estimates(
    shared_sites: Sequence[dict[str, Any]], cell: Cell, replicates: int, seed: int
) -> dict[str, Any]:
    """The two marginals of the joint cell, used only as a benchmark.

    `position_factor` is the enrichment of the bin over both orientations.
    `orientation_factor` is the sense fraction over the whole host family,
    measured against parity, since exchangeability within host makes any
    host-relative orientation expectation exactly 1 by construction.

    The product of the two is reported as a benchmark for the verdict and is
    never presented as an effect size of the joint cell.
    """
    hosts = [
        h
        for h in build_host_table(shared_sites).values()
        if opportunity_fraction(h["host_len"], cell.bin_lo, cell.bin_hi) > 0
    ]
    if not hosts:
        return {"status": "not_estimable"}

    n_in_bin = 0
    n_sense = 0
    n_events = 0
    expected_in_bin = 0.0
    for host in hosts:
        fraction = opportunity_fraction(host["host_len"], cell.bin_lo, cell.bin_hi)
        expected_in_bin += len(host["events"]) * fraction
        for offset, orientation, _s in host["events"]:
            n_events += 1
            if orientation == SENSE:
                n_sense += 1
            if cell.bin_lo <= offset <= cell.bin_hi:
                n_in_bin += 1
    position_factor = (n_in_bin / expected_in_bin) if expected_in_bin > 0 else None
    orientation_factor = (n_sense / n_events) / 0.5 if n_events else None
    return {
        "status": "estimated",
        "position_factor_bin_enrichment": position_factor,
        "orientation_factor_sense_vs_parity": orientation_factor,
        "n_events": n_events,
        "n_in_bin": n_in_bin,
        "n_sense": n_sense,
        "benchmark_product": (
            position_factor * orientation_factor
            if position_factor is not None and orientation_factor is not None
            else None
        ),
    }


# --------------------------------------------------------------------------
# reporting
# --------------------------------------------------------------------------


def opportunity_diagnostic(
    sites: Sequence[dict[str, Any]], rmsk_path: Path, host_family: str = "Alu"
) -> dict[str, Any]:
    """Measure, rather than assume, what a repeat-mask opportunity would be.

    Two facts are established: a strict repeat mask is degenerate (the host is
    itself an annotated repeat, so masking deletes the whole opportunity), and
    once the host's own annotation is excluded, how much of the host interval
    other repeats intrude on.
    """
    hosts = build_host_table(sites, host_family)
    wanted: dict[str, list[tuple[int, int, str]]] = collections.defaultdict(list)
    for host in hosts.values():
        key = host["key"]
        wanted[key[0]].append((int(key[2]), int(key[3]), key[1]))
    if not Path(rmsk_path).exists():
        return {
            "status": "rmsk_unavailable",
            "rmsk_path": str(rmsk_path),
            "consequence": (
                "opportunity could not be measured; the host-interval null is used "
                "and the mappability-scaled sensitivity is reported as unavailable"
            ),
        }
    annotations = common.read_rmsk_intervals(Path(rmsk_path), dict(wanted))
    fractions: list[float] = []
    for host in hosts.values():
        key = host["key"]
        nearby = annotations.get(key[0], [])
        fractions.append(
            common.residual_mask_fraction(
                (int(key[2]), int(key[3]), key[1]),
                [a for a in nearby if a[0] < int(key[3]) and a[1] > int(key[2])],
            )
        )
    if not fractions:
        return {"status": "no_hosts_measured"}
    arr = np.array(fractions)
    return {
        "status": "measured",
        "rmsk_path": str(rmsk_path),
        "hosts_measured": len(fractions),
        "residual_mask_fraction_median": float(np.median(arr)),
        "residual_mask_fraction_max": float(arr.max()),
        "hosts_with_zero_opportunity_from_other_repeats": int(np.count_nonzero(arr <= 0.001)),
        "hosts_fully_masked_by_other_repeats": int(np.count_nonzero(arr >= 0.999)),
        "verdict": (
            "repeat_mask_opportunity_identical_to_host_interval"
            if float(arr.max()) <= 0.001
            else "repeat_mask_opportunity_differs_from_host_interval"
        ),
        "consequence": (
            "excluding the host's own RepeatMasker annotation, no other repeat "
            "intrudes on any measured host interior, so the mappability-scaled "
            "variant collapses onto the host-interval null on this data; a strict "
            "repeat mask would be degenerate because the host is itself an "
            "annotated repeat. Recorded as a deviation, not assumed."
        ),
    }


def write_cells_csv(path: Path, cells: Sequence[dict[str, Any]]) -> None:
    fields = [
        "cell_id", "role", "label", "host_family", "bin_lo", "bin_hi",
        "statistic", "hosts_total", "hosts_with_opportunity",
        "hosts_structural_zero_excluded", "sites_in_scope",
        "carrier_observations_in_scope", "observed", "expected_analytic",
        "expected_null_mean", "effect_size_observed_over_expected",
        "bootstrap_ci_low", "bootstrap_ci_high", "replicates", "monte_carlo_p",
        "holm_adjusted_p", "private_exploratory_observed",
        "private_exploratory_expected", "private_exploratory_sites",
    ]
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for cell in cells:
            row = dict(cell)
            row["bin_lo"] = cell["position_bin"][0]
            row["bin_hi"] = cell["position_bin"][1]
            writer.writerow(row)


def _fmt(value: Any, spec: str = ".3f") -> str:
    if value is None:
        return "n/a"
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return "n/a"
    return format(value, spec)


def write_report(path: Path, report: dict[str, Any]) -> None:
    lines: list[str] = []
    lines.append("# Joint nested-signature enrichment (multi-sample cohort)")
    lines.append("")
    lines.append(
        "Generated by `scripts/joint_enrichment.py`. Only the pre-specified cells "
        "listed in that script are tested; the full family x position x "
        "orientation table is not scanned."
    )
    lines.append("")

    lines.append("## Headline")
    lines.append("")
    primary = report["cells"][0]
    lines.append(
        f"1. **Primary cell** ({primary['label']}): observed {primary['observed']} "
        f"against an expected { _fmt(primary['expected_null_mean'], '.2f') } "
        f"(analytic { _fmt(primary['expected_analytic'], '.2f') }), effect size "
        f"{_fmt(primary['effect_size_observed_over_expected'], '.2f')}x "
        f"(host-clustered 95% CI "
        f"{_fmt(primary['bootstrap_ci_low'], '.2f')} to "
        f"{_fmt(primary['bootstrap_ci_high'], '.2f')}), Monte Carlo p = "
        f"{_fmt(primary['monte_carlo_p'], '.4g')} over "
        f"{primary['replicates']:,} replicates."
    )
    lines.append(
        f"2. **Multiplication test.** {report['multiplication_test']['verdict']}. "
        f"Benchmarks: position factor "
        f"{_fmt(report['multiplication_test'].get('position_factor'), '.2f')}x, "
        f"orientation factor "
        f"{_fmt(report['multiplication_test'].get('orientation_factor'), '.2f')}x. "
        "These two numbers exist only as the comparison benchmark and are not an "
        "effect size; no multiplied figure is reported anywhere in this document."
    )
    lines.append(
        "3. **Secondary cells, Holm-adjusted within the Alu family:** "
        + "; ".join(
            f"{c['cell_id']} p={_fmt(c.get('holm_adjusted_p'), '.4g')}"
            for c in report["cells"][1:]
        )
        + "."
    )
    lines.append("")

    lines.append("## Opportunity null")
    lines.append("")
    opp = report["opportunity_diagnostic"]
    if opp.get("status") == "measured":
        lines.append(
            f"Measured on {opp['hosts_measured']} Alu host elements. The residual "
            f"mask fraction from repeats *other* than the host's own annotation has "
            f"median {_fmt(opp['residual_mask_fraction_median'], '.4f')} and maximum "
            f"{_fmt(opp['residual_mask_fraction_max'], '.4f')}; "
            f"{opp['hosts_with_zero_opportunity_from_other_repeats']} of "
            f"{opp['hosts_measured']} hosts have full opportunity, and "
            f"{opp['hosts_fully_masked_by_other_repeats']} are fully masked. "
            f"Verdict: `{opp['verdict']}`."
        )
    else:
        lines.append(f"Status: `{opp.get('status')}`. {opp.get('consequence', '')}")
    lines.append("")
    lines.append(
        "The frozen plan ranks the host-interval null as primary and a "
        "mappability-scaled null as sensitivity only, because a short-read track "
        "is not this cohort's callability. No mappability or gap track exists in "
        "this workspace, so the primary null is used and the substitute above is "
        "what the sensitivity variant reduces to. This is recorded as a deviation."
    )
    lines.append("")

    lines.append("## Cells")
    lines.append("")
    lines.append(
        "| cell | hosts w/ opportunity | structural-zero hosts | sites | carriers | "
        "observed | expected | effect | 95% CI | MC p | Holm p |"
    )
    lines.append("|---|---|---|---|---|---|---|---|---|---|---|")
    for cell in report["cells"]:
        unstable = cell.get("effect_size_unstable_near_zero_expected")
        effect_text = (
            "unstable*" if unstable
            else f"{_fmt(cell['effect_size_observed_over_expected'], '.2f')}x"
        )
        lines.append(
            f"| `{cell['cell_id']}` | {cell['hosts_with_opportunity']} "
            f"| {cell['hosts_structural_zero_excluded']} "
            f"| {cell['sites_in_scope']} "
            f"| {cell['carrier_observations_in_scope']} "
            f"| {cell['observed']} "
            f"| {_fmt(cell['expected_null_mean'], '.2f')} "
            f"| {effect_text} "
            f"| {_fmt(cell['bootstrap_ci_low'], '.2f')}–{_fmt(cell['bootstrap_ci_high'], '.2f')} "
            f"| {_fmt(cell['monte_carlo_p'], '.4g')} "
            f"| {_fmt(cell.get('holm_adjusted_p'), '.4g')} |"
        )
    if any(c.get("effect_size_unstable_near_zero_expected") for c in report["cells"]):
        lines.append("")
        lines.append(
            "\\* the expected count for this cell is within 1 of zero, so the "
            "observed/expected ratio and its interval are unstable and should not "
            "be quoted. The difference (observed minus expected) is the "
            "interpretable quantity for these cells."
        )
    lines.append("")

    lines.append("## Private sites (exploratory, no p-value)")
    lines.append("")
    lines.append(
        "Per the plan's estimand statement, private (single-carrier) sites are a "
        "separate exploratory analysis and no primary-test p-value is computed "
        "over them. Counts are given for transparency only."
    )
    lines.append("")
    lines.append("| cell | private sites | observed | expected |")
    lines.append("|---|---|---|---|")
    for cell in report["cells"]:
        lines.append(
            f"| `{cell['cell_id']}` | {cell['private_exploratory_sites']} "
            f"| {cell['private_exploratory_observed']} "
            f"| {_fmt(cell['private_exploratory_expected'], '.2f')} |"
        )
    lines.append("")

    lines.append("## Method")
    lines.append("")
    for item in report["method_notes"]:
        lines.append(f"- {item}")
    lines.append("")

    lines.append("## Deviations")
    lines.append("")
    for item in report["deviations"]:
        lines.append(f"- {item}")
    lines.append("")
    path.write_text("\n".join(lines) + "\n")


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--unique-sites", type=Path, default=common.DEFAULT_UNIQUE_SITES)
    p.add_argument("--callset-dir", type=Path, default=common.DEFAULT_CALLSET_DIR)
    p.add_argument("--rmsk", type=Path, default=DEFAULT_RMSK)
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--replicates", type=int, default=DEFAULT_REPLICATES)
    p.add_argument("--bootstrap", type=int, default=DEFAULT_BOOTSTRAP)
    p.add_argument("--seed", type=int, default=DEFAULT_SEED)
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    args.outdir.mkdir(parents=True, exist_ok=True)

    sites, load_report = common.load_unique_sites(args.unique_sites)
    # An empty cohort is a broken input, not a null result. Both this script and
    # `recurrence_test.py` would otherwise go on to compute a null over nothing
    # and report it as "no effect". The loader stays permissive because an
    # all-hostless table is a legitimate thing to hand it.
    if not sites:
        raise SystemExit(
            f"{args.unique_sites} yielded no usable sites: "
            f"{json.dumps(load_report)}. Refusing to report an empty cohort as a "
            f"null result."
        )

    print(f"unique sites in scope: {len(sites)}", flush=True)
    samples = sorted({s for site in sites for s in site["carriers"]})
    callsets = common.load_callsets(args.callset_dir, samples)
    sites = common.attach_call_details(sites, callsets)
    join = common.verify_join(sites)
    print(f"[join] {join['verdict']} ({join['calls_recovered']} calls)", flush=True)
    if join["verdict"] != "join_consistent_with_dedup_output":
        raise SystemExit(
            "per-call join disagrees with the dedup output; refusing to continue: "
            f"{json.dumps(join)}"
        )

    census = common.nested_census_for_sites(sites)
    print(f"NESTED census: {census}", flush=True)

    shared = [s for s in sites if not s["private"]]
    private = [s for s in sites if s["private"]]
    print(f"shared sites: {len(shared)}  private sites: {len(private)}", flush=True)

    cells: list[dict[str, Any]] = []
    for cell in CELLS:
        result = run_cell(
            cell, shared, private, args.replicates, args.bootstrap, args.seed
        )
        cells.append(result)
        if result["expected_null_mean"] is not None and result["expected_analytic"]:
            # The simulated mean and the closed form are independent routes to the
            # same expectation. A large gap means one of them is wrong, so the gap
            # is reported rather than only the prettier of the two.
            result["expectation_agreement"] = abs(
                result["expected_null_mean"] - result["expected_analytic"]
            )
        else:
            result["expectation_agreement"] = None
        print(
            f"  {cell.cell_id}: obs={result['observed']} "
            f"exp={_fmt(result['expected_null_mean'], '.2f')} "
            f"p={_fmt(result['monte_carlo_p'], '.4g')}",
            flush=True,
        )

    secondary_ps = [c.get("monte_carlo_p") for c in cells if c["role"] == "secondary"]
    adjusted = holm(secondary_ps)
    it = iter(adjusted)
    for cell in cells:
        if cell["role"] == "secondary":
            cell["holm_adjusted_p"] = next(it)

    primary = cells[0]
    primary_cell = CELLS[0]
    factors = per_factor_estimates(shared, primary_cell, args.replicates, args.seed)
    obs_ratio = primary["effect_size_observed_over_expected"]
    bench = factors.get("benchmark_product")
    # The dedicated interaction cells are the direct test of the same question.
    # Reporting the benchmark comparison without them beside it invites reading a
    # small numerical excess as a detected interaction.
    interaction_ps = {
        c["cell_id"]: c.get("holm_adjusted_p")
        for c in cells
        if c.get("statistic", "").startswith("n_sense_in_bin -")
    }
    interaction_null = [
        k for k, v in interaction_ps.items() if v is not None and v > 0.05
    ]
    if obs_ratio is None or bench is None:
        verdict = "not estimable: the joint estimate or a benchmark factor is undefined"
        exceeds = None
    else:
        exceeds = bool(obs_ratio > bench)
        base = (
            "the joint enrichment EXCEEDS what the two per-factor estimates "
            "predict"
            if exceeds
            else "the joint enrichment does NOT exceed what the two per-factor "
            "estimates predict, so the joint cell shows no super-multiplicative "
            "structure"
        )
        if exceeds and interaction_null:
            base += (
                ". That excess should NOT be read as a detected interaction: the "
                "excess is small, and the pre-specified interaction cells that test "
                "the same question directly do not reject "
                f"({', '.join(sorted(interaction_null))} after Holm adjustment)"
            )
        elif exceeds:
            base += (
                ". The pre-specified interaction cells are consistent with this, "
                "though they are not independent evidence of it"
            )
        verdict = base

    opp = opportunity_diagnostic(shared, args.rmsk)

    report: dict[str, Any] = {
        "analysis": "joint_nested_signature_enrichment",
        "plan": "docs/NESTED_ALU_ANALYSIS_PLAN.md (frozen 2026-10-01); this analysis extends it",
        "inputs": {
            "unique_sites": str(args.unique_sites),
            "callset_dir": str(args.callset_dir),
            "replicates": args.replicates,
            "bootstrap_draws": args.bootstrap,
            "seed": args.seed,
        },
        "load_report": load_report,
        "join_verification": join,
        "nested_enum_census": census,
        "nested_enum_note": (
            "the four canonical values are unnested, nested_sense, "
            "nested_antisense, nested_unknown. These callsets populate only the "
            "coarser unnested/nested pair, so the unlabeled value is kept as its "
            "own label and folded into neither orientation class. The orientation "
            "cells use the independent ORIENT field, not the NESTED label."
        ),
        "cohort_split": {
            "shared_sites": len(shared),
            "private_sites": len(private),
            "note": (
                "plan estimand statement: positional analyses use unique sites; "
                "private sites are a separate exploratory analysis and no "
                "primary-test p-value is computed over them"
            ),
        },
        "cells": cells,
        "multiplication_test": {
            "primary_cell": primary["cell_id"],
            "position_factor": factors.get("position_factor_bin_enrichment"),
            "orientation_factor": factors.get("orientation_factor_sense_vs_parity"),
            "benchmark_product": bench,
            "observed_joint_effect_size": obs_ratio,
            "joint_exceeds_benchmark": exceeds,
            "excess_over_benchmark": (
                obs_ratio / bench if obs_ratio and bench else None
            ),
            "interaction_cell_holm_p": interaction_ps,
            "verdict": verdict,
            "note": (
                "the two factors are the marginals of the joint cell and exist "
                "only as the benchmark for this comparison; no multiplied effect "
                "size is reported"
            ),
        },
        "opportunity_diagnostic": opp,
        "method_notes": [
            "host-stratified conditional randomization; no log-linear model fitted",
            "within each host the observed event count is held fixed and breakpoints "
            "are redrawn uniformly from that host's opportunity",
            "interaction cells additionally permute sense/antisense within each host, "
            "preserving that host's orientation totals",
            "position bins with zero opportunity in a host are excluded from the null "
            "for that host and never pseudocounted",
            f"Monte Carlo p = (1 + #{{sim >= obs}}) / (B + 1) with B = {args.replicates:,}",
            "effect size is observed / expected with a host-clustered percentile "
            "bootstrap CI, resampling whole hosts because events in one host are "
            "not independent",
            "GT and GQ are never read: the callset reader parses CHROM/POS/ID/INFO "
            "only and never touches the FORMAT column",
            "counting is restricted to chr1-22 and chrX",
        ],
        "deviations": [
            "no mappability or gap track exists in this workspace; the host-interval "
            "opportunity null (the plan's primary) is used and the RepeatMasker "
            "substitute for the mappability variant is reported as measured above",
            "the NESTED field in these callsets populates two of the four canonical "
            "values; the unlabeled value is reported, never folded",
            "site counts and carrier counts are reported side by side and never summed",
        ],
    }
    (args.outdir / "joint_enrichment.json").write_text(
        json.dumps(report, indent=2, sort_keys=True, default=str) + "\n"
    )
    write_cells_csv(args.outdir / "joint_cells.csv", cells)
    write_report(args.outdir / "joint_enrichment.md", report)
    print(
        f"wrote joint_enrichment.md / .json and joint_cells.csv in {args.outdir}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
