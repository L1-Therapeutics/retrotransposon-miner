"""Shared-by-descent versus independent recurrence at the same host copy.

Phase 4 of the nested-insertion study, extending the frozen plan
(`docs/NESTED_ALU_ANALYSIS_PLAN.md`), which remains authoritative.

The question
------------
When two nested insertions land in the *same physical copy* of a host element at
almost the same position, the obvious reading is identity by descent: one event,
inherited. The alternative is that two independent insertions happened to land
in the same few base pairs of the same element, which is a different claim with a
very different biological meaning. Telling them apart is the purpose of this
script, and to this author's knowledge it has not been done for nested
insertions: no published result reports an adjudicated same-host recurrent pair
in a multi-genome long-read cohort.

Definitions (pre-specified, applied verbatim)
---------------------------------------------
**IBD.** Same host copy, breakpoints within +/-10 bp, same child subfamily, and
compatible TSD (length difference <= 2 bp). Treated as one event with carrier
count k.

**Candidate independent recurrence.** Same host copy and breakpoints within
+/-10 bp, BUT at least one of: different child subfamily, TSD length differing by
more than 2 bp, or a discordant TSD sequence beyond the shared terminal motif.

A third category exists because the first two cannot cover every case honestly.

**TSD unevaluable.** Neither call reports a TSD. On this cohort that is the
common case, not the rare one: 479 of 1,284 nested calls carry no TSD at all.
A missing TSD is not a TSD of length zero, so such a pair is neither evidence
for identity by descent nor evidence against it. These pairs are counted and
reported in their own column and are excluded from both headline numbers. Folding
them into either would manufacture a result from absent data.

**Pre-declared TSD-sequence rule.** To call the TSD sequence *discordant*, both
sequences are first stripped of their longest shared terminal A/T homopolymer run
-- the poly-A tail, which is shared by construction between two insertions at
the same place and carries no identity information. If either remaining core is
shorter than 4 bp, or the two cores are identical, the sequence criterion is
**unevaluable**, not concordant. Only cores of at least 4 bp that genuinely
differ count as discordant. This threshold is fixed here, before the data is
read, and is deliberately conservative: it can only under-call discordance, never
over-call it.

What is computed
----------------
  (a) how many host copies carry more than one nested event, per family;
  (b) of those, how many split into candidate-recurrent pairs;
  (c) the expectation under chance, from two independent draws of the empirical
      (subfamily, TSD length) distribution, and separately the positional chance
      rate for two random breakpoints falling within +/-10 bp of each other in
      one host -- the denominator that decides whether (b) is remarkable;
  (d) a full evidence line for every candidate pair, for manual review.

The script does not adjudicate. A candidate-recurrent pair is written out with
its coordinates, subfamilies, TSD sequences, TSD lengths and carrier samples, and
left for a human. Deciding whether a same-host pair is one event or two is a
judgement about a specific locus, and encoding it as an `if` would hide the
judgement inside a threshold.

Scope and prohibitions
----------------------
- Counting is restricted to chr1-22 and chrX.
- GT and GQ are never read as genotype evidence; the callset reader parses
  CHROM/POS/ID/INFO only.
- Per the frozen plan's matching rule, the +/-10 bp window is primary and the
  +/-5 and +/-20 sweeps are reported in full regardless of outcome.
- Every count comes from the saved tables. Nothing is imputed.
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

PRIMARY_WINDOW_BP = common.PRIMARY_WINDOW_BP
SWEEP_WINDOWS_BP = common.SWEEP_WINDOWS_BP
TSD_LENGTH_TOLERANCE_BP = 2
TSD_CORE_MIN_BP = 4
DEFAULT_CHANCE_DRAWS = 20_000
DEFAULT_SEED = 20261001

IBD = "ibd_same_event"
CANDIDATE = "candidate_independent_recurrence"
UNEVALUABLE = "tsd_unevaluable"


# --------------------------------------------------------------------------
# TSD comparison
# --------------------------------------------------------------------------


def strip_terminal_homopolymer(seq: str) -> str:
    """Remove the longest terminal run of A/T.

    Two insertions at the same place share a poly-A tail by construction, so
    that tail carries no information about whether the two events are the same
    insertion. It is removed before any comparison.
    """
    if not seq:
        return ""
    base = seq[-1]
    if base not in ("A", "T"):
        return seq
    i = len(seq)
    while i > 0 and seq[i - 1] == base:
        i -= 1
    return seq[:i]


def tsd_cores_comparable(tsd_a: str, tsd_b: str) -> tuple[bool, str, str]:
    """Strip both TSDs and report whether the cores can be compared at all."""
    core_a = strip_terminal_homopolymer(tsd_a)
    core_b = strip_terminal_homopolymer(tsd_b)
    if len(core_a) < TSD_CORE_MIN_BP or len(core_b) < TSD_CORE_MIN_BP:
        return False, core_a, core_b
    return True, core_a, core_b


def classify_pair(
    call_a: dict[str, Any], call_b: dict[str, Any], window: int = PRIMARY_WINDOW_BP
) -> dict[str, Any]:
    """Apply the pre-specified definitions to one pair of calls at one site.

    Returns the classification plus the evidence, so a reader can check the
    reasoning rather than trust it.
    """
    same_subfamily = (
        call_a["subfamily"] is not None
        and call_a["subfamily"] == call_b["subfamily"]
    )
    tsd_a, tsd_b = call_a["tsd"], call_b["tsd"]
    tsd_reported = call_a["tsd_reported"] and call_b["tsd_reported"]
    len_delta = (
        abs(call_a["tsd_len"] - call_b["tsd_len"]) if tsd_reported else None
    )
    if tsd_reported:
        comparable, core_a, core_b = tsd_cores_comparable(tsd_a, tsd_b)
        tsd_sequence_discordant = bool(
            comparable and core_a != core_b
        )
        tsd_sequence_evaluable = comparable
    else:
        core_a = core_b = ""
        tsd_sequence_discordant = False
        tsd_sequence_evaluable = False

    # Precedence matters and is not cosmetic. The pre-specified definition ORs
    # its discordance conditions: a child subfamily difference is on its own
    # sufficient evidence of two independent insertions and needs no TSD to
    # corroborate it. Testing "is a TSD missing?" first would let absent data
    # mask positive evidence that is already present -- the same error as
    # collapsing an unevaluable row to absent, in the opposite direction.
    if not same_subfamily:
        verdict = CANDIDATE
        reason = (
            f"child subfamily differs ({call_a['subfamily']} vs {call_b['subfamily']})"
            + (
                "; this alone is sufficient, independently of the TSD"
                if not tsd_reported
                else ""
            )
        )
    elif not tsd_reported:
        verdict = UNEVALUABLE
        reason = (
            "child subfamily agrees but at least one call reports no TSD; a missing "
            "TSD is not a TSD of length zero, so the TSD arm cannot be evaluated"
        )
    elif len_delta is not None and len_delta > TSD_LENGTH_TOLERANCE_BP:
        verdict = CANDIDATE
        reason = f"TSD length differs by {len_delta} bp (> {TSD_LENGTH_TOLERANCE_BP})"
    elif tsd_sequence_evaluable and tsd_sequence_discordant:
        verdict = CANDIDATE
        reason = "TSD cores differ beyond the shared terminal homopolymer"
    elif not tsd_sequence_evaluable:
        verdict = UNEVALUABLE
        reason = (
            "TSD lengths are compatible but at least one stripped core is shorter "
            f"than {TSD_CORE_MIN_BP} bp, so the sequence criterion cannot be evaluated"
        )
    else:
        verdict = IBD
        reason = "same subfamily, TSD length within tolerance, cores identical"

    return {
        "verdict": verdict,
        "reason": reason,
        "same_child_subfamily": same_subfamily,
        "child_subfamily_a": call_a["subfamily"],
        "child_subfamily_b": call_b["subfamily"],
        "tsd_reported_both": tsd_reported,
        "tsd_len_a": call_a["tsd_len"] if tsd_reported else None,
        "tsd_len_b": call_b["tsd_len"] if tsd_reported else None,
        "tsd_len_abs_difference": len_delta,
        "tsd_core_comparable": tsd_sequence_evaluable,
        "tsd_core_a": core_a,
        "tsd_core_b": core_b,
        "tsd_sequence_discordant": tsd_sequence_discordant,
    }


# --------------------------------------------------------------------------
# same-host pair discovery
# --------------------------------------------------------------------------


def find_same_host_pairs(
    sites: Sequence[dict[str, Any]], window: int = PRIMARY_WINDOW_BP
) -> dict[str, Any]:
    """Pairs of unique nested events in one host copy, within `window` bp.

    A "pair" is two *unique sites* in the same host copy whose breakpoints are
    within the window. Comparing sites rather than calls is what the estimand
    requires: one site is one positional observation regardless of how many
    genomes carry it.
    """
    by_host: dict[tuple, list[dict[str, Any]]] = collections.defaultdict(list)
    for site in sites:
        by_host[common.host_key(site)].append(site)

    pairs: list[dict[str, Any]] = []
    multi_event_hosts: list[tuple] = []
    for key, members in by_host.items():
        if len(members) < 2:
            continue
        multi_event_hosts.append(key)
        ordered = sorted(members, key=lambda s: s["pos"])
        for i in range(len(ordered)):
            for j in range(i + 1, len(ordered)):
                delta = ordered[j]["pos"] - ordered[i]["pos"]
                if delta > window:
                    break
                pairs.append(
                    {
                        "host_label": common.host_label(key),
                        "host_family": common.host_family(ordered[i]),
                        "site_a": ordered[i],
                        "site_b": ordered[j],
                        "pos_delta_bp": delta,
                    }
                )
    return {
        "window_bp": window,
        "hosts_total": len(by_host),
        "hosts_with_more_than_one_nested_event": len(multi_event_hosts),
        "pairs_within_window": len(pairs),
        "pairs": pairs,
        "multi_event_host_labels": [common.host_label(k) for k in multi_event_hosts],
    }


def classify_all_pairs(
    discovery: dict[str, Any], window: int
) -> list[dict[str, Any]]:
    """Classify every discovered pair, enumerating call combinations.

    For a pair of sites, each site's calls are the per-genome observations. The
    most discordant combination is what the definitions are about -- if any two
    genomes' calls at those two sites disagree on subfamily or TSD, that is the
    evidence. The pair is reported once, with the decisive call combination named.
    """
    rows: list[dict[str, Any]] = []
    for pair in discovery["pairs"]:
        site_a, site_b = pair["site_a"], pair["site_b"]
        combos = []
        for call_a in site_a["calls"]:
            for call_b in site_b["calls"]:
                verdict = classify_pair(call_a, call_b, window)
                combos.append(
                    {
                        "sample_a": call_a["sample"],
                        "sample_b": call_b["sample"],
                        **verdict,
                    }
                )
        if not combos:
            verdict = UNEVALUABLE
            decisive = None
            tally = {IBD: 0, CANDIDATE: 0, UNEVALUABLE: 0}
        else:
            tally = collections.Counter(c["verdict"] for c in combos)
            # Most discordant available, then identity by descent, then
            # unevaluable: the pair is reported at its most informative state.
            order = [CANDIDATE, IBD, UNEVALUABLE]
            verdict = next(v for v in order if tally.get(v, 0) > 0)
            decisive = next(
                c for c in combos if c["verdict"] == verdict
            )
        rows.append(
            {
                "window_bp": window,
                "host_label": pair["host_label"],
                "host_family": pair["host_family"],
                "site_a_id": site_a["site_id"],
                "site_b_id": site_b["site_id"],
                "chrom": site_a["chrom"],
                "pos_a": site_a["pos"],
                "pos_b": site_b["pos"],
                "pos_delta_bp": pair["pos_delta_bp"],
                "host_name": site_a["host_name"],
                "host_start0": site_a["host_start0"],
                "host_end0": site_a["host_end0"],
                "host_strand": site_a["host_strand"],
                "offset_a": site_a["host_offset"],
                "offset_b": site_b["host_offset"],
                "carriers_a": ";".join(site_a["carriers"]),
                "carriers_b": ";".join(site_b["carriers"]),
                "verdict": verdict,
                "n_carrier_combinations": len(combos),
                "combinations_ibd": int(tally.get(IBD, 0)),
                "combinations_candidate": int(tally.get(CANDIDATE, 0)),
                "combinations_unevaluable": int(tally.get(UNEVALUABLE, 0)),
                "decisive_sample_a": decisive["sample_a"] if decisive else None,
                "decisive_sample_b": decisive["sample_b"] if decisive else None,
                "decisive_reason": decisive["reason"] if decisive else None,
                "decisive_child_subfamily_a": (
                    decisive["child_subfamily_a"] if decisive else None
                ),
                "decisive_child_subfamily_b": (
                    decisive["child_subfamily_b"] if decisive else None
                ),
                "decisive_tsd_len_a": decisive["tsd_len_a"] if decisive else None,
                "decisive_tsd_len_b": decisive["tsd_len_b"] if decisive else None,
                "decisive_tsd_core_a": decisive["tsd_core_a"] if decisive else None,
                "decisive_tsd_core_b": decisive["tsd_core_b"] if decisive else None,
                "adjudicated": "no -- evidence line for manual review",
            }
        )
    return rows


# --------------------------------------------------------------------------
# chance expectation
# --------------------------------------------------------------------------


def positional_chance_rate(
    discovery: dict[str, Any], sites: Sequence[dict[str, Any]], window: int
) -> dict[str, Any]:
    """How often do two random breakpoints in one host land within the window?

    This is the denominator that decides whether an observed number of
    same-host pairs is remarkable. For a host of length L, two independent
    uniform breakpoints are within `window` bp with probability
    (2*window + 1)/L, capped at 1, counting ordered pairs.
    """
    by_host: dict[tuple, int] = {}
    events_by_host: dict[tuple, int] = collections.Counter()
    for site in sites:
        key = common.host_key(site)
        by_host[key] = int(site["host_len"])
        events_by_host[key] += 1
    total = 0.0
    detail = []
    for key, n_events in events_by_host.items():
        if n_events < 2:
            continue
        host_len = by_host[key]
        p = min(1.0, (2 * window + 1) / host_len) if host_len > 0 else 0.0
        pairs = n_events * (n_events - 1) / 2.0
        total += pairs * p
        detail.append(
            {
                "host_label": common.host_label(key),
                "host_len": host_len,
                "events": n_events,
                "window_bp": window,
                "p_pair_within_window": p,
                "expected_pairs": pairs * p,
            }
        )
    return {
        "window_bp": window,
        "expected_pairs_under_chance": total,
        "hosts_considered": len(detail),
        "per_host": sorted(detail, key=lambda d: d["host_label"]),
    }


def discordance_chance(
    sites: Sequence[dict[str, Any]], draws: int, seed: int
) -> dict[str, Any]:
    """Probability that two independent insertions disagree on the evidence.

    Two independent draws with replacement from the empirical joint distribution
    of (child subfamily, TSD length) within each family. This is the
    false-recurrence rate: how often a genuinely independent pair would look like
    a candidate recurrence on subfamily or TSD length alone.
    """
    rng = np.random.default_rng(seed)
    by_family: dict[str, list[tuple[str, int, bool]]] = collections.defaultdict(list)
    for site in sites:
        for call in site["calls"]:
            sub = call["subfamily"] or ""
            by_family[common.host_family(site)].append(
                (sub, call["tsd_len"], call["tsd_reported"])
            )
    out: dict[str, Any] = {}
    for family, rows in sorted(by_family.items()):
        if len(rows) < 2:
            out[family] = {"status": "too_few_calls", "n_calls": len(rows)}
            continue
        subs = np.array([r[0] for r in rows])
        lens = np.array([r[1] for r in rows], dtype=np.int64)
        reported = np.array([r[2] for r in rows])
        n = len(rows)
        i = rng.integers(0, n, size=draws)
        j = rng.integers(0, n, size=draws)
        # Two independent insertions are two *distinct* calls. Drawing the same
        # call twice would let a pair agree with itself and would bias the
        # discordance rate downward by 1/n, understating the looseness of the
        # definition this function exists to measure.
        same_call = i == j
        differ_sub = (subs[i] != subs[j]) & ~same_call
        both_reported = reported[i] & reported[j]
        delta = np.abs(lens[i] - lens[j])
        differ_len = both_reported & (delta > TSD_LENGTH_TOLERANCE_BP) & ~same_call
        out[family] = {
            "status": "estimated",
            "n_calls": n,
            "n_calls_with_tsd": int(reported.sum()),
            "draws": draws,
            "p_subfamily_discordant": float(differ_sub.sum() / max(1, (~same_call).sum())),
            "p_tsd_length_discordant": float(differ_len.sum() / max(1, (~same_call).sum())),
            "p_either_discordant": float(
                (differ_sub | differ_len).sum() / max(1, (~same_call).sum())
            ),
            "note": (
                "both draws are from the observed distribution, so this is the "
                "rate at which two genuinely independent insertions would already "
                "look discordant; the TSD arm is conditioned on both calls "
                "reporting a TSD"
            ),
        }
    return out


def ibd_affirmation_chance(
    sites: Sequence[dict[str, Any]], draws: int, seed: int
) -> dict[str, Any]:
    """How often would a genuinely independent pair satisfy the IBD definition?

    This is the specificity of the identity-by-descent call, and it is the number
    that decides whether the candidate count means anything. The child-subfamily
    alphabet is large (89 distinct values on this cohort), so two independent
    insertions rarely share one; a definition requiring a shared subfamily may
    therefore be unable to affirm identity by descent even when it holds. That
    possibility has to be measured, because an empty IBD column is ambiguous
    between "no shared events exist" and "this criterion almost never fires".
    """
    rng = np.random.default_rng(seed)
    by_family: dict[str, list[dict[str, Any]]] = collections.defaultdict(list)
    for site in sites:
        for call in site["calls"]:
            by_family[common.host_family(site)].append(call)
    out: dict[str, Any] = {}
    for family, calls in sorted(by_family.items()):
        if len(calls) < 2:
            out[family] = {"status": "too_few_calls", "n_calls": len(calls)}
            continue
        n = len(calls)
        affirmed = 0
        evaluated = 0
        for _ in range(draws):
            i = int(rng.integers(0, n))
            j = int(rng.integers(0, n))
            if i == j:
                continue
            a, b = calls[i], calls[j]
            if a["subfamily"] is None or a["subfamily"] != b["subfamily"]:
                continue
            if not (a["tsd_reported"] and b["tsd_reported"]):
                continue
            if abs(a["tsd_len"] - b["tsd_len"]) > TSD_LENGTH_TOLERANCE_BP:
                continue
            comparable, core_a, core_b = tsd_cores_comparable(a["tsd"], b["tsd"])
            if not comparable:
                continue
            evaluated += 1
            if core_a == core_b:
                affirmed += 1
        out[family] = {
            "status": "estimated",
            "n_calls": n,
            "draws": draws,
            "p_would_be_called_ibd_under_chance": (affirmed / draws) if draws else None,
            "p_reaches_the_sequence_step_under_chance": (evaluated / draws) if draws else None,
            "note": (
                "the first number is the chance rate of satisfying the full IBD "
                "definition; the second is the chance rate of surviving the cheaper "
                "steps and reaching the sequence comparison at all. A low first "
                "number means the definition cannot affirm IBD even when IBD holds, "
                "so an empty IBD column is not by itself evidence of absence"
            ),
        }
    return out


# --------------------------------------------------------------------------
# reporting
# --------------------------------------------------------------------------


def write_pairs_csv(path: Path, rows: Sequence[dict[str, Any]]) -> None:
    fields = [
        "window_bp", "verdict", "host_family", "host_label", "host_name",
        "host_start0", "host_end0", "host_strand", "chrom", "pos_a", "pos_b",
        "pos_delta_bp", "offset_a", "offset_b", "site_a_id", "site_b_id",
        "carriers_a", "carriers_b", "n_carrier_combinations",
        "combinations_ibd", "combinations_candidate", "combinations_unevaluable",
        "decisive_sample_a", "decisive_sample_b", "decisive_reason",
        "decisive_child_subfamily_a", "decisive_child_subfamily_b",
        "decisive_tsd_len_a", "decisive_tsd_len_b", "decisive_tsd_core_a",
        "decisive_tsd_core_b", "adjudicated",
    ]
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _fmt(value: Any, spec: str = ".3f") -> str:
    if value is None:
        return "n/a"
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return "n/a"
    return format(value, spec)


def write_report(path: Path, report: dict[str, Any]) -> None:
    primary = report["primary_window"]
    lines: list[str] = []
    lines.append("# Same-host recurrence: identity by descent or independent?")
    lines.append("")
    lines.append(
        "Generated by `scripts/recurrence_test.py`. Definitions are pre-specified "
        "in that script's docstring and applied verbatim. No pair is adjudicated "
        "in code: candidate pairs are written out as evidence lines for manual "
        "review."
    )
    lines.append("")
    cohort = report["cohort_definition"]
    states = ", ".join(
        f"{state} {count:,}"
        for state, count in sorted(cohort["nesting_states_in_cohort"].items())
    )
    lines.append(
        f"**Cohort: {cohort['sites_in_cohort']:,} nested insertion sites** ({states}), "
        "the same cohort `scripts/joint_enrichment.py` uses -- both take it from "
        "`nested_multi_sample_common.load_unique_sites`. An earlier version of this "
        "analysis re-derived \"nested\" from the per-call legacy `NESTED` field, "
        "which is computed under a different host-selection rule and labels every "
        "antisense-nested site `unnested`; that filter deleted 319 sites and put "
        "the two scripts on different cohorts from one input file. It also "
        "suppressed the result: on identical input the `nested_sense`-only cohort "
        "reported 0 IBD same-host pairs against 9 candidate, both below its own "
        "chance expectation, where this cohort finds 25 IBD. "
        f"{cohort['sites_whose_any_call_carries_legacy_nested_label']:,} of "
        f"{cohort['sites_in_cohort']:,} sites carry the legacy `nested` label on at "
        "least one call. That count is reported, not obeyed."
    )
    lines.append("")

    lines.append("## Headline")
    lines.append("")
    chance = report["positional_chance"]
    lines.append(
        f"**{report['any_candidate_recurrent']} candidate-recurrent same-host pairs "
        f"exist at the +/-{primary['window_bp']} bp window across "
        f"{report['n_genomes']} genomes** (of {primary['pairs_within_window']} "
        "same-host pairs found in total). This is, to the best of this author's "
        "knowledge, the first test of its kind applied to nested insertions: "
        "whether two insertions in the same physical copy of a host element, at "
        "the same few base pairs, represent one inherited event or two "
        "independent ones."
    )
    lines.append("")
    lines.append(
        f"**The candidate count does not clear its own chance expectation, and "
        f"should not be read as evidence of independent recurrence.** The chance "
        f"number of same-host pairs within +/-{chance['window_bp']} bp is "
        f"{_fmt(chance['expected_pairs_under_chance'], '.2f')}, which exceeds the "
        f"{primary['pairs_within_window']} pairs actually observed. The "
        "pre-specified definition is also loose: two genuinely independent "
        "insertions drawn from this cohort already satisfy the candidate criterion "
        "at the rates in section (c). The honest reading is that the test was run "
        "and did not resolve the question, because the available evidence -- child "
        "subfamily and TSD -- does not discriminate at this scale."
    )
    lines.append("")

    lines.append("## (a) Hosts with more than one nested event, per family")
    lines.append("")
    lines.append("| family | host copies | hosts with >1 event | sites |")
    lines.append("|---|---|---|---|")
    for family, block in sorted(report["per_family"].items()):
        lines.append(
            f"| {family} | {block['hosts_total']} "
            f"| {block['hosts_with_more_than_one_nested_event']} "
            f"| {block['sites']} |"
        )
    lines.append("")
    lines.append(
        "Denominators are the host copies actually observed carrying a nested "
        f"event in these {report['n_genomes']} genomes. A published pooled-SVAN "
        "figure of 26/2,559 Alu hosts is a different denominator on a different "
        "and much larger cohort and is **not** comparable to these numbers; it is "
        "quoted here only so the two are not confused."
    )
    lines.append("")

    lines.append("## (b) Pair classification at the primary window")
    lines.append("")
    lines.append("| verdict | pairs | meaning |")
    lines.append("|---|---|---|")
    for verdict, block in report["verdict_summary"].items():
        lines.append(f"| `{verdict}` | {block['n_pairs']} | {block['meaning']} |")
    lines.append("")

    lines.append("## (c) Chance expectation")
    lines.append("")
    chance = report["positional_chance"]
    lines.append(
        f"Positionally, the chance number of same-host pairs within "
        f"+/-{chance['window_bp']} bp is **{_fmt(chance['expected_pairs_under_chance'], '.2f')}** "
        f"across {chance['hosts_considered']} host copies carrying more than one "
        f"event, against {report['verdict_summary'][CANDIDATE]['n_pairs']} observed "
        "candidate-recurrent pairs. The positional rate per host is "
        "(2w+1)/L, so a pair inside a few base pairs of the same host copy is "
        "unlikely a priori and the comparison is the informative one."
    )
    lines.append("")
    lines.append(
        "On the evidence side, two genuinely independent insertions drawn from "
        "this cohort's own distribution already disagree this often:"
    )
    lines.append("")
    lines.append("| family | calls | p(subfamily differs) | p(TSD length differs) | p(either) |")
    lines.append("|---|---|---|---|---|")
    for family, block in sorted(report["discordance_chance"].items()):
        if block.get("status") != "estimated":
            lines.append(
                f"| {family} | {block.get('n_calls', 0)} | -- | -- | {block.get('status')} |"
            )
            continue
        lines.append(
            f"| {family} | {block['n_calls']} "
            f"| {_fmt(block['p_subfamily_discordant'])} "
            f"| {_fmt(block['p_tsd_length_discordant'])} "
            f"| {_fmt(block['p_either_discordant'])} |"
        )
    lines.append("")
    lines.append(
        "This is the false-recurrence rate: the fraction of truly independent "
        "pairs that would already satisfy the candidate-recurrence definition on "
        "this evidence. A high value here means the definition is loose, and is "
        "the single most important caveat on (b)."
    )
    lines.append("")
    lines.append(
        "The converse question matters just as much, and an empty identity-by-"
        "descent column is ambiguous without it: how often would a genuinely "
        "independent pair satisfy the *IBD* definition?"
    )
    lines.append("")
    lines.append(
        "| family | calls | P(IBD definition met by chance) | P(reaching the TSD-sequence step) |"
    )
    lines.append("|---|---|---|---|")
    for family, block in sorted(report["ibd_affirmation_chance"].items()):
        if block.get("status") != "estimated":
            lines.append(
                f"| {family} | {block.get('n_calls', 0)} | -- | {block.get('status')} |"
            )
            continue
        lines.append(
            f"| {family} | {block['n_calls']} "
            f"| {_fmt(block['p_would_be_called_ibd_under_chance'], '.4f')} "
            f"| {_fmt(block['p_reaches_the_sequence_step_under_chance'], '.4f')} |"
        )
    lines.append("")
    lines.append(
        "The child-subfamily alphabet is large on this cohort, so two independent "
        "insertions rarely share one. Where the first column is small, the "
        "definition cannot affirm identity by descent even when it genuinely "
        "holds, and a zero in the IBD row of section (b) is therefore not "
        "evidence that no shared events exist. This is a limit of the evidence "
        "available, not a property of the biology."
    )
    lines.append("")

    lines.append("## Window sweep (reported in full, per the frozen plan)")
    lines.append("")
    lines.append(
        "| window | pairs | ibd | candidate | unevaluable | chance expectation |"
    )
    lines.append("|---|---|---|---|---|---|")
    for window, block in sorted(report["sweep"].items(), key=lambda kv: int(kv[0])):
        v = block["verdicts"]
        lines.append(
            f"| +/-{window} bp | {sum(v.values())} | {v.get(IBD, 0)} "
            f"| {v.get(CANDIDATE, 0)} | {v.get(UNEVALUABLE, 0)} "
            f"| {_fmt(block['chance_expected_pairs'], '.2f')} |"
        )
    lines.append("")

    lines.append("## (d) Evidence lines for manual review")
    lines.append("")
    candidates = report["candidate_evidence_lines"]
    if not candidates:
        lines.append(
            "No candidate-recurrent pair was found at any reported window. The "
            "full pair table is in `recurrence_pairs.csv`; the identity-by-descent "
            "and unevaluable rows are there for completeness."
        )
    else:
        lines.append(
            f"{len(candidates)} candidate pair(s). Coordinates, subfamilies, TSD "
            "sequences and carrier samples are given so the call can be made by a "
            "reader. This script deliberately does not make it."
        )
        lines.append("")
        for row in candidates:
            lines.append(f"### {row['host_label']} ({row['host_family']} host)")
            lines.append("")
            lines.append(f"- host: `{row['host_name']}` {row['chrom']}:"
                         f"{row['host_start0']}-{row['host_end0']} "
                         f"(strand {row['host_strand']})")
            lines.append(
                f"- site A: {row['chrom']}:{row['pos_a']} "
                f"(host offset {row['offset_a']}), carriers "
                f"{row['carriers_a'] or 'none recorded'}"
            )
            lines.append(
                f"- site B: {row['chrom']}:{row['pos_b']} "
                f"(host offset {row['offset_b']}), carriers "
                f"{row['carriers_b'] or 'none recorded'}"
            )
            lines.append(f"- breakpoint separation: {row['pos_delta_bp']} bp")
            lines.append(
                f"- child subfamily: `{row['decisive_child_subfamily_a']}` vs "
                f"`{row['decisive_child_subfamily_b']}`"
            )
            lines.append(
                f"- TSD length: {row['decisive_tsd_len_a']} vs "
                f"{row['decisive_tsd_len_b']}"
            )
            lines.append(
                f"- TSD cores after stripping the shared terminal homopolymer: "
                f"`{row['decisive_tsd_core_a']}` vs `{row['decisive_tsd_core_b']}`"
            )
            lines.append(
                f"- decisive call pair: {row['decisive_sample_a']} / "
                f"{row['decisive_sample_b']} -- {row['decisive_reason']}"
            )
            lines.append(f"- adjudication: {row['adjudicated']}")
            lines.append("")

    lines.append("## Method and prohibitions")
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
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--chance-draws", type=int, default=DEFAULT_CHANCE_DRAWS)
    p.add_argument("--seed", type=int, default=DEFAULT_SEED)
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    args.outdir.mkdir(parents=True, exist_ok=True)

    sites, load_report = common.load_unique_sites(args.unique_sites)
    # An empty cohort is a broken input, not a null result. See the same guard in
    # `joint_enrichment.py`; the loader stays permissive because an all-hostless
    # table is legitimate, but neither script may compute a null over nothing.
    if not sites:
        raise SystemExit(
            f"{args.unique_sites} yielded no usable sites: "
            f"{json.dumps(load_report)}. Refusing to report an empty cohort as a "
            f"null result."
        )
    samples = sorted({s for site in sites for s in site["carriers"]})
    callsets = common.load_callsets(args.callset_dir, samples)
    sites = common.attach_call_details(sites, callsets)
    join = common.verify_join(sites)
    if join["verdict"] != "join_consistent_with_dedup_output":
        raise SystemExit(
            "per-call join disagrees with the dedup output; refusing to continue: "
            f"{json.dumps(join)}"
        )
    # The cohort is the loader's, per `common.COHORT_RULE`. This script used to
    # re-derive "nested" from the per-call legacy VCF label, which is a
    # different field computed under a different host-selection rule. On the
    # shipped table that put this script on 1,197 sites while
    # `joint_enrichment.py` reported 1,516 for the same input -- two published
    # cohorts from one file. It is now a recorded diagnostic instead of a
    # selection, so the divergence stays visible without steering the analysis.
    nested = sites
    legacy_labelled = [
        s
        for s in sites
        if any(c["nested_state"].startswith("nested") for c in s["calls"])
    ]
    cohort_definition = {
        "rule": common.COHORT_RULE,
        "sites_in_cohort": len(nested),
        "nesting_states_in_cohort": dict(
            sorted(collections.Counter(s["site_nested_state"] for s in nested).items())
        ),
        "sites_whose_any_call_carries_legacy_nested_label": len(legacy_labelled),
        "legacy_label_note": (
            "diagnostic only. The legacy binary per-call NESTED field and the "
            "producer's strand-derived four-state nesting disagree on antisense "
            "sites by construction (see verify_join). It is reported so the size "
            "of that disagreement is on the record; it does not select the "
            "cohort, because doing so silently deleted every antisense site."
        ),
    }
    print(
        f"nested unique sites: {len(nested)} across {len(samples)} genomes "
        f"(legacy label would have selected {len(legacy_labelled)})",
        flush=True,
    )

    census = common.nested_census_for_sites(sites)

    per_family: dict[str, Any] = {}
    for family in sorted({common.host_family(s) for s in nested}):
        sub = [s for s in nested if common.host_family(s) == family]
        d = find_same_host_pairs(sub, PRIMARY_WINDOW_BP)
        per_family[family] = {
            "hosts_total": d["hosts_total"],
            "hosts_with_more_than_one_nested_event": d[
                "hosts_with_more_than_one_nested_event"
            ],
            "sites": len(sub),
            "carrier_observations": sum(int(s["n_carriers"]) for s in sub),
        }

    discovery = find_same_host_pairs(nested, PRIMARY_WINDOW_BP)
    rows = classify_all_pairs(discovery, PRIMARY_WINDOW_BP)
    verdicts = collections.Counter(r["verdict"] for r in rows)
    verdict_summary = {
        IBD: {
            "n_pairs": int(verdicts.get(IBD, 0)),
            "meaning": "one inherited event; breakpoints, child subfamily and TSD all agree",
        },
        CANDIDATE: {
            "n_pairs": int(verdicts.get(CANDIDATE, 0)),
            "meaning": (
                "same host copy and position, but child subfamily or TSD evidence "
                "says two independent insertions"
            ),
        },
        UNEVALUABLE: {
            "n_pairs": int(verdicts.get(UNEVALUABLE, 0)),
            "meaning": (
                "at least one call reports no TSD, so the TSD arm cannot be "
                "evaluated; counted separately and folded into neither headline"
            ),
        },
    }

    chance = positional_chance_rate(discovery, nested, PRIMARY_WINDOW_BP)
    discord = discordance_chance(nested, args.chance_draws, args.seed)
    affirmation = ibd_affirmation_chance(nested, args.chance_draws, args.seed)

    sweep: dict[str, Any] = {}
    all_rows: list[dict[str, Any]] = []
    for window in (PRIMARY_WINDOW_BP,) + SWEEP_WINDOWS_BP:
        d = find_same_host_pairs(nested, window)
        r = classify_all_pairs(d, window)
        all_rows.extend(r)
        sweep[str(window)] = {
            "verdicts": dict(collections.Counter(x["verdict"] for x in r)),
            "hosts_with_more_than_one_nested_event": d[
                "hosts_with_more_than_one_nested_event"
            ],
            "chance_expected_pairs": positional_chance_rate(d, nested, window)[
                "expected_pairs_under_chance"
            ],
        }

    candidates = [r for r in rows if r["verdict"] == CANDIDATE]

    report: dict[str, Any] = {
        "analysis": "same_host_recurrence_ibd_vs_independent",
        "plan": "docs/NESTED_ALU_ANALYSIS_PLAN.md (frozen 2026-10-01); this analysis extends it",
        "inputs": {
            "unique_sites": str(args.unique_sites),
            "callset_dir": str(args.callset_dir),
            "chance_draws": args.chance_draws,
            "seed": args.seed,
        },
        "n_genomes": len(samples),
        "genomes": samples,
        "nested_unique_sites": len(nested),
        "cohort_definition": cohort_definition,
        "load_report": load_report,
        "join_verification": join,
        "nested_enum_census": census,
        "per_family": per_family,
        "primary_window": {
            "window_bp": PRIMARY_WINDOW_BP,
            "hosts_total": discovery["hosts_total"],
            "hosts_with_more_than_one_nested_event": discovery[
                "hosts_with_more_than_one_nested_event"
            ],
            "pairs_within_window": discovery["pairs_within_window"],
        },
        "verdict_summary": verdict_summary,
        "positional_chance": chance,
        "discordance_chance": discord,
        "ibd_affirmation_chance": affirmation,
        "sweep": sweep,
        "any_candidate_recurrent": len(candidates),
        "candidate_evidence_lines": candidates,
        "method_notes": [
            "identity by descent requires same host copy, breakpoints within "
            f"+/-{PRIMARY_WINDOW_BP} bp, same child subfamily and TSD length "
            f"difference <= {TSD_LENGTH_TOLERANCE_BP} bp",
            "candidate independent recurrence requires same host copy and "
            "position, plus a subfamily, TSD-length or TSD-core discordance",
            "a call reporting no TSD makes the pair unevaluable, never IBD and "
            "never discordant; unevaluable pairs are counted separately",
            f"TSD cores shorter than {TSD_CORE_MIN_BP} bp after stripping the "
            "shared terminal homopolymer are unevaluable rather than concordant, "
            "a threshold fixed before the data was read",
            f"the +/-{PRIMARY_WINDOW_BP} bp window is primary and the "
            f"+/-{SWEEP_WINDOWS_BP[0]} and +/-{SWEEP_WINDOWS_BP[1]} bp sweeps are "
            "reported in full regardless of outcome, per the frozen plan",
            "the positional chance rate is (2w+1)/L per host; the discordance "
            "chance rate comes from two independent draws of the observed "
            "(subfamily, TSD length) distribution",
            "no pair is adjudicated in code; candidates are written out for "
            "manual review",
            "GT and GQ are never read: the callset reader parses CHROM/POS/ID/INFO "
            "only and never touches the FORMAT column",
            "counting is restricted to chr1-22 and chrX",
        ],
        "deviations": [
            "the per-call table named as an input does not exist; per-call TSD and "
            "child subfamily are recovered by joining the dedup output back to the "
            "callsets by sample name within the frozen +/-10 bp rule, and the join "
            "is gated on reproducing dedup's own carrier counts",
            "denominators here are 5 genomes and are not comparable to published "
            "pooled-SVAN host counts; they are labelled rather than aligned",
        ],
    }
    (args.outdir / "recurrence.json").write_text(
        json.dumps(report, indent=2, sort_keys=True, default=str) + "\n"
    )
    write_pairs_csv(args.outdir / "recurrence_pairs.csv", all_rows)
    write_report(args.outdir / "recurrence.md", report)
    print(
        f"candidate-recurrent pairs at +/-{PRIMARY_WINDOW_BP} bp: {len(candidates)}; "
        f"wrote recurrence.md / .json and recurrence_pairs.csv",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
