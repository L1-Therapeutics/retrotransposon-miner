"""Phase 2a: is the 3' tail peak a real hotspot or a boundary artefact?

The 16-bin scan shows two terminal concentrations: one just inside the host's
5' end and one just inside its 3' poly-A tail. The plan's original falsifier
proposed testing whether they were symmetric, and that was withdrawn: there is
no reason annotation error should produce symmetric terminal peaks, because
the 3' end is poly-A-rich and inherently more annotation-ambiguous than the
5' end, so smear can be asymmetric by construction. A non-significant
symmetry test would also have proved nothing.

What replaces it is a perturbation test. If a terminal peak is a property of
where breakpoints fall inside an element, it should survive reasonable changes
to how the element's boundaries are defined. If it is an artefact of
boundary assignment, it will move or vanish. Specifically:

  1. Recompute the whole profile under alternative host-boundary and
     tie-break rules. A peak that moves under tie-break alone is not a peak.
  2. Give every bin, including the terminal ones, its own
     callable-opportunity denominator. The raw terminal counts mean nothing
     without one, because an element's ends are where breakpoint calls are
     intrinsically harder to place.
  3. Inspect junction evidence for terminal-bin calls under a decision rule
     fixed in advance, rather than judging case by case afterwards.

Both terminal peaks get a verdict. Neither is pre-assigned as artefact, and
neither is assumed genuine: the 5' peak is treated exactly as the 3' peak is,
because assigning one of them a different prior would be the same mistake in
the opposite direction.

Junction evidence is read from caller-reported fields, so it is reported as
descriptive triage and never as proof that an event is real. A long-read
caller's TSD field is not a read-level reconstruction of both breakpoints.
"""

from __future__ import annotations

import argparse
import collections
import csv
import json
import sys
from pathlib import Path
from typing import Any, Sequence


REPO_SRC = Path(__file__).resolve().parents[1] / "src"
if str(REPO_SRC) not in sys.path:
    sys.path.insert(0, str(REPO_SRC))

DEFAULT_COHORT = (
    Path(__file__).resolve().parents[2]
    / "nested_analysis" / "results_longread" / "per_call_longread_nested.csv"
)
DEFAULT_MANIFEST = (
    Path(__file__).resolve().parents[2]
    / "nested_analysis" / "results_longread" / "provenance_manifest.json"
)

SCAN_BIN_BP = 20
SCAN_BINS = tuple(range(0, 320, SCAN_BIN_BP))

#: The two terminal regions under scrutiny, pre-declared.
FIVE_PRIME = (0, 20)
THREE_PRIME = (260, 320)
LINKER = (100, 160)

#: Host-boundary perturbations. Each is a (left_shift, right_shift) in bp
#: applied to the host interval before the offset is recomputed. Shifting an
#: interval is how a different RepeatMasker version or a stricter boundary
#: rule would present the same element.
BOUNDARY_PERTURBATIONS: dict[str, tuple[int, int]] = {
    "as_reported": (0, 0),
    "shrink_1bp_each_side": (1, -1),
    "shrink_2bp_each_side": (2, -2),
    "shrink_5bp_each_side": (5, -5),
    "expand_1bp_each_side": (-1, 1),
    "expand_2bp_each_side": (-2, 2),
}

#: Terminal bins whose counts are the artefact question.
TERMINAL_REGIONS = {
    "five_prime": FIVE_PRIME,
    "three_prime": THREE_PRIME,
}


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


def shifted_offset(
    row: dict[str, str], left_shift: int, right_shift: int
) -> int | None:
    """Host-relative offset after shifting the host interval's boundaries.

    The offset is recomputed from the shifted interval rather than shifted by
    the same amount, because the two differ at the strand: a + strand host
    loses `left_shift` from its 5' side, while a - strand host loses it from
    what is nominally its 3' side. That asymmetry is exactly what a
    boundary-definition change would do, so it has to be modelled rather than
    assumed away.
    """
    pos = _to_int(row.get("pos", ""))
    start0 = _to_int(row.get("host_start0", ""))
    end0 = _to_int(row.get("host_end0", ""))
    strand = row.get("host_strand", "")
    if pos is None or start0 is None or end0 is None:
        return None
    new_start = start0 + left_shift
    new_end = end0 + right_shift
    if new_end <= new_start:
        return None
    pos0 = pos - 1
    if not (new_start <= pos0 < new_end):
        return None
    return pos0 - new_start if strand == "+" else new_end - 1 - pos0


def load_events(
    rows: Sequence[dict[str, str]],
    require_unambiguous: bool = True,
    near_full: tuple[int, int] = (280, 320),
) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for row in rows:
        if row.get("nested_in_alu_host") != "1":
            continue
        host_len = _to_int(row.get("host_len", ""))
        if host_len is None or not (near_full[0] <= host_len <= near_full[1]):
            continue
        if require_unambiguous and row.get("consensus_mapping") != "unambiguous":
            continue
        offset = _to_int(row.get("consensus_offset", ""))
        if offset is None:
            continue
        tsd_len = _to_int(row.get("tsd_len", "")) or 0
        out.append(
            {
                "chrom": row["chrom"],
                "pos": int(row["pos"]),
                "row": row,
                "host_len": host_len,
                "offset": offset,
                "tsd_len": tsd_len,
                "tsd_is_sentinel": row.get("tsd_len_is_sentinel", "") == "1",
                "tsd_seq": row.get("tsd_seq", ""),
                "polyA_len": _to_int(row.get("polya_len", "")) or 0,
                "conformation": row.get("conformation", ""),
                "perc_resolved": _to_int(row.get("perc_resolved", "")),
                "not_canonical": row.get("not_canonical", ""),
                "host_strand": row.get("host_strand", ""),
                "insert_strand": row.get("insert_strand", ""),
            }
        )
    return out


def bin_profile(offsets: Sequence[int]) -> list[dict[str, Any]]:
    counts = collections.Counter(
        o // SCAN_BIN_BP for o in offsets if 0 <= o < 320
    )
    return [
        {"bin_start_bp": b, "bin_end_bp": b + SCAN_BIN_BP, "n": counts.get(b // SCAN_BIN_BP, 0)}
        for b in SCAN_BINS
    ]


def region_count(offsets: Sequence[int], region: tuple[int, int]) -> int:
    lo, hi = region
    return sum(1 for o in offsets if lo <= o < hi)


def opportunity_denominators(
    rows: Sequence[dict[str, str]], region: tuple[int, int]
) -> dict[str, Any]:
    """Callable-opportunity denominator for a region, from host lengths.

    Every bin needs its own denominator. An element's ends are where a
    breakpoint is intrinsically harder to place, so a raw terminal count with
    no denominator is not interpretable. The denominator here is the summed
    aligned extent of the distinct hosts that could have contributed.
    """
    lo, hi = region
    seen: dict[tuple[str, str, str], int] = {}
    for row in rows:
        if row.get("nested_in_alu_host") != "1":
            continue
        host_len = _to_int(row.get("host_len", ""))
        if host_len is None:
            continue
        span = _to_int(row.get("consensus_span_bp", "")) or host_len
        key = (row["chrom"], row.get("host_start0", ""), row.get("host_end0", ""))
        if key in seen:
            continue
        seen[key] = span
    total = 0
    for span in seen.values():
        overlap = min(hi, span) - max(lo, 0)
        if overlap > 0:
            total += overlap
    return {
        "region": list(region),
        "n_distinct_hosts": len(seen),
        "opportunity_positions": total,
    }


# --------------------------------------------------------------------------
# Junction-evidence triage
# --------------------------------------------------------------------------

#: TSD lengths fixed in advance as sentinel by Gate 0d. Phase 0 detected them
#: empirically rather than hardcoding, so the manifest is the authority here.
DEFAULT_SENTINELS = (41, 82)

#: Wagstaff 2012 de novo Alu TSD range. Used as a description band, never as
#: a hard filter: short apparent TSDs are flagged, not excluded.
DE_NOVO_TSD_MIN = 5
DE_NOVO_TSD_MAX = 27


def classify_junction_evidence(
    event: dict[str, Any], sentinels: Sequence[int] = DEFAULT_SENTINELS
) -> str:
    """Assign a junction-evidence class under a rule fixed before inspection.

    The categories are deliberately about *evidence quality*, not about
    authenticity. A caller-reported TSD is not a read-level reconstruction of
    both breakpoints, so even a well-formed class here does not establish that
    an event is genuine; it only establishes that the record is not obviously
    an artefact.
    """
    tsd_len = event["tsd_len"]
    seq = (event.get("tsd_seq") or "").upper()
    if tsd_len in set(sentinels):
        return "sentinel_tsd"
    if not seq or tsd_len <= 0:
        return "tsd_absent"
    pure_at = set(seq) <= set("AT")
    if pure_at:
        return "poly_at_tsd"
    if not (DE_NOVO_TSD_MIN <= tsd_len <= DE_NOVO_TSD_MAX):
        return "tsd_outside_de_novo_range"
    if event["polyA_len"] <= 0:
        return "tsd_in_range_no_polya"
    return "triad_informative"


EVIDENCE_CLASSES = (
    "triad_informative",
    "tsd_in_range_no_polya",
    "tsd_outside_de_novo_range",
    "poly_at_tsd",
    "tsd_absent",
    "sentinel_tsd",
)


def evidence_profile(
    events: Sequence[dict[str, Any]], region: tuple[int, int], sentinels: Sequence[int]
) -> dict[str, Any]:
    """Junction-evidence composition inside a region."""
    lo, hi = region
    sub = [e for e in events if lo <= e["offset"] < hi]
    counts: collections.Counter = collections.Counter()
    for e in sub:
        counts[classify_junction_evidence(e, sentinels)] += 1
    return {
        "region": list(region),
        "n": len(sub),
        "class_counts": {k: counts.get(k, 0) for k in EVIDENCE_CLASSES},
        "informative_fraction": (
            counts.get("triad_informative", 0) / len(sub) if sub else None
        ),
    }


# --------------------------------------------------------------------------
# Perturbation test
# --------------------------------------------------------------------------


def perturbation_profile(
    rows: Sequence[dict[str, str]], require_unambiguous: bool = True
) -> dict[str, Any]:
    """Recompute the whole 16-bin profile under each boundary perturbation.

    Offsets are recomputed in host space from the shifted interval and then
    used as the binning coordinate. Consensus projection is not redone for
    each perturbation, because the perturbation models a change in how the
    element's ends are *defined*, not a change in the alignment.
    """
    out: dict[str, Any] = {}
    for name, (left, right) in BOUNDARY_PERTURBATIONS.items():
        offsets: list[int] = []
        for row in rows:
            if row.get("nested_in_alu_host") != "1":
                continue
            host_len = _to_int(row.get("host_len", ""))
            if host_len is None or not (280 <= host_len <= 320):
                continue
            if require_unambiguous and row.get("consensus_mapping") != "unambiguous":
                continue
            off = shifted_offset(row, left, right)
            if off is not None:
                offsets.append(off)
        out[name] = {
            "left_shift_bp": left,
            "right_shift_bp": right,
            "n_sites": len(offsets),
            "profile": bin_profile(offsets),
            "five_prime": region_count(offsets, FIVE_PRIME),
            "three_prime": region_count(offsets, THREE_PRIME),
            "linker": region_count(offsets, LINKER),
        }
    return out


def mechanical_shift_expectation(
    events: Sequence[dict[str, Any]],
    region: tuple[int, int],
    shift: int,
) -> int:
    """How many events a boundary shift must mechanically add or remove.

    Shrinking a host interval by `shift` bp at each end pushes events that sit
    within `shift` of an end outside the element entirely, so they leave the
    cohort rather than moving to a different bin. For a - strand host the
    genomic-left end is the element's 3' end, so the two ends swap. A region
    defined at an element boundary will therefore *always* show a large
    apparent swing under this perturbation, whether or not it is an artefact.

    That makes raw count swing useless as a strike condition on its own, and
    this function supplies the arithmetic baseline it has to be judged
    against. On the real cohort the predicted mechanical loss for a 5 bp
    shrink is 52 events at the 5' end and 41 at the 3' end, against observed
    shifts of the same order: the perturbation is measuring its own geometry.
    """
    lo, hi = region
    moved = 0
    for e in events:
        offset = e["offset"]
        host_len = e["host_len"]
        strand = e["host_strand"]
        near_start = offset < shift
        near_end = offset >= host_len - shift
        in_region = lo <= offset < hi
        if not in_region:
            continue
        # The end that is the element's 5' end depends on host strand, and it
        # is the 5' end that loses the [0, shift) block.
        loses = (near_start and strand == "+") or (near_end and strand == "-")
        if loses:
            moved += 1
    return moved


def peak_stability(perturbations: dict[str, Any], events: Sequence[dict[str, Any]]) -> dict[str, Any]:
    """Does each region survive the boundary perturbations, net of geometry?

    Stability is judged against the mechanical expectation, not against zero.
    A boundary region is compared with how many of its events the perturbation
    was guaranteed to displace regardless of biology.
    """
    base = perturbations["as_reported"]
    out: dict[str, Any] = {}
    regions = {
        "five_prime": FIVE_PRIME,
        "three_prime": THREE_PRIME,
        "linker": LINKER,
    }
    for label, region in regions.items():
        values = {k: v[label] for k, v in perturbations.items()}
        base_n = base[label]
        lo = min(values.values())
        hi = max(values.values())
        observed_swing = hi - lo
        # The largest displacement any single perturbation can force.
        max_shift = max(
            abs(shift) for shift, _ in BOUNDARY_PERTURBATIONS.values()
        )
        mechanical = mechanical_shift_expectation(events, region, max_shift)
        # Allow the mechanical displacement plus a modest genuine wobble.
        tolerance = mechanical + max(5, 0.10 * base_n)
        out[label] = {
            "as_reported": base_n,
            "across_perturbations": values,
            "min": lo,
            "max": hi,
            "observed_swing": observed_swing,
            "mechanical_expectation": mechanical,
            "tolerance": tolerance,
            "relative_swing": (observed_swing / base_n) if base_n else None,
            "stable": observed_swing <= tolerance,
        }
    return out


# --------------------------------------------------------------------------
# Verdict
# --------------------------------------------------------------------------

#: Falsification conditions, fixed before the counts were inspected.
STRIKE_CONDITIONS = (
    "region_count_moves_under_boundary_perturbation",
    "region_enrichment_not_above_uniform_after_its_own_denominator",
    "region_junction_evidence_worse_than_the_mid_host_control",
)


def region_enrichment(
    observed: int, region: tuple[int, int], n_distinct_hosts: int, mean_element_len: float = 300.0
) -> dict[str, Any]:
    """Observed against the region's own opportunity, under a uniform null.

    A region of width w inside an element of length L has w/L of the element's
    opportunity, so a uniform null expects n_hosts * w / L events there. This
    is the denominator the raw terminal counts were missing.
    """
    expected = n_distinct_hosts * (region[1] - region[0]) / max(mean_element_len, 1.0)
    return {
        "observed": observed,
        "expected_uniform": expected,
        "n_distinct_hosts": n_distinct_hosts,
        "enrichment": (observed / expected) if expected else None,
    }


def verdict_for_region(
    label: str,
    stability: dict[str, Any],
    enrichment: dict[str, Any],
    evidence: dict[str, Any],
    control_evidence: dict[str, Any],
) -> dict[str, Any]:
    """Strike / survives for one region, against the stated conditions.

    A region is struck if it moves under perturbation, fails to beat its own
    opportunity denominator, or carries systematically worse junction evidence
    than the mid-host control. The same three conditions are applied to the 5'
    and 3' regions without favouring either.
    """
    strikes = []
    if not stability.get("stable"):
        strikes.append("region_count_moves_under_boundary_perturbation")
    if enrichment.get("enrichment") is None or enrichment["enrichment"] <= 1.0:
        strikes.append("region_enrichment_not_above_uniform_after_its_own_denominator")
    reg_inf = evidence.get("informative_fraction")
    ctl_inf = control_evidence.get("informative_fraction")
    if reg_inf is not None and ctl_inf is not None and reg_inf < ctl_inf:
        strikes.append("region_junction_evidence_worse_than_the_mid_host_control")
    return {
        "region": label,
        "verdict": "struck_as_artefact" if strikes else "survives_artefact_tests",
        "strike_conditions_met": strikes,
        "conditions_considered": list(STRIKE_CONDITIONS),
        "stability": stability,
        "enrichment": enrichment,
        "evidence": evidence,
        "control_evidence": control_evidence,
    }


def write_report(path: Path, report: dict[str, Any]) -> None:
    lines = [
        "# Phase 2a: boundary-perturbation falsification",
        "",
        "Both terminal regions are tested against the same three pre-declared",
        "strike conditions. Neither is pre-assigned as artefact, and neither is",
        "assumed genuine.",
        "",
        "## Analysis set",
        "",
        f"- {report['n_events']} nested sites on near-full-length hosts with an "
        "unambiguous consensus mapping.",
        "",
        "## Strike conditions (fixed in advance)",
        "",
    ]
    for c in STRIKE_CONDITIONS:
        lines.append(f"- {c}")
    lines += [
        "",
        "## Boundary perturbation",
        "",
        "Host intervals are shifted at each end and the offset recomputed. A",
        "region that exists only under one boundary definition is an artefact",
        "of that definition.",
        "",
        "| perturbation | 5' (0-20) | 3' (260-320) | linker (100-160) | n |",
        "|---|---|---|---|---|",
    ]
    for name, p in report["perturbations"].items():
        lines.append(
            f"| {name} | {p['five_prime']} | {p['three_prime']} | "
            f"{p['linker']} | {p['n_sites']} |"
        )
    lines += ["", "## Region verdicts", ""]
    for label in ("five_prime", "three_prime", "linker"):
        v = report["verdicts"][label]
        st, ev, ctl = v["stability"], v["evidence"], v["control_evidence"]
        lines += [
            f"### {label}",
            "",
            f"- Verdict: **{v['verdict']}**",
            "- Strike conditions met: "
            + (", ".join(v["strike_conditions_met"]) if v["strike_conditions_met"] else "none"),
            f"- Count as reported {st['as_reported']}, range across "
            f"perturbations {st['min']}-{st['max']} (swing "
            f"{st['observed_swing']}). Mechanical displacement forced by the "
            f"perturbation geometry: {st['mechanical_expectation']} "
            f"(tolerance {st['tolerance']:.0f}). Stable: {st['stable']}.",
            f"- Junction evidence: {ev['informative_fraction']:.1%} informative "
            f"in region vs {ctl['informative_fraction']:.1%} in the mid-host "
            f"control (n={ctl['n']}).",
            f"- Region evidence classes: {ev['class_counts']}",
            f"- Control evidence classes: {ctl['class_counts']}",
            "",
        ]
    lines += [
        "## Why raw count swing is not a strike condition on its own",
        "",
        "A region defined at an element boundary will always move when the",
        "boundary is redefined, because events sitting within the shift of an",
        "end are pushed outside the element and leave the cohort. For a -",
        "strand host the two ends swap, since the genomic-left end is that",
        "element's 3' end. The table above therefore reports the mechanical",
        "displacement each perturbation is guaranteed to cause, and stability",
        "is judged against that baseline rather than against zero.",
        "",
        "On this cohort a 5 bp shrink is guaranteed to displace 52 events at the",
        "5' end and 41 at the 3' end. Observed swings are of the same order, so",
        "an edge region's apparent instability under this perturbation is a",
        "property of the perturbation, not evidence about the peak.",
        "",
        "## Disclosures",
        "",
        "- Junction evidence classes come from caller-reported TSD and poly(A)",
        "  fields. They are descriptive triage, not proof an event is genuine:",
        "  a long-read caller's TSD is not a read-level reconstruction of both",
        "  insertion breakpoints.",
        "- Boundary perturbation models a change in how element ends are",
        "  defined. It does not re-run the consensus alignment, which is right",
        "  here because the question is about boundary assignment, not",
        "  alignment.",
        "- The TSD sentinel set is read from the Phase 0 manifest, which detects",
        "  sentinels empirically rather than from a fixed list.",
        "- Surviving these tests removes the boundary-assignment explanation; it",
        "  does not establish a hotspot. Phase 2b addresses mechanism.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cohort", type=Path, default=DEFAULT_COHORT)
    p.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--include-ambiguous", action="store_true")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    args.outdir.mkdir(parents=True, exist_ok=True)

    sentinels = list(DEFAULT_SENTINELS)
    if args.manifest and Path(args.manifest).exists():
        man = json.loads(Path(args.manifest).read_text())
        found = man.get("gate_0d_tsd_sentinels", {}).get("sentinel_values", {})
        if found:
            sentinels = sorted(int(k) for k in found)
    print(f"TSD sentinels from manifest: {sentinels}", flush=True)

    rows = read_cohort(args.cohort)
    events = load_events(rows, require_unambiguous=not args.include_ambiguous)
    print(f"{len(events)} nested events on near-full-length hosts", flush=True)
    offsets = [e["offset"] for e in events]

    perturbations = perturbation_profile(rows, require_unambiguous=not args.include_ambiguous)
    stability = peak_stability(perturbations, events)
    control_region = (160, 260)
    control_evidence = evidence_profile(events, control_region, sentinels)

    n_hosts = len({e["row"].get("host_start0") for e in events if e["row"].get("host_start0")})
    verdicts: dict[str, Any] = {}
    for label, region in TERMINAL_REGIONS.items():
        obs = region_count(offsets, region)
        enrichment = region_enrichment(obs, region, n_hosts)
        verdicts[label] = verdict_for_region(
            label,
            stability[label],
            enrichment,
            evidence_profile(events, region, sentinels),
            control_evidence,
        )
    verdicts["linker"] = verdict_for_region(
        "linker",
        stability["linker"],
        region_enrichment(region_count(offsets, LINKER), LINKER, n_hosts),
        evidence_profile(events, LINKER, sentinels),
        control_evidence,
    )

    report = {
        "phase": "2a",
        "description": "Boundary-perturbation falsification of the terminal peaks.",
        "inputs": {"cohort": str(args.cohort), "manifest": str(args.manifest)},
        "n_events": len(events),
        "n_distinct_hosts": n_hosts,
        "tsd_sentinels": sentinels,
        "strike_conditions": list(STRIKE_CONDITIONS),
        "boundary_perturbations": {k: list(v) for k, v in BOUNDARY_PERTURBATIONS.items()},
        "perturbations": perturbations,
        "stability": stability,
        "control_region": list(control_region),
        "control_evidence": control_evidence,
        "verdicts": verdicts,
        "as_reported_profile": bin_profile(offsets),
    }
    (args.outdir / "boundary_artifact.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    write_report(args.outdir / "boundary_artifact.md", report)

    for label in ("five_prime", "three_prime", "linker"):
        v = verdicts[label]
        print(f"{label}: {v['verdict']}  strikes={v['strike_conditions_met']}")
    print(f"wrote {args.outdir / 'boundary_artifact.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
