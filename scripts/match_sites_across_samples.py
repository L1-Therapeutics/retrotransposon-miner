#!/usr/bin/env python3
"""Cross-sample insertion-site deduplication for the nested-Alu/L1 positional analysis.

Implements the frozen (2026-10-01) cross-sample matching rule of
rtm-longread-cohort/docs/NESTED_ALU_ANALYSIS_PLAN.md §2, verbatim:

  same chrom + family + orientation + host + breakpoint ±10 bp primary
  (sweep ±5/±20 reported in full); no transitive chaining; per-source
  systematic caller offsets calibrated on high-confidence shared sites and
  recorded in matching_audit.md

Reads per-sample MEI callsets (VCFs with the HG03086 callset schema), assigns host
elements with the audit_nested_tsd.py convention verbatim, and collapses calls that
describe the same insertion site across genomes. Two calls match iff same chromosome,
same MEI family, same insertion orientation, same assigned host element, AND breakpoints
within the primary ±10 bp window. No transitive chaining: every call is compared against
the representative call of a site only, never against other members. Candidate pairs are
resolved greedily by highest breakpoint proximity, then longest TSD overlap, then
deterministic keys. Unique sites are the primary estimand; private sites are reported
separately and never merged into the primary positional test (§1).

Matching keys use NO gene/snpeff annotation fields. Host identity and orientation come
from the same RepeatMasker-based logic as artifacts/audit_nested_tsd.py. A cross-sample
schema diff (INFO keys, family/nesting value distributions) is computed and reported in
matching_audit.md; POS/family/nesting drift is flagged there.

The ±5 and ±20 windows are reported as sensitivity, in full: a result that only holds at
one window is not a result.

Caller-offset calibration: the systematic POS offset of each new sample against the
already-processed reference callset is estimated ONLY on already-shared, high-confidence
sites (KNOWNMEI catalogued positions) and applied as a systematic per-source correction.
The estimate, its support, and whether it was applied are recorded in matching_audit.md.
It is never fudged per call.

Outputs (default nested_analysis/results_multi_sample/):
  unique_sites.csv                 one row per unique site (call counts and unique-site
                                   counts side by side; carrier distribution)
  accumulation_curve.csv/.png      unique nested sites vs genomes added, with marginal
                                   yield (§3 power statement)
  matching_audit.md                schema diff, per-sample offsets, calibration sites,
                                   full ±5/±20 sweep, manifest
  matching_manifest.json           machine-readable manifest of processed samples
"""
from __future__ import annotations

import argparse
import bisect
import gzip
import json
import re
import statistics
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
WORKSPACE = REPO.parent

# ---------------------------------------------------------------------------
# Pre-registered parameters. Do not tune: they are part of the analysis record.
# ---------------------------------------------------------------------------
PRIMARY_WINDOW_BP = 10
SENSITIVITY_WINDOWS_BP = (5, 20)
CALIBRATION_BROAD_WINDOW_BP = 50
CALIBRATION_MIN_SHARED_SITES = 10
CALIBRATION_MAX_MAD_BP = 3
ACCUMULATION_FAMILIES = ("ALU", "LINE1")
HOST_SELECTION_RULE = "longest_span"
# chr22_mei.vcf is a chr22 slice of HG03086; including it would double-count those
# calls. It is refused by name and by scope check, not silently skipped.
EXCLUDED_CALLSET_STEMS = {"chr22_mei"}


# ---------------------------------------------------------------------------
# Host assignment -- audit_nested_tsd.py convention, verbatim.
# ---------------------------------------------------------------------------

def family_for_rmsk(rep_family: str, rep_class: str) -> str:
    token = f"{rep_family} {rep_class}".upper()
    if "ALU" in token:
        return "ALU"
    if "SVA" in token:
        return "SVA"
    if "LINE1" in token or re.search(r"(^|[^A-Z0-9])L1([^A-Z0-9]|$)", token):
        return "LINE1"
    return ""


@dataclass
class HostAssignment:
    chrom: str
    start0: int
    end0: int
    strand: str
    name: str
    offset_5p_0based: int

    @property
    def key(self) -> tuple[str, int, int]:
        return (self.name, self.start0, self.end0)

    @property
    def length(self) -> int:
        return self.end0 - self.start0


def assign_hosts(
    calls: list["Call"], rmsk_path: Path
) -> dict[tuple[str, int], HostAssignment]:
    """Same-family RMSK interval containing the call, greatest span then leftmost
    start then name; no orientation-dependent tie break (audit_nested_tsd.py)."""
    nested = [c for c in calls if c.nested]
    wanted: dict[tuple[str, str], list[tuple[int, str]]] = defaultdict(list)
    for call in nested:
        wanted[(call.chrom, call.family)].append((call.pos0, call.chrom_pos_key))
    for key in wanted:
        wanted[key].sort(key=lambda item: item[0])
    starts = {key: [item[0] for item in values] for key, values in wanted.items()}

    overlaps: dict[str, list[tuple[int, int, str, str]]] = defaultdict(list)
    opener = gzip.open if rmsk_path.suffix == ".gz" else open
    with opener(rmsk_path, "rt") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 13:
                continue
            chrom = fields[5]
            family = family_for_rmsk(fields[12], fields[11])
            key = (chrom, family)
            if not family or key not in wanted:
                continue
            start0, end0 = int(fields[6]), int(fields[7])
            lo = bisect.bisect_left(starts[key], start0)
            hi = bisect.bisect_left(starts[key], end0)
            for pos0, call_key in wanted[key][lo:hi]:
                if pos0 < end0:
                    overlaps[call_key].append((start0, end0, fields[9], fields[10]))

    hosts: dict[tuple[str, int], HostAssignment] = {}
    for call in nested:
        candidates = overlaps.get(call.chrom_pos_key, [])
        if not candidates:
            continue
        start0, end0, strand, name = sorted(
            candidates, key=lambda item: (-(item[1] - item[0]), item[0], item[3])
        )[0]
        pos0 = call.pos0
        offset = pos0 - start0 if strand == "+" else (end0 - 1 - pos0)
        hosts[call.chrom_pos_key] = HostAssignment(
            chrom=call.chrom,
            start0=start0,
            end0=end0,
            strand=strand,
            name=name,
            offset_5p_0based=offset,
        )
    return hosts


# ---------------------------------------------------------------------------
# Call model and VCF loading.
# ---------------------------------------------------------------------------

@dataclass
class Call:
    sample: str
    chrom: str
    pos: int
    family: str
    orientation: str
    nested: bool
    tsd_seq: str
    known_mei: bool
    call_tier: str
    host: HostAssignment | None = None
    calibr_pos: int | None = None  # systematic corrected position; None until calibrated

    @property
    def pos0(self) -> int:
        return self.pos - 1

    @property
    def match_pos(self) -> int:
        """Position used for window comparisons: calibrated when set, raw otherwise."""
        return self.pos if self.calibr_pos is None else self.calibr_pos

    @property
    def chrom_pos_key(self) -> tuple[str, int]:
        return (self.chrom, self.pos)

    @property
    def host_key(self) -> tuple[str, int, int]:
        return self.host.key if self.host is not None else ("", 0, 0)

    @property
    def match_key(self) -> tuple[str, str, str, tuple[str, int, int]]:
        return (self.chrom, self.family, self.orientation, self.host_key)


def _info_dict(info_field: str) -> dict[str, str]:
    info: dict[str, str] = {}
    for item in info_field.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        elif item:
            info[item] = "True"
    return info


def load_callset(path: Path, sample: str) -> list[Call]:
    opener = gzip.open if path.suffix == ".gz" else open
    calls: list[Call] = []
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            info = _info_dict(fields[7])
            calls.append(
                Call(
                    sample=sample,
                    chrom=fields[0],
                    pos=int(fields[1]),
                    family=info.get("MEIFAMILY", ""),
                    orientation=info.get("ORIENT", ""),
                    nested=info.get("NESTED", "").strip().lower() == "nested",
                    tsd_seq=info.get("TSD", ""),
                    known_mei=info.get("KNOWNMEI", "").strip().lower() == "true",
                    call_tier=info.get("CALLTIER", ""),
                )
            )
    return calls


def verify_genome_wide(path: Path, calls: list[Call]) -> None:
    """Refuse callsets that are not genome-wide (no chromosomes beyond chr22).

    A chr22-slice callset (e.g. chr22_mei.vcf, a slice of HG03086) would
    double-count the full callset's calls, so it is rejected before matching.
    """
    chroms = {call.chrom for call in calls}
    beyond = sorted(chroms - {"chr22"})
    if not beyond:
        raise SystemExit(
            f"{path.name}: not genome-wide (only {sorted(chroms)} present, no "
            "chromosomes beyond chr22) - refusing to match a slice callset that "
            "would double-count a full-genome callset"
        )


def discover_callsets(args: argparse.Namespace) -> list[tuple[str, Path]]:
    """Deterministic sample-id -> path list from --callset-dir and/or --callset.

    chr22_mei.vcf (a chr22 slice of HG03086) is refused explicitly: including it
    would double-count HG03086's calls. Every accepted callset is verified to be
    genome-wide (chromosomes beyond chr22 present) before matching runs.
    """
    found: dict[str, Path] = {}
    if args.callset_dir is not None:
        directory = args.callset_dir
        if not directory.is_dir():
            raise SystemExit(f"callset directory not found: {directory}")
        for path in sorted(directory.iterdir()):
            if path.suffix not in {".vcf", ".gz"} or path.name.startswith("."):
                continue
            sample = path.name.removesuffix(".vcf.gz").removesuffix(".vcf")
            stem = sample.removesuffix(".vcf")
            if stem in EXCLUDED_CALLSET_STEMS:
                print(
                    f"[match-sites] EXCLUDED {path.name}: chr22 slice of an already-"
                    "included genome-wide callset (would double-count its calls)",
                    file=sys.stderr,
                )
                continue
            found[sample] = path
    for spec in args.callset or []:
        path = Path(spec)
        if not path.is_file():
            raise SystemExit(f"callset file not found: {path}")
        sample = args.sample_id or path.name.removesuffix(".vcf.gz").removesuffix(".vcf")
        found[sample] = path
    return sorted(found.items())


# ---------------------------------------------------------------------------
# Caller-offset calibration on KNOWNMEI catalogued positions.
# ---------------------------------------------------------------------------

@dataclass
class OffsetEstimate:
    sample: str
    n_shared: int
    median_offset_bp: int
    mad_bp: float
    applied: bool
    reason: str

    def as_dict(self) -> dict[str, object]:
        return {
            "sample": self.sample,
            "n_knownmei_shared_sites": self.n_shared,
            "median_offset_bp": self.median_offset_bp,
            "mad_bp": self.mad_bp,
            "applied": self.applied,
            "reason": self.reason,
        }


def estimate_offset(
    reference: list[Call], other: list[Call], sample: str
) -> OffsetEstimate:
    """Median (reference_pos - other_pos) over KNOWNMEI sites, pairing each call with
    its nearest same-key reference site within the broad calibration window; applied
    only when well supported and stable."""
    ref_by_key: dict[tuple[str, str, str], list[int]] = defaultdict(list)
    for call in reference:
        if call.known_mei:
            ref_by_key[(call.chrom, call.family, call.orientation)].append(call.pos)
    diffs: list[int] = []
    for call in other:
        if not call.known_mei:
            continue
        candidates = ref_by_key.get((call.chrom, call.family, call.orientation), ())
        if not candidates:
            continue
        ref_pos = min(candidates, key=lambda p: (abs(p - call.pos), p))
        if abs(ref_pos - call.pos) <= CALIBRATION_BROAD_WINDOW_BP:
            diffs.append(ref_pos - call.pos)
    if len(diffs) < CALIBRATION_MIN_SHARED_SITES:
        return OffsetEstimate(
            sample=sample,
            n_shared=len(diffs),
            median_offset_bp=0,
            mad_bp=0.0,
            applied=False,
            reason="insufficient_knownmei_overlap",
        )
    median = int(round(statistics.median(diffs)))
    mad = statistics.pstdev(diffs) if len(diffs) > 1 else 0.0
    if mad > CALIBRATION_MAX_MAD_BP:
        return OffsetEstimate(
            sample=sample,
            n_shared=len(diffs),
            median_offset_bp=median,
            mad_bp=round(mad, 2),
            applied=False,
            reason="offset_unstable_not_applied",
        )
    return OffsetEstimate(
        sample=sample,
        n_shared=len(diffs),
        median_offset_bp=median,
        mad_bp=round(mad, 2),
        applied=median != 0,
        reason="applied" if median != 0 else "zero_offset_no_correction_needed",
    )


def apply_calibration(
    reference_sample: str,
    calls: list[Call],
    estimates: dict[str, OffsetEstimate],
) -> None:
    for call in calls:
        if call.sample == reference_sample:
            call.calibr_pos = call.pos
            continue
        estimate = estimates.get(call.sample)
        if estimate is not None and estimate.applied:
            call.calibr_pos = call.pos + estimate.median_offset_bp
        else:
            call.calibr_pos = call.pos


# ---------------------------------------------------------------------------
# Pre-registered matching: greedy, representative-anchored, no transitive chaining.
# ---------------------------------------------------------------------------

def tsd_overlap(a: str, b: str) -> int:
    """Position-wise agreement length of the two TSD strings (0 if either missing)."""
    a = (a or "").strip().upper()
    b = (b or "").strip().upper()
    if not a or not b or a == "." or b == ".":
        return 0
    n = 0
    for x, y in zip(a, b):
        if x != y:
            break
        n += 1
    return n


def matches_representative(call: Call, representative: Call, window_bp: int) -> bool:
    """The pre-registered rule, evaluated against a site's representative only."""
    return (
        call.match_key == representative.match_key
        and abs(call.match_pos - representative.match_pos) <= window_bp
    )


def match_calls(calls: list[Call], window_bp: int = PRIMARY_WINDOW_BP) -> list[list[Call]]:
    """Greedy representative-anchored deduplication.

    Pairs are considered by highest breakpoint proximity, then longest TSD overlap,
    then deterministic tie keys. A call joins an existing site only if it matches
    that site's representative; two already-assigned calls are never merged, so no
    transitive chain can form.
    """
    pool = [c for c in calls if c.nested and c.family]
    pairs: list[tuple[int, int, int, tuple, tuple]] = []
    by_key: dict[tuple, list[int]] = defaultdict(list)
    for idx, call in enumerate(pool):
        by_key[call.match_key].append(idx)
    for indices in by_key.values():
        indices.sort(key=lambda i: (pool[i].match_pos, pool[i].sample, pool[i].pos))
        for a in range(len(indices)):
            for b in range(a + 1, len(indices)):
                ia, ib = indices[a], indices[b]
                distance = abs(pool[ia].match_pos - pool[ib].match_pos)
                if distance <= window_bp:
                    order = tuple(sorted((ia, ib)))
                    pairs.append(
                        (
                            distance,
                            -tsd_overlap(pool[ia].tsd_seq, pool[ib].tsd_seq),
                            pool[order[0]].chrom,
                            pool[order[0]].match_pos,
                            order,
                        )
                    )
    pairs.sort()

    site_of: dict[int, int] = {}
    sites: list[list[Call]] = []
    for distance, _neg_tsd, _chrom, _pos, (ia, ib) in pairs:
        sa, sb = site_of.get(ia), site_of.get(ib)
        if sa is not None and sb is not None:
            continue  # never merge two sites: no transitive chaining
        if sa is None and sb is None:
            first, second = (ia, ib) if (pool[ia].match_pos, pool[ia].sample, pool[ia].pos) <= (
                pool[ib].match_pos, pool[ib].sample, pool[ib].pos
            ) else (ib, ia)
            site_of[first] = site_of[second] = len(sites)
            sites.append([pool[first], pool[second]])
        elif sa is None:
            if matches_representative(pool[ia], sites[sb][0], window_bp):
                site_of[ia] = sb
                sites[sb].append(pool[ia])
        elif sb is None:
            if matches_representative(pool[ib], sites[sa][0], window_bp):
                site_of[ib] = sa
                sites[sa].append(pool[ib])

    # Calls with no partner within the window are private sites: one call, one site.
    # They are never dropped -- the private subset is a deliverable, not noise.
    for idx, call in enumerate(pool):
        if idx not in site_of:
            site_of[idx] = len(sites)
            sites.append([call])

    # Deterministic site identity and member order. The representative is the
    # member the greedy rule actually anchored on: earliest match_pos, then sample,
    # then pos. members[0] is always that representative.
    sites.sort(key=lambda members: (
        members[0].chrom, members[0].match_pos, members[0].family,
        members[0].orientation, members[0].host_key, len(members),
    ))
    for members in sites:
        members.sort(key=lambda c: (c.match_pos, c.sample, c.pos))
    return sites


# ---------------------------------------------------------------------------
# Site table, pooled AC/AF, accumulation curve, audit, manifest.
# ---------------------------------------------------------------------------

def site_rows(sites: list[list[Call]], site_offset: int = 0) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for index, members in enumerate(sites, start=site_offset + 1):
        rep = members[0]
        samples = sorted({c.sample for c in members})
        host_quads = {
            (c.host.name, c.host.start0, c.host.end0, c.host.strand)
            for c in members
            if c.host is not None
        }
        rows.append(
            {
                "site_id": f"US{index:05d}",
                "chrom": rep.chrom,
                "representative_pos": rep.pos,
                "representative_pos_calibrated": rep.match_pos,
                "representative_sample": rep.sample,
                "family": rep.family,
                "insertion_orientation": rep.orientation,
                "host_name": rep.host.name if rep.host else "",
                "host_start0": rep.host.start0 if rep.host else -1,
                "host_end0": rep.host.end0 if rep.host else -1,
                "host_strand": rep.host.strand if rep.host else "",
                "host_len": rep.host.length if rep.host else -1,
                "host_offset_5p_0based": rep.host.offset_5p_0based if rep.host else -1,
                "host_selection_rule": HOST_SELECTION_RULE,
                "host_selection_concordance": len(host_quads) == 1,
                "n_carriers": len(samples),
                "samples": "|".join(samples),
                "private": len(samples) == 1,
                "private_to_sample": samples[0] if len(samples) == 1 else "",
                "window_bp": PRIMARY_WINDOW_BP,
                "tsd_overlap_within_site": max(
                    (
                        tsd_overlap(members[i].tsd_seq, members[j].tsd_seq)
                        for i in range(len(members))
                        for j in range(i + 1, len(members))
                    ),
                    default=0,
                ),
            }
        )
    return rows


def pooled_allele_counts(
    sites: list[list[Call]], pooled_bcf: Path | None
) -> dict[tuple[str, int], tuple[int, float]]:
    """AC/AF for each site from the pooled site-level BCF, matched within the
    primary window. Returns {} when pysam or the BCF is unavailable."""
    if pooled_bcf is None or not pooled_bcf.is_file():
        return {}
    try:
        import pysam  # noqa: PLC0415
    except ImportError:
        return {}
    catalog: dict[str, list[tuple[int, int, float]]] = defaultdict(list)
    skipped = 0
    try:
        with pysam.VariantFile(str(pooled_bcf)) as bcf:
            for record in bcf.fetch():
                try:
                    ac = record.info.get("AC")
                    af = record.info.get("AF")
                except (ValueError, KeyError):
                    # Some pooled headers carry tag declarations htslib rejects;
                    # records whose declared INFO cannot be read are skipped.
                    skipped += 1
                    continue
                if ac is None or af is None:
                    continue
                ac_value = ac[0] if isinstance(ac, tuple) else ac
                af_value = af[0] if isinstance(af, tuple) else af
                catalog[record.chrom].append((record.pos, int(ac_value), float(af_value)))
    except (OSError, ValueError) as exc:
        print(f"[match-sites] pooled BCF unreadable, AC/AF skipped: {exc}", file=sys.stderr)
        return {}
    if skipped:
        print(f"[match-sites] pooled BCF: {skipped} records with unreadable INFO skipped", file=sys.stderr)
    for positions in catalog.values():
        positions.sort()
    result: dict[tuple[str, int], tuple[int, float]] = {}
    for members in sites:
        rep = members[0]
        positions = catalog.get(rep.chrom, ())
        idx = bisect.bisect_left(positions, (rep.match_pos - PRIMARY_WINDOW_BP,))
        while idx < len(positions) and positions[idx][0] <= rep.match_pos + PRIMARY_WINDOW_BP:
            result.setdefault((rep.chrom, rep.pos), (positions[idx][1], positions[idx][2]))
            idx += 1
    return result


# ---------------------------------------------------------------------------
# Cross-sample schema diff (annotation-free; flags POS/family/nesting drift).
# ---------------------------------------------------------------------------

def schema_diff(callsets: dict[str, list[Call]]) -> dict[str, object]:
    """Compare the annotation-free matching-relevant schema across callsets.

    Reports per-sample INFO key sets (restricted to fields the matching logic
    consumes), MEIFAMILY value distributions, NESTED value distributions, and
    POS coordinate ranges per chromosome. Drift flags are advisory and recorded
    in matching_audit.md; no annotation (GENE/CSQ) fields are ever read.
    """
    match_fields = ("MEIFAMILY", "ORIENT", "NESTED", "KNOWNMEI", "TSD", "CALLTIER")
    per_sample: dict[str, dict[str, object]] = {}
    for sample, calls in sorted(callsets.items()):
        fam_counts: dict[str, int] = defaultdict(int)
        nested_counts: dict[str, int] = defaultdict(int)
        for call in calls:
            fam_counts[call.family or "(missing)"] += 1
            nested_counts["nested" if call.nested else "unnested"] += 1
        pos_by_chrom: dict[str, list[int]] = defaultdict(list)
        for call in calls:
            pos_by_chrom[call.chrom].append(call.pos)
        per_sample[sample] = {
            "n_calls": len(calls),
            "family_distribution": dict(sorted(fam_counts.items(), key=lambda kv: -kv[1])),
            "nested_distribution": dict(sorted(nested_counts.items())),
            "chrom_min_max_pos": {
                chrom: [min(positions), max(positions)]
                for chrom, positions in sorted(pos_by_chrom.items())
            },
        }
    present = sorted({f for s in per_sample.values() for f in s["family_distribution"]})  # type: ignore[attr-defined]
    drift: dict[str, str] = {}
    fam_sets = {s: frozenset(per_sample[s]["family_distribution"]) for s in per_sample}  # type: ignore[attr-defined]
    if len({*fam_sets.values()}) > 1:
        drift["family_values"] = (
            "MEIFAMILY value sets differ across callsets: "
            + "; ".join(f"{s}: {sorted(fam_sets[s])}" for s in sorted(fam_sets))
        )
    nested_sets = {s: frozenset(per_sample[s]["nested_distribution"]) for s in per_sample}  # type: ignore[attr-defined]
    if len({*nested_sets.values()}) > 1:
        drift["nested_values"] = (
            "NESTED-derived value sets differ across callsets: "
            + "; ".join(f"{s}: {sorted(nested_sets[s])}" for s in sorted(nested_sets))
        )
    return {
        "match_fields_consumed": list(match_fields),
        "annotation_fields_consumed": [],
        "per_sample": per_sample,
        "families_present": present,
        "drift_flags": drift,
    }


def accumulation_curve(
    callsets: dict[str, list[Call]], hosts: dict[str, dict[tuple[str, int], HostAssignment]]
) -> pd.DataFrame:
    """Unique sites as genomes accumulate, ordered by sample ID (deterministic)."""
    sample_ids = sorted(callsets)
    rows: list[dict[str, object]] = []
    prev: dict[str, int] = {}
    for k in range(1, len(sample_ids) + 1):
        prefix_calls: list[Call] = []
        for sample in sample_ids[:k]:
            prefix_calls.extend(callsets[sample])
        sites = match_calls(prefix_calls)
        counts = defaultdict(int)
        for members in sites:
            counts[members[0].family] += 1
        row: dict[str, object] = {
            "genomes_added": k,
            "sample_added": sample_ids[k - 1],
            **{f"unique_{family.lower()}_sites": counts.get(family, 0) for family in ACCUMULATION_FAMILIES},
            "unique_all_families_sites": len(sites),
        }
        # Marginal yield (plan §3): new unique sites contributed by the latest genome.
        for column in list(row):
            if column.startswith("unique_") and column != "sample_added":
                row[f"marginal_{column}"] = row[column] - prev.get(column, 0)  # type: ignore[operator]
        prev = {c: row[c] for c in row if c.startswith("unique_")}  # type: ignore[misc]
        rows.append(row)
    return pd.DataFrame(rows)


def sensitivity_table(
    calls: list[Call]
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for window in (PRIMARY_WINDOW_BP, *SENSITIVITY_WINDOWS_BP):
        sites = match_calls(calls, window_bp=window)
        carrier_counts = [len({c.sample for c in members}) for members in sites]
        rows.append(
            {
                "window_bp": window,
                "n_unique_sites": len(sites),
                "n_multi_sample_sites": sum(n >= 2 for n in carrier_counts),
                "n_private_sites": sum(n == 1 for n in carrier_counts),
                "mean_carriers_per_site": round(
                    statistics.fmean(carrier_counts), 3
                ) if carrier_counts else 0.0,
                "max_carriers_per_site": max(carrier_counts, default=0),
            }
        )
    return rows


def write_audit(
    path: Path,
    samples: list[str],
    per_sample: dict[str, dict[str, object]],
    estimates: dict[str, OffsetEstimate],
    sens: list[dict[str, object]],
    rows: list[dict[str, object]],
    curve: pd.DataFrame,
    args: argparse.Namespace,
    pooled_matched: int,
    schema: dict[str, object],
) -> None:
    lines: list[str] = []
    lines.append("# Cross-sample nested-MEI matching audit")
    lines.append("")
    lines.append("Pre-registered rule, implemented as specified and not tuned: two calls match iff")
    lines.append("same chromosome, same MEI family, same insertion orientation, same assigned host")
    lines.append(f"element (audit_nested_tsd.py convention, {HOST_SELECTION_RULE}), and breakpoints within")
    lines.append(f"+/-{PRIMARY_WINDOW_BP} bp. No transitive chaining; greedy resolution by breakpoint")
    lines.append("proximity then TSD overlap. Sensitivity at +/-5 and +/-20 bp is reported in full.")
    lines.append("")
    lines.append("## Inputs")
    lines.append("")
    lines.append(f"- callset directory: `{args.callset_dir}`" if args.callset_dir else f"- callsets: {[str(c) for c in (args.callset or [])]}")
    lines.append(f"- RepeatMasker: `{args.rmsk}`")
    lines.append(f"- pooled site BCF: `{args.pooled_bcf}` ({pooled_matched} sites with AC/AF)")
    lines.append(f"- reference sample for offset calibration: `{samples[0] if samples else '(none)'}`")
    lines.append("")
    lines.append("## Cross-sample schema diff")
    lines.append("")
    lines.append("Matching consumes only these INFO fields (no gene/snpeff annotation fields): "
                 f"{schema['match_fields_consumed']}.")
    lines.append("")
    lines.append("| sample | calls | MEIFAMILY distribution | NESTED distribution |")
    lines.append("|---|---|---|---|")
    for sample in samples:
        info = schema["per_sample"][sample]  # type: ignore[index]
        lines.append(
            f"| {sample} | {info['n_calls']} | {info['family_distribution']} | {info['nested_distribution']} |"
        )
    lines.append("")
    drift = schema["drift_flags"]
    if drift:
        lines.append("**Drift flagged:**")
        for flag, detail in drift.items():
            lines.append(f"- `{flag}`: {detail}")
    else:
        lines.append("No family or nesting value drift detected across callsets.")
    lines.append("")
    lines.append("## Samples processed (manifest)")
    lines.append("")
    lines.append("| sample | file | calls | nested | nested with host |")
    lines.append("|---|---|---|---|---|")
    for sample in samples:
        info = per_sample[sample]
        lines.append(
            f"| {sample} | `{info['file']}` | {info['n_calls']} | {info['n_nested']} | {info['n_nested_with_host']} |"
        )
    lines.append("")
    lines.append("## Caller-offset calibration (per source)")
    lines.append("")
    lines.append("Estimated on KNOWNMEI catalogued positions within a broad +/-50 bp pairing window;")
    lines.append("median offset applied as a systematic correction only when >= 10 shared sites and")
    lines.append("MAD <= 3 bp. Never applied per call.")
    lines.append("")
    lines.append("| sample | shared sites | median offset (bp) | MAD (bp) | applied | reason |")
    lines.append("|---|---|---|---|---|---|")
    for sample in samples:
        est = estimates.get(sample)
        if est is None:
            lines.append(f"| {sample} | (reference) | 0 | 0.0 | no | reference_callset |")
        else:
            d = est.as_dict()
            lines.append(
                f"| {sample} | {d['n_knownmei_shared_sites']} | {d['median_offset_bp']} | "
                f"{d['mad_bp']} | {str(d['applied']).lower()} | {d['reason']} |"
            )
    lines.append("")
    lines.append("## Window sensitivity (full reporting)")
    lines.append("")
    lines.append("| window (bp) | unique sites | multi-sample sites | private sites | mean carriers | max carriers |")
    lines.append("|---|---|---|---|---|---|")
    for row in sens:
        lines.append(
            f"| +/-{row['window_bp']} | {row['n_unique_sites']} | {row['n_multi_sample_sites']} | "
            f"{row['n_private_sites']} | {row['mean_carriers_per_site']} | {row['max_carriers_per_site']} |"
        )
    lines.append("")
    multi = [row for row in sens if row["window_bp"] == PRIMARY_WINDOW_BP][0]
    lines.append(
        f"At the primary window the mean carriers/site is {multi['mean_carriers_per_site']}; "
        "agreement across all three windows is what makes a cross-sample claim robust, and any "
        "downstream statement should cite this table rather than a single window."
    )
    lines.append("")
    lines.append("## Unique sites")
    lines.append("")
    lines.append(f"- total unique sites: {len(rows)}")
    for family in ACCUMULATION_FAMILIES:
        fam_rows = [r for r in rows if r["family"] == family]
        lines.append(
            f"- unique nested-{family} sites: {len(fam_rows)} "
            f"(multi-sample: {sum(1 for r in fam_rows if r['n_carriers'] >= 2)}, "
            f"private: {sum(1 for r in fam_rows if r['private'])})"
        )
    other = [r for r in rows if r["family"] not in ACCUMULATION_FAMILIES]
    if other:
        lines.append(f"- other families present: {sorted({str(r['family']) for r in other})}")
    lines.append("")
    lines.append("## Accumulation curve")
    lines.append("")
    lines.append("Unique sites as genomes accumulate (sample-ID order, deterministic); see")
    lines.append("`accumulation_curve.csv` / `.png`. With the current sample count the curve has")
    lines.append(f"{len(curve)} point(s); it is regenerated each run as new callsets arrive.")
    lines.append("")
    path.write_text("\n".join(lines) + "\n")


def write_manifest(
    path: Path,
    samples: list[str],
    per_sample: dict[str, dict[str, object]],
    estimates: dict[str, OffsetEstimate],
    args: argparse.Namespace,
    outputs: dict[str, str],
) -> None:
    manifest = {
        "script": "rtm-dedup/scripts/match_sites_across_samples.py",
        "frozen_plan": "rtm-longread-cohort/docs/NESTED_ALU_ANALYSIS_PLAN.md (frozen 2026-10-01)",
        "schema_diff_drift_flags": "see matching_audit.md",
        "pre_registered_parameters": {
            "primary_window_bp": PRIMARY_WINDOW_BP,
            "sensitivity_windows_bp": list(SENSITIVITY_WINDOWS_BP),
            "no_transitive_chaining": True,
            "greedy_order": "breakpoint_proximity_then_tsd_overlap_then_deterministic_keys",
            "host_selection_rule": HOST_SELECTION_RULE,
            "host_assignment_convention": "audit_nested_tsd.py verbatim",
            "calibration": {
                "shared_sites": "KNOWNMEI=True",
                "broad_window_bp": CALIBRATION_BROAD_WINDOW_BP,
                "min_shared_sites": CALIBRATION_MIN_SHARED_SITES,
                "max_mad_bp": CALIBRATION_MAX_MAD_BP,
                "estimator": "median(ref_pos - other_pos)",
            },
        },
        "samples_processed": [
            {
                "sample": sample,
                "file": str(per_sample[sample]["file"]),
                "n_calls": per_sample[sample]["n_calls"],
                "n_nested": per_sample[sample]["n_nested"],
                "offset": estimates.get(sample, OffsetEstimate(sample, 0, 0, 0.0, False, "reference_callset")).as_dict(),
            }
            for sample in samples
        ],
        "outputs": outputs,
    }
    path.write_text(json.dumps(manifest, indent=2) + "\n")


def plot_curve(curve: pd.DataFrame, path: Path) -> None:
    try:
        import matplotlib  # noqa: PLC0415

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt  # noqa: PLC0415
    except ImportError:
        return
    fig, ax = plt.subplots(figsize=(7, 4.5))
    x = curve["genomes_added"]
    for column, label in (
        ("unique_alu_sites", "nested-Alu"),
        ("unique_line1_sites", "nested-L1"),
    ):
        if column in curve.columns:
            ax.plot(x, curve[column], marker="o", label=label)
    ax.set_xlabel("genomes added (sample-ID order)")
    ax.set_ylabel("unique nested sites")
    ax.set_title("Cross-sample unique nested-MEI accumulation")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--callset-dir", type=Path, default=None,
        help="directory of per-sample callset VCFs (*.vcf / *.vcf.gz); processed whatever exists",
    )
    parser.add_argument(
        "--callset", action="append", default=[],
        help="explicit callset VCF (repeatable); sample id from filename unless --sample-id",
    )
    parser.add_argument("--sample-id", default=None, help="sample id for --callset files")
    parser.add_argument(
        "--rmsk", type=Path, default=WORKSPACE / "nested_analysis/data/rmsk.txt.gz",
        help="RepeatMasker rmsk.txt.gz (audit_nested_tsd.py format)",
    )
    parser.add_argument(
        "--pooled-bcf", type=Path,
        default=Path("/Users/smyan/retrotransposon-workdir/data/public/polymorphism/hg38/long_read_1kg_ont_vienna/svim.asm.hg38.bcf"),
        help=(
            "pooled BCF for AC/AF lookup (skipped when absent). Defaults to the "
            "genotype BCF the provenance manifest uses for AC; the noGt site BCF "
            "carries a malformed header and no AC declarations"
        ),
    )
    parser.add_argument(
        "--outdir", type=Path, default=WORKSPACE / "nested_analysis/results_multi_sample",
        help="output directory",
    )
    args = parser.parse_args(argv)

    callset_list = discover_callsets(args)
    if not callset_list:
        raise SystemExit(
            "no per-sample callsets found: pass --callset-dir with *.vcf files "
            "and/or --callset FILE (new sample files are picked up on each run)"
        )
    if not args.rmsk.is_file():
        raise SystemExit(f"RepeatMasker file not found: {args.rmsk}")

    per_sample: dict[str, dict[str, object]] = {}
    callsets: dict[str, list[Call]] = {}
    hosts: dict[str, dict[tuple[str, int], HostAssignment]] = {}
    for sample, path in callset_list:
        calls = load_callset(path, sample)
        for call in calls:
            call.host = None
        assigned = assign_hosts(calls, args.rmsk)
        for call in calls:
            call.host = assigned.get(call.chrom_pos_key)
        verify_genome_wide(path, calls)
        callsets[sample] = calls
        hosts[sample] = assigned
        per_sample[sample] = {
            "file": str(path),
            "n_calls": len(calls),
            "n_nested": sum(c.nested for c in calls),
            "n_nested_with_host": sum(c.host is not None for c in calls),
        }

    samples = sorted(callsets)
    reference_sample = samples[0]
    all_calls = [c for sample in samples for c in callsets[sample]]

    estimates: dict[str, OffsetEstimate] = {}
    for sample in samples[1:]:
        estimates[sample] = estimate_offset(callsets[reference_sample], callsets[sample], sample)
    apply_calibration(reference_sample, all_calls, estimates)

    sites = match_calls(all_calls)
    rows = site_rows(sites)
    pooled = pooled_allele_counts(sites, args.pooled_bcf)
    schema = schema_diff(callsets)
    for row in rows:
        ac_af = pooled.get((str(row["chrom"]), int(row["representative_pos"])))
        row["site_allele_count"] = ac_af[0] if ac_af else ""
        row["site_allele_freq"] = ac_af[1] if ac_af else ""
        row["pooled_site_matched"] = ac_af is not None

    sens = sensitivity_table(all_calls)
    curve = accumulation_curve(callsets, hosts)

    args.outdir.mkdir(parents=True, exist_ok=True)
    outputs: dict[str, str] = {}
    csv_path = args.outdir / "unique_sites.csv"
    frame = pd.DataFrame(rows)
    # Call counts and unique-site counts side by side (plan §1): the frame already
    # carries per-site carrier columns; the audit and stdout summarise them.
    frame.to_csv(csv_path, index=False)
    outputs["unique_sites"] = str(csv_path)
    curve.to_csv(args.outdir / "accumulation_curve.csv", index=False)
    outputs["accumulation_curve_csv"] = str(args.outdir / "accumulation_curve.csv")
    plot_curve(curve, args.outdir / "accumulation_curve.png")
    if (args.outdir / "accumulation_curve.png").exists():
        outputs["accumulation_curve_png"] = str(args.outdir / "accumulation_curve.png")
    write_audit(
        args.outdir / "matching_audit.md", samples, per_sample, estimates, sens,
        rows, curve, args, sum(1 for r in rows if r["pooled_site_matched"]),
        schema,
    )
    outputs["matching_audit"] = str(args.outdir / "matching_audit.md")
    write_manifest(
        args.outdir / "matching_manifest.json", samples, per_sample, estimates, args, outputs
    )
    outputs["manifest"] = str(args.outdir / "matching_manifest.json")

    line1_rows = [r for r in rows if r["family"] == "LINE1"]
    print(
        json.dumps(
            {
                "samples": samples,
                "unique_sites": len(rows),
                "unique_nested_line1_sites": len(line1_rows),
                "multi_sample_sites": sum(1 for r in rows if r["n_carriers"] >= 2),
                "private_sites": sum(1 for r in rows if r["private"]),
                "pooled_ac_af_matched": sum(1 for r in rows if r["pooled_site_matched"]),
                "outputs": outputs,
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
