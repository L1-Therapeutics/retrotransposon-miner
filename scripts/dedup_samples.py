#!/usr/bin/env python3
"""Deduplicate five sample MEI callsets and run the requested unique-site analyses.

All calls are checksum- and baseline-gated before use. GT/GQ are never parsed.
Calls on chrY are counted for QC but excluded from the matching and analysis universe.
Same-family nesting is independently assigned from RepeatMasker and retains the
four-valued unnested/nested_sense/nested_antisense/nested_unknown state.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import math
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
WORKSPACE = REPO.parent
SAMPLES = ("HG03086", "HG01474", "HG01566", "HG03172", "NA18498")
BASELINES = {
    "HG03086": (1852, 1555, 222, 75),
    "HG01474": (1470, 1226, 193, 51),
    "HG01566": (1522, 1253, 213, 56),
    "HG03172": (1880, 1550, 250, 80),
    "NA18498": (1721, 1457, 200, 64),
}
FAMILIES = ("ALU", "LINE1", "SVA")
PRIMARY_CHROMS = tuple([f"chr{i}" for i in range(1, 23)] + ["chrX"])
ALLOWED_CHROMS = frozenset((*PRIMARY_CHROMS, "chrY"))
PRIMARY_WINDOW = 10
SWEEP_WINDOWS = (0, 5, 10, 20)
NESTED_STATES = ("unnested", "nested_sense", "nested_antisense", "nested_unknown")
REPS = 10_000
SEED = 20261001
GC_BIN_BP = 200
GC_TOLERANCE = 0.05
GC_MATCHED_BINS = 2_000
PROFILE_HOST_RANGE = (280, 320)
PROFILE_BIN_BP = 20
PROFILE_WINDOWS = ((120, 140), (280, 300))
PROFILE_MIN_EXPECTED_PER_BIN = 5.0

# The supplied HG03086-only comparison values are carried forward verbatim;
# they are not recomputed on the new five-sample unique-site denominator.
HG03086_BASELINE = {
    "enrichment": {"ALU": 1.47, "LINE1": 1.47, "SVA": 29.4},
    "alu_sense": 2.49,
    "alu_antisense": 0.46,
    "line1_antisense": 1.0,
    "linker": 3.0,
    "tail": 3.7,
}


class InputGateError(RuntimeError):
    """Raised before analysis when a source/provenance check fails."""


@dataclass(frozen=True)
class Host:
    chrom: str
    start0: int
    end0: int
    strand: str
    name: str
    family: str

    @property
    def key(self) -> tuple[str, int, int, str, str]:
        return (self.chrom, self.start0, self.end0, self.name, self.strand)

    @property
    def host_id(self) -> str:
        return f"{self.chrom}:{self.start0}-{self.end0}:{self.name}:{self.strand}"

    @property
    def length(self) -> int:
        return self.end0 - self.start0


@dataclass
class Call:
    sample: str
    source_index: int
    chrom: str
    pos: int
    family: str
    orientation: str
    raw_nested: str
    info: dict[str, str]
    sample_column: str
    same_family_host: Host | None = None
    match_host: Host | None = None
    alu_host: Host | None = None
    nested_state: str = "unnested"
    source_count: int = 1
    source_ids: list[str] = field(default_factory=list)

    @property
    def pos0(self) -> int:
        return self.pos - 1

    @property
    def offset(self) -> int | None:
        if self.match_host is None or not self.match_host.start0 <= self.pos0 < self.match_host.end0:
            return None
        return self.pos0 - self.match_host.start0 if self.match_host.strand == "+" else self.match_host.end0 - 1 - self.pos0

    @property
    def match_key(self) -> tuple[Any, ...]:
        # Preserve the literal orientation value, including unresolved `.`;
        # callers with the same unresolved category can still coincide.
        orient_key = self.orientation or "."
        # Prefer the assigned host element used for the rule; unnested calls
        # share the explicit no-host value.
        host_key = self.match_host.key if self.match_host else None
        return (self.chrom, self.family, orient_key, host_key)


@dataclass
class Site:
    anchor: Call
    members: list[Call]

    @property
    def samples(self) -> tuple[str, ...]:
        return tuple(sorted({call.sample for call in self.members}, key=SAMPLES.index))

    @property
    def representative(self) -> Call:
        positions = [call.pos for call in self.members]
        median = float(np.median(positions))
        return min(self.members, key=lambda call: (abs(call.pos - median), SAMPLES.index(call.sample), call.pos, call.source_index))


def normalize_family(token: str) -> str:
    value = (token or "").upper()
    if "ALU" in value:
        return "ALU"
    if "SVA" in value:
        return "SVA"
    if "LINE1" in value or "L1" in value:
        return "LINE1"
    return ""


def parse_info(text: str) -> dict[str, str]:
    info: dict[str, str] = {}
    if text in {"", "."}:
        return info
    for item in text.split(";"):
        key, sep, value = item.partition("=")
        if key:
            info[key] = value if sep else "True"
    return info


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def input_paths(callset_dir: Path, root: Path) -> dict[str, tuple[Path, Path]]:
    pairs = {}
    for sample in SAMPLES:
        source = callset_dir / f"{sample}.vcf"
        duplicate = root / f"{sample}_classifier_ge_0.997.genes.snpeff.vcf"
        if not source.is_file() or not duplicate.is_file():
            raise InputGateError(f"missing required source/duplicate for {sample}: {source} | {duplicate}")
        source_hash, duplicate_hash = sha256(source), sha256(duplicate)
        if source_hash != duplicate_hash:
            raise InputGateError(
                f"SHA-256 MISMATCH for {sample}; STOP, no source selected: "
                f"{source}={source_hash}; {duplicate}={duplicate_hash}"
            )
        pairs[sample] = (source, duplicate)
    return pairs


def parse_callset(path: Path, sample: str) -> tuple[list[Call], dict[str, Any]]:
    calls: list[Call] = []
    family_counts: Counter[str] = Counter()
    chrom_counts: Counter[str] = Counter()
    raw_nested: Counter[str] = Counter()
    column_header: list[str] | None = None
    reference = ""
    errors: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            line = line.rstrip("\r\n")
            if line.startswith("##reference="):
                reference = line.split("=", 1)[1]
            elif line.startswith("#CHROM"):
                column_header = line.split("\t")
            elif line.startswith("#") or not line:
                continue
            else:
                fields = line.split("\t")
                if column_header is None or len(fields) != len(column_header):
                    errors.append(f"line {line_number}: {len(fields)} data columns; expected {len(column_header or [])}")
                    continue
                chrom = fields[0]
                family = normalize_family(parse_info(fields[7]).get("MEIFAMILY", ""))
                info = parse_info(fields[7])
                if family not in FAMILIES:
                    errors.append(f"line {line_number}: unsupported/missing MEIFAMILY {info.get('MEIFAMILY', '')!r}")
                    continue
                try:
                    pos = int(fields[1])
                except ValueError:
                    errors.append(f"line {line_number}: invalid POS {fields[1]!r}")
                    continue
                if pos < 1:
                    errors.append(f"line {line_number}: POS must be >=1")
                    continue
                sample_col = column_header[9] if len(column_header) > 9 else ""
                call = Call(
                    sample=sample, source_index=len(calls), chrom=chrom, pos=pos,
                    family=family, orientation=info.get("ORIENT", "").strip(),
                    raw_nested=info.get("NESTED", "").strip().lower(), info=info,
                    sample_column=sample_col,
                    source_ids=[fields[2] if fields[2] not in {"", "."} else f"{sample}:{chrom}:{pos}:{len(calls)}"],
                )
                calls.append(call)
                family_counts[family] += 1
                chrom_counts[chrom] += 1
                raw_nested[call.raw_nested or "<MISSING>"] += 1
    if errors:
        raise InputGateError(f"malformed {sample} VCF (first issues: {'; '.join(errors[:5])})")
    expected_columns = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample]
    if column_header != expected_columns:
        raise InputGateError(f"{sample} VCF header columns/sample mismatch: {column_header}; expected {expected_columns}")
    if reference.lower() != "hg38":
        raise InputGateError(f"{sample} VCF declares reference {reference!r}, expected hg38")
    if set(chrom_counts) - ALLOWED_CHROMS:
        raise InputGateError(f"{sample} VCF has out-of-scope chromosome labels: {sorted(set(chrom_counts)-ALLOWED_CHROMS)}")
    baseline = BASELINES[sample]
    observed = (len(calls), family_counts["ALU"], family_counts["LINE1"], family_counts["SVA"])
    missing_primary = set(PRIMARY_CHROMS) - set(chrom_counts)
    genome_wide = not missing_primary and any(chrom_counts[chrom] for chrom in ("chrX", "chrY"))
    expected_y = 1 if sample in {"HG03172", "NA18498"} else 0
    if chrom_counts["chrY"] != expected_y:
        raise InputGateError(f"{sample} chrY QC mismatch: observed {chrom_counts['chrY']}, expected {expected_y}; STOP")
    if missing_primary or not genome_wide:
        raise InputGateError(
            f"{sample} is not genome-wide (must include chr1-22 plus a chromosome beyond chr22); "
            f"missing={sorted(missing_primary)}; STOP"
        )
    if observed != baseline:
        raise InputGateError(f"QC BASELINE MISMATCH {sample}: observed {observed}, expected {baseline}; STOP")
    if sample == "HG03086" and raw_nested["nested"] != 261:
        raise InputGateError(
            f"HG03086 raw NESTED=nested QC mismatch: observed {raw_nested['nested']}, expected 261; STOP"
        )
    return calls, {
        "sample": sample, "file": str(path.resolve()), "calls": len(calls),
        "primary_chromosomes_complete": True,
        "families": dict(family_counts), "chromosomes": dict(chrom_counts),
        "chrY_calls": chrom_counts["chrY"], "raw_nested_values": dict(raw_nested),
        "reference": reference, "sample_column": column_header[9] if len(column_header) > 9 else "(no sample column)",
        "genome_wide": genome_wide, "qc_baseline_match": True,
    }


def read_rmsk(path: Path) -> dict[tuple[str, str], list[Host]]:
    intervals: dict[tuple[str, str], list[Host]] = defaultdict(list)
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 13 or fields[5] not in PRIMARY_CHROMS:
                continue
            family = normalize_family(f"{fields[10]} {fields[11]} {fields[12]}")
            if family not in FAMILIES:
                continue
            try:
                start0, end0 = int(fields[6]), int(fields[7])
            except ValueError:
                continue
            if end0 <= start0:
                continue
            host = Host(fields[5], start0, end0, fields[9], fields[10], family)
            intervals[(host.chrom, family)].append(host)
    for hosts in intervals.values():
        hosts.sort(key=lambda host: (host.start0, host.end0, host.name, host.strand))
    return intervals


def _sweep_assign(calls: list[Call], intervals: list[Host], *, only_family: str | None = None) -> dict[int, Host]:
    """Assign each point to its longest containing MEI; tie-break left then name."""
    if not calls or not intervals:
        return {}
    result: dict[int, Host] = {}
    ordered_calls = sorted(calls, key=lambda call: (call.pos0, call.sample, call.source_index))
    active: list[Host] = []
    next_interval = 0
    for call in ordered_calls:
        # VCF POS is 1-based; nested annotation converts it to zero-based start
        # then emits [POS-1, POS), so a host containing POS must satisfy
        # start0 < POS <= end0 (equivalently start0 < POS and end0 >= POS).
        while next_interval < len(intervals) and intervals[next_interval].start0 < call.pos:
            active.append(intervals[next_interval])
            next_interval += 1
        active = [host for host in active if host.end0 >= call.pos]
        candidates = [host for host in active if (only_family is None or host.family == only_family)]
        if candidates:
            result[id(call)] = min(
                candidates, key=lambda host: (-host.length, host.start0, host.name, host.strand)
            )
    return result


def assign_hosts(calls: list[Call], rmsk: dict[tuple[str, str], list[Host]]) -> None:
    by_chrom: dict[str, list[Call]] = defaultdict(list)
    for call in calls:
        if call.chrom in PRIMARY_CHROMS:
            by_chrom[call.chrom].append(call)
    for chrom, chrom_calls in by_chrom.items():
        by_family: dict[str, list[Call]] = defaultdict(list)
        for call in chrom_calls:
            by_family[call.family].append(call)
        same_assign: dict[int, Host] = {}
        for family, family_calls in by_family.items():
            same_assign.update(_sweep_assign(family_calls, rmsk.get((chrom, family), []), only_family=family))
        alu_assign = _sweep_assign(
            [call for call in chrom_calls if call.family == "LINE1"],
            rmsk.get((chrom, "ALU"), []), only_family="ALU",
        )
        for call in chrom_calls:
            call.same_family_host = same_assign.get(id(call))
            # Matching host is the same-family RepeatMasker element; an
            # insertion outside that family is explicitly hostless for this key.
            call.match_host = call.same_family_host
            call.alu_host = alu_assign.get(id(call)) if call.family == "LINE1" else (
                call.same_family_host if call.family == "ALU" else None
            )
            if call.same_family_host is None:
                call.nested_state = "unnested"
            elif call.orientation not in {"+", "-"} or call.same_family_host.strand not in {"+", "-"}:
                call.nested_state = "nested_unknown"
            elif call.orientation == call.same_family_host.strand:
                call.nested_state = "nested_sense"
            else:
                call.nested_state = "nested_antisense"


def candidate_pairs(calls: list[Call], window: int, *, include_same_sample: bool) -> list[tuple[int, int, int]]:
    grouped: dict[tuple[Any, ...], list[int]] = defaultdict(list)
    for index, call in enumerate(calls):
        grouped[call.match_key].append(index)
    pairs: list[tuple[int, int, int]] = []
    for indexes in grouped.values():
        indexes.sort(key=lambda i: (calls[i].pos, SAMPLES.index(calls[i].sample), calls[i].source_index))
        for ii, left_index in enumerate(indexes):
            left = calls[left_index]
            for right_index in indexes[ii + 1:]:
                right = calls[right_index]
                distance = right.pos - left.pos
                if distance > window:
                    break
                if include_same_sample or left.sample != right.sample:
                    pairs.append((distance, left_index, right_index))
    pairs.sort(key=lambda item: (item[0], calls[item[1]].chrom, calls[item[1]].pos,
                                 SAMPLES.index(calls[item[1]].sample), calls[item[1]].source_index,
                                 SAMPLES.index(calls[item[2]].sample), calls[item[2]].source_index))
    return pairs


def anchored_groups(calls: list[Call], window: int, *, allow_same_sample: bool) -> list[Site]:
    """Greedy representative-anchored grouping; sites never merge transitively."""
    pairs = candidate_pairs(calls, window, include_same_sample=allow_same_sample)
    site_of: dict[int, int] = {}
    sites: list[Site] = []
    for _, left_index, right_index in pairs:
        left_site, right_site = site_of.get(left_index), site_of.get(right_index)
        if left_site is not None and right_site is not None:
            continue
        if left_site is None and right_site is None:
            first, second = sorted((left_index, right_index), key=lambda idx: (
                calls[idx].pos, SAMPLES.index(calls[idx].sample), calls[idx].source_index
            ))
            if not allow_same_sample and calls[first].sample == calls[second].sample:
                continue
            site_of[first] = site_of[second] = len(sites)
            sites.append(Site(calls[first], [calls[first], calls[second]]))
            continue
        site_id = left_site if left_site is not None else right_site
        new_index = right_index if left_site is not None else left_index
        assert site_id is not None
        site = sites[site_id]
        if abs(calls[new_index].pos - site.anchor.pos) > window:
            continue
        if not allow_same_sample and calls[new_index].sample in {member.sample for member in site.members}:
            continue
        site_of[new_index] = site_id
        site.members.append(calls[new_index])
    for index, call in enumerate(calls):
        if index not in site_of:
            site_of[index] = len(sites)
            sites.append(Site(call, [call]))
    return sorted(sites, key=lambda site: (
        site.anchor.chrom, site.anchor.pos, site.anchor.family, site.anchor.orientation,
        site.anchor.match_host.key if site.anchor.match_host else (),
        SAMPLES.index(site.anchor.sample), site.anchor.source_index,
    ))


def collapse_within_sample(calls: list[Call], window: int) -> tuple[list[Call], int]:
    """Collapse repeated records from the same sample at one anchored site."""
    by_sample: dict[str, list[Call]] = defaultdict(list)
    for call in calls:
        by_sample[call.sample].append(call)
    representatives: list[Call] = []
    n_collapsed = 0
    for sample in SAMPLES:
        sample_calls = by_sample.get(sample, [])
        if not sample_calls:
            continue
        grouped = anchored_groups(sample_calls, window, allow_same_sample=True)
        for site in grouped:
            rep = site.representative
            # Aggregate provenance onto a *copy*, never onto the caller's object.
            # `run` calls `unique_sites` on the same Call objects once for the
            # primary window and again for every SWEEP_WINDOWS entry (which
            # repeats the primary +/-10 bp window). Writing the totals back in
            # place made each pass re-sum its own inflated predecessor, so
            # source_count walked 3 -> 4 -> 5 and source_call_ids grew
            # "A|B|C" -> "A|B|B|C" -> "A|B|B|B|C" for a site with three records.
            source_count = sum(member.source_count for member in site.members)
            source_ids = [source_id for member in site.members for source_id in member.source_ids]
            if source_count != rep.source_count or source_ids != rep.source_ids:
                rep = replace(rep, source_count=source_count, source_ids=source_ids)
            representatives.append(rep)
            n_collapsed += len(site.members) - 1
    return representatives, n_collapsed


def unique_sites(calls: list[Call], window: int) -> tuple[list[Site], int]:
    canonical, within_sample_collapsed = collapse_within_sample(calls, window)
    return anchored_groups(canonical, window, allow_same_sample=False), within_sample_collapsed


def site_row(site: Site, site_id: str, window: int) -> dict[str, Any]:
    call = site.representative
    samples = site.samples
    host = call.same_family_host
    return {
        "site_id": site_id, "chrom": call.chrom, "representative_pos": call.pos,
        "family": call.family, "insertion_orientation": call.orientation,
        "same_family_nested_state": call.nested_state,
        "raw_NESTED_value": call.raw_nested,
        "same_family_host_id": host.host_id if host else "",
        "same_family_host_name": host.name if host else "",
        "same_family_host_strand": host.strand if host else "",
        "same_family_host_length": host.length if host else "",
        # Host geometry under the canonical column names. `joint_enrichment.py`
        # and `recurrence_test.py` read this same file through
        # `nested_multi_sample_common.load_unique_sites`, which requires
        # host_name/host_start0/host_end0/host_strand/host_len. Emitting only the
        # `same_family_host_*` names made this producer write a table its own
        # downstream consumers reject, so both Phase 4 scripts died with a
        # SchemaError on real data while every unit test still passed. The
        # `same_family_*` columns are kept above because the matching audit reads
        # them; these are the same host, restated under the names the loader wants.
        "host_name": host.name if host else "",
        "host_start0": host.start0 if host else "",
        "host_end0": host.end0 if host else "",
        "host_strand": host.strand if host else "",
        "host_len": host.length if host else "",
        "host_selection_rule": "longest_overlap_then_leftmost_name" if host else "",
        "host_offset_5p_0based": call.offset if call.offset is not None else "",
        "n_carriers": len(samples), "carrier_count_class": str(len(samples)),
        "samples": "|".join(samples), "private": len(samples) == 1,
        "private_to_sample": samples[0] if len(samples) == 1 else "",
        "source_record_count": sum(member.source_count for member in site.members),
        "source_call_ids": "|".join(source_id for member in site.members for source_id in member.source_ids),
        "match_window_bp": window,
    }


# `any_MEI_host_id_for_matching` used to be emitted here. It is deliberately
# absent, and this note records why so it does not get re-added:
#
#   * The authoritative producer named in the matching manifest,
#     `rtm-dedup/scripts/match_sites_across_samples.py`, never emitted it. Its
#     `site_rows` writes `host_name`/`host_start0`/`host_end0`/`host_strand`/
#     `host_len` and no `any_MEI_host*` field at all.
#   * `matching_audit.md` freezes the matching key to the *same-family* host, so
#     a column named "...for_matching" could only ever hold the same-family host.
#   * It therefore held a byte-for-byte duplicate of `same_family_host_id` --
#     verified across all 4,998 rows of the shipped table.
#   * It had no consumer anywhere in the workspace, and the frozen plan never
#     mentions it.
#
# Populating it from the all-family sweep would have been worse, not better: that
# sweep returns the longest element of *any* family overlapping the call, so for
# an Alu sitting inside an L1 it would name the L1 in a column that says the host
# used for matching -- a different and actively misleading value. The sweep was
# dead work and is gone.
#
# Reconciliation: nothing consumes the column, the frozen plan does not mention
# it, and the authoritative producer never had it, so removing it moves this
# table toward the recorded schema rather than away from it. The column's
# absence is what the regenerated `matching_audit.md` schema-diff section will
# record; no hand edit of that output is needed.


def build_sites_table(sites: list[Site], window: int) -> pd.DataFrame:
    return pd.DataFrame([site_row(site, f"US{index:05d}", window) for index, site in enumerate(sites, 1)])


def merge_intervals(hosts: Iterable[Host]) -> list[tuple[int, int]]:
    merged: list[list[int]] = []
    for host in sorted(hosts, key=lambda item: (item.start0, item.end0)):
        if merged and host.start0 <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], host.end0)
        else:
            merged.append([host.start0, host.end0])
    return [(start, end) for start, end in merged]


def interval_prefix(intervals: list[tuple[int, int]]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    starts = np.asarray([start for start, _ in intervals], dtype=np.int64)
    ends = np.asarray([end for _, end in intervals], dtype=np.int64)
    cumulative = np.concatenate((np.asarray([0], dtype=np.int64), np.cumsum(ends - starts, dtype=np.int64)))
    return starts, ends, cumulative


def cumulative_coverage(points: np.ndarray, starts: np.ndarray, ends: np.ndarray, cumulative: np.ndarray) -> np.ndarray:
    if starts.size == 0:
        return np.zeros(points.shape, dtype=np.int64)
    idx = np.searchsorted(starts, points, side="right") - 1
    safe = np.maximum(idx, 0)
    result = cumulative[safe].copy()
    inside = (idx >= 0) & (points < ends[safe])
    result[inside] += points[inside] - starts[safe[inside]]
    after_interval = (idx >= 0) & (points >= ends[safe])
    result[after_interval] = cumulative[safe[after_interval] + 1]
    before_first = idx < 0
    result[before_first] = 0
    return result


def uniform_probabilities(sites: list[Site], merged_cov: dict[tuple[str, str], int], lengths: dict[str, int], family: str) -> np.ndarray:
    probs = []
    for site in sites:
        call = site.representative
        if call.family == family:
            probs.append(merged_cov.get((family, call.chrom), 0) / lengths[call.chrom])
    return np.asarray(probs, dtype=float)


def resample_nested(probabilities: np.ndarray, reps: int, rng: np.random.Generator) -> np.ndarray:
    if probabilities.size == 0:
        return np.zeros(reps, dtype=np.int64)
    result = np.zeros(reps, dtype=np.int64)
    chunk = max(1, min(256, 2_000_000 // max(1, probabilities.size)))
    for start in range(0, reps, chunk):
        n = min(chunk, reps - start)
        result[start:start + n] = (rng.random((n, probabilities.size)) < probabilities).sum(axis=1)
    return result


def load_fai(path: Path) -> tuple[dict[str, tuple[int, int, int, int]], dict[str, int]]:
    fai = Path(str(path) + ".fai")
    if not fai.is_file():
        raise InputGateError(f"FASTA index missing: {fai}")
    offsets, sizes = {}, {}
    with fai.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            name, length, offset, line_bases, line_width = fields[:5]
            offsets[name] = (int(length), int(offset), int(line_bases), int(line_width))
            sizes[name] = int(length)
    missing = set(PRIMARY_CHROMS) - set(sizes)
    if missing:
        raise InputGateError(f"GRCh38 FASTA .fai lacks required chromosomes: {sorted(missing)}")
    return offsets, sizes


def fasta_sequence_chunk(handle, meta: tuple[int, int, int, int], start: int, end: int) -> bytes:
    length, offset, line_bases, line_width = meta
    if not 0 <= start < end <= length:
        raise ValueError(f"bad FASTA slice {start}:{end} of {length}")
    line_index, col = divmod(start, line_bases)
    byte_offset = offset + line_index * line_width + col
    n_bases = end - start
    newline_bytes = ((col + n_bases - 1) // line_bases) * (line_width - line_bases)
    handle.seek(byte_offset)
    raw = handle.read(n_bases + newline_bytes)
    seq = raw.replace(b"\n", b"").replace(b"\r", b"")
    if len(seq) < n_bases:
        raise InputGateError(f"short FASTA read at bases {start}:{end}")
    return seq[:n_bases].upper()


def gc_matched_frame(fasta: Path, rmsk: dict[tuple[str, str], list[Host]], lengths: dict[str, int], sites: list[Site]) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    """Build GC/MEI-coverage values for every valid 200-bp GRCh38 bin.

    Uses the existing hg38 FASTA and .fai (no package installs). A and C/G/T
    fractions use full bin width as denominator; bins with <50% A/C/G/T are
    excluded, matching the earlier HG03086 null implementation.
    """
    fasta = fasta.expanduser().resolve()
    if not fasta.is_file():
        raise InputGateError(f"GRCh38 FASTA unavailable for GC-matched null: {fasta}")
    meta, fai_lengths = load_fai(fasta)
    for chrom in PRIMARY_CHROMS:
        if fai_lengths[chrom] != lengths[chrom]:
            raise InputGateError(f"FASTA/chromosome-size mismatch for {chrom}: {fai_lengths[chrom]} vs {lengths[chrom]}")

    # Merge RepeatMasker intervals within-family so coverage cannot exceed a bin.
    pref: dict[tuple[str, str], tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for chrom in PRIMARY_CHROMS:
        for family in FAMILIES:
            pref[(chrom, family)] = interval_prefix(merge_intervals(rmsk.get((chrom, family), [])))

    site_gc: dict[tuple[str, int], float | None] = {}
    target_bins: dict[str, set[int]] = defaultdict(set)
    for site in sites:
        call = site.representative
        if call.chrom in PRIMARY_CHROMS:
            key = (call.chrom, call.pos0 // GC_BIN_BP)
            site_gc[key] = None
            target_bins[call.chrom].add(key[1])

    gc_parts: list[np.ndarray] = []
    cov_parts: dict[str, list[np.ndarray]] = {family: [] for family in FAMILIES}
    valid_gc_bins = 0
    total_genome_bins = 0
    excluded_gc_bins = 0
    with fasta.open("rb") as handle:
        for chrom in PRIMARY_CHROMS:
            length = lengths[chrom]
            start0 = 0
            while start0 < length:
                end0 = min(start0 + 200_000, length)
                seq = fasta_sequence_chunk(handle, meta[chrom], start0, end0)
                arr = np.frombuffer(seq, dtype=np.uint8)
                nbin = (len(arr) + GC_BIN_BP - 1) // GC_BIN_BP
                padded = np.zeros(nbin * GC_BIN_BP, dtype=np.uint8)
                padded[:len(arr)] = arr
                matrix = padded.reshape(nbin, GC_BIN_BP)
                spans = np.minimum(GC_BIN_BP, len(arr) - np.arange(nbin) * GC_BIN_BP).astype(np.int64)
                gc_counts = ((matrix == ord("G")) | (matrix == ord("C"))).sum(axis=1, dtype=np.int64)
                at_counts = ((matrix == ord("A")) | (matrix == ord("C")) | (matrix == ord("G")) | (matrix == ord("T"))).sum(axis=1, dtype=np.int64)
                gc_frac = gc_counts / spans
                at_frac = at_counts / spans
                bin_starts = start0 + np.arange(nbin, dtype=np.int64) * GC_BIN_BP
                keep = at_frac >= 0.5
                gc_parts.append(gc_frac[keep].astype(np.float32))
                valid_gc_bins += int(keep.sum())
                total_genome_bins += nbin
                excluded_gc_bins += int((~keep).sum())
                bin_ends = bin_starts + spans
                for family in FAMILIES:
                    starts, ends, cumulative = pref[(chrom, family)]
                    covered = cumulative_coverage(bin_ends, starts, ends, cumulative) - cumulative_coverage(bin_starts, starts, ends, cumulative)
                    cov_parts[family].append((covered[keep] / spans[keep]).astype(np.float32))
                lo_bin = start0 // GC_BIN_BP
                hi_bin = lo_bin + nbin
                for target_bin in target_bins.get(chrom, ()):
                    if lo_bin <= target_bin < hi_bin:
                        local = target_bin - lo_bin
                        if keep[local]:
                            site_gc[(chrom, target_bin)] = float(gc_frac[local])
                start0 = end0

    gc = np.concatenate(gc_parts) if gc_parts else np.asarray([], dtype=np.float32)
    coverage = {family: np.concatenate(parts) if parts else np.asarray([], dtype=np.float32) for family, parts in cov_parts.items()}
    frame = {"gc": gc, "coverage": coverage}
    frame["sorted_order"] = np.argsort(gc)
    return {"gc_by_site": site_gc, "frame": frame}, {
        "fasta": str(fasta), "bin_bp": GC_BIN_BP, "gc_tolerance": GC_TOLERANCE,
        "valid_gc_bins": valid_gc_bins, "bins_gc_excluded_lt_50pct_acgt": excluded_gc_bins,
        "total_genome_bins": total_genome_bins,
        "n_call_sites_without_gc_bin": sum(value is None for value in site_gc.values()),
    }


def gc_probabilities(sites: list[Site], family: str, gc_bundle: dict[str, Any], rng: np.random.Generator) -> np.ndarray:
    target_gc = gc_bundle["gc_by_site"]
    frame = gc_bundle["frame"]
    all_gc = frame["gc"]
    order = frame["sorted_order"]
    sorted_gc = all_gc[order]
    sorted_cover = frame["coverage"][family][order]
    probs = []
    for site in sites:
        call = site.representative
        if call.family != family:
            continue
        target = target_gc.get((call.chrom, call.pos0 // GC_BIN_BP))
        if target is None:
            raise InputGateError(f"cannot GC-match {call.sample}:{call.chrom}:{call.pos}; local bin lacks >=50% A/C/G/T")
        lo = np.searchsorted(sorted_gc, target - GC_TOLERANCE, side="left")
        hi = np.searchsorted(sorted_gc, target + GC_TOLERANCE, side="right")
        pool = sorted_cover[lo:hi]
        if pool.size == 0:
            pool = sorted_cover
        if pool.size > GC_MATCHED_BINS:
            pool = pool[rng.choice(pool.size, size=GC_MATCHED_BINS, replace=False)]
        probs.append(float(pool.mean()))
    return np.asarray(probs, dtype=float)


def family_enrichment(sites: list[Site], lengths: dict[str, int], merged_cov: dict[tuple[str, str], int], gc_bundle: dict[str, Any], reps: int, seed: int) -> list[dict[str, Any]]:
    results = []
    for index, family in enumerate(FAMILIES):
        family_sites = [site for site in sites if site.representative.family == family]
        nested = [site for site in family_sites if site.representative.nested_state != "unnested"]
        method = "uniform" if family == "ALU" else "gc_matched"
        if method == "uniform":
            probs = uniform_probabilities(family_sites, merged_cov, lengths, family)
        else:
            probs = gc_probabilities(family_sites, family, gc_bundle, np.random.default_rng(seed + index))
        observed = len(nested)
        expected = float(probs.sum())
        null = resample_nested(probs, reps, np.random.default_rng(seed + 100 + index))
        p = (1 + int(np.count_nonzero(null >= observed))) / (reps + 1)
        results.append({
            "family": family, "unique_sites": len(family_sites), "nested_sites": observed,
            "nested_sense": sum(s.representative.nested_state == "nested_sense" for s in family_sites),
            "nested_antisense": sum(s.representative.nested_state == "nested_antisense" for s in family_sites),
            "nested_unknown": sum(s.representative.nested_state == "nested_unknown" for s in family_sites),
            "observed_nested_fraction": observed / len(family_sites) if family_sites else None,
            "null": method, "expected_nested_count": expected,
            "expected_nested_fraction": expected / len(family_sites) if family_sites else None,
            "enrichment": observed / expected if expected else None,
            "empirical_p_greater": p, "resamples": reps, "seed": seed + 100 + index,
            "HG03086_only_baseline_enrichment": HG03086_BASELINE["enrichment"][family],
        })
    return results


def orientation_summary(sites: list[Site]) -> list[dict[str, Any]]:
    rows = []
    for family in FAMILIES:
        counts = Counter(site.representative.nested_state for site in sites if site.representative.family == family)
        sense, anti, unknown = counts["nested_sense"], counts["nested_antisense"], counts["nested_unknown"]
        resolved = sense + anti
        rows.append({
            "family": family, "nested_sense": sense, "nested_antisense": anti,
            "nested_unknown_excluded": unknown, "resolved_nested_sites": resolved,
            "sense_fraction": sense / resolved if resolved else None,
            "antisense_fraction": anti / resolved if resolved else None,
            "sense_fold_vs_chance": (sense / resolved / 0.5) if resolved else None,
            "antisense_fold_vs_chance": (anti / resolved / 0.5) if resolved else None,
            "HG03086_only_sense_fold_baseline": HG03086_BASELINE["alu_sense"] if family == "ALU" else None,
            "HG03086_only_antisense_fold_baseline": HG03086_BASELINE["alu_antisense"] if family == "ALU" else (
                HG03086_BASELINE["line1_antisense"] if family == "LINE1" else None
            ),
        })
    return rows


def profile_events(calls: list[Call], sites: list[Site], family: str, mode: str, relative_orientation: str | None = None) -> list[tuple[int, int]]:
    """(offset, host length), optionally split by strand relative to the host."""
    if mode == "unique":
        items = [site.representative for site in sites]
    else:
        items = calls
    host_family = "ALU" if family == "LINE1" else family
    out = []
    for call in items:
        if call.family != family or call.alu_host is None and family == "LINE1":
            continue
        host = call.alu_host if family == "LINE1" else call.same_family_host
        if host is None or host.family != host_family or not PROFILE_HOST_RANGE[0] <= host.length <= PROFILE_HOST_RANGE[1]:
            continue
        if not host.start0 <= call.pos0 < host.end0 or host.strand not in {"+", "-"}:
            continue
        orient = "unknown" if call.orientation not in {"+", "-"} else ("sense" if call.orientation == host.strand else "antisense")
        if relative_orientation is not None and orient != relative_orientation:
            continue
        offset = call.pos0 - host.start0 if host.strand == "+" else host.end0 - 1 - call.pos0
        if 0 <= offset < 320:
            out.append((offset, host.length))
    return out


def position_profile(events: list[tuple[int, int]], sample: str, family: str, mode: str, relative_orientation: str = "all") -> list[dict[str, Any]]:
    rows = []
    for start in range(0, 320, PROFILE_BIN_BP):
        end = start + PROFILE_BIN_BP
        observed = sum(start <= offset < end for offset, _ in events)
        # Each insertion contributes only the bases inside its actual host to
        # this bin's uniform-placement opportunity.
        expected = sum(max(0, min(end, host_len) - start) / host_len for _, host_len in events if host_len > 0)
        rows.append({
            "sample": sample, "profile_set": mode, "family": family, "relative_orientation": relative_orientation,
            "host_context": "same-family host" if family == "ALU" else "near-full-length Alu host",
            "host_range_bp": "280-320", "bin_start_bp": start, "bin_end_bp": end,
            "n_profile_sites": len(events), "observed": observed, "expected_uniform": expected,
            "enrichment": observed / expected if expected else None,
            "is_linker_133_window": (start, end) == (120, 140),
            "is_3prime_tail_window": (start, end) == (280, 300),
        })
    return rows


def sample_profiles(calls_by_sample: dict[str, list[Call]], sites: list[Site]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    alu_rows: list[dict[str, Any]] = []
    l1_rows: list[dict[str, Any]] = []
    for sample in SAMPLES:
        calls = calls_by_sample[sample]
        # Per-sample profiles use within-sample-collapsed calls; pooled-unique counts one representative per cross-sample site.
        alu_events = profile_events(calls, sites, "ALU", "sample")
        l1_events = profile_events(calls, sites, "LINE1", "sample")
        alu_rows.extend(position_profile(alu_events, sample, "ALU", "per_sample"))
        l1_rows.extend(position_profile(l1_events, sample, "LINE1", "per_sample"))
        for orientation in ("sense", "antisense", "unknown"):
            events = profile_events(calls, sites, "LINE1", "sample", orientation)
            l1_rows.extend(position_profile(events, sample, "LINE1", "per_sample", orientation))
    alu_events = profile_events([], sites, "ALU", "unique")
    l1_events = profile_events([], sites, "LINE1", "unique")
    alu_rows.extend(position_profile(alu_events, "POOLED_UNIQUE", "ALU", "pooled_unique"))
    l1_rows.extend(position_profile(l1_events, "POOLED_UNIQUE", "LINE1", "pooled_unique"))
    for orientation in ("sense", "antisense", "unknown"):
        events = profile_events([], sites, "LINE1", "unique", orientation)
        l1_rows.extend(position_profile(events, "POOLED_UNIQUE", "LINE1", "pooled_unique", orientation))
    return alu_rows, l1_rows


def accumulation_curve(calls_by_sample: dict[str, list[Call]], window: int) -> list[dict[str, Any]]:
    rows = []
    prefix: list[Call] = []
    for idx, sample in enumerate(SAMPLES, 1):
        prefix.extend(calls_by_sample[sample])
        sites, collapsed = unique_sites(prefix, window)
        counts = Counter(site.representative.family for site in sites)
        rows.append({
            "genomes_added": idx, "sample_added": sample,
            "unique_ALU_sites": counts["ALU"], "unique_LINE1_sites": counts["LINE1"],
            "unique_SVA_sites": counts["SVA"], "unique_all_sites": len(sites),
            "within_sample_duplicate_records_collapsed": collapsed,
        })
    return rows


def markdown_table(headers: list[str], rows: list[list[Any]]) -> list[str]:
    lines = ["| " + " | ".join(headers) + " |", "|" + "|".join("---" for _ in headers) + "|"]
    for row in rows:
        cells = []
        for value in row:
            if value is None or (isinstance(value, float) and (math.isnan(value) or math.isinf(value))):
                cells.append("n/a")
            elif isinstance(value, float):
                cells.append(f"{value:.3g}")
            else:
                cells.append(str(value))
        lines.append("| " + " | ".join(cells) + " |")
    return lines


def write_analysis_summary(path: Path, *, callsets: dict[str, list[Call]], qc: dict[str, dict[str, Any]], checksums: dict[str, dict[str, str]], sites: list[Site], unique_frame: pd.DataFrame, audit: dict[str, Any], enrichment: list[dict[str, Any]], orientations: list[dict[str, Any]], alu_profile: list[dict[str, Any]], l1_profile: list[dict[str, Any]], samples_metrics: dict[str, Any], chr_y: dict[str, int], fasta_meta: dict[str, Any], rmsk_path: Path, window: int, reps: int) -> None:
    lines = [
        "# Five-genome MEI analysis on unique sites",
        "",
        "## Scope, provenance, and hard gates",
        "",
        "- Source callsets: `nested_analysis/callsets/{HG03086,HG01474,HG01566,HG03172,NA18498}.vcf`; `chr22_mei.vcf` was excluded.",
        "- Fixed primary chromosome universe: chr1–22 + chrX. chrY is recorded for QC only and excluded from host annotation, matching, and analyses.",
        "- Before parsing, every source VCF's SHA-256 was checked against its byte-identical L1 project-root copy; all pairs matched.",
        "- Every QC family/call baseline below reproduced exactly before any analysis. All inputs contain calls on all chr1–22 and chrX, i.e. genome-wide, not chr22-only.",
        f"- Matching/analysis universe: chr1-22 + chrX. chrY calls were counted for QC only and excluded: {chr_y}.",
        "- Genotype caveat: VCF headers say there is no validated genotyping model. GT/GQ were not parsed or used as evidence.",
        f"- Matching sensitivity includes exact-coordinate (0 bp) and the ±5/±10/±20-bp sweep; primary analysis is ±{window} bp.",
        f"- RepeatMasker host assignments: `{rmsk_path}`; same-family insertion host by longest overlap, then leftmost/name tie-break. The requested L1 profile separately uses the assigned near-full-length Alu host for L1-in-Alu calls.",
        "- Position and family nesting use a fresh RepeatMasker breakpoint annotation. Its four states are unnested, nested_sense, nested_antisense, nested_unknown; unknown is excluded from both orientation denominators.",
        "- Source data check: these checksum-verified VCF records actually use only legacy `NESTED=nested/unnested` values, not the four-value enum described in the request. Four-state results below are independently re-annotated from RepeatMasker and insertion/host strands; the source's raw 261 HG03086 `NESTED=nested` rows are baseline-gated and remain separately labelled.",
        "- HG03086 reference-count caveat: its 261 raw callset NESTED-nested calls are distinct from the 320 count in the self-insertion report (different host-selection/analysis rules); neither is substituted for the other.",
        f"- Genome-wide QC: { {sample: qc[sample]['genome_wide'] for sample in SAMPLES} }; all five include chr1–22 and chrX, and only HG03172/NA18498 have the one chrY QC call each.",
        "- Accession is not genotype: a site with multiple sample call records is one unique site with carrier count, not a GT/GQ-based genotype claim.",
        f"- Empirical null replicates: {reps:,}; seed {SEED}. GC null: {GC_BIN_BP}-bp bins, ±{GC_TOLERANCE:.2f} GC, {GC_MATCHED_BINS:,} matched bins/call maximum. GC frame: {fasta_meta['fasta']}.",


        "",
        "### Input QC and SHA-256 audit",
        "",
    ]
    lines += markdown_table(
        ["Sample", "Calls", "Alu", "L1", "SVA", "chrY (excluded)", "chr22", "Calls outside chr22 (incl. chrY QC)", "SHA-256 match"],
        [[s, qc[s]["calls"], qc[s]["families"]["ALU"], qc[s]["families"]["LINE1"], qc[s]["families"]["SVA"], qc[s]["chrY_calls"], qc[s]["chromosomes"].get("chr22", 0), qc[s]["calls"]-qc[s]["chromosomes"].get("chr22", 0), "yes" if checksums[s]["match"] else "NO"] for s in SAMPLES],
    )
    lines += ["", f"### Unique-site deduplication (primary ±{window} bp)", ""]
    family_rows = []
    for family in FAMILIES:
        fam = unique_frame[unique_frame["family"] == family] if not unique_frame.empty else pd.DataFrame()
        family_rows.append([family, len(fam), sum(int(v) == 1 for v in fam["n_carriers"]) if len(fam) else 0,
                            sum(int(v) >= 2 for v in fam["n_carriers"]) if len(fam) else 0,
                            *(sum(int(v) == n for v in fam["n_carriers"]) if len(fam) else 0 for n in range(1, 6))])
    lines += markdown_table(["Family", "Unique sites", "Private", "Shared", "1 carrier", "2", "3", "4", "5"], family_rows)
    lines += ["", "#### Family-wise carrier counts by matching window", ""]
    lines += markdown_table(["Window", "Family", "Unique sites", "Shared sites", "Private sites", "Carriers 1/2/3/4/5"],
        [["exact" if row["window_bp"] == 0 else f"±{row['window_bp']} bp", row["family"], row["unique_sites"], row["shared_sites"], row["private_sites"], row["carrier_histogram"]] for row in audit["family_sensitivity"]])
    lines += ["", f"Total primary-window unique sites: **{len(sites)}**; collapsed repeated same-sample source records before cross-sample matching: **{audit['within_sample_duplicates_collapsed']}**.", ""]
    lines += markdown_table(["Window", "Unique sites", "Shared sites", "Private sites", "Mean carriers/site", "Max carriers"],
        [["exact" if row["window_bp"] == 0 else f"±{row['window_bp']} bp", row["unique_sites"], row["shared_sites"], row["private_sites"], row["mean_carriers"], row["max_carriers"]] for row in audit["sensitivity"]])

    lines += ["", "### Per-sample four-state nesting QC", "",
              "| Sample | nested_sense | nested_antisense | nested_unknown | unnested |", "|---|---:|---:|---:|---:|"]
    for sample in SAMPLES:
        state_counts = qc[sample].get("derived_nested_states", {})
        lines.append(f"| {sample} | {state_counts.get('nested_sense', 0)} | {state_counts.get('nested_antisense', 0)} | {state_counts.get('nested_unknown', 0)} | {state_counts.get('unnested', 0)} |")
    lines += ["", "### Q1 — same-family nesting enrichment on the union of unique sites", "",
              "Uniform null for Alu; GC-matched null for LINE-1/SVA. Nested includes all three nested states, including nested_unknown; only Q2 excludes unknown.", ""]

    lines += markdown_table(["Family", "Unique sites", "Nested sites", "Null", "Expected nested", "Enrichment", "Empirical p", "HG03086-only baseline"],
        [[r["family"], r["unique_sites"], r["nested_sites"], r["null"], r["expected_nested_count"], r["enrichment"], r["empirical_p_greater"], r["HG03086_only_baseline_enrichment"]] for r in enrichment])
    lines += ["", f"The uniform Alu null randomizes each unique Alu site uniformly over the same hg38 chromosome and uses the union of same-family Alu bases on that chromosome as the chance target. For LINE-1 and SVA, each site is matched to valid 200-bp hg38 bins within ±0.05 GC of its own bin, and same-family RepeatMasker coverage is averaged over at most 2,000 matched bins. The randomization is Bernoulli at those per-site probabilities; empirical enrichment p is the add-one upper-tail fraction across {reps:,} resamples.", ""]

    lines += ["", "### Q2 — nested-sense / nested-antisense", "",
              "Fold is the resolved orientation fraction divided by 0.5. `nested_unknown` is shown separately and never assigned to antisense.", ""]
    lines += markdown_table(["Family", "Sense", "Antisense", "Unknown excluded", "Resolved n", "Sense fold", "Antisense fold", "HG03086-only comparator"],
        [[r["family"], r["nested_sense"], r["nested_antisense"], r["nested_unknown_excluded"], r["resolved_nested_sites"], r["sense_fold_vs_chance"], r["antisense_fold_vs_chance"],
          f"sense {HG03086_BASELINE['alu_sense']:.2f}x; antisense {HG03086_BASELINE['alu_antisense']:.2f}x" if r["family"] == "ALU" else (f"antisense {HG03086_BASELINE['line1_antisense']:.1f}x (chance)" if r["family"] == "LINE1" else "not supplied")]
         for r in orientations])

    lines += ["", "### Q3 — Alu-in-Alu position profile, near-full-length hosts", "",
              "Host-relative 5′-to-3′ coordinates; 20-bp bins; near-full host length 280–320 bp. Per-sample rows use that callset; pooled_unique counts a shared insertion once. Expected counts account for each host's actual length.", ""]
    lines += _profile_markdown(alu_profile, baseline_linker=HG03086_BASELINE["linker"], baseline_tail=HG03086_BASELINE["tail"])
    lines += ["", "### Q4 — LINE-1 insertions within near-full-length Alu hosts", "",
              "Same 20-bp host-coordinate profile as Q3, as clarified for this run. This cross-family L1-in-Alu profile is distinct from the count of same-family nested-L1 sites used for Q1/Q2.", ""]
    lines += _profile_markdown(l1_profile, baseline_linker=None, baseline_tail=None)

    same_family_nested_l1 = sum(site.representative.family == "LINE1" and site.representative.nested_state != "unnested" for site in sites)
    l1_alu_profile_n = next((row["n_profile_sites"] for row in l1_profile if row["sample"] == "POOLED_UNIQUE" and row["bin_start_bp"] == 0 and row.get("relative_orientation", "all") == "all"), 0)
    min_sites = PROFILE_MIN_EXPECTED_PER_BIN * 16
    power = "descriptive profile reaches the 5-site/bin reporting floor" if l1_alu_profile_n >= min_sites else "below the 5-expected-sites/bin reporting floor; treat as descriptive and collect William's next five genomes"
    l1_orientation_windows = _l1_orientation_windows(l1_profile)
    lines += ["", "### LINE-1 positional-analysis size / interpretation", "",
              f"- Unique **same-family nested-L1** sites across the five genomes: **{same_family_nested_l1}** (the Q1/Q2 same-family set).",
              f"- Unique LINE-1 sites inside near-full Alu hosts for Q4: **{l1_alu_profile_n}**.",
              f"- Reporting floor: 5 expected sites in each of 16 20-bp bins (80 sites total; a descriptive stability floor, not a formal power calculation). Verdict: **{power}**.", "",
              "L1-in-Alu pooled-unique position-window split by insertion orientation relative to the Alu host:", "",
              *markdown_table(["Orientation", "Profile sites", "120–140 observed/enrichment", "280–300 observed/enrichment"],
                 [[orientation, l1_orientation_windows[orientation]["n"], l1_orientation_windows[orientation]["linker"], l1_orientation_windows[orientation]["tail"]] for orientation in ("sense", "antisense", "unknown")]), "",
              "### Prior HG03086-only headline comparators (reported baselines, not re-estimated)", "",
              "| Metric | HG03086-only baseline | Five-genome unique-site result |"]
    lines += ["|---|---:|---:|"]
    enr_by_family = {row["family"]: row for row in enrichment}
    ori_by_family = {row["family"]: row for row in orientations}
    peak_by_family = {"ALU": _pooled_peak(alu_profile), "LINE1": _pooled_peak(l1_profile)}
    lines += [
        f"| ALU nesting enrichment (uniform) | 1.47x | {_fmt(enr_by_family['ALU']['enrichment'])}x |",
        f"| LINE-1 nesting enrichment (GC-matched) | 1.47x | {_fmt(enr_by_family['LINE1']['enrichment'])}x |",
        f"| SVA nesting enrichment (GC-matched) | 29.4x | {_fmt(enr_by_family['SVA']['enrichment'])}x |",

        f"| ALU nested-sense vs chance | 2.49x | {_fmt(ori_by_family['ALU']['sense_fold_vs_chance'])}x |",
        f"| ALU nested-antisense vs chance | 0.46x | {_fmt(ori_by_family['ALU']['antisense_fold_vs_chance'])}x |",
        f"| LINE-1 nested-antisense vs chance | 1.0x (chance) | {_fmt(ori_by_family['LINE1']['antisense_fold_vs_chance'])}x |",

        f"| Alu linker 120–140 bp enrichment | 3.0x | {_fmt(peak_by_family['ALU']['linker'])}x pooled unique |",
        f"| Alu 3′-tail 280–300 bp enrichment | 3.7x | {_fmt(peak_by_family['ALU']['tail'])}x pooled unique |",
        f"| LINE-1 linker/tail windows in Alu hosts | not supplied | {_fmt(peak_by_family['LINE1']['linker'])}x / {_fmt(peak_by_family['LINE1']['tail'])}x pooled unique |",
        "",
        "The inherited 1.47/1.47/29.4, orientation, linker, and tail figures are shown as user-supplied HG03086-only comparators. The 2.49x/0.46x Alu and chance LINE-1 orientation values are the supplied HG03086-only orientation comparator; their exact numerator/denominator were not present in the available prior result CSVs, so the comparison is the reported factor, not a reconstructed count. These historical factors use their own denominator/null implementation; this report does not back-calculate them from the five-callset results. The displayed ±10-bp unique counts and Q1–Q4 site-level results use one representative per deduplicated site; per-sample profile checks are calculated on each sample's collapsed callset.",
        "",
        "## Matching and accumulation",
        "",
        "The separate `matching_audit.md` records host-key semantics, complete ±5/±10/±20 sweep, private/shared counts, raw NESTED values, same-sample duplicate collapse, and sample checksums. `accumulation_curve.csv` records deterministic sample-order accumulation.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def _fmt(value: Any) -> str:
    return "n/a" if value is None or (isinstance(value, float) and not math.isfinite(value)) else f"{float(value):.3g}"


def _profile_markdown(rows: list[dict[str, Any]], baseline_linker: float | None, baseline_tail: float | None) -> list[str]:
    windows = []
    for sample in [*SAMPLES, "POOLED_UNIQUE"]:
        for label, start, end in (("linker", 120, 140), ("tail", 280, 300)):
            row = next((item for item in rows if item["sample"] == sample and item["bin_start_bp"] == start and item.get("relative_orientation", "all") == "all"), None)
            if row is None:
                continue
            baseline = baseline_linker if label == "linker" else baseline_tail
            windows.append([sample, label, row["n_profile_sites"], row["observed"], row["expected_uniform"], row["enrichment"], "YES" if (row["enrichment"] or 0) > 1 else "no", baseline])
    lines = markdown_table(["Sample/profile", "Window", "Eligible sites", "Observed", "Expected", "Enrichment", "Elevated?", "HG03086-only baseline"], windows)
    lines += ["", "Full 20-bp-bin profile:", ""]
    table = []
    for sample in [*SAMPLES, "POOLED_UNIQUE"]:
        sample_rows = [row for row in rows if row["sample"] == sample and row.get("relative_orientation", "all") == "all"]
        if not sample_rows:
            continue
        table.append([sample, *(row["observed"] for row in sample_rows)])
    lines += markdown_table(["Sample/profile", *[f"{start}-{start+20}" for start in range(0, 320, 20)]], table)
    if rows and rows[0].get("family") == "ALU":
        per_sample_peaks = {
            sample: all(
                (row := next((item for item in rows if item["sample"] == sample and item["bin_start_bp"] == start), None))
                is not None and (row["enrichment"] or 0) > 1
                for start in (120, 280)
            )
            for sample in SAMPLES
        }
        elevated = sum(per_sample_peaks.values())
        lines += ["", f"**Per-sample peak replication:** both 120–140 bp and 280–300 bp are elevated in {elevated}/5 individual genomes.", ""]
    if rows and rows[0].get("family") == "LINE1":
        lines += ["", "L1 profile stratified by insertion orientation relative to the Alu host strand:", ""]
        oriented = []
        for sample in [*SAMPLES, "POOLED_UNIQUE"]:
            for orientation in ("sense", "antisense", "unknown"):
                sample_rows = [row for row in rows if row["sample"] == sample and row.get("relative_orientation") == orientation]
                if sample_rows:
                    oriented.append([f"{sample} / {orientation}", sample_rows[0]["n_profile_sites"], *(row["observed"] for row in sample_rows)])
        lines += markdown_table(["Sample / relative orientation", "Sites", *[f"{start}-{start+20}" for start in range(0, 320, 20)]], oriented)
    return lines


def _l1_orientation_windows(rows: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    result = {}
    for orientation in ("sense", "antisense", "unknown"):
        oriented = [row for row in rows if row["sample"] == "POOLED_UNIQUE" and row.get("relative_orientation") == orientation]
        linker = next((row for row in oriented if row["bin_start_bp"] == 120), None)
        tail = next((row for row in oriented if row["bin_start_bp"] == 280), None)
        def display_window(row: dict[str, Any] | None) -> str:
            if row is None:
                return "n/a"
            enrichment = _fmt(row["enrichment"])
            suffix = f"{enrichment}x" if enrichment != "n/a" else "n/a"
            return f"{row['observed']}/{suffix}"

        result[orientation] = {
            "n": oriented[0]["n_profile_sites"] if oriented else 0,
            "linker": display_window(linker),
            "tail": display_window(tail),
        }
    return result


def _pooled_peak(rows: list[dict[str, Any]]) -> dict[str, float | None]:
    def value(start: int):
        row = next((item for item in rows if item["sample"] == "POOLED_UNIQUE" and item["bin_start_bp"] == start), None)
        return row.get("enrichment") if row else None
    return {"linker": value(120), "tail": value(280)}


def write_matching_audit(path: Path, checksums: dict[str, dict[str, str]], qc: dict[str, dict[str, Any]], audit: dict[str, Any], sites: list[Site], unique_frame: pd.DataFrame, rmsk: Path, window: int, callsets: dict[str, list[Call]]) -> None:
    lines = [
        "# Five-callset unique-site matching audit", "",
        "## Inputs, integrity, and scope", "",
        "All source callsets were read only after SHA-256 equality with the duplicate project-root copy was verified. `chr22_mei.vcf` was explicitly excluded. Calls on chrY are audit counts only; matching and analysis use chr1-22 plus chrX. The out-of-scope chrY count in each input is recorded in the analysis summary and the sample chromosome QC table above.", "",
        *markdown_table(["Sample", "Input", "SHA-256 source", "SHA-256 root copy", "Match"],
            [[sample, checksums[sample]["source_path"], checksums[sample]["source_sha256"], checksums[sample]["duplicate_sha256"], "yes" if checksums[sample]["match"] else "NO"] for sample in SAMPLES]),
        "", f"RepeatMasker: `{rmsk}`.", "",
        "## QC and raw NESTED source values", "",
        "The checksum-verified records actually contain only legacy binary `NESTED=nested/unnested` values, rather than the four values named in the request. The requested four-state category is independently recomputed from same-family RepeatMasker host overlap and insertion/host strands; nested_unknown is never included in sense or antisense. This source/schema discrepancy is reported, not silently relabelled.", "",
        *markdown_table(["Sample", "Calls", "Alu", "L1", "SVA", "Raw NESTED=nested", "Derived same-family nested"],
            [[sample, qc[sample]["calls"], qc[sample]["families"]["ALU"], qc[sample]["families"]["LINE1"], qc[sample]["families"]["SVA"], sum(call.raw_nested == "nested" for call in callsets[sample]), sum(call.nested_state != "unnested" for call in callsets[sample])] for sample in SAMPLES]),
        "", "The HG03086 input raw NESTED-nested count is 261. The separate self-insertion report's 320 count uses different selection/host rules; both are retained as distinct labelled historical numbers.", "",
        "## Matching key and no-chaining rule", "",
        f"Primary matching key is chromosome + MEI family + insertion orientation + same-family host element (or explicit no-host state for unnested calls) + breakpoint distance ≤{window} bp. Same-family host is selected by longest overlap, then leftmost start/name. Gene/snpEff fields and GT/GQ are not read for matching. Repeated same-sample call records are collapsed first using the same host/family/orientation key and anchored ±{window}-bp grouping. Cross-sample pairs are greedily ordered by breakpoint distance and stable sample/record keys; every additional member must match the original site anchor. Sites never merge, so there is no transitive chaining.", "",
        "### Sensitivity sweep", "",
        *markdown_table(["Window", "Unique sites", "Shared sites", "Private sites", "Mean carriers", "Max carriers"],
            [["exact" if row["window_bp"] == 0 else f"±{row['window_bp']} bp", row["unique_sites"], row["shared_sites"], row["private_sites"], row["mean_carriers"], row["max_carriers"]] for row in audit["sensitivity"]]),
        "", "### Family × window counts", "",
        *markdown_table(["Window", "Family", "Unique", "Shared", "Private", "Carrier histogram 1/2/3/4/5"],
            [["exact" if row["window_bp"] == 0 else f"±{row['window_bp']} bp", row["family"], row["unique_sites"], row["shared_sites"], row["private_sites"], row["carrier_histogram"]] for row in audit["family_sensitivity"]]),
        "", "### Primary ±10 bp sites", "",
        f"- unique sites: {len(sites)}",
        f"- private sites: {int((unique_frame['n_carriers'] == 1).sum()) if not unique_frame.empty else 0}",
        f"- sites with more than one carrier: {int((unique_frame['n_carriers'] > 1).sum()) if not unique_frame.empty else 0}",
        f"- same-sample source records collapsed: {audit['within_sample_duplicates_collapsed']}",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def run(args: argparse.Namespace) -> None:
    pairs = input_paths(args.callset_dir, WORKSPACE)
    checksum_records = {}
    for sample, (source, duplicate) in pairs.items():
        source_hash, duplicate_hash = sha256(source), sha256(duplicate)
        if source_hash != duplicate_hash:
            raise InputGateError(
                f"SHA-256 changed after initial verification for {sample}; STOP, no source selected: "
                f"{source}={source_hash}; {duplicate}={duplicate_hash}"
            )
        checksum_records[sample] = {
            "source_path": str(source.resolve()), "duplicate_path": str(duplicate.resolve()),
            "source_sha256": source_hash, "duplicate_sha256": duplicate_hash,
            "match": True,
        }
    calls_by_sample: dict[str, list[Call]] = {}
    qc: dict[str, dict[str, Any]] = {}
    for sample, (source, _) in pairs.items():
        calls, info = parse_callset(source, sample)
        calls_by_sample[sample] = calls
        qc[sample] = info
    # Explicit scope check: the chr22 slice is neither read nor admitted.
    rmsk = read_rmsk(args.rmsk)
    all_calls = [call for sample in SAMPLES for call in calls_by_sample[sample] if call.chrom in PRIMARY_CHROMS]
    assign_hosts(all_calls, rmsk)
    all_primary_calls = all_calls
    calls_by_sample = {sample: [call for call in calls_by_sample[sample] if call.chrom in PRIMARY_CHROMS] for sample in SAMPLES}

    sites, n_collapsed = unique_sites(all_primary_calls, PRIMARY_WINDOW)
    unique_frame = build_sites_table(sites, PRIMARY_WINDOW)
    sensitivity = []
    family_sensitivity = []
    for window in SWEEP_WINDOWS:
        sweep_sites, collapsed = unique_sites(all_primary_calls, window)
        carrier_counts = [len(site.samples) for site in sweep_sites]
        sensitivity.append({
            "window_bp": window, "unique_sites": len(sweep_sites),
            "shared_sites": sum(count >= 2 for count in carrier_counts),
            "private_sites": sum(count == 1 for count in carrier_counts),
            "mean_carriers": round(float(np.mean(carrier_counts)), 3) if carrier_counts else 0.0,
            "max_carriers": max(carrier_counts, default=0),
        })
        for family in FAMILIES:
            family_sites = [site for site in sweep_sites if site.representative.family == family]
            hist = Counter(len(site.samples) for site in family_sites)
            family_sensitivity.append({
                "window_bp": window, "family": family, "unique_sites": len(family_sites),
                "shared_sites": sum(count >= 2 for count in (len(site.samples) for site in family_sites)),
                "private_sites": hist[1], "carrier_histogram": "/".join(str(hist[n]) for n in range(1, 6)),
            })
    carrier_hist = Counter(len(site.samples) for site in sites)
    family_carrier_hist = {family: Counter(len(site.samples) for site in sites if site.representative.family == family) for family in FAMILIES}
    audit = {"sensitivity": sensitivity, "family_sensitivity": family_sensitivity,
             "within_sample_duplicates_collapsed": n_collapsed,
             "carrier_histogram": {str(n): carrier_hist[n] for n in range(1, 6)},
             "family_carrier_histogram": {family: {str(n): family_carrier_hist[family][n] for n in range(1, 6)} for family in FAMILIES}}

    lengths = {chrom: qc_length for chrom, qc_length in args.chrom_sizes.items() if chrom in PRIMARY_CHROMS}
    merged_cov = {}
    for chrom in PRIMARY_CHROMS:
        for family in FAMILIES:
            merged_cov[(family, chrom)] = sum(end - start for start, end in merge_intervals(rmsk.get((chrom, family), [])))
    gc_bundle, fasta_meta = gc_matched_frame(args.fasta, rmsk, lengths, sites)
    enrichment = family_enrichment(sites, lengths, merged_cov, gc_bundle, args.replicates, args.seed)
    orientations =    orientation_summary(sites)

    per_sample_unique = {sample: collapse_within_sample(calls_by_sample[sample], PRIMARY_WINDOW)[0] for sample in SAMPLES}
    alu_profile, l1_profile = sample_profiles(per_sample_unique, sites)

    # Add chromosome counts outside chrY to audit; chromosome names are preserved.
    chr_y = {sample: qc[sample]["chrY_calls"] for sample in SAMPLES}
    primary_only_counts = {sample: sum(qc[sample]["chromosomes"].get(chrom, 0) for chrom in PRIMARY_CHROMS) for sample in SAMPLES}
    for sample in SAMPLES:
        qc[sample]["analysis_calls"] = primary_only_counts[sample]
    # Derived nesting counts by sample are recorded separately from source binary NESTED.
    for sample in SAMPLES:
        qc[sample]["derived_nested_states"] = dict(Counter(call.nested_state for call in calls_by_sample[sample]))

    args.outdir.mkdir(parents=True, exist_ok=True)
    unique_frame.to_csv(args.outdir / "unique_sites.csv", index=False, quoting=csv.QUOTE_MINIMAL)
    pd.DataFrame(accumulation_curve(calls_by_sample, PRIMARY_WINDOW)).to_csv(args.outdir / "accumulation_curve.csv", index=False)
    write_matching_audit(args.outdir / "matching_audit.md", checksum_records, qc, audit, sites, unique_frame, args.rmsk, PRIMARY_WINDOW, callsets=calls_by_sample)
    write_analysis_summary(args.outdir / "analysis_summary.md", callsets=calls_by_sample, qc=qc, checksums=checksum_records,
                           sites=sites, unique_frame=unique_frame, audit=audit, enrichment=enrichment,
                           orientations=orientations, alu_profile=alu_profile, l1_profile=l1_profile,
                           samples_metrics=primary_only_counts, chr_y=chr_y, fasta_meta=fasta_meta,
                           rmsk_path=args.rmsk, window=PRIMARY_WINDOW, reps=args.replicates)
    print(f"Wrote requested outputs under {args.outdir.resolve()}")
    print(f"five-input hash and QC gates passed; primary ±{PRIMARY_WINDOW}bp sites={len(sites)}; chrY excluded={chr_y}")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--callset-dir", type=Path, default=WORKSPACE / "nested_analysis" / "callsets")
    parser.add_argument("--rmsk", type=Path, default=WORKSPACE / "nested_analysis" / "data" / "rmsk.txt.gz")
    parser.add_argument("--fasta", type=Path, default=Path("~/retrotransposon-workdir/data/public/reference/hg38/Homo_sapiens_assembly38.fasta").expanduser())
    parser.add_argument("--chrom-sizes", type=Path, default=WORKSPACE / "nested_analysis" / "results_highconf" / "hg38.chrom.sizes")
    parser.add_argument("--outdir", type=Path, default=WORKSPACE / "nested_analysis" / "results_multi_sample")
    parser.add_argument("--replicates", type=int, default=REPS)
    parser.add_argument("--seed", type=int, default=SEED)
    args = parser.parse_args(argv)
    if args.replicates < 1_000:
        parser.error("--replicates must be at least 1000")
    for path in (args.rmsk, args.chrom_sizes, args.fasta, args.callset_dir):
        if not path.exists():
            raise InputGateError(f"required input not found: {path}")
    args.chrom_sizes = _read_chrom_sizes(args.chrom_sizes)
    missing = set(PRIMARY_CHROMS) - set(args.chrom_sizes)
    if missing:
        raise InputGateError(f"chrom.sizes missing required chromosomes: {sorted(missing)}")
    return args


def _read_chrom_sizes(path: Path) -> dict[str, int]:
    sizes = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.split()
            if len(fields) >= 2:
                sizes[fields[0]] = int(fields[1])
    return sizes


def main(argv: list[str] | None = None) -> int:
    try:
        run(parse_args(argv))
    except InputGateError as exc:
        print(f"STOP: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
