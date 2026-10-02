#!/usr/bin/env python3
"""Ten-genome multi-sample MEI analysis: dedup, layers, enrichment, profiles.

William Brandler's standing instruction is "look at only unique variants not
shared between samples", so Layer 1 (private sites, carrier count == 1) leads
the report.  Layer 2 is the dedup union and Layer 3 the shared set.

Scientific rules enforced here (see the task specification):

* GT/GQ are never parsed -- the VCF headers state there is no validated
  genotyping model, so they are not evidence.
* chrY is counted for QC but excluded from matching and analysis.
* The binary VCF NESTED field is never trusted for the four-state enum; it is
  recomputed from RepeatMasker host overlap plus orientation, and ``unknown`` is
  never folded into ``antisense``.
* Gene / snpEff annotation fields are never used as matching keys.
* No transitive chaining: every additional member must match the site anchor.
* Exactly one call per sample per site, so one carrier contributes once.

Matching key: chromosome + family + orientation + host element + breakpoint
within +/-10 bp (primary), with an exact / +/-5 / +/-20 sweep also reported.
"""
from __future__ import annotations

import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

# Shared plumbing comes from the project's existing modules: `dedup_samples`
# supplies the host model, RepeatMasker parsing and .fai FASTA access, while
# `mei_reference_opportunity` supplies the reference-opportunity and consensus
# machinery.  Neither is duplicated here.
import dedup_samples as engine  # noqa: E402
import mei_reference_opportunity as refopp  # noqa: E402

SAMPLES: tuple[str, ...] = (
    "HG03086",
    "HG01474",
    "HG01566",
    "HG03172",
    "NA18498",
    "HG00171",
    "HG01058",
    "NA18939",
    "NA19017",
    "NA20845",
)
SAMPLE_INDEX = {s: i for i, s in enumerate(SAMPLES)}
ORIGINAL_FIVE: tuple[str, ...] = SAMPLES[:5]

ALLOWED_CHROMS = frozenset((*engine.PRIMARY_CHROMS, "chrY"))
PRIMARY_WINDOW = 10
SWEEP_WINDOWS = (0, 5, 10, 20)
GC_BIN_BP = 200
GC_TOLERANCE = 0.05
GC_MATCHED_BINS = 2_000
REPS = 10_000
SEED = 20261001
ALU_HOST_RANGE = (280, 320)
PROFILE_BIN_BP = 20
LINKER_WINDOW = (120, 140)
TAIL_WINDOW = (280, 300)
L1_BIN_BP = 500
MIN_EXPECTED_PER_BIN = 5.0

WORKSPACE = Path("/Users/smyan/Desktop/Research/Research_Projects/L1")


class GateError(RuntimeError):
    """Raised when a verification gate fails; the run must stop there."""


@dataclass
class Call:
    sample: str
    source_index: int
    chrom: str
    pos: int
    family: str
    orientation: str
    raw_nested: str
    info: dict[str, str] = field(default_factory=dict)
    same_family_host: engine.Host | None = None
    match_host: engine.Host | None = None
    nested_state: str = "unnested"
    source_count: int = 1
    source_ids: list[str] = field(default_factory=list)

    @property
    def pos0(self) -> int:
        return self.pos - 1

    @property
    def match_key(self) -> tuple[Any, ...]:
        host_key = self.match_host.key if self.match_host else None
        return (self.chrom, self.family, self.orientation or ".", host_key)

    @property
    def offset(self) -> int | None:
        host = self.match_host
        if host is None or not host.start0 <= self.pos0 < host.end0:
            return None
        return (
            self.pos0 - host.start0
            if host.strand == "+"
            else host.end0 - 1 - self.pos0
        )


@dataclass
class Site:
    anchor: Call
    members: list[Call]

    @property
    def samples(self) -> tuple[str, ...]:
        return tuple(
            sorted({c.sample for c in self.members}, key=lambda s: SAMPLE_INDEX[s])
        )

    @property
    def representative(self) -> Call:
        median = float(np.median([c.pos for c in self.members]))
        return min(
            self.members,
            key=lambda c: (
                abs(c.pos - median),
                SAMPLE_INDEX[c.sample],
                c.pos,
                c.source_index,
            ),
        )


# --------------------------------------------------------------------------
# parsing
# --------------------------------------------------------------------------
def parse_info(text: str) -> dict[str, str]:
    info: dict[str, str] = {}
    if text in {"", "."}:
        return info
    for item in text.split(";"):
        key, sep, value = item.partition("=")
        if key:
            info[key] = value if sep else "True"
    return info


def parse_callset(path: Path, sample: str) -> tuple[list[Call], dict[str, Any]]:
    """Parse one callset.  Gene/snpeff fields and GT/GQ are never read."""
    calls: list[Call] = []
    fam = Counter()
    chrom = Counter()
    raw_nested = Counter()
    errors: list[str] = []
    reference = ""
    header: list[str] | None = None
    with path.open(encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, 1):
            line = line.rstrip("\r\n")
            if line.startswith("##reference="):
                reference = line.split("=", 1)[1].strip()
                continue
            if line.startswith("#CHROM"):
                header = line.split("\t")
                continue
            if line.startswith("#") or not line:
                continue
            fields = line.split("\t")
            if header is None or len(fields) != len(header):
                errors.append(f"line {line_no}: {len(fields)} columns")
                continue
            # only CHROM/POS/ID/INFO are touched; never GENE/snpEff/GT/GQ
            info = parse_info(fields[7])
            family = engine.normalize_family(info.get("MEIFAMILY", ""))
            if family not in engine.FAMILIES:
                errors.append(f"line {line_no}: bad MEIFAMILY {info.get('MEIFAMILY')!r}")
                continue
            try:
                pos = int(fields[1])
            except ValueError:
                errors.append(f"line {line_no}: bad POS {fields[1]!r}")
                continue
            if pos < 1:
                errors.append(f"line {line_no}: POS < 1")
                continue
            calls.append(
                Call(
                    sample=sample,
                    source_index=len(calls),
                    chrom=fields[0],
                    pos=pos,
                    family=family,
                    orientation=info.get("ORIENT", "").strip(),
                    raw_nested=info.get("NESTED", "").strip().lower(),
                    info=info,
                    source_ids=[
                        fields[2]
                        if fields[2] not in {"", "."}
                        else f"{sample}:{fields[0]}:{pos}:{len(calls)}"
                    ],
                )
            )
            fam[family] += 1
            chrom[fields[0]] += 1
            raw_nested[info.get("NESTED", "").strip().lower() or "<MISSING>"] += 1
    if errors:
        raise GateError(f"{sample}: malformed VCF ({'; '.join(errors[:5])})")
    if reference and reference.lower() != "hg38":
        raise GateError(f"{sample}: reference {reference!r}, expected hg38")
    stray = set(chrom) - ALLOWED_CHROMS
    if stray:
        raise GateError(f"{sample}: out-of-scope chromosomes {sorted(stray)}")
    return calls, {
        "sample": sample,
        "calls": len(calls),
        "families": dict(fam),
        "chrY_calls": chrom.get("chrY", 0),
        "raw_nested_values": dict(raw_nested),
        "primary_only_calls": sum(
            chrom.get(c, 0) for c in engine.PRIMARY_CHROMS
        ),
    }


# --------------------------------------------------------------------------
# host assignment (longest containing MEI; tie-break left, then name)
# --------------------------------------------------------------------------
def assign_hosts(
    calls: Sequence[Call], hosts: dict[tuple[str, str], list[engine.Host]]
) -> None:
    """Assign each call the longest same-family host containing its position.

    A sweep over hosts sorted by start, with an active set pruned by end, so
    this stays near-linear instead of scanning every host for every call.
    """
    by_chrom: dict[str, list[Call]] = defaultdict(list)
    for call in calls:
        if call.chrom in engine.PRIMARY_CHROMS:
            by_chrom[call.chrom].append(call)
    for chrom, chrom_calls in by_chrom.items():
        for family in engine.FAMILIES:
            fam_calls = [c for c in chrom_calls if c.family == family]
            if not fam_calls:
                continue
            fam_hosts = sorted(
                hosts.get((chrom, family), []),
                key=lambda h: (h.start0, h.end0, h.name, h.strand),
            )
            fam_calls.sort(key=lambda c: c.pos)
            active: list[engine.Host] = []
            next_host = 0
            for call in fam_calls:
                # VCF POS is 1-based; a host contains POS when start0 < POS <= end0
                while next_host < len(fam_hosts) and fam_hosts[next_host].start0 < call.pos:
                    active.append(fam_hosts[next_host])
                    next_host += 1
                active = [h for h in active if h.end0 >= call.pos]
                if not active:
                    call.same_family_host = None
                    call.match_host = None
                    call.nested_state = "unnested"
                    continue
                best = min(
                    active, key=lambda h: (-h.length, h.start0, h.name, h.strand)
                )
                call.same_family_host = best
                call.match_host = best
                if call.orientation not in {"+", "-"} or best.strand not in {"+", "-"}:
                    call.nested_state = "nested_unknown"
                elif call.orientation == best.strand:
                    call.nested_state = "nested_sense"
                else:
                    call.nested_state = "nested_antisense"


# --------------------------------------------------------------------------
# dedup (frozen rules)
# --------------------------------------------------------------------------
def candidate_pairs(
    calls: Sequence[Call], window: int, *, include_same_sample: bool
) -> list[tuple[int, int, int]]:
    grouped: dict[tuple[Any, ...], list[int]] = defaultdict(list)
    for index, call in enumerate(calls):
        grouped[call.match_key].append(index)
    pairs: list[tuple[int, int, int]] = []
    for indexes in grouped.values():
        indexes.sort(
            key=lambda i: (calls[i].pos, SAMPLE_INDEX[calls[i].sample], calls[i].source_index)
        )
        for position, left_index in enumerate(indexes):
            left = calls[left_index]
            for right_index in indexes[position + 1 :]:
                right = calls[right_index]
                distance = right.pos - left.pos
                if distance > window:
                    break
                if include_same_sample or left.sample != right.sample:
                    pairs.append((distance, left_index, right_index))
    pairs.sort(
        key=lambda item: (
            item[0],
            calls[item[1]].chrom,
            calls[item[1]].pos,
            SAMPLE_INDEX[calls[item[1]].sample],
            calls[item[1]].source_index,
            SAMPLE_INDEX[calls[item[2]].sample],
            calls[item[2]].source_index,
        )
    )
    return pairs


def anchored_groups(
    calls: Sequence[Call], window: int, *, allow_same_sample: bool
) -> list[Site]:
    """Greedy representative-anchored grouping; sites never merge transitively."""
    pairs = candidate_pairs(calls, window, include_same_sample=allow_same_sample)
    site_of: dict[int, int] = {}
    sites: list[Site] = []
    for _, left_index, right_index in pairs:
        left_site = site_of.get(left_index)
        right_site = site_of.get(right_index)
        if left_site is not None and right_site is not None:
            continue
        if left_site is None and right_site is None:
            first, second = sorted(
                (left_index, right_index),
                key=lambda idx: (
                    calls[idx].pos,
                    SAMPLE_INDEX[calls[idx].sample],
                    calls[idx].source_index,
                ),
            )
            if not allow_same_sample and calls[first].sample == calls[second].sample:
                continue
            site_of[first] = site_of[second] = len(sites)
            sites.append(Site(calls[first], [calls[first], calls[second]]))
            continue
        site_id = left_site if left_site is not None else right_site
        new_index = right_index if left_site is not None else left_index
        site = sites[site_id]
        # no transitive chaining: every member must match the original anchor
        if abs(calls[new_index].pos - site.anchor.pos) > window:
            continue
        if not allow_same_sample and calls[new_index].sample in {
            member.sample for member in site.members
        }:
            continue
        site_of[new_index] = site_id
        site.members.append(calls[new_index])
    for index, call in enumerate(calls):
        if index not in site_of:
            site_of[index] = len(sites)
            sites.append(Site(call, [call]))
    return sorted(
        sites,
        key=lambda site: (
            site.anchor.chrom,
            site.anchor.pos,
            site.anchor.family,
            site.anchor.orientation,
            site.anchor.match_host.key if site.anchor.match_host else (),
            SAMPLE_INDEX[site.anchor.sample],
            site.anchor.source_index,
        ),
    )


def collapse_within_sample(
    calls: Sequence[Call], window: int
) -> tuple[list[Call], int]:
    """Collapse repeated same-sample records at one anchored site.

    Aggregates provenance onto a *copy*; the caller's Call objects are never
    mutated, because the same objects are reused across the window sweep and
    accumulation curve.
    """
    by_sample: dict[str, list[Call]] = defaultdict(list)
    for call in calls:
        by_sample[call.sample].append(call)
    representatives: list[Call] = []
    collapsed = 0
    for sample in SAMPLES:
        sample_calls = by_sample.get(sample)
        if not sample_calls:
            continue
        for site in anchored_groups(sample_calls, window, allow_same_sample=True):
            rep = site.representative
            source_count = sum(m.source_count for m in site.members)
            source_ids = [i for m in site.members for i in m.source_ids]
            if source_count != rep.source_count or source_ids != rep.source_ids:
                rep = replace(rep, source_count=source_count, source_ids=source_ids)
            representatives.append(rep)
            collapsed += len(site.members) - 1
    return representatives, collapsed


def unique_sites(calls: Sequence[Call], window: int) -> tuple[list[Site], int]:
    canonical, collapsed = collapse_within_sample(calls, window)
    return anchored_groups(canonical, window, allow_same_sample=False), collapsed


def site_row(site: Site, site_id: str, window: int) -> dict[str, Any]:
    call = site.representative
    samples = site.samples
    host = call.same_family_host
    return {
        "site_id": site_id,
        "chrom": call.chrom,
        "representative_pos": call.pos,
        "family": call.family,
        "insertion_orientation": call.orientation,
        "same_family_nested_state": call.nested_state,
        "raw_NESTED_value": call.raw_nested,
        "host_id": host.host_id if host else "",
        "host_name": host.name if host else "",
        "host_start0": host.start0 if host else "",
        "host_end0": host.end0 if host else "",
        "host_strand": host.strand if host else "",
        "host_len": host.length if host else "",
        "host_offset_5p_0based": call.offset if call.offset is not None else "",
        "n_carriers": len(samples),
        "carrier_count_class": str(len(samples)),
        "samples": "|".join(samples),
        "presence_bitmap": "".join(
            "1" if s in set(samples) else "0" for s in SAMPLES
        ),
        "private": len(samples) == 1,
        "private_to_sample": samples[0] if len(samples) == 1 else "",
        "source_record_count": sum(m.source_count for m in site.members),
        "source_call_ids": "|".join(i for m in site.members for i in m.source_ids),
        "match_window_bp": window,
    }


def build_sites_table(sites: Sequence[Site], window: int) -> pd.DataFrame:
    return pd.DataFrame(
        [site_row(site, f"US{index:05d}", window) for index, site in enumerate(sites, 1)]
    )
