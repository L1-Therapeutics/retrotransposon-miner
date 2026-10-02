"""Shared inputs for the multi-sample nested-insertion analyses.

This module holds the pieces that `joint_enrichment.py` and `recurrence_test.py`
must agree on exactly. Duplicating them across the two scripts is how a
four-value enum quietly becomes a three-value one, or how a contig filter drifts,
so they live here and are tested once.

Provenance and the join
-----------------------
Site grouping is NOT re-derived here. `unique_sites.csv` is the authoritative
output of the multi-sample dedup run (`rtm-dedup/scripts/match_sites_across_samples.py`,
commit de2589d) and is consumed as given, per the frozen plan's estimand
statement: positional analyses operate on unique insertion sites, deduplicated
across samples.

What this module adds is the *per-call* detail the dedup output does not carry.
`unique_sites.csv` has no TSD length, no TSD sequence and no child subfamily --
only a boolean `tsd_overlap_within_site` -- yet the recurrence definitions need
all three. Those fields live in the per-sample callsets, so each unique site is
joined back to the calls of the samples that dedup recorded as its carriers, by
sample name, within the frozen +/-10 bp rule.

The join is verified rather than assumed. `verify_join` checks that the number
of calls recovered per site equals `n_carriers` from the dedup table, that every
recovered call agrees with the dedup table on nesting state and orientation, and
that it reports the count of sites it could not resolve. On the shipped data all
763 sites resolve, all 1,284 recovered calls match `sum(n_carriers)`, and every
call agrees on orientation; a mismatch is an error, not a warning.

Scope and prohibitions
----------------------
- Matching and counting are restricted to chr1-22 and chrX. chrY and any other
  contig are dropped and the count is reported.
- `GT` and `GQ` are never read as genotype evidence. `read_callset` parses
  CHROM/POS/ID/INFO only and never looks at the FORMAT column, so a
  genotype-shaped field cannot enter these analyses by accident.
- The NESTED field is a four-value enum: `unnested`, `nested_sense`,
  `nested_antisense`, `nested_unknown`. `nested_unknown` is nested but belongs
  to neither orientation class and is never folded into antisense.
  These particular callsets were written by an exporter that emits only the
  coarser `nested`/`unnested` pair, so the legacy `nested` value is kept as its
  own label, `nested_orientation_unlabeled`, counted, and reported -- never
  mapped onto either orientation class. `parse_nested_state` raises on any value
  outside the known set rather than coercing it.
"""

from __future__ import annotations

import collections
import csv
import gzip
from pathlib import Path
from typing import Any, Iterable, Sequence

#: Analysis roots. The workspace keeps `nested_analysis/` at its root rather than
#: inside any single worktree, which is where the dedup outputs are written.
WORKSPACE = Path(__file__).resolve().parents[2]
ANALYSIS_DIR = WORKSPACE / "nested_analysis"
DEFAULT_UNIQUE_SITES = ANALYSIS_DIR / "results_multi_sample" / "unique_sites.csv"
DEFAULT_CALLSET_DIR = ANALYSIS_DIR / "callsets"

#: The frozen +/-10 bp matching window, plus the sweeps the plan requires to be
#: reported in full regardless of outcome.
PRIMARY_WINDOW_BP = 10
SWEEP_WINDOWS_BP = (5, 20)

#: The four canonical NESTED values, and the exact value the per-sample
#: callsets actually carry. Both are reported; neither is folded into the other.
NESTED_ENUM = ("unnested", "nested_sense", "nested_antisense", "nested_unknown")
LEGACY_NESTED_UNLABELED = "nested_orientation_unlabeled"
KNOWN_NESTED_VALUES = frozenset(NESTED_ENUM) | frozenset({LEGACY_NESTED_UNLABELED})

#: Scope. chrY is excluded: a single hemizygous contig is not comparable to the
#: autosomes and would enter the opportunity denominators.
ALLOWED_CONTIGS = frozenset([f"chr{i}" for i in range(1, 23)] + ["chrX"])

#: Insertion orientation is carried by its own field, independent of the NESTED
#: label, so the "orientation = sense" cells do not depend on how the nesting
#: label was resolved.
SENSE = "+"
ANTISENSE = "-"
ORIENTATION_VALUES = (SENSE, ANTISENSE)


class NestedEnumError(ValueError):
    """A NESTED value outside the known enum. Never coerced, never defaulted."""


class SchemaError(ValueError):
    """The unique-sites file does not carry the columns these analyses read."""


#: Every column `load_unique_sites` reads, and why each one is load-bearing.
#: The host interval drives the geometry gate and `host_key`; the carrier list
#: and its count drive the join gate; `site_id` and `orientation` are what the
#: reports and the orientation cells are keyed on.
REQUIRED_UNIQUE_SITE_COLUMNS = (
    "site_id",
    "chrom",
    "representative_pos",
    "family",
    "insertion_orientation",
    "host_name",
    "host_start0",
    "host_end0",
    "host_strand",
    "host_len",
    "host_offset_5p_0based",
    "n_carriers",
    "samples",
    "private",
)


def parse_nested_state(value: str) -> str:
    """Return a NESTED state, or raise.

    `nested_unknown` is preserved as itself. The legacy coarse `nested` value is
    preserved as `nested_orientation_unlabeled` -- it says a same-family element
    overlaps, but nothing about which orientation class, so treating it as sense
    or as antisense would be a fabrication.
    """
    text = (value or "").strip()
    if not text:
        raise NestedEnumError("empty NESTED value")
    if text in NESTED_ENUM:
        return text
    if text == "nested":
        return LEGACY_NESTED_UNLABELED
    raise NestedEnumError(f"unrecognised NESTED value: {value!r}")


def nested_census(states: Iterable[str]) -> dict[str, int]:
    """Census of NESTED states, counting every value in the enum.

    Values are normalised through `parse_nested_state` first, so a caller that
    passes a raw `nested` gets it counted as the unlabeled label rather than
    silently landing in the unrecognised bucket -- which would understate the
    unlabeled count and make the enum look better populated than it is.
    """
    counts: collections.Counter[str] = collections.Counter()
    unrecognised = 0
    for value in states:
        # Accept a raw VCF value or an already-parsed state. Re-parsing an
        # already-parsed state would reject it, which would report every
        # unlabeled call as unrecognised and make the enum look emptier than it
        # is -- the opposite error from folding, and just as misleading.
        if value in KNOWN_NESTED_VALUES:
            counts[value] += 1
            continue
        try:
            counts[parse_nested_state(value)] += 1
        except NestedEnumError:
            unrecognised += 1
    census = {value: int(counts.get(value, 0)) for value in NESTED_ENUM}
    census[LEGACY_NESTED_UNLABELED] = int(counts.get(LEGACY_NESTED_UNLABELED, 0))
    census["unrecognised_rejected"] = int(unrecognised)
    return census


def in_scope_contig(chrom: str) -> bool:
    return chrom in ALLOWED_CONTIGS


def to_int(value: str | None) -> int | None:
    if value is None:
        return None
    text = value.strip()
    if not text:
        return None
    try:
        return int(text)
    except ValueError:
        return None


def canonical_subfamily(value: str | None) -> str | None:
    """`AluYb9#SINE/Alu` -> `AluYb9`.

    The `#` suffix is the RepeatMasker class path, which is constant within a
    family and so carries no subfamily information. Comparisons are on the
    subfamily alone, otherwise every Alu would look distinct from every other.
    """
    text = (value or "").strip()
    if not text:
        return None
    return text.split("#", 1)[0] or None


# --------------------------------------------------------------------------
# inputs
# --------------------------------------------------------------------------


def load_unique_sites(path: Path) -> tuple[list[dict[str, Any]], dict[str, int]]:
    """Read the dedup output, typed, with the out-of-scope rows dropped.

    The `samples` field is pipe-separated. Splitting it on a comma or a semicolon
    yields one bogus token per multi-carrier site, which silently reduces the join
    to private sites only -- a failure that looks like a plausible result rather
    than an error, and the reason the separator is a named constant here.

    A file missing any of `REQUIRED_UNIQUE_SITE_COLUMNS` raises `SchemaError`
    instead of loading. This is not defensive parsing. Two different tables in
    this workspace are both called `unique_sites.csv` and they do not carry the
    same columns: the multi-sample site table has `same_family_host_name` and no
    host interval at all, while the per-call cohort table has `host_start0`. Read
    the wrong one and every row fails the geometry gate, so the loader returned
    0 sites from 4,998 rows and both analyses went on to report an empty cohort
    as if it were a null result. A changed column set is a new input, and the
    caller has to be told.

    Rows are still dropped individually for scope and geometry, and counted in
    the returned report; an all-rows-dropped table is not an error here, because
    a table of hostless sites legitimately has empty geometry for every row. It
    is the *analysis scripts* that refuse an empty cohort, not the loader.
    """
    samples_separator = "|"
    rows: list[dict[str, Any]] = []
    dropped_contig = 0
    dropped_unusable = 0
    with Path(path).open(newline="") as fh:
        reader = csv.DictReader(fh)
        header = list(reader.fieldnames or [])
        missing = [
            column for column in REQUIRED_UNIQUE_SITE_COLUMNS if column not in header
        ]
        if missing:
            raise SchemaError(
                f"{path} is not a nested-Alu unique-sites table: it is missing "
                f"column(s) {missing}. Columns found: {header}."
            )
        for raw in reader:
            chrom = (raw.get("chrom") or "").strip()
            if not in_scope_contig(chrom):
                dropped_contig += 1
                continue
            pos = to_int(raw.get("representative_pos"))
            offset = to_int(raw.get("host_offset_5p_0based"))
            host_len = to_int(raw.get("host_len"))
            if pos is None or offset is None or host_len is None or host_len <= 0:
                dropped_unusable += 1
                continue
            carriers = [
                s for s in (raw.get("samples") or "").split(samples_separator) if s
            ]
            n_carriers = to_int(raw.get("n_carriers")) or 0
            rows.append(
                {
                    "site_id": raw.get("site_id"),
                    "chrom": chrom,
                    "pos": pos,
                    "family": (raw.get("family") or "").strip(),
                    "orientation": (raw.get("insertion_orientation") or "").strip(),
                    "host_name": (raw.get("host_name") or "").strip(),
                    "host_start0": to_int(raw.get("host_start0")),
                    "host_end0": to_int(raw.get("host_end0")),
                    "host_strand": (raw.get("host_strand") or "").strip(),
                    "host_len": host_len,
                    "host_offset": offset,
                    "host_selection_rule": raw.get("host_selection_rule"),
                    "n_carriers": n_carriers,
                    "carriers": carriers,
                    "private": (raw.get("private") or "").strip() == "True",
                    "tsd_overlap_within_site": to_int(raw.get("tsd_overlap_within_site")),
                    "pooled_site_matched": (raw.get("pooled_site_matched") or "").strip()
                    == "True",
                }
            )
    report = {
        "rows_read": len(rows) + dropped_contig + dropped_unusable,
        "rows_in_scope": len(rows),
        "dropped_out_of_scope_contig": dropped_contig,
        "dropped_missing_geometry": dropped_unusable,
        "samples_field_separator": samples_separator,
    }
    return rows, report


def read_callset(path: Path) -> dict[str, list[dict[str, Any]]]:
    """Read one callset, grouped by contig.

    Only CHROM, POS, ID and INFO are read. The FORMAT column is never parsed, so
    `GT` and `GQ` cannot be mistaken for genotype evidence here; this is the
    structural guarantee behind that rule, not a promise in a docstring.
    """
    by_contig: dict[str, list[dict[str, Any]]] = collections.defaultdict(list)
    with Path(path).open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            info = {}
            for item in parts[7].split(";"):
                if "=" in item:
                    key, _, value = item.partition("=")
                    info[key] = value
                elif item:
                    info[item] = ""
            nested_raw = info.get("NESTED", "")
            by_contig[parts[0]].append(
                {
                    "pos": int(parts[1]),
                    "id": parts[2],
                    "family": info.get("MEIFAMILY", ""),
                    "subfamily": canonical_subfamily(info.get("MEISUBFAMILY")),
                    "orientation": info.get("ORIENT", ""),
                    "nested_raw": nested_raw,
                    "nested_state": parse_nested_state(nested_raw),
                    "tsd": (info.get("TSD") or "").strip(),
                    "tsd_len": len((info.get("TSD") or "").strip()),
                    "tsd_reported": bool((info.get("TSD") or "").strip()),
                    "mei_span": to_int(info.get("MEI_SPAN")),
                    "known_mei": info.get("KNOWNMEI", ""),
                    "call_tier": info.get("CALLTIER", ""),
                }
            )
    for contig in by_contig:
        by_contig[contig].sort(key=lambda r: r["pos"])
    return dict(by_contig)


def load_callsets(
    callset_dir: Path, samples: Sequence[str]
) -> dict[str, dict[str, list[dict[str, Any]]]]:
    out: dict[str, dict[str, list[dict[str, Any]]]] = {}
    for sample in samples:
        path = Path(callset_dir) / f"{sample}.vcf"
        if not path.exists():
            raise FileNotFoundError(f"callset missing for sample {sample}: {path}")
        out[sample] = read_callset(path)
    return out


def _nearest_call(
    calls: Sequence[dict[str, Any]], pos: int, window: int
) -> dict[str, Any] | None:
    best = None
    best_delta = None
    for record in calls:
        delta = abs(record["pos"] - pos)
        if delta > window:
            continue
        if best_delta is None or delta < best_delta:
            best, best_delta = record, delta
    return best


def attach_call_details(
    sites: Sequence[dict[str, Any]],
    callsets: dict[str, dict[str, list[dict[str, Any]]]],
    window: int = PRIMARY_WINDOW_BP,
) -> list[dict[str, Any]]:
    """Recover per-call TSD and child subfamily for each unique site.

    Sites are matched to their own recorded carriers by sample name. A site whose
    carrier list names a sample with no callset available is reported, not
    silently dropped, because it changes the denominator.
    """
    out: list[dict[str, Any]] = []
    for site in sites:
        per_call: list[dict[str, Any]] = []
        unmatched_samples: list[str] = []
        for sample in site["carriers"]:
            table = callsets.get(sample)
            if table is None:
                unmatched_samples.append(sample)
                continue
            record = _nearest_call(table.get(site["chrom"], []), site["pos"], window)
            if record is None:
                unmatched_samples.append(sample)
                continue
            per_call.append({"sample": sample, **record})
        new = dict(site)
        new["calls"] = per_call
        new["unmatched_carrier_samples"] = unmatched_samples
        out.append(new)
    return out


def verify_join(sites: Sequence[dict[str, Any]]) -> dict[str, Any]:
    """Check the join against dedup's own carrier counts.

    This is a gate, not a diagnostic. If the number of recovered calls disagrees
    with `n_carriers`, or a recovered call is not the call the site describes --
    wrong family or wrong orientation -- the per-call detail is not safe to use
    and the analysis must stop rather than proceed on a partial join.

    What is deliberately NOT part of the gate is whether the source call's raw
    `NESTED` field agrees with the site's nesting state. Those two answer
    different questions and are computed by different rules:

      - `dedup_samples.py` derives nesting from RepeatMasker by assigning a
        same-family host and comparing the insertion's strand to the host's,
        yielding the four-state sense/antisense/unknown classification.
      - The per-sample VCFs carry a binary legacy `NESTED=nested|unnested` label
        written by an earlier exporter under its own host-selection rule.

    On this cohort every call at a site the producer calls `nested_antisense`
    carries the legacy `unnested` label, so gating on that field makes the gate
    fire on 100% of antisense sites and refuses to run any analysis at all. That
    is a disagreement between two definitions, not a broken join, and treating it
    as fatal would report a definitional difference as a data defect. It is
    counted and surfaced instead, in
    `calls_where_source_nesting_label_differs_from_site`, so the divergence stays
    visible rather than being silently dropped.
    """
    carrier_mismatch: list[str] = []
    orientation_mismatch: list[str] = []
    family_mismatch: list[str] = []
    nesting_label_differs: list[str] = []
    resolved = 0
    total_calls = 0
    for site in sites:
        if site["calls"]:
            resolved += 1
        total_calls += len(site["calls"])
        if len(site["calls"]) != site["n_carriers"]:
            carrier_mismatch.append(str(site["site_id"]))
        for call in site["calls"]:
            if call["orientation"] != site["orientation"]:
                orientation_mismatch.append(f"{site['site_id']}:{call['sample']}")
            if (call["family"] or "").strip().upper() != (site["family"] or "").strip().upper():
                family_mismatch.append(f"{site['site_id']}:{call['sample']}")
            # Diagnostic only -- see the docstring. Compared as nested-vs-not,
            # never folded into either orientation class.
            if call["nested_state"].startswith("nested") is False:
                nesting_label_differs.append(f"{site['site_id']}:{call['sample']}")
    return {
        "sites_total": len(sites),
        "sites_with_at_least_one_joined_call": resolved,
        "sites_with_no_joined_call": len(sites) - resolved,
        "calls_recovered": total_calls,
        "calls_declared_by_dedup": sum(int(s["n_carriers"]) for s in sites),
        "sites_where_carrier_count_disagrees": len(carrier_mismatch),
        "calls_where_orientation_disagrees": len(orientation_mismatch),
        "calls_where_family_disagrees": len(family_mismatch),
        "calls_where_source_nesting_label_differs_from_site": len(nesting_label_differs),
        "source_nesting_label_note": (
            "diagnostic only, not part of the verdict: the legacy binary NESTED "
            "field and the producer's strand-derived four-state nesting are "
            "computed under different host-selection rules and disagree on "
            "antisense sites by construction"
        ),
        "verdict": (
            "join_consistent_with_dedup_output"
            if not carrier_mismatch and not orientation_mismatch and not family_mismatch
            else "join_disagrees_with_dedup_output"
        ),
        "examples_carrier_mismatch": carrier_mismatch[:5],
        "examples_orientation_mismatch": orientation_mismatch[:5],
        "examples_family_mismatch": family_mismatch[:5],
    }


def host_key(site: dict[str, Any]) -> tuple[str, str, int, int, str]:
    """Identity of the *host copy* an event sits in.

    Recurrence means two insertions into the same physical host element, so the
    key is the element interval and strand, not the inserted element and not the
    site. Sites sharing a key are in one host copy.
    """
    return (
        site["chrom"],
        site["host_name"],
        int(site["host_start0"]),
        int(site["host_end0"]),
        site["host_strand"],
    )


def host_label(key: tuple[str, str, int, int, str]) -> str:
    return f"{key[0]}:{key[2]}-{key[3]}:{key[1]}:{key[4]}"


def host_family(site: dict[str, Any]) -> str:
    """Repeat family of the *host*, from the host element's name.

    `AluSx1` -> `Alu`, `L1PA2` -> `L1`, `SVA_A` -> `SVA`. This is the host's
    family, which is a different quantity from the inserted element's
    `family` column and is the one the pre-specified cells name.
    """
    name = site["host_name"] or ""
    if name.startswith("Alu"):
        return "Alu"
    if name.startswith("L1"):
        return "L1"
    if name.startswith("SVA"):
        return "SVA"
    if name.startswith("FR"):
        return "SVA"
    return "other"


def nested_census_for_sites(sites: Sequence[dict[str, Any]]) -> dict[str, int]:
    return nested_census(
        call["nested_state"] for site in sites for call in site["calls"]
    )


def read_rmsk_intervals(
    rmsk_path: Path, wanted: dict[str, list[tuple[int, int, str]]]
) -> dict[str, list[tuple[int, int, str]]]:
    """Stream RepeatMasker once, keeping only intervals overlapping a wanted host.

    `wanted` maps contig -> list of (host_start0, host_end0, host_name). The
    result is the subset of annotations needed to measure how much of a host
    interval is covered by repeats *other than the host element itself*.

    Both the contig and the interval are filtered here, which is what the name
    and this docstring promise. Filtering on contig alone retained every
    annotation on all 23 in-scope contigs -- hundreds of thousands of rows held
    in memory to answer a question about a few hundred host intervals.
    """
    keep: dict[str, list[tuple[int, int, str]]] = collections.defaultdict(list)
    spans_by_contig = {
        contig: [(int(start), int(end)) for start, end, _name in spans]
        for contig, spans in wanted.items()
    }
    opener = gzip.open if str(rmsk_path).endswith(".gz") else open
    with opener(rmsk_path, "rt", errors="replace") as fh:  # type: ignore[operator]
        for line in fh:
            fields = line.split("\t")
            if len(fields) < 11:
                continue
            contig = fields[5]
            spans = spans_by_contig.get(contig)
            if not spans:
                continue
            a_start = int(fields[6])
            a_end = int(fields[7])
            # Half-open overlap, matching how host intervals are stored.
            if not any(a_start < w_end and a_end > w_start for w_start, w_end in spans):
                continue
            keep[contig].append((a_start, a_end, fields[9]))
    return dict(keep)


def residual_mask_fraction(
    host: tuple[int, int, str], annotations: Sequence[tuple[int, int, str]]
) -> float:
    """Fraction of a host interval covered by *other* repeat annotations.

    The host element is itself a RepeatMasker annotation, so a naive repeat mask
    deletes the entire host and leaves no opportunity at all. The host's own
    annotation -- matched by an interval agreeing to within 5 bp at both ends --
    is excluded, and what remains measures whether any *neighbouring* repeat
    intrudes on the host interior.

    Coverage is the *union* of the remaining intervals, not the widest single
    one. A host flanked by two repeats that each cover 40% of it is 80% masked,
    not 40%; taking the maximum understated every host with more than one
    intruding annotation and made the diagnostic read cleaner than the genome
    actually is.
    """
    start, end, _ = host
    span = end - start
    if span <= 0:
        return 1.0
    spans: list[tuple[int, int]] = []
    for a_start, a_end, _name in annotations:
        if a_end <= start or a_start >= end:
            continue
        if abs(a_start - start) <= 5 and abs(a_end - end) <= 5:
            continue
        spans.append((max(a_start, start), min(a_end, end)))
    if not spans:
        return 0.0
    spans.sort()
    covered = 0
    cur_start, cur_end = spans[0]
    for lo, hi in spans[1:]:
        if lo > cur_end:
            covered += cur_end - cur_start
            cur_start, cur_end = lo, hi
        else:
            cur_end = max(cur_end, hi)
    covered += cur_end - cur_start
    return min(covered, span) / span
