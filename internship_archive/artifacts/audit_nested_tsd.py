#!/usr/bin/env python3
"""Audit saved HG03086 TSD/poly(A) call fields and nested-host coordinates.

This deliberately does NOT call a consensus TSD a validated TSD: sequence access to
both reference flanks and original reads is needed to reconstruct breakpoint boundaries.
The script only emits QC triage fields and repeat-coordinate reconstructions.
"""
from __future__ import annotations

import bisect
import csv
import gzip
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
VCF = ROOT / "HG03086_classifier_ge_0.997.genes.snpeff (1).vcf"
PER_CALL = ROOT / "artifacts/selfins_per_call.csv"
RMSK = ROOT / "nested_analysis/data/rmsk.txt.gz"
OUT = ROOT / "artifacts"


def family_for_rmsk(rep_family: str, rep_class: str) -> str:
    token = f"{rep_family} {rep_class}".upper()
    if "ALU" in token:
        return "ALU"
    if "SVA" in token:
        return "SVA"
    if "LINE1" in token or re.search(r"(^|[^A-Z0-9])L1([^A-Z0-9]|$)", token):
        return "LINE1"
    return ""


def read_vcf() -> dict[tuple[str, int], dict[str, str]]:
    records = {}
    with VCF.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            info = {}
            for item in fields[7].split(";"):
                if "=" in item:
                    key, value = item.split("=", 1)
                    info[key] = value
            key = (fields[0], int(fields[1]))
            if key in records:
                raise ValueError(f"duplicate VCF coordinate {key}")
            records[key] = info
    return records


def read_saved_calls() -> list[dict[str, str]]:
    with PER_CALL.open(newline="") as handle:
        return list(csv.DictReader(handle))


def support_counts(raw: str) -> dict[str, str]:
    result = {}
    for entry in raw.split("|"):
        if ":" in entry:
            key, value = entry.split(":", 1)
            result[key] = value
    return result


def tsd_triage(seq: str, family: str) -> tuple[int, float | None, float | None, str]:
    if not seq or seq == ".":
        return 0, None, None, "missing"
    seq = seq.upper()
    n = len(seq)
    at_fraction = (seq.count("A") + seq.count("T")) / n if n else None
    dominant = max(Counter(seq).values()) / n if n else None
    flags = []
    if family == "ALU":
        flags.append("within_5_27bp_Wagstaff_de_novo_range_only" if 5 <= n <= 27
                     else "outside_5_27bp_Wagstaff_de_novo_range_review")
    elif family == "LINE1":
        if 9 <= n <= 11:
            flags.append("9_11bp_L1_short_TSD_false_positive_caution")
        else:
            flags.append("length_not_by_itself_evidence_for_L1_TSD")
        if n < 9 or n > 60:
            flags.append("outside_9_60bp_TSDfinder_search_or_review_limit")
    elif family == "SVA":
        flags.append("within_4_18bp_Damert_chr19_range_only" if 4 <= n <= 18
                     else "outside_4_18bp_Damert_chr19_range_review")
        if n == 42:
            flags.append("Damert_chr19_had_one_42bp_exception")
    else:
        flags.append("no_family_specific_range_assigned")
    if at_fraction is not None and at_fraction >= 0.80:
        flags.append("AT_rich_low_complexity_review")
    if dominant is not None and dominant >= 0.80:
        flags.append("single_base_low_complexity_review")
    return n, at_fraction, dominant, ";".join(flags)


def reconstruct_hosts(calls: list[dict[str, str]], vcf: dict) -> list[dict]:
    nested = [r for r in calls if r["nested"].lower() == "true"]
    wanted = defaultdict(list)
    for row in nested:
        key = (row["chrom"], row["MEIFAMILY"])
        wanted[key].append((int(row["pos"]) - 1, row))
    for key in wanted:
        wanted[key].sort(key=lambda item: item[0])
    starts = {key: [item[0] for item in values] for key, values in wanted.items()}
    overlaps: dict[tuple[str, int], list[tuple[int, int, str, str, str]]] = defaultdict(list)

    # Stream RepeatMasker once and retain only intervals that actually contain one of
    # the 320 saved same-family nested call coordinates.
    with gzip.open(RMSK, "rt") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 13:
                continue
            chrom, family = fields[5], family_for_rmsk(fields[12], fields[11])
            key = (chrom, family)
            if not family or key not in wanted:
                continue
            start0, end0 = int(fields[6]), int(fields[7])
            lo = bisect.bisect_left(starts[key], start0)
            hi = bisect.bisect_left(starts[key], end0)
            for pos0, row in wanted[key][lo:hi]:
                if pos0 < end0:
                    name = fields[10]
                    strand = fields[9]
                    overlaps[(row["chrom"], int(row["pos"]))].append(
                        (start0, end0, strand, name, family)
                    )

    result = []
    for row in calls:
        chrom, pos = row["chrom"], int(row["pos"])
        info = vcf[(chrom, pos)]
        tsd = info.get("TSD", "")
        tsd_len, at_frac, dominant, triage = tsd_triage(tsd, row["MEIFAMILY"])
        support = support_counts(info.get("SUPPORT", ""))
        rec = {
            "chrom": chrom,
            "pos": pos,
            "family": row["MEIFAMILY"],
            "subfamily": row["MEISUBFAMILY"],
            "insertion_orientation": row["ORIENT"],
            "nested_same_family": row["nested"],
            "saved_host_name": row["ref_name"],
            "saved_host_strand": row["ref_strand"],
            "saved_host_length": row["ref_len"],
            "tsd_sequence_from_vcf": tsd,
            "tsd_length_bp": tsd_len,
            "tsd_AT_fraction": "" if at_frac is None else round(at_frac, 4),
            "tsd_dominant_base_fraction": "" if dominant is None else round(dominant, 4),
            "tsd_QC_triage_not_validation": triage,
            "polyA_min_bp_classifier_field": info.get("POLYA_MIN_BP", ""),
            "nested_value_vcf": info.get("NESTED", ""),
            "call_tier": info.get("CALLTIER", ""),
            "classifier_gold_score": info.get("L1TXGOLDSCORE", ""),
            "insertion_score": info.get("INSERTIONSCORE", ""),
            "known_mei": info.get("KNOWNMEI", ""),
            "known_source": info.get("KNOWN_SRC", ""),
            "sample_status": info.get("SAMPLE_STATUS", ""),
            "SR_L": support.get("SR_L", ""),
            "SR_R": support.get("SR_R", ""),
            "DPE_L": support.get("DPE_L", ""),
            "DPE_R": support.get("DPE_R", ""),
            "MEI_MAPPED": support.get("MEI_MAPPED", ""),
            "polyA_MAPPED": support.get("polyA_MAPPED", ""),
            "host_start0_reconstructed": "",
            "host_end0_reconstructed": "",
            "host_strand_reconstructed": "",
            "host_name_reconstructed": "",
            "host_position_5prime_0based": "",
            "host_reconstruction_matches_saved": "",
        }
        if row["nested"].lower() == "true":
            candidates = overlaps.get((chrom, pos), [])
            saved_name = row["ref_name"].split("#")[0]
            named = [item for item in candidates if item[3] == saved_name]
            pool = named or candidates
            # The saved analysis resolves overlapping hosts by greatest span, then
            # leftmost start; no orientation-dependent tie break.
            if pool:
                start0, end0, strand, name, fam = sorted(
                    pool, key=lambda item: (-(item[1] - item[0]), item[0], item[3])
                )[0]
                pos0 = pos - 1
                offset = pos0 - start0 if strand == "+" else (end0 - 1 - pos0)
                rec.update({
                    "host_start0_reconstructed": start0,
                    "host_end0_reconstructed": end0,
                    "host_strand_reconstructed": strand,
                    "host_name_reconstructed": name,
                    "host_position_5prime_0based": offset,
                    "host_reconstruction_matches_saved": str(
                        name == saved_name and strand == row["ref_strand"]
                        and (end0 - start0) == int(float(row["ref_len"]))
                    ),
                })
            else:
                rec["host_reconstruction_matches_saved"] = "no_matching_rmsk_interval"
        result.append(rec)
    return result


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    calls = read_saved_calls()
    vcf = read_vcf()
    coords = {(r["chrom"], int(r["pos"])) for r in calls}
    if len(calls) != 1852 or len(vcf) != 1852 or coords != set(vcf):
        raise ValueError(f"expected exact 1,852-call join; csv={len(calls)} vcf={len(vcf)}")
    audit = reconstruct_hosts(calls, vcf)
    write_csv(OUT / "selfins_tsd_call_qc.csv", audit)

    summary = {
        "inputs": {
            "vcf": VCF.name,
            "per_call": "artifacts/selfins_per_call.csv",
            "repeatmasker": "nested_analysis/data/rmsk.txt.gz",
            "join_key": "exact (chrom, POS); all 1,852 coordinates must match",
        },
        "warnings": [
            "VCF TSD string and POLYA_MIN_BP are classifier/caller fields, not independently reconstructed breakpoint evidence.",
            "No FASTA or BAM/CRAM is available at this workspace root; do not calculate a TPRT triad or call a TSD authentic from sequence length alone.",
            "A length inside a literature range is a QC descriptor only; it is not proof of a duplication at both junctions.",
        ],
        "overall": {},
        "by_family": {},
        "nested_alu_near_full_host_position_profile": [],
        "sva_nested_calls": [],
    }
    present = [r for r in audit if r["tsd_length_bp"] > 0]
    summary["overall"] = {
        "n_calls": len(audit),
        "n_tsd_strings": len(present),
        "n_tsd_missing": len(audit) - len(present),
        "TSD_length_counts_all_calls": dict(sorted(Counter(r["tsd_length_bp"] for r in present).items())),
        "n_polyA_min_missing": sum(not r["polyA_min_bp_classifier_field"] for r in audit),
        "n_nested": sum(r["nested_same_family"] == "True" for r in audit),
    }
    for family in ("ALU", "LINE1", "SVA"):
        subset = [r for r in audit if r["family"] == family]
        nested = [r for r in subset if r["nested_same_family"] == "True"]
        family_present = [r for r in subset if r["tsd_length_bp"] > 0]
        summary["by_family"][family] = {
            "n_calls": len(subset),
            "n_nested": len(nested),
            "n_tsd_strings": len(family_present),
            "n_tsd_missing": len(subset) - len(family_present),
            "n_nested_with_tsd_string": sum(r["tsd_length_bp"] > 0 for r in nested),
            "nested_tsd_length_counts": dict(sorted(Counter(r["tsd_length_bp"] for r in nested if r["tsd_length_bp"] > 0).items())),
            "n_nested_flag_AT_rich_or_single_base_low_complexity": sum(
                "low_complexity_review" in r["tsd_QC_triage_not_validation"] for r in nested
            ),
        }
    alu = [r for r in audit if r["family"] == "ALU" and r["nested_same_family"] == "True"
           and r["host_position_5prime_0based"] != ""
           and 280 <= int(float(r["saved_host_length"])) <= 320]
    for start in range(0, 320, 20):
        bin_rows = [r for r in alu if start <= int(r["host_position_5prime_0based"]) < start + 20]
        lens = [r["tsd_length_bp"] for r in bin_rows if r["tsd_length_bp"] > 0]
        polya = [int(r["polyA_min_bp_classifier_field"]) for r in bin_rows if r["polyA_min_bp_classifier_field"].isdigit()]
        summary["nested_alu_near_full_host_position_profile"].append({
            "bin_start_bp": start,
            "bin_end_bp": start + 20,
            "n_calls": len(bin_rows),
            "n_with_tsd_string": len(lens),
            "median_tsd_string_length": None if not lens else sorted(lens)[len(lens) // 2],
            "n_with_polyA_min": len(polya),
            "median_polyA_min_bp_classifier_field": None if not polya else sorted(polya)[len(polya) // 2],
        })
    summary["nested_alu_near_full_host_position_profile"].append({
        "note": "Reconstructed RMSK/coordinate profile is separate from the saved 29/36 profile; coordinate and host-selection conventions can shift counts. The host position is 0-based from its 5-prime end; it is not a base-resolution Levy alignment."
    })
    summary["position_profile_reconstruction"] = {
        "n_near_full_host_nested_alu": len(alu),
        "n_at_120_140bp": sum(120 <= int(r["host_position_5prime_0based"]) < 140 for r in alu),
        "n_at_280_300bp": sum(280 <= int(r["host_position_5prime_0based"]) < 300 for r in alu),
        "n_reconstructed_host_positions": sum(bool(r["host_position_5prime_0based"]) for r in alu),
    }
    summary["host_reconstruction"] = {
        "n_nested_calls": sum(r["nested_same_family"] == "True" for r in audit),
        "n_reconstructed_matches_saved": sum(r["host_reconstruction_matches_saved"] == "True" for r in audit),
        "n_no_matching_rmsk_interval": sum(r["host_reconstruction_matches_saved"] == "no_matching_rmsk_interval" for r in audit),
        "n_rmsk_host_name_or_length_different": sum(r["host_reconstruction_matches_saved"] == "False" for r in audit),
        "reconstruction_note": "Six saved nested calls had no matching GRCh38 hg38 RMSK same-family interval at the normalized coordinate; two had a same-family interval but a different selected RepeatMasker record. Likely repeat track version or call-boundary differences; these are not silently recoded.",
    }
    summary["sva_nested_calls"] = [
        {k: r[k] for k in (
            "chrom", "pos", "subfamily", "insertion_orientation", "tsd_sequence_from_vcf",
            "tsd_length_bp", "tsd_AT_fraction", "tsd_QC_triage_not_validation",
            "polyA_min_bp_classifier_field", "classifier_gold_score", "known_source",
            "sample_status", "SR_L", "SR_R", "DPE_L", "DPE_R", "MEI_MAPPED", "polyA_MAPPED",
            "host_name_reconstructed", "host_strand_reconstructed", "host_position_5prime_0based",
        )}
        for r in audit if r["family"] == "SVA" and r["nested_same_family"] == "True"
    ]
    with (OUT / "selfins_tsd_audit_summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")
    print(json.dumps({
        "calls": summary["overall"]["n_calls"],
        "tsd_strings": summary["overall"]["n_tsd_strings"],
        "missing_tsd": summary["overall"]["n_tsd_missing"],
        "polyA_min_missing": summary["overall"]["n_polyA_min_missing"],
        "nested": summary["overall"]["n_nested"],
        "host_reconstruction": summary["host_reconstruction"],
        "sva_nested": len(summary["sva_nested_calls"]),
        "outputs": ["artifacts/selfins_tsd_call_qc.csv", "artifacts/selfins_tsd_audit_summary.json"],
    }, indent=2))


if __name__ == "__main__":
    main()
