"""Sample-aware overlap of RTM MEI calls against MELT and ONT genotype callsets.

Explore-first: extract per-sample carriers, measure recall, sweep a simple
score/tier cutoff, and list RTM calls that are novel to both truth sets.
"""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

MEI_FAMILIES = frozenset({"ALU", "LINE1", "SVA"})
_SCORE_COLUMNS = (
    "insertion_model_score",
    "read_support_heuristic_score",
    "coherence_score",
)
_TIER_RANK = {
    "high_conf_two_sided": 4,
    "provisional_one_sided": 3,
    "mei_with_complex": 2,
    "complex_ins": 1,
    "none": 0,
}


@dataclass(frozen=True)
class Variant:
    chrom: str
    pos: int
    end: int
    variant_id: str
    family: str
    source: str
    genotype: str = ""
    svtype: str = ""
    svlen: int = -1
    strand: str = ""
    extra: dict[str, Any] = field(default_factory=dict)

    @property
    def key(self) -> str:
        return self.variant_id or f"{self.chrom}:{self.pos}"


def normalize_chrom(chrom: str) -> str:
    text = str(chrom or "").strip()
    if not text:
        return ""
    if text.lower().startswith("chr"):
        return "chr" + text[3:]
    return f"chr{text}"


def normalize_mei_family(token: str) -> str:
    text = (token or "").upper()
    if "ALU" in text:
        return "ALU"
    if "SVA" in text:
        return "SVA"
    if "LINE1" in text or "L1" in text:
        return "LINE1"
    return ""


def normalize_strand(token: str) -> str:
    text = (token or "").strip().upper()
    if text in {"+", "PLUS", "POS", "POSITIVE", "FWD", "FORWARD"}:
        return "+"
    if text in {"-", "MINUS", "NEG", "NEGATIVE", "REV", "REVERSE"}:
        return "-"
    return ""


def strand_from_meinfo(meinfo: str) -> str:
    """MELT MEINFO is NAME,START,END,POLARITY."""
    parts = [p.strip() for p in str(meinfo or "").split(",")]
    if len(parts) >= 4:
        return normalize_strand(parts[3])
    return ""


def is_mei_like_text(*parts: str) -> bool:
    text = " ".join(str(p or "") for p in parts).upper()
    return any(marker in text for marker in ("ALU", "SVA", "LINE", "L1", "MEI", "INS:ME"))


def is_carrier_alleles(allele_indices: Iterable[int | None] | None) -> bool:
    if allele_indices is None:
        return False
    return any(a is not None and int(a) > 0 for a in allele_indices)


def is_carrier_gt_string(gt: str) -> bool:
    token = (gt or "").split(":", 1)[0].replace("|", "/")
    if token in {"", ".", "./.", ".|."}:
        return False
    return any(part not in {"0", "."} for part in token.split("/") if part != "")


def padded_interval(pos: int, end: int, pad_bp: int) -> tuple[int, int]:
    lo = min(int(pos), int(end))
    hi = max(int(pos), int(end))
    pad = max(0, int(pad_bp))
    start = max(0, lo - pad)
    stop = max(start + 1, hi + pad)
    return start, stop


def intervals_overlap(left: Variant, right: Variant, pad_bp: int) -> bool:
    if normalize_chrom(left.chrom) != normalize_chrom(right.chrom):
        return False
    a0, a1 = padded_interval(left.pos, left.end, pad_bp)
    b0, b1 = padded_interval(right.pos, right.end, 0)
    return a0 < b1 and b0 < a1


def match_variants(
    queries: list[Variant],
    references: list[Variant],
    *,
    pad_bp: int,
    require_family: bool = False,
) -> dict[str, list[Variant]]:
    """Map each query key to overlapping reference variants."""
    hits: dict[str, list[Variant]] = {q.key: [] for q in queries}
    if not queries or not references:
        return hits
    by_chrom: dict[str, list[Variant]] = {}
    for ref in references:
        by_chrom.setdefault(normalize_chrom(ref.chrom), []).append(ref)
    for query in queries:
        cand = by_chrom.get(normalize_chrom(query.chrom), [])
        for ref in cand:
            if require_family and query.family and ref.family and query.family != ref.family:
                continue
            if intervals_overlap(query, ref, pad_bp):
                hits[query.key].append(ref)
    return hits


def overlap_metrics(
    truth: list[Variant],
    calls: list[Variant],
    *,
    pad_bp: int,
    require_family: bool = False,
) -> dict[str, float | int]:
    truth_hits = match_variants(truth, calls, pad_bp=pad_bp, require_family=require_family)
    call_hits = match_variants(calls, truth, pad_bp=pad_bp, require_family=require_family)
    n_truth = len(truth)
    n_calls = len(calls)
    tp_truth = sum(1 for v in truth if truth_hits[v.key])
    tp_calls = sum(1 for v in calls if call_hits[v.key])
    return {
        "n_truth": n_truth,
        "n_calls": n_calls,
        "truth_recovered": tp_truth,
        "recall": (tp_truth / n_truth) if n_truth else 0.0,
        "calls_overlapping_truth": tp_calls,
        "novel_calls": n_calls - tp_calls,
        "precision_vs_truth": (tp_calls / n_calls) if n_calls else 0.0,
    }


def rtm_breakpoint(row: dict[str, Any]) -> int:
    for col in (
        "consensus_insertion_breakpoint_pos",
        "insertion_breakpoint_pos",
        "window_start",
    ):
        raw = row.get(col)
        try:
            pos = int(float(raw))
        except (TypeError, ValueError):
            continue
        if pos > 0:
            return pos
    start = _as_int(row.get("window_start"), 0)
    end = _as_int(row.get("window_end"), start)
    if start > 0:
        return (start + max(start, end)) // 2
    return 0


def rtm_family(row: dict[str, Any]) -> str:
    for col in ("consensus_mei_family", "event_family", "mei_family"):
        fam = normalize_mei_family(str(row.get(col, "") or ""))
        if fam:
            return fam
    return ""


def load_rtm_calls(path: Path, chrom: str | None = None) -> tuple[list[Variant], list[dict[str, Any]]]:
    import pandas as pd

    df = pd.read_csv(path, sep="\t", low_memory=False)
    want = normalize_chrom(chrom) if chrom else ""
    rows: list[dict[str, Any]] = []
    variants: list[Variant] = []
    for i, rec in enumerate(df.to_dict(orient="records")):
        row_chrom = normalize_chrom(str(rec.get("chrom", "") or ""))
        if want and row_chrom != want:
            continue
        pos = rtm_breakpoint(rec)
        if pos <= 0:
            continue
        end = _as_int(rec.get("window_end"), pos)
        vid = str(rec.get("locus_id") or rec.get("window_id") or f"rtm_{row_chrom}_{pos}_{i}")
        var = Variant(
            chrom=row_chrom,
            pos=pos,
            end=max(pos, end),
            variant_id=vid,
            family=rtm_family(rec),
            source="rtm",
            extra={
                "insertion_call_tier": str(rec.get("insertion_call_tier") or "none"),
                "insertion_model_score": _as_float(rec.get("insertion_model_score")),
                "read_support_heuristic_score": _as_float(rec.get("read_support_heuristic_score")),
                "coherence_score": _as_float(rec.get("coherence_score")),
                "known_mei_polymorphism": rec.get("known_mei_polymorphism"),
            },
        )
        variants.append(var)
        rec["_rtm_key"] = var.key
        rec["_rtm_pos"] = var.pos
        rec["_rtm_family"] = var.family
        rows.append(rec)
    return variants, rows


def load_melt_sample_meis(
    vcf_path: Path,
    sample: str,
    chrom: str,
) -> list[Variant]:
    want = normalize_chrom(chrom)
    out: list[Variant] = []
    fmt = "%CHROM\t%POS\t%END\t%ID\t%ALT\t%INFO/SVTYPE\t%INFO/SVLEN\t%INFO/MEINFO\t[%GT]\n"
    for row in _query_vcf_rows(vcf_path, chrom=want, sample=sample, fmt=fmt):
        chrom_s, pos_s, end_s, vid, alt, svtype, svlen_s, meinfo, gt = _pad_row(row, 9)
        if not is_carrier_gt_string(gt):
            continue
        if not is_mei_like_text(vid, alt, svtype, meinfo):
            continue
        pos = _as_int(pos_s, 0)
        if pos <= 0:
            continue
        end = _as_int(end_s, pos)
        out.append(
            Variant(
                chrom=normalize_chrom(chrom_s),
                pos=pos,
                end=max(pos, end),
                variant_id=vid if vid not in {"", "."} else f"melt_{chrom_s}_{pos}",
                family=normalize_mei_family(f"{vid} {alt} {meinfo}"),
                source="melt",
                genotype=gt.split(":")[0],
                svtype=svtype or alt,
                svlen=_as_int(svlen_s, -1),
                strand=strand_from_meinfo(meinfo),
            )
        )
    return out


def load_ont_sample_meis(
    svim_path: Path,
    svan_path: Path | None,
    sample: str,
    chrom: str,
) -> list[Variant]:
    want = normalize_chrom(chrom)
    anns = load_svan_index(svan_path, want) if svan_path is not None else {}
    out: list[Variant] = []
    fmt = "%CHROM\t%POS\t%END\t%ID\t%INFO/SVTYPE\t%INFO/SVLEN\t[%GT]\n"
    for row in _query_vcf_rows(svim_path, chrom=want, sample=sample, fmt=fmt):
        chrom_s, pos_s, end_s, vid, svtype, svlen_s, gt = _pad_row(row, 7)
        if not is_carrier_gt_string(gt):
            continue
        pos = _as_int(pos_s, 0)
        if pos <= 0:
            continue
        ann = anns.get(vid) or anns.get(f"{want}:{pos}")
        family = normalize_mei_family(str((ann or {}).get("fam_n", "")))
        if family not in MEI_FAMILIES:
            continue
        itype = str((ann or {}).get("itype_n", "") or "")
        ins_len = _as_int((ann or {}).get("ins_len"), _as_int(svlen_s, -1))
        end = _as_int(end_s, pos)
        out.append(
            Variant(
                chrom=normalize_chrom(chrom_s),
                pos=pos,
                end=max(pos, end),
                variant_id=vid if vid not in {"", "."} else f"ont_{chrom_s}_{pos}",
                family=family,
                source="ont",
                genotype=gt.split(":")[0],
                svtype=svtype or itype or "INS",
                svlen=ins_len,
                strand=normalize_strand(str((ann or {}).get("strand", ""))),
                extra={"itype": itype, "svan_id": str((ann or {}).get("id", ""))},
            )
        )
    return out


def load_svan_index(svan_path: Path, chrom: str) -> dict[str, dict[str, Any]]:
    want = normalize_chrom(chrom)
    index: dict[str, dict[str, Any]] = {}
    fmt = "%CHROM\t%POS\t%ID\t%INFO/FAM_N\t%INFO/ITYPE_N\t%INFO/INS_LEN\t%INFO/STRAND\n"
    for row in _query_vcf_rows(svan_path, chrom=want, sample=None, fmt=fmt):
        chrom_s, pos_s, vid, fam_n, itype_n, ins_len, strand = _pad_row(row, 7)
        pos = _as_int(pos_s, 0)
        payload = {
            "id": vid,
            "fam_n": fam_n,
            "itype_n": itype_n,
            "ins_len": _as_int(ins_len, -1),
            "strand": strand,
            "pos": pos,
        }
        if vid not in {"", "."}:
            index[vid] = payload
        if pos > 0:
            index[f"{normalize_chrom(chrom_s)}:{pos}"] = payload
    return index


def _query_vcf_rows(
    vcf_path: Path,
    *,
    chrom: str,
    sample: str | None,
    fmt: str,
) -> list[list[str]]:
    bcftools = shutil.which("bcftools")
    if bcftools is None:
        raise RuntimeError("bcftools is required to read genotype callsets")
    aliases = [chrom]
    if chrom.startswith("chr"):
        aliases.append(chrom[3:])
    else:
        aliases.append(f"chr{chrom}")
    last_err = ""
    for region in aliases:
        cmd = [bcftools, "query", "-r", region, "-f", fmt, str(vcf_path)]
        if sample:
            cmd = [bcftools, "query", "-r", region, "-s", sample, "-f", fmt, str(vcf_path)]
        proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if proc.returncode == 0:
            return [line.rstrip("\n").split("\t") for line in proc.stdout.splitlines() if line.strip()]
        last_err = (proc.stderr or proc.stdout or "").strip()
    raise RuntimeError(f"bcftools query failed for {vcf_path} {chrom}: {last_err[:800]}")


def _pad_row(row: list[str], n: int) -> list[str]:
    if len(row) >= n:
        return row[:n]
    return row + [""] * (n - len(row))


def label_rtm_calls(
    calls: list[Variant],
    melt: list[Variant],
    ont: list[Variant],
    *,
    pad_bp: int,
    require_family: bool = False,
) -> list[dict[str, Any]]:
    melt_hits = match_variants(calls, melt, pad_bp=pad_bp, require_family=require_family)
    ont_hits = match_variants(calls, ont, pad_bp=pad_bp, require_family=require_family)
    labeled = []
    for call in calls:
        melt_ids = [v.variant_id for v in melt_hits[call.key]]
        ont_ids = [v.variant_id for v in ont_hits[call.key]]
        sources = []
        if melt_ids:
            sources.append("melt")
        if ont_ids:
            sources.append("ont")
        labeled.append(
            {
                "chrom": call.chrom,
                "pos": call.pos,
                "end": call.end,
                "rtm_id": call.variant_id,
                "family": call.family,
                "insertion_call_tier": call.extra.get("insertion_call_tier", "none"),
                "insertion_model_score": call.extra.get("insertion_model_score"),
                "read_support_heuristic_score": call.extra.get("read_support_heuristic_score"),
                "coherence_score": call.extra.get("coherence_score"),
                "overlap_melt": bool(melt_ids),
                "overlap_ont": bool(ont_ids),
                "overlap_any_truth": bool(sources),
                "truth_source": ",".join(sources),
                "melt_ids": ",".join(melt_ids),
                "ont_ids": ",".join(ont_ids),
                "novel": not sources,
            }
        )
    return labeled


def sweep_cutoffs(
    labeled_calls: list[dict[str, Any]],
    melt: list[Variant],
    ont: list[Variant],
    *,
    pad_bp: int,
    require_family: bool = False,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    tiers = ["all", "not_none", "provisional_or_higher", "high_conf_two_sided"]
    for tier in tiers:
        subset = [c for c in labeled_calls if _tier_keeps(str(c.get("insertion_call_tier", "none")), tier)]
        rows.append(_sweep_row("insertion_call_tier", tier, subset, melt, ont, pad_bp, require_family))
    for col in _SCORE_COLUMNS:
        values = sorted({_as_float(c.get(col)) for c in labeled_calls if c.get(col) is not None})
        if not values:
            continue
        # A few representative cutoffs plus 0.
        probes = [0.0, 0.25, 0.4, 0.5, 0.6, 0.75]
        for thr in probes:
            subset = [c for c in labeled_calls if _as_float(c.get(col), -1.0) >= thr]
            rows.append(_sweep_row(col, f">={thr}", subset, melt, ont, pad_bp, require_family))
    return rows


def suggest_heuristic_cutoff(
    sweep_rows: list[dict[str, Any]],
    *,
    min_ont_recall: float = 0.80,
) -> dict[str, Any] | None:
    """Prefer a cutoff that keeps ONT recall high while cutting novel calls."""
    ranked = []
    for row in sweep_rows:
        if row["n_calls"] <= 0:
            continue
        ont_ok = float(row["ont_recall"]) >= min_ont_recall
        ranked.append((not ont_ok, int(row["novel_calls"]), -float(row["ont_recall"]), row))
    if not ranked:
        return None
    ranked.sort(key=lambda x: (x[0], x[1], x[2]))
    best = dict(ranked[0][3])
    best["meets_min_ont_recall"] = float(best["ont_recall"]) >= min_ont_recall
    return best


def _sweep_row(
    cutoff_name: str,
    cutoff_value: str,
    subset_rows: list[dict[str, Any]],
    melt: list[Variant],
    ont: list[Variant],
    pad_bp: int,
    require_family: bool,
) -> dict[str, Any]:
    calls = [
        Variant(
            chrom=str(r["chrom"]),
            pos=int(r["pos"]),
            end=int(r["end"]),
            variant_id=str(r["rtm_id"]),
            family=str(r.get("family") or ""),
            source="rtm",
        )
        for r in subset_rows
    ]
    melt_m = overlap_metrics(melt, calls, pad_bp=pad_bp, require_family=require_family)
    ont_m = overlap_metrics(ont, calls, pad_bp=pad_bp, require_family=require_family)
    both = melt + [v for v in ont if v.key not in {x.key for x in melt}]
    both_m = overlap_metrics(both, calls, pad_bp=pad_bp, require_family=require_family)
    return {
        "cutoff": cutoff_name,
        "value": cutoff_value,
        "n_calls": len(calls),
        "melt_recall": melt_m["recall"],
        "melt_recovered": melt_m["truth_recovered"],
        "melt_truth": melt_m["n_truth"],
        "ont_recall": ont_m["recall"],
        "ont_recovered": ont_m["truth_recovered"],
        "ont_truth": ont_m["n_truth"],
        "union_recall": both_m["recall"],
        "novel_calls": both_m["novel_calls"],
        "precision_vs_union": both_m["precision_vs_truth"],
    }


def _tier_keeps(tier: str, rule: str) -> bool:
    rank = _TIER_RANK.get(tier, 0)
    if rule == "all":
        return True
    if rule == "not_none":
        return rank > 0
    if rule == "provisional_or_higher":
        return rank >= 3
    if rule == "high_conf_two_sided":
        return rank >= 4
    return True


def _as_int(value: Any, default: int = 0) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def _as_float(value: Any, default: float | None = None) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan") if default is None else default
