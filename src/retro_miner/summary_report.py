"""Scientific Summary Report Generator for MEI Callset Quality Control.

Parses annotated VCF v4.3 callset files and produces structured JSON payloads
plus publication-ready Markdown reports summarizing genotype distributions,
subfamily breakdowns, TSD statistics, transduction events, and TPRT motif scores.

Literature Anchor: Gardner et al. (2017) Nucleic Acids Res 45(12):e108.
"""

from __future__ import annotations

import json
import re
from dataclasses import asdict, dataclass
from pathlib import Path


@dataclass(frozen=True)
class CallsetSummary:
    total_candidates: int
    genotype_counts: dict[str, int]
    subfamily_breakdown: dict[str, int]
    mean_vaf: float
    tsds_detected_count: int
    poly_a_tails_detected_count: int
    transductions_detected_count: int
    canonical_tprt_count: int


_VCF_INFO_RE = re.compile(r"([^=;]+)(?:=([^;]*))?")


def _parse_info(info: str) -> dict[str, str]:
    out: dict[str, str] = {}
    if not info or info == ".":
        return out
    for m in _VCF_INFO_RE.finditer(info):
        key = m.group(1)
        val = m.group(2) if m.group(2) is not None else ""
        out[key] = val
    return out


def _parse_format(format_str: str, sample_str: str) -> dict[str, str]:
    fields = format_str.split(":")
    values = sample_str.split(":")
    out: dict[str, str] = {}
    for f, v in zip(fields, values, strict=False):
        out[f] = v
    return out


def generate_scientific_summary_report(
    vcf_path: str | Path,
    output_md_path: str | Path,
) -> CallsetSummary:
    """Parse an annotated MEI VCF and generate a scientific summary report.

    Writes both a Markdown report and a companion JSON file whose path is
    derived from *output_md_path* by replacing the ``.md`` suffix with
    ``.json``.

    Args:
        vcf_path: Path to input VCF v4.3 file.
        output_md_path: Path to write the Markdown summary report.

    Returns:
        CallsetSummary with aggregated callset metrics.
    """
    vcf_path = Path(vcf_path)
    output_md_path = Path(output_md_path)
    output_md_path.parent.mkdir(parents=True, exist_ok=True)

    genotype_counts: dict[str, int] = {}
    subfamily_breakdown: dict[str, int] = {}
    vafs: list[float] = []
    tsds = 0
    poly_a = 0
    transductions = 0
    tprt = 0

    with vcf_path.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 8:
                continue
            info = _parse_info(parts[7])
            fmt = parts[8] if len(parts) > 8 else ""
            sample = parts[9] if len(parts) > 9 else ""

            subfamily = str(info.get("SUBFAM", info.get("MEI_TYPE", "UNKNOWN")))
            subfamily_breakdown[subfamily] = subfamily_breakdown.get(subfamily, 0) + 1

            tsd_seq = str(info.get("TSD", ""))
            if tsd_seq and tsd_seq != "NONE" and tsd_seq != ".":
                tsds += 1

            polya = str(info.get("POLYA", "0"))
            if polya == "1":
                poly_a += 1

            tr_type = str(info.get("TRANSDUCTION_TYPE", "NONE"))
            if tr_type != "NONE":
                transductions += 1

            tprt_score = float(info.get("TPRT_MOTIF_SCORE", 0.0))
            if tprt_score >= 0.7:
                tprt += 1

            if fmt and sample:
                sample_fields = _parse_format(fmt, sample)
                gt = str(sample_fields.get("GT", "./."))
                genotype_counts[gt] = genotype_counts.get(gt, 0) + 1
                vaf_str = str(sample_fields.get("VAF", "0.0"))
                try:
                    vafs.append(float(vaf_str))
                except ValueError:
                    pass

    mean_vaf = sum(vafs) / len(vafs) if vafs else 0.0
    summary = CallsetSummary(
        total_candidates=sum(genotype_counts.values()),
        genotype_counts=dict(sorted(genotype_counts.items())),
        subfamily_breakdown=dict(sorted(subfamily_breakdown.items())),
        mean_vaf=round(mean_vaf, 4),
        tsds_detected_count=tsds,
        poly_a_tails_detected_count=poly_a,
        transductions_detected_count=transductions,
        canonical_tprt_count=tprt,
    )

    md = _build_markdown(summary)
    output_md_path.write_text(md, encoding="utf-8")

    json_path = output_md_path.with_suffix(".json")
    json_path.write_text(json.dumps(asdict(summary), indent=2), encoding="utf-8")

    return summary


def _build_markdown(summary: CallsetSummary) -> str:
    lines = [
        "# MEI Callset Scientific Summary Report",
        "",
        "## Callset Overview",
        "",
        f"- **Total Candidates:** {summary.total_candidates}",
        f"- **Mean VAF:** {summary.mean_vaf:.4f}",
        "",
        "## Genotype Distribution",
        "",
        "| Genotype | Count |",
        "|----------|-------|",
    ]
    for gt, count in sorted(summary.genotype_counts.items()):
        lines.append(f"| {gt} | {count} |")
    lines.extend(
        [
            "",
            "## Subfamily Breakdown",
            "",
            "| Subfamily | Count |",
            "|-----------|-------|",
        ]
    )
    for sub, count in sorted(summary.subfamily_breakdown.items()):
        lines.append(f"| {sub} | {count} |")
    lines.extend(
        [
            "",
            "## Mechanistic Annotations",
            "",
            f"- **TSDs Detected:** {summary.tsds_detected_count}",
            f"- **Poly(A) Tails Detected:** {summary.poly_a_tails_detected_count}",
            f"- **3' Transductions Detected:** {summary.transductions_detected_count}",
            f"- **Canonical TPRT Motifs (score ≥ 0.7):** {summary.canonical_tprt_count}",
            "",
        ]
    )
    return "\n".join(lines)


def load_summary_json(json_path: str | Path) -> CallsetSummary:
    """Load a CallsetSummary from a JSON file produced by ``generate_scientific_summary_report``."""
    payload = json.loads(Path(json_path).read_text(encoding="utf-8"))
    return CallsetSummary(**payload)
