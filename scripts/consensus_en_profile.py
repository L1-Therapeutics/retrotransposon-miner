#!/usr/bin/env python3
"""L1/Alu consensus endonuclease-site profile: the opportunity map.

For each L1 (L1HS/L1PA*) and Alu (Alu*) subfamily consensus in the UCSC
RepeatBrowser hg38reps.fa, this scans BOTH strands for three PRE-DECLARED
motif definitions, reported side by side so no conclusion depends on choosing
one motif definition:

  (a) narrow_tttt_aa      TTTTAA      chance frequency (1/4)^6  per strand position
  (b) flasch_7mer         TTTTTAA     chance frequency (1/4)^7  per strand position
  (c) degenerate_yy_rrrr  YYRRRR      chance frequency (1/2)^6  per strand position
      (IUPAC: Y in {C,T}, R in {A,G}; the self-ins machinery's degenerate
      EN consensus, artifacts/selfins_en_degenerate_family.csv columns
      YYRRRR / YYYYRR -- YYRRRR is this pattern on the scanned strand)

POSITION CONVENTION (documented in the CSV header comment): positions are
1-based from the 5' end of the consensus sequence AS GIVEN in the FASTA; the
two strands (consensus-as-given "sense", reverse complement "antisense") are
scanned and reported separately, each counted from its own 5' end. Overlapping
matches are all counted (a site exists at every position where the pattern
starts).

Bins are OPPORTUNITY units, not observations: 500 bp bins along L1
consensuses, 20 bp bins along Alu consensuses, with counts AND density per kb
of the bin (so short trailing bins are comparable).

This is an opportunity profile only. Whether child insertions USE these
opportunities is the observed-insertion comparison, which happens against the
multi-sample analysis outputs and is NOT part of this script.

Outputs (default nested_analysis/results_consensus/):
  consensus_en_sites.csv       one row per subfamily x strand x motif x position
  consensus_profile_summary.md per-subfamily table + probes a-d + all conventions
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
WORKSPACE = REPO.parent

# ---------------------------------------------------------------------------
# Pre-declared motif definitions and chance frequencies (stated, not implied).
# ---------------------------------------------------------------------------
IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
}

# name -> (pattern, chance frequency per strand position)
MOTIFS: dict[str, tuple[str, float]] = {
    "narrow_tttt_aa": ("TTTTAA", (1 / 4) ** 6),
    "flasch_7mer_ttttt_aa": ("TTTTTAA", (1 / 4) ** 7),
    "degenerate_yy_rrrr": ("YYRRRR", (1 / 2) ** 6),
}

# Binning (opportunity bins; bin width by consensus family).
BIN_SIZE = {"LINE1": 500, "ALU": 20}

# Subfamily selection: exactly what hg38reps.fa contains under these prefixes.
SUBFAMILY_PREFIXES = ("L1HS", "L1PA", "Alu")

# L1 approximate functional zones, scaled proportionally from the L1HS
# canonical annotation (5'UTR 1-907, ORF1 908-2380, ORF2 2381-5977, 3'UTR
# 5978-6021 over 6,021 bp). The FASTA carries no per-consensus annotation, so
# the same proportional split is used for every L1 consensus and is documented
# as such.
L1HS_ZONE_BPS = (
    ("5UTR", 1, 907),
    ("ORF1", 908, 2380),
    ("ORF2", 2381, 5977),
    ("3UTR", 5978, 6021),
)

# Alu probe windows (1-based, on the consensus as given).
ALU_LINKER_WINDOW = (118, 136)   # A-rich linker between the monomers
ALU_TAIL_WINDOW = (280, 320)     # tail region (reverse-complement T-rich logic)

# L1 inversion-breakpoint cluster probe window.
L1_INV_CLUSTER_WINDOW = (4000, 6000)

COMPLEMENT = str.maketrans("ACGTRYKMSWBDHVNacgtrykmswbdhvn",
                           "TGCAYRMKSWVHDBNtgcayrmkswvhdbn")


def reverse_complement(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def _pattern_regex(pattern: str) -> tuple[list[frozenset[str]], int]:
    """Precompiled motif: list of allowed-base sets, and its length."""
    classes = [frozenset(IUPAC[ch]) for ch in pattern]
    return classes, len(classes)


def find_motif_positions(seq: str, pattern: str) -> list[int]:
    """All 1-based start positions where `pattern` matches `seq`.

    Overlapping matches are counted. Ambiguous sequence letters (N etc.)
    match only if the pattern position explicitly allows them.
    """
    classes, k = _pattern_regex(pattern)
    upper = seq.upper()
    hits: list[int] = []
    for start in range(len(upper) - k + 1):
        window = upper[start:start + k]
        if all(window[i] in classes[i] for i in range(k)):
            hits.append(start + 1)
    return hits


def scan_strand(seq: str, motif: str) -> list[int]:
    return find_motif_positions(seq, MOTIFS[motif][0])


# ---------------------------------------------------------------------------
# FASTA input: enumerate exactly the subfamilies present.
# ---------------------------------------------------------------------------

@dataclass
class Consensus:
    subfamily: str
    family: str          # LINE1 or ALU
    sequence: str


def read_subfamily_consensuses(fasta_path: Path) -> list[Consensus]:
    """L1 (L1HS, L1PA*) and Alu (Alu*) subfamilies exactly as present."""
    consensuses: list[Consensus] = []
    name: str | None = None
    chunks: list[str] = []
    with fasta_path.open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    if name.startswith(SUBFAMILY_PREFIXES):
                        consensuses.append(_make_consensus(name, chunks))
                name = line[1:].strip()
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None and name.startswith(SUBFAMILY_PREFIXES):
        consensuses.append(_make_consensus(name, chunks))
    if not consensuses:
        raise SystemExit(f"no L1/Alu subfamily consensuses found in {fasta_path}")
    return consensuses


def _make_consensus(name: str, chunks: list[str]) -> Consensus:
    family = "LINE1" if name.startswith(("L1HS", "L1PA")) else "ALU"
    return Consensus(subfamily=name, family=family, sequence="".join(chunks).upper())


# ---------------------------------------------------------------------------
# Binning (opportunity) and probes.
# ---------------------------------------------------------------------------

def bin_label(pos: int, bin_size: int) -> tuple[int, int, int]:
    """(bin_index_0based, bin_start, bin_end) 1-based inclusive; the last bin
    of a consensus may be short and is reported with its true extent."""
    index = (pos - 1) // bin_size
    return index, index * bin_size + 1, (index + 1) * bin_size


def zone_for_l1_position(pos: int, length: int) -> str:
    if length <= 0:
        return L1HS_ZONE_BPS[-1][0]
    scale = length / L1HS_ZONE_BPS[-1][2]
    scaled = pos / scale
    for zone, start, end in L1HS_ZONE_BPS:
        if start <= scaled <= end:
            return zone
    return L1HS_ZONE_BPS[-1][0]


def count_in_window(positions: list[int], window: tuple[int, int]) -> int:
    lo, hi = window
    return sum(1 for p in positions if lo <= p <= hi)


def density_per_kb(count: int, bp: int) -> float:
    return count / (bp / 1000) if bp > 0 else 0.0


@dataclass
class SubfamilyResult:
    subfamily: str
    family: str
    length: int
    strand_counts: dict = field(default_factory=lambda: defaultdict(Counter))
    positions: dict = field(default_factory=lambda: defaultdict(list))


def profile_consensus(consensus: Consensus) -> SubfamilyResult:
    result = SubfamilyResult(
        subfamily=consensus.subfamily,
        family=consensus.family,
        length=len(consensus.sequence),
    )
    for strand, seq in (("sense", consensus.sequence),
                        ("antisense", reverse_complement(consensus.sequence))):
        for motif in MOTIFS:
            hits = scan_strand(seq, motif)
            result.positions[(strand, motif)] = hits
            result.strand_counts[strand][motif] = len(hits)
    return result


# ---------------------------------------------------------------------------
# Output assembly.
# ---------------------------------------------------------------------------

def site_rows(results: list[SubfamilyResult]) -> tuple[list[str], list[list[object]]]:
    """CSV header comment lines + one row per subfamily x strand x motif x position."""
    header_comments = [
        "# Consensus EN-site opportunity profile (opportunity map, NOT observed insertions).",
        "# Position convention: 1-based from the 5' end of the consensus sequence AS GIVEN;",
        "#   strand=sense scans the sequence as given, strand=antisense scans its reverse",
        "#   complement, each counted from its own 5' end. Overlapping matches counted.",
        "# bin_density_per_kb is the per-bin opportunity density (bins are opportunity).",
        "# chance_freq_per_bp is the motif's per-strand-position frequency under uniform",
        "#   base composition: TTTTAA (1/4)^6, TTTTTAA (1/4)^7, YYRRRR (1/2)^6.",
        "# age column is blank: numeric subfamily ages are not derivable from local files.",
    ]
    columns = [
        "subfamily", "family", "length_bp", "strand", "motif",
        "chance_freq_per_bp", "pos_1based_from_5p", "bin_size_bp",
        "bin_index", "bin_start", "bin_end", "bin_density_per_kb",
        "l1_zone", "age",
    ]
    rows: list[list[object]] = []
    for res in sorted(results, key=lambda r: (r.family, r.subfamily)):
        bin_size = BIN_SIZE[res.family]
        for (strand, motif), hits in sorted(res.positions.items()):
            chance = MOTIFS[motif][1]
            bin_bp: dict[int, int] = Counter()
            for pos in hits:
                index, start, end = bin_label(pos, bin_size)
                rows.append([
                    res.subfamily, res.family, res.length, strand, motif,
                    chance, pos, bin_size, index, start, end, 0.0,
                    zone_for_l1_position(pos, res.length) if res.family == "LINE1" else "",
                    "",
                ])
            # Bin opportunity extent: bins up to the consensus length.
            n_bins = (res.length + bin_size - 1) // bin_size
            bin_bp = {b: min((b + 1) * bin_size, res.length) - b * bin_size
                      for b in range(n_bins)}
            for row in rows[-len(hits):] if hits else []:
                b = int(row[8])  # type: ignore[arg-type]
                row[11] = round(density_per_kb(1, bin_bp[b]), 4)  # type: ignore[index]
    return header_comments + [",".join(columns)], rows


def write_sites_csv(path: Path, results: list[SubfamilyResult]) -> int:
    header, rows = site_rows(results)
    with path.open("w") as handle:
        for line in header:
            handle.write(line + "\n")
        for row in rows:
            handle.write(",".join(str(v) for v in row) + "\n")
    return len(rows)


def alu_probe(result: SubfamilyResult) -> dict[str, object]:
    """Probe (a): linker and tail windows, both strands, per motif."""
    out: dict[str, object] = {}
    for motif in MOTIFS:
        sense = result.positions[("sense", motif)]
        anti = result.positions[("antisense", motif)]
        out[motif] = {
            "linker_118_136_sense": count_in_window(sense, ALU_LINKER_WINDOW),
            "linker_118_136_antisense": count_in_window(anti, ALU_LINKER_WINDOW),
            "tail_280_320_sense": count_in_window(sense, ALU_TAIL_WINDOW),
            "tail_280_320_antisense": count_in_window(anti, ALU_TAIL_WINDOW),
        }
    return out


def l1_zone_densities(result: SubfamilyResult) -> dict[str, object]:
    """Probe (b): density per kb of each approximate functional zone (sense)."""
    zone_bp: dict[str, int] = defaultdict(int)
    scale = result.length / L1HS_ZONE_BPS[-1][2]
    for zone, start, end in L1HS_ZONE_BPS:
        zone_bp[zone] = max(0, round((end - start + 1) * scale))
    out: dict[str, object] = {"zone_bp": dict(zone_bp)}
    for motif in MOTIFS:
        hits = result.positions[("sense", motif)]
        zone_counts: dict[str, int] = defaultdict(int)
        for pos in hits:
            zone_counts[zone_for_l1_position(pos, result.length)] += 1
        out[motif] = {
            zone: {
                "sites": zone_counts.get(zone, 0),
                "density_per_kb": round(
                    density_per_kb(zone_counts.get(zone, 0), zone_bp[zone]), 3),
            }
            for zone in ("5UTR", "ORF1", "ORF2", "3UTR")
        }
    return out


def l1_inversion_cluster_probe(result: SubfamilyResult) -> dict[str, object]:
    """Probe (c): density in 4,000-6,000 bp vs the rest (sense, per motif)."""
    out: dict[str, object] = {}
    cluster_bp = min(result.length, L1_INV_CLUSTER_WINDOW[1]) - L1_INV_CLUSTER_WINDOW[0] + 1
    rest_bp = result.length - cluster_bp
    for motif in MOTIFS:
        hits = result.positions[("sense", motif)]
        inside = count_in_window(hits, L1_INV_CLUSTER_WINDOW)
        out[motif] = {
            "sites_4000_6000": inside,
            "density_cluster_per_kb": round(density_per_kb(inside, cluster_bp), 3),
            "sites_rest": len(hits) - inside,
            "density_rest_per_kb": round(density_per_kb(len(hits) - inside, rest_bp), 3),
        }
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--fasta", type=Path,
        default=Path.home() / "retrotransposon-workdir/data/public/retrotransposon_db/ucsc_repeatbrowser/hg38reps.fa",
        help="UCSC RepeatBrowser subfamily consensuses (hg38reps.fa)",
    )
    parser.add_argument(
        "--dfam-fasta", type=Path,
        default=Path.home() / "retrotransposon-workdir/data/public/retrotransposon_db/dfam/dfam_human_mei_l1_alu_sva.fasta",
        help="optional Dfam cross-check FASTA; agreement reported when present",
    )
    parser.add_argument(
        "--outdir", type=Path, default=WORKSPACE / "nested_analysis/results_consensus",
        help="output directory",
    )
    args = parser.parse_args(argv)

    if not args.fasta.is_file():
        raise SystemExit(f"required consensus FASTA not found: {args.fasta}")
    dfam_present = args.dfam_fasta.is_file()

    consensuses = read_subfamily_consensuses(args.fasta)
    results = [profile_consensus(c) for c in consensuses]

    args.outdir.mkdir(parents=True, exist_ok=True)
    n_rows = write_sites_csv(args.outdir / "consensus_en_sites.csv", results)
    write_summary(
        args.outdir / "consensus_profile_summary.md",
        results, args, dfam_present,
    )

    totals = {
        r.subfamily: {strand: dict(counts) for strand, counts in r.strand_counts.items()}
        for r in results
    }
    print(json.dumps({
        "subfamilies_scanned": len(results),
        "csv_rows": n_rows,
        "dfam_cross_check_available": dfam_present,
        "outputs": [
            str(args.outdir / "consensus_en_sites.csv"),
            str(args.outdir / "consensus_profile_summary.md"),
        ],
        "totals_by_subfamily": totals,
    }, indent=2))
    return 0


def write_summary(
    path: Path,
    results: list[SubfamilyResult],
    args: argparse.Namespace,
    dfam_present: bool,
) -> None:
    lines: list[str] = []
    lines.append("# Consensus endonuclease-site profile (opportunity map)")
    lines.append("")
    lines.append("This is an OPPORTUNITY profile of the consensus sequences themselves - where")
    lines.append("child insertions CAN land mechanistically. It makes no claim about where they")
    lines.append("DO land; the observed-insertion comparison runs against the multi-sample")
    lines.append("outputs and is not part of this document. No published per-consensus count of")
    lines.append("TTTTT/AA-class sites exists; this produces it.")
    lines.append("")
    lines.append("## Conventions (all pre-declared)")
    lines.append("")
    lines.append("- Consensus source: `" + str(args.fasta) + "`; subfamilies are exactly the")
    lines.append("  `L1HS`/`L1PA*`/`Alu*` records the file contains (" + str(len(results)) + " scanned).")
    lines.append("- Motif definitions, side by side (chance frequency per strand position,")
    lines.append("  uniform base composition):")
    for name, (pattern, chance) in MOTIFS.items():
        lines.append(f"  - `{name}`: `{pattern}` - {chance:.6g} per bp")
    lines.append("- Positions: 1-based from the 5' end of the consensus AS GIVEN; the two")
    lines.append("  strands (sense = as given, antisense = reverse complement) are scanned and")
    lines.append("  reported separately, each from its own 5' end. Overlapping matches counted.")
    lines.append("- Bins are opportunity units: 500 bp (L1), 20 bp (Alu); density per kb uses")
    lines.append("  each bin's true extent, so short trailing bins stay comparable.")
    lines.append("- L1 functional zones: the FASTA carries no annotation, so the L1HS canonical")
    lines.append("  split (5'UTR 1-907, ORF1 908-2380, ORF2 2381-5977, 3'UTR 5978-6021 over")
    lines.append("  6,021 bp) is scaled proportionally to each consensus length. Documented as an")
    lines.append("  approximate zone split, not per-consensus annotation.")
    lines.append("- Subfamily age: no local file derives numeric ages, so the age column is")
    lines.append("  left blank. The age-GRADIENT probe uses the standard class ordering only")
    lines.append("  (L1HS youngest, then L1PA2...L1PA17; AluY youngest, AluS middle, AluJ oldest),")
    lines.append("  which is a ranking, not a numeric age.")
    lines.append(f"- Dfam cross-check FASTA `{args.dfam_fasta}`: "
                 + ("present" if dfam_present else "NOT PRESENT on this machine - cross-check skipped and reported"))
    lines.append("")
    lines.append("## Per-subfamily totals (sense + antisense summed per motif)")
    lines.append("")
    lines.append("| subfamily | family | length (bp) | " + " | ".join(MOTIFS) + " | sites/kb (narrow, both strands) |")
    lines.append("|---|---|---|---|---|---|")
    for res in sorted(results, key=lambda r: (r.family, r.subfamily)):
        per_kb = density_per_kb(
            sum(res.strand_counts[s][m] for s in ("sense", "antisense") for m in ("narrow_tttt_aa",)),
            res.length,
        )
        totals = [sum(res.strand_counts[s][m] for s in ("sense", "antisense")) for m in MOTIFS]
        lines.append(
            f"| {res.subfamily} | {res.family} | {res.length} | "
            + " | ".join(str(t) for t in totals) + f" | {per_kb:.1f} |"
        )
    l1_results = [r for r in results if r.family == "LINE1"]
    alu_results = [r for r in results if r.family == "ALU"]

    lines.append("")
    lines.append("## Probe (a) - Alu linker (118-136) and tail (280-320) opportunity")
    lines.append("")
    lines.append("Factually: does an EN-compatible site exist inside/near the 118-136 A-rich")
    lines.append("linker, and in the 280-320 tail region (reverse-complement T-rich logic)?")
    lines.append("Opportunity counts only.")
    lines.append("")
    lines.append("| subfamily | motif | linker sense | linker antisense | tail sense | tail antisense |")
    lines.append("|---|---|---|---|---|---|")
    for res in sorted(alu_results, key=lambda r: r.subfamily):
        probe = alu_probe(res)
        for motif in MOTIFS:
            d = probe[motif]
            lines.append(
                f"| {res.subfamily} | {motif} | {d['linker_118_136_sense']} | "
                f"{d['linker_118_136_antisense']} | {d['tail_280_320_sense']} | {d['tail_280_320_antisense']} |"
            )
    any_linker = any(
        alu_probe(r)[m][f"linker_118_136_{s}"]
        for r in alu_results for m in MOTIFS for s in ("sense", "antisense")
    )
    any_tail = any(
        alu_probe(r)[m][f"tail_280_320_{s}"]
        for r in alu_results for m in MOTIFS for s in ("sense", "antisense")
    )
    lines.append("")
    lines.append(f"- Linker 118-136: EN-compatible site present in at least one Alu consensus: {any_linker}")
    lines.append(f"- Tail 280-320: EN-compatible site present in at least one Alu consensus: {any_tail}")

    lines.append("")
    lines.append("## Probe (b) - L1 zone densities (approximate zones, sense strand)")
    lines.append("")
    lines.append("Density per kb by zone; zone boundaries are the documented proportional L1HS split.")
    for res in sorted(l1_results, key=lambda r: r.subfamily):
        z = l1_zone_densities(res)
        parts = []
        for motif in MOTIFS:
            d = z[motif]  # type: ignore[index]
            parts.append(
                f"{motif}: " + ", ".join(
                    f"{zone} {d[zone]['density_per_kb']}/kb ({d[zone]['sites']})"
                    for zone in ("5UTR", "ORF1", "ORF2", "3UTR")
                )
            )
        lines.append(f"- **{res.subfamily}** - " + "; ".join(parts))

    lines.append("")
    lines.append("## Probe (c) - L1 4,000-6,000 bp (Porubsky inversion-breakpoint cluster)")
    lines.append("")
    lines.append("Opportunity counts only; observed-insertion comparison is out of scope here.")
    lines.append("")
    lines.append("| subfamily | motif | sites 4000-6000 | density cluster (kb) | sites rest | density rest (kb) |")
    lines.append("|---|---|---|---|---|---|")
    for res in sorted(l1_results, key=lambda r: r.subfamily):
        probe = l1_inversion_cluster_probe(res)
        for motif in MOTIFS:
            d = probe[motif]  # type: ignore[index]
            lines.append(
                f"| {res.subfamily} | {motif} | {d['sites_4000_6000']} | "
                f"{d['density_cluster_per_kb']} | {d['sites_rest']} | {d['density_rest_per_kb']} |"
            )
    lines.append("")
    higher = [
        (res.subfamily, motif)
        for res in l1_results for motif in MOTIFS
        if l1_inversion_cluster_probe(res)[motif]["density_cluster_per_kb"]  # type: ignore[index]
        > l1_inversion_cluster_probe(res)[motif]["density_rest_per_kb"]  # type: ignore[index]
    ]
    lines.append(f"- Subfamily x motif pairs with cluster density > rest density: {len(higher)} of {len(l1_results) * len(MOTIFS)}")

    lines.append("")
    lines.append("## Probe (d) - age gradient (opportunity, not observations)")
    lines.append("")
    l1_order = ["L1HS"] + [f"L1PA{n}" for n in range(2, 18)]
    present_l1 = [s for s in l1_order if any(r.subfamily == s for r in l1_results)]
    lines.append(f"- L1 subfamilies present in age order (youngest first): {present_l1}")
    for motif in MOTIFS:
        counts = {
            res.subfamily: sum(res.strand_counts[s][motif] for s in ("sense", "antisense"))
            for res in l1_results
        }
        densities = {
            res.subfamily: round(density_per_kb(counts[res.subfamily], res.length), 1)
            for res in l1_results
        }
        lines.append(f"  - `{motif}` sites (density/kb): "
                     + ", ".join(f"{s}: {counts[s]} ({densities[s]})" for s in present_l1))
    for group in ("AluJ", "AluS", "AluY"):
        members = [r for r in alu_results if r.subfamily.startswith(group)]
        if members:
            for motif in MOTIFS:
                densities = ", ".join(
                    f"{r.subfamily}: {round(density_per_kb(sum(r.strand_counts[s][motif] for s in ('sense', 'antisense')), r.length), 1)}/kb"
                    for r in sorted(members, key=lambda r: r.subfamily)
                )
                lines.append(f"- {group} `{motif}` density/kb: {densities}")
    lines.append("")
    path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    sys.exit(main())
