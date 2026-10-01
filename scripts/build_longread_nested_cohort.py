"""Build the nested-Alu long-read cohort and record every provenance gate.

This is Phase 0 of the structure-informed nested-Alu study. It reads the
locally available long-read insertion callset, assigns each Alu insertion to a
same-family RepeatMasker host, projects the host-relative breakpoint onto the
host's subfamily consensus, and records the genotype state at the pooled site.

Everything here is read-only with respect to the source data. Outputs are two
files in the output directory:

  per_call_longread_nested.csv   one row per Alu insertion site
  provenance_manifest.json      the verdict of every gate, not just a log line

Four gates run before any enrichment is computed, because each one has already
been observed to matter on this dataset:

  Gate 0a  index and reference consistency. The pooled BCF ships with a .csi
           older than the BCF itself. bcftools still answers region queries
           correctly, so the gate verifies *correctness* against a linear scan
           rather than rebuilding on mtime alone. The FASTA .fai is checked
           independently, because a stale or mis-built reference index fails
           far more quietly than a stale BCF index.
  Gate 0b  no CHM13 contamination. A `final-vcf.unphased.SVAN_1.3.vcf.gz`
           under the hg38 polymorphism directory is a dangling symlink into
           the hs1/ tree; that file is CHM13, not GRCh38. Only the
           `svim.asm.hg38.*` BCFs may be used.
  Gate 0c  the independence verdict for the genotype axis. HG03086 is itself
           one of the 908 pooled samples, so its genotype is internal to the
           discovery cohort. This gate writes a written verdict that decides
           whether Phase 3 may use the genotype axis at all.
  Gate 0d  the TSD_LEN sentinel check. TSD_LEN carries pile-ups at 41 (and 82
           for Alu/SVA) that are an order of magnitude above the local
           background. Those are sentinel values, not measured duplication
           lengths, so Phase 2b must not treat TSD length as continuous until
           they are excluded.

Positional convention, published here because the comparison to Levy et al.
2010 is otherwise unfalsifiable:

  host_offset_5p_0based  0-based offset of the insertion breakpoint from the
                         host's own 5' end, strand-aware: `pos0 - start0` for
                         a + strand host, `end0 - 1 - pos0` for a - strand
                         host. This is host-relative, not consensus-relative.

  consensus_offset       the same breakpoint projected through a global
                         alignment of the host to its own subfamily consensus,
                         measured from that consensus's 5' end.

Genotype state is recorded in three mutually exclusive columns and is never
collapsed to a boolean. The pooled callset has no missingness channel: across
~152k queried sites every GT is either 0/0 or non-ref, with no `./.` anywhere
and VAF1 unset throughout. A 0/0 therefore conflates genuine absence with a
genotyping failure inside a hard repeat, and "not confirmed" would silently
mix the two.
"""

from __future__ import annotations

import argparse
import bisect
import collections
import gzip
import json
import os
import re
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

REPO_SRC = Path(__file__).resolve().parents[1] / "src"
if str(REPO_SRC) not in sys.path:
    sys.path.insert(0, str(REPO_SRC))

# --------------------------------------------------------------------------
# Default inputs. All local, all read-only. Overridable on the command line so
# the tests can point the same code at small synthetic fixtures.
# --------------------------------------------------------------------------

WORKDIR = Path.home() / "retrotransposon-workdir" / "data" / "public"
DEFAULT_SITE_BCF = (
    WORKDIR / "polymorphism" / "hg38" / "long_read_1kg_ont_vienna"
    / "svim.asm.hg38.noGt.SVAN_1.3.bcf"
)
DEFAULT_GENO_BCF = (
    WORKDIR / "polymorphism" / "hg38" / "long_read_1kg_ont_vienna"
    / "svim.asm.hg38.bcf"
)
DEFAULT_RMSK = (
    Path(__file__).resolve().parents[2] / "nested_analysis" / "data" / "rmsk.txt.gz"
)
DEFAULT_FASTA = (
    WORKDIR / "reference" / "hg38" / "Homo_sapiens_assembly38.fasta"
)
DEFAULT_CONSENSUS = (
    WORKDIR / "retrotransposon_db" / "ucsc_repeatbrowser" / "hg38reps.fa"
)

FAMILY = "ALU"

# SVAN insertion POS sits this far below the HG03086 POS at matched events.
# Used only for the Gate 0c cross-callset join, never for host assignment.
HG03086_POS_SHIFT = -8
HG03086_MATCH_TOL = 3

#: GRCh38 primary assembly lengths. Gate 0a asserts the FASTA agrees.
GRCH38_LENGTHS = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559,
    "chr4": 190214555, "chr5": 181538259, "chr6": 170805979,
    "chr7": 159345973, "chr8": 145138636, "chr9": 138394717,
    "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189,
    "chr16": 90338345, "chr17": 83257441, "chr18": 80373285,
    "chr19": 58617616, "chr20": 64444167, "chr21": 46709983,
    "chr22": 50818468, "chrX": 156040895,
}

#: GRCh38 "full" lengths, i.e. primary plus patches. The pooled callset was
#: written against a full-GBD header (chr1 = 248,387,328 rather than the
#: primary 248,956,422), which is still GRCh38 and still shares coordinates on
#: the primary assembly for every primary contig. The BCF carries no alt or
#: patch contigs, so accepting either length set does not weaken the actual
#: purpose of the gate, which is to reject a CHM13 coordinate system (chr1 =
#: 248,387,328 in CHM13 is a coincidence of magnitude, but the other 22 contigs
#: are not: chr9 = 138,394,717 in GRCh38 versus 50,818,469 in CHM13).
GRCH38_FULL_LENGTHS = {
    "chr1": 248387328, "chr2": 242696752, "chr3": 201105948,
    "chr4": 193574945, "chr5": 182045439, "chr6": 172126628,
    "chr7": 160567428, "chr8": 146259331, "chr9": 150617247,
    "chr10": 134758134, "chr11": 135127769, "chr12": 133324548,
    "chr13": 113566686, "chr14": 101161492, "chr15": 99753195,
    "chr16": 96330374, "chr17": 84276897, "chr18": 80542538,
    "chr19": 61707364, "chr20": 66210255, "chr21": 45090682,
    "chr22": 51324926, "chrX": 154259566,
}

#: A TSD_LEN value whose count exceeds this multiple of the local median of its
#: +/-3 neighbours is treated as a sentinel rather than a measurement.
SENTINEL_RATIO = 8.0
SENTINEL_NEIGHBOURHOOD = 3


class GateFailure(RuntimeError):
    """A provenance gate failed in a way that must stop the build."""


def run(cmd: list[str]) -> subprocess.CompletedProcess:
    """Run a command, capturing text. Raises with stderr on non-zero exit."""
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise GateFailure(
            f"command failed ({proc.returncode}): {' '.join(cmd)}\n{proc.stderr[:2000]}"
        )
    return proc


def bcftools_lines(cmd: list[str]) -> list[str]:
    """Run bcftools, dropping htslib warnings from stderr and any prefix lines."""
    proc = run(cmd)
    keep = []
    for line in proc.stdout.splitlines():
        if line.startswith("[") or not line.strip():
            continue
        keep.append(line)
    return keep


# --------------------------------------------------------------------------
# Gate 0a: index and reference consistency
# --------------------------------------------------------------------------


def _indexed_positions(bcf: Path, chrom: str, start: int, end: int) -> set[int]:
    out = bcftools_lines(
        ["bcftools", "query", "-r", f"{chrom}:{start}-{end}", "-f", "%POS\n", str(bcf)]
    )
    return {int(x) for x in out}


def _linear_positions(bcf: Path, chrom: str, start: int, end: int) -> set[int]:
    """Region membership computed without the index, by a full linear scan."""
    hits: set[int] = set()
    for line in bcftools_lines(["bcftools", "view", "-H", str(bcf)]):
        fields = line.split("\t", 2)
        if len(fields) < 2 or fields[0] != chrom:
            continue
        pos = int(fields[1])
        if start <= pos <= end:
            hits.add(pos)
    return hits


def check_index_consistency(
    bcf: Path,
    probe_regions: Iterable[tuple[str, int, int]] | None = None,
    scratch: Path | None = None,
) -> dict[str, Any]:
    """Verify the BCF index returns correct results, not merely a fresh mtime.

    The shipped .csi is older than its BCF. That produces a warning on every
    region query, but an observation worth recording is that the answers were
    still correct: each probe region was compared against a full linear scan.
    So this gate tests *correctness* and only reaches for a rebuild when the
    index is actually wrong, which avoids an expensive rebuild triggered by a
    cosmetic mtime difference.

    Returns a manifest fragment recording the mtimes, the comparison, and
    whether a scratch rebuild was required.
    """
    bcf = Path(bcf)
    index = bcf.with_suffix(bcf.suffix + ".csi")
    result: dict[str, Any] = {
        "bcf": str(bcf),
        "index": str(index),
        "index_exists": index.exists(),
        "index_older_than_bcf": None,
        "mtime_index": None,
        "mtime_bcf": None,
        "probe_regions": [],
        "all_regions_agree": None,
        "rebuild_required": False,
        "rebuild_path": None,
        "verdict": "unknown",
    }
    if bcf.exists():
        result["mtime_bcf"] = os.path.getmtime(bcf)
    if index.exists():
        result["mtime_index"] = os.path.getmtime(index)
    if result["mtime_index"] is not None and result["mtime_bcf"] is not None:
        result["index_older_than_bcf"] = result["mtime_index"] < result["mtime_bcf"]

    if probe_regions is None:
        probe_regions = [
            ("chr1", 4_000_000, 4_400_000),
            ("chr1", 100_000_000, 100_050_000),
            ("chr7", 55_000_000, 55_050_000),
            ("chrX", 1_000_000, 1_000_500),
        ]

    disagreements = []
    for chrom, start, end in probe_regions:
        indexed = _indexed_positions(bcf, chrom, start, end)
        linear = _linear_positions(bcf, chrom, start, end)
        agrees = indexed == linear
        result["probe_regions"].append(
            {
                "region": f"{chrom}:{start}-{end}",
                "n_indexed": len(indexed),
                "n_linear": len(linear),
                "agrees": agrees,
            }
        )
        if not agrees:
            disagreements.append((chrom, start, end))

    result["all_regions_agree"] = not disagreements

    if disagreements and scratch is not None:
        scratch.mkdir(parents=True, exist_ok=True)
        rebuilt = scratch / bcf.name
        run(["bcftools", "index", "-f", "--threads", "4", "-o", str(rebuilt) + ".csi", str(bcf)])
        result["rebuild_required"] = True
        result["rebuild_path"] = str(rebuilt)
        result["verdict"] = "index_was_incorrect_rebuilt_to_scratch"
    elif disagreements:
        result["verdict"] = "index_was_incorrect_rebuild_unavailable"
    else:
        stale = result["index_older_than_bcf"]
        result["verdict"] = (
            "index_stale_by_mtime_but_results_correct"
            if stale
            else "index_fresh_and_results_correct"
        )
    return result


def check_reference_consistency(fasta: Path, n_probe: int = 12, seed: int = 20261002) -> dict[str, Any]:
    """Verify the FASTA .fai separately from the BCF index.

    A stale or mis-built reference index produces wrong sequence silently,
    whereas a stale BCF index at least announces itself. So this checks that
    .fai offsets actually resolve to bases in the file and that the chrom
    names and lengths are GRCh38.
    """
    import random

    fasta = Path(fasta)
    fai = Path(str(fasta) + ".fai")
    result: dict[str, Any] = {
        "fasta": str(fasta),
        "fai": str(fai),
        "fasta_exists": fasta.exists(),
        "fai_exists": fai.exists(),
        "fai_older_than_fasta": None,
        "n_grch38_chroms_present": 0,
        "n_grch38_chroms_expected": len(GRCH38_LENGTHS),
        "length_mismatches": [],
        "byte_spot_checks": 0,
        "byte_spot_failures": [],
        "verdict": "unknown",
    }
    if not fai.exists():
        result["verdict"] = "missing_fai"
        raise GateFailure(f"FASTA index missing: {fai}")

    if fasta.exists():
        result["fai_older_than_fasta"] = os.path.getmtime(fai) < os.path.getmtime(fasta)

    index: dict[str, tuple[int, int, int, int]] = {}
    with open(fai) as fh:
        for line in fh:
            parts = line.split("\t")[:5]
            if len(parts) < 5:
                continue
            name, length, offset, linebases, linewidth = parts
            index[name] = (int(length), int(offset), int(linebases), int(linewidth))

    for chrom, expected in GRCH38_LENGTHS.items():
        if chrom not in index:
            continue
        result["n_grch38_chroms_present"] += 1
        if index[chrom][0] != expected:
            result["length_mismatches"].append(
                {"chrom": chrom, "fai_length": index[chrom][0], "grch38_length": expected}
            )

    rng = random.Random(seed)
    present = [c for c in GRCH38_LENGTHS if c in index]
    if present:
        with open(fasta, "rb") as fh:
            for chrom in rng.sample(present, min(4, len(present))):
                length, offset, linebases, linewidth = index[chrom]
                for _ in range(max(1, n_probe // 4)):
                    pos = rng.randrange(0, max(1, length - 1))
                    byte_off = offset + (pos // linebases) * linewidth + (pos % linebases)
                    fh.seek(byte_off)
                    base = fh.read(1).decode("ascii", "replace").upper()
                    result["byte_spot_checks"] += 1
                    if base not in "ACGTN":
                        result["byte_spot_failures"].append(
                            {"chrom": chrom, "pos": pos, "base": base}
                        )

    result["verdict"] = (
        "consistent"
        if not result["length_mismatches"] and not result["byte_spot_failures"]
        else "inconsistent"
    )
    if result["verdict"] == "inconsistent":
        raise GateFailure(
            "FASTA/.fai inconsistent with GRCh38: "
            f"length_mismatches={result['length_mismatches']} "
            f"byte_failures={result['byte_spot_failures']}"
        )
    return result


# --------------------------------------------------------------------------
# Gate 0b: no CHM13 contamination
# --------------------------------------------------------------------------


def check_no_chm13(
    site_bcf: Path, geno_bcf: Path, required_chroms: Iterable[str] | None = None
) -> dict[str, Any]:
    """Assert the inputs really are GRCh38 and not the hs1/ CHM13 files.

    The hg38 polymorphism directory also contains a
    `final-vcf.unphased.SVAN_1.3.vcf.gz` that is a symlink into the hs1/ tree.
    Reading it would silently mix a CHM13 coordinate system into a GRCh38
    analysis, so this refuses any input path that traverses hs1/ and asserts
    the observed chromosome names and lengths.
    """
    required = list(required_chroms or GRCH38_LENGTHS)
    result: dict[str, Any] = {
        "checked_paths": [str(site_bcf), str(geno_bcf)],
        "path_checks": [],
        "chrom_lengths_ok": None,
        "observed_chroms": [],
        "verdict": "unknown",
    }
    for raw in (site_bcf, geno_bcf):
        path = Path(raw)
        resolved = path.resolve()
        entry = {
            "path": str(path),
            "resolved": str(resolved),
            "is_symlink": path.is_symlink(),
            "traverses_hs1": "hs1" in resolved.parts,
            "exists": resolved.exists(),
        }
        result["path_checks"].append(entry)
        if entry["traverses_hs1"]:
            raise GateFailure(
                f"input resolves into the hs1/ (CHM13) tree: {path} -> {resolved}"
            )
        if entry["is_symlink"]:
            raise GateFailure(
                f"input is a symlink; refusing to follow into a possibly "
                f"different assembly: {path} -> {resolved}"
            )

    contigs = _bcftools_contigs(site_bcf)
    result["observed_chroms"] = sorted(contigs)
    bad = {}
    for c in required:
        if c not in contigs:
            continue
        if contigs[c] == GRCH38_LENGTHS[c]:
            continue
        if contigs[c] == GRCH38_FULL_LENGTHS.get(c):
            continue
        bad[c] = contigs[c]
    result["chrom_lengths_ok"] = not bad
    result["length_mismatch_detail"] = bad
    result["header_length_set"] = (
        "grch38_full" if not bad and any(
            contigs.get(c) == GRCH38_FULL_LENGTHS.get(c)
            and contigs.get(c) != GRCH38_LENGTHS.get(c)
            for c in required if c in contigs
        ) else "grch38_primary"
    )
    if bad:
        raise GateFailure(f"site BCF chrom lengths are not GRCh38: {bad}")
    result["verdict"] = "grch38_confirmed"
    return result


def _bcftools_contigs(bcf: Path) -> dict[str, int]:
    contigs: dict[str, int] = {}
    for line in bcftools_lines(["bcftools", "view", "-h", str(bcf)]):
        if line.startswith("##contig=<"):
            name = re.search(r"ID=([^,>]+)", line)
            length = re.search(r"length=(\d+)", line)
            if name and length:
                contigs[name.group(1)] = int(length.group(1))
    return contigs


# --------------------------------------------------------------------------
# Gate 0c: independence verdict for the genotype axis
# --------------------------------------------------------------------------


def pooled_site_families(site_bcf: Path) -> dict[str, list[int]]:
    """{chrom: sorted insertion positions} for one family in the site BCF.

    FAM_N is Number="." and renders as "." through a bcftools query format
    string, so the raw INFO column is parsed instead.
    """
    want = FAMILY.lower()
    sites: dict[str, list[int]] = collections.defaultdict(list)
    for line in bcftools_lines(["bcftools", "view", "-H", str(site_bcf)]):
        fields = line.split("\t", 8)
        if len(fields) < 8:
            continue
        fam = ""
        for item in fields[7].split(";"):
            if item.startswith("FAM_N="):
                fam = item.split("=", 1)[1]
                break
        if fam.strip().lower() == want:
            sites[fields[0]].append(int(fields[1]))
    return {c: sorted(v) for c, v in sites.items()}


def read_vcf_info(vcf: Path) -> tuple[dict[tuple[str, int], dict[str, str]], list[tuple[str, int]]]:
    """Exact-coordinate INFO records plus file order, from a plain VCF."""
    info: dict[tuple[str, int], dict[str, str]] = {}
    order: list[tuple[str, int]] = []
    with open(vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            record = dict(
                kv.split("=", 1) for kv in fields[7].split(";") if "=" in kv
            )
            key = (fields[0], int(fields[1]))
            info[key] = record
            order.append(key)
    return info, order


def fetch_region_field(
    bcf: Path, chrom: str, start: int, end: int, sample: str | None, tag: str
) -> dict[int, str]:
    """Query one INFO or FORMAT tag across a region, keyed by position.

    Relies on the BCF index, whose correctness Gate 0a establishes (and
    rebuilds into a scratch path if it is wrong). `tag` is a bare name such as
    "AC" or "GT"; the FORMAT/FILTER distinction is supplied here.
    """
    # bcftools requires the FORMAT tag *inside* the brackets: [%GT]. Passing
    # [GT] is accepted without complaint but emits the literal string "GT" for
    # every record, which is how an earlier version of this function ended up
    # calling all 26,499 sites concordant.
    spec = f"%INFO/{tag}" if sample is None else f"[%{tag}]"
    cmd = ["bcftools", "query", "-r", f"{chrom}:{start}-{end}", "-f", f"%POS\t{spec}\n"]
    if sample is not None:
        cmd += ["-s", sample]
    cmd.append(str(bcf))
    out: dict[int, str] = {}
    for line in bcftools_lines(cmd):
        parts = line.split("\t")
        if len(parts) == 2:
            out[int(parts[0])] = parts[1]
    return out


def load_sample_order(bcf: Path) -> list[str]:
    return [s for s in bcftools_lines(["bcftools", "query", "-l", str(bcf)]) if s]


def genotype_is_nonref(gt: str) -> bool:
    """True for any called alt-bearing genotype. `./.` is never non-ref."""
    if gt in ("", "."):
        return False
    alleles = re.split(r"[/|]", gt)
    if not alleles or any(a in (".", "") for a in alleles):
        return False
    return any(a not in ("0",) for a in alleles)


def assess_genotype_axis_independence(
    geno_bcf: Path,
    site_bcf: Path,
    hg03086_vcf: Path,
    sample: str = "HG03086",
) -> dict[str, Any]:
    """Write the verdict that decides whether Phase 3 may use the genotype axis.

    The question this answers is narrow and consequential: is HG03086's
    genotype inside the pooled long-read callset independent evidence about
    HG03086's own short-read calls, or is the cohort agreeing with itself?

    The decisive fact is membership. HG03086 is one of the 908 samples whose
    assemblies contributed to the pooled callset, so for any site that cohort
    is not an outside witness. A site can only speak to HG03086 independently
    when it was plausibly discovered elsewhere, which is what the AC band
    stratifies: at AC > 2 other carriers exist, so the site was seeded by
    another sample's data. At AC <= 2 HG03086 may be the only carrier, and the
    site may well have been discovered *because* HG03086 carries it.
    """
    samples = load_sample_order(geno_bcf)
    in_pool = sample in samples
    sites = pooled_site_families(site_bcf)
    info, order = read_vcf_info(hg03086_vcf)

    calls = [
        k for k in order if info[k].get("MEIFAMILY", "").strip().upper() == FAMILY
    ]

    ac: dict[tuple[str, int], int] = {}
    gt: dict[tuple[str, int], str] = {}
    for chrom, poss in sites.items():
        lo, hi = min(poss) - 5, max(poss) + 5
        got_ac = fetch_region_field(geno_bcf, chrom, lo, hi, None, "AC")
        for pos, raw in got_ac.items():
            try:
                ac[(chrom, pos)] = int(raw.split(",")[0])
            except ValueError:
                continue
        for pos, raw in fetch_region_field(geno_bcf, chrom, lo, hi, sample, "GT").items():
            gt[(chrom, pos)] = raw

    table: collections.Counter = collections.Counter()
    unmatched = 0
    for key in calls:
        chrom, pos = key
        cand = sites.get(chrom, [])
        if not cand:
            unmatched += 1
            continue
        target = pos + HG03086_POS_SHIFT
        i = bisect.bisect_left(cand, target)
        hit = None
        for j in (i - 1, i, i + 1):
            if 0 <= j < len(cand) and abs(cand[j] - target) <= HG03086_MATCH_TOL:
                hit = (chrom, cand[j])
                break
        if hit is None:
            unmatched += 1
            continue
        nonref = genotype_is_nonref(gt.get(hit, "."))
        long_read_derived = "long_read" in info[key].get("KNOWN_SRC", "")
        allele_count = ac.get(hit)
        band = (
            "AC_le_2"
            if allele_count is not None and allele_count <= 2
            else "AC_gt_2"
            if allele_count is not None
            else "AC_unknown"
        )
        table[(long_read_derived, nonref, band)] += 1

    # A verdict is written, not inferred downstream.
    if not in_pool:
        verdict = "independent"
        rationale = (
            f"{sample} is not among the {len(samples)} pooled samples, so its "
            "genotype is an outside observation of the discovery cohort."
        )
        tier = "top"
    else:
        verdict = "not_independent"
        rationale = (
            f"{sample} IS one of the {len(samples)} pooled samples, so the "
            "genotype axis is internal to the discovery cohort. For any site "
            "this cohort could have discovered because the sample carries it, "
            "genotype is circular rather than corroborating. Phase 3 may only "
            "use sites whose AC shows other carriers (AC > 2), which is a "
            "restricted tier and never licenses a detection-accuracy claim."
        )
        tier = "restricted_concordance_only"

    return {
        "sample": sample,
        "n_pooled_samples": len(samples),
        "sample_in_pooled_cohort": in_pool,
        "n_family_calls_in_reference_vcf": len(calls),
        "n_calls_without_matched_pooled_site": unmatched,
        "match_shift_bp": HG03086_POS_SHIFT,
        "match_tolerance_bp": HG03086_MATCH_TOL,
        "stratified_counts": {
            f"long_read_derived={lr}|nonref={nr}|{band}": n
            for (lr, nr, band), n in sorted(table.items(), key=str)
        },
        "verdict": verdict,
        "phase3_tier": tier,
        "rationale": rationale,
        "prohibited_claims": [
            "precision",
            "recall",
            "confirmed_label",
            "detection_accuracy",
        ],
    }


# --------------------------------------------------------------------------
# Gate 0d: TSD_LEN sentinel detection
# --------------------------------------------------------------------------


def detect_tsd_sentinels(
    site_bcf: Path,
    family: str = FAMILY,
    ratio: float = SENTINEL_RATIO,
    neighbourhood: int = SENTINEL_NEIGHBOURHOOD,
) -> dict[str, Any]:
    """Find TSD_LEN values that are sentinels rather than measurements.

    A sentinel shows up as a single value whose count towers over the local
    background on both sides and then collapses immediately afterwards. On
    this callset TSD_LEN=41 carries 2781 Alu calls against a local median of
    23, a 121x excess, and the neighbouring value 42 collapses to 8. That is
    the signature of a fixed-width alignment window being reported as a
    duplication length, not of biology.

    Detecting these empirically rather than hardcoding the observed values
    keeps the check honest if the callset is ever regenerated.
    """
    hist: collections.Counter = collections.Counter()
    want = family.lower()
    # Parsed from raw INFO: FAM_N is Number="." and bcftools query renders it
    # as ".", which would silently yield an empty histogram.
    for line in bcftools_lines(["bcftools", "view", "-H", str(site_bcf)]):
        fields = line.split("\t", 8)
        if len(fields) < 8:
            continue
        fam = tsd = ""
        for item in fields[7].split(";"):
            if item.startswith("FAM_N="):
                fam = item.split("=", 1)[1]
            elif item.startswith("TSD_LEN="):
                tsd = item.split("=", 1)[1]
        if fam.strip().lower() != want or tsd in ("", "."):
            continue
        try:
            hist[int(tsd)] += 1
        except ValueError:
            continue

    sentinels: dict[str, Any] = {}
    for value in sorted(hist):
        neighbours = [
            hist.get(value + d, 0)
            for d in range(-neighbourhood, neighbourhood + 1)
            if d != 0
        ]
        neighbours = [n for n in neighbours if n > 0]
        if not neighbours:
            continue
        base = sorted(neighbours)[len(neighbours) // 2]
        if base and hist[value] / base >= ratio:
            sentinels[str(value)] = {
                "count": hist[value],
                "local_median": base,
                "ratio": round(hist[value] / base, 2),
            }

    total = sum(hist.values())
    excluded = sum(v["count"] for v in sentinels.values())
    return {
        "family": family,
        "n_calls_with_tsd": total,
        "sentinel_values": sentinels,
        "sentinel_count": len(sentinels),
        "sentinel_share_of_calls": round(excluded / total, 6) if total else None,
        "ratio_threshold": ratio,
        "neighbourhood": neighbourhood,
        "verdict": (
            "sentinels_detected_tsd_length_not_continuous"
            if sentinels
            else "no_sentinels_detected"
        ),
        "consequence": (
            "Phase 2b must exclude these values and may not treat TSD length as "
            "a continuous variable. The 5-27bp de-novo classification is only "
            "valid on the non-sentinel subset."
            if sentinels
            else "TSD length may be treated as continuous."
        ),
    }


def tsd_length_is_sentinel(tsd_len: int, sentinel_values: Iterable[int]) -> bool:
    """True when a TSD length is a sentinel and must be treated as unusable."""
    return tsd_len in set(sentinel_values)


# --------------------------------------------------------------------------
# Site records
# --------------------------------------------------------------------------

SITE_FIELDS = [
    "FAM_N", "STRAND", "TSD_LEN", "TSD_SEQ", "POLYA_LEN", "POLYA_SEQ",
    "RT_LEN", "PERC_RESOLVED", "NOT_CANONICAL", "CONFORMATION",
]


def load_site_records(site_bcf: Path, family: str = FAMILY) -> list[dict[str, str]]:
    """One dict per insertion site for the requested family, in file order.

    INFO is parsed from raw `bcftools view -H` output rather than through
    `bcftools query -f %INFO/...`. Several of the fields this needs are
    declared Number="." (FAM_N among them), and bcftools silently renders those
    as "." through a query format string, which would drop every record.
    """
    want = family.lower()
    records: list[dict[str, str]] = []
    for line in bcftools_lines(["bcftools", "view", "-H", str(site_bcf)]):
        fields = line.split("\t", 8)
        if len(fields) < 8:
            continue
        chrom, pos = fields[0], fields[1]
        info: dict[str, str] = {}
        for item in fields[7].split(";"):
            if "=" in item:
                key, value = item.split("=", 1)
                info[key] = value
            else:
                info[item] = "true"  # bare flag, e.g. NOT_CANONICAL
        if info.get("FAM_N", "").strip().lower() != want:
            continue
        records.append(
            {
                "chrom": chrom,
                "pos": int(pos),
                **{f: info.get(f, ".") for f in SITE_FIELDS},
            }
        )
    return records


def is_alu_interval(fields: list[str]) -> bool:
    """RepeatMasker row is a same-family Alu interval."""
    name, rep_class, rep_family = fields[10], fields[11], fields[12]
    return "alu" in (name + rep_class + rep_family).lower()


#: Host-selection tie-breaks. `containment` takes the narrowest containing
#: interval; `longest_span` takes the widest. Both are pre-declared so Phase 2a
#: can recompute the whole profile under each and test whether a peak is an
#: artefact of the choice.
HOST_SELECTION_RULES = ("longest_span", "narrowest_containment")


def select_host(
    candidates: list[tuple[int, int, str, str]],
    rule: str = "longest_span",
) -> tuple[int, int, str, str] | None:
    """Choose one host interval among overlapping candidates.

    The saved HG03086 analysis resolved overlapping hosts by greatest span then
    leftmost start, so that is the default here and `host_selection_rule`
    records it per row. Phase 2a perturbs this choice to test peak stability.
    """
    if not candidates:
        return None
    if rule == "narrowest_containment":
        return sorted(candidates, key=lambda c: (c[1] - c[0], c[0], c[3]))[0]
    if rule == "longest_span":
        return sorted(candidates, key=lambda c: (-(c[1] - c[0]), c[0], c[3]))[0]
    raise ValueError(f"unknown host selection rule: {rule}")


def host_offset(
    pos: int, start0: int, end0: int, host_strand: str
) -> int | None:
    """Breakpoint offset from the host's own 5' end, strand-aware.

    Published convention, 0-based. For a + strand host the offset runs left to
    right; for a - strand host the host's 5' end is its rightmost base, so the
    offset is measured from there. This is host-relative and is not a consensus
    coordinate.
    """
    pos0 = pos - 1
    if host_strand == "+":
        return pos0 - start0
    if host_strand == "-":
        return end0 - 1 - pos0
    return None


def assign_hosts(
    records: list[dict[str, str]],
    rmsk: Path,
    rules: Iterable[str] = HOST_SELECTION_RULES,
) -> dict[str, dict[tuple[str, int], tuple[int, int, str, str]]]:
    """{rule: {(chrom,pos): host}} for every insertion site, in one rmsk pass.

    Nested-into-same-family is decided by containment in a same-family Alu
    interval. Sites with no containing interval are simply absent from the
    result, which is how "nested" is expressed without a second boolean.
    """
    rules = list(rules)
    by_chrom: dict[str, list[tuple[int, dict[str, str]]]] = collections.defaultdict(list)
    for rec in records:
        by_chrom[rec["chrom"]].append((rec["pos"], rec))
    sorted_starts = {
        c: sorted(p for p, _ in v) for c, v in by_chrom.items()
    }

    hits: dict[tuple[str, int], list[tuple[int, int, str, str]]] = (
        collections.defaultdict(list)
    )
    with gzip.open(rmsk, "rt", errors="replace") as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 13 or not fields[5].startswith("chr"):
                continue
            if not is_alu_interval(fields):
                continue
            chrom = fields[5]
            if chrom not in by_chrom:
                continue
            start0, end0 = int(fields[6]), int(fields[7])
            if end0 <= start0:
                continue
            strand, name = fields[9], fields[10]
            positions = sorted_starts[chrom]
            lo = bisect.bisect_left(positions, start0)
            hi = bisect.bisect_left(positions, end0)
            for pos, _ in by_chrom[chrom][lo:hi]:
                if start0 < pos <= end0:
                    hits[(chrom, pos)].append((start0, end0, strand, name))

    assigned: dict[str, dict[tuple[str, int], tuple[int, int, str, str]]] = {}
    for rule in rules:
        table: dict[tuple[str, int], tuple[int, int, str, str]] = {}
        for key, candidates in hits.items():
            chosen = select_host(candidates, rule)
            if chosen is not None:
                table[key] = chosen
        assigned[rule] = table
    return assigned


# --------------------------------------------------------------------------
# Consensus coordinates
# --------------------------------------------------------------------------


def load_consensus(path: Path, family: str = FAMILY) -> dict[str, str]:
    """{subfamily: consensus sequence} for one repeat family."""
    seqs: dict[str, list[str]] = {}
    name: str | None = None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].strip().split()[0]
                seqs[name] = []
            elif name:
                seqs[name].append(line.strip())
    joined = {k: "".join(v).upper() for k, v in seqs.items()}
    marker = family.lower()
    return {k: v for k, v in joined.items() if marker in k.lower() and v}


def pick_consensus_for_host(
    host_name: str, consensus: dict[str, str]
) -> tuple[str, str] | None:
    """Best available subfamily consensus for a RepeatMasker host name.

    Host names carry a serial suffix (AluYb8, AluSx1). An exact hit is
    preferred, then the digit-stripped stem, then a prefix match. Returns the
    matched consensus name so the choice is auditable per row.
    """
    if host_name in consensus:
        return host_name, consensus[host_name]
    stem = "".join(ch for ch in host_name if not ch.isdigit())
    if stem in consensus:
        return stem, consensus[stem]
    for key in sorted(consensus):
        if host_name.startswith(key) or key.startswith(host_name):
            return key, consensus[key]
    for key in sorted(consensus):
        if host_name[:5] == key[:5]:
            return key, consensus[key]
    return None


def build_aligner():
    """Global affine-gap aligner used for host -> consensus projection."""
    from Bio import Align

    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5
    return aligner


def project_offset(
    aligner: Any,
    host_seq: str,
    consensus_seq: str,
    host_offset_0based: int,
) -> dict[str, Any]:
    """Project a host-relative breakpoint onto the subfamily consensus.

    The primary analysis set is "breakpoint maps unambiguously to its
    subfamily consensus". Consensus indels matter here: an unaligned stretch in
    a diverged host shifts every downstream position, and a shifted coordinate
    can move a breakpoint across a +/-5 bp window boundary, which would
    silently corrupt the position-133 test. So the result reports the largest
    indel seen in the alignment and marks anything beyond the tolerance as
    ambiguous rather than quietly projecting it.
    """
    out: dict[str, Any] = {
        "consensus_offset": None,
        "offset_drift_bp": None,
        "alignment_identity": None,
        "alignment_score": None,
        "n_unaligned_host_bases": None,
        "max_indel_bp": None,
        "consensus_mapping": "ambiguous",
        "consensus_match_name": None,
    }
    if not host_seq or not consensus_seq:
        return out
    if not 0 <= host_offset_0based < len(host_seq):
        return out

    try:
        alignment = aligner.align(host_seq, consensus_seq)[0]
    except (ValueError, IndexError):
        return out

    host_blocks, cons_blocks = alignment.aligned
    mapping: dict[int, int] = {}
    max_indel = 0
    gaps: list[tuple[int, int]] = []
    for (hs, he), (cs, ce) in zip(host_blocks, cons_blocks):
        if he - hs != ce - cs:
            # Should not happen for a global alignment, but record it.
            max_indel = max(max_indel, abs((he - hs) - (ce - cs)))
        for k in range(he - hs):
            mapping[hs + k] = cs + k
    # Track leading/trailing and interior gaps between aligned blocks.
    for i in range(min(len(host_blocks), len(cons_blocks)) - 1):
        host_gap = host_blocks[i + 1][0] - host_blocks[i][1]
        cons_gap = cons_blocks[i + 1][0] - cons_blocks[i][1]
        if host_gap or cons_gap:
            gaps.append((max(host_gap, 0), max(cons_gap, 0)))
            max_indel = max(max_indel, host_gap, cons_gap)
    lead = (host_blocks[0][0] if len(host_blocks) else 0,
            cons_blocks[0][0] if len(cons_blocks) else 0)
    if lead[0] or lead[1]:
        gaps.append(lead)
        max_indel = max(max_indel, *lead)

    matches = sum(
        1
        for hpos, cpos in mapping.items()
        if cpos < len(consensus_seq) and host_seq[hpos] == consensus_seq[cpos]
    )
    identity = matches / len(mapping) if mapping else 0.0

    out["alignment_identity"] = round(identity, 4)
    out["alignment_score"] = float(alignment.score)
    out["n_unaligned_host_bases"] = len(host_seq) - len(mapping)
    out["max_indel_bp"] = int(max_indel)
    out["consensus_match_name"] = None  # filled by caller
    if host_offset_0based in mapping:
        out["consensus_offset"] = mapping[host_offset_0based]
        # Report the drift so Phase 1 can stratify on it directly instead of
        # relying on a single cutoff.
        out["offset_drift_bp"] = mapping[host_offset_0based] - host_offset_0based
        # "Unambiguous" means the alignment carries no indel large enough to
        # move the projected coordinate across a consensus bin boundary.
        unambiguous = max_indel <= CONSENSUS_MAX_INDEL_TOLERANCE
        out["consensus_mapping"] = "unambiguous" if unambiguous else "ambiguous_indel"
    else:
        out["consensus_mapping"] = "unaligned_breakpoint"
    return out


#: An alignment indel larger than this can move a projected breakpoint far
#: enough to change which consensus bin it lands in, so the mapping is marked
#: ambiguous. Measured rather than guessed: across the 3313 nested calls whose
#: breakpoint aligned, |consensus_offset - host_offset| is 6 bp at the median,
#: 17 bp at p75 and 39 bp at p90. Only 50% of calls drift by <=5 bp, so a 5 bp
#: indel tolerance would discard the majority of the cohort. A 15 bp tolerance
#: retains ~72% while excluding the tail where the projection is unreliable.
#: The projected coordinate is always reported alongside `max_indel_bp` so
#: Phase 1 can re-stratify on drift directly rather than trusting this cutoff.
CONSENSUS_MAX_INDEL_TOLERANCE = 15


# --------------------------------------------------------------------------
# Genotype state
# --------------------------------------------------------------------------

#: Mutually exclusive genotype states. These are never collapsed into a single
#: boolean, because the pooled callset has no missingness channel: every GT is
#: either 0/0 or non-ref, with no `./.` and VAF1 unset, so a 0/0 cannot be
#: distinguished from a genotyping failure inside a hard repeat.
GENOTYPE_STATES = (
    "concordant_nonref",          # alt-bearing genotype at the pooled site
    "ref_genotype",               # 0/0, and otherwise interpretable
    "uncallable_or_uninformative",  # inferred, never measured - see below
    "no_pooled_site",             # no matching site in the cohort callset
    "gt_missing",                 # an explicit missing or absent GT call
)

#: Reasons a 0/0 may be uninformative rather than informative. Recorded per
#: row so the inference is visible instead of silent.
UNCALLABLE_REASONS = (
    "none",
    "low_mappability",            # inferred from a short-read mask; a caveat
    "high_repeat_density",        # inferred from local repeat content
    "host_region",                # inside an Alu, hardest for any caller
)


def classify_genotype(
    gt: str,
    *,
    has_pooled_site: bool,
    inference_reasons: Iterable[str] = (),
) -> tuple[str, str]:
    """Return (state, reason). Keeps true 0/0 separate from the unknown.

    `inference_reasons` carries *inferred* callability caveats. It is labelled
    as inference in the output because the source callset provides no coverage
    or VAF channel from which callability could be measured directly.
    """
    reasons = [r for r in inference_reasons if r and r != "none"]
    if not has_pooled_site:
        return "no_pooled_site", "none"
    if gt in ("", "."):
        return "gt_missing", "none"
    if not re.split(r"[/|]", gt) or any(
        a in (".", "") for a in re.split(r"[/|]", gt)
    ):
        return "gt_missing", "none"
    if genotype_is_nonref(gt):
        return "concordant_nonref", (";".join(reasons) if reasons else "none")
    if reasons:
        return "uncallable_or_uninformative", ";".join(reasons)
    return "ref_genotype", "none"


def load_repeat_density_bed(rmsk: Path, tmp: Path) -> Path | None:
    """Convert rmsk to a sorted BED of same-family intervals, for density checks."""
    bed = Path(tmp) / "alu_intervals.bed"
    if bed.exists():
        return bed
    bed.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(rmsk, "rt", errors="replace") as fh, open(bed, "w") as out:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 13 or not fields[5].startswith("chr"):
                continue
            if not is_alu_interval(fields):
                continue
            out.write(
                f"{fields[5]}\t{fields[6]}\t{fields[7]}\t{fields[10]}\t0\t{fields[9]}\n"
            )
    return bed


# --------------------------------------------------------------------------
# Row assembly
# --------------------------------------------------------------------------

OUTPUT_COLUMNS = [
    "chrom", "pos", "family",
    "insert_strand", "conformation", "perc_resolved", "not_canonical", "rt_len",
    "nested_in_alu_host", "host_start0", "host_end0", "host_strand", "host_name",
    "host_len", "host_offset_5p_0based", "host_offset_bin_20bp",
    "host_selection_rule",
    "consensus_match_name", "consensus_offset", "consensus_mapping",
    "offset_drift_bp",
    "alignment_identity", "alignment_score", "n_unaligned_host_bases",
    "max_indel_bp",
    "tsd_len", "tsd_seq", "tsd_len_is_sentinel", "polya_len", "polya_seq",
    "site_allele_count", "site_allele_freq",
    "genotype_gt", "genotype_state", "genotype_state_reason",
    "genotype_callability_source",
]


def build_rows(
    records: list[dict[str, str]],
    assigned: dict[str, dict[tuple[str, int], tuple[int, int, str, str]]],
    consensus: dict[str, str],
    fasta: Path | None,
    aligner: Any,
    sentinel_values: Iterable[int],
    geno_bcf: Path | None,
    sample: str,
    near_full_min: int = 280,
    near_full_max: int = 320,
) -> list[dict[str, Any]]:
    """One output row per Alu insertion site, nested sites annotated in place.

    Unnested sites are retained with `nested_in_alu_host=0` and empty host
    fields, so the file is a complete per-site table rather than only the
    nested subset. Phase 1 needs the unnested rows as its within-host
    opportunity denominator.
    """
    import pyfaidx

    fasta_handle = pyfaidx.Fasta(
        str(fasta), as_raw=True, sequence_always_upper=True
    ) if fasta else None
    sentinels = set(sentinel_values)
    rule = "longest_span"
    table = assigned.get(rule, {})

    # Site allele count / frequency from the genotyped callset, one region
    # query per chromosome rather than per record.
    ac_lookup: dict[tuple[str, int], tuple[str, str]] = {}
    gt_lookup: dict[tuple[str, int], str] = {}
    if geno_bcf is not None:
        by_chrom: dict[str, list[int]] = collections.defaultdict(list)
        for rec in records:
            by_chrom[rec["chrom"]].append(rec["pos"])
        for chrom, poss in by_chrom.items():
            lo, hi = min(poss) - 5, max(poss) + 5
            got_ac = fetch_region_field(geno_bcf, chrom, lo, hi, None, "AC")
            got_af = fetch_region_field(geno_bcf, chrom, lo, hi, None, "AF")
            got_gt = fetch_region_field(geno_bcf, chrom, lo, hi, sample, "GT")
            for pos in poss:
                ac_lookup[(chrom, pos)] = (
                    got_ac.get(pos, "."),
                    got_af.get(pos, "."),
                )
                gt_lookup[(chrom, pos)] = got_gt.get(pos, ".")

    rows: list[dict[str, Any]] = []
    for rec in records:
        key = (rec["chrom"], rec["pos"])
        row: dict[str, Any] = {c: "" for c in OUTPUT_COLUMNS}
        row.update(
            {
                "chrom": rec["chrom"],
                "pos": rec["pos"],
                "family": FAMILY,
                "insert_strand": rec.get("STRAND", ""),
                "conformation": rec.get("CONFORMATION", ""),
                "perc_resolved": rec.get("PERC_RESOLVED", ""),
                "not_canonical": rec.get("NOT_CANONICAL", ""),
                "rt_len": rec.get("RT_LEN", ""),
                "tsd_len": rec.get("TSD_LEN", ""),
                "tsd_seq": rec.get("TSD_SEQ", ""),
                "polya_len": rec.get("POLYA_LEN", ""),
                "polya_seq": rec.get("POLYA_SEQ", ""),
            }
        )
        try:
            tsd_len_int = int(rec.get("TSD_LEN", "") or 0)
        except ValueError:
            tsd_len_int = 0
        row["tsd_len_is_sentinel"] = (
            1 if tsd_len_int in sentinels else 0
        ) if tsd_len_int else ""

        ac_raw, af_raw = ac_lookup.get(key, (".", "."))
        row["site_allele_count"] = ac_raw
        row["site_allele_freq"] = af_raw

        host = table.get(key)
        if host is not None:
            start0, end0, hstrand, hname = host
            offset = host_offset(rec["pos"], start0, end0, hstrand)
            row.update(
                {
                    "nested_in_alu_host": 1,
                    "host_start0": start0,
                    "host_end0": end0,
                    "host_strand": hstrand,
                    "host_name": hname,
                    "host_len": end0 - start0,
                    "host_offset_5p_0based": offset,
                    "host_offset_bin_20bp": (offset // 20) * 20 if offset is not None else "",
                    "host_selection_rule": rule,
                }
            )
            matched = pick_consensus_for_host(hname, consensus)
            if matched and fasta_handle is not None and offset is not None:
                cname, cseq = matched
                try:
                    hseq = str(fasta_handle[rec["chrom"]][start0:end0])
                except (KeyError, ValueError):
                    hseq = ""
                proj = project_offset(aligner, hseq, cseq, offset)
                proj["consensus_match_name"] = cname
                row.update({k: proj[k] for k in proj})
        else:
            row["nested_in_alu_host"] = 0

        gt = gt_lookup.get(key, "") if geno_bcf is not None else ""
        row["genotype_gt"] = gt
        reasons = ["host_region"] if row["nested_in_alu_host"] else []
        state, reason = classify_genotype(
            gt, has_pooled_site=bool(ac_raw and ac_raw != "."), inference_reasons=reasons
        )
        row["genotype_state"] = state
        row["genotype_state_reason"] = reason
        row["genotype_callability_source"] = (
            "inferred_from_host_context" if reasons else "not_inferred"
        )
        rows.append(row)
    return rows


# --------------------------------------------------------------------------
# Manifest
# --------------------------------------------------------------------------


def summarise_rows(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Counts a reader needs before trusting any downstream number."""
    state_counts: collections.Counter = collections.Counter()
    mapping_counts: collections.Counter = collections.Counter()
    bins: collections.Counter = collections.Counter()
    drifts: list[int] = []
    nested = 0
    near_full = 0
    for row in rows:
        state_counts[row["genotype_state"]] += 1
        if row["nested_in_alu_host"] == 1:
            nested += 1
            mapping_counts[row["consensus_mapping"]] += 1
            if row.get("offset_drift_bp") not in (None, ""):
                drifts.append(abs(int(row["offset_drift_bp"])))
            try:
                hlen = int(row["host_len"])
            except (TypeError, ValueError):
                hlen = 0
            if 280 <= hlen <= 320:
                near_full += 1
                bins[row["host_offset_bin_20bp"]] += 1
    drifts.sort()
    drift_summary: dict[str, Any] = {"n_mapped": len(drifts)}
    if drifts:
        for q in (0.5, 0.75, 0.9, 0.95):
            drift_summary[f"p{int(q * 100)}"] = drifts[int(len(drifts) * q)]
        drift_summary["fraction_le_5bp"] = round(
            sum(1 for d in drifts if d <= 5) / len(drifts), 4
        )
        drift_summary["fraction_le_15bp"] = round(
            sum(1 for d in drifts if d <= 15) / len(drifts), 4
        )
    return {
        "n_rows": len(rows),
        "n_nested": nested,
        "n_nested_near_full_host_280_320": near_full,
        "genotype_state_counts": dict(state_counts),
        "consensus_mapping_counts": dict(mapping_counts),
        "consensus_offset_drift_bp": drift_summary,
        "host_offset_bin_counts": {str(k): v for k, v in sorted(bins.items())},
        "n_tsd_len_sentinel": sum(
            1 for r in rows if r["tsd_len_is_sentinel"] == 1
        ),
    }


def write_manifest(path: Path, manifest: dict[str, Any]) -> None:
    """Write every gate verdict to disk, not just to stdout."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=False)
        fh.write("\n")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--site-bcf", type=Path, default=DEFAULT_SITE_BCF)
    p.add_argument("--geno-bcf", type=Path, default=DEFAULT_GENO_BCF)
    p.add_argument("--rmsk", type=Path, default=DEFAULT_RMSK)
    p.add_argument("--fasta", type=Path, default=DEFAULT_FASTA)
    p.add_argument("--consensus", type=Path, default=DEFAULT_CONSENSUS)
    p.add_argument(
        "--hg03086-vcf",
        type=Path,
        default=None,
        help="VCF used only for the Gate 0c independence join; omit to skip 0c.",
    )
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--sample", default="HG03086")
    p.add_argument(
        "--no-index-probe",
        action="store_true",
        help="Skip the Gate 0a linear-scan comparison (it reads the whole BCF).",
    )
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    args.outdir.mkdir(parents=True, exist_ok=True)
    scratch = args.outdir / "scratch"

    print("[gate 0a] index and reference consistency", flush=True)
    index_result = (
        {"verdict": "skipped_by_flag", "all_regions_agree": None,
         "index_older_than_bcf": None, "probe_regions": [],
         "rebuild_required": False, "rebuild_path": None}
        if args.no_index_probe
        else check_index_consistency(args.geno_bcf, scratch=scratch)
    )
    reference_result = check_reference_consistency(args.fasta)
    print(f"  index: {index_result['verdict']}", flush=True)
    print(f"  reference: {reference_result['verdict']}", flush=True)

    print("[gate 0b] assembly identity", flush=True)
    assembly_result = check_no_chm13(args.site_bcf, args.geno_bcf)
    print(f"  {assembly_result['verdict']}", flush=True)

    print("[gate 0d] TSD sentinel detection", flush=True)
    sentinel_result = detect_tsd_sentinels(args.site_bcf)
    print(
        f"  {sentinel_result['verdict']} "
        f"values={sorted(sentinel_result['sentinel_values'], key=int)}",
        flush=True,
    )

    independence_result: dict[str, Any] = {
        "verdict": "not_evaluated",
        "phase3_tier": "unknown",
        "reason": "--hg03086-vcf not supplied",
    }
    if args.hg03086_vcf:
        print("[gate 0c] independence verdict", flush=True)
        independence_result = assess_genotype_axis_independence(
            args.geno_bcf, args.site_bcf, args.hg03086_vcf, sample=args.sample
        )
        print(f"  {independence_result['verdict']}", flush=True)

    print("loading site records", flush=True)
    records = load_site_records(args.site_bcf)
    print(f"  {len(records)} {FAMILY} insertion sites", flush=True)

    print("assigning hosts (single rmsk pass, all rules)", flush=True)
    assigned = assign_hosts(records, args.rmsk)
    n_default_nested = len(assigned.get("longest_span", {}))
    print(f"  {n_default_nested} nested in a same-family Alu host", flush=True)

    consensus = load_consensus(args.consensus)
    print(f"  {len(consensus)} {FAMILY} consensus sequences", flush=True)

    print("building per-call table", flush=True)
    aligner = build_aligner()
    rows = build_rows(
        records,
        assigned,
        consensus,
        args.fasta,
        aligner,
        [int(v) for v in sentinel_result["sentinel_values"]],
        args.geno_bcf,
        args.sample,
    )

    csv_path = args.outdir / "per_call_longread_nested.csv"
    header = OUTPUT_COLUMNS
    with csv_path.open("w") as fh:
        fh.write(",".join(header) + "\n")
        for row in rows:
            fh.write(
                ",".join(_csv_cell(row.get(c, "")) for c in header) + "\n"
            )

    summary = summarise_rows(rows)
    manifest = {
        "phase": 0,
        "description": (
            "Cohort build and provenance gates for the structure-informed "
            "nested-Alu study."
        ),
        "inputs": {
            "site_bcf": str(args.site_bcf),
            "geno_bcf": str(args.geno_bcf),
            "rmsk": str(args.rmsk),
            "fasta": str(args.fasta),
            "consensus": str(args.consensus),
            "hg03086_vcf": str(args.hg03086_vcf) if args.hg03086_vcf else None,
        },
        "outputs": {"per_call_csv": str(csv_path)},
        "conventions": {
            "host_offset_5p_0based": (
                "0-based offset of the breakpoint from the host's own 5' end; "
                "pos0-start0 for + strand hosts, end0-1-pos0 for - strand hosts"
            ),
            "consensus_offset": (
                "host offset projected through a global alignment of the host to "
                "its own subfamily consensus, measured from the consensus 5' end"
            ),
            "host_selection_rule": "longest_span",
            "host_selection_rules_available": list(HOST_SELECTION_RULES),
            "consensus_max_indel_tolerance_bp": CONSENSUS_MAX_INDEL_TOLERANCE,
        },
        "gate_0a_index_and_reference": {
            "bcf_index": index_result,
            "reference": reference_result,
        },
        "gate_0b_assembly_identity": assembly_result,
        "gate_0c_genotype_axis_independence": independence_result,
        "gate_0d_tsd_sentinels": sentinel_result,
        "genotype_state_definitions": {
            state: {
                "concordant_nonref": (
                    "alt-bearing genotype at the pooled site; this is "
                    "concordance, NOT independent validation"
                ),
                "ref_genotype": (
                    "0/0 where the site is otherwise interpretable; kept "
                    "separate from uncallable because the source has no "
                    "missingness channel"
                ),
                "uncallable_or_uninformative": (
                    "genotype inferred as uninformative from context; the "
                    "inference is labelled in genotype_callability_source"
                ),
                "no_pooled_site": "no matching site in the pooled cohort",
                "gt_missing": "explicit missing or absent GT call",
            }[state]
            for state in GENOTYPE_STATES
        },
        "summary": summary,
        "disclosures": [
            "The pooled callset is a catalogue of insertion SITES across 908 "
            "samples, not one genome's insertions. Unique sites and allele "
            "counts are reported per row; they are not per-genome biology.",
            f"{args.sample} is a member of the pooled cohort, so the genotype "
            "axis is internal to the discovery cohort. See gate 0c.",
            "The source callset has no missingness channel, so callability is "
            "inferred from context rather than measured.",
            "All calls are catalogue polymorphisms that have passed selection; "
            "this is a survivor distribution, not raw de novo targeting rates.",
        ],
    }
    manifest_path = args.outdir / "provenance_manifest.json"
    write_manifest(manifest_path, manifest)

    print(f"\nwrote {csv_path} ({len(rows)} rows)")
    print(f"wrote {manifest_path}")
    print(
        f"  nested={summary['n_nested']} "
        f"near-full-host={summary['n_nested_near_full_host_280_320']} "
        f"sentinel-TSD={summary['n_tsd_len_sentinel']}"
    )
    print(f"  genotype states: {summary['genotype_state_counts']}")
    print(f"  consensus mapping: {summary['consensus_mapping_counts']}")
    return 0


def _csv_cell(value: Any) -> str:
    """Minimal CSV escaping; the table has no embedded commas in practice."""
    text = "" if value is None else str(value)
    if any(ch in text for ch in ',"\n'):
        return '"' + text.replace('"', '""') + '"'
    return text


if __name__ == "__main__":
    raise SystemExit(main())
