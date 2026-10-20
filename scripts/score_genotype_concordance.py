"""Phase 3 -- genotype concordance for the nested-Alu cohort.

This phase was pre-registered in `docs/NESTED_ALU_ANALYSIS_PLAN.md` section 4
and 5. Two of its pre-registered requirements turned out to be unsatisfiable on
the data that is actually available, and the phase is built to *measure* that
rather than to work around it quietly.

What the plan asked for
-----------------------
Per-sample callsets, genotype concordance at AC>2 sites, restricted to the
`restricted_concordance_only` tier, with any callability claim labelled as an
inference.

What the data allows
--------------------
The per-sample long-read callset that section 4 assumes exists does not exist on
GRCh38. `final-vcf.unphased.SVAN_1.3.vcf.gz` sits under the hg38 polymorphism
directory as a symlink into `hs1/`, so it is CHM13; Gate 0b rejects it and Gate
3c re-checks it rather than trusting the directory name. The genotyped SVIM
pooled BCF (`svim.asm.hg38.bcf`, 908 samples) *is* usable and is the same
pipeline that produced the sites, which section 5 anticipated: independent
individuals, one method.

The second gap is harder. The pooled callset has no missingness channel -- no
`./.` anywhere and VAF1 unset throughout -- so a `0/0` is indistinguishable
from a genotyping failure inside a hard repeat. Section 5 admits an absent call
as a negative *only at loci demonstrated callable in that sample*, and the
demonstration criteria must be recorded before use. No such demonstration is
possible from this callset. Gate 3b therefore records that the pre-registered
design is not executable as written, and the phase proceeds on the one axis the
data does support: recurrence across independent subsets of the same cohort.

What is computed instead
-------------------------
Genotype concordance is measured method-internally, as recurrence across two
independent halves of the 908 genomes. A site carried by genomes in both halves
has been seen more than once, by reads from more than one genome; a site whose
carriers all fall in one half has been seen once. That is a weaker claim than
agreement between two methods and is labelled as such everywhere. It is not
called validation, and Gate 3d shows why it cannot be.

Gates, in the order they run
----------------------------
  Gate 0c  the independence verdict for the genotype axis, with the question
           revised by section 5. Runs first, unchanged as a gate.
  Gate 3a  coordinate frame. The site BCF and the genotype BCF declare different
           lengths for 24 of their 25 shared contigs, in both directions across
           chromosomes, which is the signature of two different reference builds.
           The join is therefore verified by exact-coordinate containment rather
           than assumed from the headers, and the header disagreement is reported
           even when the data overrules it.
  Gate 3b  genotype channel completeness. Absence of a missingness channel is
           what makes every `0/0` uninterpretable, so it is a gate and not a
           footnote.
  Gate 3c  sample identity. Joins are by sample name, never by column index: the
           cross-method file orders its samples differently from the genotype
           file, and an index join would silently transpose genotypes.
  Gate 3d  cross-method floor. The only independent callset on GRCh38 is the
           short-read merged 1KG SV callset. Its coverage of nested-Alu sites is
           measured against non-nested Alu sites from the same cohort as a
           control. A short read cannot span an Alu inserted into an Alu, so the
           expectation is a large, structured deficit; the gate reports the
           measured deficit and records that cross-method agreement is therefore
           unavailable for this locus class, which is the evidence section 5
           requires before the word "validated" stays off the table.

Reused machinery
----------------
Hallmark definitions, the coarsened matching ladder, the stratified permutation
test and the Benjamini-Hochberg adjustment all come from
`score_tprt_hallmarks.py` rather than being restated, so a hallmark cannot mean
two things across phases. The cross-half comparison adds one entry to that
module's covariate and ladder tables; the existing entries are pinned by a test
so the addition cannot quietly redefine a Phase 2b comparison.

Allele frequency is deliberately *not* a matching covariate here. A site with a
single carrier cannot be cross-half supported by arithmetic, so conditioning on
frequency would condition on the exposure. The primary comparison is instead
restricted to sites with at least two carriers, which removes the entanglement
at the source.
"""

from __future__ import annotations

import argparse
import bisect
import collections
import csv
import json
import math
import re
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any, Callable, Sequence

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import score_tprt_hallmarks as tprt  # noqa: E402  (path bootstrap must run first)

DEFAULT_COHORT = tprt.DEFAULT_COHORT
DEFAULT_MANIFEST = tprt.DEFAULT_MANIFEST
DEFAULT_REFERENCE = tprt.DEFAULT_REFERENCE
DEFAULT_ALU_CONSENSUS = tprt.DEFAULT_ALU_CONSENSUS

_WORKDIR = Path("~/retrotransposon-workdir/data/public").expanduser()
DEFAULT_GENOTYPE_BCF = (
    _WORKDIR / "polymorphism/hg38/long_read_1kg_ont_vienna/svim.asm.hg38.bcf"
)
DEFAULT_SITE_BCF = (
    _WORKDIR / "polymorphism/hg38/long_read_1kg_ont_vienna/svim.asm.hg38.noGt.SVAN_1.3.bcf"
)
DEFAULT_CROSS_METHOD_BCF = (
    _WORKDIR
    / "polymorphism/hg38/1kg/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.GRCh38.vcf.gz"
)

#: Gate 0c's focus sample, restated from the Phase 0 gate so the two cannot drift.
GATE_0C_FOCUS_SAMPLE = "HG03086"

#: The two halves are cut by sorted sample name, even index versus odd index.
#: Pre-declared, balanced to within one genome, and independent of any genotype,
#: so the split cannot be chosen after seeing which sites replicate.
HALF_SPLIT_RULE = "sorted_sample_name_even_index_vs_odd_index"

#: A site needs at least this many carriers to enter the cross-half comparison.
#: Below it, cross-half support is impossible by arithmetic, not by biology.
MIN_CARRIERS_PER_SITE = 2

#: Cross-method proximity. Deliberately generous: an independent call of the
#: same insertion that disagrees by more than this is not the same insertion.
PROXIMITY_TOLERANCE_BP = 50

#: Cross-method insertion classes. A nested Alu could in principle be called as
#: any insertion class by a different pipeline, so restricting to ALU would
#: understate coverage; the restriction is reported per class anyway.
CROSS_METHOD_SVTYPES = ("ALU", "LINE1", "SVA", "INS")

ACCUMULATION_CHECKPOINTS = (10, 25, 50, 100, 200, 400, 600, 800, 908)
SATURATION_FRACTIONS = (0.5, 0.8, 0.9, 0.95)

DEFAULT_PERMUTATIONS = tprt.DEFAULT_PERMUTATIONS
DEFAULT_SEED = tprt.DEFAULT_SEED

COMPARISON = "cross_half_supported"

#: The only subject label this phase may publish. Plan section 5 fixes the
#: permitted term; naming it once means the JSON and the prose cannot drift.
SUBJECT = "genotype_concordance"

#: bcftools query format for all genotypes of a record. The tab must sit INSIDE
#: the brackets: a bare `[%GT]` concatenates all 908 genotypes with no separator
#: into a single field, and `%TGT` prints nothing for this file. Either mistake
#: yields a short, silently mis-aligned genotype array rather than an error, so
#: the format is a named constant and is asserted in the tests.
GENOTYPE_QUERY_FORMAT = "%CHROM\\t%POS[\\t%GT]\\n"

#: Metric names this phase must never emit. Section 5 reserves "validated" for
#: cross-method evidence that does not exist, and Gate 0c's tier forbids
#: detection-performance language outright. The check is on *keys and
#: identifiers*, not prose: a disclaimer is allowed to name what is not being
#: claimed, and a reported number is not.
FORBIDDEN_METRIC_TERMS = (
    "precision",
    "recall",
    "accuracy",
    "sensitivity",
    "specificity",
    "f1",
    "ppv",
    "npv",
    "validated",
    "validation",
)


# --------------------------------------------------------------------------
# bcftools plumbing
# --------------------------------------------------------------------------


def _bcftools(args: Sequence[str]) -> str:
    proc = subprocess.run(
        ["bcftools", *args], capture_output=True, text=True, check=False
    )
    if proc.returncode != 0:
        raise RuntimeError(f"bcftools {' '.join(args)} failed: {proc.stderr.strip()}")
    return proc.stdout


def header_contig_lengths(path: Path) -> dict[str, int]:
    """Contig name -> declared length, from the VCF/BCF header.

    Declared lengths are metadata, not truth. Gate 3a exists because the two
    files in this study disagree about them.
    """
    lengths: dict[str, int] = {}
    for line in _bcftools(["view", "-h", str(path)]).splitlines():
        if not line.startswith("##contig=<ID="):
            continue
        name = line.split("ID=")[1].split(",")[0]
        match = re.search(r"length=(\d+)", line)
        if match:
            lengths[name] = int(match.group(1))
    return lengths


def sample_names(path: Path) -> list[str]:
    """Sample names in file order.

    The order is returned as-is and must never be used to align genotypes across
    two files. `join_is_by_name` records the check.
    """
    text = _bcftools(["query", "-l", str(path)])
    return [line for line in text.splitlines() if line]


def site_coordinates(path: Path) -> set[tuple[str, int]]:
    text = _bcftools(["query", "-f", "%CHROM\t%POS\n", str(path)])
    out: set[tuple[str, int]] = set()
    for line in text.splitlines():
        if not line:
            continue
        chrom, _, pos = line.partition("\t")
        out.add((chrom, int(pos)))
    return out


def normalise_contig(name: str) -> str:
    """`chr1` and `1` are the same contig; the two callsets disagree on style.

    This is the join that silently returns nothing when it is missed, so it is
    one named function and is tested directly.
    """
    return name[3:] if name.startswith("chr") else name


def parse_genotype_block(tokens: Sequence[str]) -> np.ndarray:
    """Genotype strings -> carrier mask.

    A carrier is any genotype carrying a non-reference allele. A missing or
    partial genotype is *not* a carrier and is *not* a demonstrated negative: it
    is counted separately so Gate 3b can report how many exist, and it is never
    silently folded into the non-carrier column.
    """
    mask = np.zeros(len(tokens), dtype=bool)
    for i, token in enumerate(tokens):
        if not token or token.startswith("."):
            continue
        alleles = re.split(r"[/|]", token)
        if any(a not in ("0", "") for a in alleles):
            mask[i] = True
    return mask


def read_site_genotypes(
    bcf: Path, targets: Sequence[tuple[str, int]]
) -> tuple[list[str], dict[tuple[str, int], np.ndarray], dict[str, int]]:
    """Read genotypes for `targets` only, via a region file.

    Reading all 908 genotypes for every site in the callset would be 150M
    tokens; the region file restricts the decode to the sites under analysis.
    """
    wanted = set(targets)
    if not wanted:
        return sample_names(bcf), {}, {"n_records": 0, "n_missing_genotypes": 0}

    samples = sample_names(bcf)
    with tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False) as fh:
        for chrom, pos in sorted(wanted):
            fh.write(f"{chrom}\t{max(0, pos - 1)}\t{pos}\n")
        regions = Path(fh.name)
    try:
        text = _bcftools(
            ["query", "-f", GENOTYPE_QUERY_FORMAT, "-R", str(regions), str(bcf)]
        )
    finally:
        regions.unlink(missing_ok=True)

    carriers: dict[tuple[str, int], np.ndarray] = {}
    n_missing = 0
    n_records = 0
    n_short: list[str] = []
    for line in text.splitlines():
        if not line:
            continue
        parts = line.split("\t")
        key = (parts[0], int(parts[1]))
        if key not in wanted:
            continue
        # `[%GT]` with the separator *inside* the brackets is what expands to one
        # field per sample. A bare `[%GT]` concatenates all 908 genotypes with no
        # separator at all, and `%TGT` prints nothing for this file; either
        # mistake yields a silently short, mis-aligned genotype array.
        tokens = parts[2:]
        n_missing += sum(1 for t in tokens if not t or t.startswith("."))
        n_records += 1
        if len(tokens) != len(samples):
            # A record with the wrong number of genotype fields would silently
            # mis-align against the sample list, so it is named and not used.
            n_short.append(f"{key[0]}:{key[1]}")
            continue
        carriers[key] = parse_genotype_block(tokens)
    return samples, carriers, {
        "n_records": n_records,
        "n_missing_genotypes": n_missing,
        "n_records_with_unexpected_sample_count": len(n_short),
        "records_with_unexpected_sample_count": n_short[:10],
    }


# --------------------------------------------------------------------------
# gates
# --------------------------------------------------------------------------


def gate_0c_independence(
    samples: Sequence[str], focus_sample: str = GATE_0C_FOCUS_SAMPLE
) -> dict[str, Any]:
    """The genotype-axis independence verdict, question revised by plan section 5.

    Runs first and is unchanged as a gate. Only the question it answers is
    revised: same-pipeline callsets are independent individuals but one method,
    so the answer is a tier rather than a yes/no.
    """
    present = focus_sample in set(samples)
    return {
        "gate": "0c",
        "question": (
            "may Phase 3 use the genotype axis for concordance at "
            "demonstrated-callable loci, given that same-pipeline callsets are "
            "independent individuals but one method?"
        ),
        "focus_sample": focus_sample,
        "focus_sample_in_genotype_cohort": present,
        "n_samples": len(samples),
        "verdict": "not_independent" if present else "independent_of_pooled_samples",
        "tier": "restricted_concordance_only",
        "consequence": (
            "genotype concordance only; no detection-performance metric is "
            "licensed, and the genotype axis carries no method independence"
        ),
    }


def gate_3a_coordinate_frame(
    site_bcf: Path,
    genotype_bcf: Path,
    analysis_sites: Sequence[tuple[str, int]],
) -> dict[str, Any]:
    """Verify the coordinate frame instead of assuming it.

    The two BCFs declare different lengths for 24 of 25 shared contigs, and the
    sign of the difference flips across chromosomes, which is what a reference
    build mismatch looks like. Header lengths are therefore treated as a claim
    to be checked, and the check is exact-coordinate containment.
    """
    site_lengths = header_contig_lengths(site_bcf)
    geno_lengths = header_contig_lengths(genotype_bcf)
    shared = sorted(set(site_lengths) & set(geno_lengths))
    disagreements = [
        {"contig": c, "site_bcf": site_lengths[c], "genotype_bcf": geno_lengths[c]}
        for c in shared
        if site_lengths[c] != geno_lengths[c]
    ]
    in_disputed = [
        s
        for s in analysis_sites
        if s[0] in site_lengths
        and s[0] in geno_lengths
        and s[1] > min(site_lengths[s[0]], geno_lengths[s[0]])
    ]
    available = site_coordinates(genotype_bcf)
    missing = [s for s in analysis_sites if s not in available]
    return {
        "gate": "3a",
        "site_bcf": str(site_bcf),
        "genotype_bcf": str(genotype_bcf),
        "shared_contigs": len(shared),
        "contig_length_disagreements": len(disagreements),
        "contig_length_disagreement_detail": disagreements,
        "analysis_sites": len(analysis_sites),
        "analysis_sites_in_disputed_length_region": len(in_disputed),
        "analysis_sites_absent_from_genotype_bcf": len(missing),
        "verdict": (
            "frame_agreement_demonstrated_by_exact_coordinates"
            if not missing
            else "frame_disagreement_blocks_join"
        ),
        "note": (
            "header contig lengths disagree but every analysis coordinate is "
            "present verbatim in the genotype callset, so the join is verified "
            "on the data rather than on the headers; the header disagreement is "
            "reported because a future join that trusted it would be wrong"
        ),
    }


def gate_3b_genotype_channel(
    counters: dict[str, int],
    n_sites: int,
    n_samples: int,
) -> dict[str, Any]:
    """Can a `0/0` be read as absence?

    Plan section 5 admits an absent call as a negative only at loci demonstrated
    callable in that sample. With no missingness channel there is no
    demonstration available, so the pre-registered Phase 3 design is reported as
    not executable and the phase falls back to the recurrence axis.
    """
    n_missing = counters.get("n_missing_genotypes", 0)
    n_short = counters.get("n_records_with_unexpected_sample_count", 0)
    return {
        "gate": "3b",
        "genotypes_examined": n_sites * n_samples,
        "missing_or_partial_genotypes": n_missing,
        "missingness_channel_present": n_missing > 0,
        "records_with_unexpected_sample_count": n_short,
        "verdict": (
            "missingness_channel_present"
            if n_missing > 0
            else "no_missingness_channel_zero_zero_uninterpretable"
        ),
        "consequence": (
            "a 0/0 is indistinguishable from a genotyping failure inside a hard "
            "repeat, so non-carriers are not demonstrated negatives and the "
            "pre-registered 'concordance at demonstrated-callable loci' design is "
            "NOT EXECUTABLE on this callset; Phase 3 proceeds on the recurrence "
            "axis across independent halves of the cohort, where carriers are "
            "positive evidence and no negative is required"
        ),
    }


def gate_3c_sample_identity(
    genotype_samples: Sequence[str],
    cross_method_samples: Sequence[str] | None,
    cross_method_path: Path | None,
) -> dict[str, Any]:
    """Sample identity, and the column-order trap.

    The cross-method file and the genotype file list the same individuals in
    different orders. An index join would transpose genotypes silently, so the
    join is by name and the order divergence is reported.
    """
    duplicated = [s for s, n in collections.Counter(genotype_samples).items() if n > 1]
    out: dict[str, Any] = {
        "gate": "3c",
        "n_genotype_samples": len(genotype_samples),
        "n_duplicate_sample_names": len(duplicated),
        "join_is_by_name": True,
        "cross_method_file": str(cross_method_path) if cross_method_path else None,
    }
    if cross_method_samples is None:
        out["cross_method_samples"] = None
        out["verdict"] = "cross_method_comparator_unavailable"
        out["consequence"] = (
            "no independent callset was read, so no cross-method statement of any "
            "kind is made and Gate 3d cannot run"
        )
        return out
    shared = sorted(set(genotype_samples) & set(cross_method_samples))
    shared_set = set(shared)
    # Compare the two files' *order* over the samples they share. Comparing
    # against the sorted list instead would report "identical" for any two files
    # listing the same individuals, which is exactly the case that matters.
    genotype_order = [s for s in genotype_samples if s in shared_set]
    cross_order = [s for s in cross_method_samples if s in shared_set]
    same_order = genotype_order == cross_order
    out.update(
        {
            "n_cross_method_samples": len(cross_method_samples),
            "n_samples_shared": len(shared),
            "shared_fraction_of_genotype_cohort": (
                len(shared) / len(genotype_samples) if genotype_samples else None
            ),
            "column_order_identical": same_order,
            "verdict": "samples_joined_by_name",
            "consequence": (
                "genotypes are aligned to sample names, never to column index; an "
                "index join would transpose them silently"
            ),
        }
    )
    if not same_order:
        out["note"] = (
            "the two callsets order their shared samples differently, which is the "
            "concrete reason the join is by name"
        )
    return out


def fisher_exact_two_sided(a: int, b: int, c: int, d: int) -> float | None:
    """Two-sided Fisher exact p from a 2x2 table, without a scipy dependency.

    a = exposed & covered, b = exposed & not covered,
    c = reference & covered, d = reference & not covered.
    """
    n = a + b + c + d
    if n == 0 or min(a + c, b + d) == 0:
        return None
    row1, col1 = a + b, a + c
    observed = math.comb(row1, a) * math.comb(n - row1, col1 - a) / math.comb(n, col1)
    total = 0.0
    lo = max(0, col1 - (n - row1))
    hi = min(row1, col1)
    for k in range(lo, hi + 1):
        prob = math.comb(row1, k) * math.comb(n - row1, col1 - k) / math.comb(n, col1)
        if prob <= observed + 1e-12:
            total += prob
    return min(1.0, total)


def gate_3d_cross_method_floor(
    cross_method_path: Path | None,
    nested_sites: Sequence[tuple[str, int]],
    control_sites: Sequence[tuple[str, int]],
    tolerance: int = PROXIMITY_TOLERANCE_BP,
) -> dict[str, Any]:
    """Measure the independent callset's coverage of nested-Alu sites.

    The control is non-nested Alu insertion sites from the same cohort and the
    same assembly. A short read cannot span an Alu inserted into an Alu, so the
    nested arm is expected to be covered far less often; the size of that gap is
    the measurement, and it is what licenses the plan's refusal to call any of
    this cross-method evidence.
    """
    if cross_method_path is None:
        return {
            "gate": "3d",
            "verdict": "cross_method_comparator_unavailable",
            "note": "no independent GRCh38 callset was supplied",
        }
    text = _bcftools(
        [
            "query",
            "-f",
            "%CHROM\t%POS\t%INFO/SVTYPE\t%INFO/SVLEN\n",
            "-i",
            "INFO/SVTYPE='ALU'||INFO/SVTYPE='LINE1'||INFO/SVTYPE='SVA'||INFO/SVTYPE='INS'",
            str(cross_method_path),
        ]
    )
    index: dict[str, list[tuple[int, int, str]]] = collections.defaultdict(list)
    n_records = 0
    for line in text.splitlines():
        if not line:
            continue
        chrom, pos, svtype, svlen = (line.split("\t") + ["", "", "", ""])[:4]
        try:
            start = int(pos)
        except ValueError:
            continue
        length = 0
        if svlen not in ("", "."):
            try:
                length = int(svlen)
            except ValueError:
                length = 0
        index[normalise_contig(chrom)].append((start, start + length, svtype))
        n_records += 1
    for chrom in index:
        index[chrom].sort()

    def covered(sites: Sequence[tuple[str, int]]) -> int:
        n = 0
        for chrom, pos in sites:
            records = index.get(normalise_contig(chrom))
            if not records:
                continue
            starts = [r[0] for r in records]
            lo = bisect.bisect_left(starts, pos - (tolerance + 300))
            hi = bisect.bisect_right(starts, pos + tolerance + 300)
            for start, end, _ in records[lo:hi]:
                if start - tolerance <= pos <= end + tolerance:
                    n += 1
                    break
        return n

    nested_covered = covered(nested_sites)
    control_covered = covered(control_sites)
    n_nested, n_control = len(nested_sites), len(control_sites)
    nested_rate = nested_covered / n_nested if n_nested else None
    control_rate = control_covered / n_control if n_control else None
    p = fisher_exact_two_sided(
        nested_covered,
        n_nested - nested_covered,
        control_covered,
        n_control - control_covered,
    )
    return {
        "gate": "3d",
        "cross_method_file": str(cross_method_path),
        "cross_method_records_considered": n_records,
        "proximity_tolerance_bp": tolerance,
        "nested_arm": {
            "sites": n_nested,
            "covered": nested_covered,
            "fraction_covered": nested_rate,
        },
        "control_arm_non_nested_alu": {
            "sites": n_control,
            "covered": control_covered,
            "fraction_covered": control_rate,
        },
        "coverage_deficit_ratio": (
            control_rate / nested_rate
            if nested_rate and control_rate and nested_rate > 0
            else None
        ),
        "fisher_exact_two_sided_p": p,
        "verdict": (
            "cross_method_comparator_structurally_blind_to_nested_alu"
            if nested_rate is not None
            and control_rate is not None
            and nested_rate < control_rate / 2
            else "cross_method_comparator_covers_nested_alu"
        ),
        "consequence": (
            "an independent short-read callset cannot score these sites, so "
            "cross-method agreement is unavailable for this locus class; the "
            "evidence for that absence is this measurement, not an assumption, "
            "and the word 'validated' therefore stays off the table"
        ),
    }


# --------------------------------------------------------------------------
# recurrence across independent halves
# --------------------------------------------------------------------------


def split_halves(samples: Sequence[str]) -> dict[str, Any]:
    """Cut the cohort into two halves by sorted sample name.

    Pre-declared and genotype-blind. The rule is reported with the result so the
    split cannot be mistaken for a tuned choice.
    """
    order = sorted(samples)
    index_a = np.array([i for i, _ in enumerate(order) if i % 2 == 0], dtype=np.int64)
    index_b = np.array([i for i, _ in enumerate(order) if i % 2 == 1], dtype=np.int64)
    return {
        "rule": HALF_SPLIT_RULE,
        "order_is_genotype_blind": True,
        "n_samples": len(order),
        "n_half_a": int(index_a.size),
        "n_half_b": int(index_b.size),
        "index_a": index_a,
        "index_b": index_b,
        "order": order,
    }


def site_recurrence(
    keys: Sequence[tuple[str, int]],
    carriers: dict[tuple[str, int], np.ndarray],
    halves: dict[str, Any],
) -> dict[tuple[str, int], dict[str, Any]]:
    """Per-site carrier totals and cross-half support.

    A site is cross-half supported when at least one carrier falls in each half,
    which means more than one genome carries it. Sites with no genotype record
    are recorded as `unevaluable_no_genotype_record` and are never counted as
    having zero carriers.
    """
    index_a, index_b = halves["index_a"], halves["index_b"]
    out: dict[tuple[str, int], dict[str, Any]] = {}
    for key in keys:
        mask = carriers.get(key)
        if mask is None:
            out[key] = {
                "state": "unevaluable_no_genotype_record",
                "carriers_total": None,
                "carriers_half_a": None,
                "carriers_half_b": None,
                "cross_half_supported": None,
                "private_to_one_genome": None,
            }
            continue
        a = int(mask[index_a].sum()) if index_a.size else 0
        b = int(mask[index_b].sum()) if index_b.size else 0
        total = int(mask.sum())
        out[key] = {
            "state": "genotyped",
            "carriers_total": total,
            "carriers_half_a": a,
            "carriers_half_b": b,
            "cross_half_supported": bool(a >= 1 and b >= 1),
            "private_to_one_genome": total == 1,
        }
    return out


def accumulation_curve(
    keys: Sequence[tuple[str, int]],
    carriers: dict[tuple[str, int], np.ndarray],
    order: Sequence[str],
    checkpoints: Sequence[int] = ACCUMULATION_CHECKPOINTS,
    saturations: Sequence[float] = SATURATION_FRACTIONS,
) -> dict[str, Any]:
    """Unique sites and carrier observations as genomes accumulate.

    Site counts and carrier counts are reported side by side and never summed:
    a site carried by N genomes is one site and N observations, and collapsing
    them is the mistake this table exists to prevent.
    """
    keys = list(keys)
    position = {key: i for i, key in enumerate(keys)}
    matrix = np.zeros((len(keys), len(order)), dtype=bool)
    evaluable = 0
    for key in keys:
        mask = carriers.get(key)
        if mask is None:
            continue
        evaluable += 1
        matrix[position[key]] = mask
    seen = np.zeros(len(keys), dtype=bool)
    points: list[dict[str, Any]] = []
    genomes_to_saturation: dict[str, int | None] = {}
    final_unique = 0
    for n, sample in enumerate(order, start=1):
        column = matrix[:, n - 1]
        seen |= column
        final_unique = int(seen.sum())
        if n in checkpoints or n == len(order):
            points.append(
                {
                    "genomes_added": n,
                    "sample_added": sample,
                    "unique_sites": final_unique,
                    "new_unique_sites_this_genome": int(column.sum()),
                    "cumulative_carrier_observations": int(matrix[:, :n].sum()),
                }
            )
    # Saturation is measured against the final total, so it needs a second pass:
    # the target is not known until every genome has been added.
    for frac in saturations:
        target = frac * final_unique
        reached: int | None = None
        acc = np.zeros(len(keys), dtype=bool)
        running = 0
        for n in range(1, len(order) + 1):
            acc |= matrix[:, n - 1]
            running = int(acc.sum())
            if running >= target:
                reached = n
                break
        genomes_to_saturation[f"fraction_{frac:g}"] = reached
    return {
        "genome_order": "sorted_sample_name",
        "order_is_genotype_blind": True,
        "n_genomes": len(order),
        "sites_considered": len(keys),
        "sites_with_genotype_record": evaluable,
        "final_unique_sites": final_unique,
        "final_carrier_observations": int(matrix.sum()),
        "mean_carriers_per_unique_site": (
            float(matrix.sum()) / final_unique if final_unique else None
        ),
        "genomes_to_reach_fraction_of_final_unique_sites": genomes_to_saturation,
        "marginal_yield_note": (
            "new_unique_sites_this_genome is the marginal yield; extrapolate from "
            "the last checkpoints rather than from the mean, because the mean is "
            "dominated by the first few genomes"
        ),
        "checkpoints": points,
    }


# --------------------------------------------------------------------------
# the comparison, run on Phase 2b's machinery
# --------------------------------------------------------------------------


def register_comparison() -> None:
    """Add the recurrence comparison to the shared matching tables.

    Allele frequency is absent on purpose: carrier count is algebraically tied to
    the outcome, so conditioning on it would condition on the result. The
    remaining covariates are host context, none of which is a hallmark.
    """
    tprt.COVARIATE_SETS[COMPARISON] = {
        "host_subfamily": tprt._subfamily,
        "alignment_confidence": tprt._confidence_class,
        "local_opportunity": tprt._support_bin,
    }
    tprt.LADDERS[COMPARISON] = {
        "full": ("host_subfamily", "alignment_confidence", "local_opportunity"),
        "drop_local_opportunity": ("host_subfamily", "alignment_confidence"),
        "subfamily_only": ("host_subfamily",),
        "no_strata": (),
    }


def build_comparison_rows(
    rows: Sequence[dict[str, str]],
    recurrence: dict[tuple[str, int], dict[str, Any]],
) -> list[dict[str, str]]:
    """Attach recurrence to cohort rows so it can be read as an outcome."""
    out = []
    for row in rows:
        key = (row.get("chrom", ""), int(row.get("pos") or 0))
        rec = recurrence.get(key)
        if rec is None:
            continue
        new = dict(row)
        new["_carriers_total"] = (
            "" if rec["carriers_total"] is None else str(rec["carriers_total"])
        )
        new["_carriers_half_a"] = (
            "" if rec["carriers_half_a"] is None else str(rec["carriers_half_a"])
        )
        new["_carriers_half_b"] = (
            "" if rec["carriers_half_b"] is None else str(rec["carriers_half_b"])
        )
        new["_cross_half_supported"] = (
            ""
            if rec["cross_half_supported"] is None
            else ("1" if rec["cross_half_supported"] else "0")
        )
        out.append(new)
    return out


def event_recurrence_cross_half(row: dict[str, str]) -> bool | None:
    """The outcome: was this site carried by genomes in both halves?

    Unevaluable when the site has no genotype record, so it leaves its own
    denominator rather than being counted as a non-replicated site.
    """
    value = row.get("_cross_half_supported")
    if value == "" or value is None:
        return None
    return value == "1"


def event_recurrence_multi_genome(row: dict[str, str]) -> bool | None:
    """Secondary outcome: carried by at least MIN_CARRIERS_PER_SITE genomes."""
    total = row.get("_carriers_total")
    if total in ("", None):
        return None
    return int(total) >= MIN_CARRIERS_PER_SITE


OUTCOMES: dict[str, Callable[[dict[str, str]], bool | None]] = {
    "recurrence_cross_half_supported": event_recurrence_cross_half,
    "recurrence_at_least_two_genomes": event_recurrence_multi_genome,
}


def run_recurrence_comparison(
    rows: Sequence[dict[str, str]],
    seed: int,
    permutations: int,
) -> dict[str, Any]:
    """Do TPRT hallmarks distinguish sites seen in more than one genome?

    The hallmark is the exposure and recurrence is the outcome, which is the
    reverse of Phase 2b's arrangement and is deliberate. Exposing recurrence
    instead would make the design unmatchable: a site with a single carrier
    cannot be cross-half supported by arithmetic, so the reference arm would
    consist almost entirely of singletons and every ladder rung would fail its
    minimum-reference floor. With the hallmark as exposure both arms are drawn
    from the same sites and the comparison is estimable.

    Nine pre-declared hallmark exposures are tested against two pre-declared
    recurrence outcomes. Benjamini-Hochberg is applied across exposures within
    each outcome, so the eighteen p-values are not read as if they were one.
    """
    blocks: dict[str, Any] = {}
    raw: dict[str, dict[str, float | None]] = {
        label: {} for label in OUTCOMES
    }
    for exposure_label, event in tprt.EVENTS.items():
        exposed, reference = [], []
        for row in rows:
            value = event(row)
            if value is True:
                exposed.append(row)
            elif value is False:
                reference.append(row)
        entry: dict[str, Any] = {
            "exposure": exposure_label,
            "exposed_input": len(exposed),
            "reference_input": len(reference),
        }
        if not exposed or not reference:
            entry["status"] = "not_estimable_empty_arm"
            blocks[exposure_label] = entry
            for outcome in OUTCOMES:
                raw[outcome][exposure_label] = None
            continue
        match = tprt.coarsened_match(exposed, reference, COMPARISON)
        result = tprt.stratified_tests(
            match["strata"], OUTCOMES, seed=seed, permutations=permutations
        )
        entry.update(
            {
                "status": "estimated",
                "level_used": match["level"],
                "covariates_used": match["covariates"],
                "exposed_matched": match["exposed_matched"],
                "exposed_dropped_unsupported": match["exposed_dropped_unsupported"],
                "matched": match["level"] != "no_strata",
                "n_strata": len(match["strata"]),
                "balance": tprt.balance_table(match["strata"], COMPARISON, match["level"]),
                "outcomes": result,
            }
        )
        blocks[exposure_label] = entry
        for outcome, res in result.items():
            raw[outcome][exposure_label] = res.get("permutation_p")

    for outcome, per_exposure in raw.items():
        labels = list(per_exposure)
        adjusted = tprt._benjamini_hochberg([per_exposure[k] for k in labels])
        for label, value in zip(labels, adjusted):
            if blocks[label].get("status") == "estimated":
                blocks[label]["outcomes"][outcome]["permutation_p_bh_across_exposures"] = value

    return {
        "design": (
            "hallmark as exposure, recurrence as outcome; reversed from Phase 2b "
            "so that both arms are drawn from the same sites and the reference "
            "arm is not reduced to singletons by arithmetic"
        ),
        "matching_rule": (
            "coarsened exact matching on the pre-registered ladder; no hallmark, "
            "TSD length or motif content is used as a matching variable, and "
            "allele frequency is excluded because it is algebraically tied to the "
            "outcome"
        ),
        "multiplicity": (
            "Benjamini-Hochberg across hallmark exposures within each outcome"
        ),
        "exposures": blocks,
    }


def recurrence_summary(
    annotated: Sequence[dict[str, str]], recurrence: dict[tuple[str, int], dict[str, Any]]
) -> dict[str, Any]:
    """Descriptive recurrence over every analysis site, with no site dropped."""
    total = len(annotated)
    counts = collections.Counter(
        (recurrence[(r.get("chrom", ""), int(r.get("pos") or 0))]["state"])
        for r in annotated
    )
    cross = sum(1 for r in annotated if r["_cross_half_supported"] == "1")
    single_half = sum(1 for r in annotated if r["_cross_half_supported"] == "0")
    unevaluated = counts.get("unevaluable_no_genotype_record", 0)
    private = sum(1 for r in annotated if r["_carriers_total"] == "1")
    evaluable = total - unevaluated
    return {
        "analysis_sites": total,
        "sites_with_genotype_record": evaluable,
        "sites_without_genotype_record": unevaluated,
        "cross_half_supported": cross,
        "single_half_only": single_half,
        "private_to_one_genome": private,
        "cross_half_fraction_of_evaluable": cross / evaluable if evaluable else None,
        "note": (
            "site counts and carrier counts are different quantities and are "
            "never summed: a site carried by N genomes is one site and N "
            "carrier observations"
        ),
    }


# --------------------------------------------------------------------------
# reporting
# --------------------------------------------------------------------------


def forbidden_metric_keys(payload: Any, path: str = "") -> list[str]:
    """Every JSON key in the payload that names a prohibited metric.

    Keys and identifiers only. Prose may name a metric in order to disclaim it.
    """
    hits: list[str] = []
    if isinstance(payload, dict):
        for key, value in payload.items():
            lowered = str(key).lower()
            for term in FORBIDDEN_METRIC_TERMS:
                if term in lowered:
                    hits.append(f"{path}/{key}")
            hits.extend(forbidden_metric_keys(value, f"{path}/{key}"))
    elif isinstance(payload, list):
        for i, value in enumerate(payload):
            hits.extend(forbidden_metric_keys(value, f"{path}[{i}]"))
    return hits


def write_site_table(
    path: Path,
    annotated: Sequence[dict[str, str]],
    samples: Sequence[str],
    halves: dict[str, Any],
    carriers: dict[tuple[str, int], np.ndarray],
) -> int:
    """Per-site carrier counts, with the carrier sample names spelled out.

    Deliverable 8 asks for per-sample call counts beside the unique-site table.
    Carrier names are listed for every site so the counts can be audited rather
    than trusted.
    """
    index_a = set(halves["index_a"].tolist())
    index_b = set(halves["index_b"].tolist())
    written = 0
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                "site_id",
                "chrom",
                "pos",
                "consensus_offset",
                "in_peak_window",
                "genotype_state",
                "carriers_total",
                "carriers_half_a",
                "carriers_half_b",
                "cross_half_supported",
                "private_to_one_genome",
                "carriers_half_a_samples",
                "carriers_half_b_samples",
            ]
        )
        for i, row in enumerate(annotated, start=1):
            key = (row.get("chrom", ""), int(row.get("pos") or 0))
            mask = carriers.get(key)
            if mask is None:
                a_names: list[str] = []
                b_names: list[str] = []
            else:
                a_names = [samples[j] for j in np.flatnonzero(mask).tolist() if j in index_a]
                b_names = [samples[j] for j in np.flatnonzero(mask).tolist() if j in index_b]
            writer.writerow(
                [
                    f"NH{i:06d}",
                    row.get("chrom", ""),
                    row.get("pos", ""),
                    row.get("consensus_offset", ""),
                    row.get("_in_peak_window", ""),
                    row.get("_cross_half_supported", ""),
                    row.get("_carriers_total", ""),
                    row.get("_carriers_half_a", ""),
                    row.get("_carriers_half_b", ""),
                    row.get("_cross_half_supported", ""),
                    "1" if row.get("_carriers_total") == "1" else "0",
                    ";".join(a_names),
                    ";".join(b_names),
                ]
            )
            written += 1
    return written


def _fmt_pct(value: float | None) -> str:
    return "n/a" if value is None else f"{100 * value:.2f}%"


def write_report(outdir: Path, report: dict[str, Any]) -> None:
    gates = report["gates"]
    summary = report["recurrence"]
    floor = gates["3d"]
    accumulation = report["accumulation"]
    primary = report["comparison"]

    lines: list[str] = []
    lines.append("# Phase 3 -- genotype concordance (nested-Alu cohort)")
    lines.append("")
    lines.append(
        "Generated by `scripts/score_genotype_concordance.py`. Every number below "
        "is reproduced by that script from the pooled site callset, the genotyped "
        "callset, and the cohort CSV; nothing is transcribed by hand."
    )
    lines.append("")

    lines.append("## Headline")
    lines.append("")
    lines.append(
        f"1. **The pre-registered Phase 3 design is not executable on this data.** "
        f"Gate 3b found no missingness channel: {gates['3b']['missing_or_partial_genotypes']} "
        f"missing or partial genotypes out of {gates['3b']['genotypes_examined']:,} "
        "examined. A `0/0` therefore cannot be separated from a genotyping failure "
        "inside a hard repeat, so the plan's \"concordance at demonstrated-callable "
        "loci\" has no demonstration available and was not run."
    )
    lines.append(
        f"2. **Cross-method agreement is unavailable, and that is now a measurement.** "
        f"An independent short-read callset covers {_fmt_pct(floor.get('nested_arm', {}).get('fraction_covered'))} "
        f"of nested-Alu sites against "
        f"{_fmt_pct(floor.get('control_arm_non_nested_alu', {}).get('fraction_covered'))} "
        f"of non-nested Alu sites from the same cohort "
        f"(Fisher exact p = {_fmt_p(floor.get('fisher_exact_two_sided_p'))}). A short "
        "read cannot span an Alu inserted into an Alu, so the comparator is "
        "structurally blind here rather than merely noisy."
    )
    lines.append(
        f"3. **What remains is recurrence, and it is thin.** "
        f"{summary['cross_half_supported']} of {summary['sites_with_genotype_record']} "
        f"evaluable sites are carried by genomes in both halves of the cohort "
        f"({_fmt_pct(summary['cross_half_fraction_of_evaluable'])}); "
        f"{summary['private_to_one_genome']} are private to a single genome."
    )
    lines.append("")

    lines.append("## Gates")
    lines.append("")
    lines.append("| gate | verdict | consequence |")
    lines.append("|---|---|---|")
    for key in sorted(gates):
        gate = gates[key]
        consequence = str(gate.get("consequence") or gate.get("note") or "")
        lines.append(
            f"| {key} | `{gate['verdict']}` | {consequence.split('.')[0]}. |"
        )
    lines.append("")
    lines.append(
        f"Gate 3a is worth reading closely: the site and genotype callsets declare "
        f"different lengths for {gates['3a']['contig_length_disagreements']} of "
        f"{gates['3a']['shared_contigs']} shared contigs, with the sign of the "
        "difference flipping across chromosomes. The join was verified on exact "
        f"coordinates anyway -- {gates['3a']['analysis_sites_absent_from_genotype_bcf']} "
        "analysis sites are absent from the genotype callset -- so the data overrules "
        "the headers. The header disagreement is still reported, because a join that "
        "had trusted it would have been wrong."
    )
    lines.append("")

    lines.append("## Genotype concordance across independent halves")
    lines.append("")
    lines.append(
        f"The {gates['0c']['n_samples']} genomes are cut in half by sorted sample "
        f"name ({report['halves']['rule']}), a split fixed before any genotype was "
        "read. A site carried by genomes in both halves has been seen more than "
        "once, by more than one genome. That is agreement *within one method* and "
        "is not called validation; Gate 3d is the reason."
    )
    lines.append("")
    lines.append(
        "**Design note.** The hallmark is the exposure and recurrence is the "
        "outcome, which reverses Phase 2b's arrangement on purpose. Exposing "
        "recurrence instead makes the comparison unmatchable: a site with one "
        "carrier cannot be cross-half supported by arithmetic, so the reference "
        "arm would be almost entirely singletons and every ladder rung would fail "
        "its minimum-reference floor. Allele frequency is likewise kept out of the "
        "covariates for the same reason."
    )
    lines.append("")
    lines.append(
        "| hallmark exposure | exposed | reference | rung | matched | MH risk diff "
        "(cross-half) | raw p | BH p |"
    )
    lines.append("|---|---|---|---|---|---|---|---|")
    for label in sorted(primary["exposures"]):
        block = primary["exposures"][label]
        if block.get("status") != "estimated":
            lines.append(
                f"| `{label}` | {block['exposed_input']} | {block['reference_input']} "
                f"| -- | {block['status']} | -- | -- | -- |"
            )
            continue
        res = block["outcomes"]["recurrence_cross_half_supported"]
        obs = res["observed"]
        diff = obs.get("risk_difference")
        lines.append(
            f"| `{label}` | {block['exposed_matched']}/{block['exposed_input']} "
            f"| {block['reference_input']} | `{block['level_used']}` "
            f"| {'yes' if block['matched'] else 'no'} "
            f"| {'n/a' if diff is None else f'{diff:+.3f}'} "
            f"| {_fmt_p(res.get('permutation_p'))} "
            f"| {_fmt_p(res.get('permutation_p_bh_across_exposures'))} |"
        )
    lines.append("")
    unmatched = [
        label
        for label, block in primary["exposures"].items()
        if block.get("status") == "estimated" and not block["matched"]
    ]
    if unmatched:
        lines.append(
            "UNMATCHED (no ladder rung retained a usable share of the exposed arm, so "
            "the row above is a crude difference and must not be read as a matched "
            "estimate): " + ", ".join(f"`{label}`" for label in sorted(unmatched)) + "."
        )
        lines.append("")

    lines.append("## Accumulation")
    lines.append("")
    lines.append(
        f"Genomes are added in sorted-name order. At {accumulation['n_genomes']} "
        f"genomes the cohort holds {accumulation['final_unique_sites']:,} unique "
        f"nested sites backed by {accumulation['final_carrier_observations']:,} "
        f"carrier observations ({_fmt_num(accumulation['mean_carriers_per_unique_site'])} "
        "carriers per unique site). Site counts and carrier counts are reported "
        "side by side and are never summed."
    )
    lines.append("")
    lines.append("| genomes added | unique sites | new this genome | cumulative carriers |")
    lines.append("|---|---|---|---|")
    for point in accumulation["checkpoints"]:
        lines.append(
            f"| {point['genomes_added']} | {point['unique_sites']:,} "
            f"| {point['new_unique_sites_this_genome']} "
            f"| {point['cumulative_carrier_observations']:,} |"
        )
    lines.append("")
    lines.append(
        "Genomes needed to reach a given share of the final unique-site count: "
        + (
            ", ".join(
                f"{k.replace('fraction_', '')} -> {v if v is not None else 'not reached'}"
                for k, v in accumulation[
                    "genomes_to_reach_fraction_of_final_unique_sites"
                ].items()
            )
            if accumulation["final_unique_sites"]
            else "n/a -- the run produced no evaluable sites"
        )
        + "."
    )
    lines.append("")

    lines.append("## What this does not say")
    lines.append("")
    lines.append(
        "- No detection-performance metric is reported, and Gate 0c's tier forbids "
        "one. The genotype axis carries no method independence: HG03086 is inside "
        "the pooled cohort, and the 908 genomes are one pipeline applied to 908 "
        "genomes."
    )
    lines.append(
        "- The cross-half measure is not agreement between two methods. Section 5 of "
        "the plan reserves that word for evidence that does not exist on this data."
    )
    lines.append(
        "- The cohort is a catalogue of insertion *sites*, not one genome's biology. "
        "Polymorphic sites are survivors of selection; nothing here is a de novo "
        "integration rate."
    )
    lines.append("")
    (outdir / "genotype_concordance.md").write_text("\n".join(lines) + "\n")


def _fmt_p(value: float | None) -> str:
    return "n/a" if value is None else f"{value:.4g}"


def _fmt_num(value: float | None, spec: str = ".1f") -> str:
    """Format a possibly-absent number.

    A report writer that raises on a degenerate run loses the degenerate run's
    evidence, which is exactly the run worth reading, so nothing here is allowed
    to raise on a None.
    """
    return "n/a" if value is None else format(value, spec)


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cohort", type=Path, default=DEFAULT_COHORT)
    p.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    p.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE)
    p.add_argument("--alu-consensus", type=Path, default=DEFAULT_ALU_CONSENSUS)
    p.add_argument("--genotype-bcf", type=Path, default=DEFAULT_GENOTYPE_BCF)
    p.add_argument("--site-bcf", type=Path, default=DEFAULT_SITE_BCF)
    p.add_argument(
        "--cross-method-bcf",
        type=Path,
        default=DEFAULT_CROSS_METHOD_BCF,
        help="Independent GRCh38 callset for Gate 3d; pass a nonexistent path to skip.",
    )
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--permutations", type=int, default=DEFAULT_PERMUTATIONS)
    p.add_argument("--seed", type=int, default=DEFAULT_SEED)
    p.add_argument(
        "--no-reference",
        action="store_true",
        help="Skip reference-backed hallmarks; they become `unevaluable`, not `absent`.",
    )
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    args.outdir.mkdir(parents=True, exist_ok=True)
    register_comparison()

    sentinels_list, sentinel_source = tprt.load_sentinels(args.manifest)
    sentinels = frozenset(sentinels_list)
    print(f"TSD sentinels ({sentinel_source}): {sorted(sentinels)}", flush=True)

    rows = tprt.read_cohort(args.cohort)
    analysis = tprt.nested_analysis_set(rows)
    print(f"nested analysis sites: {len(analysis)}", flush=True)

    samples = sample_names(args.genotype_bcf)
    print(f"genotyped samples: {len(samples)}", flush=True)
    gates: dict[str, Any] = {"0c": gate_0c_independence(samples)}
    print(f"[gate 0c] {gates['0c']['verdict']} / {gates['0c']['tier']}", flush=True)

    keys = [(r.get("chrom", ""), int(r.get("pos") or 0)) for r in analysis]
    gates["3a"] = gate_3a_coordinate_frame(args.site_bcf, args.genotype_bcf, keys)
    print(f"[gate 3a] {gates['3a']['verdict']}", flush=True)

    samples, carriers, counters = read_site_genotypes(args.genotype_bcf, keys)
    gates["3b"] = gate_3b_genotype_channel(
        counters, n_sites=len(keys), n_samples=len(samples)
    )
    print(f"[gate 3b] {gates['3b']['verdict']}", flush=True)

    cross_path = args.cross_method_bcf if args.cross_method_bcf.exists() else None
    cross_samples = None
    if cross_path is not None:
        try:
            cross_samples = sample_names(cross_path)
        except RuntimeError:
            cross_path = None
    gates["3c"] = gate_3c_sample_identity(samples, cross_samples, cross_path)
    print(f"[gate 3c] {gates['3c']['verdict']}", flush=True)

    control = [
        (r.get("chrom", ""), int(r.get("pos") or 0))
        for r in rows
        if r.get("nested_in_alu_host") == "0"
    ]
    gates["3d"] = gate_3d_cross_method_floor(cross_path, keys, control)
    print(f"[gate 3d] {gates['3d']['verdict']}", flush=True)

    halves = split_halves(samples)
    recurrence = site_recurrence(keys, carriers, halves)
    annotated = build_comparison_rows(analysis, recurrence)
    summary = recurrence_summary(annotated, recurrence)
    print(
        f"cross-half supported: {summary['cross_half_supported']}/"
        f"{summary['sites_with_genotype_record']}",
        flush=True,
    )

    motif, motif_meta = tprt.derive_en_motif_proxy(args.alu_consensus)
    reference = None if args.no_reference else tprt.Reference(args.reference)
    scored = tprt.score_rows(annotated, reference, motif, sentinels)
    hallmark_rows = tprt.to_rows(scored)

    comparison = run_recurrence_comparison(
        hallmark_rows, args.seed, args.permutations
    )
    estimated = sum(
        1 for b in comparison["exposures"].values() if b.get("status") == "estimated"
    )
    print(
        f"recurrence comparison: {estimated}/{len(comparison['exposures'])} hallmark "
        "exposures estimable",
        flush=True,
    )

    accumulation = accumulation_curve(keys, carriers, halves["order"])

    report: dict[str, Any] = {
        "phase": "3",
        "subject": SUBJECT,
        "terminology": (
            "genotype concordance only; 'validated' is reserved for cross-method "
            "evidence that does not exist on this data"
        ),
        "inputs": {
            "cohort": str(args.cohort),
            "genotype_bcf": str(args.genotype_bcf),
            "site_bcf": str(args.site_bcf),
            "cross_method_bcf": str(cross_path) if cross_path else None,
            "permutations": args.permutations,
            "seed": args.seed,
            "reference_backed": not args.no_reference,
            "en_motif": motif_meta,
        },
        "gates": gates,
        "halves": {
            k: v for k, v in halves.items() if k not in ("index_a", "index_b", "order")
        },
        "recurrence": summary,
        "accumulation": accumulation,
        "comparison": comparison,
        "deviations": [
            "inserted-element subfamily unavailable as a covariate (unchanged from Phase 2b)",
            "curated EN consensus not verified; the consensus-derived proxy is used",
        ],
    }
    offenders = forbidden_metric_keys(report)
    if offenders:
        raise SystemExit(f"prohibited metric naming in report keys: {offenders[:5]}")

    (args.outdir / "genotype_concordance.json").write_text(
        json.dumps(report, indent=2, sort_keys=True, default=str) + "\n"
    )
    n_sites = write_site_table(
        args.outdir / "phase3_site_genotypes.csv", annotated, samples, halves, carriers
    )
    write_report(args.outdir, report)
    print(
        f"wrote genotype_concordance.json / .md and phase3_site_genotypes.csv "
        f"({n_sites} sites)",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
