"""Phase 2b: do the TPRT hallmarks distinguish nested-Alu calls at all?

Phases 1 and 2a established a positional fact and removed one artefact
explanation. Neither says anything about mechanism. A catalogue enrichment can
be real and still be explained by something other than target-primed reverse
transcription (TPRT) of the inserted element by host machinery. This phase asks
the mechanistic question: do calls carrying more hallmarks look different from
matched calls carrying fewer, in a way a mechanism would produce?

Three hallmarks, each scored independently:

  1. TSD, checked against the reference rather than taken on trust. The caller
     reports `TSD_SEQ` and `TSD_LEN`. A reported TSD is not a resolved TSD, so
     it is compared base-for-base against the reference sequence the insertion
     point implies, and the number and *positions* of the disagreements are
     kept. This arm turns out to be a measurement of caller/reference
     consistency rather than of mechanism, and the report says so with numbers
     attached instead of quietly scoring it as if it were evidence.
  2. A 3' poly(A) tail, from the caller-reported `POLYA_SEQ`/`POLYA_LEN`.
  3. An EN motif at the integration nick, strand-compatible. The Alu's 5' end
     carries the endonuclease recognition motif, so the motif must sit at the
     side of the nick the Alu's 5' end maps to. A plus-strand insert puts it on
     the reference plus strand to the left of the nick; a minus-strand insert
     puts its reverse complement on the plus strand to the right. A motif found
     only on the wrong side is `motif_incompatible`, not `motif_absent`.

Four rules constrain the scoring, each fixed before any result was read:

  * `unevaluable` is not `absent`. Missing sequence, an unmapped breakpoint or
    a sentinel TSD length make a component *unevaluable*, and unevaluable rows
    leave the denominator of every rate they feed. Coding them `absent` would
    silently manufacture evidence against the mechanism.
  * TSD length is never treated as continuous. Gate 0d found sentinel pile-ups
    at 41 bp (120.9x local median) and 82 bp (21.3x), 10.8% of calls. Those
    values are excluded from every TSD-length statement and their TSD component
    is `unevaluable`. Length is used only as a category: de novo range or longer.
  * Matching never uses a hallmark. Two comparisons are run, each matched on
    covariates chosen before the outcomes were seen. Matching on TSD length,
    poly(A) or motif content would be matching on the answer.
  * The pure-A/T TSD fraction is reported only against the sequence- and
    callability-matched control, because a TSD that is pure A/T is partly a
    statement about where the sequence was, not only about what happened.

The reference cannot reconstruct a TSD from both junctions, and it is worth
being blunt about why. A TSD is a duplication of a pre-existing target, so in
the *sample* allele the target appears twice, once on each side of the inserted
element. In the *reference*, which lacks the inserted element, it appears
exactly once. There is no second copy in the reference to compare against, and
no sample reads are available here, so a true two-junction reconstruction is
not possible from the inputs this project has. What is possible is a
consistency check: does the reported TSD sequence agree with the reference at
the insertion point. Measured on this cohort, that check mostly fails -- see
`tsd_verification` in the report, which is reported as a finding rather than
folded into a score.

The EN motif is a *consensus-derived proxy*: the 5' terminal window of the Alu
repeat consensus, matched with a small mismatch budget. A curated EN consensus
could not be verified from a primary source here, so the proxy is labelled as
one everywhere it appears and `--en-motif-seq` overrides it.

That proxy carries a confound which must be stated before any of its numbers are
read. Its sequence is Alu sequence, so a site inside an Alu is inside sequence
that is more likely to match it, for reasons of composition alone. The
nested-vs-unnested comparison is confounded by construction -- the exposed arm
is inside Alus by definition and the reference arm is not -- so its EN-motif
result is reported as confounded and is not evidence about mechanism. The
peak-vs-rest-of-host comparison does not have this problem: both arms lie
inside the same Alu hosts, so composition is held roughly constant by the
design rather than by an adjustment.

Because of that, a second, non-Alu-derived probe is reported alongside it: the
L1 ORF2p pentamer 5'-TTTT/A-3' adjacent to the nick. It is a different sequence
with a different bias, and agreement or disagreement between the two probes is
more informative than either alone.

Both comparisons use coarsened exact matching. Strata are built from the
covariate set, and if the finest rung leaves exposed sites without reference
support the rung is coarsened along a fixed ladder. The rung used and the
exposed sites lost to it are both reported; a comparison that silently discards
most of the peak is worse than one that says it did.

Inference is exact-conditional rather than asymptotic. Within each stratum the
two group sizes are held fixed and labels are permuted, which is the standard
randomisation null and needs no distributional assumption; the statistic is the
Mantel-Haenszel pooled risk difference. No scipy: the permutation is vectorised
over rows with numpy.

The primary comparison is peak window against the rest of the host. The
mid-host region used as a control in Phase 2a holds only 185 rows, which cannot
support any stratification finer than a single covariate; it is kept as a
descriptive arm, and the peak-vs-rest comparison carries the inference.

This phase measures whether hallmark content *separates* the groups. It does
not establish that separation is causal, and it licenses no detection-accuracy
claim.
"""

from __future__ import annotations

import argparse
import collections
import csv
import json
import sys
from pathlib import Path
from typing import Any, Callable, Sequence

import numpy as np

REPO_SRC = Path(__file__).resolve().parents[1] / "src"
if str(REPO_SRC) not in sys.path:
    sys.path.insert(0, str(REPO_SRC))

RESULTS_DIR = Path(__file__).resolve().parents[2] / "nested_analysis" / "results_longread"
DEFAULT_COHORT = RESULTS_DIR / "per_call_longread_nested.csv"
DEFAULT_MANIFEST = RESULTS_DIR / "provenance_manifest.json"
DEFAULT_REFERENCE = (
    Path("~/retrotransposon-workdir/data/public/reference/hg38/Homo_sapiens_assembly38.fasta").expanduser()
)
DEFAULT_ALU_CONSENSUS = (
    Path(
        "~/retrotransposon-workdir/data/public/retrotransposon_db/ucsc_repeatbrowser/hg38reps.fa"
    ).expanduser()
)

#: The pre-registered peak window, restated from Phase 1 so a Phase 2b run
#: cannot silently drift onto a different peak.
PEAK_WINDOW = (128, 139)  # consensus 133 +/- 5, half-open
#: Mid-host control. Chosen in Phase 2a and not re-chosen here.
MID_HOST_CONTROL = (160, 260)

NEAR_FULL_HOST = (280, 320)

#: Gate 0d sentinel values, used only as a fallback when no manifest is given.
FALLBACK_SENTINELS = (41, 82)

#: de novo TSD range asserted by Phase 2a, valid only off the sentinel set.
DE_NOVO_TSD_RANGE = (5, 27)
POLYA_MIN_LEN = 8
#: Mismatch budget for calling a reported TSD corroborated by the reference.
TSD_MAX_MISMATCH = 2
#: Secondary poly(A) threshold: the primary one is satisfied by nearly every
#: call, so it cannot separate anything and a second cut is reported too.
POLYA_SECONDARY_MIN_LEN = 20

EN_MOTIF_LENGTH_BP = 15
EN_MOTIF_MAX_MISMATCH = 2
EN_MOTIF_FLANK_BP = 40

#: L1 ORF2p nicks the bottom strand at 5'-TTTT/A-3'. This probe is *not*
#: derived from Alu sequence, so it does not share the composition confound of
#: the Alu-5'-window motif above, and the two probes are reported side by side
#: for exactly that reason.
L1_EN_PENTAMER = "TTTTA"
L1_EN_SIDE_BP = 5

MIN_REFERENCE_PER_STRATUM = 20
MIN_EXPOSED_RETAINED_FRACTION = 0.50
DEFAULT_PERMUTATIONS = 10_000
DEFAULT_SEED = 20261001

EVENT_TRUE = np.int8(1)
EVENT_FALSE = np.int8(0)
EVENT_UNEVALUABLE = np.int8(-1)


# --------------------------------------------------------------------------
# inputs
# --------------------------------------------------------------------------


def read_cohort(path: Path) -> list[dict[str, str]]:
    with Path(path).open(newline="") as fh:
        return list(csv.DictReader(fh))


def _to_int(value: str | None) -> int | None:
    if value is None or value in ("", "."):
        return None
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def _to_float(value: str | None) -> float | None:
    if value is None or value in ("", "."):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def load_sentinels(manifest: Path | None) -> tuple[list[int], str]:
    if manifest and Path(manifest).exists():
        man = json.loads(Path(manifest).read_text())
        found = man.get("gate_0d_tsd_sentinels", {}).get("sentinel_values", {})
        if found:
            return sorted(int(k) for k in found), "manifest"
    return list(FALLBACK_SENTINELS), "fallback"


def read_fasta(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    name: str | None = None
    chunks: list[str] = []
    with Path(path).open() as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    out[name] = "".join(chunks)
                name = line[1:].strip().split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        out[name] = "".join(chunks)
    return out


def derive_en_motif_proxy(consensus_path: Path) -> tuple[str, dict[str, Any]]:
    """Take the 5' terminal window of the Alu consensus as the EN motif proxy.

    The proxy is the modal base at each of the first EN_MOTIF_LENGTH_BP columns
    across every Alu consensus, not one chosen subfamily, so a single divergent
    AluJ consensus cannot define the motif on its own.
    """
    seqs = {k: v for k, v in read_fasta(consensus_path).items() if k.startswith("Alu")}
    usable = [s for s in seqs.values() if len(s) >= EN_MOTIF_LENGTH_BP + 10]
    if not usable:
        raise SystemExit(f"no usable Alu consensus entries in {consensus_path}")
    motif_parts: list[str] = []
    support: list[float] = []
    for i in range(EN_MOTIF_LENGTH_BP):
        col = collections.Counter(s[i].upper() for s in usable)
        base, n = col.most_common(1)[0]
        motif_parts.append(base)
        support.append(n / sum(col.values()))
    return "".join(motif_parts), {
        "n_alu_consensus_entries": len(usable),
        "length_bp": EN_MOTIF_LENGTH_BP,
        "motif": "".join(motif_parts),
        "max_mismatch": EN_MOTIF_MAX_MISMATCH,
        "flank_bp": EN_MOTIF_FLANK_BP,
        "per_column_support_min": min(support),
        "per_column_support_median": sorted(support)[len(support) // 2],
        "status": "consensus_derived_proxy_not_curated_en_consensus",
    }


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


class Reference:
    """Thin pyfaidx wrapper that reports unevaluable rather than raising."""

    def __init__(self, path: Path) -> None:
        import pyfaidx

        self._fasta = pyfaidx.Fasta(str(path), as_raw=True, sequence_always_upper=True)
        self._cache: dict[tuple[str, int, int], str | None] = {}

    def fetch(self, chrom: str, start0: int, end0: int) -> str | None:
        """Half-open 0-based fetch, or None if the region is not readable."""
        if start0 < 0 or end0 <= start0:
            return None
        key = (chrom, start0, end0)
        if key not in self._cache:
            self._cache[key] = self._fetch_uncached(chrom, start0, end0)
        return self._cache[key]

    def _fetch_uncached(self, chrom: str, start0: int, end0: int) -> str | None:
        try:
            rec = self._fasta[chrom]
        except (KeyError, ValueError, TypeError):
            return None
        if end0 > len(rec):
            return None
        try:
            return str(rec[start0:end0])
        except (KeyError, IndexError, ValueError):
            return None


# --------------------------------------------------------------------------
# hallmark components
# --------------------------------------------------------------------------


def _hamming(a: str, b: str) -> int:
    if len(a) != len(b):
        return max(len(a), len(b))
    return sum(1 for x, y in zip(a, b) if x != y)


def tsd_reference_check(
    row: dict[str, str],
    reference: Reference | None,
    sentinels: frozenset[int],
) -> dict[str, Any]:
    """Compare the caller-reported TSD sequence against the reference at POS.

    `TSD_SEQ` is read back from the caller's own reads, and the reference does
    not contain the duplicated copy, so this can only ever corroborate, never
    reconstruct. The mismatch *positions* are kept because their distribution
    is the evidence: mismatches bunched at the terminal bases say the caller's
    TSD ran into the neighbouring sequence, which is a different failure from
    mismatches scattered along its length.
    """
    tsd_len = _to_int(row.get("tsd_len"))
    tsd_seq = (row.get("tsd_seq") or "").strip().upper()
    pos = _to_int(row.get("pos"))
    out: dict[str, Any] = {
        "tsd_len": tsd_len,
        "tsd_seq": tsd_seq,
        "state": "unevaluable",
        "reason": None,
        "n_mismatches": None,
        "exact": None,
        "mismatch_offsets_from_3prime_end": None,
    }
    if tsd_len is None or not tsd_seq:
        out["reason"] = "no_caller_tsd"
        return out
    if tsd_len in sentinels:
        out["reason"] = "sentinel_tsd_length"
        return out
    if reference is None or pos is None:
        out["reason"] = "no_reference_available"
        return out
    if len(tsd_seq) != tsd_len:
        out["reason"] = "tsd_seq_length_disagrees_with_tsd_len"
        return out
    cut = pos - 1
    window = reference.fetch(row["chrom"], cut, cut + tsd_len)
    if window is None:
        out["reason"] = "insertion_point_unreadable"
        return out
    mismatches = [i for i, (x, y) in enumerate(zip(window, tsd_seq)) if x != y]
    out["n_mismatches"] = len(mismatches)
    out["exact"] = not mismatches
    out["mismatch_offsets_from_3prime_end"] = [tsd_len - 1 - i for i in mismatches]
    out["state"] = (
        "reference_verifiable"
        if len(mismatches) <= TSD_MAX_MISMATCH
        else "reference_inconsistent"
    )
    return out


#: Retention of the caller's reported TSD, checked against the reference.
#: Backwards-compatible alias for the pre-rename name.
tsd_junction_state = tsd_reference_check


def tsd_length_class(tsd_len: int | None, sentinels: frozenset[int]) -> str:
    """Categorical only. Gate 0d forbids treating TSD length as continuous."""
    if tsd_len is None:
        return "unevaluable_no_length"
    if tsd_len in sentinels:
        return "sentinel_pileup"
    lo, hi = DE_NOVO_TSD_RANGE
    return "de_novo_range_5_27" if lo <= tsd_len <= hi else "longer_than_de_novo_28_plus"


def polya_state(row: dict[str, str]) -> dict[str, Any]:
    seq = (row.get("polya_seq") or "").strip().upper()
    declared = _to_int(row.get("polya_len"))
    if not seq:
        return {"state": "unevaluable", "reason": "no_caller_polya_field", "polya_len": None}
    if declared is None:
        return {"state": "unevaluable", "reason": "polya_len_missing", "polya_len": None}
    body = seq[:declared] if 0 < declared <= len(seq) else seq
    return {
        "state": "present" if len(body) >= POLYA_MIN_LEN else "absent",
        "reason": None,
        "polya_len": declared,
        "polya_a_fraction": body.count("A") / len(body) if body else 0.0,
    }


def en_motif_state(
    row: dict[str, str],
    reference: Reference | None,
    motif: str,
    flank: int = EN_MOTIF_FLANK_BP,
    max_mismatch: int = EN_MOTIF_MAX_MISMATCH,
) -> dict[str, Any]:
    """Search for the EN motif proxy at the nick, on the side the Alu's 5' end maps to.

    The Alu's 5' end carries the motif, so strand predicts the side: a plus
    insert puts the motif on the reference plus strand to the *left* of the
    insertion point, a minus insert puts its reverse complement on the plus
    strand to the *right*. Finding it only on the other side is a real
    incompatibility and is reported as such rather than as absence. A hit that
    straddles the nick is assigned by its start position, which is stated in the
    disclosures because it is the one place the two windows overlap.
    """
    strand = (row.get("insert_strand") or "").strip()
    out: dict[str, Any] = {
        "state": "unevaluable",
        "reason": None,
        "insert_strand": strand,
        "expected_side": None,
        "hit_side": None,
        "mismatches": None,
    }
    if strand not in ("+", "-"):
        out["reason"] = "insert_strand_unresolved"
        return out
    pos = _to_int(row.get("pos"))
    if reference is None or pos is None:
        out["reason"] = "no_reference_available"
        return out
    cut = pos - 1
    lo = cut - flank - len(motif)
    window = reference.fetch(row["chrom"], lo, cut + flank + len(motif))
    if window is None:
        out["reason"] = "insertion_point_unreadable"
        return out
    nick_index = flank + len(motif)
    out["expected_side"] = "left" if strand == "+" else "right"

    def _side(start: int) -> str:
        return "left" if start < nick_index else "right"

    compatible: tuple[int, int] | None = None
    incompatible: tuple[int, int] | None = None
    for pattern in (motif, revcomp(motif)):
        hit = _best_hit(window, pattern, max_mismatch)
        if hit is None:
            continue
        start, mism = hit
        if _side(start) == out["expected_side"]:
            if compatible is None or mism < compatible[1]:
                compatible = (start, mism)
        elif incompatible is None or mism < incompatible[1]:
            incompatible = (start, mism)
    if compatible is not None:
        out["state"] = "motif_compatible"
        out["hit_side"] = out["expected_side"]
        out["mismatches"] = compatible[1]
    elif incompatible is not None:
        out["state"] = "motif_incompatible"
        out["hit_side"] = "left" if out["expected_side"] == "right" else "right"
        out["mismatches"] = incompatible[1]
    else:
        out["state"] = "motif_absent"
    return out


def l1_en_nick_state(
    row: dict[str, str],
    reference: Reference | None,
    side_bp: int = L1_EN_SIDE_BP,
) -> dict[str, Any]:
    """Look for the L1 endonuclease 5'-TTTT/A-3' pentamer adjacent to the nick.

    Reported on both sides of the insertion point without claiming which side
    the inserted element's 5' end maps to: that mapping depends on the strand
    convention of the caller, and guessing it would bake an assumption into a
    probe whose whole purpose is to be convention-free. The per-side states are
    kept so a downstream reader can impose whichever convention they trust.
    """
    out: dict[str, Any] = {
        "state": "unevaluable",
        "reason": None,
        "five_prime_side": None,
        "three_prime_side": None,
    }
    pos = _to_int(row.get("pos"))
    if reference is None or pos is None:
        out["reason"] = "no_reference_available"
        return out
    cut = pos - 1
    five = reference.fetch(row["chrom"], cut - side_bp, cut)
    three = reference.fetch(row["chrom"], cut, cut + side_bp)
    if five is None or three is None:
        out["reason"] = "insertion_point_unreadable"
        return out
    out["five_prime_side"] = five == L1_EN_PENTAMER
    out["three_prime_side"] = three == L1_EN_PENTAMER
    if five == L1_EN_PENTAMER and three == L1_EN_PENTAMER:
        out["state"] = "tttt_a_both_sides"
    elif five == L1_EN_PENTAMER:
        out["state"] = "tttt_a_five_prime_only"
    elif three == L1_EN_PENTAMER:
        out["state"] = "tttt_a_three_prime_only"
    else:
        out["state"] = "no_tttt_a_at_nick"
    return out


def _best_hit(window: str, motif: str, max_mismatch: int) -> tuple[int, int] | None:
    best: tuple[int, int] | None = None
    for i in range(len(window) - len(motif) + 1):
        mism = _hamming(window[i : i + len(motif)], motif)
        if mism <= max_mismatch and (best is None or mism < best[1]):
            best = (i, mism)
    return best


def score_triad(components: dict[str, dict[str, Any]]) -> dict[str, Any]:
    """Combine the three components without ever collapsing unevaluable to absent.

    `n_observed / n_evaluable` is the rate to compare between groups. The
    boolean `triad_complete` is only true when all three components are
    evaluable, so a row with one measurable hallmark cannot pad the numerator
    of a three-hallmark count.
    """
    tsd, polya, en = components["tsd"], components["polya"], components["en"]
    evaluable = {
        "tsd": tsd["state"] in ("reference_verifiable", "reference_inconsistent"),
        "polya": polya["state"] in ("present", "absent"),
        "en": en["state"] in ("motif_compatible", "motif_incompatible", "motif_absent"),
    }
    observed = {
        "tsd": tsd["state"] == "reference_verifiable",
        "polya": polya["state"] == "present",
        "en": en["state"] == "motif_compatible",
    }
    n_eval = sum(evaluable.values())
    n_obs = sum(observed[k] for k in observed if evaluable[k])
    return {
        "evaluable": evaluable,
        "observed": observed,
        "n_evaluable": n_eval,
        "n_observed": n_obs,
        "triad_complete": n_eval == 3,
        "triad_complete_observed": n_eval == 3 and n_obs == 3,
        "hallmarks_per_evaluable": (n_obs / n_eval) if n_eval else None,
    }


# --------------------------------------------------------------------------
# analysis sets
# --------------------------------------------------------------------------


def nested_analysis_set(rows: Sequence[dict[str, str]]) -> list[dict[str, str]]:
    """The Phase 1/2a analysis set: nested, near-full host, unambiguous offset.

    Reused verbatim so a Phase 2b number is never computed on a different
    denominator than the Phase 1 number printed beside it.
    """
    out = []
    for row in rows:
        if row.get("nested_in_alu_host") != "1":
            continue
        host_len = _to_int(row.get("host_len"))
        if host_len is None or not (NEAR_FULL_HOST[0] <= host_len <= NEAR_FULL_HOST[1]):
            continue
        if row.get("consensus_mapping") != "unambiguous":
            continue
        if _to_int(row.get("consensus_offset")) is None:
            continue
        out.append(row)
    return out


def split_by_offset(
    rows: Sequence[dict[str, str]]
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    peak, control = [], []
    for row in rows:
        off = _to_int(row.get("consensus_offset"))
        if off is None:
            continue
        if PEAK_WINDOW[0] <= off < PEAK_WINDOW[1]:
            peak.append(row)
        elif MID_HOST_CONTROL[0] <= off < MID_HOST_CONTROL[1]:
            control.append(row)
    return peak, control


# --------------------------------------------------------------------------
# matching
# --------------------------------------------------------------------------


def _subfamily(row: dict[str, str]) -> str:
    name = row.get("host_name") or row.get("consensus_match_name") or ""
    return name.split("_")[0] if name else "unknown"


def _confidence_class(row: dict[str, str]) -> str:
    ident = _to_float(row.get("alignment_identity"))
    if ident is None:
        return "unknown"
    if ident < 0.55:
        return "low"
    return "mid" if ident < 0.80 else "high"


def _truncation_class(row: dict[str, str]) -> str:
    return "truncated" if (row.get("conformation") or "").startswith("TRUN") else "full_or_other"


def _af_bin(row: dict[str, str]) -> str:
    af = _to_float(row.get("site_allele_freq"))
    if af is None:
        return "unknown"
    if af < 0.001:
        return "lt_0.001"
    if af < 0.01:
        return "0.001_0.01"
    return "0.01_0.1" if af < 0.1 else "ge_0.1"


def _support_bin(row: dict[str, str]) -> str:
    pct = _to_float(row.get("perc_resolved"))
    if pct is None:
        return "unknown"
    if pct >= 99.0:
        return "ge_99"
    return "90_99" if pct >= 90.0 else "lt_90"


#: Covariates per comparison. None of these is a hallmark, a TSD length or a
#: motif: they are properties of the call and of its host context that a
#: mechanical explanation could also produce. Matching on any hallmark would be
#: matching on the outcome.
COVARIATE_SETS: dict[str, dict[str, Callable[[dict[str, str]], str]]] = {
    "peak_vs_rest_of_host": {
        "host_subfamily": _subfamily,
        "alignment_confidence": _confidence_class,
        "local_opportunity": lambda r: f"{r.get('host_strand', '?')}|{_support_bin(r)}",
    },
    "nested_vs_unnested": {
        "truncation": _truncation_class,
        "allele_frequency": _af_bin,
        "local_opportunity": _support_bin,
    },
}

#: The coarsening ladder, finest first, per comparison.
LADDERS: dict[str, dict[str, tuple[str, ...]]] = {
    "peak_vs_rest_of_host": {
        "full": ("host_subfamily", "alignment_confidence", "local_opportunity"),
        "drop_local_opportunity": ("host_subfamily", "alignment_confidence"),
        "drop_alignment_confidence": ("host_subfamily", "local_opportunity"),
        "subfamily_only": ("host_subfamily",),
        "no_strata": (),
    },
    "nested_vs_unnested": {
        "full": ("truncation", "allele_frequency", "local_opportunity"),
        "drop_local_opportunity": ("truncation", "allele_frequency"),
        "drop_allele_frequency": ("truncation", "local_opportunity"),
        "single_covariate": ("local_opportunity",),
        "no_strata": (),
    },
}


def stratum_key(row: dict[str, str], comparison: str, level: str) -> tuple[str, ...]:
    fns = COVARIATE_SETS[comparison]
    return tuple(fns[c](row) for c in LADDERS[comparison][level])


def coarsened_match(
    exposed: Sequence[dict[str, str]],
    reference: Sequence[dict[str, str]],
    comparison: str,
) -> dict[str, Any]:
    """Pick the finest ladder rung that retains enough exposed sites.

    The rung chosen and the exposed rows lost to it are reported together. A
    rung is usable when it keeps at least MIN_EXPOSED_RETAINED_FRACTION of the
    exposed rows; the last rung keeps everything but strata on nothing, which is
    the honest floor and is labelled `no_strata` when it is used.
    """
    ladder = LADDERS[comparison]
    attempts = []
    for level, covariates in ladder.items():
        index: dict[tuple[str, ...], list[dict[str, str]]] = collections.defaultdict(list)
        for row in reference:
            index[stratum_key(row, comparison, level)].append(row)
        kept = [
            r
            for r in exposed
            if len(index.get(stratum_key(r, comparison, level), [])) >= MIN_REFERENCE_PER_STRATUM
        ]
        frac = len(kept) / len(exposed) if exposed else 0.0
        attempts.append(
            {
                "level": level,
                "covariates": list(covariates),
                "exposed_kept": len(kept),
                "exposed_dropped": len(exposed) - len(kept),
                "exposed_fraction_retained": frac,
            }
        )
        if frac >= MIN_EXPOSED_RETAINED_FRACTION:
            return _pack_match(kept, exposed, reference, comparison, level, attempts)

    level = list(ladder)[-1]
    return _pack_match(list(exposed), exposed, reference, comparison, level, attempts)


def _pack_match(
    kept: list[dict[str, str]],
    exposed: Sequence[dict[str, str]],
    reference: Sequence[dict[str, str]],
    comparison: str,
    level: str,
    attempts: list[dict[str, Any]],
) -> dict[str, Any]:
    index: dict[tuple[str, ...], list[dict[str, str]]] = collections.defaultdict(list)
    for row in reference:
        index[stratum_key(row, comparison, level)].append(row)
    keys = sorted({stratum_key(r, comparison, level) for r in kept})
    return {
        "level": level,
        "covariates": list(LADDERS[comparison][level]),
        "exposed_total": len(exposed),
        "exposed_matched": len(kept),
        "exposed_dropped_unsupported": len(exposed) - len(kept),
        "strata": [
            {
                "key": list(k),
                "exposed": [r for r in kept if stratum_key(r, comparison, level) == k],
                "reference": index.get(k, []),
            }
            for k in keys
        ],
        "attempts": attempts,
    }


def balance_table(
    strata: Sequence[dict[str, Any]], comparison: str, level: str
) -> list[dict[str, Any]]:
    """Standardised mean difference per covariate, on matched strata only.

    A covariate whose levels never co-occur in a stratum is reported as `not
    estimable`, which is itself the reason that covariate can be dropped from
    the ladder.
    """
    out = []
    n_e = sum(len(st["exposed"]) for st in strata)
    n_r = sum(len(st["reference"]) for st in strata)
    for name in LADDERS[comparison][level]:
        fn = COVARIATE_SETS[comparison][name]
        per_level: dict[str, list[int]] = collections.defaultdict(lambda: [0, 0])
        for st in strata:
            for row in st["exposed"]:
                per_level[fn(row)][0] += 1
            for row in st["reference"]:
                per_level[fn(row)][1] += 1
        levels = []
        for lvl, (n_exp, n_ref) in sorted(per_level.items()):
            p_e = n_exp / n_e if n_e else 0.0
            p_r = n_ref / n_r if n_r else 0.0
            denom = ((p_e * (1 - p_e) + p_r * (1 - p_r)) / 2) ** 0.5
            smd = (p_e - p_r) / denom if denom else 0.0
            levels.append(
                {"level_value": lvl, "exposed": n_exp, "reference": n_ref, "smd": smd}
            )
        estimable = [lv for lv in levels if lv["exposed"] and lv["reference"]]
        out.append(
            {
                "covariate": name,
                "levels": levels,
                "max_abs_smd": max((abs(lv["smd"]) for lv in estimable), default=None),
                "estimable": bool(estimable),
            }
        )
    return out


# --------------------------------------------------------------------------
# exact-conditional inference
# --------------------------------------------------------------------------


def _status(event: Callable[[dict[str, str]], bool | None], row: dict[str, str]) -> np.int8:
    v = event(row)
    if v is True:
        return EVENT_TRUE
    return EVENT_FALSE if v is False else EVENT_UNEVALUABLE


def _pooled_risk_difference(tables: Sequence[tuple[int, int, int, int]]) -> float | None:
    """Mantel-Haenszel risk difference: stratum risk differences weighted by stratum size.

    Weighting by the stratum total N_i and dividing by sum(N_i) makes a
    single-stratum cohort reduce to its own crude difference. Leaving the
    weight off silently shrinks the statistic by a factor of N_i, which would
    have made every effect here look indistinguishable from zero.
    """
    num = den = 0.0
    for a, b, c, d in tables:
        n_e, n_r = a + b, c + d
        if n_e == 0 or n_r == 0:
            continue
        n_i = n_e + n_r
        num += n_i * (a / n_e - c / n_r)
        den += n_i
    return (num / den) if den else None


def _pooled_odds_ratio(tables: Sequence[tuple[int, int, int, int]]) -> float | None:
    num = den = 0.0
    for a, b, c, d in tables:
        n = a + b + c + d
        if n == 0 or 0 in (a, b, c, d):
            continue
        num += a * d / n
        den += b * c / n
    return (num / den) if den else None


def _benjamini_hochberg(pvalues: Sequence[float | None]) -> list[float | None]:
    """BH-adjusted p-values, with None carried through untouched."""
    idx = [i for i, p in enumerate(pvalues) if p is not None]
    out: list[float | None] = list(pvalues)
    if not idx:
        return out
    m = len(idx)
    ordered = sorted(idx, key=lambda i: pvalues[i])
    running = 1.0
    for rank in range(m, 0, -1):
        i = ordered[rank - 1]
        running = min(running, pvalues[i] * m / rank)
        out[i] = running
    return out


def stratified_tests(
    strata: Sequence[dict[str, Any]],
    events: dict[str, Callable[[dict[str, str]], bool | None]],
    seed: int = DEFAULT_SEED,
    permutations: int = DEFAULT_PERMUTATIONS,
) -> dict[str, dict[str, Any]]:
    """Stratified risk differences with an exact-conditional randomisation p.

    Group sizes are held fixed inside each stratum and labels are permuted, so
    the null is the sharp one and needs no asymptotic approximation. Rows whose
    event is None are unevaluable: they leave their own denominator rather than
    being counted as non-events.

    All outcomes share one permutation per replicate. Each outcome's marginal
    null is unchanged by that -- the permuted labels do not depend on the
    outcome -- but the replicate count drops by the number of outcomes and the
    p-values become positively correlated, which is what shared random numbers
    are for. Outcomes are Benjamini-Hochberg adjusted within a comparison,
    since they are eight hypotheses tested on the same rows and an unadjusted
    set of eight p-values invites the reader to cherry-pick.
    """
    labels = list(events)
    packed = []
    for st in strata:
        exposed = [[_status(events[lab], r) for r in st["exposed"]] for lab in labels]
        reference = [[_status(events[lab], r) for r in st["reference"]] for lab in labels]
        n_e = len(st["exposed"])
        status = np.array([e + r for e, r in zip(exposed, reference)], dtype=np.int8)
        counts = []
        for j in range(len(labels)):
            a = sum(1 for v in exposed[j] if v == EVENT_TRUE)
            b = sum(1 for v in exposed[j] if v == EVENT_FALSE)
            c = sum(1 for v in reference[j] if v == EVENT_TRUE)
            d = sum(1 for v in reference[j] if v == EVENT_FALSE)
            counts.append(
                {
                    "a": a, "b": b, "c": c, "d": d,
                    "unevaluable_exposed": n_e - a - b,
                    "unevaluable_reference": len(st["reference"]) - c - d,
                }
            )
        packed.append({"key": st["key"], "status": status, "n_exposed": n_e, "counts": counts})

    tables_by_outcome = [
        [(p["counts"][j]["a"], p["counts"][j]["b"], p["counts"][j]["c"], p["counts"][j]["d"])
         for p in packed]
        for j in range(len(labels))
    ]
    observed_rd = [_pooled_risk_difference(t) for t in tables_by_outcome]
    observed_or = [_pooled_odds_ratio(t) for t in tables_by_outcome]

    n_out = len(labels)
    null = np.zeros((n_out, permutations), dtype=np.float64)
    has_null = np.zeros(n_out, dtype=bool)
    rng = np.random.default_rng(seed)
    for i in range(permutations):
        tables: list[list[tuple[int, int, int, int]]] = [[] for _ in range(n_out)]
        for p in packed:
            n_e, size = p["n_exposed"], p["status"].shape[1]
            for j in range(n_out):
                cn = p["counts"][j]
                if n_e <= 0 or n_e >= size:
                    tables[j].append((cn["a"], cn["b"], cn["c"], cn["d"]))
                    continue
                chosen = np.argpartition(rng.random(size), n_e)[:n_e]
                sel = p["status"][j][chosen]
                a = int(np.count_nonzero(sel == EVENT_TRUE))
                b = int(np.count_nonzero(sel == EVENT_FALSE))
                un = int(np.count_nonzero(sel == EVENT_UNEVALUABLE))
                # The stratum's totals are fixed; only the split moves. The new
                # reference count is the stratum total minus what the exposed
                # group took, NOT the observed reference count minus that, which
                # would subtract the events the exposed group never had.
                tables[j].append(
                    (a, n_e - a - un, cn["a"] + cn["c"] - a, cn["b"] + cn["d"] - b)
                )
        for j in range(n_out):
            v = _pooled_risk_difference(tables[j])
            if v is not None:
                null[j, i] = v
                has_null[j] = True

    raw_ps: list[float | None] = []
    for j in range(n_out):
        rd = observed_rd[j]
        if rd is None or not has_null[j]:
            raw_ps.append(None)
            continue
        extreme = int(np.count_nonzero(np.abs(null[j]) >= abs(rd) - 1e-12))
        raw_ps.append((extreme + 1) / (permutations + 1))
    adjusted = _benjamini_hochberg(raw_ps)

    out: dict[str, dict[str, Any]] = {}
    for j, lab in enumerate(labels):
        tot_a = sum(p["counts"][j]["a"] for p in packed)
        tot_b = sum(p["counts"][j]["b"] for p in packed)
        tot_c = sum(p["counts"][j]["c"] for p in packed)
        tot_d = sum(p["counts"][j]["d"] for p in packed)
        exposed_rate = tot_a / (tot_a + tot_b) if tot_a + tot_b else None
        reference_rate = tot_c / (tot_c + tot_d) if tot_c + tot_d else None
        out[lab] = {
            "observed": {
                "exposed_events": tot_a,
                "exposed_nonevents": tot_b,
                "reference_events": tot_c,
                "reference_nonevents": tot_d,
                "exposed_unevaluable": sum(p["counts"][j]["unevaluable_exposed"] for p in packed),
                "reference_unevaluable": sum(p["counts"][j]["unevaluable_reference"] for p in packed),
                "exposed_rate": exposed_rate,
                "reference_rate": reference_rate,
                "crude_difference": (
                    exposed_rate - reference_rate
                    if exposed_rate is not None and reference_rate is not None
                    else None
                ),
                "risk_difference": observed_rd[j],
                "odds_ratio": observed_or[j],
                "power_note": (
                    "fewer than 10 events in the exposed arm: underpowered, and the point "
                    "estimate should not be read as an effect size"
                    if min(tot_a, tot_c) < 10
                    else None
                ),
                "note": (
                    "crude_difference pools matched rows and ignores strata; "
                    "risk_difference is the Mantel-Haenszel stratified estimate. When they "
                    "disagree in sign, the stratification is doing the work and the crude "
                    "number is the misleading one."
                ),
            },
            "permutation_p": raw_ps[j],
            "permutation_p_bh": adjusted[j],
            "permutations": permutations,
            "seed": seed,
            "n_strata_contributing": len(packed),
            "strata": [
                {
                    "key": p["key"],
                    "exposed_total": p["n_exposed"],
                    "reference_total": int(p["status"].shape[1]) - p["n_exposed"],
                    **{k: p["counts"][j][k] for k in
                       ("a", "b", "c", "d", "unevaluable_exposed", "unevaluable_reference")},
                }
                for p in packed
            ],
        }
    return out


# --------------------------------------------------------------------------
# TSD structure (the plan's Phase 2b deliverable)
# --------------------------------------------------------------------------


def tsd_verification(scored: Sequence[dict[str, Any]]) -> dict[str, Any]:
    """How well the caller's reported TSD survives contact with the reference.

    Reported because it is the measurement that decides whether the TSD arm of
    the triad means anything. Two numbers matter: how many calls match exactly,
    and *where* the mismatches sit when they do not. Mismatches bunched at the
    terminal bases are a different defect from mismatches scattered along the
    sequence, and conflating them would hide which one is happening.
    """
    n_eval = n_exact = 0
    n_all_terminal = 0
    mm_hist: collections.Counter = collections.Counter()
    offset_hist: collections.Counter = collections.Counter()
    for s in scored:
        tsd = s["components"]["tsd"]
        if tsd["state"] not in ("reference_verifiable", "reference_inconsistent"):
            continue
        n_eval += 1
        if tsd["exact"]:
            n_exact += 1
        mm_hist[tsd["n_mismatches"]] += 1
        offsets = tsd["mismatch_offsets_from_3prime_end"] or []
        for off in offsets:
            offset_hist[min(off, 6)] += 1
        if offsets and all(off <= 2 for off in offsets):
            n_all_terminal += 1
    total_offsets = sum(offset_hist.values()) or 1
    terminal = sum(v for k, v in offset_hist.items() if k <= 2)
    near_allowance = sum(v for k, v in mm_hist.items() if k == TSD_MAX_MISMATCH)
    return {
        "n_evaluable": n_eval,
        "n_exact_match": n_exact,
        "exact_match_fraction": n_exact / n_eval if n_eval else None,
        "n_at_mismatch_allowance": near_allowance,
        "share_at_mismatch_allowance": near_allowance / n_eval if n_eval else None,
        "n_more_than_allowance": sum(v for k, v in mm_hist.items() if k > TSD_MAX_MISMATCH),
        "mismatch_count_histogram": dict(sorted(mm_hist.items())),
        "mismatch_offset_from_3prime_end_histogram": dict(sorted(offset_hist.items())),
        "share_of_mismatches_in_last_3_bases": terminal / total_offsets,
        "n_calls_with_all_mismatches_in_last_3_bases": n_all_terminal,
        "share_of_calls_whose_mismatches_are_all_terminal": (
            n_all_terminal / n_eval if n_eval else None
        ),
        "max_mismatch_allowed": TSD_MAX_MISMATCH,
        "finding": (
            "A true two-junction TSD reconstruction is impossible from these inputs: "
            "the reference lacks the duplicated copy and no sample reads are available. "
            "This is a one-sided consistency check, and the numbers above say how far "
            "it gets."
        ),
    }


def _quantiles(values: Sequence[int]) -> dict[str, Any]:
    s = sorted(values)
    if not s:
        return {"n": 0}

    def q(p: float) -> int:
        return s[min(len(s) - 1, int(p * len(s)))]

    return {"n": len(s), "min": s[0], "p25": q(0.25), "median": q(0.5),
            "p75": q(0.75), "p90": q(0.9), "p99": q(0.99), "max": s[-1]}


def tsd_structure(rows: Sequence[dict[str, str]], sentinels: frozenset[int]) -> dict[str, Any]:
    """Measure the sentinel-free TSD distribution and stratify drift by TSD class.

    Gate 0d's consequence is that TSD length is not a continuous variable and
    the sentinels may not be used. This reports the distribution twice -- once
    with them, once without -- because the difference is the point. Drift is
    stratified rather than assumed away: a TSD class with systematically larger
    consensus-vs-host offset drift is *less* trustworthy, not more enriched.
    """
    lens_all = [
        n for n in (_to_int(r.get("tsd_len")) for r in rows) if n is not None
    ]
    lens_clean = [n for n in lens_all if n not in sentinels]

    drift_by_class: dict[str, list[int]] = collections.defaultdict(list)
    for row in rows:
        drift = _to_int(row.get("offset_drift_bp"))
        if drift is not None:
            drift_by_class[tsd_length_class(_to_int(row.get("tsd_len")), sentinels)].append(
                abs(drift)
            )

    return {
        "sentinels": sorted(sentinels),
        "with_sentinels": _quantiles(lens_all),
        "sentinel_excluded": _quantiles(lens_clean),
        "n_sentinel_rows_excluded": len(lens_all) - len(lens_clean),
        "share_sentinel_excluded": (
            (len(lens_all) - len(lens_clean)) / len(lens_all) if lens_all else None
        ),
        "absolute_offset_drift_by_tsd_class": {
            cls: _quantiles(v) for cls, v in sorted(drift_by_class.items())
        },
    }


def pure_at_tsd_fraction(rows: Sequence[dict[str, str]], sentinels: frozenset[int]) -> dict[str, Any]:
    n_at = n_eval = 0
    for row in rows:
        n = _to_int(row.get("tsd_len"))
        seq = (row.get("tsd_seq") or "").strip().upper()
        if n is None or not seq or n in sentinels or len(seq) != n:
            continue
        n_eval += 1
        if set(seq) <= {"A", "T"}:
            n_at += 1
    return {
        "n_evaluable": n_eval,
        "n_pure_at": n_at,
        "fraction": n_at / n_eval if n_eval else None,
        "interpretation_limit": (
            "Only interpretable against a sequence- and callability-matched control; "
            "a pure-A/T TSD is partly a statement about the local sequence."
        ),
    }


# --------------------------------------------------------------------------
# scoring pipeline
# --------------------------------------------------------------------------


def score_rows(
    rows: Sequence[dict[str, str]],
    reference: Reference | None,
    motif: str,
    sentinels: frozenset[int],
) -> list[dict[str, Any]]:
    """Score every row and return it with its hallmark state attached.

    Unscorable components are recorded with their reason, so an `unevaluable`
    row can always be traced back to *why* it was unevaluable.
    """
    scored = []
    for row in rows:
        tsd = tsd_reference_check(row, reference, sentinels)
        components = {
            "tsd": tsd,
            "polya": polya_state(row),
            "en": en_motif_state(row, reference, motif),
            "l1_en": l1_en_nick_state(row, reference),
        }
        scored.append(
            {
                "row": row,
                "components": components,
                "triad": score_triad(components),
                "tsd_length_class": tsd_length_class(tsd["tsd_len"], sentinels),
                "pure_at_tsd": bool(tsd["tsd_seq"] and set(tsd["tsd_seq"].upper()) <= {"A", "T"}),
            }
        )
    return scored


def event_tsd_verifiable(row: dict[str, str]) -> bool | None:
    return {"reference_verifiable": True, "reference_inconsistent": False}.get(
        row.get("_tsd_state"), None
    )


def event_tsd_exact(row: dict[str, str]) -> bool | None:
    state = row.get("_tsd_state")
    if state == "reference_verifiable":
        return row.get("_tsd_exact") == "1"
    if state == "reference_inconsistent":
        return False
    return None


def event_tsd_de_novo_range(row: dict[str, str]) -> bool | None:
    cls = row.get("_tsd_length_class")
    if cls == "de_novo_range_5_27":
        return True
    if cls == "longer_than_de_novo_28_plus":
        return False
    return None


def event_polya_present(row: dict[str, str]) -> bool | None:
    return {"present": True, "absent": False}.get(row.get("_polya_state"), None)


def event_polya_long(row: dict[str, str]) -> bool | None:
    n = _to_int(row.get("polya_len"))
    if row.get("_polya_state") == "unevaluable" or n is None:
        return None
    return n >= POLYA_SECONDARY_MIN_LEN


def event_en_compatible(row: dict[str, str]) -> bool | None:
    return {"motif_compatible": True, "motif_incompatible": False, "motif_absent": False}.get(
        row.get("_en_state"), None
    )


def event_en_any_side(row: dict[str, str]) -> bool | None:
    return {"motif_compatible": True, "motif_incompatible": True, "motif_absent": False}.get(
        row.get("_en_state"), None
    )


def event_l1_en(row: dict[str, str]) -> bool | None:
    return {
        "tttt_a_both_sides": True,
        "tttt_a_five_prime_only": True,
        "tttt_a_three_prime_only": True,
        "no_tttt_a_at_nick": False,
    }.get(row.get("_l1_en_state"), None)


def event_triad_observed(row: dict[str, str]) -> bool | None:
    if row.get("_triad_all_observed") == "1":
        return True
    return False if row.get("_triad_complete") == "1" else None


EVENTS: dict[str, Callable[[dict[str, str]], bool | None]] = {
    "tsd_reference_verifiable": event_tsd_verifiable,
    "tsd_exact_match": event_tsd_exact,
    "tsd_in_de_novo_range_5_27": event_tsd_de_novo_range,
    "polya_present": event_polya_present,
    "polya_len_ge_20": event_polya_long,
    "en_motif_compatible": event_en_compatible,
    "en_motif_any_side": event_en_any_side,
    "l1_en_tttt_a_either_side": event_l1_en,
    "triad_all_three_observed": event_triad_observed,
}


def run_comparison(
    name: str,
    exposed: Sequence[dict[str, str]],
    reference: Sequence[dict[str, str]],
    seed: int,
    permutations: int,
) -> dict[str, Any]:
    match = coarsened_match(exposed, reference, name)
    out: dict[str, Any] = {
        "comparison": name,
        "covariates_available": sorted(COVARIATE_SETS[name]),
        "matching_rule": (
            "coarsened exact matching on a fixed ladder; no hallmark, TSD length or "
            "motif content is used as a matching variable"
        ),
        "exposed_input": len(exposed),
        "reference_input": len(reference),
        "level_used": match["level"],
        "covariates_used": match["covariates"],
        "exposed_matched": match["exposed_matched"],
        "exposed_dropped_unsupported": match["exposed_dropped_unsupported"],
        "matched": match["level"] != "no_strata",
        "interpretation": (
            "coarsened-matched estimate"
            if match["level"] != "no_strata"
            else "UNMATCHED: no ladder rung retained a usable share of the exposed arm, "
            "so this is a crude difference reported for transparency, not a matched "
            "estimate, and must not be read as one"
        ),
        "ladder_attempts": match["attempts"],
        "n_strata": len(match["strata"]),
        "balance": balance_table(match["strata"], name, match["level"]),
        "outcomes": {},
    }
    for label, res in stratified_tests(
        match["strata"], EVENTS, seed=seed, permutations=permutations
    ).items():
        out["outcomes"][label] = res
    return out


def to_rows(scored: Sequence[dict[str, Any]]) -> list[dict[str, str]]:
    """Flatten scored rows back into cohort-shaped dicts the matchers can read."""
    out = []
    for s in scored:
        row = dict(s["row"])
        row["_tsd_state"] = s["components"]["tsd"]["state"]
        row["_tsd_exact"] = "1" if s["components"]["tsd"]["exact"] else "0"
        row["_tsd_n_mismatches"] = s["components"]["tsd"]["n_mismatches"]
        row["_tsd_length_class"] = s["tsd_length_class"]
        row["_polya_state"] = s["components"]["polya"]["state"]
        row["_en_state"] = s["components"]["en"]["state"]
        row["_l1_en_state"] = s["components"]["l1_en"]["state"]
        row["_triad_complete"] = "1" if s["triad"]["triad_complete"] else "0"
        row["_triad_all_observed"] = "1" if s["triad"]["triad_complete_observed"] else "0"
        row["_in_peak_window"] = (
            "1"
            if PEAK_WINDOW[0] <= (_to_int(row.get("consensus_offset")) or -1) < PEAK_WINDOW[1]
            else "0"
        )
        out.append(row)
    return out


def component_summary(scored: Sequence[dict[str, Any]]) -> dict[str, Any]:
    rows = [s["row"] for s in scored]
    return {
        "n_rows": len(rows),
        "tsd_states": dict(collections.Counter(s["components"]["tsd"]["state"] for s in scored)),
        "tsd_reasons": dict(
            collections.Counter(
                s["components"]["tsd"]["reason"] for s in scored if s["components"]["tsd"]["reason"]
            )
        ),
        "polya_states": dict(collections.Counter(s["components"]["polya"]["state"] for s in scored)),
        "en_states": dict(collections.Counter(s["components"]["en"]["state"] for s in scored)),
        "l1_en_states": dict(collections.Counter(s["components"]["l1_en"]["state"] for s in scored)),
        "en_reasons": dict(
            collections.Counter(
                s["components"]["en"]["reason"] for s in scored if s["components"]["en"]["reason"]
            )
        ),
        "tsd_length_classes": dict(collections.Counter(s["tsd_length_class"] for s in scored)),
        "triad_complete": sum(1 for s in scored if s["triad"]["triad_complete"]),
        "triad_complete_observed": sum(1 for s in scored if s["triad"]["triad_complete_observed"]),
    }


# --------------------------------------------------------------------------
# reporting
# --------------------------------------------------------------------------


def _fmt(value: float | None, digits: int = 4) -> str:
    return "n/a" if value is None else f"{value:.{digits}f}"


def _pct(value: float | None) -> str:
    return "n/a" if value is None else f"{value:.1%}"


def _pct(value: float | None) -> str:
    return "n/a" if value is None else f"{value:.1%}"


def write_report(path: Path, report: dict[str, Any]) -> None:
    en = report["en_motif"]
    tsd = report["tsd_structure"]
    L: list[str] = [
        "# Phase 2b: TPRT hallmarks in matched comparisons",
        "",
        "Halmarks are scored per call and compared only after coarsened matching.",
        "`unevaluable` is never counted as `absent`; it leaves the denominator.",
        "",
        "## Set sizes",
        "",
    ]
    for k, v in report["set_sizes"].items():
        L.append(f"- {k}: {v}")
    L += [
        "",
        "## EN motif",
        "",
        f"- Status: **{en['status']}**",
        f"- Sequence `{en['motif']}` ({en['length_bp']} bp, <= {en['max_mismatch']} mismatches, "
        f"flank {en['flank_bp']} bp, min per-column support "
        f"{_fmt(en['per_column_support_min'], 3)}, "
        f"derived from {en['n_alu_consensus_entries']} Alu consensus entries)",
        "",
        "## TSD structure",
        "",
        f"- Sentinels excluded: {tsd['sentinels']} — {tsd['n_sentinel_rows_excluded']} rows "
        f"({tsd['share_sentinel_excluded']:.1%} of calls carrying a TSD length).",
        f"- With sentinels: median {tsd['with_sentinels']['median']}, p90 "
        f"{tsd['with_sentinels']['p90']}, max {tsd['with_sentinels']['max']}, n "
        f"{tsd['with_sentinels']['n']}.",
        f"- Sentinel-excluded: median {tsd['sentinel_excluded']['median']}, p90 "
        f"{tsd['sentinel_excluded']['p90']}, max {tsd['sentinel_excluded']['max']}, n "
        f"{tsd['sentinel_excluded']['n']}.",
        "",
        "Absolute consensus-vs-host offset drift by TSD class. A class with larger drift",
        "is less positionally trustworthy, not more enriched:",
        "",
        "| TSD class | n | median drift | p90 |",
        "|---|---|---|---|",
    ]
    for cls, prof in tsd["absolute_offset_drift_by_tsd_class"].items():
        if prof.get("n"):
            L.append(f"| {cls} | {prof['n']} | {prof['median']} | {prof['p90']} |")

    tv = report["tsd_verification"]["analysis_set"]
    L += [
        "",
        "## TSD verification (this is a finding, not a score)",
        "",
        f"- Evaluable calls: {tv['n_evaluable']}. Exact reference match: "
        f"{tv['n_exact_match']} ({_pct(tv['exact_match_fraction'])}).",
        f"- Exactly at the {tv['max_mismatch_allowed']}-mismatch allowance: "
        f"{tv['n_at_mismatch_allowance']} ({_pct(tv['share_at_mismatch_allowance'])}); "
        f"beyond it: {tv['n_more_than_allowance']}.",
        f"- All of a call's mismatches lie in the last 3 bases of the reported TSD for "
        f"{tv['n_calls_with_all_mismatches_in_last_3_bases']} calls "
        f"({_pct(tv['share_of_calls_whose_mismatches_are_all_terminal'])}); across all "
        f"individual mismatches the share in the last 3 bases is "
        f"{_pct(tv['share_of_mismatches_in_last_3_bases'])}.",
        f"- {tv['finding']}",
        "",
        "Read those numbers before drawing a conclusion from them. The spike at exactly",
        "the allowance, together with the terminal clustering among those calls, is",
        "consistent with a reported TSD whose last bases ran into neighbouring sequence;",
        "the calls beyond the allowance disagree with the reference along their whole",
        "length and are a different failure. The full histogram is in the JSON. Either",
        "way the TSD arm measures whether the caller's field agrees with the reference,",
        "not whether a target site was duplicated, so it is **not** read as evidence for",
        "or against TPRT.",
    ]

    for name in ("peak_vs_rest_of_host", "nested_vs_unnested"):
        c = report["comparisons"][name]
        L += [
            "",
            f"## Comparison: {name}",
            "",
            f"- Exposed offered {c['exposed_input']}, matched {c['exposed_matched']}, dropped "
            f"for want of reference support {c['exposed_dropped_unsupported']}.",
            f"- Ladder rung used: **{c['level_used']}** on covariates "
            f"{c['covariates_used'] or ['none — no strata']}, {c['n_strata']} strata.",
            f"- Status: **{c['interpretation']}**",
            "",
            "| outcome | exposed | reference | crude diff | MH risk diff | OR | exact p (BH) | unevaluable E/R |",
            "|---|---|---|---|---|---|---|---|",
        ]
        for label, o in c["outcomes"].items():
            obs = o["observed"]
            pcell = f"{_fmt(o['permutation_p'])} (BH {_fmt(o['permutation_p_bh'])})"
            L.append(
                f"| {label} | {_fmt(obs['exposed_rate'])} | {_fmt(obs['reference_rate'])} | "
                f"{_fmt(obs['crude_difference'])} | "
                f"{_fmt(obs['risk_difference'])} | {_fmt(obs['odds_ratio'])} | {pcell} | "
                f"{obs['exposed_unevaluable']}/{obs['reference_unevaluable']} |"
            )
            if obs.get("power_note"):
                L.append(f"| _{label}_ | | | | | | | {obs['power_note']} |")
        L += [
            "",
            "Covariate balance on the matched strata. `not estimable` means the covariate",
            "never varied within a stratum, which is the reason it can be dropped from the",
            "ladder rather than a missing measurement:",
            "",
            "| covariate | max abs SMD |",
            "|---|---|",
        ]
        for b in c["balance"]:
            v = b["max_abs_smd"]
            L.append(f"| {b['covariate']} | {'not estimable' if v is None else f'{v:.3f}'} |")

    L += [
        "",
        "## Read this before reading the EN-motif rows",
        "",
        "The Alu-5'-window motif is Alu sequence. A site inside an Alu sits inside",
        "sequence that matches it more often, for composition alone. The",
        "`nested_vs_unnested` comparison is confounded by construction and its",
        "`en_motif_*` rows must not be read as mechanism. `peak_vs_rest_of_host` does",
        "not share that problem: both arms are inside the same hosts. The",
        "`l1_en_tttt_a_either_side` rows use a probe that is not Alu-derived.",
    ]
    pa = report["pure_at_tsd"]["analysis_set"]
    L += [
        "",
        "## Pure-A/T TSD fraction",
        "",
        f"- Analysis set: {pa['n_pure_at']} of {pa['n_evaluable']} evaluable "
        f"({pa['fraction']:.1%}). Peak {report['pure_at_tsd']['peak']['fraction']:.1%}, "
        f"mid-host control {report['pure_at_tsd']['mid_host_control']['fraction']:.1%}, "
        f"unnested {report['pure_at_tsd']['unnested']['fraction']:.1%}.",
        f"- Limit: {pa['interpretation_limit']}",
        "",
        "## Deviations",
        "",
    ]
    for d in report["deviations"]:
        L.append(f"- **{d['item']}** ({d['status']}): {d['detail']}")
    L += ["", "## Disclosures", ""]
    for d in report["disclosures"]:
        L.append(f"- {d}")
    path.write_text("\n".join(L) + "\n", encoding="utf-8")


# --------------------------------------------------------------------------
# entry point
# --------------------------------------------------------------------------


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cohort", type=Path, default=DEFAULT_COHORT)
    p.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    p.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE)
    p.add_argument("--alu-consensus", type=Path, default=DEFAULT_ALU_CONSENSUS)
    p.add_argument(
        "--en-motif-seq",
        type=str,
        default=None,
        help="Override the consensus-derived EN motif proxy with a curated sequence.",
    )
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--permutations", type=int, default=DEFAULT_PERMUTATIONS)
    p.add_argument("--seed", type=int, default=DEFAULT_SEED)
    p.add_argument(
        "--no-reference",
        action="store_true",
        help="Skip reference-backed components (TSD junctions, EN motif). Everything "
        "reference-dependent becomes `unevaluable` rather than `absent`.",
    )
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    args.outdir.mkdir(parents=True, exist_ok=True)

    sentinels_list, sentinel_source = load_sentinels(args.manifest)
    sentinels = frozenset(sentinels_list)
    print(f"TSD sentinels ({sentinel_source}): {sorted(sentinels)}", flush=True)

    if args.en_motif_seq:
        motif = args.en_motif_seq.strip().upper()
        motif_meta = {
            "motif": motif,
            "length_bp": len(motif),
            "max_mismatch": EN_MOTIF_MAX_MISMATCH,
            "flank_bp": EN_MOTIF_FLANK_BP,
            "per_column_support_min": None,
            "per_column_support_median": None,
            "n_alu_consensus_entries": None,
            "status": "user_supplied",
        }
    else:
        motif, motif_meta = derive_en_motif_proxy(args.alu_consensus)

    reference = None if args.no_reference else Reference(args.reference)
    if reference is None:
        print("no reference: TSD-junction and EN-motif components become unevaluable", flush=True)

    rows = read_cohort(args.cohort)
    analysis = nested_analysis_set(rows)
    peak_rows, mid_rows = split_by_offset(analysis)
    # An empty analysis set is a broken input, not a null result. Scoring an
    # empty list and running both comparisons over empty arms produces a
    # well-formed report full of zeros and "no difference detected", which reads
    # as a decisive negative. The usual cause is a schema break: a gate column
    # missing from the cohort table makes `.get()` return None, every row fails
    # the gate, and the cohort empties silently. `test_position_enrichment.py`
    # and both Phase 4 scripts already refuse this case; this was the one
    # consumer that did not.
    if not analysis:
        raise SystemExit(
            f"{args.cohort} yielded an empty analysis set "
            f"({len(rows)} rows read; gate is nested_in_alu_host == 1, "
            f"{NEAR_FULL_HOST[0]}-{NEAR_FULL_HOST[1]} bp host, unambiguous "
            f"consensus mapping, numeric consensus_offset). Refusing to report "
            f"an empty cohort as a null result."
        )
    unnested_rows = [r for r in rows if r.get("nested_in_alu_host") != "1"]
    print(
        f"analysis set {len(analysis)} | peak {len(peak_rows)} | mid-host control "
        f"{len(mid_rows)} | unnested arm {len(unnested_rows)}",
        flush=True,
    )

    scored_analysis = score_rows(analysis, reference, motif, sentinels)
    scored_peak = score_rows(peak_rows, reference, motif, sentinels)
    scored_mid = score_rows(mid_rows, reference, motif, sentinels)
    scored_unnested = score_rows(unnested_rows, reference, motif, sentinels)

    peak = to_rows(scored_peak)
    mid = to_rows(scored_mid)
    nested_all = to_rows(scored_analysis)
    unnested = to_rows(scored_unnested)
    rest_of_host = [r for r in nested_all if r.get("_in_peak_window") != "1"]

    comparisons = {
        "peak_vs_rest_of_host": run_comparison(
            "peak_vs_rest_of_host", peak, rest_of_host, args.seed, args.permutations
        ),
        "nested_vs_unnested": run_comparison(
            "nested_vs_unnested", nested_all, unnested, args.seed, args.permutations
        ),
    }

    report = {
        "phase": "2b",
        "description": (
            "TPRT hallmark triad scored per call, compared only under coarsened matching "
            "on covariates that are not themselves hallmarks."
        ),
        "inputs": {
            "cohort": str(args.cohort),
            "manifest": str(args.manifest),
            "reference": None if args.no_reference else str(args.reference),
            "alu_consensus": None if args.en_motif_seq else str(args.alu_consensus),
            "permutations": args.permutations,
            "seed": args.seed,
        },
        "tsd_sentinels": sorted(sentinels),
        "tsd_sentinel_source": sentinel_source,
        "en_motif": motif_meta,
        "de_novo_tsd_range": list(DE_NOVO_TSD_RANGE),
        "polya_min_len": POLYA_MIN_LEN,
        "peak_window": list(PEAK_WINDOW),
        "mid_host_control": list(MID_HOST_CONTROL),
        "set_sizes": {
            "cohort_rows": len(rows),
            "analysis_set": len(nested_all),
            "peak": len(peak),
            "mid_host_control": len(mid),
            "unnested": len(unnested),
        },
        "component_summary": {
            "analysis_set": component_summary(scored_analysis),
            "peak": component_summary(scored_peak),
            "mid_host_control": component_summary(scored_mid),
            "unnested": component_summary(scored_unnested),
        },
        "tsd_structure": tsd_structure(nested_all, sentinels),
        "tsd_verification": {
            "analysis_set": tsd_verification(scored_analysis),
            "peak": tsd_verification(scored_peak),
            "mid_host_control": tsd_verification(scored_mid),
            "unnested": tsd_verification(scored_unnested),
        },
        "pure_at_tsd": {
            "analysis_set": pure_at_tsd_fraction(nested_all, sentinels),
            "peak": pure_at_tsd_fraction(peak, sentinels),
            "mid_host_control": pure_at_tsd_fraction(mid, sentinels),
            "unnested": pure_at_tsd_fraction(unnested, sentinels),
        },
        "comparisons": comparisons,
        "deviations": [
            {
                "item": "inserted-element subfamily as a matching covariate",
                "status": "unavailable",
                "detail": (
                    "The nested-vs-unnested comparison is meant to match on inserted "
                    "subfamily. The cohort records only the *host* subfamily "
                    "(host_name), and only for nested rows; non-nested rows carry no "
                    "host columns at all, so host subfamily cannot stand in for it. "
                    "The covariate is therefore absent from COVARIATE_SETS and the "
                    "comparison is weaker for it, rather than silently substituted."
                ),
            },
            {
                "item": "curated EN consensus",
                "status": "not_verified",
                "detail": (
                    "No primary source for a curated Alu EN consensus could be read in "
                    "this environment, so the motif is derived from the Alu repeat "
                    "consensus 5' window and labelled a proxy. Supply --en-motif-seq to "
                    "replace it."
                ),
            },
        ],
        "disclosures": [
            "The Alu-5'-window EN motif is a consensus-derived proxy whose sequence is "
            "Alu sequence. In nested_vs_unnested the exposed arm lies inside Alus by "
            "construction, so that comparison's en_motif_* rows are confounded by local "
            "composition and are not evidence about mechanism. peak_vs_rest_of_host does "
            "not share this confound: both arms lie inside the same hosts.",
            "The L1 5'-TTTT/A-3' probe is reported on both sides of the nick without "
            "asserting which side the inserted element's 5' end maps to, because that "
            "depends on the caller's strand convention.",
            "A sequence-compatibility probe cannot show that L1 ORF2p nicked at that base.",
            "TSD reconstruction compares the caller-reported TSD sequence against both "
            "reference junctions. It is a consistency check on the reported call, not a "
            "read-level assembly of the insertion.",
            "poly(A) length is a caller-reported field with its own resolution floor; "
            "rows without the field are `unevaluable`, not `absent`.",
            "An EN hit that straddles the nick is assigned by its start position, so the "
            "left and right windows overlap by len(motif)-1 bases.",
            "Hallmark content separating two matched arms is not evidence that the "
            "hallmark caused the difference; both arms come from one caller.",
            "Every arm shares one caller's conventions. Separation is a statement about "
            "calls, not about biology.",
            "No precision, recall or detection-accuracy claim is made or licensed.",
        ],
    }

    (args.outdir / "tprt_hallmarks.json").write_text(
        json.dumps(report, indent=2, default=str) + "\n", encoding="utf-8"
    )
    write_report(args.outdir / "tprt_hallmarks.md", report)

    for name, c in comparisons.items():
        print(f"\n{name}: rung={c['level_used']} matched {c['exposed_matched']}/{c['exposed_input']}")
        for label, o in c["outcomes"].items():
            obs = o["observed"]
            print(
                f"  {label}: {obs['exposed_rate']} vs {obs['reference_rate']} "
                f"crude={obs['crude_difference']} MH_RD={obs['risk_difference']} "
                f"p={o['permutation_p']} (BH {o['permutation_p_bh']})"
            )
    print(f"\nwrote {args.outdir / 'tprt_hallmarks.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
