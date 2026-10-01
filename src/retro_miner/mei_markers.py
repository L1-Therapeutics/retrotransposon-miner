"""Structural markers for MEI insertion calls, with literature-anchored constants.

Every constant below is a structural signature from the primary literature —
measurable, testable features of an insertion event.  None of them is a
clinical biomarker; the clinical layer (plasma ORF1p, cfDNA hypomethylation)
is deliberately out of scope here because those are population statistics over
assays, not per-call features a structural pipeline can compute.

Provenance discipline mirrors the rest of the repo: no helper is claimed to
"validate" anything.  Marker functions *score* evidence; callers must still
reconstruct breakpoints from reads before any per-call claim is published.

Tiers
-----
Tier 1  features already computed by the caller pipeline (prioritization and
        null-model use), kept here with exact literature numbers for reference.
Tier 2  genuinely new per-call markers this module adds: 5' microhomology
        route tags, internal-priming (false-polyA) detection, complex-event
        flags, and null-side context features.

References
----------
Flasch et al. 2019    Cell 176:928.        doi:10.1016/j.cell.2019.02.002
Wagstaff et al. 2012  PLoS Genet 8:e1002842.
Zingler et al. 2005   Genome Res 15:780.   doi:10.1101/gr.3421505
Srikanta et al. 2009  J Mol Biol 385:814.  doi:10.1016/j.jmb.2008.10.094
An et al. 2006        Nucleic Acids Res 34:1699 (internal priming).
Gilbert et al. 2005   Nat Genet 37:1166 (inversion/deletion structures).
Goodier et al. 2000 / Pickeral et al. 2000 / Szak et al. 2003
                      (L1 3' transduction rate surveys).
Damert et al. 2009    Genome Res 19:1992 (SVA 5' transduction).
Cost et al. 2002      Nat Genet 32:114 (L1 microhomology context).
Jurka 1997            PNAS 94:1872 (original 5'-TTTT/AA-3' consensus).
Levin 2025            Genetics 231:iyaf202 (truncation review).
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from enum import Enum
from functools import lru_cache
from importlib import resources

__all__ = [
    "L1_EN_CONSENSUS_MOTIF",
    "EN_NICK_OFFSET",
    "EN_MOTIF_7MERS_TOTAL",
    "EN_MOTIF_7MERS_97PCT",
    "EN_MOTIF_RANK_TTTTTAA",
    "EN_MOTIF_CANONICAL_PREF_PCT",
    "EN_MOTIF_C2345_FRACTION",
    "TSD_MIN_BP",
    "TSD_MAX_BP",
    "TSD_MEAN_BP",
    "TRIAD_FRACTION",
    "MH_MAX_LEN",
    "MH_TRUNCATED_L1_RATE",
    "MH_TRUNCATED_L1_COUNTS",
    "MH_FULL_LENGTH_L1_RATE",
    "MH_RATE_DESCRIPTION",
    "MH_ALU_NOT_MICROHOMOLOGY_DEPENDENT",
    "INTERNAL_PRIMING_MIN_TSD",
    "L1_TRANSDUCTION_RATE_RANGE",
    "SVA_TRANSDUCTION_RATE",
    "SVA_MAST2_GROUP_FRACTION_SVA_F",
    "ALU_COMPLEX_EVENT_RATE",
    "L1_TRUNCATED_INVERSION_RATE",
    "NULL_GC_WINDOW_50BP",
    "ComplexEventKind",
    "EN_EMPIRICAL_CONSENSUS_MOTIF",
    "EN_MOTIF_WEIGHTS",
    "en_motif_site_rank",
    "en_motif_site_weight",
    "en_motif_scan_table",
    "en_motif_scan_table_corrected",
    "flasch_weights_provenance",
    "ComplexEventFlags",
    "InternalPrimingSignature",
    "RouteEvidence",
    "RouteClassification",
    "classify_route",
    "scan_en_motif_strand",
    "en_motif_opportunity",
    "tsd_len_in_range",
    "microhomology_at_5p_junction",
    "internal_priming_signature",
    "flag_complex_event",
    "at_richness_footprint",
    "AluHotspotPosition",
    "ALU_HOTSPOT_POSITIONS",
    "levy_hotspot_likelihood",
    "revcomp",
]

# ---------------------------------------------------------------------------
# Tier 1 — EN target motif (Flasch et al. 2019; Jurka 1997; Levy et al. 2010)
# ---------------------------------------------------------------------------
#: The original 1997 consensus of Jurka, still quoted by Flasch et al. 2019
#: when citing Feng 1996: 5'-TTTT/AA-3', "/" marking the nick.
L1_EN_CONSENSUS_MOTIF = "TTTT/AA"

#: The empirical consensus from Flasch et al. 2019 logo plots (Fig. 2C):
#: 5'-TTTTT/AA-3' — five T's before the nick, not four.  This is the correct
#: nick offset for structural work on de novo events and for the 7-mer
#: weight table, whose sites all span the nick between offsets 5 and 6.
EN_EMPIRICAL_CONSENSUS_MOTIF = "TTTTT/AA"

#: Offset of the nick within the 7-bp empirical motif window (0-based):
#: bases 0..4 lie 5' of the nick ("TTTTT"), bases 5..6 are the "AA" side.
EN_NICK_OFFSET = 5

_EN_MOTIF_6MER = "TTTTAA"
_EN_MOTIF_WINDOW = len(_EN_MOTIF_6MER)

#: Flasch et al. 2019: 6% of all possible 7-mers (743 of 12,288) account for
#: 97% of de novo L1 insertions; the canonical TTTTT/AA 7-mer itself is only
#: 9.6% by normalized preference (rank 21), and 45% of observed 7-mers carry
#: a C in positions 2-5.
EN_MOTIF_7MERS_TOTAL = 12_288
EN_MOTIF_7MERS_97PCT = 743
EN_MOTIF_RANK_TTTTTAA = 21
EN_MOTIF_CANONICAL_PREF_PCT = 9.6
EN_MOTIF_C2345_FRACTION = 0.45

#: Published top-set weight tables (Supplemental Dataset 4).  Weights are
#: relative to the table's top site; the corrected model normalizes for
#: genomic site frequency, which is why canonical TTTTT/AA ranks 21st there.
#: Only the top 21 (corrected) / 20 (uncorrected) of the 743 used sites are
#: published — see the JSON's description field for what absence means.
_EN_WEIGHTS_RESOURCE = "data/flasch2019_7mer_weights.json"


@lru_cache(maxsize=1)
def _load_flasch_weights() -> dict[str, dict[str, float]]:
    """Load the embedded Flasch 2019 weight tables (cached)."""
    text = resources.files("retro_miner").joinpath(_EN_WEIGHTS_RESOURCE).read_text(
        encoding="utf-8"
    )
    payload = json.loads(text)
    return {
        name: dict(table)
        for name, table in payload.items()
        if name in {"corrected", "uncorrected"} and isinstance(table, dict)
    }


#: Module-level handle on the embedded tables (loaded once at import).
EN_MOTIF_WEIGHTS: dict[str, dict[str, float]] = _load_flasch_weights()


def flasch_weights_provenance() -> dict[str, object]:
    """Provenance block for the embedded weight tables (for reports/exports)."""
    text = resources.files("retro_miner").joinpath(_EN_WEIGHTS_RESOURCE).read_text(
        encoding="utf-8"
    )
    payload = json.loads(text)
    return {
        key: payload[key]
        for key in (
            "source_paper",
            "source_file",
            "url",
            "pmcid",
            "sha256",
            "retrieved",
            "description",
        )
    }


def en_motif_site_weight(site: str, *, corrected: bool = True) -> float | None:
    """Flasch model weight for a 7-mer site, or None if outside the published
    top set.  The corrected table is the default: it is the one that matches
    the paper's normalized-preference statement (TTTTT/AA rank 21, 0.2186).
    Absence is "not in the published top 21", not zero.
    """
    table = EN_MOTIF_WEIGHTS["corrected" if corrected else "uncorrected"]
    return table.get(site.upper())


def en_motif_site_rank(site: str, *, corrected: bool = True) -> int | None:
    """1-based rank of a 7-mer site within the published top set, or None."""
    table = EN_MOTIF_WEIGHTS["corrected" if corrected else "uncorrected"]
    ordered = sorted(table.items(), key=lambda kv: -kv[1])
    for idx, (candidate, _weight) in enumerate(ordered, start=1):
        if candidate == site.upper():
            return idx
    return None


def en_motif_scan_table(*, corrected: bool = True) -> tuple[str, ...]:
    """7-mer sites for the calibration tier of :func:`en_motif_score`.

    The tier system exists because the tables cover only the published top
    set: *table* = these 21 sites, weight-scaled; *fallback* = canonical
    TTTTAA hexamer (always computable); *none* otherwise.  A tier-aware
    caller must not compare a table-tier score against a fallback-tier score
    as if both were on one scale.
    """
    table = EN_MOTIF_WEIGHTS["corrected" if corrected else "uncorrected"]
    return tuple(sorted(table))


def en_motif_scan_table_corrected() -> tuple[str, ...]:
    """Corrected-model sites (see :func:`en_motif_scan_table`)."""
    return en_motif_scan_table(corrected=True)

# ---------------------------------------------------------------------------
# Tier 1 — TPRT triad (Wagstaff et al. 2012, de novo tagged Alu in HeLa)
# ---------------------------------------------------------------------------
#: Inclusive range of target-site-duplication lengths at fully characterized
#: de novo Alu events ("TSD 5-27 bp, mean 14").
TSD_MIN_BP = 5
TSD_MAX_BP = 27
TSD_MEAN_BP = 14

#: Fraction of fully characterized de novo Alu carrying the *joint* triad:
#: target-site duplication + 3' poly(A) tail + L1 EN-like target site.
TRIAD_FRACTION = 0.96

# ---------------------------------------------------------------------------
# Tier 2a — 5' junction microhomology (Zingler et al. 2005; Cost et al. 2002)
# ---------------------------------------------------------------------------
#: 1-12 nt overlap between target-flank end and inserted-element start.
MH_MAX_LEN = 12

#: 5'-truncated L1 insertions with 5'-junction microhomology: 1,503 of 2,427
#: (62%) — a significant excess over chance.  Route tag, not mechanism proof.
MH_TRUNCATED_L1_RATE = 0.62
MH_TRUNCATED_L1_COUNTS = (1_503, 2_427)

#: Full-length L1 5' junctions show microhomology at ~24%, i.e. ≈ chance —
#: the discriminator between truncated and full-length L1 structures.
MH_FULL_LENGTH_L1_RATE = 0.24

#: One-line route-tag interpretation for reports.  Deliberately NOT a
#: twin-priming fingerprint claim.
MH_RATE_DESCRIPTION = (
    "microhomology at a 5' junction says 'L1-type truncated/alternative route', "
    "its absence with intact structure says 'clean Alu-type TPRT'"
)

#: Zingler et al. 2005 negative result kept for reports: 10,062 standard Alus
#: showed only a minor microhomology shift — Alu integration is not
#: microhomology-dependent.
MH_ALU_NOT_MICROHOMOLOGY_DEPENDENT = True

# ---------------------------------------------------------------------------
# Tier 2b — Internal priming / false-polyA detection (Srikanta et al. 2009)
# ---------------------------------------------------------------------------
#: Srikanta-classified internal-priming events retain TSDs of >= 6 bp, so a
#: TSD shorter than 6 bp at an A-rich insertion is ambiguous with tail polyA.
INTERNAL_PRIMING_MIN_TSD = 6

#: Srikanta et al. 2009 verified 20 events (6 human-specific: 2 Alu, 4 L1;
#: 14 hominin-shared), validated against outgroup-primate empty sites.
INTERNAL_PRIMING_VERIFIED_EVENTS = 20
INTERNAL_PRIMING_HUMAN_SPECIFIC = (2, 4)  # (Alu, L1)

# ---------------------------------------------------------------------------
# Tier 2c — Transduced flanking sequence (source tags)
# ---------------------------------------------------------------------------
#: L1 3' transduction rates: 15/66 (23%, mean 207 nt) in Goodier's selected
#: set; 8.8-15% genome-wide (Pickeral 2000; Szak 2003); ~9% by stricter TSD
#: criteria.  The tuple is the genome-wide survey range.
L1_TRANSDUCTION_RATE_RANGE = (0.088, 0.15)

#: SVA 5' transduction: 220/2,398 (9.2%) of TSD-bearing SVAs, ~8% overall
#: (Damert et al. 2009).  MAST2 group = >= 32% of SVA-F members (78 members).
SVA_TRANSDUCTION_RATE = 0.08
SVA_MAST2_GROUP_FRACTION_SVA_F = 0.32

# ---------------------------------------------------------------------------
# Tier 2d — Complex-event flags (Wagstaff 2012; Gilbert 2005)
# ---------------------------------------------------------------------------
#: Fraction of de novo Alu lacking normal TPRT hallmarks: 8/226 = 3.5%
#: (6 chimeric/recombination-like, 2 EN-site-and-repeat-negative).  The paper
#: does not define a "two-tailed" class — do not round to "~4% two-tailed".
ALU_COMPLEX_EVENT_RATE = 8 / 226

#: Gilbert et al. 2005: 19 of 94 5'-truncated events carried inversion/deletion
#: structures.
L1_TRUNCATED_INVERSION_RATE = 19 / 94


class ComplexEventKind(str, Enum):
    """Structurally distinct complex-event classes.

    Each names a *structural* observation, not a mechanism.  MELT annotates
    5' inversions and orientation as routine call features; no published
    method uses inversion signatures as standalone call validation, so these
    flags stay descriptive.
    """

    CHIMERIC = "chimeric"
    RECOMBINATION_LIKE = "recombination_like"
    EN_AND_REPEAT_NEGATIVE = "en_and_repeat_negative"
    INVERSION = "inversion"
    DELETION = "deletion"


@dataclass
class ComplexEventFlags:
    """Raw complex-event observations for one insertion call.

    ``None`` means "not computed"; the consumer decides what that means.
    ``inversion_fraction`` is the fraction of the inserted body aligning
    antisense to its consensus (the MELT-style 5' inversion observation).
    """

    inserted_len_bp: int | None = None
    expected_len_bp: int | None = None
    inversion_fraction: float | None = None
    chimeric: bool = False
    recombination_like: bool = False
    en_site_negative: bool = False
    repeat_negative: bool = False
    notes: str = ""
    kinds: frozenset[ComplexEventKind] = field(default_factory=frozenset)

    def compute_kinds(self) -> frozenset[ComplexEventKind]:
        """Derive flag kinds from the raw observations."""
        kinds: set[ComplexEventKind] = set()
        if self.inversion_fraction is not None and self.inversion_fraction > 0:
            kinds.add(ComplexEventKind.INVERSION)
        if (
            self.inserted_len_bp is not None
            and self.expected_len_bp is not None
            and self.inserted_len_bp < 0.3 * self.expected_len_bp
        ):
            kinds.add(ComplexEventKind.DELETION)
        if self.chimeric:
            kinds.add(ComplexEventKind.CHIMERIC)
        if self.recombination_like:
            kinds.add(ComplexEventKind.RECOMBINATION_LIKE)
        if self.en_site_negative and self.repeat_negative:
            kinds.add(ComplexEventKind.EN_AND_REPEAT_NEGATIVE)
        return frozenset(kinds)


# ---------------------------------------------------------------------------
# Tier 3 — Null-side context (never positive evidence)
# ---------------------------------------------------------------------------
#: Gasior et al. 2006: 50-bp window at L1 cut sites is 32% GC vs 41% genome,
#: but 20-kb flanks are GC-neutral (41% vs 41%) with higher Alu content
#: nearby.  The AT-richness is a 50-bp footprint, not a regional preference —
#: keep the window small or the signal washes out.
NULL_GC_WINDOW_50BP = 50

#: Sultana et al. 2019: nucleosome depletion at L1 target sites is largely
#: explained by AT-richness of the motif itself — use in the null, not as a
#: feature.
NULL_NUCLEOSOME_DEPLETION_MEDIATED_BY_AT = True

#: Flasch et al. 2019 directions (no effect sizes in text): replication-timing
#: effects are cell-type dependent (late-replicating enriched in NPCs and
#: PA-1, depleted in hESCs); EN-competent L1 cDNA integrates preferentially
#: on leading-strand templates; in FANCD2-deficient cells EN-deficient L1
#: flips to lagging strand.  Mechanism context, not detectors — the module
#: intentionally exposes no scoring function for either.
NULL_REPLICATION_TIMING_CELL_TYPE_DEPENDENT = True
NULL_LEADING_STRAND_BIAS_DIRECTION = "leading_strand"

#: Raiz et al. 2011: de novo SVA insertions prefer G+C-rich regions (direction
#: verified; no percentage in text) — null-side only.
NULL_SVA_GC_PREFERENCE_DIRECTION_ONLY = True

# ---------------------------------------------------------------------------
# Alu internal hotspots (Levy et al. 2010; saved HG03086 profile; Nummi 2025)
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class AluHotspotPosition:
    """A named Alu-internal insertion hotspot region.

    Coordinates are on the Alu consensus axis, 0-based half-open, 5'-to-3'.
    """

    name: str
    start: int
    end: int
    source: str
    note: str = ""


#: Levy et al. 2010: prominent Alu-into-Alu hotspot just after the A-rich
#: linker at consensus coordinate ~133; AluJo linker at 118-136.  The saved
#: HG03086 profile peaks at 120-140 bp (29 calls, 3.0x uniform) and 280-300
#: bp (36, 3.7x); Nummi et al. 2025 independently reports 100-155 and 260-310
#: on long-read germline data.
ALU_HOTSPOT_POSITIONS: tuple[AluHotspotPosition, ...] = (
    AluHotspotPosition(
        name="a_rich_linker",
        start=118,
        end=140,
        source=(
            "Levy et al. 2010 (linker 118-136, hotspot ~133); "
            "HG03086 120-140 bp peak 3.0x; Nummi et al. 2025 100-155"
        ),
        note="coordinate-compatible reproduction at 20-bp resolution, not base-level",
    ),
    AluHotspotPosition(
        name="3p_a_tail",
        start=280,
        end=300,
        source=(
            "HG03086 280-300 bp peak 3.7x; Nummi et al. 2025 260-310 region"
        ),
        note=(
            "Levy's positional analysis excluded consensus ends; exploratory, "
            "TPRT-hallmark validation pending"
        ),
    ),
)

#: Taper length (bp) outside a hotspot interval over which the likelihood
#: decays linearly to 0.
HOTSPOT_FALLOFF_BP = 10


def levy_hotspot_likelihood(pos: float) -> float:
    """Score proximity of a host-relative position to named Alu hotspots.

    1.0 inside any hotspot interval, linear falloff to 0 over
    :data:`HOTSPOT_FALLOFF_BP` bp outside, taking the maximum across
    intervals.  A convenience score for prioritization — never a p-value.
    """
    likelihood = 0.0
    for hs in ALU_HOTSPOT_POSITIONS:
        if hs.start <= pos < hs.end:
            return 1.0
        dist = _edge_distance(pos, hs.start, hs.end)
        likelihood = max(likelihood, max(0.0, 1.0 - dist / HOTSPOT_FALLOFF_BP))
    return likelihood


def _edge_distance(pos: float, start: int, end: int) -> float:
    """Distance from *pos* to the nearest edge of [start, end); 0 inside."""
    if pos < start:
        return start - pos
    if pos >= end:
        return pos - end + 1
    return 0.0


# ---------------------------------------------------------------------------
# Route classification — combining the markers
# ---------------------------------------------------------------------------
_ROUTES = (
    "clean_tprt",
    "l1_truncated_route",
    "internal_priming_suspect",
    "complex",
    "unknown",
)


@dataclass
class RouteEvidence:
    """All marker evidence for one insertion call, in one place.

    Built by the caller (read-backed breakpoint reconstruction first), then
    classified by :func:`classify_route`.  ``None`` always means "not
    computed" so downstream code can distinguish "no signal" from "no data".
    """

    family: str = ""
    tsd_len_bp: int | None = None
    tsd_sequence: str = ""
    en_motif_present: bool | None = None
    poly_a_tail_present: bool | None = None
    poly_a_len_bp: int | None = None
    element_3p_truncated: bool | None = None
    element_5p_truncated: bool | None = None
    microhomology_bp: int | None = None
    inserted_len_bp: int | None = None
    complex_kinds: frozenset[ComplexEventKind] = field(default_factory=frozenset)

    def to_dict(self) -> dict[str, object]:
        """Flat dict for CSV/JSON serialization (kinds sorted for stable output)."""
        return {
            "family": self.family,
            "tsd_len_bp": self.tsd_len_bp,
            "tsd_sequence": self.tsd_sequence,
            "en_motif_present": self.en_motif_present,
            "poly_a_tail_present": self.poly_a_tail_present,
            "poly_a_len_bp": self.poly_a_len_bp,
            "element_3p_truncated": self.element_3p_truncated,
            "element_5p_truncated": self.element_5p_truncated,
            "microhomology_bp": self.microhomology_bp,
            "inserted_len_bp": self.inserted_len_bp,
            "complex_kinds": (
                ",".join(sorted(k.value for k in self.complex_kinds)) or None
            ),
        }


@dataclass
class RouteClassification:
    """Route verdict plus the reasons that produced it.

    ``unknown`` is a first-class outcome: a call with insufficient evidence is
    never silently folded into a route claim (provenance discipline).
    """

    route: str
    reasons: tuple[str, ...]

    def __post_init__(self) -> None:
        if self.route not in _ROUTES:
            raise ValueError(f"unknown route {self.route!r}")


def classify_route(ev: RouteEvidence) -> RouteClassification:
    """Classify the likely insertion route from marker evidence.

    Priority ladder (first match wins):

    1. internal-priming suspect — full Srikanta signature;
    2. complex — any complex-event kind flagged;
    3. L1-type truncated route — family LINE1 + 5' truncation + (5'
       microhomology OR missing polyA tail); Zingler's 62% vs 24%
       discriminator;
    4. clean TPRT — TSD in the Wagstaff range + EN motif + polyA tail + no
       microhomology + no 5' truncation;
    5. unknown otherwise.

    None of these verdicts substitutes for read-backed breakpoint
    reconstruction; ``unknown`` is the honest default.
    """
    reasons: list[str] = []

    ip = internal_priming_signature(
        tsd_len_bp=ev.tsd_len_bp,
        poly_a_tail_present=ev.poly_a_tail_present,
        element_3p_truncated=ev.element_3p_truncated,
        en_motif_present=ev.en_motif_present,
    )
    if ip.is_suspect:
        reasons.extend(ip.reasons)
        return RouteClassification("internal_priming_suspect", tuple(reasons))

    if ev.complex_kinds:
        reasons.append(
            "complex_event_kinds=" + ",".join(sorted(k.value for k in ev.complex_kinds))
        )
        return RouteClassification("complex", tuple(reasons))

    if ev.family.upper() == "LINE1" and ev.element_5p_truncated is True:
        mh = ev.microhomology_bp or 0
        if mh > 0:
            reasons.append(
                f"5' microhomology ({mh} bp) at truncated L1 — Zingler 2005 route tag "
                f"(62% of truncated L1 vs 24% full-length)"
            )
            return RouteClassification("l1_truncated_route", tuple(reasons))
        if ev.poly_a_tail_present is False:
            reasons.append("truncated L1 without 5' microhomology and without polyA tail")
            return RouteClassification("l1_truncated_route", tuple(reasons))

    if (
        ev.tsd_len_bp is not None
        and TSD_MIN_BP <= ev.tsd_len_bp <= TSD_MAX_BP
        and ev.en_motif_present is True
        and ev.poly_a_tail_present is True
        and not (ev.microhomology_bp or 0) > 0
        and ev.element_5p_truncated is False
    ):
        reasons.append("TSD + polyA + EN motif triad present, no 5' microhomology")
        return RouteClassification("clean_tprt", tuple(reasons))

    return RouteClassification("unknown", tuple(reasons))


# ---------------------------------------------------------------------------
# Motif scanning
# ---------------------------------------------------------------------------

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    """Reverse complement of an ACGTN string (case-preserving per base)."""
    return seq.translate(_COMPLEMENT)[::-1]


def scan_en_motif_strand(seq: str) -> list[int]:
    """0-based start positions of canonical EN motif 6-mers in *seq*.

    Scans for 5'-TTTT/AA-3' written as the contiguous 6-mer ``TTTTAA`` on
    the given strand (nick between offsets 3 and 4 of the hexamer).  For the
    nick-offset convention of the empirical 7-mer consensus (TTTTT/AA), see
    :data:`EN_NICK_OFFSET`.  Callers scanning a full duplex should run this
    on both the sequence and its reverse complement.

    For calibrated per-site preference, use :func:`en_motif_score` with the
    Flasch 2019 weight table; this hexamer scan stays the structural
    detection fallback that remains computable on every sequence.
    """
    if not seq:
        return []
    s = seq.upper()
    matches: list[int] = []
    for i in range(max(0, len(s) - _EN_MOTIF_WINDOW + 1)):
        if s[i : i + _EN_MOTIF_WINDOW] == _EN_MOTIF_6MER:
            matches.append(i)
    return matches


def en_motif_score(
    seq: str, *, corrected: bool = True
) -> tuple[float | None, str]:
    """Best Flasch-calibrated site weight in *seq* plus its evidence tier.

    Returns ``(weight, tier)`` where tier is one of:

    * ``"table"`` — a 7-mer window matches a published top-set site; weight
      is that site's model weight (0.2186..1.0 in the corrected model).
    * ``"fallback"`` — no table site, but the canonical TTTTAA hexamer is
      present; weight is None.  Flasch's Fig. 2E: 45% of observed sites are
      exactly this class (a single C in positions 2-5), so the fallback is a
      real part of the spectrum — just unranked by the published table.
    * ``"none"`` — no EN-like site found; weight is None.

    Only same-tier comparisons are calibrated.  Sites outside the published
    top 21 of 743 are unverifiable, not zero-preference.
    """
    if not seq:
        return (None, "none")
    s = seq.upper()
    table = en_motif_scan_table(corrected=corrected)
    best: float | None = None
    for site in table:
        start = 0
        while True:
            idx = s.find(site, start)
            if idx < 0:
                break
            weight = en_motif_site_weight(site, corrected=corrected)
            if weight is not None and (best is None or weight > best):
                best = weight
            start = idx + 1
    if best is not None:
        return (best, "table")
    if scan_en_motif_strand(s):
        return (None, "fallback")
    return (None, "none")


def en_motif_opportunity(seq: str) -> float:
    """Fraction of 6-mer windows in *seq* matching the canonical EN motif.

    An *opportunity* measure for null-model construction, not an enrichment
    statistic.  Null models must stay composition-aware (see Tier 3 notes).
    For a Flasch-style weighted null, sample sites with the weight tables
    (:func:`en_motif_scan_table`) instead of thresholding on this hexamer.
    """
    windows = max(0, len(seq) - _EN_MOTIF_WINDOW + 1)
    if windows == 0:
        return 0.0
    return len(scan_en_motif_strand(seq)) / windows


def tsd_len_in_range(n: int | None) -> bool:
    """True if *n* is within the Wagstaff de novo TSD range [5, 27] bp."""
    return n is not None and TSD_MIN_BP <= n <= TSD_MAX_BP


# ---------------------------------------------------------------------------
# 5' junction microhomology (Zingler route tag)
# ---------------------------------------------------------------------------

def microhomology_at_5p_junction(
    flank_seq: str,
    element_seq: str,
    *,
    max_len: int = MH_MAX_LEN,
) -> int:
    """Longest overlap (bp) between the end of the genomic flank and the start
    of the inserted element.

    *flank_seq* is the reference target sequence immediately 5' of the
    insertion (already oriented 5'-3' along the insertion); *element_seq*
    is the inserted body in the same orientation.  0 means clean junction;
    1-12 bp is the Zingler route tag.  NOT a twin-priming fingerprint —
    see :data:`MH_RATE_DESCRIPTION`.
    """
    if not flank_seq or not element_seq:
        return 0
    f = flank_seq.upper()
    e = element_seq.upper()
    lim = min(max_len, len(f), len(e))
    for k in range(lim, 0, -1):
        if f[-k:] == e[:k]:
            return k
    return 0


# ---------------------------------------------------------------------------
# Internal priming (Srikanta et al. 2009 signature)
# ---------------------------------------------------------------------------

@dataclass
class InternalPrimingSignature:
    """Result of the internal-priming test for one call."""

    is_suspect: bool
    reasons: tuple[str, ...]


def internal_priming_signature(
    *,
    tsd_len_bp: int | None,
    poly_a_tail_present: bool | None,
    element_3p_truncated: bool | None,
    en_motif_present: bool | None,
) -> InternalPrimingSignature:
    """Detect the Srikanta et al. 2009 internal-priming architecture.

    Signature: 3'-truncated element, **no** poly(A) tail, TSD >= 6 bp
    retained, no canonical EN-site preference.  The features are individually
    common — only the joint architecture is the test, so every component must
    be computed (``None`` blocks the signature rather than implying absence).
    """
    required: tuple[tuple[str, bool], ...] = (
        ("3'-truncated element", element_3p_truncated is True),
        ("no poly(A) tail", poly_a_tail_present is False),
        (
            f"TSD >= {INTERNAL_PRIMING_MIN_TSD} bp",
            tsd_len_bp is not None and tsd_len_bp >= INTERNAL_PRIMING_MIN_TSD,
        ),
        ("no canonical EN site", en_motif_present is False),
    )
    reasons = [label for label, ok in required if ok]
    is_suspect = all(ok for _, ok in required)
    if is_suspect:
        reasons = ["Srikanta internal-priming signature complete"]
    elif reasons:
        reasons.append(
            "signature incomplete — features are individually common, "
            "joint architecture is the test"
        )
    return InternalPrimingSignature(is_suspect=is_suspect, reasons=tuple(reasons))


# ---------------------------------------------------------------------------
# Complex-event flags
# ---------------------------------------------------------------------------

def flag_complex_event(
    *,
    inserted_len_bp: int | None = None,
    expected_len_bp: int | None = None,
    inversion_fraction: float | None = None,
    chimeric: bool = False,
    recombination_like: bool = False,
    en_site_negative: bool = False,
    repeat_negative: bool = False,
) -> ComplexEventFlags:
    """Build :class:`ComplexEventFlags` and derive ``kinds`` in one step."""
    flags = ComplexEventFlags(
        inserted_len_bp=inserted_len_bp,
        expected_len_bp=expected_len_bp,
        inversion_fraction=inversion_fraction,
        chimeric=chimeric,
        recombination_like=recombination_like,
        en_site_negative=en_site_negative,
        repeat_negative=repeat_negative,
    )
    flags.kinds = flags.compute_kinds()
    return flags


# ---------------------------------------------------------------------------
# Null-side context
# ---------------------------------------------------------------------------

def at_richness_footprint(window_seq: str) -> float:
    """AT fraction of a *small* window at the nick.

    Gasior et al. 2006: the AT-richness is a 50-bp footprint, not a regional
    preference — 20-kb flanks are GC-neutral.  Pass a window around the cut
    site (:data:`NULL_GC_WINDOW_50BP`); larger windows wash the signal out.
    """
    s = (window_seq or "").upper()
    n = len(s)
    if n == 0:
        return 0.0
    n_at = sum(1 for c in s if c in {"A", "T"})
    return n_at / n
