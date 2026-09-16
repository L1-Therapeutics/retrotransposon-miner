"""Somatic vs. Germline Odds-Ratio and Panel-of-Normals (PON) Filter.

Discriminates low-frequency somatic mobile element insertions (MEIs) from
germline inheritance or mapping/panel noise using a Fisher exact log
odds-ratio against a matched panel-of-normals (PON).

Somatic MEIs in focal brain mosaicism and tumor genomes exhibit low variant
allele frequency (VAF <= 0.15), whereas germline MEIs sit near diploid
expectation (VAF ~ 0.50 or 1.00).  When the same alt-supporting reads appear
across the control panel, the call is more plausibly a germline polymorphism
or sequencing artifact; when the panel is clean but the sample is enriched,
the call is a de-novo somatic candidate.

The 2x2 contingency matrix layout is::

                        Sample Alt   Sample Ref
    Sample locus            k_alt         k_ref
    Panel/Control           c_alt         c_ref

with Fisher's exact test two-tailed p-value from ``scipy.stats.fisher_exact``.

Literature anchors: Evrony et al. (2012) Cell 151(6):1243-1256; Iskow et al.
(2010) Cell 141(7):1253-1261.
"""

from __future__ import annotations

from dataclasses import dataclass

from scipy.stats import fisher_exact

SOMATIC_CANDIDATE = "SOMATIC_CANDIDATE"
GERMLINE_PASS = "GERMLINE_PASS"
CONTROL_PANEL_ARTIFACT = "CONTROL_PANEL_ARTIFACT"

_MIN_SOMATIC_ALT_READS = 3
_SOMATIC_P_VALUE_THRESHOLD = 0.01
_CONTROL_ARTIFACT_ALT_THRESHOLD = 2


@dataclass(frozen=True)
class SomaticFilterResult:
    is_somatic: bool
    odds_ratio: float
    p_value: float
    classification: str


def evaluate_somatic_odds(
    k_alt: int,
    k_ref: int,
    c_alt: int = 0,
    c_ref: int = 30,
) -> SomaticFilterResult:
    """Evaluate whether a low-VAF MEI call is somatic vs. germline or PON noise.

    Args:
        k_alt: Alt-supporting reads in the target sample (split, discordant, clip).
        k_ref: Spanning reference reads in the target sample.
        c_alt: Alt-supporting reads across the matched panel-of-normals.
        c_ref: Spanning reference reads across the matched panel-of-normals.

    Returns:
        SomaticFilterResult with odds ratio, Fisher's exact two-tailed p-value,
        and a SOMATIC_CANDIDATE / GERMLINE_PASS / CONTROL_PANEL_ARTIFACT label.
    """
    k_alt = max(0, int(k_alt))
    k_ref = max(0, int(k_ref))
    c_alt = max(0, int(c_alt))
    c_ref = max(0, int(c_ref))

    if k_alt + k_ref == 0 and c_alt + c_ref == 0:
        return SomaticFilterResult(False, 0.0, 1.0, GERMLINE_PASS)

    odds_ratio = (k_alt * c_ref) / max(1, k_ref * c_alt)

    _, p_value = fisher_exact([[k_alt, k_ref], [c_alt, c_ref]], alternative="two-sided")
    p_value = float(p_value)

    if c_alt >= _CONTROL_ARTIFACT_ALT_THRESHOLD:
        classification = CONTROL_PANEL_ARTIFACT
    elif c_alt == 0 and k_alt >= _MIN_SOMATIC_ALT_READS and p_value <= _SOMATIC_P_VALUE_THRESHOLD:
        classification = SOMATIC_CANDIDATE
    else:
        classification = GERMLINE_PASS

    return SomaticFilterResult(
        is_somatic=classification == SOMATIC_CANDIDATE,
        odds_ratio=round(odds_ratio, 4),
        p_value=p_value,
        classification=classification,
    )
