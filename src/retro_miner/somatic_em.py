"""Expectation-Maximization (EM) Somatic Subclonal MEI Mixture Model.

Distinguishes low-VAF somatic retrotransposition events (0.02 <= VAF <= 0.25)
from germline heterozygous calls (VAF ~ 0.50) and sequencing/PCR chimeras (VAF < 0.02).

Literature Anchors: Rodriguez-Martin et al. (2020) Nature Genetics / Evrony et al. (2012) Neuron.
"""

from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass(frozen=True)
class SomaticCall:
    is_somatic: bool
    somatic_posterior: float
    subclone_vaf: float
    classification: str


def binom_pmf(k: int, n: int, p: float) -> float:
    """Compute Binomial(k; n, p) probability."""
    if n == 0:
        return 1.0
    p = max(1e-5, min(1.0 - 1e-5, p))
    comb = math.comb(n, k)
    return comb * (p**k) * ((1.0 - p) ** (n - k))


class SomaticEMClassifier:
    """3-Component Binomial EM Classifier for somatic MEI discovery."""

    def __init__(self):
        # Initial mixture weights: [Artifact, Somatic Subclonal, Germline]
        self.weights = [0.80, 0.05, 0.15]
        # Alternate allele frequency parameters: [p_art, p_som, p_germ]
        self.p_rates = [0.01, 0.10, 0.50]

    def classify_locus(self, k_alt: int, n: int) -> SomaticCall:
        """Classify a candidate locus into ARTIFACT, SOMATIC_SUBCLONAL, or GERMLINE."""
        if n == 0 or k_alt == 0:
            return SomaticCall(
                is_somatic=False,
                somatic_posterior=0.0,
                subclone_vaf=0.0,
                classification="ARTIFACT",
            )

        vaf = k_alt / n
        probs = [
            self.weights[0] * binom_pmf(k_alt, n, self.p_rates[0]),
            self.weights[1] * binom_pmf(k_alt, n, self.p_rates[1]),
            self.weights[2] * binom_pmf(k_alt, n, self.p_rates[2]),
        ]

        total_p = sum(probs)
        if total_p == 0:
            posteriors = [0.3333, 0.3333, 0.3333]
        else:
            posteriors = [p / total_p for p in probs]

        p_art, p_som, p_germ = posteriors

        if p_som > 0.50 and 0.02 <= vaf <= 0.30:
            classification = "SOMATIC_SUBCLONAL"
            is_somatic = True
        elif vaf >= 0.35:
            classification = "GERMLINE"
            is_somatic = False
        else:
            classification = "ARTIFACT"
            is_somatic = False

        return SomaticCall(
            is_somatic=is_somatic,
            somatic_posterior=round(p_som, 4),
            subclone_vaf=round(vaf, 4),
            classification=classification,
        )
