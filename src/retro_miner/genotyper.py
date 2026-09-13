"""Bayesian MEI Genotyper & Allele Balance Model.

Calculates diploid genotype likelihoods (0/0, 0/1, 1/1), Phred-scaled genotype
quality (GQ), and Variant Allele Frequency (VAF) from alt evidence reads (k_alt)
and spanning concordant reference reads (k_ref).

Literature Anchor: Li (2011) Bioinformatics / Garrison & Marth (2012) FreeBayes model.
"""

from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass(frozen=True)
class GenotypeCall:
    genotype: str
    vaf: float
    log_likelihood_ratio: float
    genotype_quality: float
    posterior_probs: dict[str, float]


_MIN_TOTAL_READS = 3


def log_binomial_pmf(k: int, n: int, p: float) -> float:
    """Calculate log2 binomial probability log2 Binomial(k; n, p)."""
    if n == 0:
        return 0.0
    if p <= 0.0:
        p = 1e-6
    if p >= 1.0:
        p = 1.0 - 1e-6

    # Log binomial coefficient using lgamma
    log_comb = math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
    log_comb_2 = log_comb / math.log(2)

    log_prob = k * math.log2(p) + (n - k) * math.log2(1.0 - p)
    return log_comb_2 + log_prob


def calculate_mei_genotype(
    k_alt: int,
    k_ref: int,
    seq_error_rate: float = 0.02,
) -> GenotypeCall:
    """Compute Bayesian posterior genotype likelihoods and VAF for an MEI locus.

    Args:
        k_alt: Number of supporting non-reference evidence reads (split, discordant, clip).
        k_ref: Number of spanning concordant reference reads across breakpoint.
        seq_error_rate: Base sequencing/mapping error rate epsilon (default 0.02).

    Returns:
        GenotypeCall containing call string ("0/0", "0/1", "1/1", or "./."), VAF, GQ.
    """
    n = k_alt + k_ref
    if n < _MIN_TOTAL_READS:
        return GenotypeCall(
            genotype="./.",
            vaf=0.0,
            log_likelihood_ratio=0.0,
            genotype_quality=0.0,
            posterior_probs={"0/0": 0.3333, "0/1": 0.3333, "1/1": 0.3333},
        )

    vaf = round(k_alt / n, 4)

    # Expected alt frequencies under diploid hypothesis
    p_g00 = seq_error_rate
    p_g01 = 0.50
    p_g11 = 1.0 - seq_error_rate

    log_l_00 = log_binomial_pmf(k_alt, n, p_g00)
    log_l_01 = log_binomial_pmf(k_alt, n, p_g01)
    log_l_11 = log_binomial_pmf(k_alt, n, p_g11)

    log_likelihoods = {"0/0": log_l_00, "0/1": log_l_01, "1/1": log_l_11}

    # Normalize log likelihoods to probabilities via log-sum-exp
    max_log = max(log_likelihoods.values())
    unnorm_probs = {g: 2.0 ** (log_l - max_log) for g, log_l in log_likelihoods.items()}
    total_prob = sum(unnorm_probs.values())
    posteriors = {g: round(prob / total_prob, 6) for g, prob in unnorm_probs.items()}

    sorted_gt = sorted(posteriors.items(), key=lambda x: x[1], reverse=True)
    top_gt, top_p = sorted_gt[0]
    runner_up_gt, runner_up_p = sorted_gt[1]

    # Phred-scaled Genotype Quality: GQ = -10 * log10(1 - P(G_max))
    error_p = max(1e-10, 1.0 - top_p)
    gq = round(min(99.0, -10.0 * math.log10(error_p)), 2)

    # Log Likelihood Ratio between top two calls
    llr = round(log_likelihoods[top_gt] - log_likelihoods[runner_up_gt], 4)

    return GenotypeCall(
        genotype=top_gt,
        vaf=vaf,
        log_likelihood_ratio=llr,
        genotype_quality=gq,
        posterior_probs=posteriors,
    )
