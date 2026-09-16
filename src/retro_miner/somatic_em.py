"""Expectation-Maximization Subclonal Somatic MEI Mixture Model.

Standard diploid genotypers (Li 2011) assume biallelic diploidy (0/0, 0/1,
1/1) and therefore fail in cancer genomes, neurodevelopmental tissue, and
bulk-tissue somatic mosaicism, where retrotransposition occurs subclonally
in a fraction of cells: a somatic insertion with 3 alt reads out of 60
total reads (VAF = 0.05) collapses onto homozygous-reference under diploid
likelihoods.

This module solves a 3-component Binomial Mixture Model by Expectation-
Maximization over three latent states -- chimeric PCR artifact
(``theta_art`` ~ 0.01), somatic subclone (``theta_som`` learnable,
biologically VAF in [0.02, 0.25]), and germline diploid (``theta_germ``
~ 0.50):

    L(pi, theta) = sum_i log( sum_k pi_k Binom(k_alt_i | n_i, theta_k) )

The E-step computes component responsibilities gamma_{i,k}; the M-step
re-estimates the mixture weights and the somatic allele fraction
``theta_som`` by gamma-weighted alt-read pooling.

Literature Anchors: Rodriguez-Martin et al. (2020) Nature Genetics
52(3):306-319; Evrony et al. (2012) Neuron 75(6):1017-1027.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from retro_miner.genotyper import log_binomial_pmf

ARTIFACT = "ARTIFACT"
SOMATIC_SUBCLONAL = "SOMATIC_SUBCLONAL"
GERMLINE = "GERMLINE"

_COMPONENTS = (ARTIFACT, SOMATIC_SUBCLONAL, GERMLINE)
_SOMATIC_INDEX = 1
_N_COMPONENTS = 3

_THETA_SOM_LOWER = 0.02
_THETA_SOM_UPPER = 0.45
_LOG_PI_FLOOR = 1e-300
_P_MIN = 1e-5

_INVALID_INPUT_MSG = "k_alt_list and n_list must be parallel non-empty sequences"
_UNFITTED_MSG = "SomaticEMClassifier must be fit before predict() is called"


@dataclass(frozen=True)
class SomaticCall:
    is_somatic: bool
    somatic_posterior: float
    subclone_vaf: float
    classification: str


def binom_pmf(k: int, n: int, p: float) -> float:
    """Compute the Binomial (k; n, p) probability mass function."""
    if n == 0:
        return 1.0
    p = max(_P_MIN, min(1.0 - _P_MIN, p))
    comb = math.comb(n, k)
    return comb * (p**k) * ((1.0 - p) ** (n - k))


class SomaticEMClassifier:
    """Fit and apply a 3-component Binomial mixture for subclonal MEI calls.

    Component order is fixed as ``(ARTIFACT, SOMATIC_SUBCLONAL, GERMLINE)``.
    ``theta_art`` and ``theta_germ`` are held at their fixed biological
    priors while ``theta_som``, the latent subclonal allele fraction, is
    learned from the data in the M-step.
    """

    def __init__(
        self,
        pi_art: float = 0.80,
        pi_som: float = 0.05,
        pi_germ: float = 0.15,
        theta_art: float = 0.01,
        theta_som: float = 0.10,
        theta_germ: float = 0.50,
    ) -> None:
        self._pi: list[float] = [float(pi_art), float(pi_som), float(pi_germ)]
        self._theta: list[float] = [float(theta_art), float(theta_som), float(theta_germ)]
        self._n_iter: int = 0
        self._converged: bool = False
        self._fitted: bool = False

    @property
    def is_fitted(self) -> bool:
        return self._fitted

    @property
    def theta_som(self) -> float:
        return self._theta[_SOMATIC_INDEX]

    @property
    def n_iter(self) -> int:
        return self._n_iter

    @property
    def converged(self) -> bool:
        return self._converged

    @property
    def pi(self) -> dict[str, float]:
        return {name: self._pi[i] for i, name in enumerate(_COMPONENTS)}

    def fit_predict(
        self,
        k_alt_list: list[int],
        n_list: list[int],
        max_iter: int = 100,
        tol: float = 1e-5,
    ) -> list[SomaticCall]:
        """Fit the mixture by EM and return one :class:`SomaticCall` per locus.

        Args:
            k_alt_list: Alternate (non-reference) supporting read counts.
            n_list: Total read depths per locus (reference + alternate).
            max_iter: Maximum number of EM iterations.
            tol: Log-likelihood convergence threshold |L(t) - L(t-1)| < tol.

        Returns:
            Parallel list of SomaticCall objects using the fitted model.
        """
        if len(k_alt_list) != len(n_list) or not k_alt_list:
            raise ValueError(_INVALID_INPUT_MSG)

        k_alt = [min(int(k), int(n)) for k, n in zip(k_alt_list, n_list, strict=True)]
        n = [max(0, int(total)) for total in n_list]

        prev_ll = math.inf
        for _iteration in range(max_iter):
            gamma = self._e_step(k_alt, n)
            log_lik = self._log_likelihood(k_alt, n)
            if abs(prev_ll - log_lik) < tol:
                self._converged = True
                break
            prev_ll = log_lik
            self._m_step(k_alt, n, gamma)
        else:
            self._converged = False

        self._n_iter = _iteration + 1
        self._fitted = True

        gamma = self._e_step(k_alt, n)
        return [self._call_from_gamma(row) for row in gamma]

    def predict(self, k_alt: int, n: int) -> SomaticCall:
        """Score a single locus against the fitted model's posterior.

        Args:
            k_alt: Alternate supporting read count for the locus.
            n: Total read depth for the locus.

        Returns:
            SomaticCall using the last fitted model parameters.

        Raises:
            RuntimeError: If the classifier has not been fit yet.
        """
        if not self._fitted:
            raise RuntimeError(_UNFITTED_MSG)
        k_alt = min(int(k_alt), int(n))
        n = max(0, int(n))
        gamma = self._e_step([k_alt], [n])[0]
        return self._call_from_gamma(gamma)

    def classify_locus(self, k_alt: int, n: int) -> SomaticCall:
        """Classify a single candidate locus with the fixed prior parameters.

        Convenience single-locus posterior classifier used to classify an
        individual locus without fitting a batch mixture (legacy API).

        Args:
            k_alt: Alternate supporting read count for the locus.
            n: Total read depth for the locus.

        Returns:
            SomaticCall flagged ARTIFACT / SOMATIC_SUBCLONAL / GERMLINE.
        """
        if n == 0 or k_alt == 0:
            return SomaticCall(
                is_somatic=False,
                somatic_posterior=0.0,
                subclone_vaf=0.0,
                classification=ARTIFACT,
            )

        vaf = k_alt / n
        probs = [
            self._pi[c] * binom_pmf(k_alt, n, self._theta[c])
            for c in range(_N_COMPONENTS)
        ]
        total_p = sum(probs)
        if total_p == 0:
            posteriors = [1 / _N_COMPONENTS] * _N_COMPONENTS
        else:
            posteriors = [p / total_p for p in probs]

        p_som = posteriors[_SOMATIC_INDEX]
        if p_som > 0.50 and 0.02 <= vaf <= 0.30:
            classification = SOMATIC_SUBCLONAL
            is_somatic = True
        elif vaf >= 0.35:
            classification = GERMLINE
            is_somatic = False
        else:
            classification = ARTIFACT
            is_somatic = False

        return SomaticCall(
            is_somatic=is_somatic,
            somatic_posterior=round(p_som, 4),
            subclone_vaf=round(vaf, 4),
            classification=classification,
        )

    def _e_step(self, k_alt: list[int], n: list[int]) -> list[list[float]]:
        """Return the responsibilities gamma_{i,k} (posterior of component k)."""
        ln2 = math.log(2.0)
        log_pi = [math.log(max(pi, _LOG_PI_FLOOR)) for pi in self._pi]
        responsibilities: list[list[float]] = []
        for k_i, total in zip(k_alt, n, strict=True):
            log_weights = [
                log_pi[c] + log_binomial_pmf(k_i, total, self._theta[c]) * ln2
                for c in range(_N_COMPONENTS)
            ]
            max_log = max(log_weights)
            log_sum = max_log + math.log(sum(math.exp(lw - max_log) for lw in log_weights))
            responsibilities.append([math.exp(lw - log_sum) for lw in log_weights])
        return responsibilities

    def _m_step(self, k_alt: list[int], n: list[int], gamma: list[list[float]]) -> None:
        """Update mixture weights and the latent somatic allele fraction."""
        n_loci = len(k_alt)
        for c in range(_N_COMPONENTS):
            self._pi[c] = sum(row[c] for row in gamma) / n_loci

        som_num = sum(
            row[_SOMATIC_INDEX] * k_i for row, k_i in zip(gamma, k_alt, strict=True)
        )
        som_den = sum(
            row[_SOMATIC_INDEX] * total for row, total in zip(gamma, n, strict=True)
        )
        if som_den > 0:
            som_rate = som_num / som_den
            self._theta[_SOMATIC_INDEX] = min(
                _THETA_SOM_UPPER, max(_THETA_SOM_LOWER, som_rate)
            )

    def _log_likelihood(self, k_alt: list[int], n: list[int]) -> float:
        """Observed-data log-likelihood sum_i log(sum_k pi_k Binom(k_i; n_i, theta_k))."""
        ln2 = math.log(2.0)
        log_pi = [math.log(max(pi, _LOG_PI_FLOOR)) for pi in self._pi]
        total = 0.0
        for k_i, total_n in zip(k_alt, n, strict=True):
            terms = [
                log_pi[c] + log_binomial_pmf(k_i, total_n, self._theta[c]) * ln2
                for c in range(_N_COMPONENTS)
            ]
            max_log = max(terms)
            total += max_log + math.log(sum(math.exp(t - max_log) for t in terms))
        return total

    def _call_from_gamma(self, gamma_row: list[float]) -> SomaticCall:
        """Build the SomaticCall for one locus from its responsibility row."""
        top_index = max(range(_N_COMPONENTS), key=lambda c: gamma_row[c])
        classification = _COMPONENTS[top_index]
        return SomaticCall(
            is_somatic=classification == SOMATIC_SUBCLONAL,
            somatic_posterior=round(gamma_row[_SOMATIC_INDEX], 6),
            subclone_vaf=round(self.theta_som, 4),
            classification=classification,
        )