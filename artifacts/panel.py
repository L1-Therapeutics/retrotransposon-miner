"""Which loci you measure, not how many.

Earlier sweeps compared a bulk average against every locus at once. Real
locus-specific studies measure a *panel* -- the published multi-cancer null
used 11 loci -- and a panel chosen without knowing where signal lives can miss
it entirely. This module separates two questions that are easily conflated:

  1. does locus-resolved measurement help at all?  (see ``compare``)
  2. does it help *when you had to choose the loci in advance*?

Selection is where an insertion caller could earn its place: divergence,
5'UTR retention and flank anchorability are all computable from reference data
alone, so a panel can in principle be enriched before any sample is sequenced.

**Whether that enrichment works depends entirely on ``age_enrichment``, and for
tumour/normal the evidence says it does not.** Lanciano et al. 2024 analysed
~12,000 Illumina 450K samples from GEO (young L1PA probes n = 695, old
n = 189) and found the methylation of young L1PAs "high and generally similar
to or higher than that of old L1PAs in most situations." Young L1PAs are
hypomethylated specifically in pluripotent stem cells, trophoblast, embryonal
carcinoma, seminoma, placenta, fetal membranes and hydatidiform moles -- early
embryonic, extra-embryonic and male germline contexts. In normal fibroblasts
and several carcinoma lines (MCF-7, HeLa-S3, HepG2) young L1PAs are relatively
*hyper*methylated.

So the default here is ``age_enrichment=1.0``: no enrichment, which is what the
array data supports for carcinomas. At that setting age-informed selection gains no advantage within noise versus
blind selection, and the panel-selection advantage this module was written to
demonstrate disappears. In the corrected 1,500-replicate comparison the hit
rates were 70.4% blind and 68.2% age-informed (about 1.7 percentage points SE
on their difference): this is not a meaningful negative effect estimate. Values above 1 are defensible
only for germ-cell/embryonal tumours, where the young-L1 hypomethylation the
assumption needs is actually observed.

One caveat on the mapping: the array study compares young-vs-old methylation
*levels within a sample*, whereas panel selection needs young-vs-old
*differential* methylation between tumour and normal. Those are related but
not identical, so this is strong evidence against the assumption rather than a
direct refutation of it.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class PoolConfig:
    """A genome-scale pool of addressable loci with age-structured signal."""

    n_pool: int = 2600            # chr22 catalogue scaled genome-wide
    signal_fraction: float = 0.10
    age_enrichment: float = 1.0   # signal odds in the youngest loci relative to
    #   the oldest. 1.0 = none, the value the 450K evidence supports for
    #   carcinomas (see module docstring). Raise it only for germ-cell or
    #   embryonal contexts; ~8.0 reproduces the falsified optimistic case this study
    #   originally assumed, and is retained only to show what it was worth.
    seed: int = 0


def make_pool(cfg: PoolConfig) -> tuple[np.ndarray, np.ndarray]:
    """Return (age_rank, is_signal) for the pool.

    ``age_rank`` is 0 (youngest, lowest divergence) to 1 (oldest). Signal
    probability decays with age, so young promoter-intact copies are the more
    likely carriers -- the assumption the selection strategy is betting on.
    """
    rng = np.random.default_rng(cfg.seed)
    age = rng.uniform(0, 1, cfg.n_pool)
    # weight decays from `age_enrichment` at age 0 to 1.0 at age 1
    w = cfg.age_enrichment ** (1.0 - age)
    p = w / w.sum() * cfg.signal_fraction * cfg.n_pool
    is_signal = rng.uniform(0, 1, cfg.n_pool) < np.clip(p, 0, 1)
    return age, is_signal


def select_panel(age: np.ndarray, k: int, *, strategy: str, seed: int = 0) -> np.ndarray:
    """Indices of ``k`` loci chosen by ``strategy``.

    ``random``   -- uniform choice, the blind baseline.
    ``youngest`` -- lowest divergence first, computable from reference alone.
    """
    rng = np.random.default_rng(seed)
    if strategy == "random":
        return rng.choice(len(age), size=k, replace=False)
    if strategy == "youngest":
        return np.argsort(age)[:k]
    raise ValueError(f"unknown strategy: {strategy}")


def panel_hit_rate(cfg: PoolConfig, k: int, strategy: str, n_rep: int = 400) -> dict[str, float]:
    """How often does a ``k``-locus panel actually contain signal loci?"""
    hits, counts = [], []
    for rep in range(n_rep):
        age, sig = make_pool(PoolConfig(**{**cfg.__dict__, "seed": rep}))
        idx = select_panel(age, k, strategy=strategy, seed=rep)
        n_sig = int(sig[idx].sum())
        counts.append(n_sig)
        hits.append(n_sig > 0)
    return {
        "strategy": strategy,
        "k": k,
        "p_any_signal": float(np.mean(hits)),
        "mean_signal_loci": float(np.mean(counts)),
    }
