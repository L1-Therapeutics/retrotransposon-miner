"""Generate per-locus L1 promoter methylation for tumour and normal samples.

The generative model mirrors how the two competing assays actually relate:
a locus-specific assay observes the per-locus vector, while the conventional
bulk-bisulfite assay observes (approximately) its mean. Simulating the vector
and averaging it is therefore the honest way to compare them -- both readouts
come from one underlying biology.

Two independent sources of tumour/normal difference are modelled, because the
literature supports both:

* a **sample-level** shift -- global hypomethylation, the effect the bulk assay
  was designed around and which it measures efficiently;
* a **locus-specific** shift confined to a subset of loci -- the effect that
  motivates locus resolution, and which the bulk mean dilutes by ~1/n_loci.

Calibration (verified against sources, 2026-09-19):
  Baba et al. 2010, Mol Cancer 9:125 (n=869 colorectal tumours, pyrosequencing)
  -- tumour LINE-1 methylation mean 61.4, SD 9.6, range 23.1-90.3, approximately
  normal apart from an excess of extreme hypomethylators (<40; 2.5%).
  Normal colonic mucosa is taken at ~76% with a tighter SD; this figure is the
  weakest link in the calibration and is exposed as a parameter.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

# --- literature-calibrated defaults (percentage points on a 0-100 scale) ---
NORMAL_MEAN = 76.0
NORMAL_SD = 4.5
TUMOUR_MEAN = 61.4   # Baba 2010
TUMOUR_SD = 9.6      # Baba 2010


@dataclass
class SimConfig:
    n_loci: int = 35              # chr22 addressable catalogue; ~2100 genome-wide
    n_normal: int = 60
    n_tumour: int = 60
    signal_fraction: float = 0.15  # fraction of loci with locus-specific shift
    locus_effect: float = 25.0     # extra methylation loss at those loci (pp)
    global_shift: float = 14.6     # NORMAL_MEAN - TUMOUR_MEAN
    sample_sd_normal: float = NORMAL_SD
    sample_sd_tumour: float = TUMOUR_SD
    locus_sd: float = 8.0          # biological variation between loci
    coverage: int = 10             # informative reads per junction (30x WGBS)
    seed: int = 0


def simulate(cfg: SimConfig) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (X, y, is_signal_locus).

    ``X`` is (n_samples, n_loci) observed methylation percentage;
    ``y`` is 1 for tumour. Read sampling is binomial at ``cfg.coverage``, so
    per-locus estimates carry realistic measurement noise -- the cost that
    partly offsets locus resolution's advantage.
    """
    rng = np.random.default_rng(cfg.seed)
    n_sig = int(round(cfg.signal_fraction * cfg.n_loci))
    is_signal = np.zeros(cfg.n_loci, dtype=bool)
    is_signal[rng.choice(cfg.n_loci, size=n_sig, replace=False)] = True

    # Per-locus baseline methylation in normal tissue.
    locus_base = np.clip(rng.normal(NORMAL_MEAN, cfg.locus_sd, cfg.n_loci), 5, 98)

    rows, labels = [], []
    for grp, n, sd in (("normal", cfg.n_normal, cfg.sample_sd_normal),
                       ("tumour", cfg.n_tumour, cfg.sample_sd_tumour)):
        for _ in range(n):
            true = locus_base.copy()
            if grp == "tumour":
                true = true - cfg.global_shift          # sample-level component
                true[is_signal] -= cfg.locus_effect     # locus-specific component
            true = true + rng.normal(0, sd)             # sample random effect
            true = np.clip(true, 1, 99)
            # read sampling: binomial with `coverage` informative reads
            obs = rng.binomial(cfg.coverage, true / 100.0) / cfg.coverage * 100.0
            rows.append(obs)
            labels.append(1 if grp == "tumour" else 0)

    return np.asarray(rows), np.asarray(labels), is_signal
