"""Reproduce the locus-vs-global study end to end."""
import json
import numpy as np
import pandas as pd
from l1meth import catalog as C, compare as K, panel as P, simulate as S


def main() -> None:
    seq = json.load(open("data/chr22_seq.json"))["dna"]
    l1 = pd.read_csv("l1_locus_catalog_chr22.csv")
    n_addr = int(l1.addressable.sum())
    print(f"addressable loci: {n_addr}  promoter CpGs: {int(l1[l1.addressable].n_promoter_cpg.sum())}")

    for f in (0.0, 0.05, 0.10, 0.20):
        aucs = [K.compare(*S.simulate(S.SimConfig(n_loci=n_addr, signal_fraction=f, seed=s))[:2], seed=s)
                for s in range(10)]
        g = np.mean([a["auc_global"] for a in aucs])
        lo = np.mean([a["auc_locus"] for a in aucs])
        print(f"  signal_fraction={f:<5} global={g:.3f}  locus={lo:.3f}  delta={lo-g:+.3f}")

    for strat in ("random", "youngest"):
        r = P.panel_hit_rate(P.PoolConfig(signal_fraction=0.10), 11, strat)
        print(f"  11-locus panel [{strat:<8}] P(contains signal)={r['p_any_signal']:.3f}"
              f"  mean signal loci={r['mean_signal_loci']:.2f}")


if __name__ == "__main__":
    main()
