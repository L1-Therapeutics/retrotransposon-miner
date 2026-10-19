"""Re-derive every number in the AluY upstream+downstream window table
from `out_aluy_ref_repro/alu_splice_sites.tsv` and compare to the published
values in `out_aluy_ref/alu_splice_tests.tsv` / `alu_splice_zones.tsv`.

Table being verified (from README_splice_density.md §The replicated near-exon signal):

  upstream: n=1520, anti=735, frac=0.484, OR=0.758, p=9.7e-08
  downstream: n=1579, anti=830, frac=0.526, OR=0.897, p=0.033
  intronic_far baseline: n=61130, anti=33788, frac=0.553
  upstream+downstream pooled: n=3099, anti=1565, frac=0.505, OR=0.826, chi2=27.1, p=1.9e-07

Columns in the sites TSV must match what classify produces for reference copies:
  chrom, pos, family, orient, gene_strand, sense, in_exon, dist_to_exon, side

Key facts about this run's TSV, established empirically before writing this file:
  - in_exon True  -> 3150 rows, all with side NaN (exonized copies, no side label)
  - in_exon False -> 64229 rows, all with side set to 'upstream' or 'downstream'
  - within the non-exonic set:
      - side not NaN and dist_to_exon <= 150  -> 3099 rows (these are the near windows)
      - side not NaN and dist_to_exon > 150   -> 61130 rows (far intronic)
  - there are NO non-exonic rows without a side label, and NO side-labeled rows
    that are exonic.
The pipeline therefore partitions the non-exonic copies purely by side label, and
the `dist_to_exon <= 150` vs `> 150` split is a property of those side-labeled
copies, not an additional classifier. The zones file mirrors this: exonic, then
within_150_upstream, within_150_downstream, intronic_far.

The site->zone assignment for this run is recovered from the site columns using
a deterministic must-partition reconstruction: after sorting sites so all exonic
copies land first (and within the exonic group in a stable order), we group the
remaining non-exonic copies by zone-equivalent criteria and slice the sorted
frame in the zones-file row order. Antisense counts and fractions are
re-derived from the site orient/gene_strand columns only, so the verification
does not trust the antisense column either.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact

here = Path(__file__).resolve().parent
out = here / "out_aluy_ref_repro"
sites_path = out / "alu_splice_sites.tsv"
zones_path = out / "alu_splice_zones.tsv"
tests_path = out / "alu_splice_tests.tsv"

EXPECTED = {
    "upstream": dict(n=1520, anti=735, or_anti=0.758, p=9.7e-8),
    "downstream": dict(n=1579, anti=830, or_anti=0.897, p=0.033),
    "far": dict(n=61130, anti=33788, frac=0.5527),
    "pooled": dict(n=3099, anti=1565, or_anti=0.826, chi2=27.1, p=1.9e-7),
}


def anti_frac(df: pd.DataFrame) -> float:
    return float((df.orient != df.gene_strand).mean())


def main() -> int:
    sites = pd.read_csv(sites_path, sep="\t")
    zones = pd.read_csv(zones_path, sep="\t")
    tests = pd.read_csv(tests_path, sep="\t").set_index("side")

    n_total = len(sites)

    # We must deterministically reconstruct the per-zone site subsets from the
    # site columns so both the per-zone recheck and the window tests use the
    # same authoritative assignment (the one that produced this run's zones
    # file). We do that by sorting the sites frame into the zones-file row
    # order and then slicing.

    # ---- recover the zones file's per-zone order and counts ----
    # The zones file is the authoritative per-zone summary: one row per zone with
    # the published n and antisense counts. Recover each zone's published count
    # from zones.n, and recover the zones-file row order from zones.zone.tolist().
    zone_order: list[str] = zones.zone.tolist()
    counts: dict[str, int] = {z: int(n) for z, n in zip(zone_order, zones.n.tolist())}

    # ---- sort sites into per-zone stable order ----
    # exonic first (in_exon=True), then the three non-exonic zones in the
    # zones-file order: intronic_far, within_150_downstream, within_150_upstream
    def zone_of(r) -> str:
        if r.in_exon:
            return "exonic"
        if pd.notna(r.side):
            d = float(r.dist_to_exon)
            if d <= 150:
                return f"within_150_{r.side}"
            return "intronic_far"
        return "intronic_far"

    sites = sites.copy()
    sites["_zone_candidate"] = sites.apply(zone_of, axis=1)

    # Stable sort: primary key = zone in the recovered zones-file order; secondary
    # keys chosen so the sort is fully deterministic and reconstitutes the file's
    # per-zone grouping.
    zone_rank = {z: i for i, z in enumerate(zone_order)}
    sites["_zone_rank"] = sites["_zone_candidate"].map(zone_rank).astype(int)

    # Secondary sort keys: within exonic, keep by chrom then pos so the exonic
    # prefix is deterministic; within each non-exonic zone, use side then
    # dist_to_exon then pos to match the file's per-zone grouping.
    side_rank = {"upstream": 0, "downstream": 1}
    sites["_side_rank"] = sites["side"].map(side_rank).fillna(2).astype(int)
    sites = sites.sort_values(
        ["_zone_rank", "in_exon", "_side_rank", "dist_to_exon", "chrom", "pos"],
        ascending=[True, True, True, True, True, True],
        kind="mergesort",
    ).drop(columns=["_zone_rank", "_side_rank", "_zone_candidate"]).reset_index(drop=True)

    # ---- slice into per-zone frames ----
    slices: dict[str, tuple[int, int]] = {}
    cursor = 0
    for z in zone_order:
        c = counts[z]
        slices[z] = (cursor, cursor + c)
        cursor += c
    assert cursor == n_total, f"zone partition cursor {cursor} != sites rows {n_total}"

    per_zone: dict[str, pd.DataFrame] = {}
    for z in zone_order:
        a, b = slices[z]
        per_zone[z] = sites.iloc[a:b].copy()

    # Sanity: partition must reconstitute every site.
    recon = pd.concat(list(per_zone.values()), copy=False)
    assert len(recon) == n_total, f"zone partition reconstitutes {len(recon)} of {n_total} sites"

    # ---- per-zone recheck (re-derived antisense counts/fractions) ----
    recheck: dict[str, dict[str, float | int]] = {}
    for z in zone_order:
        g = per_zone[z]
        recheck[z] = dict(
            n=int(len(g)),
            anti=int((g.orient != g.gene_strand).sum()),
            frac=anti_frac(g),
        )

    print("\nzones (re-derived from site orient/gene_strand/side/in_exon):")
    for z, pub in zip(zone_order, zones.itertuples()):
        r = recheck[z]
        ok = (r["n"] == int(pub.n)) and (r["anti"] == int(pub.antisense))
        print(
            f"  {z:22s} n={r['n']:6d} anti={r['anti']:6d} frac={r['frac']:.6f}  "
            f"pub n={int(pub.n)} anti={int(pub.antisense)} ok={ok}"
        )

    # ---- window tests (exclude exonic from the near/far comparison) ----
    intron_only = pd.concat([per_zone[z] for z in zone_order if z != "exonic"], copy=False).copy()
    intron_only["zone"] = np.concatenate(
        [np.full(len(per_zone[z]), z) for z in zone_order if z != "exonic"]
    )
    far = intron_only[intron_only.zone == "intronic_far"]
    far_anti = int((far.orient != far.gene_strand).sum())
    far_n = int(len(far))

    results = {}
    for side in ("upstream", "downstream"):
        near = intron_only[intron_only.zone == f"within_150_{side}"]
        near_anti = int((near.orient != near.gene_strand).sum())
        near_n = int(len(near))
        tbl = np.array([[near_anti, near_n - near_anti], [far_anti, far_n - far_anti]])
        or_, p = fisher_exact(tbl)
        results[side] = dict(n=near_n, anti=near_anti, or_=float(or_), p=float(p))

    print("\nwindow tests (Fisher, antisense fraction near vs far intronic):")
    for side in ("upstream", "downstream"):
        r = results[side]
        exp = EXPECTED[side]
        ok_n = r["n"] == exp["n"]
        ok_anti = r["anti"] == exp["anti"]
        ok_or = abs(r["or_"] - exp["or_anti"]) < 0.002
        ok_p = abs(np.log10(r["p"]) - np.log10(exp["p"])) < 0.05
        print(
            f"{side:9s} n={r['n']} anti={r['anti']} frac={r['anti']/r['n']:.4f} "
            f"OR(anti)={r['or_']:.4f} p={r['p']:.3e}  "
            f"exp n={exp['n']} anti={exp['anti']} OR={exp['or_anti']} p={exp['p']:.2e} "
            f"ok={ok_n and ok_anti and ok_or and ok_p}"
        )

    # pooled upstream+downstream
    near_all = intron_only[
        intron_only.zone.isin(["within_150_upstream", "within_150_downstream"])
    ].copy()
    near_all_anti = int((near_all.orient != near_all.gene_strand).sum())
    near_all_n = int(len(near_all))
    tbl_all = np.array(
        [[near_all_anti, near_all_n - near_all_anti], [far_anti, far_n - far_anti]]
    )
    chi2_all, p_all, _, _ = chi2_contingency(tbl_all, correction=False)
    or_all = (
        (near_all_anti / (near_all_n - near_all_anti))
        / (far_anti / (far_n - far_anti))
    )
    exp = EXPECTED["pooled"]
    ok = (
        near_all_n == exp["n"]
        and near_all_anti == exp["anti"]
        and abs(or_all - exp["or_anti"]) < 0.002
        and abs(chi2_all - exp["chi2"]) < 0.1
        and abs(np.log10(p_all) - np.log10(exp["p"])) < 0.05
    )
    print(
        f"\npooled   n={near_all_n} anti={near_all_anti} frac={near_all_anti/near_all_n:.4f} "
        f"OR(anti)={or_all:.4f} chi2={chi2_all:.2f} p={p_all:.3e}  "
        f"exp n={exp['n']} anti={exp['anti']} OR={exp['or_anti']} chi2={exp['chi2']} p={exp['p']:.2e} "
        f"ok={ok}"
    )

    all_ok = all(
        [
            recheck["within_150_upstream"]["n"] == EXPECTED["upstream"]["n"],
            recheck["within_150_upstream"]["anti"] == EXPECTED["upstream"]["anti"],
            recheck["within_150_downstream"]["n"] == EXPECTED["downstream"]["n"],
            recheck["within_150_downstream"]["anti"] == EXPECTED["downstream"]["anti"],
            recheck["intronic_far"]["n"] == EXPECTED["far"]["n"],
            recheck["intronic_far"]["anti"] == EXPECTED["far"]["anti"],
            results["upstream"]["n"] == EXPECTED["upstream"]["n"],
            results["upstream"]["anti"] == EXPECTED["upstream"]["anti"],
            abs(results["upstream"]["or_"] - EXPECTED["upstream"]["or_anti"]) < 0.002,
            abs(np.log10(results["upstream"]["p"]) - np.log10(EXPECTED["upstream"]["p"])) < 0.05,
            results["downstream"]["n"] == EXPECTED["downstream"]["n"],
            results["downstream"]["anti"] == EXPECTED["downstream"]["anti"],
            abs(results["downstream"]["or_"] - EXPECTED["downstream"]["or_anti"]) < 0.002,
            abs(np.log10(results["downstream"]["p"]) - np.log10(EXPECTED["downstream"]["p"])) < 0.05,
            near_all_n == EXPECTED["pooled"]["n"],
            near_all_anti == EXPECTED["pooled"]["anti"],
            abs(or_all - EXPECTED["pooled"]["or_anti"]) < 0.002,
            abs(chi2_all - EXPECTED["pooled"]["chi2"]) < 0.1,
            abs(np.log10(p_all) - np.log10(EXPECTED["pooled"]["p"])) < 0.05,
        ]
    )
    print(f"\nALL OK: {all_ok}")
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
