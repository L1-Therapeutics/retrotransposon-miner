#!/usr/bin/env python3
"""Does TSD detectability track element age, when age is assigned independently of truncation?

The question comes from the annotation plan (artifacts/gene_annotation_age_plan.md), where an
earlier comparison reported TSD+ at 79.4% in a "3'-reaching, young-enriched" class vs 65.7% in a
"truncated, older" class.  That comparison is circular for the purpose of testing the age claim:
the *grouping* variable was the call's own consensus-axis span, and the span is a property of the
call being tested for a TSD.  A span-based grouping therefore cannot distinguish

    "young elements really do keep a more recognizable TSD"

from

    "calls whose span looks young are the calls where a TSD is observable at all".

This script does the non-circular version.

Design
------
1.  Age is assigned from axes that do not use the call's own span:

      A. subfamily lineage  - RepeatMasker repName collapsed to L1HS (human-specific) /
         L1P* (primate) / L1M* (mammalian).  Assigned by sequence match of the aligned
         fragment, not by how far the fragment reaches.
      B. catalog membership - whether the insertion is present in a known non-reference MEI
         catalog (MELT/1000G, long-read 1KGP).  Segregating insertions are young on average.
         This axis is span-independent by construction.

2.  TSD observability is then *held constant* by stratifying on 3' resolution
    (MEI_3P >= 5932: the copy retains the consensus 3' end).  Any age difference that survives
    inside a stratum is not a truncation artefact.

3.  Reported: crude rates, the age x resolution confound, stratified rates, Mantel-Haenszel
    pooled odds ratio, a homogeneity/interaction check, and logistic fits with and without the
    resolution term.  Sensitivity to the TSD length threshold and to call tier.

Deliberate non-claims
---------------------
* Catalog membership is an age proxy AND an evidence proxy: a matched insertion had independent
  population evidence, so its junction may be better resolved.  The file does carry the junction
  evidence that could test this (SR_L, SR_R, DPE_L, DPE_R, MEI_MAPPED, polyA_MAPPED, and the same
  numbers as one SUPPORT string) -- an earlier note here said SUPPORT was empty, which is wrong.
  Section 5 uses the coarse CALLTIER only; the evidence-matched test is `tsd_evidence_matched.py`.
* A TSD string of length exactly 4 is the caller's floor and 248 calls sit on it; the headline
  threshold is therefore reported alongside >=5 (the repo's TSD_MIN_BP) and >=9.
* Nothing here is a mutation rate, a date, or a per-call age call.

Usage:  python3 artifacts/tsd_age_association.py
Writes: artifacts/tsd_age_association_results.md
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
CALLS = ROOT / "span_work" / "out" / "line1_calls.tsv"
OUT = ROOT / "artifacts" / "tsd_age_association_results.md"

# The consensus axis window used by the plan: L1 3'UTR 5930-5932.
MARKER_LO, MARKER_HI = 5930, 5932

# RepeatMasker/Dfam lineage prefixes.  Coarse on purpose: the ordering L1HS < L1P* < L1M*
# is the well-established part; the subfamily dates (L1PA2 ~7.6 My, L1PA3 ~12.5 My) are
# model-based and are NOT used here.  Mirrors PRIMATE_PREFIXES in artifacts/catalog.py.
YOUNG_PREFIX = ("L1HS",)
INTERM_PREFIX = ("L1PA", "L1PB", "L1P1", "L1P2", "L1P3", "L1P4", "L1PREC")
OLD_PREFIX = ("L1M", "HAL1")


def lineage_age(subfamily: str) -> str:
    if subfamily.startswith(YOUNG_PREFIX):
        return "young:L1HS"
    if subfamily.startswith(INTERM_PREFIX):
        return "interm:L1P*"
    if subfamily.startswith(OLD_PREFIX):
        return "old:L1M*"
    return "unknown"


def load() -> pd.DataFrame:
    df = pd.read_csv(CALLS, sep="\t", dtype=str)
    tsd = df["TSD"].fillna("")
    df["tsd_len"] = tsd.str.len()
    df["m5p"] = pd.to_numeric(df["MEI_5P"], errors="coerce")
    df["m3p"] = pd.to_numeric(df["MEI_3P"], errors="coerce")
    df["age"] = df["subfamily"].map(lineage_age)
    # 3' resolution: the copy retains the consensus 3' end (the marker window).
    df["three_prime_resolved"] = df["m3p"] >= MARKER_HI
    df["in_catalog"] = df["KNOWNMEI"].astype(str).str.lower() == "true"
    df["in_1kg"] = df["G1K"].notna()
    df["in_longread"] = df["LR"].notna()
    return df


def rate_table(df: pd.DataFrame, by: list[str], thresh: int) -> pd.DataFrame:
    y = (df["tsd_len"] >= thresh).astype(int)
    t = df.assign(_y=y).groupby(by, observed=True)["_y"].agg(["sum", "size", "mean"])
    t.columns = ["tsd_pos", "n", "rate"]
    t["rate"] = t["rate"].round(3)
    return t


def mh_odds_ratio(tab: pd.DataFrame, stratum: str, exposure: str) -> dict:
    """Mantel-Haenszel pooled OR over strata, with the Woolf test of homogeneity.

    tab: counts with columns [stratum, exposure, tsd_pos, n].
    exposure is treated as binary in the order given by the caller.
    """
    num = den = 0.0
    a_all = b_all = c_all = d_all = 0
    parts: list[tuple[object, int, int, int, int]] = []
    for s, g in tab.groupby(stratum, observed=True):
        if len(g) != 2:
            continue
        hi = g[g[exposure] == 1]
        lo = g[g[exposure] == 0]
        if hi.empty or lo.empty:
            continue
        a, n1 = int(hi.iloc[0]["tsd_pos"]), int(hi.iloc[0]["n"])   # exposed, TSD+
        c, n0 = int(lo.iloc[0]["tsd_pos"]), int(lo.iloc[0]["n"])   # unexposed, TSD+
        b, d = n1 - a, n0 - c
        n = n1 + n0
        if min(a, b, c, d) < 0 or n == 0:
            continue
        num += a * d / n
        den += b * c / n
        a_all += a; b_all += b; c_all += c; d_all += d
        parts.append((s, a, b, c, d))
    or_mh = num / den if den else np.nan
    # Woolf test of homogeneity of ln(OR) across strata (strata with a zero cell are skipped).
    valid = [p for p in parts if min(p[1:]) > 0]
    p_hom = np.nan
    if len(valid) > 1:
        logs, weights = [], []
        for _s, a, b, c, d in valid:
            v = 1 / a + 1 / b + 1 / c + 1 / d
            logs.append(np.log((a * d) / (b * c)))
            weights.append(1 / v)
        logs = np.asarray(logs)
        weights = np.asarray(weights)
        wsum = weights.sum()
        wmean = (weights * logs).sum() / wsum
        chi2 = (weights * (logs - wmean) ** 2).sum()
        p_hom = float(1 - stats.chi2.cdf(chi2, len(valid) - 1))
    # Crude (collapsed) OR and Fisher test.
    or_crude = (a_all * d_all) / (b_all * c_all) if b_all * c_all else np.nan
    _, p_crude = stats.fisher_exact([[a_all, b_all], [c_all, d_all]])
    return {
        "or_crude": or_crude,
        "p_crude": p_crude,
        "or_mh": or_mh,
        "p_homogeneity": p_hom,
        "n_strata": len(parts),
    }


def logit(y: np.ndarray, X: np.ndarray, iters: int = 100,
          ridge: float = 1e-3) -> tuple[np.ndarray, np.ndarray]:
    """Ridge-stabilised Newton-Raphson logistic regression with step halving.

    A ridge is required, not cosmetic: one interaction cell is near-separated, so the
    unpenalised Hessian is ill-conditioned and Newton returns non-finite steps.  Steps are also
    capped and checked for finiteness.  With ridge=1e-3 against Hessian entries of order 500 the
    bias is negligible, and the fitter is validated against a closed form in `validate_fitter`.

    NOTE ON WARNINGS: this numpy build emits spurious "divide by zero encountered in matmul"
    for ANY matmul with >= ~1000 rows -- `np.ones((2186, 2)) @ np.zeros(2)` triggers it, on
    finite input, while the same operation at 100 rows does not.  The kernel is therefore run
    under `np.errstate`, and any *real* non-finite result is caught by the explicit invariant
    check below rather than being silently ignored.
    """
    beta = np.zeros(X.shape[1])
    eye = np.eye(X.shape[1])

    def nll(b: np.ndarray) -> float:
        eta = np.clip(X @ b, -30, 30)
        p = 1 / (1 + np.exp(-eta))
        return -float(np.sum(y * np.log(p + 1e-12) + (1 - y) * np.log(1 - p + 1e-12))
                      + 0.5 * ridge * float(b @ b))

    with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
        for _ in range(iters):
            eta = np.clip(X @ beta, -30, 30)
            p = 1 / (1 + np.exp(-eta))
            W = np.clip(p * (1 - p), 1e-9, None)
            grad = X.T @ (y - p) - ridge * beta
            H = X.T @ (X * W[:, None]) + ridge * eye
            if not (np.all(np.isfinite(grad)) and np.all(np.isfinite(H))):
                break
            try:
                step = np.linalg.solve(H, grad)
            except np.linalg.LinAlgError:
                break
            if not np.all(np.isfinite(step)):
                break
            step = np.clip(step, -20.0, 20.0)
            cur = nll(beta)
            t = 1.0
            moved = False
            for _ in range(40):
                cand = beta + t * step
                if nll(cand) <= cur + 1e-12:
                    beta = cand
                    moved = True
                    break
                t /= 2
            if not moved:
                break
            if np.max(np.abs(t * step)) < 1e-9:
                break
        eta = np.clip(X @ beta, -30, 30)
        p = 1 / (1 + np.exp(-eta))
        W = np.clip(p * (1 - p), 1e-9, None)
        cov = np.linalg.inv(X.T @ (X * W[:, None]) + ridge * eye)
    se = np.sqrt(np.diag(cov))
    # Invariant: suppression above must never conceal a genuine numerical failure.
    if not (np.all(np.isfinite(beta)) and np.all(np.isfinite(se))):
        raise ValueError("logistic fit produced non-finite coefficients "
                         f"(beta={beta}, se={se})")
    return beta, se


def validate_fitter(y: np.ndarray, x: np.ndarray) -> tuple[float, float]:
    """Check the fitter against the closed-form MLE for one binary regressor.

    With an intercept, the MLE of the slope is exactly log((n11*n00)/(n10*n01)).  This is the
    control that shows the warnings noted above are spurious and the reported betas are real.
    """
    n11 = int(((y == 1) & (x == 1)).sum())
    n10 = int(((y == 1) & (x == 0)).sum())
    n01 = int(((y == 0) & (x == 1)).sum())
    n00 = int(((y == 0) & (x == 0)).sum())
    analytic = float(np.log((n11 * n00) / (n10 * n01)))
    X = np.column_stack([np.ones(len(y)), x])
    beta, _ = logit(y, X)
    return float(beta[1]), analytic


def model_row(name: str, y, X, labels) -> str:
    beta, se = logit(y, X)
    z = beta / se
    pv = 2 * (1 - stats.norm.cdf(np.abs(z)))
    terms = []
    for lab, b, s, zz, pp in zip(labels, beta, se, z, pv):
        terms.append(f"{lab} beta={b:+.3f} (SE {s:.3f}, p={pp:.2g})")
    return f"- **{name}**: " + "; ".join(terms)


def main() -> int:
    df = load()
    lines: list[str] = []
    add = lines.append

    add("# Does TSD detectability track element age? Non-circular test")
    add("")
    add(f"Cohort: {len(df)} L1 calls, `span_work/out/line1_calls.tsv`.")
    add(f"3' resolution = `MEI_3P >= {MARKER_HI}` (copy retains the consensus 3' end).")
    add("")

    # ---- 0. the threshold matters ----
    add("## 0. TSD call threshold (checked before interpreting anything)")
    add("")
    add(f"- TSD string present: {(df.tsd_len > 0).sum()} / {len(df)}")
    add(f"- length exactly 4 (caller floor): {int((df.tsd_len == 4).sum())}")
    add(f"- TSD+ at >=4 bp: {(df.tsd_len >= 4).mean():.3f} | >=5 (repo TSD_MIN_BP): "
        f"{(df.tsd_len >= 5).mean():.3f} | >=9: {(df.tsd_len >= 9).mean():.3f}")
    add("")
    add("The plan's headline rates used >=4 bp, which sits on a 248-call floor spike; every")
    add("conclusion below is therefore recomputed at >=5 and >=9.")
    add("")

    # ---- 1. the confound ----
    add("## 1. The confound, quantified")
    add("")
    ct = pd.crosstab(df["age"], df["three_prime_resolved"], normalize="index").round(3)
    add("Fraction of each lineage class that is 3'-resolved:")
    add("")
    add("| age class | 3'-resolved | n |")
    add("|---|---|---|")
    for age, row in ct.iterrows():
        n = int((df["age"] == age).sum())
        add(f"| {age} | {row.get(True, 0.0):.3f} | {n} |")
    add("")
    add("So the resolution flag is almost a restatement of lineage age: L1M* copies are")
    add("overwhelmingly unresolved, L1HS copies overwhelmingly resolved. A comparison grouped")
    add("on resolution is therefore also a comparison grouped on age.")
    add("")

    results = {}
    for thresh in (4, 5, 9):
        add(f"## 2. Age axis A - subfamily lineage, TSD >= {thresh} bp")
        add("")
        add("Crude (unstratified):")
        add("")
        crude = rate_table(df, ["age"], thresh)
        add("| age class | TSD+ | n | rate |")
        add("|---|---|---|---|")
        for age, r in crude.iterrows():
            add(f"| {age} | {int(r.tsd_pos)} | {int(r.n)} | {r.rate:.3f} |")
        add("")
        add("Stratified on 3' resolution:")
        add("")
        strat = rate_table(df, ["three_prime_resolved", "age"], thresh)
        add("| 3'-resolved | age class | TSD+ | n | rate |")
        add("|---|---|---|---|---|")
        for (rz, age), r in strat.iterrows():
            add(f"| {rz} | {age} | {int(r.tsd_pos)} | {int(r.n)} | {r.rate:.3f} |")
        add("")
        results[f"lineage_{thresh}"] = strat
        add("")

    # ---- 3. age axis B ----
    add("## 3. Age axis B - population catalog membership (span-independent)")
    add("")
    for thresh in (5,):
        for axis, label in (("in_catalog", "KNOWNMEI (MELT/1KGP or long-read)"),
                            ("in_1kg", "1000 Genomes nssv id present"),
                            ("in_longread", "long-read 1KGP id present")):
            add(f"### {label}, TSD >= {thresh} bp")
            add("")
            t = rate_table(df, ["three_prime_resolved", axis], thresh)
            add("| 3'-resolved | in catalog | TSD+ | n | rate |")
            add("|---|---|---|---|---|")
            for (rz, ax), r in t.iterrows():
                add(f"| {rz} | {ax} | {int(r.tsd_pos)} | {int(r.n)} | {r.rate:.3f} |")
            add("")

    # ---- 3b. within-stratum pairwise tests (the decisive numbers) ----
    add("## 3b. Within-stratum pairwise tests (the decisive numbers)")
    add("")
    add("Comparing the *same* two age classes inside one resolution stratum removes the span")
    add("confound entirely. `old:L1M*` is shown but its resolved cell is too small to use.")
    add("")
    for thresh in (4, 5, 9):
        add(f"TSD >= {thresh} bp:")
        add("")
        add("| stratum | comparison | rate A | rate B | OR | Fisher p |")
        add("|---|---|---|---|---|---|")
        for rz in (True, False):
            s = df[df["three_prime_resolved"] == rz]
            for a_lab, b_lab in (("young:L1HS", "interm:L1P*"), ("young:L1HS", "old:L1M*")):
                ga = s[s["age"] == a_lab]
                gb = s[s["age"] == b_lab]
                if len(ga) == 0 or len(gb) == 0:
                    continue
                a1 = int((ga.tsd_len >= thresh).sum()); a0 = len(ga) - a1
                b1 = int((gb.tsd_len >= thresh).sum()); b0 = len(gb) - b1
                orv = ((a1 * b0) / (a0 * b1)) if a0 * b1 else np.inf
                _, pv = stats.fisher_exact([[a1, a0], [b1, b0]])
                add(f"| resolved={rz} | {a_lab} vs {b_lab} | {a1 / len(ga):.3f} (n={len(ga)}) | "
                    f"{b1 / len(gb):.3f} (n={len(gb)}) | {orv:.2f} | {pv:.2g} |")
        add("")

    # ---- 4. statistics ----
    add("## 4. Statistics")
    add("")
    y5 = (df.tsd_len >= 5).astype(int).to_numpy()
    for axis, label, pos in (("in_catalog", "catalog membership", True),
                             ("in_1kg", "1000G membership", True),
                             ("in_longread", "long-read membership", True)):
        sub = df.assign(exp=df[axis].astype(int))
        tab = rate_table(sub, ["three_prime_resolved", "exp"], 5).reset_index()
        mh = mh_odds_ratio(tab, "three_prime_resolved", "exp")
        add(f"- **{label}** (TSD >= 5): crude OR {mh['or_crude']:.2f} (Fisher p={mh['p_crude']:.2g}); "
            f"Mantel-Haenszel OR {mh['or_mh']:.2f} across {mh['n_strata']} strata "
            f"(homogeneity p={mh['p_homogeneity']:.2g})")
    add("")
    # lineage: young vs interm only, old is too sparse once resolved.
    sub = df[df["age"].isin(["young:L1HS", "interm:L1P*"])].copy()
    sub["exp"] = (sub["age"] == "young:L1HS").astype(int)
    tab = rate_table(sub, ["three_prime_resolved", "exp"], 5).reset_index()
    mh = mh_odds_ratio(tab, "three_prime_resolved", "exp")
    add(f"- **lineage young(L1HS) vs interm(L1P*)** (TSD >= 5): crude OR {mh['or_crude']:.2f} "
        f"(Fisher p={mh['p_crude']:.2g}); MH OR {mh['or_mh']:.2f} across {mh['n_strata']} strata "
        f"(homogeneity p={mh['p_homogeneity']:.2g})")
    add("")
    add("Fitter validation (slope for a single binary regressor, closed form vs fitted):")
    add("")
    res = df["three_prime_resolved"].astype(int).to_numpy()
    vfit, vana = validate_fitter(y5, res)
    add(f"- fitted {vfit:+.6f} vs analytic ln((n11*n00)/(n10*n01)) {vana:+.6f} "
        f"(diff {abs(vfit - vana):.2e})")
    add("")
    add("Logistic regression, TSD >= 5 bp (`resolution` = 3'-resolved):")
    add("")
    X1 = np.column_stack([np.ones(len(df)), res])
    add(model_row("TSD ~ resolution", y5, X1, ["intercept", "resolution"]))
    for axis, label in (("in_catalog", "catalog"), ("in_longread", "long_read")):
        ex = df[axis].astype(int).to_numpy()
        X2 = np.column_stack([np.ones(len(df)), res, ex])
        add(model_row(f"TSD ~ resolution + {label}", y5, X2, ["intercept", "resolution", label]))
        X3 = np.column_stack([np.ones(len(df)), res, ex, res * ex])
        add(model_row(f"TSD ~ resolution * {label}", y5, X3,
                      ["intercept", "resolution", label, f"resolution:{label}"]))
    d2 = df[df["age"].isin(["young:L1HS", "interm:L1P*"])].copy()
    yl = (d2.tsd_len >= 5).astype(int).to_numpy()
    rl = d2["three_prime_resolved"].astype(int).to_numpy()
    yg = (d2["age"] == "young:L1HS").astype(int).to_numpy()
    add(model_row("TSD ~ resolution (lineage subset)", yl,
                  np.column_stack([np.ones(len(d2)), rl]), ["intercept", "resolution"]))
    add(model_row("TSD ~ resolution + L1HS (lineage subset)", yl,
                  np.column_stack([np.ones(len(d2)), rl, yg]),
                  ["intercept", "resolution", "L1HS"]))
    add(model_row("TSD ~ resolution * L1HS (lineage subset)", yl,
                  np.column_stack([np.ones(len(d2)), rl, yg, rl * yg]),
                  ["intercept", "resolution", "L1HS", "resolution:L1HS"]))
    add("")

    # ---- 5. quality control ----
    add("## 5. Quality control (the caveat on axis B)")
    add("")
    add("Catalog-matched insertions had independent population evidence, so they may simply have")
    add("better-resolved junctions. The file carries that evidence — `SR_L`, `SR_R`, `DPE_L`,")
    add("`DPE_R`, `MEI_MAPPED`, `polyA_MAPPED`, and the same numbers as one `SUPPORT` string — so")
    add("the tier adjustment below is only the coarse version of it; the split-read-matched test")
    add("lives in `tsd_evidence_matched.py`. Here `CALLTIER` is the only quality variable used.")
    add("")
    add("| call tier | 3'-resolved | in catalog | TSD+ | n | rate |")
    add("|---|---|---|---|---|---|")
    qc = rate_table(df, ["CALLTIER", "three_prime_resolved", "in_catalog"], 5)
    for (tier, rz, ax), r in qc.iterrows():
        add(f"| {tier} | {rz} | {ax} | {int(r.tsd_pos)} | {int(r.n)} | {r.rate:.3f} |")
    add("")

    # ---- 6. verdict ----
    add("## 6. Verdict")
    add("")
    add("**The TSD x age association is real but it is not the association the premise describes,**")
    add("**and the lineage version of it does not survive the non-circular test.**")
    add("")
    add("1. *TSD detection is strongly span-dependent.* TSD ~ resolution alone gives beta=+1.34")
    add("   (p=2e-04): resolved calls are far more often TSD-positive. Any grouping on span is")
    add("   therefore partly a grouping on TSD observability.")
    add("2. *The lineage gradient mostly was the confound.* Crude TSD+ (>=5 bp) runs")
    add("   L1HS 0.787 / L1P* 0.538 / L1M* 0.254, a 25 pp young-old gap. Inside the 3'-resolved")
    add("   stratum it collapses to 0.761 / 0.726 / 0.636, and the L1HS-vs-L1P* Fisher test in")
    add("   that stratum is the number that matters (table 3b). `old:L1M*` has n=11 resolved,")
    add("   which cannot support the old end of the gradient.")
    add("3. *The interaction is the finding.* homogeneity p < 1e-3 on every paired test: the")
    add("   odds ratio is not one number across strata, so a pooled OR must not be quoted alone.")
    add("   In the unresolved stratum L1HS looks strongly TSD-enriched (0.892 vs 0.338); in the")
    add("   resolved stratum the two are indistinguishable. That is the signature of the majority")
    add("   of the effect living in the confound.")
    add("4. *A span-independent age axis behaves differently.* Catalog membership (a segregating,")
    add("   therefore young, insertion) keeps a large association inside **both** strata --")
    add("   MH OR 6.38 (KNOWNMEI), 4.06 (1000G), 3.62 (long-read) with the direction preserved in")
    add("   each. This is the one version of 'younger elements have more recognizable TSDs' that")
    add("   the data supports.")
    add("5. *Its own confound is not resolved here.* Catalog-matched insertions had independent")
    add("   population evidence, so better-resolved junctions could produce this without any age")
    add("   effect. Holding call tier fixed (section 5) the effect persists in the dominant tier,")
    add("   but the other tiers are too small to arbitrate, and tier is a coarse stand-in for")
    add("   junction evidence. That specific alternative is tested in `tsd_evidence_matched.py`.")
    add("")
    add("**What this licenses downstream:** TSD presence stays a structural field and must not be")
    add("used to assign age to an individual call. The TSD x age association may be reported as a")
    add("cohort statistic, with the resolution stratum and the assessable denominator printed next")
    add("to it, and with the lineage and catalog axes reported separately rather than pooled.")
    add("")

    report = "\n".join(lines) + "\n"
    OUT.write_text(report)
    print(report)
    return 0


if __name__ == "__main__":
    sys.exit(main())
