#!/usr/bin/env python3
"""Is the surviving catalog-membership TSD effect just junction evidence quality?

Follow-up to `tsd_age_association.py`.  That analysis left one live alternative: catalog
membership (KNOWNMEI / 1000G / long-read) kept a large TSD association inside both 3'
strata (MH OR 6.38 / 4.06 / 3.62), and its own caveat was that a catalog-matched insertion had
independent population evidence, so its *junction may simply be better resolved* than a
singleton call's.  If that is the whole story, the catalog axis is an evidence-quality proxy
and not an age axis at all.  This script tries to kill the effect with evidence.

The file carries the evidence as first-class numeric columns -- the note in
`tsd_age_association.py` that "SUPPORT is empty in this file" is wrong, and this script depends
on that correction:

    SR_L, SR_R          split-read support flanking the left / right breakpoint
    DPE_L, DPE_R        discordant paired-end support
    MEI_MAPPED          element bases/reads mapped
    polyA_MAPPED        poly(A) bases/reads mapped
    L1TXGOLDSCORE, INSERTIONSCORE, PEAKDEPTHZ, SPLITCLUSTERZ, CALLTIER
    (SUPPORT is the same numbers as one `key:value|...` string.)

Primary evidence score: `min_sr = min(SR_L, SR_R)` -- both breakpoints observed by split
reads.  A call a population catalog resolved should have min_sr >= 1 far more often than a
singleton, which is exactly the proposed mechanism.

Six instruments, because one matched test is easy to fool:

  T1  the evidence imbalance the matching has to remove
  T2  dose-response: does evidence quality itself predict TSD presence?
  T3  evidence-matched pooling: strata = (resolution x min_sr bucket)
  T4  continuous adjustment: logit(TSD ~ resolution + catalog [+ evidence terms])
  T5  the evidence-saturated subset: min_sr >= 1, >= 2, sr_total >= 5
  T6  symmetric control: the same machinery applied to the lineage axis, to show the matching
      does not mechanically destroy every association

Pooling note: the point estimate is Mantel-Haenszel (reused from `tsd_age_association.py`);
the confidence interval is Woolf's, computed from the same strata, so the interval is
self-consistent rather than borrowed from a different estimator.  Strata with a zero cell are
excluded from the interval and counted in the table, never silently dropped from the count.

Verdict text is selected by the computed numbers, not written in advance; the script prints
which branch its own instruments chose.

Usage:  python3 artifacts/tsd_evidence_matched.py                        # LINE1, from the TSV
        python3 artifacts/tsd_evidence_matched.py --callsets nested_analysis/callsets
                                                                        # every family, from VCFs
Writes: artifacts/tsd_evidence_matched_results.md                       # LINE1 route
        artifacts/tsd_evidence_matched_families_results.md              # --callsets route
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import tsd_age_association as base  # noqa: E402  (validated MH / rate_table / logit helpers)

ROOT = Path(__file__).resolve().parent.parent
CALLS = ROOT / "span_work" / "out" / "line1_calls.tsv"
OUT = ROOT / "artifacts" / "tsd_evidence_matched_results.md"
OUT_FAMILIES = ROOT / "artifacts" / "tsd_evidence_matched_families_results.md"

MARKER_HI = 5932          # 3' resolution: the copy retains the consensus 3' end
TSD_MIN_BP = 5            # repo TSD_MIN_BP; 4 is the caller floor spike
EVIDENCE_COLS = ["SR_L", "SR_R", "DPE_L", "DPE_R", "MEI_MAPPED", "polyA_MAPPED"]
BUCKETS = ("sr0", "sr1", "sr2plus")


def load() -> pd.DataFrame:
    df = pd.read_csv(CALLS, sep="\t", dtype=str)
    df = df[df["MEIFAMILY"].astype(str) == "LINE1"].copy()
    df["tsd_len"] = df["TSD"].fillna("").str.len()
    df["three_prime_resolved"] = pd.to_numeric(df["MEI_3P"], errors="coerce") >= MARKER_HI
    df["age"] = df["subfamily"].map(base.lineage_age)
    for c in EVIDENCE_COLS:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).astype(int)
    df["min_sr"] = np.minimum(df["SR_L"], df["SR_R"])
    df["sr_total"] = df["SR_L"] + df["SR_R"]
    df["min_sr_bucket"] = np.where(df["min_sr"] == 0, "sr0",
                                   np.where(df["min_sr"] == 1, "sr1", "sr2plus"))
    df["in_catalog"] = df["KNOWNMEI"].astype(str).str.lower() == "true"
    df["in_1kg"] = df["G1K"].notna()
    df["in_longread"] = df["LR"].notna()
    return df


def strata_pairs(df: pd.DataFrame, exposure: str, thresh: int,
                 cols: list[str]) -> tuple[list[tuple[int, int, int, int]], list[dict]]:
    """(a, b, c, d) per stratum = (exposed TSD+, exposed TSD-, unexposed TSD+, unexposed TSD-)."""
    pairs: list[tuple[int, int, int, int]] = []
    rows: list[dict] = []
    for key, g in df.groupby(cols, observed=True):
        key = key if isinstance(key, tuple) else (key,)
        e = g[g[exposure].astype(int) == 1]
        u = g[g[exposure].astype(int) == 0]
        if e.empty or u.empty:
            continue
        a = int((e["tsd_len"] >= thresh).sum()); b = len(e) - a
        c = int((u["tsd_len"] >= thresh).sum()); d = len(u) - c
        _, p = stats.fisher_exact([[a, b], [c, d]])
        pairs.append((a, b, c, d))
        rows.append({"stratum": " | ".join(str(x) for x in key),
                     "exposed": f"{a}/{a + b} ({a / max(a + b, 1):.3f})",
                     "unexposed": f"{c}/{c + d} ({c / max(c + d, 1):.3f})",
                     "or": ((a * d) / (b * c)) if b * c else np.inf, "p": p})
    return pairs, rows


def woolf(pairs: list[tuple[int, int, int, int]]) -> tuple[float, float, float, int]:
    """Woolf pooled OR with the matching Woolf interval; zero-cell strata are excluded.

    Var(ln OR_Woolf) = 1 / sum(w_i), with w_i = 1 / (1/a_i + 1/b_i + 1/c_i + 1/d_i).
    For a single stratum this reduces to the textbook 1/a + 1/b + 1/c + 1/d.
    """
    num = den = 0.0
    usable = 0
    for a, b, c, d in pairs:
        if min(a, b, c, d) <= 0:
            continue
        w = 1.0 / (1 / a + 1 / b + 1 / c + 1 / d)
        num += w * np.log((a * d) / (b * c))
        den += w
        usable += 1
    if usable == 0 or den <= 0:
        return float("nan"), float("nan"), float("nan"), 0
    or_w = float(np.exp(num / den))
    se = float(np.sqrt(1.0 / den))
    return or_w, float(np.exp(np.log(or_w) - 1.96 * se)), float(np.exp(np.log(or_w) + 1.96 * se)), usable


def rows_to_tab(df: pd.DataFrame, exposure: str, thresh: int, cols: list[str]) -> pd.DataFrame:
    """[(stratum, exp, tsd_pos, n)] in the shape `tsd_age_association.mh_odds_ratio` expects."""
    d = df.assign(exp=df[exposure].astype(int),
                  stratum=df[cols].astype(str).agg(" | ".join, axis=1),
                  y=(df["tsd_len"] >= thresh).astype(int))
    t = d.groupby(["stratum", "exp"], observed=True)["y"].agg(["sum", "size"]).reset_index()
    t.columns = ["stratum", "exp", "tsd_pos", "n"]
    return t


def pooled(df: pd.DataFrame, exposure: str, thresh: int, cols: list[str]) -> dict:
    """MH point estimate + Woolf interval over the same strata, with a direction tally."""
    pairs, rows = strata_pairs(df, exposure, thresh, cols)
    if not pairs:
        return {}
    or_w, lo, hi, usable = woolf(pairs)
    or_mh = base.mh_odds_ratio(rows_to_tab(df, exposure, thresh, cols),
                               "stratum", "exp").get("or_mh", float("nan"))
    same = sum(1 for a, b, c, d in pairs
               if (a / max(a + b, 1)) > (c / max(c + d, 1)))
    return {"or_mh": or_mh, "or_w": or_w, "lo": lo, "hi": hi, "usable": usable,
            "rows": rows, "n_strata": len(pairs), "direction_same": same,
            "zero_cell_strata": sum(1 for p in pairs if min(p) == 0)}


# --------------------------------------------------------------------------- #
# family-wide route: LINE1 is not special, so the same design runs per family   #
# --------------------------------------------------------------------------- #

# Family-specific lineage prefixes.  LINE1 keeps the labels the parent report already
# published; ALU and SVA get the analogous coarse split.  The lineages are NOT
# comparable across families -- "young:AluY" and "young:L1HS" are two different
# claims, and T6 is a control on each family, not a comparison between them.
FAMILY_LINEAGE: dict[str, tuple[tuple[tuple[str, ...], str], ...]] = {
    "LINE1": (
        (("L1HS",), "young:L1HS"),
        (("L1PA", "L1PB", "L1P1", "L1P2", "L1P3", "L1P4", "L1PREC"), "interm:L1P*"),
        (("L1M", "HAL1"), "old:L1M*"),
    ),
    "ALU": (
        (("AluY",), "young:AluY"),
        (("AluS",), "interm:AluS"),
        (("AluJ",), "old:AluJ"),
    ),
    "SVA": (
        (("SVA_F",), "young:SVA_F"),
        (("SVA_E", "SVA_D"), "interm:SVA_D/E"),
        (("SVA_A", "SVA_B", "SVA_C"), "old:SVA_A/B/C"),
    ),
}


def lineage_age_family(family: str, subfamily: str) -> str:
    """Coarse lineage class for one family's subfamily name."""
    for prefixes, label in FAMILY_LINEAGE.get(family, ()):
        if subfamily.startswith(prefixes):
            return label
    return "unknown"


def load_family_wide(callsets: Path) -> pd.DataFrame:
    """Rebuild the family-wide cohort exactly as span_work/line1_span_support.py builds it.

    Reusing that loader is deliberate: `span_work/out/line1_calls.tsv` is nothing but its
    LINE1 filter, so the sentinel rule and the SUPPORT parsing stay identical to the
    published cohort.  The LINE1 subset of what this returns must still number 2186 calls;
    `family_wide_main` asserts that before reporting anything.
    """
    sys.path.insert(0, str(ROOT / "span_work"))
    import line1_span_support as spansup  # cohort loader, reused not reimplemented

    usable = spansup.annotate(spansup.load_cohort(callsets))
    df = usable.copy()
    df["family"] = df["MEIFAMILY"].astype(str)
    df["tsd_len"] = df["TSD"].fillna("").str.len()
    df["three_prime_resolved"] = pd.to_numeric(df["MEI_3P"], errors="coerce") >= MARKER_HI
    df["age"] = [lineage_age_family(f, s) for f, s in
                 zip(df["family"], df["subfamily"].fillna("").astype(str))]
    df["in_catalog"] = df["KNOWNMEI"].astype(str).str.lower() == "true"
    df["in_1kg"] = df["G1K"].notna()
    df["in_longread"] = df["LR"].notna()
    for c in EVIDENCE_COLS:
        df[c] = pd.to_numeric(df.get(c), errors="coerce").fillna(0).astype(int)
    df["min_sr"] = np.minimum(df["SR_L"], df["SR_R"])
    df["sr_total"] = df["SR_L"] + df["SR_R"]
    df["min_sr_bucket"] = np.where(df["min_sr"] == 0, "sr0",
                                   np.where(df["min_sr"] == 1, "sr1", "sr2plus"))
    return df


def _sat_subset(df: pd.DataFrame, sel) -> dict | None:
    s = df[sel]
    e, u = s[s["in_catalog"]], s[~s["in_catalog"]]
    if e.empty or u.empty:
        return None
    a = int((e["tsd_len"] >= TSD_MIN_BP).sum()); b = len(e) - a
    c = int((u["tsd_len"] >= TSD_MIN_BP).sum()); d = len(u) - c
    _, p = stats.fisher_exact([[a, b], [c, d]])
    return {"or": (a * d) / (b * c) if b * c else np.inf, "p": p,
            "a": a, "b": b, "c": c, "d": d}


def family_summary(df: pd.DataFrame, fam: str) -> dict:
    """The instruments, on one family.  Numbers only; `family_report` writes the text."""
    res_cols, mat_cols = ["three_prime_resolved"], ["three_prime_resolved", "min_sr_bucket"]
    young, interm = FAMILY_LINEAGE[fam][0][1], FAMILY_LINEAGE[fam][1][1]
    un = pooled(df, "in_catalog", TSD_MIN_BP, res_cols) or {}
    mt = pooled(df, "in_catalog", TSD_MIN_BP, mat_cols) or {}
    tier = pooled(df, "in_catalog", TSD_MIN_BP,
                  ["three_prime_resolved", "min_sr_bucket", "CALLTIER"]) or {}
    sat = {lab: v for lab, v in (("min_sr >= 1", _sat_subset(df, df["min_sr"] >= 1)),
                                 ("min_sr >= 2", _sat_subset(df, df["min_sr"] >= 2)),
                                 ("sr_total >= 5", _sat_subset(df, df["sr_total"] >= 5)))
           if v is not None}
    lin = df[df["age"].isin([young, interm])].copy()
    lin["is_young"] = (lin["age"] == young).astype(int)
    lin_ok = lin["is_young"].nunique() == 2
    un_lin = (pooled(lin, "is_young", TSD_MIN_BP, res_cols) or {}) if lin_ok else {}
    mt_lin = (pooled(lin, "is_young", TSD_MIN_BP, mat_cols) or {}) if lin_ok else {}
    # T4: how far the catalog coefficient falls once junction evidence is in the model.
    res = df["three_prime_resolved"].astype(int).to_numpy(dtype=float)
    cat = df["in_catalog"].astype(int).to_numpy(dtype=float)
    one = np.ones(len(df))
    y5 = (df["tsd_len"] >= TSD_MIN_BP).astype(int).to_numpy()
    sr_l = np.log1p(df["SR_L"].to_numpy(dtype=float))
    sr_r = np.log1p(df["SR_R"].to_numpy(dtype=float))
    # A small family can be separated, and base.logit raises on non-finite coefficients.
    # Report that family's adjustment as unavailable rather than aborting the whole run.
    try:
        b_adj = base.logit(y5, np.column_stack([one, res, cat]))[0][2]
        b_full = base.logit(y5, np.column_stack([one, res, cat, sr_l, sr_r]))[0][2]
        fit_ok = True
    except (ValueError, np.linalg.LinAlgError):
        b_adj = b_full = float("nan")
        fit_ok = False
    return {"family": fam, "n": int(len(df)), "young": young, "interm": interm,
            "or_unmatched": un.get("or_mh", np.nan), "or_matched": mt.get("or_mh", np.nan),
            "lo": mt.get("lo", np.nan), "hi": mt.get("hi", np.nan),
            "strata": mt.get("n_strata", 0), "zero": mt.get("zero_cell_strata", 0),
            "rows": mt.get("rows", []), "tier_or": tier.get("or_mh", np.nan),
            "sat": sat, "lin_unmatched": un_lin.get("or_mh", np.nan),
            "lin_matched": mt_lin.get("or_mh", np.nan), "b_adj": b_adj, "b_full": b_full,
            "fit_ok": fit_ok}


def family_report(df: pd.DataFrame, fam: str) -> tuple[list[str], dict]:
    s = family_summary(df, fam)
    lines: list[str] = []
    add = lines.append
    add(f"## {fam} -- {s['n']} calls")
    add("")
    add("| min(SR_L,SR_R) | not in catalog | in catalog |")
    add("|---|---|---|")
    for b in BUCKETS:
        g = df[df["min_sr_bucket"] == b]
        add(f"| {b} | {int((~g['in_catalog']).sum())} | {int(g['in_catalog'].sum())} |")
    add("")
    add("| axis | unmatched OR | evidence-matched OR (MH) | 95% CI (Woolf) | strata | zero-cell strata |")
    add("|---|---|---|---|---|---|")
    for ax, label in (("in_catalog", "catalog (KNOWNMEI)"), ("in_1kg", "1000G id"),
                      ("in_longread", "long-read id")):
        un = pooled(df, ax, TSD_MIN_BP, ["three_prime_resolved"]) or {}
        mt = pooled(df, ax, TSD_MIN_BP, ["three_prime_resolved", "min_sr_bucket"]) or {}
        add(f"| {label} | {un.get('or_mh', np.nan):.2f} | {mt.get('or_mh', np.nan):.2f} | "
            f"{mt.get('lo', np.nan):.2f} - {mt.get('hi', np.nan):.2f} | {mt.get('n_strata', 0)} | "
            f"{mt.get('zero_cell_strata', 0)} |")
    add("")
    add("| stratum | not-in-catalog TSD+ | in-catalog TSD+ | OR | Fisher p |")
    add("|---|---|---|---|---|")
    for r in s["rows"]:
        add(f"| {r['stratum']} | {r['unexposed']} | {r['exposed']} | {r['or']:.2f} | {r['p']:.3g} |")
    add("")
    if s["sat"]:
        add("| subset | not in catalog TSD+ | in catalog TSD+ | OR | Fisher p |")
        add("|---|---|---|---|---|")
        for lab, v in s["sat"].items():
            add(f"| {lab} | {v['c']}/{v['c'] + v['d']} | {v['a']}/{v['a'] + v['b']} | "
                f"{v['or']:.2f} | {v['p']:.3g} |")
        add("")
    add(f"- call-tier-matched catalog OR: {s['tier_or']:.2f}")
    s["attenuation"] = (100 * (1 - s["b_full"] / s["b_adj"])
                        if s["fit_ok"] and s["b_adj"] else float("nan"))
    if s["fit_ok"]:
        add(f"- catalog coefficient: unadjusted {s['b_adj']:+.3f}, after split-read terms "
            f"{s['b_full']:+.3f} (log-OR attenuation {s['attenuation']:.1f}%)")
    else:
        add("- catalog coefficient: the evidence-adjusted logistic fit did not converge on this "
            "family, so the continuous adjustment is unavailable here and T3/T5 stand alone.")
    add(f"- lineage control ({s['young']} vs {s['interm']}): matched MH OR "
        f"{s['lin_matched']:.2f} vs unmatched {s['lin_unmatched']:.2f}")
    add("")
    or_mt, lo = s["or_matched"], s["lo"]
    sat_or = s["sat"].get("min_sr >= 1", {}).get("or", np.nan)
    ci_excludes_1 = bool(np.isfinite(lo) and lo > 1.0)
    if not np.isfinite(or_mt) or s["strata"] == 0:
        s["verdict"] = "not-assessable"
        add(f"**Cannot be assessed on {fam}.** No evidence-matched stratum holds both catalog")
        add("states, so the comparison this test is built on does not exist here. The unmatched")
        add("number must not be read as having survived.")
    elif ci_excludes_1 and or_mt >= 1.5 and np.isfinite(sat_or) and sat_or >= 1.5:
        s["verdict"] = "survives"
        add("**Survives.** The catalog association stays large, its interval excludes 1 after")
        add("matching on junction evidence, and it holds inside the evidence-saturated subset, so")
        add("better-resolved junctions do not account for it. That does not make the catalog axis")
        add("an age axis: catalog membership still confounds age with ascertainment.")
    elif (not ci_excludes_1) or or_mt < 1.5:
        s["verdict"] = "does-not-survive"
        add("**Does not survive.** Held fixed, junction evidence takes the catalog association to")
        add("1.5 or below with an interval including 1, so it is substantially an evidence-quality")
        add("effect: catalog-matched insertions have better-resolved junctions, and those are where")
        add("a TSD is observable at all.")
    else:
        s["verdict"] = "instruments-disagree"
        add("**Instruments disagree -- reported as such.** The matched interval and the")
        add("evidence-saturated subset do not agree, so at this sample size the design cannot")
        add("separate them and no single OR should be quoted.")
    add("")
    return lines, s


def family_wide_main(callsets: Path, families: list[str]) -> int:
    df = load_family_wide(callsets)
    lines: list[str] = []
    add = lines.append
    add("# Is the catalog TSD effect just junction evidence quality? -- all families")
    add("")
    add(f"Cohort: {len(df)} calls rebuilt from `{callsets}` via `span_work/line1_span_support.py` "
        "(same sentinel rule as the published LINE1 export).")
    add(f"TSD+ = TSD length >= {TSD_MIN_BP} bp; 3' resolution = `MEI_3P >= {MARKER_HI}`; "
        "strata = (3' resolution x min(SR_L,SR_R) bucket).")
    add("")
    summaries: list[dict] = []
    parts: list[str] = []
    for fam in families:
        if fam not in FAMILY_LINEAGE:
            add(f"- {fam}: no lineage rule defined for this family -- skipped.")
            continue
        sub = df[df["family"] == fam]
        if sub.empty:
            add(f"- {fam}: no calls in this cohort -- skipped.")
            continue
        fam_lines, s = family_report(sub, fam)
        summaries.append(s)
        parts.append("\n".join(fam_lines))
    add("")
    add("## Cross-family comparison")
    add("")
    add("| family | n | catalog unmatched OR | evidence-matched OR (MH) | 95% CI | strata | saturated OR (min_sr>=1) | verdict |")
    add("|---|---|---|---|---|---|---|---|")
    for s in summaries:
        add(f"| {s['family']} | {s['n']} | {s['or_unmatched']:.2f} | {s['or_matched']:.2f} | "
            f"{s['lo']:.2f} - {s['hi']:.2f} | {s['strata']} | "
            f"{s['sat'].get('min_sr >= 1', {}).get('or', np.nan):.2f} | {s['verdict']} |")
    add("")
    add("An evidence-quality artefact should move every family the same way; a real effect need not.")
    add("")
    add("The LINE1 control below is the published cohort reproduced through this route. If it does")
    add("not land on 2186 calls, the input directory is not the one the published numbers came from.")
    add("")
    line1 = df[df["family"] == "LINE1"]
    line1_s = next((s for s in summaries if s["family"] == "LINE1"), None)
    add(f"- LINE1 calls rebuilt: {len(line1)} (published cohort: 2186)")
    if line1_s is not None:
        add(f"- LINE1 catalog MH OR rebuilt: {line1_s['or_unmatched']:.2f} (published: 6.38)")
    add("")
    report = "\n".join(lines) + "\n" + "\n".join(parts)
    OUT_FAMILIES.write_text(report)
    print(report)
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description="evidence-matched TSD test")
    ap.add_argument("--callsets", type=Path, default=None,
                    help="directory of per-sample call VCFs: runs the design for every family "
                         "(LINE1, ALU, SVA) instead of the LINE1-only TSV route")
    ap.add_argument("--families", default="LINE1,ALU,SVA",
                    help="comma-separated families for the --callsets route")
    args = ap.parse_args()
    if args.callsets is not None:
        return family_wide_main(args.callsets,
                                [f.strip() for f in args.families.split(",") if f.strip()])

    df = load()
    lines: list[str] = []
    add = lines.append
    y5 = (df["tsd_len"] >= TSD_MIN_BP).astype(int).to_numpy()

    add("# Is the catalog TSD effect just junction evidence quality?")
    add("")
    add(f"Cohort: {len(df)} L1 calls, `span_work/out/line1_calls.tsv`. TSD+ = TSD length >= "
        f"{TSD_MIN_BP} bp; 3' resolution = `MEI_3P >= {MARKER_HI}`.")
    add("")
    add("`min_sr = min(SR_L, SR_R)` is the primary evidence score: both breakpoints observed by")
    add("split reads. `SR_L`, `SR_R`, `DPE_L`, `DPE_R`, `MEI_MAPPED`, `polyA_MAPPED` are read as")
    add("numeric columns; `SUPPORT` carries the same numbers as one string.")
    add("")

    # ---- T1 the imbalance the matching has to remove ----------------------
    add("## T1. The evidence imbalance the matching has to remove")
    add("")
    add("| min(SR_L,SR_R) | not in catalog | in catalog |")
    add("|---|---|---|")
    for b in BUCKETS:
        e = df[df["min_sr_bucket"] == b]
        add(f"| {b} | {int((~e['in_catalog']).sum())} | {int(e['in_catalog'].sum())} |")
    add("")
    for ax in ("in_catalog", "in_1kg", "in_longread"):
        if df[ax].nunique() == 2:
            add(f"- `{ax}`: median min_sr "
                f"{df[~df[ax]]['min_sr'].median():.0f} (n={int((~df[ax]).sum())}) vs "
                f"{df[df[ax]]['min_sr'].median():.0f} (n={int(df[ax].sum())}); share with "
                f"min_sr>=1 {(df[~df[ax]]['min_sr'] >= 1).mean():.3f} vs "
                f"{(df[df[ax]]['min_sr'] >= 1).mean():.3f}")
    add("")

    # ---- T2 does evidence itself predict TSD? -----------------------------
    add("## T2. Does evidence quality itself predict TSD presence?")
    add("")
    add("| 3'-resolved | min(SR_L,SR_R) | TSD+ | n | rate |")
    add("|---|---|---|---|---|")
    for rz in (True, False):
        for b in BUCKETS:
            g = df[(df["three_prime_resolved"] == rz) & (df["min_sr_bucket"] == b)]
            if len(g):
                add(f"| {rz} | {b} | {int((g['tsd_len'] >= TSD_MIN_BP).sum())} | {len(g)} | "
                    f"{(g['tsd_len'] >= TSD_MIN_BP).mean():.3f} |")
    add("")

    # ---- T3 evidence-matched pooling --------------------------------------
    add("## T3. Evidence-matched pooling (the decisive instrument)")
    add("")
    add("Strata = (3' resolution x min_sr bucket). The unmatched column is the number under test.")
    add("")
    add("| axis | unmatched OR | evidence-matched OR (MH) | 95% CI (Woolf) | strata | same direction | zero-cell strata |")
    add("|---|---|---|---|---|---|---|")
    res_cols = ["three_prime_resolved"]
    mat_cols = ["three_prime_resolved", "min_sr_bucket"]
    out = {}
    for ax, label in (("in_catalog", "catalog (KNOWNMEI)"),
                      ("in_1kg", "1000G id"), ("in_longread", "long-read id")):
        un = pooled(df, ax, TSD_MIN_BP, res_cols)
        mt = pooled(df, ax, TSD_MIN_BP, mat_cols)
        out[ax] = {"unmatched": un, "matched": mt}
        add(f"| {label} | {un['or_mh']:.2f} | {mt['or_mh']:.2f} | {mt['lo']:.2f} - {mt['hi']:.2f} | "
            f"{mt['n_strata']} | {mt['direction_same']}/{mt['n_strata']} | {mt['zero_cell_strata']} |")
    add("")
    add("Per-stratum detail, catalog axis:")
    add("")
    add("| stratum | not-in-catalog TSD+ | in-catalog TSD+ | OR | Fisher p |")
    add("|---|---|---|---|---|")
    for r in out["in_catalog"]["matched"]["rows"]:
        add(f"| {r['stratum']} | {r['unexposed']} | {r['exposed']} | {r['or']:.2f} | {r['p']:.3g} |")
    add("")
    add("Secondary matching, adding call tier (finer strata, smaller cells):")
    add("")
    add("| axis | matched OR (MH) | 95% CI | strata | zero-cell strata |")
    add("|---|---|---|---|---|")
    for ax, label in (("in_catalog", "catalog (KNOWNMEI)"), ("in_longread", "long-read id")):
        t3 = pooled(df, ax, TSD_MIN_BP, ["three_prime_resolved", "min_sr_bucket", "CALLTIER"])
        add(f"| {label} | {t3['or_mh']:.2f} | {t3['lo']:.2f} - {t3['hi']:.2f} | {t3['n_strata']} | "
            f"{t3['zero_cell_strata']} |")
    add("")

    # ---- T4 continuous adjustment ----------------------------------------
    add("## T4. Continuous adjustment: does the catalog term survive evidence?")
    add("")
    res = df["three_prime_resolved"].astype(int).to_numpy(dtype=float)
    cat = df["in_catalog"].astype(int).to_numpy(dtype=float)
    sr_l = np.log1p(df["SR_L"].to_numpy(dtype=float))
    sr_r = np.log1p(df["SR_R"].to_numpy(dtype=float))
    pa = np.log1p(df["polyA_MAPPED"].to_numpy(dtype=float))
    one = np.ones(len(df))
    models = [
        ("TSD ~ resolution", np.column_stack([one, res]), ["intercept", "resolution"]),
        ("TSD ~ resolution + catalog", np.column_stack([one, res, cat]),
         ["intercept", "resolution", "catalog"]),
        ("TSD ~ resolution + catalog + split-read", np.column_stack([one, res, cat, sr_l, sr_r]),
         ["intercept", "resolution", "catalog", "log_SR_L", "log_SR_R"]),
        ("TSD ~ resolution + catalog + split-read + polyA",
         np.column_stack([one, res, cat, sr_l, sr_r, pa]),
         ["intercept", "resolution", "catalog", "log_SR_L", "log_SR_R", "log_polyA"]),
        ("TSD ~ resolution * catalog", np.column_stack([one, res, cat, res * cat]),
         ["intercept", "resolution", "catalog", "resolution:catalog"]),
    ]
    betas = {}
    for name, X, labels in models:
        add(base.model_row(name, y5, X, labels))
        beta, _ = base.logit(y5, X)
        betas[name] = dict(zip(labels, beta))
    add("")
    b_adj = betas["TSD ~ resolution + catalog"]["catalog"]
    b_full = betas["TSD ~ resolution + catalog + split-read"]["catalog"]
    b_full2 = betas["TSD ~ resolution + catalog + split-read + polyA"]["catalog"]
    add(f"Catalog coefficient: unadjusted **{b_adj:+.3f}**, after split-read terms "
        f"**{b_full:+.3f}** (log-OR attenuation {100 * (1 - b_full / b_adj):.1f}%), after "
        f"split-read and poly(A) **{b_full2:+.3f}**.")
    add("")

    # ---- T5 evidence-saturated subset ------------------------------------
    add("## T5. The evidence-saturated subset")
    add("")
    add("Restricting to calls whose junction evidence is already good removes the confound by")
    add("construction, at the cost of sample size.")
    add("")
    add("| subset | not in catalog TSD+ | in catalog TSD+ | OR | Fisher p |")
    add("|---|---|---|---|---|")
    saturated = {}
    for lab, sel in (("min_sr >= 1", df["min_sr"] >= 1),
                     ("min_sr >= 2", df["min_sr"] >= 2),
                     ("sr_total >= 5", df["sr_total"] >= 5)):
        s = df[sel]
        e = s[s["in_catalog"]]
        u = s[~s["in_catalog"]]
        if len(e) and len(u):
            a = int((e["tsd_len"] >= TSD_MIN_BP).sum()); b = len(e) - a
            c = int((u["tsd_len"] >= TSD_MIN_BP).sum()); d = len(u) - c
            orv = ((a * d) / (b * c)) if b * c else np.inf
            _, p = stats.fisher_exact([[a, b], [c, d]])
            saturated[lab] = {"or": orv, "p": p, "e": (a, a + b), "u": (c, c + d)}
            add(f"| {lab} | {c}/{c + d} ({c / max(c + d, 1):.3f}) | {a}/{a + b} "
                f"({a / max(a + b, 1):.3f}) | {orv:.2f} | {p:.3g} |")
    add("")

    # ---- T6 symmetric control --------------------------------------------
    add("## T6. Symmetric control -- the same machinery on the lineage axis")
    add("")
    add("Matching must not be able to erase any association it is pointed at, or T3 proves nothing.")
    add("")
    lin = df[df["age"].isin(["young:L1HS", "interm:L1P*"])].copy()
    lin["is_young"] = (lin["age"] == "young:L1HS").astype(int)
    un_lin = pooled(lin, "is_young", TSD_MIN_BP, res_cols)
    mt_lin = pooled(lin, "is_young", TSD_MIN_BP, mat_cols)
    add("| axis | unmatched OR | evidence-matched OR (MH) | strata | same direction |")
    add("|---|---|---|---|---|")
    add(f"| L1HS vs L1P* | {un_lin['or_mh']:.2f} | {mt_lin['or_mh']:.2f} | {mt_lin['n_strata']} | "
        f"{mt_lin['direction_same']}/{mt_lin['n_strata']} |")
    add("")

    # ---- T7 verdict, selected by the numbers above ------------------------
    add("## T7. Verdict (selected by the numbers above, not written in advance)")
    add("")
    c = out["in_catalog"]
    or_un = c["unmatched"]["or_mh"]
    or_mt = c["matched"]["or_mh"]
    lo, hi = c["matched"]["lo"], c["matched"]["hi"]
    sat = saturated.get("min_sr >= 1", {}).get("or", float("nan"))
    log_atten = 1 - np.log(or_mt) / np.log(or_un) if or_un > 1 and or_mt > 0 else float("nan")
    ci_excludes_1 = np.isfinite(lo) and lo > 1.0
    add(f"- catalog effect, unmatched: MH OR {or_un:.2f}")
    add(f"- catalog effect, evidence-matched: MH OR {or_mt:.2f} "
        f"(Woolf 95% CI {lo:.2f}-{hi:.2f}; {c['matched']['usable']} strata usable for the interval, "
        f"{c['matched']['zero_cell_strata']} with a zero cell)")
    add(f"- log-OR attenuation under matching: {100 * log_atten:.1f}%")
    add(f"- evidence-saturated subset (min_sr >= 1): OR {sat:.2f}")
    add(f"- lineage control: matched MH OR {mt_lin['or_mh']:.2f} vs unmatched {un_lin['or_mh']:.2f}")
    add("")
    if ci_excludes_1 and or_mt >= 1.5 and np.isfinite(sat) and sat >= 1.5:
        add("**Survives.** The catalog association stays large and its interval excludes 1 after")
        add("matching on junction evidence, and it also holds inside the evidence-saturated subset,")
        add("so better-resolved junctions do not account for it. On this cohort the catalog axis is")
        add("not merely an evidence-quality proxy -- which does not make it an age axis: catalog")
        add("membership still confounds age with ascertainment.")
    elif not ci_excludes_1 or or_mt < 1.5:
        add("**Does not survive.** Once junction evidence is held fixed the catalog association")
        add("falls to or below 1.5 with an interval that includes 1, so the catalog-membership")
        add("effect reported in `tsd_age_association.py` (section 4, MH OR 6.38) is substantially an")
        add("evidence-quality effect: catalog-matched insertions have better-resolved junctions, and")
        add("better-resolved junctions are where a TSD is observable at all. The catalog axis should")
        add("then not be read as an age axis.")
    else:
        add("**Instruments disagree -- reported as such.** The matched interval and the")
        add("evidence-saturated subset do not agree on whether the effect survives, and at this")
        add("cohort size that means the design cannot separate them. The numbers above are the")
        add("bound; no single OR should be quoted from this run.")
    add("")

    report = "\n".join(lines) + "\n"
    OUT.write_text(report)
    print(report)
    return 0


if __name__ == "__main__":
    sys.exit(main())
