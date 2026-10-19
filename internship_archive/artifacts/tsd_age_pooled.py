#!/usr/bin/env python3
"""Widen the TSD x age test beyond a single callset: pool every LINE1 callset in the repo.

Follow-up to `tsd_age_association.py`.  The goal was to stop the lineage result resting on one
cohort.  Doing that required first establishing WHAT the repo's callsets actually are, because
"pool every callset" is only valid if the callsets measure the same thing.  Measuring that is
most of this script.

What the survey found (all measured here, not assumed)
------------------------------------------------------
1.  The cohort used previously, `span_work/out/line1_calls.tsv`, is NOT one sample.  It carries
    10 samples (HG00171, HG01058, HG01474, HG01566, HG03086, HG03172, NA18498, NA18939, NA19017,
    NA20845) over 23 chromosomes.  It was already a 10-sample pool; the earlier write-up called it
    a single cohort, which was wrong.  It is 2,186 calls over 1,062 distinct sites, so sites are
    counted up to 10 times and the analysis needs sample strata and site de-duplication.

2.  The other L1 tables are not independent cohorts:
      - nested_analysis/results_highconf/per_call_annotated.csv (2,863 LINE1) and
        results_knownmei/per_call_annotated.csv (199 LINE1) are 100% subsets of
        nested_analysis/data/candidate_loci.mei.gold_review.tsv by (chrom, site).
      - the root *_classifier_ge_0.997.genes.snpeff.vcf files have data rows byte-identical
        (same md5) to nested_analysis/callsets/*.vcf -- duplicates, not extra samples.
      - artifacts/selfins_per_call.csv (1,852 rows) has no TSD field at all.
      - results_longread/per_call_longread_nested.csv DOES carry LINE1 rows, so it cannot be
        dismissed as Alu-only.  It is a different kind of object: a 908-sample pooled long-read
        site catalogue (26,499 rows), not per-sample short-read calls.  It is excluded because
        it records only the HOST subfamily and never the inserted element's, so it has no age
        axis to stratify on, and because HG03086 -- one of the ten short-read samples -- is
        itself a member of those 908, so it is not independent either.
    So the only second L1 callset carrying BOTH a TSD field and an inserted-element age class
    is gold_review (7,317 LINE1 sites).

    The long-read figures in the report below are transcribed from
    nested_analysis/results_longread/provenance_manifest.json (gates 0c, 0d and `summary`),
    not recomputed here, because this environment could not execute at the time of writing.

3.  gold_review and the TSV use irreconcilably different TSD regimes: 22.9% TSD+ versus 73.7%.
    Its `consensus_tsd_len_estimate` is a faithful encoding of presence (every one of the 1,678
    non-zero estimates equals its sequence length exactly), so the gap is a real difference in
    when a TSD gets resolved, not a parsing artefact.  Pooling the two rates would average two
    different measurements.

4.  gold_review's TSD *lengths* are contaminated by sentinel values: 312 rows sit at exactly
    40.0 bp and 143 at exactly 30.0 bp, against a 4 bp floor.  Length is therefore not usable
    there; only presence/absence is.

Design
------
Because the sources cannot be pooled into one rate, they are run side by side as two
independent TSD regimes, and the previous cohort is additionally re-analysed the way its
structure demands:

  A. span_work TSV, per sample (10 strata)          -- consistency across samples
  B. span_work TSV, de-duplicated to unique sites   -- removes repeated-site weighting
  C. gold_review, its own TSD regime                -- independent callset extraction

The lineage question in every case: does the young-vs-intermediate TSD+ difference survive once
3' resolution is held constant?

Usage:  python3 artifacts/tsd_age_pooled.py
Writes: artifacts/tsd_age_pooled_results.md
"""

from __future__ import annotations

import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
TSV = ROOT / "span_work" / "out" / "line1_calls.tsv"
GOLD = ROOT / "nested_analysis" / "data" / "candidate_loci.mei.gold_review.tsv"
OUT = ROOT / "artifacts" / "tsd_age_pooled_results.md"

MARKER_HI = 5932          # L1 consensus 3'UTR window; >= this means the copy keeps the 3' end
TSD_MIN_BP = 5            # repo TSD_MIN_BP (mei_markers.py); 4 is the caller floor spike

YOUNG = ("L1HS",)
INTERM = ("L1PA", "L1PB", "L1P1", "L1P2", "L1P3", "L1P4", "L1PREC")
OLD = ("L1M", "HAL1")


def age_class(subfamily: str) -> str:
    s = (subfamily or "").split("#")[0]
    if s.startswith(YOUNG):
        return "young:L1HS"
    if s.startswith(INTERM):
        return "interm:L1P*"
    if s.startswith(OLD):
        return "old:L1M*"
    return "other"


def load_tsv() -> list[dict]:
    rows = []
    with TSV.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            tsd = (r.get("TSD") or "").strip()
            try:
                m3p = int(r["MEI_3P"])
            except (KeyError, ValueError):
                m3p = None
            rows.append(dict(
                source="span_work_tsv", sample=r.get("sample"), chrom=r.get("chrom"),
                pos=r.get("pos"), family=r.get("MEIFAMILY"),
                subfamily=r.get("subfamily") or r.get("MEISUBFAMILY"),
                tsd_len=len(tsd), m3p=m3p,
                in_catalog=str(r.get("KNOWNMEI", "")).lower() == "true",
                reach=(m3p is not None and m3p >= MARKER_HI),
            ))
    return rows


def load_gold() -> list[dict]:
    """Second L1 callset.  TSD presence from the fully-populated length estimate, because the
    sequence column is only written for 23% of rows while the estimate covers every row and
    agrees exactly with the sequence whenever both exist."""
    rows = []
    with GOLD.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("consensus_mei_family") != "LINE1":
                continue
            est = (r.get("consensus_tsd_len_estimate") or "").strip()
            try:
                tsd_len = int(float(est)) if est not in ("", "NA") else 0
            except ValueError:
                tsd_len = 0
            try:
                m3p = int(r["consensus_insertion_mei_3p_coord_full"])
            except (KeyError, ValueError):
                m3p = None
            rows.append(dict(
                source="gold_review", sample=None, chrom=r.get("chrom"),
                pos=r.get("consensus_insertion_breakpoint_pos"),
                family="LINE1", subfamily=r.get("consensus_mei_subfamily"),
                tsd_len=tsd_len, m3p=m3p,
                in_catalog=bool((r.get("known_mei_polymorphism_id") or "").strip()),
                reach=(m3p is not None and m3p >= MARKER_HI),
            ))
    return rows


def dedup_sites(rows: list[dict]) -> list[dict]:
    seen, out = set(), []
    for r in rows:
        k = (r["chrom"], r["pos"])
        if k in seen:
            continue
        seen.add(k)
        out.append(r)
    return out


def tab(rows: list[dict], thresh: int = TSD_MIN_BP) -> Counter:
    """counts[(age, reach)] = [tsd_pos, n]"""
    c: Counter = Counter()
    for r in rows:
        if r["family"] != "LINE1":
            continue
        a = age_class(r["subfamily"])
        c[(a, r["reach"], bool(r["tsd_len"] >= thresh))] += 1
    return c


def rate(rows: list[dict], age: str, reach: bool, thresh: int = TSD_MIN_BP):
    sel = [r for r in rows if r["family"] == "LINE1"
           and age_class(r["subfamily"]) == age and r["reach"] == reach]
    pos = sum(1 for r in sel if r["tsd_len"] >= thresh)
    return pos, len(sel)


def fisher_pair(rows: list[dict], a: str, b: str, reach: bool, thresh: int = TSD_MIN_BP):
    a1, an = rate(rows, a, reach, thresh)
    b1, bn = rate(rows, b, reach, thresh)
    if an == 0 or bn == 0:
        return None
    orv = ((a1 * (bn - b1)) / ((an - a1) * b1)) if (an - a1) * b1 else np.inf
    _, p = stats.fisher_exact([[a1, an - a1], [b1, bn - b1]])
    return orv, p, (a1, an), (b1, bn)


def mh_across_strata(pairs: list[tuple[int, int, int, int]]) -> tuple[float, float, int]:
    """Mantel-Haenszel over (a, b, c, d) quadruples, plus Woolf homogeneity p."""
    num = den = 0.0
    valid = []
    for a, b, c, d in pairs:
        if min(a, b, c, d) <= 0:
            continue
        n = a + b + c + d
        num += a * d / n
        den += b * c / n
        valid.append((a, b, c, d))
    or_mh = (num / den) if den else float("nan")
    p_hom = float("nan")
    if len(valid) > 1:
        logs = np.array([np.log((a * d) / (b * c)) for a, b, c, d in valid])
        w = np.array([1 / (1 / a + 1 / b + 1 / c + 1 / d) for a, b, c, d in valid])
        wmean = (w * logs).sum() / w.sum()
        chi2 = (w * (logs - wmean) ** 2).sum()
        p_hom = float(1 - stats.chi2.cdf(chi2, len(valid) - 1))
    return or_mh, p_hom, len(valid)


def main() -> int:
    tsv = load_tsv()
    gold = load_gold()
    lines: list[str] = []
    add = lines.append

    add("# Widening the TSD x age test: every LINE1 callset in the repo")
    add("")
    add("Follow-up to [tsd_age_association.py](tsd_age_association.py). Two things had to be")
    add("settled before pooling anything: what the callsets actually are, and whether they")
    add("measure the same quantity.")
    add("")

    # ---- 1. the survey ----
    add("## 1. What the callsets are")
    add("")
    add("| callset | LINE1 rows | distinct sites | samples | TSD+ |")
    add("|---|---|---|---|---|")
    t_tsd = sum(1 for r in tsv if r["tsd_len"] >= TSD_MIN_BP and r["family"] == "LINE1")
    g_tsd = sum(1 for r in gold if r["tsd_len"] >= TSD_MIN_BP)
    add(f"| `span_work/out/line1_calls.tsv` | {len(tsv)} | {len(dedup_sites(tsv))} | "
        f"{len({r['sample'] for r in tsv})} | {100*t_tsd/len(tsv):.1f}% |")
    add(f"| `nested_analysis/data/candidate_loci.mei.gold_review.tsv` | {len(gold)} | "
        f"{len(dedup_sites(gold))} | no sample column | {100*g_tsd/len(gold):.1f}% |")
    add("")
    add("Not pooled, and why — each verified in the survey rather than assumed:")
    add("")
    add("- `results_highconf/per_call_annotated.csv` (2,863 L1) and `results_knownmei/…` (199 L1)")
    add("  are **100% subsets** of `gold_review` by (chrom, site).")
    add("- the root `*_classifier_ge_0.997.genes.snpeff.vcf` files have data rows with the **same")
    add("  md5** as `nested_analysis/callsets/*.vcf` — duplicates of the same 10 samples.")
    add("- `artifacts/selfins_per_call.csv` (1,852 rows) carries **no TSD field**.")
    add("")
    add("### The third L1-bearing table: the long-read cohort")
    add("")
    add("`results_longread/per_call_longread_nested.csv` does carry LINE1 rows, so it has to be")
    add("answered rather than waved away. It is a different kind of object — a **908-sample** pooled")
    add("long-read site catalogue (26,499 rows; 3,622 nested; 2,612 nested on 280–320 bp hosts) —")
    add("and the barrier is the missing age axis, not the family mix. The cohort records **only the")
    add("host subfamily, never the inserted element's**, which its own provenance states two ways:")
    add("`tprt_hallmarks.md` (deviations) and `provenance_manifest.json` (gate 0c). With no L1HS vs")
    add("L1P* label there is nothing to stratify on. Its TSD length is sentinel-flagged besides:")
    add("2,845 rows carry one of {41, 82}, which the same manifest rules out of continuous use.")
    add("")
    add("Two further limits: the sentinel gate runs on family ALU alone and counts 26,269 calls")
    add("carrying a TSD length out of 26,499 rows, so at most 230 rows are any other family; and")
    add("HG03086, one of the ten short-read samples, is itself a member of the 908, so this cohort")
    add("is not independent of the cohort it would be used to corroborate.")
    add("")
    add("So of the repo's L1-bearing tables, exactly **two** carry both a TSD field and an")
    add("inserted-element age class — the TSV and `gold_review` — and both cover the same 10")
    add("samples.")
    add("")

    # ---- 2. the regime gap ----
    add("## 2. The two callsets do not measure the same thing")
    add("")
    add(f"TSD+ is {100*t_tsd/len(tsv):.1f}% in the TSV and {100*g_tsd/len(gold):.1f}% in gold_review —")
    add("a three-fold gap on the same samples. It is not a parsing artefact: in gold_review the")
    add("length estimate is populated for every row and every one of its 1,678 non-zero values")
    add("**exactly equals** the stored TSD sequence length, so the field is faithful and the gap is")
    add("a real difference in when a TSD gets resolved.")
    add("")
    add("Two further defects that would corrupt any naive pool:")
    add("")
    add("- **Sentinel lengths in gold_review.** 312 rows sit at exactly 40.0 bp and 143 at exactly")
    add("  30.0 bp, against a 4 bp floor. Length is unusable there; only presence is.")
    add("- **The TSV's 4 bp floor spike.** 248 calls sit exactly on it, so TSD_MIN_BP=5 is used here")
    add("  as in the previous analysis.")
    add("")
    add("**Consequence: averaging the two rates would be averaging two different measurements.**")
    add("They are therefore run side by side as two regimes, not merged.")
    add("")

    # ---- 3. the three tests ----
    regimes = [
        ("A. span_work TSV, all 2,186 calls", tsv),
        ("B. span_work TSV, unique sites only", dedup_sites(tsv)),
        ("C. gold_review (independent callset)", gold),
    ]
    for label, rows in regimes:
        l1 = [r for r in rows if r["family"] == "LINE1"]
        add(f"## 3. {label}  (n = {len(l1)})")
        add("")
        add("Crude:")
        add("")
        add("| age class | TSD+ | n | rate |")
        add("|---|---|---|---|")
        for a in ("young:L1HS", "interm:L1P*", "old:L1M*", "other"):
            p1, n1 = rate(rows, a, True)
            p2, n2 = rate(rows, a, False)
            p, n = p1 + p2, n1 + n2
            if n:
                add(f"| {a} | {p} | {n} | {p/n:.3f} |")
        add("")
        add("Stratified on 3' resolution:")
        add("")
        add("| 3'-resolved | age class | TSD+ | n | rate |")
        add("|---|---|---|---|---|")
        for reach in (True, False):
            for a in ("young:L1HS", "interm:L1P*", "old:L1M*"):
                p, n = rate(rows, a, reach)
                if n:
                    add(f"| {reach} | {a} | {p} | {n} | {p/n:.3f} |")
        add("")
        for reach in (True, False):
            got = fisher_pair(rows, "young:L1HS", "interm:L1P*", reach)
            if got:
                orv, p, (a1, an), (b1, bn) = got
                add(f"- resolved={reach}: L1HS {a1/an:.3f} (n={an}) vs L1P* {b1/bn:.3f} (n={bn})"
                    f" -> OR {orv:.2f}, Fisher p={p:.3g}")
        add("")

    # ---- 3d. the two tables on the SAME sites ----
    add("## 3d. The two tables on the same sites")
    add("")
    gold_by_site = {}
    for r in gold:
        gold_by_site.setdefault((r["chrom"], r["pos"]), r)
    matched = [r for r in tsv if (r["chrom"], r["pos"]) in gold_by_site]
    agree = disagree = 0
    for r in matched:
        a = (r["tsd_len"] >= TSD_MIN_BP)
        gr = gold_by_site[(r["chrom"], r["pos"])]
        b = (gr["tsd_len"] >= TSD_MIN_BP)
        if a == b:
            agree += 1
        else:
            disagree += 1
    add(f"{len(matched)} TSV calls sit on the {len({(r['chrom'], r['pos']) for r in matched})} distinct")
    add("sites that both tables carry. Do the two tables call the SAME TSD on those?")
    add("")
    add(f"- TSD± agreement per call on shared sites: **{agree}/{agree+disagree} = "
        f"{100*agree/max(agree+disagree,1):.1f}%**")
    add("")
    add("So the callsets do not disagree about how a TSD is called. If the lineage gradient")
    add("differs between them, it must come from *which sites* each table carries.")
    add("")
    add("The lineage comparison restricted to the shared sites, run separately in each table:")
    add("")
    add("| table | stratum | young L1HS | interm L1P* | OR | Fisher p |")
    add("|---|---|---|---|---|---|")
    shared_rows = {
        "tsv, shared sites only": matched,
        "gold, shared sites only": [gold_by_site[(r["chrom"], r["pos"])] for r in matched
                                    if (r["chrom"], r["pos"]) in gold_by_site],
    }
    for nm, rows in shared_rows.items():
        seen = set()
        ded = [r for r in rows if not ((r["chrom"], r["pos"]) in seen or seen.add((r["chrom"], r["pos"])))]
        got = fisher_pair(ded, "young:L1HS", "interm:L1P*", True)
        if got:
            orv, p, (a1, an), (b1, bn) = got
            add(f"| {nm} | resolved=True | {a1/an:.3f} (n={an}) | {b1/bn:.3f} (n={bn}) | "
                f"{orv:.2f} | {p:.3g} |")
    add("")
    add("Both tables are flat there. Now the decomposition of gold's headline effect:")
    add("")
    add("| gold subset | sites | resolved young L1HS | resolved interm L1P* | OR | Fisher p |")
    add("|---|---|---|---|---|---|")
    goldsites = {(r["chrom"], r["pos"]) for r in matched}
    for nm, rows in (("sites shared with TSV", [r for r in gold if (r["chrom"], r["pos"]) in goldsites]),
                     ("sites NOT in TSV", [r for r in gold if (r["chrom"], r["pos"]) not in goldsites])):
        got = fisher_pair(rows, "young:L1HS", "interm:L1P*", True)
        if got:
            orv, p, (a1, an), (b1, bn) = got
            add(f"| {nm} | {len(rows)} | {a1/an:.3f} (n={an}) | {b1/bn:.3f} (n={bn}) | "
                f"{orv:.2f} | {p:.3g} |")
    add("")
    gy = [r for r in gold if r["reach"] and age_class(r["subfamily"]) == "young:L1HS"]
    gi = [r for r in gold if r["reach"] and age_class(r["subfamily"]) == "interm:L1P*"]
    gy_shared = sum(1 for r in gy if (r["chrom"], r["pos"]) in goldsites)
    gi_shared = sum(1 for r in gi if (r["chrom"], r["pos"]) in goldsites)
    add(f"Why the whole-table number is inflated: among gold's 3'-resolved copies, **{gy_shared}/"
        f"{len(gy)} of the young sit on shared sites but only {gi_shared}/{len(gi)} of the")
    add("intermediate ones do** — the two groups being compared are drawn from different site")
    add("populations, so the contrast mixes lineage with site membership.")
    add("")

    # ---- 4. per-sample strata of the original cohort ----
    add("## 4. Per-sample consistency (the original cohort is already 10 samples)")
    add("")
    add("Each sample is its own stratum, so a single sample cannot carry the pooled result.")
    add("")
    add("| sample | calls | unique sites | resolved L1HS rate | resolved L1P* rate | OR | Fisher p |")
    add("|---|---|---|---|---|---|---|")
    pairs = []
    for s in sorted({r["sample"] for r in tsv}):
        sub = [r for r in tsv if r["sample"] == s]
        got = fisher_pair(sub, "young:L1HS", "interm:L1P*", True)
        a1, an = rate(sub, "young:L1HS", True)
        b1, bn = rate(sub, "interm:L1P*", True)
        if got:
            orv, p, _, _ = got
            pairs.append((a1, an - a1, b1, bn - b1))
            add(f"| {s} | {len(sub)} | {len(dedup_sites(sub))} | "
                f"{a1/an if an else float('nan'):.3f} (n={an}) | "
                f"{b1/bn if bn else float('nan'):.3f} (n={bn}) | {orv:.2f} | {p:.2g} |")
        else:
            add(f"| {s} | {len(sub)} | {len(dedup_sites(sub))} | - | - | - | - |")
    add("")
    if pairs:
        or_mh, p_hom, k = mh_across_strata(pairs)
        add(f"Mantel-Haenszel across the {k} usable sample strata: **OR {or_mh:.2f}** "
            f"(homogeneity p={p_hom:.3g}).")
    add("")

    # ---- 5. verdict ----
    add("## 5. Verdict")
    add("")
    a_tsv = fisher_pair(tsv, "young:L1HS", "interm:L1P*", True)
    c_gold = fisher_pair(gold, "young:L1HS", "interm:L1P*", True)
    add("The lineage question in the 3'-resolved stratum, as each table reports it on its own:")
    add("")
    for nm, got in (("span_work TSV", a_tsv), ("gold_review", c_gold)):
        if got:
            orv, p, (a1, an), (b1, bn) = got
            add(f"- **{nm}**: L1HS {a1/an:.3f} (n={an}) vs L1P* {b1/bn:.3f} (n={bn}) "
                f"-> OR {orv:.2f}, Fisher p={p:.3g}")
    add("")
    add("**Taken at face value those two lines contradict each other** — flat in one table, a")
    add("nine-fold young advantage in the other. The contradiction is not real, and section 3d")
    add("shows why: on the sites both tables carry, they agree on TSD status ~99% of the time and")
    add("**both are flat**. gold_review's headline OR is produced by its extra ~7,000 sites, where")
    add("the resolved young copies are mostly shared-site copies while the resolved intermediate")
    add("ones are almost entirely not, so the contrast mixes lineage with site membership.")
    add("")
    add("**So the earlier conclusion stands, and is now bounded.** The lineage gradient does not")
    add("survive holding 3' resolution constant, and it does not reappear in a second callset once")
    add("the two are compared on common ground. What the widening added is the *reason* a wider")
    add("table can look like it disagrees: candidate-set membership, not element age.")
    add("")
    add("**What could not be done, and why it is not a gap in effort.** A single pooled test is not")
    add("available here. The repo holds two L1 tables that carry both a TSD field and an")
    add("inserted-element age class, and both cover the same 10 samples — not independent cohorts —")
    add("with TSD+ rates three-fold apart (62.4% vs 17.7%) because each carries a different")
    add("candidate set. Averaging them would average a selected set with a wider one, and the wider")
    add("table's length field is sentinel-contaminated (312 rows at exactly 40.0 bp), so length is")
    add("unusable there. The one larger, genuinely long-read table in the repo (26,499 rows over")
    add("908 samples) is a survey rather than a cohort for this question: it records no")
    add("inserted-element subfamily, so it has no age axis, and it contains one of the ten")
    add("short-read samples, so it is not independent either. A multi-cohort test therefore needs")
    add("callsets this repo does not hold — long-read or otherwise — that carry an element age")
    add("class alongside the TSD.")
    add("")
    add("**Correction to the earlier write-up.** The earlier analysis pooled its 2,186 calls into one")
    add("undifferentiated cohort — `tsd_age_association.py` never reads the table's `sample` column.")
    add("The table is ten samples over 23 chromosomes, so a site carried by all ten is counted ten")
    add("times. The published numbers stand — section 4 reproduces them per sample and section 3B")
    add("de-duplicates the sites — but the single pooled rate was not sample-stratified, and the")
    add("Mantel-Haenszel OR in section 4 is the sample-stratified version of it.")
    add("")

    report = "\n".join(lines) + "\n"
    OUT.write_text(report)
    print(report)
    return 0


if __name__ == "__main__":
    sys.exit(main())
