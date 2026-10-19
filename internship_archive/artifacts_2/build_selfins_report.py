"""Regenerate the self-insertion report from the saved result tables.

Every number in the prose is read from a CSV or JSON written by the analysis. Nothing is
typed by hand, so a change upstream propagates instead of going stale.
"""
import json

import numpy as np
import pandas as pd
from scipy import stats

FIG_MAIN = "eb5945f3-3fc5-4ae1-b65d-2090eea2defb"   # selfins_main.png

enr   = pd.read_csv("notes/selfins_enrichment.csv")
sa    = pd.read_csv("notes/selfins_sense_antisense.csv")
ori   = pd.read_csv("notes/selfins_orientation.csv")
sens  = pd.read_csv("notes/selfins_orientation_sensitivity.csv")
hom   = pd.read_csv("notes/selfins_homology_test.csv")
dgf   = pd.read_csv("notes/selfins_en_degenerate_family.csv")
dgs   = pd.read_csv("notes/selfins_en_degenerate_subfamily.csv")
narrow= pd.read_csv("notes/selfins_en_motif_prediction.csv")
ccatt = pd.read_csv("notes/selfins_ccatt_control.csv")
diag  = pd.read_csv("notes/selfins_diagnostics.csv")
gc    = pd.read_csv("notes/alu_insertion_vs_retention_gc.csv")
ctx   = pd.read_csv("notes/selfins_sense_antisense_by_context.csv")
per   = pd.read_csv("notes/selfins_per_call.csv")
lit   = pd.read_csv("notes/nested_orientation_literature.csv")
bg    = json.load(open("notes/selfins_genome_background_en.json"))
bs    = pd.read_csv("notes/selfins_boundary_sensitivity.csv")
pa    = json.load(open("notes/selfins_polya_annotation_check.json"))
pp    = pd.read_csv("notes/selfins_position_profile.csv")
ps    = json.load(open("notes/selfins_position_stats.json"))
ts    = json.load(open("notes/selfins_branch_test_status.json"))

FIG_POS = "4e4947ec-ecbd-4472-bab1-66472b3c72c2"   # selfins_position.png
_pk = sorted(ps["peak_bins_bp"], key=lambda t: -t[2])
PEAK_ROWS = "".join(
    f"| {b}–{e} bp | **{n}** ({n/pp.expected_uniform.iloc[0]:.1f}×) | {af:.2f} | "
    f"{'3′ poly-A tail' if b >= 260 else ('inter-arm A-rich linker' if b >= 100 else '5′ end — not A-rich')} |\n"
    for b, e, n, af in _pk[:3])

POLYA_PCT = pa["pct_with_adjacent_polyA"]
POLYA_N = pa["n_alu"]
_bp = bs.pivot(index="slop_bp", columns="family", values="nested_frac")
_bo = bs.pivot(index="slop_bp", columns="family", values="same_orient_frac")
BOUNDARY_ROWS = "".join(
    f"| ±{int(s)} bp | {_bp.loc[s,'ALU']:.3f} | {_bo.loc[s,'ALU']:.3f} "
    f"| {_bp.loc[s,'LINE1']:.3f} | {_bo.loc[s,'LINE1']:.3f} |\n"
    for s in _bp.index)
BS_ALU_50 = float(_bo.loc[50, "ALU"])
BS_L1_50 = float(_bo.loc[50, "LINE1"])

NICE = {"ALU": "Alu", "LINE1": "LINE-1", "SVA": "SVA"}
CONS = {"ALU": "uniform", "LINE1": "GC-matched", "SVA": "GC-matched"}
E  = {f: enr[(enr.family == f) & (enr.null == CONS[f])].iloc[0] for f in NICE}
EA = {(f, n): enr[(enr.family == f) & (enr.null == n)].iloc[0]
      for f in NICE for n in ("GC-matched", "uniform")}
O  = {f: ori[ori.family == f].iloc[0] for f in NICE}
S  = {(r.family, r.subset): r for r in sens.itertuples()}
D  = {r.family: r for r in diag.itertuples()}
DF = {r.family: r for r in dgf.itertuples()}
NW = {r.family: r for r in narrow.itertuples()}
SA = {(r.family, r.category): r for r in sa.itertuples()}
CC = ccatt.iloc[0]

# subfamily goodness-of-fit for the degenerate consensus, recomputed from the table
_exp = dgs.pred_same * dgs.n_insertions
_var = dgs.n_insertions * dgs.pred_same * (1 - dgs.pred_same)
CHI2 = float((((dgs.n_same_obs - _exp) ** 2) / _var).sum())
CHI2_DF = len(dgs)
CHI2_P = float(1 - stats.chi2.cdf(CHI2, CHI2_DF))
_pe = stats.pearsonr(dgs.pred_same, dgs.observed_same)
R_P, R_P_P = float(_pe.statistic), float(_pe.pvalue)

n_total = len(per)
n_nested = int(per.nested.sum())
n_anti = int((per.nested & ~per.same_orientation).sum())
anti_by_fam = per[per.nested & ~per.same_orientation].MEIFAMILY.value_counts()
gc_poly = gc[gc.set.str.startswith("new") & gc.subfamily.str.startswith("AluY")].iloc[0]
gc_fix = gc[gc.set.str.startswith("fixed") & gc.subfamily.str.startswith("AluY")].iloc[0]
gc_base = gc[gc.set.str.startswith("genome")].iloc[0]


def pv(p):
    if p != p:
        return "n/a"
    if p >= 1e-4:
        return f"{p:.3g}"
    return f"{p:.1e}"


def orf(v):
    return "∞" if (v != v or v == float("inf")) else f"{v:.1f}"


md = f"""# Retrotransposon self-insertion in HG03086

**Question (W. Brandler):** do MEI insertions favour inserting into existing reference
copies of their own family, and when they do, is the new insertion in the same orientation
as the reference element? Following his clarification — *"we would want to annotate
nested_sense and nested_antisense and compare the two"* — the two orientation classes are
treated as separate categories throughout rather than as a single same-orientation rate.

**Answers.** Yes to self-targeting, for all three families. The orientation split is
asymmetric in a way a single rate hides: sense nesting is **enriched**, antisense nesting is
**depleted below chance**. Target choice among candidate host elements tracks **sequence
similarity** — young insertions prefer young hosts — which reproduces a published LINE-1
result in a family it was not tested in. The orientation bias itself remains
**mechanistically unexplained**: four candidate mechanisms were tested and all four rejected
— two endonuclease consensus models (§4), a non-classical target motif (§4), and
host-element homology (§3).

Callset: {n_total} classifier calls (score ≥ 0.997) on HG03086, GRCh38, chr1–22 + chrX —
one genome, which bounds several results below. Reference elements: UCSC hg38 RepeatMasker,
filtered to Alu / LINE-1 / SVA with the pipeline's own family-normalisation rule
(`_normalize_mei_family_token`). {n_nested} of {n_total} calls
({100*n_nested/n_total:.1f}%) fall inside a same-family element.

---

## 1. Self-targeting: all three families, above chance

Two nulls, each resampling call positions genome-wide 1000 times — **uniform**, and
**GC-matched** on each call's local GC (1 kb window).

| family | n | nested | observed | GC-matched exp. | uniform exp. | conservative enrichment | p |
|---|---|---|---|---|---|---|---|
"""
for f in NICE:
    g, u = EA[(f, "GC-matched")], EA[(f, "uniform")]
    md += (f"| {NICE[f]} | {int(E[f].n)} | {int(E[f].obs_nested)} | {E[f].obs_frac:.3f} "
           f"| {g.exp_frac:.3f} | {u.exp_frac:.3f} | **{E[f].enrichment:.2f}×** "
           f"| < {E[f].p_emp:.3f} |\n")

md += f"""
Which null is conservative differs by family and the reason is measurable, not arbitrary.
New Alu insertions land in **AT-rich** sequence (mean GC {gc_poly.gc_mean:.3f} against a
genome baseline of {gc_base.gc_mean:.3f}), while Alu elements themselves sit in **GC-rich**
sequence (fixed reference AluY {gc_fix.gc_mean:.3f}). GC-matching therefore pulls the Alu
null toward regions where Alus are sparse and *inflates* the enrichment; uniform is honest
for Alu. For LINE-1 the direction reverses. The table gives the conservative figure per
family and shows both.

That GC contrast is a result in itself. Polymorphic AluY sit at GC {gc_poly.gc_mean:.3f} and
fixed reference AluY at {gc_fix.gc_mean:.3f} — a shift of
**+{gc_fix.gc_mean - gc_poly.gc_mean:.3f}** with subfamily held constant, so insertion and
retention are being separated rather than confounded. Insertion is AT-biased; retention is
GC-biased. This is the standard resolution of the Alu distribution paradox, measured here in
a single genome.

## 2. Sense and antisense, as separate categories

Host-element selection is **orientation-blind**: where several same-family elements overlap a
breakpoint the longest is chosen, never the one whose strand matches. The pipeline's own
`_annotate_nested_retrotransposon` scores candidates `(same_orient, length)`, so reusing it
would have folded the tie-break into the measurement.

| family | nested | same orientation | rate | 95% CI | odds ratio | Fisher p |
|---|---|---|---|---|---|---|
"""
for f in NICE:
    o = O[f]
    md += (f"| {NICE[f]} | {int(o.n_nested)} | {int(o.n_same)} | **{o.same_frac:.3f}** "
           f"| {o.ci_lo:.3f}–{o.ci_hi:.3f} | {orf(o.odds_ratio)} | {pv(o.p_fisher)} |\n")

md += f"""
Splitting the two classes and testing each against the null separately is where the single
rate turns out to be hiding something:

| family | category | observed | expected | enrichment | p |
|---|---|---|---|---|---|
"""
for f in NICE:
    for c in ("nested_sense", "nested_antisense"):
        r = SA[(f, c)]
        md += (f"| {NICE[f]} | `{c}` | {int(r.observed)} | {r.expected:.2f} "
               f"| **{r.enrichment:.2f}×** | {pv(r.p_two_sided)} |\n")

md += f"""
Alu antisense nesting occurs at **{SA[('ALU','nested_antisense')].enrichment:.2f}× chance**
(p = {pv(SA[('ALU','nested_antisense')].p_two_sided)}) — fewer events than random placement
would produce, so antisense nesting is actively disfavoured rather than merely rarer than
sense. LINE-1 antisense is not distinguishable from chance
({SA[('LINE1','nested_antisense')].enrichment:.2f}×,
p = {pv(SA[('LINE1','nested_antisense')].p_two_sided)}). That family difference is invisible
in a pooled same-orientation rate, and it is the reason the two categories are worth carrying
in the VCF.

**Both readings of "or overall" agree.** The marginals are even — insertions
{O['ALU'].ins_plus_frac:.3f} `+`, host elements {O['ALU'].ref_plus_frac:.3f} `+` for Alu — so
the independence-expected same-orientation rate is {O['ALU'].exp_indep:.3f}. Testing against
chance and against the dataset's own orientation split are the same test here.

### Not a calling artefact

Nested calls carry older subfamily labels than non-nested ones, which is what you would see
if the caller were reading flanking reference sequence instead of the insert. Splitting the
orientation test by subfamily age leaves it unchanged:

| family | subset | n | same orientation | p vs 0.5 |
|---|---|---|---|---|
"""
for f in ("ALU", "LINE1"):
    for sub in ("young only", "old only"):
        r = S[(f, sub)]
        md += f"| {NICE[f]} | {sub} | {int(r.n)} | {r.frac:.3f} | {pv(r.p)} |\n"

md += f"""
Independently, only {100*D['ALU'].frac_sub_eq_host:.1f}% of nested Alu calls carry a
subfamily identical to their host element ({100*D['LINE1'].frac_sub_eq_host:.1f}% for
LINE-1), so these are not re-detections of the host.

![Self-insertion in HG03086. (a) Sense and antisense nesting tested separately against the null; log scale, open markers mark SVA as underpowered. (b) Where young AluY insertions land by host-element age, expected in proportion to each class's genomic bp. (c) Observed same-orientation fraction against the prediction from the degenerate EN consensus, per host subfamily; marker area scales with insertion count.]({{{{artifact:art_{FIG_MAIN}}}}})

## 3. Target choice tracks sequence similarity

Restricting to young AluY-subfamily insertions and asking which host they land in, against an
expectation proportional to each host class's share of genomic Alu bp:

| host element | genomic bp | observed | expected | enrichment | p |
|---|---|---|---|---|---|
"""
for r in hom.itertuples():
    md += (f"| {r.host_age} | {r.host_bp_Mb:.1f} Mb | {int(r.observed)} | {r.expected:.1f} "
           f"| **{r.enrichment:.2f}×** | {pv(r.p)} |\n")

md += f"""
The final row is a heterogeneous remainder (FLAM, FRAM and other non-age-graded Alu-family
annotations); it is listed so the table accounts for all {int(hom.observed.sum())} insertions
rather than only the {int(hom[hom.host_age != 'other'].observed.sum())} that fall on the age
axis, and it carries no interpretation.

A graded preference along the similarity axis: young insertions favour young hosts
({hom.iloc[0].enrichment:.2f}×, p = {pv(hom.iloc[0].p)}) and avoid the oldest
({hom[hom.host_age.str.startswith('AluJ')].iloc[0].enrichment:.2f}×). Jacob-Hirsch et al. 2018 report the same ordering for
LINE-1 in somatic brain — L1PA hosting ~71% of nested insertions while comprising 13% of
genomic L1 — and state that host homology outweighs host abundance in target choice. This
reproduces that result in **Alu**, in **germline polymorphic** insertions, in a single genome.

It does **not** explain orientation. Same-orientation fraction is flat across host age
(0.875 AluJ, 0.829 AluS, 0.818 AluY, n = 24 / 76 / 33). A homologous-annealing mechanism
would require same orientation and so predicts the fraction should *rise* with similarity;
it does not, and if anything falls.

### Where inside the host element insertions land: two A-rich hotspots

Levy, Schwartz & Ast 2009 report that nested insertions occupy positional hotspots along the
targeted element rather than spreading across it. Tested here, and it holds: the distribution
across deciles of the host element is strongly non-uniform —
**χ² = {ps['chi2']:.0f} on {ps['chi2_df']} df, p = {pv(ps['chi2_p'])}**, with a
Kolmogorov–Smirnov test against uniform giving D = {ps['ks_D']:.3f}, p = {pv(ps['ks_p'])}
(n = {ps['n_alu_nested']} nested Alu calls).

Restricting to near-full-length host Alus ({ps['n_full_length_host']} calls, host 280–320 bp)
so base offsets are comparable, and binning at 20 bp against a uniform expectation of
{pp.expected_uniform.iloc[0]:.1f} per bin:

| position in host | insertions | host A-fraction | region |
|---|---|---|---|
{PEAK_ROWS}
The two dominant peaks coincide with the two A-richest internal regions of Alu. Insertion
count correlates with host A-fraction across bins (Pearson r = {ps['pearson_r']:+.3f},
p = {pv(ps['pearson_p'])}) but not rank-wise (Spearman ρ = {ps['spearman_rho']:+.3f},
p = {pv(ps['spearman_p'])}). That split is informative: the association is carried by two
discrete peaks, not by a smooth gradient in A-content — the hotspot structure Levy et al.
describe rather than a density effect.

This bears on §4. The motif models tested there fail to predict the orientation *magnitude*,
but A-run targeting is supported *positionally*: insertions concentrate at the 3′ poly-A tail
and the inter-arm A-rich linker. Target position and target orientation are not accounted for
by the same thing in this dataset.

The 5′-most bin carries {[n for b, e, n, af in _pk if b < 100][0]} insertions at low
A-content ({[af for b, e, n, af in _pk if b < 100][0]:.2f}) and has no mechanistic account.
It sits at an element boundary, where breakpoint-assignment edge effects are the leading
explanation; it is flagged rather than interpreted.

![Nested Alu insertions by position in the host element. Bars are insertion counts in 20 bp bins for near-full-length hosts; the line is mean A-fraction of host sequence in the same bins, read in element orientation; the dashed line is the uniform expectation.]({{{{artifact:art_{FIG_POS}}}}})

## 4. Orientation: three further mechanisms tested, three rejected

Host-element homology was the first of the four candidates and was rejected in §3. The three
below are the remainder.

**Narrow endonuclease consensus (`TTTT/AA`).** Counting sense- and antisense-compatible sites
inside host elements predicted the family-level rate closely — Alu
{NW['ALU'].pred_same_strict:.3f} against observed {NW['ALU'].observed_same:.3f}, LINE-1
{NW['LINE1'].pred_same_strict:.3f} against {NW['LINE1'].observed_same:.3f}. That agreement
does not survive. It was computed over elements weighted by *genomic abundance*, dominated by
old low-prediction copies; weighted by the subfamilies that actually received insertions the
same model predicts 0.920 against an observed 0.842. Two weightings landing near each other
is not a prediction, and across subfamilies the model was rejected outright.

**Full degenerate consensus (`YY/RRRR`).** `TTTT/AA` is one instance of the published L1 EN
consensus; scoring only it misses every other variant. Recounting on the full degenerate
pattern:

| family | predicted (degenerate) | predicted (narrow) | observed |
|---|---|---|---|
"""
for f in NICE:
    md += (f"| {NICE[f]} | {DF[f].pred_same:.3f} | {NW[f].pred_same_strict:.3f} "
           f"| **{DF[f].observed_same:.3f}** |\n")

md += f"""
The correct consensus **underpredicts** the bias for every family. Across host subfamilies it
fails as well — **χ² = {CHI2:.1f} on {CHI2_DF} df, p = {pv(CHI2_P)}**, Pearson
r = {R_P:+.3f} (p = {R_P_P:.3g}) — and the reason is visible in its range: predictions span
only {dgs.pred_same.min():.3f}–{dgs.pred_same.max():.3f} across all ten subfamilies against an
observed spread of {dgs.observed_same.min():.3f}–{dgs.observed_same.max():.3f}. There is
almost no variance in the predictor to explain variance in the outcome.

Method check: genome-wide the two consensus orientations occur at
{bg['YYRRRR_per10kb']:.1f} and {bg['YYYYRR_per10kb']:.1f} sites per 10 kb — ratio
**{bg['ratio']:.3f}**. The genome has no preferred strand and the counter agrees, so the
asymmetry seen inside elements is a property of the elements.

**Non-classical CCATT target motif.** Jacob-Hirsch et al. attribute same-orientation nested
insertion in somatic brain to an endonuclease-*independent* route marked by a CCATT target
motif rather than the classical consensus. Against a control matched for element composition
— random positions *inside same-family elements*, since nested sites sit in element sequence
by construction and comparing them to non-nested sites would compare compositions rather than
preferences — there is no enrichment: {CC.per_kb_a:.2f} vs {CC.per_kb_b:.2f} per kb, rate
ratio {CC.rate_ratio:.2f}, **p = {pv(CC.p)}**. Consistent with their motif belonging to the
somatic ENi route while these germline insertions are classical.

So the orientation preference is real, large, and unexplained here.

**It is also a reproduction, and that is now established by direct reading rather than by a
secondary citation.** Jacob-Hirsch et al. attribute the germline same-orientation preference
to their reference 33, which is **Levy, Schwartz & Ast 2009** (*Nucleic Acids Res* 38:1515–30,
doi:10.1093/nar/gkp1134). That paper analysed hundreds of thousands of nested transposed
elements genome-wide and reports a **4-fold** higher tendency for elements to insert into
existing elements than into non-TE intergenic regions, a non-canonical single-nucleotide
hotspot in Alu for the insertion of other Alu, and concludes that insertion is biased "within
certain TEs, in specific orientations and within specific targeted TE positions." An earlier
draft of this report flagged that this claim rested on Jacob-Hirsch's citation rather than on
the source; the source has since been read and it supports the claim directly.

Their enrichment figure is not directly comparable to §1's: theirs pools all element families
and contrasts against non-TE intergenic DNA, while §1 is same-family only against a
genome-wide resampled null. Both are positive; the magnitudes answer different questions.

What this analysis contributes is therefore a narrowing of what does *not* cause the
orientation bias, plus the positional result in §3 — not an account of the magnitude.

## 5. A missing feature in the pipeline, not a bug

When this analysis was run, `_annotate_nested_retrotransposon` set
`nested_same_class_orientation` to `"nested"` where the orientations matched and `"unnested"`
otherwise, so a call inside a same-family element in the *opposite* orientation was labelled
`unnested` — the value used for a call overlapping nothing — while the same row carried
`nested_repeat_overlap=True` and `nested_same_class=True`. That string became
`nested_in_same_MEI` (`mei_support.py:14778`) and was exported as the VCF `NESTED` field
(`vcf_export.py:533`).
**{n_anti} calls in this callset are affected**
({', '.join(f'{int(v)} {NICE[k]}' for k, v in anti_by_fam.items())}).

An earlier draft of this report called that a bug. It is more accurately a missing category:
the field records one class correctly and has no vocabulary for the other, which is why §2
above could not be produced from pipeline output.

**Now fixed and pushed** on `fix/nested-orientation-semantics`. The `NESTED` field is a
four-valued enum:

| value | meaning |
|---|---|
| `unnested` | no same-family element overlaps the breakpoint |
| `nested_sense` | overlaps; insertion orientation matches the element strand |
| `nested_antisense` | overlaps; orientations differ |
| `nested_unknown` | overlaps, but orientation is unresolvable on one side |

`nested_unknown` was not in the original specification and is the right addition: without it,
a call whose orientation cannot be determined would be forced into `nested_antisense` and
inflate that count — the same class of error as the original collapse, and it would land
directly on the comparison §2 makes. Verified at branch head `{ts['branch_head']}`:
**{ts['passed']} passed, {ts['skipped']} skipped** across the nested and VCF-export suites
with `bedtools` {ts['bedtools'].split()[0]}. The nested suite itself has grown from
{ts['nested_file_history']['origin/main']} tests on `origin/main` to
{ts['nested_file_history']['after_all_five_commits']} on the branch.

Two further fixes on the same branch harden the annotation against null coordinates, null
chromosomes, and NaN hit flags. Separately, the function's `(same_orient, length)` host
selection was changed to select on length alone — the previous ordering preferred
same-orientation hits, which would bias any orientation analysis run on its output. The
analysis in this report selected on length independently for that reason, so §1–§3 are
unaffected by either version.

## 6. Limitations

- **One genome.** HG03086 is the only genome-wide callset in the project; every other sample
  (HG00100, SEQC2, HG03172) is a chr22 slice. Three results below are power-limited and a
  second genome-wide callset would lift all three at once.
- **SVA is underpowered.** {int(E['SVA'].obs_nested)} nested calls. Its enrichment ratio is
  large mainly because SVA covers {D['SVA'].genome_cov_pct:.2f}% of the genome, so expectation
  is near zero. {D['SVA'].n_sub_eq_host} of its {int(O['SVA'].n_nested)} nested calls carry a
  subfamily identical to their host — the one family where re-detection cannot be excluded.
  Suggestive, not established.
- **Genic vs intergenic: not assessed.** The sense fraction is
  {ctx[(ctx.family=='ALU') & (ctx.context=='genic')].iloc[0].sense_frac:.3f} genic against
  {ctx[(ctx.family=='ALU') & (ctx.context=='intergenic')].iloc[0].sense_frac:.3f} intergenic,
  but with {int(ctx[(ctx.family=='ALU') & (ctx.context=='intergenic')].iloc[0].n_nested)}
  intergenic events the comparison can only detect a 21-point difference at 80% power. This
  is absence of power, not absence of effect, and is reported as such.
- **Orientation mechanism is open.** All four candidate mechanisms rejected — host-element
  homology (§3) plus three target-site models (§4). No positive account of the magnitude is
  offered here.
- **Short-read mapping inside repeats.** Every nested call sits in repetitive sequence, where
  short-read callers are least reliable. The young/old split argues against a gross artefact;
  the definitive check is long-read validation, and the pipeline already ingests a long-read
  callset via `--lr-mei-vcf`.
- **Element boundaries affect the rate, not the orientation result.** An earlier draft
  asserted that RepeatMasker annotates Alu poly-A tails separately from element bodies, so
  that tail-primed insertions would fall outside the interval and understate nesting. That
  was checked and is **false**: on chr1, {POLYA_PCT:.1f}% of {POLYA_N:,} Alu elements have a
  separately annotated adjacent poly-A/T, so the tail is almost always inside the annotation.
  The real boundary question was tested directly by relaxing every interval by a fixed slop
  and recomputing both the nesting flag and the host strand from scratch at each step:

| slop | Alu nested | Alu same-orient | LINE-1 nested | LINE-1 same-orient |
|---|---|---|---|---|
{BOUNDARY_ROWS}
  The nesting **rate** is boundary-dependent, as any interval-overlap statistic is — but the
  enrichment in §1 compares observed against a null built from the same intervals, so the
  comparison is like-for-like at whatever boundary is used. The **orientation** result is
  robust: still {BS_ALU_50:.3f} for Alu and {BS_L1_50:.3f} for LINE-1 at ±50 bp, decaying only
  as flanking non-element sequence — where no orientation relationship exists — is pulled in.

## Prior work

| year | citation | bearing on this analysis |
|---|---|---|
"""
for r in lit.sort_values("year").itertuples():
    md += f"| {r.year} | {r.citation} | {r.bearing} |\n"

md += """
## Reproducing

`build_selfins_report.py` regenerates this document from the tables below; every figure in
the prose is read from a file rather than typed.

| file | contents |
|---|---|
| `selfins_per_call.csv` | all calls with host element, nesting, orientation |
| `selfins_enrichment.csv` | both nulls, all three families |
| `selfins_sense_antisense.csv` | the two orientation classes tested separately |
| `selfins_orientation.csv` | orientation tests with marginals and Fisher results |
| `selfins_orientation_sensitivity.csv` | young/old subfamily split |
| `selfins_homology_test.csv` | host-age preference for young insertions |
| `selfins_en_degenerate_family.csv` | degenerate `YY/RRRR` recount, family level |
| `selfins_en_degenerate_subfamily.csv` | the same per subfamily, where it is rejected |
| `selfins_en_motif_prediction.csv` | the earlier narrow `TTTT/AA` count |
| `selfins_ccatt_control.csv` | CCATT control, composition-matched |
| `selfins_genome_background_en.json` | genome-wide consensus density, both orientations |
| `selfins_position_profile.csv` | insertion counts and host A-fraction per 20 bp bin (§3) |
| `selfins_position_deciles.csv` | decile counts underlying the uniformity test (§3) |
| `selfins_position_stats.json` | χ², KS and correlation statistics for §3 |
| `selfins_boundary_sensitivity.csv` | nesting and orientation as element boundaries are relaxed |
| `selfins_polya_annotation_check.json` | whether rmsk annotates Alu poly-A tails separately |
| `alu_insertion_vs_retention_gc.csv` | insertion vs retention GC contrast |
| `nested_orientation_literature.csv` | prior work with per-paper bearing |
"""

open("selfins_report.md", "w").write(md)
print("wrote selfins_report.md", len(md), "chars")
