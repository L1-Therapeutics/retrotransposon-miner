"""Regenerate the self-insertion report from the saved result tables.

Core quantitative summaries are rendered from CSV/JSON outputs. Literature context and
interpretation are curated in this template and should be reviewed when upstream results change.
"""
import json
import re
import sys
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.platypus import Image, Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle

import numpy as np
import pandas as pd
from scipy import stats

BASE = Path(__file__).resolve().parent
DATA_DIR = BASE / "notes" if (BASE / "notes").is_dir() else BASE
sys.path.insert(0, str(BASE))
from retraction_guard import validate_sources  # noqa: E402

REPORT_PATH = BASE / "selfins_report.md"
PDF_PATH = BASE / "selfins_report.pdf"

FIG_MAIN = "eb5945f3-3fc5-4ae1-b65d-2090eea2defb"   # selfins_main.png

enr   = pd.read_csv(DATA_DIR / "selfins_enrichment.csv")
sa    = pd.read_csv(DATA_DIR / "selfins_sense_antisense.csv")
ori   = pd.read_csv(DATA_DIR / "selfins_orientation.csv")
sens  = pd.read_csv(DATA_DIR / "selfins_orientation_sensitivity.csv")
hom   = pd.read_csv(DATA_DIR / "selfins_homology_test.csv")
dgf   = pd.read_csv(DATA_DIR / "selfins_en_degenerate_family.csv")
dgs   = pd.read_csv(DATA_DIR / "selfins_en_degenerate_subfamily.csv")
narrow= pd.read_csv(DATA_DIR / "selfins_en_motif_prediction.csv")
ccatt = pd.read_csv(DATA_DIR / "selfins_ccatt_control.csv")
diag  = pd.read_csv(DATA_DIR / "selfins_diagnostics.csv")
gc    = pd.read_csv(DATA_DIR / "alu_insertion_vs_retention_gc.csv")
ctx   = pd.read_csv(DATA_DIR / "selfins_sense_antisense_by_context.csv")
per   = pd.read_csv(DATA_DIR / "selfins_per_call.csv")
lit   = pd.read_csv(DATA_DIR / "nested_orientation_literature.csv")
bg    = json.loads((DATA_DIR / "selfins_genome_background_en.json").read_text())
bs    = pd.read_csv(DATA_DIR / "selfins_boundary_sensitivity.csv")
pa    = json.loads((DATA_DIR / "selfins_polya_annotation_check.json").read_text())
pp    = pd.read_csv(DATA_DIR / "selfins_position_profile.csv")
ps    = json.loads((DATA_DIR / "selfins_position_stats.json").read_text())
ts    = json.loads((DATA_DIR / "selfins_branch_test_status.json").read_text())

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

**Answers.** Same-family nesting exceeds the selected placement nulls for Alu and LINE-1.
For SVA the result is 9 observed vs 0.306 expected (n=9; repeat-aware callability untested),
so it remains provisional. Nummi et al. 2025 report SVA self-target preference, providing
systematic precedent but not validation of these nine calls. The nested set has a strong
same-orientation tendency for Alu and LINE-1; within-nested tests and category-vs-placement
tests are distinct. Among AluY insertions, young-host targeting is above an abundance-based
expectation (1.86×), an age-associated pattern compatible with—but not proof of—homology-
mediated targeting. The orientation tendency remains unexplained by these analyses: tested
motif/opportunity models do not explain its magnitude, and no enrichment was detected for
CCATT under the selected control. These results do not assign or exclude a route.

**Scope note.** This is the table-driven report and retains exploratory subfamily and motif
sections. See the [integrated evidence report](NESTED_MEI_INTEGRATED_REPORT_2026-10-02.md)
for the current literature-reconciled synthesis and the [validation input guide](SELFINS_VALIDATION_INPUT_GUIDE.md)
for the blocked sequence/read-backed checks. The figures, tables, and quantitative summaries below
are rendered into both this Markdown and PDF by this script.

Callset: {n_total} classifier calls (score ≥ 0.997) on HG03086, GRCh38, chr1–22 + chrX —
one genome, which bounds several results below. Reference elements: UCSC hg38 RepeatMasker,
filtered to Alu / LINE-1 / SVA with the pipeline's own family-normalisation rule
(`_normalize_mei_family_token`). {n_nested} of {n_total} calls
({100*n_nested/n_total:.1f}%) fall inside a same-family element.

---

## 1. Same-family nesting against family-specific selected nulls

Two nulls, each resampling call positions genome-wide 1,000 times — **uniform**, and
**GC-matched** on each call's local GC (1 kb window). Empirical p-values at 0.000999 are
at the simulation floor, not precise tail-probability estimates.

| family | n | nested | observed | GC-matched exp. | uniform exp. | conservative enrichment | p |
|---|---|---|---|---|---|---|---|
"""
for f in NICE:
    g, u = EA[(f, "GC-matched")], EA[(f, "uniform")]
    md += (f"| {NICE[f]} | {int(E[f].n)} | {int(E[f].obs_nested)} | {E[f].obs_frac:.3f} "
           f"| {g.exp_frac:.3f} | {u.exp_frac:.3f} | **{E[f].enrichment:.2f}×** "
           f"| < {E[f].p_emp:.3f} |\n")

md += f"""
The selected null differs by family. GC matching yields an Alu expectation of
{EA[('ALU','GC-matched')].exp_nested:.1f} vs {EA[('ALU','uniform')].exp_nested:.1f} under uniform
({EA[('ALU','GC-matched')].enrichment:.2f}× vs {EA[('ALU','uniform')].enrichment:.2f}×), and raises expectations for LINE-1 and SVA relative to uniform. These are model-dependent
placement comparisons; neither null is proven universally correct. Polymorphic AluY in this callset have lower mean GC than fixed reference AluY ({gc_poly.gc_mean:.3f} vs {gc_fix.gc_mean:.3f}).
This descriptive contrast does not by itself separate integration preference from sequence
composition, ascertainment, or survival/selection.

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
Against the orientation-blind placement simulations, Alu antisense nesting is below the
selected null ({SA[('ALU','nested_antisense')].enrichment:.2f}×, p = {pv(SA[('ALU','nested_antisense')].p_two_sided)}), whereas LINE-1 antisense is compatible with its null
({SA[('LINE1','nested_antisense')].enrichment:.2f}×, p = {pv(SA[('LINE1','nested_antisense')].p_two_sided)}). These are category-vs-placement tests, not a direct sense-vs-antisense
contrast among nested calls. The conditional same-orientation fractions and their binomial
tests are reported separately above; Alu's marginal-based expectation ({O['ALU'].exp_indep:.3f}) is close to,
but not conceptually identical with, the 0.5 binomial baseline.

**The chance baseline is close to the marginal-based expectation, but the tests are distinct.** For Alu, new insertions are {O['ALU'].ins_plus_frac:.3f} `+` and reference hosts are {O['ALU'].ref_plus_frac:.3f} `+`, giving an independence-expected same-orientation rate of {O['ALU'].exp_indep:.3f}. The saved conditional exact-binomial test uses 0.5; the saved Fisher test uses the actual orientation marginals. These are similar here, not the same model in general.

### A limited subfamily negative control

Nested calls carry different subfamily-age distributions from non-nested calls, which argues
against one simple host-sequence misassignment explanation. It is not a general exclusion of
calling or mapping artefacts; the orientation split by age is a robustness check, not read-level
validation:

| family | subset | n | same orientation | p vs 0.5 |
|---|---|---|---|---|
"""
for f in ("ALU", "LINE1"):
    for sub in ("young only", "old only"):
        r = S[(f, sub)]
        md += f"| {NICE[f]} | {sub} | {int(r.n)} | {r.frac:.3f} | {pv(r.p)} |\n"

md += f"""
Only {100*D['ALU'].frac_sub_eq_host:.1f}% of nested Alu calls carry a subfamily label identical to their host element
({100*D['LINE1'].frac_sub_eq_host:.1f}% for LINE-1). This argues against one simple same-label host re-detection explanation;
it is not a general exclusion of calling or mapping artefacts and does not independently
validate the calls.

![Figure 1. Self-insertion in HG03086: (a) sense and antisense nesting tested separately against the null; log scale, open markers mark SVA as underpowered. (b) where young AluY insertions land by host-element age, expected in proportion to each class's genomic bp. (c) observed same-orientation fraction against the degenerate EN consensus, per host subfamily; marker area scales with insertion count.](selfins_main.png)

## 3. AluY host-age pattern (sequence similarity is not isolated)

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

The abundance-based table shows an age-associated distribution: young AluY insertions favour young hosts ({hom.iloc[0].enrichment:.2f}×, p = {pv(hom.iloc[0].p)}), while the oldest AluJ class is below expectation ({hom[hom.host_age.str.startswith('AluJ')].iloc[0].enrichment:.2f}×, p = {pv(hom[hom.host_age.str.startswith('AluJ')].iloc[0].p)}). Host age is only a proxy for sequence similarity and also covaries with composition, genomic location, and survival; this result does not by itself demonstrate homology-directed integration. Jacob-Hirsch et al. 2018 describe homology/abundance effects for a distinct somatic-brain LINE-1 class. The direction is qualitatively compatible, but the cohorts and mechanisms differ, so this is not a direct replication.

This age-associated pattern does not explain orientation: same-orientation fractions are
similar across host-age strata (0.875 AluJ, 0.829 AluS, 0.818 AluY, n = 24 / 76 / 33).
A simple model in which host homology drives the orientation preference would predict a
stronger same-orientation fraction in more similar hosts; these data do not show that trend.

### Where inside the host element insertions land: two A-rich peaks

Levy, Schwartz & Ast (2009 online; 2010 issue, *Nucleic Acids Research* 38:1515–1530, doi:10.1093/nar/gkp1134; Results p. 1518, Figure 2) report a prominent Alu-into-Alu insertion hotspot directly after the A-rich linker at targeted-Alu consensus coordinate 133. Their Figure 2 describes the AluJo consensus and labels the internal A-rich linker at 118–136; Figure 1 defines insertion position as the coordinate just upstream of the insertion. Under a 5′-aligned consensus-to-host convention, coordinate 133 falls in our 120–140 bp bin, with the breakpoint around the 133/134 boundary—not an exact bin edge. The paper does not explicitly say that 133 is counted from the first consensus nucleotide, so that mapping is an inference. Figure 2's sense-oriented Alu self-insertion display uses 15,759 Alu self-insertion events. Our 20-bp profile is positionally compatible, not a nucleotide-resolution replication.

The HG03086 distribution across deciles of the host element is strongly non-uniform —
**χ² = {ps['chi2']:.0f} on {ps['chi2_df']} df, p = {pv(ps['chi2_p'])}**, with a
Kolmogorov–Smirnov test against uniform giving D = {ps['ks_D']:.3f}, p = {pv(ps['ks_p'])}
(n = {ps['n_alu_nested']} nested Alu calls).

Restricting to near-full-length host Alus ({ps['n_full_length_host']} calls, host 280–320 bp)
so base offsets are comparable, and binning at 20 bp against a uniform expectation of
{pp.expected_uniform.iloc[0]:.1f} per bin:

| position in host | insertions | host A-fraction | region |
|---|---|---|---|
{PEAK_ROWS}
The separate 280–300 bp peak (36, 3.7×) is near the terminal poly(A) tail under this host convention. Levy et al. excluded target-element consensus ends from positional tests because of RepeatMasker/short-flank identification limitations, so their data do not directly adjudicate our terminal-tail cluster. This is a coverage blind spot, not evidence by itself that the cluster is novel or genuine. Nummi et al. 2025 independently report germline Alu target peaks at 100–155 and 260–310 bp, broadly covering both of our peak windows, but their long-read cohort and positional method differ.

Insertion count correlates with host A-fraction across bins (Pearson r = {ps['pearson_r']:+.3f}, p = {pv(ps['pearson_p'])}) but not rank-wise (Spearman ρ = {ps['spearman_rho']:+.3f}, p = {pv(ps['spearman_p'])}): the association is carried by discrete peaks, not a smooth gradient. Neither positional peak validates its individual insertions as authentic TPRT events.

This bears on §4. The motif models tested there fail to predict the orientation *magnitude*. A-rich target positions coincide with the two larger positional peaks, but this positional association does not explain the orientation result.

The 5′-most bin carries {[n for b, e, n, af in _pk if b < 100][0]} insertions at low
A-content ({[af for b, e, n, af in _pk if b < 100][0]:.2f}) and has no mechanistic account.
It sits at an element boundary, where breakpoint-assignment edge effects are the leading
explanation; it is flagged rather than interpreted.

![Figure 2. Nested Alu insertion positions within the host element. Bars show insertion counts in 20 bp bins for near-full-length hosts; line shows mean host A-fraction in element orientation; dashed line is the uniform expectation.](selfins_position.png)

## 4. Orientation: tested models do not explain the observed tendency

The age-associated host-choice result in §3 does not account for orientation. Three
sequence/opportunity hypotheses are evaluated below; none explains the observed magnitude.

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
The tested full degenerate consensus predicts a lower same-orientation fraction than
observed for each family. Across host subfamilies its predictions have little range relative
to the observed fractions. This model does not account for the observed magnitude; it does
not establish that endonuclease opportunity is biologically irrelevant.

Method check: genome-wide the two consensus orientations occur at
{bg['YYRRRR_per10kb']:.1f} and {bg['YYYYRR_per10kb']:.1f} sites per 10 kb — ratio **{bg['ratio']:.3f}**. This background check finds no global
strand skew in the motif counts; it does not by itself explain the nested-call orientation
pattern.

**CCATT target motif (literature context, not route assignment).** Jacob-Hirsch et al. report a
CCATT-associated EN-independent route for a specific somatic-brain nested L1 class. Against a
control matched for element composition—random positions *inside same-family elements*, since
nested sites sit in element sequence by construction—no enrichment was detected in this
callset under this control: {CC.per_kb_a:.2f} vs {CC.per_kb_b:.2f} per kb, rate ratio {CC.rate_ratio:.2f}, **p = {pv(CC.p)}**. This does not
establish that these calls used the classical route; route assignment needs independent
junction, truncation, and target-motif evidence.

The observed orientation tendency is pronounced in these nested calls, but remains unexplained by the tested sequence/opportunity models. These tests do not establish that the corresponding biological mechanisms are absent.

**Levy et al.'s “~4-fold” figure has a different unit than this report's enrichment.** The “~4-fold” figure from Levy et al. 2010 is not the ratio of total insertions into TE sequence to insertions into non-TE sequence. It is the ratio of enriched to depleted insertion-type categories: 516 overrepresented versus 116 underrepresented insertion-type categories out of 36,478 tested (young × old × orientation), i.e. ~4.45:1 by category. Their per-type null is random insertion into unique intergenic DNA scaled to each old TE class’s inferred ancestral coverage. Any total-event ratio must be computed from the underlying event counts, which the paper does not report as a single number.

That category-count result is not directly comparable to §1's family-specific observed-vs-resampled event counts; the two answer different statistical questions. Levy et al.'s exact-insertion database is a broader set than the young-into-old unique-intergenic analysis.

These tests identify hypotheses not supported by the current measurements; they do not rule
out the corresponding biological mechanisms. The orientation tendency remains unexplained by
this analysis.

## 5. Pipeline annotation context (separate from validation of this callset)

Earlier pipeline output collapsed same-family opposite-orientation overlaps into the `unnested`
label, so those records could not be separated from genuinely non-nested calls using that field
alone. In the analysis-table reconstruction, {n_anti} same-family opposite-orientation calls were identified ({', '.join(f'{int(v)} {NICE[k]}' for k, v in anti_by_fam.items())}). This describes a pipeline-annotation limitation and is not
read-level validation of the current callset.

The analysis reconstructed the two orientation classes independently of that collapsed field.
A separate pipeline branch documented in the prior project handoff added a four-valued `NESTED`
enum (branch/test details are not independently reverified in this report):

| value | meaning |
|---|---|
| `unnested` | no same-family element overlaps the breakpoint |
| `nested_sense` | overlaps; insertion orientation matches the element strand |
| `nested_antisense` | overlaps; orientations differ |
| `nested_unknown` | overlaps, but orientation is unresolvable on one side |

`nested_unknown` retains records whose host or insertion orientation cannot be resolved.
For this report's analysis, host selection was reconstructed without preferring matching
strands; raw reads were not independently validated.

## 6. Limitations

- **One genome.** HG03086 is the only genome-wide callset in the project; other samples are
  chr22 slices. These are catalogue-pattern estimates for one individual, not population-wide
  integration rates.
- **SVA is provisional.** {int(E['SVA'].obs_nested)} nested calls ({int(E['SVA'].obs_nested)} observed vs {E['SVA'].exp_nested:.3f} expected under the selected
  GC-matched null). Repeat-aware callability/mappability was not modeled, so the 29.4× ratio
  is not a settled effect estimate. Nummi et al. 2025 report systematic SVA self-target
  preference, which is precedent—not validation of these nine calls.
- **Genic vs intergenic: not assessed.** The sense fraction is
  {ctx[(ctx.family=='ALU') & (ctx.context=='genic')].iloc[0].sense_frac:.3f} genic against
  {ctx[(ctx.family=='ALU') & (ctx.context=='intergenic')].iloc[0].sense_frac:.3f} intergenic,
  but with {int(ctx[(ctx.family=='ALU') & (ctx.context=='intergenic')].iloc[0].n_nested)}
  intergenic events the comparison can only detect a 21-point difference at 80% power. This
  is absence of power, not absence of effect, and is reported as such.
- **Orientation mechanism is open.** The age-associated host-choice result does not explain
  orientation. Tested motif/opportunity models do not account for its magnitude; this does not
  prove the corresponding biological mechanisms absent.
- **Read-backed validation is incomplete.** Repetitive loci challenge short-read evidence.
  The current checkout lacks HG03086 BAM/CRAM and the full reference sequence needed to
  reconstruct junctions; known-source and read-support annotations do not substitute for
  independent per-locus validation.
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
    display_year = "2009/2010" if str(r.year) == "2009" else str(r.year)
    md += f"| {display_year} | {r.citation} | {r.bearing} |\n"

md += """
## Reproducing

Run `python artifacts/build_selfins_report.py` from the repository root to regenerate this
Markdown and `selfins_report.pdf` from the same report text, tables, and figures. Quantitative
summaries are read from saved results; literature context and interpretation are curated in
the generator template. PDF rendering requires ReportLab (`python -m pip install -r
artifacts/requirements-pdf.txt`).

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

validate_sources(
    {
        REPORT_PATH: md,
        DATA_DIR / "nested_orientation_literature.csv": (DATA_DIR / "nested_orientation_literature.csv").read_text(),
        DATA_DIR / "alu_mechanism_literature.csv": (DATA_DIR / "alu_mechanism_literature.csv").read_text(),
        DATA_DIR / "panel.py": (DATA_DIR / "panel.py").read_text(),
        DATA_DIR / "l1_locus_catalog_chr22.csv": (DATA_DIR / "l1_locus_catalog_chr22.csv").read_text(),
    }
)
REPORT_PATH.write_text(md, encoding="utf-8")


def _markdown_structure(markdown: str) -> dict[str, object]:
    """Summarize visible Markdown structure for a parity check against PDF flowables."""
    lines = markdown.splitlines()
    sections = [line[3:].strip() for line in lines if line.startswith("## ")]
    figure_count = sum(line.lstrip().startswith("![") for line in lines)
    table_rows = []
    index = 0
    while index < len(lines):
        if not lines[index].strip().startswith("|"):
            index += 1
            continue
        rows = []
        while index < len(lines) and lines[index].strip().startswith("|"):
            cells = [cell.strip() for cell in lines[index].strip().strip("|").split("|")]
            if not all(re.fullmatch(r":?-{3,}:?", cell or "-") for cell in cells):
                rows.append(cells)
            index += 1
        table_rows.append(len(rows))
    return {"sections": sections, "figures": figure_count, "table_rows": table_rows}


def render_pdf(markdown: str, pdf_path: Path) -> dict[str, object]:
    """Render Markdown and verify its sections, figures, and tables remain represented."""
    expected = _markdown_structure(markdown)
    rendered_sections = []
    rendered_figures = 0
    rendered_table_rows = []
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="ReportTitle", parent=styles["Title"], fontSize=20, leading=24, spaceAfter=14))
    styles.add(ParagraphStyle(name="Section", parent=styles["Heading1"], fontSize=15, leading=18, spaceBefore=14, spaceAfter=8))
    styles.add(ParagraphStyle(name="Subsection", parent=styles["Heading2"], fontSize=11, leading=14, spaceBefore=10, spaceAfter=6))
    styles.add(ParagraphStyle(name="BodySmall", parent=styles["BodyText"], fontSize=8.5, leading=11, spaceAfter=5))
    styles.add(ParagraphStyle(name="TableCell", parent=styles["BodyText"], fontSize=7, leading=9))
    doc = SimpleDocTemplate(
        str(pdf_path), pagesize=letter, rightMargin=0.62 * inch, leftMargin=0.62 * inch,
        topMargin=0.58 * inch, bottomMargin=0.62 * inch, title="Retrotransposon self-insertion in HG03086",
        author="", subject="Synchronized Markdown report generated from saved analysis tables",
    )
    story = []
    lines = markdown.splitlines()
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            story.append(Spacer(1, 3))
            i += 1
            continue
        if line == "---":
            story.append(Spacer(1, 5))
            i += 1
            continue
        if line.startswith("!["):
            rendered_figures += 1
            match = re.match(r"!\[(.*?)\]\((.*?)\)", line)
            if match:
                image_path = (DATA_DIR / match.group(2)).resolve()
                if not image_path.is_relative_to(DATA_DIR.resolve()):
                    raise ValueError(f"figure path must stay under {DATA_DIR}: {match.group(2)}")
                if image_path.exists():
                    image = Image(str(image_path), width=6.8 * inch, height=3.5 * inch, kind="proportional", hAlign="CENTER")
                    caption = Paragraph(_inline_format(match.group(1)), styles["BodySmall"])
                    story.extend([Spacer(1, 5), image, caption, Spacer(1, 8)])
                else:
                    raise FileNotFoundError(f"missing report figure: {image_path}")
            i += 1
            continue
        if line.startswith("|"):
            rows = []
            while i < len(lines) and lines[i].strip().startswith("|"):
                row = lines[i].strip().strip("|")
                cells = [cell.strip() for cell in row.split("|")]
                if not all(re.fullmatch(r":?-{3,}:?", cell or "-") for cell in cells):
                    rows.append([Paragraph(_inline_format(cell), styles["TableCell"]) for cell in cells])
                i += 1
            if rows:
                rendered_table_rows.append(len(rows))
                ncols = max(len(row) for row in rows)
                rows = [row + [Paragraph("", styles["TableCell"])] * (ncols - len(row)) for row in rows]
                widths = [6.76 * inch / ncols] * ncols
                table = Table(rows, colWidths=widths, repeatRows=1, hAlign="LEFT")
                table.setStyle(TableStyle([
                    ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#e7edf3")),
                    ("TEXTCOLOR", (0, 0), (-1, 0), colors.HexColor("#17324d")),
                    ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#9aa8b5")),
                    ("VALIGN", (0, 0), (-1, -1), "TOP"),
                    ("LEFTPADDING", (0, 0), (-1, -1), 4), ("RIGHTPADDING", (0, 0), (-1, -1), 4),
                    ("TOPPADDING", (0, 0), (-1, -1), 3), ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
                ]))
                story.extend([table, Spacer(1, 7)])
            continue
        if line.startswith("# "):
            story.append(Paragraph(_inline_format(line[2:]), styles["ReportTitle"]))
        elif line.startswith("## "):
            rendered_sections.append(line[3:].strip())
            story.append(Paragraph(_inline_format(line[3:]), styles["Section"]))
        elif line.startswith("### "):
            story.append(Paragraph(_inline_format(line[4:]), styles["Subsection"]))
        else:
            paragraph_lines = [line]
            i += 1
            while i < len(lines) and lines[i].strip() and not lines[i].lstrip().startswith(("#", "|", "![", "---")):
                paragraph_lines.append(lines[i].strip())
                i += 1
            story.append(Paragraph(_inline_format(" ".join(paragraph_lines)), styles["BodySmall"]))
            continue
        i += 1
    actual = {"sections": rendered_sections, "figures": rendered_figures, "table_rows": rendered_table_rows}
    if actual != expected:
        raise RuntimeError(f"Markdown/PDF structure mismatch: Markdown={expected!r}; PDF flowables={actual!r}")
    doc.build(story, onFirstPage=_page_footer, onLaterPages=_page_footer)
    return {**actual, "pages": doc.page}


def _xml_escape(text: str) -> str:
    from xml.sax.saxutils import escape
    return escape(text)


def _inline_format(text: str) -> str:
    """Escape report text and translate a small, safe subset of Markdown inline syntax."""
    from xml.sax.saxutils import escape

    tokens: list[str] = []

    def preserve_tag(match: re.Match[str]) -> str:
        tokens.append(match.group(0))
        return f"REPORTMARKUPTOKEN{len(tokens) - 1}END"

    text = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", r"\1 (\2)", text)
    text = re.sub(r"`([^`]+)`", lambda m: f"<font name=\"Courier\">{escape(m.group(1))}</font>", text)
    text = re.sub(r"\*\*(.+?)\*\*", lambda m: f"<b>{m.group(1)}</b>", text)
    text = re.sub(r"\*(.+?)\*", lambda m: f"<i>{m.group(1)}</i>", text)
    text = re.sub(r"<b>.*?</b>|<i>.*?</i>", preserve_tag, text)
    escaped = escape(text)
    for index, tag in enumerate(tokens):
        escaped = escaped.replace(f"REPORTMARKUPTOKEN{index}END", tag)
    return escaped


def _page_footer(canvas, document) -> None:
    canvas.saveState()
    canvas.setFont("Helvetica", 8)
    canvas.setFillColor(colors.HexColor("#687887"))
    canvas.drawString(0.62 * inch, 0.35 * inch, "HG03086 self-insertion report · generated from selfins_report.md")
    canvas.drawRightString(letter[0] - 0.62 * inch, 0.35 * inch, f"Page {document.page}")
    canvas.restoreState()


pdf_structure = render_pdf(md, PDF_PATH)
print("wrote", REPORT_PATH, len(md), "chars")
print("wrote", PDF_PATH)
print(
    "PDF parity:",
    f"{len(pdf_structure['sections'])} sections in Markdown order",
    f"{pdf_structure['figures']} figures",
    f"{len(pdf_structure['table_rows'])} tables with data-row counts {pdf_structure['table_rows']}",
    f"{pdf_structure['pages']} pages",
)
