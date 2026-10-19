
import pandas as pd, json, numpy as np
parts = []
W = parts.append

W(f"""# Mobile element insertions in HG03086: orientation, element depletion, and gene constraint

**Sample:** HG03086, single genome, Illumina short-read only
**Input:** `HG03086_classifier_ge_0.997.genes.snpeff.vcf` — {int(SCOPE['total classifier calls']):,} classifier calls at score ≥ 0.997
**Assembly:** GRCh38 (`##reference=hg38`)
**Annotation:** snpEff (`rtm annotate-genes`), as merged in PR #51

---

## Summary

Three analyses requested by W. Brandler, all reproducing known population-genetic
expectations from a single genome.

1. **Sense-oriented intronic insertions are depleted.** Sense fraction {comb.sense_frac:.3f}
   ({int(comb.sense)} sense / {int(comb.antisense)} antisense, n = {int(comb.n)}),
   Holm-adjusted p = {comb.p_holm:.4f}. Direction matches Medstrand et al. 2002.
2. **Coding exons are ~8× depleted, and it is not base composition.**
   O/E = {gc2.loc['coding_exon','obs_over_exp']:.3f} against a GC-matched null
   (Holm p = {gc2.loc['coding_exon','p_holm']:.3f}).
3. **Promoter depletion is entirely base composition.** O/E = {un2.loc['promoter','obs_over_exp']:.2f}
   under a uniform null but exactly {gc2.loc['promoter','obs_over_exp']:.2f} under a GC-matched null.
   **This is the most important methodological result in the report** — the uniform null
   alone would have produced a false finding.
4. **Intronic depletion is concentrated in the most constrained genes.** The lowest LOEUF
   decile carries {int(b0.obs)} insertions against {b0.exp:.0f} expected
   (O/E = {b0.oe:.2f}, 95% CI {b0.lo:.2f}–{b0.hi:.2f}); the other nine deciles sit at 1.0.

Two annotation bugs were found in the input and are reported in §6.

---

## 1. Data

{scope_table}

Of the {int(SCOPE['intronic calls (Task 1 denominator)']):,} intronic calls,
{int(SCOPE['  of which multi-gene'])} overlap more than one gene. `CSQ_TERMS` is a *distinct set*
of consequence terms rather than a per-gene aligned list, so for those records it is not
recoverable which gene received the intronic hit. **The primary analysis therefore uses the
{int(SCOPE['  of which single-gene'])} single-gene intronic calls**, where `GENEID` and
`GENE_STRAND` align 1:1 and the assignment is unambiguous.

---

## 2. Task 1 — Intronic insertions are depleted in the sense orientation

### Method
Orientation relative to the host gene was derived as `sense` when the insertion's `ORIENT`
equals the gene's `GENE_STRAND`, `antisense` otherwise — the annotator does not emit this
field (§6). Each stratum was tested against a 50:50 null with a two-sided exact binomial
test; the four tests (three families plus combined) were Holm-corrected.

### Results

{t1_table}

{fig1}

**Internal control.** Insertions annotated `upstream_gene_variant` or
`downstream_gene_variant` — beside genes rather than inside them, where transcriptional
interference cannot operate — give a sense fraction of {fl.sense_frac:.3f}
({int(fl.sense)}/{int(fl.antisense)}, n = {int(fl.n)}), sitting on the null as expected.

**This control is weaker than it looks.** Its 95% CI is {fl.ci_low:.2f}–{fl.ci_high:.2f}, wide
enough to contain both 0.50 and the intronic estimate of {comb.sense_frac:.3f}. The direct
intronic-vs-flanking contrast is not significant (Fisher exact, two-sided, p = 0.41). The
control is *consistent with* selection being confined to introns; it does not demonstrate it.

**Power.** Computed for a one-sample binomial against p₀ = 0.50, α = 0.05 two-sided:

{power_table}

Alu and the combined test are well powered. **LINE-1 and SVA are not.** The non-significant
SVA result is consistent with anything down to a sense fraction of {sva_detect:.2f} and must
be reported as underpowered, not as absence of an effect.

**Sensitivity.** Adding the {n_ext} multi-gene intronic calls whose genes all share one strand
(orientation unambiguous even though the specific gene is not) gives n = {n_ext_tot} and a
combined sense fraction of {ext_comb:.3f} — unchanged in direction or magnitude.

### Interpretation
The direction reproduces Medstrand et al. 2002, who found L1 and LTR elements
*"particularly those in the same orientation as the nearest gene, are significantly
underrepresented within genes"*. The mechanism is described by Han, Szak & Boeke 2004:
sense-oriented L1 inside an intron impairs transcriptional elongation, an effect they
attribute to the A-rich bias of the L1 sense strand and show to be length-dependent.

**A caveat that limits how far this interpretation can be carried.** Alu accounts for
{int(prim[prim.family=='ALU'].n.iloc[0])} of the {int(comb.n)} intronic calls — 86% of the
signal. Han's mechanism is specific to L1 ORF sequence and explicitly length-dependent
(*"longer ORF2 sequences led to lower concentrations of RNA"*; tandemised 5'UTR, which lacks
the A-rich bias, *"had no effect"*). Alu is ~300 bp, has no ORFs, and lacks that bias. The
mechanism therefore does not obviously extend to the elements that dominate this dataset.
The orientation bias in Alu is real in these data but unexplained by the cited mechanism.

---

## 3. Task 2 — Element depletion by simulation

### Method
Observed and expected counts were produced by the **same** interval model, so the comparison
is not confounded by differences between two annotation systems. Element intervals
(CDS, 5'UTR, 3'UTR, promoter, gene body) were built from Ensembl GRCh38 release 114,
protein-coding transcripts only, merged per class, and assigned by priority
coding exon > 5'UTR > 3'UTR > promoter > intron > intergenic. Promoter was defined as
TSS −2,000/+200 bp; this was a choice, not a specification from the brief.

The accessible genome (non-N sequence on chr1–22, X) is {acc_gb:.3f} Gb. All
{int(SCOPE['total classifier calls']):,} observed calls fall inside it.

Two nulls, 1,000 replicates each, {int(SCOPE['total classifier calls']):,} positions per replicate:

- **uniform** — positions drawn uniformly over accessible sequence.
- **GC-matched** — positions drawn so the distribution of local GC (1 kb window) matches the
  observed calls. Required because these elements land in AT-rich DNA: observed calls have
  mean local GC {obs_gc:.4f} against a genome-wide {gen_gc:.4f}, and functional elements
  differ systematically in GC. Without this control, GC-rich element classes appear depleted
  for reasons that have nothing to do with selection.

**Deviation from the original plan, for review.** The plan was to push each simulated
replicate through `rtm annotate-genes`. At William's measured 59 s per 1,852 calls this is
~16 CPU-hours for 1,000 replicates. Using one interval model for both sides instead took
19 s for 2,000 replicates and removes the observed-vs-expected annotation mismatch entirely.
The trade-off is that the model is mine rather than snpEff's; §3.1 validates it.

### 3.1 Validation of the interval model

| check | result |
|---|---|
| CDS as a fraction of accessible genome | {cds_pct:.2f}% (published range 1.2–1.5%) |
| analytic vs simulated uniform expectation, coding exon | 22.4 vs {un2.loc['coding_exon','expected']:.1f} |
| coding-exon calls found by this model | 2 — the same two snpEff flagged as inframe insertions |
| snpEff `intron_variant` calls inside any-biotype gene body | 1,203 / 1,204 |

The last row explains a discrepancy worth stating plainly: this model counts
{int(OBS_CLS['intron'])} intronic calls where snpEff reports
{int(SCOPE['intronic calls (Task 1 denominator)']):,}. snpEff annotates against all
transcript biotypes; this model is restricted to protein-coding genes. The 451-call gap is
lncRNA and pseudogene introns, not an error — but it means "intronic" denotes different sets
in §2 and §3, and the two sections' counts should not be compared directly.

### Results

{t2_table}

{fig2}

**Coding exons.** The only class that survives correction under both nulls. Two insertions
observed where {gc2.loc['coding_exon','expected']:.1f} are expected under composition
matching — an {1/gc2.loc['coding_exon','obs_over_exp']:.0f}-fold depletion,
Holm p = {gc2.loc['coding_exon','p_holm']:.3f}.

**Promoters.** {int(un2.loc['promoter','observed'])} observed against
{un2.loc['promoter','expected']:.1f} expected under the uniform null (O/E
{un2.loc['promoter','obs_over_exp']:.2f}) — but {gc2.loc['promoter','expected']:.1f} expected
under GC matching, giving O/E {gc2.loc['promoter','obs_over_exp']:.2f}. Promoters are
GC-rich and these elements are AT-biased; that accounts for the entire apparent depletion.

**UTRs.** Both are depleted in point estimate and neither survives correction. 5'UTR has
n = {int(gc2.loc['5UTR','observed'])} and no power at this sample size.

---

## 4. Task 3 — Coding-exon disruptions and gene constraint

### 4.1 The requested table

Across the whole genome, **two** calls disrupt a coding exon.

{cod_table}

Constraint scores are gnomAD v4.1, MANE Select transcript — the release served by the
gnomAD website. v2.1.1 values are retained in the CSV as a sensitivity column; the two
releases disagree on individual genes (v4.1 has roughly five times the exomes) but classify
both of these genes as LoF-tolerant.

**n = 2 is the result, not a shortfall.** Coding-exon insertions are under strong purifying
selection (§3), so a single genome is expected to carry almost none. Two observed against
{gc2.loc['coding_exon','expected']:.1f} expected is the finding.

### 4.2 Extension: constraint across all intronic insertions

Because n = 2 supports no statistical claim, the same biological question — are genes that
cannot tolerate damage protected from insertions? — was asked of all intronic calls, where
there is power. Poisson regression of per-gene insertion count on constraint, with an offset
for exposure. Two offsets were used: log intronic length, and log expected count from the
GC-matched null (which absorbs length *and* composition).

**Length correction is mandatory here.** Constrained genes have a median intronic length of
61.2 kb against 28.0 kb for tolerant genes. Uncorrected, that difference alone produces
apparent enrichment.

{t3_table}

{fig3}

All four parameterisations point the same way; two clear p = 0.05. The binned view is
sharper than the regression: the depletion is **not a gradient**. The most constrained LOEUF
decile sits at O/E {b0.oe:.2f} with its interval excluding 1.0, and the remaining nine
deciles scatter around 1.0. That is a threshold effect, and it is not what a smooth
dose-response model would predict.

Both coding-exon-disrupting genes fall far into the tolerant tail
(*OVCH2* LOEUF {ovch2:.2f}, *ZNF404* LOEUF {znf:.2f}) — consistent with, but in no way
demonstrated by, n = 2.

---

## 5. What these three results say together

Three analytically independent cuts of one genome point the same way: purifying selection
against insertions, strongest where the consequence is greatest. Sense-oriented intronic
insertions are depleted relative to antisense; coding exons are depleted ~8× relative to
composition-matched expectation; and intronic insertions are depleted specifically in the
genes least tolerant of loss of function.

None of this is novel — all three are established population-genetic findings, and
reproducing them from a single Illumina genome was the point of the exercise. The value is
that the pipeline recovers them.

---

## 6. Bugs found in the annotated VCF

**6.1 `MEI_GENE_ORIENT` is absent.** The output carries `GENE`, `GENEID`, `GENE_STRAND`,
`CSQ`, `CSQ_TERMS`, `CSQ_NTX` but not the derived insertion-vs-gene orientation, which had
existed on the VEP backend. Every analysis in §2 had to re-derive it. Recommend re-adding it
to the snpEff path as one comma-separated entry per gene, index-aligned with `GENEID`.

**6.2 Gene-keyed INFO lists fall out of register.** In four records `GENE` has fewer entries
than `GENEID` and `GENE_STRAND`, so code that zips them positionally silently mis-assigns:

| locus | GENE | GENEID | GENE_STRAND |
|---|---|---|---|
| chr5:77395411 | 1 entry | 2 | 2 |
| chr7:92432101 | 1 entry | 2 | 2 |
| chr16:2362374 | 2 entries | 3 | 3 |
| chr19:53287456 | 2 entries | 3 | 3 |

Recommend one entry per gene in `GENEID` order, substituting the Ensembl ID when a symbol is
missing, with a test asserting equal lengths across all records.

---

## 7. Limitations

1. **One genome.** No population-level inference is possible; every result is a
   reproduction of a known effect, not new evidence for it.
2. **Classifier threshold 0.997.** Only high-confidence calls are visible. All absolute
   counts are lower bounds, and any threshold-correlated bias propagates into the O/E ratios.
3. **The 50:50 orientation null is an assumption.** It is the field-standard null but
   assumes integration is orientation-agnostic; this was not independently verified here.
4. **The promoter window is a choice.** TSS −2,000/+200 bp. snpEff's `upstream_gene_variant`
   uses a 5 kb window and is not equivalent; results for this class are definition-dependent.
5. **GC matching may not be the right composition control.** It matches 1 kb local GC.
   L1 endonuclease has a degenerate TTTT/AA-like target preference that GC content only
   approximates; a motif-anchored null would be stricter and was not run.
6. **Multi-gene calls are excluded from Task 1 primary.** 348 of 1,204 intronic calls.
   The sensitivity analysis recovers {n_ext} of them without changing the result.
7. **Two definitions of "intronic" appear in this report** (§3.1) — protein-coding only in
   §3 and §4, all biotypes in §2.
8. **Assumption diagnostics were not run** on the Poisson models (overdispersion, influence).
   Not assessed.

---

## 8. References

1. Medstrand P, van de Lagemaat LN, Mager DL (2002). Retroelement distributions in the human
   genome: variations associated with age and proximity to genes. *Genome Res* 12:1483–1495.
   doi:10.1101/gr.388902
2. Han JS, Szak ST, Boeke JD (2004). Transcriptional disruption by the L1 retrotransposon and
   implications for mammalian transcriptomes. *Nature* 429:268–274. doi:10.1038/nature02536
3. Belancio VP, Hedges DJ, Deininger P (2006). LINE-1 RNA splicing and influences on mammalian
   gene expression. *Nucleic Acids Res* 34:1512–1521. doi:10.1093/nar/gkl027
4. Cingolani P, et al. (2012). A program for annotating and predicting the effects of single
   nucleotide polymorphisms, SnpEff. *Fly* 6:80–92. doi:10.4161/fly.19695
5. Lek M, et al. (2016). Analysis of protein-coding genetic variation in 60,706 humans.
   *Nature* 536:285–291. doi:10.1038/nature19057
6. Karczewski KJ, et al. (2020). The mutational constraint spectrum quantified from variation
   in 141,456 humans. *Nature* 581:434–443. doi:10.1038/s41586-020-2308-7
7. Stewart C, et al. (2011). A comprehensive map of mobile element insertion polymorphisms in
   humans. *PLoS Genet* 7:e1002236. doi:10.1371/journal.pgen.1002236

**Data sources.** Ensembl GRCh38 release 114 GTF; UCSC hg38.2bit; gnomAD v4.1 and v2.1.1
constraint metrics (AWS Open Data mirror).
""")
REPORT = "".join(parts)
