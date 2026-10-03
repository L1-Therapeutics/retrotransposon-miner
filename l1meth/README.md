# Locus-specific L1 promoter methylation for tumour/normal discrimination

A feasibility and design study. **No bisulfite data was analysed.** What is
real here is the locus catalogue (reference coordinates and CpG positions) and
the reference-computable selection features; what is simulated is the
methylation itself, calibrated to published distributions. Treat the AUC
numbers as a statement about *study design*, not as a validated classifier.

## Why this question

Bulk bisulfite assays report one LINE-1 methylation number averaged over
~500,000 copies. The motivation for locus resolution is that an average cannot
localise anything. But a published multi-cancer study of L1 promoter
methylation at 11 loci found **no significant clinical advantage of
locus-specific markers over global markers** for separating tumour from normal
tissue. Any proposal for locus resolution has to engage that result rather
than ignore it.

## What was built

| module | role |
|---|---|
| `catalog.py` | addressable L1 loci + promoter CpG coordinates, from RepeatMasker + reference sequence |
| `simulate.py` | per-locus methylation for tumour/normal, calibrated to published distributions |
| `compare.py` | bulk average vs per-locus vector, same estimator and folds |
| `panel.py`   | panel selection: blind vs reference-informed |

Run `python run_study.py`; tests in `tests/test_l1meth.py` (26 passing).

## Definition of an addressable locus

A locus is addressable when all of the following hold. Each condition removed
loci that a naive filter would have kept, and two of them were found only by
inspecting real data.

1. **Retains 5'UTR sequence.** The promoter CpG island *is* the silencing
   switch; a 3'-truncated copy has nothing to measure.
2. **Subfamily consensus >= 5 kb.** *Found the hard way, twice.* Some ancient
   subfamilies are fragment-derived: 16 of the 125 L1 subfamilies on chr22 have
   consensus models under 5 kb, and two are shorter than the 900 bp promoter
   window itself -- **L1ME4b (892 bp) and L1ME4c (830 bp)**, with L1ME3G
   (924 bp) barely longer. For those the "first 900 bp" is effectively the
   whole model, a naive `cons_start < 900` test passes for their copies, and
   they dominated the first catalogue built -- 278 loci, mostly degraded
   fragments with no promoter at all.

   Consensus length itself was initially computed wrong. UCSC orders the three
   consensus-coordinate fields **differently by strand**: for `-` rows
   `repStart` holds the negative remainder and `repLeft` the consensus start,
   the reverse of `+` rows. Applying the `+` formula to every row adds a
   coordinate to a length, which inflated L1ME4b's 892 bp model to values up
   to 11,799 bp and yielded 179 distinct "consensus lengths" for one
   subfamily. Consensus length is a property of the **subfamily**, not of a
   copy, so `subfamily_consensus_length()` reduces rows by mode (the median
   lands on a real value only by luck). Correcting this *added* 9 addressable
   loci, 34 -> 43: the corrupted column had been spuriously excluding primate
   copies.
3. **Primate (L1P) lineage.** The modern human-type 5'UTR promoter arose
   through repeated 5'UTR replacement; older L1M 5' ends are a different
   architecture and are not comparable CpG-for-CpG.
4. **>=150 bp of contiguous non-repeat, non-N flank.** A read must reach
   unique sequence to be assigned to this copy. Assembly gaps are runs of N and
   are not mappable -- masking them removed one locus whose apparent 50 kb
   flank was centromeric gap.
5. **>=3 measurable promoter CpGs.**

**Result on chr22: 43 addressable loci carrying 512 promoter CpG sites**, out
of 11,526 L1 copies. Scaling naively by chromosome length gives ~2,600
genome-wide, the same order as the "thousands of L1 copies" profiled by
bs-ATLAS-seq -- an independent check that the definition is not absurd.

## Findings

**1. Locus resolution is not free.** With only a sample-level (global)
hypomethylation shift and no locus-specific component, the per-locus
representation scores *worse* than the bulk average (AUC 0.838 vs 0.897): the
extra parameters buy nothing and cost variance. The published null is
reproducible from this alone.

**2. It wins as soon as any locus-specific biology exists.** At >=5% of loci
carrying a tumour-specific change, per-locus overtakes the bulk average
(0.959 vs 0.932); the advantage peaks near 10% (0.984 vs 0.943) and closes
again once nearly all loci move together, where the bulk mean captures
everything.

**3. Panel size does not explain the published null -- panel *selection*
does.** A blind 11-locus panel drawn from a 2,100-locus pool with 10% signal
has only a **70.4% chance of containing even one signal locus** (1,500
replicates). Roughly three studies in ten would report a null regardless of
whether the biology is real.

**4. Reference-computable selection could fix that -- but not for the outcome
this study targets. FALSIFIED against published data; retained here because
the failure is the result.**

The original claim was: order candidate loci by divergence-from-consensus,
take the youngest, and the hit rate rises from 70.7% to 94.8% at an 11-locus
panel (+24.1 pp; 1.11 -> 2.66 signal loci). That holds *only* under an assumed
8x enrichment of signal in young loci, which was motivated but never measured.

Checking it against Lanciano et al. 2024 -- ~12,000 Illumina 450K samples from
GEO, young L1PA probes n = 695, old n = 189 -- the assumption does not survive
for carcinomas. Young L1PA methylation is "high and generally similar to or
higher than that of old L1PAs in most situations." Young-L1 hypomethylation is
specific to pluripotent stem cells, trophoblast, embryonal carcinoma,
seminoma, placenta, fetal membranes and hydatidiform moles; in normal
fibroblasts and several carcinoma lines (MCF-7, HeLa-S3, HepG2) young L1PAs
are relatively *hyper*methylated.

With `age_enrichment=1.0` -- now the default -- age-informed selection gains
**nothing**: 68.2% vs 70.4% blind at an 11-locus panel (-2.2 pp, against a
sampling SE of ~1.7 pp at 1,500 replicates, so not distinguishable from zero;
1.10 vs 1.11 signal loci). The advantage was the assumption, not the data.

What survives: findings 1-3 are untouched, since none of them uses
`age_enrichment`. And the falsification is specific rather than total -- the
young-L1 hypomethylation the strategy needs *is* observed in germ-cell and
embryonal tumours. A locus-resolved study of **seminoma or embryonal
carcinoma** would be betting on an effect that has been measured, which
tumour-versus-normal in carcinoma is not.

One caveat in the other direction: the array study compares young-vs-old
methylation *levels within a sample*, while panel selection needs young-vs-old
*differential* methylation between tumour and normal. Strong evidence against
the assumption, not a direct refutation.

## Calibration and its weakest link

| parameter | value | source |
|---|---|---|
| tumour L1 methylation mean / SD | 61.4 / 9.6 | Baba et al. 2010, *Mol Cancer* 9:125 (n=869 CRC, pyrosequencing) |
| normal colonic mucosa mean | ~76 | **weakest link** -- typical pyrosequencing figure, not verified to a specific source in this study; exposed as a parameter |
| informative reads per junction | 10 | arithmetic: 30x WGBS, 150 bp reads, >=50 bp each side of the junction => coverage x 0.33 |
| age enrichment of signal | **1.0 (none)** | was an 8x assumption; **falsified for carcinomas** by Lanciano et al. 2024 (~12,000 450K samples). Raise only for germ-cell/embryonal contexts |

## What to do next, in order

1. **Stop proposing age-informed panels for carcinoma tumour/normal.** The
   evidence is against it and a reviewer who knows this literature will say so.
2. **If the tumour/normal framing is kept**, the honest pitch is findings 1-3:
   locus resolution helps when locus-specific biology exists, the published
   null is explained by panel selection rather than panel size, and a blind
   11-locus panel misses signal ~30% of the time. None of that needs the
   falsified assumption.
3. **If young-L1 biology is the point**, switch tissue: germ-cell and
   embryonal tumours are where the effect is measured.
4. **Either way, replace the simulated matrix with real data.** The same 450K
   GEO datasets Lanciano used are public, and `compare.py` takes a matrix and
   labels regardless of where they came from.

## Limitations

- **Simulation, not measurement.** No real WGBS was analysed. The natural next
  step is to run the catalogue against public bisulfite data (ENCODE, TCGA,
  GEO) and replace the simulated matrix with a measured one.
- Finding 4 is falsified for carcinomas and the default now reflects that.
  The module retains the `age_enrichment` parameter so the optimistic case can
  be reproduced, and two tests pin both behaviours.
- Genic context is annotated in the catalogue but **not** used for selection:
  Lanciano et al. link intronic L1 methylation to host transcription *within*
  a sample, which is not evidence that intronic loci carry more tumour-normal
  differential signal. Replacing one unmeasured selection assumption with
  another would repeat the mistake finding 4 just made.
- "Non-repeat, non-N flank" is a proxy for mappability. True uniqueness needs
  k-mer analysis against the whole genome; these counts are an upper bound.
- chr22 only, and chr22 is GC-rich and L1-poor, so the genome-wide scaling
  likely *under*counts.
- Sweeps use n_loci=43, the corrected catalogue size. Re-running at the
  earlier (wrong) size of 35 shifts the AUCs by <0.02 and changes no
  conclusion.
