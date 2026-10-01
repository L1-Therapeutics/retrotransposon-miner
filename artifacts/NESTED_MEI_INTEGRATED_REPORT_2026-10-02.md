# Nested mobile-element insertions in HG03086 — integrated evidence report

**Prepared 2 October 2026.** This report integrates the saved HG03086 analysis, a read-only audit of the full 1,852-call VCF against the saved per-call table and hg38 RepeatMasker track, and checks of primary papers. The report does not combine the separate 6,939-candidate `nested_analysis/results_highconf` cohort with the 1,852 classifier-call cohort. The Claude export's manifest lists an 835,393,456-byte `hg38.2bit` and 50,818,620-byte `chr22_seq.json`, but the exported JSON contains metadata/transcript references, not these file payloads; neither file is present in this checkout. Exact restoration/data requirements are in [the validation input guide](SELFINS_VALIDATION_INPUT_GUIDE.md).

**Scope note:** the updated [legacy report](selfins_report.md) and its generator retain several older exploratory sections/claims not reproduced here; use this integrated report for the current conservative synthesis. Its PDF was not regenerated.

## Executive summary

- In one genome-wide HG03086 callset (1,555 Alu, 222 LINE-1, 75 SVA), 239 Alu, 72 LINE-1 and 9 SVA insertions were annotated as nested in their own-family reference elements. The reported enrichment against the selected resampling nulls is 1.47× for Alu (239 vs 162.5 expected), 1.47× for LINE-1 (72 vs 48.9 expected), and 29.4× for SVA (9 vs 0.306 expected). The SVA numerator is only nine; this is a candidate signal, not an established HG03086 effect, and the reference callability of SVA was not tested. Nummi et al. 2025 reports SVA self-target preference, so this is not a first systematic precedent.
- Conditional on nesting, same orientation was observed for 202/239 Alu, 53/72 LINE-1 and 9/9 SVA. The separate sense/antisense resampling table reports 202 Alu sense vs 81.27 expected (2.49×) and 37 antisense vs 81.27 (0.46×); LINE-1 53 vs 24.43 (2.17×) and 19 vs 24.43 (0.78×); SVA 9 vs 0.153 and 0 vs 0.153. Those are comparisons of each class against an orientation-blind placement null, not a direct test of sense vs antisense within nested events. The within-nested sense fractions and their binomial tests are provided separately in [selfins_orientation.csv](selfins_orientation.csv).
- The saved Alu host-position profile has peaks at 120–140 bp (29 calls, 3.0× uniform bin expectation) and 280–300 bp (36, 3.7×), with an additional 0–20 bp peak (13). The 120–140 peak is positionally compatible with Levy et al.'s Alu–Alu insertion hotspot near consensus coordinate 133; 133 lies in this 20-bp bin under the stated consensus-to-host 5′ alignment. This is a coarse coordinate correspondence, not nucleotide-resolution replication.
- The 280–300 cluster is at the host Alu's terminal A-rich tail under this report's host-coordinate convention. It is outside Levy et al.'s target-end positional analysis, which excluded target consensus ends because of repeat-mapping limitations. It is therefore an open, exploratory result—not an established novel hotspot and not yet a tested TPRT-hallmark-positive cluster.
- A new exact-coordinate join shows that the VCF provides TSD strings for 1,626/1,852 records and a `POLYA_MIN_BP` classifier field for every record. A basic TSD-string histogram is bimodal-looking, but it is **not** an independent TSD validation: the VCF TSD sequence is a caller-supplied candidate, not a reconstruction from both insertion junctions. The reported `POLYA_MIN_BP` is likewise classifier evidence, not a directly measured inserted 3′ tail. Current evidence does not support a triad-frequency estimate for this callset.
- Nummi et al. (2025) is the strongest newly checked comparative context: from nanopore data they report germline TEs enriched in repeats; germline SVA insertions preferentially target reference SVA, and their Alu insertions peak in two sequence intervals overlapping this callset's two positional clusters. That is independent support for repeat-target preference and Alu positional structure, but their cohort, ascertainment, reference, and positional method differ; it does not validate our tail-peak calls one-by-one or remove our callability concerns.

## 1. Cohort and existing analyses

The classifier callset is HG03086, GRCh38, chromosomes 1–22 and X, classifier score ≥0.997. The count is 1,852 (1,555 Alu, 222 LINE-1, 75 SVA); calls are polymorphic catalogue events according to the handoff. The input VCF is `HG03086_classifier_ge_0.997.genes.snpeff (1).vcf`; the saved per-call table is [selfins_per_call.csv](selfins_per_call.csv).

Other completed analyses requested by William are documented in [HG03086_MEI_report.md](HG03086_MEI_report.md): intronic orientation relative to genes (382/856 sense; combined Holm-corrected exact-binomial p=0.0074), depletion of coding-exon insertions (2 observed vs 16.6 GC-matched expected; Holm p=0.01199), and promoter depletion disappearing under GC matching (33 vs 42.9 uniform expected but 33 vs 33 GC-matched expected). Both coding-exon calls are in low-pLI/high-LOEUF genes in the reported gnomAD version; see [task3_coding_exon_table.csv](task3_coding_exon_table.csv). These analyses concern functional genomic annotation and must not be confused with target-repeat nesting.

### Nested enrichment results and correct null language

| New insertion family | Calls | Same-family nested | Selected null | Expected nested calls | Enrichment | Empirical p |
|---|---:|---:|---|---:|---:|---:|
| Alu | 1,555 | 239 | uniform | 162.548 | 1.47× | 0.000999 |
| LINE-1 | 222 | 72 | GC-matched | 48.866 | 1.47× | 0.000999 |
| SVA | 75 | 9 | GC-matched | 0.306 | 29.4× | 0.000999 |

Each null used 1,000 replicates. The selected null differs by family: GC-matching moves Alu expected nesting down to 149.636 (which would make the estimate 1.60×), while for LINE-1 and SVA it moves expected nesting up relative to uniform. The family-specific model selection is recorded in [selfins_enrichment.csv](selfins_enrichment.csv); the manuscript should not describe the three ratios as if all had the same null. With 1,000 replicates, p=0.000999 is the empirical floor, not a precise tail probability.

SVA must be reported as **9 observed vs 0.306 expected**, not as “29×” alone. The current GC-matched placement null has no repeat-aware SVA-callability or mappability component. We therefore do not know whether its expectation is biased by which reference SVAs are detectable by the short-read pipeline. No SVA-specific inference should be called robust until that is tested.

### Nested orientation: denominators kept separate

| Family | Nested total | Same orientation (within nested) | Fraction (95% CI) | Exact test vs 0.5 | Fisher against orientation marginals |
|---|---:|---:|---:|---:|---:|
| Alu | 239 | 202 | 0.845 (0.793–0.889) | 1.06e-28 | 1.37e-28 |
| LINE-1 | 72 | 53 | 0.736 (0.619–0.833) | 7.56e-05 | 0.000202 |
| SVA | 9 | 9 | 1.000 (0.664–1.000) | 0.00391 | 0.00794 |

Alu's insertion and reference plus-strand marginals are both close to 50%, so its independent-orientation expectation is 0.500. For LINE-1 the expected same-orientation rate from the callset marginals is 0.522, and for SVA it is 0.506. Consequently the null is approximately one-half but is not mathematically identical for each family. The tests above are conditional on the 320 nested calls; the enrichment tests below compare sense/antisense categories separately to random-placement expectations and are not a direct within-nested contrast.

| Family | Category | Observed | Expected under selected orientation-blind placement null | O/E | p (saved table) |
|---|---|---:|---:|---:|---:|
| Alu | nested sense | 202 | 81.274 | 2.49× | 3.39e-29 |
| Alu | nested antisense | 37 | 81.274 | 0.46× | 6.16e-08 |
| LINE-1 | nested sense | 53 | 24.433 | 2.17× | 7.56e-07 |
| LINE-1 | nested antisense | 19 | 24.433 | 0.78× | 0.318 |
| SVA | nested sense | 9 | 0.153 | 58.8× | 2.21e-13 |
| SVA | nested antisense | 0 | 0.153 | 0× | 1.0 |

These category-vs-placement comparisons can establish deviation from the specific simulation null, but their p-values do not substitute for a direct paired/conditional comparison of sense against antisense among nested calls. That direct test is available from the orientation counts; do not combine the two nulls into one claim. The saved rates are [selfins_sense_antisense.csv](selfins_sense_antisense.csv) and [selfins_orientation.csv](selfins_orientation.csv).

## 2. Position along Alu hosts and literature comparison

The saved profile is based on 154 near-full-length host Alus (280–320 bp) and binned in 20-bp windows. The expected uniform count is 9.62/bin. It records 29 calls in 120–140 bp (host A-fraction 0.554), 36 in 280–300 bp (A-fraction 0.792), and 13 in 0–20 bp (A-fraction 0.125). Across all 239 nested Alus, positional distribution is nonuniform (χ²=103.13, df=9, p=3.65e-18; KS D=0.161, p=7.02e-6). Source: [selfins_position_profile.csv](selfins_position_profile.csv) and [selfins_position_stats.json](selfins_position_stats.json).

Levy, Schwartz & Ast (2009 online; 2010 issue, *Nucleic Acids Research* 38:1515–1530, doi:10.1093/nar/gkp1134; Results p. 1518, Figure 2) report a prominent Alu-into-Alu insertion hotspot directly after the A-rich linker at targeted-Alu consensus coordinate 133. Their Figure 2 describes the AluJo consensus and places the internal A-rich linker at 118–136; Figure 1 defines insertion position as the coordinate just upstream of the insertion. It is appropriate to say that under a 5′-aligned consensus-to-host convention, coordinate 133 falls in our 120–140 bp bin. The insertion breakpoint is around the 133/134 boundary, not an exact “133-base bin” edge. The paper does not explicitly state that 133 is counted from the first consensus nucleotide, so this 5′ coordinate mapping is an alignment inference. Our current bin is coarse, with no consensus-aligned breakpoint sequence; the finding is a coordinate-compatible reproduction, not base-level remeasurement. Figure 2's sense-oriented Alu self-insertion display uses 15,759 Alu self-insertion events.

### Levy's “~4-fold” — exact unit of analysis

The “~4-fold” figure from Levy et al. 2010 is not the ratio of total insertions into TE sequence to insertions into non-TE sequence. It is the ratio of enriched to depleted insertion-type categories: 516 overrepresented versus 116 underrepresented out of 36,478 tested (young × old × orientation), i.e. ~4.45:1 by category. Their per-type null is random insertion into unique intergenic DNA scaled to each old TE class’s inferred ancestral coverage. Any total-event ratio must be computed from the underlying event counts, which the paper does not report as a single number.

Levy's positional analysis excluded consensus ends due to RepeatMasker/short-flank identification limitations. This means our 280–300bp tail cluster falls outside that particular positional analysis; it is a coverage blind spot, **not evidence by itself that the cluster is novel or genuine**. The older paper's 4-fold category ratio and our family-specific placement enrichment are not directly comparable.

### Nummi et al. 2025, independent long-read study

Primary source: Nummi et al., *Mobile DNA* 16:20 (2025), doi:10.1186/s13100-025-00357-w. They analyzed long-read whole-genome data from 56 colorectal cancers and 112 uterine leiomyomas (plus normal samples), reporting 1,495 somatic L1 insertions and 7,265 distinct germline TE polymorphisms in the study. Their Table 1's target-hallmark comparison includes only 897 germline L1s, a subset of that total, with 569/897 (63.4%) meeting their EN-site criterion; for somatic L1 it is 1,048/1,495 (70.1%). They report TSDs for 65.4% of somatic L1 and 69.5% of germline L1, and median L1 insertion lengths of 394 vs 1,240 bp. Those metrics are L1-specific and not comparable to our mixed Alu/L1/SVA classifier fields.

Their Figure 5 and Results report germline insertions enriched in repetitive DNA; germline insertions preferentially target their own repeat type, with the strongest preference reported for SVA insertions into reference SVAs. This makes our 9 SVA nested calls biologically plausible but does not establish a first observation; Nummi is now a systematic precedent. Their Alu target-position plot shows peaks across positions 100–155 and 260–310, with AT peaks at 115–145 and 280–310. That independently supports two positional territories corresponding to our 120–140 and 280–300 bins. Do not compare event counts or effect sizes as if cohorts and methods were matched. Their 26% short-read miss figure refers to long-read somatic insertions not detected by short-read data, not to a repeat-overlap rate.

Nummi explicitly cautions that germline polymorphisms have passed through evolutionary selection. These are survivor distributions, not raw de novo targeting rates. This applies to our own 1000 Genomes catalogue calls too. Their study makes the broad self-targeting and two-region Alu positional observations less novel; our contribution is the single-genome classifier analysis, orientation decomposition and the independent historic Levy-position compatibility check.

## 3. Tail-cluster and TSD/poly(A) audit — what is computed and what is not

Reproducible read-only audit: [audit_nested_tsd.py](audit_nested_tsd.py). Outputs: [selfins_tsd_audit_summary.json](selfins_tsd_audit_summary.json), [selfins_tsd_call_qc.csv](selfins_tsd_call_qc.csv). The script exact-joins VCF and saved per-call rows on `(chrom, POS)` and rejects any mismatch. It streams hg38 RepeatMasker to reconstruct saved nested host intervals and 5′-relative offsets. No source records are altered.

**Join totals:** 1,852/1,852 coordinates match; 1,626 VCF TSD strings are non-empty and 226 are missing; `POLYA_MIN_BP` is present on all calls. The VCF string length histogram has a prominent 14–17 bp central mode (906 calls) and a short 9–11 bp mode (207 calls), with long reported values through 40 bp. The lengths have caller boundary ambiguity and possible low-complexity/poly(A/T) matches; therefore this is a histogram of *reported TSD strings*, not of validated target-site duplications.

The Szak et al. 2002 TSDfinder figure (16,266 genomic L1 TSD predictions) reports a bimodal length distribution and notes that 9–11 nt predictions may contain false positives; their method searched young, 3′-intact L1s and assigned TSDs from reference sequence. Zingler et al. 2005 reproduced peaks at 9 and 15 nt; after filtering noncanonical L1 structures, the short 9-nt peak dropped sharply and the main peak shifted to 16 nt. This is an appropriate QC caution for short reported strings, not a universal exclusion rule for Alu/SVA or our caller.

The audit's Alu near-full-host profile reconstructed from GRCh38 RepeatMasker using normalized VCF POS gives 152 calls rather than the saved profile's 154; 30 at 120–140 and 33 at 280–300 (rather than saved 29 and 36). This mismatch is disclosed and the original saved statistics are preserved. The current code cannot resolve whether discrepancies come from coordinate/host-selection conventions or track versions. In these reconstructed bins, 25/30 vs 6/33 calls have VCF TSD strings; their median reported lengths are 16 and 15 bp. `POLYA_MIN_BP` medians are 84 and 91 respectively. This is descriptive only; it cannot test TPRT triad frequency or distinguish genuine tail peak insertions from artifact.

Of 320 nested calls, the RepeatMasker reconstruction matches the saved host name/strand/length for 312; 6 have no matching interval at the normalized coordinate and 2 have a same-family interval but a different selected record. These 8 are retained and explicitly flagged in the per-call audit. The per-call CSV additionally contains the VCF's call tier, score, known-source and SR/DPE read-support annotations for manual triage.

SVA per-call triage is especially limited: 9 nested SVA calls; 7 have reported TSD strings, of lengths 4, 9, 11, 19, 21 (2 calls), and 28 bp; 2 are missing. Seven calls have nonzero classifier `POLYA_MIN_BP`, but this field is not a sequence-level tail reconstruction. The 9 rows and their available support values are in the JSON summary and CSV. These are not 9 individually validated TPRT insertions.

### Why the “TPRT triad now fully specified” claim is not test-ready

Three kinds of evidence would have to be independently reconstructed for each event: (1) target-site duplication from both original flanks and correct insertion boundaries, (2) the insert's 3′ A-rich sequence, and (3) a strand-compatible L1 EN target motif at the integration nick. The current VCF stores a single `TSD` candidate and minimum polyA classifier feature, not the underlying target sequence or separate flanking junctions; it lacks the motif sequence. It also lacks read alignments needed to decide whether either junction is a split/truncated, complex, or boundary-ambiguous call. The Claude export shows the full hg38 2bit was downloaded in a different prior workspace, but it does not embed the binary; the present checkout has no reference FASTA/2bit, HG03086 BAM/CRAM, or callable/mappability track. The Nummi/Wagstaff rates are useful literature context but are not a calibration transferable to these inherited classifier records. We do not call the 96% Wagstaff de novo Alu hallmark estimate a validation threshold for this cohort.

Wagstaff et al. 2012 report 226 de novo tagged Alu insertions in HeLa; 96% of fully characterized events carried the *combined* hallmark set (direct repeats, 3′ oligo-dA and L1 EN-like target). That is an experimental cohort and a combined hallmark definition, not an individual-call classifier sensitivity/specificity benchmark. Nummi 2025 provides more recent long-read germline and somatic distributions, but again with their pipeline's event reconstruction and filters. Neither justifies inferring authenticity from a TSD length within range plus a `POLYA_MIN_BP` integer.

## 4. Literature on mechanisms and comparability

**Jacob-Hirsch et al. 2018** (*Cell Research* 28:187–203, doi:10.1038/cr.2018.8) analyzed whole-genome sequencing of 20 brain samples and 80 non-brain samples (not single neurons). They report most somatic brain insertions nested in existing repeats, and a specific class of nested somatic L1 insertions with same-host orientation, both-end truncation, and CCATT-rich targets consistent with an EN-independent route. This is not a numeric comparator for our inherited mixed-family events. Our previous CCATT test did not find enrichment vs a composition-matched same-family element control (33 motif occurrences in 16 kb vs 574 in 320 kb; rate ratio 1.15, p=0.445; [selfins_ccatt_control.csv](selfins_ccatt_control.csv)). That supports “no detected enrichment under this test,” not a proven classical mechanism.

**Lavi & Carmel 2018** (*RNA Biology* 15:715–725, doi:10.1080/15476286.2018.1429880) analyze Alu-derived transcript cleavage/polyadenylation sites, not genomic retrotransposition insertion breakpoints. Their data show apparent tail-adjacent expression-site peaks at 281/283 and 300 strongly reduced or erased in APADB, consistent with internal-priming artifacts in transcript-end annotation; their linker-related position 119 itself is complicated by A-rich sequence and cleavage-site assignment. This is a caution about expression-derived coordinates, not evidence that a DNA insertion hotspot is artifact. It was published in 2018, after Levy's insertion-site study, so it is not the positional comparator Levy should have excluded.

**Damert et al. 2009** concern SVA 5′-flanking sequence transduction via transcription/splicing and subsequent retrotransposition; this is not a nested genomic insertion of one SVA into a reference SVA. But “no precedent” should not be asserted for SVA–SVA target nesting, because Nummi et al. 2025 now report SVA self-target preference directly.

## 5. Conclusions safe to send to William

1. In the HG03086 classifier catalogue, same-family nesting exceeds the selected genome-placement expectations for Alu and LINE-1 (both 1.47×); SVA is 9 vs 0.306 expected and must remain provisional pending repeat-aware detectability control.
2. The nested set shows same-orientation excess for Alu and LINE-1, with SVA too small for a stable mechanistic claim. Report the conditional fractions and the simulated class-vs-placement tests as distinct analyses.
3. The 120–140 bp Alu peak is compatible with Levy et al.'s hotspot at consensus position ~133 at our 20-bp resolution. Nummi 2025 reports a similar Alu region at 100–155 bp. Our 280–300 bp tail peak is independently compatible with Nummi's 260–310 region, but evidence of its event-level authenticity is pending.
4. The 1,626/1,852 reported TSD strings and `POLYA_MIN_BP` values support a reproducible **QC audit**, not the requested triad validation. No sequence-backed TPRT triad percentages or CCATT-vs-classical route assignment should be published for these calls from the current files.
5. Treat every result as the pattern of catalogue polymorphisms surviving demographic and purifying selection, not as an unfiltered integration-rate experiment.

## References checked

- Levy A, Schwartz S, Ast G. 2009 online / 2010 issue. *Nucleic Acids Res* 38:1515–1530. doi:10.1093/nar/gkp1134. [PMC full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC2836564/).
- Szak ST et al. 2002. *Genome Biol* 3:research0052. doi:10.1186/gb-2002-3-10-research0052. [PMC full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC134481/).
- Zingler N et al. 2005. *Genome Res* 15:780–789. doi:10.1101/gr.3421505. [PMC full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC1142468/).
- Wagstaff BJ et al. 2012. *PLoS Genet* 8:e1002842. doi:10.1371/journal.pgen.1002842. [PMC full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC3415434/).
- Jacob-Hirsch J et al. 2018. *Cell Res* 28:187–203. doi:10.1038/cr.2018.8. [PMC full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC5799824/).
- Lavi E, Carmel L. 2018. *RNA Biol* 15:715–725. doi:10.1080/15476286.2018.1429880. [PMC full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC6152444/).
- Damert A et al. 2009. *Genome Res* 19:1992–2008. doi:10.1101/gr.093435.109. [PMC full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC2775593/).
- Hancks DC et al. 2009. *Genome Res* 19:1983–1991. doi:10.1101/gr.093153.109. [PMC full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC2775590/).
- Nummi P et al. 2025. *Mobile DNA* 16:20. doi:10.1186/s13100-025-00357-w. [PMC full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC12016303/).
