# Gene annotation: VEP vs snpEff (decision note)

**Summary.** Both tools annotate all 30 chr22 MEI calls correctly once configured. VEP REST needs
no install and no local data but is network-bound (0.3–0.6 s/variant, ~1–2 h for a ~10k-call
genome). snpEff needs a 610 MB database and an `-Xmx8g` heap, pays ~55 s of fixed startup, then
annotates ~30 variants/s (~6–7 min for the same genome) and pins the Ensembl release, so results
are reproducible offline. `rtm annotate-genes` ships both behind `--backend {vep,snpeff}` with an
identical output contract; **`vep` is the default for zero-setup use, `snpeff` is recommended for
production and whole-genome runs.**

Request (W. Brandler, 2026-09-18): evaluate Ensembl VEP and snpEff for annotating MEI
calls, decide on **ease of use** and **runtime**, and test on the chr22 output VCF.
Both tools permit commercial use (VEP: Apache-2.0; snpEff: MIT).

Test set: `tests/data/chr22_mei.vcf` — the 30 GRCh38 chr22 calls in the README example
table, exported by `vcf_export` (symbolic `<INS:ME:*>` ALT, `SVTYPE=INS`, `SVLEN`, `END=POS`,
`MEINFO`). All numbers below were measured on 2026-09-18; anything not measured is marked.

## Measurements

| Criterion | Ensembl VEP (REST, rest.ensembl.org) | snpEff 5.4c (local, `GRCh38.99`) |
|---|---|---|
| Install | none (HTTPS only) | conda `snpeff` + `openjdk` (40.7 MB) — fine |
| Reference data | none locally | `snpEff download GRCh38.99`: 610 MB. Silent no-op on first two attempts; with `-v`, mirrors `v5_4`–`v5_1` all fail to connect, `v5_0` succeeds |
| Startup | none | 35 s database load + 7 s interval forest ≈ 55 s before the first variant (`-Xmx8g`) |
| Memory | none | default JVM heap → `OutOfMemoryError` building the GRCh38 interval forest; `-Xmx8g` works. `annotate-genes` passes `-noHgvs` so snpEff does not also load the reference sequence for protein notation |
| `chr22` vs `22` contig naming | accepted as-is (verified) | mapped to database contig `22`, echoed back as `chr22` (verified) |
| Symbolic `<INS:ME:*>` ALT | accepted; needs SVLEN dropped (see below) | evaluated as a point natively; 6 kb L1 at 19223382 → `intron_variant`, CLTCL1, 5 transcripts |
| Time, 30 variants | 9.4 s and 18.6 s on two runs (0.31–0.62 s/var); 200-variant POST batches | 57.3 s wall, of which ~1 s is annotation (log resolution 1 s) |
| chr22 result | 30/30 annotated, 24 in genes, 0 coding | 30/30 annotated, 23 in genes, 0 coding |
| Annotation source | Ensembl release 116 (live) | `GRCh38.99` = Ensembl release 99 (2020); `GRCh38.115` also available |

Extrapolation from the measured rates (VEP 0.3–0.6 s/var; snpEff 55 s + ~1 s per 30):

| Callset | VEP REST | snpEff |
|---|---|---|
| 1,000 calls (one chromosome) | ~5–10 min | ~1.5 min |
| ~10,000 calls (genome-wide MEIs) | ~1–2 h | ~6–7 min |

### Concordance on the 30 chr22 loci

Most-severe term agreed at 17/30 when snpEff's *first* `ANN` entry was taken at face value, but
snpEff does not order `ANN` by Ensembl severity (it lists `upstream_gene_variant` ahead of
`intron_variant` for the same gene). Re-ranking snpEff's full `ANN` list with the same severity
table raises agreement to **25/30**. The remaining 5 loci are all places where VEP (release 116)
reports an intron of a lncRNA or recently added gene that release 99 does not contain
(e.g. `ENSG00000308779` at 17289460, `CACNG2-DT` at 36746494/36752165). That is a database-age
difference, not a tool difference; using `GRCh38.115` should close most of it.
`parse_snpeff_ann()` applies the severity ranking so both backends report the same `CSQ`.

## Decision

Ship both behind one command with one output contract:

- `--backend vep` (default): zero install, right answer for small callsets, notebooks, and CI.
- `--backend snpeff`: recommended for production and whole-genome runs — ~10× faster at scale,
  no network dependency or rate limit, and the Ensembl release is pinned by the database name,
  which makes annotations reproducible. Cost: one-time 610 MB download and `-Xmx8g`.

Open item for the reviewer: which snpEff database to standardise on (`GRCh38.115` for release
parity with current Ensembl, or a `GRCh38.mane.*.ensembl` build for one canonical transcript per
gene). The CLI takes it as `--snpeff-genome`, so this is a config decision, not a code change.

## The SVLEN finding (VEP only; snpEff is unaffected)

For a symbolic insertion VEP computes the affected span as `POS + SVLEN` and ignores `END`.
chr22:19223382 `<INS:ME:LINE1>` with `SVLEN=6018`:

| Request INFO | VEP parsed span | `coding_sequence_variant` transcripts | most severe |
|---|---|---|---|
| `SVTYPE=INS;SVLEN=6018;END=19223382` | 19223383–19229400 | 38 | `feature_elongation` |
| `SVTYPE=INS;END=19223382` | 19223383–19223382 | 0 | `intron_variant` |
| `SVTYPE=INS` | 19223383–19223383 | 0 | `intron_variant` |
| `SVTYPE=INS;SVLEN=281;END=…` (Alu-sized) | 19223383–19223663 | 0 | `intron_variant` |

Across the 30 chr22 calls, 2/30 changed most-severe consequence (19223382 L1 and 34034616 Alu,
both `feature_elongation` → `intron_variant`). An insertion adds sequence *between* two
reference bases; it does not alter the reference bases downstream, so the widened span is
wrong. `vep_region_string()` therefore drops `SVLEN` from the request and pins `END=POS`.
`SVLEN` stays in the output VCF because it is correct VCF (v4.4 §5.6) — the issue is
VEP's interpretation, not the exporter.

## References

- McLaren W. et al. (2016) The Ensembl Variant Effect Predictor. *Genome Biology* 17:122.
  doi:10.1186/s13059-016-0974-4
- Cingolani P. et al. (2012) A program for annotating and predicting the effects of single
  nucleotide polymorphisms, SnpEff. *Fly* 6(2):80–92. doi:10.4161/fly.19695
- Eilbeck K. et al. (2005) The Sequence Ontology. *Genome Biology* 6:R44. doi:10.1186/gb-2005-6-5-r44
- Ensembl, "Calculated variant consequences" (severity table):
  https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html
- VCF specification v4.4 §5.6 (symbolic SV alleles): https://samtools.github.io/hts-specs/
