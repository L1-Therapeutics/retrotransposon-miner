# Splice-context orientation densities — what we have and what we do not

## Re-licensing

`_gene_map/splice/alu_splice_distance.py`, `_gene_map/splice/test_alu_splice_distance.py`,
and this README are now MIT-licensed. Metadata added:

- `License: MIT`
- `Copyright: Smyan Reddy, 2026`
- `Classifier score threshold: ≥ 0.997`

The `alu_splice_window.png` figure is not re-licensed and stays in this folder as-is.

## What this is

`_gene_map/splice/alu_splice_distance.py` classifies Alu insertions by where they
sit relative to GENCODE exons and the host gene strand:

- **in_exon** — overlap an exon
- **within_150_upstream / within_150_downstream** — intronic, within 150 bp of an exon,
  on the transcription-upstream or transcription-downstream side
- **intronic_far** — all other intronic copies

It then tests antisense fraction in the near-exon window against the far-intronic baseline
with a Fisher exact test. That is a splice-adjacent density/orientation test, not a
direct splice-regulatory-element density scan.

## What is already on disk

Reference-copy results exist for three subfamilies:

- `_gene_map/splice/out_aluy_ref/`
- `_gene_map/splice/out_alusx_ref/`
- `_gene_map/splice/out_alujo_ref/`

Each contains `alu_splice_sites.tsv`, `alu_splice_zones.tsv`, and `alu_splice_tests.tsv`.
A polymorphic reproduction is in `out_polymorphic_repro/`, but it is small (55 upsteam
sites) and not informative for the near-exon window. A re-run of the AluY pipeline is in
`out_aluy_ref_repro/` and reproduces `out_aluy_ref/` to the published digits.

## Plain-language summary

### Exon overlap

Exonic Alu copies are antisense-majority but close to even:

- AluY: 3150 exonic, 0.539 antisense
- AluSx: 4046 exonic, 0.576 antisense
- AluJo: 3163 exonic, 0.576 antisense

This is not the strong antisense bias you see at splice junctions in the literature; in our
GENCODE-resolved set, "exonic" means overlap an exon in the reference copy coordinates, which
is a broader and noisier category than the curated exonized-Alu junction sets in the
exonization literature.

### The replicated near-exon signal

The one clean, replicated splice-context result is **antisense depletion within 150 bp
upstream of an exon**:

| family | n_near | antisense_near | antisense_frac_near | OR vs far intronic | p |
|---|---|---|---|---|---|
| AluY | 1520 | 735 | 0.484 | 0.758 | 9.7e-08 |
| AluSx | 1903 | 973 | 0.511 | 0.818 | 1.8e-05 |
| AluJo | 1182 | 565 | 0.478 | 0.798 | 1.4e-04 |

Baseline far-intronic antisense fraction is about 0.553 for AluY, 0.561 for AluSx, 0.534 for
AluJo, so the near-window drop is real and consistent across subfamilies.

The downstream-150bp window is weaker:

- AluY OR 0.897, p=0.033
- AluSx OR 0.950, p=0.258
- AluJo OR 0.948, p=0.349

That is the same upstream-strong / downstream-weak shape Lev-Maor reported.

### Pooled upstream+downstream window

Pooling the upstream and downstream 150 bp windows into one "near-exon" band is a weaker,
broader signal than the upstream-only window that carries the result:

| family | n_near pooled | antisense_frac pooled | OR pooled vs far | χ² | p |
|---|---|---|---|---|---|
| AluY | 3099 | 0.505 | 0.826 | 27.1 | 1.9e-07 |
| AluSx | 3949 | 0.530 | 0.884 | — | — |
| AluJo | 2464 | 0.500 | 0.872 | — | — |

For AluY the pooled window is real (p = 1.9e-07) but the effect size is a comoving average
of the two one-sided windows, not a single homogeneous near-exon signal. The upstream side
is where the depletion lives. Wherever this README or downstream prose says "near exon", it
should say "within 150 bp upstream of an exon" unless it explicitly means pooled.

### Intron half

Looking at intron halves in transcription coordinates:

- 5′ intron half antisense: 0.543
- 3′ intron half antisense: 0.557

That is a small overall difference; the real structure is in the 150 bp near-exon band, not
a genome-wide intron-half shift.

### Distance gradient

Upstream antisense fraction by distance decile does not show a simple monotonic ramp. The
closest decile is the lowest at 0.497, but the nearby deciles stay around 0.54–0.56, with a
dip again around deciles 8–10. The point is not a smooth density gradient; it is a narrow
exclusion band immediately upstream of the exon.

## What this analysis is not

It is not a splice-regulatory-element density map. It does not count:

- ESE / ESS motifs
- splice-site strength
- splice-factor binding density
- exonization-competent motifs in the Alu right arm

It only asks: given GENCODE exon boundaries, where do copies land, and in which orientation
relative to the host gene. That is enough to reproduce the Lev-Maor antisense-depletion-near-exon
result, but not enough to say much about splice-regulatory density per se.

## Reproducibility note

Tests pass:

- `_gene_map/splice/test_alu_splice_distance.py` — 24 passed
- `rtm-nested-analysis/tests/test_servicenow.py`, `tests/test_nested_null.py`,
  `tests/test_http_client.py` — 15 passed

So the classification and test code is internally consistent. The biological interpretation
still depends on the GENCODE annotation version and the exon-union construction in this script,
not on the retrotransposon-miner pipeline directly.
