---
name: real-locus-test-fixtures
description: >-
  Build pytest fixtures from real BAM/CRAM snippets and extract tables for a
  single MEI locus. Use when adding or fixing locus-specific tests, breakpoint
  pickers, sentinel sites, nssv/gold variants, support-count or gold rules, or
  any test that would otherwise invent cartoon split/DPE piles.
---

# Real locus test fixtures

Do not invent read piles for a named locus. Extract the real BAM/CRAM snippet
and the extract tables that the pipeline actually saw.

When a feature is developed on sentinel examples, replace the cartoon tests
with real sentinel-locus fixtures. Do not keep a parallel hand-built DataFrame
of clips and DPE for the same rule. This is also
`.cursor/rules/real-locus-fixtures.mdc` (always apply).

## When this applies

- A test claims to cover a catalog site (`nssv*`, sentinel, gold row, IGV locus)
- Breakpoint / clustering / support-count / gold-stage logic is being changed
- A cartoon unit test passed but the chr22 (or other) rerun still missed the site
- Orientation, two-sided support, parked DPE, polyA-side, or clip-to-flank rules

## Do not

- Invent a pile of 4 clips and 17 DPE and call that a named locus
- Keep cartoon scoring-rule tests once a real sentinel demonstrates the rule
- Treat “this is just a scoring rule” as permission to skip extraction

Cartoon frames are only allowed when no real locus yet shows the behavior.
Extract one as soon as a sentinel, gold row, or IGV site exists.

## Required fixture

Put files under `tests/fixtures/loci/<catalog_id>/`:

- `manifest.json` — chrom, expected breakpoint, discovery window, BAM name,
  `plot_rank`, `expect_gold` when the test is a gold/support rule
- `candidate.tsv` — picker input window (merged discovery window, not the published BP)
- `gold_locus.tsv` — the `candidate_loci.mei.tsv` row for this window (gold tests)
- `split_evidence.{disease,control}.parquet`
- `discordant_evidence.{disease,control}.parquet` with annotate-time `mei_hit` flags
- `supporting_reads_detail.mei.tsv` slice
- `<sample>.<chrom>_<fetchStart>_<fetchEnd>.bam` + `.bai` (region plus mates)

`.gitignore` ignores `*.bam` globally. Keep this exception:

```
!tests/fixtures/**/*.bam
!tests/fixtures/**/*.bai
```

## How to cut the snippet

From an existing results dir and the source CRAM, inside `rtm-miner`:

```bash
python scripts/extract_locus_fixture.py \
  --results-dir "$RESULTS" \
  --cram "$CRAM" \
  --reference-fasta "$FASTA" \
  --chrom chr22 \
  --breakpoint 49879732 \
  --window-start 49878612 \
  --window-end 49880399 \
  --catalog-id nssv14073986 \
  --outdir tests/fixtures/loci/nssv14073986
```

Fetch window = discovery window ± 800 bp. Include discordant mates. Sort + index.

Extract parquets often lack annotate-time `mei_hit`. Copy MEI flags from
`supporting_reads_detail.mei.tsv` onto the discordant rows by `read_name`.

Copy the matching `candidate_loci.mei.tsv` row to `gold_locus.tsv` when the
test assigns gold or counts flank support.

## How to write the test

Load the fixture and run the real function (`_choose_window_breakpoints`,
`_genomic_flank_evidence_table`, `_assign_gold_stage`, extract, etc.).
Assert the catalog breakpoint or the gold keep/drop. Also assert the BAM still
contains the junction read named in the detail table.

Map each scoring assertion onto a real sentinel that has that geometry
(one-sided pile, parked DPE, MEI+polyA, two-sided keep). Do not also keep a
synthetic DataFrame that restates the same rule.

## Check the fixture actually fails the old bug

Replay the fixture through the current picker or gold assigner **before**
calling the test done. If the fixture publishes the same wrong coordinate or
the same one-sided gold as chr22, the test is good. If it publishes the
catalog site / gold pass while chr22 does not, the fixture is still a cartoon
— missing leftover clips, window tags, or DPE flags.
