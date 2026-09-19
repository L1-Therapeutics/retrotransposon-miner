---
name: real-locus-test-fixtures
description: >-
  Build pytest fixtures from real BAM/CRAM snippets and extract tables for a
  single MEI locus. Use when adding or fixing locus-specific tests, breakpoint
  pickers, sentinel sites, nssv/gold variants, or any test that would otherwise
  invent cartoon split/DPE piles.
---

# Real locus test fixtures

Do not invent read piles for a named locus. Extract the real BAM/CRAM snippet
and the extract tables that the pipeline actually saw.

## When this applies

- A test claims to cover a catalog site (`nssv*`, sentinel, gold row, IGV locus)
- Breakpoint / clustering / support-count logic is being changed
- A cartoon unit test passed but the chr22 (or other) rerun still missed the site

## Required fixture

Put files under `tests/fixtures/loci/<catalog_id>/`:

- `manifest.json` — chrom, expected breakpoint, discovery window, BAM name
- `candidate.tsv` — picker input window (merged discovery window, not the published BP)
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

## How to write the test

Load the fixture and run the real function (`_choose_window_breakpoints`, extract,
etc.). Assert the catalog breakpoint. Also assert the BAM still contains the
junction read named in the detail table.

Do not replace this with a hand-built DataFrame of 4 clips and 17 DPE unless
that test is explicitly about a scoring rule, not a named locus.

## Check the fixture actually fails the old bug

Replay the fixture through the current picker **before** calling the test done.
If the fixture publishes the same wrong coordinate as chr22, the test is good.
If it publishes the catalog site while chr22 does not, the fixture is still a
cartoon — missing leftover clips, window tags, or DPE flags.
