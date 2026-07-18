# retrotransposon-miner

Short-read mobile element insertion (MEI) detection and annotation pipeline for LINE-1, Alu, and SVA.

`retrotransposon-miner` searches next-generation sequencing data for retrotransposon (mobile element insertion, MEI) events in the human genome.
Retrotransposons are virus-like elements that can activate under stress and are implicated in disease biology.

The pipeline detects multiple MEI classes and outputs a candidate insertion table annotated with evidence and context useful for triage.

## What This Tool Can Do

- Detect mobile element insertion (MEI) candidates from short-read data using split-read and discordant paired-end evidence.
- Support paired disease/control (tumor-normal) workflows and control-focused/germline style analyses.
- Call major retrotransposon classes (`LINE-1`, `Alu`, `SVA`) and report family/subfamily assignments.
- Annotate candidates with:
  - estimated insertion coordinates,
  - polyA/polyT support,
  - target site duplication (TSD) length/sequence when resolvable,
  - overlap with known variant resources (short-read and long-read sets),
  - optional local assembly-derived features.
- Auto-generate review snapshots in IGV (Integrative Genomics Viewer) as PNG files for manual QC.
- Run cleanly in Linux environments with included bootstrap/validation scripts and an Amazon Elastic Compute Cloud (EC2)-first workflow.

## How It Compares to Other Tools

The table below summarizes `retrotransposon-miner` against commonly used tools (`xTea`, `mobster`, `MELT`, `RetroNet`, `TraFiC`, `TotalReCall`, `MEIba`).  
`retrotransposon-miner` feature claims are based on this repository; other-tool columns are high-level, publicly documented capability summaries and may vary by version/workflow.

Legend: `✅` yes, `❌` no, `➖` limited/partial/not definitive.

| Feature | ![L1 Therapeutics](https://avatars.githubusercontent.com/l1-therapeutics?s=40) | xTea | mobster | MELT | RetroNet | TraFiC | TotalReCall | MEIba |
|---|---|---|---|---|---|---|---|---|
| Split-read support | ✅ | ✅ | ✅ | ✅ | ✅ | ➖ | ➖ | ➖ |
| Discordant paired-end support | ✅ | ✅ | ✅ | ✅ | ➖ | ✅ | ✅ | ➖ |
| Germline analysis | ✅ | ✅ | ✅ | ✅ | ✅ | ➖ | ➖ | ✅ |
| Paired disease/control analysis | ✅ | ✅ | ➖ | ➖ | ➖ | ✅ | ✅ | ➖ |
| Alu / SVA / LINE-1 | ✅ | ✅ | ✅ | ✅ | ✅ | ➖ | ➖ | ➖ |
| hg19 / GRCh38 / hs1 support | ✅ | ➖ | ➖ | ➖ | ➖ | ➖ | ➖ | ➖ |
| Target site duplication (TSD) detection | ✅ | ➖ | ➖ | ➖ | ➖ | ➖ | ➖ | ➖ |
| PolyA/polyT characterization | ✅ | ✅ | ✅ | ➖ | ✅ | ✅ | ✅ | ➖ |
| Annotation of known variant catalogs | ✅ | ➖ | ➖ | ✅ | ➖ | ➖ | ➖ | ➖ |
| Optional local assembly | ✅ | ✅ | ❌ | ❌ | ❌ | ❌ | ❌ | ❌ |
| Cloud support | ✅ | ➖ | ➖ | ➖ | ➖ | ➖ | ➖ | ➖ |
| Linux environment management scripts | ✅ | ➖ | ➖ | ➖ | ➖ | ➖ | ➖ | ➖ |
| Integrative Genomics Viewer (IGV) + JupyterLab workflow | ✅ | ❌ | ❌ | ❌ | ❌ | ❌ | ❌ | ❌ |

## Current Limitations

- Designed primarily for Amazon Web Services (AWS) machines today; relatively straightforward to adapt to Google Cloud Platform (GCP), Azure, or local Linux. For larger runs, `r6i.4xlarge` or greater is recommended, and local assembly/candidate processing support parallel execution.
- Artificial intelligence/machine learning (AI/ML) genotyping confidence models are still under active development.
- Reverse-transcribed pseudogene insertion support is not yet added.
- Support for species other than *Homo sapiens* (for example, *Mus musculus*) is not yet implemented.
- Long-read native calling is not yet supported.
- Single-cell sequencing data is not yet supported.
- Other short-read platforms (for example, Ultima Genomics) are not yet validated/supported.
- Local assembly is parallelized but still compute-expensive; it is optional and not recommended by default for routine runs.

## MEI Consensus Remapping (Why `bwa mem`)

Split-read soft-clips and discordant clipped ends are remapped to the Dfam Alu / LINE-1 / SVA consensus FASTA with:

```bash
bwa mem -t <N> -k 10 -T 10 <mei.fasta> <queries.fa>
```

`<N>` comes from `annotate-mei-support --bwa-threads` (CLI default **1**). The pipeline wrapper sets it automatically: `nproc` for single-chrom / serial chroms, and **1** under multi-chrom concurrency (`--chr_concurrency > 1`). Override anytime with `--bwa-threads N` on the CLI or wrapper.

Before remap, longest A/T homopolymers (≥8 bp) are trimmed from the query so polyA+tip clips are scored on the tip. Hits then pass a length-aware gate: short trimmed queries (≤30 bp) need qcov ≥ 0.80 and pid ≥ 0.90; longer queries may tip-align (pid ≥ 0.90 and alnlen ≥ 20) without requiring full-query coverage. Primary and coord hit fields share this single trimmed+gated alignment (no second bwa pass). SAM parsing keeps primary alignments only and scores qcov against the full CIGAR query length (so hard-clipped supplementaries cannot fake qcov=1).

Local assembly contig-to-MEI alignment still uses `minimap2` (longer sequences).

Panel hit coordinates (e.g. `L1HS_5end`, `L1HS_3end`) are projected onto a shared full-length consensus axis at annotate time using a one-time prep table (`mei_fragment_to_full_coords.tsv`). That table is built with sensitive `bwa mem -a` of each Dfam fragment onto the full-consensus panel (prefer `{subfamily}_full`, else family canonical). Gold unions never mix Alu and LINE-1 axes (or different `*_full` targets).

### Why not default `minimap2 -x sr`?

Default short-read minimap2 uses ~21 bp seeds. Perfect 20–30 bp Alu / L1 / SVA tips therefore often produce **no hit**, so they never enter `SR_*` / `MEI_MAPPED` even when the clip is an exact consensus substring. That under-counts real MEI split support relative to polyA (which is sequence-rule based from ~8 bp).

### Benchmark (reproducible)

`scripts/benchmark_mei_aligners.py` samples random substrings (and reverse complements) from every Alu, LINE-1, and SVA entry in `dfam_human_mei_l1_alu_sva.fasta`, then realigns them with each tool.

Recent run on this machine (840 queries = 40 tips × 3 families × 7 lengths `{20,22,25,30,40,60,100}` × ~50% RC; plus a 5000 × 40 bp timing batch):

| method | overall hit% | 20 bp Alu / L1 / SVA | 5000×40 bp wall time | notes |
|---|---:|---|---:|---|
| `minimap2 -x sr` (old default) | 38% | 0% / 0% / 0% | **0.045 s** (fastest) | misses all ≤30 bp tips |
| `minimap2 -x sr -k10 -w5` | 52% | 0% / 0% / 0% | 0.27 s | helps ≥30 bp only |
| `bwa mem` (defaults) | 56% | 0% / 0% / 0% | 0.19 s | seed `k=19` too long for 20 bp |
| **`bwa mem -k10 -T10` (chosen)** | **99.9%** | **100% / 100% / 100%** | 0.25 s | best short-tip recall; stays ~100% at 100 bp |
| `bwa aln -l10` | 98% | 100% / 98% / 100% | **0.07 s** | great at 20 bp; drops on longer L1/SVA (~85–95% at 100 bp) |
| `bowtie2 --very-sensitive-local` | 83% | 0% / 0% / 0% | 0.17 s | good from ~22–30 bp; misses 20 bp |

**Runtime:** on this small MEI database, `bwa aln` is slightly faster than `bwa mem -k10 -T10`, and both are far cheaper than annotation I/O. We still prefer **`bwa mem -k10 -T10`** because (1) it matches `aln` on 20 bp tips for Alu **and** LINE-1 **and** SVA, (2) it does not lose long L1/SVA fragments the way `aln` does, and (3) it is a single-pass aligner (no `aln`+`samse` staging). Bowtie2 is competitive for ≥22 bp but not for the 20 bp SR gate.

Re-run:

```bash
conda activate rtm-miner || micromamba activate rtm-miner
PYTHONPATH=src python scripts/benchmark_mei_aligners.py \
  --mei-fasta "${RTM_PUBLIC_DATA_DIR:-$HOME/retrotransposon-workdir/data/public}/retrotransposon_db/dfam/dfam_human_mei_l1_alu_sva.fasta"
```

## Example Variant Calls (GRCh38)

The table below lists candidate insertion calls from a tumor/normal chr22 run.  
In the example table, `SR` denotes split-read evidence and `DPE` denotes discordant paired-end evidence.

### Example output from sample tumor/normal data

Gold-tier calls from the SEQC2 tumor/normal chr22 annotate run with `bwa mem` MEI remap + polyA/T trim/gate, family-safe full-axis projection, and BWA-built fragment→full map (top 25 of n=1335 by `read_support_heuristic_score`; top 100 family mix ≈ 85% Alu / 14% LINE1 / 1% SVA).

| chrom | consensus_insertion_breakpoint_pos | window_start | window_end | control_supporting_reads | disease_supporting_reads | sample_status_label | consensus_tsd_seq | consensus_poly_at_min_bp | consensus_mei_family | consensus_mei_subfamily | known_mei_polymorphism_id | known_mei_polymorphism_source | consensus_insertion_orientation | nested_in_same_MEI | consensus_insertion_mei_span_full | consensus_insertion_mei_5p_coord_full | consensus_insertion_mei_3p_coord_full |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| chr22 | 49029650 | 49029645 | 49029656 | SR_L=11,SR_R=0,DPE_L=2,DPE_R=152,MEI_MAPPED=154,polyA_MAPPED=46,VNTR_MAPPED=0,polyA_side=L | SR_L=32,SR_R=0,DPE_L=8,DPE_R=373,MEI_MAPPED=372,polyA_MAPPED=92,VNTR_MAPPED=1,polyA_side=L | shared | AAGAAAACTCCT | 50 | SVA | SVA_D#Retroposon/SVA | nssv14064350 | melt_1kg | + | unnested | 1380 | 1 | 1380 |
| chr22 | 37529127 | 37529108 | 37529146 | SR_L=0,SR_R=2,DPE_L=151,DPE_R=2,MEI_MAPPED=66,polyA_MAPPED=47,polyA_side=L | SR_L=0,SR_R=6,DPE_L=365,DPE_R=16,MEI_MAPPED=171,polyA_MAPPED=86,polyA_side=L | shared | GAAGCGGAGGTTGCAGTGAGCCGAGATTGCGCCACTGCA | 90 | ALU | AluYb8#SINE/Alu |  |  | + | nested | 318 | 1 | 318 |
| chr22 | 20075438 | 20075432 | 20075444 | SR_L=18,SR_R=23,DPE_L=128,DPE_R=19,MEI_MAPPED=156,polyA_MAPPED=51,polyA_side=L | SR_L=8,SR_R=17,DPE_L=89,DPE_R=16,MEI_MAPPED=93,polyA_MAPPED=40,polyA_side=R | shared | AGATTTCTTTTCT | 39 | ALU | AluYk12#SINE/Alu |  |  | - | unnested | 311 | 1 | 311 |
| chr22 | 50495046 | 50494596 | 50495497 | SR_L=2,SR_R=0,DPE_L=59,DPE_R=82,MEI_MAPPED=71,polyA_MAPPED=34,polyA_side=R | SR_L=8,SR_R=1,DPE_L=102,DPE_R=174,MEI_MAPPED=141,polyA_MAPPED=49,polyA_side=R | shared |  | 14 | ALU | AluYa5#SINE/Alu |  |  | - | unnested | 310 | 2 | 311 |
| chr22 | 45595784 | 45595639 | 45595930 | SR_L=0,SR_R=0,DPE_L=384,DPE_R=2,MEI_MAPPED=125,polyA_MAPPED=106,polyA_side=L | SR_L=0,SR_R=0,DPE_L=452,DPE_R=0,MEI_MAPPED=130,polyA_MAPPED=102,polyA_side=L | shared |  | 91 | ALU | AluJb_short_#SINE/Alu |  |  | - | unnested | 293 | 10 | 302 |
| chr22 | 31355872 | 31355858 | 31355887 | SR_L=12,SR_R=0,DPE_L=7,DPE_R=148,MEI_MAPPED=110,polyA_MAPPED=6 | SR_L=15,SR_R=0,DPE_L=19,DPE_R=273,MEI_MAPPED=187,polyA_MAPPED=16,polyA_side=R | shared | CCGCCTCGGCTTCCCAAAGTGCTGGGATTA | 71 | ALU | AluY_short_#SINE/Alu |  |  | - | nested | 311 | 1 | 311 |
| chr22 | 34034616 | 34034610 | 34034623 | SR_L=14,SR_R=10,DPE_L=62,DPE_R=1,MEI_MAPPED=73,polyA_MAPPED=14,polyA_side=R | SR_L=22,SR_R=27,DPE_L=92,DPE_R=2,MEI_MAPPED=119,polyA_MAPPED=32,polyA_side=R | shared | CAAATGGAACTTTT | 64 | ALU | AluYb8#SINE/Alu | nssv14071620 | melt_1kg | - | unnested | 309 | 1 | 309 |
| chr22 | 45166725 | 45166719 | 45166731 | SR_L=11,SR_R=8,DPE_L=28,DPE_R=3,MEI_MAPPED=30,polyA_MAPPED=14,polyA_side=L | SR_L=47,SR_R=16,DPE_L=66,DPE_R=3,MEI_MAPPED=86,polyA_MAPPED=54,polyA_side=L | shared | AAAGAATTATGTC | 64 | ALU | AluYc#SINE/Alu | g1k:nssv14054938;lr:chr22-45651200-INS->s904290<s909202>s904291-125 | melt_1kg,long_read_1kg_ont_vienna | + | unnested | 30 | 270 | 299 |
| chr22 | 33132520 | 33132513 | 33132527 | SR_L=28,SR_R=7,DPE_L=44,DPE_R=7,MEI_MAPPED=69,polyA_MAPPED=31,polyA_side=L | SR_L=19,SR_R=7,DPE_L=76,DPE_R=31,MEI_MAPPED=109,polyA_MAPPED=23,polyA_side=L | shared | AAAAGTCATTATTAG | 56 | ALU | AluYg6#SINE/Alu | nssv14075885 | melt_1kg | + | unnested | 311 | 1 | 311 |
| chr22 | 41050312 | 41050276 | 41050348 | SR_L=3,SR_R=0,DPE_L=1,DPE_R=85,MEI_MAPPED=64,polyA_MAPPED=27,polyA_side=R | SR_L=1,SR_R=0,DPE_L=4,DPE_R=166,MEI_MAPPED=124,polyA_MAPPED=44,polyA_side=R | shared |  | 34 | ALU | AluSz#SINE/Alu |  |  | + | unnested | 312 | 1 | 312 |
| chr22 | 45872288 | 45872282 | 45872294 | SR_L=10,SR_R=0,DPE_L=30,DPE_R=8,MEI_MAPPED=24,polyA_MAPPED=35,polyA_side=L | SR_L=31,SR_R=11,DPE_L=41,DPE_R=35,MEI_MAPPED=65,polyA_MAPPED=66,polyA_side=L | shared | AAAAAAAAAAAAA | 41 | ALU | AluSg#SINE/Alu |  |  | + | nested | 27 | 280 | 306 |
| chr22 | 28994190 | 28994177 | 28994203 | SR_L=0,SR_R=0,DPE_L=13,DPE_R=5,MEI_MAPPED=34,polyA_MAPPED=28,polyA_side=R | SR_L=0,SR_R=28,DPE_L=41,DPE_R=10,MEI_MAPPED=69,polyA_MAPPED=66,polyA_side=R | shared | TTTTTTTTTTTTTTTTTTTTTTTTTTT | 47 | ALU | AluYb8#SINE/Alu |  |  | - | nested | 43 | 270 | 312 |
| chr22 | 47908556 | 47908553 | 47908560 | SR_L=0,SR_R=2,DPE_L=2,DPE_R=2,MEI_MAPPED=16,polyA_MAPPED=18,polyA_side=R | SR_L=1,SR_R=35,DPE_L=5,DPE_R=26,MEI_MAPPED=58,polyA_MAPPED=76,polyA_side=R | disease_only | TTTTTTTT | 28 | ALU | AluSz#SINE/Alu |  |  | - | unnested | 284 | 28 | 311 |
| chr22 | 38557968 | 38557957 | 38557980 | SR_L=4,SR_R=7,DPE_L=11,DPE_R=15,MEI_MAPPED=24,polyA_MAPPED=23,polyA_side=R | SR_L=7,SR_R=33,DPE_L=23,DPE_R=25,MEI_MAPPED=62,polyA_MAPPED=54,polyA_side=R | shared | TTTTTTTTTTTTTTTTTTTTTTTT | 43 | ALU | AluYb8#SINE/Alu |  |  | - | nested | 56 | 262 | 317 |
| chr22 | 19223382 | 19223373 | 19223390 | SR_L=17,SR_R=14,DPE_L=60,DPE_R=27,MEI_MAPPED=108,polyA_MAPPED=17,polyA_side=L | SR_L=0,SR_R=0,DPE_L=0,DPE_R=0,MEI_MAPPED=0,polyA_MAPPED=0 | control_only | AAAAACCACCTATGCTGG | 66 | LINE1 | L1HS_5end#LINE/L1 | g1k:nssv14064681;lr:chr22-19600083-INS->s899391<s914453>s899392-6059 | melt_1kg,long_read_1kg_ont_vienna | + | unnested | 6018 | 1 | 6018 |
| chr22 | 47882294 | 47882287 | 47882300 | SR_L=0,SR_R=5,DPE_L=1,DPE_R=7,MEI_MAPPED=12,polyA_MAPPED=13,polyA_side=R | SR_L=0,SR_R=38,DPE_L=5,DPE_R=27,MEI_MAPPED=58,polyA_MAPPED=67,polyA_side=R | shared | TTTTTTTTTTTTTT | 32 | ALU | AluY#SINE/Alu |  |  | - | nested | 41 | 266 | 306 |
| chr22 | 33878142 | 33878130 | 33878154 | SR_L=2,SR_R=0,DPE_L=7,DPE_R=21,MEI_MAPPED=38,polyA_MAPPED=32,polyA_side=R | SR_L=18,SR_R=20,DPE_L=11,DPE_R=32,MEI_MAPPED=62,polyA_MAPPED=50,polyA_side=R | shared | TTTTTTTTTTTTTTTTTTTTTTTTT | 56 | ALU | AluYb9#SINE/Alu |  |  | - | nested | 45 | 272 | 316 |
| chr22 | 42859704 | 42859696 | 42859712 | SR_L=17,SR_R=0,DPE_L=14,DPE_R=12,MEI_MAPPED=29,polyA_MAPPED=33,polyA_side=L | SR_L=36,SR_R=0,DPE_L=14,DPE_R=15,MEI_MAPPED=57,polyA_MAPPED=64,polyA_side=L | shared |  | 38 | ALU | AluSq4#SINE/Alu |  |  | + | unnested | 40 | 272 | 311 |
| chr22 | 47896792 | 47896640 | 47896943 | SR_L=0,SR_R=0,DPE_L=3,DPE_R=3,MEI_MAPPED=3,polyA_MAPPED=43,polyA_side=L | SR_L=20,SR_R=6,DPE_L=3,DPE_R=14,MEI_MAPPED=44,polyA_MAPPED=133,polyA_side=L | disease_only |  | 35 | ALU | 7SLRNA#SINE/Alu |  |  | + | unnested | 49 | 263 | 311 |
| chr22 | 34539961 | 34539493 | 34540429 | SR_L=15,SR_R=0,DPE_L=2,DPE_R=23,MEI_MAPPED=27,polyA_MAPPED=22,polyA_side=L | SR_L=30,SR_R=0,DPE_L=11,DPE_R=54,MEI_MAPPED=61,polyA_MAPPED=49,polyA_side=L | shared |  | 47 | ALU | AluJo#SINE/Alu |  |  | + | unnested | 33 | 280 | 312 |
| chr22 | 17224410 | 17224401 | 17224418 | SR_L=7,SR_R=16,DPE_L=31,DPE_R=50,MEI_MAPPED=82,polyA_MAPPED=27,polyA_side=R | SR_L=0,SR_R=0,DPE_L=0,DPE_R=0,MEI_MAPPED=0,polyA_MAPPED=0 | control_only | AACAAGTGCTAATAATTT | 68 | ALU | AluYb8#SINE/Alu | g1k:nssv14074719;lr:chr22-17900865-INS->s898731>s907592>s898732-334 | melt_1kg,long_read_1kg_ont_vienna | - | unnested | 291 | 1 | 291 |
| chr22 | 17567662 | 17567655 | 17567669 | SR_L=1,SR_R=13,DPE_L=33,DPE_R=44,MEI_MAPPED=72,polyA_MAPPED=7,polyA_side=R | SR_L=2,SR_R=25,DPE_L=57,DPE_R=80,MEI_MAPPED=122,polyA_MAPPED=9,polyA_side=R | shared | TATCCTTGCTTTTAT | 61 | ALU | AluYb8#SINE/Alu | chr22-18235412-INS->s898803>s907604>s907605>s907606>s898804-358 | long_read_1kg_ont_vienna | - | unnested | 318 | 1 | 318 |
| chr22 | 28371561 | 28371510 | 28371612 | SR_L=1,SR_R=9,DPE_L=16,DPE_R=10,MEI_MAPPED=29,polyA_MAPPED=60,polyA_side=R | SR_L=2,SR_R=0,DPE_L=32,DPE_R=25,MEI_MAPPED=44,polyA_MAPPED=116,polyA_side=L | shared |  | 34 | ALU | AluYi6#SINE/Alu |  |  | + | nested | 110 | 201 | 310 |
| chr22 | 29236892 | 29236685 | 29237100 | SR_L=0,SR_R=8,DPE_L=37,DPE_R=18,MEI_MAPPED=44,polyA_MAPPED=18,polyA_side=R | SR_L=0,SR_R=22,DPE_L=50,DPE_R=35,MEI_MAPPED=82,polyA_MAPPED=22,polyA_side=R | shared |  | 71 | LINE1 | L1MCa_5end#LINE/L1 |  |  | - | unnested | 1483 | 96 | 1578 |
| chr22 | 41639674 | 41639668 | 41639680 | SR_L=6,SR_R=5,DPE_L=15,DPE_R=2,MEI_MAPPED=35,polyA_MAPPED=44,polyA_side=R | SR_L=0,SR_R=1,DPE_L=31,DPE_R=3,MEI_MAPPED=54,polyA_MAPPED=77,polyA_side=R | shared | TTTTTTTTTTTTT | 42 | ALU | FLAM_C#SINE/Alu |  |  | - | nested | 0 | -1 | -1 |


## Examples

See [`docs/EXAMPLES.md`](docs/EXAMPLES.md) for additional annotated IGV review snapshots and read-architecture plots.

### Illumina

![Illumina chr22 retrotransposon insertion example](docs/examples/retrotransposon.gif)

The gif shows screenshots from random sections of chromosome 22 in a healthy individual. Grey bars represent unmutated DNA, and colors indicate either a mutation or errors in sequencing. The final screenshot shows a barcode-like signature indicating a retrotransposon insertion at one location. This insertion was not previously reported in this individual in published studies using the same data.

### SVA insertion (GRCh38 chr22:49029650)

<img src="docs/examples/grch38_sva_read_arch_chr22_49029238_49029720.png" alt="GRCh38 chr22 SVA read architecture" width="1470" />

Top-ranked shared SVA (`SVA_D`, `nssv14064350`) with two-sided SR/DPE support, TSD `AAGAAAACTCCT`, and polyA/VNTR rescue counts in the support string. Full-consensus span is ~1–1386.

### Alu insertion (GRCh38 chr22:31355872)

<img src="docs/examples/grch38_alu_read_arch_chr22_31355856_31355889.png" alt="GRCh38 chr22 Alu read architecture" width="1470" />

Top-ranked shared Alu (`AluYh7`) with a long TSD (`GCCCGCCTCGGCTTCCCAAAGTGCTGGGATTACA`) and near-full consensus coverage (~1–299).

### LINE-1 insertion (GRCh38 chr22:22131981)

<img src="docs/examples/grch38_line1_read_arch_chr22_22131552_22132407.png" alt="GRCh38 chr22 LINE-1 read architecture" width="1470" />

Known control-only LINE-1 (`L1HS`, `nssv14066334`) with split-read and discordant paired-end support. Panel `L1HS_5end` / `L1HS_3end` hits are projected onto the shared full-length L1 axis (near-full ~3–6018), not min/max’d on short fragment references.

## Getting Started on Amazon EC2 (Elastic Compute Cloud)

For whole-genome runs, use at least `r6i.4xlarge`.

The EC2 helper script (`scripts/ec2_jlab.sh`) works with **any existing EC2 instance** in your AWS account. Instance IDs and names are **not hardcoded in the repository**; each user binds their own instance locally to `.ec2-instance.env` (gitignored).

Run `./scripts/ec2_jlab.sh help` for the full command list.

### Quick Start (Bring Your Own EC2)

From your local machine:

```bash
git clone https://github.com/<org>/retrotransposon-miner.git
cd retrotransposon-miner
chmod +x scripts/ec2_jlab.sh
./scripts/ec2_jlab.sh list-instances
./scripts/ec2_jlab.sh use <instance-id-or-name>
./scripts/ec2_jlab.sh up
ssh retro-ec2
```

One-shot bind, start, and SSH config:

```bash
./scripts/ec2_jlab.sh up <instance-id-or-name>
./scripts/ec2_jlab.sh connect
```

After the instance is running:
- Secure Shell (SSH): `ssh retro-ec2`
- JupyterLab (after `start-jlab` + `start-tunnel`): `http://127.0.0.1:8890/lab?token=<printed-token>`

### Create a New EC2 Instance

Use `bootstrap` only when you want the script to provision a new instance (key pair, security group, Elastic IP, JupyterLab):

```bash
./scripts/ec2_jlab.sh bootstrap
```

`start-instance`, `stop-instance`, and `reboot-instance` operate on the **bound** instance only and do not create new instances.

### EC2 CLI Reference

| Command | Description |
|---|---|
| `list-instances` | List EC2 instances in the configured region |
| `use <instance-id-or-name>` | Bind an instance (by ID or `Name` tag) and update SSH config |
| `up [instance-id-or-name]` | Bind (optional), start, refresh SSH config |
| `connect [instance-id-or-name]` | `up` + SSH |
| `status` | Show bound instance state |
| `start-instance` | Start the bound instance |
| `stop-instance` | Stop the bound instance |
| `reboot-instance` | Reboot the bound instance |
| `bootstrap` | Create and configure a new EC2 instance |
| `start-jlab` / `stop-jlab` / `start-tunnel` | JupyterLab lifecycle |
| `help` | Show usage |

Optional environment variables:

- `REGION` — AWS region (default: from `aws configure`, else `us-east-1`)
- `INSTANCE_ID` / `INSTANCE_NAME` — override bound instance without editing `.ec2-instance.env`
- `HOST_ALIAS` — SSH config alias (default: `retro-ec2`)
- `SSH_USER` — SSH login user (auto-detected from AMI if unset; e.g. `ec2-user`, `ubuntu`)
- `KEY_PATH` — path to PEM for the instance key pair
- `INSTANCE_TYPE` — instance type for `bootstrap` only (default: `r6i.4xlarge`)

Binding is saved to `.ec2-instance.env` in the repo checkout. Rebind anytime with `use`.

### Prerequisites

Local tools:
- `aws` command-line interface (CLI) v2 (`aws configure` complete)
- `ssh`
- `curl`
- `git`

Identity and Access Management (IAM) permissions:

For bring-your-own-EC2 (`use`, `up`, `connect`, lifecycle commands):
- `ec2:Describe*`
- `ec2:StartInstances`
- `ec2:StopInstances`
- `ec2:RebootInstances`
- `ec2:AuthorizeSecurityGroupIngress` (SSH proxy refreshes your current IP on connect)

Additional permissions for `bootstrap` (new instance provisioning):
- `ec2:RunInstances`
- `ec2:CreateTags`
- `ec2:CreateKeyPair`
- `ec2:CreateSecurityGroup`
- `ec2:AllocateAddress`
- `ec2:AssociateAddress`
- `ec2:DescribeAddresses`
- `ec2:DescribeVpcs`
- `ec2:DescribeSubnets`
- `ssm:GetParameter` (Amazon Linux AMI lookup)
- `iam:PassRole` (if attaching an instance profile)

### What `scripts/ec2_jlab.sh` Does

- Binds to any existing EC2 instance by ID or `Name` tag (`use`).
- Saves the binding locally in `.ec2-instance.env` (not committed to git).
- Starts, stops, and reboots the bound instance without creating new ones.
- Writes SSH aliases (`retro-ec2`, `jlab`) into local `~/.ssh/config`.
- Refreshes SSH security group ingress for your current public IP on connect.
- Optionally creates a new instance (`bootstrap`), key pair, security group, and Elastic IP.
- Starts JupyterLab remotely and tunnels it locally.

### Quickstart Runs (chr22)

Use the main workflow wrapper:

- `scripts/run_candidate_discovery_and_annotation.sh`

Important: these quickstart commands do not download reference/public inputs automatically.

Step 0: download public/reference data first using the provided script:

```bash
conda activate rtm-miner || micromamba activate rtm-miner
python3 scripts/download_public_data.py \
  --references hg38 \
  --categories test_bam \
  --outdir "${RTM_PUBLIC_DATA_DIR:-$HOME/retrotransposon-workdir/data/public}"
```

SEQC2 chr22 test BAMs are sliced to chr22 **plus interchrom discordant mates** (needed for
discordant-mate MEI consensus remapping). Re-download with `--force` after updating the downloader.

If you plan to run both GRCh38 and hs1 workflows:

```bash
conda activate rtm-miner || micromamba activate rtm-miner
python3 scripts/download_public_data.py \
  --references hg38 hs1 \
  --outdir "${RTM_PUBLIC_DATA_DIR:-$HOME/retrotransposon-workdir/data/public}"
```

Tumor/normal chr22 quickstart (SEQC2 public test pair):

Local assembly is **off by default** (faster; sufficient for `bwa mem` mate/clip MEI
remapping and `supporting_reads_detail.mei.tsv`). Pass `--local-assembly` when you need
`asm_*` breakpoint/TSD fields from per-locus SPAdes.

Empirical gold gating (`--empirical-stage`) is also **off by default** (expensive BAM
depth/MAPQ/NM null sampling with little callset impact). Pass `--empirical-stage` only
when you want that extra filter.

Gold also requires `MEI_MAPPED>=3` in disease **or** control (silver loci with only 1–2
MEI-mapped reads stay silver). This keeps review/IGV/read-architecture plots focused on
better-supported calls.

```bash
RUN_IN_ENV=1 bash scripts/run_candidate_discovery_and_annotation.sh \
  --reference-build hg38 \
  --disease-bam "${RTM_PUBLIC_DATA_DIR:-$HOME/retrotransposon-workdir/data/public}/test_data/seqc2/chr22/disease.chr22.hg38.bam" \
  --control-bam "${RTM_PUBLIC_DATA_DIR:-$HOME/retrotransposon-workdir/data/public}/test_data/seqc2/chr22/control.chr22.hg38.bam" \
  --disease-mate-bam "${RTM_PUBLIC_DATA_DIR:-$HOME/retrotransposon-workdir/data/public}/test_data/seqc2/chr22/disease.chr22.hg38.bam" \
  --control-mate-bam "${RTM_PUBLIC_DATA_DIR:-$HOME/retrotransposon-workdir/data/public}/test_data/seqc2/chr22/control.chr22.hg38.bam" \
  --mei-fasta "${RTM_PUBLIC_DATA_DIR:-$HOME/retrotransposon-workdir/data/public}/retrotransposon_db/dfam/dfam_human_mei_l1_alu_sva.fasta" \
  --chr chr22 \
  --outdir "${RTM_RESULTS_DIR:-$HOME/retrotransposon-workdir/results}/quickstart_seqc2_chr22"
```

For SEQC2 chr22 BAMs built with `include_discordant_mates`, point `--disease-mate-bam` /
`--control-mate-bam` at the same files as `--disease-bam` / `--control-bam`.

Re-run annotation only (reuse existing `split_evidence.*` / `candidate_loci.tsv` after code or
BAM changes that affect MEI consensus remapping). Local assembly and empirical stage stay off
unless you pass `--local-assembly` / `--empirical-stage`:

```bash
RUN_IN_ENV=1 bash scripts/run_candidate_discovery_and_annotation.sh \
  --reference-build hg38 \
  --annotate-only \
  --disease-bam "${RTM_PUBLIC_DATA_DIR:-$HOME/retrotransposon-workdir/data/public}/test_data/seqc2/chr22/disease.chr22.hg38.bam" \
  --control-bam "${RTM_PUBLIC_DATA_DIR:-$HOME/retrotransposon-workdir/data/public}/test_data/seqc2/chr22/control.chr22.hg38.bam" \
  --disease-mate-bam "${RTM_PUBLIC_DATA_DIR:-$HOME/retrotransposon-workdir/data/public}/test_data/seqc2/chr22/disease.chr22.hg38.bam" \
  --control-mate-bam "${RTM_PUBLIC_DATA_DIR:-$HOME/retrotransposon-workdir/data/public}/test_data/seqc2/chr22/control.chr22.hg38.bam" \
  --mei-fasta "${RTM_PUBLIC_DATA_DIR:-$HOME/retrotransposon-workdir/data/public}/retrotransposon_db/dfam/dfam_human_mei_l1_alu_sva.fasta" \
  --chr chr22 \
  --outdir "${RTM_RESULTS_DIR:-$HOME/retrotransposon-workdir/results}/mei_step1_hg38_chr22"
```

Pipeline outputs include `supporting_reads_detail.mei.tsv` (per-read anchor/mate MEI coords for
architecture plots). Annotate writes gold-only plots under
`<out>.read_architecture/` by default (`--no-read-architecture-plots` to skip).

Single locus or batch from an existing gold review table:

```bash
# one locus (tables loaded once via cache)
python scripts/plot_locus_read_architecture.py \
  --gold-review-tsv "${RTM_RESULTS_DIR:-$HOME/retrotransposon-workdir/results}/quickstart_seqc2_chr22/candidate_loci.mei.gold_review.tsv" \
  --chrom chr22 --pos 49029650 --sample disease \
  --out-png "${RTM_RESULTS_DIR:-$HOME/retrotransposon-workdir/results}/quickstart_seqc2_chr22/plots/read_arch_chr22_49029650.png"

# all gold loci (same load-once path)
python scripts/plot_locus_read_architecture.py \
  --gold-review-tsv "${RTM_RESULTS_DIR:-$HOME/retrotransposon-workdir/results}/quickstart_seqc2_chr22/candidate_loci.mei.gold_review.tsv" \
  --all-gold \
  --out-dir "${RTM_RESULTS_DIR:-$HOME/retrotransposon-workdir/results}/quickstart_seqc2_chr22/read_architecture"
```

HG0001-style germline/control chr22 quickstart (replace with your BAM path):

```bash
bash scripts/run_candidate_discovery_and_annotation.sh \
  --reference-build hg38 \
  --disease-bam "/path/to/HG0001.chr22.hg38.bam" \
  --control-bam "/path/to/HG0001.chr22.hg38.bam" \
  --mei-fasta "${RTM_PUBLIC_DATA_DIR:-$HOME/retrotransposon-workdir/data/public}/retrotransposon_db/dfam/dfam_human_mei_l1_alu_sva.fasta" \
  --chr chr22 \
  --outdir "${RTM_RESULTS_DIR:-$HOME/retrotransposon-workdir/results}/quickstart_hg0001_chr22"
```

### Sync Repository on Virtual Machine (VM)

```bash
cd ~
git clone https://github.com/<org>/retrotransposon-miner.git
cd retrotransposon-miner
```

### Install Environment on Virtual Machine (VM)

```bash
bash scripts/bootstrap_env.sh
bash scripts/install_ucsc_tools.sh
conda activate rtm-miner || micromamba activate rtm-miner
bash scripts/validate_environment.sh
```

If needed:

```bash
eval "$($HOME/.local/bin/micromamba shell hook -s bash)"
micromamba activate rtm-miner
```

### Download Public Data

GRCh38:

```bash
conda activate rtm-miner || micromamba activate rtm-miner
python3 scripts/download_public_data.py \
  --references hg38 \
  --outdir "${RTM_PUBLIC_DATA_DIR:-$HOME/retrotransposon-workdir/data/public}"
```

hs1:

```bash
conda activate rtm-miner || micromamba activate rtm-miner
python3 scripts/download_public_data.py \
  --references hs1 \
  --outdir "${RTM_PUBLIC_DATA_DIR:-$HOME/retrotransposon-workdir/data/public}"
```

Both:

```bash
conda activate rtm-miner || micromamba activate rtm-miner
python3 scripts/download_public_data.py \
  --references hg38 hs1 \
  --outdir "${RTM_PUBLIC_DATA_DIR:-$HOME/retrotransposon-workdir/data/public}"
```

### Connect From Cursor

1. Open command palette.
2. Run `Remote-SSH: Connect to Host...`
3. Select `retro-ec2`.
4. Open `~/retrotransposon-miner`.

### Notes

- Designed for headless Linux execution with optional Integrative Genomics Viewer (IGV) snapshot generation.
- Instance bindings are local (`.ec2-instance.env`); nothing instance-specific is committed to git.
- If your public IP changes, SSH connect refreshes security group ingress automatically via `ec2_jlab.sh`.
- Set `KEY_PATH` if your PEM is not in `~/.ssh/` or the default search paths.
- For production use, review security hardening, key lifecycle, and cost controls.

## License

This project is licensed under the Apache License 2.0.

- Full text: [`LICENSE`](LICENSE)
- SPDX identifier: `Apache-2.0`

## Contributing

Contributions are welcome and encouraged.

- Contribution guide: [`CONTRIBUTING.md`](CONTRIBUTING.md)
- Community standards: [`CODE_OF_CONDUCT.md`](CODE_OF_CONDUCT.md)
- Contact: open a GitHub issue/discussion first, or email `william [at] l1tx [dot] com`.

If you submit code, please include clear validation steps and update documentation when behavior changes.

## Repository Metadata

- Repository description: `Retrotransposon mobile element insertion (MEI) caller for short-read whole-genome sequencing (WGS) with split-read + discordant-pair evidence, candidate ranking, MEI annotation, and Integrative Genomics Viewer (IGV) snapshot review workflows.`
- GitHub topics: `retrotransposon`, `mobile-element-insertion`, `mei`, `line1`, `alu`, `sva`, `genomics`, `bioinformatics`, `structural-variation`, `nextflow`, `igv`, `jupyterlab`, `aws`, `ec2`.
- Search keywords: `retrotransposon detection`, `mobile element insertion calling`, `LINE-1 insertion`, `Alu insertion`, `SVA insertion`, `short-read MEI pipeline`, `tumor normal MEI`, `germline MEI`, `IGV MEI review`.
