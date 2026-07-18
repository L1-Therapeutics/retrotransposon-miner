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

## Example Variant Calls (GRCh38)

The table below lists candidate insertion calls from a tumor/normal chr22 run.  
In the example table, `SR` denotes split-read evidence and `DPE` denotes discordant paired-end evidence.

### Example output from sample tumor/normal data

Gold-tier calls from the SEQC2 tumor/normal chr22 annotate run with polyA-trimmed MEI consensus remap, polyA-only (not MEI_MAPPED) junction clips, and polyA/T TSD filtering (top 30 of n=1042 by review rank; mix 22 Alu / 7 LINE1 / 1 SVA).

| chrom | consensus_insertion_breakpoint_pos | window_start | window_end | control_supporting_reads | disease_supporting_reads | sample_status_label | consensus_tsd_seq | consensus_poly_at_min_bp | consensus_mei_family | consensus_mei_subfamily | known_mei_polymorphism_id | known_mei_polymorphism_source | consensus_insertion_orientation | nested_in_same_MEI | consensus_insertion_mei_span_full | consensus_insertion_mei_5p_coord_full | consensus_insertion_mei_3p_coord_full |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| chr22 | 49029650 | 49029645 | 49029656 | SR_L=0,SR_R=0,DPE_L=2,DPE_R=152,MEI_MAPPED=127,polyA_MAPPED=46,VNTR_MAPPED=0,polyA_side=L | SR_L=0,SR_R=0,DPE_L=8,DPE_R=373,MEI_MAPPED=306,polyA_MAPPED=92,VNTR_MAPPED=1,polyA_side=L | shared | AAGAAAACTCCT | 50 | SVA | SVA_D#Retroposon/SVA | nssv14064350 | melt_1kg | - | unnested | 1366 | 1 | 1366 |
| chr22 | 31355872 | 31355858 | 31355887 | SR_L=12,SR_R=0,DPE_L=7,DPE_R=148,MEI_MAPPED=108,polyA_MAPPED=6 | SR_L=15,SR_R=0,DPE_L=19,DPE_R=273,MEI_MAPPED=188,polyA_MAPPED=16,polyA_side=R | shared | CCGCCTCGGCTTCCCAAAGTGCTGGGATTA | 71 | ALU | AluY_short_#SINE/Alu |  |  | - | nested | 281 | 1 | 281 |
| chr22 | 37529127 | 37529108 | 37529146 | SR_L=0,SR_R=1,DPE_L=151,DPE_R=2,MEI_MAPPED=66,polyA_MAPPED=47,polyA_side=L | SR_L=0,SR_R=3,DPE_L=365,DPE_R=16,MEI_MAPPED=167,polyA_MAPPED=86,polyA_side=L | shared | GAAGCGGAGGTTGCAGTGAGCCGAGATTGCGCCACTGCA | 90 | ALU | AluYb8#SINE/Alu |  |  | + | nested | 288 | 1 | 288 |
| chr22 | 17567662 | 17567655 | 17567669 | SR_L=6,SR_R=6,DPE_L=33,DPE_R=44,MEI_MAPPED=71,polyA_MAPPED=7,polyA_side=R | SR_L=23,SR_R=17,DPE_L=57,DPE_R=80,MEI_MAPPED=139,polyA_MAPPED=9,polyA_side=R | shared | TATCCTTGCTTTTAT | 61 | ALU | AluYb8#SINE/Alu | chr22-18235412-INS->s898803>s907604>s907605>s907606>s898804-358 | long_read_1kg_ont_vienna | - | unnested | 288 | 1 | 288 |
| chr22 | 50495066 | 50494596 | 50495537 | SR_L=2,SR_R=0,DPE_L=59,DPE_R=82,MEI_MAPPED=71,polyA_MAPPED=34,polyA_side=R | SR_L=8,SR_R=1,DPE_L=102,DPE_R=174,MEI_MAPPED=138,polyA_MAPPED=49,polyA_side=R | shared |  | 14 | ALU | AluYa5#SINE/Alu |  |  | - | unnested | 280 | 2 | 281 |
| chr22 | 45595784 | 45595639 | 45595930 | SR_L=0,SR_R=0,DPE_L=384,DPE_R=2,MEI_MAPPED=119,polyA_MAPPED=106,polyA_side=L | SR_L=0,SR_R=0,DPE_L=452,DPE_R=0,MEI_MAPPED=128,polyA_MAPPED=102,polyA_side=L | shared |  | 91 | ALU | AluJb_short_#SINE/Alu |  |  | - | unnested | 282 | 1 | 282 |
| chr22 | 41050312 | 41050276 | 41050348 | SR_L=3,SR_R=0,DPE_L=1,DPE_R=85,MEI_MAPPED=65,polyA_MAPPED=27,polyA_side=R | SR_L=1,SR_R=0,DPE_L=4,DPE_R=166,MEI_MAPPED=125,polyA_MAPPED=44,polyA_side=R | shared |  | 34 | ALU | AluSq#SINE/Alu |  |  | - | nested | 283 | 1 | 283 |
| chr22 | 20075438 | 20075432 | 20075444 | SR_L=12,SR_R=8,DPE_L=128,DPE_R=19,MEI_MAPPED=111,polyA_MAPPED=51,polyA_side=L | SR_L=7,SR_R=3,DPE_L=89,DPE_R=16,MEI_MAPPED=65,polyA_MAPPED=40,polyA_side=R | shared | AGATTTCTTTTCT | 39 | ALU | AluYk12#SINE/Alu |  |  | + | unnested | 281 | 1 | 281 |
| chr22 | 40007330 | 40007328 | 40007332 | SR_L=0,SR_R=0,DPE_L=76,DPE_R=12,MEI_MAPPED=68,polyA_MAPPED=1,polyA_side=R | SR_L=4,SR_R=0,DPE_L=98,DPE_R=25,MEI_MAPPED=98,polyA_MAPPED=3,polyA_side=R | shared | CTCCT | 114 | ALU | AluYb8#SINE/Alu |  |  | - | unnested | 288 | 1 | 288 |
| chr22 | 19223382 | 19223373 | 19223390 | SR_L=0,SR_R=14,DPE_L=60,DPE_R=27,MEI_MAPPED=91,polyA_MAPPED=17,polyA_side=L | SR_L=0,SR_R=0,DPE_L=0,DPE_R=0,MEI_MAPPED=0,polyA_MAPPED=0 | control_only | AAAAACCACCTATGCTGG | 66 | LINE1 | L1HS_5end#LINE/L1 | g1k:nssv14064681;lr:chr22-19600083-INS->s899391<s914453>s899392-6059 | melt_1kg,long_read_1kg_ont_vienna | + | unnested | 6018 | 1 | 6018 |
| chr22 | 34034616 | 34034610 | 34034623 | SR_L=14,SR_R=2,DPE_L=62,DPE_R=1,MEI_MAPPED=61,polyA_MAPPED=14,polyA_side=R | SR_L=22,SR_R=6,DPE_L=92,DPE_R=2,MEI_MAPPED=90,polyA_MAPPED=32,polyA_side=R | shared | CAAATGGAACTTTT | 64 | ALU | AluYb8#SINE/Alu | nssv14071620 | melt_1kg | - | unnested | 288 | 1 | 288 |
| chr22 | 50351083 | 50350826 | 50351340 | SR_L=0,SR_R=0,DPE_L=121,DPE_R=18,MEI_MAPPED=50,polyA_MAPPED=12,polyA_side=R | SR_L=0,SR_R=0,DPE_L=218,DPE_R=16,MEI_MAPPED=85,polyA_MAPPED=38,polyA_side=R | shared |  | 38 | ALU | AluSc5#SINE/Alu |  |  | - | nested | 264 | 1 | 264 |
| chr22 | 33132520 | 33132513 | 33132527 | SR_L=0,SR_R=7,DPE_L=44,DPE_R=7,MEI_MAPPED=37,polyA_MAPPED=31,polyA_side=L | SR_L=1,SR_R=7,DPE_L=76,DPE_R=31,MEI_MAPPED=85,polyA_MAPPED=23,polyA_side=L | shared | AAAAGTCATTATTAG | 56 | ALU | AluYg6#SINE/Alu | nssv14075885 | melt_1kg | + | unnested | 281 | 1 | 281 |
| chr22 | 36746494 | 36746492 | 36746495 | SR_L=0,SR_R=2,DPE_L=7,DPE_R=77,MEI_MAPPED=44,polyA_MAPPED=2 | SR_L=1,SR_R=11,DPE_L=11,DPE_R=157,MEI_MAPPED=85,polyA_MAPPED=2,polyA_side=R | shared | CTCC | 48 | ALU | AluSz#SINE/Alu |  |  | - | unnested | 221 | 10 | 230 |
| chr22 | 29236892 | 29236685 | 29237100 | SR_L=0,SR_R=6,DPE_L=37,DPE_R=18,MEI_MAPPED=43,polyA_MAPPED=18,polyA_side=R | SR_L=0,SR_R=24,DPE_L=50,DPE_R=35,MEI_MAPPED=83,polyA_MAPPED=22,polyA_side=R | shared |  | 71 | LINE1 | L1MCa_5end#LINE/L1 |  |  | - | unnested | 1483 | 96 | 1578 |
| chr22 | 20521112 | 20521097 | 20521126 | SR_L=0,SR_R=0,DPE_L=39,DPE_R=19,MEI_MAPPED=50,polyA_MAPPED=19,polyA_side=R | SR_L=0,SR_R=1,DPE_L=55,DPE_R=33,MEI_MAPPED=78,polyA_MAPPED=20,polyA_side=R | shared |  | 28 | ALU | AluYa5#SINE/Alu |  |  | + | unnested | 281 | 1 | 281 |
| chr22 | 19919244 | 19919236 | 19919251 | SR_L=0,SR_R=5,DPE_L=21,DPE_R=117,MEI_MAPPED=78,polyA_MAPPED=16,polyA_side=L | SR_L=0,SR_R=5,DPE_L=17,DPE_R=97,MEI_MAPPED=63,polyA_MAPPED=3,polyA_side=L | shared | CCCAGGCTGGAGTGCA | 71 | ALU | AluSp#SINE/Alu | nssv14053291 | melt_1kg | - | nested | 273 | 1 | 273 |
| chr22 | 20595738 | 20595139 | 20596336 | SR_L=0,SR_R=9,DPE_L=3,DPE_R=43,MEI_MAPPED=30,polyA_MAPPED=0 | SR_L=0,SR_R=26,DPE_L=5,DPE_R=77,MEI_MAPPED=78,polyA_MAPPED=0 | shared |  | 75 | ALU | AluSx3#SINE/Alu |  |  | - | nested | 67 | 69 | 135 |
| chr22 | 41835230 | 41835229 | 41835232 | SR_L=27,SR_R=18,DPE_L=19,DPE_R=2,MEI_MAPPED=46,polyA_MAPPED=5,polyA_side=L | SR_L=50,SR_R=6,DPE_L=43,DPE_R=6,MEI_MAPPED=70,polyA_MAPPED=10,polyA_side=L | shared | ATAG | 84 | ALU | AluJo#SINE/Alu |  |  | - | nested | 238 | 44 | 281 |
| chr22 | 36752165 | 36751652 | 36752678 | SR_L=0,SR_R=0,DPE_L=66,DPE_R=8,MEI_MAPPED=29,polyA_MAPPED=3,polyA_side=R | SR_L=1,SR_R=0,DPE_L=145,DPE_R=18,MEI_MAPPED=68,polyA_MAPPED=6,polyA_side=R | shared |  | 14 | ALU | AluSz#SINE/Alu |  |  | - | unnested | 210 | 10 | 219 |
| chr22 | 42818644 | 42818163 | 42819124 | SR_L=0,SR_R=0,DPE_L=51,DPE_R=25,MEI_MAPPED=67,polyA_MAPPED=19,polyA_side=R | SR_L=0,SR_R=0,DPE_L=0,DPE_R=0,MEI_MAPPED=2,polyA_MAPPED=0 | control_only |  | 26 | ALU | AluYa5#SINE/Alu | chr22-43299733-INS->s903600<s909229>s903601-318 | long_read_1kg_ont_vienna | + | unnested | 281 | 1 | 281 |
| chr22 | 29239707 | 29239402 | 29240012 | SR_L=1,SR_R=0,DPE_L=44,DPE_R=4,MEI_MAPPED=34,polyA_MAPPED=1,polyA_side=L | SR_L=4,SR_R=1,DPE_L=72,DPE_R=16,MEI_MAPPED=65,polyA_MAPPED=5,polyA_side=R | shared |  | 59 | LINE1 | L1MCa_5end#LINE/L1 |  |  | - | nested | 364 | 96 | 459 |
| chr22 | 17224410 | 17224401 | 17224418 | SR_L=7,SR_R=0,DPE_L=31,DPE_R=50,MEI_MAPPED=63,polyA_MAPPED=27,polyA_side=R | SR_L=0,SR_R=0,DPE_L=0,DPE_R=0,MEI_MAPPED=1,polyA_MAPPED=0 | control_only | AACAAGTGCTAATAATTT | 68 | ALU | AluYb8#SINE/Alu | g1k:nssv14074719;lr:chr22-17900865-INS->s898731>s907592>s898732-334 | melt_1kg,long_read_1kg_ont_vienna | - | unnested | 288 | 1 | 288 |
| chr22 | 42705164 | 42704870 | 42705458 | SR_L=0,SR_R=2,DPE_L=6,DPE_R=71,MEI_MAPPED=35,polyA_MAPPED=2 | SR_L=0,SR_R=5,DPE_L=21,DPE_R=121,MEI_MAPPED=63,polyA_MAPPED=6,polyA_side=L | shared |  | 48 | ALU | AluYm1#SINE/Alu |  |  | + | unnested | 280 | 2 | 281 |
| chr22 | 33124372 | 33124353 | 33124392 | SR_L=1,SR_R=14,DPE_L=9,DPE_R=19,MEI_MAPPED=32,polyA_MAPPED=1,polyA_side=R | SR_L=15,SR_R=17,DPE_L=12,DPE_R=47,MEI_MAPPED=63,polyA_MAPPED=0 | shared | GAAAGAAGGAAGGAAGGAAGGAAGGAAGGAAGGGAGGAAG | 16 | LINE1 | L1M1_5end#LINE/L1 |  |  | + | unnested | 1462 | 72 | 1533 |
| chr22 | 49760480 | 49760479 | 49760482 | SR_L=27,SR_R=10,DPE_L=3,DPE_R=7,MEI_MAPPED=39,polyA_MAPPED=0 | SR_L=27,SR_R=27,DPE_L=2,DPE_R=10,MEI_MAPPED=63,polyA_MAPPED=1,polyA_side=R | shared | CATG | 81 | LINE1 | L1M2a1_5end#LINE/L1 |  |  | + | nested | 189 | 149 | 337 |
| chr22 | 31380162 | 31379890 | 31380433 | SR_L=0,SR_R=9,DPE_L=10,DPE_R=30,MEI_MAPPED=26,polyA_MAPPED=1,polyA_side=L | SR_L=0,SR_R=22,DPE_L=17,DPE_R=63,MEI_MAPPED=61,polyA_MAPPED=0 | shared |  | 66 | LINE1 | L1PREC2_orf2#LINE/L1 |  |  | - | unnested | 2154 | 2070 | 4223 |
| chr22 | 23938127 | 23937906 | 23938348 | SR_L=0,SR_R=0,DPE_L=44,DPE_R=20,MEI_MAPPED=21,polyA_MAPPED=33,polyA_side=R | SR_L=0,SR_R=0,DPE_L=79,DPE_R=50,MEI_MAPPED=60,polyA_MAPPED=60,polyA_side=R | shared |  | 70 | ALU | AluYc#SINE/Alu |  |  | + | nested | 262 | 8 | 269 |
| chr22 | 41051286 | 41050603 | 41051968 | SR_L=0,SR_R=0,DPE_L=67,DPE_R=0,MEI_MAPPED=36,polyA_MAPPED=19,polyA_side=R | SR_L=0,SR_R=0,DPE_L=96,DPE_R=4,MEI_MAPPED=59,polyA_MAPPED=30,polyA_side=R | shared |  | 31 | ALU | AluSp#SINE/Alu |  |  | - | unnested | 283 | 1 | 283 |
| chr22 | 17289460 | 17289460 | 17289460 | SR_L=0,SR_R=11,DPE_L=0,DPE_R=38,MEI_MAPPED=45,polyA_MAPPED=0 | SR_L=0,SR_R=13,DPE_L=3,DPE_R=60,MEI_MAPPED=59,polyA_MAPPED=0 | shared |  | 74 | LINE1 | L1P2_5end#LINE/L1 |  |  | - | unnested | 499 | 1056 | 1554 |

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

## MEI consensus remapping

Soft-clips and discordant clipped ends remap to the Dfam Alu/LINE-1/SVA panel with `bwa mem -k10 -T10` (`--bwa-threads`; wrapper auto: `nproc` single-chrom, `1` under multi-chrom concurrency). Queries are polyA/T-trimmed (≥8 bp) before align; consensus targets are also terminal-polyA-trimmed (sidecar `*.nopolya.fa`, and prep writes body-only Dfam/panel FASTAs) so clips cannot map onto the A-tail. Junction clips that are themselves polyA/T count as `polyA_MAPPED` only — never also `MEI_MAPPED`/SR. Short tips (≤30 bp) need qcov≥0.80 and pid≥0.90, longer clips need pid≥0.90 and alnlen≥20. Panel fragment hits project onto one family-consistent `*_full` axis via `mei_fragment_to_full_coords.tsv` (prep: `bwa mem -a` on trimmed sequences). Assembly contig-to-MEI still uses minimap2. Benchmark: `scripts/benchmark_mei_aligners.py`.

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
