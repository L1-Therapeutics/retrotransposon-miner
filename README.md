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

- Designed primarily for Amazon Web Services (AWS) machines today; relatively straightforward to adapt to Google Cloud Platform (GCP), Azure, or local Linux. A single full genome wants about 64 vCPU and 256 GiB of memory (see the EC2 section). Local assembly and candidate processing support parallel execution.
- Artificial intelligence/machine learning (AI/ML) genotyping confidence models are still under active development.
- Genotyping is not supported: VCF output leaves `GT`/`GQ` blank. Calls are pooled disease-vs-control read support, not per-individual diploid genotypes, and a simple alt/ref read-ratio model is not valid here — split/discordant "alt" evidence and proper-pair "ref" evidence are structurally different read populations, and ploidy varies with chromosome (chrX), somatic copy-number context, and mosaicism. Reliable MEI genotyping would need haplotype-resolved/pangenome references or long reads.
- Reverse-transcribed pseudogene insertion support is not yet added.
- Support for species other than *Homo sapiens* (for example, *Mus musculus*) is not yet implemented.
- Long-read native calling is not yet supported.
- Single-cell sequencing data is not yet supported.
- Other short-read platforms (for example, Ultima Genomics) are not yet validated/supported.
- Local assembly is parallelized but still compute-expensive; it is optional and not recommended by default for routine runs.

## Example Variant Calls (GRCh38)

HG03086 chromosome 22, classifier `gold_score` >= 0.997 (21 calls).

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG03086
chr22	17567662	L1TX-chr22-17567662-ALU	N	<INS:ME:ALU>	.	PASS	SVTYPE=INS;END=17567662;WINDOWSTART=17567655;WINDOWEND=17567669;MEIFAMILY=ALU;MEISUBFAMILY=AluYb8#SINE/Alu;TSD=TATCCTTGCTTTTAT;POLYA_MIN_BP=119;ORIENT=-;NESTED=unnested;MEI_SPAN=288;MEI_5P=1;MEI_3P=288;LR=chr22-18235412-INS->s898803>s907604>s907605>s907606>s898804-358;KNOWN_SRC=long_read_1kg_ont_vienna;SAMPLE_STATUS=shared;SUPPORT=SR_L:4|SR_R:3|DPE_L:35|DPE_R:9|MEI_MAPPED:47|polyA_MAPPED:17|polyA_side:R;L1TXGOLDSCORE=1.000;L1TXRANKPCT=98;INSERTIONSCORE=0.378;PEAKDEPTHZ=0.0804;SPLITCLUSTERZ=4.03;CALLTIER=none;KNOWNMEI=True	GT:GQ	./.:.
chr22	17577694	L1TX-chr22-17577694-ALU	N	<INS:ME:ALU>	.	PASS	SVTYPE=INS;END=17577694;WINDOWSTART=17577685;WINDOWEND=17577704;MEIFAMILY=ALU;MEISUBFAMILY=AluSc5#SINE/Alu;TSD=ATTCTCCTGCCTCAGCCTCC;POLYA_MIN_BP=94;ORIENT=-;NESTED=nested;MEI_SPAN=266;MEI_5P=10;MEI_3P=275;SAMPLE_STATUS=shared;SUPPORT=SR_L:0|SR_R:0|DPE_L:39|DPE_R:23|MEI_MAPPED:37|polyA_MAPPED:10|polyA_side:R;L1TXGOLDSCORE=0.9970;L1TXRANKPCT=90;INSERTIONSCORE=0.379;PEAKDEPTHZ=-0.00109;SPLITCLUSTERZ=2.03;CALLTIER=none;KNOWNMEI=False	GT:GQ	./.:.
chr22	17653864	L1TX-chr22-17653864-SVA	N	<INS:ME:SVA>	.	PASS	SVTYPE=INS;END=17653864;WINDOWSTART=17653857;WINDOWEND=17653870;MEIFAMILY=SVA;MEISUBFAMILY=SVA_F#Retroposon/SVA;TSD=AAAAATTGTTTATC;POLYA_MIN_BP=71;ORIENT=+;NESTED=unnested;MEI_SPAN=693;MEI_5P=670;MEI_3P=1362;G1K=nssv14065494;LR=chr22-18321622-INS->s898844<s916635<s916634<s916633<s916632<s916631>s898845-1712;KNOWN_SRC=melt_1kg|long_read_1kg_ont_vienna;SAMPLE_STATUS=shared;SUPPORT=SR_L:0|SR_R:6|DPE_L:11|DPE_R:18|MEI_MAPPED:35|polyA_MAPPED:12|VNTR_MAPPED:0|polyA_side:L;L1TXGOLDSCORE=0.9994;L1TXRANKPCT=91;INSERTIONSCORE=0.584;PEAKDEPTHZ=-0.0337;SPLITCLUSTERZ=5.3;CALLTIER=none;KNOWNMEI=True	GT:GQ	./.:.
chr22	19919244	L1TX-chr22-19919244-ALU	N	<INS:ME:ALU>	.	PASS	SVTYPE=INS;END=19919244;WINDOWSTART=19919236;WINDOWEND=19919251;MEIFAMILY=ALU;MEISUBFAMILY=AluSp#SINE/Alu;TSD=CCCAGGCTGGAGTGCA;POLYA_MIN_BP=69;ORIENT=-;NESTED=nested;MEI_SPAN=282;MEI_5P=1;MEI_3P=282;G1K=nssv14053291;KNOWN_SRC=melt_1kg;SAMPLE_STATUS=shared;SUPPORT=SR_L:0|SR_R:6|DPE_L:25|DPE_R:44|MEI_MAPPED:68|polyA_MAPPED:8|polyA_side:L;L1TXGOLDSCORE=1.000;L1TXRANKPCT=97;INSERTIONSCORE=0.422;PEAKDEPTHZ=0.0397;SPLITCLUSTERZ=11.7;CALLTIER=none;KNOWNMEI=True	GT:GQ	./.:.
chr22	20075438	L1TX-chr22-20075438-ALU	N	<INS:ME:ALU>	.	PASS	SVTYPE=INS;END=20075438;WINDOWSTART=20075432;WINDOWEND=20075444;MEIFAMILY=ALU;MEISUBFAMILY=AluYk12#SINE/Alu;POLYA_MIN_BP=109;ORIENT=-;NESTED=unnested;MEI_SPAN=281;MEI_5P=1;MEI_3P=281;SAMPLE_STATUS=shared;SUPPORT=SR_L:4|SR_R:14|DPE_L:3|DPE_R:99|MEI_MAPPED:93|polyA_MAPPED:52|polyA_side:R;L1TXGOLDSCORE=0.9999;L1TXRANKPCT=94;INSERTIONSCORE=0.459;PEAKDEPTHZ=0.602;SPLITCLUSTERZ=10.8;CALLTIER=none;KNOWNMEI=False	GT:GQ	./.:.
chr22	22131981	L1TX-chr22-22131981-LINE1	N	<INS:ME:LINE1>	.	PASS	SVTYPE=INS;END=22131981;WINDOWSTART=22131976;WINDOWEND=22131986;MEIFAMILY=LINE1;MEISUBFAMILY=L1HS_5end#LINE/L1;TSD=GCATATTTCTT;POLYA_MIN_BP=82;ORIENT=-;NESTED=unnested;MEI_SPAN=6015;MEI_5P=4;MEI_3P=6018;G1K=nssv14066334;KNOWN_SRC=melt_1kg;SAMPLE_STATUS=shared;SUPPORT=SR_L:4|SR_R:0|DPE_L:15|DPE_R:19|MEI_MAPPED:38|polyA_MAPPED:10|polyA_side:R;L1TXGOLDSCORE=1.000;L1TXRANKPCT=99;INSERTIONSCORE=0.557;PEAKDEPTHZ=-0.0255;SPLITCLUSTERZ=6.35;CALLTIER=none;KNOWNMEI=True	GT:GQ	./.:.
chr22	23928268	L1TX-chr22-23928268-ALU	N	<INS:ME:ALU>	.	PASS	SVTYPE=INS;END=23928268;WINDOWSTART=23928262;WINDOWEND=23928275;MEIFAMILY=ALU;MEISUBFAMILY=AluYb8#SINE/Alu;TSD=AAGAGATGGACTGA;POLYA_MIN_BP=79;ORIENT=+;NESTED=unnested;MEI_SPAN=288;MEI_5P=1;MEI_3P=288;G1K=nssv14081199;LR=chr22-24372840-INS->s900600>s911750>s911751>s911752>s900601-338;KNOWN_SRC=melt_1kg|long_read_1kg_ont_vienna;SAMPLE_STATUS=shared;SUPPORT=SR_L:0|SR_R:10|DPE_L:13|DPE_R:19|MEI_MAPPED:38|polyA_MAPPED:7|polyA_side:L;L1TXGOLDSCORE=1.000;L1TXRANKPCT=98;INSERTIONSCORE=0.618;PEAKDEPTHZ=-0.00109;SPLITCLUSTERZ=9.37;CALLTIER=none;KNOWNMEI=True	GT:GQ	./.:.
chr22	26557759	L1TX-chr22-26557759-ALU	N	<INS:ME:ALU>	.	PASS	SVTYPE=INS;END=26557759;WINDOWSTART=26557752;WINDOWEND=26557766;MEIFAMILY=ALU;MEISUBFAMILY=AluSx4#SINE/Alu;TSD=AGAAGAGAGATGGGG;POLYA_MIN_BP=100;ORIENT=+;NESTED=unnested;MEI_SPAN=69;MEI_5P=2;MEI_3P=70;LR=chr22-27020130-INS->s900992>s907497>s900993-322;KNOWN_SRC=long_read_1kg_ont_vienna;SAMPLE_STATUS=shared;SUPPORT=SR_L:2|SR_R:8|DPE_L:11|DPE_R:23|MEI_MAPPED:43|polyA_MAPPED:8|polyA_side:L;L1TXGOLDSCORE=0.9999;L1TXRANKPCT=95;INSERTIONSCORE=0.447;PEAKDEPTHZ=0.0234;SPLITCLUSTERZ=8.53;CALLTIER=none;KNOWNMEI=True	GT:GQ	./.:.
chr22	27154836	L1TX-chr22-27154836-ALU	N	<INS:ME:ALU>	.	PASS	SVTYPE=INS;END=27154836;WINDOWSTART=27154829;WINDOWEND=27154842;MEIFAMILY=ALU;MEISUBFAMILY=AluYb8#SINE/Alu;TSD=AAGAATAGACACGT;POLYA_MIN_BP=82;ORIENT=+;NESTED=unnested;MEI_SPAN=288;MEI_5P=1;MEI_3P=288;G1K=nssv14073646;LR=chr22-27615738-INS->s901049<s911832<s911831<s911830>s901050-326;KNOWN_SRC=melt_1kg|long_read_1kg_ont_vienna;SAMPLE_STATUS=shared;SUPPORT=SR_L:0|SR_R:5|DPE_L:27|DPE_R:22|MEI_MAPPED:53|polyA_MAPPED:9|polyA_side:L;L1TXGOLDSCORE=1.000;L1TXRANKPCT=96;INSERTIONSCORE=0.514;PEAKDEPTHZ=-0.0255;SPLITCLUSTERZ=4.77;CALLTIER=none;KNOWNMEI=True	GT:GQ	./.:.
chr22	28042280	L1TX-chr22-28042280-ALU	N	<INS:ME:ALU>	.	PASS	SVTYPE=INS;END=28042280;WINDOWSTART=28042272;WINDOWEND=28042288;MEIFAMILY=ALU;MEISUBFAMILY=AluY_short_#SINE/Alu;TSD=CATGCACACGTATTTTT;POLYA_MIN_BP=67;ORIENT=-;NESTED=unnested;MEI_SPAN=281;MEI_5P=1;MEI_3P=281;G1K=nssv14075695;LR=chr22-28503683-INS->s901149>s907506>s901150-322;KNOWN_SRC=melt_1kg|long_read_1kg_ont_vienna;SAMPLE_STATUS=shared;SUPPORT=SR_L:6|SR_R:0|DPE_L:20|DPE_R:13|MEI_MAPPED:33|polyA_MAPPED:4;L1TXGOLDSCORE=0.9996;L1TXRANKPCT=92;INSERTIONSCORE=0.593;PEAKDEPTHZ=-0.0337;SPLITCLUSTERZ=2.63;CALLTIER=none;KNOWNMEI=True	GT:GQ	./.:.
chr22	28059921	L1TX-chr22-28059921-LINE1	N	<INS:ME:LINE1>	.	PASS	SVTYPE=INS;END=28059921;WINDOWSTART=28059914;WINDOWEND=28059928;MEIFAMILY=LINE1;MEISUBFAMILY=L1HS_5end#LINE/L1;TSD=TACCCATTTATTTTC;POLYA_MIN_BP=76;ORIENT=-;NESTED=nested;MEI_SPAN=6007;MEI_5P=12;MEI_3P=6018;LR=chr22-28521335-INS->s901150<s919771>s901151-6079;KNOWN_SRC=long_read_1kg_ont_vienna;SAMPLE_STATUS=shared;SUPPORT=SR_L:4|SR_R:0|DPE_L:11|DPE_R:28|MEI_MAPPED:43|polyA_MAPPED:15|polyA_side:R;L1TXGOLDSCORE=1.000;L1TXRANKPCT=98;INSERTIONSCORE=0.599;PEAKDEPTHZ=-0.0337;SPLITCLUSTERZ=8.68;CALLTIER=none;KNOWNMEI=True	GT:GQ	./.:.
chr22	31355872	L1TX-chr22-31355872-ALU	N	<INS:ME:ALU>	.	PASS	SVTYPE=INS;END=31355872;WINDOWSTART=31355858;WINDOWEND=31355887;MEIFAMILY=ALU;MEISUBFAMILY=AluSx3#SINE/Alu;TSD=CCGCCTCGGCTTCCCAAAGTGCTGGGATTA;POLYA_MIN_BP=82;ORIENT=+;NESTED=unnested;MEI_SPAN=46;MEI_5P=58;MEI_3P=103;SAMPLE_STATUS=shared;SUPPORT=SR_L:8|SR_R:1|DPE_L:30|DPE_R:31|MEI_MAPPED:63|polyA_MAPPED:0;L1TXGOLDSCORE=0.9977;L1TXRANKPCT=90;INSERTIONSCORE=0.399;PEAKDEPTHZ=0.0478;SPLITCLUSTERZ=11.6;CALLTIER=none;KNOWNMEI=False	GT:GQ	./.:.
chr22	35673763	L1TX-chr22-35673763-ALU	N	<INS:ME:ALU>	.	PASS	SVTYPE=INS;END=35673763;WINDOWSTART=35673757;WINDOWEND=35673769;MEIFAMILY=ALU;MEISUBFAMILY=AluYa5#SINE/Alu;TSD=AAAAAGGGGAGGC;POLYA_MIN_BP=49;ORIENT=+;NESTED=unnested;MEI_SPAN=281;MEI_5P=1;MEI_3P=281;G1K=nssv14074982;LR=chr22-36133492-INS->s902223>s916221>s916222>s916223>s902224-332;KNOWN_SRC=melt_1kg|long_read_1kg_ont_vienna;SAMPLE_STATUS=shared;SUPPORT=SR_L:1|SR_R:7|DPE_L:18|DPE_R:9|MEI_MAPPED:28|polyA_MAPPED:6|polyA_side:L;L1TXGOLDSCORE=0.9989;L1TXRANKPCT=91;INSERTIONSCORE=0.646;PEAKDEPTHZ=-0.0337;SPLITCLUSTERZ=7.55;CALLTIER=none;KNOWNMEI=True	GT:GQ	./.:.
chr22	37529132	L1TX-chr22-37529132-ALU	N	<INS:ME:ALU>	.	PASS	SVTYPE=INS;END=37529132;WINDOWSTART=37529118;WINDOWEND=37529146;MEIFAMILY=ALU;MEISUBFAMILY=AluYb8#SINE/Alu;TSD=TTGCAGTGAGCCGAGATTGCGCCACTGCA;POLYA_MIN_BP=111;ORIENT=+;NESTED=nested;MEI_SPAN=266;MEI_5P=1;MEI_3P=266;SAMPLE_STATUS=shared;SUPPORT=SR_L:2|SR_R:4|DPE_L:42|DPE_R:14|MEI_MAPPED:42|polyA_MAPPED:15|polyA_side:R;L1TXGOLDSCORE=0.9994;L1TXRANKPCT=92;INSERTIONSCORE=0.57;PEAKDEPTHZ=-0.0255;SPLITCLUSTERZ=13.8;CALLTIER=none;KNOWNMEI=False	GT:GQ	./.:.
chr22	41050312	L1TX-chr22-41050312-ALU	N	<INS:ME:ALU>	.	PASS	SVTYPE=INS;END=41050312;WINDOWSTART=41050276;WINDOWEND=41050348;MEIFAMILY=ALU;MEISUBFAMILY=AluSp#SINE/Alu;POLYA_MIN_BP=88;ORIENT=-;NESTED=nested;MEI_SPAN=278;MEI_5P=4;MEI_3P=281;SAMPLE_STATUS=shared;SUPPORT=SR_L:3|SR_R:1|DPE_L:29|DPE_R:37|MEI_MAPPED:62|polyA_MAPPED:21|polyA_side:R;L1TXGOLDSCORE=0.9997;L1TXRANKPCT=92;INSERTIONSCORE=0.411;PEAKDEPTHZ=0.121;SPLITCLUSTERZ=3.99;CALLTIER=none;KNOWNMEI=False	GT:GQ	./.:.
chr22	41300515	L1TX-chr22-41300515-ALU	N	<INS:ME:ALU>	.	PASS	SVTYPE=INS;END=41300515;WINDOWSTART=41300511;WINDOWEND=41300519;MEIFAMILY=ALU;MEISUBFAMILY=AluSq2#SINE/Alu;TSD=GCCCAGGCT;POLYA_MIN_BP=65;ORIENT=-;NESTED=nested;MEI_SPAN=282;MEI_5P=1;MEI_3P=282;SAMPLE_STATUS=shared;SUPPORT=SR_L:1|SR_R:0|DPE_L:2|DPE_R:57|MEI_MAPPED:35|polyA_MAPPED:6|polyA_side:L;L1TXGOLDSCORE=0.9993;L1TXRANKPCT=91;INSERTIONSCORE=0.43;PEAKDEPTHZ=0.0234;SPLITCLUSTERZ=8.14;CALLTIER=none;KNOWNMEI=False	GT:GQ	./.:.
chr22	41739555	L1TX-chr22-41739555-ALU	N	<INS:ME:ALU>	.	PASS	SVTYPE=INS;END=41739555;WINDOWSTART=41739549;WINDOWEND=41739561;MEIFAMILY=ALU;MEISUBFAMILY=AluYb8#SINE/Alu;TSD=AAAAATTGAAACA;POLYA_MIN_BP=100;ORIENT=+;NESTED=unnested;MEI_SPAN=288;MEI_5P=1;MEI_3P=288;G1K=nssv14070718;LR=chr22-42218501-INS->s903362>s913963>s913964>s913965>s903363-325;KNOWN_SRC=melt_1kg|long_read_1kg_ont_vienna;SAMPLE_STATUS=shared;SUPPORT=SR_L:1|SR_R:4|DPE_L:13|DPE_R:21|MEI_MAPPED:35|polyA_MAPPED:9|polyA_side:L;L1TXGOLDSCORE=0.9999;L1TXRANKPCT=96;INSERTIONSCORE=0.343;PEAKDEPTHZ=0.0641;SPLITCLUSTERZ=4.0;CALLTIER=none;KNOWNMEI=True	GT:GQ	./.:.
chr22	46105398	L1TX-chr22-46105398-ALU	N	<INS:ME:ALU>	.	PASS	SVTYPE=INS;END=46105398;WINDOWSTART=46105391;WINDOWEND=46105406;MEIFAMILY=ALU;MEISUBFAMILY=AluYd8#SINE/Alu;TSD=AGTGTGTGCTTTTTCT;POLYA_MIN_BP=87;ORIENT=-;NESTED=unnested;MEI_SPAN=269;MEI_5P=1;MEI_3P=269;G1K=nssv14071372;LR=chr22-46590080-INS->s904510<s908699>s904511-316;KNOWN_SRC=melt_1kg|long_read_1kg_ont_vienna;SAMPLE_STATUS=shared;SUPPORT=SR_L:15|SR_R:2|DPE_L:3|DPE_R:33|MEI_MAPPED:49|polyA_MAPPED:14|polyA_side:R;L1TXGOLDSCORE=1.000;L1TXRANKPCT=98;INSERTIONSCORE=0.549;PEAKDEPTHZ=0.154;SPLITCLUSTERZ=11.8;CALLTIER=none;KNOWNMEI=True	GT:GQ	./.:.
chr22	46843468	L1TX-chr22-46843468-ALU	N	<INS:ME:ALU>	.	PASS	SVTYPE=INS;END=46843468;WINDOWSTART=46843462;WINDOWEND=46843475;MEIFAMILY=ALU;MEISUBFAMILY=AluYb8#SINE/Alu;TSD=CATTCACTGTTATT;POLYA_MIN_BP=53;ORIENT=-;NESTED=unnested;MEI_SPAN=288;MEI_5P=1;MEI_3P=288;G1K=nssv14066827;LR=chr22-47331882-INS->s904797<s909378<s909377<s909376>s904798-324;KNOWN_SRC=melt_1kg|long_read_1kg_ont_vienna;SAMPLE_STATUS=shared;SUPPORT=SR_L:8|SR_R:1|DPE_L:13|DPE_R:13|MEI_MAPPED:30|polyA_MAPPED:11|polyA_side:R;L1TXGOLDSCORE=1.000;L1TXRANKPCT=99;INSERTIONSCORE=0.643;PEAKDEPTHZ=0.0478;SPLITCLUSTERZ=10.4;CALLTIER=none;KNOWNMEI=True	GT:GQ	./.:.
chr22	49029650	L1TX-chr22-49029650-SVA	N	<INS:ME:SVA>	.	PASS	SVTYPE=INS;END=49029650;WINDOWSTART=49029645;WINDOWEND=49029656;MEIFAMILY=SVA;MEISUBFAMILY=SVA_D#Retroposon/SVA;TSD=AAGAAAACTCCT;POLYA_MIN_BP=116;ORIENT=-;NESTED=unnested;MEI_SPAN=1366;MEI_5P=1;MEI_3P=1366;G1K=nssv14064350;KNOWN_SRC=melt_1kg;SAMPLE_STATUS=shared;SUPPORT=SR_L:0|SR_R:0|DPE_L:57|DPE_R:9|MEI_MAPPED:61|polyA_MAPPED:27|VNTR_MAPPED:0|polyA_side:L;L1TXGOLDSCORE=1.000;L1TXRANKPCT=99;INSERTIONSCORE=0.478;PEAKDEPTHZ=0.056;SPLITCLUSTERZ=13.3;CALLTIER=none;KNOWNMEI=True	GT:GQ	./.:.
chr22	50495454	L1TX-chr22-50495454-ALU	N	<INS:ME:ALU>	.	PASS	SVTYPE=INS;END=50495454;WINDOWSTART=50495444;WINDOWEND=50495464;MEIFAMILY=ALU;MEISUBFAMILY=AluYa5#SINE/Alu;TSD=TCTCCTGCCTCACCCTCCCCA;POLYA_MIN_BP=94;ORIENT=-;NESTED=nested;MEI_SPAN=270;MEI_5P=12;MEI_3P=281;SAMPLE_STATUS=shared;SUPPORT=SR_L:3|SR_R:2|DPE_L:56|DPE_R:6|MEI_MAPPED:41|polyA_MAPPED:16|polyA_side:L;L1TXGOLDSCORE=0.9984;L1TXRANKPCT=90;INSERTIONSCORE=0.441;PEAKDEPTHZ=0.00706;SPLITCLUSTERZ=-0.228;CALLTIER=none;KNOWNMEI=False	GT:GQ	./.:.
```

## Examples

`docs/examples/HG03086_chr22_gold_review.vcf` is the raw gold-review output
for this chromosome (410 calls).
`docs/examples/HG03086_chr22_classifier_ge_0.997.vcf` keeps calls scored by a
classifier trained on 1000 Genomes; that model is described in
[pull request #69](https://github.com/L1-Therapeutics/retrotransposon-miner/pull/69).

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

Soft-clips and discordant clipped ends remap to the Dfam Alu/LINE-1/SVA panel with `bwa mem -k10 -T10` (`--bwa-threads`; wrapper auto: `nproc` single-chrom, `1` under multi-chrom concurrency). Nearby same-chrom clipped DPE mates remap the clip only; interchrom mates and same-chrom mates ≥500 kb away remap the reference-aligned body so a mate sitting on a reference-copy Alu/L1/SVA still counts as MEI identity. Coherent same-chrom deletion bridges are removed separately before MEI support is counted. Annotate writes fetched mate sequences to `discordant_mate_cache.*.parquet` beside the extract so later remaps skip CRAM. Queries are polyA/T-trimmed (≥8 bp) before align; consensus targets are also terminal-polyA-trimmed (sidecar `*.nopolya.fa`, and prep writes body-only Dfam/panel FASTAs) so clips cannot map onto the A-tail. Junction clips that are themselves polyA/T count as `polyA_MAPPED` only — never also `MEI_MAPPED`/SR. Short tips (≤30 bp) need qcov≥0.80 and pid≥0.90, longer clips need pid≥0.90 and alnlen≥20. Panel fragment hits project onto one family-consistent `*_full` axis via `mei_fragment_to_full_coords.tsv` (prep: `bwa mem -a` on trimmed sequences). Assembly contig-to-MEI still uses minimap2. Benchmark: `scripts/benchmark_mei_aligners.py`.

## Getting Started on Amazon EC2 (Elastic Compute Cloud)

### Whole-genome machine and disk

One 30× short-read genome (disease and control pointed at the same alignment) fits on a **64 vCPU / 256 GiB** machine with a **200 GB gp3** root volume. The example used here is an on-demand `m7i.16xlarge` in `us-east-1`: 64 vCPU, 256 GiB, up to 20 Gbps to EBS. The volume is 200 GB, 16000 IOPS, and 1000 MB/s throughput (gp3's throughput cap; 1000 MB/s requires at least 4000 IOPS).

```bash
INSTANCE_TYPE=m7i.16xlarge \
ROOT_VOLUME_GB=200 \
ROOT_VOLUME_IOPS=16000 \
ROOT_VOLUME_THROUGHPUT_MB=1000 \
S3_BUCKET=s3://<your-bucket> \
./scripts/ec2_jlab.sh bootstrap
```

Convert the CRAM to one coordinate-sorted BAM before the run (`samtools view -@ 24 -b -T ref.fa -o sample.bam sample.cram`, then `samtools index -@ 24`). Pass that BAM as both `--disease-bam` and `--control-bam`. Extract, mate fetch, peak-depth, and IGV then read it in place. IGV converts a path only when it ends in `.cram`. A 30× CRAM of about 15 GB becomes a BAM of about 40 GB. With the reference (~12 GB) and per-chromosome tables, 200 GB still has room to keep the CRAM until the BAM is indexed.

Run a full genome with `--chr all --chr_concurrency 16`. `--chr all` expands to chrX, chrY, then chr1 through chr22, so chromosome X starts in the first 16 slots. Chromosome X is about the length of chromosome 8. Each chromosome is mostly one thread, so extra cores past the chromosome count do not shorten the longest chromosome. Sixteen leaves memory for the jobs and for caching the BAM; 24 is the chromosome count and the useful ceiling. A larger instance does not finish faster than chromosome 1.

One 30× genome on this machine took **1 hour 45 minutes** after the BAM was indexed. Chromosome 2 was the longest chromosome and set that time. The 24-thread CRAM-to-BAM conversion took 4 minutes, so the job from the start of conversion was **1 hour 49 minutes**.

A 64-vCPU on-demand instance uses the whole default standard-family vCPU quota (64) on a new account. Stop other A/C/D/H/I/M/R/T/Z instances before launch.

A disease and normal pair (two different alignments) wants the same 64 vCPU / 256 GiB shape and a **300 GB** gp3 volume at the same 16000 IOPS and 1000 MB/s. Budget two ~40 GB BAMs plus the reference and two evidence tables. Use `--chr_concurrency 12` so both BAMs can stay cached next to the chromosome jobs.

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

Use `bootstrap` only when you want the script to provision a new instance (key pair, security group, Elastic IP, JupyterLab). On a shared AWS account, each IAM user gets their own key pair (`retrotransposon-miner-<region>-<iam-user>`). If `~/.ssh/id_ed25519.pub` or `id_rsa.pub` exists, that public key is imported — bootstrap does not reuse another user’s PEM.

`bootstrap` launches an **on-demand** `m7i.4xlarge` by default (16 vCPU / 64 GiB). A full genome uses `INSTANCE_TYPE=m7i.16xlarge` with the disk settings in the section above. `SPOT=1` switches to cheaper Spot (one-time, stop-on-interruption: reclaim keeps the EBS disk and does not start the instance again). Bring it back with `start-instance`. Launch does **not** pin an AZ; AWS places the instance in a default-VPC zone that has capacity. `SUBNET_ID` pins a subnet (and therefore an AZ). `start-instance` cannot change AZ; if start fails for capacity, retry later or `bootstrap` a new VM.

```bash
S3_BUCKET=s3://<your-bucket> ./scripts/ec2_jlab.sh bootstrap
SPOT=1 S3_BUCKET=s3://<your-bucket> ./scripts/ec2_jlab.sh bootstrap   # cheaper Spot
```

That uses your **local** AWS CLI profile only on your laptop, to:

1. Create/reuse an IAM role + instance profile scoped to that bucket.
2. Attach the profile to the new instance (the VM assumes the role via instance metadata; **do not copy** `~/.aws` keys onto the instance).
3. Write `RTM_S3_CACHE=s3://<your-bucket>/public` on the instance.

For an instance that is already running:

```bash
S3_BUCKET=s3://<your-bucket> ./scripts/ec2_jlab.sh attach-s3
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
| `attach-s3` | Grant the bound instance IAM access to `S3_BUCKET` |
| `install-my-key` | Push this laptop’s public key onto the bound instance (Instance Connect) |
| `start-jlab` / `stop-jlab` / `start-tunnel` | JupyterLab lifecycle |
| `help` | Show usage |

Optional environment variables:

- `REGION` — AWS region (default: from `aws configure`, else `us-east-1`)
- `INSTANCE_ID` / `INSTANCE_NAME` — override bound instance without editing `.ec2-instance.env`
- `HOST_ALIAS` — SSH config alias (default: `retro-ec2`)
- `SSH_USER` — SSH login user (auto-detected from AMI if unset; e.g. `ec2-user`, `ubuntu`)
- `KEY_PATH` — path to the private key for the instance (PEM or `~/.ssh/id_ed25519`)
- `KEY_NAME` / `KEY_OWNER` — override the per-user EC2 key pair name (default: `retrotransposon-miner-<region>-<iam-user>`)
- `INSTANCE_TYPE` — instance type for `bootstrap` only (default: `m7i.4xlarge`)
- `SPOT` — `0` (default) launches on-demand; `1` launches a Spot instance (`bootstrap` only)
- `SUBNET_ID` — pin `bootstrap` to one subnet/AZ. Default: omit subnet so AWS chooses an AZ with capacity.
- `ROOT_VOLUME_GB` — root EBS size for `bootstrap` only (default: `200`)
- `ROOT_VOLUME_IOPS` / `ROOT_VOLUME_THROUGHPUT_MB` — gp3 IOPS and MB/s for `bootstrap` only (default: `4000` / `1000`). AWS requires throughput ≤ 0.25 × IOPS; 1000 MB/s needs at least 4000 IOPS. The gp3 baseline (3000 / 125) is the usual WGS stage bottleneck.
- `S3_BUCKET` — bucket the instance may read/write (example: `s3://<your-bucket>`); creates/reuses an instance profile, does not copy local keys
- `S3_CACHE_PREFIX` — object prefix for public-data cache (default: `s3://<bucket>/public`)
- `IAM_INSTANCE_PROFILE` / `IAM_ROLE_NAME` — override the default `ec2-retrotransposon-s3-profile` / `ec2-retrotransposon-s3-role`

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
- `ec2:CreateKeyPair` / `ec2:ImportKeyPair` (per-IAM-user key; imports `~/.ssh/id_ed25519.pub` when present)
- `ec2-instance-connect:SendSSHPublicKey` (break-glass: `./scripts/ec2_jlab.sh install-my-key`)
- `ec2:CreateSecurityGroup`
- `ec2:AllocateAddress`
- `ec2:AssociateAddress`
- `ec2:DescribeAddresses`
- `ec2:DescribeVpcs`
- `ec2:DescribeSubnets`
- `ec2:AssociateIamInstanceProfile` / `ec2:ReplaceIamInstanceProfileAssociation` (when `S3_BUCKET` is set)
- `ssm:GetParameter` (Amazon Linux AMI lookup)
- `iam:PassRole` / `iam:GetRole` / `iam:GetInstanceProfile` (reuse the existing instance role; do not copy laptop keys onto the VM)
- `iam:CreateRole` / `iam:CreateInstanceProfile` / `iam:AddRoleToInstanceProfile` / `iam:PutRolePolicy` (first-time admin setup only; skipped when the role already exists)
- `s3:CreateBucket` / `s3:HeadBucket` / `s3:ListBucket` (optional; only if the cache bucket does not exist yet)

### What `scripts/ec2_jlab.sh` Does

- Binds to any existing EC2 instance by ID or `Name` tag (`use`).
- Saves the binding locally in `.ec2-instance.env` (not committed to git).
- Starts, stops, and reboots the bound instance without creating new ones.
- Writes SSH aliases (`retro-ec2`, `jlab`) into local `~/.ssh/config`.
- Refreshes SSH security group ingress for your current public IP on connect.
- Optionally creates a new instance (`bootstrap`), a **per-IAM-user** key pair (or imports your laptop `id_ed25519.pub`), security group, and Elastic IP. New instances are **on-demand** by default (`SPOT=1` for cheaper Spot). Does not reuse another user's PEM.
- `install-my-key` pushes your laptop public key via EC2 Instance Connect and appends it to `authorized_keys`.
- Optionally attaches an IAM instance profile for a user-specified `S3_BUCKET` and caches public data under `s3://<bucket>/public`.
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

### Export to VCF

Convert the annotated candidate-loci table into standard VCF v4.3 so the
callset can be fed to downstream tools (`bcftools`, IGV, annotation
pipelines):

```bash
python -m retro_miner.cli export-vcf \
  --in-tsv "${RTM_RESULTS_DIR:-$HOME/retrotransposon-workdir/results}/quickstart_seqc2_chr22/candidate_loci.mei.tsv" \
  --out-vcf "${RTM_RESULTS_DIR:-$HOME/retrotransposon-workdir/results}/quickstart_seqc2_chr22/candidate_loci.mei.vcf" \
  --sample-name seqc2_tumor_normal
```

The same command accepts a genome-wide gold review table. A classifier-ranked
table has no breakpoint column (`mei_family` is `Alu`/`L1`/`SVA`, and
`gold_score` is the classifier probability). Join it to the gold table and
keep a score cutoff:

```bash
python -m retro_miner.cli export-vcf \
  --in-tsv hg03086_gold_by_classifier_score.tsv \
  --breakpoint-tsv candidate_loci.mei.gold_review.tsv \
  --min-score 0.997 \
  --out-vcf HG03086.classifier_ge_0.997.vcf
```

`L1TXGOLDSCORE` is the classifier probability at 4 significant figures.
`L1TXRANKPCT` is the percentile of `classifier_rank` (100 is best, nearest
integer). `L1TXGOLDRANKPCT` is written only when the table has no classifier
rank. `--chr all` writes `candidate_loci.mei.gold_review.vcf` after the genome
gold table is aggregated. `##reference` is the run's `--reference-build`
value, read from `pipeline_params.env`, and `##assembly` names that build
once (for example `GRCh38`). Contig lines stay `##contig=<ID=chr1>`. A
classifier export inherits that header from the `--breakpoint-tsv` gold
table when the classifier file is outside the run directory. The VCF ID is
`L1TX-<chrom>-<pos>-<family>`. Overlapping catalog accessions are
`G1K` (the 1000 Genomes MELT id, which is the dbVar nssv accession) and
`LR` instead of the ID column. When disease and control support are the
same string, the VCF writes `SUPPORT` and omits the duplicate disease
field. Distinct disease and control samples keep `CTRL_SUPPORT` and
`DISEASE_SUPPORT`. `SVLEN` is
left unset so Ensembl VEP does not treat the insertion as a reference span;
the element length stays in `MEI_SPAN`, and `END` equals `POS`.

All 21 chromosome 22 classifier calls are listed under Example Variant Calls. The header and both chromosome 22 files are linked from Examples.

Records are coordinate-sorted, so the output can be compressed and indexed
directly:

```bash
bgzip -c candidate_loci.mei.vcf > candidate_loci.mei.vcf.gz
bcftools index candidate_loci.mei.vcf.gz
bcftools view -r chr22:19000000-32000000 candidate_loci.mei.vcf.gz
```

Insertions are emitted as symbolic ALT alleles (`<INS:ME:ALU>`,
`<INS:ME:LINE1>`, `<INS:ME:SVA>`), following the 1000 Genomes MEI VCF
convention; `REF` is `N` because `POS` marks the insertion breakpoint
rather than a called reference base. MEI family/subfamily, TSD sequence,
poly-A tail length, orientation, nesting status, full-length MEI
span/coordinates, known-polymorphism cross-references, and the raw
per-cohort supporting-read evidence strings are carried in `INFO`.

**Genotype fields are intentionally left blank** (`GT=./.`, `GQ=.`). This
table reports pooled disease-vs-control read support across a cohort
comparison, not per-individual diploid genotypes, so there is no genotype
to estimate. See "Current Limitations" below.

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
  --outdir "${RTM_PUBLIC_DATA_DIR:-$HOME/retrotransposon-workdir/data/public}" \
  --s3-cache-prefix "${RTM_S3_CACHE:-s3://l1tx-data/public}"
```

Test BAMs default to a **chr22 slice** (plus discordant mates). To park the entire SEQC2 WGS pair on S3 without filling the instance disk, stream NCBI HTTP through the instance (`curl | aws s3 cp -`) — the object never lands on the local volume:

```bash
python3 scripts/download_public_data.py \
  --references hg38 \
  --dataset-ids seqc2_disease_bam seqc2_control_bam \
  --test-bam-mode full \
  --s3-cache-prefix "${RTM_S3_CACHE:-s3://l1tx-data/public}" \
  --skip-postprocess

# Later chrom slices reuse that S3 object (no NCBI full-BAM scan):
python3 scripts/download_public_data.py \
  --references hg38 \
  --dataset-ids seqc2_disease_bam seqc2_control_bam \
  --test-bam-mode slice \
  --test-bam-chrom chr1 \
  --s3-cache-prefix "${RTM_S3_CACHE:-s3://l1tx-data/public}" \
  --skip-postprocess
```

Chromosome slices pull a local `.bai` (from S3 or the NCBI sidecar) and use `samtools view -X` for the region plus discordant-mate windows. They do not stream the whole BAM with `samtools view -N`.

`--test-bam-mode full` never writes the ~200 GiB BAMs to local disk. Add `--slice-after-full` only when you also want a local `--test-bam-chrom` slice.

S3 copies (cache sync, BAM/BAI staging, S3→S3) use 64 concurrent 64 MiB multipart parts instead of the AWS CLI 10×8 MiB default. Turn concurrency down with `RTM_S3_MAX_CONCURRENCY` (process-local AWS config; `~/.aws/config` is not rewritten). Staging a WGS pair to EBS still needs volume throughput: `bootstrap` gp3 defaults to 1000 MB/s (4000 IOPS). The 125 MB/s gp3 baseline will hold a 200 GiB pair at ~25 min even with the fast copier.

The candidate pipeline **does** stage remote (`s3://` or `http(s)://`) disease/control BAMs to `${RTM_WORKDIR}/data/bam_stage` when the run is multiple chromosomes, `--chr all`, or `--chr_concurrency > 1`. It skips the copy when the local file already matches the remote size, and refuses to start if free disk is below BAM size plus headroom. Single-chromosome runs keep streaming. Override the dest with `--bam-stage-dir` / `RTM_BAM_STAGE_DIR`, or disable with `--no-bam-stage` / `RTM_BAM_STAGE=0`.

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
- Set `KEY_PATH` if your private key is not `~/.ssh/id_ed25519` or `~/.ssh/<key-name>.pem`.
- Locked out of an instance you created? `./scripts/ec2_jlab.sh use <instance-id>` then `./scripts/ec2_jlab.sh install-my-key` (needs `ec2-instance-connect:SendSSHPublicKey`). Do not share another user's PEM.
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
