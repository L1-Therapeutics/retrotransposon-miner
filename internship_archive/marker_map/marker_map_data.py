"""
marker_map_data.py — The Retrotransposon Marker Map: single source of truth.

Schema
  NODES: (id, label, section, kind, evidence, source, views)
  EDGES: (src, dst, tier, label, views)
  kind      : chain | marker | resource | gap | clinical | structure
  evidence  : demonstrated | associated | proposed | unknown | gap
  views     : subsets the node/edge appears in ('spine' and/or 'full')
  tier      : edge evidence tier, styled in the renderer
Verification flags inside labels: [V]=page-verified primary text, [S]=snippet-level,
[C]=computed from published percentages, [W]=our own data. No number without a source.
"""

SECTIONS = {
    "SPINE": ("CAUSAL SPINE",                       "#E8E8E8"),
    "A":     ("A  TARGETING (where it lands)",      "#D6EAF8"),
    "B":     ("B  CREDIBILITY (is the call real)",  "#D5F5E3"),
    "C":     ("C  ELEMENT STATE (is it armed)",     "#FCF3CF"),
    "D":     ("D  MEDICAL & DISEASE",               "#FDE2E4"),
    "E":     ("E  RESOURCES (known sites)",         "#E8DAEF"),
    "G":     ("G  OUR DATA",                        "#FDEBD0"),
    "X":     ("X  VERIFIED GAPS (hexagons)",        "#F5B7B1"),
}

KIND_SHAPE = {"chain": "box", "marker": "ellipse", "resource": "folder",
              "gap": "hexagon", "clinical": "box3d", "structure": "diamond"}

EVIDENCE_STYLE = {
    "demonstrated": ("solid",  "black"),
    "associated":   ("dashed", "#1F618D"),
    "proposed":     ("dotted", "#7D6608"),
    "unknown":      ("dotted", "gray"),
    "gap":          ("dashed", "#C0392B"),
}

SPINE_IDS = {"meth", "rna", "orf2p", "en", "tprt", "position", "junction", "lineage", "reseal", "clinical"}

# =====================================================================================
#  NODES
# =====================================================================================
NODES = [
    # ---------------- SPINE: the causal chain ----------------
    ("meth",     "CpG methylation +\nH3K9me3 / SETDB1\n[permissive, not sufficient]", "SPINE", "chain", "demonstrated", "Lanciano 2024 NAR; Karlsson 2026", ["full", "spine"]),
    ("rna",      "L1/Alu RNA\n(61/379 loci express\nafter demethylation)",            "SPINE", "chain", "demonstrated", "Lanciano 2024 Fig 6/S8", ["full", "spine"]),
    ("orf2p",    "ORF1p+ORF2p RNP\n(RNA<->ORF1p\nrho = 0.55-0.79)",                   "SPINE", "chain", "demonstrated", "McKerrow 2022 CPTAC", ["full", "spine"]),
    ("en",       "ORF2p-EN finds target:\nTTTTT/AA-class 7-mer\n(6% of 7-mers = 97% of sites)", "SPINE", "chain", "demonstrated", "Flasch 2019 Cell; Feng 1996 Cell", ["full", "spine"]),
    ("tprt",     "TPRT: nick1 -> 3'OH primer;\nnick2 +15-16 bp downstream\n(TSD = gap - A-rich)", "SPINE", "chain", "demonstrated", "Levy 2010 NAR p.1518; Thawani 2024 Nature 8UW3; Ghanim 2025 Science", ["full", "spine"]),
    ("position", "INSERTION POSITION:\nlinker 133 (validated) |\n3' tail 280-300 (ours, validating)", "SPINE", "chain", "unknown", "Levy 2010 Fig 2; OUR Q4", ["full", "spine"]),
    ("junction", "JUNCTION ARCHITECTURE:\nTSD + 3'A-tail + route tags",               "SPINE", "chain", "demonstrated", "Wagstaff 2012 (~96% triad)", ["full", "spine"]),
    ("lineage",  "LINEAGE TAG:\ntransduced flank -> source element",                   "SPINE", "chain", "demonstrated", "Goodier 2000; Damert 2009; Tubio 2014", ["full", "spine"]),
    ("reseal",   "Methylation re-seal\n(YY1 proposed;\nno timing clock measured)",     "SPINE", "chain", "proposed", "Lanciano 2024", ["full", "spine"]),
    ("clinical", "CLINICAL READOUT:\nplasma ORF1p; cfDNA hypometh;\nAlu quantity/fragmentation", "SPINE", "chain", "associated", "Taylor 2023; Nagai 2017", ["full", "spine"]),

    # ---------------- A. TARGETING ----------------
    # A1. the motif
    ("a_motif7",    "EN 7-mer consensus\nTTTTT/AA: 9.6% raw use;\nrank 21 after opportunity norm", "A", "marker", "demonstrated", "Flasch 2019 (64,973 calls, 4 cell lines)", ["full"]),
    ("a_motifC",    "45% of sites carry 1 C\nwithin motif positions 2-5\n(degenerate, not a code)", "A", "marker", "demonstrated", "Flasch 2019", ["full"]),
    ("a_motifopp",  "743/12,288 7-mers (6%)\naccount for 97% of insertions;\noccupy ~23% (~750 Mb) of genome", "A", "marker", "demonstrated", "Flasch 2019", ["full"]),
    ("a_nick2",     "Second nick ~15-16 bp down\n(TNTN/AA class)",                   "A", "marker", "demonstrated", "Levy 2010 p.1518", ["full"]),
    ("a_polytsd",   "polyA/TSD inverse relation\nSpearman rho = -0.21",              "A", "marker", "demonstrated", "Levy 2010", ["full"]),
    # A2. DNA-shape / structural targeting (cryo-EM era)
    ("a_ssDNA",     "Upstream ssDNA positions\nduplex in EN active site",            "A", "structure", "demonstrated", "Thawani 2024 (8UW3)", ["full"]),
    ("a_dnashape",  "EN recognition driven by\nDNA conformation, not just bases\n(compression, not APE1-bend)", "A", "structure", "demonstrated", "Miller 2021 NAR (7N94)", ["full"]),
    ("a_wedge",     "CTD zinc-finger wedges\nbetween target strands;\nORF2p bends/unzips target", "A", "structure", "demonstrated", "Ghanim 2025 (9HDQ/R)", ["full"]),
    ("a_stag",      "In-vitro top-strand cuts at\n+7/+11/+19 nt from nick1\n(-> variable TSD lengths)", "A", "structure", "demonstrated", "Ghanim 2025", ["full"]),
    ("a_ends",      "EN cleaves near DNA ends\n& at mismatches (40x faster\non mismatch substrate, in vitro)", "A", "marker", "associated", "Miller 2025 NAR [S]", ["full"]),
    # A3. composition windows
    ("a_gc50",      "50-bp cut window: 32% GC\nvs 41% genome\n(AT-rich FOOTPRINT)",  "A", "marker", "demonstrated", "Gasior 2006 (n=103)", ["full"]),
    ("a_gc20kb",    "20-kb flanks: 41% GC = genome\n(Alu content 13.0% vs 10.6%)\nregional GC-NEUTRAL", "A", "marker", "demonstrated", "Gasior 2006", ["full"]),
    ("a_alu1mb",    "Alu 1-Mb windows w/ inserts:\nGC 0.4125 vs 0.3984;\ngene content 0.3027 vs 0.2164", "A", "marker", "demonstrated", "Wagstaff 2012 (n=226)", ["full"]),
    ("a_svagc",     "SVA de novo: G+C-rich\nregional preference (no % in text)",     "A", "marker", "associated", "Raiz 2011 J Virol [S]", ["full"]),
    # A4. genic targeting
    ("a_genic_alu", "Alu de novo: 57.7% in genes (130/226);\nof those 96.2% intronic (125/130)", "A", "marker", "demonstrated", "Wagstaff 2012", ["full"]),
    ("a_genic_l1",  "L1 de novo: intronic 32.5-41.4%\nbut NOT enriched vs seq-matched null;\nseveral cell types depleted", "A", "marker", "demonstrated", "Flasch 2019", ["full"]),
    ("a_genic_sva", "SVA ref copies: 51% intronic\n(2,618/5,107); gene-dense clustering", "A", "marker", "demonstrated", "Chu 2023 NAR", ["full"]),
    ("a_exonprox",  "Exon-proximity exclusion zones:\nSINE 100 bp, LINE 50 bp,\nnear-full LINE 2,000 bp", "A", "marker", "demonstrated", "Zhang 2011 [V]", ["full"]),
    ("a_intronic",  "Sela 2007: intronic L1 65.7%\nantisense (185,428 vs 96,718);\nAlu 66% intronic, no SINE bias", "A", "marker", "demonstrated", "Sela 2007 [V]", ["full"]),
    ("a_lavialu",   "Alu by position: intronic 54.9% ant;\nat splice junctions 86.1% ant;\nat gene ends 77.7% SENSE", "A", "marker", "demonstrated", "Lavi & Carmel 2018 (579,747 Alus)", ["full"]),
    ("a_levmaor",   "Antisense-Alu depletion within\n150 bp upstream of exons\n(p=0.000012; downstream p=0.056)", "A", "marker", "demonstrated", "Lev-Maor 2008", ["full"]),
    # A5. positional hotspots
    ("a_133",       "LINKER HOTSPOT pos 133\ndirectly after A-rich linker\n(AluJo linker 118-136/133)", "A", "marker", "demonstrated", "Levy 2010 Fig 2 p.1518 [V]", ["full"]),
    ("a_133ours",   "OUR replication:\n120-140 bp bin, 3.0x, n=29\nA-frac 0.55; chi2=103 p=4e-18", "A", "marker", "demonstrated", "OUR HG03086 Q4 [W]", ["full"]),
    ("a_tail",      "3' POLY-A TAIL PEAK\n280-300 bp: 3.7x, n=36, A-frac 0.79\nfield blind spot (Levy excluded ends)", "A", "marker", "unknown", "OUR Q4 [W]; under validation", ["full"]),
    ("a_5p",        "5'-end peak: n=13, A-frac 0.13\n(edge-effect candidate)",       "A", "marker", "unknown", "OUR Q4 [W]", ["full"]),
    ("a_sense",     "Same orientation ~85%\nflat across host age (0.82-0.88)",       "A", "marker", "demonstrated", "OUR Q2/Q3 [W]", ["full"]),
    ("a_antis",     "Alu antisense DEPLETED 0.46x (p=6e-8);\nL1 antisense at chance (p=0.32)", "A", "marker", "demonstrated", "OUR Q2 [W]", ["full"]),
    ("a_homol",     "Homology gradient: AluY children\n-> AluY hosts 1.86x, AluS 0.88x,\nAluJ 0.72x (graded by similarity)", "A", "marker", "demonstrated", "OUR Q3 [W]; Jacob-Hirsch 2018 L1", ["full"]),
    ("a_ccatt",     "CCATT motif NOT enriched\nat our breakpoints (classical germline\nroute, not somatic brain route)", "A", "marker", "demonstrated", "OUR [W]; Jacob-Hirsch 2018 (brain: CCATT+)", ["full"]),
    # A6. replication & chromatin context
    ("a_repl",      "Replication timing: cell-type\ndependent (NPC/PA-1 late+, hESC late-)\nno single direction", "A", "marker", "demonstrated", "Flasch 2019 Fig 7", ["full"]),
    ("a_fork",      "Leading-strand template bias\n(EN-dependent; flips in FANCD2-/-)", "A", "marker", "demonstrated", "Flasch 2019", ["full"]),
    ("a_chrom",     "Chromatin states: no enrichment\n(15 states); transcribed regions\ndepleted (32.6% HeLa, 19.4% PA-1)", "A", "marker", "demonstrated", "Flasch 2019", ["full"]),
    ("a_nuc",       "Nucleosome depletion at targets\nlargely = AT-richness of motif", "A", "marker", "associated", "Sultana 2019 Mol Cell [S]", ["full"]),
    ("a_consensus", "Alu core 280-282 bp (AluSx 1-282);\nfull element ~300 bp + variable tail;\nleft monomer ~1-120, right ~134-282", "A", "marker", "demonstrated", "Price 2004; Konkel 2015; Lavi & Carmel 2018 [V]", ["full"]),

    # ---------------- B. CREDIBILITY ----------------
    # B1. triad
    ("b_triad",     "TPRT TRIAD (TSD + 3'A-tail\n+ EN motif): ~96% of de novo Alu", "B", "marker", "demonstrated", "Wagstaff 2012 PLoS Genet", ["full"]),
    ("b_tsd_alu",   "TSD Alu de novo: 5-27 bp,\nmean 14.0 +/- 3.0 (n=226)",          "B", "marker", "demonstrated", "Wagstaff 2012", ["full"]),
    ("b_tsd_l1",    "TSD L1 genomic: bimodal 9/15 bp\n(16,266 calls) -> 16 bp after\nfiltering (short = chance matches)", "B", "marker", "demonstrated", "Szak 2002; Zingler 2005", ["full"]),
    ("b_tsd_l1fl",  "TSD full-length L1: 14-16 bp peak",                             "B", "marker", "demonstrated", "Kojima 2010 Mobile DNA [V]", ["full"]),
    ("b_tsd_sva",   "TSD SVA: 4-20 bp\n(no published histogram - gap)",             "B", "marker", "demonstrated", "Hancks 2012; Ostertag 2003", ["full"]),
    ("b_tsd_cat",   "TSD presence in catalogs:\nAlu 55.6% (9,756/17,543);\nL1 33.0%; SVA 32.4%", "B", "marker", "demonstrated", "Niu 2022 [V] (older 1KGP set)", ["full"]),
    ("b_tsd_lr",    "Long-read hallmark-bearing:\n96.4% (29,603/30,718)",            "B", "marker", "demonstrated", "Schloissnig preprint [S]", ["full"]),
    ("b_tsdshort",  "Short-TSD caution (9-11 bp)\nis L1-specific; imperfect longer\ndups miscall as short perfect TSDs", "B", "marker", "demonstrated", "Szak 2002; Kojima 2010", ["full"]),
    # B2. junction route tags
    ("b_micro",     "5' junction microhomology:\nL1-trunc 62% (1,503/2,427);\nL1-full 24%; ALU ~NONE (10,062)", "B", "marker", "demonstrated", "Zingler 2005 [V]", ["full"]),
    ("b_inv",       "5' inversion/twin priming:\n~25% recent genomic L1;\n6/38 de novo (Symer); ALU = UNKNOWN", "B", "marker", "demonstrated", "Ostertag 2001 GR 11:2059 [V]; Symer 2002 [V]", ["full"]),
    ("b_invporub",  "Polymorphic L1: 26% (330/1,271)\ninversions; breakpoints cluster\nat L1 bases 4,000-6,000 (88%)", "B", "marker", "demonstrated", "Porubsky 2022 [V]", ["full"]),
    ("b_ip",        "INTERNAL PRIMING (alt route):\n20 verified events (2 Alu, 4 L1 human);\n3'-trunc, no polyA, TSD>=6", "B", "marker", "demonstrated", "Srikanta 2009 AJHG [V]", ["full"]),
    ("b_del",       "Target-site DELETION:\n13/100 de novo L1 (2 bp - >14 kb);\ncell-culture max 71 kb", "B", "marker", "demonstrated", "Gilbert 2005 MCB [V]; Han 2005", ["full"]),
    ("b_filler",    "Filler DNA: 6 single-base + 7\n(4-20 bp) per 100 de novo L1;\nextranucs: 14.6% trunc-L1, 11.7% trunc-Alu", "B", "marker", "demonstrated", "Gilbert 2005 [V]; Kojima 2010 [V]", ["full"]),
    ("b_trans_l1",  "L1 3' TRANSDUCTION: 23% (15/66,\nmean 207 nt) selected set; 8.8-15%\ngenomic; 24% of somatic", "B", "marker", "demonstrated", "Goodier 2000 HMG 9:653 [V]; Pickeral 2000; Tubio 2014", ["full"]),
    # Damert 2009 as published.  A measured boundary is recorded for this node in
    # render_marker_map_v5.py (the correction ledger): our own calls cannot express a
    # 5'-transduced copy at all and name no source element, and the one SVA 5' junction
    # with reads gives clips of 2-100 bp with nothing inside them to subtract.  Written
    # by sva_5prime_transduction_probe.py, the SVA arm of l1_readthrough_probe.py's test.
    ("b_trans_sva", "SVA 5' TRANSDUCTION: 220/2,398 (9.2%);\nlengths 14-2,161 bp, mean 223;\nMAST2 group >=32% of SVA-F", "B", "marker", "demonstrated", "Damert 2009 GR 19:1992 [V]", ["full"]),
    ("b_rearr",     "Rearranged inserts: Alu 8/226 (3.5%)\n= 6 chimeric + 2 EN-negative;\nL1: 94/100 truncated, 19 inv/del", "B", "marker", "demonstrated", "Wagstaff 2012 [V]; Gilbert 2005 [V]", ["full"]),
    ("b_trunc",     "Truncation asymmetry:\nL1 genomic 85% truncated;\nALU only 8/226 (3.5%) truncated", "B", "marker", "demonstrated", "Zingler 2005; Wagstaff 2012", ["full"]),
    ("b_recomb",    "Post-insertion recombination:\n73 L1-del (~447 kb, median 3,239);\n492 Alu-del (~400 kb)", "B", "marker", "demonstrated", "Han 2008 PNAS [V]; Sen 2006 AJHG [V]", ["full"]),
    ("b_source",   "Source competence: 90 intact -> 82 cloned\n-> 40 active (49%); 6 hot = 84%\nof all activity", "B", "marker", "demonstrated", "Brouha 2003 [V]; Beck 2010 [V]", ["full"]),
    # B3. detection (the EMPTY fingerprint niche)
    ("b_callers",   "MEI callers = read evidence:\nalu-detect 97%/85% sim; Mobster 13bp\nmedian TSD; TANGRAM; MELT; assembly FDR<5%", "B", "marker", "demonstrated", "David 2013; Thung 2014; Lee 2014; Gardner 2017; Wildschutte 2015", ["full"]),
    ("b_niche",     "NO CALLER BENCHMARKS the Levy\nfingerprint (hotspot/orientation/triad)\nas ranking features -> OUR PHASE 3", "B", "gap", "gap", "verified absence [V]", ["full"]),

    # ---------------- C. ELEMENT STATE ----------------
    ("c_chain1",    "Locus methylation -> expression:\nweak rho; permissive NOT sufficient\n(only 61/379 express post-demethylation)", "C", "marker", "demonstrated", "Lanciano 2024 Fig 6+S8 [V]", ["full"]),
    ("c_chain2",    "RNA -> ORF1p rho 0.55-0.79;\nORF1p -> somatic insertions\nrho 0.38 (p=0.004, endometrial)", "C", "marker", "demonstrated", "McKerrow 2022 [V]", ["full"]),
    ("c_tubio",     "Source-locus hypomethylation\n<-> source activity during\ntumor evolution (290 pairs, 12 types)", "C", "marker", "associated", "Tubio 2014 Science [V]", ["full"]),
    ("c_freeman",   "Expressed L1 loci = hypometh +\nopen chromatin + active marks\n('transitional' loci identified)", "C", "marker", "demonstrated", "Freeman 2022 NAR [V]", ["full"]),
    ("c_byun",      "Young-Alu epigenetics: Yb8 >=8x\nH3K9me3 vs AluJ; monoallelic Ya5\nmethylates locally, represses host gene", "C", "marker", "demonstrated", "Byun 2012 PNAS [V]", ["full"]),
    ("c_setdb1",    "SETDB1/H3K9me3 perturbation\nactivates young L1s AND loses\npromoter methylation (hNPC)", "C", "marker", "demonstrated", "Karlsson 2026 [V]", ["full"]),
    ("c_youngctx",  "Young-L1PA hypomethylation ONLY in:\nPSC, trophoblast, placenta, seminoma,\nEC cells - NOT most carcinomas", "C", "marker", "demonstrated", "Lanciano 2024 (~12k GEO 450K) [V]", ["full"]),
    ("c_pgc",       "Developmental resets: PGC\n16% -> 4.5% CpG; >90% of escapee\nregions TE-linked; ICM ~30-35%", "C", "marker", "demonstrated", "Lee & Surani 2024 [V]", ["full"]),
    ("c_clock",     "Insertion -> methylation timing:\nNO CLOCK MEASURED (YY1 only proposed)", "C", "gap", "gap", "verified absence [V]", ["full"]),
    ("c_bulk",      "Bulk methylation -> insertion burden:\nNO GENERAL RULE ESTABLISHED", "C", "gap", "gap", "verified absence [V]", ["full"]),

    # ---------------- D. MEDICAL & DISEASE ----------------
    ("d_tissue",    "Tissue ORF1p IHC: 482/1,027 (47%)\nof cancers immunoreactive",  "D", "clinical", "demonstrated", "Rodic 2014 [V]", ["full"]),
    ("d_plasma",    "Plasma ORF1p (Simoa): 58% CRC (n=101),\n71% ovarian (n=145) @ 98% spec;\n~99% of 400+ controls undetectable", "D", "clinical", "associated", "Taylor 2023 [V]", ["full"]),
    ("d_multi",     "Ovarian multianalyte:\nCA125+HE4 alone AUC 0.94/59%;\n+ORF1p AUC 0.98/91% (the ORF1p gain)", "D", "clinical", "associated", "Taylor 2023 Fig 3D [V]", ["full"]),
    ("d_cfdna",     "cfDNA L1 hypomethylation: CRC\nAUC 0.81; stage I/II AUC 0.79,\n63.2% sens / 90% spec (114+53)", "D", "clinical", "associated", "Nagai 2017 [V]", ["full"]),
    ("d_alucf",     "Alu cfDNA meta (13 studies, 807 CRC):\n81% sens / 90% spec / AUC 0.92;\nprognosis HR 2.33 (1.80-3.03)", "D", "clinical", "associated", "Taghizadeh 2024 [V]", ["full"]),
    ("d_integrity", "Alu247/115 integrity: advanced\nlung AUC 0.88; CRC chemo-response\nchange AUC 0.952 (small single-center)", "D", "clinical", "associated", "Chudasama 2018; Eskander 2022 [V]", ["full"]),
    ("d_erv",       "9-ERV ccRCC signature (trial-linked):\nOS HR 0.49/0.41/0.55 across cohorts", "D", "clinical", "associated", "Zhou 2022 [V]", ["full"]),
    ("d_fda",       "FDA-CLEARED TE-BASED TEST:\nNONE EXISTS",                       "D", "gap", "gap", "EDRN list checked [V]", ["full"]),
    ("d_germline",  "Disease-causing germline MEIs:\n2012 tally 96 events (25 L1/60 Alu/\n7 SVA); 2017: 124 diseases; F8 exon 14:\n2/240 hemophilia A (first, 1988)", "D", "marker", "associated", "Hancks 2012 [S]; Kazazian 1988 [S]", ["full"]),
    ("d_xdp",       "SVA disease: X-linked dystonia-\nparkinsonism (TAF1 founder SVA;\nhexameric repeat expansion -> onset)", "D", "marker", "associated", "Bragg 2017 PNAS [S]", ["full"]),
    ("d_cancerdrv", "Pan-cancer somatic retrotransposition:\n19,166 events, 38 subtypes, 35% of\n2,954 genomes; CDKN2A 5.3-Mb del\nexample; APC ~1% CRC", "D", "marker", "demonstrated", "Rodriguez-Martin 2020 [V]; Cajuso 2019 [S]", ["full"]),
    ("d_neuron",    "Neuronal L1: somatic insertions REAL\nbut rate disputed (<0.6 vs 13.7 per\nneuron); no agreed per-neuron count", "D", "marker", "associated", "Evrony 2012 [V]; Upton 2015 [V]; Zhao 2019 [V]", ["full"]),

    # ---------------- E. RESOURCES ----------------
    ("e_counts",    "Genome content: L1 ~17% (~520k copies);\nAlu ~10-11% (1.1M; 1,238,897 in one\nhg38 analysis); SVA ~0.15% (5,107:\n1,927 FL + 3,180 trunc); pseudogenes ~7,819", "E", "marker", "demonstrated", "Han 2005; Lavi&Carmel 2018; Chu 2023; Zhang 2003 [V]", ["full"]),
    ("e_l1hs",      "L1Hs: 1,642 annotated (Chu); 321 >6kb;\nL1Base: 146 full-length INTACT L1s\n(13,418 >=4.5kb non-intact); 80-100\ncompetent; 6 hot = 84% of activity", "E", "marker", "demonstrated", "Chu 2023; Penzkofer L1Base 2 [V]; Brouha 2003 [V]", ["full"]),
    ("e_ages",      "Ages: AluJ ~65-40 Mya, AluS ~40-25,\nAluY youngest; L1PA1 ~3.1 Mya (2.3-3.9);\nrate: 1 de novo L1 per ~62.5 births\n(95% CI 31-154; prior range 20-200)", "E", "marker", "demonstrated", "Bailey 2003; Khan 2006 [V]; Feusier 2019 [V]", ["full"]),
    ("e_freq",      "Frequency spectrum: 70.7% of MEIs\nAF <0.1%; >30% singletons per type;\nAlu/L1/SVA spectra broadly SIMILAR", "E", "marker", "demonstrated", "Niu 2022 [V]; Stewart 2011; Gardner 2017", ["full"]),
    ("e_ta",        "L1 diagnostics: Ta vs pre-Ta = ACA/ACG\nat 3'UTR 5930-5932; Ta-0/Ta-1 = ORF2\n5536/5539; Beck 5'UTR bp-74 G-del", "E", "marker", "demonstrated", "Beck 2010 [V]; Shin 2018 [V]; Myers [S]", ["full"]),
    ("e_taours",    "TA ASSESSABILITY (ours): 58.4% of\n2,186 L1 calls reach the 3'UTR window\n5930-5932; 79.1% of L1HS+L1PA copies\ndo; ORF2 5536/5539 (Ta-0 vs Ta-1)\nreached by only 51.0% -> Ta level is\nOFTEN UNRESOLVABLE (not_assessable)\nand must print as its own denominator", "E", "marker", "demonstrated", "OUR span_work/out/line1_calls.tsv [W]", ["full"]),
    ("e_structs",   "Structures: ORF1p trimer 2YKP;\nORF2p+RNA 8UW3; ORF2p core 8SXT\n(tower/wrist domains); EN-DNA 7N94;\nTPRT complex 9HDP/9HDR; D145A, D702/3", "E", "structure", "demonstrated", "Khazina 2011; Thawani; Baldwin; Miller; Ghanim [V]", ["full"]),
    ("e_hmeid",     "HMEID: 36,699 MEI, 5,675 genomes\n(26,553 Alu / 7,353 L1 / 2,667 SVA)", "E", "resource", "demonstrated", "bigdata.ibp.ac.cn/HMEID", ["full"]),
    ("e_1kgp",      "1000G MEI: >20,000 polymorphic\n(pilot 7,380); DGVa nstd144",       "E", "resource", "demonstrated", "IGSR / DGVa", ["full"]),
    ("e_hgsvc",     "HGSVC: 9,453 resolved MEIs\n(7,738 Alu / 1,175 L1 / 540 SVA), 64 haplotypes", "E", "resource", "demonstrated", "Ebert 2021; hgsvc.org", ["full"]),
    ("e_lr1019",    "1,019-genome long-read: 31,302 MEI\n(23,212 Alu / 4,851 L1 / 3,239 SVA);\n96.4% hallmark-bearing", "E", "resource", "demonstrated", "Schloissnig 2025 [V]", ["full"]),
    ("e_logsdon",   "Logsdon 2025: 12,919 MEIs vs T2T-CHM13\n(65 individuals, 130 haplotypes, 39% T2T)", "E", "resource", "demonstrated", "Logsdon 2025 [V]", ["full"]),
    ("e_gnomad",    "gnomAD-SV: MEI subset of 1.19M SVs\n(63,046 genomes); DGVa nstd166", "E", "resource", "associated", "gnomAD", ["full"]),
    ("e_nesteddb",  "NESTED-MEI DATABASE\n(host-relative coordinates):",       "E", "gap", "gap", "DOES NOT EXIST ANYWHERE [V]", ["full"]),

    # ---------------- G. OUR DATA ----------------
    ("g_hg",        "HG03086 (short-read): 1,852 MEI calls\n- Q1: same-family nesting 1.47x/1.47x/29.4x\n- Q2: sense 2.49x, antisense 0.46x (Alu)\n- Q3: homology gradient 1.86x/0.88x/0.72x\n- Q4: linker 3.0x + tail 3.7x peaks", "G", "marker", "demonstrated", "selfins_*.csv [W]", ["full"]),
    ("g_lr",        "Long-read SVAN cohort (pooled 908):\n~2,586 nested-Alu hosts, ~400 tail-peak\nTSD/polyA raw fields; TSD_LEN=41\nsentinel GATE unresolved; HG03086 in-cohort", "G", "marker", "demonstrated", "svim.asm.hg38.bcf [W]", ["full"]),
    ("g_meth",      "l1meth: chr22 locus methylation\n11,526 repeats, 43 addressable;\nage_enrichment default corrected 8x -> 1x\n(finding 4 falsified by Lanciano)", "G", "marker", "demonstrated", "l1_locus_catalog_chr22.csv [W]", ["full"]),
    ("g_plan",      "Plan v2 (approved): Phase 0 gates ->\nPhase 1 pos-133 (3 nulls) -> Phase 2a\nboundary-perturbation falsifier ->\nPhase 2b triad (gated) -> Phase 3 concordance (contingent)", "G", "marker", "demonstrated", "rtm-nested-analysis plan [W]", ["full"]),
    ("g_william",   "Incoming: William processing extra\ngenomes (~2 h each); power counted in\nUNIQUE SITES (dedup pre-registered:\nchrom+family+orient+host, +-10 bp)", "G", "marker", "demonstrated", "agreed 2026-09 session [W]", ["full"]),
    ("g_tsdage",    "TSD x AGE, NON-CIRCULAR (ours): crude\nL1HS .787 / L1P* .538 / L1M* .254\nCOLLAPSES to .761 / .726 / .636 once\n3' resolution is held constant; resolved\nL1HS vs L1P* OR 1.20, Fisher p=0.18\n(>=5bp) -> the lineage trend does NOT\nsurvive, so TSD presence cannot assign\nage. A span-independent axis does\nsurvive: catalog membership MH OR 6.38\n(KNOWNMEI) / 4.06 (1KG) / 3.62 (long-\nread), direction kept in BOTH strata,\nbut it is confounded with evidence\navailability (SUPPORT empty in file)", "G", "marker", "demonstrated", "artifacts/tsd_age_association.py [W]", ["full"]),

    # ---------------- X. VERIFIED GAPS ----------------
    ("x_tsdcomp",   "TSD base-composition per position,\nper family: NEVER PUBLISHED\n(our TSD_SEQ could be first ever)", "X", "gap", "gap", "verified absence [V]", ["full"]),
    ("x_tsdagedecay", "TSD DECAY WITH ELEMENT AGE:\nNO per-age-class presence rate and NO\ndecay curve published. Szak 2002's\n'more recognizable in younger elements'\nrests on the 88% (+-7) mean identity of\na TSD-POSITIVE set; Zingler 2005's 90%\nis composition within 10,034 TSD-positive\ncopies. NEITHER has an all-copies\ndenominator per age class, so the premise\nmust not be cited as a measured effect", "X", "gap", "gap", "verified absence [V]", ["full"]),
    ("x_tss",       "Distance-to-TSS distribution\nfor germline insertions: NOT PUBLISHED", "X", "gap", "gap", "verified absence [V]", ["full"]),
    ("x_intronhalf","5' vs 3' intron-half preference:\nNOT QUANTIFIED anywhere",             "X", "gap", "gap", "verified absence [V]", ["full"]),
    ("x_splice",    "Splice-site density (per kb,\nsubfamily x orientation): NOT PUBLISHED\n(Sorek/Haesler = consensus-level only)", "X", "gap", "gap", "verified absence [V]", ["full"]),
    ("x_aluinv",    "Alu 5'-inversion frequency:\nNOT PUBLISHED (Wagstaff silent;\nZingler excluded inversions)", "X", "gap", "gap", "verified absence [V]", ["full"]),
    ("x_svatsd",    "SVA TSD histogram: NOT PUBLISHED\n(4-20 bp range only)",         "X", "gap", "gap", "verified absence [V]", ["full"]),
    ("x_fda",       "FDA-cleared TE-based test: NONE",                               "X", "gap", "gap", "EDRN checked [V]", ["full"]),
    ("x_trial",     "Dedicated L1/TE-antigen clinical trial:\nNONE VERIFIED 2023-2025\n(lamivudine CRC = off-label, pre-2023)", "X", "gap", "gap", "verified absence [S]", ["full"]),
]

# =====================================================================================
#  EDGES   (src, dst, tier, label, views)
# =====================================================================================
EDGES = [
    # ----- spine -----
    ("meth", "rna", "demonstrated", "permissive", ["full", "spine"]),
    ("rna", "orf2p", "demonstrated", "", ["full", "spine"]),
    ("orf2p", "en", "demonstrated", "", ["full", "spine"]),
    ("en", "tprt", "demonstrated", "", ["full", "spine"]),
    ("tprt", "position", "unknown", "", ["full", "spine"]),
    ("tprt", "junction", "demonstrated", "", ["full", "spine"]),
    ("junction", "lineage", "demonstrated", "", ["full", "spine"]),
    ("lineage", "reseal", "proposed", "", ["full", "spine"]),
    ("reseal", "clinical", "associated", "", ["full", "spine"]),
    ("meth", "clinical", "associated", "state readout", ["full", "spine"]),

    # ----- A: targeting -> spine anchors -----
    ("a_motif7", "en", "demonstrated", "the motif", ["full"]),
    ("a_motifC", "a_motif7", "demonstrated", "", ["full"]),
    ("a_motifopp", "a_motif7", "demonstrated", "opportunity", ["full"]),
    ("a_nick2", "tprt", "demonstrated", "", ["full"]),
    ("a_polytsd", "tprt", "demonstrated", "", ["full"]),
    ("a_ssDNA", "tprt", "demonstrated", "architecture", ["full"]),
    ("a_dnashape", "en", "demonstrated", "DNA shape", ["full"]),
    ("a_wedge", "tprt", "demonstrated", "", ["full"]),
    ("a_stag", "junction", "demonstrated", "TSD variability", ["full"]),
    ("a_ends", "en", "associated", "", ["full"]),
    ("a_gc50", "en", "demonstrated", "local", ["full"]),
    ("a_gc20kb", "a_gc50", "demonstrated", "contrast", ["full"]),
    ("a_alu1mb", "a_gc50", "demonstrated", "", ["full"]),
    ("a_svagc", "en", "associated", "family flip", ["full"]),
    ("a_genic_alu", "position", "associated", "", ["full"]),
    ("a_genic_l1", "position", "associated", "null-controlled", ["full"]),
    ("a_genic_sva", "position", "associated", "", ["full"]),
    ("a_exonprox", "position", "demonstrated", "selection", ["full"]),
    ("a_intronic", "position", "demonstrated", "", ["full"]),
    ("a_lavialu", "position", "demonstrated", "positional", ["full"]),
    ("a_levmaor", "a_lavialu", "demonstrated", "mechanism", ["full"]),
    ("a_133", "position", "demonstrated", "PUBLISHED", ["full"]),
    ("a_133ours", "a_133", "demonstrated", "REPLICATED [W]", ["full"]),
    ("a_tail", "position", "unknown", "NOVEL [W]", ["full"]),
    ("a_5p", "position", "unknown", "", ["full"]),
    ("a_sense", "position", "demonstrated", "", ["full"]),
    ("a_antis", "a_sense", "demonstrated", "family split", ["full"]),
    ("a_homol", "a_sense", "demonstrated", "does NOT explain", ["full"]),
    ("a_ccatt", "a_sense", "demonstrated", "route", ["full"]),
    ("a_repl", "en", "demonstrated", "context", ["full"]),
    ("a_fork", "en", "associated", "", ["full"]),
    ("a_chrom", "en", "demonstrated", "no pref", ["full"]),
    ("a_nuc", "a_chrom", "associated", "", ["full"]),
    ("a_consensus", "a_133", "demonstrated", "landmarks", ["full"]),

    # ----- B: credibility -> junction/orf2p -----
    ("b_triad", "junction", "demonstrated", "yardstick", ["full"]),
    ("b_tsd_alu", "b_triad", "demonstrated", "", ["full"]),
    ("b_tsd_l1", "b_triad", "demonstrated", "", ["full"]),
    ("b_tsd_l1fl", "b_tsd_l1", "demonstrated", "", ["full"]),
    ("b_tsd_sva", "b_triad", "demonstrated", "", ["full"]),
    ("b_tsd_cat", "b_triad", "demonstrated", "catalog rates", ["full"]),
    ("b_tsd_lr", "b_triad", "demonstrated", "long-read", ["full"]),
    ("b_tsdshort", "b_tsd_l1", "demonstrated", "caution", ["full"]),
    ("b_micro", "junction", "demonstrated", "L1-route tag", ["full"]),
    ("b_inv", "junction", "demonstrated", "L1-route tag", ["full"]),
    ("b_invporub", "b_inv", "demonstrated", "", ["full"]),
    ("b_ip", "junction", "demonstrated", "alt route", ["full"]),
    ("b_del", "junction", "demonstrated", "", ["full"]),
    ("b_filler", "junction", "demonstrated", "", ["full"]),
    ("b_rearr", "junction", "demonstrated", "complex", ["full"]),
    ("b_trunc", "junction", "demonstrated", "family asymmetry", ["full"]),
    ("b_trans_l1", "lineage", "demonstrated", "parent tag", ["full"]),
    ("b_trans_sva", "lineage", "demonstrated", "parent tag", ["full"]),
    ("b_recomb", "junction", "demonstrated", "aftermath", ["full"]),
    ("b_source", "orf2p", "demonstrated", "6 hot = 84%", ["full"]),
    ("b_callers", "junction", "demonstrated", "read evidence", ["full"]),
    ("b_niche", "b_callers", "gap", "unclaimed", ["full"]),

    # ----- C: state -> meth/rna/orf2p -----
    ("c_chain1", "meth", "demonstrated", "link 1", ["full"]),
    ("c_chain2", "rna", "demonstrated", "link 2-3", ["full"]),
    ("c_tubio", "meth", "associated", "source-level", ["full"]),
    ("c_freeman", "meth", "demonstrated", "", ["full"]),
    ("c_byun", "meth", "demonstrated", "young-Alu", ["full"]),
    ("c_setdb1", "meth", "demonstrated", "causal", ["full"]),
    ("c_youngctx", "meth", "demonstrated", "context map", ["full"]),
    ("c_pgc", "meth", "demonstrated", "reset", ["full"]),
    ("c_clock", "reseal", "gap", "empty", ["full"]),
    ("c_bulk", "meth", "gap", "empty", ["full"]),

    # ----- D: medical -> clinical -----
    ("d_tissue", "clinical", "demonstrated", "47%", ["full"]),
    ("d_plasma", "clinical", "associated", "", ["full"]),
    ("d_multi", "d_plasma", "associated", "the ORF1p gain", ["full"]),
    ("d_cfdna", "clinical", "associated", "", ["full"]),
    ("d_alucf", "clinical", "associated", "", ["full"]),
    ("d_integrity", "clinical", "associated", "", ["full"]),
    ("d_erv", "clinical", "associated", "", ["full"]),
    ("d_fda", "clinical", "gap", "empty", ["full"]),
    ("d_germline", "clinical", "associated", "germline disease", ["full"]),
    ("d_xdp", "d_germline", "associated", "SVA", ["full"]),
    ("d_cancerdrv", "clinical", "demonstrated", "somatic", ["full"]),
    ("d_neuron", "clinical", "associated", "disputed rate", ["full"]),

    # ----- E: resources -> spine -----
    ("e_counts", "position", "demonstrated", "scale", ["full"]),
    ("e_l1hs", "orf2p", "demonstrated", "sources", ["full"]),
    ("e_ages", "position", "demonstrated", "timescale", ["full"]),
    ("e_freq", "e_counts", "demonstrated", "rare-skewed", ["full"]),
    ("e_ta", "e_l1hs", "demonstrated", "lineage marks", ["full"]),
    ("e_taours", "e_ta", "demonstrated", "assessable fraction [W]", ["full"]),
    ("e_structs", "en", "demonstrated", "mechanism", ["full"]),
    ("e_hmeid", "position", "associated", "known sites", ["full"]),
    ("e_1kgp", "position", "associated", "known sites", ["full"]),
    ("e_hgsvc", "junction", "associated", "resolved", ["full"]),
    ("e_lr1019", "junction", "associated", "hallmarks", ["full"]),
    ("e_logsdon", "position", "associated", "T2T", ["full"]),
    ("e_gnomad", "position", "associated", "freq", ["full"]),
    ("e_nesteddb", "position", "gap", "missing", ["full"]),

    # ----- G: our data -> spine -----
    ("g_hg", "position", "demonstrated", "measured [W]", ["full"]),
    ("g_lr", "junction", "demonstrated", "measured [W]", ["full"]),
    ("g_meth", "meth", "demonstrated", "measured [W]", ["full"]),
    ("g_plan", "g_lr", "demonstrated", "pre-registered", ["full"]),
    ("g_william", "g_plan", "demonstrated", "incoming", ["full"]),
    ("g_tsdage", "b_tsd_l1", "demonstrated", "age confound resolved [W]", ["full"]),
    ("g_tsdage", "junction", "demonstrated", "measured [W]", ["full"]),
    ("g_tsdage", "e_taours", "demonstrated", "same cohort [W]", ["full"]),

    # ----- X: gaps -> their spine homes -----
    ("x_tsdcomp", "junction", "gap", "EMPTY", ["full"]),
    ("x_tsdagedecay", "junction", "gap", "EMPTY", ["full"]),
    ("x_tss", "position", "gap", "EMPTY", ["full"]),
    ("x_intronhalf", "position", "gap", "EMPTY", ["full"]),
    ("x_splice", "junction", "gap", "EMPTY", ["full"]),
    ("x_aluinv", "junction", "gap", "EMPTY", ["full"]),
    ("x_svatsd", "junction", "gap", "EMPTY", ["full"]),
    ("x_fda", "clinical", "gap", "EMPTY", ["full"]),
    ("x_trial", "clinical", "gap", "EMPTY", ["full"]),
]
