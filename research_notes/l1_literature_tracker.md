# LINE-1 literature tracker

Running log of what's been read/found and what's still open. Update as you go —
this is meant to accumulate, not get rewritten each time.

## Status of PRs (for context, not science)
- `feat/gene-annotation` -> PR #51, open, awaiting William's review. No further
  coding until he responds.
- `feat/ec2-ami-validation` -> in progress (Kilo Code, background).

---

## Open question 1: Is locus-specific L1 activation causal for chromatin state, or a consequence?

**Your answer so far: both.** Supported, but with a scope caveat — read critically.

- Chromatin -> L1 (the direction you already knew): DNA methylation / H3K9me3 /
  KRAB-ZFP repression silences L1; loss of that repression (aging, cancer
  hypomethylation) de-represses it. Well established, not the interesting half.
- L1 -> chromatin (causal, the newer half): LINE1 RNA itself (the transcript,
  not the DNA copy) organizes chromatin.
  - Percharde et al. 2018, *Cell* — LINE1 RNA partners with Nucleolin;
    required for ESC self-renewal and pre-implantation progression.
  - Jachowicz et al. 2017, *Nat Genet* — LINE-1 activation after fertilization
    regulates global chromatin accessibility in the early mouse embryo.
  - Wei et al. 2022, *Science* (FTO-LINE1 axis) — LINE1 RNA promotes open
    chromatin at LINE1-containing loci via histone-modifier recruitment;
    removing the RNA closes chromatin and installs repressive marks. Cleanest
    causal (RNA -> chromatin) mechanism found so far.
  - Zhang et al. 2024, *Dev Cell* — LINE1 RNA + PRC2 maintain nucleolar
    architecture / 3D genome organization in human ESCs.

**Caveat (unresolved): all four papers above are ESC / early-embryo systems.**
None yet shown in cancer or aging tissue, which is the context the biomarker
pitch actually needs. Next read: does any paper make this causal claim outside
stem cells/embryos? If not, "both" is true in principle but not yet shown where
you want to apply it.

---

## Open question 2: L1 gets altered post-expression — are we (the tool) missing this?

**Yes, and it's real** — but it lives in a different data type than what
retrotransposon-miner reads (DNA/WGS). Can't be added as a quick feature;
would need RNA-seq or a functional assay, not an insertion caller.

Mechanisms matching your three guesses, with the actual literature:
- "RNAi-like" -> piRNA pathway (Piwi proteins) degrades L1 transcripts +
  methylates L1 promoters de novo. **Germline-specific** — don't generalize to
  somatic/cancer/aging contexts.
- "NMD" -> MOV10 (RNA helicase) forms a complex with UPF1 (core NMD factor),
  reduces L1 RNA half-life, restricts retrotransposition. Guan et al. 2023 —
  loss of MOV10 increases L1 retrotransposition in **somatic** tissue too.
  This one is actually relevant to disease contexts.
- "Antiviral proteins" -> ZAP (ZC3HAV1) and APOBEC3 both restrict L1 RNA/cDNA;
  co-localize with L1 RNPs in stress granules along with G3BP1, TIA1, PABP.
  RNase L also degrades L1 RNA directly.

**Open question to answer before proposing a feature**: what assay would
actually show MOV10/ZAP/APOBEC3 activity on a specific L1 copy? (RNA-seq for
transcript abundance/half-life; proteomics or IP for protein binding.) Until
that's answered, this is a "different tool" idea, not a "this tool is missing
a column" idea.

---

## Open question 3: Is "L1 as a biomarker" novel?

**No — 15+ year old field, and results are mixed, not a clean win.**

- LINE-1 makes up ~17% of the genome; bulk methylation level is used as a
  surrogate for genome-wide hypomethylation. Hypomethylation linked to worse
  prognosis across GI, breast, lung, liver, prostate, oropharyngeal cancers
  (various single-cohort studies, mostly Ogino-lab lineage).
- But: 2024 meta-analysis, 19 studies / 8169 CRC patients -> **no significant
  pooled association** with overall survival, despite earlier single-cohort
  enthusiasm. A ~800-patient prostate cohort found no association with
  high-risk disease or recurrence either. A dedicated review states direct
  evidence L1 activation can be used as a cancer biomarker is still limited.

**Critical distinction**: all of the above measures ONE genome-averaged
number across ~500,000 L1 copies (bisulfite consensus assay) — not any
individual locus. That's a different, much cruder measurement than what
retrotransposon-miner does (individual, locus-specific insertion calls).

**New lead, UNVERIFIED, from a side-chat search — check before citing
anywhere**: claim that the genome-averaged methylation number is a weak
discriminator in a large CRC cohort (~869 tumors, R^2=0.084, AUC<0.63), and
that "target the CpG sites specifically" splits into two different unmet
needs: (a) locus-specific methylation detection, which needs bisulfite or
long-read methylation calling — a different pipeline than the alignment-based
insertion calling the tool does today — versus (b) actual epigenetic editing
(dCas9-DNMT3A etc.), which is still research-stage. **Action: find and read
the actual 869-tumor study before repeating these numbers to anyone.**

---

## Where this leaves the "novel capability" pitch

The genuinely open door isn't "L1 as a biomarker" (crowded, mixed results) —
it's specifically **locus-specific** methylation or activation state, which
almost nobody has measured because it needs a different assay than the bulk
consensus method everyone else uses. That's also exactly the assay gap from
Question 2. Both open questions point at the same missing thing: something
that reads L1 at single-locus resolution, at the RNA or methylation level,
not the DNA-insertion level. Worth stating explicitly next time you talk to
William, once you've done the reading to back it.

---

## Open question 4: L1 5'UTR evolution — new subfamilies, and can we detect them?

**Your idea, and it holds up.** The 5'UTR is simultaneously (a) the internal
RNA Pol II promoter, (b) the CpG island that methylation acts on, and (c) the
fastest-evolving part of the element. Those aren't three facts — they're one
piece of DNA doing three jobs, which is why "target the CpG sites" and
"silence the promoter" are the same act.

Verified:
- Internal Pol II promoter lives in the 5'UTR; required for autonomous L1
  transcription and for initiating retrotransposition.
- 5'UTR contains sense promoter, antisense promoter, TF binding sites, CpG
  sites and CpG islands, all overlapping (Khan/Lee analyses of 443 full-length
  human-specific L1s).
- The L1 promoter has a HIGHER substitution frequency than the L1 average —
  it is the fastest-changing regulatory region in the element.
- New subfamilies arise by 5'UTR REPLACEMENT: >=8 episodes in ~70 My of
  primate evolution; most recent ~40 MYA giving rise to L1PA8, whose
  architecture persists through L1PA7 -> L1PA1.
- Mechanism/why: new 5'UTRs let emergent subfamilies escape host suppression.

**The canonical arms-race example (memorize this one):** L1PA6-L1PA3 evolved a
ZNF93 binding motif in the 5'UTR -> recruits ZNF93 -> KAP1-mediated silencing.
A 129-bp DELETION in the 5'UTR removing that motif let a subset of
L1PA3/L1PA2/L1PA1 escape ZNF93. Host evolves a repressor; L1 deletes the
binding site. The escape event is a visible sequence change.

**Animal L1 (your "is animal L1 evolving?"):** yes, and more so than human.
Human = ONE contemporarily active subfamily (L1HS). Mouse = FOUR active in the
last 1 My (A_I, Tf_I, Tf_II, Gf_I). Mouse 5'UTRs are structurally different:
tandem repeated "monomers" plus a tether domain, not a single 5'UTR. YY1
activates the mouse Tf subfamily via monomer binding sites (but not Gf).

### Tool finding (checked in code, 2026-09-19)

retrotransposon-miner **cannot detect a novel subfamily, by construction**:
- Subfamily = best alignment to a fixed **Dfam panel** of known consensus
  sequences (L1HS_5end, L1HS_3end, orf2, ...).
- There is NO "matches nothing well" / unclassified category. A novel 5'UTR is
  forced into the nearest known bucket and reported with confidence.
- `mei_support.py:299` strips `_5end`/`_3end`/`_orf2` suffixes, collapsing the
  5'-end-specific signal into a bare family name before it reaches output.

**Candidate capability (NOT started — PR #51 must land first):** an
unmatched-5'UTR / novel-subfamily flag. Report per-insertion 5'-end alignment
identity+coverage against the panel, and flag insertions whose 5' end aligns
poorly to every known member. Would also surface 5'UTR deletions of the
ZNF93-escape type. Cheap to compute — the 5end panel alignments already exist,
they're just discarded.

---

## Corrections logged (ideas that didn't survive contact)

- "Hypomethylation **caused by** L1" — inverted. L1 copies are the SUBSTRATE
  where genome-wide hypomethylation is measured, not its cause. Real L1->
  methylation direction exists but with opposite sign: an L1 insertion can
  nucleate silencing of a neighbouring gene (local HYPERmethylation).
  Confidence moderate — read the L1-mediated gene silencing literature before
  citing this.
- "Locate exact area AND disease from L1" — first half is the genuine gap
  (bulk bisulfite averages ~500k copies into one number; you cannot localize
  an average). Second half is a much bigger claim: inferring WHICH disease
  needs labeled cohorts and a trained classifier, and that is exactly where
  the existing biomarker literature keeps failing to replicate.

---

## Open question 5: radiation, mouse strains, and environmental effects on L1 methylation

**Terminology conflation to fix first:** "cell lineage," "mouse strain," and "L1
subfamily" are three independent axes, not synonyms. The paper below uses all
three at once, which is probably where the mixing happened.

### Primary source: Miousse et al. 2017, Int J Mol Sci 18:1430
DOI 10.3390/ijms18071430 (full text fetched and read, 2026-09-19)

Design: 5'UTR LINE-1 methylation in HSC -> HPC -> MNC (hematopoietic CELL
LINEAGE) in two mouse STRAINs (C57BL/6J radioresistant vs CBA/J radiosensitive/
cancer-prone), basal and after ionizing radiation (0.1, 1 Gy, 137Cs).

Findings, precise:
- Basal methylation anti-correlates with SUBFAMILY age: elements <1 Myr old
  carry ~3x the methylation of elements >1 Myr old. Holds in every cell type
  tested -> age effect, not a lineage-stage effect.
- Radiation effect was STRAIN-specific: only CBA/J (radiosensitive) lost HSC
  methylation post-IR; C57BL/6J did not respond at all.
- Ruled out the obvious explanation: Dnmt1, Dnmt3a, Uhrf1 (the actual
  methylation-maintenance enzymes) did NOT differ between strains. Whatever
  drives CBA/J's sensitivity, it isn't "fewer methylation enzymes."
- SUBFAMILY response was non-monotonic across strain x cell-lineage-stage:
  L1MdA_I (youngest) 2.2x HIGHER in C57BL/6J than CBA/J in HSCs; but
  L1MdTf_I and L1MdFanc_I ("ancient") 4.1x/3.6x HIGHER in CBA/J than C57BL/6J
  in HPCs, and 18.3x/4.6x(trend) higher in MNCs. Same two strains, opposite
  direction, depending on subfamily AND lineage stage. No single number
  captures this -- per-subfamily measurement is necessary, not just nice to
  have, to see the effect at all.
- At 2 months post-IR (when hypomethylation was observed in CBA/J HSCs): NO
  residual DNA damage, NO elevated ROS, NO cell-cycle change detectable in
  those same cells. The acute damage signals were gone; the methylation
  change was not. Mechanism connecting the two is UNRESOLVED.
- Authors' own hedge, explicitly stated: epigenetic change is presented as a
  *potential* driving force of radiation carcinogenesis; direct causal link
  not demonstrated. SAME causality-direction trap as Idea 1 (chromatin
  thread) -- this time the paper's own authors flag it, not just us.

### Broader literature (multiple sources, web search 2026-09-19)
- Type AND dose of radiation both matter (gamma, X-ray, protons, 56Fe heavy
  ions all tested, effects differ).
- Effect is BIPHASIC, not one-directional: mouse heart showed L1
  hypomethylation at day 7 post-IR, flipping to hypermethylation (of L1 +
  ERV2 + SINE B1 + satellite DNA) by day 90. Direction depends on timepoint.
- Methylation and expression DECOUPLE over time: one rat model showed
  methylation recovered to baseline by 12-24 wks but LINE-1 protein stayed
  elevated regardless. The mark and its consequence stopped moving together.
- "More prevalent" (your phrase) is actually two different measurable things,
  both documented separately: (a) more transcription/translation of EXISTING
  copies (protein level increase), vs (b) more NEW copies via retrotransposition
  (one culture study: up to 4-fold increase under gamma IR, scaling with
  gammaH2AX, the standard double-strand-break marker). Decide which you mean;
  they need different assays.
- REAL HUMAN data exists: industrial radiographers (occupational radiation
  exposure) show LINE-1 hypomethylation associated with radiation-induced
  genomic instability. Cleaner exposure variable (measured cumulative dose)
  than the noisy cancer-outcome studies in Q3's biomarker literature.
- UNREAD lead: paper titled on "densely ionizing radiation" affecting
  "selective" LINE-1 elements (Environ Res 150:470) -- title suggests
  element-specific response, consistent with Miousse, not yet read/verified.

### Open mechanism question (not yet answered by field or by us)
By the time CBA/J mice showed hypomethylation (2 months post-IR), every
marker of the original radiation damage was gone. What connects "radiation
happened" to "a methylation change that outlasts every marker of the
radiation itself"? Two hypotheses, different experiments needed:
 (1) change installed during the acute damage window, stable for months after
     (hit-and-run: only ever see the aftermath in a single 2-month snapshot)
 (2) change doesn't start until well after damage clears -> not damage-driven
     at all, something slower/indirect (population shift? delayed stress
     pathway?)
A single 2-month timepoint can't distinguish these -- would need a time
course. Posed to [user] 2026-09-19, answer pending.

### Connection to the subfamily-divergence pitch already sent to William
NOT sending a follow-up yet -- he hasn't responded to PR #51 or the first
pitch. But logging the supporting data point: the L1MdA_I / L1MdTf_I /
L1MdFanc_I strain x lineage-stage divergence above is a concrete case where
subfamily identity determines the direction of a biological response, not
just its magnitude. That's a stronger example than anything in the pitch as
sent. Good ammunition for a follow-up conversation once he's replied once.

### Follow-up: timing question answered, verified against the actual data (2026-09-19)

[User]'s answer: hypothesis 1 (installed during the acute damage window, stable
after). Checked against the paper's actual section structure before accepting
or rejecting this.

**Verified fact:** Miousse et al. 2017 sections 2.4/2.5/2.6 (all irradiated-mouse
methylation data) are ALL titled "Long-Term Effects..." There is no
intermediate timepoint between baseline and 2 months post-IR anywhere in this
paper. The study CANNOT distinguish "installed early, stable for 2 months"
from "installed late, near the 2-month mark" -- the data to tell them apart
doesn't exist in this dataset.

**Reasoning errors flagged (not just the conclusion -- the argument for it):**
1. "Methylation is the hallmark of everything, so it cannot happen afterwards"
   -- category error. Prevalence (how common a feature is) says nothing about
   timing (when it arose). Named explicitly so it doesn't recur.
2. "If it happened afterwards there would be differences" -- incomplete as
   stated; doesn't specify differences in what, expected when. Actually
   describes the experiment needed (a time-course), not evidence already in
   hand.
3. Internal tension: [user]'s own proposed chain (radiation -> genomic
   instability -> hypomethylation, "hypomethylation was in response to...
   genomic instability") makes hypomethylation DOWNSTREAM of instability,
   which undercuts "happens during [the radiation exposure]" as stated --
   "during vs after" isn't well-defined once there are 3 links without saying
   during/after which transition.
4. Independent evidence against "that's just how it is" as a general rule:
   the mouse HEART study (already in tracker, Q5) measured 2 timepoints (day
   7, day 90) and found the direction FLIPPED. The one case in this
   literature with >2 timepoints checked shows time-dependence. Reason to
   expect the HSC study might too, if anyone had looked.

**Conclusion: genuinely unresolved, not "probably hypothesis 2" either.**
Correct next step (not yet done): a time-course design, not a 2-point design.

**Open sub-question posed back:** "environmental and genetic point of view"
already IS this study's design (strain=genetic, radiation=environmental).
Asked [user] to name ONE specific variable on either axis not already tested
(age at exposure? sex? chronic low-dose vs the acute doses used? other
environmental stressor?) -- answer pending.

---

## Open question 6: is the divergence pitch actually novel? (prior-art check, 2026-09-19)

Ran this AFTER the pitch was already sent to William. Result: the sent message
survives, but my framing of it to [user] was overclaimed. Logging both.

### CORRECTION to my own earlier claim
I told [user]: "reporting divergence per called insertion is not what MEI
callers do -- they emit labels." **That was wrong.**
- **MELT emits a DIFF field** = per-insertion differences from consensus.
  MELT author: the only way to get consensus sequences is to write a script
  to reconstruct them from the MEI and DIFF fields.
- **MELT has an ambiguity bucket** -- users reference "L1Ambig" for
  insertions whose subfamily can't be resolved. So a "matches-nothing"
  category is not unheard-of in the field.

### Why the SENT message is still defensible
[User]'s draft claimed *retrotransposon-miner* has no matches-nothing category
and forces novel subfamilies into the nearest bucket. That's a claim about
THIS tool, verified from its source, and it stands. The draft did not claim
the field lacks this. No correction needed to William.

### Likely William question to be ready for: "how is this different from MELT's DIFF?"
What survives as genuinely distinct:
1. MELT's DIFF is raw reconstruction data -- you script your own analysis from
   it. Not a score, not ranked, not phylogenetically anchored.
2. xTea uses divergence as a FILTER, not an output: candidates in a
   same-family repetitive region are REMOVED if the reference repeat's
   divergence from RepeatMasker consensus is below a threshold. Note this is
   the divergence of the *reference element at the locus* (a mappability
   guard) -- NOT the divergence of the *inserted* element. Different
   measurement.
3. **5'UTR-restricted divergence: no prior art found.** Whole-element
   divergence is diluted by ORF2 (conserved). Still the strongest novelty
   claim.
4. Subfamily -> suppressor-evasion -> functional consequence linkage: still
   nothing found.

### Point that STRENGTHENS the pitch
Published review names the gap directly: "a major limitation of the existing
de novo TE insertion callers is that they still only focus on detecting
insertions, like other SV callers, but do not provide TE-specific annotation,
such as TE family or target site duplication." retrotransposon-miner ALREADY
emits family + TSD -- so it is ahead of that critique, and divergence is the
natural next field in exactly the direction the review says tools are
lacking. Good framing for a reply to William.

---

## Scaffolding dispatched to Kilo (2026-09-19)

Decision: build the **validation harness, not the feature**, on a branch off
origin/main while PR #51 is in review. Rationale -- if William says yes the
validation is ready; if he says "do it differently" the harness still works
because it measures rather than implements; and it's small enough not to
swamp review. Explicitly NOT endorsing "build it and he'll approve anyway" --
volume of finished code makes reviewers less likely to engage, not more.

### Key enabling finding (verified in source)
The repo ALREADY ingests the UCSC RepeatMasker rmsk table via the existing
`--rmsk-table` option (used today only for nested-insertion annotation).
`_write_rmsk_mei_bed` (mei_support.py ~13687) splits each line and reads
parts[5,6,7,9,10,11,12] -- **walking straight past parts[2] = milliDiv**,
RepeatMasker's divergence-from-consensus in units of 1/1000. `grep -rn
"milliDiv" src/` returns NOTHING. The field's standard divergence measure,
and our ground truth, is already in a file the repo already opens, one column
index away from code that already runs.

Harness scope: parse milliDiv -> pure `panel_divergence_from_nm(nm, aln_len)`
helper -> `scripts/validate_divergence_vs_rmsk.py` correlating NM-derived
divergence against milliDiv/1000 -> tests on synthetic fixtures. <400 lines.
Branch `feat/divergence-harness` off origin/main, separate worktree.
Prompt saved: notes/kilo_prompt_divergence_harness.md

### Repo hazard noted
Main checkout is in DETACHED HEAD with uncommitted modifications including
DELETIONS (callset_eval.py, compare_sample_mei_callsets.py) -- the other Kilo
session mid-work. Prompt instructs the new agent to never touch that dir.

---

## Open question 7: 5-hydroxymethylcytosine (5hmC) / oxBS-seq as a new axis (2026-09-19)

Source: pasted attachment from [user] (provenance asked, not yet confirmed --
reads like a side-chat or parallel-session synthesis, already using
verified/inference/guess framing). Independently checked every load-bearing
claim rather than trusting the self-applied "verified" labels.

### Confirmed exactly as stated
- oxBS-seq mechanism: oxidize 5hmC -> 5fC -> bisulfite converts 5fC to
  uracil; only true 5mC still reads as C; 5hmC inferred by comparing to a
  parallel plain-BS-seq run. Matches the method's own paper essentially
  verbatim.
- 5hmC generated from 5mC by TET family oxidation.
- 5hmC is a stable mark, not just a transient demethylation intermediate.

### Confirmed with a correction
- Brain enrichment: sources consistently say ~10-FOLD (vs ESCs, vs other
  neurons/cell types), tighter than the pasted "5 to 10x" range. Use ~10-fold
  going forward.
- Alzheimer's and Huntington's disease links: solidly confirmed (HD: "genome-
  wide loss of 5-hmC is a novel epigenetic feature of Huntington's disease,"
  Wang et al 2013 Hum Mol Genet; AD: multiple sources incl. a 1079-brain
  genome-wide study).
- Cerebral ischemia and FXTAS links: NOT independently confirmed by this
  check (didn't surface in search). Not claiming they're wrong -- just
  unconfirmed by me. Don't cite as settled until sourced.
- "5hmC persists for months without turnover" -- general stable-mark concept
  confirmed, the specific "months" quantifier not pinned to a source.

### Correction to something I (assistant) told [user] previously
Pasted text said "I haven't found a paper applying oxBS-seq specifically to
LINE-1 loci" and flagged that as an open gap. **Found one:** the FOUNDING
oxBS-seq paper itself (Booth et al. 2012, Science -- the paper that invented
the method) reports high 5hmC in CGIs and in "long interspersed nuclear
elements" in mouse ES cells. Caveats: mouse ESCs not disease tissue; "LINEs"
as a class (L1-dominated but not stated L1-specific); almost certainly
genome-wide/regional resolution, not individual-copy. So: the specific
combination has precedent from the very first method paper, but the
individual-locus multi-mapping problem below is probably still unsolved by
it.

### The problem 5hmC/oxBS-seq does NOT solve (already logged, reconfirmed)
LINE-1 has ~500,000 nearly-identical copies; short bisulfite reads (oxBS or
plain BS) can't be uniquely mapped to one specific copy. Adding the 5hmC axis
adds a real new DIMENSION (methylated vs hydroxymethylated) but does not
bypass this locus-resolution problem -- same issue as flagged for plain
bisulfite in Q2/Q3.

### Idea (untested, mine, labeled as inference): does the tool's own
locus-disambiguation logic generalize?
retrotransposon-miner already resolves specific-insertion identity via
breakpoint position + flanking UNIQUE sequence + TSD -- it doesn't need to
map uniquely inside the repeat body. A bisulfite read spanning the boundary
between unique flank and a specific L1 copy could in principle be anchored
the same way. Would NOT help reads that fall entirely inside the repeat body.
Untested; needs a targeted literature check (does "amplicon/flank-anchored
locus-specific L1 methylation" already exist as a technique?) before treating
as a real lead. Not yet searched.

### Correctly-reasoned parts of the pasted text (no correction needed)
- Mouse->human monomer-structure caveat: consistent with, and correctly
  cross-references, this tracker's own Q4 finding.
- "Humans won't stay at one subfamily" flagged as unfalsifiable at any
  human-relevant timescale: consistent with Q4's verified rate (>=8 episodes
  / ~70 My primate evolution). Good instinct not to build on it.

---

## MEASURED DATA: chr22 milliDiv validation (2026-09-19, ours)

First real measurement of the session. Source: UCSC REST API, hg38 rmsk track,
chr22 (79,521 repeat rows; data timestamp 2022-10-18).
Artifacts: l1_divergence_chr22.png, l1_subfamily_divergence_chr22.csv

**Schema claim CONFIRMED empirically.** UCSC rmsk field order is
`bin, swScore, milliDiv, milliDel, milliIns, genoName, genoStart, genoEnd,
genoLeft, strand, repName, repClass, repFamily, ...` -- so milliDiv IS index 2,
and the repo's `_write_rmsk_mei_bed` reads exactly [5,6,7,9,10,11,12], skipping
it. Verified from data, not from memory.

**Result 1 -- divergence tracks subfamily age almost perfectly.**
Spearman rho = 0.995, p = 3.9e-12 over the age-ordered chain L1HS -> L1PA16
(13 subfamilies, 722 copies). Median divergence L1HS 0.011 -> L1PA16 0.147,
a 13-fold span. This validates the premise of the whole divergence pitch on
real data: milliDiv is a usable age/phylogeny proxy, empirically, on our own
test chromosome.
- chr22 carries 11,526 L1 copies = 4.42 Mb = 8.7% of the chromosome.
- Only 8 L1HS copies on chr22; 4 are >=5 kb. Young+intact is RARE -- which is
  exactly the population a novel-subfamily flag would be hunting in.

**Result 2 -- 5'UTR retention collapses with age.** (strand-corrected; see
caveat below)
L1HS 62.5% retain 5'UTR -> L1PA10/L1PA15/L1PA16 0%. Young (L1HS-L1PA4) mean
28.2% vs old (L1PA13+) mean 3.0%, a ~9-fold drop. Median element length
L1HS 3593 bp vs ~300 bp for old subfamilies (classic 5' truncation).
**Implication for the pitch: 5'UTR-restricted divergence is SELF-FILTERING.**
Only young elements still have a 5'UTR to measure, and young elements are
precisely the ones that could be novel/active. The measurement excludes the
ancient-fossil population automatically -- no tolerance-widening needed, which
is the opposite of the "extend the variability it accepts" approach rejected
earlier.

**Methods caveat (caught and fixed mid-analysis):** UCSC rmsk stores
repStart/repEnd/repLeft REVERSED for '-' strand elements. The naive
(strand-blind) calculation gave 59.5% vs 47.3% for young vs old 5'UTR
retention -- nearly flat, and wrong. Strand-corrected values are 28.2% vs
3.0%. Anyone redoing this must handle the strand encoding.

---

## Open question 8: the 869-tumor CRC claim -- side chat INVERTED the paper (2026-09-19)

Tracker previously carried, flagged do-not-cite: "one 869-tumor CRC study:
R^2=0.084, AUC<0.63" presented as evidence that averaged L1 methylation is a
**weak discriminator**. Found the paper. **The interpretation is backwards.**

Paper: Baba, Huttenhower, Nosho, ... Ogino (2010) *Mol Cancer* 9:125,
"Epigenomic diversity of colorectal cancer indicated by LINE-1 methylation in
a database of 869 tumors." DOI 10.1186/1476-4598-9-125

**What the model actually does:** LINE-1 methylation is the OUTCOME variable,
not the predictor. They regress L1 methylation ON ~20 clinical/molecular
features (age, sex, BMI, family history, smoking, tumor location, stage,
grade, CIMP, MSI, TP53, CDKN1A, CTNNB1, PTGS2, FASN, KRAS/BRAF/PIK3CA
mutations). The AUC measures "an ability of a regression model to diagnose
binary LINE-1 hypomethylation" -- i.e. how well everything else PREDICTS L1
methylation.

**So a low R^2 means the opposite of what was claimed.** ~8% of L1-methylation
variance explained by every other measured tumor feature => L1 methylation is
LARGELY INDEPENDENT of known tumor classifications. That is the paper's point
and its title: L1 methylation *indicates* epigenomic diversity -- it is a
distinct axis, not a redundant one. Low R^2 here is evidence FOR L1 carrying
unique information, not against it.
(Note: the numeric values R^2=0.084 / AUC<0.63 were not themselves re-verified
in this pass -- only the study design and what the statistics measure. Don't
quote the digits without checking the paper directly.)

Other real findings from it: L1 methylation ranged 23.1-90.3 (mean 61.4,
median 62.3, SD 9.6), approximately normal EXCEPT an excess of extreme
hypomethylators (<40; n=22, 2.5%, far more than normal-distribution
expectation), and those extreme hypomethylators were significantly associated
with YOUNGER patients (p=0.0058).

**Process lesson (second side-chat correction this session):** both pasted
side-chat syntheses have now needed correction on a load-bearing point. The
5hmC one was directionally fine with number/gap errors; this one inverted a
paper's conclusion. Side-chat output is agent-authored -- verify before it
enters the tracker as fact, same as anything else.

---

## Open question 9: the flank-anchored idea EXISTS -- bs-ATLAS-seq (2026-09-19)

Last turn I proposed, labeled as my own untested inference: anchor bisulfite
reads to a specific L1 copy using the unique flanking sequence, the same way
retrotransposon-miner identifies insertions. **This is published.**

**Lanciano et al. 2024, *Cell Genomics*** -- "Locus-level L1 DNA methylation
profiling reveals the epigenetic and transcriptional interplay between L1s and
their integration sites." Method: **bs-ATLAS-seq** (bisulfite ATLAS
sequencing), built on ATLAS, a suppression-PCR L1 insertion-mapping assay.
- Queries methylation of thousands of L1 copies AND their flanking DNA.
- Interrogates **the first 15 CpG sites of the L1 5'UTR promoter** -- exactly
  the CpGs this session identified as the silencing switch.
- Asymmetric paired-end sequencing yields genomic location AND methylation per
  locus -- i.e. the flank-anchoring mechanism, as guessed.
- Applied across 12 human cell lines; combined with Nanopore long reads for
  locus- and allele-specific resolution.
So: the reasoning was sound, the idea is not novel. Published Feb 2024.

**Findings from it that matter to us:**
- Young L1s hypermethylated in all cell lines EXCEPT embryonal lineages.
  (Consistent with Q1's ESC-specific causal literature.)
- Only 61 of 379 L1HS loci were "reactivable" vs 318 non-reactivable --
  locus identity, not just subfamily, determines reactivation potential.
- **Confirms a claim I earlier rated only moderate-confidence:** "the presence
  of an L1 element can dictate the methylation status of the adjacent genomic
  DNA," and "the DNA methylation status of L1 upstream sequences mirrors that
  of the L1 promoter itself." My correction to [user]'s Idea 1 (L1 -> local
  hypermethylation of neighbours) is now VERIFIED, not inferred.

**COUNTER-EVIDENCE the pitch must survive:** a separate study (inter- and
intra-locus L1 promoter methylation heterogeneity across bladder/colon/
pancreas/prostate/stomach cancers, 11 loci by pyrosequencing + NGS bisulfite)
reports that multivariate logistic regression "did not reveal a significant
clinical advantage of locus-specific methylation markers over global
methylation markers in distinguishing tumors from normal tissues." Locus-
specific is NOT automatically better for tumor-vs-normal discrimination.
Caveats: only 11 loci, pyrosequencing-era, and tumor-vs-normal is one
narrow endpoint. But this must be acknowledged, not ignored.
That study did also find specific CpGs within the L1 consensus are
preferentially hypomethylated ("seeds for hypomethylation") and that those
CpGs more likely face the histones in the nucleosome.

### THE ACTUAL REMAINING GAP (synthesis of Q7 + Q9)
bs-ATLAS-seq's own authors state the limitation explicitly: "5-methylcytosine
(5mC) and 5-hydroxymethylcytosine (5hmC) are both protected from
bisulfite-induced deamination, thus **bs-ATLAS-seq cannot discriminate between
these two DNA modifications**."

So the two threads meet exactly here:
 - bs-ATLAS-seq = locus resolution, CANNOT separate 5mC from 5hmC.
 - oxBS-seq     = separates 5mC from 5hmC, NO locus resolution for L1.
 - Nobody has combined them. The oxidation step is a pre-treatment; in
   principle it could sit upstream of the bs-ATLAS-seq library prep.
Closest prior work found: a 2017 mouse study doing methylation AND
hydroxymethylation on active L1 subfamilies -- but by subfamily-specific
pyrosequencing (TfI, A, GfII), i.e. subfamily resolution, not locus, and mouse
not human.
**Status: genuinely open as far as this search goes. NOT yet verified as
unattempted -- needs a dedicated search before anyone claims it.** This is a
wet-lab method proposal, not something retrotransposon-miner can do; it would
be a collaboration/assay idea, not a code change.

---

## MEASURED DATA 2: is locus-specific L1 methylation computationally reachable? (2026-09-19, ours)

Trigger: [user] asked whether there's a computational route to bisulfite work
without a wet lab.

**Category correction first:** bisulfite treatment is a CHEMICAL reaction
(sodium bisulfite deaminating unmethylated cytosines in a tube). There is no
computational substitute for generating it. What IS purely computational:
ANALYSING bisulfite data that already exists. Public WGBS is abundant (ENCODE,
TCGA, GEO) -- that is the real computational path, and it beats simulating
data, because it's real.

**Feasibility test run on chr22 rmsk (already downloaded).** Proxy for
"can a short read pin this L1 copy to its specific locus": how much
NON-REPEAT sequence flanks each copy? A read must reach into unique sequence
to be anchorable.

| threshold | % of chr22 L1 copies |
|---|---|
| any non-repeat flank | 84.8% |
| >=150 bp (one short read) | 40.3% |
| >=500 bp (amplicon-able) | 19.1% |

Per subfamily, >=150bp flank: L1HS 62.5%, L1PA5 66.2%, L1PA4 57.8%,
L1PA2 56.7% ... falling to L1PA15 8.9%, L1PA11 26.1%, L1PA10 29.7%.
**Young subfamilies are MORE anchorable than old ones** -- same self-filtering
direction as the 5'UTR-retention result. The elements you care about are the
ones you can actually address.

Genome-wide (naive x61 scaling from chr22, which is GC-rich and L1-poor so
this UNDERcounts): ~700k L1 copies, ~284k with >=150bp unique flank.

**Read-depth arithmetic** (junction-spanning reads with >=50bp on both sides
= coverage x 0.33 for 150bp reads): 10x WGBS -> ~3 reads/junction;
30x -> ~10; 60x -> ~20. At standard 30x, ~10 informative reads per junction is
enough for a rough per-locus methylation fraction. Arithmetic, not a
simulation -- treat as an order-of-magnitude feasibility estimate only.

**Caveat on the whole analysis:** "non-repeat flank" is a PROXY for unique
sequence. True mappability needs k-mer uniqueness against the whole genome,
not just absence of a RepeatMasker annotation. Numbers above are an upper
bound on anchorability.

### CORRECTION TO MY OWN FRAMING (assistant error, logged)
Last turn I asked [user] whether a novel-subfamily flag was "a detector for
something too rare to ever fire," citing 8 L1HS copies on chr22. **That used
the wrong denominator and the question was partly malformed:**
- A NOVEL subfamily is by definition NOT in the reference annotation. Counting
  reference L1HS copies cannot bound how often a novel-subfamily flag fires.
- The flag would evaluate every CALLED (non-reference) insertion, a different
  population from reference copies.
- For the METHYLATION project specifically the denominator is ~284k
  anchorable L1 loci genome-wide, not 8. [User]'s "it's not even that rare"
  is CORRECT for that project; my rarity worry does not apply to it.
The rarity question remains genuinely open ONLY for the narrow case of
"how often would a novel-subfamily flag actually fire on real callsets" --
which is an empirical question for the Kilo harness, not something reference
counts answer.

### Still unanswered and blocking: predictor of WHAT?
[User] proposes making this "a predictor." A predictor requires a named
outcome variable + labeled training data. Not yet specified. Until it is,
this is a measurement, not a predictor. Posed back 2026-09-19.

### Speculation flagged, not logged as science
"ten years into the future, who knows" / "it could have many impacts" -- not
arguments, same category as the earlier "humans won't stay at one subfamily"
hedge. Fine as motivation; cannot appear in a pitch to William.

---

## BUILT: l1meth -- locus-vs-global feasibility study (2026-09-19)

Outcome chosen by [user]: tumour vs normal tissue. Built standalone in the
workspace (NOT in the repo -- Kilo is working there). 24 tests passing.
Artifacts: l1meth/{catalog,simulate,compare,panel}.py, README.md,
tests/test_l1meth.py, run_study.py, l1_locus_catalog_chr22.csv,
l1_locus_vs_global.png

**Framing, stated up front in the README:** no bisulfite data was analysed.
The catalogue is real (reference coordinates, real CpG positions); the
methylation is simulated from published distributions. The AUCs are a
statement about STUDY DESIGN, not a validated classifier.

### Real catalogue: 34 addressable loci on chr22, 418 promoter CpG sites
(out of 11,526 L1 copies). Naive genome scaling ~2,100 -- same order as the
"thousands of L1 copies" bs-ATLAS-seq profiles, an independent sanity check.

**TWO methodological traps found by working on real data, both fixed:**
1. *Short-consensus subfamilies.* L1ME4b consensus = 892 bp, L1ME4c = 830 bp,
   L1ME3G = 924 bp -- SHORTER than the 900 bp promoter window. A naive
   `cons_start < 900` test passes for every copy of those subfamilies. First
   catalogue built came out as 278 loci dominated by L1ME4b/L1ME4c/L1MEc --
   degraded ancient fragments with no promoter at all. Fixed with a
   >=5 kb consensus-length guard plus an L1P (primate lineage) restriction,
   since the human-type 5'UTR promoter postdates the L1M subfamilies.
   NOTE: this ALSO corrects the earlier "284,000 anchorable loci genome-wide"
   figure in this tracker -- that number required only a >=150bp flank and
   counted 3' fragments with no promoter. The promoter-bearing number is
   ~2,100, two orders of magnitude smaller.
2. *Assembly gaps.* One locus reported a 50 kb flank that was centromeric N
   sequence. N runs are not mappable. Recomputed flanks as contiguous non-N,
   non-repeat sequence outward from each junction; 15 loci changed, 1 dropped
   from the addressable set (35 -> 34). Aggregate anchorable fraction
   unchanged at 40.4%.

### Findings (simulated, calibrated to Baba 2010)
1. **Locus resolution is not free.** With ONLY a global hypomethylation shift
   and no locus-specific component, per-locus scores WORSE than the bulk
   average: AUC 0.853 vs 0.908. The published null is reproducible from this
   alone -- the extra parameters cost variance and buy nothing.
2. **It wins as soon as locus-specific biology exists.** At >=5% of loci
   affected: 0.957 vs 0.910 (delta +0.046).
3. **Panel SIZE does not explain the published null; panel SELECTION does.**
   A blind 11-locus panel from a 2,100-locus pool with 10% signal has only a
   **72.7% chance of containing even one signal locus**. ~1 study in 4 reports
   a null regardless of whether the biology is real.
4. **Reference-computable features fix that before sequencing anything.**
   Ordering by divergence-from-consensus and taking the youngest raises the
   hit rate 72.7% -> 96.0% at the SAME panel size, capturing 2.3x more signal
   loci (2.71 vs 1.14). This is the concrete link from the divergence pitch
   (PR #51 follow-on) to the methylation project: the tool's output becomes
   the panel-selection input.

### Load-bearing assumptions (named in README, not buried)
- Normal colonic mucosa ~76% methylation -- weakest calibration link, not
  traced to a specific source; exposed as a parameter.
- Age enrichment of signal (8x youngest vs oldest) -- an ASSUMPTION motivated
  by Lanciano 2024 (only 61/379 L1HS loci reactivable), magnitude not
  measured. Finding 4 depends entirely on it; setting age_enrichment=1.0
  makes the advantage vanish, and that case is covered by a test.
- "Non-repeat non-N flank" is a proxy for mappability; true uniqueness needs
  genome-wide k-mer analysis. Counts are an upper bound.
- chr22 only, and chr22 is GC-rich/L1-poor so scaling likely UNDERcounts.

### Next step (not done)
Replace the simulated matrix with real public WGBS (ENCODE/TCGA/GEO) run
against this catalogue. That converts the whole thing from a design study
into a measurement. Nothing else in the study needs to change -- compare.py
takes a matrix and labels either way.

### CORRECTION to the entry above (review catch, 2026-09-19)

The "L1ME4b consensus = 892 bp" claim above was stated as a subfamily fact but
derived from a MIS-COMPUTED column, and the catalogue numbers were wrong.

**The bug:** consensus length was computed as `repEnd + |repLeft|` for every
row. UCSC orders the three consensus-coordinate fields differently by strand --
for '-' rows `repStart` carries the negative remainder and `repLeft` the
consensus START. Applying the '+' formula to '-' rows adds a coordinate to a
length. L1ME4b's single 892 bp model came out as 179 distinct values spanning
87-11,799 bp; the median happened to land near the true value, so the claim
read correctly by luck.

**Fixed:** `consensus_length()` branches on strand; `subfamily_consensus_length()`
reduces rows to one value per subfamily by MODE (consensus length is a property
of the subfamily model, not of a copy). Two regression tests added (26 total).

**What changed in the results:**
- Addressable loci **34 -> 43**; promoter CpG sites **418 -> 512**. The
  correction ADDED loci -- the corrupted column had been spuriously excluding
  primate copies. Genome-wide scaling ~2,100 -> ~2,600.
- Zero-signal penalty: locus 0.838 vs global 0.897 (was 0.853/0.908).
  Five-percent advantage: 0.959 vs 0.932 (was 0.957/0.910).
- **No conclusion changed.** Locus resolution still costs accuracy when no
  locus-specific biology exists, still wins from ~5% signal, and the panel-
  selection result (72.7% -> 96.0%) is untouched since it never used cons_len.

**Also corrected:** L1ME3G (924 bp) is NOT shorter than the 900 bp promoter
window as previously written -- it is barely longer. Only L1ME4b (892) and
L1ME4c (830) fall below it. 16 of 125 subfamilies are below the 5 kb
full-length guard.

---

## Open question 10: the age-enrichment assumption is FALSIFIED for carcinomas (2026-09-19)

Chased the question posed last turn -- "how would you test whether young loci
really are the differentially-methylated ones, using data that already
exists?" -- by fetching the bs-ATLAS-seq paper's full text. **They already ran
that exact test, and the answer is no.**

Source: Lanciano et al. 2024, *Cell Genomics*, DOI 10.1016/j.xgen.2024.100498
(full text fetched and read 2026-09-19).

**Their design:** all publicly available GEO Illumina 450K array datasets
(~12,000 samples), young L1PA probes (L1HS-L1PA3-short, n=695) vs old
(L1PA3-long-L1PA8, n=189), two-sided Wilcoxon rank-sum. This is the
already-exists-data test, done at scale.

**Their result, verbatim in substance:** methylation of young L1PAs is "high
and generally similar to or higher than that of old L1PAs in most
situations." Young L1PAs are hypomethylated specifically in pluripotent stem
cells (ESC/iPSC), trophoblast, embryonal carcinoma, seminoma, placenta, fetal
membranes and hydatidiform moles -- early embryogenesis, extra-embryonic
tissue, male germline tumours. In normal fibroblasts AND several carcinoma
lines (MCF-7, HeLa-S3, HepG2) young L1PAs are relatively HYPERmethylated.
The abstract states it directly: youngest primate L1 families are
hypomethylated in pluripotent stem cells and placenta "but not in most
tumors."

**Consequence for our study:** finding 4 (age-informed panel selection) was
built on an assumed 8x enrichment of signal in young loci. Measured impact of
removing it: age-informed selection goes from +23.3 points over blind
(96.0% vs 72.7%) to **-3.6 points (69.7% vs 73.3%)**. The advantage was
entirely the assumption.

**Fixed in code, not just prose:** `PoolConfig.age_enrichment` default changed
8.0 -> 1.0; module docstring carries the evidence; two tests pin both
behaviours (27 passing). README finding 4 rewritten as FALSIFIED with the
failure retained as the result. Figure panel B now plots all three curves so
the gap between assumption and evidence is visible.

**What survives:** findings 1-3 never used age_enrichment and are untouched --
locus resolution costs accuracy when no locus-specific biology exists
(0.829 vs 0.893), wins from ~5% signal, and the published null is explained
by panel SELECTION (blind 11-locus panel misses signal ~27% of the time)
rather than panel size.

**Specific redirect, not a dead end:** the young-L1 hypomethylation the
strategy needs IS measured -- in germ-cell and embryonal tumours. A
locus-resolved study of seminoma or embryonal carcinoma bets on an observed
effect; carcinoma tumour-vs-normal does not.

**Caveat in the other direction (do not overstate the falsification):** the
450K analysis compares young-vs-old methylation LEVELS WITHIN a sample, while
panel selection needs young-vs-old DIFFERENTIAL methylation between tumour and
normal. Related but not identical -- strong evidence against, not a direct
refutation.

### Other findings from the same paper, logged
- L1 methylation state propagates to the proximal region **up to 300 bp** --
  a concrete distance for the neighbour effect verified in Q9.
- **L1 hypomethylation alone is typically insufficient to trigger L1
  expression**, due to redundant silencing pathways. Relevant to the earlier
  "altered L1" thread: demethylation is necessary, not sufficient.
- TF binding at hypomethylated L1s: YY1, and ESR1 in MCF-7 (ESR1 binds ~25% of
  L1 loci there; knocking it down limits L1 expression).
- Their panel: 12 human cell lines; ~312 full-length L1HS per cell line.
  Public portal: https://L1methdb.ircan.org

### Our own new measurement this turn
Genic context added to the catalogue from UCSC ncbiRefSeqCurated chr22
(1,960 transcripts, 667 genes). Of 43 addressable loci: 27 intergenic,
15 intronic, 1 exonic. Addressable loci are LESS intronic than L1 overall
(34.9% vs 52.6%) -- anchorability requires non-repeat flanks. Ten of the 15
intronic ones sit in TTC28. NOTE: TTC28 intron 1 is, from memory, a known
somatic L1 insertion hotspot in colorectal cancer -- NOT verified this
session, do not repeat without checking.
Deliberately NOT used as a selection feature: Lanciano link intronic L1
methylation to host transcription WITHIN a sample, which is not evidence of
more tumour-normal differential signal. Replacing one unmeasured selection
assumption with another would repeat the mistake finding 4 just made.

### CORRECTION to Q10 (review catch, 2026-09-19) -- two numeric errors

**1. Panel-selection deltas were quoted from low-replication runs, and one
from a stale pool size.** The "-0.8 pp (1.17 vs 1.10)" figure predates the
n_pool 2100 -> 2600 change; "-3.6 pp" came from a 300-replicate run. Recomputed
at 1,500 replicates with the shipped defaults:

| case | blind | age-informed | delta |
|---|---|---|---|
| assumed 8x enrichment | 70.7% | 94.8% | +24.1 pp (1.11 -> 2.66 signal loci) |
| evidence-supported 1x | 70.4% | 68.2% | -2.2 pp (1.11 -> 1.10 signal loci) |

Sampling SE on a difference of proportions at 1,500 reps is ~1.7 pp, so the
1x delta is **not distinguishable from zero**. The honest claim is "gains
nothing," NOT a specific negative number -- quoting -0.8 or -3.6 pp as a
result overstates the precision. Conclusion unchanged: the advantage was the
assumption. Blind 11-locus panels miss signal ~30% of the time (was written
as 27%).

**Process note:** the first recomputation in-kernel still reported the 8x
behaviour because the edited `panel.py` had not been reloaded into the running
kernel -- the module-level default had changed on disk but not in memory.
Caught by asserting `PoolConfig().age_enrichment == 1.0` before trusting the
output. Worth doing routinely after editing a module mid-session.

**2. TTC28 locus count was wrong.** I wrote "ten of the 15 intronic
addressable loci sit in TTC28"; counted from the saved catalogue it is
**7**. Full breakdown of host genes among the 15 intronic addressable loci:
TTC28 7, GAB4 3, LINC02885 1, LOC339685 1, LINC01399 1, SYN3 1, LARGE1 1.
The "ten" came from eyeballing a `head(10)` display rather than counting.
The TTC28-as-somatic-L1-hotspot recollection remains UNVERIFIED -- still do
not repeat it without checking.
