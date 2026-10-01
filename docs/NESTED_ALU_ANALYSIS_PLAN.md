# Nested-Alu analysis plan — structure-informed study

**Plan consolidated and cross-sample design pre-registered 2026-10-01, before any
new per-sample genome was processed.** The matching rules in the next section
were frozen at that date and may not be revised after cross-sample results are
seen.

Provenance of this document: it consolidates commitments that already exist in
committed, executable form — the Phase 0 gates and positional conventions in
`scripts/build_longread_nested_cohort.py` (module docstring, the authoritative
statement), the Phase 1 pre-registration in
`scripts/test_position_enrichment.py` (module docstring, the authoritative
statement), and the Phase 1 result report
(`nested_analysis/results_longread/position_enrichment.md`). Phase 2a/2b and
Phase 3 commitments are stated here for the first time as a single plan. Where
this document and a script docstring disagree, the script docstring wins and
this file is corrected — it is the plan of record for *cross-sample* decisions
only.

**Version-control provenance.** The plan first entered version control in
commit `387eff3` (2026-10-01), bundled with the concurrent boundary-perturbation
falsification work because both staging trees shared one index. The
pre-registration date stands regardless of that bundling: the cross-sample
design decisions below were written and frozen **2026-10-01, before any new
per-sample genome was processed**, and the matching rule was frozen before
new-sample data existed. This note only records how the commit came to be
shared; it does not alter the freeze.

---

## Phase 0 — cohort build (complete)

`build_longread_nested_cohort.py` builds the nested-Alu cohort from the pooled
long-read BCF behind four executable provenance gates, each of which writes its
verdict to `provenance_manifest.json`:

- **Gate 0a** — BCF `.csi` verified against a linear scan (correct despite
  mtime), FASTA `.fai` checked independently.
- **Gate 0b** — no CHM13 contamination; only `svim.asm.hg38.*` BCFs accepted.
- **Gate 0c** — independence verdict for the genotype axis (HG03086 is inside
  the pooled cohort). Runs **first**, before any per-sample data is processed;
  its verdict question is revised in Phase 3 below, the gate itself is not.
- **Gate 0d** — `TSD_LEN` sentinel pile-ups excluded from continuous use;
  sentinel-free distribution measured rather than assumed.

Positional conventions (`host_offset_5p_0based`, `consensus_offset`) and the
three-column genotype state are defined in the Phase 0 docstring and are not
restated here.

---

## Cross-sample design (pre-registered 2026-10-01)

**Freeze statement.** Everything in this section was fixed before the new
per-sample genomes were processed. Any deviation must be reported as a
deviation with its date and reason, not silently absorbed into the analysis.

### 1. Estimand statement

The positional analyses — Phase 1 (pos-133) and Phase 2a (tail cluster) —
operate on **unique insertion sites**, deduplicated across samples per the
matching rule in §2. This matches the accumulated-site estimand of Levy,
Schwartz & Ast 2010: a catalogue of sites, each counted once, with carrier
information carried separately.

Private/single-sample sites are a **separate exploratory analysis** (recent,
unselected insertions). They are never merged into the primary positional
test, and no primary-test p-value is computed over them.

Per-sample call counts and unique-site counts are always reported side by
side, in every table, including the deliverables in §Deliverables.

### 2. Matching rule (verbatim, frozen)

> same chrom + family + orientation + host + breakpoint ±10 bp primary
> (sweep ±5/±20 reported in full); no transitive chaining; per-source
> systematic caller offsets calibrated on high-confidence shared sites and
> recorded in matching_audit.md

Non-normative clarification (the rule above is the frozen text): "per-source"
means per callset source — a systematic offset between two callsets' breakpoint
conventions is calibrated, applied, and recorded in `matching_audit.md` before
cross-callset matching runs. The operational definition of that offset
(estimator, calibration-set criteria, application point) is fixed in
`matching_audit.md` before any new-sample matching is run. Calibration uses
high-confidence shared sites from **already-processed** data, so the
calibration is not tuned on the samples whose sites will be deduplicated. The
±5 and ±20 sweeps are reported in full regardless of outcome. **No transitive
chaining**: two sites match a third only by directly satisfying the rule, not
by being linked through each other.

### 3. Power statement (revision)

The earlier projection of ~200 nested L1s at 8 genomes is **retired**. Genome
count is decided by the measured **marginal yield of unique nested sites** —
how many *new* unique nested sites each additional genome contributes — read
off an accumulation curve, not by multiplying per-genome call averages.
Genome-count decisions follow the curve: additional genomes are collected only
while the measured marginal yield of unique nested sites supports the test
being powered, and the stopping threshold is written down before the next
batch is unblinded for any positional test.

Shared polymorphisms add **carrier counts, not independent positional
observations**. A site carried by N samples contributes one observation to a
positional test and N to a carrier count; the two are never summed or
averaged into a single quantity.

### 4. Cohort assignment

| cohort | used for |
|---|---|
| pooled BCF (catalogued sites) | catalogued-site positional analyses: Phase 1 pos-133, Phase 2a tail |
| per-sample callsets | recurrence, carrier counts, private-site analysis, and revised Phase 3 |

When a pooled-cohort positional stratum needs carrier weighting, carrier counts
come from per-sample callsets under the §2 matching rule — not from the pooled
BCF's genotype columns alone, whose `0/0` conflates genuine absence with
genotyping failure (documented in the Phase 0 docstring). Pooled-BCF carrier
annotation remains valid only as a descriptive stratum of the pooled cohort
itself.

### 5. Phase 3 revision (genotype concordance)

- Multiple per-sample callsets from the same pipeline are **independent
  individuals but NOT independent methods**. Method independence would require
  a different caller, reference, or assay; N same-pipeline callsets are one
  method applied to N genomes.
- Terminology remains **"genotype concordance"**. It is method-internal
  agreement, not validation; the word "validated" is reserved for cross-method
  evidence that does not currently exist.
- An absent call in another sample counts as a negative **only at loci
  demonstrated callable in that sample**. Callability is an *inference* and is
  labeled as such in every table that uses it; the demonstration criteria are
  recorded before use and never read off the same absence being tested.
- **Gate 0c** still runs first, unchanged as a gate. Its verdict question is
  updated to the revised Phase 3: *"may Phase 3 use the genotype axis for
  concordance at demonstrated-callable loci, given that same-pipeline callsets
  are independent individuals but one method?"* — replacing the earlier binary
  use/not-use question.

---

## Phase 1 — position-133 enrichment (pre-registered; complete on the pooled cohort)

Pre-registered in `scripts/test_position_enrichment.py`; results in
`nested_analysis/results_longread/position_enrichment.md`. Commitments:
primary window consensus 133 ± 5 bp with a pre-declared slop series;
host-interval opportunity null as primary, mappability-scaled null as
sensitivity only (a short-read track is not this cohort's callability), EN-motif
conditioning as Null C; 16-bin scan corrected as a scan with max-statistic
permutation FWER; host-clustered uncertainty with leave-one-host-out;
site and carrier counts never collapsed. Under the §1 estimand, the cross-sample
rerun operates on unique sites per §2.

## Phase 2a — tail-cluster hallmark audit

The 280–300 bp tail cluster and pre-specified comparison bins are audited blind
against TPRT hallmarks (TSD reconstruction from both junctions, 3′ A-rich
sequence, strand-compatible EN motif at the nick). Comparison bins are fixed
before the audit; calls with favorable hallmarks are never selected into the
analysis set. Under the §1 estimand the audit operates on unique sites per §2.

## Phase 2b — TSD structure

`TSD_LEN` sentinels are excluded per Gate 0d; the sentinel-free distribution is
measured, and offset drift is stratified rather than assumed away. TSD length is
not treated as continuous until sentinels are excluded.

## Phase 3 — genotype concordance (revised §5; complete)

`scripts/score_genotype_concordance.py`; results in
`nested_analysis/results_longread/genotype_concordance.md` / `.json`, with the
per-site table in `phase3_site_genotypes.csv`.

Two pre-registered requirements were measured and found unsatisfiable, and both
are recorded as gate verdicts rather than worked around:

- **The per-sample callset §4 assumes does not exist on GRCh38.** The
  `final-vcf.unphased.SVAN_1.3.vcf.gz` under the hg38 polymorphism directory is
  a symlink into `hs1/` and is therefore CHM13. The genotyped SVIM pooled BCF
  (`svim.asm.hg38.bcf`, 908 samples) is the usable substrate, and §5
  anticipated it: independent individuals, one method.
- **The pre-registered concordance design is not executable.** The callset has no
  missingness channel, so a `0/0` cannot be separated from a genotyping failure
  and §5's "demonstrated callable" requirement cannot be met (Gate 3b).

Phase 3 therefore measures **recurrence across two independent halves of the 908
genomes**, split by sorted sample name, with the hallmark as exposure and
recurrence as outcome. That ordering is a deliberate reversal of Phase 2b: a
site with one carrier cannot be cross-half supported by arithmetic, so exposing
recurrence reduces the reference arm to singletons and every ladder rung fails
its minimum-reference floor. Allele frequency is excluded from the covariates
for the same reason.

New gates, all executable and all in the JSON:

- **Gate 3a** — the site and genotype callsets declare different lengths for 24
  of 25 shared contigs, with the sign flipping across chromosomes. The join is
  verified by exact-coordinate containment (0 of 1,335 analysis sites absent)
  rather than trusted from the headers, and the header disagreement is reported
  because a join that had trusted it would have been wrong.
- **Gate 3b** — no missingness channel; 0 missing genotypes of 1,212,180.
- **Gate 3c** — samples joined by name. The cross-method callset orders its
  samples differently, so an index join would transpose genotypes silently.
- **Gate 3d** — the cross-method floor. An independent short-read 1KG callset
  covers 0.90% of nested-Alu sites against 41.33% of non-nested Alu sites from
  the same cohort (Fisher exact p = 5.8e-12). A short read cannot span an Alu
  inserted into an Alu, so the comparator is structurally blind rather than
  noisy. This is the measured evidence for §5's refusal to use the word
  "validated", and it replaces an assumption with a number.

Result: **no TPRT hallmark predicts cross-genome recurrence.** Seven of nine
pre-declared hallmark exposures were estimable; all matched. The one nominal
signal, long poly(A) (≥20 nt) being *less* recurrent, is MH RD −0.122 at raw
p = 0.0073 and **BH p = 0.051** across exposures, so it does not survive
multiplicity. This is consistent with Phase 2b's finding that the triad does not
separate nested from non-nested calls.

Accumulation (deliverable 9) over the 908 genomes: 1,328 unique nested sites
backed by 103,196 carrier observations, 77.7 carriers per unique site; half the
final unique-site count is reached at 253 genomes, 90% at 687. Site counts and
carrier counts are reported side by side and never summed.

---

## Deliverables

| # | deliverable | produced by | cohort consumed |
|---|---|---|---|
| 1 | `per_call_longread_nested.csv` — one row per Alu insertion site | Phase 0 | pooled BCF |
| 2 | `provenance_manifest.json` — verdict of every gate | Phase 0 (Gate 0c updated per §5) | pooled BCF + per-sample (0c) |
| 3 | `position_enrichment.md` / `.json` — Phase 1 primary + sensitivities | Phase 1 | pooled BCF (unique sites) |
| 4 | `unique_sites` table — deduplicated sites per the §2 matching rule, per-sample call counts side by side | cross-sample dedup (§1–§2); *produced at n=2 callsets by the sibling `rtm-nested-analysis` project, not yet at 908* | pooled BCF + per-sample callsets |
| 5 | `matching_audit.md` — per-source caller offsets, calibration sites, full ±5/±20 sweep | §2 calibration (before new-sample matching); *same n=2 caveat* | already-processed shared sites |
| 6 | tail-cluster TPRT-hallmark audit (tail + comparison bins, blind) | Phase 2a | pooled BCF (unique sites) |
| 7 | sentinel-excluded TSD distribution + drift stratification | Phase 2b | pooled BCF |
| 8 | recurrence and carrier-count table (per-sample call counts side by side with unique sites) | Phase 3 — `phase3_site_genotypes.csv`, carrier names listed per site | genotyped BCF (908) |
| 9 | **accumulation curves** — unique nested sites vs genomes added, with marginal yield; basis for genome-count decisions (§3) | Phase 3 — `genotype_concordance.md` §Accumulation | genotyped BCF (908) |
| 10 | genotype concordance report (callability inference labeled) + private-site exploratory analysis | Phase 3 — `genotype_concordance.md` / `.json` | genotyped BCF (908) |

---

## Disclosures

- The cohort is a catalogue of insertion **sites**, not one genome's biology.
  Catalogue polymorphisms have passed selection; these are survivor
  distributions, not raw de novo integration rates.
- This is an enrichment/concordance study. It is not a detection-accuracy
  result and licenses no precision/recall claim.
- Positions are consensus-relative through host-to-consensus alignment;
  projection is exact only for the unambiguous subset.
- Same-pipeline per-sample callsets are independent individuals but one method;
  concordance is method-internal agreement, not validation.
- **Shared polymorphism dedup is pre-registered; matching window sensitivity
  reported in full.** Deviations from the §2 rule are reported as deviations,
  with dates.
- **Output-directory collision (operational).** `nested_analysis/` sits at the
  workspace root rather than inside any one worktree, and both this worktree and
  the sibling `rtm-nested-analysis` project write into
  `nested_analysis/results_longread/`. Filenames do not currently overlap, but
  the two projects are one directory apart with no coordination, so a future
  filename reuse would silently overwrite the other project's result. Each
  producer should own a subdirectory.
- **Phase 0 is byte-reproducible.** Re-running `build_longread_nested_cohort.py`
  from scratch reproduced `per_call_longread_nested.csv` with an identical MD5,
  so cohort differences between runs are not a source of variation here.

---

## Phase 4 — joint signature and same-host recurrence (post-freeze extension)

**This section was added after the 2026-10-01 freeze and changes nothing above
it.** Everything above remains as frozen; this records work that extends the plan
rather than revising it. Scripts: `scripts/joint_enrichment.py`,
`scripts/recurrence_test.py`, shared inputs in
`scripts/nested_multi_sample_common.py`. Results in
`nested_analysis/results_phase4/`.

### Substrate

Consumes `nested_analysis/results_multi_sample/unique_sites.csv` (763 unique
sites, 5 genomes) as the authoritative site grouping. The per-call table named as
an input **does not exist**; per-call TSD and child subfamily are recovered by
joining the dedup output back to the callsets by sample name within the frozen
±10 bp rule. That join is a **gate**, not a diagnostic: it passes only if it
reproduces dedup's own `n_carriers` for every site and agrees on nesting state and
orientation. On the shipped data all 763 sites resolve, all 1,284 recovered calls
match `sum(n_carriers)`, and there are zero disagreements.

### Opportunity null — a measured deviation

No mappability or gap track exists in this workspace. The RepeatMasker substitute
was measured rather than assumed. A strict repeat mask is **degenerate**: the host
element is itself a RepeatMasker annotation, so masking repeats deletes the entire
opportunity. Once the host's own annotation is excluded, repeats *other* than the
host intrude on **0.000** of the host interval for **all** 216 measured Alu hosts.
The mappability-scaled variant therefore collapses onto the host-interval null —
which is the plan's own designated primary — and is recorded as a deviation.

### Joint enrichment

Pre-specified cells only; the ~1.2M-cell table is not scanned. Host-stratified
conditional randomization, no log-linear model, ≥10,000 replicates, MC
p = (1 + #{sim ≥ obs})/(B+1), structural-zero bins excluded rather than
pseudocounted.

| cell | observed | expected | effect | 95% CI | MC p | Holm p |
|---|---|---|---|---|---|---|
| `{Alu, 120-140, sense}` (primary) | 22 | 7.97 | 2.76× | 1.85–4.00 | 9.999e-05 | n/a |
| `{Alu, 280-300, sense}` | 14 | 4.92 | 2.85× | 1.58–4.33 | 9.999e-04 | 0.003 |
| interaction 120-140 | 5 | −0.43 | unstable\* | — | 0.100 | 0.200 |
| interaction 280-300 | −3 | 0.58 | unstable\* | — | 0.917 | 0.917 |

\* expected within 1 of zero, so the ratio and its interval are not quotable.

The primary cell **replicates the position-133 signal on an independent
5-genome cohort**. Per §1, private sites are a separate exploratory analysis and
no primary p-value is computed over them; they are reported for transparency
(324 private sites, observed 31, expected 11.65).

**Multiplication test.** The joint estimate (2.76×) exceeds the per-factor
benchmark (2.48×) by 1.11×. No multiplied figure is reported. The excess is
small and the pre-specified interaction cells — which test the same question
directly — do not reject after Holm adjustment, so the excess must not be read as
a detected interaction.

### Same-host recurrence

Definitions applied verbatim; **no pair is adjudicated in code**. One finding
about the *definitions* came out of this and is recorded because it changes how
the counts may be read:

- The criterion ORs its discordance conditions, so a child-subfamily difference
  is sufficient on its own. An earlier revision tested "is a TSD missing?" first,
  which let absent data mask positive evidence already present. Fixed, with a
  regression test.
- **479 of 1,284 nested calls report no TSD at all.** A missing TSD is not a
  TSD of length zero, so those pairs are `tsd_unevaluable` and enter neither
  headline count.
- At ±10 bp: **2 candidate-recurrent pairs, 1 unevaluable, 0 IBD** among 3
  same-host pairs. Per family: Alu 37/514 host copies carry >1 event, L1 5/161,
  SVA 2/27, other 1/11.
- **The candidate count does not clear its chance expectation.** The positional
  chance count is 4.31 against 3 observed. Worse, the definition is loose: two
  genuinely independent insertions drawn from this cohort already satisfy the
  candidate criterion with probability 0.949 (Alu). The converse matters as much —
  two independent pairs satisfy the *IBD* definition with probability ≈0, because
  the cohort carries 89 distinct child subfamilies. An empty IBD column is
  therefore **not** evidence that no shared events exist.
- Denominators are 5 genomes and are **not comparable** to the published pooled
  SVAN 26/2,559 Alu-host figure; both are labelled rather than aligned.

The honest reading: the test was run, candidate pairs exist, and it did not
resolve the question, because child subfamily and TSD do not discriminate at this
scale.
