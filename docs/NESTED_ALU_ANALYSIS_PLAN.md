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

## Phase 3 — genotype concordance (revised §5)

Per-sample callsets only, under the §4 assignment and the §5 revision.

---

## Deliverables

| # | deliverable | produced by | cohort consumed |
|---|---|---|---|
| 1 | `per_call_longread_nested.csv` — one row per Alu insertion site | Phase 0 | pooled BCF |
| 2 | `provenance_manifest.json` — verdict of every gate | Phase 0 (Gate 0c updated per §5) | pooled BCF + per-sample (0c) |
| 3 | `position_enrichment.md` / `.json` — Phase 1 primary + sensitivities | Phase 1 | pooled BCF (unique sites) |
| 4 | `unique_sites` table — deduplicated sites per the §2 matching rule, per-sample call counts side by side | cross-sample dedup (§1–§2) | pooled BCF + per-sample callsets |
| 5 | `matching_audit.md` — per-source caller offsets, calibration sites, full ±5/±20 sweep | §2 calibration (before new-sample matching) | already-processed shared sites |
| 6 | tail-cluster TPRT-hallmark audit (tail + comparison bins, blind) | Phase 2a | pooled BCF (unique sites) |
| 7 | sentinel-excluded TSD distribution + drift stratification | Phase 2b | pooled BCF |
| 8 | recurrence and carrier-count table (per-sample call counts side by side with unique sites) | §4 assignment | per-sample callsets |
| 9 | **accumulation curves** — unique nested sites vs genomes added, with marginal yield; basis for genome-count decisions (§3) | §3 procedure | per-sample callsets + pooled BCF |
| 10 | genotype concordance report (callability inference labeled) + private-site exploratory analysis | revised Phase 3 | per-sample callsets |

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
