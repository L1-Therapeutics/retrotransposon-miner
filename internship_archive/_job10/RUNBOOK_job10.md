# Job 10 — exact commands

Rehearsed end to end on **your** interpreter (`_verify/.venv-verify/bin/python`, Python 3.15.0,
pandas 3.0.6, twobitreader 4.0.2). Wall time on your machine: **5 min 16 s**. All 21 verification
checks passed. Every expected value below is copied from that run.

This job reproduces the result that was gated — the poly(A) carrier replication — plus the LINE-1
probe and its family controls.

---

# Block 1 — the one command

```bash
bash ~/Desktop/Research/Research_Projects/L1/_job10/run_job10.sh
```

Leave it running for about five and a half minutes. Expect this near the start:

```
[0/5] interpreter
      python 3.15.0 | pandas 3.0.6 | twobitreader ok
[1/5] inputs present
[2/5] checksums
[3/5] tests  (expect: 50 passed)
50 passed in 3.04s
```

Then the verification block. **All 21 lines must say `OK`:**

```
  OK  AluY copies in genes                   67379   (expected 67379)
  OK  AluY carriers                          11923   (expected 11923)
  OK  AluY adjusted OR                       0.9496   (expected 0.9496)
  OK  AluY canonical hexamers on consensus   0   (expected 0)
  OK  AluSx copies in genes                  75839   (expected 75839)
  OK  AluSx carriers                         17795   (expected 17795)
  OK  AluSx adjusted OR                      0.9379   (expected 0.9379)
  OK  AluJo carriers                         18598   (expected 18598)
  OK  AluJo adjusted OR                      0.9988   (expected 0.9988)
  OK  L1PA3-7 downstream OR                  1.561   (expected 1.561)
  OK  AluSx upstream OR                      0.809   (expected 0.809)
  OK  MIR upstream OR                        0.889   (expected 0.889)
  OK  L2 upstream OR                         1.563   (expected 1.563)
  ...
  ALL CHECKS PASSED
```

If any line says `BAD`, stop and send me the block.

**What you have reproduced.** The poly(A) carrier asymmetry in three Alu subfamilies: AluY
OR 0.950, AluSx 0.938 (stronger), AluJo 0.999 (absent). Plus the clean negative — zero canonical
AATAAA/ATTAAA on any of the three consensuses. Plus the LINE-1 near-exon result and three family
controls. These are now your numbers.

---

# Block 2 — the git work, if you have not done it yet

Unchanged from Job 9, repeated here so you do not have to go back. **Note the `cd`.**

```bash
cd ~/Desktop/Research/Research_Projects/L1/rtm-mei-gene-orient
git status
```

Expect `On branch feat/mei-gene-orient-field`, `nothing to commit, working tree clean`.

```bash
git checkout -b fix/gene-register-with-tests 3be6c01
PYTHONPATH=src ~/Desktop/Research/Research_Projects/L1/_verify/.venv-verify/bin/python \
  -m pytest tests/test_gene_annotation.py -q
```

Expect `24 passed, 1 skipped`. **Do not run `pytest tests/`** — that venv lacks `click` and dies
with 36 collection errors.

```bash
git push -u origin fix/gene-register-with-tests
```

Open the PR from `fix/gene-register-with-tests`, say it supersedes PR #81 (same fix, more tests),
then close #81. This branch deliberately leaves out the `MEI_GENE_ORIENT` commit William called
redundant.

---

# Block 3 — send William the replication

This is now unblocked. Suggested message:

> Following up on the Alu orientation work. I reproduced two things.
>
> First, a published result: antisense Alu is depleted just upstream of exons, which is the
> position where an antisense Alu's right arm can be spliced in as a new exon (Lev-Maor 2008).
> On 67,379 fixed AluY copies inside single-strand genes, antisense fraction is 0.482 within
> 150 bp upstream of an exon against 0.553 elsewhere in introns, OR 0.753, p = 5.6 x 10^-8.
> It holds in all six Alu subfamilies I tested and in every divergence quartile.
>
> Second, something of our own. The Alu consensus carries no canonical poly(A) signal at all —
> zero AATAAA or ATTAAA on either strand, in AluY, AluSx and AluJo alike — so any signal has to be
> created per copy. Copies that do carry one are more often antisense to their host gene:
> AluY OR 0.950 (p=0.014), AluSx 0.938 (p=0.0003). In AluJo it disappears (0.999, p=0.96),
> which I think is drift decoupling present-day sequence from selective history at that age.
>
> This also resolves something that looked contradictory in our data. Across the whole gene, Alu
> is depleted in the SENSE orientation. Within 150 bp of an exon it is depleted in the ANTISENSE
> orientation. Both are real and act over different distances — termination selects against sense
> across the gene body, exonization selects against antisense next to an exon.

**One caveat to include if he asks about AluJo:** its consensus rests on only 90 copies with 35
low-support columns, because the full-length filter removes 99.9% of AluJo. The per-copy result
(51,437 copies) is unaffected, but the AluJo *consensus* claim is thin.

---

# Block 4 — the splice question for William

`splicing_annotation_proposal.md` section 7, still unsent. The short version: `gene_annotation.py`
disables snpEff splice annotation (his own commit 17a7462) for two reasons — a 3 GB OOM and
"an insertion is not a splice-site variant" — and the 3 GB figure predates the current
`--snpeff-xmx` default of 8g.

---

# What I am doing next

1. Restricting the Alu analyses to introns proper using the exon union.
2. Replacing the 20-stratum Mantel-Haenszel with a logistic model on the continuous covariates.
3. The L2 anomaly: L2 is antisense-ENRICHED upstream of exons (OR 1.563, p = 1e-5), the exact
   opposite of Alu. Nothing in the map predicts it.
