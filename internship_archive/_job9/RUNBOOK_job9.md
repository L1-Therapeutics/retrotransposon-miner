# Job 9 — exact commands

Everything below was **rehearsed on your own interpreter** (`_verify/.venv-verify/bin/python`,
Python 3.15.0, pandas 3.0.6, numpy 2.5.4, scipy 1.18.1, pytest 9.1.1). The expected output is
copied from that run, not predicted. If your output differs by even one digit, send it to me.

## The four shell errors you have hit before, and how this avoids them

| What happened | Cause | How this job avoids it |
|---|---|---|
| `zsh: parse error near '&'` | multi-line Python pasted into zsh | everything is a `.sh` or `.py` file, invoked by path |
| `command not found: python` | no bare `python` on macOS | the script calls `_verify/.venv-verify/bin/python` by full path |
| `fatal: not a git repository` | git run from the L1 root | Block 3 tells you which folder to `cd` into first |
| `zsh: event not found` | `!` inside double quotes | no `!` in any command below |

---

# Block 1 — the one command (about 80 seconds)

```bash
bash ~/Desktop/Research/Research_Projects/L1/_job9/run_job9.sh
```

It checks the interpreter, checks the inputs exist, verifies the four scripts against
`CHECKSUMS.txt`, runs 38 tests, runs both analyses, and prints a verification block.

**Expect this at the top:**

```
[0/4] interpreter
      python 3.15.0 | pandas 3.0.6 | numpy 2.5.4 | scipy 1.18.1
[1/4] inputs present
      rmsk    148M
      gencode  89M
[2/4] checksums
alu_exon_proximity.py: OK
exon_distance_annotate.py: OK
test_alu_exon_proximity.py: OK
test_exon_distance_annotate.py: OK
[3/4] tests  (expect: 38 passed)
38 passed in 2.41s
```

**Expect this at the end — every line must say `OK`:**

```
  OK  copies_total                     114541   (expected 114541)
  OK  copies_in_single_strand_genes    67379   (expected 67379)
  OK  overall_sense_fraction           0.4501   (expected 0.4501)
  OK  upstream n_near                  1508   (expected 1508)
  OK  upstream antisense_frac          0.482095   (expected 0.482095)
  OK  upstream odds_ratio              0.753271   (expected 0.753271)
  OK  upstream p                       5.614667e-08   (expected 5.614667e-08)
  OK  EXON_SIDE=.                      784   (expected 784)
  OK  EXON_SIDE=upstream               537   (expected 537)
  OK  EXON_SIDE=downstream             531   (expected 531)

  ALL CHECKS PASSED
```

If any line says `BAD`, stop and send me that block. Do not use the numbers.

**What you have just reproduced.** The replication of Lev-Maor 2008: antisense Alu is depleted
within 150 bp upstream of an exon, 0.482 against 0.553 elsewhere in introns, p = 5.6 × 10⁻⁸ on
1,508 copies. And a working annotator that adds `EXON_DIST` and `EXON_SIDE` to a VCF. Those are now
your numbers, not mine.

---

# Block 2 — look at what you produced (optional, 10 seconds)

```bash
cd ~/Desktop/Research/Research_Projects/L1/_job9/out
column -t -s$'\t' aluY/zones.tsv
```

```
zone                   n      antisense  antisense_frac  median_milli_div
intronic_far           61130  33788      0.552724        69.0
overlaps_exon          3181   1714       0.538824        68.0
within_150_downstream  1560   821        0.526282        67.0
within_150_upstream    1508   727        0.482095        68.0
```

Then the new INFO fields on a real record:

```bash
grep -v '^#' HG03086.exondist.vcf | head -3 | cut -f1,2,8 | tr ';' '\n' | grep -E 'EXON_|^chr'
```

---

# Block 3 — the git work that is still open (5 minutes)

This is the one thing that has been waiting for several sessions. **Note the `cd` — git does not
work from the L1 root.**

```bash
cd ~/Desktop/Research/Research_Projects/L1/rtm-mei-gene-orient
git status
```

Expect `On branch feat/mei-gene-orient-field` and `nothing to commit, working tree clean`.
If it is not clean, stop and send me the output.

```bash
PYTHONPATH=src ~/Desktop/Research/Research_Projects/L1/_verify/.venv-verify/bin/python \
  -m pytest tests/test_gene_annotation.py -q
```

Expect `28 passed, 1 skipped`. The skip needs `RTM_LIVE_SNPEFF=1` and is correct.

**Do not run `pytest tests/` on this venv.** I tried it: it dies with
`36 errors during collection`, starting with `ModuleNotFoundError: No module named 'click'`.
The verify venv has pandas/numpy/scipy but not the pipeline's own dependencies, so only the
pure-Python test modules collect there. That is a venv limitation, not a problem with your branch.

```bash
git push -u origin feat/mei-gene-orient-field
```

Then open the pull request on GitHub. **Before you do, read the note in step 4 below** — one commit
on that branch may need dropping.

---

# Block 4 — the decision you have to make before opening that PR

The branch has three commits:

| Commit | What it does | Status |
|---|---|---|
| 9351c36 | fixes the gene-symbol register bug | **wanted** — identical to PR #81's fix |
| 3be6c01 | adds the shared-symbol tests | **wanted** |
| a1aec2d | adds the `MEI_GENE_ORIENT` INFO field | **William called this redundant** |

He said `MEI_GENE_ORIENT` can be derived from `ORIENT` and `GENE_STRAND`, and you agreed to drop
it. So open the PR with the first two commits only:

```bash
cd ~/Desktop/Research/Research_Projects/L1/rtm-mei-gene-orient
git checkout -b fix/gene-register-with-tests 3be6c01
PYTHONPATH=src ~/Desktop/Research/Research_Projects/L1/_verify/.venv-verify/bin/python \
  -m pytest tests/test_gene_annotation.py -q
```

Expect `24 passed, 1 skipped` here — four fewer than the full branch, because the four
`MEI_GENE_ORIENT` tests live on the commit you are leaving off. I verified both counts by
extracting each commit and running them on your interpreter.

```bash
git push -u origin fix/gene-register-with-tests
```

That gives a branch with the register fix plus its tests and nothing he has objected to. Open the
PR from `fix/gene-register-with-tests`, and in the description say it supersedes PR #81 — same
fix, more tests — then close #81.

I also confirmed at that commit that `MEI_GENE_ORIENT` does not appear anywhere in
`gene_annotation.py`, and that the register fix (`gid_to_symbol`) does — so the split lands
exactly the content you want and none you do not.

Keep `feat/mei-gene-orient-field` locally. Do not delete it; if William ever wants the field,
commit a1aec2d is there. Note there is also a *different* remote branch `feat/mei-gene-orient`
(no `-field`) at 8e82aea — do not push over it.

---

# Block 5 — what to send William

Two things, in this order.

**First, the replication.** It is a reproduction of a published result, which is the framing he
asked for:

> I reproduced Lev-Maor 2008's finding that antisense Alu is depleted just upstream of exons —
> the position where an antisense Alu's right arm can be spliced in as a new exon. On 67,379 fixed
> AluY copies inside single-strand genes, antisense fraction is 0.482 within 150 bp upstream of an
> exon against 0.553 elsewhere in introns (OR 0.753, p = 5.6 × 10⁻⁸). Downstream is much weaker
> (0.526, p = 0.039), which is the same upstream/downstream asymmetry they reported. It holds in
> all six Alu subfamilies I tested, from AluJo to AluY, and in every divergence quartile.
>
> It also resolves something that looked contradictory in our own data. Across the whole gene,
> Alu is depleted in the SENSE orientation. Within 150 bp of an exon, it is depleted in the
> ANTISENSE orientation. Both are real; they act over different distances. Termination selects
> against sense across the gene body, exonization selects against antisense next to an exon.

**Second, the splice question.** `splicing_annotation_proposal.md` section 7 has it drafted. The
short version is that `gene_annotation.py` disables snpEff's splice annotation (your commit
17a7462) for two reasons — a 3 GB OOM, and "an insertion is not a splice-site variant" — and the
3 GB figure predates the current `--snpeff-xmx` default of 8g.

---

# What I am doing next, so you do not duplicate it

1. The AluSx/AluJo replication of the **poly(A) carrier** result — still the gate on sending that
   analysis to William. The exon-proximity sweep is done; the poly(A) one is not.
2. Extending exon proximity to LINE-1 using older subfamilies, since L1HS puts only 7 copies in
   the window and cannot be tested.
3. Restricting the Alu analyses to introns proper using the exon union.
