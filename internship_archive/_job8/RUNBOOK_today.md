# Exact commands — today

## Read this first: four errors from your last session, and what each one was

Your paste shows four distinct failures. None is a problem with the analysis; all four are
shell-level, and all four are now impossible if you use the scripts in Block 3.

| What you saw | Cause | Fix |
|---|---|---|
| `zsh: parse error near '&'` | Python pasted straight into zsh. `&` is zsh's background operator, so it never reached Python. | Python goes in a `.py` file, run as a file. See Block 3b, never a paste. |
| `zsh: command not found: python` | The venv was not activated, and macOS has no bare `python`. | Either `source _verify/.venv-verify/bin/activate` first, or call `_verify/.venv-verify/bin/python` by full path. Block 3's script does the second, so activation is not needed. |
| `fatal: not a git repository` | `git grep` run from the L1 folder, which is not a repo — the clones inside it are. | `cd _verify` (or a `retrotransposon-miner` clone) before any `git` command. |
| `zsh: event not found: \s` | A `!` inside double quotes triggers zsh history expansion. | Single-quote any pattern containing `!`, or avoid `!` in `git grep` patterns. |

A general rule that would have prevented three of them: when a command spans more than one
line or contains `&`, `!`, or quotes inside quotes, put it in a file and run the file.

Everything is verified: I ran every block below on your machine's own interpreter
(`_verify/.venv-verify/bin/python`, Python 3.15.0, pandas 3.0.6, numpy 2.5.4, scipy 1.18.1,
twobitreader 4.0.2, pytest 9.1.1) and the expected output is what it actually printed, not
what I think it should print.

Run the blocks in order. After each, compare with "Expect". If it differs, stop and paste it.

Shell setup, once per terminal:

```bash
cd ~/Desktop/Research/Research_Projects/L1
source _verify/.venv-verify/bin/activate
```

Note the path: your L1 folder moved to a single `Desktop`, not the doubled one.

---

## Block 0 — diagnostics (2 min, run these first)

```bash
cd ~/Desktop/Research/Research_Projects/L1
pwd
ls -d vcf ref gencode nested_analysis artifacts _verify _job8
wc -c ref/hg38.2bit
python3 -c "import sys; print(sys.version.split()[0])"
```

Expect: the path ends `/Desktop/Research/Research_Projects/L1` with ONE `Desktop`; all seven
directories listed with no error; `wc -c` prints `835393456`; any python3 version.

```bash
cd _verify
git -c safe.directory='*' status --short
git -c safe.directory='*' branch --show-current
git -c safe.directory='*' fetch origin
git -c safe.directory='*' log --oneline -1 origin/main
cd ..
```

Expect: `status --short` prints nothing. Current branch `analysis/line1-span-reach`. Last line
starts `05fce2a`. If `git` prints nothing at all and returns non-zero, that is the sandbox
config problem, not yours — ignore it, your own shell is fine.

---

## Block 1 — the two PRs that are already pushed and have no PR (15 min)

Nothing to run. Both branches are on the remote; I checked `git ls-remote`. The highest PR
number is 82, and neither of these has one. William asked for the histograms one on Monday.

1. Open https://github.com/L1-Therapeutics/retrotransposon-miner/pull/new/analysis/insertion-size-histograms
   Description: paste `PR_insertion_size_histograms.md`.
2. Open https://github.com/L1-Therapeutics/retrotransposon-miner/pull/new/analysis/line1-span-reach
   Description: paste `PR_line1_span_reach.md`, with one edit first —

   Find: `that is 11.6%; against a 281 bp Alu`
   Replace: `that is 11.6% per junction; against a 281 bp Alu`

   The 11.6% is the reach from one edge. Without "per junction" a reader doubles it by accident.

Then check both appeared:

```bash
cd ~/Desktop/Research/Research_Projects/L1/_verify
git -c safe.directory='*' ls-remote origin 'refs/pull/*/head' | awk '{print $2}' | \
  sed 's|refs/pull/||;s|/head||' | sort -n | tail -3
cd ..
```

Expect: `82`, `83`, `84` — the two new numbers. Before you open them it prints only up to `82`.

---

## Block 2 — PRs #81 and #82 need review, not work (5 min)

I checked both: `136edcd` (#81, the GENE/GENEID register fix) merges into current `main` with
zero conflicts, and `main`'s one newer commit `05fce2a` touches only `scripts/` and two
unrelated tests, so #81 does not need a rebase even though it is behind.

So there is nothing to fix. Message William:

> #81 (gene-annotation register fix) and #82 (gene-orientation analysis) are both ready for
> review. #81 is behind main by one commit but merges cleanly — main's newer commit only
> touches scripts/. Happy to rebase if you'd rather see it on top of main.

Confirm the merge-clean claim yourself if you want:

```bash
cd ~/Desktop/Research/Research_Projects/L1/_verify
git -c safe.directory='*' merge-tree --write-tree origin/main 136edcd >/dev/null 2>&1; echo "rc=$?"
cd ..
```

Expect: `rc=0`. Non-zero means a conflict appeared since I checked.

---

## Block 3 — reproduce the Alu poly(A) analysis (one command, 45 s of compute)

Everything is in `L1/_job8/`. One script does the whole block, calls the venv interpreter by
full path, and fails loudly if an input is missing. You do not need to activate anything.

```bash
bash ~/Desktop/Research/Research_Projects/L1/_job8/run_block3.sh
```

It prints five sections. Expect, in order:

1. `python 3.15.0 | pandas 3.0.6` then `inputs present`.
2. Four `: OK` lines from `shasum -c`.
3. `copied l1_strand_scan.py from artifacts/`.
4. `17 passed`.
5. The analysis table, which must match to the digit:

```
AluY consensus 311 bp from 4679 copies (mean column support 0.9538)
A-rich regions (>=0.6 A over 20 bp): [(112, 139), (274, 311)] = 66 bp

test                                       sense  anti  ratio   mono p  markov p
--------------------------------------------------------------------------------
full_consensus                                 2     0   None   0.2559    0.4473
full_consensus_canonical_only                  0     0   None   1.0000    1.0000
a_rich_regions_removed                         0     0   None   1.0000    1.0000
a_rich_regions_removed_canonical_only          0     0   None   1.0000    1.0000

per-copy test: 67,379 AluY copies inside single-strand genes, overall sense fraction 0.4501
hexamer location     carriers   sense     non  crude OR  adj OR          95% CI        p
----------------------------------------------------------------------------------------
any_position           11,923  0.4295  0.4546     0.903   0.950     0.911-0.989   0.0135
a_tail_only             7,955  0.4265  0.4533     0.897   0.954     0.909-1.002   0.0639
element_body_only       4,560  0.4287  0.4517     0.911   0.933     0.878-0.992   0.0271
```

6. Then the verification section:

```
strand	pos_1based	motif	in_a_rich_region
+	116	ACTAAA	1
+	122	AATACA	1
consensus length 311 | AATAAA 0 | ATTAAA 0
tail TCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
copies 67379 | sense fraction 0.4501
orient     antisense  sense
n_hex_any
False          30248  25208
True            6802   5121
```

`AATAAA 0 | ATTAAA 0` is the first result: the AluY consensus carries no canonical poly(A)
signal, and its tail is `TC` plus exactly 30 A's, which cannot contain one because the motif
needs a T.

When it all matches, tell me "job 8 reproduced". Do not send these numbers to William before
then — see "What NOT to do today".

### Block 3b — the LINE-1 span counts (the snippet that failed in zsh)

```bash
cd ~/Desktop/Research/Research_Projects/L1
_verify/.venv-verify/bin/python _job8/check_spans.py
```

Expect exactly:

```
VCF files                     10
all calls                     16610
LINE-1 non-sentinel           2186
  of those, span >= 5500      336
  of those, MEI_5P <= 100     328
  of those, MEI_3P >= 6000    330
```

Run from the L1 folder, not from `_job8`. From the wrong folder it exits 1 with
`no VCFs found - run this from the L1 folder, not from _job8`.

These counts are consistent with William's two-cluster rule but do not test it: `MEI_5P` and
`MEI_3P` are what that rule outputs, so the VCF cannot separate reads that genuinely reach
consensus position 1 from a span the rule filled in.

## Block 4 — reproduce the exon number so the blog figure is yours (10 min, 3 min of compute)

This one is still outstanding from before. The 0.58 in your technical post and in the PR #82
description is my run, not yours.

```bash
cd ~/Desktop/Research/Research_Projects/L1/_verify
git -c safe.directory='*' checkout analysis/gene-orientation
git -c safe.directory='*' log --oneline -1
```

Expect the last line starts `87b1fbb`.

```bash
cd ~/Desktop/Research/Research_Projects/L1
PYTHONPATH=_verify/src _verify/.venv-verify/bin/python _verify/scripts/exon_depletion_null.py \
  vcf gencode/gencode.v49.annotation.gtf.gz ref/hg38.2bit \
  hg38_chrom_sizes.json hg38_gaps.json --out-dir exon_run
```

(Full path to the interpreter, so this works whether or not the venv is activated — bare
`python` is what failed for you last time.)

Expect, in the printed table, the `ALL` row: `gc_matched_within_genes_250bp` obs_over_exp
`0.577`, and the 50 bp row `0.597`.

```bash
cd _verify
git -c safe.directory='*' checkout analysis/line1-span-reach
cd ..
```

Go back to the span branch so the tree is where you left it.

---

## What NOT to do today

- Do not send William the Alu poly(A) result. Adjusted odds ratio 0.950 against a baseline
  sense fraction of 0.45 is a small effect, half the crude association was age and length
  confounding, and non-carriers are themselves strongly sense-depleted — so this is at most a
  contributing term, not the explanation. It needs the AluSx/AluJo replication first.
- Do not quote the exon 0.58 anywhere until Block 4 matches.
- Do not `git add -A` in any clone. Named files only.
- Do not commit from `retrotransposon-miner` or any `rtm-*` worktree. Only `_verify`.

## Still open from earlier, unchanged

- Rotate the AWS key and Azure password that were pasted into the chat log.
- Move the three employment PDFs out of the top level of the L1 folder.
