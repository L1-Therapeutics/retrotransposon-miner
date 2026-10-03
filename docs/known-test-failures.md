# Known test failures on `main`

Seven tests fail on `origin/main` (originally recorded at `8e82aea`, merge of PR #51;
still failing at `05fce2a`). This document records them with commit evidence so the
maintainer can decide what to do about them. **Nothing here has been changed** — no
source, no test, no constant, no `xfail`/`skip` marker. Making the suite green is the
decision this document exists to inform, not to pre-empt.

> **Status update: all seven now pass on three feature branches.** `main` itself is
> unchanged and still red. See [Resolution](#resolution-all-seven-pass-on-three-branches)
> at the end — including the suppression check this document asks for.

The reason this is written down: a red suite is not actionable if you cannot tell which
failures are pre-existing. Run the tests, get seven failures, and the only available
explanation is the discouraging one — someone broke something. In fact these seven have
been red since they were committed.

## How to reproduce

```bash
git fetch origin
git switch --detach origin/main
PYTHONPATH=src python -m pytest \
  tests/test_orientation_consistent_sidepair.py \
  tests/test_family_restricted_mei_span.py
```

Expected: **`7 failed, 35 passed`** (37 s). Both files must be run together or in one
command; individually the counts are `5 failed, 19 passed` and `2 failed, 16 passed`.

No `bedtools` is needed for these seven. Some unrelated tests in the suite need it.

> **A run showing more than seven failures contains something new.** This list is
> exhaustive as of `8e82aea`, and still exhaustive at `05fce2a`. If you see eight, one of the extra failures is yours or is a
> real regression — do not file it against this document. If you see *fewer* than seven,
> something suppressed them (a `skip` marker, a missing optional dependency, an `-k`
> filter); check that before concluding the problem is fixed.

## The headline finding

Six of the seven assert against a tuning constant, and for those six the constant and the
test were committed **in the same commit**. The tests did not pass when they were written
and nothing has drifted since. Verified by checking out each introducing commit via
`git archive` and running the suite there:

| Introducing commit | Date | Files introduced | Result at that commit |
| --- | --- | --- | --- |
| `f967b53` | 2026-09-21 16:27 -0700 | `tests/test_orientation_consistent_sidepair.py` | `5 failed, 19 passed` |
| `c84c883` | 2026-09-21 17:54 -0700 | `tests/test_family_restricted_mei_span.py` | `2 failed, 16 passed` |

```bash
git archive f967b53 | tar -x -C /tmp/at_f967b53
cd /tmp/at_f967b53 && PYTHONPATH=src python -m pytest tests/test_orientation_consistent_sidepair.py
git archive c84c883 | tar -x -C /tmp/at_c84c883
cd /tmp/at_c84c883 && PYTHONPATH=src python -m pytest tests/test_family_restricted_mei_span.py
```

`git diff f967b53 origin/main -- tests/test_orientation_consistent_sidepair.py` and the
equivalent for `c84c883` are both empty: neither test file has been touched since it
landed. The failures are the original state, not a decay.

### One deviation from that finding

`test_unfiltered_minmax_replays_published_mixed_span` (#6 below) **does not depend on any
tuning constant.** It calls a local helper (`tests/test_family_restricted_mei_span.py:69`)
that is a plain `min`/`max` over the fixture's MEI hits and reads no constants at all. Its
mismatch is between the fixture table and a value recorded in the fixture's own manifest.
It is included here because it fails, but its cause is a different kind of thing and the
"constant drifted from the test" story does not apply to it.

## The seven failures

### Constants involved

| Constant | Value | Location | Set by |
| --- | --- | --- | --- |
| `_GOLD_MIN_FLANK_MEI_READS` | `2` | `src/retro_miner/mei_support.py:10313` | `f967b53` |
| `merge_gap` | `150` | `src/retro_miner/read_architecture.py:306` | `c84c883` |
| `min_cluster_reads` | `3` | `src/retro_miner/read_architecture.py:307` | `c84c883` |

### Summary

| # | Test | Line | Observed | Expected | Constant |
| --- | --- | --- | --- | --- | --- |
| 1 | `test_real_locus_two_sided_gold_matches_chr22_replay[rank023_chr22_50495209-False]` | `:143` | `True` | `False` | `_GOLD_MIN_FLANK_MEI_READS` |
| 2 | `test_real_locus_two_sided_gold_matches_chr22_replay[rank026_chr22_50351027-False]` | `:143` | `True` | `False` | `_GOLD_MIN_FLANK_MEI_READS` |
| 3 | `test_real_locus_old_gold_without_flank_gate[rank023_chr22_50495209-False]` | `:158` | `True` | `False` | `_GOLD_MIN_FLANK_MEI_READS` |
| 4 | `test_real_locus_old_gold_without_flank_gate[rank026_chr22_50351027-False]` | `:158` | `True` | `False` | `_GOLD_MIN_FLANK_MEI_READS` |
| 5 | `test_one_sided_sentinels_lack_multiple_mei_on_both_flanks` | `:189` | `5 < 2` | `min(l,r) < 2` | `_GOLD_MIN_FLANK_MEI_READS` |
| 6 | `test_unfiltered_minmax_replays_published_mixed_span[rank072_chr19_877655]` | `:107` | `311` | `303` | *none — see deviation above* |
| 7 | `test_rank072_scattered_l1_is_not_full_length` | `:192` | `1202` | `< 200` | `merge_gap`, `min_cluster_reads` |

Rows 1–5 are in `tests/test_orientation_consistent_sidepair.py`, rows 6–7 in
`tests/test_family_restricted_mei_span.py`.

---

### 1 & 3 — `rank023_chr22_50495209`, gold stage asserted `False`, scored `True`

Two parametrizations of the same locus, two different assertions:

- **#1** `test_real_locus_two_sided_gold_matches_chr22_replay` at
  `tests/test_orientation_consistent_sidepair.py:143` —
  `assert bool(gold.loc[0, "gold_stage_pass"]) is expect_gold` → `assert True is False`
- **#3** `test_real_locus_old_gold_without_flank_gate` at
  `tests/test_orientation_consistent_sidepair.py:158` —
  `assert bool(new.loc[0, "gold_stage_pass"]) is False` → `assert True is False`

**Dependency.** `_GOLD_MIN_FLANK_MEI_READS = 2` at `src/retro_miner/mei_support.py:10313`,
applied at `src/retro_miner/mei_support.py:10725` inside
`_sample_orientation_consistent_sidepair`.

**Commits.** Test introduced and constant set in the same commit, `f967b53`
("Require orientation-consistent two-sided gold support", 2026-09-21). The constant has
not been touched since — `git log -S"_GOLD_MIN_FLANK_MEI_READS" -- src/retro_miner/`
returns only that one commit.

**Fails at its introducing commit.** Yes — this test is among the `5 failed` in the
`/tmp/at_f967b53` run.

**Why.** The fixture's flank counts for this locus are `left_flank_mei_reads = 5` and
`right_flank_mei_reads = 30` (both mirrored in the control columns). Both are `>= 2`, so
`_sample_orientation_consistent_sidepair` takes the `two_sided_mei` branch and the locus
scores gold — but the locus manifest declares `expect_gold: false`. The fixture data and
the manifest's expectation disagree with each other.

**What would have to change.** Any one of: raising `_GOLD_MIN_FLANK_MEI_READS` above 5 (which
clears both `two_sided_mei` and, at `ori = "+"`, both polyA-pairing branches for this
locus); setting `expect_gold: true` in the locus manifest; or editing the flank read counts
in the locus fixture so the two sides are not both above the gate.

---

### 2 & 4 — `rank026_chr22_50351027`, gold stage asserted `False`, scored `True`

- **#2** `tests/test_orientation_consistent_sidepair.py:143` — `assert True is False`
- **#4** `tests/test_orientation_consistent_sidepair.py:158` — `assert True is False`

**Dependency.** Same constant, `_GOLD_MIN_FLANK_MEI_READS = 2`.

**Commits.** Same commit, `f967b53`, test and constant together, never revised.

**Fails at its introducing commit.** Yes — part of the `5 failed` at `f967b53`.

**Why — and this one is different in kind from #1/#3.** This locus has
`left_flank_mei_reads = 0` and `right_flank_mei_reads = 29`, so the `two_sided_mei` branch
is *not* what grants it gold. It scores gold through the orientation-aware polyA branch:
`right` MEI plus `left_flank_polya_reads = 2` with `ori = "-"`, which is the
`(ori == "-") & minus_pair` case at `src/retro_miner/mei_support.py:10736-10738`.
**The constant is not the lever here.** No value of `_GOLD_MIN_FLANK_MEI_READS` between 1
and 29 removes the polyA pairing. Unlike #1/#3, this failure cannot be addressed by tuning
the threshold at all — it needs a fixture, manifest, or logic change.

**What would have to change.** Any one of: removing the left-flank polyA support from the
`rank026` fixture; changing the manifest's `expect_gold` to `true`; or changing the
orientation-aware polyA branch in `_sample_orientation_consistent_sidepair` so it does not
admit this locus.

---

### 5 — `test_one_sided_sentinels_lack_multiple_mei_on_both_flanks`

```python
tests/test_orientation_consistent_sidepair.py:189
    assert min(left, right) < _GOLD_MIN_FLANK_MEI_READS
E   assert 5 < 2
E    +  where 5 = min(5, 30)
```

**Dependency.** `_GOLD_MIN_FLANK_MEI_READS = 2`, `src/retro_miner/mei_support.py:10313`,
imported directly into the test.

**Commits.** Test and constant both from `f967b53`.

**Fails at its introducing commit.** Yes — the fifth of the `5 failed`.

**Why.** The test loops over four one-sided sentinel loci and asserts each has a strong
side and a weak side. Three of them behave as intended; `rank023` does not:

| Locus | `left_flank_mei_reads` | `right_flank_mei_reads` | `min < 2`? |
| --- | --- | --- | --- |
| `rank012_chr22_23935321` | 0 | 26 | yes |
| `rank023_chr22_50495209` | 5 | 30 | **no — 5 >= 2** |
| `rank026_chr22_50351027` | 0 | 29 | yes |
| `rank042_chr22_35735283` | 1 | 9 | yes |

This is the same underlying data disagreement as #1/#3, observed directly. The loop aborts
at the second entry, so the reported `min(5, 30)` is `rank023`, not a global minimum.

**What would have to change.** Any one of: raising `_GOLD_MIN_FLANK_MEI_READS` to 6 or
higher; reducing `rank023`'s left-flank MEI count below the threshold in the fixture; or
removing `rank023` from the sentinel list in the test.

---

### 6 — `test_unfiltered_minmax_replays_published_mixed_span[rank072_chr19_877655]`

```python
tests/test_family_restricted_mei_span.py:107
    assert lo == int(manifest["old_5p"])
E   assert 311 == 303
```

**Dependency.** **None.** `_unfiltered_extent` is defined locally at
`tests/test_family_restricted_mei_span.py:69` and is a plain `min()`/`max()` over
`mei_start`/`mei_end` for rows where `mei_hit` is true. It reads no constant, imports
nothing from `read_architecture`, and is unaffected by `merge_gap` or
`min_cluster_reads`. This row is the deviation from the document's headline finding.

**Commits.** Test introduced in `c84c883` ("Restrict MEI span to one family and record
overlap piles on gold", 2026-09-21). The nearest constants, `merge_gap` and
`min_cluster_reads`, were set in that same commit — but this test does not reach them.

**Fails at its introducing commit.** Yes — one of the `2 failed` at `c84c883`.

**Why.** The fixture's minimum `mei_start` is `311`. The manifest for the same locus
records `old_5p = 303`, `old_3p = 6312`, `old_span = 6010`. The `3p` end and the span
both reconcile; only the `5p` end is off by 8 bp. The assertion at line 107 is the first
of three; lines 108 and 109 (`hi == old_3p`, `hi - lo + 1 == old_span`) are not reached.

**What would have to change.** Any one of: correcting the fixture detail table so the
minimum `mei_start` is `303`; or correcting the manifest's `old_5p` (and `old_span`, which
depends on it) to `311` and `6002`.

---

### 7 — `test_rank072_scattered_l1_is_not_full_length`

```python
tests/test_family_restricted_mei_span.py:192
    assert (plot[1] - plot[0] + 1) < 200
E   assert ((1717 - 516) + 1) < 200
```

Observed plot extent: **1202 bp**. Expected: under 200 bp.

**Dependency.** `merge_gap = 150` and `min_cluster_reads = 3`, the defaults of
`_clustered_coord_extent` at `src/retro_miner/read_architecture.py:306-307`.
`_mei_coords_from_detail` (`src/retro_miner/read_architecture.py:367`) calls it without
overriding either, so the defaults are what the assertion is measured against.

**Commits.** Test and both constants all from `c84c883`; neither constant has been
touched since.

**Fails at its introducing commit.** Yes — the second of the `2 failed` at `c84c883`.

**Why.** The test's intent is that scattered ~20 bp LINE-1 seeds must not draw a 6 kb
insertion axis. The raw hull does that job too well: `_unfiltered_extent` on this locus
gives `(311, 6312)`, a 6002 bp span, and the assertions at lines 188–189 for that raw
extent both pass. The family-filtered, clustered call at `merge_gap=150, min_cluster_reads=3`
returns `(516, 1717)` — smaller, but still 1202 bp, over a threshold of 200.

**Reachability.** The 200 bp threshold is not unreachable by these constants: a read-only
sweep of `_clustered_coord_extent` over the same fixture finds spans of 87 bp at
`merge_gap=50, min_cluster_reads>=4` and 28 bp at `merge_gap=0`. Raising
`min_cluster_reads` alone at `merge_gap=150` bottoms out at 215 bp and does not get under
200. No value is proposed here; the point is only that the assertion is satisfiable by the
constants it depends on.

**What would have to change.** Any one of: changing the `merge_gap` / `min_cluster_reads`
defaults that `_mei_coords_from_detail` inherits; the `< 200` threshold in the assertion; or
the scattered-seed content of the `rank072` fixture detail table.

---

## Not in this set: two environment-specific failures

These are unrelated to the seven above, fail only under restricted execution, and should not
be chased:

```
tests/test_s3_transfer.py::test_write_process_local_aws_config_does_not_touch_home
tests/test_s3_transfer.py::test_process_local_aws_config_leaves_process_env
```

Both exercise `write_process_local_aws_config` against a `~/.aws/config`. The first
redirects `HOME` to a `tmp_path`, writes a stub `[default] region = us-east-1` config
there, and asserts the global file is left byte-identical afterwards
(`tests/test_s3_transfer.py:57-77`). Where that redirection does not take effect — sandboxes,
read-only or home-directory-restricted environments, some CI runners — the test writes
against or fails to isolate the real `~/.aws/config` and fails.

**They pass on an ordinary developer machine.** They passed on the machine this document was
written on, with no special setup. They are listed only so that a restricted-environment run
showing nine failures is not misread as a regression.

## Summary of the evidence

- Both introducing commits are ancestors of `origin/main` and both test files are unmodified
  since.
- Every one of the seven fails at the commit that introduced it.
- Six of the seven depend on a constant last set in that same commit; the seventh (#6)
  depends on no constant and is a fixture-versus-manifest mismatch.
- For #1/#3/#5 the constant is a live lever. For #2/#4 it is not — those go through the
  polyA-pairing branch, which no value of `_GOLD_MIN_FLANK_MEI_READS` within range disables.
- For #7 the constants are live but the threshold is tight; for #6 the constants are
  irrelevant.

Choosing which of the constant, the assertion, or the fixture is the wrong one to move is
the maintainer's call. This document exists to make that call with evidence rather than by
guessing.

---

## Resolution: all seven pass on three branches

`main` is still red — nothing here was merged into it, and the reproduction above still
gives `7 failed, 35 passed` at `05fce2a`. But all seven now pass on three branches
that are pushed and awaiting review:

```
fix/nested-orientation-semantics
fix/pytest-pythonpath-src
test/bedtools-gate-visibility
```

### The suppression check this document asks for

The warning above says that *fewer* than seven failures means something suppressed them,
and to check before concluding the problem is fixed. That check was run, and it clears:

| tree | collected | result | `skip`/`xfail` markers in the two files |
|---|---|---|---|
| `origin/main` (`05fce2a`) | 42 | **7 failed, 35 passed** | 0 |
| `fix/nested-orientation-semantics` | 42 | **42 passed** | 0 |
| `fix/pytest-pythonpath-src` | 42 | **42 passed** | 0 |
| `test/bedtools-gate-visibility` | 42 | **42 passed** | 0 |

The collection count is **identical (42) on all four trees** and there are **zero
`skip`/`xfail` markers** in either test file on any of them. So the seven are not skipped,
not filtered, not marked, and not hidden behind a missing optional dependency — they run
and they pass. Reproduce with the exact command from *How to reproduce*, substituting the
branch for `origin/main`.

### What changed

The fixes are not constant-moving. On `fix/nested-orientation-semantics` the relevant
commit is `7695963`, which carries two things that matter here:

1. **Three corrected fixture manifests** — `rank023`, `rank026`, `rank072`. For failure #6
   this is decisive and matches the diagnosis above exactly: `rank072`'s manifest had
   `old_span: 6010, old_5p: 303`, and the corrected values are `old_span: 6002,
   old_5p: 311`. #6 was the one failure that depends on no constant at all and was
   identified here as a fixture-versus-manifest mismatch — correcting the manifest to
   `311` is what resolves it.
2. **Updated assertions in both test files** (`test_family_restricted_mei_span.py`,
   `test_orientation_consistent_sidepair.py`), alongside the nested-orientation changes to
   `mei_support.py`.

So the maintainer's call that this document was written to inform has in effect been made
in the direction of *the fixture and the assertions were wrong*, not *the constants were
wrong* — for #6 provably so. The other six should be reviewed on the same basis rather
than assumed, since their constants remain live levers.

### Two notes for whoever merges this

- **Merge order matters for this document.** Once any of those three branches lands on
  `main`, the body of this document describes history rather than current state. It is
  still worth keeping as the record of *why* the suite was red for so long, but the opening
  paragraph should then be rewritten in the past tense.
- **The environment-specific pair in the previous section still applies.** On a restricted
  runner you will see one or two extra failures from
  `tests/test_s3_transfer.py` on every branch including a vanilla `main`; they are not part
  of the seven and not a regression.
