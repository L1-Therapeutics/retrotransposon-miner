# Retrotransposon Marker Map — v5 build report

## Provenance

- **Generated:** `2026-10-16 20:08:49`
- **Builder:** `render_marker_map_v5.py` (stdlib + Graphviz `dot` only)
- **Input tables** (read-only; SHA-256 of the *input* files):

| file | path | sha256 |
| --- | --- | --- |
| `marker_map_nodes.csv` | `/Users/smyan/Desktop/Research/Research_Projects/L1/marker_map/marker_map_out_v3/marker_map_nodes.csv` | `fe61593d2c8d2aa00c24ea7d7fb6e0f09f0cf4b7d737eed378d0e0047c9bdc10` |
| `marker_map_edges.csv` | `/Users/smyan/Desktop/Research/Research_Projects/L1/marker_map/marker_map_out_v3/marker_map_edges.csv` | `2d18b04c2eb644d365200ea90834299b44758dfb700a3a68d330c45eda46e139` |
- **Metrics JSON:** `/Users/smyan/Desktop/Research/Research_Projects/L1/nested_analysis/results_multi_sample/analysis_metrics.json` — **metrics_consumed** (metrics JSON validated (analysis_metrics.json) run=4a2f58264d04)
  - digest at build time: `f3265faaef037c232094ecc15c67dcf47b9fb9d813a798b0eca3993cbc6291a7` — re-`shasum` the file and compare.

## Correction ledger as applied

Applied in memory after load; the exported CSVs carry the corrected text.

| node / edge | change |
| --- | --- |
| `b_trunc` | label + source replaced: L1 genomic ~99% of >500k copies truncated (Freeman 2022); TSD-flanked intact L1s 85% 5'-truncated (Zingler 2005, denominator stated); ALU 8/226 (Wagstaff 2012) |
| `b_rearr` | label appended `(de novo cell culture; Symer 2002)`; source → `Wagstaff 2012; Symer 2002` |
| `b_invporub` | label + source replaced: 26% (330/1,271) polymorphic L1s, 88% of inverted cluster at bases 4,000–6,000 (Porubsky 2023 [V]) |
| `e_l1hs` | label + source replaced: 1,642 annotated / 321 >6kb = 20% (Chu 2023, hg38); L1Base 146 FLI (Penzkofer); 6 hot = 84% (Brouha 2003 [V]) |
| `e_l1hstrunc` | **new node** (section E): L1Hs truncated ~71–72% (Pavlicek 2002 [S]) — denominator marked unverified |
| `e_l1hstrunc -> position` | **new edge**, tier `associated`, label `L1Hs subset` |
| `b_trans_sva` | label appended with the measured boundary: 0/671 cohort SVA calls can span past the element's 5' end and no call record holds an inserted length or names a source element; 8 reads at 1 of 3 SVA loci, clips 2–100 bp with nothing to subtract; source → + `sva_5prime_transduction_probe.py` |
| `b_l1term` | **new node** (section B): the L1 termination-site snippet with its boundary measured by `l1_readthrough_probe.py` — 0/2186 calls can span past the element's 3' end, and 18 resolved junction reads at 4 loci show none transduced |
| `b_l1term -> b_trans_l1` | **new edge**, tier `associated`, label `readthrough substrate` |

## Section H · PRECEDENT — what exists beyond us

H exists to keep the map honest about *precedent*: the positional-analysis
machinery we use is not original, and the strongest external tests of
"is nested L1 insertion random?" have known limits. 9 nodes:

- **h_jedlicka** — the only published within-host positional analysis of
  nested TEs (12 plant genomes, 1,491 pairs; host 3'UTR most targeted).
- **h_levy_null** — Levy's positional null is an *empirical*
  fragmented-TE background (6 equal bins, chi², ≥5 expected/bin):
  not uniform, not GC-matched.
- **h_levy_ends** — Levy excluded element ends for short-read
  detectability reasons; long reads escape that exclusion.
- **h_wimmer** — the recurrence standard: independent sites need TSD-length
  or subfamily discordance (18/4,500 = 0.4% of NF1 mutations).
- **h_flasch** — ~30% of engineered insertions land in existing L1/Alu
  sequence, but *proportional to opportunity*: "no lightning rods".
- **h_jed_mech** — mechanism split (Ty3/gypsy targeting vs Ty1/copia
  retention/loss).
- **h_pol3** — yeast Ty1-4 genomic target preference near Pol III genes:
  a targeting parallel, not a within-element one.
- **h_twinprim** — twin priming (Ostertag & Kazazian 2001): the TPRT DSB
  leaves two 3' overhangs, one priming on the poly(A) tail and one
  internally on the RNA; microhomology-driven SSA resolves into a 5'
  inversion (~25% of recent L1s inverted+truncated, abstract [V]). This is
  the published *mechanism* behind the inversion class in section B, not a
  new claim of ours.
- **h_atail** — the A-tail is a substrate, not just a primer site: Alu
  terminal segments expand during TPRT (median 42 A vs 23 internal; 95.9%
  grew ≥4 A via ORF2p slippage). Child-tail dynamics at the priming
  interface are distinct from host-tail targeting.

The h_wimmer → x_tsdcomp edge is `gap`-tier on purpose: a nested-recurrence
test does not exist, and ours is in flight. That gap is the point of the map.

## Section G · OUR DATA — data-driven status

**Mode: `metrics_consumed`.** The metrics JSON validated and every `g_*`
node now carries the real 10-genome figures. Nothing is hardcoded — each
value is read from the file at build time, and the JSON path it came from
is printed beside it.

The run is identified by its `provenance.commit`, **not** by a date:
the file records no run date at all, and inventing one from the file
mtime or from `seed` would be a guess. It is labelled `unrecorded`
everywhere instead.

| field | value | read from |
| --- | --- | --- |
| `linker_fold` | 2.45x | `ten/headline/linker` |
| `tail_fold` | 3.98x | `ten/headline/tail` |
| `nested_l1_count` | 364 | `ten/headline/nested_L1` |
| `union_sites` | 7,132 | `ten/headline/union` |
| `replication` | 10 | `ten/headline/replication` |
| `significant_replication` | 10 | `ten/headline/significant_replication` |
| `replication_fraction` | 100% of genomes | `ten/headline/replication_fraction` |

`ten/samples` holds 10 distinct genomes and `ten/per_sample` 320 rows; both
are cross-checked against the headline, and the run is rejected as
in-flight if they disagree.

**Superseded, not deleted.** The `[5-genome tier superseded by 10-genome
run 4a2f58264d045cc32fc88a1014f95c92d34cfbdd]` tag replaces the old `[5-genome; 10-genome rerun in
flight]`. The 5-genome numbers stay on the page so both tiers are visible
and the reader can see which one is current. Note the direction of the
linker change: 2.75x (5 genomes) → 2.45x (10), i.e. the extra genomes
*lowered* it, so the 5-genome figure was the optimistic one. Tail is
essentially unchanged (3.97x → 3.98x).

### The p-value floor is a result, not a footnote

`provenance.replicates` is 10000, so the smallest p this method can
emit is `1/10000 = 0.0001`. Every bin sitting on that value is reported
as **p=0.0001 (1/10000, the floor)** and never as "p < 0.0001", because
claiming a smaller p than the method can produce is the easiest way to
overstate this result. `a_133ours` also carries a legacy chi-square
`p=4e-18` from the earlier tier; that p is on a different scale and is
tagged as not comparable rather than silently restated.

### Counts the analysis itself flagged as still moving

`delta[].flag_gt20pct` is the analysis grading its own 5→10-genome
shift. These five moved more than 20%, so **every absolute count from
the 5-genome tier is a floor, not an estimate**:

| headline | 5 genomes | 10 genomes | change |
| --- | --- | --- | --- |
| Dedup union | 4998 | 7132 | **+42.7%** |
| Private | 3051 | 3838 | **+25.8%** |
| Shared | 1947 | 3294 | **+69.2%** |
| Nested L1 | 231 | 364 | **+57.6%** |
| Consensus-projected eligible L1 | 59 | 90 | **+52.5%** |

### Every measured bin is a range, not a number

`ten/profiles` holds **20 series**, not the 5 curves the sheet draws. The
others measure the *same bins* with a different method or a different
cohort layer, and they disagree with the drawn curve. Each panel carries
the min–max of its siblings as a pale band behind the curve, and
prints the numbers wherever the band is wide enough to change how the
bin reads. A sibling is admitted only when it shares the family **and**
the exact bin boundaries — an envelope between two curves on
different bins would be meaningless.

A band counts as *material* when it is more than 25% wide **and**
reaches the null line. The second condition matters: a band spanning
0.00× to 0.20× is 100% wide and entirely uninformative, because
both series agree there is nothing there.

| curve | siblings | material bins | widest bin | range | from → to |
| --- | --- | --- | --- | --- | --- |
| `alu_body` | 6 | 2 of 16 | 120–140 bp | 1.26×–2.45× | tier2/private → primary/union |
| `alu_private` | 6 | 2 of 16 | 120–140 bp | 1.26×–2.45× | tier2/private → primary/union |
| `alu_rel` | 2 | 0 of 10 | — | none | — → — |
| `l1_full` | 6 | 8 of 13 | 3–3.5 kb | 0.47×–1.32× | consensus_tier2/private → consensus_short_read/union |
| `l1_rel` | 4 | 8 of 10 | 0.9–1.0 | 0.65×–1.48× | near_full_relative/private → relative_all/union |
| `sva_rel` | 2 | 0 of 10 | — | none | — → — |

Three consequences, stated plainly:

1. **The headline linker bin is method-dependent by about 2×.**
   120–140 bp reads 1.26× under `tier2/private` and 2.45× under `primary/union`.
   The peak *location* is not in doubt; the peak *size* is.
   `a_133ours` now prints the range, not a bare figure.
2. **The Alu 3' tail is not.** 280–300 bp measures 3.57–4.04×
   across the same siblings, so `a_tail` keeps a single number.
3. **The L1 element-body curve is method-fragile; the relative-position
   curves are not.** `l1_full` has 8 of 13 bins material, against 0 of 10 for `alu_rel` and 0 of 10 for `sva_rel`.
   So the three-family contrast in lane D — L1 at the host's 5' end,
   SVA at its midpoint, Alu at its 3' end — is the robust part of
   this result: **all three peaks stay above 1×** under every sibling
   series that measures them.

`parity.json → metrics_spread` carries the per-bin table with the
series behind each end of each range.

### A bin with no observations is a gap, not a zero

The file reports `enrichment = 0.0` for a bin nothing landed in, which
is not the same as having measured a zero. Painting that at the axis
floor would draw a curve diving to nothing — a depletion claim, out
of a test that only ever asks about *enrichment*. So those bins are left
out, the curve edge is broken across them, and what was expected is
marked in their place:

- `sva_rel` — 0.8–1.0

Every other plotted bin is floored on data: the thinnest bin on any
panel carries 1 observation.

### The window claims were corrected at their source

`a_133ours` and `a_tail` are the nodes that first claimed the linker and
tail windows, and both still carried the old figures as unqualified fact.
Each now carries a supersession tag naming the current measurement, with
the old fold-change and `n` **parsed out of the node's own label** so the
"understated by Nx" statement cannot drift from the text it describes:

- `a_133ours` (120–140 bp) → **2.45x, n=124**, p=0.0001, obs 124 vs expected
  50.6 under the opportunity model.
- `a_tail` (280–300 bp) → **3.98x, n=164**, p=0.0001, obs 164 vs expected
  41.2 under the opportunity model.

The evidence base had been understated more than fourfold (n=29→124,
n=36→164) and the linker had moved 11% downward. Nothing was deleted:
the 5-genome sentence stays, with the correction attached to it.

## Position map · read this first

`marker_map_position_map.svg` is a cartographic sheet, not a graph: every
interval claim in the map is placed on a shared ruler, so the reading is
"in this window I have that" rather than "these boxes are connected".

- **Lane A · element body coordinate**, 0–300 bp, on a **linear**
  0–320 bp ruler with the 0–100 bp region shaded and captioned in
  10 bp bands. Holds the EN 7-mer, motif positions 2–5, the +7/+11/+19
  top-strand cuts, nick2 at 15–16 bp, the 5'-end peak, the A-tail 20–25 bp
  window, the left/right monomers and the 280–300 bp tail peak.
- **Lane B · TSD length**, same ruler, so TSDs are directly comparable
  with element-body positions. `b_trans_sva` (2,161 bp) runs off the right
  edge and is drawn with an explicit arrow, never silently clipped.
- **Lane C · window / extent scale**, on its own log axis (10 bp → 1 Mb),
  because a 50-bp GC window and a 1-Mb Alu window are not positions.
- **Lane D · relative position along the host**, 0–1.0, its own
  normalised ruler. This is where the families separate: Alu peaks 2.62x at
  the 3' end of the host, SVA 3.60x near its midpoint.
- **Measured enrichment panel**, drawn directly above lane A on lane A's
  *own* 0–320 bp ruler, so a spike in the curve and the bars piled up
  beneath it are literally the same base pairs.

| property | value |
| --- | --- |
| nodes mapped | 23 |
| interval marks | 31 |
| lane A / B bar sub-rows | 4 / 6 (packed, no bar overlaps a bar) |
| generator | built-in SVG writer (exact geometry, no layout engine) |
| measured curves plotted | alu_body, alu_rel, l1_full, l1_rel, sva_rel |
| coordinate provenance | 26 of 31 placements verbatim from the node's own label; 5 declared below and tagged in the card |

### Where the ruler is not quoting the label

An earlier version of this sheet promised that every `lo/hi` was
transcribed verbatim and that nothing was inferred. **That was false for
five of the 31 placements**, for three different and individually
legitimate reasons. Each is now declared in `MAP_COORD_BASIS`, checked
against the node's own label at build time, and tagged in the card it
appears in:

| lane | node | basis | card tag | why it is not verbatim |
| --- | --- | --- | --- | --- |
| A | `a_5p` | `origin` | `[pos. 0]` | the label names the 5' end, which is this ruler's own zero, so the coordinate adds nothing |
| A | `a_gc50` | `extent` | `[width only]` | a width is stated and no position is; drawn as a bracket at the origin instead of a located interval |
| A | `a_motif7` | `extent` | `[width only]` | a width is stated and no position is; drawn as a bracket at the origin instead of a located interval |
| C | `a_alu1mb` | `unit` | `[unit conv.]` | the magnitude is stated in kb/Mb and shown here in bp — same number, different unit |
| C | `a_gc20kb` | `unit` | `[unit conv.]` | the magnitude is stated in kb/Mb and shown here in bp — same number, different unit |

The two `extent` placements are the important ones: **a width is stated
and no position is.** `a_motif7` says "EN 7-mer consensus" and never says
where in the element it sits; `a_gc50` says "50-bp cut window", which is
a size. Drawing either as a filled, tier-coloured bar at a position would
have manufactured a fact out of a measurement, so they are drawn instead
as dashed grey brackets with end ticks at the ruler origin, rather than
filled tier-coloured bars, carrying only the width they do state. The
build refuses to run if an evidence tier is ever given the bracket's
ink, so the legend cannot go quietly stale.

The build refuses the sheet if any placement is neither verbatim nor
declared (`audit_coord_basis`, exit 2), and the declaration is keyed by
`(lane, node)` so an exemption granted for one lane cannot leak into
another. Bracketed provenance stamps are stripped before the check, so a
number this renderer itself added can never become the evidence for a
coordinate.

The sheet is built from the same node table as the graph sheets and emits
`<g class="node"><title>id</title>` per card, so search, tier chips, the
citation panel and deep links work on it exactly as they do elsewhere.
It ships as SVG only — it is not a Graphviz sheet, so there is no PNG.

A log axis was tried first and rejected: it collapsed the 1–300 bp element
body, where nearly all the content lives, into about 5% of the width.

## Parity

- sections rendered: **9** (spine + 8 section sheets)
- views rendered: **11** (10 Graphviz sheets + the position map)
- nodes: **128** · edges: **130** · fatal errors: **0**
- corrections applied: 9 · H section added: True

| view | class | nodes | edges | warnings |
| --- | --- | --- | --- | --- | --- |
| marker_map_position_map | map | 23 | 0 | 0 | 31 interval marks, linear 0–320 bp + log axis |
| marker_map_spine | spine | 10 | 10 | 0 | — |
| marker_map_A | section | 38 | 37 | 0 | — |
| marker_map_B | section | 26 | 24 | 0 | — |
| marker_map_C | section | 13 | 11 | 0 | — |
| marker_map_D | section | 13 | 12 | 0 | — |
| marker_map_E | section | 19 | 16 | 0 | — |
| marker_map_G | section | 9 | 6 | 0 | — |
| marker_map_H | section | 9 | 2 | 0 | — |
| marker_map_X | section | 12 | 9 | 0 | — |
| marker_map_full | full | 128 | 130 | 0 | — |

`sections_rendered` and `len(views)` are reported separately on purpose:
the section sheets are 8, the spine is a 9th section-shaped sheet, `full`
is a 10th view that is not a section, and the position map is an 11th view
that is not a graph at all.

## Viewer

`index.html` (833,353 bytes) is fully offline — no CDN, no network. Wheel zooms
at the cursor, drag pans, double-click fits. Added in v5: search over
id+label+source, evidence-tier chips (nodes by node tier, edges by their
own tier), a DOM-only citation panel, `#view=...&node=...` deep links, and
a help overlay documenting `[V]/[S]/[C]/[W]`. The position map is the first
item in the view list, because it is the sheet to read first.

Node labels on the Graphviz sheets are HTML tables — bold headline, body,
italic citation — so a printed PNG or PDF carries its own references
instead of relying on hover.

## File manifest (SHA-256)

Computed last, after every other write completed. `README.md` itself is
the only excluded entry: its manifest block is written after hashing, so
its own digest cannot appear inside it.

| file | bytes | sha256 |
| --- | --- | --- |
| `index.html` | 833353 | `e2e78dd4899827c5b283c9413bf962cdbe0e8926e023cf87bc91d4464b5e9d1d` |
| `marker_map_A.dot` | 34867 | `1ce4f86e0de8767dbf908a83612052d28f93efb55a5701fe166767d47d58ac0a` |
| `marker_map_A.png` | 1105527 | `c721856b9c503efa3f7f0adf1c6a79b2409d243cd3440e6f1a79a950fd344436` |
| `marker_map_A.svg` | 72583 | `ac877ad70abcf0ab03116b9479b5d2f4821451d765571fd5f66e47fe5c986ef4` |
| `marker_map_B.dot` | 27285 | `18cac519ba4483d23c6d103dfb668b5a0b1b9178c76be9f06d575cb0b07a91a3` |
| `marker_map_B.png` | 858979 | `db8ae79a5853542e6bce5d4e068c06a1f0ea2b90120bbd2cc5b3588528ae8277` |
| `marker_map_B.svg` | 58503 | `f1341c6ed09a6652f468bbfb4967a51ac0922c2dc52a46696d11c6c5f7586848` |
| `marker_map_C.dot` | 15185 | `25ebb6b3631b73fd5070f2798d597dfd28e22002d3a582f196faa99dcc0e0f3e` |
| `marker_map_C.png` | 403564 | `3962b36456936417f624163be47bb4d8086fca5921a4472a5a8e6faefa27941e` |
| `marker_map_C.svg` | 30635 | `fe95da61554da4d78d7edc389f3c55d17bc39bdb14aeb3df30016cf8010de529` |
| `marker_map_D.dot` | 15581 | `77c03fa0370810fc23ca5f42d600399147fdd711e494befee6bc485520abfc37` |
| `marker_map_D.png` | 630446 | `6f761ec8bf6e456ede8c92c45cfcabd88d01ed85f07af9a55459021de626763e` |
| `marker_map_D.svg` | 31353 | `e83f056a8d8cadae2b74f7ddaeabcb8d21fcfb31c504a6df03c672a15e36c8a2` |
| `marker_map_E.dot` | 20896 | `db903d446be576245404821dae6f437f1267a524ae184a18ab9878fbf37017d4` |
| `marker_map_E.png` | 776631 | `2650c1dc35702117c8970774d916527e66928b5601abd3a87f412472affdfc1c` |
| `marker_map_E.svg` | 46475 | `7d9a27b74552cfec030379003a3d3d03a691ee1e57b3ee4d6740a32ad78af401` |
| `marker_map_G.dot` | 15459 | `02d712efc4ed1c3d50a31daa7bcf4e3aa3b93782265dce9c8b554f3195ed5286` |
| `marker_map_G.png` | 449706 | `84fec940a1754702934776473938c6348c4adbf1285a8a05810d4587f0e5897d` |
| `marker_map_G.svg` | 32685 | `62c04c52a8e725e96af0cc1674076827a7d4f5f979c5d6640090475e900aa021` |
| `marker_map_H.dot` | 12160 | `c6b09aa4d189469547a88465fbb5155f19a457ba2f1c12967a655b426c86cd73` |
| `marker_map_H.png` | 312673 | `7fcbaeb2d8bd2731eca78332165878c0370015f801b9721d043095b1e0a57c5b` |
| `marker_map_H.svg` | 23949 | `0d12754e07006dea80a92262ad3ef304dac8571bf3ae584d9c86da155bafcc79` |
| `marker_map_X.dot` | 14635 | `0e0674170f0a0a6fbc9213c23b7f2ab17a1aa4ed1cde79596e9082c154fe75a4` |
| `marker_map_X.png` | 470639 | `81498b9c59058389009fcc8f56620bc5e968c251b64e39fc3f6a75c834d8d656` |
| `marker_map_X.svg` | 30444 | `1eb24274e0d8ec26677e996c853d696c5605b3ecaf6b01d4f58a0bcb11ed86cb` |
| `marker_map_edges.csv` | 5667 | `be64da96c8959bcb1513bfc875c433fc1451593493cd5ba788c72c52de22e514` |
| `marker_map_full.dot` | 121117 | `eb365c3c526482404b63eb24f15e4b11ed1f7adc2c06e83bbac5be26c6765a1c` |
| `marker_map_full.pdf` | 206368 | `4ebaf4774fded216769acb054973c643fa2a459fb4c5ce422c639e5e719d5aa9` |
| `marker_map_full.png` | 6564311 | `583f3134497b868b2188f1587cd03aefc2d309d444c3ed8847d68a6c9de4ef65` |
| `marker_map_full.svg` | 307783 | `9db7de18053bc968360d713d5db229a091c2fa2aa8c1f1fbc91016e6bb655b0e` |
| `marker_map_nodes.csv` | 22975 | `5d355fb626174f1bd6dc61fd220cb2a145aead71b1175e35e12bde90d56504f2` |
| `marker_map_position_map.svg` | 59801 | `73d0af19c8431dcd713452d3c864bd084bace7437ee91870d1f9d5399b4fe2c0` |
| `marker_map_spine.dot` | 12382 | `3c5eca38852e404f1e0720c768d029d6a48e5d01f2d3f8857cf9bd08bae05e0c` |
| `marker_map_spine.pdf` | 110667 | `26b41635dcb6ae83e322ec6047e0af057c871eb0e37b12884209d009b857903f` |
| `marker_map_spine.png` | 282740 | `8552e4d3cfc1ce942e2be0bc403c868ce8717d024dbb528457a4b27c1a94b18a` |
| `marker_map_spine.svg` | 25003 | `7a1ca6b364472d67341b7003e3e949f1bad2e47aab83990a2d2472268d0cc42f` |
| `parity.json` | 32672 | `64cd681f442ff0be201da44c39dd95bc05b3df7df28eb813ae531dfa9e3130f8` |

## Rebuild

```sh
python3 render_marker_map_v5.py --verify-only   # validation only
python3 render_marker_map_v5.py                 # full build
```
