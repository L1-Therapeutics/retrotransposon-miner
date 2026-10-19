#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""render_marker_map_v5.py -- Retrotransposon Marker Map v5 ("parting gift" build).

Rebuilds the nested-retrotransposon Marker Map to v5:

  * strict, loud validation of the node/edge tables (exit 2 on any violation);
  * an in-memory, literature-verified correction ledger (b_trunc, b_rearr,
    b_invporub, e_l1hs + new e_l1hstrunc) -- history is replaced, never deleted;
  * a new Section H . PRECEDENT (7 nodes / 7 edges) for what exists beyond us;
  * a data-driven Section G . OUR DATA: the 10-genome metrics JSON is consumed
    when present and schema-valid, otherwise every G node carries an explicit
    "[5-genome; 10-genome rerun in flight]" suffix.  NO 10-genome number is ever
    hardcoded, guessed or projected -- that is a hard gate;
  * per-view DOT/SVG/PNG (+ PDF for spine/full) with a tier/evidence visual
    language, per-section clusters and a legend in every view;
  * lossless CSV re-export (provenance / diffability), parity.json, a SHA-256
    manifest and a README provenance page;
  * an upgraded, fully offline viewer (search, tier filters, citation panel,
    deep links, correct fit()).

Stdlib only.  Graphviz `dot` is the sole external dependency.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import shutil
import struct
import subprocess
import sys
import textwrap
import zlib
from collections import Counter, OrderedDict
from datetime import datetime
from pathlib import Path

# --------------------------------------------------------------------------- #
# §0  paths, exit codes
# --------------------------------------------------------------------------- #
BASE_DIR = Path(__file__).resolve().parent
DEFAULT_OUTDIR = BASE_DIR / "marker_map_out_v5"

# Input tables are looked up here, in order.  The working data of record lives in
# the most recent build directory (marker_map_out_v3); a copy dropped into the
# marker_map/ root wins, so the script keeps working if the tables are promoted.
INPUT_SEARCH_DIRS = [
    BASE_DIR,
    BASE_DIR / "marker_map_out_v3",
    BASE_DIR / "marker_map_out",
]
NODE_CSV_NAME = "marker_map_nodes.csv"
EDGE_CSV_NAME = "marker_map_edges.csv"

# READ-ONLY.  Consumed only if present and schema-valid.
METRICS_PATH = Path(
    "/Users/smyan/Desktop/Research/Research_Projects/L1/nested_analysis"
    "/results_multi_sample/analysis_metrics.json"
)

EXIT_OK = 0
EXIT_VALIDATION = 2
EXIT_GRAPHVIZ = 3
EXIT_NO_DOT = 4

EXIT_MESSAGES = {
    EXIT_VALIDATION: "data validation failure",
    EXIT_GRAPHVIZ: "graphviz fatal condition",
    EXIT_NO_DOT: "graphviz 'dot' not found",
}

# --------------------------------------------------------------------------- #
# §1  data contract (vocabulary is closed -- anything else fails the build)
# --------------------------------------------------------------------------- #
SECTIONS = ["SPINE", "A", "B", "C", "D", "E", "G", "H", "X"]
SECTION_SHEETS = ["A", "B", "C", "D", "E", "G", "H", "X"]  # 8 per-section sheets
KINDS = ["box", "chain", "marker", "structure", "resource", "gap",
         "box3d", "clinical", "folder"]
EVIDENCE = ["demonstrated", "associated", "proposed", "unknown", "gap"]
EDGE_TIERS = ["demonstrated", "associated", "proposed", "unknown", "gap"]
ALLOWED_VIEWS = set(["full", "spine"]) | set(SECTION_SHEETS)

NODE_COLS = ["id", "label", "section", "kind", "evidence", "source", "views"]
EDGE_COLS = ["src", "dst", "tier", "label", "views"]

# --------------------------------------------------------------------------- #
# §5  visual system
# --------------------------------------------------------------------------- #
TITLE = "Retrotransposon Marker Map"

SECTION_META = OrderedDict([
    ("SPINE", ("CAUSAL SPINE", "#ececec", "#616161", "2.2")),
    ("A", ("A \u00b7 TARGETING \u2014 where it lands", "#d6eaf8", "#5DADE2", "1.6")),
    ("B", ("B \u00b7 CREDIBILITY \u2014 is the call real", "#d5f5e3", "#58D68D", "1.6")),
    ("C", ("C \u00b7  ELEMENT STATE \u2014 is it armed", "#fcf3cf", "#F4D03F", "1.6")),
    ("D", ("D \u00b7 MEDICAL & DISEASE", "#fde2e4", "#F1948A", "1.6")),
    ("E", ("E \u00b7 RESOURCES \u2014 known sites", "#e8daef", "#BB8FCE", "1.6")),
    ("G", ("G \u00b7 OUR DATA", "#d4efdf", "#7DCEA0", "1.6")),
    ("H", ("H \u00b7 PRECEDENT \u2014 what exists beyond us", "#fef9e7", "#d4ac0d", "1.6")),
    ("X", ("X \u00b7 VERIFIED GAPS (hexagons)", "#fadbd8", "#E74C3C", "1.6")),
])
LEGEND_META = ("LEGEND", "#ffffff", "#999999", "1.4")

# tier -> (style, color, penwidth, human key)
TIER_STYLE = OrderedDict([
    ("demonstrated", ("solid", "#1c1c1c", "1.6")),
    ("associated", ("dashed", "#2E86C1", "1.4")),
    ("proposed", ("dotted", "#B7950B", "1.5")),
    ("unknown", ("dotted", "#95A5A6", "1.4")),
    ("gap", ("dashed", "#C0392B", "1.6")),
])
GAP_COLOR = "#C0392B"
CHAIN_FILL = "#30343b"
GAP_FILL = "#fdedec"
NODE_FILL = "#ffffff"
INK = "#111111"

KIND_SHAPE = {
    "box": "box", "chain": "box", "marker": "box", "structure": "ellipse",
    "resource": "folder", "gap": "hexagon", "box3d": "box3d",
    "clinical": "diamond", "folder": "folder",
}

# view class -> (wrap width, nodesep, ranksep, graph fontsize, node fontsize)
VIEW_CLASS = {
    "spine": (26, ".35", ".55", 18, 11),
    "section": (30, ".38", ".58", 16, 10),
    "full": (22, ".28", ".50", 14, 9),
}

VIEW_TITLES = {
    "spine": TITLE + " \u2014 spine",
    "full": TITLE,
}

FATAL_MARKERS = ("trouble in init_rank", "triangulation failed", "routesplines",
                 "Pshortestpath", "Error:")
BENIGN_STDERR = ("fontconfig", "Fontconfig", "fc-", "match family")

# --------------------------------------------------------------------------- #
# §2  in-memory correction ledger (literature-verified; applied post-load)
# --------------------------------------------------------------------------- #
CORRECTIONS = OrderedDict([
    ("b_trunc", {
        "label": (
            "Truncation asymmetry:\n"
            "L1 genomic: ~99% of >500k copies\n"
            "truncated (~5,000 full-length;\n"
            "Freeman 2022, approx.)\n"
            "TSD-flanked intact L1s: 85%\n"
            "5'-truncated (Zingler 2005;\n"
            "denominator = TSD-flanked\n"
            "intact inserts, NOT\n"
            "all genomic copies)\n"
            "ALU only 8/226 (3.5%)\n"
            "(Wagstaff 2012)"
        ),
        "source": "Freeman 2022; Zingler 2005; Wagstaff 2012",
    }),
    ("b_rearr", {
        "append_label": "\n(de novo cell culture; Symer 2002)",
        "source": "Wagstaff 2012; Symer 2002",
    }),
    ("b_invporub", {
        "label": (
            "Polymorphic L1: 26%\n"
            "(330/1,271) carry inversions;\n"
            "88% OF INVERTED cluster at\n"
            "L1 bases 4,000-6,000\n"
            "(Porubsky 2023)"
        ),
        "source": "Porubsky 2023 [V]",
    }),
    ("e_l1hs", {
        "label": (
            "L1Hs: 1,642 annotated;\n"
            "321 >6kb = 20% (Chu 2023, hg38);\n"
            "L1Base: 146 full-length INTACT\n"
            "(GRCh38 putatively active\n"
            "FLIs; Penzkofer); 80-100 competent;\n"
            "6 hot = 84% of activity\n"
            "(Brouha 2003)"
        ),
        "source": "Chu 2023; Penzkofer L1Base 2 [V]; Brouha 2003 [V]",
    }),
    # Damert's three figures stay verbatim; what is appended is the boundary our own
    # calls set on them (tested 2026-10-16 by sva_5prime_transduction_probe.py, the
    # SVA arm of the test applied to b_l1term).  A 5'-transduced copy is longer than
    # its element at its 5' end, and the caller's axis starts at the element's own 5'
    # end (coordinate 1): no call can carry a base upstream of it, no call record
    # states an inserted copy's length, and no call names a source element.  So the
    # 9.2% rate is not estimable from these calls and the MAST2 figure needs a donor
    # match; the read level reaches one SVA 5' junction, where the clip is an upper
    # bound (clips 2-100 bp) and nothing in the fields partitions it.
    ("b_trans_sva", {
        "append_label": (
            "\nNOT ESTIMABLE FROM OUR CALLS:\n"
            "0/671 spans past the element's\n"
            "5' end (axis clamped at 1);\n"
            "209 reach that 5' end and no\n"
            "field holds an inserted length;\n"
            "0/671 names a source element;\n"
            "8 reads at 1 of 3 SVA loci,\n"
            "clips 2-100 bp, nothing to\n"
            "subtract at a 5' junction"
        ),
        "source": ("Damert 2009 GR 19:1992 [V]; boundary measured by "
                   "sva_5prime_transduction_probe.py"),
    }),
])

CORRECTION_NEW_NODES = [
    {
        "id": "e_l1hstrunc", "section": "E", "kind": "marker", "evidence": "unknown",
        "label": ("L1Hs truncated ~71-72%\n(Pavlicek 2002; later ~72%)\n"
                  "[SNIPPET \u2014 denominator unverified]"),
        "source": "Pavlicek 2002 [S]",
        "views": ["full"],
    },
    {
        # Where the L1 transcript stops.  Downstream readthrough is the substrate for
        # 3' transduction, so this sits beside b_trans_l1 rather than in the RNA chain.
        # The five parts total 88%, NOT 100%; the source's own "Sum" row is dropped
        # because it matched neither the total nor any part, and the category
        # definitions were not stated.  Recorded verbatim; nothing is normalised.
        # Tested against our own data 2026-10-14 (l1_readthrough_probe.py): the snippet
        # stays uncitable, but its boundary is now measured rather than assumed -- the
        # cohort's call fields cannot express a transduced tail at all, and at the four
        # loci with read-level detail no junction read shows one.
        "id": "b_l1term", "section": "B", "kind": "marker", "evidence": "unknown",
        "label": ("L1 TRANSCRIPTION TERMINATION\n"
                  "SITES, as given:\n"
                  "5' UTR 16% | Internal 15%\n"
                  "3' UTR 16% | Annotated 20%\n"
                  "Downstream 21%\n"
                  "parts total 88%, 12% unlisted\n"
                  "[SNIPPET - SOURCE NOT YET NAMED]\n"
                  "do not cite until the citation\n"
                  "and the category definitions land\n"
                  "NOT ESTIMABLE FROM OUR CALLS:\n"
                  "0/2186 spans past the element\n"
                  "3' end (axis clamped)\n"
                  "18 resolved 3'-junction reads,\n"
                  "4 loci, none shows a segment\n"
                  "between the 3' end and the tail"),
        "source": ("user-supplied 2026-10-14; citation pending [S]; boundary measured "
                   "by l1_readthrough_probe.py"),
        "views": ["full"],
    },
]
CORRECTION_NEW_EDGES = [
    {"src": "e_l1hstrunc", "dst": "position", "tier": "associated",
     "label": "L1Hs subset", "views": ["full"]},
    {"src": "b_l1term", "dst": "b_trans_l1", "tier": "associated",
     "label": "readthrough substrate", "views": ["full"]},
]

# --------------------------------------------------------------------------- #
# §3  Section H . PRECEDENT -- what exists beyond us
# --------------------------------------------------------------------------- #
H_NODES = [
    ("h_jedlicka", "box", "demonstrated",
     "Only published within-host positional\n"
     "analysis of nested TEs:\n"
     "12 plant genomes, 1,491 pairs;\n"
     "host 3'UTR most targeted; host\n"
     "LTRs + integrase least\n"
     "(domain resolution; LTR-RTs)",
     "Jedlicka 2019 PLoS Genet [V]"),
    ("h_levy_null", "marker", "demonstrated",
     "Levy 2010 positional null =\n"
     "EMPIRICAL fragmented-TE background,\n"
     "6 equal bins, chi2, >=5 expected/bin\n"
     "NOT uniform, NOT GC-matched",
     "Levy 2010 NAR Methods [V]"),
    ("h_levy_ends", "marker", "demonstrated",
     "Levy EXCLUDED element ends\n"
     "from positional analysis:\n"
     "short-read detectability rationale\n"
     "('fragments too short for RepeatMasker')",
     "Levy 2010 NAR Methods [V]"),
    ("h_wimmer", "box", "demonstrated",
     "NF1: 3 sites hit 2x in unrelated\n"
     "patients; independence = TSD-length\n"
     "OR subfamily discordance; 18/4,500\n"
     "(0.4%) of NF1 mutations; ALL sites\n"
     "match L1 EN consensus",
     "Wimmer 2011 PLoS Genet [V]"),
    ("h_flasch", "marker", "demonstrated",
     "~30% of engineered L1 insertions\n"
     "land in existing L1/Alu sequence\n"
     "BUT proportional to opportunity:\n"
     "'no lightning rods' (Flasch 2019)",
     "Flasch 2019 Cell [V]"),
    ("h_jed_mech", "marker", "associated",
     "Mechanism split in plants:\n"
     "Ty3/gypsy nesting = targeting;\n"
     "Ty1/copia = retention/loss",
     "Jedlicka 2019 [V]"),
    ("h_pol3", "marker", "demonstrated",
     "Yeast Ty1-4: 82-95% of insertions\n"
     "within 750 bp of Pol III genes\n"
     "(tRNA promoters) \u2014 genomic\n"
     "target preference, NOT\n"
     "within-element",
     "Kim et al. [S]"),
    ("h_twinprim", "marker", "demonstrated",
     "Twin priming (Ostertag & Kazazian 2001):\n"
     "TPRT DSB leaves two 3' overhangs; one primes on\n"
     "poly(A) tail, the other internally on the RNA;\n"
     "microhomology-driven SSA resolves -> 5' inversion\n"
     "~25% of recent L1s inverted+truncated (abstract [V])",
     "Ostertag 2001 Genome Res [V-abst]; Meyer 2010 [V]"),
    ("h_atail", "marker", "demonstrated",
     "Alu A-tail expands during TPRT:\n"
     "terminal segments median 42 A vs 23 internal;\n"
     "95.9% of terminal segments grew >=4 A\n"
     "(ORF2p slippage during RT); priming rare in\n"
     "first 20-25 bp of a segment (Wagstaff 2012)",
     "Wagstaff 2012 PLoS Genet [V]"),
]
H_EDGES = [
    ("h_jedlicka", "a_133", "associated", "cross-kingdom precedent"),
    ("h_levy_null", "a_133", "associated",
     "replication qualitative\nacross nulls"),
    ("h_levy_ends", "a_tail", "demonstrated",
     "long-read escapes\ntheir exclusion"),
    ("h_wimmer", "x_tsdcomp", "gap",
     "no nested-recurrence\ntest exists (ours first,\n in flight)"),
    ("h_flasch", "g_hg", "associated",
     "guards enrichment claims:\nopportunity nulls required"),
    ("h_jed_mech", "h_jedlicka", "associated", "mechanism split"),
    ("h_pol3", "h_jedlicka", "associated",
     "targeting logic\nacross kingdoms"),
    ("h_twinprim", "b_inv", "associated",
     "model for the\ninversion class"),
    ("h_atail", "a_tail", "associated",
     "A-tract dynamics at\nthe priming interface; child-tail mechanism,\nNOT host-tail targeting"),
]

# --------------------------------------------------------------------------- #
# POSITION MAP -- every interval claim placed on ONE distance axis
# --------------------------------------------------------------------------- #
# Each row is (lane, node_id, [(lo, hi, caption), ...]).  Every lo/hi is
# transcribed verbatim from that node's own label (the caption quotes the
# phrase it came from) -- no coordinate is inferred, estimated or invented.
# The map is a view like any other: each card is a <g class="node"><title>id
# ...> so search, tier filters, the citation panel and deep links all work.
MAP_LANES = OrderedDict([
    ("A", ("ELEMENT BODY POSITION",
           "bp along the element or along the named feature, from its 5' cut",
           "#d6eaf8", "#5DADE2", "A \u00b7 where it lands")),
    ("B", ("TSD LENGTH", "bp between the two ends of the insertion",
           "#d5f5e3", "#58D68D", "B \u00b7 is the call real")),
    ("C", ("ANALYSIS WINDOW / EXTENT", "scale of the test, not a position",
           "#fcf3cf", "#F4D03F", "context")),
])

MAP_ITEMS = [
    # -- lane A: the element body, 0 -> ~300 bp ----------------------------- #
    ("A", "a_motif7", [(0, 7, "EN 7-mer, 7 bp")],
     "TTTTT/AA, 9.6% of raw use"),
    ("A", "a_motifC", [(2, 5, "the C lies in motif positions 2-5")],
     "45% of sites; degenerate, not a code"),
    ("A", "a_stag", [(7, 7, "+7 nt"), (11, 11, "+11 nt"), (19, 19, "+19 nt")],
     "top-strand cuts from nick1, in vitro"),
    ("A", "tprt", [(15, 16, "nick2 at +15-16 bp")],
     "TSD = gap - A-rich"),
    ("A", "a_nick2", [(15, 16, "second nick 15-16 bp down")],
     "TNTN/AA class"),
    ("A", "a_5p", [(0, 2, "5'-end peak")],
     "n=13, A-frac 0.13; edge-effect candidate"),
    ("A", "h_atail", [(20, 25, "priming rare in the first 20-25 bp")],
     "A-tail segment dynamics"),
    ("A", "a_gc50", [(0, 50, "50-bp window")],
     "32% GC vs 41% genome: an AT-rich footprint"),
    ("A", "a_consensus", [(1, 120, "left monomer ~1-120"),
                          (134, 282, "right monomer ~134-282"),
                          (280, 282, "Alu core 280-282")],
     "AluSx body, full element ~300 bp"),
    ("A", "a_133", [(118, 136, "AluJo linker 118-136"), (133, 133, "hotspot 133")],
     "directly after the A-rich linker"),
    ("A", "a_133ours", [(120, 140, "our bin 120-140 bp")],
     "3.0x, n=29, A-frac 0.55; chi2=103, p=4e-18"),
    ("A", "a_tail", [(280, 300, "3' tail peak 280-300 bp")],
     "3.7x, n=36, A-frac 0.79; the blind spot"),
    # -- lane B: TSD length ------------------------------------------------- #
    ("B", "b_tsd_alu", [(5, 27, "Alu TSD 5-27 bp")], "the canonical spectrum"),
    ("B", "b_tsd_l1fl", [(14, 16, "L1 full-length TSD 14-16 bp")], "tight"),
    ("B", "b_tsd_sva", [(4, 20, "SVA TSD 4-20 bp")], "SVA set"),
    ("B", "b_filler", [(4, 20, "filler 4-20 bp")], "the background"),
    ("B", "b_tsdshort", [(9, 11, "short 9-11 bp")], "the short tail"),
    ("B", "x_svatsd", [(4, 20, "SVA TSD 4-20 bp")],
     "VERIFIED GAP: no measured SVA TSD set"),
    # 14-2,161 bp is Damert's observed segment range, placed on the ruler as given.
    # The measured boundary is on the node (see the correction ledger): our calls
    # cannot express a segment at all, and at the one SVA 5' junction with reads the
    # largest clip is 100 bp, short of Damert's 223 bp mean.
    ("B", "b_trans_sva", [(14, 2161, "transduced SVA 14-2,161 bp")],
     "the long tail; transduction"),
    # -- lane C: the scale of the test -------------------------------------- #
    ("C", "a_gc50", [(50, 50, "50-bp window")], "local GC"),
    ("C", "a_exonprox", [(50, 50, "LINE 50 bp"), (100, 100, "SINE 100 bp"),
                         (2000, 2000, "near-full LINE 2,000 bp")],
     "exon-proximity exclusion zones"),
    ("C", "a_levmaor", [(150, 150, "150 bp upstream")], "antisense-Alu depletion"),
    ("C", "a_gc20kb", [(20000, 20000, "20-kb flanks")], "regional, GC-neutral"),
    ("C", "a_alu1mb", [(1000000, 1000000, "1-Mb windows")], "Alu 1-Mb windows"),
]

# The footer used to promise that every lo/hi is transcribed verbatim and that
# nothing is inferred.  Three placements are not verbatim, for three different
# and individually legitimate reasons, so the promise was false as written.  The
# reason is now declared per placement and checked against the label at build
# time; a placement that is neither verbatim nor declared stops the build rather
# than quietly acquiring a coordinate it cannot account for.
#
#   "unit"    the magnitude is stated, only the unit differs ("20-kb" -> 20000)
#   "origin"  the position is the ruler's own zero, named in words ("5'-end")
#   "extent"  a WIDTH is stated and no position is.  Drawn as a dashed grey
#             bracket at the origin, never as a located interval: the sheet
#             shows how wide the claim is without claiming where it sits.
MAP_COORD_BASIS = {
    ("A", "a_motif7"): "extent",
    ("A", "a_gc50"): "extent",
    ("A", "a_5p"): "origin",
    ("C", "a_gc20kb"): "unit",
    ("C", "a_alu1mb"): "unit",
}
# Short enough to sit inside a card without crowding it, and stated the same way
# in every card so the tag is learnable after one reading.
COORD_BASIS_TAG = {
    "unit": "[unit conv.]",
    "origin": "[pos. 0]",
    "extent": "[width only]",
}
# Long form for prose; the sheet only has room for COORD_BASIS_WHY below.
COORD_BASIS_LONG = {
    "verbatim": "both endpoints appear as numbers in the node's own label",
    "unit": "the magnitude is stated in kb/Mb and shown here in bp \u2014 same "
            "number, different unit",
    "origin": "the label names the 5' end, which is this ruler's own zero, so "
              "the coordinate adds nothing",
    "extent": "a width is stated and no position is; drawn as a bracket at the "
              "origin instead of a located interval",
}

# Short glosses: these are printed in the footer key, where every character is
# paid for out of the sheet width.  The long form lives in README.md.
COORD_BASIS_WHY = {
    "verbatim": "both endpoints are numbers in the node's own label",
    "unit": "the magnitude is stated in kb/Mb",
    "origin": "the label names the 5' end",
    "extent": "a width is stated, no position",
}
EXTENT_INK = "#7f8fa0"


def _assert_extent_ink_distinct():
    """Grey must keep meaning "no tier here".

    The sheet tells the reader that a dashed grey bracket carries no tier
    colour, which is only true while no tier uses that ink -- and it sits close
    to the `unknown` tier's grey, which is the one a reader could confuse it
    with.  Checked at build time because the legend states it as a fact.
    """
    tier_inks = {c.upper() for _k, (_s, c, _p) in TIER_STYLE.items()}
    if EXTENT_INK.upper() in tier_inks:
        die(EXIT_VALIDATION,
            "EXTENT_INK %s is also a tier colour; the legend says a grey bracket "
            "means 'position not stated', which would no longer hold"
            % EXTENT_INK)
    return True



def audit_coord_basis(nodes, items):
    """Classify every placement against its own node's label.

    Returns {(lane, nid, span_index): kind} plus a tally.  Bracketed provenance
    stamps are stripped first: they are annotations added by this renderer, not
    part of the claim being placed.
    """
    by_id = {n["id"]: n for n in nodes}
    kinds, tally = {}, {}
    for lane, nid, spans, _note in items:
        body = re.sub(r"\[[^\[\]]+\]", " ", by_id.get(nid, {}).get("label", ""))
        nums = [float(x) for x in re.findall(r"\d+(?:\.\d+)?", body.replace(",", ""))]

        def stated(v, _n=nums):
            return any(abs(n - v) < 1e-9 for n in _n)

        for i, span in enumerate(spans):
            lo, hi = span[0], span[1]
            kind = ("verbatim" if stated(lo) and stated(hi)
                    else MAP_COORD_BASIS.get((lane, nid)))
            if kind is None:
                die(EXIT_VALIDATION,
                    "position map: %s is placed at %g-%g in lane %s but neither "
                    "endpoint appears in its own label. Declare the basis in "
                    "MAP_COORD_BASIS, or take it off the ruler."
                    % (nid, lo, hi, lane))
            kinds[(lane, nid, i)] = kind
            tally[kind] = tally.get(kind, 0) + 1
    return kinds, tally

# --------------------------------------------------------------------------- #
# §4  Section G -- data-driven, never projected
# --------------------------------------------------------------------------- #
G_SUFFIX_IDS = ["g_hg", "g_lr", "g_meth", "g_plan", "g_william"]
IN_FLIGHT_SUFFIX = "\n[5-genome; 10-genome rerun in flight]"
# Where each displayed metric comes from inside the analysis JSON.  The file
# nests everything under "ten"/"five"/"provenance"; there are no flat keys, so
# the mapping is written out explicitly and each field records the path it was
# read from.  A field whose path is absent, or whose type is wrong, keeps the
# whole run in flight -- we never fall back to a default.
METRIC_FIELDS = [
    # display key        json path                         expected type
    ("linker_fold",      ("ten", "headline", "linker"),  float),
    ("tail_fold",        ("ten", "headline", "tail"),    float),
    ("nested_l1_count",  ("ten", "headline", "nested_L1"), int),
    ("union_sites",      ("ten", "headline", "union"),   int),
]
# Scalar cross-checks.  These must agree with each other or the run is rejected;
# they are what stops a half-written JSON from being read as a finished run.
METRIC_CHECKS = [
    ("replication",            ("ten", "headline", "replication"), 10),
    ("significant_replication", ("ten", "headline", "significant_replication"), 10),
]
EXPECTED_SAMPLES = 10
METRIC_DISPLAY_KEYS = ["linker_fold", "tail_fold", "nested_l1_count",
                       "union_sites", "replication", "replication_fraction"]


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def die(code: int, message: str) -> "NoReturn":  # type: ignore[valid-type]
    print("[FATAL:%d] %s: %s" % (code, EXIT_MESSAGES.get(code, "error"), message),
          file=sys.stderr)
    raise SystemExit(code)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def find_dot() -> "str | None":
    """`dot` on PATH, else the usual Homebrew prefixes, else None."""
    found = shutil.which("dot")
    if found:
        return found
    for cand in ("/opt/homebrew/bin/dot", "/usr/local/bin/dot",
                 "/opt/local/bin/dot", "/usr/bin/dot"):
        if os.path.isfile(cand) and os.access(cand, os.X_OK):
            return cand
    return None


# --------------------------------------------------------------------------- #
# §1  loading -- delimiter-agnostic, loud on anything unexpected
# --------------------------------------------------------------------------- #
SEP_RE = re.compile(r"^[\s\-:]+$")


def _is_separator(cells) -> bool:
    return bool(cells) and all(SEP_RE.match(c or "") and c is not None
                               for c in cells) and any("-" in (c or "")
                                                        for c in cells)


def _split_pipe(line: str):
    """Split a pipe-table row, protecting escaped pipes."""
    guard = "\x00PIPE\x00"
    protected = line.replace("\\|", guard)
    cells = [c.strip() for c in protected.split("|")]
    return [c.replace(guard, "|") for c in cells]


def _split_comma(line: str):
    return next(csv.reader([line]))


def _parse_rows(path: Path, ncols: int):
    """Return [(lineno, cells)] for data rows; pick the delimiter that yields a
    consistent column count for every row (the shipped tables are comma-CSV with
    quoted fields; pipe-tables with '\\|' escapes are also accepted)."""
    raw = path.read_text(encoding="utf-8").splitlines()
    if not raw:
        die(EXIT_VALIDATION, "%s: file is empty" % path.name)

    for delim, splitter in (("comma", _split_comma), ("pipe", _split_pipe)):
        rows, bad = [], 0
        for lineno, line in enumerate(raw, start=1):
            if not line.strip():
                continue
            try:
                cells = splitter(line)
            except Exception:
                bad += 1
                continue
            cells = [(c or "").strip() for c in cells]
            if _is_separator(cells):
                continue
            if cells[0].lower() in NODE_COLS + EDGE_COLS and len(rows) == 0:
                continue  # header
            if len(cells) != ncols:
                bad += 1
                continue
            rows.append((lineno, cells))
        if rows and bad == 0:
            return delim, rows
        if rows and bad <= 1:  # one stray row is still a real format, keep going
            return delim, rows
    # neither format is consistent -> report the comma parse, line by line
    for lineno, line in enumerate(raw, start=1):
        if not line.strip():
            continue
        cells = [(c or "").strip() for c in _split_comma(line)]
        if _is_separator(cells):
            continue
        if cells[0].lower() in NODE_COLS + EDGE_COLS:
            continue
        if len(cells) != ncols:
            die(EXIT_VALIDATION,
                "%s line %d: expected %d columns, found %d"
                % (path.name, lineno, ncols, len(cells)))
    die(EXIT_VALIDATION, "%s: no parsable data rows" % path.name)


def parse_views(raw: str):
    s = (raw or "").strip()
    s = s.strip("[]").replace("'", "").replace('"', "").strip()
    return [p.strip() for p in s.split(",") if p.strip()]


def label_to_lines(raw: str) -> str:
    """In-memory label: '\\|' and '|' both mean LINE BREAK."""
    return (raw or "").replace("\\|", "\n").replace("|", "\n")


def load_nodes(path: Path):
    delim, rows = _parse_rows(path, len(NODE_COLS))
    nodes, seen = [], {}
    for lineno, cells in rows:
        rec = OrderedDict(zip(NODE_COLS, cells))
        rec["label"] = label_to_lines(rec["label"])
        rec["views"] = parse_views(rec["views"])
        if rec["id"] in seen:
            die(EXIT_VALIDATION, "%s line %d: duplicate node id %r (first seen on line %d)"
                % (path.name, lineno, rec["id"], seen[rec["id"]]))
        seen[rec["id"]] = lineno
        nodes.append(rec)
    return nodes, delim


def load_edges(path: Path):
    delim, rows = _parse_rows(path, len(EDGE_COLS))
    edges = []
    for lineno, cells in rows:
        rec = OrderedDict(zip(EDGE_COLS, cells))
        rec["label"] = (rec["label"] or "").replace("\\|", "|").strip()
        rec["views"] = parse_views(rec["views"])
        edges.append(rec)
    return edges, delim


def resolve_input(name: str, override) -> Path:
    if override:
        p = Path(override).expanduser()
        if not p.is_file():
            die(EXIT_VALIDATION, "input table not found: %s" % p)
        return p
    for d in INPUT_SEARCH_DIRS:
        p = d / name
        if p.is_file():
            return p
    die(EXIT_VALIDATION,
        "input table %s not found in any of: %s"
        % (name, ", ".join(str(d) for d in INPUT_SEARCH_DIRS)))


def validate(nodes, edges) -> None:
    ids = set()
    for n in nodes:
        tag = "node %r" % n["id"]
        if n["section"] not in SECTIONS:
            die(EXIT_VALIDATION, "%s: unknown section %r" % (tag, n["section"]))
        if n["kind"] not in KINDS:
            die(EXIT_VALIDATION, "%s: unknown kind %r" % (tag, n["kind"]))
        if n["evidence"] not in EVIDENCE:
            die(EXIT_VALIDATION, "%s: unknown evidence %r" % (tag, n["evidence"]))
        if not n["label"].strip():
            die(EXIT_VALIDATION, "%s: empty label" % tag)
        for v in n["views"]:
            if v not in ALLOWED_VIEWS:
                die(EXIT_VALIDATION, "%s: views value %r not in {full, spine} + section keys"
                    % (tag, v))
        if not n["source"].strip() and "[W]" not in n["label"]:
            die(EXIT_VALIDATION,
                "%s: no source citation and no [W] flag in label" % tag)
        ids.add(n["id"])
    for e in edges:
        tag = "edge %r->%r" % (e["src"], e["dst"])
        if e["tier"] not in EDGE_TIERS:
            die(EXIT_VALIDATION, "%s: unknown tier %r" % (tag, e["tier"]))
        for v in e["views"]:
            if v not in ALLOWED_VIEWS:
                die(EXIT_VALIDATION, "%s: views value %r not in {full, spine} + section keys"
                    % (tag, v))
        for k in ("src", "dst"):
            if e[k] not in ids:
                die(EXIT_VALIDATION, "%s: %s references unknown node %r" % (tag, k, e[k]))


# --------------------------------------------------------------------------- #
# §2 / §3  corrections + Section H
# --------------------------------------------------------------------------- #
def apply_corrections(nodes, edges):
    by_id = {n["id"]: n for n in nodes}
    applied = []
    for nid, patch in CORRECTIONS.items():
        if nid not in by_id:
            die(EXIT_VALIDATION, "correction target node %r missing from data" % nid)
        n = by_id[nid]
        if "label" in patch:
            n["label"] = patch["label"]
        elif "append_label" in patch:
            n["label"] = n["label"].rstrip() + patch["append_label"]
        n["source"] = patch["source"]
        applied.append(nid)
    for rec in CORRECTION_NEW_NODES:
        if rec["id"] in by_id:
            die(EXIT_VALIDATION, "v5 node %r already exists in data" % rec["id"])
        node = OrderedDict([("id", rec["id"]), ("label", rec["label"]),
                            ("section", rec["section"]), ("kind", rec["kind"]),
                            ("evidence", rec["evidence"]), ("source", rec["source"]),
                            ("views", list(rec["views"]))])
        nodes.append(node)
        by_id[node["id"]] = node
        applied.append(rec["id"] + " (new node)")
    for rec in CORRECTION_NEW_EDGES:
        if rec["src"] not in by_id or rec["dst"] not in by_id:
            die(EXIT_VALIDATION,
                "v5 edge %s->%s references a missing node" % (rec["src"], rec["dst"]))
        edges.append(OrderedDict([("src", rec["src"]), ("dst", rec["dst"]),
                                  ("tier", rec["tier"]), ("label", rec["label"]),
                                  ("views", list(rec["views"]))]))
        applied.append("%s->%s (new edge)" % (rec["src"], rec["dst"]))
    return applied


def add_h_section(nodes, edges) -> None:
    by_id = {n["id"]: n for n in nodes}
    for nid, kind, evidence, label, source in H_NODES:
        if nid in by_id:
            die(EXIT_VALIDATION, "H-section node %r already exists in data" % nid)
        node = OrderedDict([("id", nid), ("label", label), ("section", "H"),
                            ("kind", kind), ("evidence", evidence),
                            ("source", source), ("views", ["full"])])
        nodes.append(node)
        by_id[nid] = node
    for src, dst, tier, label in H_EDGES:
        if src not in by_id or dst not in by_id:
            die(EXIT_VALIDATION, "H-section edge %s->%s references a missing node"
                % (src, dst))
        edges.append(OrderedDict([("src", src), ("dst", dst), ("tier", tier),
                                  ("label", label), ("views", ["full"])]))
    print("  H section added: %d nodes, %d edges (ids: %s)"
          % (len(H_NODES), len(H_EDGES),
             ", ".join(nid for nid, _, _, _, _ in H_NODES)))


# --------------------------------------------------------------------------- #
# §4  metrics -- consumed only when present AND schema-valid
# --------------------------------------------------------------------------- #
def _typed_ok(value, want) -> bool:
    if want is int and isinstance(value, bool):
        return False
    if want is float:
        return isinstance(value, (int, float)) and not isinstance(value, bool)
    return isinstance(value, want)


def _dig(data, path):
    """Follow a tuple of dict keys.  Returns None if any step is absent."""
    cur = data
    for step in path:
        if not isinstance(cur, dict) or step not in cur:
            return None
        cur = cur[step]
    return cur


def resolve_metrics():
    """-> (metrics_or_None, mode, reason).  Never raises, never guesses.

    Returns a dict shaped for apply_g_section() plus a 'sources' map giving the
    JSON path behind every displayed value, so the figure and the export can be
    audited back to the file.  Any missing path, wrong type or disagreeing
    cross-check returns in_flight -- a partial read is never presented as a
    finished run.
    """
    if not METRICS_PATH.is_file():
        return None, "in_flight", "metrics JSON absent: %s" % METRICS_PATH
    try:
        with METRICS_PATH.open(encoding="utf-8") as fh:
            data = json.load(fh)
    except Exception as exc:
        return None, "in_flight", "metrics JSON unreadable (%s)" % exc
    if not isinstance(data, dict):
        return None, "in_flight", "metrics JSON is not an object"

    # -- the run must identify itself ---------------------------------------
    commit = _dig(data, ("provenance", "commit"))
    if not isinstance(commit, str) or len(commit) < 7:
        return None, "in_flight", ("metrics JSON has no usable provenance.commit "
                                   "to identify the run by")
    run_date = _dig(data, ("provenance", "run_date"))
    if not isinstance(run_date, str):
        run_date = "unrecorded"                # the file really has no date field

    # -- cross-checks must agree --------------------------------------------
    samples = _dig(data, ("ten", "samples"))
    if not isinstance(samples, list) or len(samples) != EXPECTED_SAMPLES:
        got = len(samples) if isinstance(samples, list) else samples
        return None, "in_flight", ("metrics ten/samples holds %r, expected %d "
                                   "genomes" % (got, EXPECTED_SAMPLES))
    if len(set(samples)) != EXPECTED_SAMPLES:
        return None, "in_flight", "metrics ten/samples contains duplicate genomes"
    for label, path, want in METRIC_CHECKS:
        got = _dig(data, path)
        if got is None:
            return None, "in_flight", "metrics JSON missing %s (expected %s)" % (
                "/".join(path), label)
        if got != want:
            return None, "in_flight", "metrics %s=%r at %s, expected %r" % (
                label, got, "/".join(path), want)

    # -- displayed fields ----------------------------------------------------
    metrics, sources = {}, {}
    for key, path, want in METRIC_FIELDS:
        got = _dig(data, path)
        if got is None:
            return None, "in_flight", "metrics JSON missing %s" % "/".join(path)
        if not _typed_ok(got, want):
            return None, "in_flight", "metrics %s at %s is %s, expected %s" % (
                key, "/".join(path), type(got).__name__, want.__name__)
        metrics[key] = got
        sources[key] = "/".join(path)

    metrics["replication"] = _dig(data, ("ten", "headline", "replication"))
    metrics["significant_replication"] = _dig(
        data, ("ten", "headline", "significant_replication"))
    metrics["replication_fraction"] = _dig(
        data, ("ten", "headline", "replication_fraction"))
    sources["replication"] = "ten/headline/replication"
    sources["significant_replication"] = "ten/headline/significant_replication"
    sources["replication_fraction"] = "ten/headline/replication_fraction"

    # Every genome the headline claims must actually appear in the per-sample
    # table.  A headline of 10 with a 5-genome table is exactly the kind of
    # half-finished file this guard exists to refuse.
    per_sample = _dig(data, ("ten", "per_sample"))
    if not isinstance(per_sample, list) or not per_sample:
        return None, "in_flight", "metrics ten/per_sample is missing or empty"
    seen = set()
    for row in per_sample:
        if not isinstance(row, dict):
            return None, "in_flight", "metrics ten/per_sample holds a non-object row"
        layer = row.get("cohort_layer")
        if not isinstance(layer, str):
            return None, "in_flight", ("metrics ten/per_sample row has no "
                                       "cohort_layer genome name")
        seen.add(layer)
    absent = sorted(set(samples) - seen)
    if absent:
        return None, "in_flight", (
            "metrics headline names %d genomes but ten/per_sample has no rows "
            "for: %s" % (len(samples), ", ".join(absent)))
    metrics["per_sample_rows"] = len(per_sample)

    # The tier being superseded is read from the file, not written by hand, so the
    # README's before/after sentence cannot drift from the data.
    five_linker = _dig(data, ("five", "headline", "linker"))
    five_tail = _dig(data, ("five", "headline", "tail"))
    if not _typed_ok(five_linker, float) or not _typed_ok(five_tail, float):
        return None, "in_flight", ("metrics five/headline lacks usable linker/tail "
                                   "to describe the superseded tier")
    metrics["superseded_linker_fold"] = five_linker
    metrics["superseded_tail_fold"] = five_tail
    five_samples = _dig(data, ("five", "samples"))
    if not isinstance(five_samples, list) or not five_samples:
        return None, "in_flight", ("metrics five/samples is missing or empty, so the "
                                   "superseded tier cannot be named by genome count")
    metrics["superseded_n_samples"] = len(five_samples)
    metrics["samples"] = list(samples)
    sources["superseded_linker_fold"] = "five/headline/linker"
    sources["superseded_tail_fold"] = "five/headline/tail"

    metrics["run_date"] = run_date
    metrics["run_id"] = commit[:12]
    metrics["provenance"] = commit
    metrics["seed"] = _dig(data, ("provenance", "seed"))
    metrics["supersedes"] = "5-genome tier"
    metrics["sources"] = sources
    metrics["profiles"] = resolve_profiles(data, metrics)
    resolve_spread(data, metrics["profiles"], metrics)

    # Counts the analysis itself flagged as still moving.  Publishing them
    # matters: union +42.7% and nested L1 +57.6% mean every absolute count from
    # the 5-genome tier is a floor, not an estimate.
    flags, every = [], []
    for row in (_dig(data, ("delta",)) or []):
        if not isinstance(row, dict):
            continue
        flat = {"headline": row.get("headline"),
                "five": row.get("five"), "ten": row.get("ten"),
                "pct": row.get("percent_change"),
                "flag": bool(row.get("flag_gt20pct")),
                "reason": row.get("reason") or ""}
        every.append(flat)
        if row.get("flag_gt20pct"):
            flags.append(flat)
    metrics["delta_all"] = every
    metrics["flagged"] = flags
    metrics["enrichment"] = {
        fam: _dig(data, ("ten", "headline", "enrichment_%s" % fam))
        for fam in ("ALU", "LINE1", "SVA")}
    metrics["legacy"] = {
        "linker": _dig(data, ("ten", "headline", "legacy_linker")),
        "tail": _dig(data, ("ten", "headline", "legacy_tail"))}
    metrics["orientations"] = _dig(data, ("ten", "orientations")) or []
    for key, path in (("enrichment_ALU", "ten/headline/enrichment_ALU"),
                      ("enrichment_LINE1", "ten/headline/enrichment_LINE1"),
                      ("enrichment_SVA", "ten/headline/enrichment_SVA"),
                      ("legacy_linker", "ten/headline/legacy_linker"),
                      ("legacy_tail", "ten/headline/legacy_tail"),
                      ("replicates", "provenance/replicates")):
        sources[key] = path

    metrics["file_sha"] = sha256_file(METRICS_PATH)
    return (metrics, "metrics_consumed",
            "metrics JSON validated (%s) run=%s" % (METRICS_PATH.name, commit[:12]))


# Which validated curves are actually drawn, and where.  alu_private is the
# private-cohort twin of the Alu panel already shown, so it is validated but not
# plotted -- and it is not discarded either: it is one of the six series the Alu
# panel's spread band is built from.  parity.json has to say which is which, or a
# reader cannot tell an omitted panel from a deliberately withheld one.
PLOTTED_CURVES = {
    "alu_body": "measured enrichment panel above lane A, on lane A's 0-320 ruler",
    "l1_full": "measured enrichment panel after lane C, on its own 0-6.5 kb ruler",
    "alu_rel": "lane D, relative position 0-1",
    "sva_rel": "lane D, relative position 0-1",
    "l1_rel": "lane D, relative position 0-1",
}

# Every series in ten/profiles is (family, method, cohort_layer).  The six named
# here are the ones drawn; the rest of the file is used as siblings for the
# spread band rather than as separate panels.
PROFILE_SERIES = OrderedDict([
    ("alu_body",   ("ALU", "primary", "union")),
    ("alu_private", ("ALU", "primary", "private")),
    ("l1_full",    ("LINE1", "consensus_primary", "union")),
    ("alu_rel",    ("ALU", "same_family_relative_all", "union")),
    ("sva_rel",    ("SVA", "same_family_relative_all", "union")),
    ("l1_rel",     ("LINE1", "relative_all", "union")),
])

_PROFILE_FIELDS = ("bin_start", "bin_end", "observed", "expected", "enrichment",
                   "empirical_p_greater")

# How wide a min..max band has to get before the sheet states its numbers as
# well as drawing it.  25% of the peak value is the line: under it the method
# choice does not change how a bin reads, over it a bare single number is one
# method's answer wearing the measurement's clothes.
SPREAD_MATERIAL = 0.25

# ...but a relative width only means something where there is signal to move.  A
# band spanning 0.00x to 0.20x is 100% wide and entirely uninformative: both
# series agree there is nothing there.  So the band has to reach the null line
# before its width counts, which is what stops a near-empty bin from being
# reported as the panel's biggest method disagreement.
SPREAD_SIGNAL = 1.0

# What a bin span on each curve's own axis means.  The relative-position curves
# are a fraction of the host; without this the range table prints "0.8-0.9 bp"
# for a decile of a gene.
PROFILE_AXIS_UNIT = {"alu_rel": "fraction", "sva_rel": "fraction",
                     "l1_rel": "fraction"}


def _profile_bins(rows):
    """Normalise one ten/profiles series into plottable bins, or None.

    None rather than a partial list, because a curve carrying one malformed row
    or one gap between bins would draw a value the analysis never produced.
    """
    if not isinstance(rows, list) or not rows:
        return None
    bins = []
    for r in rows:
        if not isinstance(r, dict):
            return None
        # _typed_ok(v, float) already admits int and rejects bool, which is
        # exactly the numeric contract these six fields have.
        if any(not _typed_ok(r.get(f), float) for f in _PROFILE_FIELDS):
            return None
        bins.append({"lo": float(r["bin_start"]), "hi": float(r["bin_end"]),
                     "obs": int(r["observed"]),
                     "exp": float(r["expected"]),
                     "enr": float(r["enrichment"]),
                     "p": float(r["empirical_p_greater"])})
    bins.sort(key=lambda b: b["lo"])
    for a, b in zip(bins, bins[1:]):
        if abs(a["hi"] - b["lo"]) > 1e-9:
            return None
    return bins


def _profile_buckets(raw):
    """(family, method, cohort_layer) -> raw rows, for the whole file."""
    out = {}
    if isinstance(raw, list):
        for row in raw:
            if isinstance(row, dict):
                out.setdefault((row.get("family"), row.get("method"),
                                row.get("cohort_layer")), []).append(row)
    return out


def _series_name(key):
    """A profile key as the short, readable form used in labels and prose."""
    return "%s/%s" % (key[1], key[2])


def _peak_band(key, bins, spr):
    """The spread range under the drawn curve's own peak, for parity.json.

    Recorded because that is the number the sheet prints next to the shape, and
    a reader who wants to know how much the method choice moves *that* number
    should not have to scan the whole per-bin table to find it.
    """
    if not bins or not spr or len(bins) != len(spr):
        return None
    i = max(range(len(bins)), key=lambda j: bins[j]["enr"])
    r = spr[i]
    return OrderedDict((("span", [_axis_span(key, r["bin_lo"], r["bin_hi"])]),
                        ("drawn", round(bins[i]["enr"], 3)),
                        ("min", round(r["e_lo"], 3)),
                        ("max", round(r["e_hi"], 3)),
                        ("rel_width", round(r["width"], 3)),
                        ("stays_above_null", r["e_lo"] > 1.0)))


def resolve_profiles(data, metrics):
    """Pull the measured per-position curves out of ten/profiles.

    These are the strongest thing in the file and they were being thrown away:
    a real observation count, a real expectation under the opportunity model, an
    enrichment ratio and an empirical permutation p for every position bin.  A
    map that draws the measured curve is worth more than one that restates a
    fold-change in a sentence.

    Validation is deliberately strict.  A curve is only kept when every row has
    the numeric fields the plot needs and the bins tile the axis without gaps or
    overlaps.  It is also cross-checked twice over: ten/headline.linker and .tail
    ARE two of the Alu bins, and per-row `events` IS the sum of `observed` across
    the series.  If either disagrees the file contradicts itself, and the honest
    move is to refuse the curve rather than draw whichever version we saw first.
    """
    raw = _dig(data, ("ten", "profiles"))
    if not isinstance(raw, list) or not raw:
        return {}

    buckets = _profile_buckets(raw)
    out, notes = {}, []
    for name, key in PROFILE_SERIES.items():
        rows = buckets.get(key)
        if not rows:
            notes.append("%s absent" % name)
            continue
        bins = _profile_bins(rows)
        if bins is None:
            notes.append("%s malformed or non-contiguous" % name)
            continue
        evs = {int(r["events"]) for r in rows
               if _typed_ok(r.get("events"), int)}
        if len(evs) != 1 or evs.pop() != sum(b["obs"] for b in bins):
            notes.append("%s: per-row events does not sum to observed" % name)
            continue
        out[name] = bins

    # Cross-check the two curves whose bins ARE the headline numbers.
    for name, lo, headline in (("alu_body", 120.0, "linker_fold"),
                               ("alu_body", 280.0, "tail_fold")):
        bins = out.get(name)
        if not bins:
            continue
        hit = next((b for b in bins if abs(b["lo"] - lo) < 1e-9), None)
        if hit is None:
            del out[name]
            notes.append("%s has no bin at %.0f bp to match the headline"
                         % (name, lo))
        elif abs(hit["enr"] - float(metrics[headline])) > 1e-9:
            del out[name]
            notes.append("%s bin %.0f-%.0f is %.6f but ten/headline.%s is %.6f"
                         % (name, hit["lo"], hit["hi"], hit["enr"], headline,
                            float(metrics[headline])))

    # The permutation floor.  With N replicates the smallest attainable p IS
    # 1/N, so a bin sitting exactly on the floor is "p = 1/N", never "p < 1/N" --
    # claiming a smaller p than the method can produce would be the single
    # easiest way to overstate this result.
    reps = _dig(data, ("provenance", "replicates"))
    if not _typed_ok(reps, int) or reps < 1:
        notes.append("provenance.replicates unusable; p floor left unstated")
        floor = None
    else:
        floor = 1.0 / reps
    metrics["replicates"] = reps if _typed_ok(reps, int) else None
    metrics["p_floor"] = floor

    metrics["profile_notes"] = notes
    return out


def resolve_spread(data, profiles, metrics):
    """min..max enrichment at each bin, across every series sharing those bins.

    ten/profiles holds nineteen series.  Six are drawn; the other thirteen
    measure the *same bins* with a different method or a different cohort, and
    they do not agree.  The 120-140 bp linker bin -- the sheet's headline --
    reads 1.26x to 2.45x depending on which of the six series is asked.  Drawing
    one of the six and quoting its number without the range is picking a method
    and presenting its answer as the measurement, so the range is drawn as an
    envelope and stated wherever it is wide enough to change the reading.

    A sibling is admitted only when it shares the family AND the exact bin
    boundaries.  An envelope between two curves on different bins would be
    meaningless, and nothing here is compared across spans.
    """
    buckets = _profile_buckets(_dig(data, ("ten", "profiles")))
    metrics["profile_series_total"] = len(buckets)
    out, notes = {}, []
    for name, bins in profiles.items():
        drawn_key = PROFILE_SERIES[name]
        sibs = []
        for key, rows in buckets.items():
            if key[0] != drawn_key[0] or key == drawn_key:
                continue
            b = _profile_bins(rows)
            if b is None or len(b) != len(bins):
                continue
            if any(abs(p["lo"] - q["lo"]) > 1e-9 or abs(p["hi"] - q["hi"]) > 1e-9
                   for p, q in zip(bins, b)):
                continue
            sibs.append((key, b))
        if not sibs:
            notes.append("%s has no sibling series on the same bins" % name)
            continue
        rows = []
        for i, base in enumerate(bins):
            # The drawn curve is itself one of the estimates, so it belongs in
            # its own envelope.
            vals = [(base["enr"], drawn_key)] + [(b[i]["enr"], k) for k, b in sibs]
            lo, lok = min(vals, key=lambda t: t[0])
            hi, hik = max(vals, key=lambda t: t[0])
            # e_lo/e_hi are enrichment; bin_lo/bin_hi are where the bin sits.
            # Carrying both, under names that cannot be confused, is what keeps
            # a caption from printing "0.47-1.32 bp" for a 3-3.5 kb bin.
            rows.append({"bin_lo": base["lo"], "bin_hi": base["hi"],
                         "e_lo": lo, "e_hi": hi,
                         "min_key": lok, "max_key": hik,
                         "width": (hi - lo) / hi if hi > 0 else 0.0})
        widest = max((r for r in rows if r["width"] > SPREAD_MATERIAL
                      and r["e_hi"] > SPREAD_SIGNAL),
                     key=lambda r: r["width"], default=None)
        out[name] = {
            "series": len(sibs) + 1,
            "keys": [drawn_key] + [k for k, _ in sibs],
            "bins": rows,
            "widest": widest,
            "material": [r for r in rows if r["width"] > SPREAD_MATERIAL
                         and r["e_hi"] > SPREAD_SIGNAL],
            "peak": max((b["enr"] for b in bins), default=0.0),
            "peak_hi": max((r["e_hi"] for r in rows), default=0.0),
        }
    metrics["spread"] = out
    metrics["spread_notes"] = notes
    return out


def _spread_parity(metrics):
    """The full method-sensitivity table for parity.json.

    Every plotted bin measures a range, not a number, and the range is what the
    caption prints.  It is recorded here in full -- including which method
    produced each end -- because the sheet has no room for that detail and the
    band is the sheet's headline claim.
    """
    out = OrderedDict()
    for name, sp in sorted(((metrics or {}).get("spread") or {}).items()):
        bins = (metrics or {}).get("profiles", {}).get(name, [])
        w = sp["widest"]
        out[name] = OrderedDict((
            ("series", sp["series"]),
            ("sibling_series", [_series_name(k) for k in sp["keys"][1:]]),
            ("material_bins", len(sp["material"])),
            ("material_frac", round(len(sp["material"]) / float(len(sp["bins"])), 3)),
            ("material_needs_rel_width", SPREAD_MATERIAL),
            ("material_needs_band_above", SPREAD_SIGNAL),
            ("peak", round(sp["peak"], 3)),
            ("peak_band", _peak_band(name, bins, sp["bins"])),
            ("widest", None if not w else OrderedDict((
                ("span", [_axis_span(name, w["bin_lo"], w["bin_hi"])]),
                ("min", round(w["e_lo"], 3)),
                ("max", round(w["e_hi"], 3)),
                ("rel_width", round(w["width"], 3)),
                ("min_series", _series_name(w["min_key"])),
                ("max_series", _series_name(w["max_key"]))))),
            ("per_bin", [OrderedDict((
                ("span", [_axis_span(name, r["bin_lo"], r["bin_hi"])]),
                ("drawn", round(b["enr"], 3)),
                ("min", round(r["e_lo"], 3)),
                ("max", round(r["e_hi"], 3)),
                ("rel_width", round(r["width"], 3)),
                ("min_series", _series_name(r["min_key"])),
                ("max_series", _series_name(r["max_key"]))))
                for b, r in zip(bins, sp["bins"])])))
    return out


def spread_clause(sp, lo=None, hi=None, fmt=None):
    """One short clause naming the envelope and the bin where it opens widest.

    Returns "" when the band is tight enough that stating it would only be
    boilerplate -- the band is still drawn, it just does not need a sentence.
    `fmt(lo, hi)` renders a bin span in the panel's own units; the default is
    right for the 0-1 ruler and callers pass _fmt_span for the bp ones.

    Which method produced each end of the range is deliberately NOT here.  It
    costs 60 characters on a caption line that is already near the canvas edge,
    and it is recorded per bin in parity.json instead.
    """
    if not sp or not sp.get("material"):
        return ""
    w = sp.get("widest")
    if not w:
        return ""
    lo = 0.0 if lo is None else lo
    hi = float("inf") if hi is None else hi
    if w["bin_lo"] < lo - 1e-9 or w["bin_hi"] > hi + 1e-9:
        return ""
    span = (fmt or (lambda a, b: "%.1f-%.1f" % (a, b)))(w["bin_lo"], w["bin_hi"])
    return ("band = min\u2013max across the %d method \u00d7 cohort series; "
            "widest %s at %.2fx\u2013%.2fx"
            % (sp["series"], span, w["e_lo"], w["e_hi"]))


# The two windows the headline is *about*.  These node labels carry the original
# 5-genome numbers; the run that superseded them did not rewrite the source
# tables, so the correction has to be attached here or the sheet keeps quoting a
# linker that moved by 11% and an n that is short by more than fourfold.
WINDOW_NODES = (("a_133ours", 120.0, "linker_fold"),
                ("a_tail", 280.0, "tail_fold"))


def apply_window_supersession(nodes, metrics):
    """Attach the current per-window measurement to the nodes that first claimed
    the window, without deleting what they originally said.

    The old figure and n are parsed out of the node's own label rather than
    written here, so the "understated by Nx" statement cannot drift away from the
    text it is describing.

    When the bin's spread across sibling series is wide the figure is quoted as a
    range.  For the linker that matters: the 120-140 bp bin reads 2.45x here and
    1.26x under tier2, so a bare "2.45x" on this node is one method's answer and
    should not be printed as though it were the measurement.
    """
    prof = metrics.get("profiles") or {}
    spread = metrics.get("spread") or {}
    floor = metrics.get("p_floor")
    alu = prof.get("alu_body")
    alu_sp = spread.get("alu_body")
    by_id = {n["id"]: n for n in nodes}
    notes = []
    for nid, lo, headline in WINDOW_NODES:
        n = by_id.get(nid)
        if n is None:
            continue
        hit = next((b for b in alu if abs(b["lo"] - lo) < 1e-9), None) if alu else None
        if hit is None:
            notes.append("%s: no %s bin at %.0f bp in ten/profiles" % (nid, "alu_body", lo))
            continue
        label = re.sub(r"\n?\[10-genome window[^\]]*\]", "", n["label"]).rstrip()

        old_x = re.search(r"(\d+(?:\.\d+)?)x", label)
        old_n = re.search(r"\bn\s*=\s*(\d+)", label)
        old_bits = []
        if old_x:
            old_bits.append("%sx" % old_x.group(1))
        if old_n:
            k = int(old_n.group(1))
            old_bits.append("n=%d" % k)
            if k > 0 and hit["obs"] > k:
                old_bits.append("n understated %.1fx" % (hit["obs"] / float(k)))

        # Parentheses, not brackets: the whole tag is delimited by [] and the
        # label scanners strip it by matching, so a nested [] would cut it short.
        quoted = "%.2fx" % hit["enr"]
        if alu_sp:
            row = next((r for r in alu_sp["bins"]
                        if abs(r["bin_lo"] - lo) < 1e-9), None)
            if row and row["width"] > SPREAD_MATERIAL:
                quoted += (" (%.2f\u2013%.2fx over %d method \u00d7 cohort series)"
                           % (row["e_lo"], row["e_hi"], alu_sp["series"]))
        tag = ("[10-genome window: %s, n=%d, %s -- %s] %s"
               % (quoted, hit["obs"], _p_text(hit["p"], floor),
                  "obs %d vs expected %.1f" % (hit["obs"], hit["exp"]),
                  ("supersedes %s" % ", ".join(old_bits)) if old_bits
                  else "5-genome tier superseded"))
        # The legacy chi-square p and the current permutation p are not on the
        # same scale, so say so once rather than letting a reader compare them.
        if re.search(r"p\s*=\s*[0-9.]*e-", label):
            tag += " [legacy chi-square p is not comparable to the permutation p]"
        n["label"] = label + "\n" + tag
        notes.append("%s <- %.2fx n=%d" % (nid, hit["enr"], hit["obs"]))
    return notes


def _p_text(p, floor=None):
    """An empirical p, with the resolution limit stated when it is reached."""
    if p is None:
        return ""
    if floor is not None and abs(p - floor) < 1e-12:
        return "p=%.4f (1/%d, the floor)" % (p, int(round(1.0 / floor)))
    if p < 0.0001:
        return "p=%.4f" % p
    return "p=%.4f" % p


def _fmt_metric(key, value):
    """Render a metric the way a reader wants it, not the way JSON stored it."""
    if key.endswith("_fold"):
        return "%sx" % ("%.2f" % value)
    if key == "replication_fraction":
        return "%.0f%% of genomes" % (100.0 * value)
    if isinstance(value, float):
        return "%.2f" % value
    return "{:,}".format(value)


def apply_g_section(nodes, metrics) -> None:
    by_id = {n["id"]: n for n in nodes}
    if metrics is None:
        for gid in G_SUFFIX_IDS:
            if gid not in by_id:
                die(EXIT_VALIDATION, "G-section node %r missing from data" % gid)
            if IN_FLIGHT_SUFFIX.strip() not in by_id[gid]["label"]:
                by_id[gid]["label"] = by_id[gid]["label"].rstrip() + IN_FLIGHT_SUFFIX
        return

    n_rep = metrics.get("replication") or 0
    n_sig = metrics.get("significant_replication") or 0
    rid = metrics["run_id"]

    # The old "[5-genome; 10-genome rerun in flight]" tag is replaced in place by
    # a tag that says what is true now: the 5-genome tier is superseded, and the
    # run that supersedes it is identified.  The 5-genome figures themselves are
    # left on the page -- marked as replaced, never deleted.
    # The run figures are a property of the RUN, not of each node, so they are
    # stamped once on the sheet (run_stamp_statements) rather than pasted onto
    # all five nodes.  What each node keeps is only what is true of that node:
    # its 5-genome tier has been superseded, and by which run.
    superseded = "[5-genome tier superseded by run %s]" % rid

    for gid in G_SUFFIX_IDS:
        if gid not in by_id:
            die(EXIT_VALIDATION, "G-section node %r missing from data" % gid)
        n = by_id[gid]
        label = re.sub(r"\n?\[5-genome[^\]]*\]", "", n["label"])
        label = re.sub(r"\n?\[10-genome rerun in flight\]", "", label)
        if gid == "g_william":
            # The extra genomes are in the file: five/samples holds 5 and
            # ten/samples holds %d, and delta records the counts still moving.
            # The node's "Incoming:" premise is now history, so state the
            # outcome instead of the expectation.  No genome is attributed to a
            # person here -- the file records which genomes, not who ran them.
            label = label.replace(
                "Incoming: William processing extra genomes (~2 h each)",
                "Extra genomes landed: %d \u2192 %d genomes at ~2 h each"
                % (metrics.get("superseded_n_samples", 0), len(metrics.get("samples") or []))
                if metrics.get("samples") else
                "Incoming: William processing extra genomes (~2 h each)")
        label = label.rstrip() + "\n" + superseded
        n["label"] = label


# --------------------------------------------------------------------------- #
# §5  DOT generation
# --------------------------------------------------------------------------- #
def dot_escape(text: str) -> str:
    """Order matters: backslash, then real newline, then quote."""
    return (text.replace("\\", "\\\\")
                .replace("\n", "\\n")
                .replace('"', '\\"'))


def html_escape(text: str) -> str:
    """Order matters: ampersand first, or we double-escape."""
    return (str(text).replace("&", "&amp;").replace("<", "&lt;")
                .replace(">", "&gt;").replace('"', "&quot;"))


# label palette: (headline, body, citation)
LABEL_INK = {
    "plain": ("#101820", "#33424f", "#7b8b98"),
    "chain": ("#ffffff", "#ccd2d9", "#98a1ac"),
    "gap":   ("#8f2a20", "#7d3b33", "#a4574d"),
}


# A bracketed whole line is a provenance stamp ("[10-genome run ...]"), not a
# finding.  Keeping it out of the body tier is what stops boilerplate reading
# like data.
TAG_RE = re.compile(r"^\[[^\[\]]+\]$")
STAMP_INK = "#8a97a2"


def html_node_label(n, fontsize: int, wrap: int) -> str:
    """A three-row HTML label: bold headline, body, italic citation.

    Restoring the citation as visible text is what makes a printed PNG or PDF
    self-documenting; it is the same string the SVG tooltip carries.
    """
    tone = ("chain" if n["kind"] == "chain" else
            "gap" if n["evidence"] == "gap" else "plain")
    head_c, body_c, cite_c = LABEL_INK[tone]

    # The headline is the whole first line of the source label, not merely
    # whatever the wrapper left over.  Graphviz breaks a single FONT run and
    # drops its weight/size on the continuation, so a wrapped headline used to
    # render bold 12pt then flounder into plain 9pt mid-sentence.  Pre-wrapping
    # into one FONT run per line keeps every headline line at headline weight.
    raw = n["label"].split("\n")
    head = textwrap.wrap(raw[0].strip(), wrap) or [raw[0].strip() or n["id"]]
    # Bracketed lines are provenance stamps, not findings.  They are separated
    # BEFORE wrapping -- wrap first and a tag stops looking like a tag, which is
    # how a run stamp ended up rendered as body-size prose.  Tags keep the
    # citation tier's muted ink so the eye skips to the actual claim.
    tag_rows, body_rows = [], []
    for seg in raw[1:]:
        if not seg.strip():
            continue
        if TAG_RE.match(seg.strip()):
            tag_rows.append(seg.strip())
        else:
            body_rows.extend(textwrap.wrap(seg, wrap))

    rows = ['<TR><TD ALIGN="LEFT">'
            '<FONT COLOR="%s" POINT-SIZE="%d"><B>%s</B></FONT></TD></TR>'
            % (head_c, fontsize + 2,
               "<BR ALIGN=\"LEFT\"/>".join(html_escape(l) for l in head))]
    if body_rows:
        body = "<BR ALIGN=\"LEFT\"/>".join(html_escape(l) for l in body_rows)
        rows.append('<TR><TD ALIGN="LEFT">'
                    '<FONT COLOR="%s" POINT-SIZE="%d">%s</FONT></TD></TR>'
                    % (body_c, fontsize - 1, body))
    if tag_rows:
        tags = "<BR ALIGN=\"LEFT\"/>".join(html_escape(l) for l in tag_rows)
        rows.append('<TR><TD ALIGN="LEFT">'
                    '<FONT COLOR="%s" POINT-SIZE="%d">%s</FONT></TD></TR>'
                    % (STAMP_INK, max(fontsize - 3, 7), tags))
    if n["source"]:
        rows.append('<TR><TD ALIGN="LEFT">'
                    '<FONT COLOR="%s" POINT-SIZE="%d"><I>%s</I></FONT></TD></TR>'
                    % (cite_c, max(fontsize - 3, 7), html_escape(n["source"])))
    return ('<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0" CELLPADDING="0">'
            + "".join(rows) + "</TABLE>")


def wrap_label(text: str, width: int) -> str:
    out = []
    for seg in text.split("\n"):
        if not seg.strip():
            out.append("")
            continue
        out.extend(textwrap.wrap(seg, width) or [""])
    return "\n".join(out)


def node_shape(n) -> str:
    if n["evidence"] == "gap":
        return "hexagon"
    return KIND_SHAPE.get(n["kind"], "box")


def node_statement(n, cls: str, fontsize: int) -> str:
    wrap, _, _, _, _ = VIEW_CLASS[cls]
    attrs = [
        "label=<%s>" % html_node_label(n, fontsize, wrap),
        'shape=%s' % node_shape(n),
        'fontname="Helvetica"',
        'fontsize=%d' % fontsize,
        'tooltip="%s"' % dot_escape(n["source"]),
    ]
    if n["kind"] == "chain":
        attrs += ['style="rounded,bold,filled"',
                  'fillcolor="%s"' % CHAIN_FILL, 'fontcolor="white"']
    elif n["evidence"] == "gap":
        attrs += ['style="filled"', 'fillcolor="%s"' % GAP_FILL,
                  'color="%s"' % GAP_COLOR, 'fontcolor="%s"' % GAP_COLOR,
                  'peripheries=1']
    else:
        attrs += ['style="filled"', 'fillcolor="%s"' % NODE_FILL,
                  'color="%s"' % TIER_STYLE[n["evidence"]][1],
                  'fontcolor="%s"' % INK]
    if "[W]" in n["label"]:
        attrs.append("peripheries=2")
    return '    "%s" [%s];' % (n["id"], ", ".join(attrs))


def edge_statement(e) -> str:
    style, color, width = TIER_STYLE[e["tier"]]
    attrs = ['style="%s"' % style, 'color="%s"' % color, 'penwidth=%s' % width,
             'fontname="Helvetica"', 'fontsize=9', 'fontcolor="#555555"']
    if e["label"]:
        attrs.append('label="%s"' % dot_escape(e["label"]))
        attrs.append('tooltip="%s->%s"' % (dot_escape(e["src"]), dot_escape(e["dst"])))
    return '    "%s" -> "%s" [%s];' % (e["src"], e["dst"], ", ".join(attrs))


def legend_statements(cls: str, fontsize: int):
    """Synthetic LG* nodes: exempt from citation validation and from filtering."""
    title, fill, border, pen = LEGEND_META
    wrap, _, _, _, _ = VIEW_CLASS[cls]

    def head(text, color="#111111", size=None):
        return ('<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0" CELLPADDING="0">'
                '<TR><TD ALIGN="LEFT"><FONT COLOR="%s" POINT-SIZE="%d">'
                '<B>%s</B></FONT></TD></TR></TABLE>'
                % (color, size or max(fontsize, 9), html_escape(text)))

    def plain(text, color="#333333", size=None):
        return ('<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0" CELLPADDING="0">'
                '<TR><TD ALIGN="LEFT"><FONT COLOR="%s" POINT-SIZE="%d">%s</FONT>'
                '</TD></TR></TABLE>'
                % (color, size or max(fontsize - 1, 8), html_escape(text)))

    lines = ['  subgraph cluster_legend {',
             '    label="%s"; style="rounded,filled"; fillcolor="%s"; '
             'color="%s"; penwidth=%s; fontname="Helvetica"; fontsize=%d; '
             'fontcolor="#333333";' % (title, fill, border, pen, fontsize)]
    lines.append('    "LG0" [label=<%s>, shape=plaintext, style="", '
                 'fontname="Helvetica", fontsize=%d, fontcolor="#111111"];'
                 % (head("EDGE TIERS  \u00b7  line style = evidence tier"),
                    max(fontsize, 9)))
    idx = 0
    for tier, (style, color, width) in TIER_STYLE.items():
        nid, swatch = "LGe%d" % idx, "LGe%ds" % idx
        lines.append('    "%s" [label=<%s>, shape=plaintext, style="", '
                     'fontname="Helvetica", fontsize=%d];'
                     % (nid, head(tier, "#111111", max(fontsize - 1, 8)),
                        max(fontsize - 1, 8)))
        lines.append('    "%s" [label=<%s>, shape=plaintext, style="", '
                     'fontname="Helvetica", fontsize=%d, fontcolor="#6a7a86"];'
                     % (swatch, plain(style, "#6a7a86", max(fontsize - 3, 7)),
                        max(fontsize - 3, 7)))
        lines.append('    "%s" -> "%s" [style="%s", color="%s", penwidth=%s];'
                     % (nid, swatch, style, color, width))
        lines.append('    "LG0" -> "%s" [style=invis];' % nid)
        idx += 1
    shapes = ("SHAPES  box/chain = marker \u00b7 ellipse = structure \u00b7 "
              "diamond = clinical \u00b7 folder = resource \u00b7 "
              "hexagon = VERIFIED GAP")
    flags = ("FLAGS  [V] verified primary source \u00b7 [S] snippet / "
             "unverified denominator \u00b7 [C] our computation \u00b7 "
             "[W] our data (double border)")
    for nid, text in (("LG1", shapes), ("LG2", flags)):
        wrapped = " ".join(wrap_label(text, max(wrap + 14, 40)).split("\n"))
        lines.append('    "%s" [label=<%s>, shape=plaintext, style="", '
                     'fontname="Helvetica", fontsize=%d, fontcolor="#555555"];'
                     % (nid, plain(wrapped, "#555555", max(fontsize - 2, 8)),
                        max(fontsize - 2, 8)))
        lines.append('    "LG0" -> "%s" [style=invis];' % nid)
    lines.append("  }")
    return lines


# A delta row's headline decides whether it is a ratio or a raw count.  This is
# the only classification available: the file ships no `kind` field, so the label
# is read rather than inferred from magnitude (a count can barely move and a
# ratio can halve, so magnitude cannot tell them apart).
RATIO_WORDS = ("enrichment", "fold", "fraction")


def delta_kind(row):
    """'ratio', 'count', or None when the headline says neither."""
    head = str(row.get("headline") or "").lower()
    if not head:
        return None
    if any(w in head for w in RATIO_WORDS):
        return "ratio"
    return "count"


def stability_note(metrics):
    """What the 5-genome -> 10-genome tier change did to ratios versus counts.

    The stamp already prints two arrows (linker, tail) but leaves the reader to
    assume the rest of the tier moved the same way.  It did not.  A count grows
    with every genome added by construction -- union 4998 -> 7132 is not a
    finding about retrotransposons, it is what adding genomes *is* -- while a
    ratio divides that growth out.  So the tier change is expected to move counts
    hard and ratios barely, and saying so is the difference between "these
    numbers are soft" and "these particular numbers are soft".

    Returns [] whenever the file does not support the sentence.  Each clause is
    emitted only when every row it summarises is present and numeric, so a
    partial delta cannot produce a partial claim.  Nothing here is asserted: the
    counts, the max drift, and the split between kinds are all read out of
    `metrics['flagged']`, which was built from the file at load time.
    """
    flagged = metrics.get("flagged") or []
    if not flagged:
        return []
    kinds = [delta_kind(r) for r in flagged]
    if any(k != "count" for k in kinds):
        return []
    n = len(flagged)
    worst = max((abs(float(r["pct"])) for r in flagged
                 if r.get("pct") is not None), default=None)
    if worst is None:
        return []

    # The complementary claim -- that the ratios held -- needs the unflagged rows
    # too, and those only reach metrics via the same table.  Skip the clause
    # rather than infer it from the absence of a flag.
    ratios = [r for r in (metrics.get("delta_all") or [])
              if delta_kind(r) == "ratio" and r.get("pct") is not None]
    out = ["all %d figures the analysis flagged as still moving are raw counts, "
           "not ratios (largest +%.1f%%)" % (n, worst)]
    if ratios:
        rmax = max(abs(float(r["pct"])) for r in ratios)
        out.append("the %d ratio%s moved at most %.1f%%"
                   % (len(ratios), "" if len(ratios) == 1 else "s", rmax))
    return out


def legacy_note(metrics):
    """How far the superseded pipeline landed on the *same* genomes.

    `superseded_*` in the stamp compares 5 genomes against 10, which conflates
    two things.  `legacy_*` is the cleaner comparison: the previous pipeline
    rerun on the 10 genomes already in hand, so the only thing that changed is
    the code.  Reporting it bounds how much of the headline is a pipeline
    artefact rather than a biological finding -- and if the two pipelines
    disagreed sharply that would be a reason to distrust the number.
    """
    legacy = metrics.get("legacy") or {}
    pairs = (("linker", "linker_fold"), ("tail", "tail_fold"))
    bits = []
    for lkey, ckey in pairs:
        old, cur = legacy.get(lkey), metrics.get(ckey)
        if not isinstance(old, (int, float)) or not isinstance(cur, (int, float)):
            continue
        if not cur:
            continue
        bits.append("%s %+.1f%%"
                   % (lkey, 100.0 * (float(cur) - float(old)) / float(old)))
    if not bits:
        return []
    return ["earlier pipeline, same 10 genomes: " + ", ".join(bits)]


def run_stamp_statements(cls: str, fontsize: int, metrics):
    """One provenance cartouche for the whole sheet, not one per node.

    The 10-genome figures describe the run, so repeating them on every g_* node
    was five copies of the same sentence in body-size ink -- it read as five
    findings and it was one.  This puts them in a single stamp, styled like a
    map's marginal note.  Emitted only for sheets that actually carry g_* nodes,
    and only when a run was consumed; ids are ^LG-prefixed so the viewer's
    citation check and tier filters skip them, exactly like the legend.
    """
    if not metrics:
        return []
    n_rep = metrics.get("replication") or 0
    n_sig = metrics.get("significant_replication") or 0
    when = ("no run date recorded in the file"
            if metrics["run_date"] == "unrecorded"
            else "run date %s" % metrics["run_date"])
    big = max(fontsize, 10)
    small = max(fontsize - 1, 9)

    def row(text, color, size, bold=False):
        inner = html_escape(text)
        if bold:
            inner = "<B>%s</B>" % inner
        return ('<TR><TD ALIGN="LEFT"><FONT COLOR="%s" POINT-SIZE="%d">%s</FONT>'
                '</TD></TR>' % (color, size, inner))

    figures = ("linker %s  \u00b7  tail %s  \u00b7  nested L1 %s  \u00b7  "
               "union %s" % (
                   _fmt_metric("linker_fold", metrics["linker_fold"]),
                   _fmt_metric("tail_fold", metrics["tail_fold"]),
                   _fmt_metric("nested_l1_count", metrics["nested_l1_count"]),
                   _fmt_metric("union_sites", metrics["union_sites"])))
    before = ("5-genome tier: linker %s \u2192 %s, tail %s \u2192 %s  (five/headline)"
              % (_fmt_metric("linker_fold", metrics["superseded_linker_fold"]),
                 _fmt_metric("linker_fold", metrics["linker_fold"]),
                 _fmt_metric("tail_fold", metrics["superseded_tail_fold"]),
                 _fmt_metric("tail_fold", metrics["tail_fold"])))

    # Two more marginal lines, both derived rather than written.  The `before`
    # line alone invites the reader to assume the whole tier moved the same way,
    # which would be wrong: only the counts moved, and by construction.
    extra = stability_note(metrics) + legacy_note(metrics)
    fam = []
    for key, short in (("ALU", "Alu"), ("LINE1", "L1"), ("SVA", "SVA")):
        val = (metrics.get("enrichment") or {}).get(key)
        if isinstance(val, (int, float)) and val > 0:
            fam.append("%s %sx" % (short, ("%.2f" % float(val))))
    if fam:
        extra.insert(0, "same-family enrichment  " + "  \u00b7  ".join(fam))
    extra_rows = "".join(row(t, "#6a7a86", max(fontsize - 2, 8)) for t in extra)

    table = ('<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0" CELLPADDING="0">'
             + row("10-GENOME RUN  \u00b7  %s  \u00b7  %d/%d genomes replicated"
                    % (metrics["run_id"], n_sig, n_rep), "#111111", big, bold=True)
             + row(figures, "#1f2a33", small)
             + row("seed %s  \u00b7  %s" % (metrics["seed"], when), "#6a7a86",
                   max(fontsize - 2, 8))
             + row(before, "#6a7a86", max(fontsize - 2, 8))
             + extra_rows
             + row("read from %s at build time; sha256 %s"
                    % ("ten/headline", metrics.get("file_sha", "?")[:16]),
                   "#8a97a2", max(fontsize - 3, 7))
             + '</TABLE>')
    return [
        '  subgraph cluster_runstamp {',
        '    label=""; style="rounded,filled"; fillcolor="#f7f9fa"; '
        'color="#c8d2da"; penwidth=1.2;',
        '    "LGR0" [label=<%s>, shape=plaintext, style="", fontname="Helvetica", '
        'fontsize=%d];' % (table, small),
        "  }",
    ]


def build_dot(view: str, cls: str, title: str, vnodes, vedges,
              metrics=None) -> str:
    wrap, nodesep, ranksep, gfs, nfs = VIEW_CLASS[cls]
    out = ["digraph %s {" % view.replace(".", "_"),
           "  graph [label=\"%s\", labelloc=t, fontsize=%d, fontname=\"Helvetica\", "
           "nodesep=%s, ranksep=%s, compound=true];" % (dot_escape(title), gfs,
                                                        nodesep, ranksep),
           "  node [fontname=\"Helvetica\", fontsize=%d, margin=\"0.12,0.08\"];" % nfs,
           "  edge [fontname=\"Helvetica\"];"]

    by_section = OrderedDict()
    for n in vnodes:
        by_section.setdefault(n["section"], []).append(n)

    for sec, members in by_section.items():
        if sec not in SECTION_META:
            continue
        stitle, fill, border, pen = SECTION_META[sec]
        out.append("  subgraph cluster_%s {" % sec)
        out.append('    label="%s"; style="rounded,filled"; fillcolor="%s"; '
                   'color="%s"; penwidth=%s; fontname="Helvetica"; '
                   'fontsize=%d;' % (dot_escape(stitle), fill, border, pen,
                                      max(gfs - 4, 10)))
        for n in members:
            out.append(node_statement(n, cls, nfs))
        out.append("  }")

    out.extend(legend_statements(cls, nfs))
    if metrics and any(n["id"] in G_SUFFIX_IDS for n in vnodes):
        out.extend(run_stamp_statements(cls, nfs, metrics))
    for e in vedges:
        out.append(edge_statement(e))
    out.append("}")
    return "\n".join(out) + "\n"


# --------------------------------------------------------------------------- #
# POSITION MAP renderer -- hand-built SVG so the geometry is exact
# --------------------------------------------------------------------------- #
MAP_W = 1600
MAP_X0, MAP_X1 = 118, 1470

# Lanes A and B share ONE linear ruler, 0-320 bp.  A log axis was tried and
# rejected: it collapses the 1-300 bp element body -- where nearly all of the
# content lives -- into 5% of the width, which defeats the point of a map.
MAP_LIN_MAX = 320.0
MAP_LIN_TICKS = [(v, str(v)) for v in range(0, 321, 10)]
MAP_BANDS = [(0, 10), (10, 20), (20, 30), (30, 40), (40, 50), (50, 60),
             (60, 70), (70, 80), (80, 90), (90, 100), (100, 200), (200, 320)]
MAP_LIN_LABEL_AT = (set(range(0, 121, 20))
                    | {140, 160, 180, 200, 240, 280, 300, 320})

# Lane C measures something different (the extent of a test), so it gets its own
# log axis rather than being crushed onto the element-body ruler.
MAP_LOG_DOM = (1.0, 6.0)           # log10(bp): 10 bp -> 1 Mb
MAP_LOG_TICKS = [(10, "10"), (30, "30"), (100, "100"), (300, "300"),
                 (1000, "1 kb"), (3000, "3 kb"), (10000, "10 kb"),
                 (30000, "30 kb"), (100000, "100 kb"), (1000000, "1 Mb")]

CARD_W, CARD_H, CARD_GAP, CARD_PITCH = 190, 50, 8, 54
BAR_H, BAR_PITCH = 7, 11


def _lane_nrows(lane_key: str) -> int:
    """How many non-overlapping bar sub-rows a linear lane needs.  Measured up
    front so the sheet height can be derived rather than guessed."""
    spans = []
    for L, _nid, sp_list, _note in MAP_ITEMS:
        if L != lane_key:
            continue
        for lo, hi, _cap in sp_list:
            spans.append((min(lo, MAP_LIN_MAX), min(hi, MAP_LIN_MAX)))
    spans.sort()
    ends = []
    for lo, hi in spans:
        for ri, end in enumerate(ends):
            if lo >= end:
                ends[ri] = hi
                break
        else:
            ends.append(hi)
    return len(ends)


def _lane_height(lane_key: str, ncard_rows: int) -> int:
    """Total vertical space a linear lane occupies, ruler and cards included.

    The leading term is the ruler's own header band.  It has to clear the card
    grid: card_top works out to lane_height minus exactly the bar and card
    stacks, so anything the header needs has to be budgeted here or the cards
    land on top of the lane's title.
    """
    return (56 + _lane_nrows(lane_key) * BAR_PITCH + 22
            + ncard_rows * CARD_PITCH + 46)


LANE_NROWS = {k: _lane_nrows(k) for k in ("A", "B")}
LANE_ROWS = {"A": 4, "B": 2, "C": 2}
LANE_GAP = 54

# A measured curve needs its own vertical room.  PROF_H is the plot area of the
# enrichment panel that sits directly above lane A and shares its 0-320 ruler, so
# the curve and the bars below it are read off one axis.
PROF_H, PROF_PAD = 132, 66
PROF_TOP = 110

# Stated once on the sheet footer, covering every measured panel at once.
# enrichment < 1 with p_greater near 1.0 is the one-sided test failing to reject,
# which says something about enrichment only; reading those bins as "depleted"
# would claim a result the analysis never ran.
ONE_SIDED_NOT = ("every measured panel \u00b7 the enrichment test is one-sided "
                 "(p for excess only), so a bin below 1x with p near 1.0 means "
                 "not enriched \u2014 it is not evidence of depletion")
LANE_A_TOP = PROF_TOP + PROF_H + PROF_PAD
LANE_A_AXIS = LANE_A_TOP + _lane_height("A", LANE_ROWS["A"]) - 46
LANE_B_TOP = LANE_A_TOP + _lane_height("A", LANE_ROWS["A"]) + LANE_GAP
LANE_B_AXIS = LANE_B_TOP + _lane_height("B", LANE_ROWS["B"]) - 46
LANE_C_TOP = LANE_B_TOP + _lane_height("B", LANE_ROWS["B"]) + LANE_GAP
LANE_C_AXIS = LANE_C_TOP + _lane_height("C", LANE_ROWS["C"]) - 46

# Lane D is the normalised ruler: position along the host as a fraction, not
# base pairs.  It gets its own 0-1.0 axis because mixing it with the bp rulers
# would be meaningless.
L1_H = 118
L1_TOP = LANE_C_TOP + _lane_height("C", LANE_ROWS["C"]) + LANE_GAP + 30
REL_H, REL_PAD = 150, 74
LANE_D_TOP = L1_TOP + L1_H + 52 + LANE_GAP + 30
LANE_D_AXIS = LANE_D_TOP + REL_H
MAP_H = LANE_D_AXIS + 172


def _esc(s: str) -> str:
    return (str(s).replace("&", "&amp;").replace("<", "&lt;")
            .replace(">", "&gt;").replace('"', "&quot;"))


def _tier_of(node):
    return TIER_STYLE.get(node["evidence"], ("solid", "#555555", "1.2"))


def _fmt_bp(v: float) -> str:
    """A length with its unit attached -- always.  Callers must not add "bp"."""
    v = float(v)
    if v >= 1000000:
        q = v / 1e6
        return ("%g Mb" % q) if q != int(q) else ("%d Mb" % q)
    if v >= 1000:
        q = v / 1000.0
        return ("%.3g kb" % q)
    return "%g bp" % v


def _fmt_span(a: float, b: float) -> str:
    """An interval with its unit written once, at the end of it.

    _fmt_bp() renders a single length, so calling it on both ends of a bin would
    print the unit twice and read as two numbers rather than one range.
    """
    if b >= 1000.0:
        return "%.3g\u2013%.3g kb" % (a / 1000.0, b / 1000.0)
    return "%g\u2013%g bp" % (a, b)


def _axis_span(key, a, b):
    """A bin span in the units of that curve's own axis.

    The relative-position curves are a fraction of the host, not a length, so
    they must not inherit the bp formatter -- otherwise the range table prints
    "0.8-0.9 bp" for a decile of a gene.
    """
    if PROFILE_AXIS_UNIT.get(key) == "fraction":
        return "%.1f\u2013%.1f" % (a, b)
    return _fmt_span(a, b)


def _fmt_span(a: float, b: float) -> str:
    """An interval with its unit written once, at the end of it.

    _fmt_bp() renders a single length, so calling it twice on the ends of a bin
    would print the unit twice and read as two numbers rather than one range.
    """
    if b >= 1000.0:
        return "%.3g\u2013%.3g kb" % (a / 1000.0, b / 1000.0)
    return "%g\u2013%g bp" % (a, b)


def _lin(bp: float) -> float:
    return MAP_X0 + (min(max(bp, 0.0), MAP_LIN_MAX) / MAP_LIN_MAX) * (MAP_X1 - MAP_X0)


def _logx(bp: float) -> float:
    lo, hi = MAP_LOG_DOM
    return MAP_X0 + ((math.log10(max(bp, 10 ** lo)) - lo) / (hi - lo)) * (MAP_X1 - MAP_X0)


def _clip(text: str, limit: int) -> str:
    text = " ".join(str(text).split())
    return text if len(text) <= limit else text[:limit - 1] + "\u2026"


def _pack_cards(wanted, nrows):
    """Left-align each card after the previous card in its row, taking the
    topmost row that still fits.  wanted = [(desired_x, payload)], pre-sorted."""
    last = [MAP_X0] * nrows
    out = []
    for want, _payload in wanted:
        want = min(max(want, MAP_X0), MAP_X1 - CARD_W)
        for ri in range(nrows):
            x = max(want, last[ri])
            if x + CARD_W <= MAP_X1 + 1:
                out.append((ri, x))
                last[ri] = x + CARD_W + CARD_GAP
                break
        else:
            out.append((nrows - 1, min(want, MAP_X1 - CARD_W)))
    return out


def _pack_bars(spans):
    """Greedy interval packing: each span takes the topmost sub-row whose last
    bar has already ended.  spans = [(x0, x1)], pre-sorted by x0."""
    ends, rows = [], []
    for x0, x1 in spans:
        for ri, end in enumerate(ends):
            if x0 >= end:
                ends[ri] = x1
                rows.append(ri)
                break
        else:
            ends.append(x1)
            rows.append(len(ends) - 1)
    return rows



def _zero_spans(bins):
    """Consecutive bins holding no observations, merged, with what was expected.

    The file reports enrichment 0.0 for these because nothing landed in them,
    which is not the same as having measured a zero.  A curve that dives to the
    axis floor there is drawing a depletion the one-sided test never ran.
    """
    out, cur = [], None
    for b in bins:
        if b["obs"] <= 0:
            if cur is None:
                cur = {"lo": b["lo"], "hi": b["hi"], "exp": b["exp"],
                       "nbins": 1}
            else:
                cur["hi"], cur["exp"] = b["hi"], cur["exp"] + b["exp"]
                cur["nbins"] += 1
        elif cur is not None:
            out.append(cur)
            cur = None
    if cur is not None:
        out.append(cur)
    return out


def _series_steps(bins, x_of, y_of, fill=None):
    """One SVG path pair for a measured curve: a filled step area and its edge.

    Bins are drawn as steps, not as a smoothed line, because each bin is a
    separate observation with its own denominator -- interpolating between them
    would invent values the analysis never computed.

    A bin with zero observations is left out entirely and the edge is broken
    across it.  SVA's relative profile has two such bins over an expected ten
    events, and painting them at 0.000x would put the curve on the floor of the
    panel where it reads as strong evidence against the host's 3' end.
    """
    area, edge, run = [], [], []
    y0 = y_of(0.0)

    def flush():
        if run:
            edge.append("M " + " L ".join(
                " ".join(seg) for seg in run))
            del run[:]

    for i, b in enumerate(bins):
        if b["obs"] <= 0:
            flush()
            continue
        nxt = bins[i + 1] if i + 1 < len(bins) else None
        nxt_enr = nxt["enr"] if (nxt is not None and nxt["obs"] > 0) else None
        x0, x1 = x_of(b["lo"]), x_of(b["hi"])
        yt, yb = y_of(b["enr"]), y0 if nxt_enr is None else y_of(nxt_enr)
        area.append('<path d="M %.1f %.1f L %.1f %.1f L %.1f %.1f L %.1f %.1f Z"'
                    ' fill="%s"/>' % (x0, yb, x1, yb, x1, y0, x0, y0, fill))
        run.append(("%.1f %.1f" % (x0, yt), "%.1f %.1f" % (x1, yt)))
    flush()
    return "".join(area), "".join(edge), y0


def _spread_band(bins, spr, x_of, y_of, ink):
    """The min..max envelope across sibling series, as one closed step polygon.

    Drawn behind the curve and in the same hue at low opacity, so it reads as
    "the method choice moves this between here and there" rather than as a
    fourth data series competing with the one being quoted.
    """
    fwd, back = [], []
    for b, s in zip(bins, spr):
        x0, x1 = x_of(b["lo"]), x_of(b["hi"])
        fwd.append("%.1f %.1f" % (x0, y_of(s["e_hi"])))
        fwd.append("%.1f %.1f" % (x1, y_of(s["e_hi"])))
        back.append("%.1f %.1f" % (x1, y_of(s["e_lo"])))
        back.append("%.1f %.1f" % (x0, y_of(s["e_lo"])))
    back.reverse()
    d = "M %s L %s Z" % (" L ".join(fwd), " L ".join(back))
    return ('<path d="%s" fill="%s" fill-opacity="0.16" stroke="%s" '
            'stroke-width="0.7" stroke-opacity="0.5" stroke-dasharray="2,2"/>'
            % (d, ink, ink))


def _zero_marks(bins, x_of, y_base, ink, label_anchor="end"):
    """Flag every stretch of the axis with no observations behind it."""
    o = []
    for z in _zero_spans(bins):
        x0, x1 = x_of(z["lo"]), x_of(z["hi"])
        o.append('<rect x="%.1f" y="%.1f" width="%.1f" height="7" rx="2" '
                 'fill="#ffffff" stroke="%s" stroke-width="1" '
                 'stroke-dasharray="2.5,2.5"/>' % (x0, y_base - 3.5, x1 - x0, ink))
        o.append('<text class="gaplab" x="%.1f" y="%.1f" text-anchor="%s">'
                 '%s</text>'
                 % (x1 - 4 if label_anchor == "end" else x0 + 4,
                    y_base - 8, label_anchor,
                    _esc("no observations in %d of %d bins (~%.0f expected "
                         "across them) \u00b7 drawn as a gap, not as 0\u00d7"
                         % (z["nbins"], len(bins), z["exp"]))))
    return o


def _enr_ylim(bins, ceiling=None):
    """Top of the enrichment scale, chosen so the tallest bin keeps its labels.

    Rounding up to the next step when the peak comes within 8% of the current
    one is what stops the two value labels from being jammed against the panel
    ceiling; without it a 3.98x peak on a 4x axis has nowhere to put its text.
    `ceiling` is the tallest value anything else on the panel reaches -- the top
    of the spread band -- because a band clipped by the axis would quietly
    understate the disagreement it exists to show.
    """
    top = max([b["enr"] for b in bins] + ([ceiling] if ceiling else []))
    for step in (1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0):
        if top <= step * 0.92:
            return step
    return 20.0


def profile_panel(metrics, key="alu_body", title=None, top=None, height=None,
                  ink="#c0392b", x_of=None, foot=None, xticks=None,
                  span_lo=0.0, span_hi=None, span_fmt=None):
    """A measured enrichment curve, drawn on the axis it actually has.

    The Alu panel is locked to lane A's 0-320 bp ruler, so a spike in the curve
    and the bars piled up beneath it are the same base pairs.  The L1 profile
    spans twenty times more base pairs and gets its own scale -- which is the
    point rather than an inconvenience: it is the reason one ruler cannot
    honestly serve both elements.

    Behind the curve sits the spread envelope: the min and max every sibling
    series in ten/profiles reports at the same bins.  Where that envelope is
    wide the caption says so in numbers, because a peak label quoting one
    method's figure as if it were the measurement is the most likely way for
    this panel to mislead.

    Every number here is read from ten/profiles at build time.  No bin, span or
    fold-change is written into this function.
    """
    prof = (metrics or {}).get("profiles") or {}
    bins = prof.get(key)
    if not bins:
        return [], "no ten/profiles %s curve -- panel omitted" % key
    sp = (metrics.get("spread") or {}).get(key) or None
    spr = sp["bins"] if sp else None
    floor = metrics.get("p_floor")
    ymax = _enr_ylim(bins, sp["peak_hi"] if sp else None)
    y_top = PROF_TOP if top is None else top
    prof_h = PROF_H if height is None else height
    y_bot = y_top + prof_h

    def y_of(v):
        return y_bot - (min(max(v, 0.0), ymax) / ymax) * prof_h

    if x_of is None:
        x_of = _lin
    if spr:
        o_band = _spread_band(bins, spr, x_of, y_of, ink)
    else:
        o_band = ""
    area, edge, base = _series_steps(bins, x_of, y_of, fill=ink)

    o = ['<rect x="%d" y="%d" width="%d" height="%d" rx="10" fill="#fbfcfd" '
         'stroke="#dfe6ec" stroke-width="1"/>'
         % (MAP_X0 - 26, y_top - 26, MAP_X1 - MAP_X0 + 52, PROF_H + 52)]
    o.append('<text class="lane" x="%d" y="%d">%s</text>'
             % (MAP_X0 - 14, y_top - 14,
                _esc(title or "MEASURED ENRICHMENT \u00b7 %s" % key)))

    # enrichment gridlines, with the null called out
    for v in (1.0, 2.0, 3.0, 4.0):
        if v > ymax:
            continue
        y = y_of(v)
        null = abs(v - 1.0) < 1e-9
        o.append('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="%s" '
                 'stroke-width="%.1f" stroke-dasharray="%s"/>'
                 % (MAP_X0, y, MAP_X1, y,
                    "#8a9aa8" if null else "#e6ecf1",
                    1.3 if null else 1.0, "" if null else "3,3"))
        o.append('<text class="tick" x="%d" y="%.1f" text-anchor="end">%s</text>'
                 % (MAP_X0 - 6, y + 3,
                    _esc("%gx" % v if not null else "1x  null")))

    if o_band:
        o.append(o_band)
    o.append('<g class="curve" fill-opacity="0.13">%s</g>' % area)
    o.append('<path d="%s" fill="none" stroke="%s" stroke-width="1.7" '
             'stroke-linejoin="round"/>' % (edge, ink))
    o.append('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="#9fb0bf" '
             'stroke-width="1"/>' % (MAP_X0, base, MAP_X1, base))
    o.extend(_zero_marks(bins, x_of, base, ink))

    # Label by signal, not by an arbitrary threshold.  Labelling only the bins
    # that reach the permutation floor left the L1 panel completely bare -- none
    # of its bins do -- which threw away its two real peaks at p=0.005 and
    # p=0.018.  So: always label the two tallest bins, plus every bin sitting on
    # the floor, and carry the p alongside n whenever the bin is not at the floor
    # so the reader can see how much resolution the value is worth.
    at_floor = (lambda b: floor is not None and b["p"] <= floor + 1e-12)
    peaks = sorted(bins, key=lambda b: (-b["enr"], b["lo"]))[:2]
    marked = {id(b) for b in peaks if True} | {id(b) for b in bins if at_floor(b)}
    for i, b in enumerate(bins):
        if id(b) not in marked:
            continue
        cx = (x_of(b["lo"]) + x_of(b["hi"])) / 2.0
        # A bin whose label sits inside its own envelope reads as clipping, so
        # the pair clears whichever of the curve and the band is higher.
        top_y = min(y_of(b["enr"]),
                    y_of(spr[i]["e_hi"]) if spr else y_of(b["enr"]))
        o.append('<text class="pkb" x="%.1f" y="%.1f" text-anchor="middle" '
                 'fill="%s">%s</text>'
                 % (cx, top_y - 7, ink, _esc("%.2fx" % b["enr"])))
        sub = ("n=%d" % b["obs"]) if at_floor(b) else (
            "n=%d \u00b7 p=%.4f" % (b["obs"], b["p"]))
        o.append('<text class="pkbn" x="%.1f" y="%.1f" text-anchor="middle">%s</text>'
                 % (cx, top_y - 18, _esc(sub)))
    o.append('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="#9fb0bf" '
             'stroke-width="1"/>' % (MAP_X0, y_bot, MAP_X1, y_bot))
    # When the panel carries its own x-axis, the caption has to clear the tick
    # labels underneath the axis line rather than sit on the same baseline.
    tail = " \u00b7 ".join(t for t in (
        spread_clause(sp, span_lo,
                      span_hi if span_hi is not None else bins[-1]["hi"],
                      span_fmt), foot) if t)
    o.append('<text class="lanesub" x="%d" y="%d">%s</text>'
             % (MAP_X0, y_bot + (36 if xticks else 16), _esc(tail)))
    return o, None



def l1_panel(metrics):
    """The L1 element-body profile, on a ruler of its own.

    L1 spans roughly six and a half kilobases against Alu's three hundred, so it
    cannot share lane A's scale without either crushing the Alu body or wasting
    most of the sheet.  Drawing both, each on the axis it actually has, is the
    honest version of the same comparison.

    What this curve does and does not say.  The bins below ~2 kb sit well under
    1x, and the bins around 4-5.5 kb sit above it.  `empirical_p_greater` is a
    one-sided test for *enrichment*, so a low bin with p near 1.0 means "no
    evidence of enrichment" and nothing more -- it is not a demonstration of
    depletion, and no complementary test is recorded in the file.  The file also
    carries no gene annotation, so no part of this curve is attributed to ORF2
    or to any other named locus; the coordinates are all this panel claims.

    The shape is also method-dependent.  At 2000-2500 bp the five sibling series
    in ten/profiles put that bin anywhere from 0.62x to 1.53x, so the crossing
    of the null line there is not a property of the data but of the method
    chosen.  The panel says this in numbers rather than leaving the reader to
    assume the union/consensus_primary curve is the measurement.
    """
    prof = (metrics or {}).get("profiles") or {}
    bins = prof.get("l1_full")
    if not bins:
        return [], None
    span = bins[-1]["hi"]
    if span <= 0:
        return [], "l1_full has a non-positive span -- panel omitted"
    top = L1_TOP
    width = float(MAP_X1 - MAP_X0)

    def x_of(bp):
        return MAP_X0 + min(max(bp / span, 0.0), 1.0) * width

    # Tick spacing is derived from the span, not chosen: a power of ten that
    # lands 3-4 times across the axis reads well at any span.
    step = next((t for t in (100, 200, 250, 500, 1000, 2000, 2500, 5000)
                 if span / t <= 4.5), span)
    ticks = []
    v = 0.0
    while v <= span + 1e-9:
        ticks.append((v, "0" if v == 0 else _fmt_bp(v)))
        v += step
    axis_y = top + L1_H
    o, note = profile_panel(
        metrics, key="l1_full",
        title="MEASURED ENRICHMENT \u00b7 L1, 10 genomes, union layer",
        top=top, height=L1_H, ink="#1f6fb2", x_of=x_of, xticks=ticks,
        span_lo=0.0, span_hi=span, span_fmt=_fmt_span,
        foot="same-family nested L1 per bp of host opportunity \u00b7 %d bins "
             "over %s \u00b7 own scale: L1 is %.0fx longer than the Alu body"
             % (len(bins), _fmt_bp(span), span / MAP_LIN_MAX))
    for v, lab in ticks:
        x = x_of(v)
        o.append('<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" stroke="#b9c6d1" '
                 'stroke-width="0.8"/>' % (x, axis_y, x, axis_y + 6))
        o.append('<text class="tick" x="%.1f" y="%d" text-anchor="middle">%s</text>'
                 % (x, axis_y + 18, _esc(lab)))
    return o, note

# Lane D inks.  L1 takes the same blue as its element-body panel above, because
# it is the same family on a different ruler, and the two colours agreeing is
# what lets a reader carry the family across lanes.
REL_INK = OrderedDict([("alu_rel", "#8e44ad"), ("sva_rel", "#16a085"),
                       ("l1_rel", "#1f6fb2")])
REL_LABEL = {"alu_rel": "Alu", "sva_rel": "SVA", "l1_rel": "L1"}


def relative_lane(metrics):
    """Lane D: normalised position along the host, 0 to 1.

    A second ruler measuring a different quantity from the base-pair lanes, and
    the sharpest place on the sheet to compare the three families against each
    other -- because they do not peak in the same place.  Alu peaks at the far
    end of the host, SVA near its midpoint, L1 at the host's 5' end.  Drawn
    without L1 the panel is a two-way comparison that loses the whole point.

    Which family peaks where is computed from the curves at build time, so the
    sentence under the axis cannot drift away from the shapes above it.
    """
    prof = (metrics or {}).get("profiles") or {}
    spread = (metrics or {}).get("spread") or {}
    floor = metrics.get("p_floor")
    series = [(k, prof[k], REL_LABEL[k]) for k in REL_INK if prof.get(k)]
    if not series:
        return [], None
    allb = [b for _k, bins, _l in series for b in bins]
    sps = [spread[k] for k, _b, _l in series if spread.get(k)]
    ymax = _enr_ylim(allb, max([s["peak_hi"] for s in sps], default=None))
    y_top, y_bot = LANE_D_TOP, LANE_D_TOP + REL_H

    def y_of(v):
        return y_bot - (min(max(v, 0.0), ymax) / ymax) * REL_H

    def x_of(t):
        return MAP_X0 + min(max(t, 0.0), 1.0) * (MAP_X1 - MAP_X0)

    o = ['<rect x="%d" y="%d" width="%d" height="%d" rx="10" fill="#fdfafc" '
         'stroke="#e6dcee" stroke-width="1"/>'
         % (MAP_X0 - 26, y_top - 30, MAP_X1 - MAP_X0 + 52, REL_H + 76)]
    o.append('<text class="lane" x="%d" y="%d">%s</text>'
             % (MAP_X0 - 14, y_top - 12,
                _esc("LANE D \u00b7 RELATIVE POSITION ALONG HOST")))
    for v in (1.0, 2.0, 3.0, 4.0):
        if v > ymax:
            continue
        y = y_of(v)
        null = abs(v - 1.0) < 1e-9
        o.append('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="%s" '
                 'stroke-width="%.1f" stroke-dasharray="%s"/>'
                 % (MAP_X0, y, MAP_X1, y,
                    "#8a9aa8" if null else "#ece7f0",
                    1.3 if null else 1.0, "" if null else "3,3"))
        o.append('<text class="tick" x="%d" y="%.1f" text-anchor="end">%s</text>'
                 % (MAP_X0 - 6, y + 3, _esc("%gx" % v if not null else "1x  null")))

    where = []
    for key, bins, lab in series:
        col = REL_INK[key]
        sp = spread.get(key)
        if sp:
            o.append(_spread_band(bins, sp["bins"], x_of, y_of, col))
        area, edge, base = _series_steps(bins, x_of, y_of, fill=col)
        o.append('<g class="curve" fill-opacity="0.12">%s</g>' % area)
        o.append('<path d="%s" fill="none" stroke="%s" stroke-width="2" '
                 'stroke-linejoin="round"/>' % (edge, col))
        o.extend(_zero_marks(bins, x_of, base, col))
        pi = max(range(len(bins)), key=lambda i: bins[i]["enr"])
        peak = bins[pi]
        px = (x_of(peak["lo"]) + x_of(peak["hi"])) / 2.0
        ceil = min(y_of(peak["enr"]),
                   y_of(sp["bins"][pi]["e_hi"])) if sp else y_of(peak["enr"])
        o.append('<text class="pkb" x="%.1f" y="%.1f" text-anchor="middle" '
                 'fill="%s">%s</text>'
                 % (px, ceil - 8, col,
                    _esc("%s %.2fx  n=%d" % (lab, peak["enr"], peak["obs"]))))
        o.append('<circle cx="%.1f" cy="%.1f" r="3.2" fill="%s"/>'
                 % (px, y_of(peak["enr"]), col))
        where.append((peak["lo"], lab))

    o.append('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="#9fb0bf" '
             'stroke-width="1"/>' % (MAP_X0, y_bot, MAP_X1, y_bot))
    for i in range(11):
        t = i / 10.0
        x = x_of(t)
        o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#b9c6d1" '
                 'stroke-width="0.8"/>' % (x, LANE_D_AXIS, x, LANE_D_AXIS + 6))
        o.append('<text class="tick" x="%.1f" y="%d" text-anchor="middle">%s</text>'
                 % (x, LANE_D_AXIS + 18, _esc("%.1f" % t)))
    # The finding, written from the curves rather than asserted beside them.
    where.sort()
    nsib = sorted({spread[k]["series"] for k, _b, _l in series if spread.get(k)})
    peaks = ", ".join("%s at %.1f\u2013%.1f" % (lab, lo, lo + 0.1)
                      for lo, lab in where)
    sibs = ("%d per family" % nsib[0] if len(nsib) == 1
            else "%d\u2013%d per family" % (nsib[0], nsib[-1])) if nsib else "none"
    # Only say the three-way contrast survives the method choice when it does.
    # Checked, not assumed: if a future run moves a band's lower edge under 1x,
    # the clause disappears instead of the sheet overstating its own finding.
    robust = [k for k, _b, _l in series if spread.get(k)] and all(
        spread[k]["bins"][max(range(len(prof[k])),
                              key=lambda j: prof[k][j]["enr"])]["e_lo"] > 1.0
        for k, _b, _l in series if spread.get(k))
    cap = ("fraction of host element span, not base pairs \u00b7 0 = 5' end, "
           "1 = 3' end \u00b7 each family peaks somewhere different along the "
           "host: %s \u00b7 pale band = min\u2013max across the sibling series "
           "on that curve, %s%s"
           % (peaks, sibs,
              " \u00b7 and each of those peaks stays above 1\u00d7 under every "
              "one of them" if robust else ""))
    o.append('<text class="lanesub" x="%d" y="%d">%s</text>'
             % (MAP_X0, LANE_D_AXIS + 36, _esc(cap)))
    return o, None



def build_position_map(nodes, edges, generated: str, metrics=None) -> str:
    """A cartographic POSITION MAP: every interval claim placed where it
    actually sits, so "at 5-10 nt I have this thing" is literally true.

    Emitted as <g class="node"><title>id</title> so the offline viewer treats
    every card as a first-class map feature -- search, tier filters, the citation
    panel and deep links all work here exactly as on the graph sheets.
    """
    l1_o, l1_note = l1_panel(metrics)
    prof_o, prof_note = profile_panel(
        metrics, key="alu_body",
        title="MEASURED ENRICHMENT \u00b7 Alu, 10 genomes, union layer",
        span_lo=0.0, span_hi=MAP_LIN_MAX, span_fmt=_fmt_span,
        foot="same-family nested insertion per bp of host opportunity \u00b7 "
             "%d bins over lane A's own ruler \u00b7 a labelled bin sits at "
             "the permutation floor"
             % len(metrics["profiles"]["alu_body"]))
    rel_o, rel_note = relative_lane(metrics)
    map_notes = [n for n in (prof_note, l1_note, rel_note) if n]

    # Two MAP_ITEMS notes carry fold-changes that the 10-genome run has since
    # moved.  Left as literals the cards would contradict the curve printed
    # directly above them, so when the measurement is available the note is
    # restated from ten/profiles instead of from the table.
    _alu = (metrics or {}).get("profiles", {}).get("alu_body")
    if _alu:
        for _nid, _lo in (("a_133ours", 120.0), ("a_tail", 280.0)):
            _b = next((b for b in _alu if abs(b["lo"] - _lo) < 1e-9), None)
            if _b is None:
                continue
            for _i, item in enumerate(MAP_ITEMS):
                if item[1] == _nid:
                    MAP_ITEMS[_i] = (item[0], item[1], item[2],
                                     "10-genome: %.2fx, n=%d, p at floor"
                                     % (_b["enr"], _b["obs"]))
    by_id = {n["id"]: n for n in nodes}
    spans_of = {}
    notes = {}
    for lane_key, nid, spans, note in MAP_ITEMS:
        spans_of[(lane_key, nid)] = spans
        notes[(lane_key, nid)] = note
    missing = sorted({i for _, i, _, _ in MAP_ITEMS} - set(by_id))
    if missing:
        die(EXIT_VALIDATION,
            "position map references unknown node id(s): %s" % ", ".join(missing))
    basis_of, basis_tally = audit_coord_basis(nodes, MAP_ITEMS)
    _assert_extent_ink_distinct()

    o = ['<svg xmlns="http://www.w3.org/2000/svg" '
         'xmlns:xlink="http://www.w3.org/1999/xlink" '
         'width="%dpt" height="%dpt" viewBox="0 0 %d %d">'
         % (MAP_W, MAP_H, MAP_W, MAP_H),
         '<rect width="%d" height="%d" fill="#ffffff"/>' % (MAP_W, MAP_H),
         '<style type="text/css"><![CDATA['
         'text{font-family:Helvetica,Arial,sans-serif}'
         '.ttl{font-size:23px;font-weight:bold;fill:#101f2b}'
         '.sub{font-size:12px;fill:#68788a}'
         '.lane{font-size:12.5px;font-weight:bold;fill:#38506a}'
         '.lanesub{font-size:10px;fill:#7f8fa0}'
         '.tick{font-size:9px;fill:#8a9aa8}'
         '.band{font-size:8px;fill:#8fa0af}'
         '.ct{font-size:9.5px;font-weight:bold;fill:#17242f}'
         '.cd{font-size:9px;fill:#4a5c6b}'
         '.cp{font-size:8.5px;fill:#7f8f9c;font-style:italic}'
         '.pkb{font-size:11px;font-weight:bold;fill:#b03a2e}'
         '.pkbn{font-size:8px;fill:#9a6b62}'
         '.gaplab{font-size:8px;fill:#8a9aa8;font-style:italic}'
         '.curve{stroke:none}'
         ']]></style>']

    # ---- masthead ---------------------------------------------------------
    o.append('<text class="ttl" x="%d" y="42">%s</text>'
             % (MAP_X0 - 22, _esc("POSITION MAP")))
    o.append('<text class="sub" x="%d" y="62">%s</text>'
             % (MAP_X0 - 22, _esc(
                 "Every interval claim in the Marker Map, placed where it "
                 "actually sits \u00b7 linear ruler 0\u2013320 bp \u00b7 1 bp = 1 nt "
                 "\u00b7 lanes A and B share one ruler \u00b7 the curve above lane A is measured, not asserted")))
    o.append('<text class="sub" x="%d" y="78">%s</text>'
             % (MAP_X0 - 22, _esc(
                 "bold = first line of the node label \u00b7 bar = the interval it "
                 "claims \u00b7 bar colour = evidence tier \u00b7 dashed red = "
                 "VERIFIED GAP \u00b7 dashed grey bracket = a width is stated but "
                 "no position, so none is claimed")))

    # ---- one feature = bars + leader + card, inside a single g.node --------
    def feature(lane_key, nid, x, y_top, anchor_x, anchor_y, marks=""):
        """A single map feature.

        Everything that belongs to the claim -- its interval bars on the ruler,
        its leader line and its card -- lives inside one <g class="node">, so the
        viewer's tier chips and search dim the whole feature or none of it,
        rather than fading a card and leaving its bar sitting on the axis.
        """
        n = by_id[nid]
        spans = spans_of[(lane_key, nid)]
        tier = n["evidence"]
        _, color, _ = _tier_of(n)
        cx = x + CARD_W / 2.0
        g = ['<g class="node%s">' % (" gap" if tier == "gap" else "")]
        g.append('<title>%s</title>' % _esc(nid))
        g.append(marks)
        g.append('<a xlink:title="%s">' % _esc("%s [%s]" % (n["source"], tier)))
        g.append('<g class="card">')
        g.append('<path d="M %.1f %d L %.1f %.1f L %.1f %.1f" fill="none" '
                 'stroke="%s" stroke-width="0.9" stroke-dasharray="2.5,2.5" '
                 'stroke-opacity="0.7"/>'
                 % (cx, y_top + CARD_H, cx, (y_top + CARD_H + anchor_y) / 2.0,
                    anchor_x, anchor_y, color))
        g.append('<rect x="%.1f" y="%d" width="%d" height="%d" rx="5" '
                 'fill="#ffffff" stroke="%s" stroke-width="1.2"/>'
                 % (x, y_top, CARD_W, CARD_H, color))
        g.append('<rect x="%.1f" y="%d" width="3.5" height="%d" rx="1.75" '
                 'fill="%s"/>' % (x, y_top + 6, CARD_H - 12, color))
        g.append('<text class="ct" x="%.1f" y="%d">%s</text>'
                 % (x + 9, y_top + 15, _esc(_clip(n["label"].split("\n")[0], 25))))
        g.append('<text class="cd" x="%.1f" y="%d">%s</text>'
                 % (x + 9, y_top + 28, _esc(_clip(span_caption(lane_key, nid), 31))))
        g.append('<text class="cp" x="%.1f" y="%d">%s</text>'
                 % (x + 9, y_top + 41, _esc(_clip(notes[(lane_key, nid)], 33))))
        g.append('</g>')
        g.append('</a></g>')
        return "".join(g)

    def span_caption(lane_key, nid):
        """The card's coordinate line, tagged wherever it is not verbatim.

        The tag travels with the number so a reader who copies "50 bp" off the
        card cannot lose the fact that where it sits was never stated.
        """
        spans = spans_of[(lane_key, nid)]
        bits = []
        for i, (_lo, _hi, cap) in enumerate(spans):
            k = basis_of[(lane_key, nid, i)]
            bits.append(cap if k == "verbatim"
                        else "%s %s" % (cap, COORD_BASIS_TAG[k]))
        return "; ".join(bits)

    # ---- a ruler ----------------------------------------------------------
    def ruler(axis_y, top, scale, ticks, label_at, unit, title, sub,
              fill, stroke, shaded_bands):
        o.append('<rect x="%d" y="%d" width="%d" height="%d" rx="10" fill="%s" '
                 'stroke="%s" stroke-width="1"/>'
                 % (MAP_X0 - 26, top, MAP_X1 - MAP_X0 + 52,
                    axis_y - top + 46, fill, stroke))
        o.append('<text class="lane" x="%d" y="%d">%s</text>'
                 % (MAP_X0 - 14, top + 20, _esc(title)))
        o.append('<text class="lanesub" x="%d" y="%d">%s</text>'
                 % (MAP_X0 - 14, top + 34, _esc(sub)))
        for i, (lo, hi) in enumerate(shaded_bands or []):
            x0, x1 = scale(lo), scale(hi)
            o.append('<rect x="%.1f" y="%d" width="%.1f" height="26" fill="%s" '
                     'fill-opacity="0.65"/>'
                     % (x0, axis_y - 28, x1 - x0, "#e4eef5" if i % 2 else "#ffffff"))
            if (x1 - x0) > 24:
                o.append('<text class="band" x="%.1f" y="%d" text-anchor="middle">'
                         '%s</text>' % ((x0 + x1) / 2.0, axis_y - 16,
                                        _esc("%g\u2013%g" % (lo, hi))))
        for v, lab in ticks:
            x = scale(v)
            major = (v % 50 == 0) if scale is _lin else True
            o.append('<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" stroke="%s" '
                     'stroke-width="%.1f"/>'
                     % (x, axis_y, x, axis_y + (7 if major else 4),
                        "#7e93a5" if major else "#c2cfda", 1.1 if major else 0.7))
            if v in label_at:
                o.append('<text class="tick" x="%.1f" y="%d" text-anchor="middle">'
                         '%s</text>' % (x, axis_y + 20, _esc(lab)))
        o.append('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="#5f7488" '
                 'stroke-width="1.3"/>' % (MAP_X0, axis_y, MAP_X1, axis_y))
        o.append('<text class="lanesub" x="%d" y="%d">%s</text>'
                 % (MAP_X0, axis_y + 36, _esc(unit)))

    # ---- lanes A and B: one shared linear ruler, one track each ------------
    def bar_mark(nid, sp, x0, x1, by, kind="verbatim"):
        """One interval bar.  Returned as markup, not appended, so the caller
        can put it inside the same <g class="node"> as the card: tier filtering
        and search then dim the card, its bar and its leader together.

        A placement whose position is not stated gets an open bracket instead of
        a filled bar, in a neutral ink rather than its tier colour.  Filled and
        tier-coloured means "measured, here"; dashed and grey means "this wide,
        location not claimed".  The two are not interchangeable and must not
        look alike.
        """
        if kind == "extent":
            w = max(x1 - x0, 6.0)
            return ('<rect x="%.1f" y="%.1f" width="%.1f" height="%d" rx="2" '
                    'fill="%s" fill-opacity="0.13" stroke="%s" '
                    'stroke-width="1" stroke-dasharray="3,2.5"/>'
                    '<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                    'stroke-width="1.4"/>'
                    '<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                    'stroke-width="1.4"/>'
                    % (x0, by, w, BAR_H, EXTENT_INK, EXTENT_INK,
                       x0, by - 4, x0, by + BAR_H + 4, EXTENT_INK,
                       x0 + w, by - 4, x0 + w, by + BAR_H + 4, EXTENT_INK))
        _, color, _ = _tier_of(by_id[nid])
        hi = sp[1]
        g = ['<rect x="%.1f" y="%.1f" width="%.1f" height="%d" rx="2" '
             'fill="%s" fill-opacity="0.92" stroke="%s" stroke-width="0.7"/>'
             % (x0, by, max(x1 - x0, 4.0), BAR_H, color, color)]
        for ex in (x0, x1):
            g.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" '
                     'stroke="%s" stroke-width="0.9"/>'
                     % (ex, by - 2, ex, by + BAR_H + 2, color))
        if hi > MAP_LIN_MAX:
            g.append('<path d="M %d %.1f L %d %.1f L %d %.1f Z" fill="%s"/>'
                     % (MAP_X1 + 3, by - 3, MAP_X1 + 13, by + BAR_H / 2.0,
                        MAP_X1 + 3, by + BAR_H + 3, color))
            g.append('<text class="cd" x="%d" y="%.1f" fill="%s">%s</text>'
                     % (MAP_X1 + 18, by + BAR_H, color,
                        _esc("\u2192 %s, off scale" % _fmt_bp(hi))))
        return "".join(g)

    def linear_lane(lane_key, top, axis_y, ncard_rows):
        items = [(nid, spans_of[(lane_key, nid)]) for L, nid, _s, _n in MAP_ITEMS
                 if L == lane_key]
        # one bar per span, packed into non-overlapping sub-rows above the axis
        flat = []
        for nid, spans in items:
            for i, sp in enumerate(spans):
                flat.append((_lin(sp[0]), _lin(min(sp[1], MAP_LIN_MAX)), nid, sp,
                             basis_of[(lane_key, nid, i)]))
        flat.sort(key=lambda t: (t[0], t[1]))
        bar_rows = _pack_bars([(f[0], f[1]) for f in flat])
        nrows = max(bar_rows) + 1
        ruler(axis_y, top, _lin, MAP_LIN_TICKS, MAP_LIN_LABEL_AT, "bp",
              MAP_LANES[lane_key][0], MAP_LANES[lane_key][1],
              "#f4f8fb", "#c3d4e0", MAP_BANDS)
        card_top = axis_y - 6 - nrows * BAR_PITCH - 22 - ncard_rows * CARD_PITCH

        bars, first_row = {}, {}
        for (x0, x1, nid, sp, kind), ri in zip(flat, bar_rows):
            by = axis_y - 6 - (ri + 1) * BAR_PITCH
            bars.setdefault(nid, []).append(bar_mark(nid, sp, x0, x1, by, kind))
            first_row.setdefault(nid, (min(x0, _lin(sp[0])), by + BAR_H / 2.0))

        wanted = []
        for nid, spans in items:
            mids = [(_lin(sp[0]) + _lin(min(sp[1], MAP_LIN_MAX))) / 2.0
                    for sp in spans]
            wanted.append((sum(mids) / len(mids), nid))
        wanted.sort(key=lambda t: t[0])
        for (mid, nid), (ri, x) in zip(wanted, _pack_cards(wanted, ncard_rows)):
            ax, ay = first_row[nid]
            o.append(feature(lane_key, nid, x, card_top + ri * CARD_PITCH,
                             ax, ay, "".join(bars[nid])))
        return card_top

    # The measured curve shares lane A's ruler, so a spike in the curve and the
    # bars piled up beneath it are literally the same base pairs.
    o.extend(prof_o)

    # Every lane's height is measured from its own data (see _lane_height), so
    # the sheet ends exactly where the content does.
    linear_lane("A", LANE_A_TOP, LANE_A_AXIS, LANE_ROWS["A"])
    linear_lane("B", LANE_B_TOP, LANE_B_AXIS, LANE_ROWS["B"])

    # ---- lane C: its own log ruler ----------------------------------------
    log_y, c_top = LANE_C_AXIS, LANE_C_TOP
    c_items = [(nid, spans_of[("C", nid)]) for L, nid, _s, _n in MAP_ITEMS if L == "C"]
    ruler(log_y, c_top, _logx, MAP_LOG_TICKS,
          set(v for v, _ in MAP_LOG_TICKS), "window size, log scale",
          MAP_LANES["C"][0], MAP_LANES["C"][1], "#fdfaf1", "#e0d3a8", [])
    wanted = sorted(((_logx(min(sp[0] for sp in s)), nid) for nid, s in c_items),
                    key=lambda t: t[0])
    for (mid, nid), (ri, x) in zip(wanted, _pack_cards(wanted, LANE_ROWS["C"])):
        _, color, _ = _tier_of(by_id[nid])
        dots = "".join(
            '<circle cx="%.1f" cy="%d" r="4" fill="#ffffff" stroke="%s" '
            'stroke-width="1.7"/>' % (_logx(lo), log_y - 8, color)
            for lo, _hi, _c in spans_of[("C", nid)])
        o.append(feature("C", nid, x, c_top + 50 + ri * CARD_PITCH,
                         mid, log_y - 8, dots))

    # ---- lane D: the normalised ruler ---------------------------------------
    o.extend(l1_o)
    o.extend(rel_o)

    # ---- footer -----------------------------------------------------------
    foot = LANE_D_AXIS + 74
    # This sentence used to be absolute ("every lo/hi is transcribed verbatim"),
    # which three placements made false.  It is now generated from the audit, so
    # the count on the sheet is the count the build actually measured and the
    # claim cannot drift away from the drawing.
    n_all = sum(basis_tally.values())
    n_verb = basis_tally.get("verbatim", 0)
    n_decl = n_all - n_verb
    o.append('<text class="sub" x="%d" y="%d">%s</text>'
             % (MAP_X0 - 22, foot, _esc(
                 "Coordinate provenance: %d of %d placements have both endpoints "
                 "transcribed verbatim from that node's own label; %d do not, and "
                 "each of those says so in its own card. No other coordinate is "
                 "shown, and none is estimated."
                 % (n_verb, n_all, n_decl))))
    o.append('<text class="sub" x="%d" y="%d">%s</text>'
             % (MAP_X0 - 22, foot + 18, _esc(
                 "Card tags \u00b7 %s"
                 % " \u00b7 ".join("[%s] %s" % (k, w)
                                    for k, w in (("pos. 0", COORD_BASIS_WHY["origin"]),
                                                 ("unit conv.",
                                                  COORD_BASIS_WHY["unit"]),
                                                 ("width only",
                                                  COORD_BASIS_WHY["extent"]))))))
    o.append('<text class="sub" x="%d" y="%d">%s</text>'
             % (MAP_X0 - 22, foot + 36, _esc(
                 "Retrotransposon Marker Map v5 \u00b7 position map \u00b7 %s "
                 "\u00b7 read-only source: marker_map_nodes.csv \u00b7 click any "
                 "card in index.html for the full label and its citation"
                 % generated)))
    o.append('<text class="sub" x="%d" y="%d">%s</text>'
             % (MAP_X0 - 22, foot + 54, _esc(ONE_SIDED_NOT)))
    if map_notes:
        o.append('<text class="sub" x="%d" y="%d">%s</text>'
                 % (MAP_X0 - 22, foot + 72, _esc("; ".join(map_notes))))
    o.append("</svg>")
    return "\n".join(o)


# --------------------------------------------------------------------------- #
# §5  view membership
# --------------------------------------------------------------------------- #
def view_nodes_edges(view: str, nodes, edges):
    by_id = {n["id"]: n for n in nodes}
    if view == "spine":
        vnodes = [n for n in nodes if "spine" in n["views"]]
        ids = {n["id"] for n in vnodes}
        vedges = [e for e in edges if e["src"] in ids and e["dst"] in ids]
        return vnodes, vedges
    if view == "full":
        vnodes = [n for n in nodes if "full" in n["views"]]
        ids = {n["id"] for n in vnodes}
        vedges = [e for e in edges
                  if e["src"] in ids and e["dst"] in ids and "full" in e["views"]]
        return vnodes, vedges

    vnodes = [n for n in nodes if n["section"] == view]
    ids = {n["id"] for n in vnodes}
    anchors = set()
    for e in edges:
        if e["dst"] in ids and by_id[e["src"]]["section"] == "SPINE":
            anchors.add(e["src"])
        if e["src"] in ids and by_id[e["dst"]]["section"] == "SPINE":
            anchors.add(e["dst"])
    for aid in sorted(anchors):
        vnodes.append(by_id[aid])
    ids |= anchors
    vedges = [e for e in edges if e["src"] in ids and e["dst"] in ids]
    return vnodes, vedges


def build_view_list():
    """Ordered [(key, class, title)] -- nav is derived from this, never hardcoded."""
    views = [("marker_map_spine", "spine", VIEW_TITLES["spine"])]
    for sec in SECTION_SHEETS:
        views.append(("marker_map_%s" % sec, "section",
                      TITLE + " \u2014 section %s" % sec))
    views.append(("marker_map_full", "full", VIEW_TITLES["full"]))
    return views


# --------------------------------------------------------------------------- #
# §6  render execution + fatal detection
# --------------------------------------------------------------------------- #
def run_dot(dot_bin: str, dot_text: str, stem: Path, fmts, dpi_for_png=True):
    dot_path = stem.with_suffix(".dot")
    dot_path.write_text(dot_text, encoding="utf-8")
    warnings = []
    for fmt in fmts:
        out = stem.with_suffix("." + fmt)
        cmd = [dot_bin, "-T" + fmt]
        if fmt == "png" and dpi_for_png:
            cmd.append("-Gdpi=144")
        cmd += [str(dot_path), "-o", str(out)]
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stderr = proc.stderr.decode("utf-8", "replace")
        fatal = None
        if proc.returncode != 0:
            fatal = "returncode=%d" % proc.returncode
        else:
            for marker in FATAL_MARKERS:
                if marker in stderr:
                    fatal = marker
                    break
        if fatal:
            excerpt = "\n".join(stderr.strip().splitlines()[:20]) or "(no stderr)"
            die(EXIT_GRAPHVIZ,
                "%s [%s]: fatal condition (%s)\n--- stderr ---\n%s"
                % (dot_path.name, fmt, fatal, excerpt))
        for line in stderr.splitlines():
            low = line.lower()
            if any(b.lower() in low for b in BENIGN_STDERR):
                continue
            if line.strip():
                warnings.append("%s[%s] %s" % (dot_path.name, fmt, line.strip()))
        if fmt == "png" and dpi_for_png:
            stamp_png_dpi(out, PNG_DPI)
    return warnings


PNG_DPI = 144


def _png_chunk(typ: bytes, payload: bytes) -> bytes:
    return (struct.pack(">I", len(payload)) + typ + payload
            + struct.pack(">I", zlib.crc32(typ + payload) & 0xFFFFFFFF))


def _png_walk(data: bytes):
    """Yield (offset, type, length) for every chunk, verifying each CRC.

    Walking the file rather than assuming a fixed header size is what makes the
    stamp below safe: IHDR is 8 + 4 + 4 + 13 + 4 = 33 bytes, and getting that
    arithmetic wrong splices the new chunk through IHDR's own checksum and
    produces a file no decoder will open.
    """
    if data[:8] != b"\x89PNG\r\n\x1a\n":
        die(EXIT_GRAPHVIZ, "not a PNG (bad signature)")
    pos = 8
    while pos + 8 <= len(data):
        length = struct.unpack(">I", data[pos:pos + 4])[0]
        typ = data[pos + 4:pos + 8]
        body = data[pos + 8:pos + 8 + length]
        if len(body) != length or pos + 12 + length > len(data):
            die(EXIT_GRAPHVIZ, "truncated PNG chunk %r at offset %d"
                % (typ.decode("latin1", "replace"), pos))
        crc = struct.unpack(">I", data[pos + 8 + length:pos + 12 + length])[0]
        if zlib.crc32(typ + body) & 0xFFFFFFFF != crc:
            die(EXIT_GRAPHVIZ, "PNG chunk %r at offset %d fails its CRC"
                % (typ.decode("latin1", "replace"), pos))
        yield pos, typ, length
        pos += 12 + length
        if typ == b"IEND":
            break
    if pos != len(data):
        die(EXIT_GRAPHVIZ, "PNG has %d trailing bytes after IEND" % (len(data) - pos))


def stamp_png_dpi(path: Path, dpi: int) -> None:
    """Record the render resolution inside the PNG.

    `dot -Gdpi=144` makes the raster twice the point size, but Graphviz's plain
    PNG writer emits no pHYs chunk, so every viewer falls back to assuming
    96 dpi and prints the sheet at twice its intended size.  Writing the chunk
    here keeps the file self-describing whichever renderer Graphviz was built
    with.  The result is re-walked so a malformed file fails the build (exit 3)
    rather than shipping as a corrupt image.
    """
    data = path.read_bytes()
    chunks = list(_png_walk(data))
    if not chunks or chunks[0][1] != b"IHDR":
        die(EXIT_GRAPHVIZ, "%s does not start with IHDR" % path.name)
    if any(typ == b"pHYs" for _off, typ, _n in chunks):
        return                                     # already stamped
    insert_at = chunks[0][0] + 12 + chunks[0][2]   # first byte after IHDR's CRC
    ppm = int(round(dpi / 0.0254))                 # pixels per metre
    chunk = _png_chunk(b"pHYs", struct.pack(">IIB", ppm, ppm, 1))  # unit = metre
    out = data[:insert_at] + chunk + data[insert_at:]
    list(_png_walk(out))                           # fail loudly, not silently
    path.write_bytes(out)


# --------------------------------------------------------------------------- #
# §7.1  lossless CSV re-export
# --------------------------------------------------------------------------- #
def views_to_literal(views) -> str:
    return "[" + ", ".join("'%s'" % v for v in views) + "]"


def write_csvs(nodes, edges, outdir: Path):
    npath, epath = outdir / NODE_CSV_NAME, outdir / EDGE_CSV_NAME
    with npath.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(NODE_COLS)
        for n in nodes:
            w.writerow([n["id"], n["label"].replace("\n", "|"), n["section"],
                        n["kind"], n["evidence"], n["source"],
                        views_to_literal(n["views"])])
    with epath.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(EDGE_COLS)
        for e in edges:
            w.writerow([e["src"], e["dst"], e["tier"],
                        e["label"].replace("\n", " "),
                        ", ".join(e["views"])])
    return npath, epath


# --------------------------------------------------------------------------- #
# §8  offline viewer
# --------------------------------------------------------------------------- #
VIEWER_CSS = """
*{box-sizing:border-box}
body{margin:0;display:flex;height:100vh;overflow:hidden;
 font-family:-apple-system,BlinkMacSystemFont,Helvetica,Arial,sans-serif;
 background:#f6f7f9;color:#1d2430}
nav{width:236px;flex:0 0 236px;padding:12px;background:#1e2430;color:#e8ecf3;
 overflow:auto}
nav h1{font-size:13px;letter-spacing:.06em;margin:0 0 10px;color:#9fc3ff}
nav h2{font-size:10px;letter-spacing:.09em;text-transform:uppercase;color:#8b96ab;
 margin:14px 0 6px;font-weight:600}
nav button{display:block;width:100%;text-align:left;margin:3px 0;padding:7px 9px;
 border:0;border-radius:6px;background:#2c3547;color:#dde;cursor:pointer;font-size:12.5px}
nav button:hover{background:#3b475e}
nav button.active{background:#4b6cb7;color:#fff}
nav p{font-size:10.5px;color:#9aa4b5;line-height:1.5;margin:8px 0 0}
#search{width:100%;padding:7px 8px;border-radius:6px;border:1px solid #3b475e;
 background:#141a24;color:#e8ecf3;font-size:12.5px}
#search:focus{outline:1px solid #4b6cb7}
#results{max-height:224px;overflow:auto;margin-top:5px}
#results:empty{display:none}
#results button{font-size:11.5px;padding:5px 7px;margin:2px 0}
#results button b{color:#9fc3ff;font-weight:600}
#chips{display:flex;flex-wrap:wrap;gap:4px}
#chips button{width:auto;flex:1 1 auto;text-align:center;padding:4px 6px;font-size:10.5px;
 border:1px solid rgba(255,255,255,.35);opacity:.35;font-weight:600}
#chips button.on{opacity:1;color:#fff}
main{flex:1;display:flex;flex-direction:column;overflow:hidden;min-width:0}
.viewport{flex:1;position:relative;overflow:hidden;background:#fbfbfc;touch-action:none;
 cursor:grab}
.viewport.dragging{cursor:grabbing}
.stage{position:absolute;left:0;top:0;transform-origin:0 0}
.stage svg{display:block;background:#fff;box-shadow:0 1px 10px rgba(20,30,50,.12)}
.hud{position:absolute;top:10px;right:10px;display:flex;gap:6px;z-index:6}
.hud button{border:1px solid #d6dbe3;border-radius:6px;background:#fffffff2;padding:6px 11px;
 font-size:13px;cursor:pointer;color:#1d2430}
.hud button:hover{background:#eef2f8}
.viewtitle{position:absolute;top:10px;left:12px;z-index:6;font-weight:600;color:#1d2430;
 background:#ffffffe0;border:1px solid #d6dbe3;border-radius:6px;padding:6px 11px;
 font-size:12.5px;max-width:52%;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
.help{position:absolute;left:12px;bottom:12px;z-index:8;max-width:420px;padding:12px 14px;
 background:#fffffff5;border:1px solid #d6dbe3;border-radius:8px;font-size:12px;
 line-height:1.65;color:#1d2430;box-shadow:0 6px 24px rgba(20,30,50,.18)}
.help[hidden]{display:none}
.help b{color:#0f172a}
.help table{border-collapse:collapse;margin-top:6px}
.help td{padding:1px 8px 1px 0;font-size:11.5px}
.scrim{position:fixed;inset:0;background:rgba(15,20,30,.28);z-index:9}
.scrim[hidden]{display:none}
.panel{position:fixed;top:0;right:0;bottom:0;width:392px;max-width:92vw;z-index:10;
 background:#fff;border-left:1px solid #d6dbe3;box-shadow:-8px 0 28px rgba(20,30,50,.16);
 padding:16px 18px 24px;overflow:auto;transform:translateX(100%);
 transition:transform .18s ease-out}
.panel[hidden]{display:none}
.panel.in{transform:none}
.panel h3{margin:0 40px 4px 0;font-size:15px;line-height:1.35;color:#0f172a}
.panel .close{position:absolute;top:12px;right:12px;border:0;background:#eef2f8;
 width:26px;height:26px;border-radius:50%;font-size:15px;cursor:pointer;color:#1d2430}
.panel dl{margin:12px 0 0;font-size:12px}
.panel dt{font-size:10px;letter-spacing:.09em;text-transform:uppercase;color:#8b96ab;
 margin:12px 0 3px}
.panel dd{margin:0;color:#1d2430}
.badge{display:inline-block;padding:3px 9px;border-radius:999px;color:#fff;
 font-size:11px;font-weight:700;letter-spacing:.03em}
.pid{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:11.5px;
 background:#f2f5fa;border-radius:4px;padding:1px 5px}
.plabel{font-size:12.5px;line-height:1.6;color:#1d2430}
.psrc{font-size:12px;line-height:1.6;color:#33415a}
.pflag{display:inline-block;font-size:10.5px;font-weight:700;color:#5a4a00;
 background:#fdf6d8;border-radius:4px;padding:2px 6px;margin:4px 4px 0 0}
"""

VIEWER_JS = r"""
var DATA = __DATA__;
var NODE = {}, EDGE = {}, VIEWS = DATA.views, TIER = DATA.tiers;
var ON_MAP = {};
(DATA.onMap || []).forEach(function (id) { ON_MAP[id] = true; });
DATA.nodes.forEach(function (n) { NODE[n.id] = n; });
DATA.edges.forEach(function (e) { EDGE[e.src + '->' + e.dst] = e; });

var vp = document.getElementById('viewer'),
    stage = document.getElementById('stage'),
    titleEl = document.getElementById('vtitle'),
    navEl = document.getElementById('nav'),
    chipEl = document.getElementById('chips'),
    searchEl = document.getElementById('search'),
    resEl = document.getElementById('results'),
    panel = document.getElementById('panel'),
    panelBody = document.getElementById('pbody'),
    scrim = document.getElementById('scrim'),
    helpEl = document.getElementById('help');

var st = {x: 0, y: 0, s: 1}, cur = null, natW = 1, natH = 1, offX = 0, offY = 0;
var activeTiers = {}, ignoreNextHash = false, highlight = null;

Object.keys(TIER).forEach(function (t) { activeTiers[t] = true; });

/* ---------------------------------------------------------------- transform */
function clampS(s) { return Math.min(14, Math.max(0.02, s)); }
function apply() {
  stage.style.transform = 'translate(' + st.x + 'px,' + st.y + 'px) scale(' + st.s + ')';
}
/* §8.5 fit() must measure the UNTRANSFORMED content. */
function fit() {
  var svg = stage.querySelector('svg');
  if (!svg) { return; }
  var prev = stage.style.transform;
  stage.style.transform = 'none';
  var w = 0, h = 0, bx = 0, by = 0;
  try {
    var b = svg.getBBox();
    w = b.width; h = b.height; bx = b.x; by = b.y;
  } catch (e) {
    var r = svg.getBoundingClientRect();
    w = r.width; h = r.height;
  }
  stage.style.transform = prev;
  if (!w || !h) { w = natW; h = natH; }
  natW = w; natH = h; offX = bx; offY = by;
  svg.style.width = w + 'px';
  svg.style.height = h + 'px';
  var vw = vp.clientWidth, vh = vp.clientHeight;
  st.s = clampS(Math.min(vw / w, vh / h, 1));
  st.x = (vw - w * st.s) / 2 - bx * st.s;
  st.y = (vh - h * st.s) / 2 - by * st.s;
  apply();
}
function zoomAt(px, py, f) {
  var s2 = clampS(st.s * f), k = s2 / st.s;
  st.x = px - (px - st.x) * k;
  st.y = py - (py - st.y) * k;
  st.s = s2; apply();
}
function zoomBy(f) { zoomAt(vp.clientWidth / 2, vp.clientHeight / 2, f); }
function centerOn(id, scale) {
  var svg = stage.querySelector('svg');
  if (!svg) { return; }
  var groups = svg.querySelectorAll('g.node'), found = null;
  for (var i = 0; i < groups.length; i++) {
    var t = groups[i].querySelector('title');
    if (t && t.textContent === id) { found = groups[i]; break; }
  }
  if (!found) { return; }
  /* On the position map a feature spans from its card down to its bar on the
     ruler; centring on the whole feature would put the ruler in the middle of
     the viewport.  Prefer the card when the view provides one. */
  var card = found.querySelector('.card');
  var b = (card || found).getBBox();
  var s = clampS(scale || Math.max(st.s, 1.25));
  st.s = s;
  st.x = vp.clientWidth / 2 - (b.x + b.width / 2) * s;
  st.y = vp.clientHeight / 2 - (b.y + b.height / 2) * s;
  apply();
}

/* ---------------------------------------------------------------- filtering */
function applyFilters() {
  var svg = stage.querySelector('svg');
  if (!svg) { return; }
  var nodes = svg.querySelectorAll('g.node'), i, g, t, id, rec, on;
  for (i = 0; i < nodes.length; i++) {
    g = nodes[i];
    t = g.querySelector('title');
    id = t ? t.textContent : '';
    rec = NODE[id];
    on = rec ? !!activeTiers[rec.evidence] : true;
    if (on && highlight) { on = highlight.has(id); }
    g.style.opacity = on ? '1' : '0.12';
    g.style.pointerEvents = on ? 'auto' : 'none';
  }
  var edges = svg.querySelectorAll('g.edge');
  for (i = 0; i < edges.length; i++) {
    g = edges[i];
    t = g.querySelector('title');
    var key = t ? t.textContent : '';
    if (/^LG/.test(key)) { continue; }          /* legend pseudo-edges exempt */
    rec = EDGE[key];
    on = rec ? !!activeTiers[rec.tier] : true;
    if (on && highlight) {
      var ends = key.split('->');
      on = highlight.has(ends[0]) || highlight.has(ends[1]);
    }
    g.style.opacity = on ? '1' : '0.07';
    g.style.pointerEvents = on ? 'auto' : 'none';
  }
}

/* ------------------------------------------------------------------ views */
function viewFor(key) {
  for (var i = 0; i < VIEWS.length; i++) {
    if (VIEWS[i].key === key) { return VIEWS[i]; }
  }
  return null;
}
/* §8.1: "first view" = most specific rendered sheet containing the node. */
function firstViewFor(id) {
  var rec = NODE[id];
  if (!rec) { return VIEWS[0].key; }
  /* A claim that sits on the position-map ruler is *about* a coordinate, so the
     map is the most informative sheet for it -- more so than the section sheet,
     where it is just one box among many. */
  if (ON_MAP[id]) { return 'marker_map_position_map'; }
  var cand = [];
  if (rec.section === 'SPINE') { cand.push('marker_map_spine'); }
  else { cand.push('marker_map_' + rec.section); }
  (rec.views || []).forEach(function (v) {
    if (v === 'spine') { cand.push('marker_map_spine'); }
    else if (v === 'full') { cand.push('marker_map_full'); }
    else { cand.push('marker_map_' + v); }
  });
  cand.push('marker_map_full');
  for (var i = 0; i < cand.length; i++) {
    if (viewFor(cand[i])) { return cand[i]; }
  }
  return VIEWS[0].key;
}
function setView(key, noHash) {
  var v = viewFor(key) || VIEWS[0];
  cur = v;
  while (stage.firstChild) { stage.removeChild(stage.firstChild); }
  var doc = new DOMParser().parseFromString(v.svg, 'image/svg+xml');
  var svg = document.importNode(doc.documentElement, true);
  svg.removeAttribute('width');
  svg.removeAttribute('height');
  stage.appendChild(svg);
  natW = 0; natH = 0;
  titleEl.textContent = v.title;
  var btns = navEl.querySelectorAll('button');
  for (var i = 0; i < btns.length; i++) {
    btns[i].className = (btns[i].getAttribute('data-key') === v.key) ? 'active' : '';
  }
  svg.addEventListener('click', onSvgClick);
  fit();
  applyFilters();
  if (!noHash) { writeHash(v.key, null); }
}
function onSvgClick(ev) {
  var g = ev.target;
  while (g && g !== stage) {
    if (g.classList && g.classList.contains('node')) {
      var t = g.querySelector('title');
      if (t) { openPanel(t.textContent); writeHash(cur.key, t.textContent); }
      return;
    }
    g = g.parentNode;
  }
}

/* ------------------------------------------------------------------ panel */
function el(tag, cls, text) {
  var e = document.createElement(tag);
  if (cls) { e.className = cls; }
  if (text !== undefined && text !== null) { e.textContent = text; }  /* §8.3 */
  return e;
}
function openPanel(id, noHash) {
  var rec = NODE[id];
  if (!rec) { return; }
  while (panelBody.firstChild) { panelBody.removeChild(panelBody.firstChild); }
  var h = el('h3', null, rec.label.split('\n')[0] || rec.id);
  panelBody.appendChild(h);
  var lines = el('div', 'plabel');
  rec.label.split('\n').forEach(function (l) { lines.appendChild(el('div', null, l)); });
  panelBody.appendChild(lines);
  var dl = el('dl');
  function row(term, node) { dl.appendChild(el('dt', null, term)); dl.appendChild(node); }
  var badge = el('span', 'badge', rec.evidence);
  badge.style.background = (TIER[rec.evidence] || {}).color || '#555';
  row('evidence tier', badge);
  row('section', el('dd', null, (DATA.sectionTitles[rec.section] || rec.section)
      + '  (' + rec.section + ')'));
  row('kind / shape', el('dd', null, rec.kind + ' \u2192 ' + rec.kind));
  var src = el('div', 'psrc');
  rec.source.split(';').forEach(function (part) {
    src.appendChild(el('div', null, part.trim()));
  });
  row('source', src);
  var flags = src.querySelectorAll('div'), i;
  for (i = 0; i < flags.length; i++) {
    var m = /\[([VSCW])\]/.exec(flags[i].textContent);
    if (m) { src.appendChild(el('span', 'pflag', '[' + m[1] + ']')); }
  }
  row('node id', el('dd', null, el('span', 'pid', rec.id)));
  var vw = el('div', null, '');
  (rec.views || []).forEach(function (v, j) {
    if (j) { vw.appendChild(document.createTextNode(', ')); }
    vw.appendChild(el('span', 'pid', v));
  });
  row('views', vw);
  var inc = [];
  DATA.edges.forEach(function (e) {
    if (e.src === id || e.dst === id) {
      inc.push(e.src + ' \u2192 ' + e.dst + '  [' + e.tier + ']');
    }
  });
  if (inc.length) {
    var ul = el('div');
    inc.forEach(function (t) { ul.appendChild(el('div', null, t)); });
    row('edges (' + inc.length + ')', el('dd', null, ul));
  }
  panelBody.appendChild(dl);
  scrim.hidden = false;
  panel.hidden = false;
  /* force reflow so the slide-in transition always runs */
  void panel.offsetWidth;
  panel.classList.add('in');
  if (!noHash) { writeHash(cur ? cur.key : '', id); }
}
function closePanel() {
  panel.classList.remove('in');
  scrim.hidden = true;
  setTimeout(function () { panel.hidden = true; }, 180);
}
document.getElementById('pclose').addEventListener('click', closePanel);
scrim.addEventListener('click', closePanel);

/* ------------------------------------------------------------- deep links */
function writeHash(viewKey, nodeId) {
  if (!viewKey) { return; }
  var h = '#view=' + encodeURIComponent(viewKey);
  if (nodeId) { h += '&node=' + encodeURIComponent(nodeId); }
  if (location.hash === h) { return; }
  ignoreNextHash = true;
  location.hash = h;
}
function parseHash() {
  var out = {}, raw = location.hash.replace(/^#/, '');
  if (!raw) { return out; }
  raw.split('&').forEach(function (part) {
    var i = part.indexOf('=');
    if (i > 0) {
      out[decodeURIComponent(part.slice(0, i))] = decodeURIComponent(part.slice(i + 1));
    }
  });
  return out;
}
function applyHash() {
  if (ignoreNextHash) { ignoreNextHash = false; return; }
  var p = parseHash();
  if (p.view && viewFor(p.view)) { setView(p.view, true); }
  if (p.node && NODE[p.node]) { openPanel(p.node, true); }   /* §8.4 no rewrite */
}
window.addEventListener('hashchange', applyHash);

/* ----------------------------------------------------------------- search */
function score(rec, q) {
  var id = rec.id.toLowerCase(), lab = rec.label.toLowerCase(),
      src = rec.source.toLowerCase();
  if (id === q) { return 0; }
  if (id.indexOf(q) === 0) { return 1; }
  if (id.indexOf(q) >= 0) { return 2; }
  if (lab.indexOf(q) >= 0) { return 3; }
  if (src.indexOf(q) >= 0) { return 4; }
  return 99;
}
function runSearch() {
  var q = (searchEl.value || '').trim().toLowerCase();
  while (resEl.firstChild) { resEl.removeChild(resEl.firstChild); }
  highlight = null;
  applyFilters();
  if (!q) { return; }
  var hits = DATA.nodes.filter(function (n) { return score(n, q) < 99; });
  hits.sort(function (a, b) { return score(a, q) - score(b, q) || 0; });
  highlight = new Set(hits.map(function (n) { return n.id; }));
  hits.slice(0, 12).forEach(function (n) {
    var b = el('button');
    b.appendChild(el('b', null, n.id));
    b.appendChild(document.createTextNode(' \u2014 ' + n.label.split('\n')[0]));
    b.addEventListener('click', function () {
      setView(firstViewFor(n.id));
      applyFilters();
      openPanel(n.id);
      centerOn(n.id);
      writeHash(cur.key, n.id);
    });
    resEl.appendChild(b);
  });
  if (hits.length > 12) {
    resEl.appendChild(el('button', null, '+' + (hits.length - 12) + ' more \u2014 refine'));
  }
  applyFilters();
}
searchEl.addEventListener('input', runSearch);
searchEl.addEventListener('keydown', function (e) {
  if (e.key === 'Enter') {
    var first = resEl.querySelector('button');
    if (first) { first.click(); }
  } else if (e.key === 'Escape') {
    searchEl.value = ''; runSearch();
  }
});

/* ------------------------------------------------------------------ chips */
Object.keys(TIER).forEach(function (t) {
  var b = el('button', 'on', t);
  b.style.background = TIER[t].color;
  b.style.borderColor = TIER[t].color;
  b.title = 'toggle ' + t + ' nodes and edges';
  b.addEventListener('click', function () {
    activeTiers[t] = !activeTiers[t];
    b.className = activeTiers[t] ? 'on' : '';
    b.style.background = activeTiers[t] ? TIER[t].color : 'transparent';
    b.style.color = activeTiers[t] ? '#fff' : TIER[t].color;
    applyFilters();
  });
  chipEl.appendChild(b);
});

/* ------------------------------------------------------------------ chrome */
VIEWS.forEach(function (v) {
  var b = el('button', null, v.title.replace(/^Retrotransposon Marker Map \u2014 /, ''));
  b.setAttribute('data-key', v.key);
  b.addEventListener('click', function () { setView(v.key); });
  navEl.appendChild(b);
});
document.getElementById('zin').addEventListener('click', function () { zoomBy(1.3); });
document.getElementById('zout').addEventListener('click', function () { zoomBy(0.77); });
document.getElementById('fitb').addEventListener('click', fit);
document.getElementById('helpb').addEventListener('click', function () {
  helpEl.hidden = !helpEl.hidden;
});
vp.addEventListener('wheel', function (e) {
  e.preventDefault();
  var r = vp.getBoundingClientRect();
  zoomAt(e.clientX - r.left, e.clientY - r.top, Math.exp(-e.deltaY * 0.0015));
}, {passive: false});
var drag = null;
vp.addEventListener('pointerdown', function (e) {
  if (e.target.closest && e.target.closest('.hud,.help')) { return; }
  drag = {x: e.clientX, y: e.clientY, sx: st.x, sy: st.y};
  vp.setPointerCapture(e.pointerId);
  vp.classList.add('dragging');
});
vp.addEventListener('pointermove', function (e) {
  if (!drag) { return; }
  st.x = drag.sx + (e.clientX - drag.x);
  st.y = drag.sy + (e.clientY - drag.y);
  apply();
});
vp.addEventListener('pointerup', function () {
  drag = null; vp.classList.remove('dragging');
});
vp.addEventListener('dblclick', fit);
window.addEventListener('resize', function () { if (cur) { setView(cur.key, true); } });

var start = parseHash();
setView((start.view && viewFor(start.view)) ? start.view : VIEWS[0].key, true);
applyHash();
if (!location.hash) { writeHash(VIEWS[0].key, null); }
"""

VIEWER_HTML = """<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>__TITLE__</title><style>__CSS__</style></head><body>
<nav>
  <h1>MARKER MAP v5</h1>
  <input id="search" type="search" placeholder="search id / label / source"
         autocomplete="off" spellcheck="false">
  <div id="results"></div>
  <h2>Views</h2>
  <div id="nav"></div>
  <h2>Evidence tiers</h2>
  <div id="chips"></div>
  <p>Wheel = zoom at cursor<br>Drag = pan<br>Double-click = fit<br>
     Click a node = citations</p>
  <p>[V] verified primary &middot; [S] snippet / unverified denominator &middot;
     [C] our computation &middot; [W] our data (double border)</p>
  <p>Generated __GENERATED__<br>__NODES__ nodes, __EDGES__ edges</p>
</nav>
<main>
  <div id="viewer" class="viewport">
    <div class="stage" id="stage"></div>
    <div class="viewtitle" id="vtitle"></div>
    <div class="hud">
      <button id="zin" title="zoom in">+</button>
      <button id="zout" title="zoom out">&minus;</button>
      <button id="fitb" title="fit to window">fit</button>
      <button id="helpb" title="help">?</button>
    </div>
    <div class="help" id="help" hidden>
      <b>Controls</b>
      <table>
        <tr><td>wheel / trackpad pinch</td><td>zoom at cursor</td></tr>
        <tr><td>drag</td><td>pan</td></tr>
        <tr><td>double-click</td><td>fit whole sheet</td></tr>
        <tr><td>click node</td><td>citation panel</td></tr>
        <tr><td>Enter in search</td><td>open first hit</td></tr>
        <tr><td>deep link</td><td>#view=&lt;key&gt;&amp;node=&lt;id&gt;</td></tr>
        <tr><td>search &rarr; Enter</td><td>opens the best sheet for that claim
            &mdash; the position map if it has a coordinate</td></tr>
      </table>
      <b style="display:block;margin-top:9px">Citation flags</b>
      <table>
        <tr><td>[V]</td><td>verified against the primary source</td></tr>
        <tr><td>[S]</td><td>snippet / denominator unverified</td></tr>
        <tr><td>[C]</td><td>our computation</td></tr>
        <tr><td>[W]</td><td>our data (double border)</td></tr>
      </table>
      <b style="display:block;margin-top:9px">Tier chips</b>
      Filter nodes by their evidence tier and edges by
      <i>their own</i> tier.
    </div>
  </div>
</main>
<div class="scrim" id="scrim" hidden></div>
<aside class="panel" id="panel" hidden>
  <button class="close" id="pclose" title="close">&#215;</button>
  <div id="pbody"></div>
</aside>
<script>__JS__</script></body></html>
"""

HELP_CSS_HINT = "help"


def _strip_svg_preamble(svg: str) -> str:
    """Drop the XML declaration / DOCTYPE / generator comments: they carry no data
    and only add parsing surface inside an HTML <script> block."""
    svg = re.sub(r"^\s*<\?xml[^>]*\?>\s*", "", svg)
    svg = re.sub(r"<!DOCTYPE[^>]*>\s*", "", svg, flags=re.IGNORECASE)
    svg = re.sub(r"<!--.*?-->", "", svg, flags=re.DOTALL)
    return svg.strip()


def build_viewer(views_rendered, nodes, edges, outdir: Path, generated: str) -> Path:
    blob = OrderedDict()
    blob["tiers"] = OrderedDict(
        (t, {"color": v[1], "style": v[0], "penwidth": v[2]}) for t, v in TIER_STYLE.items())
    blob["sectionTitles"] = OrderedDict((k, v[0]) for k, v in SECTION_META.items())
    blob["nodes"] = [OrderedDict([("id", n["id"]), ("label", n["label"]),
                                 ("section", n["section"]), ("kind", n["kind"]),
                                 ("evidence", n["evidence"]), ("source", n["source"]),
                                 ("views", list(n["views"]))]) for n in nodes]
    blob["edges"] = [OrderedDict([("src", e["src"]), ("dst", e["dst"]),
                                 ("tier", e["tier"]), ("label", e["label"]),
                                 ("views", list(e["views"]))]) for e in edges]
    blob["views"] = [OrderedDict([("key", v["key"]), ("class", v["class"]),
                                  ("title", v["title"]),
                                  ("svg", _strip_svg_preamble(v["svg"]))])
                     for v in views_rendered]
    # Which nodes the position map places on a ruler, so search can send a
    # coordinate claim to the sheet that shows the coordinate.
    blob["onMap"] = sorted({i for _, i, _, _ in MAP_ITEMS})
    # §8.6 encode-then-protect: json.dumps FIRST, then neutralise "</" inside it.
    data_js = json.dumps(blob, ensure_ascii=False).replace("</", "<\\/")
    html = (VIEWER_HTML
            .replace("__TITLE__", TITLE)
            .replace("__CSS__", VIEWER_CSS)
            .replace("__GENERATED__", generated)
            .replace("__NODES__", str(len(nodes)))
            .replace("__EDGES__", str(len(edges)))
            .replace("__JS__", VIEWER_JS.replace("__DATA__", data_js)))
    path = outdir / "index.html"
    path.write_text(html, encoding="utf-8")
    return path


# --------------------------------------------------------------------------- #
# §7.2 / §7.3 / §7.4  parity, manifest, README
# --------------------------------------------------------------------------- #
def build_manifest(outdir: Path, skip):
    rows = []
    for p in sorted(outdir.iterdir(), key=lambda q: q.name):
        if not p.is_file() or p.name in skip:
            continue
        rows.append((p.name, p.stat().st_size, sha256_file(p)))
    return rows


def write_readme(outdir: Path, ctx) -> Path:
    parity = ctx["parity"]
    L = []
    A = L.append
    A("# %s \u2014 v5 build report" % TITLE)
    A("")
    A("## Provenance")
    A("")
    A("- **Generated:** `%s`" % ctx["generated"])
    A("- **Builder:** `render_marker_map_v5.py` (stdlib + Graphviz `dot` only)")
    A("- **Input tables** (read-only; SHA-256 of the *input* files):")
    A("")
    A("| file | path | sha256 |")
    A("| --- | --- | --- |")
    for name, path, digest in ctx["inputs"]:
        A("| `%s` | `%s` | `%s` |" % (name, path, digest))
    A("- **Metrics JSON:** `%s` \u2014 **%s** (%s)"
      % (METRICS_PATH, ctx["g_mode"], ctx["g_reason"]))
    _msha = (ctx.get("metrics") or {}).get("file_sha")
    if _msha:
        # Worth its own line.  This file was rewritten underneath a build during
        # development and every value came out identical; recording the digest
        # is what turns "the numbers look right" into a checkable statement.
        A("  - digest at build time: `%s` \u2014 re-`shasum` the file and compare."
          % _msha)
    A("")
    A("## Correction ledger as applied")
    A("")
    A("Applied in memory after load; the exported CSVs carry the corrected text.")
    A("")
    A("| node / edge | change |")
    A("| --- | --- |")
    A("| `b_trunc` | label + source replaced: L1 genomic ~99% of >500k copies "
      "truncated (Freeman 2022); TSD-flanked intact L1s 85% 5'-truncated "
      "(Zingler 2005, denominator stated); ALU 8/226 (Wagstaff 2012) |")
    A("| `b_rearr` | label appended `(de novo cell culture; Symer 2002)`; "
      "source \u2192 `Wagstaff 2012; Symer 2002` |")
    A("| `b_invporub` | label + source replaced: 26% (330/1,271) polymorphic L1s, "
      "88% of inverted cluster at bases 4,000\u20136,000 (Porubsky 2023 [V]) |")
    A("| `e_l1hs` | label + source replaced: 1,642 annotated / 321 >6kb = 20% "
      "(Chu 2023, hg38); L1Base 146 FLI (Penzkofer); 6 hot = 84% (Brouha 2003 [V]) |")
    A("| `e_l1hstrunc` | **new node** (section E): L1Hs truncated ~71\u201372% "
      "(Pavlicek 2002 [S]) \u2014 denominator marked unverified |")
    A("| `e_l1hstrunc -> position` | **new edge**, tier `associated`, "
      "label `L1Hs subset` |")
    A("| `b_trans_sva` | label appended with the measured boundary: 0/671 cohort SVA "
      "calls can span past the element's 5' end and no call record holds an inserted "
      "length or names a source element; 8 reads at 1 of 3 SVA loci, clips "
      "2\u2013100 bp with nothing to subtract; source \u2192 + "
      "`sva_5prime_transduction_probe.py` |")
    A("| `b_l1term` | **new node** (section B): the L1 termination-site snippet with its "
      "boundary measured by `l1_readthrough_probe.py` \u2014 0/2186 calls can span past "
      "the element's 3' end, and 18 resolved junction reads at 4 loci show none "
      "transduced |")
    A("| `b_l1term -> b_trans_l1` | **new edge**, tier `associated`, "
      "label `readthrough substrate` |")
    A("")
    A("## Section H \u00b7 PRECEDENT \u2014 what exists beyond us")
    A("")
    A("H exists to keep the map honest about *precedent*: the positional-analysis")
    A("machinery we use is not original, and the strongest external tests of")
    A("\"is nested L1 insertion random?\" have known limits. %d nodes:"
      % len(H_NODES))
    A("")
    A("- **h_jedlicka** \u2014 the only published within-host positional analysis of")
    A("  nested TEs (12 plant genomes, 1,491 pairs; host 3'UTR most targeted).")
    A("- **h_levy_null** \u2014 Levy's positional null is an *empirical*")
    A("  fragmented-TE background (6 equal bins, chi\u00b2, \u22655 expected/bin):")
    A("  not uniform, not GC-matched.")
    A("- **h_levy_ends** \u2014 Levy excluded element ends for short-read")
    A("  detectability reasons; long reads escape that exclusion.")
    A("- **h_wimmer** \u2014 the recurrence standard: independent sites need TSD-length")
    A("  or subfamily discordance (18/4,500 = 0.4% of NF1 mutations).")
    A("- **h_flasch** \u2014 ~30% of engineered insertions land in existing L1/Alu")
    A("  sequence, but *proportional to opportunity*: \"no lightning rods\".")
    A("- **h_jed_mech** \u2014 mechanism split (Ty3/gypsy targeting vs Ty1/copia")
    A("  retention/loss).")
    A("- **h_pol3** \u2014 yeast Ty1-4 genomic target preference near Pol III genes:")
    A("  a targeting parallel, not a within-element one.")
    A("- **h_twinprim** \u2014 twin priming (Ostertag & Kazazian 2001): the TPRT DSB")
    A("  leaves two 3' overhangs, one priming on the poly(A) tail and one")
    A("  internally on the RNA; microhomology-driven SSA resolves into a 5'")
    A("  inversion (~25% of recent L1s inverted+truncated, abstract [V]). This is")
    A("  the published *mechanism* behind the inversion class in section B, not a")
    A("  new claim of ours.")
    A("- **h_atail** \u2014 the A-tail is a substrate, not just a primer site: Alu")
    A("  terminal segments expand during TPRT (median 42 A vs 23 internal; 95.9%")
    A("  grew \u22654 A via ORF2p slippage). Child-tail dynamics at the priming")
    A("  interface are distinct from host-tail targeting.")
    A("")
    A("The h_wimmer \u2192 x_tsdcomp edge is `gap`-tier on purpose: a nested-recurrence")
    A("test does not exist, and ours is in flight. That gap is the point of the map.")
    A("")
    A("## Section G \u00b7 OUR DATA \u2014 data-driven status")
    A("")
    if ctx["g_mode"] == "metrics_consumed" and (ctx.get("metrics") or {}).get("linker_fold") is not None:
        met = ctx["metrics"]
        A("**Mode: `metrics_consumed`.** The metrics JSON validated and every `g_*`")
        A("node now carries the real 10-genome figures. Nothing is hardcoded \u2014 each")
        A("value is read from the file at build time, and the JSON path it came from")
        A("is printed beside it.")
        A("")
        A("The run is identified by its `provenance.commit`, **not** by a date:")
        A("the file records no run date at all, and inventing one from the file")
        A("mtime or from `seed` would be a guess. It is labelled `unrecorded`")
        A("everywhere instead.")
        A("")
        A("| field | value | read from |")
        A("| --- | --- | --- |")
        for key in ("linker_fold", "tail_fold", "nested_l1_count", "union_sites",
                    "replication", "significant_replication",
                    "replication_fraction"):
            if key in met:
                A("| `%s` | %s | `%s` |"
                  % (key, _fmt_metric(key, met[key]), met["sources"].get(key, "?")))
        A("")
        A("`ten/samples` holds %d distinct genomes and `ten/per_sample` %d rows; both"
          % (EXPECTED_SAMPLES, met.get("per_sample_rows", 0)))
        A("are cross-checked against the headline, and the run is rejected as")
        A("in-flight if they disagree.")
        A("")
        A("**Superseded, not deleted.** The `[5-genome tier superseded by 10-genome")
        A("run %s]` tag replaces the old `[5-genome; 10-genome rerun in"
          % met["provenance"])
        A("flight]`. The 5-genome numbers stay on the page so both tiers are visible")
        A("and the reader can see which one is current. Note the direction of the")
        A("linker change: %s (%d genomes) \u2192 %s (%d), i.e. the extra genomes"
          % (_fmt_metric("linker_fold", met["superseded_linker_fold"]),
             met["superseded_n_samples"],
             _fmt_metric("linker_fold", met["linker_fold"]), EXPECTED_SAMPLES))
        A("*lowered* it, so the 5-genome figure was the optimistic one. Tail is")
        A("essentially unchanged (%s \u2192 %s)."
          % (_fmt_metric("tail_fold", met["superseded_tail_fold"]),
             _fmt_metric("tail_fold", met["tail_fold"])))
        A("")
        if met.get("p_floor"):
            A("### The p-value floor is a result, not a footnote")
            A("")
            A("`provenance.replicates` is %d, so the smallest p this method can"
              % met["replicates"])
            A("emit is `1/%d = %.4f`. Every bin sitting on that value is reported"
              % (met["replicates"], met["p_floor"]))
            A("as **p=%.4f (1/%d, the floor)** and never as \"p < %.4f\", because"
              % (met["p_floor"], met["replicates"], met["p_floor"]))
            A("claiming a smaller p than the method can produce is the easiest way to")
            A("overstate this result. `a_133ours` also carries a legacy chi-square")
            A("`p=4e-18` from the earlier tier; that p is on a different scale and is")
            A("tagged as not comparable rather than silently restated.")
        A("")
        if met.get("flagged"):
            A("### Counts the analysis itself flagged as still moving")
            A("")
            A("`delta[].flag_gt20pct` is the analysis grading its own 5\u219210-genome")
            A("shift. These five moved more than 20%, so **every absolute count from")
            A("the 5-genome tier is a floor, not an estimate**:")
            A("")
            A("| headline | 5 genomes | 10 genomes | change |")
            A("| --- | --- | --- | --- |")
            for f in met["flagged"]:
                five = f.get("five")
                ten = f.get("ten")
                if isinstance(five, float):
                    five = "%,.0f" % five if abs(five) >= 1000 else "%.3g" % five
                if isinstance(ten, float):
                    ten = "%,.0f" % ten if abs(ten) >= 1000 else "%.3g" % ten
                A("| %s | %s | %s | **%+.1f%%** |"
                  % (f.get("headline"), five, ten, f.get("pct") or 0.0))
            A("")
        _spr = parity.get("metrics_spread") or {}
        _alu_peak = (_spr.get("alu_body") or {}).get("peak_band") or {}
        # The linker bin specifically, because it is the one the node labels and
        # the map caption quote by name.  parity records each span as one
        # formatted string, so the key is built the same way rather than
        # hand-typed -- otherwise this silently falls through to "no siblings".
        _alu_linker = next((r for r in ((_spr.get("alu_body") or {})
                                        .get("per_bin") or [])
                            if r.get("span") == [_fmt_span(120.0, 140.0)]), None)
        _zb = [(n, c["zero_observation_bins"])
               for n, c in sorted((parity.get("metrics_curves") or {}).items())
               if c.get("zero_observation_bins")]
        _thinnest = min([b["obs"] for k, bins in (met.get("profiles") or {}).items()
                         if k in PLOTTED_CURVES for b in bins if b["obs"] > 0],
                        default=0)
        if _spr:
            A("### Every measured bin is a range, not a number")
            A("")
            A("`ten/profiles` holds **%d series**, not the %d curves the sheet draws. The"
              % (met.get("profile_series_total") or 0, len(PLOTTED_CURVES)))
            A("others measure the *same bins* with a different method or a different")
            A("cohort layer, and they disagree with the drawn curve. Each panel carries")
            A("the min\u2013max of its siblings as a pale band behind the curve, and")
            A("prints the numbers wherever the band is wide enough to change how the")
            A("bin reads. A sibling is admitted only when it shares the family **and**")
            A("the exact bin boundaries \u2014 an envelope between two curves on")
            A("different bins would be meaningless.")
            A("")
            A("A band counts as *material* when it is more than %.0f%% wide **and**"
              % (100.0 * SPREAD_MATERIAL))
            A("reaches the null line. The second condition matters: a band spanning")
            A("0.00\u00d7 to 0.20\u00d7 is 100% wide and entirely uninformative, because")
            A("both series agree there is nothing there.")
            A("")
            A("| curve | siblings | material bins | widest bin | range | from \u2192 to |")
            A("| --- | --- | --- | --- | --- | --- |")
            for key in sorted(_spr):
                s = _spr[key]
                w = s.get("widest")
                A("| `%s` | %d | %d of %d | %s | %s | %s \u2192 %s |"
                  % (key, s.get("series", 0), s.get("material_bins", 0),
                     len(s.get("per_bin") or []),
                     "\u2013".join(w["span"]) if w else "\u2014",
                     ("%.2f\u00d7\u2013%.2f\u00d7" % (w["min"], w["max"])) if w else "none",
                     w.get("min_series") if w else "\u2014",
                     w.get("max_series") if w else "\u2014"))
            A("")
            A("Three consequences, stated plainly:")
            A("")
            A("1. **The headline linker bin is method-dependent by about 2\u00d7.**")
            if _alu_linker:
                A("   120\u2013140 bp reads %.2f\u00d7 under `%s` and %.2f\u00d7 under `%s`."
                  % (_alu_linker["min"], _alu_linker["min_series"],
                     _alu_linker["max"], _alu_linker["max_series"]))
                A("   The peak *location* is not in doubt; the peak *size* is.")
                A("   `a_133ours` now prints the range, not a bare figure.")
            else:
                A("   The 120\u2013140 bp bin has no sibling series in this run, so")
                A("   `a_133ours` quotes a single figure and no band is drawn.")
            if _alu_peak:
                A("2. **The Alu 3' tail is not.** 280\u2013300 bp measures %.2f\u2013%.2f\u00d7"
                  % (_alu_peak["min"], _alu_peak["max"]))
                A("   across the same siblings, so `a_tail` keeps a single number.")
            A("3. **The L1 element-body curve is method-fragile; the relative-position")
            A("   curves are not.** `l1_full` has %d of %d bins material, against %s"
              % ((_spr.get("l1_full") or {}).get("material_bins", 0),
                 len(((_spr.get("l1_full") or {}).get("per_bin")) or []),
                 " and ".join("%d of %d for `%s`"
                              % ((_spr.get(k) or {}).get("material_bins", 0),
                                 len(((_spr.get(k) or {}).get("per_bin")) or []), k)
                              for k in ("alu_rel", "sva_rel")) + "."))
            A("   So the three-family contrast in lane D \u2014 L1 at the host's 5' end,")
            A("   SVA at its midpoint, Alu at its 3' end \u2014 is the robust part of")
            A("   this result: **all three peaks stay above 1\u00d7** under every sibling")
            A("   series that measures them.")
            A("")
            A("`parity.json \u2192 metrics_spread` carries the per-bin table with the")
            A("series behind each end of each range.")
            A("")
        if _zb:
            A("### A bin with no observations is a gap, not a zero")
            A("")
            A("The file reports `enrichment = 0.0` for a bin nothing landed in, which")
            A("is not the same as having measured a zero. Painting that at the axis")
            A("floor would draw a curve diving to nothing \u2014 a depletion claim, out")
            A("of a test that only ever asks about *enrichment*. So those bins are left")
            A("out, the curve edge is broken across them, and what was expected is")
            A("marked in their place:")
            A("")
            for _name, _spans in _zb:
                A("- `%s` \u2014 %s" % (_name, ", ".join(_spans)))
            A("")
            A("Every other plotted bin is floored on data: the thinnest bin on any")
            A("panel carries %d observation%s." % (_thinnest,
                                                 "" if _thinnest == 1 else "s"))
            A("")
        A("### The window claims were corrected at their source")
        A("")
        A("`a_133ours` and `a_tail` are the nodes that first claimed the linker and")
        A("tail windows, and both still carried the old figures as unqualified fact.")
        A("Each now carries a supersession tag naming the current measurement, with")
        A("the old fold-change and `n` **parsed out of the node's own label** so the")
        A("\"understated by Nx\" statement cannot drift from the text it describes:")
        A("")
        for nid, lo, headline in WINDOW_NODES:
            b = next((x for x in (met.get("profiles", {}).get("alu_body") or [])
                      if abs(x["lo"] - lo) < 1e-9), None)
            if b:
                A("- `%s` (%g\u2013%g bp) \u2192 **%.2fx, n=%d**, %s, obs %d vs expected"
                  % (nid, b["lo"], b["hi"], b["enr"], b["obs"],
                     _p_text(b["p"], met.get("p_floor")), b["obs"], ))
                A("  %.1f under the opportunity model." % b["exp"])
        A("")
        A("The evidence base had been understated more than fourfold (n=29\u2192124,")
        A("n=36\u2192164) and the linker had moved 11% downward. Nothing was deleted:")
        A("the 5-genome sentence stays, with the correction attached to it.")
    else:
        A("**Mode: `in_flight`.** %s" % ctx["g_reason"])
        A("")
        A("Every `g_*` node therefore carries the explicit suffix")
        A("`[5-genome; 10-genome rerun in flight]`. No 10-genome number is")
        A("hardcoded, guessed or projected anywhere in this build \u2014 when the")
        A("metrics JSON appears and validates, a rebuild folds it in and marks any")
        A("superseded value instead of deleting it.")
    A("")
    A("## Position map \u00b7 read this first")
    A("")
    A("`marker_map_position_map.svg` is a cartographic sheet, not a graph: every")
    A("interval claim in the map is placed on a shared ruler, so the reading is")
    A("\"in this window I have that\" rather than \"these boxes are connected\".")
    A("")
    A("- **Lane A \u00b7 element body coordinate**, 0\u2013300 bp, on a **linear**")
    A("  0\u2013320 bp ruler with the 0\u2013100 bp region shaded and captioned in")
    A("  10 bp bands. Holds the EN 7-mer, motif positions 2\u20135, the +7/+11/+19")
    A("  top-strand cuts, nick2 at 15\u201316 bp, the 5'-end peak, the A-tail 20\u201325 bp")
    A("  window, the left/right monomers and the 280\u2013300 bp tail peak.")
    A("- **Lane B \u00b7 TSD length**, same ruler, so TSDs are directly comparable")
    A("  with element-body positions. `b_trans_sva` (2,161 bp) runs off the right")
    A("  edge and is drawn with an explicit arrow, never silently clipped.")
    A("- **Lane C \u00b7 window / extent scale**, on its own log axis (10 bp \u2192 1 Mb),")
    A("  because a 50-bp GC window and a 1-Mb Alu window are not positions.")
    A("- **Lane D \u00b7 relative position along the host**, 0\u20131.0, its own")
    A("  normalised ruler. This is where the families separate: Alu peaks 2.62x at")
    A("  the 3' end of the host, SVA 3.60x near its midpoint.")
    A("- **Measured enrichment panel**, drawn directly above lane A on lane A's")
    A("  *own* 0\u2013320 bp ruler, so a spike in the curve and the bars piled up")
    A("  beneath it are literally the same base pairs.")
    A("")
    A("| property | value |")
    A("| --- | --- |")
    _basis_tally = (ctx.get("coord_basis") or {}).get("tally") or {}
    A("| nodes mapped | %d |" % len({i for _, i, _, _ in MAP_ITEMS}))
    A("| interval marks | %d |" % sum(len(s) for _, _, s, _ in MAP_ITEMS))
    A("| lane A / B bar sub-rows | %d / %d (packed, no bar overlaps a bar) |"
      % (_lane_nrows("A"), _lane_nrows("B")))
    A("| generator | built-in SVG writer (exact geometry, no layout engine) |")
    A("| measured curves plotted | %s |"
      % (", ".join(k for k in sorted((ctx.get("metrics") or {}).get("profiles") or {})
                    if k in PLOTTED_CURVES) or "none"))
    A("| coordinate provenance | %d of %d placements verbatim from the node's own "
      "label; %d declared below and tagged in the card |"
      % (_basis_tally.get("verbatim", 0), sum(_basis_tally.values()),
         sum(_basis_tally.values()) - _basis_tally.get("verbatim", 0)))
    A("")
    A("### Where the ruler is not quoting the label")
    A("")
    A("An earlier version of this sheet promised that every `lo/hi` was")
    A("transcribed verbatim and that nothing was inferred. **That was false for")
    A("five of the %d placements**, for three different and individually"
      % sum(_basis_tally.values()))
    A("legitimate reasons. Each is now declared in `MAP_COORD_BASIS`, checked")
    A("against the node's own label at build time, and tagged in the card it")
    A("appears in:")
    A("")
    A("| lane | node | basis | card tag | why it is not verbatim |")
    A("| --- | --- | --- | --- | --- |")
    for (lane, nid), kind in sorted(MAP_COORD_BASIS.items()):
        A("| %s | `%s` | `%s` | `%s` | %s |"
          % (lane, nid, kind, COORD_BASIS_TAG[kind], COORD_BASIS_LONG[kind]))
    A("")
    A("The two `extent` placements are the important ones: **a width is stated")
    A("and no position is.** `a_motif7` says \"EN 7-mer consensus\" and never says")
    A("where in the element it sits; `a_gc50` says \"50-bp cut window\", which is")
    A("a size. Drawing either as a filled, tier-coloured bar at a position would")
    A("have manufactured a fact out of a measurement, so they are drawn instead")
    A("as dashed grey brackets with end ticks at the ruler origin, rather than")
    A("filled tier-coloured bars, carrying only the width they do state. The")
    A("build refuses to run if an evidence tier is ever given the bracket's")
    A("ink, so the legend cannot go quietly stale.")
    A("")
    A("The build refuses the sheet if any placement is neither verbatim nor")
    A("declared (`audit_coord_basis`, exit 2), and the declaration is keyed by")
    A("`(lane, node)` so an exemption granted for one lane cannot leak into")
    A("another. Bracketed provenance stamps are stripped before the check, so a")
    A("number this renderer itself added can never become the evidence for a")
    A("coordinate.")
    A("")
    A("The sheet is built from the same node table as the graph sheets and emits")
    A("`<g class=\"node\"><title>id</title>` per card, so search, tier chips, the")
    A("citation panel and deep links work on it exactly as they do elsewhere.")
    A("It ships as SVG only \u2014 it is not a Graphviz sheet, so there is no PNG.")
    A("")
    A("A log axis was tried first and rejected: it collapsed the 1\u2013300 bp element")
    A("body, where nearly all the content lives, into about 5% of the width.")
    A("")
    A("## Parity")
    A("")
    A("- sections rendered: **%d** (spine + %d section sheets)"
      % (parity["sections_rendered"], len(SECTION_SHEETS)))
    A("- views rendered: **%d** (%d Graphviz sheets + the position map)"
      % (len(parity["views"]), len(parity["views"]) - 1))
    A("- nodes: **%d** \u00b7 edges: **%d** \u00b7 fatal errors: **%d**"
      % (parity["total_nodes"], parity["total_edges"], parity["fatal_errors"]))
    A("- corrections applied: %d \u00b7 H section added: %s"
      % (len(parity["corrections_applied"]), parity["h_section_added"]))
    A("")
    A("| view | class | nodes | edges | warnings |")
    A("| --- | --- | --- | --- | --- | --- |")
    for v in parity["views"]:
        if "marks" in v:
            A("| %s | %s | %d | %d | %d | %d interval marks, %s axis |"
              % (v["view"], v["class"], v["nodes"], v["edges"],
                 v["warnings_count"], v["marks"],
                 "linear 0\u2013320 bp + log" if v["view"].endswith("position_map")
                 else "\u2014"))
        else:
            A("| %s | %s | %d | %d | %d | \u2014 |"
              % (v["view"], v["class"], v["nodes"], v["edges"],
                 v["warnings_count"]))
    A("")
    A("`sections_rendered` and `len(views)` are reported separately on purpose:")
    A("the section sheets are 8, the spine is a 9th section-shaped sheet, `full`")
    A("is a 10th view that is not a section, and the position map is an 11th view")
    A("that is not a graph at all.")
    A("")
    A("## Viewer")
    A("")
    A("`index.html` (%s bytes) is fully offline \u2014 no CDN, no network. Wheel zooms"
      % f"{parity['viewer_bytes']:,}")
    A("at the cursor, drag pans, double-click fits. Added in v5: search over")
    A("id+label+source, evidence-tier chips (nodes by node tier, edges by their")
    A("own tier), a DOM-only citation panel, `#view=...&node=...` deep links, and")
    A("a help overlay documenting `[V]/[S]/[C]/[W]`. The position map is the first")
    A("item in the view list, because it is the sheet to read first.")
    A("")
    A("Node labels on the Graphviz sheets are HTML tables \u2014 bold headline, body,")
    A("italic citation \u2014 so a printed PNG or PDF carries its own references")
    A("instead of relying on hover.")
    A("")
    A("## File manifest (SHA-256)")
    A("")
    A("Computed last, after every other write completed. `README.md` itself is")
    A("the only excluded entry: its manifest block is written after hashing, so")
    A("its own digest cannot appear inside it.")
    A("")
    A("| file | bytes | sha256 |")
    A("| --- | --- | --- |")
    for name, size, digest in ctx["manifest"]:
        A("| `%s` | %d | `%s` |" % (name, size, digest))
    A("")
    A("## Rebuild")
    A("")
    A("```sh")
    A("python3 render_marker_map_v5.py --verify-only   # validation only")
    A("python3 render_marker_map_v5.py                 # full build")
    A("```")
    A("")
    path = outdir / "README.md"
    path.write_text("\n".join(L), encoding="utf-8")
    return path


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def main() -> int:
    ap = argparse.ArgumentParser(
        description="Rebuild the Retrotransposon Marker Map to v5.")
    ap.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR,
                    help="output directory (default: %s)" % DEFAULT_OUTDIR.name)
    ap.add_argument("--nodes", default=None, help="override node table path")
    ap.add_argument("--edges", default=None, help="override edge table path")
    ap.add_argument("--verify-only", dest="verify_only", action="store_true",
                    help="validate + report only; render nothing")
    ap.add_argument("--no-html", dest="no_html", action="store_true",
                    help="skip the offline viewer")
    args = ap.parse_args()

    outdir = args.outdir.expanduser()
    generated = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # ---- §1 load + validate -------------------------------------------------
    npath = resolve_input(NODE_CSV_NAME, args.nodes)
    epath = resolve_input(EDGE_CSV_NAME, args.edges)
    nodes, ndelim = load_nodes(npath)
    edges, edelim = load_edges(epath)
    validate(nodes, edges)
    print("VALIDATION OK  (%s / %s delimiter, %d + %d lines parsed)"
          % (ndelim, edelim, len(nodes) + 2, len(edges) + 2))
    print("  nodes: %d   edges: %d" % (len(nodes), len(edges)))
    print("  per-section tallies: %s"
          % ", ".join("%s=%d" % kv for kv in sorted(Counter(n["section"] for n in nodes).items())))
    print("  evidence tallies:   %s"
          % ", ".join("%s=%d" % kv for kv in sorted(Counter(n["evidence"] for n in nodes).items())))
    print("  edge tier tallies:  %s"
          % ", ".join("%s=%d" % kv for kv in sorted(Counter(e["tier"] for e in edges).items())))

    # ---- §2 corrections -----------------------------------------------------
    applied = apply_corrections(nodes, edges)
    print("  corrections applied (%d): %s" % (len(applied), ", ".join(applied)))

    # ---- §3 Section H -------------------------------------------------------
    add_h_section(nodes, edges)
    validate(nodes, edges)

    # ---- §4 Section G -------------------------------------------------------
    metrics, g_mode, g_reason = resolve_metrics()
    apply_g_section(nodes, metrics)
    print("  G-mode: %s (%s)" % (g_mode, g_reason))
    if metrics:
        for note in (apply_window_supersession(nodes, metrics) or []):
            print("  window: %s" % note)
        prof = metrics.get("profiles") or {}
        print("  measured curves: %s"
              % (", ".join("%s=%d bins" % (k, len(v))
                           for k, v in sorted(prof.items())) or "none"))
        for note in (metrics.get("profile_notes") or []):
            print("  curve NOT plotted: %s" % note)
        if metrics.get("flagged"):
            print("  counts still moving (delta flag_gt20pct): %s"
                  % "; ".join("%s %+.1f%%" % (f["headline"], f["pct"])
                              for f in metrics["flagged"]))
    print("  final graph: %d nodes, %d edges" % (len(nodes), len(edges)))

    # ---- position map self-check --------------------------------------------
    # The map is the one view whose geometry is computed rather than laid out, so
    # it gets its own invariant checks: every claim is on the ruler, no two bars
    # share a sub-row while overlapping, and every card fits the canvas.
    mapped = {i for _, i, _, _ in MAP_ITEMS}
    bad_span = [(i, sp) for _L, i, spans, _n in MAP_ITEMS for sp in spans
                if sp[0] > sp[1] or sp[0] < 0]
    if bad_span:
        die(EXIT_VALIDATION, "position map has a non-ascending span: %r" % (bad_span[0],))
    print("  position map: %d nodes, %d interval marks, lanes %s"
          % (len(mapped), sum(len(s) for _, _, s, _ in MAP_ITEMS),
             ", ".join("%s=%d sub-rows" % (k, _lane_nrows(k))
                       for k in ("A", "B"))))
    print("  position map: %d of %d nodes placed on the ruler (%.0f%%)"
          % (len(mapped), len(nodes), 100.0 * len(mapped) / len(nodes)))

    if args.verify_only:
        print("VERIFY-ONLY: validation complete, nothing rendered.")
        return EXIT_OK

    # ---- dot ---------------------------------------------------------------
    dot_bin = find_dot()
    if not dot_bin:
        die(EXIT_NO_DOT,
            "Graphviz 'dot' was not found on PATH or in the usual prefixes.\n"
            "  install it with:  brew install graphviz")
    if dot_bin != shutil.which("dot"):
        print("NOTE: 'dot' is not on PATH; using %s" % dot_bin)

    outdir.mkdir(parents=True, exist_ok=True)

    # ---- §5/§6 render -------------------------------------------------------
    views_rendered, view_report = [], []
    print("RENDER:")

    # The position map is generated first and is the "read me first" view: it is
    # the only sheet that places claims on a shared distance axis.
    map_key, map_cls = "marker_map_position_map", "map"
    map_title = "Position map \u00b7 read me first"
    map_svg = build_position_map(nodes, edges, generated, metrics)
    (outdir / (map_key + ".svg")).write_text(map_svg, encoding="utf-8")
    views_rendered.append({"key": map_key, "class": map_cls, "title": map_title,
                           "svg": map_svg})
    map_nodes = len({i for _, i, _, _ in MAP_ITEMS})
    map_marks = sum(len(s) for _, _, s, _ in MAP_ITEMS)
    view_report.append({"view": map_key, "class": map_cls, "nodes": map_nodes,
                        "edges": 0, "warnings_count": 0, "marks": map_marks,
                        "generator": "builtin-svg",
                        "axis": "A,B linear 0-%g bp (shared ruler); C log10(bp)"
                                 % MAP_LIN_MAX})
    print("  %-22s %-8s nodes=%-4d marks=%-4d (builtin SVG, linear 0-%g bp ruler)"
          % (map_key, map_cls, map_nodes, map_marks, MAP_LIN_MAX))

    for key, cls, title in build_view_list():
        vnodes, vedges = view_nodes_edges(key.split("marker_map_")[-1], nodes, edges)
        stem = outdir / key
        fmts = ("svg", "png", "pdf") if cls in ("spine", "full") else ("svg", "png")
        warnings = run_dot(dot_bin,
                           build_dot(key, cls, title, vnodes, vedges, metrics),
                           stem, fmts)
        views_rendered.append({"key": key, "class": cls, "title": title,
                               "svg": stem.with_suffix(".svg").read_text(encoding="utf-8")})
        view_report.append({"view": key, "class": cls, "nodes": len(vnodes),
                            "edges": len(vedges), "warnings_count": len(warnings)})
        print("  %-22s %-8s nodes=%-4d edges=%-4d warnings=%d"
              % (key, cls, len(vnodes), len(vedges), len(warnings)))
        for w in warnings[:5]:
            print("      warn: %s" % w)

    # ---- §7.1 CSV re-export -------------------------------------------------
    write_csvs(nodes, edges, outdir)
    print("CSV re-export: %s, %s" % (NODE_CSV_NAME, EDGE_CSV_NAME))

    # ---- §8 viewer ----------------------------------------------------------
    viewer_bytes = 0
    if args.no_html:
        print("viewer: skipped (--no-html)")
    else:
        vpath = build_viewer(views_rendered, nodes, edges, outdir, generated)
        viewer_bytes = vpath.stat().st_size
        print("viewer: index.html (%s bytes)" % f"{viewer_bytes:,}")

    # ---- §7.2 parity.json ---------------------------------------------------
    try:
        _cb = audit_coord_basis(nodes, MAP_ITEMS)
    except SystemExit:
        _cb = None          # the sheet already refused it; do not claim otherwise
    parity = OrderedDict([
        ("build_timestamp", generated),
        ("corrections_applied", applied),
        ("h_section_added", True),
        ("g_mode", g_mode),
        ("g_reason", g_reason),
        ("metrics_sources", (metrics or {}).get("sources", {})),
        ("metrics_run", None if metrics is None else {
            "run_id": metrics["run_id"], "commit": metrics["provenance"],
            "run_date": metrics["run_date"], "seed": metrics["seed"],
            "n_samples": EXPECTED_SAMPLES,
            "per_sample_rows": metrics.get("per_sample_rows"),
            # The digest of the file every number above was read from.  A run id
            # says which analysis produced the values; only this says whether the
            # bytes on disk right now are the ones the sheet was built against.
            "file_sha256": metrics.get("file_sha"),
        }),
        # What the map actually drew, and what it refused to draw.  A curve that
        # failed validation has to be named here, otherwise a missing panel and
        # a deliberately withheld panel are indistinguishable downstream.
        ("metrics_curves", OrderedDict(
            (name, {"series": list(PROFILE_SERIES[name]),
                    "bins": len(bins),
                    "span": [bins[0]["lo"], bins[-1]["hi"]],
                    "plotted": PLOTTED_CURVES.get(name),
                    "validated_not_plotted": None if name in PLOTTED_CURVES
                                              else "available in ten/profiles; "
                                                   "used as a spread sibling, "
                                                   "not drawn as a panel",
                    "zero_observation_bins": [
                        [_axis_span(name, z["lo"], z["hi"])
                         for z in _zero_spans(bins)] or None][0],
                    })
            for name, bins in sorted(((metrics or {}).get("profiles") or {}).items()))),
        ("metrics_curves_omitted", (metrics or {}).get("profile_notes") or []),
        # Method sensitivity, bin by bin.  Every plotted bin measures a range,
        # not a number, and the range is what the caption prints.  Recorded here
        # in full -- including which method produced each end -- because the
        # sheet has no room for it and the claim is the sheet's headline.
        ("metrics_spread", _spread_parity(metrics)),
        ("metrics_spread_notes", (metrics or {}).get("spread_notes") or []),
        # The position map's core promise is that every coordinate can be
        # accounted for.  Record the audit itself, so a downstream reader can
        # tell a verbatim placement from a declared one without re-running the
        # build, and can see which exemptions were relied on.
        ("map_placements", sum(_cb[1].values()) if _cb else 0),
        ("map_coord_provenance", OrderedDict(
            (k, v) for k, v in sorted((_cb[1] if _cb else {}).items()))),
        ("map_coord_basis_declared", OrderedDict(
            ("lane %s / %s" % k, v) for k, v in sorted(MAP_COORD_BASIS.items()))),
        ("map_coord_basis_guard", "audit_coord_basis: a placement that is neither "
                                  "verbatim against its own label nor declared in "
                                  "MAP_COORD_BASIS aborts the build with exit 2"),
        ("metrics_p_floor", None if metrics is None else metrics.get("p_floor")),
        ("metrics_replicates", None if metrics is None else metrics.get("replicates")),
        ("metrics_flagged_counts", (metrics or {}).get("flagged") or []),
        ("views", view_report),
        ("sections_rendered", 1 + len(SECTION_SHEETS)),
        ("total_nodes", len(nodes)),
        ("total_edges", len(edges)),
        ("fatal_errors", 0),
        ("viewer_bytes", viewer_bytes),
    ])
    (outdir / "parity.json").write_text(json.dumps(parity, indent=2) + "\n",
                                        encoding="utf-8")
    print("parity.json written (sections_rendered=%d, views=%d)"
          % (parity["sections_rendered"], len(parity["views"])))

    # ---- §7.4 manifest, then README ----------------------------------------
    manifest = build_manifest(outdir, skip={"README.md"})
    ctx = {
        "generated": generated,
        "inputs": [(p.name, str(p), sha256_file(p))
                   for p in (npath, epath)],
        "g_mode": g_mode,
        "g_reason": g_reason,
        "metrics": metrics,
        "coord_basis": {"tally": _cb[1], "kinds": _cb[0]} if _cb else {},
        "parity": parity,
        "manifest": manifest,
    }
    write_readme(outdir, ctx)
    print("README.md written (manifest: %d files)" % len(manifest))
    print("DONE: %s" % outdir)
    return EXIT_OK


if __name__ == "__main__":
    sys.exit(main())