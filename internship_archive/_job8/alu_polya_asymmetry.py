"""Does the AluY consensus carry sense-strand poly(A) signals beyond composition?

Background
----------
Intronic Alu insertions are depleted in the sense orientation relative to their
host gene (Medstrand 2002; reproduced in this project at a sense fraction of
0.457, n=733 intronic Alu calls in HG03086). Han, Szak & Boeke 2004 explain the
effect for LINE-1 via transcriptional elongation, but their mechanism is
specific to L1 ORF sequence and length-dependent, so it does not cover Alu.

Conley & Jordan 2012 supply a candidate that does: TE sequences provide
transcription termination sites, those sites are more often sense-oriented, and
Alu is the single largest contributor. Roy-Engel 2005 locates the signal in the
A-rich tails of both SINEs and LINEs.

That last point is the problem this script exists to handle. AATAAA is an A-rich
hexamer, so an A-rich strand carries more of them *for composition reasons
alone*. The same test on L1.3 found its 8.1x sense-strand skew was entirely
composition (shuffle null 8.12x, empirical p=0.99). So a raw count proves
nothing.

Design
------
Three nested questions, each stricter than the last:

1. COUNT. Does the sense strand carry more poly(A) hexamers than a
   composition-matched null predicts? Two nulls: mononucleotide shuffle
   (composition only) and a first-order Markov chain (composition plus
   dinucleotide structure).

2. POSITION. Alu has two A-rich regions that will generate hexamers by
   composition: the internal A-rich linker between the two arms, and the 3'
   A-tail. Where do the observed hits actually sit?

3. SURVIVAL. Masking those two regions, does a sense-strand excess remain? If
   yes, Alu carries poly(A) signals in its body and the termination mechanism is
   general. If no, the orientation bias traces to the A-tail specifically, which
   is a narrower and more testable claim.

The consensus is built from the genome rather than downloaded: ungapped,
full-length, low-divergence AluY copies from the UCSC RepeatMasker track,
stacked on consensus coordinates and called by column majority.
"""
from __future__ import annotations

import argparse
import collections
import gzip
import re
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
try:
    from l1_strand_scan import (POLYA_HEXAMERS, POLYA_STRONG, fit_markov1,
                                revcomp, scan_polya, shuffle_mono,
                                simulate_markov1, strand_counts)
except ModuleNotFoundError as exc:  # pragma: no cover
    raise SystemExit(
        "l1_strand_scan.py must sit beside this script (it is in L1/artifacts/)"
    ) from exc

MAIN_CHROMS = {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]}
CONSENSUS_LEN = 311


# --- consensus construction ------------------------------------------------

def select_copies(rmsk_path, rep_name="AluY", max_milli_div=60,
                  min_rep_end=300, max_cons_start=2):
    """Ungapped, full-length, low-divergence copies, on consensus coordinates.

    UCSC orders the three consensus-coordinate columns differently by strand:
    for '-' rows `repLeft` holds the consensus start and `repStart` the negative
    remainder, the reverse of '+' rows. Getting this backwards inflates lengths
    (it cost 9 loci in an earlier analysis in this project), so the strand
    branch below is load-bearing.

    Indel-free copies only (milliDel == milliIns == 0) and genomic length equal
    to consensus length, so the stack needs no alignment.
    """
    kept = []
    opener = gzip.open if str(rmsk_path).endswith(".gz") else open
    with opener(rmsk_path, "rt") as handle:
        for line in handle:
            f = line.rstrip("\n").split("\t")
            if len(f) < 16 or f[10] != rep_name or f[5] not in MAIN_CHROMS:
                continue
            milli_div, milli_del, milli_ins = int(f[2]), int(f[3]), int(f[4])
            if milli_div > max_milli_div or milli_del or milli_ins:
                continue
            strand = f[9]
            rep_start, rep_end, rep_left = int(f[13]), int(f[14]), int(f[15])
            cons_start = rep_start if strand == "+" else rep_left
            if cons_start > max_cons_start or rep_end < min_rep_end:
                continue
            start, end = int(f[6]), int(f[7])
            if (end - start) != (rep_end - cons_start + 1):
                continue
            kept.append((f[5], start, end, strand, cons_start, rep_end))
    return kept


def select_all_copies(rmsk_path, rep_name="AluY", chroms=MAIN_CHROMS):
    """Every annotated copy of `rep_name`, with divergence. Part 2 wants all of
    them, not just the ungapped full-length subset the consensus needs."""
    out = []
    opener = gzip.open if str(rmsk_path).endswith(".gz") else open
    with opener(rmsk_path, "rt") as handle:
        for line in handle:
            f = line.rstrip("\n").split("\t")
            if len(f) < 16 or f[10] != rep_name or f[5] not in chroms:
                continue
            out.append((f[5], int(f[6]), int(f[7]), f[9], int(f[2])))
    return out


def build_consensus(copies, twobit_path, cons_len=CONSENSUS_LEN):
    import twobitreader
    genome = twobitreader.TwoBitFile(str(twobit_path))
    columns = [collections.Counter() for _ in range(cons_len)]
    used = 0
    for chrom, start, end, strand, cons_start, _rep_end in copies:
        seq = genome[chrom][start:end].upper()
        if strand == "-":
            seq = revcomp(seq)
        offset = cons_start - 1
        for i, base in enumerate(seq):
            j = offset + i
            if j < cons_len and base in "ACGT":
                columns[j][base] += 1
        used += 1
    consensus = "".join(max(c, key=c.get) if c else "N" for c in columns)
    depth = np.array([sum(c.values()) for c in columns])
    support = np.array([(max(c.values()) / sum(c.values())) if c else 0.0
                        for c in columns])
    return consensus, {"copies_used": used,
                       "min_depth": int(depth.min()),
                       "median_depth": int(np.median(depth)),
                       "mean_column_support": round(float(support.mean()), 4),
                       "columns_below_0.6_support": int((support < 0.6).sum())}


# --- A-rich regions --------------------------------------------------------

def a_rich_regions(seq, window=20, min_a_frac=0.6):
    """1-based inclusive spans whose A-fraction reaches `min_a_frac`.

    Found from the sequence rather than hardcoded, so the regions are a property
    of the consensus and move if the consensus does.
    """
    n = len(seq)
    flagged = np.zeros(n, dtype=bool)
    for i in range(n - window + 1):
        if seq[i:i + window].count("A") / window >= min_a_frac:
            flagged[i:i + window] = True
    spans, start = [], None
    for i, hot in enumerate(flagged):
        if hot and start is None:
            start = i
        elif not hot and start is not None:
            spans.append((start + 1, i))
            start = None
    if start is not None:
        spans.append((start + 1, n))
    return spans


def mask_spans(seq, spans):
    """Delete the spans. Deleting rather than N-filling keeps the null's
    composition equal to the masked sequence's own composition."""
    keep = np.ones(len(seq), dtype=bool)
    for lo, hi in spans:
        keep[lo - 1:hi] = False
    return "".join(b for b, k in zip(seq, keep) if k)


# --- the test ---------------------------------------------------------------

def null_distribution(seq, n_rep, rng, kind, canonical_only=False):
    sense = np.empty(n_rep, dtype=int)
    anti = np.empty(n_rep, dtype=int)
    if kind == "markov":
        states, p0, T = fit_markov1(seq)
    for i in range(n_rep):
        sim = (shuffle_mono(seq, rng) if kind == "mono"
               else simulate_markov1(len(seq), states, p0, T, rng))
        sense[i], anti[i] = strand_counts(sim, canonical_only)
    return sense, anti


def test_sequence(seq, label, n_rep, seed, canonical_only=False):
    obs_sense, obs_anti = strand_counts(seq, canonical_only)
    out = {"label": label, "length_bp": len(seq),
           "a_fraction": round(seq.count("A") / len(seq), 4),
           "canonical_only": canonical_only,
           "observed_sense": int(obs_sense), "observed_antisense": int(obs_anti),
           "observed_ratio": (round(obs_sense / obs_anti, 3) if obs_anti
                              else None),
           "n_replicates": n_rep, "nulls": {}}
    for kind in ("mono", "markov"):
        rng = np.random.default_rng(seed)   # same stream per null, reproducible
        ns, na = null_distribution(seq, n_rep, rng, kind, canonical_only)
        # one-sided: is the SENSE count higher than composition predicts?
        p_sense = (1 + int((ns >= obs_sense).sum())) / (n_rep + 1)
        ratios = ns / np.maximum(na, 1)
        obs_ratio = obs_sense / max(obs_anti, 1)
        p_ratio = (1 + int((ratios >= obs_ratio).sum())) / (n_rep + 1)
        out["nulls"][kind] = {
            "expected_sense_mean": round(float(ns.mean()), 2),
            "expected_sense_p95": float(np.percentile(ns, 95)),
            "expected_ratio_median": round(float(np.median(ratios)), 3),
            "p_sense_count_one_sided": round(p_sense, 5),
            "p_ratio_one_sided": round(p_ratio, 5),
        }
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("rmsk", type=Path, help="UCSC rmsk table (.txt.gz)")
    ap.add_argument("twobit", type=Path, help="hg38.2bit")
    ap.add_argument("--rep-name", default="AluY")
    ap.add_argument("--replicates", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=17)
    ap.add_argument("--a-window", type=int, default=20)
    ap.add_argument("--a-frac", type=float, default=0.6)
    ap.add_argument("--gencode", type=Path, default=None,
                    help="GENCODE GTF(.gz); enables the per-copy test (part 2)")
    ap.add_argument("--tail-start", type=int, default=270,
                    help="1-based consensus position where the 3' A-tail begins")
    ap.add_argument("--out-dir", type=Path, default=Path("alu_polya_out"))
    a = ap.parse_args(argv)

    a.out_dir.mkdir(parents=True, exist_ok=True)
    copies = select_copies(a.rmsk, a.rep_name)
    if not copies:
        raise SystemExit(f"no usable {a.rep_name} copies found in {a.rmsk}")
    consensus, qc = build_consensus(copies, a.twobit)
    spans = a_rich_regions(consensus, a.a_window, a.a_frac)
    masked = mask_spans(consensus, spans)

    report = {
        "rep_name": a.rep_name,
        "consensus": consensus,
        "consensus_qc": qc,
        "hexamers_scanned": POLYA_HEXAMERS,
        "strong_hexamers": sorted(POLYA_STRONG),
        "a_rich_regions_1based": [list(s) for s in spans],
        "a_rich_bp": sum(hi - lo + 1 for lo, hi in spans),
        "tests": [
            test_sequence(consensus, "full_consensus", a.replicates, a.seed),
            test_sequence(consensus, "full_consensus_canonical_only",
                          a.replicates, a.seed, canonical_only=True),
            test_sequence(masked, "a_rich_regions_removed", a.replicates, a.seed),
            test_sequence(masked, "a_rich_regions_removed_canonical_only",
                          a.replicates, a.seed, canonical_only=True),
        ],
    }

    # per-hit table, so position claims are auditable
    lines = ["strand\tpos_1based\tmotif\tin_a_rich_region"]
    for h in sorted(scan_polya(consensus), key=lambda x: (x.l1_pos, x.strand)):
        inside = any(lo <= h.l1_pos <= hi for lo, hi in spans)
        lines.append(f"{h.strand}\t{h.l1_pos}\t{h.motif}\t{int(inside)}")
    (a.out_dir / "polya_hits.tsv").write_text("\n".join(lines) + "\n")
    print(f"{a.rep_name} consensus {len(consensus)} bp from {qc['copies_used']} copies "
          f"(mean column support {qc['mean_column_support']})")
    print(f"A-rich regions (>={a.a_frac} A over {a.a_window} bp): {spans} "
          f"= {report['a_rich_bp']} bp")
    print()
    hdr = f"{'test':<42}{'sense':>6}{'anti':>6}{'ratio':>7}{'mono p':>9}{'markov p':>10}"
    print(hdr); print("-" * len(hdr))
    for t in report["tests"]:
        print(f"{t['label']:<42}{t['observed_sense']:>6}{t['observed_antisense']:>6}"
              f"{str(t['observed_ratio']):>7}"
              f"{t['nulls']['mono']['p_sense_count_one_sided']:>9.4f}"
              f"{t['nulls']['markov']['p_sense_count_one_sided']:>10.4f}")
    if a.gencode:
        gene_idx = load_gene_strands(a.gencode)
        all_copies = select_all_copies(a.rmsk, a.rep_name)
        per_copy, _df = per_copy_analysis(all_copies, a.twobit, gene_idx,
                                         a.out_dir, a.tail_start)
        report["per_copy"] = per_copy
        print()
        print(f"per-copy test: {per_copy['n_copies_in_single_strand_genes']:,} "
              f"{a.rep_name} copies inside single-strand genes, "
              f"overall sense fraction {per_copy['overall_sense_fraction']}")
        h = (f"{'hexamer location':<20}{'carriers':>9}{'sense':>8}{'non':>8}"
             f"{'crude OR':>10}{'adj OR':>8}{'95% CI':>16}{'p':>9}")
        print(h); print("-" * len(h))
        for label, e in per_copy["exposures"].items():
            adj = e["adjusted"] or {}
            ci = (f"{adj['ci95'][0]:.3f}-{adj['ci95'][1]:.3f}" if adj else "-")
            print(f"{label:<20}{e['n_carriers']:>9,}"
                  f"{e['carrier_sense_fraction']:>8.4f}{e['noncarrier_sense_fraction']:>8.4f}"
                  f"{e['crude_or']:>10.3f}{adj.get('or_mh', float('nan')):>8.3f}"
                  f"{ci:>16}{adj.get('p', float('nan')):>9.4f}")

    (a.out_dir / "alu_polya_asymmetry.json").write_text(
        json.dumps(report, indent=2) + "\n")
    print(f"\nwrote {a.out_dir}/alu_polya_asymmetry.json and polya_hits.tsv")
    return report




# ---------------------------------------------------------------------------
# Part 2: per-copy test in the reference genome
#
# The consensus test above answers "does the element, as a species, carry a
# signal". It does not, because AluY's consensus has no canonical hexamer on
# either strand and its 3' tail is a pure A-run, which cannot contain AATAAA.
# Roy-Engel 2005 says the signals are *created* in the A-rich tails, i.e. by
# mutation in individual copies. So the test has to be per copy.
#
# Design: among fixed reference AluY copies inside genes, compare the
# sense-orientation fraction of copies that carry a canonical hexamer on their
# own sense strand against copies that do not. A carrier can only terminate host
# transcription when it is sense-oriented, so if termination drives the
# orientation bias, carriers should be more depleted in sense than non-carriers.
# This is a within-Alu contrast: everything about being an Alu is held fixed.
#
# Confounding is the whole problem here. Carriers are older and longer than
# non-carriers (both raise the chance of a hexamer existing at all, and age also
# raises exposure to selection), so the crude comparison is not interpretable.
# Estimates below are Mantel-Haenszel, stratified on divergence and length.
# ---------------------------------------------------------------------------

STRONG_HEXAMERS = ("AATAAA", "ATTAAA")


def load_gene_strands(gtf_path, chroms=MAIN_CHROMS):
    """Per-chromosome, per-strand sorted gene-body intervals from a GENCODE GTF."""
    per = {}
    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    with opener(gtf_path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.split("\t", 9)
            if f[2] != "gene" or f[0] not in chroms:
                continue
            per.setdefault((f[0], f[6]), []).append((int(f[3]), int(f[4])))
    out = {}
    for key, spans in per.items():
        spans.sort()
        starts = np.array([s for s, _ in spans])
        ends = np.array([e for _, e in spans])
        out[key] = (starts, ends)
    return out


def in_gene(gene_idx, chrom, pos, strand):
    key = (chrom, strand)
    if key not in gene_idx:
        return False
    starts, ends = gene_idx[key]
    i = int(np.searchsorted(starts, pos, side="right"))
    return bool(i and (ends[:i] >= pos).any())


def count_hexamers(seq, lo=1, hi=None):
    """Canonical hexamers on `seq` as given, within 1-based [lo, hi]."""
    hi = len(seq) if hi is None else hi
    n = 0
    for hexamer in STRONG_HEXAMERS:
        for m in re.finditer(f"(?={hexamer})", seq):
            if lo <= m.start() + 1 <= hi:
                n += 1
    return n


def mantel_haenszel(table_iter):
    """MH odds ratio, RBG standard error, and the MH chi-square p-value.

    `table_iter` yields (a, b, c, d) = (exposed-sense, exposed-antisense,
    unexposed-sense, unexposed-antisense) per stratum.
    """
    from scipy.stats import chi2 as _chi2
    num = den = 0.0
    sa = sb = sc = 0.0
    obs = exp = var = 0.0
    strata = 0
    for a, b, c, d in table_iter:
        n = a + b + c + d
        if min(a, b, c, d) < 5 or n < 20:
            continue
        R, S = a * d / n, b * c / n
        num += R
        den += S
        sa += ((a + d) * R) / n
        sb += (((a + d) * S) + ((b + c) * R)) / n
        sc += ((b + c) * S) / n
        obs += a
        exp += (a + b) * (a + c) / n
        var += (a + b) * (c + d) * (a + c) * (b + d) / (n * n * (n - 1))
        strata += 1
    if not strata or den == 0 or var == 0:
        return None
    or_mh = num / den
    se = np.sqrt(sa / (2 * num ** 2) + sb / (2 * num * den) + sc / (2 * den ** 2))
    chi_sq = (abs(obs - exp) - 0.5) ** 2 / var
    return {"or_mh": round(float(or_mh), 4),
            "ci95": [round(float(np.exp(np.log(or_mh) - 1.96 * se)), 4),
                     round(float(np.exp(np.log(or_mh) + 1.96 * se)), 4)],
            "chi2": round(float(chi_sq), 3),
            "p": float(_chi2.sf(chi_sq, 1)),
            "strata": strata}

def per_copy_analysis(copies_all, twobit_path, gene_idx, out_dir,
                      tail_start=270):
    """Carrier-vs-non-carrier orientation test over reference copies in genes."""
    import twobitreader
    genome = twobitreader.TwoBitFile(str(twobit_path))
    cache = {}
    recs = []
    for chrom, start, end, strand, milli_div in copies_all:
        if chrom not in cache:
            cache[chrom] = genome[chrom]
        mid = (start + end) // 2 + 1
        plus = in_gene(gene_idx, chrom, mid, "+")
        minus = in_gene(gene_idx, chrom, mid, "-")
        if plus == minus:           # intergenic, or both strands: not assignable
            continue
        gene_strand = "+" if plus else "-"
        seq = cache[chrom][start:end].upper()
        if strand == "-":
            seq = revcomp(seq)
        recs.append({
            "chrom": chrom, "start": start, "end": end,
            "elem_strand": strand, "gene_strand": gene_strand,
            "orient": "sense" if strand == gene_strand else "antisense",
            "milli_div": milli_div, "length_bp": end - start,
            "n_hex_any": count_hexamers(seq),
            "n_hex_tail": count_hexamers(seq, tail_start),
            "n_hex_body": count_hexamers(seq, 1, tail_start - 1),
        })
    if not recs:
        raise SystemExit("no genic copies assignable to a single gene strand")

    import pandas as pd
    df = pd.DataFrame(recs)
    df["sense"] = df.orient.eq("sense")
    df["div_bin"] = pd.cut(df.milli_div, [-1, 40, 60, 80, 100, 10 ** 6],
                           labels=["<=40", "41-60", "61-80", "81-100", ">100"])
    df["len_bin"] = pd.cut(df.length_bp, [0, 290, 300, 310, 10 ** 6],
                           labels=["<=290", "291-300", "301-310", ">310"])
    df.to_csv(out_dir / "alu_per_copy.csv.gz", index=False, compression="gzip")

    results = {"n_copies_in_single_strand_genes": int(len(df)),
               "overall_sense_fraction": round(float(df.sense.mean()), 4),
               "covariate_means_by_carrier": {}, "exposures": {}}
    for col, label in [("n_hex_any", "any_position"),
                       ("n_hex_tail", "a_tail_only"),
                       ("n_hex_body", "element_body_only")]:
        exposed = df[col] > 0
        results["covariate_means_by_carrier"][label] = {
            "carrier_milli_div": round(float(df.loc[exposed, "milli_div"].mean()), 2),
            "noncarrier_milli_div": round(float(df.loc[~exposed, "milli_div"].mean()), 2),
            "carrier_length_bp": round(float(df.loc[exposed, "length_bp"].mean()), 2),
            "noncarrier_length_bp": round(float(df.loc[~exposed, "length_bp"].mean()), 2),
        }

        def tables():
            for _, g in df.groupby(["div_bin", "len_bin"], observed=True):
                e = g[col] > 0
                yield (int((e & g.sense).sum()), int((e & ~g.sense).sum()),
                       int((~e & g.sense).sum()), int((~e & ~g.sense).sum()))

        results["exposures"][label] = {
            "n_carriers": int(exposed.sum()),
            "carrier_sense_fraction": round(float(df.loc[exposed, "sense"].mean()), 4),
            "noncarrier_sense_fraction": round(float(df.loc[~exposed, "sense"].mean()), 4),
            "crude_or": round(float(
                (df.loc[exposed, "sense"].sum() * (~df.loc[~exposed, "sense"]).sum()) /
                max((~df.loc[exposed, "sense"]).sum() * df.loc[~exposed, "sense"].sum(), 1)), 4),
            "adjusted": mantel_haenszel(tables()),
        }
    return results, df


if __name__ == "__main__":
    main()
