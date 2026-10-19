"""Test whether the Jacob-Hirsch 2018 CCATT target motif is enriched at our nested
insertion breakpoints, relative to a control matched for element composition.

Design note, which is the whole point of the script: nested insertion sites sit inside
Alu/LINE-1/SVA elements by construction, so their flanking sequence IS element sequence.
Comparing nested against non-nested insertion sites would therefore compare element
composition against genomic composition, not target preference against target preference.
The correct control for a nested site is a RANDOM POSITION INSIDE A SAME-FAMILY ELEMENT.
Both controls are computed; the matched one is the one to read.

Motifs, all counted on both strands unless stated:
  CCATT / AATGG   - Jacob-Hirsch 2018 non-classical nested-insertion target motif
  TTTTAA / TTAAAA - classical L1 endonuclease consensus (TTTT/AA)
"""
import gzip
import json

import numpy as np
import pandas as pd
import py2bit
from scipy import stats

W = 30                      # breakpoint +/- W
N_CTRL_PER_CALL = 20        # matched control positions drawn per nested call
COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


def normalize_family(rep_name, rep_class, rep_family):
    t = f"{rep_name} {rep_class} {rep_family}".upper()
    if "ALU" in t:
        return "ALU"
    if "SVA" in t:
        return "SVA"
    if "LINE1" in t or "L1" in t:
        return "LINE1"
    return ""


def load_rmsk(path, chroms):
    rm = pd.read_csv(
        path, sep="\t", header=None, usecols=[5, 6, 7, 9, 10, 11, 12],
        names=["chrom", "start0", "end0", "strand", "repName", "repClass", "repFamily"],
        dtype={"chrom": str, "strand": str, "repName": str, "repClass": str, "repFamily": str},
        compression="gzip", low_memory=False,
    )
    rm = rm[rm.chrom.isin(chroms)].copy()
    rm["family"] = [normalize_family(a, b, c)
                    for a, b, c in zip(rm.repName, rm.repClass, rm.repFamily)]
    return rm[(rm.family != "") & (rm.end0 > rm.start0)].copy()


def seq_at(tb, clen, chrom, pos, w=W):
    s = tb.sequence(chrom, max(0, pos - w), min(clen[chrom], pos + w)).upper()
    return s if "N" not in s else s.replace("N", "")


def motif_counts(seqs):
    """Return per-sequence counts for each motif, both strands pooled."""
    out = {}
    for name, fwd, rev in [("CCATT", "CCATT", "AATGG"),
                           ("EN", "TTTTAA", "TTAAAA")]:
        out[name] = np.array([s.count(fwd) + s.count(rev) for s in seqs])
    out["bp"] = np.array([len(s) for s in seqs])
    return out


def rate_test(a_counts, a_bp, b_counts, b_bp, label_a, label_b):
    """Poisson rate-ratio test between two sets of motif counts over unequal total bp."""
    ka, kb = int(a_counts.sum()), int(b_counts.sum())
    na, nb = int(a_bp.sum()), int(b_bp.sum())
    ra, rbv = ka / na, kb / nb
    # conditional binomial test: given ka+kb events, expected share of a is na/(na+nb)
    p = stats.binomtest(ka, ka + kb, na / (na + nb)).pvalue if (ka + kb) else np.nan
    return dict(set_a=label_a, set_b=label_b, k_a=ka, bp_a=na, per_kb_a=1e3 * ra,
                k_b=kb, bp_b=nb, per_kb_b=1e3 * rbv,
                rate_ratio=(ra / rbv) if rbv else np.nan, p=p)


def main(per_call_csv, rmsk_path, twobit_path, out_prefix="notes/ccatt"):
    per = pd.read_csv(per_call_csv)
    tb = py2bit.open(twobit_path)
    allc = tb.chroms()
    chroms = sorted(set(per.chrom.unique()))
    clen = {c: allc[c] for c in chroms}
    rm = load_rmsk(rmsk_path, set(chroms))

    nested = per[per.nested].reset_index(drop=True)
    unnested = per[~per.nested].reset_index(drop=True)

    seq_nested = [seq_at(tb, clen, r.chrom, int(r.pos)) for r in nested.itertuples()]
    seq_unnested = [seq_at(tb, clen, r.chrom, int(r.pos)) for r in unnested.itertuples()]

    # matched control: random positions inside same-family elements, same family mix
    rng = np.random.default_rng(17)
    grp = {k: g[["chrom", "start0", "end0"]].to_numpy()
           for k, g in rm.groupby("family", sort=False)}
    seq_ctrl, ctrl_fam = [], []
    for fam, n in nested.MEIFAMILY.value_counts().items():
        iv = grp.get(fam)
        if iv is None:
            continue
        pick = rng.integers(0, len(iv), n * N_CTRL_PER_CALL)
        for j in pick:
            c, s0, e0 = iv[j]
            if int(e0) - int(s0) < 5:
                continue
            p = int(rng.integers(int(s0), int(e0)))
            seq_ctrl.append(seq_at(tb, clen, c, p))
            ctrl_fam.append(fam)

    mn, mu, mc = motif_counts(seq_nested), motif_counts(seq_unnested), motif_counts(seq_ctrl)
    rows = []
    for motif in ("CCATT", "EN"):
        rows.append({**rate_test(mn[motif], mn["bp"], mc[motif], mc["bp"],
                                 "nested insertion sites",
                                 "random positions inside same-family elements"),
                     "motif": motif, "comparison": "MATCHED (read this one)"})
        rows.append({**rate_test(mn[motif], mn["bp"], mu[motif], mu["bp"],
                                 "nested insertion sites", "non-nested insertion sites"),
                     "motif": motif, "comparison": "unmatched (composition-confounded)"})
    res = pd.DataFrame(rows)[["motif", "comparison", "set_a", "set_b", "k_a", "per_kb_a",
                              "k_b", "per_kb_b", "rate_ratio", "p"]]
    res.to_csv(f"{out_prefix}_enrichment.csv", index=False)

    # orientation-resolved: is CCATT oriented relative to the INSERTION's own strand?
    orient_rows = []
    for cat, mask in [("nested_sense", nested.same_orientation),
                      ("nested_antisense", ~nested.same_orientation)]:
        sub = nested[mask]
        if not len(sub):
            continue
        fwd = rev = bp = 0
        for r in sub.itertuples():
            s = seq_at(tb, clen, r.chrom, int(r.pos))
            if r.ORIENT == "-":
                s = rc(s)                     # read in the INSERTION's orientation
            fwd += s.count("CCATT")
            rev += s.count("AATGG")
            bp += len(s)
        orient_rows.append(dict(category=cat, n=len(sub), bp=bp,
                                CCATT_insertion_sense=fwd, CCATT_insertion_antisense=rev,
                                per_kb_sense=1e3 * fwd / bp, per_kb_antisense=1e3 * rev / bp,
                                strand_ratio=fwd / rev if rev else np.nan))
    ori = pd.DataFrame(orient_rows)
    ori.to_csv(f"{out_prefix}_orientation.csv", index=False)

    json.dump({"W": W, "n_nested": len(nested), "n_unnested": len(unnested),
               "n_control": len(seq_ctrl), "control_per_call": N_CTRL_PER_CALL},
              open(f"{out_prefix}_meta.json", "w"), indent=1)
    return res, ori


if __name__ == "__main__":
    import sys
    res, ori = main(sys.argv[1], "data/rmsk.txt.gz", "data/hg38.2bit")
    print(res.to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    print()
    print(ori.to_string(index=False, float_format=lambda v: f"{v:.4g}"))
