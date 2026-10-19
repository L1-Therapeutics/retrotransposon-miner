#!/bin/bash
# Job 10 — the poly(A) subfamily replication (the gate on sending that result to
# William) and the LINE-1 downstream probe with its family controls.
# Run:  bash ~/Desktop/Research/Research_Projects/L1/_job10/run_job10.sh
set -u

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
L="$(cd "$HERE/.." && pwd)"
PY="$L/_verify/.venv-verify/bin/python"
OUT="$HERE/out"

fail() { echo; echo "STOPPED: $1"; exit 1; }

echo "=============================================="
echo " Job 10  —  subfamily replication + LINE-1 probe"
echo " Takes about 5.5 minutes. Safe to leave running."
echo "=============================================="
echo

[ -x "$PY" ] || fail "no interpreter at $PY"
echo "[0/5] interpreter"
"$PY" -c "import sys,pandas,numpy,scipy,twobitreader;print('      python',sys.version.split()[0],'| pandas',pandas.__version__,'| twobitreader ok')" \
  || fail "the interpreter is missing a package (pandas/numpy/scipy/twobitreader)"

RMSK="$L/nested_analysis/data/rmsk.txt.gz"
GTF="$L/gencode/gencode.v49.annotation.gtf.gz"
TWOBIT="$L/ref/hg38.2bit"
for f in "$RMSK" "$GTF" "$TWOBIT"; do [ -f "$f" ] || fail "missing input: $f"; done
echo "[1/5] inputs present"

echo "[2/5] checksums"
( cd "$HERE" && shasum -a 256 -c CHECKSUMS.txt ) || fail "a script does not match CHECKSUMS.txt"

echo "[3/5] tests  (expect: 50 passed)"
( cd "$HERE" && "$PY" -m pytest test_alu_polya_asymmetry.py test_alu_exon_proximity.py \
    test_line1_downstream_probe.py -q ) || fail "tests failed — send me the output, do not use any numbers"

mkdir -p "$OUT"
echo
echo "[4/5] poly(A) carrier replication — 3 subfamilies, about 3 minutes"
echo "      NOTE the --max-milli-div values: the default 60 keeps 31% of AluY but"
echo "      under 1% of AluSx/AluJo, so each subfamily gets its own 31st percentile."
"$PY" "$HERE/alu_polya_asymmetry.py" "$RMSK" "$TWOBIT" --gencode "$GTF" \
      --rep-name AluY  --max-milli-div 60  --tail-start 270 --out-dir "$OUT/AluY"  > "$OUT/AluY.log"  2>&1 || fail "AluY run failed"
echo "      AluY done"
"$PY" "$HERE/alu_polya_asymmetry.py" "$RMSK" "$TWOBIT" --gencode "$GTF" \
      --rep-name AluSx --max-milli-div 105 --out-dir "$OUT/AluSx" > "$OUT/AluSx.log" 2>&1 || fail "AluSx run failed"
echo "      AluSx done"
"$PY" "$HERE/alu_polya_asymmetry.py" "$RMSK" "$TWOBIT" --gencode "$GTF" \
      --rep-name AluJo --max-milli-div 151 --out-dir "$OUT/AluJo" > "$OUT/AluJo.log" 2>&1 || fail "AluJo run failed"
echo "      AluJo done"

echo
echo "[5/5] LINE-1 downstream probe and family controls, about 2.5 minutes"
"$PY" "$HERE/line1_downstream_probe.py" "$RMSK" "$GTF" --out-dir "$OUT/L1PA3_7" > "$OUT/L1PA3_7.log" 2>&1 || fail "LINE-1 probe failed"
echo "      L1PA3-7 done"
for N in AluSx MIR L2; do
  "$PY" "$HERE/line1_downstream_probe.py" "$RMSK" "$GTF" --rep-name "$N" --out-dir "$OUT/ctrl_$N" > "$OUT/ctrl_$N.log" 2>&1 \
    || fail "control $N failed"
  echo "      control $N done"
done

echo
echo "=============================================="
echo " VERIFICATION — every line must say OK"
echo "=============================================="
"$PY" - "$OUT" <<'EOF'
import csv, json, pathlib, sys
out = pathlib.Path(sys.argv[1]); ok = True
def check(name, got, want, tol=None):
    global ok
    good = (abs(got - want) <= tol) if tol is not None else (got == want)
    if not good: ok = False
    print(f"  {'OK ' if good else 'BAD'} {name:38s} {got}   (expected {want})")

exp = {"AluY":  (67379, 0.4501, 11923, 0.9496),
       "AluSx": (75839, 0.4397, 17795, 0.9379),
       "AluJo": (51437, 0.4646, 18598, 0.9988)}
for rep, (n_genes, sense, n_carr, or_mh) in exp.items():
    j = json.load(open(out / rep / "alu_polya_asymmetry.json"))
    pc = j["per_copy"]; e = pc["exposures"]["any_position"]
    canon = [t for t in j["tests"] if t["label"] == "full_consensus_canonical_only"][0]
    check(f"{rep} copies in genes", pc["n_copies_in_single_strand_genes"], n_genes)
    check(f"{rep} sense fraction", pc["overall_sense_fraction"], sense, 1e-9)
    check(f"{rep} carriers", e["n_carriers"], n_carr)
    check(f"{rep} adjusted OR", e["adjusted"]["or_mh"], or_mh, 1e-4)
    check(f"{rep} canonical hexamers on consensus",
          canon["observed_sense"] + canon["observed_antisense"], 0)

rows = list(csv.DictReader(open(out / "L1PA3_7" / "overall.tsv"), delimiter="\t"))
d = next(r for r in rows if r["side"] == "downstream")
check("L1PA3-7 downstream n", int(d["n_near"]), 309)
check("L1PA3-7 downstream OR", float(d["OR"]), 1.561, 1e-3)
u = next(r for r in rows if r["side"] == "upstream")
check("L1PA3-7 upstream OR", float(u["OR"]), 1.010, 1e-3)

for fam, want_up in [("AluSx", 0.809), ("MIR", 0.889), ("L2", 1.563)]:
    rows = list(csv.DictReader(open(out / f"ctrl_{fam}" / "overall.tsv"), delimiter="\t"))
    up = next(r for r in rows if r["side"] == "upstream")
    check(f"{fam} upstream OR", float(up["OR"]), want_up, 1e-3)
print()
print("  ALL CHECKS PASSED" if ok else "  SOMETHING DID NOT MATCH — send me this block")
EOF
echo
echo "Outputs are in $OUT"
