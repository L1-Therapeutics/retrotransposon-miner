#!/bin/bash
# Job 9 — exon-proximity replication and the EXON_DIST/EXON_SIDE annotator.
# Run from anywhere:  bash ~/Desktop/Research/Research_Projects/L1/_job9/run_job9.sh
set -u

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
L="$(cd "$HERE/.." && pwd)"
PY="$L/_verify/.venv-verify/bin/python"
OUT="$L/_job9/out"

fail() { echo; echo "STOPPED: $1"; exit 1; }

echo "=============================================="
echo " Job 9  —  splicing / exon proximity"
echo " L1 folder: $L"
echo "=============================================="
echo

# ---- 0. environment ---------------------------------------------------
[ -x "$PY" ] || fail "no interpreter at $PY
  The venv lives at L1/_verify/.venv-verify. If you moved or rebuilt it,
  tell me and I will re-point this script."
echo "[0/4] interpreter"
"$PY" -c "import sys,pandas,numpy,scipy;print('      python',sys.version.split()[0],'| pandas',pandas.__version__,'| numpy',numpy.__version__,'| scipy',scipy.__version__)" \
  || fail "the interpreter cannot import pandas/numpy/scipy"

# ---- 1. inputs --------------------------------------------------------
RMSK="$L/nested_analysis/data/rmsk.txt.gz"
GTF="$L/gencode/gencode.v49.annotation.gtf.gz"
VCF="$L/vcf/HG03086_classifier_ge_0.997.genes.snpeff.vcf"
for f in "$RMSK" "$GTF" "$VCF"; do
  [ -f "$f" ] || fail "missing input: $f"
done
echo "[1/4] inputs present"
echo "      rmsk    $(du -h "$RMSK" | cut -f1)"
echo "      gencode $(du -h "$GTF"  | cut -f1)"

# ---- 2. checksums -----------------------------------------------------
echo "[2/4] checksums"
( cd "$HERE" && shasum -a 256 -c CHECKSUMS.txt ) || fail "a script does not match CHECKSUMS.txt.
  You may have an older copy. Re-copy the four .py files from _job9 and retry."

# ---- 3. tests ---------------------------------------------------------
echo "[3/4] tests  (expect: 38 passed)"
( cd "$HERE" && "$PY" -m pytest test_alu_exon_proximity.py test_exon_distance_annotate.py -q ) \
  || fail "tests failed — do not trust any numbers below. Send me the output."

# ---- 4. analyses ------------------------------------------------------
mkdir -p "$OUT"
echo
echo "[4/4] analysis A — AluY exon proximity  (~35 s)"
"$PY" "$HERE/alu_exon_proximity.py" "$RMSK" "$GTF" --out-dir "$OUT/aluY" || fail "analysis A failed"

echo
echo "[4/4] analysis B — annotate HG03086 with EXON_DIST / EXON_SIDE  (~35 s)"
"$PY" "$HERE/exon_distance_annotate.py" "$VCF" "$GTF" "$OUT/HG03086.exondist.vcf" || fail "analysis B failed"

# ---- verification -----------------------------------------------------
echo
echo "=============================================="
echo " VERIFICATION — these must match exactly"
echo "=============================================="
"$PY" - "$OUT" <<'EOF'
import json, sys, pathlib
out = pathlib.Path(sys.argv[1])
s = json.load(open(out / "aluY" / "summary.json"))
exp = {"copies_total": 114541, "copies_in_single_strand_genes": 67379, "overall_sense_fraction": 0.4501}
ok = True
for k, v in exp.items():
    got = s[k]
    flag = "OK " if got == v else "BAD"
    if got != v:
        ok = False
    print(f"  {flag} {k:32s} {got}   (expected {v})")

import csv
rows = list(csv.DictReader(open(out / "aluY" / "tests.tsv"), delimiter="\t"))
up = next(r for r in rows if r["side"] == "upstream")
flag = "OK " if up["n_near"] == "1508" else "BAD"
if flag == "BAD":
    ok = False
print(f"  {flag} {'upstream n_near':32s} {up['n_near']}   (expected 1508)")
for name, key, want in [("upstream antisense_frac", "antisense_frac_near", 0.482095),
                        ("upstream odds_ratio", "odds_ratio", 0.753271)]:
    got = float(up[key])
    flag = "OK " if abs(got - want) < 1e-6 else "BAD"
    if flag == "BAD":
        ok = False
    print(f"  {flag} {name:32s} {got:.6f}   (expected {want:.6f})")
p = float(up["p"])
flag = "OK " if 5.5e-8 < p < 5.7e-8 else "BAD"
if flag == "BAD":
    ok = False
print(f"  {flag} {'upstream p':32s} {p:.6e}   (expected 5.614667e-08)")

side = {}
for line in open(out / "HG03086.exondist.vcf"):
    if line.startswith("#"):
        continue
    for f in line.split("\t")[7].split(";"):
        if f.startswith("EXON_SIDE="):
            side[f[10:]] = side.get(f[10:], 0) + 1
exp_side = {".": 784, "upstream": 537, "downstream": 531}
for k, v in exp_side.items():
    got = side.get(k, 0)
    flag = "OK " if got == v else "BAD"
    if got != v:
        ok = False
    print(f"  {flag} {'EXON_SIDE=' + k:32s} {got}   (expected {v})")
print()
print("  ALL CHECKS PASSED" if ok else "  SOMETHING DID NOT MATCH — send me this block")
EOF

echo
echo "Outputs are in $OUT"
echo "Done."
