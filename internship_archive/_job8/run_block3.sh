#!/bin/bash
# Reproduce the Alu poly(A) analysis. Run from anywhere:
#     bash ~/Desktop/Research/Research_Projects/L1/_job8/run_block3.sh
# Does not need an activated venv - it calls the venv's interpreter directly.
set -euo pipefail

L="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PY="$L/_verify/.venv-verify/bin/python"
JOB="$L/_job8"

echo "== 0. environment =="
test -x "$PY" || { echo "FAIL: no interpreter at $PY"; exit 1; }
"$PY" -c "import sys,pandas,numpy,scipy,twobitreader; print('python',sys.version.split()[0],'| pandas',pandas.__version__)"
for f in "$L/nested_analysis/data/rmsk.txt.gz" "$L/ref/hg38.2bit" "$L/gencode/gencode.v49.annotation.gtf.gz"; do
  test -f "$f" || { echo "FAIL: missing $f"; exit 1; }
done
echo "inputs present"

echo
echo "== 1. file integrity =="
cd "$JOB"
shasum -a 256 -c CHECKSUMS.txt

echo
echo "== 2. the scanner this script imports =="
cp -f "$L/artifacts/l1_strand_scan.py" "$JOB/l1_strand_scan.py"
echo "copied l1_strand_scan.py from artifacts/"

echo
echo "== 3. tests =="
PYTHONPATH=. "$PY" -m pytest test_alu_polya_asymmetry.py -q

echo
echo "== 4. analysis (about 45 s) =="
PYTHONPATH=. "$PY" alu_polya_asymmetry.py \
  "$L/nested_analysis/data/rmsk.txt.gz" \
  "$L/ref/hg38.2bit" \
  --gencode "$L/gencode/gencode.v49.annotation.gtf.gz" \
  --replicates 2000 --out-dir out

echo
echo "== 5. verification =="
cat out/polya_hits.tsv
PYTHONPATH=. "$PY" - <<'EOF'
import json
import pandas as pd
c = json.load(open("out/alu_polya_asymmetry.json"))["consensus"]
print(f"consensus length {len(c)} | AATAAA {c.count('AATAAA')} | ATTAAA {c.count('ATTAAA')}")
print(f"tail {c[-32:]}")
d = pd.read_csv("out/alu_per_copy.csv.gz")
print(f"copies {len(d)} | sense fraction {d.orient.eq('sense').mean():.4f}")
print(d.groupby(d.n_hex_any > 0).orient.value_counts().unstack().to_string())
EOF
echo
echo "DONE. Every number above should match RUNBOOK_today.md Block 3."
