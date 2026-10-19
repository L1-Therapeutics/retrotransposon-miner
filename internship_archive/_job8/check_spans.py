"""Count LINE-1 near-full-length calls and how many reach the consensus ends.

Run from the L1 folder:  python _job8/check_spans.py
(This is a Python file on purpose. Pasting Python into zsh gives
 'parse error near &' because zsh reads & as a background operator.)
"""
import glob
import sys

rows = []
files = sorted(glob.glob("vcf/*_classifier_ge_0.997.genes.snpeff.vcf"))
if not files:
    sys.exit("no VCFs found - run this from the L1 folder, not from _job8")
for path in files:
    with open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            info = dict(kv.split("=", 1)
                        for kv in line.split("\t")[7].split(";") if "=" in kv)
            rows.append((info.get("MEIFAMILY"),
                         int(float(info.get("MEI_SPAN") or 0)),
                         int(float(info.get("MEI_5P") or -1)),
                         int(float(info.get("MEI_3P") or -1))))

line1 = [r for r in rows if r[0] == "LINE1" and r[1] > 0 and r[2] >= 0]
big = [r for r in line1 if r[1] >= 5500]
print(f"VCF files                     {len(files)}")
print(f"all calls                     {len(rows)}")
print(f"LINE-1 non-sentinel           {len(line1)}")
print(f"  of those, span >= 5500      {len(big)}")
print(f"  of those, MEI_5P <= 100     {sum(1 for r in big if r[2] <= 100)}")
print(f"  of those, MEI_3P >= 6000    {sum(1 for r in big if r[3] >= 6000)}")
