"""Compare your reproduction of the gene-annotation analysis with the reference run.

Usage (from the L1 folder):
    python check_gene_numbers.py gene_run1/gene_annotation_report.json gene_run2/gene_gtf_report.json

Stdlib only. Prints PASS/FAIL per number and exits non-zero if any FAIL.
Expected values come from the reference run (GENCODE v49, 10 genomes, merge 50 bp).
"""
import json
import sys

r1 = json.load(open(sys.argv[1]))
r2 = json.load(open(sys.argv[2]))
a, l1 = r1["annotation"], r1["line1_strand_asymmetry"]
sr, ora = r2["strand_resolution"], r2["orientation_agreement_calls"]

checks = [
    ("unique sites", r1["unique_sites"], 6410, 0),
    ("calls where GENE and GENEID differ in length", a["gene_vs_geneid_length_mismatch"], 27, 0),
    ("calls where GENE_STRAND and GENEID differ in length", a["strand_vs_geneid_length_mismatch"], 0, 0),
    ("intronic calls (snpEff CSQ)", a["intronic_calls"], 11008, 0),
    ("  ... with an up/downstream term", a["intronic_with_updown_term"], 3679, 0),
    ("  ... mixed-strand gene list", a["intronic_mixed_strand"], 2591, 0),
    ("  ... mixed-strand AND up/downstream term", a["intronic_mixed_strand_with_updown_term"], 1817, 0),
    ("LINE-1 sites analysed (snpEff set)", r1["analysed_sites_by_family"]["LINE1"], 452, 0),
    ("LINE-1 sense in + genes", l1["plus_sense"], 76, 0),
    ("LINE-1 + genes total", l1["plus_n"], 239, 0),
    ("LINE-1 sense in - genes", l1["minus_sense"], 102, 0),
    ("LINE-1 - genes total", l1["minus_n"], 213, 0),
    ("gene IDs found in GENCODE v49", r2["gene_ids_found_in_gtf"], 5427, 0),
    ("intronic sites, GENCODE overlap", sr["intronic_sites"], 4240, 0),
    ("single-strand in both, strand agrees", sr["strand_agrees"], 3232, 0),
    ("snpEff mixed but GENCODE single-strand", sr["snpeff_mixed_but_gtf_single"], 496, 0),
    ("snpEff single but GENCODE mixed", sr["snpeff_single_but_gtf_mixed"], 0, 0),
    ("sites analysed on GENCODE set", r2["analysed_sites"]["n"], 3461, 0),
    ("ORIENT vs poly(A) agreement, Alu", ora["ALU"]["agreement"], 0.9099, 5e-5),
    ("ORIENT vs poly(A) agreement, LINE-1", ora["LINE1"]["agreement"], 0.7052, 5e-5),
    ("ORIENT vs poly(A) agreement, SVA", ora["SVA"]["agreement"], 0.7641, 5e-5),
    ("LINE-1 sense fraction (ORIENT)", r2["sense_by_ORIENT"]["LINE1"]["sense_frac"], 0.398, 5e-4),
    ("LINE-1 asymmetry p, ORIENT", r2["sense_by_ORIENT"]["LINE1"]["asymmetry_fisher_p"], 0.0010, 5e-5),
    ("LINE-1 asymmetry p, poly(A) orientation", r2["sense_by_polyA"]["LINE1"]["asymmetry_fisher_p"], 0.1371, 5e-5),
    ("LINE-1 asymmetry p, measures agree", r2["sense_by_ORIENT_where_measures_agree"]["LINE1"]["asymmetry_fisher_p"], 0.0170, 5e-4),
]
bad = 0
for label, got, want, tol in checks:
    ok = abs(got - want) <= tol
    bad += not ok
    print(f"{'PASS' if ok else 'FAIL'}  {label:55s} got {got:<10.5g} want {want:<10.5g}")
print(f"\n{len(checks) - bad}/{len(checks)} match")
sys.exit(1 if bad else 0)
