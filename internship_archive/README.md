# Internship working code (July - October 2026)

Verbatim copies of the analysis code written during the L1 Therapeutics internship, kept
exactly as run. NOT FOR MERGE AS-IS: this is an archive, so the same script appears in
more than one place where a delivered package carried its own copy.

Layout mirrors the working folder:

| Directory | What it is |
|---|---|
| `_job8/`, `_job9/`, `_job10/` | Self-contained delivered packages: scripts, tests, a `run_*.sh` runner, a runbook and `CHECKSUMS.txt`. `_job10` holds the latest versions of the exon-proximity, poly(A)-carrier and LINE-1 downstream scripts. |
| `_gene_map/` | The same analyses, one folder each (splice, exon_prox, exon_dist, polya_subfamily, carrier_model, l1_downstream), plus `verify_aluY_table.py`. |
| `artifacts/`, `artifacts_2/` | Earlier nested-insertion, TSD, L1-length and report-building scripts. `artifacts_2` is an older snapshot. |
| `marker_map/` | Scripts that render the evidence map. |
| `report/` | Figure scripts for the report. |

Notes
- The packages are frozen deliveries. Each `CHECKSUMS.txt` also lists input and output
  files that are not in the repository, so `shasum -c` reports those as missing.
- A few scripts contain a home-directory path as a default. Override it on the command line.
- Related branches: `proposal/splice-exonization-annotation` (a design change that needs a
  decision) and `analysis/l1-transduction-readthrough` (three more standalone scripts).
