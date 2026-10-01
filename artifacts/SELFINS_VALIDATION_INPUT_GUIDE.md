# HG03086 nested-MEI validation: exact input and action guide

**Status:** the saved analyses and per-call QC are complete. The Claude export confirms that `hg38.2bit` (835,393,456 bytes) and `chr22_seq.json` (50,818,620 bytes) existed as artifacts in the earlier session, but the export stores their metadata and conversation references—not their binary/file contents. The current checkout and `artifacts.zip` do not contain those payloads. The HG03086 BAM/CRAM and its index are also absent. Therefore the event-level TPRT triad, de novo junction TSD reconstruction, and read-backed validation cannot be finished from the files currently on disk.

## Smallest useful handoff

Please provide **either** (A) access to the original local analysis workspace/data directory, **or** (B) the following files in a directory you can share/read locally:

1. **Reference:** GRCh38/hg38 reference sequence matching the callset, preferably `hg38.2bit` (the earlier Claude artifact reports 835,393,456 bytes), or a primary-assembly FASTA plus `.fai`. The exact source/release and sequence checksum should be recorded. `chr22_seq.json` alone only covers chr22 and cannot support the genome-wide analysis.
2. **Read evidence:** the HG03086 BAM plus `.bai`, or CRAM plus `.crai` **and the exact reference used to decode it**. Confirm that reads come from the same sample/build as the VCF. If controlled access prevents sharing the alignment, provide a secure local environment where the existing QC script can run, or export per-locus evidence (details below).
3. **Candidate calls / insertion metadata:** use the existing [HG03086 VCF](../HG03086_classifier_ge_0.997.genes.snpeff%20%281%29.vcf) and [per-call QC table](selfins_tsd_call_qc.csv). Confirm their coordinate convention (VCF POS is 1-based; reference intervals used by scripts should be 0-based half-open) and provide any caller output that stores both breakpoint/junction coordinates, assembled inserted sequence, or candidate left/right target flanks. The current VCF TSD string and `POLYA_MIN_BP` are not substitutes for those data.
4. **Callability/opportunity background:** provide the exact pipeline's genome-wide callable/mappable MEI-target intervals or enough inputs/configuration to build them (read length, aligner and version, mapping-quality threshold, depth/coverage constraints, repeat-aware masks, and any candidate-detection filters). For the SVA null, we especially need the fraction of reference SVA sequence that the same pipeline could detect under comparable insertion conditions. A generic mappability track alone is not an MEI-calling detection model.
5. **Versioning:** attach the HG03086 pipeline/caller version and command/configuration used to create the classifier VCF, plus checksums for each supplied file. A copy of the VCF already exists here, so only resend if its provenance/version differs.

## If full BAM/CRAM cannot be shared

Run locally in the controlled environment and return one record per call with:

- call ID and chr/1-based POS;
- left and right breakpoint coordinates and the exact reference-oriented flanking sequences (at least 100 bp each, longer if available);
- split-read and discordant-pair names/counts or sanitized evidence summaries supporting **each** junction;
- local assembled alternate haplotype spanning both junctions and the inserted sequence, if the caller/assembler can produce it;
- coverage, mapping qualities, strand/orientation, and whether the site is uniquely/ambiguously mappable;
- the callable/opportunity status for this genomic interval and the matching reference-repeat annotation.

Do not include sample identifiers or raw reads in an unsecured handoff. Coordinates and limited sequence flanks can themselves be sensitive genomic data; use William's approved secure channel. If sequence cannot leave that environment, run the validation there and return aggregated statistics plus per-call pass/fail reasons and checksums.

## Analyses to run once inputs are available

1. **Rebuild event structure:** reconcile VCF calls against read-backed left/right breakpoints and assembled alleles. Classify resolved, ambiguous, one-sided, and unsupported calls; retain unresolved calls in the denominator transparently.
2. **Validate TSDs:** independently infer duplicated target sequence at both insertion junctions, with orientation and microhomology accounted for. Compare reconstructed sequence and length with the VCF's `TSD` field; flag low-complexity matches and boundary ambiguity. Report separate distributions by family and call-support tier. Do not count a plausible length alone as a TSD.
3. **Validate the poly(A) tail:** read/assembly-backed inserted 3′ sequence, in insertion orientation; report tail length and sequence-composition rule. Treat `POLYA_MIN_BP` only as classifier evidence.
4. **Score endonuclease motifs:** score strand-compatible L1 EN motif opportunity at the reconstructed nick/breakpoint, with sequence composition and motif opportunity in callable matched controls. Report the exact motif and matching logic, not a generic “EN motif present” label.
5. **Tail-peak cluster:** apply the same validation blind to the 280–300 bp bin and to prespecified comparison bins, then compare support fractions and TPRT hallmark completeness. Do not select only calls that already have favorable TSD/poly(A) values.
6. **SVA repeat-aware null:** rerun the same-family nesting null using detectable/callable reference SVA opportunities, with the calling pipeline's real filters. Show the uncorrected and corrected expected counts side by side; with nine observed events, include uncertainty and avoid treating one point estimate as settled.
7. **Orientation/mechanism:** only interpret the CCATT vs classical-route contrast after the event-level motif and truncation evidence is scored; the existing CCATT null test is “no enrichment detected under this control,” not a route assignment.

## Current deliverables to inspect first

- [Integrated evidence report](NESTED_MEI_INTEGRATED_REPORT_2026-10-02.md)
- [Per-call TSD/poly(A)/support QC](selfins_tsd_call_qc.csv)
- [Audit summary, including all nine nested SVA calls](selfins_tsd_audit_summary.json)
- [Reproducible audit script](audit_nested_tsd.py)
- [Earlier full analysis report](selfins_report.md)

**Important distinction:** the Claude export's artifact manifest lists the 2bit file, but that does not mean the actual binary is in this checkout. Its transcript confirms successful downloads into a different Claude workspace, while current filesystem searches find neither `hg38.2bit` nor the 50 MB `chr22_seq.json`, and the archived `artifacts.zip` omits them. The exact next step is to restore/copy those files from that original workspace if still available, or provide the reference and controlled-access read data as above. No genome-wide read-backed validation has been silently inferred from the manifest.