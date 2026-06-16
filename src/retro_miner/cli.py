from __future__ import annotations

from pathlib import Path

import click

from retro_miner.candidate_loci import build_candidate_loci
from retro_miner.evidence_extract import (
    ExtractionSummary,
    extract_discordant_evidence,
    extract_split_evidence,
)
from retro_miner.mei_support import annotate_candidate_loci_with_mei


@click.group()
def cli() -> None:
    """Retrotransposon miner command-line interface."""


@cli.command("check-env")
def check_env() -> None:
    """Print basic environment status."""
    click.echo("retrotransposon-miner CLI is installed.")
    click.echo("Run scripts/validate_environment.sh for full validation.")


@cli.command("extract-split-evidence")
@click.option("--tumor-bam", type=click.Path(exists=True, dir_okay=False, path_type=Path), required=True)
@click.option("--normal-bam", type=click.Path(exists=True, dir_okay=False, path_type=Path), required=True)
@click.option("--outdir", type=click.Path(file_okay=False, path_type=Path), required=True)
@click.option("--region", default="chr22", show_default=True, help="Region to scan, e.g. chr22 or chr1:1-1000000.")
@click.option(
    "--regions",
    default=None,
    help="Optional comma-separated region/chromosome list (overrides --region), e.g. chr15,chr16,chr17.",
)
@click.option("--min-mapq", default=20, show_default=True, type=int, help="Minimum split-read anchor mapping quality.")
@click.option(
    "--min-mapq-discordant",
    default=0,
    show_default=True,
    type=int,
    help="Minimum discordant-read anchor mapping quality.",
)
@click.option("--min-clip-len", default=20, show_default=True, type=int, help="Minimum soft-clip length.")
@click.option(
    "--poly-tail-rescue-min-clip-len",
    default=8,
    show_default=True,
    type=int,
    help="Minimum soft-clip length to rescue reads with strong polyA/polyT tails at breakpoint.",
)
@click.option(
    "--poly-tail-rescue-min-run",
    default=8,
    show_default=True,
    type=int,
    help="Minimum longest homopolymer run of A/T in clip sequence for poly-tail rescue.",
)
@click.option(
    "--poly-tail-rescue-min-frac",
    default=0.8,
    show_default=True,
    type=float,
    help="Minimum A/T fraction in clip sequence for poly-tail rescue.",
)
@click.option(
    "--with-discordant/--no-discordant",
    default=True,
    show_default=True,
    help="Also extract discordant paired-end evidence in the same run.",
)
@click.option(
    "--discordant-quantile",
    default=0.995,
    show_default=True,
    type=float,
    help="Quantile used to define large-insert discordant thresholds.",
)
@click.option(
    "--discordant-min-abs-tlen",
    default=1000,
    show_default=True,
    type=int,
    help="Minimum absolute template length allowed for large-insert discordant evidence.",
)
@click.option(
    "--discordant-poly-tail-rescue-window-bases",
    default=25,
    show_default=True,
    type=int,
    help="Number of terminal bases to scan on each read end for discordant-anchor polyA/polyT rescue.",
)
@click.option(
    "--discordant-poly-tail-rescue-min-run",
    default=10,
    show_default=True,
    type=int,
    help="Minimum A/T homopolymer run in breakpoint-proximal discordant-anchor sequence.",
)
@click.option(
    "--discordant-poly-tail-rescue-min-frac",
    default=0.8,
    show_default=True,
    type=float,
    help="Minimum A/T fraction in breakpoint-proximal discordant-anchor sequence.",
)
@click.option(
    "--discordant-poly-tail-rescue-min-abs-tlen",
    default=500,
    show_default=True,
    type=int,
    help="Minimum |TLEN| structural context for discordant-anchor polyA/polyT rescue.",
)
def extract_split_evidence_cmd(
    tumor_bam: Path,
    normal_bam: Path,
    outdir: Path,
    region: str,
    regions: str | None,
    min_mapq: int,
    min_mapq_discordant: int,
    min_clip_len: int,
    poly_tail_rescue_min_clip_len: int,
    poly_tail_rescue_min_run: int,
    poly_tail_rescue_min_frac: float,
    with_discordant: bool,
    discordant_quantile: float,
    discordant_min_abs_tlen: int,
    discordant_poly_tail_rescue_window_bases: int,
    discordant_poly_tail_rescue_min_run: int,
    discordant_poly_tail_rescue_min_frac: float,
    discordant_poly_tail_rescue_min_abs_tlen: int,
) -> None:
    """Extract split-read MEI evidence and optional discordant evidence."""
    outdir.mkdir(parents=True, exist_ok=True)
    region_list = [r.strip() for r in regions.split(",")] if regions else [region]
    region_list = [r for r in region_list if r]
    if not region_list:
        raise click.ClickException("No valid regions provided via --region/--regions.")

    split_summaries: list[ExtractionSummary] = []
    discordant_summaries: dict[str, ExtractionSummary] = {}
    for sample, bam in (("tumor", tumor_bam), ("normal", normal_bam)):
        click.echo(f"[extract] sample={sample} bam={bam} regions={','.join(region_list)}")
        split_summary = extract_split_evidence(
            bam_path=bam,
            sample_name=sample,
            outdir=outdir,
            regions=region_list,
            min_mapq=min_mapq,
            min_clip_len=min_clip_len,
            poly_tail_rescue_min_clip_len=poly_tail_rescue_min_clip_len,
            poly_tail_rescue_min_run=poly_tail_rescue_min_run,
            poly_tail_rescue_min_frac=poly_tail_rescue_min_frac,
        )
        split_summaries.append(split_summary)
        click.echo(
            f"[done] sample={sample} scanned={split_summary.total_reads_scanned} "
            f"passing={split_summary.passing_reads} split_rows={split_summary.split_evidence_rows}"
        )
        if with_discordant:
            click.echo(f"[extract-discordant] sample={sample} regions={','.join(region_list)}")
            discordant_summary = extract_discordant_evidence(
                bam_path=bam,
                sample_name=sample,
                outdir=outdir,
                regions=region_list,
                min_mapq=min_mapq_discordant,
                insert_quantile=discordant_quantile,
                min_abs_tlen=discordant_min_abs_tlen,
                poly_tail_rescue_window_bases=discordant_poly_tail_rescue_window_bases,
                poly_tail_rescue_min_run=discordant_poly_tail_rescue_min_run,
                poly_tail_rescue_min_frac=discordant_poly_tail_rescue_min_frac,
                poly_tail_rescue_min_abs_tlen=discordant_poly_tail_rescue_min_abs_tlen,
            )
            discordant_summaries[sample] = discordant_summary
            click.echo(
                f"[done-discordant] sample={sample} scanned={discordant_summary.total_reads_scanned} "
                f"passing={discordant_summary.passing_reads} discordant_rows={discordant_summary.discordant_evidence_rows} "
                f"insert_threshold={discordant_summary.insert_size_threshold}"
            )

    summary_path = outdir / "split_evidence.summary.tsv"
    with summary_path.open("w", encoding="utf-8") as handle:
        handle.write(
            "sample\ttotal_reads_scanned\tpassing_reads\tsplit_evidence_rows\tdiscordant_evidence_rows\tinsert_size_threshold\n"
        )
        for s in split_summaries:
            d = discordant_summaries.get(s.sample)
            discordant_rows = d.discordant_evidence_rows if d else 0
            insert_threshold = d.insert_size_threshold if d else 0
            handle.write(
                f"{s.sample}\t{s.total_reads_scanned}\t{s.passing_reads}\t{s.split_evidence_rows}\t"
                f"{discordant_rows}\t{insert_threshold}\n"
            )
    click.echo(f"[summary] {summary_path}")


@cli.command("build-candidate-loci")
@click.option(
    "--evidence-dir",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    required=True,
    help="Directory containing split/discordant evidence outputs.",
)
@click.option(
    "--outdir",
    type=click.Path(file_okay=False, path_type=Path),
    default=None,
    help="Output directory for candidate loci table (defaults to evidence-dir).",
)
@click.option("--window-size", type=int, default=200, show_default=True, help="Window size in bp for locus binning.")
@click.option(
    "--split-cluster-bp",
    type=int,
    default=30,
    show_default=True,
    help="Breakpoint clustering distance for split-read evidence (bp).",
)
@click.option(
    "--discordant-cluster-bp",
    type=int,
    default=400,
    show_default=True,
    help="Breakpoint clustering distance for discordant evidence (bp).",
)
@click.option(
    "--max-locus-span-bp",
    type=int,
    default=2000,
    show_default=True,
    help="Maximum merged locus span to prevent over-merging chained discordant intervals.",
)
@click.option("--pseudocount", type=float, default=1.0, show_default=True, help="Pseudocount for enrichment ratio.")
@click.option(
    "--segdup-bed",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Optional BED for segmental duplications; overlaps are flagged in output.",
)
@click.option(
    "--segdup-min-fraction",
    type=float,
    default=0.1,
    show_default=True,
    help="Minimum fraction of candidate window required for segdup overlap flag.",
)
@click.option(
    "--mappability-bedgraph",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Optional BEDGRAPH for mappability; low-score overlaps are flagged in output.",
)
@click.option(
    "--mappability-low-threshold",
    type=float,
    default=0.5,
    show_default=True,
    help="Values below this mappability score are flagged as low mappability.",
)
@click.option(
    "--mappability-min-fraction",
    type=float,
    default=0.5,
    show_default=True,
    help="Minimum fraction of candidate window required for low-mappability overlap flag.",
)
@click.option(
    "--giab-highconf-bed",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Optional GIAB high-confidence BED; windows outside these intervals are flagged.",
)
@click.option(
    "--gap-bed",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Optional gap/mask BED (or UCSC table interpreted as BED columns); overlaps are flagged.",
)
@click.option(
    "--gap-min-fraction",
    type=float,
    default=0.1,
    show_default=True,
    help="Minimum fraction of candidate window required for gap overlap flag.",
)
@click.option(
    "--encode-blacklist-bed",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Optional ENCODE blacklist BED; overlaps are flagged.",
)
@click.option(
    "--encode-blacklist-min-fraction",
    type=float,
    default=0.1,
    show_default=True,
    help="Minimum fraction of candidate window required for ENCODE blacklist overlap flag.",
)
def build_candidate_loci_cmd(
    evidence_dir: Path,
    outdir: Path | None,
    window_size: int,
    split_cluster_bp: int,
    discordant_cluster_bp: int,
    max_locus_span_bp: int,
    pseudocount: float,
    segdup_bed: Path | None,
    segdup_min_fraction: float,
    mappability_bedgraph: Path | None,
    mappability_low_threshold: float,
    mappability_min_fraction: float,
    giab_highconf_bed: Path | None,
    gap_bed: Path | None,
    gap_min_fraction: float,
    encode_blacklist_bed: Path | None,
    encode_blacklist_min_fraction: float,
) -> None:
    """Aggregate split/discordant evidence into tumor-vs-normal candidate loci."""
    target_outdir = outdir if outdir is not None else evidence_dir
    tsv_path = build_candidate_loci(
        evidence_dir=evidence_dir,
        outdir=target_outdir,
        window_size=window_size,
        split_cluster_bp=split_cluster_bp,
        discordant_cluster_bp=discordant_cluster_bp,
        max_locus_span_bp=max_locus_span_bp,
        pseudocount=pseudocount,
        segdup_bed=segdup_bed,
        segdup_min_fraction=segdup_min_fraction,
        low_mappability_bedgraph=mappability_bedgraph,
        low_mappability_threshold=mappability_low_threshold,
        low_mappability_min_fraction=mappability_min_fraction,
        giab_highconf_bed=giab_highconf_bed,
        gap_bed=gap_bed,
        gap_min_fraction=gap_min_fraction,
        encode_blacklist_bed=encode_blacklist_bed,
        encode_blacklist_min_fraction=encode_blacklist_min_fraction,
    )
    click.echo(f"[candidate-loci] {tsv_path}")


@cli.command("annotate-mei-support")
@click.option(
    "--evidence-dir",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    required=True,
    help="Directory containing split evidence outputs.",
)
@click.option(
    "--candidate-loci",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=True,
    help="Candidate loci TSV from build-candidate-loci.",
)
@click.option(
    "--mei-fasta",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=True,
    help="MEI reference FASTA (family/subfamily headers retained).",
)
@click.option(
    "--out-tsv",
    type=click.Path(dir_okay=False, path_type=Path),
    required=True,
    help="Output TSV path for candidate loci with MEI support annotation.",
)
@click.option(
    "--reference-fasta",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Optional reference FASTA for TSD sequence extraction and breakpoint-context annotation.",
)
@click.option(
    "--tumor-bam-depth",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Optional tumor BAM path for local BAM-depth normalization on consistent/junk-clean loci.",
)
@click.option(
    "--normal-bam-depth",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Optional normal BAM path for local BAM-depth normalization on consistent/junk-clean loci.",
)
@click.option(
    "--rmsk-table",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Optional UCSC RepeatMasker rmsk table (e.g. rmsk.txt.gz) for nested-insertion annotation.",
)
@click.option(
    "--g1k-mei-vcf",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Optional 1000G/MELT MEI VCF(.gz) for overlap annotation.",
)
@click.option(
    "--g1k-split-padding-bp",
    type=int,
    default=200,
    show_default=True,
    help="Padding for split-resolved breakpoint-centered overlap queries (bp).",
)
@click.option(
    "--g1k-dpe-padding-min-bp",
    type=int,
    default=200,
    show_default=True,
    help="Minimum padding for discordant-dominant overlap queries (bp).",
)
@click.option(
    "--g1k-dpe-padding-max-bp",
    type=int,
    default=200,
    show_default=True,
    help="Maximum padding cap for discordant-dominant overlap queries (bp).",
)
@click.option(
    "--g1k-dpe-padding-tlen-factor",
    type=float,
    default=0.0,
    show_default=True,
    help="Scale factor applied to discordant TLEN mean when deriving discordant-dominant query padding.",
)
def annotate_mei_support_cmd(
    evidence_dir: Path,
    candidate_loci: Path,
    mei_fasta: Path,
    out_tsv: Path,
    reference_fasta: Path | None,
    tumor_bam_depth: Path | None,
    normal_bam_depth: Path | None,
    rmsk_table: Path | None,
    g1k_mei_vcf: Path | None,
    g1k_split_padding_bp: int,
    g1k_dpe_padding_min_bp: int,
    g1k_dpe_padding_max_bp: int,
    g1k_dpe_padding_tlen_factor: float,
) -> None:
    """Annotate candidate loci with MEI family/subfamily support and insertion span estimates."""
    click.echo("[mei-annotate] starting minimap2 clip-to-MEI alignment and locus annotation")
    if (tumor_bam_depth is None) ^ (normal_bam_depth is None):
        raise click.ClickException("Provide both --tumor-bam-depth and --normal-bam-depth, or neither.")
    out_path = annotate_candidate_loci_with_mei(
        evidence_dir=evidence_dir,
        candidate_loci_path=candidate_loci,
        mei_fasta=mei_fasta,
        out_path=out_tsv,
        reference_fasta=reference_fasta,
        tumor_bam_path=tumor_bam_depth,
        normal_bam_path=normal_bam_depth,
        rmsk_table_path=rmsk_table,
        g1k_mei_vcf=g1k_mei_vcf,
        g1k_split_padding_bp=g1k_split_padding_bp,
        g1k_dpe_padding_min_bp=g1k_dpe_padding_min_bp,
        g1k_dpe_padding_max_bp=g1k_dpe_padding_max_bp,
        g1k_dpe_padding_tlen_factor=g1k_dpe_padding_tlen_factor,
    )
    click.echo(f"[mei-annotate] done {out_path}")


if __name__ == "__main__":
    cli()
