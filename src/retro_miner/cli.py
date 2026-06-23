from __future__ import annotations

from pathlib import Path
import time

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
@click.option("--disease-bam", type=click.Path(exists=True, dir_okay=False, path_type=Path), required=True)
@click.option("--control-bam", type=click.Path(exists=True, dir_okay=False, path_type=Path), required=True)
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
@click.option(
    "--require-strong-discordant-reason/--allow-weak-discordant-only",
    default=True,
    show_default=True,
    help=(
        "Require at least one strong discordant reason "
        "(mate_unmapped/interchrom/large_insert/poly_tail_anchor_rescue). "
        "Weak-only same_strand/improper_pair reads are filtered by default."
    ),
)
def extract_split_evidence_cmd(
    disease_bam: Path,
    control_bam: Path,
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
    require_strong_discordant_reason: bool,
) -> None:
    """Extract split-read MEI evidence and optional discordant evidence."""
    cmd_t0 = time.monotonic()
    outdir.mkdir(parents=True, exist_ok=True)
    region_list = [r.strip() for r in regions.split(",")] if regions else [region]
    region_list = [r for r in region_list if r]
    if not region_list:
        raise click.ClickException("No valid regions provided via --region/--regions.")

    split_summaries: list[ExtractionSummary] = []
    discordant_summaries: dict[str, ExtractionSummary] = {}
    for sample, bam in (("disease", disease_bam), ("control", control_bam)):
        sample_t0 = time.monotonic()
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
            f"passing={split_summary.passing_reads} split_rows={split_summary.split_evidence_rows} "
            f"elapsed={time.monotonic() - sample_t0:.1f}s"
        )
        if with_discordant:
            disc_t0 = time.monotonic()
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
                require_strong_discordant_reason=require_strong_discordant_reason,
            )
            discordant_summaries[sample] = discordant_summary
            click.echo(
                f"[done-discordant] sample={sample} scanned={discordant_summary.total_reads_scanned} "
                f"passing={discordant_summary.passing_reads} discordant_rows={discordant_summary.discordant_evidence_rows} "
                f"insert_threshold={discordant_summary.insert_size_threshold} "
                f"weak_only_filtered={discordant_summary.weak_only_discordant_filtered_rows} "
                f"elapsed={time.monotonic() - disc_t0:.1f}s"
            )

    summary_path = outdir / "split_evidence.summary.tsv"
    with summary_path.open("w", encoding="utf-8") as handle:
        handle.write(
            "sample\ttotal_reads_scanned\tpassing_reads\tsplit_evidence_rows\tdiscordant_evidence_rows\tinsert_size_threshold\tweak_only_discordant_filtered_rows\n"
        )
        for s in split_summaries:
            d = discordant_summaries.get(s.sample)
            discordant_rows = d.discordant_evidence_rows if d else 0
            insert_threshold = d.insert_size_threshold if d else 0
            weak_only_filtered = d.weak_only_discordant_filtered_rows if d else 0
            handle.write(
                f"{s.sample}\t{s.total_reads_scanned}\t{s.passing_reads}\t{s.split_evidence_rows}\t"
                f"{discordant_rows}\t{insert_threshold}\t{weak_only_filtered}\n"
            )
    click.echo(f"[summary] {summary_path}")
    click.echo(f"[extract] total_elapsed={time.monotonic() - cmd_t0:.1f}s")


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
    """Aggregate split/discordant evidence into disease-vs-control candidate loci."""
    t0 = time.monotonic()
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
    click.echo(f"[candidate-loci] {tsv_path} elapsed={time.monotonic() - t0:.1f}s")


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
    "--mei-full-fasta",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help=(
        "Optional full-consensus MEI FASTA for coordinate normalization "
        "(recommended: full-length L1HS/AluY/SVA canonical references). "
        "If omitted, an auto-selected canonical subset is generated from --mei-fasta."
    ),
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
    "--disease-bam-depth",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Optional disease BAM path for local BAM-depth controlization on consistent/junk-clean loci.",
)
@click.option(
    "--control-bam-depth",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Optional control BAM path for local BAM-depth controlization on consistent/junk-clean loci.",
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
    "--lr-mei-vcf",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Optional 1000G ONT Vienna long-read SVAN MEI VCF(.gz) for overlap annotation.",
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
@click.option(
    "--empirical-stage/--no-empirical-stage",
    default=False,
    show_default=True,
    help="Enable empirical BAM/context outlier scoring stage.",
)
@click.option(
    "--empirical-random-windows",
    type=int,
    default=1000,
    show_default=True,
    help="Number of random windows sampled per chromosome when scope=chromosome (or total when scope=genome).",
)
@click.option(
    "--empirical-random-scope",
    type=click.Choice(["chromosome", "genome"], case_sensitive=False),
    default="chromosome",
    show_default=True,
    help="Scope for random-window sampling used in empirical context scoring.",
)
@click.option(
    "--empirical-random-seed",
    type=int,
    default=13,
    show_default=True,
    help="Random seed for reproducible empirical window sampling.",
)
@click.option(
    "--empirical-highconf-bed",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Optional BED restricting empirical random windows to high-confidence intervals.",
)
@click.option(
    "--empirical-exclude-merged-bed",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=Path("data/public/annotation/hg38/junk/junk_exclusion_merged.bed"),
    show_default=True,
    help="Pre-merged junk exclusion BED for empirical random-window sampling.",
)
@click.option(
    "--empirical-exclude-segdup-bed",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=Path("data/public/annotation/hg38/segdup/genomicSuperDups.bed"),
    show_default=True,
    help="Optional segdup BED to exclude from empirical random-window sampling.",
)
@click.option(
    "--empirical-exclude-mappability-bedgraph",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=Path("data/public/annotation/hg38/mappability/k100.Umap.MultiTrackMappability.low_lt0.5.bed"),
    show_default=True,
    help="Optional low-mappability mask (BED preferred; bedGraph also accepted) for empirical sampling exclusion.",
)
@click.option(
    "--empirical-exclude-mappability-threshold",
    type=float,
    default=0.5,
    show_default=True,
    help="Mappability scores below this threshold are excluded from empirical sampling.",
)
@click.option(
    "--empirical-exclude-gap-bed",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=Path("data/public/annotation/hg38/masks/gap.txt.gz"),
    show_default=True,
    help="Optional gap BED to exclude from empirical random-window sampling.",
)
@click.option(
    "--empirical-exclude-blacklist-bed",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=Path("data/public/annotation/hg38/blacklist/ENCFF356LFX.bed.gz"),
    show_default=True,
    help="Optional blacklist BED to exclude from empirical random-window sampling.",
)
@click.option(
    "--empirical-cache-dir",
    type=click.Path(file_okay=False, path_type=Path),
    default=None,
    help="Optional cache directory for empirical random-window metrics (defaults to out-tsv directory/empirical_cache).",
)
@click.option(
    "--igv-plots/--no-igv-plots",
    default=True,
    show_default=True,
    help="Generate IGV snapshot PNGs for gold-review variants (requires reference + disease/control BAMs).",
)
@click.option(
    "--igv-top-n",
    type=int,
    default=0,
    show_default=True,
    help="Maximum number of prioritized gold-review variants to snapshot with IGV (<=0 means all).",
)
@click.option(
    "--igv-snapshot-dir",
    type=click.Path(file_okay=False, path_type=Path),
    default=None,
    help="Output directory for IGV snapshots (default: <out-tsv-stem>.gold_review.igv).",
)
@click.option(
    "--igv-launcher",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Path to igv.sh (default: search PATH and CONDA_PREFIX/bin).",
)
@click.option(
    "--igv-gold-only/--igv-all-tiers",
    default=True,
    show_default=True,
    help="Restrict IGV snapshots to gold-tier variants only.",
)
@click.option(
    "--igv-panel-height-min",
    type=int,
    default=250,
    show_default=True,
    help="Minimum IGV maxPanelHeight per snapshot (pixels).",
)
@click.option(
    "--igv-panel-height-max",
    type=int,
    default=8000,
    show_default=True,
    help="Maximum IGV maxPanelHeight per snapshot (pixels).",
)
@click.option(
    "--igv-timeout-sec",
    type=int,
    default=None,
    help="Optional timeout for the IGV batch run.",
)
@click.option(
    "--local-assembly/--no-local-assembly",
    default=False,
    show_default=True,
    help="Run per-silver-locus local assembly (SPAdes) for disease/control read stacks.",
)
@click.option(
    "--assembly-cache-dir",
    type=click.Path(file_okay=False, path_type=Path),
    default=None,
    help="Optional cache directory for local assembly artifacts (default: output dir/assembly_cache).",
)
@click.option(
    "--assembly-interval-pad-bp",
    type=int,
    default=250,
    show_default=True,
    help="Local-assembly interval half-padding around breakpoint/window center.",
)
@click.option(
    "--assembly-retry-pad-bp",
    type=int,
    default=600,
    show_default=True,
    help="Retry interval padding if initial local assembly fails.",
)
@click.option(
    "--assembly-max-reads-per-sample",
    type=int,
    default=600,
    show_default=True,
    help="Cap extracted reads per sample per locus for local assembly.",
)
@click.option(
    "--assembly-spades-threads",
    type=int,
    default=1,
    show_default=True,
    help="SPAdes threads per locus assembly.",
)
@click.option(
    "--assembly-spades-memory-gb",
    type=int,
    default=8,
    show_default=True,
    help="SPAdes memory limit (GB) per locus assembly.",
)
@click.option(
    "--assembly-minimap2-threads",
    type=int,
    default=1,
    show_default=True,
    help="minimap2 threads for per-locus contig alignments.",
)
@click.option(
    "--assembly-locus-workers",
    type=int,
    default=0,
    show_default=True,
    help="Number of loci to process in parallel during local assembly (0 = auto-size from CPU/RAM).",
)
@click.option(
    "--assembly-reuse-cache-only/--no-assembly-reuse-cache-only",
    default=False,
    show_default=True,
    help="Only reuse existing assembled cache entries; do not run SPAdes on cache miss.",
)
def annotate_mei_support_cmd(
    evidence_dir: Path,
    candidate_loci: Path,
    mei_fasta: Path,
    mei_full_fasta: Path | None,
    out_tsv: Path,
    reference_fasta: Path | None,
    disease_bam_depth: Path | None,
    control_bam_depth: Path | None,
    rmsk_table: Path | None,
    g1k_mei_vcf: Path | None,
    lr_mei_vcf: Path | None,
    g1k_split_padding_bp: int,
    g1k_dpe_padding_min_bp: int,
    g1k_dpe_padding_max_bp: int,
    g1k_dpe_padding_tlen_factor: float,
    empirical_stage: bool,
    empirical_random_windows: int,
    empirical_random_scope: str,
    empirical_random_seed: int,
    empirical_highconf_bed: Path | None,
    empirical_exclude_merged_bed: Path | None,
    empirical_exclude_segdup_bed: Path | None,
    empirical_exclude_mappability_bedgraph: Path | None,
    empirical_exclude_mappability_threshold: float,
    empirical_exclude_gap_bed: Path | None,
    empirical_exclude_blacklist_bed: Path | None,
    empirical_cache_dir: Path | None,
    igv_plots: bool,
    igv_top_n: int,
    igv_snapshot_dir: Path | None,
    igv_launcher: Path | None,
    igv_gold_only: bool,
    igv_panel_height_min: int,
    igv_panel_height_max: int,
    igv_timeout_sec: int | None,
    local_assembly: bool,
    assembly_cache_dir: Path | None,
    assembly_interval_pad_bp: int,
    assembly_retry_pad_bp: int,
    assembly_max_reads_per_sample: int,
    assembly_spades_threads: int,
    assembly_spades_memory_gb: int,
    assembly_minimap2_threads: int,
    assembly_locus_workers: int,
    assembly_reuse_cache_only: bool,
) -> None:
    """Annotate candidate loci with MEI family/subfamily support and insertion span estimates."""
    t0 = time.monotonic()
    click.echo("[mei-annotate] starting minimap2 clip-to-MEI alignment and locus annotation")
    if (disease_bam_depth is None) ^ (control_bam_depth is None):
        raise click.ClickException("Provide both --disease-bam-depth and --control-bam-depth, or neither.")
    out_path = annotate_candidate_loci_with_mei(
        evidence_dir=evidence_dir,
        candidate_loci_path=candidate_loci,
        mei_fasta=mei_fasta,
        out_path=out_tsv,
        reference_fasta=reference_fasta,
        disease_bam_path=disease_bam_depth,
        control_bam_path=control_bam_depth,
        rmsk_table_path=rmsk_table,
        g1k_mei_vcf=g1k_mei_vcf,
        lr_mei_vcf=lr_mei_vcf,
        g1k_split_padding_bp=g1k_split_padding_bp,
        g1k_dpe_padding_min_bp=g1k_dpe_padding_min_bp,
        g1k_dpe_padding_max_bp=g1k_dpe_padding_max_bp,
        g1k_dpe_padding_tlen_factor=g1k_dpe_padding_tlen_factor,
        empirical_stage=empirical_stage,
        empirical_random_windows=empirical_random_windows,
        empirical_random_scope=empirical_random_scope.lower(),
        empirical_random_seed=empirical_random_seed,
        empirical_highconf_bed=empirical_highconf_bed,
        empirical_exclude_merged_bed=empirical_exclude_merged_bed,
        empirical_exclude_segdup_bed=empirical_exclude_segdup_bed,
        empirical_exclude_mappability_bedgraph=empirical_exclude_mappability_bedgraph,
        empirical_exclude_mappability_threshold=empirical_exclude_mappability_threshold,
        empirical_exclude_gap_bed=empirical_exclude_gap_bed,
        empirical_exclude_blacklist_bed=empirical_exclude_blacklist_bed,
        empirical_cache_dir=empirical_cache_dir,
        igv_plots=igv_plots,
        igv_top_n=igv_top_n,
        igv_snapshot_dir=igv_snapshot_dir,
        igv_launcher=igv_launcher,
        igv_gold_only=igv_gold_only,
        igv_panel_height_min=igv_panel_height_min,
        igv_panel_height_max=igv_panel_height_max,
        igv_timeout_sec=igv_timeout_sec,
        local_assembly=local_assembly,
        assembly_cache_dir=assembly_cache_dir,
        assembly_interval_pad_bp=assembly_interval_pad_bp,
        assembly_retry_pad_bp=assembly_retry_pad_bp,
        assembly_max_reads_per_sample=assembly_max_reads_per_sample,
        assembly_spades_threads=assembly_spades_threads,
        assembly_spades_memory_gb=assembly_spades_memory_gb,
        assembly_minimap2_threads=assembly_minimap2_threads,
        assembly_locus_workers=assembly_locus_workers,
        assembly_reuse_cache_only=assembly_reuse_cache_only,
        mei_full_fasta=mei_full_fasta,
    )
    click.echo(f"[mei-annotate] done {out_path} elapsed={time.monotonic() - t0:.1f}s")


if __name__ == "__main__":
    cli()
