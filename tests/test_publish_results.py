"""Genome gold and plot upload selects the right files and bucket."""

from __future__ import annotations

from pathlib import Path

import pytest

from retro_miner.publish_results import (
    GENOME_GOLD_NAME,
    GENOME_VCF_NAME,
    IGV_DIR_NAME,
    READ_ARCH_DIR_NAME,
    bucket_uri,
    publish_outdir,
    resolve_bucket,
    results_destination,
    sync_argv,
)


def test_results_go_to_the_bucket_root_not_the_public_prefix() -> None:
    assert bucket_uri("s3://l1tx-data") == "s3://l1tx-data"
    assert bucket_uri("s3://l1tx-data/public") == "s3://l1tx-data"
    assert bucket_uri("l1tx-data") == "s3://l1tx-data"
    assert results_destination("s3://l1tx-data/public", "HG03086_germline_wgs") == (
        "s3://l1tx-data/results/HG03086_germline_wgs"
    )
    with pytest.raises(ValueError, match="empty"):
        bucket_uri("s3://")
    with pytest.raises(ValueError, match="single path"):
        results_destination("s3://l1tx-data", "a/b")


def test_resolve_bucket_prefers_flag_then_env(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("RTM_S3_BUCKET", "s3://from-env")
    monkeypatch.setenv("S3_BUCKET", "s3://other")
    assert resolve_bucket("s3://flag") == "s3://flag"
    assert resolve_bucket("") == "s3://from-env"
    monkeypatch.delenv("RTM_S3_BUCKET")
    assert resolve_bucket(None) == "s3://other"
    monkeypatch.delenv("S3_BUCKET")
    with pytest.raises(ValueError, match="RTM_S3_BUCKET"):
        resolve_bucket("")


def test_sync_keeps_root_gold_and_plots_only(tmp_path: Path) -> None:
    sample = tmp_path / "HG03086_germline_wgs"
    (sample / "chr1").mkdir(parents=True)
    (sample / GENOME_GOLD_NAME).write_text("genome\n", encoding="utf-8")
    (sample / "chr1" / GENOME_GOLD_NAME).write_text("per-chrom\n", encoding="utf-8")
    (sample / ".rtm_bam_stage.env").write_text("BAM=x\n", encoding="utf-8")
    igv = sample / "chr1" / IGV_DIR_NAME
    igv.mkdir()
    (igv / "rank001.png").write_bytes(b"png")
    arch = sample / "chr1" / READ_ARCH_DIR_NAME
    arch.mkdir()
    (arch / "rank001_read_arch.png").write_bytes(b"png")
    single_igv = sample / IGV_DIR_NAME
    single_igv.mkdir()
    (single_igv / "igv_snapshot_index.tsv").write_text("i\n", encoding="utf-8")

    dest = results_destination("s3://l1tx-data", sample.name)
    argv = sync_argv(sample, dest)
    assert argv[:4] == ["s3", "sync", str(sample), dest]
    assert "--exclude" in argv and "*" in argv
    assert GENOME_GOLD_NAME in argv
    assert GENOME_VCF_NAME in argv
    assert f"*/{IGV_DIR_NAME}/*" in argv
    assert f"*/{READ_ARCH_DIR_NAME}/*" in argv
    assert f"{IGV_DIR_NAME}/*" in argv
    assert "chr1/candidate_loci.mei.gold_review.tsv" not in argv
    assert ".rtm_bam_stage.env" not in argv

    seen: list[list[str]] = []
    published = publish_outdir(sample, "s3://l1tx-data/public", run_aws=seen.append)
    assert published == dest
    assert seen == [argv]
