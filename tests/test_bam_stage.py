"""Unit tests for remote BAM staging predicates and dest-path logic (no live S3)."""

from __future__ import annotations

from pathlib import Path

import pytest

from retro_miner.bam_stage import (
    alignment_basename,
    apply_bam_stage,
    default_bam_stage_dir,
    index_sidecar_uris,
    is_remote_alignment_uri,
    local_copy_is_complete,
    plan_bam_stage,
    plan_staged_dests,
    required_free_bytes,
    should_stage_remote_bams,
    staged_alignment_dest,
    write_env_file,
)
from retro_miner.s3_transfer import DEFAULT_S3_MAX_CONCURRENCY


def test_is_remote_alignment_uri() -> None:
    assert is_remote_alignment_uri("s3://bucket/key.bam"), "s3:// URI should be remote"
    assert is_remote_alignment_uri("https://example.com/x.bam"), "https:// URI should be remote"
    assert is_remote_alignment_uri("http://example.com/x.bam"), "http:// URI should be remote"
    assert not is_remote_alignment_uri("/tmp/local.bam"), "local path should not be remote"
    assert not is_remote_alignment_uri("file:///tmp/local.bam"), "file:// URI should not be remote"
    assert not is_remote_alignment_uri(""), "empty string should not be remote"
    assert not is_remote_alignment_uri(None), "None should not be remote"


def test_should_stage_remote_bams_multi_chrom_and_all() -> None:
    assert should_stage_remote_bams(remote_bam_present=True, chromosome_count=6) is True, (
        "multi-chrom remote BAMs should be staged"
    )
    assert should_stage_remote_bams(remote_bam_present=True, chromosome_count=1, chr_all=True) is True, (
        "chr_all remote BAM should be staged"
    )
    assert should_stage_remote_bams(
        remote_bam_present=True, chromosome_count=1, chr_concurrency=6
    ) is True, "chr_concurrency should enable staging"
    assert should_stage_remote_bams(remote_bam_present=True, chromosome_count=1) is False, (
        "single-chrom remote BAM without concurrency should not be staged"
    )
    assert should_stage_remote_bams(remote_bam_present=False, chromosome_count=24) is False, (
        "local BAMs should not be staged"
    )
    assert should_stage_remote_bams(
        remote_bam_present=True, chromosome_count=24, enabled=False
    ) is False, "disabled staging should not stage"


def test_staged_alignment_dest_uses_basename() -> None:
    dest = staged_alignment_dest(
        "/work/data/bam_stage",
        "s3://l1tx-data/public/test_data/full/seqc2_disease_bam/WGS_EA_T_1.bwa.dedup.bam",
    )
    assert dest == Path("/work/data/bam_stage/WGS_EA_T_1.bwa.dedup.bam"), (
        f"expected basename dest, got {dest}"
    )
    http_dest = staged_alignment_dest(
        "/work/data/bam_stage",
        "https://example.com/path/control.final.cram",
    )
    assert http_dest == Path("/work/data/bam_stage/control.final.cram"), (
        f"expected basename dest for https, got {http_dest}"
    )


def test_plan_staged_dests_hash_suffix_on_basename_collision() -> None:
    stage = Path("/work/data/bam_stage")
    a = "s3://bucket-a/tumor/sample.bam"
    b = "s3://bucket-b/normal/sample.bam"
    planned = plan_staged_dests(stage, [a, a, b])
    assert planned[a] != planned[b], "colliding basenames should get different destinations"
    assert planned[a].name.startswith("sample."), (
        f"collision dest should start with 'sample.', got {planned[a].name}"
    )
    assert planned[a].name.endswith(".bam"), (
        f"collision dest should end with .bam, got {planned[a].name}"
    )
    assert planned[b].name.endswith(".bam"), (
        f"collision dest should end with .bam, got {planned[b].name}"
    )
    same = plan_staged_dests(stage, [a, a])
    assert same[a] == stage / "sample.bam", (
        f"duplicate URI should map to single dest, got {same[a]}"
    )


def test_alignment_basename_and_index_candidates() -> None:
    bam = "s3://l1tx-data/public/foo/WGS_EA_T_1.bwa.dedup.bam"
    assert alignment_basename(bam) == "WGS_EA_T_1.bwa.dedup.bam", (
        f"expected basename WGS_EA_T_1.bwa.dedup.bam, got {alignment_basename(bam)}"
    )
    cands = index_sidecar_uris(bam)
    assert bam + ".bai" in cands, "expected .bai candidate"
    assert bam[:-4] + ".bai" in cands, "expected .bai candidate without extension"
    assert bam + ".csi" in cands, "expected .csi candidate"
    cram = "https://example.com/x.cram"
    assert cram + ".crai" in index_sidecar_uris(cram), "expected .crai candidate for CRAM"


def test_local_copy_is_complete(tmp_path: Path) -> None:
    dest = tmp_path / "x.bam"
    dest.write_bytes(b"abc")
    assert local_copy_is_complete(dest, 3) is True, "3-byte file should be complete"
    assert local_copy_is_complete(dest, 4) is False, "3-byte file should be incomplete for size 4"
    assert local_copy_is_complete(dest, None) is False, "None size should be incomplete"
    assert local_copy_is_complete(tmp_path / "missing.bam", 3) is False, (
        "missing file should be incomplete"
    )


def test_required_free_bytes_adds_headroom() -> None:
    assert required_free_bytes([]) == 0, "empty list should require 0 bytes"
    need = required_free_bytes([1000], headroom_ratio=0.10, min_headroom_bytes=50)
    assert need == 1000 + max(100, 50), (
        f"expected 1100, got {need}"
    )
    big = required_free_bytes([100], headroom_ratio=0.10, min_headroom_bytes=50)
    assert big == 150, f"expected 150, got {big}"


def test_plan_bam_stage_rewrites_only_when_staging(tmp_path: Path) -> None:
    disease = "s3://bucket/d.bam"
    control = "s3://bucket/c.bam"
    plan = plan_bam_stage(
        disease_bam=disease,
        control_bam=control,
        disease_mate_bam=disease,
        control_mate_bam=control,
        stage_dir=tmp_path,
        chromosome_count=6,
    )
    assert plan.should_stage is True, "multi-chrom remote BAMs should be staged"
    assert plan.rewritten["DISEASE_BAM"] == str(tmp_path / "d.bam"), (
        f"expected staged disease BAM at {tmp_path / 'd.bam'}, got {plan.rewritten['DISEASE_BAM']}"
    )
    assert plan.rewritten["CONTROL_BAM"] == str(tmp_path / "c.bam"), (
        f"expected staged control BAM at {tmp_path / 'c.bam'}, got {plan.rewritten['CONTROL_BAM']}"
    )
    assert plan.rewritten["DISEASE_MATE_BAM"] == plan.rewritten["DISEASE_BAM"], (
        "disease mate BAM should reuse disease BAM path"
    )
    single = plan_bam_stage(
        disease_bam=disease,
        control_bam=control,
        stage_dir=tmp_path,
        chromosome_count=1,
        chr_concurrency=1,
    )
    assert single.should_stage is False, "single-chrom with concurrency=1 should not stage"
    assert single.rewritten["DISEASE_BAM"] == disease, (
        "non-staged plan should preserve original disease BAM URI"
    )


def test_default_bam_stage_dir_env(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    monkeypatch.delenv("RTM_BAM_STAGE_DIR", raising=False)
    monkeypatch.setenv("RTM_WORKDIR", str(tmp_path))
    assert default_bam_stage_dir() == tmp_path / "data" / "bam_stage", (
        "default stage dir should be under RTM_WORKDIR/data/bam_stage"
    )
    monkeypatch.setenv("RTM_BAM_STAGE_DIR", str(tmp_path / "custom"))
    assert default_bam_stage_dir() == tmp_path / "custom", (
        "RTM_BAM_STAGE_DIR should override default"
    )


def test_apply_bam_stage_skips_matching_local_copy(tmp_path: Path) -> None:
    disease_uri = "s3://bucket/d.bam"
    control_uri = "s3://bucket/c.bam"
    dest_d = tmp_path / "d.bam"
    dest_c = tmp_path / "c.bam"
    dest_d.write_bytes(b"D" * 10)
    dest_c.write_bytes(b"C" * 8)
    (tmp_path / "d.bam.bai").write_bytes(b"idxd")
    (tmp_path / "c.bam.bai").write_bytes(b"idxc")
    copies: list[str] = []
    plan = plan_bam_stage(
        disease_bam=disease_uri,
        control_bam=control_uri,
        stage_dir=tmp_path,
        chromosome_count=2,
    )
    rewritten = apply_bam_stage(
        plan,
        head_size=lambda uri: 10 if uri.endswith("d.bam") else 8,
        copy_object=lambda src, dest: copies.append(src),
        disk_free=lambda _p: 10**12,
        log=lambda _m: None,
    )
    assert copies == [], "complete local copies should not trigger staging"
    assert rewritten["DISEASE_BAM"] == str(dest_d), (
        f"expected existing local copy at {dest_d}, got {rewritten['DISEASE_BAM']}"
    )


def test_apply_bam_stage_copies_bam_and_index(tmp_path: Path) -> None:
    copies: list[str] = []

    def copy_fn(src: str, dest: Path) -> None:
        copies.append(src)
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(b"idx" if src.endswith((".bai", ".csi")) else b"x" * 5)

    plan = plan_bam_stage(
        disease_bam="s3://bucket/d.bam",
        control_bam="s3://bucket/c.bam",
        stage_dir=tmp_path,
        chromosome_count=2,
    )
    rewritten = apply_bam_stage(
        plan,
        head_size=lambda uri: 5 if uri.endswith(".bam") else 3,
        copy_object=copy_fn,
        disk_free=lambda _p: 10**12,
        log=lambda _m: None,
    )
    assert "s3://bucket/d.bam" in copies, "expected disease BAM to be staged"
    assert "s3://bucket/c.bam" in copies, "expected control BAM to be staged"
    assert any(c.endswith(".bai") or c.endswith(".csi") for c in copies), (
        "expected at least one index file to be staged"
    )
    assert Path(rewritten["DISEASE_BAM"]).is_file(), (
        f"staged disease BAM should exist at {rewritten['DISEASE_BAM']}"
    )
    assert (tmp_path / "d.bam.bai").is_file() or (tmp_path / "d.bai").is_file(), (
        "staged index file should exist"
    )


def test_apply_bam_stage_refuses_low_disk(tmp_path: Path) -> None:
    plan = plan_bam_stage(
        disease_bam="s3://bucket/d.bam",
        control_bam="s3://bucket/c.bam",
        stage_dir=tmp_path,
        chromosome_count=2,
    )
    with pytest.raises(RuntimeError, match="not enough free space"):
        apply_bam_stage(
            plan,
            head_size=lambda _uri: 1000,
            copy_object=lambda _s, _d: None,
            disk_free=lambda _p: 10,
            log=lambda _m: None,
        )


def test_write_env_file_quotes_paths(tmp_path: Path) -> None:
    env = tmp_path / "bam.env"
    write_env_file(env, {"DISEASE_BAM": "/tmp/a.bam", "CONTROL_BAM": "/tmp/b.bam"})
    text = env.read_text(encoding="utf-8")
    assert "DISEASE_BAM='/tmp/a.bam'" in text, "DISEASE_BAM should be single-quoted"
    assert "CONTROL_BAM='/tmp/b.bam'" in text, "CONTROL_BAM should be single-quoted"


def test_default_copy_uses_fast_s3_helper(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    seen: list[tuple[str, str]] = []

    def _fake_download(uri: str, dest) -> None:
        seen.append((uri, str(dest)))
        Path(dest).write_bytes(b"x")

    monkeypatch.setattr("retro_miner.bam_stage.download_s3_uri", _fake_download)
    from retro_miner.bam_stage import _default_copy

    dest = tmp_path / "x.bam"
    _default_copy("s3://bucket/x.bam", dest)
    assert seen == [("s3://bucket/x.bam", str(dest))], (
        f"expected single S3 download, got {seen}"
    )
    assert DEFAULT_S3_MAX_CONCURRENCY == 64, (
        f"expected S3 max concurrency 64, got {DEFAULT_S3_MAX_CONCURRENCY}"
    )
