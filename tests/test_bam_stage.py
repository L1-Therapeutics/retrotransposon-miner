"""Remote staging regression tests: fake transfers, no cloud credentials or calls."""
from __future__ import annotations

import json
import multiprocessing
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from retro_miner.bam_stage import (
    RemoteMetadata,
    _default_copy,
    _default_head_metadata,
    _default_head_size,
    _stage_lock,
    alignment_basename,
    apply_bam_stage,
    bam_stage_enabled,
    default_bam_stage_dir,
    index_sidecar_dests,
    index_sidecar_uris,
    is_remote_alignment_uri,
    local_copy_is_complete,
    main,
    plan_bam_stage,
    plan_staged_dests,
    required_free_bytes,
    should_stage_remote_bams,
    staged_alignment_dest,
    write_env_file,
)

URI = "s3://bucket/d.bam"


def _plan(stage: Path, uri: str = URI):
    return plan_bam_stage(disease_bam=uri, control_bam="local.bam", stage_dir=stage, chromosome_count=2)


def _metadata(uri: str) -> RemoteMetadata:
    index = uri.endswith((".bai", ".csi", ".crai"))
    return RemoteMetadata(5, etag='"index-v1"' if index else '"bam-v1"')


def _copy(_uri: str, dest: Path) -> None:
    dest.write_bytes(b"reads")


def _apply(stage: Path, **kwargs) -> Path:
    options = dict(head_metadata=_metadata, copy_object=_copy, disk_free=lambda _p: 10**12, log=lambda _m: None)
    options.update(kwargs)
    return Path(apply_bam_stage(_plan(stage), **options)["DISEASE_BAM"])


def test_remote_uri_predicates_and_staging_policy() -> None:
    for uri in (URI, "http://example.com/x.bam", "https://example.com/x.cram"):
        assert is_remote_alignment_uri(uri)
    for uri in (None, "", "/tmp/local.bam", "file:///tmp/local.bam"):
        assert not is_remote_alignment_uri(uri)
    assert should_stage_remote_bams(remote_bam_present=True, chromosome_count=2)
    assert should_stage_remote_bams(remote_bam_present=True, chromosome_count=1, chr_all=True)
    assert should_stage_remote_bams(remote_bam_present=True, chromosome_count=1, chr_concurrency=2)
    assert not should_stage_remote_bams(remote_bam_present=True, chromosome_count=1)
    assert not should_stage_remote_bams(remote_bam_present=False, chromosome_count=24)
    assert not should_stage_remote_bams(remote_bam_present=True, chromosome_count=24, enabled=False)
    assert not bam_stage_enabled("off")
    assert bam_stage_enabled("1")


def test_destinations_isolate_sources_across_plans(tmp_path: Path) -> None:
    a, b = "s3://tumor/sample.bam", "s3://normal/sample.bam"
    both = plan_staged_dests(tmp_path, [a, a, b])
    assert both[a] != both[b]
    assert both[a] == plan_staged_dests(tmp_path, [a])[a]
    assert both[b] == plan_staged_dests(tmp_path, [b])[b]
    assert both[a].name == both[b].name == "sample.bam"
    assert both[a].is_relative_to(tmp_path / "objects")
    assert not tmp_path.joinpath("objects").exists()  # planning is pure
    assert staged_alignment_dest(tmp_path, a, colliding_basenames={"sample.bam"}) == both[a]


def test_alignment_basename_and_index_candidates() -> None:
    assert alignment_basename(URI) == "d.bam"
    assert index_sidecar_uris(URI) == [URI + ".bai", "s3://bucket/d.bai", URI + ".csi", "s3://bucket/d.csi"]
    assert index_sidecar_uris("s3://bucket/x.cram") == ["s3://bucket/x.cram.crai", "s3://bucket/x.crai"]
    assert index_sidecar_uris("s3://bucket/sample?part.bam")[0] == "s3://bucket/sample?part.bam.bai"


@pytest.mark.parametrize("suffix,index", [(".bam", ".bai"), (".cram", ".crai")])
def test_http_index_candidates_preserve_query(suffix: str, index: str) -> None:
    uri = f"https://example.org/sample{suffix}?token=secret#download"
    assert index_sidecar_uris(uri)[:2] == [
        f"https://example.org/sample{suffix}{index}?token=secret#download",
        f"https://example.org/sample{index}?token=secret#download",
    ]


def test_local_copy_and_disk_budget(tmp_path: Path) -> None:
    dest = tmp_path / "x.bam"
    dest.write_bytes(b"abc")
    assert local_copy_is_complete(dest, 3)
    for size in (None, -1, 4):
        assert not local_copy_is_complete(dest, size)
    assert not local_copy_is_complete(tmp_path / "missing", 3)
    assert required_free_bytes([]) == 0
    assert required_free_bytes([1000], headroom_ratio=0.1, min_headroom_bytes=50) == 1100
    assert required_free_bytes([100], headroom_ratio=0.1, min_headroom_bytes=50) == 150


def test_plan_normalizes_remote_whitespace_and_keeps_local_paths(tmp_path: Path) -> None:
    plan = _plan(tmp_path, "  " + URI + "  ")
    assert list(plan.dest_by_uri) == [URI]
    assert plan.rewritten["DISEASE_BAM"] == str(plan.dest_by_uri[URI])
    assert plan.rewritten["CONTROL_BAM"] == "local.bam"
    single = plan_bam_stage(disease_bam=URI, control_bam="local.bam", chromosome_count=1)
    assert not single.should_stage
    assert apply_bam_stage(single, head_metadata=lambda _u: pytest.fail("must not HEAD")) == single.rewritten


def test_default_stage_dir(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    monkeypatch.delenv("RTM_BAM_STAGE_DIR", raising=False)
    monkeypatch.setenv("RTM_WORKDIR", str(tmp_path))
    assert default_bam_stage_dir() == tmp_path / "data" / "bam_stage"
    monkeypatch.setenv("RTM_BAM_STAGE_DIR", str(tmp_path / "custom"))
    assert default_bam_stage_dir() == tmp_path / "custom"


def test_verified_cache_reused_without_copy_or_free_space(tmp_path: Path) -> None:
    first = _apply(tmp_path)
    second = _apply(tmp_path, copy_object=lambda _u, _p: pytest.fail("cache must be reused"), disk_free=lambda _p: 0)
    assert first == second
    manifest = json.loads((first.parent / "manifest.json").read_text())
    assert manifest["remote"]["etag"] == '"bam-v1"'
    assert manifest["index_remote"]["etag"] == '"index-v1"'
    assert first.read_bytes() == b"reads"


@pytest.mark.parametrize("suffix", [".bam", ".cram"])
def test_same_basename_datasets_have_distinct_pairs(tmp_path: Path, suffix: str) -> None:
    sources = [f"s3://tumor/sample{suffix}", f"s3://normal/sample{suffix}"]
    paths = []
    for uri in sources:
        plan = _plan(tmp_path, uri)
        result = apply_bam_stage(plan, head_metadata=_metadata, copy_object=_copy,
                                 disk_free=lambda _p: 10**12, log=lambda _m: None)
        dest = Path(result["DISEASE_BAM"])
        assert any(p.exists() for p in index_sidecar_dests(dest))
        paths.append(dest)
    assert paths[0] != paths[1]
    assert all(p.read_bytes() == b"reads" for p in paths)


@pytest.mark.parametrize("field", ["etag", "version_id", "last_modified"])
def test_same_size_remote_change_preserves_previous_generation(tmp_path: Path, field: str) -> None:
    first = _apply(tmp_path)

    def changed(uri: str) -> RemoteMetadata:
        if uri == URI:
            values = dict(size=5, etag='"bam-v1"', version_id="", last_modified="")
            values[field] = "new-value"
            return RemoteMetadata(**values)
        return _metadata(uri)

    second = _apply(tmp_path, head_metadata=changed, copy_object=lambda _u, p: p.write_bytes(b"newer"))
    assert second != first
    assert first.read_bytes() == b"reads"
    assert second.read_bytes() == b"newer"


def test_changed_index_gets_new_pair_and_can_be_reused(tmp_path: Path) -> None:
    first = _apply(tmp_path)
    head = lambda uri: RemoteMetadata(5, etag='"index-v2"') if uri != URI else _metadata(uri)
    second = _apply(tmp_path, head_metadata=head)
    assert first != second
    assert first.exists()
    assert _apply(tmp_path, head_metadata=head, copy_object=lambda _u, _p: pytest.fail("reuse")) == second


@pytest.mark.parametrize("damage", ["missing", "truncated", "same_size", "bad_manifest", "path_traversal"])
def test_invalid_cache_is_not_trusted_or_overwritten(tmp_path: Path, damage: str) -> None:
    first = _apply(tmp_path)
    manifest_path = first.parent / "manifest.json"
    if damage == "missing":
        manifest_path.unlink()
    elif damage == "bad_manifest":
        manifest_path.write_text("not JSON")
    elif damage == "path_traversal":
        manifest = json.loads(manifest_path.read_text())
        manifest["index_name"] = "../../outside.bai"
        manifest_path.write_text(json.dumps(manifest))
    else:
        first.write_bytes(b"x" if damage == "truncated" else b"other")
    before = first.read_bytes()
    second = _apply(tmp_path)
    assert second != first
    assert first.read_bytes() == before
    assert second.read_bytes() == b"reads"
    assert _apply(tmp_path, copy_object=lambda _u, _p: pytest.fail("replacement should be reusable")) == second


@pytest.mark.parametrize("validator", ["", 'W/"weak"'])
def test_absent_or_weak_validators_never_reuse(tmp_path: Path, validator: str) -> None:
    head = lambda _u: RemoteMetadata(5, etag=validator, last_modified="yesterday")
    first = _apply(tmp_path, head_metadata=head)
    second = _apply(tmp_path, head_metadata=head)
    assert first != second
    assert first.exists() and second.exists()


def test_size_only_provider_is_supported_but_not_reusable(tmp_path: Path) -> None:
    def run():
        return apply_bam_stage(_plan(tmp_path), head_size=lambda _u: 5, copy_object=_copy,
                               disk_free=lambda _p: 10**12, log=lambda _m: None)
    assert run()["DISEASE_BAM"] != run()["DISEASE_BAM"]


def test_legacy_cache_is_preserved_but_never_adopted(tmp_path: Path) -> None:
    legacy = tmp_path / "d.bam"
    legacy.write_bytes(b"wrong")
    (tmp_path / "d.bam.bai").write_bytes(b"index")
    dest = _apply(tmp_path)
    assert dest != legacy
    assert legacy.read_bytes() == b"wrong"
    assert dest.read_bytes() == b"reads"


@pytest.mark.parametrize("failure", ["exception", "short_copy", "index_failure", "metadata_change"])
def test_failed_refresh_never_publishes_partial_generation(tmp_path: Path, failure: str) -> None:
    first = _apply(tmp_path)
    bam_heads = 0

    def head(uri: str) -> RemoteMetadata:
        nonlocal bam_heads
        if uri == URI:
            bam_heads += 1
            etag = '"v2"' if bam_heads == 1 or failure != "metadata_change" else '"v3"'
            return RemoteMetadata(5, etag=etag)
        return _metadata(uri)

    def copy(uri: str, dest: Path) -> None:
        if uri == URI:
            dest.write_bytes(b"x" if failure == "short_copy" else b"newer")
            if failure == "exception":
                raise RuntimeError("interrupted")
        else:
            dest.write_bytes(b"partial")
            if failure == "index_failure":
                raise RuntimeError("index unavailable")
            dest.write_bytes(b"index")

    with pytest.raises(RuntimeError, match="interrupted|size mismatch|Could not stage|changed during"):
        _apply(tmp_path, head_metadata=head, copy_object=copy)
    assert first.read_bytes() == b"reads"
    assert list(tmp_path.rglob("manifest.json")) == [first.parent / "manifest.json"]
    assert not list(tmp_path.rglob(".rtm-stage-*"))


@pytest.mark.parametrize("first_failure", ["empty", "exception", "wrong_size"])
def test_index_fallback_uses_matching_local_name(tmp_path: Path, first_failure: str) -> None:
    def copy(uri: str, dest: Path) -> None:
        if uri.endswith(".bam.bai"):
            dest.write_bytes(b"" if first_failure == "empty" else b"partial")
            if first_failure == "exception":
                raise RuntimeError("unavailable")
        else:
            dest.write_bytes(b"reads")
    dest = _apply(tmp_path, copy_object=copy)
    assert dest.with_suffix(".bai").read_bytes() == b"reads"
    assert not Path(str(dest) + ".bai").exists()


@pytest.mark.parametrize("size", [0, -1, None])
def test_invalid_size_fails_before_copy(tmp_path: Path, size: int | None) -> None:
    with pytest.raises(RuntimeError, match="positive size"):
        apply_bam_stage(_plan(tmp_path), head_size=lambda _u: size,
                        copy_object=lambda _u, _p: pytest.fail("invalid size must not copy"))


def test_low_disk_and_conflicting_providers_fail_before_copy(tmp_path: Path) -> None:
    with pytest.raises(RuntimeError, match="not enough free space"):
        _apply(tmp_path, disk_free=lambda _p: 10)
    with pytest.raises(ValueError, match="not both"):
        apply_bam_stage(_plan(tmp_path), head_size=lambda _u: 5, head_metadata=_metadata)


def test_metadata_parser_handles_s3_versions(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    seen = []

    def run(cmd, **kwargs):
        seen.append((cmd, kwargs))
        return subprocess.CompletedProcess(cmd, 0, json.dumps({
            "ContentLength": 5, "ETag": '"multipart-2"', "VersionId": "v42", "LastModified": "today",
        }), "")
    monkeypatch.setattr("retro_miner.bam_stage.shutil.which", lambda _name: "/fake/aws")
    monkeypatch.setattr("retro_miner.bam_stage.subprocess.run", run)
    assert _default_head_metadata(URI) == RemoteMetadata(5, '"multipart-2"', "today", "v42")
    assert _default_head_size(URI) == 5
    assert seen[0][1]["timeout"] == 60
    assert "--key" in seen[0][0]


@pytest.mark.parametrize("response", ["null", "[]", "invalid", '{"ContentLength":"bad"}'])
def test_malformed_s3_head_is_not_trusted(monkeypatch: pytest.MonkeyPatch, response: str) -> None:
    monkeypatch.setattr("retro_miner.bam_stage.shutil.which", lambda _name: "/fake/aws")
    monkeypatch.setattr("retro_miner.bam_stage.subprocess.run",
                        lambda *a, **kw: subprocess.CompletedProcess(a, 0, response, ""))
    assert _default_head_metadata(URI) is None


def test_http_head_uses_identity_encoding_and_validators(monkeypatch: pytest.MonkeyPatch) -> None:
    response = MagicMock()
    response.__enter__.return_value.headers = {"Content-Length": "5", "ETag": '"v1"', "Last-Modified": "today"}
    seen = []

    def open_url(request, **kwargs):
        seen.append(request)
        return response
    monkeypatch.setattr("retro_miner.bam_stage.urllib.request.urlopen", open_url)
    assert _default_head_metadata("https://example.org/d.bam") == RemoteMetadata(5, '"v1"', "today")
    assert seen[0].get_method() == "HEAD"
    assert seen[0].get_header("Accept-encoding") == "identity"


def test_default_copy_routes_to_s3_helper(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    seen = []
    monkeypatch.setattr("retro_miner.bam_stage.download_s3_uri", lambda u, p: seen.append((u, p)))
    _default_copy(URI, tmp_path / "d.bam")
    assert seen == [(URI, tmp_path / "d.bam")]


def test_reserved_manifest_filename_rejected(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="manifest.json"):
        staged_alignment_dest(tmp_path, "https://example.org/manifest.json")


def test_http_copy_errors_do_not_expose_signed_tokens(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr("retro_miner.bam_stage.shutil.which", lambda _n: "/fake/curl")

    def fail(_cmd):
        raise RuntimeError("curl failed: https://example.org/d.bam?token=VERY_SECRET")
    monkeypatch.setattr("retro_miner.bam_stage._run_cmd", fail)
    with pytest.raises(RuntimeError, match="HTTP alignment/index transfer failed") as error:
        _default_copy("https://example.org/d.bam?token=VERY_SECRET", tmp_path / "d.bam")
    assert "VERY_SECRET" not in str(error.value)
    assert error.value.__suppress_context__


def test_signed_urls_never_appear_in_manifest_or_logs(tmp_path: Path) -> None:
    uri = "https://example.org/d.bam?token=VERY_SECRET#fragment"
    logs = []
    result = apply_bam_stage(_plan(tmp_path, uri), head_metadata=_metadata, copy_object=_copy,
                             disk_free=lambda _p: 10**12, log=logs.append)
    manifest = (Path(result["DISEASE_BAM"]).parent / "manifest.json").read_text()
    assert "VERY_SECRET" not in manifest + "\n".join(logs)


def _process_stage(stage, queue, entered=None, release=None):
    """Spawn-safe fake transfer worker for real cross-process lock coverage."""
    copies = 0

    def copy(uri, dest):
        nonlocal copies
        copies += 1
        if entered is not None and uri == URI:
            entered.set()
            if not release.wait(10):
                raise RuntimeError("test release timed out")
        dest.write_bytes(b"reads")
    try:
        result = apply_bam_stage(_plan(Path(stage)), head_metadata=_metadata, copy_object=copy,
                                 disk_free=lambda _p: 10**12, log=lambda _m: None, lock_timeout=5)
        queue.put(("ok", result["DISEASE_BAM"], copies))
    except Exception as exc:  # noqa: BLE001 - report child errors to the parent
        queue.put(("error", str(exc), copies))


def test_concurrent_processes_share_one_complete_generation(tmp_path: Path) -> None:
    ctx = multiprocessing.get_context("spawn")
    queue, entered, release = ctx.Queue(), ctx.Event(), ctx.Event()
    first = ctx.Process(target=_process_stage, args=(str(tmp_path), queue, entered, release))
    second = ctx.Process(target=_process_stage, args=(str(tmp_path), queue))
    first.start()
    try:
        assert entered.wait(10)
        second.start()
        release.set()
        results = [queue.get(timeout=15), queue.get(timeout=15)]
        assert all(r[0] == "ok" for r in results), results
        assert results[0][1] == results[1][1]
        assert sorted(r[2] for r in results) == [0, 2]
    finally:
        release.set()
        for child in (first, second):
            if child.pid is not None:
                child.join(10)
                if child.is_alive():
                    child.terminate()
                    child.join(5)
        queue.close()
    assert first.exitcode == second.exitcode == 0


def test_lock_timeout_and_release(tmp_path: Path) -> None:
    def acquire():
        with _stage_lock(tmp_path, 0.1):
            pytest.fail("lock must not be acquired while held")
    with _stage_lock(tmp_path, 0):
        with ThreadPoolExecutor(max_workers=1) as pool:
            with pytest.raises(RuntimeError, match="Timed out"):
                pool.submit(acquire).result(timeout=5)
    with _stage_lock(tmp_path, 0):
        pass


def _hold_lock(stage, entered):
    import time
    with _stage_lock(Path(stage), 1):
        entered.set()
        time.sleep(30)


def test_kernel_releases_lock_after_process_crash(tmp_path: Path) -> None:
    ctx = multiprocessing.get_context("spawn")
    entered = ctx.Event()
    child = ctx.Process(target=_hold_lock, args=(str(tmp_path), entered))
    child.start()
    try:
        assert entered.wait(10)
    finally:
        child.terminate()
        child.join(5)
    with _stage_lock(tmp_path, 1):
        pass


@pytest.mark.parametrize("timeout", [-1, float("inf"), float("nan")])
def test_invalid_lock_timeout(tmp_path: Path, timeout: float) -> None:
    with pytest.raises(ValueError, match="finite"):
        _apply(tmp_path, lock_timeout=timeout)


def test_env_file_quotes_values_and_publishes_atomically(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    dest = tmp_path / "bam.env"
    dest.write_text("original\n")
    original_replace = Path.replace
    seen = []

    def replace(path, target):
        assert dest.read_text() == "original\n"
        seen.append(path)
        return original_replace(path, target)
    monkeypatch.setattr(Path, "replace", replace)
    write_env_file(dest, {"DISEASE_BAM": "/tmp/it's a.bam", "CONTROL_BAM": "$(no-execution)"})
    assert "DISEASE_BAM='/tmp/it'\"'\"'s a.bam'" in dest.read_text()
    assert "CONTROL_BAM='$(no-execution)'" in dest.read_text()
    assert len(seen) == 1 and not seen[0].exists()


def test_orphaned_partial_generation_is_never_adopted(tmp_path: Path) -> None:
    source = staged_alignment_dest(tmp_path, URI)
    orphan = source.parent / ".rtm-stage-crashed"
    orphan.mkdir(parents=True)
    (orphan / source.name).write_bytes(b"reads")
    (orphan / (source.name + ".bai")).write_bytes(b"index")
    dest = _apply(tmp_path)
    assert dest.parent != orphan
    assert (orphan / source.name).read_bytes() == b"reads"
    assert (dest.parent / "manifest.json").is_file()


def test_generation_published_as_one_complete_directory(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    original_rename = Path.rename
    published = []

    def rename(path, target):
        assert not target.exists()
        assert sorted(p.name for p in path.iterdir()) == ["d.bam", "d.bam.bai", "manifest.json"]
        published.append(target)
        return original_rename(path, target)
    monkeypatch.setattr(Path, "rename", rename)
    dest = _apply(tmp_path)
    assert published == [dest.parent]


def test_env_publication_failure_preserves_previous_output(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    dest = tmp_path / "bam.env"
    dest.write_text("original\n")

    def fail(_path, _target):
        raise OSError("rename failed")
    monkeypatch.setattr(Path, "replace", fail)
    with pytest.raises(OSError, match="rename failed"):
        write_env_file(dest, {"DISEASE_BAM": "new.bam"})
    assert dest.read_text() == "original\n"
    assert sorted(p.name for p in tmp_path.iterdir()) == ["bam.env"]


def test_failed_apply_does_not_replace_previous_environment(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    env = tmp_path / "bam.env"
    env.write_text("original\n")

    def fail(*_a, **_kw):
        raise RuntimeError("transfer failed")
    monkeypatch.setattr("retro_miner.bam_stage.apply_bam_stage", fail)
    with pytest.raises(RuntimeError, match="transfer failed"):
        main(["--disease-bam", URI, "--control-bam", "local.bam", "--stage-dir", str(tmp_path),
              "--chr-count", "2", "--out-env", str(env), "--apply"])
    assert env.read_text() == "original\n"


def test_cli_will_not_publish_unverified_planned_paths(tmp_path: Path) -> None:
    env = tmp_path / "bam.env"
    env.write_text("original\n")
    with pytest.raises(ValueError, match="--apply is required"):
        main(["--disease-bam", URI, "--control-bam", "local.bam", "--stage-dir", str(tmp_path),
              "--chr-count", "2", "--out-env", str(env)])
    assert env.read_text() == "original\n"
