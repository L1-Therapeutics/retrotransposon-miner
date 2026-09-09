"""S3 cache URI helpers for scripts/download_public_data.py."""

from __future__ import annotations

import importlib.util
import shutil
import sys
from pathlib import Path

import pytest

_SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "download_public_data.py"


def _load_download_public_data():
    spec = importlib.util.spec_from_file_location("download_public_data", _SCRIPT)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def dl():
    return _load_download_public_data()


def test_normalize_s3_prefix_empty(dl) -> None:
    assert dl.normalize_s3_prefix("") == ""
    assert dl.normalize_s3_prefix(None) == ""
    assert dl.normalize_s3_prefix("   ") == ""


def test_normalize_s3_prefix_strips_slash_and_scheme(dl) -> None:
    assert dl.normalize_s3_prefix("s3://l1tx-data/") == "s3://l1tx-data"
    assert dl.normalize_s3_prefix("s3://l1tx-data/public/") == "s3://l1tx-data/public"
    assert dl.normalize_s3_prefix("l1tx-data/public") == "s3://l1tx-data/public"


def test_normalize_test_bam_mode(dl) -> None:
    assert dl.normalize_test_bam_mode("chr22") == "slice"
    assert dl.normalize_test_bam_mode("slice") == "slice"
    assert dl.normalize_test_bam_mode("full") == "full"
    assert dl.normalize_test_bam_mode("entire") == "full"
    with pytest.raises(ValueError):
        dl.normalize_test_bam_mode("both")


def test_split_s3_uri_and_full_alignment_key(dl) -> None:
    assert dl.split_s3_uri("s3://l1tx-data/public/foo.bam") == ("l1tx-data", "public/foo.bam")
    uri = dl.full_alignment_s3_uri(
        "s3://l1tx-data/public/",
        "seqc2_disease_bam",
        "https://example.com/path/WGS_EA_T_1.bwa.dedup.bam",
    )
    assert uri == "s3://l1tx-data/public/test_data/full/seqc2_disease_bam/WGS_EA_T_1.bwa.dedup.bam"


def test_local_slice_target_replaces_chr22(dl) -> None:
    ds = dl.Dataset(
        dataset_id="seqc2_disease_bam",
        category="test_bam",
        description="",
        source="NCBI",
        url="https://example.com/x.bam",
        target_path="test_data/seqc2/chr22/disease.chr22.hg38.bam",
        region="chr22",
    )
    assert dl._local_slice_target(ds, "chr1") == "test_data/seqc2/chr1/disease.chr1.hg38.bam"


def test_merge_interval_windows_and_samtools_region(dl) -> None:
    merged = dl.merge_interval_windows(
        [("chr2", 100, 200), ("chr2", 250, 300), ("chr3", 10, 20)],
        merge_gap=100,
    )
    assert ("chr2", 100, 300) in merged
    assert ("chr3", 10, 20) in merged
    split = dl.merge_interval_windows([("chr2", 100, 200), ("chr2", 500, 600)], merge_gap=100)
    assert split == [("chr2", 100, 200), ("chr2", 500, 600)]
    assert dl.samtools_region_1based("chr1", 0, 100) == "chr1:1-100"
    assert dl.samtools_region_1based("chrX", 999, 1200) == "chrX:1000-1200"


def test_http_index_sidecar_urls(dl) -> None:
    bam = "https://example.com/foo.bam"
    assert dl._http_index_sidecar_urls(bam)[0] == "https://example.com/foo.bai"
    assert "https://example.com/foo.bam.bai" in dl._http_index_sidecar_urls(bam)
    cram = "https://example.com/bar.final.cram"
    assert dl._http_index_sidecar_urls(cram)[0] == "https://example.com/bar.final.cram.crai"
    assert "https://example.com/bar.final.crai" in dl._http_index_sidecar_urls(cram)
    assert dl._sidecar_url(bam) == "https://example.com/foo.bai"


def test_index_candidates(dl, tmp_path: Path) -> None:
    bam = tmp_path / "sample.bam"
    assert dl._local_index_candidates(bam)[0] == Path(str(bam) + ".bai")
    s3_bam = "s3://l1tx-data/public/test_data/full/seqc2_disease_bam/WGS_EA_T_1.bwa.dedup.bam"
    cands = dl._s3_index_candidates(s3_bam)
    assert s3_bam + ".bai" in cands
    assert s3_bam[:-4] + ".bai" in cands


def test_remote_slice_without_index_does_not_name_scan(dl, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    calls: list[list[str]] = []

    def _spy(cmd: list[str], required: bool = True):
        calls.append(cmd)
        raise RuntimeError("samtools should not run without an index")

    monkeypatch.setattr(dl, "_run_cmd", _spy)
    with pytest.raises(RuntimeError, match="index"):
        dl._slice_remote_alignment_with_mates(
            "https://example.com/missing.bam",
            "chr1",
            tmp_path / "out.bam",
            threads=1,
            force=True,
        )
    assert calls == []
    assert not any("-N" in part for cmd in calls for part in cmd)


def _have_pysam_and_samtools() -> bool:
    if shutil.which("samtools") is None:
        return False
    try:
        import pysam  # noqa: F401
    except ImportError:
        return False
    return True


def test_index_aware_slice_recovers_discordant_mates(dl, tmp_path: Path) -> None:
    if not _have_pysam_and_samtools():
        pytest.skip("pysam and samtools are required for mate-fetch integration")
    import pysam

    raw = tmp_path / "raw.bam"
    sorted_bam = tmp_path / "sorted.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "SQ": [{"SN": "chrA", "LN": 10_000}, {"SN": "chrB", "LN": 2_000_000}],
    }
    seq = "A" * 50

    def _seg(bam, qname: str, chrom: str, pos: int, mate_chrom: str, mate_pos: int, is_read1: bool) -> None:
        a = pysam.AlignedSegment()
        a.query_name = qname
        a.query_sequence = seq
        a.flag = 1 | (64 if is_read1 else 128)
        a.reference_id = bam.get_tid(chrom)
        a.reference_start = pos
        a.next_reference_id = bam.get_tid(mate_chrom)
        a.next_reference_start = mate_pos
        a.mapping_quality = 60
        a.cigarstring = f"{len(seq)}M"
        a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
        bam.write(a)

    with pysam.AlignmentFile(str(raw), "wb", header=header) as bam:
        _seg(bam, "disc", "chrA", 100, "chrB", 5_000, True)
        _seg(bam, "disc", "chrB", 5_000, "chrA", 100, False)
        _seg(bam, "junk", "chrB", 1_000_000, "chrB", 1_000_100, True)

    ok, msg = dl._run_cmd(["samtools", "sort", "-o", str(sorted_bam), str(raw)], required=True)
    assert ok, msg
    ok, msg = dl._run_cmd(["samtools", "index", str(sorted_bam)], required=True)
    assert ok, msg

    out = tmp_path / "slice.bam"
    result = dl._slice_remote_alignment_with_mates(
        str(sorted_bam),
        "chrA",
        out,
        threads=1,
        force=True,
    )
    assert result["mate_fetch"] == "index"
    assert result["mate_qnames"] == 1
    assert result["mate_windows"] == 1
    names: set[str] = set()
    chroms: set[str] = set()
    with pysam.AlignmentFile(str(out), "rb") as bam:
        for read in bam:
            names.add(read.query_name)
            chroms.add(bam.get_reference_name(read.reference_id))
    assert names == {"disc"}
    assert chroms == {"chrA", "chrB"}
