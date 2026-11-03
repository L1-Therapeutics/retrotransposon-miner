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


def test_s3_uri_from_url_virtual_host_and_path(dl) -> None:
    cram = "https://1000genomes.s3.amazonaws.com/1000G_2504_high_coverage/data/ERR3240117/HG00100.final.cram"
    assert dl.s3_uri_from_url(cram) == (
        "s3://1000genomes/1000G_2504_high_coverage/data/ERR3240117/HG00100.final.cram"
    )
    assert dl.s3_uri_from_url(
        "https://s3.amazonaws.com/1000genomes/1000G_2504_high_coverage/data/ERR3240117/HG00100.final.cram.crai"
    ) == "s3://1000genomes/1000G_2504_high_coverage/data/ERR3240117/HG00100.final.cram.crai"
    assert dl.s3_uri_from_url(
        "https://1000genomes.s3.us-east-1.amazonaws.com/foo/bar.cram"
    ) == "s3://1000genomes/foo/bar.cram"
    assert dl.s3_uri_from_url("s3://1000genomes/foo/bar.cram") == "s3://1000genomes/foo/bar.cram"
    assert dl.s3_uri_from_url("https://ftp-trace.ncbi.nlm.nih.gov/foo.bam") is None


def test_mirror_prefers_s3_copy_for_s3_hosted_http(dl, monkeypatch: pytest.MonkeyPatch) -> None:
    copies: list[tuple[str, str]] = []
    streams: list[str] = []
    ds = dl.Dataset(
        dataset_id="hg00100_shortread_highcov_cram",
        category="test_bam",
        description="",
        source="1000 Genomes",
        url="https://1000genomes.s3.amazonaws.com/1000G_2504_high_coverage/data/ERR3240117/HG00100.final.cram",
        target_path="test_data/1kg_hg00100/chr22/hg00100.shortread.chr22.hg38.bam",
    )
    monkeypatch.setattr(dl, "_s3_head_size", lambda _uri: None)
    monkeypatch.setattr(
        dl,
        "_s3_copy",
        lambda src, dst: copies.append((src, dst))
        or {"status": "s3_copied", "src": src, "s3": dst, "bytes": 1, "seconds": 0.1},
    )
    monkeypatch.setattr(
        dl,
        "_http_stream_to_s3",
        lambda url, dest: streams.append(url) or {"status": "streamed_to_s3", "url": url, "s3": dest},
    )
    result = dl._mirror_full_alignment_to_s3(ds, "s3://l1tx-data/public", force=False)
    assert streams == []
    assert copies[0][0] == (
        "s3://1000genomes/1000G_2504_high_coverage/data/ERR3240117/HG00100.final.cram"
    )
    assert result["status"] == "s3_copied"


def test_s3_copy_uses_shared_high_concurrency_helper(dl, monkeypatch: pytest.MonkeyPatch) -> None:
    calls: list[tuple[str, str]] = []
    monkeypatch.setattr(dl, "copy_s3_uri", lambda src, dst: calls.append((src, dst)))
    monkeypatch.setattr(dl, "_s3_head_size", lambda _uri: 123)
    result = dl._s3_copy("s3://src/a.bam", "s3://dst/a.bam")
    assert calls == [("s3://src/a.bam", "s3://dst/a.bam")]
    assert result["status"] == "s3_copied"
    assert result["bytes"] == 123


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


def _write_discordant_pair_bam(dl, tmp_path: Path):
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
    return sorted_bam


def _alignment_keys(bam_path: Path) -> list[tuple[str, str, int, int]]:
    import pysam

    keys: list[tuple[str, str, int, int]] = []
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam:
            keys.append(
                (
                    str(read.query_name),
                    str(bam.get_reference_name(read.reference_id)),
                    int(read.flag),
                    int(read.reference_start),
                )
            )
    return keys


def test_index_aware_slice_recovers_discordant_mates(dl, tmp_path: Path) -> None:
    if not _have_pysam_and_samtools():
        pytest.skip("pysam and samtools are required for mate-fetch integration")

    sorted_bam = _write_discordant_pair_bam(dl, tmp_path)
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
    keys = _alignment_keys(out)
    names = {q for q, _c, _f, _p in keys}
    chroms = {c for _q, c, _f, _p in keys}
    assert names == {"disc"}
    assert chroms == {"chrA", "chrB"}
    assert len(keys) == 2
    assert len(set(keys)) == len(keys)
    qname_chrom = {(q, c) for q, c, _f, _p in keys}
    assert qname_chrom == {("disc", "chrA"), ("disc", "chrB")}
    assert len(keys) == len(qname_chrom)


def test_index_aware_slice_collapses_region_anchor_leaked_in_mate_window(
    dl, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Mate-window view that also dumps the region-chrom end must not 2× chrA."""
    if not _have_pysam_and_samtools():
        pytest.skip("pysam and samtools are required for mate-fetch integration")

    sorted_bam = _write_discordant_pair_bam(dl, tmp_path)
    orig = dl._samtools_view_indexed

    def _leak_other_end(url, index_path, region, out_bam, threads):
        if ":" not in str(region):
            return orig(url, index_path, region, out_bam, threads)
        import pysam

        with pysam.AlignmentFile(str(url), "rb") as src, pysam.AlignmentFile(
            str(out_bam), "wb", template=src
        ) as dst:
            for read in src:
                if read.query_name == "disc":
                    dst.write(read)

    monkeypatch.setattr(dl, "_samtools_view_indexed", _leak_other_end)
    out = tmp_path / "slice.bam"
    result = dl._slice_remote_alignment_with_mates(
        str(sorted_bam),
        "chrA",
        out,
        threads=1,
        force=True,
    )
    keys = _alignment_keys(out)
    qname_chrom = [(q, c) for q, c, _f, _p in keys]
    assert qname_chrom.count(("disc", "chrA")) == 1
    assert qname_chrom.count(("disc", "chrB")) == 1
    assert len(keys) == len(set(keys))
    assert "duplicate_alignments_dropped" in result


HG38_MEI_GENOTYPE_IDS = {
    "melt_1kg_phase3_sv_genotypes_vcf",
    "melt_1kg_phase3_sv_genotypes_vcf_tbi",
    "lr_1kg_ont_vienna_svim_asm_hg38_bcf",
    "lr_1kg_ont_vienna_svim_asm_hg38_bcf_csi",
    "lr_1kg_ont_vienna_svim_asm_hg38_svan_bcf",
    "lr_1kg_ont_vienna_svim_asm_hg38_svan_bcf_csi",
}


def test_hg38_download_includes_1kg_mei_genotype_callsets(dl) -> None:
    selected = dl._select_dataset_ids_for_references(("hg38",))
    missing = sorted(HG38_MEI_GENOTYPE_IDS - selected)
    assert missing == []
    assert "mm39_reference_fasta" not in selected
    assert "mm10_reference_fasta" not in selected

    cfg = Path(__file__).resolve().parents[1] / "resources" / "public_datasets.yaml"
    catalog = {ds.dataset_id: ds for ds in dl._load_config(cfg)}
    for ds_id in HG38_MEI_GENOTYPE_IDS:
        assert ds_id in catalog, ds_id
        assert catalog[ds_id].category == "mei_polymorphism"
        assert catalog[ds_id].required is True


HAZEN2016_ACCESSIONS = {
    "SRR3184315",
    "SRR3184320",
    "SRR3184323",
    "SRR3184308",
    "SRR3184309",
    "SRR3184313",
    "SRR3184375",
    "SRR3184314",
    "SRR3184307",
    "SRR3184310",
    "SRR3184327",
}
SUGO2025_ACCESSIONS = {f"DRR6979{n:02d}" for n in range(8, 31)}


def _catalog(dl):
    cfg = Path(__file__).resolve().parents[1] / "resources" / "public_datasets.yaml"
    return {ds.dataset_id: ds for ds in dl._load_config(cfg)}


def test_default_references_stay_human(dl) -> None:
    assert dl._resolve_selected_references(None) == ("hg19", "hg38", "hs1")
    assert dl._resolve_selected_references([]) == ("hg19", "hg38", "hs1")


def test_mm39_selection_adds_mouse_tracks_and_skips_human_smoke_bams(dl) -> None:
    selected = dl._select_dataset_ids_for_references(("mm39",))
    required = {
        "mm39_reference_fasta",
        "mm39_repeatmasker_rmsk",
        "gencode_vm39_gtf",
        "ferraj2023_merged_sv_vcf_zip",
        "mgp_v5_sv_insertions_bed",
        "mgp_v5_sv_insertions_bed_tbi",
        "dfam_human_families_embl",
        "mm10_segmental_duplications",
        "mm10_blacklist",
        "mm10_to_mm39_chain",
        "mm9_to_mm39_chain",
        "mm39_mappability_100bp_bw",
        "mm39_gap",
        "mm39_simple_repeat",
    }
    assert required <= selected
    assert "mm10_reference_fasta" not in selected
    assert "seqc2_control_bam" not in selected
    assert "gnomad_v41_sv_non_neuro_bb" not in selected
    assert "hazen2016_b2_r1" not in selected


def test_mouse_test_fastqs_match_hazen_and_sugo_libraries(dl) -> None:
    catalog = _catalog(dl)
    fastqs = [ds for ds in catalog.values() if ds.category == "test_fastq"]
    assert len(fastqs) == 68
    accessions = set()
    for ds in fastqs:
        assert ds.target_path.startswith("test_data/mouse_wgs/")
        accession = Path(ds.target_path).name.split("_", 1)[0]
        accessions.add(accession)
    assert accessions == HAZEN2016_ACCESSIONS | SUGO2025_ACCESSIONS
    assert catalog["hazen2016_b2_r1"].url.endswith("/SRR318/005/SRR3184315/SRR3184315_1.fastq.gz")
    assert "/DRR697/DRR697908/DRR697908_1.fastq.gz" in catalog["sugo2025_a1_r1"].url
    assert "Hazen2016/SRR3184308_1.fastq.gz" in catalog["hazen2016_c1_r1"].target_path

    ids = {ds.dataset_id for ds in fastqs}
    selected = dl._select_dataset_ids_for_references(
        ("mm39",),
        include_mouse_test_fastqs=True,
        mouse_test_fastq_ids=ids,
    )
    assert ids <= selected
    without = dl._select_dataset_ids_for_references(("mm39",))
    assert ids.isdisjoint(without)


def test_s3_sync_excludes_mouse_wgs(dl) -> None:
    source = Path(dl.__file__).read_text(encoding="utf-8")
    assert '"test_data/mouse_wgs/*"' in source
    assert '"*/test_data/mouse_wgs/*"' in source


def test_write_mgp_insertion_beds_splits_retroseq(dl, tmp_path: Path) -> None:
    import gzip

    src = tmp_path / "mgp.bed.gz"
    with gzip.open(src, "wt", encoding="utf-8") as handle:
        handle.write("#chr\tstart\tend\ttype\n")
        handle.write("1\t100\t300\tINS\t1:100-300;INS;NA;1/1:NA:NA:retroseq\n")
        handle.write("Y\t5\t15\tINS\t0/0\n")
    all_bed = tmp_path / "all.bed"
    mei_bed = tmp_path / "mei.bed"
    result = dl._write_mgp_insertion_beds(src, all_bed, mei_bed, force=True)
    assert result["rows"] == 2
    assert result["mei_rows"] == 1
    assert all_bed.read_text(encoding="utf-8").splitlines() == [
        "chr1\t100\t300\tINS",
        "chrY\t5\t15\tINS",
    ]
    assert mei_bed.read_text(encoding="utf-8") == "chr1\t100\t300\tINS\n"


def test_export_mouse_mei_subset_keeps_l1_sine_ltr(dl, tmp_path: Path) -> None:
    curated = tmp_path / "curated.fa"
    curated.write_text(
        ">L1Md_A#LINE/L1\nACGTAAAAAAAAAAAA\n"
        ">B1_Mus1#SINE/Alu\nTTTTTTTTTTGGGG\n"
        ">IAP#LTR/ERVK\nAAAAACCCCCTTTTT\n"
        ">B2_Mm1a#SINE/B2\nATGC\n"
        ">MER20#DNA/hAT\nGGGGAAAA\n",
        encoding="utf-8",
    )
    out = tmp_path / "mouse.fa"
    stats = dl._export_mouse_mei_subset(curated, out)
    text = out.read_text(encoding="utf-8")
    assert stats["mei_families"] == 4
    assert "L1Md_A#LINE/L1" in text
    assert "B1_Mus1#SINE/Alu" in text
    assert "IAP#LTR/ERVK" in text
    assert "MER20" not in text
    assert "AAAAAAAAAAAA" not in text
    assert text.index("B1_Mus1") < text.index("GGGG")
