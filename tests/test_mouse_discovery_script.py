"""The discovery wrapper accepts mm39 and runs on a mouse chromosome.

The public catalogs are not a MELT VCF. RetroSeq is a BED, and the Ferraj
zip's AppleDouble entries are not the insertion VCF.

The reference sequence and the BAM are a 20 kb slice of real mm10 chr1 from
nf-core/test-datasets, branch tfactivity (Lee et al., Nat Commun 2021, GEO
test subset). Source files:

https://github.com/nf-core/test-datasets/blob/tfactivity/reference/chr1.fa.gz
https://github.com/nf-core/test-datasets/blob/tfactivity/bams/p6/WT_p6_CONTROL_R1.mLb.clN.sorted.bam

Coordinates on the original chromosome are chr1:13270000-13289999. They are
shifted to 1-20000 on a contig still named chr1. This is GRCm38 sequence, not
GRCm39, and the BAM is a ChIP-seq control, not a whole-genome MEI library.
"""

from __future__ import annotations

import gzip
import importlib.util
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
WRAPPER = ROOT / "scripts" / "run_candidate_discovery_and_annotation.sh"
_DOWNLOAD = ROOT / "scripts" / "download_public_data.py"
_NFCORE = ROOT / "tests" / "fixtures" / "mouse" / "nfcore_tfactivity"
_NFCORE_FA = _NFCORE / "mm10_chr1_13270000_13290000.fa"
_NFCORE_BAM = _NFCORE / "mm10_chr1_13270000_13290000.bam"
# First bases of nf-core mm10 chr1:13270000-13289999. Guards against swapping
# the slice back for a synthetic repeat.
_NFCORE_PREFIX = "GGTGTTTCTAAATTAGCAAGGGTTTTCTTAGACTTACAGAAAATTTTCCCTTCCTGAT"


def _load_download_public_data():
    spec = importlib.util.spec_from_file_location("download_public_data_mouse", _DOWNLOAD)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def test_select_ferraj_vcf_skips_appledouble() -> None:
    dl = _load_download_public_data()
    paths = [
        Path("extracted/__MACOSX/._variants_freeze5_sv_insdel_sym.vcf.gz"),
        Path("extracted/merged_SV_VCF/variants_freeze5_sv_insdel_alt.vcf.gz"),
        Path("extracted/merged_SV_VCF/variants_freeze5_sv_insdel_sym.vcf.gz"),
    ]
    chosen = dl._select_ferraj_vcf(paths)
    assert chosen.name == "variants_freeze5_sv_insdel_sym.vcf.gz"


def test_retroseq_bed_keeps_distinct_known_mei_ids(tmp_path: Path) -> None:
    from retro_miner.mei_support import _build_known_mei_bed_from_intervals

    src = tmp_path / "mgpv5.SV_insertions.retroseq.bed"
    src.write_text("chr19\t150\t180\tINS\nchr19\t400\t420\tINS\n", encoding="utf-8")
    out = tmp_path / "known.bed"
    kept = _build_known_mei_bed_from_intervals(src, out)
    assert kept == 2
    rows = [line.split("\t") for line in out.read_text(encoding="utf-8").splitlines()]
    assert all(len(row) == 9 for row in rows)
    assert rows[0][3] == "chr19:150-180"
    assert rows[1][3] == "chr19:400-420"
    assert rows[0][3] != rows[1][3]


def _require_tools() -> None:
    missing = [name for name in ("samtools", "bedtools", "bwa") if shutil.which(name) is None]
    if missing:
        pytest.skip("discovery script needs " + ", ".join(missing))


def _write_gz(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(text)


def _mouse_public_tree(root: Path) -> None:
    """Smallest mm39-layout files the discovery script opens.

    The FASTA sequence is the nf-core mm10 chr1 slice, written where the mm39
    default path expects ``reference/mm39/mm39.fa``.
    """
    ref = root / "reference" / "mm39"
    ref.mkdir(parents=True)
    fasta = _NFCORE_FA.read_text(encoding="utf-8")
    if not fasta.startswith(">chr1 ") or _NFCORE_PREFIX not in fasta:
        raise AssertionError("nf-core mouse fixture is not the expected chr1 slice")
    (ref / "mm39.fa").write_text(fasta, encoding="utf-8")
    subprocess.run(["samtools", "faidx", str(ref / "mm39.fa")], check=True)

    # UCSC rmsk: chrom is column 6, then start, end, genoLeft, strand, name, class, family.
    _write_gz(
        root / "annotation" / "mm39" / "repeats" / "rmsk.txt.gz",
        "608\t0\t0\t0\t0\tchr1\t100\t400\t0\t+\tL1MdA_VI\tLINE\tL1\n",
    )
    _write_gz(
        root / "annotation" / "mm39" / "masks" / "gap.txt.gz",
        "585\tchr1\t0\t50\t1\tN\t50\ttelomere\tno\n",
    )
    ann = root / "annotation" / "mm39"
    (ann / "blacklist").mkdir(parents=True)
    (ann / "blacklist" / "blacklist.mm39.bed").write_text("chr1\t0\t10\n", encoding="utf-8")
    (ann / "segdup").mkdir()
    (ann / "segdup" / "genomicSuperDups.bed").write_text(
        "chr1\t20\t30\tchr1:20\n", encoding="utf-8"
    )
    (ann / "mappability").mkdir()
    (ann / "mappability" / "k100.Umap.MultiTrackMappability.bedGraph").write_text(
        "chr1\t0\t20000\t1\n", encoding="utf-8"
    )
    (ann / "mappability" / "k100.Umap.MultiTrackMappability.low_lt0.5.bed").write_text(
        "", encoding="utf-8"
    )
    (ann / "junk").mkdir()
    (ann / "junk" / "junk_exclusion_merged.bed").write_text("chr1\t0\t10\n", encoding="utf-8")
    _write_gz(
        ann / "gencode" / "gencode.vM39.annotation.gtf.gz",
        'chr1\tHAVANA\tgene\t100\t200\t.\t+\t.\tgene_id "ENSMUSG00000102693.2";\n',
    )
    mei = root / "retrotransposon_db" / "dfam"
    mei.mkdir(parents=True)
    (mei / "dfam_mouse_mei_l1_sine_erv.fasta").write_text(
        ">L1Md_A#LINE/L1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n"
        ">B1_Mus1#SINE/Alu\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n"
        ">B2_Mm1a#SINE/B2\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n"
        ">IAPEz-int#LTR/ERVK\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n"
        ">ETnERV3-int#LTR/ERV1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n",
        encoding="utf-8",
    )
    mgp = root / "polymorphism" / "mm39" / "mgp"
    mgp.mkdir(parents=True)
    (mgp / "mgpv5.SV_insertions.retroseq.bed").write_text(
        "chr1\t150\t180\tINS\n", encoding="utf-8"
    )


def _nfcore_bam(path: Path) -> None:
    """Copy the sliced nf-core BAM so the test does not rewrite the fixture index."""
    path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(_NFCORE_BAM, path)
    shutil.copyfile(Path(str(_NFCORE_BAM) + ".bai"), Path(str(path) + ".bai"))


def test_discovery_script_runs_on_mm39_with_nfcore_chr1(tmp_path: Path) -> None:
    _require_tools()
    public = tmp_path / "public"
    _mouse_public_tree(public)
    bam = tmp_path / "nfcore_chr1.bam"
    _nfcore_bam(bam)
    outdir = tmp_path / "results"
    proc = subprocess.run(
        [
            "bash",
            str(WRAPPER),
            "--reference-build",
            "mm39",
            "--chr",
            "chr1",
            "--disease-bam",
            str(bam),
            "--control-bam",
            str(bam),
            "--outdir",
            str(outdir),
            "--bwa-threads",
            "1",
            "--no-s3-results",
            "--no-local-assembly",
            "--no-empirical-stage",
        ],
        cwd=ROOT,
        env={
            **os.environ,
            "RTM_PUBLIC_DATA_DIR": str(public),
            "RTM_WORKDIR": str(tmp_path / "work"),
            "RTM_RESULTS_DIR": str(tmp_path / "results-default"),
            "PYTHONPATH": str(ROOT / "src"),
        },
        capture_output=True,
        text=True,
        timeout=180,
    )
    assert proc.returncode == 0, proc.stderr[-4000:] + proc.stdout[-2000:]
    params = (outdir / "pipeline_params.env").read_text(encoding="utf-8")
    assert "reference_build=mm39\n" in params
    assert "reference/mm39/mm39.fa" in params
    assert "gencode.vM39.annotation.gtf.gz" in params
    assert "stage=extract-split-evidence done region=chr1" in proc.stdout
    assert "stage=annotate-mei-support done region=chr1" in proc.stdout
    assert (outdir / "candidate_loci.mei.tsv").is_file()
    assert "unsupported --reference-build" not in proc.stderr
    reference = (public / "reference" / "mm39" / "mm39.fa").read_text(encoding="utf-8")
    assert _NFCORE_PREFIX in reference
    assert "ACGT" * 20 not in reference
