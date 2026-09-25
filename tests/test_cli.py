"""CLI integration tests using Click's CliRunner (no real BAMs required).

All tests that exercise command bodies create minimal temp files to satisfy
``click.Path(exists=True)`` validators; no pysam / samtools calls are made
because the tests hit error-handling paths before any BAM access.
"""

from __future__ import annotations

import pytest
from click.testing import CliRunner

from retro_miner.cli import cli


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


@pytest.fixture()
def runner() -> CliRunner:
    return CliRunner()


# ---------------------------------------------------------------------------
# check-env
# ---------------------------------------------------------------------------


def test_check_env_exits_ok(runner: CliRunner) -> None:
    result = runner.invoke(cli, ["check-env"])
    assert result.exit_code == 0


def test_check_env_output(runner: CliRunner) -> None:
    result = runner.invoke(cli, ["check-env"])
    assert "retrotransposon-miner CLI is installed" in result.output
    assert "validate_environment.sh" in result.output


def test_check_env_help(runner: CliRunner) -> None:
    result = runner.invoke(cli, ["check-env", "--help"])
    assert result.exit_code == 0


# ---------------------------------------------------------------------------
# --help for remaining commands (smoke-test option parsing)
# ---------------------------------------------------------------------------


def test_extract_split_evidence_help(runner: CliRunner) -> None:
    result = runner.invoke(cli, ["extract-split-evidence", "--help"])
    assert result.exit_code == 0
    assert "--disease-bam" in result.output
    assert "--region" in result.output


def test_build_candidate_loci_help(runner: CliRunner) -> None:
    result = runner.invoke(cli, ["build-candidate-loci", "--help"])
    assert result.exit_code == 0
    assert "--evidence-dir" in result.output


def test_annotate_mei_support_help(runner: CliRunner) -> None:
    result = runner.invoke(cli, ["annotate-mei-support", "--help"])
    assert result.exit_code == 0
    assert "--mei-fasta" in result.output


# ---------------------------------------------------------------------------
# extract-split-evidence: region validation
#
# The validation fires at the top of the command body (before any BAM access),
# so we create minimal 1-byte files to satisfy click.Path(exists=True) and
# pass an empty/blank outdir path.
# ---------------------------------------------------------------------------


def _make_temp_bams(tmp_path):
    """Create two minimal placeholder files and a temp outdir."""
    disease_bam = tmp_path / "disease.bam"
    control_bam = tmp_path / "control.bam"
    disease_bam.write_bytes(b"\x00")
    control_bam.write_bytes(b"\x00")
    outdir = tmp_path / "out"
    return str(disease_bam), str(control_bam), str(outdir)


def test_regions_comma_only_raises(runner: CliRunner, tmp_path) -> None:
    """--regions ',' produces an empty region list → ClickException."""
    d, c, o = _make_temp_bams(tmp_path)
    result = runner.invoke(
        cli,
        [
            "extract-split-evidence",
            "--disease-bam", d,
            "--control-bam", c,
            "--outdir", o,
            "--regions", ",",
        ],
    )
    assert result.exit_code != 0
    assert "No valid regions" in result.output


def test_region_empty_string_raises(runner: CliRunner, tmp_path) -> None:
    """--region '' with no --regions produces an empty list → ClickException."""
    d, c, o = _make_temp_bams(tmp_path)
    result = runner.invoke(
        cli,
        [
            "extract-split-evidence",
            "--disease-bam", d,
            "--control-bam", c,
            "--outdir", o,
            "--region", "",
        ],
    )
    assert result.exit_code != 0
    assert "No valid regions" in result.output


def test_regions_whitespace_only_raises(runner: CliRunner, tmp_path) -> None:
    """--regions '  ,  , ' strips to empty → ClickException."""
    d, c, o = _make_temp_bams(tmp_path)
    result = runner.invoke(
        cli,
        [
            "extract-split-evidence",
            "--disease-bam", d,
            "--control-bam", c,
            "--outdir", o,
            "--regions", "  ,  , ",
        ],
    )
    assert result.exit_code != 0
    assert "No valid regions" in result.output


# ---------------------------------------------------------------------------
# annotate-mei-support: --disease-bam-depth XOR --control-bam-depth
# ---------------------------------------------------------------------------


def test_annotate_bam_depth_xor_raises(runner: CliRunner, tmp_path) -> None:
    """Providing --disease-bam-depth without --control-bam-depth should error."""
    evidence_dir = tmp_path / "evidence"
    evidence_dir.mkdir()
    candidate_loci = tmp_path / "loci.tsv"
    candidate_loci.write_text("")
    mei_fasta = tmp_path / "mei.fa"
    mei_fasta.write_text(">MEI\nACGT\n")
    bam_depth = tmp_path / "disease.depth.bam"
    bam_depth.write_bytes(b"\x00")
    out_tsv = tmp_path / "out.tsv"

    result = runner.invoke(
        cli,
        [
            "annotate-mei-support",
            "--evidence-dir", str(evidence_dir),
            "--candidate-loci", str(candidate_loci),
            "--mei-fasta", str(mei_fasta),
            "--out-tsv", str(out_tsv),
            "--disease-bam-depth", str(bam_depth),
            # intentionally omit --control-bam-depth
        ],
    )
    assert result.exit_code != 0
    assert "disease-bam-depth" in result.output or "control-bam-depth" in result.output


# ---------------------------------------------------------------------------
# annotate-genes
# ---------------------------------------------------------------------------


def test_annotate_genes_help_lists_options(runner):
    result = runner.invoke(cli, ["annotate-genes", "--help"])
    assert result.exit_code == 0
    for opt in ("--vcf", "--out", "--tsv", "--snpeff-genome", "--snpeff-bin", "--snpeff-config", "--snpeff-xmx"):
        assert opt in result.output
    assert "--backend" not in result.output
    assert "--batch-size" not in result.output


def test_annotate_genes_requires_existing_vcf(runner, tmp_path):
    result = runner.invoke(cli, ["annotate-genes", "--vcf", str(tmp_path / "nope.vcf"), "--out", str(tmp_path / "o.vcf")])
    assert result.exit_code != 0
    assert "does not exist" in result.output


def test_annotate_genes_offline(runner, tmp_path, monkeypatch):
    import functools
    import subprocess
    from pathlib import Path

    from retro_miner import gene_annotation as ga

    canned = (Path(__file__).parent / "data" / "chr22_mei.snpeff.vcf").read_text()
    seen = {}

    def fake_run(cmd, **kw):
        seen["cmd"] = cmd
        return subprocess.CompletedProcess(cmd, 0, stdout=canned, stderr="")

    monkeypatch.setattr(ga, "annotate_vcf", functools.partial(ga.annotate_vcf, run=fake_run))
    src = Path(__file__).parent / "data" / "chr22_mei.vcf"
    out = tmp_path / "out.vcf"
    result = runner.invoke(cli, ["annotate-genes", "--snpeff-genome", "GRCh38.115",
                                 "--snpeff-xmx", "4g", "--vcf", str(src), "--out", str(out)])
    assert result.exit_code == 0, result.output
    assert "records=30 annotated=30 with_gene=26 unmatched=0" in result.output
    assert seen["cmd"][:12] == [
        "snpEff", "-Xmx4g", "-noStats", "-noHgvs",
        "-spliceSiteSize", "0", "-spliceRegionExonSize", "0",
        "-spliceRegionIntronMin", "0", "-spliceRegionIntronMax", "0",
    ] and "GRCh38.115" in seen["cmd"]
    assert "GENE=CLTCL1" in out.read_text()
