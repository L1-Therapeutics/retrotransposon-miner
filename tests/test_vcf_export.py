import re
import tempfile
from pathlib import Path

import pysam
import pytest

from retro_miner.genotyper import GenotypeCall
from retro_miner.subfamily_voter import SubfamilyCall
from retro_miner.tsd_refiner import TSDResult
from retro_miner.vcf_export import write_mei_vcf

VCF_LINE_RE = re.compile(
    r"^(\S+)\t(\d+)\t(MEI_\S+)\t([A-ZN.])\t(<INS:MEI:[^>]+>)\t([.\d]+)\t"
    r"(PASS|LowQual)\t([^\t]+)\t(GT:GQ:VAF:AD)\t(\S+)$"
)


def _typed_record() -> dict:
    return {
        "chrom": "chr1",
        "pos": 15000,
        "ref_base": "A",
        "family": "L1",
        "genotype_call": GenotypeCall(
            genotype="0/1",
            vaf=0.50,
            log_likelihood_ratio=2.3134,
            genotype_quality=72.0,
            posterior_probs={"0/0": 0.0, "0/1": 1.0, "1/1": 0.0},
        ),
        "tsd_result": TSDResult(
            chrom="chr1",
            pos=15000,
            tsd_seq="ATTGCAG",
            tsd_length=7,
            tsd_confidence_score=0.85,
            polyA_tail_detected=True,
            entropy=1.4,
            method="tsd_flush_match",
        ),
        "subfamily_call": SubfamilyCall(
            family="L1",
            top_subfamily="L1HS",
            log_likelihood_ratio=3.5,
            posterior_prob=0.99,
            matching_kmers=4,
        ),
        "k_ref": 12,
        "k_alt": 13,
    }


def test_typed_record_header_and_body_formatting():
    with tempfile.TemporaryDirectory() as tmpdir:
        vcf_file = Path(tmpdir) / "test_output.vcf"
        write_mei_vcf([_typed_record()], vcf_file, sample_name="NA12878")

        assert vcf_file.exists()
        content = vcf_file.read_text()

        assert "##fileformat=VCFv4.3" in content
        assert "##contig=<ID=chr1>" in content
        assert "##INFO=<ID=SUBFAM" in content
        assert "##INFO=<ID=SUPPORT" in content
        assert "#CHROM\tPOS\tID\tREF\tALT" in content
        assert "NA12878" in content

        assert "chr1\t15000\tMEI_chr1_15000\tA\t<INS:MEI:L1HS>\t72.0\tPASS" in content
        assert "TSD=ATTGCAG" in content
        assert "TSDLEN=7" in content
        assert "POLYA=1" in content
        assert "SUBFAM=L1HS" in content
        assert "SUPPORT=13" in content
        assert "0/1:72.0:0.50:12,13" in content


def test_pysam_roundtrip_typed_record():
    with tempfile.TemporaryDirectory() as tmpdir:
        vcf_file = Path(tmpdir) / "typed.vcf"
        write_mei_vcf([_typed_record()], vcf_file, sample_name="NA12878")

        with pysam.VariantFile(str(vcf_file)) as reader:
            record = next(iter(reader))
            assert record.contig == "chr1"
            assert record.pos == 15000
            assert record.ref == "A"
            assert str(record.alts[0]) == "<INS:MEI:L1HS>"
            assert record.info["SVTYPE"] == "INS"
            assert record.info["MEI_TYPE"] == "L1HS"
            assert record.info["SUBFAM"] == "L1HS"
            assert record.info["SUPPORT"] == 13
            assert record.info["TSD"] == "ATTGCAG"
            assert record.info["TSDLEN"] == 7
            assert record.info["POLYA"] == 1
            sample = record.samples["NA12878"]
            assert sample["GT"] == (0, 1)
            assert sample["GQ"] == pytest.approx(72.0)
            assert sample["VAF"] == pytest.approx(0.50)
            assert sample["AD"] == (12, 13)


def test_flat_dict_record_backward_compatibility():
    record = {
        "chrom": "chr1",
        "pos": 15000,
        "ref_base": "A",
        "family": "L1",
        "subfamily": "L1HS",
        "genotype": "0/1",
        "genotype_quality": 75.0,
        "vaf": 0.48,
        "k_ref": 12,
        "k_alt": 11,
        "tsd_seq": "ATTGCAG",
        "tsd_len": 7,
        "poly_a_detected": True,
        "mei_llr": 4.25,
    }

    with tempfile.TemporaryDirectory() as tmpdir:
        vcf_file = Path(tmpdir) / "flat.vcf"
        write_mei_vcf([record], vcf_file)

        content = vcf_file.read_text()
        assert "chr1\t15000\tMEI_chr1_15000\tA\t<INS:MEI:L1HS>" in content
        assert "PASS" in content
        assert "0/1:75.0:0.48:12,11" in content
        assert "MEI_LLR=4.25" in content

        with pysam.VariantFile(str(vcf_file)) as reader:
            parsed = next(iter(reader))
            assert parsed.info["SUPPORT"] == 11
            assert parsed.samples["SAMPLE"]["AD"] == (12, 11)
            assert parsed.samples["SAMPLE"]["GQ"] == pytest.approx(75.0)


def test_low_gq_record_flagged_lowqual():
    record = {
        "chrom": "chr2",
        "pos": 20000,
        "genotype": "0/1",
        "genotype_quality": 12.5,
        "vaf": 0.15,
        "k_ref": 20,
        "k_alt": 3,
    }

    with tempfile.TemporaryDirectory() as tmpdir:
        vcf_file = Path(tmpdir) / "low_qual.vcf"
        write_mei_vcf([record], vcf_file)

        content = vcf_file.read_text()
        assert "chr2\t20000\tMEI_chr2_20000\tN\t<INS:MEI:L1HS>" in content
        assert "LowQual" in content
        assert "\tPASS\t" not in content


def test_high_gq_low_support_flagged_lowqual():
    record = {
        "chrom": "chr3",
        "pos": 30000,
        "genotype": "0/1",
        "genotype_quality": 55.0,
        "vaf": 0.95,
        "k_ref": 40,
        "k_alt": 2,
    }

    with tempfile.TemporaryDirectory() as tmpdir:
        vcf_file = Path(tmpdir) / "low_support.vcf"
        write_mei_vcf([record], vcf_file)

        content = vcf_file.read_text()
        assert "SUPPORT=2" in content
        assert "LowQual" in content
        assert "\tPASS\t" not in content


def test_no_poly_a_tail_exported_as_zero():
    record = dict(_typed_record())
    record["tsd_result"] = TSDResult(
        chrom="chr1",
        pos=15000,
        tsd_seq="GGCC",
        tsd_length=4,
        tsd_confidence_score=0.60,
        polyA_tail_detected=False,
        entropy=2.0,
        method="blunt",
    )

    with tempfile.TemporaryDirectory() as tmpdir:
        vcf_file = Path(tmpdir) / "no_polya.vcf"
        write_mei_vcf([record], vcf_file)

        assert "POLYA=0" in vcf_file.read_text()
        assert "POLYA=1" not in vcf_file.read_text()


def test_alt_definitions_match_used_families():
    with tempfile.TemporaryDirectory() as tmpdir:
        vcf_file = Path(tmpdir) / "l1_only.vcf"
        write_mei_vcf([_typed_record()], vcf_file)

        content = vcf_file.read_text()
        assert "##ALT=<ID=INS:MEI:L1HS" in content
        assert "##ALT=<ID=INS:MEI:ALU" not in content
        assert "##ALT=<ID=INS:MEI:SVA" not in content


def test_none_optional_fields_flat_record_exports_with_defaults():
    record = {
        "chrom": "chr1",
        "pos": 15000,
        "ref_base": None,
        "family": "L1",
        "subfamily": None,
        "genotype": None,
        "genotype_quality": None,
        "vaf": None,
        "k_ref": None,
        "k_alt": None,
        "tsd_seq": None,
        "tsd_len": None,
        "poly_a_detected": None,
        "mei_llr": None,
    }

    with tempfile.TemporaryDirectory() as tmpdir:
        vcf_file = Path(tmpdir) / "none_flat.vcf"
        write_mei_vcf([record], vcf_file)

        content = vcf_file.read_text()
        assert "chr1\t15000\tMEI_chr1_15000\tN\t<INS:MEI:L1HS>" in content
        assert "0/1:30.0:0.50:10,10" in content
        assert "MEI_LLR=3.00" in content
        assert "SUBFAM=L1HS" in content
        assert "TSDLEN=0" in content
        assert "POLYA=0" in content

        with pysam.VariantFile(str(vcf_file)) as reader:
            parsed = next(iter(reader))
            assert parsed.info["SUBFAM"] == "L1HS"
            assert parsed.info["MEI_LLR"] == pytest.approx(3.0)
            assert parsed.info["SUPPORT"] == 0
            sample = parsed.samples["SAMPLE"]
            assert sample["GT"] == (0, 1)
            assert sample["GQ"] == pytest.approx(30.0)
            assert sample["VAF"] == pytest.approx(0.50)
            assert sample["AD"] == (10, 10)


def test_explicit_zero_flat_fields_preserved():
    record = {
        "chrom": "chr8",
        "pos": 80000,
        "ref_base": "C",
        "family": "L1",
        "subfamily": "L1HS",
        "genotype": "0/1",
        "genotype_quality": 0.0,
        "vaf": 0.0,
        "k_ref": 0,
        "k_alt": 3,
        "tsd_seq": "GG",
        "tsd_len": 0,
        "poly_a_detected": False,
        "mei_llr": 0.0,
    }

    with tempfile.TemporaryDirectory() as tmpdir:
        vcf_file = Path(tmpdir) / "zero_flat.vcf"
        write_mei_vcf([record], vcf_file)

        content = vcf_file.read_text()
        assert "chr8\t80000\tMEI_chr8_80000\tC\t<INS:MEI:L1HS>" in content
        assert ".\tLowQual\t" in content
        assert "TSDLEN=0" in content
        assert "POLYA=0" in content
        assert "MEI_LLR=0.00" in content
        assert "0/1:0.0:0.00:0,3" in content

        with pysam.VariantFile(str(vcf_file)) as reader:
            parsed = next(iter(reader))
            assert parsed.info["MEI_LLR"] == pytest.approx(0.0)
            assert parsed.info["POLYA"] == 0
            assert parsed.info["SUPPORT"] == 3
            sample = parsed.samples["SAMPLE"]
            assert sample["GQ"] == pytest.approx(0.0)
            assert sample["VAF"] == pytest.approx(0.0)
            assert sample["AD"] == (0, 3)


def test_record_lines_conform_to_vcf_syntax():
    records = [
        _typed_record(),
        {
            "chrom": "chr2",
            "pos": 20000,
            "genotype": "0/1",
            "genotype_quality": 12.5,
            "vaf": 0.15,
            "k_ref": 20,
            "k_alt": 3,
        },
        {
            "chrom": "chr3",
            "pos": 30000,
            "genotype": "1/1",
            "genotype_quality": 55.0,
            "vaf": 0.95,
            "k_ref": 1,
            "k_alt": 45,
        },
    ]

    with tempfile.TemporaryDirectory() as tmpdir:
        vcf_file = Path(tmpdir) / "syntax.vcf"
        write_mei_vcf(records, vcf_file)

        data_lines = [ln for ln in vcf_file.read_text().splitlines() if not ln.startswith("#")]
        assert len(data_lines) == 3
        for line in data_lines:
            match = VCF_LINE_RE.match(line)
            assert match, f"VCF line failed structural regex: {line}"
            assert len(line.split("\t")) == 10