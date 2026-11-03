"""ENA paths and Mouse Genomes Project row parsing for the mouse catalog."""

from __future__ import annotations

import pytest

from retro_miner.mouse_public_data import (
    ena_fastq_subdir,
    ena_fastq_url,
    mgp_insertion_bed_fields,
    mouse_ucsc_chrom,
)


def test_ena_fastq_url_ten_character_srr() -> None:
    url = ena_fastq_url("SRR3184315", 1)
    assert url == (
        "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR318/005/SRR3184315/SRR3184315_1.fastq.gz"
    )
    assert ena_fastq_subdir("SRR3184315") == "SRR318/005"


def test_ena_fastq_url_nine_character_drr_has_no_extra_subdir() -> None:
    url = ena_fastq_url("DRR697908", 2)
    assert url == "https://ftp.sra.ebi.ac.uk/vol1/fastq/DRR697/DRR697908/DRR697908_2.fastq.gz"
    assert ena_fastq_subdir("DRR697908") == "DRR697"


def test_ena_fastq_url_rejects_bad_mate() -> None:
    with pytest.raises(ValueError):
        ena_fastq_url("SRR3184315", 3)


def test_mouse_ucsc_chrom() -> None:
    assert mouse_ucsc_chrom("1") == "chr1"
    assert mouse_ucsc_chrom("chr19") == "chr19"
    assert mouse_ucsc_chrom("MT") == "chrM"
    assert mouse_ucsc_chrom("X") == "chrX"
    assert mouse_ucsc_chrom("#chr") == ""
    assert mouse_ucsc_chrom("GL456210.1") == ""


def test_mgp_insertion_bed_fields_marks_retroseq() -> None:
    row = [
        "1",
        "3021907",
        "3022107",
        "INS",
        "A/J",
        "1:3021907-3022107;INS;NA;1/1:NA:NA:retroseq",
    ]
    parsed = mgp_insertion_bed_fields(row)
    assert parsed == ("chr1", 3021907, 3022107, "INS", True)


def test_mgp_insertion_bed_fields_skips_header_and_non_mei() -> None:
    assert mgp_insertion_bed_fields(["#chr", "start", "end", "type"]) is None
    parsed = mgp_insertion_bed_fields(["2", "10", "30", "INS", "0/0"])
    assert parsed == ("chr2", 10, 30, "INS", False)
