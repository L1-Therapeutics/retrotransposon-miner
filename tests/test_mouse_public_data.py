"""ENA paths and Mouse Genomes Project row parsing for the mouse catalog."""

from __future__ import annotations

import pytest

from retro_miner.mouse_public_data import (
    ENA_MATCHING_QUAL_DEFLINE,
    ENA_MATCHING_SEQ_DEFLINE,
    BWA_MEM3_BIN,
    BWA_MEM3_INDEX_SUFFIXES,
    BWA_MEM3_MIMALLOC_PURGE_DELAY_MS,
    bwa_mem3_env,
    bwa_mem_command,
    ena_fastq_subdir,
    ena_fastq_url,
    fasterq_dump_command,
    mgp_insertion_bed_fields,
    mouse_ucsc_chrom,
    ncbi_sra_s3_uri,
    s3_cp_command,
    s3_single_file_config_text,
    scratch_bytes_for_sra,
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


def test_ncbi_sra_object_is_the_free_public_run() -> None:
    assert ncbi_sra_s3_uri("SRR3184315") == "s3://sra-pub-run-odp/sra/SRR3184315/SRR3184315"
    assert ncbi_sra_s3_uri(" DRR697908 ") == "s3://sra-pub-run-odp/sra/DRR697908/DRR697908"
    with pytest.raises(ValueError):
        ncbi_sra_s3_uri("SRR3184315_1")


def test_s3_copy_targets_one_file_at_disk_speed() -> None:
    config = s3_single_file_config_text()
    assert "preferred_transfer_client = crt" in config
    assert "target_bandwidth = 10Gb/s" in config
    cmd = s3_cp_command(
        "s3://sra-pub-run-odp/sra/SRR3184315/SRR3184315",
        "/tmp/SRR3184315",
        anonymous=True,
    )
    assert cmd[-1] == "--no-sign-request"
    signed = s3_cp_command("/tmp/a.bam", "s3://l1tx-data/mouse/align/mm39/a.bam", anonymous=False)
    assert "--no-sign-request" not in signed


def test_dump_and_align_commands_match_the_ena_reads() -> None:
    dump = fasterq_dump_command(
        "/tmp/SRR3184315",
        "/tmp/fq",
        "/tmp/scratch",
        16,
        total_ram_bytes=64 * 1024**3,
    )
    assert ENA_MATCHING_SEQ_DEFLINE in dump
    assert ENA_MATCHING_QUAL_DEFLINE in dump
    assert dump[dump.index("--threads") + 1] == "16"
    assert dump[dump.index("--mem") + 1] == "2G"
    assert dump[dump.index("--bufsize") + 1] == "8M"
    wider = fasterq_dump_command(
        "/tmp/SRR3184315",
        "/tmp/fq",
        "/tmp/scratch",
        16,
        total_ram_bytes=128 * 1024**3,
    )
    assert wider[wider.index("--mem") + 1] == "4G"
    assert wider[wider.index("--bufsize") + 1] == "16M"
    # MemTotal on the 16-thread align host, in kB from /proc/meminfo.
    this_host = fasterq_dump_command(
        "/tmp/SRR3184315",
        "/tmp/fq",
        "/tmp/scratch",
        16,
        total_ram_bytes=64776600 * 1024,
    )
    assert this_host[this_host.index("--mem") + 1] == "1976M"
    assert this_host[this_host.index("--bufsize") + 1] == "7M"
    align = bwa_mem_command("/ref/mm39-mem3/mm39", "/tmp/r1.fq", "/tmp/r2.fq", 14, "hazen2016_b2")
    assert align[:4] == [BWA_MEM3_BIN, "mem", "-t", "14"]
    assert BWA_MEM3_INDEX_SUFFIXES == (".amb", ".ann", ".bwt.2bit.64", ".pac")
    read_group = align[align.index("-R") + 1]
    assert read_group == "@RG\\tID:hazen2016_b2\\tSM:hazen2016_b2\\tPL:ILLUMINA"
    assert "\t" not in read_group
    assert align[-2:] == ["/tmp/r1.fq", "/tmp/r2.fq"]
    assert bwa_mem3_env({"PATH": "/usr/bin"}) == {
        "PATH": "/usr/bin",
        "MIMALLOC_PURGE_DELAY": BWA_MEM3_MIMALLOC_PURGE_DELAY_MS,
    }
    assert BWA_MEM3_MIMALLOC_PURGE_DELAY_MS == "100"


def test_scratch_reserves_dump_temp_and_sort_room() -> None:
    assert scratch_bytes_for_sra(1000) == 9000 + 80 * 1024**3
    with pytest.raises(ValueError):
        scratch_bytes_for_sra(-1)


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
