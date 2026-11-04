"""Mouse public-data helpers for the downloader.

Pure functions only: ENA FASTQ paths, UCSC chromosome names, and the
Mouse Genomes Project insertion table. The dataset catalog lives in
``resources/public_datasets.yaml``.
"""

from __future__ import annotations

import os

MOUSE_REFERENCES = ("mm10", "mm39")

# Active mouse ME classes in Dfam/RepeatMasker notation. B1 is filed under
# SINE/Alu. LTR covers IAP (ERVK), ETn (ERV1), and MaLR, which are the
# non-classical insertion classes the human L1/Alu/SVA subset does not have.
MOUSE_MEI_CLASS_RE = r"#(?:LINE/L1|SINE/|LTR/)"


def ena_fastq_subdir(accession: str) -> str:
    """ENA FASTQ directory under ``/vol1/fastq``.

    9-character accessions (DRR697908) use the 6-character prefix only.
    10-character accessions (SRR3184315) add a zero-padded last digit
    (``SRR318/005``).
    """
    acc = (accession or "").strip()
    if len(acc) < 6:
        raise ValueError(f"ENA accession is too short: {accession!r}")
    n = len(acc)
    if n <= 9:
        return acc[:6]
    if n == 10:
        return f"{acc[:6]}/00{acc[-1]}"
    if n == 11:
        return f"{acc[:6]}/0{acc[-2:]}"
    return f"{acc[:6]}/{acc[-3:]}"


def ena_fastq_url(accession: str, mate: int) -> str:
    """HTTPS URL for one mate of an ENA-submitted paired FASTQ."""
    if mate not in (1, 2):
        raise ValueError(f"ENA FASTQ mate must be 1 or 2, got {mate}")
    acc = accession.strip()
    sub = ena_fastq_subdir(acc)
    return f"https://ftp.sra.ebi.ac.uk/vol1/fastq/{sub}/{acc}/{acc}_{mate}.fastq.gz"


# Free public SRA runs in us-east-1. The FASTQs themselves are not in this bucket.
NCBI_SRA_BUCKET = "sra-pub-run-odp"

# gp3 on the mouse instance tops out at 1000 MiB/s. The CRT client sizes its
# one-file transfer from this target. 10 Gb/s is enough to fill that disk.
S3_SINGLE_FILE_TARGET_BANDWIDTH = "10Gb/s"

# ENA FASTQ records use this defline. Checked against SRR3184315 and DRR697908:
# the sequence, quality, and the BAM QNAME (the token before the first space)
# match a dump of the public SRA object.
ENA_MATCHING_SEQ_DEFLINE = "@$ac.$si $sn/$ri"
ENA_MATCHING_QUAL_DEFLINE = "+"

# fasterq-dump keeps a temporary copy of the uncompressed FASTQ beside the
# final files, and the SRA has to stay on disk while it reads. A WGS run
# expands to about four times the archive size per copy.
FASTERQ_SRA_SCRATCH_FACTOR = 9
BAM_SORT_HEADROOM_BYTES = 80 * 1024**3


def validate_run_accession(accession: str) -> str:
    """Return a bare SRA run accession."""
    acc = (accession or "").strip()
    if len(acc) < 9 or not acc.startswith(("SRR", "ERR", "DRR")) or not acc[3:].isdigit():
        raise ValueError(f"Not an SRA run accession: {accession!r}")
    return acc


def ncbi_sra_s3_uri(accession: str) -> str:
    """URI of the free ``.sra`` object for one run."""
    acc = validate_run_accession(accession)
    return f"s3://{NCBI_SRA_BUCKET}/sra/{acc}/{acc}"


def scratch_bytes_for_sra(sra_size: int) -> int:
    """Bytes of free disk to hold one run through dump and sort."""
    if sra_size < 0:
        raise ValueError(f"SRA size must be >= 0, got {sra_size}")
    return sra_size * FASTERQ_SRA_SCRATCH_FACTOR + BAM_SORT_HEADROOM_BYTES


def s3_single_file_config_text() -> str:
    """AWS CLI config that pushes one object as fast as the instance disk allows."""
    return (
        "[default]\n"
        "region = us-east-1\n"
        "s3 =\n"
        "    preferred_transfer_client = crt\n"
        f"    target_bandwidth = {S3_SINGLE_FILE_TARGET_BANDWIDTH}\n"
        "    multipart_chunksize = 64MB\n"
        "    max_concurrent_requests = 20\n"
    )


def s3_cp_command(src: str, dest: str, *, anonymous: bool) -> list[str]:
    """Copy one object. Anonymous is for the public SRA bucket."""
    cmd = ["aws", "s3", "cp", src, dest]
    if anonymous:
        cmd.append("--no-sign-request")
    return cmd


def fasterq_dump_command(
    sra_path: str,
    outdir: str,
    tempdir: str,
    threads: int,
) -> list[str]:
    """Dump paired FASTQ with the ENA read names."""
    return [
        "fasterq-dump",
        "--split-files",
        "--force",
        "--progress",
        "--seq-defline",
        ENA_MATCHING_SEQ_DEFLINE,
        "--qual-defline",
        ENA_MATCHING_QUAL_DEFLINE,
        "--threads",
        str(threads),
        "--temp",
        tempdir,
        "--outdir",
        outdir,
        sra_path,
    ]


# Pinned to the build under test. The mm39 index for this binary is
# s3://l1tx-data/mouse/public/reference/mm39-mem3/ and is not a bwa index.
BWA_MEM3_BIN = "bwa-mem3"
BWA_MEM3_VERSION = "0.16.0"
BWA_MEM3_INDEX_SUFFIXES = (".amb", ".ann", ".bwt.2bit.64", ".pac")
# bwa-mem3 ships mimalloc with purge disabled (-1 = never return pages).
# On a 64 GB align host the resident set grows for the whole run and the
# OOM killer cuts the SAM stream. 100 ms is mimalloc's usual purge delay.
BWA_MEM3_MIMALLOC_PURGE_DELAY_MS = "100"


def bwa_mem_command(
    reference: str,
    read1: str,
    read2: str,
    threads: int,
    sample: str,
) -> list[str]:
    """Align a pair with ``bwa-mem3``. The read group does not change coordinates.

    ``bwa-mem3`` rejects a real tab in ``-R`` and wants the two-character ``\\t`` escape.
    ``reference`` is the mem3 index prefix (``mm39``), not ``mm39.fa``.
    """
    read_group = f"@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA"
    return [
        BWA_MEM3_BIN,
        "mem",
        "-t",
        str(threads),
        "-R",
        read_group,
        reference,
        read1,
        read2,
    ]


def bwa_mem3_env(base: dict[str, str] | None = None) -> dict[str, str]:
    """Environment for a ``bwa-mem3`` subprocess, with mimalloc purge enabled."""
    env = dict(os.environ if base is None else base)
    env["MIMALLOC_PURGE_DELAY"] = BWA_MEM3_MIMALLOC_PURGE_DELAY_MS
    return env


def mouse_ucsc_chrom(token: str) -> str:
    """Map an Ensembl or UCSC mouse sequence name to a UCSC ``chr*`` name.

    Returns an empty string for headers and non-chromosomal tokens so callers
    can skip them. ``MT`` becomes ``chrM``.
    """
    text = (token or "").strip()
    if not text or text.startswith("#"):
        return ""
    if text.startswith("chr"):
        return text
    if text in {"MT", "M", "Mt"}:
        return "chrM"
    if text in {"X", "Y"} or text.isdigit():
        return "chr" + text
    return ""


def mgp_insertion_bed_fields(row: list[str]) -> tuple[str, int, int, str, bool] | None:
    """Parse one Mouse Genomes Project v5 insertion row.

    The release uses Ensembl chromosome names (``1``, not ``chr1``) on GRCm38.
    A row is treated as an MEI when any column mentions RetroSeq, which is the
    mobile-element caller in that table. Coordinates are kept as published.
    """
    if len(row) < 4 or row[0].startswith("#"):
        return None
    chrom = mouse_ucsc_chrom(row[0])
    if not chrom:
        return None
    try:
        start = int(row[1])
        end = int(row[2])
    except ValueError:
        return None
    if end <= start:
        return None
    kind = row[3] or "."
    mei = "retroseq" in "\t".join(row).lower()
    return chrom, start, end, kind, mei
