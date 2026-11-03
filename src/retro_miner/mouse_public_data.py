"""Mouse public-data helpers for the downloader.

Pure functions only: ENA FASTQ paths, UCSC chromosome names, and the
Mouse Genomes Project insertion table. The dataset catalog lives in
``resources/public_datasets.yaml``.
"""

from __future__ import annotations

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
