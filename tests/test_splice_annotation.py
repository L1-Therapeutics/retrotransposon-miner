"""Unit tests for transcript-aware splice disruption, GTF splice-distance calculation,
isoform/MANE prioritization, and retrotransposon exonization / truncation modeling.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from retro_miner import gene_annotation as ga


def _write_synthetic_gtf(path: Path) -> Path:
    """Write a synthetic GTF with a + strand gene (2 exons) and a - strand gene (2 exons).

    Plus-strand gene ``GENE_POS`` (``ENSG00000100001``) on ``chr1``:
      * Canonical / MANE Select protein-coding transcript ``ENST00000100001.1``:
        - Exon 1: 1000..1100
        - Intron 1/1: 1101..1999 (donor at 1101..1102, acceptor at 1998..1999)
        - Exon 2: 2000..2100
      * Retained-intron isoform ``ENST00000100002.1`` spanning 1000..2100 as a single exon.

    Minus-strand gene ``GENE_NEG`` (``ENSG00000200001``) on ``chr1``:
      * Canonical / MANE Select protein-coding transcript ``ENST00000200001.1``:
        - Genomic lower exon (Exon 2 in 5'->3' transcription): 5000..5100
        - Intron 1/1: 5101..5999 (donor is near 6000, acceptor is near 5100)
        - Genomic upper exon (Exon 1 in 5'->3' transcription): 6000..6100
    """
    lines = [
        'chr1\thavana\tgene\t1000\t2100\t.\t+\t.\tgene_id "ENSG00000100001.1"; gene_name "GENE_POS"; gene_biotype "protein_coding";',
        'chr1\thavana\ttranscript\t1000\t2100\t.\t+\t.\tgene_id "ENSG00000100001.1"; transcript_id "ENST00000100001.1"; gene_name "GENE_POS"; transcript_biotype "protein_coding"; tag "MANE_Select"; tag "Ensembl_canonical";',
        'chr1\thavana\texon\t1000\t1100\t.\t+\t.\tgene_id "ENSG00000100001.1"; transcript_id "ENST00000100001.1"; exon_number "1"; transcript_biotype "protein_coding"; tag "MANE_Select";',
        'chr1\thavana\texon\t2000\t2100\t.\t+\t.\tgene_id "ENSG00000100001.1"; transcript_id "ENST00000100001.1"; exon_number "2"; transcript_biotype "protein_coding"; tag "MANE_Select";',
        'chr1\thavana\ttranscript\t1000\t2100\t.\t+\t.\tgene_id "ENSG00000100001.1"; transcript_id "ENST00000100002.1"; gene_name "GENE_POS"; transcript_biotype "retained_intron";',
        'chr1\thavana\texon\t1000\t2100\t.\t+\t.\tgene_id "ENSG00000100001.1"; transcript_id "ENST00000100002.1"; exon_number "1"; transcript_biotype "retained_intron";',
        'chr1\thavana\tgene\t5000\t6100\t.\t-\t.\tgene_id "ENSG00000200001.1"; gene_name "GENE_NEG"; gene_biotype "protein_coding";',
        'chr1\thavana\ttranscript\t5000\t6100\t.\t-\t.\tgene_id "ENSG00000200001.1"; transcript_id "ENST00000200001.1"; gene_name "GENE_NEG"; transcript_biotype "protein_coding"; tag "MANE_Select";',
        'chr1\thavana\texon\t5000\t5100\t.\t-\t.\tgene_id "ENSG00000200001.1"; transcript_id "ENST00000200001.1"; exon_number "2"; transcript_biotype "protein_coding"; tag "MANE_Select";',
        'chr1\thavana\texon\t6000\t6100\t.\t-\t.\tgene_id "ENSG00000200001.1"; transcript_id "ENST00000200001.1"; exon_number "1"; transcript_biotype "protein_coding"; tag "MANE_Select";',
    ]
    path.write_text("\n".join(lines) + "\n")
    return path


def _make_ann_entry(
    allele: str,
    annotation: str,
    impact: str,
    gene_name: str,
    gene_id: str,
    tx_id: str,
    biotype: str = "protein_coding",
    rank: str = "1/1",
) -> str:
    cols = [
        allele,
        annotation,
        impact,
        gene_name,
        gene_id,
        "transcript",
        tx_id,
        biotype,
        rank,
        "", "", "", "", "", "", "",
    ]
    return "|".join(cols)


def _run_canned_annotation(
    records: list[ga.VcfRecord],
    ann_by_pos: dict[int, str],
    gtf_path: Path | None,
) -> list[ga.VcfRecord]:
    headers = [
        "##fileformat=VCFv4.3",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
    ]
    out_lines = list(headers)
    for r in records:
        info_copy = dict(r.info)
        if r.pos in ann_by_pos:
            info_copy["ANN"] = ann_by_pos[r.pos]
        rec_copy = ga.VcfRecord(
            r.chrom, r.pos, r.id, r.ref, r.alt, r.qual, r.filter, info_copy, list(r.rest)
        )
        out_lines.append(rec_copy.to_line())
    stdout_text = "\n".join(out_lines) + "\n"

    def fake_run(cmd, **kw):
        return subprocess.CompletedProcess(cmd, 0, stdout=stdout_text, stderr="")

    annotated, _ = ga.annotate_records(
        records,
        headers,
        "GRCh38.99",
        gene_gtf=gtf_path,
        run=fake_run,
    )
    return annotated


# ---------------------------------------------------------------------------
# 1. Core Splice Disruption (+1 donor, -2 acceptor, plus and minus strands)
# ---------------------------------------------------------------------------


def test_core_splice_donor_and_acceptor_disruption(tmp_path: Path) -> None:
    gtf = _write_synthetic_gtf(tmp_path / "test.gtf")
    # Intron 1 of GENE_POS (+ strand) spans 1101..1999:
    #   +1 donor coordinate: 1101 (pos - 1100 == 1)
    #   -2 acceptor coordinate: 1998 (2000 - pos == 2)
    rec_donor = ga.VcfRecord(
        chrom="chr1", pos=1101, id="donor_ins", ref="N", alt="<INS:ME:ALU>",
        qual=".", filter="PASS",
        info={"SVTYPE": "INS", "END": "1101", "MEIFAMILY": "ALU", "ORIENT": "-"},
        rest=["GT:GQ", "./.:."],
    )
    rec_acceptor = ga.VcfRecord(
        chrom="chr1", pos=1998, id="acceptor_ins", ref="N", alt="<INS:ME:ALU>",
        qual=".", filter="PASS",
        info={"SVTYPE": "INS", "END": "1998", "MEIFAMILY": "ALU", "ORIENT": "-"},
        rest=["GT:GQ", "./.:."],
    )
    ann_map = {
        1101: _make_ann_entry(
            "<INS:ME:ALU>", "splice_donor_variant&intron_variant", "HIGH",
            "GENE_POS", "ENSG00000100001", "ENST00000100001.1", "protein_coding", "1/1",
        ),
        1998: _make_ann_entry(
            "<INS:ME:ALU>", "splice_acceptor_variant&intron_variant", "HIGH",
            "GENE_POS", "ENSG00000100001", "ENST00000100001.1", "protein_coding", "1/1",
        ),
    }
    annotated = _run_canned_annotation([rec_donor, rec_acceptor], ann_map, gtf)
    donor_out, acceptor_out = annotated[0], annotated[1]

    assert donor_out.info["CSQ"] == "splice_donor_variant"
    assert donor_out.info["SPLICE_SIDE"] == "donor"
    assert donor_out.info["SPLICE_DIST"] == "1"
    assert donor_out.info["DIST_TO_SPLICE_DONOR"] == "1"
    assert donor_out.info["DIST_TO_SPLICE_ACCEPTOR"] == "899"
    assert "SPLICE_REGION" in donor_out.info and donor_out.info["SPLICE_REGION"] is None

    assert acceptor_out.info["CSQ"] == "splice_acceptor_variant"
    assert acceptor_out.info["SPLICE_SIDE"] == "acceptor"
    assert acceptor_out.info["SPLICE_DIST"] == "2"
    assert acceptor_out.info["DIST_TO_SPLICE_DONOR"] == "898"
    assert acceptor_out.info["DIST_TO_SPLICE_ACCEPTOR"] == "2"
    assert "SPLICE_REGION" in acceptor_out.info and acceptor_out.info["SPLICE_REGION"] is None


def test_minus_strand_donor_and_acceptor_coordinates(tmp_path: Path) -> None:
    """Verify minus-strand reversal: donor is at upper genomic coordinate (5999), acceptor at lower (5102)."""
    gtf = _write_synthetic_gtf(tmp_path / "test.gtf")
    # For GENE_NEG (- strand), exon 1 is 6000..6100 and exon 2 is 5000..5100.
    # +1 of intron (donor) is at 5999 (6000 - 5999 == 1).
    # -2 of intron (acceptor) is at 5102 (5102 - 5100 == 2).
    rec_donor = ga.VcfRecord(
        chrom="chr1", pos=5999, id="neg_donor", ref="N", alt="<INS:ME:ALU>",
        qual=".", filter="PASS",
        info={"SVTYPE": "INS", "END": "5999", "MEIFAMILY": "ALU", "ORIENT": "+"},
        rest=["GT:GQ", "./.:."],
    )
    rec_acceptor = ga.VcfRecord(
        chrom="chr1", pos=5102, id="neg_acceptor", ref="N", alt="<INS:ME:ALU>",
        qual=".", filter="PASS",
        info={"SVTYPE": "INS", "END": "5102", "MEIFAMILY": "ALU", "ORIENT": "+"},
        rest=["GT:GQ", "./.:."],
    )
    ann_map = {
        5999: _make_ann_entry(
            "<INS:ME:ALU>", "intron_variant", "MODIFIER",
            "GENE_NEG", "ENSG00000200001", "ENST00000200001.1", "protein_coding", "1/1",
        ),
        5102: _make_ann_entry(
            "<INS:ME:ALU>", "intron_variant", "MODIFIER",
            "GENE_NEG", "ENSG00000200001", "ENST00000200001.1", "protein_coding", "1/1",
        ),
    }
    annotated = _run_canned_annotation([rec_donor, rec_acceptor], ann_map, gtf)
    donor_out, acceptor_out = annotated[0], annotated[1]

    assert donor_out.info["CSQ"] == "splice_donor_variant"
    assert donor_out.info["SPLICE_SIDE"] == "donor"
    assert donor_out.info["SPLICE_DIST"] == "1"
    assert "SPLICE_REGION" in donor_out.info

    assert acceptor_out.info["CSQ"] == "splice_acceptor_variant"
    assert acceptor_out.info["SPLICE_SIDE"] == "acceptor"
    assert acceptor_out.info["SPLICE_DIST"] == "2"
    assert "SPLICE_REGION" in acceptor_out.info


# ---------------------------------------------------------------------------
# 2. Splice Region Proximity (+5 of an intron)
# ---------------------------------------------------------------------------


def test_splice_region_proximity_plus_5(tmp_path: Path) -> None:
    gtf = _write_synthetic_gtf(tmp_path / "test.gtf")
    # Coordinate +5 of intron 1 (1100 + 5 = 1105)
    rec_region = ga.VcfRecord(
        chrom="chr1", pos=1105, id="region_ins", ref="N", alt="<INS:ME:ALU>",
        qual=".", filter="PASS",
        info={"SVTYPE": "INS", "END": "1105", "MEIFAMILY": "ALU", "ORIENT": "-"},
        rest=["GT:GQ", "./.:."],
    )
    ann_map = {
        1105: _make_ann_entry(
            "<INS:ME:ALU>", "splice_region_variant&intron_variant", "LOW",
            "GENE_POS", "ENSG00000100001", "ENST00000100001.1", "protein_coding", "1/1",
        ),
    }
    annotated = _run_canned_annotation([rec_region], ann_map, gtf)
    out = annotated[0]

    assert out.info["CSQ"] == "splice_region_variant"
    assert out.info["SPLICE_SIDE"] == "donor"
    assert out.info["SPLICE_DIST"] == "5"
    assert "SPLICE_REGION" in out.info and out.info["SPLICE_REGION"] is None


# ---------------------------------------------------------------------------
# 3. Deep Intronic vs. Exonization (+150 bp antisense vs sense Alu, plus LINE-1)
# ---------------------------------------------------------------------------


def test_intronic_alu_exonization_risk_antisense_vs_sense(tmp_path: Path) -> None:
    gtf = _write_synthetic_gtf(tmp_path / "test.gtf")
    # Coordinate +150 bp into intron 1 of GENE_POS (+ strand): 1100 + 150 = 1250
    rec_antisense = ga.VcfRecord(
        chrom="chr1", pos=1250, id="alu_as", ref="N", alt="<INS:ME:ALU>",
        qual=".", filter="PASS",
        info={"SVTYPE": "INS", "END": "1250", "MEIFAMILY": "ALU", "ORIENT": "-"},
        rest=["GT:GQ", "./.:."],
    )
    rec_sense = ga.VcfRecord(
        chrom="chr1", pos=1250, id="alu_sense", ref="N", alt="<INS:ME:ALU>",
        qual=".", filter="PASS",
        info={"SVTYPE": "INS", "END": "1250", "MEIFAMILY": "ALU", "ORIENT": "+"},
        rest=["GT:GQ", "./.:."],
    )
    ann_entry = _make_ann_entry(
        "<INS:ME:ALU>", "intron_variant", "MODIFIER",
        "GENE_POS", "ENSG00000100001", "ENST00000100001.1", "protein_coding", "1/1",
    )
    out_as = _run_canned_annotation([rec_antisense], {1250: ann_entry}, gtf)[0]
    out_s = _run_canned_annotation([rec_sense], {1250: ann_entry}, gtf)[0]

    assert out_as.info["CSQ"] == "intron_variant"
    assert out_as.info["SPLICE_DIST"] == "150"
    assert "SPLICE_REGION" not in out_as.info
    assert out_as.info["MEI_GENE_ORIENT"] == "antisense"
    assert out_as.info["MEI_EXONIZATION_RISK"] == "HIGH"
    assert out_as.info["ALU_EXONIZATION_RISK"] == "HIGH"

    assert out_s.info["CSQ"] == "intron_variant"
    assert out_s.info["SPLICE_DIST"] == "150"
    assert "SPLICE_REGION" not in out_s.info
    assert out_s.info["MEI_GENE_ORIENT"] == "sense"
    assert out_s.info["MEI_EXONIZATION_RISK"] == "NONE"
    assert out_s.info["ALU_EXONIZATION_RISK"] == "NONE"


def test_deep_intronic_alu_moderates_risk_beyond_500bp(tmp_path: Path) -> None:
    gtf = tmp_path / "deep.gtf"
    gtf.write_text(
        'chr1\thavana\tgene\t1000\t5000\t.\t+\t.\tgene_id "ENSG00000300001"; gene_name "DEEPG";\n'
        'chr1\thavana\texon\t1000\t1100\t.\t+\t.\tgene_id "ENSG00000300001"; transcript_id "ENST00000300001.1"; exon_number "1"; transcript_biotype "protein_coding"; tag "MANE_Select";\n'
        'chr1\thavana\texon\t4900\t5000\t.\t+\t.\tgene_id "ENSG00000300001"; transcript_id "ENST00000300001.1"; exon_number "2"; transcript_biotype "protein_coding"; tag "MANE_Select";\n'
    )
    # Breakpoint at 2000 -> 900 bp from donor (1100), 2900 bp from acceptor (4900)
    rec = ga.VcfRecord(
        chrom="chr1", pos=2000, id="alu_deep", ref="N", alt="<INS:ME:ALU>",
        qual=".", filter="PASS",
        info={"SVTYPE": "INS", "END": "2000", "MEIFAMILY": "ALU", "ORIENT": "-"},
        rest=["GT:GQ", "./.:."],
    )
    ann = _make_ann_entry(
        "<INS:ME:ALU>", "intron_variant", "MODIFIER",
        "DEEPG", "ENSG00000300001", "ENST00000300001.1", "protein_coding", "1/1",
    )
    out = _run_canned_annotation([rec], {2000: ann}, gtf)[0]
    assert out.info["SPLICE_DIST"] == "900"
    assert out.info["MEI_EXONIZATION_RISK"] == "MODERATE"
    assert out.info["ALU_EXONIZATION_RISK"] == "MODERATE"


def test_sense_line1_truncation_risk(tmp_path: Path) -> None:
    gtf = _write_synthetic_gtf(tmp_path / "test.gtf")
    rec_l1_sense = ga.VcfRecord(
        chrom="chr1", pos=1500, id="l1_sense", ref="N", alt="<INS:ME:LINE1>",
        qual=".", filter="PASS",
        info={
            "SVTYPE": "INS", "END": "1500", "MEIFAMILY": "LINE1",
            "ORIENT": "+", "MEI_5P": "1", "MEI_3P": "6018", "POLYA_MIN_BP": "45",
        },
        rest=["GT:GQ", "./.:."],
    )
    rec_l1_antisense = ga.VcfRecord(
        chrom="chr1", pos=1500, id="l1_antisense", ref="N", alt="<INS:ME:LINE1>",
        qual=".", filter="PASS",
        info={
            "SVTYPE": "INS", "END": "1500", "MEIFAMILY": "LINE1",
            "ORIENT": "-", "MEI_5P": "1", "MEI_3P": "6018", "POLYA_MIN_BP": "45",
        },
        rest=["GT:GQ", "./.:."],
    )
    ann = _make_ann_entry(
        "<INS:ME:LINE1>", "intron_variant", "MODIFIER",
        "GENE_POS", "ENSG00000100001", "ENST00000100001.1", "protein_coding", "1/1",
    )
    out_sense = _run_canned_annotation([rec_l1_sense], {1500: ann}, gtf)[0]
    out_as = _run_canned_annotation([rec_l1_antisense], {1500: ann}, gtf)[0]

    assert out_sense.info["MEI_GENE_ORIENT"] == "sense"
    assert out_sense.info["L1_TRUNCATION_RISK"] == "HIGH"
    assert out_sense.info["MEI_EXONIZATION_RISK"] == "HIGH"

    assert out_as.info["MEI_GENE_ORIENT"] == "antisense"
    assert out_as.info["L1_TRUNCATION_RISK"] == "NONE"
    assert out_as.info["MEI_EXONIZATION_RISK"] == "NONE"


def test_sva_exon_trapping_and_cryptic_splicing_sense_and_antisense(tmp_path: Path) -> None:
    """Verify SVA cryptic splice activation and exon trapping modeling in both
    sense (5' Alu-like/VNTR cryptic 3'SA + SINE-R AATAAA poly(A) trapping) and
    antisense (inverted poly(A) polypyrimidine tract + VNTR/Alu-like cryptic splice pairs)
    orientations, including deep-intronic moderation.
    """
    gtf = _write_synthetic_gtf(tmp_path / "test.gtf")
    # Near-exon intronic locus (pos=1250 -> SPLICE_DIST=150 <= 500 bp in GENE_POS on + strand)
    rec_sva_sense = ga.VcfRecord(
        chrom="chr1", pos=1250, id="sva_sense_full", ref="N", alt="<INS:ME:SVA>",
        qual=".", filter="PASS",
        info={
            "SVTYPE": "INS", "END": "1250", "MEIFAMILY": "SVA",
            "ORIENT": "+", "MEI_5P": "1", "MEI_3P": "1366", "MEI_SPAN": "1366", "POLYA_MIN_BP": "50",
        },
        rest=["GT:GQ", "./.:."],
    )
    rec_sva_antisense_near = ga.VcfRecord(
        chrom="chr1", pos=1250, id="sva_antisense_near", ref="N", alt="<INS:ME:SVA>",
        qual=".", filter="PASS",
        info={
            "SVTYPE": "INS", "END": "1250", "MEIFAMILY": "SVA",
            "ORIENT": "-", "MEI_5P": "1", "MEI_3P": "1366", "MEI_SPAN": "1366", "POLYA_MIN_BP": "50",
        },
        rest=["GT:GQ", "./.:."],
    )
    ann_near = _make_ann_entry(
        "<INS:ME:SVA>", "intron_variant", "MODIFIER",
        "GENE_POS", "ENSG00000100001", "ENST00000100001.1", "protein_coding", "1/1",
    )
    out_sense = _run_canned_annotation([rec_sva_sense], {1250: ann_near}, gtf)[0]
    out_as_near = _run_canned_annotation([rec_sva_antisense_near], {1250: ann_near}, gtf)[0]

    assert out_sense.info["MEI_GENE_ORIENT"] == "sense"
    assert out_sense.info["SVA_EXON_TRAPPING_RISK"] == "HIGH"
    assert out_sense.info["MEI_EXONIZATION_RISK"] == "HIGH"

    assert out_as_near.info["MEI_GENE_ORIENT"] == "antisense"
    assert out_as_near.info["SVA_EXON_TRAPPING_RISK"] == "HIGH"
    assert out_as_near.info["MEI_EXONIZATION_RISK"] == "HIGH"

    # Deep intronic GTF (>500 bp from splice junctions)
    deep_gtf = tmp_path / "deep_sva.gtf"
    deep_gtf.write_text(
        'chr1\thavana\tgene\t1000\t5000\t.\t+\t.\tgene_id "ENSG00000300001"; gene_name "DEEPG";\n'
        'chr1\thavana\texon\t1000\t1100\t.\t+\t.\tgene_id "ENSG00000300001"; transcript_id "ENST00000300001.1"; exon_number "1"; transcript_biotype "protein_coding"; tag "MANE_Select";\n'
        'chr1\thavana\texon\t4900\t5000\t.\t+\t.\tgene_id "ENSG00000300001"; transcript_id "ENST00000300001.1"; exon_number "2"; transcript_biotype "protein_coding"; tag "MANE_Select";\n'
    )
    ann_deep = _make_ann_entry(
        "<INS:ME:SVA>", "intron_variant", "MODIFIER",
        "DEEPG", "ENSG00000300001", "ENST00000300001.1", "protein_coding", "1/1",
    )
    # Antisense SVA deep in intron (pos=2500 -> SPLICE_DIST=1400 > 500) -> MODERATE
    rec_sva_as_deep = ga.VcfRecord(
        chrom="chr1", pos=2500, id="sva_as_deep", ref="N", alt="<INS:ME:SVA>",
        qual=".", filter="PASS",
        info={
            "SVTYPE": "INS", "END": "2500", "MEIFAMILY": "SVA",
            "ORIENT": "-", "MEI_5P": "1", "MEI_3P": "1366", "MEI_SPAN": "1366",
        },
        rest=["GT:GQ", "./.:."],
    )
    # Severely 5'-truncated and 3'-truncated sense SVA deep in intron -> MODERATE
    rec_sva_sense_trunc_deep = ga.VcfRecord(
        chrom="chr1", pos=2500, id="sva_sense_trunc_deep", ref="N", alt="<INS:ME:SVA>",
        qual=".", filter="PASS",
        info={
            "SVTYPE": "INS", "END": "2500", "MEIFAMILY": "SVA",
            "ORIENT": "+", "MEI_5P": "700", "MEI_3P": "950", "MEI_SPAN": "250", "POLYA_MIN_BP": "0",
        },
        rest=["GT:GQ", "./.:."],
    )
    out_as_deep = _run_canned_annotation([rec_sva_as_deep], {2500: ann_deep}, deep_gtf)[0]
    out_sense_trunc = _run_canned_annotation([rec_sva_sense_trunc_deep], {2500: ann_deep}, deep_gtf)[0]

    assert out_as_deep.info["SPLICE_DIST"] == "1400"
    assert out_as_deep.info["SVA_EXON_TRAPPING_RISK"] == "MODERATE"
    assert out_as_deep.info["MEI_EXONIZATION_RISK"] == "MODERATE"

    assert out_sense_trunc.info["SPLICE_DIST"] == "1400"
    assert out_sense_trunc.info["SVA_EXON_TRAPPING_RISK"] == "MODERATE"
    assert out_sense_trunc.info["MEI_EXONIZATION_RISK"] == "MODERATE"


# ---------------------------------------------------------------------------
# 4. Isoform Awareness & MANE Select / Protein-Coding Prioritization
# ---------------------------------------------------------------------------


def test_isoform_awareness_protein_coding_takes_precedence_over_retained_intron(tmp_path: Path) -> None:
    """When an insertion overlaps both a protein_coding MANE transcript (intron_variant)
    and a retained_intron transcript (non_coding_transcript_exon_variant), the
    protein_coding / MANE transcript's consequence takes precedence in CSQ while
    TX_ID and TX_BIOTYPE retain both isoforms.
    """
    gtf = _write_synthetic_gtf(tmp_path / "test.gtf")
    rec = ga.VcfRecord(
        chrom="chr1", pos=1400, id="isoform_ins", ref="N", alt="<INS:ME:ALU>",
        qual=".", filter="PASS",
        info={"SVTYPE": "INS", "END": "1400", "MEIFAMILY": "ALU", "ORIENT": "-"},
        rest=["GT:GQ", "./.:."],
    )
    # Note: non_coding_transcript_exon_variant ranks ahead of intron_variant in raw
    # Ensembl severity, but belongs to the retained_intron isoform (ENST00000100002.1),
    # whereas ENST00000100001.1 is protein_coding and MANE_Select.
    ann = ",".join([
        _make_ann_entry(
            "<INS:ME:ALU>", "non_coding_transcript_exon_variant", "MODIFIER",
            "GENE_POS", "ENSG00000100001", "ENST00000100002.1", "retained_intron", "1/1",
        ),
        _make_ann_entry(
            "<INS:ME:ALU>", "intron_variant", "MODIFIER",
            "GENE_POS", "ENSG00000100001", "ENST00000100001.1", "protein_coding", "1/1",
        ),
    ])
    out = _run_canned_annotation([rec], {1400: ann}, gtf)[0]

    assert out.info["CSQ"] == "intron_variant"
    assert out.info["TX_ID"] == "ENST00000100002.1,ENST00000100001.1"
    assert out.info["TX_BIOTYPE"] == "retained_intron,protein_coding"
    assert out.info["TX_RANK"] == "1/1"
    assert out.info["INTRON_RANK"] == "1/1"
    assert "EXON_RANK" not in out.info
    # Because the MANE protein_coding transcript is prioritized at pos=1400,
    # SPLICE_DIST is 300 bp (to donor at 1100) rather than 0 (exonic in retained_intron).
    assert out.info["SPLICE_DIST"] == "300"
    assert out.info["SPLICE_SIDE"] == "donor"


def test_vcf_4_3_symbolic_insertion_end_equals_pos_and_no_empty_fields(tmp_path: Path) -> None:
    gtf = _write_synthetic_gtf(tmp_path / "test.gtf")
    rec = ga.VcfRecord(
        chrom="chr1", pos=1250, id="vcf_spec", ref="N", alt="<INS:ME:ALU>",
        qual=".", filter="PASS",
        info={"SVTYPE": "INS", "SVLEN": "288", "END": "999999", "MEIFAMILY": "ALU", "ORIENT": "+"},
        rest=["GT:GQ", "./.:."],
    )
    ann = _make_ann_entry(
        "<INS:ME:ALU>", "intron_variant", "MODIFIER",
        "GENE_POS", "ENSG00000100001", "ENST00000100001.1", "protein_coding", "1/1",
    )
    out = _run_canned_annotation([rec], {1250: ann}, gtf)[0]
    assert out.info["END"] == "1250"
    assert out.info["SVLEN"] == "288"
    line = out.to_line()
    info_str = line.split("\t")[7]
    for token in info_str.split(";"):
        if "=" in token:
            k, v = token.split("=", 1)
            assert k != "" and v != "" and v != "."


# ---------------------------------------------------------------------------
# 5. Sequence-Level Cryptic 5' Donor & 3' Acceptor Motif Scanner
# ---------------------------------------------------------------------------


def test_cryptic_splice_motif_scanner_scores_and_pseudoexon_pairing() -> None:
    """Verify 9-mer 5' donor scoring, 23-mer 3' acceptor + PPT scoring,
    antisense reverse-complementation, and internal 40-350 bp pseudo-exon pairing.
    """
    # Canonical 9-mer donor CAG|GTAAGT scores 1.0; non-GT scores 0.0
    assert ga.score_donor_9mer("CAGGTAAGT") == pytest.approx(1.0, abs=0.01)
    assert ga.score_donor_9mer("CAGGCAAGT") == 0.0

    # Canonical 23-mer acceptor with 17-nt poly-T tract + C + AG + GTA scores > 0.90
    strong_acc = "TTTTTTTTTTTTTTTTTCAGGTA"
    assert len(strong_acc) == 23
    assert ga.score_acceptor_23mer(strong_acc) >= 0.90
    assert ga.score_acceptor_23mer("AAAAAAAAAAAAAAAAAAAGGTA") < 0.50

    # Construct a synthetic transcribed RNA carrying 1 strong 3' acceptor followed
    # 80 bp downstream by 1 strong 5' donor (a complete 80-bp pseudo-exon).
    # Exon starts right after AG of strong_acc (3 exonic bases "GTA" already in strong_acc),
    # then 74 neutral "AC" bases + "CAG" (total 3 + 74 + 3 = 80 bp exon) + "GTAAGT" + flank.
    transcribed_pseudoexon = (
        "GCGCGCGCGC"
        + strong_acc
        + ("AC" * 37)
        + "CAGGTAAGT"
        + "GCGCGCGCGCGCGC"
    )
    scan_sense = ga.scan_cryptic_splice_motifs(transcribed_pseudoexon, orient_to_gene="sense")
    assert scan_sense.acceptor_count >= 1
    assert scan_sense.donor_count >= 1
    assert scan_sense.pseudoexon_pairs >= 1
    assert scan_sense.max_score >= 0.90

    # When the insertion is stored in element orientation and integrates ANTISENSE to the
    # host gene, reverse-complementing the element sequence recovers the identical RNA motifs.
    element_seq_for_antisense = ga.reverse_complement_dna(transcribed_pseudoexon)
    scan_antisense = ga.scan_cryptic_splice_motifs(element_seq_for_antisense, orient_to_gene="antisense")
    assert scan_antisense.acceptor_count == scan_sense.acceptor_count
    assert scan_antisense.donor_count == scan_sense.donor_count
    assert scan_antisense.pseudoexon_pairs == scan_sense.pseudoexon_pairs
    assert scan_antisense.max_score == pytest.approx(scan_sense.max_score, abs=1e-4)


def test_sequence_level_scanner_refines_exonization_risk_in_vcf(tmp_path: Path) -> None:
    """Verify that CONSENSUS_SEQ motif scanning refines ALU/SVA/LINE1 exonization risk
    and populates CRYPTIC_DONOR_COUNT, CRYPTIC_ACCEPTOR_COUNT, CRYPTIC_PSEUDOEXON_PAIRS,
    and CRYPTIC_SPLICE_MAX_SCORE.
    """
    deep_gtf = tmp_path / "deep_motif.gtf"
    deep_gtf.write_text(
        'chr1\thavana\tgene\t1000\t5000\t.\t+\t.\tgene_id "ENSG00000300001"; gene_name "DEEPG";\n'
        'chr1\thavana\texon\t1000\t1100\t.\t+\t.\tgene_id "ENSG00000300001"; transcript_id "ENST00000300001.1"; exon_number "1"; transcript_biotype "protein_coding"; tag "MANE_Select";\n'
        'chr1\thavana\texon\t4900\t5000\t.\t+\t.\tgene_id "ENSG00000300001"; transcript_id "ENST00000300001.1"; exon_number "2"; transcript_biotype "protein_coding"; tag "MANE_Select";\n'
    )
    ann_deep = _make_ann_entry(
        "<INS:ME:ALU>", "intron_variant", "MODIFIER",
        "DEEPG", "ENSG00000300001", "ENST00000300001.1", "protein_coding", "1/1",
    )

    # 1. Deep intronic antisense Alu (pos=2500, SPLICE_DIST=1400 > 500) WITH a strong
    #    assembled pseudo-exon cassette on the transcribed strand -> upgraded from MODERATE to HIGH
    transcribed_cassette = "GCGCGCGCGC" + "TTTTTTTTTTTTTTTTTCAGGTA" + ("AC" * 37) + "CAGGTAAGT" + "GCGCGCGCGC"
    alu_element_seq = ga.reverse_complement_dna(transcribed_cassette)
    rec_alu_upgraded = ga.VcfRecord(
        chrom="chr1", pos=2500, id="alu_seq_high", ref="N", alt="<INS:ME:ALU>",
        qual=".", filter="PASS",
        info={
            "SVTYPE": "INS", "END": "2500", "MEIFAMILY": "ALU",
            "ORIENT": "-", "CONSENSUS_SEQ": alu_element_seq,
        },
        rest=["GT:GQ", "./.:."],
    )

    # 2. Deep intronic antisense Alu whose sequence is mutated/devoid of any AG/GT splice
    #    dinucleotides -> downgraded to NONE
    devoid_seq = "AC" * 60
    rec_alu_downgraded = ga.VcfRecord(
        chrom="chr1", pos=2500, id="alu_seq_none", ref="N", alt="<INS:ME:ALU>",
        qual=".", filter="PASS",
        info={
            "SVTYPE": "INS", "END": "2500", "MEIFAMILY": "ALU",
            "ORIENT": "-", "CONSENSUS_SEQ": devoid_seq,
        },
        rest=["GT:GQ", "./.:."],
    )

    out_high = _run_canned_annotation([rec_alu_upgraded], {2500: ann_deep}, deep_gtf)[0]
    out_none = _run_canned_annotation([rec_alu_downgraded], {2500: ann_deep}, deep_gtf)[0]

    assert out_high.info["ALU_EXONIZATION_RISK"] == "HIGH"
    assert out_high.info["MEI_EXONIZATION_RISK"] == "HIGH"
    assert int(out_high.info["CRYPTIC_ACCEPTOR_COUNT"] or "0") >= 1
    assert int(out_high.info["CRYPTIC_DONOR_COUNT"] or "0") >= 1
    assert int(out_high.info["CRYPTIC_PSEUDOEXON_PAIRS"] or "0") >= 1
    assert float(out_high.info["CRYPTIC_SPLICE_MAX_SCORE"] or "0") >= 0.90

    assert out_none.info["ALU_EXONIZATION_RISK"] == "NONE"
    assert out_none.info["MEI_EXONIZATION_RISK"] == "NONE"
    assert out_none.info["CRYPTIC_ACCEPTOR_COUNT"] == "0"
    assert out_none.info["CRYPTIC_DONOR_COUNT"] == "0"
    assert out_none.info["CRYPTIC_PSEUDOEXON_PAIRS"] == "0"


def test_panel_consensus_span_slicing_fallback_when_consensus_seq_absent(tmp_path: Path) -> None:
    """When a call lacks CONSENSUS_SEQ, verify that (MEISUBFAMILY, MEI_5P, MEI_3P)
    is sliced from the Dfam/Repbase consensus panel and scanned for cryptic splice motifs.
    """
    deep_gtf = tmp_path / "deep_panel.gtf"
    deep_gtf.write_text(
        'chr1\thavana\tgene\t1000\t5000\t.\t+\t.\tgene_id "ENSG00000300001"; gene_name "DEEPG";\n'
        'chr1\thavana\texon\t1000\t1100\t.\t+\t.\tgene_id "ENSG00000300001"; transcript_id "ENST00000300001.1"; exon_number "1"; transcript_biotype "protein_coding"; tag "MANE_Select";\n'
        'chr1\thavana\texon\t4900\t5000\t.\t+\t.\tgene_id "ENSG00000300001"; transcript_id "ENST00000300001.1"; exon_number "2"; transcript_biotype "protein_coding"; tag "MANE_Select";\n'
    )
    ann_deep = _make_ann_entry(
        "<INS:ME:ALU>", "intron_variant", "MODIFIER",
        "DEEPG", "ENSG00000300001", "ENST00000300001.1", "protein_coding", "1/1",
    )

    # Build a synthetic full-consensus panel entry (300 bp):
    # Bases 1..150 are neutral "AC"*75 (no GT/AG splice sites),
    # Bases 151..300 contain a strong antisense-Alu pseudo-exon cassette.
    transcribed_cassette = (
        "GCGCGCGCGC"
        + "TTTTTTTTTTTTTTTTTCAGGTA"
        + ("AC" * 37)
        + "CAGGTAAGT"
        + ("GC" * 17)
    )
    assert len(transcribed_cassette) == 150
    active_element_half = ga.reverse_complement_dna(transcribed_cassette)
    neutral_half = "AC" * 75
    full_alu_consensus = neutral_half + active_element_half
    assert len(full_alu_consensus) == 300

    panel_fa = tmp_path / "mei_full_canonical.panel.fa"
    panel_fa.write_text(
        f">AluY_full#SINE/Alu|source=AluY\n{full_alu_consensus}\n",
        encoding="utf-8",
    )

    # Call 1: Full-length AluY_short_#SINE/Alu (MEI_5P=1, MEI_3P=300) without CONSENSUS_SEQ
    # Slices bases 1..300 -> includes active cassette -> upgrades deep intronic to HIGH
    rec_full_span = ga.VcfRecord(
        chrom="chr1", pos=2500, id="alu_panel_full", ref="N", alt="<INS:ME:ALU>",
        qual=".", filter="PASS",
        info={
            "SVTYPE": "INS", "END": "2500", "MEIFAMILY": "ALU",
            "MEISUBFAMILY": "AluY_short_#SINE/Alu",
            "MEI_5P": "1", "MEI_3P": "300", "ORIENT": "-",
        },
        rest=["GT:GQ", "./.:."],
    )

    # Call 2: Truncated AluY_short_#SINE/Alu (MEI_5P=1, MEI_3P=140) without CONSENSUS_SEQ
    # Slices only bases 1..140 (the neutral half, missing positions 151..300) -> 0 acceptors -> NONE
    rec_trunc_span = ga.VcfRecord(
        chrom="chr1", pos=2500, id="alu_panel_trunc", ref="N", alt="<INS:ME:ALU>",
        qual=".", filter="PASS",
        info={
            "SVTYPE": "INS", "END": "2500", "MEIFAMILY": "ALU",
            "MEISUBFAMILY": "AluY_short_#SINE/Alu",
            "MEI_5P": "1", "MEI_3P": "140", "ORIENT": "-",
        },
        rest=["GT:GQ", "./.:."],
    )

    headers = [
        "##fileformat=VCFv4.3",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
    ]
    out_lines = list(headers)
    for r in (rec_full_span, rec_trunc_span):
        ic = dict(r.info)
        ic["ANN"] = ann_deep
        out_lines.append(
            ga.VcfRecord(r.chrom, r.pos, r.id, r.ref, r.alt, r.qual, r.filter, ic, list(r.rest)).to_line()
        )
    stdout_text = "\n".join(out_lines) + "\n"

    def fake_run(cmd, **kw):
        return subprocess.CompletedProcess(cmd, 0, stdout=stdout_text, stderr="")

    annotated, _ = ga.annotate_records(
        [rec_full_span, rec_trunc_span],
        headers,
        "GRCh38.99",
        gene_gtf=deep_gtf,
        panel_fasta=panel_fa,
        run=fake_run,
    )
    full_out, trunc_out = annotated[0], annotated[1]

    assert full_out.info["ALU_EXONIZATION_RISK"] == "HIGH"
    assert full_out.info["MEI_EXONIZATION_RISK"] == "HIGH"
    assert int(full_out.info["CRYPTIC_ACCEPTOR_COUNT"] or "0") >= 1
    assert int(full_out.info["CRYPTIC_DONOR_COUNT"] or "0") >= 1
    assert int(full_out.info["CRYPTIC_PSEUDOEXON_PAIRS"] or "0") >= 1

    assert trunc_out.info["ALU_EXONIZATION_RISK"] == "NONE"
    assert trunc_out.info["MEI_EXONIZATION_RISK"] == "NONE"
    assert trunc_out.info["CRYPTIC_ACCEPTOR_COUNT"] == "0"
    assert trunc_out.info["CRYPTIC_DONOR_COUNT"] == "0"
    assert trunc_out.info["CRYPTIC_PSEUDOEXON_PAIRS"] == "0"


# ---------------------------------------------------------------------------
# 6. Twin-Priming 5' Inversion Detection & Chimeric Junction Splice Scanning
# ---------------------------------------------------------------------------


def test_twin_priming_5p_inversion_chimeric_junction_donor_and_risk_upgrade(tmp_path: Path) -> None:
    """Verify twin-priming 5' inversion detection and chimeric junction scanning:
    1. Two segments that lack any 5' donor motif individually form a canonical
       9-mer 5' donor (CAG|GTAAGT) directly across the RC(5'_inv) + 3'_noninv
       inversion junction.
    2. An antisense LINE-1 insertion carrying a twin-priming 5' inversion flips
       its 5' segment into sense orientation on the pre-mRNA and creates chimeric
       junction splice sites, upgrading L1_TRUNCATION_RISK and MEI_EXONIZATION_RISK
       from NONE to HIGH and setting TWIN_PRIMING_INV, INV_JUNCTION_SPLICE_COUNT,
       and INV_JUNCTION_SPLICE_MAX_SCORE.
    """
    # Construct two segments devoid of GT dinucleotides individually:
    # We want RC(inv_5p_seg) to end with "CAG" -> inv_5p_seg must start with RC("CAG") = "CTG".
    # And we want noninv_3p_seg to start with "GTAAGT" + neutral bases -> chimeric junction is "CAG|GTAAGT".
    inv_5p_seg = "CTG" + ("AC" * 25)      # RC(inv_5p_seg) = ("GT"*25) ... wait: RC("AC") is "GT"! Use "AA" or "CC" instead!
    inv_5p_seg = "CTG" + ("CC" * 25)      # RC is ("GG" * 25) + "CAG" (0 GT, 0 AG)
    noninv_3p_seg = "GTAAGT" + ("CC" * 25)

    # Confirm inv_5p_seg alone has 0 donors, and noninv_3p_seg alone has 0 high-scoring donors (since no exonic CAG precedes GTAAGT)
    assert ga.scan_cryptic_splice_motifs(ga.reverse_complement_dna(inv_5p_seg), orient_to_gene="sense").donor_count == 0
    assert ga.scan_cryptic_splice_motifs(noninv_3p_seg, orient_to_gene="sense").donor_count == 0

    info_sense_inv = {
        "MEIFAMILY": "LINE1",
        "TWIN_PRIMING_INV": None,
        "INV_5P_SEQ": inv_5p_seg,
        "NONINV_3P_SEQ": noninv_3p_seg,
    }
    junc_res = ga.scan_twin_priming_inversion_junction(info_sense_inv, family="LINE1", orient_to_gene="sense")
    assert junc_res.is_inverted is True
    assert junc_res.junction_donor_count == 1
    assert junc_res.junction_total_count == 1
    assert junc_res.junction_max_score == pytest.approx(1.0, abs=0.01)

    # Now test full VCF annotation of an antisense LINE-1 insertion with twin-priming 5' inversion
    # (via ORIENT="-+" and panel coordinates) vs a normal uninverted antisense LINE-1 insertion.
    gtf = _write_synthetic_gtf(tmp_path / "test_inv.gtf")
    panel_fa = tmp_path / "l1_sva_inv_panel.fa"
    # Build a 200-bp panel sequence where bases 1..100 (5' segment) have a strong sense donor
    # ("CAGGTAAGT") that becomes active on the host gene's + strand when the 3' end integrates on "-"
    # and the 5' segment is inverted by twin priming!
    seg_5p_with_sense_cassette = (
        "CCCCCCCCCC"
        + "TTTTTTTTTTTTTTTTTCAGGTA"
        + ("AC" * 25)
        + "CAGGTAAGT"
        + "CCCCCCCC"
    )
    assert len(seg_5p_with_sense_cassette) == 100
    seg_3p_neutral = "CC" * 50
    panel_fa.write_text(
        f">L1HS_full#LINE/L1|source=L1HS\n{seg_5p_with_sense_cassette + seg_3p_neutral}\n"
        f">SVA_F_full#Retroposon/SVA|source=SVA_F\n{inv_5p_seg + noninv_3p_seg}\n",
        encoding="utf-8",
    )

    # Record 1: Antisense L1 (3' segment on '-', gene on '+') with twin-priming 5' inversion
    rec_l1_twin = ga.VcfRecord(
        chrom="chr1", pos=1500, id="l1_twin_as", ref="N", alt="<INS:ME:LINE1>",
        qual=".", filter="PASS",
        info={
            "SVTYPE": "INS", "END": "1500", "MEIFAMILY": "LINE1",
            "MEISUBFAMILY": "L1HS_5end#LINE/L1",
            "ORIENT": "+-",  # 5' segment +, 3' poly(A) segment - (antisense to + strand gene)
            "INV_5P_START": "1", "INV_5P_END": "100",
            "NONINV_5P_START": "101", "NONINV_3P_END": "200",
        },
        rest=["GT:GQ", "./.:."],
    )
    # Record 2: Sense SVA on + strand gene with twin-priming 5' inversion forming a chimeric junction donor
    rec_sva_twin = ga.VcfRecord(
        chrom="chr1", pos=1500, id="sva_twin_sense", ref="N", alt="<INS:ME:SVA>",
        qual=".", filter="PASS",
        info={
            "SVTYPE": "INS", "END": "1500", "MEIFAMILY": "SVA",
            "MEISUBFAMILY": "SVA_F#Retroposon/SVA",
            "ORIENT": "-+",  # 5' segment -, 3' poly(A) segment + (sense to + strand gene)
            "INV_5P_START": "1", "INV_5P_END": "53",
            "NONINV_5P_START": "54", "NONINV_3P_END": "109",
        },
        rest=["GT:GQ", "./.:."],
    )

    headers = [
        "##fileformat=VCFv4.3",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
    ]
    ann_l1 = _make_ann_entry(
        "<INS:ME:LINE1>", "intron_variant", "MODIFIER",
        "GENE_POS", "ENSG00000100001", "ENST00000100001.1", "protein_coding", "1/1",
    )
    ann_sva = _make_ann_entry(
        "<INS:ME:SVA>", "intron_variant", "MODIFIER",
        "GENE_POS", "ENSG00000100001", "ENST00000100001.1", "protein_coding", "1/1",
    )
    out_lines = list(headers)
    for r, a in ((rec_l1_twin, ann_l1), (rec_sva_twin, ann_sva)):
        ic = dict(r.info)
        ic["ANN"] = a
        out_lines.append(
            ga.VcfRecord(r.chrom, r.pos, r.id, r.ref, r.alt, r.qual, r.filter, ic, list(r.rest)).to_line()
        )
    stdout_text = "\n".join(out_lines) + "\n"

    def fake_run(cmd, **kw):
        return subprocess.CompletedProcess(cmd, 0, stdout=stdout_text, stderr="")

    annotated, _ = ga.annotate_records(
        [rec_l1_twin, rec_sva_twin],
        headers,
        "GRCh38.99",
        gene_gtf=gtf,
        panel_fasta=panel_fa,
        run=fake_run,
    )
    l1_out, sva_out = annotated[0], annotated[1]

    assert l1_out.info["MEI_GENE_ORIENT"] == "antisense"
    assert "TWIN_PRIMING_INV" in l1_out.info and l1_out.info["TWIN_PRIMING_INV"] is None
    assert l1_out.info["L1_TRUNCATION_RISK"] == "HIGH"
    assert l1_out.info["MEI_EXONIZATION_RISK"] == "HIGH"
    assert int(l1_out.info["CRYPTIC_PSEUDOEXON_PAIRS"] or "0") >= 1

    assert sva_out.info["MEI_GENE_ORIENT"] == "sense"
    assert "TWIN_PRIMING_INV" in sva_out.info and sva_out.info["TWIN_PRIMING_INV"] is None
    assert sva_out.info["SVA_EXON_TRAPPING_RISK"] == "HIGH"
    assert sva_out.info["MEI_EXONIZATION_RISK"] == "HIGH"
    assert int(sva_out.info["INV_JUNCTION_SPLICE_COUNT"] or "0") >= 1
    assert float(sva_out.info["INV_JUNCTION_SPLICE_MAX_SCORE"] or "0") == pytest.approx(1.0, abs=0.01)



