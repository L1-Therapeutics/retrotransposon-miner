"""Tests for vcf_export, built against real rows from the README gold-tier example table."""

import shutil
import subprocess
from pathlib import Path

import pytest

from retro_miner.vcf_export import (
    build_vcf_record,
    export_vcf,
    export_vcf_from_tsv,
    VCF_HEADER_LINES,
    VCF_COLUMN_HEADER,
)

# Row 1: real SVA record from README's gold-tier chr22 example table.
SVA_ROW = {
    "chrom": "chr22",
    "consensus_insertion_breakpoint_pos": "49029650",
    "window_start": "49029645",
    "window_end": "49029656",
    "control_supporting_reads": "SR_L=0,SR_R=0,DPE_L=2,DPE_R=152,MEI_MAPPED=127,polyA_MAPPED=46,VNTR_MAPPED=0,polyA_side=L",
    "disease_supporting_reads": "SR_L=0,SR_R=0,DPE_L=8,DPE_R=373,MEI_MAPPED=306,polyA_MAPPED=92,VNTR_MAPPED=1,polyA_side=L",
    "sample_status_label": "shared",
    "consensus_tsd_seq": "AAGAAAACTCCT",
    "consensus_poly_at_min_bp": "50",
    "consensus_mei_family": "SVA",
    "consensus_mei_subfamily": "SVA_D#Retroposon/SVA",
    "known_mei_polymorphism_id": "nssv14064350",
    "known_mei_polymorphism_source": "melt_1kg",
    "consensus_insertion_orientation": "-",
    "nested_in_same_MEI": "unnested",
    "consensus_insertion_mei_span_full": "1366",
    "consensus_insertion_mei_5p_coord_full": "1",
    "consensus_insertion_mei_3p_coord_full": "1366",
}

# Row 2: real LINE1 record, control_only, with a multi-source known-id field
# (contains a literal ';' in the source list -- this is exactly the kind of
# value that must be sanitized rather than break VCF syntax).
LINE1_ROW = {
    "chrom": "chr22",
    "consensus_insertion_breakpoint_pos": "19223382",
    "window_start": "19223373",
    "window_end": "19223390",
    "control_supporting_reads": "SR_L=0,SR_R=14,DPE_L=60,DPE_R=27,MEI_MAPPED=91,polyA_MAPPED=17,polyA_side=L",
    "disease_supporting_reads": "SR_L=0,SR_R=0,DPE_L=0,DPE_R=0,MEI_MAPPED=0,polyA_MAPPED=0",
    "sample_status_label": "control_only",
    "consensus_tsd_seq": "AAAAACCACCTATGCTGG",
    "consensus_poly_at_min_bp": "66",
    "consensus_mei_family": "LINE1",
    "consensus_mei_subfamily": "L1HS_5end#LINE/L1",
    "known_mei_polymorphism_id": "g1k:nssv14064681;lr:chr22-19600083-INS->s899391<s914453>s899392-6059",
    "known_mei_polymorphism_source": "melt_1kg,long_read_1kg_ont_vienna",
    "consensus_insertion_orientation": "+",
    "nested_in_same_MEI": "unnested",
    "consensus_insertion_mei_span_full": "6018",
    "consensus_insertion_mei_5p_coord_full": "1",
    "consensus_insertion_mei_3p_coord_full": "6018",
}

# Row 3: real ALU row with several blank optional fields (no known_mei_polymorphism_id/source).
BLANK_OPTIONAL_ROW = {
    "chrom": "chr22",
    "consensus_insertion_breakpoint_pos": "31355872",
    "window_start": "31355858",
    "window_end": "31355887",
    "control_supporting_reads": "SR_L=12,SR_R=0,DPE_L=7,DPE_R=148,MEI_MAPPED=108,polyA_MAPPED=6",
    "disease_supporting_reads": "SR_L=15,SR_R=0,DPE_L=19,DPE_R=273,MEI_MAPPED=188,polyA_MAPPED=16,polyA_side=R",
    "sample_status_label": "shared",
    "consensus_tsd_seq": "CCGCCTCGGCTTCCCAAAGTGCTGGGATTA",
    "consensus_poly_at_min_bp": "71",
    "consensus_mei_family": "ALU",
    "consensus_mei_subfamily": "AluY_short_#SINE/Alu",
    "known_mei_polymorphism_id": "",
    "known_mei_polymorphism_source": "",
    "consensus_insertion_orientation": "-",
    "nested_in_same_MEI": "nested",
    "consensus_insertion_mei_span_full": "281",
    "consensus_insertion_mei_5p_coord_full": "1",
    "consensus_insertion_mei_3p_coord_full": "281",
}


class TestBuildVcfRecord:
    def test_field_count_matches_column_header(self):
        rec = build_vcf_record(SVA_ROW)
        assert len(rec.split("\t")) == len(VCF_COLUMN_HEADER.split("\t"))

    def test_chrom_pos_from_row(self):
        rec = build_vcf_record(SVA_ROW)
        fields = rec.split("\t")
        assert fields[0] == "chr22"
        assert fields[1] == "49029650"

    def test_alt_allele_maps_known_family(self):
        assert build_vcf_record(SVA_ROW).split("\t")[4] == "<INS:ME:SVA>"
        assert build_vcf_record(LINE1_ROW).split("\t")[4] == "<INS:ME:LINE1>"
        assert build_vcf_record(BLANK_OPTIONAL_ROW).split("\t")[4] == "<INS:ME:ALU>"

    def test_unknown_family_falls_back_to_generic_symbolic_allele(self):
        row = dict(SVA_ROW, consensus_mei_family="UNKNOWN_NEW_FAMILY")
        assert build_vcf_record(row).split("\t")[4] == "<INS:ME>"

    def test_id_field_uses_known_polymorphism_id_when_present(self):
        assert build_vcf_record(SVA_ROW).split("\t")[2] == "nssv14064350"

    def test_id_field_is_dot_when_absent(self):
        assert build_vcf_record(BLANK_OPTIONAL_ROW).split("\t")[2] == "."

    def test_genotype_is_always_blank_missing_not_fabricated(self):
        for row in (SVA_ROW, LINE1_ROW, BLANK_OPTIONAL_ROW):
            fields = build_vcf_record(row).split("\t")
            fmt, sample_val = fields[8], fields[9]
            assert fmt == "GT:GQ"
            assert sample_val == "./.:."

    def test_qual_and_filter_are_missing_not_fabricated(self):
        fields = build_vcf_record(SVA_ROW).split("\t")
        assert fields[5] == "."  # QUAL
        assert fields[6] == "."  # FILTER

    def test_semicolon_in_known_id_is_sanitized_not_left_raw(self):
        rec = build_vcf_record(LINE1_ROW)
        info_field = rec.split("\t")[7]
        # A raw, un-sanitized ';' here would corrupt VCF INFO key=value parsing.
        assert "KNOWN_ID=g1k:nssv14064681_lr:chr22-19600083-INS->s899391<s914453>s899392-6059" in info_field

    def test_comma_separated_evidence_string_is_sanitized(self):
        # Both ',' and internal '=' must be escaped: VCF INFO reserves '='
        # for the top-level KEY=VALUE split, so an inner "SR_L=0" cannot be
        # left with a raw '=' either, or a parser would split there too.
        rec = build_vcf_record(SVA_ROW)
        info_field = rec.split("\t")[7]
        assert "CTRL_SUPPORT=SR_L:0|SR_R:0|DPE_L:2|DPE_R:152|MEI_MAPPED:127|polyA_MAPPED:46|VNTR_MAPPED:0|polyA_side:L" in info_field

    def test_blank_optional_fields_are_omitted_from_info_not_padded(self):
        rec = build_vcf_record(BLANK_OPTIONAL_ROW)
        info_field = rec.split("\t")[7]
        assert "KNOWN_ID=" not in info_field
        assert "KNOWN_SRC=" not in info_field

    def test_no_info_field_ever_contains_a_raw_comma(self):
        # A raw unescaped comma inside a Number=1 INFO value is invalid VCF.
        for row in (SVA_ROW, LINE1_ROW, BLANK_OPTIONAL_ROW):
            info_field = build_vcf_record(row).split("\t")[7]
            for kv in info_field.split(";"):
                if "=" in kv:
                    _, value = kv.split("=", 1)
                    assert "," not in value

    def test_missing_breakpoint_pos_becomes_dot_not_zero(self):
        row = dict(SVA_ROW, consensus_insertion_breakpoint_pos="")
        assert build_vcf_record(row).split("\t")[1] == "."


class TestExportVcf:
    def test_writes_header_and_one_line_per_row(self, tmp_path: Path):
        out_path = tmp_path / "out.vcf"
        n = export_vcf([SVA_ROW, LINE1_ROW, BLANK_OPTIONAL_ROW], out_path)
        assert n == 3
        text = out_path.read_text()
        lines = text.splitlines()
        assert lines[: len(VCF_HEADER_LINES)] == VCF_HEADER_LINES
        assert lines[len(VCF_HEADER_LINES)] == VCF_COLUMN_HEADER.format(sample="SAMPLE")
        data_lines = lines[len(VCF_HEADER_LINES) + 1 :]
        assert len(data_lines) == 3

    def test_custom_sample_name_in_column_header(self, tmp_path: Path):
        out_path = tmp_path / "out.vcf"
        export_vcf([SVA_ROW], out_path, sample_name="tumor_normal_pool")
        text = out_path.read_text()
        assert "\ttumor_normal_pool" in text

    def test_empty_rows_still_produces_valid_header_only_vcf(self, tmp_path: Path):
        out_path = tmp_path / "out.vcf"
        n = export_vcf([], out_path)
        assert n == 0
        text = out_path.read_text()
        assert text.splitlines()[-1] == VCF_COLUMN_HEADER.format(sample="SAMPLE")


class TestExportVcfFromTsv:
    def test_round_trips_a_real_tsv(self, tmp_path: Path):
        tsv_path = tmp_path / "candidate_loci.annotated.tsv"
        fieldnames = list(SVA_ROW.keys())
        with open(tsv_path, "w") as fh:
            fh.write("\t".join(fieldnames) + "\n")
            fh.write("\t".join(SVA_ROW[k] for k in fieldnames) + "\n")
            fh.write("\t".join(LINE1_ROW[k] for k in fieldnames) + "\n")

        out_path = tmp_path / "out.vcf"
        n = export_vcf_from_tsv(tsv_path, out_path)
        assert n == 2
        text = out_path.read_text()
        assert "chr22\t49029650" in text
        assert "chr22\t19223382" in text

    def test_missing_input_file_raises(self, tmp_path: Path):
        with pytest.raises(FileNotFoundError):
            export_vcf_from_tsv(tmp_path / "does_not_exist.tsv", tmp_path / "out.vcf")


class TestCoordinateSorting:
    """VCF requires contig-then-position order; the upstream TSV is sorted by
    enrichment_ratio instead, which makes the file unindexable."""

    def _rows_out_of_order(self):
        return [
            dict(SVA_ROW, chrom="chr22", consensus_insertion_breakpoint_pos="49029650"),
            dict(SVA_ROW, chrom="chr22", consensus_insertion_breakpoint_pos="31355872"),
            dict(SVA_ROW, chrom="chr2", consensus_insertion_breakpoint_pos="17567662"),
            dict(SVA_ROW, chrom="chr1", consensus_insertion_breakpoint_pos="900000"),
        ]

    def test_records_are_coordinate_sorted_by_default(self, tmp_path: Path):
        out_path = tmp_path / "out.vcf"
        export_vcf(self._rows_out_of_order(), out_path)
        data = [
            line.split("\t")[:2]
            for line in out_path.read_text().splitlines()
            if not line.startswith("#")
        ]
        assert data == [
            ["chr1", "900000"],
            ["chr2", "17567662"],
            ["chr22", "31355872"],
            ["chr22", "49029650"],
        ]

    def test_contig_order_follows_header_not_lexicographic(self, tmp_path: Path):
        # Lexicographically "chr22" < "chr2"[sic] is false, but naive string
        # sorting would place chr10/chr2/chr22 wrongly relative to chr3.
        rows = [
            dict(SVA_ROW, chrom="chr3", consensus_insertion_breakpoint_pos="100"),
            dict(SVA_ROW, chrom="chr10", consensus_insertion_breakpoint_pos="100"),
            dict(SVA_ROW, chrom="chr2", consensus_insertion_breakpoint_pos="100"),
        ]
        out_path = tmp_path / "out.vcf"
        export_vcf(rows, out_path)
        chroms = [
            line.split("\t")[0]
            for line in out_path.read_text().splitlines()
            if not line.startswith("#")
        ]
        assert chroms == ["chr2", "chr3", "chr10"]

    def test_unknown_contig_sorts_last_and_is_not_dropped(self, tmp_path: Path):
        rows = [
            dict(SVA_ROW, chrom="chrUn_decoy1", consensus_insertion_breakpoint_pos="500"),
            dict(SVA_ROW, chrom="chr1", consensus_insertion_breakpoint_pos="500"),
        ]
        out_path = tmp_path / "out.vcf"
        n = export_vcf(rows, out_path)
        assert n == 2  # nothing dropped
        chroms = [
            line.split("\t")[0]
            for line in out_path.read_text().splitlines()
            if not line.startswith("#")
        ]
        assert chroms == ["chr1", "chrUn_decoy1"]

    def test_sort_false_preserves_input_order(self, tmp_path: Path):
        out_path = tmp_path / "out.vcf"
        export_vcf(self._rows_out_of_order(), out_path, sort=False)
        data = [
            line.split("\t")[:2]
            for line in out_path.read_text().splitlines()
            if not line.startswith("#")
        ]
        assert data[0] == ["chr22", "49029650"]

    def test_output_is_indexable_by_bcftools(self, tmp_path: Path):
        """The actual downstream requirement: bcftools must be able to index it."""
        bcftools = shutil.which("bcftools")
        bgzip = shutil.which("bgzip")
        if not bcftools or not bgzip:
            pytest.skip("bcftools/bgzip not installed")

        vcf_path = tmp_path / "out.vcf"
        export_vcf(self._rows_out_of_order(), vcf_path)

        gz_path = tmp_path / "out.vcf.gz"
        with open(gz_path, "wb") as fh:
            subprocess.run([bgzip, "-c", str(vcf_path)], stdout=fh, check=True)

        result = subprocess.run(
            [bcftools, "index", str(gz_path)], capture_output=True, text=True
        )
        assert result.returncode == 0, f"bcftools index failed: {result.stderr}"


class TestPysamRoundTrip:
    """Parse our own output with a real VCF library, not just string checks."""

    def test_pysam_parses_header_and_records_cleanly(self, tmp_path: Path):
        pysam = pytest.importorskip("pysam")
        out_path = tmp_path / "out.vcf"
        export_vcf([SVA_ROW, LINE1_ROW, BLANK_OPTIONAL_ROW], out_path)

        vf = pysam.VariantFile(str(out_path))
        assert list(vf.header.samples) == ["SAMPLE"]
        assert vf.header.formats["GQ"].type == "Integer"
        assert "chr22" in vf.header.contigs

        records = list(vf)
        assert len(records) == 3
        # Output is coordinate-sorted, so look the record up by position
        # rather than assuming the input order is preserved.
        sva = next(r for r in records if r.pos == 49029650)
        assert sva.chrom == "chr22"
        assert sva.id == "nssv14064350"
        assert sva.alts == ("<INS:ME:SVA>",)
        assert sva.samples[0]["GT"] == (None, None)  # ./. -- blank, not fabricated
        assert sva.samples[0]["GQ"] is None  # . -- blank, not fabricated
        assert [r.pos for r in records] == sorted(r.pos for r in records)
