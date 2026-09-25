"""Tests for vcf_export, built against real rows from the README gold-tier example table."""

import csv
import shutil
import subprocess
from pathlib import Path

import pytest

from retro_miner.vcf_export import (
    build_vcf_record,
    export_vcf,
    export_vcf_from_tsv,
    format_sigfigs,
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

    def test_id_is_l1tx_chrom_pos_family(self):
        assert build_vcf_record(SVA_ROW).split("\t")[2] == "L1TX-chr22-49029650-SVA"
        assert build_vcf_record(LINE1_ROW).split("\t")[2] == "L1TX-chr22-19223382-LINE1"
        assert build_vcf_record(BLANK_OPTIONAL_ROW).split("\t")[2] == "L1TX-chr22-31355872-ALU"

    def test_id_is_dot_when_breakpoint_is_missing(self):
        row = dict(SVA_ROW, consensus_insertion_breakpoint_pos="")
        assert build_vcf_record(row).split("\t")[2] == "."

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

    def test_catalog_ids_are_split_into_their_own_info_fields(self):
        info_field = build_vcf_record(LINE1_ROW).split("\t")[7]
        assert "NSSV=" not in info_field
        assert "G1K=nssv14064681" in info_field
        assert "LR=chr22-19600083-INS->s899391<s914453>s899392-6059" in info_field
        assert "KNOWN_ID=" not in info_field
        for kv in info_field.split(";"):
            if "=" in kv:
                _, value = kv.split("=", 1)
                assert ";" not in value

        sva_info = build_vcf_record(SVA_ROW).split("\t")[7]
        assert "NSSV=" not in sva_info
        assert "G1K=nssv14064350" in sva_info
        assert "LR=" not in sva_info
        sva_keys = [kv.split("=", 1)[0] for kv in sva_info.split(";") if "=" in kv]
        assert "CTRL_SUPPORT" in sva_keys
        assert "DISEASE_SUPPORT" in sva_keys
        assert "SUPPORT" not in sva_keys

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
        assert "G1K=" not in info_field
        assert "LR=" not in info_field
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
        assert "SVLEN" not in records[0].info
        assert records[0].info["SVTYPE"] == "INS"
        data_lines = [ln for ln in out_path.read_text().splitlines() if not ln.startswith("#")]
        for line in data_lines:
            pos = line.split("\t")[1]
            assert f"END={pos}" in line.split("\t")[7]
        # Output is coordinate-sorted, so look the record up by position
        # rather than assuming the input order is preserved.
        sva = next(r for r in records if r.pos == 49029650)
        assert sva.chrom == "chr22"
        assert sva.id == "L1TX-chr22-49029650-SVA"
        assert sva.alts == ("<INS:ME:SVA>",)
        assert sva.samples[0]["GT"] == (None, None)  # ./. -- blank, not fabricated
        assert sva.samples[0]["GQ"] is None  # . -- blank, not fabricated
        assert [r.pos for r in records] == sorted(r.pos for r in records)


FIXTURE_DIR = Path(__file__).resolve().parent / "fixtures" / "vcf"
CLASSIFIER_SLICE = FIXTURE_DIR / "hg03086_classifier_slice.tsv"
GOLD_SLICE = FIXTURE_DIR / "hg03086_gold_slice.tsv"


def test_gold_review_export_has_no_rank_pct_without_classifier(tmp_path: Path):
    """Real gold tables have no classifier_rank (and no gold_rank); emit neither pct."""
    pysam = pytest.importorskip("pysam")
    out_path = tmp_path / "gold.vcf"
    n = export_vcf_from_tsv(GOLD_SLICE, out_path, sample_name="HG03086")
    assert n == 4
    text = out_path.read_text()
    assert "L1TXRANKPCT=" not in text
    assert "L1TXGOLDRANKPCT=" not in text
    assert "L1TX-chr4-66240893-ALU" in text
    assert "L1TX-chr4-102422592-LINE1" in text
    assert "L1TX-chr10-3041566-SVA" in text
    records = list(pysam.VariantFile(str(out_path)))
    by_id = {rec.id: rec for rec in records}
    alu = by_id["L1TX-chr4-66240893-ALU"]
    assert alu.alts == ("<INS:ME:ALU>",)
    assert "NSSV" not in alu.info
    assert alu.info["G1K"] == "nssv14044437"
    assert "CTRL_SUPPORT" not in alu.info
    assert "DISEASE_SUPPORT" not in alu.info
    assert "SUPPORT" in alu.info
    line1 = by_id["L1TX-chr4-102422592-LINE1"]
    assert line1.alts == ("<INS:ME:LINE1>",)
    assert "NSSV" not in line1.info
    assert line1.info["G1K"] == "nssv14080750"
    assert str(line1.info["LR"]).startswith("chr4-105736355-INS->")
    assert "SUPPORT" in line1.info
    sva = by_id["L1TX-chr10-3041566-SVA"]
    assert sva.alts == ("<INS:ME:SVA>",)
    assert "NSSV" not in sva.info
    assert sva.info["G1K"] == "nssv14066958"
    assert "SUPPORT" in sva.info


class TestHg03086Tables:
    """Rows sliced from the HG03086 gold review and classifier ranking."""

    def test_classifier_without_breakpoint_table_is_rejected(self, tmp_path: Path):
        with pytest.raises(ValueError, match="consensus_insertion_breakpoint_pos"):
            export_vcf_from_tsv(CLASSIFIER_SLICE, tmp_path / "out.vcf")

    def test_high_score_classifier_rows_use_gold_breakpoints(self, tmp_path: Path):
        pysam = pytest.importorskip("pysam")
        out_path = tmp_path / "out.vcf"
        n = export_vcf_from_tsv(
            CLASSIFIER_SLICE,
            out_path,
            min_score=0.997,
            breakpoint_tsv=GOLD_SLICE,
        )
        assert n == 3

        gold = {
            (row["chrom"], row["window_start"], row["window_end"]): row
            for row in csv.DictReader(GOLD_SLICE.open(), delimiter="\t")
        }
        classifier_rows = list(csv.DictReader(CLASSIFIER_SLICE.open(), delimiter="\t"))
        kept = [row for row in classifier_rows if float(row["gold_score"]) >= 0.997]
        assert len(kept) == 3
        rank_n = max(int(row["classifier_rank"]) for row in classifier_rows)

        vf = pysam.VariantFile(str(out_path))
        assert list(vf.header.samples) == ["HG03086"]
        records = list(vf)
        for line in out_path.read_text().splitlines():
            if line.startswith("#"):
                continue
            pos = line.split("\t")[1]
            assert f"END={pos};" in line or line.split("\t")[7].endswith(f"END={pos}")
        assert len(records) == 3
        by_pos = {rec.pos: rec for rec in records}
        line_by_pos = {
            int(line.split("\t")[1]): line.split("\t")
            for line in out_path.read_text().splitlines()
            if not line.startswith("#")
        }
        family_token = {"Alu": "ALU", "L1": "LINE1", "SVA": "SVA"}
        saw_distinct_breakpoint = False
        for row in kept:
            locus = gold[(row["chrom"], row["window_start"], row["window_end"])]
            pos = int(locus["consensus_insertion_breakpoint_pos"])
            if pos != int(row["window_start"]):
                saw_distinct_breakpoint = True
            rec = by_pos[pos]
            fields = line_by_pos[pos]
            info = dict(part.split("=", 1) for part in fields[7].split(";") if "=" in part)
            assert rec.chrom == row["chrom"]
            assert rec.info["SVTYPE"] == "INS"
            assert "SVLEN" not in rec.info
            assert list(rec.filter) == ["PASS"]
            assert fields[2] == f"L1TX-{row['chrom']}-{pos}-{family_token[row['mei_family']]}"
            assert info["L1TXGOLDSCORE"] == format_sigfigs(row["gold_score"])
            pct = int(round(100.0 * (rank_n - int(row["classifier_rank"]) + 1) / rank_n))
            assert rec.info["L1TXRANKPCT"] == pct
            assert "L1TXGOLDRANKPCT" not in info
            assert "CLASSIFIERRANK" not in info
            assert "GOLDSCORE" not in info
            assert "GOLDRANK" not in info
            assert rec.info["INSERTIONSCORE"] == pytest.approx(float(locus["insertion_model_score"]))
            expected = {"Alu": "<INS:ME:ALU>", "L1": "<INS:ME:LINE1>", "SVA": "<INS:ME:SVA>"}
            assert rec.alts == (expected[row["mei_family"]],)
        line1 = line_by_pos[102422592]
        line1_info = dict(part.split("=", 1) for part in line1[7].split(";") if "=" in part)
        assert "NSSV" not in line1_info
        assert line1_info["G1K"] == "nssv14080750"
        assert "SUPPORT" in line1_info
        assert "CTRL_SUPPORT" not in line1_info
        assert "DISEASE_SUPPORT" not in line1_info
        assert line1_info["LR"].startswith("chr4-105736355-INS->")
        assert saw_distinct_breakpoint

    def test_classifier_export_inherits_reference_from_gold_table(self, tmp_path: Path):
        """Classifier file sits outside the run dir; ##reference comes from the gold table."""
        ranked = tmp_path / "classifier" / CLASSIFIER_SLICE.name
        run = tmp_path / "run"
        gold = run / GOLD_SLICE.name
        ranked.parent.mkdir()
        run.mkdir()
        ranked.write_text(CLASSIFIER_SLICE.read_text(), encoding="utf-8")
        gold.write_text(GOLD_SLICE.read_text(), encoding="utf-8")
        (run / "pipeline_params.env").write_text(
            "reference_build=hg38\n"
            "reference_fasta=/data/reference/hg38/Homo_sapiens_assembly38.fasta\n",
            encoding="utf-8",
        )
        out_path = tmp_path / "classifier.vcf"
        n = export_vcf_from_tsv(
            ranked,
            out_path,
            min_score=0.997,
            breakpoint_tsv=gold,
        )
        assert n == 3
        text = out_path.read_text(encoding="utf-8")
        assert "##reference=hg38\n" in text
        assert "##assembly=GRCh38\n" in text
        assert ",assembly=" not in text
        assert "L1TX-chr4-66240893-ALU" in text
        assert "L1TX-chr4-102422592-LINE1" in text
        assert "L1TX-chr10-3041566-SVA" in text

    def test_gold_review_slice_exports_without_classifier_columns(self, tmp_path: Path):
        pysam = pytest.importorskip("pysam")
        out_path = tmp_path / "gold.vcf"
        n = export_vcf_from_tsv(GOLD_SLICE, out_path, sample_name="HG03086")
        assert n == 4
        records = list(pysam.VariantFile(str(out_path)))
        assert len(records) == 4
        assert all("INSERTIONSCORE" in rec.info for rec in records)
        assert all("L1TXGOLDSCORE" not in rec.info for rec in records)
        assert all(rec.id.startswith("L1TX-") for rec in records)
        for line in out_path.read_text().splitlines():
            if line.startswith("#"):
                continue
            pos = line.split("\t")[1]
            assert f"END={pos}" in line.split("\t")[7]
            assert "SVLEN=" not in line

    def test_unmatched_classifier_window_raises(self, tmp_path: Path):
        bad = tmp_path / "bad.tsv"
        bad.write_text(
            CLASSIFIER_SLICE.read_text().splitlines()[0]
            + "\n"
            + "1\tHG03086\tchr1\t1\tAlu\tpositive\tTrue\tTrue\t0.999\t1\t10\n"
        )
        with pytest.raises(ValueError, match="did not match the breakpoint table"):
            export_vcf_from_tsv(bad, tmp_path / "out.vcf", breakpoint_tsv=GOLD_SLICE)
