import pytest
import tempfile
from pathlib import Path
from scripts.compare_callset_parity import compare_parity

def create_mock_vcf(path: Path, records: list):
    with open(path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for r in records:
            f.write(f"{r[0]}\t{r[1]}\t.\t{r[2]}\t{r[3]}\t60\tPASS\t.\n")

def test_exact_parity_match():
    with tempfile.TemporaryDirectory() as tmpdir:
        b_path = Path(tmpdir) / "base.vcf"
        c_path = Path(tmpdir) / "cand.vcf"
        recs = [("chr1", "1000", "A", "<INS:MEI>"), ("chr2", "5000", "C", "<INS:MEI>")]
        create_mock_vcf(b_path, recs)
        create_mock_vcf(c_path, recs)
        
        res = compare_parity(b_path, c_path)
        assert res["exact_match"] is True
        assert res["concordance_pct"] == 100.0

def test_divergence_detection():
    with tempfile.TemporaryDirectory() as tmpdir:
        b_path = Path(tmpdir) / "base.vcf"
        c_path = Path(tmpdir) / "cand.vcf"
        create_mock_vcf(b_path, [("chr1", "1000", "A", "<INS:MEI>")])
        create_mock_vcf(c_path, [("chr1", "1000", "A", "<INS:MEI>"), ("chr3", "200", "G", "<INS:MEI>")])
        
        res = compare_parity(b_path, c_path)
        assert res["exact_match"] is False
        assert res["extra_in_candidate"] == 1
