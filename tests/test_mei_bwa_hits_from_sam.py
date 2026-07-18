"""Parse SAM hits from bwa mem MEI remapping."""

from __future__ import annotations

from retro_miner.mei_support import _best_hits_from_sam, _cigar_alignment_spans


def test_cigar_spans():
    q, r, aln = _cigar_alignment_spans("5S20M3S")
    assert q == 20
    assert r == 20
    assert aln == 20


def test_best_hits_from_sam_picks_mapped():
    sam = "\n".join(
        [
            "@SQ\tSN:AluY#SINE/Alu\tLN:311",
            "q1\t0\tAluY#SINE/Alu\t1\t0\t20M\t*\t0\t0\tGGCCGGGCGCGGTGGCTCAC\t*",
            "q1\t256\tAluY#SINE/Alu\t10\t0\t10M\t*\t0\t0\tGGCCGGGCGC\t*",
            "q2\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t*",
        ]
    )
    hits = _best_hits_from_sam(sam)
    assert list(hits["qname"]) == ["q1"]
    row = hits.iloc[0]
    assert row.target.startswith("AluY")
    assert int(row.target_start) == 1
    assert int(row.alnlen) == 20
    assert row.family == "ALU"
    assert float(row.qcov) == 1.0


def test_best_hits_skips_supplementary_hardclip_fake_qcov():
    # Primary: 13 bp tip of a 55 bp clip (true qcov ~0.24).
    # Supplementary hard-clip would look like qcov=1.0 if SEQ is only 13 bp.
    sam = "\n".join(
        [
            "@SQ\tSN:L1HS#LINE/L1\tLN:6000",
            "q1\t16\tL1HS#LINE/L1\t1360\t0\t38S13M4S\t*\t0\t0\t"
            + ("A" * 38)
            + ("G" * 13)
            + ("T" * 4)
            + "\t*\tNM:i:0\tAS:i:13",
            "q1\t2064\tL1HS#LINE/L1\t1360\t0\t26H13M16H\t*\t0\t0\t"
            + ("G" * 13)
            + "\t*\tNM:i:0\tAS:i:13",
        ]
    )
    hits = _best_hits_from_sam(sam)
    assert len(hits) == 1
    row = hits.iloc[0]
    assert int(row.alnlen) == 13
    assert abs(float(row.qcov) - (13 / 55)) < 1e-6
