import gzip

import numpy as np
import pandas as pd
import pytest

from carrier_logistic import fit, flag_exonic, load_exons, prepare, summarise

GTF = """\
chr1\tH\tgene\t1\t9000\t.\t+\t.\tgene_id "G1";
chr1\tH\texon\t1001\t1100\t.\t+\t.\tgene_id "G1";
chr1\tH\texon\t1051\t1200\t.\t+\t.\tgene_id "G1";
chr1\tH\texon\t5001\t5200\t.\t+\t.\tgene_id "G1";
"""


@pytest.fixture
def exons(tmp_path):
    p = tmp_path / "t.gtf.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(GTF)
    return load_exons(str(p))


def test_load_exons_merges_overlapping_entries(exons):
    assert exons["chr1"].tolist() == [[1000, 1200], [5000, 5200]]


def _copies(spans):
    return pd.DataFrame({"chrom": ["chr1"] * len(spans),
                         "start": [s for s, _ in spans],
                         "end": [e for _, e in spans]})


def test_flag_exonic_detects_overlap_and_misses(exons):
    df = _copies([(500, 600), (1150, 1250), (1000, 1200), (3000, 3100), (5100, 5300)])
    assert flag_exonic(df, exons).tolist() == [False, True, True, False, True]


def test_flag_exonic_on_a_chromosome_with_no_exons(exons):
    df = pd.DataFrame({"chrom": ["chrZ"], "start": [100], "end": [200]})
    assert flag_exonic(df, exons).tolist() == [False]


def test_an_element_abutting_an_exon_does_not_count_as_overlapping(exons):
    # exon is [1000, 1200); an element ending exactly at 1000 does not overlap
    assert flag_exonic(_copies([(900, 1000)]), exons).tolist() == [False]


def test_prepare_builds_carrier_flags_and_scaled_covariates():
    df = pd.DataFrame({"sense": [True, False], "n_hex_any": [0, 2], "n_hex_body": [0, 1],
                       "milli_div": [50, 100], "length_bp": [300, 280]})
    out = prepare(df, "AluY", np.array([False, True]))
    assert out.sense.tolist() == [1, 0]
    assert out.carrier_any.tolist() == [0, 1]
    assert out.div10.tolist() == [5.0, 10.0]
    assert out.div10_sq.tolist() == [25.0, 100.0]
    assert out.len10.tolist() == [30.0, 28.0]
    assert out.exonic.tolist() == [0, 1]


def _synthetic(n=8000, carrier_log_or=-0.3, seed=0):
    """Carriers get a known log-odds shift; divergence is a real confounder."""
    rng = np.random.default_rng(seed)
    div10 = rng.uniform(2, 20, n)
    carrier = (rng.random(n) < 1 / (1 + np.exp(-(0.15 * div10 - 2.0)))).astype(int)
    lin = 0.2 - 0.05 * div10 + carrier_log_or * carrier
    sense = (rng.random(n) < 1 / (1 + np.exp(-lin))).astype(int)
    return pd.DataFrame({"sense": sense, "carrier_any": carrier, "div10": div10,
                         "div10_sq": div10 ** 2, "len10": rng.normal(29, 1, n)})


def test_model_recovers_a_planted_carrier_effect():
    d = _synthetic(carrier_log_or=-0.3)
    r = summarise(fit(d, "carrier_any"), "carrier_any", "synthetic", len(d))
    assert r["ci_lo"] < np.exp(-0.3) < r["ci_hi"]
    assert r["p"] < 0.01


def test_model_reports_no_effect_when_none_was_planted():
    d = _synthetic(carrier_log_or=0.0, seed=3)
    r = summarise(fit(d, "carrier_any"), "carrier_any", "synthetic", len(d))
    assert r["ci_lo"] < 1.0 < r["ci_hi"]
    assert r["p"] > 0.05


def test_ignoring_the_confounder_biases_the_estimate_but_the_model_does_not():
    """Divergence drives both carrier status and sense, so the crude OR is wrong."""
    d = _synthetic(carrier_log_or=-0.3, seed=7)
    crude = sm_crude(d)
    adj = np.exp(fit(d, "carrier_any").params["carrier_any"])
    truth = np.exp(-0.3)
    assert abs(adj - truth) < abs(crude - truth)


def sm_crude(d):
    a = ((d.carrier_any == 1) & (d.sense == 0)).sum()
    b = ((d.carrier_any == 1) & (d.sense == 1)).sum()
    c = ((d.carrier_any == 0) & (d.sense == 0)).sum()
    e = ((d.carrier_any == 0) & (d.sense == 1)).sum()
    return (b * c) / (a * e)


def test_summarise_returns_an_odds_ratio_not_a_log_odds():
    d = _synthetic(carrier_log_or=-0.5)
    r = summarise(fit(d, "carrier_any"), "carrier_any", "x", len(d))
    assert 0 < r["OR"] < 1 and r["ci_lo"] < r["OR"] < r["ci_hi"]
