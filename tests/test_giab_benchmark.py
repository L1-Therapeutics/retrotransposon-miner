from benchmarks.benchmark_giab_truthset import evaluate_benchmark_parity


def test_perfect_parity_benchmark():
    truth = [("chr1", 1000, "L1HS"), ("chr1", 5000, "ALU")]
    cand = [("chr1", 1005, "L1HS"), ("chr1", 4995, "ALU")] # Within 50 bp window

    res = evaluate_benchmark_parity(truth, cand, window_bp=50)
    assert res["precision"] == 1.0
    assert res["recall"] == 1.0
    assert res["f1_score"] == 1.0

def test_imperfect_calls_benchmark():
    truth = [("chr1", 1000, "L1HS"), ("chr1", 5000, "ALU")]
    cand = [("chr1", 1000, "L1HS"), ("chr2", 9000, "SPURIOUS")] # 1 TP, 1 FP, 1 FN

    res = evaluate_benchmark_parity(truth, cand, window_bp=50)
    assert res["true_positives"] == 1
    assert res["false_positives"] == 1
    assert res["false_negatives"] == 1
    assert res["precision"] == 0.5
    assert res["recall"] == 0.5
