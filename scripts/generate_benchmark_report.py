from pathlib import Path
import pandas as pd
import sys

def build_report(baseline_path, perf_path, out_md):
    b_df = pd.read_csv(baseline_path, sep='\t')
    p_df = pd.read_csv(perf_path, sep='\t')
    
    report = f"""# Pipeline Benchmark Performance Summary

| Metric | Baseline (`main`) | Optimized (`perf`) | Parity |
| :--- | :--- | :--- | :--- |
| **Total Candidates** | {len(b_df):,} | {len(p_df):,} | {'✅ Match' if len(b_df) == len(p_df) else '❌ Mismatch'} |
| **Mean MEI Score** | {b_df['mei_score'].mean():.3f} | {p_df['mei_score'].mean():.3f} | {'✅ Match' if abs(b_df['mei_score'].mean() - p_df['mei_score'].mean()) < 1e-4 else '❌ Mismatch'} |
"""
    Path(out_md).write_text(report)

if __name__ == "__main__":
    if len(sys.argv) > 3:
        build_report(sys.argv[1], sys.argv[2], sys.argv[3])