from __future__ import annotations

import pickle

import numpy as np
import pandas as pd

import retro_miner.mei_support as mei_support
from retro_miner.local_assembly import (
    _build_locus_row_payloads,
    _interval_from_row,
    _window_locus_id_from_row,
)
from retro_miner.mei_support import (
    ClipAlignmentSummary,
    _candidate_subset_for_sample,
    _empty_mei_remap_result,
    _run_named_process_jobs,
    _sample_has_assigned_rows,
)



def _dummy_named_worker(*, sample: str, payload: list[int]) -> dict[str, object]:
    return {"sample": sample, "total": sum(payload)}



def test_build_locus_row_payloads_preserves_window_helpers_and_pickles():
    frame = pd.DataFrame(
        [
            {
                "chrom": "chr1",
                "window_start": np.int64(100),
                "window_end": np.int64(160),
                "insertion_breakpoint_pos": np.int64(125),
                "support_fraction": np.float64(0.75),
                "is_pass": np.bool_(True),
                "maybe_missing": pd.NA,
                "tags": ["split", np.int64(2)],
                "meta": {"depth": np.int64(7)},
            }
        ]
    )

    payload = _build_locus_row_payloads(frame)[0]
    series = frame.iloc[0]

    assert type(payload["window_start"]) is int
    assert type(payload["support_fraction"]) is float
    assert type(payload["is_pass"]) is bool
    assert payload["maybe_missing"] is None
    assert payload["tags"] == ["split", 2]
    assert payload["meta"] == {"depth": 7}
    pickle.dumps(payload)
    assert _window_locus_id_from_row(payload) == _window_locus_id_from_row(series)
    assert _interval_from_row(payload, 25) == _interval_from_row(series, 25)



def test_sample_has_assigned_rows_requires_nonempty_frame():
    empty = pd.DataFrame(columns=["chrom", "window_start", "window_end"])
    nonempty = pd.DataFrame([{"chrom": "chr1", "window_start": 100, "window_end": 200}])

    assert not _sample_has_assigned_rows(None, empty)
    assert _sample_has_assigned_rows(empty, nonempty)



def test_candidate_subset_for_sample_filters_loci_in_candidate_order():
    candidates = pd.DataFrame(
        [
            {"chrom": "chr1", "window_start": 10, "window_end": 20, "score": 1},
            {"chrom": "chr2", "window_start": 30, "window_end": 40, "score": 2},
            {"chrom": "chr3", "window_start": 50, "window_end": 60, "score": 3},
        ]
    )
    split_df = pd.DataFrame(
        [
            {"chrom": "chr3", "window_start": 50, "window_end": 60},
            {"chrom": "chr1", "window_start": 10, "window_end": 20},
        ]
    )
    discordant_df = pd.DataFrame(columns=["chrom", "window_start", "window_end"])

    subset = _candidate_subset_for_sample(candidates, split_df, discordant_df)

    assert subset[["chrom", "window_start", "window_end", "score"]].to_dict(orient="records") == [
        {"chrom": "chr1", "window_start": 10, "window_end": 20, "score": 1},
        {"chrom": "chr3", "window_start": 50, "window_end": 60, "score": 3},
    ]



def test_empty_mei_remap_result_returns_zero_summaries():
    result = _empty_mei_remap_result("disease")

    assert result["split_hits"].empty
    assert result["disc_hits"].empty
    assert result["split_summary"] == ClipAlignmentSummary(sample="disease", clip_count=0, paf_hits=0)
    assert result["disc_summary"] == ClipAlignmentSummary(sample="disease", clip_count=0, paf_hits=0)
    assert result["disc_mate_summary"] == ClipAlignmentSummary(sample="disease_mate", clip_count=0, paf_hits=0)



def test_filtered_named_jobs_skip_process_pool_when_only_one_sample_has_rows(monkeypatch):
    disease_split = pd.DataFrame(columns=["chrom", "window_start", "window_end"])
    disease_discordant = pd.DataFrame(columns=["chrom", "window_start", "window_end"])
    control_split = pd.DataFrame([{"chrom": "chr2", "window_start": 30, "window_end": 40}])
    control_discordant = pd.DataFrame(columns=["chrom", "window_start", "window_end"])
    jobs: dict[str, dict[str, object]] = {}
    if _sample_has_assigned_rows(disease_split, disease_discordant):
        jobs["disease"] = {"sample": "disease", "payload": [1, 2, 3]}
    if _sample_has_assigned_rows(control_split, control_discordant):
        jobs["control"] = {"sample": "control", "payload": [4, 5]}

    class FailPool:
        def __init__(self, *args, **kwargs):
            raise AssertionError("ProcessPoolExecutor should not be constructed")

    monkeypatch.setattr(mei_support, "ProcessPoolExecutor", FailPool)

    result = _run_named_process_jobs(_dummy_named_worker, jobs, max_workers=2)

    assert result == {"control": {"sample": "control", "total": 9}}
