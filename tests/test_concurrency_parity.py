from __future__ import annotations

from retro_miner.local_assembly import _run_ordered_process_jobs
from retro_miner.mei_support import _run_named_process_jobs


def _dummy_locus_worker(*, idx: int, total_loci: int, row_data: dict[str, object], marker: str):
    return (
        {
            "idx": idx,
            "total_loci": total_loci,
            "chrom": row_data["chrom"],
            "window_start": row_data["window_start"],
            "window_end": row_data["window_end"],
            "marker": marker,
        },
        f"log-{idx}-{row_data['chrom']}",
    )



def _dummy_named_worker(*, sample: str, payload: list[int], multiplier: int = 1):
    return {
        "sample": sample,
        "payload": [value * multiplier for value in payload],
        "total": sum(payload) * multiplier,
    }



def test_run_ordered_process_jobs_matches_serial_output():
    jobs = [
        {
            "idx": 1,
            "total_loci": 3,
            "row_data": {"chrom": "chr1", "window_start": 100, "window_end": 200},
            "marker": "alpha",
        },
        {
            "idx": 2,
            "total_loci": 3,
            "row_data": {"chrom": "chr2", "window_start": 300, "window_end": 450},
            "marker": "beta",
        },
        {
            "idx": 3,
            "total_loci": 3,
            "row_data": {"chrom": "chr3", "window_start": 500, "window_end": 650},
            "marker": "gamma",
        },
    ]

    serial = _run_ordered_process_jobs(_dummy_locus_worker, jobs, max_workers=1)
    parallel = _run_ordered_process_jobs(_dummy_locus_worker, jobs, max_workers=2)

    assert parallel == serial
    assert [result[0]["idx"] for result in parallel] == [1, 2, 3]



def test_run_ordered_process_jobs_fast_path_single_worker():
    jobs = [
        {
            "idx": 1,
            "total_loci": 1,
            "row_data": {"chrom": "chr7", "window_start": 10, "window_end": 20},
            "marker": "solo",
        }
    ]

    result = _run_ordered_process_jobs(_dummy_locus_worker, jobs, max_workers=1)

    assert result == [
        (
            {
                "idx": 1,
                "total_loci": 1,
                "chrom": "chr7",
                "window_start": 10,
                "window_end": 20,
                "marker": "solo",
            },
            "log-1-chr7",
        )
    ]



def test_run_named_process_jobs_matches_serial_output():
    jobs = {
        "disease": {"sample": "disease", "payload": [1, 2, 3], "multiplier": 2},
        "control": {"sample": "control", "payload": [4, 5], "multiplier": 3},
    }

    serial = _run_named_process_jobs(_dummy_named_worker, jobs, max_workers=1)
    parallel = _run_named_process_jobs(_dummy_named_worker, jobs, max_workers=2)

    assert parallel == serial
    assert set(parallel) == {"disease", "control"}
