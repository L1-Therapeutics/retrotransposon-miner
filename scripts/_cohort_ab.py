"""A/B the cohort rule on identical inputs. Scratch driver, not an analysis.

The regenerated Phase 4 artifacts differ from the previously published ones in
two ways at once: the cohort rule changed, *and* the underlying table grew from
5 genomes to 10. Comparing those two artifact sets measures both at once and
attributes neither. This driver holds the input fixed and varies only the cohort
predicate, so the difference is attributable to the rule alone.

    python _cohort_ab.py /tmp/ab_nested              # the settled rule
    python _cohort_ab.py /tmp/ab_sense_only sense_only  # the withdrawn rule

The sense-only arm works by patching `common.is_nested_site` in a generated
wrapper before the analysis module is run. Both scripts do
`import nested_multi_sample_common as common`, which resolves to the already
patched module object, so the override reaches them without either script
knowing about it -- the same thing that let the two definitions drift apart.
"""

from __future__ import annotations

import json
import subprocess
import sys
import textwrap
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))

import nested_multi_sample_common as common  # noqa: E402

US = common.DEFAULT_UNIQUE_SITES
CALLSETS = common.DEFAULT_CALLSET_DIR

_WRAPPER = textwrap.dedent(
    f"""
    import runpy, sys
    sys.path.insert(0, {str(SCRIPTS)!r})
    import nested_multi_sample_common as common
    def is_nested_site(state):
        return (state or "").strip() == "nested_sense"
    common.is_nested_site = is_nested_site
    sys.argv = sys.argv[1:]
    runpy.run_path(sys.argv[0], run_name="__main__")
    """
)


def run(outdir: Path, rule: str, replicates: str = "10000", draws: str = "20000") -> dict:
    outdir.mkdir(parents=True, exist_ok=True)
    if rule == "nested":
        entry = sys.executable
        prefix: list[str] = []
    elif rule == "sense_only":
        wrapper = outdir / "_sense_only_wrapper.py"
        wrapper.write_text(_WRAPPER)
        entry = sys.executable
        prefix = [str(wrapper)]
    else:
        raise SystemExit(f"unknown rule: {rule}")

    subprocess.run(
        [entry]
        + prefix
        + [
            str(SCRIPTS / "joint_enrichment.py"),
            "--unique-sites", str(US),
            "--callset-dir", str(CALLSETS),
            "--outdir", str(outdir),
            "--replicates", replicates,
            "--bootstrap", "2000",
            "--seed", "20261001",
        ],
        check=True,
        capture_output=True,
    )
    subprocess.run(
        [entry]
        + prefix
        + [
            str(SCRIPTS / "recurrence_test.py"),
            "--unique-sites", str(US),
            "--callset-dir", str(CALLSETS),
            "--outdir", str(outdir),
            "--chance-draws", draws,
            "--seed", "20261001",
        ],
        check=True,
        capture_output=True,
    )
    joint = json.loads((outdir / "joint_enrichment.json").read_text())
    rec = json.loads((outdir / "recurrence.json").read_text())
    return {
        "rule": rule,
        "sites_in_cohort": joint["cohort_definition"]["sites_in_cohort"],
        "joint_cells": {
            c["cell_id"]: [
                c["observed"],
                round(c["expected_null_mean"], 2),
                c["monte_carlo_p"],
            ]
            for c in joint["cells"]
        },
        "joint_position_factor": round(
            joint["multiplication_test"].get("position_factor") or 0, 3
        ),
        "joint_orientation_factor": round(
            joint["multiplication_test"].get("orientation_factor") or 0, 3
        ),
        "recurrence": {
            "sites": rec["nested_unique_sites"],
            "candidate": rec["any_candidate_recurrent"],
            "ibd": rec["verdict_summary"]["ibd_same_event"]["n_pairs"],
            "unevaluable": rec["verdict_summary"]["tsd_unevaluable"]["n_pairs"],
            "chance_expected": round(
                rec["positional_chance"]["expected_pairs_under_chance"], 2
            ),
            "hosts_considered": rec["positional_chance"]["hosts_considered"],
            "per_family_total_multi_site": sum(
                b["hosts_with_more_than_one_nested_event"]
                for b in rec["per_family"].values()
            ),
        },
    }


if __name__ == "__main__":
    out = Path(sys.argv[1])
    rule = sys.argv[2] if len(sys.argv) > 2 else "nested"
    print(json.dumps(run(out, rule), indent=2, sort_keys=True))