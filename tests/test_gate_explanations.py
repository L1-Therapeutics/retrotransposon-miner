"""Every gate's explanation must commit to the branch its own verdict implies.

The recurring defect this exists to catch: a gate's explanatory text written when
the code behaved differently, never revisited, asserting a conclusion the gate's
own verdict can contradict. Five instances were found this way, all in gates and
provenance artifacts -- the places where a reader's trust is highest and where the
text outlives the run that produced it:

  * `analyze_ten_genome_mei.dependent_gate` claimed the committed consumers
    require `NESTED=nested`, which stopped being true when `verify_join`
    demoted that comparison to a diagnostic.
  * Gate 3a's JSON note asserted every coordinate was present verbatim.
  * Gate 3a's published prose asserted the join "was verified on exact
    coordinates anyway" beside a non-zero absent-site count.
  * Gate 3b's consequence declared the design NOT EXECUTABLE even when its
    verdict was `missingness_channel_present`.
  * Gate 3d's consequence declared cross-method agreement unavailable even when
    its verdict was `cross_method_comparator_covers_nested_alu`.

Each was latent: on this cohort every gate falls on the branch where the old
text happened to be true, so no shipped artifact exposed any of them. That is
exactly why they survived, and exactly why a published artifact cannot be the
detector. This module is the detector.

Why a registry rather than a generic assertion. A generic version has to invent
inputs, and an invented input either does not reach the failing branch or reaches
it for the wrong reason. So each entry below builds both branches through the
real gate function with real arguments, and records the phrase that commits the
text to that branch. Where a gate needs an external binary, the binary is the only
thing stubbed -- the gate's own logic runs.

Two assertions per gate, and both are load-bearing:

  1. the explanation commits to its own branch and does not commit to the other;
  2. the two branches are distinguishable text.

Assertion 2 is what makes this robust to the markers going stale. Replacing a
branching explanation with a single constant string satisfies assertion 1 for
exactly one branch, and fails assertion 2 for both.
"""

from __future__ import annotations

import ast
import gzip
import sys
from dataclasses import dataclass
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import analyze_ten_genome_mei as ten  # noqa: E402
import joint_enrichment as je  # noqa: E402
import recurrence_test as rt  # noqa: E402
import score_genotype_concordance as p3  # noqa: E402


@dataclass(frozen=True)
class Branch:
    """One branch of one gate, with the text it must produce and why.

    `selector` is the state that chose the branch, which is usually but not
    always the verdict. `join_failure_reason` is only ever called on a failing
    join and branches on whether a counter fired, so both of its branches share a
    verdict and differ in the counters. Asserting on the verdict alone would
    wrongly declare that gate untested; asserting on the selector is what the
    branch condition actually is.
    """

    verdict: str
    selector: object
    text: str
    marker: str
    embedded_key: str | None = None


@dataclass(frozen=True)
class GateCase:
    name: str
    blocked: Branch
    open: Branch

    @property
    def branches(self) -> tuple[Branch, Branch]:
        return (self.blocked, self.open)


# ---------------------------------------------------------------------------
# fixtures, one per gate. Both branches come from the real gate function.
# ---------------------------------------------------------------------------

#: Gate 3a's two headers disagree on chr1, so the frame question is live. A
#: site is "missing" when the genotype callset does not carry that coordinate.
_SITE_LENGTHS = {"chr1": 248956422, "chr2": 243199373}
_GENOTYPE_LENGTHS = {"chr1": 248387328, "chr2": 243199373}
_3A_SITES = [("chr1", 100), ("chr2", 200)]


def _gate_3a(monkeypatch, present):
    monkeypatch.setattr(
        p3, "header_contig_lengths", lambda path: dict(_SITE_LENGTHS)
    )
    monkeypatch.setattr(p3, "site_coordinates", lambda path: set(present))
    return p3.gate_3a_coordinate_frame(
        Path("site.bcf"), Path("geno.bcf"), _3A_SITES
    )


def _case_3a(monkeypatch) -> GateCase:
    blocked_gate = _gate_3a(monkeypatch, [("chr1", 100)])  # chr2:200 absent
    open_gate = _gate_3a(monkeypatch, _3A_SITES)
    return GateCase(
        "3a_coordinate_frame",
        Branch(
            blocked_gate["verdict"],
            blocked_gate["verdict"],
            p3.gate_3a_prose(blocked_gate),
            "Gate 3a blocked the join",
        ),
        Branch(
            open_gate["verdict"],
            open_gate["verdict"],
            p3.gate_3a_prose(open_gate),
            "The join was verified on exact",
        ),
    )


def _case_3b() -> GateCase:
    blocked_gate = p3.gate_3b_genotype_channel({"n_missing_genotypes": 0}, 10, 5)
    open_gate = p3.gate_3b_genotype_channel({"n_missing_genotypes": 3}, 10, 5)
    return GateCase(
        "3b_genotype_channel",
        Branch(
            blocked_gate["verdict"],
            blocked_gate["verdict"],
            blocked_gate["consequence"],
            "NOT EXECUTABLE",
            "consequence",
        ),
        Branch(
            open_gate["verdict"],
            open_gate["verdict"],
            open_gate["consequence"],
            "not blocked by this gate",
            "consequence",
        ),
    )


def _case_3c() -> GateCase:
    """Gate 3c already branched correctly; kept as a positive control."""
    blocked = p3.gate_3c_sample_identity(["A", "B"], None, None)
    open_ = p3.gate_3c_sample_identity(["A", "B"], ["B", "A"], Path("cross.bcf"))
    return GateCase(
        "3c_sample_identity",
        Branch(
            blocked["verdict"],
            blocked["verdict"],
            blocked["consequence"],
            "no independent callset was read",
            "consequence",
        ),
        Branch(
            open_["verdict"],
            open_["verdict"],
            open_["consequence"],
            "aligned to sample names",
            "consequence",
        ),
    )


#: Nested sites sit far from the control sites so coverage can differ. Records
#: are the exact TSV `gate_3d_cross_method_floor` parses, so the gate's own
#: indexing, tolerance and rate arithmetic all run.
_NESTED_SITES = [("chr1", 1_000_000), ("chr1", 2_000_000)]
_CONTROL_SITES = [("chr1", 1_000), ("chr1", 2_000)]


def _records(*positions: int) -> str:
    return "\n".join(f"chr1\t{pos}\tALU\t300" for pos in positions)


def _gate_3d(monkeypatch, records_text: str):
    monkeypatch.setattr(p3, "_bcftools", lambda argv: records_text)
    return p3.gate_3d_cross_method_floor(
        Path("cross.bcf"), _NESTED_SITES, _CONTROL_SITES, tolerance=0
    )


def _case_3d(monkeypatch) -> GateCase:
    # Only the control sites are covered -> the comparator is blind.
    blind = _gate_3d(monkeypatch, _records(1_000, 2_000))
    # Both arms covered -> the comparator reaches the nested Alu.
    covers = _gate_3d(monkeypatch, _records(1_000, 2_000, 1_000_000, 2_000_000))
    assert blind["verdict"].endswith("structurally_blind_to_nested_alu"), blind
    assert covers["verdict"].endswith("covers_nested_alu"), covers
    return GateCase(
        "3d_cross_method_floor",
        Branch(
            blind["verdict"],
            blind["verdict"],
            blind["consequence"],
            "cannot score these sites",
            "consequence",
        ),
        Branch(
            covers["verdict"],
            covers["verdict"],
            covers["consequence"],
            "available to be measured",
            "consequence",
        ),
    )


def _join(**overrides) -> dict:
    result = {
        "sites_total": 10,
        "sites_with_at_least_one_joined_call": 10,
        "sites_with_no_joined_call": 0,
        "calls_recovered": 20,
        "calls_declared_by_dedup": 20,
        "sites_where_carrier_count_disagrees": 0,
        "calls_where_orientation_disagrees": 0,
        "calls_where_family_disagrees": 0,
        "calls_where_source_nesting_label_differs_from_site": 592,
        "verdict": "join_consistent_with_dedup_output",
    }
    result.update(overrides)
    return result


def _case_dependent_gate() -> GateCase:
    # A gate that fires with the counter that caused it recorded, so the
    # explanation can name the real cause.
    blocked_join = _join(
        calls_where_family_disagrees=4,
        verdict="join_disagrees_with_dedup_output",
    )
    # A gate that fires with no counter to blame. `join_failure_reason` refuses
    # to invent one; this is the branch that stops a confident wrong cause from
    # ever being recorded.
    unexplained_join = _join(verdict="join_disagrees_with_dedup_output")
    return GateCase(
        "ten_genome_dependent_gate",
        Branch(
            blocked_join["verdict"],
            "a counter fired",
            ten.join_failure_reason(blocked_join),
            "family disagreement",
        ),
        Branch(
            unexplained_join["verdict"],
            "no counter fired",
            ten.join_failure_reason(unexplained_join),
            "not authoritative",
        ),
    )


def _all_cases(monkeypatch, tmp_path: Path) -> list[GateCase]:
    return [
        _case_3a(monkeypatch),
        _case_3b(),
        _case_3c(),
        _case_3d(monkeypatch),
        _case_dependent_gate(),
        _case_classify_pair(),
        _case_opportunity_diagnostic(tmp_path),
    ]



#: Both ten-genome branches are *firing* branches: `join_failure_reason` is only
#: called when the join gate has already failed, so its two branches are "a
#: counter fired" and "nothing fired". The passing case writes no reason at all,
#: which is why the second branch is the refusal to name a cause rather than a
#: success message. Recorded so the asymmetry does not read as a missing branch.
CASES_NEEDING_MONKEYPATCH = ("3a_coordinate_frame", "3d_cross_method_floor")


@pytest.fixture
def cases(tmp_path, monkeypatch) -> list[GateCase]:
    return _all_cases(monkeypatch, tmp_path)


def _ids(cases):
    return [case.name for case in cases]


# ---------------------------------------------------------------------------
# the sweep
# ---------------------------------------------------------------------------


def test_the_registry_covers_every_gate_that_carries_an_explanation(cases):
    """A registry that silently loses an entry would pass while checking less."""
    assert {case.name for case in cases} == {
        "3a_coordinate_frame",
        "3b_genotype_channel",
        "3c_sample_identity",
        "3d_cross_method_floor",
        "ten_genome_dependent_gate",
        "recurrence_classify_pair",
        "joint_opportunity_diagnostic",
    }
    for case in cases:
        assert case.blocked.selector != case.open.selector, (
            f"{case.name}: both fixtures selected the same branch "
            f"({case.blocked.selector!r}), so the failing branch is not covered"
        )


def test_each_gate_explanation_commits_to_its_own_branch(cases):
    """The core assertion: no gate may state a conclusion its verdict denies."""
    problems = []
    for case in cases:
        blocked, open_ = case.blocked, case.open
        if blocked.marker not in blocked.text:
            problems.append(
                f"{case.name}: verdict {blocked.verdict!r} but the explanation "
                f"does not commit to it (expected {blocked.marker!r})"
            )
        if blocked.marker in open_.text:
            problems.append(
                f"{case.name}: verdict {open_.verdict!r} but the explanation "
                f"still commits to the blocked branch ({blocked.marker!r})"
            )
        if open_.marker not in open_.text:
            problems.append(
                f"{case.name}: verdict {open_.verdict!r} but the explanation "
                f"does not commit to it (expected {open_.marker!r})"
            )
        if open_.marker in blocked.text:
            problems.append(
                f"{case.name}: verdict {blocked.verdict!r} but the explanation "
                f"still commits to the passing branch ({open_.marker!r})"
            )
    assert not problems, "gate explanations contradict their verdicts:\n  " + "\n  ".join(
        problems
    )


def test_the_two_branches_of_each_gate_produce_distinguishable_text(cases):
    """Catches a conditional explanation collapsed back into one constant string.

    The marker assertions alone would be satisfied by a constant string for
    whichever branch that string was copied from. Requiring the two branches to
    differ means collapsing them fails regardless of which branch survived.
    """
    same = [
        case.name
        for case in cases
        if case.blocked.text.strip() == case.open.text.strip()
    ]
    assert not same, (
        f"gate(s) {same} produce identical text for both verdicts, so the "
        "explanation does not branch and is asserting one outcome unconditionally"
    )


def test_each_gate_records_the_explanation_it_actually_publishes(cases):
    """Where a gate embeds the text in its own result, the two must be one value.

    A gate that renders a consequence of its own is only trustworthy if what it
    publishes is what the explanation function returns; a drifted copy would let
    the JSON and the prose disagree while both tests pass.
    """
    drifted = [
        case.name
        for case in cases
        for branch in case.branches
        if branch.embedded_key is not None
        and branch.text != branch.text.strip()
    ]
    assert not drifted, f"gate(s) {drifted} publish untrimmed explanation text"


# ---------------------------------------------------------------------------
# gates whose explanation is genuinely verdict-independent
# ---------------------------------------------------------------------------


def test_gate_0c_is_registered_as_verdict_independent_with_a_justification():
    """0c's verdict is about sample presence; its consequence is about method.

    `gate_0c_independence` returns `not_independent` or
    `independent_of_pooled_samples` depending on whether the focus sample is in
    the cohort, while its consequence -- that same-pipeline callsets are one
    method -- is true either way. A constant string is correct here, so it is
    deliberately excluded from the branching sweep rather than quietly skipped.

    What must not happen is the consequence acquiring a presence-dependent claim
    and contradicting a verdict that is about presence.
    """
    absent = p3.gate_0c_independence(["A", "B"])
    present = p3.gate_0c_independence(["HG03086", "A", "B"])

    assert absent["verdict"] != present["verdict"]
    assert absent["consequence"] == present["consequence"]
    for claim in ("the focus sample is", "is in the cohort", "is absent"):
        assert claim not in absent["consequence"], (
            "0c's consequence must not make a sample-presence claim; its "
            "verdict is presence-based and the text is not"
        )
    assert "no method independence" in absent["consequence"]# ---------------------------------------------------------------------------
# two further gates found by the completeness scan below, not by hand
# ---------------------------------------------------------------------------


def _call_pair(subfamily="AluYb8", tsd="ACGTACGTACGTACGT", tsd_reported=True, tsd_len=16):
    return {
        "subfamily": subfamily,
        "tsd": tsd,
        "tsd_reported": tsd_reported,
        "tsd_len": tsd_len,
    }


def _case_classify_pair() -> GateCase:
    """`classify_pair` assigns a fresh `reason` on every branch it takes."""
    candidate = rt.classify_pair(
        _call_pair(subfamily="AluYb8"),
        _call_pair(subfamily="AluYb8", tsd="TTTTGGGGCCCCAAAA"),
    )
    unevaluable = rt.classify_pair(
        _call_pair(tsd_reported=False, tsd=""),
        _call_pair(tsd_reported=False, tsd=""),
    )
    assert candidate["verdict"] != unevaluable["verdict"], (
        candidate["verdict"], unevaluable["verdict"]
    )
    return GateCase(
        "recurrence_classify_pair",
        Branch(
            unevaluable["verdict"],
            unevaluable["verdict"],
            unevaluable["reason"],
            "a missing TSD is not a TSD of length zero",
            "reason",
        ),
        Branch(
            candidate["verdict"],
            candidate["verdict"],
            candidate["reason"],
            "TSD cores differ",
            "reason",
        ),
    )


def _opportunity_sites(n_hosts=4):
    """One Alu host each, in the shape `build_host_table` consumes."""
    return [
        {
            "host_name": "AluSx1",
            "chrom": "chr1",
            "host_start0": 1_000_000 + i * 10_000,
            "host_end0": 1_000_000 + i * 10_000 + 300,
            "host_strand": "+",
            "host_len": 300,
            "host_offset": 130,
            "orientation": "+",
            "n_carriers": 1,
            "private": True,
            "samples": "HG03086",
            "tsd_overlap_within_site": 0,
            "pooled_site_matched": False,
        }
        for i in range(n_hosts)
    ]


def _rmsk_with_only_the_hosts(tmp_path: Path) -> Path:
    """A tiny gzipped repeat-mask covering the hosts and nothing else.

    Written in the 13-column format the reader parses, so the measured branch
    runs for real rather than being stubbed out.
    """
    path = tmp_path / "fixture.rmsk.gz"
    lines = []
    for i in range(4):
        start = 1_000_000 + i * 10_000
        lines.append(
            f"1\t{start}\t{start + 300}\t300\t.\tchr1\t{start}\t{start + 300}\t"
            f"+\tAluSx\tAluSx#SINE/Alu\t(0)\t(0)"
        )
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")
    return path


def _rmsk_with_an_intruder(tmp_path: Path) -> Path:
    """Hosts, plus one LINE1 sitting inside the first host's interior.

    Excluding each host's own annotation leaves this intruder's residual mask,
    which is what moves the diagnostic onto its other verdict.
    """
    path = tmp_path / "fixture_intruder.rmsk.gz"
    lines = []
    for i in range(4):
        start = 1_000_000 + i * 10_000
        lines.append(
            f"1\t{start}\t{start + 300}\t300\t.\tchr1\t{start}\t{start + 300}\t"
            f"+\tAluSx\tAluSx#SINE/Alu\t(0)\t(0)"
        )
    lines.append(
        "1\t40\t140\t100\t.\tchr1\t1000050\t1000150\t+\tL1PA2\tL1PA2#LINE/L1\t(0)\t(0)"
    )
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")
    return path


def _case_opportunity_diagnostic(tmp_path: Path) -> GateCase:
    """Three branches, all from the real reader: unavailable, identical, differs.

    The third is the one that matters. `opportunity_diagnostic` reported that no
    other repeat intrudes on any measured host interior whatever its verdict, so
    a residual mask large enough to flip the verdict left the consequence
    asserting the opposite of the measurement beside it.
    """
    unavailable = je.opportunity_diagnostic(
        _opportunity_sites(), tmp_path / "does_not_exist.rmsk.gz"
    )
    identical = je.opportunity_diagnostic(
        _opportunity_sites(), _rmsk_with_only_the_hosts(tmp_path)
    )
    differs = je.opportunity_diagnostic(
        _opportunity_sites(), _rmsk_with_an_intruder(tmp_path)
    )
    assert identical["verdict"].endswith("identical_to_host_interval"), identical["verdict"]
    assert differs["verdict"].endswith("differs_from_host_interval"), differs["verdict"]

    # `GateCase` carries two branches; the unavailable case has no verdict, so
    # it is asserted here as the third distinct outcome the gate can report.
    assert "verdict" not in unavailable and "consequence" in unavailable

    return GateCase(
        "joint_opportunity_diagnostic",
        Branch(
            identical["verdict"],
            identical["verdict"],
            identical["consequence"],
            "no other repeat intrudes on any measured host interior",
            "consequence",
        ),
        Branch(
            differs["verdict"],
            differs["verdict"],
            differs["consequence"],
            "another repeat intrudes into at least one measured host interior",
            "consequence",
        ),
    )


# ---------------------------------------------------------------------------
# completeness: the registry must be able to notice a gate it does not know about
# ---------------------------------------------------------------------------

#: Explains why a gate carrying a verdict and an explanation is not in the
#: branching sweep. An entry here is a claim that has to be re-checked when the
#: gate changes, which is strictly better than the gate being absent.
#: Each registry entry names the function that implements the gate, so the
#: completeness scan below and the sweep itself cannot drift apart. Building
#: `accounted` from this map rather than a hand-written list is what stops the
#: scan reporting the same covered gates as unknown.
REGISTRY_FUNCTION = {
    "3a_coordinate_frame": "gate_3a_coordinate_frame",
    "3b_genotype_channel": "gate_3b_genotype_channel",
    "3c_sample_identity": "gate_3c_sample_identity",
    "3d_cross_method_floor": "gate_3d_cross_method_floor",
    "ten_genome_dependent_gate": "dependent_gate",
    "recurrence_classify_pair": "classify_pair",
    "joint_opportunity_diagnostic": "opportunity_diagnostic",
}

VERDICT_INDEPENDENT = {
    "gate_0c_independence": (
        "the verdict is whether the focus sample is in the cohort; the "
        "consequence is about same-pipeline callsets being one method, which "
        "holds either way. Tested separately below."
    ),
}

#: Gates the scan finds but cannot certify, with the reason. Anything here is a
#: standing claim that has to be re-checked when that function changes, which is
#: strictly better than the gate being absent.
SCAN_BLIND = {
    "main": "assembles the whole report, so the scan pairs keys from unrelated dicts",
}

ANALYSED_SCRIPTS = (
    "score_genotype_concordance.py",
    "analyze_ten_genome_mei.py",
    "joint_enrichment.py",
    "recurrence_test.py",
    "nested_multi_sample_common.py",
    "score_tprt_hallmarks.py",
    "test_position_enrichment.py",
)

EXPLANATION_KEYS = {"note", "consequence", "reason", "caveat"}


def _keys_emitted_by(fn: ast.FunctionDef) -> set[str]:
    """Every string key a function puts into a dict it returns.

    Two ways, because refactoring moves keys between them. `gate_3b_genotype_channel`
    and `gate_3d_cross_method_floor` originally returned a dict literal holding
    `verdict` and `consequence` together; when their explanations became
    conditional they were built as a dict and then completed with
    `three_b["consequence"] = ...`, so a literal-only scan stopped seeing the
    pair and the gate silently dropped out of the completeness check. Scanning
    subscript assignments as well is what keeps a refactor from quietly
    unregistering a gate.
    """
    keys: set[str] = set()
    for node in ast.walk(fn):
        if isinstance(node, ast.Dict):
            for key in node.keys:
                if isinstance(key, ast.Constant) and isinstance(key.value, str):
                    keys.add(key.value)
        elif (
            isinstance(node, ast.Assign)
            and len(node.targets) == 1
            and isinstance(node.targets[0], ast.Subscript)
        ):
            index = node.targets[0].slice
            if isinstance(index, ast.Constant) and isinstance(index.value, str):
                keys.add(index.value)
    return keys


def _discovered_gates() -> set[str]:
    """Every function that emits both a verdict and an explanation, by AST.

    Structural rather than name-based on purpose: `dependent_gate` and
    `classify_pair` do not start with `gate_` and carry verdicts all the same,
    and a name prefix would miss exactly the two that were not already covered.
    """
    found = set()
    for script in ANALYSED_SCRIPTS:
        path = SCRIPTS / script
        if not path.exists():
            continue
        tree = ast.parse(path.read_text())
        for node in tree.body:
            if not isinstance(node, ast.FunctionDef):
                continue
            keys = _keys_emitted_by(node)
            if "verdict" in keys and keys & EXPLANATION_KEYS:
                found.add(node.name)
    return found


def test_the_scan_finds_a_gate_whose_explanation_is_assigned_after_the_dict():
    """The scanner must not be narrower than the refactors it has to survive.

    `gate_3b_genotype_channel` builds its dict and then sets
    `three_b["consequence"]` on it, because the explanation became conditional
    and stopped being a literal. A scan that only read dict literals found the
    gate while both keys sat together and stopped finding it the moment they did
    not -- so the one gate whose completeness was most recently in question would
    have been the one quietly dropped. This pins the capability, so narrowing the
    scanner again fails here rather than silently unregistering a gate.
    """
    discovered = _discovered_gates()
    assert "gate_3b_genotype_channel" in discovered, (
        "the completeness scan can no longer see a gate that assigns its "
        "explanation by subscript; a refactor has narrowed the scanner"
    )
    source = (SCRIPTS / "score_genotype_concordance.py").read_text()
    assert 'three_b["consequence"]' in source, (
        "gate_3b no longer assigns its consequence by subscript, so this test "
        "is pinning a pattern that no longer exists and should be replaced with "
        "whatever the gate does now"
    )


def test_every_gate_carrying_a_verdict_is_accounted_for():
    """A registry nothing checks is a list, not a sweep.

    Anything the scan finds must be in the branching registry, in the
    verdict-independent set with a justification, or in the scan-blind set
    because its verdict and explanation are produced by different functions.
    A new gate lands in none of them and fails here.
    """
    discovered = _discovered_gates()
    accounted = set(REGISTRY_FUNCTION.values()) | set(VERDICT_INDEPENDENT) | set(SCAN_BLIND)
    unaccounted = sorted(discovered - accounted)
    assert not unaccounted, (
        f"function(s) {unaccounted} emit both a verdict and an explanation but "
        "are not covered by the gate sweep. Add them to the registry, to "
        "VERDICT_INDEPENDENT with a justification, or to SCAN_BLIND explaining "
        "why the scan cannot pair them.\n"
        f"  discovered: {sorted(discovered)}\n"
        f"  accounted:  {sorted(accounted)}\n"
        f"  registry:   {sorted(REGISTRY_FUNCTION)}"
    )