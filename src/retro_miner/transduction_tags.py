"""L1/SVA transduction source tagging: link MEI calls to their source element.

Retrotransposition read-through transduces flank DNA along with the element
body.  A *shared unique transduced segment* links multiple insertions to one
active source element: Damert et al. 2009 built the SVA 5' transduction
groups (220/2,398 TSD'd SVAs, 93 source loci, MAST2 group = 78 members) with
exactly this logic — "junctions between the 5'-ends of the truncated Alu-like
region and the 5' transductions are identical, and 5'-transduced sequences
are overlapping and can be traced back to the same source locus".  Pickeral
et al. 2000 and Goodier et al. 2000 established the L1 3' equivalent.

Architecture (mirrors :mod:`retro_miner.mei_markers` provenance discipline):

* **Source-flank extraction** (:class:`SourceElement`, :func:`extract_source_flank`)
  — strand-correct: an L1 on the '+' strand reads through its 3' end into
  *downstream* reference; an SVA on '+' reads through its 5' end from
  *upstream* reference; '-' strands mirror both.
* **Seed index** (:class:`TransductionSeedIndex`) — k-mer seeds from unique
  source-flank sequence, k chosen so a shared segment has many independent
  seed hits (Damert criterion 1: shared junction; criterion 2: overlapping
  transduced sequence).
* **Matching** (:class:`TransducedSegment`, :func:`match_call_to_source`,
  :func:`annotate_callset`) — a call is tagged only when the evidence clears
  the seed floor AND the best source is unambiguous; ties are "ambiguous"
  because guessing a source would be exactly the unvalidated claim the rest
  of the pipeline refuses to make.
* **Group linking** (:func:`build_transduction_groups`) — union-find over
  shared-source links so multi-member clusters emerge, MAST2-style.

Repetitive-sequence discipline: flanks are screened before seeding
(:func:`uniqueness_score`); low-complexity/short flanks never enter the
index, because a repeat-derived seed would link unrelated calls — a
false "source" assignment is worse than no assignment.

References
----------
Damert et al. 2009   Genome Res 19:1992. doi:10.1101/gr.093435.109
Pickeral et al. 2000 Genome Res 10:1188 (L1 3' transduction, 8.8%).
Goodier et al. 2000  Genome Res 10:1178 (L1 3' transduction, mean ~207 nt).
Szak et al. 2003     J Mol Biol 327:521 (stricter TSD-verified ~9%).
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from typing import Callable, Iterable, Mapping, Sequence

__all__ = [
    "DEFAULT_KMER",
    "DEFAULT_MIN_SEED_HITS",
    "DEFAULT_MIN_SEED_SITES",
    "DEFAULT_MIN_UNIQUE_FRACTION",
    "DEFAULT_L1_FLANK_MAX_BP",
    "DEFAULT_SVA_FLANK_MAX_BP",
    "SequenceFetcher",
    "SourceElement",
    "TransducedSegment",
    "TransductionGroup",
    "uniqueness_score",
    "extract_source_flank",
    "build_source_index",
    "match_call_to_source",
    "build_transduction_groups",
    "annotate_callset",
]

# --- defaults (each with its justification) --------------------------------

#: Seed k.  A 25-mer is long enough that random 25-mer matches between
#: unrelated genomic loci are rare (expected ≪1 across a callset), yet short
#: enough that a ≥100 bp shared transduced segment yields ~75 independent
#: seed windows — so a genuine link produces dozens of hits, not one.
DEFAULT_KMER = 25

#: Minimum independent seed-site hits (distinct k-mer positions, not raw
#: occurrences) before a call/source match is considered a candidate link.
DEFAULT_MIN_SEED_HITS = 4

#: Minimum distinct seed k-mers from one source hit by a call.  Requiring
#: *distinct* seeds defeats a single low-complexity stretch faking depth.
DEFAULT_MIN_SEED_SITES = 4

#: Flank bases screened before entering the seed index: fraction of
#: maximally-unique seeds a flank must support.  Damert screened transduced
#: sequence against repeats before assigning sources; this is the same gate.
DEFAULT_MIN_UNIQUE_FRACTION = 0.4

#: L1 3' transductions observed 22-1,388 bp (mean ~207 nt; Goodier 2000 /
#: Szak 2003); cap the extraction window well above the mean.
DEFAULT_L1_FLANK_MAX_BP = 1_500

#: SVA 5' transductions observed 14-2,161 bp (mean 223; Damert 2009 Fig. 4).
DEFAULT_SVA_FLANK_MAX_BP = 2_300


#: Protocol for reference-sequence access.  Coordinates are 1-based
#: *inclusive* on the forward strand — the convention pyfaidx users get from
#: ``genome[chrom][a:b].seq`` after converting, and unambiguous to test.
SequenceFetcher = Callable[[str, int, int], str]

# 0-based half-open helper used throughout: [start0, end0).
_Complement = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _revcomp(seq: str) -> str:
    return seq.translate(_Complement)[::-1]


def _bases(seq: str) -> str:
    return "".join(ch for ch in seq.upper() if ch in "ACGT")


def _window_kmers(seq: str, k: int) -> list[str]:
    """k-mers of pure ACGT windows (N runs break the window, not the frame)."""
    kmers: list[str] = []
    run: list[str] = []
    for ch in seq.upper():
        if ch in "ACGT":
            run.append(ch)
        else:
            if len(run) >= k:
                kmers.extend("".join(run[i : i + k]) for i in range(len(run) - k + 1))
            run = []
    if len(run) >= k:
        kmers.extend("".join(run[i : i + k]) for i in range(len(run) - k + 1))
    return kmers


def uniqueness_score(seq: str, k: int = DEFAULT_KMER) -> float:
    """Fraction of k-mer windows in *seq* that are unique 25-mers (no repeated
    25-mer within the sequence itself).  Cheap intra-sequence screen; real
    uniqueness is enforced downstream by the seed-floor mechanics (a call
    matches a repetitive source only if the *call's own* unique segment
    overlaps the flank).
    """
    kmers = _window_kmers(seq, k)
    if not kmers:
        return 0.0
    counts = Counter(kmers)
    return sum(1 for n in counts.values() if n == 1) / len(kmers)


# --- source element model ----------------------------------------------------

#: Direction the transduced flank extends from the element body.
FLANK_3P = "3p"  # L1: read-through past the poly(A) site
FLANK_5P = "5p"  # SVA: upstream sequence under an external promoter


@dataclass(frozen=True)
class SourceElement:
    """A candidate active source element with its transducible flank.

    ``element_start0``/``element_end0`` are 0-based half-open reference
    coordinates of the element body on the *forward* strand; ``strand`` is
    the element's own orientation.  ``flank_start0``/``flank_end0`` bound the
    extracted transducible flank in reference coordinates (already strand-
    oriented by :func:`extract_source_flank`).
    """

    source_id: str
    family: str  # "LINE1" or "SVA"
    chrom: str
    strand: str  # "+" or "-"
    element_start0: int
    element_end0: int
    flank_start0: int
    flank_end0: int
    flank_seq: str  # reference sequence, 5'->3' along the element's sense
    flank_dir: str  # FLANK_3P or FLANK_5P
    uniqueness: float = 0.0

    @property
    def flank_len(self) -> int:
        return len(self.flank_seq)


def extract_source_flank(
    source_id: str,
    family: str,
    chrom: str,
    strand: str,
    element_start0: int,
    element_end0: int,
    fetch: SequenceFetcher,
    *,
    flank_max_bp: int | None = None,
    min_unique_fraction: float = DEFAULT_MIN_UNIQUE_FRACTION,
) -> SourceElement | None:
    """Extract the transducible flank for one candidate source element.

    L1 reads through its 3' end (downstream flank on the element's sense
    strand); SVA reads through its 5' end (upstream flank).  Strand is the
    *element's* orientation: for '-' elements the transduced flank is the
    reverse complement of the forward-strand reference beyond the element
    boundary.  Returns None when the flank is unusable (too short, too
    repetitive, or boundary coordinates out of range) — an unusable source is
    silently absent from the index rather than guessed at.
    """
    fam = family.upper()
    if fam == "LINE1":
        dirn = FLANK_3P
        default_max = DEFAULT_L1_FLANK_MAX_BP
    elif fam == "SVA":
        dirn = FLANK_5P
        default_max = DEFAULT_SVA_FLANK_MAX_BP
    else:
        return None
    max_bp = flank_max_bp if flank_max_bp is not None else default_max

    if strand == "+":
        if dirn == FLANK_3P:
            f_start, f_end = element_end0, element_end0 + max_bp
        else:
            f_start, f_end = max(0, element_start0 - max_bp), element_start0
    elif strand == "-":
        if dirn == FLANK_3P:
            f_start, f_end = max(0, element_start0 - max_bp), element_start0
        else:
            f_start, f_end = element_end0, element_end0 + max_bp
    else:
        return None
    if f_end <= f_start:
        return None

    raw = fetch(chrom, f_start + 1, f_end)  # protocol is 1-based inclusive
    if not raw:
        return None
    raw = _bases(raw)
    flank = raw if strand == "+" else _revcomp(raw)

    # Uniqueness gate: screen the *sense-oriented* flank before indexing.
    frac = uniqueness_score(flank)
    if frac < min_unique_fraction:
        return None
    return SourceElement(
        source_id=source_id,
        family=fam,
        chrom=chrom,
        strand=strand,
        element_start0=element_start0,
        element_end0=element_end0,
        flank_start0=f_start,
        flank_end0=f_end,
        flank_seq=flank,
        flank_dir=dirn,
        uniqueness=frac,
    )


# --- seed index ----------------------------------------------------------------


@dataclass
class TransductionSeedIndex:
    """k-mer -> source-element membership index over unique source flanks.

    Each k-mer maps to a set of ``(source_id, offset_in_flank)`` pairs, so a
    call's matching segment can be located within the source flank, not just
    named.
    """

    k: int
    _seeds: dict[str, set[tuple[str, int]]] = field(default_factory=dict)
    _sources: dict[str, SourceElement] = field(default_factory=dict)

    def add(self, element: SourceElement) -> None:
        self._sources[element.source_id] = element
        for offset, kmer in enumerate(_window_kmers(element.flank_seq, self.k)):
            self._seeds.setdefault(kmer, set()).add((element.source_id, offset))

    def __contains__(self, source_id: str) -> bool:
        return source_id in self._sources

    def sources(self) -> Mapping[str, SourceElement]:
        return self._sources

    def __len__(self) -> int:
        return len(self._sources)

    def lookup(self, seq: str) -> dict[str, list[tuple[str, int, int]]]:
        """Return {source_id: [(kmer, call_offset, flank_offset), ...]}.

        ``call_offset`` is where the seed hit inside *seq*; ``flank_offset``
        is where that k-mer sits in the source flank — together they locate
        the shared segment on both sides of the junction.
        """
        hits: dict[str, list[tuple[str, int, int]]] = {}
        for call_offset, kmer in enumerate(_window_kmers(seq, self.k)):
            for source_id, flank_offset in self._seeds.get(kmer, ()):
                hits.setdefault(source_id, []).append((kmer, call_offset, flank_offset))
        return hits


def build_source_index(
    elements: Iterable[SourceElement],
    *,
    k: int = DEFAULT_KMER,
) -> TransductionSeedIndex:
    """Build a seed index from extracted source flanks."""
    index = TransductionSeedIndex(k=k)
    for element in elements:
        index.add(element)
    return index


# --- call matching ---------------------------------------------------------------


@dataclass(frozen=True)
class TransducedSegment:
    """One call's transduced-segment match to a source element.

    ``hit_offsets`` are positions in the *call's* 5'-to-3' sequence where
    seed k-mers from the source flank were found; their span approximates
    the transduced segment's extent inside the call.
    """

    source_id: str
    family: str
    n_seed_sites: int  # distinct k-mers shared with the source flank
    n_hits: int  # total hits (a repeated k-mer counts per occurrence)
    call_span0: int  # first hit offset in the call sequence
    call_span1: int  # one past the last hit offset
    flank_span0: int  # matching region inside the source flank
    flank_span1: int
    unique: bool  # single unambiguous best source


def match_call_to_source(
    call_id: str,
    family: str,
    call_seq: str,
    index: TransductionSeedIndex,
    *,
    min_seed_sites: int = DEFAULT_MIN_SEED_SITES,
    min_hits: int = DEFAULT_MIN_SEED_HITS,
) -> TransducedSegment | None:
    """Tag one call with its feeding source element, or None.

    ``call_seq`` is the call's non-element sequence 5'->3' in the call's own
    orientation (e.g. a locally assembled insertion allele with the element
    body removed, or the extra sequence the caller recovered beyond the
    consensus alignment).  A source must be the same family and clear both
    seed floors; if two sources tie at the maximum, the result is None — a
    guessed source is worse than no source (Damert only assigned the 207 of
    220 transductions whose sequences mapped back to a localized source
    locus, and left the rest unassigned).
    """
    if not call_seq:
        return None
    fam = family.upper()
    hits = index.lookup(call_seq)
    scored: list[tuple[str, int, int, list[int], list[int]]] = []
    for source_id, kmer_hits in hits.items():
        element = index.sources()[source_id]
        if element.family != fam:
            continue
        sites = {kmer for kmer, _c, _f in kmer_hits}
        call_offsets = [c for _k, c, _f in kmer_hits]
        flank_offsets = [f for _k, _c, f in kmer_hits]
        scored.append((source_id, len(sites), len(kmer_hits), call_offsets, flank_offsets))
    if not scored:
        return None
    scored.sort(key=lambda item: (-item[1], -item[2], item[0]))
    best = scored[0]
    if best[1] < min_seed_sites or best[2] < min_hits:
        return None
    if len(scored) > 1 and scored[1][1] == best[1] and scored[1][2] == best[2]:
        # Tie at the top: refuse to name a source.
        return None
    source_id, n_sites, n_hits, call_offsets, flank_offsets = best
    call_offsets_sorted = sorted(call_offsets)
    flank_offsets_sorted = sorted(flank_offsets)
    return TransducedSegment(
        source_id=source_id,
        family=fam,
        n_seed_sites=n_sites,
        n_hits=n_hits,
        call_span0=call_offsets_sorted[0],
        call_span1=call_offsets_sorted[-1] + index.k,
        flank_span0=flank_offsets_sorted[0],
        flank_span1=flank_offsets_sorted[-1] + index.k,
        unique=True,
    )


# --- group linking ---------------------------------------------------------------


@dataclass
class TransductionGroup:
    """A union of calls linked to one source element (or transitively)."""

    group_id: str
    source_ids: frozenset[str]
    call_ids: frozenset[str]

    @property
    def size(self) -> int:
        return len(self.call_ids)


def build_transduction_groups(
    assignments: Sequence[tuple[str, str]],
) -> list[TransductionGroup]:
    """Union-find over (call_id, source_id) assignments.

    Two calls land in one group when they share a source element, which is
    Damert's group definition; transitive chains (A->S1, B->S2 with S1,S2
    linked elsewhere) are out of scope here — a secondary-transduction chain
    is a sequence-level claim this module does not make from seeds alone.
    """
    parent: dict[str, str] = {}

    def find(x: str) -> str:
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: str, b: str) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for call_id, source_id in assignments:
        union(f"call:{call_id}", f"src:{source_id}")

    groups: dict[str, tuple[set[str], set[str]]] = {}
    for node in list(parent):
        root = find(node)
        calls, srcs = groups.setdefault(root, (set(), set()))
        if node.startswith("call:"):
            calls.add(node[5:])
        else:
            srcs.add(node[4:])
    out: list[TransductionGroup] = []
    for i, (root, (calls, srcs)) in enumerate(sorted(groups.items()), start=1):
        out.append(
            TransductionGroup(
                group_id=f"TG{i:04d}",
                source_ids=frozenset(srcs),
                call_ids=frozenset(calls),
            )
        )
    return out


# --- callset entry point ---------------------------------------------------------


def annotate_callset(
    calls: Sequence[Mapping[str, object]],
    index: TransductionSeedIndex,
    *,
    call_id_key: str = "call_id",
    family_key: str = "family",
    seq_key: str = "extra_seq",
) -> list[TransducedSegment | None]:
    """Tag every call in *calls* (order preserved) with its transduction match.

    Each call mapping needs ``call_id``, ``family``, and ``extra_seq`` (the
    call's non-element 5'->3' sequence).  Unmatched calls get None; the
    caller decides how to serialize that (absent column vs "unassigned").
    """
    out: list[TransducedSegment | None] = []
    for call in calls:
        out.append(
            match_call_to_source(
                str(call.get(call_id_key, "")),
                str(call.get(family_key, "")),
                str(call.get(seq_key, "")),
                index,
            )
        )
    return out
