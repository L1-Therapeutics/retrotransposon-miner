"""Per-call transduction detection: name the donor, or refuse to.

Why this module exists
----------------------
:mod:`retro_miner.transduction_tags` answers *which source element feeds a set of
calls* once each call's transduced segment is known; it stops short of deciding
whether a given call's junction sequence is a transduced segment at all.  The
marker catalog carries the population rates (L1 3' 8.8-15%, SVA 5' 9.2%) and
records the gap explicitly: "no per-call detector -- needs a source-element
index".  This module is that detector, and this module builds that index.

The distinction it exists to make.  An insertion that is longer than its element
is only evidence of read-through if the *extra* sequence came from somewhere
else.  Three cases look identical in a length field and are not the same claim:

  * **partnered** -- the extra sequence matches a transducible flank of some
    *other* reference element of the same family.  The donor can be named, and
    the transduced length measured.
  * **orphan** -- the extra sequence is absent from the reference altogether, so
    no donor can be named.  It is still a transduction candidate (the donor may
    be a non-reference or somatically rearranged locus, or the assembly may
    contain errors), and it is never upgraded to a named source.
  * **overhang** -- the extra sequence matches the reference *at the insertion's
    own locus*.  That is a mis-anchored alignment, not a transduction.

Damert et al. 2009 localized a source locus for only 207 of 220 SVA 5'
transductions and left the rest unassigned; Pickeral 2000 and Goodier 2000
established the L1 3' equivalent (Goodier's mean transduced length ~207 nt,
observed range 22-1,388 nt; Szak 2003 tightened the TSD-verified rate to ~9%).
The refusal to name a source on a tie, on a short segment, or on low complexity
is inherited from that practice and from ``transduction_tags``, not invented
here.

Gates, each one a place this measurement could lie
--------------------------------------------------
* **Length floor.**  A segment under ``min_segment_bp`` cannot carry enough
  independent k-mer seeds to outrun chance, so it is ``unevaluable`` -- never
  ``not_transduction``, which would claim absence from silence.
* **Poly(A/T) extension.**  A pure A/T run at the junction is the signature of a
  transcript that terminated at its own signal, so it is positive evidence of a
  *non*-transduced junction and reports as ``not_transduction``.  Where the
  segment's own terminal bases are A/T the boundary between segment and tail is
  not knowable from sequence alone, so the whole terminal run is treated as tail
  and every reported transduced length is a lower bound.
* **Low complexity.**  Repeat-derived sequence matches many donors equally and
  is screened out before donor search, because a false donor is worse than none.
* **Both orientations.**  A stored junction sequence may be the reverse
  complement of the element sense direction; donor search therefore runs on the
  segment and its reverse complement and reports which one matched.
* **Tie refusal.**  Two donors with equal seed evidence produce
  ``transduction_source_ambiguous``: a candidate transduction with no named
  source.
* **Matched-span specificity.**  The seed floors are counts, and a count can be
  met inside a very short span: at ``k = 25`` a 28 bp shared piece yields exactly
  the 4 seed windows ``DEFAULT_MIN_SEED_SITES`` asks for.  Such a piece is not
  specific to one copy -- the one this gate was calibrated against is carried by
  46 genome-wide 3'-competent L1 flanks -- so naming the top-scoring carrier
  reports which copies happened to be indexed, not which element transduced.
  A match shorter than ``min_matched_bp`` therefore produces
  ``transduction_source_ambiguous`` with reason
  ``matched_span_below_specificity_floor``: still a candidate transduction, with
  no nameable donor.

Reference access is a :class:`~retro_miner.transduction_tags.SequenceFetcher`
(1-based inclusive), so the detector is pure logic over a sequence source.  Two
adapters ship here: :class:`TwoBitReference` (UCSC ``.2bit``, because that is
what this project keeps on disk and no 2bit library is installed) and
:func:`fasta_fetcher` (pyfaidx).

References
----------
Damert et al. 2009   Genome Res 19:1992. doi:10.1101/gr.093435.109
Pickeral et al. 2000 Genome Res 10:1188 (L1 3' transduction, 8.8%).
Goodier et al. 2000  Genome Res 10:1178 (mean transduced length ~207 nt).
Szak et al. 2003     J Mol Biol 327:521 (TSD-verified ~9%).
"""

from __future__ import annotations

import argparse
import gzip
import struct
from dataclasses import dataclass
from pathlib import Path
from typing import IO, Iterable, Iterator, Mapping, Sequence, cast

from retro_miner.transduction_tags import (
    DEFAULT_KMER,
    DEFAULT_MIN_SEED_HITS,
    DEFAULT_MIN_SEED_SITES,
    DEFAULT_MIN_UNIQUE_FRACTION,
    SequenceFetcher,
    SourceElement,
    TransductionSeedIndex,
    build_source_index,
    extract_source_flank,
    uniqueness_score,
)

__all__ = [
    "OUTCOMES",
    "DEFAULT_MIN_SEGMENT_BP",
    "DEFAULT_MIN_TAIL_RUN",
    "DEFAULT_MIN_MATCHED_BP",
    "DEFAULT_END_TOLERANCE_BP",
    "DEFAULT_MIN_MODEL_COVERAGE",
    "PARTNERED",
    "SOURCE_AMBIGUOUS",
    "ORPHAN_CANDIDATE",
    "NOT_TRANSDUCTION",
    "UNEVALUABLE",
    "TwoBitReference",
    "fasta_fetcher",
    "rmsk_rows",
    "consensus_coords",
    "candidate_source_elements",
    "build_annotation_index",
    "strip_poly_at_ends",
    "find_donors",
    "DonorCandidate",
    "TransductionCall",
    "detect_junction_transduction",
    "annotate_transduction",
]

# --- outcome vocabulary ------------------------------------------------------

PARTNERED = "partnered_transduction"
SOURCE_AMBIGUOUS = "transduction_source_ambiguous"
ORPHAN_CANDIDATE = "orphan_transduction_candidate"
NOT_TRANSDUCTION = "not_transduction"
UNEVALUABLE = "unevaluable"

#: Closed vocabulary, in the order the detector tests for them.
OUTCOMES = (PARTNERED, SOURCE_AMBIGUOUS, ORPHAN_CANDIDATE, NOT_TRANSDUCTION, UNEVALUABLE)

#: Below this, no segment can carry enough independent seeds to be tested; the
#: question is unanswered rather than answered negatively.
DEFAULT_MIN_SEGMENT_BP = 30

#: Shortest A/T run counted as a poly(A/T) tail at a junction.
DEFAULT_MIN_TAIL_RUN = 8

#: How close to a consensus model's end a RepeatMasker copy must reach to count
#: as carrying that end.  A read-through source needs its 3' end present; a
#: flank-utilizing source needs its 5' end present.
DEFAULT_END_TOLERANCE_BP = 150

#: Fraction of the consensus model a candidate source copy must cover.  A
#: fragment-derived row is not a competent source element.
DEFAULT_MIN_MODEL_COVERAGE = 0.95

#: Shortest seed-supported span that can single out one donor.  Set at two k-mer
#: lengths because a span of one k-mer's width carries one seed window whatever
#: the seed count says, and the seed floors are satisfied by counts alone: on real
#: hg38 a 28 bp span cleared ``min_seed_sites`` (4 windows at k = 25) while being
#: carried by 46 genome-wide 3'-competent L1 flanks.  A named donor is worth less
#: than none, so a span under this floor is refused rather than resolved by rank.
DEFAULT_MIN_MATCHED_BP = 2 * DEFAULT_KMER

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def _bases(seq: str) -> str:
    return "".join(ch for ch in (seq or "").upper() if ch in "ACGT")


# --- reference sequence access ----------------------------------------------


class TwoBitReference:
    """Read UCSC ``.2bit`` sequence; callable as a ``SequenceFetcher``.

    Format per Kent's specification: a 16-byte header (signature
    ``0x1A412743`` little-endian, version, sequence count, reserved), then one
    index entry per sequence (1-byte name length, name, 4-byte offset), then per
    sequence a record holding the DNA size, the N-block and mask-block tables,
    and the 2-bit packed bases.  Packed codes are T=00, C=01, A=10, G=11, most
    significant base first within each byte.  N blocks override the packed bases
    and are the only source of ``N`` in the output.

    Mask blocks (soft-masking) are parsed past but not applied: this class
    returns uppercase ACGTN, so downstream code cannot mistake case for evidence.
    Coordinates on the call protocol are 1-based inclusive.
    """

    _SIGNATURE = 0x1A412743

    def __init__(self, path: str | Path) -> None:
        self.path = Path(path)
        self._index: dict[str, tuple[int, int]] | None = None
        self._handle: IO[bytes] | None = None
        self._records: dict[str, tuple[int, int, tuple[tuple[int, int], ...]]] = {}

    def _open(self) -> IO[bytes]:
        """One handle, kept open: base-by-base fetching is the normal access pattern."""
        if self._handle is None or self._handle.closed:
            self._handle = self.path.open("rb")
        return self._handle

    def close(self) -> None:
        if self._handle is not None and not self._handle.closed:
            self._handle.close()
        self._handle = None

    def __enter__(self) -> "TwoBitReference":
        return self

    def __exit__(self, *_exc: object) -> None:
        self.close()

    def _load_index(self) -> dict[str, tuple[int, int]]:
        if self._index is not None:
            return self._index
        index: dict[str, tuple[int, int]] = {}
        with self.path.open("rb") as handle:   # index is read once; records reuse _open()
            signature, version, count, _reserved = struct.unpack("<IIII", handle.read(16))
            if signature != self._SIGNATURE:
                raise ValueError(f"{self.path}: not a 2bit file (signature {signature:#x})")
            if version != 0:
                raise ValueError(f"{self.path}: unsupported 2bit version {version}")
            # Read the whole index before seeking: each entry's record offset is
            # followed by the *next* entry's name, so seeking to a record mid-loop
            # would land the next read inside that record (it did -- the real hg38
            # file caught this, and a one-sequence test file cannot).
            entries: list[tuple[str, int]] = []
            for _ in range(count):
                name_len = handle.read(1)[0]
                name = handle.read(name_len).decode("ascii")
                offset = struct.unpack("<I", handle.read(4))[0]
                entries.append((name, offset))
            for name, offset in entries:
                handle.seek(offset)
                dna_size = struct.unpack("<I", handle.read(4))[0]
                index[name] = (offset, int(dna_size))
        self._index = index
        return index

    def chroms(self) -> tuple[str, ...]:
        return tuple(self._load_index())

    def _resolve_name(self, chrom: str) -> str | None:
        index = self._load_index()
        for candidate in (chrom, f"chr{chrom}", chrom[3:] if chrom.startswith("chr") else chrom):
            if candidate in index:
                return candidate
        return None

    def _layout(self, name: str) -> tuple[int, int, tuple[tuple[int, int], ...]]:
        """Return ``(packed_offset, dna_size, n_blocks)`` for a chromosome, cached."""
        cached = self._records.get(name)
        if cached is not None:
            return cached
        offset, dna_size = self._load_index()[name]
        handle = self._open()
        handle.seek(offset)
        handle.read(4)  # dna size, already known from the index
        n_block_count = struct.unpack("<I", handle.read(4))[0]
        raw = struct.unpack(f"<{2 * n_block_count}I", handle.read(8 * n_block_count)) \
            if n_block_count else ()
        # The table is two arrays -- all `nBlockStarts`, then all `nBlockSizes` --
        # not interleaved (start, size) pairs.  Pairing the halves against each
        # other feeds starts in as sizes, which turns each block's *end* into the
        # next block's *size*: the real hg38 file then reports megabase-scale N
        # runs that swallow real sequence (chr22:45.59 Mb went all-N, which is
        # where the fixture reads align).  A one-block test file parses either
        # way, so only a multi-block table catches it.
        n_blocks = tuple(zip(raw[:n_block_count], raw[n_block_count:]))
        mask_count = struct.unpack("<I", handle.read(4))[0]
        handle.read(8 * mask_count)  # mask starts + sizes, parsed past
        handle.read(4)  # reserved
        layout = (handle.tell(), int(dna_size), n_blocks)
        self._records[name] = layout
        return layout

    def fetch(self, chrom: str, start: int, end: int) -> str:
        """Return reference bases for 1-based inclusive ``[start, end]``."""
        name = self._resolve_name(chrom)
        if name is None or start > end:
            return ""
        packed_offset, dna_size, n_blocks = self._layout(name)
        start = max(1, int(start))
        end = min(int(end), dna_size)
        if end < start:
            return ""
        n_start = start - 1  # 0-based inclusive
        n_end = end          # 0-based exclusive
        handle = self._open()
        first_byte = n_start // 4
        last_byte = (n_end - 1) // 4
        handle.seek(packed_offset + first_byte)
        packed = handle.read(last_byte - first_byte + 1)

        codes = ("T", "C", "A", "G")
        bases: list[str] = []
        # Decode the whole byte-aligned window, then trim to the request.
        for i, byte in enumerate(packed):
            for shift in (6, 4, 2, 0):
                bases.append(codes[(byte >> shift) & 0b11])
        lead = n_start - 4 * first_byte
        seq = bases[lead : lead + (n_end - n_start)]

        for block_start, block_size in n_blocks:
            block_end = block_start + block_size
            lo = max(block_start, n_start)
            hi = min(block_end, n_end)
            for pos in range(lo, hi):
                seq[pos - n_start] = "N"
        return "".join(seq)

    def __call__(self, chrom: str, start: int, end: int) -> str:
        return self.fetch(chrom, start, end)


def fasta_fetcher(path: str | Path) -> SequenceFetcher:
    """pyfaidx-backed ``SequenceFetcher``; the adapter ``transduction_tags`` documents."""
    import pyfaidx  # lazy: the module stays importable without it

    fasta = pyfaidx.Fasta(str(path), as_raw=True, sequence_always_upper=True)

    def fetch(chrom: str, start: int, end: int) -> str:
        try:
            return str(fasta[chrom][start - 1 : end])
        except KeyError:
            if chrom.startswith("chr"):
                return str(fasta[chrom[3:]][start - 1 : end])
            return str(fasta[f"chr{chrom}"][start - 1 : end])

    return fetch


# --- source-element index from an annotation ---------------------------------

#: UCSC rmsk table columns, with the leading ``bin`` a real table carries.
RMSK_COLUMNS = (
    "bin", "swScore", "milliDiv", "milliDel", "milliIns", "genoName", "genoStart",
    "genoEnd", "genoLeft", "strand", "repName", "repClass", "repFamily", "repStart",
    "repEnd", "repLeft", "id",
)


def rmsk_rows(path: str | Path) -> Iterator[dict[str, str]]:
    """Yield RepeatMasker rows as dicts, handling the optional leading ``bin``."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(str(path), "rt", encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            # A real table has no header; skip one if it is there.  The names are
            # matched against the leading column, never a numeric bin.
            if len(fields) == len(RMSK_COLUMNS) and fields[0] == RMSK_COLUMNS[0]:
                continue
            if len(fields) == len(RMSK_COLUMNS) - 1 and fields[0] == RMSK_COLUMNS[1]:
                continue
            if len(fields) == len(RMSK_COLUMNS):
                yield dict(zip(RMSK_COLUMNS, fields))
            elif len(fields) == len(RMSK_COLUMNS) - 1:
                yield dict(zip(RMSK_COLUMNS[1:], fields))


def family_of_repclass(rep_class: str, rep_family: str) -> str:
    """Map RepeatMasker class/family onto the two families with transducible flanks."""
    rep_class, rep_family = rep_class.strip(), rep_family.strip()
    if rep_class == "LINE" and rep_family == "L1":
        return "LINE1"
    if rep_class == "Retroposon" and rep_family == "SVA":
        return "SVA"
    return ""


def consensus_coords(row: Mapping[str, str]) -> tuple[int, int, int]:
    """Return ``(copy_start, copy_end, model_len)`` in consensus coordinates.

    UCSC stores ``repStart``/``repEnd``/``repLeft`` in a strand-dependent order:
    on a '+' row ``repStart`` is the consensus start and ``repLeft`` the
    remainder; on a '-' row the first and third are swapped, so ``repStart``
    carries the remainder and ``repLeft`` the consensus start.  Reading a '-'
    row with the '+' formula adds a coordinate to a length and inflates the
    model, which is the bug ``artifacts/catalog.py`` documents for the chr22
    catalogue.  Either way ``repEnd`` is the copy's far consensus coordinate, so
    model length is ``repEnd + |remainder|``.
    """
    rep_start = int(row["repStart"])
    rep_end = int(row["repEnd"])
    rep_left = int(row["repLeft"])
    if str(row["strand"]) == "-":
        copy_start, remainder = rep_left, rep_start
    else:
        copy_start, remainder = rep_start, rep_left
    model_len = rep_end + abs(remainder)
    return copy_start, rep_end, model_len


def _carries_transducible_end(family: str, copy_start: int, copy_end: int, model_len: int,
                             end_tolerance_bp: int) -> bool:
    """Does this copy carry the end a transducing transcript must include?

    An L1 transduces its 3' flank, so the copy needs the *consensus 3' end*; an
    SVA transduces its 5' flank, so the copy needs the *consensus 5' end*.  The
    tolerance absorbs RepeatMasker's alignment ends not landing exactly on the
    model boundary.
    """
    if family == "LINE1":
        return copy_end >= model_len - end_tolerance_bp
    if family == "SVA":
        return copy_start <= end_tolerance_bp
    return False


def candidate_source_elements(
    rows: Iterable[Mapping[str, str]],
    fetch: SequenceFetcher,
    *,
    families: Sequence[str] = ("LINE1", "SVA"),
    end_tolerance_bp: int = DEFAULT_END_TOLERANCE_BP,
    min_model_coverage: float = DEFAULT_MIN_MODEL_COVERAGE,
    min_unique_fraction: float = DEFAULT_MIN_UNIQUE_FRACTION,
    max_elements: int | None = None,
) -> list[SourceElement]:
    """Turn annotation rows into indexable candidate source elements.

    A row becomes a source only when it is a same-family MEI, full enough to be
    a competent source (``min_model_coverage`` of its consensus model), and
    carries the transducible end (see
    :func:`_carries_transducible_end`).  Flank extraction, strand handling and
    the k-mer uniqueness gate are ``transduction_tags``' own; a flank too
    repetitive to seed is dropped there rather than guessed at here.
    """
    wanted = {fam.upper() for fam in families}
    elements: list[SourceElement] = []
    for row in rows:
        family = family_of_repclass(str(row.get("repClass", "")), str(row.get("repFamily", "")))
        if family not in wanted:
            continue
        try:
            copy_start, copy_end, model_len = consensus_coords(row)
        except (KeyError, ValueError):
            continue
        if model_len <= 0:
            continue
        coverage = (copy_end - copy_start + 1) / float(model_len)
        if coverage < min_model_coverage:
            continue
        if not _carries_transducible_end(family, copy_start, copy_end, model_len,
                                         end_tolerance_bp):
            continue
        chrom = str(row["genoName"])
        start0, end0 = int(row["genoStart"]), int(row["genoEnd"])
        if end0 <= start0:
            continue
        element = extract_source_flank(
            source_id=f"{chrom}:{start0}-{end0}:{row.get('repName', '')}",
            family=family,
            chrom=chrom,
            strand=str(row["strand"]),
            element_start0=start0,
            element_end0=end0,
            fetch=fetch,
            min_unique_fraction=min_unique_fraction,
        )
        if element is not None:
            elements.append(element)
            if max_elements is not None and len(elements) >= max_elements:
                break
    return elements


def build_annotation_index(
    annotation_path: str | Path,
    fetch: SequenceFetcher,
    *,
    families: Sequence[str] = ("LINE1", "SVA"),
    k: int = DEFAULT_KMER,
    **kwargs: object,
) -> tuple[TransductionSeedIndex, list[SourceElement]]:
    """Build the source-element index the catalog asked for.

    Returns the index and the elements that went into it, so a caller can report
    how many candidate sources a given annotation yielded rather than only
    whether a call matched one.
    """
    elements = candidate_source_elements(rmsk_rows(annotation_path), fetch, families=families,
                                         **kwargs)  # type: ignore[arg-type]
    return build_source_index(elements, k=k), elements


# --- per-call detection ------------------------------------------------------


def strip_poly_at_ends(seq: str, min_run: int = DEFAULT_MIN_TAIL_RUN) -> tuple[str, int, str]:
    """Strip leading/trailing A/T homopolymers; return ``(core, run_bp, sides)``.

    ``sides`` is ``"5p"``, ``"3p"``, ``"5p,3p"`` or ``""``.  A junction read
    presents the tail at the 3' end of the copy; a read stored in the reverse
    orientation presents it at the 5' end of the stored string, so both ends are
    checked.  ``run_bp`` is the longest run removed, which is the tail length the
    evidence supports.  A run shorter than ``min_run`` is left in place: trimming
    an incidental A or two would move the segment boundary on no evidence.
    """
    s = _bases(seq)
    if not s:
        return "", 0, ""
    sides: list[str] = []
    longest = 0
    head = 0
    while head < len(s) and s[head] in "AT":
        head += 1
    if head >= min_run:
        sides.append("5p")
        longest = max(longest, head)
        s = s[head:]
    tail = len(s)
    while tail > 0 and s[tail - 1] in "AT":
        tail -= 1
    run = len(s) - tail
    if run >= min_run:
        sides.append("3p")
        longest = max(longest, run)
        s = s[:tail]
    return s, longest, ",".join(sides)


def _dominant_fraction(seq: str) -> float:
    if not seq:
        return 1.0
    counts = {base: seq.count(base) for base in "ACGT"}
    return max(counts.values()) / float(len(seq))


@dataclass(frozen=True)
class DonorCandidate:
    """One candidate source element a junction segment matches.

    ``matched_bp`` is the seed-supported span inside the oriented segment, which
    is a *lower* bound on the transduced length: seeds only exist where windows
    matched, so genuine transduced bases beyond the outermost seed are not
    counted.  ``seg_start``/``seg_end`` locate that span in the segment as it was
    searched (i.e. in ``orientation``), and ``donor_start0``/``donor_end0``
    locate the matched region in the donor's transducible flank.
    """

    source_id: str
    family: str
    chrom: str
    strand: str
    flank_dir: str
    orientation: str  # "sense" (segment as given) or "antisense" (its revcomp)
    n_seed_sites: int
    n_hits: int
    seg_start: int
    seg_end: int
    donor_start0: int
    donor_end0: int

    @property
    def matched_bp(self) -> int:
        return self.seg_end - self.seg_start


def find_donors(
    segment: str,
    family: str,
    index: TransductionSeedIndex,
    *,
    min_seed_sites: int = DEFAULT_MIN_SEED_SITES,
    min_hits: int = DEFAULT_MIN_SEED_HITS,
) -> list[DonorCandidate]:
    """Rank donor candidates for a segment, searching both orientations.

    Returns candidates clearing both seed floors, best first.  Ties are *not*
    resolved here: the caller must refuse to name a source when two candidates
    carry equal evidence, which is Damert's rule.
    """
    fam = family.upper()
    best_per_source: dict[str, DonorCandidate] = {}
    for orientation, search_seq in (("sense", segment), ("antisense", _revcomp(segment))):
        if not search_seq:
            continue
        for source_id, kmer_hits in index.lookup(search_seq).items():
            element = index.sources()[source_id]
            if element.family != fam:
                continue
            sites = {kmer for kmer, _c, _f in kmer_hits}
            if len(sites) < min_seed_sites or len(kmer_hits) < min_hits:
                continue
            call_offsets = [c for _k, c, _f in kmer_hits]
            flank_offsets = [f for _k, _c, f in kmer_hits]
            k = index.k
            candidate = DonorCandidate(
                source_id=source_id,
                family=element.family,
                chrom=element.chrom,
                strand=element.strand,
                flank_dir=element.flank_dir,
                orientation=orientation,
                n_seed_sites=len(sites),
                n_hits=len(kmer_hits),
                seg_start=min(call_offsets),
                seg_end=max(call_offsets) + k,
                donor_start0=element.flank_start0 + min(flank_offsets),
                donor_end0=element.flank_start0 + max(flank_offsets) + k,
            )
            previous = best_per_source.get(source_id)
            if previous is None or (candidate.n_seed_sites, candidate.n_hits) > (
                    previous.n_seed_sites, previous.n_hits):
                best_per_source[source_id] = candidate
    return sorted(best_per_source.values(),
                  key=lambda c: (-c.n_seed_sites, -c.n_hits, c.source_id, c.orientation))


def _is_tie(ranked: Sequence[DonorCandidate]) -> bool:
    if len(ranked) < 2:
        return False
    return (ranked[0].n_seed_sites, ranked[0].n_hits) == (ranked[1].n_seed_sites,
                                                          ranked[1].n_hits)


@dataclass(frozen=True)
class TransductionCall:
    """One call's transduction verdict.  Frozen, and every field is evidence."""

    call_id: str
    family: str
    outcome: str
    reason: str
    junction_bp: int
    segment_bp: int
    tail_run_bp: int
    tail_side: str
    transduced_bp: int | None
    observed_excess_bp: int | None
    unexplained_excess_bp: int | None
    donor: DonorCandidate | None
    n_donor_candidates: int

    @property
    def is_transduction(self) -> bool:
        """True for partnered and ambiguous outcomes; an orphan stays a candidate."""
        return self.outcome in (PARTNERED, SOURCE_AMBIGUOUS)

    def to_dict(self) -> dict[str, object]:
        donor = self.donor
        return {
            "call_id": self.call_id,
            "family": self.family,
            "outcome": self.outcome,
            "reason": self.reason,
            "junction_bp": self.junction_bp,
            "segment_bp": self.segment_bp,
            "tail_run_bp": self.tail_run_bp,
            "tail_side": self.tail_side,
            "transduced_bp": self.transduced_bp,
            "observed_excess_bp": self.observed_excess_bp,
            "unexplained_excess_bp": self.unexplained_excess_bp,
            "n_donor_candidates": self.n_donor_candidates,
            "donor_source_id": donor.source_id if donor else "",
            "donor_chrom": donor.chrom if donor else "",
            "donor_strand": donor.strand if donor else "",
            "donor_flank_dir": donor.flank_dir if donor else "",
            "donor_matched_start0": donor.donor_start0 if donor else "",
            "donor_matched_end0": donor.donor_end0 if donor else "",
            "donor_orientation": donor.orientation if donor else "",
            "donor_seed_sites": donor.n_seed_sites if donor else "",
            "donor_seed_hits": donor.n_hits if donor else "",
        }


def detect_junction_transduction(
    call_id: str,
    family: str,
    junction_seq: str,
    index: TransductionSeedIndex | None,
    *,
    observed_excess_bp: int | None = None,
    local_reference_seq: str = "",
    min_segment_bp: int = DEFAULT_MIN_SEGMENT_BP,
    min_tail_run: int = DEFAULT_MIN_TAIL_RUN,
    min_seed_sites: int = DEFAULT_MIN_SEED_SITES,
    min_hits: int = DEFAULT_MIN_SEED_HITS,
    min_matched_bp: int = DEFAULT_MIN_MATCHED_BP,
    max_dominant_fraction: float = 0.9,
    min_unique_fraction: float = DEFAULT_MIN_UNIQUE_FRACTION,
) -> TransductionCall:
    """Decide whether a junction sequence is a transduced segment.

    ``junction_seq`` is the call's non-element sequence at the junction, 5'->3'
    in the orientation it was recovered (either orientation is accepted).
    ``observed_excess_bp`` is the length by which the call exceeds its element,
    when a caller can measure one; it is recorded and, for partnered calls,
    compared against the seed-supported transduced length so an excess larger
    than its explanation stays visible.  ``local_reference_seq`` is the reference
    sequence immediately beyond the insertion site: a segment matching *that* is
    an overhang, not a transduction.  ``min_matched_bp`` is the shortest span a
    donor may be named from; see :data:`DEFAULT_MIN_MATCHED_BP`.
    """
    junction = _bases(junction_seq)
    core, tail_run, tail_side = strip_poly_at_ends(junction, min_run=min_tail_run)

    def call(outcome: str, reason: str, *, donor: DonorCandidate | None = None,
             n_candidates: int = 0, transduced_bp: int | None = None,
             segment_bp: int | None = None) -> TransductionCall:
        seg_bp = len(core) if segment_bp is None else segment_bp
        unexplained: int | None = None
        if observed_excess_bp is not None and transduced_bp is not None:
            unexplained = max(int(observed_excess_bp) - int(transduced_bp), 0)
        return TransductionCall(
            call_id=call_id, family=family.upper(), outcome=outcome, reason=reason,
            junction_bp=len(junction), segment_bp=seg_bp, tail_run_bp=tail_run,
            tail_side=tail_side, transduced_bp=transduced_bp,
            observed_excess_bp=observed_excess_bp, unexplained_excess_bp=unexplained,
            donor=donor, n_donor_candidates=n_candidates,
        )

    if not junction:
        return call(UNEVALUABLE, "no_junction_sequence")
    if index is None or len(index) == 0:
        return call(UNEVALUABLE, "no_source_index")
    if tail_run > 0 and len(core) == 0:
        # The whole junction is an A/T run: the transcript was polyadenylated
        # where its own signal told it to stop.  Positive evidence, not silence.
        return call(NOT_TRANSDUCTION, "polyA_extension_only", transduced_bp=0)
    if len(core) < min_segment_bp:
        return call(UNEVALUABLE, "segment_below_floor"
                    if tail_run == 0 else "segment_below_floor_after_tail_strip")
    if _dominant_fraction(core) > max_dominant_fraction:
        return call(NOT_TRANSDUCTION, "low_complexity_segment", transduced_bp=0)
    if uniqueness_score(core) < min_unique_fraction:
        # A repeat-derived segment matches many donors equally; it is screened out
        # rather than allowed to name one of them.
        return call(NOT_TRANSDUCTION, "repeat_like_segment", transduced_bp=0)

    ranked = find_donors(core, family, index, min_seed_sites=min_seed_sites,
                         min_hits=min_hits)
    local = _bases(local_reference_seq)
    if not ranked:
        # A segment drawn from the reference next to the insertion is a
        # mis-anchored overhang, not a transduction.  Checked only when no donor
        # matched, so a genuine donor is never overruled by local similarity.
        if local and _local_containment(core, local, index.k) >= 0.5:
            return call(NOT_TRANSDUCTION, "overhang_matches_local_reference", transduced_bp=0)
        return call(ORPHAN_CANDIDATE, "no_reference_match", n_candidates=0,
                    transduced_bp=len(core))
    specific = [candidate for candidate in ranked if candidate.matched_bp >= min_matched_bp]
    if not specific:
        # Matched, but only over a span several unrelated copies also carry: the
        # candidates are reported as a count, and the span as the lower bound it
        # is, without turning whichever one ranked first into a named source.
        return call(SOURCE_AMBIGUOUS, "matched_span_below_specificity_floor",
                    n_candidates=len(ranked),
                    transduced_bp=max(candidate.matched_bp for candidate in ranked))
    ranked = specific
    if _is_tie(ranked):
        return call(SOURCE_AMBIGUOUS, "donor_tie", n_candidates=len(ranked),
                    transduced_bp=max(ranked[0].matched_bp, ranked[1].matched_bp))
    best = ranked[0]
    return call(PARTNERED, "donor_match", donor=best, n_candidates=len(ranked),
                transduced_bp=best.matched_bp)


def _local_containment(segment: str, local: str, k: int) -> float:
    """Fraction of the segment's k-mers that occur in the local reference.

    The overhang question is containment, not similarity: a mis-anchored segment
    *is* the reference continuation, so nearly all of its k-mers appear in the
    local sequence, while an unrelated 25-mer appears there essentially never.
    """
    if len(segment) < k or len(local) < k:
        return 0.0
    seg_kmers = {segment[i : i + k] for i in range(len(segment) - k + 1)}
    local_kmers = {local[i : i + k] for i in range(len(local) - k + 1)}
    if not seg_kmers:
        return 0.0
    return len(seg_kmers & local_kmers) / float(len(seg_kmers))


def annotate_transduction(
    calls: Iterable[Mapping[str, object]],
    index: TransductionSeedIndex | None,
    *,
    call_id_key: str = "call_id",
    family_key: str = "family",
    junction_key: str = "junction_seq",
    excess_key: str = "observed_excess_bp",
    local_key: str = "local_reference_seq",
    **kwargs: object,
) -> list[TransductionCall]:
    """Run :func:`detect_junction_transduction` over a table of calls, in order.

    Missing optional fields are passed through as absent (``None`` / ``""``), so
    a call with no measurable excess or no local reference is still evaluated
    rather than skipped.
    """
    out: list[TransductionCall] = []
    for row in calls:
        excess = row.get(excess_key)
        try:
            excess_bp: int | None = None if excess in (None, "") else int(float(str(excess)))
        except (TypeError, ValueError):
            excess_bp = None
        out.append(
            detect_junction_transduction(
                str(row.get(call_id_key, "")),
                str(row.get(family_key, "")),
                str(row.get(junction_key, "")),
                index,
                observed_excess_bp=excess_bp,
                local_reference_seq=str(row.get(local_key, "") or ""),
                **kwargs,  # type: ignore[arg-type]
            )
        )
    return out


# --- CLI ---------------------------------------------------------------------


def _open_reference(path: str) -> SequenceFetcher:
    if str(path).endswith(".2bit"):
        return TwoBitReference(path)
    return fasta_fetcher(path)


def main(argv: Sequence[str] | None = None) -> int:
    """Build a source index from an annotation and annotate a calls table."""
    import pandas as pd

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--rmsk", required=True, help="RepeatMasker table (plain or .gz)")
    parser.add_argument("--reference", required=True, help="reference FASTA or .2bit")
    parser.add_argument("--calls", required=True,
                        help="TSV with call_id, family, junction_seq and optional "
                             "observed_excess_bp, local_reference_seq")
    parser.add_argument("--out", required=True, help="output TSV")
    parser.add_argument("--families", default="LINE1,SVA")
    parser.add_argument("--min-model-coverage", type=float, default=DEFAULT_MIN_MODEL_COVERAGE)
    parser.add_argument("--min-segment-bp", type=int, default=DEFAULT_MIN_SEGMENT_BP)
    parser.add_argument("--min-matched-bp", type=int, default=DEFAULT_MIN_MATCHED_BP)
    parser.add_argument("--end-tolerance-bp", type=int, default=DEFAULT_END_TOLERANCE_BP)
    parser.add_argument("--max-elements", type=int, default=None)
    args = parser.parse_args(argv)

    fetch = _open_reference(args.reference)
    index, elements = build_annotation_index(
        args.rmsk, fetch,
        families=tuple(part for part in args.families.split(",") if part),
        min_model_coverage=args.min_model_coverage,
        end_tolerance_bp=args.end_tolerance_bp,
        max_elements=args.max_elements,
    )
    calls = pd.read_csv(args.calls, sep="\t", low_memory=False)
    results = annotate_transduction(
        cast("Iterable[Mapping[str, object]]", calls.to_dict("records")),
        index,
        min_segment_bp=args.min_segment_bp,
        min_matched_bp=args.min_matched_bp,
    )
    frame = pd.DataFrame([r.to_dict() for r in results])
    frame.to_csv(args.out, sep="\t", index=False)

    print(f"source elements indexed: {len(elements)} (families {args.families})")
    print(f"calls annotated: {len(frame)}")
    if not frame.empty:
        for outcome, count in frame["outcome"].value_counts().items():
            print(f"  {outcome}: {count}")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
