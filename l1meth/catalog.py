"""Build a catalogue of L1 loci addressable for locus-specific methylation.

A locus is *addressable* when (a) it retains 5'UTR sequence -- the promoter CpG
island that actually controls silencing -- and (b) it carries enough adjacent
non-repeat sequence that a short read crossing the element/genome junction can
be assigned to this copy rather than to one of the ~500,000 others.

Coordinates are real: RepeatMasker annotations and reference sequence, not
consensus models.
"""

from __future__ import annotations

from dataclasses import dataclass, field

# The L1 internal promoter occupies roughly the first 900 bp of the element.
UTR_LEN = 900
# bs-ATLAS-seq (Lanciano et al. 2024, Cell Genomics) interrogates the first 15
# promoter CpGs; we mirror that window so results are comparable.
N_PROMOTER_CPG = 15

# A "first 900 bp = promoter" rule is only meaningful when the subfamily's
# consensus model spans a full-length element (~6 kb). Some ancient subfamilies
# are fragment-derived: on chr22, 16 of 125 L1 subfamilies have consensus
# models under 5 kb, and two -- L1ME4b (892 bp) and L1ME4c (830 bp) -- are
# shorter than the promoter window itself, with L1ME3G (924 bp) barely longer.
# For those, "the first 900 bp" is effectively the whole model, so the 5'UTR
# test passes trivially and floods the catalogue with degraded elements that
# have no promoter at all.
MIN_CONSENSUS_LEN = 5000

# The modern human-type 5'UTR promoter belongs to the primate (L1P) lineage.
# Older L1M subfamilies predate the relevant 5'UTR replacements (>=8 episodes
# over ~70 My of primate evolution), so their 5' ends are a different promoter
# architecture and are not comparable CpG-for-CpG.
PRIMATE_PREFIXES = ("L1HS", "L1PA", "L1PB", "L1P")

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def consensus_length(rep_start: int, rep_end: int, rep_left: int, strand: str) -> int:
    """Length of the subfamily consensus model implied by one RepeatMasker row.

    UCSC stores the three consensus-coordinate fields in a *strand-dependent*
    order. For '+' rows they are (start, end, -remaining); for '-' rows the
    first and third are swapped, so ``rep_start`` carries the negative
    remainder and ``rep_left`` the consensus start. Applying the '+' formula
    to a '-' row therefore adds a coordinate to a length and inflates the
    result -- on chr22 that turned a 892 bp model into values up to 11,799 bp
    and produced 179 distinct "consensus lengths" for a single subfamily.

    Use :func:`subfamily_consensus_length` to reduce many rows to one value;
    consensus length is a property of the subfamily, not of an individual copy.
    """
    return rep_end + abs(rep_left if strand == "+" else rep_start)


def subfamily_consensus_length(rows) -> dict[str, int]:
    """Map subfamily -> consensus length, taking the modal value per subfamily.

    Even strand-corrected, individual rows disagree (RepeatMasker splits some
    elements across rows and consensus models are versioned), so the mode is
    the robust estimator; the median is not -- it lands on a real value only
    by luck.
    """
    seen: dict[str, dict[int, int]] = {}
    for r in rows:
        n = consensus_length(
            int(r["repStart"]), int(r["repEnd"]), int(r["repLeft"]), str(r["strand"])
        )
        seen.setdefault(str(r["repName"]), {})
        seen[str(r["repName"])][n] = seen[str(r["repName"])].get(n, 0) + 1
    return {fam: max(counts.items(), key=lambda kv: (kv[1], kv[0]))[0]
            for fam, counts in seen.items()}


@dataclass
class L1Locus:
    """One L1 copy and the promoter CpGs that can be measured on it."""

    chrom: str
    start: int          # 0-based genomic start of the element
    end: int            # genomic end (exclusive)
    strand: str         # '+' or '-'
    subfamily: str
    divergence: float   # RepeatMasker milliDiv / 1000
    cons_start: int     # consensus coordinate where this copy begins
    flank_bp: int       # largest adjacent non-repeat run (anchorability proxy)
    cons_len: int = 0   # length of this subfamily's consensus model
    cpg_genomic: list[int] = field(default_factory=list)  # CpG genomic starts,
    #                                                       ordered 5'->3' along
    #                                                       the ELEMENT

    @property
    def n_cpg(self) -> int:
        return len(self.cpg_genomic)

    @property
    def full_length_consensus(self) -> bool:
        """Is the subfamily's consensus long enough for a 5'UTR to be defined?"""
        return self.cons_len >= MIN_CONSENSUS_LEN

    @property
    def primate_lineage(self) -> bool:
        return self.subfamily.startswith(PRIMATE_PREFIXES)

    @property
    def has_utr(self) -> bool:
        return self.full_length_consensus and 0 <= self.cons_start < UTR_LEN

    def addressable(
        self, min_flank: int = 150, min_cpg: int = 3, primate_only: bool = True
    ) -> bool:
        if primate_only and not self.primate_lineage:
            return False
        return self.has_utr and self.flank_bp >= min_flank and self.n_cpg >= min_cpg

    @property
    def locus_id(self) -> str:
        return f"{self.chrom}:{self.start}-{self.end}:{self.strand}:{self.subfamily}"


def utr_window(start: int, end: int, strand: str, cons_start: int) -> tuple[int, int]:
    """Genomic interval covering the element's retained 5'UTR portion.

    Consensus coordinates run 5'->3' along the element, which for a '-' strand
    copy is right-to-left on the genome. ``cons_start`` is where this (often
    5'-truncated) copy picks up the consensus.
    """
    if cons_start >= UTR_LEN:
        return (0, 0)
    span = UTR_LEN - cons_start
    if strand == "+":
        return (start, min(start + span, end))
    return (max(end - span, start), end)


def promoter_cpgs(seq: str, start: int, end: int, strand: str, cons_start: int) -> list[int]:
    """Genomic starts of promoter CpGs, ordered 5'->3' along the element.

    A CpG is palindromic across strands, so the dinucleotide is found on the
    forward sequence either way; only the *ordering* depends on strand.
    """
    w0, w1 = utr_window(start, end, strand, cons_start)
    if w1 <= w0:
        return []
    window = seq[w0:w1].upper()
    hits = [w0 + i for i in range(len(window) - 1) if window[i] == "C" and window[i + 1] == "G"]
    if strand == "-":
        hits.reverse()
    return hits[:N_PROMOTER_CPG]


def build_catalog(rmsk_rows, seq: str, chrom: str = "chr22") -> list[L1Locus]:
    """Assemble loci from RepeatMasker rows plus reference sequence.

    ``rmsk_rows`` yields dicts with keys: genoStart, genoEnd, strand, repName,
    milliDiv, cons_start, flank_bp, cons_len.
    """
    out: list[L1Locus] = []
    for r in rmsk_rows:
        locus = L1Locus(
            chrom=chrom,
            start=int(r["genoStart"]),
            end=int(r["genoEnd"]),
            strand=str(r["strand"]),
            subfamily=str(r["repName"]),
            divergence=float(r["milliDiv"]) / 1000.0,
            cons_start=int(r["cons_start"]),
            flank_bp=int(r["flank_bp"]),
            cons_len=int(r.get("cons_len", 0)),
        )
        if locus.has_utr:
            locus.cpg_genomic = promoter_cpgs(
                seq, locus.start, locus.end, locus.strand, locus.cons_start
            )
        out.append(locus)
    return out
