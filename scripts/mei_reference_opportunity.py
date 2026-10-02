"""Reference-coordinate opportunity for the ten-genome MEI analysis.

Enrichment is conditional on reference-genome opportunity; reference host spans
are assumed present in the cohort. This is not measured sample carriage.
"""
from __future__ import annotations

import bisect
import gzip
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from Bio import Align, SeqIO

CAVEAT = "conditional on reference-genome opportunity; reference host spans are assumed present in the cohort"
TIER2 = "excluding known-MEI-overlapping hosts"


def merge(spans):
    out = []
    for lo, hi in sorted(spans):
        if hi <= lo:
            continue
        if out and lo <= out[-1][1]:
            out[-1] = (out[-1][0], max(hi, out[-1][1]))
        else:
            out.append((lo, hi))
    return out


class SpanMask(list):
    """Disjoint interval union with cached ends for repeated host queries."""
    def __init__(self, spans):
        super().__init__(merge(spans))
        self.ends = [hi for _, hi in self]


def subtract(spans, masks):
    """Subtract a sorted union of half-open masks from half-open spans."""
    masks = masks if isinstance(masks, SpanMask) else SpanMask(masks)
    ends = masks.ends
    out = []
    for lo, hi in spans:
        cursor = lo
        j = bisect.bisect_right(ends, lo)
        while j < len(masks) and masks[j][0] < hi:
            a, b = masks[j]
            if a > cursor:
                out.append((cursor, min(a, hi)))
            cursor = max(cursor, b)
            if cursor >= hi:
                break
            j += 1
        if cursor < hi:
            out.append((cursor, hi))
    return out


def overlap(spans, lo, hi):
    """Total overlap; spans must already be a disjoint sorted union."""
    ends = spans.ends if isinstance(spans, SpanMask) else [b for _, b in spans]
    j = bisect.bisect_right(ends, lo)
    total = 0
    while j < len(spans) and spans[j][0] < hi:
        a, b = spans[j]
        total += max(0, min(b, hi) - max(a, lo))
        j += 1
    return total


def runs(mask, offset=0):
    edges = np.diff(np.r_[False, mask, False].astype(np.int8))
    return [(int(a + offset), int(b + offset)) for a, b in zip(np.flatnonzero(edges == 1), np.flatnonzero(edges == -1))]


def n_gaps(engine, fasta, lengths):
    """Exact N/n runs from the indexed hg38 FASTA, not a mappability proxy."""
    meta, _ = engine.load_fai(fasta)
    out = {}
    with fasta.open('rb') as fh:
        for chrom in engine.PRIMARY_CHROMS:
            spans = []
            for start in range(0, lengths[chrom], 2_000_000):
                seq = engine.fasta_sequence_chunk(fh, meta[chrom], start, min(start + 2_000_000, lengths[chrom]))
                spans.extend(runs(np.frombuffer(seq, dtype=np.uint8) == ord('N'), start))
            out[chrom] = SpanMask(spans)
    return out


def known_exclusions(hosts, calls):
    """Exclude hosts containing a known-MEI insertion breakpoint, independently of nesting.

    These are symbolic insertions: END=POS and they consume no reference span.
    Their overlap interval is therefore [POS-1, POS), not WINDOWSTART/WINDOWEND.
    """
    points = defaultdict(set)
    for call in calls:
        if call.info.get('KNOWNMEI', '').lower() == 'true':
            points[call.chrom].add(call.pos0)
    points = {c: sorted(p) for c, p in points.items()}
    excluded = set()
    for hs in hosts.values():
        for h in hs:
            p = points.get(h.chrom, [])
            i = bisect.bisect_left(p, h.start0)
            if i < len(p) and p[i] < h.end0:
                excluded.add(h.key)
    return excluded


def host_segments(host, gaps, masks=None):
    spans = subtract([(host.start0, host.end0)], gaps.get(host.chrom, []))
    if masks:
        spans = subtract(spans, masks.get(host.chrom, []))
    return [(a-host.start0, b-host.start0) if host.strand == '+' else
            (host.end0-b, host.end0-a) for a, b in spans]


def bin_opportunity(hosts, gaps, edges, excluded=frozenset(), masks=None, relative=False):
    """Host-coordinate base opportunity across the FULL reference host universe.

    Each physical annotation contributes once, independently of observed calls.
    Relative bins may have fractional bp because edges divide host length.
    """
    result = np.zeros(len(edges)-1, dtype=float)
    seen = set()
    for host in hosts:
        if host.key in seen or host.key in excluded or host.strand not in {'+', '-'}:
            continue
        seen.add(host.key)
        spans = host_segments(host, gaps, masks)
        for i, (a, b) in enumerate(zip(edges, edges[1:])):
            lo, hi = (a*host.length, b*host.length) if relative else (a, b)
            result[i] += sum(max(0, min(hi, y)-max(lo, x)) for x, y in spans)
    return result


def consensus_metadata(path, engine, rmsk):
    """Exact subfamily matches and UCSC RepeatMasker aligned consensus extents.

    UCSC rmsk repStart/repEnd use zero-based half-open repeat coordinates.
    Minus-strand rows store repLeft/repEnd/repStart in columns 13/14/15.
    Extents (not hypothetical full-length consensuses) define the primary
    consensus opportunity, as in the project's aligned-extent Phase 1 null.
    """
    sequences = {r.id: str(r.seq).upper() for r in SeqIO.parse(path, 'fasta')}
    wanted = {h.key for hs in rmsk.values() for h in hs if h.family == 'LINE1'}
    metadata = {}
    census = Counter()
    chroms = rmsk_chroms(rmsk)
    with gzip.open(engine, 'rt') as fh:
        for line in fh:
            f = line.rstrip().split('\t')
            if len(f) < 17 or f[5] not in chroms:
                continue
            key = (f[5], int(f[6]), int(f[7]), f[10], f[9])
            if key not in wanted or key in metadata:
                continue
            name = f[10]
            if name not in sequences:
                census['no_exact_subfamily_consensus'] += 1
                continue
            lo = int(f[13] if f[9] == '+' else f[15])
            hi = int(f[14])
            if not 0 <= lo < hi <= len(sequences[name]):
                census['invalid_consensus_extent'] += 1
                continue
            metadata[key] = (name, lo, hi)
            census['eligible_reference_annotations'] += 1
    return sequences, metadata, dict(census)


def rmsk_chroms(rmsk):
    return {c for c, _ in rmsk}


def aligned_mapping(host_sequence, consensus_sequence):
    """Affine-gap global mapping in host 5'-to-3' coordinates.

    Target is the host's own RMSK-aligned subfamily slice: terminal truncation
    does not count as an internal indel. Score convention matches Phase 0.
    """
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(host_sequence, consensus_sequence)[0]
    mapping = np.full(len(host_sequence), -1, dtype=np.int32)
    hb, cb = alignment.aligned
    for (a, b), (c, d) in zip(hb, cb):
        if b-a != d-c:
            raise ValueError('unequal ungapped alignment blocks')
        mapping[a:b] = np.arange(c, d)
    valid = np.flatnonzero(mapping >= 0)
    identity = sum(host_sequence[i] == consensus_sequence[mapping[i]] for i in valid)/len(valid) if len(valid) else 0
    indel = 0
    for i in range(len(hb)-1):
        indel = max(indel, int(hb[i+1][0]-hb[i][1]), int(cb[i+1][0]-cb[i][1]))
    return mapping, {'alignment_identity': identity, 'max_internal_indel_bp': indel,
                     'alignment_score': float(alignment.score)}


def consensus_opportunity(hosts, gaps, metadata, edges, excluded=frozenset(), masks=None):
    """RMSK aligned-extent opportunity; zero bins are left zero, not pseudocounted.

    Extents are already aligned to the named subfamily. If gaps/masks intersect
    an annotation, allocate only its remaining reference fraction to its aligned
    extent; this explicit extent-level approximation is reported. Non-event
    hosts are NOT selected by child calls or sample carriage.
    """
    out = np.zeros(len(edges)-1, dtype=float)
    masked = 0
    seen = set()
    for host in hosts:
        if host.key in excluded or host.key in seen or host.key not in metadata:
            continue
        seen.add(host.key)
        _, lo, hi = metadata[host.key]
        segments = host_segments(host, gaps, masks)
        available = sum(b-a for a, b in segments)
        if available < host.length:
            masked += 1
        # Map callable segment edges affinely onto the aligned extent. This is
        # the aligned-extent null, not a claim of an exact per-base RMSK CIGAR.
        for a, b in segments:
            x = lo + a/host.length*(hi-lo)
            y = lo + b/host.length*(hi-lo)
            first = max(0, int(np.searchsorted(edges, x, side='right')-1))
            last = min(len(out)-1, int(np.searchsorted(edges, y, side='left')))
            for i in range(first, last+1):
                out[i] += max(0, min(y, edges[i+1])-max(x, edges[i]))
    return out, masked


def read_low_mappability(path, chroms):
    out = defaultdict(list)
    if not path or not Path(path).is_file():
        return None
    with Path(path).open() as fh:
        for line in fh:
            f = line.split('\t')
            if len(f) >= 3 and f[0] in chroms:
                out[f[0]].append((int(f[1]), int(f[2])))
    return {c: SpanMask(v) for c, v in out.items()}
