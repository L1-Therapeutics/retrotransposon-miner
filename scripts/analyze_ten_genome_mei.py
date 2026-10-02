#!/usr/bin/env python3
"""Ten-genome MEI analysis; production must run from a committed detached snapshot.

All enrichment is conditional on reference-genome opportunity; reference host
spans are assumed present in the cohort, not measured as universally carried.
The existing dedup module supplies the frozen grouping algorithm unchanged.
"""
from __future__ import annotations

import argparse
import copy
import csv
import importlib.util
import json
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

import mei_reference_opportunity as opportunity

ROOT = Path('/Users/smyan/Desktop/Research/Research_Projects/L1')
SAMPLES = ('HG03086', 'HG01474', 'HG01566', 'HG03172', 'NA18498',
           'HG00171', 'HG01058', 'NA18939', 'NA19017', 'NA20845')
ORIGINAL = SAMPLES[:5]
FAMILIES = ('ALU', 'LINE1', 'SVA')
STATES = ('unnested', 'nested_sense', 'nested_antisense', 'nested_unknown')
SEED = 20261001
CAVEAT = opportunity.CAVEAT
INFO_FIELDS = {'MEIFAMILY', 'ORIENT', 'NESTED', 'KNOWNMEI', 'MEISUBFAMILY', 'TSD'}


def load_engine(directory):
    path = Path(directory)/'dedup_samples.py'
    spec = importlib.util.spec_from_file_location('ten_genome_frozen_dedup', path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    module.SAMPLES = SAMPLES  # runtime manifest only; no source or rules modified
    return module


def read_calls(engine, directory):
    """Read canonical, already-ingested inputs; never parse GT/GQ or gene INFO."""
    calls, qc = {}, []
    for sample in SAMPLES:
        path = directory/f'{sample}.vcf'
        rows, chroms, families, raw = [], Counter(), Counter(), Counter()
        with path.open() as fh:
            for line in fh:
                if line.startswith('#') or not line.strip():
                    continue
                f = line.rstrip().split('\t')
                info = {k: v for token in f[7].split(';')
                        for k, sep, v in [token.partition('=')] if sep and k in INFO_FIELDS}
                family = engine.normalize_family(info.get('MEIFAMILY', ''))
                if family not in FAMILIES:
                    raise engine.InputGateError(f'unsupported family: {sample}:{f[0]}:{f[1]}')
                rows.append(engine.Call(sample, len(rows), f[0], int(f[1]), family,
                                        info.get('ORIENT', ''), info.get('NESTED', ''), info, sample,
                                        source_ids=[f'{sample}:{f[2]}']))
                chroms[f[0]] += 1
                families[family] += 1
                raw[info.get('NESTED', '')] += 1
        observed = (len(rows), *(families[f] for f in FAMILIES))
        if sample in engine.BASELINES and observed != engine.BASELINES[sample]:
            raise engine.InputGateError(f'baseline mismatch {sample}: {observed}')
        if set(engine.PRIMARY_CHROMS)-set(chroms):
            raise engine.InputGateError(f'non-genome-wide {sample}')
        if sample == 'HG03086' and raw['nested'] != 261:
            raise engine.InputGateError('HG03086 raw nested gate !=261')
        if set(raw)-{'nested', 'unnested'}:
            raise engine.InputGateError(f'invalid binary NESTED: {sample}')
        qc.append({'sample': sample, 'calls': len(rows), **dict(families),
                   'chrY': chroms['chrY'], 'chromosomes': dict(chroms),
                   'raw_nested': raw['nested'], 'genome_wide': True, 'schema': 'PASS (ingest verified)',
                   'sha256': engine.sha256(path), 'note': '2011 calls: accepted +5.8% deviation' if sample == 'NA19017'
                   else ('SVA 3.85–4.73%: accepted cohort-consistent deviation' if sample not in ORIGINAL else '')})
        calls[sample] = [c for c in rows if c.chrom in engine.PRIMARY_CHROMS]
    return calls, qc


def dedup(engine, calls, window=10):
    """Copies isolate committed collapse provenance mutations across repeated runs."""
    sites, collapsed = engine.unique_sites(copy.deepcopy(calls), window)
    for site in sites:
        if len(site.members) != len(site.samples) or not 1 <= len(site.samples) <= 10:
            raise engine.InputGateError('one-carrier-per-sample gate failed')
        if any(abs(c.pos-site.anchor.pos) > window for c in site.members):
            raise engine.InputGateError('anchored non-transitive window gate failed')
    return sites, collapsed


def rows_for_sites(engine, sites):
    out = []
    for i, site in enumerate(sites, 1):
        row = engine.site_row(site, f'US{i:05d}', 10)
        host = site.representative.same_family_host
        row.update({'host_name': host.name if host else '', 'host_start0': host.start0 if host else '',
                    'host_end0': host.end0 if host else '', 'host_strand': host.strand if host else '',
                    'host_len': host.length if host else '',
                    'host_selection_rule': 'longest_containing_then_leftmost_name',
                    'presence_bitmap': ''.join('1' if s in site.samples else '0' for s in SAMPLES),
                    'opportunity_caveat': CAVEAT})
        row.update({f'present_{s}': int(s in site.samples) for s in SAMPLES})
        out.append(row)
    return out


def write_csv(path, rows):
    if not rows:
        raise ValueError(f'cannot infer empty table schema: {path}')
    with path.open('w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def coverage_data(engine, hosts, gaps, excluded):
    """Disjoint covered reference spans for original/Tier1/Tier2 probabilities."""
    result = {}
    for chrom in engine.PRIMARY_CHROMS:
        for family in FAMILIES:
            hs = hosts.get((chrom, family), [])
            original = engine.merge_intervals(hs)
            primary = opportunity.subtract(original, gaps[chrom])
            secondary = opportunity.subtract(engine.merge_intervals(h for h in hs if h.key not in excluded), gaps[chrom])
            result[chrom, family] = (original, primary, secondary)
    return result


def gc_frames(engine, fasta, lengths, gaps, coverage, targets):
    """200bp GC controls: prior tolerance/pool/filter conventions, N-masked bases.

    GC fraction retains the rtm-gc-null full-bin convention. Host coverage is
    divided by non-N opportunity for Tier1/2; legacy full-span controls are kept
    solely to reconstruct the previous five-genome headlines.
    """
    meta, _ = engine.load_fai(fasta)
    gc_parts, cov = [], {tier: {f: [] for f in FAMILIES} for tier in ('legacy', 'primary', 'tier2')}
    site_gc, stats = {}, Counter()
    target_bins = defaultdict(set)
    for call in targets:
        target_bins[call.chrom].add(call.pos0//200)
    with fasta.open('rb') as fh:
        for chrom in engine.PRIMARY_CHROMS:
            prefixes = {(family, tier): engine.interval_prefix(coverage[chrom, family][j])
                        for family in FAMILIES for j, tier in enumerate(cov)}
            for start in range(0, lengths[chrom], 200_000):
                end = min(start+200_000, lengths[chrom])
                seq = engine.fasta_sequence_chunk(fh, meta[chrom], start, end)
                arr = np.frombuffer(seq, dtype=np.uint8)
                nb = (len(arr)+199)//200
                padded = np.zeros(nb*200, dtype=np.uint8)
                padded[:len(arr)] = arr
                matrix = padded.reshape(nb, 200)
                spans = np.minimum(200, len(arr)-np.arange(nb)*200)
                gc = ((matrix == 71)|(matrix == 67)).sum(axis=1)/spans
                acgt = ((matrix == 65)|(matrix == 67)|(matrix == 71)|(matrix == 84)).sum(axis=1)
                non_n = spans-(matrix == 78).sum(axis=1)
                keep = (acgt/spans >= .5) & (non_n > 0)
                gc_parts.append(gc[keep].astype(np.float32))
                stats['total_bins'] += nb
                stats['kept_bins'] += int(keep.sum())
                starts = start+np.arange(nb)*200
                ends = starts+spans
                for (family, tier), prefix in prefixes.items():
                    covered = engine.cumulative_coverage(ends, *prefix)-engine.cumulative_coverage(starts, *prefix)
                    denominator = spans if tier == 'legacy' else non_n
                    cov[tier][family].append((covered[keep]/denominator[keep]).astype(np.float32))
                for target in target_bins[chrom]:
                    local = target-start//200
                    if 0 <= local < nb:
                        site_gc[chrom, target] = float(gc[local]) if keep[local] else None
            print(f'GC frame {chrom}', flush=True)
    gc = np.concatenate(gc_parts)
    order = np.argsort(gc)
    bundles = {tier: {'gc_by_site': site_gc, 'frame': {'gc': gc, 'sorted_order': order,
                 'coverage': {f: np.concatenate(parts) for f, parts in cs.items()}}} for tier, cs in cov.items()}
    return bundles, dict(stats)


def family_enrichment(engine, sites, lengths, coverage, bundles, gaps, excluded, reps, tier):
    result = []
    cov_index = {'legacy': 0, 'primary': 1, 'tier2': 2}[tier]
    uniform_by_chrom = {chrom: sum(b-a for a,b in coverage[chrom, 'ALU'][cov_index]) /
                        (lengths[chrom] if tier == 'legacy' else lengths[chrom]-sum(b-a for a,b in gaps[chrom]))
                        for chrom in engine.PRIMARY_CHROMS}
    for j, family in enumerate(FAMILIES):
        # Tier2 is a new eligible-host universe: remove observations assigned to
        # excluded hosts as well as removing those hosts from opportunity.
        sub = [s for s in sites if s.representative.family == family and
               (tier != 'tier2' or not s.representative.same_family_host or s.representative.same_family_host.key not in excluded)]
        obs = sum(s.representative.same_family_host is not None for s in sub)
        if family == 'ALU':
            probabilities = np.asarray([uniform_by_chrom[s.representative.chrom] for s in sub])
        else:
            probabilities = engine.gc_probabilities(sub, family, bundles[tier], np.random.default_rng(SEED+j))
        expected = float(probabilities.sum())
        null = engine.resample_nested(probabilities, reps, np.random.default_rng(SEED+100+j))
        result.append({'family': family, 'sites': len(sub), 'nested': obs,
                       'expected': expected, 'enrichment': obs/expected if expected else None,
                       'p': (1+int((null >= obs).sum()))/(reps+1), 'resamples': reps,
                       'null': 'uniform' if family == 'ALU' else 'GC-matched 200bp ±0.05; max2000',
                       'tier': tier, 'opportunity_caveat': CAVEAT})
    return result


def profile_table(offsets, edges, bases, reps, label, family, method):
    """Aggregate reference-base null; structural-zero bins never enter the null."""
    observed = np.histogram(offsets, bins=edges)[0]
    total = float(sum(bases))
    if offsets and total <= 0:
        raise ValueError(f'no opportunity for observed events: {method}')
    probabilities = np.asarray(bases)/total if total else np.zeros(len(bases))
    active = probabilities > 0
    if any(observed[~active]):
        raise ValueError('observed event in zero-opportunity bin')
    null = np.zeros((reps, len(bases)), dtype=int)
    if any(active):
        null[:, active] = np.random.default_rng(SEED).multinomial(len(offsets), probabilities[active]/probabilities[active].sum(), size=reps)
    return [{'cohort_layer': label, 'family': family, 'method': method,
             'bin_start': float(a), 'bin_end': float(b), 'observed': int(observed[i]),
             'opportunity_bp': float(bases[i]), 'expected': len(offsets)*float(probabilities[i]),
             'observed_per_opportunity_bp': float(observed[i]/bases[i]) if bases[i] else None,
             'enrichment': float(observed[i]/(len(offsets)*probabilities[i])) if len(offsets)*probabilities[i] else None,
             'empirical_p_greater': (1+int((null[:, i] >= observed[i]).sum()))/(reps+1) if active[i] else None,
             'zero_opportunity_excluded': not bool(active[i]), 'events': len(offsets), 'opportunity_caveat': CAVEAT}
            for i, (a, b) in enumerate(zip(edges, edges[1:]))]


def unmasked_calls(calls, masks):
    if not masks:
        return calls
    return [c for c in calls if not opportunity.overlap(masks.get(c.chrom, []), c.pos0, c.pos0+1)]


def alu_offsets(calls, excluded=frozenset()):
    return [c.offset for c in calls if c.family == 'ALU' and c.same_family_host is not None
            and 280 <= c.same_family_host.length <= 320 and c.same_family_host.key not in excluded
            and c.same_family_host.strand in {'+', '-'}]


def l1_event_projection(engine, calls, fasta, sequences, metadata):
    """Align event-bearing hosts once, then map all contained child breakpoints."""
    import pysam
    fa = pysam.FastaFile(str(fasta))
    cache, rows = {}, []
    complement = str.maketrans('ACGTN', 'TGCAN')
    for call in calls:
        h = call.same_family_host
        if call.family != 'LINE1' or h is None:
            continue
        row = {'sample': call.sample, 'source_index': call.source_index, 'chrom': call.chrom,
               'pos': call.pos, 'host_id': h.host_id, 'host_subfamily': h.name,
               'host_offset': call.offset, 'consensus_offset': None, 'status': 'no_exact_consensus_extent',
               'alignment_identity': None, 'max_internal_indel_bp': None, 'alignment_score': None,
               'opportunity_caveat': CAVEAT}
        if h.key in metadata and h.strand in {'+', '-'}:
            name, lo, hi = metadata[h.key]
            if h.key not in cache:
                seq = fa.fetch(h.chrom, h.start0, h.end0).upper()
                if h.strand == '-':
                    seq = seq.translate(complement)[::-1]
                mapping, metrics = opportunity.aligned_mapping(seq, sequences[name][lo:hi])
                cache[h.key] = (mapping, metrics, lo)
            mapping, metrics, lo = cache[h.key]
            row.update(metrics)
            coord = int(mapping[call.offset]) if call.offset is not None else -1
            row['status'] = ('low_identity' if metrics['alignment_identity'] < .65 else
                             'large_internal_indel' if metrics['max_internal_indel_bp'] > 500 else
                             'unaligned_breakpoint' if coord < 0 else 'eligible')
            if row['status'] == 'eligible':
                row['consensus_offset'] = coord+lo
        rows.append(row)
    fa.close()
    return {(r['sample'], r['source_index']): r for r in rows}, rows


def layer_counts(sites):
    layers = {'private': [s for s in sites if len(s.samples) == 1], 'union': sites}
    rows = []
    for family in FAMILIES:
        sub = [s for s in sites if s.representative.family == family]
        counts = Counter(len(s.samples) for s in sub)
        row = {'family': family, 'private': counts[1], 'union': len(sub), 'shared': sum(v for k, v in counts.items() if k >= 2),
               'carrier_histogram': {str(i): counts[i] for i in range(1, 11)}}
        for label, members in layers.items():
            row[f'{label}_states'] = dict(Counter(s.representative.nested_state for s in members if s.representative.family == family))
        rows.append(row)
    return rows


def cohort_analysis(engine, samples, callsets, lengths, coverage, bundles, gaps, excluded, alu_opps,
                    l1_opps, projections, hosts, reps, masks=None, other_opps=None):
    calls = [c for s in samples for c in callsets[s]]
    sites, collapsed = dedup(engine, calls)
    private = [s for s in sites if len(s.samples) == 1]
    enrich = {label: {tier: family_enrichment(engine, ss, lengths, coverage, bundles, gaps, excluded, reps, tier)
                      for tier in ('legacy', 'primary', 'tier2')} for label, ss in [('private', private), ('union', sites)]}
    profiles = []
    for label, ss in [('private', private), ('union', sites)]:
        reps_calls = [s.representative for s in ss]
        for tier in ('primary', 'tier2', 'short_read'):
            exc = excluded if tier == 'tier2' else frozenset()
            events = alu_offsets(unmasked_calls(reps_calls, masks) if tier == 'short_read' else reps_calls, exc)
            profiles.extend(profile_table(events, np.arange(0, 321, 20), alu_opps[tier], reps, label, 'ALU', tier))
        # Private/union descriptive same-family relative profiles also cover
        # SVA (no canonical full-length coordinate window specified for SVA).
        for family in ('ALU', 'SVA'):
            eligible = [c for c in reps_calls if c.family == family and c.same_family_host is not None
                        and c.offset is not None and c.same_family_host.strand in {'+', '-'}]
            bases = other_opps[family]
            profiles.extend(profile_table([c.offset/c.same_family_host.length for c in eligible], np.linspace(0, 1, 11),
                                          bases, reps, label, family, 'same_family_relative_all'))
        for method in ('consensus_primary', 'consensus_tier2', 'consensus_short_read', 'relative_all', 'near_full_relative'):
            selected = [c for c in reps_calls if c.family == 'LINE1' and c.same_family_host is not None]
            if method.startswith('consensus'):
                selected = [c for c in selected if projections[c.sample, c.source_index]['status'] == 'eligible' and
                            (method != 'consensus_tier2' or c.same_family_host.key not in excluded)]
                if method == 'consensus_short_read':
                    selected = unmasked_calls(selected, masks)
                offsets = [projections[c.sample, c.source_index]['consensus_offset'] for c in selected]
                edges, bases = l1_opps[method]
            else:
                if method == 'near_full_relative':
                    selected = [c for c in selected if 5500 <= c.same_family_host.length <= 7000]
                selected = [c for c in selected if c.offset is not None and c.same_family_host.strand in {'+', '-'}]
                offsets = [c.offset/c.same_family_host.length for c in selected]
                edges, bases = l1_opps[method]
            profiles.extend(profile_table(offsets, edges, bases, reps, label, 'LINE1', method))
    per_sample = []
    for sample in samples:
        sample_sites, _ = dedup(engine, callsets[sample])
        for tier in ('primary', 'tier2'):
            exc = excluded if tier == 'tier2' else frozenset()
            rows = profile_table(alu_offsets([s.representative for s in sample_sites], exc), np.arange(0, 321, 20),
                                 alu_opps[tier], reps, sample, 'ALU', tier)
            per_sample.extend(rows)
    primary_alu = [r for r in profiles if r['family'] == 'ALU' and r['cohort_layer'] == 'union' and r['method'] == 'primary']
    linker = next(r for r in primary_alu if r['bin_start'] == 120)
    tail = next(r for r in primary_alu if r['bin_start'] == 280)
    replication = sum(all(next(r for r in per_sample if r['cohort_layer'] == s and r['method'] == 'primary' and r['bin_start'] == b)['enrichment'] > 1
                          for b in (120, 280)) for s in samples)
    significant = sum(all(next(r for r in per_sample if r['cohort_layer'] == s and r['method'] == 'primary' and r['bin_start'] == b)['empirical_p_greater'] < .025
                          for b in (120, 280)) for s in samples)
    l1bins = [r for r in profiles if r['family'] == 'LINE1' and r['cohort_layer'] == 'union' and r['method'] == 'consensus_primary']
    legacy_profile = engine.position_profile(engine.profile_events([], sites, 'ALU', 'unique'), 'union', 'ALU', 'unique')
    headline = {'union': len(sites), 'private': len(private), 'shared': len(sites)-len(private),
                'nested_L1': sum(s.representative.family == 'LINE1' and s.representative.same_family_host is not None for s in sites),
                'linker': linker['enrichment'], 'tail': tail['enrichment'], 'replication': replication,
                'replication_fraction': replication/len(samples), 'significant_replication': significant,
                'l1_projected': l1bins[0]['events'],
                'legacy_linker': next(r['enrichment'] for r in legacy_profile if r['bin_start_bp'] == 120),
                'legacy_tail': next(r['enrichment'] for r in legacy_profile if r['bin_start_bp'] == 280)}
    headline.update({f'enrichment_{r["family"]}': r['enrichment'] for r in enrich['union']['primary']})
    return {'samples': samples, 'sites': sites, 'collapsed': collapsed, 'layers': layer_counts(sites), 'enrichment': enrich,
            'profiles': profiles, 'per_sample': per_sample, 'headline': headline,
            'orientations': engine.orientation_summary(sites)}


def delta_rows(five, ten):
    labels = {'enrichment_ALU': 'Alu same-family enrichment', 'enrichment_LINE1': 'L1 same-family enrichment',
              'enrichment_SVA': 'SVA same-family enrichment', 'linker': 'Alu linker enrichment', 'tail': 'Alu tail enrichment',
              'replication_fraction': 'Both windows elevated: genome fraction', 'union': 'Dedup union', 'private': 'Private',
              'shared': 'Shared', 'nested_L1': 'Nested L1', 'l1_projected': 'Consensus-projected eligible L1'}
    out = []
    for key, label in labels.items():
        a, b = five['headline'][key], ten['headline'][key]
        change = (b/a-1)*100 if a else None
        reason = ('Additional genomes reveal more distinct sites and carrier combinations.' if key in {'union','private','shared','nested_L1','l1_projected'} else
                  'Additional genomes change the breakpoint/host mix under the same reference null.')
        out.append({'headline': label, 'five': a, 'ten': b, 'percent_change': change,
                    'flag_gt20pct': change is not None and abs(change) > 20,
                    'reason': reason if change is not None and abs(change) > 20 else ''})
    return out


def fmt(x):
    if x is None:
        return 'NA'
    return f'{x:.4g}' if isinstance(x, float) else str(x)


def md_table(headers, rows):
    return ['| '+' | '.join(headers)+' |', '| '+' | '.join(['---']*len(headers))+' |'] + [
        '| '+' | '.join(fmt(v).replace('|', ',') for v in row)+' |' for row in rows]


def dependent_gate(scripts, outdir, callsets):
    """Probe committed consumers unchanged; never turn an incompatible join into zero events."""
    spec = importlib.util.spec_from_file_location('ten_genome_committed_common', scripts/'nested_multi_sample_common.py')
    common = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(common)
    sites, load = common.load_unique_sites(outdir/'unique_sites.csv')
    samples = sorted({s for site in sites for s in site['carriers']})
    attached = common.attach_call_details(sites, common.load_callsets(callsets, samples))
    join = common.verify_join(attached)
    result = {'provenance': 'committed-code snapshot; no dirty-tree fallback', 'load': load, 'join': join,
              'status': 'incompatible' if join['verdict'] != 'join_consistent_with_dedup_output' else 'ready',
              'opportunity_caveat': CAVEAT}
    if result['status'] == 'incompatible':
        result['reason'] = ('Committed consumers require every recovered source call to have raw NESTED=nested; '
                            'the producer independently assigns all same-family RMSK overlaps, including raw unnested antisense calls. '
                            'This is a nesting-definition/input-contract mismatch, not zero recurrence. Consumers were not patched.')
    else:
        # Use existing CLI unmodified, temporary output subdirectory owned here.
        dest = outdir/'phase4_run'
        dest.mkdir(exist_ok=False)
        for script in ('joint_enrichment.py', 'recurrence_test.py'):
            command = [sys.executable, str(scripts/script), '--unique-sites', str(outdir/'unique_sites.csv'),
                       '--callset-dir', str(callsets), '--outdir', str(dest)]
            if script == 'joint_enrichment.py':
                command += ['--rmsk', str(ROOT/'nested_analysis/data/rmsk.txt.gz')]
            completed = subprocess.run(command, capture_output=True, text=True)
            result[script] = {'returncode': completed.returncode, 'stdout': completed.stdout, 'stderr': completed.stderr}
            if completed.returncode:
                result['status'] = 'failed'
                break
        for path in dest.iterdir():
            path.rename(outdir/f'phase4_{path.name}')
        dest.rmdir()
    (outdir/'phase4_status.json').write_text(json.dumps(result, indent=2)+'\n')
    return result


def provenance(scripts, engine, args):
    sha = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=scripts.parent, text=True).strip()
    hashes = {p.name: engine.sha256(p) for p in [scripts/'analyze_ten_genome_mei.py', scripts/'mei_reference_opportunity.py',
                scripts/'dedup_samples.py', scripts/'nested_multi_sample_common.py', scripts/'joint_enrichment.py', scripts/'recurrence_test.py']}
    return {'commit': sha, 'scripts_sha256': hashes, 'python': sys.version, 'python_executable': sys.executable,
            'replicates': args.replicates, 'seed': SEED, 'gap_source': str(args.fasta)+' (exact N/n runs)',
            'opportunity_caveat': CAVEAT}


def write_report(outdir, five, ten, qc, deltas, provenance_data, diagnostics, carriage, dependent):
    lines = ['# Ten-genome MEI analysis — private sites lead', '',
             '> William: "look at only unique variants not shared between samples."', '',
             f'All enrichment values are **{CAVEAT}**. This is a reference-coordinate null, not measured carriage.', '',
             '## QC and ingestion', '']
    lines += md_table(['Sample','Calls','Alu','L1','SVA','chrY','Genome-wide','Schema','Note'],
                      [[q['sample'],q['calls'],q['ALU'],q['LINE1'],q['SVA'],q['chrY'],'PASS',q['schema'],q['note']] for q in qc])
    lines += ['', 'Folder/ZIP/canonical checksums and normalized headers passed during ingestion. Existing copies were preserved; five manifest lines appended.',
              '', '## Layer 1 — PRIVATE sites (primary interpretation)', '']
    lines += md_table(['Family','Private','Unnested','Sense','Antisense','Unknown'],
                      [[r['family'],r['private'],*[r['private_states'].get(s,0) for s in STATES]] for r in ten['layers']])
    lines += ['', 'Private means carrier count exactly one in this ten-genome panel, not a claim of population rarity or de novo origin.',
              'Private orientation and position results are primary; union results are secondary context. Shared sites receive counts/histograms only.', '',
              '### Private same-family enrichment', '']
    lines += md_table(['Family','Observed nested','Expected','Enrichment','MC p','Sensitivity enrichment'],
                      [[a['family'],a['nested'],a['expected'],a['enrichment'],a['p'],b['enrichment']]
                       for a,b in zip(ten['enrichment']['private']['primary'],ten['enrichment']['private']['tier2'])])
    lines += ['', '## Computed original-five versus ten delta', '']
    lines += md_table(['Headline','Five (this code)','Ten (same code)','Change %','>20% flag / reason'],
                      [[r['headline'],r['five'],r['ten'],r['percent_change'],r['reason']] for r in deltas])
    lines += ['', f"Both-window descriptive elevation: {five['headline']['replication']}/5 versus {ten['headline']['replication']}/10; "
              f"both-window p<0.025: {five['headline']['significant_replication']}/5 versus {ten['headline']['significant_replication']}/10 (two windows, no across-genome adjustment).",
              'These are reference-opportunity-conditioned per-sample results, not replication at demonstrated-carried/callable hosts.',
              '', 'Historical five-genome numbers (1.67/1.51/12.8; linker 2.76; tail 3.86; 5/5; 4998/3051/1947; nested L1 231) used a different position denominator. '
              'The primary delta above uses the current reference-universe null in BOTH columns; the historical event-host-conditioned calculation is reproduced separately below.', '']
    lines += md_table(['Legacy reconstructed five headline','Value'],
                      [[r['family']+' nesting',r['enrichment']] for r in five['enrichment']['union']['legacy']] +
                      [['Linker (event-host-conditioned)',five['headline']['legacy_linker']],['Tail (event-host-conditioned)',five['headline']['legacy_tail']]])
    lines += ['', '## Layer 2 — dedup union / Layer 3 — shared counts only','']
    lines += md_table(['Family','Union','Private','Shared','Carrier histogram 1–10'],
                      [[r['family'],r['union'],r['private'],r['shared'],r['carrier_histogram']] for r in ten['layers']])
    lines += ['', '### Union same-family enrichment (secondary)', '']
    lines += md_table(['Family','Observed nested','Expected','Enrichment','MC p','Sensitivity enrichment'],
                      [[a['family'],a['nested'],a['expected'],a['enrichment'],a['p'],b['enrichment']]
                       for a,b in zip(ten['enrichment']['union']['primary'],ten['enrichment']['union']['tier2'])])
    lines += ['', '## Alu positions — near-full 280–320bp hosts; half-open 20bp bins','']
    other = [r for r in ten['profiles'] if r['method']=='same_family_relative_all']
    lines += md_table(['Layer','Family','Relative decile','Observed','Reference bp','O/E'],
                      [[r['cohort_layer'],r['family'],f"{r['bin_start']:g}–{r['bin_end']:g}",r['observed'],r['opportunity_bp'],r['enrichment']] for r in other])
    lines += ['', 'All host-relative views include truncation; they are secondary to the specified Alu and L1 primary coordinates.', '']
    alu = [r for r in ten['profiles'] if r['family']=='ALU' and r['method']=='primary']
    lines += md_table(['Layer','Bin','Observed','Reference bp','Expected','Enrichment','MC p'],
                      [[r['cohort_layer'],f"{r['bin_start']:g}–{r['bin_end']:g}",r['observed'],r['opportunity_bp'],r['expected'],r['enrichment'],r['empirical_p_greater']] for r in alu])
    lines += ['', '### Per-genome linker / tail', '']
    lines += md_table(['Sample','Window','Observed','Expected','Enrichment','MC p'],
                      [[r['cohort_layer'],'linker' if r['bin_start']==120 else 'tail',r['observed'],r['expected'],r['enrichment'],r['empirical_p_greater']]
                       for r in ten['per_sample'] if r['method']=='primary' and r['bin_start'] in (120,280)])
    lines += ['', '## L1-into-L1 — primary consensus-projected profile', '',
              'Host sequences are oriented 5′→3′ on the HOST strand and globally aligned to the exact subfamily consensus slice named by RepeatMasker. '
              'No nearest-name consensus substitution is used. Identity <0.65, internal indel >500bp, or unaligned breakpoints are excluded and audited.',
              'Opportunity uses ALL eligible reference L1 annotations, not event-bearing hosts: their RepeatMasker aligned subfamily extents, gap-masked. '
              'The repeat-coordinate extents are measured alignments, but exact per-base RMSK CIGAR is unavailable; gap/mappability subsegments are allocated affinely within each aligned extent. '
              'This aligned-extent approximation is explicit, not a claim of exact basewise opportunity.', '']
    l1 = [r for r in ten['profiles'] if r['family']=='LINE1' and r['method']=='consensus_primary']
    lines += md_table(['Layer','Consensus bin bp','Observed','Opportunity bp','Observed/opportunity','Expected','O/E','MC p'],
                      [[r['cohort_layer'],f"{r['bin_start']:g}–{r['bin_end']:g}",r['observed'],r['opportunity_bp'],r['observed_per_opportunity_bp'],r['expected'],r['enrichment'],r['empirical_p_greater']] for r in l1])
    for label in ('private','union'):
        sub = [r for r in l1 if r['cohort_layer']==label and not r['zero_opportunity_excluded']]
        peak = max(sub,key=lambda r: r['enrichment'] or 0) if sub else None
        lines += ['', f"{label}: {sub[0]['events'] if sub else 0} eligible projected L1-into-L1 sites; largest descriptive ratio "
                  f"{fmt(peak['enrichment']) if peak else 'NA'}x at {fmt(peak['bin_start']) if peak else 'NA'}–{fmt(peak['bin_end']) if peak else 'NA'}bp. "
                  f"{sum(r['expected']<5 for r in sub)}/{len(sub)} active bins have expected counts <5. Sparse bins and multiple comparisons limit shape/domain claims; "
                  'per-bin p-values are exploratory, not a formal powered domain test. Zero-opportunity bins are excluded, never pseudocounted.']
    lines += ['', '### Secondary L1 position sensitivities','',
              'Relative-position deciles include truncated hosts, so truncation changes which consensus regions are represented and biases this view. '
              'Near-full-only sensitivity uses 5500–7000bp reference hosts; this does not establish intact ORFs or activity.', '']
    sens = [r for r in ten['profiles'] if r['family']=='LINE1' and r['method']!='consensus_primary']
    lines += md_table(['Layer','Method','Bin','Observed','Expected','O/E'],
                      [[r['cohort_layer'],r['method'],f"{r['bin_start']:g}–{r['bin_end']:g}",r['observed'],r['expected'],r['enrichment']] for r in sens])
    lines += ['', '### Nested L1 states and child-family composition of L1 hosts', '']
    for layer in ten['layers']:
        if layer['family']=='LINE1':
            lines.append(f"- Private states: {layer['private_states']}; union states: {layer['union_states']}.")
    lines += [f"- Child-family composition / L1 child contexts: {diagnostics['child_contexts']}.",
              f"- L1-in-Alu is secondary ONLY: {diagnostics['l1_in_alu']} union L1 child sites overlap an Alu host; this never substitutes for L1-into-L1.",
              '- Literature novelty and Jedlicka 2019 plant-LTR domain precedent are operator-supplied context, not an independently verified exhaustive literature claim.',
              '', '## Opportunity specification and descriptive host evidence', '',
              f"- Exact gap source: {provenance_data['gap_source']}; primary masks N/n only, not short-read mappability.",
              f"- {CAVEAT}.",
              f"- Tier2 is labelled **{opportunity.TIER2}**: exclude any reference host containing the breakpoint of a KNOWNMEI=True symbolic insertion in these canonical callsets. "
              'Symbolic insertions consume no reference span; overlap uses [POS−1,POS). The remainder is NOT established fixed or universally carried; bias direction is unknown.',
              '- A parent host observed in a sample’s own calls provides descriptive positive context only. Using those counts as opportunity denominators would condition on the calls being tested and is circular.',
              '- Absence of a child call never establishes host absence or a callable negative. A genome lacking a host is ineligible, not a negative; host absence is unmeasured here.',
              f"- Short-read sensitivity: {diagnostics['mappability_path']}; mask Umap k100<0.5 only for that sensitivity, never the primary null.",
              f"- Reference diagnostics: {diagnostics['opportunity']}", '',
              '### Hosts seen as parent elements — descriptive, NOT carriage denominators','']
    lines += md_table(['Family','Window','Host id','Samples with parent evidence','Number of genomes'],
                      [[r['family'],r['window'],r['host_id'],r['samples'],r['n_genomes']] for r in carriage])
    lines += ['', '## Dependent joint-signature / recurrence analyses', '', f"Status: {dependent['status']}. {dependent.get('reason','')}",
              f"Committed-code probe: {dependent.get('join',{})}. No incompatible result is reported as zero.", '',
              '## William’s three answers (private-first; union separately labelled)', '']
    private_enrich = ten['enrichment']['private']['primary'][0]
    union_enrich = ten['enrichment']['union']['primary'][0]
    al = next(r for r in ten['layers'] if r['family']=='ALU')
    ps = al['private_states']
    resolved = ps.get('nested_sense',0)+ps.get('nested_antisense',0)
    private_bins = {r['bin_start']:r for r in alu if r['cohort_layer']=='private'}
    lines += [f"1. Into Alu? Private Alu sites show {fmt(private_enrich['enrichment'])}x same-family nesting enrichment "
              f"(MC p={fmt(private_enrich['p'])}; union secondary {fmt(union_enrich['enrichment'])}x), {CAVEAT}.",
              f"2. Same orientation? {ps.get('nested_sense',0)}/{resolved} resolved private nested Alu sites are sense "
              f"({100*ps.get('nested_sense',0)/resolved if resolved else 0:.1f}%); unknown={ps.get('nested_unknown',0)} is excluded, not antisense.",
              f"3. At the polyA? Private Alu tail 280–300bp has {private_bins[280]['observed']} sites and {fmt(private_bins[280]['enrichment'])}x "
              f"versus linker 120–140bp {private_bins[120]['observed']} sites and {fmt(private_bins[120]['enrichment'])}x, {CAVEAT}; the tail window is a positional proxy, not a measured host polyA sequence.",
              '', '## Standing caveats', '',
              '- GT/GQ are never parsed or used: headers explicitly lack validated MEI genotyping. Gene/snpEff fields are not analysis inputs or matching keys.',
              '- chrY counts appear in QC only; no chrY matching, enrichment, or profiles. chr22_mei.vcf is an excluded HG03086 slice.',
              '- Source NESTED is binary. Independent RepeatMasker overlap + strand derives unnested/nested_sense/nested_antisense/nested_unknown; unknown never becomes antisense.',
              '- HG03086 raw 261 is reproduced; the self-insertion report’s 320 uses different host-selection rules and remains a distinct historical count. Fresh reannotation is not required to equal either.',
              '- Frozen grouping: chromosome+family+literal insertion orientation+same-family host; ±10bp anchor window, no transitive chaining, one carrier per sample. Exact/±5/±20 sweep persisted.',
              '- Per-site union representatives count once; per-genome tables collapse within-sample records. Every rerun copies call objects to avoid provenance accumulation.',
              '- Reference opportunity is an assumption, not individualized carriage or detection probability. Shared/private classifications depend on this ten-genome panel and detection.',
              '', '## Provenance', '', '```json', json.dumps(provenance_data,indent=2), '```','']
    (outdir/'analysis_summary.md').write_text('\n'.join(lines))


def run(args):
    scripts = Path(__file__).resolve().parent
    engine = load_engine(scripts)
    calls, qc = read_calls(engine, args.callset_dir)
    all_calls = [c for s in SAMPLES for c in calls[s]]
    print('Reading reference hosts', flush=True)
    hosts = engine.read_rmsk(args.rmsk)
    engine.assign_hosts(all_calls, hosts)
    _, lengths = engine.load_fai(args.fasta)
    lengths = {c:lengths[c] for c in engine.PRIMARY_CHROMS}
    print('Deriving N gaps from FASTA', flush=True)
    gaps = opportunity.n_gaps(engine, args.fasta, lengths)
    excluded = opportunity.known_exclusions(hosts, all_calls)
    cov = coverage_data(engine, hosts, gaps, excluded)
    print('Building GC controls', flush=True)
    bundles, gc_stats = gc_frames(engine, args.fasta, lengths, gaps, cov, all_calls)
    print('Reading consensus extents', flush=True)
    seqs, metadata, consensus_stats = opportunity.consensus_metadata(args.consensus, args.rmsk, hosts)
    l1hosts = [h for (c,f),hs in hosts.items() if f=='LINE1' for h in hs]
    aluhosts = [h for (c,f),hs in hosts.items() if f=='ALU' for h in hs if 280<=h.length<=320]
    print('Reading short-read sensitivity mask', flush=True)
    masks = opportunity.read_low_mappability(args.low_mappability, engine.PRIMARY_CHROMS)
    alu_opps = {'primary': opportunity.bin_opportunity(aluhosts,gaps,np.arange(0,321,20)),
                'tier2': opportunity.bin_opportunity(aluhosts,gaps,np.arange(0,321,20),excluded),
                'short_read': opportunity.bin_opportunity(aluhosts,gaps,np.arange(0,321,20),masks=masks)}
    maximum = max((hi for _,lo,hi in metadata.values()),default=6000)
    edges = np.arange(0, ((maximum+499)//500)*500+1,500)
    l1_opps, masked_hosts = {}, {}
    for name, ex, mask in [('consensus_primary',frozenset(),None),('consensus_tier2',excluded,None),('consensus_short_read',frozenset(),masks)]:
        bases, nmask = opportunity.consensus_opportunity(l1hosts,gaps,metadata,edges,ex,mask)
        l1_opps[name] = (edges,bases)
        masked_hosts[name] = nmask
    l1_opps['relative_all'] = (np.linspace(0,1,11),opportunity.bin_opportunity(l1hosts,gaps,np.linspace(0,1,11),relative=True))
    l1_opps['near_full_relative'] = (np.linspace(0,1,11),opportunity.bin_opportunity([h for h in l1hosts if 5500<=h.length<=7000],gaps,np.linspace(0,1,11),relative=True))
    other_opps = {family: opportunity.bin_opportunity([h for (chrom,f),hs in hosts.items() if f==family for h in hs],
                  gaps,np.linspace(0,1,11),relative=True) for family in ('ALU','SVA')}
    print('Aligning L1 event hosts', flush=True)
    projections, projection_rows = l1_event_projection(engine, all_calls, args.fasta, seqs, metadata)
    print('Analyzing original five and ten with identical code', flush=True)
    five = cohort_analysis(engine,ORIGINAL,calls,lengths,cov,bundles,gaps,excluded,alu_opps,l1_opps,projections,hosts,args.replicates,masks,other_opps)
    ten = cohort_analysis(engine,SAMPLES,calls,lengths,cov,bundles,gaps,excluded,alu_opps,l1_opps,projections,hosts,args.replicates,masks,other_opps)
    print('Five/ten analysis completed; checking registered matching counts', flush=True)
    # Registered count gates are matching checks, not expectations of the new null.
    if (five['headline']['union'],five['headline']['private'],five['headline']['shared'],five['headline']['nested_L1']) != (4998,3051,1947,231):
        raise engine.InputGateError(f'five-genome matching reproduction failed: {five["headline"]}')
    sweep = []
    for w in (0,5,10,20):
        ss, collapsed = dedup(engine, all_calls,w)
        sweep.append({'window':w,'union':len(ss),'private':sum(len(s.samples)==1 for s in ss),
                      'shared':sum(len(s.samples)>=2 for s in ss),'within_sample_collapsed':collapsed})
    site_rows = rows_for_sites(engine,ten['sites'])
    private_rows = [r for r in site_rows if r['n_carriers']==1]
    args.outdir.mkdir(parents=True,exist_ok=True)
    write_csv(args.outdir/'unique_sites.csv',site_rows)
    write_csv(args.outdir/'private_sites.csv',private_rows)
    per_call = []
    for sample in SAMPLES:
        ss,_ = dedup(engine,calls[sample])
        for row in rows_for_sites(engine,ss):
            row['sample'] = sample
            per_call.append(row)
    write_csv(args.outdir/'per_sample_calls.csv',per_call)
    write_csv(args.outdir/'position_profiles.csv',ten['profiles']+ten['per_sample'])
    write_csv(args.outdir/'l1_consensus_projection.csv',projection_rows)
    accumulation = []
    for i,sample in enumerate(SAMPLES,1):
        ss,_ = dedup(engine,[c for s in SAMPLES[:i] for c in calls[s]])
        row={'samples_included':i,'sample_added':sample,'unique_sites':len(ss)}
        row.update({f:sum(s.representative.family==f for s in ss) for f in FAMILIES})
        row['opportunity_caveat']=CAVEAT
        accumulation.append(row)
    write_csv(args.outdir/'accumulation_curve.csv',accumulation)
    carriage_groups = defaultdict(set)
    for call in all_calls:
        h=call.same_family_host
        if h and call.offset is not None:
            for lo,hi,label in [(120,140,'linker'),(280,300,'tail')]:
                if call.family=='ALU' and 280<=h.length<=320 and lo<=call.offset<hi:
                    carriage_groups[call.family,label,h.host_id].add(call.sample)
            if call.family=='LINE1':
                carriage_groups[call.family,'any L1-host position',h.host_id].add(call.sample)
    carriage = [{'family':f,'window':w,'host_id':h,'samples':','.join(s for s in SAMPLES if s in ss),
                 'n_genomes':len(ss),'opportunity_caveat':CAVEAT} for (f,w,h),ss in sorted(carriage_groups.items())]
    write_csv(args.outdir/'host_parent_evidence.csv',carriage)
    contexts=Counter()
    unique_contexts=Counter()
    representatives=[s.representative for s in ten['sites']]
    # Child composition of all L1 hosts uses independent ANY-child overlap.
    for chrom in engine.PRIMARY_CHROMS:
        cc=[c for c in all_calls if c.chrom==chrom]
        assigned=engine._sweep_assign(cc,hosts.get((chrom,'LINE1'),[]))
        for c in cc:
            if id(c) in assigned:
                contexts[c.family]+=1
        unique_calls=[c for c in representatives if c.chrom==chrom]
        unique_assigned=engine._sweep_assign(unique_calls,hosts.get((chrom,'LINE1'),[]))
        for c in unique_calls:
            if id(c) in unique_assigned:
                unique_contexts[c.family]+=1
    lc=[s.representative for s in ten['sites'] if s.representative.family=='LINE1']
    ctx={'per_call_children_in_L1_hosts':dict(contexts),'union_children_in_L1_hosts':dict(unique_contexts),
         'union_L1_child_context':dict(Counter('L1 host' if c.same_family_host else 'Alu host' if c.alu_host else 'neither L1 nor Alu' for c in lc))}
    diagnostics={'opportunity':{'N_gap_bp':{c:sum(b-a for a,b in gs) for c,gs in gaps.items()},
                 'excluded_known_mei_hosts':len(excluded),'gc':gc_stats,'consensus':consensus_stats,
                 'consensus_masked_host_annotations':masked_hosts,'projection_status':dict(Counter(r['status'] for r in projection_rows))},
                 'mappability_path':str(args.low_mappability) if masks is not None else 'not available',
                 'child_contexts':ctx,'l1_in_alu':sum(c.alu_host is not None for c in lc)}
    prov=provenance(scripts,engine,args)
    deltas=delta_rows(five,ten)
    # Persist full five rerun numeric substrate in the task-owned audit JSON.
    audit={'provenance':prov,'qc':qc,'sweep':sweep,'five':{k:v for k,v in five.items() if k!='sites'},
           'ten':{k:v for k,v in ten.items() if k!='sites'},'delta':deltas,'diagnostics':diagnostics}
    (args.outdir/'analysis_metrics.json').write_text(json.dumps(audit,indent=2,default=lambda x: x.item() if isinstance(x,np.generic) else str(x))+'\n')
    lines=['# Matching audit','',CAVEAT,'','Frozen chromosome+family+orientation+same-family-host ±10bp; no transitive chaining; one carrier per sample.','']
    lines+=md_table(['Window bp','Union','Private','Shared','Same-sample collapse'],[[r[k] for k in ('window','union','private','shared','within_sample_collapsed')] for r in sweep])
    lines+=['','Exclusions: chrY (HG03172=1, NA18498=1, all others=0); chr22_mei.vcf slice excluded; unknown remains separate.',
            'Call objects copied before every grouping pass; matching rules are committed dedup code, configured with the ten-sample manifest only.',
            'Schema extended additively with numeric host geometry and ten-sample bitmaps; source IDs prefixed by sample.','',
            '## Every matched site decision','']
    lines+=md_table(['Site','Anchor breakpoint','Representative breakpoint','Carriers','Host'],
                   [[r['site_id'],s.anchor.pos,r['representative_pos'],r['samples'],r['same_family_host_id'] or 'no host'] for r,s in zip(site_rows,ten['sites'])])
    (args.outdir/'matching_audit.md').write_text('\n'.join(lines)+'\n')
    dependent=dependent_gate(scripts,args.outdir,args.callset_dir)
    write_report(args.outdir,five,ten,qc,deltas,prov,diagnostics,carriage,dependent)
    print(json.dumps({'five':five['headline'],'ten':ten['headline'],'layers':ten['layers'],'dependent':dependent},indent=2),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--callset-dir',type=Path,default=ROOT/'nested_analysis/callsets')
    p.add_argument('--outdir',type=Path,default=ROOT/'nested_analysis/results_multi_sample')
    p.add_argument('--rmsk',type=Path,default=ROOT/'nested_analysis/data/rmsk.txt.gz')
    p.add_argument('--fasta',type=Path,default=Path('/Users/smyan/retrotransposon-workdir/data/public/reference/hg38/Homo_sapiens_assembly38.fasta'))
    p.add_argument('--consensus',type=Path,default=Path('/Users/smyan/retrotransposon-workdir/data/public/retrotransposon_db/ucsc_repeatbrowser/hg38reps.fa'))
    p.add_argument('--low-mappability',type=Path,default=Path('/Users/smyan/retrotransposon-workdir/data/public/annotation/hg38/mappability/k100.Umap.MultiTrackMappability.low_lt0.5.bed'))
    p.add_argument('--replicates',type=int,default=10000)
    args=p.parse_args()
    if args.replicates<1000:
        p.error('at least 1000 resamples required')
    run(args)


if __name__=='__main__':
    main()
