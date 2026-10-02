"""Synthetic gates for the provenance-locked ten-genome extension."""
import gzip
import sys
from pathlib import Path

import numpy as np
import pytest

SCRIPTS = Path(__file__).resolve().parents[1]/'scripts'
sys.path.insert(0,str(SCRIPTS))
import analyze_ten_genome_mei as m  # noqa: E402
import mei_reference_opportunity as o  # noqa: E402


def engine():
    # Unit-test grouping using the committed snapshot, never the concurrent
    # dirty dedup. In a detached production checkout this falls back to itself.
    snapshot = Path('/tmp/rtm_head_10g/scripts')
    return m.load_engine(snapshot if snapshot.exists() else SCRIPTS)


def call(e, sample, pos, host=None, orient='+'):
    return e.Call(sample,pos,'chr1',pos,'ALU',orient,'unnested',{},sample,
                  same_family_host=host,match_host=host,source_ids=[f'{sample}:{pos}'])


def test_gap_subtraction_is_half_open_and_union_safe():
    assert o.subtract([(10,30)],[(0,12),(20,25),(24,28),(30,40)]) == [(12,20),(28,30)]
    assert o.subtract([(10,30)],[(10,30)]) == []
    assert o.subtract([(10,30)],[(0,10),(30,40)]) == [(10,30)]


def test_minus_host_opportunity_orients_from_host_five_prime():
    e=engine()
    h=e.Host('chr1',100,400,'-','AluY','ALU')
    assert o.host_segments(h,{'chr1':[(100,120)]}) == [(0,280)]
    bases=o.bin_opportunity([h],{'chr1':[(100,120)]},[0,20,280,300])
    assert bases.tolist()==[20,260,0]


def test_opportunity_is_reference_universe_not_event_hosts():
    e=engine()
    hosts=[e.Host('chr1',100,400,'+','AluY','ALU'),e.Host('chr1',500,800,'+','AluY','ALU')]
    assert o.bin_opportunity(hosts,{'chr1':[]},[0,120,140,300]).tolist()==[240,40,320]
    # Duplicate annotation is one physical host opportunity.
    assert np.array_equal(o.bin_opportunity(hosts,{'chr1':[]},[0,300]),o.bin_opportunity(hosts+hosts,{'chr1':[]},[0,300]))


def test_known_mei_exclusion_only_uses_known_breakpoint():
    e=engine()
    h=e.Host('chr1',100,400,'+','AluY','ALU')
    a=call(e,'HG03086',101,h)
    a.info={'KNOWNMEI':'True','GT':'0/0','GENE':'ignored'}
    assert o.known_exclusions({('chr1','ALU'):[h]},[a])=={h.key}
    a.pos=100
    assert o.known_exclusions({('chr1','ALU'):[h]},[a])==set()
    a.pos=400
    assert o.known_exclusions({('chr1','ALU'):[h]},[a])=={h.key}
    a.info['KNOWNMEI']='False'
    assert o.known_exclusions({('chr1','ALU'):[h]},[a])==set()


def test_ten_carrier_bitmap_and_private_layers():
    e=engine()
    cs=[call(e,s,200) for s in m.SAMPLES]
    sites,_=m.dedup(e,cs)
    row=m.rows_for_sites(e,sites)[0]
    assert row['n_carriers']==10 and row['presence_bitmap']=='1111111111'
    counts=m.layer_counts(sites)[0]
    assert (counts['private'],counts['union'],counts['shared'])==(0,1,1)
    assert counts['carrier_histogram']['10']==1


def test_repeated_grouping_does_not_mutate_provenance():
    e=engine()
    cs=[call(e,'HG03086',200),call(e,'HG03086',205),call(e,'HG01474',203)]
    results=[m.dedup(e,cs,10) for _ in range(3)]
    assert all(sum(c.source_count for c in ss[0].members)==3 for ss,_ in results)
    assert all(c.source_count==1 and len(c.source_ids)==1 for c in cs)


def test_frozen_no_chaining_and_orientation_host_keys():
    e=engine()
    ss,_=m.dedup(e,[call(e,'HG03086',100),call(e,'HG01474',110),call(e,'HG00171',120)])
    assert sorted(len(s.samples) for s in ss)==[1,2]
    ss,_=m.dedup(e,[call(e,'HG03086',100),call(e,'HG00171',100,orient='-')])
    assert len(ss)==2


def test_zero_opportunity_is_excluded_without_pseudocount():
    rows=m.profile_table([10,30],[0,20,40,60],[20,20,0],1000,'private','LINE1','consensus_primary')
    assert rows[-1]['zero_opportunity_excluded']
    assert rows[-1]['expected']==0 and rows[-1]['enrichment'] is None and rows[-1]['empirical_p_greater'] is None
    with pytest.raises(ValueError,match='zero-opportunity'):
        m.profile_table([50],[0,20,40,60],[20,20,0],1000,'private','LINE1','consensus_primary')


def test_real_alignment_projects_indel_and_truncated_slice():
    seq='ACGTCGATCGGATCCTAGCTAGGCTAACGT'
    host=seq[:10]+'AAA'+seq[10:]
    mapping,metrics=o.aligned_mapping(host,seq)
    assert mapping[20]==17
    assert metrics['max_internal_indel_bp']==3
    assert metrics['alignment_identity']==1


def test_repeatmasker_consensus_coordinates_are_zero_based(tmp_path):
    e=engine()
    plus=e.Host('chr1',100,110,'+','L1TEST','LINE1')
    minus=e.Host('chr1',200,210,'-','L1TEST','LINE1')
    rmsk=tmp_path/'rmsk.gz'
    with gzip.open(rmsk,'wt') as fh:
        fh.write('\t'.join(['0','0','0','0','0','chr1','100','110','0','+','L1TEST','LINE','L1','5','15','-5','1'])+'\n')
        fh.write('\t'.join(['0','0','0','0','0','chr1','200','210','0','-','L1TEST','LINE','L1','-5','15','5','2'])+'\n')
    cons=tmp_path/'cons.fa'
    cons.write_text('>L1TEST\n'+'A'*20+'\n')
    seq,meta,stats=o.consensus_metadata(cons,rmsk,{('chr1','LINE1'):[plus,minus]})
    assert meta[plus.key]==meta[minus.key]==('L1TEST',5,15)
    bases,_=o.consensus_opportunity([plus,minus],{'chr1':[]},meta,[0,5,10,15,20])
    assert bases.tolist()==[0,10,10,0]


def test_delta_computed_flags_more_than_twenty_percent():
    values={'enrichment_ALU':1,'enrichment_LINE1':1,'enrichment_SVA':1,'linker':1,'tail':1,
            'replication_fraction':1,'union':100,'private':60,'shared':40,'nested_L1':10,'l1_projected':10}
    delta=m.delta_rows({'headline':values},{'headline':{**values,'union':130}})
    row=next(r for r in delta if r['headline']=='Dedup union')
    assert row['flag_gt20pct'] and row['percent_change']==pytest.approx(30)


def test_n_gaps_exact_runs_across_chunk_edges(tmp_path):
    class E:
        PRIMARY_CHROMS=('chr1',)
        @staticmethod
        def load_fai(path):return {'chr1':None}, {'chr1':10}
        @staticmethod
        def fasta_sequence_chunk(fh,meta,start,end):return b'ACNNNTNNAC'[start:end]
    (tmp_path/'ref').write_bytes(b'ACNNNTNNAC')
    assert o.n_gaps(E,tmp_path/'ref',{'chr1':10}) == {'chr1':[(2,5),(6,8)]}


def test_unmasked_events_use_same_short_read_mask():
    e=engine()
    cs=[call(e,'HG03086',101),call(e,'HG00171',201)]
    assert m.unmasked_calls(cs,{'chr1':o.SpanMask([(100,101)])})==[cs[1]]
