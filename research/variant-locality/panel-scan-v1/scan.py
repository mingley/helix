#!/usr/bin/env python3
"""Outcome-free exact feature scan; uses only panel record id and sequence."""
import argparse, collections, csv, hashlib, itertools, json, pathlib, statistics
BASES = 'ACGT'
COMP = str.maketrans(BASES, 'TGCA')
EXPECTED_PANEL_SHA = 'fce6b2a1aaf5aa4d530a2e298dccb840bc7c594a05b6831fdb0c6c0f9e2a9be3'
def rc(s): return s.translate(COMP)[::-1]
def canon(s): return min(s,rc(s))
KMERS = sorted({canon(''.join(x)) for x in itertools.product(BASES,repeat=3)})
KINDEX = {k:i for i,k in enumerate(KMERS)}
def features(s):
    v=[0]*33
    for i in range(len(s)-2): v[KINDEX[canon(s[i:i+3])]]+=1
    v[32]=sum(c in 'GC' for c in s)
    return tuple(v)
def delta(s,a):
    t=s[:2]+a+s[3:]
    x,y=features(s),features(t)
    return tuple(v-u for u,v in zip(x,y))
def quantile(a,q):
    t=(len(a)-1)*q; lo=int(t); hi=min(lo+1,len(a)-1)
    return a[lo]+(a[hi]-a[lo])*(t-lo)
def digest(b): return hashlib.sha256(b).hexdigest()
def write_json(path,x): path.write_text(json.dumps(x,indent=2,sort_keys=True)+'\n')
def write_csv(path,headers,rows):
    with path.open('w',newline='') as f:
        w=csv.writer(f);w.writerow(headers);w.writerows(rows)
def main():
    p=argparse.ArgumentParser();p.add_argument('--panel',required=True);p.add_argument('--out',required=True);a=p.parse_args()
    out=pathlib.Path(a.out);out.mkdir(parents=True,exist_ok=True)
    if (out/'summary.json').exists(): raise SystemExit('refusing existing completed output')
    script_bytes=pathlib.Path(__file__).read_bytes()
    (out/'scan.py').write_bytes(script_bytes)
    raw=pathlib.Path(a.panel).read_bytes();assert digest(raw)==EXPECTED_PANEL_SHA
    source=json.loads(raw)
    # Discard all other fields immediately; numerical outcomes are never selected or used.
    records=[(r['id'],r['sequence']) for r in source['records']];del source
    assert len(records)==55276 and len({k for k,s in records})==55276
    assert all(len(s)==200 and set(s)<=set(BASES) for k,s in records)
    templates={}
    for z in itertools.product(BASES,repeat=5):
        s=''.join(z)
        for alt in BASES:
            if alt!=s[2]: templates[(s,alt)]=delta(s,alt)
    assert len(templates)==3072 and len(KMERS)==32
    zero={key for key,d in templates.items() if not any(d)}
    theoretical={(s,alt) for s,alt in templates if alt==s[2].translate(COMP) and s[0]==s[4].translate(COMP) and s[1]==s[3].translate(COMP)}
    assert zero==theoretical and len(zero)==64
    checks=collections.Counter()
    for (s,alt),d in templates.items():
        swap=(s[:2]+alt+s[3:],s[2]);reverse=(rc(s),alt.translate(COMP))
        assert templates[swap]==tuple(-v for v in d);checks['allele_swap_negates_exact_delta']+=1
        assert templates[reverse]==d;checks['reverse_complement_preserves_exact_delta']+=1
        assert (swap in zero)==((s,alt) in zero);checks['allele_swap_preserves_zero_status']+=1
        assert (reverse in zero)==((s,alt) in zero);checks['reverse_complement_preserves_zero_status']+=1
    zero_alternates={s:sum((s,b) in zero for b in BASES if b!=s[2]) for s,_ in templates}
    assert set(zero_alternates.values())=={0,1}
    totals=collections.Counter();first={};per_element=[];hist=collections.Counter();sequence_hashes=set();rc_identities=set()
    for element_id,s in records:
        windows=collections.Counter(s[i-2:i+3] for i in range(2,198))
        assert sum(windows.values())==196
        blind=sum(n*zero_alternates[k] for k,n in windows.items())
        hist[blind]+=1;per_element.append((element_id,blind,588,blind/588))
        totals.update(windows)
        sequence_hashes.add(digest(s.encode()));rc_identities.add(digest(min(s,rc(s)).encode()))
        for i in range(2,198):
            k=s[i-2:i+3]
            if k not in first:first[k]=(element_id,i,s)
    assert sum(totals.values())==55276*196
    class_rows=[];example_rows=[];observed_signatures=set();observed_classes=0;direct_checks=0
    for (s,alt),d in sorted(templates.items()):
        n=totals[s];example=first.get(s)
        if n:
            observed_classes+=1;observed_signatures.add(d)
            element_id,i,full=example
            changed=full[:i]+alt+full[i+1:]
            assert tuple(v-u for u,v in zip(features(full),features(changed)))==d
            direct_checks+=1
        else: element_id,i,full='',-1,''
        class_rows.append([s,s[2],alt,n,int(not any(d)),*d,rc(s),alt.translate(COMP),s[:2]+alt+s[3:],s[2],element_id,i,i+1 if i>=0 else -1])
        if not any(d) and n:
            assert features(full)==features(full[:i]+alt+full[i+1:])
            example_rows.append([element_id,i,i+1,s[2],alt,s,rc(s),digest(full.encode()),n])
    blind_total=sum(totals[s] for s,alt in zero)
    opportunity_total=sum(r[3] for r in class_rows)
    assert opportunity_total==32502288 and sum(x[1] for x in per_element)==blind_total
    counts=sorted(x[1] for x in per_element)
    write_csv(out/'opportunity-classes.csv',['reference_5mer','reference_base','alternate_base','panel_slot_count','exact_zero_delta',*[k+'_count_delta' for k in KMERS],'GC_numerator_delta','RC_reference_5mer','RC_alternate_base','swapped_reference_5mer','swapped_alternate_base','example_element_id','example_core_position_0based','example_core_position_1based'],class_rows)
    write_csv(out/'per-element.csv',['element_id','exact_zero_delta_interior_SNVs','all_interior_SNVs','exact_zero_fraction'],per_element)
    write_csv(out/'zero-feature-examples.csv',['element_id','core_position_0based','core_position_1based','reference_base','alternate_base','reference_5mer','RC_reference_5mer','reference_200nt_sha256','panel_class_opportunity_count'],example_rows)
    checks['direct_full_200nt_mutation_feature_recomputations']=direct_checks
    checks['exact_zero_class_characterization_all_directed_templates']=len(templates)
    checks['per_element_196_positions_times_3_alternates']=len(records)
    write_json(out/'checks.json',{'passed':True,'checks':dict(checks),'feature_vector':'32 exact reverse-complement-canonical overlapping 3mer counts followed by exact GC numerator; normalization fixed at198 and200 respectively for200nt sequence.','excluded_inputs':'No observations, reporter activities, source effects, model coefficients, fitted predictions or clinical labels are selected or used.'})
    write_json(out/'summary.json',{
      'schema_version':1,'purpose':'Outcome-free exact feature equivalence scan, not a biological neutrality test.',
      'source':{'path':'datasets/context-panel/panel.json','bytes':len(raw),'sha256':digest(raw),'selected_fields':['records[].id','records[].sequence']},
      'panel_slots':len(records),'length_nt':200,'unique_exact_sequence_identities':len(sequence_hashes),'unique_RC_sequence_identities':len(rc_identities),
      'position_definition':{'zero_based_inclusive':[2,197],'one_based_inclusive':[3,198],'positions_per_element':196,'alternates_per_position':3},
      'interior_SNV_opportunities':opportunity_total,'exact_zero_feature_delta_opportunities':blind_total,'exact_zero_feature_delta_fraction':blind_total/opportunity_total,
      'exact_oriented_5mer_plus_alternate_classes':len(templates),'observed_oriented_classes':observed_classes,'distinct_exact_delta_vectors':len(set(templates.values())),'observed_distinct_exact_delta_vectors':len(observed_signatures),'exact_zero_oriented_classes':len(zero),
      'per_element_exact_zero_count':{'min':min(counts),'q025_type7':quantile(counts,.025),'q25_type7':quantile(counts,.25),'median':statistics.median(counts),'mean':statistics.mean(counts),'q75_type7':quantile(counts,.75),'q975_type7':quantile(counts,.975),'max':max(counts),'elements_with_any':sum(x>0 for x in counts),'elements_with_none':sum(x==0 for x in counts),'histogram':dict(sorted(hist.items()))},
      'exact_zero_characterization':'All64andonly64classes have alternate=complement(reference_center), outerleft=complement(outerright), and innerleft=complement(innerright). Equivalently the mutant5mer equals the reverse complement of the reference5mer. The GCnumerator and canonical3mer spectrum are unchanged.',
      'limits':['Slots count every pre-existing assay-panel record. Duplicate or RC-identical sequences are not deduplicated and no genomic-unique count is claimed.','Exact feature equality forces equal outputs only for deterministic predictors using solely this feature vector. No measured functional effect, clinical consequence or biological neutrality is inferred.','Only interior single-nucleotide substitutions are enumerated. Boundary positions, indels, multi-edit haplotypes, positional or longer-context features are outside this scan.','Reported quantiles describe this fixed sequence panel; they are not inferential intervals.'],
      'implementation_sha256':digest(pathlib.Path(__file__).read_bytes())})
    artifact_names=['scan.py','summary.json','checks.json','opportunity-classes.csv','per-element.csv','zero-feature-examples.csv']
    files={name:{'bytes':(out/name).stat().st_size,'sha256':digest((out/name).read_bytes())} for name in artifact_names}
    write_json(out/'manifest.json',{'schema_version':1,'files':files})
    print(json.dumps({'panel_slots':len(records),'opportunities':opportunity_total,'zero_opportunities':blind_total,'zero_fraction':blind_total/opportunity_total,'observed_classes':observed_classes,'distinct_delta_vectors':len(observed_signatures),'per_element':{k:v for k,v in json.loads((out/'summary.json').read_text())['per_element_exact_zero_count'].items() if k!='histogram'},'checks':dict(checks)},indent=2))
if __name__=='__main__':main()
