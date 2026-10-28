#!/usr/bin/env python3
"""Outcome-free exhaustive pair geometry for fixed canonical trimer models.

Run from repository root: python3 research/variant-locality/pair-geometry-v1/implementation.py --out /tmp/pair-geometry-replay
Output directory must not exist. No experimental response is read.
"""
import argparse
import collections
import csv
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path

BASES = 'ACGT'
COMP = str.maketrans('ACGT', 'TGCA')
def rc(s): return s.translate(COMP)[::-1]
WORDS = [''.join(t) for t in itertools.product(BASES, repeat=3)]
CANON = sorted({min(w, rc(w)) for w in WORDS})
LOOKUP = {w:i for i,w in enumerate(CANON)}
def canon(w): return min(w, rc(w))
def all_counts(seq):
    v = [0] * 33
    for i in range(len(seq) - 2): v[LOOKUP[canon(seq[i:i+3])]] += 1
    v[-1] = sum(b in 'GC' for b in seq)
    return v

def substitute(seq, index, b): return seq[:index] + b + seq[index+1:]
def full_interaction(ref,p,q,a,b):
    A = substitute(ref,p,a)
    B = substitute(ref,q,b)
    AB = substitute(A,q,b)
    vectors = [all_counts(z) for z in (AB,ref,A,B)]
    return tuple(w+x-y-z for w,x,y,z in zip(*vectors))
def local_interaction(ref,p,q,a,b):
    # Direct per-window inclusion-exclusion of the at most two shared windows.
    out = [0]*33
    A = substitute(ref,p,a); B = substitute(ref,q,b); AB = substitute(A,q,b)
    for i in range(len(ref)-2):
        if i <= p <= i+2 and i <= q <= i+2:
            for seq, sign in ((AB,1),(ref,1),(A,-1),(B,-1)):
                out[LOOKUP[canon(seq[i:i+3])]] += sign
    return tuple(out)
def digest(path): return {'bytes':path.stat().st_size,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--out',required=True,type=Path); args=parser.parse_args()
    if args.out.exists(): raise SystemExit('Refusing to overwrite existing results')
    source = Path('research/variant-locality/model-envelope-v1/saved-fits.json')
    fits=json.loads(source.read_text())
    assert fits['feature_names'] == ['rc3mer_'+s for s in CANON] + ['gc_fraction']
    assert len(fits['raw_slopes']) == 5
    models=fits['models']
    # Treat each saved binary64 coefficient as its exact rational value.
    slopes = {m:[sum((Fraction.from_float(float(f[m][i])) for f in fits['raw_slopes']),Fraction())/5/198 for i in range(32)] for m in models}
    rows=[]; hist={}; max_info={}; checks=collections.Counter()
    for d in (1,2):
        size=5-d # d=1: four bases, d=2: three
        p=2-d; q=2
        vectors=collections.Counter(); l1=collections.Counter(); mags={m:Fraction() for m in models}; extrema={m:[] for m in models}
        for chars in itertools.product(BASES, repeat=size):
            ref=''.join(chars)
            for a in BASES:
                if a == ref[p]: continue
                for b in BASES:
                    if b == ref[q]: continue
                    v=local_interaction(ref,p,q,a,b)
                    assert v==full_interaction(ref,p,q,a,b)
                    assert not v[-1] and sum(v[:32])==0
                    A=substitute(ref,p,a); B=substitute(ref,q,b); AB=substitute(A,q,b)
                    # Reverse orientation and both edit directions preserve E;
                    # reversing either one changes its sign.
                    assert v==local_interaction(rc(ref),size-1-q,size-1-p,b.translate(COMP),a.translate(COMP))
                    assert v==local_interaction(AB,p,q,ref[p],ref[q])
                    assert tuple(-x for x in v)==local_interaction(A,p,q,ref[p],b)
                    checks['local_full_rc_and_direction_classes']+=1
                    # All external two-base left/right flanks; full count implementation.
                    for lf in map(''.join,itertools.product(BASES,repeat=2)):
                        for rf in map(''.join,itertools.product(BASES,repeat=2)):
                            assert v==full_interaction(lf+ref+rf,p+2,q+2,a,b)
                            checks['full_count_flank_comparisons']+=1
                    # One literal 200-nt full sequence per class.
                    reference='A'*97+ref+'C'*(103-size)
                    assert len(reference)==200
                    assert v==full_interaction(reference,p+97,q+97,a,b)
                    checks['full_200nt_comparisons']+=1
                    vectors[v]+=1; l1[sum(abs(x) for x in v)]+=1
                    values={m:sum((slopes[m][i]*v[i] for i in range(32)),Fraction()) for m in models}
                    for m,val in values.items():
                        if abs(val)>mags[m]: mags[m]=abs(val); extrema[m]=[]
                        if abs(val)==mags[m]: extrema[m].append({'reference_context':ref,'alt_first':a,'alt_second':b,'signed_value':float(val)})
                    rows.append([d,ref,p,q,a,b,AB,','.join(map(str,v)),sum(map(abs,v)),*[format(float(values[m]),'.17g') for m in models]])
        hist[str(d)]={'minimal_oriented_classes':sum(vectors.values()),'distinct_interaction_vectors':len(vectors),'zero_feature_interactions':vectors[(0,)*33], 'count_l1_histogram':dict(sorted(l1.items())),'vector_multiplicity_histogram':dict(sorted(collections.Counter(vectors.values()).items()))}
        max_info[str(d)]={m:{'absolute_max':float(mags[m]),'exact_max_numerator':str(mags[m].numerator),'exact_max_denominator':str(mags[m].denominator),'ties':extrema[m]} for m in models}
    assert hist['1']['minimal_oriented_classes']==2304 and hist['2']['minimal_oriented_classes']==576
    assert hist['1']['zero_feature_interactions']==hist['2']['zero_feature_interactions']==0
    args.out.mkdir(parents=True)
    with (args.out/'all-pairs.csv').open('w',newline='') as f:
        w=csv.writer(f,lineterminator='\n'); w.writerow(['distance_nt','reference_minimal_context','first_index','second_index','alt_first','alt_second','double_context','interaction_count_vector','count_l1',*models]); w.writerows(rows)
    summary={'description':'Exhaustive oriented minimal shared-window contexts for interior pairs, fixed-200nt features. Input edits are real substitutions at both sites.', 'prediction_units':'Interaction in original predicted common reporter activity units. Fractions use exact saved binary64 raw slopes, mean over five fits, trimer denominator 198; no refitting.','canonical_trimers':CANON,'source_model':{'path':str(source),**digest(source)},'outcomes_read':False,'models_fitted':False,'model_numbers_are_not_biological_effect_bounds':True,'pair_summary':hist,'prediction_extrema':max_info,'checks':dict(checks),'caveats':['The nonzero result concerns the feature vector, not necessarily a particular fitted prediction; exact cancellation with any coefficients is possible.','The zero-feature result is verified exhaustively for interior edits; boundary pairs may have fewer shared windows.','Maximum predictions concern these fixed fitted models, not biological interactions or measured reporter expression.','Fixed affine transformations preserve the feature mixed difference up to a linear map; GC and intercept cancel.']}
    (args.out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    (args.out/'implementation.py').write_bytes(Path(__file__).read_bytes())
    (args.out/'manifest.json').write_text(json.dumps({'files':{p.name:digest(p) for p in sorted(args.out.iterdir())}},indent=2)+'\n')
    print(json.dumps({'pair_summary':hist,'extrema':{d:{m:x['absolute_max'] for m,x in vals.items()} for d,vals in max_info.items()},'checks':dict(checks)},indent=2))
if __name__=='__main__': main()
