#!/usr/bin/env python3
"""Independent outcome-free component replay. No numeric reporter file is accepted."""
from __future__ import annotations
import argparse,collections,csv,gzip,hashlib,json,pathlib,re,shutil,subprocess,sys,tempfile
sys.dont_write_bytecode=True
ROOT=pathlib.Path(__file__).resolve().parents[3]
HERE=pathlib.Path(__file__).resolve().parent
NEW=ROOT/'research/variant-transfer/source/design-mapping-v2/all-reference-validated-pairs.jsonl.gz'
EXPECTED_NEW='e37f6e4b621159fd432fb89d180a317bde145863e8bd706ba8f315c6f2531979'
PRIOR=ROOT/'research/variant-transfer/prior-exposure-inventory/inventory-v1'
EXPECTED_PRIOR='861d838582ce1f0d52d846e6d20ee5d6c0f28c059db257032c9cff7a1fae6faf'
CONTROL=ROOT/'research/variant-transfer/source/S30-technical-control-matches.json'
PROTOCOL=ROOT/'experiments/variant-transfer-v1.json'
EXPECTED_PROTOCOL='811930c2d3a5ee763014ea611acda4d2c84404b7d5ae399c7843c62d0a828149'
def hashfile(p):
 with pathlib.Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def require(c,m):
 if not c:raise ValueError(m)
def loads(s):
 def obj(pairs):
  d={}
  for k,v in pairs:require(k not in d,'duplicate JSON key');d[k]=v
  return d
 return json.loads(s,object_pairs_hook=obj,parse_constant=lambda x:(_ for _ in ()).throw(ValueError(x)))
def load(p):return loads(pathlib.Path(p).read_text())
def lines(p):
 with gzip.open(p,'rt') as f:
  for line in f:yield loads(line)
def dump(p,a):pathlib.Path(p).write_text(json.dumps(a,sort_keys=True,indent=2,allow_nan=False)+'\n')
def zlines(p,rows):
 with pathlib.Path(p).open('xb') as f:
  with gzip.GzipFile(filename='',mode='wb',fileobj=f,mtime=0) as z:
   for row in rows:z.write((json.dumps(row,sort_keys=True,separators=(',',':'))+'\n').encode())
def literal(s):return isinstance(s,str) and len(s)==200 and not(set(s)-set('ACGT'))
def reverse(s):return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
def prepare(out):
 require(not out.exists(),'output exists; choose a new directory');require(hashfile(NEW)==EXPECTED_NEW,'new DNA SHA mismatch');require(hashfile(PRIOR/'manifest.json')==EXPECTED_PRIOR,'prior manifest SHA mismatch');require(hashfile(PROTOCOL)==EXPECTED_PROTOCOL,'semantic protocol changed')
 manifest=load(PRIOR/'manifest.json')
 for name,item in manifest['files'].items():
  require(pathlib.Path(name).name==name,'unsafe prior path');p=PRIOR/name
  require(not p.is_symlink() and p.stat().st_size==item['bytes'] and hashfile(p)==item['sha256'],'prior inventory mismatch:'+name)
 controls=load(CONTROL);control_ids=sorted({r['id'] for r in controls['matched_external_variants']});require(len(control_ids)==41,'unexpected technical control ID count')
 out.mkdir(parents=True);scratch=pathlib.Path(tempfile.mkdtemp(prefix='helix-independent-groups-'))
 metadata=[];newpath=scratch/'new.tsv';ordered=sorted(lines(NEW),key=lambda q:q['id']);ids=set()
 with newpath.open('w') as f:
  for row in ordered:
   require(set(row)<={'id','source_row_1based','reasons','source_hg19','whole_insert_chain_candidates','hg38','sequence_ref','sequence_alt','gc_preserving','parent_id'},'outcome/unexpected new input field')
   identity=row['id'];require(identity not in ids and not any(c in identity for c in '\t\r\n'),'duplicate/unsafe variant ID');ids.add(identity)
   ref,alt=row['sequence_ref'],row['sequence_alt'];require(literal(ref) and literal(alt),'new DNA is not literal200nt');require([i for i,(a,b) in enumerate(zip(ref,alt)) if a!=b]==[99],'expected centered SNV')
   g=row['hg38'];require(re.fullmatch(r'chr(?:[1-9]|1[0-9]|2[0-2]|X|Y)',g['chromosome']) is not None,'nonprimary chromosome');require(type(g['start']) is int and type(g['end']) is int and 0<=g['start'] and g['end']-g['start']==200,'invalid hg38 interval')
   parent=row.get('parent_id',identity);require(not any(c in parent for c in '\t\r\n'),'unsafe parent ID')
   f.write('\t'.join([identity,g['chromosome'],str(g['start']),str(g['end']),parent,ref,alt])+'\n');metadata.append({'id':identity,'hg38':g})
 require(len(metadata)==332719,'full pre-exclusion population changed');del ordered
 oldpath=scratch/'old.tsv';guards=[];oldcounts={}
 with oldpath.open('w') as f:
  for name in ['agarwal-cores.jsonl.gz','kircher-deletion-windows200.jsonl.gz','hbb-auxiliary-windows200.jsonl.gz','kircher-reference-windows200.jsonl.gz']:
   count=0
   for row in lines(PRIOR/name):
    s=row['sequence'];require(literal(s),'prior DNA length/type');require(hashlib.sha256(s.encode()).hexdigest()==row['sequence_sha256'],'prior sequence SHA');require(hashlib.sha256(min(s,reverse(s)).encode()).hexdigest()==row['canonical_rc_sha256'],'prior RC SHA');require(not any(c in row['id'] for c in '\t\r\n'),'unsafe old ID')
    factor=name=='kircher-reference-windows200.jsonl.gz';mutations=[]
    if factor:
     positions=set()
     for p,alts in row['observed_snv_alts']:
      require(type(p) is int and 0<=p<200 and p not in positions and isinstance(alts,str) and alts and not(set(alts)-set('ACGT')) and s[p] not in alts,'invalid prior mutation factorization');positions.add(p);mutations.append(str(p)+':'+alts)
    f.write('\t'.join(['factor' if factor else 'explicit',row['id'],s,','.join(mutations)])+'\n');count+=1
    if name=='agarwal-cores.jsonl.gz':guards.append(('explored_Agarwal_insert',row))
   oldcounts[name]=count
 for name,category in [('kircher-loci.json','explored_Kircher_locus'),('hbb-auxiliary-loci.json','explored_HBB_locus')]:guards.extend((category,q) for q in load(PRIOR/name)['loci'])
 guardpath=scratch/'guards.tsv'
 with guardpath.open('w') as f:
  for category,g in guards:
   require(g['assembly']=='GRCh38' and g['coordinate_convention']=='zero_based_half_open','prior coordinate/build mismatch')
   f.write('\t'.join([category,g['chromosome'],str(g['start']),str(g['end']),g['id']])+'\n')
 controlpath=scratch/'controls.tsv';controlpath.write_text('\n'.join(control_ids)+'\n')
 executable=scratch/'independent-grouping';flags=['-O3','-std=c++20','-Wall','-Wextra','-pedantic'];compilation=subprocess.run(['g++',*flags,str(HERE/'group.cpp'),'-o',str(executable)],capture_output=True,text=True)
 (out/'compiler-output.txt').write_text(compilation.stdout+compilation.stderr);require(compilation.returncode==0,'C++ compilation failed; see compiler-output.txt')
 shutil.copyfile(HERE/'group.cpp',out/'implementation.cpp');shutil.copyfile(HERE/'run.py',out/'orchestration.py')
 compiler=subprocess.check_output(['g++','--version'],text=True)
 bindings={'schema_version':1,'phase':'prepared_before_independent_execution_no_response_file','protocol':{'path':str(PROTOCOL.relative_to(ROOT)),'sha256':EXPECTED_PROTOCOL},'new_DNA':{'path':str(NEW.relative_to(ROOT)),'sha256':EXPECTED_NEW,'count':len(metadata)},'prior_inventory':{'path':str(PRIOR.relative_to(ROOT)),'manifest_sha256':EXPECTED_PRIOR,'counts':oldcounts},'technical_controls':{'path':str(CONTROL.relative_to(ROOT)),'sha256':hashfile(CONTROL),'ids':control_ids,'scope':'Reported separately from prior-exposure flags; primary grouping source-cohort gate remains separate.'},'algorithm':{'grouping':'Exact source parent; floor(start/1Mb) bin keyed by chromosome; separate strict-positive-interval overlap; any 200nt allele/either strand Hamming<=20. Complete new-new components precede old-exposure exclusion.','new_new_index':'21 disjoint pieces, sorted flat C++ index, forward stored alleles and both query strands; verify every cross-component candidate by seven-word packed Hamming. Skip only already-connected components.','old_explicit':'21 pieces for <=20 mismatches; every retained allele and both strands screened.','old_factorized':'22 pieces for <=21 reference mismatches, exact allowed observed-SNV repair at21; <=20 reference matches directly.','genomic_guard':'Merge each collection/chromosome old closed interval expanded by10000bp; test closed intersection with each new half-open interval endpoints, reproducing inclusive gap<=10000 convention.','pruning':'Whole already-exposed components need no further old-hit enumeration. Old coordinates never union new components, and old nodes never receive1Mb bins.','candidate_budget':1000000000,'group_ID':'SHA256 of newline-joined lexicographically sorted variant IDs without trailingnewline'},'compiler':{'version':compiler,'flags':flags,'binary_sha256':hashfile(executable)},'scratch_inputs':{p.name:{'sha256':hashfile(p),'bytes':p.stat().st_size} for p in [newpath,oldpath,guardpath,controlpath]},'new_numeric_outcome_access':False}
 dump(out/'bindings.json',bindings);zlines(out/'new-metadata-index.jsonl.gz',metadata)
 return scratch,executable,newpath,oldpath,guardpath,controlpath

def execute(out):
 scratch,binary,new,old,guard,controls=prepare(out)
 try:
  with (out/'execution.log').open('w') as log:
   code=subprocess.call([str(binary),str(new),str(old),str(guard),str(controls),str(scratch/'partition.tsv'),str(scratch/'witnesses.tsv')],stdout=log,stderr=subprocess.STDOUT)
  require(code==0,'independent accelerator failed; preserve execution.log')
  with (scratch/'partition.tsv').open() as source:records=list(csv.DictReader(source,delimiter='\t'))
  require(len(records)==332719,'incomplete independent partition')
  members=collections.defaultdict(list)
  for row in records:members[row['root_index']].append(row['id'])
  groupids={root:hashlib.sha256('\n'.join(sorted(ids)).encode()).hexdigest() for root,ids in members.items()}
  ledger=[];components={}
  for row in records:
   group=groupids[row['root_index']];retained=row['retained_prior']=='1';reasons=[] if not row['reasons'] else row['reasons'].split(',');item={'id':row['id'],'component':group,'retained_after_prior_exposure':retained,'component_exclusion_reasons':reasons,'technical_control_component':row['technical_control_component']=='1'};ledger.append(item)
   if group not in components:components[group]={'component':group,'variant_ids':sorted(members[row['root_index']]),'size':len(members[row['root_index']]),'retained_after_prior_exposure':retained,'exclusion_reasons':reasons,'technical_control_component':item['technical_control_component']}
   else:require(components[group]['retained_after_prior_exposure']==retained and components[group]['exclusion_reasons']==reasons and components[group]['technical_control_component']==item['technical_control_component'],'inconsistent component flags')
  zlines(out/'variant-ledger.jsonl.gz',ledger);dump(out/'components.json',{'schema_version':1,'components':sorted(components.values(),key=lambda q:q['component'])})
  with (scratch/'witnesses.tsv').open('rb') as source,(out/'witnesses.tsv.gz').open('xb') as f:
   with gzip.GzipFile(filename='',mode='wb',fileobj=f,mtime=0) as target:shutil.copyfileobj(source,target)
  summary={'schema_version':1,'independent_component_reconstruction_complete':True,'new_numeric_outcome_access':False,'new_variants':len(ledger),'components':len(components),'retained_prior_variants':sum(q['retained_after_prior_exposure'] for q in ledger),'retained_prior_components':sum(q['retained_after_prior_exposure'] for q in components.values()),'technical_control_components':sum(q['technical_control_component'] for q in components.values()),'technical_control_still_retained_variants':sum(q['retained_after_prior_exposure'] and q['technical_control_component'] for q in ledger),'source_cohort_gate_passed':False,'comparison_to_primary':'pending separate command','retained_ID_order_sha256':hashlib.sha256('\n'.join(q['id'] for q in ledger if q['retained_after_prior_exposure']).encode()).hexdigest()};dump(out/'summary.json',summary)
 except Exception as error:
  dump(out/'failure.json',{'type':type(error).__name__,'error':str(error),'scratch_path':str(scratch),'new_numeric_outcome_access':False});raise
 finally:
  dump(out/'manifest.json',{'schema_version':1,'files':{p.name:{'bytes':p.stat().st_size,'sha256':hashfile(p)} for p in sorted(out.iterdir()) if p.name!='manifest.json' and p.is_file()}})
 print(json.dumps(summary,sort_keys=True))

def verify_manifest(directory):
 manifest=load(directory/'manifest.json');expected=set(manifest['files'])|{'manifest.json'}
 require(expected=={p.name for p in directory.iterdir() if p.is_file()},'unexpected/missing archived file')
 for name,item in manifest['files'].items():
  path=directory/name;require(pathlib.Path(name).name==name and not path.is_symlink(),'unsafe archived manifest path')
  require(path.stat().st_size==item['bytes'] and hashfile(path)==item['sha256'],'archive member size/SHA mismatch:'+name)
 return hashfile(directory/'manifest.json')

def compare(primary,independent,out):
 require(not out.exists(),'comparison output exists');verify_manifest(primary);verify_manifest(independent)
 require(load(primary/'summary.json')['grouping_and_prior_exposure_gates_passed'] is True,'primary grouping incomplete')
 require(load(independent/'summary.json')['independent_component_reconstruction_complete'] is True,'independent grouping incomplete')
 p=load(primary/'components.json')['components'];q=load(independent/'components.json')['components'];pm={v['component']:v for v in p};qm={v['component']:v for v in q};checks=[]
 def check(name,ok,detail=None):checks.append({'check':name,'passed':bool(ok),'detail':detail})
 primary_inputs=load(primary/'inputs.json');independent_inputs=load(independent/'bindings.json')
 check('same_new_DNA_input',primary_inputs['new_pairs']['sha256']==independent_inputs['new_DNA']['sha256'])
 check('same_prior_inventory',primary_inputs['prior_inventory']['manifest_sha256']==independent_inputs['prior_inventory']['manifest_sha256'])
 check('component_key_set',set(pm)==set(qm),{'primary':len(pm),'independent':len(qm)})
 for key in sorted(set(pm)|set(qm)):
  if key not in pm or key not in qm:continue
  for field in ['variant_ids','size','retained_after_prior_exposure','exclusion_reasons']:check('component_'+key+'_'+field,pm[key][field]==qm[key][field])
 a=list(lines(primary/'variant-ledger.jsonl.gz'));b=list(lines(independent/'variant-ledger.jsonl.gz'));am={r['id']:r for r in a};bm={r['id']:r for r in b};check('all_variant_IDs',set(am)==set(bm))
 mismatches=[]
 for identity in sorted(set(am)&set(bm)):
  if any(am[identity][field]!=bm[identity][field] for field in ['component','retained_after_prior_exposure','component_exclusion_reasons']):mismatches.append(identity)
 check('every_variant_partition_exclusion_and_reasons',not mismatches,{'mismatch_count':len(mismatches),'ids':mismatches})
 out.mkdir(parents=True);dump(out/'comparison.json',{'schema_version':1,'matched':all(c['passed'] for c in checks),'checks_passed':sum(c['passed'] for c in checks),'failed_checks':[c for c in checks if not c['passed']],'checks':checks,'bindings':{'primary_manifest_sha256':hashfile(primary/'manifest.json'),'independent_manifest_sha256':hashfile(independent/'manifest.json')},'scope':'Full outcome-free new-component and prior-exposure reconstruction. Excluded-component witness catalogs need not enumerate every old hit; exact component reason classes are compared. S30 technical-control flags retained separately until source-cohort binding.','new_numeric_outcome_access':False});shutil.copyfile(HERE/'run.py',out/'implementation.py');dump(out/'manifest.json',{'files':{p.name:{'bytes':p.stat().st_size,'sha256':hashfile(p)} for p in sorted(out.iterdir())}})
 print(json.dumps({'matched':all(c['passed'] for c in checks),'checks':len(checks),'failures':sum(not c['passed'] for c in checks)}));require(all(c['passed'] for c in checks),'independent grouping differs; retained full failures')

def main():
 p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True);r=sub.add_parser('run');r.add_argument('--out',type=pathlib.Path,required=True);c=sub.add_parser('compare');c.add_argument('--primary',type=pathlib.Path,required=True);c.add_argument('--independent',type=pathlib.Path,required=True);c.add_argument('--out',type=pathlib.Path,required=True);a=p.parse_args();execute(a.out) if a.command=='run' else compare(a.primary,a.independent,a.out)
if __name__=='__main__':main()
