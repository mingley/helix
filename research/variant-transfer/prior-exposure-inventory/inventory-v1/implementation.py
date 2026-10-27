#!/usr/bin/env python3
"""Freeze old explored DNA/locus identities without reading new-study outcomes.

Old CSV identity fields and original author allele fields are read; numeric old
score/count columns are ignored. No Siraj file is accepted by this program.
"""
from __future__ import annotations
import argparse,collections,csv,gzip,hashlib,io,json,pathlib,re,shutil,sys
sys.dont_write_bytecode=True
ROOT=pathlib.Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'scripts'))
from twobit_reference import TwoBitReader
REFERENCE_SHA='1f67aaa17a77b327738fe750ab37430a85dddf73c7d9189385ad1839259acec0'
PANEL_SHA='fce6b2a1aaf5aa4d530a2e298dccb840bc7c594a05b6831fdb0c6c0f9e2a9be3'
PATTERN=re.compile(r'n\.(\d+)(?:([ACGT])>([ACGT])|=)')
COMPLEMENT=str.maketrans('ACGT','TGCA')
AUTHOR={'00000015-a-1':'GRCh38_F9.tsv','00000016-a-1':'GRCh38_FOXE1.tsv','00000017-a-1':'GRCh38_GP1BA.tsv','00000018-a-1':'GRCh38_HBB.tsv','00000019-a-1':'GRCh38_HBG1.tsv','00000020-a-1':'GRCh38_HNF4A.tsv','00000023-a-1':'GRCh38_LDLR.tsv','00000023-a-2':'GRCh38_LDLR.2.tsv','00000024-a-1':'GRCh38_MSMB.tsv','00000027-a-1':'GRCh38_PKLR-24h.tsv','00000027-b-1':'GRCh38_PKLR-48h.tsv','00000031-a-1':'GRCh38_TERT-HEK.tsv','00000031-b-1':'GRCh38_TERT-GBM.tsv','00000031-c-1':'GRCh38_TERT-GAa.tsv','00000031-d-1':'GRCh38_TERT-GSc.tsv'}

def require(ok,message):
 if not ok:raise ValueError(message)
def sha(data):return hashlib.sha256(data).hexdigest()
def filehash(path):
 h=hashlib.sha256()
 with pathlib.Path(path).open('rb') as f:
  for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
 return h.hexdigest()
def rc(s):return s.translate(COMPLEMENT)[::-1]
def dna(s):require(bool(s) and not(set(s)-set('ACGT')),'invalid literal DNA');return s
def seqmeta(s):
 dna(s);return {'sequence':s,'length':len(s),'sequence_sha256':sha(s.encode()),'canonical_rc_sha256':sha(min(s,rc(s)).encode())}
def read_fasta(path):return dna(''.join(line.strip() for line in pathlib.Path(path).read_text().splitlines() if not line.startswith('>')))
def write_json(path,value):path.write_text(json.dumps(value,sort_keys=True,indent=2)+'\n')
def write_lines(path,rows):
 count=0
 with pathlib.Path(path).open('xb') as stream:
  with gzip.GzipFile(filename='',mode='wb',fileobj=stream,mtime=0,compresslevel=6) as z:
   for row in rows:
    z.write((json.dumps(row,sort_keys=True,separators=(',',':'))+'\n').encode());count+=1
 return count

def build(output,reference_path,author_cache):
 output=pathlib.Path(output)
 require(not output.exists(),'output exists; choose a new directory')
 require(filehash(reference_path)==REFERENCE_SHA,'locked hg38 reference SHA mismatch')
 inputs={}
 def load(rel):
  path=ROOT/rel;inputs[rel]={'bytes':path.stat().st_size,'sha256':filehash(path)};return json.loads(path.read_text())
 def bind(rel,expected=None):
  path=ROOT/rel;digest=filehash(path)
  if expected is not None:require(digest==expected,'old source SHA mismatch: '+rel)
  inputs[rel]={'bytes':path.stat().st_size,'sha256':digest}
  return path
 # Bind every retained old raw source, including raw source/header/QC manifests.
 panel_manifest=load('datasets/context-panel/manifest.json')
 for name,expected in panel_manifest['files'].items():bind('datasets/context-panel/'+name,expected)
 panel=json.loads(bind('datasets/context-panel/panel.json',PANEL_SHA).read_text())
 grouping=load('research/results/context-activity-v3/core/grouping.json')
 activity_manifest=load('research/results/context-activity-v3/core/artifacts.json')
 require(grouping['panel_sha256']==PANEL_SHA,'old grouping DNA binding mismatch')
 training_folds=collections.defaultdict(list)
 for fold in grouping['folds']:
  for identifier in fold['train_ids']:training_folds[identifier].append(fold['fold'])
 cohort_manifest=load('datasets/regulatory-cohort/manifest.json')
 for name,expected in cohort_manifest['files'].items():bind('datasets/regulatory-cohort/'+name,expected)
 selection=load('datasets/regulatory-cohort/metadata-selection.json')
 require(selection['assays']==[{k:v for k,v in a.items() if k!='scores_file'} for a in cohort_manifest['assays']],'cohort metadata/source manifest disagreement')
 acquisitions=load('research/variant-transfer/prior-exposure-inventory/author-acquisition.json')
 author_receipts={q['name']:q for q in acquisitions['files']}
 metadata=load('research/forensic-sources/osf-grch38-files.json')
 osf={q['attributes']['name']:q['attributes'] for q in metadata['data']}
 reference_receipts=[];agarwal=[];loci={};edits={};row_bindings=[];author_unmapped=[]
 with TwoBitReader(reference_path) as reference:
  for record in panel['records']:
   s=dna(record['sequence']);strand=record['source_metadata']['reference_match_strand']
   require(len(s)==200,'old panel core length differs')
   genome=reference.fetch(record['chromosome'],record['start'],record['end'],strand=strand)
   require(genome==s,'old panel no longer equals locked reference: '+record['id'])
   agarwal.append({'id':record['id'],'collection':'Agarwal2025','assembly':'GRCh38','chromosome':record['chromosome'],'start':record['start'],'end':record['end'],'coordinate_convention':'zero_based_half_open','reference_strand':strand,'prior_roles':['explored_specificity','explored_published_activity'],'old_training_folds':training_folds[record['id']],**seqmeta(s)})
  require(len(agarwal)==55276 and len({q['id'] for q in agarwal})==55276,'old panel coverage changed')
  require(set(training_folds)=={q['id'] for q in agarwal},'five-fold training union is not the entire old panel')
  for assay in cohort_manifest['assays']:
   element=assay['element_id'];mapping=assay['mapping'];s=dna(assay['target_sequence']);short=assay['assay_id'].split(':')[-1]
   require(mapping['assembly']=='GRCh38','old promoter is not declared GRCh38')
   require(mapping['end']-mapping['start']==len(s),'old target length/coordinate mismatch')
   genome=reference.fetch('chr'+mapping['chromosome'],mapping['start'],mapping['end'],strand=mapping['target_reference_alignment_strand'])
   require(genome==s,'old promoter fails exact full-reference validation: '+element)
   source=load('datasets/regulatory-cohort/'+assay['metadata_file'])
   require(source['extraMetadata']['reference']=='hg38','old promoter source build declaration changed')
   require((int(source['extraMetadata']['start'])-1,int(source['extraMetadata']['end']))==(mapping['start'],mapping['end']),'old source one-based inclusive conversion mismatch')
   require(source['targetGenes'][0]['targetSequence']['sequence']==s,'source target DNA mismatch')
   if element not in loci:
    loci[element]={'id':'Kircher2019:'+element,'collection':'Kircher2019','element_id':element,'assembly':'GRCh38','chromosome':'chr'+mapping['chromosome'],'start':mapping['start'],'end':mapping['end'],'coordinate_convention':'zero_based_half_open','reference_strand':mapping['target_reference_alignment_strand'],'source_coordinate_convention':'one_based_inclusive','genome_wide_uniqueness_established':False,'assays':[],**seqmeta(s)}
   locus=loci[element];require(locus['sequence']==s,'same old locus has inconsistent DNA')
   locus['assays'].append({'assay_id':assay['assay_id'],'prior_role':assay['role'],'cell_context':assay['cell_context'],'metadata_path':'datasets/regulatory-cohort/'+assay['metadata_file'],'scores_path':'datasets/regulatory-cohort/'+assay['scores_file']})
   rawname=AUTHOR[short];rawpath=pathlib.Path(author_cache)/rawname;expected=author_receipts[rawname]
   require(filehash(rawpath)==expected['sha256']==osf[rawname]['extra']['hashes']['sha256'],'old author table SHA mismatch')
   require(rawpath.stat().st_size==expected['bytes'],'old author table size mismatch')
   original={};omitted=[]
   for line,row in enumerate(csv.DictReader(rawpath.open(),delimiter='\t'),2):
    if row['Position']=='NA' and all(row[k]=='NA' for k in ['Chromosome','Ref','Alt']):
     omitted.append({'source_line_1based':line,'reason':'All genomic identity fields NA in original author file; not present in previously explored MaveDB table.'});continue
    pos=int(row['Position'])-1;offset=pos-mapping['start']
    require(row['Chromosome']==mapping['chromosome'] and 0<=offset<len(s),'old author coordinate outside exact promoter')
    require(row['Ref']==s[offset],'old author reference allele fails exact DNA')
    alt=row['Alt'];require(alt=='-' or (alt in 'ACGT' and len(alt)==1 and alt!=row['Ref']),'unsupported old author allele')
    key=offset,alt;require(key not in original,'duplicate old author molecular identity key')
    original[key]={'source_line_1based':line,'reference':row['Ref'],'alternate':alt,'genomic_position_1based':pos+1}
   with bind('datasets/regulatory-cohort/'+assay['scores_file']).open() as f:
    score_rows=list(csv.DictReader(f))
   require(len(score_rows)==len(original)==assay['expected_source_rows'],'old source row identity coverage mismatch')
   seen=set();counts=collections.Counter()
   for row in score_rows:
    match=PATTERN.fullmatch(row['hgvs_nt']);require(match is not None,'unsupported old HGVS identity')
    offset=int(match[1])-1;require(0<=offset<len(s),'old HGVS coordinate outside target')
    alt=match[3] or '-';key=offset,alt
    require(key in original and key not in seen,'old author/MaveDB identity not bijective');seen.add(key)
    if match[2]:require(match[2]==s[offset],'old MaveDB reference allele mismatch')
    typ='single_nucleotide_deletion' if alt=='-' else 'single_nucleotide_substitution'
    counts[typ]+=1;edited=s[:offset]+('' if alt=='-' else alt)+s[offset+1:]
    edit_key=(element,offset,alt)
    if edit_key not in edits:
     edits[edit_key]={'id':f'Kircher2019:{element}:{offset}:{alt}','locus_id':locus['id'],'element_id':element,'assembly':'GRCh38','chromosome':locus['chromosome'],'parent_start':mapping['start'],'parent_end':mapping['end'],'parent_reference_strand':mapping['target_reference_alignment_strand'],'offset_0based':offset,'reference':s[offset],'alternate':alt,'edit_type':typ,'genomic_position_0based':mapping['start']+offset,'coordinate_note':'Genomic position refers to the parent interval; edited sequence length differs for deletions. It is not a contiguous reference interval.','sources':[],**seqmeta(edited)}
    edits[edit_key]['sources'].append({'assay_id':assay['assay_id'],'mavedb_accession':row['accession'],'mavedb_literal_hgvs':row['hgvs_nt'],'author_file':rawname,'author_sha256':expected['sha256'],'author_source_line_1based':original[key]['source_line_1based']})
   row_bindings.append({'assay_id':assay['assay_id'],'element_id':element,'source_rows':len(score_rows),'author_mapped_rows':len(original),'molecular_identity_bijection':True,'representation_counts':dict(counts),'author_name':rawname,'author_sha256':expected['sha256'],'numeric_outcomes_compared':False,'note':'GP1BA author filename maps to source chr22 GP1BB by exact target and alleles; no gene-name-only mapping.' if element=='GP1BB' else None})
   if omitted:author_unmapped.append({'assay_id':assay['assay_id'],'author_file':rawname,'rows':omitted})
  aux=[]
  for path,fastaname,start_key,end_key,role in [('datasets/hbb-reference/manifest.json','genomic.fasta','sequence_start','sequence_end','explored_HBB_gene_reference_and_mechanism_work'),('datasets/clinvar-hbb/manifest.json','reference.fasta','reference_start','reference_end','explored_ClinVar_HBB_classification_reference')]:
   a=load(path);base=str(pathlib.PurePosixPath(path).parent);fasta=bind(base+'/'+fastaname,a['files'][fastaname]);s=read_fasta(fasta)
   require(a['assembly']=='GRCh38' and len(s)==a[end_key]-a[start_key],'HBB auxiliary reference interval mismatch')
   require(reference.fetch('chr11',a[start_key],a[end_key])==s,'HBB auxiliary reference fails locked hg38 comparison')
   aux.append({'id':base,'collection':'HBB_reference','assembly':'GRCh38','chromosome':'chr11','start':a[start_key],'end':a[end_key],'coordinate_convention':'zero_based_half_open','reference_strand':'+','prior_roles':[role],**seqmeta(s)})
 # No numeric source outcome fields are emitted. Every old identity is kept, without QC.
 output.mkdir(parents=True)
 write_lines(output/'agarwal-cores.jsonl.gz',sorted(agarwal,key=lambda q:q['id']))
 write_json(output/'kircher-loci.json',{'coordinate_convention':'zero_based_half_open','loci':sorted(loci.values(),key=lambda q:q['id'])})
 write_json(output/'hbb-auxiliary-loci.json',{'loci':aux})
 native_edits=sorted(edits.values(),key=lambda q:q['id'])
 write_lines(output/'kircher-native-edits.jsonl.gz',native_edits)
 reference_windows=[];ref_window_sequences=set();snv_by_element=collections.defaultdict(lambda:collections.defaultdict(set))
 for e in native_edits:
  if e['edit_type']=='single_nucleotide_substitution':snv_by_element[e['element_id']][e['offset_0based']].add(e['alternate'])
 for locus in sorted(loci.values(),key=lambda q:q['id']):
  require(locus['reference_strand']=='+','window coordinate conversion currently requires validated plus-strand targets')
  for offset in range(max(0,len(locus['sequence'])-199)):
   s=locus['sequence'][offset:offset+200];ref_window_sequences.add(s);ref_window_sequences.add(rc(s))
   mutations=[[pos-offset,''.join(sorted(alts))] for pos,alts in sorted(snv_by_element[locus['element_id']].items()) if offset<=pos<offset+200]
   reference_windows.append({'id':locus['id']+f':window200:{offset}','locus_id':locus['id'],'element_id':locus['element_id'],'assembly':'GRCh38','chromosome':locus['chromosome'],'start':locus['start']+offset,'end':locus['start']+offset+200,'reference_strand':'+','parent_offset_0based':offset,'derived_window_not_native_construct':True,'observed_snv_alts':mutations,**seqmeta(s)})
 write_lines(output/'kircher-reference-windows200.jsonl.gz',reference_windows)
 # Only deletion-created windows not already present in any native reference are
 # needed beyond the reference+single-substitution factorization. Deduplicate DNA,
 # retaining every parent/edit/window occurrence for an audit trail.
 deletion_windows={};all_deletion_occurrences=0
 for edit in native_edits:
  if edit['edit_type']!='single_nucleotide_deletion':continue
  s=edit['sequence']
  for offset in range(max(0,len(s)-199)):
   window=s[offset:offset+200]
   if window in ref_window_sequences:continue
   all_deletion_occurrences+=1;canonical=min(window,rc(window));key=sha(canonical.encode())
   if key not in deletion_windows:deletion_windows[key]={'id':'Kircher2019:deletion_window200:'+key,'derived_window_not_native_construct':True,'coordinate_note':'Derived edited-DNA substring; use parent_locus_id for genomic guard. No contiguous reference mapping is asserted.','occurrences':[],**seqmeta(canonical)}
   deletion_windows[key]['occurrences'].append({'native_edit_id':edit['id'],'parent_locus_id':edit['locus_id'],'edited_sequence_offset_0based':offset,'canonical_strand':'+' if window==canonical else '-'})
 write_lines(output/'kircher-deletion-windows200.jsonl.gz',sorted(deletion_windows.values(),key=lambda q:q['id']))
 # Auxiliary HBB windows represent explored reference material, not reporter assays.
 aux_windows=[]
 for locus in aux:
  for offset in range(max(0,len(locus['sequence'])-199)):
   s=locus['sequence'][offset:offset+200]
   aux_windows.append({'id':locus['id']+f':window200:{offset}','parent_locus_id':locus['id'],'assembly':'GRCh38','chromosome':'chr11','start':locus['start']+offset,'end':locus['start']+offset+200,'reference_strand':'+','derived_window_not_native_construct':True,**seqmeta(s)})
 write_lines(output/'hbb-auxiliary-windows200.jsonl.gz',aux_windows)
 write_json(output/'molecular-identity-bindings.json',{'bindings':row_bindings,'unmapped_author_rows_not_in_explored_mavedb':author_unmapped,'numeric_old_outcome_fields_read':False,'note':'CSV parsers read whole source rows as text, but only accession/HGVS or chromosome/position/reference/alternate fields are inspected; scores/counts are neither converted nor emitted. No Siraj outcome source is opened.'})
 write_json(output/'inputs.json',{'reference':{'path':str(reference_path),'sha256':REFERENCE_SHA,'bytes':pathlib.Path(reference_path).stat().st_size,'assembly':'GRCh38'},'files':inputs,'old_author_tables':acquisitions['files']})
 summary={'schema_version':1,'complete_for_declared_prior_exposure_scope':True,'scope':'All 55,276 old Agarwal cores, all 15 explored Kircher assays at ten unique promoter loci (nine primary non-HBB plus HBB), every source-resolved old SNV/deletion identity irrespective of QC, and explored HBB genomic/ClinVar reference intervals. This is a DNA/locus exposure inventory, not a new biological cohort.','assembly':'GRCh38','coordinate_convention':'zero_based_half_open','reference_sha256':REFERENCE_SHA,'agarwal_cores':len(agarwal),'agarwal_five_fold_training_union':len(training_folds),'kircher_assays':len(row_bindings),'kircher_unique_loci':len(loci),'kircher_explored_mavedb_rows':sum(q['source_rows'] for q in row_bindings),'kircher_snv_source_rows':sum(q['representation_counts']['single_nucleotide_substitution'] for q in row_bindings),'kircher_deletion_source_rows':sum(q['representation_counts']['single_nucleotide_deletion'] for q in row_bindings),'kircher_distinct_native_edits':len(native_edits),'kircher_distinct_snv_edits':sum(e['edit_type']=='single_nucleotide_substitution' for e in native_edits),'kircher_distinct_deletion_edits':sum(e['edit_type']=='single_nucleotide_deletion' for e in native_edits),'kircher_reference_windows200':len(reference_windows),'kircher_unique_new_deletion_windows200':len(deletion_windows),'kircher_new_deletion_window_occurrences':all_deletion_occurrences,'hbb_auxiliary_loci':len(aux),'hbb_auxiliary_windows200':len(aux_windows),'all_full_reference_checks_passed':True,'all_old_molecular_identity_bijections_passed':True,'new_study_outcome_access':False,'window_semantics':'Full native constructs are retained at original lengths. The 200nt windows are derived exclusion-only DNA substrings, not claimed native reporter constructs. HBB187nt and its186nt deletion constructs have no200nt windows; their locus is still excluded by the genomic guard.','snv_homology_factorization':'For a query200nt sequence, match each reference window on both strands. <=20 mismatches qualifies directly. At21 mismatches, an exposed single substitution can qualify iff its documented position/alternate matches the query at one mismatch. >21 mismatches cannot be rescued by one substitution. This exactly represents >=90% homology to every exposed single-SNV200nt window without materializing about one million copies. All distinct deletion-created windows are explicit.','limits':['No new source target or response eligibility is chosen here; apply the separate frozen protocol to this inventory.','Exact reference checks establish interval/strand consistency, not genome-wide uniqueness or biological causality.','DNA-only author/MaveDB correspondence does not assert numerical outcome equality.','Seven all-NA identity rows in the original MSMB file were not in the old explored MaveDB table and are explicitly recorded rather than assigned DNA.','GRCh38_GP1BA.tsv is the author filename for the exact chr22 GP1BB source target here; the mapping is established from full DNA and every allele, not the filename.']}
 write_json(output/'summary.json',summary)
 shutil.copyfile(pathlib.Path(__file__),output/'implementation.py')
 shutil.copyfile(ROOT/'scripts/twobit_reference.py',output/'twobit_reference.py')
 files={p.name:{'bytes':p.stat().st_size,'sha256':filehash(p)} for p in sorted(output.iterdir()) if p.is_file()}
 write_json(output/'manifest.json',{'schema_version':1,'files':files,'summary_sha256':files['summary.json']['sha256']})
 print(json.dumps(summary,sort_keys=True))

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=pathlib.Path,required=True);p.add_argument('--reference',type=pathlib.Path,default=pathlib.Path('/workspace/.tools/reference/hg38.2bit'));p.add_argument('--author-cache',type=pathlib.Path,default=pathlib.Path('/workspace/.tools/prior-exposure-author-tables'));a=p.parse_args();build(a.out,a.reference,a.author_cache)
if __name__=='__main__':main()
