import collections,gzip,json,pathlib,hashlib,re,openpyxl,datetime
p=pathlib.Path('/tmp/helix-variant-source')
f=p/'siraj2026-S2-oligo-sequences.xlsx'
old=json.load(open('/workspace/helix/datasets/context-panel/panel.json'))['records']
rc=lambda s:s.translate(str.maketrans('ACGT','TGCA'))[::-1]
oldseq=collections.defaultdict(list)
for r in old:
 oldseq[min(r['sequence'],rc(r['sequence']))].append(r['id'])
wb=openpyxl.load_workbook(f,read_only=True,data_only=False);sheet=wb['S2.Variant Oligo Seqs']
assert [c.value for c in next(sheet.iter_rows(min_row=2,max_row=2))]==['Variant','Allele 1 Oligo','Allele 2 Oligo']
count=collections.Counter();lens=collections.Counter();diffpos=collections.Counter();seen=collections.Counter();issues=[];overlap=[];bad_examples=[]
out=p/'siraj2026-S2-DNA.tsv.gz'
with gzip.GzipFile(filename=str(out),mode='wb',mtime=0) as g:
 g.write(b'row_1based\tvariant\tallele_1_oligo\tallele_2_oligo\n')
 for nr,row in enumerate(sheet.iter_rows(min_row=3),3):
  vals=[c.value for c in row];count['stored_rows']+=1
  if len(vals)!=3 or not all(isinstance(v,str) for v in vals):
   count['non_text_or_missing']+=1;issues.append({'row':nr,'types':[c.data_type for c in row]});continue
  vid,a,b=vals;seen[vid]+=1;g.write(('\t'.join([str(nr),vid,a,b])+'\n').encode());lens[f'{len(a)},{len(b)}']+=1
  if set(a)|set(b)<=set('ACGT'):count['ACGT_only_pair']+=1
  else:count['non_ACGT_pair']+=1
  parts=vid.split(':');is_snv=len(parts)==4 and re.fullmatch(r'chr(?:[1-9]|1\d|2[0-2]|X|Y|M)',parts[0]) and parts[1].isdigit() and len(parts[2])==len(parts[3])==1 and parts[2] in 'ACGT' and parts[3] in 'ACGT' and parts[2]!=parts[3]
  if is_snv:
   count['nominal_snv']+=1
   if len(a)==len(b)==200:
    d=[i for i,(x,y) in enumerate(zip(a,b)) if x!=y]
    if len(d)==1:
     count['nominal_snv_200nt_single_difference']+=1;diffpos[d[0]]+=1
     if a[d[0]]==parts[2] and b[d[0]]==parts[3]:count['single_difference_matches_ID_ref_alt']+=1
     elif len(bad_examples)<15:bad_examples.append({'row':nr,'id':vid,'diffpos':d,'base1':a[d[0]],'base2':b[d[0]]})
    else:count['nominal_snv_200nt_not_single_difference']+=1
   else:count['nominal_snv_not_200nt_pair']+=1
  else:count['not_nominal_primary_chromosome_snv']+=1
  for allele,s in [(1,a),(2,b)]:
   if set(s)<=set('ACGT'):
    match=oldseq.get(min(s,rc(s)))
    if match:overlap.append({'source_row':nr,'variant':vid,'allele':allele,'old_panel_ids':match})
wb.close()
summary={'timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'old_panel_sha256':hashlib.sha256(pathlib.Path('/workspace/helix/datasets/context-panel/panel.json').read_bytes()).hexdigest(),'purpose':'DNA-only pre-outcome source inventory, no target values acquired or used; counts do not select final population','counts':dict(count),'unique_ids':len(seen),'duplicate_ids':{k:v for k,v in seen.items() if v>1},'length_pairs':dict(lens),'SNV_single_difference_zero_based_positions':dict(diffpos),'mismatch_examples':bad_examples,'nontext_issues':issues,'exact_or_RC_overlap_rows':overlap,'DNA_export_sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'source_build_status':'likely hg19 based on Methods liftover and code; must validate before coordinate comparison'}
(p/'siraj2026-S2-DNA-inventory.json').write_text(json.dumps(summary,indent=2));print(json.dumps({k:v for k,v in summary.items() if k not in ['exact_or_RC_overlap_rows','duplicate_ids','nontext_issues']},indent=2));print('overlaprecords',len(overlap),'duplicateIDs',len(summary['duplicate_ids']))
