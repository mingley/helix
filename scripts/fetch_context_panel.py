#!/usr/bin/env python3
"""Strict Agarwal source import with separate failed-activity and ESS endpoints.

No network is used here: an independently checked retrieval receipt names the
untouched workbooks. The caller must bind an announced frozen protocol hash.
"""
from __future__ import annotations
import argparse
import collections
import hashlib
import gzip
import json
import math
from pathlib import Path
import posixpath
import re
import tempfile
import xml.etree.ElementTree as ET
import zipfile

from twobit_reference import TwoBitReader

NS = {'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
REL = 'http://schemas.openxmlformats.org/officeDocument/2006/relationships'
CTX = ['HepG2','K562','WTC11']
S7 = 'Table-S7-joint-design.xlsx'
S9 = 'Table-S9-folds-performance.xlsx'
SEQ_HEADER = "230nt sequence (15nt 5' adaptor - 200nt element - 15nt 3' adaptor)"
DESIGN_HEADERS = ['name','category','chr.hg19','start.hg19','stop.hg19','str.hg19','chr.hg38','start.hg38','stop.hg38','str.hg38',SEQ_HEADER]
ESS_HEADERS = ['Fold','Sequence ID'] + [f'{cell}, Element Specificity Score' for cell in CTX]
REFERENCE = {'size_bytes':835393456,'md5':'dcc3ea27079aa6dc3f9deccd7275e0f8','sha256':'1f67aaa17a77b327738fe750ab37430a85dddf73c7d9189385ad1839259acec0'}


def digest(raw): return hashlib.sha256(raw).hexdigest()


def hashes(path):
    sha,md5=hashlib.sha256(),hashlib.md5()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):
            sha.update(block);md5.update(block)
    return {'sha256':sha.hexdigest(),'md5':md5.hexdigest(),'size_bytes':Path(path).stat().st_size}


def write_new(path,obj):
    with Path(path).open('x') as stream:
        json.dump(obj,stream,indent=2,allow_nan=False);stream.write('\n')


def verified_inputs(directory,protocol,expected):
    root=Path(directory).resolve();raw_root=root/'raw'
    protocol_raw=Path(protocol).read_bytes()
    if digest(protocol_raw)!=expected:raise ValueError('frozen protocol SHA256 changed')
    specification=json.loads(protocol_raw)
    if specification.get('status')!='frozen':raise ValueError('protocol is not frozen')
    receipt_raw=(raw_root/'retrieval.json').read_bytes();receipt=json.loads(receipt_raw)
    license_path=Path('research/context-source-alternatives/zenodo10558183-metadata.json')
    license_bytes=license_path.read_bytes();license_metadata=json.loads(license_bytes)
    if str(license_metadata['id'])!='10558183' or license_metadata['metadata']['license']['id']!='cc-by-4.0':
        raise ValueError('source Zenodo record does not declare the expected data grant')
    retained_license=raw_root/'zenodo-record-metadata.json'
    if retained_license.exists():
        if retained_license.read_bytes()!=license_bytes:raise ValueError('retained source-license metadata differs')
    else:
        with retained_license.open('xb') as stream:stream.write(license_bytes)
    # V2 is deliberately a new interpretation of raw workbooks acquired under V1.
    receipts=[]
    allowed={entry['key']:entry for entry in [specification['source']['design_workbook'],specification['source']['folds_workbook']]}
    for item in receipt['files']:
        if item['original_filename'] not in allowed:raise ValueError('receipt names an unpermitted source')
        path=(raw_root/item['file']).resolve()
        if path.parent!=raw_root.resolve():raise ValueError('source path escapes raw directory')
        actual=hashes(path);declared=allowed[item['original_filename']]
        if actual['size_bytes']!=declared['size_bytes'] or actual['md5']!=declared['md5']:
            raise ValueError('source size/published MD5 mismatch')
        if actual['sha256']!=item['sha256'] or actual['md5']!=item['observed_md5'] or not item['md5_verified']:
            raise ValueError('source retrieval receipt does not match bytes')
        receipts.append(dict(item,independently_rechecked=actual))
    if {item['file'] for item in receipts}!={S7,S9}:raise ValueError('both bounded workbooks are required')
    return root,specification,dict(protocol_sha256=expected,raw_acquisition_protocol_sha256=receipt['protocol_sha256'],retrieval_receipt_sha256=digest(receipt_raw),source_license_metadata_sha256=digest(license_bytes),source_license='CC-BY-4.0 declared by the source Zenodo record, covering its listed workbooks; retain source attribution',sources=receipts)


class Workbook:
    """Small read-only OpenXML parser; ZIP CRC checks remain enabled."""
    def __init__(self,path):
        self.zip=zipfile.ZipFile(path)
        if sum(i.file_size for i in self.zip.infolist())>300_000_000:raise ValueError('workbook uncompressed size exceeds budget')
        for info in self.zip.infolist():
            if info.flag_bits & 1:raise ValueError('encrypted workbook unsupported')
        self.shared=[]
        if 'xl/sharedStrings.xml' in self.zip.namelist():
            tree=ET.fromstring(self.zip.read('xl/sharedStrings.xml'))
            self.shared=[''.join(n.itertext()) for n in tree.findall('m:si',NS)]
        book=ET.fromstring(self.zip.read('xl/workbook.xml'))
        relationships=ET.fromstring(self.zip.read('xl/_rels/workbook.xml.rels'))
        targets={node.attrib['Id']:node.attrib['Target'] for node in relationships if not node.attrib.get('TargetMode')=='External'}
        self.sheets={}
        for sheet in book.findall('m:sheets/m:sheet',NS):
            target=targets[sheet.attrib['{'+REL+'}id']]
            path=posixpath.normpath(target.lstrip('/') if target.startswith('/') else 'xl/'+target)
            if not path.startswith('xl/') or '..' in path.split('/'):raise ValueError('workbook relationship escapes archive')
            self.sheets[sheet.attrib['name']]=path
    def rows(self,sheet,columns=None,stop_after=None):
        if sheet not in self.sheets:raise ValueError('expected source sheet is absent')
        with self.zip.open(self.sheets[sheet]) as stream:
            for event,node in ET.iterparse(stream,events=('end',)):
                if node.tag!='{'+NS['m']+'}row':continue
                row_number=int(node.attrib['r']);cells={}
                for cell in node.findall('m:c',NS):
                    column=re.match(r'[A-Z]+',cell.attrib['r']).group()
                    if columns is not None and column not in columns:continue
                    kind=cell.attrib.get('t','n');value=cell.findtext('m:v',None,NS)
                    if kind=='s':value=self.shared[int(value)] if value is not None else None
                    elif kind=='inlineStr':value=''.join(cell.find('m:is',NS).itertext())
                    elif kind=='e':value={'excel_error':value}
                    cells[column]=value
                yield row_number,cells
                node.clear()
                if stop_after is not None and row_number>=stop_after:return
    def close(self):self.zip.close()
    def __enter__(self):return self
    def __exit__(self,*_):self.close()


def row_values(cells,count):return [cells.get(chr(65+i)) for i in range(count)]

def header_check(path,sheet,row,expected):
    with Workbook(path) as workbook:
        found=next((cells for n,cells in workbook.rows(sheet,stop_after=row) if n==row),None)
    if found is None or row_values(found,len(expected))!=expected:raise ValueError('source header differs from frozen schema binding')


def source_schema(root):
    raw=root/'raw'
    header_check(raw/S7,'joint library',2,DESIGN_HEADERS)
    header_check(raw/S9,'Fold assignment, joint library',1,ESS_HEADERS)
    return {'schema_version':1,'workbooks':{
        S7:{'role':'design_and_coordinate_metadata_only','sheet':'joint library','header_row':2,'headers':DESIGN_HEADERS},
        S9:{'role':'published_element_specificity_score_not_raw_reporter_activity','sheet':'Fold assignment, joint library','header_row':1,'headers':ESS_HEADERS}},
        'v1_required_endpoint':'mean-of-three-median-normalized-log2-rna-dna-activity',
        'v1_semantic_endpoint_acceptance':False,
        'v1_failure_reasons':['S7 contains design and coordinate metadata without observed activity columns.','S9 labels each cell column Element Specificity Score, explicitly forbidden as a V1 observed activity label.','CNN prediction output cannot fill missing measurements.'],
        'no_ess_values_parsed_for_v1':True,
        'core_extraction':{'explicit_source_header':SEQ_HEADER,'assayed_oligo_length':230,'left_adapter_bases':15,'genomic_core_bases':200,'right_adapter_bases':15,'python_slice':'sequence[15:215]'},
        'coordinate_binding':'Declared chr.hg38/start.hg38/stop.hg38; verify each 200nt core against independent hg38 with the unmodified declared interval. No coordinate searching or silent liftover.',
        'parent_binding':'No separate parent field in S7; preserve stable design IDs and let ROOT-owned grouping establish documented unions.'}


def integer(value):
    if value is None or isinstance(value,dict):raise ValueError('missing_integer')
    number=float(value)
    if not math.isfinite(number) or number<0 or number!=int(number):raise ValueError('invalid_integer')
    return int(number)


def chromosome(value):
    if value is None or isinstance(value,dict) or value=='NA':raise ValueError('missing_chromosome')
    if re.fullmatch(r'\d+(?:\.0+)?',value):name=str(integer(value))
    else:name=value.removeprefix('chr')
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*',name):raise ValueError('invalid_declared_chromosome')
    return 'chr'+name


def validate_design(root,reference):
    observed=hashes(reference)
    if observed!=REFERENCE:raise ValueError('independent hg38 reference size/hash mismatch')
    records=[];excluded=[];reasons=collections.Counter();strands=collections.Counter();ids=set();all_rows=0;declared_disagreements=0
    with Workbook(root/'raw'/S7) as workbook,TwoBitReader(reference) as genome:
        for number,cells in workbook.rows('joint library'):
            if number<=2:continue
            all_rows+=1;row=row_values(cells,11);identifier=row[0]
            metadata={'source_row_1based':number,'id':identifier,'category':row[1],'source_hg38_fields':row[6:10]}
            try:
                if not isinstance(identifier,str) or not identifier or '\n' in identifier:raise ValueError('invalid_source_id')
                if identifier in ids:raise ValueError('duplicate_source_id')
                ids.add(identifier)
                full=row[10]
                if not isinstance(full,str) or len(full)!=230 or set(full)-set('ACGT'):raise ValueError('invalid_assayed_oligo')
                chrom=chromosome(row[6]);start=integer(row[7]);end=integer(row[8])
                if end-start!=200:raise ValueError('genomic_core_coordinate_length_mismatch')
                core=full[15:215]
                forward=genome.fetch(chrom,start,end)
                reverse=forward.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
                matches=[strand for strand,seq in [('+',forward),('-',reverse)] if seq==core]
                if len(matches)==0:raise ValueError('deposited_core_does_not_match_declared_hg38_interval')
                if len(matches)>1:raise ValueError('ambiguous_palindromic_reference_orientation')
                strand=matches[0]
                if row[9] not in ('+','-'):raise ValueError('missing_or_invalid_declared_strand')
                if row[9]!=strand:
                    declared_disagreements+=1
                    raise ValueError('declared_strand_disagrees_with_exact_reference')
                strands[strand]+=1
                records.append({'id':identifier,'sequence':core,'chromosome':chrom,'start':start,'end':end,'parent_id':identifier,
                    'source_metadata':dict(metadata,assayed_oligo=full,left_adapter=full[:15],right_adapter=full[215:],assayed_oligo_sha256=digest(full.encode()),core_sha256=digest(core.encode()),source_declared_strand=row[9],reference_match_strand=strand,declared_strand_matches=(row[9]==strand),parent_policy='No separate documented parent column; stable ID used as self-parent without inventing cross-ID relationships.')})
            except (ValueError,KeyError) as error:
                reason=str(error);reasons[reason]+=1;excluded.append(dict(metadata,reason=reason))
    return records,excluded,{'schema_version':1,'source_design_rows':all_rows,'reference_checked_design_records':len(records),'exact_forward_matches':strands['+'],'reverse_complement_matches':strands['-'],'excluded_counts':dict(sorted(reasons.items())),'declared_strand_disagreements':declared_disagreements,'reference_sha256':observed['sha256'],'reference_md5':observed['md5'],'reference_size_bytes':observed['size_bytes'],'reference_verified':True,'source_coordinate_convention':'zero-based half-open supported by exact200nt matches at the unshifted declared intervals','endpoint_status':'design_only_no_observations'}


def design_only(args):
    root,spec,receipt=verified_inputs(args.directory,args.protocol,args.expected_protocol_sha256)
    schema=source_schema(root)
    if spec.get('protocol_id')!='context-prediction-v1':raise ValueError('design-only V1 receipt must retain original endpoint identity')
    records,excluded,reference=validate_design(root,args.reference)
    write_new(root/'v1-source-schema-binding.json',dict(schema,input_receipt=receipt))
    write_new(root/'v1-blocked-acquisition.json',{'schema_version':1,'protocol_id':spec['protocol_id'],'protocol_sha256':receipt['protocol_sha256'],'status':'blocked_missing_required_observed_activity_endpoint','fitting_performed':False,'source_schema_binding_sha256':digest((root/'v1-source-schema-binding.json').read_bytes()),'reasons':schema['v1_failure_reasons']})
    write_new(root/'design-catalog.json',{'schema_version':1,'endpoint':'design_only_no_observations','records':records})
    write_new(root/'design-exclusions.json',{'schema_version':1,'records':excluded})
    write_new(root/'design-reference-validation.json',reference)
    print(json.dumps({'v1_status':'blocked_missing_required_endpoint','design_rows':reference['source_design_rows'],'valid_reference_designs':len(records),'excluded_designs':len(excluded),'no_ess_values_parsed':True}))


def write_compact(path,obj):
    with Path(path).open('x') as stream:
        if isinstance(obj,dict) and isinstance(obj.get('records'),list):
            metadata={key:value for key,value in obj.items() if key!='records'}
            stream.write(json.dumps(metadata,separators=(',',':'),allow_nan=False)[:-1]+',"records":[\n')
            for index,record in enumerate(obj['records']):
                if index:stream.write(',\n')
                stream.write(json.dumps(record,separators=(',',':'),allow_nan=False))
            stream.write('\n]}\n')
        else:
            json.dump(obj,stream,separators=(',',':'),allow_nan=False);stream.write('\n')


def write_gzip_new(path,obj):
    payload=(json.dumps(obj,separators=(',',':'),allow_nan=False)+'\n').encode()
    with Path(path).open('xb') as stream:stream.write(gzip.compress(payload,mtime=0))


def source_number(value):
    if value is None:return None,'missing'
    if isinstance(value,dict):return None,'excel_error'
    if not isinstance(value,str):return None,'invalid_encoding'
    if value.strip().upper() in ('','NA','N/A','NULL','.'):return None,'missing'
    try: result=float(value)
    except ValueError:return None,'non_numeric'
    if not math.isfinite(result):return None,'non_finite'
    return result,'finite'


def build_ess(args):
    root,spec,receipt=verified_inputs(args.directory,args.protocol,args.expected_protocol_sha256)
    if spec.get('protocol_id')!='context-specificity-v2' or spec.get('frozen_before_ess_numeric_value_inspection_and_model_fitting') is not True:
        raise ValueError('ESS requires the separately frozen context-specificity-v2 protocol')
    if spec['semantic_endpoint_contract']['endpoint_id']!='published-three-context-element-specificity-score':
        raise ValueError('unexpected semantic target; never reinterpret ESS as raw activity')
    if (root/'panel.json').exists() or (root/'manifest.json').exists():raise ValueError('output exists; preserve the previous panel')
    if spec['schema_binding']['design']['all_header_columns']!=DESIGN_HEADERS or spec['schema_binding']['scores']['context_columns']!={c:ESS_HEADERS[2+i] for i,c in enumerate(CTX)}:
        raise ValueError('frozen schema headers differ from explicit source mapping')
    schema=source_schema(root)
    # The new gate is persisted before reading any numeric ESS cell.
    gate=dict(schema_version=1,protocol_id=spec['protocol_id'],protocol_sha256=receipt['protocol_sha256'],raw_acquisition_protocol_sha256=receipt['raw_acquisition_protocol_sha256'],raw_source_receipt_sha256=receipt['retrieval_receipt_sha256'],frozen_before_ess_numeric_inspection=True,ess_interpretation='Published source ESS unchanged; source whole-panel preprocessing is not training-only.')
    write_new(root/'v2-ess-inspection-gate.json',gate)
    designs,design_exclusions,reference_report=validate_design(root,args.reference)
    design_by_id={r['id']:r for r in designs}
    invalid_design_by_id={r['id']:r for r in design_exclusions if isinstance(r.get('id'),str)}
    known_design_ids=set(design_by_id)|set(invalid_design_by_id)
    included=[];excluded=[];context_missing=collections.Counter();source_rows=0;source_ids=set();blank_rows=[];score_statuses={c:collections.Counter() for c in CTX};metadata=[]
    with Workbook(root/'raw'/S9) as workbook:
        for row,cells in workbook.rows('Fold assignment, joint library'):
            if row==1:continue
            source_rows+=1;identifier=cells.get('B')
            if not isinstance(identifier,str) or not identifier:
                blank_rows.append(row);excluded.append({'source':'S9','source_row_1based':row,'id':identifier,'reason':'blank_source_id'});continue
            if identifier in source_ids:raise ValueError('duplicate source element observation')
            source_ids.add(identifier)
            observations={};raw_values={};bad={}
            for i,cell in enumerate(CTX):
                raw=cells.get(chr(67+i));raw_values[cell]=raw;value,status=source_number(raw);score_statuses[cell][status]+=1
                if status!='finite':bad[cell]=status;context_missing[cell]+=1
                else:observations[cell]=value
            if identifier not in known_design_ids:
                excluded.append({'source':'S9','source_row_1based':row,'id':identifier,'reason':'no_exact_source_design_id','source_ess_raw':raw_values});continue
            if identifier not in design_by_id:
                excluded.append({'source':'S7+S9','source_row_1based':row,'id':identifier,'reason':invalid_design_by_id[identifier]['reason'],'design':invalid_design_by_id[identifier],'source_ess_raw':raw_values});continue
            if bad:
                excluded.append({'source':'S9','source_row_1based':row,'id':identifier,'reason':'incomplete_or_nonfinite_context','contexts':bad,'source_ess_raw':raw_values});continue
            design=design_by_id[identifier];full_metadata=design['source_metadata']
            full=full_metadata['assayed_oligo'];core=design['sequence'];left=full_metadata['left_adapter'];right=full_metadata['right_adapter']
            metadata.append(dict(full_metadata,id=identifier,source_s9_row_1based=row,source_ess_raw=raw_values,left_adapter_sha256=digest(left.encode()),right_adapter_sha256=digest(right.encode())))
            included.append({k:design[k] for k in ('id','sequence','chromosome','start','end','parent_id')}|{'observations':observations,'source_metadata':{'s7_row_1based':full_metadata['source_row_1based'],'s9_row_1based':row,'reference_match_strand':full_metadata['reference_match_strand'],'source_declared_strand':full_metadata['source_declared_strand'],'source_design_metadata_record':identifier,'core_sha256':digest(core.encode())}})
    for identifier in sorted(known_design_ids-source_ids):
        excluded.append({'source':'S7','id':identifier,'reason':'no_exact_source_ess_id','design_source_row_1based':design_by_id[identifier]['source_metadata']['source_row_1based'] if identifier in design_by_id else invalid_design_by_id[identifier]['source_row_1based']})
    if not included:raise ValueError('no eligible complete three-context reference-validated records')
    included.sort(key=lambda r:r['id'])
    ref_counts=collections.Counter(r['source_metadata']['reference_match_strand'] for r in included)
    exclusion_counts=collections.Counter(r['reason'] for r in excluded)
    final_reference={'schema_version':1,'records_checked':len(included),'exact_forward_matches':ref_counts['+'],'reverse_complement_matches':ref_counts['-'],'excluded_counts':dict(sorted(exclusion_counts.items())),'reference_sha256':REFERENCE['sha256'],'reference_md5':REFERENCE['md5'],'reference_verified':True,'reference_size_bytes':REFERENCE['size_bytes'],'source_design_validation':reference_report,'match_policy':'Every admitted core200 matches independenthg38 at the exact unshifted declared coordinate and strand; no lookup rescue or liftover.'}
    inclusion={'schema_version':1,'endpoint':'published_element_specificity_score','source_s7_design_rows':reference_report['source_design_rows'],'source_s9_physical_data_rows':source_rows,'source_s9_nonempty_ids':len(source_ids),'source_s9_dimension_not_record_count':'A1:E59108','source_s9_implicit_empty_rows':59107-source_rows,'explicit_blank_id_rows':blank_rows,'admitted_records':len(included),'admitted_context_observations':len(included)*3,'excluded_counts':dict(sorted(exclusion_counts.items())),'missing_context_counts':dict(context_missing),'numeric_status_counts':{c:dict(score_statuses[c]) for c in CTX},'all_design_exclusions':design_exclusions,'excluded_records':excluded,'no_score_based_exclusions':True,'no_observation_imputation_or_target_recomputation':True,'id_join':'Exact case-sensitive S9 SequenceID equalsS7 name; no suffix/prefix/fuzzy transformations.'}
    binding={'schema_version':1,'protocol_id':spec['protocol_id'],'protocol_sha256':receipt['protocol_sha256'],'endpoint':'published_element_specificity_score','semantic_endpoint_acceptance':True,'required_workbooks':schema['workbooks'],'schema_binding':spec['schema_binding'],'source_semantic_evidence':spec['semantic_endpoint_contract'],'source_roles':{'S7':'design_sequence_coordinates','S9':'published_ESS_not_rawactivity','source_CNN_predictions':'not_consumed_not_empirical'},'source_version_caveat':spec['source']['source_version_caveat'],'no_new_normalization':True,'previous_v1_status':'blocked_missing_required_observed_activity_endpoint'}
    panel={'schema_version':1,'endpoint':'published_element_specificity_score','contexts':CTX,'records':included}
    write_new(root/'source-schema-binding-v2.json',binding)
    write_compact(root/'panel.json',panel)
    write_gzip_new(root/'source-design-metadata-v2.json.gz',{'schema_version':1,'role':'Full230ntassayedDNA,explicit15ntadapters,200ntcorehashes,strands,sourceRowsandunchangedESSnumericrepresentations','parent_policy':'No separate documented sourceparentfield; stable sourceID isselfparent; ROOT-owned grouping mayestablishadditionalrelationships.','records':metadata})
    write_new(root/'inclusion-report.json',inclusion)
    write_new(root/'reference-validation.json',final_reference)
    protocol_raw=Path(args.protocol).read_bytes()
    with (root/'frozen-protocol-v2.json').open('xb') as stream:stream.write(protocol_raw)
    files={str(p.relative_to(root)):digest(p.read_bytes()) for p in sorted(root.rglob('*')) if p.is_file() and p.name!='manifest.json'}
    manifest={'schema_version':1,'protocol_id':spec['protocol_id'],'protocol_sha256':receipt['protocol_sha256'],'endpoint':'published_element_specificity_score','contexts':CTX,'record_count':len(included),'files':files,'raw_acquisition':receipt,'source_version':spec['source']['source_version_caveat'],'source_endpoint':spec['semantic_endpoint_contract'],'source_selection':spec['source'].get('selection_caveat'),'reference':final_reference,'grouping_status':'ROOT-owned separate frozen grouping; no grouping or model fitting done by importer.','previous_failed_protocol':'context-prediction-v1','limitations':['PublishedESS includes whole-source-panel normalization and cell centering; these are not raw activities or train-only target preprocessing.','Three knowncellcontexts in one selectedsourcejointpanel; no unseen-cell or genomewide generalization claim.','Referenceverification is not empiricalfunctionalvalidation.']}
    write_new(root/'manifest.json',manifest)
    print(json.dumps({'endpoint':panel['endpoint'],'admitted_records':len(included),'admitted_context_observations':len(included)*3,'excluded_counts':dict(sorted(exclusion_counts.items())),'manifest_sha256':digest((root/'manifest.json').read_bytes())}))


def verify_panel(args):
    root=Path(args.directory).resolve();manifest=json.loads((root/'manifest.json').read_text())
    if manifest['endpoint']!='published_element_specificity_score' or manifest['protocol_id']!='context-specificity-v2':raise ValueError('unexpected panel endpoint')
    for name,expected in manifest['files'].items():
        path=(root/name).resolve()
        if not path.is_relative_to(root) or not path.is_file() or digest(path.read_bytes())!=expected:raise ValueError('panel source/artifact commitment mismatch: '+name)
    specification=(root/'frozen-protocol-v2.json').read_bytes()
    if digest(specification)!=manifest['protocol_sha256']:raise ValueError('retained protocol hash mismatch')
    panel=json.loads((root/'panel.json').read_text())
    if panel['schema_version']!=1 or panel['endpoint']!=manifest['endpoint'] or panel['contexts']!=CTX or len(panel['records'])!=manifest['record_count']:raise ValueError('panel identity/count mismatch')
    by_id={record['id']:record for record in panel['records']}
    if len(by_id)!=len(panel['records']):raise ValueError('duplicate admitted element ID')
    for record in by_id.values():
        if len(record['sequence'])!=200 or set(record['sequence'])-set('ACGT') or record['end']-record['start']!=200 or record['start']<0:raise ValueError('invalid admitted genomic core')
        if set(record['observations'])!=set(CTX) or any(isinstance(v,bool) or not isinstance(v,(float,int)) or not math.isfinite(v) for v in record['observations'].values()):raise ValueError('incomplete/nonfinite admitted ESS vector')
    source_schema(root)
    seen_scores=set();seen_design=set()
    with Workbook(root/'raw'/S9) as workbook:
        for row,cells in workbook.rows('Fold assignment, joint library'):
            identifier=cells.get('B')
            if row==1 or identifier not in by_id:continue
            if identifier in seen_scores:raise ValueError('duplicate source ESS observation')
            seen_scores.add(identifier)
            for i,cell in enumerate(CTX):
                value,status=source_number(cells.get(chr(67+i)))
                if status!='finite' or value!=by_id[identifier]['observations'][cell]:raise ValueError('admitted ESS differs from exact source cell')
    with Workbook(root/'raw'/S7) as workbook:
        for row,cells in workbook.rows('joint library'):
            identifier=cells.get('A')
            if row<=2 or identifier not in by_id:continue
            if identifier in seen_design:raise ValueError('duplicate source design')
            seen_design.add(identifier);record=by_id[identifier];full=cells.get('K')
            if not isinstance(full,str) or len(full)!=230 or full[15:215]!=record['sequence'] or chromosome(cells.get('G'))!=record['chromosome'] or integer(cells.get('H'))!=record['start'] or integer(cells.get('I'))!=record['end'] or cells.get('J')!=record['source_metadata']['source_declared_strand']:raise ValueError('admitted design differs from exact source fields')
    if seen_scores!=set(by_id) or seen_design!=set(by_id):raise ValueError('admitted source-ID coverage incomplete')
    reference_checked=False
    if args.reference:
        if hashes(args.reference)!=REFERENCE:raise ValueError('independent reference checksum mismatch')
        with TwoBitReader(args.reference) as genome:
            for record in by_id.values():
                if genome.fetch(record['chromosome'],record['start'],record['end'],strand=record['source_metadata']['source_declared_strand'])!=record['sequence']:raise ValueError('admitted core fails independent reference replay')
        reference_checked=True
    reference=json.loads((root/'reference-validation.json').read_text())
    if reference['records_checked']!=len(by_id) or reference['exact_forward_matches']+reference['reverse_complement_matches']!=len(by_id):raise ValueError('reference receipt coverage mismatch')
    print(json.dumps({'verified_offline':True,'endpoint':panel['endpoint'],'records':len(by_id),'source_ess_vectors_rechecked':len(seen_scores),'source_designs_rechecked':len(seen_design),'full_reference_replayed':reference_checked,'manifest_sha256':digest((root/'manifest.json').read_bytes())}))


def self_test():
    with tempfile.TemporaryDirectory() as directory:
        path=Path(directory)/'test.xlsx'
        with zipfile.ZipFile(path,'w') as archive:
            archive.writestr('xl/workbook.xml','<workbook xmlns="'+NS['m']+'" xmlns:r="'+REL+'"><sheets><sheet name="example" sheetId="1" r:id="rId1"/></sheets></workbook>')
            archive.writestr('xl/_rels/workbook.xml.rels','<Relationships><Relationship Id="rId1" Target="worksheets/sheet1.xml"/></Relationships>')
            archive.writestr('xl/sharedStrings.xml','<sst xmlns="'+NS['m']+'"><si><t>ID</t></si><si><r><t>AC</t></r><r><t>GT</t></r></si></sst>')
            archive.writestr('xl/worksheets/sheet1.xml','<worksheet xmlns="'+NS['m']+'"><sheetData><row r="1"><c r="A1" t="s"><v>0</v></c><c r="B1" t="inlineStr"><is><t>DNA</t></is></c></row><row r="2"><c r="A2" t="n"><v>12</v></c><c r="B2" t="s"><v>1</v></c><c r="C2" t="e"><v>#N/A</v></c></row></sheetData></worksheet>')
        with Workbook(path) as workbook:
            rows=list(workbook.rows('example'));assert rows[0][1]=={'A':'ID','B':'DNA'};assert rows[1][1]['B']=='ACGT';assert rows[1][1]['C']=={'excel_error':'#N/A'}
        with Workbook(path) as workbook:assert list(workbook.rows('example',columns={'A'},stop_after=1))==[(1,{'A':'ID'})]
    assert chromosome('1')=='chr1' and integer('200.0')==200
    print('Strict XLSX metadata/type handling checks passed.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);commands=parser.add_subparsers(dest='command',required=True)
    p=commands.add_parser('design-only');p.add_argument('directory');p.add_argument('--protocol',required=True);p.add_argument('--expected-protocol-sha256',required=True);p.add_argument('--reference',default='/workspace/.tools/reference/hg38.2bit')
    p=commands.add_parser('build-ess');p.add_argument('directory');p.add_argument('--protocol',required=True);p.add_argument('--expected-protocol-sha256',required=True);p.add_argument('--reference',default='/workspace/.tools/reference/hg38.2bit')
    p=commands.add_parser('verify');p.add_argument('directory');p.add_argument('--reference',help='Optional full independent hg38 replay; serialized reference receipt alone is not called a replay')
    commands.add_parser('self-test');args=parser.parse_args()
    if args.command=='design-only':design_only(args)
    elif args.command=='build-ess':build_ess(args)
    elif args.command=='verify':verify_panel(args)
    else:self_test()
