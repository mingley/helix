#!/usr/bin/env python3
"""Acquire bounded old Kircher author tables for DNA identity only; no new-study outcomes."""
import argparse,concurrent.futures,datetime,hashlib,json,pathlib,urllib.request
ROOT=pathlib.Path(__file__).resolve().parents[3]
CACHE=pathlib.Path('/workspace/.tools/prior-exposure-author-tables')
NAMES={'GRCh38_GP1BA.tsv','GRCh38_F9.tsv','GRCh38_FOXE1.tsv','GRCh38_HBB.tsv','GRCh38_HBG1.tsv','GRCh38_HNF4A.tsv','GRCh38_LDLR.2.tsv','GRCh38_LDLR.tsv','GRCh38_MSMB.tsv','GRCh38_PKLR-24h.tsv','GRCh38_PKLR-48h.tsv','GRCh38_TERT-GBM.tsv','GRCh38_TERT-GAa.tsv','GRCh38_TERT-GSc.tsv','GRCh38_TERT-HEK.tsv'}
def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--cache',type=pathlib.Path,default=CACHE)
 parser.add_argument('--receipt',type=pathlib.Path,default=pathlib.Path(__file__).with_name('author-acquisition.json'))
 args=parser.parse_args();cache=args.cache;dest=args.receipt
 if dest.exists():raise ValueError('receipt exists; choose a fresh receipt path')
 cache.mkdir(exist_ok=True,parents=True)
 snapshot=ROOT/'research/forensic-sources/osf-grch38-files.json'
 metadata=json.loads(snapshot.read_text())
 rows=[r for r in metadata['data'] if r['attributes']['name'] in NAMES]
 if {r['attributes']['name'] for r in rows}!=NAMES:raise ValueError('Incomplete old source inventory')
 def fetch(r):
  a=r['attributes'];name=a['name'];expected=a['extra']['hashes']['sha256'];dst=cache/name;url=r['links']['download']
  if dst.exists():raw=dst.read_bytes();status='verified_existing_cache';content_type=None;final_url=url
  elif name=='GRCh38_HBB.tsv':raw=(ROOT/'research/forensic-sources/GRCh38_HBB.tsv').read_bytes();status='verified_retained_old_source';content_type=None;final_url=url
  else:
   req=urllib.request.Request(url,headers={'User-Agent':'Helix-source-identity/1.0'})
   with urllib.request.urlopen(req,timeout=40) as response:
    raw=response.read(1000001);status=response.status;content_type=response.headers.get('Content-Type');final_url=response.url
    if len(raw)>1000000:raise ValueError('old table exceeds1MB cap')
  sha=hashlib.sha256(raw).hexdigest()
  if len(raw)!=a['size'] or sha!=expected:raise ValueError('old source checksum/size mismatch:'+name)
  if not raw.startswith(b'Chromosome\tPosition\tRef\tAlt\t'):raise ValueError('unexpected old table schema')
  if not dst.exists():dst.write_bytes(raw)
  return {'name':name,'url':url,'final_url':final_url,'bytes':len(raw),'sha256':sha,'osf_version':a['current_version'],'status':status,'content_type':content_type,'cache_path':str(dst),'role':'Old explored-study molecular identity forensics only. Numeric outcome columns not inspected or used.'}
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:result=list(pool.map(fetch,rows))
 out={'schema_version':1,'retrieved_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'upstream_metadata_path':str(snapshot.relative_to(ROOT)),'upstream_metadata_sha256':hashlib.sha256(snapshot.read_bytes()).hexdigest(),'files':sorted(result,key=lambda a:a['name']),'license_note':'Unchanged old author OSF deposit, not relabeled CC0. Original producer-use statement and article citation retained under research/forensic-sources. Raw cache is outside Git; derived DNA identity facts retain this receipt.','new_study_exposure':'No Siraj S4 outcome file requested or read.'}
 dest.write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({'verified_tables':len(result),'total_bytes':sum(r['bytes'] for r in result),'receipt':str(dest)}))
if __name__=='__main__':main()
