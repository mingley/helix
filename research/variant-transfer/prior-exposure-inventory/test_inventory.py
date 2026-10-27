"""Independent controls for fixed-DNA inventory integrity and window factorization."""
import collections,gzip,hashlib,itertools,json,pathlib,unittest
ROOT=pathlib.Path(__file__).resolve().parent
BUNDLE=ROOT/'inventory-v1'
def rows(name):
 with gzip.open(BUNDLE/name,'rt') as stream:return [json.loads(line) for line in stream]
def hamming(a,b):return sum(x!=y for x,y in zip(a,b))
def reverse(s):return s.translate(str.maketrans('ACGT','TGCA'))[::-1]

class ExactFactorizationTests(unittest.TestCase):
 def test_triangle_shortcut_equals_all_exposed_mutants_exhaustively(self):
  reference='ACGTA';exposed={0:{'T'},2:{'A','C'},4:{'C'}};threshold=1
  old=[reference]
  for pos,alts in exposed.items():
   old.extend(reference[:pos]+alt+reference[pos+1:] for alt in alts)
  for chars in itertools.product('ACGT',repeat=len(reference)):
   query=''.join(chars);distance=hamming(query,reference)
   factorized=distance<=threshold or (distance==threshold+1 and any(query[pos] in alts and query[pos]!=reference[pos] for pos,alts in exposed.items()))
   explicit=any(hamming(query,sequence)<=threshold for sequence in old)
   self.assertEqual(factorized,explicit,query)
 def test_reverse_complement_must_reverse_and_complement_edit(self):
  reference='AACCG';mutant='AATCG';query=reverse(mutant)
  self.assertEqual(hamming(query,reverse(reference)),1)
  # The forward C>T edit at index2 becomes G>A at reverse index2.
  self.assertEqual(query[2],'A');self.assertEqual(reverse(reference)[2],'G')
 def test_deletion_windows_need_explicit_representation(self):
  reference='AACGTT';deleted=reference[:2]+reference[3:]
  query=deleted[1:5]
  self.assertNotIn(query,[reference[i:i+4] for i in range(3)])
  self.assertEqual(query,'AGTT')

class FrozenInventoryTests(unittest.TestCase):
 def test_all_bytes_match_retained_manifest(self):
  manifest=json.loads((BUNDLE/'manifest.json').read_text())
  expected=set(manifest['files'])|{'manifest.json'}
  self.assertEqual(expected,{p.name for p in BUNDLE.iterdir()})
  for name,item in manifest['files'].items():
   payload=(BUNDLE/name).read_bytes()
   self.assertEqual(len(payload),item['bytes'],name)
   self.assertEqual(hashlib.sha256(payload).hexdigest(),item['sha256'],name)
 def test_factorized_source_alleles_cover_every_native_snv_window(self):
  edits=rows('kircher-native-edits.jsonl.gz');observed=collections.defaultdict(lambda:collections.defaultdict(set))
  for e in edits:
   if e['edit_type']=='single_nucleotide_substitution':observed[e['element_id']][e['offset_0based']].add(e['alternate'])
  for window in rows('kircher-reference-windows200.jsonl.gz'):
   start=window['parent_offset_0based'];actual={int(p):set(alts) for p,alts in window['observed_snv_alts']}
   expected={p-start:alts for p,alts in observed[window['element_id']].items() if 0<=p-start<200}
   self.assertEqual(actual,expected,window['id'])
 def test_native_edits_are_exact_single_edits_of_source_reference(self):
  loci={q['id']:q for q in json.loads((BUNDLE/'kircher-loci.json').read_text())['loci']}
  for e in rows('kircher-native-edits.jsonl.gz'):
   reference=loci[e['locus_id']]['sequence'];p=e['offset_0based'];self.assertEqual(reference[p],e['reference'])
   replacement='' if e['edit_type']=='single_nucleotide_deletion' else e['alternate']
   self.assertEqual(e['sequence'],reference[:p]+replacement+reference[p+1:])
   self.assertEqual(hashlib.sha256(e['sequence'].encode()).hexdigest(),e['sequence_sha256'])
   self.assertTrue(e['sources'])
 def test_counts_reconcile_to_source_exposure_without_qc_filter(self):
  summary=json.loads((BUNDLE/'summary.json').read_text());bindings=json.loads((BUNDLE/'molecular-identity-bindings.json').read_text())['bindings']
  self.assertEqual(summary['kircher_explored_mavedb_rows'],18468)
  self.assertEqual(sum(b['source_rows'] for b in bindings),18468)
  self.assertEqual(sum(b['representation_counts']['single_nucleotide_substitution'] for b in bindings),15699)
  self.assertEqual(sum(b['representation_counts']['single_nucleotide_deletion'] for b in bindings),2769)
  self.assertTrue(all(b['molecular_identity_bijection'] and not b['numeric_outcomes_compared'] for b in bindings))
  self.assertEqual(summary['agarwal_five_fold_training_union'],55276)

if __name__=='__main__':unittest.main()
