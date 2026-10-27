"""Small executable counterexamples for the independent C++ grouping implementation."""
import csv,pathlib,subprocess,tempfile,unittest
HERE=pathlib.Path(__file__).resolve().parent

def rev(s):return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
def mutate(s,positions,base):
 a=list(s)
 for p in positions:a[p]=base
 return ''.join(a)
def row(identity,chrom,start,sequence):return [identity,chrom,str(start),str(start+200),identity,sequence,mutate(sequence,[99],'G' if sequence[99]!='G' else 'A')]

class ExecutableCounterexamples(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  cls.temp=tempfile.TemporaryDirectory(prefix='helix-grouping-controls-');cls.binary=pathlib.Path(cls.temp.name)/'groups'
  subprocess.run(['g++','-O3','-std=c++20','-Wall','-Wextra','-pedantic',str(HERE/'group.cpp'),'-o',str(cls.binary)],check=True,capture_output=True)
 @classmethod
 def tearDownClass(cls):cls.temp.cleanup()
 def run_case(self,new,old=(),guards=(),controls=()):
  with tempfile.TemporaryDirectory(dir=self.temp.name) as directory:
   p=pathlib.Path(directory)
   for name,values in [('new.tsv',new),('old.tsv',old),('guards.tsv',guards)]:
    (p/name).write_text(''.join('\t'.join(v)+'\n' for v in values))
   (p/'controls.tsv').write_text('\n'.join(controls))
   result=subprocess.run([str(self.binary),*(str(p/n) for n in ['new.tsv','old.tsv','guards.tsv','controls.tsv','partition.tsv','events.tsv'])],capture_output=True,text=True)
   self.assertEqual(result.returncode,0,result.stderr)
   with (p/'partition.tsv').open() as source:
    return {v['id']:v for v in csv.DictReader(source,delimiter='\t')}
 def test_start_bin_does_not_union_a_crossing_interval_with_a_separated_next_bin(self):
  q=self.run_case([row('a','chr1',999950,'A'*200),row('b','chr1',1000200,'C'*200)])
  self.assertNotEqual(q['a']['root_index'],q['b']['root_index'])
 def test_adjacency_across_start_bin_boundary_is_not_positive_overlap(self):
  q=self.run_case([row('a','chr1',999950,'A'*200),row('b','chr1',1000150,'C'*200)])
  self.assertNotEqual(q['a']['root_index'],q['b']['root_index'])
 def test_actual_overlap_joins_across_start_bin_boundary(self):
  q=self.run_case([row('a','chr1',999950,'A'*200),row('b','chr1',1000149,'C'*200)])
  self.assertEqual(q['a']['root_index'],q['b']['root_index'])
 def test_distant_same_start_bin_still_joins(self):
  q=self.run_case([row('a','chr1',1,'A'*200),row('b','chr1',900000,'C'*200)])
  self.assertEqual(q['a']['root_index'],q['b']['root_index'])
 def test_chromosome_is_part_of_bin_key(self):
  q=self.run_case([row('a','chr1',1,'A'*200),row('b','chr2',1,'C'*200)])
  self.assertNotEqual(q['a']['root_index'],q['b']['root_index'])
 def test_twenty_distributed_mismatches_found(self):
  changed=mutate('A'*200,[j*200//21 for j in range(20)],'C')
  q=self.run_case([row('a','chr1',1,'A'*200),row('b','chr2',1,changed)])
  self.assertEqual(q['a']['root_index'],q['b']['root_index'])
 def test_twenty_one_mismatches_do_not_qualify_as_new_homology(self):
  changed=mutate('A'*200,[j*200//21 for j in range(21)],'C')
  q=self.run_case([row('a','chr1',1,'A'*200),row('b','chr2',1,changed)])
  self.assertNotEqual(q['a']['root_index'],q['b']['root_index'])
 def test_reverse_strand_homology_joins(self):
  changed=mutate('A'*200,list(range(10)),'C')
  q=self.run_case([row('a','chr1',1,'A'*200),row('b','chr2',1,rev(changed))])
  self.assertEqual(q['a']['root_index'],q['b']['root_index'])
 def test_guard_boundary_is_inclusive(self):
  guard=[['explored_Agarwal_insert','chr1','0','200','old']]
  q=self.run_case([row('a','chr1',10200,'A'*200)],guards=guard);self.assertEqual(q['a']['retained_prior'],'0')
  q=self.run_case([row('a','chr1',10201,'A'*200)],guards=guard);self.assertEqual(q['a']['retained_prior'],'1')
 def test_whole_component_exclusion_propagates(self):
  q=self.run_case([row('a','chr1',1,'A'*200),row('b','chr1',900000,'C'*200)],guards=[['explored_Kircher_locus','chr1','0','200','old']])
  self.assertEqual([q[x]['retained_prior'] for x in ['a','b']],['0','0'])
 def test_factorized_twenty_one_mismatch_repair_uses_twenty_two_pieces(self):
  changed=mutate('A'*200,[j*200//22 for j in range(21)],'C')
  q=self.run_case([row('a','chr1',1,changed)],old=[['factor','old','A'*200,'0:C']]);self.assertEqual(q['a']['retained_prior'],'0')
  q=self.run_case([row('a','chr1',1,changed)],old=[['factor','old','A'*200,'0:T']]);self.assertEqual(q['a']['retained_prior'],'1')
 def test_explicit_old_hits_exclude_without_joining_new_components(self):
  q=self.run_case([row('a','chr1',1,'A'*200),row('b','chr2',1,'C'*200)],old=[['explicit','oldA','A'*200],['explicit','oldC','C'*200]])
  self.assertNotEqual(q['a']['root_index'],q['b']['root_index']);self.assertEqual(q['a']['retained_prior'],'0');self.assertEqual(q['b']['retained_prior'],'0')
 def test_technical_control_flag_follows_full_component(self):
  q=self.run_case([row('a','chr1',1,'A'*200),row('b','chr1',900000,'C'*200)],controls=['a'])
  self.assertEqual([q[x]['technical_control_component'] for x in ['a','b']],['1','1']);self.assertEqual(q['a']['retained_prior'],'1')
 def test_source_parent_can_join_different_bins(self):
  a=row('a','chr1',1,'A'*200);b=row('b','chr2',1,'C'*200);a[4]=b[4]='source-parent'
  q=self.run_case([a,b]);self.assertEqual(q['a']['root_index'],q['b']['root_index'])

if __name__=='__main__':unittest.main()
