/- GC fraction and reverse-complement-canonical 3-mer counts.

For a sequence of at least four bases, each interior base occurs in three
complete overlapping 3-mers.  The first and last bases occur once each, and
the second and penultimate bases occur twice each.  The theorems below prove
3*GC = weightedTrimerGC + 2*GC(first) + GC(second) + GC(penultimate)
        + 2*GC(last).

The proof is exact, uses integers and finite lists, and is outcome-free.  The
same identity first holds for arbitrary integer base weights.  We then define
DNA, reverse complementation, trimer windows, and any choice of one member of
a trimer/reverse-complement pair.  Reverse complementation preserves the GC
weight.  Equality of all canonical-trimer counts and the boundary GC mass
therefore forces equal GC counts.  A six-base counterexample shows that GC is
not globally redundant when the boundary changes.

The model's fixed-length feature denominators and frozen affine preprocessing
must be incorporated separately into weights.  The score algebra here proves
integer equalities, not the behavior of floating-point execution or empirical
biological effects. No user-declared axioms, sorry, native_decide, or external
packages are used. Lean's `#print axioms` reports `propext` and `Quot.sound`
for the main identity; the count-to-permutation result also uses
`Classical.choice` from the standard library.
-/
import Std
set_option autoImplicit false

namespace HelixGCTrimerIdentity

def total : List Int → Int
  | [] => 0
  | x :: xs => x + total xs

def tripleTotal : List Int → Int
  | [] => 0
  | x :: xs =>
    match xs with
    | y :: z :: _ => x + y + z + tripleTotal xs
    | _ => 0

theorem tripleTotal_step (x y z : Int) (xs : List Int) :
    tripleTotal (x :: y :: z :: xs) = x + y + z + tripleTotal (y :: z :: xs) := by
  rfl

theorem boundary_identity (middle : List Int) :
    ∀ (a b c d : Int),
      3 * total (a :: b :: (middle ++ [c, d])) =
        tripleTotal (a :: b :: (middle ++ [c, d])) + 2*a + b + c + 2*d := by
  induction middle with
  | nil =>
    intro a b c d
    simp only [List.nil_append, total, tripleTotal]
    omega
  | cons x xs ih =>
    intro a b c d
    have h := ih b x c d
    simp only [List.cons_append] at h ⊢
    rw [tripleTotal_step]
    simp only [total] at h ⊢
    omega

-- Every list of length at least two has a penultimate and final element.
theorem split_last_two {α : Type} (xs : List α) (h : 2 ≤ xs.length) :
    ∃ (middle : List α) (c d : α), xs = middle ++ [c, d] := by
  induction xs with
  | nil => simp at h
  | cons x xs ih =>
    cases xs with
    | nil => simp at h
    | cons y ys =>
      cases ys with
      | nil => exact ⟨[], x, y, rfl⟩
      | cons z zs =>
        obtain ⟨middle, c, d, eq⟩ := ih (by simp)
        exact ⟨x :: middle, c, d, by simp [eq]⟩

theorem every_long_list (xs : List Int) (h : 4 ≤ xs.length) :
    ∃ (a b c d : Int) (middle : List Int),
      xs = a :: b :: (middle ++ [c, d]) ∧
      3 * total xs = tripleTotal xs + 2*a + b + c + 2*d := by
  cases xs with
  | nil => simp at h
  | cons a xs =>
    cases xs with
    | nil => simp at h
    | cons b rest =>
      have hRest : 2 ≤ rest.length := by simp at h; omega
      obtain ⟨middle, c, d, eq⟩ := split_last_two rest hRest
      subst rest
      exact ⟨a, b, c, d, middle, rfl, boundary_identity middle a b c d⟩

theorem linear_reparameterization (g t edge a b : Int)
    (h : 3*g = t + edge) :
    3*(a*t + b*g) = (3*a + b)*t + b*edge := by
  calc
    3*(a*t + b*g) = 3*a*t + b*(3*g) := by
      simp only [Int.mul_add]
      ac_rfl
    _ = 3*a*t + b*(t+edge) := by rw [h]
    _ = (3*a+b)*t + b*edge := by
      simp only [Int.mul_add, Int.add_mul]
      ac_rfl

inductive Base where
  | A | C | G | T
  deriving DecidableEq, Repr

def gc : Base → Int
  | .A => 0
  | .C => 1
  | .G => 1
  | .T => 0

def complement : Base → Base
  | .A => .T
  | .C => .G
  | .G => .C
  | .T => .A

theorem gc_complement (b : Base) : gc (complement b) = gc b := by
  cases b <;> rfl

structure Trimer where
  first : Base
  second : Base
  third : Base
  deriving DecidableEq, Repr

def reverseComplement (t : Trimer) : Trimer :=
  ⟨complement t.third, complement t.second, complement t.first⟩

def trimerGC (t : Trimer) : Int := gc t.first + gc t.second + gc t.third

theorem trimer_gc_reverse_complement (t : Trimer) :
    trimerGC (reverseComplement t) = trimerGC t := by
  cases t with
  | mk a b c =>
    simp only [trimerGC, reverseComplement, gc_complement]
    omega

-- This holds for lexicographic RC canonicalization and any other pair choice.
theorem canonical_trimer_gc (t canonical : Trimer)
    (choice : canonical = t ∨ canonical = reverseComplement t) :
    trimerGC canonical = trimerGC t := by
  rcases choice with same | reverse
  · rw [same]
  · rw [reverse, trimer_gc_reverse_complement]

def genomicGC (dna : List Base) : Int := total (dna.map gc)
def slidingGC (dna : List Base) : Int := tripleTotal (dna.map gc)

theorem dna_boundary_identity (middle : List Base) (a b c d : Base) :
    3 * genomicGC (a :: b :: (middle ++ [c,d])) =
      slidingGC (a :: b :: (middle ++ [c,d])) + 2*gc a + gc b + gc c + 2*gc d := by
  simpa only [genomicGC, slidingGC, List.map_cons, List.map_append, List.map_nil]
    using boundary_identity (middle.map gc) (gc a) (gc b) (gc c) (gc d)

theorem equal_trimer_and_edges_imply_equal_gc
    (m1 m2 : List Base) (a b c d : Base)
    (same : slidingGC (a :: b :: (m1 ++ [c,d])) =
      slidingGC (a :: b :: (m2 ++ [c,d])) ) :
    genomicGC (a :: b :: (m1 ++ [c,d])) =
      genomicGC (a :: b :: (m2 ++ [c,d])) := by
  have h1 := dna_boundary_identity m1 a b c d
  have h2 := dna_boundary_identity m2 a b c d
  omega

def windows : List Base → List Trimer
  | [] => []
  | x :: xs =>
    match xs with
    | y :: z :: _ => ⟨x, y, z⟩ :: windows xs
    | _ => []

def sumTrimerGC : List Trimer → Int
  | [] => 0
  | x :: xs => trimerGC x + sumTrimerGC xs

theorem windows_sum_is_sliding_gc (dna : List Base) :
    sumTrimerGC (windows dna) = slidingGC dna := by
  induction dna with
  | nil => rfl
  | cons a tail ih =>
    cases tail with
    | nil => rfl
    | cons b rest =>
      cases rest with
      | nil => rfl
      | cons c rest =>
        change trimerGC ⟨a,b,c⟩ + sumTrimerGC (windows (b :: c :: rest)) =
          gc a + gc b + gc c + slidingGC (b :: c :: rest)
        simp only [trimerGC]
        rw [ih]

theorem canonical_windows_preserve_gc
    (xs : List Trimer) (canonical : Trimer → Trimer)
    (choice : ∀ t, canonical t = t ∨ canonical t = reverseComplement t) :
    sumTrimerGC (xs.map canonical) = sumTrimerGC xs := by
  induction xs with
  | nil => rfl
  | cons t rest ih =>
    simp only [List.map_cons, sumTrimerGC]
    rw [canonical_trimer_gc t (canonical t) (choice t), ih]

theorem canonical_dna_boundary_identity
    (middle : List Base) (a b c d : Base) (canonical : Trimer → Trimer)
    (choice : ∀ t, canonical t = t ∨ canonical t = reverseComplement t) :
    3 * genomicGC (a :: b :: (middle ++ [c,d])) =
      sumTrimerGC ((windows (a :: b :: (middle ++ [c,d]))).map canonical) +
        2*gc a + gc b + gc c + 2*gc d := by
  rw [canonical_windows_preserve_gc _ canonical choice, windows_sum_is_sliding_gc]
  exact dna_boundary_identity middle a b c d

-- Equal counts imply a permutation, and hence equal weighted trimer sums.
theorem sum_trimer_gc_perm {xs ys : List Trimer} (permutation : xs.Perm ys) :
    sumTrimerGC xs = sumTrimerGC ys := by
  induction permutation with
  | nil => rfl
  | cons x _ ih => simp only [sumTrimerGC, ih]
  | swap x y rest => simp only [sumTrimerGC]; omega
  | trans _ _ ih1 ih2 => exact Eq.trans ih1 ih2

theorem equal_canonical_counts_and_boundary_imply_equal_gc
    (m1 m2 : List Base) (a b c d : Base) (canonical : Trimer → Trimer)
    (choice : ∀ t, canonical t = t ∨ canonical t = reverseComplement t)
    (counts : ∀ t, List.count t ((windows (a :: b :: (m1 ++ [c,d]))).map canonical) =
      List.count t ((windows (a :: b :: (m2 ++ [c,d]))).map canonical)) :
    genomicGC (a :: b :: (m1 ++ [c,d])) = genomicGC (a :: b :: (m2 ++ [c,d])) := by
  have eqSum := sum_trimer_gc_perm (List.perm_iff_count.mpr counts)
  have h1 := canonical_dna_boundary_identity m1 a b c d canonical choice
  have h2 := canonical_dna_boundary_identity m2 a b c d canonical choice
  omega

theorem canonical_all_long_dna (dna : List Base) (h : 4 ≤ dna.length)
    (canonical : Trimer → Trimer)
    (choice : ∀ t, canonical t = t ∨ canonical t = reverseComplement t) :
    ∃ (a b c d : Base) (middle : List Base),
      dna = a :: b :: (middle ++ [c,d]) ∧
      3 * genomicGC dna = sumTrimerGC ((windows dna).map canonical) +
        2*gc a + gc b + gc c + 2*gc d := by
  cases dna with
  | nil => simp at h
  | cons a dna =>
    cases dna with
    | nil => simp at h
    | cons b rest =>
      have hRest : 2 ≤ rest.length := by simp at h; omega
      obtain ⟨middle, c, d, eq⟩ := split_last_two rest hRest
      subst rest
      exact ⟨a, b, c, d, middle, rfl, canonical_dna_boundary_identity middle a b c d canonical choice⟩

theorem equal_canonical_counts_and_edge_mass_imply_equal_gc
    (m1 m2 : List Base) (a1 b1 c1 d1 a2 b2 c2 d2 : Base)
    (canonical : Trimer → Trimer)
    (choice : ∀ t, canonical t = t ∨ canonical t = reverseComplement t)
    (counts : ∀ t, List.count t ((windows (a1 :: b1 :: (m1 ++ [c1,d1]))).map canonical) =
      List.count t ((windows (a2 :: b2 :: (m2 ++ [c2,d2]))).map canonical))
    (sameEdge : 2*gc a1 + gc b1 + gc c1 + 2*gc d1 =
      2*gc a2 + gc b2 + gc c2 + 2*gc d2) :
    genomicGC (a1 :: b1 :: (m1 ++ [c1,d1])) = genomicGC (a2 :: b2 :: (m2 ++ [c2,d2])) := by
  have eqSum := sum_trimer_gc_perm (List.perm_iff_count.mpr counts)
  have h1 := canonical_dna_boundary_identity m1 a1 b1 c1 d1 canonical choice
  have h2 := canonical_dna_boundary_identity m2 a2 b2 c2 d2 canonical choice
  omega

-- The trimer score here may be arbitrary; t is its GC-weighted trimer sum.
theorem general_score_reparameterization (trimerScore g t edge coefficient : Int)
    (h : 3*g = t + edge) :
    3*(trimerScore + coefficient*g) =
      3*trimerScore + coefficient*t + coefficient*edge := by
  calc
    3*(trimerScore + coefficient*g) =
        3*trimerScore + coefficient*(3*g) := by
      simp only [Int.mul_add]
      ac_rfl
    _ = 3*trimerScore + coefficient*(t+edge) := by rw [h]
    _ = 3*trimerScore + coefficient*t + coefficient*edge := by
      simp only [Int.mul_add]
      omega

def edgeExample1 : List Base := [.A, .A, .A, .C, .A, .A]
def edgeExample2 : List Base := [.A, .C, .A, .A, .A, .C]

theorem identical_trimer_counts_can_have_different_gc (canonical : Trimer → Trimer) :
    (∀ t, List.count t ((windows edgeExample1).map canonical) =
      List.count t ((windows edgeExample2).map canonical)) ∧
    genomicGC edgeExample1 ≠ genomicGC edgeExample2 := by
  have permutation : (windows edgeExample1).Perm (windows edgeExample2) := by decide
  constructor
  · intro t
    exact (permutation.map canonical).count_eq t
  · decide

end HelixGCTrimerIdentity
