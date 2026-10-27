/- Exact locality of fixed linear DNA-trimer + base-composition predictors.
   Integer weights encode rational weights after common-denominator scaling.
   For fixed length, per-feature frequency denominators and frozen training
   centering/scaling can be absorbed into the intercept and weights. That
   algebraic bridge is described in the study; this file proves the resulting
   exact finite-prefix score identities, not Rust/IEEE-754 rounding behavior.
   GC is the special baseWeights = gcWeights coefficient case.
   This proves an architecture restriction, not a property of genomic biology.
   No axioms, sorry, native_decide, source outcomes, or external packages. -/
import Std
set_option autoImplicit false

namespace HelixVariantLocality

inductive Base where
  | A | C | G | T
  deriving DecidableEq, Repr

abbrev DNA := Nat → Base

structure Trimer where
  first : Base
  second : Base
  third : Base
  deriving DecidableEq, Repr

def mutate (dna : DNA) (position : Nat) (replacement : Base) : DNA :=
  fun index => if index = position then replacement else dna index

def window (dna : DNA) (start : Nat) : Trimer :=
  ⟨dna start, dna (start + 1), dna (start + 2)⟩

def inWindow (start position : Nat) : Prop :=
  start ≤ position ∧ position ≤ start + 2

instance (start position : Nat) : Decidable (inWindow start position) :=
  inferInstanceAs (Decidable (start ≤ position ∧ position ≤ start + 2))

theorem three_possible_affected_starts (start position : Nat)
    (inside : inWindow start position) :
    start = position ∨ start + 1 = position ∨ start + 2 = position := by
  unfold inWindow at inside
  omega

theorem window_unchanged_outside (dna : DNA) (start position : Nat)
    (replacement : Base) (outside : ¬ inWindow start position) :
    window (mutate dna position replacement) start = window dna start := by
  have h0 : start ≠ position := by unfold inWindow at outside; omega
  have h1 : start + 1 ≠ position := by unfold inWindow at outside; omega
  have h2 : start + 2 ≠ position := by unfold inWindow at outside; omega
  simp [window, mutate, h0, h1, h2]

theorem mutate_commute (dna : DNA) (p q : Nat) (bp bq : Base)
    (different : p ≠ q) :
    mutate (mutate dna p bp) q bq = mutate (mutate dna q bq) p bp := by
  funext index
  by_cases hp : index = p
  · subst index
    simp [mutate, different]
  · by_cases hq : index = q
    · subst index
      simp [mutate, hp]
    · simp [mutate, hp, hq]

theorem separated_positions_share_no_window (start p q : Nat)
    (separation : p + 3 ≤ q) :
    ¬ (inWindow start p ∧ inWindow start q) := by
  unfold inWindow
  omega

def sumAt : List Nat → (Nat → Int) → Int
  | [], _ => 0
  | index :: rest, value => value index + sumAt rest value

theorem sumAt_pair_identity (indices : List Nat) (a b c d : Nat → Int)
    (pointwise : ∀ index, a index + b index = c index + d index) :
    sumAt indices a + sumAt indices b = sumAt indices c + sumAt indices d := by
  induction indices with
  | nil => rfl
  | cons index rest ih =>
    simp only [sumAt]
    have atIndex := pointwise index
    omega

theorem sumAt_sub (indices : List Nat) (a b : Nat → Int) :
    sumAt indices a - sumAt indices b = sumAt indices (fun index => a index - b index) := by
  induction indices with
  | nil => rfl
  | cons index rest ih => simp only [sumAt]; omega

theorem sumAt_filter_zeros (indices : List Nat) (value : Nat → Int)
    (keep : Nat → Bool) (outside : ∀ index, keep index = false → value index = 0) :
    sumAt indices value = sumAt (indices.filter keep) value := by
  induction indices with
  | nil => rfl
  | cons index rest ih =>
    cases retained : keep index with
    | false => simp [sumAt, List.filter_cons, retained, outside index retained, ih]
    | true => simp [sumAt, List.filter_cons, retained, ih]

theorem sumAt_append (left right : List Nat) (value : Nat → Int) :
    sumAt (left ++ right) value = sumAt left value + sumAt right value := by
  induction left with
  | nil => simp [sumAt]
  | cons index rest ih => simp only [List.cons_append, sumAt, ih]; omega

theorem sumAt_range_single (length position : Nat) (value : Int) :
    sumAt (List.range length) (fun index => if index = position then value else 0) =
      if position < length then value else 0 := by
  induction length with
  | zero => simp [sumAt]
  | succ length ih =>
    rw [List.range_succ, sumAt_append]
    simp only [sumAt, Int.add_zero, ih]
    by_cases inside : position < length
    · have ne : length ≠ position := by omega
      have larger : position < length + 1 := by omega
      simp [inside, ne, larger]
    · by_cases equal : length = position
      · subst position; simp
      · have beyond : ¬ position < length + 1 := by omega
        simp [inside, equal, beyond]

def trimerScore (weights : Trimer → Int) (dna : DNA) (start : Nat) : Int :=
  weights (window dna start)

def baseScore (weights : Base → Int) (dna : DNA) (position : Nat) : Int :=
  weights (dna position)

def affectedStarts (length position : Nat) : List Nat :=
  (List.range (length - 2)).filter (fun start => decide (inWindow start position))

-- This is a full finite-window sum identity, not an assumed additive SNV model.
-- Every surviving summand starts at p, p-1 or p-2 when that start is valid.
theorem trimer_mutation_delta_is_local (length position : Nat) (weights : Trimer → Int)
    (dna : DNA) (replacement : Base) :
    sumAt (List.range (length - 2)) (trimerScore weights (mutate dna position replacement)) -
      sumAt (List.range (length - 2)) (trimerScore weights dna) =
    sumAt (affectedStarts length position) (fun start =>
      trimerScore weights (mutate dna position replacement) start -
        trimerScore weights dna start) := by
  rw [sumAt_sub]
  apply sumAt_filter_zeros
  intro start outside
  have notInside : ¬ inWindow start position := by simpa using outside
  simp [trimerScore, window_unchanged_outside dna start position replacement notInside]

theorem base_mutation_delta (length position : Nat) (weights : Base → Int)
    (dna : DNA) (replacement : Base) :
    sumAt (List.range length) (baseScore weights (mutate dna position replacement)) -
      sumAt (List.range length) (baseScore weights dna) =
    if position < length then weights replacement - weights (dna position) else 0 := by
  rw [sumAt_sub]
  have pointwise : (fun index => baseScore weights (mutate dna position replacement) index -
      baseScore weights dna index) =
      (fun index => if index = position then weights replacement - weights (dna position) else 0) := by
    funext index
    by_cases equal : index = position
    · subst index; simp [baseScore, mutate]
    · simp [baseScore, mutate, equal]
  rw [pointwise, sumAt_range_single]

def gcWeights (coefficient : Int) : Base → Int
  | .A => 0
  | .C => coefficient
  | .G => coefficient
  | .T => 0

-- Any trimer weight function includes reverse-complement canonical weights:
-- choose weights(t) = coefficient(canonical(t)). The proof needs no orientation.
theorem trimer_pair_identity (weights : Trimer → Int) (dna : DNA)
    (start p q : Nat) (bp bq : Base) (separation : p + 3 ≤ q) :
    trimerScore weights (mutate (mutate dna p bp) q bq) start +
      trimerScore weights dna start =
    trimerScore weights (mutate dna p bp) start +
      trimerScore weights (mutate dna q bq) start := by
  unfold trimerScore
  by_cases qOutside : ¬ inWindow start q
  · rw [window_unchanged_outside _ start q bq qOutside,
        window_unchanged_outside dna start q bq qOutside]
  · have pOutside : ¬ inWindow start p := by
      intro pInside
      exact separated_positions_share_no_window start p q separation
        ⟨pInside, by unfold inWindow at qOutside ⊢; omega⟩
    have different : p ≠ q := by omega
    rw [mutate_commute dna p q bp bq different,
        window_unchanged_outside _ start p bp pOutside,
        window_unchanged_outside dna start p bp pOutside]
    omega

theorem base_pair_identity (weights : Base → Int) (dna : DNA)
    (index p q : Nat) (bp bq : Base) (different : p ≠ q) :
    baseScore weights (mutate (mutate dna p bp) q bq) index +
      baseScore weights dna index =
    baseScore weights (mutate dna p bp) index +
      baseScore weights (mutate dna q bq) index := by
  by_cases hp : index = p
  · subst index
    simp [baseScore, mutate, different]
  · by_cases hq : index = q
    · subst index
      simp [baseScore, mutate, hp]
      omega
    · simp [baseScore, mutate, hp, hq]

def score (length : Nat) (intercept : Int) (trimerWeights : Trimer → Int)
    (baseWeights : Base → Int) (dna : DNA) : Int :=
  intercept + sumAt (List.range (length - 2)) (trimerScore trimerWeights dna) +
    sumAt (List.range length) (baseScore baseWeights dna)

theorem single_mutation_score_delta_is_local (length position : Nat) (intercept : Int)
    (trimerWeights : Trimer → Int) (baseWeights : Base → Int)
    (dna : DNA) (replacement : Base) :
    score length intercept trimerWeights baseWeights (mutate dna position replacement) -
      score length intercept trimerWeights baseWeights dna =
    sumAt (affectedStarts length position) (fun start =>
      trimerScore trimerWeights (mutate dna position replacement) start -
        trimerScore trimerWeights dna start) +
      (if position < length then baseWeights replacement - baseWeights (dna position) else 0) := by
  have trimerDelta := trimer_mutation_delta_is_local length position trimerWeights dna replacement
  have baseDelta := base_mutation_delta length position baseWeights dna replacement
  unfold score
  omega

-- The affected windows extend at most two bases on either side of the edit.
-- The inequality form also handles positions zero and one without underflow.
def sameNeighborhood (left right : DNA) (position : Nat) : Prop :=
  ∀ index, position ≤ index + 2 → index ≤ position + 2 → left index = right index

theorem same_neighborhood_gives_same_affected_window (left right : DNA)
    (position start : Nat) (agreement : sameNeighborhood left right position)
    (inside : inWindow start position) : window left start = window right start := by
  unfold inWindow at inside
  have h0 : left start = right start := agreement start (by omega) (by omega)
  have h1 : left (start + 1) = right (start + 1) :=
    agreement (start + 1) (by omega) (by omega)
  have h2 : left (start + 2) = right (start + 2) :=
    agreement (start + 2) (by omega) (by omega)
  simp [window, h0, h1, h2]

theorem matched_mutations_preserve_neighborhood (left right : DNA)
    (position : Nat) (replacement : Base) (agreement : sameNeighborhood left right position) :
    sameNeighborhood (mutate left position replacement) (mutate right position replacement) position := by
  intro index lower upper
  by_cases equal : index = position
  · simp [mutate, equal]
  · simp [mutate, equal, agreement index lower upper]

theorem same_neighborhood_gives_same_snv_effect (length position : Nat) (intercept : Int)
    (trimerWeights : Trimer → Int) (baseWeights : Base → Int) (left right : DNA)
    (replacement : Base) (agreement : sameNeighborhood left right position) :
    score length intercept trimerWeights baseWeights (mutate left position replacement) -
      score length intercept trimerWeights baseWeights left =
    score length intercept trimerWeights baseWeights (mutate right position replacement) -
      score length intercept trimerWeights baseWeights right := by
  rw [single_mutation_score_delta_is_local, single_mutation_score_delta_is_local]
  have sameTrimerDeltas : (fun start =>
      trimerScore trimerWeights (mutate left position replacement) start -
        trimerScore trimerWeights left start) = (fun start =>
      trimerScore trimerWeights (mutate right position replacement) start -
        trimerScore trimerWeights right start) := by
    funext start
    unfold trimerScore
    by_cases inside : inWindow start position
    · rw [same_neighborhood_gives_same_affected_window left right position start agreement inside,
          same_neighborhood_gives_same_affected_window _ _ position start
            (matched_mutations_preserve_neighborhood left right position replacement agreement) inside]
    · rw [window_unchanged_outside left start position replacement inside,
          window_unchanged_outside right start position replacement inside]
      omega
  have atPosition : left position = right position := agreement position (by omega) (by omega)
  rw [sameTrimerDeltas, atPosition]

-- Fixed sequence length and fixed trained weights are essential. Feature
-- centering, scaling and fixed denominators can be absorbed into these weights.
theorem distant_substitutions_are_additive (length : Nat) (intercept : Int)
    (trimerWeights : Trimer → Int) (baseWeights : Base → Int) (dna : DNA)
    (p q : Nat) (bp bq : Base) (separation : p + 3 ≤ q) :
    score length intercept trimerWeights baseWeights (mutate (mutate dna p bp) q bq) +
      score length intercept trimerWeights baseWeights dna =
    score length intercept trimerWeights baseWeights (mutate dna p bp) +
      score length intercept trimerWeights baseWeights (mutate dna q bq) := by
  have trimerIdentity := sumAt_pair_identity (List.range (length - 2))
    (trimerScore trimerWeights (mutate (mutate dna p bp) q bq))
    (trimerScore trimerWeights dna)
    (trimerScore trimerWeights (mutate dna p bp))
    (trimerScore trimerWeights (mutate dna q bq))
    (fun start => trimer_pair_identity trimerWeights dna start p q bp bq separation)
  have different : p ≠ q := by omega
  have baseIdentity := sumAt_pair_identity (List.range length)
    (baseScore baseWeights (mutate (mutate dna p bp) q bq))
    (baseScore baseWeights dna)
    (baseScore baseWeights (mutate dna p bp))
    (baseScore baseWeights (mutate dna q bq))
    (fun index => base_pair_identity baseWeights dna index p q bp bq different)
  unfold score
  omega

def epistasis (length : Nat) (intercept : Int) (trimerWeights : Trimer → Int)
    (baseWeights : Base → Int) (dna : DNA) (p q : Nat) (bp bq : Base) : Int :=
  score length intercept trimerWeights baseWeights (mutate (mutate dna p bp) q bq) -
    score length intercept trimerWeights baseWeights (mutate dna p bp) -
    score length intercept trimerWeights baseWeights (mutate dna q bq) +
    score length intercept trimerWeights baseWeights dna

theorem distant_epistasis_is_zero (length : Nat) (intercept : Int)
    (trimerWeights : Trimer → Int) (baseWeights : Base → Int) (dna : DNA)
    (p q : Nat) (bp bq : Base) (separation : p + 3 ≤ q) :
    epistasis length intercept trimerWeights baseWeights dna p q bp bq = 0 := by
  have identity := distant_substitutions_are_additive length intercept
    trimerWeights baseWeights dna p q bp bq separation
  unfold epistasis
  omega

-- Distance two can share a window: this architecture can have local
-- epistasis even though the distant interaction is structurally forced to zero.
def allA : DNA := fun _ => .A

def complement : Base → Base
  | .A => .T
  | .C => .G
  | .G => .C
  | .T => .A

def reverseComplement (trimer : Trimer) : Trimer :=
  ⟨complement trimer.third, complement trimer.second, complement trimer.first⟩

def cacWeight (trimer : Trimer) : Int :=
  if trimer = ⟨.C, .A, .C⟩ ∨ trimer = ⟨.G, .T, .G⟩ then 1 else 0

theorem cac_weight_is_reverse_complement_invariant (trimer : Trimer) :
    cacWeight (reverseComplement trimer) = cacWeight trimer := by
  cases trimer with
  | mk first second third =>
    cases first <;> cases second <;> cases third <;> decide

theorem distance_two_can_have_nonzero_epistasis :
    epistasis 3 0 cacWeight (gcWeights 0) allA 0 2 .C .C = 1 := by
  decide

end HelixVariantLocality
