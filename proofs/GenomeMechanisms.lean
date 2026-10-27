/- Definitions and finite arithmetic for the HBB reference research run.
   These theorems do not verify the Rust implementation or biological premises. -/
import Std
set_option autoImplicit false

namespace HelixMechanisms

inductive Base where
  | A | C | G | T
  deriving DecidableEq

def complement : Base → Base
  | .A => .T
  | .T => .A
  | .C => .G
  | .G => .C

theorem complement_involution (base : Base) :
    complement (complement base) = base := by
  cases base <;> rfl

def reverseComplement (sequence : List Base) : List Base :=
  sequence.reverse.map complement

theorem reverse_complement_involution (sequence : List Base) :
    reverseComplement (reverseComplement sequence) = sequence := by
  simp [reverseComplement, List.map_reverse, List.map_map, Function.comp_def,
    complement_involution]

theorem one_based_roundtrip (offset : Nat) :
    (offset + 1) - 1 = offset := by omega

theorem minus_strand_offset_in_range (position endOffset : Nat)
    (inside : position < endOffset) :
    endOffset - 1 - position < endOffset := by omega

theorem category_only_cannot_fit_opposed_labels (predict : Nat → Bool)
    (category : Nat) :
    ¬ (predict category = true ∧ predict category = false) := by
  intro h
  have impossible : true = false := h.1.symm.trans h.2
  cases impossible

-- HBB NM_000518.5 reference quantities; provenance is checked by Rust, while
-- Lean checks only the displayed arithmetic.
theorem hbb_spliced_length : (263 + 223 + 142 : Nat) = 628 := by decide
theorem hbb_coding_length : (444 : Nat) = 3 * (147 + 1) := by decide
theorem hbb_locus_saturation : (1608 * 3 : Nat) = 4824 := by decide

end HelixMechanisms
