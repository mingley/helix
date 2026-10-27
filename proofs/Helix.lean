/- The kernel checks implication and arithmetic, not empirical biological truth.
   Scientific premises and bridge rules remain explicit assumptions. No axioms,
   sorry, native_decide, external packages, or untrusted agent-generated tactics. -/
set_option autoImplicit false

namespace Helix

theorem implication_transitive (a b c : Prop)
    (ab : a → b) (bc : b → c) (ha : a) : c :=
  bc (ab ha)

def Contradicts (sameVariant sameContext oppositeConclusion : Prop) : Prop :=
  sameVariant ∧ sameContext ∧ oppositeConclusion

theorem contradiction_needs_context
    (v c o : Prop) (h : Contradicts v c o) : c :=
  h.2.1

theorem no_context_no_contradiction
    (v c o : Prop) (hc : ¬c) : ¬Contradicts v c o :=
  fun h => hc h.2.1

end Helix
