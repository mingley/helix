/- Small method contracts only: no assertion about biological truth, statistical
   calibration, fitted coefficients, or correctness of the Rust/Python programs. -/
import Std
set_option autoImplicit false
namespace HelixRegulatory

theorem guarded_supports_disjoint (trainCenter testCenter : Int)
    (before : trainCenter + 20 < testCenter) :
    trainCenter + 10 < testCenter - 10 := by omega

theorem guarded_supports_disjoint_after (trainCenter testCenter : Int)
    (after : testCenter + 20 < trainCenter) :
    testCenter + 10 < trainCenter - 10 := by omega

theorem identical_input_has_one_prediction
    {Input Output : Type} (predict : Input → Output) (x : Input)
    (a b : Output) (ha : predict x = a) (hb : predict x = b) : a = b :=
  ha.symm.trans hb

-- With 99 controls, a plus-one tail is >=1/100. Its six-comparison Bonferroni
-- value cannot be below 1/20. This is arithmetic, not a calibrated p-value claim.
theorem control_resolution (exceedances : Nat) :
    (100 : Nat) < 20 * (6 * (1 + exceedances)) := by omega

end HelixRegulatory

namespace Helix.Generated
theorem conditional_conclusion (p0 : Prop) (p1 : Prop) (h0 : p0) (r0 : p0 → p1) : p1 := (r0 h0)

theorem probability_0 : (15563 : Nat) ≤ 15699 ∧ (0 : Nat) < 15699 := by decide

end Helix.Generated
