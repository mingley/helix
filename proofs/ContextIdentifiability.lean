/- Exact identifiability statements for an IDEAL three-coordinate row-centering
   map. This does not assert which normalization produced a published workbook.
   Std-only: integer arithmetic represents three times the centered row, without
   truncated division. Rational finite panels admit a common-denominator encoding.
   No axioms, sorry, native_decide, source outcomes, or external packages. -/
import Std
set_option autoImplicit false

namespace HelixContextIdentifiability

structure Row where
  first : Int
  second : Int
  third : Int
  deriving DecidableEq

@[ext] theorem Row.ext {x y : Row}
    (first : x.first = y.first) (second : x.second = y.second)
    (third : x.third = y.third) : x = y := by
  cases x
  cases y
  simp_all

def zeroRow : Row := ⟨0, 0, 0⟩

def shift (x : Row) (common : Int) : Row :=
  ⟨x.first + common, x.second + common, x.third + common⟩

def rowSum (x : Row) : Int := x.first + x.second + x.third

-- This equals 3 * (x - mean(x) * 1) over rationals, not integer division by 3.
def center3 (x : Row) : Row :=
  ⟨3 * x.first - rowSum x,
   3 * x.second - rowSum x,
   3 * x.third - rowSum x⟩

theorem common_shift_invariant (x : Row) (common : Int) :
    center3 (shift x common) = center3 x := by
  apply Row.ext <;> simp only [center3, shift, rowSum] <;> omega

theorem centered_row_sum_zero (x : Row) : rowSum (center3 x) = 0 := by
  simp only [center3, rowSum]
  omega

theorem common_shift_changes_sum (x : Row) (common : Int) :
    rowSum (shift x common) = rowSum x + 3 * common := by
  simp only [rowSum, shift]
  omega

theorem common_shift_is_distinct (x : Row) (common : Int) (nonzero : common ≠ 0) :
    shift x common ≠ x := by
  intro equal
  have first := congrArg Row.first equal
  simp only [shift] at first
  omega

-- The full equivalence class, not only one direction of invariance.
theorem equal_specificity_iff_common_shift (x y : Row) :
    center3 y = center3 x ↔ ∃ common : Int, y = shift x common := by
  constructor
  · intro equal
    have first := congrArg Row.first equal
    have second := congrArg Row.second equal
    have third := congrArg Row.third equal
    simp only [center3, rowSum] at first second third
    refine ⟨y.first - x.first, ?_⟩
    apply Row.ext <;> simp only [shift] <;> omega
  · rintro ⟨common, rfl⟩
    exact common_shift_invariant x common

def fromTwoContrasts (firstVsThird secondVsThird : Int) : Row :=
  ⟨2 * firstVsThird - secondVsThird,
   2 * secondVsThird - firstVsThird,
   -firstVsThird - secondVsThird⟩

theorem determined_by_two_contrasts (x : Row) :
    center3 x = fromTwoContrasts (x.first - x.third) (x.second - x.third) := by
  apply Row.ext <;> simp only [center3, rowSum, fromTwoContrasts] <;> omega

theorem first_contrast_preserved_up_to_three (x : Row) :
    (center3 x).first - (center3 x).third = 3 * (x.first - x.third) := by
  simp only [center3, rowSum]
  omega

theorem second_contrast_preserved_up_to_three (x : Row) :
    (center3 x).second - (center3 x).third = 3 * (x.second - x.third) := by
  simp only [center3, rowSum]
  omega

-- The arbitrary per-element function need not be constant across the panel.
def centerPanel {Element : Type} (panel : Element → Row) : Element → Row :=
  fun element => center3 (panel element)

def shiftPanel {Element : Type} (panel : Element → Row)
    (commonActivity : Element → Int) : Element → Row :=
  fun element => shift (panel element) (commonActivity element)

theorem arbitrary_common_activity_vector_is_unobservable
    {Element : Type} (panel : Element → Row) (commonActivity : Element → Int) :
    centerPanel (shiftPanel panel commonActivity) = centerPanel panel := by
  funext element
  exact common_shift_invariant (panel element) (commonActivity element)

theorem arbitrary_common_activity_vector_changes_total
    {Element : Type} (panel : Element → Row) (commonActivity : Element → Int)
    (element : Element) :
    rowSum (shiftPanel panel commonActivity element) =
      rowSum (panel element) + 3 * commonActivity element := by
  exact common_shift_changes_sum (panel element) (commonActivity element)

-- This includes any fixed cell rescaling after centering, even an operation
-- whose output no longer has an unweighted zero sum. No rank computation or
-- inference about the actual source's operation order is assumed.
theorem downstream_postprocessing_cannot_restore_common_activity
    {Output : Type} (postprocess : Row → Output) (x : Row) (common : Int) :
    postprocess (center3 (shift x common)) = postprocess (center3 x) := by
  rw [common_shift_invariant]

def addRows (x delta : Row) : Row :=
  ⟨x.first + delta.first, x.second + delta.second, x.third + delta.third⟩

-- Division-cleared fixed affine normalization: normalized_c = gain_c*x_c+offset_c.
-- Rational fixed gains/offsets can be encoded by a common positive denominator.
def cellAffine (x gains offsets : Row) : Row :=
  ⟨gains.first * x.first + offsets.first,
   gains.second * x.second + offsets.second,
   gains.third * x.third + offsets.third⟩

-- An upstream raw shift is hidden when it produces the SAME normalized shift
-- in every cell. This is a condition to check, not a claim that unequal gains
-- preserve a literal common shift in raw units.
theorem fixed_cell_affine_hidden_direction
    (x delta gains offsets : Row) (normalizedCommon : Int)
    (first : gains.first * delta.first = normalizedCommon)
    (second : gains.second * delta.second = normalizedCommon)
    (third : gains.third * delta.third = normalizedCommon) :
    center3 (cellAffine (addRows x delta) gains offsets) =
      center3 (cellAffine x gains offsets) := by
  have transformed : cellAffine (addRows x delta) gains offsets =
      shift (cellAffine x gains offsets) normalizedCommon := by
    apply Row.ext
    · simp only [cellAffine, addRows, shift, Int.mul_add]
      rw [first]
      omega
    · simp only [cellAffine, addRows, shift, Int.mul_add]
      rw [second]
      omega
    · simp only [cellAffine, addRows, shift, Int.mul_add]
      rw [third]
      omega
  rw [transformed, common_shift_invariant]

theorem equal_cell_gains_preserve_common_raw_shift
    (x offsets : Row) (gain common : Int) :
    center3 (cellAffine (shift x common) ⟨gain, gain, gain⟩ offsets) =
      center3 (cellAffine x ⟨gain, gain, gain⟩ offsets) := by
  exact fixed_cell_affine_hidden_direction x ⟨common, common, common⟩
    ⟨gain, gain, gain⟩ offsets (gain * common) rfl rfl rfl

-- Unequal preprocessing gains can destroy common-raw-shift invariance.
-- x=0, gains=(1,2,3), offsets=0: a raw shift of 1 yields (-3,0,3) after center3.
theorem unequal_cell_gains_need_not_preserve_common_raw_shift :
    center3 (cellAffine (shift zeroRow 1) ⟨1, 2, 3⟩ zeroRow) ≠
      center3 (cellAffine zeroRow ⟨1, 2, 3⟩ zeroRow) := by
  decide

theorem zero_specificity_does_not_imply_zero_activity
    (common : Int) (nonzero : common ≠ 0) :
    center3 (shift zeroRow common) = zeroRow ∧ shift zeroRow common ≠ zeroRow := by
  constructor
  · rw [common_shift_invariant]
    decide
  · exact common_shift_is_distinct zeroRow common nonzero

theorem arbitrarily_large_common_activity_difference (bound : Nat) :
    ∃ x y : Row, center3 y = center3 x ∧ y.first - x.first > (bound : Int) := by
  refine ⟨zeroRow, shift zeroRow ((bound : Int) + 1), ?_, ?_⟩
  · exact common_shift_invariant zeroRow ((bound : Int) + 1)
  · simp only [shift, zeroRow]
    omega

theorem exact_three_cell_illustration :
    center3 ⟨2, 4, 3⟩ = ⟨-3, 3, 0⟩ ∧
      center3 ⟨7, 9, 8⟩ = ⟨-3, 3, 0⟩ := by
  decide

theorem no_exact_raw_activity_decoder :
    ¬ ∃ decode : Row → Row, ∀ x : Row, decode (center3 x) = x := by
  rintro ⟨decode, reconstruct⟩
  have atZero := reconstruct zeroRow
  have atShift := reconstruct (shift zeroRow 1)
  rw [common_shift_invariant] at atShift
  have impossible : shift zeroRow 1 = zeroRow := atShift.symm.trans atZero
  have first := congrArg Row.first impossible
  simp only [shift, zeroRow] at first
  omega

/- A stronger finite witness: both COMPLETE panels satisfy the same column
   first/second moments. Dividing either panel by the same positive real sample
   SD sqrt(14/3), or population SD sqrt(14/4), yields zero-mean/unit-variance
   columns. Those square-root/variance operations are described in the note;
   all underlying integer moment and specificity equalities are checked here. -/
def momentPanelX : List Row :=
  [⟨2, 1, 0⟩, ⟨0, 2, 1⟩, ⟨1, 0, 2⟩, ⟨-3, -3, -3⟩]

def momentPanelY : List Row :=
  [⟨0, -1, -2⟩, ⟨-2, 0, -1⟩, ⟨-1, -2, 0⟩, ⟨3, 3, 3⟩]

def squaredRow (x : Row) : Row :=
  ⟨x.first * x.first, x.second * x.second, x.third * x.third⟩

def columnSums (panel : List Row) : Row := panel.foldl addRows zeroRow

def columnSquareSums (panel : List Row) : Row := columnSums (panel.map squaredRow)

def fourElementMomentClass (panel : List Row) : Prop :=
  panel.length = 4 ∧ columnSums panel = zeroRow ∧
    columnSquareSums panel = ⟨14, 14, 14⟩

theorem moment_panel_x_valid : fourElementMomentClass momentPanelX := by
  unfold fourElementMomentClass
  decide

theorem moment_panel_y_valid : fourElementMomentClass momentPanelY := by
  unfold fourElementMomentClass
  decide

theorem moment_panels_differ_by_element_common_activity :
    momentPanelY = List.zipWith shift momentPanelX [-2, -2, -2, 6] := by decide

theorem moment_panels_have_distinct_common_activity :
    momentPanelX.map rowSum = [3, 3, 3, -9] ∧
      momentPanelY.map rowSum = [-3, -3, -3, 9] := by decide

theorem moment_panels_identical_specificity :
    momentPanelX.map center3 = momentPanelY.map center3 := by decide

theorem moment_panels_equal_after_common_specificity_postprocessing
    {Output : Type} (postprocess : List Row → Output) :
    postprocess (momentPanelX.map center3) = postprocess (momentPanelY.map center3) := by
  exact congrArg postprocess moment_panels_identical_specificity

theorem moment_panels_specificity_is_nonconstant :
    momentPanelX.map center3 =
      [⟨3, 0, -3⟩, ⟨-3, 3, 0⟩, ⟨0, -3, 3⟩, zeroRow] := by decide

theorem moment_panels_distinct : momentPanelX ≠ momentPanelY := by decide

theorem column_moment_constraints_do_not_identify_input :
    ∃ x y : List Row, fourElementMomentClass x ∧ fourElementMomentClass y ∧
      x.map center3 = y.map center3 ∧ x ≠ y := by
  exact ⟨momentPanelX, momentPanelY, moment_panel_x_valid, moment_panel_y_valid,
    moment_panels_identical_specificity, moment_panels_distinct⟩

theorem no_exact_decoder_even_with_column_moment_constraints :
    ¬ ∃ decode : List Row → List Row, ∀ panel : List Row,
      fourElementMomentClass panel → decode (panel.map center3) = panel := by
  rintro ⟨decode, reconstruct⟩
  have x := reconstruct momentPanelX moment_panel_x_valid
  have y := reconstruct momentPanelY moment_panel_y_valid
  rw [← moment_panels_identical_specificity] at y
  exact moment_panels_distinct (x.symm.trans y)

end HelixContextIdentifiability
