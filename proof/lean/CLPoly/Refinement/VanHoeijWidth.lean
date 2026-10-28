import CLPoly.Generated.StrictVanHoeijParameters
import Mathlib.Data.Int.Order.Basic

set_option autoImplicit false

namespace Refinement.StrictRecombine
open Generated.StrictVanHoeijParameters
open Generated.StrictRecombine

theorem int32_max_toInt (left right : Int32) :
    (max left right).toInt = max left.toInt right.toInt := by
  change (if left ≤ right then right else left).toInt = _
  split
  · rename_i h
    have hh := Int32.le_iff_toInt_le.mp h
    omega
  · rename_i h
    have hh : ¬ left.toInt ≤ right.toInt := by
      simpa only [Int32.le_iff_toInt_le] using h
    omega

theorem int32_min_toInt (left right : Int32) :
    (min left right).toInt = min left.toInt right.toInt := by
  change (if left ≤ right then left else right).toInt = _
  split
  · rename_i h
    have hh := Int32.le_iff_toInt_le.mp h
    omega
  · rename_i h
    have hh : ¬ left.toInt ≤ right.toInt := by
      simpa only [Int32.le_iff_toInt_le] using h
    omega

theorem cld_poly_width_toInt (poly : SparsePolyZZ) (hnonempty : poly.isEmpty = false)
    (hfits : poly[0]!.1.deg + 1 < 2 ^ 31) :
    ((get_deg poly).toInt32 + 1).toInt = ((poly[0]!.1.deg + 1 : Nat) : Int) := by
  unfold get_deg
  rw [hnonempty]
  simp only [Bool.false_eq_true, ↓reduceIte]
  change ((Int64.ofNat poly[0]!.1.deg).toInt32 + Int32.ofNat 1).toInt = _
  rw [Int64.toInt32_ofNat', ← Int32.ofNat_add]
  exact Int32.toInt_ofNat_of_lt hfits

theorem cldWidthLoop_eq_natural (cld : Array SparsePolyZZ)
    (hfits : ∀ i (hi : i < cld.size), cld[i].isEmpty = false →
      cld[i][0]!.1.deg + 1 < 2 ^ 31)
    (index maximum : Nat) (width : Int32)
    (hwidth : width.toInt = (maximum : Int)) (hmaximum : maximum < 2 ^ 31) :
    (cldWidthLoop cld index width).toInt = (cldSpiralWidthLoop cld index maximum : Int) ∧
      cldSpiralWidthLoop cld index maximum < 2 ^ 31 := by
  induction hmeasure : cld.size - index using Nat.strong_induction_on
      generalizing index maximum width with
  | h measure ih =>
    rw [cldWidthLoop, cldSpiralWidthLoop]
    by_cases hi : index < cld.size
    · simp only [dif_pos hi]
      cases he : cld[index].isEmpty with
      | true =>
        simp only [he, ↓reduceIte, Nat.max_zero]
        exact ih _ (by omega) (index + 1) maximum width hwidth hmaximum rfl
      | false =>
        simp only [he, Bool.false_eq_true, ↓reduceIte]
        have hp := hfits index hi he
        apply ih (cld.size - (index + 1)) (by omega) (index + 1)
          (max maximum (cld[index][0]!.1.deg + 1))
          (max width ((get_deg cld[index]).toInt32 + 1))
        · rw [int32_max_toInt, hwidth, cld_poly_width_toInt cld[index] he hp]
          omega
        · omega
        · rfl
    · simp only [dif_neg hi]
      exact ⟨hwidth, hmaximum⟩

theorem cldWidth_eq_natural (cld : Array SparsePolyZZ)
    (hfits : ∀ i (hi : i < cld.size), cld[i].isEmpty = false →
      cld[i][0]!.1.deg + 1 < 2 ^ 31) :
    (cldWidth cld).toInt = (cldSpiralWidth cld : Int) ∧
      cldSpiralWidth cld < 2 ^ 31 := by
  exact cldWidthLoop_eq_natural cld hfits 0 0 0 (by decide) (by decide)

/-- The source initial target satisfies the column loop's erased range premise,
even before a nonnegative-width invariant is used. -/
theorem initialCldTarget_bound (width : Int32) :
    (min (2 : Int32) width).toInt ≤ width.toInt := by
  rw [int32_min_toInt]
  omega

theorem cldWidth_initialTarget_nonnegative (cld : Array SparsePolyZZ)
    (hfits : ∀ i (hi : i < cld.size), cld[i].isEmpty = false →
      cld[i][0]!.1.deg + 1 < 2 ^ 31) :
    0 ≤ (min (2 : Int32) (cldWidth cld)).toInt := by
  rw [int32_min_toInt, (cldWidth_eq_natural cld hfits).1]
  change 0 ≤ min (2 : Int) (cldSpiralWidth cld : Int)
  omega

end Refinement.StrictRecombine
