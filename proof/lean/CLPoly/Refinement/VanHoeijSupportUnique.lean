import CLPoly.Refinement.VanHoeijRecoveryHeight
import CLPoly.Refinement.VanHoeijGroupCoverage

set_option autoImplicit false

namespace Refinement.StrictRecombine

/-- At the same explicit separation precision, an exact modular support
determines at most one positive-leading primitive divisor of the source. -/
theorem positive_primitive_divisor_unique_of_support
    (source left right : Polynomial Int) (count : Nat) (divisor : Int)
    (hf : source ≠ 0) (hn : 0 < source.natDegree) (hd : 0 ≤ divisor)
    (hl : left ∣ source) (hr : right ∣ source)
    (hlp : left.IsPrimitive) (hrp : right.IsPrimitive)
    (hlc : 0 < left.leadingCoeff) (hrc : 0 < right.leadingCoeff)
    (modulus : Nat) (hm : 0 < modulus)
    (hprecision : separationResultantThreshold source count divisor < (modulus : Int))
    (support : Polynomial (ZMod modulus))
    (hleft : left.map (Int.castRingHom (ZMod modulus)) =
      Polynomial.C (left.leadingCoeff : ZMod modulus) * support)
    (hright : right.map (Int.castRingHom (ZMod modulus)) =
      Polynomial.C (right.leadingCoeff : ZMod modulus) * support) : left = right := by
  have hcross : Polynomial.C right.leadingCoeff * left =
      Polynomial.C left.leadingCoeff * right := by
    apply symmetric_recovery_closed_left _ _ modulus hm
    · simp only [Polynomial.map_mul, Polynomial.map_C, hleft, hright]
      change Polynomial.C (right.leadingCoeff : ZMod modulus) *
        (Polynomial.C (left.leadingCoeff : ZMod modulus) * support) =
        Polynomial.C (left.leadingCoeff : ZMod modulus) *
        (Polynomial.C (right.leadingCoeff : ZMod modulus) * support)
      ring
    · intro degree
      exact (separation_threshold_covers_scaled_recovery source left right count divisor
        hf hn hd hl hr modulus hprecision degree).le
    · exact separation_threshold_covers_scaled_recovery source right left count divisor
        hf hn hd hr hl modulus hprecision
  have hcoeff := congrArg Polynomial.content hcross
  rw [Polynomial.content_C_mul, Polynomial.content_C_mul,
    hlp.content_eq_one, hrp.content_eq_one, mul_one, mul_one,
    Int.normalize_of_nonneg hlc.le, Int.normalize_of_nonneg hrc.le] at hcoeff
  rw [hcoeff] at hcross
  exact mul_left_cancel₀ (Polynomial.C_ne_zero.mpr (ne_of_gt hlc)) hcross

end Refinement.StrictRecombine
