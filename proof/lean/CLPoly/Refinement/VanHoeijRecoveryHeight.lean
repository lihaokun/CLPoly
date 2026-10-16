import CLPoly.Refinement.VanHoeijSeparationHeight

set_option autoImplicit false

namespace Refinement.StrictRecombine

/-- The separation threshold also covers nonmonic recovery's scaled target
at every stage whose factor and complement divide the original input. -/
theorem separation_threshold_covers_scaled_recovery
    (source factor complement : Polynomial Int) (count : Nat) (divisor : Int)
    (hf : source ≠ 0) (hn : 0 < source.natDegree) (hd : 0 ≤ divisor)
    (hfactor : factor ∣ source) (hcomplement : complement ∣ source)
    (modulus : Nat) (hm : separationResultantThreshold source count divisor < (modulus : Int)) :
    ∀ degree, ((Polynomial.C complement.leadingCoeff * factor).coeff degree).natAbs * 2 < modulus := by
  let B := separationPolynomialHeight source count divisor
  have hB : 1 ≤ B := separationPolynomialHeight_ge_one source count divisor hd
  have hthreshold : 2 * B ^ 2 ≤ separationResultantThreshold source count divisor := by
    have h := resultant_coefficient_bound_mono 2 (2 * source.natDegree) B (by omega) hB
    simpa [separationResultantThreshold, B] using h
  intro degree
  have hl : |complement.leadingCoeff| ≤ B :=
    source_factor_coeff_le_separationPolynomialHeight source complement count divisor hd hf hcomplement
      complement.natDegree
  have hc : |factor.coeff degree| ≤ B :=
    source_factor_coeff_le_separationPolynomialHeight source factor count divisor hd hf hfactor degree
  have hprod : |(Polynomial.C complement.leadingCoeff * factor).coeff degree| ≤ B ^ 2 := by
    rw [Polynomial.coeff_C_mul, abs_mul, pow_two]
    exact mul_le_mul hl hc (abs_nonneg _) (by omega)
  have hlt : 2 * |(Polynomial.C complement.leadingCoeff * factor).coeff degree| < (modulus : Int) :=
    (mul_le_mul_of_nonneg_left hprod (by norm_num)).trans_lt (hthreshold.trans_lt hm)
  rw [Int.abs_eq_natAbs] at hlt
  have hlt' : (((Polynomial.C complement.leadingCoeff * factor).coeff degree).natAbs : Int) * 2 <
      (modulus : Int) := by linarith
  exact_mod_cast hlt'

end Refinement.StrictRecombine
