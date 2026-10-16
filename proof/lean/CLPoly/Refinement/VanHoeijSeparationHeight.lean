import CLPoly.Refinement.VanHoeijUniformBounds
import CLPoly.Algorithm.VanHoeijSeparation

set_option autoImplicit false

namespace Refinement.StrictRecombine

def sourceL1 (f : Polynomial Int) : Nat :=
  ∑ k ∈ Finset.range (f.natDegree + 1), (f.coeff k).natAbs

def sourceFactorHeight (f : Polynomial Int) : Int :=
  ((Nat.choose f.natDegree (f.natDegree / 2) * sourceL1 f : Nat) : Int)

def separationPolynomialHeight (f : Polynomial Int) (count : Nat) (divisor : Int) : Int :=
  1 + sourceFactorHeight f + separationResidualHeight count f.natDegree divisor +
    separationCoordinateBound count f.natDegree * (f.natDegree : Int) * (sourceL1 f : Int)

def separationResultantThreshold (f : Polynomial Int) (count : Nat) (divisor : Int) : Int :=
  (((2 * f.natDegree).factorial : Nat) : Int) *
    separationPolynomialHeight f count divisor ^ (2 * f.natDegree)

theorem source_coeff_natAbs_le (f : Polynomial Int) (k : Nat) :
    (f.coeff k).natAbs ≤ sourceL1 f := by
  by_cases hk : k < f.natDegree + 1
  · exact Finset.single_le_sum (fun i _ => Nat.zero_le ((f.coeff i).natAbs)) (Finset.mem_range.mpr hk)
  · rw [Polynomial.coeff_eq_zero_of_natDegree_lt (by omega)]
    simp

theorem source_derivative_coeff_abs_le (f : Polynomial Int) (k : Nat) :
    |(Polynomial.derivative f).coeff k| ≤ (f.natDegree : Int) * (sourceL1 f : Int) := by
  rw [Polynomial.coeff_derivative, abs_mul]
  by_cases hk : k + 1 ≤ f.natDegree
  · have hc : |f.coeff (k + 1)| ≤ (sourceL1 f : Int) := by
      rw [Int.abs_eq_natAbs]
      exact_mod_cast source_coeff_natAbs_le f (k + 1)
    have hn : |((k : Int) + 1)| ≤ (f.natDegree : Int) := by
      rw [abs_of_nonneg (by positivity)]
      exact_mod_cast hk
    exact (mul_le_mul hc hn (abs_nonneg _) (by positivity)).trans_eq (by ring)
  · rw [Polynomial.coeff_eq_zero_of_natDegree_lt (by omega)]
    simp only [abs_zero, zero_mul]
    positivity

theorem source_factor_coeff_abs_le (f factor : Polynomial Int) (hf : f ≠ 0)
    (hdiv : factor ∣ f) (k : Nat) : |factor.coeff k| ≤ sourceFactorHeight f := by
  rw [Int.abs_eq_natAbs]
  unfold sourceFactorHeight sourceL1
  exact_mod_cast mignotte_bound f factor hf hdiv k

theorem separationPolynomialHeight_ge_one (f : Polynomial Int) (count : Nat)
    (divisor : Int) (hd : 0 ≤ divisor) : 1 ≤ separationPolynomialHeight f count divisor := by
  have hK := separationCoordinateBound_nonneg count f.natDegree
  have hS : 0 ≤ separationResidualHeight count f.natDegree divisor := by
    unfold separationResidualHeight
    positivity
  have hA : 0 ≤ sourceFactorHeight f := by unfold sourceFactorHeight; positivity
  have hD : 0 ≤ separationCoordinateBound count f.natDegree * (f.natDegree : Int) * (sourceL1 f : Int) := by
    positivity
  unfold separationPolynomialHeight
  omega

theorem residual_difference_coeff_bound (f H : Polynomial Int) (count : Nat)
    (divisor scalar : Int) (_hd : 0 ≤ divisor)
    (he : |scalar| ≤ separationCoordinateBound count f.natDegree)
    (hH : ∀ k, |H.coeff k| ≤ separationResidualHeight count f.natDegree divisor)
    (k : Nat) :
    |(H - Polynomial.C scalar * Polynomial.derivative f).coeff k| ≤
      separationPolynomialHeight f count divisor := by
  rw [Polynomial.coeff_sub, Polynomial.coeff_C_mul]
  have ht : |scalar * (Polynomial.derivative f).coeff k| ≤
      separationCoordinateBound count f.natDegree * (f.natDegree : Int) * (sourceL1 f : Int) := by
    rw [abs_mul]
    simpa only [mul_assoc] using mul_le_mul he (source_derivative_coeff_abs_le f k)
      (abs_nonneg _) (separationCoordinateBound_nonneg count f.natDegree)
  have hb := (abs_sub (H.coeff k) (scalar * (Polynomial.derivative f).coeff k)).trans
    (add_le_add (hH k) ht)
  have hA : 0 ≤ sourceFactorHeight f := by unfold sourceFactorHeight; positivity
  unfold separationPolynomialHeight
  omega

theorem source_factor_coeff_le_separationPolynomialHeight
    (f factor : Polynomial Int) (count : Nat) (divisor : Int) (hd : 0 ≤ divisor)
    (hf : f ≠ 0) (hdiv : factor ∣ f) (k : Nat) :
    |factor.coeff k| ≤ separationPolynomialHeight f count divisor := by
  have hb := source_factor_coeff_abs_le f factor hf hdiv k
  have hK := separationCoordinateBound_nonneg count f.natDegree
  have hS : 0 ≤ separationResidualHeight count f.natDegree divisor := by
    unfold separationResidualHeight
    positivity
  have hD : 0 ≤ separationCoordinateBound count f.natDegree * (f.natDegree : Int) * (sourceL1 f : Int) := by
    positivity
  unfold separationPolynomialHeight
  omega

theorem resultant_coefficient_bound_mono (degree limit : Nat) (B : Int)
    (hd : degree ≤ limit) (hB : 1 ≤ B) :
    (degree.factorial : Int) * B ^ degree ≤ (limit.factorial : Int) * B ^ limit := by
  have hf : (degree.factorial : Int) ≤ limit.factorial := by exact_mod_cast Nat.factorial_le hd
  exact mul_le_mul hf (pow_le_pow_right₀ hB hd) (by positivity) (by positivity)

end Refinement.StrictRecombine
