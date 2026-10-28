import CLPoly.Generated.StrictVanHoeijParameters
import Mathlib.Data.Int.Order.Basic
import Mathlib.Tactic

set_option autoImplicit false

namespace Refinement.StrictRecombine
open Generated.StrictVanHoeijParameters

theorem divisorLoop_positive_and_bound (height divisor : ZZ) (positive : 0 < divisor) :
    0 < divisorLoop height divisor positive ∧
      height ≤ divisorLoop height divisor positive * 32 := by
  fun_induction divisorLoop height divisor positive with
  | case1 divisor positive hlt ih =>
    exact ih
  | case2 divisor positive hlt =>
    dsimp only [ZZ] at *
    exact ⟨positive, by omega⟩

/-- The bit-size policy keeps an input-fixed lower and upper divisor bound. -/
theorem truncationDivisor_bounds (count : Int32) (height modulus : ZZ) :
    divisorLoop height 1 (by decide) ≤ truncationDivisor count height modulus ∧
      truncationDivisor count height modulus ≤
        divisorLoop height 1 (by decide) * divisorLoop height 1 (by decide) := by
  dsimp only [ZZ] at *
  have hp := (divisorLoop_positive_and_bound height 1 (by decide)).1
  have hs : divisorLoop height 1 (by decide) ≤
      divisorLoop height 1 (by decide) * divisorLoop height 1 (by decide) := by
    dsimp only [ZZ] at hp
    nlinarith
  unfold truncationDivisor
  dsimp only
  split
  · exact ⟨le_min hs (le_max_left _ _), min_le_left _ _⟩
  · exact ⟨le_refl _, hs⟩

theorem truncationDivisor_positive (count : Int32) (height modulus : ZZ) :
    0 < truncationDivisor count height modulus :=
  lt_of_lt_of_le (divisorLoop_positive_and_bound height 1 (by decide)).1
    (truncationDivisor_bounds count height modulus).1

/-- The bounded bit-size policy preserves the base height inequality. -/
theorem truncationDivisor_height_bound (count : Int32) (height modulus : ZZ) :
    height ≤ truncationDivisor count height modulus * 32 := by
  dsimp only [ZZ] at *
  have hb := (divisorLoop_positive_and_bound height 1 (by decide)).2
  exact hb.trans (mul_le_mul_of_nonneg_right
    (truncationDivisor_bounds count height modulus).1 (by decide))

theorem ceil_div_le_of_le_mul (height divisor ceiling : Int)
    (hd : 0 < divisor) (hh : height ≤ divisor * ceiling) :
    (height + divisor - 1) / divisor ≤ ceiling := by
  have heq := Int.mul_ediv_add_emod (height + divisor - 1) divisor
  have hr := Int.emod_nonneg (height + divisor - 1) (ne_of_gt hd)
  by_contra h
  have hq : ceiling + 1 ≤ (height + divisor - 1) / divisor := by omega
  nlinarith

theorem selectParameters_spec (count : Int32) (height originalModulus : ZZ)
    (parameters : Parameters)
    (hrun : selectParameters count height originalModulus = some parameters) :
    height < originalModulus ∧
      0 < parameters.divisor ∧
      parameters.divisor = truncationDivisor count height originalModulus ∧
      parameters.modulus = originalModulus / parameters.divisor ∧
      parameters.residual = (height + parameters.divisor - 1) / parameters.divisor +
        2 * count.toInt ∧
      1 < parameters.modulus ∧
      8 * count.toInt * parameters.residual < parameters.modulus := by
  dsimp only [ZZ] at *
  unfold selectParameters at hrun
  split at hrun
  next hm => contradiction
  next hm =>
    dsimp only at hrun
    split at hrun
    next hsmall => contradiction
    next hsmall =>
      split at hrun
      next hweak => contradiction
      next hweak =>
        have he := Option.some.inj hrun
        rw [← he]
        dsimp only [ZZ] at *
        exact ⟨by omega, truncationDivisor_positive count height originalModulus,
          rfl, rfl, rfl, by omega, by omega⟩

/-- The actual selected residual cutoff is uniformly bounded in the
factor count, despite the modulus-dependent bit-size truncation policy. -/
theorem selectParameters_residual_bounds (count : Int32) (height originalModulus : ZZ)
    (parameters : Parameters) (hcount : 0 ≤ count.toInt) (hheight : 0 ≤ height)
    (hrun : selectParameters count height originalModulus = some parameters) :
    0 ≤ parameters.residual ∧ parameters.residual ≤ 32 + 2 * count.toInt := by
  dsimp only [ZZ] at *
  obtain ⟨_, hd, heq, _, hresidual, _, _⟩ :=
    selectParameters_spec count height originalModulus parameters hrun
  have hb := truncationDivisor_height_bound count height originalModulus
  rw [← heq] at hb
  dsimp only [ZZ] at hd hresidual hb
  have hceil := ceil_div_le_of_le_mul height parameters.divisor 32 hd (by simpa [mul_comm] using hb)
  have hnonneg : (0 : Int) ≤ (height + parameters.divisor - 1) / parameters.divisor :=
    Int.ediv_nonneg (by omega) (le_of_lt hd)
  rw [hresidual]
  constructor <;> omega

/-- With a bound on the number of added columns, the source's squared
short-vector cutoff is uniformly bounded by input parameters. -/
theorem selectParameters_squared_cutoff_le
    (count : Int32) (height originalModulus : ZZ) (parameters : Parameters)
    (hcount : 0 ≤ count.toInt) (hheight : 0 ≤ height)
    (hrun : selectParameters count height originalModulus = some parameters)
    (added degree : Nat) (hadded : added ≤ degree) :
    count.toInt * 32 * 32 + (added : Int) * parameters.residual * parameters.residual ≤
      count.toInt * 32 * 32 + (degree : Int) * (32 + 2 * count.toInt) ^ 2 := by
  dsimp only [ZZ] at *
  obtain ⟨hr0, hr⟩ := selectParameters_residual_bounds count height originalModulus
    parameters hcount hheight hrun
  have hs : parameters.residual ^ 2 ≤ (32 + 2 * count.toInt) ^ 2 := by nlinarith
  have hc : (added : Int) ≤ degree := by exact_mod_cast hadded
  have hb := mul_le_mul hc hs (sq_nonneg _) (by positivity : (0 : Int) ≤ degree)
  nlinarith only [hb]

end Refinement.StrictRecombine
