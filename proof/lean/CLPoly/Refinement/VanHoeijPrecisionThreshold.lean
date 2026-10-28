import CLPoly.Refinement.VanHoeijSeparationHeight
import Mathlib.Algebra.Order.Archimedean.Basic

set_option autoImplicit false

namespace Refinement.StrictRecombine
open Generated.StrictVanHoeijParameters

def inputDivisorLimit (height : Int) : Int :=
  divisorLoop height 1 (by decide) * divisorLoop height 1 (by decide)

def parameterSuccessThreshold (count : Int32) (height : Int) : Int :=
  max height (inputDivisorLimit height * (8 * count.toInt * (32 + 2 * count.toInt) + 2))

theorem inputDivisorLimit_positive (height : Int) : 0 < inputDivisorLimit height := by
  have h := (divisorLoop_positive_and_bound height 1 (by decide)).1
  exact mul_pos h h

/-- Parameter selection eventually succeeds at an input-fixed modulus bound,
even though the fast bit-size truncation policy still varies below its cap. -/
theorem selectParameters_complete_of_large_modulus
    (count : Int32) (height modulus : Int) (hc : 0 ≤ count.toInt)
    (_hh : 0 ≤ height) (hm : parameterSuccessThreshold count height < modulus) :
    ∃ parameters, selectParameters count height modulus = some parameters := by
  have hheight : height < modulus := (le_max_left _ _).trans_lt hm
  have hlarge : inputDivisorLimit height * (8 * count.toInt * (32 + 2 * count.toInt) + 2) < modulus :=
    (le_max_right _ _).trans_lt hm
  let D := truncationDivisor count height modulus
  have hd : 0 < D := truncationDivisor_positive count height modulus
  have hb : D ≤ inputDivisorLimit height := (truncationDivisor_bounds count height modulus).2
  have hu : 0 ≤ 8 * count.toInt * (32 + 2 * count.toInt) + 2 := by positivity
  have hdu : D * (8 * count.toInt * (32 + 2 * count.toInt) + 2) < modulus :=
    (mul_le_mul_of_nonneg_right hb hu).trans_lt hlarge
  have hquot : 8 * count.toInt * (32 + 2 * count.toInt) + 2 ≤ modulus / D := by
    have he := Int.mul_ediv_add_emod modulus D
    have hr := Int.emod_lt_of_pos modulus hd
    by_contra hn
    have hq : modulus / D ≤ 8 * count.toInt * (32 + 2 * count.toInt) + 1 := by omega
    nlinarith
  have hres : (height + D - 1) / D + 2 * count.toInt ≤ 32 + 2 * count.toInt := by
    have hhD := truncationDivisor_height_bound count height modulus
    have hceil := ceil_div_le_of_le_mul height D 32 hd hhD
    omega
  have hqpos : 1 < modulus / D := by
    have hterm : 0 ≤ 8 * count.toInt * (32 + 2 * count.toInt) := by positivity
    omega
  have hstrong : 8 * count.toInt * ((height + D - 1) / D + 2 * count.toInt) < modulus / D := by
    have he := mul_le_mul_of_nonneg_left hres (show 0 ≤ 8 * count.toInt by positivity)
    omega
  refine ⟨⟨D, modulus / D, (height + D - 1) / D + 2 * count.toInt⟩, ?_⟩
  unfold selectParameters
  rw [if_neg (not_le_of_gt hheight)]
  change (if modulus / D ≤ 1 then none else
    if modulus / D ≤ 8 * count.toInt * ((height + D - 1) / D + 2 * count.toInt) then none
    else some (⟨D, modulus / D, (height + D - 1) / D + 2 * count.toInt⟩ : Parameters)) = _
  rw [if_neg (not_le_of_gt hqpos), if_neg (not_le_of_gt hstrong)]

theorem separationPolynomialHeight_mono (source : Polynomial Int) (count : Nat)
    (left right : Int) (h : left ≤ right) :
    separationPolynomialHeight source count left ≤ separationPolynomialHeight source count right := by
  have hK := separationCoordinateBound_nonneg count source.natDegree
  have hres : separationResidualHeight count source.natDegree left ≤
      separationResidualHeight count source.natDegree right := by
    exact mul_le_mul_of_nonneg_right
      (mul_le_mul_of_nonneg_right h (by positivity)) hK
  unfold separationPolynomialHeight
  omega

theorem separationResultantThreshold_mono (source : Polynomial Int) (count : Nat)
    (left right : Int) (hl : 0 ≤ left) (h : left ≤ right) :
    separationResultantThreshold source count left ≤ separationResultantThreshold source count right := by
  have hb := separationPolynomialHeight_mono source count left right h
  have hn : 0 ≤ separationPolynomialHeight source count left :=
    le_trans (by decide) (separationPolynomialHeight_ge_one source count left hl)
  exact mul_le_mul_of_nonneg_left (pow_le_pow_left₀ hn hb _) (by positivity)

/-- This bound depends on the input height, count and source, never on the
current modulus or its chosen truncation divisor. -/
def vanHoeijInputThreshold (source : Polynomial Int) (count : Int32) (height : Int) : Int :=
  max (parameterSuccessThreshold count height)
    (max (separationResultantThreshold source count.toInt.toNat (inputDivisorLimit height))
      (2 * separationCoordinateBound count.toInt.toNat source.natDegree))

theorem selected_parameters_above_input_threshold
    (source : Polynomial Int) (count : Int32) (height modulus : Int)
    (hc : 0 ≤ count.toInt) (hh : 0 ≤ height)
    (hm : vanHoeijInputThreshold source count height < modulus) :
    ∃ parameters, selectParameters count height modulus = some parameters ∧
      separationResultantThreshold source count.toInt.toNat parameters.divisor < modulus ∧
      2 * separationCoordinateBound count.toInt.toNat source.natDegree < modulus := by
  have hp : parameterSuccessThreshold count height < modulus := (le_max_left _ _).trans_lt hm
  have hmax : max (separationResultantThreshold source count.toInt.toNat (inputDivisorLimit height))
      (2 * separationCoordinateBound count.toInt.toNat source.natDegree) < modulus :=
    (le_max_right _ _).trans_lt hm
  obtain ⟨parameters, hs⟩ := selectParameters_complete_of_large_modulus count height modulus hc hh hp
  obtain ⟨_, hd, he, _, _, _, _⟩ := selectParameters_spec count height modulus parameters hs
  have hb : parameters.divisor ≤ inputDivisorLimit height := by
    rw [he]; exact (truncationDivisor_bounds count height modulus).2
  refine ⟨parameters, hs, ?_, (le_max_right _ _).trans_lt hmax⟩
  exact (separationResultantThreshold_mono source count.toInt.toNat parameters.divisor
    (inputDivisorLimit height) hd.le hb).trans_lt ((le_max_left _ _).trans_lt hmax)

/-- Sufficient prime-power precision exists for the actual bounded-divisor
selection policy. This does not assert that a bounded machine retry loop
already reaches the required exponent. -/
theorem prime_power_eventually_sufficient_parameters
    (source : Polynomial Int) (count : Int32) (height : Int) (prime : Nat)
    (hp : 1 < prime) (hc : 0 ≤ count.toInt) (hh : 0 ≤ height) :
    ∃ exponent : Nat, 0 < exponent ∧ ∀ later, exponent ≤ later →
      ∃ parameters, selectParameters count height ((prime : Int) ^ later) = some parameters ∧
        separationResultantThreshold source count.toInt.toNat parameters.divisor < (prime : Int) ^ later ∧
        2 * separationCoordinateBound count.toInt.toNat source.natDegree < (prime : Int) ^ later := by
  have hp' : (1 : Int) < prime := by exact_mod_cast hp
  obtain ⟨exponent, he⟩ := pow_unbounded_of_one_lt
    (max (vanHoeijInputThreshold source count height) 1) hp'
  have hexponent : 0 < exponent := by
    by_contra h
    have hz : exponent = 0 := by omega
    simp only [hz, pow_zero] at he
    exact (not_lt_of_ge (le_max_right _ _)) he
  refine ⟨exponent, hexponent, ?_⟩
  intro later hl
  have hb : vanHoeijInputThreshold source count height < (prime : Int) ^ later :=
    ((le_max_left _ _).trans_lt he).trans_le (pow_le_pow_right₀ hp'.le hl)
  exact selected_parameters_above_input_threshold source count height _ hc hh hb

end Refinement.StrictRecombine
