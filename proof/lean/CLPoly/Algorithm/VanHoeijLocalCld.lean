import Mathlib.Algebra.Polynomial.Derivative
import Mathlib.Algebra.Polynomial.Div
import Mathlib.RingTheory.Coprime.Basic
import Mathlib.Tactic

set_option autoImplicit false

namespace CLPoly.VanHoeij

open Polynomial

/-- The diagonal CLD agrees with f' modulo its own monic factor, over a
general commutative ring (in particular over a prime-power residue ring). -/
theorem monic_factor_dvd_derivative_sub_cld {R : Type*} [CommRing R]
    (f q : Polynomial R) (hq : q.Monic) (hdiv : q ∣ f) :
    q ∣ derivative f - (f /ₘ q) * derivative q := by
  obtain ⟨g, rfl⟩ := hdiv
  rw [mul_divByMonic_cancel_left g hq, derivative_mul]
  refine ⟨derivative g, ?_⟩
  ring

/-- Every off-diagonal CLD is divisible by the selected local factor;
coprimality cancels the other monic divisor from its quotient. -/
theorem coprime_factor_dvd_cld {R : Type*} [CommRing R]
    (f q other : Polynomial R) (ho : other.Monic)
    (hqf : q ∣ f) (hof : other ∣ f) (hcop : IsCoprime q other) :
    q ∣ (f /ₘ other) * derivative other := by
  obtain ⟨g, hfg⟩ := hof
  rw [hfg, mul_divByMonic_cancel_left g ho]
  have hg : q ∣ g := hcop.dvd_of_dvd_mul_left (hfg ▸ hqf)
  exact dvd_mul_of_dvd_left hg _

/-- Local divisibility of an arbitrary weighted CLD sum. Only input
monic divisors and their pairwise coprimality are used. -/
theorem local_factor_dvd_weighted_cld_residual {R ι : Type*} [CommRing R]
    (support : Finset ι) (f : Polynomial R) (factors : ι → Polynomial R)
    (weights : ι → R) (selected : ι) (hselected : selected ∈ support)
    (hmonic : ∀ i ∈ support, (factors i).Monic)
    (hdiv : ∀ i ∈ support, factors i ∣ f)
    (hcop : ∀ i ∈ support, i ≠ selected → IsCoprime (factors selected) (factors i)) :
    factors selected ∣
      (∑ i ∈ support, C (weights i) * ((f /ₘ factors i) * derivative (factors i))) -
        C (weights selected) * derivative f := by
  classical
  have hterms : ∀ i ∈ support, factors selected ∣
      C (weights i) * ((f /ₘ factors i) * derivative (factors i)) -
        (if i = selected then C (weights selected) * derivative f else 0) := by
    intro i hi
    by_cases he : i = selected
    · subst i
      simp only [↓reduceIte]
      have hd := dvd_neg.mpr (monic_factor_dvd_derivative_sub_cld f (factors selected)
        (hmonic selected hselected) (hdiv selected hselected))
      convert dvd_mul_of_dvd_right hd (C (weights selected)) using 1 <;> ring
    · simp only [he, ↓reduceIte, sub_zero]
      exact dvd_mul_of_dvd_right (coprime_factor_dvd_cld f (factors selected) (factors i)
        (hmonic i hi) (hdiv selected hselected) (hdiv i hi) (hcop i hi he)) _
  have hs := Finset.dvd_sum hterms
  simpa only [Finset.sum_sub_distrib, Finset.sum_ite_eq', if_pos hselected] using hs

end CLPoly.VanHoeij
