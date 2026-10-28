import CLPoly.Algorithm.Recombine
import Mathlib.Data.Nat.Choose.Vandermonde
import Mathlib.Data.Nat.Choose.Bounds

set_option autoImplicit false

open Polynomial

namespace CLPoly.VanHoeij

/-- The CLD coefficient bound before replacing the binomial coefficient by
a power of two.  This applies to nonmonic factors as well. -/
theorem norm_coeff_mul_derivative_le
    (g h : Polynomial ℂ) (d : Nat) (hg : g.natDegree = d + 1) (j : Nat) :
    ‖(h * derivative g).coeff j‖ ≤
      ((d + 1 : Nat) : ℝ) * ((h.natDegree + d).choose j : ℝ) *
        (h * g).mahlerMeasure := by
  have hderiv (k : Nat) : ‖(derivative g).coeff k‖ ≤
      (((d + 1) * d.choose k : Nat) : ℝ) * g.mahlerMeasure := by
    rw [coeff_derivative, norm_mul]
    have hcast : ‖((k : ℂ) + 1)‖ = ((k + 1 : Nat) : ℝ) := by
      norm_cast
    rw [hcast]
    calc
      ‖g.coeff (k + 1)‖ * ((k + 1 : Nat) : ℝ) ≤
          ((g.natDegree.choose (k + 1) : ℝ) * g.mahlerMeasure) *
            ((k + 1 : Nat) : ℝ) :=
        mul_le_mul_of_nonneg_right
          (norm_coeff_le_choose_mul_mahlerMeasure (k + 1) g) (by positivity)
      _ = ((((d + 1).choose (k + 1) * (k + 1) : Nat) : ℝ)) *
          g.mahlerMeasure := by rw [hg]; push_cast; ring
      _ = _ := by rw [← Nat.add_one_mul_choose_eq]
  have hv : (∑ ij ∈ Finset.antidiagonal j,
      (h.natDegree.choose ij.1 : ℝ) * (d.choose ij.2 : ℝ)) =
      ((h.natDegree + d).choose j : ℝ) := by
    exact_mod_cast (Nat.add_choose_eq h.natDegree d j).symm
  rw [coeff_mul]
  calc
    ‖∑ ij ∈ Finset.antidiagonal j, h.coeff ij.1 * (derivative g).coeff ij.2‖ ≤
        ∑ ij ∈ Finset.antidiagonal j,
          ‖h.coeff ij.1‖ * ‖(derivative g).coeff ij.2‖ := by
      simpa only [norm_mul] using norm_sum_le
        (Finset.antidiagonal j) (fun ij => h.coeff ij.1 * (derivative g).coeff ij.2)
    _ ≤ ∑ ij ∈ Finset.antidiagonal j,
        ((h.natDegree.choose ij.1 : ℝ) * h.mahlerMeasure) *
          (((d + 1) * d.choose ij.2 : Nat) : ℝ) * g.mahlerMeasure := by
      apply Finset.sum_le_sum
      intro ij _
      simpa only [mul_assoc] using mul_le_mul
        (norm_coeff_le_choose_mul_mahlerMeasure ij.1 h) (hderiv ij.2)
        (norm_nonneg _)
        (mul_nonneg (Nat.cast_nonneg _) (mahlerMeasure_nonneg h))
    _ = ∑ ij ∈ Finset.antidiagonal j,
        (((d + 1 : Nat) : ℝ) * (h.mahlerMeasure * g.mahlerMeasure)) *
          ((h.natDegree.choose ij.1 : ℝ) * (d.choose ij.2 : ℝ)) := by
      apply Finset.sum_congr rfl
      intro ij _
      push_cast
      ring
    _ = (((d + 1 : Nat) : ℝ) * (h.mahlerMeasure * g.mahlerMeasure)) *
        ((h.natDegree + d).choose j : ℝ) := by rw [← Finset.mul_sum, hv]
    _ = _ := by rw [mahlerMeasure_mul]; ring

/-- A uniform CLD bound expressed only in the original polynomial. -/
theorem norm_coeff_mul_derivative_le_l1
    (g h : Polynomial ℂ) (hg : g ≠ 0) (hh : h ≠ 0)
    (hdegree : 0 < g.natDegree) (j : Nat) :
    ‖(h * derivative g).coeff j‖ ≤
      ((h * g).natDegree : ℝ) * (2 : ℝ) ^ ((h * g).natDegree - 1) *
        (h * g).sum (fun _ coefficient => ‖coefficient‖) := by
  obtain ⟨d, hd⟩ := Nat.exists_eq_succ_of_ne_zero (Nat.ne_of_gt hdegree)
  have hd' : g.natDegree = d + 1 := hd
  have hn : (h * g).natDegree = h.natDegree + (d + 1) := by
    rw [natDegree_mul hh hg, hd']
  have hexponent : (h * g).natDegree - 1 = h.natDegree + d := by omega
  have hchoose : ((h.natDegree + d).choose j : ℝ) ≤
      (2 : ℝ) ^ ((h * g).natDegree - 1) := by
    rw [hexponent]
    exact_mod_cast Nat.choose_le_two_pow (h.natDegree + d) j
  have hdeg : ((d + 1 : Nat) : ℝ) ≤ ((h * g).natDegree : ℝ) := by
    exact_mod_cast (show d + 1 ≤ (h * g).natDegree by omega)
  have hnumeric := mul_le_mul hdeg hchoose (by positivity) (by positivity)
  calc
    ‖(h * derivative g).coeff j‖ ≤
        ((d + 1 : Nat) : ℝ) * ((h.natDegree + d).choose j : ℝ) *
          (h * g).mahlerMeasure := norm_coeff_mul_derivative_le g h d hd' j
    _ ≤ ((h * g).natDegree : ℝ) * (2 : ℝ) ^ ((h * g).natDegree - 1) *
          (h * g).mahlerMeasure :=
      mul_le_mul_of_nonneg_right hnumeric (mahlerMeasure_nonneg _)
    _ ≤ _ := mul_le_mul_of_nonneg_left (mahlerMeasure_le_sum_norm_coeff _)
      (by positivity)

/-- Integer form of exactly the global coefficient bound computed by
`__vanhoeij_cld_bound`.  It does not assume irreducibility of an output. -/
theorem integer_cld_coeff_bound
    (g h : Polynomial Int) (hg : g ≠ 0) (hh : h ≠ 0)
    (hdegree : 0 < g.natDegree) (j : Nat) :
    ((h * derivative g).coeff j).natAbs ≤
      (h * g).natDegree * 2 ^ ((h * g).natDegree - 1) *
        (Finset.range ((h * g).natDegree + 1)).sum
          (fun k => ((h * g).coeff k).natAbs) := by
  let cast : Int →+* ℂ := Int.castRingHom ℂ
  have hinjective : Function.Injective cast := Int.cast_injective
  have hbound := norm_coeff_mul_derivative_le_l1 (g.map cast) (h.map cast)
    ((Polynomial.map_ne_zero_iff hinjective).2 hg)
    ((Polynomial.map_ne_zero_iff hinjective).2 hh)
    (by simpa only [natDegree_map_eq_of_injective hinjective] using hdegree) j
  simp only [derivative_map, ← Polynomial.map_mul,
    natDegree_map_eq_of_injective hinjective] at hbound
  have hnorm : ‖((h * derivative g).map cast).coeff j‖ =
      (((h * derivative g).coeff j).natAbs : ℝ) := by
    rw [coeff_map]
    change ‖(((h * derivative g).coeff j : Int) : ℂ)‖ = _
    rw [Complex.norm_intCast, Nat.cast_natAbs, Int.cast_abs]
  have hsum : ((h * g).map cast).sum (fun _ coefficient => ‖coefficient‖) =
      (((Finset.range ((h * g).natDegree + 1)).sum
        (fun k => ((h * g).coeff k).natAbs) : Nat) : ℝ) := by
    rw [((h * g).map cast).sum_over_range (fun _ => by simp),
      natDegree_map_eq_of_injective hinjective]
    push_cast
    apply Finset.sum_congr rfl
    intro k _
    rw [coeff_map]
    change ‖(((h * g).coeff k : Int) : ℂ)‖ = _
    rw [Complex.norm_intCast, Nat.cast_natAbs, Int.cast_abs]
  rw [hnorm, hsum] at hbound
  exact_mod_cast hbound

end CLPoly.VanHoeij
