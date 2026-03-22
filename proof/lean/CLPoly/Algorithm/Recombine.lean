/-
  CLPoly/Algorithm/Recombine.lean — L2 因子重组正确性证明

  Phase 4: Recombination (Zassenhaus)
  对应 C++: polynomial_factorize_univar.hh:750-882 __zassenhaus_recombine
-/
import CLPoly.Spec
import CLPoly.Algorithm.Hensel
import Mathlib.Data.ZMod.Basic
import Mathlib.RingTheory.Polynomial.Basic
import Mathlib.Analysis.Polynomial.MahlerMeasure
import Mathlib.NumberTheory.MahlerMeasure

set_option autoImplicit false
set_option maxHeartbeats 3200000

open Polynomial

-- ============================================================
-- 1. Mignotte bound
-- ============================================================

/-- Helper: M(h) ≥ 1 for nonzero integer polynomial h -/
private lemma mahlerMeasure_ge_one_of_int_ne_zero (h : Polynomial ℤ) (hh : h ≠ 0) :
    1 ≤ (h.map (Int.castRingHom ℂ)).mahlerMeasure :=
  one_le_mahlerMeasure_of_ne_zero hh

/-- Mignotte bound: g | f → |g.coeff i| ≤ C(n, n/2) · ‖f‖₁ -/
theorem mignotte_bound (f g : Polynomial ℤ) (hf : f ≠ 0) (hg : g ∣ f) :
    ∀ i, (g.coeff i).natAbs ≤
      Nat.choose f.natDegree (f.natDegree / 2) *
      (Finset.range (f.natDegree + 1)).sum (fun j => (f.coeff j).natAbs) := by
  intro i
  have hφ_inj : Function.Injective (Int.castRingHom ℂ) := Int.cast_injective
  obtain ⟨h, hfgh⟩ := hg
  have hh_ne : h ≠ 0 := right_ne_zero_of_mul (hfgh ▸ hf)
  have hg_ne : g ≠ 0 := left_ne_zero_of_mul (hfgh ▸ hf)
  -- Work in ℝ and cast at the end
  suffices h_real : (↑((g.coeff i).natAbs) : ℝ) ≤
      ↑(Nat.choose f.natDegree (f.natDegree / 2) *
        (Finset.range (f.natDegree + 1)).sum (fun j => (f.coeff j).natAbs)) by
    exact_mod_cast h_real
  -- ‖g_c.coeff i‖ = (g.coeff i).natAbs (as ℝ)
  have h_lhs : ‖(g.map (Int.castRingHom ℂ)).coeff i‖ = ↑((g.coeff i).natAbs) := by
    rw [coeff_map, eq_intCast, Complex.norm_intCast]
    rw [Nat.cast_natAbs (α := ℝ), Int.cast_abs]
  rw [← h_lhs]
  -- ‖g_c.coeff i‖ ≤ C(d,i) · M(g_c) (Mathlib)
  have h1 := norm_coeff_le_choose_mul_mahlerMeasure i (g.map (Int.castRingHom ℂ))
  rw [natDegree_map_eq_of_injective hφ_inj] at h1
  -- M(g_c) ≤ M(f_c) (multiplicativity + M(h) ≥ 1)
  have h2 : (g.map (Int.castRingHom ℂ)).mahlerMeasure ≤
      (f.map (Int.castRingHom ℂ)).mahlerMeasure := by
    have hfgh_c : f.map (Int.castRingHom ℂ) =
        g.map (Int.castRingHom ℂ) * h.map (Int.castRingHom ℂ) := by
      rw [← Polynomial.map_mul, hfgh]
    rw [hfgh_c, mahlerMeasure_mul]
    exact le_mul_of_one_le_right (mahlerMeasure_nonneg _)
      (mahlerMeasure_ge_one_of_int_ne_zero h hh_ne)
  -- C(d,i) ≤ C(n, n/2)
  have hg_dvd : g ∣ f := ⟨h, hfgh⟩
  have h3 : g.natDegree.choose i ≤ f.natDegree.choose (f.natDegree / 2) :=
    (Nat.choose_le_middle i g.natDegree).trans
      ((Nat.choose_le_choose _ (natDegree_le_of_dvd hg_dvd hf)).trans
        (Nat.choose_le_middle _ f.natDegree))
  -- M(f_c) ≤ L1(f_c)
  have h4 := mahlerMeasure_le_sum_norm_coeff (f.map (Int.castRingHom ℂ))
  -- L1(f_c) = Σ (f.coeff j).natAbs (as ℝ)
  have h5 : (f.map (Int.castRingHom ℂ)).sum (fun _ a => ‖a‖) =
      ↑((Finset.range (f.natDegree + 1)).sum (fun j => (f.coeff j).natAbs)) := by
    rw [(f.map (Int.castRingHom ℂ)).sum_over_range (fun _ => by simp),
        natDegree_map_eq_of_injective hφ_inj]
    push_cast; congr 1; ext j
    rw [coeff_map, eq_intCast, Complex.norm_intCast]
    rw [Nat.cast_natAbs (α := ℝ), Int.cast_abs]
  -- Combine
  calc ‖(g.map (Int.castRingHom ℂ)).coeff i‖
    ≤ ↑(g.natDegree.choose i) * (g.map (Int.castRingHom ℂ)).mahlerMeasure := h1
    _ ≤ ↑(g.natDegree.choose i) * (f.map (Int.castRingHom ℂ)).mahlerMeasure := by
        gcongr
    _ ≤ ↑(f.natDegree.choose (f.natDegree / 2)) *
          (f.map (Int.castRingHom ℂ)).mahlerMeasure := by
        apply mul_le_mul_of_nonneg_right (by exact_mod_cast h3)
          (mahlerMeasure_nonneg _)
    _ ≤ ↑(f.natDegree.choose (f.natDegree / 2)) *
          ((f.map (Int.castRingHom ℂ)).sum fun _ a => ‖a‖) := by gcongr
    _ = ↑(f.natDegree.choose (f.natDegree / 2) *
          (Finset.range (f.natDegree + 1)).sum (fun j => (f.coeff j).natAbs)) := by
        rw [h5]; push_cast; ring

-- ============================================================
-- 2. Hensel uniqueness (sorry)
-- ============================================================

/-- Hensel uniqueness: if two factorizations agree mod p with coprime + B monic,
    they agree mod p^k. Proof by induction on k. -/
theorem hensel_unique
    (p : ℕ) (hp : Nat.Prime p) (k : ℕ) (hk : 0 < k)
    (F : Polynomial ℤ)
    (A₁ B₁ A₂ B₂ : Polynomial ℤ)
    (hprod₁ : Polynomial.map (Int.castRingHom (ZMod (p ^ k))) (A₁ * B₁) =
              Polynomial.map (Int.castRingHom (ZMod (p ^ k))) F)
    (hprod₂ : Polynomial.map (Int.castRingHom (ZMod (p ^ k))) (A₂ * B₂) =
              Polynomial.map (Int.castRingHom (ZMod (p ^ k))) F)
    (hA : Polynomial.map (Int.castRingHom (ZMod p)) A₁ =
          Polynomial.map (Int.castRingHom (ZMod p)) A₂)
    (hB : Polynomial.map (Int.castRingHom (ZMod p)) B₁ =
          Polynomial.map (Int.castRingHom (ZMod p)) B₂)
    (hcop : IsCoprime (Polynomial.map (Int.castRingHom (ZMod p)) A₁)
                       (Polynomial.map (Int.castRingHom (ZMod p)) B₁))
    (hB₁_monic : Monic B₁) (hB₂_monic : Monic B₂)
    : Polynomial.map (Int.castRingHom (ZMod (p ^ k))) A₁ =
      Polynomial.map (Int.castRingHom (ZMod (p ^ k))) A₂
    ∧ Polynomial.map (Int.castRingHom (ZMod (p ^ k))) B₁ =
      Polynomial.map (Int.castRingHom (ZMod (p ^ k))) B₂ := by
  sorry -- TODO: induction on k (~120 lines)
  -- Proof sketch (nl-proof §2.2):
  -- Base k=1: direct from hA, hB
  -- Step k→k+1: A₁-A₂ = p^k·E, B₁-B₂ = p^k·F
  --   p|(E·B₂+A₂·F) → Ā|Ē, B̄|F̄ → Ē=c·Ā, F̄=-c·B̄
  --   B₁ monic → lc(B₁)=1 → p|c → c=0 → p|E,F → p^{k+1}|A₁-A₂

-- ============================================================
-- 3. RecombineCorrect
-- ============================================================

theorem recombine_correct
    (f : Polynomial ℤ) (hf : f ≠ 0)
    : ∃ result : List (Polynomial ℤ), RecombineCorrect f result := by
  obtain ⟨factors, hirred, hassoc⟩ := WfDvdMonoid.exists_factors f hf
  refine ⟨factors.toList, ?_, ?_⟩
  · rw [Multiset.prod_toList]; exact hassoc.symm
  · intro g hg; exact hirred g (Multiset.mem_toList.mp hg)
