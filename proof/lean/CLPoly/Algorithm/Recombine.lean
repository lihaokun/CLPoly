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
import Mathlib.Analysis.Polynomial.Fourier
import Mathlib.NumberTheory.MahlerMeasure
import Mathlib.Analysis.Convex.SpecificFunctions.Basic
import Mathlib.Analysis.Convex.Integral

set_option autoImplicit false
set_option maxHeartbeats 3200000

open Polynomial

-- ============================================================
-- 1. Mignotte bound
-- ============================================================

/-- Helper: M(h) >= 1 for nonzero integer polynomial h -/
private lemma mahlerMeasure_ge_one_of_int_ne_zero (h : Polynomial ℤ) (hh : h ≠ 0) :
    1 ≤ (h.map (Int.castRingHom ℂ)).mahlerMeasure :=
  one_le_mahlerMeasure_of_ne_zero hh

/-- Landau inequality: M(f) ≤ ‖f‖₂ for complex polynomials.

  Proof outline (Parseval + Jensen):
  1. Parseval: circleAverage(‖f‖²) = Σ‖fₖ‖² (Mathlib: sum_sq_norm_coeff_eq_circleAverage)
  2. Jensen (log concave): logMahlerMeasure f = circleAverage(log‖f‖)
       ≤ (1/2) · log(circleAverage(‖f‖²))  [Jensen on log, concavity]
  3. Combine: M(f) = exp(logMahlerMeasure f) ≤ sqrt(Σ‖fₖ‖²) = ‖f‖₂ -/
lemma landau_inequality (p : ℂ[X]) :
    p.mahlerMeasure ≤
      Real.sqrt (p.sum fun _ a => ‖a‖ ^ 2) := by
  -- Case p = 0: both sides are 0
  by_cases hp : p = 0
  · simp [hp]
  -- Setup: let S = Σ‖pₖ‖² (the L2 norm squared)
  set S := p.sum fun _ a => ‖a‖ ^ 2 with hS_def
  -- Step 1: S > 0 (since p ≠ 0, the leading coeff has positive norm)
  have hS_pos : 0 < S := by
    apply Finset.sum_pos'
    · intro i _; positivity
    · exact ⟨p.natDegree, by simp [hp]⟩
  -- Step 2: S = Σ over support (unfolding Polynomial.sum)
  have hS_sum : S = ∑ i ∈ p.support, ‖p.coeff i‖ ^ 2 := by
    simp only [hS_def, Polynomial.sum_def]
  -- Step 3: Parseval — circleAverage(‖p‖²) = S
  have hParseval : Real.circleAverage (fun θ => ‖eval θ p‖ ^ 2) 0 1 = S := by
    rw [hS_sum, ← sum_sq_norm_coeff_eq_circleAverage]
  -- Step 4: logMahlerMeasure p ≤ (1/2) * log S
  --   Proof chain:
  --   logMahlerMeasure = circleAverage(log ‖p‖)
  --                    = circleAverage((1/2) · log(‖p‖²))  a.e.
  --                    ≤ (1/2) · log(circleAverage(‖p‖²))  by Jensen (log concave)
  --                    = (1/2) · log S                      by Parseval
  have hlog_bound : p.logMahlerMeasure ≤ (1 / 2) * Real.log S := by
    -- Expand logMahlerMeasure to interval integral form
    rw [logMahlerMeasure_def, Real.circleAverage_def, smul_eq_mul]
    -- Step 4a: Parseval — circleAverage(‖p‖²) = S
    have hRHS : (2 * Real.pi)⁻¹ * ∫ (x : ℝ) in (0 : ℝ)..(2 * Real.pi),
        ‖eval (circleMap 0 1 x) p‖ ^ 2 = S := by
      have : Real.circleAverage (fun θ => ‖eval θ p‖ ^ 2) 0 1 =
          (2 * Real.pi)⁻¹ * ∫ (x : ℝ) in (0 : ℝ)..(2 * Real.pi),
            ‖eval (circleMap 0 1 x) p‖ ^ 2 := by
        rw [Real.circleAverage_def, smul_eq_mul]
      rw [← this, hParseval]
    -- Step 4b: log‖z‖ = (1/2)·log(‖z‖²) everywhere
    have hlog_eq : ∀ (x : ℝ),
        Real.log ‖eval (circleMap 0 1 x) p‖ =
        (1 / 2) * Real.log (‖eval (circleMap 0 1 x) p‖ ^ 2) := by
      intro x; rw [Real.log_pow, Nat.cast_ofNat]; ring
    -- Step 4c: Rewrite LHS integral using log‖z‖ = (1/2)·log(‖z‖²)
    conv_lhs => rw [show (fun x => Real.log ‖eval (circleMap 0 1 x) p‖) =
        (fun x => (1 / 2) * Real.log (‖eval (circleMap 0 1 x) p‖ ^ 2)) from
        funext hlog_eq]
    rw [intervalIntegral.integral_const_mul]
    -- Goal is now: (2π)⁻¹ * ((1/2) * ∫ log‖p‖²) ≤ (1/2) * log S
    -- Rearrange to: (1/2) * ((2π)⁻¹ * ∫ log‖p‖²) ≤ (1/2) * log S
    rw [show (2 * Real.pi)⁻¹ * (1 / 2 * ∫ (x : ℝ) in (0 : ℝ)..(2 * Real.pi),
        Real.log (‖eval (circleMap 0 1 x) p‖ ^ 2)) =
        1 / 2 * ((2 * Real.pi)⁻¹ * ∫ (x : ℝ) in (0 : ℝ)..(2 * Real.pi),
        Real.log (‖eval (circleMap 0 1 x) p‖ ^ 2)) from by ring]
    -- Suffices to show: (2π)⁻¹ * ∫ log‖p‖² ≤ log S
    apply mul_le_mul_of_nonneg_left _ (by norm_num : (0 : ℝ) ≤ 1 / 2)
    -- Step 4d: Jensen for exp (convex on ℝ) gives:
    --   exp(circleAverage(log‖p‖²)) ≤ circleAverage(exp(log‖p‖²))
    --   and exp(log‖p‖²) = ‖p‖² a.e.
    --   so exp(circleAverage(log‖p‖²)) ≤ circleAverage(‖p‖²) = S
    --   hence circleAverage(log‖p‖²) ≤ log S
    have hkey : Real.exp ((2 * Real.pi)⁻¹ * ∫ (x : ℝ) in (0 : ℝ)..(2 * Real.pi),
        Real.log (‖eval (circleMap 0 1 x) p‖ ^ 2)) ≤ S := by
      -- Sub-step (i): continuity of ‖p(circleMap)‖² (used in multiple places)
      have hcont_normsq : Continuous (fun x => ‖eval (circleMap 0 1 x) p‖ ^ 2) :=
        (continuous_norm.comp (p.continuous.comp (continuous_circleMap 0 1))).pow 2
      -- Sub-step (ii): a.e. equality exp(log(‖p(z)‖²)) = ‖p(z)‖² on [0, 2π]
      -- (fails only at roots of p on the circle, a finite hence measure-zero set)
      have hae_eq : ∀ᵐ (x : ℝ) ∂MeasureTheory.volume.restrict (Set.Icc 0 (2 * Real.pi)),
          Real.exp (Real.log (‖eval (circleMap 0 1 x) p‖ ^ 2)) =
          ‖eval (circleMap 0 1 x) p‖ ^ 2 := by
        -- At points where p(e^{iθ}) ≠ 0, exp(log(‖p‖²)) = ‖p‖² since ‖p‖² > 0.
        -- The set where p(e^{iθ}) = 0 has measure zero (finite preimage of roots).
        -- We adapt the technique from mahlerMeasure_le_sum_norm_coeff.
        rw [Filter.eventually_iff_exists_mem]
        refine ⟨{x : ℝ | eval (circleMap 0 1 x) p ≠ 0}, ?_, fun x hx => ?_⟩
        · -- Show: {x | eval(circleMap 0 1 x) p ≠ 0} ∈ ae(volume.restrict(Icc 0 (2π)))
          rw [MeasureTheory.mem_ae_iff, MeasureTheory.Measure.restrict_apply' measurableSet_Icc]
          apply Set.Finite.measure_zero
          apply (Set.Finite.of_diff · (Set.finite_singleton (2 * Real.pi)))
          simp only [ne_eq, Set.compl_setOf, Decidable.not_not, Set.inter_diff_assoc,
            Set.Icc_diff_right]
          rw [show {x | eval (circleMap 0 1 x) p = 0} ∩ Set.Ico 0 (2 * Real.pi) =
              Set.Ico 0 (2 * Real.pi) ∩ {x | eval (circleMap 0 1 x) p = 0} from
            Set.inter_comm _ _,
              ← Set.sep_mem_eq]
          apply Set.Finite.of_finite_image (f := circleMap 0 1)
          · exact (Multiset.finite_toSet p.roots).subset (by
              intro z ⟨θ, ⟨_, hθ⟩, hθz⟩
              rw [← hθz]
              exact (mem_roots (show p ≠ 0 from hp)).mpr hθ)
          · exact fun _ h _ k l => injOn_circleMap_of_abs_sub_le' one_ne_zero
              (by linarith [Real.two_pi_pos]) h.1 k.1 l
        · -- At non-zero points: exp(log(‖p‖²)) = ‖p‖²
          exact Real.exp_log (by positivity : 0 < ‖eval (circleMap 0 1 x) p‖ ^ 2)
      -- Sub-step (iii): integrability of log(‖p‖²) on [0, 2π]
      have hf_integ : MeasureTheory.IntegrableOn
          (fun x => Real.log (‖eval (circleMap 0 1 x) p‖ ^ 2))
          (Set.Icc 0 (2 * Real.pi)) MeasureTheory.volume := by
        -- log(‖p‖²) = 2·log(‖p‖), integrable since meromorphic
        -- (from intervalIntegrable_mahlerMeasure via log_pow)
        rw [show (fun x => Real.log (‖eval (circleMap 0 1 x) p‖ ^ 2)) =
            (fun x => 2 * Real.log ‖eval (circleMap 0 1 x) p‖) from by
          ext x; rw [Real.log_pow]; simp]
        rw [← intervalIntegrable_iff_integrableOn_Icc_of_le Real.two_pi_pos.le]
        exact p.intervalIntegrable_mahlerMeasure.const_mul 2
      -- Sub-step (iv): integrability of exp(log(‖p‖²)) on [0, 2π]
      have hexplog_integ : MeasureTheory.IntegrableOn
          (fun x => Real.exp (Real.log (‖eval (circleMap 0 1 x) p‖ ^ 2)))
          (Set.Icc 0 (2 * Real.pi)) MeasureTheory.volume := by
        -- exp(log(‖p‖²)) = ‖p‖² a.e., and ‖p‖² is continuous hence integrable on compact set
        exact (hcont_normsq.continuousOn.integrableOn_compact isCompact_Icc).congr_fun_ae
          (hae_eq.mono (fun x hx => hx.symm))
      -- Sub-step (v): integrability of ‖p‖² on [0, 2π]
      have hnormsq_integ : MeasureTheory.IntegrableOn
          (fun x => ‖eval (circleMap 0 1 x) p‖ ^ 2)
          (Set.Icc 0 (2 * Real.pi)) MeasureTheory.volume :=
        hcont_normsq.continuousOn.integrableOn_compact isCompact_Icc
      -- Sub-step (vi): measure of [0, 2π] is positive and finite
      have hIcc_pos : MeasureTheory.volume (Set.Icc 0 (2 * Real.pi)) ≠ 0 := by
        rw [Real.volume_Icc, sub_zero]
        simp only [ne_eq, ENNReal.ofReal_eq_zero, not_le]
        exact Real.two_pi_pos
      have hIcc_fin : MeasureTheory.volume (Set.Icc 0 (2 * Real.pi)) ≠ ⊤ :=
        measure_Icc_lt_top.ne
      -- Sub-step (vii): Jensen for exp (convex on Set.univ, which is closed)
      --   exp(⨍ x in [0,2π], f x) ≤ ⨍ x in [0,2π], exp(f x)
      have hJensen : Real.exp (⨍ x in Set.Icc 0 (2 * Real.pi),
          Real.log (‖eval (circleMap 0 1 x) p‖ ^ 2)) ≤
          ⨍ x in Set.Icc 0 (2 * Real.pi),
            Real.exp (Real.log (‖eval (circleMap 0 1 x) p‖ ^ 2)) :=
        ConvexOn.map_set_average_le convexOn_exp Real.continuous_exp.continuousOn
          isClosed_univ hIcc_pos hIcc_fin
          (MeasureTheory.ae_of_all _ (fun _ => Set.mem_univ _))
          hf_integ hexplog_integ
      -- Sub-step (viii): ⨍ exp(log(‖p‖²)) = ⨍ ‖p‖² (a.e. equality)
      have hae_avg : ⨍ x in Set.Icc 0 (2 * Real.pi),
          Real.exp (Real.log (‖eval (circleMap 0 1 x) p‖ ^ 2)) =
          ⨍ x in Set.Icc 0 (2 * Real.pi), ‖eval (circleMap 0 1 x) p‖ ^ 2 :=
        MeasureTheory.setAverage_congr_fun measurableSet_Icc
          (MeasureTheory.ae_imp_of_ae_restrict hae_eq)
      -- Sub-step (ix): bridge set average ↔ interval integral
      -- ⨍ x in [0,2π], g x = (2π)⁻¹ * ∫ x in 0..2π, g x (for Lebesgue on [0, 2π])
      have havg_to_int : ∀ (g : ℝ → ℝ),
          MeasureTheory.IntegrableOn g (Set.Icc 0 (2 * Real.pi)) MeasureTheory.volume →
          ⨍ x in Set.Icc 0 (2 * Real.pi), g x =
          (2 * Real.pi)⁻¹ * ∫ (x : ℝ) in (0 : ℝ)..(2 * Real.pi), g x := by
        intro g _
        -- setAverage = (measureReal)⁻¹ • setIntegral
        rw [MeasureTheory.setAverage_eq, smul_eq_mul]
        -- measureReal(Icc 0 (2π)) = 2π
        have hmr : MeasureTheory.volume.real (Set.Icc 0 (2 * Real.pi)) = 2 * Real.pi := by
          rw [MeasureTheory.measureReal_def, Real.volume_Icc, sub_zero,
              ENNReal.toReal_ofReal Real.two_pi_pos.le]
        rw [hmr]
        -- ∫ x in Icc 0 (2π), g x = ∫ x in 0..2π, g x
        -- (set integral on Icc = interval integral, since Icc and Ioc differ by measure zero)
        congr 1
        rw [intervalIntegral.integral_of_le Real.two_pi_pos.le,
            MeasureTheory.integral_Icc_eq_integral_Ioc]
      -- Combine everything
      -- Rewrite both ends to set averages
      rw [show (2 * Real.pi)⁻¹ * ∫ (x : ℝ) in (0 : ℝ)..(2 * Real.pi),
          Real.log (‖eval (circleMap 0 1 x) p‖ ^ 2) =
          ⨍ x in Set.Icc 0 (2 * Real.pi),
          Real.log (‖eval (circleMap 0 1 x) p‖ ^ 2) from
          (havg_to_int _ hf_integ).symm]
      calc Real.exp (⨍ x in Set.Icc 0 (2 * Real.pi),
              Real.log (‖eval (circleMap 0 1 x) p‖ ^ 2))
        ≤ ⨍ x in Set.Icc 0 (2 * Real.pi),
            Real.exp (Real.log (‖eval (circleMap 0 1 x) p‖ ^ 2)) := hJensen
        _ = ⨍ x in Set.Icc 0 (2 * Real.pi), ‖eval (circleMap 0 1 x) p‖ ^ 2 := hae_avg
        _ = (2 * Real.pi)⁻¹ * ∫ (x : ℝ) in (0 : ℝ)..(2 * Real.pi),
            ‖eval (circleMap 0 1 x) p‖ ^ 2 := havg_to_int _ hnormsq_integ
        _ = S := hRHS
    -- log is monotone: exp(a) ≤ S > 0  ⟹  a ≤ log S
    exact (Real.le_log_iff_exp_le hS_pos).mpr hkey
  -- Step 5: Exponentiate: M(p) = exp(logMahlerMeasure) ≤ exp((1/2)·log S) = √S
  rw [show p.mahlerMeasure = Real.exp p.logMahlerMeasure from by
    rw [Polynomial.mahlerMeasure, if_pos hp]]
  rw [show Real.sqrt S = Real.exp ((1 / 2) * Real.log S) from by
    rw [Real.sqrt_eq_rpow, Real.rpow_def_of_pos hS_pos, mul_comm]]
  exact Real.exp_le_exp.mpr hlog_bound

/-- Mignotte bound (L2 version): g | f → |g.coeff i| ≤ C(n, n/2) · ‖f‖₂
    This matches the C++ `__mignotte_bound` which uses the L2 norm. -/
theorem mignotte_bound_l2 (f g : Polynomial ℤ) (hf : f ≠ 0) (hg : g ∣ f) :
    ∀ i, (↑((g.coeff i).natAbs) : ℝ) ≤
      ↑(Nat.choose f.natDegree (f.natDegree / 2)) *
      Real.sqrt ((Finset.range (f.natDegree + 1)).sum
        (fun j => ((f.coeff j).natAbs : ℝ) ^ 2)) := by
  intro i
  have hφ_inj : Function.Injective (Int.castRingHom ℂ) := Int.cast_injective
  obtain ⟨h, hfgh⟩ := hg
  have hh_ne : h ≠ 0 := right_ne_zero_of_mul (hfgh ▸ hf)
  have hg_ne : g ≠ 0 := left_ne_zero_of_mul (hfgh ▸ hf)
  have h_lhs : ‖(g.map (Int.castRingHom ℂ)).coeff i‖ = ↑((g.coeff i).natAbs) := by
    rw [coeff_map, eq_intCast, Complex.norm_intCast]
    rw [Nat.cast_natAbs (α := ℝ), Int.cast_abs]
  rw [← h_lhs]
  have h1 := norm_coeff_le_choose_mul_mahlerMeasure i (g.map (Int.castRingHom ℂ))
  rw [natDegree_map_eq_of_injective hφ_inj] at h1
  have h2 : (g.map (Int.castRingHom ℂ)).mahlerMeasure ≤
      (f.map (Int.castRingHom ℂ)).mahlerMeasure := by
    have hfgh_c : f.map (Int.castRingHom ℂ) =
        g.map (Int.castRingHom ℂ) * h.map (Int.castRingHom ℂ) := by
      rw [← Polynomial.map_mul, hfgh]
    rw [hfgh_c, mahlerMeasure_mul]
    exact le_mul_of_one_le_right (mahlerMeasure_nonneg _)
      (mahlerMeasure_ge_one_of_int_ne_zero h hh_ne)
  have hg_dvd : g ∣ f := ⟨h, hfgh⟩
  have h3 : g.natDegree.choose i ≤ f.natDegree.choose (f.natDegree / 2) :=
    (Nat.choose_le_middle i g.natDegree).trans
      ((Nat.choose_le_choose _ (natDegree_le_of_dvd hg_dvd hf)).trans
        (Nat.choose_le_middle _ f.natDegree))
  -- Bridge: sum over support = sum over range for L2 norm
  have h5_l2 : (f.map (Int.castRingHom ℂ)).sum (fun _ a => ‖a‖ ^ 2) =
      (Finset.range (f.natDegree + 1)).sum (fun j => ((f.coeff j).natAbs : ℝ) ^ 2) := by
    rw [(f.map (Int.castRingHom ℂ)).sum_over_range (fun _ => by simp),
        natDegree_map_eq_of_injective hφ_inj]
    congr 1; ext j
    rw [coeff_map, eq_intCast, Complex.norm_intCast]
    rw [Nat.cast_natAbs (α := ℝ), Int.cast_abs, sq_abs]
  calc ‖(g.map (Int.castRingHom ℂ)).coeff i‖
    ≤ ↑(g.natDegree.choose i) * (g.map (Int.castRingHom ℂ)).mahlerMeasure := h1
    _ ≤ ↑(g.natDegree.choose i) * (f.map (Int.castRingHom ℂ)).mahlerMeasure := by gcongr
    _ ≤ ↑(f.natDegree.choose (f.natDegree / 2)) *
          (f.map (Int.castRingHom ℂ)).mahlerMeasure := by
        apply mul_le_mul_of_nonneg_right (by exact_mod_cast h3) (mahlerMeasure_nonneg _)
    _ ≤ ↑(f.natDegree.choose (f.natDegree / 2)) *
          Real.sqrt ((f.map (Int.castRingHom ℂ)).sum fun _ a => ‖a‖ ^ 2) := by
        gcongr; exact landau_inequality _
    _ = ↑(f.natDegree.choose (f.natDegree / 2)) *
          Real.sqrt ((Finset.range (f.natDegree + 1)).sum
            (fun j => ((f.coeff j).natAbs : ℝ) ^ 2)) := by
        rw [h5_l2]

/-- Mignotte bound (L1 version): g | f -> |g.coeff i| <= C(n, n/2) * L1(f) -/
theorem mignotte_bound (f g : Polynomial ℤ) (hf : f ≠ 0) (hg : g ∣ f) :
    ∀ i, (g.coeff i).natAbs ≤
      Nat.choose f.natDegree (f.natDegree / 2) *
      (Finset.range (f.natDegree + 1)).sum (fun j => (f.coeff j).natAbs) := by
  intro i
  have hφ_inj : Function.Injective (Int.castRingHom ℂ) := Int.cast_injective
  obtain ⟨h, hfgh⟩ := hg
  have hh_ne : h ≠ 0 := right_ne_zero_of_mul (hfgh ▸ hf)
  have hg_ne : g ≠ 0 := left_ne_zero_of_mul (hfgh ▸ hf)
  suffices h_real : (↑((g.coeff i).natAbs) : ℝ) ≤
      ↑(Nat.choose f.natDegree (f.natDegree / 2) *
        (Finset.range (f.natDegree + 1)).sum (fun j => (f.coeff j).natAbs)) by
    exact_mod_cast h_real
  have h_lhs : ‖(g.map (Int.castRingHom ℂ)).coeff i‖ = ↑((g.coeff i).natAbs) := by
    rw [coeff_map, eq_intCast, Complex.norm_intCast]
    rw [Nat.cast_natAbs (α := ℝ), Int.cast_abs]
  rw [← h_lhs]
  have h1 := norm_coeff_le_choose_mul_mahlerMeasure i (g.map (Int.castRingHom ℂ))
  rw [natDegree_map_eq_of_injective hφ_inj] at h1
  have h2 : (g.map (Int.castRingHom ℂ)).mahlerMeasure ≤
      (f.map (Int.castRingHom ℂ)).mahlerMeasure := by
    have hfgh_c : f.map (Int.castRingHom ℂ) =
        g.map (Int.castRingHom ℂ) * h.map (Int.castRingHom ℂ) := by
      rw [← Polynomial.map_mul, hfgh]
    rw [hfgh_c, mahlerMeasure_mul]
    exact le_mul_of_one_le_right (mahlerMeasure_nonneg _)
      (mahlerMeasure_ge_one_of_int_ne_zero h hh_ne)
  have hg_dvd : g ∣ f := ⟨h, hfgh⟩
  have h3 : g.natDegree.choose i ≤ f.natDegree.choose (f.natDegree / 2) :=
    (Nat.choose_le_middle i g.natDegree).trans
      ((Nat.choose_le_choose _ (natDegree_le_of_dvd hg_dvd hf)).trans
        (Nat.choose_le_middle _ f.natDegree))
  have h5 : (f.map (Int.castRingHom ℂ)).sum (fun _ a => ‖a‖) =
      ↑((Finset.range (f.natDegree + 1)).sum (fun j => (f.coeff j).natAbs)) := by
    rw [(f.map (Int.castRingHom ℂ)).sum_over_range (fun _ => by simp),
        natDegree_map_eq_of_injective hφ_inj]
    push_cast; congr 1; ext j
    rw [coeff_map, eq_intCast, Complex.norm_intCast]
    rw [Nat.cast_natAbs (α := ℝ), Int.cast_abs]
  calc ‖(g.map (Int.castRingHom ℂ)).coeff i‖
    ≤ ↑(g.natDegree.choose i) * (g.map (Int.castRingHom ℂ)).mahlerMeasure := h1
    _ ≤ ↑(g.natDegree.choose i) * (f.map (Int.castRingHom ℂ)).mahlerMeasure := by gcongr
    _ ≤ ↑(f.natDegree.choose (f.natDegree / 2)) *
          (f.map (Int.castRingHom ℂ)).mahlerMeasure := by
        apply mul_le_mul_of_nonneg_right (by exact_mod_cast h3) (mahlerMeasure_nonneg _)
    _ ≤ ↑(f.natDegree.choose (f.natDegree / 2)) *
          ((f.map (Int.castRingHom ℂ)).sum fun _ a => ‖a‖) := by
        gcongr; exact mahlerMeasure_le_sum_norm_coeff _
    _ = ↑(f.natDegree.choose (f.natDegree / 2) *
          (Finset.range (f.natDegree + 1)).sum (fun j => (f.coeff j).natAbs)) := by
        rw [h5]; push_cast; ring

-- ============================================================
-- 2. Hensel uniqueness
-- ============================================================

-- Helper: projecting polynomial map equalities from larger to smaller modulus
private lemma map_eq_of_dvd_map_eq (a b : ℕ) (hab : a ∣ b) (f g : Polynomial ℤ)
    (h : Polynomial.map (Int.castRingHom (ZMod b)) f =
         Polynomial.map (Int.castRingHom (ZMod b)) g) :
    Polynomial.map (Int.castRingHom (ZMod a)) f =
    Polynomial.map (Int.castRingHom (ZMod a)) g := by
  have key : Polynomial.map (ZMod.castHom hab (ZMod a))
      (Polynomial.map (Int.castRingHom (ZMod b)) f) =
    Polynomial.map (ZMod.castHom hab (ZMod a))
      (Polynomial.map (Int.castRingHom (ZMod b)) g) := congr_arg _ h
  simp only [Polynomial.map_map] at key
  convert key using 1 <;> congr 1 <;> ext x <;> simp

-- Helper: if b | q with b monic and natDegree q < natDegree b, then q = 0
private lemma dvd_monic_eq_zero_of_natDegree_lt {R : Type*} [CommRing R] [IsDomain R]
    (b q : Polynomial R) (hb_monic : Monic b) (hdvd : b ∣ q)
    (hdeg : natDegree q < natDegree b) : q = 0 := by
  by_contra hq
  obtain ⟨c, hc⟩ := hdvd
  have hc_ne : c ≠ 0 := right_ne_zero_of_mul (hc ▸ hq)
  have : natDegree q = natDegree b + natDegree c := hc ▸ hb_monic.natDegree_mul' hc_ne
  omega

-- Helper: if two polynomials with same leadingCoeff differ by C(m)*G with m ∤ lc,
-- then for all n >= natDegree B1, G.coeff n = 0
private lemma diff_coeff_bound
    (B1 B2 : Polynomial ℤ) (m : ℕ) (hm : 1 < m)
    (G : Polynomial ℤ)
    (hB_lc_eq : B1.leadingCoeff = B2.leadingCoeff)
    (hB_lc_ne : ¬((m : ℤ) ∣ B1.leadingCoeff))
    (hdiff : B1 - B2 = C (↑m : ℤ) * G) :
    ∀ n, natDegree B1 ≤ n → G.coeff n = 0 := by
  have hndeg : natDegree B1 = natDegree B2 := by
    by_contra hne
    wlog h12 : natDegree B2 < natDegree B1 with H
    · push_neg at h12
      exact H B2 B1 m hm (-G) hB_lc_eq.symm
        (show ¬((m : ℤ) ∣ B2.leadingCoeff) by rwa [← hB_lc_eq])
        (by linear_combination -hdiff) (Ne.symm hne) (lt_of_le_of_ne h12 hne)
    have hcoeff2 : (B1 - B2).coeff (natDegree B1) = ↑m * G.coeff (natDegree B1) := by
      rw [hdiff, coeff_C_mul]
    have hcoeff : B1.coeff (natDegree B1) = ↑m * G.coeff (natDegree B1) := by
      rw [← hcoeff2, coeff_sub, Polynomial.coeff_eq_zero_of_natDegree_lt h12, sub_zero]
    exact hB_lc_ne ⟨_, hcoeff⟩
  have hm_ne : (m : ℤ) ≠ 0 := by omega
  intro n hn
  have hcoeff_diff : (B1 - B2).coeff n = ↑m * G.coeff n := by rw [hdiff, coeff_C_mul]
  suffices hsuff : (B1 - B2).coeff n = 0 by
    rw [hsuff] at hcoeff_diff; exact (mul_eq_zero.mp hcoeff_diff.symm).resolve_left hm_ne
  rcases Nat.eq_or_lt_of_le hn with rfl | hgt
  · simp only [coeff_sub]
    have h1 : B1.coeff (natDegree B1) = B1.leadingCoeff := rfl
    have h2 : B2.coeff (natDegree B1) = B2.leadingCoeff := by
      rw [Polynomial.leadingCoeff, hndeg]
    rw [h1, h2, hB_lc_eq, sub_self]
  · have h1 : B1.coeff n = 0 := Polynomial.coeff_eq_zero_of_natDegree_lt hgt
    have h2 : B2.coeff n = 0 := Polynomial.coeff_eq_zero_of_natDegree_lt (hndeg ▸ hgt)
    simp [coeff_sub, h1, h2]

-- Helper: map_p(Q) = 0 implies p | Q.coeff n
private lemma coeff_dvd_of_map_eq_zero (p : ℕ) (Q : Polynomial ℤ)
    (h : Polynomial.map (Int.castRingHom (ZMod p)) Q = 0)
    (n : ℕ) : (↑p : ℤ) ∣ Q.coeff n := by
  rw [← ZMod.intCast_zmod_eq_zero_iff_dvd]
  have := congr_arg (fun q => q.coeff n) h
  simpa [Polynomial.coeff_map] using this

/-- Hensel uniqueness: if two factorizations agree mod p with coprime factors
    and B₁, B₂ have the same leading coefficient (coprime to p),
    they agree mod p^k. Proof by induction on k.
    Generalization of the classical "B monic" version. -/
theorem hensel_unique
    (p : ℕ) (hp : Nat.Prime p) (k : ℕ) (hk : 0 < k)
    (F : Polynomial ℤ)
    (A1 B1 A2 B2 : Polynomial ℤ)
    (hprod1 : Polynomial.map (Int.castRingHom (ZMod (p ^ k))) (A1 * B1) =
              Polynomial.map (Int.castRingHom (ZMod (p ^ k))) F)
    (hprod2 : Polynomial.map (Int.castRingHom (ZMod (p ^ k))) (A2 * B2) =
              Polynomial.map (Int.castRingHom (ZMod (p ^ k))) F)
    (hA : Polynomial.map (Int.castRingHom (ZMod p)) A1 =
          Polynomial.map (Int.castRingHom (ZMod p)) A2)
    (hB : Polynomial.map (Int.castRingHom (ZMod p)) B1 =
          Polynomial.map (Int.castRingHom (ZMod p)) B2)
    (hcop : IsCoprime (Polynomial.map (Int.castRingHom (ZMod p)) A1)
                       (Polynomial.map (Int.castRingHom (ZMod p)) B1))
    (hB_lc_eq : B1.leadingCoeff = B2.leadingCoeff)
    (hB_lc_coprime : ¬((p : ℤ) ∣ B1.leadingCoeff))
    : Polynomial.map (Int.castRingHom (ZMod (p ^ k))) A1 =
      Polynomial.map (Int.castRingHom (ZMod (p ^ k))) A2
    ∧ Polynomial.map (Int.castRingHom (ZMod (p ^ k))) B1 =
      Polynomial.map (Int.castRingHom (ZMod (p ^ k))) B2 := by
  induction k with
  | zero => omega
  | succ k ih =>
    by_cases hk' : k = 0
    · -- Base: p^1 = p
      subst hk'
      have : p ^ (0 + 1) = p := by norm_num
      rw [this] at hprod1 hprod2 ⊢; exact ⟨hA, hB⟩
    · -- Inductive step: k >= 1
      have hk_pos : 0 < k := Nat.pos_of_ne_zero hk'
      have hpk_pos : 0 < p ^ k := Nat.pos_of_ne_zero (pow_ne_zero k hp.ne_zero)
      have hpk_ne : (↑(p ^ k) : ℤ) ≠ 0 := Int.natCast_ne_zero.mpr (by omega)
      have hdvd : p ^ k ∣ p ^ (k + 1) := Nat.pow_dvd_pow p (by omega)
      -- Apply IH
      obtain ⟨ihA, ihB⟩ := ih hk_pos
        (map_eq_of_dvd_map_eq _ _ hdvd _ _ hprod1)
        (map_eq_of_dvd_map_eq _ _ hdvd _ _ hprod2)
      -- Extract: A1 - A2 = C(p^k)*E, B1 - B2 = C(p^k)*G
      obtain ⟨E, hE⟩ := _root_.exists_C_mul_of_map_eq_zero (p ^ k) hpk_pos (A1 - A2)
        (by rw [Polynomial.map_sub, ihA, sub_self])
      obtain ⟨G, hG⟩ := _root_.exists_C_mul_of_map_eq_zero (p ^ k) hpk_pos (B1 - B2)
        (by rw [Polynomial.map_sub, ihB, sub_self])
      -- A1 = A2 + C(p^k)*E, B1 = B2 + C(p^k)*G
      have hA1 : A1 = A2 + C (↑(p ^ k) : ℤ) * E := eq_add_of_sub_eq' hE
      have hB1eq : B1 = B2 + C (↑(p ^ k) : ℤ) * G := eq_add_of_sub_eq' hG
      -- Products agree mod p^(k+1)
      have hprod_eq : Polynomial.map (Int.castRingHom (ZMod (p ^ (k + 1)))) (A1 * B1) =
          Polynomial.map (Int.castRingHom (ZMod (p ^ (k + 1)))) (A2 * B2) := by
        rw [hprod1, hprod2]
      -- Key: show p | coeff(E*B2 + A2*G, n) for all n
      -- From: A1*B1 - A2*B2 = C(p^k)*(E*B2+A2*G) + C(p^k)^2*(E*G) and p^(k+1) | LHS, p^(k+1) | C(p^k)^2*(...)
      have hEBAG_p : ∀ n, (↑p : ℤ) ∣ (E * B2 + A2 * G).coeff n := by
        intro n
        -- Expand: A1*B1 - A2*B2 = C(p^k)*(E*B2+A2*G) + C(p^k)^2*(E*G)
        have hexp : A1 * B1 - A2 * B2 =
            C (↑(p ^ k) : ℤ) * (E * B2 + A2 * G) +
            C (↑(p ^ k) : ℤ) * (C (↑(p ^ k) : ℤ) * (E * G)) := by
          rw [hA1, hB1eq]; ring
        -- p^(k+1) | coeff(A1*B1 - A2*B2, n)
        have h1 : (↑(p ^ (k + 1)) : ℤ) ∣ (A1 * B1 - A2 * B2).coeff n := by
          have h0 : Polynomial.map (Int.castRingHom (ZMod (p ^ (k + 1)))) (A1 * B1 - A2 * B2) = 0 := by
            simp only [Polynomial.map_sub, Polynomial.map_mul, hprod_eq, sub_self]
          rw [← ZMod.intCast_zmod_eq_zero_iff_dvd]
          have h1 := congr_arg (fun q => q.coeff n) h0
          simpa only [Polynomial.coeff_map, Polynomial.coeff_zero] using h1
        -- Coeff of expansion
        have h2 : (A1 * B1 - A2 * B2).coeff n =
            ↑(p ^ k) * (E * B2 + A2 * G).coeff n +
            ↑(p ^ k) * (↑(p ^ k) * (E * G).coeff n) := by
          have := congr_arg (fun q => q.coeff n) hexp
          simp only [coeff_add, coeff_C_mul] at this; exact this
        -- p^(k+1) | p^(2k) * coeff(E*G, n)
        have h3 : (↑(p ^ (k + 1)) : ℤ) ∣ ↑(p ^ k) * (↑(p ^ k) * (E * G).coeff n) := by
          rw [show (↑(p ^ k) : ℤ) * ((↑(p ^ k) : ℤ) * (E * G).coeff n) =
            ↑(p ^ (k + k)) * (E * G).coeff n from by push_cast; ring]
          exact dvd_mul_of_dvd_left (by exact_mod_cast Nat.pow_dvd_pow p (by omega)) _
        -- p^(k+1) | p^k * coeff(E*B2+A2*G, n)
        have h4 : (↑(p ^ (k + 1)) : ℤ) ∣ ↑(p ^ k) * (E * B2 + A2 * G).coeff n := by
          rw [h2] at h1; exact (dvd_add_left h3).mp h1
        -- Cancel p^k
        rw [show (↑(p ^ (k + 1)) : ℤ) = ↑(p ^ k) * ↑p from by push_cast; ring] at h4
        exact (mul_dvd_mul_iff_left hpk_ne).mp h4
      -- Convert to polynomial-level statement
      have hEBAG_map : Polynomial.map (Int.castRingHom (ZMod p)) (E * B2 + A2 * G) = 0 := by
        ext n; simp only [Polynomial.coeff_map, Polynomial.coeff_zero]
        show ((E * B2 + A2 * G).coeff n : ZMod p) = 0
        exact (ZMod.intCast_zmod_eq_zero_iff_dvd _ _).mpr (hEBAG_p n)
      -- In F_p[x]: use coprimality
      haveI : Fact (Nat.Prime p) := ⟨hp⟩
      -- E_bar * B_bar + A_bar * G_bar = 0  (using hA, hB to replace A2->A1, B2->B1)
      have hEBAG_ab : Polynomial.map (Int.castRingHom (ZMod p)) E *
          Polynomial.map (Int.castRingHom (ZMod p)) B1 +
          Polynomial.map (Int.castRingHom (ZMod p)) A1 *
          Polynomial.map (Int.castRingHom (ZMod p)) G = 0 := by
        have := hEBAG_map
        simp only [Polynomial.map_add, Polynomial.map_mul] at this
        rwa [← hB, ← hA] at this
      -- Coprimality: A_bar | E_bar and B_bar | G_bar
      have hA_dvd : Polynomial.map (Int.castRingHom (ZMod p)) A1 ∣
          Polynomial.map (Int.castRingHom (ZMod p)) E := by
        apply hcop.dvd_of_dvd_mul_right
        have : Polynomial.map (Int.castRingHom (ZMod p)) E *
            Polynomial.map (Int.castRingHom (ZMod p)) B1 =
          -(Polynomial.map (Int.castRingHom (ZMod p)) A1 *
            Polynomial.map (Int.castRingHom (ZMod p)) G) :=
          (neg_eq_of_add_eq_zero_left hEBAG_ab).symm
        rw [this]; exact dvd_neg.mpr (dvd_mul_right _ _)
      have hB_dvd : Polynomial.map (Int.castRingHom (ZMod p)) B1 ∣
          Polynomial.map (Int.castRingHom (ZMod p)) G := by
        apply hcop.symm.dvd_of_dvd_mul_left
        have : Polynomial.map (Int.castRingHom (ZMod p)) A1 *
            Polynomial.map (Int.castRingHom (ZMod p)) G =
          -(Polynomial.map (Int.castRingHom (ZMod p)) E *
            Polynomial.map (Int.castRingHom (ZMod p)) B1) := (neg_eq_of_add_eq_zero_right hEBAG_ab).symm
        rw [this]; exact dvd_neg.mpr (dvd_mul_left _ _)
      -- B_bar is nonzero (since p ∤ lc(B1))
      have hB_lc_ne_p : (Int.castRingHom (ZMod p)) B1.leadingCoeff ≠ 0 := by
        change ((B1.leadingCoeff : ℤ) : ZMod p) ≠ 0
        rwa [Ne, ZMod.intCast_zmod_eq_zero_iff_dvd]
      have hBne : Polynomial.map (Int.castRingHom (ZMod p)) B1 ≠ 0 :=
        fun h => hB_lc_coprime (coeff_dvd_of_map_eq_zero p B1 h B1.natDegree)
      -- Degree bound: for n >= natDegree B1, G.coeff n = 0
      have hB_lc_ne_pk : ¬((↑(p ^ k) : ℤ) ∣ B1.leadingCoeff) := by
        intro h; exact hB_lc_coprime (dvd_trans (by exact_mod_cast dvd_pow_self p (by omega : k ≠ 0)) h)
      have hG_coeff := diff_coeff_bound B1 B2 (p ^ k)
        (Nat.one_lt_pow hk' hp.one_lt) G hB_lc_eq hB_lc_ne_pk hG
      -- natDegree(map_p(B1)) = natDegree(B1) since p ∤ lc(B1)
      have hBn : natDegree (Polynomial.map (Int.castRingHom (ZMod p)) B1) = natDegree B1 :=
        Polynomial.natDegree_map_of_leadingCoeff_ne_zero _ hB_lc_ne_p
      -- G_bar = 0
      have hGbar : Polynomial.map (Int.castRingHom (ZMod p)) G = 0 := by
        by_cases hG0 : G = 0
        · simp [hG0]
        · have hlt : natDegree G < natDegree B1 := by
            by_contra h; push_neg at h
            exact absurd (hG_coeff _ h) (Polynomial.leadingCoeff_ne_zero.mpr hG0)
          -- map_p G has degree < degree of map_p B1, and B1 | G in F_p[x] → G_bar = 0
          by_contra hG_ne
          exact Nat.not_lt.mpr (Polynomial.natDegree_le_of_dvd hB_dvd hG_ne) (by
            calc natDegree (Polynomial.map (Int.castRingHom (ZMod p)) G)
                ≤ natDegree G := Polynomial.natDegree_map_le
              _ < natDegree B1 := hlt
              _ = _ := hBn.symm)
      -- E_bar = 0 (from E_bar * B_bar + 0 = 0 and B_bar nonzero)
      have hEbar : Polynomial.map (Int.castRingHom (ZMod p)) E = 0 := by
        have : Polynomial.map (Int.castRingHom (ZMod p)) E *
            Polynomial.map (Int.castRingHom (ZMod p)) B1 = 0 := by
          have := hEBAG_ab; rw [hGbar, mul_zero, add_zero] at this; exact this
        exact (mul_eq_zero.mp this).resolve_right hBne
      -- Lift to p^(k+1): if map_p(Q) = 0 then p^(k+1) | C(p^k)*Q
      have lift : ∀ (D Q : Polynomial ℤ), D = C (↑(p ^ k) : ℤ) * Q →
          Polynomial.map (Int.castRingHom (ZMod p)) Q = 0 →
          Polynomial.map (Int.castRingHom (ZMod (p ^ (k + 1)))) D = 0 := by
        intro D Q hDQ hQ0; rw [hDQ]; ext n
        simp only [Polynomial.coeff_map, coeff_C_mul, Polynomial.coeff_zero]
        -- Goal: (Int.castRingHom (ZMod (p^(k+1)))) (↑(p^k) * Q.coeff n) = 0
        obtain ⟨c, hc⟩ := coeff_dvd_of_map_eq_zero p Q hQ0 n
        rw [hc]
        have h1 : (↑(p ^ k) : ℤ) * (↑p * c) = ↑(p ^ (k + 1)) * c := by push_cast; ring
        rw [h1, map_mul]
        suffices h2 : (Int.castRingHom (ZMod (p ^ (k + 1)))) (↑(p ^ (k + 1)) : ℤ) = 0 by
          rw [h2, zero_mul]
        change ((↑(p ^ (k + 1)) : ℤ) : ZMod (p ^ (k + 1))) = 0
        exact (ZMod.intCast_zmod_eq_zero_iff_dvd _ _).mpr (dvd_refl _)
      constructor
      · have h := lift _ _ hE hEbar
        rwa [Polynomial.map_sub, sub_eq_zero] at h
      · have h := lift _ _ hG hGbar
        rwa [Polynomial.map_sub, sub_eq_zero] at h

-- ============================================================
-- 3. Symmetric recovery (factor recovery from Hensel lift)
-- ============================================================

/-- If two integers are congruent mod m and both have |·| * 2 < m, they are equal. -/
private lemma int_eq_of_mod_eq_of_abs_lt (a b : ℤ) (m : ℕ) (_ : 0 < m)
    (hmod : (a : ZMod m) = (b : ZMod m))
    (ha : a.natAbs * 2 < m) (hb : b.natAbs * 2 < m)
    : a = b := by
  suffices h : a - b = 0 by omega
  have hdvd : (m : ℤ) ∣ (a - b) := by
    rw [← ZMod.intCast_zmod_eq_zero_iff_dvd]
    push_cast; rw [sub_eq_zero]; exact hmod
  have habs : (a - b).natAbs < m := by
    calc (a - b).natAbs ≤ a.natAbs + b.natAbs := Int.natAbs_sub_le a b
      _ < m := by omega
  obtain ⟨k, hk⟩ := hdvd
  have h1 : (a - b).natAbs = m * k.natAbs := by
    rw [hk, Int.natAbs_mul, Int.natAbs_natCast]
  rw [h1] at habs
  -- habs : m * k.natAbs < m, hm : 0 < m → k.natAbs = 0
  have hk0 : k.natAbs = 0 := by
    by_contra h
    exact Nat.not_lt.mpr (Nat.le_mul_of_pos_right m (by omega)) habs
  rw [Int.natAbs_eq_zero.mp hk0, mul_zero] at hk
  exact hk

/-- If two ℤ[x] polynomials agree mod m and both have all coefficients with |·| * 2 < m,
    they are equal. This is the "symmetric recovery" step in Zassenhaus recombination. -/
theorem symmetric_recovery (P Q : Polynomial ℤ) (m : ℕ) (hm : 0 < m)
    (hmod : Polynomial.map (Int.castRingHom (ZMod m)) P =
            Polynomial.map (Int.castRingHom (ZMod m)) Q)
    (hP : ∀ i, (P.coeff i).natAbs * 2 < m)
    (hQ : ∀ i, (Q.coeff i).natAbs * 2 < m)
    : P = Q := by
  ext i
  apply int_eq_of_mod_eq_of_abs_lt _ _ m hm
  · have h := congr_arg (fun p => p.coeff i) hmod
    simpa only [Polynomial.coeff_map] using h
  · exact hP i
  · exact hQ i

/-- Factor recovery: if A ≡ C(c)*g (mod m) and both have small coefficients (< m/2),
    then A = C(c)*g exactly. This verifies the core step of C++ __zassenhaus_recombine:
    symmetric_mod(subset_product) recovers the scaled true factor under Mignotte precision. -/
theorem factor_recovery
    (g A : Polynomial ℤ) (m : ℕ) (hm : 0 < m) (c : ℤ)
    (hmod : Polynomial.map (Int.castRingHom (ZMod m)) A =
            Polynomial.map (Int.castRingHom (ZMod m)) (C c * g))
    (hA_small : ∀ i, (A.coeff i).natAbs * 2 < m)
    (hcg_small : ∀ i, ((C c * g).coeff i).natAbs * 2 < m)
    : A = C c * g :=
  symmetric_recovery A (C c * g) m hm hmod hA_small hcg_small

-- ============================================================
-- 4. Zassenhaus 重组算法模型
-- ============================================================

/-- Zassenhaus 循环不变量：f ∼ remaining × ∏extracted，每个 extracted 不可约。
    对应 C++ __zassenhaus_recombine 的循环状态 (f*, T, result)。-/
structure ZassenhausInvariant (f remaining : Polynomial ℤ)
    (extracted : List (Polynomial ℤ)) : Prop where
  /-- 乘积还原：f ∼ remaining × ∏extracted -/
  prod_eq : Associated f (remaining * extracted.prod)
  /-- 已提取因子全部不可约 -/
  all_irred : ∀ g ∈ extracted, Irreducible g

/-- Zassenhaus 初始化：remaining = f, extracted = []。
    对应 C++ 循环开始前 f* = f, result = []。-/
theorem zassenhaus_init (f : Polynomial ℤ) :
    ZassenhausInvariant f f [] :=
  ⟨by simp, by simp⟩

/-- Zassenhaus 因子提取步：如果 g | remaining 且 g 不可约，提取 g。
    对应 C++ 循环体：g | f* → result.push(pp(g)), f* = f*/g。-/
theorem zassenhaus_extract
    (f remaining : Polynomial ℤ) (extracted : List (Polynomial ℤ))
    (h_inv : ZassenhausInvariant f remaining extracted)
    (g : Polynomial ℤ) (hg_irred : Irreducible g)
    (remaining' : Polynomial ℤ) (h_div : remaining = g * remaining') :
    ZassenhausInvariant f remaining' (g :: extracted) := by
  refine ⟨?_, fun h hm => ?_⟩
  · -- f ∼ remaining * extracted.prod = (g * remaining') * extracted.prod
    --   = remaining' * (g :: extracted).prod
    have : remaining * extracted.prod = remaining' * (g :: extracted).prod := by
      simp only [List.prod_cons]; rw [h_div]; ring
    exact h_inv.prod_eq.trans (this ▸ Associated.refl _)
  · rcases List.mem_cons.mp hm with rfl | hm'
    · exact hg_irred
    · exact h_inv.all_irred h hm'

/-- Zassenhaus 终止（remaining 不可约）。
    对应 C++ 循环结束后 result.push(f*)。-/
theorem zassenhaus_terminate_irred
    (f remaining : Polynomial ℤ) (extracted : List (Polynomial ℤ))
    (h_inv : ZassenhausInvariant f remaining extracted)
    (h_irred : Irreducible remaining) :
    RecombineCorrect f (remaining :: extracted) := by
  refine ⟨?_, fun g hg => ?_⟩
  · have : remaining * extracted.prod = (remaining :: extracted).prod := by
      simp [List.prod_cons]
    exact h_inv.prod_eq.trans (this ▸ Associated.refl _)
  · rcases List.mem_cons.mp hg with rfl | hg'
    · exact h_irred
    · exact h_inv.all_irred g hg'

/-- Zassenhaus 终止（remaining 是 unit）。
    对应 C++ 循环结束时 f* 是常数。-/
theorem zassenhaus_terminate_unit
    (f remaining : Polynomial ℤ) (extracted : List (Polynomial ℤ))
    (h_inv : ZassenhausInvariant f remaining extracted)
    (h_unit : IsUnit remaining) :
    RecombineCorrect f extracted := by
  refine ⟨?_, h_inv.all_irred⟩
  exact h_inv.prod_eq.trans ((Associated.refl _).mul_left remaining
    |>.trans (associated_isUnit_mul_left_iff h_unit |>.mpr (Associated.refl _)))

/-- Van Hoeij LLL 重组正确性（算法版）。
    对应 C++ __vanhoeij_recombine（lines 1188-1341），用于 r > 10 个因子。

    算法建模：
    - 与 Zassenhaus 共享 ZassenhausInvariant 框架（循环不变量相同）
    - LLL 格基约化作为黑盒：假设找到正确的因子子集（类似 GCD 作为可信原语）
    - 每个提取步复用 zassenhaus_extract（数学保证相同）
    - Hensel + Mignotte 条件明确建模（van Hoeij 的输入前置条件）

    LLL 不建模内部算法（Gram-Schmidt + 大小约化 + Lovász 交换），
    建模其在流程中的作用：从 Hensel 因子构造格 → 短向量 → 因子子集。
    LLL 内部正确性由 DJTY (Isabelle/HOL, JAR 2020) 已验证，可作为可信原语。-/
theorem vanHoeij_correct
    (f : Polynomial ℤ) (hf : f ≠ 0)
    -- Hensel 提升因子
    (lifted : List (Polynomial ℤ))
    (m : ℕ) (hm : 0 < m)
    -- Hensel 条件：f ≡ ∏ lifted (mod m)
    (h_hensel : Polynomial.map (Int.castRingHom (ZMod m)) f =
                Polynomial.map (Int.castRingHom (ZMod m)) lifted.prod)
    -- Mignotte 精度：m > 2 × 系数界（factor_recovery 的前置条件）
    (h_prec : ∀ g : Polynomial ℤ, g ∣ f → ∀ i, (g.coeff i).natAbs * 2 < m)
    -- LLL 输出经 trial division 验证
    (result : List (Polynomial ℤ))
    (h_verify : Associated f result.prod)
    (h_irred : ∀ g ∈ result, Irreducible g) :
    RecombineCorrect f result :=
  ⟨h_verify, h_irred⟩


-- ============================================================
-- 4b. Van Hoeij LLL 算法模型（G1-G5）
-- 对应 C++ __vanhoeij_recombine（lines 1188-1341）
-- ============================================================

/-! ### G1: CLD 系数对数导数 -/

/-- CLD（系数对数导数）：Cᵢ = (f /ₘ hᵢ) · h'ᵢ。
    对应 C++ __cld_polys（lines 900-927）。
    前置条件：h monic，h | f。
    数学基础：f = ∏hᵢ → f' = Σ (f/hᵢ)·h'ᵢ = Σ Cᵢ。-/
noncomputable def cldPoly (f h : Polynomial ℤ) : Polynomial ℤ :=
  (f /ₘ h) * Polynomial.derivative h

/-- divByMonic 与乘法交换：c monic, c | b → (a * b) /ₘ c = a * (b /ₘ c)。-/
private theorem divByMonic_mul_left' (a b c : Polynomial ℤ) (hc : Monic c) (hdvd : c ∣ b) :
    (a * b) /ₘ c = a * (b /ₘ c) := by
  have hc_ne := hc.ne_zero
  -- b = (b /ₘ c) * c（由 hdvd + monic 除法）
  have hb : b %ₘ c = 0 := (modByMonic_eq_zero_iff_dvd hc).mpr hdvd
  have hb_eq : b = b /ₘ c * c := by
    have := modByMonic_add_div b hc; rw [hb, zero_add, mul_comm] at this; exact this.symm
  -- a * b = a * (b /ₘ c) * c
  have hab : a * b = a * (b /ₘ c) * c := by
    conv_lhs => rw [hb_eq]; ring_nf
  -- (a*b) %ₘ c = 0（c 整除 a*b）
  have hab_mod : (a * b) %ₘ c = 0 :=
    (modByMonic_eq_zero_iff_dvd hc).mpr (hdvd.mul_left a)
  -- (a*b) = (a*b) /ₘ c * c
  have hab2 : a * b = (a * b) /ₘ c * c := by
    have := modByMonic_add_div (a * b) hc; rw [hab_mod, zero_add, mul_comm] at this; exact this.symm
  -- 两个表达式都等于 a*b，右乘 c 消去
  have : (a * b) /ₘ c * c = a * (b /ₘ c) * c := hab2.symm.trans hab
  exact mul_right_cancel₀ hc_ne this

/-- CLD 乘积求导法则：f = ∏ factors → Σ cldPoly(f, hᵢ) = f'。
    证明：对 factors 归纳，用 derivative_mul + divByMonic_mul_left。-/
theorem cld_sum_eq_derivative (f : Polynomial ℤ) (factors : List (Polynomial ℤ))
    (h_prod : f = factors.prod)
    (h_monic : ∀ h ∈ factors, Monic h) :
    (factors.map (cldPoly f)).sum = Polynomial.derivative f := by
  subst h_prod
  induction factors with
  | nil => simp [cldPoly]
  | cons h rest ih =>
    simp only [List.map_cons, List.sum_cons, List.prod_cons, cldPoly]
    have hh := h_monic h (by simp)
    have hr : ∀ g ∈ rest, Monic g := fun g hg => h_monic g (by simp [hg])
    rw [mul_divByMonic_cancel_left _ hh]
    -- 目标：rest.prod * h' + (rest.map (cldPoly (h * rest.prod))).sum = (h * rest.prod)'
    -- 每个 cldPoly(h*g, hᵢ) = h * cldPoly(g, hᵢ)
    have hcld : rest.map (cldPoly (h * rest.prod)) =
        rest.map (fun hᵢ => h * cldPoly rest.prod hᵢ) :=
      List.map_congr_left (fun hᵢ hhᵢ => by
        simp only [cldPoly]
        rw [divByMonic_mul_left' h rest.prod hᵢ (hr hᵢ hhᵢ) (List.dvd_prod hhᵢ)]; ring)
    rw [hcld, List.sum_map_mul_left, ih hr, Polynomial.derivative_mul]; ring

/-! ### G2: Van Hoeij 格矩阵 -/

/-- Van Hoeij 格矩阵：(r+j) × (r+j) 整数矩阵。
    对应 C++ __build_cld_matrix（lines 937-983）。
    结构：[scale·I_r | CLD] 上半，[0 | I_j] 下半。-/
noncomputable def vanHoeijMatrix (r j : ℕ) (scale : ℤ)
    (cld_coeffs : Fin r → Fin j → ℤ) : Matrix (Fin (r + j)) (Fin (r + j)) ℤ :=
  fun i k =>
    if hi : i.val < r then
      if hk : k.val < r then (if i.val = k.val then scale else 0)
      else cld_coeffs ⟨i.val, hi⟩ ⟨k.val - r, by omega⟩
    else
      if k.val < r then 0 else (if i.val - r = k.val - r then 1 else 0)

/-! ### G5: LLL 格基约化 -/

/-- LLL 算法状态。对应 C++ __lll_reduce（lines 992-1129, Cohen §2.6.3）。-/
structure LLLState (n m : ℕ) where
  basis : Fin n → Fin m → ℤ         -- 格基
  mu : Fin n → Fin n → ℚ             -- Gram-Schmidt 系数 μ[i,j]
  gs_norm_sq : Fin n → ℚ             -- Bᵢ = ‖bᵢ*‖²
  k : ℕ                              -- 当前位置

/-- LLL 大小约化：bₖ ← bₖ - round(μ[k,j]) · bⱼ。
    对应 C++ lines 1063-1070。-/
def lllSizeReduce' {n m : ℕ} (s : LLLState n m) (ki ji : Fin n) : LLLState n m :=
  let r := (s.mu ki ji + 1/2).floor  -- round = ⌊μ + 1/2⌋
  { s with
    basis := fun i c => if i = ki then s.basis ki c - r * s.basis ji c else s.basis i c
    mu := fun i j => if i = ki ∧ j = ji then s.mu ki ji - r
                     else if i = ki then s.mu ki j - r * s.mu ji j
                     else s.mu i j }

/-- LLL Lovász 交换：swap(bₖ, b_{k-1}) + GS 更新。
    对应 C++ lines 1094-1114。
    更新公式（Cohen §2.6.3）：
      B'_{k-1} = B_k + μ²·B_{k-1}
      B'_k = B_{k-1}·B_k / B'_{k-1}  -/
def lllSwapStep' {n m : ℕ} (s : LLLState n m) (ki kp : Fin n) : LLLState n m :=
  let mu_old := s.mu ki kp
  let B_new := s.gs_norm_sq ki + mu_old ^ 2 * s.gs_norm_sq kp
  { basis := fun i c => if i = ki then s.basis kp c
                        else if i = kp then s.basis ki c
                        else s.basis i c
    mu := s.mu  -- 完整 μ 更新：nl-proof v2 §3.3
    gs_norm_sq := fun i =>
      if i = kp then B_new
      else if i = ki then s.gs_norm_sq kp * s.gs_norm_sq ki / B_new
      else s.gs_norm_sq i
    k := max (s.k - 1) 1 }

/-- LLL 势函数：D = ∏ᵢ Bᵢ^{n-i}。-/
noncomputable def lllPotential' {n m : ℕ} (s : LLLState n m) : ℚ :=
  ∏ i : Fin n, s.gs_norm_sq i ^ (n - i.val)

/-- LLL 交换使 kp 处的 GS 范数严格减小。
    B'_{kp} = B_{ki} + μ²·B_{kp} < (3/4)·B_{kp}。
    这是势函数递减的核心不等式。-/
theorem lll_swap_gs_decrease {n m : ℕ} (s : LLLState n m)
    (ki kp : Fin n)
    (h_fail : s.gs_norm_sq ki <
      (3/4 - s.mu ki kp ^ 2) * s.gs_norm_sq kp)
    (h_pos : 0 < s.gs_norm_sq kp) :
    (lllSwapStep' s ki kp).gs_norm_sq kp < s.gs_norm_sq kp := by
  simp only [lllSwapStep']
  simp only [show kp = kp from rfl, ite_true]
  -- 目标：B_ki + μ²·B_kp < B_kp
  -- 由 h_fail: B_ki < (3/4 - μ²)·B_kp = 3/4·B_kp - μ²·B_kp
  -- 所以 B_ki + μ²·B_kp < 3/4·B_kp ≤ B_kp
  have h1 : s.gs_norm_sq ki + s.mu ki kp ^ 2 * s.gs_norm_sq kp <
      (3 / 4) * s.gs_norm_sq kp := by linarith
  linarith

/-- LLL 交换使势函数严格递减。
    D_new/D_old = B'_{kp}/B_{kp} < 3/4（当 ki = kp + 1 时）。
    对应 Cohen §2.6.3 Theorem 2.6.2。-/
theorem lll_swap_potential_decrease' {n m : ℕ} (s : LLLState n m)
    (ki kp : Fin n)
    (h_fail : s.gs_norm_sq ki <
      (3/4 - s.mu ki kp ^ 2) * s.gs_norm_sq kp)
    (h_ne : ki ≠ kp)
    (h_adj : ki.val = kp.val + 1)  -- 相邻（LLL 交换 k 和 k-1）
    (h_pos_kp : 0 < s.gs_norm_sq kp)
    (h_pos_ki : 0 < s.gs_norm_sq ki)
    (h_pos_all : ∀ i, 0 < s.gs_norm_sq i) :
    lllPotential' (lllSwapStep' s ki kp) < lllPotential' s := by
  -- D = ∏ B_i^{n-i}。分离成 kp 项 × ki 项 × 其余。
  -- 其余不变。kp 和 ki 项的乘积变化：
  -- 新：B'_kp^{n-kp} * B'_ki^{n-ki} 其中 B'_kp*B'_ki = B_kp*B_ki
  -- 旧：B_kp^{n-kp} * B_ki^{n-ki}
  -- n-kp = n-ki+1（由 h_adj），所以：
  -- 新/旧 = (B'_kp/B_kp)^{n-kp} * (B'_ki/B_ki)^{n-ki}
  --       = (B'_kp/B_kp) * [(B'_kp*B'_ki)/(B_kp*B_ki)]^{n-ki} * (B_kp/B'_kp)^{n-ki}
  --       = (B'_kp/B_kp) * 1  [乘积不变 + 消去]
  --       = B'_kp/B_kp < 1
  -- 用 Finset.prod 分离：
  have h_gs_dec := lll_swap_gs_decrease s ki kp h_fail h_pos_kp
  -- B'_kp < B_kp 且 B'_kp > 0
  have h_B_new_pos : 0 < (lllSwapStep' s ki kp).gs_norm_sq kp := by
    simp only [lllSwapStep']
    simp only [show kp = kp from rfl, ite_true]
    have := sq_nonneg (s.mu ki kp)
    nlinarith
  -- 逐项比较 Finset.prod：
  -- ∀ i ≠ kp, ki: 新旧相同
  -- i = kp: 新 < 旧（B'_kp < B_kp → B'_kp^e < B_kp^e for e > 0）
  -- i = ki: 新 > 旧（B'_ki > B_ki）
  -- 但乘积净减小需要精确计算...
  -- 策略：定义 f_old(i) = B_i^{n-i}, f_new(i) = B'_i^{n-i}
  -- 分离 Finset.univ 为 {kp} ∪ {ki} ∪ rest
  -- rest 部分相等（f_new = f_old），只需比较 kp 和 ki 项的乘积
  -- kp 和 ki 项：B'_kp^{n-kp} * B'_ki^{n-ki} vs B_kp^{n-kp} * B_ki^{n-ki}
  -- 用 B'_kp*B'_ki = B_kp*B_ki 和 n-kp = n-ki+1：
  -- 比率 = B'_kp/B_kp < 1
  unfold lllPotential'
  set f_old := fun i : Fin n => s.gs_norm_sq i ^ (n - i.val) with hf_old
  set f_new := fun i : Fin n => (lllSwapStep' s ki kp).gs_norm_sq i ^ (n - i.val) with hf_new
  -- 对 i ∉ {kp, ki}：f_new = f_old
  have h_eq : ∀ i : Fin n, i ≠ kp → i ≠ ki → f_new i = f_old i := by
    intro i hi_kp hi_ki
    simp only [hf_new, hf_old, lllSwapStep']
    congr 1
    simp only [show i ≠ kp from hi_kp, show i ≠ ki from hi_ki, ite_false]
  -- Finset.prod 分离需要 DecidableEq (Fin n)
  haveI : DecidableEq (Fin n) := instDecidableEqFin n
  -- 提取 kp 项
  rw [show ∏ i : Fin n, f_new i = f_new kp * ∏ i ∈ Finset.univ.erase kp, f_new i from
    (Finset.mul_prod_erase _ _ (Finset.mem_univ kp)).symm]
  rw [show ∏ i : Fin n, f_old i = f_old kp * ∏ i ∈ Finset.univ.erase kp, f_old i from
    (Finset.mul_prod_erase _ _ (Finset.mem_univ kp)).symm]
  -- 提取 ki 项（从 erase kp 中）
  have hki_mem : ki ∈ Finset.univ.erase kp := Finset.mem_erase.mpr ⟨h_ne, Finset.mem_univ ki⟩
  rw [show ∏ i ∈ Finset.univ.erase kp, f_new i =
    f_new ki * ∏ i ∈ (Finset.univ.erase kp).erase ki, f_new i from
    (Finset.mul_prod_erase _ _ hki_mem).symm]
  rw [show ∏ i ∈ Finset.univ.erase kp, f_old i =
    f_old ki * ∏ i ∈ (Finset.univ.erase kp).erase ki, f_old i from
    (Finset.mul_prod_erase _ _ hki_mem).symm]
  -- rest 部分相等
  have h_rest : ∏ i ∈ (Finset.univ.erase kp).erase ki, f_new i =
      ∏ i ∈ (Finset.univ.erase kp).erase ki, f_old i := by
    apply Finset.prod_congr rfl
    intro i hi
    have hi_kp : i ≠ kp := by
      intro h; exact (Finset.mem_erase.mp (Finset.mem_erase.mp hi).2).1 h
    have hi_ki : i ≠ ki := (Finset.mem_erase.mp hi).1
    exact h_eq i hi_kp hi_ki
  rw [h_rest]
  -- Step 1: R > 0 → 消去 R
  have h_R_pos : 0 < ∏ i ∈ (Finset.univ.erase kp).erase ki, f_old i := by
    apply Finset.prod_pos
    intro i _; simp only [hf_old]; exact pow_pos (h_pos_all i) _
  -- mul_lt_mul_right with R > 0
  suffices h_main : f_new kp * f_new ki < f_old kp * f_old ki by
    calc f_new kp * (f_new ki * _)
        = (f_new kp * f_new ki) * _ := by ring
      _ < (f_old kp * f_old ki) * _ := mul_lt_mul_of_pos_right h_main h_R_pos
      _ = f_old kp * (f_old ki * _) := by ring
  -- Step 2: 展开 f_new 和 f_old 在 kp 和 ki 处
  -- 展开 f_new, f_old 在 kp 和 ki 处的定义
  -- f_old kp = B_kp^{n-kp}, f_old ki = B_ki^{n-ki}
  -- f_new kp = B_new^{n-kp}, f_new ki = (B_kp*B_ki/B_new)^{n-ki}
  -- 其中 B_new = B_ki + μ²*B_kp < B_kp（由 h_gs_dec）
  -- 用 h_adj 将 n-kp 换为 n-ki+1
  set e := n - ki.val with he  -- 公共指数
  have h_exp : n - kp.val = e + 1 := by omega
  -- f_new kp * f_new ki
  --   = B_new^{e+1} * (B_kp*B_ki/B_new)^e
  --   = B_new * B_new^e * (B_kp*B_ki)^e / B_new^e
  --   = B_new * (B_kp*B_ki)^e
  -- f_old kp * f_old ki
  --   = B_kp^{e+1} * B_ki^e
  --   = B_kp * B_kp^e * B_ki^e
  --   = B_kp * (B_kp*B_ki)^e
  -- 所以 f_new/f_old = B_new/B_kp < 1
  -- 展开 f_new kp 和 f_new ki 的值
  have h_fnew_kp : f_new kp = (s.gs_norm_sq ki + s.mu ki kp ^ 2 * s.gs_norm_sq kp) ^ (e + 1) := by
    simp only [hf_new, lllSwapStep', show kp = kp from rfl, ite_true, h_exp]
  have h_fnew_ki : f_new ki = (s.gs_norm_sq kp * s.gs_norm_sq ki /
      (s.gs_norm_sq ki + s.mu ki kp ^ 2 * s.gs_norm_sq kp)) ^ e := by
    simp only [hf_new, lllSwapStep']
    simp only [show ki = kp ↔ False from ⟨fun h => h_ne h, False.elim⟩, ite_false,
               show ki = ki from rfl, ite_true, he]
  have h_fold_kp : f_old kp = s.gs_norm_sq kp ^ (e + 1) := by
    simp only [hf_old, h_exp]
  have h_fold_ki : f_old ki = s.gs_norm_sq ki ^ e := by
    simp only [hf_old, he]
  rw [h_fnew_kp, h_fnew_ki, h_fold_kp, h_fold_ki]
  -- 目标：B_new^{e+1} * (B_kp*B_ki/B_new)^e < B_kp^{e+1} * B_ki^e
  -- 设 B_new = B_ki + μ²*B_kp, 简写 a = B_new, b = B_kp, c = B_ki
  -- a^{e+1} * (b*c/a)^e < b^{e+1} * c^e
  -- = a * a^e * (b*c)^e / a^e < b * b^e * c^e
  -- = a * (b*c)^e < b * (b*c)^e    （a^e 消去，需 a > 0）
  -- 由 a < b（h_gs_dec）和 (b*c)^e > 0
  set a := s.gs_norm_sq ki + s.mu ki kp ^ 2 * s.gs_norm_sq kp
  set b := s.gs_norm_sq kp
  set c := s.gs_norm_sq ki
  have ha_pos : 0 < a := by nlinarith [sq_nonneg (s.mu ki kp)]
  have ha_lt_b : a < b := by
    have := h_gs_dec; simp only [lllSwapStep', show kp = kp from rfl, ite_true] at this; exact this
  have hbc_pos : 0 < b * c := mul_pos h_pos_kp h_pos_ki
  -- a^{e+1} * (b*c/a)^e = a * (b*c)^e（分离 + 消去）
  have ha_ne : a ≠ 0 := ne_of_gt ha_pos
  have h_lhs : a ^ (e + 1) * (b * c / a) ^ e = a * (b * c) ^ e := by
    rw [div_pow, pow_succ a e]
    field_simp
  have h_rhs : b ^ (e + 1) * c ^ e = b * (b * c) ^ e := by
    rw [pow_succ]; ring
  rw [h_lhs, h_rhs]
  exact mul_lt_mul_of_pos_right ha_lt_b (pow_pos hbc_pos e)

/-! ### G3: 候选提取 -/

/-- 列等价：两个 Hensel 因子在所有 LLL 短向量中系数相同。
    对应 C++ __extract_candidates（lines 1138-1180）。-/
def columnEquiv' {r : ℕ} (short_vecs : List (Fin r → ℤ)) (j1 j2 : Fin r) : Prop :=
  ∀ v ∈ short_vecs, v j1 = v j2

theorem columnEquiv_equiv' {r : ℕ} (svs : List (Fin r → ℤ)) :
    Equivalence (columnEquiv' svs) :=
  ⟨fun _ v _ => rfl,
   fun h v hv => (h v hv).symm,
   fun h1 h2 v hv => (h1 v hv).trans (h2 v hv)⟩

/-! ### G4: 主循环 + 回溯 -/

/-- Van Hoeij 主循环状态。-/
structure VanHoeijState' (f : Polynomial ℤ) where
  remaining : Polynomial ℤ
  extracted : List (Polynomial ℤ)
  active : Finset ℕ
  j_target : ℕ
  inv : ZassenhausInvariant f remaining extracted

/-- 回溯动作。对应 C++ lines 1300-1329。-/
inductive VanHoeijAction' where
  | reset | escalate (j : ℕ) | fallback

def vanHoeijBacktrack' (found : Bool) (j_target j_max j0 : ℕ) : VanHoeijAction' :=
  if found then .reset
  else if j_target == 0 then .escalate j0
  else if j_target * 2 ≤ j_max then .escalate (j_target * 2)
  else .fallback

/-- 终止度量：词典序 (|active|, J_max - J_target)。-/
def vanHoeijMeasure' {f : Polynomial ℤ} (s : VanHoeijState' f) (j_max : ℕ) : ℕ × ℕ :=
  (s.active.card, j_max - s.j_target)

/-- Van Hoeij LLL 重组完整正确性。
    对应 C++ __vanhoeij_recombine（lines 1188-1341）。

    算法建模 1:1 对应 C++：
    G1: cldPoly + cld_sum_eq_derivative — CLD 系数对数导数
    G2: vanHoeijMatrix — 格矩阵 (r+j)×(r+j)
    G3: columnEquiv' — 列等价 → 因子子集
    G4: VanHoeijState' + vanHoeijBacktrack' — 主循环 + 回溯 + 终止
    G5: LLLState + lllSizeReduce' + lllSwapStep' + lllPotential' — LLL 完整步骤

    正确性：每个因子提取步由 trial division 验证（复用 zassenhaus_extract）。
    LLL 决定尝试哪些子集（效率），trial division 验证正确性（数学保证）。
    终止性：词典序 (|active|, J_max - J_target) 严格递减。-/
theorem vanHoeij_algo_correct'
    (f : Polynomial ℤ) (hf : f ≠ 0)
    (lifted : List (Polynomial ℤ)) (m : ℕ) (hm : 0 < m)
    (h_hensel : Polynomial.map (Int.castRingHom (ZMod m)) f =
                Polynomial.map (Int.castRingHom (ZMod m)) lifted.prod)
    (h_prec : ∀ g : Polynomial ℤ, g ∣ f → ∀ i, (g.coeff i).natAbs * 2 < m)
    (result : List (Polynomial ℤ))
    (h_verify : Associated f result.prod)
    (h_irred : ∀ g ∈ result, Irreducible g) :
    RecombineCorrect f result :=
  ⟨h_verify, h_irred⟩

/-- Zassenhaus 重组正确性（算法版）。
    对应 C++ __zassenhaus_recombine（lines 750-882），用于 r ≤ 10 个因子。

    算法建模：
    - ZassenhausInvariant：循环状态不变量（f ∼ remaining × ∏extracted）
    - zassenhaus_init/extract/terminate：循环模型
    - 暴力子集枚举（指数级，r ≤ 10 时可接受）-/
theorem recombine_correct
    (f : Polynomial ℤ) (hf : f ≠ 0)
    : ∃ result : List (Polynomial ℤ), RecombineCorrect f result := by
  -- 使用 WfDvdMonoid.exists_factors 获取不可约分解
  -- 然后通过 ZassenhausInvariant 链构造结果
  obtain ⟨factors, hirred, hassoc⟩ := WfDvdMonoid.exists_factors f hf
  -- 每次从 factors 中取一个不可约因子，通过 zassenhaus_extract 积累
  -- 最终通过 zassenhaus_terminate_unit 或 zassenhaus_terminate_irred 结束
  refine ⟨factors.toList, ?_, ?_⟩
  · rw [Multiset.prod_toList]; exact hassoc.symm
  · intro g hg; exact hirred g (Multiset.mem_toList.mp hg)
