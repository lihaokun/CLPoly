/-
  CLPoly/Pipeline/FactorZpInstantiate.lean — 端到端 Zp 因式分解定理

  两个版本：
  1. factor_Zp_instantiate：需要 splits_fn（EDF 去随机化参数）
  2. factor_Zp_instantiate_unconditional：无条件版（EDF 用 UFD 替代）
-/
import CLPoly.Pipeline.FactorZp
import CLPoly.Algorithm.SquarefreeZp
import CLPoly.Algorithm.DDF
import CLPoly.Algorithm.EDF

set_option autoImplicit false
set_option maxHeartbeats 800000

open Polynomial

variable {p : ℕ} [hp : Fact (Nat.Prime p)]

/-- 端到端 Zp[x] 因式分解（带 splits_fn 参数版）。 -/
theorem factor_Zp_instantiate
    (f : Polynomial (ZMod p)) (hf : f ≠ 0)
    (splits_fn : Polynomial (ZMod p) → ℕ → List (Polynomial (ZMod p)))
    (hsplits : ∀ (g : Polynomial (ZMod p)) (d : ℕ),
        Monic g → Squarefree g →
        (∀ q, Irreducible q → q ∣ g → q.natDegree = d) →
        ∀ q ∈ edf g d (splits_fn g d), q.natDegree ≤ d)
    : ∃ (lc : ZMod p) (factors : List (Polynomial (ZMod p) × ℕ)),
        FactorZpCorrect f lc factors := by
  exact factor_Zp_correct f hf
    sqfZp
    (fun g hg => sqf_correct g hg)
    ddf
    (fun g hm hsq => ddf_correct g hm hsq)
    (fun g d => if g.natDegree = 0 then [] else edf g d (splits_fn g d))
    (fun g d hm hsq hfactors => by
      simp only
      split
      · rename_i hg_deg
        constructor
        · simp only [List.prod_nil]
          have hg_eq_one : g = 1 := by
            have hg_eq := Polynomial.eq_C_of_natDegree_eq_zero hg_deg
            have : g.coeff 0 = 1 := by
              have := hm.leadingCoeff; rw [Polynomial.leadingCoeff, hg_deg] at this; exact this
            rw [hg_eq, this, map_one]
          rw [hg_eq_one]
        · intro q hq; simp at hq
      · rename_i hg_deg
        have hg_pos : 0 < g.natDegree := Nat.pos_of_ne_zero hg_deg
        have hg_ne : g ≠ 0 := Monic.ne_zero hm
        have hg_nu : ¬IsUnit g := fun hu => absurd (natDegree_eq_zero_of_isUnit hu) hg_deg
        obtain ⟨q₀, hq₀_irr, hq₀_dvd⟩ := WfDvdMonoid.exists_irreducible_factor hg_nu hg_ne
        have hd_pos : 0 < d := by
          have := hfactors q₀ hq₀_irr hq₀_dvd; rw [← this]
          exact Irreducible.natDegree_pos hq₀_irr
        exact edf_correct g d hm hsq hg_pos hd_pos hfactors
          (splits_fn g d) (hsplits g d hm hsq hfactors))

/-- 端到端 Zp[x] 因式分解（无条件版）。
    不依赖外部 splits_fn 假设。EDF 步骤由 F_p[x] UFD 直接提供。 -/
theorem factor_Zp_instantiate_unconditional
    (f : Polynomial (ZMod p)) (hf : f ≠ 0)
    : ∃ (lc : ZMod p) (factors : List (Polynomial (ZMod p) × ℕ)),
        FactorZpCorrect f lc factors := by
  classical
  -- EDF function: use edf_correct_unconditional via Classical.choose
  -- (bypasses Cantor-Zassenhaus algorithm; uses UFD normalizedFactors internally)
  let edf_unconditional : Polynomial (ZMod p) → ℕ → List (Polynomial (ZMod p)) :=
    fun g d =>
      if hg_deg : g.natDegree = 0 then []
      else
        -- All preconditions are decidable via Classical
        if hpre : g.Monic ∧ Squarefree g ∧ 0 < d ∧
            (∀ q : Polynomial (ZMod p), Irreducible q → q ∣ g → q.natDegree = d) then
          (edf_correct_unconditional g d hpre.1 hpre.2.1
            (Nat.pos_of_ne_zero hg_deg) hpre.2.2.1 hpre.2.2.2).choose
        else [g]
  exact factor_Zp_correct f hf
    sqfZp (fun g hg => sqf_correct g hg)
    ddf (fun g hm hsq => ddf_correct g hm hsq)
    edf_unconditional
    (fun g d hm hsq hfactors => by
      simp only [edf_unconditional]
      split
      · -- natDegree = 0: unit
        rename_i hg_deg
        constructor
        · simp only [List.prod_nil]
          have hg_eq_one : g = 1 := by
            have hg_eq := Polynomial.eq_C_of_natDegree_eq_zero hg_deg
            have : g.coeff 0 = 1 := by
              have := hm.leadingCoeff; rw [Polynomial.leadingCoeff, hg_deg] at this; exact this
            rw [hg_eq, this, map_one]
          rw [hg_eq_one]
        · intro q hq; simp at hq
      · -- natDegree > 0
        rename_i hg_deg
        have hg_pos := Nat.pos_of_ne_zero hg_deg
        have hg_nu : ¬IsUnit g := fun hu => absurd (natDegree_eq_zero_of_isUnit hu) hg_deg
        obtain ⟨q₀, hq₀_irr, hq₀_dvd⟩ :=
          WfDvdMonoid.exists_irreducible_factor hg_nu (Monic.ne_zero hm)
        have hd_pos : 0 < d := by
          rw [← hfactors q₀ hq₀_irr hq₀_dvd]; exact Irreducible.natDegree_pos hq₀_irr
        rw [dif_pos ⟨hm, hsq, hd_pos, hfactors⟩]
        exact (edf_correct_unconditional g d hm hsq hg_pos hd_pos hfactors).choose_spec)
