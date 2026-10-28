/-
  CLPoly/Pipeline/FactorZpInstantiate.lean — 端到端 Zp 因式分解定理

  只保留显式 splits_fn 及其逐步正确性假设的参数化 L2 定理。
  曾经使用 `Classical.choose`/UFD witness 冒充 EDF 实现的无条件版本已删除。
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
    (sqf_correct f hf)
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
