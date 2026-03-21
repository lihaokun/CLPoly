/-
  CLPoly/Pipeline/FactorZpInstantiate.lean — 端到端 Zp 因式分解定理

  实例化 Pipeline 框架（FactorZp.lean）的参数化证明，
  用 L2 算法实现（sqfZp, ddf, edf）替代抽象子过程。

  结果：FactorZpCorrect f lc factors（完整不可约分解）
  条件：f ≠ 0 + EDF splits 充分
-/
import CLPoly.Pipeline.FactorZp
import CLPoly.Algorithm.SquarefreeZp
import CLPoly.Algorithm.DDF
import CLPoly.Algorithm.EDF

set_option autoImplicit false
set_option maxHeartbeats 800000

open Polynomial

variable {p : ℕ} [hp : Fact (Nat.Prime p)]

/-- 端到端 Zp[x] 因式分解：用 L2 算法实例化 Pipeline 框架。
    sqfZp + ddf + edf 组合满足 FactorZpCorrect。

    唯一的外部条件：EDF 的 splits 函数对每个等度多项式提供充分的随机元素。
    这是概率性条件——在实际运行中以压倒性概率成立。 -/
theorem factor_Zp_instantiate
    (f : Polynomial (ZMod p)) (hf : f ≠ 0)
    -- EDF 的 splits 提供者：对每个多项式和度数，给出一个随机元素列表
    (splits_fn : Polynomial (ZMod p) → ℕ → List (Polynomial (ZMod p)))
    -- splits 充分性条件：对每个 DDF 输出的等度多项式，edf 结果的每个元素 natDegree ≤ d
    (hsplits : ∀ (g : Polynomial (ZMod p)) (d : ℕ),
        Monic g → Squarefree g →
        (∀ q, Irreducible q → q ∣ g → q.natDegree = d) →
        ∀ q ∈ edf g d (splits_fn g d), q.natDegree ≤ d)
    : ∃ (lc : ZMod p) (factors : List (Polynomial (ZMod p) × ℕ)),
        FactorZpCorrect f lc factors := by
  exact factor_Zp_correct f hf
    -- sqf = sqfZp
    sqfZp
    (fun g hg => sqf_correct g hg)
    -- ddf = ddf
    ddf
    (fun g hm hsq => ddf_correct g hm hsq)
    -- edf: wrap to handle natDegree = 0 case (return [] for units)
    (fun g d => if g.natDegree = 0 then [] else edf g d (splits_fn g d))
    (fun g d hm hsq hfactors => by
      simp only
      split
      · -- g.natDegree = 0: g is unit (monic + deg 0 = 1). No irred factors. Return [].
        rename_i hg_deg
        constructor
        · -- Associated g [].prod = Associated g 1. g = 1 (monic + deg 0).
          simp only [List.prod_nil]
          have hg_eq_one : g = 1 := by
            have hg_eq := Polynomial.eq_C_of_natDegree_eq_zero hg_deg
            have : g.coeff 0 = 1 := by
              have := hm.leadingCoeff; rw [Polynomial.leadingCoeff, hg_deg] at this; exact this
            rw [hg_eq, this, map_one]
          rw [hg_eq_one]
        · intro q hq; simp at hq
      · -- g.natDegree > 0
        rename_i hg_deg
        have hg_pos : 0 < g.natDegree := Nat.pos_of_ne_zero hg_deg
        have hg_ne : g ≠ 0 := Monic.ne_zero hm
        have hg_nu : ¬IsUnit g := fun hu => absurd (natDegree_eq_zero_of_isUnit hu) hg_deg
        obtain ⟨q₀, hq₀_irr, hq₀_dvd⟩ := WfDvdMonoid.exists_irreducible_factor hg_nu hg_ne
        have hd_pos : 0 < d := by
          have := hfactors q₀ hq₀_irr hq₀_dvd; rw [← this]
          exact Irreducible.natDegree_pos hq₀_irr
        exact edf_correct g d hm hsq hg_pos hd_pos hfactors
          (splits_fn g d) (hsplits g d hm hsq hfactors))
