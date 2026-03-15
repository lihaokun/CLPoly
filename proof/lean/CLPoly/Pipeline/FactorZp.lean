/-
  CLPoly/Pipeline/FactorZp.lean — Zp[x] 顶层正确性骨架

  Phase 1：陈述顶层定理 + sorry，锁定接口。
  定理含义：若 SQF、DDF、EDF 三个子过程各自正确，
  则它们的组合产生完整的 Zp[x] 不可约分解。

  证明策略（Phase 3 填充）：
  1. SQF 将 f 分解为 ∏ gᵢ^eᵢ，每个 gᵢ 首一无平方
  2. 对每个 gᵢ，DDF 产出 [(gd_j, d_j)]，按度分组
  3. 对每组 (gd_j, d_j)，EDF 拆出不可约因子
  4. 收集所有因子，附上重数 eᵢ 和 leading coefficient
-/
import CLPoly.Spec

set_option autoImplicit false

open Polynomial

variable {p : ℕ} [Fact (Nat.Prime p)]

/-- Zp[x] 因式分解的顶层正确性：
    假设 SQF、DDF、EDF 各自正确，则组合结果是完整因式分解 -/
theorem factor_Zp_correct
    (f : Polynomial (ZMod p)) (hf : f ≠ 0)
    -- 假设各子过程存在且正确
    (sqf : Polynomial (ZMod p) → List (Polynomial (ZMod p) × ℕ))
    (hsqf : ∀ g, g ≠ 0 → SquarefreeDecomp g (sqf g))
    (ddf : Polynomial (ZMod p) → List (Polynomial (ZMod p) × ℕ))
    (hddf : ∀ g, Monic g → Squarefree g → DDFCorrect g (ddf g))
    (edf : Polynomial (ZMod p) → ℕ → List (Polynomial (ZMod p)))
    (hedf : ∀ g d, Monic g → Squarefree g →
            (∀ q, Irreducible q → q ∣ g → Polynomial.natDegree q = d) →
            EDFCorrect g d (edf g d))
    : ∃ (lc : ZMod p) (factors : List (Polynomial (ZMod p) × ℕ)),
        FactorZpCorrect f lc factors := by
  sorry
