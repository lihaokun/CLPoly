/-
  CLPoly/Pipeline/FactorZZInstantiate.lean — 端到端 Z[x] 因式分解定理

  实例化 Pipeline 框架（FactorZZ.lean）的参数化证明，
  用 L2 算法实现替代抽象子过程。
-/
import CLPoly.Pipeline.FactorZZ
import CLPoly.Pipeline.FactorZpInstantiate
import CLPoly.Algorithm.Hensel
import CLPoly.Algorithm.Recombine

set_option autoImplicit false

open Polynomial

/-- 端到端 Z[x] 因式分解：任何非零整数多项式有不可约分解。
    RecombineCorrect ≡ FactorZZCorrect（定义相同），
    recombine_correct 直接给出 Z[x] UFD 分解。 -/
theorem factor_ZZ_instantiate
    (f : Polynomial ℤ) (hf : f ≠ 0)
    : ∃ result : List (Polynomial ℤ), FactorZZCorrect f result := by
  -- RecombineCorrect ≡ FactorZZCorrect (definitionally equal)
  exact recombine_correct f hf
