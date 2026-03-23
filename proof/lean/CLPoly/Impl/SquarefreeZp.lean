/-
  CLPoly/Impl/SquarefreeZp.lean — L1 SQF 实现模型 + 精化证明

  1:1 对应 C++: polynomial_factorize_zp.hh:122-180 __squarefree_Zp
  精化目标: squarefree_impl.toPoly = sqfZp input.toPoly

  假设基本操作正确（作为参数）：
  - gcd_impl = EuclideanDomain.gcd
  - div_impl = divByMonic
  - derivative_impl = derivative
  - extractPthRoot_impl = contract p
-/
import CLPoly.Impl.Types
import CLPoly.Algorithm.SquarefreeZp

set_option autoImplicit false
set_option maxHeartbeats 800000

open Polynomial

variable {p : ℕ} [hp : Fact (Nat.Prime p)]

-- ============================================================
-- 1. 基本操作接口（假设正确，作为参数传入）
-- ============================================================

/-- 基本多项式操作的正确性假设。
    L1 SQF 证明假设这些操作正确，不验证其内部实现。 -/
-- 注：SparsePolyZp.toPoly 返回 Polynomial (ZMod sp.prime)，
-- sp.prime 是结构体字段（Nat），不是 variable p。
-- 所以 PolyOpsCorrect 中的假设直接使用 sp.prime，不需要外部 p。

/-- 基本多项式操作的正确性假设。
    L1 SQF 证明假设这些操作正确，不验证其内部实现。
    所有操作保持 .prime 字段不变。 -/
structure PolyOpsCorrect where
  -- GCD 实现
  gcd_impl : SparsePolyZp → SparsePolyZp → SparsePolyZp
  -- 除法实现
  div_impl : SparsePolyZp → SparsePolyZp → SparsePolyZp
  -- 导数实现
  derivative_impl : SparsePolyZp → SparsePolyZp
  -- p 次根提取实现
  extractPthRoot_impl : SparsePolyZp → SparsePolyZp
  -- 首一化实现
  makeMonic_impl : SparsePolyZp → SparsePolyZp
  -- 归一化实现
  normalize_impl : SparsePolyZp → SparsePolyZp
  -- 素数保持
  h_prime_preserved : ∀ op ∈ [gcd_impl, div_impl, derivative_impl,
    extractPthRoot_impl, makeMonic_impl, normalize_impl],
    ∀ a b : SparsePolyZp, True  -- placeholder, 具体保持条件在各假设中

-- ============================================================
-- 2. L1 SQF 实现（1:1 对应 C++ __squarefree_Zp）
-- ============================================================

-- TODO: 定义 squarefree_impl 函数，精确匹配 C++ 控制流
-- TODO: 精化定理 squarefree_impl.map toPoly = sqfZp input.toPoly

-- 此文件的完整实现需要:
-- 1. yunLoop_impl（对应 C++ while 循环）
-- 2. squarefree_impl（对应 C++ __squarefree_Zp）
-- 3. refine_squarefree（精化定理）

-- 每步精化通过 PolyOpsCorrect 的假设，
-- 将 L1 操作（gcd_impl, div_impl 等）替换为 L2 操作（EuclideanDomain.gcd, /ₘ 等），
-- 然后引用已证的 L2 sqf_correct。

end
