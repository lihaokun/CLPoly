/-
  E2: 定理 2.1 桥接路径验证
  Phase 0 实验 — 确认从 Mathlib 的有限域理论到 DDF 正确性定理的桥接可行性

  定理 2.1：有限域上 squarefree 多项式可以分解为不同次数不可约因子的乘积
  X^{p^d} - X = ∏_{k|d} (所有 k 次不可约多项式的乘积)
-/
import Mathlib.FieldTheory.Finite.Basic
import Mathlib.FieldTheory.Finite.GaloisField
import Mathlib.Algebra.Squarefree.Basic
import Mathlib.FieldTheory.Separable
import Mathlib.RingTheory.Polynomial.Basic
import Mathlib.Algebra.Polynomial.FieldDivision
import Mathlib.RingTheory.UniqueFactorizationDomain.NormalizedFactors

open Polynomial FiniteField

noncomputable section

variable (p : ℕ) [hp : Fact (Nat.Prime p)]

-- ============================================================
-- 步骤 1：确认 roots 定理的精确签名
-- ============================================================

#check @roots_X_pow_card_sub_X
-- 预期：X^{|F|} - X 的根 = F 中所有元素

-- ============================================================
-- 步骤 2：X^{p^d} - X 的 Separable 性质
-- ============================================================

-- 思路：derivative(X^{p^d} - X) = p^d · X^{p^d-1} - 1 = -1（char p 下 p^d = 0）
-- gcd(X^{p^d} - X, -1) = 1 → Separable

-- 先检查 derivative 相关引理
#check @Polynomial.derivative_sub
#check @Polynomial.derivative_X_pow
#check @Polynomial.derivative_X
#check @CharP.cast_eq_zero

-- 尝试证明 X^{p^d} - X 是 Separable 的
example (d : ℕ) (hd : 0 < d) :
    Separable (X ^ (p ^ d) - X : Polynomial (ZMod p)) := by
  sorry -- TODO: 尝试 simp + ring 或手动展开

-- ============================================================
-- 步骤 3：UFD 中 Squarefree 与 normalizedFactors 的关系
-- ============================================================

-- 在 UFD 中，squarefree ↔ normalizedFactors 无重复
#check @UniqueFactorizationMonoid.squarefree_iff_nodup_normalizedFactors
-- 完整名称是 UniqueFactorizationMonoid.squarefree_iff_nodup_normalizedFactors

-- Polynomial (ZMod p) 是 UFD
example : UniqueFactorizationMonoid (Polynomial (ZMod p)) := inferInstance

-- ============================================================
-- 步骤 4：不可约多项式与有限域的关系
-- ============================================================

-- GaloisField 存在且有唯一的 p^n 阶有限域
#check @GaloisField
#check @GaloisField.card

-- 不可约多项式 g ∈ (ZMod p)[X], deg g = k
-- X^{|K|} - X 在 K 上完全分裂
#check @FiniteField.splits_X_pow_card_sub_X

-- K 是 X^{p^n} - X 的分裂域
#check @FiniteField.isSplittingField_sub

-- 域嵌入：finrank K | finrank L → ∃ K →ₐ[F] L
#check @nonempty_algHom_of_finrank_dvd
-- 双向条件
#check @nonempty_algHom_iff_finrank_dvd

-- GaloisField p n 的 finrank = n
#check @GaloisField.finrank

-- ============================================================
-- 步骤 5：整除性谓词方法（替代 normalizedFactors）
-- ============================================================

-- 我们的 DDF 正确性可以用纯整除谓词表达：
-- ∀ g, Irreducible g → natDeg g = d → g ∣ gd
-- 不需要枚举所有不可约因子，只需整除关系

-- 这个方向需要的关键引理：
-- Irreducible g ∧ deg g = d → g ∣ (X^{p^d} - X)

-- 检查是否有现成引理
-- 核心定理：X^{p^d} - X 整除关系
-- Polynomial.dvd_sub_mod 可能不存在，用 sorry 占位
-- 关键：Irreducible g ∧ deg g = d → g ∣ (X^{p^d} - X)
-- 这需要：g 的分裂域嵌入 F_{p^d} + 根的刻画

-- ============================================================
-- 步骤 6：Frobenius 与幂映射
-- ============================================================

#check @FiniteField.expand_card  -- X^{p^n} = expand p^n X?
-- 或 Frobenius 的多项式版本
#check @Polynomial.expand

end
