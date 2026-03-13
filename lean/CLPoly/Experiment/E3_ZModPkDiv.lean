/-
  E3: ZMod (p^k) 上的多项式除法
  Phase 0 实验 — 确认 Hensel 提升所需的基础设施

  关键问题：divByMonic 在非整环 ZMod (p^k) 上是否正常工作？
  castHom 层间投影是否可用？
-/
import Mathlib.RingTheory.Polynomial.Basic
import Mathlib.Algebra.Polynomial.Div
import Mathlib.Data.ZMod.Basic
import Mathlib.Data.ZMod.ValMinAbs

open Polynomial

-- ============================================================
-- 1. ZMod (p^k) 的代数结构
-- ============================================================

-- CommRing（不是 IsDomain）
example : CommRing (ZMod (7 ^ 3)) := inferInstance
example : CommRing (ZMod 343) := inferInstance

-- 确认不是 IsDomain（取消注释应该报错）
-- example : IsDomain (ZMod (7 ^ 3)) := inferInstance

-- ============================================================
-- 2. 首一多项式除法（divByMonic / modByMonic）
-- ============================================================

noncomputable section

-- 构造 ZMod 343 上的多项式
def f343 : Polynomial (ZMod 343) := X ^ 3 + C 2 * X + C 1
def g343 : Polynomial (ZMod 343) := X ^ 2 + C 3  -- 首一

-- divByMonic 应该工作（只需 Ring，不需 IsDomain）
def q343 := f343 /ₘ g343
def r343 := f343 %ₘ g343

-- 恒等式：p = q * d + r 当 d 为 Monic
-- 先验证 g343 是首一的
example : Monic g343 := by
  unfold g343
  exact monic_X_pow_add_C 3 (by norm_num)

-- 除法恒等式：modByMonic_add_div 给出 f %ₘ g + g * (f /ₘ g) = f
example : Monic g343 → f343 %ₘ g343 + g343 * (f343 /ₘ g343) = f343 := by
  intro hm
  exact modByMonic_add_div f343 hm

-- ============================================================
-- 3. castHom 层间投影
-- ============================================================

-- ZMod 343 →+* ZMod 7（343 = 7^3, 7 | 343）
def proj : ZMod 343 →+* ZMod 7 :=
  ZMod.castHom (show 7 ∣ 343 by norm_num) (ZMod 7)

-- Polynomial.map 配合 castHom
def f7_from_343 := Polynomial.map proj f343

-- ZMod 343 →+* ZMod 49（49 | 343）
def proj2 : ZMod 343 →+* ZMod 49 :=
  ZMod.castHom (show 49 ∣ 343 by norm_num) (ZMod 49)

-- ============================================================
-- 4. 单位与可逆性
-- ============================================================

-- 在 ZMod 343 中，与 343 互素的元素是单位
#check @ZMod.unitOfCoprime

-- 3 和 343 互素 → 3 是单位
example : IsUnit (3 : ZMod 343) := by
  sorry -- 尝试 ZMod.unitOfCoprime 或 decide

-- ============================================================
-- 5. 自然提升：ℤ → ZMod n
-- ============================================================

#check @ZMod.intCast_zmod_eq_zero_iff_dvd
#check @ZMod.val
#check @ZMod.valMinAbs  -- 对称模表示

-- ZMod.valMinAbs 的签名
-- 预期：返回 ℤ，范围在 [-n/2, n/2]
example : ZMod.valMinAbs (5 : ZMod 7) = -2 := by decide

-- ============================================================
-- 6. Polynomial.map 保持 Monic
-- ============================================================

#check @Polynomial.Monic.map

-- map proj g343 应该仍然是 Monic
example : Monic (Polynomial.map proj g343) := by
  exact Monic.map proj (by unfold g343; exact monic_X_pow_add_C 3 (by norm_num))

end
