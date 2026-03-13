/-
  E1: Mathlib Polynomial (ZMod p) API 可用性验证
  Phase 0 实验 — 确认关键 typeclass 实例和基本操作
-/
import Mathlib.RingTheory.Polynomial.Basic
import Mathlib.Algebra.Polynomial.FieldDivision
import Mathlib.Data.ZMod.Basic
import Mathlib.FieldTheory.Finite.Basic
import Mathlib.FieldTheory.Separable
import Mathlib.Algebra.Squarefree.Basic

open Polynomial

-- ============================================================
-- 1. Typeclass 实例验证
-- ============================================================

-- 需要告诉 Lean 7 是素数
instance : Fact (Nat.Prime 7) := ⟨by decide⟩

-- ZMod p 是域（p 素数）
example : Field (ZMod 7) := inferInstance

-- Polynomial (ZMod 7) 是 EuclideanDomain（noncomputable：依赖 Field 的除法）
noncomputable example : EuclideanDomain (Polynomial (ZMod 7)) := inferInstance

-- GCDMonoid（noncomputable：依赖 EuclideanDomain）
noncomputable example : GCDMonoid (Polynomial (ZMod 7)) := inferInstance

-- DecidableEq
example : DecidableEq (Polynomial (ZMod 7)) := inferInstance

-- ============================================================
-- 2. 基本多项式构造与运算
-- ============================================================

noncomputable section

variable (p : ℕ) [hp : Fact (Nat.Prime p)]

-- 类型别名
abbrev ZpPoly := Polynomial (ZMod p)

-- 构造具体多项式
def testPoly : Polynomial (ZMod 7) := X ^ 3 + C 2 * X + C 1
def testDiv  : Polynomial (ZMod 7) := X ^ 2 + C 3

-- divByMonic / modByMonic
def testQ := testPoly /ₘ testDiv
def testR := testPoly %ₘ testDiv

end

-- ============================================================
-- 3. gcd 归一化行为
-- ============================================================

-- GCDMonoid.gcd 用于域上多项式
#check @GCDMonoid.gcd
-- EuclideanDomain.gcd 也可用
#check @EuclideanDomain.gcd

-- gcd 结果的 normalize 行为
#check @normalize_gcd

-- ============================================================
-- 4. Monic API
-- ============================================================

-- Monic 是 Prop
#check @Polynomial.Monic  -- Polynomial R → Prop

-- divByMonic 签名：不需要 Monic 证明
#check @Polynomial.divByMonic  -- Polynomial R → Polynomial R → Polynomial R

-- 关键恒等式
#check @Polynomial.modByMonic_add_div

-- ============================================================
-- 5. Squarefree / Separable
-- ============================================================

#check @Polynomial.Separable.squarefree
#check @Polynomial.separable_def

-- ============================================================
-- 6. Irreducible
-- ============================================================

#check @Irreducible
-- Irreducible p ↔ ¬IsUnit p ∧ ∀ a b, p = a * b → IsUnit a ∨ IsUnit b
