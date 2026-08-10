/-
  L3 数学层：Model.lean 的 ZZ / Nat 运算的形式化 spec 与证明

  Generated.Corpus 只 import Model.lean，不依赖此文件 → B2B / native 编译不受影响。
  此文件 import Mathlib 用于 ring / linarith / push_cast 等数学策略。
-/

import CLPoly.Model
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Push
import Mathlib.Data.ZMod.Basic

namespace CLPoly.Math

open Nat

-- ============================================================
-- Bezout 等式：Nat.extGcd 的核心 spec
-- 定理：a * x + b * y = g, 其中 (g, x, y) := Nat.extGcd a b
-- ============================================================

theorem Nat.extGcd_bezout (a b : Nat) :
    let (g, x, y) := Nat.extGcd a b
    (a : Int) * x + (b : Int) * y = (g : Int) := by
  induction b using Nat.strongRecOn generalizing a with
  | _ b ih =>
    rw [Nat.extGcd]
    by_cases hb : b = 0
    · simp [hb]
    · simp only [hb, dite_false]
      have hbpos : 0 < b := Nat.pos_of_ne_zero hb
      have hr_lt : a % b < b := Nat.mod_lt a hbpos
      have ih' := ih (a % b) hr_lt b
      rcases hres : Nat.extGcd b (a % b) with ⟨g, x', y'⟩
      simp only [hres] at ih'
      -- 数学骨架（在 Int）：
      --   设 q = a / b, r = a % b。
      --   IH: b * x' + r * y' = g
      --   目标: a * y' + b * (x' - q * y') = g
      --   由 a = b * q + r → a - b * q = r：
      --     a * y' + b * (x' - q * y')
      --     = a * y' + b * x' - b * q * y'
      --     = b * x' + (a - b * q) * y'
      --     = b * x' + r * y'
      --     = g
      have hcast : ((b : Int)) * ((a / b : Nat) : Int) + ((a % b : Nat) : Int) = (a : Int) := by
        have hn : b * (a / b) + a % b = a := Nat.div_add_mod a b
        exact_mod_cast hn
      -- 把 b * (a/b) 提取为单一变量给 linarith
      set bq : Int := (b : Int) * ((a / b : Nat) : Int) with hbq
      set r : Int := ((a % b : Nat) : Int) with hr_def
      -- hcast : bq + r = a
      have hr_eq : (a : Int) - bq = r := by linarith
      calc (a : Int) * y' + (b : Int) * (x' - ((a / b : Nat) : Int) * y')
          = (b : Int) * x' + ((a : Int) - bq) * y' := by
              show _ = (b : Int) * x' + ((a : Int) - (b : Int) * ((a/b : Nat) : Int)) * y'
              ring
        _ = (b : Int) * x' + r * y' := by rw [hr_eq]
        _ = (g : Int) := ih'

/-- A successful call to the concrete GMP-style `ZZ.invert` implementation
returns an actual multiplicative inverse in `ZMod m`.  This unfolds the
generated runtime implementation and derives the result from
`Nat.extGcd_bezout`; success is not treated as an abstract assumption. -/
theorem ZZ.invert_success_mul_eq_one (a : Int) (m : Nat)
    (hsuccess : (ZZ.invert 0 a (m : Int)).1 = true) :
    (a : ZMod m) * ((ZZ.invert 0 a (m : Int)).2 : ZMod m) = 1 := by
  simp only [ZZ.invert, ZZ.invertImpl] at hsuccess ⊢
  by_cases hm : m ≤ 1
  · simp [hm] at hsuccess
  · simp [hm] at hsuccess ⊢
    let am : Int := a % (m : Int)
    by_cases ham : am = 0
    · simp [am, ham] at hsuccess
    · simp [am, ham] at hsuccess ⊢
      rcases hext : Nat.extGcd am.natAbs m with ⟨d, x, y⟩
      have hextA : Nat.extGcd (a % (m : Int)).natAbs m = (d, x, y) := by
        simpa [am] using hext
      by_cases hd : d = 1
      · simp [hd]
        have hbez := Nat.extGcd_bezout am.natAbs m
        simp only [hext] at hbez
        have hmInt : (m : Int) ≠ 0 := by omega
        have hamNonneg : 0 ≤ am := by
          dsimp [am]
          exact Int.emod_nonneg a hmInt
        have hnatAbs : (am.natAbs : Int) = am :=
          Int.natAbs_of_nonneg hamNonneg
        have hbezCast := congrArg (fun z : Int => (z : ZMod m)) hbez
        simp [hnatAbs, hd] at hbezCast
        simpa [am] using hbezCast
      · simp [hextA, hd] at hsuccess

-- ============================================================
-- 数值验证（使用 Bezout 测试已实现的 extGcd）
-- ============================================================

example : (Nat.extGcd 12 18).1 = 6 := by
  -- 直接 #reduce / 计算; 用 native_decide 因 def 是 well-founded 不一定 reduce
  native_decide

example :
    let (g, x, y) := Nat.extGcd 12 18
    (12 : Int) * x + (18 : Int) * y = (g : Int) := by
  -- 由 Bezout 自动得到（定理 Nat.extGcd_bezout 的应用）
  exact Nat.extGcd_bezout 12 18

end CLPoly.Math
