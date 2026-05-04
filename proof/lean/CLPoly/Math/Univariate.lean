/-
  L3 数学层：SparsePolyZp 算术的 Mathlib Polynomial bridge

  设计：
    SparsePolyZp（Array, 可执行）── toPoly ──> Polynomial (ZMod p) （Mathlib）

  通过 toPoly 同态把 ring 公理外包给 Mathlib 的 Polynomial 实例。
  完整 ring 公理需要：
    1. toPoly_add, toPoly_mul, toPoly_neg：homomorphism
    2. WellFormed 在算术下保持
    3. toPoly_inj_canonical：canonical SparsePolyZp 唯一对应 Polynomial

  本 commit 范围（Phase 2A.4a）：
    - Bridge 类型定义（Zp.toZMod, WellFormed, toPoly）
    - Empty / Zero 平凡情形证明
    - listSum 辅助定义 + 简单展开引理

  后续 commit（Phase 2A.4b）将补：
    - listSum_mergeAdd（核心：mergeAdd 与 list 单项式和的关系）
       难点：Zp.toZMod_add 需要 UInt64 不溢出预条件 (2p < 2^64)
    - toPoly_add, toPoly_mul, toPoly_neg
    - 最终 ring 公理（add_comm/assoc, mul_comm/assoc, distrib, ...）
-/

import CLPoly.Model
import Mathlib.Algebra.Polynomial.Basic
import Mathlib.Data.ZMod.Basic
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Linarith

namespace CLPoly.Math

open Polynomial

-- ============================================================
-- §1. Zp ↔ ZMod p bridge
-- ============================================================

-- Zp 值（带运行时 prime 字段）→ ZMod p（类型层 prime）
-- 仅当 z.prime.toNat = p 时数学语义正确（由 WellFormed 保证）
def Zp.toZMod (p : Nat) (z : Zp) : ZMod p :=
  (z.val.toNat : ZMod p)

@[simp] theorem Zp.toZMod_zero (p : Nat) :
    Zp.toZMod p { val := 0, prime := 1 } = 0 := by
  simp [Zp.toZMod]

-- ============================================================
-- §2. SparsePolyZp 良好形式（WellFormed）
-- ============================================================

-- 所有 Zp 元素都有相同的 prime（与 C++ 一致：static class-level prime）
def SparsePolyZp.WellFormed (p : Nat) (f : SparsePolyZp) : Prop :=
  ∀ x ∈ f, x.snd.prime.toNat = p

theorem SparsePolyZp.WellFormed.empty (p : Nat) :
    SparsePolyZp.WellFormed p (#[] : SparsePolyZp) := by
  intro x hx
  simp at hx

-- ============================================================
-- §3. listSum：把单项式列表加和为多项式（递归实现，便于归纳）
-- ============================================================

noncomputable def listSum (p : Nat) : List (UMonomial × Zp) → Polynomial (ZMod p)
  | [] => 0
  | (m, c) :: rest => Polynomial.monomial m.deg (Zp.toZMod p c) + listSum p rest

@[simp] theorem listSum_nil (p : Nat) : listSum p [] = 0 := rfl

@[simp] theorem listSum_cons (p : Nat) (m : UMonomial) (c : Zp) (rest : List _) :
    listSum p ((m, c) :: rest) =
      Polynomial.monomial m.deg (Zp.toZMod p c) + listSum p rest := rfl

-- 拼接的列表加和等于各自加和的和
theorem listSum_append (p : Nat) (xs ys : List (UMonomial × Zp)) :
    listSum p (xs ++ ys) = listSum p xs + listSum p ys := by
  induction xs with
  | nil => simp
  | cons x xs ih =>
    rcases x with ⟨m, c⟩
    simp [listSum, ih, add_assoc]

-- ============================================================
-- §4. SparsePolyZp.toPoly：array → Polynomial 桥
-- ============================================================

noncomputable def SparsePolyZp.toPoly (p : Nat) (f : SparsePolyZp) : Polynomial (ZMod p) :=
  listSum p f.toList

@[simp] theorem SparsePolyZp.toPoly_empty (p : Nat) :
    SparsePolyZp.toPoly p (#[] : SparsePolyZp) = 0 := by
  simp [SparsePolyZp.toPoly, listSum]

-- ============================================================
-- §5. WellFormed 在简单算术下保持（先证 negImpl，加乘留 2A.4b）
-- ============================================================

theorem SparsePolyZp.WellFormed.neg (p : Nat) (f : SparsePolyZp)
    (hf : SparsePolyZp.WellFormed p f) :
    SparsePolyZp.WellFormed p (SparsePolyZp.negImpl f) := by
  intro x hx
  simp only [SparsePolyZp.negImpl, Array.mem_map] at hx
  rcases hx with ⟨y, hy_mem, hy_eq⟩
  -- y ∈ f, x = (y.fst, -y.snd)
  -- 需要 x.snd.prime.toNat = p
  -- x.snd = -y.snd，Neg Zp 的实现为 ⟨(prime - val) % prime, prime⟩
  -- 即 (-y.snd).prime = y.snd.prime，故 (-y.snd).prime.toNat = y.snd.prime.toNat = p
  have h_prime : x.snd.prime = y.snd.prime := by
    rw [← hy_eq]
    rfl
  rw [h_prime]
  exact hf y hy_mem

-- ============================================================
-- §6. 数值验证（Mathlib decidable 通过）
-- ============================================================

-- Zp 7 的 0 → ZMod 7 的 0
example : Zp.toZMod 7 ⟨0, 7⟩ = 0 := by decide

-- Zp 7 的 3 → ZMod 7 的 3
example : Zp.toZMod 7 ⟨3, 7⟩ = (3 : ZMod 7) := by decide

-- empty SparsePolyZp 是 WellFormed
example : SparsePolyZp.WellFormed 7 #[] := SparsePolyZp.WellFormed.empty 7

end CLPoly.Math
