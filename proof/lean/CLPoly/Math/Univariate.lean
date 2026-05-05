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

-- Zp 加法对应 ZMod 加法（前提：a.prime = b.prime = p 且 UInt64 不溢出）
-- 不溢出条件：a.val + b.val 在 UInt64 内不 wrap around
theorem Zp.toZMod_add (p : Nat) (a b : Zp)
    (ha : a.prime.toNat = p) (_hb : b.prime.toNat = p)
    (hno : a.val.toNat + b.val.toNat < UInt64.size) :
    Zp.toZMod p (a + b) = Zp.toZMod p a + Zp.toZMod p b := by
  -- Zp 加法定义：(a + b).val = (a.val + b.val) % a.prime
  show Zp.toZMod p ⟨(a.val + b.val) % a.prime, a.prime⟩ = _
  unfold Zp.toZMod
  -- 步骤 1: ((a.val + b.val) % a.prime).toNat = (a.val.toNat + b.val.toNat) % p
  have step1 : ((a.val + b.val) % a.prime).toNat = (a.val.toNat + b.val.toNat) % p := by
    rw [UInt64.toNat_mod, UInt64.toNat_add, Nat.mod_eq_of_lt hno, ha]
  rw [step1]
  -- 步骤 2: ZMod.natCast_mod 消 mod，然后 Nat.cast_add 拆开
  rw [ZMod.natCast_mod]
  push_cast
  rfl

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
-- §6. listSum_mergeAdd 核心引理
-- 证明 mergeAdd 算法保持单项式总和（在 Polynomial 环中）
-- 需要 WellFormed (prime 一致) + Reduced (val < prime) 预条件
-- 全局：2p ≤ UInt64.size 防溢出
-- ============================================================

def Zp.Reduced (p : Nat) (z : Zp) : Prop :=
  z.prime.toNat = p ∧ z.val.toNat < p

def SparsePolyZp.AllReduced (p : Nat) (xs : List (UMonomial × Zp)) : Prop :=
  ∀ x ∈ xs, Zp.Reduced p x.snd

-- 把 AllReduced 假设放进 conclusion，以便 mergeAdd.induct 自动生成的 IH 包含
theorem listSum_mergeAdd (p : Nat)
    (h2p : 2 * p ≤ UInt64.size)
    (xs ys : List (UMonomial × Zp)) :
    SparsePolyZp.AllReduced p xs → SparsePolyZp.AllReduced p ys →
    listSum p (SparsePolyZp.mergeAdd xs ys) = listSum p xs + listSum p ys := by
  induction xs, ys using SparsePolyZp.mergeAdd.induct with
  | case1 ys =>
    intro _ _
    rw [SparsePolyZp.mergeAdd]
    simp [listSum_nil]
  | case2 f fs =>
    intro _ _
    rw [SparsePolyZp.mergeAdd]
    simp [listSum_nil]
  | case3 f fs g gs hgt ih =>
    intro hxs hys
    have hxs' : SparsePolyZp.AllReduced p fs := fun x hx => hxs x (List.mem_cons_of_mem _ hx)
    rw [SparsePolyZp.mergeAdd]
    simp only [hgt, ↓reduceIte]
    rw [listSum_cons, ih hxs' hys]
    rcases f with ⟨mf, cf⟩
    rcases g with ⟨mg, cg⟩
    simp only [listSum_cons]
    ring
  | case4 f fs g gs hngt hlt ih =>
    intro hxs hys
    have hys' : SparsePolyZp.AllReduced p gs := fun x hx => hys x (List.mem_cons_of_mem _ hx)
    rw [SparsePolyZp.mergeAdd]
    simp only [hngt, ↓reduceIte, hlt]
    rw [listSum_cons, ih hxs hys']
    rcases f with ⟨mf, cf⟩
    rcases g with ⟨mg, cg⟩
    simp only [listSum_cons]
    ring
  | case5 f fs g gs hngt hnlt s heq_zero ih =>
    -- 此 case：f.deg = g.deg, s = f.snd + g.snd（let-bound）, s.val = 0（和为零）
    intro hxs hys
    have hxs' : SparsePolyZp.AllReduced p fs := fun x hx => hxs x (List.mem_cons_of_mem _ hx)
    have hys' : SparsePolyZp.AllReduced p gs := fun x hx => hys x (List.mem_cons_of_mem _ hx)
    have hf_red : Zp.Reduced p f.snd := hxs f List.mem_cons_self
    have hg_red : Zp.Reduced p g.snd := hys g List.mem_cons_self
    have hdeg : f.fst.deg = g.fst.deg := by
      rcases lt_trichotomy f.fst.deg g.fst.deg with h | h | h
      · exact absurd h hnlt
      · exact h
      · exact absurd h hngt
    have hno : f.snd.val.toNat + g.snd.val.toNat < UInt64.size := by
      have h1 := hf_red.2
      have h2 := hg_red.2
      omega
    -- LHS = listSum p (mergeAdd (f::fs) (g::gs))
    -- 由 case5 hypothesis（s.val = 0）：mergeAdd 落到 mergeAdd fs gs 分支
    have h_eq : SparsePolyZp.mergeAdd (f :: fs) (g :: gs) = SparsePolyZp.mergeAdd fs gs := by
      conv_lhs => rw [SparsePolyZp.mergeAdd]
      rw [if_neg hngt, if_neg hnlt]
      exact if_pos heq_zero
    rw [h_eq, ih hxs' hys']
    rcases f with ⟨mf, cf⟩
    rcases g with ⟨mg, cg⟩
    simp only [listSum_cons]
    -- 用 toZMod_add: 由 s.val = 0 推 cf.toZMod + cg.toZMod = 0
    have h_sum_zero : Zp.toZMod p cf + Zp.toZMod p cg = 0 := by
      rw [← Zp.toZMod_add p cf cg hf_red.1 hg_red.1 hno]
      unfold Zp.toZMod
      rw [show (cf + cg).val = 0 from heq_zero]
      simp
    have hmono_zero : Polynomial.monomial mf.deg (Zp.toZMod p cf) +
                      Polynomial.monomial mg.deg (Zp.toZMod p cg) = 0 := by
      rw [show mg.deg = mf.deg from hdeg.symm]
      rw [← Polynomial.monomial_add]
      rw [h_sum_zero]
      simp
    -- 目标重排：listSum fs + listSum gs = (mono cf + listSum fs) + (mono cg + listSum gs)
    have : Polynomial.monomial mf.deg (Zp.toZMod p cf) + listSum p fs +
           (Polynomial.monomial mg.deg (Zp.toZMod p cg) + listSum p gs) =
           (Polynomial.monomial mf.deg (Zp.toZMod p cf) +
            Polynomial.monomial mg.deg (Zp.toZMod p cg)) +
           (listSum p fs + listSum p gs) := by ring
    rw [this, hmono_zero, zero_add]
  | case6 f fs g gs hngt hnlt s heq_nonzero ih =>
    -- 此 case：f.deg = g.deg, s.val ≠ 0（和非零）
    intro hxs hys
    have hxs' : SparsePolyZp.AllReduced p fs := fun x hx => hxs x (List.mem_cons_of_mem _ hx)
    have hys' : SparsePolyZp.AllReduced p gs := fun x hx => hys x (List.mem_cons_of_mem _ hx)
    have hf_red : Zp.Reduced p f.snd := hxs f List.mem_cons_self
    have hg_red : Zp.Reduced p g.snd := hys g List.mem_cons_self
    have hdeg : f.fst.deg = g.fst.deg := by
      rcases lt_trichotomy f.fst.deg g.fst.deg with h | h | h
      · exact absurd h hnlt
      · exact h
      · exact absurd h hngt
    have hno : f.snd.val.toNat + g.snd.val.toNat < UInt64.size := by
      have h1 := hf_red.2
      have h2 := hg_red.2
      omega
    have h_eq : SparsePolyZp.mergeAdd (f :: fs) (g :: gs) =
        (f.fst, f.snd + g.snd) :: SparsePolyZp.mergeAdd fs gs := by
      conv_lhs => rw [SparsePolyZp.mergeAdd]
      rw [if_neg hngt, if_neg hnlt]
      exact if_neg heq_nonzero
    rw [h_eq]
    rw [listSum_cons, ih hxs' hys']
    rcases f with ⟨mf, cf⟩
    rcases g with ⟨mg, cg⟩
    simp only [listSum_cons]
    rw [Zp.toZMod_add p cf cg hf_red.1 hg_red.1 hno]
    rw [Polynomial.monomial_add]
    rw [show mf.deg = mg.deg from hdeg]
    ring

-- ============================================================
-- §6b. Zp.toZMod 在 neg / mul 下的 homomorphism
-- ============================================================

-- Zp 取负对应 ZMod 取负（前提：a.prime = p，a.val ≤ a.prime 不下溢）
theorem Zp.toZMod_neg (p : Nat) (a : Zp) (ha : a.prime.toNat = p)
    (hred : a.val.toNat < p) :
    Zp.toZMod p (-a) = -Zp.toZMod p a := by
  show Zp.toZMod p ⟨(a.prime - a.val) % a.prime, a.prime⟩ = -Zp.toZMod p a
  unfold Zp.toZMod
  -- (a.prime - a.val).toNat = p - a.val.toNat（由 a.val ≤ a.prime 不下溢）
  have h_le : a.val ≤ a.prime := by
    rw [show a.val ≤ a.prime ↔ a.val.toNat ≤ a.prime.toNat from UInt64.le_iff_toNat_le]
    omega
  have h_sub : (a.prime - a.val).toNat = a.prime.toNat - a.val.toNat :=
    UInt64.toNat_sub_of_le _ _ h_le
  rw [UInt64.toNat_mod, h_sub, ha]
  -- 目标：((p - a.val.toNat) % p : ZMod p) = -↑a.val.toNat
  rw [ZMod.natCast_mod]
  -- (p - n : Nat) % p = p - n when n < p, else 0; 直接计算
  -- ↑(p - a.val.toNat) = ↑p - ↑a.val.toNat (in ZMod p)
  -- ↑p = 0
  rw [Nat.cast_sub (by omega : a.val.toNat ≤ p)]
  rw [ZMod.natCast_self, zero_sub]

-- Zp 乘法对应 ZMod 乘法（前提：a.prime = b.prime = p，无溢出 p^2 ≤ UInt64.size）
theorem Zp.toZMod_mul (p : Nat) (a b : Zp)
    (ha : a.prime.toNat = p) (_hb : b.prime.toNat = p)
    (hno : a.val.toNat * b.val.toNat < UInt64.size) :
    Zp.toZMod p (a * b) = Zp.toZMod p a * Zp.toZMod p b := by
  show Zp.toZMod p ⟨(a.val * b.val) % a.prime, a.prime⟩ = _
  unfold Zp.toZMod
  rw [UInt64.toNat_mod, UInt64.toNat_mul, Nat.mod_eq_of_lt hno, ha]
  rw [ZMod.natCast_mod]
  push_cast
  rfl

-- ============================================================
-- §7. SparsePolyZp.WellFormed_arr：Array 版本 + toPoly_add/_neg/_sub/_mul
-- ============================================================

-- Array 版本（包装 toList）
def SparsePolyZp.WellFormed_arr (p : Nat) (f : SparsePolyZp) : Prop :=
  SparsePolyZp.AllReduced p f.toList

-- toPoly_add: 由 listSum_mergeAdd 直接推
theorem SparsePolyZp.toPoly_add (p : Nat) (h2p : 2 * p ≤ UInt64.size)
    (f g : SparsePolyZp)
    (hf : SparsePolyZp.WellFormed_arr p f) (hg : SparsePolyZp.WellFormed_arr p g) :
    SparsePolyZp.toPoly p (f + g) = SparsePolyZp.toPoly p f + SparsePolyZp.toPoly p g := by
  show listSum p (SparsePolyZp.addImpl f g).toList = listSum p f.toList + listSum p g.toList
  have h_toList : (SparsePolyZp.addImpl f g).toList = SparsePolyZp.mergeAdd f.toList g.toList := by
    unfold SparsePolyZp.addImpl
    simp
  rw [h_toList]
  exact listSum_mergeAdd p h2p f.toList g.toList hf hg

-- listSum_map_neg: 取负在列表上的同态
theorem listSum_map_neg (p : Nat) (xs : List (UMonomial × Zp))
    (hxs : SparsePolyZp.AllReduced p xs) :
    listSum p (xs.map (fun (m, c) => (m, -c))) = -listSum p xs := by
  induction xs with
  | nil => simp
  | cons x xs ih =>
    rcases x with ⟨m, c⟩
    have hc_red : Zp.Reduced p c := hxs (m, c) List.mem_cons_self
    have hxs' : SparsePolyZp.AllReduced p xs := fun y hy => hxs y (List.mem_cons_of_mem _ hy)
    simp only [List.map_cons, listSum_cons]
    rw [ih hxs', Zp.toZMod_neg p c hc_red.1 hc_red.2]
    rw [Polynomial.monomial_neg]
    ring

-- toPoly_neg: 取负的 polynomial 同态
theorem SparsePolyZp.toPoly_neg (p : Nat) (f : SparsePolyZp)
    (hf : SparsePolyZp.WellFormed_arr p f) :
    SparsePolyZp.toPoly p (-f) = -SparsePolyZp.toPoly p f := by
  show listSum p (SparsePolyZp.negImpl f).toList = -listSum p f.toList
  unfold SparsePolyZp.negImpl
  -- (f.map g).toList = f.toList.map g
  rw [Array.toList_map]
  exact listSum_map_neg p f.toList hf

-- toPoly_sub: 由 toPoly_add + toPoly_neg 推
-- WellFormed_arr 在 negImpl 下需保持（neg 不影响 prime）
theorem SparsePolyZp.WellFormed_arr.neg (p : Nat) (f : SparsePolyZp)
    (hf : SparsePolyZp.WellFormed_arr p f) :
    SparsePolyZp.WellFormed_arr p (-f) := by
  intro x hx
  show Zp.Reduced p x.snd
  -- (-f) = SparsePolyZp.negImpl f = f.map (fun (m, c) => (m, -c))
  -- toList 后用 List.mem_map 解构成员关系
  simp only [SparsePolyZp.WellFormed_arr, Neg.neg, SparsePolyZp.negImpl,
             Array.toList_map, List.mem_map] at hx
  rcases hx with ⟨y, hy_mem, hy_eq⟩
  have hy_red : Zp.Reduced p y.snd := hf y hy_mem
  rw [← hy_eq]
  refine ⟨?_, ?_⟩
  · show y.snd.prime.toNat = p
    exact hy_red.1
  · show ((y.snd.prime - y.snd.val) % y.snd.prime).toNat < p
    rw [UInt64.toNat_mod]
    have hp_eq : y.snd.prime.toNat = p := hy_red.1
    have hp_pos : 0 < p := by
      have := hy_red.2  -- y.snd.val.toNat < p
      omega
    rw [hp_eq]
    exact Nat.mod_lt _ hp_pos

theorem SparsePolyZp.toPoly_sub (p : Nat) (h2p : 2 * p ≤ UInt64.size)
    (f g : SparsePolyZp)
    (hf : SparsePolyZp.WellFormed_arr p f) (hg : SparsePolyZp.WellFormed_arr p g) :
    SparsePolyZp.toPoly p (f - g) = SparsePolyZp.toPoly p f - SparsePolyZp.toPoly p g := by
  -- f - g = subImpl f g = addImpl f (negImpl g)
  show listSum p (SparsePolyZp.subImpl f g).toList = listSum p f.toList - listSum p g.toList
  unfold SparsePolyZp.subImpl
  -- 转 toPoly 后用 add + neg
  have h1 : listSum p (SparsePolyZp.addImpl f (SparsePolyZp.negImpl g)).toList =
      listSum p f.toList + listSum p (SparsePolyZp.negImpl g).toList := by
    have hg_neg : SparsePolyZp.WellFormed_arr p (-g) :=
      SparsePolyZp.WellFormed_arr.neg p g hg
    have h_toList : (SparsePolyZp.addImpl f (SparsePolyZp.negImpl g)).toList =
        SparsePolyZp.mergeAdd f.toList (SparsePolyZp.negImpl g).toList := by
      unfold SparsePolyZp.addImpl; simp
    rw [h_toList]
    exact listSum_mergeAdd p h2p f.toList (SparsePolyZp.negImpl g).toList hf hg_neg
  rw [h1]
  -- listSum p (negImpl g).toList = -listSum p g.toList
  unfold SparsePolyZp.negImpl
  rw [Array.toList_map]
  rw [listSum_map_neg p g.toList hg]
  ring

-- ============================================================
-- §7c. WellFormed 在加减下保持
-- ============================================================

-- mergeAdd 保持 AllReduced（核心：归纳 + case5 输出新 Zp 仍是 Reduced）
theorem mergeAdd_AllReduced (p : Nat) (xs ys : List (UMonomial × Zp)) :
    SparsePolyZp.AllReduced p xs → SparsePolyZp.AllReduced p ys →
    SparsePolyZp.AllReduced p (SparsePolyZp.mergeAdd xs ys) := by
  induction xs, ys using SparsePolyZp.mergeAdd.induct with
  | case1 ys =>
    intro _ hys
    rw [SparsePolyZp.mergeAdd]
    exact hys
  | case2 f fs =>
    intro hfs _
    rw [SparsePolyZp.mergeAdd]
    exact hfs
  | case3 f fs g gs hgt ih =>
    intro hxs hys
    have hxs' : SparsePolyZp.AllReduced p fs :=
      fun y hy => hxs y (List.mem_cons_of_mem _ hy)
    have h_eq : SparsePolyZp.mergeAdd (f :: fs) (g :: gs) =
        f :: SparsePolyZp.mergeAdd fs (g :: gs) := by
      conv_lhs => rw [SparsePolyZp.mergeAdd]
      rw [if_pos hgt]
    rw [h_eq]
    intro x hx
    simp only [List.mem_cons] at hx
    rcases hx with hx_eq | hx_in
    · rw [hx_eq]; exact hxs f List.mem_cons_self
    · exact ih hxs' hys x hx_in
  | case4 f fs g gs hngt hlt ih =>
    intro hxs hys
    have hys' : SparsePolyZp.AllReduced p gs :=
      fun y hy => hys y (List.mem_cons_of_mem _ hy)
    have h_eq : SparsePolyZp.mergeAdd (f :: fs) (g :: gs) =
        g :: SparsePolyZp.mergeAdd (f :: fs) gs := by
      conv_lhs => rw [SparsePolyZp.mergeAdd]
      rw [if_neg hngt, if_pos hlt]
    rw [h_eq]
    intro x hx
    simp only [List.mem_cons] at hx
    rcases hx with hx_eq | hx_in
    · rw [hx_eq]; exact hys g List.mem_cons_self
    · exact ih hxs hys' x hx_in
  | case5 f fs g gs hngt hnlt s heq_zero ih =>
    intro hxs hys
    have hxs' : SparsePolyZp.AllReduced p fs :=
      fun y hy => hxs y (List.mem_cons_of_mem _ hy)
    have hys' : SparsePolyZp.AllReduced p gs :=
      fun y hy => hys y (List.mem_cons_of_mem _ hy)
    have h_eq : SparsePolyZp.mergeAdd (f :: fs) (g :: gs) = SparsePolyZp.mergeAdd fs gs := by
      conv_lhs => rw [SparsePolyZp.mergeAdd]
      rw [if_neg hngt, if_neg hnlt]; exact if_pos heq_zero
    rw [h_eq]
    exact ih hxs' hys'
  | case6 f fs g gs hngt hnlt s heq_nonzero ih =>
    intro hxs hys
    have hxs' : SparsePolyZp.AllReduced p fs :=
      fun y hy => hxs y (List.mem_cons_of_mem _ hy)
    have hys' : SparsePolyZp.AllReduced p gs :=
      fun y hy => hys y (List.mem_cons_of_mem _ hy)
    have hf_red : Zp.Reduced p f.snd := hxs f List.mem_cons_self
    have h_eq : SparsePolyZp.mergeAdd (f :: fs) (g :: gs) =
        (f.fst, f.snd + g.snd) :: SparsePolyZp.mergeAdd fs gs := by
      conv_lhs => rw [SparsePolyZp.mergeAdd]
      rw [if_neg hngt, if_neg hnlt]; exact if_neg heq_nonzero
    rw [h_eq]
    intro x hx
    simp only [List.mem_cons] at hx
    rcases hx with hx_eq | hx_in
    · rw [hx_eq]
      -- 新元素 (f.fst, f.snd + g.snd) 是 Reduced p
      refine ⟨?_, ?_⟩
      · show f.snd.prime.toNat = p
        exact hf_red.1
      · show ((f.snd.val + g.snd.val) % f.snd.prime).toNat < p
        rw [UInt64.toNat_mod, hf_red.1]
        have : 0 < p := lt_of_le_of_lt (Nat.zero_le _) hf_red.2
        exact Nat.mod_lt _ this
    · exact ih hxs' hys' x hx_in

theorem SparsePolyZp.WellFormed_arr.add (p : Nat) (f g : SparsePolyZp)
    (hf : SparsePolyZp.WellFormed_arr p f) (hg : SparsePolyZp.WellFormed_arr p g) :
    SparsePolyZp.WellFormed_arr p (f + g) := by
  -- f + g = addImpl f g = (mergeAdd f.toList g.toList).toArray
  -- (f + g).toList = mergeAdd f.toList g.toList
  show SparsePolyZp.AllReduced p (SparsePolyZp.addImpl f g).toList
  unfold SparsePolyZp.addImpl
  -- (List.toArray l).toList = l 由 simp 自动得
  have h : (SparsePolyZp.mergeAdd f.toList g.toList).toArray.toList =
           SparsePolyZp.mergeAdd f.toList g.toList := by simp
  rw [h]
  exact mergeAdd_AllReduced p f.toList g.toList hf hg

theorem SparsePolyZp.WellFormed_arr.sub (p : Nat) (f g : SparsePolyZp)
    (hf : SparsePolyZp.WellFormed_arr p f) (hg : SparsePolyZp.WellFormed_arr p g) :
    SparsePolyZp.WellFormed_arr p (f - g) := by
  -- f - g = addImpl f (negImpl g)
  show SparsePolyZp.AllReduced p (SparsePolyZp.subImpl f g).toList
  unfold SparsePolyZp.subImpl
  show SparsePolyZp.AllReduced p (SparsePolyZp.addImpl f (SparsePolyZp.negImpl g)).toList
  -- 由 WellFormed_arr.neg 得 negImpl g 是 WF；再用 add 保持性
  have hg_neg : SparsePolyZp.WellFormed_arr p (-g) :=
    SparsePolyZp.WellFormed_arr.neg p g hg
  -- (-g) = negImpl g 定义上相等
  have h_eq : (-g : SparsePolyZp) = SparsePolyZp.negImpl g := rfl
  rw [h_eq] at hg_neg
  -- 此时 hg_neg : WellFormed_arr p (negImpl g)
  -- 转 WellFormed.add
  show SparsePolyZp.AllReduced p
      ((SparsePolyZp.addImpl f (SparsePolyZp.negImpl g)) : SparsePolyZp).toList
  -- f + (negImpl g) 当成 WellFormed.add 应用
  have h_add : SparsePolyZp.WellFormed_arr p (f + SparsePolyZp.negImpl g) :=
    SparsePolyZp.WellFormed_arr.add p f (SparsePolyZp.negImpl g) hf hg_neg
  -- f + (negImpl g) = addImpl f (negImpl g)
  exact h_add

-- ============================================================
-- §7d. Ring 公理（加法部分）via toPoly bridge
-- 这些定理把 SparsePolyZp 加法的 ring 公理 reduce 到 Polynomial 已有的 ring 公理。
-- 形式：toPoly p (LHS) = toPoly p (RHS)（不直接 Array 等式；后者需 toPoly_inj）
-- ============================================================

-- 加法交换律
theorem SparsePolyZp.add_comm_via_toPoly (p : Nat) (h2p : 2 * p ≤ UInt64.size)
    (f g : SparsePolyZp)
    (hf : SparsePolyZp.WellFormed_arr p f) (hg : SparsePolyZp.WellFormed_arr p g) :
    SparsePolyZp.toPoly p (f + g) = SparsePolyZp.toPoly p (g + f) := by
  rw [SparsePolyZp.toPoly_add p h2p f g hf hg, SparsePolyZp.toPoly_add p h2p g f hg hf]
  ring

-- 加法结合律
theorem SparsePolyZp.add_assoc_via_toPoly (p : Nat) (h2p : 2 * p ≤ UInt64.size)
    (f g h : SparsePolyZp)
    (hf : SparsePolyZp.WellFormed_arr p f)
    (hg : SparsePolyZp.WellFormed_arr p g)
    (hh : SparsePolyZp.WellFormed_arr p h) :
    SparsePolyZp.toPoly p ((f + g) + h) = SparsePolyZp.toPoly p (f + (g + h)) := by
  have hfg := SparsePolyZp.WellFormed_arr.add p f g hf hg
  have hgh := SparsePolyZp.WellFormed_arr.add p g h hg hh
  rw [SparsePolyZp.toPoly_add p h2p _ _ hfg hh,
      SparsePolyZp.toPoly_add p h2p f g hf hg,
      SparsePolyZp.toPoly_add p h2p f _ hf hgh,
      SparsePolyZp.toPoly_add p h2p g h hg hh]
  ring

-- 零元（左/右）：toPoly p (0 + f) = toPoly p f, toPoly p (f + 0) = toPoly p f
-- 注意：0 在 SparsePolyZp 是 #[]（OfNat instance），WellFormed_arr p 0 trivially holds.
theorem SparsePolyZp.zero_add_via_toPoly (p : Nat) (h2p : 2 * p ≤ UInt64.size)
    (f : SparsePolyZp) (hf : SparsePolyZp.WellFormed_arr p f) :
    SparsePolyZp.toPoly p (0 + f) = SparsePolyZp.toPoly p f := by
  have h0 : SparsePolyZp.WellFormed_arr p (0 : SparsePolyZp) := by
    intro x hx
    simp [show (0 : SparsePolyZp) = #[] from rfl] at hx
  rw [SparsePolyZp.toPoly_add p h2p _ _ h0 hf]
  show SparsePolyZp.toPoly p (#[] : SparsePolyZp) + _ = _
  rw [SparsePolyZp.toPoly_empty]
  ring

-- 加法逆元：toPoly p (f - f) = 0
theorem SparsePolyZp.sub_self_via_toPoly (p : Nat) (h2p : 2 * p ≤ UInt64.size)
    (f : SparsePolyZp) (hf : SparsePolyZp.WellFormed_arr p f) :
    SparsePolyZp.toPoly p (f - f) = 0 := by
  rw [SparsePolyZp.toPoly_sub p h2p f f hf hf]
  ring

-- ============================================================
-- §8. 数值验证（Mathlib decidable 通过）
-- ============================================================
--
-- 余下工作（未在本 commit 完成；难度递增）：
-- - listSum_filterMap_scale + toPoly_scaleByMonomial（标量乘法同态）
-- - WellFormed_arr.scaleByMonomial 保持性
-- - toPoly_mul（foldl 分配律 — 最复杂）
-- - WellFormed_arr.mul 保持性
-- - 乘法相关 ring 公理（mul_comm/_assoc, distrib）via toPoly bridge
-- - toPoly_inj_canonical（用于反推 Array 等式 spec）
-- ============================================================

-- Zp 7 的 0 → ZMod 7 的 0
example : Zp.toZMod 7 ⟨0, 7⟩ = 0 := by decide

-- Zp 7 的 3 → ZMod 7 的 3
example : Zp.toZMod 7 ⟨3, 7⟩ = (3 : ZMod 7) := by decide

-- empty SparsePolyZp 是 WellFormed
example : SparsePolyZp.WellFormed 7 #[] := SparsePolyZp.WellFormed.empty 7

end CLPoly.Math
