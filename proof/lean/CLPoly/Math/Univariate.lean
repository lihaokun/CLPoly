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
import Mathlib.Algebra.Polynomial.Coeff
import Mathlib.Algebra.Polynomial.Degree.Definitions
import Mathlib.Data.List.Chain
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
  -- Add Zp 实例（Nat-based 防大素数溢出）：
  -- (a + b).val = ((a.val.toNat + b.val.toNat) % a.prime.toNat).toUInt64
  show Zp.toZMod p ⟨((a.val.toNat + b.val.toNat) % a.prime.toNat).toUInt64, a.prime⟩ = _
  unfold Zp.toZMod
  -- 关键 round-trip: (((a+b)%p).toUInt64).toNat = (a+b)%p（因 (a+b)%p ≤ a+b < UInt64.size）
  have h_mod_le : (a.val.toNat + b.val.toNat) % a.prime.toNat ≤
                  a.val.toNat + b.val.toNat := Nat.mod_le _ _
  have h_mod_lt : (a.val.toNat + b.val.toNat) % a.prime.toNat < UInt64.size :=
    Nat.lt_of_le_of_lt h_mod_le hno
  have step1 : (((a.val.toNat + b.val.toNat) % a.prime.toNat).toUInt64).toNat =
               (a.val.toNat + b.val.toNat) % p := by
    change (OfNat.ofNat _ : UInt64).toNat = _
    rw [UInt64.toNat_ofNat, Nat.mod_eq_of_lt h_mod_lt, ha]
  rw [step1, ZMod.natCast_mod]
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
  -- Mul Zp 实例（Nat-based）：(a*b).val = ((a.val.toNat * b.val.toNat) % a.prime.toNat).toUInt64
  show Zp.toZMod p ⟨((a.val.toNat * b.val.toNat) % a.prime.toNat).toUInt64, a.prime⟩ = _
  unfold Zp.toZMod
  have h_mod_le : (a.val.toNat * b.val.toNat) % a.prime.toNat ≤
                  a.val.toNat * b.val.toNat := Nat.mod_le _ _
  have h_mod_lt : (a.val.toNat * b.val.toNat) % a.prime.toNat < UInt64.size :=
    Nat.lt_of_le_of_lt h_mod_le hno
  have step1 : (((a.val.toNat * b.val.toNat) % a.prime.toNat).toUInt64).toNat =
               (a.val.toNat * b.val.toNat) % p := by
    change (OfNat.ofNat _ : UInt64).toNat = _
    rw [UInt64.toNat_ofNat, Nat.mod_eq_of_lt h_mod_lt, ha]
  rw [step1, ZMod.natCast_mod]
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
      -- (f.snd + g.snd).val = ((f.snd.val.toNat + g.snd.val.toNat) % f.snd.prime.toNat).toUInt64
      refine ⟨?_, ?_⟩
      · show f.snd.prime.toNat = p
        exact hf_red.1
      · -- 目标:((((f.snd.val.toNat + g.snd.val.toNat) % f.snd.prime.toNat).toUInt64)).toNat < p
        show (((f.snd.val.toNat + g.snd.val.toNat) % f.snd.prime.toNat).toUInt64).toNat < p
        have hp_pos : 0 < p := lt_of_le_of_lt (Nat.zero_le _) hf_red.2
        have h_mod_lt_prime : (f.snd.val.toNat + g.snd.val.toNat) % f.snd.prime.toNat <
                              f.snd.prime.toNat := by
          apply Nat.mod_lt; rw [hf_red.1]; exact hp_pos
        have h_mod_lt_size : (f.snd.val.toNat + g.snd.val.toNat) % f.snd.prime.toNat <
                             UInt64.size :=
          lt_of_lt_of_le h_mod_lt_prime (Nat.le_of_lt f.snd.prime.toNat_lt)
        change (OfNat.ofNat _ : UInt64).toNat < p
        rw [UInt64.toNat_ofNat, Nat.mod_eq_of_lt h_mod_lt_size, hf_red.1]
        rw [hf_red.1] at h_mod_lt_prime
        exact h_mod_lt_prime
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
-- §7e. scaleByMonomial 同态 + WellFormed (Phase 2A.4c Stage 1)
-- ============================================================

-- Lemma 1: scaleByMonomial 内层 filterMap 在 list 上对应 polynomial 单项式乘
-- 每项 (mf, cf) 映射到 monomial(m+mf, c*cf)；零项被 filterMap 丢掉
theorem listSum_filterMap_scale (p : Nat) (h_p2 : p * p ≤ UInt64.size)
    (m : UMonomial) (c : Zp) (hc : Zp.Reduced p c)
    (xs : List (UMonomial × Zp)) (hxs : SparsePolyZp.AllReduced p xs) :
    listSum p (xs.filterMap (fun term : UMonomial × Zp =>
      if (c * term.snd).val = 0 then none
      else some (UMonomial.mk (m.deg + term.fst.deg), c * term.snd)))
    = Polynomial.monomial m.deg (Zp.toZMod p c) * listSum p xs := by
  induction xs with
  | nil =>
    simp [listSum_nil, mul_zero]
  | cons x rest ih =>
    rcases x with ⟨mf, cf⟩
    have hcf_red : Zp.Reduced p cf := hxs (mf, cf) List.mem_cons_self
    have hxs' : SparsePolyZp.AllReduced p rest := fun y hy => hxs y (List.mem_cons_of_mem _ hy)
    -- c.val * cf.val < p * p ≤ UInt64.size，故 Zp.toZMod_mul 前提满足
    have h_no_overflow : c.val.toNat * cf.val.toNat < UInt64.size := by
      have h1 : c.val.toNat < p := hc.2
      have h2 : cf.val.toNat < p := hcf_red.2
      have hp_pos : 0 < p := Nat.zero_lt_of_lt h1
      calc c.val.toNat * cf.val.toNat
          < p * p := Nat.mul_lt_mul_of_lt_of_le h1 (Nat.le_of_lt h2) hp_pos
        _ ≤ UInt64.size := h_p2
    have h_toZMod_mul : Zp.toZMod p (c * cf) = Zp.toZMod p c * Zp.toZMod p cf :=
      Zp.toZMod_mul p c cf hc.1 hcf_red.1 h_no_overflow
    -- monomial 乘法等式：monomial(m+mf, toZMod(c*cf)) = monomial m (toZMod c) * monomial mf (toZMod cf)
    have h_mono_eq : Polynomial.monomial (m.deg + mf.deg) (Zp.toZMod p (c * cf))
                    = Polynomial.monomial m.deg (Zp.toZMod p c) *
                      Polynomial.monomial mf.deg (Zp.toZMod p cf) := by
      rw [h_toZMod_mul, Polynomial.monomial_mul_monomial]
    by_cases h_zero : (c * cf).val = 0
    · -- filterMap 在该位 drop
      simp only [List.filterMap_cons, h_zero, ↓reduceIte]
      rw [ih hxs', listSum_cons, mul_add]
      -- 证 monomial m.deg (toZMod c) * monomial mf.deg (toZMod cf) = 0
      have h_toZMod_zero : Zp.toZMod p (c * cf) = 0 := by
        unfold Zp.toZMod
        rw [show (c * cf).val.toNat = 0 from by rw [h_zero]; rfl]
        simp
      have h_mono_zero :
          Polynomial.monomial m.deg (Zp.toZMod p c) *
          Polynomial.monomial mf.deg (Zp.toZMod p cf) = 0 := by
        rw [← h_mono_eq, h_toZMod_zero]; exact Polynomial.monomial_zero_right _
      rw [h_mono_zero, zero_add]
    · -- filterMap 在该位 keep
      simp only [List.filterMap_cons, h_zero, ↓reduceIte]
      rw [listSum_cons, ih hxs']
      rw [listSum_cons, mul_add, h_mono_eq]

-- Lemma 2: toPoly_scaleByMonomial — 整 SparsePoly 对应 polynomial 单项式乘
theorem SparsePolyZp.toPoly_scaleByMonomial (p : Nat) (h_p2 : p * p ≤ UInt64.size)
    (m : UMonomial) (c : Zp) (hc : Zp.Reduced p c)
    (f : SparsePolyZp) (hf : SparsePolyZp.WellFormed_arr p f) :
    SparsePolyZp.toPoly p (SparsePolyZp.scaleByMonomial m c f)
    = Polynomial.monomial m.deg (Zp.toZMod p c) * SparsePolyZp.toPoly p f := by
  unfold SparsePolyZp.scaleByMonomial
  by_cases h_c_zero : c.val = 0
  · -- c.val = 0：scaleByMonomial 返 #[]；同时 toZMod c = 0 → RHS = 0 * _ = 0
    simp only [h_c_zero, ↓reduceIte]
    show SparsePolyZp.toPoly p (#[] : SparsePolyZp) = _
    rw [SparsePolyZp.toPoly_empty p]
    have h_toZMod_c_zero : Zp.toZMod p c = 0 := by
      unfold Zp.toZMod
      rw [show c.val.toNat = 0 from by rw [h_c_zero]; rfl]
      simp
    rw [h_toZMod_c_zero, Polynomial.monomial_zero_right, zero_mul]
  · -- c.val ≠ 0：套 Lemma 1
    simp only [h_c_zero, ↓reduceIte]
    show listSum p (Array.filterMap _ f).toList = _
    rw [Array.toList_filterMap]
    -- 现在目标 listSum p (f.toList.filterMap _) = monomial _ _ * toPoly p f
    -- toPoly p f = listSum p f.toList
    show listSum p (f.toList.filterMap _) = _ * listSum p f.toList
    exact listSum_filterMap_scale p h_p2 m c hc f.toList hf

-- Lemma 3: WellFormed_arr 在 scaleByMonomial 下闭合
theorem SparsePolyZp.WellFormed_arr.scaleByMonomial (p : Nat)
    (h_p2 : p * p ≤ UInt64.size)
    (m : UMonomial) (c : Zp) (hc : Zp.Reduced p c)
    (f : SparsePolyZp) (hf : SparsePolyZp.WellFormed_arr p f) :
    SparsePolyZp.WellFormed_arr p (SparsePolyZp.scaleByMonomial m c f) := by
  unfold SparsePolyZp.scaleByMonomial
  by_cases h_c_zero : c.val = 0
  · -- c.val = 0：scaleByMonomial 返 #[]，平凡 WellFormed
    simp only [h_c_zero, ↓reduceIte]
    intro x hx
    simp at hx
  · -- c.val ≠ 0：每项 = (m+mf, c*cy) 中 c*cy 的 prime = c.prime；val < c.prime
    simp only [h_c_zero, ↓reduceIte]
    intro x hx
    rw [Array.toList_filterMap, List.mem_filterMap] at hx
    rcases hx with ⟨y, hy_mem, hy_some⟩
    rcases y with ⟨my, cy⟩
    have hcy_red : Zp.Reduced p cy := hf (my, cy) hy_mem
    -- hy_some 形如 `(if ... then none else some (..., c*cy)) = some x`
    simp only at hy_some
    by_cases h_prod_zero : (c * cy).val = 0
    · simp [h_prod_zero] at hy_some
    · simp only [h_prod_zero, ↓reduceIte, Option.some.injEq] at hy_some
      -- hy_some : ({ deg := m.deg + my.deg }, c * cy) = x
      subst hy_some
      -- 目标：Zp.Reduced p (c * cy)
      refine ⟨?_, ?_⟩
      · -- (c*cy).prime.toNat = p：Mul Zp 实例返 a.prime
        show c.prime.toNat = p
        exact hc.1
      · -- (c*cy).val.toNat < p
        -- Mul Zp 实例: (c*cy).val = ((c.val.toNat * cy.val.toNat) % c.prime.toNat).toUInt64
        show (((c.val.toNat * cy.val.toNat) % c.prime.toNat).toUInt64).toNat < p
        have hp_pos : 0 < p := Nat.zero_lt_of_lt hc.2
        have h_no : c.val.toNat * cy.val.toNat < UInt64.size := by
          calc c.val.toNat * cy.val.toNat
              < p * p := Nat.mul_lt_mul_of_lt_of_le hc.2 (Nat.le_of_lt hcy_red.2) hp_pos
            _ ≤ UInt64.size := h_p2
        have h_mod_lt_prime : (c.val.toNat * cy.val.toNat) % c.prime.toNat <
                              c.prime.toNat := by
          apply Nat.mod_lt; rw [hc.1]; exact hp_pos
        have h_mod_lt_size : (c.val.toNat * cy.val.toNat) % c.prime.toNat <
                             UInt64.size :=
          Nat.lt_of_le_of_lt (Nat.mod_le _ _) h_no
        change (OfNat.ofNat _ : UInt64).toNat < p
        rw [UInt64.toNat_ofNat, Nat.mod_eq_of_lt h_mod_lt_size, hc.1]
        rw [hc.1] at h_mod_lt_prime
        exact h_mod_lt_prime

-- ============================================================
-- §7f. toPoly_mul + WellFormed_arr.mul (Phase 2A.4c Stage 2)
-- ============================================================

-- Helper: 对 List.foldl 的不变量
-- foldl 累积 acc：每步 acc' = acc + scaleByMonomial (mf, cf) g
-- 不变量：toPoly p (foldl ... acc xs) = toPoly p acc + listSum p xs * toPoly p g
theorem toPoly_foldl_mulStep (p : Nat)
    (h_2p : 2 * p ≤ UInt64.size) (h_p2 : p * p ≤ UInt64.size)
    (g : SparsePolyZp) (hg : SparsePolyZp.WellFormed_arr p g)
    (xs : List (UMonomial × Zp)) :
    ∀ acc : SparsePolyZp,
    SparsePolyZp.WellFormed_arr p acc → SparsePolyZp.AllReduced p xs →
    SparsePolyZp.toPoly p (List.foldl (fun a (t : UMonomial × Zp) =>
        SparsePolyZp.addImpl a (SparsePolyZp.scaleByMonomial t.fst t.snd g)) acc xs)
    = SparsePolyZp.toPoly p acc + listSum p xs * SparsePolyZp.toPoly p g := by
  induction xs with
  | nil =>
    intro acc _ _
    simp [listSum_nil, zero_mul, add_zero]
  | cons x rest ih =>
    intro acc hacc hxs
    rcases x with ⟨mf, cf⟩
    have hcf_red : Zp.Reduced p cf := hxs (mf, cf) List.mem_cons_self
    have hxs' : SparsePolyZp.AllReduced p rest :=
      fun y hy => hxs y (List.mem_cons_of_mem _ hy)
    have h_scale_wf : SparsePolyZp.WellFormed_arr p
        (SparsePolyZp.scaleByMonomial mf cf g) :=
      SparsePolyZp.WellFormed_arr.scaleByMonomial p h_p2 mf cf hcf_red g hg
    have hacc' : SparsePolyZp.WellFormed_arr p
        (SparsePolyZp.addImpl acc (SparsePolyZp.scaleByMonomial mf cf g)) :=
      SparsePolyZp.WellFormed_arr.add p acc _ hacc h_scale_wf
    rw [List.foldl_cons]
    rw [ih _ hacc' hxs']
    -- 目标: toPoly p (addImpl acc (scale mf cf g)) + listSum rest * toPoly g
    --     = toPoly p acc + listSum ((mf,cf) :: rest) * toPoly g
    show SparsePolyZp.toPoly p (acc + SparsePolyZp.scaleByMonomial mf cf g) +
         listSum p rest * SparsePolyZp.toPoly p g = _
    rw [SparsePolyZp.toPoly_add p h_2p acc _ hacc h_scale_wf]
    rw [SparsePolyZp.toPoly_scaleByMonomial p h_p2 mf cf hcf_red g hg]
    rw [listSum_cons, add_mul]
    ring

-- Lemma 4 (核心): toPoly_mul — SparsePolyZp 乘法的 polynomial 同态
theorem SparsePolyZp.toPoly_mul (p : Nat)
    (h_2p : 2 * p ≤ UInt64.size) (h_p2 : p * p ≤ UInt64.size)
    (f g : SparsePolyZp)
    (hf : SparsePolyZp.WellFormed_arr p f) (hg : SparsePolyZp.WellFormed_arr p g) :
    SparsePolyZp.toPoly p (f * g) = SparsePolyZp.toPoly p f * SparsePolyZp.toPoly p g := by
  show SparsePolyZp.toPoly p (SparsePolyZp.mulImpl f g) = _
  unfold SparsePolyZp.mulImpl
  rw [← Array.foldl_toList]
  have h_empty_wf : SparsePolyZp.WellFormed_arr p (#[] : SparsePolyZp) := by
    intro x hx
    simp at hx
  rw [toPoly_foldl_mulStep p h_2p h_p2 g hg f.toList #[] h_empty_wf hf]
  rw [SparsePolyZp.toPoly_empty p, zero_add]
  -- 目标: listSum p f.toList * toPoly p g = toPoly p f * toPoly p g
  show listSum p f.toList * SparsePolyZp.toPoly p g = SparsePolyZp.toPoly p f * _
  rfl

-- Lemma 5: WellFormed_arr 在 * 下闭合
theorem SparsePolyZp.WellFormed_arr.mul (p : Nat)
    (h_p2 : p * p ≤ UInt64.size)
    (f g : SparsePolyZp)
    (hf : SparsePolyZp.WellFormed_arr p f) (hg : SparsePolyZp.WellFormed_arr p g) :
    SparsePolyZp.WellFormed_arr p (f * g) := by
  show SparsePolyZp.WellFormed_arr p (SparsePolyZp.mulImpl f g)
  unfold SparsePolyZp.mulImpl
  rw [← Array.foldl_toList]
  have h_empty_wf : SparsePolyZp.WellFormed_arr p (#[] : SparsePolyZp) := by
    intro x hx; simp at hx
  -- 内层 helper：xs 优先归纳，acc/hxs 留到 IH 内层
  suffices h : ∀ xs : List (UMonomial × Zp), ∀ acc : SparsePolyZp,
      SparsePolyZp.WellFormed_arr p acc → SparsePolyZp.AllReduced p xs →
      SparsePolyZp.WellFormed_arr p (List.foldl
        (fun a (t : UMonomial × Zp) =>
          SparsePolyZp.addImpl a (SparsePolyZp.scaleByMonomial t.fst t.snd g)) acc xs)
    from h f.toList #[] h_empty_wf hf
  intro xs
  induction xs with
  | nil =>
    intro acc hacc _; simp [List.foldl_nil]; exact hacc
  | cons x rest ih =>
    intro acc hacc hxs
    rcases x with ⟨mf, cf⟩
    have hcf_red : Zp.Reduced p cf := hxs (mf, cf) List.mem_cons_self
    have hxs' : SparsePolyZp.AllReduced p rest :=
      fun y hy => hxs y (List.mem_cons_of_mem _ hy)
    have h_scale_wf : SparsePolyZp.WellFormed_arr p
        (SparsePolyZp.scaleByMonomial mf cf g) :=
      SparsePolyZp.WellFormed_arr.scaleByMonomial p h_p2 mf cf hcf_red g hg
    have hacc' : SparsePolyZp.WellFormed_arr p
        (SparsePolyZp.addImpl acc (SparsePolyZp.scaleByMonomial mf cf g)) :=
      SparsePolyZp.WellFormed_arr.add p acc _ hacc h_scale_wf
    rw [List.foldl_cons]
    exact ih _ hacc' hxs'

-- ============================================================
-- §7g. 乘法 ring 公理（toPoly-level, Phase 2A.4c Stage 3a）
-- ============================================================
-- 通过 toPoly bridge 把 Mathlib Polynomial 的 ring 公理外包到 SparsePolyZp。
-- 注：这些是 toPoly-level 等式（toPoly p (f*g) = toPoly p (g*f) 等），
-- 不是 Array-level 等式。Array-level 等式需 toPoly_inj_canonical（Stage 3b）。

theorem SparsePolyZp.mul_comm_via_toPoly (p : Nat)
    (h_2p : 2 * p ≤ UInt64.size) (h_p2 : p * p ≤ UInt64.size)
    (f g : SparsePolyZp)
    (hf : SparsePolyZp.WellFormed_arr p f) (hg : SparsePolyZp.WellFormed_arr p g) :
    SparsePolyZp.toPoly p (f * g) = SparsePolyZp.toPoly p (g * f) := by
  rw [SparsePolyZp.toPoly_mul p h_2p h_p2 f g hf hg]
  rw [SparsePolyZp.toPoly_mul p h_2p h_p2 g f hg hf]
  ring

theorem SparsePolyZp.mul_assoc_via_toPoly (p : Nat)
    (h_2p : 2 * p ≤ UInt64.size) (h_p2 : p * p ≤ UInt64.size)
    (f g h : SparsePolyZp)
    (hf : SparsePolyZp.WellFormed_arr p f) (hg : SparsePolyZp.WellFormed_arr p g)
    (hh : SparsePolyZp.WellFormed_arr p h) :
    SparsePolyZp.toPoly p ((f * g) * h) = SparsePolyZp.toPoly p (f * (g * h)) := by
  have hfg := SparsePolyZp.WellFormed_arr.mul p h_p2 f g hf hg
  have hgh := SparsePolyZp.WellFormed_arr.mul p h_p2 g h hg hh
  rw [SparsePolyZp.toPoly_mul p h_2p h_p2 _ h hfg hh,
      SparsePolyZp.toPoly_mul p h_2p h_p2 f g hf hg,
      SparsePolyZp.toPoly_mul p h_2p h_p2 f _ hf hgh,
      SparsePolyZp.toPoly_mul p h_2p h_p2 g h hg hh]
  ring

theorem SparsePolyZp.left_distrib_via_toPoly (p : Nat)
    (h_2p : 2 * p ≤ UInt64.size) (h_p2 : p * p ≤ UInt64.size)
    (f g h : SparsePolyZp)
    (hf : SparsePolyZp.WellFormed_arr p f) (hg : SparsePolyZp.WellFormed_arr p g)
    (hh : SparsePolyZp.WellFormed_arr p h) :
    SparsePolyZp.toPoly p (f * (g + h)) =
    SparsePolyZp.toPoly p (f * g) + SparsePolyZp.toPoly p (f * h) := by
  have hgh := SparsePolyZp.WellFormed_arr.add p g h hg hh
  rw [SparsePolyZp.toPoly_mul p h_2p h_p2 f _ hf hgh,
      SparsePolyZp.toPoly_add p h_2p g h hg hh,
      SparsePolyZp.toPoly_mul p h_2p h_p2 f g hf hg,
      SparsePolyZp.toPoly_mul p h_2p h_p2 f h hf hh]
  ring

theorem SparsePolyZp.right_distrib_via_toPoly (p : Nat)
    (h_2p : 2 * p ≤ UInt64.size) (h_p2 : p * p ≤ UInt64.size)
    (f g h : SparsePolyZp)
    (hf : SparsePolyZp.WellFormed_arr p f) (hg : SparsePolyZp.WellFormed_arr p g)
    (hh : SparsePolyZp.WellFormed_arr p h) :
    SparsePolyZp.toPoly p ((f + g) * h) =
    SparsePolyZp.toPoly p (f * h) + SparsePolyZp.toPoly p (g * h) := by
  have hfg := SparsePolyZp.WellFormed_arr.add p f g hf hg
  rw [SparsePolyZp.toPoly_mul p h_2p h_p2 _ h hfg hh,
      SparsePolyZp.toPoly_add p h_2p f g hf hg,
      SparsePolyZp.toPoly_mul p h_2p h_p2 f h hf hh,
      SparsePolyZp.toPoly_mul p h_2p h_p2 g h hg hh]
  ring

-- 单位元：toPoly p (1 * f) = toPoly p f
-- 注：SparsePolyZp 没有显式 1 元素；常数 1 多项式 ⟨(0, 1), ...⟩ 需 Reduced
-- 此条留 Stage 3b（与 Canonical 一起）

-- ============================================================
-- §7h. Canonical 谓词 (Phase 2A.4c Stage 3b 起步)
-- ============================================================

-- Canonical: 严格 canonical 形式
-- (1) WellFormed_arr：所有 Zp 元素 .prime.toNat = p, .val.toNat < p
-- (2) Chain'：deg 严格降序（用 sorted 降序保 leading 项位置确定）
-- (3) 无零系数
def SparsePolyZp.Canonical (p : Nat) (f : SparsePolyZp) : Prop :=
  SparsePolyZp.WellFormed_arr p f ∧
  List.IsChain (fun a b : UMonomial × Zp => a.fst.deg > b.fst.deg) f.toList ∧
  ∀ x ∈ f.toList, x.snd.val ≠ 0

-- Helper: list 中所有项 deg < n → listSum p xs.coeff n = 0
theorem listSum_coeff_zero_of_all_lt (p : Nat) (n : Nat)
    (xs : List (UMonomial × Zp))
    (h : ∀ x ∈ xs, x.fst.deg < n) :
    (listSum p xs).coeff n = 0 := by
  induction xs with
  | nil => simp [listSum_nil]
  | cons x rest ih =>
    have hx : x.fst.deg < n := h x List.mem_cons_self
    have hrest : ∀ y ∈ rest, y.fst.deg < n :=
      fun y hy => h y (List.mem_cons_of_mem _ hy)
    rcases x with ⟨mx, cx⟩
    rw [listSum_cons, Polynomial.coeff_add, Polynomial.coeff_monomial]
    rw [if_neg (Nat.ne_of_lt hx)]
    rw [ih hrest]; ring

-- Helper: toZMod 在 Reduced 上 injective
theorem Zp.toZMod_inj_of_reduced (p : Nat) (a b : Zp)
    (ha : Zp.Reduced p a) (hb : Zp.Reduced p b)
    (heq : Zp.toZMod p a = Zp.toZMod p b) :
    a = b := by
  -- toZMod a = (a.val.toNat : ZMod p)；用 ZMod.val 反推
  have h_val_eq : a.val.toNat = b.val.toNat := by
    have ha_lt : a.val.toNat < p := ha.2
    have hb_lt : b.val.toNat < p := hb.2
    have ha_cast : ((a.val.toNat : ZMod p)).val = a.val.toNat := by
      rw [ZMod.val_natCast]; exact Nat.mod_eq_of_lt ha_lt
    have hb_cast : ((b.val.toNat : ZMod p)).val = b.val.toNat := by
      rw [ZMod.val_natCast]; exact Nat.mod_eq_of_lt hb_lt
    unfold Zp.toZMod at heq
    rw [← ha_cast, ← hb_cast, heq]
  have h_val : a.val = b.val := UInt64.toNat.inj h_val_eq
  have h_prime : a.prime = b.prime :=
    UInt64.toNat.inj (ha.1.trans hb.1.symm)
  rcases a with ⟨va, pa⟩
  rcases b with ⟨vb, pb⟩
  simp only at h_val h_prime
  rw [h_val, h_prime]

-- Helper: Chain' 严格降序 → head 之后所有 deg < head.deg
theorem chain_gt_all_after_head
    (head : UMonomial × Zp) (rest : List (UMonomial × Zp))
    (h_chain : List.IsChain (fun a b : UMonomial × Zp => a.fst.deg > b.fst.deg)
                 (head :: rest)) :
    ∀ x ∈ rest, x.fst.deg < head.fst.deg := by
  induction rest generalizing head with
  | nil =>
    intro x hx; simp at hx
  | cons r₀ rs ih =>
    intro x hx
    -- h_chain: Chain' R (head :: r₀ :: rs)
    rw [List.isChain_cons_cons] at h_chain
    obtain ⟨h_head_r₀, h_chain_r₀_rs⟩ := h_chain
    -- h_head_r₀ : head.fst.deg > r₀.fst.deg
    rcases List.mem_cons.mp hx with rfl | hx_in_rs
    · exact h_head_r₀
    · -- x ∈ rs；用 IH for (r₀ :: rs) on r₀
      have := ih r₀ h_chain_r₀_rs x hx_in_rs
      -- this : x.fst.deg < r₀.fst.deg；trans with h_head_r₀
      exact Nat.lt_trans this h_head_r₀

-- Helper: 排序列表 leading coeff
theorem listSum_coeff_at_head (p : Nat)
    (head : UMonomial × Zp) (rest : List (UMonomial × Zp))
    (h_chain : List.IsChain (fun a b : UMonomial × Zp => a.fst.deg > b.fst.deg)
                 (head :: rest)) :
    (listSum p (head :: rest)).coeff head.fst.deg = Zp.toZMod p head.snd := by
  rcases head with ⟨mh, ch⟩
  rw [listSum_cons, Polynomial.coeff_add, Polynomial.coeff_monomial]
  rw [if_pos rfl]
  have h_rest_lt : ∀ x ∈ rest, x.fst.deg < mh.deg :=
    chain_gt_all_after_head ⟨mh, ch⟩ rest h_chain
  rw [listSum_coeff_zero_of_all_lt p mh.deg rest h_rest_lt]; ring

-- Helper: Reduced p c + c.val ≠ 0 → toZMod p c ≠ 0
-- 因为 c.val < p 且 c.val ≠ 0 → (c.val.toNat : ZMod p) 不被 p 整除
theorem Zp.toZMod_ne_zero_of_val_ne_zero (p : Nat) (c : Zp)
    (hc : Zp.Reduced p c) (hval : c.val ≠ 0) :
    Zp.toZMod p c ≠ 0 := by
  intro heq
  -- toZMod p c = (c.val.toNat : ZMod p) = 0 → p ∣ c.val.toNat
  -- 但 c.val.toNat < p 且 c.val ≠ 0 → c.val.toNat ≠ 0 → ¬(p ∣ c.val.toNat)
  unfold Zp.toZMod at heq
  -- 用 ZMod.val: (c.val.toNat : ZMod p).val = c.val.toNat % p = c.val.toNat（因 < p）
  have h_cast : ((c.val.toNat : ZMod p)).val = c.val.toNat := by
    rw [ZMod.val_natCast]; exact Nat.mod_eq_of_lt hc.2
  rw [heq, ZMod.val_zero] at h_cast
  -- h_cast: 0 = c.val.toNat
  apply hval
  exact UInt64.toNat.inj h_cast.symm

-- List-level helper: canonical 列表（含 sorted、reduced、no-zero）的 listSum 单射
theorem listSum_inj_canonical (p : Nat)
    (xs ys : List (UMonomial × Zp))
    (hx_red : SparsePolyZp.AllReduced p xs)
    (hy_red : SparsePolyZp.AllReduced p ys)
    (hx_sorted : List.IsChain (fun a b : UMonomial × Zp => a.fst.deg > b.fst.deg) xs)
    (hy_sorted : List.IsChain (fun a b : UMonomial × Zp => a.fst.deg > b.fst.deg) ys)
    (hx_nz : ∀ x ∈ xs, x.snd.val ≠ 0)
    (hy_nz : ∀ y ∈ ys, y.snd.val ≠ 0)
    (heq : listSum p xs = listSum p ys) :
    xs = ys := by
  induction xs generalizing ys with
  | nil =>
    -- xs = []：listSum = 0；ys 必为 []
    rcases ys with _ | ⟨⟨mh, ch⟩, rest⟩
    · rfl
    · exfalso
      rw [listSum_nil] at heq
      have hch_red : Zp.Reduced p ch := hy_red (mh, ch) List.mem_cons_self
      have hch_nz : ch.val ≠ 0 := hy_nz (mh, ch) List.mem_cons_self
      have h_ne : Zp.toZMod p ch ≠ 0 :=
        Zp.toZMod_ne_zero_of_val_ne_zero p ch hch_red hch_nz
      have h_coeff : (listSum p ((mh, ch) :: rest)).coeff mh.deg = Zp.toZMod p ch :=
        listSum_coeff_at_head p (mh, ch) rest hy_sorted
      rw [← heq] at h_coeff
      simp at h_coeff
      exact h_ne h_coeff.symm
  | cons xf rest_f ih =>
    rcases xf with ⟨mf, cf⟩
    rcases ys with _ | ⟨⟨mg, cg⟩, rest_g⟩
    · -- ys = []：对称矛盾
      exfalso
      rw [listSum_nil] at heq
      have hcf_red : Zp.Reduced p cf := hx_red (mf, cf) List.mem_cons_self
      have hcf_nz : cf.val ≠ 0 := hx_nz (mf, cf) List.mem_cons_self
      have h_ne : Zp.toZMod p cf ≠ 0 :=
        Zp.toZMod_ne_zero_of_val_ne_zero p cf hcf_red hcf_nz
      have h_coeff : (listSum p ((mf, cf) :: rest_f)).coeff mf.deg = Zp.toZMod p cf :=
        listSum_coeff_at_head p (mf, cf) rest_f hx_sorted
      rw [heq] at h_coeff
      simp at h_coeff
      exact h_ne h_coeff.symm
    · -- 两边非空：比 leading deg + leading coef
      have hcf_red : Zp.Reduced p cf := hx_red (mf, cf) List.mem_cons_self
      have hcg_red : Zp.Reduced p cg := hy_red (mg, cg) List.mem_cons_self
      have hcf_nz : cf.val ≠ 0 := hx_nz (mf, cf) List.mem_cons_self
      have hcg_nz : cg.val ≠ 0 := hy_nz (mg, cg) List.mem_cons_self
      have hf_ne : Zp.toZMod p cf ≠ 0 :=
        Zp.toZMod_ne_zero_of_val_ne_zero p cf hcf_red hcf_nz
      have hg_ne : Zp.toZMod p cg ≠ 0 :=
        Zp.toZMod_ne_zero_of_val_ne_zero p cg hcg_red hcg_nz
      have h_deg_eq : mf.deg = mg.deg := by
        rcases lt_trichotomy mf.deg mg.deg with hlt | heq_d | hgt
        · exfalso
          have h_f_zero : (listSum p ((mf, cf) :: rest_f)).coeff mg.deg = 0 := by
            apply listSum_coeff_zero_of_all_lt
            intro x hx
            rcases List.mem_cons.mp hx with rfl | hx_in_rest
            · exact hlt
            · have := chain_gt_all_after_head (mf, cf) rest_f hx_sorted x hx_in_rest
              exact Nat.lt_trans this hlt
          have h_g_lead :
              (listSum p ((mg, cg) :: rest_g)).coeff mg.deg = Zp.toZMod p cg :=
            listSum_coeff_at_head p (mg, cg) rest_g hy_sorted
          rw [heq] at h_f_zero
          rw [h_g_lead] at h_f_zero
          exact hg_ne h_f_zero
        · exact heq_d
        · exfalso
          have h_g_zero : (listSum p ((mg, cg) :: rest_g)).coeff mf.deg = 0 := by
            apply listSum_coeff_zero_of_all_lt
            intro x hx
            rcases List.mem_cons.mp hx with rfl | hx_in_rest
            · exact hgt
            · have := chain_gt_all_after_head (mg, cg) rest_g hy_sorted x hx_in_rest
              exact Nat.lt_trans this hgt
          have h_f_lead :
              (listSum p ((mf, cf) :: rest_f)).coeff mf.deg = Zp.toZMod p cf :=
            listSum_coeff_at_head p (mf, cf) rest_f hx_sorted
          rw [← heq] at h_g_zero
          rw [h_f_lead] at h_g_zero
          exact hf_ne h_g_zero
      have h_cf_eq_cg : cf = cg := by
        have h_f_lead : (listSum p ((mf, cf) :: rest_f)).coeff mf.deg = Zp.toZMod p cf :=
          listSum_coeff_at_head p (mf, cf) rest_f hx_sorted
        have h_g_lead : (listSum p ((mg, cg) :: rest_g)).coeff mg.deg = Zp.toZMod p cg :=
          listSum_coeff_at_head p (mg, cg) rest_g hy_sorted
        rw [heq, h_deg_eq] at h_f_lead
        rw [h_g_lead] at h_f_lead
        exact Zp.toZMod_inj_of_reduced p cf cg hcf_red hcg_red h_f_lead.symm
      have h_mf_eq_mg : mf = mg := by
        rcases mf with ⟨df⟩; rcases mg with ⟨dg⟩
        simp at h_deg_eq
        exact congrArg UMonomial.mk h_deg_eq
      have h_rest_listSum_eq : listSum p rest_f = listSum p rest_g := by
        have heq' := heq
        rw [listSum_cons, listSum_cons] at heq'
        rw [h_mf_eq_mg, h_cf_eq_cg] at heq'
        exact add_left_cancel heq'
      have hx_red' : SparsePolyZp.AllReduced p rest_f :=
        fun y hy => hx_red y (List.mem_cons_of_mem _ hy)
      have hy_red' : SparsePolyZp.AllReduced p rest_g :=
        fun y hy => hy_red y (List.mem_cons_of_mem _ hy)
      have hx_sorted' : List.IsChain (fun a b : UMonomial × Zp => a.fst.deg > b.fst.deg) rest_f :=
        (List.isChain_cons.mp hx_sorted).2
      have hy_sorted' : List.IsChain (fun a b : UMonomial × Zp => a.fst.deg > b.fst.deg) rest_g :=
        (List.isChain_cons.mp hy_sorted).2
      have hx_nz' : ∀ x ∈ rest_f, x.snd.val ≠ 0 :=
        fun y hy => hx_nz y (List.mem_cons_of_mem _ hy)
      have hy_nz' : ∀ y ∈ rest_g, y.snd.val ≠ 0 :=
        fun y hy => hy_nz y (List.mem_cons_of_mem _ hy)
      have h_rest_eq := ih rest_g hx_red' hy_red' hx_sorted' hy_sorted'
                        hx_nz' hy_nz' h_rest_listSum_eq
      rw [h_mf_eq_mg, h_cf_eq_cg, h_rest_eq]

-- 核心定理：toPoly_inj_canonical (Array level)
theorem SparsePolyZp.toPoly_inj_canonical (p : Nat)
    (f g : SparsePolyZp)
    (hf : SparsePolyZp.Canonical p f) (hg : SparsePolyZp.Canonical p g)
    (heq : SparsePolyZp.toPoly p f = SparsePolyZp.toPoly p g) :
    f = g := by
  suffices h : f.toList = g.toList by
    cases f; cases g; congr 1
  obtain ⟨hf_wf, hf_sorted, hf_nonzero⟩ := hf
  obtain ⟨hg_wf, hg_sorted, hg_nonzero⟩ := hg
  unfold SparsePolyZp.toPoly at heq
  exact listSum_inj_canonical p _ _ hf_wf hg_wf hf_sorted hg_sorted
    hf_nonzero hg_nonzero heq

-- ============================================================
-- §8. 数值验证（Mathlib decidable 通过）
-- ============================================================
--
-- 余下工作（Phase 2A.4c Stage 3b 续）：
-- - listSum_coeff_zero_of_all_lt, listSum_coeff_at_head (coeff helpers)
-- - Zp.toZMod_inj_of_reduced (toZMod 单射)
-- - toPoly_inj_canonical（核心难点）
-- - Canonical.add / Canonical.mul 保持性
-- - Array-level ring 公理（f * g = g * f 等，via toPoly_inj）
-- ============================================================

-- Zp 7 的 0 → ZMod 7 的 0
example : Zp.toZMod 7 ⟨0, 7⟩ = 0 := by decide

-- Zp 7 的 3 → ZMod 7 的 3
example : Zp.toZMod 7 ⟨3, 7⟩ = (3 : ZMod 7) := by decide

-- empty SparsePolyZp 是 WellFormed
example : SparsePolyZp.WellFormed 7 #[] := SparsePolyZp.WellFormed.empty 7

end CLPoly.Math
