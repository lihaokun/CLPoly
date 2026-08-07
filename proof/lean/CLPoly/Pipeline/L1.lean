/-
  CLPoly/Pipeline/L1.lean — L1 翻译代码的管线实例化

  用 L1 `_ir` 函数（来自 Corpus.lean，翻译自 C++）实例化 factor_Zp_correct。
  精化定理填补后，此文件直接验证翻译出的 C++ 代码。

  当前不导出算法级 L1 正确性；各阶段只有在严格 generated→L2
  证明闭合后才能进入此文件。
-/

import CLPoly.Pipeline.FactorZp
import CLPoly.Pipeline.FactorZZ
import CLPoly.Pipeline.FactorZZInstantiate
import CLPoly.Algorithm.EDF
import CLPoly.Algorithm.SquarefreeZp
import CLPoly.Algorithm.DDF
import CLPoly.Refinement.Basic
import CLPoly.Refinement.SquarefreeZp
import CLPoly.Refinement.DDF
import CLPoly.Refinement.EDF
import CLPoly.Refinement.Hensel
import CLPoly.Refinement.Recombine
import CLPoly.Algorithm.Hensel
import CLPoly.Algorithm.Recombine
import Mathlib.Algebra.Polynomial.Degree.Support

set_option autoImplicit false

open Polynomial
open CLPoly.Math
open Refinement
open Classical

noncomputable section L1Pipeline

variable {p : ℕ} [hp : Fact (Nat.Prime p)]

-- ============================================================
-- §1. Bridge: Polynomial (ZMod p) ←→ SparsePolyZp
-- ============================================================

/-- 将 Polynomial (ZMod p) 转换为 SparsePolyZp（Array 表示）。
    利用 Polynomial.support（非零系数的度集合）构造项列表。
    noncomputable 因为 Polynomial 基于 Finsupp（非可计算）。 -/
noncomputable def toSparsePolyZp (f : (ZMod p)[X]) : SparsePolyZp :=
  (f.support.sort (· ≥ ·)).map (fun n =>
    ({deg := n}, Zp.ofUInt64 (UInt64.ofNat (ZMod.val (f.coeff n))) (UInt64.ofNat p))
  ) |>.toArray

/-- Integer-polynomial counterpart used by the generated Hensel entry point. -/
noncomputable def toSparsePolyZZ (f : Polynomial ℤ) : SparsePolyZZ :=
  (f.support.sort (· ≥ ·)).map (fun n =>
    ({deg := n}, f.coeff n)
  ) |>.toArray

private lemma sortedListB_of_pairwise : ∀ l : List (UMonomial × Zp),
    l.Pairwise (fun a b => b.1.deg < a.1.deg) → SparsePolyZp.sortedListB l = true := by
  intro l
  induction l with
  | nil => simp [SparsePolyZp.sortedListB]
  | cons a rest ih =>
      intro hs
      rw [SparsePolyZp.sortedListB_iff]
      exact ⟨(List.pairwise_cons.mp hs).1, ih (List.pairwise_cons.mp hs).2⟩

private lemma pairwise_gt_of_pairwise_ge_nodup : ∀ l : List Nat,
    l.Pairwise (· ≥ ·) → l.Nodup → l.Pairwise (· > ·) := by
  intro l
  induction l with
  | nil => simp
  | cons a rest ih =>
      intro hge hnd
      rw [List.pairwise_cons] at hge ⊢
      rw [List.nodup_cons] at hnd
      refine ⟨?_, ih hge.2 hnd.2⟩
      intro b hb
      have hab := hge.1 b hb
      have hne : a ≠ b := fun he => hnd.1 (he ▸ hb)
      omega

lemma toSparsePolyZp_sorted (f : (ZMod p)[X]) :
    SparsePolyZp.Sorted (toSparsePolyZp f) := by
  unfold toSparsePolyZp SparsePolyZp.Sorted
  simp only [List.toList_toArray]
  apply sortedListB_of_pairwise
  rw [List.pairwise_map]
  simpa only using pairwise_gt_of_pairwise_ge_nodup (f.support.sort (· ≥ ·))
    (Finset.sort_sorted f.support (· ≥ ·)) (Finset.sort_nodup f.support (· ≥ ·))

lemma toSparsePolyZp_wellFormed (f : (ZMod p)[X]) (hp_size : 2 * p ≤ UInt64.size) :
    SparsePolyZp.WellFormed p (toSparsePolyZp f) := by
  unfold toSparsePolyZp SparsePolyZp.WellFormed
  intro x hx
  have hx_list : x ∈ ((f.support.sort (· ≥ ·)).map (fun n =>
      ({deg := n}, Zp.ofUInt64 (UInt64.ofNat (ZMod.val (f.coeff n))) (UInt64.ofNat p)))) :=
    Array.mem_toArray.mp hx
  rcases List.mem_map.mp hx_list with ⟨n, hn, rfl⟩
  have hp_prime : Nat.Prime p := hp.out
  have hp_pos : 0 < p := Nat.Prime.pos hp_prime
  have hp_lt : p < UInt64.size := by nlinarith
  unfold Zp.ofUInt64
  simp [hp_lt]

lemma toSparsePolyZp_allReduced (f : (ZMod p)[X]) (hp_size : 2 * p ≤ UInt64.size) :
    SparsePolyZp.AllReduced p (toSparsePolyZp f).toList := by
  unfold toSparsePolyZp
  intro x hx
  rcases List.mem_map.mp hx with ⟨n, hn, rfl⟩
  have hp_prime : Nat.Prime p := hp.out
  have hp_pos : 0 < p := Nat.Prime.pos hp_prime
  have hp_lt : p < UInt64.size := by nlinarith
  have hval_lt_p : ZMod.val (f.coeff n) < p := ZMod.val_lt (f.coeff n)
  have hval_lt_64 : ZMod.val (f.coeff n) < UInt64.size := lt_trans hval_lt_p hp_lt
  unfold Zp.ofUInt64
  refine ⟨by simp [hp_lt], ?_⟩
  have h1 : (UInt64.ofNat (ZMod.val (f.coeff n))).toNat = ZMod.val (f.coeff n) := by simp [hval_lt_64]
  have h2 : (UInt64.ofNat p).toNat = p := by simp [hp_lt]
  have h_mod : (UInt64.ofNat (ZMod.val (f.coeff n)) % UInt64.ofNat p).toNat = ZMod.val (f.coeff n) % p := by
    simp [h1, h2]
  have h_mod_eq : ZMod.val (f.coeff n) % p = ZMod.val (f.coeff n) := Nat.mod_eq_of_lt hval_lt_p
  simp [h_mod, h_mod_eq, hval_lt_p]

lemma toSparsePolyZp_val_nonzero (f : (ZMod p)[X]) (hp_size : 2 * p ≤ UInt64.size) :
    ∀ x ∈ (toSparsePolyZp f).toList, x.snd.val.toNat ≠ 0 := by
  unfold toSparsePolyZp
  intro x hx
  -- x is in the array built from f.support; all support elements have non-zero coefficients
  have hx_list : x ∈ ((f.support.sort (· ≥ ·)).map (fun n =>
    ({deg := n}, Zp.ofUInt64 (UInt64.ofNat (ZMod.val (f.coeff n))) (UInt64.ofNat p)))).toArray.toList := hx
  rw [List.toList_toArray] at hx_list
  rcases List.mem_map.mp hx_list with ⟨n, hn, rfl⟩
  -- n ∈ f.support, so f.coeff n ≠ 0, so ZMod.val (f.coeff n) ≠ 0
  have hn' : n ∈ f.support := Finset.mem_sort (· ≥ ·) |>.mp hn
  have h_coeff_ne_zero : f.coeff n ≠ 0 := Polynomial.mem_support_iff.mp hn'
  have h_val_ne_zero : ZMod.val (f.coeff n) ≠ 0 :=
    mt (ZMod.val_eq_zero _).mp h_coeff_ne_zero
  -- Zp.ofUInt64 (UInt64.ofNat v) p has val = v % p (UInt64)；由 v < p < 2^64 化简为 v ≠ 0
  have hp_prime : Nat.Prime p := hp.out
  haveI : NeZero p := ⟨hp_prime.ne_zero⟩
  have hv_lt : ZMod.val (f.coeff n) < p := ZMod.val_lt (f.coeff n)
  have hp_lt_sz : p < UInt64.size := by
    have h2 : (2 : ℕ) ≤ p := hp_prime.two_le
    omega
  have hv_lt_sz : ZMod.val (f.coeff n) < UInt64.size := lt_trans hv_lt hp_lt_sz
  unfold Zp.ofUInt64
  simp only [UInt64.toNat_mod]
  have h1 : (UInt64.ofNat (ZMod.val (f.coeff n))).toNat = ZMod.val (f.coeff n) := by simp [hv_lt_sz]
  have h2 : (UInt64.ofNat p).toNat = p := by simp [hp_lt_sz]
  rw [h1, h2, Nat.mod_eq_of_lt hv_lt]
  exact h_val_ne_zero

lemma toSparsePolyZp_toPoly (f : (ZMod p)[X]) (hp_size : 2 * p ≤ UInt64.size) :
    SparsePolyZp.toPoly p (toSparsePolyZp f) = f := by
  unfold toSparsePolyZp SparsePolyZp.toPoly
  have hp_prime : Nat.Prime p := hp.out
  have hp_pos : 0 < p := Nat.Prime.pos hp_prime
  have hp_lt : p < UInt64.size := by nlinarith
  have hval_lt (n : ℕ) : ZMod.val (f.coeff n) < UInt64.size := by
    by_cases hn : f.coeff n = 0
    · have hzero : ZMod.val (f.coeff n) = 0 := by simpa [hn] using ZMod.val_zero
      have hpos : 0 < UInt64.size := by decide
      simpa [hzero]
    · have h_support : n ∈ f.support := Polynomial.mem_support_iff.mpr hn
      have h_lt_p : ZMod.val (f.coeff n) < p := ZMod.val_lt (f.coeff n)
      exact lt_trans h_lt_p hp_lt
  have h_coeff (n : ℕ) : Zp.toZMod p (Zp.ofUInt64 (UInt64.ofNat (ZMod.val (f.coeff n))) (UInt64.ofNat p)) =
      f.coeff n := by
    unfold Zp.toZMod Zp.ofUInt64
    have hp_lt' : (UInt64.ofNat p).toNat = p := by simp [hp_lt]
    have h_mod_toNat : (UInt64.ofNat (ZMod.val (f.coeff n)) % UInt64.ofNat p).toNat = ZMod.val (f.coeff n) := by
      have h1 : (UInt64.ofNat (ZMod.val (f.coeff n))).toNat = ZMod.val (f.coeff n) := by simp [hval_lt n]
      calc
        (UInt64.ofNat (ZMod.val (f.coeff n)) % UInt64.ofNat p).toNat
            = (UInt64.ofNat (ZMod.val (f.coeff n))).toNat % (UInt64.ofNat p).toNat := by simp
        _ = ZMod.val (f.coeff n) % p := by simp [h1, hp_lt']
        _ = ZMod.val (f.coeff n) := Nat.mod_eq_of_lt (ZMod.val_lt (f.coeff n))
    simp [h_mod_toNat, ZMod.natCast_zmod_val]
  have h_nodup : (f.support.sort (· ≥ ·)).Nodup :=
    Finset.sort_nodup (s := f.support) (r := (· ≥ ·))
  have hlistSum (l : List ℕ) (hn : l.Nodup) : listSum p (l.map (fun n =>
      ({deg := n}, Zp.ofUInt64 (UInt64.ofNat (ZMod.val (f.coeff n))) (UInt64.ofNat p)))) =
      ∑ n ∈ l.toFinset, monomial n (f.coeff n) := by
    induction l with
    | nil => simp
    | cons x xs ih =>
      cases hn with
      | cons hx_not_mem hxs_nodup =>
        have hx_not_mem' : x ∉ xs := λ hx_mem => hx_not_mem x hx_mem rfl
        have hx_not_mem_fs : x ∉ xs.toFinset :=
          (mt ((List.mem_toFinset (l := xs)).mp)) hx_not_mem'
        calc
          listSum p (((x :: xs).map (fun n =>
              ({deg := n}, Zp.ofUInt64 (UInt64.ofNat (ZMod.val (f.coeff n))) (UInt64.ofNat p))))
            ) = monomial x (f.coeff x) + listSum p (xs.map (fun n =>
              ({deg := n}, Zp.ofUInt64 (UInt64.ofNat (ZMod.val (f.coeff n))) (UInt64.ofNat p)))) := by
            simp [listSum, h_coeff x]
          _ = monomial x (f.coeff x) + ∑ n ∈ xs.toFinset, monomial n (f.coeff n) := by
            rw [ih hxs_nodup]
          _ = ∑ n ∈ insert x xs.toFinset, monomial n (f.coeff n) := by
            rw [Finset.sum_insert hx_not_mem_fs]
          _ = ∑ n ∈ (x :: xs).toFinset, monomial n (f.coeff n) := by simp
  calc
    listSum p ((f.support.sort (· ≥ ·)).map (fun n =>
        ({deg := n}, Zp.ofUInt64 (UInt64.ofNat (ZMod.val (f.coeff n))) (UInt64.ofNat p)))) =
      ∑ n ∈ ((f.support.sort (· ≥ ·)).toFinset), monomial n (f.coeff n) :=
      hlistSum (f.support.sort (· ≥ ·)) h_nodup
    _ = ∑ n ∈ f.support, monomial n (f.coeff n) := by simp
    _ = f := (Polynomial.as_sum_support f).symm

/-- toSparsePolyZp 各项度数不超过 natDegree（支集元素 ≤ natDegree）。 -/
lemma toSparsePolyZp_deg_le (f : (ZMod p)[X]) :
    ∀ x ∈ (toSparsePolyZp f).toList, x.1.deg ≤ f.natDegree := by
  intro x hx
  unfold toSparsePolyZp at hx
  simp only [List.mem_map] at hx
  obtain ⟨n, hn_mem, hn_eq⟩ := hx
  rw [← hn_eq]
  exact Polynomial.le_natDegree_of_mem_supp n ((Finset.mem_sort _).mp hn_mem)

-- No L1 algorithm wrapper is exported here until it targets a total,
-- cpp2lean-generated execution tree with a direct refinement theorem.

/- DDF, EDF, Hensel, recombination, and end-to-end L1 exports are intentionally
omitted until their cpp2lean-generated entries have direct refinement proofs.
No L2 algorithm or existence witness is used as an implementation fallback. -/

end L1Pipeline
