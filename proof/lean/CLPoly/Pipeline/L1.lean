/-
  CLPoly/Pipeline/L1.lean — L1 翻译代码的管线实例化

  用 L1 `_ir` 函数（来自 Corpus.lean，翻译自 C++）实例化 factor_Zp_correct。
  精化定理填补后，此文件直接验证翻译出的 C++ 代码。

  当前已闭合的 L1 覆盖：SQF (`__squarefree_Zp_ir`)
  DDF/EDF 仍有独立精化债务；EDF 暂用 L2 `edf_correct_unconditional`。
-/

import CLPoly.Pipeline.FactorZp
import CLPoly.Pipeline.FactorZZ
import CLPoly.Pipeline.FactorZZInstantiate
import CLPoly.Algorithm.EDF
import CLPoly.Algorithm.SquarefreeZp
import CLPoly.Algorithm.DDF
import CLPoly.Generated.Corpus
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

-- ============================================================
-- §2. L1 Wrappers
-- ============================================================

/-- L1 sqf wrapper — 使用 `__squarefree_Zp_ir_safe`（含常数多项式保护）。
     正确性由 `__squarefree_Zp_ir_refines` 保证。 -/
noncomputable def sqfZp_l1 (hp_size : 2 * p ≤ UInt64.size) (f : (ZMod p)[X])
    : List ((ZMod p)[X] × ℕ) :=
  toPolyList (__squarefree_Zp_ir_safe p (toSparsePolyZp f)) p

theorem sqfZp_l1_correct (hp_size : 2 * p ≤ UInt64.size) (hp2 : p * p ≤ UInt64.size)
    (f : (ZMod p)[X]) (hf : f ≠ 0) (hdeg : f.natDegree * p < 2 ^ 64) :
    SquarefreeDecomp f (sqfZp_l1 hp_size f) := by
  unfold sqfZp_l1
  have hwf : SparsePolyZp.WellFormed p (toSparsePolyZp f) :=
    toSparsePolyZp_wellFormed f hp_size
  have hred : SparsePolyZp.AllReduced p (toSparsePolyZp f).toList :=
    toSparsePolyZp_allReduced f hp_size
  have hppos : 0 < p := (Fact.out (p := Nat.Prime p)).pos
  -- 系数 val < p、度数 ≤ natDegree ⇒ 乘积 ≤ natDegree·p < 2^64（无 UInt64 溢出）。
  have h_no_overflow : ∀ x ∈ (toSparsePolyZp f).toList, x.2.val.toNat * x.1.deg < 2 ^ 64 := by
    intro x hx
    have hval : x.2.val.toNat < p := (hred x hx).2
    have hdle : x.1.deg ≤ f.natDegree := toSparsePolyZp_deg_le f x hx
    have hle : x.2.val.toNat * x.1.deg ≤ f.natDegree * p := by
      rw [Nat.mul_comm f.natDegree p]; exact Nat.mul_le_mul (Nat.le_of_lt hval) hdle
    exact Nat.lt_of_le_of_lt hle hdeg
  have h_deg_bound : ∀ x ∈ (toSparsePolyZp f).toList, x.1.deg < 2 ^ 64 := by
    intro x hx
    have hdle : x.1.deg ≤ f.natDegree := toSparsePolyZp_deg_le f x hx
    have hnp : f.natDegree ≤ f.natDegree * p := Nat.le_mul_of_pos_right _ hppos
    omega
  have h_p_deg_bound : ∀ x ∈ (toSparsePolyZp f).toList,
      p * x.1.deg < 2 ^ 64 := by
    intro x hx
    have hdle := toSparsePolyZp_deg_le f x hx
    have hle : p * x.1.deg ≤ f.natDegree * p := by
      rw [Nat.mul_comm p x.1.deg]
      exact Nat.mul_le_mul_right p hdle
    exact lt_of_le_of_lt hle hdeg
  have hf_signed : f.natDegree < 2 ^ 63 := by
    have hp_ge_two : 2 ≤ p := (Fact.out (p := Nat.Prime p)).two_le
    have htwop : f.natDegree * 2 ≤ f.natDegree * p :=
      Nat.mul_le_mul_left f.natDegree hp_ge_two
    have hlt : f.natDegree * 2 < 2 ^ 64 := lt_of_le_of_lt htwop hdeg
    norm_num [pow_succ] at hlt ⊢
    omega
  have h_signed_deg_bound : ∀ x ∈ (toSparsePolyZp f).toList, x.1.deg < 2 ^ 63 := by
    intro x hx
    exact lt_of_le_of_lt (toSparsePolyZp_deg_le f x hx) hf_signed
  have h_refines : toPolyList (__squarefree_Zp_ir_safe p (toSparsePolyZp f)) p = sqfZp f := by
    have h := __squarefree_Zp_ir_refines p (toSparsePolyZp f) (toSparsePolyZp_sorted f)
      hwf hred hp_size hp2 h_no_overflow h_p_deg_bound h_deg_bound h_signed_deg_bound
      (toSparsePolyZp_val_nonzero f hp_size)
    rw [toSparsePolyZp_toPoly f hp_size] at h
    exact h
  rw [h_refines]
  exact sqf_correct f hf

/-- L1 ddf wrapper — 使用 `__ddf_Zp_ir`（Corpus.lean，翻译自 C++）。
    正确性由 `__ddf_Zp_ir_refines` 保证。 -/
noncomputable def ddf_l1 (hp_size : 2 * p ≤ UInt64.size) (f : (ZMod p)[X])
    : List ((ZMod p)[X] × ℕ) :=
  toPolyList (Generated.__ddf_Zp_ir (toSparsePolyZp f)) p

theorem ddf_l1_correct (hp_size : 2 * p ≤ UInt64.size) (f : (ZMod p)[X]) (hm : Monic f)
    (hsq : Squarefree f) : DDFCorrect f (ddf_l1 hp_size f) := by
  unfold ddf_l1
  have hwf : SparsePolyZp.WellFormed p (toSparsePolyZp f) :=
    toSparsePolyZp_wellFormed f hp_size
  have hred : SparsePolyZp.AllReduced p (toSparsePolyZp f).toList :=
    toSparsePolyZp_allReduced f hp_size
  have h_refines : toPolyList (Generated.__ddf_Zp_ir (toSparsePolyZp f)) p = ddf f := by
    have h := __ddf_Zp_ir_refines p (toSparsePolyZp f) hwf hred hp_size
    simpa [toSparsePolyZp_toPoly f hp_size] using h
  rw [h_refines]
  exact ddf_correct f hm hsq

-- ============================================================
-- §3. EDF（使用 __edf_Zp_ir，C++ 翻译）
--     正确性由 `__edf_Zp_ir_refines` + `edf_correct_unconditional` 保证。
-- ============================================================

/-- L1 edf wrapper — 使用 `__edf_Zp_ir`（Corpus.lean，翻译自 C++）。
    当 `__edf_Zp_ir_refines` 填补后，将比对 C++ 输出与 L2 `edf` 的一致性。 -/
noncomputable def edf_l1 (hp_size : 2 * p ≤ UInt64.size) (g : (ZMod p)[X]) (d : ℕ)
    : List ((ZMod p)[X]) :=
  if hg_deg : g.natDegree = 0 then []
  else
    if hpre : Monic g ∧ Squarefree g ∧ 0 < d ∧
        (∀ q : (ZMod p)[X], Irreducible q → q ∣ g → q.natDegree = d) then
      (edf_correct_unconditional g d hpre.1 hpre.2.1
        (Nat.pos_of_ne_zero hg_deg) hpre.2.2.1 hpre.2.2.2).choose
    else [g]

theorem edf_l1_correct (hp_size : 2 * p ≤ UInt64.size) (g : (ZMod p)[X]) (d : ℕ)
    (hm : Monic g) (hsq : Squarefree g)
    (hdeg : ∀ q, Irreducible q → q ∣ g → q.natDegree = d) :
    EDFCorrect g d (edf_l1 hp_size g d) := by
  by_cases hg_deg : g.natDegree = 0
  · have htarget : edf_l1 hp_size g d = [] := by
      unfold edf_l1
      rw [dif_pos hg_deg]
    rw [htarget]
    have hg_eq_one : g = 1 := by
      have hg_eq := Polynomial.eq_C_of_natDegree_eq_zero hg_deg
      have h_lc : g.coeff 0 = 1 := by
        have := hm.leadingCoeff
        rw [Polynomial.leadingCoeff, hg_deg] at this
        exact this
      rw [hg_eq, h_lc, map_one]
    subst hg_eq_one
    exact ⟨Associated.refl 1, by intro q hq; simp at hq⟩
  · have hg_pos : 0 < g.natDegree := Nat.pos_of_ne_zero hg_deg
    have hg_nu : ¬IsUnit g := fun hu =>
      absurd (natDegree_eq_zero_of_isUnit hu) hg_deg
    obtain ⟨q₀, hq₀_irr, hq₀_dvd⟩ :=
      WfDvdMonoid.exists_irreducible_factor hg_nu (Monic.ne_zero hm)
    have hd_pos : 0 < d := by
      rw [← hdeg q₀ hq₀_irr hq₀_dvd]
      have hq_nu : ¬ IsUnit q₀ := hq₀_irr.1
      by_contra hzero
      have hzero' : q₀.natDegree = 0 := Nat.eq_zero_of_not_pos hzero
      have hconst : q₀ = Polynomial.C (q₀.coeff 0) :=
        Polynomial.eq_C_of_natDegree_eq_zero hzero'
      have hc_ne_zero : q₀.coeff 0 ≠ 0 := by
        intro hc0
        apply hq₀_irr.ne_zero
        rw [hconst, hc0, Polynomial.C_0]
      have h_unit : IsUnit q₀ := by
        rw [hconst]
        refine Polynomial.isUnit_C.mpr ?_
        refine ⟨⟨q₀.coeff 0, (q₀.coeff 0)⁻¹, ?_, ?_⟩, rfl⟩
        · field_simp [hc_ne_zero]
        · field_simp [hc_ne_zero]
      exact hq_nu h_unit
    have htarget : edf_l1 hp_size g d = (edf_correct_unconditional g d hm hsq hg_pos hd_pos hdeg).choose := by
      unfold edf_l1
      rw [dif_neg hg_deg]
      rw [dif_pos ⟨hm, hsq, hd_pos, hdeg⟩]
    rw [htarget]
    exact (edf_correct_unconditional g d hm hsq hg_pos hd_pos hdeg).choose_spec

-- ============================================================
-- §4. Zp 因式分解端到端定理（L1 包装）
-- ============================================================

/-- 使用 L1 翻译代码的 Zp 因式分解（需要 hp_size 硬件约束）。 -/
theorem factor_Zp_l1 (hp_size : 2 * p ≤ UInt64.size) (hp2 : p * p ≤ UInt64.size) (f : (ZMod p)[X]) (hf : f ≠ 0)
    (hdeg : f.natDegree * p < 2 ^ 64) :
    ∃ (lc : ZMod p) (factors : List ((ZMod p)[X] × ℕ)),
      FactorZpCorrect f lc factors :=
  factor_Zp_correct f hf
    (sqfZp_l1 hp_size) (sqfZp_l1_correct hp_size hp2 f hf hdeg)
    (ddf_l1 hp_size) (ddf_l1_correct hp_size)
    (edf_l1 hp_size) (edf_l1_correct hp_size)

/-- L1 Zp 因式分解函数 — 适配 factor_ZZ_correct 的接口（返回 `(lc, factors)` 对）。
    大度数（溢出）分支返回占位 (0, [])；正确性仅在 `g.natDegree * p < 2^64` 时成立。 -/
noncomputable def factor_zp_l1_func (hp_size : 2 * p ≤ UInt64.size) (hp2 : p * p ≤ UInt64.size) (g : Polynomial (ZMod p))
    : ZMod p × List ((ZMod p)[X] × ℕ) :=
  if hg : g = 0 then (0, [])
  else if hb : g.natDegree * p < 2 ^ 64 then
    have h := factor_Zp_l1 hp_size hp2 g hg hb
    (h.choose, h.choose_spec.choose)
  else (0, [])

lemma factor_zp_l1_func_correct (hp_size : 2 * p ≤ UInt64.size) (hp2 : p * p ≤ UInt64.size)
    (g : Polynomial (ZMod p)) (hg : g ≠ 0) (hb : g.natDegree * p < 2 ^ 64) :
    FactorZpCorrect g (factor_zp_l1_func hp_size hp2 g).1 (factor_zp_l1_func hp_size hp2 g).2 := by
  unfold factor_zp_l1_func
  rw [dif_neg hg, dif_pos hb]
  exact (factor_Zp_l1 hp_size hp2 g hg hb).choose_spec.choose_spec

-- ============================================================
-- §5. Hensel 提升（L1 包装 — TODO: 替换为 __hensel_lift_upoly_ir 精化版本）
-- ============================================================

/-- L1 Hensel 提升 — 使用 `__hensel_lift_upoly_ir`（Corpus.lean，翻译自 C++）。
    正确性由 `__hensel_lift_upoly_ir_refines` 保证。
    当前暂用 L2 `hensel_lift`（待精化定理填补后改回 C++ 路径）。 -/
noncomputable def hensel_l1 (hp_size : 2 * p ≤ UInt64.size) (k : ℕ) (hk : 0 < k) (f : Polynomial ℤ)
    (facs_p : List (Polynomial (ZMod p))) : List (Polynomial (ZMod (p ^ k))) :=
  hensel_lift p k hk f facs_p

lemma hensel_l1_correct (hp_size : 2 * p ≤ UInt64.size) (k : ℕ) (hk : 0 < k) (f : Polynomial ℤ)
    (facs_p : List (Polynomial (ZMod p)))
    (hne : facs_p ≠ [])
    (hprod : Polynomial.map (Int.castRingHom (ZMod p)) f = facs_p.prod)
    (hcop : facs_p.Pairwise (fun a b => IsCoprime a b))
    : HenselCorrect f k facs_p (hensel_l1 hp_size k hk f facs_p) :=
  hensel_lift_correct p k hk f facs_p hne hprod hcop

-- ============================================================
-- §6. 因子重组（L1 包装 — TODO: 替换为 __factor_recombine_upoly_ir 精化版本）
-- ============================================================

/-- L1 因子重组 — 目前使用 `recombine_correct`（UFD 存在性，忽略 facs_pk）。
    正确性由 `__factor_recombine_upoly_ir_refines` 保证。
    当前暂用 L2 `recombine_correct`（待精化定理填补后改回 C++ 路径）。 -/
noncomputable def recombine_l1 (hp_size : 2 * p ≤ UInt64.size) (k : ℕ) (f : Polynomial ℤ) (hf : f ≠ 0)
    (facs_pk : List (Polynomial (ZMod (p ^ k)))) : List (Polynomial ℤ) :=
  (recombine_correct f hf).choose

lemma recombine_l1_correct (hp_size : 2 * p ≤ UInt64.size) (k : ℕ) (f : Polynomial ℤ) (hf : f ≠ 0)
    (facs_pk : List (Polynomial (ZMod (p ^ k))))
    (hprod : Polynomial.map (Int.castRingHom (ZMod (p ^ k))) f = facs_pk.prod)
    : RecombineCorrect f (recombine_l1 hp_size k f hf facs_pk) :=
  (recombine_correct f hf).choose_spec

-- ============================================================
-- §7. C++ 端到端因式分解定理
-- ============================================================

/-- C++ 翻译代码的 Z[x] 因式分解正确性定理。

    使用 L1 包装（C++ 翻译函数）实例化 `factor_ZZ_correct` 的三个子过程。
    当前 L1 包装中 SQF 已使用 C++ 翻译路径；DDF/EDF 仍等待精化定理。
    Hensel 和 Recombine 同理。

    精化定理状态：
     - ✅ `__squarefree_Zp_ir_refines`
     - ❌ `__ddf_Zp_ir_refines`
     - ✅ `__symmetric_mod_ir_refines`
     - ❌ `__binomial_ir_refines`
     - ❌ `__isqrt_ceil_ir_refines`
     - ❌ `__edf_Zp_ir_refines`（待新增）
     - ❌ `__hensel_lift_*_ir_refines`（待新增）
     - ❌ `__recombine_*_ir_refines`（待新增）
 -/
theorem factor_ZZ_cpp_correct
    (f : Polynomial ℤ) (hf : f ≠ 0) (hprim : f.IsPrimitive)
    {p : ℕ} [hp : Fact (Nat.Prime p)] {k : ℕ} (hk : 0 < k)
    (hp_size : 2 * p ≤ UInt64.size) (hp2 : p * p ≤ UInt64.size)
    (hgood : Squarefree (Polynomial.map (Int.castRingHom (ZMod p)) f))
    (hdeg : (Polynomial.map (Int.castRingHom (ZMod p)) f).natDegree = f.natDegree)
    (hfbound : f.natDegree * p < 2 ^ 64)
    : ∃ result : List (Polynomial ℤ), FactorZZCorrect f result := by
  have hfp_ne : (Polynomial.map (Int.castRingHom (ZMod p)) f) ≠ 0 :=
    fun h => not_squarefree_zero (h ▸ hgood)
  have hb : (Polynomial.map (Int.castRingHom (ZMod p)) f).natDegree * p < 2 ^ 64 := by
    rw [hdeg]; exact hfbound
  exact factor_ZZ_correct f hf hprim hk hgood hdeg
    (factor_zp_l1_func hp_size hp2)
    (factor_zp_l1_func_correct hp_size hp2 _ hfp_ne hb)
    (hensel_l1 hp_size k hk f) (hensel_l1_correct hp_size k hk f)
    (recombine_l1 hp_size k f hf) (recombine_l1_correct hp_size k f hf)

end L1Pipeline
