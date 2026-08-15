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
import CLPoly.Refinement.Generated
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

-- ============================================================
-- §2. Strict generated stage contracts consumed by the pipeline
-- ============================================================

/-- Pipeline boundary for the exact generated C++ `__squarefree_Zp` entry.
The proof source remains the centralized generated refinement contract. -/
theorem strictSQFStage
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (physical : Refinement.StrictSquarefreeZp.YunRawGCDWorkspaceProvider
      this hcfg)
    (source : SparsePolyZp)
    (hcanonical : SparsePolyZp.Canonical this._p.toNat source)
    (hmonic : (SparsePolyZp.toPoly this._p.toNat source).Monic)
    (hnonempty : 0 < source.size)
    (hpositive : 0 < (SparsePolyZp.toPoly this._p.toNat source).natDegree)
    (hbound : CLPoly.Impl.StrictPolynomialGCDRefinement.sparseDenseLength
      source ≤ 2 ^ 63) :
    ∃ factors,
      Generated.StrictSquarefreeZp.__squarefree_Zp_raw_ir
          (Refinement.StrictSquarefreeGenerated.strictSQFRawOps
            this hcfg physical)
          source (fun _ => ⟨hcanonical, hmonic, hnonempty, hpositive,
            hbound⟩) = .ok factors ∧
      toPolyList factors this._p.toNat =
        sqfZp (SparsePolyZp.toPoly this._p.toNat source) ∧
      ∀ item ∈ factors.toList,
        SparsePolyZp.Canonical this._p.toNat item.1 := by
  exact Refinement.__squarefree_Zp_raw_ir_refines_sqfZp this hcfg physical
    source hcanonical hmonic hnonempty hpositive hbound

/-- Exact representation and L2 facts required to pass one concrete SQF
output factor to the generated DDF entry. -/
structure StrictSQFFactorDDFReady (this : DenseUPolyZp)
    (factor : SparsePolyZp) : Prop where
  primeMatches : factor[0]!.2.prime = this._p
  canonical : SparsePolyZp.Canonical this._p.toNat factor
  degreeBound : ∀ term ∈ factor.toList, term.1.deg < 2 ^ 62
  monic : (SparsePolyZp.toPoly this._p.toNat factor).Monic
  squarefree : Squarefree (SparsePolyZp.toPoly this._p.toNat factor)

private theorem sqfZp_factor_natDegree_le
    {q : Nat} [Fact (Nat.Prime q)] (f : Polynomial (ZMod q))
    (hf : f ≠ 0) (item : Polynomial (ZMod q) × Nat)
    (hitem : item ∈ sqfZp f) : item.1.natDegree ≤ f.natDegree := by
  have hcorrect := sqf_correct f hf
  have hexponent : 1 ≤ item.2 := hcorrect.2.2.1 item hitem
  have hpowMem : item.1 ^ item.2 ∈
      (sqfZp f).map (fun factor => factor.1 ^ factor.2) :=
    List.mem_map.mpr ⟨item, hitem, rfl⟩
  have hfactorDvdPow : item.1 ∣ item.1 ^ item.2 :=
    dvd_pow_self item.1 (Nat.one_le_iff_ne_zero.mp hexponent)
  have hfactorDvd : item.1 ∣ f := hfactorDvdPow.trans
    ((List.dvd_prod hpowMem).trans hcorrect.1.symm.dvd)
  exact Polynomial.natDegree_le_of_dvd hfactorDvd hf

/-- The concrete result array produced by strict SQF carries every condition
needed by strict DDF.  No property is reconstructed by changing the result:
the proof follows membership in the actual generated array. -/
theorem strictSQFStage_preparesDDF
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (physical : Refinement.StrictSquarefreeZp.YunRawGCDWorkspaceProvider
      this hcfg)
    (source : SparsePolyZp)
    (hcanonical : SparsePolyZp.Canonical this._p.toNat source)
    (hmonic : (SparsePolyZp.toPoly this._p.toNat source).Monic)
    (hnonempty : 0 < source.size)
    (hpositive : 0 < (SparsePolyZp.toPoly this._p.toNat source).natDegree)
    (hbound : CLPoly.Impl.StrictPolynomialGCDRefinement.sparseDenseLength
      source ≤ 2 ^ 63)
    (hsourceDegree :
      (SparsePolyZp.toPoly this._p.toNat source).natDegree < 2 ^ 62) :
    ∃ factors,
      Generated.StrictSquarefreeZp.__squarefree_Zp_raw_ir
          (Refinement.StrictSquarefreeGenerated.strictSQFRawOps
            this hcfg physical)
          source (fun _ => ⟨hcanonical, hmonic, hnonempty, hpositive,
            hbound⟩) = .ok factors ∧
      toPolyList factors this._p.toNat =
        sqfZp (SparsePolyZp.toPoly this._p.toNat source) ∧
      ∀ item ∈ factors.toList,
        StrictSQFFactorDDFReady this item.1 := by
  rcases strictSQFStage this hcfg physical source hcanonical hmonic
      hnonempty hpositive hbound with
    ⟨factors, hrun, hsemantic, hresultCanonical⟩
  refine ⟨factors, hrun, hsemantic, ?_⟩
  intro item hitem
  have hdecoded :
      (SparsePolyZp.toPoly this._p.toNat item.1, item.2.toNat) ∈
        sqfZp (SparsePolyZp.toPoly this._p.toNat source) := by
    rw [← hsemantic]
    unfold toPolyList
    simp only [Array.toList_map]
    exact List.mem_map.mpr ⟨item, hitem, rfl⟩
  have hproperties :=
    (sqf_correct (SparsePolyZp.toPoly this._p.toNat source)
      hmonic.ne_zero).2.1 _ hdecoded
  have hitemCanonical := hresultCanonical item hitem
  have hitemNonempty :=
    Refinement.StrictSquarefreeZp.sparsePolyZp_size_pos_of_toPoly_ne_zero
      this._p.toNat item.1 hproperties.2.ne_zero
  have hheadMem : item.1[0] ∈ item.1.toList := by
    simpa using Array.getElem_mem item.1 0 hitemNonempty
  have hhead : item.1[0]! = item.1[0] := by
    rw [getElem!_def, getElem?_def]
    simp [hitemNonempty]
  have hprimeNat : item.1[0]!.2.prime.toNat = this._p.toNat := by
    rw [hhead]
    exact (hitemCanonical.1 item.1[0] hheadMem).1
  refine ⟨UInt64.toNat_inj.mp hprimeNat, hitemCanonical, ?_,
    hproperties.2, hproperties.1⟩
  intro term hterm
  have htermLe := Refinement.StrictDDF.canonical_term_degree_le_natDegree
    this._p.toNat item.1 hitemCanonical term hterm
  have hfactorLe := sqfZp_factor_natDegree_le
    (SparsePolyZp.toPoly this._p.toNat source) hmonic.ne_zero
    (SparsePolyZp.toPoly this._p.toNat item.1, item.2.toNat) hdecoded
  exact lt_of_le_of_lt (htermLe.trans hfactorLe) hsourceDegree

/-- Pipeline boundary for the exact generated C++ `__ddf_Zp` entry. -/
theorem strictDDFStage
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (providers : Refinement.StrictDDF.DDFRawProviders this)
    (f : SparsePolyZp) (hfPrime : f[0]!.2.prime = this._p)
    (hfCanonical : SparsePolyZp.Canonical this._p.toNat f)
    (hfDegree : ∀ term ∈ f.toList, term.1.deg < 2 ^ 62)
    (hfMonic : (SparsePolyZp.toPoly this._p.toNat f).Monic)
    (hfSquarefree : Squarefree (SparsePolyZp.toPoly this._p.toNat f)) :
    ∃ output,
      Generated.StrictDDF.__ddf_Zp_raw_ir
          (Refinement.StrictDDF.strictDDFRawOps this providers
            (SparsePolyZp.toPoly this._p.toNat f)) f
          (fun _ => Refinement.StrictDDF.DDFLoopInvariant.initial this f
            f[0]!.2.prime hfPrime hfCanonical hfDegree hfMonic
            hfSquarefree) = .ok output ∧
      Refinement.StrictDDF.ddfResultToL2 this._p.toNat output =
        ddf (SparsePolyZp.toPoly this._p.toNat f) ∧
      (∀ item ∈ output.toList,
        SparsePolyZp.Canonical this._p.toNat item.1) ∧
      (∀ item ∈ output.toList,
        0 < (SparsePolyZp.toPoly this._p.toNat item.1).natDegree) ∧
      ∀ item ∈ output.toList, 0 < item.2.toNat := by
  exact Refinement.__ddf_Zp_raw_ir_refines_ddf this providers f hfPrime
    hfCanonical hfDegree hfMonic hfSquarefree

/-- Genuine SQF-to-DDF execution composition.  Every concrete sparse factor
returned by the generated SQF run is passed unchanged to a generated DDF run;
the existential output is the result of that raw execution. -/
theorem strictSQFStage_runsDDF
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (physical : Refinement.StrictSquarefreeZp.YunRawGCDWorkspaceProvider
      this hcfg)
    (providers : Refinement.StrictDDF.DDFRawProviders this)
    (source : SparsePolyZp)
    (hcanonical : SparsePolyZp.Canonical this._p.toNat source)
    (hmonic : (SparsePolyZp.toPoly this._p.toNat source).Monic)
    (hnonempty : 0 < source.size)
    (hpositive : 0 < (SparsePolyZp.toPoly this._p.toNat source).natDegree)
    (hbound : CLPoly.Impl.StrictPolynomialGCDRefinement.sparseDenseLength
      source ≤ 2 ^ 63)
    (hsourceDegree :
      (SparsePolyZp.toPoly this._p.toNat source).natDegree < 2 ^ 62) :
    ∃ factors,
      Generated.StrictSquarefreeZp.__squarefree_Zp_raw_ir
          (Refinement.StrictSquarefreeGenerated.strictSQFRawOps
            this hcfg physical)
          source (fun _ => ⟨hcanonical, hmonic, hnonempty, hpositive,
            hbound⟩) = .ok factors ∧
      toPolyList factors this._p.toNat =
        sqfZp (SparsePolyZp.toPoly this._p.toNat source) ∧
      ∀ item ∈ factors.toList,
        ∃ (hfactorReady : StrictSQFFactorDDFReady this item.1)
          (output : Array (SparsePolyZp × UInt64)),
          Generated.StrictDDF.__ddf_Zp_raw_ir
              (Refinement.StrictDDF.strictDDFRawOps this providers
                (SparsePolyZp.toPoly this._p.toNat item.1)) item.1
              (fun _ => Refinement.StrictDDF.DDFLoopInvariant.initial
                this item.1 item.1[0]!.2.prime
                hfactorReady.primeMatches hfactorReady.canonical
                hfactorReady.degreeBound hfactorReady.monic
                hfactorReady.squarefree) = .ok output ∧
          Refinement.StrictDDF.ddfResultToL2 this._p.toNat output =
            ddf (SparsePolyZp.toPoly this._p.toNat item.1) ∧
          (∀ component ∈ output.toList,
            SparsePolyZp.Canonical this._p.toNat component.1) ∧
          (∀ component ∈ output.toList,
            0 < (SparsePolyZp.toPoly this._p.toNat component.1).natDegree) ∧
          ∀ component ∈ output.toList, 0 < component.2.toNat := by
  rcases strictSQFStage_preparesDDF this hcfg physical source hcanonical
      hmonic hnonempty hpositive hbound hsourceDegree with
    ⟨factors, hsqfRun, hsqfSemantic, hready⟩
  refine ⟨factors, hsqfRun, hsqfSemantic, ?_⟩
  intro item hitem
  have hfactorReady := hready item hitem
  rcases strictDDFStage this providers item.1 hfactorReady.primeMatches
      hfactorReady.canonical hfactorReady.degreeBound hfactorReady.monic
      hfactorReady.squarefree with
    ⟨output, hrun, hsemantic, hcanonicalOutput, hdegreePositive,
      hindexPositive⟩
  exact ⟨hfactorReady, output, hrun, hsemantic, hcanonicalOutput,
    hdegreePositive, hindexPositive⟩

/-- A generated DDF execution prepares every concrete output component for
the generated EDF entry.  The equal-degree property comes from DDF semantics;
all representation facts remain attached to the actual raw output array. -/
theorem strictDDFStage_preparesEDF
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (providers : Refinement.StrictDDF.DDFRawProviders this)
    (f : SparsePolyZp) (hready : StrictSQFFactorDDFReady this f) :
    ∃ output,
      Generated.StrictDDF.__ddf_Zp_raw_ir
          (Refinement.StrictDDF.strictDDFRawOps this providers
            (SparsePolyZp.toPoly this._p.toNat f)) f
          (fun _ => Refinement.StrictDDF.DDFLoopInvariant.initial this f
            f[0]!.2.prime hready.primeMatches hready.canonical
            hready.degreeBound hready.monic hready.squarefree) = .ok output ∧
      Refinement.StrictDDF.ddfResultToL2 this._p.toNat output =
        ddf (SparsePolyZp.toPoly this._p.toNat f) ∧
      ∀ item ∈ output.toList,
        Refinement.StrictEDF.EDFEntryInvariant this item.1 item.2 := by
  rcases strictDDFStage this providers f hready.primeMatches
      hready.canonical hready.degreeBound hready.monic hready.squarefree with
    ⟨output, hrun, hsemantic, hcanonical, hdegreePositive,
      hindexPositive⟩
  refine ⟨output, hrun, hsemantic, ?_⟩
  intro item hitem
  have hdecoded :
      (SparsePolyZp.toPoly this._p.toNat item.1, item.2.toNat) ∈
        ddf (SparsePolyZp.toPoly this._p.toNat f) := by
    rw [← hsemantic]
    unfold Refinement.StrictDDF.ddfResultToL2
    exact List.mem_map.mpr ⟨item, hitem, rfl⟩
  have hcorrect := ddf_correct (SparsePolyZp.toPoly this._p.toNat f)
    hready.monic hready.squarefree
  have hcomponentDvd := hcorrect.1 _ hdecoded
  have hcomponentMonic := hcorrect.2.2.2.2 _ hdecoded
  have hcomponentCanonical := hcanonical item hitem
  have hcomponentPositive := hdegreePositive item hitem
  have hcomponentNonempty :=
    Refinement.StrictSquarefreeZp.sparsePolyZp_size_pos_of_toPoly_ne_zero
      this._p.toNat item.1 hcomponentMonic.ne_zero
  have hheadMem : item.1[0] ∈ item.1.toList := by
    simpa using Array.getElem_mem item.1 0 hcomponentNonempty
  have hhead : item.1[0]! = item.1[0] := by
    rw [getElem!_def, getElem?_def]
    simp [hcomponentNonempty]
  have hprimeNat : item.1[0]!.2.prime.toNat = this._p.toNat := by
    rw [hhead]
    exact (hcomponentCanonical.1 item.1[0] hheadMem).1
  have hfDegree := Refinement.StrictDDF.canonical_natDegree_lt_of_terms_lt
    this._p.toNat f hready.canonical hready.monic.ne_zero (2 ^ 62)
    hready.degreeBound
  have hcomponentDegree :
      (SparsePolyZp.toPoly this._p.toNat item.1).natDegree < 2 ^ 62 :=
    lt_of_le_of_lt
      (Polynomial.natDegree_le_of_dvd hcomponentDvd hready.monic.ne_zero)
      hfDegree
  refine ⟨hcomponentCanonical, ?_, ?_, hcomponentMonic,
    hcomponentPositive, hindexPositive item hitem,
    hready.squarefree.squarefree_of_dvd hcomponentDvd, ?_⟩
  · intro _
    exact UInt64.toNat_inj.mp hprimeNat
  · intro term hterm
    exact lt_of_le_of_lt
      (Refinement.StrictDDF.canonical_term_degree_le_natDegree
        this._p.toNat item.1 hcomponentCanonical term hterm)
      hcomponentDegree
  · intro q hq hqDvd
    exact hcorrect.2.1 _ hdecoded q hq hqDvd

/-- Pipeline boundary for the exact generated C++ `__edf_Zp` entry,
including the concrete RNG transition. -/
theorem strictEDFStage
    {State : Type} (engine : Generated.StrictEDF.RandomEngine State)
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (providers : Refinement.StrictDDF.DDFRawProviders this)
    (termination : Generated.StrictEDF.EDFTermination
      (Refinement.StrictEDF.strictEDFRawOps engine this providers))
    (result : Array SparsePolyZp) (f : SparsePolyZp) (d : UInt64)
    (rng : State)
    (hinvariant : Refinement.StrictEDF.EDFEntryInvariant this f d) :
    ∃ output rng' factors,
      Generated.StrictEDF.__edf_Zp_raw_ir
          (Refinement.StrictEDF.strictEDFRawOps engine this providers)
          (Refinement.StrictEDF.strictEDFSplitLaw engine this providers)
          termination result f d rng hinvariant = .ok (output, rng') ∧
      Refinement.StrictEDF.edfResultToL2 this._p.toNat output =
        Refinement.StrictEDF.edfResultToL2 this._p.toNat result ++ factors ∧
      EDFCorrect (SparsePolyZp.toPoly this._p.toNat f) d.toNat factors := by
  exact Refinement.__edf_Zp_raw_ir_refines_edf engine this providers
    termination result f d rng hinvariant

/-- Genuine DDF-to-EDF execution composition.  Each concrete DDF component
is fed unchanged to the generated EDF entry, for an arbitrary incoming RNG
state; the outgoing RNG state is the one returned by that raw execution. -/
theorem strictDDFStage_runsEDF
    {State : Type} (engine : Generated.StrictEDF.RandomEngine State)
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (providers : Refinement.StrictDDF.DDFRawProviders this)
    (termination : Generated.StrictEDF.EDFTermination
      (Refinement.StrictEDF.strictEDFRawOps engine this providers))
    (f : SparsePolyZp) (hready : StrictSQFFactorDDFReady this f) :
    ∃ ddfOutput,
      Generated.StrictDDF.__ddf_Zp_raw_ir
          (Refinement.StrictDDF.strictDDFRawOps this providers
            (SparsePolyZp.toPoly this._p.toNat f)) f
          (fun _ => Refinement.StrictDDF.DDFLoopInvariant.initial this f
            f[0]!.2.prime hready.primeMatches hready.canonical
            hready.degreeBound hready.monic hready.squarefree) =
        .ok ddfOutput ∧
      Refinement.StrictDDF.ddfResultToL2 this._p.toNat ddfOutput =
        ddf (SparsePolyZp.toPoly this._p.toNat f) ∧
      ∀ component ∈ ddfOutput.toList, ∀ rng : State,
        ∃ (hinvariant : Refinement.StrictEDF.EDFEntryInvariant this
              component.1 component.2)
          (edfOutput : Array SparsePolyZp) (rng' : State)
          (factors : List (Polynomial (ZMod this._p.toNat))),
          Generated.StrictEDF.__edf_Zp_raw_ir
              (Refinement.StrictEDF.strictEDFRawOps engine this providers)
              (Refinement.StrictEDF.strictEDFSplitLaw engine this providers)
              termination #[] component.1 component.2 rng hinvariant =
            .ok (edfOutput, rng') ∧
          Refinement.StrictEDF.edfResultToL2 this._p.toNat edfOutput =
            factors ∧
          EDFCorrect (SparsePolyZp.toPoly this._p.toNat component.1)
            component.2.toNat factors := by
  rcases strictDDFStage_preparesEDF this providers f hready with
    ⟨ddfOutput, hddfRun, hddfSemantic, hreadyEDF⟩
  refine ⟨ddfOutput, hddfRun, hddfSemantic, ?_⟩
  intro component hcomponent rng
  have hinvariant := hreadyEDF component hcomponent
  rcases strictEDFStage engine this providers termination #[] component.1
      component.2 rng hinvariant with
    ⟨edfOutput, rng', factors, hedfRun, hedfSemantic, hcorrect⟩
  refine ⟨hinvariant, edfOutput, rng', factors, hedfRun, ?_, hcorrect⟩
  simpa using hedfSemantic

/-- Pipeline boundary for the exact generated C++ quadratic
`__hensel_step` entry. -/
theorem strictHenselStepStage
    (node : HenselNode) (f : SparsePolyZZ) (m : Nat)
    (hinvariant : Refinement.StrictHensel.HenselStepRefinementInvariant
      Refinement.StrictHensel.concreteDivmodTermination node f m) :
    ∃ output,
      Generated.StrictHensel.__hensel_step_raw_ir
          (Refinement.StrictHensel.strictHenselRawOps
            Refinement.StrictHensel.concreteDivmodTermination)
          node f (m : Int) = .ok output ∧
      Refinement.StrictHensel.HenselStepCorrect f m node output := by
  exact Refinement.__hensel_step_raw_ir_refines node f m hinvariant

/-- Pipeline boundary for the exact generated C++ `__hensel_tree_build`
entry.  Its concrete output supplies both the initial algebraic root invariant
and the whole-tree topology required by the later raw lift/extraction stages. -/
theorem strictHenselTreeBuildStage
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (h2p : 2 * this._p.toNat ≤ UInt64.size)
    (hp2 : this._p.toNat * this._p.toNat ≤ UInt64.size)
    (mulProvider : Refinement.StrictDDF.RawMulWorkspaceProvider this)
    (factors : Array SparsePolyZp)
    (hfactors : ∀ factor ∈ factors.toList,
      SparsePolyZp.Canonical this._p.toNat factor)
    (hfactorsNonempty : ∀ factor ∈ factors.toList, 0 < factor.size)
    (hpairwise : ∀ i j (hi : i < factors.size) (hj : j < factors.size),
      i < j → IsCoprime
        (SparsePolyZp.toPoly this._p.toNat (getElem factors i hi))
        (SparsePolyZp.toPoly this._p.toNat (getElem factors j hj)))
    (htwo : 2 ≤ factors.size)
    (hfitsInt32 : Refinement.StrictHensel.henselTreeInternalNodeCount
      0 factors.size < 2 ^ 31) :
    let tree := Refinement.StrictHensel.henselTreeBuildTopology
      0 factors.size 0
    ∃ output,
      Generated.StrictHensel.__hensel_tree_build_raw_ir
          (Refinement.StrictHensel.strictHenselTreeBuildRawOps
            this mulProvider) factors this._p = .ok output ∧
      ∃ hroot : 0 < output.size,
      output.size = tree.nodeCount ∧ tree.rootIndex = 0 ∧
      Refinement.StrictHensel.liftChildMatches
        (getElem output 0 hroot).left
        (match tree with | .node _ left _ => left) ∧
      Refinement.StrictHensel.liftChildMatches
        (getElem output 0 hroot).right
        (match tree with | .node _ _ right => right) ∧
      Refinement.StrictHensel.HenselExtractInvariant tree output ∧
      Refinement.StrictHensel.HenselTreeSemanticBuildCertificate
        this._p.toNat factors 0 0 factors.size tree output ∧
      Refinement.StrictHensel.HenselTreeNodeInitialInvariant
        this._p.toNat factors 0 factors.size (getElem output 0 hroot) ∧
      Refinement.StrictHensel.HenselArrayCanonical output := by
  exact Refinement.__hensel_tree_build_raw_ir_refines this hcfg h2p hp2
    mulProvider factors hfactors hfactorsNonempty hpairwise htwo hfitsInt32

/-- Pipeline boundary for the complete generated C++
`__hensel_lift_upoly` entry.  It exposes the actual strict L1 execution and
the full target/adjust/build/lift/extract/normalize L2 trace. -/
theorem strictHenselLiftUpolyStage
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (h2p : 2 * this._p.toNat ≤ UInt64.size)
    (hp2 : this._p.toNat * this._p.toNat ≤ UInt64.size)
    (mulProvider : Refinement.StrictDDF.RawMulWorkspaceProvider this)
    (f : SparsePolyZZ) (factors : Array SparsePolyZp) (aTarget : Int32)
    (hinvariant : Refinement.StrictHensel.HenselLiftEntryInvariant this
      Refinement.StrictHensel.concreteDivmodTermination mulProvider f factors
      aTarget) :
    ∃ output,
      Generated.StrictHensel.__hensel_lift_upoly_raw_ir
          (Refinement.StrictHensel.strictHenselRawOps
            Refinement.StrictHensel.concreteDivmodTermination)
          (Refinement.StrictHensel.strictHenselTreeBuildRawOps
            this mulProvider)
          f factors this._p aTarget
            (Nat.Prime.two_le (Fact.out : Nat.Prime this._p.toNat)) =
              .ok output ∧
      Refinement.StrictHensel.HenselLiftEntryCorrect
        Refinement.StrictHensel.concreteDivmodTermination f factors this._p
        aTarget output := by
  exact Refinement.__hensel_lift_upoly_raw_ir_refines this hcfg h2p hp2
    mulProvider f factors aTarget hinvariant

/-- End-to-end pipeline boundary for the original generated C++ `__factor_Zp`
entry.  This is a direct re-export of the centralized generated contract, not
an L2-only reconstruction from `FactorZpInstantiate`. -/
theorem strictFactorZpStage
    {State : Type} (engine : Generated.StrictEDF.RandomEngine State)
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (sqfPhysical : Refinement.StrictSquarefreeZp.YunRawGCDWorkspaceProvider
      this hcfg)
    (providers : Refinement.StrictDDF.DDFRawProviders this)
    (termination : Generated.StrictEDF.EDFTermination
      (Refinement.StrictEDF.strictEDFRawOps engine this providers))
    (sort : Refinement.StrictFactorZp.SortByDegreeProvider)
    (initialRng : State) (f : SparsePolyZp)
    (hcanonical : SparsePolyZp.Canonical this._p.toNat f)
    (hnonempty : 0 < f.size)
    (hdegreePositive :
      0 < (SparsePolyZp.toPoly this._p.toNat f).natDegree)
    (hdegreeBound :
      (SparsePolyZp.toPoly this._p.toNat f).natDegree < 2 ^ 62) :
    let ops : Generated.StrictFactorZp.FactorZpRawOps State := {
      makeMonic := Refinement.StrictSquarefreeZp.upolyMakeMonicIR this
      squarefree := Refinement.StrictFactorZp.strictSQFCall
        this hcfg sqfPhysical
      ddf := Refinement.StrictFactorZp.strictDDFCall this providers
      edf := Refinement.StrictFactorZp.strictEDFCall
        engine this providers termination
      sortByDegree := sort.run }
    ∃ lc output,
      Generated.StrictFactorZp.__factor_Zp_raw_ir ops initialRng f =
        .ok (lc, output) ∧
      FactorZpCorrect (SparsePolyZp.toPoly this._p.toNat f)
        (Zp.toZMod this._p.toNat lc)
        (Refinement.StrictFactorZp.factorResultToL2
          this._p.toNat output) := by
  exact Refinement.__factor_Zp_raw_ir_refines_FactorZpCorrect engine this
    hcfg sqfPhysical providers termination sort initialRng f hcanonical
    hnonempty hdegreePositive hdegreeBound

/- The ZZ end-to-end wrapper must compose the concrete prime-selection,
Zp-factorization, Hensel and recombination executions. -/

end L1Pipeline
