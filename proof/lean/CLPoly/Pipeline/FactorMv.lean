/-
  CLPoly/Pipeline/FactorMv.lean — 多变量 Z[x₁,...,xₙ] 因式分解

  1. mv_factor_instantiate：直接 UFD
  2. mv_factor_correct：参数化框架（SQF + Wang → 构造 f 的不可约分解）
-/
import CLPoly.Spec
import Mathlib.RingTheory.Polynomial.UniqueFactorization

set_option autoImplicit false
set_option maxHeartbeats 1600000

open MvPolynomial

-- ============================================================
-- 1. 直接版（UFD）
-- ============================================================

theorem mv_factor_instantiate {σ : Type*}
    (f : MvPolynomial σ ℤ) (hf : f ≠ 0)
    : ∃ result : List (MvPolynomial σ ℤ), MvFactorCorrect f result := by
  obtain ⟨factors, hirred, hassoc⟩ := WfDvdMonoid.exists_factors f hf
  refine ⟨factors.toList, ?_, fun g hg => hirred g (Multiset.mem_toList.mp hg)⟩
  rw [Multiset.prod_toList]; exact hassoc.symm

-- ============================================================
-- 2. 辅助引理
-- ============================================================

private lemma associated_prod_of_forall₂ {α : Type*} [CommMonoid α]
    {l₁ l₂ : List α} (h : List.Forall₂ Associated l₁ l₂) :
    Associated l₁.prod l₂.prod := by
  induction h with
  | nil => exact Associated.refl 1
  | cons h _ ih => simp only [List.prod_cons]; exact h.mul_mul ih

private lemma prod_flatMap_replicate {α : Type*} [CommMonoid α]
    (l : List (α × ℕ)) (g : α → List α) :
    (l.flatMap (fun pr => (List.replicate pr.2 (g pr.1)).flatten)).prod =
    (l.map (fun pr => (g pr.1).prod ^ pr.2)).prod := by
  induction l with
  | nil => simp
  | cons pr rest ih =>
    simp only [List.flatMap_cons, List.map_cons, List.prod_cons, List.prod_append]
    congr 1
    · induction pr.2 with
      | zero => simp
      | succ n ihn =>
        simp only [List.replicate_succ, List.flatten_cons, List.prod_append, ihn, pow_succ,
                   mul_comm]

-- ============================================================
-- 3. 参数化框架
-- ============================================================

/-- 多变量因式分解参数化框架。

    从 SQF 和 Wang/单变量的输出**构造** f 的不可约分解：
    - SQF: f ~ ∏ gₖ^mₖ（每个 gₖ squarefree）
    - Wang: 每个 gₖ ~ ∏ hₖⱼ（每个 hₖⱼ irreducible）
    - result = flatten(replicate mₖ [hₖ₁,...] for each k)
    - f ~ result.prod 且每个元素 irreducible -/
theorem mv_factor_correct {n : ℕ}
    (f : MvPolynomial (Fin n) ℤ) (_hf : f ≠ 0)
    -- 子过程 1: SQF 分解
    (sqf : List (MvPolynomial (Fin n) ℤ × ℕ))
    (hsqf_prod : Associated f (sqf.map (fun pr => pr.1 ^ pr.2)).prod)
    (_hsqf_props : ∀ pr ∈ sqf, Squarefree pr.1 ∧ pr.1 ≠ 0 ∧ pr.2 ≥ 1)
    -- 子过程 2: 每个 squarefree 分量的不可约分解
    (factor_sqfree : MvPolynomial (Fin n) ℤ → List (MvPolynomial (Fin n) ℤ))
    (hfactor_prod : ∀ pr ∈ sqf, Associated pr.1 (factor_sqfree pr.1).prod)
    (hfactor_irred : ∀ pr ∈ sqf, ∀ g ∈ factor_sqfree pr.1, Irreducible g)
    : ∃ result : List (MvPolynomial (Fin n) ℤ), MvFactorCorrect f result := by
  -- 构造 result
  set result := sqf.flatMap (fun pr =>
    (List.replicate pr.2 (factor_sqfree pr.1)).flatten) with hresult_def
  refine ⟨result, ?_, ?_⟩
  · -- Associated f result.prod
    rw [hresult_def, prod_flatMap_replicate]
    -- Goal: Associated f (sqf.map (fun pr => (factor_sqfree pr.1).prod ^ pr.2)).prod
    -- 由 hsqf_prod 和每个分量的 Associated 组合
    exact hsqf_prod.trans (associated_prod_of_forall₂ (by
      rw [List.forall₂_map_right_iff, List.forall₂_map_left_iff]
      exact List.forall₂_same.mpr (fun pr hpr =>
        (hfactor_prod pr hpr).pow_pow)))
  · -- ∀ g ∈ result, Irreducible g
    intro g hg
    simp only [hresult_def, List.mem_flatMap, List.mem_flatten, List.mem_replicate] at hg
    obtain ⟨pr, hpr_mem, l, ⟨_, hl_eq⟩, hg_mem⟩ := hg
    rw [hl_eq] at hg_mem
    exact hfactor_irred pr hpr_mem g hg_mem
