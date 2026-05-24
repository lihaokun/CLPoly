/-
  CLPoly/Pipeline/L1.lean — L1 翻译代码的管线实例化

  用 L1 `_ir` 函数（来自 Corpus.lean，翻译自 C++）实例化 factor_Zp_correct。
  精化定理填补后，此文件直接验证翻译出的 C++ 代码。

  当前 L1 覆盖：sqf (__squarefree_Zp_ir), ddf (__ddf_Zp_ir)
  EDF 暂用 L2 edf_correct_unconditional（同 FactorZpInstantiate）
-/

import CLPoly.Pipeline.FactorZp
import CLPoly.Algorithm.EDF
import CLPoly.Generated.Corpus
import CLPoly.Refinement.Basic

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
  (f.support.sort (· ≤ ·)).map (fun n =>
    ({deg := n}, Zp.ofUInt64 (UInt64.ofNat (ZMod.val (f.coeff n))) (UInt64.ofNat p))
  ) |>.toArray

lemma toSparsePolyZp_wellFormed (f : (ZMod p)[X]) :
    SparsePolyZp.WellFormed p (toSparsePolyZp f) := by
  sorry

lemma toSparsePolyZp_allReduced (f : (ZMod p)[X]) :
    SparsePolyZp.AllReduced p (toSparsePolyZp f).toList := by
  sorry

lemma toSparsePolyZp_toPoly (f : (ZMod p)[X]) :
    SparsePolyZp.toPoly p (toSparsePolyZp f) = f := by
  sorry

-- ============================================================
-- §2. L1 Wrappers
-- ============================================================

/-- L1 sqf wrapper — 使用 __squarefree_Zp_ir（Corpus.lean，翻译自 C++） -/
noncomputable def sqfZp_l1 (f : (ZMod p)[X]) : List ((ZMod p)[X] × ℕ) :=
  toPolyList (Generated.__squarefree_Zp_ir (toSparsePolyZp f)) p

theorem sqfZp_l1_correct (f : (ZMod p)[X]) (hf : f ≠ 0) :
    SquarefreeDecomp f (sqfZp_l1 f) := by
  sorry

/-- L1 ddf wrapper — 使用 __ddf_Zp_ir（Corpus.lean，翻译自 C++） -/
noncomputable def ddf_l1 (f : (ZMod p)[X]) : List ((ZMod p)[X] × ℕ) :=
  toPolyList (Generated.__ddf_Zp_ir (toSparsePolyZp f)) p

theorem ddf_l1_correct (f : (ZMod p)[X]) (hm : Monic f) (hsq : Squarefree f) :
    DDFCorrect f (ddf_l1 f) := by
  sorry

-- ============================================================
-- §3. EDF（暂用 L2 edf_correct_unconditional，同 FactorZpInstantiate）
-- ============================================================

/-- L1 edf wrapper — 目前使用 L2 edf_correct_unconditional（Classical.choose）。
    TODO: 替换为 __edf_Zp_ir 精化版本。 -/
noncomputable def edf_l1 (g : (ZMod p)[X]) (d : ℕ) : List ((ZMod p)[X]) :=
  if hg_deg : g.natDegree = 0 then []
  else
    if hpre : Monic g ∧ Squarefree g ∧ 0 < d ∧
        (∀ q : (ZMod p)[X], Irreducible q → q ∣ g → q.natDegree = d) then
      (edf_correct_unconditional g d hpre.1 hpre.2.1
        (Nat.pos_of_ne_zero hg_deg) hpre.2.2.1 hpre.2.2.2).choose
    else [g]

theorem edf_l1_correct (g : (ZMod p)[X]) (d : ℕ)
    (hm : Monic g) (hsq : Squarefree g)
    (hdeg : ∀ q, Irreducible q → q ∣ g → q.natDegree = d) :
    EDFCorrect g d (edf_l1 g d) := by
  by_cases hg_deg : g.natDegree = 0
  · have htarget : edf_l1 g d = [] := by
      unfold edf_l1
      have hzero : g.natDegree = 0 := hg_deg
      rw [dif_pos hzero]
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
      exact Irreducible.natDegree_pos hq₀_irr
    have htarget : edf_l1 g d = (edf_correct_unconditional g d hm hsq hg_pos hd_pos hdeg).choose := by
      unfold edf_l1
      rw [dif_neg hg_deg]
      rw [dif_pos ⟨hm, hsq, hd_pos, hdeg⟩]
    rw [htarget]
    exact (edf_correct_unconditional g d hm hsq hg_pos hd_pos hdeg).choose_spec

-- ============================================================
-- §4. 端到端定理
-- ============================================================

/-- 使用 L1 翻译代码的 Zp 因式分解（精化定理填补后即验证 C++）。 -/
theorem factor_Zp_l1 (f : (ZMod p)[X]) (hf : f ≠ 0) :
    ∃ (lc : ZMod p) (factors : List ((ZMod p)[X] × ℕ)),
      FactorZpCorrect f lc factors :=
  factor_Zp_correct f hf
    sqfZp_l1 sqfZp_l1_correct
    ddf_l1 ddf_l1_correct
    edf_l1 edf_l1_correct

end L1Pipeline
