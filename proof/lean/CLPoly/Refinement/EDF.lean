/-
  CLPoly/Refinement/EDF.lean — EDF 的总化精化边界

  生成的 `__edf_Zp_ir` 包含无界随机重试，对任意 RNG 不可能证明
  决定性终止（例如 xorshift 的零种子是不动点）。因此本模块不对
  partial 函数声称错误的无条件计算等式，而是提供一个总函数入口，
  并证明它在 DDF 产生的 EDF 前置条件下满足 `EDFCorrect`。
-/
import CLPoly.Algorithm.EDF
import CLPoly.Generated.EDFSafe
import CLPoly.Refinement.Basic
import CLPoly.Math.Univariate

set_option autoImplicit false

open Polynomial
open CLPoly.Math

namespace Refinement

variable {p : ℕ} [hp : Fact (Nat.Prime p)]

/-- EDF 的总化入口。在满足等度分解前置条件时，取 L2 存在性
    定理给出的规范结果；其余分支使函数对所有输入有定义。 -/
noncomputable def edfCertified (g : Polynomial (ZMod p)) (d : ℕ) :
    List (Polynomial (ZMod p)) := by
  classical
  exact
    if hg0 : g.natDegree = 0 then []
    else if hpre : Monic g ∧ Squarefree g ∧ 0 < d ∧
        (∀ q : Polynomial (ZMod p), Irreducible q → q ∣ g → q.natDegree = d) then
      (edf_correct_unconditional g d hpre.1 hpre.2.1
        (Nat.pos_of_ne_zero hg0) hpre.2.2.1 hpre.2.2.2).choose
    else
      [g]

private lemma equalDegree_pos (g : Polynomial (ZMod p)) (d : ℕ)
    (hm : Monic g) (hg0 : g.natDegree ≠ 0)
    (hdeg : ∀ q, Irreducible q → q ∣ g → q.natDegree = d) : 0 < d := by
  have hg_nu : ¬ IsUnit g := fun hu =>
    hg0 (natDegree_eq_zero_of_isUnit hu)
  obtain ⟨q, hq_irr, hq_dvd⟩ :=
    WfDvdMonoid.exists_irreducible_factor hg_nu (Monic.ne_zero hm)
  rw [← hdeg q hq_irr hq_dvd]
  exact (Polynomial.natDegree_pos_iff_degree_pos.mpr hq_irr.degree_pos)

/-- `edfZpSafe` 在 DDF 保证的首一、无平方、等度输入上满足 EDF 规约。 -/
theorem edfCertified_correct (g : Polynomial (ZMod p)) (d : ℕ)
    (hm : Monic g) (hsq : Squarefree g)
    (hdeg : ∀ q, Irreducible q → q ∣ g → q.natDegree = d) :
    EDFCorrect g d (edfCertified g d) := by
  by_cases hg0 : g.natDegree = 0
  · have hg_one : g = 1 := by
      have hg_const := Polynomial.eq_C_of_natDegree_eq_zero hg0
      have hc : g.coeff 0 = 1 := by
        have := hm.leadingCoeff
        rw [Polynomial.leadingCoeff, hg0] at this
        exact this
      rw [hg_const, hc, map_one]
    subst g
    simp [edfCertified, EDFCorrect]
  · have hd : 0 < d := equalDegree_pos g d hm hg0 hdeg
    unfold edfCertified
    rw [dif_neg hg0, dif_pos ⟨hm, hsq, hd, hdeg⟩]
    exact (edf_correct_unconditional g d hm hsq
      (Nat.pos_of_ne_zero hg0) hd hdeg).choose_spec

/-- Map a successful bounded L1 execution back to the L2 polynomial model. -/
noncomputable def edfTryResultToPoly (p : ℕ) (result : Array SparsePolyZp × Rng) :
    List (Polynomial (ZMod p)) :=
  (result.1.map (SparsePolyZp.toPoly p)).toList

/-- Execute the fuel-bounded generated EDF path. A successful result is used
only when it satisfies the L2 contract; exhaustion or a rejected result takes
the certified total fallback. -/
noncomputable def edfZpSafe (g : Polynomial (ZMod p)) (f : SparsePolyZp)
    (d : ℕ) (fuel : ℕ) (rng : Rng) : List (Polynomial (ZMod p)) := by
  classical
  exact match Generated.edfZpTrySafe fuel #[] f d.toUInt64 rng with
    | none => edfCertified g d
    | some candidate =>
      let result := edfTryResultToPoly p candidate
      if EDFCorrect g d result then result else edfCertified g d

/-- The bounded L1 path plus certified fallback always meets the EDF contract
on the equal-degree inputs supplied by DDF. -/
theorem edfZpSafe_correct (g : Polynomial (ZMod p)) (f : SparsePolyZp)
    (d fuel : ℕ) (rng : Rng)
    (hm : Monic g) (hsq : Squarefree g)
    (hdeg : ∀ q, Irreducible q → q ∣ g → q.natDegree = d) :
    EDFCorrect g d (edfZpSafe g f d fuel rng) := by
  classical
  unfold edfZpSafe
  split
  · exact edfCertified_correct g d hm hsq hdeg
  · dsimp only
    split
    · assumption
    · exact edfCertified_correct g d hm hsq hdeg

/-- A bounded generated result that passes the EDF contract is the observable
result of the safe entry point (the fallback is not taken). -/
theorem edfZpSafe_of_try_correct (g : Polynomial (ZMod p)) (f : SparsePolyZp)
    (d fuel : ℕ) (rng rng' : Rng) (out : Array SparsePolyZp)
    (htry : Generated.edfZpTrySafe fuel #[] f d.toUInt64 rng = some (out, rng'))
    (hout : EDFCorrect g d
      ((out.map (SparsePolyZp.toPoly p)).toList)) :
    edfZpSafe g f d fuel rng = (out.map (SparsePolyZp.toPoly p)).toList := by
  classical
  unfold edfZpSafe
  rw [htry]
  dsimp [edfTryResultToPoly]
  exact if_pos hout

/-- Fuel exhaustion selects the certified fallback. -/
theorem edfZpSafe_of_try_exhausted (g : Polynomial (ZMod p)) (f : SparsePolyZp)
    (d fuel : ℕ) (rng : Rng)
    (htry : Generated.edfZpTrySafe fuel #[] f d.toUInt64 rng = none) :
    edfZpSafe g f d fuel rng = edfCertified g d := by
  classical
  unfold edfZpSafe
  rw [htry]

end Refinement
