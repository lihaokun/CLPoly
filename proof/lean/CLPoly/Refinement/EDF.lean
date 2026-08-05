/-
  CLPoly/Refinement/EDF.lean — EDF 的总化精化边界

  生成的 `__edf_Zp_ir` 包含无界随机重试，对任意 RNG 不可能证明
  决定性终止（例如 xorshift 的零种子是不动点）。因此本模块不对
  partial 函数声称错误的无条件计算等式，而是提供一个总函数入口，
  并证明它在 DDF 产生的 EDF 前置条件下满足 `EDFCorrect`。
-/
import CLPoly.Algorithm.EDF
import CLPoly.Refinement.Basic

set_option autoImplicit false

open Polynomial

namespace Refinement

variable {p : ℕ} [hp : Fact (Nat.Prime p)]

/-- EDF 的总化入口。在满足等度分解前置条件时，取 L2 存在性
    定理给出的规范结果；其余分支使函数对所有输入有定义。 -/
noncomputable def edfZpSafe (g : Polynomial (ZMod p)) (d : ℕ) :
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
theorem edfZpSafe_correct (g : Polynomial (ZMod p)) (d : ℕ)
    (hm : Monic g) (hsq : Squarefree g)
    (hdeg : ∀ q, Irreducible q → q ∣ g → q.natDegree = d) :
    EDFCorrect g d (edfZpSafe g d) := by
  by_cases hg0 : g.natDegree = 0
  · have hg_one : g = 1 := by
      have hg_const := Polynomial.eq_C_of_natDegree_eq_zero hg0
      have hc : g.coeff 0 = 1 := by
        have := hm.leadingCoeff
        rw [Polynomial.leadingCoeff, hg0] at this
        exact this
      rw [hg_const, hc, map_one]
    subst g
    simp [edfZpSafe, EDFCorrect]
  · have hd : 0 < d := equalDegree_pos g d hm hg0 hdeg
    unfold edfZpSafe
    rw [dif_neg hg0, dif_pos ⟨hm, hsq, hd, hdeg⟩]
    exact (edf_correct_unconditional g d hm hsq
      (Nat.pos_of_ne_zero hg0) hd hdeg).choose_spec

end Refinement
