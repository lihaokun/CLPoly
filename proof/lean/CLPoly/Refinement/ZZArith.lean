import CLPoly.Model
import CLPoly.Model
import CLPoly.Generated.Corpus
import CLPoly.Refinement.Basic
import CLPoly.Math.Univariate
import Mathlib.Tactic

set_option autoImplicit false

open Polynomial
open CLPoly.Math

namespace Refinement

variable {p : ℕ} [hp : Fact (Nat.Prime p)]

theorem __symmetric_mod_ir_refines (p : ℕ) [hp : Fact (Nat.Prime p)]
    (a : ZZ) (m : ZZ)
    :   Generated.__symmetric_mod_ir a m = ZZ.symmetricMod a m :=
by
  unfold Generated.__symmetric_mod_ir ZZ.symmetricMod
  dsimp
  have h_fdiv_r : ZZ.fdiv_r 0 a m = a.fmod m := rfl
  have h_fdiv_q' : ZZ.fdiv_q 0 m 2 = Int.fdiv m 2 := by
    unfold ZZ.fdiv_q; simp
  rw [h_fdiv_r, h_fdiv_q']
  rw [Int.fdiv_eq_ediv_of_nonneg m (by norm_num : 0 ≤ (2 : Int))]
  set r := a.fmod m with hr
  by_cases h : r * 2 ≤ m
  · have h_not_gt : ¬ (r > m / 2) := by omega
    simp [h, h_not_gt]
  · have h_gt : r > m / 2 := by omega
    simp [h, h_gt]

/- Binomial and integer-square-root refinement declarations were removed until
   direct generated execution proofs are available; no placeholder theorem is
   exported from this module. -/

end Refinement
