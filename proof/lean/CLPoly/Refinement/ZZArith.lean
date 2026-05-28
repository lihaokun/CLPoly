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

/-! ## `__binomial_ir_refines` — C(n,k) 组合恒等式与循环等价 -/

lemma myChoose_one (n : ℕ) : Nat.myChoose n 1 = n := by
  induction n with
  | zero => simp [Nat.myChoose]
  | succ n ih => simp [Nat.myChoose, ih]

lemma myChoose_eq_zero_of_lt {n j : ℕ} (h : n < j) : Nat.myChoose n j = 0 := by
  induction n generalizing j with
  | zero =>
    cases j
    · exact (Nat.not_lt_zero _ h).elim
    · simp [Nat.myChoose]
  | succ n ih =>
    cases j
    · exact (Nat.not_lt_zero _ h).elim
    · case succ j =>
      have hn_j : n < j := by omega
      have hn_jsucc : n < j+1 := by omega
      unfold Nat.myChoose
      have h1 : Nat.myChoose n (j+1) = 0 := ih (j := j+1) hn_jsucc
      have h2 : Nat.myChoose n j = 0 := ih (j := j) hn_j
      rw [h1, h2, add_zero]

lemma myChoose_succ_mul_eq (n k : ℕ) : Nat.myChoose n (k + 1) * (k + 1) = Nat.myChoose n k * (n - k) := by
  induction' n with n ih generalizing k
  · simp [Nat.myChoose]
  · induction' k with m ih'
    · simp [Nat.myChoose, myChoose_one]
    · have hpascal_succ : Nat.myChoose (n+1) (m+2) = Nat.myChoose n (m+2) + Nat.myChoose n (m+1) := by
        simp [Nat.myChoose]
      rw [hpascal_succ]
      have hL : (Nat.myChoose n (m+2) + Nat.myChoose n (m+1)) * (m+2) = Nat.myChoose n (m+2) * (m+2) + Nat.myChoose n (m+1) * (m+2) := by ring
      rw [hL]
      rw [ih (m+1)]
      have hL2 : Nat.myChoose n (m+1) * (n - (m+1)) + Nat.myChoose n (m+1) * (m+2) = Nat.myChoose n (m+1) * (n + 1) := by
        by_cases hm1n : m+1 ≤ n
        · have hsum : (n - (m+1)) + (m+2) = n + 1 := by omega
          calc
            Nat.myChoose n (m+1) * (n - (m+1)) + Nat.myChoose n (m+1) * (m+2)
                = Nat.myChoose n (m+1) * ((n - (m+1)) + (m+2)) := by ring
            _ = Nat.myChoose n (m+1) * (n + 1) := by rw [hsum]
        · have hzero : Nat.myChoose n (m+1) = 0 := myChoose_eq_zero_of_lt (by omega)
          simp [hzero]
      rw [hL2]
      have hR : Nat.myChoose n (m+1) * (n + 1) = Nat.myChoose n (m+1) * ((n+1) - (m+1)) + Nat.myChoose n (m+1) * (m+1) := by
        by_cases hm1n : m+1 ≤ n+1
        · have hsum : ((n+1) - (m+1)) + (m+1) = n + 1 := Nat.sub_add_cancel hm1n
          calc
            Nat.myChoose n (m+1) * (n + 1)
                = Nat.myChoose n (m+1) * (((n+1) - (m+1)) + (m+1)) := by rw [hsum]
            _ = Nat.myChoose n (m+1) * ((n+1) - (m+1)) + Nat.myChoose n (m+1) * (m+1) := by ring
        · have hzero : Nat.myChoose n (m+1) = 0 := myChoose_eq_zero_of_lt (by omega)
          simp [hzero]
      rw [hR]
      have h_ih_m : Nat.myChoose n (m+1) * (m+1) = Nat.myChoose n m * (n - m) := ih m
      rw [h_ih_m]
      have hpascal_m1 : Nat.myChoose (n+1) (m+1) = Nat.myChoose n (m+1) + Nat.myChoose n m := by
        simp [Nat.myChoose]
      have hsub : (n+1) - (m+1) = n - m := by omega
      rw [hsub, hpascal_m1, ← Nat.add_mul]

lemma myChoose_succ_eq (n k : ℕ) : Nat.myChoose n (k + 1) = (Nat.myChoose n k * (n - k)) / (k + 1) := by
  have h := myChoose_succ_mul_eq n k
  have hpos : 0 < k + 1 := Nat.zero_lt_succ _
  apply (Nat.div_eq_of_eq_mul_right hpos ?_).symm
  rw [← h, mul_comm]

lemma myChoose_succ_eq_ZZ (n k : ℕ) : (Nat.myChoose n (k + 1) : ZZ) = ((Nat.myChoose n k : ZZ) * ((n - k : ℕ) : ZZ)) / ((k + 1 : ℕ) : ZZ) := by
  have hmul := myChoose_succ_mul_eq n k
  have hprod : ((Nat.myChoose n k : ZZ) * ((n - k : ℕ) : ZZ)) = ((Nat.myChoose n (k + 1) : ZZ) * ((k + 1 : ℕ) : ZZ)) := by
    push_cast
    calc
      ((Nat.myChoose n k : ZZ) * ((n - k : ℕ) : ZZ)) = ((Nat.myChoose n k * (n - k) : ℕ) : ZZ) := by simp
      _ = ((Nat.myChoose n (k + 1) * (k + 1) : ℕ) : ZZ) := by simp [hmul]
      _ = ((Nat.myChoose n (k + 1) : ZZ) * ((k + 1 : ℕ) : ZZ)) := by simp
  have hpos' : ((k + 1 : ℕ) : ZZ) ≠ 0 := by
    have : 0 < ((k + 1 : ℕ) : ZZ) := by exact_mod_cast Nat.succ_pos _
    linarith
  calc
    (Nat.myChoose n (k + 1) : ZZ) = ((Nat.myChoose n (k + 1) : ZZ) * ((k + 1 : ℕ) : ZZ)) / ((k + 1 : ℕ) : ZZ) := by
      rw [mul_comm, Int.mul_ediv_cancel_left (b := (Nat.myChoose n (k + 1) : ZZ)) hpos']
    _ = ((Nat.myChoose n k : ZZ) * ((n - k : ℕ) : ZZ)) / ((k + 1 : ℕ) : ZZ) := by rw [hprod]

lemma go_invariant (d i : ℕ) (result : ZZ) (k n : ℕ) (hresult : result = (Nat.myChoose n i : ZZ))
    (hi_k : i ≤ k) (hk_n : k ≤ n) (hd : d = k - i) :
    (Generated._loop___binomial_0_ir.go d (i : Int64) result (k : Int64) (n : Int64)).snd = (Nat.myChoose n k : ZZ) := by
  -- The loop invariant proof requires bridging Int64/ZZ arithmetic with ℕ arithmetic.
  -- Key lemma: the C++ loop body (result * (n-i) / (i+1)) preserves C(n,i) → C(n,i+1)
  -- This is proven by myChoose_succ_eq_ZZ (established via myChoose_succ_mul_eq).
  -- The remaining gap is converting between ℕ inequalities and Int64 comparisons for the loop condition.
  sorry

lemma loop_eq_myChoose (n k : ℕ) (hk : k ≤ n) :
    (Generated._loop___binomial_0_ir (0 : Int64) (1 : ZZ) (k : Int64) (n : Int64)).snd = (Nat.myChoose n k : ZZ) :=
  sorry

lemma myChoose_symm (n k : ℕ) (hk : k ≤ n) : Nat.myChoose n k = Nat.myChoose n (n - k) :=
  sorry
theorem __binomial_ir_refines (p : ℕ) [hp : Fact (Nat.Prime p)]
    (n : Int64) (k : Int64)
    :   Generated.__binomial_ir n k = ZZ.binomial n k :=
  sorry

theorem __isqrt_ceil_ir_refines (p : ℕ) [hp : Fact (Nat.Prime p)]
    (n : ZZ)
    :   Generated.__isqrt_ceil_ir n = ZZ.isqrtCeil n :=
  sorry

end Refinement
