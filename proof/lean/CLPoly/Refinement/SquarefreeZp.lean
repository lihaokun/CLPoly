/-
  CLPoly/Refinement/SquarefreeZp.lean — __squarefree_Zp_ir（C++） ≃ sqfZp（L2）

  结构对应：
    C++（SparsePolyZp）                   L2（Polynomial (ZMod p)）
    ──────────────────────────────────────────────────────
    derivative f = ∅                     derivative f = 0
      → __extract_pth_root_ir f           → Polynomial.contract p f
      → __upoly_make_monic_ir             → normalize
      → __squarefree_Zp_ir (递归)         → sqfZp (递归)
      → scale exponents × p               → .map (pr.1, pr.2 * p)

    derivative f ≠ ∅                     derivative f ≠ 0
      → polynomial_GCD f f'              → gcd f f'
      → pair_vec_div f / gcd             → f /ₘ gcd
      → normalization                    → normalize
      → yunLoop                          → yunLoop
      → c_rem (rest)                     → c_rem
        → __extract_pth_root_ir c_rem      → Polynomial.contract p c_rem
        → __upoly_make_monic_ir             → normalize
        → __squarefree_Zp_ir (递归)         → sqfZp (递归)
        → scale exponents × p               → .map (pr.1, pr.2 * p)

  证明策略：强归纳于 f.natDegree。每一步用 dsimp 展开 partial def，
  用 L2 已证明的引理（yunLoop 正确性、natDegree 递减、contract p 等）。
-/

import CLPoly.Algorithm.SquarefreeZp
import CLPoly.Model
import CLPoly.Generated.Corpus
import CLPoly.Refinement.Basic
import CLPoly.Refinement.ZpArith
import CLPoly.Math.Univariate

set_option autoImplicit false

open Polynomial
open CLPoly.Math

namespace Refinement

variable {p : ℕ} [hp : Fact (Nat.Prime p)]

-- ============================================================
-- §1. 辅助引理：SparsePolyZp 操作 ↔ Polynomial (ZMod p) 操作
-- ============================================================

/-- 辅助引理：extGcdAux 的结果系数关于 (old_s, s) 线性，且 g 整除 old_r 和 r。
    存在 x, y, g 使 extGcdAux old_r r old_s s = (g, x*old_s + y*s) 且
    x*old_r + y*r = g 且 g ∣ old_r 且 g ∣ r。
    证明：强归纳于 r.natAbs，扩展 extGcdAux_linearity。 -/
private lemma extGcdAux_linearity_gcd (old_r r old_s s : Int) :
    ∃ (x y g : Int), Zp.extGcdAux old_r r old_s s = (g, x * old_s + y * s) ∧ x * old_r + y * r = g ∧ g ∣ old_r ∧ g ∣ r := by
  have h_nat : ∀ (n : ℕ), (∀ (m : ℕ), m < n → ∀ (old_r' r' old_s' s' : Int), r'.natAbs = m →
    ∃ (x y g : Int), Zp.extGcdAux old_r' r' old_s' s' = (g, x * old_s' + y * s') ∧ x * old_r' + y * r' = g ∧ g ∣ old_r' ∧ g ∣ r') →
    ∀ (old_r r old_s s : Int), r.natAbs = n → ∃ (x y g : Int), Zp.extGcdAux old_r r old_s s = (g, x * old_s + y * s) ∧ x * old_r + y * r = g ∧ g ∣ old_r ∧ g ∣ r := by
    intro n IH old_r r old_s s hn
    by_cases hzero : (r : Int) = 0
    · refine ⟨1, 0, old_r, ?_, ?_, ?_, ?_⟩
      · simp [Zp.extGcdAux, hzero]
      · simp
      · exact dvd_refl old_r
      · subst hzero; simp
    · have hB_ne_zero : r ≠ 0 := hzero
      have h_sub : Zp.extGcdAux old_r r old_s s = Zp.extGcdAux r (old_r % r) s (old_s - (old_r / r) * s) := by
        have hzero_beq : ¬ (r == 0) := by
          intro h
          apply hzero
          simpa using h
        rw [Zp.extGcdAux, dif_neg hzero_beq]
        have h_eq : old_r - (old_r / r) * r = old_r % r := by
          linarith [Int.ediv_add_emod old_r r]
        rw [← h_eq]
      rw [h_sub]
      let q := old_r / r
      let r' := old_r % r
      have h_decr : r'.natAbs < r.natAbs := by
        have h_nonneg : 0 ≤ r' := Int.emod_nonneg old_r hB_ne_zero
        have h_mod_lt : r' < (r.natAbs : Int) := Int.emod_lt old_r hB_ne_zero
        have h_natAbs_lt : r'.natAbs < r.natAbs := by
          have h_eq : (r'.natAbs : Int) = r' := Int.natAbs_of_nonneg h_nonneg
          have h' : (r'.natAbs : Int) < (r.natAbs : Int) := by rw [h_eq]; exact h_mod_lt
          exact (Int.ofNat_lt.mp h')
        exact h_natAbs_lt
      have h_r'_lt_n : r'.natAbs < n := by rw [← hn]; exact h_decr
      rcases IH r'.natAbs h_r'_lt_n r r' s (old_s - q * s) rfl with ⟨x', y', g', h_eq_inner, h_bezout_inner, h_gcd1, h_gcd2⟩
      -- h_bezout_inner: x'*r + y'*r' = g'
      -- h_gcd1: g' ∣ r
      -- h_gcd2: g' ∣ r'
      let x := y'
      let y := x' - q * y'
      let g := g'
      have h_eq_outer : Zp.extGcdAux r r' s (old_s - q * s) = (g, x * old_s + y * s) := by
        rw [h_eq_inner]
        dsimp [x, y, g]
        ring
      have h_bezout_outer : x * old_r + y * r = g := by
        dsimp [x, y, g]
        have h_r'_eq : r' = old_r - q * r := by
          linarith [Int.ediv_add_emod old_r r]
        rw [h_r'_eq] at h_bezout_inner
        linarith
      -- g ∣ old_r: 因为 g ∣ r 且 g ∣ r' = old_r - q*r，所以 g ∣ old_r
      have h_gcd_old_r : g' ∣ old_r := by
        have h_gcd_r' : g' ∣ r' := h_gcd2
        have h_gcd_r : g' ∣ r := h_gcd1
        have h_r'_eq : r' = old_r - q * r := by
          linarith [Int.ediv_add_emod old_r r]
        have h_gcd_r'_expr : g' ∣ old_r - q * r := by
          rw [h_r'_eq] at h_gcd2
          exact h_gcd2
        have : g' ∣ old_r := by
          have : g' ∣ q * r := by
            -- g' ∣ r → g' ∣ q*r
            exact h_gcd_r.mul_left q
          -- g' ∣ (old_r - q*r) + q*r = old_r
          have : g' ∣ (old_r - q * r) + q * r := dvd_add h_gcd_r'_expr this
          simpa [sub_add_cancel] using this
        exact this
      exact ⟨x, y, g, h_eq_outer, h_bezout_outer, h_gcd_old_r, h_gcd1⟩
  have h_all : ∀ (n : ℕ), ∀ (old_r r old_s s : Int), r.natAbs = n →
    ∃ (x y g : Int), Zp.extGcdAux old_r r old_s s = (g, x * old_s + y * s) ∧ x * old_r + y * r = g ∧ g ∣ old_r ∧ g ∣ r := by
    intro n
    induction n using Nat.strong_induction_on with
    | h m IH_all =>
      intro old_r r old_s s hn
      apply h_nat m (λ k hk old_r' r' old_s' s' h_eq => ?_) old_r r old_s s hn
      exact IH_all k (Nat.lt_of_lt_of_le hk (by omega)) old_r' r' old_s' s' h_eq
  apply h_all r.natAbs old_r r old_s s rfl

/-- Zp.extGcdAux 的 Bezout 恒等式：对于 (g, s') = Zp.extGcdAux A B 0 1，
    存在整数 u 使 s' * B + u * A = g。 -/
private lemma extGcdAux_bezout (A B : Int) :
    let (g, s') := Zp.extGcdAux A B 0 1
    ∃ (u : Int), s' * B + u * A = g := by
  rcases extGcdAux_linearity_gcd A B 0 1 with ⟨x, y, g, h_eq, h_bezout, h_gcdA, h_gcdB⟩
  have h_coeff : (Zp.extGcdAux A B 0 1).2 = y := by
    calc
      (Zp.extGcdAux A B 0 1).2 = (g, x * (0 : Int) + y * (1 : Int)).2 := by rw [h_eq]
      _ = x * (0 : Int) + y * (1 : Int) := rfl
      _ = y := by ring
  have h_g_val : (Zp.extGcdAux A B 0 1).1 = g := by
    calc
      (Zp.extGcdAux A B 0 1).1 = (g, x * (0 : Int) + y * (1 : Int)).1 := by rw [h_eq]
      _ = g := rfl
  refine ⟨x, ?_⟩
  rw [h_coeff, h_g_val]
  linarith

/-- extGcdAux 对非负输入返回非负值的第一个分量。
    证明：强归纳于 B（第二个参数），利用 extGcdAux_linearity_gcd 连通不同 (old_s, s) 的 .1。 -/
private lemma extGcdAux_fst_independent (old_r r old_s s : Int) :
    (Zp.extGcdAux old_r r old_s s).1 = (Zp.extGcdAux old_r r 0 1).1 := by
  have h_nat : ∀ (n : ℕ), (∀ m < n, ∀ (old_r' r' old_s' s' : Int), r'.natAbs = m →
    (Zp.extGcdAux old_r' r' old_s' s').1 = (Zp.extGcdAux old_r' r' 0 1).1) →
    ∀ (old_r r old_s s : Int), r.natAbs = n → (Zp.extGcdAux old_r r old_s s).1 = (Zp.extGcdAux old_r r 0 1).1 := by
    intro n IH old_r r old_s s hn
    by_cases hzero : r == 0
    · simp [Zp.extGcdAux, hzero]
    · have hzero_beq : ¬ (r == 0) := hzero
      have h_sub1 : Zp.extGcdAux old_r r old_s s = Zp.extGcdAux r (old_r - (old_r / r) * r) s (old_s - (old_r / r) * s) := by
        rw [Zp.extGcdAux, dif_neg hzero_beq]
      have h_sub2 : Zp.extGcdAux old_r r 0 1 = Zp.extGcdAux r (old_r - (old_r / r) * r) 1 (-(old_r / r)) := by
        rw [Zp.extGcdAux, dif_neg hzero_beq]
        dsimp; simp
      have h_eq : old_r - (old_r / r) * r = old_r % r := by
        linarith [Int.ediv_add_emod old_r r]
      have h_decr : (old_r % r).natAbs < r.natAbs := by
        have h_nonneg : 0 ≤ old_r % r := Int.emod_nonneg old_r (by
          intro h; apply hzero; simp [h])
        have h_mod_lt : old_r % r < (r.natAbs : Int) := Int.emod_lt old_r (by
          intro h; apply hzero; simp [h])
        have h_natAbs_lt : (old_r % r).natAbs < r.natAbs := by
          have h_eq' : ((old_r % r).natAbs : Int) = old_r % r := Int.natAbs_of_nonneg h_nonneg
          have h' : ((old_r % r).natAbs : Int) < (r.natAbs : Int) := by rw [h_eq']; exact h_mod_lt
          exact (Int.ofNat_lt.mp h')
        exact h_natAbs_lt
      have h_lt_n : (old_r % r).natAbs < n := by rw [← hn]; exact h_decr
      have h_IH1 : (Zp.extGcdAux r (old_r % r) s (old_s - (old_r / r) * s)).1 = (Zp.extGcdAux r (old_r % r) 0 1).1 :=
        IH (old_r % r).natAbs h_lt_n r (old_r % r) s (old_s - (old_r / r) * s) rfl
      have h_IH2 : (Zp.extGcdAux r (old_r % r) 1 (-(old_r / r))).1 = (Zp.extGcdAux r (old_r % r) 0 1).1 :=
        IH (old_r % r).natAbs h_lt_n r (old_r % r) (1 : Int) (-(old_r / r)) rfl
      calc
        (Zp.extGcdAux old_r r old_s s).1 = (Zp.extGcdAux r (old_r - (old_r / r) * r) s (old_s - (old_r / r) * s)).1 := by rw [h_sub1]
        _ = (Zp.extGcdAux r (old_r % r) s (old_s - (old_r / r) * s)).1 := by rw [h_eq]
        _ = (Zp.extGcdAux r (old_r % r) 0 1).1 := h_IH1
        _ = (Zp.extGcdAux r (old_r % r) 1 (-(old_r / r))).1 := by symm; exact h_IH2
        _ = (Zp.extGcdAux r (old_r - (old_r / r) * r) 1 (-(old_r / r))).1 := by rw [← h_eq]
        _ = (Zp.extGcdAux old_r r 0 1).1 := by rw [← h_sub2]
  have h_all : ∀ (n : ℕ), ∀ (old_r r old_s s : Int), r.natAbs = n →
    (Zp.extGcdAux old_r r old_s s).1 = (Zp.extGcdAux old_r r 0 1).1 := by
    intro n
    induction n using Nat.strong_induction_on with
    | h m IH_all =>
      intro old_r r old_s s hm
      apply h_nat m (λ k hk old_r' r' old_s' s' hk_eq => IH_all k hk old_r' r' old_s' s' hk_eq) old_r r old_s s hm
  apply h_all r.natAbs old_r r old_s s rfl

private lemma extGcdAux_gcd_nonneg (A B : ℕ) : 0 ≤ (Zp.extGcdAux (A : Int) (B : Int) 0 1).1 := by
  revert A
  induction B using Nat.strong_induction_on with
  | h B IH =>
    intro A
    by_cases hB : (B : ℕ) = 0
    · subst hB; simp [Zp.extGcdAux]
    · have hB_int_ne_zero : (B : Int) ≠ 0 := by exact_mod_cast hB
      have hzero_beq : ¬ (B : Int) == 0 := by
        intro h; apply hB_int_ne_zero; simpa using h
      have h_sub : Zp.extGcdAux (A : Int) (B : Int) 0 1 =
          Zp.extGcdAux (B : Int) ((A : Int) % (B : Int)) (1 : Int) (-((A : Int) / (B : Int))) := by
        rw [Zp.extGcdAux, dif_neg hzero_beq]; dsimp
        have h_eq_rem : (A : Int) - ((A : Int) / (B : Int)) * (B : Int) = (A : Int) % (B : Int) := by
          linarith [Int.ediv_add_emod (A : Int) (B : Int)]
        simp [h_eq_rem]
      rw [h_sub]
      rw [extGcdAux_fst_independent (B : Int) ((A : Int) % (B : Int)) (1 : Int) (-((A : Int) / (B : Int)))]
      have h_nonneg : 0 ≤ (A : Int) % (B : Int) :=
        Int.emod_nonneg (A : Int) hB_int_ne_zero
      have h_lt_B : ((A : Int) % (B : Int)).natAbs < B := by
        have h_mod_lt : (A : Int) % (B : Int) < (B : Int) :=
          Int.emod_lt (A : Int) hB_int_ne_zero
        have h_natAbs_lt : ((A : Int) % (B : Int)).natAbs < (B : Int).natAbs := by
          have h_eq : (((A : Int) % (B : Int)).natAbs : Int) = (A : Int) % (B : Int) :=
            Int.natAbs_of_nonneg h_nonneg
          have h' : (((A : Int) % (B : Int)).natAbs : Int) < ((B : Int).natAbs : Int) := by
            rw [h_eq]; exact h_mod_lt
          exact (Int.ofNat_lt.mp h')
        have : (B : Int).natAbs = B := by simp
        rw [this] at h_natAbs_lt
        exact h_natAbs_lt
      -- IH generalizes over A, so we can apply it with A' = B
      have h_IH : ∀ (m : ℕ), m < B → ∀ (A' : ℕ), 0 ≤ (Zp.extGcdAux (A' : Int) (m : Int) 0 1).1 := IH
      have h_result := h_IH ((A : Int) % (B : Int)).natAbs h_lt_B B
      -- h_result: 0 ≤ (extGcdAux (B : Int) (((A:Int) % (B:Int)).natAbs : Int) 0 1).1
      -- but we need: 0 ≤ (extGcdAux (B : Int) ((A:Int) % (B:Int)) 0 1).1
      simpa [Int.natAbs_of_nonneg h_nonneg] using h_result

/-- 对于非零且 AllReduced 的 Zp 元素 a，在 ZMod p 中有 (a * a.inv).toZMod p = 1。
    证明：extGcdAux_bezout → gcd=1 → modInv = (s%p).toNat → 模算术。 -/
private lemma Zp_toZMod_inv_mul_self (a : Zp) (hred_a : Zp.Reduced p a) (hval_nonzero : a.val.toNat ≠ 0)
    (hp_lt_U64 : p < UInt64.size) : Zp.toZMod p (a * a.inv) = (1 : ZMod p) := by
  rcases hred_a with ⟨h_prime_eq, h_val_lt⟩
  have hp_prime : Nat.Prime p := hp.out
  have hp_pos : 0 < p := Nat.Prime.pos hp_prime
  -- Step 1: Bezout identity from extGcdAux
  rcases extGcdAux_bezout (p : Int) (a.val.toNat : Int) with ⟨u, h_bezout⟩
  rcases extGcdAux_linearity_gcd (p : Int) (a.val.toNat : Int) 0 1 with ⟨x, y, g', h_eq, h_bezout2, h_gcd_p, h_gcd_a⟩
  have h_g_val : (Zp.extGcdAux (p : Int) (a.val.toNat : Int) 0 1).1 = g' := by simpa [h_eq]
  have h_s_val : (Zp.extGcdAux (p : Int) (a.val.toNat : Int) 0 1).2 = y := by rw [h_eq]; simp
  -- Step 2: gcd=1 via primality and bounds
  have h_g'_nonneg : 0 ≤ g' := by
    have h_nonneg : 0 ≤ (Zp.extGcdAux (p : Int) (a.val.toNat : Int) 0 1).1 :=
      extGcdAux_gcd_nonneg p (a.val.toNat)
    rw [← h_g_val]; exact h_nonneg
  have h_g'_eq_nat : (g'.toNat : ℤ) = g' := Int.toNat_of_nonneg h_g'_nonneg
  have h_gcd_nat_dvd_p : g'.toNat ∣ p := by
    have h_dvd_int : (g'.toNat : ℤ) ∣ (p : ℤ) := by rw [h_g'_eq_nat]; exact h_gcd_p
    exact (Int.ofNat_dvd (m := g'.toNat) (n := p)).mp h_dvd_int
  have h_gcd_nat_dvd_a : g'.toNat ∣ a.val.toNat := by
    have h_dvd_int : (g'.toNat : ℤ) ∣ (a.val.toNat : ℤ) := by rw [h_g'_eq_nat]; exact h_gcd_a
    exact (Int.ofNat_dvd (m := g'.toNat) (n := a.val.toNat)).mp h_dvd_int
  have h_gcd_eq_one_nat : g'.toNat = 1 := by
    rcases hp_prime.eq_one_or_self_of_dvd g'.toNat h_gcd_nat_dvd_p with (h | h)
    · exact h
    · have h_a_dvd : p ∣ a.val.toNat := by rw [← h]; exact h_gcd_nat_dvd_a
      have h_le : p ≤ a.val.toNat := Nat.le_of_dvd (by omega) h_a_dvd
      omega
  have h_gcd_eq_one : g' = 1 := by
    rw [← h_g'_eq_nat, h_gcd_eq_one_nat, Nat.cast_one]
  -- Step 3: y * a.val + u * p = 1
  have h_bezout_one : y * (a.val.toNat : Int) + u * (p : Int) = 1 := by
    simpa [h_g_val, h_s_val, h_gcd_eq_one] using h_bezout
  -- Step 4: Let r = (y % p).toNat. Then r * a.val ≡ 1 (mod p)
  have h_emod_nonneg : 0 ≤ y % (p : ℤ) := Int.emod_nonneg y (by exact mod_cast hp_pos.ne')
  have h_emod_lt : y % (p : ℤ) < (p : ℤ) := Int.emod_lt y (by exact mod_cast hp_pos.ne')
  set r := (y % (p : ℤ)).toNat with hr
  have h_r_int : (r : ℤ) = y % (p : ℤ) := by
    rw [hr, Int.toNat_of_nonneg h_emod_nonneg]
  -- From y = (y/p)*p + r, rewrite h_bezout_one to express r * a.val % p = 1 % p
  have h_y_eq : y = (y / (p : ℤ)) * (p : ℤ) + y % (p : ℤ) := by
    simpa [mul_comm] using (Int.ediv_add_emod y (p : ℤ)).symm
  rw [h_y_eq] at h_bezout_one
  have h_mod_int : (r : ℤ) * (a.val.toNat : ℤ) % (p : ℤ) = (1 : ℤ) % (p : ℤ) := by
    -- ((y/p)*p + r) * a.val + u * p = 1
    -- → r * a.val + p * (...) = 1  →  r * a.val = 1 - p * (...)
    -- → r * a.val % p = 1 % p
    have h_bezout' : (r : ℤ) * (a.val.toNat : ℤ) = 1 - (u + (y / (p : ℤ)) * (a.val.toNat : ℤ)) * (p : ℤ) := by
      rw [← h_r_int] at h_bezout_one
      nlinarith
    rw [h_bezout']
    -- Goal: (1 - K * p) % p = 1 % p
    set K := u + (y / (p : ℤ)) * (a.val.toNat : ℤ) with hK
    have h_mod_neg_zero : (-(K * (p : ℤ))) % (p : ℤ) = 0 := by
      calc
        (-(K * (p : ℤ))) % (p : ℤ) = ((-K) * (p : ℤ)) % (p : ℤ) := by ring
        _ = 0 := Int.emod_eq_zero_of_dvd ⟨-K, by ring⟩
    calc
      (1 - K * (p : ℤ)) % (p : ℤ) = ((1 : ℤ) % (p : ℤ) + (-(K * (p : ℤ))) % (p : ℤ)) % (p : ℤ) := by
        simp [Int.sub_eq_add_neg, Int.add_emod]
      _ = ((1 : ℤ) % (p : ℤ) + 0) % (p : ℤ) := by rw [h_mod_neg_zero]
      _ = (1 : ℤ) % (p : ℤ) := by simp
  -- Step 5: Convert to ℕ: r * a.val.toNat % p = 1
  have h_mod_nat : r * a.val.toNat % p = 1 := by
    have h_eq_int : ((r * a.val.toNat % p : ℕ) : ℤ) = ((1 : ℕ) : ℤ) := by
      calc
        ((r * a.val.toNat % p : ℕ) : ℤ) = ((r : ℤ) * (a.val.toNat : ℤ)) % (p : ℤ) := by simp
        _ = (1 : ℤ) % (p : ℤ) := h_mod_int
        _ = (1 : ℤ) := by
          have h1_lt_p : (1 : ℤ) < (p : ℤ) := by exact mod_cast (Nat.Prime.one_lt hp_prime)
          simp [Int.emod_eq_of_lt (by omega) h1_lt_p]
    exact_mod_cast h_eq_int
  -- Step 6: Connect modInv: (Zp.modInv a.val a.prime).toNat = r
  have h_modInv_toNat : (Zp.modInv a.val a.prime).toNat = r := by
    -- modInv defined with `let` bindings. Expand via dsimp, then substitute y.
    unfold Zp.modInv
    have hzero_val : a.val ≠ 0 := by
      intro h; apply hval_nonzero; simp [h]
    have hzero_bool : (a.val == 0) = false := by simp [hzero_val]
    rw [hzero_bool]
    -- Expand `let (_, s) := extGcdAux ...` and subsequent `let` bindings
    dsimp
    -- Now the goal has (extGcdAux ...).2, replace it with y (by h_s_val)
    have h_ext2 : (Zp.extGcdAux (a.prime.toNat : Int) (a.val.toNat : Int) 0 1).2 = y := by
      rw [h_prime_eq]; exact h_s_val
    rw [h_ext2]
    -- Goal: ((if y % (a.prime.toNat : Int) < 0 then ... else ...).toNat.toUInt64).toNat = r
    -- Since y % p ≥ 0 and a.prime.toNat = p, the `if` is false
    have h_nonneg : 0 ≤ y % (a.prime.toNat : Int) := by rw [h_prime_eq]; exact h_emod_nonneg
    split_ifs with h_lt
    · exfalso; linarith
    · calc
        (y % (a.prime.toNat : Int)).toNat.toUInt64.toNat = (y % (a.prime.toNat : Int)).toNat := by
          rw [h_prime_eq]
          have h_bound : (y % (p : ℤ)).toNat < UInt64.size := by
            have h_lt_p : (y % (p : ℤ)).toNat < p := by
              have h_lt_int : y % (p : ℤ) < (p : ℤ) := h_emod_lt
              have h_nonneg_int : 0 ≤ y % (p : ℤ) := h_emod_nonneg
              omega
            have hp_ge_2 : 2 ≤ p := Nat.Prime.two_le hp_prime
            have hp_lt_U64' : p < UInt64.size := hp_lt_U64
            omega
          simpa [h_bound]
        _ = (y % (p : ℤ)).toNat := by
          rw [show (a.prime.toNat : Int) = (p : Int) from by exact_mod_cast h_prime_eq]
        _ = r := by rw [← hr]
  have h_mul_val : (a * a.inv).val.toNat = 1 := by
    have h_mul_def : (a * a.inv).val = ((a.val.toNat * a.inv.val.toNat) % a.prime.toNat).toUInt64 := rfl
    rw [h_mul_def]
    simp [Zp.inv, h_modInv_toNat, h_prime_eq, h_mod_nat, mul_comm]
  -- Step 8: toZMod converts value modulo p
  calc
    Zp.toZMod p (a * a.inv) = ((a * a.inv).val.toNat : ZMod p) := rfl
    _ = (1 : ZMod p) := by rw [h_mul_val]; simp

/-- 非零约化系数的实现层逆元确实令原始 `val` 等于 1。
这是 `divmod` 合法分支和 GCD 递归不变量所需的表示层版本。 -/
private lemma Zp_inv_mul_val_eq_one (a : Zp) (hred_a : Zp.Reduced p a)
    (hval_nonzero : a.val.toNat ≠ 0) (hp_lt_U64 : p < UInt64.size) :
    (a.inv * a).val = 1 := by
  have hp_pos : 0 < p := hp.out.pos
  have ha_prime : a.prime.toNat = p := hred_a.1
  have hainv_prime : a.inv.prime.toNat = p := by simpa [Zp.inv] using ha_prime
  have hto : Zp.toZMod p (a * a.inv) = (1 : ZMod p) :=
    Zp_toZMod_inv_mul_self a hred_a hval_nonzero hp_lt_U64
  have hmulred := SparsePolyZp.Zp_mul_reduced a a.inv (by simpa [ha_prime] using hp_pos)
  have hval_lt : (a * a.inv).val.toNat < p := by
    rw [UInt64.lt_iff_toNat_lt] at hmulred
    simpa [ha_prime] using hmulred.1
  have hcast : ((a * a.inv).val.toNat : ZMod p) = 1 := by simpa [Zp.toZMod] using hto
  have hval_nat : (a * a.inv).val.toNat = 1 := by
    letI : NeZero p := ⟨hp.out.ne_zero⟩
    letI : Fact (1 < p) := ⟨hp.out.one_lt⟩
    have hv := congrArg ZMod.val hcast
    calc
      (a * a.inv).val.toNat = ZMod.val ((a * a.inv).val.toNat : ZMod p) := by
        rw [ZMod.val_natCast, Nat.mod_eq_of_lt hval_lt]
      _ = ZMod.val (1 : ZMod p) := hv
      _ = 1 := ZMod.val_one p
  have hraw : (a * a.inv).val = 1 := UInt64.toNat_inj.mp (by simpa using hval_nat)
  change (((a.inv.val.toNat * a.val.toNat) % a.inv.prime.toNat).toUInt64) = 1
  change (((a.val.toNat * a.inv.val.toNat) % a.prime.toNat).toUInt64) = 1 at hraw
  simpa [Zp.inv, Nat.mul_comm] using hraw

/-- `modInv a q` 的返回值恒落在 `[0, q)`。由定义构造直接得出（无需 Bezout / 逆元正确性）：
    `a = 0` 时返回 0；否则返回 `s % q`（`Int.emod` 值域 `[0, q)`），再经 toUInt64 往返保值。 -/
private lemma modInv_val_lt (a q : UInt64) (hq : 1 < q.toNat) (hq_size : q.toNat ≤ UInt64.size) :
    (Zp.modInv a q).toNat < q.toNat := by
  unfold Zp.modInv
  by_cases hb : (a == 0) = true
  · rw [if_pos hb]; simp; omega
  · rw [if_neg hb]
    dsimp only
    set s : Int := (Zp.extGcdAux (q.toNat : Int) (a.toNat : Int) 0 1).2 with hs
    have hne : (q.toNat : Int) ≠ 0 := by
      have : q.toNat ≠ 0 := by omega
      exact_mod_cast this
    have hpos : (0 : Int) < (q.toNat : Int) := by
      have : 0 < q.toNat := by omega
      exact_mod_cast this
    have h_nonneg : 0 ≤ s % (q.toNat : Int) := Int.emod_nonneg s hne
    have h_lt : s % (q.toNat : Int) < (q.toNat : Int) := Int.emod_lt_of_pos s hpos
    rw [if_neg (by linarith)]
    have h_toNat_lt : (s % (q.toNat : Int)).toNat < q.toNat := by omega
    have h_bound : (s % (q.toNat : Int)).toNat < UInt64.size :=
      lt_of_lt_of_le h_toNat_lt hq_size
    calc ((s % (q.toNat : Int)).toNat.toUInt64).toNat
        = (s % (q.toNat : Int)).toNat := by simp [h_bound]
      _ < q.toNat := h_toNat_lt

/-- 非零 Zp 元素的逆仍然非零：a.val ≠ 0 → (a.inv).val ≠ 0。 -/
private lemma Zp_inv_val_nonzero (a : Zp) (hred_a : Zp.Reduced p a) (hval_nonzero : a.val.toNat ≠ 0)
    (hp_size : 2 * p ≤ UInt64.size) : (a.inv).val.toNat ≠ 0 := by
  have ha_prime : a.prime.toNat = p := hred_a.1
  have ha_inv_prime : (a.inv).prime.toNat = p := by
    simpa [Zp.inv] using ha_prime
  have hp_lt_U64 : p < UInt64.size := by
    nlinarith
  have h_mul_one : Zp.toZMod p (a * a.inv) = (1 : ZMod p) :=
    Zp_toZMod_inv_mul_self a hred_a hval_nonzero hp_lt_U64
  have h_mul_split : Zp.toZMod p (a * a.inv) = Zp.toZMod p a * Zp.toZMod p (a.inv) :=
    Zp.toZMod_mul_weak a a.inv ha_prime ha_inv_prime hp_size
  rw [h_mul_split] at h_mul_one
  intro hzero
  have hzero_toZMod : Zp.toZMod p (a.inv) = 0 := by
    simp [Zp.toZMod, hzero]
  rw [hzero_toZMod] at h_mul_one
  simp at h_mul_one

/-- Zp 乘法保持 val ≠ 0（ZMod p 是整环）。 -/
private lemma Zp_mul_val_nonzero {a b : Zp} (hred_a : Zp.Reduced p a) (hred_b : Zp.Reduced p b)
    (ha : a.val.toNat ≠ 0) (hb : b.val.toNat ≠ 0) (hp_size : 2 * p ≤ UInt64.size) :
    (a * b).val.toNat ≠ 0 := by
  have ha_prime : a.prime.toNat = p := hred_a.1
  have ha_lt : a.val.toNat < p := hred_a.2
  have hb_lt : b.val.toNat < p := hred_b.2
  have hp_prime : Nat.Prime p := hp.out
  have hp_lt_U64 : p < UInt64.size := by nlinarith
  have h_val : (a * b).val.toNat = (a.val.toNat * b.val.toNat) % p := by
    have h_mod_lt : (a.val.toNat * b.val.toNat) % p < UInt64.size :=
      calc
        (a.val.toNat * b.val.toNat) % p < p := Nat.mod_lt (a.val.toNat * b.val.toNat) (Nat.Prime.pos hp_prime)
        _ < UInt64.size := hp_lt_U64
    have h_mul_val : (a * b).val = ((a.val.toNat * b.val.toNat) % a.prime.toNat).toUInt64 := rfl
    calc
      (a * b).val.toNat = (((a.val.toNat * b.val.toNat) % a.prime.toNat).toUInt64).toNat := by rw [h_mul_val]
      _ = (a.val.toNat * b.val.toNat) % a.prime.toNat := by
        -- .toUInt64.toNat = id for values < UInt64.size
        have h_mod_lt' : (a.val.toNat * b.val.toNat) % a.prime.toNat < UInt64.size := by
          rw [ha_prime]; exact h_mod_lt
        simp [h_mod_lt']
      _ = (a.val.toNat * b.val.toNat) % p := by rw [ha_prime]
  intro hzero
  rw [h_val] at hzero
  have h_dvd : p ∣ a.val.toNat * b.val.toNat := Nat.dvd_of_mod_eq_zero hzero
  rcases hp_prime.dvd_mul.mp h_dvd with (hp_a | hp_b)
  · exact ha (Nat.eq_zero_of_dvd_of_lt hp_a ha_lt)
  · exact hb (Nat.eq_zero_of_dvd_of_lt hp_b hb_lt)

private lemma UInt64_mul_mod_toZMod_eq (a b m : UInt64) (p' : ℕ)
    (hm : m.toNat = p') (h_no_ov : a.toNat * b.toNat < 2 ^ 64) :
    ((a * b % m).toNat : ZMod p') = (a.toNat : ZMod p') * (b.toNat : ZMod p') := by
  have h1 : ((a * b % m).toNat : ZMod p') = ((a * b).toNat : ZMod p') := by
    calc
      ((a * b % m).toNat : ZMod p') = (((a * b).toNat % m.toNat : ℕ) : ZMod p') := by
        simp [UInt64.toNat_mod]
      _ = (((a * b).toNat % p' : ℕ) : ZMod p') := by rw [hm]
      _ = ((a * b).toNat : ZMod p') := by
        simpa using (ZMod.natCast_mod ((a * b).toNat) p')
  have h2 : ((a * b).toNat : ZMod p') = (a.toNat : ZMod p') * (b.toNat : ZMod p') := by
    have h_mul := UInt64.toNat_mul a b
    have h_mul' : ((a * b).toNat : ZMod p') = ((a.toNat * b.toNat % 2 ^ 64 : ℕ) : ZMod p') := by
      simpa using congrArg (fun (n : Nat) => (n : ZMod p')) h_mul
    calc
      ((a * b).toNat : ZMod p') = ((a.toNat * b.toNat % 2 ^ 64 : ℕ) : ZMod p') := h_mul'
      _ = ((a.toNat * b.toNat : ℕ) : ZMod p') := by
        rw [Nat.mod_eq_of_lt h_no_ov]
      _ = (a.toNat : ZMod p') * (b.toNat : ZMod p') := by simp
  calc
    ((a * b % m).toNat : ZMod p') = ((a * b).toNat : ZMod p') := h1
    _ = (a.toNat : ZMod p') * (b.toNat : ZMod p') := h2

-- 辅助引理：listSum 层面的求导对应（不涉及 Array，只在 List 上操作）
private lemma listSum_derivative_eq (p : ℕ) (xs : List (UMonomial × Zp))
    (hprime : ∀ x ∈ xs, x.2.prime.toNat = p)
    (h_no_overflow : ∀ x ∈ xs, x.2.val.toNat * x.1.deg < 2 ^ 64)
    (h_deg_bound : ∀ x ∈ xs, x.1.deg < 2 ^ 64) :
    listSum p (xs.filterMap (fun (x : UMonomial × Zp) =>
      let m := x.1; let c := x.2
      if m.deg = 0 then none
      else
        let new_val := c.val * m.deg.toUInt64 % c.prime
        if new_val = 0 then none
        else some (⟨m.deg - 1⟩, ⟨new_val, c.prime⟩)
    )) = Polynomial.derivative (listSum p xs) := by
  let F : UMonomial × Zp → Option (UMonomial × Zp) :=
    fun (x : UMonomial × Zp) =>
      let m := x.1; let c := x.2
      if m.deg = 0 then none
      else
        if c.val * m.deg.toUInt64 % c.prime = 0 then none
        else some (⟨m.deg - 1⟩, ⟨c.val * m.deg.toUInt64 % c.prime, c.prime⟩)
  induction xs with
  | nil => simp [listSum, listSum_nil, F, Polynomial.derivative_zero]
  | cons x xs ih =>
    rcases x with ⟨m, c⟩
    have hprime_xs : ∀ x ∈ xs, x.2.prime.toNat = p :=
      fun y hy => hprime y (List.mem_cons_of_mem _ hy)
    have h_no_overflow_xs : ∀ x ∈ xs, x.2.val.toNat * x.1.deg < 2 ^ 64 :=
      fun y hy => h_no_overflow y (List.mem_cons_of_mem _ hy)
    have h_deg_bound_xs : ∀ x ∈ xs, x.1.deg < 2 ^ 64 :=
      fun y hy => h_deg_bound y (List.mem_cons_of_mem _ hy)
    have hp_eq : c.prime.toNat = p := hprime (m, c) (by simp)
    have h_no_overflow_term : c.val.toNat * m.deg < 2 ^ 64 := h_no_overflow (m, c) (by simp)
    have h_deg_bound_term : m.deg < 2 ^ 64 := h_deg_bound (m, c) (by simp)
    have h_degU64_toNat : (m.deg.toUInt64).toNat = m.deg := by
      have : m.deg % 2 ^ 64 = m.deg := Nat.mod_eq_of_lt h_deg_bound_term
      simpa [UInt64.toNat_ofNat] using this
    have h_no_overflow_uint64 : c.val.toNat * (m.deg.toUInt64).toNat < 2 ^ 64 := by
      rw [h_degU64_toNat]
      exact h_no_overflow_term
    have h_ih := ih hprime_xs h_no_overflow_xs h_deg_bound_xs
    by_cases hdeg : m.deg = 0
    · have h_deriv_zero' : Polynomial.derivative (monomial m.deg (Zp.toZMod p c)) = 0 := by
        simp [Polynomial.derivative_monomial, hdeg, Zp.toZMod]
      have h1 : listSum p (List.filterMap F ((m,c) :: xs)) = listSum p (List.filterMap F xs) := by
        unfold F; simp [List.filterMap, hdeg]
      have h2 : listSum p (List.filterMap F xs) = Polynomial.derivative (listSum p xs) := h_ih
      rw [h1, h2, listSum_cons, Polynomial.derivative_add, h_deriv_zero', zero_add]
    · let new_val : UInt64 := c.val * m.deg.toUInt64 % c.prime
      by_cases hnew : new_val = 0
      · dsimp [new_val] at hnew
        have hzero' : (Zp.toZMod p c) * (m.deg : ZMod p) = 0 := by
          calc
            (Zp.toZMod p c) * (m.deg : ZMod p) = (c.val.toNat : ZMod p) * (m.deg : ZMod p) := by
              simp [Zp.toZMod]
            _ = (c.val.toNat : ZMod p) * ((m.deg.toUInt64).toNat : ZMod p) := by
              simp [h_degU64_toNat]
            _ = ((c.val * m.deg.toUInt64 % c.prime).toNat : ZMod p) :=
              (UInt64_mul_mod_toZMod_eq c.val m.deg.toUInt64 c.prime p hp_eq h_no_overflow_uint64).symm
            _ = 0 := by simp [hnew]
        have h1 : listSum p (List.filterMap F ((m,c) :: xs)) = listSum p (List.filterMap F xs) := by
          unfold F; simp [List.filterMap, hdeg, hnew]
        have h2 : listSum p (List.filterMap F xs) = Polynomial.derivative (listSum p xs) := h_ih
        rw [h1, h2, listSum_cons, Polynomial.derivative_add, Polynomial.derivative_monomial, hzero', monomial_zero_right, zero_add]
      · dsimp [new_val] at hnew
        have hcoeff' : ((new_val : UInt64).toNat : ZMod p) =
            (Zp.toZMod p c) * (m.deg : ZMod p) := by
          calc
            ((new_val : UInt64).toNat : ZMod p) =
                ((c.val * m.deg.toUInt64 % c.prime).toNat : ZMod p) := rfl
            _ = (c.val.toNat : ZMod p) * ((m.deg.toUInt64).toNat : ZMod p) :=
              UInt64_mul_mod_toZMod_eq c.val m.deg.toUInt64 c.prime p hp_eq h_no_overflow_uint64
            _ = (c.val.toNat : ZMod p) * (m.deg : ZMod p) := by simp [h_degU64_toNat]
            _ = (Zp.toZMod p c) * (m.deg : ZMod p) := by simp [Zp.toZMod]
        have h1 : listSum p (List.filterMap F ((m,c) :: xs)) =
            monomial (m.deg - 1) ((new_val.toNat : ZMod p)) + listSum p (List.filterMap F xs) := by
          calc
            listSum p (List.filterMap F ((m,c) :: xs))
                = listSum p (((⟨m.deg - 1⟩, ⟨c.val * m.deg.toUInt64 % c.prime, c.prime⟩) : UMonomial × Zp) :: List.filterMap F xs) := by
              simp [List.filterMap, F, hdeg, hnew]
            _ = monomial (m.deg - 1) (Zp.toZMod p ⟨c.val * m.deg.toUInt64 % c.prime, c.prime⟩) + listSum p (List.filterMap F xs) := by
              rfl
            _ = monomial (m.deg - 1) (((c.val * m.deg.toUInt64 % c.prime).toNat : ZMod p)) + listSum p (List.filterMap F xs) := by
              simp [Zp.toZMod]
            _ = monomial (m.deg - 1) ((new_val.toNat : ZMod p)) + listSum p (List.filterMap F xs) := rfl
        have h2 : monomial (m.deg - 1) ((new_val.toNat : ZMod p)) + listSum p (List.filterMap F xs) =
            monomial (m.deg - 1) ((Zp.toZMod p c) * (m.deg : ZMod p)) + Polynomial.derivative (listSum p xs) := by
          rw [hcoeff', h_ih]
        rw [h1, h2, listSum_cons, Polynomial.derivative_add, Polynomial.derivative_monomial]

/-- derivative 的 toPoly 对应：SparsePolyZp 的求导 ↔ Polynomial 的求导 -/
private lemma derivative_toPoly_eq (f : SparsePolyZp) (hwf : SparsePolyZp.WellFormed p f)
    (hp_size : 2 * p ≤ UInt64.size)
    (h_no_overflow : ∀ x ∈ f.toList, x.2.val.toNat * x.1.deg < 2 ^ 64)
    (h_deg_bound : ∀ x ∈ f.toList, x.1.deg < 2 ^ 64)
    : SparsePolyZp.toPoly p (SparsePolyZp.derivative f) = Polynomial.derivative (SparsePolyZp.toPoly p f) := by
  -- UInt64 乘法取模到 ZMod 的转换
  have UInt64_mul_mod_toZMod_eq (a b m : UInt64) (p' : ℕ)
    (hm : m.toNat = p') (h_no_ov : a.toNat * b.toNat < 2 ^ 64) :
    ((a * b % m).toNat : ZMod p') = (a.toNat : ZMod p') * (b.toNat : ZMod p') := by
    have h1 : ((a * b % m).toNat : ZMod p') = ((a * b).toNat : ZMod p') := by
      calc
        ((a * b % m).toNat : ZMod p') = (((a * b).toNat % m.toNat : ℕ) : ZMod p') := by
          simp [UInt64.toNat_mod]
        _ = (((a * b).toNat % p' : ℕ) : ZMod p') := by rw [hm]
        _ = ((a * b).toNat : ZMod p') := by
          simpa using (ZMod.natCast_mod ((a * b).toNat) p')
    have h2 : ((a * b).toNat : ZMod p') = (a.toNat : ZMod p') * (b.toNat : ZMod p') := by
      have h_mul := UInt64.toNat_mul a b
      have h_mul' : ((a * b).toNat : ZMod p') = ((a.toNat * b.toNat % 2 ^ 64 : ℕ) : ZMod p') := by
        simpa using congrArg (fun (n : Nat) => (n : ZMod p')) h_mul
      calc
        ((a * b).toNat : ZMod p') = ((a.toNat * b.toNat % 2 ^ 64 : ℕ) : ZMod p') := h_mul'
        _ = ((a.toNat * b.toNat : ℕ) : ZMod p') := by
          rw [Nat.mod_eq_of_lt h_no_ov]
        _ = (a.toNat : ZMod p') * (b.toNat : ZMod p') := by simp
    calc
      ((a * b % m).toNat : ZMod p') = ((a * b).toNat : ZMod p') := h1
      _ = (a.toNat : ZMod p') * (b.toNat : ZMod p') := h2
  unfold SparsePolyZp.toPoly SparsePolyZp.derivative
  by_cases h_empty : f.isEmpty
  · have h_empty_list : f.toList = [] := by
      have : f.size = 0 := (Array.isEmpty_iff_size_eq_zero.mp h_empty)
      simpa [Array.length_toList]
    simp [h_empty, h_empty_list, listSum, listSum_nil, Polynomial.derivative_zero]
  · have hprime0 : (f[0]!.snd.prime).toNat = p := by
      have hpos : 0 < f.size := by
        rw [Array.isEmpty_iff_size_eq_zero] at h_empty
        omega
      have hmem : f[0] ∈ f := Array.getElem_mem (i := 0) (h := hpos) (xs := f)
      have hmem' : f[0]! ∈ f := by simpa [hpos] using hmem
      exact hwf (f[0]!) hmem'
    have hprime_all : ∀ x ∈ f.toList, x.2.prime.toNat = p := by
      intro x hx
      have hx' : x ∈ f := (Array.mem_def.mpr hx)
      exact hwf x hx'
    have hprime_eq : ∀ x ∈ f.toList, x.2.prime = f[0]!.snd.prime := by
      intro x hx
      have hcalc : x.2.prime.toNat = (f[0]!.snd.prime).toNat := by
        calc
          x.2.prime.toNat = p := hprime_all x hx
          _ = (f[0]!.snd.prime).toNat := by symm; exact hprime0
      exact ((UInt64.toNat_inj (a := x.2.prime) (b := f[0]!.snd.prime)).mp hcalc)
    have h_filter_eq : List.filterMap (fun (x : UMonomial × Zp) =>
        let m := x.1; let c := x.2
        if m.deg = 0 then none
        else
          let new_val := c.val * m.deg.toUInt64 % f[0]!.snd.prime
          if new_val = 0 then none
          else some ((⟨m.deg - 1⟩ : UMonomial), (⟨new_val, f[0]!.snd.prime⟩ : Zp))
      ) f.toList = List.filterMap (fun (x : UMonomial × Zp) =>
        let m := x.1; let c := x.2
        if m.deg = 0 then none
        else
          let new_val := c.val * m.deg.toUInt64 % c.prime
          if new_val = 0 then none
          else some ((⟨m.deg - 1⟩ : UMonomial), (⟨new_val, c.prime⟩ : Zp))
      ) f.toList := by
      refine List.filterMap_congr ?_
      intro x hx
      rcases x with ⟨m, c⟩
      have h_eq_prime : c.prime = f[0]!.snd.prime := by
        simpa using hprime_eq (m,c) hx
      simp [h_eq_prime]
    have h_target : listSum p (List.filterMap (fun (x : UMonomial × Zp) =>
        let m := x.1; let c := x.2
        if m.deg = 0 then none
        else
          let new_val := c.val * m.deg.toUInt64 % f[0]!.snd.prime
          if new_val = 0 then none
          else some (⟨m.deg - 1⟩, ⟨new_val, f[0]!.snd.prime⟩)
      ) f.toList) = Polynomial.derivative (listSum p f.toList) := by
      rw [h_filter_eq]
      exact listSum_derivative_eq p f.toList hprime_all h_no_overflow h_deg_bound
    simpa [Array.toList_filterMap, h_empty]

/-- Helper lemma: in a non-empty SparsePolyZp, the first element (via r[0]!) is in its toList. -/
private lemma mem_getFirst_toList (r : SparsePolyZp) (hpos : 0 < r.size) : (r[0]! ∈ r.toList) := by
  have hlen : r.toList.length = r.size := by simp
  have hpos' : 0 < r.toList.length := by rw [hlen]; exact hpos
  have hmem_list : r.toList.get ⟨0, hpos'⟩ ∈ r.toList := List.get_mem _ _
  have h_eq : r.toList.get ⟨0, hpos'⟩ = r[0] := by
    simp
  have h_eq' : r[0]! = r[0] := by
    rw [getElem!_def, getElem?_def]
    simp [hpos]
  have hmem : r[0] ∈ r.toList := by
    rw [← h_eq]
    exact hmem_list
  rw [h_eq']
  exact hmem

lemma drop_eq_get_cons_general {α : Type} (l : List α) (n : ℕ) (h : n < l.length) :
    l.drop n = (l.get ⟨n, h⟩) :: l.drop (n + 1) := by
  revert l
  induction n with
  | zero =>
    intro l h
    cases l
    · simp at h
    · simp
  | succ n ih =>
    intro l h
    cases l with
    | nil => simp at h
    | cons hd tl =>
      have h' : n < tl.length := by
        simpa [List.length_cons] using h
      have h_ih := ih tl h'
      simpa using h_ih

/-- drop n xs = xs[n] :: drop (n+1) xs when n < length xs。 -/
lemma drop_eq_get_cons' (l : List (UMonomial × Zp)) (n : ℕ) (h : n < l.length) :
    l.drop n = (l.get ⟨n, h⟩) :: l.drop (n + 1) := by
  revert l
  induction n with
  | zero =>
    intro l h
    cases l with
    | nil => simp at h
    | cons hd tl => simp
  | succ n ih =>
    intro l h
    cases l with
    | nil => simp at h
    | cons hd tl =>
      have h' : n < tl.length := by
        simpa [List.length_cons] using h
      have h_ih := ih tl h'
      simpa using h_ih

/-- toPoly 在 Array.push 下的表现。 -/
lemma toPoly_push (g : SparsePolyZp) (m : UMonomial) (c : Zp) :
    SparsePolyZp.toPoly p (g.push (m, c)) = SparsePolyZp.toPoly p g + Polynomial.monomial m.deg (Zp.toZMod p c) := by
  unfold SparsePolyZp.toPoly
  simp [Array.push, listSum_append, listSum, add_assoc]

/-- listSum of a single element。 -/
lemma listSum_singleton (m : UMonomial) (c : Zp) : listSum p [(m, c)] = Polynomial.monomial m.deg (Zp.toZMod p c) := by
  simp [listSum]

/-- Boolean `ReducedB` 转为数学侧逐项约化性质。 -/
private lemma allReduced_of_reducedB (f : SparsePolyZp) (pm : UInt64)
    (hpm : pm.toNat = p) (hred : SparsePolyZp.ReducedB f pm) :
    SparsePolyZp.AllReduced p f.toList := by
  have aux : ∀ xs : List (UMonomial × Zp),
      SparsePolyZp.reducedListB pm xs = true → SparsePolyZp.AllReduced p xs := by
    intro xs
    induction xs with
    | nil => simp [SparsePolyZp.AllReduced]
    | cons a rest ih =>
      intro h
      rcases (SparsePolyZp.reducedListB_cons pm a rest).mp h with ⟨ha, hrest⟩
      intro x hx
      rcases List.mem_cons.mp hx with rfl | hx
      · constructor
        · rw [ha.2, hpm]
        · rw [UInt64.lt_iff_toNat_lt] at ha
          simpa [hpm] using ha.1
      · exact ih hrest x hx
  exact aux f.toList hred

/-- SQF/GCD 所需的稀疏表示规范形。空数组也满足该谓词。 -/
private def CanonicalRep (p : Nat) (f : SparsePolyZp) : Prop :=
  SparsePolyZp.Sorted f ∧ SparsePolyZp.NonZeroB f ∧ SparsePolyZp.AllReduced p f.toList

private lemma nonzeroB_of_val_nonzero (f : SparsePolyZp)
    (h : ∀ x ∈ f.toList, x.2.val ≠ 0) : SparsePolyZp.NonZeroB f := by
  unfold SparsePolyZp.NonZeroB
  have aux : ∀ l : List (UMonomial × Zp),
      (∀ x ∈ l, x.2.val ≠ 0) → SparsePolyZp.nonzeroListB l = true := by
    intro l
    induction l with
    | nil => simp [SparsePolyZp.nonzeroListB]
    | cons a rest ih =>
        intro hl
        rw [SparsePolyZp.nonzeroListB_cons]
        exact ⟨hl a (by simp), ih (fun x hx => hl x (List.mem_cons_of_mem a hx))⟩
  exact aux f.toList h

private lemma val_nonzero_of_nonzeroB (f : SparsePolyZp)
    (h : SparsePolyZp.NonZeroB f) : ∀ x ∈ f.toList, x.2.val ≠ 0 := by
  unfold SparsePolyZp.NonZeroB at h
  have aux : ∀ l : List (UMonomial × Zp),
      SparsePolyZp.nonzeroListB l = true → ∀ x ∈ l, x.2.val ≠ 0 := by
    intro l
    induction l with
    | nil => simp
    | cons a rest ih =>
        intro hl x hx
        rw [SparsePolyZp.nonzeroListB_cons] at hl
        rcases List.mem_cons.mp hx with rfl | hx
        · exact hl.1
        · exact ih hl.2 x hx
  exact aux f.toList h

private lemma sortedListB_filterMap_pred
    (G : (UMonomial × Zp) → Option (UMonomial × Zp))
    (hG : ∀ x y, G x = some y → y.1.deg + 1 = x.1.deg) :
    ∀ l : List (UMonomial × Zp), SparsePolyZp.sortedListB l = true →
      SparsePolyZp.sortedListB (l.filterMap G) = true := by
  intro l
  induction l with
  | nil => simp [SparsePolyZp.sortedListB]
  | cons a rest ih =>
      intro hs
      cases hGa : G a with
      | none =>
          rw [List.filterMap_cons_none hGa]
          exact ih ((SparsePolyZp.sortedListB_iff a rest).mp hs).2
      | some b =>
          rw [List.filterMap_cons_some hGa, SparsePolyZp.sortedListB_iff]
          have hst := (SparsePolyZp.sortedListB_iff a rest).mp hs
          refine ⟨?_, ih hst.2⟩
          intro z hz
          rcases List.mem_filterMap.mp hz with ⟨x, hx, hxz⟩
          have hxa := hst.1 x hx
          have hba := hG a b hGa
          have hzx := hG x z hxz
          omega

/-- 稀疏求导在非空规范输入上保持 GCD 所需的规范形。 -/
private lemma derivative_canonical (f : SparsePolyZp) (hf : CanonicalRep p f)
    (hfne : ¬f.isEmpty) : CanonicalRep p (SparsePolyZp.derivative f) := by
  let D : (UMonomial × Zp) → Option (UMonomial × Zp) := fun x =>
    if x.1.deg = 0 then none
    else
      let v := x.2.val * x.1.deg.toUInt64 % f[0]!.2.prime
      if v = 0 then none else some (⟨x.1.deg - 1⟩, ⟨v, f[0]!.2.prime⟩)
  have hlist : (SparsePolyZp.derivative f).toList = f.toList.filterMap D := by
    simp [SparsePolyZp.derivative, hfne, D, Array.toList_filterMap]
  have hfrontmem : f[0]! ∈ f.toList := mem_getFirst_toList f (by
    rw [Array.isEmpty_iff_size_eq_zero] at hfne
    omega)
  have hprime : f[0]!.2.prime.toNat = p := (hf.2.2 f[0]! hfrontmem).1
  have hsorted : SparsePolyZp.Sorted (SparsePolyZp.derivative f) := by
    unfold SparsePolyZp.Sorted
    rw [hlist]
    apply sortedListB_filterMap_pred D
    · intro x y hxy
      by_cases hd : x.1.deg = 0
      · simp [D, hd] at hxy
      let v := x.2.val * x.1.deg.toUInt64 % f[0]!.2.prime
      by_cases hv : v = 0
      · simp [D, hd, v, hv] at hxy
      simp [D, hd, v, hv] at hxy
      have hydeg : y.1.deg = x.1.deg - 1 := by
        simpa using congrArg (fun z : UMonomial × Zp => z.1.deg) hxy.symm
      rw [hydeg]
      omega
    · simpa [SparsePolyZp.Sorted] using hf.1
  have hnonzero : SparsePolyZp.NonZeroB (SparsePolyZp.derivative f) := by
    apply nonzeroB_of_val_nonzero
    intro z hz
    rw [hlist] at hz
    rcases List.mem_filterMap.mp hz with ⟨x, hx, hxz⟩
    simp only [D] at hxz
    split at hxz <;> try contradiction
    split at hxz <;> try contradiction
    injection hxz with hxz
    simpa [← hxz] using ‹_›
  have hreduced : SparsePolyZp.AllReduced p (SparsePolyZp.derivative f).toList := by
    intro z hz
    rw [hlist] at hz
    rcases List.mem_filterMap.mp hz with ⟨x, hx, hxz⟩
    simp only [D] at hxz
    split at hxz <;> try contradiction
    split at hxz <;> try contradiction
    injection hxz with hxz
    rw [← hxz]
    constructor
    · exact hprime
    · rw [← hprime]
      rw [← UInt64.lt_iff_toNat_lt]
      exact UInt64.mod_lt (x.2.val * x.1.deg.toUInt64) (by
        rw [UInt64.lt_iff_toNat_lt, hprime]
        exact hp.out.pos)
  exact ⟨hsorted, hnonzero, hreduced⟩

private lemma derivative_deg_bound_of (B : Nat) (f : SparsePolyZp)
    (hdeg : ∀ x ∈ f.toList, x.1.deg < B) :
    ∀ x ∈ (SparsePolyZp.derivative f).toList, x.1.deg < B := by
  intro y hy
  unfold SparsePolyZp.derivative at hy
  split at hy
  · simp at hy
  · rw [Array.toList_filterMap] at hy
    rcases List.mem_filterMap.mp hy with ⟨x, hx, hxy⟩
    simp only at hxy
    split at hxy <;> try contradiction
    split at hxy <;> try contradiction
    injection hxy with hxy
    have hydeg : y.1.deg = x.1.deg - 1 := by
      simpa using congrArg (fun z : UMonomial × Zp => z.1.deg) hxy.symm
    rw [hydeg]
    exact lt_of_le_of_lt (Nat.sub_le _ _) (hdeg x hx)

/-- 若稀疏求导在 `ZMod p` 上为空，则每个存储项的次数都被 `p` 整除。 -/
private lemma degree_dvd_of_derivative_empty (f : SparsePolyZp)
    (hf : CanonicalRep p f) (hfne : ¬f.isEmpty)
    (hderiv : SparsePolyZp.derivative f = SparsePolyZp.empty)
    (h_no_overflow : ∀ x ∈ f.toList, x.2.val.toNat * x.1.deg < 2 ^ 64)
    (h_deg_bound : ∀ x ∈ f.toList, x.1.deg < 2 ^ 64) :
    ∀ x ∈ f.toList, p ∣ x.1.deg := by
  let D : (UMonomial × Zp) → Option (UMonomial × Zp) := fun x =>
    if x.1.deg = 0 then none
    else
      let v := x.2.val * x.1.deg.toUInt64 % f[0]!.2.prime
      if v = 0 then none else some (⟨x.1.deg - 1⟩, ⟨v, f[0]!.2.prime⟩)
  have hfilter : f.toList.filterMap D = [] := by
    have h := congrArg Array.toList hderiv
    simpa [SparsePolyZp.derivative, hfne, SparsePolyZp.empty, D] using h
  have hfrontmem : f[0]! ∈ f.toList := mem_getFirst_toList f (by
    rw [Array.isEmpty_iff_size_eq_zero] at hfne
    omega)
  have hprime : f[0]!.2.prime.toNat = p := (hf.2.2 f[0]! hfrontmem).1
  intro x hx
  by_cases hd0 : x.1.deg = 0
  · simp [hd0]
  have hv0 : x.2.val * x.1.deg.toUInt64 % f[0]!.2.prime = 0 := by
    by_contra hv
    have hm : (⟨x.1.deg - 1⟩,
        ⟨x.2.val * x.1.deg.toUInt64 % f[0]!.2.prime,
          f[0]!.2.prime⟩) ∈ f.toList.filterMap D := by
      apply List.mem_filterMap.mpr
      refine ⟨x, hx, ?_⟩
      simp [D, hd0, hv]
    rw [hfilter] at hm
    simp at hm
  have hdeg_round : x.1.deg.toUInt64.toNat = x.1.deg := by
    have hm := Nat.mod_eq_of_lt (h_deg_bound x hx)
    simpa [UInt64.toNat_ofNat] using hm
  have hmul_round : (x.2.val * x.1.deg.toUInt64).toNat =
      x.2.val.toNat * x.1.deg := by
    rw [UInt64.toNat_mul, hdeg_round,
      Nat.mod_eq_of_lt (h_no_overflow x hx)]
  have hmod0 : (x.2.val.toNat * x.1.deg) % p = 0 := by
    have h := congrArg UInt64.toNat hv0
    simpa [UInt64.toNat_mod, hmul_round, hprime] using h
  have hdvd_mul : p ∣ x.2.val.toNat * x.1.deg := Nat.dvd_of_mod_eq_zero hmod0
  rcases hp.out.dvd_mul.mp hdvd_mul with hval | hdeg
  · have hval_lt := (hf.2.2 x hx).2
    have hval0 := Nat.eq_zero_of_dvd_of_lt hval hval_lt
    exact False.elim ((val_nonzero_of_nonzeroB f hf.2.1 x hx)
      (UInt64.toNat_inj.mp (by simpa using hval0)))
  · exact hdeg

private lemma isChain_of_sortedListB : ∀ xs : List (UMonomial × Zp),
    SparsePolyZp.sortedListB xs = true →
      List.IsChain (fun a b : UMonomial × Zp => a.fst.deg > b.fst.deg) xs := by
  intro xs
  induction xs with
  | nil => simp
  | cons a rest ih =>
      intro hs
      rw [SparsePolyZp.sortedListB_iff] at hs
      cases rest with
      | nil => simp
      | cons b rest =>
          rw [List.isChain_cons]
          refine ⟨?_, ih hs.2⟩
          intro y hy
          simp at hy
          subst y
          exact hs.1 b (by simp)

private lemma listSum_natDegree_le_of_all_le (xs : List (UMonomial × Zp)) (d : Nat)
    (h : ∀ x ∈ xs, x.1.deg ≤ d) : (listSum p xs).natDegree ≤ d := by
  induction xs with
  | nil => simp
  | cons a rest ih =>
      rcases a with ⟨m, c⟩
      rw [listSum_cons]
      refine (Polynomial.natDegree_add_le _ _).trans ?_
      simp only [max_le_iff]
      exact ⟨(Polynomial.natDegree_monomial_le (R := ZMod p) c).trans (h (m, c) (by simp)),
        ih (fun x hx => h x (List.mem_cons_of_mem _ hx))⟩

/-- 规范稀疏表示的数组首项确实给出数学多项式的次数和首项系数。 -/
private lemma toPoly_head_data (f : SparsePolyZp) (hf : CanonicalRep p f)
    (hne : ¬f.isEmpty) :
    (SparsePolyZp.toPoly p f).natDegree = f[0]!.fst.deg ∧
      (SparsePolyZp.toPoly p f).coeff f[0]!.fst.deg = Zp.toZMod p f[0]!.snd := by
  have hpos : 0 < f.size := by
    rw [Array.isEmpty_iff_size_eq_zero] at hne
    omega
  have hlist := SparsePolyZp.toList_cons_of_ne_empty f hne
  have hchain : List.IsChain (fun a b : UMonomial × Zp => a.fst.deg > b.fst.deg)
      f.toList := isChain_of_sortedListB f.toList hf.1
  rw [hlist] at hchain
  have hcoeff : (SparsePolyZp.toPoly p f).coeff f[0]!.fst.deg =
      Zp.toZMod p f[0]!.snd := by
    unfold SparsePolyZp.toPoly
    rw [hlist]
    exact listSum_coeff_at_head p f[0]! _ hchain
  have hmem : f[0]! ∈ f.toList := mem_getFirst_toList f hpos
  have hheadnz : f[0]!.snd.val ≠ 0 := by
    have hnz := hf.2.1
    unfold SparsePolyZp.NonZeroB at hnz
    rw [hlist, SparsePolyZp.nonzeroListB_cons] at hnz
    exact hnz.1
  have hcoef_ne : Zp.toZMod p f[0]!.snd ≠ 0 :=
    Zp.toZMod_ne_zero_of_val_ne_zero p f[0]!.snd (hf.2.2 f[0]! hmem)
      hheadnz
  have hlower : f[0]!.fst.deg ≤ (SparsePolyZp.toPoly p f).natDegree :=
    Polynomial.le_natDegree_of_ne_zero (hcoeff.trans_ne hcoef_ne)
  have hall : ∀ x ∈ f.toList, x.1.deg ≤ f[0]!.fst.deg := by
    intro x hx
    rcases List.mem_cons.mp (hlist ▸ hx) with hx0 | hxr
    · simpa [hx0]
    · exact Nat.le_of_lt (chain_gt_all_after_head f[0]! _ hchain x hxr)
  have hupper := listSum_natDegree_le_of_all_le (p := p) f.toList f[0]!.fst.deg hall
  exact ⟨Nat.le_antisymm hupper hlower, hcoeff⟩

private lemma reducedB_of_allReduced (f : SparsePolyZp) (pm : UInt64)
    (hpm : pm.toNat = p) (hred : SparsePolyZp.AllReduced p f.toList) :
    SparsePolyZp.ReducedB f pm := by
  have aux : ∀ xs : List (UMonomial × Zp), SparsePolyZp.AllReduced p xs →
      SparsePolyZp.reducedListB pm xs = true := by
    intro xs
    induction xs with
    | nil => intro _; rfl
    | cons a rest ih =>
      intro hs
      rw [SparsePolyZp.reducedListB_cons]
      have ha := hs a (by simp)
      have hrest : SparsePolyZp.AllReduced p rest :=
        fun x hx => hs x (List.mem_cons_of_mem _ hx)
      constructor
      · constructor
        · rw [UInt64.lt_iff_toNat_lt, hpm]
          exact ha.2
        · apply UInt64.toNat_inj.mp
          rw [ha.1, hpm]
      · exact ih hrest
  exact aux f.toList hred

private lemma head_val_nonzero_of_nonZeroB (f : SparsePolyZp) (hne : ¬f.isEmpty)
    (hnz : SparsePolyZp.NonZeroB f) : f[0]!.snd.val.toNat ≠ 0 := by
  have hlist : f.toList = f[0]! :: f.toList.tail := SparsePolyZp.toList_cons_of_ne_empty f hne
  have hhead := ((SparsePolyZp.nonzeroListB_cons f[0]! f.toList.tail).mp (hlist ▸ hnz)).1
  intro hz
  apply hhead
  exact UInt64.toNat_inj.mp (by simpa using hz)

private lemma Zp_inv_reduced (a : Zp) (ha : Zp.Reduced p a)
    (hp_lt : p < UInt64.size) : Zp.Reduced p a.inv := by
  have hp_one : 1 < p := hp.out.one_lt
  refine ⟨by simpa [Zp.inv] using ha.1, ?_⟩
  calc
    a.inv.val.toNat = (Zp.modInv a.val a.prime).toNat := rfl
    _ < a.prime.toNat := modInv_val_lt a.val a.prime
      (by simpa [ha.1] using hp_one) (by simpa [ha.1] using hp_lt.le)
    _ = p := ha.1

private lemma sortedListB_append_singleton (xs : List (UMonomial × Zp))
    (z : UMonomial × Zp) (hs : SparsePolyZp.sortedListB xs = true)
    (hz : ∀ x ∈ xs, z.1.deg < x.1.deg) :
    SparsePolyZp.sortedListB (xs ++ [z]) = true := by
  induction xs with
  | nil => simp [SparsePolyZp.sortedListB]
  | cons a rest ih =>
      change SparsePolyZp.sortedListB (a :: (rest ++ [z])) = true
      rw [SparsePolyZp.sortedListB_iff] at hs ⊢
      constructor
      · intro x hx
        rcases List.mem_append.mp hx with hx | hx
        · exact hs.1 x hx
        · simp at hx
          subst x
          exact hz a (by simp)
      · exact ih hs.2 (fun x hx => hz x (List.mem_cons_of_mem _ hx))

private lemma nonzeroListB_append_singleton (xs : List (UMonomial × Zp))
    (z : UMonomial × Zp) (hs : SparsePolyZp.nonzeroListB xs = true)
    (hz : z.2.val ≠ 0) : SparsePolyZp.nonzeroListB (xs ++ [z]) = true := by
  induction xs with
  | nil => simpa [SparsePolyZp.nonzeroListB] using hz
  | cons a rest ih =>
      change SparsePolyZp.nonzeroListB (a :: (rest ++ [z])) = true
      rw [SparsePolyZp.nonzeroListB_cons] at hs ⊢
      exact ⟨hs.1, ih hs.2⟩

private lemma val_ne_zero_of_nonzeroListB (xs : List (UMonomial × Zp))
    (hs : SparsePolyZp.nonzeroListB xs = true) : ∀ x ∈ xs, x.2.val ≠ 0 := by
  induction xs with
  | nil => simp
  | cons a rest ih =>
      rw [SparsePolyZp.nonzeroListB_cons] at hs
      intro x hx
      rcases List.mem_cons.mp hx with rfl | hx
      · exact hs.1
      · exact ih hs.2 x hx

private lemma nonzeroListB_filterMap_of_output
    (G : (UMonomial × Zp) → Option (UMonomial × Zp))
    (hG : ∀ x y, G x = some y → y.2.val ≠ 0) :
    ∀ xs : List (UMonomial × Zp),
      SparsePolyZp.nonzeroListB (xs.filterMap G) = true := by
  intro xs
  induction xs with
  | nil => rfl
  | cons a rest ih =>
      cases hGa : G a with
      | none => simpa [List.filterMap_cons_none hGa] using ih
      | some b =>
          rw [List.filterMap_cons_some hGa, SparsePolyZp.nonzeroListB_cons]
          exact ⟨hG a b hGa, ih⟩

/-- 两个规范形多项式构成 `divmod` 的合法分支输入。 -/
private lemma divmod_valid_of_canonical (f g : SparsePolyZp) (hg_ne : ¬g.isEmpty)
    (hf : CanonicalRep p f) (hg : CanonicalRep p g) (hp_lt : p < UInt64.size) :
    ¬ g.isEmpty ∧ ((g[0]!.snd.inv * (g[0]!).snd).val = 1) ∧
      SparsePolyZp.Sorted g ∧ SparsePolyZp.Sorted f ∧
      0 < (g[0]!.snd.prime).toNat ∧ SparsePolyZp.ReducedB g (g[0]!.snd.prime) ∧
      SparsePolyZp.ReducedB f (g[0]!.snd.prime) ∧ SparsePolyZp.NonZeroB f := by
  have hmem : g[0]! ∈ g.toList := by
    have hpos : 0 < g.size := by
      rw [Array.isEmpty_iff_size_eq_zero] at hg_ne
      omega
    exact mem_getFirst_toList g hpos
  have hgred : Zp.Reduced p g[0]!.snd := hg.2.2 g[0]! hmem
  have hpm : g[0]!.snd.prime.toNat = p := hgred.1
  have hval_ne : g[0]!.snd.val.toNat ≠ 0 :=
    head_val_nonzero_of_nonZeroB g hg_ne hg.2.1
  exact ⟨hg_ne, Zp_inv_mul_val_eq_one g[0]!.snd hgred hval_ne hp_lt,
    hg.1, hf.1, by rw [hpm]; exact hp.out.pos,
    reducedB_of_allReduced g _ hpm hg.2.2,
    reducedB_of_allReduced f _ hpm hf.2.2, hf.2.1⟩

/-- 实际 `divmodAux` 保持多项式长除法恒等式。

证明直接使用实现定义生成的良基归纳原理；循环不变量是
`F = q * g + r`。递归步把首项消去项同时加入商、从余式中减去，
再由 `toPoly_sub`/`toPoly_mul` 化为环恒等式。 -/
private theorem divmodAux_toPoly_identity (g : SparsePolyZp) (dg : Nat)
    (lc_g_inv : Zp) (pm : UInt64) (hq : 0 < pm.toNat)
    (hg_ne : ¬ g.isEmpty) (hg_red : SparsePolyZp.ReducedB g pm)
    (h_dg : (g[0]!).fst.deg = dg) (hlp : lc_g_inv.prime = pm)
    (h_lc : (lc_g_inv * (g[0]!).snd).val = 1) (h_sorted_g : SparsePolyZp.Sorted g)
    (hpm : pm.toNat = p) (h2p : 2 * p ≤ UInt64.size) (hp2 : p * p ≤ UInt64.size) :
    ∀ (q r : SparsePolyZp) (h_sorted_r : SparsePolyZp.Sorted r)
      (hr_red : SparsePolyZp.ReducedB r pm) (hr_nz : SparsePolyZp.NonZeroB r),
      SparsePolyZp.AllReduced p q.toList →
      let out := SparsePolyZp.divmodAux g dg lc_g_inv pm hq hg_ne hg_red h_dg hlp h_lc
        h_sorted_g q r h_sorted_r hr_red hr_nz
      SparsePolyZp.toPoly p q * SparsePolyZp.toPoly p g + SparsePolyZp.toPoly p r =
        SparsePolyZp.toPoly p out.1 * SparsePolyZp.toPoly p g + SparsePolyZp.toPoly p out.2 := by
  let motive : ∀ (q r : SparsePolyZp), SparsePolyZp.Sorted r →
      SparsePolyZp.ReducedB r pm → SparsePolyZp.NonZeroB r → Prop :=
    fun q r h_sorted_r hr_red hr_nz =>
      SparsePolyZp.AllReduced p q.toList →
      let out := SparsePolyZp.divmodAux g dg lc_g_inv pm hq hg_ne hg_red h_dg hlp h_lc
        h_sorted_g q r h_sorted_r hr_red hr_nz
      SparsePolyZp.toPoly p q * SparsePolyZp.toPoly p g + SparsePolyZp.toPoly p r =
        SparsePolyZp.toPoly p out.1 * SparsePolyZp.toPoly p g + SparsePolyZp.toPoly p out.2
  apply SparsePolyZp.divmodAux.induct g dg lc_g_inv pm hq hg_ne hg_red h_dg hlp h_lc
    h_sorted_g motive
  · intro q r hsr hrr hrnz hr hqred
    simp [SparsePolyZp.divmodAux, hr]
  · intro q r hsr hrr hrnz hr
    dsimp only
    intro hd hqred
    simp [SparsePolyZp.divmodAux, hr, hd]
  · intro q r hsr hrr hrnz hr
    dsimp only
    intro hnotlt hsr' hcoeff_prime htgred htgnz hrr' hrnz' ih hqred
    let coeff : Zp := r[0]!.snd * lc_g_inv
    let d : Nat := r[0]!.fst.deg - dg
    let term : SparsePolyZp := #[(⟨d⟩, coeff)]
    let r' : SparsePolyZp := r - term * g
    let q' : SparsePolyZp := q.push (⟨d⟩, coeff)
    have hcoeff_red : Zp.Reduced p coeff := by
      have hrlist : r.toList = r[0]! :: r.toList.tail := SparsePolyZp.toList_cons_of_ne_empty r hr
      have hrhead := ((SparsePolyZp.reducedListB_cons pm r[0]! r.toList.tail).mp
        (hrlist ▸ hrr)).1
      have hcprime : coeff.prime = pm := by simpa [coeff] using hcoeff_prime
      have hcval : coeff.val < pm := by
        dsimp [coeff]
        have hrprime : r[0]!.snd.prime = pm := hrhead.2
        have hmul := SparsePolyZp.Zp_mul_reduced r[0]!.snd lc_g_inv (by simpa [hrprime] using hq)
        simpa [hrprime] using hmul.1
      constructor
      · rw [hcprime, hpm]
      · rw [UInt64.lt_iff_toNat_lt] at hcval
        simpa [hpm] using hcval
    have hterm_red : SparsePolyZp.AllReduced p term.toList := by
      intro x hx
      have hx' : x = (⟨d⟩, coeff) := by simpa [term] using hx
      subst x
      exact hcoeff_red
    have hg_all : SparsePolyZp.AllReduced p g.toList :=
      allReduced_of_reducedB g pm hpm hg_red
    have htgred_all : SparsePolyZp.AllReduced p (term * g).toList :=
      allReduced_of_reducedB (term * g) pm hpm (by simpa [term, coeff] using htgred)
    have hrr_all : SparsePolyZp.AllReduced p r.toList :=
      allReduced_of_reducedB r pm hpm hrr
    have hrr'_all : SparsePolyZp.AllReduced p r'.toList :=
      allReduced_of_reducedB r' pm hpm (by simpa [r', term, coeff, d] using hrr')
    have hq'_red : SparsePolyZp.AllReduced p q'.toList := by
      intro x hx
      have hx' : x ∈ q.toList ∨ x = (⟨d⟩, coeff) := by simpa [q'] using hx
      rcases hx' with hx | rfl
      · exact hqred x hx
      · exact hcoeff_red
    have hqpoly : SparsePolyZp.toPoly p q' =
        SparsePolyZp.toPoly p q + SparsePolyZp.toPoly p term := by
      rw [toPoly_push]
      simp [term, SparsePolyZp.toPoly, listSum]
    have hrpoly : SparsePolyZp.toPoly p r' =
        SparsePolyZp.toPoly p r - SparsePolyZp.toPoly p (term * g) := by
      exact SparsePolyZp.toPoly_sub p h2p r (term * g) hrr_all htgred_all
    have htgpoly : SparsePolyZp.toPoly p (term * g) =
        SparsePolyZp.toPoly p term * SparsePolyZp.toPoly p g := by
      exact SparsePolyZp.toPoly_mul p h2p hp2 term g hterm_red hg_all
    have hih := ih hq'_red
    unfold SparsePolyZp.divmodAux
    simp only [hr, hnotlt, dif_neg, ↓reduceIte]
    change SparsePolyZp.toPoly p q * SparsePolyZp.toPoly p g + SparsePolyZp.toPoly p r = _
    calc
      SparsePolyZp.toPoly p q * SparsePolyZp.toPoly p g + SparsePolyZp.toPoly p r
          = (SparsePolyZp.toPoly p q + SparsePolyZp.toPoly p term) *
              SparsePolyZp.toPoly p g +
            (SparsePolyZp.toPoly p r - SparsePolyZp.toPoly p term * SparsePolyZp.toPoly p g) := by ring
      _ = SparsePolyZp.toPoly p q' * SparsePolyZp.toPoly p g + SparsePolyZp.toPoly p r' := by
        rw [hqpoly, hrpoly, htgpoly]
      _ = _ := hih

/-- 合法规范形输入上，真实 `SparsePolyZp.divmod` 满足 `f = q*g+r`。 -/
private theorem divmod_toPoly_identity (f g : SparsePolyZp)
    (hvalid : ¬ g.isEmpty ∧ ((g[0]!.snd.inv * (g[0]!).snd).val = 1) ∧
      SparsePolyZp.Sorted g ∧ SparsePolyZp.Sorted f ∧
      0 < (g[0]!.snd.prime).toNat ∧ SparsePolyZp.ReducedB g (g[0]!.snd.prime) ∧
      SparsePolyZp.ReducedB f (g[0]!.snd.prime) ∧ SparsePolyZp.NonZeroB f)
    (hpm : (g[0]!.snd.prime).toNat = p)
    (h2p : 2 * p ≤ UInt64.size) (hp2 : p * p ≤ UInt64.size) :
    SparsePolyZp.toPoly p f =
      SparsePolyZp.toPoly p (SparsePolyZp.divmod f g).1 * SparsePolyZp.toPoly p g +
      SparsePolyZp.toPoly p (SparsePolyZp.divmod f g).2 := by
  have hqred : SparsePolyZp.AllReduced p (#[] : SparsePolyZp).toList := by
    intro x hx
    simp at hx
  have hinv := divmodAux_toPoly_identity g (g[0]!.fst.deg) (g[0]!.snd.inv)
    (g[0]!.snd.prime) hvalid.2.2.2.2.1 hvalid.1 hvalid.2.2.2.2.2.1 rfl rfl
    hvalid.2.1 hvalid.2.2.1 hpm h2p hp2 #[] f hvalid.2.2.2.1
    hvalid.2.2.2.2.2.2.1 hvalid.2.2.2.2.2.2.2 hqred
  unfold SparsePolyZp.divmod
  rw [dif_pos hvalid]
  simpa using hinv

/-- `normalization` 只删除零系数项，不改变表示的数学多项式。 -/
private lemma normalization_toPoly (f : SparsePolyZp) :
    SparsePolyZp.toPoly p (SparsePolyZp.normalization f) = SparsePolyZp.toPoly p f := by
  unfold SparsePolyZp.normalization SparsePolyZp.toPoly
  rw [Array.toList_filter]
  induction f.toList with
  | nil => simp [listSum]
  | cons a rest ih =>
    rcases a with ⟨m, c⟩
    by_cases hc : c.val = 0
    · simpa [hc, listSum, Zp.toZMod] using ih
    · simp [hc, listSum, ih]

private lemma normalization_eq_of_nonZeroB (f : SparsePolyZp)
    (hf : SparsePolyZp.NonZeroB f) : SparsePolyZp.normalization f = f := by
  apply Array.ext'
  unfold SparsePolyZp.normalization
  rw [Array.toList_filter]
  apply List.filter_eq_self.mpr
  intro x hx
  have hv := val_ne_zero_of_nonzeroListB f.toList hf x hx
  simpa [bne_iff_ne] using hv

private lemma normalization_canonical (f : SparsePolyZp) (hf : CanonicalRep p f) :
    CanonicalRep p (SparsePolyZp.normalization f) := by
  rw [normalization_eq_of_nonZeroB f hf.2.1]
  exact hf

private lemma scalarMul_canonical (c : Zp) (f : SparsePolyZp)
    (hf : CanonicalRep p f) (hcprime : c.prime.toNat = p) :
    CanonicalRep p (SparsePolyZp.scalarMul c f) := by
  let G : (UMonomial × Zp) → Option (UMonomial × Zp) := fun x =>
    let new_val := x.2 * c
    if new_val.val = 0 then none else some (x.1, new_val)
  have hlist : (SparsePolyZp.scalarMul c f).toList = f.toList.filterMap G := by
    unfold SparsePolyZp.scalarMul
    rw [Array.toList_filterMap]
  have hdeg : ∀ x y, G x = some y → y.1.deg = 0 + x.1.deg := by
    intro x y hxy
    dsimp [G] at hxy
    split at hxy
    · contradiction
    · injection hxy with hxy
      subst y
      simp
  have hs : SparsePolyZp.Sorted (SparsePolyZp.scalarMul c f) := by
    unfold SparsePolyZp.Sorted
    rw [hlist]
    exact (SparsePolyZp.sortedListB_filterMapShift 0 G hdeg f.toList).2 hf.1
  have hnz : SparsePolyZp.NonZeroB (SparsePolyZp.scalarMul c f) := by
    unfold SparsePolyZp.NonZeroB
    rw [hlist]
    apply nonzeroListB_filterMap_of_output G
    intro x y hxy
    dsimp [G] at hxy
    split at hxy
    · contradiction
    · rename_i hz
      injection hxy with hxy
      subst y
      exact hz
  have hred : SparsePolyZp.AllReduced p (SparsePolyZp.scalarMul c f).toList := by
    rw [hlist]
    intro y hy
    rcases List.mem_filterMap.mp hy with ⟨x, hx, hxy⟩
    dsimp [G] at hxy
    split at hxy
    · contradiction
    · rename_i hz
      injection hxy with hxy
      subst y
      have hxred := hf.2.2 x hx
      have hm := SparsePolyZp.Zp_mul_reduced x.2 c (by rw [hxred.1]; exact hp.out.pos)
      exact ⟨by
        change x.2.prime.toNat = p
        exact hxred.1, by
        rw [UInt64.lt_iff_toNat_lt] at hm
        simpa [hxred.1] using hm.1⟩
  exact ⟨hs, hnz, hred⟩

private lemma makeMonic_canonical (f : SparsePolyZp) (hf : CanonicalRep p f) :
    CanonicalRep p (SparsePolyZp.makeMonic f) := by
  unfold SparsePolyZp.makeMonic
  split_ifs with hempty
  · exact hf
  · have hpos : 0 < f.size := by
      rw [Array.isEmpty_iff_size_eq_zero] at hempty
      omega
    have hprime : f[0]!.snd.inv.prime.toNat = p := by
      have hmem := mem_getFirst_toList f hpos
      simpa [Zp.inv] using (hf.2.2 f[0]! hmem).1
    exact scalarMul_canonical f[0]!.snd.inv f hf hprime

private lemma wellFormed_of_canonical (f : SparsePolyZp) (hf : CanonicalRep p f) :
    SparsePolyZp.WellFormed p f := by
  intro x hx
  exact (hf.2.2 x (Array.mem_def.mp hx)).1

private lemma squarefreeMeasure_eq (f : SparsePolyZp) (hf : CanonicalRep p f) :
    Generated.squarefreeMeasure f =
      if f.isEmpty then 0 else (SparsePolyZp.toPoly p f).natDegree + 1 := by
  unfold Generated.squarefreeMeasure
  split_ifs with hempty
  · rfl
  · rw [(toPoly_head_data (p := p) f hf hempty).1]

private lemma get_deg_pos_iff (f : SparsePolyZp) (hf : CanonicalRep p f)
    (hdeg : ∀ x ∈ f.toList, x.1.deg < 2 ^ 63) :
    (get_deg f > (0 : Int64)) ↔ 0 < (SparsePolyZp.toPoly p f).natDegree := by
  by_cases hempty : f.isEmpty
  · have hfarr : f = #[] := Array.eq_empty_of_size_eq_zero
      (Array.isEmpty_iff_size_eq_zero.mp hempty)
    simp [get_deg, hempty, hfarr, SparsePolyZp.toPoly_empty]
  · have hpos : 0 < f.size := by
      rw [Array.isEmpty_iff_size_eq_zero] at hempty
      omega
    have hmem := mem_getFirst_toList f hpos
    have hheadlt : f[0]!.fst.deg < 2 ^ 63 := hdeg f[0]! hmem
    rw [(toPoly_head_data (p := p) f hf hempty).1]
    unfold get_deg
    rw [if_neg hempty]
    change ((0 : Int64) < f[0]!.fst.deg.toUInt64.toInt64) ↔ _
    rw [Int64.lt_iff_toInt_lt]
    change (0 : Int) < (Int64.ofNat f[0]!.fst.deg).toInt ↔ _
    rw [Int64.toInt_ofNat_of_lt hheadlt]
    exact_mod_cast (show 0 < f[0]!.fst.deg ↔ 0 < f[0]!.fst.deg from Iff.rfl)

private lemma toPoly_ne_zero_of_canonical_nonempty (f : SparsePolyZp)
    (hf : CanonicalRep p f) (hne : ¬f.isEmpty) : SparsePolyZp.toPoly p f ≠ 0 := by
  have hdata := toPoly_head_data (p := p) f hf hne
  have hpos : 0 < f.size := by
    rw [Array.isEmpty_iff_size_eq_zero] at hne
    omega
  have hmem := mem_getFirst_toList f hpos
  have hcoeffne := Zp.toZMod_ne_zero_of_val_ne_zero p f[0]!.snd
    (hf.2.2 f[0]! hmem) (by
      intro hz
      exact head_val_nonzero_of_nonZeroB f hne hf.2.1 (by simp [hz]))
  intro hzero
  have hc := hdata.2
  rw [hzero] at hc
  simp at hc
  exact hcoeffne hc.symm

private lemma nonempty_of_toPoly_ne_zero (f : SparsePolyZp)
    (hne : SparsePolyZp.toPoly p f ≠ 0) : ¬f.isEmpty := by
  intro hempty
  have hfarr : f = #[] := Array.eq_empty_of_size_eq_zero
    (Array.isEmpty_iff_size_eq_zero.mp hempty)
  apply hne
  simp [hfarr, SparsePolyZp.toPoly_empty]

private lemma canonical_deg_bound_of_natDegree_lt (f : SparsePolyZp)
    (hf : CanonicalRep p f) {B : Nat}
    (hB : (SparsePolyZp.toPoly p f).natDegree < B) :
    ∀ x ∈ f.toList, x.1.deg < B := by
  by_cases hempty : f.isEmpty
  · have hfarr : f = #[] := Array.eq_empty_of_size_eq_zero
      (Array.isEmpty_iff_size_eq_zero.mp hempty)
    simp [hfarr]
  · have hdata := toPoly_head_data (p := p) f hf hempty
    have hlist := SparsePolyZp.toList_cons_of_ne_empty f hempty
    have hs := (SparsePolyZp.sortedListB_iff f[0]! f.toList.tail).mp (hlist ▸ hf.1)
    intro x hx
    rcases List.mem_cons.mp (hlist ▸ hx) with rfl | hxt
    · simpa [hdata.1] using hB
    · exact lt_trans (hs.1 x hxt) (by simpa [hdata.1] using hB)

private lemma squarefree_loop_cond_iff (f : SparsePolyZp) (hf : CanonicalRep p f)
    (hdeg : ∀ x ∈ f.toList, x.1.deg < 2 ^ 63) :
    ((!f.isEmpty) && decide (get_deg f > (0 : Int64))) = true ↔
      0 < (SparsePolyZp.toPoly p f).natDegree := by
  by_cases hempty : f.isEmpty
  · have hfarr : f = #[] := Array.eq_empty_of_size_eq_zero
      (Array.isEmpty_iff_size_eq_zero.mp hempty)
    simp [hempty, hfarr, SparsePolyZp.toPoly_empty]
  · simp [hempty, (get_deg_pos_iff (p := p) f hf hdeg)]

private lemma normalize_ne_zero_iff_local {R : Type*} [CommMonoidWithZero R]
    [NormalizationMonoid R] {a : R} : normalize a ≠ 0 ↔ a ≠ 0 := by
  rw [ne_eq, normalize_eq_zero, ne_eq]

private lemma natDegree_normalize_eq_local (f : Polynomial (ZMod p)) :
    (normalize f).natDegree = f.natDegree := by
  rcases eq_or_ne f 0 with rfl | hne
  · simp
  · have hne' : normalize f ≠ 0 := normalize_ne_zero_iff_local.mpr hne
    have hdeg := degree_eq_degree_of_associated (normalize_associated f)
    rw [Polynomial.degree_eq_natDegree hne, Polynomial.degree_eq_natDegree hne'] at hdeg
    exact Nat.cast_injective hdeg

private lemma divByMonic_ne_zero_local (f g : Polynomial (ZMod p))
    (hg : Monic g) (hdvd : g ∣ f) (hf : f ≠ 0) : f /ₘ g ≠ 0 := by
  intro h
  have hid := modByMonic_add_div f g
  rw [(modByMonic_eq_zero_iff_dvd hg).mpr hdvd, h, mul_zero, zero_add] at hid
  exact hf hid.symm

private lemma normalize_eq_of_monic_local (f : Polynomial (ZMod p))
    (hf : Monic f) : normalize f = f := by
  symm
  exact Polynomial.eq_of_monic_of_associated hf
    (Polynomial.monic_normalize hf.ne_zero) (normalize_associated f).symm

/-- `divmodAux` 的余式保留实现层三个规范形不变量。 -/
private theorem divmodAux_snd_invariants (g : SparsePolyZp) (dg : Nat)
    (lc_g_inv : Zp) (pm : UInt64) (hq : 0 < pm.toNat)
    (hg_ne : ¬ g.isEmpty) (hg_red : SparsePolyZp.ReducedB g pm)
    (h_dg : (g[0]!).fst.deg = dg) (hlp : lc_g_inv.prime = pm)
    (h_lc : (lc_g_inv * (g[0]!).snd).val = 1) (h_sorted_g : SparsePolyZp.Sorted g) :
    ∀ (q r : SparsePolyZp) (h_sorted_r : SparsePolyZp.Sorted r)
      (hr_red : SparsePolyZp.ReducedB r pm) (hr_nz : SparsePolyZp.NonZeroB r),
      let out := SparsePolyZp.divmodAux g dg lc_g_inv pm hq hg_ne hg_red h_dg hlp h_lc
        h_sorted_g q r h_sorted_r hr_red hr_nz
      SparsePolyZp.Sorted out.2 ∧ SparsePolyZp.ReducedB out.2 pm ∧
        SparsePolyZp.NonZeroB out.2 := by
  let motive : ∀ (q r : SparsePolyZp), SparsePolyZp.Sorted r →
      SparsePolyZp.ReducedB r pm → SparsePolyZp.NonZeroB r → Prop :=
    fun q r hsr hrr hrnz =>
      let out := SparsePolyZp.divmodAux g dg lc_g_inv pm hq hg_ne hg_red h_dg hlp h_lc
        h_sorted_g q r hsr hrr hrnz
      SparsePolyZp.Sorted out.2 ∧ SparsePolyZp.ReducedB out.2 pm ∧ SparsePolyZp.NonZeroB out.2
  apply SparsePolyZp.divmodAux.induct g dg lc_g_inv pm hq hg_ne hg_red h_dg hlp h_lc
    h_sorted_g motive
  · intro q r hsr hrr hrnz hr
    dsimp only [motive]
    rw [SparsePolyZp.divmodAux]
    simp only [hr, dif_pos]
    exact ⟨hsr, hrr, hrnz⟩
  · intro q r hsr hrr hrnz hr
    dsimp only
    intro hd
    dsimp only [motive]
    rw [SparsePolyZp.divmodAux]
    simp only [hr, hd, dif_neg, dif_pos]
    exact ⟨hsr, hrr, hrnz⟩
  · intro q r hsr hrr hrnz hr
    dsimp only
    intro hnotlt hsr' hcoeff_prime htgred htgnz hrr' hrnz' ih
    dsimp only [motive] at ih ⊢
    rw [SparsePolyZp.divmodAux]
    simp only [hr, hnotlt, dif_neg]
    exact ih

/-- `divmodAux` 从空商开始时也保持商的规范形。额外的 `orderInv` 记录
每次追加的商次数严格递减。 -/
private theorem divmodAux_fst_invariants (g : SparsePolyZp) (dg : Nat)
    (lc_g_inv : Zp) (pm : UInt64) (hq : 0 < pm.toNat)
    (hg_ne : ¬ g.isEmpty) (hg_red : SparsePolyZp.ReducedB g pm)
    (h_dg : (g[0]!).fst.deg = dg) (hlp : lc_g_inv.prime = pm)
    (h_lc : (lc_g_inv * (g[0]!).snd).val = 1) (h_sorted_g : SparsePolyZp.Sorted g)
    (hpm : pm.toNat = p) (hlc_red : Zp.Reduced p lc_g_inv)
    (h2p : 2 * p ≤ UInt64.size) :
    ∀ (q r : SparsePolyZp) (h_sorted_r : SparsePolyZp.Sorted r)
      (hr_red : SparsePolyZp.ReducedB r pm) (hr_nz : SparsePolyZp.NonZeroB r),
      SparsePolyZp.Sorted q → SparsePolyZp.ReducedB q pm → SparsePolyZp.NonZeroB q →
      (∀ x ∈ q.toList, ¬r.isEmpty → dg ≤ r[0]!.fst.deg →
        r[0]!.fst.deg - dg < x.1.deg) →
      let out := SparsePolyZp.divmodAux g dg lc_g_inv pm hq hg_ne hg_red h_dg hlp h_lc
        h_sorted_g q r h_sorted_r hr_red hr_nz
      SparsePolyZp.Sorted out.1 ∧ SparsePolyZp.ReducedB out.1 pm ∧
        SparsePolyZp.NonZeroB out.1 := by
  let motive : ∀ (q r : SparsePolyZp), SparsePolyZp.Sorted r →
      SparsePolyZp.ReducedB r pm → SparsePolyZp.NonZeroB r → Prop :=
    fun q r hsr hrr hrnz =>
      SparsePolyZp.Sorted q → SparsePolyZp.ReducedB q pm → SparsePolyZp.NonZeroB q →
      (∀ x ∈ q.toList, ¬r.isEmpty → dg ≤ r[0]!.fst.deg →
        r[0]!.fst.deg - dg < x.1.deg) →
      let out := SparsePolyZp.divmodAux g dg lc_g_inv pm hq hg_ne hg_red h_dg hlp h_lc
        h_sorted_g q r hsr hrr hrnz
      SparsePolyZp.Sorted out.1 ∧ SparsePolyZp.ReducedB out.1 pm ∧
        SparsePolyZp.NonZeroB out.1
  apply SparsePolyZp.divmodAux.induct g dg lc_g_inv pm hq hg_ne hg_red h_dg hlp h_lc
    h_sorted_g motive
  · intro q r hsr hrr hrnz hr hqs hqr hqnz hord
    dsimp only [motive]
    rw [SparsePolyZp.divmodAux]
    simp only [hr, dif_pos]
    exact ⟨hqs, hqr, hqnz⟩
  · intro q r hsr hrr hrnz hr
    dsimp only
    intro hd hqs hqr hqnz hord
    dsimp only [motive]
    rw [SparsePolyZp.divmodAux]
    simp only [hr, hd, dif_neg, dif_pos]
    exact ⟨hqs, hqr, hqnz⟩
  · intro q r hsr hrr hrnz hr
    dsimp only
    intro hnotlt hsr' hcoeff_prime htgred htgnz hrr' hrnz' ih
      hqs hqr hqnz hord
    dsimp only [motive] at ih ⊢
    let coeff : Zp := r[0]!.snd * lc_g_inv
    let d : Nat := r[0]!.fst.deg - dg
    let term : SparsePolyZp := #[(⟨d⟩, coeff)]
    let r' : SparsePolyZp := r - term * g
    let q' : SparsePolyZp := q.push (⟨d⟩, coeff)
    have hrlist := SparsePolyZp.toList_cons_of_ne_empty r hr
    have hrhead_pm : r[0]!.snd.prime = pm :=
      ((SparsePolyZp.reducedListB_cons pm r[0]! r.toList.tail).mp (hrlist ▸ hrr)).1.2
    have hrhead_red : Zp.Reduced p r[0]!.snd := by
      have hv := ((SparsePolyZp.reducedListB_cons pm r[0]! r.toList.tail).mp
        (hrlist ▸ hrr)).1.1
      exact ⟨by simpa [hrhead_pm] using hpm, by
        rw [UInt64.lt_iff_toNat_lt] at hv
        simpa [hrhead_pm, hpm] using hv⟩
    have hrhead_nz : r[0]!.snd.val.toNat ≠ 0 :=
      head_val_nonzero_of_nonZeroB r hr hrnz
    have hcoeff_red : Zp.Reduced p coeff := by
      have hm := SparsePolyZp.Zp_mul_reduced r[0]!.snd lc_g_inv (by simpa [hrhead_pm] using hq)
      exact ⟨by
        change r[0]!.snd.prime.toNat = p
        rw [hrhead_pm, hpm], by
        rw [UInt64.lt_iff_toNat_lt] at hm
        simpa [coeff, hrhead_pm, hpm] using hm.1⟩
    have hlc_nz : lc_g_inv.val.toNat ≠ 0 := by
      intro hz
      change (((lc_g_inv.val.toNat * g[0]!.snd.val.toNat) %
        lc_g_inv.prime.toNat).toUInt64) = 1 at h_lc
      rw [hz] at h_lc
      simp at h_lc
    have hcoeff_nz_nat : coeff.val.toNat ≠ 0 :=
      Zp_mul_val_nonzero hrhead_red hlc_red hrhead_nz hlc_nz h2p
    have hcoeff_nz : coeff.val ≠ 0 := by
      intro hz
      exact hcoeff_nz_nat (by simp [hz])
    have hqs' : SparsePolyZp.Sorted q' := by
      unfold SparsePolyZp.Sorted
      rw [show q'.toList = q.toList ++ [(⟨d⟩, coeff)] by simp [q']]
      exact sortedListB_append_singleton q.toList (⟨d⟩, coeff) hqs
        (fun x hx => hord x hx hr (Nat.le_of_not_lt hnotlt))
    have hqr' : SparsePolyZp.ReducedB q' pm := by
      apply reducedB_of_allReduced q' pm hpm
      intro x hx
      have hx' : x ∈ q.toList ∨ x = (⟨d⟩, coeff) := by simpa [q'] using hx
      rcases hx' with hx | rfl
      · exact allReduced_of_reducedB q pm hpm hqr x hx
      · exact hcoeff_red
    have hqnz' : SparsePolyZp.NonZeroB q' := by
      unfold SparsePolyZp.NonZeroB
      rw [show q'.toList = q.toList ++ [(⟨d⟩, coeff)] by simp [q']]
      exact nonzeroListB_append_singleton q.toList (⟨d⟩, coeff) hqnz hcoeff_nz
    have hdrop : ∀ z ∈ r'.toList, z.fst.deg < r[0]!.fst.deg := by
      simpa [r', term, coeff, d] using
        (SparsePolyZp.divmod_step_drop g r dg lc_g_inv pm hq hg_ne h_sorted_g hg_red
          h_dg hlp h_lc hr hsr hrr hrnz (Nat.le_of_not_lt hnotlt))
    have hord' : ∀ x ∈ q'.toList, ¬r'.isEmpty → dg ≤ r'[0]!.fst.deg →
        r'[0]!.fst.deg - dg < x.1.deg := by
      intro x hx hr'ne hdg'
      have hr'pos : 0 < r'.size := by
        rw [Array.isEmpty_iff_size_eq_zero] at hr'ne
        omega
      have hlead_drop := hdrop r'[0]! (mem_getFirst_toList r' hr'pos)
      have hx' : x ∈ q.toList ∨ x = (⟨d⟩, coeff) := by simpa [q'] using hx
      rcases hx' with hxold | rfl
      · have hold := hord x hxold hr (Nat.le_of_not_lt hnotlt)
        omega
      · change r'[0]!.fst.deg - dg < d
        simp only [d]
        omega
    rw [SparsePolyZp.divmodAux]
    simp only [hr, hnotlt, dif_neg]
    exact ih hqs' hqr' hqnz' hord'

/-- `divmod` 合法分支的余式规范形和首项次数界。 -/
private theorem divmod_snd_invariants_and_deg (f g : SparsePolyZp)
    (hvalid : ¬ g.isEmpty ∧ ((g[0]!.snd.inv * (g[0]!).snd).val = 1) ∧
      SparsePolyZp.Sorted g ∧ SparsePolyZp.Sorted f ∧
      0 < (g[0]!.snd.prime).toNat ∧ SparsePolyZp.ReducedB g (g[0]!.snd.prime) ∧
      SparsePolyZp.ReducedB f (g[0]!.snd.prime) ∧ SparsePolyZp.NonZeroB f) :
    let r := (SparsePolyZp.divmod f g).2
    SparsePolyZp.Sorted r ∧ SparsePolyZp.ReducedB r (g[0]!.snd.prime) ∧
      SparsePolyZp.NonZeroB r ∧
      (r.isEmpty = true ∨ r[0]!.fst.deg < g[0]!.fst.deg) := by
  have hinv := divmodAux_snd_invariants g (g[0]!.fst.deg) (g[0]!.snd.inv)
    (g[0]!.snd.prime) hvalid.2.2.2.2.1 hvalid.1 hvalid.2.2.2.2.2.1 rfl rfl
    hvalid.2.1 hvalid.2.2.1 #[] f hvalid.2.2.2.1
    hvalid.2.2.2.2.2.2.1 hvalid.2.2.2.2.2.2.2
  have hdeg := SparsePolyZp.divmodAux_snd_deg_lt g (g[0]!.fst.deg) (g[0]!.snd.inv)
    (g[0]!.snd.prime) hvalid.2.2.2.2.1 hvalid.1 hvalid.2.2.2.2.2.1 rfl rfl
    hvalid.2.1 hvalid.2.2.1 #[] f hvalid.2.2.2.1
    hvalid.2.2.2.2.2.2.1 hvalid.2.2.2.2.2.2.2
  unfold SparsePolyZp.divmod
  rw [dif_pos hvalid]
  exact ⟨hinv.1, hinv.2.1, hinv.2.2, hdeg⟩

private lemma divmod_snd_canonical (f g : SparsePolyZp)
    (hf : CanonicalRep p f) (hg : CanonicalRep p g) (hg_ne : ¬g.isEmpty)
    (hp_lt : p < UInt64.size) : CanonicalRep p (SparsePolyZp.divmod f g).2 := by
  have hvalid := divmod_valid_of_canonical f g hg_ne hf hg hp_lt
  have hi := divmod_snd_invariants_and_deg f g hvalid
  have hmem : g[0]! ∈ g.toList := by
    have hpos : 0 < g.size := by
      rw [Array.isEmpty_iff_size_eq_zero] at hg_ne
      omega
    exact mem_getFirst_toList g hpos
  have hpm : g[0]!.snd.prime.toNat = p := (hg.2.2 g[0]! hmem).1
  exact ⟨hi.1, hi.2.2.1, allReduced_of_reducedB _ _ hpm hi.2.1⟩

private lemma divmod_fst_canonical (f g : SparsePolyZp)
    (hf : CanonicalRep p f) (hg : CanonicalRep p g) (hg_ne : ¬g.isEmpty)
    (h2p : 2 * p ≤ UInt64.size) : CanonicalRep p (SparsePolyZp.divmod f g).1 := by
  have hp_lt : p < UInt64.size := by nlinarith
  have hvalid := divmod_valid_of_canonical f g hg_ne hf hg hp_lt
  have hpos : 0 < g.size := by
    rw [Array.isEmpty_iff_size_eq_zero] at hg_ne
    omega
  have hmem : g[0]! ∈ g.toList := mem_getFirst_toList g hpos
  have hgred : Zp.Reduced p g[0]!.snd := hg.2.2 g[0]! hmem
  have hpm : g[0]!.snd.prime.toNat = p := hgred.1
  have hlcred : Zp.Reduced p g[0]!.snd.inv := Zp_inv_reduced g[0]!.snd hgred hp_lt
  have hi := divmodAux_fst_invariants g (g[0]!.fst.deg) (g[0]!.snd.inv)
    (g[0]!.snd.prime) hvalid.2.2.2.2.1 hvalid.1 hvalid.2.2.2.2.2.1 rfl rfl
    hvalid.2.1 hvalid.2.2.1 hpm hlcred h2p #[] f hvalid.2.2.2.1
    hvalid.2.2.2.2.2.2.1 hvalid.2.2.2.2.2.2.2
    (by rfl) (by rfl) (by rfl) (by simp)
  unfold SparsePolyZp.divmod
  rw [dif_pos hvalid]
  exact ⟨hi.1, hi.2.2, allReduced_of_reducedB _ _ hpm hi.2.1⟩

/-- 实现层 Euclid 循环与数学侧 GCD 相伴，并保持规范形。 -/
private theorem gcdAux_refines (h2p : 2 * p ≤ UInt64.size) (hp2 : p * p ≤ UInt64.size) :
    ∀ (f g : SparsePolyZp), CanonicalRep p f → CanonicalRep p g →
      CanonicalRep p (SparsePolyZp.gcdAux f g) ∧
      Associated (SparsePolyZp.toPoly p (SparsePolyZp.gcdAux f g))
        (EuclideanDomain.gcd (SparsePolyZp.toPoly p f) (SparsePolyZp.toPoly p g)) := by
  have hp_lt : p < UInt64.size := by nlinarith
  let motive : SparsePolyZp → SparsePolyZp → Prop := fun f g =>
    CanonicalRep p f → CanonicalRep p g →
      CanonicalRep p (SparsePolyZp.gcdAux f g) ∧
      Associated (SparsePolyZp.toPoly p (SparsePolyZp.gcdAux f g))
        (EuclideanDomain.gcd (SparsePolyZp.toPoly p f) (SparsePolyZp.toPoly p g))
  apply SparsePolyZp.gcdAux.induct motive
  · intro f g hg_empty hf hg
    have hg_arr : g = #[] := Array.eq_empty_of_size_eq_zero
      (Array.isEmpty_iff_size_eq_zero.mp hg_empty)
    rw [SparsePolyZp.gcdAux]
    simp only [hg_empty, dif_pos]
    refine ⟨hf, ?_⟩
    rw [hg_arr, SparsePolyZp.toPoly_empty]
    exact Associated.of_eq (EuclideanDomain.gcd_zero_right (SparsePolyZp.toPoly p f)).symm
  · intro f g hg_ne
    dsimp only
    intro hr_empty hf hg
    have hvalid := divmod_valid_of_canonical f g hg_ne hf hg hp_lt
    have hpm : g[0]!.snd.prime.toNat = p := by
      have hpos : 0 < g.size := by
        rw [Array.isEmpty_iff_size_eq_zero] at hg_ne
        omega
      exact (hg.2.2 g[0]! (mem_getFirst_toList g hpos)).1
    have hid := divmod_toPoly_identity f g hvalid hpm h2p hp2
    have hr_arr : (SparsePolyZp.divmod f g).2 = #[] := Array.eq_empty_of_size_eq_zero
      (Array.isEmpty_iff_size_eq_zero.mp hr_empty)
    have hdvd : SparsePolyZp.toPoly p g ∣ SparsePolyZp.toPoly p f := by
      refine ⟨SparsePolyZp.toPoly p (SparsePolyZp.divmod f g).1, ?_⟩
      rw [hid, hr_arr, SparsePolyZp.toPoly_empty, add_zero, mul_comm]
    rw [SparsePolyZp.gcdAux]
    simp only [hg_ne, hr_empty, dif_neg, dif_pos]
    exact ⟨hg, associated_of_dvd_dvd
      (EuclideanDomain.dvd_gcd hdvd dvd_rfl)
      (EuclideanDomain.gcd_dvd_right _ _)⟩
  · intro f g hg_ne
    dsimp only
    intro hr_ne hdeg ih hf hg
    dsimp only [motive] at ih ⊢
    let r := (SparsePolyZp.divmod f g).2
    have hvalid := divmod_valid_of_canonical f g hg_ne hf hg hp_lt
    have hrcan : CanonicalRep p r := divmod_snd_canonical f g hf hg hg_ne hp_lt
    have hih := ih hg hrcan
    have hpm : g[0]!.snd.prime.toNat = p := by
      have hpos : 0 < g.size := by
        rw [Array.isEmpty_iff_size_eq_zero] at hg_ne
        omega
      exact (hg.2.2 g[0]! (mem_getFirst_toList g hpos)).1
    have hid := divmod_toPoly_identity f g hvalid hpm h2p hp2
    let F := SparsePolyZp.toPoly p f
    let G := SparsePolyZp.toPoly p g
    let R := SparsePolyZp.toPoly p r
    let Q := SparsePolyZp.toPoly p (SparsePolyZp.divmod f g).1
    have hdiv : F = Q * G + R := by simpa [F, G, R, Q, r] using hid
    have hfg_dvd_r : EuclideanDomain.gcd F G ∣ R := by
      have hsub : R = F - Q * G := by rw [hdiv]; ring
      rw [hsub]
      exact dvd_sub (EuclideanDomain.gcd_dvd_left F G)
        (dvd_mul_of_dvd_right (EuclideanDomain.gcd_dvd_right F G) Q)
    have hgr_dvd_f : EuclideanDomain.gcd G R ∣ F := by
      rw [hdiv]
      exact dvd_add (dvd_mul_of_dvd_right (EuclideanDomain.gcd_dvd_left G R) Q)
        (EuclideanDomain.gcd_dvd_right G R)
    have hgcd_assoc : Associated (EuclideanDomain.gcd G R) (EuclideanDomain.gcd F G) :=
      associated_of_dvd_dvd
        (EuclideanDomain.dvd_gcd hgr_dvd_f (EuclideanDomain.gcd_dvd_left G R))
        (EuclideanDomain.dvd_gcd (EuclideanDomain.gcd_dvd_right F G) hfg_dvd_r)
    rw [SparsePolyZp.gcdAux]
    simp only [hg_ne, hr_ne, hdeg, dif_neg]
    exact ⟨hih.1, hih.2.trans (by simpa [F, G, R] using hgcd_assoc)⟩
  · intro f g hg_ne
    dsimp only
    intro hr_ne hnotdeg hf hg
    have hvalid := divmod_valid_of_canonical f g hg_ne hf hg hp_lt
    have hi := divmod_snd_invariants_and_deg f g hvalid
    rcases hi.2.2.2 with hempty | hdeg
    · exact (hr_ne hempty).elim
    · exact (hnotdeg hdeg).elim

private theorem makeMonic_toPoly_monic (h2p : 2 * p ≤ UInt64.size)
    (f : SparsePolyZp) (hf : CanonicalRep p f) (hne : ¬f.isEmpty) :
    Monic (SparsePolyZp.toPoly p (SparsePolyZp.makeMonic f)) := by
  have hdata := toPoly_head_data (p := p) f hf hne
  have hpos : 0 < f.size := by
    rw [Array.isEmpty_iff_size_eq_zero] at hne
    omega
  have hmem : f[0]! ∈ f.toList := mem_getFirst_toList f hpos
  have hlcred : Zp.Reduced p f[0]!.snd := hf.2.2 f[0]! hmem
  have hlcnz : f[0]!.snd.val.toNat ≠ 0 :=
    head_val_nonzero_of_nonZeroB f hne hf.2.1
  have hinvprime : f[0]!.snd.inv.prime.toNat = p := by
    simpa [Zp.inv] using hlcred.1
  unfold SparsePolyZp.makeMonic
  rw [if_neg hne, toPoly_scalarMul f[0]!.snd.inv f hf.2.2 hinvprime h2p]
  rw [Polynomial.Monic.def, Polynomial.leadingCoeff_mul]
  simp only [Polynomial.leadingCoeff_C]
  rw [Polynomial.leadingCoeff, hdata.1, hdata.2]
  have hp_lt : p < UInt64.size := by nlinarith
  have hone := Zp_toZMod_inv_mul_self f[0]!.snd hlcred hlcnz hp_lt
  have hsplit := Zp.toZMod_mul_weak f[0]!.snd f[0]!.snd.inv hlcred.1 hinvprime h2p
  rw [hsplit] at hone
  simpa [mul_comm] using hone

private theorem makeMonic_toPoly_eq_normalize (h2p : 2 * p ≤ UInt64.size)
    (f : SparsePolyZp) (hf : CanonicalRep p f) (hne : ¬f.isEmpty) :
    SparsePolyZp.toPoly p (SparsePolyZp.makeMonic f) =
      normalize (SparsePolyZp.toPoly p f) := by
  have houtmonic := makeMonic_toPoly_monic (p := p) h2p f hf hne
  have hdata := toPoly_head_data (p := p) f hf hne
  have hfinne : SparsePolyZp.toPoly p f ≠ 0 := by
    intro hz
    have hpos : 0 < f.size := by
      rw [Array.isEmpty_iff_size_eq_zero] at hne
      omega
    have hmem := mem_getFirst_toList f hpos
    have hcoeff_ne := Zp.toZMod_ne_zero_of_val_ne_zero p f[0]!.snd
      (hf.2.2 f[0]! hmem) (by
        intro hzv
        exact head_val_nonzero_of_nonZeroB f hne hf.2.1 (by simp [hzv]))
    have hc := hdata.2
    rw [hz] at hc
    simp at hc
    exact hcoeff_ne hc.symm
  have hnormmonic : Monic (normalize (SparsePolyZp.toPoly p f)) :=
    Polynomial.monic_normalize hfinne
  have hassoc : Associated (SparsePolyZp.toPoly p (SparsePolyZp.makeMonic f))
      (SparsePolyZp.toPoly p f) := by
    unfold SparsePolyZp.makeMonic
    rw [if_neg hne]
    have hpos : 0 < f.size := by
      rw [Array.isEmpty_iff_size_eq_zero] at hne
      omega
    have hmem := mem_getFirst_toList f hpos
    have hlcred := hf.2.2 f[0]! hmem
    have hinvprime : f[0]!.snd.inv.prime.toNat = p := by
      simpa [Zp.inv] using hlcred.1
    rw [toPoly_scalarMul f[0]!.snd.inv f hf.2.2 hinvprime h2p]
    have hcne : Zp.toZMod p f[0]!.snd.inv ≠ 0 := by
      have hp_lt : p < UInt64.size := by nlinarith
      have hlcnz := head_val_nonzero_of_nonZeroB f hne hf.2.1
      have hone := Zp_toZMod_inv_mul_self f[0]!.snd hlcred hlcnz hp_lt
      intro hz
      have hsplit := Zp.toZMod_mul_weak f[0]!.snd f[0]!.snd.inv
        hlcred.1 hinvprime h2p
      rw [hsplit, hz] at hone
      simp at hone
    have hcunit : IsUnit (Polynomial.C (Zp.toZMod p f[0]!.snd.inv)) :=
      Polynomial.isUnit_C.mpr ((isUnit_iff_ne_zero).mpr hcne)
    obtain ⟨u, hu⟩ := hcunit
    rw [← hu]
    simpa using ((show Associated (↑u) 1 from unit_associated_one).mul_right
      (SparsePolyZp.toPoly p f))
  exact Polynomial.eq_of_monic_of_associated houtmonic hnormmonic
    (hassoc.trans (normalize_associated _).symm)

/-- `polynomial_GCD` 的首一化只乘以非零常数，因此与数学 GCD 相伴。 -/
private theorem polynomial_GCD_refines (h2p : 2 * p ≤ UInt64.size)
    (hp2 : p * p ≤ UInt64.size) (f g : SparsePolyZp)
    (hf : CanonicalRep p f) (hg : CanonicalRep p g) :
    Associated (SparsePolyZp.toPoly p (polynomial_GCD f g))
      (EuclideanDomain.gcd (SparsePolyZp.toPoly p f) (SparsePolyZp.toPoly p g)) := by
  let d := SparsePolyZp.gcdAux f g
  have hd := gcdAux_refines (p := p) h2p hp2 f g hf hg
  have hdcan : CanonicalRep p d := by simpa [d] using hd.1
  have hdassoc : Associated (SparsePolyZp.toPoly p d)
      (EuclideanDomain.gcd (SparsePolyZp.toPoly p f) (SparsePolyZp.toPoly p g)) := by
    simpa [d] using hd.2
  change Associated (SparsePolyZp.toPoly p (SparsePolyZp.makeMonic d)) _
  unfold SparsePolyZp.makeMonic
  split_ifs with hempty
  · exact hdassoc
  · have hpos : 0 < d.size := by
      rw [Array.isEmpty_iff_size_eq_zero] at hempty
      omega
    have hmem : d[0]! ∈ d.toList := mem_getFirst_toList d hpos
    have hlcred : Zp.Reduced p d[0]!.snd := hdcan.2.2 d[0]! hmem
    have hlcnz : d[0]!.snd.val.toNat ≠ 0 := by
      exact head_val_nonzero_of_nonZeroB d hempty hdcan.2.1
    have hinvprime : d[0]!.snd.inv.prime.toNat = p := by
      simpa [Zp.inv] using hlcred.1
    rw [toPoly_scalarMul d[0]!.snd.inv d hdcan.2.2 hinvprime h2p]
    have hcne : Zp.toZMod p d[0]!.snd.inv ≠ 0 := by
      have hp_lt : p < UInt64.size := by nlinarith
      have hone := Zp_toZMod_inv_mul_self d[0]!.snd hlcred hlcnz hp_lt
      intro hzero
      have hprime : d[0]!.snd.prime.toNat = p := hlcred.1
      have hsplit := Zp.toZMod_mul_weak d[0]!.snd d[0]!.snd.inv
        hprime hinvprime h2p
      rw [hsplit, hzero] at hone
      simp at hone
    have hcunit : IsUnit (Polynomial.C (Zp.toZMod p d[0]!.snd.inv)) :=
      Polynomial.isUnit_C.mpr ((isUnit_iff_ne_zero).mpr hcne)
    obtain ⟨u, hu⟩ := hcunit
    have hscale : Associated
        (Polynomial.C (Zp.toZMod p d[0]!.snd.inv) * SparsePolyZp.toPoly p d)
        (SparsePolyZp.toPoly p d) := by
      rw [← hu]
      simpa using ((show Associated (↑u) 1 from unit_associated_one).mul_right
        (SparsePolyZp.toPoly p d))
    exact hscale.trans hdassoc

private theorem polynomial_GCD_canonical (h2p : 2 * p ≤ UInt64.size)
    (hp2 : p * p ≤ UInt64.size) (f g : SparsePolyZp)
    (hf : CanonicalRep p f) (hg : CanonicalRep p g) :
    CanonicalRep p (polynomial_GCD f g) := by
  have hd := (gcdAux_refines (p := p) h2p hp2 f g hf hg).1
  change CanonicalRep p (SparsePolyZp.makeMonic (SparsePolyZp.gcdAux f g))
  exact makeMonic_canonical _ hd

/-- 非零情形下，实现 GCD 的具体首一代表就是数学侧的 `normalize gcd`。 -/
private theorem polynomial_GCD_toPoly_eq_normalize (h2p : 2 * p ≤ UInt64.size)
    (hp2 : p * p ≤ UInt64.size) (f g : SparsePolyZp)
    (hf : CanonicalRep p f) (hg : CanonicalRep p g)
    (hne : ¬(SparsePolyZp.gcdAux f g).isEmpty) :
    SparsePolyZp.toPoly p (polynomial_GCD f g) =
      normalize (EuclideanDomain.gcd (SparsePolyZp.toPoly p f) (SparsePolyZp.toPoly p g)) := by
  let d := SparsePolyZp.gcdAux f g
  have hd := gcdAux_refines (p := p) h2p hp2 f g hf hg
  have hdcan : CanonicalRep p d := by simpa [d] using hd.1
  have hdne : ¬d.isEmpty := by simpa [d] using hne
  have houtmonic : Monic (SparsePolyZp.toPoly p (polynomial_GCD f g)) := by
    change Monic (SparsePolyZp.toPoly p (SparsePolyZp.makeMonic d))
    exact makeMonic_toPoly_monic (p := p) h2p d hdcan hdne
  have hassoc := polynomial_GCD_refines (p := p) h2p hp2 f g hf hg
  have houtne : SparsePolyZp.toPoly p (polynomial_GCD f g) ≠ 0 :=
    houtmonic.ne_zero
  have hgcdne : EuclideanDomain.gcd (SparsePolyZp.toPoly p f)
      (SparsePolyZp.toPoly p g) ≠ 0 := hassoc.ne_zero_iff.mp houtne
  have hnormmonic : Monic (normalize (EuclideanDomain.gcd
      (SparsePolyZp.toPoly p f) (SparsePolyZp.toPoly p g))) :=
    Polynomial.monic_normalize hgcdne
  exact Polynomial.eq_of_monic_of_associated houtmonic hnormmonic
    (hassoc.trans (normalize_associated _).symm)

private theorem polynomial_GCD_step_data (h2p : 2 * p ≤ UInt64.size)
    (hp2 : p * p ≤ UInt64.size) (f g : SparsePolyZp)
    (hf : CanonicalRep p f) (hg : CanonicalRep p g)
    (hfpoly : SparsePolyZp.toPoly p f ≠ 0) :
    let y := polynomial_GCD f g
    CanonicalRep p y ∧ ¬y.isEmpty ∧
      SparsePolyZp.toPoly p y = normalize (EuclideanDomain.gcd
        (SparsePolyZp.toPoly p f) (SparsePolyZp.toPoly p g)) ∧
      Monic (SparsePolyZp.toPoly p y) := by
  let d := SparsePolyZp.gcdAux f g
  let y := polynomial_GCD f g
  have hd := gcdAux_refines (p := p) h2p hp2 f g hf hg
  have hgcdne : EuclideanDomain.gcd (SparsePolyZp.toPoly p f)
      (SparsePolyZp.toPoly p g) ≠ 0 := by
    intro hz
    have hdiv := EuclideanDomain.gcd_dvd_left (SparsePolyZp.toPoly p f)
      (SparsePolyZp.toPoly p g)
    rw [hz] at hdiv
    exact hfpoly (zero_dvd_iff.mp hdiv)
  have hdpolyne : SparsePolyZp.toPoly p d ≠ 0 := by
    exact (by simpa [d] using hd.2.ne_zero_iff.mpr hgcdne)
  have hdne : ¬d.isEmpty := nonempty_of_toPoly_ne_zero (p := p) d hdpolyne
  have hycan : CanonicalRep p y := by
    simpa [y] using polynomial_GCD_canonical (p := p) h2p hp2 f g hf hg
  have hyne : ¬y.isEmpty := by
    apply nonempty_of_toPoly_ne_zero (p := p)
    simpa [y] using (polynomial_GCD_refines (p := p) h2p hp2 f g hf hg).ne_zero_iff.mpr hgcdne
  have hyeq : SparsePolyZp.toPoly p y = normalize (EuclideanDomain.gcd
      (SparsePolyZp.toPoly p f) (SparsePolyZp.toPoly p g)) := by
    simpa [y, d] using polynomial_GCD_toPoly_eq_normalize (p := p) h2p hp2 f g hf hg hdne
  exact ⟨hycan, hyne, hyeq, hyeq ▸ Polynomial.monic_normalize hgcdne⟩

/-- 对规范输入及首一除数，实际 `divmod` 的商等于 mathlib 的首一长除商。 -/
private theorem divmod_fst_toPoly_eq_divByMonic (h2p : 2 * p ≤ UInt64.size)
    (hp2 : p * p ≤ UInt64.size) (f g : SparsePolyZp)
    (hf : CanonicalRep p f) (hg : CanonicalRep p g) (hg_ne : ¬g.isEmpty)
    (hg_monic : Monic (SparsePolyZp.toPoly p g))
    (hdvd : SparsePolyZp.toPoly p g ∣ SparsePolyZp.toPoly p f) :
    SparsePolyZp.toPoly p (SparsePolyZp.divmod f g).1 =
      SparsePolyZp.toPoly p f /ₘ SparsePolyZp.toPoly p g := by
  have hp_lt : p < UInt64.size := by nlinarith
  have hvalid := divmod_valid_of_canonical f g hg_ne hf hg hp_lt
  have hpos : 0 < g.size := by
    rw [Array.isEmpty_iff_size_eq_zero] at hg_ne
    omega
  have hpm : g[0]!.snd.prime.toNat = p :=
    (hg.2.2 g[0]! (mem_getFirst_toList g hpos)).1
  have hid := divmod_toPoly_identity f g hvalid hpm h2p hp2
  let q := (SparsePolyZp.divmod f g).1
  let r := (SparsePolyZp.divmod f g).2
  let F := SparsePolyZp.toPoly p f
  let G := SparsePolyZp.toPoly p g
  let Q := SparsePolyZp.toPoly p q
  let R := SparsePolyZp.toPoly p r
  have hiden : F = Q * G + R := by simpa [F, G, Q, R, q, r] using hid
  have hGdvdR : G ∣ R := by
    rcases hdvd with ⟨k, hk⟩
    refine ⟨k - Q, ?_⟩
    change F = G * k at hk
    calc
      R = F - Q * G := by rw [hiden]; ring
      _ = G * k - Q * G := by rw [hk]
      _ = G * (k - Q) := by ring
  have hrzero : R = 0 := by
    have hi := divmod_snd_invariants_and_deg f g hvalid
    rcases hi.2.2.2 with hr_empty | hr_deg
    · have hr_arr : r = #[] := Array.eq_empty_of_size_eq_zero
        (Array.isEmpty_iff_size_eq_zero.mp hr_empty)
      simp [R, hr_arr, SparsePolyZp.toPoly_empty]
    · by_cases hr_ne : ¬r.isEmpty
      · have hrcan : CanonicalRep p r := divmod_snd_canonical f g hf hg hg_ne hp_lt
        have hrdata := toPoly_head_data (p := p) r hrcan hr_ne
        have hgdata := toPoly_head_data (p := p) g hg hg_ne
        apply Polynomial.eq_zero_of_dvd_of_natDegree_lt hGdvdR
        simpa [R, G, r, hrdata.1, hgdata.1] using hr_deg
      · have hr_empty' : r.isEmpty := Decidable.not_not.mp hr_ne
        have hr_arr : r = #[] := Array.eq_empty_of_size_eq_zero
          (Array.isEmpty_iff_size_eq_zero.mp hr_empty')
        simp [R, hr_arr, SparsePolyZp.toPoly_empty]
  have hmul : G * Q = F := by
    rw [hiden, hrzero, add_zero, mul_comm]
  have huniq := div_modByMonic_unique Q 0 hg_monic
    ⟨by simpa [hmul], by
      simpa only [Polynomial.degree_zero, bot_lt_iff_ne_bot, Polynomial.degree_ne_bot]
        using hg_monic.ne_zero⟩
  simpa [Q, F, G, q] using huniq.1.symm

set_option maxHeartbeats 0 in
private theorem yunLoop_ir_refines (h2p : 2 * p ≤ UInt64.size)
    (hp2 : p * p ≤ UInt64.size) :
    ∀ (i : UInt64) (w c : SparsePolyZp) (result : Array (SparsePolyZp × UInt64)),
      ∀ (hw : CanonicalRep p w) (hc : CanonicalRep p c)
      (hdw : ∀ x ∈ w.toList, x.1.deg < 2 ^ 63)
      (hdc : ∀ x ∈ c.toList, x.1.deg < 2 ^ 63)
      (hcne : SparsePolyZp.toPoly p c ≠ 0)
      (hcmonic : Monic (SparsePolyZp.toPoly p c)),
      i.toNat + (SparsePolyZp.toPoly p w).natDegree +
          (SparsePolyZp.toPoly p c).natDegree + 2 < 2 ^ 64 →
      let out := Generated._loop___squarefree_Zp_1_ir_def i w c result
      CanonicalRep p out.2.1 ∧
      (∀ x ∈ out.2.1.toList, x.1.deg < 2 ^ 63) ∧
      SparsePolyZp.toPoly p out.2.1 ≠ 0 ∧
      Monic (SparsePolyZp.toPoly p out.2.1) ∧
      (toPolyList out.2.2 p, SparsePolyZp.toPoly p out.2.1) =
        yunLoop (SparsePolyZp.toPoly p w) (SparsePolyZp.toPoly p c) i.toNat
          (toPolyList result p) hcne := by
  let motive : UInt64 → SparsePolyZp → SparsePolyZp →
      Array (SparsePolyZp × UInt64) → Prop := fun i w c result =>
    ∀ (hw : CanonicalRep p w) (hc : CanonicalRep p c)
    (hdw : ∀ x ∈ w.toList, x.1.deg < 2 ^ 63)
    (hdc : ∀ x ∈ c.toList, x.1.deg < 2 ^ 63)
    (hcne : SparsePolyZp.toPoly p c ≠ 0)
    (hcmonic : Monic (SparsePolyZp.toPoly p c)),
    i.toNat + (SparsePolyZp.toPoly p w).natDegree +
        (SparsePolyZp.toPoly p c).natDegree + 2 < 2 ^ 64 →
    let out := Generated._loop___squarefree_Zp_1_ir_def i w c result
    CanonicalRep p out.2.1 ∧
    (∀ x ∈ out.2.1.toList, x.1.deg < 2 ^ 63) ∧
    SparsePolyZp.toPoly p out.2.1 ≠ 0 ∧
    Monic (SparsePolyZp.toPoly p out.2.1) ∧
    (toPolyList out.2.2 p, SparsePolyZp.toPoly p out.2.1) =
      yunLoop (SparsePolyZp.toPoly p w) (SparsePolyZp.toPoly p c) i.toNat
        (toPolyList result p) hcne
  apply Generated._loop___squarefree_Zp_1_ir_def.induct motive
  · intro i w c result hcond
    dsimp only
    intro hdec ih hw hc hdw hdc hcne hcmonic hibound
    dsimp only [motive] at ih ⊢
    have hwpos : 0 < (SparsePolyZp.toPoly p w).natDegree :=
      (squarefree_loop_cond_iff (p := p) w hw hdw).mp hcond
    have hwne : SparsePolyZp.toPoly p w ≠ 0 := by
      intro hz
      simp [hz] at hwpos
    let W := SparsePolyZp.toPoly p w
    let C := SparsePolyZp.toPoly p c
    let y := polynomial_GCD w c
    have hydata := polynomial_GCD_step_data (p := p) h2p hp2 w c hw hc hwne
    have hycan : CanonicalRep p y := by simpa [y] using hydata.1
    have hyne : ¬y.isEmpty := by simpa [y] using hydata.2.1
    have hYeq : SparsePolyZp.toPoly p y = normalize (EuclideanDomain.gcd W C) := by
      simpa [y, W, C] using hydata.2.2.1
    have hYmonic : Monic (SparsePolyZp.toPoly p y) := by simpa [y] using hydata.2.2.2
    have hYdvdW : SparsePolyZp.toPoly p y ∣ W := by
      rw [hYeq]
      exact normalize_dvd_iff.mpr (EuclideanDomain.gcd_dvd_left W C)
    have hYdvdC : SparsePolyZp.toPoly p y ∣ C := by
      rw [hYeq]
      exact normalize_dvd_iff.mpr (EuclideanDomain.gcd_dvd_right W C)
    let qw := (SparsePolyZp.divmod w y).1
    let qc := (SparsePolyZp.divmod c y).1
    have hqwcan : CanonicalRep p qw := by
      exact divmod_fst_canonical (p := p) w y hw hycan hyne h2p
    have hqccan : CanonicalRep p qc := by
      exact divmod_fst_canonical (p := p) c y hc hycan hyne h2p
    have hqwpoly : SparsePolyZp.toPoly p qw = W /ₘ SparsePolyZp.toPoly p y := by
      simpa [qw, W] using divmod_fst_toPoly_eq_divByMonic (p := p) h2p hp2
        w y hw hycan hyne hYmonic hYdvdW
    have hqcpoly : SparsePolyZp.toPoly p qc = C /ₘ SparsePolyZp.toPoly p y := by
      simpa [qc, C] using divmod_fst_toPoly_eq_divByMonic (p := p) h2p hp2
        c y hc hycan hyne hYmonic hYdvdC
    have hqwpoly_ne : SparsePolyZp.toPoly p qw ≠ 0 := by
      rw [hqwpoly]
      exact divByMonic_ne_zero_local W _ hYmonic hYdvdW hwne
    have hqcpoly_ne : SparsePolyZp.toPoly p qc ≠ 0 := by
      rw [hqcpoly]
      exact divByMonic_ne_zero_local C _ hYmonic hYdvdC (by simpa [C] using hcne)
    have hqwne : ¬qw.isEmpty := nonempty_of_toPoly_ne_zero (p := p) qw hqwpoly_ne
    have hqcne : ¬qc.isEmpty := nonempty_of_toPoly_ne_zero (p := p) qc hqcpoly_ne
    have hnormqw : SparsePolyZp.normalization qw = qw :=
      normalization_eq_of_nonZeroB qw hqwcan.2.1
    have hnormqc : SparsePolyZp.normalization qc = qc :=
      normalization_eq_of_nonZeroB qc hqccan.2.1
    have hw_nd_lt : W.natDegree < 2 ^ 63 := by
      have hwarrne := nonempty_of_toPoly_ne_zero (p := p) w hwne
      rw [(toPoly_head_data (p := p) w hw hwarrne).1]
      exact hdw w[0]! (mem_getFirst_toList w (by
        rw [Array.isEmpty_iff_size_eq_zero] at hwarrne
        omega))
    have hc_nd_lt : C.natDegree < 2 ^ 63 := by
      have hcarrne := nonempty_of_toPoly_ne_zero (p := p) c (by simpa [C] using hcne)
      rw [(toPoly_head_data (p := p) c hc hcarrne).1]
      exact hdc c[0]! (mem_getFirst_toList c (by
        rw [Array.isEmpty_iff_size_eq_zero] at hcarrne
        omega))
    have hYnd_leW : (SparsePolyZp.toPoly p y).natDegree ≤ W.natDegree :=
      Polynomial.natDegree_le_of_dvd hYdvdW hwne
    have hqwnd_le : (SparsePolyZp.toPoly p qw).natDegree ≤ W.natDegree := by
      rw [hqwpoly, Polynomial.natDegree_divByMonic W hYmonic]
      omega
    have hqcnd_le : (SparsePolyZp.toPoly p qc).natDegree ≤ C.natDegree := by
      rw [hqcpoly, Polynomial.natDegree_divByMonic C hYmonic]
      omega
    have hdy : ∀ x ∈ y.toList, x.1.deg < 2 ^ 63 :=
      canonical_deg_bound_of_natDegree_lt (p := p) y hycan (lt_of_le_of_lt hYnd_leW hw_nd_lt)
    have hdqw : ∀ x ∈ qw.toList, x.1.deg < 2 ^ 63 :=
      canonical_deg_bound_of_natDegree_lt (p := p) qw hqwcan (lt_of_le_of_lt hqwnd_le hw_nd_lt)
    have hdqc : ∀ x ∈ qc.toList, x.1.deg < 2 ^ 63 :=
      canonical_deg_bound_of_natDegree_lt (p := p) qc hqccan (lt_of_le_of_lt hqcnd_le hc_nd_lt)
    have hc_mul : SparsePolyZp.toPoly p y * SparsePolyZp.toPoly p qc = C := by
      have hid := modByMonic_add_div C (SparsePolyZp.toPoly p y)
      rw [(modByMonic_eq_zero_iff_dvd hYmonic).mpr hYdvdC, zero_add] at hid
      simpa [hqcpoly] using hid
    have hdegC : (SparsePolyZp.toPoly p y).natDegree +
        (SparsePolyZp.toPoly p qc).natDegree = C.natDegree := by
      have hm := Polynomial.natDegree_mul
        (toPoly_ne_zero_of_canonical_nonempty (p := p) y hycan hyne) hqcpoly_ne
      rw [hc_mul] at hm
      exact hm.symm
    have hmeasure : Generated.squarefreeMeasure y + Generated.squarefreeMeasure qc <
        Generated.squarefreeMeasure w + Generated.squarefreeMeasure c := by
      rw [squarefreeMeasure_eq (p := p) y hycan,
        squarefreeMeasure_eq (p := p) qc hqccan,
        squarefreeMeasure_eq (p := p) w hw,
        squarefreeMeasure_eq (p := p) c hc]
      simp only [if_neg hyne, if_neg hqcne]
      have hwarrne := nonempty_of_toPoly_ne_zero (p := p) w hwne
      have hcarrne := nonempty_of_toPoly_ne_zero (p := p) c (by simpa [C] using hcne)
      rw [if_neg hwarrne, if_neg hcarrne]
      dsimp [W, C] at hdegC hwpos ⊢
      omega
    have hdec_eq : Generated.squarefreeMeasure y +
        Generated.squarefreeMeasure (SparsePolyZp.normalization qc) <
        Generated.squarefreeMeasure w + Generated.squarefreeMeasure c := by
      simpa [hnormqc] using hmeasure
    have hiadd : (i + 1).toNat = i.toNat + 1 := by
      rw [UInt64.toNat_add]
      have hi_lt : i.toNat + 1 < 2 ^ 64 := by omega
      exact Nat.mod_eq_of_lt hi_lt
    have hibound' : (i + 1).toNat + (SparsePolyZp.toPoly p y).natDegree +
        (SparsePolyZp.toPoly p qc).natDegree + 2 < 2 ^ 64 := by
      rw [hiadd]
      dsimp [W, C] at hdegC hibound
      omega
    have hdivinst : HasPolyDivmod.polyDivmod (α := SparsePolyZp) = SparsePolyZp.divmod := rfl
    have hnormqw_impl : (SparsePolyZp.divmod w (polynomial_GCD w c)).1.normalization =
        (SparsePolyZp.divmod w (polynomial_GCD w c)).1 := by
      simpa [qw, y] using hnormqw
    have hqcmonic : Monic (SparsePolyZp.toPoly p qc) := by
      rw [Polynomial.Monic.def]
      have hlc := congrArg Polynomial.leadingCoeff hc_mul
      rw [Polynomial.leadingCoeff_mul, hYmonic.leadingCoeff, hcmonic.leadingCoeff,
        one_mul] at hlc
      exact hlc
    have hih_raw := ih hycan (normalization_canonical qc hqccan) hdy
      (by simpa [y, qc, pair_vec_div, hdivinst, hnormqc] using hdqc)
      (by simpa [y, qc, pair_vec_div, hdivinst, normalization_toPoly] using hqcpoly_ne)
      (by simpa [normalization_toPoly] using hqcmonic)
      (by simpa [y, qc, pair_vec_div, hdivinst, normalization_toPoly] using hibound')
    have hdec_impl : Generated.squarefreeMeasure (polynomial_GCD w c) +
        Generated.squarefreeMeasure
          (HasPolyDivmod.polyDivmod c (polynomial_GCD w c)).1.normalization <
        Generated.squarefreeMeasure w + Generated.squarefreeMeasure c := by
      simpa [y, qc, hdivinst, hnormqc] using hmeasure
    rw [Generated._loop___squarefree_Zp_1_ir_def.eq_1, if_pos hcond]
    simp only [pair_vec_div, id_eq]
    rw [dif_pos hdec_impl]
    refine ⟨?_, ?_, ?_, ?_, ?_⟩
    · simpa [y, qc, pair_vec_div, hdivinst, hnormqc] using hih_raw.1
    · simpa [y, qc, pair_vec_div, hdivinst, hnormqc] using hih_raw.2.1
    · simpa [y, qc, pair_vec_div, hdivinst, normalization_toPoly] using hih_raw.2.2.1
    · simpa [y, qc, pair_vec_div, hdivinst, normalization_toPoly] using hih_raw.2.2.2.1
    rw [yunLoop]
    rw [dif_neg hwpos.ne']
    have hYeq0 : SparsePolyZp.toPoly p (polynomial_GCD w c) =
        normalize (EuclideanDomain.gcd (SparsePolyZp.toPoly p w)
          (SparsePolyZp.toPoly p c)) := by
      simpa [y, W, C] using hYeq
    have hqcarg : SparsePolyZp.toPoly p
          (SparsePolyZp.divmod c (polynomial_GCD w c)).1.normalization =
        normalize (SparsePolyZp.toPoly p c /ₘ
          normalize (EuclideanDomain.gcd (SparsePolyZp.toPoly p w)
            (SparsePolyZp.toPoly p c))) := by
      calc
        SparsePolyZp.toPoly p (SparsePolyZp.divmod c (polynomial_GCD w c)).1.normalization =
            SparsePolyZp.toPoly p qc := by
              simpa [qc, y] using (normalization_toPoly (p := p) qc)
        _ = C /ₘ SparsePolyZp.toPoly p y := hqcpoly
        _ = normalize (SparsePolyZp.toPoly p qc) := by
          rw [← hqcpoly]
          exact (normalize_eq_of_monic_local _ hqcmonic).symm
        _ = _ := by rw [hqcpoly, hYeq]
    have himod : (i.toNat + 1) % 2 ^ 64 = i.toNat + 1 := by
      have hi_lt : i.toNat + 1 < 2 ^ 64 := by omega
      exact Nat.mod_eq_of_lt hi_lt
    have himod' : (i.toNat + 1) % 18446744073709551616 = i.toNat + 1 := by
      simpa using himod
    -- The implementation and L2 both choose whether to append the normalized quotient.
    by_cases hzpos : 0 < (SparsePolyZp.toPoly p qw).natDegree
    · have hzcond : ((! (SparsePolyZp.normalization qw).isEmpty) &&
          decide (get_deg (SparsePolyZp.normalization qw) > (0 : Int64))) = true := by
        rw [hnormqw]
        exact (squarefree_loop_cond_iff (p := p) qw hqwcan hdqw).mpr hzpos
      have hzcond_impl :
          ((!Array.isEmpty (HasPolyDivmod.polyDivmod w (polynomial_GCD w c)).1.normalization) &&
            decide (get_deg (HasPolyDivmod.polyDivmod w (polynomial_GCD w c)).1.normalization >
              (0 : Int64))) = true := by
        simpa [qw, y, hdivinst] using hzcond
      have hzcond_raw :
          ((!Array.isEmpty (pair_vec_div SparsePolyZp.empty w (polynomial_GCD w c) w.comp).normalization) &&
            decide (get_deg (pair_vec_div SparsePolyZp.empty w (polynomial_GCD w c) w.comp).normalization >
              (0 : Int64))) = true := by
        simpa [qw, y] using hzcond
      simp [hdivinst, qw, y, hzcond, hnormqw]
      have hzmake : SparsePolyZp.toPoly p (Generated.__upoly_make_monic_ir_def qw).snd =
          normalize (SparsePolyZp.toPoly p qw) := by
        calc
          SparsePolyZp.toPoly p (Generated.__upoly_make_monic_ir_def qw).snd =
              SparsePolyZp.toPoly p (Generated.__upoly_make_monic_ir qw).snd := rfl
          _ = SparsePolyZp.toPoly p (SparsePolyZp.makeMonic qw) :=
            __upoly_make_monic_ir_refines p qw (wellFormed_of_canonical qw hqwcan)
              hqwcan.2.2 h2p
          _ = normalize (SparsePolyZp.toPoly p qw) :=
            makeMonic_toPoly_eq_normalize (p := p) h2p qw hqwcan hqwne
      have hpush : toPolyList (result.push ((Generated.__upoly_make_monic_ir_def qw).snd, i)) p =
          toPolyList result p ++ [(normalize (SparsePolyZp.toPoly p qw), i.toNat)] := by
        simp [toPolyList, hzmake]
      have hih_pos := hih_raw.2.2.2.2
      simp [hzcond_raw, pair_vec_div, hdivinst, hnormqw_impl, hnormqc] at hih_pos
      rw [hih_pos]
      have hzprop : ¬(SparsePolyZp.divmod w (polynomial_GCD w c)).1 = #[] ∧
          0 < get_deg (SparsePolyZp.divmod w (polynomial_GCD w c)).1 := by
        constructor
        · intro he
          apply hqwne
          simp [qw, y, he]
        · exact (get_deg_pos_iff (p := p) qw hqwcan hdqw).mpr hzpos
      have hacc : toPolyList
          (if ¬(SparsePolyZp.divmod w (polynomial_GCD w c)).1 = #[] ∧
              0 < get_deg (SparsePolyZp.divmod w (polynomial_GCD w c)).1 then
            result.push ((Generated.__upoly_make_monic_ir_def
              (SparsePolyZp.divmod w (polynomial_GCD w c)).1).snd, i)
          else result) p =
          toPolyList result p ++ [(normalize (SparsePolyZp.toPoly p w /ₘ
            normalize (EuclideanDomain.gcd (SparsePolyZp.toPoly p w)
              (SparsePolyZp.toPoly p c))), i.toNat)] := by
        rw [if_pos hzprop]
        simpa [qw, y, hqwpoly, hYeq] using hpush
      have hzpos_math : 0 < (normalize (SparsePolyZp.toPoly p w /ₘ
          normalize (EuclideanDomain.gcd (SparsePolyZp.toPoly p w)
            (SparsePolyZp.toPoly p c)))).natDegree := by
        rw [natDegree_normalize_eq_local]
        rw [← hYeq0]
        simpa [W] using (hqwpoly ▸ hzpos)
      simp only [hYeq0, hqcarg, hacc, if_pos hzpos_math]
      rw [himod']
    · have hzcond : ¬((! (SparsePolyZp.normalization qw).isEmpty) &&
          decide (get_deg (SparsePolyZp.normalization qw) > (0 : Int64))) = true := by
        rw [hnormqw]
        exact fun h => hzpos ((squarefree_loop_cond_iff (p := p) qw hqwcan hdqw).mp h)
      have hzcond_impl : ¬
          ((!Array.isEmpty (HasPolyDivmod.polyDivmod w (polynomial_GCD w c)).1.normalization) &&
            decide (get_deg (HasPolyDivmod.polyDivmod w (polynomial_GCD w c)).1.normalization >
              (0 : Int64))) = true := by
        simpa [qw, y, hdivinst] using hzcond
      have hzcond_raw : ¬
          ((!Array.isEmpty (pair_vec_div SparsePolyZp.empty w (polynomial_GCD w c) w.comp).normalization) &&
            decide (get_deg (pair_vec_div SparsePolyZp.empty w (polynomial_GCD w c) w.comp).normalization >
              (0 : Int64))) = true := by
        simpa [qw, y] using hzcond
      simp [hdivinst, qw, y, hzcond, hnormqw]
      have hih_neg := hih_raw.2.2.2.2
      simp [hzcond_raw, pair_vec_div, hdivinst, hnormqw_impl, hnormqc] at hih_neg
      rw [hih_neg]
      have hz0 : (W /ₘ normalize (EuclideanDomain.gcd W C)).natDegree = 0 := by
        simpa [hYeq, hqwpoly] using Nat.eq_zero_of_not_pos hzpos
      have hz0' : (SparsePolyZp.toPoly p w /ₘ
          normalize (EuclideanDomain.gcd (SparsePolyZp.toPoly p w)
            (SparsePolyZp.toPoly p c))).natDegree = 0 := by
        simpa [W, C] using hz0
      have hznpos : ¬0 < (SparsePolyZp.toPoly p w /ₘ
          normalize (EuclideanDomain.gcd (SparsePolyZp.toPoly p w)
            (SparsePolyZp.toPoly p c))).natDegree := by
        rw [hz0']
        omega
      have hzprop : ¬(¬(SparsePolyZp.divmod w (polynomial_GCD w c)).1 = #[] ∧
          0 < get_deg (SparsePolyZp.divmod w (polynomial_GCD w c)).1) := by
        simpa [qw, y, hnormqw_impl] using hzcond
      have hacc : toPolyList
          (if ¬(SparsePolyZp.divmod w (polynomial_GCD w c)).1 = #[] ∧
              0 < get_deg (SparsePolyZp.divmod w (polynomial_GCD w c)).1 then
            result.push ((Generated.__upoly_make_monic_ir_def
              (SparsePolyZp.divmod w (polynomial_GCD w c)).1).snd, i)
          else result) p = toPolyList result p := by
        rw [if_neg hzprop]
      simp only [hYeq0, hqcarg, hacc, natDegree_normalize_eq_local,
        if_neg hznpos]
      rw [himod']
  · intro i w c result hcond
    dsimp only
    intro hnotdec hw hc hdw hdc hcne hcmonic hibound
    dsimp only [motive] at ⊢
    have hwpos := (squarefree_loop_cond_iff (p := p) w hw hdw).mp hcond
    exfalso
    -- The same algebraic degree argument as the recursive case makes the guard true.
    let y := polynomial_GCD w c
    have hwne : SparsePolyZp.toPoly p w ≠ 0 := by intro hz; simp [hz] at hwpos
    have hydata := polynomial_GCD_step_data (p := p) h2p hp2 w c hw hc hwne
    have hycan : CanonicalRep p y := by simpa [y] using hydata.1
    have hyne : ¬y.isEmpty := by simpa [y] using hydata.2.1
    have hYmonic : Monic (SparsePolyZp.toPoly p y) := by simpa [y] using hydata.2.2.2
    let qc := (SparsePolyZp.divmod c y).1
    have hYdvdC : SparsePolyZp.toPoly p y ∣ SparsePolyZp.toPoly p c := by
      rw [hydata.2.2.1]
      exact normalize_dvd_iff.mpr (EuclideanDomain.gcd_dvd_right _ _)
    have hqccan := divmod_fst_canonical (p := p) c y hc hycan hyne h2p
    have hqcpoly := divmod_fst_toPoly_eq_divByMonic (p := p) h2p hp2 c y hc hycan
      hyne hYmonic hYdvdC
    have hqcpoly_ne : SparsePolyZp.toPoly p qc ≠ 0 := by
      simpa [qc, hqcpoly] using divByMonic_ne_zero_local (SparsePolyZp.toPoly p c)
        (SparsePolyZp.toPoly p y) hYmonic hYdvdC hcne
    have hqcne := nonempty_of_toPoly_ne_zero (p := p) qc hqcpoly_ne
    have hc_mul : SparsePolyZp.toPoly p y * SparsePolyZp.toPoly p qc = SparsePolyZp.toPoly p c := by
      have hid := modByMonic_add_div (SparsePolyZp.toPoly p c) (SparsePolyZp.toPoly p y)
      rw [(modByMonic_eq_zero_iff_dvd hYmonic).mpr hYdvdC, zero_add] at hid
      simpa [qc, hqcpoly] using hid
    have hdegC := Polynomial.natDegree_mul
      (toPoly_ne_zero_of_canonical_nonempty (p := p) y hycan hyne) hqcpoly_ne
    rw [hc_mul] at hdegC
    have hmeasure : Generated.squarefreeMeasure y + Generated.squarefreeMeasure qc <
        Generated.squarefreeMeasure w + Generated.squarefreeMeasure c := by
      rw [squarefreeMeasure_eq (p := p) y hycan, squarefreeMeasure_eq (p := p) qc hqccan,
        squarefreeMeasure_eq (p := p) w hw, squarefreeMeasure_eq (p := p) c hc]
      have hwarrne := nonempty_of_toPoly_ne_zero (p := p) w hwne
      have hcarrne := nonempty_of_toPoly_ne_zero (p := p) c hcne
      simp only [if_neg hyne, if_neg hqcne, if_neg hwarrne, if_neg hcarrne]
      omega
    apply hnotdec
    have hdivinst : HasPolyDivmod.polyDivmod (α := SparsePolyZp) = SparsePolyZp.divmod := rfl
    simpa [y, qc, pair_vec_div, hdivinst,
      normalization_eq_of_nonZeroB qc hqccan.2.1] using hmeasure
  · intro i w c result hcond hw hc hdw hdc hcne hcmonic hibound
    dsimp only [motive]
    have hwzero : (SparsePolyZp.toPoly p w).natDegree = 0 := by
      by_contra hn
      apply hcond
      exact (squarefree_loop_cond_iff (p := p) w hw hdw).mpr (Nat.pos_of_ne_zero hn)
    rw [Generated._loop___squarefree_Zp_1_ir_def.eq_1, if_neg hcond]
    unfold yunLoop
    refine ⟨hc, hdc, hcne, hcmonic, ?_⟩
    rw [dif_pos hwzero]
/-- __upoly_make_monic_ir 保持 AllReduced（需 hp_size 保证 UInt64 roundtrip）。 -/
lemma upoly_make_monic_allReduced (f : SparsePolyZp) (hred : SparsePolyZp.AllReduced p f.toList)
    (hp_size : 2 * p ≤ UInt64.size) :
    SparsePolyZp.AllReduced p (Generated.__upoly_make_monic_ir f).snd.toList := by
  intro x hx
  simp [Generated.__upoly_make_monic_ir, Generated.__upoly_make_monic_ir_def] at hx
  split_ifs at hx with h
  · simp at hx; exact hred x hx.val
  · simp at hx
    have hloop : (Generated._loop___upoly_make_monic_0_ir 0 f (Zp.inv ((SparsePolyZp.front! f).snd))).snd.toList =
        (f.toList).map (fun (m, x) => (m, x * (Zp.inv ((SparsePolyZp.front! f).snd)))) := by
      have := loop_result_toList f (Zp.inv ((SparsePolyZp.front! f).snd)) 0 (Nat.zero_le _)
      simpa [List.take, List.drop] using this
    have hx_val : x ∈ (Generated._loop___upoly_make_monic_0_ir 0 f (Zp.inv ((SparsePolyZp.front! f).snd))).snd.toList :=
      (show x ∈ ((Generated._loop___upoly_make_monic_0_ir 0 f (Zp.inv ((SparsePolyZp.front! f).snd))).snd : SparsePolyZp) from hx).val
    have hx_mem : x ∈ (f.toList).map (fun (m, x) => (m, x * (Zp.inv ((SparsePolyZp.front! f).snd)))) := by
      rw [← hloop]; exact hx_val
    rcases List.mem_map.mp hx_mem with ⟨y, hy, rfl⟩
    have hy_red : Zp.Reduced p y.2 := hred y hy
    have hp_pos : 0 < p := Nat.Prime.pos (hp.out)
    have hp_le_size : p ≤ UInt64.size := by
      have : 2 * p ≤ UInt64.size := hp_size
      have hp_pos : p > 0 := Nat.Prime.pos (hp.out)
      have : p ≤ 2 * p := by nlinarith
      exact Nat.le_trans this hp_size
    have h_mod_lt_uint64 : (y.2.val.toNat * (Zp.inv ((SparsePolyZp.front! f).snd)).val.toNat) % y.2.prime.toNat < UInt64.size := by
      have h_mod_lt_p : (y.2.val.toNat * (Zp.inv ((SparsePolyZp.front! f).snd)).val.toNat) % y.2.prime.toNat < y.2.prime.toNat :=
        Nat.mod_lt _ (by
          have : y.2.prime.toNat = p := hy_red.1
          rw [this]
          exact Nat.Prime.pos (hp.out))
      have h_p_eq : y.2.prime.toNat = p := hy_red.1
      have h_mod_lt_p' : (y.2.val.toNat * (Zp.inv ((SparsePolyZp.front! f).snd)).val.toNat) % p < p := by
        simpa [h_p_eq] using h_mod_lt_p
      have h_temp : (y.2.val.toNat * (Zp.inv ((SparsePolyZp.front! f).snd)).val.toNat) % p < UInt64.size :=
        Nat.lt_of_lt_of_le h_mod_lt_p' hp_le_size
      simpa [hy_red.1] using h_temp
    have h_val : (y.2 * (Zp.inv ((SparsePolyZp.front! f).snd))).val.toNat =
        (y.2.val.toNat * (Zp.inv ((SparsePolyZp.front! f).snd)).val.toNat) % y.2.prime.toNat := by
      calc
        (y.2 * (Zp.inv ((SparsePolyZp.front! f).snd))).val.toNat
            = ((y.2.val.toNat * (Zp.inv ((SparsePolyZp.front! f).snd)).val.toNat) % y.2.prime.toNat).toUInt64.toNat := rfl
        _ = (y.2.val.toNat * (Zp.inv ((SparsePolyZp.front! f).snd)).val.toNat) % y.2.prime.toNat := by
          simp [UInt64.toNat_ofNat, h_mod_lt_uint64]
    have : (y.2 * (Zp.inv ((SparsePolyZp.front! f).snd))).prime.toNat = p := by
      calc
        (y.2 * (Zp.inv ((SparsePolyZp.front! f).snd))).prime.toNat = y.2.prime.toNat := rfl
        _ = p := hy_red.1
    have hval_lt_p : (y.2 * (Zp.inv ((SparsePolyZp.front! f).snd))).val.toNat < p := by
      rw [h_val, hy_red.1]
      exact Nat.mod_lt _ hp_pos
    exact ⟨this, hval_lt_p⟩

/-- __upoly_make_monic_ir 保持 WellFormed。 -/
lemma upoly_make_monic_wellFormed (f : SparsePolyZp) (hwf : SparsePolyZp.WellFormed p f) :
    SparsePolyZp.WellFormed p (Generated.__upoly_make_monic_ir f).snd := by
  intro x hx
  simp [Generated.__upoly_make_monic_ir, Generated.__upoly_make_monic_ir_def] at hx
  split_ifs at hx with h
  · simp at hx; exact hwf x hx
  · simp at hx
    have hloop : (Generated._loop___upoly_make_monic_0_ir 0 f (Zp.inv ((SparsePolyZp.front! f).snd))).snd.toList =
        (f.toList).map (fun (m, x) => (m, x * (Zp.inv ((SparsePolyZp.front! f).snd)))) := by
      have := loop_result_toList f (Zp.inv ((SparsePolyZp.front! f).snd)) 0 (Nat.zero_le _)
      simpa [List.take, List.drop] using this
    have hx_val : x ∈ (Generated._loop___upoly_make_monic_0_ir 0 f (Zp.inv ((SparsePolyZp.front! f).snd))).snd.toList :=
      (show x ∈ ((Generated._loop___upoly_make_monic_0_ir 0 f (Zp.inv ((SparsePolyZp.front! f).snd))).snd : SparsePolyZp) from hx).val
    have hx_mem : x ∈ (f.toList).map (fun (m, x) => (m, x * (Zp.inv ((SparsePolyZp.front! f).snd)))) := by
      rw [← hloop]; exact hx_val
    rcases List.mem_map.mp hx_mem with ⟨y, hy, rfl⟩
    exact hwf y (Array.Mem.mk hy)

/-- `__upoly_make_monic_ir` 只改变系数，保持严格降序的稀疏表示。 -/
lemma upoly_make_monic_sorted (f : SparsePolyZp) (hs : SparsePolyZp.Sorted f) :
    SparsePolyZp.Sorted (Generated.__upoly_make_monic_ir f).snd := by
  unfold SparsePolyZp.Sorted
  simp [Generated.__upoly_make_monic_ir, Generated.__upoly_make_monic_ir_def]
  split_ifs with h
  · exact hs
  · have hloop :
        (Generated._loop___upoly_make_monic_0_ir 0 f
          (Zp.inv ((SparsePolyZp.front! f).snd))).snd.toList =
        f.toList.map (fun (m, x) =>
          (m, x * Zp.inv ((SparsePolyZp.front! f).snd))) := by
        have hres := loop_result_toList f
          (Zp.inv ((SparsePolyZp.front! f).snd)) 0 (Nat.zero_le _)
        simpa [List.take, List.drop] using hres
    rw [hloop]
    exact (SparsePolyZp.sortedListB_map_fst _
      (by rintro ⟨m, c⟩; rfl) f.toList).trans hs

/-- `__upoly_make_monic_ir` 不改变项的次数，因而保持任意次数上界。 -/
lemma upoly_make_monic_deg_bound_of (B : Nat) (f : SparsePolyZp)
    (h_deg_bound : ∀ x ∈ f.toList, x.1.deg < B) :
    ∀ x ∈ (Generated.__upoly_make_monic_ir f).snd.toList, x.1.deg < B := by
  intro x hx
  simp [Generated.__upoly_make_monic_ir, Generated.__upoly_make_monic_ir_def] at hx ⊢
  split_ifs at hx with h
  · simp at hx; exact h_deg_bound x hx.val
  · simp at hx
    have hloop : (Generated._loop___upoly_make_monic_0_ir 0 f (Zp.inv ((SparsePolyZp.front! f).snd))).snd.toList =
        (f.toList).map (fun (m, x) => (m, x * (Zp.inv ((SparsePolyZp.front! f).snd)))) := by
      have := loop_result_toList f (Zp.inv ((SparsePolyZp.front! f).snd)) 0 (Nat.zero_le _)
      simpa [List.take, List.drop] using this
    have hx_val : x ∈ (Generated._loop___upoly_make_monic_0_ir 0 f (Zp.inv ((SparsePolyZp.front! f).snd))).snd.toList :=
      (show x ∈ ((Generated._loop___upoly_make_monic_0_ir 0 f (Zp.inv ((SparsePolyZp.front! f).snd))).snd : SparsePolyZp) from hx).val
    have hx_mem : x ∈ (f.toList).map (fun (m, x) => (m, x * (Zp.inv ((SparsePolyZp.front! f).snd)))) := by
      rw [← hloop]; exact hx_val
    rcases List.mem_map.mp hx_mem with ⟨y, hy, rfl⟩
    exact h_deg_bound y hy

/-- `__upoly_make_monic_ir` 保持 UInt64 次数上界。 -/
lemma upoly_make_monic_deg_bound (f : SparsePolyZp)
    (h_deg_bound : ∀ x ∈ f.toList, x.1.deg < 2 ^ 64) :
    ∀ x ∈ (Generated.__upoly_make_monic_ir f).snd.toList, x.1.deg < 2 ^ 64 :=
  upoly_make_monic_deg_bound_of (2 ^ 64) f h_deg_bound

lemma upoly_make_monic_p_deg_bound (f : SparsePolyZp)
    (h : ∀ x ∈ f.toList, p * x.1.deg < 2 ^ 64) :
    ∀ x ∈ (Generated.__upoly_make_monic_ir f).snd.toList,
      p * x.1.deg < 2 ^ 64 := by
  intro x hx
  simp [Generated.__upoly_make_monic_ir, Generated.__upoly_make_monic_ir_def] at hx ⊢
  split_ifs at hx with hc
  · exact h x hx.val
  · simp at hx
    have hloop :
        (Generated._loop___upoly_make_monic_0_ir 0 f
          (Zp.inv ((SparsePolyZp.front! f).snd))).snd.toList =
        f.toList.map (fun (m, y) =>
          (m, y * Zp.inv ((SparsePolyZp.front! f).snd))) := by
      have hr := loop_result_toList f (Zp.inv ((SparsePolyZp.front! f).snd))
        0 (Nat.zero_le _)
      simpa [List.take, List.drop] using hr
    have hxlist : x ∈
        (Generated._loop___upoly_make_monic_0_ir 0 f
          (Zp.inv ((SparsePolyZp.front! f).snd))).snd.toList := hx.val
    rw [hloop] at hxlist
    rcases List.mem_map.mp hxlist with ⟨y, hy, rfl⟩
    exact h y hy

lemma upoly_make_monic_val_nonzero (f : SparsePolyZp)
    (hred : SparsePolyZp.AllReduced p f.toList)
    (hval : ∀ x ∈ f.toList, x.2.val.toNat ≠ 0)
    (hfne : ¬f.isEmpty) (h2p : 2 * p ≤ UInt64.size) :
    ∀ x ∈ (Generated.__upoly_make_monic_ir f).snd.toList,
      x.2.val.toNat ≠ 0 := by
  have hfrontmem : SparsePolyZp.front! f ∈ f.toList := mem_getFirst_toList f (by
    rw [Array.isEmpty_iff_size_eq_zero] at hfne
    omega)
  have hredlc := hred (SparsePolyZp.front! f) hfrontmem
  have hvallc := hval (SparsePolyZp.front! f) hfrontmem
  have hredinv : Zp.Reduced p (SparsePolyZp.front! f).snd.inv := by
    have hprime := hredlc.1
    constructor
    · simpa [Zp.inv] using hprime
    · calc
        (SparsePolyZp.front! f).snd.inv.val.toNat =
            (Zp.modInv (SparsePolyZp.front! f).snd.val
              (SparsePolyZp.front! f).snd.prime).toNat := rfl
        _ < (SparsePolyZp.front! f).snd.prime.toNat :=
          modInv_val_lt _ _ (by rw [hprime]; exact hp.out.one_lt)
            (by rw [hprime]; omega)
        _ = p := hprime
  have hinvval : (SparsePolyZp.front! f).snd.inv.val.toNat ≠ 0 :=
    Zp_inv_val_nonzero (SparsePolyZp.front! f).snd hredlc hvallc h2p
  intro x hx
  simp [Generated.__upoly_make_monic_ir, Generated.__upoly_make_monic_ir_def] at hx
  split_ifs at hx with h
  · exact hval x hx.val
  · simp at hx
    let inv := (SparsePolyZp.front! f).snd.inv
    have hloop : (Generated._loop___upoly_make_monic_0_ir 0 f inv).snd.toList =
        f.toList.map (fun (m, y) => (m, y * inv)) := by
      have hr := loop_result_toList f inv 0 (Nat.zero_le _)
      simpa [List.take, List.drop] using hr
    have hxlist : x ∈ (Generated._loop___upoly_make_monic_0_ir 0 f inv).snd.toList := hx.val
    rw [hloop] at hxlist
    rcases List.mem_map.mp hxlist with ⟨y, hy, rfl⟩
    exact Zp_mul_val_nonzero (hred y hy) (by simpa [inv] using hredinv)
      (hval y hy) (by simpa [inv] using hinvval) h2p

/-- _loop___extract_pth_root_0_ir_def 的 toList = acc.toList ++ (drop idx f.toList).map Φ。 -/
private lemma loop_extract_toList (f : SparsePolyZp) (acc : SparsePolyZp) (idx : Nat) (p_1 : UInt64) (hidx : idx ≤ Array.size f) :
    (Generated._loop___extract_pth_root_0_ir_def idx acc f p_1).snd.toList =
    acc.toList ++ (List.drop idx (f.toList)).map (λ term : UMonomial × Zp =>
      (UMonomial.mk ((term.1.deg.toUInt64 / p_1).toInt64), term.2)) := by
  have h_wf : WellFounded (λ (a b : SparsePolyZp × Nat) => Array.size a.1 - a.2 < Array.size b.1 - b.2) :=
    (measure (λ (p : SparsePolyZp × Nat) => Array.size p.1 - p.2)).wf
  refine h_wf.induction (f, idx) (C := λ p => ∀ (acc : SparsePolyZp), p.2 ≤ Array.size p.1 →
    (Generated._loop___extract_pth_root_0_ir_def p.2 acc p.1 p_1).snd.toList =
    acc.toList ++ (List.drop p.2 (p.1.toList)).map (λ term => (UMonomial.mk ((term.1.deg.toUInt64 / p_1).toInt64), term.2))) ?_ acc hidx
  intro p ih acc hidx
  rw [Generated._loop___extract_pth_root_0_ir_def.eq_1]
  by_cases h : p.2 < Array.size p.1
  · simp [h]
    have hx_len : p.2 < p.1.toList.length := by simpa using h
    have hx_drop : List.drop p.2 (p.1.toList) = (p.1[p.2]) :: List.drop (p.2 + 1) (p.1.toList) := by
      have htemp := drop_eq_get_cons (p.1.toList) p.2 hx_len
      calc
        List.drop p.2 (p.1.toList)
            = (p.1.toList).get ⟨p.2, hx_len⟩ :: List.drop (p.2 + 1) (p.1.toList) := htemp
        _ = p.1[p.2] :: List.drop (p.2 + 1) (p.1.toList) := by simp
    have h_measure : Array.size p.1 - (p.2 + 1) < Array.size p.1 - p.2 := by omega
    have h_idx_succ : p.2 + 1 ≤ Array.size p.1 := by omega
    have h_ih := ih (p.1, p.2 + 1) h_measure
      (Array.push acc (Prod.mk (UMonomial.mk (((p.1[p.2]).1.deg.toUInt64 / p_1).toInt64)) (p.1[p.2]).2))
      h_idx_succ
    have h_target : (Generated._loop___extract_pth_root_0_ir_def (p.2 + 1)
        (Array.push acc (Prod.mk (UMonomial.mk (((p.1[p.2]).1.deg.toUInt64 / p_1).toInt64)) (p.1[p.2]).2))
        p.1 p_1).snd.toList = acc.toList ++ List.drop p.2 ((p.1.toList).map (λ term => (UMonomial.mk ((term.1.deg.toUInt64 / p_1).toInt64), term.2))) :=
      calc
        (Generated._loop___extract_pth_root_0_ir_def (p.2 + 1)
            (Array.push acc (Prod.mk (UMonomial.mk (((p.1[p.2]).1.deg.toUInt64 / p_1).toInt64)) (p.1[p.2]).2))
            p.1 p_1).snd.toList
            = (Array.push acc (Prod.mk (UMonomial.mk (((p.1[p.2]).1.deg.toUInt64 / p_1).toInt64)) (p.1[p.2]).2)).toList ++
              (List.drop (p.2 + 1) (p.1.toList)).map (λ term => (UMonomial.mk ((term.1.deg.toUInt64 / p_1).toInt64), term.2)) :=
          h_ih
        _ = (acc.toList ++ [Prod.mk (UMonomial.mk (((p.1[p.2]).1.deg.toUInt64 / p_1).toInt64)) (p.1[p.2]).2]) ++
              (List.drop (p.2 + 1) (p.1.toList)).map (λ term => (UMonomial.mk ((term.1.deg.toUInt64 / p_1).toInt64), term.2)) := by simp
        _ = acc.toList ++ ([Prod.mk (UMonomial.mk (((p.1[p.2]).1.deg.toUInt64 / p_1).toInt64)) (p.1[p.2]).2] ++
              (List.drop (p.2 + 1) (p.1.toList)).map (λ term => (UMonomial.mk ((term.1.deg.toUInt64 / p_1).toInt64), term.2))) := by
          simp [List.append_assoc]
        _ = acc.toList ++ (((p.1[p.2]) :: List.drop (p.2 + 1) (p.1.toList)).map (λ term => (UMonomial.mk ((term.1.deg.toUInt64 / p_1).toInt64), term.2))) := by simp
        _ = acc.toList ++ (List.drop p.2 (p.1.toList)).map (λ term => (UMonomial.mk ((term.1.deg.toUInt64 / p_1).toInt64), term.2)) := by
          have h_map_eq : ((p.1[p.2]) :: List.drop (p.2 + 1) (p.1.toList)).map (λ term => (UMonomial.mk ((term.1.deg.toUInt64 / p_1).toInt64), term.2))
                        = (List.drop p.2 (p.1.toList)).map (λ term => (UMonomial.mk ((term.1.deg.toUInt64 / p_1).toInt64), term.2)) :=
            congrArg (·.map (λ term => (UMonomial.mk ((term.1.deg.toUInt64 / p_1).toInt64), term.2))) hx_drop.symm
          rw [h_map_eq]
        _ = acc.toList ++ List.drop p.2 ((p.1.toList).map (λ term => (UMonomial.mk ((term.1.deg.toUInt64 / p_1).toInt64), term.2))) := by
          simp [List.map_drop]
    exact h_target
  · -- i ≥ size, loop returns acc directly, List.drop is empty
    have h_drop_empty : List.drop p.2 (p.1.toList) = [] :=
      List.drop_eq_nil_of_le (by
        have : p.1.toList.length = Array.size p.1 := by simp
        rw [this]
        omega)
    simp [Generated._loop___extract_pth_root_0_ir_def, h, h_drop_empty]

/-- __extract_pth_root_ir 保持 WellFormed。 -/
lemma extract_pth_root_wellFormed (g : SparsePolyZp) (hwf : SparsePolyZp.WellFormed p g) :
    SparsePolyZp.WellFormed p (Generated.__extract_pth_root_ir g) := by
  intro x hx
  have h_loop : (Generated.__extract_pth_root_ir g).toList =
      g.toList.map (λ term : UMonomial × Zp => (UMonomial.mk ((term.1.deg.toUInt64 / (SparsePolyZp.front! g).snd.prime).toInt64), term.2)) := by
    unfold Generated.__extract_pth_root_ir Generated.__extract_pth_root_ir_def
    have h_loop' := loop_extract_toList g SparsePolyZp.empty 0 (SparsePolyZp.front! g).snd.prime (by simp)
    simpa [SparsePolyZp.empty, List.drop] using h_loop'
  have hx_list : x ∈ (Generated.__extract_pth_root_ir g).toList := by
    simpa using hx
  -- Apply List.mem_map directly via `apply`
  have hx_map : x ∈ (g.toList.map (λ term : UMonomial × Zp => (UMonomial.mk ((term.1.deg.toUInt64 / (SparsePolyZp.front! g).snd.prime).toInt64), term.2))) := by
    -- Rewrite using h_loop in the membership
    rw [← h_loop]
    exact hx_list
  -- Now use mem_map to decompose
  simp only [List.mem_map] at hx_map
  obtain ⟨y, hy, h_eq⟩ := hx_map
  rw [← h_eq]
  have hy_arr : y ∈ g := by simpa using hy
  exact hwf y hy_arr

/-- __extract_pth_root_ir 保持 AllReduced。 -/
lemma extract_pth_root_allReduced (g : SparsePolyZp) (hred : SparsePolyZp.AllReduced p g.toList) :
    SparsePolyZp.AllReduced p (Generated.__extract_pth_root_ir g).toList := by
  have h_loop : (Generated.__extract_pth_root_ir g).toList =
      g.toList.map (λ term : UMonomial × Zp => (UMonomial.mk ((term.1.deg.toUInt64 / (SparsePolyZp.front! g).snd.prime).toInt64), term.2)) := by
    unfold Generated.__extract_pth_root_ir Generated.__extract_pth_root_ir_def
    have h_loop' := loop_extract_toList g SparsePolyZp.empty 0 (SparsePolyZp.front! g).snd.prime (by simp)
    simpa [SparsePolyZp.empty, List.drop] using h_loop'
  rw [h_loop]
  intro x hx
  rcases List.mem_map.mp hx with ⟨y, hy, rfl⟩
  exact hred y hy

/-- __extract_pth_root_ir 保持系数的 val ≠ 0 性质（无零系数项）。 -/
lemma extract_pth_root_val_nonzero (g : SparsePolyZp)
    (h_val_nonzero : ∀ x ∈ g.toList, x.snd.val.toNat ≠ 0) :
    ∀ x ∈ (Generated.__extract_pth_root_ir g).toList, x.snd.val.toNat ≠ 0 := by
  have h_loop : (Generated.__extract_pth_root_ir g).toList =
      g.toList.map (λ term : UMonomial × Zp => (UMonomial.mk ((term.1.deg.toUInt64 / (SparsePolyZp.front! g).snd.prime).toInt64), term.2)) := by
    unfold Generated.__extract_pth_root_ir Generated.__extract_pth_root_ir_def
    have h_loop' := loop_extract_toList g SparsePolyZp.empty 0 (SparsePolyZp.front! g).snd.prime (by simp)
    simpa [SparsePolyZp.empty, List.drop] using h_loop'
  rw [h_loop]
  intro x hx
  rcases List.mem_map.mp hx with ⟨y, hy, rfl⟩
  -- The new term's .snd is y.2, which is unchanged by extract_pth_root_ir
  simpa using h_val_nonzero y hy

/-- __extract_pth_root_ir 保持 degree bound。 -/
lemma extract_pth_root_deg_bound (g : SparsePolyZp) (h_deg_bound : ∀ x ∈ g.toList, x.1.deg < 2 ^ 64) :
    ∀ x ∈ (Generated.__extract_pth_root_ir g).toList, x.1.deg < 2 ^ 64 := by
  have h_loop : (Generated.__extract_pth_root_ir g).toList =
      g.toList.map (λ term : UMonomial × Zp => (UMonomial.mk ((term.1.deg.toUInt64 / (SparsePolyZp.front! g).snd.prime).toInt64), term.2)) := by
    unfold Generated.__extract_pth_root_ir Generated.__extract_pth_root_ir_def
    have h_loop' := loop_extract_toList g SparsePolyZp.empty 0 (SparsePolyZp.front! g).snd.prime (by simp)
    simpa [SparsePolyZp.empty, List.drop] using h_loop'
  rw [h_loop]
  intro x hx
  rcases List.mem_map.mp hx with ⟨y, hy, rfl⟩
  have hdeg : y.1.deg < 2 ^ 64 := h_deg_bound y hy
  have hU64_div (a b : UInt64) : (a / b).toNat = a.toNat / b.toNat := by
    calc
      (a / b).toNat = (a.toBitVec.udiv b.toBitVec).toNat := rfl
      _ = a.toBitVec.toNat / b.toBitVec.toNat := by
        unfold BitVec.udiv; simp
      _ = a.toNat / b.toNat := by simp
  have hdiv_le : (UMonomial.mk ((y.1.deg.toUInt64 / (SparsePolyZp.front! g).snd.prime).toInt64)).deg ≤ y.1.deg := by
    -- UMonomial.mk with Int64 argument has deg = (toInt64).toNatClampNeg. When non-negative and < 2^63,
    -- this equals the UInt64.toNat of the original value.
    -- For y.1.deg < 2^64 and division by prime ≥ 2, this holds. Needs a dedicated lemma.
    have h_deg_lt : (UMonomial.mk ((y.1.deg.toUInt64 / (SparsePolyZp.front! g).snd.prime).toInt64)).deg ≤
        (y.1.deg.toUInt64 / (SparsePolyZp.front! g).snd.prime).toNat := by
      calc
        (UMonomial.mk ((y.1.deg.toUInt64 / (SparsePolyZp.front! g).snd.prime).toInt64)).deg
            = ((y.1.deg.toUInt64 / (SparsePolyZp.front! g).snd.prime).toInt64).toNatClampNeg := rfl
        _ ≤ (y.1.deg.toUInt64 / (SparsePolyZp.front! g).snd.prime).toNat :=
          UInt64_toInt64_toNatClampNeg_le_toNat _
    have hval_le : (y.1.deg.toUInt64 / (SparsePolyZp.front! g).snd.prime).toNat ≤ y.1.deg := by
      calc
        (y.1.deg.toUInt64 / (SparsePolyZp.front! g).snd.prime).toNat
            = (y.1.deg.toUInt64).toNat / ((SparsePolyZp.front! g).snd.prime).toNat := by rw [hU64_div]
        _ = y.1.deg / ((SparsePolyZp.front! g).snd.prime).toNat := by
          have h_degU64_toNat : (y.1.deg.toUInt64).toNat = y.1.deg := by
            have : y.1.deg % 2 ^ 64 = y.1.deg := Nat.mod_eq_of_lt hdeg
            simpa [UInt64.toNat_ofNat] using this
          rw [h_degU64_toNat]
        _ ≤ y.1.deg := Nat.div_le_self _ _
    omega
    -- Original calc kept below for reference
    /- calc
      (UMonomial.mk ((y.1.deg.toUInt64 / (SparsePolyZp.front! g).snd.prime).toInt64)).deg
          = ((y.1.deg.toUInt64 / (SparsePolyZp.front! g).snd.prime).toInt64 : ℕ) := rfl
      _ = ((y.1.deg.toUInt64 / (SparsePolyZp.front! g).snd.prime).toNat) := by simp
      _ = (y.1.deg.toUInt64).toNat / ((SparsePolyZp.front! g).snd.prime).toNat := by rw [hU64_div]
      _ = y.1.deg / ((SparsePolyZp.front! g).snd.prime).toNat := by
        have h_degU64_toNat : (y.1.deg.toUInt64).toNat = y.1.deg := by
          have : y.1.deg % 2 ^ 64 = y.1.deg := Nat.mod_eq_of_lt hdeg
          simpa [UInt64.toNat_ofNat] using this
        rw [h_degU64_toNat]
      _ ≤ y.1.deg := Nat.div_le_self _ _
    -/
  have h_deg_bound' : (UMonomial.mk ((y.1.deg.toUInt64 / (SparsePolyZp.front! g).snd.prime).toInt64)).deg < 2 ^ 64 := by
    calc
      (UMonomial.mk ((y.1.deg.toUInt64 / (SparsePolyZp.front! g).snd.prime).toInt64)).deg
          ≤ y.1.deg := hdiv_le
      _ < 2 ^ 64 := hdeg
  exact h_deg_bound'

/-- `__extract_pth_root_ir` 保持 `get_deg` 的有符号 64 位安全范围。 -/
lemma extract_pth_root_signed_deg_bound (g : SparsePolyZp)
    (h_deg_bound : ∀ x ∈ g.toList, x.1.deg < 2 ^ 63) :
    ∀ x ∈ (Generated.__extract_pth_root_ir g).toList, x.1.deg < 2 ^ 63 := by
  have h_loop : (Generated.__extract_pth_root_ir g).toList =
      g.toList.map (λ term : UMonomial × Zp =>
        (UMonomial.mk ((term.1.deg.toUInt64 /
          (SparsePolyZp.front! g).snd.prime).toInt64), term.2)) := by
    unfold Generated.__extract_pth_root_ir Generated.__extract_pth_root_ir_def
    have h := loop_extract_toList g SparsePolyZp.empty 0
      (SparsePolyZp.front! g).snd.prime (by simp)
    simpa [SparsePolyZp.empty, List.drop] using h
  rw [h_loop]
  intro x hx
  rcases List.mem_map.mp hx with ⟨y, hy, rfl⟩
  have hy63 := h_deg_bound y hy
  have hy64 : y.1.deg < 2 ^ 64 := lt_trans hy63 (by norm_num)
  have hU64_div (a b : UInt64) : (a / b).toNat = a.toNat / b.toNat := by
    calc
      (a / b).toNat = (a.toBitVec.udiv b.toBitVec).toNat := rfl
      _ = a.toBitVec.toNat / b.toBitVec.toNat := by
        unfold BitVec.udiv
        simp
      _ = a.toNat / b.toNat := by simp
  have hdiv_le :
      (UMonomial.mk ((y.1.deg.toUInt64 /
        (SparsePolyZp.front! g).snd.prime).toInt64)).deg ≤ y.1.deg := by
    have hclamp :
        (UMonomial.mk ((y.1.deg.toUInt64 /
          (SparsePolyZp.front! g).snd.prime).toInt64)).deg ≤
          (y.1.deg.toUInt64 / (SparsePolyZp.front! g).snd.prime).toNat := by
      calc
        (UMonomial.mk ((y.1.deg.toUInt64 /
            (SparsePolyZp.front! g).snd.prime).toInt64)).deg =
            ((y.1.deg.toUInt64 /
              (SparsePolyZp.front! g).snd.prime).toInt64).toNatClampNeg := rfl
        _ ≤ (y.1.deg.toUInt64 /
              (SparsePolyZp.front! g).snd.prime).toNat :=
          UInt64_toInt64_toNatClampNeg_le_toNat _
    calc
      _ ≤ (y.1.deg.toUInt64 /
          (SparsePolyZp.front! g).snd.prime).toNat := hclamp
      _ = (y.1.deg.toUInt64).toNat /
          (SparsePolyZp.front! g).snd.prime.toNat := by
        rw [hU64_div]
      _ = y.1.deg / (SparsePolyZp.front! g).snd.prime.toNat := by
        have hround : (y.1.deg.toUInt64).toNat = y.1.deg := by
          have hm : y.1.deg % 2 ^ 64 = y.1.deg := Nat.mod_eq_of_lt hy64
          simpa [UInt64.toNat_ofNat] using hm
        rw [hround]
      _ ≤ y.1.deg := Nat.div_le_self _ _
  exact lt_of_le_of_lt hdiv_le hy63

private theorem UInt64_toNat_div (a b : UInt64) : (a / b).toNat = a.toNat / b.toNat := by
  calc
    (a / b).toNat = (a.toBitVec.udiv b.toBitVec).toNat := rfl
    _ = a.toBitVec.toNat / b.toBitVec.toNat := by
      unfold BitVec.udiv; simp
    _ = a.toNat / b.toNat := by simp

/-- __extract_pth_root_ir 保持 no_overflow，且因 degree 被 p 除，还有 p * deg < 2^64。 -/
lemma extract_pth_root_no_overflow (g : SparsePolyZp) (hwf : SparsePolyZp.WellFormed p g) (h_no_overflow : ∀ x ∈ g.toList, x.2.val.toNat * x.1.deg < 2 ^ 64)
    (h_deg_bound : ∀ x ∈ g.toList, x.1.deg < 2 ^ 64) :
    (∀ x ∈ (Generated.__extract_pth_root_ir g).toList, x.2.val.toNat * x.1.deg < 2 ^ 64) ∧
    (∀ x ∈ (Generated.__extract_pth_root_ir g).toList, p * x.1.deg < 2 ^ 64) := by
  have h_loop : (Generated.__extract_pth_root_ir g).toList =
      g.toList.map (λ term : UMonomial × Zp => (UMonomial.mk ((term.1.deg.toUInt64 / (SparsePolyZp.front! g).snd.prime).toInt64), term.2)) := by
    unfold Generated.__extract_pth_root_ir Generated.__extract_pth_root_ir_def
    have h_loop' := loop_extract_toList g SparsePolyZp.empty 0 (SparsePolyZp.front! g).snd.prime (by simp)
    simpa [SparsePolyZp.empty, List.drop] using h_loop'
  constructor
  · intro x hx
    rw [h_loop] at hx
    rcases List.mem_map.mp hx with ⟨y, hy, rfl⟩
    have h_ov : y.2.val.toNat * y.1.deg < 2 ^ 64 := h_no_overflow y hy
    have hdiv_le : (UMonomial.mk ((y.1.deg.toUInt64 / (SparsePolyZp.front! g).snd.prime).toInt64)).deg ≤ y.1.deg := by
      have hdeg : y.1.deg < 2 ^ 64 := h_deg_bound y hy
      have hp_ge_2 : 2 ≤ p := Nat.Prime.two_le hp.out
      have hbound : (y.1.deg.toUInt64 / (SparsePolyZp.front! g).snd.prime).toNat < 2 ^ 63 := by
        rw [UInt64_toNat_div]
        have h_degU64_toNat : (y.1.deg.toUInt64).toNat = y.1.deg := by
          have : y.1.deg % 2 ^ 64 = y.1.deg := Nat.mod_eq_of_lt hdeg
          simpa [UInt64.toNat_ofNat] using this
        rw [h_degU64_toNat]
        have hp_div : ((SparsePolyZp.front! g).snd.prime).toNat = p := by
          have hpos : 0 < g.size := by
            have : (g.toList).length > 0 := List.length_pos_of_mem hy
            have hsize : g.size = (g.toList).length := by simp
            omega
          have hmem_front : SparsePolyZp.front! g ∈ g.toList := mem_getFirst_toList g hpos
          exact hwf (SparsePolyZp.front! g) (Array.mem_def.mpr hmem_front)
        rw [hp_div]
        apply (Nat.div_lt_iff_lt_mul (Nat.Prime.pos hp.out)).mpr
        calc
          y.1.deg < 2 ^ 64 := hdeg
          _ = 2 ^ 63 * 2 := by ring
          _ ≤ 2 ^ 63 * p := Nat.mul_le_mul_left _ hp_ge_2
      have h_deg_eq : (UMonomial.mk ((y.1.deg.toUInt64 / (SparsePolyZp.front! g).snd.prime).toInt64)).deg =
          (y.1.deg.toUInt64 / (SparsePolyZp.front! g).snd.prime).toNat := by
        calc
          (UMonomial.mk ((y.1.deg.toUInt64 / (SparsePolyZp.front! g).snd.prime).toInt64)).deg
              = ((y.1.deg.toUInt64 / (SparsePolyZp.front! g).snd.prime).toInt64).toNatClampNeg := rfl
          _ = (y.1.deg.toUInt64 / (SparsePolyZp.front! g).snd.prime).toNat :=
            UInt64_toInt64_toNatClampNeg_eq_toNat_of_lt hbound
      rw [h_deg_eq]
      calc
        (y.1.deg.toUInt64 / (SparsePolyZp.front! g).snd.prime).toNat
            = (y.1.deg.toUInt64).toNat / ((SparsePolyZp.front! g).snd.prime).toNat := by rw [UInt64_toNat_div]
        _ = y.1.deg / ((SparsePolyZp.front! g).snd.prime).toNat := by
          have h_degU64_toNat : (y.1.deg.toUInt64).toNat = y.1.deg := by
            have : y.1.deg % 2 ^ 64 = y.1.deg := Nat.mod_eq_of_lt (h_deg_bound y hy)
            simpa [UInt64.toNat_ofNat] using this
          rw [h_degU64_toNat]
        _ ≤ y.1.deg := Nat.div_le_self _ _
    have h_div_ov : y.2.val.toNat * (UMonomial.mk ((y.1.deg.toUInt64 / (SparsePolyZp.front! g).snd.prime).toInt64)).deg < 2 ^ 64 := by
      calc
        y.2.val.toNat * (UMonomial.mk ((y.1.deg.toUInt64 / (SparsePolyZp.front! g).snd.prime).toInt64)).deg
            ≤ y.2.val.toNat * y.1.deg := Nat.mul_le_mul_left (y.2.val.toNat) hdiv_le
        _ < 2 ^ 64 := h_ov
    simpa using h_div_ov
  · intro x hx
    rw [h_loop] at hx
    rcases List.mem_map.mp hx with ⟨y, hy, rfl⟩
    have hdeg : y.1.deg < 2 ^ 64 := h_deg_bound y hy
    have h_wf_y : y.2.prime.toNat = p := hwf y (Array.mem_def.mpr hy)
    have hprime_front : ((SparsePolyZp.front! g).snd.prime).toNat = p := by
      have hpos : 0 < g.size := by
        have : (g.toList).length > 0 := List.length_pos_of_mem hy
        have hsize : g.size = (g.toList).length := by simp
        omega
      have hmem_front : SparsePolyZp.front! g ∈ g.toList := mem_getFirst_toList g hpos
      exact hwf (SparsePolyZp.front! g) (Array.mem_def.mpr hmem_front)
    have h_deg_eq : (UMonomial.mk ((y.1.deg.toUInt64 / (SparsePolyZp.front! g).snd.prime).toInt64)).deg =
        (y.1.deg.toUInt64 / (SparsePolyZp.front! g).snd.prime).toNat := by
      have hp_ge_2 : 2 ≤ p := Nat.Prime.two_le hp.out
      have hbound : (y.1.deg.toUInt64 / (SparsePolyZp.front! g).snd.prime).toNat < 2 ^ 63 := by
        rw [UInt64_toNat_div]
        have h_degU64_toNat : (y.1.deg.toUInt64).toNat = y.1.deg := by
          have : y.1.deg % 2 ^ 64 = y.1.deg := Nat.mod_eq_of_lt hdeg
          simpa [UInt64.toNat_ofNat] using this
        rw [h_degU64_toNat, hprime_front]
        apply (Nat.div_lt_iff_lt_mul (Nat.Prime.pos hp.out)).mpr
        calc
          y.1.deg < 2 ^ 64 := hdeg
          _ = 2 ^ 63 * 2 := by ring
          _ ≤ 2 ^ 63 * p := Nat.mul_le_mul_left _ hp_ge_2
      calc
        (UMonomial.mk ((y.1.deg.toUInt64 / (SparsePolyZp.front! g).snd.prime).toInt64)).deg
            = ((y.1.deg.toUInt64 / (SparsePolyZp.front! g).snd.prime).toInt64).toNatClampNeg := rfl
        _ = (y.1.deg.toUInt64 / (SparsePolyZp.front! g).snd.prime).toNat :=
          UInt64_toInt64_toNatClampNeg_eq_toNat_of_lt hbound
    have h_p_deg : p * (UMonomial.mk ((y.1.deg.toUInt64 / (SparsePolyZp.front! g).snd.prime).toInt64)).deg < 2 ^ 64 := by
      calc
        p * (UMonomial.mk ((y.1.deg.toUInt64 / (SparsePolyZp.front! g).snd.prime).toInt64)).deg
            = p * ((y.1.deg.toUInt64 / (SparsePolyZp.front! g).snd.prime).toNat) := by rw [h_deg_eq]
        _ = p * ((y.1.deg.toUInt64).toNat / ((SparsePolyZp.front! g).snd.prime).toNat) := by rw [UInt64_toNat_div]
        _ = p * (y.1.deg / ((SparsePolyZp.front! g).snd.prime).toNat) := by
          have h_degU64_toNat : (y.1.deg.toUInt64).toNat = y.1.deg := by
            have : y.1.deg % 2 ^ 64 = y.1.deg := Nat.mod_eq_of_lt hdeg
            simpa [UInt64.toNat_ofNat] using this
          rw [h_degU64_toNat]
        _ = p * (y.1.deg / p) := by rw [hprime_front]
        _ ≤ y.1.deg := Nat.mul_div_le _ _
        _ < 2 ^ 64 := hdeg
    simpa

/-- __upoly_make_monic_ir 保持 no_overflow（需 p * deg < 2^64 保证缩放后的 val * deg < 2^64）。 -/
lemma upoly_make_monic_no_overflow (f : SparsePolyZp) (hred : SparsePolyZp.AllReduced p f.toList)
    (h_no_overflow : ∀ x ∈ f.toList, x.2.val.toNat * x.1.deg < 2 ^ 64)
    (h_p_mul_deg_lt : ∀ x ∈ f.toList, p * x.1.deg < 2 ^ 64) (hp_size : 2 * p ≤ UInt64.size) :
    ∀ x ∈ (Generated.__upoly_make_monic_ir f).snd.toList, x.2.val.toNat * x.1.deg < 2 ^ 64 := by
  intro x hx
  simp [Generated.__upoly_make_monic_ir, Generated.__upoly_make_monic_ir_def] at hx ⊢
  split_ifs at hx with h
  · simp at hx; exact h_no_overflow x hx.val
  · simp at hx
    have hloop : (Generated._loop___upoly_make_monic_0_ir 0 f (Zp.inv ((SparsePolyZp.front! f).snd))).snd.toList =
        (f.toList).map (fun (m, x) => (m, x * (Zp.inv ((SparsePolyZp.front! f).snd)))) := by
      have := loop_result_toList f (Zp.inv ((SparsePolyZp.front! f).snd)) 0 (Nat.zero_le _)
      simpa [List.take, List.drop] using this
    have hx_val : x ∈ (Generated._loop___upoly_make_monic_0_ir 0 f (Zp.inv ((SparsePolyZp.front! f).snd))).snd.toList :=
      (show x ∈ ((Generated._loop___upoly_make_monic_0_ir 0 f (Zp.inv ((SparsePolyZp.front! f).snd))).snd : SparsePolyZp) from hx).val
    have hx_mem : x ∈ (f.toList).map (fun (m, x) => (m, x * (Zp.inv ((SparsePolyZp.front! f).snd)))) := by
      rw [← hloop]; exact hx_val
    rcases List.mem_map.mp hx_mem with ⟨y, hy, rfl⟩
    have hy_red : Zp.Reduced p y.2 := hred y hy
    have hp_pos : 0 < p := Nat.Prime.pos (hp.out)
    have hp_le_size : p ≤ UInt64.size := by
      have : 2 * p ≤ UInt64.size := hp_size
      have : p ≤ 2 * p := by omega
      omega
    have h_new_val_lt_p : (y.2 * (Zp.inv ((SparsePolyZp.front! f).snd))).val.toNat < p := by
      have hp_le_size : p ≤ UInt64.size := by
        have : 2 * p ≤ UInt64.size := hp_size
        have : p ≤ 2 * p := by omega
        omega
      have h_mod_lt_uint64 : (y.2.val.toNat * (Zp.inv ((SparsePolyZp.front! f).snd)).val.toNat) % y.2.prime.toNat < UInt64.size := by
        have h_mod_lt_p : (y.2.val.toNat * (Zp.inv ((SparsePolyZp.front! f).snd)).val.toNat) % y.2.prime.toNat < y.2.prime.toNat :=
          Nat.mod_lt _ (by
            have : y.2.prime.toNat = p := hy_red.1
            rw [this]; exact hp_pos)
        have h_p_eq : y.2.prime.toNat = p := hy_red.1
        have h_mod_lt_p' : (y.2.val.toNat * (Zp.inv ((SparsePolyZp.front! f).snd)).val.toNat) % p < p := by
          simpa [h_p_eq] using h_mod_lt_p
        exact calc
          (y.2.val.toNat * (Zp.inv ((SparsePolyZp.front! f).snd)).val.toNat) % y.2.prime.toNat
              = (y.2.val.toNat * (Zp.inv ((SparsePolyZp.front! f).snd)).val.toNat) % p := by simp [hy_red.1]
          _ < p := h_mod_lt_p'
          _ ≤ UInt64.size := hp_le_size
      have h_val : (y.2 * (Zp.inv ((SparsePolyZp.front! f).snd))).val.toNat =
          (y.2.val.toNat * (Zp.inv ((SparsePolyZp.front! f).snd)).val.toNat) % y.2.prime.toNat := by
        calc
          (y.2 * (Zp.inv ((SparsePolyZp.front! f).snd))).val.toNat
              = ((y.2.val.toNat * (Zp.inv ((SparsePolyZp.front! f).snd)).val.toNat) % y.2.prime.toNat).toUInt64.toNat := rfl
          _ = (y.2.val.toNat * (Zp.inv ((SparsePolyZp.front! f).snd)).val.toNat) % y.2.prime.toNat := by
            simp [UInt64.toNat_ofNat, h_mod_lt_uint64]
      rw [h_val, hy_red.1]
      exact Nat.mod_lt _ hp_pos
    have hp_mul_deg_lt : p * y.1.deg < 2 ^ 64 := h_p_mul_deg_lt y hy
    nlinarith

private lemma Φ_deg_div (x : UMonomial × Zp) (p_1 : UInt64) (hp_1_eq_p : p_1.toNat = p) (hdeg : (x.1.deg : ℕ) < 2 ^ 64) :
    ((UMonomial.mk ((x.1.deg.toUInt64 / p_1).toInt64)).deg : ℕ) = (x.1.deg : ℕ) / p := by
  have hp_ge_2 : 2 ≤ p := Nat.Prime.two_le hp.out
  have hbound : (x.1.deg.toUInt64 / p_1).toNat < 2 ^ 63 := by
    rw [UInt64_toNat_div]
    have h_degU64_toNat : (x.1.deg.toUInt64).toNat = x.1.deg := by
      have : x.1.deg % 2 ^ 64 = x.1.deg := Nat.mod_eq_of_lt hdeg
      simpa [UInt64.toNat_ofNat] using this
    rw [h_degU64_toNat, hp_1_eq_p]
    apply (Nat.div_lt_iff_lt_mul (Nat.Prime.pos hp.out)).mpr
    calc
      x.1.deg < 2 ^ 64 := hdeg
      _ = 2 ^ 63 * 2 := by ring
      _ ≤ 2 ^ 63 * p := Nat.mul_le_mul_left _ hp_ge_2
  calc
    ((UMonomial.mk ((x.1.deg.toUInt64 / p_1).toInt64)).deg : ℕ)
        = ((x.1.deg.toUInt64 / p_1).toInt64).toNatClampNeg := rfl
    _ = (x.1.deg.toUInt64 / p_1).toNat := UInt64_toInt64_toNatClampNeg_eq_toNat_of_lt hbound
    _ = (x.1.deg.toUInt64).toNat / p_1.toNat := by rw [UInt64_toNat_div]
    _ = x.1.deg / p_1.toNat := by
      have h_degU64_toNat : (x.1.deg.toUInt64).toNat = x.1.deg := by
        have : x.1.deg % 2 ^ 64 = x.1.deg := Nat.mod_eq_of_lt hdeg
        simpa [UInt64.toNat_ofNat] using this
      rw [h_degU64_toNat]
    _ = (x.1.deg : ℕ) / p := by rw [hp_1_eq_p]

private lemma sortedListB_map_of_strict
    (F : UMonomial × Zp → UMonomial × Zp) :
    ∀ l : List (UMonomial × Zp), SparsePolyZp.sortedListB l = true →
      (∀ a ∈ l, ∀ b ∈ l, b.1.deg < a.1.deg → (F b).1.deg < (F a).1.deg) →
      SparsePolyZp.sortedListB (l.map F) = true := by
  intro l
  induction l with
  | nil => simp [SparsePolyZp.sortedListB]
  | cons a rest ih =>
      intro hs hstrict
      rw [List.map_cons, SparsePolyZp.sortedListB_iff]
      have hst := (SparsePolyZp.sortedListB_iff a rest).mp hs
      constructor
      · intro z hz
        rcases List.mem_map.mp hz with ⟨b, hb, rfl⟩
        exact hstrict a (by simp) b (List.mem_cons_of_mem a hb) (hst.1 b hb)
      · exact ih hst.2 (by
          intro x hx y hy hxy
          exact hstrict x (List.mem_cons_of_mem a hx) y
            (List.mem_cons_of_mem a hy) hxy)

/-- `derivative = 0` 分支中，p 次根保持稀疏规范形和严格次数顺序。 -/
private lemma extract_pth_root_canonical (g : SparsePolyZp) (hg : CanonicalRep p g)
    (hgne : ¬g.isEmpty)
    (hderiv : SparsePolyZp.derivative g = SparsePolyZp.empty)
    (h_no_overflow : ∀ x ∈ g.toList, x.2.val.toNat * x.1.deg < 2 ^ 64)
    (h_deg_bound : ∀ x ∈ g.toList, x.1.deg < 2 ^ 64) :
    CanonicalRep p (Generated.__extract_pth_root_ir g) := by
  let p_1 := (SparsePolyZp.front! g).snd.prime
  let F : UMonomial × Zp → UMonomial × Zp := fun x =>
    (UMonomial.mk ((x.1.deg.toUInt64 / p_1).toInt64), x.2)
  have hlist : (Generated.__extract_pth_root_ir g).toList = g.toList.map F := by
    unfold Generated.__extract_pth_root_ir Generated.__extract_pth_root_ir_def
    have h := loop_extract_toList g SparsePolyZp.empty 0 p_1 (by simp)
    simpa [SparsePolyZp.empty, List.drop, F, p_1] using h
  have hfrontmem : g[0]! ∈ g.toList := mem_getFirst_toList g (by
    rw [Array.isEmpty_iff_size_eq_zero] at hgne
    omega)
  have hp1 : p_1.toNat = p := by
    simpa [p_1] using (hg.2.2 g[0]! hfrontmem).1
  have hdvd := degree_dvd_of_derivative_empty (p := p) g hg hgne hderiv
    h_no_overflow h_deg_bound
  have hsorted : SparsePolyZp.Sorted (Generated.__extract_pth_root_ir g) := by
    unfold SparsePolyZp.Sorted
    rw [hlist]
    apply sortedListB_map_of_strict F g.toList hg.1
    intro a ha b hb hba
    rw [show (F b).1.deg = b.1.deg / p by
          exact Φ_deg_div (p := p) b p_1 hp1 (h_deg_bound b hb),
      show (F a).1.deg = a.1.deg / p by
          exact Φ_deg_div (p := p) a p_1 hp1 (h_deg_bound a ha)]
    apply (Nat.div_lt_iff_lt_mul hp.out.pos).mpr
    rw [Nat.div_mul_cancel (hdvd a ha)]
    exact hba
  have hnonzero : SparsePolyZp.NonZeroB (Generated.__extract_pth_root_ir g) := by
    apply nonzeroB_of_val_nonzero
    intro x hx
    have hxnat := extract_pth_root_val_nonzero g
      (fun y hy => by
        exact fun hz => (val_nonzero_of_nonzeroB g hg.2.1 y hy)
          (UInt64.toNat_inj.mp (by simpa using hz))) x hx
    exact fun hz => hxnat (hz ▸ rfl)
  exact ⟨hsorted, hnonzero, extract_pth_root_allReduced g hg.2.2⟩

-- ============================================================
-- §2b. 辅助引理：Φ 映射下的系数对应
-- ============================================================

private lemma coeff_listSum_eq (xs : List (UMonomial × Zp)) (k : ℕ) :
    coeff (listSum p xs) k = ((xs.filter (λ x => (x.1.deg : ℕ) = k)).map (λ x => Zp.toZMod p x.2)).sum := by
  induction xs with
  | nil => simp [listSum]
  | cons x xs ih =>
    rcases x with ⟨m, c⟩
    simp [listSum, coeff_add, coeff_monomial, ih]
    by_cases h : (m.deg : ℕ) = k
    · simp [h]
    · simp [h]

private lemma filter_map_Φ_eq (xs : List (UMonomial × Zp)) (p_1 : UInt64) (k : ℕ) (Φ : UMonomial × Zp → UMonomial × Zp) :
    ((xs.map Φ).filter (λ x => (x.1.deg : ℕ) = k)) = (xs.filter (λ x => ((Φ x).1.deg : ℕ) = k)).map Φ := by
  induction xs with
  | nil => simp
  | cons x xs ih =>
    simp
    by_cases h : ((Φ x).1.deg : ℕ) = k
    · simp [h]; rw [ih]
    · simp [h]; exact ih

private lemma partition_sum_eq (xs : List (UMonomial × Zp)) (p_1 : UInt64) (hp_1_eq_p : p_1.toNat = p) (k : ℕ) (Φ : UMonomial × Zp → UMonomial × Zp)
    (h_Φ_div : ∀ x ∈ xs, ((Φ x).1.deg : ℕ) = (x.1.deg : ℕ) / p) :
    ((xs.filter (λ x => ((Φ x).1.deg : ℕ) = k)).map (λ x => Zp.toZMod p x.2)).sum =
    ((xs.filter (λ x => (x.1.deg : ℕ) = k * p)).map (λ x => Zp.toZMod p x.2)).sum +
    ((xs.filter (λ x => ((Φ x).1.deg : ℕ) = k ∧ (x.1.deg : ℕ) ≠ k * p)).map (λ x => Zp.toZMod p x.2)).sum := by
  induction xs with
  | nil => simp
  | cons a xs ih =>
    have ha_div : ((Φ a).1.deg : ℕ) = (a.1.deg : ℕ) / p := h_Φ_div a (by simp)
    simp
    by_cases h1 : ((Φ a).1.deg : ℕ) = k
    · by_cases h2 : (a.1.deg : ℕ) = k * p
      · have ih' := ih (λ x hx => h_Φ_div x (List.mem_cons_of_mem _ hx))
        simp [h1, h2, List.filter_cons, not_true_eq_false, List.map_cons, List.sum_cons, add_assoc, ih']
      · have ih' := ih (λ x hx => h_Φ_div x (List.mem_cons_of_mem _ hx))
        -- h1:True, h2:False. LHS: a included; RHS1: a not; RHS2: a included.
        -- a' + L = R1 + (a' + R2) where L = R1 + R2 → a' + R1 + R2 = R1 + a' + R2
        simp [h1, h2, List.filter_cons, not_true_eq_false, List.map_cons, List.sum_cons, ih']
        ring
    · have ih' := ih (λ x hx => h_Φ_div x (List.mem_cons_of_mem _ hx))
      by_cases h2 : (a.1.deg : ℕ) = k * p
      · -- h1: False, h2: True → impossible by ha_div
        have h_contra : ((Φ a).1.deg : ℕ) = k := by
          rw [ha_div, h2]
          have hp_pos : 0 < p := Nat.Prime.pos hp.out
          calc
            (k * p) / p = (p * k) / p := by rw [mul_comm]
            _ = k := Nat.mul_div_cancel_left k hp_pos
        exfalso; exact h1 h_contra
      · simp [h1, h2, List.filter_cons, not_true_eq_false, List.map_cons, List.sum_cons, ih']

private lemma extra_union_eq (xs : List (UMonomial × Zp)) (p_1 : UInt64) (hp_1_eq_p : p_1.toNat = p) (k : ℕ) (Φ : UMonomial × Zp → UMonomial × Zp) (hp_pos : 0 < p)
    (h_Φ_div : ∀ x ∈ xs, ((Φ x).1.deg : ℕ) = (x.1.deg : ℕ) / p) :
    List.sum (List.map (λ x : UMonomial × Zp => Zp.toZMod p x.2)
      (xs.filter (λ x => ((Φ x).1.deg : ℕ) = k ∧ (x.1.deg : ℕ) ≠ k * p))) =
    (Finset.Ico 1 p).sum (λ r : ℕ => List.sum (List.map (λ x : UMonomial × Zp => Zp.toZMod p x.2)
      (xs.filter (λ x : UMonomial × Zp => (x.1.deg : ℕ) = k * p + r)))) := by
  induction xs with
  | nil => simp
  | cons a xs ih =>
    have ha_div : ((Φ a).1.deg : ℕ) = (a.1.deg : ℕ) / p := h_Φ_div a (by simp)
    have h_Φ_div_xs : ∀ x ∈ xs, ((Φ x).1.deg : ℕ) = (x.1.deg : ℕ) / p := by
      intro x hx; exact h_Φ_div x (by simp [hx])
    by_cases hcond : ((Φ a).1.deg : ℕ) = k ∧ (a.1.deg : ℕ) ≠ k * p
    · have hΦ := hcond.1
      have hne := hcond.2
      rw [ha_div] at hΦ
      -- hΦ: a.deg / p = k, so a.deg = k*p + r where r = a.deg % p, 0 < r < p
      set r := (a.1.deg : ℕ) % p with hr_def
      have hr_pos : 0 < r := by
        by_contra hzero
        have hr_zero : r = 0 := by omega
        apply hne
        have h_total : (a.1.deg : ℕ) = k * p + r := by
          calc
            (a.1.deg : ℕ) = p * ((a.1.deg : ℕ) / p) + (a.1.deg : ℕ) % p :=
              (Nat.div_add_mod (a.1.deg : ℕ) p).symm
            _ = p * k + r := by rw [hΦ, hr_def]
            _ = k * p + r := by ring
        calc
          (a.1.deg : ℕ) = k * p + r := h_total
          _ = k * p := by simp [hr_zero]
      have hr_lt_p : r < p := Nat.mod_lt _ hp_pos
      have hr_mem : r ∈ Finset.Ico 1 p := Finset.mem_Ico.mpr ⟨hr_pos, hr_lt_p⟩
      have h_a_eq : (a.1.deg : ℕ) = k * p + r := by
        calc
          (a.1.deg : ℕ) = p * ((a.1.deg : ℕ) / p) + (a.1.deg : ℕ) % p :=
            (Nat.div_add_mod (a.1.deg : ℕ) p).symm
          _ = p * k + r := by rw [hΦ, hr_def]
          _ = k * p + r := by ring
      -- Uniqueness: a.deg = k*p + r' implies r' = r
      have h_unique : ∀ r', (a.1.deg : ℕ) = k * p + r' → r' = r := by
        intro r' heq; omega
      -- Compute RHS(a::xs) = a.2.toZMod p + RHS(xs) as a separate lemma
      have h_rhs_decomp : (Finset.Ico 1 p).sum (λ r' : ℕ => List.sum (List.map (λ x : UMonomial × Zp => Zp.toZMod p x.2)
          ((a :: xs).filter (λ x : UMonomial × Zp => (x.1.deg : ℕ) = k * p + r')))) =
          Zp.toZMod p a.2 +
          (Finset.Ico 1 p).sum (λ r' : ℕ => List.sum (List.map (λ x : UMonomial × Zp => Zp.toZMod p x.2)
            (xs.filter (λ x : UMonomial × Zp => (x.1.deg : ℕ) = k * p + r')))) := by
        calc
          (Finset.Ico 1 p).sum (λ r' : ℕ => List.sum (List.map (λ x : UMonomial × Zp => Zp.toZMod p x.2)
            ((a :: xs).filter (λ x : UMonomial × Zp => (x.1.deg : ℕ) = k * p + r')))) =
              (Finset.Ico 1 p).sum (λ r' : ℕ =>
                (if (a.1.deg : ℕ) = k * p + r' then Zp.toZMod p a.2 else 0) +
                List.sum (List.map (λ x : UMonomial × Zp => Zp.toZMod p x.2)
                  (xs.filter (λ x : UMonomial × Zp => (x.1.deg : ℕ) = k * p + r')))) := by
            refine Finset.sum_congr rfl (λ r' hr' => ?_)
            by_cases ha_eq' : (a.1.deg : ℕ) = k * p + r'
            · simp [ha_eq']
            · simp [ha_eq']
          _ = ((Finset.Ico 1 p).sum (λ r' : ℕ => if (a.1.deg : ℕ) = k * p + r' then Zp.toZMod p a.2 else 0)) +
              (Finset.Ico 1 p).sum (λ r' : ℕ => List.sum (List.map (λ x : UMonomial × Zp => Zp.toZMod p x.2)
                (xs.filter (λ x : UMonomial × Zp => (x.1.deg : ℕ) = k * p + r')))) := by
            rw [Finset.sum_add_distrib]
          _ = Zp.toZMod p a.2 +
              (Finset.Ico 1 p).sum (λ r' : ℕ => List.sum (List.map (λ x : UMonomial × Zp => Zp.toZMod p x.2)
                (xs.filter (λ x : UMonomial × Zp => (x.1.deg : ℕ) = k * p + r')))) := by
            have h_sum_ite : (Finset.Ico 1 p).sum (λ r' : ℕ => if (a.1.deg : ℕ) = k * p + r' then Zp.toZMod p a.2 else 0) =
                Zp.toZMod p a.2 := by
              -- Only r' = r satisfies the condition, and r ∈ Ico 1 p
              have h_fn_eq : (λ r' : ℕ => if (a.1.deg : ℕ) = k * p + r' then Zp.toZMod p a.2 else 0) =
                             (λ r' : ℕ => if r' = r then Zp.toZMod p a.2 else 0) := by
                ext r'
                by_cases h_eq : (a.1.deg : ℕ) = k * p + r'
                · have hr_eq : r' = r := h_unique r' h_eq
                  simp [h_eq, hr_eq]
                · have hr_ne : r' ≠ r := by
                    intro hsub; apply h_eq; rw [hsub]; exact h_a_eq
                  simp [h_eq, hr_ne]
              rw [h_fn_eq]
              -- Use Finset.sum_ite_eq to compute the sum
              calc
                (Finset.Ico 1 p).sum (λ r' : ℕ => if r' = r then Zp.toZMod p a.2 else 0)
                    = (Finset.Ico 1 p).sum (λ r' : ℕ => if r = r' then Zp.toZMod p a.2 else 0) := by
                  refine Finset.sum_congr rfl (λ r' hr' => ?_)
                  simp [eq_comm]
                _ = (if r ∈ Finset.Ico 1 p then Zp.toZMod p a.2 else 0) := by
                  rw [Finset.sum_ite_eq (Finset.Ico 1 p) r (λ _ => Zp.toZMod p a.2)]
                _ = Zp.toZMod p a.2 := by simp [hr_mem]
            rw [h_sum_ite]
      -- Now prove the main equality using the decomposition
      calc
        List.sum (List.map (λ x : UMonomial × Zp => Zp.toZMod p x.2)
          ((a :: xs).filter (λ x => ((Φ x).1.deg : ℕ) = k ∧ (x.1.deg : ℕ) ≠ k * p))) =
            Zp.toZMod p a.2 +
            List.sum (List.map (λ x : UMonomial × Zp => Zp.toZMod p x.2)
              (xs.filter (λ x => ((Φ x).1.deg : ℕ) = k ∧ (x.1.deg : ℕ) ≠ k * p))) := by
          simp [hcond]
        _ = Zp.toZMod p a.2 +
            (Finset.Ico 1 p).sum (λ r' : ℕ => List.sum (List.map (λ x : UMonomial × Zp => Zp.toZMod p x.2)
              (xs.filter (λ x : UMonomial × Zp => (x.1.deg : ℕ) = k * p + r')))) := by
          rw [ih h_Φ_div_xs]
        _ = (Finset.Ico 1 p).sum (λ r' : ℕ => List.sum (List.map (λ x : UMonomial × Zp => Zp.toZMod p x.2)
              ((a :: xs).filter (λ x : UMonomial × Zp => (x.1.deg : ℕ) = k * p + r')))) := by
          rw [h_rhs_decomp]
    · -- a does not satisfy the condition: ¬(Φ(a).deg = k ∧ a.deg ≠ k*p)
      -- So either Φ(a).deg ≠ k, or a.deg = k*p
      have h_no_match : ∀ r' ∈ Finset.Ico 1 p, (a.1.deg : ℕ) ≠ k * p + r' := by
        intro r' hr'
        by_cases hΦ_eq' : ((Φ a).1.deg : ℕ) = k
        · -- a.deg/p = k (from ha_div). From hcond: a.deg = k*p
          have hΦ_eq_orig : ((Φ a).1.deg : ℕ) = k := hΦ_eq'
          rw [ha_div] at hΦ_eq'
          by_cases h_ne_kp : (a.1.deg : ℕ) ≠ k * p
          · exfalso; apply hcond; exact ⟨hΦ_eq_orig, h_ne_kp⟩
          · -- a.deg = k*p. Then a.deg ≠ k*p + r' since r' ≥ 1
            have h_kp : (a.1.deg : ℕ) = k * p := by omega
            rw [h_kp]
            have hr'_pos : 1 ≤ r' := (Finset.mem_Ico.mp hr').1
            omega
        · -- Φ(a).deg ≠ k, so a.deg/p ≠ k
          rw [ha_div] at hΦ_eq'
          intro heq
          apply hΦ_eq'
          rw [heq]
          have hr'_lt_p : r' < p := (Finset.mem_Ico.mp hr').2
          have h_lt : k * p + r' < (k + 1) * p :=
            calc
              k * p + r' < k * p + p := Nat.add_lt_add_left hr'_lt_p _
              _ = (k + 1) * p := by ring
          exact Nat.div_eq_of_lt_le (Nat.le_add_right _ _) h_lt
      calc
        List.sum (List.map (λ x : UMonomial × Zp => Zp.toZMod p x.2)
          ((a :: xs).filter (λ x => ((Φ x).1.deg : ℕ) = k ∧ (x.1.deg : ℕ) ≠ k * p))) =
            List.sum (List.map (λ x : UMonomial × Zp => Zp.toZMod p x.2)
              (xs.filter (λ x => ((Φ x).1.deg : ℕ) = k ∧ (x.1.deg : ℕ) ≠ k * p))) := by
          simp [hcond]
        _ = (Finset.Ico 1 p).sum (λ r' : ℕ => List.sum (List.map (λ x : UMonomial × Zp => Zp.toZMod p x.2)
              (xs.filter (λ x : UMonomial × Zp => (x.1.deg : ℕ) = k * p + r')))) := by
          rw [ih h_Φ_div_xs]
        _ = (Finset.Ico 1 p).sum (λ r' : ℕ => List.sum (List.map (λ x : UMonomial × Zp => Zp.toZMod p x.2)
              ((a :: xs).filter (λ x : UMonomial × Zp => (x.1.deg : ℕ) = k * p + r')))) := by
          refine Finset.sum_congr rfl (λ r' hr' => ?_)
          simp [h_no_match r' hr']
private lemma coeff_listSum_map_Φ_eq_coeff_listSum_mul (xs : List (UMonomial × Zp)) (p_1 : UInt64) (hp_1_eq_p : p_1.toNat = p) (k : ℕ)
    (hcoeff : ∀ n, ¬ p ∣ n → coeff (listSum p xs) n = 0) (hp_pos : 0 < p)
    (hdeg_bound : ∀ x ∈ xs, (x.1.deg : ℕ) < 2 ^ 64) :
    coeff (listSum p (xs.map (λ (m,c) => (UMonomial.mk ((m.deg.toUInt64 / p_1).toInt64), c)))) k =
    coeff (listSum p xs) (k * p) := by
  let Φ : UMonomial × Zp → UMonomial × Zp := λ (m,c) => (UMonomial.mk ((m.deg.toUInt64 / p_1).toInt64), c)
  have h_Φ_div_on_xs : ∀ x ∈ xs, ((Φ x).1.deg : ℕ) = (x.1.deg : ℕ) / p := by
    intro x hx
    dsimp [Φ]
    apply Φ_deg_div x p_1 hp_1_eq_p
    exact hdeg_bound x hx
  rw [coeff_listSum_eq (xs.map Φ) k, coeff_listSum_eq xs (k * p)]
  rw [filter_map_Φ_eq xs p_1 k Φ, List.map_map]
  have h_extra : ((xs.filter (λ x => ((Φ x).1.deg : ℕ) = k ∧ (x.1.deg : ℕ) ≠ k * p)).map (λ x => Zp.toZMod p x.2)).sum = 0 := by
    have h_extra_union : ((xs.filter (λ x => ((Φ x).1.deg : ℕ) = k ∧ (x.1.deg : ℕ) ≠ k * p)).map (λ x => Zp.toZMod p x.2)).sum =
        (Finset.Ico 1 p).sum (λ r : ℕ => ((xs.filter (λ x : UMonomial × Zp => (x.1.deg : ℕ) = k * p + r)).map (λ x => Zp.toZMod p x.2)).sum) := by
      apply extra_union_eq xs p_1 hp_1_eq_p k Φ hp_pos h_Φ_div_on_xs
    rw [h_extra_union]
    apply Finset.sum_eq_zero
    intro r hr
    have hr_pos : 0 < r := (Finset.mem_Ico.mp hr).1
    have hr_lt_p : r < p := (Finset.mem_Ico.mp hr).2
    have h_not_dvd : ¬ p ∣ k * p + r := by
      intro h_dvd
      have h_p_dvd_kp : p ∣ k * p := Nat.dvd_mul_left p k
      have h_p_dvd_r : p ∣ r := by
        -- from h_dvd: p ∣ k * p + r and h_p_dvd_kp: p ∣ k * p
        have := (Nat.dvd_add_right h_p_dvd_kp).mp h_dvd
        exact this
      have h_r_zero : r = 0 := Nat.eq_zero_of_dvd_of_lt h_p_dvd_r hr_lt_p
      omega
    have h_zero : coeff (listSum p xs) (k * p + r) = 0 := hcoeff (k * p + r) h_not_dvd
    rw [coeff_listSum_eq xs (k * p + r)] at h_zero
    simpa using h_zero
  calc
    ((xs.filter (λ x => ((Φ x).1.deg : ℕ) = k)).map (λ x => Zp.toZMod p x.2)).sum
        = ((xs.filter (λ x => (x.1.deg : ℕ) = k * p)).map (λ x => Zp.toZMod p x.2)).sum +
          ((xs.filter (λ x => ((Φ x).1.deg : ℕ) = k ∧ (x.1.deg : ℕ) ≠ k * p)).map (λ x => Zp.toZMod p x.2)).sum := by
      apply partition_sum_eq xs p_1 hp_1_eq_p k Φ h_Φ_div_on_xs
    _ = ((xs.filter (λ x => (x.1.deg : ℕ) = k * p)).map (λ x => Zp.toZMod p x.2)).sum + 0 := by rw [h_extra]
    _ = ((xs.filter (λ x => (x.1.deg : ℕ) = k * p)).map (λ x => Zp.toZMod p x.2)).sum := by simp

/-- __extract_pth_root_ir 的 toPoly 对应 contract p。 -/
lemma extract_pth_root_toPoly_eq (g : SparsePolyZp) (h_deriv0 : SparsePolyZp.derivative g = SparsePolyZp.empty)
    (hwf : SparsePolyZp.WellFormed p g) (hred : SparsePolyZp.AllReduced p g.toList)
    (hp_size : 2 * p ≤ UInt64.size) (h_no_overflow : ∀ x ∈ g.toList, x.2.val.toNat * x.1.deg < 2 ^ 64)
    (h_deg_bound : ∀ x ∈ g.toList, x.1.deg < 2 ^ 64) :
    SparsePolyZp.toPoly p (Generated.__extract_pth_root_ir g) = Polynomial.contract p (SparsePolyZp.toPoly p g) := by
  by_cases hsize_zero : Array.size g = 0
  · -- Empty case: g = #[], both sides evaluate to 0. In practice this case is never reached
    -- (the main theorem ensures natDegree > 0 → g non-empty), but the lemma is stated
    -- unconditionally, so we discharge it directly.
    have hg_empty : g = #[] := Array.eq_empty_of_size_eq_zero hsize_zero
    subst hg_empty
    have hp_ne_zero : p ≠ 0 := Nat.Prime.ne_zero hp.out
    -- LHS: __extract_pth_root_ir #[] 的 toList 为空（循环遍历空数组），故 toPoly = 0
    have hL : SparsePolyZp.toPoly p (Generated.__extract_pth_root_ir (#[] : SparsePolyZp)) = 0 := by
      have h_tl : (Generated.__extract_pth_root_ir (#[] : SparsePolyZp)).toList = [] := by
        unfold Generated.__extract_pth_root_ir Generated.__extract_pth_root_ir_def
        have h_loop' := loop_extract_toList (#[] : SparsePolyZp) SparsePolyZp.empty 0
          (SparsePolyZp.front! (#[] : SparsePolyZp)).snd.prime (by simp)
        simpa [SparsePolyZp.empty, List.drop] using h_loop'
      show listSum p (Generated.__extract_pth_root_ir (#[] : SparsePolyZp)).toList = 0
      rw [h_tl]; simp [listSum]
    -- RHS: contract p 0 = 0（逐系数：coeff (contract p 0) n = coeff 0 (n*p) = 0）
    have hR : Polynomial.contract p (SparsePolyZp.toPoly p (#[] : SparsePolyZp)) = 0 := by
      rw [SparsePolyZp.toPoly_empty]
      apply Polynomial.ext; intro n
      rw [Polynomial.coeff_contract hp_ne_zero]; simp
    rw [hL, hR]
  · have hsize_pos : 0 < Array.size g := Nat.pos_of_ne_zero hsize_zero
    let p_1 : UInt64 := (SparsePolyZp.front! g).snd.prime
    have hp_1_eq_p : p_1.toNat = p := by
      have hfront_wf : (SparsePolyZp.front! g).snd.prime.toNat = p :=
        hwf (SparsePolyZp.front! g) (by
          apply Array.Mem.mk
          simpa [SparsePolyZp.front!] using mem_getFirst_toList g hsize_pos)
      simp [p_1, hfront_wf]
    have hp_ne_zero : p ≠ 0 := Nat.Prime.ne_zero (hp.out)
    have hp_pos : 0 < p := Nat.Prime.pos (hp.out)
    have h_deriv_poly : Polynomial.derivative (SparsePolyZp.toPoly p g) = 0 := by
      rw [← derivative_toPoly_eq g hwf hp_size h_no_overflow h_deg_bound, h_deriv0]
      exact SparsePolyZp.toPoly_empty p
    have h_coeff_zero_not_dvd : ∀ n, ¬ p ∣ n → coeff (SparsePolyZp.toPoly p g) n = 0 := by
      intro n hpn
      have h_expand : SparsePolyZp.toPoly p g =
          Polynomial.expand (ZMod p) p (Polynomial.contract p (SparsePolyZp.toPoly p g)) :=
        (Polynomial.expand_contract p h_deriv_poly hp_ne_zero).symm
      rw [h_expand, Polynomial.coeff_expand hp_pos, if_neg hpn]
    have h_coeff_eq : ∀ n, coeff (SparsePolyZp.toPoly p (Generated.__extract_pth_root_ir g)) n =
        coeff (Polynomial.contract p (SparsePolyZp.toPoly p g)) n := by
      intro n
      have h_coeff_contract : coeff (Polynomial.contract p (SparsePolyZp.toPoly p g)) n = coeff (SparsePolyZp.toPoly p g) (n * p) :=
        Polynomial.coeff_contract hp_ne_zero (SparsePolyZp.toPoly p g) n
      rw [h_coeff_contract]
      have h_toPoly_listSum (f : SparsePolyZp) : SparsePolyZp.toPoly p f = listSum p f.toList := by rfl
      rw [h_toPoly_listSum (Generated.__extract_pth_root_ir g), h_toPoly_listSum g]
      have h_loop_toList : (Generated.__extract_pth_root_ir g).toList =
          g.toList.map (fun (term : UMonomial × Zp) =>
            (UMonomial.mk ((term.1.deg.toUInt64 / p_1).toInt64), term.2)) := by
        unfold Generated.__extract_pth_root_ir Generated.__extract_pth_root_ir_def
        have h_loop' := loop_extract_toList g SparsePolyZp.empty 0 p_1 (by simp)
        simpa [SparsePolyZp.empty, List.drop] using h_loop'
      rw [h_loop_toList]
      have h_coeff_listSum_zero_not_dvd : ∀ n, ¬ p ∣ n → coeff (listSum p g.toList) n = 0 := by
        intro n' hpn'
        have h_temp : SparsePolyZp.toPoly p g = listSum p g.toList := by rfl
        rw [← h_temp]
        exact h_coeff_zero_not_dvd n' hpn'
      exact coeff_listSum_map_Φ_eq_coeff_listSum_mul (g.toList) p_1 hp_1_eq_p n
        h_coeff_listSum_zero_not_dvd hp_pos h_deg_bound
    haveI : Fact (Nat.Prime p) := ⟨hp.out⟩
    apply Polynomial.ext
    exact h_coeff_eq

-- ============================================================
-- §3. 主定理：__squarefree_Zp_ir ≃ sqfZp
-- ============================================================

/-- sqfZp 关于标量乘法的不变性：对于非零常数 c，sqfZp(c·f) = sqfZp(f)。
    证明使用强归纳于 f.natDegree。derivative=0 分支通过 contract 线性性处理；
    derivative≠0 分支需 normalize 的单位不变性（admit）。 -/
private lemma sqfZp_smul (c : ZMod p) (hc : c ≠ 0) (f : Polynomial (ZMod p)) : sqfZp (C c * f) = sqfZp f := by
  refine Nat.strongRecOn (motive := λ n => ∀ (g : Polynomial (ZMod p)), g.natDegree = n → sqfZp (C c * g) = sqfZp g)
    f.natDegree ?_ f rfl
  intro k IH g hk
  by_cases hzero : g.natDegree = 0
  · -- Both constant: sqfZp returns []
    have hzero_cg' : (C c * g).natDegree = 0 :=
      calc
        (C c * g).natDegree = g.natDegree := Polynomial.natDegree_C_mul hc
        _ = 0 := hzero
    unfold sqfZp
    rw [dif_pos hzero_cg', dif_pos hzero]
  · by_cases hderiv : Polynomial.derivative g = 0
    · -- derivative = 0
      have h_deriv_cg : Polynomial.derivative (C c * g) = 0 := by simp [hderiv]
      have hzero_cg : (C c * g).natDegree ≠ 0 := by
        rw [Polynomial.natDegree_C_mul hc]
        exact hzero
      have h_contract : Polynomial.contract p (C c * g) = C c * Polynomial.contract p g := by
        have hp_ne_zero : p ≠ 0 := Nat.Prime.ne_zero hp.out
        ext n
        simp [Polynomial.coeff_contract hp_ne_zero, Polynomial.coeff_C_mul]
      have h_deg_lt : (Polynomial.contract p g).natDegree < k := by
        have h_expand : Polynomial.expand (ZMod p) p (Polynomial.contract p g) = g :=
          Polynomial.expand_contract p hderiv (Nat.Prime.ne_zero hp.out)
        have h_nd_expand : (Polynomial.expand (ZMod p) p (Polynomial.contract p g)).natDegree =
            p * (Polynomial.contract p g).natDegree := by
          rw [Polynomial.natDegree_expand, mul_comm]
        rw [h_expand] at h_nd_expand
        rw [hk] at h_nd_expand
        by_cases hzero_contract : (Polynomial.contract p g).natDegree = 0
        · rw [hzero_contract, mul_zero] at h_nd_expand
          rw [← hk] at h_nd_expand
          exact (hzero h_nd_expand).elim
        · have hp_ge_2 : 2 ≤ p := Nat.Prime.two_le hp.out
          have hpos : 1 ≤ (Polynomial.contract p g).natDegree := by omega
          nlinarith
      unfold sqfZp
      rw [dif_neg hzero_cg, dif_pos h_deriv_cg, dif_neg hzero, dif_pos hderiv]
      dsimp
      rw [h_contract]
      rw [IH (Polynomial.contract p g).natDegree h_deg_lt (Polynomial.contract p g) rfl]
    · -- derivative ≠ 0: Yun algorithm。C c 是单位（c ≠ 0 于域 ZMod p），
      -- 故内部的 gcd / divByMonic / normalize 均在单位乘法下不变 ⇒ c、w 相等 ⇒ 整个结果相等。
      have hCc_unit : IsUnit (C c : Polynomial (ZMod p)) := isUnit_C.mpr hc.isUnit
      have h_deriv_cg : Polynomial.derivative (C c * g) = C c * Polynomial.derivative g :=
        Polynomial.derivative_C_mul c g
      have hderiv_cg_ne : Polynomial.derivative (C c * g) ≠ 0 := by
        rw [h_deriv_cg]; intro h
        rcases mul_eq_zero.mp h with h1 | h2
        · exact hc (Polynomial.C_eq_zero.mp h1)
        · exact hderiv h2
      have hzero_cg : (C c * g).natDegree ≠ 0 := by
        rw [Polynomial.natDegree_C_mul hc]; exact hzero
      -- (1) c 相等：gcd 在单位乘法下相伴，normalize 后相等
      have hg_assoc : Associated (C c * g) g := associated_unit_mul_left g (C c) hCc_unit
      have hg'_assoc : Associated (C c * Polynomial.derivative g) (Polynomial.derivative g) :=
        associated_unit_mul_left (Polynomial.derivative g) (C c) hCc_unit
      have hgcd_assoc : Associated (EuclideanDomain.gcd (C c * g) (Polynomial.derivative (C c * g)))
          (EuclideanDomain.gcd g (Polynomial.derivative g)) := by
        rw [h_deriv_cg]
        apply associated_of_dvd_dvd
        · exact EuclideanDomain.dvd_gcd
            (dvd_trans (EuclideanDomain.gcd_dvd_left _ _) hg_assoc.dvd)
            (dvd_trans (EuclideanDomain.gcd_dvd_right _ _) hg'_assoc.dvd)
        · exact EuclideanDomain.dvd_gcd
            (dvd_trans (EuclideanDomain.gcd_dvd_left _ _) hg_assoc.symm.dvd)
            (dvd_trans (EuclideanDomain.gcd_dvd_right _ _) hg'_assoc.symm.dvd)
      have hcc_eq : normalize (EuclideanDomain.gcd (C c * g) (Polynomial.derivative (C c * g)))
          = normalize (EuclideanDomain.gcd g (Polynomial.derivative g)) :=
        normalize_eq_normalize_iff_associated.mpr hgcd_assoc
      -- normalized gcd 非零且 monic
      have hgcd_g_ne : EuclideanDomain.gcd g (Polynomial.derivative g) ≠ 0 := by
        intro h
        have hdvd : EuclideanDomain.gcd g (Polynomial.derivative g) ∣ g :=
          EuclideanDomain.gcd_dvd_left g (Polynomial.derivative g)
        rw [h] at hdvd
        exact hzero (by rw [zero_dvd_iff.mp hdvd]; simp)
      have hcc_monic : Monic (normalize (EuclideanDomain.gcd g (Polynomial.derivative g))) :=
        monic_normalize hgcd_g_ne
      -- (2) w 相等：(C c * g) /ₘ cc = C c * (g /ₘ cc)，再由单位相伴 normalize 相等
      have h_div_eq : (C c * g) /ₘ normalize (EuclideanDomain.gcd g (Polynomial.derivative g))
          = C c * (g /ₘ normalize (EuclideanDomain.gcd g (Polynomial.derivative g))) := by
        refine (div_modByMonic_unique
          (C c * (g /ₘ normalize (EuclideanDomain.gcd g (Polynomial.derivative g))))
          (C c * (g %ₘ normalize (EuclideanDomain.gcd g (Polynomial.derivative g))))
          hcc_monic ⟨?_, ?_⟩).1
        · calc C c * (g %ₘ normalize (EuclideanDomain.gcd g (Polynomial.derivative g)))
                + normalize (EuclideanDomain.gcd g (Polynomial.derivative g))
                  * (C c * (g /ₘ normalize (EuclideanDomain.gcd g (Polynomial.derivative g))))
              = C c * (g %ₘ normalize (EuclideanDomain.gcd g (Polynomial.derivative g))
                + normalize (EuclideanDomain.gcd g (Polynomial.derivative g))
                  * (g /ₘ normalize (EuclideanDomain.gcd g (Polynomial.derivative g)))) := by ring
            _ = C c * g := by
                rw [modByMonic_add_div g (normalize (EuclideanDomain.gcd g (Polynomial.derivative g)))]
        · have hdeg : (C c * (g %ₘ normalize (EuclideanDomain.gcd g (Polynomial.derivative g)))).degree
              = (g %ₘ normalize (EuclideanDomain.gcd g (Polynomial.derivative g))).degree := by
            rw [Polynomial.degree_mul, Polynomial.degree_C hc, zero_add]
          rw [hdeg]; exact degree_modByMonic_lt g hcc_monic
      have hw_eq : normalize ((C c * g) /ₘ normalize (EuclideanDomain.gcd g (Polynomial.derivative g)))
          = normalize (g /ₘ normalize (EuclideanDomain.gcd g (Polynomial.derivative g))) := by
        rw [h_div_eq]
        exact normalize_eq_normalize_iff_associated.mpr
          (associated_unit_mul_left _ (C c) hCc_unit)
      -- 展开两侧的 else-else 分支，用 hcc_eq、hw_eq 对齐
      unfold sqfZp
      rw [dif_neg hzero_cg, dif_neg hderiv_cg_ne, dif_neg hzero, dif_neg hderiv]
      dsimp only
      -- CC1(=gcd of C c*g) → CC2(=gcd of g) 出现在 yunLoop 的依赖证明位，需用 simp 的
      -- proof-irrelevant congruence 重写（rw 会破坏 motive）
      simp only [hcc_eq, hw_eq]

private lemma squarefree_loop2_result (arr result : Array (SparsePolyZp × UInt64))
    (p_1 : UInt64) :
    (Generated._loop___squarefree_Zp_2_ir_def 0 arr result p_1).2.2.toList =
      result.toList ++ arr.toList.map (fun (g, e) => (g, e * p_1)) := by
  let F : SparsePolyZp × UInt64 → SparsePolyZp × UInt64 :=
    fun (g, e) => (g, e * p_1)
  have hw : WellFounded (fun (a b : Array (SparsePolyZp × UInt64) × Nat) =>
      a.1.size - a.2 < b.1.size - b.2) :=
    (measure (fun q : Array (SparsePolyZp × UInt64) × Nat => q.1.size - q.2)).wf
  have hinv := hw.induction (arr, 0) (C := fun q =>
      ∀ r : Array (SparsePolyZp × UInt64),
        (Generated._loop___squarefree_Zp_2_ir_def q.2 q.1 r p_1).2.2.toList =
          r.toList ++ (List.drop q.2 q.1.toList).map F) (by
    intro q ih r
    rw [Generated._loop___squarefree_Zp_2_ir_def.eq_1]
    by_cases hidx : q.2 < q.1.size
    · simp only [hidx, ↓reduceIte]
      have hset : q.1.set! q.2 q.1[q.2]! = q.1 := by
        rw [Array.set!_eq_setIfInBounds]
        apply Array.ext
        · simp
        · intro j hj1 hj2
          rw [Array.getElem_setIfInBounds hj2]
          by_cases heq : q.2 = j
          · subst j
            simp [getElem!_def, getElem?_def, hidx]
          · simp [heq]
      rw [hset]
      simp only [id_eq]
      have hmeasure : q.1.size - (q.2 + 1) < q.1.size - q.2 := by omega
      rw [ih (q.1, q.2 + 1) hmeasure
        (r.push ((q.1[q.2]!).1, (q.1[q.2]!).2 * p_1))]
      have hlen : q.2 < q.1.toList.length := by simpa using hidx
      have hdrop : List.drop q.2 q.1.toList =
          q.1[q.2]! :: List.drop (q.2 + 1) q.1.toList := by
        have hd := drop_eq_get_cons_general q.1.toList q.2 hlen
        have hget : q.1.toList.get ⟨q.2, hlen⟩ = q.1[q.2]! := by
          calc
            q.1.toList.get ⟨q.2, hlen⟩ = q.1[q.2] := by simp
            _ = q.1[q.2]! := by
              rw [getElem!_def, getElem?_def]
              simp [hidx]
        rw [hget] at hd
        exact hd
      rw [hdrop]
      simp [F, List.append_assoc]
    · simp only [hidx, ↓reduceIte]
      have hlen : q.1.toList.length ≤ q.2 := by
        simpa using (by omega : q.1.size ≤ q.2)
      simp [List.drop_eq_nil_of_le hlen])
  simpa [F] using hinv result

/-- safe wrapper：常数多项式直接返回 #[]，避免 C++ 模型不终止的问题 -/
noncomputable def __squarefree_Zp_ir_safe (p : ℕ) [hp : Fact (Nat.Prime p)] (f : SparsePolyZp) : Array (SparsePolyZp × UInt64) :=
  if (SparsePolyZp.toPoly p f).natDegree = 0 then #[]
  else Generated.__squarefree_Zp_ir f

/-- listSum 的 natDegree 被各项度数的公共上界所界（Branch A 指数界用）。 -/
private lemma listSum_natDegree_lt (p : ℕ) (xs : List (UMonomial × Zp)) (B : ℕ) (hB : 0 < B)
    (h : ∀ x ∈ xs, x.1.deg < B) : (listSum p xs).natDegree < B := by
  induction xs with
  | nil => simpa using hB
  | cons a rest ih =>
    rcases a with ⟨m, c⟩
    rw [listSum_cons]
    have h1 : (Polynomial.monomial m.deg (Zp.toZMod p c)).natDegree < B :=
      lt_of_le_of_lt (Polynomial.natDegree_monomial_le _) (h (m, c) (by simp))
    have h2 : (listSum p rest).natDegree < B := ih (fun x hx => h x (List.mem_cons_of_mem _ hx))
    exact lt_of_le_of_lt (Polynomial.natDegree_add_le _ _) (by simp only [max_lt_iff]; exact ⟨h1, h2⟩)

set_option maxHeartbeats 0 in
/-- 主正确性定理（使用 safe wrapper） -/
theorem __squarefree_Zp_ir_refines (p : ℕ) [hp : Fact (Nat.Prime p)]
    (f : SparsePolyZp)
    (hsorted_f : SparsePolyZp.Sorted f)
    (hwf_f : SparsePolyZp.WellFormed p f)
    (hred_f : SparsePolyZp.AllReduced p f.toList)
    (hp_size : 2 * p ≤ UInt64.size)
    (hp2 : p * p ≤ UInt64.size)
    (h_no_overflow : ∀ x ∈ f.toList, x.2.val.toNat * x.1.deg < 2 ^ 64)
    (h_p_deg_bound : ∀ x ∈ f.toList, p * x.1.deg < 2 ^ 64)
    (h_deg_bound : ∀ x ∈ f.toList, x.1.deg < 2 ^ 64)
    (h_signed_deg_bound : ∀ x ∈ f.toList, x.1.deg < 2 ^ 63)
    (h_val_nonzero_f : ∀ x ∈ f.toList, x.snd.val.toNat ≠ 0)
    :   toPolyList (__squarefree_Zp_ir_safe p f) p = sqfZp (SparsePolyZp.toPoly p f) := by
  unfold __squarefree_Zp_ir_safe
  by_cases h_deg0 : (SparsePolyZp.toPoly p f).natDegree = 0
  · simp [sqfZp, h_deg0, toPolyList_empty]
  · have h_deg_pos : (SparsePolyZp.toPoly p f).natDegree > 0 := by omega
    -- 强归纳于 natDegree，使用 __squarefree_Zp_ir_safe 避免 partial def 在常数时不停机
    have h_main : toPolyList (__squarefree_Zp_ir_safe p f) p = sqfZp (SparsePolyZp.toPoly p f) := by
      suffices ∀ n, ∀ (g : SparsePolyZp),
          (SparsePolyZp.toPoly p g).natDegree = n →
          SparsePolyZp.Sorted g →
          SparsePolyZp.WellFormed p g →
          SparsePolyZp.AllReduced p g.toList →
          2 * p ≤ UInt64.size →
          p * p ≤ UInt64.size →
          (∀ x ∈ g.toList, x.2.val.toNat * x.1.deg < 2 ^ 64) →
          (∀ x ∈ g.toList, p * x.1.deg < 2 ^ 64) →
          (∀ x ∈ g.toList, x.1.deg < 2 ^ 64) →
          (∀ x ∈ g.toList, x.1.deg < 2 ^ 63) →
          (∀ x ∈ g.toList, x.snd.val.toNat ≠ 0) →
          toPolyList (__squarefree_Zp_ir_safe p g) p = sqfZp (SparsePolyZp.toPoly p g) from
        this (SparsePolyZp.toPoly p f).natDegree f rfl hsorted_f hwf_f hred_f hp_size hp2
          h_no_overflow h_p_deg_bound h_deg_bound h_signed_deg_bound h_val_nonzero_f
      intro n
      induction n using Nat.strongRecOn with
      | ind n ih =>
        intro g h_deg_eq hsorted_g hwf_g hred_g hp_size hp2 h_no_overflow_g h_p_deg_bound_g
          h_deg_bound_g h_signed_deg_bound_g h_val_nonzero_g
        unfold __squarefree_Zp_ir_safe
        by_cases h_deg0_g : (SparsePolyZp.toPoly p g).natDegree = 0
        · simp [sqfZp, h_deg0_g, toPolyList_empty]
        · by_cases h_deriv0 : SparsePolyZp.derivative g = (SparsePolyZp.empty : SparsePolyZp)
          · -- Branch A: derivative = 0 → p-th root
            have hp_prime : Nat.Prime p := hp.out
            have hp_ne_zero : p ≠ 0 := Nat.Prime.ne_zero hp_prime
            have h_deriv_poly : Polynomial.derivative (SparsePolyZp.toPoly p g) = 0 := by
              calc
                Polynomial.derivative (SparsePolyZp.toPoly p g) = SparsePolyZp.toPoly p (SparsePolyZp.derivative g) := by
                  symm; exact derivative_toPoly_eq g hwf_g hp_size h_no_overflow_g h_deg_bound_g
                _ = SparsePolyZp.toPoly p (SparsePolyZp.empty) := by rw [h_deriv0]
                _ = 0 := SparsePolyZp.toPoly_empty p
            have h_eq_nd : (Polynomial.contract p (SparsePolyZp.toPoly p g)).natDegree * p =
                (SparsePolyZp.toPoly p g).natDegree := by
              have hexp : Polynomial.expand (ZMod p) p (Polynomial.contract p (SparsePolyZp.toPoly p g)) =
                  SparsePolyZp.toPoly p g :=
                Polynomial.expand_contract p h_deriv_poly hp_ne_zero
              calc
                (Polynomial.contract p (SparsePolyZp.toPoly p g)).natDegree * p
                    = (Polynomial.expand (ZMod p) p (Polynomial.contract p (SparsePolyZp.toPoly p g))).natDegree := by
                  rw [Polynomial.natDegree_expand, mul_comm]
                _ = (SparsePolyZp.toPoly p g).natDegree := by rw [hexp]
            have h_contract_pos : 0 < (Polynomial.contract p (SparsePolyZp.toPoly p g)).natDegree := by
              by_contra hzero
              have hzero' : (Polynomial.contract p (SparsePolyZp.toPoly p g)).natDegree = 0 := by omega
              rw [hzero', zero_mul] at h_eq_nd
              omega
            have h_contract_deg_lt : (Polynomial.contract p (SparsePolyZp.toPoly p g)).natDegree <
                (SparsePolyZp.toPoly p g).natDegree := by
              have hp_ge_2 : 2 ≤ p := Nat.Prime.two_le hp_prime
              by_cases hpos : (Polynomial.contract p (SparsePolyZp.toPoly p g)).natDegree = 0
              · omega
              · have hpos_mul : (Polynomial.contract p (SparsePolyZp.toPoly p g)).natDegree * 1 <
                    (Polynomial.contract p (SparsePolyZp.toPoly p g)).natDegree * p :=
                  mul_lt_mul_of_pos_left (by omega) (Nat.pos_of_ne_zero hpos)
                simp at hpos_mul
                have : (SparsePolyZp.toPoly p g).natDegree
                    = (Polynomial.contract p (SparsePolyZp.toPoly p g)).natDegree * p := by
                  symm; exact h_eq_nd
                rw [this]
                exact hpos_mul
            have h_extract_eq : SparsePolyZp.toPoly p (Generated.__extract_pth_root_ir g) =
                Polynomial.contract p (SparsePolyZp.toPoly p g) :=
              extract_pth_root_toPoly_eq g h_deriv0 hwf_g hred_g hp_size h_no_overflow_g h_deg_bound_g
            let g_1 := Generated.__extract_pth_root_ir g
            let g_2 := (Generated.__upoly_make_monic_ir g_1).snd
            have hg_nonzeroB : SparsePolyZp.NonZeroB g := by
              apply nonzeroB_of_val_nonzero
              intro x hx hv
              exact h_val_nonzero_g x hx (by simp [hv])
            have hgcan : CanonicalRep p g := ⟨hsorted_g, hg_nonzeroB, hred_g⟩
            have hgpolyne : SparsePolyZp.toPoly p g ≠ 0 := by
              intro hz
              simp [hz] at h_deg0_g
            have hgne : ¬g.isEmpty := nonempty_of_toPoly_ne_zero (p := p) g hgpolyne
            have hg1can : CanonicalRep p g_1 := by
              simpa [g_1] using extract_pth_root_canonical (p := p) g hgcan hgne h_deriv0
                h_no_overflow_g h_deg_bound_g
            have h_toPoly_g1 : SparsePolyZp.toPoly p g_1 = Polynomial.contract p (SparsePolyZp.toPoly p g) := h_extract_eq
            have h_wf_g1 : SparsePolyZp.WellFormed p g_1 := extract_pth_root_wellFormed g hwf_g
            have h_red_g1 : SparsePolyZp.AllReduced p g_1.toList := extract_pth_root_allReduced g hred_g
            have h_db_g1 : ∀ x ∈ g_1.toList, x.1.deg < 2 ^ 64 := extract_pth_root_deg_bound g h_deg_bound_g
            have h_sdb_g1 : ∀ x ∈ g_1.toList, x.1.deg < 2 ^ 63 :=
              extract_pth_root_signed_deg_bound g h_signed_deg_bound_g
            rcases @extract_pth_root_no_overflow p hp g hwf_g h_no_overflow_g h_deg_bound_g with ⟨h_nov_g1, h_p_mul_g1⟩
            have h_deg_g2_lt_n : (SparsePolyZp.toPoly p g_2).natDegree < n := by
              have h_nd_le : (SparsePolyZp.toPoly p g_2).natDegree ≤ (SparsePolyZp.toPoly p g_1).natDegree := by
                have h_toPoly_g2_eq : SparsePolyZp.toPoly p g_2 = SparsePolyZp.toPoly p (SparsePolyZp.makeMonic g_1) :=
                  __upoly_make_monic_ir_refines p g_1 h_wf_g1 h_red_g1 hp_size
                rw [h_toPoly_g2_eq]
                unfold SparsePolyZp.makeMonic
                split_ifs with h_empty
                · simp
                · have hfront_prime : ((SparsePolyZp.front! g_1).snd.inv).prime.toNat = p := by
                    have hfront_wf_g1 : ((SparsePolyZp.front! g_1).snd.prime).toNat = p :=
                      h_wf_g1 (SparsePolyZp.front! g_1) (by
                        apply Array.Mem.mk
                        have hsize_pos_g1 : 0 < Array.size g_1 := by
                          by_contra! hzero_sz
                          apply h_empty
                          apply Array.isEmpty_iff_size_eq_zero.mpr
                          omega
                        simpa [SparsePolyZp.front!] using mem_getFirst_toList g_1 hsize_pos_g1)
                    simp [Zp.inv, hfront_wf_g1]
                  have h_scalar_eq : SparsePolyZp.toPoly p (SparsePolyZp.scalarMul ((SparsePolyZp.front! g_1).snd.inv) g_1) =
                      Polynomial.C (((SparsePolyZp.front! g_1).snd.inv).toZMod p) * SparsePolyZp.toPoly p g_1 :=
                    toPoly_scalarMul ((SparsePolyZp.front! g_1).snd.inv) g_1 h_red_g1 hfront_prime hp_size
                  have h_scalar_eq' : SparsePolyZp.toPoly p (SparsePolyZp.scalarMul (g_1[0]!.snd.inv) g_1) =
                      Polynomial.C (((g_1[0]!).snd.inv).toZMod p) * SparsePolyZp.toPoly p g_1 := by
                    simpa using h_scalar_eq
                  rw [h_scalar_eq']
                  calc
                    (Polynomial.C (((SparsePolyZp.front! g_1).snd.inv).toZMod p) * SparsePolyZp.toPoly p g_1).natDegree
                        ≤ (Polynomial.C (((SparsePolyZp.front! g_1).snd.inv).toZMod p)).natDegree + (SparsePolyZp.toPoly p g_1).natDegree :=
                      Polynomial.natDegree_mul_le
                    _ = (SparsePolyZp.toPoly p g_1).natDegree := by simp
              have h_nd_g1_lt_n : (SparsePolyZp.toPoly p g_1).natDegree < n := by
                rw [h_toPoly_g1]
                rw [h_deg_eq] at h_contract_deg_lt
                exact h_contract_deg_lt
              omega
            have h_wf_g2 : SparsePolyZp.WellFormed p g_2 := upoly_make_monic_wellFormed g_1 h_wf_g1
            have h_sorted_g2 : SparsePolyZp.Sorted g_2 := upoly_make_monic_sorted g_1 hg1can.1
            have h_red_g2 : SparsePolyZp.AllReduced p g_2.toList := upoly_make_monic_allReduced g_1 h_red_g1 hp_size
            have h_db_g2 : ∀ x ∈ g_2.toList, x.1.deg < 2 ^ 64 := upoly_make_monic_deg_bound g_1 h_db_g1
            have h_sdb_g2 : ∀ x ∈ g_2.toList, x.1.deg < 2 ^ 63 :=
              upoly_make_monic_deg_bound_of (2 ^ 63) g_1 h_sdb_g1
            have h_nov_g2 : ∀ x ∈ g_2.toList, x.2.val.toNat * x.1.deg < 2 ^ 64 :=
              upoly_make_monic_no_overflow g_1 h_red_g1 h_nov_g1 h_p_mul_g1 hp_size
            have h_pdb_g2 : ∀ x ∈ g_2.toList, p * x.1.deg < 2 ^ 64 :=
              upoly_make_monic_p_deg_bound g_1 h_p_mul_g1
            have hmem_front_g1 : (SparsePolyZp.front! g_1) ∈ g_1.toList := by
              have hsize_pos : 0 < Array.size g_1 := by
                by_contra! hzero
                have h_empty : g_1 = #[] := Array.eq_empty_of_size_eq_zero (by omega)
                rw [h_empty] at h_toPoly_g1
                have : SparsePolyZp.toPoly p #[] = (0 : Polynomial (ZMod p)) := by simp
                rw [this] at h_toPoly_g1
                have hzero_nd : Polynomial.natDegree (0 : Polynomial (ZMod p)) = 0 := by simp
                rw [← h_toPoly_g1, hzero_nd] at h_contract_pos
                omega
              exact mem_getFirst_toList g_1 hsize_pos
            have h_val_nonzero_g1 : ∀ x ∈ g_1.toList, x.snd.val.toNat ≠ 0 :=
              extract_pth_root_val_nonzero g h_val_nonzero_g
            have h_lc_val_nonzero : ((SparsePolyZp.front! g_1).snd).val.toNat ≠ 0 :=
              h_val_nonzero_g1 (SparsePolyZp.front! g_1) hmem_front_g1
            have h_red_lc : Zp.Reduced p ((SparsePolyZp.front! g_1).snd) :=
              h_red_g1 (SparsePolyZp.front! g_1) hmem_front_g1
            have h_red_inv : Zp.Reduced p ((SparsePolyZp.front! g_1).snd.inv) := by
              -- inv preserves prime (= p); val < p follows from modInv returning value in [0, p)
              have ha_prime : ((SparsePolyZp.front! g_1).snd).prime.toNat = p := h_red_lc.1
              have h_val_lt : ((SparsePolyZp.front! g_1).snd.inv).val.toNat < p := by
                -- modInv returns r where r = s % p, so 0 ≤ r < p
                calc ((SparsePolyZp.front! g_1).snd.inv).val.toNat
                    = (Zp.modInv ((SparsePolyZp.front! g_1).snd).val
                        ((SparsePolyZp.front! g_1).snd).prime).toNat := rfl
                  _ < ((SparsePolyZp.front! g_1).snd).prime.toNat :=
                      modInv_val_lt _ _ (by rw [ha_prime]; exact hp.out.one_lt)
                        (by rw [ha_prime]; omega)
                  _ = p := ha_prime
              exact ⟨by simpa [Zp.inv] using ha_prime, h_val_lt⟩
            have h_inv_val_nonzero : ((SparsePolyZp.front! g_1).snd.inv).val.toNat ≠ 0 :=
              Zp_inv_val_nonzero ((SparsePolyZp.front! g_1).snd) h_red_lc h_lc_val_nonzero hp_size
            have h_val_nonzero_g2 : ∀ x ∈ g_2.toList, x.snd.val.toNat ≠ 0 := by
              intro x hx
              set lc_inv := ((SparsePolyZp.front! g_1).snd.inv) with hlc_inv
              -- Convert List membership to Array membership (as in upoly_make_monic_wellFormed)
              have hx_arr : x ∈ g_2 := Array.mem_def.mpr hx
              dsimp [g_2] at hx_arr
              -- hx_arr: x ∈ (Generated.__upoly_make_monic_ir g_1).snd
              simp [Generated.__upoly_make_monic_ir, Generated.__upoly_make_monic_ir_def] at hx_arr
              split_ifs at hx_arr with h
              · -- lc.val == 1: g_2 = g_1, so x ∈ g_1
                simp at hx_arr
                -- hx_arr: x ∈ g_1 (Array membership), convert to List
                have hx_list : x ∈ g_1.toList := (Array.mem_def.mp hx_arr)
                exact h_val_nonzero_g1 x hx_list
              · simp at hx_arr
                -- hx_arr: x ∈ ((loop ...).snd : SparsePolyZp)
                have hx_list : x ∈ (Generated._loop___upoly_make_monic_0_ir 0 g_1 lc_inv).snd.toList :=
                  (show x ∈ ((Generated._loop___upoly_make_monic_0_ir 0 g_1 lc_inv).snd : SparsePolyZp) from hx_arr).val
                have hloop : (Generated._loop___upoly_make_monic_0_ir 0 g_1 lc_inv).snd.toList =
                    (g_1.toList).map (fun (m, y) => (m, y * lc_inv)) := by
                  have := loop_result_toList g_1 lc_inv 0 (Nat.zero_le _)
                  simpa [List.take, List.drop] using this
                rw [hloop] at hx_list
                rcases List.mem_map.mp hx_list with ⟨y, hy, rfl⟩
                exact Zp_mul_val_nonzero (h_red_g1 y hy) h_red_inv (h_val_nonzero_g1 y hy) h_inv_val_nonzero hp_size
            have h_ih_g2 : toPolyList (__squarefree_Zp_ir_safe p g_2) p = sqfZp (SparsePolyZp.toPoly p g_2) :=
              ih (SparsePolyZp.toPoly p g_2).natDegree h_deg_g2_lt_n g_2 rfl
                h_sorted_g2 h_wf_g2 h_red_g2 hp_size hp2 h_nov_g2 h_pdb_g2 h_db_g2 h_sdb_g2
                h_val_nonzero_g2
            let p_1 : UInt64 := (SparsePolyZp.front! g).snd.prime
            have hp_1_eq_p : p_1.toNat = p := by
              have hsize_pos : 0 < Array.size g := by
                by_contra! hzero
                have hzero_sz : Array.size g = 0 := by omega
                have h_empty : g = #[] := Array.eq_empty_of_size_eq_zero hzero_sz
                have : SparsePolyZp.toPoly p g = 0 := by
                  rw [h_empty]; exact SparsePolyZp.toPoly_empty p
                have hzero_nd : (SparsePolyZp.toPoly p g).natDegree = 0 := by simp [this]
                omega
              have hfront_wf : (SparsePolyZp.front! g).snd.prime.toNat = p :=
                hwf_g (SparsePolyZp.front! g) (by
                  apply Array.Mem.mk
                  have hsize_pos' : 0 < Array.size g := hsize_pos
                  simpa [SparsePolyZp.front!] using mem_getFirst_toList g hsize_pos')
              simp [p_1, hfront_wf]
            have h_loop_lemma (arr : Array (SparsePolyZp × UInt64)) :
                (Generated._loop___squarefree_Zp_0_ir_def 0 arr #[] p_1).2.2 = arr.map (fun (g_h, e) => (g_h, e * p_1)) := by
              set f : (SparsePolyZp × UInt64) → (SparsePolyZp × UInt64) := fun (g_h, e) => (g_h, e * p_1)
              have h_wf : WellFounded (λ (a b : Array (SparsePolyZp × UInt64) × Nat) =>
                  a.1.size - a.2 < b.1.size - b.2) :=
                (measure (λ (p : Array (SparsePolyZp × UInt64) × Nat) => p.1.size - p.2)).wf
              have h_toList : (Generated._loop___squarefree_Zp_0_ir_def 0 arr #[] p_1).2.2.toList =
                  (arr.map f).toList := by
                have h_inv := h_wf.induction (arr, 0) (C := λ p => ∀ (result : Array (SparsePolyZp × UInt64)),
                    (Generated._loop___squarefree_Zp_0_ir_def p.2 p.1 result p_1).2.2.toList =
                    result.toList ++ (List.drop p.2 (p.1.toList)).map f) (by
                  intro p ih result
                  unfold Generated._loop___squarefree_Zp_0_ir_def
                  split_ifs with h
                  · -- After Corpus fix, no `let` bindings; LHS = h_ih LHS directly
                    have hx_len : p.2 < p.1.toList.length := by simpa using h
                    have h_measure : Array.size p.1 - (p.2 + 1) < Array.size p.1 - p.2 := by omega
                    have h_ih := ih (p.1, p.2 + 1) h_measure
                      (Array.push result ((p.1[p.2]!).1, (p.1[p.2]!).2 * p_1))
                    rw [h_ih]
                    -- Goal: push.toList ++ (drop (p.2+1)).map f = result.toList ++ (drop p.2).map f
                    -- Step 1: expand push.toList
                    calc
                      (Array.push result ((p.1[p.2]!).1, (p.1[p.2]!).2 * p_1)).toList ++
                          (List.drop (p.2 + 1) (p.1.toList)).map f
                          = (result.toList ++ [((p.1[p.2]!).1, (p.1[p.2]!).2 * p_1)]) ++
                              (List.drop (p.2 + 1) (p.1.toList)).map f := by
                        simp
                      _ = result.toList ++ ([((p.1[p.2]!).1, (p.1[p.2]!).2 * p_1)] ++
                              (List.drop (p.2 + 1) (p.1.toList)).map f) := by simp
                      _ = result.toList ++ ((p.1[p.2]! :: List.drop (p.2 + 1) (p.1.toList)).map f) := by
                        simp [f]
                      _ = result.toList ++ (List.drop p.2 (p.1.toList)).map f := by
                        congr 1
                        -- Prove: (arr[n]! :: drop (n+1)).map f = (drop n).map f
                        -- which follows from: drop n = arr[n]! :: drop (n+1)
                        have h_drop_eq : List.drop p.2 (p.1.toList) = (p.1[p.2]!) :: List.drop (p.2 + 1) (p.1.toList) := by
                          -- Use the general list lemma, then relate List.get to Array.get!
                          have h_temp := drop_eq_get_cons_general (p.1.toList) p.2 hx_len
                          -- h_temp: drop = l.get :: drop
                          -- p.1[p.2]! is the same element as l.get
                          have h_get : (p.1.toList).get ⟨p.2, hx_len⟩ = p.1[p.2]! := by
                            -- Exact same pattern as mem_getFirst_toList (line 673-677)
                            have hpos : p.2 < Array.size p.1 := by
                              simpa [Array.length_toList] using hx_len
                            have h_eq : (p.1.toList).get ⟨p.2, hx_len⟩ = p.1[p.2] := by simp
                            have h_eq' : p.1[p.2]! = p.1[p.2] := by
                              rw [getElem!_def, getElem?_def]
                              simp [hpos]
                            calc
                              (p.1.toList).get ⟨p.2, hx_len⟩ = p.1[p.2] := h_eq
                              _ = p.1[p.2]! := h_eq'.symm
                          rw [h_get] at h_temp
                          exact h_temp
                        rw [h_drop_eq]
                  · -- ¬ p.2 < Array.size p.1: loop returns result
                    dsimp
                    have hlen : p.1.toList.length ≤ p.2 := by
                      simpa using (by omega : Array.size p.1 ≤ p.2)
                    have hdrop : List.drop p.2 (p.1.toList) = [] := drop_all_of_le (p.1.toList) p.2 hlen
                    simp [hdrop]) #[]
                simpa [f] using h_inv
              -- Array 相等由 toList 相等直接得出（Array.ext' : toList 相等 → Array 相等）
              exact Array.ext' h_toList

            have h_loop_eq : (Generated._loop___squarefree_Zp_0_ir_def 0
                (Generated.__squarefree_Zp_ir g_2) #[] p_1).2.2 =
                (Generated.__squarefree_Zp_ir g_2).map (fun (g_h, e) => (g_h, e * p_1)) :=
              h_loop_lemma (Generated.__squarefree_Zp_ir g_2)
            have h_sqfz_g : sqfZp (SparsePolyZp.toPoly p g) = (sqfZp (Polynomial.contract p (SparsePolyZp.toPoly p g))).map
                (fun (h, e) => (h, e * p)) := by
              rw [sqfZp]
              split_ifs with h1
              · exfalso
                have : (SparsePolyZp.toPoly p g).natDegree = n := h_deg_eq
                have hpos : n > 0 := by
                  rw [← this]
                  omega
                omega
              · simp
            have h_toPoly_g2_eq : SparsePolyZp.toPoly p g_2 = SparsePolyZp.toPoly p (SparsePolyZp.makeMonic g_1) :=
              __upoly_make_monic_ir_refines p g_1 h_wf_g1 h_red_g1 hp_size
            have h_makeMonic_eq : SparsePolyZp.toPoly p (SparsePolyZp.makeMonic g_1) =
                Polynomial.C (((SparsePolyZp.front! g_1).snd.inv).toZMod p) * SparsePolyZp.toPoly p g_1 := by
              unfold SparsePolyZp.makeMonic
              split_ifs with hempty
              · have h_empty' : SparsePolyZp.toPoly p g_1 = 0 := by
                  have h_arr_empty : g_1 = #[] :=
                    Array.eq_empty_of_size_eq_zero (by
                      simpa [Array.isEmpty_iff_size_eq_zero] using hempty)
                  simpa [h_arr_empty, SparsePolyZp.toPoly_empty]
                simp [h_empty']
              · have hfront_prime_safe : ((SparsePolyZp.front! g_1).snd.inv).prime.toNat = p := by
                  have hfront_in_g1 : (SparsePolyZp.front! g_1) ∈ g_1.toList := by
                    have hsize_pos_g1 : 0 < Array.size g_1 := by
                      by_contra! hzero_sz
                      have hsz_zero : Array.size g_1 = 0 := by omega
                      have h_empty_arr : g_1 = #[] := Array.eq_empty_of_size_eq_zero hsz_zero
                      have h_not_empty : ¬ g_1.isEmpty := hempty
                      have : g_1.isEmpty := by
                        rw [h_empty_arr]; simp
                      exact h_not_empty this
                    simpa [SparsePolyZp.front!] using mem_getFirst_toList g_1 hsize_pos_g1
                  have hfront_in_g1_arr : (SparsePolyZp.front! g_1) ∈ g_1 := by
                    simpa using hfront_in_g1
                  have hprime_wf : (SparsePolyZp.front! g_1).snd.prime.toNat = p :=
                    h_wf_g1 (SparsePolyZp.front! g_1) hfront_in_g1_arr
                  have hprime_wf_inv : ((SparsePolyZp.front! g_1).snd.inv).prime.toNat = p := by
                    simpa [Zp.inv] using hprime_wf
                  exact hprime_wf_inv
                simpa using toPoly_scalarMul ((SparsePolyZp.front! g_1).snd.inv) g_1 h_red_g1 hfront_prime_safe hp_size
            have h_lc_inv_nonzero : ((SparsePolyZp.front! g_1).snd.inv).toZMod p ≠ 0 := by
              set a := (SparsePolyZp.front! g_1).snd with ha
              have hp_lt_U64 : p < UInt64.size := by
                have : 2 * p ≤ UInt64.size := hp_size; nlinarith
              have ha_prime : a.prime.toNat = p := by
                have hfront_in_arr : (SparsePolyZp.front! g_1) ∈ g_1 :=
                  Array.mem_def.mpr hmem_front_g1
                have hwf_front := h_wf_g1 (SparsePolyZp.front! g_1) hfront_in_arr
                simpa [ha] using hwf_front
              have ha_inv_prime : (a.inv).prime.toNat = p := by
                simpa [Zp.inv, ha] using ha_prime
              have h_red_a : Zp.Reduced p a := by
                have := h_red_g1 (SparsePolyZp.front! g_1) hmem_front_g1
                simpa [ha] using this
              have hval_nonzero : a.val.toNat ≠ 0 :=
                h_val_nonzero_g1 (SparsePolyZp.front! g_1) hmem_front_g1
              have h_mul_one : Zp.toZMod p (a * a.inv) = (1 : ZMod p) :=
                Zp_toZMod_inv_mul_self a h_red_a hval_nonzero hp_lt_U64
              have h_mul_split : Zp.toZMod p (a * a.inv) = Zp.toZMod p a * Zp.toZMod p (a.inv) :=
                Zp.toZMod_mul_weak a a.inv ha_prime ha_inv_prime hp_size
              rw [h_mul_split] at h_mul_one
              intro hzero
              rw [hzero] at h_mul_one
              simp at h_mul_one
            have h_nd_g2_pos : (SparsePolyZp.toPoly p g_2).natDegree > 0 := by
              have h_nd_g1_pos : (SparsePolyZp.toPoly p g_1).natDegree > 0 := by
                rw [h_toPoly_g1]
                exact h_contract_pos
              rw [h_toPoly_g2_eq, h_makeMonic_eq]
              have h_nd_mul : (Polynomial.C (((SparsePolyZp.front! g_1).snd.inv).toZMod p) * SparsePolyZp.toPoly p g_1).natDegree = (SparsePolyZp.toPoly p g_1).natDegree :=
                Polynomial.natDegree_C_mul h_lc_inv_nonzero
              rw [h_nd_mul]
              exact h_nd_g1_pos
            have h_safe_eq : __squarefree_Zp_ir_safe p g_2 = Generated.__squarefree_Zp_ir g_2 := by
              unfold __squarefree_Zp_ir_safe
              simp [h_nd_g2_pos.ne']
            have hg2can : CanonicalRep p g_2 := by
              refine ⟨h_sorted_g2, ?_, h_red_g2⟩
              apply nonzeroB_of_val_nonzero
              intro x hx hv
              exact h_val_nonzero_g2 x hx (by simp [hv])
            have hg2ne : ¬g_2.isEmpty :=
              nonempty_of_toPoly_ne_zero (p := p) g_2 (by
                intro hz
                simp [hz] at h_nd_g2_pos)
            have hmeasure_g2 : Generated.squarefreeMeasure g_2 <
                Generated.squarefreeMeasure g := by
              rw [squarefreeMeasure_eq (p := p) g_2 hg2can,
                squarefreeMeasure_eq (p := p) g hgcan,
                if_neg hg2ne, if_neg hgne]
              omega
            have hmeasure_g2_def := hmeasure_g2
            dsimp [g_2, g_1, Generated.__extract_pth_root_ir,
              Generated.__upoly_make_monic_ir] at hmeasure_g2_def
            have h_result : toPolyList (Generated.__squarefree_Zp_ir g) p = sqfZp (SparsePolyZp.toPoly p g) := by
              have h_step1 : toPolyList (Generated.__squarefree_Zp_ir g) p = toPolyList ((Generated._loop___squarefree_Zp_0_ir_def 0
                  (__squarefree_Zp_ir_safe p g_2) #[] p_1).2.2) p := by
                have h_def_eq : Generated.__squarefree_Zp_ir g = (Generated._loop___squarefree_Zp_0_ir_def 0
                    (Generated.__squarefree_Zp_ir g_2) #[] p_1).2.2 := by
                  dsimp [Generated.__squarefree_Zp_ir]
                  rw [Generated.__squarefree_Zp_ir_def.eq_1 g]
                  split_ifs with h_cond
                  · simp [hmeasure_g2_def, Generated.__squarefree_Zp_ir,
                      Generated.__extract_pth_root_ir,
                      Generated.__upoly_make_monic_ir, g_2, g_1, p_1]
                  · exfalso
                    apply h_cond
                    simpa [h_deriv0]
                calc
                  toPolyList (Generated.__squarefree_Zp_ir g) p
                      = toPolyList ((Generated._loop___squarefree_Zp_0_ir_def 0 (Generated.__squarefree_Zp_ir g_2) #[] p_1).2.2) p := by rw [h_def_eq]
                  _ = toPolyList ((Generated._loop___squarefree_Zp_0_ir_def 0 (__squarefree_Zp_ir_safe p g_2) #[] p_1).2.2) p := by rw [h_safe_eq]
              have h_calc : toPolyList ((Generated._loop___squarefree_Zp_0_ir_def 0
                  (__squarefree_Zp_ir_safe p g_2) #[] p_1).2.2) p = sqfZp (SparsePolyZp.toPoly p g) := by
                rw [← h_safe_eq] at h_loop_eq
                rw [h_loop_eq]
                -- toPolyList (arr.map ...) p = (toPolyList arr p).map ... (need hp_1_eq_p, no overflow bound)
                have h_toPolyList_specific :
                    toPolyList ((__squarefree_Zp_ir_safe p g_2).map (fun (g_h, e) => (g_h, e * p_1))) p =
                    (toPolyList (__squarefree_Zp_ir_safe p g_2) p).map (fun (h, e) => (h, e * p)) := by
                  unfold toPolyList
                  simp only [Array.map_map, Function.comp, Array.toList_map, List.map_map]
                  apply List.map_congr_left
                  intro x hx
                  -- x ∈ (__squarefree_Zp_ir_safe p g_2).toList；证指数界 x.2.toNat * p < 2^64
                  have hmem_poly : (SparsePolyZp.toPoly p x.1, x.2.toNat)
                      ∈ sqfZp (SparsePolyZp.toPoly p g_2) := by
                    rw [← h_ih_g2]; unfold toPolyList
                    rw [Array.toList_map, List.mem_map]; exact ⟨x, hx, rfl⟩
                  have he_le : x.2.toNat ≤ (SparsePolyZp.toPoly p g_2).natDegree :=
                    sqfZp_exponent_le_natDegree _ _ hmem_poly
                  have hg2_deg : (SparsePolyZp.toPoly p g_2).natDegree * p
                      = (SparsePolyZp.toPoly p g).natDegree := by
                    rw [h_toPoly_g2_eq, h_makeMonic_eq,
                      Polynomial.natDegree_C_mul h_lc_inv_nonzero, h_toPoly_g1]; exact h_eq_nd
                  have hn_lt : (SparsePolyZp.toPoly p g).natDegree < 2 ^ 64 :=
                    listSum_natDegree_lt p g.toList (2 ^ 64) (by norm_num) h_deg_bound_g
                  have hbound : x.2.toNat * p < 2 ^ 64 :=
                    lt_of_le_of_lt (le_trans (Nat.mul_le_mul_right p he_le) (le_of_eq hg2_deg)) hn_lt
                  simp only [Function.comp_apply, Prod.mk.injEq, true_and]
                  rw [UInt64.toNat_mul, hp_1_eq_p, Nat.mod_eq_of_lt hbound]
                rw [h_toPolyList_specific]
                rw [h_ih_g2]
                -- sqfZp invariant under scaling by units (makeMonic preserves sqfZp)
                have h_sqfz_smul : sqfZp (SparsePolyZp.toPoly p g_2) = sqfZp (SparsePolyZp.toPoly p g_1) := by
                  have h_g2_eq_C_mul : SparsePolyZp.toPoly p g_2 = Polynomial.C (((SparsePolyZp.front! g_1).snd.inv).toZMod p) * SparsePolyZp.toPoly p g_1 := by
                    rw [h_toPoly_g2_eq, h_makeMonic_eq]
                  rw [h_g2_eq_C_mul]
                  exact sqfZp_smul (c := ((SparsePolyZp.front! g_1).snd.inv).toZMod p) (hc := h_lc_inv_nonzero) (f := SparsePolyZp.toPoly p g_1)
                rw [h_sqfz_smul]
                rw [h_toPoly_g1]
                rw [h_sqfz_g]
              rw [h_step1, h_calc]
            -- Branch A 的结果
            simp [h_deg0_g]
            exact h_result
          · -- Branch B: derivative ≠ 0 → Yun algorithm
            have hg_nonzeroB : SparsePolyZp.NonZeroB g := by
              apply nonzeroB_of_val_nonzero
              intro x hx hv
              exact h_val_nonzero_g x hx (by simp [hv])
            have hgcan : CanonicalRep p g := ⟨hsorted_g, hg_nonzeroB, hred_g⟩
            have hgpolyne : SparsePolyZp.toPoly p g ≠ 0 := by
              intro hz
              simp [hz] at h_deg0_g
            have hgne : ¬g.isEmpty := nonempty_of_toPoly_ne_zero (p := p) g hgpolyne
            let fd := SparsePolyZp.derivative g
            have hfdcan : CanonicalRep p fd := by
              simpa [fd] using derivative_canonical (p := p) g hgcan hgne
            have hfdne : ¬fd.isEmpty := by simpa [fd] using h_deriv0
            have hfdpolyne : SparsePolyZp.toPoly p fd ≠ 0 :=
              toPoly_ne_zero_of_canonical_nonempty (p := p) fd hfdcan hfdne
            have hfd_signed : ∀ x ∈ fd.toList, x.1.deg < 2 ^ 63 := by
              simpa [fd] using derivative_deg_bound_of (2 ^ 63) g h_signed_deg_bound_g
            let c1 := polynomial_GCD g fd
            have hc1data := polynomial_GCD_step_data (p := p) hp_size hp2 g fd hgcan hfdcan hgpolyne
            have hc1can : CanonicalRep p c1 := by simpa [c1] using hc1data.1
            have hc1ne : ¬c1.isEmpty := by simpa [c1] using hc1data.2.1
            have hc1poly : SparsePolyZp.toPoly p c1 = normalize
                (EuclideanDomain.gcd (SparsePolyZp.toPoly p g)
                  (SparsePolyZp.toPoly p fd)) := by
              simpa [c1] using hc1data.2.2.1
            have hc1monic : Monic (SparsePolyZp.toPoly p c1) := by
              simpa [c1] using hc1data.2.2.2
            have hc1polyne : SparsePolyZp.toPoly p c1 ≠ 0 :=
              toPoly_ne_zero_of_canonical_nonempty (p := p) c1 hc1can hc1ne
            have hc1dvd : SparsePolyZp.toPoly p c1 ∣ SparsePolyZp.toPoly p g := by
              rw [hc1poly]
              exact normalize_dvd_iff.mpr (EuclideanDomain.gcd_dvd_left _ _)
            let wraw := (SparsePolyZp.divmod g c1).1
            have hwrawcan : CanonicalRep p wraw := by
              exact divmod_fst_canonical (p := p) g c1 hgcan hc1can hc1ne hp_size
            have hwrawpoly : SparsePolyZp.toPoly p wraw =
                SparsePolyZp.toPoly p g /ₘ SparsePolyZp.toPoly p c1 := by
              simpa [wraw] using divmod_fst_toPoly_eq_divByMonic (p := p) hp_size hp2
                g c1 hgcan hc1can hc1ne hc1monic hc1dvd
            have hwrawpolyne : SparsePolyZp.toPoly p wraw ≠ 0 := by
              rw [hwrawpoly]
              exact divByMonic_ne_zero_local _ _ hc1monic hc1dvd hgpolyne
            have hwrawne : ¬wraw.isEmpty :=
              nonempty_of_toPoly_ne_zero (p := p) wraw hwrawpolyne
            have hnormw : SparsePolyZp.normalization wraw = wraw :=
              normalization_eq_of_nonZeroB wraw hwrawcan.2.1
            let w3 := SparsePolyZp.normalization wraw
            have hw3can : CanonicalRep p w3 := normalization_canonical wraw hwrawcan
            have hw3poly : SparsePolyZp.toPoly p w3 = SparsePolyZp.toPoly p wraw := by
              simpa [w3] using normalization_toPoly (p := p) wraw
            have hgnd63 : (SparsePolyZp.toPoly p g).natDegree < 2 ^ 63 := by
              rw [(toPoly_head_data (p := p) g hgcan hgne).1]
              exact h_signed_deg_bound_g g[0]! (mem_getFirst_toList g (by
                rw [Array.isEmpty_iff_size_eq_zero] at hgne
                omega))
            have hc1ndle : (SparsePolyZp.toPoly p c1).natDegree ≤
                (SparsePolyZp.toPoly p g).natDegree :=
              Polynomial.natDegree_le_of_dvd hc1dvd hgpolyne
            have hwrawndle : (SparsePolyZp.toPoly p wraw).natDegree ≤
                (SparsePolyZp.toPoly p g).natDegree := by
              rw [hwrawpoly, Polynomial.natDegree_divByMonic _ hc1monic]
              omega
            have hc1signed : ∀ x ∈ c1.toList, x.1.deg < 2 ^ 63 :=
              canonical_deg_bound_of_natDegree_lt (p := p) c1 hc1can
                (lt_of_le_of_lt hc1ndle hgnd63)
            have hwrawsigned : ∀ x ∈ wraw.toList, x.1.deg < 2 ^ 63 :=
              canonical_deg_bound_of_natDegree_lt (p := p) wraw hwrawcan
                (lt_of_le_of_lt hwrawndle hgnd63)
            have hw3signed : ∀ x ∈ w3.toList, x.1.deg < 2 ^ 63 := by
              simpa [w3, hnormw] using hwrawsigned
            have hmul : SparsePolyZp.toPoly p c1 * SparsePolyZp.toPoly p wraw =
                SparsePolyZp.toPoly p g := by
              have hid := modByMonic_add_div (SparsePolyZp.toPoly p g)
                (SparsePolyZp.toPoly p c1)
              rw [(modByMonic_eq_zero_iff_dvd hc1monic).mpr hc1dvd, zero_add] at hid
              simpa [hwrawpoly] using hid
            have hdegsum : (SparsePolyZp.toPoly p w3).natDegree +
                (SparsePolyZp.toPoly p c1).natDegree =
                (SparsePolyZp.toPoly p g).natDegree := by
              have hm := Polynomial.natDegree_mul hc1polyne hwrawpolyne
              rw [hmul] at hm
              rw [hw3poly]
              omega
            have hloopbound : (1 : UInt64).toNat +
                (SparsePolyZp.toPoly p w3).natDegree +
                (SparsePolyZp.toPoly p c1).natDegree + 2 < 2 ^ 64 := by
              simp only [UInt64.toNat_ofNat]
              omega
            have hloop := yunLoop_ir_refines (p := p) hp_size hp2 (1 : UInt64) w3 c1 #[]
              hw3can hc1can hw3signed hc1signed hc1polyne hc1monic hloopbound
            have hfdpoly_eq : SparsePolyZp.toPoly p fd =
                Polynomial.derivative (SparsePolyZp.toPoly p g) := by
              simpa [fd] using (derivative_toPoly_eq g hwf_g hp_size
                h_no_overflow_g h_deg_bound_g)
            have hc1poly_math : SparsePolyZp.toPoly p c1 = normalize
                (EuclideanDomain.gcd (SparsePolyZp.toPoly p g)
                  (Polynomial.derivative (SparsePolyZp.toPoly p g))) := by
              simpa [hfdpoly_eq] using hc1poly
            let loopout := Generated._loop___squarefree_Zp_1_ir_def
              (1 : UInt64) w3 c1 #[]
            let c2 := loopout.2.1
            let result4 := loopout.2.2
            have hc2can : CanonicalRep p c2 := by
              simpa [loopout, c2] using hloop.1
            have hc2signed : ∀ x ∈ c2.toList, x.1.deg < 2 ^ 63 := by
              simpa [loopout, c2] using hloop.2.1
            have hc2polyne : SparsePolyZp.toPoly p c2 ≠ 0 := by
              simpa [loopout, c2] using hloop.2.2.1
            have hc2monic : Monic (SparsePolyZp.toPoly p c2) := by
              simpa [loopout, c2] using hloop.2.2.2.1
            have hloop_eq :
                (toPolyList result4 p, SparsePolyZp.toPoly p c2) =
                  yunLoop (SparsePolyZp.toPoly p w3)
                    (SparsePolyZp.toPoly p c1) 1 [] hc1polyne := by
              simpa [loopout, c2, result4, toPolyList_empty] using hloop.2.2.2.2
            have hloop_eq_norm :
                (toPolyList result4 p, SparsePolyZp.toPoly p c2) =
                  yunLoop (normalize (SparsePolyZp.toPoly p w3))
                    (SparsePolyZp.toPoly p c1) 1 [] hc1polyne := by
              rw [yunLoop_normalize_left]
              exact hloop_eq
            have hnorm_w3 : normalize (SparsePolyZp.toPoly p w3) =
                normalize (SparsePolyZp.toPoly p g /ₘ
                  normalize (EuclideanDomain.gcd (SparsePolyZp.toPoly p g)
                    (Polynomial.derivative (SparsePolyZp.toPoly p g)))) := by
              rw [hw3poly, hwrawpoly, hc1poly_math]
            have hc_math_ne : normalize
                (EuclideanDomain.gcd (SparsePolyZp.toPoly p g)
                  (Polynomial.derivative (SparsePolyZp.toPoly p g))) ≠ 0 := by
              rw [← hc1poly_math]
              exact hc1polyne
            have hloop_eq_math :
                (toPolyList result4 p, SparsePolyZp.toPoly p c2) =
                  yunLoop
                    (normalize (SparsePolyZp.toPoly p g /ₘ
                      normalize (EuclideanDomain.gcd (SparsePolyZp.toPoly p g)
                        (Polynomial.derivative (SparsePolyZp.toPoly p g)))))
                    (normalize (EuclideanDomain.gcd (SparsePolyZp.toPoly p g)
                      (Polynomial.derivative (SparsePolyZp.toPoly p g))))
                    1 [] hc_math_ne := by
              simpa [hnorm_w3, hc1poly_math] using hloop_eq_norm
            have hcrem_deriv : Polynomial.derivative (SparsePolyZp.toPoly p c2) = 0 := by
              have hm := yunLoop_sqf_remainder_derivative_zero
                (p := p) (SparsePolyZp.toPoly p g) hgpolyne hc_math_ne
              have hc2eq := congrArg Prod.snd hloop_eq_math
              simp only [Prod.snd] at hc2eq
              rw [hc2eq]
              exact hm
            have hc2ne : ¬c2.isEmpty :=
              nonempty_of_toPoly_ne_zero (p := p) c2 hc2polyne
            have hc2ndle : (SparsePolyZp.toPoly p c2).natDegree ≤
                (SparsePolyZp.toPoly p g).natDegree := by
              have hrem := congrArg (fun z => z.2.natDegree) hloop_eq_math
              dsimp only at hrem
              rw [hrem]
              exact le_trans (yunLoop_c_natDegree_le _ _ 1 [] hc_math_ne)
                (Polynomial.natDegree_le_of_dvd
                  (normalize_dvd_iff.mpr (EuclideanDomain.gcd_dvd_left _ _)) hgpolyne)
            have hc2pdeg : ∀ x ∈ c2.toList, p * x.1.deg < 2 ^ 64 := by
              intro x hx
              have hxle : x.1.deg ≤ (SparsePolyZp.toPoly p c2).natDegree := by
                have hlt := canonical_deg_bound_of_natDegree_lt (p := p) c2 hc2can
                  (Nat.lt_succ_self (SparsePolyZp.toPoly p c2).natDegree) x hx
                omega
              have hghead_p : p * (SparsePolyZp.toPoly p g).natDegree < 2 ^ 64 := by
                rw [(toPoly_head_data (p := p) g hgcan hgne).1]
                exact h_p_deg_bound_g g[0]! (mem_getFirst_toList g (by
                  rw [Array.isEmpty_iff_size_eq_zero] at hgne
                  omega))
              exact lt_of_le_of_lt (Nat.mul_le_mul_left p (le_trans hxle hc2ndle)) hghead_p
            have hc2deg64 : ∀ x ∈ c2.toList, x.1.deg < 2 ^ 64 := by
              intro x hx
              exact lt_trans (hc2signed x hx) (by norm_num)
            have hc2nov : ∀ x ∈ c2.toList,
                x.2.val.toNat * x.1.deg < 2 ^ 64 := by
              intro x hx
              exact lt_of_le_of_lt (Nat.mul_le_mul
                (Nat.le_of_lt (hc2can.2.2 x hx).2) (Nat.le_refl _)) (hc2pdeg x hx)
            have hc2derivpoly : SparsePolyZp.toPoly p (SparsePolyZp.derivative c2) = 0 := by
              rw [derivative_toPoly_eq c2 (wellFormed_of_canonical c2 hc2can)
                hp_size hc2nov hc2deg64]
              exact hcrem_deriv
            have hc2derivcan : CanonicalRep p (SparsePolyZp.derivative c2) :=
              derivative_canonical (p := p) c2 hc2can hc2ne
            have hc2deriv0 : SparsePolyZp.derivative c2 = SparsePolyZp.empty := by
              by_contra hn
              have hnempty : ¬(SparsePolyZp.derivative c2).isEmpty := by
                intro he
                apply hn
                exact Array.eq_empty_of_size_eq_zero
                  (Array.isEmpty_iff_size_eq_zero.mp he)
              exact (toPoly_ne_zero_of_canonical_nonempty (p := p)
                (SparsePolyZp.derivative c2) hc2derivcan hnempty) hc2derivpoly
            let g3 := Generated.__extract_pth_root_ir c2
            let g4 := (Generated.__upoly_make_monic_ir g3).snd
            have hg3poly : SparsePolyZp.toPoly p g3 =
                Polynomial.contract p (SparsePolyZp.toPoly p c2) := by
              simpa [g3] using extract_pth_root_toPoly_eq c2 hc2deriv0
                (wellFormed_of_canonical c2 hc2can) hc2can.2.2 hp_size hc2nov hc2deg64
            have hg3can : CanonicalRep p g3 := by
              simpa [g3] using extract_pth_root_canonical (p := p) c2 hc2can hc2ne
                hc2deriv0 hc2nov hc2deg64
            have hg3polyne : SparsePolyZp.toPoly p g3 ≠ 0 := by
              rw [hg3poly]
              intro hz
              have hexp := Polynomial.expand_contract p hcrem_deriv hp.out.ne_zero
              rw [hz, map_zero] at hexp
              exact hc2polyne hexp.symm
            have hg3ne : ¬g3.isEmpty :=
              nonempty_of_toPoly_ne_zero (p := p) g3 hg3polyne
            have hg4poly : SparsePolyZp.toPoly p g4 =
                normalize (SparsePolyZp.toPoly p g3) := by
              calc
                SparsePolyZp.toPoly p g4 =
                    SparsePolyZp.toPoly p (SparsePolyZp.makeMonic g3) := by
                  simpa [g4] using __upoly_make_monic_ir_refines p g3
                    (wellFormed_of_canonical g3 hg3can) hg3can.2.2 hp_size
                _ = normalize (SparsePolyZp.toPoly p g3) :=
                  makeMonic_toPoly_eq_normalize (p := p) hp_size g3 hg3can hg3ne
            have hg4polyne : SparsePolyZp.toPoly p g4 ≠ 0 := by
              rw [hg4poly, normalize_ne_zero_iff_local]
              exact hg3polyne
            have hg3deg64 := extract_pth_root_deg_bound c2 hc2deg64
            have hg3signed := extract_pth_root_signed_deg_bound c2 hc2signed
            have hg3val : ∀ x ∈ g3.toList, x.2.val.toNat ≠ 0 := by
              simpa [g3] using extract_pth_root_val_nonzero c2 (fun x hx => by
                exact fun hz => (val_nonzero_of_nonzeroB c2 hc2can.2.1 x hx)
                  (UInt64.toNat_inj.mp (by simpa using hz)))
            have hg3wf := extract_pth_root_wellFormed (p := p) c2
              (wellFormed_of_canonical c2 hc2can)
            have hg3red := extract_pth_root_allReduced (p := p) c2 hc2can.2.2
            have hg3ov := extract_pth_root_no_overflow (p := p) c2
              (wellFormed_of_canonical c2 hc2can) hc2nov hc2deg64
            have hg4sorted : SparsePolyZp.Sorted g4 := upoly_make_monic_sorted g3 hg3can.1
            have hg4wf : SparsePolyZp.WellFormed p g4 := upoly_make_monic_wellFormed g3 hg3wf
            have hg4red : SparsePolyZp.AllReduced p g4.toList :=
              upoly_make_monic_allReduced g3 hg3red hp_size
            have hg4deg64 := upoly_make_monic_deg_bound g3 hg3deg64
            have hg4signed := upoly_make_monic_deg_bound_of (2 ^ 63) g3 hg3signed
            have hg4ov := upoly_make_monic_no_overflow g3 hg3red hg3ov.1 hg3ov.2 hp_size
            have hg4pdeg := upoly_make_monic_p_deg_bound (p := p) g3 hg3ov.2
            have hg4val := upoly_make_monic_val_nonzero (p := p) g3 hg3red hg3val hg3ne hp_size
            have hc2_contract_deg :
                (Polynomial.contract p (SparsePolyZp.toPoly p c2)).natDegree * p =
                  (SparsePolyZp.toPoly p c2).natDegree := by
              have hexp := Polynomial.expand_contract p hcrem_deriv hp.out.ne_zero
              calc
                (Polynomial.contract p (SparsePolyZp.toPoly p c2)).natDegree * p =
                    (Polynomial.expand (ZMod p) p
                      (Polynomial.contract p (SparsePolyZp.toPoly p c2))).natDegree := by
                  rw [Polynomial.natDegree_expand, mul_comm]
                _ = (SparsePolyZp.toPoly p c2).natDegree := by rw [hexp]
            have hg4deg_lt : (SparsePolyZp.toPoly p g4).natDegree < n := by
              rw [hg4poly, natDegree_normalize_eq_local, hg3poly]
              rw [← h_deg_eq]
              by_cases hz : (Polynomial.contract p
                  (SparsePolyZp.toPoly p c2)).natDegree = 0
              · rw [hz]
                omega
              · have hlt : (Polynomial.contract p
                    (SparsePolyZp.toPoly p c2)).natDegree <
                    (Polynomial.contract p
                      (SparsePolyZp.toPoly p c2)).natDegree * p :=
                  by
                    have hm := Nat.mul_lt_mul_of_pos_left hp.out.one_lt
                      (Nat.pos_of_ne_zero hz)
                    simpa using hm
                rw [hc2_contract_deg] at hlt
                omega
            have hih_g4 : toPolyList (__squarefree_Zp_ir_safe p g4) p =
                sqfZp (SparsePolyZp.toPoly p g4) :=
              ih (SparsePolyZp.toPoly p g4).natDegree hg4deg_lt g4 rfl
                hg4sorted hg4wf hg4red hp_size hp2 hg4ov hg4pdeg hg4deg64
                hg4signed hg4val
            have hpolyderiv_ne : Polynomial.derivative (SparsePolyZp.toPoly p g) ≠ 0 := by
              rw [← hfdpoly_eq]
              exact hfdpolyne
            have himpl_cond_iff :
                ((!c2.isEmpty) && decide (get_deg c2 > (0 : Int64))) = true ↔
                  0 < (SparsePolyZp.toPoly p c2).natDegree :=
              squarefree_loop_cond_iff (p := p) c2 hc2can hc2signed
            by_cases hc2pos : 0 < (SparsePolyZp.toPoly p c2).natDegree
            · have himpl_true : ((!c2.isEmpty) &&
                  decide (get_deg c2 > (0 : Int64))) = true :=
                himpl_cond_iff.mpr hc2pos
              have hg4pos : 0 < (SparsePolyZp.toPoly p g4).natDegree := by
                rw [hg4poly, natDegree_normalize_eq_local, hg3poly]
                by_contra hz
                have hz0 : (Polynomial.contract p
                    (SparsePolyZp.toPoly p c2)).natDegree = 0 := by omega
                rw [hz0, zero_mul] at hc2_contract_deg
                omega
              have hsafe_g4 : __squarefree_Zp_ir_safe p g4 =
                  Generated.__squarefree_Zp_ir g4 := by
                unfold __squarefree_Zp_ir_safe
                rw [if_neg hg4pos.ne']
              have hg4can : CanonicalRep p g4 := by
                refine ⟨hg4sorted, ?_, hg4red⟩
                apply nonzeroB_of_val_nonzero
                intro x hx hv
                exact hg4val x hx (by simp [hv])
              have hg4ne : ¬g4.isEmpty :=
                nonempty_of_toPoly_ne_zero (p := p) g4 hg4polyne
              have hmeasure_g4 : Generated.squarefreeMeasure g4 <
                  Generated.squarefreeMeasure g := by
                rw [squarefreeMeasure_eq (p := p) g4 hg4can,
                  squarefreeMeasure_eq (p := p) g hgcan,
                  if_neg hg4ne, if_neg hgne]
                omega
              let p_1 : UInt64 := (SparsePolyZp.front! g).snd.prime
              have hp1 : p_1.toNat = p := by
                simpa [p_1] using (hgcan.2.2 g[0]!
                  (mem_getFirst_toList g (by
                    rw [Array.isEmpty_iff_size_eq_zero] at hgne
                    omega))).1
              let loop2out := Generated._loop___squarefree_Zp_2_ir_def 0
                (Generated.__squarefree_Zp_ir g4) result4 p_1
              have hloop2list : loop2out.2.2.toList = result4.toList ++
                  (Generated.__squarefree_Zp_ir g4).toList.map
                    (fun (q, e) => (q, e * p_1)) := by
                simpa [loop2out] using squarefree_loop2_result
                  (Generated.__squarefree_Zp_ir g4) result4 p_1
              have hg4deg_mul : (SparsePolyZp.toPoly p g4).natDegree * p =
                  (SparsePolyZp.toPoly p c2).natDegree := by
                rw [hg4poly, natDegree_normalize_eq_local, hg3poly]
                exact hc2_contract_deg
              have hloop2poly : toPolyList loop2out.2.2 p =
                  toPolyList result4 p ++
                    (toPolyList (__squarefree_Zp_ir_safe p g4) p).map
                      (fun (q, e) => (q, e * p)) := by
                unfold toPolyList
                simp only [Array.toList_map]
                rw [hloop2list, List.map_append, List.map_map, hsafe_g4]
                congr 1
                simp only [List.map_map]
                apply List.map_congr_left
                intro x hx
                have hxpoly : (SparsePolyZp.toPoly p x.1, x.2.toNat) ∈
                    sqfZp (SparsePolyZp.toPoly p g4) := by
                  rw [← hih_g4, hsafe_g4]
                  unfold toPolyList
                  rw [Array.toList_map]
                  exact List.mem_map.mpr ⟨x, hx, rfl⟩
                have hexp := sqfZp_exponent_le_natDegree _ _ hxpoly
                have hbound : x.2.toNat * p < 2 ^ 64 := by
                  have hle : x.2.toNat * p ≤
                      (SparsePolyZp.toPoly p g4).natDegree * p :=
                    Nat.mul_le_mul_right p hexp
                  rw [hg4deg_mul] at hle
                  have hc2nd63 : (SparsePolyZp.toPoly p c2).natDegree < 2 ^ 63 :=
                    lt_of_le_of_lt hc2ndle hgnd63
                  exact lt_of_le_of_lt hle (lt_trans hc2nd63 (by norm_num))
                simp only [Function.comp_apply, Prod.mk.injEq, true_and]
                rw [UInt64.toNat_mul, hp1, Nat.mod_eq_of_lt hbound]
              have hsqf_g4 : sqfZp (SparsePolyZp.toPoly p g4) =
                  sqfZp (Polynomial.contract p (SparsePolyZp.toPoly p c2)) := by
                rw [hg4poly, hg3poly, normalize_apply, mul_comm]
                exact sqfZp_smul
                  (c := (↑(normUnit (Polynomial.contract p
                    (SparsePolyZp.toPoly p c2)).leadingCoeff) : ZMod p))
                  (by simp) _
              have hgenerated : Generated.__squarefree_Zp_ir g = loop2out.2.2 := by
                have hderiv_not_empty : ¬Array.isEmpty (_root_.derivative g) = true := by
                  intro he
                  apply h_deriv0
                  simpa using Array.eq_empty_of_size_eq_zero
                    (Array.isEmpty_iff_size_eq_zero.mp he)
                dsimp [Generated.__squarefree_Zp_ir]
                rw [Generated.__squarefree_Zp_ir_def.eq_1, if_neg hderiv_not_empty]
                dsimp only
                have hdivinst : HasPolyDivmod.polyDivmod (α := SparsePolyZp) =
                    SparsePolyZp.divmod := rfl
                have hguard : Generated.squarefreeMeasure
                    (Generated.__upoly_make_monic_ir_def
                      (Generated.__extract_pth_root_ir_def
                        (Generated._loop___squarefree_Zp_1_ir_def 1
                          (pair_vec_div SparsePolyZp.empty g
                            (polynomial_GCD g (_root_.derivative g)) g.comp).normalization
                          (polynomial_GCD g (_root_.derivative g)) #[]).2.1)).2 <
                    Generated.squarefreeMeasure g := by
                  simpa [g4, g3, c2, loopout, w3, wraw, c1, fd,
                    pair_vec_div, hdivinst, Generated.__extract_pth_root_ir,
                    Generated.__upoly_make_monic_ir] using hmeasure_g4
                rw [dif_pos hguard]
                change (if ((!c2.isEmpty) && decide (get_deg c2 > (0 : Int64))) then
                  loop2out.2.2 else result4) = loop2out.2.2
                rw [if_pos himpl_true]
              have hmath : sqfZp (SparsePolyZp.toPoly p g) =
                  toPolyList loop2out.2.2 p := by
                rw [sqfZp, dif_neg h_deg0_g, dif_neg hpolyderiv_ne]
                dsimp only
                simp only [← hloop_eq_math]
                rw [if_pos hc2pos]
                rw [hloop2poly, hih_g4, hsqf_g4]
              simp [h_deg0_g, hgenerated, hmath]
            · have himpl_false : ¬(((!c2.isEmpty) &&
                  decide (get_deg c2 > (0 : Int64))) = true) :=
                fun h => hc2pos (himpl_cond_iff.mp h)
              have hprop_false : ¬(¬c2 = #[] ∧ 0 < get_deg c2) := by
                intro h
                exact hc2pos ((get_deg_pos_iff (p := p) c2 hc2can hc2signed).mp h.2)
              have hgenerated : Generated.__squarefree_Zp_ir g = result4 := by
                have hderiv_not_empty : ¬Array.isEmpty (_root_.derivative g) = true := by
                  intro he
                  apply h_deriv0
                  simpa using Array.eq_empty_of_size_eq_zero
                    (Array.isEmpty_iff_size_eq_zero.mp he)
                dsimp [Generated.__squarefree_Zp_ir]
                rw [Generated.__squarefree_Zp_ir_def.eq_1]
                rw [if_neg hderiv_not_empty]
                dsimp only
                change (if ((!c2.isEmpty) && decide (get_deg c2 > (0 : Int64))) then
                  _ else result4) = result4
                rw [if_neg himpl_false]
              have hmath : sqfZp (SparsePolyZp.toPoly p g) = toPolyList result4 p := by
                rw [sqfZp, dif_neg h_deg0_g, dif_neg hpolyderiv_ne]
                dsimp only
                simp only [← hloop_eq_math]
                simp [hc2pos]
              simp [h_deg0_g, hgenerated, hmath]
    unfold __squarefree_Zp_ir_safe at h_main
    simpa using h_main

end Refinement
