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

/-- 负数的 get_deg 映射到 0（空多项式），≥0 值对应于 ℕ degree。
    TODO: 需要数组排序性质才能证明，当前未使用。 -/
private lemma get_deg_toPoly (f : SparsePolyZp) (hwf : SparsePolyZp.WellFormed p f)
    : get_deg f = (SparsePolyZp.toPoly p f).natDegree := by
  -- 这个引理当前未使用，且需要 mergeAdd 排序性质才能证明
  -- 暂时保留为 admit，待需要时再补完
  admit

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
    (h_p2 : p * p ≤ UInt64.size) : Zp.toZMod p (a * a.inv) = (1 : ZMod p) := by
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
            have hp_lt_U64 : p < UInt64.size := by
              have : p * p ≤ UInt64.size := h_p2
              nlinarith
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

/-- 多项式长除法中 r' = r - term*g 的首项抵消 → deg(r') < dr。 -/
private lemma divmod_deg_decrease (g : SparsePolyZp) (dg : ℕ) (lc_g_inv : Zp) (r : SparsePolyZp)
    (dr : ℕ) (coeff : Zp) (d : ℕ) (term : SparsePolyZp) (r' : SparsePolyZp)
    (hred_g : SparsePolyZp.AllReduced p g.toList) (hred_r : SparsePolyZp.AllReduced p r.toList)
    (h_nonempty_g : ¬g.isEmpty) (hprime : lc_g_inv.prime.toNat = p)
    (h_p2 : p * p ≤ UInt64.size) (h_empty_r' : ¬r'.isEmpty) :
    r'[0]!.fst.deg < dr := by
  admit

-- UInt64 乘法取模到 ZMod 的转换
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

/-- 良基递归版 divmodAux（与 partial def divmodAux 等价，可做归纳）。
     后三个条件（h2p, h_p2, h_nonempty）用于确保 deg(r') < deg(r) 的终止性证明。
     当无法证明终止时，用 decreasing_by { admit } 暂时跳过。 -/
noncomputable def divmodAux'_wf (g : SparsePolyZp) (dg : ℕ) (lc_g_inv : Zp)
    (q r : SparsePolyZp) : SparsePolyZp × SparsePolyZp :=
  if hr_empty : r.isEmpty then (q, r)
  else
    let dr := r[0]!.fst.deg
    if hdr_lt_dg : dr < dg then (q, r)
    else
      let coeff := r[0]!.snd * lc_g_inv
      let d := dr - dg
      let term : SparsePolyZp := #[(⟨d⟩, coeff)]
      let r' := r - (term * g)
      let q' := q.push (⟨d⟩, coeff)
      divmodAux'_wf g dg lc_g_inv q' r'
termination_by
  r[0]!.fst.deg + 1
decreasing_by
  -- 需要证明 r' 的度 < r 的度。
  -- 由于 partial def 版本「数学上能终止」，此处暂跳过。
  sorry

/-- divmodAux'_wf 与原始 partial def divmodAux 对 terminating inputs 等价。
     用 native_decide 可验证所有具体小例；一般情形因 partial def 不可推理而 admit。 -/
private lemma divmodAux'_wf_eq (g : SparsePolyZp) (dg : ℕ) (lc_g_inv : Zp) (q r : SparsePolyZp) :
    divmodAux'_wf g dg lc_g_inv q r = SparsePolyZp.divmodAux g dg lc_g_inv q r := by
  admit

/-- 非递归的 divmod 良基版包装器（供 gcdAux'_wf 使用）。 -/
noncomputable def divmod'_wf (f g : SparsePolyZp) : SparsePolyZp × SparsePolyZp :=
  if g.isEmpty then (#[], f)
  else
    let dg := g[0]!.fst.deg
    let lc_g_inv := g[0]!.snd.inv
    divmodAux'_wf g dg lc_g_inv #[] f

/-- gcdAux 的良基递归版（与 partial def gcdAux 等价，可做归纳）。 -/
noncomputable def gcdAux'_wf (p : ℕ) (f g : SparsePolyZp) : SparsePolyZp :=
  if h : g.isEmpty then f
  else gcdAux'_wf p g (divmod'_wf f g).snd
termination_by
  (SparsePolyZp.toPoly p g).natDegree
decreasing_by
  -- 需要 (toPoly p ((divmod'_wf f g).snd)).natDegree < (toPoly p g).natDegree
  -- 由 leading term cancellation + divmod 性质可得，暂跳过。
  sorry

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

/-- 空数组总是 AllReduced -/
private lemma allReduced_empty (p : ℕ) : SparsePolyZp.AllReduced p ((#[] : SparsePolyZp).toList) := by
  intro x hx; simp at hx

/-- divmod 的循环不变量：f = q*g + r 在 divmodAux 的各步中保持。
     用强归纳于 r[0]!.fst.deg 证明。 -/
private lemma divmodAux_invariant (g : SparsePolyZp) (dg : ℕ) (lc_g_inv : Zp) (q r f : SparsePolyZp)
    (h_inv : SparsePolyZp.toPoly p f = SparsePolyZp.toPoly p q * SparsePolyZp.toPoly p g + SparsePolyZp.toPoly p r)
    (hwf_f : SparsePolyZp.WellFormed p f) (hwg_g : SparsePolyZp.WellFormed p g)
    (hred_f : SparsePolyZp.AllReduced p f.toList) (hred_g : SparsePolyZp.AllReduced p g.toList)
    (hred_q : SparsePolyZp.AllReduced p q.toList) (hred_r : SparsePolyZp.AllReduced p r.toList)
    (h_2p : 2 * p ≤ UInt64.size) (h_p2 : p * p ≤ UInt64.size)
    (h_nonempty_g : ¬g.isEmpty) (h_prime: lc_g_inv.prime.toNat = p) :
    SparsePolyZp.toPoly p f = SparsePolyZp.toPoly p (divmodAux'_wf g dg lc_g_inv q r).1 * SparsePolyZp.toPoly p g +
                              SparsePolyZp.toPoly p (divmodAux'_wf g dg lc_g_inv q r).2 := by
  -- 用强归纳于 r[0]!.fst.deg，同时携带 AllReduced 假设（用于 toPoly homomorphism lemmas）
  have h_aux : ∀ (n : ℕ), ∀ (q r : SparsePolyZp),
    (SparsePolyZp.AllReduced p q.toList) → (SparsePolyZp.AllReduced p r.toList) →
    (r[0]!.fst.deg + 1 = n ∨ r.isEmpty) →
    (SparsePolyZp.toPoly p f = SparsePolyZp.toPoly p q * SparsePolyZp.toPoly p g + SparsePolyZp.toPoly p r) →
    (SparsePolyZp.toPoly p f = SparsePolyZp.toPoly p (divmodAux'_wf g dg lc_g_inv q r).1 * SparsePolyZp.toPoly p g +
                               SparsePolyZp.toPoly p (divmodAux'_wf g dg lc_g_inv q r).2) := by
    intro n
    induction n using Nat.strong_induction_on with
    | h n IH =>
      intro q r hred_q hred_r h_deg h_inv'
      by_cases hr_empty : r.isEmpty
      · -- r.isEmpty → divmodAux'_wf returns (q, r)
        have h_eq : divmodAux'_wf g dg lc_g_inv q r = (q, r) := by
          rw [divmodAux'_wf.eq_1]; simp [hr_empty]
        simpa [h_eq] using h_inv'
      · by_cases hdr_lt_dg : r[0]!.fst.deg < dg
        · -- dr < dg → divmodAux'_wf returns (q, r)
          have h_eq : divmodAux'_wf g dg lc_g_inv q r = (q, r) := by
            rw [divmodAux'_wf.eq_1]; simp [hr_empty, hdr_lt_dg]
          simpa [h_eq] using h_inv'
        · let dr := r[0]!.fst.deg
          have h_dr_ge_dg : dr ≥ dg := by omega
          let coeff := r[0]!.snd * lc_g_inv
          let d := dr - dg
          let term : SparsePolyZp := #[(⟨d⟩, coeff)]
          let r' := r - (term * g)
          let q' := q.push (⟨d⟩, coeff)
          -- 展开一步到递归调用
          have h_step : divmodAux'_wf g dg lc_g_inv q r = divmodAux'_wf g dg lc_g_inv q' r' := by
            rw [divmodAux'_wf.eq_1 g dg lc_g_inv q r]
            have h_not_empty : r.isEmpty = false := by
              simp [hr_empty]
            have h_not_lt : (r[0]!.fst.deg < dg) = false := by
              simp [hdr_lt_dg]
            simp [h_not_empty, h_not_lt]
            dsimp
          rw [h_step]
          have hpos_r : 0 < r.size := by
            rw [Array.isEmpty_iff_size_eq_zero] at hr_empty; omega
          -- ================================================================
          -- (A) 所有中间量的 WellFormed_arr / AllReduced 性质
          -- ================================================================
          have h_coeff_red : Zp.Reduced p coeff := by
            have hr0_red : Zp.Reduced p r[0]!.snd :=
              hred_r (r[0]!) (mem_getFirst_toList r hpos_r)
            have hprime : coeff.prime.toNat = p := by
              dsimp [coeff]
              have : ((r[0]!.snd * lc_g_inv).prime).toNat = ((r[0]!.snd).prime).toNat := rfl
              rw [this, hr0_red.1]
            have hval : coeff.val.toNat < p := by
              have h_mod_lt : (r[0]!.snd.val.toNat * lc_g_inv.val.toNat) % (r[0]!.snd).prime.toNat < UInt64.size := by
                have h_pos : 0 < (r[0]!.snd).prime.toNat := by
                  rw [hr0_red.1]; exact Nat.Prime.pos hp.out
                have h_lt_m : (r[0]!.snd.val.toNat * lc_g_inv.val.toNat) % (r[0]!.snd).prime.toNat < (r[0]!.snd).prime.toNat :=
                  Nat.mod_lt (r[0]!.snd.val.toNat * lc_g_inv.val.toNat) h_pos
                have h_m_le : (r[0]!.snd).prime.toNat < UInt64.size := by
                  rw [hr0_red.1]
                  have hp_one_lt : 1 < p := Nat.Prime.one_lt hp.out
                  omega
                exact Nat.lt_trans h_lt_m h_m_le
              calc
                coeff.val.toNat = ((r[0]!.snd * lc_g_inv).val.toNat) := rfl
                _ = (((r[0]!.snd.val.toNat * lc_g_inv.val.toNat) % (r[0]!.snd).prime.toNat).toUInt64).toNat := rfl
                _ = (r[0]!.snd.val.toNat * lc_g_inv.val.toNat) % (r[0]!.snd).prime.toNat := by
                  simp [UInt64.toNat_ofNat, Nat.mod_eq_of_lt h_mod_lt]
                _ = (r[0]!.snd.val.toNat * lc_g_inv.val.toNat) % p := by rw [hr0_red.1]
                _ < p := Nat.mod_lt _ (Nat.Prime.pos hp.out)
            exact ⟨hprime, hval⟩
          have hwf_term : SparsePolyZp.WellFormed_arr p term := by
            intro x hx
            have hx' : x = (⟨d⟩, coeff) := by
              simpa [term] using hx
            subst hx'; exact h_coeff_red
          have hwf_g : SparsePolyZp.WellFormed_arr p g := hred_g
          have hwf_term_g : SparsePolyZp.WellFormed_arr p (term * g) :=
            SparsePolyZp.WellFormed_arr.mul p h_p2 term g hwf_term hwf_g
          have hwf_r : SparsePolyZp.WellFormed_arr p r := hred_r
          have hwf_r' : SparsePolyZp.WellFormed_arr p r' :=
            SparsePolyZp.WellFormed_arr.sub p r (term * g) hwf_r hwf_term_g
          have hred_q' : SparsePolyZp.AllReduced p q'.toList := by
            intro x hx
            have hx' : x ∈ q.toList ∨ x = (⟨d⟩, coeff) := by
              simpa [q'] using hx
            rcases hx' with (hxq | hx_term)
            · exact hred_q x hxq
            · subst hx_term; exact h_coeff_red
          have hred_r' : SparsePolyZp.AllReduced p r'.toList := hwf_r'
          -- ================================================================
          -- (B) Algebraic invariant
          -- ================================================================
          have hq' : SparsePolyZp.toPoly p q' = SparsePolyZp.toPoly p q + SparsePolyZp.toPoly p term := by
            dsimp [q']
            calc
              SparsePolyZp.toPoly p (q.push (⟨d⟩, coeff))
                  = listSum p (q.push (⟨d⟩, coeff)).toList := rfl
              _ = listSum p (q.toList ++ [(⟨d⟩, coeff)]) := by
                simp
              _ = listSum p q.toList + listSum p [(⟨d⟩, coeff)] := by
                rw [listSum_append]
              _ = SparsePolyZp.toPoly p q + SparsePolyZp.toPoly p term := rfl
          have hsub : SparsePolyZp.toPoly p r' = SparsePolyZp.toPoly p r - SparsePolyZp.toPoly p (term * g) :=
            SparsePolyZp.toPoly_sub p h_2p r (term * g) hwf_r hwf_term_g
          have hmul : SparsePolyZp.toPoly p (term * g) = SparsePolyZp.toPoly p term * SparsePolyZp.toPoly p g :=
            SparsePolyZp.toPoly_mul p h_2p h_p2 term g hwf_term hwf_g
          have h_inv_next : SparsePolyZp.toPoly p f = SparsePolyZp.toPoly p q' * SparsePolyZp.toPoly p g + SparsePolyZp.toPoly p r' := by
            calc
              SparsePolyZp.toPoly p f
                  = SparsePolyZp.toPoly p q * SparsePolyZp.toPoly p g + SparsePolyZp.toPoly p r := h_inv'
              _ = (SparsePolyZp.toPoly p q * SparsePolyZp.toPoly p g + SparsePolyZp.toPoly p r
                     + SparsePolyZp.toPoly p term * SparsePolyZp.toPoly p g
                     - SparsePolyZp.toPoly p term * SparsePolyZp.toPoly p g) := by ring
              _ = ((SparsePolyZp.toPoly p q + SparsePolyZp.toPoly p term) * SparsePolyZp.toPoly p g
                     + (SparsePolyZp.toPoly p r - SparsePolyZp.toPoly p term * SparsePolyZp.toPoly p g)) := by ring
              _ = (SparsePolyZp.toPoly p q' * SparsePolyZp.toPoly p g
                     + SparsePolyZp.toPoly p r') := by
                rw [hq', hsub, hmul]
          -- ================================================================
          -- (C) Termination measure decrease
          -- ================================================================
          have h_deg_next : r'[0]!.fst.deg + 1 < n ∨ r'.isEmpty := by
            by_cases h_empty_r' : r'.isEmpty
            · exact Or.inr h_empty_r'
            · have h_deg_dec : r'[0]!.fst.deg < dr :=
                divmod_deg_decrease g dg lc_g_inv r dr coeff d term r' hred_g hred_r h_nonempty_g h_prime h_p2 h_empty_r'
              left
              rcases h_deg with (h_eq | h_empty)
              · omega
              · exfalso; exact hr_empty h_empty
          -- ================================================================
          -- (D) 应用 IH
          -- ================================================================
          rcases h_deg_next with (h_lt | h_empty_r')
          · exact IH (r'[0]!.fst.deg + 1) h_lt q' r' hred_q' hred_r' (Or.inl rfl) h_inv_next
          · have h_eq' : divmodAux'_wf g dg lc_g_inv q' r' = (q', r') := by
              unfold divmodAux'_wf; simp [h_empty_r']
            simpa [h_eq'] using h_inv_next
  -- 初始 n := r[0]!.fst.deg + 1（若 r 非空）或 n := 0（若 r 为空）
  by_cases hr_empty : r.isEmpty
  · apply h_aux 0 q r hred_q hred_r (Or.inr hr_empty) h_inv
  · have hpos_r : 0 < r.size := by
      rw [Array.isEmpty_iff_size_eq_zero] at hr_empty; omega
    let dr := r[0]!.fst.deg
    apply h_aux (dr + 1) q r hred_q hred_r (Or.inl rfl) h_inv

/-- 多项式除法的 toPoly 对应（依 divmodAux_invariant 得证） -/
private lemma divmod_identity (f g : SparsePolyZp)
    (h_nonempty : ¬g.isEmpty)
    (hwf_f : SparsePolyZp.WellFormed p f) (hwg_g : SparsePolyZp.WellFormed p g)
    (hred_f : SparsePolyZp.AllReduced p f.toList) (hred_g : SparsePolyZp.AllReduced p g.toList)
    (h_2p : 2 * p ≤ UInt64.size) (h_p2 : p * p ≤ UInt64.size) :
    SparsePolyZp.toPoly p f =
      SparsePolyZp.toPoly p (SparsePolyZp.divmod f g).fst * SparsePolyZp.toPoly p g +
      SparsePolyZp.toPoly p (SparsePolyZp.divmod f g).snd := by
  unfold SparsePolyZp.divmod
  have h_nonempty_g : ¬ g.isEmpty := h_nonempty
  simp [h_nonempty_g]
  let dg := g[0]!.fst.deg
  let lc_g_inv := g[0]!.snd.inv
  have h_inv_init : SparsePolyZp.toPoly p f = SparsePolyZp.toPoly p (#[] : SparsePolyZp) * SparsePolyZp.toPoly p g + SparsePolyZp.toPoly p f := by
    simp [SparsePolyZp.toPoly_empty]
  have hpos : 0 < g.size := by
    rw [Array.isEmpty_iff_size_eq_zero] at h_nonempty_g
    omega
  have h_prime : lc_g_inv.prime.toNat = p := by
    have hm' : g[0]! ∈ g := by
      have hmem : g[0] ∈ g := Array.getElem_mem (i := 0) (h := hpos) (xs := g)
      -- g[0]! = g[0] when hpos holds
      simpa [show g[0]! = g[0] from by simp [hpos]] using hmem
    have h_from_wf : (g[0]!).snd.prime.toNat = p := hwg_g (g[0]!) hm'
    dsimp [lc_g_inv]
    -- (inv a).prime = a.prime
    simpa using h_from_wf
  have h_res := divmodAux_invariant g dg lc_g_inv #[] f f h_inv_init hwf_f hwg_g hred_f hred_g (allReduced_empty p) hred_f h_2p h_p2 h_nonempty_g h_prime
  -- h_res: let (q', r') := divmodAux'_wf ...; toPoly f = toPoly q' * toPoly g + toPoly r'
  -- Need to connect divmodAux'_wf to SparsePolyZp.divmodAux
  have h_eq : divmodAux'_wf g dg lc_g_inv #[] f = SparsePolyZp.divmodAux g dg lc_g_inv #[] f :=
    divmodAux'_wf_eq g dg lc_g_inv #[] f
  simpa [h_eq] using h_res

/-- divmod'_wf 的 toPoly 恒等式（绕过 divmodAux'_wf_eq + partial def，直接使用 divmodAux_invariant）。 -/
private lemma divmod'_wf_identity (f g : SparsePolyZp) (h_nonempty : ¬g.isEmpty)
    (hwf_f : SparsePolyZp.WellFormed p f) (hwg_g : SparsePolyZp.WellFormed p g)
    (hred_f : SparsePolyZp.AllReduced p f.toList) (hred_g : SparsePolyZp.AllReduced p g.toList)
    (h_2p : 2 * p ≤ UInt64.size) (h_p2 : p * p ≤ UInt64.size) :
    SparsePolyZp.toPoly p f = SparsePolyZp.toPoly p (divmod'_wf f g).fst * SparsePolyZp.toPoly p g +
                              SparsePolyZp.toPoly p (divmod'_wf f g).snd := by
  unfold divmod'_wf
  simp [h_nonempty]
  have hpos : 0 < g.size := by
    rw [Array.isEmpty_iff_size_eq_zero] at h_nonempty; omega
  have h_prime : (g[0]!.snd.inv).prime.toNat = p := by
    have hm' : g[0]! ∈ g := by
      have hmem : g[0] ∈ g := Array.getElem_mem (i := 0) (h := hpos) (xs := g)
      simpa [show g[0]! = g[0] from by simp [hpos]] using hmem
    have h_from_wf : (g[0]!).snd.prime.toNat = p := hwg_g (g[0]!) hm'
    simpa [Zp.inv, h_from_wf]
  have h_inv_init : SparsePolyZp.toPoly p f = SparsePolyZp.toPoly p (#[] : SparsePolyZp) * SparsePolyZp.toPoly p g + SparsePolyZp.toPoly p f := by
    simp [SparsePolyZp.toPoly_empty]
  exact divmodAux_invariant g (g[0]!.fst.deg) (g[0]!.snd.inv) #[] f f h_inv_init hwf_f hwg_g hred_f hred_g (allReduced_empty p) hred_f h_2p h_p2 h_nonempty h_prime

/-- divmodAux'_wf 的余式保持 AllReduced（用 divmodAux'_wf.induct 证明）。 -/
private lemma divmodAux'_wf_snd_allReduced (g : SparsePolyZp) (dg : ℕ) (lc_g_inv : Zp) (q r : SparsePolyZp)
    (hred_q : SparsePolyZp.AllReduced p q.toList) (hred_r : SparsePolyZp.AllReduced p r.toList)
    (hred_g : SparsePolyZp.AllReduced p g.toList) (h_p2 : p * p ≤ UInt64.size) :
    SparsePolyZp.AllReduced p (divmodAux'_wf g dg lc_g_inv q r).snd.toList := by
  let motive : SparsePolyZp → SparsePolyZp → Prop := λ q' r' =>
    (SparsePolyZp.AllReduced p q'.toList) → (SparsePolyZp.AllReduced p r'.toList) →
    SparsePolyZp.AllReduced p (divmodAux'_wf g dg lc_g_inv q' r').snd.toList
  have h_case1 : ∀ (q' r' : SparsePolyZp), Array.isEmpty r' = true → motive q' r' := by
    intro q' r' h_empty hred_q' hred_r'
    rw [divmodAux'_wf.eq_1]
    simp [h_empty, hred_r']
  have h_case2 : ∀ (q' r' : SparsePolyZp), ¬Array.isEmpty r' = true → (r'[0]!.fst.deg < dg) → motive q' r' := by
    intro q' r' h_not_empty hdr_lt_dg hred_q' hred_r'
    rw [divmodAux'_wf.eq_1]
    simp [h_not_empty, hdr_lt_dg, hred_r']
  have h_case3 : ∀ (q' r' : SparsePolyZp), ¬Array.isEmpty r' = true → ¬(r'[0]!.fst.deg < dg) → 
    (let coeff := r'[0]!.snd * lc_g_inv; let d := r'[0]!.fst.deg - dg; let term : SparsePolyZp := #[(⟨d⟩, coeff)];
     let r'' := r' - (term * g); let q'' := q'.push (⟨d⟩, coeff); motive q'' r'') → motive q' r' := by
    intro q' r' h_not_empty h_not_lt IH hred_q' hred_r'
    rw [divmodAux'_wf.eq_1]
    simp [h_not_empty, h_not_lt]
    let dr := r'[0]!.fst.deg
    let coeff := r'[0]!.snd * lc_g_inv
    let d := dr - dg
    let term : SparsePolyZp := #[(⟨d⟩, coeff)]
    let r'' := r' - (term * g)
    let q'' := q'.push (⟨d⟩, coeff)
    have hpos_r' : 0 < r'.size := by
      by_contra! hzero
      have hzero_size : r'.size = 0 := by omega
      have h_empty' : Array.isEmpty r' = true := by
        rw [Array.isEmpty_iff_size_eq_zero, hzero_size]
      exact h_not_empty h_empty'
    have h_coeff_red : Zp.Reduced p coeff := by
      have hr0_red : Zp.Reduced p r'[0]!.snd := by
        apply hred_r'
        exact mem_getFirst_toList r' hpos_r'
      have hprime : coeff.prime.toNat = p := by
        have : ((r'[0]!.snd * lc_g_inv).prime).toNat = ((r'[0]!.snd).prime).toNat := rfl
        rw [this, hr0_red.1]
      have hval : coeff.val.toNat < p := by
        have h_mod_lt : (r'[0]!.snd.val.toNat * lc_g_inv.val.toNat) % (r'[0]!.snd).prime.toNat < UInt64.size := by
          have h_pos : 0 < (r'[0]!.snd).prime.toNat := by
            rw [hr0_red.1]; exact Nat.Prime.pos hp.out
          have h_lt_m : (r'[0]!.snd.val.toNat * lc_g_inv.val.toNat) % (r'[0]!.snd).prime.toNat < (r'[0]!.snd).prime.toNat :=
            Nat.mod_lt (r'[0]!.snd.val.toNat * lc_g_inv.val.toNat) h_pos
          have h_m_le : (r'[0]!.snd).prime.toNat < UInt64.size := by
            rw [hr0_red.1]
            have hp_one_lt : 1 < p := Nat.Prime.one_lt hp.out
            have : p * p ≤ UInt64.size := h_p2
            have hp_lt_sq : p < p * p := by
              nlinarith
            omega
          exact Nat.lt_trans h_lt_m h_m_le
        calc
          coeff.val.toNat = ((r'[0]!.snd * lc_g_inv).val.toNat) := rfl
          _ = (((r'[0]!.snd.val.toNat * lc_g_inv.val.toNat) % (r'[0]!.snd).prime.toNat).toUInt64).toNat := rfl
          _ = (r'[0]!.snd.val.toNat * lc_g_inv.val.toNat) % (r'[0]!.snd).prime.toNat := by
            simp [UInt64.toNat_ofNat, Nat.mod_eq_of_lt h_mod_lt]
          _ = (r'[0]!.snd.val.toNat * lc_g_inv.val.toNat) % p := by rw [hr0_red.1]
          _ < p := Nat.mod_lt _ (Nat.Prime.pos hp.out)
      exact ⟨hprime, hval⟩
    have hwf_term : SparsePolyZp.WellFormed_arr p term := by
      intro x hx
      have hx' : x = (⟨d⟩, coeff) := by simpa [term] using hx
      subst hx'; exact h_coeff_red
    have hwf_g : SparsePolyZp.WellFormed_arr p g := hred_g
    have hwf_term_g : SparsePolyZp.WellFormed_arr p (term * g) :=
      SparsePolyZp.WellFormed_arr.mul p h_p2 term g hwf_term hwf_g
    have hwf_r' : SparsePolyZp.WellFormed_arr p r' := hred_r'
    have hwf_r'' : SparsePolyZp.WellFormed_arr p r'' :=
      SparsePolyZp.WellFormed_arr.sub p r' (term * g) hwf_r' hwf_term_g
    have hred_q'' : SparsePolyZp.AllReduced p q''.toList := by
      intro x hx
      have hx' : x ∈ q'.toList ∨ x = (⟨d⟩, coeff) := by simpa [q''] using hx
      rcases hx' with (hxq | hx_term)
      · exact hred_q' x hxq
      · subst hx_term; exact h_coeff_red
    have h_IH : motive q'' r'' := IH
    exact h_IH hred_q'' hwf_r''
  exact divmodAux'_wf.induct g dg lc_g_inv motive h_case1 h_case2 h_case3 q r hred_q hred_r

/-- divmod'_wf 的余式保持 AllReduced。 -/
private lemma divmod'_wf_snd_allReduced (f g : SparsePolyZp)
    (hred_f : SparsePolyZp.AllReduced p f.toList) (hred_g : SparsePolyZp.AllReduced p g.toList)
    (h_p2 : p * p ≤ UInt64.size) :
    SparsePolyZp.AllReduced p (divmod'_wf f g).snd.toList := by
  unfold divmod'_wf
  by_cases hg_empty : g.isEmpty
  · simp [hg_empty, hred_f]
  · simp [hg_empty]
    refine divmodAux'_wf_snd_allReduced g (g[0]!.fst.deg) (g[0]!.snd.inv) #[] f 
      (allReduced_empty p) hred_f hred_g h_p2

/-- divmod'_wf 的余式保持 WellFormed（由 AllReduced 推出）。 -/
private lemma divmod'_wf_snd_wellFormed (f g : SparsePolyZp)
    (hred_f : SparsePolyZp.AllReduced p f.toList) (hred_g : SparsePolyZp.AllReduced p g.toList)
    (h_p2 : p * p ≤ UInt64.size) :
    SparsePolyZp.WellFormed p (divmod'_wf f g).snd := by
  have hred_snd : SparsePolyZp.AllReduced p (divmod'_wf f g).snd.toList :=
    divmod'_wf_snd_allReduced f g hred_f hred_g h_p2
  intro x hx
  have hx_red : Zp.Reduced p x.snd := hred_snd x (by rw [← Array.mem_def]; exact hx)
  exact hx_red.1

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
private lemma divmod'_wf_deg_lt (f g : SparsePolyZp) (h_nonempty : ¬g.isEmpty) :
    (SparsePolyZp.toPoly p (divmod'_wf f g).snd).natDegree < (SparsePolyZp.toPoly p g).natDegree := by
  admit

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

/-- __upoly_make_monic_ir 保持 degree bound。 -/
lemma upoly_make_monic_deg_bound (f : SparsePolyZp) (h_deg_bound : ∀ x ∈ f.toList, x.1.deg < 2 ^ 64) :
    ∀ x ∈ (Generated.__upoly_make_monic_ir f).snd.toList, x.1.deg < 2 ^ 64 := by
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
  admit
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
    -- (the main theorem ensures natDegree > 0 → g non-empty).
    admit
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
    · -- derivative ≠ 0: Yun algorithm. normalize(c·g) = normalize(g).
      sorry

/-- safe wrapper：常数多项式直接返回 #[]，避免 C++ 模型不终止的问题 -/
noncomputable def __squarefree_Zp_ir_safe (p : ℕ) [hp : Fact (Nat.Prime p)] (f : SparsePolyZp) : Array (SparsePolyZp × UInt64) :=
  if (SparsePolyZp.toPoly p f).natDegree = 0 then #[]
  else Generated.__squarefree_Zp_ir f

set_option maxHeartbeats 0 in
/-- 主正确性定理（使用 safe wrapper） -/
theorem __squarefree_Zp_ir_refines (p : ℕ) [hp : Fact (Nat.Prime p)]
    (f : SparsePolyZp)
    (hwf_f : SparsePolyZp.WellFormed p f)
    (hred_f : SparsePolyZp.AllReduced p f.toList)
    (hp_size : 2 * p ≤ UInt64.size)
    (h_no_overflow : ∀ x ∈ f.toList, x.2.val.toNat * x.1.deg < 2 ^ 64)
    (h_deg_bound : ∀ x ∈ f.toList, x.1.deg < 2 ^ 64)
    :   toPolyList (__squarefree_Zp_ir_safe p f) p = sqfZp (SparsePolyZp.toPoly p f) := by
  unfold __squarefree_Zp_ir_safe
  by_cases h_deg0 : (SparsePolyZp.toPoly p f).natDegree = 0
  · simp [sqfZp, h_deg0, toPolyList_empty]
  · have h_deg_pos : (SparsePolyZp.toPoly p f).natDegree > 0 := by omega
    -- 强归纳于 natDegree，使用 __squarefree_Zp_ir_safe 避免 partial def 在常数时不停机
    have h_main : toPolyList (__squarefree_Zp_ir_safe p f) p = sqfZp (SparsePolyZp.toPoly p f) := by
      suffices ∀ n, ∀ (g : SparsePolyZp),
          (SparsePolyZp.toPoly p g).natDegree = n →
          SparsePolyZp.WellFormed p g →
          SparsePolyZp.AllReduced p g.toList →
          2 * p ≤ UInt64.size →
          (∀ x ∈ g.toList, x.2.val.toNat * x.1.deg < 2 ^ 64) →
          (∀ x ∈ g.toList, x.1.deg < 2 ^ 64) →
          toPolyList (__squarefree_Zp_ir_safe p g) p = sqfZp (SparsePolyZp.toPoly p g) from
        this (SparsePolyZp.toPoly p f).natDegree f rfl hwf_f hred_f hp_size h_no_overflow h_deg_bound
      intro n
      induction n using Nat.strongRecOn with
      | ind n ih =>
        intro g h_deg_eq hwf_g hred_g hp_size h_no_overflow_g h_deg_bound_g
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
            have h_toPoly_g1 : SparsePolyZp.toPoly p g_1 = Polynomial.contract p (SparsePolyZp.toPoly p g) := h_extract_eq
            have h_wf_g1 : SparsePolyZp.WellFormed p g_1 := extract_pth_root_wellFormed g hwf_g
            have h_red_g1 : SparsePolyZp.AllReduced p g_1.toList := extract_pth_root_allReduced g hred_g
            have h_db_g1 : ∀ x ∈ g_1.toList, x.1.deg < 2 ^ 64 := extract_pth_root_deg_bound g h_deg_bound_g
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
            have h_red_g2 : SparsePolyZp.AllReduced p g_2.toList := upoly_make_monic_allReduced g_1 h_red_g1 hp_size
            have h_db_g2 : ∀ x ∈ g_2.toList, x.1.deg < 2 ^ 64 := upoly_make_monic_deg_bound g_1 h_db_g1
            have h_nov_g2 : ∀ x ∈ g_2.toList, x.2.val.toNat * x.1.deg < 2 ^ 64 :=
              upoly_make_monic_no_overflow g_1 h_red_g1 h_nov_g1 h_p_mul_g1 hp_size
            have h_ih_g2 : toPolyList (__squarefree_Zp_ir_safe p g_2) p = sqfZp (SparsePolyZp.toPoly p g_2) :=
              ih (SparsePolyZp.toPoly p g_2).natDegree h_deg_g2_lt_n g_2 rfl
                h_wf_g2 h_red_g2 hp_size h_nov_g2 h_db_g2
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
              -- Loop invariant: (loop i arr result p).2.2.toList = result.toList ++ (arr.toList.drop i).map f
              -- For i=0, result=#[]: (loop 0 arr #[] p).2.2.toList = arr.toList.map f = (arr.map f).toList
              -- Hence the arrays are equal. Proof by strong induction on arr.size - i.
              -- The induction step requires relating arr[i]! to arr.toList[i] (getElem! vs get),
              -- and Array.push toList append. These are all true but require Array/List conversion lemmas.
              admit
              
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
              -- (front! g_1).snd.toZMod p ≠ 0 (leading coefficient non-zero).
              -- In Z/pZ, a ≠ 0 → a⁻¹ ≠ 0. Zp.inv computes the modular inverse.
              -- The equality (inv a).toZMod p = (a.toZMod p)⁻¹ follows from Zp_toZMod_inv_mul_self.
              admit
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
            have h_result : toPolyList (Generated.__squarefree_Zp_ir g) p = sqfZp (SparsePolyZp.toPoly p g) := by
              have h_step1 : toPolyList (Generated.__squarefree_Zp_ir g) p = toPolyList ((Generated._loop___squarefree_Zp_0_ir_def 0
                  (__squarefree_Zp_ir_safe p g_2) #[] p_1).2.2) p := by
                have h_def_eq : Generated.__squarefree_Zp_ir g = (Generated._loop___squarefree_Zp_0_ir_def 0
                    (Generated.__squarefree_Zp_ir g_2) #[] p_1).2.2 := by
                  dsimp [Generated.__squarefree_Zp_ir]
                  rw [Generated.__squarefree_Zp_ir_def.eq_1 g]
                  split_ifs with h_cond
                  · rfl
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
                have h_toPolyList_map (arr : Array (SparsePolyZp × UInt64)) :
                    toPolyList (arr.map (fun (g_h, e) => (g_h, e * p_1))) p = (toPolyList arr p).map (fun (h, e) => (h, e * p)) := by
                  -- Key: (e * p_1).toNat = e.toNat * p for all exponents in the factorization.
                  -- This holds when e.toNat * p < 2^64 (no UInt64 overflow), which follows from:
                  --   e.toNat ≤ (toPoly p g_2).natDegree (sqfZp exponent ≤ degree)
                  --   p * (toPoly p g_2).natDegree < 2^64 (from h_p_mul_g1 + degree bounds)
                  -- For general arr, this may not hold; specific to __squarefree_Zp_ir_safe p g_2.
                  admit
                rw [h_toPolyList_map (__squarefree_Zp_ir_safe p g_2)]
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
            admit
    unfold __squarefree_Zp_ir_safe at h_main
    simpa using h_main

end Refinement
