/-
  CLPoly/Math/FiniteFieldFact.lean — L3 有限域数学定理

  Phase 2 数学基石：DDF/EDF 正确性所需的纯数学定理。

  T2.1: X^{p^d} - X 的 Separable 性              ✓ 无 sorry
  T2.2: 不可约多项式整除 X^{p^d} - X             ✓ 无 sorry
  T2.2': T2.2 的逆命题                            ✓ 无 sorry
  T2.3: gcd(X^{p^d}-X, f) 的不可约因子刻画       ✓ 无 sorry
  T2.4: EDF 三分性（Euler 判据）                  ✓ 无 sorry
  T2.5: 半幂根计数                                ✓ 无 sorry
-/
import Mathlib.FieldTheory.Separable
import Mathlib.FieldTheory.Finite.Basic
import Mathlib.RingTheory.AdjoinRoot
import Mathlib.Data.ZMod.Basic
import Mathlib.Algebra.CharP.Lemmas

set_option autoImplicit false
set_option maxHeartbeats 400000

open Polynomial

variable {p : ℕ} [Fact (Nat.Prime p)]

-- ============================================================
-- T2.1: X^{p^d} - X 的 Separable 性
-- ============================================================

/-- X^(p^d) - X 在 ZMod p 上是 separable 的。
    证明：导数 = -1（char p 下），-1 是单位。 -/
theorem X_pow_sub_X_separable (d : ℕ) (hd : 0 < d) :
    Separable (X ^ (p ^ d) - X : Polynomial (ZMod p)) := by
  rw [separable_def, derivative_sub, derivative_X_pow, derivative_X]
  have hp : (↑(p ^ d) : ZMod p) = 0 := by
    rw [Nat.cast_pow, CharP.cast_eq_zero, zero_pow (by omega : d ≠ 0)]
  rw [hp, map_zero, zero_mul, zero_sub]
  exact ⟨0, -1, by ring⟩

-- ============================================================
-- T2.2: 不可约多项式整除 X^{p^d} - X
-- ============================================================

/-- 不可约多项式 g 整除 X^{p^d} - X，当 deg(g) | d。
    证明：在 AdjoinRoot g 中，root g 满足 Frobenius 迭代。 -/
theorem irreducible_dvd_X_pow_sub_X
    (g : Polynomial (ZMod p)) (hg : Irreducible g)
    (d : ℕ) (hdvd : g.natDegree ∣ d) :
    g ∣ (X ^ (p ^ d) - X : Polynomial (ZMod p)) := by
  haveI : Fact (Irreducible g) := ⟨hg⟩
  letI : Fintype (AdjoinRoot g) :=
    Fintype.ofEquiv _ (AdjoinRoot.powerBasis hg.ne_zero).basis.equivFun.toEquiv.symm
  rw [← AdjoinRoot.mk_eq_zero, ← AdjoinRoot.aeval_eq]
  simp only [map_sub, map_pow, aeval_X]
  rw [sub_eq_zero]
  have hcard : Fintype.card (AdjoinRoot g) = p ^ g.natDegree := by
    rw [Fintype.card_congr (AdjoinRoot.powerBasis hg.ne_zero).basis.equivFun.toEquiv,
        Fintype.card_fun, Fintype.card_fin, ZMod.card]
    rfl
  have hexp : p ^ d = Fintype.card (AdjoinRoot g) ^ (d / g.natDegree) := by
    rw [hcard, ← pow_mul, Nat.mul_div_cancel' hdvd]
  rw [hexp]
  exact FiniteField.pow_card_pow _ _

-- ============================================================
-- 辅助引理：AdjoinRoot 中所有元素满足 Frobenius 不动点
-- ============================================================

/-- 若 g 不可约且 g | X^{p^d} - X，则 AdjoinRoot g 中每个元素满足 x^{p^d} = x。
    提取自 T2.2' 的 Step A+B，供 T2.4 复用。 -/
lemma all_pow_eq_of_dvd_X_pow_sub_X
    (g : Polynomial (ZMod p)) (hg : Irreducible g)
    (d : ℕ) (hdvd : g ∣ (X ^ (p ^ d) - X : Polynomial (ZMod p))) :
    ∀ (x : AdjoinRoot g), x ^ (p ^ d) = x := by
  haveI : Fact (Irreducible g) := ⟨hg⟩
  haveI : CharP (AdjoinRoot g) p :=
    charP_of_injective_algebraMap (algebraMap (ZMod p) (AdjoinRoot g)).injective p
  letI : Fintype (AdjoinRoot g) :=
    Fintype.ofEquiv _ (AdjoinRoot.powerBasis hg.ne_zero).basis.equivFun.toEquiv.symm
  have hα : (AdjoinRoot.root g) ^ (p ^ d) = AdjoinRoot.root g := by
    have h0 := AdjoinRoot.mk_eq_zero.mpr hdvd
    simp only [map_sub, map_pow, AdjoinRoot.mk_X] at h0
    rwa [sub_eq_zero] at h0
  let φ : AdjoinRoot g →ₐ[ZMod p] AdjoinRoot g :=
    { iterateFrobenius (AdjoinRoot g) p d with
      commutes' := fun r => by
        change (algebraMap (ZMod p) (AdjoinRoot g) r) ^ (p ^ d) =
          algebraMap (ZMod p) (AdjoinRoot g) r
        rw [← map_pow]; congr 1
        have := FiniteField.pow_card_pow d r
        rwa [ZMod.card p] at this }
  have hφ_apply : ∀ x, φ x = x ^ (p ^ d) := fun _ => rfl
  have hφ_eq_id : φ = AlgHom.id (ZMod p) (AdjoinRoot g) := by
    apply AdjoinRoot.algHom_ext
    simp [hφ_apply, hα]
  intro x
  have := congr_fun (congr_arg DFunLike.coe hφ_eq_id) x
  simp [hφ_apply] at this
  exact this

-- ============================================================
-- T2.2': 逆命题 — q | X^{p^d} - X → deg q | d
-- ============================================================

/-- T2.2 的逆命题：不可约 g 整除 X^{p^d} - X → deg(g) | d。
    证明：由辅助引理得 ∀ x, x^{p^d} = x，然后反证法 + 根计数矛盾。 -/
theorem dvd_X_pow_sub_X_imp_natDegree_dvd
    (g : Polynomial (ZMod p)) (hg : Irreducible g)
    (d : ℕ) (h : g ∣ (X ^ (p ^ d) - X : Polynomial (ZMod p))) :
    g.natDegree ∣ d := by
  haveI : Fact (Irreducible g) := ⟨hg⟩
  letI : Fintype (AdjoinRoot g) :=
    Fintype.ofEquiv _ (AdjoinRoot.powerBasis hg.ne_zero).basis.equivFun.toEquiv.symm
  set k := g.natDegree
  have hk_pos : 0 < k := Irreducible.natDegree_pos hg
  have hp : Nat.Prime p := Fact.out
  have hcard : Fintype.card (AdjoinRoot g) = p ^ k := by
    rw [Fintype.card_congr (AdjoinRoot.powerBasis hg.ne_zero).basis.equivFun.toEquiv,
        Fintype.card_fun, Fintype.card_fin, ZMod.card]; rfl
  have hall := all_pow_eq_of_dvd_X_pow_sub_X g hg d h
  -- Step C: k | d by contradiction (root counting)
  by_contra hnd
  set r := d % k
  have hr_pos : 0 < r := Nat.pos_of_ne_zero (fun h => hnd (Nat.dvd_of_mod_eq_zero h))
  have hr_lt : r < k := Nat.mod_lt d hk_pos
  -- ∀ x, x^{p^r} = x
  have hall_r : ∀ (x : AdjoinRoot g), x ^ (p ^ r) = x := by
    intro x
    have h2 : x ^ ((p ^ k) ^ (d / k)) = x := by
      rw [← hcard]; exact FiniteField.pow_card_pow (d / k) x
    have : x ^ (p ^ d) = x ^ (p ^ r) := by
      rw [show p ^ d = (p ^ k) ^ (d / k) * p ^ r from by
        rw [← pow_mul, ← pow_add]; congr 1; exact (Nat.div_add_mod d k).symm]
      rw [pow_mul, h2]
    rw [← this]; exact hall x
  -- X^{p^r} - X has degree p^r but p^k > p^r roots → contradiction
  set f := (X ^ (p ^ r) - X : Polynomial (AdjoinRoot g))
  have hf_deg : f.natDegree = p ^ r := by
    have : natDegree (X : Polynomial (AdjoinRoot g)) <
        natDegree ((X : Polynomial (AdjoinRoot g)) ^ (p ^ r)) := by
      rw [natDegree_X, natDegree_X_pow]; exact Nat.one_lt_pow hr_pos.ne' hp.one_lt
    simp only [f, natDegree_sub_eq_left_of_natDegree_lt this, natDegree_X_pow]
  have hf_ne : f ≠ 0 := by
    intro h; rw [h, natDegree_zero] at hf_deg
    exact absurd hf_deg.symm (pow_pos hp.pos r).ne'
  have hroots : ∀ (x : AdjoinRoot g), x ∈ f.roots := by
    intro x; rw [mem_roots hf_ne, IsRoot, eval_sub, eval_pow, eval_X, sub_eq_zero]; exact hall_r x
  exact absurd (Nat.pow_lt_pow_right hp.one_lt hr_lt) (not_lt.mpr (
    calc p ^ k = Fintype.card (AdjoinRoot g) := hcard.symm
      _ = Finset.univ.card := Finset.card_univ.symm
      _ ≤ f.roots.toFinset.card :=
          Finset.card_le_card (fun x _ => Multiset.mem_toFinset.mpr (hroots x))
      _ ≤ f.roots.card := Multiset.toFinset_card_le f.roots
      _ ≤ f.natDegree := card_roots' f
      _ = p ^ r := hf_deg))

-- ============================================================
-- T2.3: gcd(X^{p^d} - X, f) 的不可约因子刻画
-- ============================================================

/-- gcd(X^{p^d} - X, f) 的不可约因子恰为 f 中度整除 d 的不可约因子。
    正方向用 T2.2'，反方向用 T2.2。 -/
theorem gcd_X_pow_sub_X_factors
    (f : Polynomial (ZMod p)) (d : ℕ)
    (q : Polynomial (ZMod p)) (hq : Irreducible q) :
    q ∣ EuclideanDomain.gcd (X ^ (p ^ d) - X) f ↔ q ∣ f ∧ q.natDegree ∣ d := by
  constructor
  · -- 正方向: q | gcd → q | f ∧ deg q | d
    intro h
    exact ⟨dvd_trans h (EuclideanDomain.gcd_dvd_right _ _),
           dvd_X_pow_sub_X_imp_natDegree_dvd q hq d
             (dvd_trans h (EuclideanDomain.gcd_dvd_left _ _))⟩
  · -- 反方向: q | f ∧ deg q | d → q | gcd
    intro ⟨hqf, hqd⟩
    exact EuclideanDomain.dvd_gcd
      (irreducible_dvd_X_pow_sub_X q hq d hqd) hqf

-- ============================================================
-- T2.4: EDF 三分性 — Euler 判据
-- ============================================================

/-- EDF 三分性：对 p 奇素数，不可约 q 满足 deg(q) | d 时，
    任意多项式 a 满足 q | a ∨ q | (a^m - 1) ∨ q | (a^m + 1)，
    其中 m = (p^d - 1) / 2。Cantor-Zassenhaus EDF 算法的核心数学依据。

    证明：在 AdjoinRoot q 中 ā^{p^d} = ā，直接因式分解
    ā^{p^d} - ā = ā · (ā^m - 1) · (ā^m + 1) = 0（利用 p 奇 → 2m+1 = p^d），
    域无零因子得三分，再由 mk_eq_zero 翻译回多项式整除。 -/
theorem edf_trichotomy
    (hp_odd : p ≠ 2)
    (q : Polynomial (ZMod p)) (hq : Irreducible q)
    (d : ℕ) (hdvd : q.natDegree ∣ d)
    (a : Polynomial (ZMod p)) :
    let m := (p ^ d - 1) / 2
    q ∣ a ∨ q ∣ (a ^ m - 1) ∨ q ∣ (a ^ m + 1) := by
  intro m
  haveI : Fact (Irreducible q) := ⟨hq⟩
  have hall := all_pow_eq_of_dvd_X_pow_sub_X q hq d
    (irreducible_dvd_X_pow_sub_X q hq d hdvd)
  set ā := AdjoinRoot.mk q a
  -- p^d = 2m + 1 (p 奇素数 → p^d 奇数)
  have hpd_eq : p ^ d = 2 * m + 1 := by
    have hp_prime : Nat.Prime p := Fact.out
    suffices Odd (p ^ d) by obtain ⟨k, hk⟩ := this; omega
    have : p % 2 = 1 := by
      rcases Nat.mod_two_eq_zero_or_one p with h | h
      · exact absurd (hp_prime.eq_one_or_self_of_dvd 2 (Nat.dvd_of_mod_eq_zero h))
          (by omega)
      · exact h
    exact Nat.odd_iff.mpr (by rw [Nat.pow_mod]; simp [this])
  -- ā · (ā^m - 1) · (ā^m + 1) = 0（ring 恒等式 + Frobenius 不动点）
  have h_triple : ā * (ā ^ m - 1) * (ā ^ m + 1) = 0 := by
    have hring : ā * (ā ^ m - 1) * (ā ^ m + 1) = ā ^ (2 * m + 1) - ā := by ring
    rw [hring, ← hpd_eq, sub_eq_zero]
    exact hall ā
  -- 域无零因子 → 三分，翻译回多项式整除
  rcases mul_eq_zero.mp h_triple with hab | hc
  · rcases mul_eq_zero.mp hab with ha | hb
    · left; rwa [← AdjoinRoot.mk_eq_zero]
    · right; left; rw [← AdjoinRoot.mk_eq_zero]
      have : AdjoinRoot.mk q (a ^ m - 1) = ā ^ m - 1 := by
        simp [ā, map_sub, map_pow, map_one]
      rw [this]; exact hb
  · right; right; rw [← AdjoinRoot.mk_eq_zero]
    have : AdjoinRoot.mk q (a ^ m + 1) = ā ^ m + 1 := by
      simp [ā, map_add, map_pow, map_one]
    rw [this]; exact hc

-- ============================================================
-- T2.5: 半幂根计数
-- ============================================================

/-- 有限域 F（奇特征）中，满足 a^{(|F|-1)/2} = 1 的元素恰有 (|F|-1)/2 个。
    Cantor-Zassenhaus EDF 概率分析的核心代数事实。

    证明：X^{q-1}-1 = (X^m - 1)(X^m + 1) 互素分解（m = (q-1)/2），
    Fermat 给 F* 的 q-1 个根分入两个因子，度数上界各 m，夹逼得各恰 m 个。 -/
theorem card_pow_half_eq_one
    {F : Type*} [Field F] [Fintype F] [DecidableEq F]
    (hF_odd : ringChar F ≠ 2) :
    (Finset.univ.filter (fun a : F => a ^ ((Fintype.card F - 1) / 2) = 1)).card
    = (Fintype.card F - 1) / 2 := by
  set q := Fintype.card F with hq_def
  set m := (q - 1) / 2 with hm_def
  have hp : Nat.Prime (ringChar F) := CharP.char_is_prime F (ringChar F)
  have hq_pos : 0 < q := Fintype.card_pos
  -- q 是奇数（char ≠ 2 → |F| = p^n, p 奇 → q 奇）
  have hq_odd : q % 2 = 1 := by
    have := FiniteField.card F (ringChar F)
    obtain ⟨n, hn_pos, hcard⟩ := this
    rw [hq_def, hcard]
    have : ringChar F % 2 = 1 := by
      rcases Nat.mod_two_eq_zero_or_one (ringChar F) with h | h
      · exact absurd (hp.eq_one_or_self_of_dvd 2 (Nat.dvd_of_mod_eq_zero h)) (by omega)
      · exact h
    rw [Nat.pow_mod]; simp [this]
  -- q ≥ 3（char ≠ 2 → ringChar ≥ 3 → q = ringChar^n ≥ 3）
  have hq_ge3 : 3 ≤ q := by
    obtain ⟨npn, _, hcard⟩ := FiniteField.card F (ringChar F)
    rw [hq_def, hcard]
    calc 3 ≤ ringChar F := by have := hp.two_le; omega
      _ = ringChar F ^ 1 := (pow_one _).symm
      _ ≤ ringChar F ^ (npn : ℕ) := Nat.pow_le_pow_right hp.pos npn.pos
  have h2m : 2 * m = q - 1 := by omega
  have hm_pos : 0 < m := by omega
  -- 1 ≠ -1（char ≠ 2）
  have hone_ne_neg : (1 : F) ≠ -1 := by
    intro h
    have h2 : (1 : F) + 1 = 0 := by
      nth_rewrite 1 [h]; exact neg_add_cancel 1
    haveI : CharP F (ringChar F) := ⟨ringChar.spec F⟩
    have h_dvd : ringChar F ∣ 2 :=
      (CharP.cast_eq_zero_iff F (ringChar F) 2).mp (by exact_mod_cast h2)
    exact absurd (Nat.le_antisymm (Nat.le_of_dvd (by omega) h_dvd) hp.two_le) hF_odd
  -- 度数辅助引理
  have hdeg_lt : natDegree (1 : Polynomial F) < natDegree ((X : Polynomial F) ^ m) := by
    rw [natDegree_one, natDegree_X_pow]; exact hm_pos
  -- A = {a^m = 1}, B = {a^m = -1}
  set A := Finset.univ.filter (fun a : F => a ^ m = 1)
  set B := Finset.univ.filter (fun a : F => a ^ m = -1)
  -- |A| ≤ m（X^m - 1 至多 m 个根）
  have hA_le : A.card ≤ m := by
    have hf_ne : (X ^ m - 1 : Polynomial F) ≠ 0 := by
      intro h
      have := natDegree_sub_eq_left_of_natDegree_lt hdeg_lt
      rw [h, natDegree_zero, natDegree_X_pow] at this; omega
    have hf_deg : (X ^ m - 1 : Polynomial F).natDegree = m := by
      rw [natDegree_sub_eq_left_of_natDegree_lt hdeg_lt, natDegree_X_pow]
    calc A.card
        ≤ (X ^ m - 1 : Polynomial F).roots.toFinset.card := by
          apply Finset.card_le_card; intro a ha
          have ha_eq := (Finset.mem_filter.mp ha).2
          rw [Multiset.mem_toFinset, mem_roots hf_ne, IsRoot,
              eval_sub, eval_pow, eval_X, eval_one, sub_eq_zero]
          exact ha_eq
      _ ≤ (X ^ m - 1 : Polynomial F).roots.card := Multiset.toFinset_card_le _
      _ ≤ (X ^ m - 1 : Polynomial F).natDegree := card_roots' _
      _ = m := hf_deg
  -- |B| ≤ m（X^m + 1 至多 m 个根）
  have hB_le : B.card ≤ m := by
    have hdeg_lt' : natDegree (C (-1 : F)) <
        natDegree ((X : Polynomial F) ^ m) := by
      rw [natDegree_C, natDegree_X_pow]; exact hm_pos
    have hg_ne : (X ^ m + 1 : Polynomial F) ≠ 0 := by
      intro h
      have hd : natDegree (X ^ m - C (-1 : F)) = m := by
        rw [natDegree_sub_eq_left_of_natDegree_lt hdeg_lt', natDegree_X_pow]
      have heq : (X ^ m - C (-1 : F) : Polynomial F) = X ^ m + 1 := by
        rw [C_neg, C_1]; ring
      rw [heq, h, natDegree_zero] at hd; omega
    have hg_deg : (X ^ m + 1 : Polynomial F).natDegree = m := by
      have : (X ^ m + 1 : Polynomial F) = X ^ m - C (-1) := by rw [C_neg, C_1]; ring
      rw [this, natDegree_sub_eq_left_of_natDegree_lt hdeg_lt', natDegree_X_pow]
    calc B.card
        ≤ (X ^ m + 1 : Polynomial F).roots.toFinset.card := by
          apply Finset.card_le_card; intro a ha
          have ha_eq := (Finset.mem_filter.mp ha).2
          rw [Multiset.mem_toFinset, mem_roots hg_ne, IsRoot,
              eval_add, eval_pow, eval_X, eval_one]
          simp [ha_eq]
      _ ≤ (X ^ m + 1 : Polynomial F).roots.card := Multiset.toFinset_card_le _
      _ ≤ (X ^ m + 1 : Polynomial F).natDegree := card_roots' _
      _ = m := hg_deg
  -- 非零元素覆盖：Fermat → (a^m)² = 1 → a^m = ±1
  have hcover : ∀ a : F, a ≠ 0 → a ∈ A ∨ a ∈ B := by
    intro a ha_ne
    have hfermat : a ^ (q - 1) = 1 := FiniteField.pow_card_sub_one_eq_one a ha_ne
    have hsq : (a ^ m) ^ 2 = 1 := by
      rw [← pow_mul, show m * 2 = q - 1 from by omega]; exact hfermat
    have hprod : (a ^ m - 1) * (a ^ m + 1) = 0 := by
      have : (a ^ m - 1) * (a ^ m + 1) = (a ^ m) ^ 2 - 1 := by ring
      rw [this, hsq, sub_self]
    rcases mul_eq_zero.mp hprod with h | h
    · left; exact Finset.mem_filter.mpr ⟨Finset.mem_univ _, sub_eq_zero.mp h⟩
    · right; exact Finset.mem_filter.mpr ⟨Finset.mem_univ _, eq_neg_of_add_eq_zero_left h⟩
  -- A ∩ B = ∅
  have hAB_disj : Disjoint A B := by
    rw [Finset.disjoint_filter]
    intro a _ h1 h2; exact hone_ne_neg (h1.symm.trans h2)
  -- |A| + |B| ≥ q - 1
  have hAB_ge : q - 1 ≤ A.card + B.card := by
    rw [← Finset.card_union_of_disjoint hAB_disj]
    calc q - 1
        = (Finset.univ.erase (0 : F)).card := by
          rw [Finset.card_erase_of_mem (Finset.mem_univ _), Finset.card_univ]
      _ ≤ (A ∪ B).card := by
          apply Finset.card_le_card; intro a ha
          rw [Finset.mem_erase] at ha
          exact (hcover a ha.1).elim (Finset.mem_union_left B) (Finset.mem_union_right A)
  -- 夹逼：|A| = m
  omega
