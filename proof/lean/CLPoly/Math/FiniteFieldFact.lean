/-
  CLPoly/Math/FiniteFieldFact.lean — L3 有限域数学定理

  Phase 2 数学基石：DDF/EDF 正确性所需的纯数学定理。

  T2.1: X^{p^d} - X 的 Separable 性              ✓ 无 sorry
  T2.2: 不可约多项式整除 X^{p^d} - X             ✓ 无 sorry
  T2.2': T2.2 的逆命题                            ✓ 无 sorry
  T2.3: gcd(X^{p^d}-X, f) 的不可约因子刻画       ✓ 无 sorry
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
-- T2.2': 逆命题 — q | X^{p^d} - X → deg q | d
-- ============================================================

/-- T2.2 的逆命题：不可约 g 整除 X^{p^d} - X → deg(g) | d。
    证明：
    1. 在 K = AdjoinRoot g 中，α^{p^d} = α
    2. 迭代 Frobenius 作为 AlgHom 固定生成元 → 固定全域 → ∀ x, x^{p^d} = x
    3. 若 k ∤ d，则 ∀ x, x^{p^r} = x（r = d mod k），但 X^{p^r}-X 的度
       p^r < p^k = |K|，根计数矛盾。 -/
theorem dvd_X_pow_sub_X_imp_natDegree_dvd
    (g : Polynomial (ZMod p)) (hg : Irreducible g)
    (d : ℕ) (h : g ∣ (X ^ (p ^ d) - X : Polynomial (ZMod p))) :
    g.natDegree ∣ d := by
  haveI : Fact (Irreducible g) := ⟨hg⟩
  haveI : Nontrivial (AdjoinRoot g) := inferInstance
  haveI : CharP (AdjoinRoot g) p :=
    charP_of_injective_algebraMap (algebraMap (ZMod p) (AdjoinRoot g)).injective p
  letI : Fintype (AdjoinRoot g) :=
    Fintype.ofEquiv _ (AdjoinRoot.powerBasis hg.ne_zero).basis.equivFun.toEquiv.symm
  set k := g.natDegree
  have hk_pos : 0 < k := Irreducible.natDegree_pos hg
  have hp : Nat.Prime p := Fact.out
  have hcard : Fintype.card (AdjoinRoot g) = p ^ k := by
    rw [Fintype.card_congr (AdjoinRoot.powerBasis hg.ne_zero).basis.equivFun.toEquiv,
        Fintype.card_fun, Fintype.card_fin, ZMod.card]; rfl
  -- Step A: α^{p^d} = α
  have hα : (AdjoinRoot.root g) ^ (p ^ d) = AdjoinRoot.root g := by
    have h0 := AdjoinRoot.mk_eq_zero.mpr h
    simp only [map_sub, map_pow, AdjoinRoot.mk_X] at h0
    rwa [sub_eq_zero] at h0
  -- Step B: ∀ x, x^{p^d} = x (Frobenius AlgHom + AdjoinRoot.algHom_ext)
  have hall : ∀ (x : AdjoinRoot g), x ^ (p ^ d) = x := by
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
