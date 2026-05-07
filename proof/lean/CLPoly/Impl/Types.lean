/-
  CLPoly/Impl/Types.lean — L1 实现模型：基础类型定义

  Phase L1.0: C++ 类型的 Lean 模型
  对应 C++: number.hh (Zp class), upolynomial.hh (upolynomial_<Zp>)

  设计决策：
  - 系数用 Nat（∈ [0, p)），不单独建模 Zp 类
  - 度数用 Nat（简化，不处理 > 2^64）
  - 素数在 SparsePolyZp 结构级共享
  - Barrett 预计算省略（不影响数学正确性）
-/
import Mathlib.Data.ZMod.Basic
import Mathlib.Algebra.Polynomial.Basic
import Mathlib.Algebra.Polynomial.Coeff
import Mathlib.Algebra.Polynomial.Degree.Defs
import Mathlib.Algebra.Polynomial.Degree.Operations
import Mathlib.Data.List.Chain

set_option autoImplicit false

open Polynomial

-- ============================================================
-- 1. Zp 算术（内联操作，对应 C++ Zp 类的运算）
-- ============================================================

/-- 模加法（对应 C++ Zp::operator+）。
    C++ 实现用 `p - a > b ? a + b : a + b - p` 避免溢出。
    L1 用 Nat 回避溢出问题，数学语义等价。 -/
def zp_add (a b p : Nat) : Nat :=
  if a + b < p then a + b else a + b - p

/-- 模减法（对应 C++ Zp::operator-）。 -/
def zp_sub (a b p : Nat) : Nat :=
  if a ≥ b then a - b else p - b + a

/-- 模乘法（对应 C++ Zp::__nmod_mul，Barrett 省略）。 -/
def zp_mul (a b p : Nat) : Nat :=
  (a * b) % p

/-- 模取负（对应 C++ Zp::operator- (unary)）。 -/
def zp_neg (a p : Nat) : Nat :=
  if a = 0 then 0 else p - a

-- ============================================================
-- 2. Zp 算术精化定理
-- ============================================================

theorem zp_add_spec (a b p : Nat) (hp : 0 < p) (ha : a < p) (hb : b < p) :
    zp_add a b p = (a + b) % p := by
  unfold zp_add; split
  · exact (Nat.mod_eq_of_lt ‹_›).symm
  · rename_i h; push_neg at h
    conv_rhs => rw [show a + b = (a + b - p) + 1 * p from by omega]
    rw [Nat.add_mul_mod_self_right, Nat.mod_eq_of_lt (by omega)]

theorem zp_sub_spec (a b p : Nat) (hp : 0 < p) (ha : a < p) (hb : b < p) :
    zp_sub a b p = (a + (p - b)) % p := by
  unfold zp_sub; split
  · rename_i h
    conv_rhs => rw [show a + (p - b) = (a - b) + 1 * p from by omega]
    rw [Nat.add_mul_mod_self_right, Nat.mod_eq_of_lt (by omega)]
  · rename_i h; push_neg at h
    rw [show p - b + a = a + (p - b) from by omega, Nat.mod_eq_of_lt (by omega)]

theorem zp_mul_spec (a b p : Nat) : zp_mul a b p = (a * b) % p := rfl

theorem zp_neg_spec (a p : Nat) (hp : 0 < p) (ha : a < p) :
    zp_neg a p = (p - a) % p := by
  simp only [zp_neg]
  split
  · rename_i h; subst h; simp [Nat.mod_self, hp]
  · rename_i h; rw [Nat.mod_eq_of_lt (by omega)]

-- ============================================================
-- 3. SparsePolyZp（稀疏多项式，1:1 对应 C++ upolynomial_<Zp>）
-- ============================================================

/-- 稀疏多项式 over Z/pZ。1:1 对应 C++ 的 upolynomial_<Zp>。
    terms: (度数, 系数值) 对，降序排列。
    不变量：度数严格降序 + 系数非零 + 系数在 [0, p) 内。 -/
structure SparsePolyZp where
  terms : List (Nat × Nat)
  prime : Nat
  h_prime : Nat.Prime prime
  h_sorted : terms.IsChain (fun a b => a.1 > b.1)
  h_nonzero : ∀ t ∈ terms, t.2 ≠ 0
  h_range : ∀ t ∈ terms, t.2 < prime

-- ============================================================
-- 4. SparsePolyZp 便利方法
-- ============================================================

/-- 多项式度数（空多项式 = 0）-/
def SparsePolyZp.deg (sp : SparsePolyZp) : Nat :=
  match sp.terms with
  | [] => 0
  | (d, _) :: _ => d

/-- 是否为空（零多项式） -/
def SparsePolyZp.isEmpty (sp : SparsePolyZp) : Bool :=
  sp.terms.isEmpty

/-- 首项系数（空多项式 = 0） -/
def SparsePolyZp.leadCoeff (sp : SparsePolyZp) : Nat :=
  match sp.terms with
  | [] => 0
  | (_, c) :: _ => c

/-- 是否 monic（首项系数 = 1） -/
def SparsePolyZp.isMonic (sp : SparsePolyZp) : Bool :=
  sp.leadCoeff = 1

-- ============================================================
-- 5. 精化关系：SparsePolyZp → Polynomial (ZMod p)
-- ============================================================

/-- 将 L1 稀疏多项式转换为 L2 数学多项式。
    逐项构建：∑ C(cᵢ) · X^{dᵢ}。 -/
noncomputable def SparsePolyZp.toPoly (sp : SparsePolyZp) :
    Polynomial (ZMod sp.prime) :=
  sp.terms.foldr
    (fun (d, c) acc => acc + Polynomial.C (c : ZMod sp.prime) * Polynomial.X ^ d)
    0

-- ============================================================
-- 6. 基本精化引理
-- ============================================================

/-- 空多项式对应零多项式 -/
theorem SparsePolyZp.toPoly_empty (sp : SparsePolyZp) (h : sp.terms = []) :
    sp.toPoly = 0 := by
  simp [toPoly, h]

-- 辅助引理：高于 n 的"严格小于 n 的所有项度数"foldr 系数为 0
private lemma foldr_at_higher_eq_zero
    {prime : ℕ} (terms : List (Nat × Nat)) (n : Nat)
    (hn : ∀ t ∈ terms, t.1 < n) :
    (terms.foldr
      (fun (dc : Nat × Nat) acc => acc + Polynomial.C ((dc.2 : ZMod prime)) * Polynomial.X ^ dc.1)
      0).coeff n = 0 := by
  induction terms with
  | nil => simp
  | cons hd tl ih =>
    have hhd : hd.1 < n := hn hd (List.mem_cons.mpr (Or.inl rfl))
    have htl : ∀ t ∈ tl, t.1 < n := fun t ht => hn t (List.mem_cons_of_mem _ ht)
    simp only [List.foldr_cons, Polynomial.coeff_add, Polynomial.coeff_C_mul,
      Polynomial.coeff_X_pow]
    rw [ih htl, if_neg (by omega : n ≠ hd.1)]
    ring

-- 辅助引理：Chain'(>) 在 head :: tail 上 → head.1 严格大于所有 tail 元素的 .1
private lemma chain'_gt_head_gt_all
    (head : Nat × Nat) (tail : List (Nat × Nat))
    (h : List.IsChain (fun a b : Nat × Nat => a.1 > b.1) (head :: tail)) :
    ∀ t ∈ tail, t.1 < head.1 := by
  induction tail generalizing head with
  | nil => intros _ ht; exact absurd ht (by simp)
  | cons mid rest ih =>
    rw [List.isChain_cons_cons] at h
    obtain ⟨hhd_mid, hchain⟩ := h
    intro t ht
    rcases List.mem_cons.mp ht with rfl | ht'
    · exact hhd_mid
    · have : t.1 < mid.1 := ih mid hchain t ht'
      omega

/-- 非空多项式非零：head 的系数在 ZMod p 上非零，故 toPoly ≠ 0 -/
theorem SparsePolyZp.toPoly_ne_zero_of_nonempty (sp : SparsePolyZp)
    (h : sp.terms ≠ []) :
    sp.toPoly ≠ 0 := by
  match hterms : sp.terms with
  | [] => exact absurd hterms h
  | head :: tail =>
    intro hzero
    have h_head_mem : head ∈ sp.terms := hterms ▸ List.mem_cons.mpr (Or.inl rfl)
    have h_head_nonzero : head.2 ≠ 0 := sp.h_nonzero head h_head_mem
    have h_head_range : head.2 < sp.prime := sp.h_range head h_head_mem
    have h_tail_lt : ∀ t ∈ tail, t.1 < head.1 := by
      have hsorted := sp.h_sorted
      rw [hterms] at hsorted
      exact chain'_gt_head_gt_all head tail hsorted
    -- 系数 head.2 在 ZMod p 中非零（因为 0 < head.2 < prime 且 prime 是素数）
    have h_cast_ne : (head.2 : ZMod sp.prime) ≠ 0 := by
      haveI : Fact (Nat.Prime sp.prime) := ⟨sp.h_prime⟩
      rw [ne_eq, ZMod.natCast_eq_zero_iff]
      intro hdvd
      exact absurd (Nat.le_of_dvd (Nat.pos_of_ne_zero h_head_nonzero) hdvd)
        (Nat.not_le_of_lt h_head_range)
    -- toPoly.coeff head.1 = head.2（在 ZMod p 上）
    have h_coeff : sp.toPoly.coeff head.1 = (head.2 : ZMod sp.prime) := by
      unfold toPoly
      rw [hterms, List.foldr_cons, Polynomial.coeff_add,
          foldr_at_higher_eq_zero tail head.1 h_tail_lt,
          Polynomial.coeff_C_mul, Polynomial.coeff_X_pow, if_pos rfl]
      ring
    rw [hzero, Polynomial.coeff_zero] at h_coeff
    exact h_cast_ne h_coeff.symm

/-- deg 精化：SparsePolyZp.deg = natDegree(toPoly) -/
theorem SparsePolyZp.deg_eq_natDegree (sp : SparsePolyZp) :
    sp.deg = sp.toPoly.natDegree := by
  match hterms : sp.terms with
  | [] =>
    simp only [SparsePolyZp.deg, hterms]
    rw [sp.toPoly_empty hterms]
    simp
  | head :: tail =>
    simp only [SparsePolyZp.deg, hterms]
    -- 显示 natDegree = head.1
    have h_head_mem : head ∈ sp.terms := hterms ▸ List.mem_cons.mpr (Or.inl rfl)
    have h_head_nonzero : head.2 ≠ 0 := sp.h_nonzero head h_head_mem
    have h_head_range : head.2 < sp.prime := sp.h_range head h_head_mem
    have h_tail_lt : ∀ t ∈ tail, t.1 < head.1 := by
      have hsorted := sp.h_sorted
      rw [hterms] at hsorted
      exact chain'_gt_head_gt_all head tail hsorted
    have h_cast_ne : (head.2 : ZMod sp.prime) ≠ 0 := by
      haveI : Fact (Nat.Prime sp.prime) := ⟨sp.h_prime⟩
      rw [ne_eq, ZMod.natCast_eq_zero_iff]
      intro hdvd
      exact absurd (Nat.le_of_dvd (Nat.pos_of_ne_zero h_head_nonzero) hdvd)
        (Nat.not_le_of_lt h_head_range)
    have h_coeff_head : sp.toPoly.coeff head.1 = (head.2 : ZMod sp.prime) := by
      unfold toPoly
      rw [hterms, List.foldr_cons, Polynomial.coeff_add,
          foldr_at_higher_eq_zero tail head.1 h_tail_lt,
          Polynomial.coeff_C_mul, Polynomial.coeff_X_pow, if_pos rfl]
      ring
    -- 上界：natDegree ≤ head.1（高于 head.1 的系数都为 0）
    have h_deg_le : sp.toPoly.degree ≤ (head.1 : WithBot ℕ) := by
      rw [Polynomial.degree_le_iff_coeff_zero]
      intro m hm
      have h_m_gt : head.1 < m := by exact_mod_cast hm
      unfold toPoly
      rw [hterms]
      apply foldr_at_higher_eq_zero
      intro t ht
      rcases ht with _ | ⟨_, ht'⟩
      · exact h_m_gt
      · exact lt_trans (h_tail_lt t ht') h_m_gt
    have h_natDeg_le : sp.toPoly.natDegree ≤ head.1 :=
      Polynomial.natDegree_le_of_degree_le h_deg_le
    -- 下界：head.1 ≤ natDegree（contrapositive of coeff_eq_zero_of_natDegree_lt）
    have h_natDeg_ge : head.1 ≤ sp.toPoly.natDegree := by
      by_contra h_lt
      push_neg at h_lt
      have h_zero : sp.toPoly.coeff head.1 = 0 :=
        Polynomial.coeff_eq_zero_of_natDegree_lt h_lt
      rw [h_coeff_head] at h_zero
      exact h_cast_ne h_zero
    omega
