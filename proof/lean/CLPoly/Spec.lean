/-
  CLPoly/Spec.lean — 子过程接口规约
  Phase 1 核心文件：定义每个子过程"做了什么"的 Prop

  这些规约是验证链的接口层：
  - Phase 1 用它们组装顶层骨架（factor_Zp_correct / factor_ZZ_correct）
  - Phase 2-4 逐个证明具体算法满足这些规约
  - 规约一旦锁定，后续只填充证明，不修改签名
-/
import Mathlib.RingTheory.Polynomial.Basic
import Mathlib.RingTheory.Polynomial.GaussLemma
import Mathlib.RingTheory.Polynomial.Content
import Mathlib.Algebra.Polynomial.FieldDivision
import Mathlib.Algebra.Polynomial.Degree.Defs
import Mathlib.Data.ZMod.Basic
import Mathlib.FieldTheory.Finite.Basic
import Mathlib.FieldTheory.Separable
import Mathlib.Algebra.Squarefree.Basic
import Mathlib.RingTheory.MvPolynomial.Basic

set_option autoImplicit false

open Polynomial

variable {p : ℕ} [Fact (Nat.Prime p)]

-- ============================================================
-- 1. Zp[x] 因式分解子过程规约
-- ============================================================

/-- SQF 规约：输入非零多项式 → 输出无平方因子分解

    result 是 (因子, 重数) 列表，满足：
    1. 乘积还原到 associate（差一个单位）
    2. 每个因子无平方且首一
    3. 重数 ≥ 1
    4. 因子两两互素 -/
def SquarefreeDecomp (f : Polynomial (ZMod p))
    (result : List (Polynomial (ZMod p) × ℕ)) : Prop :=
  -- 1. 乘积还原（到 associate）
  Associated f (result.map (fun pr => pr.1 ^ pr.2)).prod
  -- 2. 每个因子无平方且首一
  ∧ (∀ pr ∈ result, Squarefree pr.1 ∧ Monic pr.1)
  -- 3. 重数 ≥ 1
  ∧ (∀ pr ∈ result, pr.2 ≥ 1)
  -- 4. 因子两两互素
  ∧ (∀ pr₁ ∈ result, ∀ pr₂ ∈ result, pr₁ ≠ pr₂ → IsCoprime pr₁.1 pr₂.1)

/-- DDF 规约：输入首一无平方多项式 → 输出按度分组的不可约因子积

    result 是 (gd, d) 列表，其中 gd 是 f 中所有 d 次不可约因子之积。满足：
    1. 每个 gd 整除 f
    2. gd 的每个不可约因子度数 = d
    3. f 的每个不可约因子都在某个 gd 中
    4. 乘积还原 -/
def DDFCorrect (f : Polynomial (ZMod p))
    (result : List (Polynomial (ZMod p) × ℕ)) : Prop :=
  -- 1. 每个 gd 整除 f
  (∀ pr ∈ result, pr.1 ∣ f)
  -- 2. gd 的每个不可约因子度数 = d
  ∧ (∀ pr ∈ result, ∀ q : Polynomial (ZMod p),
      Irreducible q → q ∣ pr.1 → Polynomial.natDegree q = pr.2)
  -- 3. 不遗漏：f 的每个不可约因子都在某个 gd 中
  ∧ (∀ q : Polynomial (ZMod p), Irreducible q → q ∣ f →
      ∃ pr ∈ result, q ∣ pr.1)
  -- 4. 乘积还原
  ∧ Associated f (result.map (fun pr => pr.1)).prod
  -- 5. 每个 gd 首一（域上 gcd 的自然性质，显式声明以衔接 EDF 前置条件）
  ∧ (∀ pr ∈ result, Monic pr.1)

/-- EDF 规约：输入等度多项式（所有不可约因子度 = d）→ 输出不可约因子列表

    result 是不可约因子列表，满足：
    1. 乘积还原
    2. 每个因子不可约、首一、度 = d -/
def EDFCorrect (g : Polynomial (ZMod p)) (d : ℕ)
    (result : List (Polynomial (ZMod p))) : Prop :=
  -- 1. 乘积还原
  Associated g result.prod
  -- 2. 每个因子不可约、首一、度 = d
  ∧ (∀ q ∈ result, Irreducible q ∧ Monic q ∧ Polynomial.natDegree q = d)

-- ============================================================
-- 2. Z[x] 因式分解子过程规约
-- ============================================================

/-- Hensel 提升规约：模 p 的因式分解 → 模 p^k 的因式分解

    给定 f ∈ Z[x] 和模 p 的因子列表，提升到模 p^k，满足：
    1. 模 p^k 下乘积 ≡ f
    2. 因子数量一致 -/
def HenselCorrect
    (f : Polynomial ℤ) (k : ℕ)
    (factors_mod_p : List (Polynomial (ZMod p)))
    (factors_mod_pk : List (Polynomial (ZMod (p ^ k)))) : Prop :=
  -- 1. 模 p^k 下乘积 ≡ f（通过 intCast 映射）
  Polynomial.map (Int.castRingHom (ZMod (p ^ k))) f = factors_mod_pk.prod
  -- 2. 因子数量一致
  ∧ factors_mod_p.length = factors_mod_pk.length

/-- 因子重组规约：Hensel 因子 → 真正的 Z[x] 不可约因子

    给定 f ∈ Z[x]（本原、无平方），输出不可约因子列表，满足：
    1. 乘积还原
    2. 每个因子不可约 -/
def RecombineCorrect
    (f : Polynomial ℤ)
    (result : List (Polynomial ℤ)) : Prop :=
  -- 1. 乘积还原
  Associated f result.prod
  -- 2. 每个因子不可约
  ∧ (∀ g ∈ result, Irreducible g)

-- ============================================================
-- 3. 顶层因式分解规约
-- ============================================================

/-- Zp[x] 完整因式分解规约：输出 leading coefficient + 不可约因子列表 -/
def FactorZpCorrect (f : Polynomial (ZMod p))
    (lc : ZMod p) (factors : List (Polynomial (ZMod p) × ℕ)) : Prop :=
  -- 1. 乘积还原
  f = Polynomial.C lc * (factors.map (fun pr => pr.1 ^ pr.2)).prod
  -- 2. 每个因子不可约、首一、重数 ≥ 1
  ∧ (∀ pr ∈ factors, Irreducible pr.1 ∧ Monic pr.1 ∧ pr.2 ≥ 1)

/-- Z[x] 完整因式分解规约 -/
def FactorZZCorrect (f : Polynomial ℤ)
    (result : List (Polynomial ℤ)) : Prop :=
  -- 1. 乘积还原
  Associated f result.prod
  -- 2. 每个因子不可约
  ∧ (∀ g ∈ result, Irreducible g)

-- ============================================================
-- 3b. 素数选取规约
-- ============================================================

/-- 素数选取结果。对应 C++ __select_prime（lines 1388-1470）。
    好的素数满足：p ∤ lc(f)，f mod p 度数保持，f mod p 无平方。
    C++ 枚举至多 3 个有效素数，选因子数最少的。-/
structure GoodPrime (f : Polynomial ℤ) (p : ℕ) : Prop where
  /-- p 是素数 -/
  prime : Nat.Prime p
  /-- p 不整除首项系数（度数保持） -/
  lc_nonzero : (f.leadingCoeff : ZMod p) ≠ 0
  /-- f mod p 无平方 -/
  sqfree : Squarefree (Polynomial.map (Int.castRingHom (ZMod p)) f)

/-- 好素数存在。对应 C++ __select_prime 的数学保证。
    f squarefree/ℤ → ∃ Bézout 系数 Sf+Tf'=D (D≠0) → p ∤ D·lc(f) 的素数是好的。-/
theorem good_prime_exists (f : Polynomial ℤ) (hf : f ≠ 0) (hd : 0 < f.natDegree)
    -- Bézout 条件：∃ D ≠ 0, S*f + T*f' = C(D)。
    -- 这是 Squarefree f 在 ℚ[x] 上 IsCoprime 清分母后的结果。
    -- C++ __select_prime 隐含使用此条件（gcd(f,f')=1 mod p 的前提）。
    (D : ℤ) (hD : D ≠ 0)
    (S T : Polynomial ℤ)
    (hbez : S * f + T * Polynomial.derivative f = Polynomial.C D) :
    ∃ p : ℕ, GoodPrime f p := by
  -- 取 p > max(|D|, |lc(f)|) 的素数
  obtain ⟨p, hp_ge, hp_prime⟩ := Nat.exists_infinite_primes
    (max D.natAbs f.leadingCoeff.natAbs + 1)
  refine ⟨p, ⟨hp_prime, ?_, ?_⟩⟩
  · -- lc(f) mod p ≠ 0：p > |lc(f)| → p ∤ lc(f)
    intro h
    have hdvd : (p : ℤ) ∣ f.leadingCoeff := (ZMod.intCast_zmod_eq_zero_iff_dvd _ _).mp
      (by exact_mod_cast h)
    have hlt : f.leadingCoeff.natAbs < p := by omega
    have hne : f.leadingCoeff ≠ 0 := Polynomial.leadingCoeff_ne_zero.mpr hf
    -- (p : ℤ) | lc(f) → p ≤ |lc(f)| → |lc(f)| < p → 矛盾
    obtain ⟨k, hk⟩ := hdvd
    have hk_ne : k ≠ 0 := by intro h; rw [h, mul_zero] at hk; exact hne hk
    have : p ≤ f.leadingCoeff.natAbs := by
      rw [hk, Int.natAbs_mul]; simp
      exact Nat.le_mul_of_pos_right _ (Int.natAbs_pos.mpr hk_ne)
    omega
  · -- f mod p squarefree：Bézout mod p → gcd = 1 → squarefree
    haveI : Fact (Nat.Prime p) := ⟨hp_prime⟩
    -- S*f + T*f' = C(D) → (S mod p)(f mod p) + (T mod p)(f mod p)' = C(D mod p)
    -- D mod p ≠ 0（p > |D|）→ gcd(f mod p, f' mod p) | C(D mod p) → gcd = 1
    -- → squarefree
    -- map Bézout 到 𝔽ₚ：(S mod p)(f mod p) + (T mod p)(f mod p)' = D mod p
    set φ := Int.castRingHom (ZMod p) with hφ
    set fp := Polynomial.map φ f
    -- D mod p ≠ 0（p > |D|）：同 lc 的论证
    have hD_ne : (φ D : ZMod p) ≠ 0 := by
      intro h; rw [hφ] at h
      have hdvd_D : (p : ℤ) ∣ D := (ZMod.intCast_zmod_eq_zero_iff_dvd _ _).mp h
      obtain ⟨k, hk⟩ := hdvd_D
      have hk_ne : k ≠ 0 := by intro h; rw [h, mul_zero] at hk; exact hD hk
      have : p ≤ D.natAbs := by
        rw [hk, Int.natAbs_mul]; simp
        exact Nat.le_mul_of_pos_right _ (Int.natAbs_pos.mpr hk_ne)
      omega
    -- map Bézout 到 𝔽ₚ + derivative_map
    have h_sep : Polynomial.Separable fp := by
      -- Separable = IsCoprime fp fp'。用 Bézout 模约化。
      -- S*f + T*f' = C(D) map 到 𝔽ₚ → map(S)*fp + map(T)*fp' = C(D mod p)
      -- C(D mod p) unit → 除以它 → coprime
      rw [Polynomial.separable_def]
      have hD_unit : IsUnit (φ D : ZMod p) := IsUnit.mk0 _ hD_ne
      -- Bézout mod p
      have key : Polynomial.map φ S * fp + Polynomial.map φ T *
          Polynomial.map φ (Polynomial.derivative f) = Polynomial.C (φ D) := by
        have := congr_arg (Polynomial.map φ) hbez
        simp only [Polynomial.map_add, Polynomial.map_mul, Polynomial.map_C] at this
        exact this
      -- fp' = map(f') = derivative(map f) = derivative fp
      rw [show Polynomial.map φ (Polynomial.derivative f) = Polynomial.derivative fp from
        Polynomial.derivative_map f φ |>.symm] at key
      -- 除以 C(D mod p)
      set Dinv := (φ D)⁻¹
      refine ⟨Polynomial.map φ S * Polynomial.C Dinv,
             Polynomial.map φ T * Polynomial.C Dinv, ?_⟩
      have h1 : Polynomial.map φ S * fp + Polynomial.map φ T * Polynomial.derivative fp =
          Polynomial.C (φ D) := key
      -- 两边右乘 C(Dinv)
      have h2 : (Polynomial.map φ S * fp + Polynomial.map φ T * Polynomial.derivative fp) *
          Polynomial.C Dinv = 1 := by
        rw [h1, ← Polynomial.C_mul, mul_inv_cancel₀ hD_ne, Polynomial.C_1]
      calc Polynomial.map φ S * Polynomial.C Dinv * fp +
          Polynomial.map φ T * Polynomial.C Dinv * Polynomial.derivative fp
          = (Polynomial.map φ S * fp + Polynomial.map φ T * Polynomial.derivative fp) *
            Polynomial.C Dinv := by ring
        _ = 1 := h2
    exact h_sep.squarefree

-- ============================================================
-- 3c. QQ[x] 因式分解规约
-- ============================================================

set_option maxHeartbeats 1600000 in
/-- Gauss 引理（ℤ → ℚ 方向）：Z[x] 不可约 → Q[x] 不可约。-/
private theorem irreducible_map_of_irreducible_zz (g : Polynomial ℤ)
    (hg : Irreducible g) (hd : 0 < g.natDegree) :
    Irreducible (Polynomial.map (Int.castRingHom ℚ) g) := by
  have hg_prim := hg.isPrimitive hd.ne'
  exact (hg_prim.irreducible_iff_irreducible_map_fraction_map (K := ℚ)).mp hg

theorem factor_QQ_of_factor_ZZ
    (f : Polynomial ℚ) (hf : f ≠ 0)
    -- ZZ 版本的因式分解正确
    (f_zz : Polynomial ℤ) (lcd : ℤ) (hlcd : lcd ≠ 0)
    (h_convert : Polynomial.map (Int.castRingHom ℚ) f_zz = Polynomial.C (lcd : ℚ) * f)
    (result_zz : List (Polynomial ℤ))
    (h_zz_correct : RecombineCorrect f_zz result_zz)
    -- C++ 输出因子均非常数（常数吸收到 content）
    (h_nonconstant : ∀ g ∈ result_zz, 0 < g.natDegree)
    -- 结论：QQ 因式分解正确
    : ∃ result_qq : List (Polynomial ℚ),
        Associated f result_qq.prod ∧ ∀ g ∈ result_qq, Irreducible g := by
  -- 将 ZZ 因子映射到 QQ
  set result_qq := result_zz.map (Polynomial.map (Int.castRingHom ℚ)) with hres
  refine ⟨result_qq, ?_, ?_⟩
  · -- Associated f result_qq.prod
    -- h_zz_correct.1 : Associated f_zz result_zz.prod
    -- map 保持 Associated → map(f_zz) ∼ map(result_zz.prod) = result_qq.prod
    -- map(f_zz) = lcd * f → f ∼ result_qq.prod (除以 lcd)
    -- h_zz_correct.1 : Associated f_zz result_zz.prod
    -- map 保持 Associated
    have h_map_assoc : Associated (Polynomial.map (Int.castRingHom ℚ) f_zz)
        (Polynomial.map (Int.castRingHom ℚ) result_zz.prod) :=
      h_zz_correct.1.map (Polynomial.mapRingHom (Int.castRingHom ℚ))
    -- map(result_zz.prod) = result_qq.prod
    have h_map_prod : Polynomial.map (Int.castRingHom ℚ) result_zz.prod = result_qq.prod := by
      simp only [hres, Polynomial.map_list_prod, List.map_map]
    -- map(f_zz) = C(lcd) * f
    rw [h_convert] at h_map_assoc
    rw [h_map_prod] at h_map_assoc
    -- C(lcd) * f ∼ result_qq.prod → f ∼ result_qq.prod（C(lcd) 是 ℚ[x] 中的 unit）
    have hlcd_unit : IsUnit (Polynomial.C (lcd : ℚ)) :=
      Polynomial.isUnit_C.mpr ⟨⟨(lcd : ℚ), (lcd : ℚ)⁻¹,
        mul_inv_cancel₀ (Int.cast_ne_zero.mpr hlcd),
        inv_mul_cancel₀ (Int.cast_ne_zero.mpr hlcd)⟩, rfl⟩
    exact (associated_isUnit_mul_left_iff hlcd_unit).mp h_map_assoc
  · -- ∀ g ∈ result_qq, Irreducible g
    -- 每个 g = map(gᵢ) where gᵢ ∈ result_zz
    -- gᵢ 在 Z[x] 不可约 (h_zz_correct.2)
    -- Gauss 引理 (IsPrimitive.irreducible_iff_irreducible_map_fraction_map):
    --   gᵢ 本原不可约 → map(gᵢ) 在 Q[x] 不可约
    intro g hg
    obtain ⟨g_zz, hg_mem, rfl⟩ := List.mem_map.mp hg
    have hg_irr := h_zz_correct.2 g_zz hg_mem
    -- gᵢ 不可约 → gᵢ 本原 (Irreducible → IsPrimitive for nonconstant)
    -- Gauss 引理：Z[x] 不可约 → Q[x] 不可约
    exact irreducible_map_of_irreducible_zz g_zz hg_irr (h_nonconstant g_zz hg_mem)

-- ============================================================
-- 4. 多变量因式分解规约
-- ============================================================

/-- Z[x₁,...,xₙ] 多变量完整因式分解规约。
    结构与 FactorZZCorrect 相同，作用于 MvPolynomial。
    对应 C++: polynomial_factorize_wang.hh __factor_multivar -/
def MvFactorCorrect {σ : Type*}
    (f : MvPolynomial σ ℤ)
    (result : List (MvPolynomial σ ℤ)) : Prop :=
  -- 1. 乘积还原
  Associated f result.prod
  -- 2. 每个因子不可约
  ∧ (∀ g ∈ result, Irreducible g)
