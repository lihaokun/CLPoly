/-
  CLPoly/Spec.lean — 子过程接口规约
  Phase 1 核心文件：定义每个子过程"做了什么"的 Prop

  这些规约是验证链的接口层：
  - Phase 1 用它们组装顶层骨架（factor_Zp_correct / factor_ZZ_correct）
  - Phase 2-4 逐个证明具体算法满足这些规约
  - 规约一旦锁定，后续只填充证明，不修改签名
-/
import Mathlib.RingTheory.Polynomial.Basic
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
