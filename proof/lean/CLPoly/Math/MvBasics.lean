/-
  CLPoly/Math/MvBasics.lean — 多变量多项式数学基础

  Wang 算法依赖的 L3 引理：
  - eval 保持 Associated
  - 多变量 content（关于第一个变量）
  - finSuccEquiv 求值分解（引用 Mathlib）
-/
import Mathlib.Algebra.MvPolynomial.Equiv
import Mathlib.RingTheory.Polynomial.Content
import Mathlib.RingTheory.Polynomial.UniqueFactorization
import Mathlib.RingTheory.UniqueFactorizationDomain.GCDMonoid
import Mathlib.Algebra.EuclideanDomain.Int

set_option autoImplicit false

open Polynomial MvPolynomial

-- ============================================================
-- 1. eval 保持 Associated
-- ============================================================

/-- 求值保持 Associated 关系。-/
lemma MvPolynomial.eval_associated {σ R : Type*} [CommSemiring R]
    (α : σ → R) {f g : MvPolynomial σ R}
    (h : Associated f g) : Associated (MvPolynomial.eval α f) (MvPolynomial.eval α g) :=
  h.map (MvPolynomial.eval α)

-- ============================================================
-- 2. NormalizedGCDMonoid instance for MvPolynomial (Fin n) ℤ
-- ============================================================

/-- MvPolynomial (Fin n) ℤ 有 NormalizedGCDMonoid（从 UFD 推导）。
    这使得 Polynomial.content 可用于 Polynomial (MvPolynomial (Fin n) ℤ)。-/
noncomputable instance mvPolyIntNormalizedGCDMonoid (n : ℕ) :
    NormalizedGCDMonoid (MvPolynomial (Fin n) ℤ) := by
  letI : NormalizationMonoid (MvPolynomial (Fin n) ℤ) :=
    UniqueFactorizationMonoid.normalizationMonoid
  exact UniqueFactorizationMonoid.toNormalizedGCDMonoid _

-- ============================================================
-- 3. 多变量 content（关于第一个变量）
-- ============================================================

/-- 多变量多项式关于第一个变量的 content。
    通过 finSuccEquiv 视为 Polynomial (MvPolynomial (Fin n) ℤ)，取 content。-/
noncomputable def mvContentX0 (n : ℕ) (f : MvPolynomial (Fin (n + 1)) ℤ) :
    MvPolynomial (Fin n) ℤ :=
  (MvPolynomial.finSuccEquiv ℤ n f).content

/-- 多变量 IsPrimitive：关于第一个变量的 content 是 unit。-/
def mvIsPrimitiveX0 (n : ℕ) (f : MvPolynomial (Fin (n + 1)) ℤ) : Prop :=
  IsUnit (mvContentX0 n f)

-- ============================================================
-- 4. content 基本性质
-- ============================================================

/-- Gauss 引理：content(f * g) = content(f) * content(g)。-/
lemma mvContentX0_mul (n : ℕ) (f g : MvPolynomial (Fin (n + 1)) ℤ) :
    mvContentX0 n (f * g) = mvContentX0 n f * mvContentX0 n g := by
  simp only [mvContentX0, map_mul, Polynomial.content_mul]
