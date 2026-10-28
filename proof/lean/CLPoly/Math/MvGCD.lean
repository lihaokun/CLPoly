/-
  CLPoly/Math/MvGCD.lean — 多变量 GCD / cont / leadcoeff 数学 spec

  Phase F-impl-A.1: 用 Mathlib `UniqueFactorizationMonoid` + `GCDMonoid` 给
  MvPolyZZ 的 GCD / cont / leadcoeff 定义 noncomputable 数学 spec。

  目的：
  - 为后续 L2 算法证明（refinement）准备数学层接口
  - 执行层（L1）由 Phase F-impl-A.2 的 FFI 提供
  - refinement 定理形态：`gcd_impl = gcd_spec`（用 FFI/手写 impl 对照数学 spec）

  设计：
  - 桥 `MvPolyZZ → MvPolynomial Variable ℤ`（σ = Variable = UInt64，σ 无限但
    每个具体 poly support 有限，Mathlib UFD 实例适用）
  - 用 Mathlib `UniqueFactorizationMonoid.toGCDMonoid`（noncomputable）拿 gcd
  - 反桥用 `toFinset` 枚举 support
-/
import Mathlib.Algebra.MvPolynomial.Basic
import Mathlib.Algebra.MvPolynomial.Equiv
import Mathlib.Algebra.EuclideanDomain.Int
import Mathlib.RingTheory.PrincipalIdealDomain
import Mathlib.RingTheory.Polynomial.UniqueFactorization
import Mathlib.RingTheory.UniqueFactorizationDomain.GCDMonoid
import CLPoly.Model

set_option autoImplicit false

open MvPolynomial

namespace CLPoly.Math.MvGCD

-- ============================================================
-- §1. Variable 类型适配：MvPolyZZ 的 Variable = UInt64
-- ============================================================

-- σ 选 Variable（= UInt64），Mathlib MvPolynomial 接受无限 σ 类型；
-- 每个具体 poly 的 support 是 Finset，所以 UFD 性质保持（用
-- `exists_finset_rename` reduction）。

abbrev MvPolyZZSpec : Type := MvPolynomial Variable ℤ

-- ============================================================
-- §2. Monomial → Finsupp 桥
-- ============================================================

/-- 把 `Monomial = Array (Variable × Int64)` 转成 Mathlib 的 `Variable →₀ ℕ`。
    负 exp（理论不该出现）按 0 截断。-/
noncomputable def monomialToFinsupp (m : Monomial) : Variable →₀ ℕ :=
  m.foldr (fun (vp : Variable × Int64) acc =>
    acc + Finsupp.single vp.fst vp.snd.toNatClampNeg) (0 : Variable →₀ ℕ)

-- ============================================================
-- §3. MvPolyZZ → MvPolyZZSpec 桥（toMathlib）
-- ============================================================

/-- 把 Model 的 `MvPolyZZ = Array (Monomial × Int)` 转成 Mathlib
    `MvPolynomial Variable ℤ`。逐项构建：∑ cᵢ · X^mᵢ。-/
noncomputable def toMathlib (f : MvPolyZZ) : MvPolyZZSpec :=
  f.foldr (fun term acc => acc + MvPolynomial.monomial (monomialToFinsupp term.fst) term.snd)
    (0 : MvPolyZZSpec)

-- ============================================================
-- §4. MvPolyZZSpec → MvPolyZZ 反桥（fromMathlib）
-- ============================================================

/-- Finsupp 反桥到 Monomial：枚举 support。 -/
noncomputable def finsuppToMonomial (s : Variable →₀ ℕ) : Monomial :=
  let l := s.support.toList
  let arr : Array (Variable × Int64) :=
    l.toArray.map (fun v => (v, (s v : Int64)))
  arr

/-- MvPolynomial 反桥：枚举 support，每项 `(monomial, coeff)`，组合成
    `MvPolyZZ`。 -/
noncomputable def fromMathlib (p : MvPolyZZSpec) : MvPolyZZ :=
  let l := p.support.toList
  let arr : MvPolyZZ :=
    l.toArray.map (fun s => (finsuppToMonomial s, p.coeff s))
  arr

-- ============================================================
-- §5. GCD spec（基于 Mathlib UFD → GCDMonoid）
-- ============================================================

/-- `MvPolyZZSpec = MvPolynomial Variable ℤ` 通过 Mathlib UFD 得到 GCDMonoid。
    noncomputable，但概念上良定义。 -/
noncomputable instance mvPolyZZSpecGCDMonoid : GCDMonoid MvPolyZZSpec :=
  UniqueFactorizationMonoid.toGCDMonoid MvPolyZZSpec

/-- 多变量 GCD 数学 spec：通过桥到 Mathlib，调 GCDMonoid.gcd，再桥回。 -/
noncomputable def gcd_spec (f g : MvPolyZZ) : MvPolyZZ :=
  fromMathlib (GCDMonoid.gcd (toMathlib f) (toMathlib g))

-- ============================================================
-- §6. cont / leadcoeff spec（w.r.t. 主变量 = 第一个变量）
-- ============================================================

/-- 多变量 cont 数学 spec：关于第一个变量。具体 Mathlib 用 `finSuccEquiv`
    视为 `Polynomial (MvPolynomial _ ℤ)` 后取 `Polynomial.content`。
    简化版：直接 GCD 系数（视为 Mathlib 多项式系数集合）。 -/
noncomputable def cont_spec (f : MvPolyZZ) : MvPolyZZ :=
  -- 通过 Mathlib 桥，对 toMathlib f 的系数集合取 gcd
  -- 暂用占位的 toMathlib + 反桥；精确 finSuccEquiv 路径见 MvBasics.lean
  fromMathlib (toMathlib f)

/-- 多变量 leadcoeff 数学 spec：关于第一个变量取首项系数（在剩余变量中的
    多项式）。 -/
noncomputable def leadcoeff_spec (f : MvPolyZZ) : MvPolyZZ :=
  -- 同 cont_spec：暂用桥+反桥占位；具体 finSuccEquiv 路径留 refinement 阶段细化
  fromMathlib (toMathlib f)

-- ============================================================
-- §7. Spec 基本性质（refinement 阶段填充）
-- ============================================================

-- 性质 `GCDMonoid.gcd (toMathlib f) (toMathlib g) = toMathlib (gcd_spec f g)`
-- 等需要先证 `fromMathlib ∘ toMathlib = id`（同构性质）。
-- refinement 阶段（Phase F-proof）填充。

end CLPoly.Math.MvGCD
