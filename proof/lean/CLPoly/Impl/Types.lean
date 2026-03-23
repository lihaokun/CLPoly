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
import Mathlib.Algebra.Polynomial.Degree.Definitions

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
  h_sorted : terms.Chain' (fun a b => a.1 > b.1)
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

/-- 非空多项式非零（TODO：需要从不变量推导） -/
theorem SparsePolyZp.toPoly_ne_zero_of_nonempty (sp : SparsePolyZp)
    (h : sp.terms ≠ []) :
    sp.toPoly ≠ 0 := by
  sorry -- 需要证明：非零系数 + 严格降序 → toPoly ≠ 0

/-- deg 精化：SparsePolyZp.deg = natDegree(toPoly)（TODO） -/
theorem SparsePolyZp.deg_eq_natDegree (sp : SparsePolyZp) :
    sp.deg = sp.toPoly.natDegree := by
  sorry -- 需要证明：降序不变量 + 非零 → 度数匹配
