# 多变量因式分解 L3 数学基石

> 状态：nl-proof v1
> 目标：Wang 算法依赖的数学基础引理

---

## 0. 依赖的 Mathlib API（调研确认）

| API | 状态 | 路径 |
|-----|------|------|
| `MvPolynomial.eval` | ✅ | `Algebra.MvPolynomial.Eval` |
| `MvPolynomial.finSuccEquiv` | ✅ | `Algebra.MvPolynomial.Equiv` |
| `eval₂_dvd` | ✅ | `Algebra.MvPolynomial.Eval` |
| `eval_eq_eval_mv_eval'` | ✅ | `Algebra.MvPolynomial.Equiv` |
| `natDegree_finSuccEquiv` | ✅ | `Algebra.MvPolynomial.Equiv` |
| `Polynomial.content` | ✅ | `RingTheory.Polynomial.Content` |
| `MvPolynomial.content` | ❌ | 需自定义 |

---

## 1. 求值保持乘积和整除

### 1.1 已有（Mathlib 直接提供）

```
eval_mul : eval α (f * g) = eval α f * eval α g     -- 环同态
eval₂_dvd : f ∣ g → eval₂ φ ψ f ∣ eval₂ φ ψ g      -- 保持整除
eval_prod : eval α (∏ fᵢ) = ∏ (eval α fᵢ)            -- 保持有限积
```

这些是 Mathlib 已有的，**不需要证明**。

### 1.2 需要证明：求值保持 Associated

```lean
/-- 求值保持 Associated 关系 -/
lemma eval_associated {σ R : Type*} [CommRing R]
    (α : σ → R) (f g : MvPolynomial σ R)
    (h : Associated f g) : Associated (eval α f) (eval α g)
```

**证明**：Associated f g → ∃ u : Rˣ, f = g * u → eval α f = eval α g * eval α u。
eval α u 是 unit（eval 是环同态，保持 unit）。✓

**Lean 路径**：`Associated.map (eval α)` 或 `h.map (MvPolynomial.eval α)`。
RingHom.map 保持 Associated（Mathlib 有 `Associated.map`）。

---

## 2. 多变量 content（关于一个变量）

### 2.1 定义

Wang 算法需要 "content of f with respect to x₁" = gcd of coefficients when f is viewed as polynomial in x₁。

通过 `finSuccEquiv`：
```
f ∈ MvPolynomial (Fin (n+1)) R  →  finSuccEquiv R n f ∈ Polynomial (MvPolynomial (Fin n) R)
```

然后用 `Polynomial.content` 对 `finSuccEquiv R n f` 取 content。

```lean
/-- 多变量多项式关于第一个变量的 content -/
noncomputable def mvContent {n : ℕ} (f : MvPolynomial (Fin (n + 1)) ℤ) :
    MvPolynomial (Fin n) ℤ :=
  (finSuccEquiv ℤ n f).content
```

### 2.2 性质

**需要确认**：`Polynomial.content` 要求系数环是 `GCDMonoid`。
`MvPolynomial (Fin n) ℤ` 是 UFD → 是 `GCDMonoid`？

Mathlib 中 `UniqueFactorizationMonoid` 继承 `GCDMonoid`（通过 `WfDvdMonoid.gcdMonoid`）。
`MvPolynomial (Fin n) ℤ` 是 UFD（§调研确认）→ 有 `GCDMonoid` 实例。✓

**但**：`GCDMonoid` 实例可能是 noncomputable（需要 `Classical.choice`）。这在证明中没问题，但需要 `noncomputable` 标注。

### 2.3 content 的关键引理

```lean
/-- Gauss 引理：content(f * g) = content(f) * content(g)（up to associated）-/
-- Mathlib: Polynomial.content_mul（需要系数环是 GCDMonoid + 域条件...）
-- 需要确认 Mathlib 中 Gauss 引理的精确条件
```

**注意**：Mathlib 的 `Polynomial.content_mul` 可能需要特定条件。需要进一步调研。

---

## 3. 求值与因子数

### 3.1 Wang 条件 (d) 的数学内容

**定理**：设 f ∈ Z[x₁,...,xₙ] 本原无平方，α = (α₂,...,αₙ)。
若 f₀ = f(x₁,α₂,...,αₙ) ∈ Z[x₁] 无平方且 lc(f,x₁)(α) ≠ 0，
则 f₀ 的不可约因子数 ≥ f 的不可约因子数。

**更强版本**（Hilbert 不可约性）：对"几乎所有" α，不可约因子数相等。

### 3.2 形式化路径

完整的 Hilbert 不可约性定理在 Mathlib 中可能没有。但我们只需要**弱方向**：

```lean
/-- 求值不减少不可约因子数（弱方向）：
    若 f = g₁ · ... · gₜ 在 MvPolynomial 中不可约分解，
    则 eval α f 在 Polynomial 中的不可约因子数 ≥ t。

    直觉：eval 是环同态，gᵢ(x₁,α) 可能不再不可约，但 eval 不会"创造"新因子。
    实际上方向相反：eval 可能合并因子（gᵢ(α) 和 gⱼ(α) 可能有公共因子）。
    所以不可约因子数可能增多（分裂）也可能减少（合并）。-/
```

**实际上**：方向分析——
- f = g₁ · g₂ 在 MvPoly 中
- f(x₁,α) = g₁(x₁,α) · g₂(x₁,α)
- g₁(x₁,α) 可能可约 → 因子数增多 ✓
- g₁(x₁,α) 和 g₂(x₁,α) 可能有公共因子 → 因子数减少 ✗

所以**求值可以增加也可以减少因子数**。Wang 条件 (d) 要求选择 α 使得因子数**恰好**等于 r（单变量因子数 = 多变量因子数）。这是一个**存在性**条件（好的 α 存在），不是对所有 α 成立。

### 3.3 务实方案

条件 (d) 的证明非常深（本质上是 Hilbert 不可约性定理）。**在 L2 中，我们将条件 (d) 作为 EvalPointGood 的前置假设**，不在 L3 证明它。

这类似于单变量 EDF 中的 `splits_fn` 假设——概率性条件作为外部参数。

---

## 4. finSuccEquiv 的求值分解

### 4.1 已有关键定理

Mathlib 中 `eval_eq_eval_mv_eval'`：

```
eval (Fin.cons y s) f = Polynomial.eval y (Polynomial.map (eval s) (finSuccEquiv R n f))
```

即：将 f ∈ MvPoly(Fin(n+1), R) 在 (y, s₁,...,sₙ) 处求值 =
先用 finSuccEquiv 转为 Poly(MvPoly(Fin n, R))，然后 map(eval s) 得到 Poly(R)，最后 eval y。

### 4.2 应用

这直接给出了 Wang 的"求值降维"步骤的数学正当性：
- f ∈ MvPoly(Fin(n+1), ℤ)
- 选 α = (α₂,...,αₙ₊₁)
- f₀ = f(x₁,α₂,...,αₙ₊₁) = Poly.eval_at_α (finSuccEquiv ℤ n f)

不需要额外证明，Mathlib 已有。✓

---

## 5. 总结：需要新证明的 L3 定理

| 定理 | 估计行数 | 说明 |
|------|---------|------|
| `eval_associated` | ~5 | Associated.map |
| `mvContent` 定义 | ~5 | finSuccEquiv + Polynomial.content |
| `mvContent_mul`（Gauss 引理） | ~20 | 需确认 Mathlib 的 Polynomial.content_mul 条件 |
| `mvIsPrimitive` 定义 | ~5 | content = 1 |
| **总计** | **~35** |

**不需要证明**（Mathlib 已有或作为假设）：
- eval 保持乘积/整除（Mathlib）
- finSuccEquiv 求值分解（Mathlib）
- 条件 (d)（作为 EvalPointGood 假设）
- Hilbert 不可约性（不在 L3 范围内）

---

## 6. 形式化计划

```lean
-- 文件：CLPoly/Math/MvPolynomialBasics.lean
import Mathlib.Algebra.MvPolynomial.Equiv
import Mathlib.RingTheory.Polynomial.Content
import Mathlib.RingTheory.Polynomial.UniqueFactorization

-- 1. eval 保持 Associated
lemma eval_associated ...

-- 2. mvContent 定义
noncomputable def mvContent ...

-- 3. mvIsPrimitive 定义
def mvIsPrimitive ...

-- 4. Gauss 引理（如果 Mathlib 条件匹配）
lemma mvContent_mul ...
```
