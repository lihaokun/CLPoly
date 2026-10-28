# 多变量 L2 模块 1：求值 + 单变量因式分解

> 状态：nl-proof v1
> 对应 C++：`__wang_core` lines 1-6（求值 + factorize(f₀)）
> 依赖：`factor_ZZ_correct`（已验证 ✅）

---

## 0. C++ 逻辑

```cpp
// __wang_core, 简化
f0 = eval(g, x1, alpha);           // f₀ = g(x₁, α₂,...,αₙ)
auto [facs, lc_f0] = factorize(f0); // 单变量因式分解
r = facs.size();
if (r <= 1) { mark_dead(x1); break; } // 不可约方向
```

## 1. Lean 模型

### 1.1 定义 eval_at_α

```lean
/-- 在 (α₂,...,αₙ₊₁) 处求值，保留 x₁。
    对应 C++ eval(g, x1, alpha)。
    实现：finSuccEquiv 转为 Poly(MvPoly)，然后 map(eval α)。-/
noncomputable def eval_at_α {n : ℕ}
    (f : MvPolynomial (Fin (n + 1)) ℤ) (α : Fin n → ℤ) : Polynomial ℤ :=
  Polynomial.map (MvPolynomial.eval α ∘ (Int.castRingHom _))
    ... -- 需要仔细处理类型
```

**实际上**：`finSuccEquiv ℤ n f` 的类型是 `Polynomial (MvPolynomial (Fin n) ℤ)`。
对其系数 apply `MvPolynomial.eval α` 得到 `Polynomial ℤ`。

```lean
noncomputable def eval_at_α {n : ℕ}
    (f : MvPolynomial (Fin (n + 1)) ℤ) (α : Fin n → ℤ) : Polynomial ℤ :=
  Polynomial.map (MvPolynomial.eval α) (MvPolynomial.finSuccEquiv ℤ n f)
```

但 `MvPolynomial.eval α : MvPolynomial (Fin n) ℤ →+* ℤ`，
而 `Polynomial.map` 需要 `R →+* S`，这里 `R = MvPolynomial (Fin n) ℤ`，`S = ℤ`。
所以 `Polynomial.map (MvPolynomial.eval α)` 把 `Polynomial (MvPolynomial (Fin n) ℤ)` 映射到 `Polynomial ℤ`。✓

### 1.2 eval_at_α 与 MvPolynomial.eval 的关系

Mathlib 已有 `eval_eq_eval_mv_eval'`：
```
MvPolynomial.eval (Fin.cons y α) f =
  Polynomial.eval y (Polynomial.map (MvPolynomial.eval α) (finSuccEquiv ℤ n f))
```

即 `MvPolynomial.eval (Fin.cons y α) f = Polynomial.eval y (eval_at_α f α)`。

这说明：`eval_at_α f α` 是"对 f 在 x₂=α₁,...,xₙ₊₁=αₙ 求值后得到的关于 x₁ 的单变量多项式"。✓

### 1.3 eval_at_α 保持乘积

```lean
lemma eval_at_α_mul {n : ℕ} (f g : MvPolynomial (Fin (n + 1)) ℤ) (α : Fin n → ℤ) :
    eval_at_α (f * g) α = eval_at_α f α * eval_at_α g α := by
  simp [eval_at_α, map_mul]
```

`finSuccEquiv` 是代数同构（保持乘法），`Polynomial.map` 是环同态（保持乘法）。两者组合保持乘法。✓

### 1.4 eval_at_α 保持 Associated

```lean
lemma eval_at_α_associated {n : ℕ} {f g : MvPolynomial (Fin (n + 1)) ℤ}
    (h : Associated f g) (α : Fin n → ℤ) :
    Associated (eval_at_α f α) (eval_at_α g α) := by
  exact h.map (Polynomial.mapRingHom (MvPolynomial.eval α) |>.comp
    (finSuccEquiv ℤ n).toRingHom)
```

---

## 2. 求值点条件

### 2.1 C++ `__select_eval_point` 的条件

```cpp
// 条件 (a): f₀ ≠ 0
// 条件 (b): f₀ squarefree
// 条件 (c): lc(f, x₁)(α) ≠ 0
// 条件 (d): f₀ 的不可约因子数 = f 的不可约因子数
```

### 2.2 Lean 规约

```lean
/-- 求值点 α 对 f 是"好的"：求值后保持必要性质。
    对应 C++ __select_eval_point 的条件检查。-/
structure EvalPointGood {n : ℕ}
    (f : MvPolynomial (Fin (n + 1)) ℤ) (α : Fin n → ℤ) : Prop where
  -- (a) 求值后非零
  eval_ne_zero : eval_at_α f α ≠ 0
  -- (b) 求值后 squarefree
  eval_sqfree : Squarefree (eval_at_α f α)
  -- (c) lc(f, x₁) 在 α 处不为零
  lc_ne_zero : MvPolynomial.eval α ((finSuccEquiv ℤ n f).leadingCoeff) ≠ 0
```

注：条件 (d)（因子数保持）不在 EvalPointGood 中——它是算法正确性的深层条件，在 L2 中作为 Wang 核心的前提假设。

### 2.3 求值点条件的满足

条件 (a)(b)(c) 对"几乎所有" α 成立（有限域上的非零条件）。
C++ 通过枚举尝试找到好的 α。L2 中我们取 α 为参数，假设 EvalPointGood 成立。

---

## 3. Wang 单步（eval + factor + 检查）

### 3.1 定理

```lean
/-- Wang 第一步：求值 + 单变量因式分解。
    对应 C++ __wang_core 的 eval + factorize(f₀) 部分。

    给定好的求值点，单变量因式分解给出 f₀ 的不可约因子。
    这些因子是后续 LC 分配和 Hensel 提升的输入。-/
theorem wang_eval_step {n : ℕ}
    (f : MvPolynomial (Fin (n + 1)) ℤ) (hf : f ≠ 0)
    (α : Fin n → ℤ) (hα : EvalPointGood f α)
    : ∃ univar_facs : List (Polynomial ℤ),
        FactorZZCorrect (eval_at_α f α) univar_facs := by
  exact factor_ZZ_instantiate (eval_at_α f α) hα.eval_ne_zero
```

### 3.2 说明

这个定理本身很简单——直接调用已验证的单变量因式分解。
它的价值在于：
1. 定义了 `eval_at_α`（桥接多变量和单变量）
2. 定义了 `EvalPointGood`（求值点条件规约）
3. 后续模块（LC 分配、Hensel）以此为起点

---

## 4. lc_x1 定义

### 4.1 定义

```lean
/-- f 关于 x₁ 的首项系数。
    对应 C++ lc(g, x₁)。-/
noncomputable def lc_x1 {n : ℕ}
    (f : MvPolynomial (Fin (n + 1)) ℤ) : MvPolynomial (Fin n) ℤ :=
  (MvPolynomial.finSuccEquiv ℤ n f).leadingCoeff
```

### 4.2 性质

```lean
-- lc(f*g, x₁) = lc(f, x₁) * lc(g, x₁)（整环上）
lemma lc_x1_mul (f g : MvPolynomial (Fin (n + 1)) ℤ)
    (hf : f ≠ 0) (hg : g ≠ 0) :
    lc_x1 (f * g) = lc_x1 f * lc_x1 g

-- lc(f, x₁)(α) = lc(f(x₁, α), x₁) = lc(eval_at_α f α)
lemma eval_lc_x1 (f : MvPolynomial (Fin (n + 1)) ℤ) (α : Fin n → ℤ)
    (h : MvPolynomial.eval α (lc_x1 f) ≠ 0) :
    (eval_at_α f α).leadingCoeff = MvPolynomial.eval α (lc_x1 f)
```

第一个：`Polynomial.leadingCoeff_mul`（需要整环 + 非零）。
第二个：`Polynomial.leadingCoeff_map_of_leadingCoeff_ne_zero`（map 不降度时 lc 保持）。

---

## 5. 形式化估计

| 内容 | 行数 |
|------|------|
| `eval_at_α` 定义 + mul/associated 引理 | ~20 |
| `EvalPointGood` 结构体 | ~10 |
| `wang_eval_step` 定理 | ~5 |
| `lc_x1` 定义 + 性质 | ~30 |
| **总计** | **~65** |

这是 Wang L2 的**第一个模块**。完成后接 LC 分配（模块 2）。
