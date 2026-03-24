# 多变量因式分解 Spec 扩展

> 状态：nl-proof v1
> 目标：扩展 Spec.lean 以支持多变量因式分解（Wang 算法）
> 对应 C++：`polynomial_factorize_wang.hh` 全部

---

## 0. 设计决策

### 0.1 多变量多项式类型

Mathlib 提供 `MvPolynomial σ R`（σ 上的多变量多项式环）。CLPoly 用迭代表示：
`Polynomial (Polynomial (... ℤ ...))` ≃ `MvPolynomial (Fin n) ℤ`。

**选择**：Spec 用 `MvPolynomial (Fin n) ℤ`（数学自然），L2 算法模型用迭代表示
`Polynomial (MvPolynomial (Fin n) ℤ)`（匹配递归结构）。

桥接：Mathlib 的 `MvPolynomial.finSuccEquiv` 提供同构：
```
MvPolynomial (Fin (n+1)) R ≃ₐ Polynomial (MvPolynomial (Fin n) R)
```

### 0.2 顶层 Spec 结构

Wang 算法的子过程需要以下规约（自顶向下）：

```
MvFactorCorrect           -- 多变量完整因式分解（顶层）
├── MvSqfreeDecomp        -- 多变量无平方分解
├── WangCorrect           -- Wang 核心：本原无平方 → 不可约因子
│   ├── EvalPointGood     -- 求值点满足条件 (a)-(d)
│   ├── LCDistribCorrect  -- LC 分配正确性
│   ├── MvHenselCorrect   -- 多变量 Hensel 提升
│   └── TrialDivCorrect   -- 试除验证
└── FactorZZCorrect       -- 单变量因式分解（已有 ✅）
```

---

## 1. 顶层规约：MvFactorCorrect

### 1.1 定义

```lean
/-- 多变量 Z[x₁,...,xₙ] 完整因式分解规约。
    与 FactorZZCorrect 结构相同，但作用于 MvPolynomial。 -/
def MvFactorCorrect {n : ℕ}
    (f : MvPolynomial (Fin n) ℤ)
    (result : List (MvPolynomial (Fin n) ℤ)) : Prop :=
  -- 1. 乘积还原
  Associated f result.prod
  -- 2. 每个因子不可约
  ∧ (∀ g ∈ result, Irreducible g)
```

### 1.2 说明

- 与 `FactorZZCorrect` 定义相同，只是类型从 `Polynomial ℤ` 变为 `MvPolynomial (Fin n) ℤ`
- `MvPolynomial (Fin n) ℤ` 是 UFD（Mathlib：`MvPolynomial.instUniqueFactorizationMonoid`），所以不可约分解存在
- n = 1 时通过 `MvPolynomial.finSuccEquiv` 退化为单变量情况

---

## 2. Wang 核心规约：WangCorrect

### 2.1 定义

```lean
/-- Wang 算法核心规约：本原无平方多变量多项式 → 不可约因子列表。
    对应 C++ __wang_core。-/
def WangCorrect {n : ℕ}
    (f : MvPolynomial (Fin (n + 1)) ℤ)
    (result : List (MvPolynomial (Fin (n + 1)) ℤ)) : Prop :=
  -- 1. 乘积还原
  Associated f result.prod
  -- 2. 每个因子不可约
  ∧ (∀ g ∈ result, Irreducible g)
```

### 2.2 说明

- 数学上与 MvFactorCorrect 相同，但 Wang 的前置条件更强（本原、无平方）
- 前置条件由外层 `__factor_multivar` 保证（先做 SQF 分解）
- `n + 1` 保证至少有一个变量（Wang 需要选主变量）

---

## 3. 求值点规约：EvalPointGood

### 3.1 定义

```lean
/-- 求值点 α = (α₂,...,αₙ₊₁) 对 f ∈ Z[x₁,...,xₙ₊₁] 满足 Wang 条件。
    对应 C++ __select_eval_point 的条件 (a)(b)(b')(d)。-/
def EvalPointGood {n : ℕ}
    (f : MvPolynomial (Fin (n + 2)) ℤ)
    (α : Fin (n + 1) → ℤ)
    (r : ℕ)  -- 期望的单变量因子数
    : Prop :=
  let f₀ := MvPolynomial.eval α (通过 finSuccEquiv 将 f 视为 x₁ 上的多项式，
             在 x₂=α₁,...,xₙ₊₂=αₙ₊₁ 处求值得到 f₀ ∈ Z[x₁])
  -- (a) f₀ ≠ 0
  f₀ ≠ 0
  -- (b) f₀ squarefree
  ∧ Squarefree f₀
  -- (c) lc(f, x₁) 在 α 处不为零
  ∧ MvPolynomial.eval α (leadingCoeff_x1 f) ≠ 0
  -- (d) f₀ 恰好有 r 个不可约因子
  ∧ (不可约因子数 f₀ = r)
```

### 3.2 说明

- 条件 (a)(b) 保证求值后的单变量多项式可因式分解
- 条件 (c) 保证 LC 不退化（LC 分配需要）
- 条件 (d) 保证因子数一致（不会"合并"或"分裂"）
- 条件 (d) 最关键也最难形式化——它要求 f 模 (x₂-α₂,...,xₙ₊₂-αₙ₊₂) 的不可约因子数等于 f 的不可约因子数

### 3.3 工程问题

- `leadingCoeff_x1` 需要从 `MvPolynomial (Fin (n+2)) ℤ` 提取关于 `x₁` 的首项系数，结果是 `MvPolynomial (Fin (n+1)) ℤ`
- 这需要 `finSuccEquiv` 将 f 转为 `Polynomial (MvPolynomial (Fin (n+1)) ℤ)`，然后取 `leadingCoeff`
- "不可约因子数" 需要 UFD 中的唯一分解定理 + 计数

---

## 4. LC 分配规约：LCDistribCorrect

### 4.1 定义

```lean
/-- LC 分配规约：将 lc(f, x₁) 分配到各因子的 LC 目标。
    对应 C++ __wang_leading_coeff。-/
def LCDistribCorrect {n : ℕ}
    (L : MvPolynomial (Fin (n + 1)) ℤ)   -- lc(f, x₁)
    (factors_univar : List (Polynomial ℤ)) -- 单变量因子 u₁,...,uᵣ
    (σ : List (MvPolynomial (Fin (n + 1)) ℤ))  -- LC 目标
    (α : Fin (n + 1) → ℤ)                 -- 求值点
    : Prop :=
  -- 1. 数量一致
  factors_univar.length = σ.length
  -- 2. 乘积条件：∏σᵢ ∼ L（up to content）
  ∧ Associated L σ.prod
  -- 3. 求值一致：σᵢ(α) = lc(uᵢ)（每个 LC 目标在求值点处等于对应单变量因子的 LC）
  ∧ List.Forall₂ (fun u s => MvPolynomial.eval α s = Polynomial.leadingCoeff u) factors_univar σ
```

### 4.2 说明

- 这是 Wang 算法最关键也最容易出错的步骤（5/11 bug 出在这里）
- C++ 实现经历了多次重写：power extraction → GCD matching → non-divisor algorithm
- 核心数学：Gauss 引理 + 多变量 UFD 中的 LC 因子分解

---

## 5. 多变量 Hensel 提升规约：MvHenselCorrect

### 5.1 定义

```lean
/-- 多变量 Hensel 提升规约。
    给定 f ∈ Z[x₁,...,xₙ] 和求值点 α，以及 f₀ = f(x₁,α₂,...,αₙ) 的因子 u₁,...,uᵣ，
    提升到 Z[x₁,...,xₙ] 的因子 G₁,...,Gᵣ，满足在 α 处求值还原到 uᵢ。
    对应 C++ __mtshl_lift。-/
def MvHenselCorrect {n : ℕ}
    (f : MvPolynomial (Fin (n + 1)) ℤ)
    (α : Fin n → ℤ)
    (factors_univar : List (Polynomial ℤ))
    (factors_mv : List (MvPolynomial (Fin (n + 1)) ℤ))
    : Prop :=
  -- 1. 乘积还原
  Associated f factors_mv.prod
  -- 2. 因子数一致
  ∧ factors_univar.length = factors_mv.length
  -- 3. 求值一致：Gᵢ(x₁,α₂,...) = uᵢ(x₁)
  ∧ List.Forall₂ (fun u G => eval_at_α G α = u) factors_univar factors_mv
```

### 5.2 说明

- 这里 `eval_at_α` 是在 x₂=α₁,...,xₙ₊₁=αₙ 处求值，保留 x₁
- C++ 用 MTSHL（稀疏 Hensel）或 dense Hensel（退化情况）
- 不变量：每步提升后 ∏Gᵢ ≡ f (mod (x₂-α₁)^k₂ · ... · (xₙ₊₁-αₙ)^kₙ₊₁)

---

## 6. Pipeline 框架：MvFactorZZ

### 6.1 定义

```lean
/-- 多变量 Z[x₁,...,xₙ] 因式分解的顶层正确性。
    组合 SQF + Wang，类似 factor_ZZ_correct。-/
theorem mv_factor_correct {n : ℕ}
    (f : MvPolynomial (Fin n) ℤ) (hf : f ≠ 0)
    -- 假设子过程正确
    (sqf : ...) (hsqf : ...)
    (wang : ...) (hwang : ...)
    (factor_univar : ...) (hfactor_univar : ...)
    : ∃ result, MvFactorCorrect f result
```

### 6.2 说明

- 结构与 `factor_ZZ_correct` 平行
- SQF 分解得到无平方分量 → 对每个分量调用 Wang → 合并结果
- Wang 内部递归调用单变量因式分解（已验证 ✅）

---

## 7. Lean 形式化注意事项

### 7.1 MvPolynomial 在 Mathlib 中的状态

需要确认的 Mathlib API：
- `MvPolynomial.finSuccEquiv`：`MvPolynomial (Fin (n+1)) R ≃ₐ Polynomial (MvPolynomial (Fin n) R)`
- `MvPolynomial.eval`：`(σ → R) → MvPolynomial σ R → R`
- `MvPolynomial` 是 UFD：`instUniqueFactorizationMonoid`（需要 R 是 UFD + σ 有限）
- `Squarefree` 在 `MvPolynomial` 上的行为

### 7.2 与单变量的桥接

通过 `finSuccEquiv`，`MvPolynomial (Fin 1) ℤ ≃ₐ Polynomial ℤ`。
所以 `MvFactorCorrect` 在 n=1 时等价于 `FactorZZCorrect`。
单变量验证成果直接可用。

### 7.3 递归结构

Wang 算法递归减少变量数：
- `__factor_multivar` 对 n 变量调用 `__wang_core`
- `__wang_core` 对 (n-1) 变量 LC 递归调用 `factorize`（变量数减 1）
- 终止条件：变量数 = 1 → 单变量因式分解（已验证 ✅）

这需要对 n 做归纳的 Pipeline 证明。

---

## 8. 形式化估计

| 规约 | 行数 |
|------|------|
| `MvFactorCorrect` | ~10 |
| `WangCorrect` | ~10 |
| `EvalPointGood` | ~20 |
| `LCDistribCorrect` | ~15 |
| `MvHenselCorrect` | ~15 |
| Pipeline 框架 `mv_factor_correct` | ~30 |
| **总计** | **~100** |

注：这只是 Spec 扩展。L2 算法模型（证明 Wang 满足这些 Spec）是后续的大工程（~1150 行）。
