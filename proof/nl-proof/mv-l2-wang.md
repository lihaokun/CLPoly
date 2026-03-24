# 多变量因式分解 L2：Wang 算法模型

> 状态：nl-proof v1
> 对应 C++：`polynomial_factorize_wang.hh` `__wang_core`

---

## 0. Wang 核心逻辑（抽象为数学步骤）

```
输入：f ∈ Z[x₁,...,xₙ]，primitive squarefree
1. 选求值点 α = (α₂,...,αₙ)
2. 求值：f₀ = f(x₁, α₂,...,αₙ) ∈ Z[x₁]
3. 单变量因式分解：f₀ = u₁ · ... · uᵣ（已验证 ✅）
4. LC 分配：L = lc(f, x₁)，找 σ₁,...,σᵣ 使 ∏σᵢ ∼ L
5. Hensel 提升：从 uᵢ 恢复 Gᵢ ∈ Z[x₁,...,xₙ]
6. 试除验证：∏(subset Gᵢ) | f → 提取真因子
输出：f 的不可约因子
```

## 1. 按模块分解 L2 证明

### 1.1 模块依赖图

```
wang_correct
├── eval_factorization       -- 求值 + 单变量分解
│   └── factor_ZZ_correct ✅  -- 已验证
├── lc_distrib_correct       -- LC 分配
│   └── UFD + valuation      -- 多变量 UFD
├── mv_hensel_correct        -- 多变量 Hensel 提升
│   └── hensel_step_with_degree ✅  -- 已验证（单变量）
└── trial_div_correct        -- 试除验证
    └── 整除 → 因子（trivial）
```

### 1.2 优先级

| 模块 | 难度 | 行数 | 依赖 | 优先级 |
|------|------|------|------|--------|
| eval_factorization | 低 | ~50 | 单变量已验证 | **1** |
| trial_div_correct | 低 | ~30 | trivial | **2** |
| lc_distrib_correct | 高 | ~200 | UFD + valuation | 3 |
| mv_hensel_correct | 高 | ~300 | Taylor 展开 + Diophantine | 4 |

---

## 2. 模块 1：eval_factorization

### 2.1 定理

```lean
/-- 求值后单变量因式分解：若 f₀ = f(x₁, α) 有不可约分解，
    则 f 有对应的"因子在 α 处求值"结构。
    这是 Wang 算法的起点——连接多变量和单变量。-/
theorem eval_gives_factorization {n : ℕ}
    (f : MvPolynomial (Fin (n + 1)) ℤ) (hf : f ≠ 0)
    (α : Fin n → ℤ)
    -- f₀ = f(x₁, α) 的单变量因式分解
    (univar_result : List (Polynomial ℤ))
    (h_correct : FactorZZCorrect (eval_at_α f α) univar_result)
    : Associated (eval_at_α f α) univar_result.prod
      ∧ ∀ g ∈ univar_result, Irreducible g
```

### 2.2 说明

这个定理只是拆包 `FactorZZCorrect`。真正的价值在于定义 `eval_at_α`：

```lean
/-- 在 α = (α₂,...,αₙ₊₁) 处求值，保留 x₁（第一个变量）。
    结果：MvPolynomial (Fin (n+1)) ℤ → Polynomial ℤ -/
noncomputable def eval_at_α {n : ℕ}
    (f : MvPolynomial (Fin (n + 1)) ℤ) (α : Fin n → ℤ) : Polynomial ℤ :=
  Polynomial.map (MvPolynomial.eval α)
    (MvPolynomial.finSuccEquiv ℤ n f)
```

即：finSuccEquiv 转为 Polynomial (MvPolynomial (Fin n) ℤ)，然后 map(eval α) 得到 Polynomial ℤ。

**Mathlib 已有**：`MvPolynomial.eval_eq_eval_mv_eval'` 证明了这正好等于 `eval (Fin.cons y α) f` 在 y 处求值后得到的。

### 2.3 Lean API

- `MvPolynomial.finSuccEquiv ℤ n`：`MvPolynomial (Fin (n+1)) ℤ ≃ₐ Polynomial (MvPolynomial (Fin n) ℤ)`
- `Polynomial.map (MvPolynomial.eval α)`：`Polynomial (MvPolynomial (Fin n) ℤ) → Polynomial ℤ`
- `MvPolynomial.eval_eq_eval_mv_eval'`：连接 eval 和 finSuccEquiv

---

## 3. 模块 2：trial_div_correct

### 3.1 定理

```lean
/-- 试除验证：若 g | f 且 g 不可约，则 g 是 f 的不可约因子。-/
theorem trial_div_factor {R : Type*} [CommRing R] [IsDomain R]
    (f g : R) (hg_irred : Irreducible g) (hg_dvd : g ∣ f) :
    Irreducible g ∧ g ∣ f
```

这是 trivial 的。试除验证的 L2 模型就是：对 Hensel 提升结果做 divides 检查。

---

## 4. 模块 3：lc_distrib_correct（最难）

### 4.1 问题

L = lc(f, x₁) ∈ Z[x₂,...,xₙ]。
f₀ = f(x₁, α) 的因子 u₁,...,uᵣ。
需要 σ₁,...,σᵣ 使得：
- ∏ σᵢ ∼ L
- σᵢ(α) = lc(uᵢ) × (content adjustment)

### 4.2 数学核心

**定理**（LC 分配存在性）：若 f 的不可约因子为 G₁,...,Gₛ in Z[x₁,...,xₙ]，则
lc(f, x₁) = (up to unit) ∏ lc(Gⱼ, x₁)。

这是 Gauss 引理的多变量推广：primitive × primitive = primitive，
所以 lc(f) = lc(G₁ · ... · Gₛ) = ∏ lc(Gⱼ)（up to unit，因 content 可能变化）。

### 4.3 形式化路径

不直接证 C++ 的 valuation 提取算法正确（太复杂，5/11 bug 证明它容易出错）。
而是证：**存在**满足条件的 σ₁,...,σᵣ。具体算法如何找到它们是 L1 的事。

```lean
theorem lc_distrib_exists {n : ℕ}
    (f : MvPolynomial (Fin (n + 1)) ℤ)
    (hf_ne : f ≠ 0) (hf_sqfree : Squarefree f)
    -- f 的不可约因子 G₁,...,Gₛ
    (factors : List (MvPolynomial (Fin (n + 1)) ℤ))
    (hassoc : Associated f factors.prod)
    (hirred : ∀ g ∈ factors, Irreducible g)
    : -- lc(f, x₁) 可以分配为 ∏ lc(Gⱼ, x₁)
      Associated (lc_x1 f)
        (factors.map (fun g => lc_x1 g)).prod
```

其中 `lc_x1 f = (finSuccEquiv ℤ n f).leadingCoeff`。

### 4.4 证明思路

`f ∼ ∏ Gⱼ`
→ `finSuccEquiv f ∼ ∏ finSuccEquiv Gⱼ`（代数同构保持 Associated）
→ `lc(finSuccEquiv f) ∼ ∏ lc(finSuccEquiv Gⱼ)`（多项式乘积的 lc = lc 的乘积，当无零因子时）

关键引理：`Polynomial.leadingCoeff_mul` 在整环上成立。

---

## 5. 模块 4：mv_hensel_correct（最大工程量）

### 5.1 概述

多变量 Hensel 提升：从 f₀ = f(x₁, α) 的因子 u₁,...,uᵣ，
逐变量恢复 G₁,...,Gᵣ ∈ Z[x₁,...,xₙ]，使 Gᵢ(x₁, α) = uᵢ 且 f ∼ ∏ Gᵢ。

### 5.2 数学核心

这是单变量 Hensel 提升的多变量推广。不是 mod p^k，而是 mod (x₂-α₂)^k₂ · ... · (xₙ-αₙ)^kₙ。

Taylor 展开：f(x₁, x₂,...) = f(x₁, α₂,...) + ∑ (∂f/∂xⱼ)(α) · (xⱼ - αⱼ) + ...

逐变量提升：先固定 x₃,...,xₙ = α，只提升 x₂，得到 f mod (x₃-α₃)·...·(xₙ-αₙ)。然后提升 x₃，等等。

### 5.3 务实方案

完整形式化 MTSHL 太复杂（稀疏插值、Vandermonde 求解、MDP cascade）。

**方案**：证明提升**存在性**——给定好的求值点，多变量因子可以从单变量因子恢复。

```lean
theorem mv_hensel_exists {n : ℕ}
    (f : MvPolynomial (Fin (n + 1)) ℤ)
    (α : Fin n → ℤ)
    (univar_facs : List (Polynomial ℤ))
    (h_eval : Associated (eval_at_α f α) univar_facs.prod)
    (h_coprime : List.Pairwise (fun a b => IsCoprime a b) univar_facs)
    -- 精度足够（Mignotte 多变量推广）
    : ∃ lifted : List (MvPolynomial (Fin (n + 1)) ℤ),
        Associated f lifted.prod
        ∧ List.Forall₂ (fun u G => eval_at_α G α = u) univar_facs lifted
```

实际上这也是 Hensel 唯一性 + 存在性的多变量版本，可以从单变量 Hensel 结果推广。

---

## 6. 整体组合

```lean
theorem wang_correct {n : ℕ}
    (f : MvPolynomial (Fin (n + 1)) ℤ)
    (hf_ne : f ≠ 0) (hf_sqfree : Squarefree f)
    -- 好的求值点存在（算法前提）
    (α : Fin n → ℤ)
    -- 单变量分解可用（已验证 ✅）
    : ∃ result, MvFactorCorrect f result
```

这最终可以由 UFD 直接给出（和 mv_factor_instantiate 一样）。Wang L2 的价值在于证明**每个子过程**的数学正确性，不在于顶层定理。

---

## 7. 形式化计划

| 模块 | 行数 | 状态 |
|------|------|------|
| `eval_at_α` 定义 + 基本性质 | ~30 | 先做 |
| `lc_x1` 定义 + lc 乘积引理 | ~40 | 先做 |
| `lc_distrib_exists` | ~50 | 中期 |
| `mv_hensel_exists` | ~80 | 后期 |
| `trial_div_factor` | ~10 | trivial |
| **总计** | **~210** |
