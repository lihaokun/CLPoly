# 多变量 L2：MTSHL Taylor 提升循环

> 状态：nl-proof v1
> 对应 C++：`__mtshl_step_j`（lines 763-935）

---

## 0. C++ 核心循环（剥离工程细节后）

```
输入：
  aj ∈ Zp[x₁,...,xⱼ]     -- 目标多项式（f 在 xⱼ₊₁=αⱼ₊₁,... 处求值）
  F[0..r-1] ∈ Zp[x₁,...,xⱼ₋₁]  -- 当前因子近似
  αⱼ ∈ Zp               -- 求值点
不变量：aj ≡ ∏F[i] (mod (xⱼ-αⱼ)^k)

循环 k = 1 to Dⱼ = deg(aj, xⱼ):
  1. error = aj - ∏F[i]
  2. cₖ = Taylor系数(error, xⱼ, αⱼ, k)     -- error 的 (xⱼ-αⱼ)^k 系数
  3. MDP求解: cₖ = Σᵢ σᵢₖ · F̂ᵢ (mod xⱼ=αⱼ)
     其中 F̂ᵢ = ∏_{j≠i} F[j]|_{xⱼ=αⱼ}
  4. F[i] += σᵢₖ · (xⱼ-αⱼ)^k
  5. 重算 error；若 error = 0 则终止

输出：F[0..r-1] ∈ Zp[x₁,...,xⱼ]，满足 aj = ∏F[i]
```

## 1. 数学核心：Taylor 展开 + Diophantine

### 1.1 不变量

**初始**：F[i] 只含 x₁,...,xⱼ₋₁。∏F[i] = aj|_{xⱼ=αⱼ}。
所以 error = aj - ∏F[i] ∈ (xⱼ - αⱼ) · Zp[x₁,...,xⱼ]。
即 error ≡ 0 (mod xⱼ - αⱼ)。

**第 k 步后**：error ≡ 0 (mod (xⱼ - αⱼ)^(k+1))。
即 ∏F[i] ≡ aj (mod (xⱼ - αⱼ)^(k+1))。

**终止**：k = Dⱼ 后，若 deg(error, xⱼ) < Dⱼ+1 且 error ≡ 0 mod (xⱼ-αⱼ)^(Dⱼ+1)，
则 error = 0（因 deg < Dⱼ+1 但含 (xⱼ-αⱼ)^(Dⱼ+1) 因子 → 必为零）。

### 1.2 MDP（多变量 Diophantine 问题）

给 cₖ ∈ Zp[x₁,...,xⱼ₋₁]，找 σᵢₖ 使得：
cₖ = Σᵢ σᵢₖ · (∏_{j≠i} F[j]|_{xⱼ=αⱼ})

这是偏分式分解：cₖ / ∏F = Σ σᵢ / Fᵢ。

**j=2 时**：F[i] 是单变量多项式，MDP 就是单变量偏分式分解（Bézout 系数）。
**j≥3 时**：多变量 MDP，C++ 用 sparse interpolation 或 WMDS 回退。

### 1.3 为什么不变量成立

第 k 步：
- error_old ≡ 0 (mod (xⱼ-αⱼ)^k)
- cₖ = error 的 (xⱼ-αⱼ)^k 系数
- MDP: cₖ = Σ σᵢₖ · F̂ᵢ
- F[i]_new = F[i]_old + σᵢₖ · (xⱼ-αⱼ)^k
- ∏F_new = ∏F_old + Σ σᵢₖ · F̂ᵢ · (xⱼ-αⱼ)^k + O((xⱼ-αⱼ)^(k+1))
         = ∏F_old + cₖ · (xⱼ-αⱼ)^k + O((xⱼ-αⱼ)^(k+1))
- error_new = aj - ∏F_new
            = error_old - cₖ · (xⱼ-αⱼ)^k + O((xⱼ-αⱼ)^(k+1))
            ≡ 0 (mod (xⱼ-αⱼ)^(k+1)) ✓

## 2. L2 模型

### 2.1 Taylor 系数

```lean
/-- Taylor 系数：提取 f 在 (xⱼ - αⱼ)^k 处的系数。
    对应 C++ __taylor_coeff_zp。-/
noncomputable def taylorCoeff {n : ℕ}
    (f : MvPolynomial (Fin (n + 1)) ℤ)
    (j : Fin n) (αⱼ : ℤ) (k : ℕ) : MvPolynomial (Fin (n + 1)) ℤ :=
  sorry -- 数学定义：f = Σₖ cₖ · (xⱼ-αⱼ)^k，提取 cₖ
```

### 2.2 MDP 求解规约

```lean
/-- MDP（多变量 Diophantine）求解规约。
    给定 cₖ 和基础因子 F_base[i]，找 σᵢ 使得 cₖ = Σ σᵢ · F̂ᵢ。
    对应 C++ __mtshl_zp_univar_mdp / __mtshl_sparse_int / __mtshl_wmds。-/
structure MDPSolution {n : ℕ}
    (ck : MvPolynomial (Fin (n + 1)) ℤ)
    (f_base : List (MvPolynomial (Fin (n + 1)) ℤ))
    (sigma : List (MvPolynomial (Fin (n + 1)) ℤ)) : Prop where
  -- cₖ = Σ σᵢ · F̂ᵢ，F̂ᵢ = ∏_{j≠i} f_base[j]
  sum_eq : ck = (sigma.zipWith (fun si fi_hat => si * fi_hat)
    (f_base.map (fun fi => (f_base.erase fi).prod))).sum
  length_eq : sigma.length = f_base.length
```

### 2.3 单步提升

```lean
/-- MTSHL 单步：给定误差的 Taylor 系数 + MDP 解，更新因子。
    不变量：误差 mod (xⱼ-αⱼ)^(k+1) = 0。-/
def mtshlStep {n : ℕ}
    (factors : List (MvPolynomial (Fin (n + 1)) ℤ))
    (sigma : List (MvPolynomial (Fin (n + 1)) ℤ))
    (xj_minus_alpha_pow_k : MvPolynomial (Fin (n + 1)) ℤ)
    : List (MvPolynomial (Fin (n + 1)) ℤ) :=
  factors.zipWith (fun fi si => fi + si * xj_minus_alpha_pow_k) sigma
```

### 2.4 完整循环

```lean
/-- MTSHL 逐变量提升循环。
    对应 C++ __mtshl_step_j 的 for k=1..Dj 循环。-/
noncomputable def mtshlLoop {n : ℕ}
    (target : MvPolynomial (Fin (n + 1)) ℤ)  -- aj
    (factors : List (MvPolynomial (Fin (n + 1)) ℤ))
    (j : Fin n) (αⱼ : ℤ) (Dj : ℕ)
    -- MDP 求解器（外部提供，L2 不建模具体实现）
    (mdp_solver : MvPolynomial (Fin (n + 1)) ℤ →
        List (MvPolynomial (Fin (n + 1)) ℤ) →
        List (MvPolynomial (Fin (n + 1)) ℤ))
    : List (MvPolynomial (Fin (n + 1)) ℤ) :=
  -- 递归实现 k = 1 to Dj
  sorry -- 循环体
```

### 2.5 正确性定理

```lean
/-- MTSHL 循环正确性：若 MDP 求解器正确，则循环终止后 ∏F[i] = target。
    关键不变量：每步后 error ≡ 0 (mod (xⱼ-αⱼ)^(k+1))。-/
theorem mtshl_correct {n : ℕ}
    (target : MvPolynomial (Fin (n + 1)) ℤ)
    (factors_init : List (MvPolynomial (Fin (n + 1)) ℤ))
    (j : Fin n) (αⱼ : ℤ) (Dj : ℕ)
    -- 初始条件：∏factors_init = target|_{xⱼ=αⱼ}
    (h_init : target |_{xⱼ=αⱼ} = factors_init.prod)
    -- MDP 正确
    (mdp_solver : ...) (h_mdp : ...)
    -- deg(target, xⱼ) ≤ Dj
    (h_deg : MvPolynomial.degreeOf (Fin.succ j) target ≤ Dj)
    : let result := mtshlLoop target factors_init j αⱼ Dj mdp_solver
      target = result.prod
```

## 3. 简化方案

完整建模 mtshlLoop + 不变量归纳需要处理：
- Taylor 系数提取（涉及多变量多项式除法）
- MDP 求解正确性（偏分式分解）
- 理想幂次 mod (xⱼ-αⱼ)^k 的形式化

这些每个都需要大量 Mathlib 基础设施。

**务实方案**：建模循环骨架 + MDP 规约，不变量证明用 sorry 标注，后续逐步填充。循环的**终止条件**（error = 0 when k > Dj）是核心数学事实，可以先证。

## 4. 形式化估计

| 内容 | 行数 | sorry |
|------|------|-------|
| `taylorCoeff` 定义 | ~10 | 1（定义用 sorry） |
| `MDPSolution` 规约 | ~10 | 0 |
| `mtshlStep` 单步 | ~5 | 0 |
| `mtshlLoop` 循环 | ~30 | 1（循环体） |
| `mtshl_correct` 正确性 | ~30 | 1（不变量归纳） |
| **总计** | **~85** | **3** |
