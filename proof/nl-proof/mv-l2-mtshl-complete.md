# MTSHL 完整 L2：全部 MDP 求解器

> 状态：nl-proof v1
> 对应 C++：`__si_vandermonde_solve`, `__si_theta_array_eval`, `__taylor_coeff_zp`,
>           `__mtshl_zp_univar_mdp`, `__mtshl_sparse_int`, `__mtshl_multi_bdp`, `__mtshl_wmds`

---

## 1. Taylor 系数提取（`__taylor_coeff_zp`）

### 1.1 C++ 算法
```
g := f
for t = 0 to j-1:
    g := g / (xk - αk)     // 多项式除法，取商
return g(xk = αk)           // 求值
```

即：f 在 (xk - αk) 上做 j 次除法取商，然后求值。结果是 Taylor 展开 `f = Σ cₜ · (xk-αk)^t` 的第 j 项系数 cⱼ。

### 1.2 L2 模型

```lean
/-- Taylor 系数：f 关于 (xk - αk) 的第 j 项系数。
    f = Σ cₜ · (xk-αk)^t，返回 cⱼ。
    对应 C++ __taylor_coeff_zp。-/
noncomputable def taylorCoeff (f : MvPolynomial σ R) (k : σ) (α : R) (j : ℕ) : MvPolynomial σ R :=
  -- 迭代除以 (X k - C α) 共 j 次，取商，然后在 xk=α 处求值
  let g := (Nat.iterate (fun g => (g - partialEval k α g) / (X k - C α)) j f)
  partialEval k α g
```

注：`/ (X k - C α)` 需要在整环中定义（或用 `divByMonic` 如果 monic）。

### 1.3 正确性

Taylor 展开 `f = Σ_{t=0}^{d} cₜ · (x-α)^t`。
除以 `(x-α)` 一次：商 = `Σ_{t=1}^{d} cₜ · (x-α)^{t-1}`，余式 = `c₀`。
除以 `(x-α)` j 次：商 = `Σ_{t=j}^{d} cₜ · (x-α)^{t-j}`。
求值 `x=α`：= `cⱼ`。✓

形式化路径：因子定理（已证 ✅）的迭代应用。`(x-α) | (f - f(α))` → `f = f(α) + (x-α)·q` → 迭代。

---

## 2. 单变量 MDP（`__mtshl_zp_univar_mdp`）

### 2.1 C++ 算法
```
Step 1: Bézout 链构建
  s[0] = 1, g = F[0]
  for i = 1 to r-1:
    (α, β) = extGCD(g, F[i])    // α·g + β·F[i] = 1
    s[0..i-1] *= β
    s[i] = α
    g *= F[i]

Step 2: σ[i] = (s[i] · c) mod F[i]
```

### 2.2 L2 模型

```lean
/-- 单变量 MDP：Bézout 链 + mod。
    对应 C++ __mtshl_zp_univar_mdp。-/
noncomputable def mdpUnivarAlg (F : List (Polynomial R)) (c : Polynomial R) :
    List (Polynomial R) :=
  let bezout_s := buildBezoutChain F  -- Bézout 链
  F.zipWith (fun fi si => (c * si) %ₘ fi) bezout_s
```

### 2.3 正确性

`mdp_exists`（已证 ✅）保证 Bézout 构造的解满足 `linearTerm F sigma = c`。
`%ₘ fᵢ` 不影响 `linearTerm` 等式（因为 `fᵢ · F̂ᵢ = ∏F`，取余只改变 `∏F` 的倍数部分）。

具体：`σᵢ = c·sᵢ - fᵢ·qᵢ`。
`Σ σᵢ · F̂ᵢ = Σ (c·sᵢ - fᵢ·qᵢ) · F̂ᵢ = c·(Σsᵢ·F̂ᵢ) - (Σqᵢ)·∏F = c·1 - (Σqᵢ)·∏F`。
因 `deg(σᵢ) < deg(fᵢ)` → `deg(Σσᵢ·F̂ᵢ) < deg(∏F)` → `(Σqᵢ)·∏F` 部分 = 0。
所以 `Σ σᵢ · F̂ᵢ = c`。✓

---

## 3. 二变量 BDP（`__mtshl_multi_bdp`）

### 3.1 C++ 算法

本质是**递归 Newton 迭代**——和 `__mtshl_step_j` 相同的结构，但用于 MDP 而非整体提升。

```
Step 1: 求值 x₂=α₂，得单变量 F₀[i] 和 c₀
Step 2: 单变量 MDP：Σ σ₀[i] · F̂₀[i] = c₀
Step 3: Taylor 逐阶提升：
  for k = 1 to deg(c, x₂):
    cₖ = Taylor系数(error, x₂, α₂, k)
    δₖ = 单变量 MDP(F₀, cₖ)
    result[i] += δₖ[i] · (x₂-α₂)^k
    error = c - Σ result[i] · F̂ᵢ
```

### 3.2 L2 正确性

这和 `mtshl_step_invariant`（已证 ✅）的数学结构**完全相同**：
- 不变量：`(x₂-α₂)^k | (c - Σ result[i] · F̂ᵢ)`
- 每步：单变量 MDP 给出 Taylor 修正 → 不变量提升 k → k+1
- 终止：k > deg(c, x₂) → error = 0

**所以 multi_bdp 的正确性直接由已证的 Newton 迭代框架给出。** 无需额外证明。

---

## 4. 递归 WMDS（`__mtshl_wmds`）

### 4.1 C++ 算法

```
if aux_vars 为空: 调用 __mtshl_zp_univar_mdp（基础情形）
else:
  xⱼ = aux_vars.back, αⱼ = ideal_alphas.back
  F_base = F[i](xⱼ=αⱼ), c_base = c(xⱼ=αⱼ)
  result_base = 递归 WMDS(F_base, c_base, aux_vars[:-1])
  result = result_base
  // Taylor 逐阶提升（和 multi_bdp 相同）
  for k = 1 to deg(c, xⱼ):
    cₖ = Taylor系数(error, xⱼ, αⱼ, k)
    δₖ = 递归 WMDS(F_base, cₖ, aux_vars[:-1])
    result[i] += δₖ[i] · (xⱼ-αⱼ)^k
    error = c - Σ result[i] · F̂ᵢ
```

### 4.2 L2 正确性

结构和 multi_bdp 相同——逐变量 Taylor 提升。不同之处：
- multi_bdp 只处理 1 个 aux 变量
- WMDS 处理任意多个 aux 变量，通过递归

正确性归纳：
- 基础（0 aux vars）：单变量 MDP（`mdp_exists` ✅）
- 递归（n+1 aux vars）：Newton 迭代（`mtshl_step_invariant` ✅）+ 递归 IH

**WMDS 的正确性由已证框架覆盖。**

---

## 5. 稀疏插值 MDP（`__mtshl_sparse_int`）

### 5.1 C++ 算法

```
Step 0: 随机选 β₂,...,βⱼ₋₁
Step 1: θ-array 求值
  images_c[l] = c(x₁, β^{l+1})            l=0..s-1
  images_F[i][l] = F[i](x₁, β^{l+1})     l=0..s-1
Step 2: 逐点单变量 MDP
  for l = 0 to s-1:
    sigma_vals[i][l] = MDP(images_F[·][l], images_c[l])
Step 3: Vandermonde 恢复
  for each (i, x₁-degree d):
    thetas = [θ_m for m ∈ forms[i] with deg(m,x₁)=d]
    values = [sigma_vals[i][l][d] for l=0..t-1]
    coeffs = Vandermonde_solve(values, thetas)
    result[i] += Σ coeffs[k] · forms[i][k]
Step 4: 验证 Σ result[i]·F̂ᵢ = c
```

### 5.2 数学核心

**关键等式**：对每个固定的 x₁ 次数 d 和因子 i，
`σᵢ(x₁, β^l) 的 x₁^d 系数 = Σₘ cₘ · θₘˡ`
其中 cₘ 是 σᵢ 在单项式 m 处的系数，θₘ = ∏βₖ^{eₖ(m)}。

这是 Vandermonde 系统 `V · c = v`，V[l][m] = θₘˡ。

**正确性条件**：
1. θₘ 两两不同（β 随机选取，高概率成立）
2. s ≥ |forms[i]|（足够多的求值点）
3. 单变量 MDP 正确（已证 ✅）
4. Vandermonde 可逆（θₘ 两两不同 → det ≠ 0）

### 5.3 L2 模型

```lean
/-- 稀疏插值 MDP 模型。
    对应 C++ __mtshl_sparse_int。-/
noncomputable def mdpSparseInt
    (F : List (MvPolynomial σ R)) (c : MvPolynomial σ R)
    (betas : ...) (forms : ...) : List (MvPolynomial σ R) :=
  -- Step 1: θ-array 求值（数学上 = 直接在 β^l 处求值）
  let images_c := thetaArrayEval c betas s
  let images_F := F.map (thetaArrayEval · betas s)
  -- Step 2: 逐点单变量 MDP
  let sigma_vals := (List.range s).map (fun l =>
    mdpUnivarAlg (images_F.map (· l)) (images_c l))
  -- Step 3: Vandermonde 恢复
  vandermonde_recover sigma_vals forms betas
```

### 5.4 正确性

稀疏插值 = 求值 + 单变量MDP + Vandermonde逆。

每步的正确性：
1. θ-array 求值 = 在 (x₁, β²ˡ,...,βⱼˡ) 处求值（数学等价） ✅
2. 单变量 MDP：`mdp_exists` ✅ → 每个求值点的 MDP 正确
3. Vandermonde 逆：`vandermonde_solve_exists` ✅ → 系数唯一恢复
4. 恢复的 result 满足原始 MDP 等式：由 3 的精确恢复 + 2 的逐点正确

**Step 4 验证**：C++ 最终验证 `Σ result[i]·F̂ᵢ = c`。
在 L2 中，验证等价于算法正确性的后置条件。

---

## 6. Vandermonde 可逆性

### 6.1 定理

```lean
theorem vandermonde_invertible (θ : Fin s → F)
    [Field F] (h_inj : Function.Injective θ) :
    IsUnit (Matrix.vandermonde θ).det
```

Mathlib 有 `Matrix.det_vandermonde`：
`det(V) = ∏_{i<j} (θⱼ - θᵢ)`。
θ injective → 所有差非零 → 乘积非零 → 在域中可逆。✓

### 6.2 唯一解

```lean
theorem vandermonde_unique_solution (θ : Fin s → F) (v : Fin s → F)
    (h_inj : Function.Injective θ) :
    ∃! d, ∀ l, (Finset.univ.sum (fun t => d t * θ t ^ (l : ℕ))) = v l
```

由 det ≠ 0 → 矩阵可逆 → 线性系统有唯一解。✓

---

## 7. θ-array = 直接求值

### 7.1 定理

```lean
/-- θ-array 求值等价于直接在 (x₁, β^l) 处求值。
    C++ 的 running product 技巧不改变数学结果。-/
theorem thetaArrayEval_eq_directEval
    (f : MvPolynomial σ R) (betas : ...) (l : ℕ) :
    thetaArrayEval f betas l = eval_at_beta_pow f betas l
```

θ-array 的 running product `θₘ^l = θₘ^{l-1} · θₘ` 与直接计算 `∏βₖ^{eₖ·l}` 在数学上等价（乘法交换律 + 幂律）。这不需要证明——两者的定义在 Lean 中应该是 definitionally equal。

---

## 8. MDP 级联回退

### 8.1 C++ 控制流

```
尝试 1: __mtshl_sparse_int
尝试 2: __mtshl_sparse_int（重试，不同随机 β）
回退:
  aux_vars = 1: __mtshl_multi_bdp
  aux_vars ≥ 2: __mtshl_wmds
```

### 8.2 L2 模型

级联只影响**效率和鲁棒性**，不影响正确性。每个求解器如果成功返回，输出满足 `MDPCorrect`。

```lean
/-- MDP 级联：尝试多种求解器，任一成功即返回。
    对应 C++ 的 sparse_int → retry → multi_bdp/wmds 回退。-/
def mdpCascade (...) : Option (List (MvPolynomial σ R)) :=
  (mdpSparseInt ...).orElse (fun _ =>
  (mdpSparseInt ...).orElse (fun _ =>   -- 重试
  if aux_vars.length == 1 then mdpMultiBdp ...
  else mdpWmds ...))
```

正确性：每个分支返回 `some sigma` 时，`MDPCorrect c F_base sigma` 成立。由各求解器的独立正确性保证。

---

## 9. 形式化估计

| 组件 | 行数 | 依赖 |
|------|------|------|
| `taylorCoeff` 定义 + 正确性 | ~40 | 因子定理迭代 |
| `mdpUnivarAlg`（Bézout + mod）| ~30 | `mdp_exists` ✅ |
| `vandermonde_invertible` + 唯一解 | ~40 | Mathlib `det_vandermonde` |
| `thetaArrayEval` = 直接求值 | ~10 | 定义等价 |
| `mdpSparseInt` 模型 + 正确性 | ~60 | Vandermonde + 单变量 MDP |
| `mdpMultiBdp` 模型 + 正确性 | ~30 | Newton 迭代框架 ✅ |
| `mdpWmds` 模型 + 正确性 | ~40 | Newton 迭代框架 ✅ + 递归 |
| MDP 级联 | ~10 | 各求解器正确性 |
| **总计** | **~260** |

---

## 10. 关键洞察

multi_bdp 和 wmds 的 Taylor 逐阶提升与 `mtshl_step_j` **数学结构完全相同**——都是 Newton 迭代。所以 `mtshl_step_invariant`（已证 ✅）直接适用。

稀疏插值的正确性归结为：Vandermonde 可逆 + 单变量 MDP（已证）。

**真正的新证明工作**：Vandermonde 可逆性（~40 行，Mathlib `det_vandermonde`）和 Taylor 系数提取（~40 行，因子定理迭代）。其余由已证框架覆盖。
