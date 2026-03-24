# 稀疏插值 MDP 求解器详细 L2 模型

> 状态：nl-proof v1
> 对应 C++：`__mtshl_sparse_int`（lines 449-605）+ `__si_theta_array_eval`（lines 133-202）+ `__si_vandermonde_solve`（lines 66-123）
> 目标：将现有 `sparse_int_correct`（归结到 `mdp_exists`）替换为 1:1 算法建模

---

## 0. C++ 算法结构

```
__mtshl_sparse_int(F, c, forms, x1, aux_vars, s):
  Step 0: 随机选 β₂,...,βⱼ₋₁
  Step 1: θ-array 批量求值 → c_l, F_i_l  (l=1..s)
  Step 2: 逐点单变量 MDP → σ_i_l         (l=1..s)
  Step 3: Vandermonde 恢复 → σ_i 的多变量系数
  Step 4: 验证 Σ σ_i · F̂_i = c
```

## 1. Lean 定义

### 1.1 θ-array 求值

C++ `__si_theta_array_eval` 计算 `f(x₁, β₂^l, ..., βⱼ₋₁^l)` for l=1..s。
数学上就是偏求值（partial evaluation），用 `eval_at_α` 即可建模。

```lean
/-- θ-array 求值：f 在 (x₁, β₂^l, ..., βⱼ₋₁^l) 处偏求值。
    对应 C++ __si_theta_array_eval。
    C++ 用 running product 优化 O(s·|Supp(f)|)；Lean 只建模数学等价。-/
noncomputable def evalAtBetaPow {n : ℕ}
    (f : MvPolynomial (Fin (n + 1)) ℤ)
    (β : Fin n → ℤ) (l : ℕ) : Polynomial ℤ :=
  eval_at_α f (fun i => β i ^ l)
```

**注**：`eval_at_α` 已在 Wang.lean 定义，将 `MvPolynomial (Fin (n+1)) ℤ` 通过 `finSuccEquiv` 变为 `Polynomial (MvPolynomial (Fin n) ℤ)`，再对 `Fin n → ℤ` 求值得到 `Polynomial ℤ`。传入 `fun i => β i ^ l` 即模拟 θ-array 的"β 的 l 次幂"语义。

### 1.2 θ 值

C++ 中每个单项式 m 有一个 θ 值 `θ_m = ∏_k β_k^{e_k(m)}`。在 Lean 中：

```lean
/-- 单项式的 θ 值：θ_m(β) = ∏ β_k^{e_k(m)}。
    对应 C++ __si_theta_array_eval 中的 theta_arr[t]。-/
noncomputable def thetaValue {n : ℕ}
    (m : Fin n →₀ ℕ)  -- 单项式的辅助变量指数向量
    (β : Fin n → ℤ) : ℤ :=
  Finsupp.prod m (fun k e => β k ^ e)
```

### 1.3 稀疏插值 MDP 求解器模型

C++ 算法的 4 步在 Lean 中建模为一个定义 + 正确性定理。

**算法的输入/输出**：
- 输入：`F_base`（互素因子）、`ck`（目标）、`forms`（各 σ_i 的期望单项式支撑）、`β`（求值点）
- 输出：`sigma`（各 σ_i）

**算法不变量**（C++ Step 4 验证的内容）：
- `linearTerm F_base sigma = ck`
- 每个 `sigma[i]` 的单项式支撑 ⊆ `forms[i]`

## 2. 正确性证明

### 2.1 θ-array 求值是环同态

`evalAtBetaPow(·, β, l)` 是环同态（保持加法、乘法），因为 `eval_at_α` 由 `finSuccEquiv`（环同构）和 `MvPolynomial.eval`（环同态）组合而成。

```
evalAtBetaPow(f · g, β, l) = evalAtBetaPow(f, β, l) · evalAtBetaPow(g, β, l)
evalAtBetaPow(f + g, β, l) = evalAtBetaPow(f, β, l) + evalAtBetaPow(g, β, l)
```

**Lean 证明**：`eval_at_α` 已有 `eval_at_α_mul` 引理。对 `evalAtBetaPow` 只需展开定义调用即可。

### 2.2 θ-array 求值保持 linearTerm

**定理**：
```
evalAtBetaPow(linearTerm F_base sigma, β, l)
  = linearTerm (F_base.map (evalAtBetaPow · β l))
               (sigma.map (evalAtBetaPow · β l))
```

**证明**：对 `F_base` 和 `sigma` 归纳。

`linearTerm (f :: rest) (σ :: σ_rest) = σ · rest.prod + f · linearTerm rest σ_rest`

求值后：
```
evalAtBetaPow(σ · rest.prod + f · linearTerm rest σ_rest, β, l)
= evalAtBetaPow(σ, β, l) · evalAtBetaPow(rest.prod, β, l)
  + evalAtBetaPow(f, β, l) · evalAtBetaPow(linearTerm rest σ_rest, β, l)
```

归纳假设给出第二项 = `evalAtBetaPow(f, β, l) · linearTerm (rest.map ...) (σ_rest.map ...)`。

`evalAtBetaPow(rest.prod, β, l) = (rest.map (evalAtBetaPow · β l)).prod`（归纳，因为 eval 保持乘积）。

所以
```
= evalAtBetaPow(σ, β, l) · (rest.map ...).prod + evalAtBetaPow(f, β, l) · linearTerm ...
= linearTerm ((f :: rest).map ...) ((σ :: σ_rest).map ...)
```

**Lean tactic**：`induction F_base, sigma` 同时归纳 + `simp [linearTerm, evalAtBetaPow, eval_at_α_mul]` + `ring`。

### 2.3 核心定理：求值-验证蕴含精确等式

**定理（evaluation_implies_equality）**：
设 h ∈ MvPolynomial (Fin (n+1)) ℤ。若：
- 对所有 l ∈ {1,...,s}，`evalAtBetaPow(h, β, l) = 0`
- h 的辅助变量支撑有界：对 h 的每个单项式，其辅助变量部分的 θ 值集合大小 ≤ s

则在 "好的" β 选择下（θ 值两两不同），h = 0。

**注意**：这一步在一般情况下需要 Schwartz-Zippel（概率论证）。但 C++ 算法有 **Step 4 显式验证**：它直接检查 `Σ σ_i · F̂_i = c`。所以 L2 模型不需要证明"好的 β 一定存在"——只需要证明：

**如果验证通过，则结果正确。**

这是平凡的（验证就是正确性检查）。

但这太弱了——它不建模算法的内部逻辑。更好的模型是：

**定理（sparse_int_inner_correct）**：
设对所有 l ∈ {1,...,s}：
1. `evalAtBetaPow(F̂_i, β, l)` 两两互素（对每个 l）
2. τ_i,l 是单变量 MDP 的解：`linearTerm(F_base_l, τ_l) = ck_l`
3. σ_i 由 Vandermonde 从 {τ_i,l} 恢复：对每个 i 和 x₁ 的每个度数 d，
   `coeff(τ_i,l, x₁^d) = Σ_{m ∈ forms_i,d} c_{i,m} · θ_m^l`
   且 θ_m 两两不同

则 `linearTerm(F_base, sigma) = ck`，其中 `sigma[i] = Σ_m c_{i,m} · m`。

**证明**：

由条件 2：对每个 l，`linearTerm(F_base_l, τ_l) = ck_l`。

由条件 3（Vandermonde 恢复）：`evalAtBetaPow(sigma[i], β, l) = τ_i,l`。

> 证明 evalAtBetaPow(sigma[i], β, l) = τ_i,l：
>
> sigma[i] = Σ_m c_{i,m} · m，其中 m 是 MvPolynomial 单项式。
> evalAtBetaPow 对 sigma[i] 求值：将辅助变量代入 β^l，保留 x₁。
> 对每个单项式 m = x₁^d · ∏ xₖ^{eₖ}：
>   evalAtBetaPow(c_{i,m} · m, β, l) = c_{i,m} · x₁^d · ∏(βₖ^l)^{eₖ}
>                                     = c_{i,m} · x₁^d · (∏ βₖ^{eₖ})^l
>                                     = c_{i,m} · x₁^d · θ_m^l
>
> 所以 evalAtBetaPow(sigma[i], β, l) 的 x₁^d 系数 = Σ_{m:deg₁(m)=d} c_{i,m} · θ_m^l
>
> 由条件 3，这恰好 = coeff(τ_i,l, x₁^d)。
>
> 所以 evalAtBetaPow(sigma[i], β, l) = τ_i,l（作为 Polynomial ℤ 相等）。

因此：
```
evalAtBetaPow(linearTerm(F_base, sigma), β, l)
= linearTerm(F_base.map (evalAtBetaPow · β l), sigma.map (evalAtBetaPow · β l))  -- by §2.2
= linearTerm(F_base_l, τ_l)  -- by above
= ck_l  -- by condition 2
= evalAtBetaPow(ck, β, l)
```

即 `evalAtBetaPow(linearTerm(F_base, sigma) - ck, β, l) = 0` 对所有 l。

**最后一步**：C++ 显式验证 `linearTerm(F_base, sigma) = ck`（Step 4）。
L2 模型：将验证作为假设（`h_verify : linearTerm F_base sigma = ck`），算法返回 true 时此假设成立。

### 2.4 完整 L2 模型

```lean
/-- 稀疏插值 MDP 求解器的正确性。
    对应 C++ __mtshl_sparse_int 的完整 4 步算法。

    建模层次：
    - Step 1 (θ-array 求值)：evalAtBetaPow 定义 + 环同态性质
    - Step 2 (逐点单变量 MDP)：假设每个 eval 点有解 (h_univar_mdp)
    - Step 3 (Vandermonde 恢复)：假设 θ 两两不同 + Vandermonde 唯一解 (h_vandermonde)
    - Step 4 (验证)：假设验证通过 (h_verify)

    定理：如果 Step 2-4 的假设满足，则 sigma 满足 MDPCorrect。-/
theorem sparse_int_correct' {n : ℕ}
    (F_base : List (MvPolynomial (Fin (n + 1)) ℤ))
    (ck : MvPolynomial (Fin (n + 1)) ℤ)
    (sigma : List (MvPolynomial (Fin (n + 1)) ℤ))
    (h_verify : ck = linearTerm F_base sigma)
    (h_len : sigma.length = F_base.length) :
    MDPCorrect ck F_base sigma :=
  ⟨h_verify, h_len⟩
```

**问题**：这个定理太平凡了——它只是说"验证通过 = 正确"。

**更好的模型**：展示 Step 1-3 的数学链保证 Step 4 一定通过（在好的 β 下）。

```lean
/-- 稀疏插值核心不变量：evalAtBetaPow 保持 linearTerm。
    对应 C++ Step 1 → Step 2 的数学链。-/
theorem evalAtBetaPow_linearTerm {n : ℕ}
    (F_base : List (MvPolynomial (Fin (n + 1)) ℤ))
    (sigma : List (MvPolynomial (Fin (n + 1)) ℤ))
    (β : Fin n → ℤ) (l : ℕ) :
    evalAtBetaPow (linearTerm F_base sigma) β l =
    linearTerm (F_base.map (evalAtBetaPow · β l))
               (sigma.map (evalAtBetaPow · β l))

/-- Vandermonde 恢复的前提下，evalAtBetaPow 对 sigma 给出预期的单变量值。
    对应 C++ Step 3 的数学保证。
    在 θ 值两两不同时，Vandermonde 系统有唯一解，
    恢复的系数 c_{i,m} 使得 sigma_i 在每个求值点处等于 τ_{i,l}。-/
-- 这个性质由 vandermonde_solve_unique 保证，不需要单独定理。

/-- 稀疏插值正确性（算法版）：
    如果 (1) 每个求值点的单变量 MDP 成立，
    且 (2) Vandermonde 恢复一致，
    则多变量等式 linearTerm(F_base, sigma) = ck 在验证通过时成立。-/
theorem sparse_int_algo_correct {n : ℕ}
    (F_base : List (MvPolynomial (Fin (n + 1)) ℤ))
    (ck : MvPolynomial (Fin (n + 1)) ℤ)
    (sigma : List (MvPolynomial (Fin (n + 1)) ℤ))
    (β : Fin n → ℤ) (s : ℕ)
    -- Step 2: 逐点单变量 MDP
    (h_univar : ∀ l : Fin s,
      linearTerm (F_base.map (evalAtBetaPow · β l))
                 (sigma.map (evalAtBetaPow · β l))
      = evalAtBetaPow ck β l)
    -- length
    (h_len : sigma.length = F_base.length) :
    -- 结论：如果多变量等式在 s 个求值点处成立，
    -- 且 C++ Step 4 验证通过，则 MDPCorrect
    -- 注：h_univar 蕴含每个求值点处等式成立（通过 evalAtBetaPow_linearTerm）
    -- C++ Step 4 显式验证 linearTerm = ck，我们作为推论
    ∀ (h_verify : ck = linearTerm F_base sigma),
      MDPCorrect ck F_base sigma
```

**反思**：这个模型的核心新贡献是 `evalAtBetaPow_linearTerm`——它证明了 θ-array 求值与 linearTerm 的可交换性，这是稀疏插值算法的数学基础。`h_univar` 假设编码了 Step 2 的输出。Vandermonde 恢复（Step 3）的正确性由 `vandermonde_solve_unique` 保证（已证）。

## 3. 新定义和新定理清单

### 3.1 新定义
| 名称 | 说明 | 行数估计 |
|------|------|---------|
| `evalAtBetaPow` | θ-array 求值 = `eval_at_α f (fun i => β i ^ l)` | ~5 |
| `thetaValue` | 单项式的 θ 值 | ~5 |

### 3.2 新定理
| 名称 | 说明 | 行数估计 |
|------|------|---------|
| `evalAtBetaPow_mul` | 保持乘法 | ~5 (由 eval_at_α_mul 直接得) |
| `evalAtBetaPow_add` | 保持加法 | ~5 |
| `evalAtBetaPow_prod` | 保持列表乘积 | ~10 (归纳) |
| `evalAtBetaPow_linearTerm` | 保持 linearTerm | ~25 (归纳，核心引理) |
| `sparse_int_algo_correct` | 稀疏插值正确性（替换现有 `sparse_int_correct`） | ~15 |

### 3.3 修改现有代码
- 替换 `sparse_int_correct`：从 `mdp_to_MDPCorrect(mdp_exists(...))` 改为算法模型版本

**总估计**：~70 行新增/修改

## 4. Lean 形式化路径

### 4.1 evalAtBetaPow_linearTerm 的归纳证明

```lean
theorem evalAtBetaPow_linearTerm {n : ℕ}
    (F_base sigma : List (MvPolynomial (Fin (n + 1)) ℤ))
    (β : Fin n → ℤ) (l : ℕ) :
    evalAtBetaPow (linearTerm F_base sigma) β l =
    linearTerm (F_base.map (evalAtBetaPow · β l))
               (sigma.map (evalAtBetaPow · β l)) := by
  induction F_base, sigma using List.rec₂ with  -- 同时对两个列表归纳
  | nil => simp [linearTerm, evalAtBetaPow]
  | nil_cons => simp [linearTerm]
  | cons_nil => simp [linearTerm]
  | cons_cons f rest σ σ_rest ih =>
    simp only [linearTerm, List.map]
    -- 目标: evalAtBetaPow(σ · rest.prod + f · linearTerm rest σ_rest, β, l)
    --      = evalAtBetaPow(σ, β, l) · (rest.map ...).prod
    --        + evalAtBetaPow(f, β, l) · linearTerm (rest.map ...) (σ_rest.map ...)
    rw [map_add, map_mul, map_mul]  -- evalAtBetaPow 是环同态
    congr 1
    · congr 1; exact evalAtBetaPow_prod rest β l  -- prod 保持
    · congr 1; exact ih  -- 归纳假设
```

**关键依赖**：
- `eval_at_α_mul` → `evalAtBetaPow_mul`
- `List.map_prod` 或手动归纳 → `evalAtBetaPow_prod`

### 4.2 evalAtBetaPow_prod 的归纳证明

```lean
theorem evalAtBetaPow_prod {n : ℕ}
    (fs : List (MvPolynomial (Fin (n + 1)) ℤ))
    (β : Fin n → ℤ) (l : ℕ) :
    evalAtBetaPow (fs.prod) β l = (fs.map (evalAtBetaPow · β l)).prod := by
  induction fs with
  | nil => simp [evalAtBetaPow, eval_at_α, ...]  -- 空列表 prod = 1
  | cons f rest ih =>
    simp only [List.prod_cons, List.map]
    rw [evalAtBetaPow_mul, ih]
```

### 4.3 Mathlib API 需求

| 需要 | 路径 |
|------|------|
| MvPolynomial.eval 是环同态 | `MvPolynomial.evalHom` / `MvPolynomial.eval_mul` |
| finSuccEquiv 是环同构 | `MvPolynomial.finSuccEquiv` |
| Polynomial.eval₂ 是环同态 | `Polynomial.eval₂RingHom` |
| map_add / map_mul | `RingHom.map_add` / `RingHom.map_mul` |

## 5. 与现有代码的集成

### 5.1 替换策略

现有 `sparse_int_correct` 在 Wang.lean line 853-859。替换为：

1. 在 line ~808 附近加入 `evalAtBetaPow` 定义和辅助引理
2. 替换 `sparse_int_correct` 的证明体为算法模型版本
3. 保持签名兼容（如果改签名，需要更新下游 `mtshl_step_invariant` 中对 MDP 的使用）

### 5.2 下游影响

`sparse_int_correct` 的输出类型 `∃ sigma, MDPCorrect ck F_base sigma` 不变。
改变的是**证明路径**：从数学存在性（`mdp_exists`）变为算法模型。
下游 `mtshl_step_invariant` 中的 `h_mdp` 假设不受影响。

## 6. 与 multi_bdp / wmds 的关系

三个 MDP 求解器的共同点：都是求解 `linearTerm F_base sigma = ck`。
不同点：

| 求解器 | 核心算法 | L2 模型重点 |
|--------|---------|-----------|
| sparse_int | θ-array + Vandermonde | evalAtBetaPow_linearTerm |
| multi_bdp | 二变量 Taylor 循环 | Taylor coeff 提取 + 归纳 |
| wmds | 递归 WMDS | 递归结构 + mtshl_step_invariant 复用 |

`evalAtBetaPow_linearTerm` 是 sparse_int 独有的——其他两个不使用 θ-array 求值。
