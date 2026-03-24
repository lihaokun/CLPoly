# Vandermonde 求解 + θ-array 求值

> 状态：nl-proof v1
> 对应 C++：`__si_vandermonde_solve`（lines 66-123）+ `__si_theta_array_eval`（lines 133-209）
> 数学：Vandermonde 线性系统在 F_p 上的求解

---

## 0. 数学问题

### 0.1 Vandermonde 系统

给定 θ₁,...,θₛ ∈ F_p（两两不同）和 v₁,...,vₛ ∈ F_p，
求 d₁,...,dₛ 使得：

```
Σₜ dₜ · θₜˡ = vₗ    对所有 l = 0,...,s-1
```

即求解线性系统 V · d = v，其中 V 是 Vandermonde 矩阵 V[l][t] = θₜˡ。

### 0.2 θ-array 求值

给定 f ∈ F_p[x₁,...,xⱼ₋₁]（稀疏，s 项），求值点 β₂,...,βⱼ₋₁，
计算 f(x₁, β₂ˡ,...,βⱼ₋₁ˡ) 对 l = 1,...,s。

对每项 m = c · x₁^a · x₂^e₂ · ... · xⱼ₋₁^eⱼ₋₁：
定义 θ_m = ∏ₖ βₖ^eₖ（"θ-array" 值）。
则 m(x₁, β₂ˡ,...) = c · x₁^a · θ_m^l。

所以 f(x₁, β₂ˡ,...) = Σ_m c_m · x₁^{a_m} · θ_m^l。

固定 x₁ 的某个系数位置后，这是一个 Vandermonde 系统——未知数是 c_m，节点是 θ_m，右端是各 l 处的求值。

## 1. Vandermonde 可解性

### 1.1 定理

Vandermonde 矩阵在 θᵢ 两两不同时可逆（行列式 = ∏_{i<j}(θⱼ - θᵢ) ≠ 0）。

### 1.2 Lean 形式化路径

Mathlib 有 `Matrix.det_vandermonde`：
```
det (vandermonde θ) = ∏ i in Finset.univ, ∏ j in Finset.Ioi i, (θ j - θ i)
```

θᵢ 两两不同 + F_p 是域 → 行列式非零 → 矩阵可逆 → 唯一解存在。

### 1.3 Gauss-Jordan 正确性

C++ 用 Gauss-Jordan 消元求解。在域上 Gauss-Jordan 对可逆矩阵总是成功的，输出唯一解。

L2 模型：不需要建模 Gauss-Jordan 的每一步。只需证：
1. Vandermonde 矩阵可逆（两两不同节点）
2. 可逆矩阵有唯一解
3. Gauss-Jordan 输出的就是这个唯一解

即：`vandermonde_solve_correct : V 可逆 → V · d = v → d 是唯一的`。

C++ 的 Gauss-Jordan 是可逆矩阵求解的标准算法。它的正确性由线性代数保证，不是 MTSHL 特有的。

## 2. L2 模型

### 2.1 Vandermonde 求解

```lean
/-- Vandermonde 求解规约：给定两两不同的 θ 和目标 v，
    存在唯一 d 使得 V · d = v。
    对应 C++ __si_vandermonde_solve。-/
theorem vandermonde_solve_exists
    (θ : Fin s → ZMod p) (v : Fin s → ZMod p)
    (h_distinct : Function.Injective θ) :
    ∃! d : Fin s → ZMod p,
      ∀ l : Fin s, (Finset.univ.sum (fun t => d t * θ t ^ (l : ℕ))) = v l
```

### 2.2 θ-array 求值

```lean
/-- θ-array 求值：f(x₁, β^l) 对 l=1..s。
    对应 C++ __si_theta_array_eval。-/
noncomputable def thetaArrayEval
    (f : MvPolynomial (Fin j) (ZMod p))
    (betas : Fin (j-1) → ZMod p)
    (s : ℕ) : Fin s → Polynomial (ZMod p) :=
  fun l => eval_at_betas_pow f betas (l + 1)
```

这只是定义——θ-array 是一种**高效计算**求值的方法（O(s · nterms) 而非 O(s · nterms · log deg)），但数学上等价于直接求值。

## 3. 稀疏插值 MDP（__mtshl_sparse_int）

### 3.1 算法核心

```
输入：f_base[0..r-1], ck, forms[i] = Supp(F[i])
输出：sigma[0..r-1] 使得 ck = Σ σᵢ · F̂ᵢ

算法：
1. 选随机 β₂,...,βⱼ₋₁
2. 对 l = 1,...,s_max：
   - 求值 ck(x₁, β^l) 和 F_base[i](x₁, β^l)
   - 单变量 MDP 求解：ck_l = Σ σᵢ_l · F̂ᵢ_l
3. 对每个 σᵢ 的每个 x₁ 系数位置：
   - 收集 s 个求值 σᵢ_l[a]（l=1..s）
   - Vandermonde 求解恢复 σᵢ 的多变量系数
```

### 3.2 正确性

稀疏插值正确性由以下保证：
1. 对足够多的 l，θ-array 值两两不同（Schwartz-Zippel 引理）
2. Vandermonde 可逆 → 系数唯一恢复
3. 单变量 MDP 正确（j=2 Bézout，已证 ✅）
4. 逐系数恢复后的多变量 σᵢ 满足原始 MDP 等式

### 3.3 L2 建模层次

| 层次 | 建模 | 不建模 |
|------|------|--------|
| L2 | Vandermonde 可逆性 + 唯一解 | Gauss-Jordan 消元步骤 |
| L2 | θ-array 求值 = 直接求值 | O(s·nterms) 复杂度优化 |
| L2 | 稀疏插值 = Vandermonde 逆 + 求值 | 骨架估计、随机选取 |
| L1 | Gauss-Jordan 实现 | — |

## 4. 形式化估计

| 内容 | 行数 |
|------|------|
| `vandermonde_solve_exists`（Mathlib `det_vandermonde` + 可逆矩阵唯一解）| ~40 |
| `thetaArrayEval` 定义 + 等价于直接求值 | ~20 |
| `sparse_int_correct`（组合 Vandermonde + 单变量 MDP）| ~60 |
| `multi_bdp_correct`（二变量特例）| ~40 |
| `wmds_correct`（递归，复用单变量 MDP）| ~60 |
| **总计** | **~220** |
