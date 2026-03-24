# MTSHL 单步不变量传播 + MDP 求解

> 状态：nl-proof v1
> 对应 C++：`__mtshl_step_j` 内层循环（lines 866-932）

---

## 0. C++ 核心循环（单步 k）

```cpp
// error = aj - ∏F[i]（当前误差）
// xk_pow = (xj - αj)^k
ck = __taylor_coeff_zp(error, xj, αj, k);  // Taylor 系数
if (ck.empty()) continue;
sigma_k = MDP_solve(F_base, ck);             // Diophantine 求解
for i = 0 to r-1:
    F[i] += sigma_k[i] * xk_pow;            // 更新因子
error = aj - ∏F[i];                          // 重算误差
```

## 1. 数学论证

### 1.1 设定

设 d = X_j - C α_j。不变量：d^k | (target - ∏F_old)。

设 error = target - ∏F_old。则 error = d^k · q 对某个 q。

Taylor 系数 c_k = partialEval j α_j q（即 q 在 x_j = α_j 处求值）。

### 1.2 MDP（Multivariate Diophantine Problem）

给定 c_k 和基础因子 F_base[i] = F_old[i]|_{x_j=α_j}（在 α_j 处求值），
找 σ_i 使得：

```
c_k = Σ_i σ_i · F̂_i  (mod 关于 x₁ 的度数约束)
```

其中 F̂_i = ∏_{j≠i} F_base[j]。

这是偏分式分解：c_k / ∏F_base = Σ σ_i / F_base[i]。

### 1.3 更新后的误差

F_new[i] = F_old[i] + σ_i · d^k

∏F_new = ∏(F_old[i] + σ_i · d^k)

展开：
∏F_new = ∏F_old + d^k · (Σ_i σ_i · ∏_{j≠i} F_old[j]) + d^{2k} · (...)

所以：
∏F_new - ∏F_old = d^k · (Σ_i σ_i · F̂_i_full) + O(d^{2k})

其中 F̂_i_full = ∏_{j≠i} F_old[j]（不是 F_base 而是完整的 F_old）。

**关键**：F̂_i_full ≡ F̂_i (mod d)（因 F_old[j] ≡ F_base[j] (mod d)）。

所以：Σ σ_i · F̂_i_full ≡ Σ σ_i · F̂_i ≡ c_k (mod d)。

error_new = target - ∏F_new
         = error - d^k · (Σ σ_i · F̂_i_full) - O(d^{2k})
         = d^k · q - d^k · (Σ σ_i · F̂_i_full) - O(d^{2k})
         = d^k · (q - Σ σ_i · F̂_i_full) - O(d^{2k})

由 MDP：c_k = Σ σ_i · F̂_i = partialEval(q)。
所以 q - Σ σ_i · F̂_i_full ≡ q - c_k ≡ q - partialEval(q) (mod d)。
由因子定理：d | (q - partialEval(q))。
所以 q - Σ σ_i · F̂_i_full = d · r 对某个 r。

error_new = d^k · d · r + O(d^{2k}) = d^{k+1} · r + O(d^{2k})。
由 2k ≥ k+1（k ≥ 1）：d^{k+1} | O(d^{2k})。
所以 d^{k+1} | error_new。✓

## 2. Lean 形式化

### 2.1 MDP 规约（加强版）

```lean
/-- MDP 求解规约：c_k = Σ σ_i · F̂_i。
    F̂_i 是去掉第 i 个因子后的乘积在 α_j 处求值。-/
structure MDPCorrect {n : ℕ}
    (ck : MvPolynomial (Fin (n + 1)) ℤ)
    (f_base : List (MvPolynomial (Fin (n + 1)) ℤ))
    (sigma : List (MvPolynomial (Fin (n + 1)) ℤ)) : Prop where
  /-- 偏分式分解等式 -/
  partial_fraction : ck = (sigma.zipWith (· * ·)
    (f_base.enum.map (fun ⟨i, _⟩ =>
      (f_base.eraseIdx i).prod))).sum
  /-- 长度一致 -/
  length_eq : sigma.length = f_base.length
```

### 2.2 乘积展开引理

```lean
/-- 核心代数引理：∏(aᵢ + bᵢ) - ∏aᵢ ≡ Σ bᵢ · ∏_{j≠i} aⱼ (mod (∏bᵢ)²)。
    更精确：∏(aᵢ + bᵢ) = ∏aᵢ + Σ bᵢ · ∏_{j≠i} aⱼ + higher order terms。
    当 bᵢ = σᵢ · d^k 时，higher order 含 d^{2k}。-/
lemma prod_add_eq_prod_add_linear_term (...)
```

### 2.3 单步定理

```lean
theorem mtshl_step_invariant_full {n : ℕ}
    (target : MvPolynomial (Fin (n + 1)) ℤ)
    (factors : List (MvPolynomial (Fin (n + 1)) ℤ))
    (j : Fin n) (αⱼ : ℤ) (k : ℕ) (hk : 1 ≤ k)
    (sigma : List (MvPolynomial (Fin (n + 1)) ℤ))
    (h_inv : MtshlInvariant target factors j αⱼ k)
    -- MDP 正确：c_k = Σ σᵢ · F̂ᵢ
    (h_mdp : MDPCorrect
        (partialEval (Fin.succ j) αⱼ (error_quot target factors j αⱼ k h_inv))
        (factors.map (partialEval (Fin.succ j) αⱼ))
        sigma)
    (h_len : sigma.length = factors.length)
    : MtshlInvariant target
        (mtshlStepUpdate factors sigma
          ((MvPolynomial.X (Fin.succ j) - MvPolynomial.C αⱼ) ^ k))
        j αⱼ (k + 1)
```

其中 `error_quot target factors j αⱼ k h_inv` 是 `(target - ∏factors) / d^k`（从整除假设提取商）。

## 3. 证明路径

1. 设 d = X_j - C α, error = target - ∏factors = d^k · q（从 h_inv）
2. MDP：partialEval q = Σ σ_i · partialEval(F̂_i)
3. 乘积展开：∏(F_i + σ_i·d^k) - ∏F_i = d^k · (Σ σ_i · F̂_i) + d^{2k} · R
4. 新误差 = d^k · q - d^k · (Σ σ_i · F̂_i) - d^{2k} · R = d^k · (q - Σ σ_i · F̂_i) - d^{2k} · R
5. q - Σ σ_i · F̂_i ≡ q - partialEval(q) (mod d)（因 F̂_i ≡ partialEval(F̂_i) mod d... 需仔细）
6. d | (q - partialEval q)（因子定理 ✅ 已证）
7. d^{k+1} | d^k · d · (...) 且 d^{k+1} | d^{2k} · R（因 2k ≥ k+1）
8. 所以 d^{k+1} | 新误差 ✓

## 4. 形式化估计

| 内容 | 行数 |
|------|------|
| `MDPCorrect` 规约 | ~15 |
| `error_quot` 定义（从整除提取商） | ~10 |
| 乘积展开引理（2 因子版本 + 归纳到 r 因子） | ~60 |
| `mtshl_step_invariant_full` 证明 | ~50 |
| Taylor 循环 `mtshlLoop` 模型 | ~30 |
| `mtshl_loop_correct` 组合定理 | ~20 |
| **总计** | **~185** |

这是最后的工程量。完成后多变量 L2 全部 0 sorry。
