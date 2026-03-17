# T2.3: gcd(X^{p^d} - X, f) 的不可约因子刻画

> 状态：已形式化，0 sorry（无需 Monic/Squarefree 假设）

---

## 定理陈述

设 p 为素数，f ∈ F_p[x]，d ∈ ℕ，q 不可约。则：
q | gcd(X^{p^d} - X, f) ↔ q | f ∧ deg(q) | d

```lean
theorem gcd_X_pow_sub_X_factors
    (f : Polynomial (ZMod p)) (d : ℕ)
    (q : Polynomial (ZMod p)) (hq : Irreducible q) :
    q ∣ EuclideanDomain.gcd (X ^ (p ^ d) - X) f ↔ q ∣ f ∧ q.natDegree ∣ d
```

注：最终形式化去掉了 Monic/Squarefree 假设——定理对任意 f 成立。

## 证明

### 反方向（←）：q | f ∧ deg q | d → q | gcd

1. q | f（给定）
2. deg q | d → q | X^{p^d} - X（T2.2，需要 q 首一；不可约多项式可乘标量变首一，不影响整除性；或者直接要求 Monic 版本）
3. q | f ∧ q | X^{p^d} - X
4. q 不可约 → q 是素元（域上多项式环是 UFD，不可约=素）
5. q | gcd(X^{p^d} - X, f)（因为 gcd 是最大公因子，q 是两者的公因子 → q | gcd；或者更精确：在 PID 中，q | a ∧ q | b → q | gcd(a,b)）

**Mathlib 路径**：
- T2.2: `irreducible_dvd_X_pow_sub_X`
- `dvd_gcd`: `a | b → a | c → a | gcd b c`（EuclideanDomain.dvd_gcd 或类似）
- 首一化：不可约 q 在域上总可以首一化（`q * C (leadingCoeff q)⁻¹` 是首一的且 Associated）

### 正方向（→）：q | gcd → q | f ∧ deg q | d

1. q | gcd → q | f（因为 gcd | f，传递性）
2. q | gcd → q | X^{p^d} - X（因为 gcd | X^{p^d} - X，传递性）
3. **需要辅助定理 T2.2'**：q | X^{p^d} - X → deg q | d

**T2.2' 是 T2.2 的逆命题，需要单独证明。**

## T2.2' 的证明：q | X^{p^d} - X → deg q | d

### 证明思路

设 K = AdjoinRoot q，card(K) = p^k（k = deg q），α = root q。

**Step A**: q | X^{p^d} - X → α^{p^d} = α（AdjoinRoot.mk_eq_zero）

**Step B**: ∀ x ∈ K, x^{p^d} = x

论证：x ↦ x^{p^d} 是 K 的 Frobenius 迭代（环自同态）。α^{p^d} = α 且对所有 a ∈ F_p 有 a^{p^d} = a。由于 K = F_p(α)（PowerBasis），K 中每个元素是 α 的 F_p 多项式。环自同态固定 F_p 和 α，故固定所有 F_p(α) = K。

**Step C**: k | d（反证法）

假设 k ∤ d。写 d = qk + r, 0 < r < k。
- ∀ x ∈ K, x^{p^k} = x（FiniteField.pow_card）
- ∀ x ∈ K, x^{p^d} = x（Step B）
- x^{p^r} = (x^{p^{qk}})^{p^r}·... 实际上：
  x^{p^d} = x^{p^{qk+r}} = (x^{p^{qk}})^{p^r}
  由 pow_card_pow: x^{p^{qk}} = x
  所以 x^{p^d} = x^{p^r}
  从 x^{p^d} = x 得 x^{p^r} = x
- ∀ x ∈ K, x^{p^r} = x，即 X^{p^r} - X 在 K 中有 p^k 个根
- 但 deg(X^{p^r} - X) = p^r < p^k（因 r < k）
- 域上度 n 多项式至多 n 个根
- 矛盾。故 r = 0，k | d。✓

### Step B 的 Lean 形式化难点

Step B（Frobenius 固定生成元则固定整个域）是这个证明中最难形式化的部分。可能的 Mathlib 路径：

1. **直接法**：用 PowerBasis 把 x 写成 ∑ aᵢ αⁱ，然后用 Frobenius 的 RingHom 性质
2. **代换法**：f^{p^d} = f(X^{p^d}) 在 char p 下成立（Polynomial.expand 相关），X^{p^d} ≡ X (mod q)，所以 f^{p^d} ≡ f (mod q)
3. **SubAlgebra 法**：Fix(φ^d) 是 K 的子代数，包含 α 和 F_p，由 PowerBasis 生成的子代数 = K

**推荐方案**：方案 2 最直接——在多项式环层面操作，避免 AdjoinRoot 内部的元素分解。

## 实际行数

- T2.3 本身: ~12 行（正方向 dvd_trans，反方向 dvd_gcd）
- T2.2'（逆命题）: ~55 行（详见 `phase2-t22p-converse.md`）
- 总计: ~67 行

## 实施回顾

1. 先实现 T2.3 框架，T2.2' 用 sorry ✓
2. 填充 T2.2'：Step B 用 AlgHom ext（比预期简洁）✓
3. 全部 0 sorry ✓
