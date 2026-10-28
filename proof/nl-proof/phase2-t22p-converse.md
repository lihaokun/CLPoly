# T2.2': 不可约多项式整除 X^{p^d} - X 的逆命题

> 状态：已形式化，0 sorry

---

## 定理陈述

设 p 为素数，g ∈ F_p[x] 不可约。若 g | X^{p^d} - X，则 deg(g) | d。

```lean
theorem dvd_X_pow_sub_X_imp_natDegree_dvd
    (g : Polynomial (ZMod p)) (hg : Irreducible g)
    (d : ℕ) (h : g ∣ (X ^ (p ^ d) - X : Polynomial (ZMod p))) :
    g.natDegree ∣ d
```

## 证明

设 K = AdjoinRoot g，k = deg(g)，α = root g，card(K) = p^k。

### Step A: α^{p^d} = α

由 g | X^{p^d} - X，在 K 中 mk g (X^{p^d} - X) = 0，即 α^{p^d} - α = 0。

**Mathlib 路径**：`AdjoinRoot.mk_eq_zero.mpr h` → `map_sub, map_pow, AdjoinRoot.mk_X` → `sub_eq_zero`

### Step B: ∀ x ∈ K, x^{p^d} = x

**关键思路**：构造迭代 Frobenius φ(x) = x^{p^d} 为 AlgHom over ZMod p，然后用 `AdjoinRoot.algHom_ext` 证明 φ = id。

1. φ 是环同态：`iterateFrobenius K p d` 给出 `K →+* K`
2. φ 是 ZMod p-代数同态：证明 `commutes'`，即 `(algebraMap r)^{p^d} = algebraMap r`
   - 因为 `r ∈ ZMod p`，`r^{p^d} = r`（`FiniteField.pow_card_pow` + `ZMod.card`）
3. φ(root g) = (root g)^{p^d} = root g = id(root g)（Step A）
4. 由 `AdjoinRoot.algHom_ext`：φ = id
5. 因此 ∀ x, x^{p^d} = x

**Mathlib 路径**：
- `iterateFrobenius (AdjoinRoot g) p d : K →+* K`（需要 `CharP K p`，由 `charP_of_injective_algebraMap` 获得）
- `AdjoinRoot.algHom_ext`：两个 AlgHom 在生成元上一致则相等
- `FiniteField.pow_card_pow d r` + `ZMod.card p`：ZMod p 中 `r^{p^d} = r`

**注**：最初考虑了三种方案（直接法/代换法/SubAlgebra），最终采用 AlgHom ext 方案——最简洁，只需约 15 行。

### Step C: k | d（反证法 + 根计数）

假设 k ∤ d。设 r = d mod k，0 < r < k。

1. ∀ x ∈ K, x^{p^k} = x（`FiniteField.pow_card_pow`，card K = p^k）
2. ∀ x ∈ K, x^{p^d} = x（Step B）
3. 由 p^d = (p^k)^{d/k} · p^r，得 x^{p^d} = (x^{(p^k)^{d/k}})^{p^r} = x^{p^r}
4. 从 x^{p^d} = x 和 x^{p^d} = x^{p^r}，得 ∀ x, x^{p^r} = x
5. X^{p^r} - X 在 K 中有 p^k 个根（每个元素都是根）
6. 但 deg(X^{p^r} - X) = p^r < p^k（因 r < k, p ≥ 2）
7. 域上度 n 多项式至多 n 个根（`Polynomial.card_roots'`）
8. p^k ≤ p^r 矛盾 p^k > p^r。故 r = 0，k | d。✓

**Mathlib 路径**：
- `Nat.div_add_mod d k`：`d = k * (d/k) + r`
- `pow_mul, pow_add`：指数算术
- `natDegree_sub_eq_left_of_natDegree_lt` + `natDegree_X_pow` + `natDegree_X`：度数计算
- `Polynomial.mem_roots` + `Polynomial.card_roots'`：根计数
- `Finset.card_le_card` + `Multiset.toFinset_card_le`：集合论不等式
- `Nat.pow_lt_pow_right`：p^k > p^r（k > r, p > 1）

## 实际行数

T2.2' 形式化约 55 行（不含注释），其中：
- Step A: 3 行
- Step B: 15 行（AlgHom 构造 + ext）
- Step C: 30 行（反证法 + 根计数）
- 前置设置: 7 行
