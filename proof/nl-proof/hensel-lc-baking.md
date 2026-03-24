# Hensel lc-baking：多因子 Hensel 提升的首项系数保持

> 状态：nl-proof v1
> 对应 C++：`polynomial_factorize_univar.hh:534-570`（lc-baking 初始化 + Hensel 提升）
> 参考：`formal-proof-univar-factorization.md` §2.1-2.3

---

## 0. 目标

证明：多因子 Hensel 提升后，H₁ 的首项系数 = lc(f)，Hᵢ (i ≥ 2) monic。

C++（line 562-568）：
```cpp
// lc-baking: 将 lc(f) 分配到 factor[0]
factors_adj[0] *= lc_mod_p;  // h̃₁ = ℓ · h̄₁ mod p
// factors[1..r] 不变（monic）
```

---

## 1. 不变量 (I3)

### 1.1 初始设定

f̄ = ℓ̄ · h̄₁ · h̄₂ · ... · h̄ᵣ（F_p[x]，h̄ᵢ monic irreducible）。

lc-baking 调整：
- h̃₁ = ℓ̄ · h̄₁（lc = ℓ̄ ≠ 0，不 monic）
- h̃ᵢ = h̄ᵢ（monic，i ≥ 2）
- 乘积：h̃₁ · h̄₂ · ... · h̄ᵣ = ℓ̄ · ∏h̄ᵢ = f̄ ✓

### 1.2 提升后不变量 (I3)

经过 Hensel 提升到 mod m = p^k：
- **(I3)** lc(H₁) = ℓ = lc(f)，lc(Hᵢ) = 1 (i ≥ 2)

### 1.3 (I3) 的证明

**关键**：`hensel_step_with_degree` 保持首项系数。

对于 2-factor 分裂 f ≡ G·H (mod m)，G 对应"含 lc" 的一侧，H 对应 monic 一侧。

hensel_step_with_degree 给出 (H5)：h'.natDegree = h.natDegree。
加上构造 h' = h + C(m)·σ_int，deg(σ_int) < natDegree(h)：
**lc(h') = lc(h)**（高次项未被修改）。

所以每步 Hensel 提升保持：
- monic 因子保持 monic（lc = 1 → lc = 1）
- 非 monic 因子保持 lc（lc = ℓ → lc = ℓ）

归纳：初始 lc(H̃₁) = ℓ，经 k 步提升后 lc(H₁) = ℓ。✓
初始 lc(H̃ᵢ) = 1，经 k 步提升后 lc(Hᵢ) = 1。✓

---

## 2. Lean 形式化

### 2.1 lc 保持引理（从 hensel_step_with_degree 推导）

```lean
/-- hensel_step_with_degree 保持 h 的首项系数 -/
theorem hensel_step_preserves_leadingCoeff
    (m : ℕ) (hm : 1 < m) (f g h : Polynomial ℤ)
    (hprod : ...) (hcop : ...) (hh_monic : ...) (hh_deg : ...)
    : ∃ g' h', [H1-H5 from hensel_step_with_degree]
      ∧ h'.leadingCoeff = h.leadingCoeff := by
  obtain ⟨g', h', h1, h2, h3, h4, h5⟩ := hensel_step_with_degree ...
  refine ⟨g', h', h1, h2, h3, h4, h5, ?_⟩
  -- h'.natDegree = h.natDegree (H5)
  -- h' = h + C(m)·σ_int, deg(C(m)·σ_int) < natDegree(h)
  -- → lc(h') = h.coeff(h.natDegree) = lc(h)
  sorry -- ~10 行：从 H5 + 构造推导
```

**问题**：当前 hensel_step_with_degree 是存在性定理（∃ g' h'），不暴露构造。
要推导 lc(h') = lc(h)，需要知道 h' 的具体构造（h' = h + C(m)·σ_int）。

**方案 A**：在 hensel_step_with_degree 内部直接加 (H6): h'.leadingCoeff = h.leadingCoeff。

**方案 B**：从 H5 + 定义推导——但存在性证明不暴露具体构造。

**选择方案 A**（最干净）：修改 hensel_step_with_degree 加入 H6。

### 2.2 修改 hensel_step_with_degree

在现有 H1-H5 基础上加：
```
∧ h'.leadingCoeff = h.leadingCoeff  -- H6
```

**证明**：在 hensel_step_with_degree 内部，h' = h + C(m)·σ_int。
- H5: h'.natDegree = h.natDegree
- 因此 h'.leadingCoeff = h'.coeff(h'.natDegree) = h'.coeff(h.natDegree)
- = h.coeff(h.natDegree) + (C(m)·σ_int).coeff(h.natDegree)
- = lc(h) + 0（因 deg(C(m)·σ_int) < natDegree(h)，该位置系数 = 0）
- = lc(h) ✓

### 2.3 多因子 lc-baking

有了 H6，多因子 lc-baking 就是归纳应用：

**定理**：设 f ≡ H₁·...·Hᵣ (mod p^k)，H₁ 的 lc = ℓ，Hᵢ monic (i ≥ 2)。
则经过进一步 Hensel 提升到 mod p^{k'}（k' > k），lc 和 monic 属性保持。

**证明**：每步 2-factor 分裂中，monic 侧（Hᵢ, i ≥ 2）经 hensel_step_with_degree H6 保持 lc = 1。
非 monic 侧（H₁ 或包含 H₁ 的乘积）保持 lc = ℓ。

（严格形式化需要对 binary tree 分裂结构做归纳，但数学上每步都是 H6 的应用。）

---

## 3. 形式化估计

| 改动 | 行数 |
|------|------|
| hensel_step_with_degree 加 H6 | ~15 |
| 多因子 lc-baking 引理（binary tree 归纳） | ~50（或 sorry） |
| **总计** | **~15-65** |

**建议**：先只加 H6 到 hensel_step_with_degree（~15 行）。多因子归纳作为 TODO。
