# EDF 分裂存在性：消除 splits_fn 假设

> 状态：nl-proof v1
> 目标：证明好的分裂元素存在，消除 edf_correct 中的 hsplits 假设

---

## 0. 当前问题

`edf_correct` 需要外部假设 `hsplits : ∀ q ∈ edf f d splits, q.natDegree ≤ d`。
这等价于假设 splits 列表中有足够好的分裂元素。

**目标**：证明对于任何有 ≥ 2 个不可约因子的等度多项式 f，**存在** T 使得
`gcd(T^((p^d-1)/2) - 1, f)` 给出非平凡分裂。

---

## 1. 数学证明

### 1.1 设定

设 p 为奇素数，f ∈ F_p[x] monic squarefree，f = g₁ · g₂ · ... · gᵣ（r ≥ 2），
每个 gᵢ monic irreducible of degree d。

### 1.2 关键事实

**事实 1**：F_p[x]/(gᵢ) ≃ F_{p^d}（有限域，阶 p^d）。

**事实 2**：F_{p^d}* 的阶为 p^d - 1。对任何 a ∈ F_{p^d}*，a^((p^d-1)/2) = ±1
（因为 (a^((p^d-1)/2))² = a^(p^d-1) = 1，在域中 x² = 1 只有 ±1 两个根）。

**事实 3**（Euler 判别法）：a^((p^d-1)/2) = 1 iff a 是 F_{p^d} 中的二次剩余，
= -1 iff a 是非二次剩余。二次剩余恰好 (p^d-1)/2 个，非二次剩余也恰好 (p^d-1)/2 个。
（F_{p^d}* 是循环群，取生成元 γ，γ^k 是二次剩余 iff k 偶。）

**事实 4**（CRT）：由于 g₁,...,gᵣ 两两互素：
F_p[x]/(f) ≃ F_p[x]/(g₁) × ... × F_p[x]/(gᵣ) ≃ F_{p^d} × ... × F_{p^d}

CRT 同构是**满射**：对任何 (a₁,...,aᵣ) ∈ ∏ F_p[x]/(gᵢ)，∃ T ∈ F_p[x] 使得
T ≡ aᵢ (mod gᵢ) 对所有 i。

### 1.3 分裂存在性证明

要证：∃ T 使得 gcd(T^((p^d-1)/2) - 1, f) 是 f 的非平凡因子。

构造：
1. 取 α ∈ F_{p^d} 为非二次剩余（存在，因为恰好一半是非二次剩余，p^d ≥ 3）。
2. 取 β = 1 ∈ F_{p^d}（二次剩余）。
3. 由 CRT，∃ T ∈ F_p[x] 使得：
   - T ≡ α (mod g₁)（α 是非二次剩余）
   - T ≡ β (mod g₂)（β = 1 是二次剩余）
   - T ≡ 任意 (mod gᵢ, i ≥ 3)（无所谓）

4. 则：
   - T^((p^d-1)/2) ≡ α^((p^d-1)/2) ≡ -1 (mod g₁)（非二次剩余 → -1）
   - T^((p^d-1)/2) ≡ β^((p^d-1)/2) ≡ 1 (mod g₂)（二次剩余 → 1）

5. 所以：
   - g₁ ∤ (T^((p^d-1)/2) - 1)（因为 ≡ -2 ≠ 0 mod g₁，p 奇 → 2 ≠ 0）
   - g₂ | (T^((p^d-1)/2) - 1)（因为 ≡ 0 mod g₂）

6. 令 h = gcd(T^((p^d-1)/2) - 1, f)：
   - g₂ | h（g₂ | 两者）
   - g₁ ∤ h（g₁ ∤ T^((p^d-1)/2) - 1）
   - 所以 h 是 f 的非平凡因子：h ≠ 1（g₂ | h → deg h ≥ d > 0）且 h ≠ f（g₁ ∤ h）✓

### 1.4 p = 2 的情况

p = 2 时 (p^d - 1)/2 = (2^d - 1)/2 不是整数。C++ 用整数除法下取整。
实际上 p = 2 时 Cantor-Zassenhaus 用不同的公式：T + T² + T⁴ + ... + T^{2^{d-1}}。
但 CLPoly 的 `__select_prime` 选择奇素数，所以 p = 2 在实践中不出现。

**处理**：L2 定理加 p ≥ 3（奇素数）条件。这匹配 C++ 实践。

---

## 2. Lean 形式化路径

### 2.1 定理签名

```lean
theorem edf_split_exists
    (f : Polynomial (ZMod p)) (d : ℕ)
    (hp_odd : p ≥ 3)
    (hf_monic : Monic f) (hf_sqfree : Squarefree f)
    (hd : 0 < d)
    (hf_factors : ∀ q, Irreducible q → q ∣ f → q.natDegree = d)
    (hf_deg : d < f.natDegree)  -- f 有 ≥ 2 个不可约因子
    : ∃ T : Polynomial (ZMod p),
        let g := normalize (EuclideanDomain.gcd (T ^ ((p ^ d - 1) / 2) - 1) f)
        0 < g.natDegree ∧ g.natDegree < f.natDegree
```

### 2.2 依赖的 Mathlib API

| 需要 | Mathlib 路径 |
|------|------------|
| F_p[x]/(g) ≃ F_{p^d} 当 g irreducible degree d | `AdjoinRoot.powerBasis`、`Fintype.card_of_finset` |
| 有限域乘法群是循环群 | `IsCyclic (ZMod p)ˣ`、`FiniteField.isCyclic` |
| 二次剩余存在 / 非二次剩余存在 | 需要：∃ a ∈ F_{p^d}, a^((p^d-1)/2) = -1 |
| CRT for polynomials | `EuclideanDomain.quotient_map_prod` 或 `Ideal.quotientInfRingEquivPiQuotient` |
| gcd 性质 | `EuclideanDomain.gcd_dvd_left/right` |

### 2.3 核心子引理

```lean
-- 1. 在 F_{p^d} 中存在非二次剩余
lemma exists_nonsquare (q : ℕ) (hq : q ≥ 3) [Fintype (ZMod q)] :
    ∃ a : (ZMod q)ˣ, a ^ ((q - 1) / 2) = -1

-- 2. CRT 构造
lemma exists_poly_with_prescribed_residues
    (g₁ g₂ : Polynomial (ZMod p))
    (hcop : IsCoprime g₁ g₂)
    (a₁ a₂ : Polynomial (ZMod p)) :
    ∃ T, T ≡ a₁ [MOD g₁] ∧ T ≡ a₂ [MOD g₂]

-- 3. 组合：构造分裂元素
lemma cantor_zassenhaus_split ...
```

### 2.4 注意事项

- `exists_nonsquare` 对 F_{p^d}（不是 F_p）。F_{p^d} 在 Lean 中表示为 `AdjoinRoot g` 或 `GaloisField p d`。
  Mathlib 有 `GaloisField p d`（`FieldTheory.Galois`），但 API 可能不完整。
  替代：直接在 `Polynomial (ZMod p)` 上工作 modulo gᵢ。

- CRT for polynomials：Mathlib 有 `EuclideanDomain` 的 CRT 但可能不在 polynomial 特化形式。
  可以用 `IsCoprime.mul_dvd`：若 g₁ | r 且 g₂ | r 且 IsCoprime g₁ g₂，则 g₁*g₂ | r。
  反过来构造 T：从 Bézout s*g₁ + t*g₂ = 1，T = a₁*t*g₂ + a₂*s*g₁。

---

## 3. 形式化估计

| 内容 | 行数 |
|------|------|
| `exists_nonsquare`（有限域非二次剩余存在） | ~30 |
| CRT 构造 | ~20 |
| `edf_split_exists` 组合 | ~30 |
| 修改 `edf` 使用 Classical.choice | ~20 |
| 更新 `edf_correct` 消除 hsplits | ~10 |
| **总计** | **~110** |

---

## 4. 替代方案：更简单的存在性

如果完整的 Cantor-Zassenhaus 证明太复杂，有一个更简单的路径：

**引理**：f 有 ≥ 2 个不可约因子 → f 不是不可约的 → ∃ 非平凡因子。

这不需要 Cantor-Zassenhaus，但只给出**因子存在**，不给出**如何用 gcd 找到**。

如果只需要消除 `hsplits`（证明 edf 输出的每个元素 natDegree ≤ d），实际上我们需要的是：
edf 递归到底时，每个叶子的 natDegree ≤ d。这等价于说：对于每个 natDegree > d 的等度多项式，edf 能成功分裂它（不会在 splits 用尽时停留在 natDegree > d 的叶子上）。

所以我们确实需要 Cantor-Zassenhaus 的存在性（或足够长的 splits 列表）。

**但**：如果修改 `edf` 为不接受 `splits` 参数，而是用 `Classical.choice` 选择分裂元素（证明好的存在），则 `edf` 总是能递归到 natDegree ≤ d，不需要外部假设。
