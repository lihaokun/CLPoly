# hensel_unique 泛化：Monic → 相同 leadingCoeff

> 状态：nl-proof v1
> 目标：将 hensel_unique 的 `Monic B1, Monic B2` 替换为更弱的条件

---

## 0. 动机

当前 `hensel_unique` 要求 B₁, B₂ monic in ℤ[x]。但在 factor_recovery 的对比对构造中，B' = h/(ℓ/ℓⱼ) 可能不 monic。泛化后可以处理 "两个 B 有相同 lc 且 p ∤ lc" 的情况。

---

## 1. 新签名

```lean
theorem hensel_unique
    (p : ℕ) (hp : Nat.Prime p) (k : ℕ) (hk : 0 < k)
    (F : Polynomial ℤ)
    (A1 B1 A2 B2 : Polynomial ℤ)
    (hprod1 : map_{p^k} (A1 * B1) = map_{p^k} F)
    (hprod2 : map_{p^k} (A2 * B2) = map_{p^k} F)
    (hA : map_p A1 = map_p A2)
    (hB : map_p B1 = map_p B2)
    (hcop : IsCoprime (map_p A1) (map_p B1))
    -- 泛化：Monic → 相同 lc + p 不整除 lc
    (hB_lc_eq : B1.leadingCoeff = B2.leadingCoeff)
    (hB_lc_coprime : ¬((p : ℤ) ∣ B1.leadingCoeff))
```

原 `Monic B1, Monic B2` 是特例：lc = 1, p ∤ 1。

---

## 2. 需要修改的辅助引理

### 2.1 `monic_diff_coeff_bound` → `lc_diff_coeff_bound`

**当前**：B₁, B₂ monic, B₁ - B₂ = C(m) * G → ∀ n ≥ natDegree(B₁), G.coeff n = 0

**泛化**：B₁.lc = B₂.lc, ¬(m | B₁.lc), B₁ - B₂ = C(m) * G → 同结论

**natDegree 相等证明的修改**：

当前（monic 版）：若 natDegree(B₂) < natDegree(B₁)：
- (B₁-B₂).coeff(natDegree B₁) = 1 - 0 = 1 = m · G.coeff(natDegree B₁)
- → m | 1, 矛盾（m > 1）✓

泛化版：若 natDegree(B₂) < natDegree(B₁)：
- (B₁-B₂).coeff(natDegree B₁) = lc(B₁) - 0 = lc(B₁) = m · G.coeff(natDegree B₁)
- → m | lc(B₁)

需要矛盾。m = p^k（在 hensel_unique 归纳中，inductive step 用 m = p^k）。
¬(p | lc(B₁)) → ¬(p^k | lc(B₁))（因 p | p^k | lc → p | lc）→ ¬(m | lc(B₁))。矛盾。✓

**修改量**：
- `hB1m.leadingCoeff` → `rfl`（用 hB_lc_eq）
- `one_ne_zero` 的矛盾 → `hB_lc_coprime` 的矛盾
- `Int.le_of_dvd one_pos` → 改用 `hB_lc_coprime` 得矛盾

**注意**：当前 `monic_diff_coeff_bound` 的矛盾用了 `Int.le_of_dvd one_pos ⟨_, hcoeff2⟩; omega`（m | 1 且 m > 1 → ⊥）。泛化后改为：m | lc(B₁) 且 ¬(m | lc(B₁)) → ⊥。

但 `monic_diff_coeff_bound` 的 `m` 参数不一定是 p^k——它是 hensel_unique 归纳步中的 p^k，但引理本身是独立的。所以需要在引理签名中接受 `¬((m : ℤ) ∣ B₁.leadingCoeff)` 作为条件。

### 2.2 `dvd_monic_eq_zero_of_natDegree_lt`

**当前**：b monic, b | q, natDegree(q) < natDegree(b) → q = 0

**需要**：b ≠ 0（不要求 monic），b | q, natDegree(q) < natDegree(b) → q = 0

**证明**：与现有相同。if q ≠ 0, then natDegree(q) = natDegree(b) + natDegree(c) ≥ natDegree(b)，矛盾。

但 `natDegree(b * c) = natDegree(b) + natDegree(c)` 需要 `b ≠ 0, c ≠ 0`（IsDomain 条件下成立）。当前证明用了 `hb_monic.natDegree_mul'` 需要 monic。

**替代路径**：直接用 `Polynomial.natDegree_le_of_dvd : p ∣ q → q ≠ 0 → natDegree p ≤ natDegree q`。

所以替换为：
```lean
by_contra hq
exact Nat.not_lt.mpr (Polynomial.natDegree_le_of_dvd hdvd hq) hdeg
```

这就不需要 monic 也不需要 `dvd_monic_eq_zero_of_natDegree_lt` 了。

### 2.3 hensel_unique 中使用 monic 的其他位置

在 hensel_unique 证明中，`hB1_monic` 还用于：
- `hBm : Monic (Polynomial.map ... B1) := hB1_monic.map _` — 给 G_bar 的度数约束
- `hBn : natDegree (map_p B1) = natDegree B1` — 从 monic 得 map 保持 natDegree

**泛化**：
- `Monic (map_p B1)`：不再有。但 `map_p B1 ≠ 0`（因 p ∤ lc(B₁) → lc(map_p B₁) ≠ 0 → map_p B₁ ≠ 0）。
- `natDegree (map_p B1) = natDegree B1`：从 `Polynomial.natDegree_map_of_leadingCoeff_ne_zero`。条件：`(Int.castRingHom (ZMod p)) B1.leadingCoeff ≠ 0`，等价于 `¬(p | lc(B₁))`，这正是我们的假设。✓

- `dvd_monic_eq_zero_of_natDegree_lt` 中 `hBm` 的使用：改用上述 `natDegree_le_of_dvd` 替代。

---

## 3. 修改计划

### 3.1 `monic_diff_coeff_bound` → `diff_coeff_bound`

```lean
private lemma diff_coeff_bound
    (B1 B2 : Polynomial ℤ) (m : ℕ) (hm : 1 < m)
    (G : Polynomial ℤ)
    (hB_lc_eq : B1.leadingCoeff = B2.leadingCoeff)
    (hB_lc_ne : ¬((m : ℤ) ∣ B1.leadingCoeff))
    (hdiff : B1 - B2 = C (↑m : ℤ) * G) :
    ∀ n, natDegree B1 ≤ n → G.coeff n = 0
```

### 3.2 hensel_unique

- 替换 `Monic B1, Monic B2` → `hB_lc_eq, hB_lc_coprime`
- `hBm` 的使用改为 `map_p B1 ≠ 0`
- `hBn` 用 `natDegree_map_of_leadingCoeff_ne_zero`
- G_bar = 0 的证明用 `natDegree_le_of_dvd` 代替 `dvd_monic_eq_zero_of_natDegree_lt`

### 3.3 保留原始版本为推论

```lean
theorem hensel_unique_monic := hensel_unique ... (by rw [hB1m.leadingCoeff, hB2m.leadingCoeff])
    (by rw [hB1m.leadingCoeff]; exact fun h => hp.one_lt.ne' (Int.eq_one_of_dvd_one ...))
```

---

## 4. 形式化估计

| 改动 | 行数变化 |
|------|---------|
| `diff_coeff_bound`（替换 `monic_diff_coeff_bound`） | ~5 行修改 |
| hensel_unique 签名 + 证明适配 | ~15 行修改 |
| **总计** | **~20 行修改**（净增很少） |
