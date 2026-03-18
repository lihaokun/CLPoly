# T2.5: 半幂根计数 — Cantor-Zassenhaus 概率基础

> 状态：已形式化，0 sorry

---

## 定理陈述

设 F 为有限域，|F| = q，q 为奇数。则 X^{(q-1)/2} - 1 在 F 中恰有 (q-1)/2 个根。

等价地：在 F* 中，恰有 (q-1)/2 个元素 a 满足 a^{(q-1)/2} = 1。

```lean
theorem card_pow_half_eq_one
    {F : Type*} [Field F] [Fintype F] [DecidableEq F]
    (hF_odd : ringChar F ≠ 2) :
    (Finset.univ.filter (fun a : Fˣ => a ^ ((Fintype.card F - 1) / 2) = 1)).card
    = (Fintype.card F - 1) / 2
```

### 用途

这是 Cantor-Zassenhaus EDF 算法概率分析的核心代数事实。

**概率分析（纸面论证，不形式化）**：设 f = g₁·...·gₖ（k ≥ 2 个不同的 d 次不可约因子），
q = p^d。取随机 a mod f，由 CRT 对应 (a₁,...,aₖ) ∈ F_q^k。
每个 aᵢ^{(q-1)/2} 独立地为 0（概率 1/q）、1（概率 (q-1)/(2q)）、-1（概率 (q-1)/(2q)）。
gcd(h-1, f) 非平凡当且仅当不是所有因子都给出相同结果。

对 k = 2, q ≥ 3：Pr[非平凡] = (q²-1)/(2q²) ≥ 4/9。
对 k ≥ 3：概率单调递增趋向 1。

## 证明

### 方案：多项式根计数（避免循环群理论）

**Step 1**: 记 m = (q-1)/2。X^{q-1} - 1 分解为两个互素因子。

在任意交换环中：X^{2m} - 1 = (X^m - 1)(X^m + 1)。
（差平方分解：a² - b² = (a-b)(a+b)，取 a = X^m, b = 1）

验证 gcd(X^m - 1, X^m + 1) = 1：
它们的差 = (X^m + 1) - (X^m - 1) = 2。
char F ≠ 2 → 2 是单位 → gcd 是 1。
所以 X^m - 1 和 X^m + 1 互素。

**Step 2**: X^{q-1} - 1 在 F 中恰有 q-1 个根。

由 Fermat 小定理（`FiniteField.pow_card_sub_one_eq_one`）：
∀ a ∈ F*, a^{q-1} = 1。
所以 F* 的全部 q-1 个元素都是 X^{q-1}-1 的根。

又 deg(X^{q-1}-1) = q-1，域上 q-1 次多项式至多 q-1 个根（`card_roots'`）。
所以恰好 q-1 个根（上界 = 下界）。

**Step 3**: 根的分配。

X^{q-1} - 1 = (X^m - 1)(X^m + 1)，两者互素。
X^{q-1}-1 的每个根 a 满足 a^m = 1 或 a^m = -1（不可能两者都满足，因 1 ≠ -1 in char ≠ 2）。

- 设 A = X^m - 1 的根集（a^m = 1），|A| ≤ m（度数上界）
- 设 B = X^m + 1 的根集（a^m = -1），|B| ≤ m（度数上界）
- A ∩ B = ∅（因 1 ≠ -1）
- A ∪ B = F*（q-1 个元素）

由 |A| + |B| = q - 1 = 2m 且 |A| ≤ m 且 |B| ≤ m：
|A| = m = (q-1)/2，|B| = m = (q-1)/2。 ✓

### Mathlib 路径

- `FiniteField.pow_card_sub_one_eq_one`：a ∈ F* → a^{|F|-1} = 1
- `Polynomial.card_roots'`：roots.card ≤ natDegree
- `sq_sub_one_eq_mul` 或 `ring`：X^{2m} - 1 = (X^m - 1)(X^m + 1)
- `IsCoprime` 的构造：从 "差 = 2 是单位" 推 "互素"
- `Finset.card_union_of_disjoint` 或直接用 `card_le` + `card_ge` 夹逼

### Lean 形式化选择

**方案 A（推荐）**：直接在 Fˣ（乘法群的单位群）上操作，用 Finset.filter。

```lean
-- 主定理
theorem card_pow_half_eq_one ... :
    (Finset.univ.filter (fun a : Fˣ => a ^ m = 1)).card = m := by
  -- 1. 上界：X^m - 1 至多 m 个根（card_roots'）
  -- 2. 下界：complement 也至多 m 个根（a^m = -1，即 a^m + 1 = 0 的根）
  -- 3. 两者之和 = |F*| = q - 1 = 2m
  -- 4. 由夹逼：两者都恰好 m 个
```

**方案 B**：用 Polynomial.roots 直接操作。

```lean
-- X^m - 1 在 F 中的根数 = m
theorem card_roots_X_pow_half_sub_one ... :
    ((X ^ m - 1 : Polynomial F).roots.toFinset).card = m
```

方案 A 更自然（EDF 中用的是元素满足条件的计数）；方案 B 更直接地对应多项式理论。
选方案 A（更贴近应用场景）。

### 技术难点

1. **Fˣ vs F 的转换**：`Fintype.card Fˣ = Fintype.card F - 1`
   Mathlib: `Fintype.card_units`。

2. **char ≠ 2 → 1 ≠ -1**：需要 `CharP.one_ne_neg_one` 或类似引理。
   或者：若 (1 : F) = -1，则 1 + 1 = 0，即 (2 : F) = 0，即 ringChar F | 2，
   与 ringChar F ≠ 2 矛盾。

3. **Nat 算术**：2 * ((q-1)/2) = q - 1 需要 q 奇数的证明。
   q = |F|，F 是有限域，char ≠ 2 → |F| = p^n（p 奇素数）→ |F| 奇数。

4. **Finset 计数的夹逼**：|A| ≤ m, |B| ≤ m, |A| + |B| = 2m → |A| = |B| = m。
   这是简单的 Nat 算术：`omega` 应该可以解决。

## 预计行数

~40-50 行：
- 辅助引理（1 ≠ -1 等）：~10 行
- 根计数上界（card_roots'）：~10 行
- 夹逼论证：~15 行
- 最终组装：~10 行
