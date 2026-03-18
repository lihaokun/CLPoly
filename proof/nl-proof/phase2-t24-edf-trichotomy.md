# T2.4: EDF 三分性 — Euler 判据的多项式版本

> 状态：已形式化，0 sorry

---

## 定理陈述

设 p 为奇素数，q ∈ F_p[x] 不可约，deg(q) | d，a ∈ F_p[x]。设 m = (p^d - 1) / 2。则：

q | a  ∨  q | (a^m - 1)  ∨  q | (a^m + 1)

```lean
theorem edf_trichotomy
    (hp_odd : p ≠ 2)
    (q : Polynomial (ZMod p)) (hq : Irreducible q)
    (d : ℕ) (hd : 0 < d) (hdvd : q.natDegree ∣ d)
    (a : Polynomial (ZMod p)) :
    let m := (p ^ d - 1) / 2
    q ∣ a ∨ q ∣ (a ^ m - 1) ∨ q ∣ (a ^ m + 1)
```

### 用途

EDF (Cantor-Zassenhaus) 算法的核心数学依据：对等度多项式 f（所有不可约因子度 = d），
取随机 a，计算 h = a^{(p^d-1)/2} mod f。则 f 的每个不可约因子 q 要么整除 a
（概率极低），要么整除 h-1 或 h+1。因此 gcd(h-1, f) 把 f 的不可约因子分成
"h ≡ 1" 和 "h ≡ -1" 两组，实现非平凡分裂。

## 辅助引理：提取 T2.2' 的 "hall" 步骤

T2.2' 的证明中，Step A+B 建立了 "∀ x ∈ AdjoinRoot g, x^{p^d} = x"。
T2.4 需要完全相同的结论（从 deg q | d 出发，先用 T2.2 得 q | X^{p^d}-X，
再推出 ∀ x, x^{p^d} = x）。建议提取为独立引理：

```lean
lemma all_pow_eq_of_dvd_X_pow_sub_X
    (g : Polynomial (ZMod p)) (hg : Irreducible g)
    (d : ℕ) (h : g ∣ (X ^ (p ^ d) - X : Polynomial (ZMod p))) :
    ∀ (x : AdjoinRoot g), x ^ (p ^ d) = x
```

**证明**：与 T2.2' 的 Step A+B 完全相同：

1. Step A: mk g (X^{p^d} - X) = 0 → root^{p^d} = root
2. Step B: 构造 iterateFrobenius 为 AlgHom，由 algHom_ext（在 root 上一致）得 φ = id
3. 因此 ∀ x, x^{p^d} = x

**重构影响**：提取后，T2.2' 可简化为：

```
Step A+B: exact all_pow_eq_of_dvd_X_pow_sub_X g hg d h  -- 一行替代原来 ~20 行
Step C: (不变)
```

需要的 Lean 实例 boilerplate（Fact, CharP, Fintype）移入引理内部。

## T2.4 证明

### Step 1: ∀ x ∈ AdjoinRoot q, x^{p^d} = x

由 T2.2 (`irreducible_dvd_X_pow_sub_X`)：deg(q) | d → q | X^{p^d} - X。
再由辅助引理 `all_pow_eq_of_dvd_X_pow_sub_X`：∀ x, x^{p^d} = x。

### Step 2: 域中的 Euler 三分

设 K = AdjoinRoot q，ā = mk q a。由 Step 1：ā^{p^d} = ā。

- 若 ā = 0：得 q | a（mk_eq_zero 反向）。
- 若 ā ≠ 0：
  - ā^{p^d} = ā → ā^{p^d - 1} = 1（两边除以 ā）
  - p 为奇素数，p^d 为奇数，p^d - 1 为正偶数
  - 设 m = (p^d - 1) / 2，则 2m = p^d - 1
  - ā^{2m} = 1，即 (ā^m)² = 1
  - 域中 y² = 1 → y = 1 ∨ y = -1（`sq_eq_one_iff_of_ne_neg_one` 或 `mul_self_eq_one_iff`）
  - 所以 ā^m = 1 ∨ ā^m = -1

### Step 3: 翻译回多项式整除

- ā = 0       ↔ mk q a = 0         ↔ q | a
- ā^m = 1     ↔ mk q (a^m - 1) = 0 ↔ q | (a^m - 1)
- ā^m = -1    ↔ mk q (a^m + 1) = 0 ↔ q | (a^m + 1)

翻译用 `AdjoinRoot.mk_eq_zero` + `map_sub`/`map_pow`/`map_add`/`map_one`。

### Mathlib 路径

- `irreducible_dvd_X_pow_sub_X`：T2.2（已证）
- `all_pow_eq_of_dvd_X_pow_sub_X`：辅助引理（新建，提取自 T2.2'）
- `AdjoinRoot.mk_eq_zero`：多项式整除 ↔ 商环中为零
- `mul_left_cancel₀`：域中 ā ≠ 0 → ā · x = ā · y → x = y（用于从 ā^{p^d} = ā 推 ā^{p^d-1} = 1）
- `sq_eq_one_iff_of_ne_neg_one` 或直接用 `mul_self_eq_one_iff`：y² = 1 → y = ±1
- `Nat.two_mul_div_two_of_even` 或 `Nat.div_add_mod`：2 * ((p^d-1)/2) = p^d - 1

### 技术难点

1. **Nat 算术**：(p^d - 1) / 2 * 2 = p^d - 1 需要 p^d 为奇数的证明。
   p 为奇素数 + d > 0 → p^d 奇数 → p^d - 1 偶数。
   Mathlib: `Nat.Prime.odd_of_ne_two` + `Odd.pow` → `Even.sub`。

2. **除以非零元素**：在域 AdjoinRoot q 中，ā ≠ 0 时从 ā^{p^d} = ā 推
   ā^{p^d-1} = 1。可以用 `mul_left_cancel₀` 或 `div_eq_one`。
   或者直接：ā^{p^d} - ā = 0 → ā · (ā^{p^d-1} - 1) = 0 → 域中 mul_eq_zero → 两种情况。

3. **y² = 1 在域中**：`sq_eq_one_iff_of_ne_neg_one` 需要 char ≠ 2（即 1 ≠ -1）。
   在 AdjoinRoot q 中，CharP K p，p ≠ 2，所以 char K ≠ 2，所以 (1 : K) ≠ -1。

## 预计行数

- `all_pow_eq_of_dvd_X_pow_sub_X` 引理：~25 行（从 T2.2' 提取）
- T2.2' 简化后：~35 行（原 ~55 行，Step A+B 缩减为 1 行调用）
- `edf_trichotomy`：~25 行
- 总计新增：~50 行，T2.2' 净减少 ~20 行
