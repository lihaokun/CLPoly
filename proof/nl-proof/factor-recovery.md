# Factor Recovery：从 Hensel 提升因子恢复真因子

> 状态：nl-proof v1
> 对应 C++：`polynomial_factorize_univar.hh:750-882` `__zassenhaus_recombine`
> 依赖：`hensel_unique`（已形式化 0 sorry），`mignotte_bound_l2`（已形式化 0 sorry）

---

## 0. 目标

证明 C++ `__zassenhaus_recombine` 的核心数学正确性：

在 Hensel 提升不变量 + Mignotte 精度下，Hensel 因子子集乘积的对称约化的原始部分 = 真因子。

**数学链**：
1. `hensel_unique` → A ≡ c·g (mod m)
2. `symmetric_recovery` → A = c·g（精确，在 ℤ[x]）
3. `primPart` → pp(A) = g

本 nl-proof 证明步骤 2（对称恢复引理），步骤 1 已有 `hensel_unique`，步骤 3 是 Gauss 引理的直接应用。

---

## 1. 整数对称恢复引理

### 1.1 定理

```
若 a ≡ b (mod m)，|a| · 2 < m，|b| · 2 < m，则 a = b。
```

### 1.2 证明

m | (a - b)。|a - b| ≤ |a| + |b| < m/2 + m/2 ≤ m。
故 |a - b| < m 且 m | (a - b) → a - b = 0。✓

### 1.3 Lean 形式化

```lean
lemma int_eq_of_mod_eq_of_abs_lt (a b : ℤ) (m : ℕ) (hm : 0 < m)
    (hmod : (a : ZMod m) = (b : ZMod m))
    (ha : a.natAbs * 2 < m) (hb : b.natAbs * 2 < m)
    : a = b := by
  suffices h : a - b = 0 by linarith [h]  -- 或 omega
  have hdvd : (m : ℤ) ∣ (a - b) := by
    rw [← ZMod.intCast_zmod_eq_zero_iff_dvd]
    push_cast; rw [sub_eq_zero]; exact hmod
  -- |a - b| < m
  have habs : (a - b).natAbs < m := by
    calc (a - b).natAbs ≤ a.natAbs + b.natAbs := Int.natAbs_sub_le a b
      _ < m := by omega  -- ha + hb: 2*|a| < m, 2*|b| < m → |a|+|b| < m
  -- m | (a-b) 且 |a-b| < m → a-b = 0
  obtain ⟨k, hk⟩ := hdvd
  by_contra hne
  have hk_ne : k ≠ 0 := by rintro rfl; simp at hne
  have : (a - b).natAbs = m * k.natAbs := by
    rw [hk, Int.natAbs_mul, Int.natAbs_natCast]
  have : 0 < k.natAbs := Int.natAbs_pos.mpr hk_ne
  omega  -- m * k.natAbs ≥ m > (a-b).natAbs，矛盾
```

**Lean 路径**：
- `ZMod.intCast_zmod_eq_zero_iff_dvd`：(a : ZMod m) = 0 ↔ m | a
- `Int.natAbs_sub_le`：|a - b| ≤ |a| + |b|
- `Int.natAbs_mul`：|a * b| = |a| * |b|
- `Int.natAbs_natCast`：|(n : ℤ)| = n
- `Int.natAbs_pos`：|k| > 0 ↔ k ≠ 0

---

## 2. 多项式对称恢复引理

### 2.1 定理

```
若 map_m(P) = map_m(Q)，且 ∀ i, |P.coeff i| · 2 < m，∀ i, |Q.coeff i| · 2 < m，
则 P = Q。
```

### 2.2 证明

逐系数应用 §1：

对任意 i，P.coeff i ≡ Q.coeff i (mod m)（从 map_m(P) = map_m(Q) 的系数相等）。
|P.coeff i| · 2 < m，|Q.coeff i| · 2 < m。
由 §1 的 `int_eq_of_mod_eq_of_abs_lt`：P.coeff i = Q.coeff i。

由多项式的 ext 引理：P = Q。✓

### 2.3 Lean 形式化

```lean
theorem symmetric_recovery (P Q : Polynomial ℤ) (m : ℕ) (hm : 0 < m)
    (hmod : Polynomial.map (Int.castRingHom (ZMod m)) P =
            Polynomial.map (Int.castRingHom (ZMod m)) Q)
    (hP : ∀ i, (P.coeff i).natAbs * 2 < m)
    (hQ : ∀ i, (Q.coeff i).natAbs * 2 < m)
    : P = Q := by
  ext i
  apply int_eq_of_mod_eq_of_abs_lt _ _ m hm
  · -- (P.coeff i : ZMod m) = (Q.coeff i : ZMod m)
    have h := congr_arg (fun p => p.coeff i) hmod
    simp only [Polynomial.coeff_map] at h
    exact_mod_cast h
  · exact hP i
  · exact hQ i
```

**Lean 路径**：
- `Polynomial.ext`：∀ i, p.coeff i = q.coeff i → p = q
- `Polynomial.coeff_map`：(map f p).coeff i = f (p.coeff i)

---

## 3. Factor Recovery 定理

### 3.1 定理

```lean
/-- 因子恢复：若 A ≡ c·g (mod m) 且两者系数都 < m/2，则 A = c·g（精确）。
    这验证了 C++ __zassenhaus_recombine 的核心数学步骤：
    symmetric_mod(subset_product) 在 Mignotte 精度下精确恢复缩放后的真因子。-/
theorem factor_recovery
    (g A : Polynomial ℤ) (m : ℕ) (hm : 0 < m) (c : ℤ)
    -- A ≡ c·g (mod m) [由 hensel_unique 提供]
    (hmod : Polynomial.map (Int.castRingHom (ZMod m)) A =
            Polynomial.map (Int.castRingHom (ZMod m)) (C c * g))
    -- A 有小系数（C++ 中 __symmetric_mod 保证）
    (hA_small : ∀ i, (A.coeff i).natAbs * 2 < m)
    -- c·g 有小系数（Mignotte bound 保证）
    (hcg_small : ∀ i, ((C c * g).coeff i).natAbs * 2 < m)
    : A = C c * g
```

### 3.2 证明

直接应用 `symmetric_recovery`。✓

### 3.3 与 C++ 算法的对应

C++ `__zassenhaus_recombine` 的核心循环：

```cpp
// line 806-859
auto trial = __symmetric_mod(__subset_product_mod(subset, lifted, m), m);
auto candidate = __upoly_primitive(trial);
if (f_star % candidate == 0) {
    result.push_back(candidate);
    f_star = f_star / candidate;
}
```

对应数学链：
1. `__subset_product_mod(S, lifted, m)` = ∏_{i∈S} Hᵢ mod m
2. `__symmetric_mod(_, m)` = sym_m(∏_{i∈S} Hᵢ) — 系数约化到 (-m/2, m/2]
3. `hensel_unique` → sym_m(∏_{i∈S} Hᵢ) ≡ (ℓ/ℓⱼ)·gⱼ (mod m)
4. `factor_recovery`（本定理）→ sym_m(∏_{i∈S} Hᵢ) = (ℓ/ℓⱼ)·gⱼ（精确）
5. `__upoly_primitive` = pp → pp((ℓ/ℓⱼ)·gⱼ) = gⱼ（因 gⱼ primitive）
6. Trial division 成功因为 gⱼ 确实整除 f

**步骤 3 的前置条件满足**：
- Hensel 因子满足 (I1)-(I4)
- 正确子集 S 存在（因子对应定理，来自 mod p 分解的唯一性）
- Mignotte 精度：m > 2ℓ · C(n,n/2) · ‖f‖₂ 保证 hcg_small
- C++ 的 `__symmetric_mod` 保证 hA_small

**步骤 5 的 primPart**：
若 A = C(c) * g，g primitive，c ≠ 0：
- content(A) = |c| · content(g) = |c|（Gauss 引理：content 是乘法的）
- primPart(A) = A / C(content(A)) = C(c/|c|) * g = ±g
- 在 Lean 中：`Associated (primPart A) g`

---

## 4. Mignotte 精度条件的验证

### 4.1 hcg_small 的满足

设 Q = C(ℓ/ℓⱼ) * g，g | f，g irreducible。

Q.coeff i = (ℓ/ℓⱼ) * g.coeff i。

|Q.coeff i| = |ℓ/ℓⱼ| * |g.coeff i| ≤ ℓ * |g.coeff i|（ℓ/ℓⱼ ≤ ℓ 因 ℓⱼ ≥ 1）。

由 Mignotte bound（`mignotte_bound_l2`）：
|g.coeff i| ≤ C(n, n/2) · ‖f‖₂

故 |Q.coeff i| ≤ ℓ · C(n, n/2) · ‖f‖₂。

精度条件 (M)：m > 2 · ℓ · C(n, n/2) · ‖f‖₂。

故 |Q.coeff i| · 2 ≤ 2ℓ · C(n, n/2) · ‖f‖₂ < m。✓

### 4.2 hA_small 的满足

C++ 的 `__symmetric_mod(poly, m)` 将每个系数映射到 (-m/2, m/2]。
所以 |A.coeff i| ≤ m/2，即 |A.coeff i| · 2 ≤ m。

（注：严格不等式 < 需要排除 |A.coeff i| = m/2 的情况。
当 m 为奇数时 m/2 不是整数所以 < 自动成立。
当 m = p^k（p 为奇素数）时 m 为奇数。
p = 2 时需要单独处理。实际上 C++ 选择 p 为奇素数。）

---

## 5. hensel_unique 的适用性

### 5.1 应用 hensel_unique 得 hmod

要得到 `hmod : map_m(A) = map_m(C(ℓ/ℓⱼ) * g)`，需要构造"对比对"满足 hensel_unique 的前置条件。

**对比对**：
- 第一对 (A, B)：Hensel 提升的结果。B monic in ℤ[x]（由 hensel_step_with_degree H6 保持）。
- 第二对 (A', B')：A' = C(ℓ/ℓⱼ) * g，B' 需要 monic in ℤ[x]。

**B' 的构造**：

设 h = f / g ∈ ℤ[x]。lc(h) = ℓ/ℓⱼ。

B' 不能简单取 h/(ℓ/ℓⱼ)（不一定在 ℤ[x] 中）。

**方案**：对 (g, h) 应用 hensel_two_factor 得到 Hensel 提升 (G', H')。
- G', H' ∈ ℤ[x]，G'*H' ≡ f (mod p^k)
- G' ≡ g (mod p)，H' ≡ h (mod p)
- H' 的 lc 由 hensel_step_with_degree H6 保持：lc(H') = lc(h) = ℓ/ℓⱼ

但 H' 不 monic（lc = ℓ/ℓⱼ ≠ 1），所以不能直接用 hensel_unique。

**需要的变体**：hensel_unique 当前要求 B monic in ℤ[x]。实际上，证明中用到 B monic 的关键点是：
- natDegree(B₁) = natDegree(B₂)（从 monic 推出）
- lc(B₁) = lc(B₂) = 1（用于推出 c = 0）

如果改为 "lc(B₁) = lc(B₂) = ℓ/ℓⱼ"（非零常数），证明仍然成立：
- 1 + p^k · lc(G) 变为 ℓ/ℓⱼ + p^k · lc(G)
- B₂ 也有 lc = ℓ/ℓⱼ → (ℓ/ℓⱼ) + p^k · lc(G) = (ℓ/ℓⱼ) → p | lc(G)

**结论**：hensel_unique 可以泛化为 "B₁, B₂ 有相同的 leadingCoeff"（不要求 = 1），证明完全相同。
但当前已有的形式化用了 `Monic B1`。

### 5.2 务实方案

**方案 A**（推荐）：factor_recovery 以 `hmod` 为假设，不在内部推导。理由：
1. hensel_unique 的泛化/适配是独立的工程工作
2. `hmod` 在 C++ 语境下成立（由 Hensel 唯一性保证）
3. 保持模块化：恢复引理和唯一性引理独立

**方案 B**：泛化 hensel_unique 支持 "相同 lc" 而非 "monic"。约 ~20 行改动。可后续做。

**选择方案 A**。factor_recovery 的 `hmod` 假设由使用者提供（来自 hensel_unique 的泛化版或特化应用）。

---

## 6. 形式化估计

| 改动 | 行数 |
|------|------|
| `int_eq_of_mod_eq_of_abs_lt`（整数版） | ~15 |
| `symmetric_recovery`（多项式版） | ~15 |
| `factor_recovery`（组合定理） | ~5 |
| **总计** | **~35** |

**依赖**（已形式化）：
- `hensel_unique`（Recombine.lean, ~130 行）
- `mignotte_bound_l2`（Recombine.lean, ~50 行）

**后续 TODO**：
- 泛化 hensel_unique 支持 "相同 lc"（~20 行改动）
- primPart 恢复引理（`Associated (primPart (C c * g)) g` for primitive g）
- 因子对应定理（mod p 不可约分解 ↔ ℤ[x] 不可约因子的子集划分）
