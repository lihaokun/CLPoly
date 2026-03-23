# Recombination 算法验证：从 UFD 捷径到 Hensel 因子恢复

> 状态：nl-proof v1
> 目标：加强 recombine_correct，使用 Hensel 因子恢复而非 UFD 存在性

---

## 0. 当前差距

当前 `recombine_correct` 用 `WfDvdMonoid.exists_factors`（UFD），完全不涉及 Hensel 因子。
C++ 用 Zassenhaus 子集枚举 + trial division + 从 Hensel 因子恢复真因子。

---

## 1. 方案

不建模 Zassenhaus 子集枚举（这是 C++ 控制流细节）。
而是证明**因子恢复定理**：在 Hensel 提升不变量 + Mignotte 精度下，真因子可从 Hensel 因子恢复。

**新定理**：
```lean
theorem recombine_from_hensel
    (f : Polynomial ℤ) (hf : f ≠ 0) (hprim : IsPrimitive f)
    (p : ℕ) (hp : Nat.Prime p) (k : ℕ) (hk : 0 < k)
    -- Hensel 因子
    (H : List (Polynomial ℤ))
    -- (I1) 乘积还原 mod p^k
    (hprod : Polynomial.map (Int.castRingHom (ZMod (p^k))) f =
             (H.map (Polynomial.map (Int.castRingHom (ZMod (p^k))))).prod)
    -- (I4) mod p 一致
    (h_mod_p : ∀ i, ...)
    -- 互素 mod p
    (hcop : ...)
    -- Mignotte 精度
    (hmig : ...)
    : ∃ result, RecombineCorrect f result
```

**但**这签名很复杂。简化方案：**保持现有 `recombine_correct`（UFD 版），另外新增 `factor_recovery` 引理证明 Hensel 因子可恢复真因子。**

---

## 2. factor_recovery 引理

### 2.1 定理

设 f ∈ ℤ[x] primitive squarefree，f = g₁·...·gₜ（不可约分解）。
设 H₁,...,Hᵣ 是 Hensel 因子（mod p^k），满足 (I1)(I4) + 互素 + Mignotte 精度。

则**每个真因子 gⱼ 可从 Hensel 因子恢复**：
存在子集 Sⱼ ⊆ {1,...,r}，使得 pp(sym_m(P_{Sⱼ})) = gⱼ。

### 2.2 证明骨架（已在 phase4-recombine.md §3 详述）

1. 真因子 gⱼ 和 Hensel 因子 Hᵢ 通过 mod p 因子对应
2. Hensel 唯一性 → ∏_{i∈Sⱼ} Hᵢ ≡ (ℓ/ℓⱼ)·gⱼ (mod m)
3. Mignotte → sym_m 精确恢复
4. pp 提取 gⱼ

### 2.3 Lean 可行性

这个引理需要大量参数和前置条件。形式化工作量约 ~200 行。
主要依赖已有的：
- `hensel_unique`
- `mignotte_bound_l2`
- Z[x] UFD（`WfDvdMonoid.exists_factors`）

---

## 3. 务实方案

与其形式化完整的 factor_recovery（~200 行，大量参数），不如：

1. **保持 `recombine_correct`（UFD 版）不变**——它数学正确，满足 RecombineCorrect spec
2. **新增 `factor_recovery_theorem`** 作为独立数学定理——证明 Hensel 因子可恢复，不用于 Pipeline
3. 这样 Pipeline 继续 0 sorry，factor_recovery 是额外的算法验证

### 3.1 factor_recovery_theorem 签名

```lean
/-- 因子恢复：在 Hensel 提升不变量 + Mignotte 精度下，
    每个不可约因子可从 Hensel 因子的子集乘积恢复。
    这验证了 C++ __zassenhaus_recombine 的数学正确性。-/
theorem factor_recovery_theorem
    (f : Polynomial ℤ) (hf : f ≠ 0)
    (p : ℕ) (hp : Nat.Prime p) (m : ℕ) (hm : 1 < m)
    -- f 在 ℤ[x] 中有不可约因子 g
    (g : Polynomial ℤ) (hg_irred : Irreducible g) (hg_dvd : g ∣ f)
    -- Hensel 因子 A, B with f ≡ A*B (mod m)，B monic
    (A B : Polynomial ℤ)
    (hprod_m : Polynomial.map (Int.castRingHom (ZMod m)) f =
               Polynomial.map (Int.castRingHom (ZMod m)) A *
               Polynomial.map (Int.castRingHom (ZMod m)) B)
    (hB_monic : Monic B)
    -- mod p: A ≡ (对应 g 的因子), B ≡ (其余)
    (hA_mod : Polynomial.map (Int.castRingHom (ZMod p)) A =
              Polynomial.map (Int.castRingHom (ZMod p)) (C (leadingCoeff f / leadingCoeff g) * g))
    -- IsCoprime mod p
    (hcop : IsCoprime (Polynomial.map (Int.castRingHom (ZMod p)) A)
                       (Polynomial.map (Int.castRingHom (ZMod p)) B))
    -- Mignotte 精度
    (hmig : ∀ i, ((leadingCoeff f / leadingCoeff g * g).coeff i).natAbs < m / 2)
    -- Hensel 唯一性给出
    : Polynomial.map (Int.castRingHom (ZMod m)) A =
      Polynomial.map (Int.castRingHom (ZMod m)) (C (leadingCoeff f / leadingCoeff g) * g) := by
  -- 由 hensel_unique：两组因子 (A,B) 和 ((ℓ/ℓⱼ)·g, (ℓⱼ/ℓ)·h*) 在 mod m 下相同
  -- 因为 mod p 一致 + B monic + 互素
  sorry
```

这就是 Hensel 唯一性的直接应用。

---

## 4. 结论

**factor_recovery_theorem** 是 hensel_unique 的包装——在具体的因式分解场景中应用唯一性。

实际上，我们已有 `hensel_unique`，它就是 factor_recovery 的核心。唯一需要补的是：
- 构造对比组 `(C(ℓ/ℓⱼ)·gⱼ, (ℓⱼ/ℓ)·hⱼ*)` 满足唯一性前置条件
- 验证 mod p 一致、B monic、互素

这大部分是参数组装，不是新数学。

**建议**：不额外写 factor_recovery_theorem。当前 `recombine_correct`（UFD）+ `hensel_unique` + `mignotte_bound_l2` 三者**组合**已经构成了 C++ 算法正确性的数学基础。差距只是没有显式组装成一个定理。

**实际影响**：对 CLAUDE.md 差距表，Recombination 项标注为 "数学基础已有（hensel_unique + mignotte），显式组装待定"。
