# Bridge lemmas: toSparsePolyZp 三个前提条件

## 做了什么

补全了 `Pipeline/L1.lean` 中 `toSparsePolyZp` 的三个 bridge lemma：

1. **`toSparsePolyZp_wellFormed`** — `SparsePolyZp.WellFormed p (toSparsePolyZp f)`
   - 每个 `Zp` 系数的 `prime` 字段 = `UInt64.ofNat p`，因此 `.prime.toNat = p`

2. **`toSparsePolyZp_allReduced`** — `SparsePolyZp.AllReduced p (toSparsePolyZp f).toList`
   - 每个 `Zp` 系数的 `val` = `UInt64.ofNat (ZMod.val (f.coeff n))`，且 `ZMod.val (f.coeff n) < p`

3. **`toSparsePolyZp_toPoly`** — `SparsePolyZp.toPoly p (toSparsePolyZp f) = f`（roundtrip）
   - 通过 `Finset.sort` → `listSum` → `Polynomial.as_sum_support` 链证明
   - 关键引理：`listSum` 对 sorted support 求和 = `∑ n ∈ f.support, monomial n (f.coeff n)`
   - 使用了 `Finset.sort_nodup`、`List.mem_toFinset`、`Polynomial.as_sum_support`

## 为什么做

三个 bridge lemma 是所有精化定理使用的前提条件。精化定理要求 `SparsePolyZp.WellFormed` 和 `AllReduced`，而 `sqfZp_l1_correct`/`ddf_l1_correct` 的证明需要 `toPoly` roundtrip。

## 关键决策

- 增加 `hp_size : 2 * p ≤ UInt64.size` 参数（与精化定理一致），用于保证 `UInt64.ofNat` 不溢出
- 使用 `nlinarith` 从 `2*p ≤ UInt64.size` 推导 `p < UInt64.size`
- 使用 `calc` 而非 `simp` 处理 `hlistSum` 的归纳步骤（避免 `simp` 无法处理 `Finset` 成员关系的困境）

## 涉及文件

- `proof/lean/CLPoly/Pipeline/L1.lean` — 新增 3 个 bridge lemma，添加 import
