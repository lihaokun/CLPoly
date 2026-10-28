# L1 管线框架完成

Date: 2026-05-24

## 做了什么

创建 `Pipeline/L1.lean`，把 L1 翻译代码（Corpus.lean 中的 `_ir` 函数）接入 `factor_Zp_correct` 管线。

## 架构

```
factor_Zp_correct (框架定理: 接受算法实现 + 正确性证明作为参数)
  └── factor_Zp_l1 (实例化: 传入 L1 wrapper)
        ├── sqfZp_l1      = toPolyList(__squarefree_Zp_ir(toSparsePolyZp f))   [sorry]
        ├── ddf_l1        = toPolyList(__ddf_Zp_ir(toSparsePolyZp f))         [sorry]
        └── edf_l1        = 暂用 L2 edf_correct_unconditional                 [已证明]
```

## 关键决策

1. **`toSparsePolyZp` 不可计算** — Polynomial (Finsupp) 本质不可计算，函数的 three `sorry` lemmas (wellFormed, allReduced, toPoly) 待填
2. **EDF 暂用 L2** — `__edf_Zp_ir` 有随机状态参数 `rng: Rng` 和输出缓冲区 `result: Array SparsePolyZp`，接口与 pipeline 期望的 `Polynomial → ℕ → List Polynomial` 差异太大，留到下一阶段处理
3. **`edf_l1_correct` 已完整证明** — 使用 `edf_correct_unconditional` + `WfDvdMonoid.exists_irreducible_factor` 证明 `0 < d` 分支可达

## 涉及文件

- `proof/lean/CLPoly/Pipeline/L1.lean` (新)
- `proof/lean/CLPoly.lean` (添加 import CLPoly.Pipeline.L1)
