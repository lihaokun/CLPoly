# 严格 Hensel 因子提取精化

日期：2026-08-11

## 做了什么

- 为 C++ `__hensel_extract_factors` 建立严格 L1 raw IR，保持先左子树、后右子树的执行顺序。
- 用 Hensel 树节点数作为良基度量，证明满足数组索引不变量时 raw IR 必然成功终止。
- 建立 `HenselExtractCorrect` 语义关系，并证明 raw 执行结果满足该关系。
- 将最终公开精化合同登记到 Pass 9，由生成器统一输出到 `CLPoly.Refinement.Generated`。

## 为什么做

完整 `__hensel_lift` 精化不仅需要提升循环，还必须覆盖提升后从树中提取因子的真实 C++ 路径。公开合同集中生成，也使已完成的原始 C++ 函数入口可直接审计。

## 关键决策

- L1 可执行定义保留在 `CLPoly.Generated.StrictHensel`，证明正文保留在 `CLPoly.Refinement.Hensel`，最终公开合同集中在 `CLPoly.Refinement.Generated`。这样既集中入口，又不混淆被证明程序与证明本身。
- 不使用 fuel、`partial def`、规格 oracle 或 L2 fallback；递归终止性完全由有限树的 `nodeCount` 给出。
- 越界仍是 raw IR 的显式错误；成功性由与 C++ 前置条件对应的 `HenselExtractInvariant` 排除，而不是修改程序语义。

## 遇到的问题与解决方式

Lean 的 equation compiler 对带嵌套 `Option` 的递归不变量生成不稳定，因此改用显式 `HenselLiftTree.rec` 定义同一结构归纳谓词，未改变算法或证明目标。

## 涉及文件

- `proof/lean/CLPoly/Generated/StrictHensel.lean`
- `proof/lean/CLPoly/Refinement/Hensel.lean`
- `proof/lean/CLPoly/Refinement/Generated.lean`
- `proof/cpp2lean_v2/class_map.py`
- `proof/cpp2lean_v2/passes/pass9_refine_gen.py`
