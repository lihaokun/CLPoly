# 严格 HGCD 矩阵不变量

日期：2026-08-08

## 做了什么

- 新增 `StrictHGCDRefinement.lean`，刻画两对多项式具有完全相同共同因子的关系。
- 证明双向二阶线性变换保持共同因子，并由此保持规范化 Euclidean GCD。
- 分别证明行列式为 `1` 和 `-1` 的 HGCD 矩阵变换保持规范化 GCD。
- 新增源码门禁，固定 C++ `_hgcd_iter`、`_hgcd_recursive`、`_gcd_hgcd`
  与 `_gcd_euclid` 的稳定 AST 哈希及真实调用结构。
- 验证严格边界、divrem、GCD 源码门禁和 Lean 构建全部通过。

## 当前边界

这一步建立的是实际 HGCD 矩阵操作必须满足的语义不变量，尚未证明 C++
原始内存中的四个矩阵多项式确实实现这些变换，因此不代表 HGCD 精化完成。
下一步需要为 `_mat_one`、`_mat_row_update`、`_mat_mul` 及其调用的原始
多项式运算建立 RawHeap 表示和逐操作精化，再接入 HGCD 的良基递归控制流。

## 生成器探测

当前转换器能够解析并降低四个 HGCD 方法的 MIR，但直接生成结果仍含未证明
占位与非良基定义：`_hgcd_iter` 为 11/2，`_hgcd_recursive` 为 30/6，
`_gcd_hgcd` 为 3/4（前者是未证明占位数，后者是非良基定义数）。主要缺口是
二级指针解引用、`sizeof`、`memcpy`、指针交换、矩阵原始操作及递归循环。
因此这些输出没有进入严格生成边界，也没有被当作证明成果。

## 涉及文件

- `proof/lean/CLPoly/Impl/StrictHGCDRefinement.lean`
- `proof/cpp2lean_v2/tests/check_strict_hgcd_source.py`
- `docs/devlog/2026-08-08-strict-hgcd-invariant.md`
