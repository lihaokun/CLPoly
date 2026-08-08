# HGCD 原始乘法 lowering

日期：2026-08-08

## 做了什么

- 固定 C++ `_classical_mul`、`_kar_mul`、`_mul` 三个方法的稳定 AST 哈希。
- 新增 schoolbook 乘法的 RawHeap 执行层。
- 内层严格使用源码的 `j_min..j_max` 点积区间，逐项读取 A/B，执行真实
  `_umul128` 和三字 `_add_carry3`。
- 外层对每个累加结果调用真实 `_lll_mod_preinv`，随后写入 C。
- 内外循环都以未处理区间为终止度量，并显式传播所有 RawHeap 故障。
- 从源码的 `j_min/j_max` 条件表达式证明每个 `A[j]` 与 `B[k-j]` 索引合法。
- 证明内层点积、外层逐系数写入及完整 schoolbook 入口在有效容量和非空输入
  下必然成功，并在所有写入后保持 RawHeap 分配布局。
- 定义与 C++ 点积循环读取完全相同 raw cells 的无界自然数求和执行，并证明
  它在有效索引下必然成功。
- 逐次组合已验证的 `_umul128/_add_carry3` 定理，证明最终三字累加器与
  “初始累加器 + 原始乘积和”在机器宽度 `2^192` 下同余。
- 由每个规范系数至多为 `p-1`，递归证明 raw 点积和不超过
  `项数 × (p-1)^2`。
- 使用项数小于 `2^64`、`p>1` 和既有 lazy-accumulation 容量不等式，证明
  点积和严格小于 `2^192`，从而把机器同余提升为普通自然数精确相等。

## 当前边界

本步只完成 schoolbook 执行层。Karatsuba 与 `_mul` 的源码身份已经锁定，
但其零填充、scratch 布局、三次递归乘法、交叉项减法和组装尚未 lowering，
因此不能称为 `_mul` 精化完成。下一步先证明 schoolbook 的终止执行、点积
内容和 L2 乘法表示，再进入 Karatsuba 的良基递归。当前精确 raw 和还需要
通过 `_lll_mod_preinv` 连接到 `ZMod p` 卷积系数。

## 涉及文件

- `proof/lean/CLPoly/Generated/StrictMul.lean`
- `proof/lean/CLPoly/Impl/StrictMulRefinement.lean`
- `proof/cpp2lean_v2/tests/check_strict_mul_source.py`
- `docs/devlog/2026-08-08-strict-mul-lowering.md`
