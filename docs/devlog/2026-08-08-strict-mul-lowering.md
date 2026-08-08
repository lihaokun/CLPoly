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
- 由精确点积容量推出三字高 limb 小于 p，直接应用已验证的生成
  `_lll_mod_preinv` 定理，证明实际归约输出等于 raw 点积和 `% p`。
- 进一步证明该输出转换到 `ZMod p` 后等于无界 raw 点积和的自然数转换。
- 定义与源码循环边界一致的 L2 系数递归和，并从 `SlicePolyRep` 的实际
  raw reads 逐项证明：无界 raw 点积和转换到 `ZMod p` 后，恰好等于同一
  区间上的 L2 系数乘积和。该桥没有调用 L2 乘法来代替 C++ 执行。
- 证明递归系数和等于 `j_min..j_max` 闭区间上的有限和；再用输入的
  `SlicePolyRep` 证明该区间外的每个标准卷积项均为零，最终得到源码点积
  区间精确等于 `Polynomial.coeff (left * right) k`。
- 对实际外层写循环证明帧性质：它只会改写 `C[k..lenC)`；任意不与这些
  cell 重合的 raw read 在完整剩余执行后保持不变。该性质将用于保持 A/B
  输入表示以及此前已经写出的 C 系数。
- 将单个源码系数的整条执行链闭合：实际 raw 点积、精确三 limb 累加、实际
  `_lll_mod_preinv`、`ZMod` 转换以及标准乘积系数全部相等。
- 从外层帧性质导出与 C 不同 allocation 的整个输入 prefix 保持不变，并
  证明 `CanonicalU64Prefix` 可沿相同 prefix 传递，为外层递归保持 A/B 的
  表示与规范剩余不变式。
- 建立外层循环的已完成输出前缀不变式：一次真实 `writeU64` 将已证明前缀
  精确扩展一个乘积系数，而剩余外层循环保持此前所有输出系数不变。

## 当前边界

本步只完成 schoolbook 执行层。Karatsuba 与 `_mul` 的源码身份已经锁定，
但其零填充、scratch 布局、三次递归乘法、交叉项减法和组装尚未 lowering，
因此不能称为 `_mul` 精化完成。下一步先证明 schoolbook 的终止执行、点积
内容和完整输出表示，再进入 Karatsuba 的良基递归。当前实际 raw 点积、模
归约、L2 区间和及标准乘法卷积系数已全部连通；还需要沿外层 raw writes
证明每个输出 cell 表示该系数，并由输入首项非零证明乘积长度规范化。

## 涉及文件

- `proof/lean/CLPoly/Generated/StrictMul.lean`
- `proof/lean/CLPoly/Impl/StrictMulRefinement.lean`
- `proof/cpp2lean_v2/tests/check_strict_mul_source.py`
- `docs/devlog/2026-08-08-strict-mul-lowering.md`
