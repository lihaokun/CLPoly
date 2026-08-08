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
- 以 `lenC-k` 为良基度量闭合完整外层递归：每轮从当前 raw heap 执行真实
  点积和写入，扩展输出前缀，并通过 allocation 不相交性保持 A/B 的
  `SlicePolyRep` 与规范剩余；退出时 C 的全部 `lenC` 个 cell 均对应乘积
  系数。
- 从完整点态输出不变式构造实际 `SlicePolyRep`：先由 raw slice 的唯一观测
  表示逐系数识别缓冲区内内容，再用卷积 antidiagonal 与 A/B 声明长度证明
  `lenA+lenB-1` 之外的所有乘积系数为零。由此得到 `_classical_mul` 的完整
  raw slice 到 L2 乘积精化定理。
- 强化输出前缀不变式，使每个实际 `_lll_mod_preinv` 结果同时携带 `< p`
  的规范剩余证明，从而得到完整输出 `CanonicalU64Prefix`。
- 在素数模数下，从输入 `RawDensePolyRep` 的真实 normalization 结果推出
  A/B 末系数非零；用卷积最高项公式证明乘积末系数非零，再直接展开实际
  `normaliseU64` 的末 cell 分支，证明输出长度保持 `lenA+lenB-1`。
- 最终定理 `classicalMul_refines` 已将真实 `_classical_mul` 从两个规范 raw
  输入严格精化为 `RawDensePolyRep (left * right)`，包括终止执行、布局、
  规范剩余、完整 L2 表示和规范长度。
- 新增 `_kar_mul` 的完整 RawHeap 生成层：按源码顺序执行半区模加、奇数
  尾项复制、P0/P1/P2 三次递归、两轮原地模减、P0 `copyU64`、间隙清零和
  交叉项原地模加。递归以 `n` 为良基度量，并从 `n≥16` 证明 `n/2` 与
  `n-n/2` 都严格小于 `n`；没有引入执行次数参数。
- 对 Karatsuba 的三类原始循环分别建立终止执行与布局保持定理：半区相加
  循环验证四次实际读取和两次写入，交叉项模减循环验证读-读-写，最终组装
  循环验证偏移 C cell 与 P1 cell 的读-读-写；全部以剩余区间良基递归。
- 对三个 helper 建立精确帧定理：半区相加只可能改写 t1/t2 的剩余同索引
  cells，模减只可能改写目标剩余区间，组装只可能改写 `C[m+i]` 的剩余区间。
  每个定理逐次应用真实 `readU64_writeU64_ne`，可用于证明输入、递归结果和
  已完成输出在后续原地更新中保持不变。
- 对半区相加 helper 的任意当前索引证明最终 raw 值：真实四次 A/B 读取经
  `nmod_add` 写入 t1/t2 后，由帧定理保持到循环结束。进一步从输入
  `SlicePolyRep` 与 `CanonicalU64Prefix` 推出两个输出分别表示对应低/高半
  L2 系数之和，并且实际机器值继续严格小于 p。
- 将当前索引结论提升到完整剩余区间：每轮 t1/t2 写入后，利用 scratch 与
  A/B allocation 不相交性逐 cell 构造 `SameU64Prefix`，重新传递输入
  `SlicePolyRep` 和规范剩余，再对 `m-(i+1)` 做良基递归。由此从 i=0 可
  得到全部 `k<m` 的 t1/t2 raw 值与 L2 半区和语义。

## 当前边界

本步只完成 schoolbook 执行层。Karatsuba 与 `_mul` 的源码身份已经锁定，
但其零填充、scratch 布局、三次递归乘法、交叉项减法和组装尚未 lowering，
因此不能称为 `_mul` 精化完成。下一步先证明 schoolbook 的终止执行、点积
内容和完整输出表示已闭合。Karatsuba 的 raw 生成执行和良基递归现已落下；
下一步证明各 helper 的终止/帧性质、scratch 分区不相交性及 Karatsuba 恒等式。
在这些语义精化和 `_mul` 分派完成前仍不能称为 `_mul` 精化完成。

## 涉及文件

- `proof/lean/CLPoly/Generated/StrictMul.lean`
- `proof/lean/CLPoly/Impl/StrictMulRefinement.lean`
- `proof/cpp2lean_v2/tests/check_strict_mul_source.py`
- `docs/devlog/2026-08-08-strict-mul-lowering.md`
