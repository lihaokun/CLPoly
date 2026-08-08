# HGCD 原始多项式加减 lowering

日期：2026-08-08

## 做了什么

- 新增 `_poly_add` 与 `_poly_sub` 的 RawHeap 执行模型。
- 共同前缀逐系数读取两个输入，并调用与 C++ 一致的溢出安全模加、模减后写入。
- 加法长尾严格保留源代码的地址相等判断和 `memcpy` 分支。
- 减法分别实现正向复制长尾和 `0 - B[i]` 长尾。
- 两个入口均执行真实 `_poly_normalise`，显式返回内存故障或最终逻辑长度。
- 所有循环按未处理长度良基终止，没有使用非良基定义。
- 新增门禁，固定 `nmod_add`、`nmod_sub`、`_poly_add`、`_poly_sub`
  四个 C++ 方法的稳定 AST 哈希。
- 证明三个逐系数循环在有效容量下必然成功，并在每次真实写入后保持整个
  RawHeap 的分配布局。
- 证明 `_poly_add`、`_poly_sub` 完整入口（包括所有别名、复制、取负和
  normalization 分支）必然成功，保持布局且返回长度不超过最大输入长度。
- 证明 C++ 溢出安全 `nmod_add/nmod_sub` 的精确自然数模等式；由输入系数
  的规范范围推出输出仍小于模数，并证明转换到 `ZMod p` 后分别等于真正的
  L2 加法与减法。
- 明确刻画 C++ 所支持的别名条件：输出与输入完全同址，或属于不同分配。
- 证明三个循环只写各自的目标区间，区间外任意地址的读值保持不变。
- 在上述别名条件下证明加减共同前缀的逐系数内容：最终输出严格来自初始
  两个输入系数所执行的真实 `nmod_add/nmod_sub`，包括原位执行。
- 证明减法长尾的 C++ 条件表达式在规范输入上仍产生规范系数，映射到
  `ZMod p` 后精确等于取负，并证明整个长尾循环在原位执行时逐系数成立。

## 当前边界

本步闭合了 HGCD 矩阵更新所调用的原始加减执行层及终止执行桥，但尚未完成其
RawHeap 输出到 L2 多项式加减的表示定理，因此不能单独称为加减精化完成。
下一步将证明别名分支、规范系数和 `SlicePolyRep` 在这些真实执行后成立，再
接到 `_mat_row_update`。

## 涉及文件

- `proof/lean/CLPoly/Generated/StrictPolyAddSub.lean`
- `proof/lean/CLPoly/Impl/StrictPolyAddSubRefinement.lean`
- `proof/cpp2lean_v2/tests/check_strict_poly_add_sub.py`
- `docs/devlog/2026-08-08-strict-poly-add-sub-lowering.md`
