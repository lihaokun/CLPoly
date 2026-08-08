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

## 当前边界

本步闭合了 HGCD 矩阵更新所调用的原始加减执行层，但尚未完成其 RawHeap
输出到 L2 多项式加减的表示定理，因此不能单独称为加减精化完成。下一步将证明
有效容量、别名分支、规范系数和 `SlicePolyRep` 在这些真实执行后成立，再接到
`_mat_row_update`。

## 涉及文件

- `proof/lean/CLPoly/Generated/StrictPolyAddSub.lean`
- `proof/cpp2lean_v2/tests/check_strict_poly_add_sub.py`
- `docs/devlog/2026-08-08-strict-poly-add-sub-lowering.md`
