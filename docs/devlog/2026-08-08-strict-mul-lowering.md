# HGCD 原始乘法 lowering

日期：2026-08-08

## 做了什么

- 固定 C++ `_classical_mul`、`_kar_mul`、`_mul` 三个方法的稳定 AST 哈希。
- 新增 schoolbook 乘法的 RawHeap 执行层。
- 内层严格使用源码的 `j_min..j_max` 点积区间，逐项读取 A/B，执行真实
  `_umul128` 和三字 `_add_carry3`。
- 外层对每个累加结果调用真实 `_lll_mod_preinv`，随后写入 C。
- 内外循环都以未处理区间为终止度量，并显式传播所有 RawHeap 故障。

## 当前边界

本步只完成 schoolbook 执行层。Karatsuba 与 `_mul` 的源码身份已经锁定，
但其零填充、scratch 布局、三次递归乘法、交叉项减法和组装尚未 lowering，
因此不能称为 `_mul` 精化完成。下一步先证明 schoolbook 的终止执行、点积
内容和 L2 乘法表示，再进入 Karatsuba 的良基递归。

## 涉及文件

- `proof/lean/CLPoly/Generated/StrictMul.lean`
- `proof/cpp2lean_v2/tests/check_strict_mul_source.py`
- `docs/devlog/2026-08-08-strict-mul-lowering.md`
