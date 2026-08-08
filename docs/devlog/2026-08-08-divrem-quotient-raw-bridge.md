# Divrem quotient 规范化 raw 桥

日期：2026-08-08

## 做了什么

- 加强 `quotientLoop_refines_polynomial`，归纳证明实际写入的 quotient 前缀全部是模数下的 canonical residue。
- 证明 remainder 循环保持 quotient 前缀，并将 canonical 性缩短到 `normaliseU64` 返回的实际 `lenQ`。
- 导出 quotient 在实际长度上的 `normaliseU64` fixed point。
- 加强 `polyDivrem_refines`、`polyDivrem_preserves_normalized_gcd` 和 `polyDivrem_next_state`，直接返回 quotient 的完整 `RawDensePolyRep`。

## 为什么做

HGCD 的下一步会把真实 `_poly_divrem` 产生的 quotient 传给 `_mat_row_update`，后者内部 `_mul` 要求规范化且 canonical 的 raw 输入。仅有 `SlicePolyRep` 不足以证明真实 C++ 调用安全和正确，因此必须在 divrem 层闭合这条 raw→safe 表示桥。

## 关键决策

- canonical 性从每轮生成的 `qi` 模乘结果 `< p` 推导，不在 HGCD 层添加 quotient 合法性假设。
- remainder 阶段通过实际写区不相交定理传递 quotient 内容。
- quotient 的最终规范化由源码执行的 `normaliseU64` 及其 fixed-point 定理给出。

## 涉及文件

- `proof/lean/CLPoly/Impl/StrictDivremRefinement.lean`
- `proof/lean/CLPoly/Impl/StrictEuclidRefinement.lean`
- `docs/devlog/2026-08-08-divrem-quotient-raw-bridge.md`
