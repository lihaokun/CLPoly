# HGCD concrete well-founded recursive call

## Date

2026-08-09

## What changed

- Added `hgcdRecursiveCallChecked`, a concrete terminating execution of the
  generated `_hgcd_recursive` body with measure `lenA`.
- Added an execution bridge proving that, on the generated ordered/decreasing
  child calls, its raw checks erase to the exact branch-admissible source body.
- Added a strong-induction refinement theorem directly for that concrete
  execution.
- Removed the abstract `HgcdRecursiveCallUnfolds` fixed-point premise and the
  old theorem family that could only refine a caller-supplied callback.
- Changed the strict source gate to require the concrete execution/theorem and
  reject reintroduction of the abstract fixed-point route.

## Why

An assumed callback satisfying the recursive source equation is not an
executable L1 semantics.  Even though the earlier semantic induction was
well-founded, leaving the callback abstract would allow a future top-level
proof to claim HGCD refinement without constructing the C++ recursion.

## Key decisions

- A failed recursion-order/decrease check yields only `RawFault.assertionFailure`.
  It never invokes an L2 algorithm or returns a specification value.
- On every physically admissible source invocation, existing first- and
  second-child bounds prove both checks true, so the checked call is exactly
  the generated body on the refinement path.
- The semantic theorem recursively invokes itself at the actual child lengths;
  no fuel or independent recursion counter is present.

## Problems and resolution

The generated source callback type does not carry its decrease proof.  The
checked call restores that proof at the raw/safe boundary, while the execution
bridge eliminates the checks using the proof-indexed generated dispatch.

## Files

- `proof/lean/CLPoly/Generated/StrictHGCD.lean`
- `proof/lean/CLPoly/Generated/StrictHGCDChecked.lean`
- `proof/lean/CLPoly/Impl/StrictHGCDRawRefinement.lean`
- `proof/lean/CLPoly/Impl/StrictHGCDFirstCallRefinement.lean`
- `proof/lean/CLPoly/Impl/StrictHGCDReconstructRefinement.lean`
- `proof/lean/CLPoly/Impl/StrictHGCDRecursiveRefinement.lean`
- `proof/lean/CLPoly/Impl/StrictHGCDContinuationRefinement.lean`
- `proof/lean/CLPoly/Impl/StrictHGCDCheckedRefinement.lean`
- `proof/cpp2lean_v2/tests/check_strict_hgcd_source.py`
- `docs/devlog/2026-08-09-hgcd-physical-continuation-packages.md`
- `docs/devlog/2026-08-09-hgcd-concrete-well-founded-call.md`

## 度量

- 耗时：约 3 小时（递归执行设计、等价桥、强归纳迁移、残留清理）
- 迭代：4 轮编译—修复循环
- Lean 新增/修改行数：约 140 行净变化；并将原 11k 行单模块按证明
  依赖拆分为可稳定生成 `.olean` 的六个模块
- 对应 C++ 行数：约 250 行（`_hgcd_recursive` 完整控制流）
- 放弃的方案：继续要求任意 `plain` 回调满足固定点方程；该方案没有
  构造可执行 L1 递归，因此从严格路径中完全移除。
