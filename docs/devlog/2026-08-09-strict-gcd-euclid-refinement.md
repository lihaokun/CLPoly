# Strict Euclid GCD raw refinement

## Date

2026-08-09

## What changed

- Proved that raw `copyU64` transports and frames normalized dense polynomial
  representations across distinct C++ allocations.
- Proved by the same well-founded divisor-length recursion that the generated
  Euclid loop returns one of its three rotating buffers and stays within the
  physical workspace capacity.
- Added the physical allocation contract for `_gcd_euclid`.
- Proved the complete raw helper refines normalized polynomial GCD, including
  both input copies, every generated `_poly_divrem`, buffer rotations, and the
  final copy to `G`.
- Pinned the new execution and refinement theorem in the strict source gate.

## Why

SQF calls polynomial GCD.  A strict connection cannot use the L2 gcd as an
implementation or stop at the Euclid loop: it must account for the actual
C++ helper's local allocations, copies, output pointer, output length, and
heap effects.

## Key decisions

- The execution function is the generated raw helper instantiated only with
  the proof of decrease derived from successful generated division.
- Workspace assumptions contain allocation validity and non-aliasing only;
  no semantic result is supplied by the caller.
- The L2 gcd occurs solely in the postcondition and is propagated through the
  division identity proved from raw execution.
- Termination uses the actual divisor length.  No fuel, `partial`, oracle, or
  alternate L2 execution path is introduced.

## Files

- `proof/lean/CLPoly/Impl/StrictEuclidRefinement.lean`
- `proof/cpp2lean_v2/tests/check_strict_hgcd_source.py`
- `docs/devlog/2026-08-09-strict-gcd-euclid-refinement.md`

## 度量

- 耗时：约 2 小时（证明草稿、raw copy/frame 引理、递归位置证明、构建）
- 迭代：3 轮编译—修复循环
- Lean 新增/修改行数：228 行
- 对应 C++ 行数：35 行（完整 `_gcd_euclid` helper）
- 放弃的方案：依赖 HGCD 文件中已有的 address-slice copy 引理；该方向会
  产生导入环。最终直接复用更底层的 region-separated raw copy 语义。
