# Strict `_gcd_hgcd` raw execution

## Date

2026-08-09

## What changed

- Added the raw result record for `_gcd_hgcd`.
- Added the exact well-founded lowering of its `while (len_j != 0)` loop.
- Added the complete raw entry: initial division, zero-remainder copy, first
  concrete HGCD call, loop divisions, small Euclid dispatch, repeated concrete
  HGCD calls, and final output length.
- Exposed all C++ workspace/vector allocations as raw pointer parameters.
- Pinned the execution and its source-length termination contract in the
  strict source gate.

## Why

The completed Euclid and concrete `_hgcd_recursive` proofs must be composed
through the actual `_gcd_hgcd` dispatch before top-level polynomial GCD or SQF
can use them.  Stopping at either helper would omit the source main loop.

## Key decisions

- The main measure is the actual source variable `len_j`.
- The termination contract is indexed by the successful concrete HGCD call;
  it cannot provide a different heap, length, or polynomial result.
- The small branch invokes the actual raw `_gcd_euclid` execution.
- No L2 operation occurs in this execution definition.

## Files

- `proof/lean/CLPoly/Generated/StrictGCDHGCD.lean`
- `proof/cpp2lean_v2/tests/check_strict_hgcd_source.py`
- `docs/devlog/2026-08-09-strict-gcd-hgcd-raw-entry.md`

## 度量

- 耗时：约 1.5 小时（源码控制流映射、良基主循环、构建）
- 迭代：3 轮编译—修复循环
- Lean 新增/修改行数：约 110 行
- 对应 C++ 行数：约 85 行（`_gcd_hgcd`）
- 放弃的方案：把整个 HGCD GCD 改写成单纯 Euclid；该方案不对应 C++
  的阈值分派和 HGCD 主循环。
