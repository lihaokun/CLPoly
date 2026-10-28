# Strict GCD Euclid raw entry

## Date

2026-08-09

## What changed

- Added an explicit result record for the raw `_gcd_euclid` helper.
- Added its source-ordered raw execution: copy both immutable inputs into the
  two local buffers, run the existing well-founded Euclid rotation, and copy
  the last live dividend into `G` while returning `len_G`.
- Extended the pinned HGCD/GCD source gate to scan the generated Euclid,
  divrem refinement, and Euclid refinement files and require the raw entry,
  termination certificate, and end-to-end loop theorem.

## Why

The object-level `dense_upoly_zp::gcd` probe still lowers vector storage to
default-valued arrays and loses the `_gcd_hgcd` output mutation.  The strict
path therefore has to connect the already verified raw helpers before an
object wrapper can be admitted.

## Key decisions

- Vector allocations are exposed as caller-provided raw regions with explicit
  validity and separation obligations.
- The wrapper reuses `euclidLoop`, whose termination measure is the actual
  divisor length and whose decrease follows the generated divrem execution.
- No object-level `partial def`, array default behavior, or specification
  implementation is introduced.

## Files

- `proof/lean/CLPoly/Generated/StrictEuclidGCD.lean`
- `proof/cpp2lean_v2/tests/check_strict_hgcd_source.py`
- `docs/devlog/2026-08-09-strict-gcd-euclid-raw-entry.md`

## 度量

- 耗时：约 1 小时（源代码核对、raw 入口建模、构建与审计门）
- 迭代：2 轮编译—修复循环
- Lean 新增/修改行数：29 行
- 对应 C++ 行数：35 行（`dense_upoly_zp::_gcd_euclid` 主体）
- 放弃的方案：直接采用对象级生成结果；该结果把 vector 降成默认数组并丢失
  `_gcd_hgcd` 的输出副作用，不能作为严格 L1 执行语义。
