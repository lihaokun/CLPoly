# Strict polynomial GCD object bridge

## Date

2026-08-09

## Natural-language proof draft

The C++ `polynomial_GCD` path converts each canonical sparse polynomial into
a dense coefficient vector, calls `dense_upoly_zp::gcd`, scales its dense
output, and converts the resulting dense vector back to descending sparse
terms.  The existing generated SQF instead resolves `polynomial_GCD` through
the hand-written `SparsePolyZp.gcd` instance, so it cannot be used as strict
L1 execution evidence.

Introduce one representation relation tying a canonical `SparsePolyZp` value
to the exact normalized raw dense buffer whose mathematical polynomial is
`SparsePolyZp.toPoly`.  Its projections give the sparse canonical invariant
and the already verified `RawDensePolyRep`; it contains no gcd result or
algorithm callback.  Introduce a corresponding dense-to-sparse output
relation, again requiring both the exact `toPoly` equality and canonical
sparse form.  The raw GCD theorem can consume the first relation immediately;
the forthcoming constructor and `to_upoly` loop proofs must establish these
relations from their actual writes/reads.

This separation prevents either conversion from being replaced by a
noncomputable polynomial encoder and prevents `SparsePolyZp.gcd` from serving
as the executable implementation.

For dense-to-sparse execution, mirror `dense_upoly_zp::to_upoly` with a
structurally recursive reverse scan.  With `remaining = n+1`, read exactly
`coeffs[n]`; skip zero, otherwise append `(n, Zp(value,p))`, then recurse with
`n`.  Validity of the complete raw slice gives validity of every actual read,
so induction on `remaining` proves the scan cannot fault.  The returned array
is produced solely by those reads and pushes; its polynomial and canonical
properties are proved in the following refinement stage.

For semantic refinement, strengthen the scan induction over an arbitrary
already-produced high-degree accumulator.  A represented prefix of length
`n+1` is its length-`n` prefix plus the monomial read at index `n`.  The zero
branch drops a zero monomial; the nonzero branch uses `Array.push` and
`listSum_append` to add that exact monomial.  Thus the final sparse `toPoly`
equals the complete raw `SlicePolyRep`, with no reconstruction oracle.

For canonicality, carry the additional scan invariant that every accumulated
term has degree at least `remaining`.  A nonzero read at index `remaining-1`
is therefore strictly below all accumulated terms; appending it preserves
strict descending order.  Raw canonicality supplies `value < p`, and the
source guard supplies `value != 0`.  Combine this proof with the semantic
proof by uniqueness of the actual successful scan, yielding the complete
`RawDenseSparseResult` for one executable output.

## What changed

- Added strict sparse/raw-dense input and output representation relations.
- Added projection and uniqueness lemmas used by the forthcoming generated
  conversion loops.
- Added the exact reverse dense-to-sparse scan and its total raw execution
  theorem.
- Proved the scan's accumulator invariant and end-to-end `toPoly` equality
  with the normalized raw input polynomial.
- Proved canonical output and combined both properties for the same actual
  `to_upoly` execution.
- Imported the completed raw HGCD-GCD refinement at the squarefree boundary.

## Why

SQF must call the actual C++ GCD path.  A precise object/raw representation
boundary is required before replacing its current typeclass GCD call.

## Files

- `proof/lean/CLPoly/Impl/StrictPolynomialGCDRefinement.lean`
- `proof/lean/CLPoly/Refinement/SquarefreeZp.lean`
- `docs/devlog/2026-08-09-strict-polynomial-gcd-object-bridge.md`

## 度量

- 耗时：约 0.5 小时（源码控制流核对、接口设计、形式化与构建）
- 迭代：1 轮编译—修复
- Lean 新增/修改行数：约 55 行
- 对应 C++ 行数：约 55 行（两个 sparse/dense 转换及 GCD 包装）
- 放弃的方案：直接证明当前 `SparsePolyZp.gcd` 正确；它不是 C++ dense
  GCD 的执行，不能用于严格 L1→L2 精化。
