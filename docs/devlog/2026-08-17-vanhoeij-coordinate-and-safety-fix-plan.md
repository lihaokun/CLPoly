# Van-Hoeij coordinate provenance and full-precision safety fix plan

## Problem

`__vanhoeij_recombine` retained an already LLL-reduced matrix between failed
CLD rounds.  The next `__lll_reduce` initialized a fresh local transform `U`,
but candidate extraction interpreted its first columns as coordinates in the
original active lifted-factor basis.  After the first nontrivial reduction
those coordinate systems need not agree.

At the top level, `__lll_factorize` could also return a reduced-cardinality
van-Hoeij result at complete Mignotte precision without executing the existing
exhaustive Zassenhaus safety net.  The current irreducibility proof covers the
full-cardinality modular case and the concrete Zassenhaus terminal case, but
does not justify accepting every such intermediate LLL grouping.

## Correction

1. Rebuild the canonical scaled diagonal lattice at every van-Hoeij retry and
   append all requested CLD columns from zero through the current target.
   Consequently each returned `U` is relative to the active lifted factors.
2. At complete precision, route every reduced-cardinality result through the
   literal existing `__zassenhaus_recombine` execution.
3. Mirror the changed source control flow in the generated strict FactorZZ L1.
4. Add native regression coverage and C++/Lean B2B coverage for the new safety
   decision.  This targeted B2B is additional to, not a substitute for, the
   final full FactorZZ B2B.

## Correctness argument

Rebuilding changes no mathematical lattice for a given target: it constructs
the same scaled factor-coordinate rows and the same ordered CLD prefix, but
prevents a local transform from being confused with a cumulative transform.
The safety route executes an already generated and proved exhaustive
recombination algorithm on the same full-precision Hensel factors, so it does
not introduce a semantic oracle or an alternative factorization witness.

## Expected cost

Failed CLD retries redo earlier LLL work.  This is a correctness-first change;
future optimization may retain the matrix only if it also composes transforms
and proves the corresponding coordinate invariant.

## Implementation and verification

The correction was implemented in the production C++ controller and mirrored
in the generated strict FactorZZ L1.  `FactorZZRawOps` now exposes the concrete
Zassenhaus operation, so the full-precision safety branch executes that
operation rather than assuming an abstract factorization result.  The
low-precision control-flow lemma was extended to account for both the direct
van-Hoeij result and the executed Zassenhaus branch.

The B2B driver build now links the two out-of-line CLPoly translation units it
actually needs (`upolynomial.cc` and `variable.cc`) instead of depending on a
prebuilt `libclpoly.so`.  This makes the source/Lean comparison reproducible on
the current platform even though the unrelated full-library build currently
fails in `realroot.hh`.

Verification performed:

- native recombination regression: 116/116 passed, including a modular-factor
  grouping case and all four safety-decision cases;
- targeted C++/Lean B2B: 4/4 passed for
  `__needs_zassenhaus_safety_net`;
- all existing effective B2B cases plus the new cases: 141 passed; the sole
  runner failure is the pre-existing `_pipeline_smoke.json`, which deliberately
  requests `_nonexistent_fn` and receives `unknown fn` from both drivers;
- strict FactorZZ generator `--check`: passed;
- `CLPoly.Generated.StrictFactorZZ`, `CLPoly.Refinement.FactorZZ`,
  `B2B.Registry`, and `B2B.Driver`: built successfully;
- B2B codec round trips, including both Bool values: passed;
- strict-boundary audit: passed (no legacy dispatch, `partial def`, fuel,
  `sorry`, or forbidden fallback interface);
- refinement axiom audit: passed with only Lean/Mathlib kernel-standard
  `propext`, `Classical.choice`, and `Quot.sound` dependencies;
- `git diff --check`: passed.

## Files

- `clpoly/polynomial_factorize_univar.hh`
- `test/test_recombine.cc`
- `proof/cpp2lean_v2/tests/build_strict_factor_zz.py`
- `proof/lean/CLPoly/Generated/StrictFactorZZ.lean`
- `proof/lean/CLPoly/Refinement/FactorZZ.lean`
- `proof/b2b/Makefile`
- `proof/b2b/registry/cpp_registry.cc`
- `proof/b2b/types/b2b_types.hh`
- `proof/b2b/types/b2b_types.cc`
- `proof/b2b/types/test_types.cc`
- `proof/b2b/vectors/__needs_zassenhaus_safety_net.json`
- `proof/lean/B2B/Types.lean`
- `proof/lean/B2B/Registry.lean`

## 度量

- 耗时：约 2.5 小时（问题复核、C++ 修正、生成同步、证明修正和回归）。
- 迭代：2 轮 Lean 编译-修复；3 轮 C++/B2B 构建-修复。
- Lean 新增/修改行数：约 64 行（含生成模板、生成产物、证明和 B2B）。
- 对应 C++ 行数：约 55 行（生产控制流、测试、codec、registry 和构建规则）。
- 放弃的方案：继续跨轮复用约化矩阵；由于需要额外组合并证明坐标变换，当前采用每轮重建规范矩阵。另放弃依赖完整 `libclpoly.so` 的 B2B 构建，改为直接链接实际需要的编译单元。
