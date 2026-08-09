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

## Main-loop refinement proof draft

For every large-loop HGCD call, use the polynomial representations produced
by the immediately preceding raw divrem call together with the physical
recursive-invocation workspace.  Apply the refinement theorem for the exact
successful `hgcdRecursiveCallChecked` execution.  Its returned invariant says
that the actual output second length is less than `lenJ / 2 + 1`.  The source
large branch establishes `hgcdRecursiveCutoff <= lenJ`; since the cutoff is
100, arithmetic gives `lenJ / 2 + 1 < lenJ`.  Transitivity therefore proves
that the exact returned `lenB` strictly decreases the source loop measure.
No alternative execution or semantic result is chosen in this argument.

The semantic loop proof then follows the source branches by well-founded
induction on `lenJ`: raw divrem changes `(G,J)` to `(J,R)` while preserving the
normalized gcd; zero remainder copies `J`; the small branch invokes the exact
raw Euclid refinement; the large branch invokes the concrete checked HGCD
refinement and recurses on its strictly smaller returned second length.

The loop-level physical providers are deliberately conditional on the raw
representations at the currently reached heap.  They provide only capacity
and non-aliasing for divrem and the nested Euclid helper.  They do not contain
quotient, remainder, gcd, output lengths, or an execution result.  A bundled
large-step lemma now obtains both gcd preservation and strict decrease from
one and the same successful checked HGCD call.

The first composition layer now has branch-local theorems for the actual
loop-head divrem, the zero-remainder copy, and the small raw-Euclid dispatch.
Together with the bundled checked-HGCD theorem, all four source branches have
semantic lemmas tied to their exact raw executions.  The remaining work is to
thread these lemmas through one well-founded loop theorem and replace the
raw entry's freely supplied decrease contract with that construction.

Fixed-capacity providers were then generalized to exact-size dynamic
providers, because checked HGCD returns new descriptor lengths in a new heap.
The dynamic contracts expose precisely the divrem slices, the nested Euclid
allocation, and the concrete recursive HGCD invocation/descendant workspaces
at each represented reached state.  The divrem, Euclid, and HGCD branch
theorems now have dynamic forms with no inherited-capacity premise.

Composition also exposed that the checked-HGCD invariant theorem was still
conditional on a successful raw return.  The first totality bridge now proves
that the exact generated recursive base arm succeeds from its physical base
workspace, for both `computeM` values.  Non-base totality remains to be built
well-foundedly from the two smaller child executions and the stored concrete
continuation executions before the main GCD loop can be called total.

The totality stack now also contains generated `_poly_add`: its common loop,
selected long-input tail copy, and final normalization are each constructed
from raw slice validity.  Using the existing total raw multiplication theorem
and this addition result, `_mat_row_update` is total in both its inactive
descriptor-swap branch and its active raw-multiply/raw-add branch.  Neither
theorem accepts an execution equality as a premise.

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
