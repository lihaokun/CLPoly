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

### Staged HGCD iteration totality proof draft

The existing loop workspace provider is intentionally useful for conditional
semantic refinement, but it is indexed by an already successful division and
two already successful row updates.  It therefore cannot establish totality.
For totality, use a staged physical provider instead.  Its first stage contains
only the source heap's division buffers and separation facts.  After the real
division returns its concrete heap and lengths, the second stage supplies the
physical workspace for row `(2,3)`.  After that real row update returns its
concrete matrix and heap, the third stage supplies the workspace for row
`(0,1)` and the four external-polynomial frame guards.  No stage supplies a
polynomial value, output length, result record, or execution equality.

Run the generated divrem from the represented operands.  Its refinement gives
the quotient and remainder representations, preserves every matrix entry, and
proves the actual `lenR < lenB`.  Run row `(2,3)` using total
`_mat_row_update`; frame the quotient and the untouched entries `(0,1)`, then
translate those representations through the returned descriptor equalities.
Run row `(0,1)` using the same total theorem.  The staged guard workspaces now
form the old conditional iteration workspace, so the existing exact-call
semantic theorem yields the next matrix transform, determinant sign, gcd, and
raw operand representations.  Thus one source iteration both succeeds and
returns the invariant required for a recursive call on the strictly smaller
real remainder length; no fuel or alternate semantic execution is involved.

This construction is now formalized by `hgcdIteration_succeeds` and lifted to
the complete generated while-loop by `hgcdIterLoop_succeeds`.  The latter is a
well-founded Lean definition whose measure is the current raw state's `lenB`;
its recursive obligation is discharged by the `lenR < lenB` fact returned by
the same concrete divrem execution.  The strict source gate, focused build,
and axiom audit all pass; the new theorems depend only on `propext`,
`Classical.choice`, and `Quot.sound`.

For the enclosing `_hgcd_iter`, first execute the real identity-matrix
initialization and the two ordered input copies using `hgcdIterInit_refines`.
Those exact writes produce the initial raw identity transform.  Feed that
returned state—not a reconstructed alternative—into the total loop theorem,
then rewrite the generated `_hgcd_iter` match with the initialization equality
and the returned loop equality.  This yields total execution and the final
raw transform/gcd invariant for the complete generated helper.

### Staged reconstruction-pair totality proof draft

The four leaf helpers used by the paired reconstruction are already total:
the two sign-selected low reconstructions execute real multiplication and
subtraction, while each shifted-high lift executes real zero fill, addition,
and normalization.  The remaining circularity is only in the old pair
provider, which returns its physical facts after all four calls have already
succeeded.  Replace it with four stages.  The first gives the B reconstruction
workspace.  Its actual returned heap frames the high-B input and supplies the
B-lift workspace.  That actual lift frames the matrix and both low inputs and
supplies the A-reconstruction workspace.  The returned A heap frames high A
and supplies the A-lift workspace.  Final prefix facts frame the completed B
output and the child matrix across the last two calls.  Each stage depends
only on prior concrete executions, never on a later success or semantic
result.  Composing the four existing total leaf theorems then constructs the
exact generated pair result and its two raw polynomial representations.

For the middle block between recursive children, run the generated divrem
directly from `HgcdRecursiveMiddleWorkspace`.  Package its actual heap and
lengths with the source pointer arithmetic for `k`, `c0`, and `d0`; divrem's
own theorem supplies quotient/remainder representations, gcd preservation,
layout preservation, and the strict remainder decrease.  No separate middle
result or success witness is accepted as input.

For one quotient-matrix entry, reuse the already total raw row-update
arithmetic with `bottom` as the multiplied entry and `top` as the in-place add
destination.  In the active branch, extract the exact multiplication and
addition executions from that call and replay them through the generated
quotient-entry constructor, whose only descriptor change is `top.len`.  In the
inactive branch the generated helper returns the input matrix immediately.
This proves totality without assuming a quotient-entry result.

For the two-column quotient application, provide the first column workspace
up front and the second only after the actual first result.  Execute and refine
column `(0,2)`, then use its returned matrix and preserved quotient to execute
and refine `(1,3)`.  The generated wrapper is exactly the second result with
its validity witness, so this staged composition proves the complete quotient
application total and semantic.

### Staged matrix-entry totality proof draft

For `_mat_mul_entry`, supply the first guarded-product workspace before any
execution.  After its actual result, frame the second product's two inputs and
supply its workspace.  After the actual second result, frame the stored first
product and supply only the physical facts for the selected tail: common
destination capacity, add aliasing, and copy disjointness.  If both products
are nonzero, run total `_poly_add`; if only the first is nonzero, return its
existing buffer; if only the second is nonzero, run the real copy; if both are
zero, return length zero.  Repackage these staged facts as the older
conditional semantic workspace only after both product executions have been
constructed, then apply the exact-call semantic theorem.  Thus no product or
tail success is assumed.

For complete `_mat_mul`, provide both the staged total entry workspace and the
existing purely spatial step frame at every concrete loop state.  Execute the
entry for index `i`, use the frame theorem on that same execution to retain
both input matrices, update only the generated output length descriptor, and
continue at `i+1`.  The source loop has exactly four entries, so recursion is
on `4-i`.  Once its total raw execution is constructed, apply the existing
exact-run matrix refinement to obtain all four product entries from that same
execution.

For the final combine block, execute the now-total two-column quotient update
on `S`.  Its actual result frames every entry of the left matrix `R` and
supplies the total four-entry multiplication workspace for `R * modifiedS`.
Run that generated matrix multiplication and package the two exact calls as
`hgcdRecursiveCombineMatrix`; the returned raw matrix is the L2 product already
proved for those same calls.

For the complete recursive finish block, first execute the four-call
reconstruction pair unconditionally, because that is the first operation in
the generated C++ tail.  Frame `R`, `S`, and the quotient across that actual
execution.  When `computeM` is false, return the reconstructed operands and
the incoming matrix exactly as the generated branch does.  When it is true,
feed the actual reconstructed heap to the total quotient/matrix-combine
workspace, execute that block, and frame both reconstructed output operands
across its actual result.  The theorem therefore constructs every success
equation before extracting its polynomial or matrix semantics; it assumes no
success of the finish block or either of its children.

For a cutoff dispatch's iterator arm, execute the complete well-founded
iterator first.  Its actual matrix and operand descriptors determine the
physical stabilization workspace.  Run the generated four-entry stage and
restore loops, then use the actual stable heap to run the generated
alias-sensitive output normalization.  Only after all three calls have been
constructed, repackage their staged workspaces as the older conditional
iterator-arm contract to reuse its length, transform, determinant, and GCD
semantic theorem.  Thus the cutoff arm becomes total without assuming the
success of the iterator, stabilization, or output copies.

For either cutoff dispatch, use a total callback contract for the large arm:
it returns an actual smaller generated execution together with the invariant
of that same result.  If the source cutoff selects the small arm, invoke the
total iterator-arm theorem; otherwise invoke that smaller-call contract.
Package the exact branch selected by `hgcdRecursiveDispatchBelow`, using
proof irrelevance only for generated validity witnesses.  This gives the
parent a concrete first or second child instead of merely explaining a child
whose success was assumed.

For the first child of a non-base recursive body, instantiate that total
dispatch at the source high-half lengths.  Use its actual heap frame to carry
the untouched low halves forward, then execute the four generated paired
reconstruction calls with the actual child matrix and high-half outputs.
Finally apply the existing proof-only reconstruction bound to those same two
executions.  The resulting strict outer decrease is therefore derived from
the actual child and reconstruction, not supplied for a preselected record.

For the early continuation, recover the full-input transform, signed
determinant, and GCD relation from that same first child and reconstruction.
Use the reconstructed stopping guard to derive the outer length invariant,
then execute the generated early-return copies (A, B, and optionally four
matrix entries).  Package that actual result as the common recursive
invariant; no successful early result is accepted as an input premise.

## Files

- `proof/lean/CLPoly/Generated/StrictGCDHGCD.lean`
- `proof/lean/CLPoly/Impl/StrictGCDHGCDRefinement.lean`
- `proof/cpp2lean_v2/tests/check_strict_hgcd_source.py`
- `docs/devlog/2026-08-09-strict-gcd-hgcd-raw-entry.md`

## 度量

- 耗时：约 1.5 小时（源码控制流映射、良基主循环、构建）
- 迭代：3 轮编译—修复循环
- Lean 新增/修改行数：约 110 行
- 对应 C++ 行数：约 85 行（`_gcd_hgcd`）
- 放弃的方案：把整个 HGCD GCD 改写成单纯 Euclid；该方案不对应 C++
  的阈值分派和 HGCD 主循环。
- 本轮增量：约 2 小时、4 轮编译—修复、Lean 约 330 行；完成 staged
  单迭代总执行及真实余式长度上的完整 iterator 良基总执行。
- 后续增量：约 2 小时、7 轮编译—修复、Lean 约 320 行；完成矩阵乘法
  单项、四项循环以及最终商更新后矩阵组合块的 staged 总执行。所有成功
  等式均由实际生成函数调用构造，旧的条件式语义工作空间只在调用成功后
  用于提取 L2 语义。
- finish 增量：约 1 小时、5 轮编译—修复、Lean 约 160 行；完成真实
  四调用重构、`computeM` 两分支和可选最终矩阵块的总执行，并在矩阵写入
  前后显式保持重构后的 A/B 表示。
- cutoff iterator 增量：约 1.5 小时、3 轮编译—修复、Lean 约 230 行；
  将良基 iterator、真实矩阵稳定化和别名敏感输出复制组成一个总执行，
  再从同一执行提取递归长度、变换、行列式与 GCD 不变式。
- dispatch 增量：约 0.5 小时、2 轮编译—修复、Lean 约 70 行；以总执行
  callback 统一 cutoff 小分支和严格变小的大分支，使父调用得到实际 child
  结果及同一次执行的完整递归不变式。
- first child 增量：约 1 小时、2 轮编译—修复、Lean 约 150 行；执行
  第一次总 dispatch、保持原始低半部、完成四调用 paired reconstruction，
  并从这两个实际执行导出严格的外层重构长度界。
- early 增量：约 1 小时、2 轮编译—修复、Lean 约 130 行；从 first child
  与重构的实际执行恢复完整输入变换/GCD，再执行 A、B 与可选矩阵复制，
  生成同一次 early result 的公共递归不变式。
