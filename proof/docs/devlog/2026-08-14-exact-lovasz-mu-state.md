# Exact generated Lovász μ state

The strict recombination refinement now identifies the concrete array value
produced by the generated C++ Lovász-swap μ update.

- `swapQQRows_ok`, `swapQQRows_size`, and `swapQQRows_get` expose the exact
  successful execution of the generated row-swap helper.
- `swapRowsArray` is a total pure representation of that bounded operation.
- `updateMuAfterSwapRow` and `updateMuAfterSwapArray` characterize the two
  coefficient corrections on every later row, with exact size and lookup
  theorems.
- `lovaszSwapMuResult_of_generated` composes the generated row swap, boundary
  coefficient write, and well-founded post-swap loop into one concrete output
  equality.

This layer does not assert an abstract LLL result or use an existence oracle;
it records the array computed by the generated C++ control flow.  The next
step uses these equations to prove preservation of the concrete
Gram--Schmidt factorization.

The follow-up lookup layer separates the intermediate corrected swap from the
tail loop and proves the complete entry formula.  In particular, it exposes
the exchanged rows, the new boundary coefficient, and both corrected
coefficients in every row strictly after `k`.  Array and row sizes are
preserved throughout, so later matrix proofs can rewrite generated storage
without any unchecked read.

The local rational algebra is now closed as well.  `lovaszLocalTransform`
defines the identity column transformation with the exact adjacent `2 × 2`
block used by the swap.  `lovaszLocalTransform_diagonal` proves that
conjugating the old diagonal norms by this block gives precisely the two
generated replacement norms, while `mul_lovaszLocalTransform_apply` exposes
the corresponding right-multiplication formula for every matrix entry.

`gsLowerPrefix_lovaszSwapMuResult_mul` now connects that algebra to generated
storage.  It proves entry by entry that the lower factor built from the exact
post-loop μ array, followed by the local column transform, equals the old
lower factor with the two adjacent basis rows exchanged.  The proof covers
the two exchanged rows and every earlier/later row separately and uses the
generated lookup equations for both corrected columns.

The same local transform is now connected to the generated norm array, and
`concreteGramSchmidt_lovaszSwap_of_above` composes both identities with the
actual integer row swap.  For every prefix containing both adjacent rows it
proves the full rational Gram matrix equality `G' = L' D' L'ᵀ`, rather than
only determinant or potential preservation.

The remaining prefix cases are now closed as well.  Prefixes strictly before
the exchanged pair are unchanged entry by entry.  The boundary prefix that
contains the new row `k - 1` but not row `k` is obtained by restricting the
proved `(k + 1)`-dimensional factorization and eliminating the final triangular
summand.  `concreteGramSchmidt_lovaszSwap` combines these cases into the full
executable Gram--Schmidt invariant for the exact generated swap state.

Finally, `concreteLLLExecutionValid_lovaszSwap` lifts this matrix identity to
the complete executable invariant: the swapped integer matrix remains square,
the exact generated μ array remains square with the same dimension, the
generated norm array keeps its size and strict positivity, and every prefix
retains the concrete Gram--Schmidt factorization.

The proof now directly inverts the generated `lllStep` failed-Lovász branch.
It follows the successful size reduction, both integer-matrix swaps, the μ-row
swap, correction, and `updateMuAfterSwapLoop`, and identifies the returned
state with the concrete swap state above.  Together with the already proved
advanced branch, `concreteLLLTermination` instantiates the generated
`LLLTermination` interface using `ConcreteLLLExecutionValid` and the concrete
determinant/index lexicographic rank.  Thus the generated `lllMainLoop` is now
justified by genuine well-founded recursion, not fuel or a termination oracle.
`concreteLLLMainLoop_preserves_execution_valid` then follows the generated
recursive execution itself and proves that every successful final state still
carries the same concrete invariant.

The generated recombination layer now exposes the complete executable wrapper
around that main loop.  `lllReduce` runs the source Gram--Schmidt
initialization, enters the well-founded LLL loop, scans every reduced row for
the concrete squared-norm bound, and deterministically orders the accepted row
indices by the source strict norm comparator.  `extractCandidates` then runs
the actual nested column-equivalence loops over the unimodular transform and
collects the resulting factor-index classes.  The former
`VanHoeijRawOps.prepareCandidates` result callback has been removed, so no
operation field can choose candidate subsets.

This audit also found that the earlier generated `VanHoeijState` omitted the
source variables `M` and `J_cur`.  That would have discarded previously added
CLD columns after an unsuccessful round.  The state now stores the lattice,
current CLD-column count, and short-vector bound.  An unsuccessful round keeps
the reduced lattice and accumulated columns; a successful extraction rebuilds
the scaled diagonal lattice and resets the column count, exactly as the C++
loop does.  The remaining proof-only `LLLExecution.inputValid` interface is
explicitly restricted to admissible source lattices, with separate preservation
obligations for CLD extension, LLL output, and lattice reset.  These obligations
cannot return executable data and replace the invalid earlier requirement that
Gram--Schmidt initialization succeed with positive norms for every arbitrary
integer matrix.

The proof layer now fixes that admissibility predicate concretely:
`ConcreteLLLInputValid` means that the raw array is square and its full integer
basis determinant is nonzero.  Any state satisfying the already established
executable Gram--Schmidt invariant yields this input condition: the determinant
of its full Gram matrix is the product of its strictly positive generated
norms, while the same Gram determinant is the square of the basis determinant.
Consequently a completed LLL round supplies the exact precondition needed when
the next C++ call reinitializes Gram--Schmidt.  `concreteLLLExecution` now
instantiates the generated wrapper modulo one sharply isolated theorem,
`LLLInitializationCorrect`, asserting correctness of the concrete generated
initialization loops on a square full-rank basis.

The initialization proof has now entered the generated code at its innermost
operation.  `dotRowsLoop_eq_Ico_sum` follows every checked array read of the
well-founded C++ dot-product loop and identifies its accumulator with the
remaining interval sum.  `dotRows_eq_fin_sum` also proves that the loop cannot
fault when the right row is at least as long as the left row and returns the
exact finite inner product used by the Gram matrix.  This removes the first raw
execution boundary needed by the μ-row and norm initialization invariants.

`gramNumeratorLoop_eq_Ico_sum` now follows the next generated loop and proves
that every successful execution subtracts exactly
`Σ l<j, μ[i,l] μ[j,l] B[l]` from the raw row inner product.
`gramNumeratorLoop_succeeds` discharges each of its six checked array accesses
from square μ storage and norm-size invariants, and
`gramNumeratorLoop_exact` packages the result as the finite sum indexed by
`Fin j`.  Thus the generated numerator is fixed by executable storage; no
abstract Gram--Schmidt coefficient is selected by the proof.

The analogous norm pass is closed too.  The three
`gramNormLoop_*` theorems follow the generated checked reads, prove totality
under the square-storage invariant, and identify the final value with
`⟨b_i,b_i⟩ - Σ j<i, μ[i,j]^2 B[j]`.  Both scalar inner loops used by C++
Gram--Schmidt initialization are therefore now concrete finite-sum
computations.

`GramStorageShape` now records the exact square μ-array and norm-array sizes
used by the generated outer loops, and `GramStorageShape.setMu` proves that the
source's nested `Array.set` preserves them.  Using the already proved dot and
numerator executions, `gramMuRowLoop_succeeds` follows every branch of the
actual `j < i` loop, proves it cannot fault on a square lattice, and shows that
the returned norm array is unchanged while the μ storage remains square.  The
next refinement step strengthens this execution theorem with the lookup
formula for every coefficient written in row `i`.

The generated outer initialization loop is now executable under exactly the
same storage conditions.  `GramStorageShape.setNorm` records that the source
norm write preserves the array dimensions.  `initializeGramSchmidtLoop_succeeds`
then follows each concrete row through `gramMuRowLoop`, the checked integer dot
product, `gramNormLoop`, the actual norm `Array.set`, and the recursive call on
the next row.  Its well-founded induction measure is the number of unprocessed
rows, `size - i`; every source step changes it to `size - (i + 1)`.  The theorem
proves no-fault execution and shape preservation only—it does not assume the
Gram--Schmidt coefficients, norm values, positivity, or LDLᵀ identity that the
next semantic invariant must derive from these concrete computations.

## 度量

- 耗时：约 1.5 小时（外层控制流组合、数组边界调试与单文件复验）
- 迭代：约 8 轮编译-修复循环
- Lean 新增/修改行数：约 64 行
- 对应 C++ 行数：约 15 行 Gram--Schmidt 外层初始化控制流
- 放弃的方案：直接在外层定理中同时闭合 LDLᵀ 语义；拆分无故障/形状与数学语义后，循环边界义务更清晰，且不会引入抽象执行接口

## Generated initialization semantic base

The semantic induction now starts from the actual generated storage rather
than an assumed Gram--Schmidt result.  Well-founded proofs for
`zeroQQRowLoop` and `zeroQQMatrixLoop` establish their exact sizes, preserve
pre-existing prefixes, and identify every appended cell as zero.  Consequently
the concrete μ array used by `initializeLLL` is proved square and pointwise
zero directly from the generated loops.

`ConcreteGramSchmidtUpTo` restricts the exact `G = L D Lᵀ` equation to the
prefix already processed by the source outer loop, while
`ProcessedGramSchmidtValid` combines that equation with exact array shape and
positivity of every completed norm.  The initial state is proved valid through
prefix one: the zero μ row makes the lower factor the one-by-one identity and
the source `dotRows matrix[0] matrix[0]` value is exactly its diagonal norm.

Positivity is not assumed.  A nonzero determinant of the full integer input
basis gives linear independence of its actual rows, hence row zero contains a
nonzero stored coefficient.  Its finite sum of integer squares is therefore
strictly positive, and casting this concrete sum to `QQ` proves positivity of
the first generated norm.  This supplies the base case required for the
row-by-row semantic extension without an existence witness or semantic oracle.

## 度量（semantic base）

- 耗时：约 2 小时（零数组执行证明、前缀不变量设计、首行正性与编译调试）
- 迭代：约 9 轮编译-修复循环
- Lean 新增/修改行数：约 205 行
- 对应 C++ 行数：约 12 行零矩阵与首行 norm 初始化
- 放弃的方案：把初始 μ 视为数学零矩阵直接化简；改为逐个证明生成循环的 size/get 语义，确保基例仍来自真实 C++ lowering

## Generated μ-row write semantics

The no-fault μ-row theorem now carries exact frame information through every
recursive source iteration.  Every row other than the current target row is
unchanged, and every target-row column preceding the loop's starting index is
unchanged.  These are proved by following the nested `Array.set` operations:
the outer write selects only row `i`, the inner write selects only column `j`,
and the recursive theorem transports both frames through all later writes.

`sourceGramCoefficient` names the closed form of exactly the scalar operations
performed by one generated iteration: the integer row dot product, the
generated `l < j` numerator subtraction, and the source zero-norm conditional
division.  It does not compute a specification result or choose an output.
`gramMuRowLoop_written_coefficient` inverts the actual generated execution at
column `j`, applies the already proved exact dot/numerator loops, and then uses
the recursive prefix frame to show that no later iteration overwrites the
cell.  Thus the final returned array's `μ[i,j]` is now identified with the
actual source computation.

## 度量（μ-row write semantics）

- 耗时：约 1.5 小时（双重数组帧、当前写入识别与依赖索引调试）
- 迭代：约 7 轮编译-修复循环
- Lean 新增/修改行数：约 145 行
- 对应 C++ 行数：约 11 行 μ-row 内层循环
- 放弃的方案：直接把整行结果声明为数学 Gram--Schmidt 系数；改为先证明每个生成写入的闭式值和后续不覆盖性质，再由逐列不变量推出 LDLᵀ

## Full-rank prefixes have positive Gram determinant

The positivity argument has now been generalized from the first generated
norm to every source row prefix.  `basisPrefixMatrixQQ` is the rectangular
matrix containing the first `r` actual integer basis rows, cast cell-by-cell
to `QQ` while retaining all source columns.  Its product with its transpose is
proved entrywise equal to the exact generated prefix Gram matrix.

The input determinant certificate first gives linear independence of all
integer rows.  Faithful base change transports this fact to `QQ`, and
restriction along the concrete `Fin r -> Fin n` row embedding proves that
every prefix remains linearly independent.  A second faithful base change to
the reals permits use of positive-definite Gram-matrix theory: the mapped
matrix `B * B^T` is positive definite, so its determinant is strictly
positive.  `Rat.cast_det` and order reflection then return the result to the
project's executable `QQ` values.  No orthogonal basis, pivot, or output is
selected existentially.

This is the positivity source needed by the next row-extension theorem: after
the generated coefficient and norm loops establish the next exact
`G = L D L^T` prefix, determinant multiplicativity will force the newly
written pivot to be positive from this theorem and the already-positive
earlier pivots.

## 度量（prefix Gram positivity）

- 耗时：约 1 小时（整基换标量、前缀限制、Gram 正定性与 `QQ`/`Real` 桥）
- 迭代：约 6 轮编译-修复循环
- Lean 新增/修改行数：约 99 行
- 对应 C++ 行数：约 8 行 LLL 满秩输入前置条件与 Gram--Schmidt norm 更新所依赖的数学边界
- 放弃的方案：只由完整整基行列式非零直接断言各 pivot 正；改为显式证明每个实际行前缀的 Gram 行列式严格为正

## Executable row-prefix recurrence

`sourceRowDot` now names the exact cast of the generated integer dot-product
loop at an explicit source row length.  `GramMuPrefixCorrect` records, for
each column already visited by `gramMuRowLoop`, the concrete LDL row equation
using only the stored `mu` and `norms` cells.  Its empty-prefix base is
immediate and contains no guessed coefficient array.

`sourceGramCoefficient_closes_column` proves the local algebraic step used by
the forthcoming loop induction.  When the previously generated norm is
nonzero, the exact source numerator/division expression written at column
`j` makes the next LDL row equation hold.  This theorem retains the source
zero-norm branch in `sourceGramCoefficient`; positivity established by the
processed-state invariant is what rules that branch out, rather than erasing
it from the generated execution.

## 度量（row-prefix recurrence）

- 耗时：约 0.5 小时（可执行逐列不变量与局部除法代数）
- 迭代：约 2 轮编译-修复循环
- Lean 新增/修改行数：约 47 行
- 对应 C++ 行数：约 6 行 Gram--Schmidt numerator/division 更新
- 下一步：沿实际 `j < i` 良基循环保持该不变量，再由 norm loop 闭合对角项

## Generated μ-row recurrence completed

The executable row invariant is now preserved by the actual generated
`gramMuRowLoop`.  `gramMuPrefixCorrect_set_next` proves the precise mutable
array step: the nested source writes replace only `mu[i][j]`, preserve every
earlier target-row cell and every non-target row, and extend the LDL row
equation by one column using the generated coefficient.

`gramMuRowLoop_prefix_correct` then follows the source recursion with measure
`i - j`.  At each recursive call it expands the real integer dot loop and
rational numerator loop, proves the written conditional quotient is exactly
`sourceGramCoefficient`, rules out the source zero-norm branch from the
already processed positive norms, and invokes the one-column update theorem.
The terminal branch restricts the accumulated prefix invariant to all
columns below `i`.  Thus the entire returned μ row is now characterized by
the generated execution, not by an assumed Gram--Schmidt result.

## 度量（complete μ row）

- 耗时：约 1.5 小时（嵌套数组写帧、逐列递归、dot/numerator 精确值复用）
- 迭代：约 6 轮编译-修复循环
- Lean 新增/修改行数：约 205 行
- 对应 C++ 行数：约 11 行完整 `j < i` μ-row 循环
- 下一步：证明实际 `gramNormLoop` 写入闭合 `(i,i)` 对角方程，并将整行方程组装为 `G = L D L^T`

## Generated diagonal norm equation

The diagonal computation is now identified with the actual generated
`gramNormLoop`.  `gramNormLoop_closes_diagonal` starts from the concrete
integer self-dot product returned by `dotRows`, applies the exact
well-founded norm-loop theorem, and reconciles every bounded regular source
read with the total reads used by the semantic invariant.  Its conclusion is
the exact diagonal LDL equation

`dot(b_i,b_i) = Σ k<i, mu[i,k]^2 * norm[k] + norm[i]`.

No positivity or Gram--Schmidt output is assumed by this equation.  The proof
only identifies the scalar returned by the generated subtraction loop.
`gsLowerNormMul_apply` additionally exposes every entry of the stored
`L * D * L^T` product as the corresponding finite sum, preparing the next
matrix-extensional row-extension proof.

## 度量（generated diagonal）

- 耗时：约 0.75 小时（norm-loop 反演、数组读一致性与矩阵乘积展开）
- 迭代：约 3 轮编译-修复循环
- Lean 新增/修改行数：约 75 行
- 对应 C++ 行数：约 7 行 `gramNormLoop` 与 norm 写回前的标量计算
- 下一步：按旧前缀、非对角新行和新对角三类矩阵元素组装 `G = L D L^T`

## One-row LDL matrix extension

The scalar equations are now assembled into the next full matrix prefix.
Exact μ/norm frame lemmas show that a generated row iteration leaves every
already processed `L` and `D` prefix entry unchanged.  The finite-sum support
lemmas then prove three concrete facts about the stored lower-triangular
matrix: an old entry receives no contribution from the new last column, a new
off-diagonal row entry reduces exactly to the generated row recurrence, and
the new diagonal reduces exactly to the generated norm recurrence.

`concreteGramSchmidtUpTo_extend_one` combines these facts extensionally.  It
handles old-prefix entries using the prior executable invariant, the new row
using `GramMuPrefixCorrect`, the new column by symmetry of the actual integer
Gram sum, and the final diagonal using the exact norm-loop equation.  The
result is `ConcreteGramSchmidtUpTo state (i + 1)`, i.e. the complete next
`G = L * D * L^T` prefix for the arrays produced by execution.  No alternate
orthogonalization, choice function, or existence certificate is involved.

## 度量（one-row LDL extension）

- 耗时：约 2 小时（前缀帧、有限和支撑、四类矩阵元素与依赖索引）
- 迭代：约 11 轮编译-修复循环
- Lean 新增/修改行数：约 280 行
- 对应 C++ 行数：约 18 行 μ/norm 两个内层循环及外层单行写回
- 下一步：将该扩展定理嵌入生成的外层初始化循环，并由前缀 Gram 行列式正性推出每个新 pivot 严格为正

## Generated pivot positivity

The row-extension invariant now also yields positivity of the newly written
Gram--Schmidt norm.  `ConcreteGramSchmidtUpTo.gram_prefix` derives the exact
prefix determinant/product identity from any processed executable prefix,
not only from a completed LLL state.  `generatedGramPivot_positive` combines
that identity at prefix `i+1` with the previously proved strict positivity of
the actual prefix Gram determinant.

Expanding `prefixNormProduct (i+1)` gives the product of all prior generated
norms times the newly stored norm.  The processed invariant makes the prior
product strictly positive, while the Gram determinant makes the total product
strictly positive; ordered-field arithmetic therefore forces the new norm to
be strictly positive.  This derives positivity after execution and does not
assume a nonzero divisor branch or a precomputed orthogonal basis.

The supporting prefix-frame and finite-sum lemmas were compiled together with
the complete one-row LDL extension, so the next outer-loop theorem can carry
both matrix semantics and positivity through each concrete source iteration.

## 度量（pivot positivity）

- 耗时：约 0.75 小时（前缀 determinant/product 桥与严格正性代数）
- 迭代：约 3 轮编译-修复循环
- Lean 新增/修改行数：约 65 行（含通用 processed-prefix determinant 引理）
- 对应 C++ 行数：不新增执行代码；闭合 norm 写回所需的满秩输入数学义务
- 下一步：沿 `size - i` 外层良基递归组合 row、diagonal、LDL 与 positivity，并删除 `LLLInitializationCorrect`

## Generated LLL initialization completed

The actual outer initialization recursion is now fully validated.
`processedGramSchmidtValid_step` combines one successful generated μ-row,
self-dot, norm loop, and norm array write.  It transports the old prefix
through the exact array frames, extends `G = L * D * L^T` by one row, derives
the new pivot's positivity from the full-rank input, and returns the complete
processed-state invariant at `i+1`.

`initializeGramSchmidtLoop_processed_valid` follows the generated outer loop
with its source measure `matrix.size - i`.  Its recursive branch invokes the
single-step theorem on the exact values returned by both inner loops; its
terminal branch proves the source index has reached the matrix size.  Hence a
successful generated outer execution carries exact LDL semantics and strict
positivity for every row.

Finally, `initializeLLL_concrete_valid` directly inverts the generated
`initializeLLL`: it follows the actual `makeInitialMatrix` result, exact row-0
dot product, seeded norm write, and generated outer loop.  The resulting
arrays satisfy `ConcreteLLLExecutionValid`; the transform value is precisely
the one returned by source execution and is not replaced or chosen by the
proof.

The former proposition `LLLInitializationCorrect` and the corresponding
argument to `concreteLLLExecution` have been deleted.  The concrete LLL
execution object now supplies its `initialized_valid` field exclusively from
`initializeLLL_concrete_valid`, eliminating the last abstract initialization
assumption.

## 度量（complete generated initialization）

- 耗时：约 2.5 小时（单步组合、外层良基递归、入口反演与首行索引换元）
- 迭代：约 12 轮编译-修复循环
- Lean 新增/修改行数：约 260 行
- 对应 C++ 行数：约 25 行完整 Gram--Schmidt/LLL 初始化控制流
- 删除的抽象边界：`LLLInitializationCorrect` Prop 参数及其注入点
- 下一步：将无参数 `concreteLLLExecution` 接入 concrete CLD/van-Hoeij 状态递归，闭合最终 recombination 合同

## Candidate-validation boundary made concrete

The candidate product path no longer exposes the intermediate proposition
alias `TrialProductStepCorrect`.  Although that proposition had already been
proved solely from the generated `multiplyNormalizeModRaw` execution, its
shape resembled an injectable correctness contract.  The replacement theorem
now quantifies directly over a successful raw execution and proves its exact
polynomial product modulo the concrete modulus.

The concrete candidate-validation dependency value has also been defined.
Both generated operation records are data-free: they contain neither a
callback nor candidate/product data, so the only possible execution is the
generated candidate availability, modular multiplication, symmetric recovery,
primitive normalization, exact division, and consumed-bit mutation chain.

## 度量（candidate validation boundary）

- 耗时：约 0.25 小时（抽象边界复查、直接定理改写与编译）
- 迭代：1 轮 Lean 编译
- Lean 净变化：删除 1 个 Prop 合同，增加 1 个无回调 concrete ops 值
- C++ 变化：无
- 下一步：证明 CLD lattice extension/reset 保持 concrete LLL 满秩输入，并构造完整 `VanHoeijRawOps`

## Van-Hoeij lattice dimension invariant exposed

The old generated `cld_extension_valid` field required every successful CLD
extension to preserve a square full-rank LLL input, but did not state the
source precondition `matrix.size = cld.size + currentColumns`.  Consequently
the field quantified over dimension-mismatched arrays for which the generated
append operation is not square; no honest concrete implementation of that
contract could exist.

The generated interface and its generator now carry the actual C++ lattice
dimension invariant.  `VanHoeijStateValid` combines concrete LLL full rank
with `matrix.size = active.size + currentColumns`.  Successful gather and CLD
executions carry exact output-size facts, CLD extension returns both full rank
and its new dimension, reset returns the exact factor-count dimension, and
LLL reduction now proves it preserves matrix size.  Both recursive branches
of `vanHoeijLoop` transport this invariant through their actual generated
outputs.

On the refinement side, the gather and complete nested CLD loops now have
direct successful-execution size theorems.  The LLL main loop also has a
well-founded execution proof that every advanced or swapped step preserves
matrix size.  Finally, `det_append_unit_lower` establishes the determinant
identity for the block-lower-triangular matrix shape produced by one CLD
column append; it is the algebraic core of the next full-rank extension proof.

## 度量（van-Hoeij dimension invariant）

- 耗时：约 1.25 小时（识别不可满足接口、生成器同步、递归 size 证明与 block determinant 引理）
- 迭代：约 5 轮生成/Lean 编译
- 生成层变化：状态新增真实 lattice dimension invariant；LLL 输出新增 size preservation
- 删除的退化风险：不再允许用过宽且事实上不可构造的 `cld_extension_valid` 掩盖维度前提
- C++ 变化：无
- 下一步：把 `appendZeroColumn`/`appendCldColumn` 的实际数组条目映射到 block determinant 引理，完成 CLD full-rank preservation

## Concrete CLD extension and lattice reset validated

The generated CLD-extension path is now proved against its exact array
execution.  The proof follows `appendZeroColumnLoop`,
`fillCldDataRowLoop`, the identity-coordinate write, and the final row push.
It identifies the resulting matrix as a unit lower block extension and applies
the block determinant theorem, so every successful `appendCldColumn` preserves
nonzero determinant.  A well-founded proof over `target - added` lifts this
fact through the complete generated `buildCldMatrixLoop`, while also retaining
the exact source dimension.

The retry reset has likewise been discharged from execution rather than an
assumption.  `InitialMatrixPrefix` tracks the exact entries written by
`setInitialDiagonalLoop`; induction on `size - index` proves that
`makeInitialMatrix size scale` returns the scaled identity matrix.  Its
determinant is nonzero whenever `scale` is nonzero.  In particular,
`resetVanHoeijLattice` uses `scale = 2 ^ vanHoeijExponent`, so its returned
matrix is a concrete full-rank LLL input with exactly one row per active
factor.

## 度量（CLD extension and reset）

- 耗时：约 2.0 小时（数组轨迹、block determinant、缩放单位阵与 codegen 检查）
- 迭代：约 14 轮 Lean 编译/构建修正
- Lean 新增：约 620 行直接执行精化与数组引理
- C++ 变化：无，因此本次无新的 C++ b2b 变更面
- 验证：`lake build CLPoly.Refinement.Recombine`、生成器 stale check、strict refinement boundary 和 `git diff --check` 通过
- 下一步：从模多项式长除的实际次数下降构造 `DivmodTermination`，然后组装无抽象字段的 `VanHoeijRawOps`

## Modular-divmod representation invariant established

The remaining `DivmodTermination.trace` field was audited before attempting
to instantiate it.  Its old universal domain includes arbitrary arrays, but
the C++ `upolynomial_` loop relies on its representation invariant: stored
terms have nonzero coefficients and strictly descending degrees.  For a
malformed divisor whose tail degree exceeds its head, the source update need
not decrease the remainder degree, so a universal trace would be an invalid
termination oracle rather than a proof of the C++ execution.

The refinement layer now begins the concrete repair from the actual sparse
representation.  `modCoeffOutput_canonical` proves that the generated
coefficient-reduction `filterMap` preserves strict degree order and removes
all zero coefficients.  `eraseLeading_canonical` and
`eraseLeading_size_lt` validate the exact zero-coefficient branch, while
`canonical_tail_degree_lt_head` extracts the strict head/tail inequality
needed for the nonzero merge branch.  These lemmas use the existing public
`SparsePolyZZCanonical` invariant rather than introducing a second weaker
notion of validity.

## 度量（divmod canonical foundation）

- 耗时：约 1.0 小时（终止接口反例审计、filterMap/erase 规范性证明）
- Lean 新增：约 80 行
- C++ 变化：无
- 验证：单文件 Lean 内核检查与 `lake build CLPoly.Refinement.Hensel` 通过
- 下一步：对 `divmodMergeLoop` 的六个实际游标分支证明严格降序保持与新首项次数下降

## Exact divmod merge preserves canonical sparse order

The nonzero-coefficient remainder update is now validated branch for branch.
`DivmodMergeCursorValid` records that the already emitted array is canonical
and that every emitted degree is above both unprocessed cursor suffixes (with
the divisor suffix shifted exactly as in C++).  Its common push theorem covers
both concrete outcomes of `pushNonzero`: retaining the array when the reduced
coefficient vanishes and appending a new nonzero term otherwise.

Three cursor transitions correspond to the source merge actions: consume the
remainder cursor, consume the shifted-divisor cursor, or consume both equal
degrees.  `divmodMergeLoop_canonical` invokes those transitions through all
six branches of the generated well-founded iterator, including each exhausted
cursor case and the terminal case.  Consequently the actual generated
`divmodRemainder` preserves `SparsePolyZZCanonical`; no list-sort operation or
mathematical polynomial remainder replaces the source execution.

## 度量（exact divmod merge canonicality）

- 耗时：约 2.0 小时（游标不变量、Bool `pushNonzero` 分支、六分支生成归纳器）
- Lean 新增：约 290 行
- C++ 变化：无
- 验证：单文件 Lean 内核检查与 `lake build CLPoly.Refinement.Hensel` 通过
- 下一步：在同一游标归纳上证明所有输出次数低于被消去的首项，构造余式次数良基 trace

## Concrete well-founded modular-divmod trace constructed

The exact merge proof now also carries an arbitrary strict degree bound over
the two live cursor suffixes.  Instantiating that bound with the current
remainder head, and using the exact signed shift identity, proves every term
emitted by `divmodRemainder` is strictly below the head that the C++ step
cancels.  This yields a source-state rank: zero for the empty remainder and
`headDegree + 1` otherwise.  Both concrete source branches strictly decrease
it: erasing a vanished leading coefficient and executing the complete
nonzero merge.

`concreteDivmodTrace` is consequently a genuine well-founded constructor for
the generated `DivmodTrace`.  It tests the exact source while condition and
coefficient, emits the corresponding `done`, `vanished`, or `subtract`
constructor, and recurses only on the actual updated remainder.  Canonicality
and machine-degree bounds are transported to the recursive state; they do not
supply a quotient, remainder, or semantic result.  The recursion uses
`termination_by divmodDegreeRank r`, with both decrease obligations discharged
by the concrete array proofs above.

## 度量（concrete modular-divmod trace）

- 耗时：约 2.0 小时（后缀 degree bound、两条 rank 下降、trace 良基构造）
- Lean 新增：约 315 行
- 度量：空余式 `0`，非空余式 `headDegree + 1`；不是 fuel 或执行上限
- C++ 变化：无
- 验证：单文件 Lean 内核检查与 `lake build CLPoly.Refinement.Hensel` 通过
- 下一步：收紧生成 `DivmodTermination` 的过宽全称域，将规范输入唯一接到 `concreteDivmodTrace`，删除可注入 trace 的抽象边界

## Concrete modular-divmod termination boundary

The generated division entry no longer asks a termination provider for a
trace on every arbitrary array.  `DivmodTermination` now exposes an explicit,
decidable source-input domain, and the raw entry faults before entering the
loop when that domain is violated.  This models the representation contract
of C++ `upolynomial_` instead of postulating termination for malformed sparse
arrays on which the source degree update need not decrease.

The sole production provider is `concreteDivmodTermination`.  Its domain is
checked by executable Boolean traversals: the divisor must be nonempty, both
arrays must have nonzero stored coefficients in strictly descending degree
order, and every degree must fit the signed 64-bit source arithmetic.  The
checker is proved equivalent to the existing canonicality and degree-bound
predicates.  Its trace field is exactly `concreteDivmodTrace`, so it contains
no semantic quotient/remainder witness and uses the previously proved
well-founded head-degree rank.

Both generated Hensel phases now transport this validity fact through their
actual raw execution.  Successful modular division implies validity, which
lets the semantic `_of_run` theorem recover the fact rather than demand an
independent oracle.  The centralized public Hensel theorems and `Pipeline/L1`
have also been specialized to the concrete provider; callers can no longer
inject an arbitrary `DivmodTermination` at those boundaries.

## 度量（concrete divmod boundary）

- 耗时：约 3.0 小时（生成接口、可执行 checker、phase 传播与 codegen 构建）
- 生成层变化：`DivmodTermination` 增加可判定输入域，raw entry 显式拒绝非法表示
- 删除的退化风险：公共 Hensel/流水线边界不再接受任意 trace provider
- C++ 变化：无，因此本次无新的 C++ b2b 变更面
- 验证：生成器 stale check、strict refinement boundary、`git diff --check`、
  `lake build CLPoly.Refinement.Hensel` 与 `lake build CLPoly.Refinement.Recombine` 通过
- 已知无关构建问题：公共层联合构建在未修改的 `SelectPrime.lean:823` 触发默认
  200k heartbeat；需以提高 heartbeat 的单文件检查复核本次公共 wrapper 修改
- 下一步：从 Hensel tree 的 canonical 不变量直接构造 concrete entry invariant，
  再接入 FactorZZ 的 Hensel/recombine raw ops

## Concrete van-Hoeij entry assembled

The generated recombination layer now includes the missing source entry
`__vanhoeij_recombine_raw_ir`.  It computes the source factor count, degree,
maximum CLD precision and initial precision, constructs the scaled diagonal
lattice through `resetVanHoeijLattice`, initializes the active indices, and
calls the existing lexicographically well-founded `vanHoeijLoop`.  After the
loop it executes the same remaining-factor append and degree sort as the C++
function via `finishZassenhaus`.

`concreteVanHoeijRawOps` fixes every dependency of that entry: concrete CLD
division, determinant-ranked LLL, proved CLD matrix extension, scaled-identity
reset, generated candidate validation, and well-founded Zassenhaus fallback.
`concreteVanHoeijTermination` is derived from the actual consumed-bit output
and reverse-erasure loop, so successful extraction decreases the active array
and cannot be supplied by a caller.

## 度量（concrete van-Hoeij entry）

- 耗时：约 0.75 小时（入口恢复、初始化证明、concrete ops 组装与 codegen）
- 生成层变化：新增完整 `__vanhoeij_recombine_raw_ir` 入口及生成器同步
- 良基度量：主循环仍为 `(active.size, precisionRank ...)` 的字典序；fallback
  使用 combination rank，LLL 使用 determinant/index rank
- C++ 变化：无，因此本次无新的 C++ b2b 变更面
- 验证：recombine 生成器 stale check、生成模块 codegen、Refinement/Recombine
  无限 heartbeat 内核检查通过
- 下一步：建立 van-Hoeij loop 的 live-product invariant并证明入口输出满足
  `RecombineCorrect`，随后接入 `FactorZZRawOps`

## Generated heuristic precision is now concrete and well-founded

The strict `FactorZZRawOps` boundary no longer accepts an injectable
`heuristicStartingPrecision` implementation.  The generated FactorZZ module
now executes the same floating-point FLINT estimate, the generated concrete
Mignotte-bound routine, and the source power loop that finds the least
exponent with `p^a > 2 |lc(f)| B_mig`.

The power loop is total under the source prime precondition `p ≥ 2`.  Its
well-founded measure is the remaining natural-number distance
`target + 1 - pa`; on the source loop branch `pa ≤ target`, multiplication by
the prime strictly increases `pa`.  Invalid non-prime inputs and an exponent
that cannot fit the source 32-bit result are reported as raw execution faults,
rather than being hidden by fuel, `partial`, or an assumed termination
callback.

## 度量（concrete heuristic precision）

- 耗时：约 0.5 小时（生成入口、乘方循环度量、机器整数边界）
- 删除的退化风险：`FactorZZRawOps` 不再允许调用者替换启发式精度执行
- 良基度量：`target + 1 - pa`，递归分支由 `p ≥ 2` 与 `pa > 0` 严格下降
- C++ 变化：无，因此本次无新的 C++ b2b 变更面
- 验证：FactorZZ 生成器 stale check、生成 Lean 模块内核检查和 diff check 通过
- 数学审计发现：严格重组目前仅闭合 trial-division 的乘积保持；最终
  `RecombineCorrect` 所需不可约性仍必须从满 Mignotte 精度、模素不可约因子和
  Zassenhaus 最小子集穷举推出，不能复用算法层接收 `h_irred` 参数的包装定理

## Concrete Hensel output exceeds its computed precision target

`HenselLiftLoopCorrect.outputM_gt_target` now proves that the modulus returned
by the generated quadratic Hensel loop is strictly larger than its input
target.  The proof follows the exact semantic trace: the terminal constructor
carries the negated source test `¬m ≤ target`, while recursive constructors
transport the terminal modulus unchanged.  No independently selected modulus
or precision certificate is involved.

The entry-level theorem
`HenselLiftEntryCorrect.outputModulus_gt_target` connects that loop fact to the
actual pair returned by `__hensel_lift_upoly_raw_ir` and to the target selected
by the generated Mignotte/explicit-precision branch.  This supplies the first
required premise for symmetric factor recovery in strict recombination.

## 度量（Hensel reached precision）

- 耗时：约 0.75 小时（循环后置条件、入口 witness 解包、整数转换）
- 新语义证据：实际 `outputM > target`，而不是仅有执行终止
- C++ 变化：无，因此本次无新的 C++ b2b 变更面
- 验证：`Hensel.lean` 以无限 heartbeat 完整内核检查通过
- 下一步：从同一个具体 Hensel trace 导出提升叶因子逐项模 p 还原到输入的关系

## Hensel tree algebra reduces through every concrete lift round

`HenselNodeReduces` and `HenselArrayReduces` now track all four algebraic
fields (`g`, `h`, `s`, `t`) of every concrete tree node.  The array relation
has proved reflexive/transitive rules, an exact `Array.set!` update rule, and
a modulus-projection rule along divisibility.

`HenselLiftRecursiveCorrect.arrayReduces` follows all four generated recursive
tree shapes.  At the current node it uses the already proved semantic result
of the actual `__hensel_step_raw_ir`; it then composes the exact left and right
recursive executions.  `HenselLiftLoopCorrect.arrayReduces_of_dvd` carries the
relation through every quadratic precision round and projects each round back
to any divisor of the initial modulus.  Taking that divisor to be the selected
prime gives the array-level fact needed to relate final lifted leaves to their
original finite-field factors.

## 度量（Hensel modular ancestry）

- 耗时：约 1.25 小时（节点/数组关系、set 更新、模数投影、递归和外层循环归纳）
- 覆盖执行：实际节点更新以及实际左/右子树调用，不是列表级替代算法
- C++ 变化：无，因此本次无新的 C++ b2b 变更面
- 验证：`Hensel.lean` 以无限 heartbeat 完整内核检查通过
- 下一步：组合初始树叶布局、最终提取遍历与 normalization，导出输出数组逐项模 p 对应

## Generated leaf extraction order is now explicit

`henselExtractedFactors` gives a pure denotation of the exact generated leaf
walk: a missing left/right child emits the current node's `g`/`h`, while a
present child recursively emits that subtree.  It therefore follows the same
left-to-right source order and reads the actual node array; it does not accept
an expected factor list.

`HenselExtractCorrect.toList_eq` proves that every concrete extraction trace
returns precisely its original prefix followed by this denotation.
`HenselExtractCertificate.extractedFactors_forall₂` then combines the builder's
constructor-shaped bounds certificate with `HenselArrayReduces`: reading the
same topology from pre-lift and post-lift arrays yields pointwise factors that
are equal modulo the selected prime.

## 度量（exact Hensel leaf extraction）

- 耗时：约 1.0 小时（提取 denotation、四分支 trace 归纳、逐项模 p 关系）
- 覆盖执行：实际 `g/h` 叶分支和实际左右递归顺序
- C++ 变化：无，因此本次无新的 C++ b2b 变更面
- 验证：`Hensel.lean` 以无限 heartbeat 完整内核检查通过
- 新发现的剩余边界：builder 当前只公开根区间乘积；还需沿实际递归 builder
  导出所有初始叶子按顺序对应 adjusted 输入因子，不能由根证书替代

## Recursive builder semantic certificate introduced

`HenselTreeSemanticBuildCertificate` is now the algebraic companion to the
existing raw topology certificate.  At every concrete node it records the
proved `HenselTreeNodeGCDInvariant` for that node's exact `[start, stop)` input
interval, and its child constructors recurse over the same midpoint halves as
the generated builder.

The certificate supports lowering its allocation-index boundary and transport
through later disjoint array updates.  Transport uses the builder's exact
`HenselTreePreservesFrom` equations, so the stored interval semantics remain
tied to the same concrete nodes.  The certificate contains neither a result
factor list nor a recombination witness.

## 度量（builder semantic certificate）

- 耗时：约 0.75 小时（递归证书、索引单调性、数组 frame 运输）
- 新证据形状：每个实际节点对应其实际输入半开区间乘积
- C++ 变化：无，因此本次无新的 C++ b2b 变更面
- 验证：`Hensel.lean` 以无限 heartbeat 完整内核检查通过
- 下一步：让现有 `strictHenselTreeBuildRecursiveRawIR_succeeds` 在四种子树分支中
  同时产出该证书，再归纳得到叶子与 adjusted 输入逐项对应

## Actual Hensel builder now emits the semantic trace

`strictHenselTreeBuildRecursiveRawIR_succeeds` now constructs a
`HenselTreeSemanticBuildCertificate` along the generated builder execution.
All four source shapes are covered: both recursive children, only the left
child, only the right child, and two direct leaves.  Recursive certificates
are transported through the actual intervening array writes using their
proved frames, so every retained node still denotes its concrete input
interval in the final builder array.

`strictHenselTreeBuildRawIR_refines_topology_root` exposes this certificate for
the canonical root topology, and `HenselLiftEntryCorrect` retains it beside
the real lift-loop, extraction, and normalization traces.  Thus downstream
proofs can recover the modular origin of each output factor without adding an
existence witness or semantic oracle.  The leaf extraction helper also now
uses explicit left/right `nodeCount` decrease lemmas and an explicit bounded
array-read bridge.

## 度量（builder semantic trace emission）

- 耗时：约 1.25 小时（四分支组合、frame 运输、入口证书接线、完整检查）
- 新证据形状：实际 builder 输出数组中每个节点保留对应输入区间语义
- 良基性：叶提取分别使用左右子树 `nodeCount` 严格下降定理，无 fuel/partial
- C++ 变化：无，因此本次无新的 C++ b2b 变更面
- 验证：`Hensel.lean` 以无限 heartbeat 完整内核检查通过
- 下一步：在 canonical topology 上归纳，证明提取叶列表与 adjusted 模因子列表
  逐项对应，再穿过 lift 与 normalization

## Ordered source interval list for leaf correspondence

`henselFactorRangeList` now denotes the exact source-array order over the same
half-open interval used by the generated Hensel product loop.  Its recursion
uses the concrete `index + 1` progression and the well-founded measure
`stop - index`.  `henselFactorRangeList_split` proves that an interval is the
ordered concatenation of its two midpoint halves, and the singleton theorem
identifies a one-element interval with the actual array read.

These lemmas provide the input side of the forthcoming pointwise `Forall₂`
theorem between modular factors and concrete builder leaves; they do not
introduce lifted outputs or expected factor values.

## 度量（ordered Hensel input intervals）

- 耗时：约 0.5 小时（良基列表 denotation、区间拆分、singleton 边界）
- C++ 变化：无，因此本次无新的 C++ b2b 变更面
- 验证：`Hensel.lean` 以无限 heartbeat 完整内核检查通过
- 下一步：组合 canonical topology 与 builder semantic certificate，得到初始叶子
  和 adjusted 输入的逐项对应

## Canonical builder leaves have exact modular origins

`HenselTreeSemanticBuildCertificate.extractedFactors_forall₂` now proves the
pointwise, ordered relation between the adjusted finite-field factor interval
and the leaves read from the actual canonical builder array.  The proof
follows all four concrete topology shapes.  Recursive sides consume the child
certificate; direct-leaf sides reduce the node's proved interval product to
the unique source array element using
`henselFactorRangeProduct_singleton`.

The theorem keeps the certificate allocation lower bound independent from
each subtree root.  This is essential because child certificates retain the
original builder boundary while their actual node indices increase.  No
expected lifted list, permutation witness, or irreducibility oracle is an
input to the theorem.

## 度量（exact modular leaf origins）

- 耗时：约 1.0 小时（四拓扑分支、singleton product、实际数组读取桥接）
- 结论强度：按源顺序逐项 `Forall₂`，不是仅证明总乘积或集合相等
- C++ 变化：无，因此本次无新的 C++ b2b 变更面
- 验证：`Hensel.lean` 以无限 heartbeat 完整内核检查通过
- 下一步：与 lift-loop 的 `HenselArrayReduces` 及精确 extraction trace 合成，
  得到 normalization 前的最终提升因子逐项模 p 来源

## Lifted extracted factors retain their modular origins

The builder semantic certificate now transports its concrete leaf reads
through `HenselArrayReduces`, using every actual node position and the proved
array-size equality.  `HenselLiftEntryCorrect.preNormalizationOrigins`
composes this transport with the canonical initial leaf-origin theorem and
`HenselExtractCorrect.toList_eq`.  Its conclusion is therefore about the
actual `extracted.toList` produced immediately before the source normalization
block.

`HenselLiftEntryCorrect` now also exposes its genuine prime-modulus dependency
as an explicit typeclass parameter.  This dependency was already required by
the stored builder certificate; making it explicit prevents later consumers
from treating a composite or unproved modulus as if it carried the same
finite-field semantics.

## 度量（pre-normalization modular origins）

- 耗时：约 0.75 小时（数组约化运输、两段 `Forall₂` 合成、extraction 接线）
- 覆盖执行：builder 初始数组 → 实际 lift loop → 实际 extraction 输出
- C++ 变化：无，因此本次无新的 C++ b2b 变更面
- 验证：`Hensel.lean` 以无限 heartbeat 完整内核检查通过
- 下一步：证明 normalization 只对首因子施加可逆单位缩放，并把模不可约性
  从 adjusted 输入运输到最终 Hensel 输出
