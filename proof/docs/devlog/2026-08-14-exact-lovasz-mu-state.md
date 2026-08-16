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

## Final Hensel normalization is pointwise unit scaling

`HenselLiftLoopCorrect.initialM_dvd_outputM` derives divisibility of the final
modulus from the concrete quadratic loop: each recursive source step changes
`m` only to `m * m`.  This makes the selected prime a divisor of the actual
normalization modulus without a separate precision assumption.

`HenselNormalizeCorrect.unitRel` follows every generated normalization branch.
The empty and already-normalized branches use unit one.  In the scaling branch,
the successful concrete `ZZ.invert` call yields an inverse in `ZMod outputM`;
the divisibility map projects it to a unit in `ZMod p`.  The proved
`__upoly_mod_coeff_raw_ir` divisor semantics then identifies the normalized
first factor with that unit scalar times the extracted first factor.  Every
other array position is proved unchanged through the actual bounded `set!`.

The entry-level origin theorem now returns this unit relation together with
the pre-normalization `Forall₂` fact from the same existential execution
trace, so downstream irreducibility transport cannot accidentally mix
different witnesses.

## 度量（Hensel normalization unit semantics）

- 耗时：约 1.25 小时（最终模数整除、逆元投影、首元素更新、入口合成）
- 结论强度：最终数组每个位置等于对应 extracted 因子的模 p 单位倍
- C++ 变化：无，因此本次无新的 C++ b2b 变更面
- 验证：`Hensel.lean` 以无限 heartbeat 完整内核检查通过
- 下一步：把 adjusted 首因子的单位关系与 Zp 输出不可约性合成，导出最终
  Hensel 数组每个模 p 像不可约

## Exact semantics of the first-factor coefficient map

`scaleZpCoeffs_toPoly` proves that the generated array `map` used by
`__hensel_adjust_first_factor_raw_ir` denotes multiplication by the concrete
field constant.  The proof traverses the actual term list, applies the
verified `Zp` multiplication semantics at every coefficient, and discharges
the machine-word product bound from `p * p ≤ UInt64.size` and the concrete
reduced-coefficient bounds.

This closes the executable scaling part of the adjustment.  Combined with the
already proved sparse normalization semantics, it will identify adjusted
factor zero as the leading-coefficient unit multiple of the actual selected
Zp factor; all other positions are unchanged by the source `set!`.

## 度量（adjust-first scaling semantics）

- 耗时：约 0.5 小时（term-map 归纳、Zp 乘法、word overflow 上界）
- C++ 变化：无，因此本次无新的 C++ b2b 变更面
- 验证：`Hensel.lean` 以无限 heartbeat 完整内核检查通过
- 下一步：结合 `GoodPrime.lc_nonzero` 和实际首项读取，证明 adjusted 因子
  保留不可约性

## First-factor adjustment is pointwise unit scaling

The concrete `Zp.ofInt` residue used by the adjustment now has proved
`Reduced` and exact `ZMod` decoding lemmas.  Using those facts,
`HenselAdjustFirstFactorCorrect.unitRel` follows the actual source read,
coefficient map, sparse normalization, and `set! 0`.  Under the selected
leading coefficient's nonzero premise, factor zero is multiplied by a field
unit and every other array position is unchanged.

Both `HenselAdjustUnitRel.irreducible` and
`HenselNormalizeUnitRel.irreducible` transport irreducibility pointwise through
their respective unit scalings.  They use the polynomial constant-unit
theorem and do not infer factor quality from a product identity.

## 度量（adjusted-factor unit and irreducibility transport）

- 耗时：约 1.0 小时（`Zp.ofInt` 解码、首位置语义、数组 frame、不可约性运输）
- 结论强度：adjusted 与最终 normalized 数组均按下标保留不可约性
- C++ 变化：无，因此本次无新的 C++ b2b 变更面
- 验证：`Hensel.lean` 以无限 heartbeat 完整内核检查通过
- 下一步：从严格 SelectPrime 的 `CandidateCorrect.quality` 构造原始 factors
  数组逐项不可约前提，并组合整条 Hensel 链

## Strict SelectPrime quality reaches adjusted Hensel factors

The new downstream module `Refinement/FactorZZ.lean` begins the genuine
cross-stage assembly without introducing an import cycle.  Its
`selectionFactors_irreducible` theorem decodes each concrete array position
from the successful `SelectionCorrect/CandidateCorrect.quality` certificate.
`selectionFactors_adjusted_irreducible` then supplies the actual selected
prime, machine-word bound, canonical sparse inputs, and source leading-term
semantics to the Hensel adjustment unit theorem.

The result states irreducibility of every concrete adjusted array element.
It neither reruns FactorZp nor replaces the selection array with an abstract
list.  A formal `lake build CLPoly.Refinement.Hensel` also exposed and removed
two stale-olean hazards: a forward helper reference and inaccessible names in
a dependent constructor pattern.  The importable Hensel artifact now builds,
and the new FactorZZ module checks against that artifact.

## 度量（SelectPrime → adjusted Hensel quality）

- 耗时：约 1.0 小时（跨模块组合、正式 olean 构建修复、逐数组元素质量）
- 新模块：`CLPoly.Refinement.FactorZZ`
- C++ 变化：无，因此本次无新的 C++ b2b 变更面
- 验证：正式 Hensel lake target 构建通过；`FactorZZ.lean` 内核检查通过
- 下一步：把 adjusted 不可约性与 entry-level leaf origins、最终 normalization
  unit relation合成，导出最终 lifted array 的逐项模 p 不可约性

## Selected irreducibility reaches the concrete Hensel output array

`henselFactorRangeList_suffix` identifies the leaf interval used by the actual
tree builder with the corresponding suffix of the adjusted source array.  In
particular, the full interval starting at zero is definitionally tied to
`adjusted.toList`; this removes the last list/array ordering ambiguity from the
entry-level origin certificate.

`selectionHenselFactors_mod_irreducible` now composes the successful concrete
SelectPrime certificate with the actual first-factor adjustment, builder leaf
order, lifted-node extraction, and final normalization.  A structural
`Forall₂` induction transports irreducibility across the exact pointwise
mod-p equalities, after which the already proved normalization unit relation
transports it to every position of the returned Hensel array.  No factor is
introduced existentially and no product-only argument is used as a substitute
for pointwise provenance.

## 度量（SelectPrime → final Hensel output quality）

- 耗时：约 1.0 小时（完整 interval/list 身份、`Forall₂` 运输、数组下标合成）
- 结论强度：实际 Hensel 返回数组的每个因子模选定素数后不可约
- C++ 变化：无，因此本次无新的 C++ b2b 变更面
- 验证：正式 `CLPoly.Refinement.Hensel` lake target 构建通过；
  `FactorZZ.lean` 内核检查通过；正式 FactorZZ target 首次重放还暴露出
  `SelectPrime.lean` 最终组合在默认 200k heartbeat 下超时，因此该大型真实
  精化模块已显式使用无限 heartbeat；随后 SelectPrime 与 FactorZZ 正式
  targets 均完整构建通过
- 下一步：把这项逐因子质量接入完整精度的真实 Zassenhaus/van Hoeij 重组，
  证明输出乘积与整数不可约性

## Concrete recombination finishing preserves factor products

The proof layer now decodes the common `finishZassenhaus` block used by both
the Zassenhaus fallback and the normal van-Hoeij exit.  The degree sort is
proved to preserve the polynomial product through the actual `mergeSort`
permutation.  The finishing theorem then records the exact remaining effect:
the live `fStar` is appended precisely when its source front degree is
positive, followed by the same concrete degree sort.

This is deliberately an exact execution fact rather than a semantic
factorization assumption.  It supplies the terminal case needed for the
source-shaped outer-loop product invariant; the constant `fStar` case will be
discharged from primitivity instead of silently dropping an arbitrary scalar.

## 度量（recombination finish product）

- 耗时：约 0.5 小时（排序 permutation、数组 push/product、分支精确化）
- C++ 变化：无，因此本次无新的 C++ b2b 变更面
- 验证：正式 `CLPoly.Refinement.Recombine` lake target 构建通过
- 下一步：以 live product 的 primitivity 为循环不变量，证明 Zassenhaus
  extraction/restart/exhaustion 的完整乘积保持

## The dropped constant remainder is proved to be a unit

The terminal product theorem is strengthened from a branch description to an
`Associated` result for the complete live product.  For a canonical sparse
integer polynomial, a concrete front degree of zero forces the backing array
to contain exactly one term: a second term would need a natural-number degree
strictly below zero.  Thus the source branch that omits `fStar` omits only a
constant polynomial.

The live-product primitivity invariant then forces that constant coefficient
to be a unit: taking polynomial content gives
`normalize coefficient * content(result product) = 1`.  The proof transports
the resulting integer unit through `Polynomial.C` and concludes association
with the emitted product.  The empty-`fStar` branch is shown impossible under
the same primitivity invariant.  No scalar is discarded by convention.

## 度量（terminal constant/unit bridge）

- 耗时：约 0.75 小时（canonical chain、常数识别、content/unit 运输）
- C++ 变化：无，因此本次无新的 C++ b2b 变更面
- 验证：正式 `CLPoly.Refinement.Recombine` lake target 构建通过
- 下一步：证明 successful Zassenhaus extraction 的 quotient primitive
  继续满足 canonical 表示，并闭合整个外层循环的不变量

## Primitive quotient normalization preserves canonical sparse form

`primitiveDivideLoop_toList` now identifies the exact backing list emitted by
the source coefficient-division loop: the existing accumulator followed by
the input suffix with every coefficient divided by the concrete divisor.
`primitiveDivideLoop_constraints` separately extracts the runtime facts that
the divisor is nonzero and divides every copied coefficient from the same
successful execution trace.

Using those facts, `primitiveRaw_canonical` proves that the actual
`__upoly_primitive` lowering preserves the strict descending degree chain and
cannot introduce a zero coefficient.  Degree order is preserved because the
loop copies each monomial unchanged; quotient nonzeroness follows by exact
division from the original canonical nonzero coefficient.  This applies
directly to the `quotientPrimitive` installed after successful Zassenhaus and
candidate-validation extraction.

## 度量（primitive quotient canonicality）

- 耗时：约 1.0 小时（精确输出列表、运行时整除约束、canonical 运输）
- C++ 变化：无，因此本次无新的 C++ b2b 变更面
- 验证：正式 `CLPoly.Refinement.Recombine` lake target 构建通过
- 下一步：把 quotient 的 canonicality、已有 product/unit-scalar 和
  primitivity 保持组合进 Zassenhaus 良基外层递归

## Exact-division remainders retain nonzero stored coefficients

The concrete sparse normalization used after each subtraction is now proved
to retain only nonzero coefficients.  The proof follows the actual grouping,
zero-filter, and merge-sort representation: merge-sort membership is carried
back through its permutation to the filtered array, whose Boolean predicate
decodes to coefficient nonzeroness.

`exactDivmodLoop_remainder_coefficients_nonzero` transports this property
through the real well-founded long-division recursion.  Terminal branches
return the current remainder unchanged; recursive branches use the concrete
`subtractScaledNormalize` result and the checked `divisionRank` decrease.
This supplies the nonzero leading coefficient needed to show that every
quotient term pushed by an exact division is itself nonzero.

## 度量（exact-division remainder nonzeroness）

- 耗时：约 0.75 小时（normalization filter、permutation membership、WF 递归）
- C++ 变化：无，因此本次无新的 C++ b2b 变更面
- 验证：正式 `CLPoly.Refinement.Recombine` lake target 构建通过
- 下一步：利用 remainder nonzero 与 rank decrease 证明 quotient 的 pushed
  degree shifts 严格递减，闭合 exact quotient canonicality

## Checked exact division returns a canonical quotient

`ExactDivmodQuotientInvariant` records the source accumulator property: the
stored quotient is canonical, and every existing quotient degree is strictly
above the next shift whenever another division step is available.  Pushing a
new term preserves the strict chain because the current shift is below the
entire existing prefix.

`exactDivmodLoop_quotient_canonical` follows every branch of the actual
well-founded long-division loop.  The pushed coefficient is nonzero by exact
division of the current nonzero remainder lead.  On recursion, the checked
`divisionRank` decrease makes the next remainder front degree—and hence its
next quotient shift—strictly smaller.  Terminal nondivisibility and degree
branches return the already canonical accumulator unchanged.

The raw-entry theorem initializes this invariant at the actual empty quotient
and proves that every successful `exactDivmodRaw` call returns a canonical
quotient.  Consequently the subsequent concrete `primitiveRaw` call can now
provide a canonical `quotientPrimitive` without any semantic replacement.

## 度量（exact quotient canonicality）

- 耗时：约 1.25 小时（accumulator invariant、push、rank/shift decrease、raw entry）
- C++ 变化：无，因此本次无新的 C++ b2b 变更面
- 验证：正式 `CLPoly.Refinement.Recombine` lake target 构建通过
- 下一步：在 Zassenhaus successful-extraction 分支组合 quotient canonical、
  primitive canonical、unit-scalar product 和 primitivity，闭合外层递归

## Canonical quotient provenance crosses attempt and combination scan

`zassenhausAttempt_extracted_quotient_trace` exposes the exact successful
`exactDivmodRaw` and following `primitiveRaw` calls that produced the quotient
returned by a successful source attempt.  The trace is recovered by following
all concrete pruning, index conversion, modular product, symmetric recovery,
primitive-factor, trial-division, empty-remainder, and quotient-primitive
branches; it introduces no alternative quotient witness.

The trace composes the exact quotient canonicality theorem with primitive-part
canonical preservation.  A second well-founded theorem transports this result
through `scanZassenhausCombinations`: rejected candidates advance via the
concrete combination rank, while the successful candidate returns the same
canonical quotient.  The Zassenhaus outer extraction branch can therefore
install its recursive `fStar` with a proved representation invariant.

## 度量（attempt/scan quotient provenance）

- 耗时：约 1.0 小时（完整 attempt trace、两段 canonical 合成、scan WF 运输）
- C++ 变化：无，因此本次无新的 C++ b2b 变更面
- 验证：正式 `CLPoly.Refinement.Recombine` lake target 构建通过
- 下一步：闭合 Zassenhaus outer loop 的 product、primitivity、canonical
  三项联合不变量

## Complete source-shaped Zassenhaus product refinement

The full generated Zassenhaus outer loop now preserves its live product up to
a unit.  The theorem recurses on exactly the source pair
`(active.size, active.size + 1 - subsetSize)`: exhausted scans advance the
subset size, while successful extraction uses the concrete removal proof to
strictly shrink the active lifted array and restart at subset size one.

In the extraction branch, whole-live-product primitivity yields primitivity of
the current `fStar`.  The scan's unit-scalar product trace, canonical returned
quotient, and pushed result factor are combined into the next live product;
content transport proves that product primitive before the recursive call.
At termination, the previously proved finishing theorem handles the positive
degree append and proves any omitted constant remainder is a unit.

`zassenhausRecombine_product_associated` closes the concrete entry, including
its zero/one lifted-factor shortcut and the full outer-loop branch.  It proves
association with the actual sorted output array, not with an abstract list of
factors.

## 度量（complete Zassenhaus product refinement）

- 耗时：约 1.5 小时（联合不变量、双分量 WF 递归、entry shortcuts）
- C++ 变化：无，因此本次无新的 C++ b2b 变更面
- 验证：正式 `CLPoly.Refinement.Recombine` lake target 构建通过
- 下一步：把该定理接入 van-Hoeij fallback，并沿 successful validation、
  precision retry 和 normal exit 闭合 van-Hoeij 主循环

## Candidate validation preserves the live canonical remainder

`validateCandidatesLoop_fStar_canonical` follows the full concrete validation
loop over candidate rows.  Empty, trivial, unavailable, rejected, and
nondividing candidates retain the current canonical `fStar`.  Every successful
candidate follows the actual modular product, symmetric recovery, primitive
factor, exact division, empty-remainder test, quotient primitive conversion,
and consumed-bit marking before recursion.

For each such extraction, exact division returns a canonical quotient and the
subsequent primitive conversion preserves it.  The well-founded candidate
index recursion therefore carries canonicality through arbitrarily many
successful candidates in one validation pass.  The entry theorem initializes
the actual replicated consumed array and exposes canonicality of the concrete
returned `fStar'` needed by the van-Hoeij state transition.

## 度量（validation live canonicality）

- 耗时：约 1.25 小时（完整 candidate loop、连续 extraction、entry composition）
- C++ 变化：无，因此本次无新的 C++ b2b 变更面
- 验证：正式 `CLPoly.Refinement.Recombine` lake target 构建通过
- 下一步：使用 validation 的 product/primitivity/canonical 三个结论闭合
  van-Hoeij 主循环的 extraction、retry、fallback 与 normal exit

## Complete well-founded van-Hoeij product refinement

`vanHoeijLoop_product_associated` now follows every branch of the generated
C++ main loop and preserves the complete live product up to a unit.  A
successful candidate-validation round uses the concrete unit-scalar product
trace, canonical quotient trace, and primitivity transport before recursively
resetting the actual lattice.  A round without extraction either advances the
bounded precision target or executes the already refined concrete Zassenhaus
fallback and appends its factors to the accumulated result.

The proof uses exactly the generated lexicographic termination measure
`(active.size, precisionRank target initial maximum)`: extraction strictly
shrinks the active array; retry strictly decreases the remaining precision
rank.  There is no fuel counter, partial definition, semantic oracle, or L2
replacement branch.  The theorem additionally returns canonicality and
primitivity of every nonempty live output remainder, which is needed to
justify the source finishing block rather than assuming that a dropped
degree-zero remainder is harmless.

`__vanhoeij_recombine_raw_ir_product_associated` lifts the loop result through
the complete generated C++ entry: reset lattice, construct the active index
array, choose the source starting precision, execute the well-founded loop,
conditionally append the remaining positive-degree factor, and sort by
degree.  Thus its actual returned sparse array has product associated to the
input primitive polynomial.

## 度量（complete van-Hoeij product refinement）

- 耗时：约 1.5 小时（双分量 WF 循环、validation 输出不变量、entry finish）
- C++ 变化：无，因此本次无新的 C++ b2b 变更面
- 验证：单文件 `lake env lean CLPoly/Refinement/Recombine.lean` 通过；
  正式 `lake build CLPoly.Refinement.Recombine` 通过（3531 jobs）
- 下一步：从 lifted Zp 因子的不可约来源与真实重组候选性质，
  闭合重组输出的 ZZ 不可约性，再接入 FactorZZ 顶层组合

## Primitive modular irreducibility bridge

`primitive_irreducible_of_irreducible_mod` supplies the non-monic reduction
criterion required by the concrete C++ recombination output.  If a primitive
integer polynomial has irreducible reduction at the selected prime and its
leading coefficient survives reduction, then it is irreducible over the
integers.

The proof analyzes an arbitrary actual product decomposition.  Survival of
the product leading coefficient forces both factor leading coefficients to
survive modulo the prime, so neither mapped factor can become constant by
degree loss.  Modular irreducibility makes one mapped factor a unit; its
preserved degree is therefore zero.  Since each divisor of the primitive
product is primitive, that constant integer factor is a unit.  No abstract
factorization witness or semantic recombination premise is introduced.

## 度量（primitive modular irreducibility）

- 耗时：约 0.5 小时（非首一 Gauss 桥、次数保持、primitive 常数单位）
- C++ 变化：无，因此本次无新的 C++ b2b 变更面
- 验证：`lake env lean CLPoly/Refinement/FactorZZ.lean` 通过；
  正式 `lake build CLPoly.Refinement.FactorZZ` 通过（3545 jobs）
- 下一步：证明真实 Zassenhaus 组合枚举无遗漏，并将 Hensel
  的模素数不可约来源输送到每个成功恢复的 primitive 块

## Concrete combination terminal-state completeness

The generated right-to-left pivot search now has an exact terminal theorem.
If it returns no pivot, every array position equals its unique maximal
fixed-size-combination digit `upper - count + position`.  The proof follows
the actual inspected-suffix recursion and establishes maximality one source
comparison at a time.

`nextCombination_false_is_final` lifts this through the complete generated
C++ helper.  Under the real size precondition, a `false` return preserves the
input array and can occur only at the final combination.  In particular, the
Zassenhaus scan cannot terminate at an interior subset.  This is one half of
the no-omission proof; the remaining half is to show that every `true` step is
the immediate successor among strictly increasing bounded arrays.

## 度量（combination terminal state）

- 耗时：约 0.5 小时（none-pivot WF 归纳、false 分支、最大组合唯一性）
- C++ 变化：无，因此本次无新的 C++ b2b 变更面
- 验证：单文件 `lake env lean CLPoly/Refinement/Recombine.lean` 通过；
  正式 `lake build CLPoly.Refinement.Recombine` 通过（3531 jobs）
- 下一步：证明 `resetCombinationSuffix` 构造严格递增的最小后缀，
  然后得到 `nextCombination` 的直接后继与 scan 全覆盖

## Exact minimal suffix of every combination successor

Two well-founded frame/value theorems now characterize the generated
`resetCombinationSuffix` recursion.  Cells before the next write position are
preserved by every later write, while each cell at or after that position is
exactly the prior reset value plus its distance.  At offset zero, every cell
right of the pivot is therefore the pivot value plus `index - pivot`.

`nextCombination_true_minimal_suffix` lifts these facts through the complete
generated C++ helper.  Every `true` result has a concrete pivot such that the
left prefix is unchanged, the pivot increases by exactly one, and the entire
right suffix is the unique consecutive minimum beginning at the new pivot.
Together with the already proved rightmost-pivot/maximal-old-suffix theorem,
this supplies both value halves of the immediate-successor argument; it does
not reason through an alternative list-combination implementation.

## 度量（combination minimal successor suffix）

- 耗时：约 0.75 小时（processed-prefix frame、suffix 精确值、true 入口组合）
- C++ 变化：无，因此本次无新的 C++ b2b 变更面
- 验证：单文件 `lake env lean CLPoly/Refinement/Recombine.lean` 通过；
  正式 `lake build CLPoly.Refinement.Recombine` 通过（3531 jobs）
- 下一步：定义固定大小严格递增组合的 lex 顺序，证明不存在
  介于当前数组与该最小后缀结果之间的合法组合

## Immediate legal successor and strict scan invariant

The concrete combination certificate no longer tracks only array size and a
per-position upper bound.  Its invariant is now the actual C++ subset
representation: fixed size, strictly increasing entries (in accumulated-gap
form), and every entry below `upper`.  The generated initial iota array
establishes this invariant, and every successful generated
`nextCombination` preserves it.

`nextCombination_true_no_legal_between` proves that no legal fixed-size
combination lies lexicographically between an input and its successful source
successor.  The proof distinguishes the concrete first differing positions.
A difference to the right of the generated pivot contradicts the old
rightmost-maximal suffix; a difference at the pivot is bounded below by the
new pivot and its minimal consecutive suffix; and a difference to the left
contradicts the unchanged prefix.  Thus the well-founded generated scan now
carries the strong legality invariant needed for a genuine no-omission
argument, rather than an unrelated list enumerator or an assumed coverage
oracle.

## 度量（immediate legal successor）

- 耗时：约 1.5 小时（合法组合不变量、位置上界、直接后继、闭包）
- C++ 变化：无，因此本次无新的 C++ b2b 变更面
- 验证：单文件 `lake env lean CLPoly/Refinement/Recombine.lean` 通过；
  正式 `lake build CLPoly.Refinement.Recombine` 通过（3531 jobs）
- 下一步：从初始最小组合和直接后继推出 generated scan 对每个合法
  固定大小候选的可达性，再连接 rejected/extracted 结果

## Generated scan exhaustion has no omitted subset

The concrete first-difference array order now has a proved trichotomy.  The
source iota state is minimal among legal fixed-size combinations, while the
state returned by a `false` successor is maximal.  These boundary facts and
the immediate-successor theorem are consumed by a well-founded induction on
the generated scan's existing rank.

`scanZassenhausCombinations_exhausted_rejects_legal` follows the actual
`zassenhausAttempt`/`nextCombination` recursion.  At the target state it
extracts the observed rejected result directly.  Before the target, a false
successor contradicts maximality; a true successor cannot jump over the
target and recurses with the generated rank decrease.  Its initial-state
corollary `scanZassenhausCombinations_exhausted_rejects_all` therefore proves
that `.exhausted` means every legal fixed-size subset was actually attempted
and rejected.  No alternative enumerator, fuel bound, or semantic coverage
premise is used.

## 度量（generated scan no omission）

- 耗时：约 1.5 小时（首差三歧性、初始/终态边界、scan 良基归纳）
- C++ 变化：无，因此本次无新的 C++ b2b 变更面
- 验证：单文件 `lake env lean CLPoly/Refinement/Recombine.lean` 通过；
  正式 `lake build CLPoly.Refinement.Recombine` 通过（3531 jobs）
- 下一步：把所有固定大小 scan 的拒绝结论穿过 Zassenhaus outer loop，
  证明 finish 留下的 primitive 块不存在任何由 Hensel 原子组成的真因子

## Hensel normalization preserves pairwise coprimality

The exact pointwise unit relation already emitted by the generated Hensel
normalization proof now transports `IsCoprime` for any two output positions.
The proof constructs the transformed Bézout identity using the inverses of
the two concrete unit scales.  Consequently the final normalization cannot
merge, duplicate, or create a common modular factor between distinct Hensel
leaves.

This is a necessary uniqueness premise for converting a nontrivial integer
divisor into a unique subset of the irreducible modular leaves.  It adds no
new certificate field or semantic assumption: both scales and their unit
proofs come from the already executed normalization relation.

## 度量（Hensel pairwise normalization）

- 耗时：约 1 小时（单位缩放 Bézout 恒等式、数组位置传递）
- C++ 变化：无，因此本次无新的 C++ b2b 变更面
- 验证：`lake env lean CLPoly/Refinement/Hensel.lean` 通过；
  正式 `lake build CLPoly.Refinement.Hensel` 通过（3328 jobs）
- 下一步：把 adjustment 的 pairwise 证书沿 leaf origins 和 lift 输出传到
  最终 Hensel factor array

## Pairwise coprimality reaches the final Hensel array

`henselFactors_mod_pairwise_coprime` now follows the actual successful Hensel
entry trace.  It rewrites the full leaf range to the adjusted input array,
uses the pointwise `Forall₂` leaf-origin relation to transport each selected
pair through the lift, and finally applies the proved normalization-unit
coprimality theorem.  Thus distinct positions of the returned Hensel factor
array remain coprime modulo the selected prime whenever the concrete
adjustment invariant supplies pairwise coprimality.

The premise is exactly the universally quantified `adjustedPairwise` field
already required by the real Hensel entry invariant; it cannot prescribe an
output or encode a factorization witness.

## 度量（end-to-end Hensel pairwise trace）

- 耗时：约 0.5 小时（Forall₂ 位置来源、normalization 组合）
- C++ 变化：无，因此本次无新的 C++ b2b 变更面
- 验证：`lake env lean CLPoly/Refinement/FactorZZ.lean` 通过；
  正式 `lake build CLPoly.Refinement.FactorZZ` 通过（3545 jobs）
- 下一步：用最终 Hensel 模因子的 irreducibility + pairwise coprimality
  证明任意非平凡整数因子唯一对应一个合法 active subset

## Every modular divisor is a concrete atom sublist product

`divisor_associated_sublist_product` proves constructively, by recursion over
the actual ordered atom list, that every divisor of its product is associated
to the product of a genuine `List.Sublist`.  At each head irreducible, the
proof either cancels that atom from both sides when it divides the candidate,
or obtains coprimality from non-divisibility and cancels the atom from the
ambient product.  The resulting sublist retains source order and occurrence
identity, including lists that might contain equal values at different
positions.

This is the algebraic half of the integer-factor-to-C++-candidate bridge.  It
uses neither `WfDvdMonoid.exists_factors` to invent a replacement output nor
an assumed subset correspondence.  The remaining half is to encode this
sublist as the strictly increasing index array consumed by the generated
attempt.

## 度量（modular divisor sublist）

- 耗时：约 0.5 小时（不可约首项二分、domain cancellation、Sublist 递归）
- C++ 变化：无，因此本次无新的 C++ b2b 变更面
- 验证：`lake env lean CLPoly/Refinement/FactorZZ.lean` 通过；
  正式 `lake build CLPoly.Refinement.FactorZZ` 通过（3545 jobs）
- 下一步：把 `List.Sublist` 证明转成合法严格递增 source indices，并证明
  generated `trialProductLoop` 对该 candidate 计算同一模乘积

## Source sublists become actual generated scan candidates

`sublist_exists_legal_source_indices` follows the concrete `List.Sublist`
derivation and constructs occurrence indices in source order.  Skipping a
source head increments every retained index; selecting it prepends zero and
increments the tail.  This deliberately avoids `indexOf`, so equal atom
values at distinct positions retain their occurrence identity.

The construction proves the same accumulated strict-gap and upper-bound
invariant used by the generated C++ combination enumerator.  Its array
wrapper `sublist_exists_legal_combination` therefore supplies a genuine
`LegalCombination source.length chosen.length`, while selecting those source
positions recovers exactly `chosen`.  Together with scan no-omission, this
turns the algebraic sublist witness into an actually attempted candidate.

## 度量（sublist-to-generated-candidate）

- 耗时：约 1 小时（occurrence-preserving index recursion、累计 gap、数组桥接）
- C++ 变化：无，因此本次无新的 C++ b2b 变更面
- 验证：单文件 `lake env lean CLPoly/Refinement/Recombine.lean` 通过；
  正式 `lake build CLPoly.Refinement.Recombine` 通过（3531 jobs）
- 下一步：证明 `combinationToInt32` 与 generated `trialProductLoop` 对该
  candidate 成功执行并计算对应 Hensel 原子乘积

## Checked candidate lowering preserves every source index

`combinationToInt32Loop_toList` now executes the generated checked lowering
loop and identifies its exact output as the already-emitted prefix followed
by the converted source suffix.  The proof traverses the real recursion and
must discharge the concrete `index < 2^31` branch at every emitted cell.

The fixed-width bridge proves that a checked `Nat.toUInt32.toInt32` is
nonnegative and round-trips through the exact `Int32.toInt64.toNat` lookup
path.  Consequently `combinationToInt32_candidate_valid` derives the
`CandidateIndicesValid` premise consumed by `trialProductLoop_refines`; its
legal-combination wrapper connects this directly to generated scan
candidates.  Overflow is not idealized away: success depends on the same
source check, and the theorem proves its exact consequence.

## 度量（checked candidate lowering）

- 耗时：约 1 小时（递归输出、Int32 符号边界、逐项 round-trip）
- C++ 变化：无，因此本次无新的 C++ b2b 变更面
- 验证：单文件 `lake env lean CLPoly/Refinement/Recombine.lean` 通过；
  正式 `lake build CLPoly.Refinement.Recombine` 通过（3531 jobs）
- 下一步：把 converted candidate 的 `SelectedProductMod` 精确改写为前一
  步恢复的 sublist product，并接到 `zassenhausAttempt` 的执行分支

## Generated trial multiplication computes the recovered sublist product

`selectedProductMod_combinationToInt32` identifies every lookup performed by
the converted `Int32` candidate with the same occurrence in the original
active Hensel array.  It combines the exact conversion output list with the
proved fixed-width round trip, then rewrites the generated semantic product
to the occurrence-sensitive `selectSourceIndices` list product.

`trialProductLoop_source_indices_refines` composes that identity with the
existing execution proof for the generated multiplication-and-reduction
loop.  Thus a successful concrete run computes the initial product times the
exact modular product of the recovered atom sublist.  There is no abstract
candidate-product callback or post-hoc witness in this bridge.

## 度量（generated candidate product bridge）

- 耗时：约 0.75 小时（转换数组逐项对应、source lookup、乘积组合）
- C++ 变化：无，因此本次无新的 C++ b2b 变更面
- 验证：单文件 `lake env lean CLPoly/Refinement/Recombine.lean` 通过；
  正式 `lake build CLPoly.Refinement.Recombine` 通过（3531 jobs）
- 下一步：进入 `zassenhausAttempt`，证明对应整数真因子的 candidate 不会
  被 leading/constant pruning 拒绝，并由 symmetric recovery 恢复该因子

## Exact execution of both Zassenhaus pruning products

`selectedLeadingProductLoop_succeeds` executes the generated leading-product
recursion for a bounded candidate whose selected active factors are nonempty.
It proves that the returned integer is the incoming accumulator times the
ordered product of the concrete coefficient stored at position zero of every
selected factor.  The out-of-bounds and empty-factor errors are discharged
at their actual source branches.

`selectedConstantProductLoop_succeeds` gives the analogous exact execution
for the generated `constantTerm` pruning product.  These results expose the
two integers tested by `zassenhausAttempt`; the next algebraic step can now
prove that the recovered divisor's leading and constant coefficients make
both remainder tests zero rather than assuming that pruning accepts it.

## 度量（Zassenhaus pruning loops）

- 耗时：约 0.75 小时（两个良基循环、越界/空数组分支、累积乘积）
- C++ 变化：无，因此本次无新的 C++ b2b 变更面
- 验证：单文件 `lake env lean CLPoly/Refinement/Recombine.lean` 通过；
  正式 `lake build CLPoly.Refinement.Recombine` 通过（3531 jobs）
- 下一步：证明这两个实际整数乘积对应候选多项式的 leading/constant
  coefficient，并由整数真因子关系推出两个 pruning remainder 均为零

## Polynomial divisors cannot be falsely rejected by boundary pruning

`polynomial_divisor_boundary_coefficients` transports an actual integer
polynomial divisibility proof to divisibility of both the leading and
constant coefficients.  The leading result uses the polynomial
`leadingCoeff` monoid homomorphism; the constant result expands the real
quotient product at coefficient zero.

`zassenhaus_prune_condition_false_of_dvd` unfolds the generated arithmetic
model `ZZ.fdiv_r = Int.fmod` and proves directly that a divisor makes the
tested remainder zero.  Its leading and constant specializations allow the
concrete recovered integer to differ by a unit, matching the primitive and
normalization steps, while still proving that neither generated rejection
condition can hold.

## 度量（boundary-pruning soundness）

- 耗时：约 0.5 小时（多项式边界系数整除、associated 传递、fmod 零余数）
- C++ 变化：无，因此本次无新的 C++ b2b 变更面
- 验证：单文件 `lake env lean CLPoly/Refinement/Recombine.lean` 通过；
  正式 `lake build CLPoly.Refinement.Recombine` 通过（3531 jobs）
- 下一步：把两个 generated loop 的具体返回值与 recovered divisor 的
  leading/constant coefficient 建立 associated 关系，完成 pruning 分支组合

## Pruning products are candidate-polynomial boundary coefficients

`leadingCoeff_list_prod` and `coeff_zero_list_prod` record the two boundary
coefficients of an ordered list product over integer polynomials.
`selectedLeadingValues_prod_eq_leadingCoeff` then transports the actual
source-array head coefficients selected by the generated candidate to the
leading coefficient of the exact occurrence-sensitive sublist product.

`selectedConstantValues_prod_eq_coeff_zero` provides the constant analogue
using the generated `constantTerm` values.  Both proofs preserve the selected
source positions rather than replacing the factor list extensionally.  The
remaining per-factor head/constant premises will be supplied from the Hensel
canonical/origin invariant before these lemmas are used by the final entry.

## 度量（candidate boundary coefficients）

- 耗时：约 0.75 小时（list product 边界系数、source occurrence lookup）
- C++ 变化：无，因此本次无新的 C++ b2b 变更面
- 验证：单文件 `lake env lean CLPoly/Refinement/Recombine.lean` 通过；
  正式 `lake build CLPoly.Refinement.Recombine` 通过（3531 jobs）
- 下一步：由最终 Hensel array 的 canonical/origin 证明每个选中因子的
  head/constant 表示定理，并组合两个 pruning accept 分支

## Canonical sparse integer heads are mathematical leading coefficients

`sparsePolyZZ_leadingCoeff_eq_head` proves directly for the integer sparse
representation that a nonempty canonical array stores the mathematical
leading coefficient in its first cell.  A dedicated integer-chain lemma
shows every tail degree is strictly below the head, and a coefficient lemma
shows the complete tail vanishes at the head degree.  Head coefficient
nonzeroness then identifies the exact polynomial degree and leading
coefficient.

This removes the per-factor head premise from the leading-pruning bridge once
the already-traced Hensel output canonicality is supplied.  It does not reuse
the finite-field result or assume an ordered semantic polynomial.

## 度量（canonical ZZ leading coefficient）

- 耗时：约 0.75 小时（整数 strict chain、tail coefficient、degree/head）
- C++ 变化：无，因此本次无新的 C++ b2b 变更面
- 验证：单文件 `lake env lean CLPoly/Refinement/Recombine.lean` 通过；
  正式 `lake build CLPoly.Refinement.Recombine` 通过（3531 jobs）
- 下一步：沿 canonical array 的实际末项证明 generated `constantTerm`
  等于 `SparsePolyZZ.toPoly.coeff 0`，再组合完整 pruning accept 控制流

## Canonical sparse integer tails are the generated constant term

`sparseZZLastConstant` mirrors the source last-element lookup on lists.
`sparseZZ_coeff_zero_eq_lastConstant` proves by strict-chain recursion that
all earlier terms have positive degree, leaving exactly the last stored term
when its degree is zero.  `sparsePolyZZ_constantTerm_eq_coeff_zero` then
identifies this list result with the generated array access at `size - 1`.

The canonical wrappers for both selected boundary products now derive their
per-factor head and constant facts internally.  Their remaining premises are
only concrete array bounds, canonicality, and (for leading access) nonempty
selected factors, all already present in the Hensel execution invariants.

## 度量（canonical ZZ constant coefficient）

- 耗时：约 1 小时（末项递归、strict-chain coeff 0、array/list last lookup）
- C++ 变化：无，因此本次无新的 C++ b2b 变更面
- 验证：单文件 `lake env lean CLPoly/Refinement/Recombine.lean` 通过；
  正式 `lake build CLPoly.Refinement.Recombine` 通过（3531 jobs）
- 下一步：从 Hensel output invariant 实例化两个 canonical boundary
  wrappers，并在 `zassenhausAttempt` 中关闭 leading/constant accept 分支

## Natural-language proof: canonical final Hensel normalization

The generated normalization block either leaves the extracted array unchanged
or replaces exactly index zero by
`modCoeffOutput (scaleCoeffs first inverse) m`.  For the changed factor, unfold
the two concrete generated traversals as a single `filterMap` over the original
term list.  Every retained term has the same monomial degree as its source, so
the strict descending chain follows from the source chain by `Pairwise.filterMap`;
the generated guard supplies the retained coefficient's nonzeroness.  This
direct argument is important because `scaleCoeffs first inverse` itself need
not be canonical when a scaled coefficient is zero.

For an arbitrary output index, split on whether it is zero.  At zero, use the
actual `HenselNormalizeCorrect.normalized` run equation and the direct
filter-map theorem above.  At every other index, `Array.set!` reads the exact
unchanged input cell.  The empty and already-one source branches are identities.
Thus pointwise canonicality is transported through the real generated branch,
without adding an output field to an execution invariant.

Nonemptiness will be derived separately from the already-proved modular origin
and unit-scaling relation: an adjusted irreducible finite-field factor is
nonzero; multiplication by the normalization unit remains nonzero; an empty
sparse array decodes to the zero polynomial, contradiction.  This avoids an
extra representation oracle and covers both normalized and untouched factors.

`modCoeffOutput_scaleCoeffs_canonical` now implements the direct traversal
argument.  `HenselNormalizeCorrect.canonical` lifts it through all three exact
normalization trace constructors and proves every output index canonical.

## 度量（Hensel final-normalization canonicality）

- 耗时：约 0.75 小时（natural-language proof、filterMap chain、dependent array index）
- 迭代：3 轮单文件编译-修复
- Lean 新增/修改行数：约 95 行（含证明草稿）
- 对应 C++ 行数：约 12 行最终 normalization 分支；C++ 变化为无，因此无新的 C++ b2b 变更面
- 放弃的方案：先证明中间 `scaleCoeffs` canonical；当缩放产生零系数时该命题为假，改为直接证明后续过滤结果
- 验证：单文件 `lake env lean CLPoly/Refinement/Hensel.lean` 通过
- 下一步：沿实际 builder/step/lift/extract trace 传递叶因子的 canonicality，
  再由 modular irreducibility 推出逐项 nonempty

## Natural-language proof: canonical generated pair-vector addition

Track the actual two source cursors and the generated output prefix.  The
cursor invariant states that the prefix is canonical and every term remaining
in either input suffix has degree strictly below every emitted prefix term.
At a source step, canonicality of an input says its cursor tail is strictly
below its current term.  The generated comparison selects the larger current
degree; therefore the other complete suffix is also below the selected term.
For equal degrees both cursors advance and `pushNonzero` either appends their
coefficient sum or emits nothing, exactly matching the C++ zero-compaction.

A generic cursor-push lemma preserves all three invariant clauses.  Functional
induction on `pairVecAddLoop` then follows its five advancing branches and its
terminal branch.  At the public `0,0,#[]` entry, the prefix clauses are
vacuous, yielding canonicality of the exact generated sparse addition result.
No semantic polynomial equality is used to infer a representation property.

`PairVecAddCursorValid`, its two concrete push lemmas, and
`pairVecAddLoop_canonical` now implement this proof.  `pairVecAdd_canonical`
specializes it to the exact empty-prefix entry invoked by generated Hensel
arithmetic.

## 度量（generated sparse addition canonicality）

- 耗时：约 1 小时（frontier invariant、五个 source branches、数组 cursor API）
- 迭代：2 轮单文件编译-修复
- Lean 新增/修改行数：约 215 行（含证明草稿）
- 对应 C++ 行数：约 25 行 `pair_vec_add` iterator merge；C++ 变化为无，因此无新的 C++ b2b 变更面
- 放弃的方案：从多项式语义等式反推 canonical；表示不唯一，无法成立
- 验证：单文件 `lake env lean CLPoly/Refinement/Hensel.lean` 通过
- 下一步：证明 generated scale/mul/sub 的 canonical preservation，组合出
  每个真实 Hensel step 的 `g/h` canonicality

## Natural-language proof: canonical generated coefficient scaling

The generated `scaleCoeffs` is the source range-for loop represented as an
array map.  Mapping retains every monomial and therefore preserves the exact
strict degree chain.  If the input coefficients and the scalar are nonzero,
integer multiplication has no zero divisors, so every mapped coefficient is
also nonzero.  This proves canonicality directly for the generated array; the
Hensel caller supplies scalar nonzeroness from its positive modulus.

## Natural-language proof: canonical generated pair-vector subtraction

Subtraction uses the same two-cursor ordering decisions and emitted-prefix
frontier as addition.  A term selected only from the right input is emitted
with its coefficient negated; integer negation preserves nonzeroness.  Equal
degrees advance both cursors and call the same `pushNonzero` compactor on the
coefficient difference.  Therefore the addition cursor invariant and generic
push lemmas apply branch-for-branch to the actual `pairVecSubLoop`, with only
the right-only coefficient proof changed from `c ≠ 0` to `-c ≠ 0`.

## Natural-language proof: canonical generated multiplication heap

The generated multiplication heap loop does not rely on product-list order.
At each nonempty frontier it computes the actual maximum pending degree,
sums every contribution at that degree, filters all of them out, and calls
`pushNonzero`.  Track that every pending degree is below every already emitted
term.  The list maximum theorem makes the chosen degree no larger than no
pending degree and at least every pending degree.  After filtering equality,
every remaining degree is therefore strictly smaller than the chosen degree.
`pushNonzero` either appends the summed coefficient or drops a zero sum, so it
preserves both canonicality and the strict frontier.  Well-founded functional
induction is exactly on the generated filtered-list length; the empty-prefix
entry makes the initial frontier condition vacuous.  This proof applies to
all generated product lists and hence to the concrete `pairVecMulProducts`.

`scaleCoeffs_canonical`, `pairVecSubLoop_canonical`, and
`pairVecMulHeapLoop_canonical` now close all three generated arithmetic
representation primitives needed by a Hensel step.  The multiplication result
is canonical even without canonical input arrays because its concrete heap
always groups a maximal degree before emitting it.

## 度量（generated scale/sub/mul canonicality）

- 耗时：约 2 小时（scale map、subtraction five-way merge、heap maximum/filter induction）
- 迭代：1 轮 scale/sub 与 5 轮 multiplication 单文件编译-修复
- Lean 新增/修改行数：约 310 行（含三段证明草稿）
- 对应 C++ 行数：约 55 行 scale/sub/multiplication heap lowering；C++ 变化为无，因此无新的 C++ b2b 变更面
- 放弃的方案：要求 multiplication product list 预排序；实际 heap 输入是 flatMap，证明改为使用每轮真实 max/filter 语义
- 验证：单文件 `lake env lean CLPoly/Refinement/Hensel.lean` 通过
- 下一步：由实际 divmod trace 取得 remainder canonicality，组合完整
  `__hensel_step_raw_ir` 的 `g/h` canonical preservation

## Natural-language proof: canonical concrete divmod remainder

Induct on the exact finite `DivmodTrace` consumed by the generated while-loop.
The done constructor returns its current remainder.  The vanished constructor
removes the current leading cell; its active guard proves the array is
nonempty, so erasing index zero preserves the strict chain and nonzero tail.
The subtract constructor replaces the remainder with the actual generated
`divmodRemainder`; its two-cursor merge canonicality theorem applies to the
current remainder and fixed divisor.  In both recursive constructors, apply
the induction hypothesis to the exact next trace.  The quotient accumulator
is irrelevant to this remainder property and is never replaced by an oracle.

`Generated.StrictHensel.DivmodTrace.remainder_canonical` now implements this
constructor induction and applies directly to the exact trace stored in a
successful generated divmod call.

## 度量（divmod trace remainder canonicality）

- 耗时：约 0.5 小时（trace induction、active/nonempty bridge）
- 迭代：1 轮单文件编译
- Lean 新增/修改行数：约 45 行（含证明草稿）
- 对应 C++ 行数：约 12 行 modular long-division while branches；C++ 变化为无，因此无新的 C++ b2b 变更面

## Successful raw Hensel operations preserve canonical sparse form

The next Hensel composition step needs representation facts about the exact
values returned by generated raw execution, rather than facts about an
independent specification-level operation.

- `divideThenReduceCoeffs` is the generated C++ coefficient traversal.  Its
  `filterMap` leaves monomial degrees unchanged, so strict descending order is
  inherited from the source array, and its guard removes every coefficient
  equal to zero.
- A successful generated multiplication, addition, or subtraction equation
  fixes the returned value to the output of the already verified concrete
  heap/merge loop.  Canonicality therefore transfers by injection of the
  actual `RawExec.ok` equation; no semantic arithmetic oracle is introduced.
- A successful generated modular-divmod entry excludes, in the same source
  branch order, an empty divisor, an invalid termination input, and a failed
  modular inverse.  These facts reconstruct the exact finite `DivmodTrace`.
  The previously proved trace invariant then gives canonicality of the actual
  remainder, and equality of successful `RawExec` results identifies it with
  the caller-visible remainder.

These bridges are intended for direct composition through the generated
Hensel factor phase.  They neither add fuel nor replace the C++ execution with
an existential quotient/remainder result.

The composition is now explicit in
`henselFactorCorrection_canonical_from_raw_runs`: it follows the generated
sequence `gh`, `difference`, `e`, `se`, concrete divmod, `tau`, and the final
`gNew`/`hNew` reductions.  The complete factor-phase theorem exports these two
canonicality facts.  `HenselStepRefinementInvariant` consequently requires
canonical target and input factors, while `HenselStepCorrect` guarantees that
the full two-phase C++ step returns canonical `g` and `h`; the second phase
changes only `s` and `t`, as recorded by its generated execution proof.

- C++ changes: none, so there is no new C++ b2b change surface in this step.

## The recursive tree builder returns a canonical node array

The local node-write theorem is now threaded through the complete generated
well-founded tree builder.

- Algebra-only pointer updates preserve `HenselNodeCanonical` because they do
  not change `g` or `h`.
- A same-size `HenselTreeArrayFrame` transfers the global array invariant by
  using the new canonical value at the selected source slot and exact equality
  at every other existing slot.
- Each source `push default` preserves the invariant: the new allocation has
  empty canonical `g`/`h` fields before the recursive call fills it.
- The left recursion consumes the concrete left-ready array; the optional
  right allocation and recursion consume the actual array returned by the
  left call.  All four source branch combinations return
  `HenselArrayCanonical output` together with their existing raw, topology,
  gcd, and semantic certificates.

The strongest tree-build theorem, its generated public contract, and the
`Pipeline/L1` boundary now expose this array certificate.  It is therefore a
direct input to the previously proved lift-loop and extraction canonicality
chain.

- C++ changes: none, so there is no new C++ b2b change surface in this step.

## Canonical fields written by the concrete tree builder

The tree builder converts finite-field range products to integer sparse arrays
before storing them in each Hensel node.  The generated
`henselTreeZpToZZIR` is now proved to preserve strict descending degree order
and nonzero coefficients directly through its array `map`.  Nonzero transfer
uses the injectivity of `UInt64.toNat` followed by the exact `Nat → Int`
embedding; it is not inferred from polynomial denotation.

`henselTreeStoreNodeRawIR_canonical` then follows the same six checked writes
as the C++ node-store helper.  Its successful raw equation identifies the
stored `g` and `h` with the two converted product-loop outputs, whose source
canonicality is available from the concrete finite-field multiplication
trace.  This is the local write fact needed by the recursive builder array
frame.

- C++ changes: none, so there is no new C++ b2b change surface in this step.

## Canonical factor fields survive recursive lifting and extraction

Canonicality is now carried above a single Hensel step without assuming that
recursive children leave unrelated slots untouched.

- `HenselArrayCanonical` quantifies over every allocated node slot.  Its
  checked-`set!` frame theorem replaces the selected slot with the canonical
  `g`/`h` returned by the concrete step and reuses the prior invariant at all
  other indices.
- `HenselLiftRecursiveCorrect.arrayCanonical` follows the exact well-founded
  tree trace.  The left and right recursive certificates consume the array
  produced by the preceding source call, so the proof matches the generated
  mutation order rather than reasoning about a mathematical replacement tree.
- `HenselLiftLoopCorrect.arrayCanonical` iterates that result along the actual
  quadratic-precision loop trace.
- `HenselExtractCorrect.outputCanonical` follows all four generated extraction
  branch shapes and proves that each concrete `push` appends only a canonical
  `g` or `h` read from the lifted node array.

The remaining input to this chain is a builder theorem establishing
`HenselArrayCanonical` for the array returned by the concrete tree builder.

- C++ changes: none, so there is no new C++ b2b change surface in this step.
- 放弃的方案：同时证明 quotient canonical；当前 Hensel `g/h` 链只消费 remainder，且 multiplication heap 对任意输入已证明 canonical
- 验证：单文件 `lake env lean CLPoly/Refinement/Hensel.lean` 通过
- 下一步：从 successful concrete divmod run 提取该 remainder 定理并组合 Hensel factor phase

## The complete generated Hensel entry returns canonical factors

The canonical certificate produced by the concrete tree builder is now
consumed by the full generated `__hensel_lift_upoly` refinement theorem.

- The actual quadratic lift-loop trace transfers `HenselArrayCanonical` from
  the builder output to the concrete lifted node array.
- The actual extraction traversal starts from the source's empty result array
  and appends only canonical `g`/`h` fields read from those lifted nodes.
- `HenselNormalizeCorrect.outputCanonical` lifts the existing pointwise
  normalization theorem to the array contract used by recombination.  It
  covers the unchanged branches and the exact source rewrite of slot zero.
- `HenselLiftEntryCorrect` now includes
  `HenselFactorArrayCanonical output.1`; the property therefore constrains the
  factors returned by the successful generated entry, rather than remaining
  an unconnected collection of helper lemmas.

No fuel, partial definition, semantic oracle, or L2 fallback is introduced.
Every certificate is transported along the intermediate arrays returned by
the corresponding strict raw stage.

- C++ changes: none, so there is no new C++ b2b change surface in this step.
- Verification: `lake env lean CLPoly/Refinement/Hensel.lean` passes with no
  error or `sorry` warning.

## The concrete primitive-content loop is an exact gcd fold

The next irreducibility obligation starts inside the generated
`primitiveRaw` used by every successful Zassenhaus candidate.  Its recursive
`contentLoop` is now identified with a source-order `List.foldl Nat.gcd` over
the remaining concrete coefficient cells.

The proof unfolds the generated well-founded loop at each array index, aligns
`Array.getElem` with the corresponding `toList.drop` head, and recurses on
`input.size - index`.  This supplies the executable gcd characterization
needed to prove that coefficient-wise exact division returns an actually
primitive polynomial; it does not assume a mathematical primitive-part
oracle.

- C++ changes: none, so there is no new C++ b2b change surface in this step.
- Verification: both `lake env lean CLPoly/Refinement/Recombine.lean` and the
  formal `lake build CLPoly.Refinement.Recombine CLPoly.Pipeline.L1` pass with
  no error or `sorry` warning.  The axiom audit reports only `propext`,
  `Classical.choice`, and `Quot.sound` (no project axiom or `sorryAx`).

## The generated primitive routine returns a primitive polynomial

The complete successful `primitiveRaw` execution is now proved to return an
`IsPrimitive` polynomial on every nonempty canonical sparse input.

- `contentLoop_dvd_iff` proves that the loop result has the exact universal
  divisibility property of the gcd of all remaining concrete coefficient
  cells.
- Canonical sparse degree order proves each stored coefficient is the
  corresponding coefficient of `SparsePolyZZ.toPoly`; conversely, divisibility
  of all stored coefficients transfers through the actual monomial sum to
  every mathematical coefficient.
- The two directions identify the nonnegative loop result with
  `Polynomial.content` by normalized gcd antisymmetry.
- The successful raw trace chooses either that gcd or its negation according
  to the actual leading-sign branch, then executes the checked coefficient
  division loop.  Taking content of its proved product equation and cancelling
  the nonzero input content yields output content one.

The proof never executes or assumes `Polynomial.primPart`; the output remains
the array returned by the generated C++ control flow.

- C++ changes: none, so there is no new C++ b2b change surface in this step.
- Verification: both `lake env lean CLPoly/Refinement/Recombine.lean` and the
  formal `lake build CLPoly.Refinement.Recombine CLPoly.Pipeline.L1` pass with
  no error or `sorry` warning.  The axiom audit reports only `propext`,
  `Classical.choice`, and `Quot.sound` (no project axiom or `sorryAx`).

## Concrete Zassenhaus candidates are canonical and primitive

The generated candidate path now carries a complete representation and
primitive-content certificate through its actual sparse operations.

- The integer normalization grouping fold is proved to retain one cell per
  degree.  A source-shaped proof of the strict `deg >` merge sort then uses
  that uniqueness to rule out the comparator's non-total equality case;
  `normalization_canonical` therefore certifies the actual C++ sort rather
  than replacing it with a mathematical sorting oracle.
- Exact list traces for `modCoeffLoop` and `symmetricModLoop` show that they
  preserve source order and omit precisely the concrete zero residues.
  Canonicality is transported through `multiplyNormalizeRaw`,
  `multiplyNormalizeModRaw`, `trialProductLoop`, and `symmetricModRaw`.
- A successful `zassenhausAttempt` returns the factor produced by the real
  symmetric-content division and the quotient produced by the real checked
  long division.  Both returned polynomials are now proved canonical and
  `IsPrimitive`; nonemptiness is derived from the successful division trace,
  not assumed.
- The same pair of certificates is propagated through the generated
  fixed-size combination scan by induction on its concrete termination rank.

No C++ files changed, so this step has no new C++ b2b change surface.

Verification: `lake env lean CLPoly/Refinement/Recombine.lean` and
`lake build CLPoly.Refinement.Recombine CLPoly.Pipeline.L1` pass with no error
or `sorry` warning.  The axiom audit reports only `propext`,
`Classical.choice`, and `Quot.sound`, with no project axiom or `sorryAx`.  The
next obligation is irreducibility: use exhaustion of all smaller source
subsets to contradict any proper factor of an extracted primitive candidate.

## Observe the concrete Hensel-modulus candidate at the selected prime

The successful Zassenhaus trace now remains exact when the Hensel modulus
`M` is reduced further at every positive divisor `p`.  This is the bridge
needed to use the selected-prime irreducibility certificates after lifting.

- The generated `% M` coefficient loop, normalized modular multiplication,
  and complete `trialProductLoop` are proved to preserve their polynomial
  value in `ZMod p` whenever `p ∣ M`.
- The generated `symmetricModLoop` and `symmetricModRaw` receive the analogous
  divisor theorem: the concrete symmetric representative modulo `M` reduces
  to the original polynomial modulo `p`.
- `zassenhausAttempt_extracted_candidate_trace` exposes the actual successful
  `combinationToInt32`, trial product, symmetric representative, and
  `primitiveRaw` result.  No intermediate value is chosen by a specification.
- Combining those exact executions proves that the returned factor, scaled by
  the content computed by the concrete primitive loop, equals modulo `p` the
  leading coefficient times the exact occurrence-sensitive selected sublist
  product.

No C++ files changed, so this step has no new C++ b2b change surface.
Verification: `lake env lean CLPoly/Refinement/Recombine.lean` and
`lake build CLPoly.Refinement.Recombine CLPoly.Pipeline.L1` pass with no error
or `sorry` warning.  The permanent axiom audit now includes the new concrete
candidate theorem and reports only `propext`, `Classical.choice`, and
`Quot.sound`.  The next step is to prove both concrete scalars are nonzero
modulo the selected prime and convert this equality to `Associated`, then
combine it with the earlier exhaustive-scan rejection theorem.

## Successful scans retain exact selected-prime certificates

The concrete successful-candidate equation is now upgraded to a modular
association certificate under the live selected-prime invariants.

- Every selected lifted factor is irreducible and hence nonzero modulo `p`,
  so the exact occurrence-sensitive candidate subproduct is nonzero.
- The current polynomial's surviving leading coefficient makes the concrete
  right-hand scalar nonzero.  The successful execution equation therefore
  forces the content scalar computed by `primitiveRaw` to be nonzero as well;
  both scalars are units in `ZMod p`, yielding `Associated` for the returned
  factor and exact selected subproduct.
- The actual exact-division and primitive-quotient trace proves that a
  successful extraction preserves nonzeroness of the quotient's leading
  coefficient modulo `p`.
- `scanZassenhausCombinations_extracted_mod_certificate` follows every
  generated rejected/next-candidate branch by the concrete combination rank
  and returns both facts for the actual successful candidate.

No C++ files changed, so this step has no new C++ b2b change surface.  The
single-file formal check passes without error or `sorry`.  The next obligation
is to carry these invariants through successful active-array removal and the
outer subset-size loop, retaining the rejection history for all smaller
candidate sizes needed to rule out proper factors.

## Active-factor erasure preserves concrete origins

The successful outer-loop transition removes the selected active positions by
executing `removeCombinationLoop` from the candidate tail toward its head.
That exact reverse-erasure recursion is now proved origin preserving: every
factor in the returned active array is a member of the input array.

The resulting `removeCombination_preserves_pointwise` theorem transports any
pointwise invariant through the actual successful removal, and will be
instantiated with irreducibility of each lifted factor modulo the selected
prime.  It neither reconstructs an abstract complement nor assumes that the
candidate has no duplicates; successful bounds checks and the concrete array
erasures provide the whole proof.

No C++ files changed, so this step has no new C++ b2b change surface.  The
single-file check passes without error or `sorry`.  The next proof combines
this origin preservation with the scan certificate and the accumulated
smaller-subset rejection history in the outer Zassenhaus loop.

## Exact symmetric reconstruction uses the correct closed boundary

The generated symmetric-remainder execution now carries the coefficient
bound required by the later precision argument.

- For positive `M`, `ZZ.symmetricMod` returns a representative satisfying
  `2 * |coefficient| ≤ M`.  The proof follows both source branches after
  identifying positive-modulus `fmod` with `emod`; it deliberately retains
  equality at the positive midpoint for even `M`.
- The exact generated loop/list trace transports this bound to every emitted
  sparse cell.  Canonicality then promotes it to every mathematical
  coefficient of `SparsePolyZZ.toPoly`, including absent zero coefficients.
- A new recovery theorem proves uniqueness when the concrete representative
  has the closed bound and the true scaled factor has the strict bound.  This
  is the boundary actually supplied by `M > 2 * |lc| * B`; it avoids the
  unsound shortcut of demanding a strict symmetric bound at even moduli.

No C++ files changed, so this step has no new C++ b2b change surface.  The
single-file check passes without error or `sorry`.  Next, the generated
binomial, L2-square norm, and integer-square-root target helpers must be tied
to the strict scaled-factor coefficient bound before recovery is applied to
the actual candidate trace.

## Generated Mignotte arithmetic executes exact square and choose folds

Two executable cores of the concrete Mignotte target are now identified.

- The generated `Array.foldl` for `__upoly_norm_l2_sq_upoly_raw_ir` equals
  the source-order integer sum of every physically stored coefficient square;
  its result is nonnegative.
- The well-founded multiplicative `__binomial_loop_raw_ir` is proved equal to
  `Nat.choose` from every reachable `(i, choose n i)` state.  At each source
  iteration, `Nat.choose_succ_right_eq` proves the exact divisibility of
  `choose(n,i) * (n-i)` by `i+1`, so the generated integer division is not
  treated as an arithmetic oracle.

No C++ files changed, so this step has no new C++ b2b change surface.  The
single-file check passes without error or `sorry`.  The next boundary proves
the actual `Int64.ofInt(degree)` and `n/2` call reaches this Nat recurrence,
then establishes the generated Newton square-root upper bound.

## Actual Mignotte machine-degree call reaches the choose recurrence

The exact machine boundary around the generated binomial helper is now
closed.  For every sparse degree below the signed 64-bit limit,
`Int64.ofInt degree` is proved to preserve the degree, the actual signed
`n / 2` instruction is proved to equal the natural floor half, and both
`toNatClampNeg` conversions recover those same natural values.  The proof
also follows the generated negative/out-of-range guard, the zero/end-point
fast path, and the concrete `min(k, n-k)` branch before invoking the already
verified well-founded multiplicative loop.  Consequently the literal call
made by `__mignotte_bound_upoly_raw_ir` returns `choose degree (degree/2)`;
there is no wraparound, truncation, clamp, or mathematical replacement of
the C++ helper hidden at this boundary.

No C++ files changed, so this step has no new C++ b2b change surface.  The
single-file check passes without error or `sorry`.  Next, the generated
Newton loop and its final correction branch must be proved to return the
least integer upper square root of the concrete stored coefficient-square
sum.

## Generated Newton square root supplies the concrete norm upper bound

The remaining executable Mignotte helpers are now composed through their
actual generated control flow.

- The well-founded `__isqrt_ceil_loop_raw_ir` is proved extensionally equal
  to `Nat.sqrt.iter` by comparing the two recurrences at every current
  estimate.  The generated loop remains the executed term; the library
  iterator supplies its already verified floor-square-root bounds.
- The source bit-length estimate `2^((bits+1)/2)` is proved large enough from
  `ZZ.lt_pow_sizeinbase_nat`.  This discharges the exact initial condition
  needed by the Newton recurrence.
- The complete generated helper follows its nonpositive branch and its final
  `root*root < n ? root+1 : root` correction, proving a nonnegative result
  whose square is at least the original integer input.
- `__mignotte_bound_upoly_raw_ir` is now exposed as the product of the exact
  central `Nat.choose` result and that exact generated Newton norm, under the
  real signed-degree bound at its stored leading term.

No C++ files changed, so this step has no new C++ b2b change surface.
`Recombine` checks without errors or `sorry`.  Next, the classical Mignotte
coefficient inequality must be instantiated for a hypothetical primitive
factor and linked to this concrete stored norm before the actual exhaustive
candidate scan can recover it.

## Actual generated Mignotte value bounds every true divisor coefficient

The mathematical coefficient bound is now connected to the exact sparse
execution value.

- A structural induction over the canonical strictly descending sparse list
  proves that the source-order sum of stored coefficient squares equals the
  sum of all mathematical coefficient squares over any containing degree
  range.  The proof handles the head degree explicitly and proves the tail
  coefficient there is zero, so duplicate-degree cancellation cannot be
  hidden in the bridge.
- The canonical head cell is proved to carry the exact mathematical
  `natDegree`.  This aligns the central `choose` parameter used by C++ with
  the degree used by `mignotte_bound_l2`.
- The real Mahler/L2 Mignotte inequality is instantiated for an arbitrary
  genuine integer divisor.  The generated Newton theorem turns its real
  square root into the concrete integer norm, and the result is cast back to
  an exact integer coefficient inequality.
- The final theorem includes the raw execution equation for
  `__mignotte_bound_upoly_raw_ir`: the very integer returned by that call is
  nonnegative and bounds every coefficient of the supplied true divisor.
  It does not choose an abstract bound or execute an L2 factorization.

No C++ files changed, so this step has no new C++ b2b change surface.  The
next step multiplies this bound by the actual leading-coefficient scale used
by `__hensel_lift`, derives the strict modulus inequality from its concrete
target loop, and feeds it into symmetric recovery for the enumerated lifted
subproduct.

## Concrete Hensel precision reaches the unique recovery interval

The generated Mignotte coefficient bound is now carried through the actual
default-precision Hensel execution.

- `hensel_output_modulus_bounds_scaled_divisor` consumes a genuine
  `HenselLiftEntryCorrect` trace.  It identifies the zero-exponent target with
  the exact generated Mignotte return value, then uses the well-founded
  lifting loop's own `outputM_gt_target` certificate.  Thus every coefficient
  of `leadingCoeff * g`, for any genuine divisor `g`, lies strictly inside
  half of the returned modulus.
- The proof rules out the positive-exponent constructor at the literal
  `Int32` zero boundary and identifies the source leading cell through its
  actual successful array lookup.  It does not accept a caller-supplied
  modulus or coefficient bound.
- `symmetricModRaw_recovers_strictly_bounded_target` combines the generated
  symmetric-remainder execution, its closed half-interval output bound, and
  congruence modulo that same modulus.  The asymmetric strict bound on the
  true target supplies uniqueness even when the modulus is even.

No C++ files changed, so this step has no new C++ b2b change surface.  The
next step instantiates the congruence premise with the actual selected Hensel
subproduct and then proves that the generated exact-division attempt accepts
the recovered true factor.

## Every true divisor has a legal candidate in the actual Hensel array

The cross-stage SelectPrime-to-Hensel completeness bridge now retains exact
array order and occurrences.

- `selectionHenselFactors_pointwise_associated` follows the generated first
  factor adjustment, lifted leaf-origin `Forall₂` trace, and final generated
  normalization.  At each concrete array index, the selected-prime factor and
  returned Hensel factor differ only by the two unit scalars actually produced
  by those stages.
- The pointwise trace is lifted to the product of the concrete arrays.  The
  selected-prime factorization certificate therefore makes the reduction of
  any genuine integer divisor divide the product of the actual normalized
  Hensel output.
- Recursive divisibility through that concrete list of irreducibles produces
  an occurrence-sensitive `List.Sublist`; `sublist_map_iff` brings the witness
  back from modular polynomials to the actual sparse output objects.
- `integer_divisor_mod_has_legal_hensel_candidate` converts that sublist to the
  exact strictly increasing natural-index array accepted by the generated
  Zassenhaus scanner.  Selecting those indices from the source array recovers
  the same sublist exactly, so no reordered set or semantic factor oracle is
  substituted at the enumeration boundary.

No C++ files changed, so this step has no new C++ b2b change surface.  Next,
the selected product's congruence at the full Hensel modulus is combined with
the concrete recovery bound, after which the generated pruning and exact
division branches can be shown to accept this legal candidate.

## Exact mutable-array frame for recursive Hensel traversal

The first structural prerequisite for the full-precision leaf-product proof
is now discharged directly from the generated recursive execution trace.

- `henselLiftTreeIndices` records, in preorder, the exact finite node indices
  visited by a generated `HenselLiftTree`; it carries no polynomial values or
  semantic outputs.
- `HenselLiftRecursiveCorrect.preserves_not_mem` inducts over the constructors
  anchored to successful `__hensel_lift_recursive_raw_ir` calls.  The only
  local mutation is the actual `set!` at the current root, and the two
  recursive hypotheses compose in source left-then-right order.
- Consequently every lookup outside the concrete traversal topology is
  exactly unchanged, not merely congruent modulo a number.  This will prove
  that child traversal cannot alter its parent and that the right traversal
  cannot invalidate the already established left-subtree product.

No C++ files changed, so this step has no new C++ b2b change surface.  Next,
canonical preorder topology disjointness is combined with this frame theorem
to prove that the actual extracted leaf product equals the round target
modulo `m^2`, then inductively modulo the returned `outputM`.

## Canonical Hensel topology indices are disjoint contiguous blocks

The generated tree builder now has the exact index-separation theorem needed
to compose mutable left and right traversals.

- `henselTreeBuildTopology_indices_nodup_bounded` follows the same midpoint
  split and four optional-child branches as `henselTreeBuildTopology`; its
  recursion decreases only by the source interval length `stop - start`.
- The preorder list `henselLiftTreeIndices` is proved duplicate-free.  When
  both children exist, the proof uses the concrete left interval ending at
  `root + 1 + leftCount` and the right interval beginning at that exact value,
  rather than assuming abstract tree disjointness.
- Every visited index lies in the half-open block
  `[root, root + henselTreeInternalNodeCount start stop)`.  Empty-child cases
  use the generated internal-node count equations, so singleton topologies
  are covered without inventing placeholder nodes.

No C++ files changed, so this step has no new C++ b2b change surface.  The
next proof combines these exact blocks with
`HenselLiftRecursiveCorrect.preserves_not_mem` to retain parent and completed
left-subtree values across the actual source-order recursive mutations.

## Actual recursive Hensel leaves multiply to the round target

The algebraic product invariant is now carried through the exact generated
left-then-right mutable traversal, rather than through an abstract Hensel
existence result.

- `henselExtractedFactors_eq_of_lookups` proves that the pure source-order
  extraction reads identically from two arrays whenever every concrete tree
  index has the same `getElem?` result.  Its recursion decreases by the actual
  generated tree `nodeCount`.
- `HenselLiftRecursiveCorrect.extractedFactors_product` inducts over the four
  constructors anchored to successful generated recursive calls.  At each
  node it uses the concrete `set!` result and the real `HenselStepCorrect`
  product at `m^2`.
- In the two-child case, topology `Nodup` splits into root exclusion and
  left/right disjointness.  The right traversal's exact frame theorem then
  transports the completed left extraction into the final array before the
  two induction products are combined.

No C++ files changed, so this step has no new C++ b2b change surface.  Next,
the canonical topology theorem instantiates this round result inside the
actual well-founded Hensel lift loop, yielding the leaf product at its returned
full recovery modulus.

## Actual Hensel loop preserves the leaf product to its returned modulus

`HenselLiftLoopCorrect.extractedFactors_product` now carries the complete
leaf-product equation through the generated quadratic-precision loop.

- The done constructor returns the supplied equation at the literal current
  modulus and unchanged node array.
- The step constructor obtains the next equation solely from the actual
  recursive traversal's `extractedFactors_product` theorem at `m^2`, rewrites
  the source loop's literal `m * m`, and passes that equation into the real
  recursive tail.
- The proof inducts over `HenselLiftLoopCorrect`, whose construction is tied
  to successful `__hensel_lift_loop_raw_ir` calls and whose termination uses
  the strictly increasing modulus measure.  It introduces no fuel and no
  alternate semantic Hensel execution.
- `HenselTreeSemanticBuildCertificate.extractedFactors_product` folds the
  builder's exact occurrence-preserving leaf-origin trace into the product of
  the actual adjusted finite-field input interval.
- `liftLoop_extractedFactors_product` then instantiates canonical topology
  `Nodup`, combines that builder product with an explicit upstream source
  factorization equation, and returns the product at the loop's concrete
  `outputM`.  The upstream equation remains visible rather than becoming a
  hidden Hensel oracle.

No C++ files changed, so this step has no new C++ b2b change surface.  Next,
extraction and normalization transport this concrete output-node equation to
the public Hensel entry result.

## Selected-prime factors seed the exact Hensel product

The explicit source-product premise of the Hensel loop has now been
discharged from the actual selected-prime result and the generated first
factor adjustment.

- `HenselAdjustFirstFactorCorrect.product_eq` unfolds the successful source
  trace, reads the concrete first factor, executes the generated coefficient
  scaling and normalization semantics, and proves the product effect of the
  literal `set! 0` write.  The scalar is exactly `lc(f) mod p`.
- `polynomial_eq_C_leadingCoeff_mul_of_associated_monic` recovers an exact
  finite-field equation from the selected candidate's associatedness.  Its
  unit is determined by the actual mapped leading coefficient and by the
  fact that every selected factor is monic.
- `selectionAdjustedFactors_product_eq_source` combines those facts, so the
  complete adjusted array multiplies exactly to the integer source reduced
  modulo the selected prime.  It does not choose a replacement factor list or
  assume a semantic factorization oracle.
- `selectionHenselFactors_preNormalization_product` supplies that equation to
  `HenselLiftEntryCorrect.preNormalizationProduct`.  Hence the factors read
  from the actual generated tree multiply to the source at the exact modulus
  returned by the well-founded quadratic lift loop.
- `HenselNormalizeCorrect.product_eq_unit_mul` records the exact unit scaling
  performed by the subsequent generated normalization block, ready for the
  recombination boundary.
- `selectionHenselFactors_normalized_product_eq_unit_mul_source` composes that
  literal normalization trace with the pre-normalized product.  The public
  Hensel output array therefore multiplies to the source modulo its actual
  returned modulus, up to precisely the unit computed by C++ normalization.

No C++ files changed, so this step has no new C++ b2b change surface.  Next,
the normalization unit and the actual exact-divmod combination scan are used
to complete the occurrence-sensitive Zassenhaus recombination argument.

## Occurrence-sensitive divisors reach the concrete scan

The selected-prime/Hensel divisor witness is now connected to the actual
generated fixed-size Zassenhaus enumerator.

- `FixedSizeScanExhausted` names the literal generated scan from its concrete
  iota initial combination while keeping the private well-founded rank
  implementation encapsulated.
- `FixedSizeScanExhausted.rejects` proves that exhaustion means every legal
  occurrence-index array was actually visited and rejected.
- `integer_divisor_candidate_rejected_of_scan_exhausted` combines this with
  the genuine divisor-to-sublist theorem.  For every integer divisor it
  retains a strictly increasing array of positions in the actual normalized
  Hensel array, its modular product association, and the concrete rejected
  execution.  Repeated equal factors remain distinct occurrences.
- `exactDivmodRaw_remainder_eq_empty_of_dvd_of_degree_lt` closes the algebraic
  end of actual sparse long division: for a concrete canonical returned
  remainder below the divisor degree, true divisibility forces that physical
  remainder array to be empty.  The remaining execution lemma must establish
  the degree bound and exclude the checked fault branches from canonical
  divisible inputs.

No C++ files changed, so this step has no new C++ b2b change surface.  The
remaining completeness boundary is now sharply localized: use the generated
precision/symmetric-recovery facts and exact long division to contradict this
rejection for a genuine bounded divisor.

## Actual exact division is complete on canonical divisible inputs

The generated sparse integer long-division recursion now has its missing
completeness theorem.

- `canonical_head_coefficient_dvd_of_poly_dvd` derives the literal C++
  leading-coefficient divisibility branch from polynomial divisibility and
  canonical nonempty sparse representations.
- `subtractScaledNormalize_divisionRank_lt` proves that the exact generated
  subtraction and normalization cancel the current leading term.  Its proof
  uses the concrete head degrees and integer quotient, then establishes the
  strict decrease of the same `divisionRank` used by generated well-founded
  recursion.  No fuel or alternate division implementation is introduced.
- `exactDivmodLoop_complete_of_dvd` inducts on that rank, discharges the
  nonempty, degree, nonzero-leading, divisibility, and decrease checks, and
  preserves divisibility through the actual recursive remainder update.
- `exactDivmodRaw_complete_of_dvd` therefore returns an actual generated
  quotient together with the physical empty sparse remainder for every
  canonical nonempty divisor that genuinely divides the canonical dividend.

No C++ files changed, so this step has no new C++ b2b change surface.  Next,
the concrete symmetric-recovery theorem supplies the recovered true divisor
to this exact-division result, ruling out rejection at its legal candidate.

## Resolve the non-monic Zassenhaus scaling exactly

The modular unit hidden by the divisor/subset `Associated` statement is now
identified algebraically instead of being carried as an arbitrary unit.

- `leading_scaled_monic_associated_divisor` proves that, for an actual
  factorization `source = divisor * quotient`, scaling the monic modular
  candidate by `source.leadingCoeff` gives exactly the reduction of
  `C quotient.leadingCoeff * divisor`.
- The quotient leading coefficient is essential: using
  `C source.leadingCoeff * divisor` would be false for a general non-monic
  integer factor and would not match the concrete C++ trial product.
- This lemma establishes the precise mod-`p` base equation required by
  `hensel_unique`.  The remaining step is to prove monicity/leading-term
  preservation for the factors returned by the actual generated Hensel tree,
  then lift this equality to its returned `p^k` modulus.

No C++ files changed, so this step has no new C++ b2b change surface.

## Preserve the generated modular-division stopping condition

The concrete modular long division used inside each Hensel correction now
exports the source-level fact that was previously lost at its semantic
boundary.

- `DivmodTrace.output_inactive` follows the actual finite generated trace to
  its `done` constructor and proves that the returned remainder is precisely
  a state where the C++ loop condition is false.
- `__upoly_divmod_mod_raw_ir_remainder_inactive_of_run` recovers that fact
  from a successful public raw execution; it inverts all checked source
  branches and uses the exact trace stored by the termination certificate.
- `__upoly_divmod_mod_raw_ir_remainder_get_deg_lt_of_run` specializes the
  stopping condition: every nonempty returned remainder has strictly smaller
  source degree than the divisor.  This is the missing executable premise for
  proving that the `h + m * remainder` Hensel update preserves its leading
  term.

No C++ files changed, so this step has no new C++ b2b change surface.  Next,
the strict degree result is propagated through the concrete scale/add/mod
operations and then through the recursive Hensel tree.

## Frame the concrete Hensel correction merge

The actual sparse addition loop now exposes the frame property needed to
carry a factor's leading term through a Hensel correction.

- `pairVecAddLoop_result_prefix` follows all six generated merge branches and
  proves that recursion only appends after the explicit result cursor.  Its
  equal-degree branch separately handles `pushNonzero` appending a sum and
  omitting an exact cancellation.
- `pairVecAddLoop_preserves_left_head` executes the first source iteration
  when every right-hand correction degree is below the concrete left head.
  It proves that the left head is emitted first, then uses the recursive frame
  theorem to show it cannot be overwritten.
- The proof explicitly reconciles generated `getElem!` reads with bounded
  `getElem` reads; no unchecked array identity is assumed.

No C++ files changed, so this step has no new C++ b2b change surface.  Next,
the modular-divmod degree stop is transported through `scaleCoeffs`, this
merge theorem, and the exact coefficient-reduction loop.

## Preserve a Hensel factor head through its concrete correction

The complete generated `h`-field correction fragment now retains a physical
coefficient-one head under the strict remainder-degree premise.

- `scaleCoeffs_degrees_below` proves directly from the generated array map
  that multiplying the correction by `m` changes no stored exponent.
- `modCoeffOutput_preserves_head` follows the generated filter-map and keeps
  a concrete head when its floor remainder is unchanged and nonzero;
  `modCoeffOutput_preserves_one_head` discharges those arithmetic checks for
  coefficient one and every modulus greater than one.
- `henselHCorrection_preserves_one_head` composes the exact scale, merge, and
  `m^2` coefficient-reduction operations in source order.  Its conclusion is
  an equality of the returned sparse array's physical list head, not only a
  congruence of decoded polynomials.

No C++ files changed, so this step has no new C++ b2b change surface.  Next,
the actual modular-divmod stopping theorem supplies the remainder-degree
premise, closing leading-term preservation for one generated Hensel step.

## Convert the concrete division stop to every sparse term

The generated modular-divmod stopping condition now supplies exactly the
termwise premise consumed by the concrete Hensel correction theorem.

- `get_deg_toInt_eq_head` identifies the signed `Int64` value returned by the
  source `get_deg` helper with the mathematical natural exponent stored in a
  bounded nonempty sparse head.
- `__upoly_divmod_mod_raw_ir_remainder_terms_below_of_run` combines that
  identity with the actual trace stop, the returned remainder's canonical
  ordering, and its machine-degree bounds.  Every physically stored
  remainder term is therefore strictly below the divisor's physical head.
- The result is derived from the exact raw call and its finite trace; no
  abstract Euclidean remainder is substituted for the generated output.

No C++ files changed, so this step has no new C++ b2b change surface.  Next,
the trace is shown to preserve machine-degree bounds to its returned
remainder, after which the theorem instantiates directly inside the generated
Hensel factor phase.

## Preserve the physical right-factor head through raw Hensel execution

The generated Hensel factor phase now has the executable bridge needed to
show that a monic right factor remains physically monic after one correction.

- `divmodLoop_remainder_degreesBound` proves by induction over the supplied
  well-founded `DivmodTrace` that every concrete division transition
  preserves the machine exponent bound.
- `__upoly_divmod_mod_raw_ir_remainder_degreesBound_of_run` inverts a
  successful generated raw call and transports that bound to its actual
  returned remainder; no abstract division result is introduced.
- `henselFactorCorrection_preserves_h_one_head_from_raw_runs` combines the
  raw multiplication canonicality, actual remainder bounds, physical degree
  stop, and the exact scale/add/mod run equations.  Its conclusion identifies
  the coefficient-one head in the returned sparse array itself.

No C++ files changed, so this step has no new C++ b2b change surface.
`lake build CLPoly.Refinement.Hensel` completed successfully (3328 jobs).
Next, this property is threaded through the generated step, recursive tree,
and extraction order so full-modulus Hensel uniqueness can be applied to the
real candidate subset.

## Carry physical monicity through the Hensel tree and precision loop

Physical right-factor monicity is now part of the strict step postcondition
and is propagated through both generated recursive control structures.

- `HenselStepCorrect` now records that every coefficient-one input `h` head
  survives the exact factor and Bézout phases.  The latter changes only `s`
  and `t`, so the physical factor-phase result is retained verbatim.
- `HasPhysicalOneHead` and `HenselArrayHOneHead` state the representation
  invariant directly on sparse arrays.  `henselArrayHOneHead_set` proves the
  exact mutable-array frame rule.
- `HenselLiftRecursiveCorrect.arrayHOneHead` follows all four finite tree
  constructors and applies the generated single-step theorem at the written
  node before each real child traversal.
- `HenselLiftLoopCorrect.arrayHOneHead` carries the invariant through every
  actual quadratic-precision iteration, deriving `1 < m` from the source
  loop's `2 ≤ m` measure and preserving it under `m := m * m`.

No C++ files changed, so this step has no new C++ b2b change surface.
`lake build CLPoly.Refinement.Hensel` completed successfully (3328 jobs).
Next, the canonical builder's initial right factors are proved physically
monic and the property is transferred through exact leaf extraction.

## Recover a literal monic head from generated field conversion

The generated tree builder's `Zp`-to-`ZZ` conversion now has the exact
representation theorem required to initialize the physical Hensel invariant.

- `henselTreeZpToZZIR_hasPhysicalOneHead_of_monic` starts from a canonical,
  nonempty sparse finite-field polynomial whose decoded polynomial is monic.
- It identifies the decoded leading coefficient with the actual array head,
  uses reduced-representative injectivity to prove that the stored `UInt64`
  value is literally one, and then follows the generated array `map` used by
  `poly_convert`.
- The conclusion is `HasPhysicalOneHead` for the concrete integer sparse
  array, rather than a modular equality or a postulated monic output.

No C++ files changed, so this step has no new C++ b2b change surface.
`lake build CLPoly.Refinement.Hensel` completed successfully (3328 jobs).
Next, monicity of the actual right-half product is established and this
conversion theorem is threaded through the recursive builder certificate.

## Prove monicity of the concrete builder interval product

`henselFactorRangeProduct_monic` now proves that every half-open product used
by the generated tree builder is monic whenever each concrete array element
read in that interval is monic.

The proof follows the same well-founded `stop - index` recursion as the
builder's mathematical denotation, performs the actual bounded array read at
each step, and composes `Polynomial.Monic.mul`.  The terminal interval is the
literal product identity.  This supplies the missing premise for applying
the generated `Zp`-to-`ZZ` physical-head theorem to each right-half product.

No C++ files changed, so this step has no new C++ b2b change surface.
The complete Hensel source file passed direct Lean checking.  Next, this
interval result is attached to each successful raw builder store and carried
in its recursive semantic certificate.

## Expose the physical monic head at the concrete builder store

`henselTreeStoreNodeRawIR_canonical_oneHead` now extracts two facts from the
same successful execution of the generated six-write tree-node store: the
stored `g`/`h` payload is canonical, and the converted right interval product
has a literal `(degree, 1)` array head.  The latter is obtained from the
actual `henselTreeZpToZZIR` output equality, not from polynomial extensionality
or a replacement conversion.

`HasPhysicalOneHead.of_algebraEq` transports this representation fact only
across the builder's proven algebraic frame.  Consequently later writes of
left/right child indices cannot silently discard it, while arbitrary changes
to `g` or `h` remain impossible to justify through this lemma.

No C++ files changed, so this step has no new C++ b2b change surface.  The
complete Hensel source passed direct Lean checking.  Next, the store result is
installed in `HenselTreeSemanticBuildCertificate` and recursively propagated
to every concrete topology node.

## Preserve selected-factor monicity after the concrete slot-zero write

The generated Hensel entry bakes the source leading coefficient into exactly
slot zero.  `HenselAdjustFirstFactorCorrect.getElem_eq_of_pos` now proves by
unfolding that actual `Array.set! 0` that every positive-index read is the
same read from the original selected-factor array.  The proof retains the
source bound used by that read; it does not weaken the claim to membership or
permutation.

`HenselAdjustFirstFactorCorrect.monic_of_pos` transports the monicity supplied
by the genuine `__factor_Zp`/SelectPrime certificate through this positional
equality.  Every right-half interval in the concrete Hensel builder starts at
a positive midpoint, so this is precisely the input needed by the interval
product theorem.

No C++ files changed, so this step has no new C++ b2b change surface.  The
complete Hensel source passed direct Lean checking.  Next, these positional
facts are supplied to the recursive raw builder and recorded at every stored
node.

## Carry physical right-factor heads through the recursive builder

`HenselTreeSemanticBuildCertificate.node` now records
`HasPhysicalOneHead value.h` for the concrete node found by its array lookup.
Both certificate transport operations retain this representation fact:
lower-bound weakening leaves the node untouched, while array-prefix
preservation carries the exact lookup equality.

The well-founded `strictHenselTreeBuildRecursiveRawIR_succeeds` proof now
requires monicity only for positive factor indices.  At each real source
recursive call it proves the right-half product monic with
`henselFactorRangeProduct_monic`; the midpoint is positive, so the adjusted
slot zero is never used by that product.  The same successful six-write store
then yields `hstoredOneHead`.  All four source branch shapes (both children,
left only, right only, neither) transport this fact through their actual
child-pointer writes and install it in the emitted semantic node.

The full Hensel entry supplies the new builder premise through the concrete
slot-zero adjustment theorem and a pointwise monicity field originating at
the selected `__factor_Zp` result.  No final Hensel output, tree, or semantic
factor is prescribed by this premise.

No C++ files changed, so this step has no new C++ b2b change surface.  The
complete Hensel source passed direct Lean checking.  Next, topology coverage
turns these per-node certificates into `HenselArrayHOneHead` for the exact
builder output.

## Cover the full concrete builder array with physical heads

`HenselTreeSemanticBuildCertificate.oneHead_of_mem` follows the exact
left/right certificate structure and proves that every index enumerated by
the source topology has a successful concrete array lookup carrying its
stored `HasPhysicalOneHead` fact.

`henselTreeBuildTopology_indices_complete` is the converse of the existing
boundedness theorem.  By well-founded recursion on `stop - start`, it proves
that the source preorder topology contains every index in the contiguous
allocation block `[root, root + internalNodeCount)`.  Its four cases match
the actual left/right allocation branches and use their exact count offsets.

The root builder contract combines completeness, the recursively proved
lookup facts, and the exact output-size equation to return
`HenselArrayHOneHead output`.  Pass 9 now generates this strengthened public
contract.  While synchronizing the template, the older generator drift that
would have reintroduced an abstract `DivmodTermination` argument was removed:
the generated step, recursive lift, and full Hensel entry continue to use the
proved `concreteDivmodTermination` implementation.

No C++ files changed, so this step has no new C++ b2b change surface.  Direct
checking of the complete Hensel source and the Pass 9 generated-file freshness
check passed.  Next, the builder-wide property is propagated through the
quadratic lift and exact leaf extraction.
