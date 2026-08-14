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
