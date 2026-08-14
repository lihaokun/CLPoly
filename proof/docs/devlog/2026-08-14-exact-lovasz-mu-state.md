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
