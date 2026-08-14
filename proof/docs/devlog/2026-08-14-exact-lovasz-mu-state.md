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
