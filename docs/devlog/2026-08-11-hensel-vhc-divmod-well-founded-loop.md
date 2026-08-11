# Close the Hensel VHC divmod execution loop

The five-argument C++ `pair_vec_div` overload used by Hensel's extended
Euclidean algorithm now has a strict quotient-and-remainder VHC loop.  The
loop executes the same frontier selection, equal-degree heap consumption,
quotient emission, reset activation, reverse reinsertion, and low-degree
remainder push as the source implementation.

Termination is well-founded on the selected frontier degree.  Each recursive
call replaces the current degree limit with the degree returned by the actual
heap/dividend selector and requires that concrete degree to be strictly
smaller.  There is no fuel parameter and no L2 division call in the executable
definition.

The new raw-to-safe bridge proves that every state reachable through the
concrete iteration bodies produces a successful result together with an exact
source-level execution trace.  Its readiness invariant supplies only safety
and degree-decrease facts for actual raw calls; it contains no expected
quotient, expected remainder, or specification oracle.

This commit does not yet claim the complete quotient/remainder semantic
theorem.  The next step is to carry the quotient product-agreement invariant
and the accumulated low-degree remainder invariant through this trace, then
instantiate the strict Hensel EEA and tree builder with the resulting concrete
division operation.

Verification:

- `lake env lean CLPoly/Refinement/Hensel.lean`: successful.
- `lake build`: successful, 3498 jobs.
- placeholder/oracle/fallback scan of strict Hensel sources: no prohibited
  occurrence; the remaining `fuel` matches are explanatory statements that
  explicitly reject fuel-based recursion.
- `git diff --check`: successful.

## Quotient semantic projection

The quotient-and-remainder iteration and its complete well-founded loop now
project definitionally onto the already verified quotient-only VHC execution.
The projection follows each successful concrete selector, consume, emit, and
reinsertion call; it does not recompute or predict a quotient.

As a result, the five-argument loop inherits the existing theorem that the
computed quotient times the divisor agrees with the dividend at every degree
at or above the divisor's leading degree.  The remaining semantic obligation
is now isolated to the concrete low-degree terms accumulated in the remainder
array.

The remainder push has also been exposed as an exact theorem about every
successful iteration.  Induction over the concrete execution trace now proves
that, from an empty initial remainder, every emitted remainder term has degree
strictly below the divisor's leading degree.  This discharges the Euclidean
degree-bound half of the remainder contract; coefficient equality with the
low-degree residual remains to be proved.

The first coefficient-level residual bridge is now closed.  At every selected
frontier below the divisor leading degree, the coefficient consumed by the
actual heap loop is proved equal to the corresponding coefficient of
`dividend - quotient * divisor`.  The proof reuses the exact heap-owner
coefficient theorem from the quotient-only VHC development and separately
proves that `quotient * leadingMonomial(divisor)` cannot contribute below the
leading degree.  It also identifies the source conditional remainder push
with precisely that consumed coefficient.
