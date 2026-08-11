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

The descending remainder invariant is now explicit.  Stored remainder terms
are tracked above the next frontier bound, while their polynomial coefficients
are tracked against the residual on the already processed interval.  A proved
interval-extension lemma combines one exact conditional push, the selected
residual coefficient, and a zero-coefficient frontier gap to extend agreement
from the old bound down to the newly selected degree.  This is the induction
step needed by the complete below-leading-degree loop proof.

The residual gap is now discharged directly from the concrete VHC state in
both control-flow branches needed by the well-founded proof.  At a selected
frontier, the generated dividend cursor and product heap independently prove
that every skipped coefficient is zero; at loop termination, exhaustion of
the cursor and heap proves the same fact below the final degree bound.  Taking
their polynomial difference yields the required zero residual coefficients.
Neither theorem assumes an expected quotient/remainder nor invokes an L2
division operation.

Verification:

- `lake env lean CLPoly/Refinement/Hensel.lean`: successful.
- `git diff --check`: successful.

## Complete remainder trace semantics

The exact quotient-and-remainder execution trace now proves coefficient
agreement for the entire final low-degree remainder.  The induction follows
the concrete selected frontier degree: iterations at or above the divisor
lead leave the remainder unchanged, while iterations below the lead identify
the consumed source coefficient with the residual and extend the processed
interval.  The terminal branch closes all degrees left below the final bound
from actual cursor and heap exhaustion.

The operational invariant supplied for reachable prefixes contains only
canonicality and the generated cursor/heap/node-chain facts.  It contains no
expected output, quotient, remainder, or L2 division result.  A final
coefficient-extensional theorem combines this low-degree remainder result
with the already proved high-degree quotient product agreement and strict
remainder degree bound to obtain
`dividend = quotient * divisor + remainder`.

## Concrete general entry execution

The fresh-output general branch of the five-argument Hensel `pair_vec_div`
overload is now executable and proved total for canonical nonempty dividends
and divisors with a genuine tail.  Each successful quotient-only source body
lifts to the quotient-and-remainder body by unfolding the identical selector,
consume, emit, and reverse-reinsertion calls; the only additional state is the
conditional source remainder push.  Strong induction then lifts the complete
well-founded VHC run while proving that its concrete quotient is unchanged.

This establishes the actual general divmod execution boundary needed by EEA,
not a wrapper around an L2 quotient or remainder.  Empty-dividend and
single-term-divisor source branches remain to be composed into the complete
EEA operation.

## Complete source branch execution

The EEA divmod operation now follows the full five-argument source branch
order.  It rejects a zero divisor, clears both fresh outputs, returns two empty
arrays for an empty dividend, executes a well-founded source-order loop for a
single-term divisor, and otherwise enters the general VHC loop.  In the
single-term branch, field coefficient division has zero coefficient remainder
exactly as the C++ `__div(Zp &, Zp &, ...)` specialization specifies; terms
whose monomials are not divisible are appended unchanged to the polynomial
remainder.

The resulting concrete function has been installed as
`strictHenselEEARawOps.divmod`, and totality is proved for canonical inputs
with a nonempty divisor.  The next proof obligation is the uniform semantic
contract and strict remainder-degree decrease across these branches, which
will provide the actual well-founded EEA termination measure.

## Single-term branch semantics

The source-order single-term loop now has an exact quotient/remainder
partition theorem.  Its quotient projection is the existing checked
single-output loop; its remainder list is precisely the original input terms
whose monomials fail `is_divexact`.  A termwise proof using the actual strict
word inverse/multiply coefficient implementation establishes the full
polynomial equation
`quotient * divisor + remainder = dividend`.

The concrete remainder is also proved canonical and every stored term has
degree strictly below the divisor monomial.  Thus the complete semantic and
termination obligations for the single-term EEA division branch are closed;
the remaining uniformization work is concentrated in the fresh general VHC
initial state and its reachable operational invariants.

Verification:

- `lake env lean CLPoly/Refinement/Hensel.lean`: successful.
- `git diff --check`: successful.
