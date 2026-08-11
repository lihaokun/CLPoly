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

## Generated terminal EEA scaling

The terminal extended-Euclidean branch no longer scales its three result
polynomials through an in-model helper. `HenselEEARawOps` now exposes the
actual raw scaling operation, and `__polynomial_GCD_eea_raw_ir` executes that
operation separately for the gcd, left Bézout coefficient, and right Bézout
coefficient, propagating every raw failure in source order.

The concrete operation maps every stored coefficient through the generated
`dense_upoly_zp_nmod_mul_ir` implementation and then performs the same sparse
normalization as the source representation. Its refinement theorem proves
reducedness of every emitted coefficient and the exact polynomial identity

`toPoly output = C (toZMod coefficient) * toPoly input`.

Consequently, the EEA execution trace records the three actual raw calls and
their returned arrays; it no longer manufactures terminal results with an L2
or model-level multiplication fallback. The remaining EEA work is to combine
these three scaling identities with the reachable Bézout invariant and prove
that the normalized terminal gcd is monic and extensionally equal to the L2
polynomial gcd.

## Unified divmod contract and concrete EEA termination

The empty-dividend, single-term-divisor, and general VHC paths are now
combined into one total contract for the complete source-order
`henselDivmodVHCIR` entry.  For canonical inputs and a nonempty divisor, the
actual raw call returns canonical quotient and remainder arrays, bounds every
remainder term strictly below the divisor lead, and proves
`quotient * divisor + remainder = dividend`.

A second execution-only theorem derives the strict remainder bound from every
successful call without assuming input canonicality or any semantic oracle.
This is exactly the universally quantified fact required by the generated EEA
termination interface.

The concrete `strictHenselEEARawOps` is therefore equipped with a genuine
well-founded measure: zero for an empty current second remainder and one plus
its leading degree otherwise.  Each successful recursive EEA transition
strictly decreases this measure because its next second remainder is the
remainder returned by the actual VHC division call.  There is no fuel counter,
expected remainder, L2 division fallback, or caller-supplied decrease proof.

Verification:

- `lake env lean CLPoly/Refinement/Hensel.lean`: successful.
- `git diff --check`: successful.
- placeholder and axiom scan of the strict Hensel generated/refinement/public
  modules: no `sorry`, `admit`, or `axiom` declaration.

## Complete remainder trace semantics

The exact quotient-and-remainder execution trace now proves coefficient
agreement for the entire final low-degree remainder.  The induction follows
the concrete selected frontier degree: iterations at or above the divisor
lead leave the remainder unchanged, while iterations below the lead identify
the consumed source coefficient with the residual and extend the processed
interval.  The terminal branch closes all degrees left below the final bound
from actual cursor and heap exhaustion.

## Concrete EEA raw-to-safe execution

The concrete VHC operation and its well-founded termination instance are now
fed into the generated `__polynomial_GCD_eea_raw_ir` recursion.  A reachable
state invariant proves that both remainder registers stay canonical and that
the previous remainder remains nonempty.  Its induction step consumes the
actual raw divmod result and the unified divmod contract; it does not choose a
quotient or remainder independently.

Consequently every nonterminal generated EEA state has an executable division
call, and every terminal state has the concrete leading term required by the
source monicization branch.  The resulting `strictHenselEEAIR_refines` theorem
is the concrete raw-to-safe execution bridge using the degree-based recursion.
It establishes execution and an exact source trace; the Bézout polynomial
semantics are intentionally the next obligation and are not claimed here.

## Reachable EEA Bézout invariant

The generated EEA assignment block now carries both algebraic identities
through every concrete reachable prefix:

- `r0 = s0 * left + t0 * right`;
- `r1 = s1 * left + t1 * right`.

For a successful source VHC result `(q, remainder)`, the proof uses its exact
division equation and the actual assignments `s2 = s0 - q*s1` and
`t2 = t0 - q*t1`.  Normalization is proved polynomial-preserving, and the
coefficient registers retain their well-formed runtime modulus.  Thus the
next remainder identity is derived by ring algebra from the concrete output,
not supplied as an EEA specification oracle.

The current polynomial bridge uses the existing bounded executable `Zp`
arithmetic lemmas (`2*p` and `p*p` within one limb).  These are arithmetic
preconditions, not expected-result fields; the tree-entry contract must
either obtain them from the selected Hensel prime or replace the legacy
bounded bridge with the newer overflow-safe `Zp` model theorem before the
public top-level theorem is emitted.

## EEA source initialization

The exact six-register initialization emitted for the source extended GCD is
now proved to establish both the execution invariant and the algebraic
invariant.  The singleton runtime polynomial created with `Zp.ofUInt64 1 p`
is proved canonical and denotes polynomial one, while the two empty arrays
denote zero.  Therefore the initial equations reduce definitionally to
`left = 1*left + 0*right` and `right = 0*left + 1*right`.

Before closing terminal monicization, the coefficient inverse boundary still
needs to be connected explicitly to the already verified generated C++
`inv_prime` implementation.  The current executable `Zp.inv` primitive is not
an L2 polynomial fallback or result oracle, but that primitive-level bridge is
required before publishing the final strict EEA contract.

## Generated terminal inverse execution

The EEA raw-operation interface now exposes the terminal inverse call instead
of evaluating `Zp.inv` invisibly inside the recursion.  The concrete Hensel
instance executes `Generated.StrictGCD.dense_upoly_zp_nmod_inv_ir`, i.e. the
generated C++ `inv_prime` control flow already verified by the strict word
arithmetic development.  The raw-to-safe trace records that exact successful
call and the inverse value used by all three source monicization loops.

For every reduced nonzero concrete leading coefficient, the new inverse
contract proves that the returned runtime coefficient is reduced modulo `p`
and that its `ZMod p` value multiplied by the actual leading coefficient is
one.  Thus terminal EEA normalization no longer relies on the unconnected
`Model.Zp.inv` implementation.

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

## Reachable general-state invariants

The general quotient-and-remainder state now has an explicit projection to
the quotient-only iteration result.  A successful double-output body projects
to exactly that structure, allowing all existing generated VHC preservation
theorems to apply to the same concrete selector/consume/emit/reinsert run.

A full operational invariant transports quotient canonicality, divisor-node
indices, state coverage, heap ownership/homogeneity/order, node denotation,
reset readiness, cursor prefixes, consumed dividend prefixes, processed
quotient prefixes, remaining-dividend bounds, active-node bounds, and
quotient readiness.  It is proved preserved by every actual double-output
body.  The fresh empty quotient/remainder and generated node array establish
the invariant initially, and induction over the concrete reachable-prefix
relation now supplies it for every state used by the semantic trace.  No field
contains an expected result or L2 division value.

## Entry-level general divmod contract

The fresh general five-argument entry now has its complete semantic contract.
A successful well-founded execution is converted directly into the exact
source trace; fresh-prefix operational invariants then prove every low-degree
remainder coefficient, while the quotient projection supplies the previously
verified high-degree product agreement.  Coefficient extensionality closes
`quotient * divisor + remainder = dividend`.

Both concrete outputs are canonical.  Quotient canonicality follows the same
generated VHC body projection.  Remainder canonicality is proved separately
along the actual trace: every source push uses the consumed strict-word
coefficient, its reducedness comes from the real modular consume operation,
the push condition supplies nonzeroness, and the selected frontier decrease
places the new term strictly below all prior remainder terms.  The final
remainder also retains the strict degree bound below the divisor lead.

Verification:

- `lake env lean CLPoly/Refinement/Hensel.lean`: successful.
- `git diff --check`: successful.
