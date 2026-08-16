# Preserve physical nonempty van-Hoeij validation results

## What changed

`exactDivmodRaw_divisor_nonempty_of_success` unfolds the generated sparse
long-division entry and proves that any successful execution on a nonempty
dividend must have taken its positive-size divisor branch.  Thus a candidate
factor accepted by literal exact division cannot be an empty sparse array.

`validateCandidatesLoop_result_nonempty` follows all branches of the generated
candidate-validation loop.  Rejected, unavailable, empty, and trivial
candidates preserve the existing result invariant.  The successful branch
uses the actual trial product, symmetric recovery, primitive call, exact
division, quotient primitive call, consumed marking, and recursive push; the
new exact-division lemma proves that the physically pushed factor is nonempty.
`validateCandidates_result_nonempty` packages the invariant for the complete
generated validation entry.

## Why

The van-Hoeij cardinality closure requires every returned polynomial to be
nonzero and nonunit before equality with the modular irreducible-atom count can
force one normalized factor per output.  This step establishes the first of
those facts at the physical sparse-array level without replacing execution by
an abstract factorization witness.

## Verification and remaining work

- `CLPoly/Refinement/Recombine.lean` compiles directly.
- The theorem is included in `CLPoly.Refinement.AxiomAudit`.
- C++ changes: none, so there is no new C++ regression or B2B change surface.
- Next: preserve canonicality for every pushed factor, derive polynomial
  nonzeroness, and use the nonempty candidate's positive modular degree to
  rule out units.

## Recover canonical and nonzero L2 factors

`validationRecoveredFactor_canonical` now packages the exact successful
validation prefix.  Starting from the current canonical nonempty `fStar`, it
proves the physical leading-coefficient singleton canonical, transports this
through the actual trial-product loop, obtains positivity from the successful
`symmetricModRaw` branch, and then follows symmetric recovery and the literal
primitive call to the returned factor.

`sparsePolyZZ_toPoly_ne_zero_of_canonical_nonempty` proves the general sparse
representation bridge by comparing the physical head with the L2 leading
coefficient.  `validationRecoveredFactor_ne_zero` combines it with the actual
successful exact-division trace, whose divisor-size check supplies physical
nonemptiness.  Thus the precise candidate factor that the generated loop can
push is nonzero as an L2 polynomial.

- C++ changes: none, so this step has no new C++ regression or B2B change
  surface.
- Next: thread this certificate through the complete result array and prove
  positive modular degree/nonunit status.

## Thread L2 nonzeroness through complete validation

`validateCandidatesLoop_result_ne_zero` strengthens the live validation
invariant to carry both the canonical current quotient and pointwise L2
nonzeroness of the physical result array.  On the successful branch it applies
`validationRecoveredFactor_ne_zero` to exactly the factor being pushed,
constructs the pushed-array invariant occurrence-sensitively, and derives the
next quotient's canonical representation from the same exact-division and
primitive-normalization traces.  Every non-success branch reuses the unchanged
physical state.

`validateCandidates_result_ne_zero` specializes the invariant to the complete
generated entry.  The final result array needed by the normalized-factor count
therefore has a direct execution-derived nonzero certificate for every member.

- C++ changes: none, so this step has no new C++ regression or B2B change
  surface.
- Next: establish nonunit status using the positive degree of nonempty products
  of the selected modular irreducible factors.

## Prove an accepted physical validation factor is nonunit

`candidateAvailableLoop_true_valid_suffix` reverses the literal checked-index
scan: a `true` result proves that every remaining `Int32` index passed both
nonnegativity and physical-array bounds.  Its entry specialization produces
the exact `CandidateIndicesValid` certificate consumed by the trial-product
refinement.

For a nonempty candidate, `selectedProductMod_not_isUnit_of_nonempty_irreducible`
uses the actual first checked array occurrence.  That occurrence denotes an
irreducible modular factor and divides the exact selected list product, so the
product cannot be a unit.  `validationRecoveredFactor_mod_eq_selected` then
follows the same trial-product, symmetric-recovery, and primitive-content
trace at the prime divisor of the Hensel modulus.

Finally, `validationRecoveredFactor_not_isUnit` proves the physical recovered
integer factor is nonunit.  Nonzero source leading coefficient and the
nonzero selected irreducible product force the actual primitive content to be
nonzero modulo the prime.  Removing the two resulting constant units gives an
association between the mapped recovered factor and exact selected product;
an integer unit would map to a modular unit, contradicting that association.

- C++ changes: none, so this step has no new C++ regression or B2B change
  surface.
- Next: thread this nonunit invariant through all successful validation pushes
  and combine it with the result-count equality theorem.

## Preserve the selected-prime leading coefficient after validation

`validationQuotient_leading_mod_ne_zero` proves that the exact primitive
quotient installed after a successful validation extraction retains a nonzero
leading coefficient modulo the selected prime.  It uses
`successfulTrialExtraction_toPoly`, hence the same returned exact quotient and
the same literal quotient-primitive call.  Taking leading coefficients and
mapping the resulting integer equality modulo the prime shows that a zero next
leading coefficient would force the prior live leading coefficient to zero.

This closes the recursive state premise needed to reapply
`validationRecoveredFactor_not_isUnit` at every later successful candidate,
rather than only at the first validation push.

- C++ changes: none, so this step has no new C++ regression or B2B change
  surface.
- Next: package canonicality, nonempty live quotient, nonzero modular leading
  coefficient, consumed-array size, and pointwise nonunit results into the
  complete generated validation-loop invariant.

## Preserve nonunit results through the generated validation loop

`validateCandidatesLoop_result_nonunit` follows the literal generated
candidate-index recursion and carries the complete physical state needed by
the final van-Hoeij cardinality argument.  At each accepted exact division it
uses the actual trial product, symmetric recovery, primitive normalization,
exact quotient, quotient normalization, and consumed-marker update.  It proves
simultaneously that the live quotient remains canonical and physically
nonempty, its leading coefficient remains nonzero modulo the selected prime,
the consumed bitmap retains the active lifted-factor length, and every factor
actually pushed into the result array is nonunit over `ZZ`.

`validateCandidates_result_nonunit` instantiates that invariant at the public
generated `validateCandidates` entry, including its literal false-filled
consumed bitmap and active-factor count.  Neither theorem replaces execution
with an abstract factorization or existence oracle; both invert and recurse
over the generated `RawExec` control flow.

- Termination measure: `candidates.size - candidateIndex`, inherited directly
  from the generated well-founded validation loop.
- C++ changes: none, so this step has no new C++ regression or B2B change
  surface.
- Next: apply the invariant to the validation call inside the generated
  van-Hoeij iteration, combine nonzero/nonunit output membership with exact
  output cardinality and associated products, and close irreducibility through
  the literal retry loop.

## Associate every recovered factor with its checked modular subset

`validationRecoveredFactor_mod_associated_selected` strengthens the existing
modular equation for one physical validation success.  It proves that mapping
the returned primitive integer factor to the selected prime is associated to
the exact `SelectedProductMod` traversed by `trialProductLoop`.  The proof
derives nonzeroness of that selected product from the actual checked indices
and pointwise modular irreducibility, then proves both scalar constants in the
generated recovery equation are units.

`validationRecoveredFactor_not_isUnit` now consumes this association rather
than reconstructing it internally.  This exposes precisely the relation needed
to turn consumed-subset cardinality into memberwise irreducibility later.

- C++ changes: none, so this step has no new C++ regression or B2B change
  surface.
- Next: preserve these modular subset associations across all validation
  pushes, together with the consumed partition, and apply the equal-cardinality
  UFD bridge on the accepted `__lll_factorize` branch.

## Preserve modular nonunits and lift the cardinality conclusion

The complete validation invariant now records three properties for every
physical result factor: it is a nonunit over `ZZ`, its reduction at the
selected prime is a nonunit, and its integer leading coefficient remains
nonzero after reduction.  The modular nonunit proof uses the newly exposed
association with the exact checked candidate product.  The leading fact comes
from `validationRecoveredFactor_leading_mod_ne_zero`, which takes leading
coefficients of the actual exact-division/primitive-quotient equation; a zero
factor leading coefficient would force the previous live leading coefficient
to zero.

`modular_irreducible_members_of_equal_length_associated_product` packages the
UFD cardinality step for a physical output array.  Given an associated product
with the irreducible modular atoms and equal physical length, the preserved
leading coefficients establish nonzero mapped outputs and the validation
invariant supplies nonunits, so every mapped output is irreducible.

Finally, `factorArrayIrreducible_of_modular` lifts those memberwise modular
irreducibility facts back to the actual primitive integer factors.  Each
physical member is primitive because it divides the primitive accumulated
product, and the preserved leading coefficient prevents degree loss modulo the
prime.

- C++ changes: none, so this step has no new C++ regression or B2B change
  surface.
- Next: preserve the global modular product and output-cardinality equations
  through `vanHoeijLoop`, then discharge the accepted first-pass branch of the
  generated `__lll_factorize_raw_ir`.

## Combine the Hensel certificate with equal-cardinality recombination

`factorArrayIrreducible_of_hensel_cardinality` is the cross-stage bridge for
an accepted `__lll_factorize` result.  It maps the actual integer recombination
product association to `ZMod p[X]`, composes it with the concrete
`LiveHenselProduct.primeProductAssociated` certificate for the same lifted
array, and obtains an associated product between physical output reductions
and the original irreducible modular atoms.

When the physical output size equals the lifted array size, the modular UFD
cardinality theorem makes every mapped output irreducible.  The theorem then
uses the primitive accumulated integer product and surviving leading
coefficients to lift that conclusion back to `FactorArrayIrreducible` for the
same returned array.

- C++ changes: none, so this step has no new C++ regression or B2B change
  surface.
- Next: derive this theorem's memberwise and cardinality premises from the
  literal `vanHoeijLoop` result and the first-pass acceptance guard.

## Account exactly for validation pushes and active removal

`markConsumedLoop_count_mono` proves that the literal generated marker loop
never clears a true bit.  Its strict companion uses the actual successful
`candidateAvailable` scan and a physical nonempty candidate to show that at
least one previously-false bit becomes true.

Threading that fact through every branch of the generated validation recursion
gives `validateCandidatesLoop_result_count_le_consumed`: the number of physical
result pushes is bounded by the increase in true consumed bits.  The public
entry specializes the initial false bitmap, so every newly returned factor is
charged to a concrete consumed position.

Independently, `removeConsumedLoop_size_add_count` follows the real reverse
erase loop and proves the exact equation between output active size and true
bits in the processed prefix.  At the public entry the whole bitmap is
processed.  Combining both sides yields
`validateRemove_result_active_size_le`: one successful van-Hoeij extraction
round cannot increase `result.size + active.size`.

- Termination measures: candidate suffix length for marking and validation;
  processed-prefix length for reverse removal.  No fuel is used.
- C++ changes: none, so this step has no new C++ regression or B2B change
  surface.
- Next: lift this per-round conservation law through the lexicographically
  well-founded `vanHoeijLoop`, including the literal Zassenhaus fallback size
  bound.

## Bound the literal Zassenhaus fallback output

`finishZassenhaus_size_le` proves that the shared terminal block either keeps
the physical result array or appends exactly one remaining positive-degree
factor; the final merge sort preserves length.

`zassenhausLoop_result_size_le` then follows the generated lexicographic
Zassenhaus recursion.  Subset exhaustion only advances the bounded subset
size.  A successful extraction uses the exact
`removeCombinationLoop_size` equation, so removing a positive candidate pays
for the one physical result push and leaves a positive active array under the
source `2 * subsetSize ≤ active.size` guard.  At termination the remaining
positive active slot pays for the optional final append.  Thus fallback output
size is bounded by `result.size + active.size` without any semantic
factorization assumption.

- Termination measure: the generated pair
  `(active.size, active.size + 1 - subsetSize)`.
- C++ changes: none, so this step has no new C++ regression or B2B change
  surface.
- Next: prove validation's source `remaining` guard leaves at least one active
  index, then lift the common size budget through `vanHoeijLoop`.

## Preserve one physical active factor through validation

`markConsumedLoop_count_le` follows the exact generated marker recursion and
bounds the number of newly true bitmap positions by the unprocessed candidate
suffix.  `validateCandidatesLoop_consumed_count_lt` then follows every branch
of the real validation loop.  Its successful-factor branch uses the source
guard `candidate.size < remaining`, charges the concrete marker update, and
recurses with the literal subtraction `remaining - candidate.size`.

At the public entry, the initial bitmap is physically all false and
`remaining = activeLifted.size`, giving a strict bound on the final true-bit
count.  Combining that bound with the exact reverse-erasure equation proves
`validateRemove_active_nonempty`: a successful extraction round cannot erase
the last active factor.  This is execution-derived state needed by the next
well-founded van-Hoeij call, not an external termination or semantic oracle.

- Termination measures: candidate suffix length for both generated recursive
  loops; no fuel is used.
- C++ changes: none, so this step has no new C++ regression or B2B change
  surface.
- Next: carry `result.size + active.size` through every constructor of the
  generated lexicographic `vanHoeijLoop` and its literal fallback.

## Carry the physical factor budget through van Hoeij

`vanHoeijLoop_finished_size_le` follows the complete generated lexicographic
recursion with the invariant `result.size + active.size ≤ budget`.  A
successful validation/removal round uses the exact conservation theorem and
the newly proved active nonemptiness before taking the smaller-active
recursive edge.  A precision retry preserves the same two physical arrays.
The fallback branch executes the literal `zassenhausLoop`, uses its physical
output bound, and accounts exactly for the subsequent source-order append.

Specializing the invariant to the generated entry initialization proves
`__vanhoeij_recombine_raw_ir_size_le`: every successful concrete C++
recombination output has at most as many slots as the lifted Hensel array.
This is the missing arithmetic fact needed to interpret the generated
`__lll_factorize` result-count guard without assuming recombination
correctness as an oracle.

- Termination measure: the generated pair
  `(active.size, precisionRank target initial maximum)`; no fuel is used.
- C++ changes: none, so this step has no new C++ regression or B2B change
  surface.
- Next: combine this upper bound with the literal `__lll_factorize` acceptance
  guard and the Hensel cardinality certificate, then lift modular
  irreducibility to the returned integer factors.

## Interpret the literal signed result-count guard

`size_toUInt32_toInt32_lt_iff` proves that the exact C++ conversion
`size_t → uint32_t → int32_t` reflects natural strict comparison below the
signed 32-bit boundary.  Combined with the concrete recombination upper bound,
`result_size_eq_of_not_machine_lt` shows that failure of the generated
`result.size() < r` test means exact equality with the Hensel factor count.
No mathematical comparison is substituted for the machine comparison.

- C++ changes: none, so this step has no new C++ regression or B2B change
  surface.
- Next: analyze the generated conjunction with `aH < aMig`; the equal-count
  acceptance branch uses modular cardinality, while every full-precision path
  must carry irreducibility through the actual van-Hoeij validation/fallback
  execution.

## Classify the generated two-pass LLL control flow

`__lll_factorize_raw_ir_low_precision_cases` unfolds the exact generated
entry after naming the concrete heuristic, first Hensel, and first van-Hoeij
results.  Under the literal low-precision condition `aH < aMig`, a successful
return is now proved to be exactly one of two source paths:

- the first result was accepted, in which case the signed count guard and the
  physical recombination upper bound force equality with the Hensel count; or
- the result is returned by the concrete second `henselLift(..., 0)` call and
  its immediately following van-Hoeij call.

`heuristic_starting_precision_first_le_second` separately follows every
branch of the generated heuristic (including its floating calculation and
overflow guard) and derives `aH ≤ aMig` from the literal final `min`.  Thus
failure of `aH < aMig` identifies a first pass already at full precision.

- C++ changes: none, so this step has no new C++ regression or B2B change
  surface.
- Next: establish the full-precision irreducibility invariant for factors
  physically appended by van-Hoeij validation, and combine it with the
  already-certified Zassenhaus fallback.

## Preserve the literal LLL basis transform end to end

`LLLTransformRel initial matrix transform` records the exact integer matrix
equation `matrix = transform * initial` on the full physical square arrays.
The proof starts from the generated diagonal initializer: at scale one its
transform is the identity matrix.  It then follows every mutation performed
by the generated C++ lowering:

- `subtractMatrixRows_preserves_transform_rel` identifies synchronized basis
  and transform row subtraction with left multiplication by the same integer
  transvection;
- `sizeReduceAt_preserves_transform_rel` and
  `extraSizeReduceLoop_preserves_transform_rel` carry that equation through
  the zero-rounding branch and the exact descending recursive loop;
- `swapMatrixRows_pair_preserves_transform_rel` proves that the two literal
  sequential row-swap calls preserve the equation;
- `lllStep_preserves_transform_rel` covers both the successful Lovasz branch
  and the swap branch, including all physically executed Gram--Schmidt update
  loops;
- `concreteLLLMainLoop_preserves_transform_rel` uses the existing determinant
  and index rank of the genuine well-founded generated loop; and
- `concreteLLLReduce_transform_rel` follows initialization, the entire main
  loop, and short-row collection to the actual generated `__lll_reduce`
  return value.

Consequently the candidate-column partition now receives a transform that is
proved to encode the returned lattice basis, rather than an unrelated array
whose dimensions merely happen to agree.  No fuel, partial definition,
existence witness, semantic LLL oracle, `sorry`, or custom axiom is used.

- Termination measure: the existing concrete LLL determinant/index rank for
  the main loop, plus literal decreasing natural indices for size reduction;
  no fuel is used.
- C++ changes: none, so this step has no new C++ regression or B2B change
  surface.
- Next: combine the transform equation with the short-row bound and the CLD
  lattice construction to prove that the generated candidate equivalence
  classes cannot merge distinct integer irreducible factors at sufficient
  precision.

## Identify generated CLD derivative and coefficient observations

The first two semantic observations used by the physical CLD lattice are now
connected directly to their generated loops:

- `derivativeZZRaw_toPoly` follows the decreasing source-array suffix of
  `derivativeZZLoop`, including the constant-term skip and zero-coefficient
  compaction branches, and proves that its concrete sparse output denotes the
  mathematical polynomial derivative;
- `sparseCoeff_eq_coeff` follows the generated linear lookup and uses strict
  canonical degree order to prove that every value inserted into a CLD matrix
  row is the corresponding mathematical polynomial coefficient.

These results prevent the later lattice argument from silently replacing the
generated sparse derivative or coefficient accessor with specification-level
operations.  Both executable loops are well-founded on the remaining array
suffix; no fuel, partial definition, semantic oracle, `sorry`, or custom axiom
is introduced.

- C++ changes: none, so this step has no new C++ regression or B2B change
  surface.
- Next: combine generated modular division, derivative, modular
  multiplication, and symmetric recovery into a pointwise semantic theorem
  for every physical element returned by `cldPolys`.

## Refine the generated CLD polynomial array

The complete physical CLD construction is now related to the mathematical
logarithmic-derivative data used by the van Hoeij argument.  In particular:

- `generatedModularQuotient_eq_divByMonic` starts from the successful return
  of generated modular sparse division, uses its physically empty remainder,
  and identifies the quotient with monic polynomial division after reduction
  modulo the requested modulus;
- `generatedCldElement_toPolyMod` follows the exact generated sequence of
  modular division, sparse differentiation, coefficient-wise reduction,
  modular multiplication, and symmetric coefficient recovery, proving that
  the returned sparse value represents the corresponding mathematical CLD
  polynomial modulo the modulus;
- `cldPolysLoop_toPolyMod` follows every recursive push performed by the
  generated array loop while maintaining the physical position-to-factor
  correspondence; and
- `cldPolys_toPolyMod` exposes the entry-point result: its output has exactly
  the active-factor cardinality and every physical array element denotes the
  matching mathematical CLD polynomial modulo the modulus.

This closes the semantic gap between the generated C++ CLD arithmetic and
the polynomial objects that will populate the lattice.  The recursive proof
uses the same well-founded remaining-active-array suffix as the generated
loop.  It introduces no fuel, partial definition, existence witness,
semantic oracle, `sorry`, or custom axiom.

- C++ changes: none, so this step has no new C++ regression or B2B change
  surface.
- Next: prove that the generated CLD data-row and matrix builders store these
  exact coefficients, then combine that matrix semantics with the concrete
  LLL transform equation and short-row bound.

## Refine generated CLD data-row entries

The first matrix-construction layer now has value semantics, not only shape
semantics:

- `fillCldDataRowLoop_entry` follows the generated indexed `Array.set` loop
  and proves that every populated coordinate is exactly the generated sparse
  coefficient lookup for the corresponding physical CLD polynomial; and
- `appendCldColumn_data_entry` follows zero-column extension, row filling,
  identity-coordinate insertion, and the final physical matrix push, proving
  that every CLD coordinate of the newly appended row contains that exact
  coefficient.

The proof uses the generated loop's well-founded measure
`cld.size - index`; it uses no fuel, partial definition, semantic matrix
oracle, `sorry`, or custom axiom.

- C++ changes: none, so this step has no new C++ regression or B2B change
  surface.
- Next: lift the one-column result through `buildCldMatrixLoop`, preserving
  old columns and identifying every spiral-selected degree in the complete
  physical lattice.
