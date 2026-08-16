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

## Refine the complete generated CLD matrix build

The single-column coefficient theorem is now lifted through the complete
generated bounded matrix builder:

- `appendCldColumn_old_entry` proves that extending the square lattice by one
  row and column preserves every pre-existing physical matrix entry;
- `CldMatrixDataRows` records the exact correspondence between each appended
  physical row, its spiral-selected polynomial degree, and every CLD factor
  coordinate;
- `appendCldColumn_preserves_data_rows` extends that invariant by one genuine
  generated column-building step;
- `buildCldMatrixLoop_data_rows` follows both generated stopping branches and
  every successful append through the well-founded `target - added` loop; and
- `buildCldMatrix_data_rows` exposes the invariant at the public generated
  matrix-builder entry point.

Thus the lattice consumed by generated LLL is no longer known merely to be a
nonsingular square matrix of the expected size: all CLD data coordinates of
all appended rows are tied to the actual generated sparse coefficient
lookups at the exact alternating low/high spiral degrees selected by C++.
There is no fuel, partial definition, semantic matrix oracle, `sorry`, or
custom axiom.

- C++ changes: none, so this step has no new C++ regression or B2B change
  surface.
- Next: combine the full matrix data invariant with `cldPolys_toPolyMod`,
  symmetric coefficient bounds, `concreteLLLReduce_transform_rel`, and the
  short-row filter to obtain the candidate-partition separation theorem.

## Identify generated candidate-column equality

The physical comparator used to form the candidate partition now has an exact
logical characterization:

- `candidateColumnsEqual_true_sound` follows every generated bounds check and
  recursive selected-row access, proving that a successful `true` result means
  the two transform columns agree on every physical short-row index; and
- `candidateColumnsEqual_true_complete` proves the converse, including that
  the generated traversal executes successfully when all selected physical
  rows and columns are valid and equal.

Both directions recurse on the literal remaining short-row suffix
`shortRows.size - index`.  No fuel, abstract equivalence oracle, `sorry`, or
custom axiom is used.

- C++ changes: none, so this step has no new C++ regression or B2B change
  surface.
- Next: lift this exact comparator meaning through class assignment and the
  outer candidate-column partition, then use the CLD/LLL separation result to
  show that each produced class corresponds to one integer factor block.

## Refine generated candidate-class assignment

The exact comparator meaning is now carried through the generated inner class
assignment loop:

- `CandidateClassRepresentativeSound` states that every physical array slot
  carrying a given class identifier has passed the generated comparator
  against that class's representative column;
- `assignCandidateClass_preserves_representative_sound` proves this invariant
  across occupied slots, unequal comparisons, successful equal comparisons,
  and the terminating branch;
- `assignCandidateClass_size` proves the physical class array is never resized;
  and
- `assignCandidateClass_member_provenance` proves that every member visible in
  the returned class either already had that identifier or was inserted only
  after an actual `candidateColumnsEqual = .ok true` execution.

The proof follows the literal `factorCount - right` well-founded recursion and
every generated error/bounds branch.  It uses no fuel, abstract partition
oracle, `sorry`, or custom axiom.

- C++ changes: none, so this step has no new C++ regression or B2B change
  surface.
- Next: establish the corresponding provenance and representative invariant
  for `partitionCandidateColumns`, starting from its concrete replicated
  `none` array.

## Prove generated outer candidate partition coverage

The outer physical candidate partition now has its first complete execution
invariant:

- `assignCandidateClass_preserves_some` proves that the inner loop never
  overwrites a class identifier already present in the physical array;
- `partitionCandidateColumns_all_assigned` carries an assigned-prefix
  invariant through the occupied-column and new-class branches of the exact
  generated outer loop;
- `partitionCandidateColumns_size` proves that neither the outer loop nor any
  nested assignment changes the class-array cardinality; and
- `partitionCandidateColumns_from_empty_all_assigned` specializes the result
  to the actual `Array.replicate factorCount none` entry state and proves that
  every factor column in the successful output has a concrete class.

Both loops use their literal decreasing suffix measures.  There is no fuel,
partial definition, abstract partition oracle, `sorry`, or custom axiom.

- C++ changes: none, so this step has no new C++ regression or B2B change
  surface.
- Next: attach representative-comparison provenance to every outer class and
  carry it through `collectCandidateClasses`/`extractCandidates`, then apply
  CLD/LLL separation to identify the resulting integer factor blocks.

## Prove generated candidate collection coverage

The assigned physical class array is now connected to the concrete candidate
arrays returned by `extractCandidates`:

- `CandidateMember` records an actual physical position and `Int32` value in
  a returned candidate class;
- `collectCandidateClasses_preserves_member` proves that all later generated
  `Array.set`/`push` iterations preserve existing members;
- `collectCandidateClasses_includes_source` follows the exact collector until
  a requested source column, proves its class lookup and append, and carries
  that member through the remaining well-founded suffix; and
- `extractCandidates_covers_columns` combines outer partition coverage with
  the literal collector to prove that, in the nonempty-short-row branch, every
  physical factor column is materialized as its exact
  `source.toUInt32.toInt32` value in some returned candidate array.

This is execution provenance, not merely a returned-size estimate.  The proof
uses the generated `classes.size - column` recursion and introduces no fuel,
partial definition, collection oracle, `sorry`, or custom axiom.

- C++ changes: none, so this step has no new C++ regression or B2B change
  surface.
- Next: prove the converse member soundness and representative equality for
  each collected class, then combine it with the CLD/LLL short-vector
  separation theorem.

## Prove generated extracted-member soundness

Candidate collection now has the converse execution guarantee as well:

- `candidateMember_set_push_provenance` analyzes the exact nested
  `Array.set`/`Array.push` mutation: a visible member was either preserved from
  the old candidate matrix or is precisely the newly appended value;
- `collectCandidateClasses_member_provenance` lifts that local fact through
  the full generated collector and returns the concrete scanned source column,
  its physical class lookup, and its exact `Int32` conversion;
- `collectCandidateClasses_from_empty_member_source` removes the impossible
  initial-member case for the actual replicated-empty entry state; and
- `extractCandidates_members_sound` proves that every member of every
  candidate returned by the nonempty-short-row branch is the conversion of a
  real source column strictly below `factorCount`.

Together with the previous coverage theorem this excludes both dropped source
columns and fabricated candidate members.  All recursion is the generated
well-founded suffix recursion; there is no fuel, partial definition, semantic
collection oracle, `sorry`, or custom axiom.

- C++ changes: none, so this step has no new C++ regression or B2B change
  surface.
- Next: carry representative-comparison certificates through the outer
  partition so members of each returned candidate class are proven equal on
  the concrete short rows, then apply CLD/LLL separation.

## Prove generated candidate-column equivalence laws

The executable column comparator now supports the equivalence reasoning needed
for whole candidate classes:

- `CandidateColumnsValid` records the actual selected-row and factor-column
  bounds required by the generated array accesses;
- `candidateColumnsEqual_refl` executes the comparator on one physical column
  and proves a successful `true` result;
- `candidateColumnsEqual_symm` derives the reversed generated execution from
  the exact pointwise meaning of a successful comparison; and
- `candidateColumnsEqual_trans` composes two successful generated comparisons
  pointwise and executes the comparator for the outer pair.

These are theorems about the literal generated comparator, not an independently
postulated equivalence relation.  They use its well-founded short-row suffix
traversal and introduce no fuel, semantic equivalence oracle, `sorry`, or
custom axiom.

- C++ changes: none, so this step has no new C++ regression or B2B change
  surface.
- Next: use these laws with class-member provenance to prove pairwise equality
  for every generated candidate class, then connect that fact to CLD/LLL
  separation.

## Prove generated candidate classes are pairwise equal on short rows

The complete physical partition and collection pipeline now carries the
generated comparator certificate through to returned candidate members:

- `assignCandidateClass_other_member_origin` proves the inner loop cannot
  create a member of any class other than the class currently being assigned;
- class-id bounds, representative soundness, and pairwise comparator equality
  are established for the literal fresh-class `Array.set` operation and
  preserved through `assignCandidateClass`;
- `partitionCandidateColumns_preserves_pairwise` carries these invariants
  through the occupied-column and fresh-class branches of the generated outer
  loop;
- `partitionCandidateColumns_from_empty_pairwise` specializes them to the
  actual replicated-`none` initial array; and
- `extractCandidates_class_members_equal` traces any two physical members of
  one returned candidate array back to real source columns below
  `factorCount`, recovers their exact `Int32` values, and proves that the
  generated `candidateColumnsEqual` execution for those source columns
  returns `.ok true`.

Thus returned candidate classes are no longer justified merely by coverage or
membership bounds: their defining equality on every concrete selected short
row is proved end to end.  All loops use their generated decreasing suffix
measures; no fuel, partial definition, partition/equivalence oracle, `sorry`,
or custom axiom is used.

- C++ changes: none, so this step has no new C++ regression or B2B change
  surface.
- Next: derive `CandidateColumnsValid` from the concrete LLL return and short
  row collector, then combine pairwise transform-column equality with the CLD
  lattice semantics and sufficient-precision separation theorem.

## Validate generated LLL candidate columns end to end

The candidate-column safety hypothesis is now discharged by the literal
generated `__lll_reduce` execution:

- `LLLTransformRel` has been strengthened to retain the physical row length
  of every transform row in addition to the exact basis equation;
- initialization obtains that shape from the generated identity matrix, while
  the literal synchronized subtraction and Lovasz row-swap branches preserve
  it through every well-founded LLL step and the complete main loop;
- `insertShortRow_indices_valid` and `collectShortRows_indices_valid` prove
  directly from the generated array loops that every returned short-row index
  names a physical reduced-matrix row; and
- `concreteLLLReduce_candidate_columns_valid` combines both results for the
  actual transform and short-row array returned by `lllReduce`, proving all
  row and factor-column accesses needed by candidate partitioning are in
  bounds.

The added row-shape field is essential: the matrix equation alone uses
`getElem!` and therefore could not rule out a physically short transform row
whose missing entries observationally default to zero.  The proof now rules
that case out from the generated execution itself.  It introduces no fuel,
partial definition, shape oracle, `sorry`, or custom axiom.

- C++ changes: none, so this step has no new C++ regression or B2B change
  surface.
- Next: combine the certified physical candidate classes with the CLD lattice
  observations and the sufficient-precision separation argument, then close
  the generated van-Hoeij result theorem.

## Preserve the concrete LLL witnesses through candidate extraction

The LLL/candidate boundary now has a composed execution theorem rather than a
pair of disconnected component facts:

- `concreteLLLReduce_extractCandidates_class_members_equal` takes the exact
  subtype returned by generated `lllReduce` and the immediately following
  successful generated `extractCandidates` call; for any two physical values
  in one returned class it recovers their real source columns and the
  successful generated comparator execution on the same transform and short
  rows; and
- `concreteLLLReduce_extractCandidates_class_members_pointwise` opens that
  comparator execution and exposes equality of the two actual integer
  transform entries at every selected physical short row.

Keeping the transform and short-row array as execution witnesses is important
because `prepareCandidates` intentionally discards both after constructing
the candidate arrays.  The next theorem can now unfold that wrapper without
inventing replacement lattice evidence.  No fuel, partial definition,
semantic candidate oracle, `sorry`, or custom axiom is introduced.

- C++ changes: none, so this step has no new C++ regression or B2B change
  surface.
- Next: preserve these witnesses through both `prepareCandidates` branches
  and connect transform-entry equality to the generated CLD matrix rows.

## Prove generated short-row collection is complete

The exact short-row collector now has the completeness direction required by
the forthcoming lattice-span argument:

- `insertShortRow_preserves_member` and
  `insertShortRow_contains_inserted` follow every recursive comparator call
  and prove that sorted insertion neither loses an old physical row index nor
  omits the newly inserted one;
- `collectShortRows_preserves_member` carries membership through all later
  generated collection iterations;
- `collectShortRows_includes_bounded_row` proves that any concrete matrix row
  whose literal generated `dotRows` result is at most the C++ bound occurs in
  the returned array; and
- `concreteLLLReduce_includes_bounded_row` specializes that fact to the
  reduced matrix and short-row array returned by the same full generated
  `lllReduce` execution.

Together with the previous index-validity theorem, short rows are now tied to
the actual reduced basis in both directions needed downstream: no fabricated
indices, and no omitted bounded basis row.  The recursion uses the generated
array suffix measures and introduces no fuel, partial definition, collection
oracle, `sorry`, or custom axiom.

- C++ changes: none, so this step has no new C++ regression or B2B change
  surface.
- Next: prove the final generated LLL state satisfies its completed
  size-reduction/Lovasz conditions, then derive that every sufficiently short
  lattice vector lies in the span of the collected rows.

## Expose the genuine LLL loop exit condition

`concreteLLLMainLoop_finished` now follows the same determinant-ranked
well-founded recursion as the generated main loop and proves that every
successful returned state satisfies `matrix.size <= k`.  The only base case
is the literal negation of the C++ `k < matrix.size` guard; recursive results
must come from an actual successful `lllStep` and its proved rank decrease.

This is the exit fact needed to turn a preserved reduced-prefix invariant into
a reduced property for every adjacent basis pair.  It is not postulated on
the LLL output and introduces no fuel, partial definition, termination oracle,
`sorry`, or custom axiom.

- C++ changes: none, so this step has no new C++ regression or B2B change
  surface.
- Next: define and preserve the exact size-reduced/Lovasz prefix processed by
  the C++ `k` cursor, then specialize it to all rows using this exit theorem.

## Preserve the generated Lovasz witness through full size reduction

The advancing branch of the generated LLL step now carries a literal
post-state certificate for the Lovasz comparison it executed:

- `sizeReduceAt_preserves_mu_entry_of_source_lt` proves directly from the
  generated `sizeReduceAt` result that reducing against a lower source column
  cannot modify any higher Gram--Schmidt coefficient;
- `extraSizeReduceLoop_preserves_mu_entry_above` follows the genuine
  descending well-founded source loop and preserves the selected higher
  coefficient through every successful physical reduction; and
- `lllStep_advanced_lovasz_at_previous` combines that preservation with the
  exact successful C++ branch test and the unchanged norm array.  Thus the
  returned advanced state satisfies the same Lovasz inequality at the row
  just processed by the C++ `k` cursor.

This closes a subtle execution gap: using the branch proposition alone would
not certify the returned state unless the later `k-2 .. 0` reductions were
proved not to change `mu[k][k-1]`.  The proof uses the generated recursion and
physical arrays, with no fuel, partial definition, semantic LLL oracle,
`sorry`, or custom axiom.

- C++ changes: none, so this step has no new C++ regression or B2B change
  surface.
- Next: preserve the accumulated reduced-prefix predicate across both the
  advancing and swapping branches, then use the proved main-loop exit guard
  to obtain the complete final LLL certificate.

## Accumulate the generated Lovasz prefix on advancing steps

The generated C++ cursor now has an explicit processed-prefix invariant:

- `sizeReduceAt_preserves_mu_row_of_ne` proves a physical reduction changes
  only its current target row;
- `extraSizeReduceLoop_preserves_mu_row_before` propagates unchanged earlier
  rows through the actual descending well-founded loop;
- `lllStep_advanced_preserves_mu_row_before` composes the adjacent reduction,
  the later reductions, and the final cursor increment; and
- `lllStep_advanced_preserves_lovaszPrefix` extends `LovaszPrefix` by exactly
  the pair checked in the successful generated branch while preserving every
  previously checked pair.

The invariant is stated on the returned physical `mu` and norm arrays, not on
an abstract L2 LLL execution.  It therefore prevents a successful comparison
from being reused after a later array update has invalidated its operands.
No fuel, partial definition, semantic LLL oracle, `sorry`, or custom axiom is
introduced.

- C++ changes: none, so this step has no new C++ regression or B2B change
  surface.
- Next: show a failed-Lovasz swap preserves the strictly earlier prefix; then
  lift the invariant through the complete well-founded main loop and its
  concrete exit condition.

## Complete the generated LLL Lovasz certificate

The accumulated Lovasz invariant now covers the other concrete branch and
the entire generated main loop:

- `lllStep_swapped_array_witness` reconstructs the exact final `mu` and norm
  arrays from the successful generated swap, correction, and
  `updateMuAfterSwapLoop` executions;
- `lllStep_swapped_preserves_lovaszPrefix` uses those physical array
  equalities to prove every pair strictly before the swapped pair is
  unchanged;
- `lllStep_preserves_lovaszPrefix` combines the advancing and swapping cases;
- `concreteLLLMainLoop_preserves_lovaszPrefix` follows the same determinant-
  ranked well-founded recursion as generated C++; and
- `concreteLLLMainLoop_lovasz_all` combines the preserved cursor prefix with
  the literal `k >= matrix.size` exit theorem, yielding the Lovasz inequality
  for every adjacent pair in the returned basis.

No abstract LLL result, alternate recursion, fuel, partial definition,
`sorry`, or custom axiom is used.  In particular, the failed comparison is
not treated as a mathematical swap oracle: all three generated row swaps and
the physical post-swap mu update must succeed before the witness exists.

- C++ changes: none, so this step has no new C++ regression or B2B change
  surface.
- Next: construct the analogous generated size-reduction prefix (including
  the exact `roundQQ` residual bound), expose both final certificates through
  `lllReduce`, and prove the required short-vector span result.

## Complete the generated LLL size-reduction certificate

The second half of the physical LLL postcondition is now proved through the
same generated execution:

- `roundQQ_eq_round` unfolds the generated `2*num + den`, `2*den`, and
  `Int.fdiv` computation and identifies it with nearest-integer rounding;
- `roundQQ_residual_abs_le_half` obtains the exact half-unit residual bound;
- `sizeReduceAt_source_abs_le_half` transfers that bound to the physical
  updated `mu[k][j]` cell;
- `extraSizeReduceLoop_abs_le_half` follows the generated descending loop and
  covers every lower source column;
- the advanced and swapped branch theorems preserve `SizeReducedPrefix` on
  the actual returned arrays; and
- `concreteLLLMainLoop_size_reduced_all` follows the determinant-ranked
  well-founded main loop and combines the accumulated prefix with its literal
  exit guard to prove `|mu[i][j]| <= 1/2` for every `j < i` in the final
  generated state.

This proof does not replace source rounding with an executable oracle: the
source quotient must run successfully and is then proved equal to the
mathematical round.  There is no fuel, partial definition, semantic LLL
oracle, `sorry`, or custom axiom.

- C++ changes: none, so this step has no new C++ regression or B2B change
  surface.
- Next: retain the internal final `mu`/norm certificate at the complete
  `lllReduce` boundary and use the genuine reduced-basis inequalities to prove
  that every sufficiently short true factor-indicator vector is represented
  by the physically collected transform rows.
