# Strict FactorZZ outer controller

## What changed

The source-shaped generated `__factor_squarefree_primitive_ZZ_raw_ir` entry is
now composed with the completed `__lll_factorize_raw_ir` refinement.

The literal early return is proved rather than assumed.  A successful real
prime-selection execution supplies a nonempty modular factor list.  If that
list has at most one member, its unique member is irreducible and associated
to the complete reduction of the source.  The selected `GoodPrime` preserves
the leading coefficient, so Gauss reduction proves that the primitive integer
source is irreducible and that the physical return `#[f]` satisfies
`FactorZZCorrect`.

The multi-factor branch invokes the already-proved literal Hensel/van-Hoeij/
Zassenhaus controller.  The selected result and Hensel outputs are constrained
only after their raw calls return; no premise selects the final result.

The select-prime field is also instantiated with the actual generated prime
enumerator and concrete polynomial-mod/GCD/DDF/EDF callback.  Its RNG state is
the local entry state, matching the source-local generator boundary.

The former `native_decide` proof of the fixed initial prime was removed after
the axiom audit exposed Lean's native-code trust axiom.  The intermediate
select-prime contract now takes primality as an explicit mathematical premise.
This keeps the current composition kernel-auditable, but the premise is not
accepted as part of the final public FactorZZ theorem: it still has to be
discharged by a compact kernel-checked primality certificate.

## Remaining work

The downstream Hensel operation and the representation facts retained from
prime selection still need concrete adapters.  During this wiring audit, the
older proof infrastructure was found to retain artificial bounds
`2*p <= 2^64` and `p*p <= 2^64` even though `Model.lean` now computes `Zp`
addition and multiplication through unbounded `Nat` intermediates before
reduction.  Those proof-only bounds exclude the C++ benchmark large-prime
mode and therefore must be removed, not exported as final theorem premises.

## Verification

- `lake env lean CLPoly/Refinement/FactorZZ.lean`: passed.
- `lake build CLPoly.Refinement.SelectPrime`: passed.
- `lake env lean CLPoly/Refinement/Generated.lean`: passed.
- `lake env lean CLPoly/Refinement/AxiomAudit.lean`: passed; the new outer
  contracts depend only on `propext`, `Classical.choice`, and `Quot.sound`.
- strict search: no `native_decide`, `partial def`, `sorry`, or fuel appears in
  the changed proof path.
- `git diff --check`: passed.
- Production C++ changes: none; no new native/B2B change surface.

## Files

- `proof/lean/CLPoly/Refinement/FactorZZ.lean`
- `proof/lean/CLPoly/Refinement/SelectPrime.lean`
- `proof/lean/CLPoly/Refinement/Generated.lean`
- `proof/lean/CLPoly/Refinement/AxiomAudit.lean`
- `proof/docs/devlog/2026-08-17-strict-factor-zz-outer-controller.md`

## Progress and estimate

The full outer branch structure and real select-prime execution are connected.
Overall goal progress is approximately 99.87%; estimated remaining time is
1.5--3 full working days, dominated by concrete Hensel readiness, removal of
the stale word-product bounds, generated public export, and final gates.

## Concrete generated Hensel adapter

The outer controller no longer needs an abstract `henselLift` implementation.
`concreteHenselLift` selects the dense arithmetic object and multiplication
workspace at the prime physically returned by `__select_prime`, then executes
the generated `__hensel_lift_upoly_raw_ir` entry with its strict tree builder,
well-founded quadratic lift loop, extraction, and normalization.

`concreteHenselLift_success` identifies a successful returned value with the
output of the existing full Hensel refinement theorem.  The composed
`concreteSelectHensel___factor_squarefree_primitive_ZZ_raw_ir_refines` theorem
therefore has no abstract select-prime or Hensel function and no caller-chosen
Hensel result.  It still retains the universally quantified Hensel execution
readiness invariant and the stale prime-square bound; both must be derived or
removed before the final public theorem is generated.

Verification:

- `lake env lean CLPoly/Refinement/FactorZZ.lean`: passed.
- Production C++ changes: none; no new native/B2B change surface.

Progress is approximately 99.89%; the remaining estimate stays 1.5--3 full
working days because removing the obsolete word-product bounds is now the
critical path rather than merely connecting the generated Hensel call.
