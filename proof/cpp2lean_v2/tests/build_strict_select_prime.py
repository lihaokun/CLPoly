"""Generate the well-founded C++ ``__select_prime`` control flow."""

from __future__ import annotations

import argparse
from pathlib import Path


OUT = Path(__file__).resolve().parents[2] / "lean/CLPoly/Generated/StrictSelectPrime.lean"


def generate() -> str:
    return r'''-- Auto-generated strict raw C++ control flow for `__select_prime`.
import CLPoly.Generated.StrictFactorZp

set_option autoImplicit false

namespace Generated.StrictSelectPrime

/-- Result of executing the concrete body for one candidate prime. -/
inductive CandidateResult (State : Type) where
  | skip (rng : State)
  | factored (fp : SparsePolyZp) (factors : Array SparsePolyZp) (rng : State)

/-- Concrete callees in one C++ candidate attempt. -/
structure CandidateRawOps (State : Type) where
  lcMod : ZZ → UInt64 → RawExec ZZ
  polynomialMod : SparsePolyZZ → UInt64 → RawExec SparsePolyZp
  derivative : SparsePolyZp → RawExec SparsePolyZp
  gcd : SparsePolyZp → SparsePolyZp → RawExec SparsePolyZp
  makeMonic : SparsePolyZp → RawExec (Zp × SparsePolyZp)
  ddf : SparsePolyZp → RawExec (Array (SparsePolyZp × UInt64))
  edf : Array SparsePolyZp → SparsePolyZp → UInt64 → State →
    RawExec (Array SparsePolyZp × State)

/-- Append one concrete EDF output array, preserving the source range-for. -/
def appendFactorsLoop (source : Array SparsePolyZp) (index : Nat)
    (result : Array SparsePolyZp) : Array SparsePolyZp :=
  if h : index < source.size then
    appendFactorsLoop source (index + 1) (result.push source[index])
  else result
termination_by source.size - index
decreasing_by omega

/-- Execute EDF for every concrete DDF component. -/
def factorComponentsLoop {State : Type} (ops : CandidateRawOps State)
    (components : Array (SparsePolyZp × UInt64)) (index : Nat)
    (rng : State) (result : Array SparsePolyZp) :
    RawExec (Array SparsePolyZp × State) :=
  if h : index < components.size then
    let component := components[index]
    match ops.edf #[] component.1 component.2 rng with
    | .error fault => .error fault
    | .ok (edfOutput, rng') =>
      factorComponentsLoop ops components (index + 1) rng'
        (appendFactorsLoop edfOutput 0 result)
  else .ok (result, rng)
termination_by components.size - index
decreasing_by omega

/-- Exact candidate body: coefficient reduction, polynomial reduction,
degree preservation, derivative/GCD squarefree test, make-monic, DDF and EDF. -/
def tryCandidateRaw {State : Type} (ops : CandidateRawOps State)
    (f : SparsePolyZZ) (degF : Int64) (lcF : ZZ) (p : UInt64)
    (rng : State) : RawExec (CandidateResult State) :=
  match ops.lcMod lcF p with
  | .error fault => .error fault
  | .ok lcMod =>
    if !ZZ.toBool lcMod then .ok (.skip rng)
    else
      match ops.polynomialMod f p with
      | .error fault => .error fault
      | .ok fp =>
        if fp.isEmpty || get_deg fp != degF then .ok (.skip rng)
        else
          match ops.derivative fp with
          | .error fault => .error fault
          | .ok fpDerivative =>
            if fpDerivative.isEmpty then .ok (.skip rng)
            else
              match ops.gcd fp fpDerivative with
              | .error fault => .error fault
              | .ok gcd =>
                if get_deg gcd > 0 then .ok (.skip rng)
                else
                  match ops.makeMonic fp with
                  | .error fault => .error fault
                  | .ok (_, monic) =>
                    match ops.ddf monic with
                    | .error fault => .error fault
                    | .ok components =>
                      match factorComponentsLoop ops components 0 rng #[] with
                      | .error fault => .error fault
                      | .ok (factors, rng') => .ok (.factored monic factors rng')

/-- Raw C++ boundaries used by prime enumeration.  `tryCandidate` is refined
separately through modular reduction, GCD, make-monic, strict DDF and strict
EDF; it is not permitted to return an L2 proposition or witness. -/
structure SelectPrimeRawOps (State : Type) where
  nextPrime : Bool → UInt64 → UInt64
  tryCandidate : SparsePolyZZ → Int64 → ZZ → UInt64 → State →
    RawExec (CandidateResult State)

/-- Termination certificate for the source prime iterator.  This is a
well-founded relation on successive machine primes, not a fuel counter. -/
structure PrimeEnumerationTermination {State : Type}
    (ops : SelectPrimeRawOps State) where
  rank : Bool → UInt64 → Nat
  next_decreases : ∀ useLargePrime p,
    rank useLargePrime (ops.nextPrime useLargePrime p) <
      rank useLargePrime p

/-- Concrete mutable variables carried by the generated prime-search loop.
This is public so the refinement proof can state its prime and best-candidate
invariants without replacing the generated control flow. -/
structure LoopState (State : Type) where
  tried : Nat
  p : UInt64
  rng : State
  bestCount : Nat
  best : PrimeSelectionResult

/-- Strict lowering of the C++ candidate loop.  The lexicographic measure
first decreases when a valid candidate is counted and otherwise follows the
well-founded machine-prime iterator. -/
def selectPrimeLoop {State : Type} (ops : SelectPrimeRawOps State)
    (termination : PrimeEnumerationTermination ops) (f : SparsePolyZZ)
    (degF : Int64) (lcF : ZZ) (useLargePrime : Bool) (maxTries : Nat)
    (state : LoopState State) : RawExec PrimeSelectionResult :=
  if hdone : maxTries ≤ state.tried then
    .ok { state.best with irreducible := false }
  else
    match ops.tryCandidate f degF lcF state.p state.rng with
    | .error fault => .error fault
    | .ok (.skip rng') =>
      selectPrimeLoop ops termination f degF lcF useLargePrime maxTries
        { state with p := ops.nextPrime useLargePrime state.p, rng := rng' }
    | .ok (.factored fp factors rng') =>
      if factors.size ≤ 1 then
        let factors' := if factors.isEmpty then factors.push fp else factors
        .ok { prime := state.p, factors := factors', irreducible := true }
      else
        let best := if factors.size < state.bestCount then
          { state.best with prime := state.p, factors := factors }
        else state.best
        let bestCount := Nat.min state.bestCount factors.size
        selectPrimeLoop ops termination f degF lcF useLargePrime maxTries
          { tried := state.tried + 1
            p := ops.nextPrime useLargePrime state.p
            rng := rng'
            bestCount := bestCount
            best := best }
termination_by (maxTries - state.tried,
  termination.rank useLargePrime state.p)
decreasing_by
  · exact Prod.Lex.right _ (termination.next_decreases _ _)
  · exact Prod.Lex.left _ _ (by omega)

/-- Strict generated entry for original C++ `__select_prime`. -/
def __select_prime_raw_ir {State : Type} (ops : SelectPrimeRawOps State)
    (termination : PrimeEnumerationTermination ops) (initialRng : State)
    (useLargePrime : Bool) (f : SparsePolyZZ) : RawExec PrimeSelectionResult :=
  if f.isEmpty || get_deg f < 2 then .error .assertionFailure
  else
    let initialPrime : UInt64 :=
      if useLargePrime then (18446744073709551615 : UInt64) - 58 else 2
    selectPrimeLoop ops termination f (get_deg f) (SparsePolyZZ.front! f).2
      useLargePrime 3
      { tried := 0
        p := initialPrime
        rng := initialRng
        bestCount := 18446744073709551615
        best := default }

end Generated.StrictSelectPrime
'''


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    source = generate()
    if args.check:
        if not OUT.exists() or OUT.read_text() != source:
            raise SystemExit(f"generated output is stale: {OUT}")
    else:
        OUT.write_text(source)


if __name__ == "__main__":
    main()
