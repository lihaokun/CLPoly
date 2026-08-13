-- Auto-generated strict raw C++ control flow for `__select_prime`.
import CLPoly.Generated.StrictFactorZp

set_option autoImplicit false

namespace Generated.StrictSelectPrime

/-- Result of executing the concrete body for one candidate prime. -/
inductive CandidateResult (State : Type) where
  | skip (rng : State)
  | factored (fp : SparsePolyZp) (factors : Array SparsePolyZp) (rng : State)

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

private structure LoopState (State : Type) where
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
