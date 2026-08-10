"""Generate strict well-founded raw control flow for C++ ``__edf_Zp``.

The source random retry loop has no unconditional termination argument.  The
generated state therefore requires a rank on concrete RNG states and proofs
that every failed retry decreases it.  This is a termination contract for the
real RNG execution, not a fuel counter or an L2 factor oracle.
"""

from __future__ import annotations

import argparse
from pathlib import Path


OUT = (
    Path(__file__).resolve().parents[2] / "lean" / "CLPoly" / "Generated" /
    "StrictEDF.lean"
)


def generate_strict_edf() -> str:
    return r'''-- Auto-generated strict raw C++ control flow for `__edf_Zp`.
import CLPoly.Model

set_option autoImplicit false
set_option maxErrors 2000

namespace Generated.StrictEDF

/-- Source RNG boundary used by the generated EDF translation.  `State` is
the concrete C++ engine state; an implementation must provide the actual
uniform draw-and-advance operation and prove its `[0, upper)` range. -/
structure RandomEngine (State : Type) where
  nextAdvance : State → UInt64 → UInt64 × State
  nextLt : ∀ state upper, 0 < upper → (nextAdvance state upper).1 < upper

def edfMeasure (f : SparsePolyZp) : Nat :=
  if f.isEmpty then 0 else f[0]!.fst.deg + 1

def certifyRawExec {α : Type} (run : RawExec α) :
    RawExec { output : α // run = .ok output } :=
  match hrun : run with
  | .error fault => .error fault
  | .ok output => .ok ⟨output, by simpa using hrun⟩

def certifyBool (value : Bool) : { result : Bool // value = result } :=
  ⟨value, rfl⟩

def __make_zp_ir (val : Int64) (p : UInt64) : Zp :=
  Zp.ofInt val.toInt p

/-- Total lowering of the generated range-for loop in
`__upoly_subtract_one`. -/
def _loop___upoly_subtract_one_0_raw_ir (idx : Nat) (found : Bool)
    (result input : SparsePolyZp) (p : UInt64) : Bool × SparsePolyZp :=
  if hidx : idx < input.size then
    let term := input[idx]
    let next :=
      if term.fst.deg == (0 : Int64) then
        let newCoefficient := term.snd - __make_zp_ir 1 p
        if newCoefficient.val != 0 then
          (true, result.push (UMonomial.mk 0, newCoefficient))
        else
          (true, result)
      else
        (found, result.push term)
    _loop___upoly_subtract_one_0_raw_ir (idx + 1) next.1 next.2 input p
  else
    (found, result)
termination_by input.size - idx
decreasing_by omega

/-- Exact generated `__upoly_subtract_one` control flow, now total. -/
def __upoly_subtract_one_raw_ir (h : SparsePolyZp) (p : UInt64) :
    RawExec SparsePolyZp :=
  let loopResult := _loop___upoly_subtract_one_0_raw_ir 0 false #[] h p
  if !loopResult.1 then
    .ok (SparsePolyZp.normalization (loopResult.2.push
      (UMonomial.mk 0, Zp.ofInt (p - 1).toInt p)))
  else
    .ok loopResult.2

/-- Well-founded lowering of the descending degree loop in
`__upoly_random`.  Coefficients are drawn from the requested modulus range;
zero coefficients are omitted exactly as in the C++ sparse constructor. -/
def _loop___upoly_random_0_raw_ir {State : Type}
    (engine : RandomEngine State) (remaining degree : Nat)
    (result : SparsePolyZp) (p : UInt64) (rng : State) :
    SparsePolyZp × State :=
  if hremaining : remaining = 0 then
    (result, rng)
  else
    let draw := engine.nextAdvance rng p
    let coefficient := draw.1
    let result' := if coefficient != 0 then
      result.push (UMonomial.mk degree, Zp.ofInt coefficient.toInt p)
    else result
    _loop___upoly_random_0_raw_ir engine (remaining - 1) (degree - 1)
      result' p draw.2
termination_by remaining
decreasing_by omega

/-- Exact total entry for the source `__upoly_random` loop. -/
def __upoly_random_raw_ir {State : Type} (engine : RandomEngine State)
    (maxDegree : Int64) (p : UInt64) (rng : State) :
    RawExec (SparsePolyZp × State) :=
  if hpositive : 0 < maxDegree then
    let count := maxDegree.toNatClampNeg
    .ok (_loop___upoly_random_0_raw_ir engine count (count - 1) #[] p rng)
  else
    .ok (#[], rng)

structure EDFRawOps (State : Type) where
  random : Int64 → UInt64 → State → RawExec (SparsePolyZp × State)
  modPoly : SparsePolyZp → SparsePolyZp → RawExec SparsePolyZp
  squareAddMod : SparsePolyZp → SparsePolyZp → SparsePolyZp →
    RawExec SparsePolyZp
  powmod : SparsePolyZp → Int → SparsePolyZp → RawExec SparsePolyZp
  gcd : SparsePolyZp → SparsePolyZp → RawExec SparsePolyZp
  exactDiv : SparsePolyZp → SparsePolyZp → RawExec SparsePolyZp
  makeMonic : SparsePolyZp → RawExec SparsePolyZp
  EntryInvariant : SparsePolyZp → UInt64 → Prop
  splitStep : ∀ (f : SparsePolyZp) (d : UInt64)
      (g hRaw gMonic hMonic : SparsePolyZp),
    EntryInvariant f d →
    exactDiv f g = .ok hRaw →
    makeMonic g = .ok gMonic →
    makeMonic hRaw.normalization = .ok hMonic →
    EntryInvariant gMonic d ∧ EntryInvariant hMonic d ∧
      edfMeasure gMonic < edfMeasure f ∧
      edfMeasure hMonic < edfMeasure f

def traceLoop {State : Type} (ops : EDFRawOps State) (d : UInt64) (f r : SparsePolyZp)
    (i : UInt64) (g : SparsePolyZp)
    (hbudget : i.toNat ≤ d.toNat ∧ d.toNat < UInt64.size) :
    RawExec SparsePolyZp :=
  if hmore : i < d then
    match certifyRawExec (ops.squareAddMod g r f) with
    | .error fault => .error fault
    | .ok hnext =>
      traceLoop ops d f r (i + 1) hnext.val (by
        have hi : i.toNat < d.toNat := by simpa using hmore
        have hiSize : i.toNat + 1 < UInt64.size :=
          Nat.lt_of_le_of_lt (Nat.succ_le_of_lt hi) hbudget.2
        constructor
        · simpa [UInt64.toNat_add, Nat.mod_eq_of_lt hiSize] using
            Nat.succ_le_of_lt hi
        · exact hbudget.2)
  else
    .ok g
termination_by d.toNat - i.toNat
decreasing_by
  have hi : i.toNat < d.toNat := by simpa using hmore
  have hiSize : i.toNat + 1 < UInt64.size :=
    Nat.lt_of_le_of_lt (Nat.succ_le_of_lt hi) hbudget.2
  simp [UInt64.toNat_add, Nat.mod_eq_of_lt hiSize]
  omega

def candidateRun {State : Type} (ops : EDFRawOps State) (f : SparsePolyZp) (d : UInt64)
    (r : SparsePolyZp) (hbudget : 0 < d.toNat ∧ d.toNat < UInt64.size) :
    RawExec SparsePolyZp :=
  if f[0]!.2.prime == 2 then
    match certifyRawExec (ops.modPoly r f) with
    | .error fault => .error fault
    | .ok hg0 =>
      match traceLoop ops d f r 1 hg0.val ⟨by simpa using hbudget.1, hbudget.2⟩ with
      | .error fault => .error fault
      | .ok trace => ops.gcd trace f
  else
    let exponent : Int := (f[0]!.2.prime.toNat ^ d.toNat - 1) / 2
    match certifyRawExec (ops.powmod r exponent f) with
    | .error fault => .error fault
    | .ok hpow =>
      match certifyRawExec (__upoly_subtract_one_raw_ir hpow.val f[0]!.2.prime) with
      | .error fault => .error fault
      | .ok hminus => ops.gcd hminus.val f

structure EDFRetryLaw {State : Type} (ops : EDFRawOps State) where
  RetryInvariant : SparsePolyZp → UInt64 → State → Prop
  retryRank : SparsePolyZp → UInt64 → State → Nat
  retryDegreePositive : ∀ (f : SparsePolyZp) (d : UInt64) (rng : State),
    RetryInvariant f d rng → 0 < d.toNat
  retryInitial : ∀ (f : SparsePolyZp) (d : UInt64) (rng : State),
    ops.EntryInvariant f d → RetryInvariant f d rng
  emptyRetry : ∀ (f : SparsePolyZp) (d : UInt64) (rng : State)
      (r : SparsePolyZp) (rngNext : State),
    RetryInvariant f d rng →
    ops.random (get_deg f) f[0]!.2.prime rng = .ok (r, rngNext) →
    r.isEmpty = true →
    RetryInvariant f d rngNext ∧
      retryRank f d rngNext < retryRank f d rng
  failedRetry : ∀ (f : SparsePolyZp) (d : UInt64) (rng : State)
      (r : SparsePolyZp) (rngNext : State) (candidate : SparsePolyZp)
      (hbudget : 0 < d.toNat ∧ d.toNat < UInt64.size),
    RetryInvariant f d rng →
    ops.random (get_deg f) f[0]!.2.prime rng = .ok (r, rngNext) →
    r.isEmpty = false →
    candidateRun ops f d r hbudget = .ok candidate →
    ¬(get_deg candidate > 0 ∧ get_deg candidate < get_deg f) →
    RetryInvariant f d rngNext ∧
      retryRank f d rngNext < retryRank f d rng

structure RetryState {State : Type} (ops : EDFRawOps State) (law : EDFRetryLaw ops)
    (f : SparsePolyZp) (d : UInt64) where
  rng : State
  valid : law.RetryInvariant f d rng

structure SplitState {State : Type} (ops : EDFRawOps State) (f : SparsePolyZp) (d : UInt64) where
  factor : SparsePolyZp
  rngBefore : State
  rng : State
  randomPoly : SparsePolyZp
  randomRun : ops.random (get_deg f) f[0]!.2.prime rngBefore = .ok (randomPoly, rng)
  candidateRun : ∃ hbudget : 0 < d.toNat ∧ d.toNat < UInt64.size,
    Generated.StrictEDF.candidateRun ops f d randomPoly hbudget = .ok factor
  nonempty : factor.isEmpty = false
  proper : get_deg factor > 0 ∧ get_deg factor < get_deg f

def retryLoop {State : Type} (ops : EDFRawOps State) (law : EDFRetryLaw ops)
    (f : SparsePolyZp) (d : UInt64)
    (state : RetryState ops law f d) : RawExec (SplitState ops f d) :=
  match certifyRawExec (ops.random (get_deg f) f[0]!.2.prime state.rng) with
  | .error fault => .error fault
  | .ok hrandom =>
    let r := hrandom.val.1
    let rngNext := hrandom.val.2
    have hrandomRun := hrandom.property
    match certifyBool r.isEmpty with
    | ⟨true, hrempty⟩ =>
      have hretry := law.emptyRetry f d state.rng r rngNext state.valid
        hrandomRun hrempty
      retryLoop ops law f d ⟨rngNext, hretry.1⟩
    | ⟨false, hrnonempty⟩ =>
      have hbudget : 0 < d.toNat ∧ d.toNat < UInt64.size :=
        ⟨law.retryDegreePositive f d state.rng state.valid, d.toNat_lt_size⟩
      match certifyRawExec (candidateRun ops f d r hbudget) with
      | .error fault => .error fault
      | .ok hcandidate =>
        let candidate := hcandidate.val
        match certifyBool (!candidate.isEmpty && get_deg candidate > 0 &&
            get_deg candidate < get_deg f) with
        | ⟨true, hproper⟩ =>
          have hp : (candidate.isEmpty = false ∧ get_deg candidate > 0) ∧
              get_deg candidate < get_deg f := by
            simpa [Bool.and_eq_true] using hproper
          .ok ⟨candidate, state.rng, rngNext, r, hrandomRun,
            ⟨hbudget, hcandidate.property⟩, hp.1.1, ⟨hp.1.2, hp.2⟩⟩
        | ⟨false, hfailed⟩ =>
          have hretry := law.failedRetry f d state.rng r rngNext candidate
            hbudget state.valid hrandomRun hrnonempty hcandidate.property (by
              intro hproper
              have hcne : candidate ≠ #[] := by
                intro hempty
                rw [hempty] at hproper
                simp [get_deg] at hproper
              have hcempty : candidate.isEmpty = false := by
                simpa [Array.isEmpty_iff]
              have : (!candidate.isEmpty && get_deg candidate > 0 &&
                  get_deg candidate < get_deg f) = true := by
                simp [hcempty, hproper]
              rw [hfailed] at this
              contradiction)
          retryLoop ops law f d ⟨rngNext, hretry.1⟩
termination_by law.retryRank f d state.rng
decreasing_by all_goals exact hretry.2

structure EDFState {State : Type} (ops : EDFRawOps State) where
  result : Array SparsePolyZp
  f : SparsePolyZp
  d : UInt64
  rng : State
  valid : ops.EntryInvariant f d

def __edf_Zp_raw_ir_state {State : Type} (ops : EDFRawOps State) (law : EDFRetryLaw ops)
    (state : EDFState ops) :
    RawExec (Array SparsePolyZp × State) :=
  if (get_deg state.f).toUInt64 == state.d then
    match certifyRawExec (ops.makeMonic state.f) with
    | .error fault => .error fault
    | .ok hmonic => .ok (state.result.push hmonic.val, state.rng)
  else if get_deg state.f ≤ 0 then
    .ok (state.result, state.rng)
  else
    match retryLoop ops law state.f state.d
        ⟨state.rng, law.retryInitial state.f state.d state.rng state.valid⟩ with
    | .error fault => .error fault
    | .ok split =>
      match certifyRawExec (ops.exactDiv state.f split.factor) with
      | .error fault => .error fault
      | .ok hhRaw =>
        match certifyRawExec (ops.makeMonic split.factor) with
        | .error fault => .error fault
        | .ok hgMonic =>
          match certifyRawExec (ops.makeMonic hhRaw.val.normalization) with
          | .error fault => .error fault
          | .ok hhMonic =>
            have hstep := ops.splitStep state.f state.d split.factor hhRaw.val
              hgMonic.val hhMonic.val state.valid hhRaw.property
              hgMonic.property hhMonic.property
            match __edf_Zp_raw_ir_state ops law
                ⟨state.result, hgMonic.val, state.d, split.rng, hstep.1⟩ with
            | .error fault => .error fault
            | .ok leftRun =>
              __edf_Zp_raw_ir_state ops law
                ⟨leftRun.1, hhMonic.val, state.d, leftRun.2, hstep.2.1⟩
termination_by edfMeasure state.f
decreasing_by
  · exact hstep.2.2.1
  · exact hstep.2.2.2

def __edf_Zp_raw_ir {State : Type} (ops : EDFRawOps State) (law : EDFRetryLaw ops)
    (result : Array SparsePolyZp)
    (f : SparsePolyZp) (d : UInt64) (rng : State)
    (hinitial : ops.EntryInvariant f d) : RawExec (Array SparsePolyZp × State) :=
  __edf_Zp_raw_ir_state ops law ⟨result, f, d, rng, hinitial⟩

end Generated.StrictEDF
'''


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    source = generate_strict_edf()
    forbidden = ("sorry", "admit", "axiom", "partial def", "fuel", "edf ")
    if any(token in source for token in forbidden):
        raise RuntimeError("strict EDF output contains placeholder or L2 fallback")
    if args.check:
        if not args.output.exists() or args.output.read_text() != source:
            raise RuntimeError(f"{args.output} is not reproducible")
        print(f"PASS: {args.output} is reproducible and placeholder-free")
        return 0
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(source)
    print(f"generated {args.output} ({source.count(chr(10))} lines)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
