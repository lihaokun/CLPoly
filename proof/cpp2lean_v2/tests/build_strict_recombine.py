"""Generate strict well-founded range loops from ``__vanhoeij_recombine``."""

from __future__ import annotations

import argparse
from pathlib import Path


OUT = Path(__file__).resolve().parents[2] / "lean/CLPoly/Generated/StrictRecombine.lean"


def generate() -> str:
    return r'''-- Auto-generated strict raw C++ control flow for recombination helpers in
-- `__vanhoeij_recombine`.
import CLPoly.Model

set_option autoImplicit false

namespace Generated.StrictRecombine

/-- Exact lowering of the range-for that maps active indices back to the
original Hensel-lifted factor array.  Invalid source indices are observable
raw faults instead of `get!` defaults. -/
def gatherActiveLoop (active : Array Int32) (lifted : Array SparsePolyZZ)
    (index : Nat) (result : Array SparsePolyZZ) :
    RawExec (Array SparsePolyZZ) :=
  if hindex : index < active.size then
    let sourceIndex := active[index]
    if hnonnegative : 0 ≤ sourceIndex then
      let sourceNat := sourceIndex.toInt64.toNat
      if hsource : sourceNat < lifted.size then
        gatherActiveLoop active lifted (index + 1)
          (result.push lifted[sourceNat])
      else .error (.outOfBounds sourceNat lifted.size)
    else .error .arithmeticDomain
  else .ok result
termination_by active.size - index
decreasing_by omega

def gatherActive (active : Array Int32) (lifted : Array SparsePolyZZ) :
    RawExec (Array SparsePolyZZ) :=
  gatherActiveLoop active lifted 0 #[]

/-- Exact lowering of the fallback range-for that moves every Zassenhaus
factor into the already accumulated van-Hoeij result. -/
def appendFallbackLoop (fallback : Array SparsePolyZZ) (index : Nat)
    (result : Array SparsePolyZZ) : RawExec (Array SparsePolyZZ) :=
  if hindex : index < fallback.size then
    appendFallbackLoop fallback (index + 1) (result.push fallback[index])
  else .ok result
termination_by fallback.size - index
decreasing_by omega

def appendFallback (fallback result : Array SparsePolyZZ) :
    RawExec (Array SparsePolyZZ) :=
  appendFallbackLoop fallback 0 result

/-- Reverse source loop removing every active entry marked consumed.  The
reverse direction is semantically relevant: erasing a larger index preserves
all still-to-be-tested smaller indices. -/
def removeConsumedLoop (consumed : Array Bool) :
    (remaining : Nat) → Array Int32 → RawExec (Array Int32)
  | 0, active => .ok active
  | remaining + 1, active =>
      let index := remaining
      if hconsumed : index < consumed.size then
        if consumed[index] then
          if hactive : index < active.size then
            removeConsumedLoop consumed remaining
              (active.eraseIdxIfInBounds index)
          else .error (.outOfBounds index active.size)
        else removeConsumedLoop consumed remaining active
      else .error (.outOfBounds index consumed.size)
termination_by remaining _ => remaining

def removeConsumed (active : Array Int32) (consumed : Array Bool) :
    RawExec (Array Int32) :=
  if _hsizes : consumed.size = active.size then
    removeConsumedLoop consumed active.size active
  else .error .assertionFailure

inductive PrecisionAction where
  | retry (target : Nat)
  | fallback
deriving DecidableEq

/-- Exact no-factor branch of the source precision schedule. -/
def nextPrecision (target initial maximum : Nat) : PrecisionAction :=
  let next := if target = 0 then initial else target * 2
  if maximum < next then .fallback else .retry next

/-- Number of remaining strict target increases before fallback.  It is
termination evidence derived from the bounded precision state, not an
execution counter passed to the source loop. -/
def precisionRank (target _initial maximum : Nat) : Nat :=
  if target = 0 then maximum + 2
  else maximum + 1 - target

theorem nextPrecision_retry_decreases (target initial maximum next : Nat)
    (hinitial : 0 < initial)
    (haction : nextPrecision target initial maximum = .retry next) :
    precisionRank next initial maximum < precisionRank target initial maximum := by
  unfold nextPrecision at haction
  split at haction <;> rename_i htarget
  · dsimp at haction
    split at haction
    next hover => contradiction
    next hfits =>
      cases PrecisionAction.retry.inj haction
      simp [precisionRank, htarget, Nat.ne_of_gt hinitial]
      omega
  · dsimp at haction
    split at haction
    next hover => contradiction
    next hfits =>
      cases PrecisionAction.retry.inj haction
      simp [precisionRank, htarget]
      omega

end Generated.StrictRecombine
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
