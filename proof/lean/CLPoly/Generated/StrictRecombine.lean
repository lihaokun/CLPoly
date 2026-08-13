-- Auto-generated strict raw C++ control flow for recombination helpers in
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

end Generated.StrictRecombine
