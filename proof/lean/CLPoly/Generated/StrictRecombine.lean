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

/-- Mutable source variables relevant to control flow.  Matrix contents and
candidate construction stay behind raw operation boundaries until their own
strict loops are refined. -/
structure VanHoeijState where
  active : Array Int32
  fStar : SparsePolyZZ
  result : Array SparsePolyZZ
  target : Nat

/-- Raw C++ callees for one van-Hoeij iteration.  `validateCandidates`
returns only concrete mutations and consumed bits; it cannot return an L2
factorization proposition. -/
structure VanHoeijRawOps where
  prepareCandidates : SparsePolyZZ → Array SparsePolyZZ → ZZ → Nat →
    RawExec (Array (Array Int32))
  validateCandidates : SparsePolyZZ → Array SparsePolyZZ → ZZ →
    Array (Array Int32) → Array SparsePolyZZ →
    RawExec (SparsePolyZZ × Array SparsePolyZZ × Array Bool)
  zassenhaus : SparsePolyZZ → Array SparsePolyZZ → ZZ →
    RawExec (Array SparsePolyZZ)

/-- Check one candidate's active-relative indices exactly as the source inner
loop does before constructing its trial product. -/
def candidateAvailableLoop (candidate : Array Int32) (consumed : Array Bool)
    (index : Nat) : RawExec Bool :=
  if hindex : index < candidate.size then
    let activeIndex := candidate[index]
    if hnonnegative : 0 ≤ activeIndex then
      let activeNat := activeIndex.toInt64.toNat
      if hactive : activeNat < consumed.size then
        if consumed[activeNat] then .ok false
        else candidateAvailableLoop candidate consumed (index + 1)
      else .error (.outOfBounds activeNat consumed.size)
    else .error .arithmeticDomain
  else .ok true
termination_by candidate.size - index
decreasing_by omega

def candidateAvailable (candidate : Array Int32) (consumed : Array Bool) :
    RawExec Bool := candidateAvailableLoop candidate consumed 0

/-- Concrete mutation performed after a successful trial division. -/
def markConsumedLoop (candidate : Array Int32) (index : Nat)
    (consumed : Array Bool) : RawExec (Array Bool) :=
  if hindex : index < candidate.size then
    let activeIndex := candidate[index]
    if hnonnegative : 0 ≤ activeIndex then
      let activeNat := activeIndex.toInt64.toNat
      if hactive : activeNat < consumed.size then
        markConsumedLoop candidate (index + 1) (consumed.set activeNat true)
      else .error (.outOfBounds activeNat consumed.size)
    else .error .arithmeticDomain
  else .ok consumed
termination_by candidate.size - index
decreasing_by omega

structure TrialProductRawOps where
  multiplyNormalizeMod : SparsePolyZZ → SparsePolyZZ → ZZ →
    RawExec SparsePolyZZ

/-- Source candidate-product loop: multiply by every selected active lifted
factor, normalize, and reduce coefficients after each multiplication. -/
def trialProductLoop (ops : TrialProductRawOps)
    (candidate : Array Int32) (activeLifted : Array SparsePolyZZ)
    (modulus : ZZ) (index : Nat) (product : SparsePolyZZ) :
    RawExec SparsePolyZZ :=
  if hindex : index < candidate.size then
    let activeIndex := candidate[index]
    if hnonnegative : 0 ≤ activeIndex then
      let activeNat := activeIndex.toInt64.toNat
      if hactive : activeNat < activeLifted.size then
        match ops.multiplyNormalizeMod product activeLifted[activeNat] modulus with
        | .error fault => .error fault
        | .ok product' => trialProductLoop ops candidate activeLifted modulus
            (index + 1) product'
      else .error (.outOfBounds activeNat activeLifted.size)
    else .error .arithmeticDomain
  else .ok product
termination_by candidate.size - index
decreasing_by omega

/-- Inner source multiplication loop for one left sparse term. -/
def multiplyRowLoop (left : UMonomial × Int) (right : SparsePolyZZ)
    (rightIndex : Nat) (terms : SparsePolyZZ) : SparsePolyZZ :=
  if hright : rightIndex < right.size then
    let rightTerm := right[rightIndex]
    multiplyRowLoop left right (rightIndex + 1)
      (terms.push (⟨left.1.deg + rightTerm.1.deg⟩,
        left.2 * rightTerm.2))
  else terms
termination_by right.size - rightIndex
decreasing_by omega

/-- Exact double range-for producing every monomial product before the C++
normalization pass merges equal degrees. -/
def multiplyTermsLoop (left right : SparsePolyZZ) (leftIndex : Nat)
    (terms : SparsePolyZZ) : SparsePolyZZ :=
  if hleft : leftIndex < left.size then
    multiplyTermsLoop left right (leftIndex + 1)
      (multiplyRowLoop left[leftIndex] right 0 terms)
  else terms
termination_by left.size - leftIndex
decreasing_by omega

def multiplyNormalizeRaw (left right : SparsePolyZZ) : RawExec SparsePolyZZ :=
  .ok (SparsePolyZZ.normalization (multiplyTermsLoop left right 0 #[]))

/-- Erased termination certificate for the successful-extraction branch.
It refers only to concrete successful raw executions and the consumed bits
they returned. -/
structure VanHoeijTermination (ops : VanHoeijRawOps) where
  extraction_decreases : ∀ (modulus : ZZ) (state : VanHoeijState)
      activeLifted candidates fStar' result'
      consumed active',
    ops.validateCandidates state.fStar activeLifted modulus candidates state.result =
      .ok (fStar', result', consumed) →
    (∃ index, ∃ hindex : index < consumed.size, consumed[index] = true) →
    removeConsumed state.active consumed = .ok active' →
    active'.size < state.active.size

def removeConsumedDecreasing (ops : VanHoeijRawOps)
    (termination : VanHoeijTermination ops) (modulus : ZZ)
    (state : VanHoeijState) (activeLifted : Array SparsePolyZZ)
    (candidates : Array (Array Int32)) (fStar' : SparsePolyZZ)
    (result' : Array SparsePolyZZ) (consumed : Array Bool)
    (hvalidate : ops.validateCandidates state.fStar activeLifted modulus
      candidates state.result = .ok (fStar', result', consumed))
    (hfound : ∃ index, ∃ hindex : index < consumed.size,
      consumed[index] = true) :
    RawExec { active' : Array Int32 // active'.size < state.active.size } :=
  match hremove : removeConsumed state.active consumed with
  | .error fault => .error fault
  | .ok active' => .ok ⟨active', termination.extraction_decreases modulus
      state activeLifted candidates fStar' result' consumed active'
      hvalidate hfound hremove⟩

/-- Source-shaped van-Hoeij main loop.  Successful extraction decreases the
active set; an unsuccessful round strictly advances bounded precision; an
overflowing precision target executes the actual Zassenhaus fallback. -/
def vanHoeijLoop (ops : VanHoeijRawOps) (termination : VanHoeijTermination ops)
    (lifted : Array SparsePolyZZ)
    (modulus : ZZ) (initial maximum : Nat) (hinitial : 0 < initial) :
    VanHoeijState → RawExec (SparsePolyZZ × Array SparsePolyZZ)
  | state =>
    if hdone : state.active.size ≤ 1 then .ok (state.fStar, state.result)
    else
      match gatherActive state.active lifted with
      | .error fault => .error fault
      | .ok activeLifted =>
        match ops.prepareCandidates state.fStar activeLifted modulus state.target with
        | .error fault => .error fault
        | .ok candidates =>
          match hvalidate : ops.validateCandidates state.fStar activeLifted modulus candidates
              state.result with
          | .error fault => .error fault
          | .ok (fStar', result', consumed) =>
            if hfound : ∃ index, ∃ hindex : index < consumed.size,
                consumed[index] = true then
              match removeConsumedDecreasing ops termination modulus state
                  activeLifted candidates fStar' result' consumed hvalidate hfound with
              | .error fault => .error fault
              | .ok activeNext =>
                vanHoeijLoop ops termination lifted modulus initial maximum hinitial
                  { active := activeNext.1, fStar := fStar', result := result',
                    target := 0 }
            else
              match hprecision : nextPrecision state.target initial maximum with
              | .retry target' =>
                vanHoeijLoop ops termination lifted modulus initial maximum hinitial
                  { state with target := target' }
              | .fallback =>
                match ops.zassenhaus state.fStar activeLifted modulus with
                | .error fault => .error fault
                | .ok fallback =>
                  match appendFallback fallback state.result with
                  | .error fault => .error fault
                  | .ok output => .ok (#[], output)
termination_by state =>
  (state.active.size, precisionRank state.target initial maximum)
decreasing_by
  · exact Prod.Lex.left _ _ activeNext.2
  · exact Prod.Lex.right _
      (nextPrecision_retry_decreases state.target initial maximum target'
        hinitial hprecision)

end Generated.StrictRecombine
