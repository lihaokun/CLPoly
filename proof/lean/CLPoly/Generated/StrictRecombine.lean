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

/-- The source `iota` constructing the first lexicographic subset
`[0, ..., count - 1]`. -/
def initialCombinationLoop (count index : Nat) (result : Array Nat) :
    Array Nat :=
  if hindex : index < count then
    initialCombinationLoop count (index + 1) (result.push index)
  else result
termination_by count - index
decreasing_by omega

def initialCombination (count : Nat) : Array Nat :=
  initialCombinationLoop count 0 #[]

/-- Exact right-to-left search in C++ `next_combination`.  `inspected`
counts positions already rejected at the right edge, avoiding signed-loop
indices while preserving the source comparisons. -/
def nextCombinationPivot (indices : Array Nat) (upper inspected : Nat) :
    Option Nat :=
  if hinspected : inspected < indices.size then
    let position := indices.size - 1 - inspected
    if indices[position] = upper - indices.size + position then
      nextCombinationPivot indices upper (inspected + 1)
    else some position
  else none
termination_by indices.size - inspected
decreasing_by omega

/-- Source suffix reset `idx[j] = idx[j-1] + 1` after incrementing the
pivot. -/
def resetCombinationSuffix (indices : Array Nat) (pivot offset : Nat) :
    Array Nat :=
  let position := pivot + 1 + offset
  if hposition : position < indices.size then
    resetCombinationSuffix
      (indices.set position (indices[pivot]'(by omega) + 1)) pivot (offset + 1)
  else indices
termination_by indices.size - (pivot + 1 + offset)
decreasing_by simp only [Array.size_set]; omega

/-- Fuel-free lowering of the complete C++ `next_combination` lambda. -/
def nextCombination (indices : Array Nat) (upper : Nat) : Bool × Array Nat :=
  if hfits : indices.size ≤ upper then
    match hpivot : nextCombinationPivot indices upper 0 with
    | none => (false, indices)
    | some pivot =>
        if hpivotBounds : pivot < indices.size then
          let incremented := indices.set pivot (indices[pivot] + 1)
          (true, resetCombinationSuffix incremented pivot 0)
        else (false, indices)
  else (false, indices)

/-- Exact C++ `__upoly_const_term`: canonical sparse inputs store a constant
term at the back when one exists. -/
def constantTerm (input : SparsePolyZZ) : ZZ :=
  if hempty : input.isEmpty then 0
  else
    let term := input[input.size - 1]'(by
      have : input.size ≠ 0 := by simpa [Array.isEmpty] using hempty
      omega)
    if term.1.deg = 0 then term.2 else 0

/-- Leading-coefficient pruning product over one active-relative subset. -/
def selectedLeadingProductLoop (candidate : Array Nat)
    (activeLifted : Array SparsePolyZZ) (index : Nat) (acc : ZZ) :
    RawExec ZZ :=
  if hindex : index < candidate.size then
    let activeIndex := candidate[index]
    if hactive : activeIndex < activeLifted.size then
      let factor := activeLifted[activeIndex]
      if hempty : factor.isEmpty then .error .assertionFailure
      else
        selectedLeadingProductLoop candidate activeLifted (index + 1)
          (acc * (factor[0]'(by
            have : factor.size ≠ 0 := by simpa [Array.isEmpty] using hempty
            omega)).2)
    else .error (.outOfBounds activeIndex activeLifted.size)
  else .ok acc
termination_by candidate.size - index
decreasing_by omega

/-- Constant-coefficient pruning product over the same concrete subset. -/
def selectedConstantProductLoop (candidate : Array Nat)
    (activeLifted : Array SparsePolyZZ) (index : Nat) (acc : ZZ) :
    RawExec ZZ :=
  if hindex : index < candidate.size then
    let activeIndex := candidate[index]
    if hactive : activeIndex < activeLifted.size then
      selectedConstantProductLoop candidate activeLifted (index + 1)
        (acc * constantTerm activeLifted[activeIndex])
    else .error (.outOfBounds activeIndex activeLifted.size)
  else .ok acc
termination_by candidate.size - index
decreasing_by omega

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

/-- Exact coefficient reduction performed after every Zassenhaus subset
product multiplication. -/
def modCoeffLoop (input : SparsePolyZZ) (modulus : ZZ) (index : Nat)
    (result : SparsePolyZZ) : RawExec SparsePolyZZ :=
  if hindex : index < input.size then
    if hmodulus : modulus ≠ 0 then
      let coefficient := input[index].2 % modulus
      modCoeffLoop input modulus (index + 1)
        (if coefficient = 0 then result
         else result.push (input[index].1, coefficient))
    else .error .arithmeticDomain
  else .ok result
termination_by input.size - index
decreasing_by omega

def multiplyNormalizeModRaw (left right : SparsePolyZZ) (modulus : ZZ) :
    RawExec SparsePolyZZ :=
  match multiplyNormalizeRaw left right with
  | .error fault => .error fault
  | .ok product => modCoeffLoop product modulus 0 #[]

/-- This operation record is deliberately data-free.  Earlier revisions
stored a semantic multiplication callback here; the strict boundary now
always executes `multiplyNormalizeModRaw`. -/
structure TrialProductRawOps where
  marker : Unit := ()

/-- Source candidate-product loop: multiply by every selected active lifted
factor, normalize, and reduce coefficients after each multiplication. -/
def trialProductLoop (_ops : TrialProductRawOps)
    (candidate : Array Int32) (activeLifted : Array SparsePolyZZ)
    (modulus : ZZ) (index : Nat) (product : SparsePolyZZ) :
    RawExec SparsePolyZZ :=
  if hindex : index < candidate.size then
    let activeIndex := candidate[index]
    if hnonnegative : 0 ≤ activeIndex then
      let activeNat := activeIndex.toInt64.toNat
      if hactive : activeNat < activeLifted.size then
        match multiplyNormalizeModRaw product activeLifted[activeNat] modulus with
        | .error fault => .error fault
        | .ok product' => trialProductLoop _ops candidate activeLifted modulus
            (index + 1) product'
      else .error (.outOfBounds activeNat activeLifted.size)
    else .error .arithmeticDomain
  else .ok product
termination_by candidate.size - index
decreasing_by omega

/-- Exact range-for in C++ `__upoly_symmetric_mod`: symmetrically reduce each
coefficient and omit concrete zero coefficients. -/
def symmetricModLoop (input : SparsePolyZZ) (modulus : ZZ) (index : Nat)
    (result : SparsePolyZZ) : RawExec SparsePolyZZ :=
  if hindex : index < input.size then
    let term := input[index]
    let coefficient := ZZ.symmetricMod term.2 modulus
    symmetricModLoop input modulus (index + 1)
      (if coefficient = 0 then result else result.push (term.1, coefficient))
  else .ok result
termination_by input.size - index
decreasing_by omega

def symmetricModRaw (input : SparsePolyZZ) (modulus : ZZ) :
    RawExec SparsePolyZZ :=
  if hmodulus : 0 < modulus then symmetricModLoop input modulus 0 #[]
  else .error .arithmeticDomain

def contentLoop (input : SparsePolyZZ) (index acc : Nat) : Nat :=
  if hindex : index < input.size then
    contentLoop input (index + 1) (Nat.gcd acc input[index].2.natAbs)
  else acc
termination_by input.size - index
decreasing_by omega

def primitiveDivideLoop (input : SparsePolyZZ) (divisor : ZZ) (index : Nat)
    (result : SparsePolyZZ) : RawExec SparsePolyZZ :=
  if hindex : index < input.size then
    if hdivisor : divisor ≠ 0 then
      if hdivides : divisor ∣ input[index].2 then
        primitiveDivideLoop input divisor (index + 1)
          (result.push (input[index].1, input[index].2 / divisor))
      else .error .arithmeticDomain
    else .error .arithmeticDomain
  else .ok result
termination_by input.size - index
decreasing_by omega

/-- Exact C++ `__upoly_primitive`: content gcd, leading-sign adjustment,
and coefficient-wise exact division. -/
def primitiveRaw (input : SparsePolyZZ) : RawExec (ZZ × SparsePolyZZ) :=
  if hempty : input.isEmpty then .ok (1, input)
  else
    let content : ZZ := contentLoop input 0 0
    let divisor := if (input[0]'(by
      have : input.size ≠ 0 := by simpa [Array.isEmpty] using hempty
      omega)).2 < 0
      then -content else content
    match primitiveDivideLoop input divisor 0 #[] with
    | .error fault => .error fault
    | .ok primitive => .ok (divisor, primitive)

def subtractScaledTermsLoop (divisor : SparsePolyZZ) (scale : ZZ)
    (degreeShift index : Nat) (terms : SparsePolyZZ) : SparsePolyZZ :=
  if hindex : index < divisor.size then
    subtractScaledTermsLoop divisor scale degreeShift (index + 1)
      (terms.push
        (⟨divisor[index].1.deg + degreeShift⟩, -(scale * divisor[index].2)))
  else terms
termination_by divisor.size - index
decreasing_by omega

def subtractScaledNormalize (remainder divisor : SparsePolyZZ) (scale : ZZ)
    (degreeShift : Nat) : SparsePolyZZ :=
  SparsePolyZZ.normalization
    (subtractScaledTermsLoop divisor scale degreeShift 0 remainder)

def divisionRank (remainder : SparsePolyZZ) : Nat :=
  if h : 0 < remainder.size then remainder[0].1.deg + 1 else 0

/-- Checked exact sparse long division used by the C++ trial-division call.
The degree guard is a runtime validation of the mathematical decrease, not a
fuel counter: a malformed/noncanonical intermediate produces a raw fault. -/
def exactDivmodLoop (divisor : SparsePolyZZ) :
    SparsePolyZZ → SparsePolyZZ → RawExec (SparsePolyZZ × SparsePolyZZ)
  | remainder, quotient =>
      if hremainder : 0 < remainder.size then
        if hdivisor : 0 < divisor.size then
          let remainderLead := remainder[0]
          let divisorLead := divisor[0]
          if hdegree : divisorLead.1.deg ≤ remainderLead.1.deg then
            if hnonzero : divisorLead.2 ≠ 0 then
              if hdivides : divisorLead.2 ∣ remainderLead.2 then
                let degreeShift := remainderLead.1.deg - divisorLead.1.deg
                let scale := remainderLead.2 / divisorLead.2
                let remainder' := subtractScaledNormalize remainder divisor scale degreeShift
                let quotient' := quotient.push (⟨degreeShift⟩, scale)
                if hdecrease : divisionRank remainder' < divisionRank remainder then
                  exactDivmodLoop divisor remainder' quotient'
                else .error .arithmeticDomain
              else .ok (quotient, remainder)
            else .error .arithmeticDomain
          else .ok (quotient, remainder)
        else .error .arithmeticDomain
      else .ok (quotient, remainder)
termination_by remainder quotient => divisionRank remainder
decreasing_by exact hdecrease

def exactDivmodRaw (dividend divisor : SparsePolyZZ) :
    RawExec (SparsePolyZZ × SparsePolyZZ) :=
  exactDivmodLoop divisor dividend #[]

/-- Checked lowering of the source `vector<int>` combination indices. -/
def combinationToInt32Loop (indices : Array Nat) (index : Nat)
    (result : Array Int32) : RawExec (Array Int32) :=
  if hindex : index < indices.size then
    if hfits : indices[index] < 2 ^ 31 then
      combinationToInt32Loop indices (index + 1)
        (result.push indices[index].toUInt32.toInt32)
    else .error .arithmeticDomain
  else .ok result
termination_by indices.size - index
decreasing_by omega

def combinationToInt32 (indices : Array Nat) : RawExec (Array Int32) :=
  combinationToInt32Loop indices 0 #[]

/-- Reverse erase of the selected active positions after one successful
Zassenhaus extraction. -/
def removeCombinationLoop (candidate : Array Nat) :
    (remaining : Nat) → Array SparsePolyZZ → RawExec (Array SparsePolyZZ)
  | 0, active => .ok active
  | remaining + 1, active =>
      let candidateIndex := remaining
      if hcand : candidateIndex < candidate.size then
        let activeIndex := candidate[candidateIndex]
        if hactive : activeIndex < active.size then
          removeCombinationLoop candidate remaining
            (active.eraseIdxIfInBounds activeIndex)
        else .error (.outOfBounds activeIndex active.size)
      else .error (.outOfBounds candidateIndex candidate.size)
termination_by remaining _ => remaining

def removeCombination (candidate : Array Nat)
    (active : Array SparsePolyZZ) : RawExec (Array SparsePolyZZ) :=
  removeCombinationLoop candidate candidate.size active

inductive ZassenhausAttemptResult where
  | rejected
  | extracted (factor quotient : SparsePolyZZ)

inductive ZassenhausScanResult where
  | exhausted
  | extracted (factor quotient : SparsePolyZZ) (candidate : Array Nat)

/-- One complete C++ Zassenhaus candidate attempt, including both scalar
pruning tests and the actual sparse exact division. -/
def zassenhausAttempt (fStar : SparsePolyZZ)
    (activeLifted : Array SparsePolyZZ) (modulus : ZZ)
    (candidate : Array Nat) : RawExec ZassenhausAttemptResult :=
  if hfstar : 0 < fStar.size then
    let leading := fStar[0].2
    match selectedLeadingProductLoop candidate activeLifted 0 leading with
    | .error fault => .error fault
    | .ok leadingProduct =>
      let recoveredLeading := ZZ.symmetricMod leadingProduct modulus
      let leadingSquare := leading * leading
      if recoveredLeading ≠ 0 ∧ ZZ.fdiv_r 0 leadingSquare recoveredLeading ≠ 0 then
        .ok .rejected
      else
        match selectedConstantProductLoop candidate activeLifted 0 leading with
        | .error fault => .error fault
        | .ok constantProduct =>
          let recoveredConstant := ZZ.symmetricMod constantProduct modulus
          let targetConstant := leading * constantTerm fStar
          if recoveredConstant ≠ 0 ∧
              ZZ.fdiv_r 0 targetConstant recoveredConstant ≠ 0 then
            .ok .rejected
          else
            match combinationToInt32 candidate with
            | .error fault => .error fault
            | .ok candidate32 =>
              match trialProductLoop ⟨()⟩ candidate32 activeLifted modulus 0
                  #[(⟨0⟩, leading)] with
              | .error fault => .error fault
              | .ok product =>
                match symmetricModRaw product modulus with
                | .error fault => .error fault
                | .ok symmetric =>
                  match primitiveRaw symmetric with
                  | .error fault => .error fault
                  | .ok (_, factor) =>
                    match exactDivmodRaw fStar factor with
                    | .error fault => .error fault
                    | .ok (quotient, remainder) =>
                      if remainder.isEmpty then
                        match primitiveRaw quotient with
                        | .error fault => .error fault
                        | .ok (_, quotientPrimitive) =>
                          .ok (.extracted factor quotientPrimitive)
                      else .ok .rejected
  else .error .assertionFailure

/-- Erased well-founded metric for the concrete lexicographic
`next_combination` execution.  The refinement layer constructs this
certificate from the generated combination arithmetic. -/
structure CombinationTermination (upper : Nat) where
  rank : Array Nat → Nat
  next_decreases : ∀ current next,
    nextCombination current upper = (true, next) → rank next < rank current

/-- Exact source do-while over all fixed-size combinations. -/
def scanZassenhausCombinations {upper : Nat}
    (termination : CombinationTermination upper)
    (fStar : SparsePolyZZ) (activeLifted : Array SparsePolyZZ)
    (modulus : ZZ) : Array Nat → RawExec ZassenhausScanResult
  | candidate =>
      match zassenhausAttempt fStar activeLifted modulus candidate with
      | .error fault => .error fault
      | .ok (.extracted factor quotient) =>
          .ok (.extracted factor quotient candidate)
      | .ok .rejected =>
          match hnext : nextCombination candidate upper with
          | (false, _) => .ok .exhausted
          | (true, next) =>
              scanZassenhausCombinations termination fStar activeLifted
                modulus next
termination_by candidate => termination.rank candidate
decreasing_by exact termination.next_decreases candidate next hnext

/-- Concrete C++ callees used inside candidate validation.  Each field returns
only computed polynomial data; no field may return a semantic proposition or
choose an L2 factorization witness. -/
structure CandidateValidationRawOps where
  product : TrialProductRawOps

/-- Exact source `for (auto& cand : candidates)` validation loop: reject empty,
trivial, or already-consumed candidates; build the modular product; recover a
primitive integer polynomial; run actual trial division; and mutate `fStar`,
`result`, and `consumed` only when the concrete remainder is empty. -/
def validateCandidatesLoop (ops : CandidateValidationRawOps)
    (candidates : Array (Array Int32)) (candidateIndex : Nat)
    (activeLifted : Array SparsePolyZZ) (modulus : ZZ) :
    SparsePolyZZ → Array SparsePolyZZ → Array Bool → Nat →
      RawExec (SparsePolyZZ × Array SparsePolyZZ × Array Bool)
  | fStar, result, consumed, remaining =>
      if hcandidates : candidateIndex < candidates.size then
        let candidate := candidates[candidateIndex]
        if hempty : candidate.isEmpty then
          validateCandidatesLoop ops candidates (candidateIndex + 1)
            activeLifted modulus fStar result consumed remaining
        else if htrivial : remaining ≤ candidate.size then
          validateCandidatesLoop ops candidates (candidateIndex + 1)
            activeLifted modulus fStar result consumed remaining
        else
          match candidateAvailable candidate consumed with
          | .error fault => .error fault
          | .ok false => validateCandidatesLoop ops candidates
              (candidateIndex + 1) activeLifted modulus fStar result consumed remaining
          | .ok true =>
            if hfstar : 0 < fStar.size then
              let initial : SparsePolyZZ := #[(⟨0⟩, fStar[0].2)]
              match trialProductLoop ops.product candidate activeLifted modulus 0 initial with
              | .error fault => .error fault
              | .ok product =>
                match symmetricModRaw product modulus with
                | .error fault => .error fault
                | .ok symmetric =>
                  match primitiveRaw symmetric with
                  | .error fault => .error fault
                  | .ok (_, factor) =>
                    match exactDivmodRaw fStar factor with
                    | .error fault => .error fault
                    | .ok (quotient, remainder) =>
                      if hremainder : remainder.isEmpty then
                        match primitiveRaw quotient with
                        | .error fault => .error fault
                        | .ok (_, quotientPrimitive) =>
                          match markConsumedLoop candidate 0 consumed with
                          | .error fault => .error fault
                          | .ok consumed' =>
                            validateCandidatesLoop ops candidates
                              (candidateIndex + 1) activeLifted modulus
                              quotientPrimitive (result.push factor) consumed'
                              (remaining - candidate.size)
                      else validateCandidatesLoop ops candidates
                        (candidateIndex + 1) activeLifted modulus
                        fStar result consumed remaining
            else .error (.outOfBounds 0 fStar.size)
      else .ok (fStar, result, consumed)
termination_by fStar result consumed remaining => candidates.size - candidateIndex
decreasing_by all_goals omega

def validateCandidates (ops : CandidateValidationRawOps)
    (fStar : SparsePolyZZ) (activeLifted : Array SparsePolyZZ) (modulus : ZZ)
    (candidates : Array (Array Int32)) (result : Array SparsePolyZZ) :
    RawExec (SparsePolyZZ × Array SparsePolyZZ × Array Bool) :=
  validateCandidatesLoop ops candidates 0 activeLifted modulus fStar result
    (Array.replicate activeLifted.size false) activeLifted.size

/-- Raw C++ callees for one van-Hoeij iteration.  Candidate validation itself
is generated above and is no longer hidden behind a semantic callback. -/
structure VanHoeijRawOps where
  prepareCandidates : SparsePolyZZ → Array SparsePolyZZ → ZZ → Nat →
    RawExec (Array (Array Int32))
  validation : CandidateValidationRawOps
  zassenhaus : SparsePolyZZ → Array SparsePolyZZ → ZZ →
    RawExec (Array SparsePolyZZ)

/-- Erased termination certificate for the successful-extraction branch.
It refers only to concrete successful raw executions and the consumed bits
they returned. -/
structure VanHoeijTermination (ops : VanHoeijRawOps) where
  extraction_decreases : ∀ (lifted : Array SparsePolyZZ) (modulus : ZZ)
      (state : VanHoeijState) activeLifted candidates fStar' result'
      consumed active',
    gatherActive state.active lifted = .ok activeLifted →
    validateCandidates ops.validation state.fStar activeLifted modulus candidates state.result =
      .ok (fStar', result', consumed) →
    (∃ index, ∃ hindex : index < consumed.size, consumed[index] = true) →
    removeConsumed state.active consumed = .ok active' →
    active'.size < state.active.size

def removeConsumedDecreasing (ops : VanHoeijRawOps)
    (termination : VanHoeijTermination ops) (lifted : Array SparsePolyZZ)
    (modulus : ZZ)
    (state : VanHoeijState) (activeLifted : Array SparsePolyZZ)
    (candidates : Array (Array Int32)) (fStar' : SparsePolyZZ)
    (result' : Array SparsePolyZZ) (consumed : Array Bool)
    (hgather : gatherActive state.active lifted = .ok activeLifted)
    (hvalidate : validateCandidates ops.validation state.fStar activeLifted modulus
      candidates state.result = .ok (fStar', result', consumed))
    (hfound : ∃ index, ∃ hindex : index < consumed.size,
      consumed[index] = true) :
    RawExec { active' : Array Int32 // active'.size < state.active.size } :=
  match hremove : removeConsumed state.active consumed with
  | .error fault => .error fault
  | .ok active' => .ok ⟨active', termination.extraction_decreases lifted modulus
      state activeLifted candidates fStar' result' consumed active'
      hgather hvalidate hfound hremove⟩

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
      match hgather : gatherActive state.active lifted with
      | .error fault => .error fault
      | .ok activeLifted =>
        match ops.prepareCandidates state.fStar activeLifted modulus state.target with
        | .error fault => .error fault
        | .ok candidates =>
          match hvalidate : validateCandidates ops.validation state.fStar activeLifted modulus candidates
              state.result with
          | .error fault => .error fault
          | .ok (fStar', result', consumed) =>
            if hfound : ∃ index, ∃ hindex : index < consumed.size,
                consumed[index] = true then
              match removeConsumedDecreasing ops termination lifted modulus state
                  activeLifted candidates fStar' result' consumed hgather hvalidate hfound with
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
