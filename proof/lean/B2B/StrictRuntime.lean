import CLPoly.Refinement.FactorZp

/-!
Executable physical runtime for strict factorization B2B.

This module is test-only.  It constructs concrete raw heaps and distinct
regions for every generated pointer boundary.  `erasedProof` is used only for
`Prop` fields that Lean removes from executable code (the formal refinement
modules prove those fields independently).  It cannot manufacture an
executable value, polynomial, RNG state, or factorization result.
-/
set_option autoImplicit false

namespace B2B.StrictRuntime

open Refinement

/-- Test-runtime inhabitant for proof-erased physical contracts.  This module
must never be imported by `CLPoly.Refinement` or `CLPoly.Pipeline`. -/
unsafe def erasedValue.{u} {α : Sort u} : α :=
  unsafeCast.{1, u} Unit.unit

def rawPtr (region : Nat) : RawPtr UInt64 :=
  { region := some region, limbOffset := 0 }

def rawWord3Ptr (region : Nat) : RawPtr Word3 :=
  { region := some region, limbOffset := 0 }

def workspaceCapacity (left right : SparsePolyZp) : Nat :=
  let width := CLPoly.Impl.StrictPolynomialGCDRefinement.sparseDenseLength left +
    CLPoly.Impl.StrictPolynomialGCDRefinement.sparseDenseLength right + 1
  512 * width * width + 4096

def workspaceHeap (left right : SparsePolyZp) : RawHeap :=
  { regions := Array.replicate 24
      (Array.replicate (workspaceCapacity left right) 0) }

def workspaceMatrix : HgcdMat :=
  { poly := #[rawPtr 16, rawPtr 17, rawPtr 18, rawPtr 19]
    len := #[0, 0, 0, 0] }

unsafe def rawMulWorkspace (this : DenseUPolyZp)
    (left right : SparsePolyZp) :
    Refinement.StrictDDF.RawMulWorkspace this left right :=
  { heap := workspaceHeap left right
    leftPtr := rawPtr 1
    rightPtr := rawPtr 2
    outputPtr := rawPtr 0
    scratchPtr := rawPtr 10
    leftValid := erasedValue
    rightValid := erasedValue
    outputValid := erasedValue
    scratchValid := erasedValue
    inputsDisjoint := erasedValue
    outputLeftDisjoint := erasedValue
    outputRightDisjoint := erasedValue
    outputScratchDisjoint := erasedValue
    scratchLeftDisjoint := erasedValue
    scratchRightDisjoint := erasedValue
    leftLengthWord := erasedValue }

unsafe def rawModWorkspace (this : DenseUPolyZp)
    (dividend divisor : SparsePolyZp) :
    Refinement.StrictDDF.RawModWorkspace this dividend divisor :=
  { heap := workspaceHeap dividend divisor
    dividendPtr := rawPtr 1
    divisorPtr := rawPtr 2
    quotientPtr := rawPtr 3
    remainderPtr := rawPtr 4
    scratchPtr := rawWord3Ptr 5
    dividendValid := erasedValue
    divisorValid := erasedValue
    quotientValid := erasedValue
    remainderValid := erasedValue
    scratchValid := erasedValue
    inputsDisjoint := erasedValue
    remainderDividendDisjoint := erasedValue
    scratchDividendDisjoint := erasedValue
    scratchDivisorDisjoint := erasedValue
    quotientDivisorDisjoint := erasedValue
    quotientScratchDisjoint := erasedValue
    remainderScratchDisjoint := erasedValue
    remainderQuotientDisjoint := erasedValue
    remainderDivisorDisjoint := erasedValue
    quotientCapacity := erasedValue }

unsafe def rawGCDWorkspace (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)] (left right : SparsePolyZp) :
    Refinement.StrictDDF.RawGCDWorkspace this left right :=
  { M := workspaceMatrix
    hM := erasedValue
    resultPtr := rawPtr 0
    leftPtr := rawPtr 1
    rightPtr := rawPtr 2
    aBuf := rawPtr 3
    bBuf := rawPtr 4
    J := rawPtr 5
    Q := rawPtr 6
    R := rawPtr 7
    W3 := rawWord3Ptr 8
    W := rawPtr 9
    scratch := rawPtr 10
    euclidQ := rawPtr 11
    euclidR := rawPtr 12
    euclidW3 := rawWord3Ptr 13
    heap := workspaceHeap left right
    hcfg := erasedValue
    hp := erasedValue
    loopDecrease := erasedValue
    leftValid := erasedValue
    rightValid := erasedValue
    leftBound := erasedValue
    rightBound := erasedValue
    inputsDisjoint := erasedValue
    readyAB := erasedValue
    readyBA := erasedValue }

unsafe def yunRawGCDWorkspace (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (left right : SparsePolyZp) :
    Refinement.StrictSquarefreeZp.YunRawGCDWorkspace this hcfg left right :=
  { M := workspaceMatrix
    hM := erasedValue
    resultPtr := rawPtr 0
    leftPtr := rawPtr 1
    rightPtr := rawPtr 2
    aBuf := rawPtr 3
    bBuf := rawPtr 4
    J := rawPtr 5
    Q := rawPtr 6
    R := rawPtr 7
    W3 := rawWord3Ptr 8
    W := rawPtr 9
    scratch := rawPtr 10
    euclidQ := rawPtr 11
    euclidR := rawPtr 12
    euclidW3 := rawWord3Ptr 13
    heap := workspaceHeap left right
    loopDecrease := erasedValue
    leftValid := erasedValue
    rightValid := erasedValue
    disjoint := erasedValue
    readyAB := erasedValue
    readyBA := erasedValue }

unsafe def rawMulProvider (this : DenseUPolyZp) :
    Refinement.StrictDDF.RawMulWorkspaceProvider this :=
  { workspace := fun left right _ _ => rawMulWorkspace this left right }

unsafe def rawModProvider (this : DenseUPolyZp) (modulus : SparsePolyZp) :
    Refinement.StrictDDF.RawModWorkspaceProvider this modulus :=
  { workspace := fun dividend => rawModWorkspace this dividend modulus }

unsafe def ddfProviders (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)] :
    Refinement.StrictDDF.DDFRawProviders this :=
  { hcfg := erasedValue
    h2p := erasedValue
    mul := rawMulProvider this
    mod := rawModProvider this
    gcd := fun left right => rawGCDWorkspace this left right }

unsafe def sqfOps (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)] :
    Generated.StrictSquarefreeZp.SQFRawOps :=
  let hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this :=
    erasedValue
  { derivative := Refinement.StrictSquarefreeZp.derivativeIR this
    extractPthRoot := Refinement.StrictSquarefreeZp.extractPthRootIR
    makeMonic := Refinement.StrictSquarefreeGenerated.makeMonicRawIR this
    gcd := Refinement.StrictSquarefreeGenerated.gcdRawIR this hcfg
      (fun left right => yunRawGCDWorkspace this hcfg left right)
    exactDiv := Refinement.StrictSquarefreeZp.pairVecDivIR this
    EntryInvariant := fun _ => True
    YunInvariant := fun _ _ _ _ _ => True
    derivativeZeroStep := erasedValue
    yunInitial := erasedValue
    yunFactorStep := erasedValue
    yunNoFactorStep := erasedValue
    remainderRootStep := erasedValue }

unsafe def ddfOps (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)] : Generated.StrictDDF.DDFRawOps :=
  let providers := ddfProviders this
  { powmod := fun base exponent modulus =>
      Refinement.StrictDDF.strictPowmodIR this base exponent modulus
        providers.mul (providers.mod modulus)
    gcd := fun left right =>
      Refinement.StrictDDF.strictDDFGCDIR this left right
        (providers.gcd left right)
    makeMonic := Refinement.StrictDDF.strictMakeMonicIR this
    exactDiv := Refinement.StrictDDF.strictExactDivIR this
    mod := fun dividend divisor =>
      Refinement.StrictDDF.strictModIR this dividend divisor
        ((providers.mod divisor).workspace dividend)
    Invariant := fun _ _ _ _ _ => True
    splitStep := erasedValue
    noSplitStep := erasedValue }

/-- Deterministic placeholder state for executions whose EDF components take
the source early-return branch and therefore never request a random word. -/
unsafe def noRetryEngine : Generated.StrictEDF.RandomEngine Nat where
  nextAdvance := fun state upper =>
    if upper = 0 then (0, state + 1) else (state.toUInt64 % upper, state + 1)
  nextLt := erasedValue

unsafe def edfOps (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)] : Generated.StrictEDF.EDFRawOps Nat :=
  let providers := ddfProviders this
  { random := Generated.StrictEDF.__upoly_random_raw_ir noRetryEngine
    modPoly := fun dividend modulus =>
      Refinement.StrictDDF.strictModIR this dividend modulus
        ((providers.mod modulus).workspace dividend)
    squareAddMod := Refinement.StrictEDF.strictSquareAddModIR this providers
    powmod := fun base exponent modulus =>
      Refinement.StrictDDF.strictPowmodIR this base exponent.toNat modulus
        providers.mul (providers.mod modulus)
    gcd := fun left right =>
      Refinement.StrictDDF.strictDDFGCDIR this left right
        (providers.gcd left right)
    exactDiv := Refinement.StrictDDF.strictExactDivIR this
    makeMonic := Refinement.StrictDDF.strictMakeMonicIR this
    EntryInvariant := fun _ _ => True }

def insertByDegree (item : SparsePolyZp × UInt64)
    (items : List (SparsePolyZp × UInt64)) :
    List (SparsePolyZp × UInt64) :=
  match items with
  | [] => [item]
  | head :: tail =>
    if get_deg item.1 < get_deg head.1 then item :: head :: tail
    else head :: insertByDegree item tail

def sortByDegree (items : Array (SparsePolyZp × UInt64)) :
    RawExec (Array (SparsePolyZp × UInt64)) :=
  .ok (items.toList.foldr insertByDegree [] |>.toArray)

unsafe def factorZpNoRetryOps (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)] :
    Generated.StrictFactorZp.FactorZpRawOps Nat :=
  let sqf := sqfOps this
  let ddf := ddfOps this
  let edf := edfOps this
  { makeMonic := Refinement.StrictSquarefreeZp.upolyMakeMonicIR this
    squarefree := fun f =>
      Generated.StrictSquarefreeZp.__squarefree_Zp_raw_ir sqf f
        (fun _ => trivial)
    ddf := fun f =>
      Generated.StrictDDF.__ddf_Zp_raw_ir ddf f (fun _ => trivial)
    edf := fun result f d rng =>
      if (get_deg f).toUInt64 == d then
        Generated.StrictEDF.__edf_Zp_raw_ir edf
          { splitStep := erasedValue }
          { retryTrace := erasedValue }
          result f d rng trivial
      else
        .error .assertionFailure
    sortByDegree := sortByDegree }

unsafe def denseContext (p : UInt64) : DenseUPolyZp :=
  let norm := CLPoly.Impl.StrictWordArithmetic.denseNorm p
  { _coeffs := #[]
    _p := p
    _norm := norm
    _ninv := Generated.StrictGCD.dense_upoly_zp___preinvert_limb_ir
      (p <<< norm) }

/-- Execute the generated strict C++ entry on inputs for which every EDF
component takes its degree-equals-`d` source branch.  The name deliberately
records that this is not yet the general random-retry executor. -/
unsafe def factorZpNoRetry (f : SparsePolyZp) :
    RawExec (Zp × Array (SparsePolyZp × UInt64)) :=
  if hf : f.isEmpty then
    .error .assertionFailure
  else
    let ctx := denseContext f[0]!.2.prime
    letI : Fact (Nat.Prime ctx._p.toNat) := erasedValue
    Generated.StrictFactorZp.__factor_Zp_raw_ir
      (@factorZpNoRetryOps ctx inferInstance) 42 f

end B2B.StrictRuntime
