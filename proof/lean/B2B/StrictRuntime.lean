import CLPoly.Refinement.FactorZZ

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

structure MT19937State where
  words : Array UInt32
  index : Nat
deriving Repr

def mtSeed (seed : UInt32) : MT19937State :=
  let words := (Array.range 623).foldl (fun (words : Array UInt32) offset =>
    let index := offset + 1
    let previous := words[offset]!
    words.push ((1812433253 : UInt32) *
      (previous ^^^ (previous >>> (30 : UInt32))) + index.toUInt32))
    #[seed]
  { words, index := 624 }

def mtTwist (state : MT19937State) : MT19937State :=
  let words : Array UInt32 := (Array.range 624).map fun index =>
    let combined := (state.words[index]! &&& (0x80000000 : UInt32)) |||
      (state.words[(index + 1) % 624]! &&& (0x7fffffff : UInt32))
    let shifted := combined >>> (1 : UInt32)
    state.words[(index + 397) % 624]! ^^^ shifted ^^^
      (if combined &&& (1 : UInt32) = 0 then 0 else 0x9908b0df)
  { words, index := 0 }

def mtNextWord (state : MT19937State) : UInt32 × MT19937State :=
  let ready := if state.index < 624 then state else mtTwist state
  let word := ready.words[ready.index]!
  let word := word ^^^ (word >>> (11 : UInt32))
  let word := word ^^^ ((word <<< (7 : UInt32)) &&& 0x9d2c5680)
  let word := word ^^^ ((word <<< (15 : UInt32)) &&& 0xefc60000)
  let word := word ^^^ (word >>> (18 : UInt32))
  (word, { ready with index := ready.index + 1 })

def lowMask (width : Nat) : UInt64 :=
  if 64 ≤ width then ~~~(0 : UInt64) else (1 <<< width.toUInt64) - 1

/-- libc++'s small-range `uniform_int_distribution<uint64_t>` path: take the
minimum low-bit window and reject values outside the requested range. -/
unsafe def uniformBelow (state : MT19937State) (upper : UInt64) :
    UInt64 × MT19937State :=
  if upper ≤ 1 then (0, state)
  else
    let width := Nat.log2 (upper.toNat - 1) + 1
    if width ≤ 32 then
      let (word, next) := mtNextWord state
      let candidate := word.toUInt64 &&& lowMask width
      if candidate < upper then (candidate, next)
      else uniformBelow next upper
    else
      -- libc++ independent_bits_engine uses two MT words for a UInt64 range.
      let firstWidth := width / 2
      let secondWidth := width - firstWidth
      let (first, state') := uniformBelow state (1 <<< firstWidth.toUInt64)
      let (second, state'') := uniformBelow state' (1 <<< secondWidth.toUInt64)
      let candidate := (first <<< secondWidth.toUInt64) ||| second
      if candidate < upper then (candidate, state'')
      else uniformBelow state'' upper

unsafe def mt19937Engine : Generated.StrictEDF.RandomEngine MT19937State where
  nextAdvance := uniformBelow
  nextLt := erasedValue

unsafe def edfOps (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)] :
    Generated.StrictEDF.EDFRawOps MT19937State :=
  let providers := ddfProviders this
  { random := Generated.StrictEDF.__upoly_random_raw_ir mt19937Engine
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

/-- Run the literal EDF retry body until it produces a proper factor, and
record that concrete run as the finite trace consumed by the generated
well-founded EDF.  There is no retry counter and no candidate oracle. -/
unsafe def buildRetryTrace
    (ops : Generated.StrictEDF.EDFRawOps MT19937State)
    (f : SparsePolyZp) (d : UInt64) (rng : MT19937State) :
    Generated.StrictEDF.RetryTrace ops f d rng :=
  match hrandom : ops.random (get_deg f) f[0]!.2.prime rng with
  | .error _ => buildRetryTrace ops f d rng
  | .ok (r, rngNext) =>
    if hempty : r.isEmpty then
      .empty rng r rngNext hrandom hempty
        (buildRetryTrace ops f d rngNext)
    else
      let hbudget : 0 < d.toNat ∧ d.toNat < UInt64.size := erasedValue
      have hnonempty : r.isEmpty = false := by
        cases hvalue : r.isEmpty <;> simp_all
      match hcandidate : Generated.StrictEDF.candidateRun ops f d r hbudget with
      | .error _ => buildRetryTrace ops f d rng
      | .ok candidate =>
        if hproper : get_deg candidate > 0 ∧ get_deg candidate < get_deg f then
          .success rng r rngNext candidate hbudget hrandom hnonempty hcandidate
            erasedValue hproper
        else
          .failed rng r rngNext candidate hbudget hrandom hnonempty hcandidate
            hproper (buildRetryTrace ops f d rngNext)

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

unsafe def factorZpRuntimeOps (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)] :
    Generated.StrictFactorZp.FactorZpRawOps MT19937State :=
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
      Generated.StrictEDF.__edf_Zp_raw_ir edf
        { splitStep := erasedValue }
        { retryTrace := fun f d rng _ => buildRetryTrace edf f d rng }
        result f d rng trivial
    sortByDegree := sortByDegree }

unsafe def denseContext (p : UInt64) : DenseUPolyZp :=
  let norm := CLPoly.Impl.StrictWordArithmetic.denseNorm p
  { _coeffs := #[]
    _p := p
    _norm := norm
    _ninv := Generated.StrictGCD.dense_upoly_zp___preinvert_limb_ir
      (p <<< norm) }

/-- Execute the generated strict C++ entry with a concrete raw workspace and
an actual finite EDF retry trace assembled from the executed retry body. -/
unsafe def factorZpRuntime (f : SparsePolyZp) :
    RawExec (Zp × Array (SparsePolyZp × UInt64)) :=
  if hf : f.isEmpty then
    .error .assertionFailure
  else
    let ctx := denseContext f[0]!.2.prime
    letI : Fact (Nat.Prime ctx._p.toNat) := erasedValue
    Generated.StrictFactorZp.__factor_Zp_raw_ir
      (@factorZpRuntimeOps ctx inferInstance) (mtSeed 42) f

unsafe def selectPrimeCandidateOps (p : UInt64) :
    Generated.StrictSelectPrime.CandidateRawOps MT19937State :=
  let ctx := denseContext p
  letI : Fact (Nat.Prime ctx._p.toNat) := erasedValue
  let ddf := ddfOps ctx
  let edf := edfOps ctx
  { lcMod := fun coefficient modulus =>
      .ok (ZZ.fdiv_r 0 coefficient modulus.toNat)
    polynomialMod := Generated.StrictPolynomialMod.polynomial_mod_raw_ir
    derivative := fun source =>
      .ok (Refinement.StrictSquarefreeZp.derivativeIR ctx source)
    gcd := fun left right =>
      Refinement.StrictDDF.strictGCDIR ctx left right
        (rawGCDWorkspace ctx left right)
    makeMonic := Refinement.StrictSquarefreeZp.upolyMakeMonicIR ctx
    ddf := fun f =>
      Generated.StrictDDF.__ddf_Zp_raw_ir ddf f (fun _ => trivial)
    edf := fun result f d rng =>
      Generated.StrictEDF.__edf_Zp_raw_ir edf
        { splitStep := erasedValue }
        { retryTrace := fun f d rng _ => buildRetryTrace edf f d rng }
        result f d rng trivial }

unsafe def trySelectPrimeCandidate (f : SparsePolyZZ) (degree : Int64)
    (leading : ZZ) (p : UInt64) (rng : MT19937State) :
    RawExec (Generated.StrictSelectPrime.CandidateResult MT19937State) :=
  Generated.StrictSelectPrime.tryCandidateRaw (selectPrimeCandidateOps p)
    f degree leading p rng

/-- Logarithmic modular exponentiation used by the executable deterministic
Miller--Rabin implementation below.  `exponent` is the source loop variant,
not an artificial execution bound. -/
def powModNatLoop (modulus accumulator base exponent : Nat) : Nat :=
  if exponent = 0 then accumulator % modulus
  else
    let accumulator' :=
      if exponent % 2 = 1 then (accumulator * base) % modulus
      else accumulator
    powModNatLoop modulus accumulator' ((base * base) % modulus) (exponent / 2)
termination_by exponent
decreasing_by omega

def powModNat (base exponent modulus : Nat) : Nat :=
  powModNatLoop modulus 1 (base % modulus) exponent

def factorTwos (value : Nat) : Nat × Nat :=
  if value ≠ 0 && value % 2 = 0 then
    let result := factorTwos (value / 2)
    (result.1, result.2 + 1)
  else (value, 0)
termination_by value
decreasing_by
  simp_all only [Bool.and_eq_true, bne_iff_ne, decide_eq_true_eq]
  omega

def millerSquareLoop (n x rounds : Nat) : Bool :=
  if rounds = 0 then false
  else
    let x' := (x * x) % n
    if x' = n - 1 then true
    else millerSquareLoop n x' (rounds - 1)
termination_by rounds
decreasing_by omega

def millerWitness (n d s base : Nat) : Bool :=
  let x := powModNat (base % n) d n
  x = 1 || x = n - 1 || millerSquareLoop n x (s - 1)

/-- Deterministic Miller--Rabin for the complete UInt64 domain, using the
standard seven-base witness set.  This computes every candidate decision;
no selected-prime table or factorization answer is embedded in the runtime. -/
def isPrime64 (n : Nat) : Bool :=
  if n < 2 then false
  else if n = 2 || n = 3 then true
  else if n % 2 = 0 then false
  else
    let decomposition := factorTwos (n - 1)
    #[2, 325, 9375, 28178, 450775, 9780504, 1795265022].all
      (fun base => n ∣ base || millerWitness n decomposition.1 decomposition.2 base)

def nextPrimeFastSearch (candidate : Nat) : RawExec UInt64 :=
  if hbound : candidate < UInt64.size then
    if isPrime64 candidate then .ok candidate.toUInt64
    else nextPrimeFastSearch (candidate + 1)
  else .error .arithmeticOverflow
termination_by UInt64.size - candidate
decreasing_by omega

def prevPrimeFastSearch (candidate : Nat) : RawExec UInt64 :=
  if isPrime64 candidate then .ok candidate.toUInt64
  else if candidate = 0 then .error .arithmeticDomain
  else prevPrimeFastSearch (candidate - 1)
termination_by candidate
decreasing_by omega

def nextPrimeFast (useLargePrime : Bool) (p : UInt64) : RawExec UInt64 :=
  if useLargePrime then
    if p.toNat ≤ 2 then .error .arithmeticDomain
    else prevPrimeFastSearch (p.toNat - 1)
  else nextPrimeFastSearch (p.toNat + 1)

unsafe def selectPrimeRuntime (f : SparsePolyZZ) (useLargePrime : Bool) :
    RawExec PrimeSelectionResult :=
  let ops : Generated.StrictSelectPrime.SelectPrimeRawOps MT19937State :=
    { nextPrime := nextPrimeFast
      tryCandidate := trySelectPrimeCandidate }
  Generated.StrictSelectPrime.__select_prime_raw_ir ops
    { rank := Generated.StrictPrimeEnumeration.rank
      next_decreases := erasedValue }
    (mtSeed 42) useLargePrime f

unsafe def henselLiftRuntime (f : SparsePolyZZ)
    (factors : Array SparsePolyZp) (p : UInt64) (target : Int32) :
    RawExec (Array SparsePolyZZ × ZZ) :=
  let ctx := denseContext p
  letI : Fact (Nat.Prime ctx._p.toNat) := erasedValue
  Generated.StrictHensel.__hensel_lift_upoly_raw_ir
    (Refinement.StrictHensel.strictHenselRawOps
      Refinement.StrictHensel.concreteDivmodTermination)
    (Refinement.StrictHensel.strictHenselTreeBuildRawOps ctx
      (rawMulProvider ctx))
    f factors p target erasedValue

unsafe def factorZZRuntimeOps : Generated.StrictFactorZZ.FactorZZRawOps :=
  { selectPrime := selectPrimeRuntime
    henselLift := henselLiftRuntime
    vanHoeijRecombine := Refinement.StrictFactorZZ.concreteVanHoeijRecombine
    zassenhausRecombine :=
      Refinement.StrictFactorZZ.concreteZassenhausRecombine }

unsafe def factorZZRuntime (useLargePrime : Bool) (f : SparsePolyZZ) :
    RawExec (Array SparsePolyZZ) :=
  Generated.StrictFactorZZ.__factor_squarefree_primitive_ZZ_raw_ir
    factorZZRuntimeOps useLargePrime f

end B2B.StrictRuntime
