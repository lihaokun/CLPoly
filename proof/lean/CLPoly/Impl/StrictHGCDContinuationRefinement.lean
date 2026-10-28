import CLPoly.Impl.StrictHGCDRecursiveRefinement

set_option autoImplicit false

namespace CLPoly.Impl.StrictHGCDRawRefinement

open Generated.StrictHGCD
open CLPoly.Impl.RawPolynomialRep
open CLPoly.Impl.StrictDivremRefinement
open CLPoly.Impl.StrictEuclidRefinement
open CLPoly.Impl.StrictPolyAddSubRefinement
open CLPoly.Impl.StrictMulRefinement
open CLPoly.Impl.StrictWordArithmetic

/-- End-to-end semantic closure of the generated non-early tail, beginning
with the exact second cutoff dispatch and ending with the exact finish
record.  The middle quotient and low/high decompositions are consumed as
representations in the dispatch input heap and are transported only by the
dispatch's physical frame. -/
theorem hgcdRecursiveSecondDispatchFinish_rawInvariant
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (bound : Nat) (recurse : HgcdRecursiveCallBelow bound)
    (M R S : HgcdMat) (hM : M.Valid) (hR : R.Valid) (hS : S.Valid)
    (computeM : Bool)
    (A B T0 lowA lowB highA highB q inputA inputB : RawPtr UInt64)
    (lenLowA lenLowB shift lenQ : Nat)
    (Q2 : RawPtr UInt64) (W3 : RawPtr Word3)
    (T1 scratch stage WNext a2 : RawPtr UInt64)
    (lenInputA lenInputB inputLength : Nat) (heap : RawHeap)
    (second : HgcdRecursiveResult) (result : HgcdRecursiveFinishResult)
    (outerA outerB left right currentA currentB remainder quotient
      lowAPoly lowBPoly :
      Polynomial (ZMod this._p.toNat))
    (firstEntries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (sgnR : Int)
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (horder : lenInputB < lenInputA) (hdecrease : lenInputA < bound)
    (iterWorkspace : HgcdRecursiveDispatchIterWorkspace this S hS highA highB
      inputA inputB lenInputA lenInputB Q2 W3 T0 T1 scratch stage heap left
      right)
    (recursiveRefines : HgcdRecursiveCallbackRefinesAt this bound recurse S hS
      highA highB inputA inputB lenInputA lenInputB WNext scratch heap horder
      hdecrease left right)
    (frame : HgcdRecursiveSecondDispatchFrameProvider this bound recurse S hS
      highA highB inputA inputB lenInputA lenInputB Q2 W3 T0 T1 scratch stage
      WNext heap horder hdecrease R hR q lowA lowB lenQ lenLowA lenLowB)
    (hRRep : HgcdMatRawDenseRep this heap R firstEntries hR)
    (hQ : RawDensePolyRep this heap q lenQ quotient)
    (hLowA : RawCanonicalPolySlice this heap lowA lenLowA lowAPoly)
    (hLowB : RawCanonicalPolySlice this heap lowB lenLowB lowBPoly)
    (hFullA : currentB = lowAPoly + Polynomial.X ^ shift * left)
    (hFullB : remainder = lowBPoly + Polynomial.X ^ shift * right)
    (hFirstTransform : CLPoly.Impl.StrictHGCDRefinement.HgcdTransform
      outerA outerB currentA currentB (firstEntries 0) (firstEntries 1)
      (firstEntries 2) (firstEntries 3))
    (hFirstDet : CLPoly.Impl.StrictHGCDRefinement.HgcdSignedDet sgnR
      (firstEntries 0) (firstEntries 1) (firstEntries 2) (firstEntries 3))
    (hDivision : currentA = quotient * currentB + remainder)
    (finishWorkspace : HgcdRecursiveFinishWorkspaceProvider this M R
      second.matrix hM hR second.valid A B T0 lowA lowB highA highB q
      lenLowA lenLowB second.lenA second.lenB shift lenQ a2 scratch second.sgn
      second.heap)
    (hsecond : hgcdRecursiveDispatchBelow this bound recurse S hS highA highB
      inputA inputB lenInputA lenInputB Q2 W3 T0 T1 scratch stage WNext heap
      horder hdecrease = .ok second)
    (hfinish : hgcdRecursiveFinish this M R second.matrix hM hR second.valid
      computeM A B T0 lowA lowB highA highB q lenLowA lenLowB second.lenA
      second.lenB shift lenQ a2 scratch sgnR second.sgn second.heap =
        .ok result)
    (hstop : result.lenB < inputLength / 2 + 1)
    (hlength : computeM = true →
      HgcdRecursiveLengthInvariant inputLength result.toResult) :
    ∃ finalA finalB combinedEntries,
      HgcdRecursiveRawInvariant this outerA outerB finalA finalB
        combinedEntries computeM A B inputLength result.toResult := by
  rcases hgcdRecursiveSecondDispatch_refines this bound recurse S hS highA
      highB inputA inputB lenInputA lenInputB Q2 W3 T0 T1 scratch stage WNext
      heap left right R hR q lowA lowB lenQ lenLowA lenLowB firstEntries
      quotient lowAPoly lowBPoly hcfg hp horder hdecrease iterWorkspace
      recursiveRefines frame hRRep hQ hLowA hLowB second hsecond with
    ⟨outputHighA, outputHighB, secondEntries, hSecond, hRAfter, hQAfter,
      hLowAAfter, hLowBAfter⟩
  have hFinish := hgcdRecursiveFinish_preserves_input this M R second.matrix
    hM hR second.valid computeM A B T0 lowA lowB highA highB q lenLowA
    lenLowB second.lenA second.lenB shift lenQ lenInputA a2 scratch sgnR
    second.sgn second.heap result firstEntries secondEntries quotient currentB
    remainder lowAPoly lowBPoly left right outputHighA outputHighB second hcfg
    hp finishWorkspace hRAfter (hSecond.matrixSemantics rfl).1 hQAfter
    hLowAAfter hLowBAfter hSecond rfl rfl (HEq.rfl) rfl rfl rfl hFullA
    hFullB hfinish
  let finalA := hgcdReconstructedLowA secondEntries lowAPoly lowBPoly
      second.sgn + Polynomial.X ^ shift * outputHighA
  let finalB := hgcdReconstructedLowB secondEntries lowAPoly lowBPoly
      second.sgn + Polynomial.X ^ shift * outputHighB
  let combinedEntries := hgcdMatProductEntry firstEntries
    (hgcdMatApplyQuotientEntries secondEntries quotient)
  have hInvariant := hgcdRecursiveRawInvariant_of_finish_execution this
    outerA outerB currentA currentB remainder quotient finalA finalB
    firstEntries secondEntries sgnR second.sgn computeM A B inputLength result
    (by simpa [finalA] using hFinish.1)
    (by simpa [finalB] using hFinish.2.1)
    hFinish.2.2.2.2.2 hFirstTransform hFirstDet hDivision
    (by simpa [finalA, finalB] using hFinish.2.2.1)
    hFinish.2.2.2.1 hFinish.2.2.2.2.1 hstop hlength
  exact ⟨finalA, finalB, combinedEntries, by
    simpa [combinedEntries] using hInvariant⟩

set_option maxHeartbeats 1200000 in
/-- Complete length invariant of the real matrix-producing non-early tail.
The proof uses the exact reconstruction and combine executions exposed by
`hgcdRecursiveFinish_exec`; every descriptor bound is transported through
their concrete heaps and returned records. -/
theorem hgcdRecursiveFinish_lengthInvariant (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (M : HgcdMat) (hM : M.Valid)
    (A B T0 lowA lowB highA highB q : RawPtr UInt64)
    (a2 scratch : RawPtr UInt64)
    (outerLength highLength m reconstructedLenA reconstructedLenB lenD k
      secondInputLength lenQ : Nat)
    (first second : HgcdRecursiveResult)
    (result : HgcdRecursiveFinishResult)
    (right entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (quotient polyLowA polyLowB polyHighA polyHighB :
      Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (firstLength : HgcdRecursiveLengthInvariant highLength first)
    (secondLength : HgcdRecursiveLengthInvariant secondInputLength second)
    (physical : HgcdRecursiveFinishWorkspaceProvider this M first.matrix
      second.matrix hM first.valid second.valid A B T0 lowA lowB highA highB
      q (Nat.min reconstructedLenB k) (Nat.min lenD k) second.lenA
      second.lenB k lenQ a2 scratch second.sgn second.heap)
    (hRRep : HgcdMatRawDenseRep this second.heap first.matrix right first.valid)
    (hSRep : HgcdMatRawDenseRep this second.heap second.matrix entries
      second.valid)
    (hQ : RawDensePolyRep this second.heap q lenQ quotient)
    (hLowA : RawCanonicalPolySlice this second.heap lowA
      (Nat.min reconstructedLenB k) polyLowA)
    (hLowB : RawCanonicalPolySlice this second.heap lowB
      (Nat.min lenD k) polyLowB)
    (hHighA : RawDensePolyRep this second.heap highA second.lenA polyHighA)
    (hHighB : RawDensePolyRep this second.heap highB second.lenB polyHighB)
    (hm : m = outerLength / 2)
    (hhigh : highLength = outerLength - m)
    (hreconstructedA : reconstructedLenA = m + first.lenA)
    (hreconstructedOrder : reconstructedLenB ≤ reconstructedLenA)
    (hreconstructedLower : m + 1 ≤ reconstructedLenB)
    (hreconstructedUpper : reconstructedLenB < outerLength)
    (hlenQ : lenQ ≤ reconstructedLenA - (reconstructedLenB - 1))
    (hk : k = 2 * m - reconstructedLenB + 1)
    (hc : secondInputLength = reconstructedLenB - k)
    (hrun : hgcdRecursiveFinish this M first.matrix second.matrix hM
      first.valid second.valid true A B T0 lowA lowB highA highB q
      (Nat.min reconstructedLenB k) (Nat.min lenD k) second.lenA second.lenB
      k lenQ a2 scratch first.sgn second.sgn second.heap = .ok result) :
    HgcdRecursiveLengthInvariant outerLength result.toResult := by
  rcases hgcdRecursiveFinish_exec this M first.matrix second.matrix hM
      first.valid second.valid true A B T0 lowA lowB highA highB q
      (Nat.min reconstructedLenB k) (Nat.min lenD k) second.lenA second.lenB
      k lenQ a2 scratch first.sgn second.sgn second.heap result hrun with
    ⟨reconstructed, hreconstruct, hlenA, _, _, htail⟩
  simp only [if_pos] at htail
  rcases htail with ⟨combined, hcombine, _, hmatrix⟩
  have hwork := physical reconstructed hreconstruct
  have hcombineWork := hwork.combine combined hcombine
  have hOperands := hgcdRecursiveFinish_operandInvariant this M first.matrix
    hM first.valid true A B T0 lowA lowB highA highB q outerLength m
    reconstructedLenB lenD k secondInputLength lenQ a2 scratch first.sgn
    second result entries polyLowA polyLowB polyHighA polyHighB hcfg hp
    secondLength hwork.reconstruct hSRep hLowA hLowB hHighA hHighB hm hk hc
    hreconstructedLower hreconstructedUpper hrun
  have hRReconstructed : HgcdMatRawDenseRep this reconstructed.heap
      first.matrix right first.valid := by
    intro i
    exact rawDensePolyRep_of_same_prefix this second.heap reconstructed.heap
      (hgcdMatPtr first.matrix first.valid i)
      (hgcdMatLen first.matrix first.valid i) (right i) hwork.afterLayout
      (hwork.rightPrefix i) (hRRep i)
  have hSReconstructed : HgcdMatRawDenseRep this reconstructed.heap
      second.matrix entries second.valid := by
    intro i
    exact rawDensePolyRep_of_same_prefix this second.heap reconstructed.heap
      (hgcdMatPtr second.matrix second.valid i)
      (hgcdMatLen second.matrix second.valid i) (entries i)
      hwork.afterLayout (hwork.secondPrefix i) (hSRep i)
  have hQReconstructed := rawDensePolyRep_of_same_prefix this second.heap
    reconstructed.heap q lenQ quotient hwork.afterLayout hwork.quotientPrefix
    hQ
  have horder' : reconstructedLenB ≤ m + first.lenA := by
    rwa [← hreconstructedA]
  have hq' : lenQ ≤ m + first.lenA - (reconstructedLenB - 1) := by
    rwa [← hreconstructedA]
  have houter : highLength + m = outerLength := by
    omega
  have hsplit : k + secondInputLength = reconstructedLenB := by
    omega
  have hleading := hgcdRecursiveFinalReconstruct_lenA_eq_of_invariant this A
    B T0 lowA lowB highA highB scratch (Nat.min reconstructedLenB k)
    (Nat.min lenD k) k secondInputLength second reconstructed entries polyLowA
    polyLowB polyHighA polyHighB hcfg hp secondLength (Nat.min_le_right _ _)
    (Nat.min_le_right _ _) hwork.reconstruct hSRep hLowA hLowB hHighA hHighB
    hreconstruct
  have hRows := hgcdRecursiveCombineMatrix_rowA_bounds this M first.matrix
    second.matrix hM first.valid second.valid q lenQ T0 a2 scratch
    reconstructed.heap combined right entries quotient outerLength highLength m
    first.lenA reconstructedLenB k secondInputLength second.lenA
    reconstructed.lenA hcfg hp hcombineWork.1 hRReconstructed hSReconstructed
    hQReconstructed houter
    (fun i => by
      fin_cases i
      · simpa [hgcdMatLen, hgcdMatLenRaw] using firstLength.row0A
      · simpa [hgcdMatLen, hgcdMatLenRaw] using firstLength.row1A
      · simpa [hgcdMatLen, hgcdMatLenRaw] using firstLength.row2A
      · simpa [hgcdMatLen, hgcdMatLenRaw] using firstLength.row3A)
    firstLength.inputBound horder' hreconstructedLower hq' hsplit
    (fun i => by
      fin_cases i
      · simpa [hgcdMatLen, hgcdMatLenRaw] using secondLength.row0A
      · simpa [hgcdMatLen, hgcdMatLenRaw] using secondLength.row1A
      · simpa [hgcdMatLen, hgcdMatLenRaw] using secondLength.row2A
      · simpa [hgcdMatLen, hgcdMatLenRaw] using secondLength.row3A)
    hleading.1 hcombine
  have hCoeff := hgcdRecursiveCombineMatrix_coeff_bounds this M first.matrix
    second.matrix hM first.valid second.valid q lenQ T0 a2 scratch
    reconstructed.heap combined right entries quotient outerLength m highLength
    first.lenA reconstructedLenA reconstructedLenB k secondInputLength hcfg hp
    hcombineWork.1 hRReconstructed hSReconstructed hQReconstructed hm hhigh
    (fun i => by
      fin_cases i
      · simpa [hgcdMatLen, hgcdMatLenRaw] using firstLength.row0A
      · simpa [hgcdMatLen, hgcdMatLenRaw] using firstLength.row1A
      · simpa [hgcdMatLen, hgcdMatLenRaw] using firstLength.row2A
      · simpa [hgcdMatLen, hgcdMatLenRaw] using firstLength.row3A)
    hreconstructedA hreconstructedOrder hreconstructedLower hlenQ hk hc
    (fun i => by
      simpa [hgcdMatLen, hgcdMatLenRaw] using secondLength.coeffBound i)
    hcombine
  rcases hRows with ⟨hCombinedValid, hCombinedRows⟩
  rcases hCoeff with ⟨hCombinedValid', hCombinedCoeff⟩
  apply hgcdRecursiveLengthInvariant_of_finish outerLength result hOperands
  · intro i
    simpa [hmatrix, hlenA, hgcdMatLen] using hCombinedRows i
  · intro i
    simpa [hmatrix, hgcdMatLen] using hCombinedCoeff i

/-- The generated recursive-HGCD base helper with matrix computation enabled
is exactly the already-refined iterator initialization prefix, modulo its
smaller return record. -/
theorem hgcdRecursiveBase_true_eq_init (M : HgcdMat)
    (A B a b : RawPtr UInt64) (lenA lenB : Nat) (heap : RawHeap) :
    hgcdRecursiveBase M true A B a b lenA lenB heap =
      (hgcdIterInit M A B A A 0 a lenA b lenB heap).map
        HgcdIterState.toRecursiveBaseResult := by
  unfold hgcdRecursiveBase hgcdIterInit
  generalize hone : dense_upoly_zp__mat_one_ir M heap = one
  cases one with
  | error fault =>
    simp only [hone, ↓reduceIte, Except.map]
  | ok pair =>
    rcases pair with ⟨heap1, matrix⟩
    generalize hcopyA : heap1.copyU64 A a lenA = copyA
    cases copyA with
    | error fault =>
      simp only [hone, hcopyA, ↓reduceIte, Except.map]
    | ok heap2 =>
      generalize hcopyB : heap2.copyU64 B b lenB = copyB
      cases copyB with
      | error fault =>
        simp only [hone, hcopyA, hcopyB, ↓reduceIte, Except.map]
      | ok heap3 =>
        simp only [hone, hcopyA, hcopyB, ↓reduceIte, Except.map,
          HgcdIterState.toRecursiveBaseResult]

/-- Semantic refinement of the exact `_hgcd_recursive` base branch when the
source requests its matrix.  The two raw copies retain their C++ order. -/
theorem hgcdRecursiveBase_true_refines (this : DenseUPolyZp)
    (M : HgcdMat) (A B a b : RawPtr UInt64) (lenA lenB : Nat)
    (heap : RawHeap) (result : HgcdRecursiveBaseResult)
    (left right : Polynomial (ZMod this._p.toNat)) (hM : M.Valid)
    (hp : 1 < this._p.toNat)
    (h0 : heap.ValidU64Slice (hgcdMatPtr M hM (0 : Fin 4)) 1)
    (h3 : heap.ValidU64Slice (hgcdMatPtr M hM (3 : Fin 4)) 1)
    (h03 : U64SlicesDisjoint (hgcdMatPtr M hM (0 : Fin 4)) 1
      (hgcdMatPtr M hM (3 : Fin 4)) 1)
    (hA : heap.ValidU64Slice A lenA) (hB : heap.ValidU64Slice B lenB)
    (hAa : U64SlicesDisjoint A lenA a lenA)
    (hBb : U64SlicesDisjoint B lenB b lenB)
    (hAb : U64SlicesDisjoint A lenA b lenB)
    (hBA : U64SlicesDisjoint B lenB A lenA)
    (h0a : U64SlicesDisjoint (hgcdMatPtr M hM (0 : Fin 4)) 1 a lenA)
    (h3a : U64SlicesDisjoint (hgcdMatPtr M hM (3 : Fin 4)) 1 a lenA)
    (h0b : U64SlicesDisjoint (hgcdMatPtr M hM (0 : Fin 4)) 1 b lenB)
    (h3b : U64SlicesDisjoint (hgcdMatPtr M hM (3 : Fin 4)) 1 b lenB)
    (hAMatrix : ∀ i : Fin 4, U64SlicesDisjoint A lenA
      (hgcdMatPtr M hM i) (identityEntryLen i))
    (hBMatrix : ∀ i : Fin 4, U64SlicesDisjoint B lenB
      (hgcdMatPtr M hM i) (identityEntryLen i))
    (hMatrixValid : ∀ i : Fin 4, heap.ValidU64Slice
      (hgcdMatPtr M hM i) (identityEntryLen i))
    (hLeft : RawDensePolyRep this heap a lenA left)
    (hRight : RawDensePolyRep this heap b lenB right)
    (horder : lenB ≤ lenA) (_hlenAPos : 0 < lenA)
    (hstop : lenB < lenA / 2 + 1)
    (hrun : hgcdRecursiveBase M true A B a b lenA lenB heap = .ok result) :
    result.lenA = lenA ∧ result.lenB = lenB ∧ result.sgn = 1 ∧
      ∃ hResultM : result.matrix.Valid,
        HgcdMatRawDenseRep this result.heap result.matrix
          (identityEntries this._p.toNat) hResultM ∧
        RawDensePolyRep this result.heap A result.lenA left ∧
        RawDensePolyRep this result.heap B result.lenB right ∧
        CLPoly.Impl.StrictHGCDRefinement.HgcdTransform left right left right
          (identityEntries this._p.toNat 0)
          (identityEntries this._p.toNat 1)
          (identityEntries this._p.toNat 2)
          (identityEntries this._p.toNat 3) ∧
        CLPoly.Impl.StrictHGCDRefinement.HgcdSignedDet result.sgn
          (identityEntries this._p.toNat 0)
          (identityEntries this._p.toNat 1)
          (identityEntries this._p.toNat 2)
          (identityEntries this._p.toNat 3) ∧
        HgcdRecursiveLengthInvariant lenA (result.toResult hResultM) := by
  rcases hgcdIterInit_refines this M A B A A 0 a lenA b lenB heap left
      right hM hp h0 h3 h03 hA hB hAa hBb hAb hBA h0a h3a h0b h3b
      hAMatrix hBMatrix hMatrixValid hLeft hRight with
    ⟨initial, hinit, hIA, hlenIA, hIB, hlenIB, _, _, _, hsgn,
      hInitialM, hMatrix,
      hARep, hBRep, hTransform, hDet⟩
  have hMatrixLengths := hgcdIterInit_matrixLengthInvariant M A B A A 0 a
    lenA b lenB heap initial hInitialM horder hinit
  have hbase : hgcdRecursiveBase M true A B a b lenA lenB heap =
      .ok initial.toRecursiveBaseResult := by
    rw [hgcdRecursiveBase_true_eq_init, hinit]
    rfl
  have heq := Except.ok.inj (hbase.symm.trans hrun)
  subst result
  refine ⟨hlenIA, hlenIB, hsgn, hInitialM, hMatrix, ?_, ?_,
    hTransform, ?_, ?_⟩
  · simpa [HgcdIterState.toRecursiveBaseResult, hIA, hlenIA] using hARep
  · simpa [HgcdIterState.toRecursiveBaseResult, hIB, hlenIB] using hBRep
  · simpa [HgcdIterState.toRecursiveBaseResult] using hDet
  · constructor
    · simpa [HgcdRecursiveBaseResult.toResult, HgcdIterState.toRecursiveBaseResult,
        hgcdMatLen, hgcdMatLenRaw] using hMatrixLengths.row0A
    · simpa [HgcdRecursiveBaseResult.toResult, HgcdIterState.toRecursiveBaseResult,
        hgcdMatLen, hgcdMatLenRaw] using hMatrixLengths.row1B
    · simpa [HgcdRecursiveBaseResult.toResult, HgcdIterState.toRecursiveBaseResult,
        hgcdMatLen, hgcdMatLenRaw] using hMatrixLengths.row2A
    · simpa [HgcdRecursiveBaseResult.toResult, HgcdIterState.toRecursiveBaseResult,
        hgcdMatLen, hgcdMatLenRaw] using hMatrixLengths.row3B
    · simpa [HgcdRecursiveBaseResult.toResult, HgcdIterState.toRecursiveBaseResult,
        hgcdMatLen, hgcdMatLenRaw] using hMatrixLengths.row1A
    · simpa [HgcdRecursiveBaseResult.toResult, HgcdIterState.toRecursiveBaseResult,
        hgcdMatLen, hgcdMatLenRaw] using hMatrixLengths.row3A
    · simp [HgcdRecursiveBaseResult.toResult,
        HgcdIterState.toRecursiveBaseResult, hlenIA]
    · simpa [HgcdRecursiveBaseResult.toResult, HgcdIterState.toRecursiveBaseResult,
        hlenIB] using hstop
    · simpa [HgcdRecursiveBaseResult.toResult, HgcdIterState.toRecursiveBaseResult,
        hlenIA] using _hlenAPos
    · simp [HgcdRecursiveBaseResult.toResult,
        HgcdIterState.toRecursiveBaseResult, hlenIA]
      omega
    · intro i
      simpa [HgcdRecursiveBaseResult.toResult, HgcdIterState.toRecursiveBaseResult,
        hgcdMatLenRaw, hgcdMatLen] using
        (hgcdIterInit_matrixCoefficientBound M A B A A 0 a lenA b lenB heap
          initial hInitialM _hlenAPos hinit i)

/-- Semantic refinement of the exact `_hgcd_recursive` base branch used by
GCD when matrix output is disabled.  No matrix call or matrix specification
is introduced; only the two source-ordered raw copies are executed. -/
theorem hgcdRecursiveBase_false_refines (this : DenseUPolyZp)
    (M : HgcdMat) (A B a b : RawPtr UInt64) (lenA lenB : Nat)
    (heap : RawHeap) (result : HgcdRecursiveBaseResult)
    (left right : Polynomial (ZMod this._p.toNat))
    (hA : heap.ValidU64Slice A lenA) (hB : heap.ValidU64Slice B lenB)
    (hAa : U64SlicesDisjoint A lenA a lenA)
    (hBb : U64SlicesDisjoint B lenB b lenB)
    (hAb : U64SlicesDisjoint A lenA b lenB)
    (hBA : U64SlicesDisjoint B lenB A lenA)
    (hLeft : RawDensePolyRep this heap a lenA left)
    (hRight : RawDensePolyRep this heap b lenB right)
    (hrun : hgcdRecursiveBase M false A B a b lenA lenB heap = .ok result) :
    result.matrix = M ∧ result.lenA = lenA ∧ result.lenB = lenB ∧
      result.sgn = 1 ∧
      RawDensePolyRep this result.heap A result.lenA left ∧
      RawDensePolyRep this result.heap B result.lenB right := by
  rcases copyU64_refines_rawDense this heap A a lenA left hA hAa hLeft with
    ⟨heap1, hcopyA, hlayout1, hA1⟩
  have hsameB1 := copyU64_preserves_prefix heap heap1 A a b lenA lenB hA
    hLeft.1 hAb hcopyA
  have hRight1 := rawDensePolyRep_of_same_prefix this heap heap1 b lenB
    right hlayout1 hsameB1 hRight
  have hB1 : heap1.ValidU64Slice B lenB := (hlayout1 B lenB).mp hB
  rcases copyU64_refines_rawDense this heap1 B b lenB right hB1 hBb
      hRight1 with ⟨heap2, hcopyB, hlayout2, hB2⟩
  have hsameA2 := copyU64_preserves_prefix heap1 heap2 B b A lenB lenA
    hB1 hRight1.1 hBA hcopyB
  have hA2 := rawDensePolyRep_of_same_prefix this heap1 heap2 A lenA left
    hlayout2 hsameA2 hA1
  have hactual : hgcdRecursiveBase M false A B a b lenA lenB heap =
      .ok (.mk heap2 M lenA lenB 1) := by
    simp [hgcdRecursiveBase, hcopyA, hcopyB]
  have heq : HgcdRecursiveBaseResult.mk heap2 M lenA lenB 1 = result :=
    Except.ok.inj (hactual.symm.trans hrun)
  subst result
  exact ⟨rfl, rfl, rfl, rfl, hA2, hB2⟩

/-- The matrix-producing base arm establishes the same semantic package
consumed by a recursive parent.  Every field comes from `_mat_one` and the
two concrete source copies executed by `hgcdRecursiveBase`. -/
theorem hgcdRecursiveBase_true_rawInvariant (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (M : HgcdMat) (A B a b : RawPtr UInt64) (lenA lenB : Nat)
    (heap : RawHeap) (result : HgcdRecursiveBaseResult)
    (left right : Polynomial (ZMod this._p.toNat)) (hM : M.Valid)
    (hp : 1 < this._p.toNat)
    (h0 : heap.ValidU64Slice (hgcdMatPtr M hM (0 : Fin 4)) 1)
    (h3 : heap.ValidU64Slice (hgcdMatPtr M hM (3 : Fin 4)) 1)
    (h03 : U64SlicesDisjoint (hgcdMatPtr M hM (0 : Fin 4)) 1
      (hgcdMatPtr M hM (3 : Fin 4)) 1)
    (hA : heap.ValidU64Slice A lenA) (hB : heap.ValidU64Slice B lenB)
    (hAa : U64SlicesDisjoint A lenA a lenA)
    (hBb : U64SlicesDisjoint B lenB b lenB)
    (hAb : U64SlicesDisjoint A lenA b lenB)
    (hBA : U64SlicesDisjoint B lenB A lenA)
    (h0a : U64SlicesDisjoint (hgcdMatPtr M hM (0 : Fin 4)) 1 a lenA)
    (h3a : U64SlicesDisjoint (hgcdMatPtr M hM (3 : Fin 4)) 1 a lenA)
    (h0b : U64SlicesDisjoint (hgcdMatPtr M hM (0 : Fin 4)) 1 b lenB)
    (h3b : U64SlicesDisjoint (hgcdMatPtr M hM (3 : Fin 4)) 1 b lenB)
    (hAMatrix : ∀ i : Fin 4, U64SlicesDisjoint A lenA
      (hgcdMatPtr M hM i) (identityEntryLen i))
    (hBMatrix : ∀ i : Fin 4, U64SlicesDisjoint B lenB
      (hgcdMatPtr M hM i) (identityEntryLen i))
    (hMatrixValid : ∀ i : Fin 4, heap.ValidU64Slice
      (hgcdMatPtr M hM i) (identityEntryLen i))
    (hLeft : RawDensePolyRep this heap a lenA left)
    (hRight : RawDensePolyRep this heap b lenB right)
    (horder : lenB ≤ lenA) (hlenAPos : 0 < lenA)
    (hstop : lenB < lenA / 2 + 1)
    (hrun : hgcdRecursiveBase M true A B a b lenA lenB heap = .ok result) :
    ∃ hResultM : result.matrix.Valid,
      HgcdRecursiveRawInvariant this left right left right
        (identityEntries this._p.toNat) true A B lenA
        (result.toResult hResultM) := by
  rcases hgcdRecursiveBase_true_refines this M A B a b lenA lenB heap result
      left right hM hp h0 h3 h03 hA hB hAa hBb hAb hBA h0a h3a h0b h3b
      hAMatrix hBMatrix hMatrixValid hLeft hRight horder hlenAPos hstop hrun with
    ⟨hlenA, hlenB, hsgn, hResultM, hMatrix, hARep, hBRep, hTransform,
      hDet, hLengths⟩
  refine ⟨hResultM, ?_⟩
  constructor
  · exact hARep
  · exact hBRep
  · intro _
    exact ⟨hMatrix, hTransform, hDet⟩
  · rfl
  · simpa [HgcdRecursiveBaseResult.toResult, hlenB] using hstop
  · intro _
    exact hLengths

/-- The matrix-disabled base arm establishes the recursive semantic package
without claiming matrix contents that the source did not compute. -/
theorem hgcdRecursiveBase_false_rawInvariant (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (M : HgcdMat) (A B a b : RawPtr UInt64) (lenA lenB : Nat)
    (heap : RawHeap) (result : HgcdRecursiveBaseResult)
    (left right : Polynomial (ZMod this._p.toNat)) (hM : M.Valid)
    (hA : heap.ValidU64Slice A lenA) (hB : heap.ValidU64Slice B lenB)
    (hAa : U64SlicesDisjoint A lenA a lenA)
    (hBb : U64SlicesDisjoint B lenB b lenB)
    (hAb : U64SlicesDisjoint A lenA b lenB)
    (hBA : U64SlicesDisjoint B lenB A lenA)
    (hLeft : RawDensePolyRep this heap a lenA left)
    (hRight : RawDensePolyRep this heap b lenB right)
    (hstop : lenB < lenA / 2 + 1)
    (hrun : hgcdRecursiveBase M false A B a b lenA lenB heap = .ok result) :
    ∃ hResultM : result.matrix.Valid,
      HgcdRecursiveRawInvariant this left right left right
        (identityEntries this._p.toNat) false A B lenA
        (result.toResult hResultM) := by
  rcases hgcdRecursiveBase_false_refines this M A B a b lenA lenB heap result
      left right hA hB hAa hBb hAb hBA hLeft hRight hrun with
    ⟨hmatrix, hlenA, hlenB, hsgn, hARep, hBRep⟩
  have hResultM : result.matrix.Valid := by simpa [hmatrix] using hM
  refine ⟨hResultM, ?_⟩
  constructor
  · exact hARep
  · exact hBRep
  · intro hfalse
    simp at hfalse
  · rfl
  · simpa [HgcdRecursiveBaseResult.toResult, hlenB] using hstop
  · intro hfalse
    simp at hfalse

/-- Uniform refinement theorem for the exact source base branch.  Splitting
on `computeM` follows the generated branch itself; the disabled arm does not
inherit any matrix claim from the enabled arm. -/
theorem hgcdRecursiveBase_rawInvariant (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (M : HgcdMat) (computeM : Bool)
    (A B a b : RawPtr UInt64) (lenA lenB : Nat)
    (heap : RawHeap) (result : HgcdRecursiveBaseResult)
    (left right : Polynomial (ZMod this._p.toNat)) (hM : M.Valid)
    (hp : 1 < this._p.toNat)
    (h0 : heap.ValidU64Slice (hgcdMatPtr M hM (0 : Fin 4)) 1)
    (h3 : heap.ValidU64Slice (hgcdMatPtr M hM (3 : Fin 4)) 1)
    (h03 : U64SlicesDisjoint (hgcdMatPtr M hM (0 : Fin 4)) 1
      (hgcdMatPtr M hM (3 : Fin 4)) 1)
    (hA : heap.ValidU64Slice A lenA) (hB : heap.ValidU64Slice B lenB)
    (hAa : U64SlicesDisjoint A lenA a lenA)
    (hBb : U64SlicesDisjoint B lenB b lenB)
    (hAb : U64SlicesDisjoint A lenA b lenB)
    (hBA : U64SlicesDisjoint B lenB A lenA)
    (h0a : U64SlicesDisjoint (hgcdMatPtr M hM (0 : Fin 4)) 1 a lenA)
    (h3a : U64SlicesDisjoint (hgcdMatPtr M hM (3 : Fin 4)) 1 a lenA)
    (h0b : U64SlicesDisjoint (hgcdMatPtr M hM (0 : Fin 4)) 1 b lenB)
    (h3b : U64SlicesDisjoint (hgcdMatPtr M hM (3 : Fin 4)) 1 b lenB)
    (hAMatrix : ∀ i : Fin 4, U64SlicesDisjoint A lenA
      (hgcdMatPtr M hM i) (identityEntryLen i))
    (hBMatrix : ∀ i : Fin 4, U64SlicesDisjoint B lenB
      (hgcdMatPtr M hM i) (identityEntryLen i))
    (hMatrixValid : ∀ i : Fin 4, heap.ValidU64Slice
      (hgcdMatPtr M hM i) (identityEntryLen i))
    (hLeft : RawDensePolyRep this heap a lenA left)
    (hRight : RawDensePolyRep this heap b lenB right)
    (horder : lenB ≤ lenA) (hlenAPos : 0 < lenA)
    (hstop : lenB < lenA / 2 + 1)
    (hrun : hgcdRecursiveBase M computeM A B a b lenA lenB heap =
      .ok result) :
    ∃ entries hResultM,
      HgcdRecursiveRawInvariant this left right left right entries computeM
        A B lenA (result.toResult hResultM) := by
  cases computeM with
  | false =>
      rcases hgcdRecursiveBase_false_rawInvariant this M A B a b lenA lenB
        heap result left right hM hA hB hAa hBb hAb hBA hLeft hRight hstop
        hrun with ⟨hResultM, hinvariant⟩
      exact ⟨identityEntries this._p.toNat, hResultM, hinvariant⟩
  | true =>
      rcases hgcdRecursiveBase_true_rawInvariant this M A B a b lenA lenB
        heap result left right hM hp h0 h3 h03 hA hB hAa hBb hAb hBA h0a
        h3a h0b h3b hAMatrix hBMatrix hMatrixValid hLeft hRight horder
        hlenAPos hstop hrun with ⟨hResultM, hinvariant⟩
      exact ⟨identityEntries this._p.toNat, hResultM, hinvariant⟩

/-- Physical facts for the exact base arm of the complete well-founded body.
They describe only generated matrix initialization and the two source-order
copies; no fact about a precomputed result is stored here. -/
structure HgcdRecursiveBaseCallWorkspace (this : DenseUPolyZp)
    (M : HgcdMat) (hM : M.Valid) (A B a b : RawPtr UInt64)
    (lenA lenB : Nat) (heap : RawHeap)
    (left right : Polynomial (ZMod this._p.toNat)) : Prop where
  valid0 : heap.ValidU64Slice (hgcdMatPtr M hM (0 : Fin 4)) 1
  valid3 : heap.ValidU64Slice (hgcdMatPtr M hM (3 : Fin 4)) 1
  disjoint03 : U64SlicesDisjoint (hgcdMatPtr M hM (0 : Fin 4)) 1
    (hgcdMatPtr M hM (3 : Fin 4)) 1
  validA : heap.ValidU64Slice A lenA
  validB : heap.ValidU64Slice B lenB
  Aa : U64SlicesDisjoint A lenA a lenA
  Bb : U64SlicesDisjoint B lenB b lenB
  Ab : U64SlicesDisjoint A lenA b lenB
  BA : U64SlicesDisjoint B lenB A lenA
  row0a : U64SlicesDisjoint (hgcdMatPtr M hM (0 : Fin 4)) 1 a lenA
  row3a : U64SlicesDisjoint (hgcdMatPtr M hM (3 : Fin 4)) 1 a lenA
  row0b : U64SlicesDisjoint (hgcdMatPtr M hM (0 : Fin 4)) 1 b lenB
  row3b : U64SlicesDisjoint (hgcdMatPtr M hM (3 : Fin 4)) 1 b lenB
  aMatrix : ∀ i : Fin 4, U64SlicesDisjoint A lenA
    (hgcdMatPtr M hM i) (identityEntryLen i)
  bMatrix : ∀ i : Fin 4, U64SlicesDisjoint B lenB
    (hgcdMatPtr M hM i) (identityEntryLen i)
  matrixValid : ∀ i : Fin 4, heap.ValidU64Slice
    (hgcdMatPtr M hM i) (identityEntryLen i)
  leftRep : RawDensePolyRep this heap a lenA left
  rightRep : RawDensePolyRep this heap b lenB right

/-- The base path of the complete strictly decreasing body establishes the
common recursive raw invariant from its actual returned record. -/
theorem hgcdRecursiveBodyBelow_base_rawInvariant (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (bound : Nat) (recurse : HgcdRecursiveCallBelow bound)
    (M : HgcdMat) (hM : M.Valid) (computeM : Bool)
    (A B a b : RawPtr UInt64) (lenA lenB : Nat)
    (W scratch : RawPtr UInt64) (heap : RawHeap)
    (left right : Polynomial (ZMod this._p.toNat))
    (hp : 1 < this._p.toNat)
    (hbound : lenA = bound) (horder : lenB < lenA)
    (workspace : HgcdRecursiveBaseCallWorkspace this M hM A B a b lenA
      lenB heap left right)
    (hlenAPos : 0 < lenA) (hstop : lenB < lenA / 2 + 1)
    (result : HgcdRecursiveResult)
    (hrun : hgcdRecursiveBodyBelow this bound recurse M hM computeM A B a b
      lenA lenB W scratch heap hbound horder
        (fun hnonbase => False.elim (hnonbase hstop))
        (fun hnonbase => False.elim (hnonbase hstop)) =
        .ok result) :
    ∃ entries,
      HgcdRecursiveRawInvariant this left right left right entries computeM
        A B lenA result := by
  rw [hgcdRecursiveBodyBelow] at hrun
  simp only [hstop, ↓reduceDIte] at hrun
  split at hrun
  next fault hbase => simp at hrun
  next base hbase =>
    rcases hgcdRecursiveBase_rawInvariant this M computeM A B a b lenA lenB
        heap base left right hM hp workspace.valid0 workspace.valid3
        workspace.disjoint03 workspace.validA workspace.validB workspace.Aa
        workspace.Bb workspace.Ab workspace.BA workspace.row0a workspace.row3a
        workspace.row0b workspace.row3b workspace.aMatrix workspace.bMatrix
        workspace.matrixValid workspace.leftRep workspace.rightRep
        (Nat.le_of_lt horder) hlenAPos hstop hbase with
      ⟨entries, hBaseValid, hInvariant⟩
    have heq : base.toResult hBaseValid = result := by
      apply HgcdRecursiveResult.ext_value
      simpa only [HgcdRecursiveResult.value, HgcdRecursiveBaseResult.toResult]
        using congrArg HgcdRecursiveResult.value (Except.ok.inj hrun)
    subst result
    exact ⟨entries, hInvariant⟩

/-- The non-base early arm of the well-founded body carries the common raw
invariant.  The first-child premise is the induction result for the actual
dispatch execution; all remaining premises describe the actual source split,
four-call reconstruction, and output-copy workspace. -/
theorem hgcdRecursiveBodyBelow_early_rawInvariant (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (bound : Nat) (recurse : HgcdRecursiveCallBelow bound)
    (M : HgcdMat) (hM : M.Valid) (computeM : Bool)
    (A B a b : RawPtr UInt64) (lenA lenB : Nat)
    (W scratch : RawPtr UInt64) (heap : RawHeap)
    (left right inputHighA inputHighB lowPolyA lowPolyB :
      Polynomial (ZMod this._p.toNat))
    (hp : 1 < this._p.toNat) (hcfg : DensePreinvConfigured this)
    (hbound : lenA = bound) (horder : lenB < lenA)
    (firstCall : HgcdRecursiveFirstCallAdmissible this bound recurse a b W
      scratch lenA lenB heap inputHighA inputHighB lowPolyA lowPolyB)
    (first : HgcdRecursiveResult)
    (reconstructed : HgcdRecursiveReconstructPairResult)
    (early : HgcdRecursiveEarlyResult)
    (hbase : ¬ lenB < lenA / 2 + 1)
    (hfirst :
      let ws := hgcdRecursiveWorkspace W lenA
      let high := hgcdRecursiveHighInput a b lenA lenB
      hgcdRecursiveDispatchBelow this bound recurse ws.R
        (hgcdRecursiveWorkspace_R_valid W lenA) ws.a3 ws.b3 high.a0 high.b0
        high.lenA0 high.lenB0 ws.q ws.W3 ws.T0 ws.T1 scratch ws.a2 ws.next
        heap (hgcdRecursiveHighInput_order a b lenA lenB horder)
          (by
            rw [← hbound]
            exact hgcdRecursiveHighInput_len_lt a b lenA lenB horder
              (by omega)) = .ok first)
    (hFullA : left = lowPolyA + Polynomial.X ^ (lenA / 2) * inputHighA)
    (hFullB : right = lowPolyB + Polynomial.X ^ (lenA / 2) * inputHighB)
    (reconstructWork :
      let ws := hgcdRecursiveWorkspace W lenA
      HgcdRecursiveReconstructPairWorkspaceProvider this ws.a2 ws.b2 ws.T0
        a b ws.a3 ws.b3 scratch (Nat.min lenA (lenA / 2))
        (Nat.min lenB (lenA / 2)) first.lenA first.lenB (lenA / 2)
        first.matrix first.valid first.sgn first.heap)
    (hreconstruct :
      let ws := hgcdRecursiveWorkspace W lenA
      hgcdRecursiveReconstructPair this ws.a2 ws.b2 ws.T0 a b ws.a3 ws.b3
        scratch (Nat.min lenA (lenA / 2)) (Nat.min lenB (lenA / 2))
        first.lenA first.lenB (lenA / 2) first.matrix first.valid first.sgn
        first.heap = .ok reconstructed)
    (hearlyGuard : reconstructed.lenB < lenA / 2 + 1)
    (earlyWork :
      let ws := hgcdRecursiveWorkspace W lenA
      HgcdEarlyReturnRefineWorkspace reconstructed.heap M first.matrix hM
        first.valid A B ws.a2 ws.b2 reconstructed.lenA reconstructed.lenB)
    (hearly :
      let ws := hgcdRecursiveWorkspace W lenA
      hgcdRecursiveEarlyReturn M first.matrix hM first.valid computeM A B
        ws.a2 ws.b2 reconstructed.lenA reconstructed.lenB first.sgn
        reconstructed.heap = .ok early)
    (result : HgcdRecursiveResult)
    (hrun :
      let providers := hgcdRecursiveFirstCall_providers this bound recurse a b
        W scratch lenA lenB heap inputHighA inputHighB lowPolyA lowPolyB hcfg
        hp horder firstCall.workspace firstCall.recursiveRefines
      hgcdRecursiveBodyBelow this bound recurse M hM computeM A B a b lenA
        lenB W scratch heap hbound horder (fun _ => providers.1)
          (fun _ => providers.2) = .ok result) :
    ∃ finalA finalB entries,
      HgcdRecursiveRawInvariant this left right finalA finalB entries computeM
        A B lenA result := by
  let ws := hgcdRecursiveWorkspace W lenA
  let high := hgcdRecursiveHighInput a b lenA lenB
  let providers := hgcdRecursiveFirstCall_providers this bound recurse a b W
    scratch lenA lenB heap inputHighA inputHighB lowPolyA lowPolyB hcfg hp
    horder firstCall.workspace firstCall.recursiveRefines
  change hgcdRecursiveBodyBelow this bound recurse M hM computeM A B a b lenA
    lenB W scratch heap hbound horder (fun _ => providers.1)
      (fun _ => providers.2) = .ok result
    at hrun
  rcases hgcdRecursiveFirstCall_refines this bound recurse a b W scratch lenA
      lenB heap inputHighA inputHighB lowPolyA lowPolyB hcfg hp firstCall first
      (hgcdRecursiveHighInput_order a b lenA lenB horder)
      (by
        rw [← hbound]
        exact hgcdRecursiveHighInput_len_lt a b lenA lenB horder (by omega))
      hfirst with
    ⟨outputHighA, outputHighB, entries, hFirst, hLowA, hLowB⟩
  have hfirstLength : HgcdRecursiveLengthInvariant high.lenA0 first :=
    providers.1 first (hgcdRecursiveHighInput_order a b lenA lenB horder)
      (by
        rw [← hbound]
        exact hgcdRecursiveHighInput_len_lt a b lenA lenB horder (by omega))
      hfirst
  have hReconstructed : HgcdFirstReconstructionInvariant lenA first
      reconstructed :=
    providers.2 first reconstructed ⟨_, _, hfirst⟩ (by
      simpa [high, hgcdRecursiveHighInput] using hfirstLength) hreconstruct
  let finalA := hgcdReconstructedLowA entries lowPolyA lowPolyB first.sgn +
    Polynomial.X ^ (lenA / 2) * outputHighA
  let finalB := hgcdReconstructedLowB entries lowPolyA lowPolyB first.sgn +
    Polynomial.X ^ (lenA / 2) * outputHighB
  rcases hgcdRecursiveReconstructPair_preserves_input this ws.a2 ws.b2 ws.T0
      a b ws.a3 ws.b3 scratch (Nat.min lenA (lenA / 2))
      (Nat.min lenB (lenA / 2)) (lenA / 2) high.lenA0 first reconstructed
      entries left right lowPolyA lowPolyB inputHighA inputHighB outputHighA
      outputHighB hcfg hp reconstructWork hFirst hLowA.toCanonicalSlice
      hLowB.toCanonicalSlice hFullA hFullB
      hreconstruct with
    ⟨hARep, hBRep, hTransform, hDet, hGcd, _⟩
  have hLength : HgcdRecursiveLengthInvariant lenA
      ⟨reconstructed.heap, first.matrix, first.valid, reconstructed.lenA,
        reconstructed.lenB, first.sgn⟩ :=
    hgcdRecursiveEarly_lengthInvariant lenA high.lenA0 (lenA / 2) first
      reconstructed reconstructed.heap first.sgn rfl (by
        simp [high, hgcdRecursiveHighInput]) hfirstLength hReconstructed
      hearlyGuard
  have hMatrixAtReconstructed : HgcdMatRawDenseRep this reconstructed.heap
      first.matrix entries first.valid :=
    hgcdRecursiveReconstructPair_preserves_matrix this ws.a2 ws.b2 ws.T0 a
      b ws.a3 ws.b3 scratch (Nat.min lenA (lenA / 2))
      (Nat.min lenB (lenA / 2)) first.lenA first.lenB (lenA / 2)
      first.matrix first.valid first.sgn first.heap reconstructed entries
      reconstructWork (hFirst.matrixSemantics rfl).1 hreconstruct
  rcases hgcdRecursiveEarlyReturn_rawInvariant this M first.matrix hM
      first.valid computeM A B ws.a2 ws.b2 reconstructed.lenA
      reconstructed.lenB lenA first.sgn left right finalA finalB entries
      reconstructed.heap early earlyWork hARep hBRep
      hMatrixAtReconstructed hTransform hDet hGcd hearlyGuard hLength
      hearly with ⟨hEarlyValid, hInvariant⟩
  rw [hgcdRecursiveBodyBelow] at hrun
  simp only [hbase, ↓reduceDIte] at hrun
  split at hrun
  next fault hdispatch =>
    simp at hrun
  next actualFirst hdispatch =>
    have hdispatchOk : (.ok actualFirst : RawExec HgcdRecursiveResult) =
        .ok first := hdispatch.symm.trans (by
          convert hfirst using 1)
    have hFirstEq : actualFirst = first :=
      Except.ok.inj hdispatchOk
    subst actualFirst
    split at hrun
    next fault hreconstructActual =>
      simp at hrun
    next actualReconstructed hreconstructActual =>
      have hreconstructOk :
          (.ok actualReconstructed : RawExec
            HgcdRecursiveReconstructPairResult) = .ok reconstructed :=
        hreconstructActual.symm.trans (by
          convert hreconstruct using 1)
      have hReconstructedEq : actualReconstructed = reconstructed :=
        Except.ok.inj hreconstructOk
      subst actualReconstructed
      simp only [hearlyGuard, ↓reduceDIte] at hrun
      split at hrun
      next fault hearlyActual =>
        simp at hrun
      next actualEarly hearlyActual =>
        have hearlyOk : (.ok actualEarly : RawExec HgcdRecursiveEarlyResult) =
            .ok early := hearlyActual.symm.trans (by
              convert hearly using 1)
        have hEarlyEq : actualEarly = early :=
          Except.ok.inj hearlyOk
        subst actualEarly
        have heq : early.toResult hEarlyValid = result := by
          apply HgcdRecursiveResult.ext_value
          simpa only [HgcdRecursiveResult.value,
              HgcdRecursiveEarlyResult.toResult] using
            congrArg HgcdRecursiveResult.value (Except.ok.inj hrun)
        subst result
        exact ⟨finalA, finalB, entries, hInvariant⟩

/-- Physical division workspace used between the two recursive HGCD calls.
Every field describes an allocation or non-aliasing fact consumed by the
actual generated `hgcdRecursiveMiddle` execution. -/
structure HgcdRecursiveMiddleWorkspace (W : RawPtr UInt64) (lenA : Nat)
    (reconstructed : HgcdRecursiveReconstructPairResult) : Prop where
  validQ :
    let ws := hgcdRecursiveWorkspace W lenA
    reconstructed.heap.ValidU64Slice ws.q
      (reconstructed.lenA - (reconstructed.lenB - 1))
  validD :
    let ws := hgcdRecursiveWorkspace W lenA
    reconstructed.heap.ValidU64Slice ws.d
      (Nat.min reconstructed.lenA (reconstructed.lenB - 1))
  validW3 :
    let ws := hgcdRecursiveWorkspace W lenA
    reconstructed.heap.ValidWord3Slice ws.W3 reconstructed.lenA
  quotientCapacity : reconstructed.lenA - (reconstructed.lenB - 1) < limbBase
  dA :
    let ws := hgcdRecursiveWorkspace W lenA
    ws.d.region ≠ ws.a2.region
  wA :
    let ws := hgcdRecursiveWorkspace W lenA
    ws.W3.region ≠ ws.a2.region
  wB :
    let ws := hgcdRecursiveWorkspace W lenA
    ws.W3.region ≠ ws.b2.region
  qB :
    let ws := hgcdRecursiveWorkspace W lenA
    ws.q.region ≠ ws.b2.region
  qW :
    let ws := hgcdRecursiveWorkspace W lenA
    ws.q.region ≠ ws.W3.region
  dW :
    let ws := hgcdRecursiveWorkspace W lenA
    ws.d.region ≠ ws.W3.region
  dQ :
    let ws := hgcdRecursiveWorkspace W lenA
    ws.d.region ≠ ws.q.region
  dB :
    let ws := hgcdRecursiveWorkspace W lenA
    ws.d.region ≠ ws.b2.region

/-- Physical data tied to the actual second recursive call and the following
finish.  It contains no claim about the recursive callback's polynomial
semantics; that claim is supplied separately by well-founded induction. -/
structure HgcdRecursiveSecondCallWorkspace (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (bound : Nat) (recurse : HgcdRecursiveCallBelow bound)
    (M : HgcdMat) (hM : M.Valid) (A B : RawPtr UInt64)
    (W scratch : RawPtr UInt64) (lenA : Nat)
    (first : HgcdRecursiveResult)
    (reconstructed : HgcdRecursiveReconstructPairResult)
    (middle : HgcdRecursiveMiddleResult) (second : HgcdRecursiveResult)
    (hsecondOrder : middle.lenD0 < middle.lenC0)
    (hsecondDecrease : middle.lenC0 < bound) : Prop where
  matrixPrefix :
    let ws := hgcdRecursiveWorkspace W lenA
    ∀ i : Fin 4, SameU64Prefix reconstructed.heap middle.heap
      (hgcdMatPtr first.matrix first.valid i)
      (hgcdMatLen first.matrix first.valid i)
  iterWorkspace : ∀ highC highD,
    RawDensePolyRep this middle.heap middle.c0 middle.lenC0 highC →
    RawDensePolyRep this middle.heap middle.d0 middle.lenD0 highD →
    let ws := hgcdRecursiveWorkspace W lenA
    HgcdRecursiveDispatchIterWorkspace this ws.S
      (hgcdRecursiveWorkspace_S_valid W lenA) ws.a3 ws.b3 middle.c0 middle.d0
      middle.lenC0 middle.lenD0 ws.a2 ws.W3 ws.T0 ws.T1 scratch ws.a2
      middle.heap highC highD
  frame :
    let ws := hgcdRecursiveWorkspace W lenA
    HgcdRecursiveSecondDispatchFrameProvider this bound recurse ws.S
      (hgcdRecursiveWorkspace_S_valid W lenA) ws.a3 ws.b3 middle.c0 middle.d0
      middle.lenC0 middle.lenD0 ws.a2 ws.W3 ws.T0 ws.T1 scratch ws.a2 ws.next
      middle.heap hsecondOrder hsecondDecrease first.matrix first.valid ws.q
      ws.b2 ws.d middle.lenQ (Nat.min reconstructed.lenB middle.k)
      (Nat.min middle.lenD middle.k)
  finish :
    let ws := hgcdRecursiveWorkspace W lenA
    HgcdRecursiveFinishWorkspaceProvider this M first.matrix second.matrix hM
      first.valid second.valid A B ws.T0 ws.b2 ws.d ws.a3 ws.b3 ws.q
      (Nat.min reconstructed.lenB middle.k) (Nat.min middle.lenD middle.k)
      second.lenA second.lenB middle.k middle.lenQ ws.a2 scratch second.sgn
      second.heap

/-- The semantic fact excluded from the second-call physical workspace.  It is
obtained from well-founded induction at `middle.lenC0 < bound`. -/
def HgcdRecursiveSecondCallbackRefines (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (bound : Nat) (recurse : HgcdRecursiveCallBelow bound)
    (W scratch : RawPtr UInt64) (lenA : Nat)
    (middle : HgcdRecursiveMiddleResult)
    (hsecondOrder : middle.lenD0 < middle.lenC0)
    (hsecondDecrease : middle.lenC0 < bound) : Prop :=
  ∀ highC highD,
    RawDensePolyRep this middle.heap middle.c0 middle.lenC0 highC →
    RawDensePolyRep this middle.heap middle.d0 middle.lenD0 highD →
    let ws := hgcdRecursiveWorkspace W lenA
    HgcdRecursiveCallbackRefinesAt this bound recurse ws.S
      (hgcdRecursiveWorkspace_S_valid W lenA) ws.a3 ws.b3 middle.c0 middle.d0
      middle.lenC0 middle.lenD0 ws.next scratch middle.heap hsecondOrder
      hsecondDecrease highC highD

/-- All data tied to the actual second recursive call and the following
finish.  Its semantic field is exactly the smaller-call induction hypothesis;
all remaining fields are bundled in `workspace`. -/
structure HgcdRecursiveSecondCallAdmissible (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (bound : Nat) (recurse : HgcdRecursiveCallBelow bound)
    (M : HgcdMat) (hM : M.Valid) (A B : RawPtr UInt64)
    (W scratch : RawPtr UInt64) (lenA : Nat)
    (first : HgcdRecursiveResult)
    (reconstructed : HgcdRecursiveReconstructPairResult)
    (middle : HgcdRecursiveMiddleResult) (second : HgcdRecursiveResult)
    (hsecondOrder : middle.lenD0 < middle.lenC0)
    (hsecondDecrease : middle.lenC0 < bound) : Prop where
  workspace : HgcdRecursiveSecondCallWorkspace this bound recurse M hM A B W
    scratch lenA first reconstructed middle second hsecondOrder
    hsecondDecrease
  recursiveRefines : HgcdRecursiveSecondCallbackRefines this bound recurse W
    scratch lenA middle hsecondOrder hsecondDecrease

/-- Attach only the well-founded child semantic theorem to an already
established physical second-call workspace. -/
theorem HgcdRecursiveSecondCallWorkspace.admissible (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (bound : Nat) (recurse : HgcdRecursiveCallBelow bound)
    (M : HgcdMat) (hM : M.Valid) (A B : RawPtr UInt64)
    (W scratch : RawPtr UInt64) (lenA : Nat)
    (first : HgcdRecursiveResult)
    (reconstructed : HgcdRecursiveReconstructPairResult)
    (middle : HgcdRecursiveMiddleResult) (second : HgcdRecursiveResult)
    (hsecondOrder : middle.lenD0 < middle.lenC0)
    (hsecondDecrease : middle.lenC0 < bound)
    (workspace : HgcdRecursiveSecondCallWorkspace this bound recurse M hM A B
      W scratch lenA first reconstructed middle second hsecondOrder
      hsecondDecrease)
    (recursiveRefines : HgcdRecursiveSecondCallbackRefines this bound recurse W
      scratch lenA middle hsecondOrder hsecondDecrease) :
    HgcdRecursiveSecondCallAdmissible this bound recurse M hM A B W scratch
      lenA first reconstructed middle second hsecondOrder hsecondDecrease :=
  ⟨workspace, recursiveRefines⟩

/-- The complete non-early arm of the well-founded body, following the
generated source from the first dispatch through reconstruction, middle
divrem, second dispatch, and finish.  All semantic values are extracted from
those successful executions; the extra premises are physical workspace and
frame obligations only. -/
theorem hgcdRecursiveBodyBelow_nonEarly_rawInvariant (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (bound : Nat) (recurse : HgcdRecursiveCallBelow bound)
    (M : HgcdMat) (hM : M.Valid) (computeM : Bool)
    (A B a b : RawPtr UInt64) (lenA lenB : Nat)
    (W scratch : RawPtr UInt64) (heap : RawHeap)
    (left right inputHighA inputHighB lowPolyA lowPolyB :
      Polynomial (ZMod this._p.toNat))
    (hp : 1 < this._p.toNat) (hcfg : DensePreinvConfigured this)
    (hbound : lenA = bound) (horder : lenB < lenA)
    (firstCall : HgcdRecursiveFirstCallAdmissible this bound recurse a b W
      scratch lenA lenB heap inputHighA inputHighB lowPolyA lowPolyB)
    (first : HgcdRecursiveResult)
    (reconstructed : HgcdRecursiveReconstructPairResult)
    (middle : HgcdRecursiveMiddleResult)
    (second : HgcdRecursiveResult)
    (finished : HgcdRecursiveFinishResult)
    (hbase : ¬ lenB < lenA / 2 + 1)
    (hfirst :
      let ws := hgcdRecursiveWorkspace W lenA
      let high := hgcdRecursiveHighInput a b lenA lenB
      hgcdRecursiveDispatchBelow this bound recurse ws.R
        (hgcdRecursiveWorkspace_R_valid W lenA) ws.a3 ws.b3 high.a0 high.b0
        high.lenA0 high.lenB0 ws.q ws.W3 ws.T0 ws.T1 scratch ws.a2 ws.next
        heap (hgcdRecursiveHighInput_order a b lenA lenB horder)
          (by
            rw [← hbound]
            exact hgcdRecursiveHighInput_len_lt a b lenA lenB horder
              (by omega)) = .ok first)
    (hFullA : left = lowPolyA + Polynomial.X ^ (lenA / 2) * inputHighA)
    (hFullB : right = lowPolyB + Polynomial.X ^ (lenA / 2) * inputHighB)
    (reconstructWork :
      let ws := hgcdRecursiveWorkspace W lenA
      HgcdRecursiveReconstructPairWorkspaceProvider this ws.a2 ws.b2 ws.T0
        a b ws.a3 ws.b3 scratch (Nat.min lenA (lenA / 2))
        (Nat.min lenB (lenA / 2)) first.lenA first.lenB (lenA / 2)
        first.matrix first.valid first.sgn first.heap)
    (hreconstruct :
      let ws := hgcdRecursiveWorkspace W lenA
      hgcdRecursiveReconstructPair this ws.a2 ws.b2 ws.T0 a b ws.a3 ws.b3
        scratch (Nat.min lenA (lenA / 2)) (Nat.min lenB (lenA / 2))
        first.lenA first.lenB (lenA / 2) first.matrix first.valid first.sgn
        first.heap = .ok reconstructed)
    (hearly : ¬ reconstructed.lenB < lenA / 2 + 1)
    (middleWork : HgcdRecursiveMiddleWorkspace W lenA reconstructed)
    (hmiddle :
      let ws := hgcdRecursiveWorkspace W lenA
      hgcdRecursiveMiddle this ws.q ws.d ws.a2 ws.b2 reconstructed.lenA
        reconstructed.lenB (lenA / 2) ws.W3 reconstructed.heap = .ok middle)
    (hsecondOrder : middle.lenD0 < middle.lenC0)
    (hsecondDecrease : middle.lenC0 < bound)
    (secondWork : HgcdRecursiveSecondCallWorkspace this bound recurse M hM A B
      W scratch lenA first reconstructed middle second hsecondOrder
      hsecondDecrease)
    (secondRecursiveRefines : HgcdRecursiveSecondCallbackRefines this bound
      recurse W scratch lenA middle hsecondOrder hsecondDecrease)
    (hsecond :
      let ws := hgcdRecursiveWorkspace W lenA
      hgcdRecursiveDispatchBelow this bound recurse ws.S
        (hgcdRecursiveWorkspace_S_valid W lenA) ws.a3 ws.b3 middle.c0
        middle.d0 middle.lenC0 middle.lenD0 ws.a2 ws.W3 ws.T0 ws.T1 scratch
        ws.a2 ws.next middle.heap hsecondOrder hsecondDecrease = .ok second)
    (hfinish :
      let ws := hgcdRecursiveWorkspace W lenA
      hgcdRecursiveFinish this M first.matrix second.matrix hM first.valid
        second.valid computeM A B ws.T0 ws.b2 ws.d ws.a3 ws.b3 ws.q
        (Nat.min reconstructed.lenB middle.k) (Nat.min middle.lenD middle.k)
        second.lenA second.lenB middle.k middle.lenQ ws.a2 scratch first.sgn
        second.sgn second.heap = .ok finished)
    (result : HgcdRecursiveResult)
    (hrun :
      let providers := hgcdRecursiveFirstCall_providers this bound recurse a b
        W scratch lenA lenB heap inputHighA inputHighB lowPolyA lowPolyB hcfg
        hp horder firstCall.workspace firstCall.recursiveRefines
      hgcdRecursiveBodyBelow this bound recurse M hM computeM A B a b lenA
        lenB W scratch heap hbound horder (fun _ => providers.1)
          (fun _ => providers.2) = .ok result) :
    ∃ finalA finalB entries,
      HgcdRecursiveRawInvariant this left right finalA finalB entries computeM
        A B lenA result := by
  let ws := hgcdRecursiveWorkspace W lenA
  let high := hgcdRecursiveHighInput a b lenA lenB
  let secondCall := secondWork.admissible this bound recurse M hM A B W scratch
    lenA first reconstructed middle second hsecondOrder hsecondDecrease
    secondRecursiveRefines
  let providers := hgcdRecursiveFirstCall_providers this bound recurse a b W
    scratch lenA lenB heap inputHighA inputHighB lowPolyA lowPolyB hcfg hp
    horder firstCall.workspace firstCall.recursiveRefines
  change hgcdRecursiveBodyBelow this bound recurse M hM computeM A B a b lenA
    lenB W scratch heap hbound horder (fun _ => providers.1)
      (fun _ => providers.2) = .ok result
    at hrun
  rcases hgcdRecursiveFirstCall_refines this bound recurse a b W scratch lenA
      lenB heap inputHighA inputHighB lowPolyA lowPolyB hcfg hp firstCall first
      (hgcdRecursiveHighInput_order a b lenA lenB horder)
      (by
        rw [← hbound]
        exact hgcdRecursiveHighInput_len_lt a b lenA lenB horder (by omega))
      hfirst with
    ⟨outputHighA, outputHighB, firstEntries, hFirst, hLowA, hLowB⟩
  have hfirstLength : HgcdRecursiveLengthInvariant high.lenA0 first :=
    providers.1 first (hgcdRecursiveHighInput_order a b lenA lenB horder)
      (by
        rw [← hbound]
        exact hgcdRecursiveHighInput_len_lt a b lenA lenB horder (by omega))
      hfirst
  have hReconstructed : HgcdFirstReconstructionInvariant lenA first
      reconstructed := providers.2 first reconstructed ⟨_, _, hfirst⟩ (by
    simpa [high, hgcdRecursiveHighInput] using hfirstLength) hreconstruct
  let currentA := hgcdReconstructedLowA firstEntries lowPolyA lowPolyB
      first.sgn + Polynomial.X ^ (lenA / 2) * outputHighA
  let currentB := hgcdReconstructedLowB firstEntries lowPolyA lowPolyB
      first.sgn + Polynomial.X ^ (lenA / 2) * outputHighB
  have hFirstReconstruct := hgcdRecursiveReconstructPair_preserves_input this
    ws.a2 ws.b2 ws.T0 a b ws.a3 ws.b3 scratch
    (Nat.min lenA (lenA / 2)) (Nat.min lenB (lenA / 2)) (lenA / 2)
    high.lenA0 first reconstructed firstEntries left right lowPolyA lowPolyB
    inputHighA inputHighB outputHighA outputHighB hcfg hp reconstructWork
    hFirst hLowA.toCanonicalSlice hLowB.toCanonicalSlice hFullA hFullB
    hreconstruct
  have hFirstMatrixAtReconstructed : HgcdMatRawDenseRep this
      reconstructed.heap first.matrix firstEntries first.valid :=
    hgcdRecursiveReconstructPair_preserves_matrix this ws.a2 ws.b2 ws.T0 a b
      ws.a3 ws.b3 scratch (Nat.min lenA (lenA / 2))
      (Nat.min lenB (lenA / 2)) first.lenA first.lenB (lenA / 2)
      first.matrix first.valid first.sgn first.heap reconstructed firstEntries
      reconstructWork (hFirst.matrixSemantics rfl).1 hreconstruct
  rcases hgcdRecursiveMiddle_refines this ws.q ws.d ws.a2 ws.b2
      reconstructed.lenA reconstructed.lenB (lenA / 2) lenA ws.W3
      reconstructed.heap currentA currentB (by omega)
      (by simpa [currentA] using hFirstReconstruct.1)
      (by simpa [currentB] using hFirstReconstruct.2.1)
      middleWork.validQ middleWork.validD middleWork.validW3
      middleWork.quotientCapacity middleWork.dA middleWork.wA middleWork.wB
      middleWork.qB middleWork.qW middleWork.dW middleWork.dQ middleWork.dB
      hcfg (Fact.out : Nat.Prime this._p.toNat) (by omega)
      (Nat.le_of_lt hReconstructed.decreases) (by omega) (by omega) with
    ⟨actualMiddle, quotient, remainder, hmiddle', _, _, _, _, _, _,
      hMiddleLayout, hQRep, hDRep, hBRep, hDivision, _, hlenQ, hlenDCap, _, hlenC,
      hlenCPos, hlenD0, hk, hcPtr, hcLen, hdPtr, hdLen⟩
  have hMiddleEq : actualMiddle = middle :=
    Except.ok.inj (hmiddle'.symm.trans hmiddle)
  subst actualMiddle
  have hMiddleSourceLayout := hgcdRecursiveMiddle_layout this ws.q ws.d
    ws.a2 ws.b2 reconstructed.lenA reconstructed.lenB (lenA / 2) ws.W3
    reconstructed.heap middle hmiddle
  have hkLtReconstructed : middle.k < reconstructed.lenB := by
    rw [hMiddleSourceLayout.2.2.2.1] at hlenCPos
    split at hlenCPos <;> omega
  have hkCapacity : middle.k ≤
      Nat.min reconstructed.lenA (reconstructed.lenB - 1) := by
    apply Nat.le_min.mpr
    constructor
    · exact (Nat.le_of_lt hkLtReconstructed).trans hReconstructed.ordered
    · omega
  have hD0Valid : middle.heap.ValidU64Slice middle.d0 middle.lenD0 := by
    rw [hMiddleSourceLayout.2.2.2.2.1,
      hMiddleSourceLayout.2.2.2.2.2]
    apply middle.heap.validU64Slice_add ws.d
      (Nat.min reconstructed.lenA (reconstructed.lenB - 1)) middle.k
      (if middle.lenD ≥ middle.k then middle.lenD - middle.k else 0)
      ((hMiddleLayout ws.d
        (Nat.min reconstructed.lenA (reconstructed.lenB - 1))).mp
          middleWork.validD)
    split <;> omega
  rcases hgcdRecursiveMiddle_split_reps this ws.q ws.d ws.a2 ws.b2
      reconstructed.lenA reconstructed.lenB (lenA / 2) ws.W3
      reconstructed.heap middle currentB remainder hBRep hDRep hD0Valid
      hlenCPos hmiddle with
    ⟨lowC, highC, lowD, highD, hLowC, hHighC, hLowD, hHighD, hSplitC,
      hSplitD⟩
  have hFirstMatrixAtMiddle : HgcdMatRawDenseRep this middle.heap first.matrix
      firstEntries first.valid := by
    intro i
    exact rawDensePolyRep_of_same_prefix this reconstructed.heap middle.heap
      (hgcdMatPtr first.matrix first.valid i)
      (hgcdMatLen first.matrix first.valid i) (firstEntries i) hMiddleLayout
      (secondCall.workspace.matrixPrefix i) (hFirstMatrixAtReconstructed i)
  rcases hgcdRecursiveSecondDispatch_refines this bound recurse ws.S
      (hgcdRecursiveWorkspace_S_valid W lenA) ws.a3 ws.b3 middle.c0 middle.d0
      middle.lenC0 middle.lenD0 ws.a2 ws.W3 ws.T0 ws.T1 scratch ws.a2
      ws.next middle.heap highC highD first.matrix first.valid ws.q ws.b2 ws.d
      middle.lenQ (Nat.min reconstructed.lenB middle.k)
      (Nat.min middle.lenD middle.k) firstEntries quotient lowC lowD hcfg hp
      hsecondOrder hsecondDecrease
      (secondCall.workspace.iterWorkspace highC highD hHighC hHighD)
      (secondCall.recursiveRefines highC highD hHighC hHighD)
      secondCall.workspace.frame
      hFirstMatrixAtMiddle hQRep hLowC hLowD second hsecond with
    ⟨finalHighA, finalHighB, secondEntries, hSecond, hRAfter, hQAfter,
      hLowCAfter, hLowDAfter⟩
  have hSecondLength : HgcdRecursiveLengthInvariant middle.lenC0 second :=
    hSecond.lengths rfl
  have hcExact : middle.lenC0 = reconstructed.lenB - middle.k := by
    rw [hcLen]
    simp [Nat.le_of_lt hkLtReconstructed]
  rcases hgcdRecursiveFinish_exec this M first.matrix second.matrix hM
      first.valid second.valid computeM A B ws.T0 ws.b2 ws.d ws.a3 ws.b3
      ws.q (Nat.min reconstructed.lenB middle.k)
      (Nat.min middle.lenD middle.k) second.lenA second.lenB middle.k
      middle.lenQ ws.a2 scratch first.sgn second.sgn second.heap finished
      hfinish with ⟨finishReconstructed, hfinishReconstruct, _, _, _, _⟩
  have hFinishWork := secondCall.workspace.finish finishReconstructed
    hfinishReconstruct
  have hFinishOperands := hgcdRecursiveFinish_operandInvariant this M
    first.matrix hM first.valid computeM A B ws.T0 ws.b2 ws.d ws.a3 ws.b3
    ws.q lenA (lenA / 2) reconstructed.lenB middle.lenD middle.k
    middle.lenC0 middle.lenQ ws.a2 scratch first.sgn second finished
    secondEntries lowC lowD finalHighA finalHighB hcfg hp hSecondLength
    hFinishWork.reconstruct (hSecond.matrixSemantics rfl).1 hLowCAfter
    hLowDAfter hSecond.aRep hSecond.bRep rfl hk hcExact (by omega)
    hReconstructed.decreases hfinish
  have hFinishStop : finished.lenB < lenA / 2 + 1 :=
    hFinishOperands.stopped
  have hFinishLength : computeM = true →
      HgcdRecursiveLengthInvariant lenA finished.toResult := by
    intro hcompute
    subst computeM
    apply hgcdRecursiveFinish_lengthInvariant this M hM A B ws.T0 ws.b2
      ws.d ws.a3 ws.b3 ws.q ws.a2 scratch lenA high.lenA0 (lenA / 2)
      reconstructed.lenA reconstructed.lenB middle.lenD middle.k
      middle.lenC0 middle.lenQ first second finished firstEntries secondEntries
      quotient lowC lowD finalHighA finalHighB hcfg hp hfirstLength
      hSecondLength secondCall.workspace.finish hRAfter
      (hSecond.matrixSemantics rfl).1
      hQAfter hLowCAfter hLowDAfter hSecond.aRep hSecond.bRep rfl
      (by simp [high, hgcdRecursiveHighInput]) hReconstructed.leadingA
      hReconstructed.ordered (by omega) hReconstructed.decreases hlenQ hk
      hcExact hfinish
  have hTail := hgcdRecursiveSecondDispatchFinish_rawInvariant this bound
    recurse M first.matrix ws.S hM first.valid
    (hgcdRecursiveWorkspace_S_valid W lenA) computeM A B ws.T0 ws.b2 ws.d
    ws.a3 ws.b3 ws.q middle.c0 middle.d0
    (Nat.min reconstructed.lenB middle.k) (Nat.min middle.lenD middle.k)
    middle.k middle.lenQ ws.a2 ws.W3 ws.T1 scratch ws.a2 ws.next ws.a2
    middle.lenC0 middle.lenD0 lenA middle.heap second finished left right
    highC highD currentA currentB remainder quotient lowC lowD firstEntries
    first.sgn hcfg hp hsecondOrder hsecondDecrease
    (secondCall.workspace.iterWorkspace highC highD hHighC hHighD)
    (secondCall.recursiveRefines highC highD hHighC hHighD)
    secondCall.workspace.frame
    hFirstMatrixAtMiddle hQRep hLowC hLowD hSplitC hSplitD
    (by simpa [currentA, currentB] using hFirstReconstruct.2.2.1)
    hFirstReconstruct.2.2.2.1 hDivision secondCall.workspace.finish hsecond
    hfinish
    hFinishStop hFinishLength
  rcases hTail with ⟨finalA, finalB, entries, hInvariant⟩
  rw [hgcdRecursiveBodyBelow] at hrun
  simp only [hbase, ↓reduceDIte] at hrun
  split at hrun
  next fault hfirstActual => simp at hrun
  next actualFirst hfirstActual =>
    have hFirstEq : actualFirst = first := Except.ok.inj
      (hfirstActual.symm.trans (by convert hfirst using 1))
    subst actualFirst
    split at hrun
    next fault hreconstructActual => simp at hrun
    next actualReconstructed hreconstructActual =>
      have hReconstructedEq : actualReconstructed = reconstructed :=
        Except.ok.inj (hreconstructActual.symm.trans (by
          convert hreconstruct using 1))
      subst actualReconstructed
      simp only [hearly, ↓reduceDIte] at hrun
      split at hrun
      next fault hmiddleActual => simp at hrun
      next actualMiddle hmiddleActual =>
        have hMiddleEq : actualMiddle = middle := Except.ok.inj
          (hmiddleActual.symm.trans (by convert hmiddle using 1))
        subst actualMiddle
        split at hrun
        next fault hsecondActual => simp at hrun
        next actualSecond hsecondActual =>
          have hSecondEq : actualSecond = second := Except.ok.inj
            (hsecondActual.symm.trans (by convert hsecond using 1))
          subst actualSecond
          split at hrun
          next fault hfinishActual => simp at hrun
          next actualFinished hfinishActual =>
            have hFinishedEq : actualFinished = finished := Except.ok.inj
              (hfinishActual.symm.trans (by convert hfinish using 1))
            subst actualFinished
            have heq : finished.toResult = result := Except.ok.inj hrun
            subst result
            exact ⟨finalA, finalB, entries, hInvariant⟩

/-- Existentially package the mathematical split and physical workspace for
the first child of a genuinely non-base invocation.  The package is never
requested by the generated base branch and contains no recursive semantics. -/
structure HgcdRecursiveNonBasePackage (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (bound : Nat) (recurse : HgcdRecursiveCallBelow bound)
    (a b W scratch : RawPtr UInt64) (lenA lenB : Nat)
    (heap : RawHeap) : Type where
  inputHighA : Polynomial (ZMod this._p.toNat)
  inputHighB : Polynomial (ZMod this._p.toNat)
  lowPolyA : Polynomial (ZMod this._p.toNat)
  lowPolyB : Polynomial (ZMod this._p.toNat)
  inputA : Polynomial (ZMod this._p.toNat)
  inputB : Polynomial (ZMod this._p.toNat)
  splitA : inputA = lowPolyA + Polynomial.X ^ (lenA / 2) * inputHighA
  splitB : inputB = lowPolyB + Polynomial.X ^ (lenA / 2) * inputHighB
  workspace : HgcdRecursiveFirstCallWorkspace this bound recurse a b W
    scratch lenA lenB heap inputHighA inputHighB lowPolyA lowPolyB

/-- The one semantic fact intentionally excluded from a non-base physical
package.  A top-level well-founded construction supplies it from the smaller
`lenA` induction result. -/
def HgcdRecursiveNonBaseCallbackRefines (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (bound : Nat) (recurse : HgcdRecursiveCallBelow bound)
    (a b W scratch : RawPtr UInt64) (lenA lenB : Nat) (heap : RawHeap)
    (package : HgcdRecursiveNonBasePackage this bound recurse a b W scratch
      lenA lenB heap) : Prop :=
  let ws := hgcdRecursiveWorkspace W lenA
  let high := hgcdRecursiveHighInput a b lenA lenB
  ∀ (hchildOrder : high.lenB0 < high.lenA0)
    (hchildDecrease : high.lenA0 < bound),
  HgcdRecursiveCallbackRefinesAt this bound recurse ws.R
    (hgcdRecursiveWorkspace_R_valid W lenA) ws.a3 ws.b3 high.a0 high.b0
    high.lenA0 high.lenB0 ws.next scratch heap hchildOrder hchildDecrease
    package.inputHighA package.inputHighB

/-- A physical non-base package becomes first-call admissible only after the
smaller-call semantic theorem has been supplied. -/
theorem HgcdRecursiveNonBasePackage.admissible (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (bound : Nat) (recurse : HgcdRecursiveCallBelow bound)
    (a b W scratch : RawPtr UInt64) (lenA lenB : Nat) (heap : RawHeap)
    (package : HgcdRecursiveNonBasePackage this bound recurse a b W scratch
      lenA lenB heap)
    (recursiveRefines : HgcdRecursiveNonBaseCallbackRefines this bound recurse
      a b W scratch lenA lenB heap package) :
    HgcdRecursiveFirstCallAdmissible this bound recurse a b W scratch lenA
      lenB heap package.inputHighA package.inputHighB package.lowPolyA
      package.lowPolyB :=
  package.workspace.admissible this bound recurse a b W scratch lenA lenB
    heap package.inputHighA package.inputHighB package.lowPolyA
    package.lowPolyB recursiveRefines

/-- Physical evidence for the successful early continuation of an actual
non-base invocation.  This package records the generated first dispatch,
reconstruction, and early-return executions, but deliberately contains no
recursive semantic theorem. -/
structure HgcdRecursiveEarlyContinuationWorkspace (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (bound : Nat) (recurse : HgcdRecursiveCallBelow bound)
    (M : HgcdMat) (hM : M.Valid) (computeM : Bool)
    (A B a b : RawPtr UInt64)
    (W scratch : RawPtr UInt64) (lenA lenB : Nat) (heap : RawHeap)
    (package : HgcdRecursiveNonBasePackage this bound recurse a b W scratch
      lenA lenB heap)
    (hbound : lenA = bound) (horder : lenB < lenA)
    (hbase : ¬ lenB < lenA / 2 + 1) : Type where
  first : HgcdRecursiveResult
  reconstructed : HgcdRecursiveReconstructPairResult
  early : HgcdRecursiveEarlyResult
  firstExec :
    let ws := hgcdRecursiveWorkspace W lenA
    let high := hgcdRecursiveHighInput a b lenA lenB
    hgcdRecursiveDispatchBelow this bound recurse ws.R
      (hgcdRecursiveWorkspace_R_valid W lenA) ws.a3 ws.b3 high.a0 high.b0
      high.lenA0 high.lenB0 ws.q ws.W3 ws.T0 ws.T1 scratch ws.a2 ws.next
      heap (hgcdRecursiveHighInput_order a b lenA lenB horder)
        (by
          rw [← hbound]
          exact hgcdRecursiveHighInput_len_lt a b lenA lenB horder
            (by omega)) = .ok first
  reconstructWorkspace :
    let ws := hgcdRecursiveWorkspace W lenA
    HgcdRecursiveReconstructPairWorkspaceProvider this ws.a2 ws.b2 ws.T0
      a b ws.a3 ws.b3 scratch (Nat.min lenA (lenA / 2))
      (Nat.min lenB (lenA / 2)) first.lenA first.lenB (lenA / 2)
      first.matrix first.valid first.sgn first.heap
  reconstructExec :
    let ws := hgcdRecursiveWorkspace W lenA
    hgcdRecursiveReconstructPair this ws.a2 ws.b2 ws.T0 a b ws.a3 ws.b3
      scratch (Nat.min lenA (lenA / 2)) (Nat.min lenB (lenA / 2))
      first.lenA first.lenB (lenA / 2) first.matrix first.valid first.sgn
      first.heap = .ok reconstructed
  guard : reconstructed.lenB < lenA / 2 + 1
  earlyWorkspace :
    let ws := hgcdRecursiveWorkspace W lenA
    HgcdEarlyReturnRefineWorkspace reconstructed.heap M first.matrix hM
      first.valid A B ws.a2 ws.b2 reconstructed.lenA reconstructed.lenB
  earlyExec :
    let ws := hgcdRecursiveWorkspace W lenA
    hgcdRecursiveEarlyReturn M first.matrix hM first.valid computeM A B ws.a2
      ws.b2 reconstructed.lenA reconstructed.lenB first.sgn
      reconstructed.heap = .ok early

/-- Physical evidence for the successful non-early continuation of an actual
non-base invocation.  Both child executions are retained, while their
polynomial semantics remain external well-founded induction hypotheses. -/
structure HgcdRecursiveNonEarlyContinuationWorkspace (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (bound : Nat) (recurse : HgcdRecursiveCallBelow bound)
    (M : HgcdMat) (hM : M.Valid) (computeM : Bool)
    (A B a b : RawPtr UInt64)
    (W scratch : RawPtr UInt64) (lenA lenB : Nat) (heap : RawHeap)
    (package : HgcdRecursiveNonBasePackage this bound recurse a b W scratch
      lenA lenB heap)
    (hbound : lenA = bound) (horder : lenB < lenA)
    (hbase : ¬ lenB < lenA / 2 + 1) : Type where
  first : HgcdRecursiveResult
  reconstructed : HgcdRecursiveReconstructPairResult
  middle : HgcdRecursiveMiddleResult
  second : HgcdRecursiveResult
  finished : HgcdRecursiveFinishResult
  firstExec :
    let ws := hgcdRecursiveWorkspace W lenA
    let high := hgcdRecursiveHighInput a b lenA lenB
    hgcdRecursiveDispatchBelow this bound recurse ws.R
      (hgcdRecursiveWorkspace_R_valid W lenA) ws.a3 ws.b3 high.a0 high.b0
      high.lenA0 high.lenB0 ws.q ws.W3 ws.T0 ws.T1 scratch ws.a2 ws.next
      heap (hgcdRecursiveHighInput_order a b lenA lenB horder)
        (by
          rw [← hbound]
          exact hgcdRecursiveHighInput_len_lt a b lenA lenB horder
            (by omega)) = .ok first
  reconstructWorkspace :
    let ws := hgcdRecursiveWorkspace W lenA
    HgcdRecursiveReconstructPairWorkspaceProvider this ws.a2 ws.b2 ws.T0
      a b ws.a3 ws.b3 scratch (Nat.min lenA (lenA / 2))
      (Nat.min lenB (lenA / 2)) first.lenA first.lenB (lenA / 2)
      first.matrix first.valid first.sgn first.heap
  reconstructExec :
    let ws := hgcdRecursiveWorkspace W lenA
    hgcdRecursiveReconstructPair this ws.a2 ws.b2 ws.T0 a b ws.a3 ws.b3
      scratch (Nat.min lenA (lenA / 2)) (Nat.min lenB (lenA / 2))
      first.lenA first.lenB (lenA / 2) first.matrix first.valid first.sgn
      first.heap = .ok reconstructed
  notEarly : ¬ reconstructed.lenB < lenA / 2 + 1
  middleWorkspace : HgcdRecursiveMiddleWorkspace W lenA reconstructed
  middleExec :
    let ws := hgcdRecursiveWorkspace W lenA
    hgcdRecursiveMiddle this ws.q ws.d ws.a2 ws.b2 reconstructed.lenA
      reconstructed.lenB (lenA / 2) ws.W3 reconstructed.heap = .ok middle
  secondOrder : middle.lenD0 < middle.lenC0
  secondDecrease : middle.lenC0 < bound
  secondWorkspace : HgcdRecursiveSecondCallWorkspace this bound recurse M hM
    A B W scratch lenA first reconstructed middle second secondOrder
    secondDecrease
  secondExec :
    let ws := hgcdRecursiveWorkspace W lenA
    hgcdRecursiveDispatchBelow this bound recurse ws.S
      (hgcdRecursiveWorkspace_S_valid W lenA) ws.a3 ws.b3 middle.c0
      middle.d0 middle.lenC0 middle.lenD0 ws.a2 ws.W3 ws.T0 ws.T1 scratch
      ws.a2 ws.next middle.heap secondOrder secondDecrease = .ok second
  finishExec :
    let ws := hgcdRecursiveWorkspace W lenA
    hgcdRecursiveFinish this M first.matrix second.matrix hM first.valid
      second.valid computeM A B ws.T0 ws.b2 ws.d ws.a3 ws.b3 ws.q
      (Nat.min reconstructed.lenB middle.k) (Nat.min middle.lenD middle.k)
      second.lenA second.lenB middle.k middle.lenQ ws.a2 scratch first.sgn
      second.sgn second.heap = .ok finished

/-- Attach the first-child well-founded induction theorem to a physical early
continuation and obtain the common raw invariant for its actual result. -/
theorem HgcdRecursiveEarlyContinuationWorkspace.rawInvariant
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (bound : Nat) (recurse : HgcdRecursiveCallBelow bound)
    (M : HgcdMat) (hM : M.Valid) (computeM : Bool)
    (A B a b : RawPtr UInt64) (lenA lenB : Nat)
    (W scratch : RawPtr UInt64) (heap : RawHeap)
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (hbound : lenA = bound) (horder : lenB < lenA)
    (hbase : ¬ lenB < lenA / 2 + 1)
    (package : HgcdRecursiveNonBasePackage this bound recurse a b W scratch
      lenA lenB heap)
    (recursiveRefines : HgcdRecursiveNonBaseCallbackRefines this bound recurse
      a b W scratch lenA lenB heap package)
    (continuation : HgcdRecursiveEarlyContinuationWorkspace this bound recurse
      M hM computeM A B a b W scratch lenA lenB heap package hbound horder
      hbase)
    (result : HgcdRecursiveResult)
    (hrun :
      let firstCall := package.admissible this bound recurse a b W scratch
        lenA lenB heap recursiveRefines
      let providers := hgcdRecursiveFirstCall_providers this bound recurse a b
        W scratch lenA lenB heap package.inputHighA package.inputHighB
        package.lowPolyA package.lowPolyB hcfg hp horder firstCall.workspace
        firstCall.recursiveRefines
      hgcdRecursiveBodyBelow this bound recurse M hM computeM A B a b lenA
        lenB W scratch heap hbound horder (fun _ => providers.1)
          (fun _ => providers.2) = .ok result) :
    ∃ finalA finalB entries,
      HgcdRecursiveRawInvariant this package.inputA package.inputB finalA
        finalB entries computeM A B lenA result := by
  let firstCall := package.admissible this bound recurse a b W scratch lenA
    lenB heap recursiveRefines
  exact hgcdRecursiveBodyBelow_early_rawInvariant this bound recurse M hM
    computeM A B a b lenA lenB W scratch heap package.inputA package.inputB
    package.inputHighA package.inputHighB package.lowPolyA package.lowPolyB hp
    hcfg hbound horder firstCall continuation.first continuation.reconstructed
    continuation.early hbase continuation.firstExec package.splitA
    package.splitB continuation.reconstructWorkspace
    continuation.reconstructExec continuation.guard continuation.earlyWorkspace
    continuation.earlyExec result hrun

/-- Attach exactly the two smaller-call induction theorems to a physical
non-early continuation and obtain the common invariant for the generated
finish result. -/
theorem HgcdRecursiveNonEarlyContinuationWorkspace.rawInvariant
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (bound : Nat) (recurse : HgcdRecursiveCallBelow bound)
    (M : HgcdMat) (hM : M.Valid) (computeM : Bool)
    (A B a b : RawPtr UInt64) (lenA lenB : Nat)
    (W scratch : RawPtr UInt64) (heap : RawHeap)
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (hbound : lenA = bound) (horder : lenB < lenA)
    (hbase : ¬ lenB < lenA / 2 + 1)
    (package : HgcdRecursiveNonBasePackage this bound recurse a b W scratch
      lenA lenB heap)
    (firstRecursiveRefines : HgcdRecursiveNonBaseCallbackRefines this bound
      recurse a b W scratch lenA lenB heap package)
    (continuation : HgcdRecursiveNonEarlyContinuationWorkspace this bound
      recurse M hM computeM A B a b W scratch lenA lenB heap package hbound
      horder hbase)
    (secondRecursiveRefines : HgcdRecursiveSecondCallbackRefines this bound
      recurse W scratch lenA continuation.middle continuation.secondOrder
      continuation.secondDecrease)
    (result : HgcdRecursiveResult)
    (hrun :
      let firstCall := package.admissible this bound recurse a b W scratch
        lenA lenB heap firstRecursiveRefines
      let providers := hgcdRecursiveFirstCall_providers this bound recurse a b
        W scratch lenA lenB heap package.inputHighA package.inputHighB
        package.lowPolyA package.lowPolyB hcfg hp horder firstCall.workspace
        firstCall.recursiveRefines
      hgcdRecursiveBodyBelow this bound recurse M hM computeM A B a b lenA
        lenB W scratch heap hbound horder (fun _ => providers.1)
          (fun _ => providers.2) = .ok result) :
    ∃ finalA finalB entries,
      HgcdRecursiveRawInvariant this package.inputA package.inputB finalA
        finalB entries computeM A B lenA result := by
  let firstCall := package.admissible this bound recurse a b W scratch lenA
    lenB heap firstRecursiveRefines
  exact hgcdRecursiveBodyBelow_nonEarly_rawInvariant this bound recurse M hM
    computeM A B a b lenA lenB W scratch heap package.inputA package.inputB
    package.inputHighA package.inputHighB package.lowPolyA package.lowPolyB hp
    hcfg hbound horder firstCall continuation.first continuation.reconstructed
    continuation.middle continuation.second continuation.finished hbase
    continuation.firstExec package.splitA package.splitB
    continuation.reconstructWorkspace continuation.reconstructExec
    continuation.notEarly continuation.middleWorkspace continuation.middleExec
    continuation.secondOrder continuation.secondDecrease
    continuation.secondWorkspace secondRecursiveRefines continuation.secondExec
    continuation.finishExec result hrun

/-- The two successful continuations of a generated non-base body.  This is a
physical sum: neither constructor contains a recursive semantic conclusion. -/
inductive HgcdRecursiveNonBaseContinuationWorkspace (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (bound : Nat) (recurse : HgcdRecursiveCallBelow bound)
    (M : HgcdMat) (hM : M.Valid) (computeM : Bool)
    (A B a b : RawPtr UInt64) (W scratch : RawPtr UInt64)
    (lenA lenB : Nat) (heap : RawHeap)
    (package : HgcdRecursiveNonBasePackage this bound recurse a b W scratch
      lenA lenB heap)
    (hbound : lenA = bound) (horder : lenB < lenA)
    (hbase : ¬ lenB < lenA / 2 + 1) : Type where
  | early
      (workspace : HgcdRecursiveEarlyContinuationWorkspace this bound recurse
        M hM computeM A B a b W scratch lenA lenB heap package hbound horder
        hbase)
  | nonEarly
      (workspace : HgcdRecursiveNonEarlyContinuationWorkspace this bound
        recurse M hM computeM A B a b W scratch lenA lenB heap package hbound
        horder hbase)

/-- The semantic obligation selected by an actual continuation: the early
constructor has no second child, while the non-early constructor receives the
strictly-smaller second-child induction theorem. -/
def HgcdRecursiveContinuationSecondRefines (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (bound : Nat) (recurse : HgcdRecursiveCallBelow bound)
    (M : HgcdMat) (hM : M.Valid) (computeM : Bool)
    (A B a b : RawPtr UInt64) (W scratch : RawPtr UInt64)
    (lenA lenB : Nat) (heap : RawHeap)
    (package : HgcdRecursiveNonBasePackage this bound recurse a b W scratch
      lenA lenB heap)
    (hbound : lenA = bound) (horder : lenB < lenA)
    (hbase : ¬ lenB < lenA / 2 + 1)
    (continuation : HgcdRecursiveNonBaseContinuationWorkspace this bound
      recurse M hM computeM A B a b W scratch lenA lenB heap package hbound
      horder hbase) : Prop :=
  match continuation with
  | .early _ => True
  | .nonEarly workspace =>
      HgcdRecursiveSecondCallbackRefines this bound recurse W scratch lenA
        workspace.middle workspace.secondOrder workspace.secondDecrease

/-- A successful physical non-base continuation plus exactly its selected
well-founded child hypotheses proves the common raw invariant. -/
theorem HgcdRecursiveNonBaseContinuationWorkspace.rawInvariant
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (bound : Nat) (recurse : HgcdRecursiveCallBelow bound)
    (M : HgcdMat) (hM : M.Valid) (computeM : Bool)
    (A B a b : RawPtr UInt64) (lenA lenB : Nat)
    (W scratch : RawPtr UInt64) (heap : RawHeap)
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (hbound : lenA = bound) (horder : lenB < lenA)
    (hbase : ¬ lenB < lenA / 2 + 1)
    (package : HgcdRecursiveNonBasePackage this bound recurse a b W scratch
      lenA lenB heap)
    (firstRecursiveRefines : HgcdRecursiveNonBaseCallbackRefines this bound
      recurse a b W scratch lenA lenB heap package)
    (continuation : HgcdRecursiveNonBaseContinuationWorkspace this bound
      recurse M hM computeM A B a b W scratch lenA lenB heap package hbound
      horder hbase)
    (secondRecursiveRefines : HgcdRecursiveContinuationSecondRefines this bound
      recurse M hM computeM A B a b W scratch lenA lenB heap package hbound
      horder hbase continuation)
    (result : HgcdRecursiveResult)
    (hrun :
      let firstCall := package.admissible this bound recurse a b W scratch
        lenA lenB heap firstRecursiveRefines
      let providers := hgcdRecursiveFirstCall_providers this bound recurse a b
        W scratch lenA lenB heap package.inputHighA package.inputHighB
        package.lowPolyA package.lowPolyB hcfg hp horder firstCall.workspace
        firstCall.recursiveRefines
      hgcdRecursiveBodyBelow this bound recurse M hM computeM A B a b lenA
        lenB W scratch heap hbound horder (fun _ => providers.1)
          (fun _ => providers.2) = .ok result) :
    ∃ finalA finalB entries,
      HgcdRecursiveRawInvariant this package.inputA package.inputB finalA
        finalB entries computeM A B lenA result := by
  cases continuation with
  | early workspace =>
      exact workspace.rawInvariant this bound recurse M hM computeM A B a b
        lenA lenB W scratch heap hcfg hp hbound horder hbase package
        firstRecursiveRefines result hrun
  | nonEarly workspace =>
      exact workspace.rawInvariant this bound recurse M hM computeM A B a b
        lenA lenB W scratch heap hcfg hp hbound horder hbase package
        firstRecursiveRefines secondRecursiveRefines result hrun

/-- Invoke the well-founded body with recursive admissibility demanded only
on the source's actual non-base branch.  Thus the base execution carries no
child workspace or unexecuted recursive premise. -/
def hgcdRecursiveBodyBranchAdmissible (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (bound : Nat) (recurse : HgcdRecursiveCallBelow bound)
    (M : HgcdMat) (hM : M.Valid) (computeM : Bool)
    (A B a b : RawPtr UInt64) (lenA lenB : Nat)
    (W scratch : RawPtr UInt64) (heap : RawHeap)
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (hbound : lenA = bound) (horder : lenB < lenA)
    (nonBase : ¬ lenB < lenA / 2 + 1 →
      HgcdRecursiveNonBasePackage this bound recurse a b W scratch lenA lenB
        heap)
    (recursiveRefines : ∀ hbase,
      HgcdRecursiveNonBaseCallbackRefines this bound recurse a b W scratch
        lenA lenB heap (nonBase hbase)) : RawExec HgcdRecursiveResult :=
  if hbase : lenB < lenA / 2 + 1 then
    hgcdRecursiveBodyBelow this bound recurse M hM computeM A B a b lenA lenB
      W scratch heap hbound horder
        (fun hnonbase => False.elim (hnonbase hbase))
        (fun hnonbase => False.elim (hnonbase hbase))
  else
    let package := nonBase hbase
    let providers := hgcdRecursiveFirstCall_providers this bound recurse a b W
      scratch lenA lenB heap package.inputHighA package.inputHighB
      package.lowPolyA package.lowPolyB hcfg hp horder
      package.workspace (recursiveRefines hbase)
    hgcdRecursiveBodyBelow this bound recurse M hM computeM A B a b lenA lenB
      W scratch heap hbound horder (fun _ => providers.1)
        (fun _ => providers.2)

/-- Uniform semantic theorem for one branch-admissible well-founded body.
The caller supplies physical evidence only for the branch that the generated
guard actually takes.  Recursive polynomial semantics are separate smaller-
length hypotheses, and successful continuation evidence is tied to the exact
body execution rather than to a preselected result. -/
theorem hgcdRecursiveBodyBranchAdmissible_rawInvariant
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (bound : Nat) (recurse : HgcdRecursiveCallBelow bound)
    (M : HgcdMat) (hM : M.Valid) (computeM : Bool)
    (A B a b : RawPtr UInt64) (lenA lenB : Nat)
    (W scratch : RawPtr UInt64) (heap : RawHeap)
    (left right : Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (hbound : lenA = bound) (horder : lenB < lenA)
    (hlenAPos : 0 < lenA)
    (baseWorkspace : lenB < lenA / 2 + 1 →
      HgcdRecursiveBaseCallWorkspace this M hM A B a b lenA lenB heap left
        right)
    (nonBase : ∀ hbase : ¬ lenB < lenA / 2 + 1,
      HgcdRecursiveNonBasePackage this bound recurse a b W scratch lenA lenB
        heap)
    (nonBaseInput : ∀ hbase,
      (nonBase hbase).inputA = left ∧ (nonBase hbase).inputB = right)
    (firstRecursiveRefines : ∀ hbase,
      HgcdRecursiveNonBaseCallbackRefines this bound recurse a b W scratch
        lenA lenB heap (nonBase hbase))
    (continuation : ∀ (hbase : ¬ lenB < lenA / 2 + 1)
      (result : HgcdRecursiveResult),
      let package := nonBase hbase
      let firstCall := package.admissible this bound recurse a b W scratch
        lenA lenB heap (firstRecursiveRefines hbase)
      let providers := hgcdRecursiveFirstCall_providers this bound recurse a b
        W scratch lenA lenB heap package.inputHighA package.inputHighB
        package.lowPolyA package.lowPolyB hcfg hp horder firstCall.workspace
        firstCall.recursiveRefines
      hgcdRecursiveBodyBelow this bound recurse M hM computeM A B a b lenA
        lenB W scratch heap hbound horder (fun _ => providers.1)
          (fun _ => providers.2) = .ok result →
      HgcdRecursiveNonBaseContinuationWorkspace this bound recurse M hM
        computeM A B a b W scratch lenA lenB heap package hbound horder hbase)
    (secondRecursiveRefines : ∀ (hbase : ¬ lenB < lenA / 2 + 1)
      (result : HgcdRecursiveResult)
      (hrun :
        let package := nonBase hbase
        let firstCall := package.admissible this bound recurse a b W scratch
          lenA lenB heap (firstRecursiveRefines hbase)
        let providers := hgcdRecursiveFirstCall_providers this bound recurse a
          b W scratch lenA lenB heap package.inputHighA package.inputHighB
          package.lowPolyA package.lowPolyB hcfg hp horder firstCall.workspace
          firstCall.recursiveRefines
        hgcdRecursiveBodyBelow this bound recurse M hM computeM A B a b lenA
          lenB W scratch heap hbound horder (fun _ => providers.1)
            (fun _ => providers.2) = .ok result),
      HgcdRecursiveContinuationSecondRefines this bound recurse M hM computeM
        A B a b W scratch lenA lenB heap (nonBase hbase) hbound horder hbase
        (continuation hbase result hrun))
    (result : HgcdRecursiveResult)
    (hrun : hgcdRecursiveBodyBranchAdmissible this bound recurse M hM computeM
      A B a b lenA lenB W scratch heap hcfg hp hbound horder nonBase
        firstRecursiveRefines = .ok result) :
    ∃ finalA finalB entries,
      HgcdRecursiveRawInvariant this left right finalA finalB entries computeM
        A B lenA result := by
  by_cases hbase : lenB < lenA / 2 + 1
  · have hrunBelow :
        hgcdRecursiveBodyBelow this bound recurse M hM computeM A B a b lenA
          lenB W scratch heap hbound horder
            (fun hnonbase => False.elim (hnonbase hbase))
            (fun hnonbase => False.elim (hnonbase hbase)) = .ok result := by
      simpa [hgcdRecursiveBodyBranchAdmissible, hbase] using hrun
    rcases hgcdRecursiveBodyBelow_base_rawInvariant this bound recurse M hM
        computeM A B a b lenA lenB W scratch heap left right hp hbound horder
        (baseWorkspace hbase) hlenAPos hbase result hrunBelow with
      ⟨entries, hInvariant⟩
    exact ⟨left, right, entries, hInvariant⟩
  · let package := nonBase hbase
    let firstRefines := firstRecursiveRefines hbase
    let firstCall := package.admissible this bound recurse a b W scratch lenA
      lenB heap firstRefines
    let providers := hgcdRecursiveFirstCall_providers this bound recurse a b W
      scratch lenA lenB heap package.inputHighA package.inputHighB
      package.lowPolyA package.lowPolyB hcfg hp horder firstCall.workspace
      firstCall.recursiveRefines
    have hrunBelow :
        hgcdRecursiveBodyBelow this bound recurse M hM computeM A B a b lenA
          lenB W scratch heap hbound horder (fun _ => providers.1)
            (fun _ => providers.2) = .ok result := by
      simpa [hgcdRecursiveBodyBranchAdmissible, hbase, package, firstRefines,
        firstCall, providers] using hrun
    let next := continuation hbase result (by
      simpa [package, firstRefines, firstCall, providers] using hrunBelow)
    have hSecond := secondRecursiveRefines hbase result (by
      simpa [package, firstRefines, firstCall, providers] using hrunBelow)
    have hResult := next.rawInvariant this bound recurse M hM computeM A B a b
      lenA lenB W scratch heap hcfg hp hbound horder hbase package firstRefines
      (by simpa [next, package] using hSecond) result hrunBelow
    rcases nonBaseInput hbase with ⟨hleft, hright⟩
    simpa [package, hleft, hright] using hResult

/-- Physical safety node for one invocation of a source-shaped recursive
callback.  It is parameterized over every possible first-child induction
proof, so the node cannot manufacture or store recursive polynomial
semantics.  Its continuation is additionally tied to the exact successful
execution obtained with that proof. -/
structure HgcdRecursiveInvocationWorkspace (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (plain : HgcdRecursiveCall) (bound : Nat)
    (M : HgcdMat) (hM : M.Valid) (computeM : Bool)
    (A B a b : RawPtr UInt64) (lenA lenB : Nat)
    (W scratch : RawPtr UInt64) (heap : RawHeap)
    (left right : Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (hbound : lenA = bound) (horder : lenB < lenA) : Type where
  positive : 0 < lenA
  base : ∀ hbase : lenB < lenA / 2 + 1,
    HgcdRecursiveBaseCallWorkspace this M hM A B a b lenA lenB heap left right
  nonBase : ∀ hbase : ¬ lenB < lenA / 2 + 1,
    HgcdRecursiveNonBasePackage this bound
      (hgcdRecursiveCallBelowOfCall bound plain) a b W scratch lenA lenB heap
  nonBaseInput : ∀ hbase,
    (nonBase hbase).inputA = left ∧ (nonBase hbase).inputB = right
  continuation : ∀ (hbase : ¬ lenB < lenA / 2 + 1)
    (firstRecursiveRefines :
      HgcdRecursiveNonBaseCallbackRefines this bound
        (hgcdRecursiveCallBelowOfCall bound plain) a b W scratch lenA lenB
        heap (nonBase hbase))
    (result : HgcdRecursiveResult),
    let package := nonBase hbase
    let firstCall := package.admissible this bound
      (hgcdRecursiveCallBelowOfCall bound plain) a b W scratch lenA lenB heap
      firstRecursiveRefines
    let providers := hgcdRecursiveFirstCall_providers this bound
      (hgcdRecursiveCallBelowOfCall bound plain) a b W scratch lenA lenB heap
      package.inputHighA package.inputHighB package.lowPolyA package.lowPolyB
      hcfg hp horder firstCall.workspace firstCall.recursiveRefines
    hgcdRecursiveBodyBelow this bound
      (hgcdRecursiveCallBelowOfCall bound plain) M hM computeM A B a b lenA
      lenB W scratch heap hbound horder (fun _ => providers.1)
        (fun _ => providers.2) = .ok result →
    HgcdRecursiveNonBaseContinuationWorkspace this bound
      (hgcdRecursiveCallBelowOfCall bound plain) M hM computeM A B a b W
      scratch lenA lenB heap package hbound horder hbase

/-- Hereditary physical safety for every strictly smaller represented call
reachable below one invocation length.  The provider is conditional on the
actual raw polynomial representations, so it grants nothing for arbitrary or
invalid pointers. -/
def HgcdRecursiveInvocationWorkspaceProviderBelow (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (plain : HgcdRecursiveCall) (hcfg : DensePreinvConfigured this)
    (hp : 1 < this._p.toNat) (limit : Nat) : Type :=
  ∀ (childLen : Nat), childLen < limit →
    ∀ (M : HgcdMat) (hM : M.Valid) (computeM : Bool)
      (A B a b : RawPtr UInt64) (lenB : Nat)
      (W scratch : RawPtr UInt64) (heap : RawHeap)
      (left right : Polynomial (ZMod this._p.toNat))
      (horder : lenB < childLen),
    RawDensePolyRep this heap a childLen left →
    RawDensePolyRep this heap b lenB right →
    HgcdRecursiveInvocationWorkspace this plain childLen M hM computeM A B a b
      childLen lenB W scratch heap left right hcfg hp rfl horder

/-- Restrict a hereditary workspace provider to a strictly smaller limit. -/
def HgcdRecursiveInvocationWorkspaceProviderBelow.mono
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (plain : HgcdRecursiveCall) (hcfg : DensePreinvConfigured this)
    (hp : 1 < this._p.toNat) (large small : Nat) (hlt : small < large)
    (provider : HgcdRecursiveInvocationWorkspaceProviderBelow this plain hcfg
      hp large) :
    HgcdRecursiveInvocationWorkspaceProviderBelow this plain hcfg hp small :=
  fun childLen hchild => provider childLen (hchild.trans hlt)

end CLPoly.Impl.StrictHGCDRawRefinement
