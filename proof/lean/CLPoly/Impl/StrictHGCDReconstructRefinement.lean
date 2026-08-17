import CLPoly.Impl.StrictHGCDFirstCallRefinement

set_option autoImplicit false

namespace CLPoly.Impl.StrictHGCDRawRefinement

open Generated.StrictHGCD
open CLPoly.Impl.RawPolynomialRep
open CLPoly.Impl.StrictDivremRefinement
open CLPoly.Impl.StrictEuclidRefinement
open CLPoly.Impl.StrictPolyAddSubRefinement
open CLPoly.Impl.StrictMulRefinement
open CLPoly.Impl.StrictWordArithmetic

/-- The first reconstructed pair already carries the complete parent length
contract when the generated early-stop guard succeeds.  Matrix descriptors
are those returned by the real first child; only the operand descriptors are
the concrete reconstruction outputs. -/
theorem hgcdRecursiveEarly_lengthInvariant
    (outerLength highLength m : Nat) (first : HgcdRecursiveResult)
    (reconstructed : HgcdRecursiveReconstructPairResult)
    (heap : RawHeap) (sgn : Int)
    (hm : m = outerLength / 2)
    (hhigh : highLength = outerLength - m)
    (hfirst : HgcdRecursiveLengthInvariant highLength first)
    (hreconstructed : HgcdFirstReconstructionInvariant outerLength first
      reconstructed)
    (hearly : reconstructed.lenB < m + 1) :
    HgcdRecursiveLengthInvariant outerLength
      ⟨heap, first.matrix, first.valid, reconstructed.lenA,
        reconstructed.lenB, sgn⟩ := by
  have houter : highLength + m = outerLength := by omega
  have hrow (i : Fin 4)
      (h : hgcdMatLenRaw first.matrix first.valid i + first.lenA ≤
        highLength + 1) :
      hgcdMatLenRaw first.matrix first.valid i + reconstructed.lenA ≤
        outerLength + 1 := by
    rw [hreconstructed.leadingA]
    omega
  have hcoeff (i : Fin 4) :
      hgcdMatLenRaw first.matrix first.valid i ≤
        outerLength - outerLength / 2 := by
    have := hfirst.coeffBound i
    omega
  exact {
    row0A := hrow 0 hfirst.row0A
    row1B := by
      have hrow1 := hrow 1 hfirst.row1A
      have hx : hgcdMatLenRaw first.matrix first.valid (1 : Fin 4) +
          reconstructed.lenB ≤ outerLength + 1 := by
        have := hreconstructed.ordered
        omega
      simpa using hx
    row2A := hrow 2 hfirst.row2A
    row3B := by
      have hrow3 := hrow 3 hfirst.row3A
      have hx : hgcdMatLenRaw first.matrix first.valid (3 : Fin 4) +
          reconstructed.lenB ≤ outerLength + 1 := by
        have := hreconstructed.ordered
        omega
      simpa using hx
    row1A := hrow 1 hfirst.row1A
    row3A := hrow 3 hfirst.row3A
    inputBound := by
      have hx : reconstructed.lenA ≤ outerLength := by
        rw [hreconstructed.leadingA]
        have hf := hfirst.inputBound
        omega
      simpa using hx
    stopped := by simpa [hm] using hearly
    positive := hreconstructed.positiveA
    aboveHalf := by
      have hx : outerLength / 2 < reconstructed.lenA := by
        rw [hreconstructed.leadingA, ← hm]
        exact Nat.lt_add_of_pos_right hfirst.positive
      simpa using hx
    coeffBound := hcoeff }

/-- Substituting the exact source split `k = 2*m-lenB2+1` into the exact
final-A reconstruction length shows that the complete recursive result stays
strictly above the outer half-length threshold. -/
theorem hgcdRecursiveFinalReconstruct_lenA_above_half
    (outerLength m reconstructedLenB k lenC0 secondLenA resultLenA : Nat)
    (hm : m = outerLength / 2)
    (hk : k = 2 * m - reconstructedLenB + 1)
    (hc : lenC0 = reconstructedLenB - k)
    (hreconstructedLower : m + 1 ≤ reconstructedLenB)
    (hreconstructedUpper : reconstructedLenB < outerLength)
    (hsecondAbove : lenC0 / 2 < secondLenA)
    (hresult : resultLenA = k + secondLenA) :
    outerLength / 2 < resultLenA := by
  omega

/-- The B output of the second reconstruction satisfies the outer HGCD stop
threshold.  This uses the source's exact `k = 2*m-lenB2+1`, the second call's
leading-A/stop pair, and the two matrix rows that physically build B. -/
theorem hgcdRecursiveFinalReconstruct_lenB_lt_half
    (outerLength m reconstructedLenB k lenC0 lenD : Nat)
    (secondLenA secondLenB s0 s2 resultLenB : Nat)
    (hm : m = outerLength / 2)
    (hk : k = 2 * m - reconstructedLenB + 1)
    (hc : lenC0 = reconstructedLenB - k)
    (hreconstructedLower : m + 1 ≤ reconstructedLenB)
    (hreconstructedUpper : reconstructedLenB < outerLength)
    (hsecondAbove : lenC0 / 2 < secondLenA)
    (hsecondStop : secondLenB < lenC0 / 2 + 1)
    (hrow0 : s0 + secondLenA ≤ lenC0 + 1)
    (hrow2 : s2 + secondLenA ≤ lenC0 + 1)
    (hresultB : resultLenB ≤ max (k + secondLenB)
      (max
        (s2 + Nat.min reconstructedLenB k - 1)
        (s0 + Nat.min lenD k - 1))) :
    resultLenB < outerLength / 2 + 1 := by
  have hminReconstructed : Nat.min reconstructedLenB k ≤ k :=
    Nat.min_le_right _ _
  have hminD : Nat.min lenD k ≤ k := Nat.min_le_right _ _
  have hkLe : k ≤ reconstructedLenB := by omega
  have hsplit : k + lenC0 = reconstructedLenB := by omega
  have hhalf : k + lenC0 / 2 ≤ m := by omega
  rw [← hm]
  apply lt_of_le_of_lt hresultB
  apply max_lt
  · omega
  · apply max_lt <;> omega

/-- Operand-length portion of the recursive contract returned by the real
non-early finish branch.  Matrix descriptor fields are supplied separately
by the concrete combine-matrix execution. -/
structure HgcdRecursiveFinishOperandInvariant (outerLength : Nat)
    (result : HgcdRecursiveFinishResult) : Prop where
  inputBoundA : result.lenA ≤ outerLength
  inputBoundB : result.lenB ≤ outerLength
  positiveA : 0 < result.lenA
  aboveHalf : outerLength / 2 < result.lenA
  stopped : result.lenB < outerLength / 2 + 1

set_option maxHeartbeats 800000 in
/-- The exact second reconstruction executed inside `hgcdRecursiveFinish`
establishes every operand-length field required by its recursive parent.
All bounds are derived from the returned second-child invariant and the
source formulas for `k`, `c0`, and the two low slices. -/
theorem hgcdRecursiveFinish_operandInvariant (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (M R : HgcdMat) (hM : M.Valid) (hR : R.Valid)
    (computeM : Bool) (A B T0 lowA lowB highA highB q : RawPtr UInt64)
    (outerLength m reconstructedLenB lenD k secondInputLength lenQ : Nat)
    (a2 scratch : RawPtr UInt64) (sgnR : Int)
    (second : HgcdRecursiveResult) (result : HgcdRecursiveFinishResult)
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (polyLowA polyLowB polyHighA polyHighB :
      Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (hinvariant : HgcdRecursiveLengthInvariant secondInputLength second)
    (physical : HgcdRecursiveReconstructPairWorkspaceProvider this A B T0
      lowA lowB highA highB scratch (Nat.min reconstructedLenB k)
      (Nat.min lenD k) second.lenA second.lenB k second.matrix second.valid
      second.sgn second.heap)
    (hMatrix : HgcdMatRawDenseRep this second.heap second.matrix entries
      second.valid)
    (hLowA : RawCanonicalPolySlice this second.heap lowA
      (Nat.min reconstructedLenB k) polyLowA)
    (hLowB : RawCanonicalPolySlice this second.heap lowB
      (Nat.min lenD k) polyLowB)
    (hHighA : RawDensePolyRep this second.heap highA second.lenA polyHighA)
    (hHighB : RawDensePolyRep this second.heap highB second.lenB polyHighB)
    (hm : m = outerLength / 2)
    (hk : k = 2 * m - reconstructedLenB + 1)
    (hc : secondInputLength = reconstructedLenB - k)
    (hreconstructedLower : m + 1 ≤ reconstructedLenB)
    (hreconstructedUpper : reconstructedLenB < outerLength)
    (hrun : hgcdRecursiveFinish this M R second.matrix hM hR second.valid
      computeM A B T0 lowA lowB highA highB q
      (Nat.min reconstructedLenB k) (Nat.min lenD k) second.lenA second.lenB
      k lenQ a2 scratch sgnR second.sgn second.heap = .ok result) :
    HgcdRecursiveFinishOperandInvariant outerLength result := by
  rcases hgcdRecursiveFinish_exec this M R second.matrix hM hR second.valid
      computeM A B T0 lowA lowB highA highB q
      (Nat.min reconstructedLenB k) (Nat.min lenD k) second.lenA second.lenB
      k lenQ a2 scratch sgnR second.sgn second.heap result hrun with
    ⟨reconstructed, hreconstruct, hlenA, hlenB, _, _⟩
  have hrefines := hgcdRecursiveReconstructPair_refines this A B T0 lowA
    lowB highA highB scratch (Nat.min reconstructedLenB k)
    (Nat.min lenD k) second.lenA second.lenB k second.matrix second.valid
    second.sgn second.heap reconstructed entries polyLowA polyLowB polyHighA
    polyHighB hcfg hp physical hMatrix hLowA hLowB hHighA hHighB
    hreconstruct
  have hleading := hgcdRecursiveFinalReconstruct_lenA_eq_of_invariant this A
    B T0 lowA lowB highA highB scratch (Nat.min reconstructedLenB k)
    (Nat.min lenD k) k secondInputLength second reconstructed entries polyLowA
    polyLowB polyHighA polyHighB hcfg hp hinvariant (Nat.min_le_right _ _)
    (Nat.min_le_right _ _) physical hMatrix hLowA hLowB hHighA hHighB
    hreconstruct
  have hlengths := hgcdRecursiveFinalReconstruct_lengths_le_input outerLength
    reconstructedLenB k secondInputLength lenD second.lenA second.lenB
    (hgcdMatLen second.matrix second.valid (0 : Fin 4))
    (hgcdMatLen second.matrix second.valid (1 : Fin 4))
    (hgcdMatLen second.matrix second.valid (2 : Fin 4))
    (hgcdMatLen second.matrix second.valid (3 : Fin 4)) reconstructed.lenA
    reconstructed.lenB (by omega) (Nat.le_of_lt hreconstructedUpper)
    hinvariant.inputBound hinvariant.stopped
    (by simpa [hgcdMatLen, hgcdMatLenRaw] using hinvariant.row0A)
    (by simpa [hgcdMatLen, hgcdMatLenRaw] using hinvariant.row1B)
    (by simpa [hgcdMatLen, hgcdMatLenRaw] using hinvariant.row2A)
    (by simpa [hgcdMatLen, hgcdMatLenRaw] using hinvariant.row3B)
    hrefines.2.2.1 hrefines.2.2.2.1
  have habove := hgcdRecursiveFinalReconstruct_lenA_above_half outerLength m
    reconstructedLenB k secondInputLength second.lenA reconstructed.lenA hm
    hk hc hreconstructedLower hreconstructedUpper hinvariant.aboveHalf
    hleading.1
  have hstop := hgcdRecursiveFinalReconstruct_lenB_lt_half outerLength m
    reconstructedLenB k secondInputLength lenD second.lenA second.lenB
    (hgcdMatLen second.matrix second.valid (0 : Fin 4))
    (hgcdMatLen second.matrix second.valid (2 : Fin 4)) reconstructed.lenB hm
    hk hc hreconstructedLower hreconstructedUpper hinvariant.aboveHalf
    hinvariant.stopped
    (by simpa [hgcdMatLen, hgcdMatLenRaw] using hinvariant.row0A)
    (by simpa [hgcdMatLen, hgcdMatLenRaw] using hinvariant.row2A)
    hrefines.2.2.2.1
  exact {
    inputBoundA := by simpa [hlenA] using hlengths.1
    inputBoundB := by simpa [hlenB] using hlengths.2
    positiveA := by simpa [hlenA] using hleading.2
    aboveHalf := by simpa [hlenA] using habove
    stopped := by simpa [hlenB] using hstop }

/-- Assemble the common recursive length contract from the finish operand
facts and the two descriptor facts returned by the real matrix block. -/
theorem hgcdRecursiveLengthInvariant_of_finish
    (outerLength : Nat) (result : HgcdRecursiveFinishResult)
    (hOperands : HgcdRecursiveFinishOperandInvariant outerLength result)
    (hRows : ∀ i : Fin 4,
      hgcdMatLen result.matrix result.valid i + result.lenA ≤ outerLength + 1)
    (hCoeff : ∀ i : Fin 4,
      hgcdMatLen result.matrix result.valid i ≤
        outerLength - outerLength / 2) :
    HgcdRecursiveLengthInvariant outerLength result.toResult := by
  have horder : result.lenB ≤ result.lenA := by
    have hstop := hOperands.stopped
    have habove := hOperands.aboveHalf
    omega
  exact {
    row0A := by simpa [HgcdRecursiveFinishResult.toResult, hgcdMatLenRaw,
      hgcdMatLen] using hRows (0 : Fin 4)
    row1B := by
      have := hRows (1 : Fin 4)
      simpa [HgcdRecursiveFinishResult.toResult, hgcdMatLenRaw,
        hgcdMatLen] using (show
          hgcdMatLen result.matrix result.valid (1 : Fin 4) + result.lenB ≤
            outerLength + 1 by omega)
    row2A := by simpa [HgcdRecursiveFinishResult.toResult, hgcdMatLenRaw,
      hgcdMatLen] using hRows (2 : Fin 4)
    row3B := by
      have := hRows (3 : Fin 4)
      simpa [HgcdRecursiveFinishResult.toResult, hgcdMatLenRaw,
        hgcdMatLen] using (show
          hgcdMatLen result.matrix result.valid (3 : Fin 4) + result.lenB ≤
            outerLength + 1 by omega)
    row1A := by simpa [HgcdRecursiveFinishResult.toResult, hgcdMatLenRaw,
      hgcdMatLen] using hRows (1 : Fin 4)
    row3A := by simpa [HgcdRecursiveFinishResult.toResult, hgcdMatLenRaw,
      hgcdMatLen] using hRows (3 : Fin 4)
    inputBound := hOperands.inputBoundA
    stopped := hOperands.stopped
    positive := hOperands.positiveA
    aboveHalf := hOperands.aboveHalf
    coeffBound := by
      intro i
      simpa [HgcdRecursiveFinishResult.toResult, hgcdMatLenRaw,
        hgcdMatLen] using hCoeff i }

/-- Purely physical obligations for the exact final matrix block.  Besides
the two existing generated-call workspaces, the frame fields state that the
quotient update does not alter any buffer of the left matrix `R`. -/
structure HgcdRecursiveCombineMatrixWorkspace (this : DenseUPolyZp)
    (R S : HgcdMat) (hR : R.Valid) (hS : S.Valid)
    (q : RawPtr UInt64) (lenQ : Nat) (T a2 scratch : RawPtr UInt64)
    (heap : RawHeap) (modified : HgcdMatQuotientResult) : Prop where
  quotient : HgcdMatApplyQuotientWorkspaceProvider this S hS q lenQ T
    scratch heap
  rightLayout : RawHeap.SameLayout heap modified.heap
  rightPrefix : ∀ i : Fin 4, SameU64Prefix heap modified.heap
    (hgcdMatPtr R hR i) (hgcdMatLen R hR i)
  multiply : HgcdMatMulLoopWorkspaceProvider this R modified.matrix hR
    modified.valid a2 scratch

def HgcdRecursiveCombineMatrixWorkspaceProvider (this : DenseUPolyZp)
    (R S : HgcdMat) (hR : R.Valid) (hS : S.Valid)
    (q : RawPtr UInt64) (lenQ : Nat) (T a2 scratch : RawPtr UInt64)
    (heap : RawHeap) : Prop :=
  ∀ modified,
    hgcdMatApplyQuotient this S hS q lenQ T scratch heap = .ok modified →
    HgcdRecursiveCombineMatrixWorkspace this R S hR hS q lenQ T a2 scratch
      heap modified

/-- End-to-end semantic refinement of the exact final C++ matrix block.
The proof consumes the real quotient update and the complete four-entry
matrix multiplication; the L2 result is obtained only after both generated
executions have succeeded. -/
theorem hgcdRecursiveCombineMatrix_refines (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (M R S : HgcdMat) (hM : M.Valid) (hR : R.Valid) (hS : S.Valid)
    (q : RawPtr UInt64) (lenQ : Nat) (T a2 scratch : RawPtr UInt64)
    (heap : RawHeap) (result : HgcdMatMulResult)
    (right entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (quotient : Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (physical : HgcdRecursiveCombineMatrixWorkspaceProvider this R S hR hS
      q lenQ T a2 scratch heap)
    (hRight : HgcdMatRawDenseRep this heap R right hR)
    (hSRep : HgcdMatRawDenseRep this heap S entries hS)
    (hQ : RawDensePolyRep this heap q lenQ quotient)
    (hrun : hgcdRecursiveCombineMatrix this M R S hM hR hS q lenQ T a2
      scratch heap = .ok result) :
    HgcdMatRawDenseRep this result.heap R right hR ∧
      ∃ hResult : result.matrix.Valid,
        HgcdMatRawDenseRep this result.heap result.matrix
          (hgcdMatProductEntry right
            (hgcdMatApplyQuotientEntries entries quotient)) hResult := by
  rcases hgcdRecursiveCombineMatrix_exec this M R S hM hR hS q lenQ T a2
      scratch heap result hrun with ⟨modified, hmodified, hmul⟩
  have hwork := physical modified hmodified
  have hModified := hgcdMatApplyQuotient_refines this S hS q lenQ T scratch
    heap modified entries quotient hcfg hp hwork.quotient hQ hSRep hmodified
  have hRightModified :
      HgcdMatRawDenseRep this modified.heap R right hR := by
    intro i
    exact rawDensePolyRep_of_same_prefix this heap modified.heap
      (hgcdMatPtr R hR i) (hgcdMatLen R hR i) (right i)
      hwork.rightLayout (hwork.rightPrefix i) (hRight i)
  have hProduct := hgcdMatMul_refines this M R modified.matrix hM hR
    modified.valid a2 scratch modified.heap result right
    (hgcdMatApplyQuotientEntries entries quotient) hcfg hp hwork.multiply
    hRightModified hModified.1 hmul
  exact ⟨hProduct.1, hProduct.2.2⟩


end CLPoly.Impl.StrictHGCDRawRefinement
