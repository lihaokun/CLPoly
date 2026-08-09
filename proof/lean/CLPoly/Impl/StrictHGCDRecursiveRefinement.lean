import CLPoly.Impl.StrictHGCDReconstructRefinement

set_option autoImplicit false

namespace CLPoly.Impl.StrictHGCDRawRefinement

open Generated.StrictHGCD
open CLPoly.Impl.RawPolynomialRep
open CLPoly.Impl.StrictDivremRefinement
open CLPoly.Impl.StrictEuclidRefinement
open CLPoly.Impl.StrictPolyAddSubRefinement
open CLPoly.Impl.StrictMulRefinement
open CLPoly.Impl.StrictWordArithmetic

/-- Length evidence for the same real quotient-update/matrix-product block.
The intermediate matrix is the concrete result returned by the generated
quotient execution, so the final product bounds remain tied to actual C++
descriptors rather than a specification-side matrix. -/
theorem hgcdRecursiveCombineMatrix_length_bounds (this : DenseUPolyZp)
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
    ∃ (modified : HgcdMatQuotientResult) (hResult : result.matrix.Valid),
      HgcdMatApplyQuotientLengthBounds S modified.matrix hS modified.valid
        lenQ ∧
      (∀ i : Fin 4,
        hgcdMatLen result.matrix hResult i ≤
          max
            (hgcdMatLen R hR ⟨2 * (i.val / 2), by omega⟩ +
              hgcdMatLen modified.matrix modified.valid
                ⟨i.val % 2, by omega⟩ - 1)
            (hgcdMatLen R hR ⟨2 * (i.val / 2) + 1, by omega⟩ +
              hgcdMatLen modified.matrix modified.valid
                ⟨2 + i.val % 2, by omega⟩ - 1)) := by
  rcases hgcdRecursiveCombineMatrix_exec this M R S hM hR hS q lenQ T a2
      scratch heap result hrun with ⟨modified, hmodified, hmul⟩
  have hwork := physical modified hmodified
  have hModified := hgcdMatApplyQuotient_refines this S hS q lenQ T scratch
    heap modified entries quotient hcfg hp hwork.quotient hQ hSRep hmodified
  have hRightModified : HgcdMatRawDenseRep this modified.heap R right hR := by
    intro i
    exact rawDensePolyRep_of_same_prefix this heap modified.heap
      (hgcdMatPtr R hR i) (hgcdMatLen R hR i) (right i)
      hwork.rightLayout (hwork.rightPrefix i) (hRight i)
  have hProduct := hgcdMatMul_refines this M R modified.matrix hM hR
    modified.valid a2 scratch modified.heap result right
    (hgcdMatApplyQuotientEntries entries quotient) hcfg hp hwork.multiply
    hRightModified hModified.1 hmul
  rcases hProduct.2.2 with ⟨hResult, hResultRep⟩
  refine ⟨modified, hResult, ?_, ?_⟩
  · exact hgcdMatApplyQuotientEntries_length_bounds this heap modified.heap
      modified.heap S modified.matrix hS modified.valid q lenQ entries
      quotient hSRep hModified.2 hModified.1
  · intro i
    exact hgcdMatProductEntry_length_le this result.heap R modified.matrix
      result.matrix hR modified.valid hResult right
      (hgcdMatApplyQuotientEntries entries quotient) hProduct.1 hProduct.2.1
      hResultRep i

/-- The source-selected second suffix has half-capacity no larger than the
excess of the reconstructed divisor over the outer half split. -/
theorem hgcdSecondInput_halfCapacity_le
    (m reconstructedLenB k secondInputLength : Nat)
    (hreconstructedLower : m + 1 ≤ reconstructedLenB)
    (hk : k = 2 * m - reconstructedLenB + 1)
    (hc : secondInputLength = reconstructedLenB - k) :
    secondInputLength - secondInputLength / 2 ≤ reconstructedLenB - m := by
  omega

/-- A first-matrix row paired with its returned A descriptor fits both the
direct second coefficient and the quotient-updated coefficient products. -/
theorem hgcdCombinedRowTerm_le
    (capacity firstLenA reconstructedLenB m lenQ r s : Nat)
    (hrA : r + firstLenA ≤ capacity + 1)
    (hreconstructedOrder : reconstructedLenB ≤ m + firstLenA)
    (hreconstructedLower : m + 1 ≤ reconstructedLenB)
    (hlenQ : lenQ ≤ m + firstLenA - (reconstructedLenB - 1))
    (hs : s ≤ reconstructedLenB - m) :
    r + s - 1 ≤ capacity ∧
      r + (lenQ + s - 1) - 1 ≤ capacity := by
  have hqsum : lenQ + reconstructedLenB ≤ m + firstLenA + 1 := by
    omega
  have hssum : s + m ≤ reconstructedLenB := by
    omega
  constructor <;> omega

/-- Both entries of one column of the exact final matrix product remain
within the outer half-length capacity.  The statement is column-parametric:
the source quotient update has the same descriptor shape for columns zero
and one, so this lemma covers all four matrix entries without duplicating
the arithmetic argument. -/
theorem hgcdRecursiveCombinedColumn_coeff_bounds
    (outerLength m highLength firstLenA reconstructedLenA reconstructedLenB
      lenQ k secondInputLength r0 r1 r2 r3 s1 s3 u1 u3 out1 out3 : Nat)
    (hm : m = outerLength / 2)
    (hhigh : highLength = outerLength - m)
    (hr0A : r0 + firstLenA ≤ highLength + 1)
    (hr1A : r1 + firstLenA ≤ highLength + 1)
    (hr2A : r2 + firstLenA ≤ highLength + 1)
    (hr3A : r3 + firstLenA ≤ highLength + 1)
    (hreconstructedA : reconstructedLenA = m + firstLenA)
    (hreconstructedOrder : reconstructedLenB ≤ reconstructedLenA)
    (hreconstructedLower : m + 1 ≤ reconstructedLenB)
    (hlenQ : lenQ ≤ reconstructedLenA - (reconstructedLenB - 1))
    (hk : k = 2 * m - reconstructedLenB + 1)
    (hc : secondInputLength = reconstructedLenB - k)
    (hs1 : s1 ≤ secondInputLength - secondInputLength / 2)
    (hs3 : s3 ≤ secondInputLength - secondInputLength / 2)
    (hu1 : u1 ≤ max s3 (lenQ + s1 - 1))
    (hu3 : u3 = s1)
    (hout1 : out1 ≤ max (r0 + u1 - 1) (r1 + u3 - 1))
    (hout3 : out3 ≤ max (r2 + u1 - 1) (r3 + u3 - 1)) :
    out1 ≤ outerLength - outerLength / 2 ∧
      out3 ≤ outerLength - outerLength / 2 := by
  have htarget : outerLength - outerLength / 2 = highLength := by
    rw [← hm]
    exact hhigh.symm
  have hcap := hgcdSecondInput_halfCapacity_le m reconstructedLenB k
    secondInputLength hreconstructedLower hk hc
  have hs1' : s1 ≤ reconstructedLenB - m := hs1.trans hcap
  have hs3' : s3 ≤ reconstructedLenB - m := hs3.trans hcap
  have horder' : reconstructedLenB ≤ m + firstLenA := by
    rwa [← hreconstructedA]
  have hq' : lenQ ≤ m + firstLenA - (reconstructedLenB - 1) := by
    rwa [← hreconstructedA]
  have boundU (r : Nat)
      (hrA : r + firstLenA ≤ highLength + 1) :
      r + u1 - 1 ≤ outerLength - outerLength / 2 := by
    rw [htarget]
    have hdirect := (hgcdCombinedRowTerm_le highLength firstLenA
      reconstructedLenB m lenQ r s3 hrA horder' (by omega) hq' hs3').1
    have hquotient := (hgcdCombinedRowTerm_le highLength firstLenA
      reconstructedLenB m lenQ r s1 hrA horder' (by omega) hq' hs1').2
    rcases le_total s3 (lenQ + s1 - 1) with hle | hle
    · have hu : u1 ≤ lenQ + s1 - 1 := by
        simpa [max_eq_right hle] using hu1
      exact (Nat.sub_le_sub_right (Nat.add_le_add_left hu r) 1).trans
        hquotient
    · have hu : u1 ≤ s3 := by
        simpa [max_eq_left hle] using hu1
      exact (Nat.sub_le_sub_right (Nat.add_le_add_left hu r) 1).trans
        hdirect
  have boundOdd (r : Nat)
      (hrA : r + firstLenA ≤ highLength + 1) :
      r + u3 - 1 ≤ outerLength - outerLength / 2 := by
    rw [htarget, hu3]
    exact (hgcdCombinedRowTerm_le highLength firstLenA reconstructedLenB m
      lenQ r s1 hrA horder' (by omega) hq' hs1').1
  constructor
  · exact hout1.trans (max_le (boundU r0 hr0A) (boundOdd r1 hr1A))
  · exact hout3.trans (max_le (boundU r2 hr2A) (boundOdd r3 hr3A))

/-- Entries 1 and 3 of the exact final matrix product remain within the
outer half-length capacity.  This specialization records the generated
matrix layout used by the recursive HGCD implementation. -/
theorem hgcdRecursiveCombined_odd_coeff_bounds
    (outerLength m highLength firstLenA reconstructedLenA reconstructedLenB
      lenQ k secondInputLength r0 r1 r2 r3 s1 s3 u1 u3 out1 out3 : Nat)
    (hm : m = outerLength / 2)
    (hhigh : highLength = outerLength - m)
    (hr0A : r0 + firstLenA ≤ highLength + 1)
    (hr1A : r1 + firstLenA ≤ highLength + 1)
    (hr2A : r2 + firstLenA ≤ highLength + 1)
    (hr3A : r3 + firstLenA ≤ highLength + 1)
    (hreconstructedA : reconstructedLenA = m + firstLenA)
    (hreconstructedOrder : reconstructedLenB ≤ reconstructedLenA)
    (hreconstructedLower : m + 1 ≤ reconstructedLenB)
    (hlenQ : lenQ ≤ reconstructedLenA - (reconstructedLenB - 1))
    (hk : k = 2 * m - reconstructedLenB + 1)
    (hc : secondInputLength = reconstructedLenB - k)
    (hs1 : s1 ≤ secondInputLength - secondInputLength / 2)
    (hs3 : s3 ≤ secondInputLength - secondInputLength / 2)
    (hu1 : u1 ≤ max s3 (lenQ + s1 - 1))
    (hu3 : u3 = s1)
    (hout1 : out1 ≤ max (r0 + u1 - 1) (r1 + u3 - 1))
    (hout3 : out3 ≤ max (r2 + u1 - 1) (r3 + u3 - 1)) :
    out1 ≤ outerLength - outerLength / 2 ∧
      out3 ≤ outerLength - outerLength / 2 := by
  exact hgcdRecursiveCombinedColumn_coeff_bounds outerLength m highLength
    firstLenA reconstructedLenA reconstructedLenB lenQ k secondInputLength
    r0 r1 r2 r3 s1 s3 u1 u3 out1 out3 hm hhigh hr0A hr1A hr2A hr3A
    hreconstructedA hreconstructedOrder hreconstructedLower hlenQ hk hc
    hs1 hs3 hu1 hu3 hout1 hout3

/-- The two column instances together give the uniform coefficient bound
required by the recursive HGCD contract for every concrete output entry. -/
theorem hgcdRecursiveCombined_all_coeff_bounds
    (outerLength m highLength firstLenA reconstructedLenA reconstructedLenB
      lenQ k secondInputLength r0 r1 r2 r3 s0 s1 s2 s3 u0 u1 u2 u3
      out0 out1 out2 out3 : Nat)
    (hm : m = outerLength / 2)
    (hhigh : highLength = outerLength - m)
    (hr0A : r0 + firstLenA ≤ highLength + 1)
    (hr1A : r1 + firstLenA ≤ highLength + 1)
    (hr2A : r2 + firstLenA ≤ highLength + 1)
    (hr3A : r3 + firstLenA ≤ highLength + 1)
    (hreconstructedA : reconstructedLenA = m + firstLenA)
    (hreconstructedOrder : reconstructedLenB ≤ reconstructedLenA)
    (hreconstructedLower : m + 1 ≤ reconstructedLenB)
    (hlenQ : lenQ ≤ reconstructedLenA - (reconstructedLenB - 1))
    (hk : k = 2 * m - reconstructedLenB + 1)
    (hc : secondInputLength = reconstructedLenB - k)
    (hs0 : s0 ≤ secondInputLength - secondInputLength / 2)
    (hs1 : s1 ≤ secondInputLength - secondInputLength / 2)
    (hs2 : s2 ≤ secondInputLength - secondInputLength / 2)
    (hs3 : s3 ≤ secondInputLength - secondInputLength / 2)
    (hu0 : u0 ≤ max s2 (lenQ + s0 - 1))
    (hu1 : u1 ≤ max s3 (lenQ + s1 - 1))
    (hu2 : u2 = s0)
    (hu3 : u3 = s1)
    (hout0 : out0 ≤ max (r0 + u0 - 1) (r1 + u2 - 1))
    (hout1 : out1 ≤ max (r0 + u1 - 1) (r1 + u3 - 1))
    (hout2 : out2 ≤ max (r2 + u0 - 1) (r3 + u2 - 1))
    (hout3 : out3 ≤ max (r2 + u1 - 1) (r3 + u3 - 1)) :
    ∀ i : Fin 4,
      [out0, out1, out2, out3][i] ≤ outerLength - outerLength / 2 := by
  have heven := hgcdRecursiveCombinedColumn_coeff_bounds outerLength m
    highLength firstLenA reconstructedLenA reconstructedLenB lenQ k
    secondInputLength r0 r1 r2 r3 s0 s2 u0 u2 out0 out2 hm hhigh hr0A
    hr1A hr2A hr3A hreconstructedA hreconstructedOrder hreconstructedLower
    hlenQ hk hc hs0 hs2 hu0 hu2 hout0 hout2
  have hodd := hgcdRecursiveCombinedColumn_coeff_bounds outerLength m
    highLength firstLenA reconstructedLenA reconstructedLenB lenQ k
    secondInputLength r0 r1 r2 r3 s1 s3 u1 u3 out1 out3 hm hhigh hr0A
    hr1A hr2A hr3A hreconstructedA hreconstructedOrder hreconstructedLower
    hlenQ hk hc hs1 hs3 hu1 hu3 hout1 hout3
  intro i
  fin_cases i <;> simp_all

/-- Sharp row/A arithmetic for one product term of the final matrix.  Unlike
the uniform half-capacity bound, this retains the exact second-child A
length and is therefore strong enough for the recursive row/A contract. -/
theorem hgcdCombinedRowTerm_add_finalA_le
    (outerLength highLength m firstLenA reconstructedLenB lenQ k
      secondInputLength secondLenA r s : Nat)
    (houter : highLength + m = outerLength)
    (hrA : r + firstLenA ≤ highLength + 1)
    (hfirstBound : firstLenA ≤ highLength)
    (hreconstructedOrder : reconstructedLenB ≤ m + firstLenA)
    (hreconstructedLower : m + 1 ≤ reconstructedLenB)
    (hlenQ : lenQ ≤ m + firstLenA - (reconstructedLenB - 1))
    (hsplit : k + secondInputLength = reconstructedLenB)
    (hsA : s + secondLenA ≤ secondInputLength + 1) :
    (r + s - 1) + (k + secondLenA) ≤ outerLength + 1 ∧
      (r + (lenQ + s - 1) - 1) + (k + secondLenA) ≤
        outerLength + 1 := by
  have hqsum : lenQ + reconstructedLenB ≤ m + firstLenA + 1 := by
    omega
  have hdirectSum : r + s + k + secondLenA ≤ outerLength + 2 := by
    omega
  have hquotientSum : r + lenQ + s + k + secondLenA ≤
      outerLength + 3 := by
    omega
  constructor
  · by_cases hsum : r + s = 0
    · simp [hsum]
      omega
    · omega
  · by_cases hinner : lenQ + s = 0
    · simp [hinner]
      omega
    · by_cases hsum : r + (lenQ + s - 1) = 0
      · simp [hsum]
        omega
      · omega

/-- Both rows of one concrete final-product column retain the sharp pairing
with the final reconstructed A descriptor. -/
theorem hgcdRecursiveCombinedColumn_rowA_bounds
    (outerLength highLength m firstLenA reconstructedLenB lenQ k
      secondInputLength secondLenA finalLenA r0 r1 r2 r3 sTop sBottom
      uTop uBottom outTop outBottom : Nat)
    (houter : highLength + m = outerLength)
    (hr0A : r0 + firstLenA ≤ highLength + 1)
    (hr1A : r1 + firstLenA ≤ highLength + 1)
    (hr2A : r2 + firstLenA ≤ highLength + 1)
    (hr3A : r3 + firstLenA ≤ highLength + 1)
    (hfirstBound : firstLenA ≤ highLength)
    (hreconstructedOrder : reconstructedLenB ≤ m + firstLenA)
    (hreconstructedLower : m + 1 ≤ reconstructedLenB)
    (hlenQ : lenQ ≤ m + firstLenA - (reconstructedLenB - 1))
    (hsplit : k + secondInputLength = reconstructedLenB)
    (hsTopA : sTop + secondLenA ≤ secondInputLength + 1)
    (hsBottomA : sBottom + secondLenA ≤ secondInputLength + 1)
    (hfinalA : finalLenA = k + secondLenA)
    (huTop : uTop ≤ max sBottom (lenQ + sTop - 1))
    (huBottom : uBottom = sTop)
    (houtTop : outTop ≤ max (r0 + uTop - 1) (r1 + uBottom - 1))
    (houtBottom : outBottom ≤
      max (r2 + uTop - 1) (r3 + uBottom - 1)) :
    outTop + finalLenA ≤ outerLength + 1 ∧
      outBottom + finalLenA ≤ outerLength + 1 := by
  have boundUpdated (r : Nat) (hrA : r + firstLenA ≤ highLength + 1) :
      (r + uTop - 1) + finalLenA ≤ outerLength + 1 := by
    rw [hfinalA]
    have hdirect := (hgcdCombinedRowTerm_add_finalA_le outerLength
      highLength m firstLenA reconstructedLenB lenQ k secondInputLength
      secondLenA r sBottom houter hrA hfirstBound hreconstructedOrder
      hreconstructedLower hlenQ hsplit hsBottomA).1
    have hquotient := (hgcdCombinedRowTerm_add_finalA_le outerLength
      highLength m firstLenA reconstructedLenB lenQ k secondInputLength
      secondLenA r sTop houter hrA hfirstBound hreconstructedOrder
      hreconstructedLower
      hlenQ hsplit hsTopA).2
    rcases le_total sBottom (lenQ + sTop - 1) with hle | hle
    · have hu : uTop ≤ lenQ + sTop - 1 := by
        simpa [max_eq_right hle] using huTop
      exact (Nat.add_le_add_right
        (Nat.sub_le_sub_right (Nat.add_le_add_left hu r) 1)
        (k + secondLenA)).trans hquotient
    · have hu : uTop ≤ sBottom := by
        simpa [max_eq_left hle] using huTop
      exact (Nat.add_le_add_right
        (Nat.sub_le_sub_right (Nat.add_le_add_left hu r) 1)
        (k + secondLenA)).trans hdirect
  have boundDirect (r : Nat) (hrA : r + firstLenA ≤ highLength + 1) :
      (r + uBottom - 1) + finalLenA ≤ outerLength + 1 := by
    rw [hfinalA, huBottom]
    exact (hgcdCombinedRowTerm_add_finalA_le outerLength highLength m
      firstLenA reconstructedLenB lenQ k secondInputLength secondLenA r
      sTop houter hrA hfirstBound hreconstructedOrder hreconstructedLower
      hlenQ hsplit hsTopA).1
  have boundOutput (out left right : Nat)
      (hout : out ≤ max left right)
      (hleft : left + finalLenA ≤ outerLength + 1)
      (hright : right + finalLenA ≤ outerLength + 1) :
      out + finalLenA ≤ outerLength + 1 := by
    rcases le_total left right with hle | hle
    · have ho : out ≤ right := by simpa [max_eq_right hle] using hout
      exact (Nat.add_le_add_right ho finalLenA).trans hright
    · have ho : out ≤ left := by simpa [max_eq_left hle] using hout
      exact (Nat.add_le_add_right ho finalLenA).trans hleft
  constructor
  · exact boundOutput outTop (r0 + uTop - 1) (r1 + uBottom - 1)
      houtTop (boundUpdated r0 hr0A) (boundDirect r1 hr1A)
  · exact boundOutput outBottom (r2 + uTop - 1) (r3 + uBottom - 1)
      houtBottom (boundUpdated r2 hr2A) (boundDirect r3 hr3A)

/-- Applying the sharp column argument twice closes the row/A pairing for
all four entries of the final matrix. -/
theorem hgcdRecursiveCombined_all_rowA_bounds
    (outerLength highLength m firstLenA reconstructedLenB lenQ k
      secondInputLength secondLenA finalLenA r0 r1 r2 r3 s0 s1 s2 s3
      u0 u1 u2 u3 out0 out1 out2 out3 : Nat)
    (houter : highLength + m = outerLength)
    (hr0A : r0 + firstLenA ≤ highLength + 1)
    (hr1A : r1 + firstLenA ≤ highLength + 1)
    (hr2A : r2 + firstLenA ≤ highLength + 1)
    (hr3A : r3 + firstLenA ≤ highLength + 1)
    (hfirstBound : firstLenA ≤ highLength)
    (hreconstructedOrder : reconstructedLenB ≤ m + firstLenA)
    (hreconstructedLower : m + 1 ≤ reconstructedLenB)
    (hlenQ : lenQ ≤ m + firstLenA - (reconstructedLenB - 1))
    (hsplit : k + secondInputLength = reconstructedLenB)
    (hs0A : s0 + secondLenA ≤ secondInputLength + 1)
    (hs1A : s1 + secondLenA ≤ secondInputLength + 1)
    (hs2A : s2 + secondLenA ≤ secondInputLength + 1)
    (hs3A : s3 + secondLenA ≤ secondInputLength + 1)
    (hfinalA : finalLenA = k + secondLenA)
    (hu0 : u0 ≤ max s2 (lenQ + s0 - 1))
    (hu1 : u1 ≤ max s3 (lenQ + s1 - 1))
    (hu2 : u2 = s0) (hu3 : u3 = s1)
    (hout0 : out0 ≤ max (r0 + u0 - 1) (r1 + u2 - 1))
    (hout1 : out1 ≤ max (r0 + u1 - 1) (r1 + u3 - 1))
    (hout2 : out2 ≤ max (r2 + u0 - 1) (r3 + u2 - 1))
    (hout3 : out3 ≤ max (r2 + u1 - 1) (r3 + u3 - 1)) :
    ∀ i : Fin 4,
      [out0, out1, out2, out3][i] + finalLenA ≤ outerLength + 1 := by
  have heven := hgcdRecursiveCombinedColumn_rowA_bounds outerLength
    highLength m firstLenA reconstructedLenB lenQ k secondInputLength
    secondLenA finalLenA r0 r1 r2 r3 s0 s2 u0 u2 out0 out2 houter hr0A
    hr1A hr2A hr3A hfirstBound hreconstructedOrder hreconstructedLower hlenQ
    hsplit hs0A hs2A hfinalA hu0 hu2 hout0 hout2
  have hodd := hgcdRecursiveCombinedColumn_rowA_bounds outerLength
    highLength m firstLenA reconstructedLenB lenQ k secondInputLength
    secondLenA finalLenA r0 r1 r2 r3 s1 s3 u1 u3 out1 out3 houter hr0A
    hr1A hr2A hr3A hfirstBound hreconstructedOrder hreconstructedLower hlenQ
    hsplit hs1A hs3A hfinalA hu1 hu3 hout1 hout3
  intro i
  fin_cases i <;> simp_all

/-- The uniform bound above is realized by the descriptors returned from the
actual generated quotient-update/matrix-product tail.  In particular, the
intermediate `modified` matrix is obtained by executing
`hgcdRecursiveCombineMatrix`; it is not chosen by the specification. -/
theorem hgcdRecursiveCombineMatrix_coeff_bounds (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (M R S : HgcdMat) (hM : M.Valid) (hR : R.Valid) (hS : S.Valid)
    (q : RawPtr UInt64) (lenQ : Nat) (T a2 scratch : RawPtr UInt64)
    (heap : RawHeap) (result : HgcdMatMulResult)
    (right entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (quotient : Polynomial (ZMod this._p.toNat))
    (outerLength m highLength firstLenA reconstructedLenA reconstructedLenB
      k secondInputLength : Nat)
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (physical : HgcdRecursiveCombineMatrixWorkspaceProvider this R S hR hS
      q lenQ T a2 scratch heap)
    (hRight : HgcdMatRawDenseRep this heap R right hR)
    (hSRep : HgcdMatRawDenseRep this heap S entries hS)
    (hQ : RawDensePolyRep this heap q lenQ quotient)
    (hm : m = outerLength / 2)
    (hhigh : highLength = outerLength - m)
    (hRRows : ∀ i : Fin 4,
      hgcdMatLen R hR i + firstLenA ≤ highLength + 1)
    (hreconstructedA : reconstructedLenA = m + firstLenA)
    (hreconstructedOrder : reconstructedLenB ≤ reconstructedLenA)
    (hreconstructedLower : m + 1 ≤ reconstructedLenB)
    (hlenQ : lenQ ≤ reconstructedLenA - (reconstructedLenB - 1))
    (hk : k = 2 * m - reconstructedLenB + 1)
    (hc : secondInputLength = reconstructedLenB - k)
    (hSBound : ∀ i : Fin 4,
      hgcdMatLen S hS i ≤ secondInputLength - secondInputLength / 2)
    (hrun : hgcdRecursiveCombineMatrix this M R S hM hR hS q lenQ T a2
      scratch heap = .ok result) :
    ∃ hResult : result.matrix.Valid,
      ∀ i : Fin 4,
        hgcdMatLen result.matrix hResult i ≤
          outerLength - outerLength / 2 := by
  rcases hgcdRecursiveCombineMatrix_length_bounds this M R S hM hR hS q
      lenQ T a2 scratch heap result right entries quotient hcfg hp physical
      hRight hSRep hQ hrun with ⟨modified, hResult, hModified, hOutput⟩
  refine ⟨hResult, ?_⟩
  have hAll := hgcdRecursiveCombined_all_coeff_bounds outerLength m highLength
    firstLenA reconstructedLenA reconstructedLenB lenQ k secondInputLength
    (hgcdMatLen R hR 0) (hgcdMatLen R hR 1)
    (hgcdMatLen R hR 2) (hgcdMatLen R hR 3)
    (hgcdMatLen S hS 0) (hgcdMatLen S hS 1)
    (hgcdMatLen S hS 2) (hgcdMatLen S hS 3)
    (hgcdMatLen modified.matrix modified.valid 0)
    (hgcdMatLen modified.matrix modified.valid 1)
    (hgcdMatLen modified.matrix modified.valid 2)
    (hgcdMatLen modified.matrix modified.valid 3)
    (hgcdMatLen result.matrix hResult 0)
    (hgcdMatLen result.matrix hResult 1)
    (hgcdMatLen result.matrix hResult 2)
    (hgcdMatLen result.matrix hResult 3) hm hhigh (hRRows 0) (hRRows 1)
    (hRRows 2) (hRRows 3) hreconstructedA hreconstructedOrder
    hreconstructedLower hlenQ hk hc (hSBound 0) (hSBound 1) (hSBound 2)
    (hSBound 3) hModified.row0 hModified.row1 hModified.row2
    hModified.row3 (by simpa using hOutput (0 : Fin 4))
    (by simpa using hOutput (1 : Fin 4))
    (by simpa using hOutput (2 : Fin 4))
    (by simpa using hOutput (3 : Fin 4))
  intro i
  fin_cases i
  · simpa using hAll (0 : Fin 4)
  · simpa using hAll (1 : Fin 4)
  · simpa using hAll (2 : Fin 4)
  · simpa using hAll (3 : Fin 4)

/-- Execution-level sharp row/A bounds for the same generated final matrix
block.  This is the descriptor half of the complete recursive length
invariant; no specification matrix appears in the statement. -/
theorem hgcdRecursiveCombineMatrix_rowA_bounds (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (M R S : HgcdMat) (hM : M.Valid) (hR : R.Valid) (hS : S.Valid)
    (q : RawPtr UInt64) (lenQ : Nat) (T a2 scratch : RawPtr UInt64)
    (heap : RawHeap) (result : HgcdMatMulResult)
    (right entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (quotient : Polynomial (ZMod this._p.toNat))
    (outerLength highLength m firstLenA reconstructedLenB k
      secondInputLength secondLenA finalLenA : Nat)
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (physical : HgcdRecursiveCombineMatrixWorkspaceProvider this R S hR hS
      q lenQ T a2 scratch heap)
    (hRight : HgcdMatRawDenseRep this heap R right hR)
    (hSRep : HgcdMatRawDenseRep this heap S entries hS)
    (hQ : RawDensePolyRep this heap q lenQ quotient)
    (houter : highLength + m = outerLength)
    (hRRows : ∀ i : Fin 4,
      hgcdMatLen R hR i + firstLenA ≤ highLength + 1)
    (hfirstBound : firstLenA ≤ highLength)
    (hreconstructedOrder : reconstructedLenB ≤ m + firstLenA)
    (hreconstructedLower : m + 1 ≤ reconstructedLenB)
    (hlenQ : lenQ ≤ m + firstLenA - (reconstructedLenB - 1))
    (hsplit : k + secondInputLength = reconstructedLenB)
    (hSRows : ∀ i : Fin 4,
      hgcdMatLen S hS i + secondLenA ≤ secondInputLength + 1)
    (hfinalA : finalLenA = k + secondLenA)
    (hrun : hgcdRecursiveCombineMatrix this M R S hM hR hS q lenQ T a2
      scratch heap = .ok result) :
    ∃ hResult : result.matrix.Valid,
      ∀ i : Fin 4,
        hgcdMatLen result.matrix hResult i + finalLenA ≤ outerLength + 1 := by
  rcases hgcdRecursiveCombineMatrix_length_bounds this M R S hM hR hS q
      lenQ T a2 scratch heap result right entries quotient hcfg hp physical
      hRight hSRep hQ hrun with ⟨modified, hResult, hModified, hOutput⟩
  refine ⟨hResult, ?_⟩
  have hAll := hgcdRecursiveCombined_all_rowA_bounds outerLength highLength m
    firstLenA reconstructedLenB lenQ k secondInputLength secondLenA finalLenA
    (hgcdMatLen R hR 0) (hgcdMatLen R hR 1)
    (hgcdMatLen R hR 2) (hgcdMatLen R hR 3)
    (hgcdMatLen S hS 0) (hgcdMatLen S hS 1)
    (hgcdMatLen S hS 2) (hgcdMatLen S hS 3)
    (hgcdMatLen modified.matrix modified.valid 0)
    (hgcdMatLen modified.matrix modified.valid 1)
    (hgcdMatLen modified.matrix modified.valid 2)
    (hgcdMatLen modified.matrix modified.valid 3)
    (hgcdMatLen result.matrix hResult 0)
    (hgcdMatLen result.matrix hResult 1)
    (hgcdMatLen result.matrix hResult 2)
    (hgcdMatLen result.matrix hResult 3) houter (hRRows 0) (hRRows 1)
    (hRRows 2) (hRRows 3) hfirstBound hreconstructedOrder
    hreconstructedLower hlenQ hsplit (hSRows 0) (hSRows 1) (hSRows 2)
    (hSRows 3) hfinalA hModified.row0 hModified.row1 hModified.row2
    hModified.row3 (by simpa using hOutput (0 : Fin 4))
    (by simpa using hOutput (1 : Fin 4))
    (by simpa using hOutput (2 : Fin 4))
    (by simpa using hOutput (3 : Fin 4))
  intro i
  fin_cases i
  · simpa using hAll (0 : Fin 4)
  · simpa using hAll (1 : Fin 4)
  · simpa using hAll (2 : Fin 4)
  · simpa using hAll (3 : Fin 4)

/-- Physical obligations that connect the real final reconstruction to the
optional quotient/matrix-product block while framing both reconstructed
output polynomials. -/
structure HgcdRecursiveFinishWorkspace (this : DenseUPolyZp)
    (M R S : HgcdMat) (hM : M.Valid) (hR : R.Valid) (hS : S.Valid)
    (A B T0 lowA lowB highA highB q : RawPtr UInt64)
    (lenLowA lenLowB lenHighA lenHighB shift lenQ : Nat)
    (a2 scratch : RawPtr UInt64) (sgnS : Int) (heap : RawHeap)
    (reconstructed : HgcdRecursiveReconstructPairResult) : Prop where
  reconstruct : HgcdRecursiveReconstructPairWorkspaceProvider this A B T0
    lowA lowB highA highB scratch lenLowA lenLowB lenHighA lenHighB shift
    S hS sgnS heap
  afterLayout : RawHeap.SameLayout heap reconstructed.heap
  rightPrefix : ∀ i : Fin 4, SameU64Prefix heap reconstructed.heap
    (hgcdMatPtr R hR i) (hgcdMatLen R hR i)
  secondPrefix : ∀ i : Fin 4, SameU64Prefix heap reconstructed.heap
    (hgcdMatPtr S hS i) (hgcdMatLen S hS i)
  quotientPrefix : SameU64Prefix heap reconstructed.heap q lenQ
  combine : ∀ (combined : HgcdMatMulResult),
    hgcdRecursiveCombineMatrix this M R S hM hR hS q lenQ T0 a2 scratch
        reconstructed.heap = .ok combined →
    HgcdRecursiveCombineMatrixWorkspaceProvider this R S hR hS q lenQ T0
      a2 scratch reconstructed.heap ∧
    RawHeap.SameLayout reconstructed.heap combined.heap ∧
    SameU64Prefix reconstructed.heap combined.heap A reconstructed.lenA ∧
    SameU64Prefix reconstructed.heap combined.heap B reconstructed.lenB

def HgcdRecursiveFinishWorkspaceProvider (this : DenseUPolyZp)
    (M R S : HgcdMat) (hM : M.Valid) (hR : R.Valid) (hS : S.Valid)
    (A B T0 lowA lowB highA highB q : RawPtr UInt64)
    (lenLowA lenLowB lenHighA lenHighB shift lenQ : Nat)
    (a2 scratch : RawPtr UInt64) (sgnS : Int) (heap : RawHeap) : Prop :=
  ∀ reconstructed,
    hgcdRecursiveReconstructPair this A B T0 lowA lowB highA highB scratch
      lenLowA lenLowB lenHighA lenHighB shift S hS sgnS heap =
        .ok reconstructed →
    HgcdRecursiveFinishWorkspace this M R S hM hR hS A B T0 lowA lowB highA
      highB q lenLowA lenLowB lenHighA lenHighB shift lenQ a2 scratch sgnS
      heap reconstructed

/-- Raw refinement of the exact final source tail.  Both reconstructed
outputs are obtained from the four generated calls and are explicitly
framed across the optional quotient update and full matrix multiplication. -/
theorem hgcdRecursiveFinish_refines (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (M R S : HgcdMat) (hM : M.Valid) (hR : R.Valid) (hS : S.Valid)
    (computeM : Bool) (A B T0 lowA lowB highA highB q : RawPtr UInt64)
    (lenLowA lenLowB lenHighA lenHighB shift lenQ : Nat)
    (a2 scratch : RawPtr UInt64) (sgnR sgnS : Int) (heap : RawHeap)
    (result : HgcdRecursiveFinishResult)
    (right entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (quotient polyLowA polyLowB polyHighA polyHighB :
      Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (physical : HgcdRecursiveFinishWorkspaceProvider this M R S hM hR hS A
      B T0 lowA lowB highA highB q lenLowA lenLowB lenHighA lenHighB shift
      lenQ a2 scratch sgnS heap)
    (hRRep : HgcdMatRawDenseRep this heap R right hR)
    (hSRep : HgcdMatRawDenseRep this heap S entries hS)
    (hQ : RawDensePolyRep this heap q lenQ quotient)
    (hLowA : RawCanonicalPolySlice this heap lowA lenLowA polyLowA)
    (hLowB : RawCanonicalPolySlice this heap lowB lenLowB polyLowB)
    (hHighA : RawDensePolyRep this heap highA lenHighA polyHighA)
    (hHighB : RawDensePolyRep this heap highB lenHighB polyHighB)
    (hrun : hgcdRecursiveFinish this M R S hM hR hS computeM A B T0 lowA
      lowB highA highB q lenLowA lenLowB lenHighA lenHighB shift lenQ a2
      scratch sgnR sgnS heap = .ok result) :
    RawDensePolyRep this result.heap A result.lenA
        (hgcdReconstructedLowA entries polyLowA polyLowB sgnS +
          Polynomial.X ^ shift * polyHighA) ∧
      RawDensePolyRep this result.heap B result.lenB
        (hgcdReconstructedLowB entries polyLowA polyLowB sgnS +
          Polynomial.X ^ shift * polyHighB) ∧
      result.sgn = -(sgnR * sgnS) ∧
      (computeM = true →
        HgcdMatRawDenseRep this result.heap result.matrix
          (hgcdMatProductEntry right
            (hgcdMatApplyQuotientEntries entries quotient)) result.valid) := by
  rcases hgcdRecursiveFinish_exec this M R S hM hR hS computeM A B T0 lowA
      lowB highA highB q lenLowA lenLowB lenHighA lenHighB shift lenQ a2
      scratch sgnR sgnS heap result hrun with
    ⟨reconstructed, hreconstruct, hlenA, hlenB, hsgn, htail⟩
  have hwork := physical reconstructed hreconstruct
  rcases hgcdRecursiveReconstructPair_refines this A B T0 lowA lowB highA
      highB scratch lenLowA lenLowB lenHighA lenHighB shift S hS sgnS heap
      reconstructed entries polyLowA polyLowB polyHighA polyHighB hcfg hp
      hwork.reconstruct hSRep hLowA hLowB hHighA hHighB hreconstruct with
    ⟨hAReconstructed, hBReconstructed, _, _, _, _⟩
  let finalA := hgcdReconstructedLowA entries polyLowA polyLowB sgnS +
    Polynomial.X ^ shift * polyHighA
  let finalB := hgcdReconstructedLowB entries polyLowA polyLowB sgnS +
    Polynomial.X ^ shift * polyHighB
  by_cases hcompute : computeM = true
  · simp [hcompute] at htail
    rcases htail with ⟨combined, hcombine, hheap, hmatrix⟩
    have hcombineWork := hwork.combine combined hcombine
    have hRReconstructed : HgcdMatRawDenseRep this reconstructed.heap R right
        hR := by
      intro i
      exact rawDensePolyRep_of_same_prefix this heap reconstructed.heap
        (hgcdMatPtr R hR i) (hgcdMatLen R hR i) (right i)
        hwork.afterLayout (hwork.rightPrefix i) (hRRep i)
    have hSReconstructed : HgcdMatRawDenseRep this reconstructed.heap S entries
        hS := by
      intro i
      exact rawDensePolyRep_of_same_prefix this heap reconstructed.heap
        (hgcdMatPtr S hS i) (hgcdMatLen S hS i) (entries i)
        hwork.afterLayout (hwork.secondPrefix i) (hSRep i)
    have hQReconstructed := rawDensePolyRep_of_same_prefix this heap
      reconstructed.heap q lenQ quotient hwork.afterLayout
      hwork.quotientPrefix hQ
    have hCombined := hgcdRecursiveCombineMatrix_refines this M R S hM hR hS
      q lenQ T0 a2 scratch reconstructed.heap combined right entries quotient
      hcfg hp hcombineWork.1 hRReconstructed hSReconstructed hQReconstructed
      hcombine
    have hACombined := rawDensePolyRep_of_same_prefix this reconstructed.heap
      combined.heap A reconstructed.lenA finalA hcombineWork.2.1
      hcombineWork.2.2.1 hAReconstructed
    have hBCombined := rawDensePolyRep_of_same_prefix this reconstructed.heap
      combined.heap B reconstructed.lenB finalB hcombineWork.2.1
      hcombineWork.2.2.2 hBReconstructed
    refine ⟨?_, ?_, hsgn, ?_⟩
    · simpa [finalA, hheap, hlenA] using hACombined
    · simpa [finalB, hheap, hlenB] using hBCombined
    · intro _
      simpa [hheap, hmatrix] using hCombined.2.2
  · have hfalse : computeM = false := by cases computeM <;> simp_all
    simp [hfalse] at htail
    refine ⟨?_, ?_, hsgn, ?_⟩
    · simpa [finalA, htail.1, hlenA] using hAReconstructed
    · simpa [finalB, htail.1, hlenB] using hBReconstructed
    · intro htrue
      simp [hfalse] at htrue

/-- The exact final source tail lifts the second recursive transform from its
high suffixes back to the complete divisor/remainder pair.  The transform is
obtained from the same four reconstruction calls that produce the returned
raw operands; the optional matrix block only frames those operands. -/
theorem hgcdRecursiveFinish_preserves_input (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (M R S : HgcdMat) (hM : M.Valid) (hR : R.Valid) (hS : S.Valid)
    (computeM : Bool) (A B T0 lowA lowB highA highB q : RawPtr UInt64)
    (lenLowA lenLowB lenHighA lenHighB shift lenQ inputLength : Nat)
    (a2 scratch : RawPtr UInt64) (sgnR sgnS : Int) (heap : RawHeap)
    (result : HgcdRecursiveFinishResult)
    (right entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (quotient fullA fullB polyLowA polyLowB polyInputHighA polyInputHighB
      polyOutputHighA polyOutputHighB : Polynomial (ZMod this._p.toNat))
    (second : HgcdRecursiveResult)
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (physical : HgcdRecursiveFinishWorkspaceProvider this M R S hM hR hS A
      B T0 lowA lowB highA highB q lenLowA lenLowB lenHighA lenHighB shift
      lenQ a2 scratch sgnS heap)
    (hRRep : HgcdMatRawDenseRep this heap R right hR)
    (hSRep : HgcdMatRawDenseRep this heap S entries hS)
    (hQ : RawDensePolyRep this heap q lenQ quotient)
    (hLowA : RawCanonicalPolySlice this heap lowA lenLowA polyLowA)
    (hLowB : RawCanonicalPolySlice this heap lowB lenLowB polyLowB)
    (hSecond : HgcdRecursiveRawInvariant this polyInputHighA polyInputHighB
      polyOutputHighA polyOutputHighB entries true highA highB inputLength
      second)
    (hSecondHeap : second.heap = heap)
    (hSecondMatrix : second.matrix = S)
    (hSecondValid : HEq second.valid hS)
    (hSecondLenA : second.lenA = lenHighA)
    (hSecondLenB : second.lenB = lenHighB)
    (hSecondSgn : second.sgn = sgnS)
    (hFullA : fullA = polyLowA + Polynomial.X ^ shift * polyInputHighA)
    (hFullB : fullB = polyLowB + Polynomial.X ^ shift * polyInputHighB)
    (hrun : hgcdRecursiveFinish this M R S hM hR hS computeM A B T0 lowA
      lowB highA highB q lenLowA lenLowB lenHighA lenHighB shift lenQ a2
      scratch sgnR sgnS heap = .ok result) :
    RawDensePolyRep this result.heap A result.lenA
        (hgcdReconstructedLowA entries polyLowA polyLowB sgnS +
          Polynomial.X ^ shift * polyOutputHighA) ∧
      RawDensePolyRep this result.heap B result.lenB
        (hgcdReconstructedLowB entries polyLowA polyLowB sgnS +
          Polynomial.X ^ shift * polyOutputHighB) ∧
      CLPoly.Impl.StrictHGCDRefinement.HgcdTransform fullA fullB
        (hgcdReconstructedLowA entries polyLowA polyLowB sgnS +
          Polynomial.X ^ shift * polyOutputHighA)
        (hgcdReconstructedLowB entries polyLowA polyLowB sgnS +
          Polynomial.X ^ shift * polyOutputHighB)
        (entries 0) (entries 1) (entries 2) (entries 3) ∧
      CLPoly.Impl.StrictHGCDRefinement.HgcdSignedDet sgnS
        (entries 0) (entries 1) (entries 2) (entries 3) ∧
      result.sgn = -(sgnR * sgnS) ∧
      (computeM = true →
        HgcdMatRawDenseRep this result.heap result.matrix
          (hgcdMatProductEntry right
            (hgcdMatApplyQuotientEntries entries quotient)) result.valid) := by
  subst heap
  subst S
  subst lenHighA
  subst lenHighB
  subst sgnS
  cases hSecondValid
  have hFinish := hgcdRecursiveFinish_refines this M R second.matrix hM hR
    second.valid computeM A B T0 lowA lowB highA highB q lenLowA lenLowB
    second.lenA second.lenB shift lenQ a2 scratch sgnR second.sgn second.heap
    result right entries quotient polyLowA polyLowB polyOutputHighA
    polyOutputHighB hcfg hp physical hRRep hSRep hQ hLowA hLowB hSecond.aRep
    hSecond.bRep hrun
  rcases hgcdRecursiveFinish_exec this M R second.matrix hM hR second.valid
      computeM A B T0 lowA lowB highA highB q lenLowA lenLowB second.lenA
      second.lenB shift lenQ a2 scratch sgnR second.sgn second.heap result hrun
      with ⟨reconstructed, hreconstruct, _, _, _, _⟩
  have hwork := physical reconstructed hreconstruct
  have hTransform := hgcdRecursiveReconstructPair_preserves_input this A B T0
    lowA lowB highA highB scratch lenLowA lenLowB shift inputLength second
    reconstructed entries fullA fullB polyLowA polyLowB polyInputHighA
    polyInputHighB polyOutputHighA polyOutputHighB hcfg hp hwork.reconstruct
    hSecond hLowA hLowB hFullA hFullB hreconstruct
  exact ⟨hFinish.1, hFinish.2.1, hTransform.2.2.1,
    hTransform.2.2.2.1, hFinish.2.2.1, hFinish.2.2.2⟩

/-- Assemble the semantic portion of the non-early recursive result after
the real middle division, second recursive transform, and finish execution
have supplied their concrete facts.  The returned matrix entries are exactly
the quotient-updated second matrix multiplied by the first matrix. -/
theorem hgcdRecursiveRawInvariant_of_finish_semantics (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (left right currentA currentB remainder quotient finalA finalB :
      Polynomial (ZMod this._p.toNat))
    (first second : Fin 4 → Polynomial (ZMod this._p.toNat))
    (sgnR sgnS : Int) (computeM : Bool)
    (A B : RawPtr UInt64) (inputLength : Nat)
    (result : HgcdRecursiveResult)
    (hARep : RawDensePolyRep this result.heap A result.lenA finalA)
    (hBRep : RawDensePolyRep this result.heap B result.lenB finalB)
    (hMatrix : computeM = true →
      HgcdMatRawDenseRep this result.heap result.matrix
        (hgcdMatProductEntry first
          (hgcdMatApplyQuotientEntries second quotient)) result.valid)
    (hFirstTransform : CLPoly.Impl.StrictHGCDRefinement.HgcdTransform
      left right currentA currentB (first 0) (first 1) (first 2) (first 3))
    (hFirstDet : CLPoly.Impl.StrictHGCDRefinement.HgcdSignedDet sgnR
      (first 0) (first 1) (first 2) (first 3))
    (hDivision : currentA = quotient * currentB + remainder)
    (hSecondTransform : CLPoly.Impl.StrictHGCDRefinement.HgcdTransform
      currentB remainder finalA finalB (second 0) (second 1) (second 2)
      (second 3))
    (hSecondDet : CLPoly.Impl.StrictHGCDRefinement.HgcdSignedDet sgnS
      (second 0) (second 1) (second 2) (second 3))
    (hsgn : result.sgn = -(sgnR * sgnS))
    (hstop : result.lenB < inputLength / 2 + 1)
    (hlengths : computeM = true →
      HgcdRecursiveLengthInvariant inputLength result) :
    HgcdRecursiveRawInvariant this left right finalA finalB
      (hgcdMatProductEntry first
        (hgcdMatApplyQuotientEntries second quotient))
      computeM A B inputLength result := by
  have hTransform := hgcdRecursiveCombined_preserves_transform left right
    currentA currentB remainder quotient finalA finalB first second
    hFirstTransform hDivision hSecondTransform
  have hDet := hgcdRecursiveCombined_preserves_signedDet sgnR sgnS quotient
    first second hFirstDet hSecondDet
  have hResultDet : CLPoly.Impl.StrictHGCDRefinement.HgcdSignedDet result.sgn
      (hgcdMatProductEntry first
        (hgcdMatApplyQuotientEntries second quotient) 0)
      (hgcdMatProductEntry first
        (hgcdMatApplyQuotientEntries second quotient) 1)
      (hgcdMatProductEntry first
        (hgcdMatApplyQuotientEntries second quotient) 2)
      (hgcdMatProductEntry first
        (hgcdMatApplyQuotientEntries second quotient) 3) := by
    rw [hsgn]
    exact hDet
  have hGcd :=
    CLPoly.Impl.StrictHGCDRefinement.normalize_gcd_eq_of_hgcd_signed_transform
      result.sgn left right finalA finalB
      (hgcdMatProductEntry first
        (hgcdMatApplyQuotientEntries second quotient) 0)
      (hgcdMatProductEntry first
        (hgcdMatApplyQuotientEntries second quotient) 1)
      (hgcdMatProductEntry first
        (hgcdMatApplyQuotientEntries second quotient) 2)
      (hgcdMatProductEntry first
        (hgcdMatApplyQuotientEntries second quotient) 3)
      hTransform hResultDet
  exact {
    aRep := hARep
    bRep := hBRep
    matrixSemantics := fun hcompute =>
      ⟨hMatrix hcompute, hTransform, hResultDet⟩
    gcdPreserved := hGcd
    stopped := hstop
    lengths := hlengths }

/-- Package the actual final-tail facts as the common recursive invariant.
The operands and matrix below are the concrete fields of one successful
`hgcdRecursiveFinish` result; the two transforms are those proved for the
first reconstruction and the second reconstruction respectively. -/
theorem hgcdRecursiveRawInvariant_of_finish_execution
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (left right currentA currentB remainder quotient finalA finalB :
      Polynomial (ZMod this._p.toNat))
    (first second : Fin 4 → Polynomial (ZMod this._p.toNat))
    (sgnR sgnS : Int) (computeM : Bool) (A B : RawPtr UInt64)
    (inputLength : Nat) (result : HgcdRecursiveFinishResult)
    (hFinishA : RawDensePolyRep this result.heap A result.lenA finalA)
    (hFinishB : RawDensePolyRep this result.heap B result.lenB finalB)
    (hFinishMatrix : computeM = true →
      HgcdMatRawDenseRep this result.heap result.matrix
        (hgcdMatProductEntry first
          (hgcdMatApplyQuotientEntries second quotient)) result.valid)
    (hFirstTransform : CLPoly.Impl.StrictHGCDRefinement.HgcdTransform
      left right currentA currentB (first 0) (first 1) (first 2) (first 3))
    (hFirstDet : CLPoly.Impl.StrictHGCDRefinement.HgcdSignedDet sgnR
      (first 0) (first 1) (first 2) (first 3))
    (hDivision : currentA = quotient * currentB + remainder)
    (hSecondTransform : CLPoly.Impl.StrictHGCDRefinement.HgcdTransform
      currentB remainder finalA finalB (second 0) (second 1) (second 2)
      (second 3))
    (hSecondDet : CLPoly.Impl.StrictHGCDRefinement.HgcdSignedDet sgnS
      (second 0) (second 1) (second 2) (second 3))
    (hsgn : result.sgn = -(sgnR * sgnS))
    (hstop : result.lenB < inputLength / 2 + 1)
    (hlength : computeM = true →
      HgcdRecursiveLengthInvariant inputLength result.toResult) :
    HgcdRecursiveRawInvariant this left right finalA finalB
      (hgcdMatProductEntry first
        (hgcdMatApplyQuotientEntries second quotient))
      computeM A B inputLength result.toResult := by
  apply hgcdRecursiveRawInvariant_of_finish_semantics this left right currentA
    currentB remainder quotient finalA finalB first second sgnR sgnS computeM
    A B inputLength result.toResult
  · simpa [HgcdRecursiveFinishResult.toResult] using hFinishA
  · simpa [HgcdRecursiveFinishResult.toResult] using hFinishB
  · intro hcompute
    simpa [HgcdRecursiveFinishResult.toResult] using hFinishMatrix hcompute
  · exact hFirstTransform
  · exact hFirstDet
  · exact hDivision
  · exact hSecondTransform
  · exact hSecondDet
  · simpa [HgcdRecursiveFinishResult.toResult] using hsgn
  · exact hstop
  · exact hlength


end CLPoly.Impl.StrictHGCDRawRefinement
