import CLPoly.Generated.StrictHGCD
import CLPoly.Impl.StrictMulDispatchRefinement
import CLPoly.Impl.StrictPolyAddSubRefinement
import CLPoly.Impl.StrictHGCDRefinement

set_option autoImplicit false

namespace CLPoly.Impl.StrictHGCDRawRefinement

open Generated.StrictHGCD
open CLPoly.Impl.RawPolynomialRep
open CLPoly.Impl.StrictDivremRefinement
open CLPoly.Impl.StrictEuclidRefinement
open CLPoly.Impl.StrictPolyAddSubRefinement
open CLPoly.Impl.StrictMulRefinement
open CLPoly.Impl.StrictWordArithmetic

/-- Bounds-safe accessors justified by the C array layout invariant. -/
def hgcdMatPtr (M : HgcdMat) (hM : M.Valid) (i : Fin 4) : RawPtr UInt64 :=
  M.poly[i.val]'(by rw [hM.1]; exact i.isLt)

def hgcdMatLen (M : HgcdMat) (hM : M.Valid) (i : Fin 4) : Nat :=
  M.len[i.val]'(by rw [hM.2]; exact i.isLt)

theorem array_getElem_eq_of_eq {α : Type} (left right : Array α)
    (h : left = right) (i : Nat) (hLeft : i < left.size)
    (hRight : i < right.size) :
    left[i]'hLeft = right[i]'hRight := by
  subst right
  rfl

/-- Raw representation of the four polynomial entries of an HGCD matrix. -/
def HgcdMatPolyRep (heap : RawHeap) (M : HgcdMat) (p : Nat)
    (entries : Fin 4 → Polynomial (ZMod p)) (hM : M.Valid) : Prop :=
  ∀ i : Fin 4, SlicePolyRep heap (hgcdMatPtr M hM i)
    (hgcdMatLen M hM i) p (entries i)

/-- Pointwise heap-frame transport for all four matrix entries. -/
theorem hgcdMatPolyRep_of_same_prefixes (before after : RawHeap)
    (M : HgcdMat) (p : Nat) (entries : Fin 4 → Polynomial (ZMod p))
    (hM : M.Valid) (hlayout : RawHeap.SameLayout before after)
    (hvalid : ∀ i : Fin 4, before.ValidU64Slice
      (hgcdMatPtr M hM i) (hgcdMatLen M hM i))
    (hsame : ∀ i : Fin 4, SameU64Prefix before after
      (hgcdMatPtr M hM i) (hgcdMatLen M hM i))
    (hrep : HgcdMatPolyRep before M p entries hM) :
    HgcdMatPolyRep after M p entries hM := by
  intro i
  exact slicePolyRep_of_same_prefix before after (hgcdMatPtr M hM i)
    (hgcdMatLen M hM i) p (entries i) (hvalid i)
    ((hlayout (hgcdMatPtr M hM i) (hgcdMatLen M hM i)).mp (hvalid i))
    (hsame i) (hrep i)

/-- A successful generated recursive memcpy preserves a represented HGCD
matrix whenever its destination is disjoint from every live matrix slice. -/
theorem copyU64_preserves_hgcdMatPolyRep (heap heap' : RawHeap)
    (dst src : RawPtr UInt64) (count : Nat) (M : HgcdMat) (p : Nat)
    (entries : Fin 4 → Polynomial (ZMod p)) (hM : M.Valid)
    (hDst : heap.ValidU64Slice dst count)
    (hSrc : heap.ValidU64Slice src count)
    (hValidMatrix : ∀ i : Fin 4, heap.ValidU64Slice
      (hgcdMatPtr M hM i) (hgcdMatLen M hM i))
    (hMatrix : ∀ i : Fin 4, U64SlicesDisjoint dst count
      (hgcdMatPtr M hM i) (hgcdMatLen M hM i))
    (hcopy : heap.copyU64 dst src count = .ok heap')
    (hrep : HgcdMatPolyRep heap M p entries hM) :
    RawHeap.SameLayout heap heap' ∧
      HgcdMatPolyRep heap' M p entries hM := by
  rcases copyU64_ok heap dst src count hDst hSrc with
    ⟨copyHeap, hcopy', hlayout⟩
  have heq : copyHeap = heap' := Except.ok.inj (hcopy'.symm.trans hcopy)
  subst copyHeap
  refine ⟨hlayout, hgcdMatPolyRep_of_same_prefixes heap heap' M p entries
    hM hlayout hValidMatrix ?_ hrep⟩
  intro i
  exact copyU64_preserves_prefix heap heap' dst src (hgcdMatPtr M hM i)
    count (hgcdMatLen M hM i) hDst hSrc (hMatrix i) hcopy

noncomputable def identityEntries (p : Nat) : Fin 4 → Polynomial (ZMod p)
  | ⟨0, _⟩ => 1
  | ⟨1, _⟩ => 0
  | ⟨2, _⟩ => 0
  | ⟨3, _⟩ => 1

def identityEntryLen : Fin 4 → Nat
  | ⟨0, _⟩ => 1
  | ⟨1, _⟩ => 0
  | ⟨2, _⟩ => 0
  | ⟨3, _⟩ => 1

/-- A normalized raw polynomial survives any execution that preserves both
its allocation layout and every cell in its declared prefix. -/
theorem rawDensePolyRep_of_same_prefix (this : DenseUPolyZp)
    (before after : RawHeap) (ptr : RawPtr UInt64) (length : Nat)
    (poly : Polynomial (ZMod this._p.toNat))
    (hlayout : RawHeap.SameLayout before after)
    (hsame : SameU64Prefix before after ptr length)
    (hrep : RawDensePolyRep this before ptr length poly) :
    RawDensePolyRep this after ptr length poly := by
  have hvalidAfter := (hlayout ptr length).mp hrep.1
  refine ⟨hvalidAfter,
    canonicalU64Prefix_of_same_prefix before after ptr length this._p
      hrep.1 hsame hrep.2.1,
    slicePolyRep_of_same_prefix before after ptr length this._p.toNat poly
      hrep.1 hvalidAfter hsame hrep.2.2.1, ?_⟩
  have hnorm := normaliseU64_eq_of_prefix_map before after ptr ptr length
    hrep.1 hsame
  exact hnorm.symm.trans hrep.2.2.2

/-- A normalized raw polynomial has zero descriptor length exactly when its
represented polynomial is zero. -/
theorem rawDensePolyRep_length_zero_iff (this : DenseUPolyZp)
    (heap : RawHeap) (ptr : RawPtr UInt64) (length : Nat)
    (poly : Polynomial (ZMod this._p.toNat))
    (hrep : RawDensePolyRep this heap ptr length poly) :
    length = 0 ↔ poly = 0 := by
  constructor
  · intro hzero
    subst length
    exact slicePolyRep_zero_length heap ptr this._p.toNat poly
      (by simpa using hrep.2.2.1)
  · intro hpoly
    by_contra hlength
    have hpositive : 0 < length := Nat.pos_of_ne_zero hlength
    have hlast := normaliseU64_poly_last_coeff_ne_zero heap ptr length
      this._p.toNat length poly hrep.1 hrep.2.2.1 hrep.2.1 hrep.2.2.2
      hlength
    rw [hpoly] at hlast
    simp at hlast

/-- For every nonzero normalized raw descriptor, its C++ length is exactly
one more than the L2 polynomial's natural degree. -/
theorem rawDensePolyRep_natDegree_add_one (this : DenseUPolyZp)
    (heap : RawHeap) (ptr : RawPtr UInt64) (length : Nat)
    (poly : Polynomial (ZMod this._p.toNat))
    (hrep : RawDensePolyRep this heap ptr length poly)
    (hpositive : 0 < length) :
    poly.natDegree + 1 = length := by
  have hdegree := normaliseU64_poly_natDegree_eq heap ptr length
    this._p.toNat length poly hrep.1 hrep.2.2.1 hrep.2.1 hrep.2.2.2
    (Nat.ne_of_gt hpositive)
  omega

/-- Convert an L2 degree bound back to the normalized C++ descriptor length.
The zero case uses the exact zero-length equivalence above. -/
theorem rawDensePolyRep_length_le_of_degree_lt (this : DenseUPolyZp)
    (heap : RawHeap) (ptr : RawPtr UInt64) (length bound : Nat)
    (poly : Polynomial (ZMod this._p.toNat))
    (hrep : RawDensePolyRep this heap ptr length poly)
    (hdegree : poly = 0 ∨ poly.natDegree < bound) :
    length ≤ bound := by
  rcases hdegree with hzero | hdegree
  · have := (rawDensePolyRep_length_zero_iff this heap ptr length poly
      hrep).mpr hzero
    omega
  · by_cases hlength : length = 0
    · omega
    · have hexact := rawDensePolyRep_natDegree_add_one this heap ptr length
        poly hrep (Nat.pos_of_ne_zero hlength)
      omega

/-- Degree bound for one product expressed solely in terms of the two exact
normalized C++ descriptor lengths. -/
theorem rawDensePolyRep_mul_zero_or_degree_lt (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (leftHeap rightHeap : RawHeap) (leftPtr rightPtr : RawPtr UInt64)
    (leftLength rightLength : Nat)
    (left right : Polynomial (ZMod this._p.toNat))
    (hLeft : RawDensePolyRep this leftHeap leftPtr leftLength left)
    (hRight : RawDensePolyRep this rightHeap rightPtr rightLength right) :
    left * right = 0 ∨
      (left * right).natDegree < leftLength + rightLength - 1 := by
  by_cases hLeftLength : leftLength = 0
  · left
    rw [(rawDensePolyRep_length_zero_iff this leftHeap leftPtr leftLength left
      hLeft).mp hLeftLength, zero_mul]
  · by_cases hRightLength : rightLength = 0
    · left
      rw [(rawDensePolyRep_length_zero_iff this rightHeap rightPtr rightLength right
        hRight).mp hRightLength, mul_zero]
    · right
      have hLeftDegree := rawDensePolyRep_natDegree_add_one this leftHeap leftPtr
        leftLength left hLeft (Nat.pos_of_ne_zero hLeftLength)
      have hRightDegree := rawDensePolyRep_natDegree_add_one this rightHeap rightPtr
        rightLength right hRight (Nat.pos_of_ne_zero hRightLength)
      have hLeftNonzero : left ≠ 0 := by
        intro hzero
        exact hLeftLength ((rawDensePolyRep_length_zero_iff this leftHeap leftPtr
          leftLength left hLeft).mpr hzero)
      have hRightNonzero : right ≠ 0 := by
        intro hzero
        exact hRightLength ((rawDensePolyRep_length_zero_iff this rightHeap rightPtr
          rightLength right hRight).mpr hzero)
      rw [Polynomial.natDegree_mul hLeftNonzero hRightNonzero]
      omega

/-- A normalized raw descriptor is either the zero polynomial or its stored
length strictly dominates its degree. -/
theorem rawDensePolyRep_zero_or_degree_lt (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (heap : RawHeap) (ptr : RawPtr UInt64) (length : Nat)
    (poly : Polynomial (ZMod this._p.toNat))
    (hrep : RawDensePolyRep this heap ptr length poly) :
    poly = 0 ∨ poly.natDegree < length := by
  by_cases hlength : length = 0
  · left
    exact (rawDensePolyRep_length_zero_iff this heap ptr length poly hrep).mp
      hlength
  · right
    have hdegree := rawDensePolyRep_natDegree_add_one this heap ptr length
      poly hrep (Nat.pos_of_ne_zero hlength)
    omega

/-- Normalization makes the descriptor length uniquely determined by the L2
polynomial, even when the two real buffers live in different heaps. -/
theorem rawDensePolyRep_length_eq (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (heap₁ heap₂ : RawHeap) (ptr₁ ptr₂ : RawPtr UInt64)
    (length₁ length₂ : Nat) (poly : Polynomial (ZMod this._p.toNat))
    (hrep₁ : RawDensePolyRep this heap₁ ptr₁ length₁ poly)
    (hrep₂ : RawDensePolyRep this heap₂ ptr₂ length₂ poly) :
    length₁ = length₂ := by
  by_cases hzero : poly = 0
  · have hzero₁ := (rawDensePolyRep_length_zero_iff this heap₁ ptr₁ length₁
      poly hrep₁).mpr hzero
    have hzero₂ := (rawDensePolyRep_length_zero_iff this heap₂ ptr₂ length₂
      poly hrep₂).mpr hzero
    omega
  · have hpos₁ : 0 < length₁ := by
      by_contra hnot
      have hlength : length₁ = 0 := by omega
      exact hzero ((rawDensePolyRep_length_zero_iff this heap₁ ptr₁ length₁
        poly hrep₁).mp hlength)
    have hpos₂ : 0 < length₂ := by
      by_contra hnot
      have hlength : length₂ = 0 := by omega
      exact hzero ((rawDensePolyRep_length_zero_iff this heap₂ ptr₂ length₂
        poly hrep₂).mp hlength)
    have hdegree₁ := rawDensePolyRep_natDegree_add_one this heap₁ ptr₁ length₁
      poly hrep₁ hpos₁
    have hdegree₂ := rawDensePolyRep_natDegree_add_one this heap₂ ptr₂ length₂
      poly hrep₂ hpos₂
    omega

/-- When the normalized low part lies strictly below the shift and the high
descriptor is nonempty, the shifted high term determines the exact normalized
output length.  This is the semantic form of the real `liftHigh` leading-limb
preservation used by recursive HGCD. -/
theorem rawDensePolyRep_add_shift_length_eq_of_lt (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (lowHeap highHeap outHeap : RawHeap)
    (lowPtr highPtr outPtr : RawPtr UInt64)
    (lowLength highLength outLength shift : Nat)
    (low high : Polynomial (ZMod this._p.toNat))
    (hLow : RawDensePolyRep this lowHeap lowPtr lowLength low)
    (hHigh : RawDensePolyRep this highHeap highPtr highLength high)
    (hOut : RawDensePolyRep this outHeap outPtr outLength
      (low + Polynomial.X ^ shift * high))
    (hLowBound : lowLength < shift + highLength)
    (hHighPos : 0 < highLength) :
    outLength = shift + highLength := by
  have hHighNonzero : high ≠ 0 := by
    intro hzero
    have := (rawDensePolyRep_length_zero_iff this highHeap highPtr highLength
      high hHigh).mpr hzero
    omega
  have hXNonzero : (Polynomial.X ^ shift :
      Polynomial (ZMod this._p.toNat)) ≠ 0 := pow_ne_zero _ Polynomial.X_ne_zero
  have hShiftNonzero : Polynomial.X ^ shift * high ≠ 0 :=
    mul_ne_zero hXNonzero hHighNonzero
  have hHighDegree := rawDensePolyRep_natDegree_add_one this highHeap highPtr
    highLength high hHigh hHighPos
  have hShiftDegree : (Polynomial.X ^ shift * high).natDegree =
      shift + high.natDegree := by
    rw [Polynomial.natDegree_mul hXNonzero hHighNonzero,
      Polynomial.natDegree_X_pow]
  have hOutNonzero : low + Polynomial.X ^ shift * high ≠ 0 := by
    intro hzero
    have heq : low = -(Polynomial.X ^ shift * high) := eq_neg_of_add_eq_zero_left
      hzero
    have hLowDegreeEq : low.natDegree =
        (Polynomial.X ^ shift * high).natDegree := by
      rw [heq, Polynomial.natDegree_neg]
    rcases rawDensePolyRep_zero_or_degree_lt this lowHeap lowPtr lowLength low
        hLow with hLowZero | hLowDegree
    · exact hShiftNonzero (by simpa [hLowZero] using hzero.symm)
    · rw [hShiftDegree] at hLowDegreeEq
      omega
  have hOutPos : 0 < outLength := by
    by_contra hnot
    have hlength : outLength = 0 := by omega
    exact hOutNonzero ((rawDensePolyRep_length_zero_iff this outHeap outPtr
      outLength (low + Polynomial.X ^ shift * high) hOut).mp hlength)
  have hOutDegree : (low + Polynomial.X ^ shift * high).natDegree =
      shift + high.natDegree := by
    rcases rawDensePolyRep_zero_or_degree_lt this lowHeap lowPtr lowLength low
        hLow with hLowZero | hLowDegree
    · simp [hLowZero, hShiftDegree]
    · rw [Polynomial.natDegree_add_eq_right_of_natDegree_lt]
      · exact hShiftDegree
      · rw [hShiftDegree]
        omega
  have hNormalized := rawDensePolyRep_natDegree_add_one this outHeap outPtr
    outLength (low + Polynomial.X ^ shift * high) hOut hOutPos
  rw [hOutDegree] at hNormalized
  omega

/-- The common disjoint-low/high specialization of
`rawDensePolyRep_add_shift_length_eq_of_lt`. -/
theorem rawDensePolyRep_add_shift_length_eq (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (lowHeap highHeap outHeap : RawHeap)
    (lowPtr highPtr outPtr : RawPtr UInt64)
    (lowLength highLength outLength shift : Nat)
    (low high : Polynomial (ZMod this._p.toNat))
    (hLow : RawDensePolyRep this lowHeap lowPtr lowLength low)
    (hHigh : RawDensePolyRep this highHeap highPtr highLength high)
    (hOut : RawDensePolyRep this outHeap outPtr outLength
      (low + Polynomial.X ^ shift * high))
    (hLowBound : lowLength ≤ shift) (hHighPos : 0 < highLength) :
    outLength = shift + highLength := by
  apply rawDensePolyRep_add_shift_length_eq_of_lt this lowHeap highHeap outHeap
    lowPtr highPtr outPtr lowLength highLength outLength shift low high hLow
    hHigh hOut
  · omega
  · exact hHighPos

/-- Exact normalized length bound for the quotient update used by the real
HGCD tail, `top + quotient * bottom`. -/
theorem rawDensePolyRep_add_mul_length_le (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (topHeap quotientHeap bottomHeap outHeap : RawHeap)
    (topPtr quotientPtr bottomPtr outPtr : RawPtr UInt64)
    (topLength quotientLength bottomLength outLength : Nat)
    (top quotient bottom : Polynomial (ZMod this._p.toNat))
    (hTop : RawDensePolyRep this topHeap topPtr topLength top)
    (hQuotient : RawDensePolyRep this quotientHeap quotientPtr quotientLength quotient)
    (hBottom : RawDensePolyRep this bottomHeap bottomPtr bottomLength bottom)
    (hOut : RawDensePolyRep this outHeap outPtr outLength
      (top + quotient * bottom)) :
    outLength ≤ max topLength (quotientLength + bottomLength - 1) := by
  have hTopDegree := rawDensePolyRep_zero_or_degree_lt this topHeap topPtr
    topLength top hTop
  have hProductDegree := rawDensePolyRep_mul_zero_or_degree_lt this quotientHeap
    bottomHeap quotientPtr bottomPtr quotientLength bottomLength quotient bottom
    hQuotient hBottom
  apply rawDensePolyRep_length_le_of_degree_lt this outHeap outPtr outLength
    (max topLength (quotientLength + bottomLength - 1))
    (top + quotient * bottom) hOut
  rcases hTopDegree with hTopZero | hTopDegree <;>
      rcases hProductDegree with hProductZero | hProductDegree
  · left
    rw [hTopZero, hProductZero, zero_add]
  · right
    rw [hTopZero, zero_add]
    exact hProductDegree.trans_le (Nat.le_max_right _ _)
  · right
    rw [hProductZero, add_zero]
    exact hTopDegree.trans_le (Nat.le_max_left _ _)
  · right
    exact (Polynomial.natDegree_add_le top (quotient * bottom)).trans_lt (by
      exact max_lt (hTopDegree.trans_le (Nat.le_max_left _ _))
        (hProductDegree.trans_le (Nat.le_max_right _ _)))

/-- Exact normalized length bound for the polynomial produced by one real
matrix-entry computation `P*Q + R*S`. -/
theorem rawDensePolyRep_sum_products_length_le (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (heap : RawHeap) (pPtr qPtr rPtr sPtr outPtr : RawPtr UInt64)
    (lenP lenQ lenR lenS outLength : Nat)
    (P Q R S : Polynomial (ZMod this._p.toNat))
    (hP : RawDensePolyRep this heap pPtr lenP P)
    (hQ : RawDensePolyRep this heap qPtr lenQ Q)
    (hR : RawDensePolyRep this heap rPtr lenR R)
    (hS : RawDensePolyRep this heap sPtr lenS S)
    (hOut : RawDensePolyRep this heap outPtr outLength (P * Q + R * S)) :
    outLength ≤ max (lenP + lenQ - 1) (lenR + lenS - 1) := by
  have hPQ := rawDensePolyRep_mul_zero_or_degree_lt this heap heap pPtr qPtr lenP
    lenQ P Q hP hQ
  have hRS := rawDensePolyRep_mul_zero_or_degree_lt this heap heap rPtr sPtr lenR
    lenS R S hR hS
  apply rawDensePolyRep_length_le_of_degree_lt this heap outPtr outLength
    (max (lenP + lenQ - 1) (lenR + lenS - 1)) (P * Q + R * S) hOut
  rcases hPQ with hPQ | hPQ <;> rcases hRS with hRS | hRS
  · left
    rw [hPQ, hRS, zero_add]
  · right
    rw [hPQ, zero_add]
    exact hRS.trans_le (Nat.le_max_right _ _)
  · right
    rw [hRS, add_zero]
    exact hPQ.trans_le (Nat.le_max_left _ _)
  · right
    exact (Polynomial.natDegree_add_le (P * Q) (R * S)).trans_lt (by
      exact max_lt (hPQ.trans_le (Nat.le_max_left _ _))
        (hRS.trans_le (Nat.le_max_right _ _)))

/-- Normalized raw representation of all four live HGCD matrix entries. -/
def HgcdMatRawDenseRep (this : DenseUPolyZp) (heap : RawHeap)
    (M : HgcdMat) (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (hM : M.Valid) : Prop :=
  ∀ i : Fin 4, RawDensePolyRep this heap (hgcdMatPtr M hM i)
    (hgcdMatLen M hM i) (entries i)

/-- Purely physical obligations of the two C++ matrix-stabilization loops.
No L2 polynomial value or expected algorithm result is stored here. -/
structure HgcdMatStabilizeWorkspace (heap : RawHeap)
    (original current : HgcdMat) (hOriginal : original.Valid)
    (hCurrent : current.Valid) (stage : RawPtr UInt64) : Prop where
  stageValid : heap.ValidU64Slice stage
    (hgcdMatStageSize current hCurrent)
  stageCurrentDisjoint : ∀ i : Fin 4,
    U64SlicesDisjoint stage (hgcdMatStageSize current hCurrent)
      (hgcdMatPtr current hCurrent i) (hgcdMatLen current hCurrent i)
  originalValid : ∀ i : Fin 4, heap.ValidU64Slice
    (hgcdMatPtr original hOriginal i) (hgcdMatLen current hCurrent i)
  originalStageDisjoint : ∀ i : Fin 4,
    U64SlicesDisjoint (hgcdMatPtr original hOriginal i)
      (hgcdMatLen current hCurrent i) stage
      (hgcdMatStageSize current hCurrent)
  originalPairwiseDisjoint : ∀ i j : Fin 4, i ≠ j →
    U64SlicesDisjoint (hgcdMatPtr original hOriginal i)
      (hgcdMatLen current hCurrent i) (hgcdMatPtr original hOriginal j)
      (hgcdMatLen current hCurrent j)

theorem validU64Slice_stage_entry (heap : RawHeap) (stage : RawPtr UInt64)
    (M : HgcdMat) (hM : M.Valid) (i : Fin 4)
    (hvalid : heap.ValidU64Slice stage (hgcdMatStageSize M hM)) :
    heap.ValidU64Slice (stage.add (hgcdMatStageOffset M hM i.val))
      (hgcdMatLen M hM i) := by
  exact heap.validU64Slice_add stage (hgcdMatStageSize M hM)
    (hgcdMatStageOffset M hM i.val) (hgcdMatLen M hM i) hvalid
    (hgcdMatStageOffset_entry_le_size M hM i.val i.isLt)

theorem u64SlicesDisjoint_stage_entry_left (stage other : RawPtr UInt64)
    (total offset length otherLength : Nat)
    (hbound : offset + length ≤ total)
    (hdisjoint : U64SlicesDisjoint stage total other otherLength) :
    U64SlicesDisjoint (stage.add offset) length other otherLength := by
  intro i hi j hj
  rcases hdisjoint (offset + i) (by omega) j hj with hregion | hoff
  · exact Or.inl hregion
  · right
    have hwidth : RawLimbWidth.width UInt64 = 1 := rfl
    simpa [RawPtr.add, hwidth, Nat.add_assoc, Nat.add_left_comm, Nat.add_comm]
      using hoff

theorem u64SlicesDisjoint_stage_entry_right (other stage : RawPtr UInt64)
    (otherLength total offset length : Nat)
    (hbound : offset + length ≤ total)
    (hdisjoint : U64SlicesDisjoint other otherLength stage total) :
    U64SlicesDisjoint other otherLength (stage.add offset) length := by
  exact u64SlicesDisjoint_symm
    (u64SlicesDisjoint_stage_entry_left stage other total offset length
      otherLength hbound (u64SlicesDisjoint_symm hdisjoint))

theorem hgcdMatStabilizeWorkspace_of_sameLayout (before after : RawHeap)
    (original current : HgcdMat) (hOriginal : original.Valid)
    (hCurrent : current.Valid) (stage : RawPtr UInt64)
    (hlayout : RawHeap.SameLayout before after)
    (hws : HgcdMatStabilizeWorkspace before original current hOriginal
      hCurrent stage) :
    HgcdMatStabilizeWorkspace after original current hOriginal hCurrent stage := by
  exact {
    stageValid := (hlayout stage (hgcdMatStageSize current hCurrent)).mp
      hws.stageValid
    stageCurrentDisjoint := hws.stageCurrentDisjoint
    originalValid := fun i =>
      (hlayout (hgcdMatPtr original hOriginal i)
        (hgcdMatLen current hCurrent i)).mp (hws.originalValid i)
    originalStageDisjoint := hws.originalStageDisjoint
    originalPairwiseDisjoint := hws.originalPairwiseDisjoint }

/-- Raw representation of the four source-order entries after the first
stabilization loop has copied them into its contiguous stage buffer. -/
def HgcdMatStageRawDenseRep (this : DenseUPolyZp) (heap : RawHeap)
    (stage : RawPtr UInt64) (M : HgcdMat)
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (hM : M.Valid) : Prop :=
  ∀ i : Fin 4, RawDensePolyRep this heap
    (stage.add (hgcdMatStageOffset M hM i.val))
    (hgcdMatLen M hM i) (entries i)

/-- Pointwise frame transport for the normalized raw matrix invariant. -/
theorem hgcdMatRawDenseRep_of_same_prefixes (this : DenseUPolyZp)
    (before after : RawHeap) (M : HgcdMat)
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (hM : M.Valid) (hlayout : RawHeap.SameLayout before after)
    (hsame : ∀ i : Fin 4, SameU64Prefix before after
      (hgcdMatPtr M hM i) (hgcdMatLen M hM i))
    (hrep : HgcdMatRawDenseRep this before M entries hM) :
    HgcdMatRawDenseRep this after M entries hM := by
  intro i
  exact rawDensePolyRep_of_same_prefix this before after
    (hgcdMatPtr M hM i) (hgcdMatLen M hM i) (entries i) hlayout
    (hsame i) (hrep i)

/-- A real recursive memcpy preserves the normalized raw matrix invariant
when its destination is disjoint from every live entry. -/
theorem copyU64_preserves_hgcdMatRawDenseRep (this : DenseUPolyZp)
    (heap heap' : RawHeap) (dst src : RawPtr UInt64) (count : Nat)
    (M : HgcdMat) (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (hM : M.Valid) (hDst : heap.ValidU64Slice dst count)
    (hSrc : heap.ValidU64Slice src count)
    (hMatrix : ∀ i : Fin 4, U64SlicesDisjoint dst count
      (hgcdMatPtr M hM i) (hgcdMatLen M hM i))
    (hcopy : heap.copyU64 dst src count = .ok heap')
    (hrep : HgcdMatRawDenseRep this heap M entries hM) :
    RawHeap.SameLayout heap heap' ∧
      HgcdMatRawDenseRep this heap' M entries hM := by
  rcases copyU64_ok heap dst src count hDst hSrc with
    ⟨copyHeap, hcopy', hlayout⟩
  have heq : copyHeap = heap' := Except.ok.inj (hcopy'.symm.trans hcopy)
  subst copyHeap
  refine ⟨hlayout, hgcdMatRawDenseRep_of_same_prefixes this heap heap' M
    entries hM hlayout ?_ hrep⟩
  intro i
  exact copyU64_preserves_prefix heap heap' dst src (hgcdMatPtr M hM i)
    count (hgcdMatLen M hM i) hDst hSrc (hMatrix i) hcopy

/-- The actual recursive RawHeap memcpy transports the complete normalized
polynomial representation to a disjoint destination. -/
theorem copyU64_refines_rawDense (this : DenseUPolyZp)
    (heap : RawHeap) (dst src : RawPtr UInt64) (length : Nat)
    (poly : Polynomial (ZMod this._p.toNat))
    (hDst : heap.ValidU64Slice dst length)
    (hdisjoint : U64SlicesDisjoint dst length src length)
    (hrep : RawDensePolyRep this heap src length poly) :
    ∃ heap', heap.copyU64 dst src length = .ok heap' ∧
      RawHeap.SameLayout heap heap' ∧
      RawDensePolyRep this heap' dst length poly := by
  rcases copyU64_refines_disjoint heap dst src length hDst hrep.1 hdisjoint
      with ⟨heap', hcopy, hlayout, hcontents⟩
  rcases copyU64_refines_slice_canonical heap dst src length
      this._p.toNat poly this._p hDst hrep.1 hdisjoint hrep.2.2.1
      hrep.2.1 with ⟨repHeap, hcopyRep, _, hslice, hcanonical⟩
  have heq : repHeap = heap' := Except.ok.inj (hcopyRep.symm.trans hcopy)
  subst repHeap
  have hnorm := normaliseU64_eq_of_prefix_map heap heap' src dst length
    hrep.1 hcontents
  exact ⟨heap', hcopy, hlayout, (hlayout dst length).mp hDst, hcanonical,
    hslice, hnorm.symm.trans hrep.2.2.2⟩

/-- A represented raw polynomial disjoint from the memcpy destination is a
frame of that actual generated copy. -/
theorem copyU64_preserves_rawDenseRep (this : DenseUPolyZp)
    (heap heap' : RawHeap) (dst src : RawPtr UInt64) (count : Nat)
    (ptr : RawPtr UInt64) (length : Nat)
    (poly : Polynomial (ZMod this._p.toNat))
    (hDst : heap.ValidU64Slice dst count)
    (hSrc : heap.ValidU64Slice src count)
    (hdisjoint : U64SlicesDisjoint dst count ptr length)
    (hcopy : heap.copyU64 dst src count = .ok heap')
    (hrep : RawDensePolyRep this heap ptr length poly) :
    RawHeap.SameLayout heap heap' ∧
      RawDensePolyRep this heap' ptr length poly := by
  rcases copyU64_ok heap dst src count hDst hSrc with
    ⟨copyHeap, hcopy', hlayout⟩
  have heq : copyHeap = heap' := Except.ok.inj (hcopy'.symm.trans hcopy)
  subst copyHeap
  exact ⟨hlayout, rawDensePolyRep_of_same_prefix this heap heap' ptr length
    poly hlayout (copyU64_preserves_prefix heap heap' dst src ptr count length
      hDst hSrc hdisjoint hcopy) hrep⟩

/-- Semantic refinement of the first actual stabilization loop.  The induction
follows its generated `i` recursion and accumulates the already copied stage
entries; the workspace contributes only validity and non-aliasing. -/
theorem hgcdMatStageLoop_refines (this : DenseUPolyZp)
    (M original : HgcdMat) (hM : M.Valid) (hOriginal : original.Valid)
    (stage : RawPtr UInt64)
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (i : Nat) (hi : i ≤ 4) (heap : RawHeap)
    (hws : HgcdMatStabilizeWorkspace heap original M hOriginal hM stage)
    (hMatrix : HgcdMatRawDenseRep this heap M entries hM)
    (hPrior : ∀ j : Fin 4, j.val < i → RawDensePolyRep this heap
      (stage.add (hgcdMatStageOffset M hM j.val))
      (hgcdMatLen M hM j) (entries j)) :
    ∃ result, hgcdMatStageLoop M hM stage i
        (hgcdMatStageOffset M hM i) heap = .ok result ∧
      result.off = hgcdMatStageSize M hM ∧
      RawHeap.SameLayout heap result.heap ∧
      HgcdMatRawDenseRep this result.heap M entries hM ∧
      HgcdMatStageRawDenseRep this result.heap stage M entries hM := by
  rw [hgcdMatStageLoop]
  split
  next hlt =>
    dsimp only
    let index : Fin 4 := ⟨i, hlt⟩
    let dst := stage.add (hgcdMatStageOffset M hM i)
    let length := hgcdMatLen M hM index
    have hDst : heap.ValidU64Slice dst length :=
      validU64Slice_stage_entry heap stage M hM index hws.stageValid
    have hSrc := hMatrix index
    have hDstSrc : U64SlicesDisjoint dst length
        (hgcdMatPtr M hM index) length := by
      exact u64SlicesDisjoint_stage_entry_left stage (hgcdMatPtr M hM index)
        (hgcdMatStageSize M hM) (hgcdMatStageOffset M hM i) length length
        (hgcdMatStageOffset_entry_le_size M hM i hlt)
        (hws.stageCurrentDisjoint index)
    rcases copyU64_refines_rawDense this heap dst (hgcdMatPtr M hM index)
        length (entries index) hDst hDstSrc hSrc with
      ⟨heap1, hcopy, hlayout, hStageI⟩
    have hMatrix1 : HgcdMatRawDenseRep this heap1 M entries hM := by
      apply (copyU64_preserves_hgcdMatRawDenseRep this heap heap1 dst
        (hgcdMatPtr M hM index) length M entries hM hDst hSrc.1 ?_
        hcopy hMatrix).2
      intro k
      exact u64SlicesDisjoint_stage_entry_left stage (hgcdMatPtr M hM k)
        (hgcdMatStageSize M hM) (hgcdMatStageOffset M hM i) length
        (hgcdMatLen M hM k)
        (hgcdMatStageOffset_entry_le_size M hM i hlt)
        (hws.stageCurrentDisjoint k)
    have hPrior1 : ∀ j : Fin 4, j.val < i + 1 →
        RawDensePolyRep this heap1
          (stage.add (hgcdMatStageOffset M hM j.val))
          (hgcdMatLen M hM j) (entries j) := by
      intro j hj
      by_cases hji : j.val = i
      · have hjeq : j = index := Fin.ext hji
        subst j
        simpa [dst, length, hgcdMatPtr, hgcdMatPtrRaw, hgcdMatLen,
          hgcdMatLenRaw] using hStageI
      · have hjlt : j.val < i := by omega
        have hOld := hPrior j hjlt
        have hBefore :
            hgcdMatStageOffset M hM j.val + hgcdMatLen M hM j ≤
              hgcdMatStageOffset M hM i := by
          exact hgcdMatStageOffset_entry_le_later M hM j.val i hjlt hlt j.isLt
        have hDisjoint : U64SlicesDisjoint dst length
            (stage.add (hgcdMatStageOffset M hM j.val))
            (hgcdMatLen M hM j) := by
          exact u64SlicesDisjoint_symm
            (u64SlicesDisjoint_add_of_le stage
              (hgcdMatStageOffset M hM j.val) (hgcdMatLen M hM j)
              (hgcdMatStageOffset M hM i) length hBefore)
        exact (copyU64_preserves_rawDenseRep this heap heap1 dst
          (hgcdMatPtr M hM index) length
          (stage.add (hgcdMatStageOffset M hM j.val))
          (hgcdMatLen M hM j) (entries j) hDst hSrc.1 hDisjoint hcopy hOld).2
    have hws1 := hgcdMatStabilizeWorkspace_of_sameLayout heap heap1 original M
      hOriginal hM stage hlayout hws
    have hcopyLoop : heap.copyU64
        (stage.add (hgcdMatStageOffset M hM i))
        (hgcdMatPtrRaw M hM index) (hgcdMatLenRaw M hM index) = .ok heap1 := by
      simpa [dst, length, hgcdMatPtr, hgcdMatPtrRaw, hgcdMatLen,
        hgcdMatLenRaw] using hcopy
    have hrec := hgcdMatStageLoop_refines this M original hM hOriginal stage
      entries (i + 1) (by omega) heap1 hws1 hMatrix1 hPrior1
    rw [hcopyLoop]
    rcases hrec with ⟨final, hrunFinal, hoffFinal, hlayoutFinal,
      hMatrixFinal, hStageFinal⟩
    refine ⟨final, ?_, hoffFinal, ?_, hMatrixFinal, hStageFinal⟩
    · simpa [index, dst, length, hgcdMatPtr, hgcdMatPtrRaw, hgcdMatLen,
        hgcdMatLenRaw, hgcdMatStageOffset_step M hM i hlt] using hrunFinal
    · intro ptr count
      exact (hlayout ptr count).trans (hlayoutFinal ptr count)
  next hstop =>
    have hi4 : i = 4 := by omega
    subst i
    refine ⟨{ heap := heap, off := hgcdMatStageOffset M hM 4 }, rfl,
      rfl, (fun _ _ => Iff.rfl), hMatrix, ?_⟩
    intro j
    exact hPrior j (by omega)
termination_by 4 - i
decreasing_by omega

/-- Every concrete copy in the staging loop is a frame for a represented
live polynomial allocated outside the staging region.  This follows the
generated loop itself; no semantic result is supplied by the caller. -/
theorem hgcdMatStageLoop_preserves_rawDenseRep (this : DenseUPolyZp)
    (M original : HgcdMat) (hM : M.Valid) (hOriginal : original.Valid)
    (stage : RawPtr UInt64)
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (i : Nat) (hi : i ≤ 4) (heap : RawHeap)
    (hws : HgcdMatStabilizeWorkspace heap original M hOriginal hM stage)
    (hMatrix : HgcdMatRawDenseRep this heap M entries hM)
    (ptr : RawPtr UInt64) (length : Nat)
    (poly : Polynomial (ZMod this._p.toNat))
    (hStagePtr : stage.region ≠ ptr.region)
    (hRep : RawDensePolyRep this heap ptr length poly)
    (result : HgcdMatStageResult)
    (hrun : hgcdMatStageLoop M hM stage i
      (hgcdMatStageOffset M hM i) heap = .ok result) :
    RawHeap.SameLayout heap result.heap ∧
      RawDensePolyRep this result.heap ptr length poly := by
  rw [hgcdMatStageLoop] at hrun
  split at hrun
  next hlt =>
    dsimp only at hrun
    split at hrun
    next fault hcopy => simp at hrun
    next heap1 hcopy =>
      let index : Fin 4 := ⟨i, hlt⟩
      let dst := stage.add (hgcdMatStageOffset M hM i)
      let count := hgcdMatLen M hM index
      have hDst : heap.ValidU64Slice dst count :=
        validU64Slice_stage_entry heap stage M hM index hws.stageValid
      have hSrc := hMatrix index
      have hLiveDisjoint : U64SlicesDisjoint dst count ptr length := by
        apply u64SlicesDisjoint_of_region_ne
        simpa [dst, RawPtr.add] using hStagePtr
      have hFrame := copyU64_preserves_rawDenseRep this heap heap1 dst
        (hgcdMatPtr M hM index) count ptr length poly hDst hSrc.1
        hLiveDisjoint (by
          simpa [dst, count, index, hgcdMatPtr, hgcdMatPtrRaw, hgcdMatLen,
            hgcdMatLenRaw] using hcopy) hRep
      have hMatrix1 : HgcdMatRawDenseRep this heap1 M entries hM := by
        apply (copyU64_preserves_hgcdMatRawDenseRep this heap heap1 dst
          (hgcdMatPtr M hM index) count M entries hM hDst hSrc.1 ?_ (by
            simpa [dst, count, index, hgcdMatPtr, hgcdMatPtrRaw, hgcdMatLen,
              hgcdMatLenRaw] using hcopy) hMatrix).2
        intro k
        exact u64SlicesDisjoint_stage_entry_left stage (hgcdMatPtr M hM k)
          (hgcdMatStageSize M hM) (hgcdMatStageOffset M hM i) count
          (hgcdMatLen M hM k)
          (hgcdMatStageOffset_entry_le_size M hM i hlt)
          (hws.stageCurrentDisjoint k)
      have hws1 := hgcdMatStabilizeWorkspace_of_sameLayout heap heap1
        original M hOriginal hM stage hFrame.1 hws
      have hrec := hgcdMatStageLoop_preserves_rawDenseRep this M original hM
        hOriginal stage entries (i + 1) (by omega) heap1 hws1 hMatrix1 ptr
        length poly hStagePtr hFrame.2 result
        (by simpa [index, dst, count, hgcdMatPtr, hgcdMatPtrRaw, hgcdMatLen,
          hgcdMatLenRaw, hgcdMatStageOffset_step M hM i hlt] using hrun)
      exact ⟨fun p n => (hFrame.1 p n).trans (hrec.1 p n), hrec.2⟩
  next hstop =>
    have heq := Except.ok.inj hrun
    subst result
    exact ⟨fun _ _ => Iff.rfl, hRep⟩
termination_by 4 - i
decreasing_by omega

/-- Entry-point refinement of the first generated stabilization loop. -/
theorem hgcdMatStageLoop_zero_refines (this : DenseUPolyZp)
    (M original : HgcdMat) (hM : M.Valid) (hOriginal : original.Valid)
    (stage : RawPtr UInt64)
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (heap : RawHeap)
    (hws : HgcdMatStabilizeWorkspace heap original M hOriginal hM stage)
    (hMatrix : HgcdMatRawDenseRep this heap M entries hM) :
    ∃ result, hgcdMatStageLoop M hM stage 0 0 heap = .ok result ∧
      result.off = hgcdMatStageSize M hM ∧
      RawHeap.SameLayout heap result.heap ∧
      HgcdMatRawDenseRep this result.heap M entries hM ∧
      HgcdMatStageRawDenseRep this result.heap stage M entries hM := by
  simpa [hgcdMatStageOffset] using
    (hgcdMatStageLoop_refines this M original hM hOriginal stage entries 0
      (by omega) heap hws hMatrix (by intro j hj; omega))

/-- Semantic refinement of every suffix of the second generated stabilization
loop.  `base` fixes the iterator-produced lengths while `current` follows the
real descriptor pointer updates performed by the loop. -/
theorem hgcdMatRestoreLoop_refines (this : DenseUPolyZp)
    (original base current : HgcdMat)
    (hOriginal : original.Valid) (hBase : base.Valid)
    (hCurrent : current.Valid)
    (stage : RawPtr UInt64)
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (i : Nat) (hi : i ≤ 4) (heap : RawHeap)
    (result : HgcdMatRestoreResult)
    (hLen : current.len = base.len)
    (hws : HgcdMatStabilizeWorkspace heap original base hOriginal hBase stage)
    (hStage : HgcdMatStageRawDenseRep this heap stage base entries hBase)
    (hPrior : ∀ j : Fin 4, j.val < i → RawDensePolyRep this heap
      (hgcdMatPtr original hOriginal j) (hgcdMatLen base hBase j)
      (entries j))
    (hrun : hgcdMatRestoreLoop original current hOriginal hCurrent stage i
      (hgcdMatStageOffset base hBase i) heap = .ok result) :
    HgcdMatStageRawDenseRep this result.heap stage base entries hBase ∧
      ∀ j : Fin 4, RawDensePolyRep this result.heap
        (hgcdMatPtr original hOriginal j) (hgcdMatLen base hBase j)
        (entries j) := by
  rw [hgcdMatRestoreLoop] at hrun
  split at hrun
  next hlt =>
    dsimp only at hrun
    split at hrun
    next fault hcopy => simp at hrun
    next heap1 hcopy =>
      let index : Fin 4 := ⟨i, hlt⟩
      let dst := hgcdMatPtr original hOriginal index
      let src := stage.add (hgcdMatStageOffset base hBase i)
      let length := hgcdMatLen base hBase index
      have hLength : hgcdMatLen current hCurrent index = length := by
        exact array_getElem_eq_of_eq current.len base.len hLen index.val
          (by rw [hCurrent.2]; exact index.isLt)
          (by rw [hBase.2]; exact index.isLt)
      have hLengthRaw : hgcdMatLenRaw current hCurrent index = length := by
        simpa [hgcdMatLenRaw, hgcdMatLen] using hLength
      have hDst : heap.ValidU64Slice dst length := hws.originalValid index
      have hSrc := hStage index
      have hDstSrc : U64SlicesDisjoint dst length src length := by
        exact u64SlicesDisjoint_stage_entry_right dst stage length
          (hgcdMatStageSize base hBase) (hgcdMatStageOffset base hBase i)
          length (hgcdMatStageOffset_entry_le_size base hBase i hlt)
          (hws.originalStageDisjoint index)
      rcases copyU64_refines_rawDense this heap dst src length (entries index)
          hDst hDstSrc hSrc with
        ⟨semanticHeap, hsemantic, hlayout, hRestoredI⟩
      have hcopyNorm : heap.copyU64 dst src length = .ok heap1 := by
        simpa [dst, src, length, index, hgcdMatPtr, hgcdMatPtrRaw,
          hLengthRaw] using hcopy
      have heq : semanticHeap = heap1 :=
        Except.ok.inj (hsemantic.symm.trans hcopyNorm)
      subst semanticHeap
      have hStage1 : HgcdMatStageRawDenseRep this heap1 stage base entries
          hBase := by
        intro k
        have hDisjoint : U64SlicesDisjoint dst length
            (stage.add (hgcdMatStageOffset base hBase k.val))
            (hgcdMatLen base hBase k) := by
          exact u64SlicesDisjoint_stage_entry_right dst stage length
            (hgcdMatStageSize base hBase)
            (hgcdMatStageOffset base hBase k.val) (hgcdMatLen base hBase k)
            (hgcdMatStageOffset_entry_le_size base hBase k.val k.isLt)
            (hws.originalStageDisjoint index)
        exact (copyU64_preserves_rawDenseRep this heap heap1 dst src length
          (stage.add (hgcdMatStageOffset base hBase k.val))
          (hgcdMatLen base hBase k) (entries k) hDst hSrc.1 hDisjoint
          hcopyNorm (hStage k)).2
      have hPrior1 : ∀ j : Fin 4, j.val < i + 1 →
          RawDensePolyRep this heap1 (hgcdMatPtr original hOriginal j)
            (hgcdMatLen base hBase j) (entries j) := by
        intro j hj
        by_cases hji : j.val = i
        · have hjeq : j = index := Fin.ext hji
          subst j
          exact hRestoredI
        · have hjlt : j.val < i := by omega
          have hDisjoint : U64SlicesDisjoint dst length
              (hgcdMatPtr original hOriginal j)
              (hgcdMatLen base hBase j) := by
            exact hws.originalPairwiseDisjoint index j (by
              intro heqIndex
              exact hji (congrArg Fin.val heqIndex).symm)
          exact (copyU64_preserves_rawDenseRep this heap heap1 dst src length
            (hgcdMatPtr original hOriginal j) (hgcdMatLen base hBase j)
            (entries j) hDst hSrc.1 hDisjoint hcopyNorm
            (hPrior j hjlt)).2
      let poly' := current.poly.set i (hgcdMatPtrRaw original hOriginal index)
        (by rw [hCurrent.1]; omega)
      let next : HgcdMat := { current with poly := poly' }
      have hNext : next.Valid := by
        exact ⟨by simp [next, poly', hCurrent.1], hCurrent.2⟩
      have hLenNext : next.len = base.len := by
        exact hLen
      have hws1 := hgcdMatStabilizeWorkspace_of_sameLayout heap heap1 original
        base hOriginal hBase stage hlayout hws
      have htail : hgcdMatRestoreLoop original next hOriginal hNext stage
          (i + 1) (hgcdMatStageOffset base hBase (i + 1)) heap1 =
            .ok result := by
        rw [hgcdMatStageOffset_step base hBase i hlt]
        simpa [next, poly', index, hLengthRaw] using hrun
      exact hgcdMatRestoreLoop_refines this original base next hOriginal hBase
        hNext stage entries (i + 1) (by omega) heap1 result hLenNext hws1
        hStage1 hPrior1 htail
  next hstop =>
    have hi4 : i = 4 := by omega
    have heq := Except.ok.inj hrun
    subst result
    constructor
    · exact hStage
    · intro j
      exact hPrior j (by omega)
termination_by 4 - i
decreasing_by omega

/-- Every concrete copy in the restore loop is a frame for a represented
live polynomial outside all four saved matrix-entry regions. -/
theorem hgcdMatRestoreLoop_preserves_rawDenseRep (this : DenseUPolyZp)
    (original base current : HgcdMat)
    (hOriginal : original.Valid) (hBase : base.Valid)
    (hCurrent : current.Valid) (stage : RawPtr UInt64)
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (i : Nat) (hi : i ≤ 4) (heap : RawHeap)
    (hLen : current.len = base.len)
    (hws : HgcdMatStabilizeWorkspace heap original base hOriginal hBase stage)
    (hStage : HgcdMatStageRawDenseRep this heap stage base entries hBase)
    (ptr : RawPtr UInt64) (length : Nat)
    (poly : Polynomial (ZMod this._p.toNat))
    (hOriginalPtr : ∀ j : Fin 4,
      (hgcdMatPtr original hOriginal j).region ≠ ptr.region)
    (hRep : RawDensePolyRep this heap ptr length poly)
    (result : HgcdMatRestoreResult)
    (hrun : hgcdMatRestoreLoop original current hOriginal hCurrent stage i
      (hgcdMatStageOffset base hBase i) heap = .ok result) :
    RawHeap.SameLayout heap result.heap ∧
      RawDensePolyRep this result.heap ptr length poly := by
  rw [hgcdMatRestoreLoop] at hrun
  split at hrun
  next hlt =>
    dsimp only at hrun
    split at hrun
    next fault hcopy => simp at hrun
    next heap1 hcopy =>
      let index : Fin 4 := ⟨i, hlt⟩
      let dst := hgcdMatPtr original hOriginal index
      let src := stage.add (hgcdMatStageOffset base hBase i)
      let count := hgcdMatLen base hBase index
      have hLength : hgcdMatLen current hCurrent index = count := by
        exact array_getElem_eq_of_eq current.len base.len hLen index.val
          (by rw [hCurrent.2]; exact index.isLt)
          (by rw [hBase.2]; exact index.isLt)
      have hLengthRaw : hgcdMatLenRaw current hCurrent index = count := by
        simpa [hgcdMatLenRaw, hgcdMatLen] using hLength
      have hDst : heap.ValidU64Slice dst count := hws.originalValid index
      have hSrc := hStage index
      have hLiveDisjoint : U64SlicesDisjoint dst count ptr length := by
        exact u64SlicesDisjoint_of_region_ne (hOriginalPtr index)
      have hcopyNorm : heap.copyU64 dst src count = .ok heap1 := by
        simpa [dst, src, count, index, hgcdMatPtr, hgcdMatPtrRaw,
          hLengthRaw] using hcopy
      have hFrame := copyU64_preserves_rawDenseRep this heap heap1 dst src
        count ptr length poly hDst hSrc.1 hLiveDisjoint hcopyNorm hRep
      have hStage1 : HgcdMatStageRawDenseRep this heap1 stage base entries
          hBase := by
        intro k
        have hDisjoint : U64SlicesDisjoint dst count
            (stage.add (hgcdMatStageOffset base hBase k.val))
            (hgcdMatLen base hBase k) := by
          exact u64SlicesDisjoint_stage_entry_right dst stage count
            (hgcdMatStageSize base hBase)
            (hgcdMatStageOffset base hBase k.val) (hgcdMatLen base hBase k)
            (hgcdMatStageOffset_entry_le_size base hBase k.val k.isLt)
            (hws.originalStageDisjoint index)
        exact (copyU64_preserves_rawDenseRep this heap heap1 dst src count
          (stage.add (hgcdMatStageOffset base hBase k.val))
          (hgcdMatLen base hBase k) (entries k) hDst hSrc.1 hDisjoint
          hcopyNorm (hStage k)).2
      let poly' := current.poly.set i
        (hgcdMatPtrRaw original hOriginal index)
        (by rw [hCurrent.1]; omega)
      let next : HgcdMat := { current with poly := poly' }
      have hNext : next.Valid := by
        exact ⟨by simp [next, poly', hCurrent.1], hCurrent.2⟩
      have hws1 := hgcdMatStabilizeWorkspace_of_sameLayout heap heap1
        original base hOriginal hBase stage hFrame.1 hws
      have hrec := hgcdMatRestoreLoop_preserves_rawDenseRep this original base next
        hOriginal hBase hNext stage entries (i + 1) (by omega) heap1
        (by exact hLen) hws1 hStage1 ptr length poly hOriginalPtr hFrame.2
        result (by
          rw [hgcdMatStageOffset_step base hBase i hlt]
          simpa [next, poly', index, hLengthRaw] using hrun)
      exact ⟨fun p n => (hFrame.1 p n).trans (hrec.1 p n), hrec.2⟩
  next hstop =>
    have heq := Except.ok.inj hrun
    subst result
    exact ⟨fun _ _ => Iff.rfl, hRep⟩
termination_by 4 - i
decreasing_by omega

/-- Entry-point refinement of the second generated loop: staged polynomial
contents are copied back to all four saved pointers, and the returned matrix
descriptor represents those same entries at the iterator-produced lengths. -/
theorem hgcdMatRestoreLoop_zero_refines (this : DenseUPolyZp)
    (original current : HgcdMat)
    (hOriginal : original.Valid) (hCurrent : current.Valid)
    (stage : RawPtr UInt64)
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (heap : RawHeap) (result : HgcdMatRestoreResult)
    (hws : HgcdMatStabilizeWorkspace heap original current hOriginal
      hCurrent stage)
    (hStage : HgcdMatStageRawDenseRep this heap stage current entries hCurrent)
    (hrun : hgcdMatRestoreLoop original current hOriginal hCurrent stage
      0 0 heap = .ok result) :
    ∃ hResult : result.matrix.Valid,
      HgcdMatStageRawDenseRep this result.heap stage current entries hCurrent ∧
      HgcdMatRawDenseRep this result.heap result.matrix entries hResult := by
  have hrefine := hgcdMatRestoreLoop_refines this original current current
    hOriginal hCurrent hCurrent stage entries 0 (by omega) heap result rfl hws
    hStage (by intro j hj; omega) (by
      simpa [hgcdMatStageOffset] using hrun)
  have hdesc := hgcdMatRestoreLoop_zero_descriptors original current hOriginal
    hCurrent stage heap result hrun
  refine ⟨hdesc.1, hrefine.1, ?_⟩
  intro j
  have hPtr : hgcdMatPtr result.matrix hdesc.1 j =
      hgcdMatPtr original hOriginal j := by
    exact array_getElem_eq_of_eq result.matrix.poly original.poly hdesc.2.1
      j.val (by rw [hdesc.1.1]; exact j.isLt)
      (by rw [hOriginal.1]; exact j.isLt)
  have hLength : hgcdMatLen result.matrix hdesc.1 j =
      hgcdMatLen current hCurrent j := by
    exact array_getElem_eq_of_eq result.matrix.len current.len hdesc.2.2
      j.val (by rw [hdesc.1.2]; exact j.isLt)
      (by rw [hCurrent.2]; exact j.isLt)
  rw [hPtr, hLength]
  exact hrefine.2 j

/-- The physical stabilization workspace is sufficient to make every suffix
of the generated restore loop terminate successfully. -/
theorem hgcdMatRestoreLoop_terminates (this : DenseUPolyZp)
    (original base current : HgcdMat)
    (hOriginal : original.Valid) (hBase : base.Valid)
    (hCurrent : current.Valid) (stage : RawPtr UInt64)
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (i : Nat) (hi : i ≤ 4) (heap : RawHeap)
    (hLen : current.len = base.len)
    (hws : HgcdMatStabilizeWorkspace heap original base hOriginal hBase stage)
    (hStage : HgcdMatStageRawDenseRep this heap stage base entries hBase) :
    ∃ result, hgcdMatRestoreLoop original current hOriginal hCurrent stage i
      (hgcdMatStageOffset base hBase i) heap = .ok result := by
  rw [hgcdMatRestoreLoop]
  split
  next hlt =>
    let index : Fin 4 := ⟨i, hlt⟩
    let dst := hgcdMatPtr original hOriginal index
    let src := stage.add (hgcdMatStageOffset base hBase i)
    let length := hgcdMatLen base hBase index
    have hLength : hgcdMatLen current hCurrent index = length := by
      exact array_getElem_eq_of_eq current.len base.len hLen index.val
        (by rw [hCurrent.2]; exact index.isLt)
        (by rw [hBase.2]; exact index.isLt)
    have hLengthRaw : hgcdMatLenRaw current hCurrent index = length := by
      simpa [hgcdMatLenRaw, hgcdMatLen] using hLength
    have hBaseLengthRaw : hgcdMatLenRaw base hBase index = length := by
      rfl
    have hDst : heap.ValidU64Slice dst length := hws.originalValid index
    have hSrc := hStage index
    have hDstSrc : U64SlicesDisjoint dst length src length := by
      exact u64SlicesDisjoint_stage_entry_right dst stage length
        (hgcdMatStageSize base hBase) (hgcdMatStageOffset base hBase i) length
        (hgcdMatStageOffset_entry_le_size base hBase i hlt)
        (hws.originalStageDisjoint index)
    rcases copyU64_refines_rawDense this heap dst src length (entries index)
        hDst hDstSrc hSrc with ⟨heap1, hcopy, hlayout, _⟩
    have hStage1 : HgcdMatStageRawDenseRep this heap1 stage base entries
        hBase := by
      intro k
      have hDisjoint : U64SlicesDisjoint dst length
          (stage.add (hgcdMatStageOffset base hBase k.val))
          (hgcdMatLen base hBase k) := by
        exact u64SlicesDisjoint_stage_entry_right dst stage length
          (hgcdMatStageSize base hBase)
          (hgcdMatStageOffset base hBase k.val) (hgcdMatLen base hBase k)
          (hgcdMatStageOffset_entry_le_size base hBase k.val k.isLt)
          (hws.originalStageDisjoint index)
      exact (copyU64_preserves_rawDenseRep this heap heap1 dst src length
        (stage.add (hgcdMatStageOffset base hBase k.val))
        (hgcdMatLen base hBase k) (entries k) hDst hSrc.1 hDisjoint hcopy
        (hStage k)).2
    let poly' := current.poly.set i (hgcdMatPtrRaw original hOriginal index)
      (by rw [hCurrent.1]; omega)
    let next : HgcdMat := { current with poly := poly' }
    have hNext : next.Valid := by
      exact ⟨by simp [next, poly', hCurrent.1], hCurrent.2⟩
    have hws1 := hgcdMatStabilizeWorkspace_of_sameLayout heap heap1 original
      base hOriginal hBase stage hlayout hws
    rcases hgcdMatRestoreLoop_terminates this original base next hOriginal hBase
        hNext stage entries (i + 1) (by omega) heap1 hLen hws1 hStage1 with
      ⟨result, htail⟩
    refine ⟨result, ?_⟩
    dsimp only
    rw [show heap.copyU64 (hgcdMatPtrRaw original hOriginal index) src
        (hgcdMatLenRaw current hCurrent index) = .ok heap1 by
      simpa [dst, length, hLengthRaw] using hcopy]
    simp only
    rw [hLengthRaw, ← hBaseLengthRaw,
      ← hgcdMatStageOffset_step base hBase i hlt]
    simpa [next, poly', index, hLengthRaw] using htail
  next hstop =>
    exact ⟨HgcdMatRestoreResult.mk heap current
      (hgcdMatStageOffset base hBase i), rfl⟩
termination_by 4 - i
decreasing_by omega

/-- Complete raw semantic refinement of the two-loop source stabilization
block, including successful termination under its physical workspace. -/
theorem hgcdMatStabilize_refines (this : DenseUPolyZp)
    (original current : HgcdMat)
    (hOriginal : original.Valid) (hCurrent : current.Valid)
    (stage : RawPtr UInt64)
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (heap : RawHeap)
    (hws : HgcdMatStabilizeWorkspace heap original current hOriginal
      hCurrent stage)
    (hMatrix : HgcdMatRawDenseRep this heap current entries hCurrent) :
    ∃ result, hgcdMatStabilize original current hOriginal hCurrent stage heap =
        .ok result ∧
      ∃ hResult : result.matrix.Valid,
        HgcdMatRawDenseRep this result.heap result.matrix entries hResult := by
  rcases hgcdMatStageLoop_zero_refines this current original hCurrent hOriginal
      stage entries heap hws hMatrix with
    ⟨staged, hstage, hoff, hlayout, _, hStage⟩
  have hws1 := hgcdMatStabilizeWorkspace_of_sameLayout heap staged.heap original
    current hOriginal hCurrent stage hlayout hws
  rcases hgcdMatRestoreLoop_terminates this original current current hOriginal
      hCurrent hCurrent stage entries 0 (by omega) staged.heap rfl hws1 hStage with
    ⟨result, hrestore⟩
  have hrestore0 : hgcdMatRestoreLoop original current hOriginal hCurrent stage
      0 0 staged.heap = .ok result := by
    simpa [hgcdMatStageOffset] using hrestore
  rcases hgcdMatRestoreLoop_zero_refines this original current hOriginal hCurrent
      stage entries staged.heap result hws1 hStage hrestore0 with
    ⟨hResult, _, hResultRep⟩
  refine ⟨result, ?_, hResult, hResultRep⟩
  simp [hgcdMatStabilize, hstage, hrestore0]

/-- The complete generated matrix stabilization is a frame for a live raw
polynomial outside both the staging allocation and the four saved matrix
entry allocations.  The result also has the same allocation layout. -/
theorem hgcdMatStabilize_preserves_rawDenseRep (this : DenseUPolyZp)
    (original current : HgcdMat)
    (hOriginal : original.Valid) (hCurrent : current.Valid)
    (stage : RawPtr UInt64)
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (heap : RawHeap)
    (hws : HgcdMatStabilizeWorkspace heap original current hOriginal
      hCurrent stage)
    (hMatrix : HgcdMatRawDenseRep this heap current entries hCurrent)
    (ptr : RawPtr UInt64) (length : Nat)
    (poly : Polynomial (ZMod this._p.toNat))
    (hStagePtr : stage.region ≠ ptr.region)
    (hOriginalPtr : ∀ j : Fin 4,
      (hgcdMatPtr original hOriginal j).region ≠ ptr.region)
    (hRep : RawDensePolyRep this heap ptr length poly)
    (result : HgcdMatRestoreResult)
    (hrun : hgcdMatStabilize original current hOriginal hCurrent stage heap =
      .ok result) :
    RawHeap.SameLayout heap result.heap ∧
      RawDensePolyRep this result.heap ptr length poly := by
  simp only [hgcdMatStabilize] at hrun
  generalize hstage : hgcdMatStageLoop current hCurrent stage 0 0 heap = staged
    at hrun
  cases staged with
  | error fault => simp at hrun
  | ok staged =>
    rcases hgcdMatStageLoop_zero_refines this current original hCurrent
        hOriginal stage entries heap hws hMatrix with
      ⟨semantic, hsemantic, _, hlayoutStage, _, hStage⟩
    have hStagedEq : semantic = staged :=
      Except.ok.inj (hsemantic.symm.trans hstage)
    subst semantic
    have hLiveStage := hgcdMatStageLoop_preserves_rawDenseRep this current
      original hCurrent hOriginal stage entries 0 (by omega) heap hws
      hMatrix ptr length poly hStagePtr hRep staged (by
        simpa [hgcdMatStageOffset] using hstage)
    have hws1 := hgcdMatStabilizeWorkspace_of_sameLayout heap staged.heap
      original current hOriginal hCurrent stage hlayoutStage hws
    have hrestore : hgcdMatRestoreLoop original current hOriginal hCurrent
        stage 0 (hgcdMatStageOffset current hCurrent 0) staged.heap =
          .ok result := by
      simpa [hgcdMatStageOffset] using hrun
    have hLiveRestore := hgcdMatRestoreLoop_preserves_rawDenseRep this
      original current current hOriginal hCurrent hCurrent stage entries 0
      (by omega) staged.heap rfl hws1 hStage ptr length poly hOriginalPtr
      hLiveStage.2 result hrestore
    exact ⟨fun p n => (hLiveStage.1 p n).trans (hLiveRestore.1 p n),
      hLiveRestore.2⟩

/-- The source's alias-protection branch really preserves both iterator
outputs: `pB` is copied away first, the `pA` representation is framed across
that write, and only then is `a3` overwritten from `pA`. -/
theorem hgcdRecursiveStoreIterOutputs_cross_refines (this : DenseUPolyZp)
    (a3 b3 pA pB : RawPtr UInt64) (lenA3 lenB3 : Nat)
    (polyA polyB : Polynomial (ZMod this._p.toNat)) (heap : RawHeap)
    (hcross : (!(pA == a3) && (pB == a3)) = true)
    (hA3 : heap.ValidU64Slice a3 lenA3)
    (hB3 : heap.ValidU64Slice b3 lenB3)
    (hARep : RawDensePolyRep this heap pA lenA3 polyA)
    (hBRep : RawDensePolyRep this heap pB lenB3 polyB)
    (hB3PB : U64SlicesDisjoint b3 lenB3 pB lenB3)
    (hB3PA : U64SlicesDisjoint b3 lenB3 pA lenA3)
    (hA3PA : U64SlicesDisjoint a3 lenA3 pA lenA3)
    (hA3B3 : U64SlicesDisjoint a3 lenA3 b3 lenB3) :
    ∃ heap', hgcdRecursiveStoreIterOutputs a3 b3 pA pB lenA3 lenB3 heap =
        .ok heap' ∧
      RawDensePolyRep this heap' a3 lenA3 polyA ∧
      RawDensePolyRep this heap' b3 lenB3 polyB := by
  rcases copyU64_refines_rawDense this heap b3 pB lenB3 polyB hB3 hB3PB
      hBRep with ⟨heap1, hcopyB, hlayout1, hB3Rep⟩
  have hARep1 := (copyU64_preserves_rawDenseRep this heap heap1 b3 pB lenB3
    pA lenA3 polyA hB3 hBRep.1 hB3PA hcopyB hARep).2
  have hA31 : heap1.ValidU64Slice a3 lenA3 :=
    (hlayout1 a3 lenA3).mp hA3
  rcases copyU64_refines_rawDense this heap1 a3 pA lenA3 polyA hA31 hA3PA
      hARep1 with ⟨heap2, hcopyA, _, hA3Rep⟩
  have hB3Rep2 := (copyU64_preserves_rawDenseRep this heap1 heap2 a3 pA lenA3
    b3 lenB3 polyB hA31 hARep1.1 hA3B3 hcopyA hB3Rep).2
  exact ⟨heap2,
    hgcdRecursiveStoreIterOutputs_cross_exec a3 b3 pA pB lenA3 lenB3
      heap heap1 heap2 hcross hcopyB hcopyA,
    hA3Rep, hB3Rep2⟩

/-- In the regular branch where both pointers need normalization, the two
source-order copies also preserve both represented iterator outputs. -/
theorem hgcdRecursiveStoreIterOutputs_regular_both_refines
    (this : DenseUPolyZp)
    (a3 b3 pA pB : RawPtr UInt64) (lenA3 lenB3 : Nat)
    (polyA polyB : Polynomial (ZMod this._p.toNat)) (heap : RawHeap)
    (hcross : (!(pA == a3) && (pB == a3)) = false)
    (hPACopy : (pA == a3) = false) (hPBCopy : (pB == b3) = false)
    (hA3 : heap.ValidU64Slice a3 lenA3)
    (hB3 : heap.ValidU64Slice b3 lenB3)
    (hARep : RawDensePolyRep this heap pA lenA3 polyA)
    (hBRep : RawDensePolyRep this heap pB lenB3 polyB)
    (hA3PA : U64SlicesDisjoint a3 lenA3 pA lenA3)
    (hA3PB : U64SlicesDisjoint a3 lenA3 pB lenB3)
    (hB3PB : U64SlicesDisjoint b3 lenB3 pB lenB3)
    (hB3A3 : U64SlicesDisjoint b3 lenB3 a3 lenA3) :
    ∃ heap', hgcdRecursiveStoreIterOutputs a3 b3 pA pB lenA3 lenB3 heap =
        .ok heap' ∧
      RawDensePolyRep this heap' a3 lenA3 polyA ∧
      RawDensePolyRep this heap' b3 lenB3 polyB := by
  rcases copyU64_refines_rawDense this heap a3 pA lenA3 polyA hA3 hA3PA
      hARep with ⟨heap1, hcopyA, hlayout1, hA3Rep⟩
  have hBRep1 := (copyU64_preserves_rawDenseRep this heap heap1 a3 pA lenA3
    pB lenB3 polyB hA3 hARep.1 hA3PB hcopyA hBRep).2
  have hB31 : heap1.ValidU64Slice b3 lenB3 :=
    (hlayout1 b3 lenB3).mp hB3
  rcases copyU64_refines_rawDense this heap1 b3 pB lenB3 polyB hB31 hB3PB
      hBRep1 with ⟨heap2, hcopyB, _, hB3Rep⟩
  have hA3Rep2 := (copyU64_preserves_rawDenseRep this heap1 heap2 b3 pB lenB3
    a3 lenA3 polyA hB31 hBRep1.1 hB3A3 hcopyB hA3Rep).2
  have hfirst : (if (pA == a3) = false then heap.copyU64 a3 pA lenA3
      else .ok heap) = .ok heap1 := by
    simp [hPACopy, hcopyA]
  have hsecond : (if (pB == b3) = false then heap1.copyU64 b3 pB lenB3
      else .ok heap1) = .ok heap2 := by
    simp [hPBCopy, hcopyB]
  exact ⟨heap2,
    hgcdRecursiveStoreIterOutputs_regular_exec a3 b3 pA pB lenA3 lenB3
      heap heap1 heap2 hcross hfirst hsecond,
    hA3Rep2, hB3Rep⟩

theorem hgcdRecursiveStoreIterOutputs_regular_skip_a_refines
    (this : DenseUPolyZp)
    (a3 b3 pA pB : RawPtr UInt64) (lenA3 lenB3 : Nat)
    (polyA polyB : Polynomial (ZMod this._p.toNat)) (heap : RawHeap)
    (hcross : (!(pA == a3) && (pB == a3)) = false)
    (hPASkip : (pA == a3) = true) (hpA : pA = a3)
    (hPBCopy : (pB == b3) = false)
    (hB3 : heap.ValidU64Slice b3 lenB3)
    (hARep : RawDensePolyRep this heap pA lenA3 polyA)
    (hBRep : RawDensePolyRep this heap pB lenB3 polyB)
    (hB3PB : U64SlicesDisjoint b3 lenB3 pB lenB3)
    (hB3A3 : U64SlicesDisjoint b3 lenB3 a3 lenA3) :
    ∃ heap', hgcdRecursiveStoreIterOutputs a3 b3 pA pB lenA3 lenB3 heap =
        .ok heap' ∧
      RawDensePolyRep this heap' a3 lenA3 polyA ∧
      RawDensePolyRep this heap' b3 lenB3 polyB := by
  have hA3Rep : RawDensePolyRep this heap a3 lenA3 polyA := by
    simpa [hpA] using hARep
  rcases copyU64_refines_rawDense this heap b3 pB lenB3 polyB hB3 hB3PB
      hBRep with ⟨heap1, hcopyB, _, hB3Rep⟩
  have hA3Rep1 := (copyU64_preserves_rawDenseRep this heap heap1 b3 pB lenB3
    a3 lenA3 polyA hB3 hBRep.1 hB3A3 hcopyB hA3Rep).2
  have hfirst : (if (pA == a3) = false then heap.copyU64 a3 pA lenA3
      else .ok heap) = .ok heap := by simp [hPASkip]
  have hsecond : (if (pB == b3) = false then heap.copyU64 b3 pB lenB3
      else .ok heap) = .ok heap1 := by simp [hPBCopy, hcopyB]
  exact ⟨heap1,
    hgcdRecursiveStoreIterOutputs_regular_exec a3 b3 pA pB lenA3 lenB3
      heap heap heap1 hcross hfirst hsecond,
    hA3Rep1, hB3Rep⟩

theorem hgcdRecursiveStoreIterOutputs_regular_skip_b_refines
    (this : DenseUPolyZp)
    (a3 b3 pA pB : RawPtr UInt64) (lenA3 lenB3 : Nat)
    (polyA polyB : Polynomial (ZMod this._p.toNat)) (heap : RawHeap)
    (hcross : (!(pA == a3) && (pB == a3)) = false)
    (hPACopy : (pA == a3) = false)
    (hPBSkip : (pB == b3) = true) (hpB : pB = b3)
    (hA3 : heap.ValidU64Slice a3 lenA3)
    (hARep : RawDensePolyRep this heap pA lenA3 polyA)
    (hBRep : RawDensePolyRep this heap pB lenB3 polyB)
    (hA3PA : U64SlicesDisjoint a3 lenA3 pA lenA3)
    (hA3B3 : U64SlicesDisjoint a3 lenA3 b3 lenB3) :
    ∃ heap', hgcdRecursiveStoreIterOutputs a3 b3 pA pB lenA3 lenB3 heap =
        .ok heap' ∧
      RawDensePolyRep this heap' a3 lenA3 polyA ∧
      RawDensePolyRep this heap' b3 lenB3 polyB := by
  have hB3Rep : RawDensePolyRep this heap b3 lenB3 polyB := by
    simpa [hpB] using hBRep
  rcases copyU64_refines_rawDense this heap a3 pA lenA3 polyA hA3 hA3PA
      hARep with ⟨heap1, hcopyA, _, hA3Rep⟩
  have hB3Rep1 := (copyU64_preserves_rawDenseRep this heap heap1 a3 pA lenA3
    b3 lenB3 polyB hA3 hARep.1 hA3B3 hcopyA hB3Rep).2
  have hfirst : (if (pA == a3) = false then heap.copyU64 a3 pA lenA3
      else .ok heap) = .ok heap1 := by simp [hPACopy, hcopyA]
  have hsecond : (if (pB == b3) = false then heap1.copyU64 b3 pB lenB3
      else .ok heap1) = .ok heap1 := by simp [hPBSkip]
  exact ⟨heap1,
    hgcdRecursiveStoreIterOutputs_regular_exec a3 b3 pA pB lenA3 lenB3
      heap heap1 heap1 hcross hfirst hsecond,
    hA3Rep, hB3Rep1⟩

theorem hgcdRecursiveStoreIterOutputs_regular_skip_both_refines
    (this : DenseUPolyZp)
    (a3 b3 pA pB : RawPtr UInt64) (lenA3 lenB3 : Nat)
    (polyA polyB : Polynomial (ZMod this._p.toNat)) (heap : RawHeap)
    (hcross : (!(pA == a3) && (pB == a3)) = false)
    (hPASkip : (pA == a3) = true) (hpA : pA = a3)
    (hPBSkip : (pB == b3) = true) (hpB : pB = b3)
    (hARep : RawDensePolyRep this heap pA lenA3 polyA)
    (hBRep : RawDensePolyRep this heap pB lenB3 polyB) :
    hgcdRecursiveStoreIterOutputs a3 b3 pA pB lenA3 lenB3 heap = .ok heap ∧
      RawDensePolyRep this heap a3 lenA3 polyA ∧
      RawDensePolyRep this heap b3 lenB3 polyB := by
  have hfirst : (if (pA == a3) = false then heap.copyU64 a3 pA lenA3
      else .ok heap) = .ok heap := by simp [hPASkip]
  have hsecond : (if (pB == b3) = false then heap.copyU64 b3 pB lenB3
      else .ok heap) = .ok heap := by simp [hPBSkip]
  exact ⟨hgcdRecursiveStoreIterOutputs_regular_exec a3 b3 pA pB lenA3 lenB3
      heap heap heap hcross hfirst hsecond,
    by simpa [hpA] using hARep, by simpa [hpB] using hBRep⟩

/-- Complete raw refinement of all pointer-comparison branches in the source
output-normalization block.  Equality bridges are facts about the concrete C++
pointer comparison; the remaining hypotheses are only capacities and aliasing. -/
theorem hgcdRecursiveStoreIterOutputs_refines (this : DenseUPolyZp)
    (a3 b3 pA pB : RawPtr UInt64) (lenA3 lenB3 : Nat)
    (polyA polyB : Polynomial (ZMod this._p.toNat)) (heap : RawHeap)
    (hPAEq : (pA == a3) = true → pA = a3)
    (hPBEq : (pB == b3) = true → pB = b3)
    (hA3 : heap.ValidU64Slice a3 lenA3)
    (hB3 : heap.ValidU64Slice b3 lenB3)
    (hARep : RawDensePolyRep this heap pA lenA3 polyA)
    (hBRep : RawDensePolyRep this heap pB lenB3 polyB)
    (hB3PB : U64SlicesDisjoint b3 lenB3 pB lenB3)
    (hB3PA : U64SlicesDisjoint b3 lenB3 pA lenA3)
    (hA3PA : U64SlicesDisjoint a3 lenA3 pA lenA3)
    (hA3PB : U64SlicesDisjoint a3 lenA3 pB lenB3)
    (hA3B3 : U64SlicesDisjoint a3 lenA3 b3 lenB3) :
    ∃ heap', hgcdRecursiveStoreIterOutputs a3 b3 pA pB lenA3 lenB3 heap =
        .ok heap' ∧
      RawDensePolyRep this heap' a3 lenA3 polyA ∧
      RawDensePolyRep this heap' b3 lenB3 polyB := by
  by_cases hcross : (!(pA == a3) && (pB == a3)) = true
  · exact hgcdRecursiveStoreIterOutputs_cross_refines this a3 b3 pA pB
      lenA3 lenB3 polyA polyB heap hcross hA3 hB3 hARep hBRep hB3PB
      hB3PA hA3PA hA3B3
  · have hcrossFalse : (!(pA == a3) && (pB == a3)) = false := by
      cases hvalue : (!(pA == a3) && (pB == a3)) <;> simp_all
    cases hpa : pA == a3 with
    | false =>
      cases hpb : pB == b3 with
      | false =>
        exact hgcdRecursiveStoreIterOutputs_regular_both_refines this a3 b3
          pA pB lenA3 lenB3 polyA polyB heap hcrossFalse hpa hpb hA3 hB3
          hARep hBRep hA3PA hA3PB hB3PB (u64SlicesDisjoint_symm hA3B3)
      | true =>
        exact hgcdRecursiveStoreIterOutputs_regular_skip_b_refines this a3 b3
          pA pB lenA3 lenB3 polyA polyB heap hcrossFalse hpa hpb (hPBEq hpb)
          hA3 hARep hBRep hA3PA hA3B3
    | true =>
      cases hpb : pB == b3 with
      | false =>
        exact hgcdRecursiveStoreIterOutputs_regular_skip_a_refines this a3 b3
          pA pB lenA3 lenB3 polyA polyB heap hcrossFalse hpa (hPAEq hpa) hpb
          hB3 hARep hBRep hB3PB (u64SlicesDisjoint_symm hA3B3)
      | true =>
        have hboth := hgcdRecursiveStoreIterOutputs_regular_skip_both_refines
          this a3 b3 pA pB lenA3 lenB3 polyA polyB heap hcrossFalse hpa
          (hPAEq hpa) hpb (hPBEq hpb) hARep hBRep
        exact ⟨heap, hboth.1, hboth.2.1, hboth.2.2⟩

/-- Matrix frame for one source-conditional memcpy. -/
theorem optionalCopy_preserves_hgcdMatRawDenseRep (this : DenseUPolyZp)
    (condition : Bool) (heap heap' : RawHeap)
    (dst src : RawPtr UInt64) (count : Nat)
    (M : HgcdMat) (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (hM : M.Valid) (hDst : heap.ValidU64Slice dst count)
    (hSrc : heap.ValidU64Slice src count)
    (hframe : ∀ j : Fin 4, U64SlicesDisjoint dst count
      (hgcdMatPtr M hM j) (hgcdMatLen M hM j))
    (hrep : HgcdMatRawDenseRep this heap M entries hM)
    (hrun : (if condition then heap.copyU64 dst src count else .ok heap) =
      .ok heap') :
    RawHeap.SameLayout heap heap' ∧
      HgcdMatRawDenseRep this heap' M entries hM := by
  cases hcondition : condition with
  | false =>
    have heq : heap = heap' := by simpa [hcondition] using hrun
    subst heap'
    exact ⟨fun _ _ => Iff.rfl, hrep⟩
  | true =>
    have hcopy : heap.copyU64 dst src count = .ok heap' := by
      simpa [hcondition] using hrun
    exact copyU64_preserves_hgcdMatRawDenseRep this heap heap' dst src count
      M entries hM hDst hSrc hframe hcopy hrep

/-- The complete alias-sensitive output normalization is a frame for every
stable matrix entry when both possible destinations are disjoint from it. -/
theorem hgcdRecursiveStoreIterOutputs_preserves_matrix
    (this : DenseUPolyZp)
    (a3 b3 pA pB : RawPtr UInt64) (lenA3 lenB3 : Nat)
    (heap heap' : RawHeap)
    (M : HgcdMat) (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (hM : M.Valid)
    (hA3 : heap.ValidU64Slice a3 lenA3)
    (hB3 : heap.ValidU64Slice b3 lenB3)
    (hPA : heap.ValidU64Slice pA lenA3)
    (hPB : heap.ValidU64Slice pB lenB3)
    (hAFrame : ∀ j : Fin 4, U64SlicesDisjoint a3 lenA3
      (hgcdMatPtr M hM j) (hgcdMatLen M hM j))
    (hBFrame : ∀ j : Fin 4, U64SlicesDisjoint b3 lenB3
      (hgcdMatPtr M hM j) (hgcdMatLen M hM j))
    (hMatrix : HgcdMatRawDenseRep this heap M entries hM)
    (hrun : hgcdRecursiveStoreIterOutputs a3 b3 pA pB lenA3 lenB3 heap =
      .ok heap') :
    RawHeap.SameLayout heap heap' ∧
      HgcdMatRawDenseRep this heap' M entries hM := by
  simp only [hgcdRecursiveStoreIterOutputs] at hrun
  split at hrun
  next hcross =>
    split at hrun
    next fault hcopyB => simp at hrun
    next heap1 hcopyB =>
      have hframe1 := copyU64_preserves_hgcdMatRawDenseRep this heap heap1
        b3 pB lenB3 M entries hM hB3 hPB hBFrame hcopyB hMatrix
      have hA31 := (hframe1.1 a3 lenA3).mp hA3
      have hPA1 := (hframe1.1 pA lenA3).mp hPA
      have hframe2 := copyU64_preserves_hgcdMatRawDenseRep this heap1 heap'
        a3 pA lenA3 M entries hM hA31 hPA1 hAFrame hrun hframe1.2
      exact ⟨fun ptr count =>
        (hframe1.1 ptr count).trans (hframe2.1 ptr count), hframe2.2⟩
  next hregular =>
    generalize hfirst :
      (if !(pA == a3) then heap.copyU64 a3 pA lenA3 else .ok heap) = first
      at hrun
    cases first with
    | error fault => simp at hrun
    | ok heap1 =>
      have hframe1 := optionalCopy_preserves_hgcdMatRawDenseRep this
        (!(pA == a3)) heap heap1 a3 pA lenA3 M entries hM hA3 hPA
        hAFrame hMatrix hfirst
      have hB31 := (hframe1.1 b3 lenB3).mp hB3
      have hPB1 := (hframe1.1 pB lenB3).mp hPB
      have hsecond :
          (if !(pB == b3) then heap1.copyU64 b3 pB lenB3 else .ok heap1) =
            .ok heap' := by
        exact hrun
      have hframe2 := optionalCopy_preserves_hgcdMatRawDenseRep this
        (!(pB == b3)) heap1 heap' b3 pB lenB3 M entries hM hB31 hPB1
        hBFrame hframe1.2 hsecond
      exact ⟨fun ptr count =>
        (hframe1.1 ptr count).trans (hframe2.1 ptr count), hframe2.2⟩

theorem slicePolyRep_zero_length_any (heap : RawHeap) (ptr : RawPtr UInt64)
    (p : Nat) : SlicePolyRep heap ptr 0 p 0 := by
  refine ⟨#[], rfl, rfl, ?_⟩
  ext degree
  rw [coeff_coeffArrayPoly]
  simp

theorem slicePolyRep_one_of_read_one (heap : RawHeap) (ptr : RawPtr UInt64)
    (p : Nat) (hvalid : heap.ValidU64Slice ptr 1)
    (hread : heap.readU64 ptr 0 = .ok 1) :
    SlicePolyRep heap ptr 1 p 1 := by
  rcases slicePolyRep_exists_unique heap ptr 1 p hvalid with
    ⟨poly, hrep, _⟩
  have heq : poly = 1 := by
    ext degree
    by_cases hd : degree = 0
    · subst degree
      rcases slicePolyRep_coeff heap ptr 1 p poly hrep 0 (by omega) with
        ⟨value, hvalue, hcoeff⟩
      have hv : value = 1 := Except.ok.inj (hvalue.symm.trans hread)
      subst value
      simpa using hcoeff
    · rw [slicePolyRep_coeff_zero_of_length_le heap ptr 1 p poly hrep
          degree (by omega)]
      simp [Polynomial.coeff_one, hd]
  simpa [heq] using hrep

/-- A zero-length C++ coefficient slice is the normalized raw
representation of the zero polynomial. -/
theorem rawDensePolyRep_zero_length (this : DenseUPolyZp)
    (heap : RawHeap) (ptr : RawPtr UInt64)
    (hvalid : heap.ValidU64Slice ptr 0) :
    RawDensePolyRep this heap ptr 0 0 := by
  refine ⟨hvalid, ?_, slicePolyRep_zero_length_any heap ptr this._p.toNat, ?_⟩
  · intro k value hk _
    omega
  · simp [RawHeap.normaliseU64]

/-- Physical buffer obligations for one guarded multiplication in the
low-half reconstruction of `_hgcd_recursive`.  The maximum operand length
matches the source's longer-first dispatch and does not encode an L2 value. -/
structure HgcdMulTermWorkspace (heap : RawHeap)
    (dst left : RawPtr UInt64) (lenLeft : Nat)
    (right : RawPtr UInt64) (lenRight : Nat)
    (scratch : RawPtr UInt64) : Prop where
  lengthFits : max lenLeft lenRight < limbBase
  dstValid : heap.ValidU64Slice dst (2 * max lenLeft lenRight - 1)
  scratchValid : heap.ValidU64Slice scratch (8 * max lenLeft lenRight)
  dstLeft : U64SlicesDisjoint dst (2 * max lenLeft lenRight - 1)
    left lenLeft
  dstRight : U64SlicesDisjoint dst (2 * max lenLeft lenRight - 1)
    right lenRight
  dstScratch : U64SlicesDisjoint dst (2 * max lenLeft lenRight - 1)
    scratch (8 * max lenLeft lenRight)
  scratchLeft : U64SlicesDisjoint scratch (8 * max lenLeft lenRight)
    left lenLeft
  scratchRight : U64SlicesDisjoint scratch (8 * max lenLeft lenRight)
    right lenRight

/-- Exact guarded multiplication for canonical fixed-length slices.  The
positive branch retains the source capacity `lenLeft + lenRight - 1`; it does
not claim that this capacity is already the normalized product length. -/
theorem hgcdRecursiveMulTerm_refines_slice (this : DenseUPolyZp)
    (dst left : RawPtr UInt64) (lenLeft : Nat)
    (right : RawPtr UInt64) (lenRight : Nat)
    (scratch : RawPtr UInt64) (heap : RawHeap)
    (leftPoly rightPoly : Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (hwork : HgcdMulTermWorkspace heap dst left lenLeft right lenRight scratch)
    (hLeft : RawCanonicalPolySlice this heap left lenLeft leftPoly)
    (hRight : RawCanonicalPolySlice this heap right lenRight rightPoly) :
    ∃ result, hgcdRecursiveMulTerm this dst left lenLeft right lenRight
        scratch heap = .ok result ∧
      RawHeap.SameLayout heap result.heap ∧
      RawCanonicalPolySlice this result.heap dst result.length
        (leftPoly * rightPoly) := by
  by_cases hLeftPos : 0 < lenLeft
  · by_cases hRightPos : 0 < lenRight
    · by_cases horder : lenRight ≤ lenLeft
      · have hmax : max lenLeft lenRight = lenLeft := Nat.max_eq_left horder
        rcases mul_refines_slice this dst left lenLeft right lenRight scratch
            heap leftPoly rightPoly hcfg hp hLeftPos hRightPos horder
            (by simpa [hmax] using hwork.lengthFits)
            (by simpa [hmax] using hwork.dstValid)
            hLeft.1 hRight.1
            (by simpa [hmax] using hwork.scratchValid)
            (by simpa [hmax] using hwork.dstLeft)
            (by simpa [hmax] using hwork.dstRight)
            (by simpa [hmax] using hwork.dstScratch)
            (by simpa [hmax] using hwork.scratchLeft)
            (by simpa [hmax] using hwork.scratchRight)
            hLeft.2.1 hRight.2.1 hLeft.2.2 hRight.2.2 with
          ⟨heap1, hrun, hlayout, hrep, hcanonical⟩
        refine ⟨{ heap := heap1, length := lenLeft + lenRight - 1 }, ?_,
          hlayout, ?_⟩
        · simp [hgcdRecursiveMulTerm, hLeftPos, hRightPos, horder, hrun]
        · exact ⟨heap1.validU64Slice_mono dst (2 * lenLeft - 1)
              (lenLeft + lenRight - 1)
              ((hlayout dst _).mp (by simpa [hmax] using hwork.dstValid))
              (by omega), hcanonical, hrep⟩
      · have hreverse : lenLeft ≤ lenRight := by omega
        have hmax : max lenLeft lenRight = lenRight :=
          Nat.max_eq_right hreverse
        rcases mul_refines_slice this dst right lenRight left lenLeft scratch
            heap rightPoly leftPoly hcfg hp hRightPos hLeftPos hreverse
            (by simpa [hmax] using hwork.lengthFits)
            (by simpa [hmax] using hwork.dstValid)
            hRight.1 hLeft.1
            (by simpa [hmax] using hwork.scratchValid)
            (by simpa [hmax] using hwork.dstRight)
            (by simpa [hmax] using hwork.dstLeft)
            (by simpa [hmax] using hwork.dstScratch)
            (by simpa [hmax] using hwork.scratchRight)
            (by simpa [hmax] using hwork.scratchLeft)
            hRight.2.1 hLeft.2.1 hRight.2.2 hLeft.2.2 with
          ⟨heap1, hrun, hlayout, hrep, hcanonical⟩
        refine ⟨{ heap := heap1, length := lenLeft + lenRight - 1 }, ?_,
          hlayout, ?_⟩
        · simp [hgcdRecursiveMulTerm, hLeftPos, hRightPos, horder, hrun,
            Nat.add_comm]
        · refine ⟨heap1.validU64Slice_mono dst (2 * lenRight - 1)
              (lenLeft + lenRight - 1)
              ((hlayout dst _).mp (by simpa [hmax] using hwork.dstValid))
              (by omega), ?_, ?_⟩
          · simpa [Nat.add_comm] using hcanonical
          · simpa [Nat.add_comm, mul_comm] using hrep
    · have hlenRight : lenRight = 0 := by omega
      subst lenRight
      have hzero : rightPoly = 0 :=
        slicePolyRep_zero_length heap right this._p.toNat rightPoly hRight.2.2
      subst rightPoly
      refine ⟨{ heap := heap, length := 0 }, ?_, (fun _ _ => Iff.rfl), ?_⟩
      · simp [hgcdRecursiveMulTerm]
      · have hz : RawCanonicalPolySlice this heap dst 0 0 := by
          refine ⟨heap.validU64Slice_mono dst
            (2 * max lenLeft 0 - 1) 0 hwork.dstValid (by omega), ?_, ?_⟩
          · intro k value hk hread
            omega
          · exact slicePolyRep_zero_length_any heap dst this._p.toNat
        simpa using hz
  · have hlenLeft : lenLeft = 0 := by omega
    subst lenLeft
    have hzero : leftPoly = 0 :=
      slicePolyRep_zero_length heap left this._p.toNat leftPoly hLeft.2.2
    subst leftPoly
    refine ⟨{ heap := heap, length := 0 }, ?_, (fun _ _ => Iff.rfl), ?_⟩
    · simp [hgcdRecursiveMulTerm]
    · have hz : RawCanonicalPolySlice this heap dst 0 0 := by
        refine ⟨heap.validU64Slice_mono dst
            (2 * max 0 lenRight - 1) 0 hwork.dstValid (by omega), ?_, ?_⟩
        · intro k value hk hread
          omega
        · exact slicePolyRep_zero_length_any heap dst this._p.toNat
      simpa using hz

/-- The exact guarded C++ multiplication block computes the polynomial
product.  In particular, a zero-length branch is justified directly from
the raw input representation. -/
theorem hgcdRecursiveMulTerm_refines (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (dst left : RawPtr UInt64) (lenLeft : Nat)
    (right : RawPtr UInt64) (lenRight : Nat)
    (scratch : RawPtr UInt64) (heap : RawHeap)
    (leftPoly rightPoly : Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (hwork : HgcdMulTermWorkspace heap dst left lenLeft right lenRight scratch)
    (hLeft : RawDensePolyRep this heap left lenLeft leftPoly)
    (hRight : RawDensePolyRep this heap right lenRight rightPoly) :
    ∃ result, hgcdRecursiveMulTerm this dst left lenLeft right lenRight
        scratch heap = .ok result ∧
      RawHeap.SameLayout heap result.heap ∧
      RawDensePolyRep this result.heap dst result.length
        (leftPoly * rightPoly) := by
  by_cases hLeftPos : 0 < lenLeft
  · by_cases hRightPos : 0 < lenRight
    · by_cases horder : lenRight ≤ lenLeft
      · have hmax : max lenLeft lenRight = lenLeft := Nat.max_eq_left horder
        rcases mul_refines_rawDense this dst left lenLeft right lenRight
            scratch heap leftPoly rightPoly hcfg hp hLeftPos hRightPos horder
            (by simpa [hmax] using hwork.lengthFits)
            (by simpa [hmax] using hwork.dstValid)
            (by simpa [hmax] using hwork.scratchValid)
            (by simpa [hmax] using hwork.dstLeft)
            (by simpa [hmax] using hwork.dstRight)
            (by simpa [hmax] using hwork.dstScratch)
            (by simpa [hmax] using hwork.scratchLeft)
            (by simpa [hmax] using hwork.scratchRight)
            hLeft hRight with ⟨heap1, hrun, hlayout, hrep⟩
        refine ⟨{ heap := heap1, length := lenLeft + lenRight - 1 }, ?_,
          hlayout, hrep⟩
        simp [hgcdRecursiveMulTerm, hLeftPos, hRightPos, horder, hrun]
      · have hreverse : lenLeft ≤ lenRight := by omega
        have hmax : max lenLeft lenRight = lenRight := Nat.max_eq_right hreverse
        rcases mul_refines_rawDense this dst right lenRight left lenLeft
            scratch heap rightPoly leftPoly hcfg hp hRightPos hLeftPos hreverse
            (by simpa [hmax] using hwork.lengthFits)
            (by simpa [hmax] using hwork.dstValid)
            (by simpa [hmax] using hwork.scratchValid)
            (by simpa [hmax] using hwork.dstRight)
            (by simpa [hmax] using hwork.dstLeft)
            (by simpa [hmax] using hwork.dstScratch)
            (by simpa [hmax] using hwork.scratchRight)
            (by simpa [hmax] using hwork.scratchLeft)
            hRight hLeft with ⟨heap1, hrun, hlayout, hrep⟩
        refine ⟨{ heap := heap1, length := lenLeft + lenRight - 1 }, ?_,
          hlayout, ?_⟩
        · simp [hgcdRecursiveMulTerm, hLeftPos, hRightPos, horder, hrun,
            Nat.add_comm]
        · simpa [Nat.add_comm, mul_comm] using hrep
    · have hlenRight : lenRight = 0 := by omega
      subst lenRight
      have hzero : rightPoly = 0 :=
        slicePolyRep_zero_length heap right this._p.toNat rightPoly
          hRight.2.2.1
      subst rightPoly
      refine ⟨{ heap := heap, length := 0 }, ?_, (fun _ _ => Iff.rfl),
        ?_⟩
      · simp [hgcdRecursiveMulTerm]
      · simpa using rawDensePolyRep_zero_length this heap dst
          (heap.validU64Slice_mono dst (2 * max lenLeft 0 - 1) 0
            hwork.dstValid (by omega))
  · have hlenLeft : lenLeft = 0 := by omega
    subst lenLeft
    have hzero : leftPoly = 0 :=
      slicePolyRep_zero_length heap left this._p.toNat leftPoly
        hLeft.2.2.1
    subst leftPoly
    refine ⟨{ heap := heap, length := 0 }, ?_, (fun _ _ => Iff.rfl),
      ?_⟩
    · simp [hgcdRecursiveMulTerm]
    · simpa using rawDensePolyRep_zero_length this heap dst
        (heap.validU64Slice_mono dst (2 * max 0 lenRight - 1) 0
          hwork.dstValid (by omega))

/-- Exact product-capacity bound for the generated guarded multiplication.
Unlike the older coarse sum bound, this retains the source's `- 1`, including
the zero-input branches where truncated subtraction still makes the bound
valid. -/
theorem hgcdRecursiveMulTerm_length_le_product (this : DenseUPolyZp)
    (dst left : RawPtr UInt64) (lenLeft : Nat)
    (right : RawPtr UInt64) (lenRight : Nat)
    (scratch : RawPtr UInt64) (heap : RawHeap) (result : HgcdMulTermResult)
    (hrun : hgcdRecursiveMulTerm this dst left lenLeft right lenRight
      scratch heap = .ok result) :
    result.length ≤ lenLeft + lenRight - 1 := by
  simp only [hgcdRecursiveMulTerm] at hrun
  split at hrun
  next hnonzero =>
    split at hrun
    next fault hmul => simp at hrun
    next heap1 hmul =>
      have heq : result =
          HgcdMulTermResult.mk heap1 (lenLeft + lenRight - 1) :=
        (Except.ok.inj hrun).symm
      subst result
      exact Nat.le_refl _
  next hzero =>
    have heq : result = HgcdMulTermResult.mk heap 0 :=
      (Except.ok.inj hrun).symm
    subst result
    exact Nat.zero_le _

/-- Frame rule for the exact guarded multiplication block.  It exposes that
the generated execution writes only its destination and scratch areas, which
is needed to sequence the two products in each reconstruction block. -/
theorem hgcdRecursiveMulTerm_preserves_guard (this : DenseUPolyZp)
    (dst left : RawPtr UInt64) (lenLeft : Nat)
    (right : RawPtr UInt64) (lenRight : Nat)
    (scratch guard : RawPtr UInt64) (guardLen : Nat)
    (heap : RawHeap) (result : HgcdMulTermResult)
    (hwork : HgcdMulTermWorkspace heap dst left lenLeft right lenRight scratch)
    (hLeft : heap.ValidU64Slice left lenLeft)
    (hRight : heap.ValidU64Slice right lenRight)
    (hDstGuard : U64SlicesDisjoint dst (2 * max lenLeft lenRight - 1)
      guard guardLen)
    (hScratchGuard : U64SlicesDisjoint scratch
      (8 * max lenLeft lenRight) guard guardLen)
    (hrun : hgcdRecursiveMulTerm this dst left lenLeft right lenRight
      scratch heap = .ok result) :
    SameU64Prefix heap result.heap guard guardLen := by
  by_cases hLeftPos : 0 < lenLeft
  · by_cases hRightPos : 0 < lenRight
    · by_cases horder : lenRight ≤ lenLeft
      · have hmax : max lenLeft lenRight = lenLeft := Nat.max_eq_left horder
        cases hmul : Generated.StrictMul.dense_upoly_zp__mul_ir this dst left lenLeft right
            lenRight scratch heap with
        | error fault =>
          simp [hgcdRecursiveMulTerm, hLeftPos, hRightPos, horder, hmul]
            at hrun
        | ok heap1 =>
          have hactual :
              (Except.ok (HgcdMulTermResult.mk heap1
                (lenLeft + lenRight - 1)) : RawExec HgcdMulTermResult) =
              Except.ok result := by
            simpa [hgcdRecursiveMulTerm, hLeftPos, hRightPos, horder, hmul]
              using hrun
          have heq : result =
              HgcdMulTermResult.mk heap1 (lenLeft + lenRight - 1) :=
            (Except.ok.inj hactual).symm
          subst result
          exact mul_preserves_prefix this dst left lenLeft right lenRight
            scratch guard guardLen heap heap1 hLeftPos hRightPos horder
            (by simpa [hmax] using hwork.dstValid) hLeft hRight
            (by simpa [hmax] using hwork.scratchValid)
            (by simpa [hmax] using hwork.scratchRight)
            (by simpa [hmax] using hDstGuard)
            (by simpa [hmax] using hScratchGuard) hmul
      · have hreverse : lenLeft ≤ lenRight := by omega
        have hmax : max lenLeft lenRight = lenRight := Nat.max_eq_right hreverse
        cases hmul : Generated.StrictMul.dense_upoly_zp__mul_ir this dst right lenRight left
            lenLeft scratch heap with
        | error fault =>
          simp [hgcdRecursiveMulTerm, hLeftPos, hRightPos, horder, hmul]
            at hrun
        | ok heap1 =>
          have hactual :
              (Except.ok (HgcdMulTermResult.mk heap1
                (lenLeft + lenRight - 1)) : RawExec HgcdMulTermResult) =
              Except.ok result := by
            simpa [hgcdRecursiveMulTerm, hLeftPos, hRightPos, horder, hmul]
              using hrun
          have heq : result =
              HgcdMulTermResult.mk heap1 (lenLeft + lenRight - 1) :=
            (Except.ok.inj hactual).symm
          subst result
          exact mul_preserves_prefix this dst right lenRight left lenLeft
            scratch guard guardLen heap heap1 hRightPos hLeftPos hreverse
            (by simpa [hmax] using hwork.dstValid) hRight hLeft
            (by simpa [hmax] using hwork.scratchValid)
            (by simpa [hmax] using hwork.scratchLeft)
            (by simpa [hmax] using hDstGuard)
            (by simpa [hmax] using hScratchGuard) hmul
    · have hlenRight : lenRight = 0 := by omega
      subst lenRight
      have hactual :
          (Except.ok (HgcdMulTermResult.mk heap 0) :
            RawExec HgcdMulTermResult) = Except.ok result := by
        simpa [hgcdRecursiveMulTerm] using hrun
      have heq : result = HgcdMulTermResult.mk heap 0 :=
        (Except.ok.inj hactual).symm
      subst result
      exact fun _ _ _ hread => hread
  · have hlenLeft : lenLeft = 0 := by omega
    subst lenLeft
    have hactual :
        (Except.ok (HgcdMulTermResult.mk heap 0) :
          RawExec HgcdMulTermResult) = Except.ok result := by
      simpa [hgcdRecursiveMulTerm] using hrun
    have heq : result = HgcdMulTermResult.mk heap 0 :=
      (Except.ok.inj hactual).symm
    subst result
    exact fun _ _ _ hread => hread

def hgcdMulCapacity (leftLength rightLength : Nat) : Nat :=
  2 * max leftLength rightLength - 1

/-- Physical obligations for the two generated products and their in-place
subtraction.  It contains capacities and aliasing only, never an expected L2
polynomial or a supplied execution result. -/
structure HgcdReconstructWorkspace (heap : RawHeap)
    (out temp left1 : RawPtr UInt64) (lenLeft1 : Nat)
    (right1 : RawPtr UInt64) (lenRight1 : Nat)
    (left2 : RawPtr UInt64) (lenLeft2 : Nat)
    (right2 : RawPtr UInt64) (lenRight2 : Nat)
    (scratch : RawPtr UInt64) : Prop where
  first : HgcdMulTermWorkspace heap out left1 lenLeft1 right1 lenRight1 scratch
  second : HgcdMulTermWorkspace heap temp left2 lenLeft2 right2 lenRight2 scratch
  subDstValid : heap.ValidU64Slice out
    (max (hgcdMulCapacity lenLeft1 lenRight1)
      (hgcdMulCapacity lenLeft2 lenRight2))
  firstDstLeft2 : U64SlicesDisjoint out
    (hgcdMulCapacity lenLeft1 lenRight1) left2 lenLeft2
  firstScratchLeft2 : U64SlicesDisjoint scratch
    (8 * max lenLeft1 lenRight1) left2 lenLeft2
  firstDstRight2 : U64SlicesDisjoint out
    (hgcdMulCapacity lenLeft1 lenRight1) right2 lenRight2
  firstScratchRight2 : U64SlicesDisjoint scratch
    (8 * max lenLeft1 lenRight1) right2 lenRight2
  secondDstFirst : U64SlicesDisjoint temp
    (hgcdMulCapacity lenLeft2 lenRight2) out
    (hgcdMulCapacity lenLeft1 lenRight1)
  secondScratchFirst : U64SlicesDisjoint scratch
    (8 * max lenLeft2 lenRight2) out
    (hgcdMulCapacity lenLeft1 lenRight1)
  subAliasTemp : ExactOrDisjoint out temp

theorem hgcdMulTermWorkspace_of_sameLayout (heap heap' : RawHeap)
    (dst left : RawPtr UInt64) (lenLeft : Nat)
    (right : RawPtr UInt64) (lenRight : Nat) (scratch : RawPtr UInt64)
    (hlayout : RawHeap.SameLayout heap heap')
    (hwork : HgcdMulTermWorkspace heap dst left lenLeft right lenRight scratch) :
    HgcdMulTermWorkspace heap' dst left lenLeft right lenRight scratch := by
  exact ⟨hwork.lengthFits,
    (hlayout dst (2 * max lenLeft lenRight - 1)).mp hwork.dstValid,
    (hlayout scratch (8 * max lenLeft lenRight)).mp hwork.scratchValid,
    hwork.dstLeft, hwork.dstRight, hwork.dstScratch,
    hwork.scratchLeft, hwork.scratchRight⟩

theorem rawCanonicalPolySlice_of_same_prefix (this : DenseUPolyZp)
    (before after : RawHeap) (ptr : RawPtr UInt64) (length : Nat)
    (poly : Polynomial (ZMod this._p.toNat))
    (hlayout : RawHeap.SameLayout before after)
    (hprefix : SameU64Prefix before after ptr length)
    (hrep : RawCanonicalPolySlice this before ptr length poly) :
    RawCanonicalPolySlice this after ptr length poly :=
  ⟨(hlayout ptr length).mp hrep.1,
    canonicalU64Prefix_of_same_prefix before after ptr length this._p
      hrep.1 hprefix hrep.2.1,
    slicePolyRep_of_same_prefix before after ptr length this._p.toNat poly
      hrep.1 ((hlayout ptr length).mp hrep.1) hprefix hrep.2.2⟩

/-- Semantic refinement of the exact `b2` low reconstruction block.  Both
products and the final sign-selected subtraction are actual generated raw
executions. -/
theorem hgcdRecursiveReconstructB_refines (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (b2 T0 r2 r0 aLo bLo scratch : RawPtr UInt64)
    (lenR2 lenR0 lenALo lenBLo : Nat) (sgn : Int) (heap : RawHeap)
    (polyR2 polyR0 polyALo polyBLo :
      Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (hwork : HgcdReconstructWorkspace heap b2 T0 r2 lenR2 aLo lenALo
      r0 lenR0 bLo lenBLo scratch)
    (hR2 : RawDensePolyRep this heap r2 lenR2 polyR2)
    (hR0 : RawDensePolyRep this heap r0 lenR0 polyR0)
    (hALo : RawCanonicalPolySlice this heap aLo lenALo polyALo)
    (hBLo : RawCanonicalPolySlice this heap bLo lenBLo polyBLo) :
    ∃ heap' length,
      hgcdRecursiveReconstructB this b2 T0 r2 r0 aLo bLo scratch
        lenR2 lenR0 lenALo lenBLo sgn heap = .ok (heap', length) ∧
      RawHeap.SameLayout heap heap' ∧
      RawDensePolyRep this heap' b2 length
        (if sgn < 0 then polyR2 * polyALo - polyR0 * polyBLo
         else polyR0 * polyBLo - polyR2 * polyALo) ∧
      length ≤ max (lenR2 + lenALo - 1) (lenR0 + lenBLo - 1) := by
  rcases hgcdRecursiveMulTerm_refines_slice this b2 r2 lenR2 aLo lenALo
      scratch heap polyR2 polyALo hcfg hp hwork.first
      hR2.toCanonicalSlice hALo with
    ⟨term1, hrun1, hlayout1, hTerm1⟩
  have hR01 : RawDensePolyRep this term1.heap r0 lenR0 polyR0 := by
    apply rawDensePolyRep_of_same_prefix this heap term1.heap r0 lenR0
      polyR0 hlayout1
    · exact hgcdRecursiveMulTerm_preserves_guard this b2 r2 lenR2 aLo
        lenALo scratch r0 lenR0 heap term1 hwork.first hR2.1 hALo.1
        hwork.firstDstLeft2 hwork.firstScratchLeft2 hrun1
    · exact hR0
  have hBLo1 : RawCanonicalPolySlice this term1.heap bLo lenBLo polyBLo := by
    apply rawCanonicalPolySlice_of_same_prefix this heap term1.heap bLo
      lenBLo polyBLo hlayout1
    · exact hgcdRecursiveMulTerm_preserves_guard this b2 r2 lenR2 aLo
        lenALo scratch bLo lenBLo heap term1 hwork.first hR2.1 hALo.1
        hwork.firstDstRight2 hwork.firstScratchRight2 hrun1
    · exact hBLo
  have hSecond1 := hgcdMulTermWorkspace_of_sameLayout heap term1.heap T0 r0
    lenR0 bLo lenBLo scratch hlayout1 hwork.second
  rcases hgcdRecursiveMulTerm_refines_slice this T0 r0 lenR0 bLo lenBLo
      scratch term1.heap polyR0 polyBLo hcfg hp hSecond1
      hR01.toCanonicalSlice hBLo1 with
    ⟨term2, hrun2, hlayout2, hTerm2⟩
  have hLen1 := hgcdRecursiveMulTerm_length_le this b2 r2 lenR2 aLo
    lenALo scratch heap term1 hrun1
  have hLen2 := hgcdRecursiveMulTerm_length_le this T0 r0 lenR0 bLo
    lenBLo scratch term1.heap term2 hrun2
  have hLen1Product := hgcdRecursiveMulTerm_length_le_product this b2 r2 lenR2 aLo
    lenALo scratch heap term1 hrun1
  have hLen2Product := hgcdRecursiveMulTerm_length_le_product this T0 r0 lenR0 bLo
    lenBLo scratch term1.heap term2 hrun2
  have hTerm1Final : RawCanonicalPolySlice this term2.heap b2 term1.length
      (polyR2 * polyALo) := by
    apply rawCanonicalPolySlice_of_same_prefix this term1.heap term2.heap b2
      term1.length (polyR2 * polyALo) hlayout2
    · apply hgcdRecursiveMulTerm_preserves_guard this T0 r0 lenR0 bLo
        lenBLo scratch b2 term1.length term1.heap term2 hSecond1 hR01.1
        hBLo1.1
      · exact u64SlicesDisjoint_mono hwork.secondDstFirst
          (Nat.le_refl _) hLen1
      · exact u64SlicesDisjoint_mono hwork.secondScratchFirst
          (Nat.le_refl _) hLen1
      · exact hrun2
    · exact hTerm1
  have hSubValid0 : heap.ValidU64Slice b2
      (max term1.length term2.length) :=
    heap.validU64Slice_mono b2
      (max (hgcdMulCapacity lenR2 lenALo)
        (hgcdMulCapacity lenR0 lenBLo))
      (max term1.length term2.length) hwork.subDstValid (by
        apply max_le
        · exact le_max_of_le_left hLen1
        · exact le_max_of_le_right hLen2)
  have hSubValid : term2.heap.ValidU64Slice b2
      (max term1.length term2.length) :=
    (hlayout2 b2 _).mp ((hlayout1 b2 _).mp hSubValid0)
  have hpWord : this._p ≠ 0 := by
    intro hzero
    have hzeroNat := congrArg UInt64.toNat hzero
    simp at hzeroNat
    omega
  by_cases hsgn : sgn < 0
  · rcases polySub_ok this b2 b2 term1.length T0 term2.length
        term2.heap hSubValid hTerm1Final.1 hTerm2.1 with
      ⟨heap3, length, hsub, hlayout3, hlength⟩
    have hrep := polySub_refines this b2 b2 term1.length T0 term2.length
      term2.heap heap3 length (polyR2 * polyALo) (polyR0 * polyBLo)
      hpWord hSubValid hTerm1Final hTerm2 (Or.inl rfl)
      hwork.subAliasTemp hsub
    refine ⟨heap3, length, ?_, ?_, ?_, ?_⟩
    · simp [hgcdRecursiveReconstructB, hrun1, hrun2, hsgn, hsub]
    · exact fun ptr count =>
        (hlayout1 ptr count).trans
          ((hlayout2 ptr count).trans (hlayout3 ptr count))
    · simpa [hsgn] using hrep
    · exact hlength.trans (max_le_max hLen1Product hLen2Product)

  · rcases polySub_ok this b2 T0 term2.length b2 term1.length
        term2.heap (by simpa [max_comm] using hSubValid) hTerm2.1
        hTerm1Final.1 with
      ⟨heap3, length, hsub, hlayout3, hlength⟩
    have hrep := polySub_refines this b2 T0 term2.length b2 term1.length
      term2.heap heap3 length (polyR0 * polyBLo) (polyR2 * polyALo)
      hpWord (by simpa [max_comm] using hSubValid) hTerm2 hTerm1Final
      hwork.subAliasTemp
      (Or.inl rfl) hsub
    refine ⟨heap3, length, ?_, ?_, ?_, ?_⟩
    · simp [hgcdRecursiveReconstructB, hrun1, hrun2, hsgn, hsub]
    · exact fun ptr count =>
        (hlayout1 ptr count).trans
          ((hlayout2 ptr count).trans (hlayout3 ptr count))
    · simpa [hsgn] using hrep
    · simpa [max_comm] using hlength.trans
        (max_le_max hLen2Product hLen1Product)

/-- Semantic refinement of the exact `a2` reconstruction block.  Its source
differs from `b2` only by reversing the sign-selected subtraction, so the
already verified physical execution is reused with an explicitly flipped
branch selector. -/
theorem hgcdRecursiveReconstructA_refines (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (a2 T0 r3 r1 aLo bLo scratch : RawPtr UInt64)
    (lenR3 lenR1 lenALo lenBLo : Nat) (sgn : Int) (heap : RawHeap)
    (polyR3 polyR1 polyALo polyBLo :
      Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (hwork : HgcdReconstructWorkspace heap a2 T0 r3 lenR3 aLo lenALo
      r1 lenR1 bLo lenBLo scratch)
    (hR3 : RawDensePolyRep this heap r3 lenR3 polyR3)
    (hR1 : RawDensePolyRep this heap r1 lenR1 polyR1)
    (hALo : RawCanonicalPolySlice this heap aLo lenALo polyALo)
    (hBLo : RawCanonicalPolySlice this heap bLo lenBLo polyBLo) :
    ∃ heap' length,
      hgcdRecursiveReconstructA this a2 T0 r3 r1 aLo bLo scratch
        lenR3 lenR1 lenALo lenBLo sgn heap = .ok (heap', length) ∧
      RawHeap.SameLayout heap heap' ∧
      RawDensePolyRep this heap' a2 length
        (if sgn < 0 then polyR1 * polyBLo - polyR3 * polyALo
         else polyR3 * polyALo - polyR1 * polyBLo) ∧
      length ≤ max (lenR3 + lenALo - 1) (lenR1 + lenBLo - 1) := by
  let flippedSign : Int := if sgn < 0 then 0 else -1
  rcases hgcdRecursiveReconstructB_refines this a2 T0 r3 r1 aLo bLo
      scratch lenR3 lenR1 lenALo lenBLo flippedSign heap polyR3 polyR1
      polyALo polyBLo hcfg hp hwork hR3 hR1 hALo hBLo with
    ⟨heap', length, hrun, hlayout, hrep, hlength⟩
  refine ⟨heap', length, ?_, hlayout, ?_, hlength⟩
  · have hfunctions :
        hgcdRecursiveReconstructA this a2 T0 r3 r1 aLo bLo scratch
            lenR3 lenR1 lenALo lenBLo sgn heap =
          hgcdRecursiveReconstructB this a2 T0 r3 r1 aLo bLo scratch
            lenR3 lenR1 lenALo lenBLo flippedSign heap := by
      by_cases hsgn : sgn < 0
      · simp [hgcdRecursiveReconstructA, hgcdRecursiveReconstructB,
          flippedSign, hsgn]
      · simp [hgcdRecursiveReconstructA, hgcdRecursiveReconstructB,
          flippedSign, hsgn]
    exact hfunctions.trans hrun
  · by_cases hsgn : sgn < 0
    · simpa [flippedSign, hsgn] using hrep
    · simpa [flippedSign, hsgn] using hrep

/-- Purely physical buffers for the generated shifted-high reconstruction.
The output covers both the existing low representation and the inserted high
half; the alias condition is the one accepted by the real in-place add. -/
structure HgcdLiftHighWorkspace (heap : RawHeap)
    (out high : RawPtr UInt64) (lowLength m highLength : Nat) : Prop where
  outValid : heap.ValidU64Slice out (max (m + highLength) lowLength)
  highValid : heap.ValidU64Slice high highLength
  zeroHighDisjoint : U64SlicesDisjoint out
    (max (m + highLength) lowLength) high highLength
  addAliasHigh : ExactOrDisjoint (out.add m) high

/-- Splitting a normalized raw polynomial at an in-range offset leaves a
normalized suffix.  This is what permits the source's in-place `_poly_add`
on `out + m` without inventing a normalized padded whole buffer. -/
theorem rawDensePolyRep_split_suffix (this : DenseUPolyZp)
    (heap : RawHeap) (ptr : RawPtr UInt64) (length m : Nat)
    (poly : Polynomial (ZMod this._p.toNat)) (hm : m ≤ length)
    (hrep : RawDensePolyRep this heap ptr length poly) :
    ∃ low high : Polynomial (ZMod this._p.toNat),
      SlicePolyRep heap ptr m this._p.toNat low ∧
      RawDensePolyRep this heap (ptr.add m) (length - m) high ∧
      poly = low + Polynomial.X ^ m * high := by
  rcases slicePolyRep_split_exists heap ptr m (length - m) this._p.toNat
      poly (by simpa [Nat.add_sub_of_le hm] using hrep.1)
      (by simpa [Nat.add_sub_of_le hm] using hrep.2.2.1) with
    ⟨low, high, hLow, hHigh, hsplit⟩
  have hcanonical := canonicalU64Prefix_split heap ptr m (length - m)
    this._p (by simpa [Nat.add_sub_of_le hm] using hrep.2.1)
  have hvalidHigh := heap.validU64Slice_add ptr length m (length - m)
    hrep.1 (by omega)
  by_cases hstrict : m < length
  · have hnormHigh : heap.normaliseU64 (ptr.add m) (length - m) =
        .ok (length - m) := by
      rcases normaliseU64_spec heap ptr length hrep.1 with
        ⟨observed, hobserved, _, _, hlast⟩
      have heq : observed = length :=
        Except.ok.inj (hobserved.symm.trans hrep.2.2.2)
      subst observed
      rcases hlast with hzero | ⟨value, hread, hne⟩
      · omega
      · have hreadHigh : heap.readU64 (ptr.add m) (length - m - 1) =
            .ok value := by
          rw [RawHeap.readU64_add]
          have hindex : m + (length - m - 1) = length - 1 := by omega
          simpa [hindex] using hread
        rw [show length - m = (length - m - 1) + 1 by omega]
        simp [RawHeap.normaliseU64, hreadHigh, hne]
    exact ⟨low, high, hLow,
      ⟨hvalidHigh, hcanonical.2, hHigh, hnormHigh⟩, hsplit⟩
  · have heq : m = length := by omega
    subst m
    have hhighZero : high = 0 :=
      slicePolyRep_zero_length heap (ptr.add length) this._p.toNat high
        (by simpa using hHigh)
    subst high
    have hzeroRep : RawDensePolyRep this heap (ptr.add length)
        (length - length) 0 := by
      simpa using rawDensePolyRep_zero_length this heap (ptr.add length)
        (by simpa using hvalidHigh)
    exact ⟨low, 0, hLow, hzeroRep, by simpa using hsplit⟩

/-- Semantic result of the conditional source `memset`: the original
normalized prefix is framed, while the enlarged physical prefix represents
the same polynomial with canonical zero limbs appended. -/
theorem hgcdLiftHigh_zero_refines (this : DenseUPolyZp)
    (out : RawPtr UInt64) (lowLength m highLength : Nat)
    (heap : RawHeap) (lowPoly : Polynomial (ZMod this._p.toNat))
    (hp : this._p ≠ 0)
    (hOut : heap.ValidU64Slice out (max (m + highLength) lowLength))
    (hLow : RawDensePolyRep this heap out lowLength lowPoly) :
    ∃ heap1,
      (if lowLength < m + highLength then
        Generated.StrictMul.mulZeroPadLoop out lowLength
          (m + highLength - lowLength) 0 heap
       else (Except.ok heap : RawExec RawHeap)) = Except.ok heap1 ∧
      RawHeap.SameLayout heap heap1 ∧
      RawDensePolyRep this heap1 out lowLength lowPoly ∧
      SlicePolyRep heap1 out (max (m + highLength) lowLength)
        this._p.toNat lowPoly ∧
      CanonicalU64Prefix heap1 out (max (m + highLength) lowLength)
        this._p := by
  by_cases hpad : lowLength < m + highLength
  · have hsum : lowLength + (m + highLength - lowLength) =
        m + highLength := by omega
    have hmax : max (m + highLength) lowLength = m + highLength :=
      Nat.max_eq_left (Nat.le_of_lt hpad)
    have hvalid : heap.ValidU64Slice out
        (lowLength + (m + highLength - lowLength)) := by
      simpa [hsum, hmax] using hOut
    rcases mulZeroPadLoop_refines out lowLength
        (m + highLength - lowLength) 0 this._p.toNat heap lowPoly this._p
        (Nat.zero_le _) hp hvalid (by simpa using hLow.2.2.1)
        (by simpa using hLow.2.1) with
      ⟨heap1, hrun, hlayout, hslice, hcanonical⟩
    have hsame := mulZeroPadLoop_preserves_before_start out lowLength
      (m + highLength - lowLength) 0 lowLength heap heap1
      (Nat.le_refl _) hvalid hrun
    have hLow1 := rawDensePolyRep_of_same_prefix this heap heap1 out
      lowLength lowPoly hlayout hsame hLow
    exact ⟨heap1, by simpa [hpad] using hrun, hlayout, hLow1,
      by simpa [hsum, hmax] using hslice,
      by simpa [hsum, hmax] using hcanonical⟩
  · have hmax : max (m + highLength) lowLength = lowLength :=
      Nat.max_eq_right (by omega)
    exact ⟨heap, by simp [hpad], fun _ _ => Iff.rfl, hLow,
      by simpa [hmax] using hLow.2.2.1,
      by simpa [hmax] using hLow.2.1⟩

/-- Re-expand a normalized raw result to the full prefix scanned by the
actual normalization call.  The discarded suffix consists of physical zero
limbs, so both the polynomial and canonical-residue properties extend. -/
theorem rawDensePolyRep_extend_to_normalise_input (this : DenseUPolyZp)
    (heap : RawHeap) (ptr : RawPtr UInt64) (fullLength resultLength : Nat)
    (poly : Polynomial (ZMod this._p.toNat)) (hp : this._p ≠ 0)
    (hFull : heap.ValidU64Slice ptr fullLength)
    (hResult : RawDensePolyRep this heap ptr resultLength poly)
    (hnorm : heap.normaliseU64 ptr fullLength = .ok resultLength) :
    SlicePolyRep heap ptr fullLength this._p.toNat poly ∧
      CanonicalU64Prefix heap ptr fullLength this._p := by
  rcases normaliseU64_spec heap ptr fullLength hFull with
    ⟨observed, hobserved, hle, hzeros, _⟩
  have heq : observed = resultLength :=
    Except.ok.inj (hobserved.symm.trans hnorm)
  subst observed
  rcases slicePolyRep_exists_unique heap ptr fullLength this._p.toNat hFull with
    ⟨actual, hActual, _⟩
  have hpoly : actual = poly := by
    ext degree
    by_cases hshort : degree < resultLength
    · rcases slicePolyRep_coeff heap ptr fullLength this._p.toNat actual
          hActual degree (by omega) with ⟨v1, hr1, hc1⟩
      rcases slicePolyRep_coeff heap ptr resultLength this._p.toNat poly
          hResult.2.2.1 degree hshort with ⟨v2, hr2, hc2⟩
      have hv : v1 = v2 := Except.ok.inj (hr1.symm.trans hr2)
      simpa [hc1, hc2, hv]
    · have hpolyZero := slicePolyRep_coeff_zero_of_length_le heap ptr
          resultLength this._p.toNat poly hResult.2.2.1 degree (by omega)
      by_cases hfull : degree < fullLength
      · rcases slicePolyRep_coeff heap ptr fullLength this._p.toNat actual
            hActual degree hfull with ⟨value, hread, hcoeff⟩
        have hreadZero := hzeros degree (by omega) hfull
        have hv : value = 0 := Except.ok.inj (hread.symm.trans hreadZero)
        subst value
        simpa [hcoeff] using hpolyZero.symm
      · rw [slicePolyRep_coeff_zero_of_length_le heap ptr fullLength
            this._p.toNat actual hActual degree (by omega), hpolyZero]
  have hSlice : SlicePolyRep heap ptr fullLength this._p.toNat poly := by
    simpa [hpoly] using hActual
  refine ⟨hSlice, ?_⟩
  intro k value hk hread
  by_cases hshort : k < resultLength
  · exact hResult.2.1 k value hshort hread
  · have hreadZero := hzeros k (by omega) hk
    have hv : value = 0 := Except.ok.inj (hread.symm.trans hreadZero)
    subst value
    have hpPos : 0 < this._p.toNat := UInt64.pos_iff_ne_zero.mpr hp
    simpa using hpPos

/-- Total raw execution bridge for zero-fill, shifted in-place addition, and
whole-buffer normalization.  Every step is the corresponding generated L1
operation; no expected L2 output is an input to the execution. -/
theorem hgcdRecursiveLiftHigh_terminates (this : DenseUPolyZp)
    (out high : RawPtr UInt64) (lowLength m highLength : Nat)
    (heap : RawHeap) (lowPoly : Polynomial (ZMod this._p.toNat))
    (hp : this._p ≠ 0)
    (hwork : HgcdLiftHighWorkspace heap out high lowLength m highLength)
    (hLow : RawDensePolyRep this heap out lowLength lowPoly) :
    ∃ result,
      hgcdRecursiveLiftHigh this out high lowLength m highLength heap =
        .ok result ∧
      RawHeap.SameLayout heap result.heap ∧
      result.length ≤ max (m + highLength) lowLength := by
  let required := m + highLength
  let fullLength := max required lowLength
  have hzero : ∃ heap1,
      (if lowLength < required then
        Generated.StrictMul.mulZeroPadLoop out lowLength
          (required - lowLength) 0 heap
       else (Except.ok heap : RawExec RawHeap)) = Except.ok heap1 ∧
      RawHeap.SameLayout heap heap1 := by
    by_cases hpad : lowLength < required
    · have hsum : lowLength + (required - lowLength) = required := by omega
      have hfull : fullLength = required := by
        simp [fullLength, Nat.max_eq_left (Nat.le_of_lt hpad)]
      have hvalidFull : heap.ValidU64Slice out fullLength := by
        simpa [fullLength, required] using hwork.outValid
      have hvalid : heap.ValidU64Slice out
          (lowLength + (required - lowLength)) := by
        simpa [hsum, hfull] using hvalidFull
      rcases mulZeroPadLoop_refines out lowLength (required - lowLength) 0
          this._p.toNat heap lowPoly this._p (Nat.zero_le _) hp hvalid
          (by simpa using hLow.2.2.1) (by simpa using hLow.2.1) with
        ⟨heap1, hrun, hlayout, _, _⟩
      exact ⟨heap1, by simpa [hpad] using hrun, hlayout⟩
    · exact ⟨heap, by simp [hpad], fun _ _ => Iff.rfl⟩
  rcases hzero with ⟨heap1, hzeroRun, hlayout1⟩
  let oldHighLength := if m ≤ lowLength then lowLength - m else 0
  have hOut1 : heap1.ValidU64Slice out fullLength :=
    (hlayout1 out fullLength).mp (by simpa [fullLength, required] using
      hwork.outValid)
  have hOutHigh1 : heap1.ValidU64Slice (out.add m) oldHighLength := by
    apply heap1.validU64Slice_add out fullLength m oldHighLength hOut1
    dsimp [oldHighLength, fullLength, required]
    split <;> omega
  have hHigh1 : heap1.ValidU64Slice high highLength :=
    (hlayout1 high highLength).mp hwork.highValid
  have hAddOut1 : heap1.ValidU64Slice (out.add m)
      (max oldHighLength highLength) := by
    apply heap1.validU64Slice_add out fullLength m
      (max oldHighLength highLength) hOut1
    dsimp [oldHighLength, fullLength, required]
    split <;> omega
  rcases polyAdd_ok this (out.add m) (out.add m) oldHighLength high
      highLength heap1 hAddOut1 hOutHigh1 hHigh1 with
    ⟨heap2, ignoredLength, hadd, hlayout2, _⟩
  have hOut2 : heap2.ValidU64Slice out fullLength :=
    (hlayout2 out fullLength).mp hOut1
  rcases normaliseU64_ok heap2 out fullLength hOut2 with
    ⟨length, hnorm, hlength⟩
  refine ⟨HgcdLiftHighResult.mk heap2 length, ?_, ?_, ?_⟩
  · simp [hgcdRecursiveLiftHigh, required, fullLength, oldHighLength,
      hzeroRun, hadd, hnorm]
  · exact fun ptr count =>
      (hlayout1 ptr count).trans (hlayout2 ptr count)
  · simpa [fullLength, required] using hlength

/-- Full semantic refinement of the generated shifted-high block.  The
temporary padded whole is treated only as a slice representation; the result
becomes normalized solely through the source's final `_poly_normalise`. -/
theorem hgcdRecursiveLiftHigh_refines (this : DenseUPolyZp)
    (out high : RawPtr UInt64) (lowLength m highLength : Nat)
    (heap : RawHeap)
    (lowPoly highPoly : Polynomial (ZMod this._p.toNat))
    (hp : this._p ≠ 0)
    (hwork : HgcdLiftHighWorkspace heap out high lowLength m highLength)
    (hLow : RawDensePolyRep this heap out lowLength lowPoly)
    (hHigh : RawDensePolyRep this heap high highLength highPoly) :
    ∃ result,
      hgcdRecursiveLiftHigh this out high lowLength m highLength heap =
        .ok result ∧
      RawHeap.SameLayout heap result.heap ∧
      RawDensePolyRep this result.heap out result.length
        (lowPoly + Polynomial.X ^ m * highPoly) := by
  let fullLength := max (m + highLength) lowLength
  let oldHighLength := if m ≤ lowLength then lowLength - m else 0
  rcases hgcdLiftHigh_zero_refines this out lowLength m highLength heap
      lowPoly hp hwork.outValid hLow with
    ⟨heap1, hzero, hlayout1, hLow1, hFullSlice1, hFullCanonical1⟩
  have hHigh1 : RawDensePolyRep this heap1 high highLength highPoly := by
    by_cases hpad : lowLength < m + highLength
    · have hrun : Generated.StrictMul.mulZeroPadLoop out lowLength
          (m + highLength - lowLength) 0 heap = .ok heap1 := by
        simpa [hpad] using hzero
      have hvalid : heap.ValidU64Slice out
          (lowLength + (m + highLength - lowLength)) := by
        have hmax : max (m + highLength) lowLength = m + highLength :=
          Nat.max_eq_left (Nat.le_of_lt hpad)
        simpa [hmax, Nat.add_sub_of_le (Nat.le_of_lt hpad)] using
          hwork.outValid
      have hsame := mulZeroPadLoop_preserves_prefix out high lowLength
        (m + highLength - lowLength) 0 highLength heap heap1 hvalid
        (by
          simpa [Nat.add_sub_of_le (Nat.le_of_lt hpad)] using
            u64SlicesDisjoint_mono hwork.zeroHighDisjoint
              (by omega) (Nat.le_refl _)) hrun
      exact rawDensePolyRep_of_same_prefix this heap heap1 high highLength
        highPoly hlayout1 hsame hHigh
    · have heq : heap1 = heap := by
        have : (Except.ok heap : RawExec RawHeap) = Except.ok heap1 := by
          simpa [hpad] using hzero
        exact (Except.ok.inj this).symm
      subst heap1
      exact hHigh
  have hparts : ∃ lowPart existingHigh : Polynomial (ZMod this._p.toNat),
      SlicePolyRep heap1 out m this._p.toNat lowPart ∧
      RawDensePolyRep this heap1 (out.add m) oldHighLength existingHigh ∧
      lowPoly = lowPart + Polynomial.X ^ m * existingHigh := by
    by_cases hm : m ≤ lowLength
    · simpa [oldHighLength, hm] using
        rawDensePolyRep_split_suffix this heap1 out lowLength m lowPoly hm hLow1
    · have hm' : lowLength < m := by omega
      have hmFull : m ≤ fullLength := by
        dsimp [fullLength]
        omega
      have hLowPrefix : SlicePolyRep heap1 out m this._p.toNat lowPoly :=
        slicePolyRep_prefix_of_coeff_zero heap1 out fullLength m
          this._p.toNat lowPoly
          (by simpa [fullLength] using
            (hlayout1 out fullLength).mp (by simpa [fullLength] using
              hwork.outValid)) hmFull
          (by simpa [fullLength] using hFullSlice1) (by
            intro degree hdegree
            exact slicePolyRep_coeff_zero_of_length_le heap1 out lowLength
              this._p.toNat lowPoly hLow1.2.2.1 degree (by omega))
      have hzeroValid : heap1.ValidU64Slice (out.add m) 0 := by
        exact heap1.validU64Slice_add out fullLength m 0
          ((hlayout1 out fullLength).mp (by simpa [fullLength] using
            hwork.outValid)) (by omega)
      refine ⟨lowPoly, 0, hLowPrefix, ?_, by simp⟩
      simpa [oldHighLength, hm] using
        rawDensePolyRep_zero_length this heap1 (out.add m) hzeroValid
  rcases hparts with ⟨lowPart, existingHigh, hLowPart, hExisting, hsplit⟩
  have hOut1 : heap1.ValidU64Slice out fullLength :=
    (hlayout1 out fullLength).mp (by simpa [fullLength] using hwork.outValid)
  have hAddOut1 : heap1.ValidU64Slice (out.add m)
      (max oldHighLength highLength) := by
    apply heap1.validU64Slice_add out fullLength m
      (max oldHighLength highLength) hOut1
    dsimp [fullLength, oldHighLength]
    split <;> omega
  rcases polyAdd_ok this (out.add m) (out.add m) oldHighLength high
      highLength heap1 hAddOut1 hExisting.1 hHigh1.1 with
    ⟨heap2, addLength, hadd, hlayout2, _⟩
  have hAdded := polyAdd_refines this (out.add m) (out.add m)
    oldHighLength high highLength heap1 heap2 addLength existingHigh highPoly
    hp hAddOut1 hExisting hHigh1 (Or.inl rfl) hwork.addAliasHigh hadd
  have hAddNorm := polyAdd_result_normalise this (out.add m) (out.add m)
    oldHighLength high highLength heap1 heap2 addLength hadd
  have hAddOut2 : heap2.ValidU64Slice (out.add m)
      (max oldHighLength highLength) :=
    (hlayout2 (out.add m) _).mp hAddOut1
  have hAddedFull := rawDensePolyRep_extend_to_normalise_input this heap2
    (out.add m) (max oldHighLength highLength) addLength
    (existingHigh + highPoly) hp hAddOut2 hAdded hAddNorm
  have hLowSame : SameU64Prefix heap1 heap2 out m :=
    polyAdd_preserves_prefix_disjoint this (out.add m) (out.add m)
      oldHighLength high out highLength m heap1 heap2 addLength hAddOut1
      hExisting.1 hHigh1.1 (by
        intro writeIndex _ readIndex _
        right
        have hwidth : RawLimbWidth.width UInt64 = 1 := rfl
        simp [RawPtr.add, hwidth]
        omega) hadd
  have hLowPart2 : SlicePolyRep heap2 out m this._p.toNat lowPart :=
    slicePolyRep_of_same_prefix heap1 heap2 out m this._p.toNat lowPart
      (heap1.validU64Slice_mono out fullLength m hOut1 (by omega))
      ((hlayout2 out m).mp
        (heap1.validU64Slice_mono out fullLength m hOut1 (by omega)))
      hLowSame hLowPart
  have hCanonicalLow1 := canonicalU64Prefix_mono heap1 out fullLength m
    this._p (by omega) (by simpa [fullLength] using hFullCanonical1)
  have hCanonicalLow2 := canonicalU64Prefix_of_same_prefix heap1 heap2 out
    m this._p (heap1.validU64Slice_mono out fullLength m hOut1 (by omega))
    hLowSame hCanonicalLow1
  have hlengthEq : m + max oldHighLength highLength = fullLength := by
    dsimp [oldHighLength, fullLength]
    split <;> omega
  have hOut2 : heap2.ValidU64Slice out fullLength :=
    (hlayout2 out fullLength).mp hOut1
  have hWholeSlice : SlicePolyRep heap2 out fullLength this._p.toNat
      (lowPoly + Polynomial.X ^ m * highPoly) := by
    have hjoin := slicePolyRep_join heap2 out m
      (max oldHighLength highLength) this._p.toNat lowPart
      (existingHigh + highPoly) (by simpa [hlengthEq] using hOut2)
      hLowPart2 hAddedFull.1
    have hpoly : lowPart + Polynomial.X ^ m *
          (existingHigh + highPoly) =
        lowPoly + Polynomial.X ^ m * highPoly := by
      rw [hsplit]
      ring
    simpa [hlengthEq, hpoly] using hjoin
  have hWholeCanonical : CanonicalU64Prefix heap2 out fullLength this._p := by
    have hjoin := canonicalU64Prefix_join heap2 out m
      (max oldHighLength highLength) this._p
      (by simpa [hlengthEq] using hOut2) hCanonicalLow2 hAddedFull.2
    simpa [hlengthEq] using hjoin
  rcases normaliseU64_ok heap2 out fullLength hOut2 with
    ⟨length, hnorm, hlength⟩
  let result := HgcdLiftHighResult.mk heap2 length
  refine ⟨result, ?_, ?_, ?_⟩
  · simp [result, hgcdRecursiveLiftHigh, fullLength, oldHighLength,
      hzero, hadd, hnorm]
  · exact fun ptr count =>
      (hlayout1 ptr count).trans (hlayout2 ptr count)
  · refine ⟨heap2.validU64Slice_mono out fullLength length hOut2 hlength,
      canonicalU64Prefix_mono heap2 out fullLength length this._p hlength
        hWholeCanonical,
      slicePolyRep_of_normaliseU64 heap2 out fullLength this._p.toNat length
        (lowPoly + Polynomial.X ^ m * highPoly) hOut2 hWholeSlice hnorm,
      normaliseU64_result_fixed heap2 out fullLength length hOut2 hnorm⟩

/-- Physical capacities needed by the optional four-entry early-return
matrix copy.  Source lengths come from `R`, exactly as in C++. -/
structure HgcdEarlyMatrixWorkspace (heap : RawHeap) (M R : HgcdMat)
    (hM : M.Valid) (hR : R.Valid) : Prop where
  targetValid : ∀ i : Fin 4, heap.ValidU64Slice
    (hgcdMatPtrRaw M hM i) (hgcdMatLenRaw R hR i)
  sourceValid : ∀ i : Fin 4, heap.ValidU64Slice
    (hgcdMatPtrRaw R hR i) (hgcdMatLenRaw R hR i)

/-- Non-aliasing obligations that make the four sequential C++ matrix copies
semantic, rather than merely executable.  Later destinations cannot overwrite
an earlier result or any still-live source entry. -/
structure HgcdEarlyMatrixRefineWorkspace (heap : RawHeap) (M R : HgcdMat)
    (hM : M.Valid) (hR : R.Valid) extends
    HgcdEarlyMatrixWorkspace heap M R hM hR where
  targetSourceDisjoint : ∀ i j : Fin 4, U64SlicesDisjoint
    (hgcdMatPtrRaw M hM i) (hgcdMatLenRaw R hR i)
    (hgcdMatPtrRaw R hR j) (hgcdMatLenRaw R hR j)
  targetPairwiseDisjoint : ∀ i j : Fin 4, i ≠ j → U64SlicesDisjoint
    (hgcdMatPtrRaw M hM i) (hgcdMatLenRaw R hR i)
    (hgcdMatPtrRaw M hM j) (hgcdMatLenRaw R hR j)

theorem hgcdEarlyMatrixRefineWorkspace_of_sameLayout
    (before after : RawHeap) (M R : HgcdMat)
    (hM : M.Valid) (hR : R.Valid)
    (hlayout : RawHeap.SameLayout before after)
    (hwork : HgcdEarlyMatrixRefineWorkspace before M R hM hR) :
    HgcdEarlyMatrixRefineWorkspace after M R hM hR := by
  exact {
    targetValid := fun i => (hlayout _ _).mp (hwork.targetValid i)
    sourceValid := fun i => (hlayout _ _).mp (hwork.sourceValid i)
    targetSourceDisjoint := hwork.targetSourceDisjoint
    targetPairwiseDisjoint := hwork.targetPairwiseDisjoint }

theorem hgcdEarlyMatrixLoop_terminates (M R : HgcdMat)
    (hM : M.Valid) (hR : R.Valid) (i : Nat) (heap : RawHeap)
    (hwork : HgcdEarlyMatrixWorkspace heap M R hM hR) :
    ∃ result, hgcdEarlyMatrixLoop M R hM hR i heap = .ok result ∧
      RawHeap.SameLayout heap result.heap ∧ result.matrix.Valid := by
  rw [hgcdEarlyMatrixLoop]
  split
  next hi =>
    let index : Fin 4 := ⟨i, hi⟩
    rcases copyU64_ok heap (hgcdMatPtrRaw M hM index)
        (hgcdMatPtrRaw R hR index) (hgcdMatLenRaw R hR index)
        (hwork.targetValid index) (hwork.sourceValid index) with
      ⟨heap1, hcopy, hlayout1⟩
    simp only [hcopy]
    let nextLen := M.len.set i (hgcdMatLenRaw R hR index)
      (by rw [hM.2]; exact hi)
    let next : HgcdMat := { M with len := nextLen }
    have hNext : next.Valid := by
      exact ⟨hM.1, by simp [next, nextLen, hM.2]⟩
    have hwork1 : HgcdEarlyMatrixWorkspace heap1 next R hNext hR := by
      constructor
      · intro j
        apply (hlayout1 (hgcdMatPtrRaw next hNext j)
          (hgcdMatLenRaw R hR j)).mp
        simpa [hgcdMatPtrRaw, next] using hwork.targetValid j
      · intro j
        exact (hlayout1 (hgcdMatPtrRaw R hR j)
          (hgcdMatLenRaw R hR j)).mp (hwork.sourceValid j)
    rcases hgcdEarlyMatrixLoop_terminates next R hNext hR (i + 1)
        heap1 hwork1 with ⟨result, hrun, hlayout2, hvalid⟩
    refine ⟨result, ?_, fun ptr count =>
      (hlayout1 ptr count).trans (hlayout2 ptr count), hvalid⟩
    simpa [index, next, nextLen, hcopy] using hrun
  next hi =>
    exact ⟨HgcdEarlyMatrixResult.mk heap M, rfl, fun _ _ => Iff.rfl, hM⟩
termination_by 4 - i
decreasing_by omega

/-- Content semantics of the actual four-copy early-return loop.  The
induction follows the generated loop and records exactly which destination
slots have already received the corresponding normalized source polynomial. -/
theorem hgcdEarlyMatrixLoop_copies (this : DenseUPolyZp)
    (M R : HgcdMat) (hM : M.Valid) (hR : R.Valid)
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (i : Nat) (heap : RawHeap) (result : HgcdEarlyMatrixResult)
    (hwork : HgcdEarlyMatrixRefineWorkspace heap M R hM hR)
    (hsource : ∀ j : Fin 4, RawDensePolyRep this heap
      (hgcdMatPtrRaw R hR j) (hgcdMatLenRaw R hR j) (entries j))
    (hdone : ∀ j : Fin 4, j.val < i → RawDensePolyRep this heap
      (hgcdMatPtrRaw M hM j) (hgcdMatLenRaw R hR j) (entries j))
    (hrun : hgcdEarlyMatrixLoop M R hM hR i heap = .ok result) :
    (∀ j : Fin 4, RawDensePolyRep this result.heap
      (hgcdMatPtrRaw R hR j) (hgcdMatLenRaw R hR j) (entries j)) ∧
    ∀ j : Fin 4, RawDensePolyRep this result.heap
      (hgcdMatPtrRaw result.matrix
        (hgcdEarlyMatrixLoop_result_valid M R hM hR i heap result hrun) j)
      (hgcdMatLenRaw R hR j) (entries j) := by
  rw [hgcdEarlyMatrixLoop] at hrun
  split at hrun
  next hi =>
    dsimp only at hrun
    split at hrun
    next fault hcopy => simp at hrun
    next heap1 hcopy =>
      let index : Fin 4 := ⟨i, hi⟩
      let dst := hgcdMatPtrRaw M hM index
      let src := hgcdMatPtrRaw R hR index
      let length := hgcdMatLenRaw R hR index
      rcases copyU64_refines_rawDense this heap dst src length (entries index)
          (hwork.targetValid index) (hwork.targetSourceDisjoint index index)
          (hsource index) with ⟨semanticHeap, hsemantic, hlayout, hcopied⟩
      have heq : semanticHeap = heap1 :=
        Except.ok.inj (hsemantic.symm.trans hcopy)
      subst semanticHeap
      have hsource1 : ∀ j : Fin 4, RawDensePolyRep this heap1
          (hgcdMatPtrRaw R hR j) (hgcdMatLenRaw R hR j) (entries j) := by
        intro j
        exact (copyU64_preserves_rawDenseRep this heap heap1 dst src length
          (hgcdMatPtrRaw R hR j) (hgcdMatLenRaw R hR j) (entries j)
          (hwork.targetValid index) (hwork.sourceValid index)
          (hwork.targetSourceDisjoint index j) hcopy (hsource j)).2
      let nextLen := M.len.set i (hgcdMatLenRaw R hR index)
        (by rw [hM.2]; exact hi)
      let next : HgcdMat := { M with len := nextLen }
      have hNext : next.Valid := by
        exact ⟨hM.1, by simp [next, nextLen, hM.2]⟩
      have hwork1 : HgcdEarlyMatrixRefineWorkspace heap1 next R hNext hR := by
        have transported := hgcdEarlyMatrixRefineWorkspace_of_sameLayout
          heap heap1 M R hM hR hlayout hwork
        exact {
          targetValid := fun j => by
            simpa [next, hgcdMatPtrRaw] using transported.targetValid j
          sourceValid := transported.sourceValid
          targetSourceDisjoint := fun j k => by
            simpa [next, hgcdMatPtrRaw] using
              transported.targetSourceDisjoint j k
          targetPairwiseDisjoint := fun j k hjk => by
            simpa [next, hgcdMatPtrRaw] using
              transported.targetPairwiseDisjoint j k hjk }
      have hdone1 : ∀ j : Fin 4, j.val < i + 1 →
          RawDensePolyRep this heap1 (hgcdMatPtrRaw next hNext j)
            (hgcdMatLenRaw R hR j) (entries j) := by
        intro j hj
        by_cases hji : j = index
        · subst j
          simpa [dst, next, hgcdMatPtrRaw] using hcopied
        · have hjlt : j.val < i := by
            have : j.val ≠ i := by
              intro hval
              exact hji (Fin.ext hval)
            omega
          have hpreserved := (copyU64_preserves_rawDenseRep this heap heap1
            dst src length (hgcdMatPtrRaw M hM j) (hgcdMatLenRaw R hR j)
            (entries j) (hwork.targetValid index) (hwork.sourceValid index)
            (hwork.targetPairwiseDisjoint index j (Ne.symm hji)) hcopy
            (hdone j hjlt)).2
          simpa [next, hgcdMatPtrRaw] using hpreserved
      exact hgcdEarlyMatrixLoop_copies this next R hNext hR entries (i + 1)
        heap1 result hwork1 hsource1 hdone1 hrun
  next hi =>
    have heq : result = HgcdEarlyMatrixResult.mk heap M :=
      (Except.ok.inj hrun).symm
    subst result
    exact ⟨hsource, fun j => hdone j (by omega)⟩
termination_by 4 - i
decreasing_by omega

/-- Entry-point refinement of the optional matrix branch: the returned
descriptor has `R`'s four lengths and its original destination pointers now
represent exactly `R`'s four L2 polynomial entries. -/
theorem hgcdEarlyMatrixLoop_zero_refines (this : DenseUPolyZp)
    (M R : HgcdMat) (hM : M.Valid) (hR : R.Valid)
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (heap : RawHeap) (result : HgcdEarlyMatrixResult)
    (hwork : HgcdEarlyMatrixRefineWorkspace heap M R hM hR)
    (hsource : HgcdMatRawDenseRep this heap R entries hR)
    (hrun : hgcdEarlyMatrixLoop M R hM hR 0 heap = .ok result) :
    ∃ hResult : result.matrix.Valid,
      result.matrix.len = R.len ∧
      HgcdMatRawDenseRep this result.heap result.matrix entries hResult := by
  have hvalid := hgcdEarlyMatrixLoop_result_valid M R hM hR 0 heap result hrun
  have hlen := hgcdEarlyMatrixLoop_lengths M R hM hR 0 heap result
    (by intro j hj; omega) hrun
  have hcopies := hgcdEarlyMatrixLoop_copies this M R hM hR entries 0 heap
    result hwork (by
      intro j
      simpa [HgcdMatRawDenseRep, hgcdMatPtr, hgcdMatLen,
        hgcdMatPtrRaw, hgcdMatLenRaw] using hsource j)
    (by intro j hj; omega) hrun
  refine ⟨hvalid, hlen, ?_⟩
  intro j
  have hLength : hgcdMatLen result.matrix hvalid j =
      hgcdMatLenRaw R hR j := by
    exact array_getElem_eq_of_eq result.matrix.len R.len hlen j.val
      (by rw [hvalid.2]; exact j.isLt)
      (by rw [hR.2]; exact j.isLt)
  rw [hLength]
  simpa [hgcdMatPtr, hgcdMatPtrRaw] using hcopies.2 j

/-- A normalized polynomial outside all four matrix destinations is framed by
the complete generated early-copy loop. -/
theorem hgcdEarlyMatrixLoop_preserves_rawDenseRep (this : DenseUPolyZp)
    (M R : HgcdMat) (hM : M.Valid) (hR : R.Valid)
    (i : Nat) (heap : RawHeap) (result : HgcdEarlyMatrixResult)
    (ptr : RawPtr UInt64) (length : Nat)
    (poly : Polynomial (ZMod this._p.toNat))
    (hwork : HgcdEarlyMatrixWorkspace heap M R hM hR)
    (hframe : ∀ j : Fin 4, U64SlicesDisjoint
      (hgcdMatPtrRaw M hM j) (hgcdMatLenRaw R hR j) ptr length)
    (hrep : RawDensePolyRep this heap ptr length poly)
    (hrun : hgcdEarlyMatrixLoop M R hM hR i heap = .ok result) :
    RawDensePolyRep this result.heap ptr length poly := by
  rw [hgcdEarlyMatrixLoop] at hrun
  split at hrun
  next hi =>
    dsimp only at hrun
    split at hrun
    next fault hcopy => simp at hrun
    next heap1 hcopy =>
      let index : Fin 4 := ⟨i, hi⟩
      have hcopyFrame := copyU64_preserves_rawDenseRep this heap heap1
        (hgcdMatPtrRaw M hM index) (hgcdMatPtrRaw R hR index)
        (hgcdMatLenRaw R hR index) ptr length poly
        (hwork.targetValid index) (hwork.sourceValid index) (hframe index)
        hcopy hrep
      have hpreserved := hcopyFrame.2
      let nextLen := M.len.set i (hgcdMatLenRaw R hR index)
        (by rw [hM.2]; exact hi)
      let next : HgcdMat := { M with len := nextLen }
      have hNext : next.Valid := by
        exact ⟨hM.1, by simp [next, nextLen, hM.2]⟩
      have hwork1 : HgcdEarlyMatrixWorkspace heap1 next R hNext hR := by
        exact {
          targetValid := fun j => by
            apply (hcopyFrame.1 _ _).mp
            simpa [next, hgcdMatPtrRaw] using hwork.targetValid j
          sourceValid := fun j => by
            apply (hcopyFrame.1 _ _).mp
            exact hwork.sourceValid j }
      have hframe1 : ∀ j : Fin 4, U64SlicesDisjoint
          (hgcdMatPtrRaw next hNext j) (hgcdMatLenRaw R hR j) ptr length := by
        intro j
        simpa [next, hgcdMatPtrRaw] using hframe j
      exact hgcdEarlyMatrixLoop_preserves_rawDenseRep this next R hNext hR
        (i + 1) heap1 result ptr length poly hwork1 hframe1 hpreserved hrun
  next hi =>
    have heq : result = HgcdEarlyMatrixResult.mk heap M :=
      (Except.ok.inj hrun).symm
    subst result
    exact hrep
termination_by 4 - i
decreasing_by omega

structure HgcdEarlyReturnWorkspace (heap : RawHeap)
    (M R : HgcdMat) (hM : M.Valid) (hR : R.Valid)
    (A B a2 b2 : RawPtr UInt64) (lenA2 lenB2 : Nat) : Prop where
  AValid : heap.ValidU64Slice A lenA2
  BValid : heap.ValidU64Slice B lenB2
  a2Valid : heap.ValidU64Slice a2 lenA2
  b2Valid : heap.ValidU64Slice b2 lenB2
  matrix : HgcdEarlyMatrixWorkspace heap M R hM hR

/-- Full separation contract for the real `_hgcd_recursive` early return.
It mentions only physical slices; expected L2 output values remain theorem
inputs and are obtained from the actual memcpy executions. -/
structure HgcdEarlyReturnRefineWorkspace (heap : RawHeap)
    (M R : HgcdMat) (hM : M.Valid) (hR : R.Valid)
    (A B a2 b2 : RawPtr UInt64) (lenA2 lenB2 : Nat) : Prop where
  AValid : heap.ValidU64Slice A lenA2
  BValid : heap.ValidU64Slice B lenB2
  matrix : HgcdEarlyMatrixRefineWorkspace heap M R hM hR
  Aa2Disjoint : U64SlicesDisjoint A lenA2 a2 lenA2
  Ab2Disjoint : U64SlicesDisjoint A lenA2 b2 lenB2
  Bb2Disjoint : U64SlicesDisjoint B lenB2 b2 lenB2
  BADisjoint : U64SlicesDisjoint B lenB2 A lenA2
  ARDisjoint : ∀ j : Fin 4, U64SlicesDisjoint A lenA2
    (hgcdMatPtrRaw R hR j) (hgcdMatLenRaw R hR j)
  BRDisjoint : ∀ j : Fin 4, U64SlicesDisjoint B lenB2
    (hgcdMatPtrRaw R hR j) (hgcdMatLenRaw R hR j)
  matrixADisjoint : ∀ j : Fin 4, U64SlicesDisjoint
    (hgcdMatPtrRaw M hM j) (hgcdMatLenRaw R hR j) A lenA2
  matrixBDisjoint : ∀ j : Fin 4, U64SlicesDisjoint
    (hgcdMatPtrRaw M hM j) (hgcdMatLenRaw R hR j) B lenB2

/-- Total execution bridge for the generated early-return branch, including
the optional matrix copy. -/
theorem hgcdRecursiveEarlyReturn_terminates (M R : HgcdMat)
    (hM : M.Valid) (hR : R.Valid) (computeM : Bool)
    (A B a2 b2 : RawPtr UInt64) (lenA2 lenB2 : Nat) (sgn : Int)
    (heap : RawHeap)
    (hwork : HgcdEarlyReturnWorkspace heap M R hM hR A B a2 b2
      lenA2 lenB2) :
    ∃ result,
      hgcdRecursiveEarlyReturn M R hM hR computeM A B a2 b2
        lenA2 lenB2 sgn heap = .ok result ∧
      RawHeap.SameLayout heap result.heap ∧
      result.matrix.Valid ∧ result.lenA = lenA2 ∧
      result.lenB = lenB2 ∧ result.sgn = sgn := by
  rcases copyU64_ok heap A a2 lenA2 hwork.AValid hwork.a2Valid with
    ⟨heap1, hcopyA, hlayout1⟩
  have hB1 := (hlayout1 B lenB2).mp hwork.BValid
  have hb21 := (hlayout1 b2 lenB2).mp hwork.b2Valid
  rcases copyU64_ok heap1 B b2 lenB2 hB1 hb21 with
    ⟨heap2, hcopyB, hlayout2⟩
  by_cases hcompute : computeM = true
  · have hmatrix2 : HgcdEarlyMatrixWorkspace heap2 M R hM hR := by
      constructor <;> intro j
      · exact (hlayout2 _ _).mp ((hlayout1 _ _).mp
          (hwork.matrix.targetValid j))
      · exact (hlayout2 _ _).mp ((hlayout1 _ _).mp
          (hwork.matrix.sourceValid j))
    rcases hgcdEarlyMatrixLoop_terminates M R hM hR 0 heap2 hmatrix2 with
      ⟨matrixResult, hmatrix, hlayout3, hvalid⟩
    refine ⟨HgcdRecursiveEarlyResult.mk matrixResult.heap
        matrixResult.matrix lenA2 lenB2 sgn, ?_, ?_, hvalid,
      rfl, rfl, rfl⟩
    · simp [hgcdRecursiveEarlyReturn, hcopyA, hcopyB, hcompute, hmatrix]
    · exact fun ptr count =>
        (hlayout1 ptr count).trans
          ((hlayout2 ptr count).trans (hlayout3 ptr count))
  · have hfalse : computeM = false := by cases computeM <;> simp_all
    refine ⟨HgcdRecursiveEarlyResult.mk heap2 M lenA2 lenB2 sgn,
      ?_, ?_, hM, rfl, rfl, rfl⟩
    · simp [hgcdRecursiveEarlyReturn, hcopyA, hcopyB, hfalse]
    · exact fun ptr count =>
        (hlayout1 ptr count).trans (hlayout2 ptr count)

/-- End-to-end semantic refinement of the generated early-return block.  It
executes both output memcpys and, when requested, the four matrix memcpys;
the returned A/B polynomials and matrix are therefore consequences of those
same RawHeap executions. -/
theorem hgcdRecursiveEarlyReturn_refines (this : DenseUPolyZp)
    (M R : HgcdMat) (hM : M.Valid) (hR : R.Valid) (computeM : Bool)
    (A B a2 b2 : RawPtr UInt64) (lenA2 lenB2 : Nat) (sgn : Int)
    (left right : Polynomial (ZMod this._p.toNat))
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (heap : RawHeap)
    (hwork : HgcdEarlyReturnRefineWorkspace heap M R hM hR
      A B a2 b2 lenA2 lenB2)
    (hLeft : RawDensePolyRep this heap a2 lenA2 left)
    (hRight : RawDensePolyRep this heap b2 lenB2 right)
    (hMatrix : HgcdMatRawDenseRep this heap R entries hR) :
    ∃ result,
      hgcdRecursiveEarlyReturn M R hM hR computeM A B a2 b2
        lenA2 lenB2 sgn heap = .ok result ∧
      RawHeap.SameLayout heap result.heap ∧
      result.lenA = lenA2 ∧ result.lenB = lenB2 ∧ result.sgn = sgn ∧
      RawDensePolyRep this result.heap A lenA2 left ∧
      RawDensePolyRep this result.heap B lenB2 right ∧
      ∃ hResult : result.matrix.Valid,
        (computeM = true →
          result.matrix.len = R.len ∧
          HgcdMatRawDenseRep this result.heap result.matrix entries hResult) ∧
        (computeM = false → result.matrix = M) := by
  rcases copyU64_refines_rawDense this heap A a2 lenA2 left hwork.AValid
      hwork.Aa2Disjoint hLeft with ⟨heap1, hcopyA, hlayout1, hA1⟩
  have hRight1 := (copyU64_preserves_rawDenseRep this heap heap1 A a2 lenA2
    b2 lenB2 right hwork.AValid hLeft.1 hwork.Ab2Disjoint hcopyA hRight).2
  have hMatrix1 : HgcdMatRawDenseRep this heap1 R entries hR := by
    intro j
    exact (copyU64_preserves_rawDenseRep this heap heap1 A a2 lenA2
      (hgcdMatPtrRaw R hR j) (hgcdMatLenRaw R hR j) (entries j)
      hwork.AValid hLeft.1 (hwork.ARDisjoint j) hcopyA (hMatrix j)).2
  have hB1 := (hlayout1 B lenB2).mp hwork.BValid
  rcases copyU64_refines_rawDense this heap1 B b2 lenB2 right hB1
      hwork.Bb2Disjoint hRight1 with ⟨heap2, hcopyB, hlayout2, hB2⟩
  have hA2 := (copyU64_preserves_rawDenseRep this heap1 heap2 B b2 lenB2
    A lenA2 left hB1 hRight1.1 hwork.BADisjoint hcopyB hA1).2
  have hMatrix2 : HgcdMatRawDenseRep this heap2 R entries hR := by
    intro j
    exact (copyU64_preserves_rawDenseRep this heap1 heap2 B b2 lenB2
      (hgcdMatPtrRaw R hR j) (hgcdMatLenRaw R hR j) (entries j)
      hB1 hRight1.1 (hwork.BRDisjoint j) hcopyB (hMatrix1 j)).2
  have hMatrixWork1 := hgcdEarlyMatrixRefineWorkspace_of_sameLayout heap
    heap1 M R hM hR hlayout1 hwork.matrix
  have hMatrixWork2 := hgcdEarlyMatrixRefineWorkspace_of_sameLayout heap1
    heap2 M R hM hR hlayout2 hMatrixWork1
  by_cases hcompute : computeM = true
  · rcases hgcdEarlyMatrixLoop_terminates M R hM hR 0 heap2
        hMatrixWork2.toHgcdEarlyMatrixWorkspace with
      ⟨matrixResult, hmatrixRun, hlayout3, hmatrixValid⟩
    rcases hgcdEarlyMatrixLoop_zero_refines this M R hM hR entries heap2
        matrixResult hMatrixWork2 hMatrix2 hmatrixRun with
      ⟨hResult, hResultLen, hResultRep⟩
    have hAResult := hgcdEarlyMatrixLoop_preserves_rawDenseRep this M R hM hR
      0 heap2 matrixResult A lenA2 left
      hMatrixWork2.toHgcdEarlyMatrixWorkspace hwork.matrixADisjoint hA2
      hmatrixRun
    have hBResult := hgcdEarlyMatrixLoop_preserves_rawDenseRep this M R hM hR
      0 heap2 matrixResult B lenB2 right
      hMatrixWork2.toHgcdEarlyMatrixWorkspace hwork.matrixBDisjoint hB2
      hmatrixRun
    let result : HgcdRecursiveEarlyResult := {
      heap := matrixResult.heap
      matrix := matrixResult.matrix
      lenA := lenA2
      lenB := lenB2
      sgn := sgn }
    refine ⟨result, ?_, ?_, rfl, rfl, rfl, hAResult, hBResult,
      hResult, ?_, ?_⟩
    · simp [hgcdRecursiveEarlyReturn, hcopyA, hcopyB, hcompute, hmatrixRun,
        result]
    · exact fun ptr count => (hlayout1 ptr count).trans
        ((hlayout2 ptr count).trans (hlayout3 ptr count))
    · intro _
      exact ⟨hResultLen, hResultRep⟩
    · intro hfalse
      simp_all
  · have hfalse : computeM = false := by cases computeM <;> simp_all
    let result : HgcdRecursiveEarlyResult := {
      heap := heap2
      matrix := M
      lenA := lenA2
      lenB := lenB2
      sgn := sgn }
    refine ⟨result, ?_, ?_, rfl, rfl, rfl, hA2, hB2, hM, ?_, ?_⟩
    · simp [hgcdRecursiveEarlyReturn, hcopyA, hcopyB, hfalse, result]
    · exact fun ptr count =>
        (hlayout1 ptr count).trans (hlayout2 ptr count)
    · intro htrue
      simp_all
    · intro _
      rfl

/-- Strict semantic bridge for the middle source `_poly_divrem`.  Its result
record also exposes the exact `k`, `c0`, and `d0` passed to the second HGCD
call, while quotient and remainder come from the same raw execution. -/
theorem hgcdRecursiveMiddle_refines (this : DenseUPolyZp)
    (q d a2 b2 : RawPtr UInt64) (lenA2 lenB2 m lenA : Nat)
    (W3 : RawPtr Word3) (heap : RawHeap)
    (dividend divisor : Polynomial (ZMod this._p.toNat))
    (hlenB : 0 < lenB2)
    (hA : RawDensePolyRep this heap a2 lenA2 dividend)
    (hB : RawDensePolyRep this heap b2 lenB2 divisor)
    (hQ : heap.ValidU64Slice q (lenA2 - (lenB2 - 1)))
    (hD : heap.ValidU64Slice d (Nat.min lenA2 (lenB2 - 1)))
    (hW3 : heap.ValidWord3Slice W3 lenA2)
    (hqCapacity : lenA2 - (lenB2 - 1) < limbBase)
    (hDA : d.region ≠ a2.region)
    (hWA : W3.region ≠ a2.region) (hWB : W3.region ≠ b2.region)
    (hQB : q.region ≠ b2.region) (hQW : q.region ≠ W3.region)
    (hDW : d.region ≠ W3.region) (hDQ : d.region ≠ q.region)
    (hDB : d.region ≠ b2.region)
    (hcfg : DensePreinvConfigured this)
    (hprime : Nat.Prime this._p.toNat)
    (hlenA : 0 < lenA) (hlenB2Bound : lenB2 ≤ lenA)
    (hm : 0 < m) (hnonearly : m + 1 ≤ lenB2) :
    ∃ result quotient remainder,
      hgcdRecursiveMiddle this q d a2 b2 lenA2 lenB2 m W3 heap =
        .ok result ∧
      SlicePolyRep result.heap q result.lenQ this._p.toNat quotient ∧
      CanonicalU64Prefix result.heap q result.lenQ this._p ∧
      result.heap.normaliseU64 q result.lenQ = .ok result.lenQ ∧
      SlicePolyRep result.heap d result.lenD this._p.toNat remainder ∧
      CanonicalU64Prefix result.heap d result.lenD this._p ∧
      result.heap.normaliseU64 d result.lenD = .ok result.lenD ∧
      RawHeap.SameLayout heap result.heap ∧
      RawDensePolyRep this result.heap q result.lenQ quotient ∧
      RawDensePolyRep this result.heap d result.lenD remainder ∧
      RawDensePolyRep this result.heap b2 lenB2 divisor ∧
      dividend = quotient * divisor + remainder ∧
      (remainder = 0 ∨ remainder.natDegree < divisor.natDegree) ∧
      result.lenQ ≤ lenA2 - (lenB2 - 1) ∧
      result.lenD ≤ Nat.min lenA2 (lenB2 - 1) ∧
      result.lenD < lenB2 ∧ result.lenC0 < lenA ∧
      0 < result.lenC0 ∧ result.lenD0 < result.lenC0 ∧
      result.k = 2 * m - lenB2 + 1 ∧
      result.c0 = b2.add result.k ∧
      result.lenC0 = (if lenB2 ≥ result.k then lenB2 - result.k else 0) ∧
      result.d0 = d.add result.k ∧
      result.lenD0 = (if result.lenD ≥ result.k then
        result.lenD - result.k else 0) := by
  rcases polyDivrem_refines this q d a2 b2 lenA2 lenB2 W3 heap
      dividend divisor hlenB hA.1 hB.1 hQ hD hW3 hA.2.1 hB.2.1
      hA.2.2.1 hB.2.2.1 hA.2.2.2 hB.2.2.2 hqCapacity hDA hWA hWB
      hQB hQW hDW hDQ hDB hcfg hprime with
    ⟨heap1, lenQ, lenD, quotient, remainder, hdiv, hQRep, hQCanonical,
      hQNorm, hDRep, hDCanonical, hDNorm, hlayout, hsameB, hidentity,
      hdegree, hlenQ, hlenD, hlt⟩
  let result : HgcdRecursiveMiddleResult := {
    heap := heap1
    lenQ := lenQ
    lenD := lenD
    k := 2 * m - lenB2 + 1
    c0 := b2.add (2 * m - lenB2 + 1)
    lenC0 := if lenB2 ≥ 2 * m - lenB2 + 1 then
      lenB2 - (2 * m - lenB2 + 1) else 0
    d0 := d.add (2 * m - lenB2 + 1)
    lenD0 := if lenD ≥ 2 * m - lenB2 + 1 then
      lenD - (2 * m - lenB2 + 1) else 0 }
  have hmiddle : hgcdRecursiveMiddle this q d a2 b2 lenA2 lenB2 m W3
      heap = .ok result := by
    simp [hgcdRecursiveMiddle, hdiv, result]
  have hQValid : heap1.ValidU64Slice q lenQ :=
    heap1.validU64Slice_mono q (lenA2 - (lenB2 - 1)) lenQ
      ((hlayout q _).mp hQ) hlenQ
  have hDValid : heap1.ValidU64Slice d lenD :=
    heap1.validU64Slice_mono d (Nat.min lenA2 (lenB2 - 1)) lenD
      ((hlayout d _).mp hD) hlenD
  have hQDense : RawDensePolyRep this heap1 q lenQ quotient :=
    ⟨hQValid, hQCanonical, hQRep, hQNorm⟩
  have hDDense : RawDensePolyRep this heap1 d lenD remainder :=
    ⟨hDValid, hDCanonical, hDRep, hDNorm⟩
  have hBDense := rawDensePolyRep_of_same_prefix this heap heap1 b2 lenB2
    divisor hlayout hsameB hB
  rcases hgcdRecursiveMiddle_second_call_bounds this q d a2 b2 lenA2 lenB2
      m lenA W3 heap result hlenA hm hlenB2Bound hnonearly hlt hmiddle with
    ⟨hcpos, hdorder, hdecrease⟩
  refine ⟨result, quotient, remainder, hmiddle, hQRep, hQCanonical, hQNorm,
    hDRep, hDCanonical, hDNorm, hlayout, hQDense, hDDense, hBDense,
    hidentity, hdegree, hlenQ, hlenD, hlt,
    hdecrease, hcpos, hdorder, rfl, rfl, rfl, rfl, rfl⟩

/-- The nonempty second-call input is obtained by slicing the actual divisor
and remainder buffers at the generated `k`.  No heap write or L2 evaluation
is used to create `c0` or `d0`; a zero remainder suffix uses only its physical
zero-length validity. -/
theorem hgcdRecursiveMiddle_suffix_reps (this : DenseUPolyZp)
    (q d a2 b2 : RawPtr UInt64) (lenA2 lenB2 m : Nat)
    (W3 : RawPtr Word3) (heap : RawHeap)
    (result : HgcdRecursiveMiddleResult)
    (divisor remainder : Polynomial (ZMod this._p.toNat))
    (hDivisor : RawDensePolyRep this result.heap b2 lenB2 divisor)
    (hRemainder : RawDensePolyRep this result.heap d result.lenD remainder)
    (hD0Valid : result.heap.ValidU64Slice result.d0 result.lenD0)
    (hC0Pos : 0 < result.lenC0)
    (hrun : hgcdRecursiveMiddle this q d a2 b2 lenA2 lenB2 m W3 heap =
      .ok result) :
    ∃ cPoly dPoly : Polynomial (ZMod this._p.toNat),
      RawDensePolyRep this result.heap result.c0 result.lenC0 cPoly ∧
      RawDensePolyRep this result.heap result.d0 result.lenD0 dPoly := by
  have hlayout := hgcdRecursiveMiddle_layout this q d a2 b2 lenA2 lenB2 m
    W3 heap result hrun
  have hkB : result.k ≤ lenB2 := by
    rw [hlayout.2.2.2.1] at hC0Pos
    split at hC0Pos
    next hk => exact hk
    next hk => simp at hC0Pos
  rcases rawDensePolyRep_split_suffix this result.heap b2 lenB2 result.k
      divisor hkB hDivisor with ⟨cLow, cPoly, hCLow, hCRep, hCSplit⟩
  have hCResult : RawDensePolyRep this result.heap result.c0 result.lenC0
      cPoly := by
    rw [hlayout.2.2.1, hlayout.2.2.2.1]
    simp [hkB]
    exact hCRep
  by_cases hkD : result.k ≤ result.lenD
  · rcases rawDensePolyRep_split_suffix this result.heap d result.lenD
        result.k remainder hkD hRemainder with
      ⟨dLow, dPoly, hDLow, hDRep, hDSplit⟩
    refine ⟨cPoly, dPoly, hCResult, ?_⟩
    rw [hlayout.2.2.2.2.1, hlayout.2.2.2.2.2]
    simp [hkD]
    exact hDRep
  · have hLenD0 : result.lenD0 = 0 := by
      rw [hlayout.2.2.2.2.2]
      simp [hkD]
    have hZero := rawDensePolyRep_zero_length this result.heap result.d0
      (by simpa [hLenD0] using hD0Valid)
    exact ⟨cPoly, 0, hCResult, by simpa [hLenD0] using hZero⟩

/-- Complete low/high decomposition consumed by the second recursive call
and final reconstruction.  Low parts are fixed-length canonical slices,
while the nonempty suffixes inherit normalization from the original divrem
outputs. -/
theorem hgcdRecursiveMiddle_split_reps (this : DenseUPolyZp)
    (q d a2 b2 : RawPtr UInt64) (lenA2 lenB2 m : Nat)
    (W3 : RawPtr Word3) (heap : RawHeap)
    (result : HgcdRecursiveMiddleResult)
    (divisor remainder : Polynomial (ZMod this._p.toNat))
    (hDivisor : RawDensePolyRep this result.heap b2 lenB2 divisor)
    (hRemainder : RawDensePolyRep this result.heap d result.lenD remainder)
    (hD0Valid : result.heap.ValidU64Slice result.d0 result.lenD0)
    (hC0Pos : 0 < result.lenC0)
    (hrun : hgcdRecursiveMiddle this q d a2 b2 lenA2 lenB2 m W3 heap =
      .ok result) :
    ∃ lowC highC lowD highD : Polynomial (ZMod this._p.toNat),
      RawCanonicalPolySlice this result.heap b2 (Nat.min lenB2 result.k)
        lowC ∧
      RawDensePolyRep this result.heap result.c0 result.lenC0 highC ∧
      RawCanonicalPolySlice this result.heap d
        (Nat.min result.lenD result.k) lowD ∧
      RawDensePolyRep this result.heap result.d0 result.lenD0 highD ∧
      divisor = lowC + Polynomial.X ^ result.k * highC ∧
      remainder = lowD + Polynomial.X ^ result.k * highD := by
  have hlayout := hgcdRecursiveMiddle_layout this q d a2 b2 lenA2 lenB2 m
    W3 heap result hrun
  have hkB : result.k ≤ lenB2 := by
    rw [hlayout.2.2.2.1] at hC0Pos
    split at hC0Pos
    next hk => exact hk
    next hk => simp at hC0Pos
  rcases rawDensePolyRep_split_suffix this result.heap b2 lenB2 result.k
      divisor hkB hDivisor with
    ⟨lowC, highC, hLowC, hHighC, hSplitC⟩
  have hLowCCanonical : RawCanonicalPolySlice this result.heap b2
      (Nat.min lenB2 result.k) lowC := by
    simp only [Nat.min_eq_right hkB]
    exact ⟨result.heap.validU64Slice_mono b2 lenB2 result.k hDivisor.1
        hkB,
      (canonicalU64Prefix_split result.heap b2 result.k
        (lenB2 - result.k) this._p (by
          simpa [Nat.add_sub_of_le hkB] using hDivisor.2.1)).1,
      hLowC⟩
  have hHighCResult : RawDensePolyRep this result.heap result.c0
      result.lenC0 highC := by
    rw [hlayout.2.2.1, hlayout.2.2.2.1]
    simp [hkB]
    exact hHighC
  by_cases hkD : result.k ≤ result.lenD
  · rcases rawDensePolyRep_split_suffix this result.heap d result.lenD
        result.k remainder hkD hRemainder with
      ⟨lowD, highD, hLowD, hHighD, hSplitD⟩
    have hLowDCanonical : RawCanonicalPolySlice this result.heap d
        (Nat.min result.lenD result.k) lowD := by
      simp only [Nat.min_eq_right hkD]
      exact ⟨result.heap.validU64Slice_mono d result.lenD result.k
          hRemainder.1 hkD,
        (canonicalU64Prefix_split result.heap d result.k
          (result.lenD - result.k) this._p (by
            simpa [Nat.add_sub_of_le hkD] using hRemainder.2.1)).1,
        hLowD⟩
    have hHighDResult : RawDensePolyRep this result.heap result.d0
        result.lenD0 highD := by
      rw [hlayout.2.2.2.2.1, hlayout.2.2.2.2.2]
      simp [hkD]
      exact hHighD
    exact ⟨lowC, highC, lowD, highD, hLowCCanonical, hHighCResult,
      hLowDCanonical, hHighDResult, hSplitC, hSplitD⟩
  · have hLenD0 : result.lenD0 = 0 := by
      rw [hlayout.2.2.2.2.2]
      simp [hkD]
    have hLowDCanonical : RawCanonicalPolySlice this result.heap d
        (Nat.min result.lenD result.k) remainder := by
      have hmin : Nat.min result.lenD result.k = result.lenD :=
        Nat.min_eq_left (by omega)
      simpa [hmin] using
        (show RawCanonicalPolySlice this result.heap d result.lenD remainder
          from ⟨hRemainder.1, hRemainder.2.1, hRemainder.2.2.1⟩)
    have hHighDZero : RawDensePolyRep this result.heap result.d0
        result.lenD0 0 := by
      simpa [hLenD0] using rawDensePolyRep_zero_length this result.heap
        result.d0 (by simpa [hLenD0] using hD0Valid)
    have hSplitD : remainder = remainder + Polynomial.X ^ result.k * 0 := by
      simp
    exact ⟨lowC, highC, remainder, 0, hLowCCanonical, hHighCResult,
      hLowDCanonical, hHighDZero, hSplitC, hSplitD⟩

/-- A readable limb `1` is the normalized raw representation of the constant
one whenever the C++ modulus has at least two residues. -/
theorem rawDensePolyRep_one_of_read_one (this : DenseUPolyZp)
    (heap : RawHeap) (ptr : RawPtr UInt64)
    (hp : 1 < this._p.toNat)
    (hvalid : heap.ValidU64Slice ptr 1)
    (hread : heap.readU64 ptr 0 = .ok 1) :
    RawDensePolyRep this heap ptr 1 1 := by
  refine ⟨hvalid, ?_, slicePolyRep_one_of_read_one heap ptr
    this._p.toNat hvalid hread, ?_⟩
  · intro k value hk hvalue
    have hk0 : k = 0 := by omega
    subst k
    have hv : value = 1 := Except.ok.inj (hvalue.symm.trans hread)
    subst value
    simpa using hp
  · simp [RawHeap.normaliseU64, hread]

theorem matOne_refines (M : HgcdMat) (heap : RawHeap) (p : Nat)
    (hM : M.Valid)
    (h0 : heap.ValidU64Slice (hgcdMatPtr M hM ⟨0, by omega⟩) 1)
    (h3 : heap.ValidU64Slice (hgcdMatPtr M hM ⟨3, by omega⟩) 1)
    (h03 : CLPoly.Impl.StrictMulRefinement.U64SlicesDisjoint
      (hgcdMatPtr M hM ⟨0, by omega⟩) 1
      (hgcdMatPtr M hM ⟨3, by omega⟩) 1) :
    ∃ heap' M', dense_upoly_zp__mat_one_ir M heap = .ok (heap', M') ∧
      RawHeap.SameLayout heap heap' ∧
      M'.poly = M.poly ∧ M'.len = #[1, 0, 0, 1] ∧
      heap'.readU64 (hgcdMatPtr M hM ⟨0, by omega⟩) 0 = .ok 1 ∧
      heap'.readU64 (hgcdMatPtr M hM ⟨3, by omega⟩) 0 = .ok 1 ∧
      ∃ hM' : M'.Valid,
        HgcdMatPolyRep heap' M' p (identityEntries p) hM' := by
  have hvalid : M.poly.size = 4 ∧ M.len.size = 4 := by
    simpa [HgcdMat.Valid] using hM
  let p0 := hgcdMatPtr M hM ⟨0, by omega⟩
  let p3 := hgcdMatPtr M hM ⟨3, by omega⟩
  rcases heap.writeU64_of_valid p0 1 0 1 (by simpa [p0] using h0)
      (by omega) with ⟨heap1, hwrite0⟩
  have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 p0 0 1 hwrite0
  have h3one := (hlayout1 p3 1).mp (by simpa [p3] using h3)
  rcases heap1.writeU64_of_valid p3 1 0 1 h3one (by omega) with
    ⟨heap2, hwrite3⟩
  have hlayout2 := RawHeap.writeU64_sameLayout heap1 heap2 p3 0 1 hwrite3
  let M' : HgcdMat := { M with len := #[1, 0, 0, 1] }
  have hwrite0' : heap.writeU64 (M.poly[0]'(by omega)) 0 1 = .ok heap1 := by
    simpa [p0, hgcdMatPtr] using hwrite0
  have hwrite3' : heap1.writeU64 (M.poly[3]'(by omega)) 0 1 = .ok heap2 := by
    simpa [p3, hgcdMatPtr] using hwrite3
  have hrun : dense_upoly_zp__mat_one_ir M heap = .ok (heap2, M') := by
    simp [dense_upoly_zp__mat_one_ir, hvalid, hwrite0', hwrite3', M']
  have hsame0 := CLPoly.Impl.StrictMulRefinement.writeU64_preserves_prefix
    heap1 heap2 p3 p0 1 1 0 1
    (by simpa [p0, p3] using
      CLPoly.Impl.StrictMulRefinement.u64SlicesDisjoint_symm h03)
    (by omega) hwrite3
  have hread0Heap1 := RawHeap.readU64_writeU64_same heap heap1 p0 0 1 hwrite0
  have hread0 : heap2.readU64 p0 0 = .ok 1 :=
    hsame0 0 1 (by omega) hread0Heap1
  have hread3 : heap2.readU64 p3 0 = .ok 1 :=
    RawHeap.readU64_writeU64_same heap1 heap2 p3 0 1 hwrite3
  have hM' : M'.Valid := by simp [M', HgcdMat.Valid, hvalid.1]
  refine ⟨heap2, M', hrun, fun ptr length =>
    (hlayout1 ptr length).trans (hlayout2 ptr length), rfl, rfl,
    (by simpa [p0]), (by simpa [p3]), hM', ?_⟩
  intro i
  fin_cases i
  · simpa [HgcdMatPolyRep, hgcdMatPtr, hgcdMatLen, M', identityEntries,
      p0] using
      slicePolyRep_one_of_read_one heap2 p0 p
        ((hlayout2 p0 1).mp ((hlayout1 p0 1).mp (by simpa [p0] using h0)))
        hread0
  · simpa [HgcdMatPolyRep, hgcdMatPtr, hgcdMatLen, M', identityEntries] using
      slicePolyRep_zero_length_any heap2
        (hgcdMatPtr M hM ⟨1, by omega⟩) p
  · simpa [HgcdMatPolyRep, hgcdMatPtr, hgcdMatLen, M', identityEntries] using
      slicePolyRep_zero_length_any heap2
        (hgcdMatPtr M hM ⟨2, by omega⟩) p
  · simpa [HgcdMatPolyRep, hgcdMatPtr, hgcdMatLen, M', identityEntries,
      p3] using
      slicePolyRep_one_of_read_one heap2 p3 p ((hlayout2 p3 1).mp h3one)
        hread3

/-- The two real writes in `_mat_one` preserve every declared prefix outside
the two matrix-entry cells. -/
theorem matOne_preserves_prefix (M : HgcdMat) (heap heap' : RawHeap)
    (M' : HgcdMat) (guard : RawPtr UInt64) (guardLen : Nat)
    (hM : M.Valid)
    (h0Guard : U64SlicesDisjoint
      (hgcdMatPtr M hM ⟨0, by omega⟩) 1 guard guardLen)
    (h3Guard : U64SlicesDisjoint
      (hgcdMatPtr M hM ⟨3, by omega⟩) 1 guard guardLen)
    (hrun : dense_upoly_zp__mat_one_ir M heap = .ok (heap', M')) :
    SameU64Prefix heap heap' guard guardLen := by
  have hvalid : M.poly.size = 4 ∧ M.len.size = 4 := by
    simpa [HgcdMat.Valid] using hM
  let p0 := hgcdMatPtr M hM ⟨0, by omega⟩
  let p3 := hgcdMatPtr M hM ⟨3, by omega⟩
  cases hwrite0 : heap.writeU64 p0 0 1 with
  | error fault =>
      have hwrite0' : heap.writeU64 (M.poly[0]'(by omega)) 0 1 =
          .error fault := by simpa [p0, hgcdMatPtr] using hwrite0
      simp [dense_upoly_zp__mat_one_ir, hvalid, hwrite0'] at hrun
  | ok heap1 =>
      have hwrite0' : heap.writeU64 (M.poly[0]'(by omega)) 0 1 =
          .ok heap1 := by simpa [p0, hgcdMatPtr] using hwrite0
      cases hwrite3 : heap1.writeU64 p3 0 1 with
      | error fault =>
          have hwrite3' : heap1.writeU64 (M.poly[3]'(by omega)) 0 1 =
              .error fault := by simpa [p3, hgcdMatPtr] using hwrite3
          simp [dense_upoly_zp__mat_one_ir, hvalid, hwrite0', hwrite3'] at hrun
      | ok heap2 =>
          have hwrite3' : heap1.writeU64 (M.poly[3]'(by omega)) 0 1 =
              .ok heap2 := by simpa [p3, hgcdMatPtr] using hwrite3
          have heq : heap' = heap2 := by
            have hrun' : (.ok (heap2, { M with len := #[1, 0, 0, 1] }) :
                RawExec (RawHeap × HgcdMat)) = .ok (heap', M') := by
              simpa [dense_upoly_zp__mat_one_ir, hvalid, hwrite0', hwrite3']
                using hrun
            exact (congrArg Prod.fst (Except.ok.inj hrun')).symm
          subst heap'
          have hsame0 := writeU64_preserves_prefix heap heap1 p0 guard 1
            guardLen 0 1 h0Guard (by omega) hwrite0
          have hsame3 := writeU64_preserves_prefix heap1 heap2 p3 guard 1
            guardLen 0 1 h3Guard (by omega) hwrite3
          exact sameU64Prefix_trans hsame0 hsame3

/-- End-to-end refinement of the exact initialization prefix of C++
`_hgcd_iter`: identity matrix construction followed by the ordered A and B
copies.  Every alias restriction below is a physical L1 memory condition. -/
theorem hgcdIterInit_refines (this : DenseUPolyZp)
    (M : HgcdMat) (A B T t : RawPtr UInt64) (lenT : Nat)
    (a : RawPtr UInt64) (lenA : Nat) (b : RawPtr UInt64) (lenB : Nat)
    (heap : RawHeap) (left right : Polynomial (ZMod this._p.toNat))
    (hM : M.Valid)
    (hp : 1 < this._p.toNat)
    (h0 : heap.ValidU64Slice (hgcdMatPtr M hM ⟨0, by omega⟩) 1)
    (h3 : heap.ValidU64Slice (hgcdMatPtr M hM ⟨3, by omega⟩) 1)
    (h03 : U64SlicesDisjoint (hgcdMatPtr M hM ⟨0, by omega⟩) 1
      (hgcdMatPtr M hM ⟨3, by omega⟩) 1)
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB)
    (hAa : U64SlicesDisjoint A lenA a lenA)
    (hBb : U64SlicesDisjoint B lenB b lenB)
    (hAb : U64SlicesDisjoint A lenA b lenB)
    (hBA : U64SlicesDisjoint B lenB A lenA)
    (h0a : U64SlicesDisjoint (hgcdMatPtr M hM ⟨0, by omega⟩) 1 a lenA)
    (h3a : U64SlicesDisjoint (hgcdMatPtr M hM ⟨3, by omega⟩) 1 a lenA)
    (h0b : U64SlicesDisjoint (hgcdMatPtr M hM ⟨0, by omega⟩) 1 b lenB)
    (h3b : U64SlicesDisjoint (hgcdMatPtr M hM ⟨3, by omega⟩) 1 b lenB)
    (hAMatrix : ∀ i : Fin 4, U64SlicesDisjoint A lenA
      (hgcdMatPtr M hM i) (identityEntryLen i))
    (hBMatrix : ∀ i : Fin 4, U64SlicesDisjoint B lenB
      (hgcdMatPtr M hM i) (identityEntryLen i))
    (hMatrixValid : ∀ i : Fin 4, heap.ValidU64Slice
      (hgcdMatPtr M hM i) (identityEntryLen i))
    (hLeft : RawDensePolyRep this heap a lenA left)
    (hRight : RawDensePolyRep this heap b lenB right) :
    ∃ initial, hgcdIterInit M A B T t lenT a lenA b lenB heap =
        .ok initial ∧
      initial.A = A ∧ initial.lenA = lenA ∧
      initial.B = B ∧ initial.lenB = lenB ∧
      initial.T = T ∧ initial.lenT = lenT ∧ initial.t = t ∧
      initial.sgn = 1 ∧
      ∃ hInitialM : initial.matrix.Valid,
        HgcdMatRawDenseRep this initial.heap initial.matrix
          (identityEntries this._p.toNat) hInitialM ∧
        RawDensePolyRep this initial.heap initial.A initial.lenA left ∧
        RawDensePolyRep this initial.heap initial.B initial.lenB right ∧
        CLPoly.Impl.StrictHGCDRefinement.HgcdTransform left right left right
          (identityEntries this._p.toNat ⟨0, by omega⟩)
          (identityEntries this._p.toNat ⟨1, by omega⟩)
          (identityEntries this._p.toNat ⟨2, by omega⟩)
          (identityEntries this._p.toNat ⟨3, by omega⟩) ∧
        CLPoly.Impl.StrictHGCDRefinement.HgcdSignedDet initial.sgn
          (identityEntries this._p.toNat ⟨0, by omega⟩)
          (identityEntries this._p.toNat ⟨1, by omega⟩)
          (identityEntries this._p.toNat ⟨2, by omega⟩)
          (identityEntries this._p.toNat ⟨3, by omega⟩) := by
  rcases matOne_refines M heap this._p.toNat hM h0 h3 h03 with
    ⟨heap1, M1, hone, hlayout1, hpoly1, hlen1, hread0, hread3,
      hM1, hMatrix1⟩
  have hsameA := matOne_preserves_prefix M heap heap1 M1 a lenA hM
    h0a h3a hone
  have hsameB := matOne_preserves_prefix M heap heap1 M1 b lenB hM
    h0b h3b hone
  have hLeft1 := rawDensePolyRep_of_same_prefix this heap heap1 a lenA left
    hlayout1 hsameA hLeft
  have hRight1 := rawDensePolyRep_of_same_prefix this heap heap1 b lenB right
    hlayout1 hsameB hRight
  have hPtr1 : ∀ i : Fin 4,
      hgcdMatPtr M1 hM1 i = hgcdMatPtr M hM i := by
    intro i
    fin_cases i <;> simp [hgcdMatPtr, hpoly1]
  have hLen1 : ∀ i : Fin 4, hgcdMatLen M1 hM1 i = identityEntryLen i := by
    intro i
    fin_cases i <;> simp [hgcdMatLen, hlen1, identityEntryLen]
  have hValidMatrix1 : ∀ i : Fin 4, heap1.ValidU64Slice
      (hgcdMatPtr M1 hM1 i) (hgcdMatLen M1 hM1 i) := by
    intro i
    rw [hPtr1 i, hLen1 i]
    exact (hlayout1 _ _).mp (hMatrixValid i)
  have hRawMatrix1 : HgcdMatRawDenseRep this heap1 M1
      (identityEntries this._p.toNat) hM1 := by
    intro i
    fin_cases i
    · simp only [identityEntries]
      rw [hLen1 ⟨0, by omega⟩]
      apply rawDensePolyRep_one_of_read_one this heap1 _ hp
      · convert hValidMatrix1 ⟨0, by omega⟩ using 1
        simpa [identityEntryLen] using (hLen1 ⟨0, by omega⟩).symm
      · rw [hPtr1 ⟨0, by omega⟩]
        exact hread0
    · simp only [identityEntries]
      rw [hLen1 ⟨1, by omega⟩]
      apply rawDensePolyRep_zero_length this heap1
      convert hValidMatrix1 ⟨1, by omega⟩ using 1
      simpa [identityEntryLen] using (hLen1 ⟨1, by omega⟩).symm
    · simp only [identityEntries]
      rw [hLen1 ⟨2, by omega⟩]
      apply rawDensePolyRep_zero_length this heap1
      convert hValidMatrix1 ⟨2, by omega⟩ using 1
      simpa [identityEntryLen] using (hLen1 ⟨2, by omega⟩).symm
    · simp only [identityEntries]
      rw [hLen1 ⟨3, by omega⟩]
      apply rawDensePolyRep_one_of_read_one this heap1 _ hp
      · convert hValidMatrix1 ⟨3, by omega⟩ using 1
        simpa [identityEntryLen] using (hLen1 ⟨3, by omega⟩).symm
      · rw [hPtr1 ⟨3, by omega⟩]
        exact hread3
  have hA1 := (hlayout1 A lenA).mp hA
  rcases copyU64_refines_rawDense this heap1 A a lenA left hA1 hAa
      hLeft1 with ⟨heap2, hcopyA, hlayout2, hA2⟩
  have hsameB2 := copyU64_preserves_prefix heap1 heap2 A a b lenA lenB
    hA1 hLeft1.1 hAb hcopyA
  have hRight2 := rawDensePolyRep_of_same_prefix this heap1 heap2 b lenB
    right hlayout2 hsameB2 hRight1
  have hAMatrix1 : ∀ i : Fin 4, U64SlicesDisjoint A lenA
      (hgcdMatPtr M1 hM1 i) (hgcdMatLen M1 hM1 i) := by
    intro i
    simpa [hPtr1 i, hLen1 i] using hAMatrix i
  rcases copyU64_preserves_hgcdMatRawDenseRep this heap1 heap2 A a lenA
      M1 (identityEntries this._p.toNat) hM1 hA1 hLeft1.1 hAMatrix1
      hcopyA hRawMatrix1 with ⟨_, hMatrix2⟩
  have hB2 := (hlayout2 B lenB).mp ((hlayout1 B lenB).mp hB)
  rcases copyU64_refines_rawDense this heap2 B b lenB right hB2 hBb
      hRight2 with ⟨heap3, hcopyB, hlayout3, hB3⟩
  have hsameA3 := copyU64_preserves_prefix heap2 heap3 B b A lenB lenA
    hB2 hRight2.1 hBA hcopyB
  have hA3 := rawDensePolyRep_of_same_prefix this heap2 heap3 A lenA left
    hlayout3 hsameA3 hA2
  have hValidMatrix2 : ∀ i : Fin 4, heap2.ValidU64Slice
      (hgcdMatPtr M1 hM1 i) (hgcdMatLen M1 hM1 i) := by
    intro i
    exact (hlayout2 _ _).mp (hValidMatrix1 i)
  have hBMatrix1 : ∀ i : Fin 4, U64SlicesDisjoint B lenB
      (hgcdMatPtr M1 hM1 i) (hgcdMatLen M1 hM1 i) := by
    intro i
    simpa [hPtr1 i, hLen1 i] using hBMatrix i
  rcases copyU64_preserves_hgcdMatRawDenseRep this heap2 heap3 B b lenB
      M1 (identityEntries this._p.toNat) hM1 hB2 hRight2.1 hBMatrix1
      hcopyB hMatrix2 with
    ⟨_, hMatrix3⟩
  let initial : HgcdIterState := {
    heap := heap3, matrix := M1, A := A, lenA := lenA, B := B,
    lenB := lenB, T := T, lenT := lenT, t := t, sgn := 1 }
  refine ⟨initial, ?_, rfl, rfl, rfl, rfl, rfl, rfl, rfl, rfl,
    hM1, ?_, ?_, ?_, ?_, ?_⟩
  · simp [hgcdIterInit, hone, hcopyA, hcopyB, initial]
  · simpa [initial] using hMatrix3
  · simpa [initial] using hA3
  · simpa [initial] using hB3
  · simp [CLPoly.Impl.StrictHGCDRefinement.HgcdTransform, identityEntries]
  · simp [initial, CLPoly.Impl.StrictHGCDRefinement.HgcdSignedDet,
      identityEntries]

/-- A stopped generated HGCD loop returns its state unchanged. -/
theorem hgcdIterLoop_stop (this : DenseUPolyZp) (m : Nat)
    (Q : RawPtr UInt64) (W3 : RawPtr Word3) (scratch : RawPtr UInt64)
    (state : HgcdIterState) (hstop : state.lenB < m + 1) :
    hgcdIterLoop this m Q W3 scratch state = .ok state := by
  rw [hgcdIterLoop]
  simp [show ¬state.lenB ≥ m + 1 by omega]

/-- Every successful nonterminal iteration is exactly one generated divrem,
the source pointer rotation, two generated row updates, and the recursive
tail call on the shorter remainder. -/
theorem hgcdIterLoop_step_shape (this : DenseUPolyZp) (m : Nat)
    (Q : RawPtr UInt64) (W3 : RawPtr Word3) (scratch : RawPtr UInt64)
    (state final : HgcdIterState) (hguard : state.lenB ≥ m + 1)
    (hrun : hgcdIterLoop this m Q W3 scratch state = .ok final) :
    ∃ heap1 lenQ lenR row23 row01,
      Generated.StrictDivrem.dense_upoly_zp__poly_divrem_ir this Q state.T
        state.A state.lenA state.B state.lenB W3 state.heap =
          .ok (heap1, lenQ, lenR) ∧
      dense_upoly_zp__mat_row_update_ir this state.matrix
        ⟨2, by omega⟩ ⟨3, by omega⟩ Q lenQ state.A state.lenT state.t
        scratch heap1 = .ok row23 ∧
      dense_upoly_zp__mat_row_update_ir this row23.matrix
        ⟨0, by omega⟩ ⟨1, by omega⟩ Q lenQ row23.T row23.lenT row23.t
        scratch row23.heap = .ok row01 ∧
      hgcdIterLoop this m Q W3 scratch {
        heap := row01.heap
        matrix := row01.matrix
        A := state.B
        lenA := state.lenB
        B := state.T
        lenB := lenR
        T := row01.T
        lenT := row01.lenT
        t := row01.t
        sgn := -state.sgn
      } = .ok final ∧ lenR < state.lenB := by
  rw [hgcdIterLoop] at hrun
  simp only [if_pos hguard] at hrun
  split at hrun
  next fault hdiv => simp at hrun
  next =>
    rename_i heap1 lenQ lenR hdiv
    split at hrun
    next fault hrow23 => simp at hrun
    next row23 hrow23 =>
      split at hrun
      next fault hrow01 => simp at hrun
      next row01 hrow01 =>
        have hlt := Generated.StrictDivrem.polyDivrem_remainder_lt this Q
          state.T state.A
          state.lenA state.B state.lenB W3 state.heap heap1 lenQ lenR hdiv
        exact ⟨heap1, lenQ, lenR, row23, row01, hdiv, hrow23,
          hrow01, hrun, hlt⟩

/-- Semantic divrem closure for one successful nonterminal HGCD iteration.
All polynomial facts below concern the exact divrem result exposed by
`hgcdIterLoop_step_shape`; no second or specification-level execution is
substituted. -/
theorem hgcdIterLoop_step_divrem_refines (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (m : Nat) (Q : RawPtr UInt64) (W3 : RawPtr Word3)
    (scratch : RawPtr UInt64) (state final : HgcdIterState)
    (dividend divisor : Polynomial (ZMod this._p.toNat))
    (hguard : state.lenB ≥ m + 1)
    (hARep : RawDensePolyRep this state.heap state.A state.lenA dividend)
    (hBRep : RawDensePolyRep this state.heap state.B state.lenB divisor)
    (hQ : state.heap.ValidU64Slice Q
      (state.lenA - (state.lenB - 1)))
    (hR : state.heap.ValidU64Slice state.T
      (Nat.min state.lenA (state.lenB - 1)))
    (hW3 : state.heap.ValidWord3Slice W3 state.lenA)
    (hqCapacity : state.lenA - (state.lenB - 1) < limbBase)
    (hRA : state.T.region ≠ state.A.region)
    (hWA : W3.region ≠ state.A.region)
    (hWB : W3.region ≠ state.B.region)
    (hQB : Q.region ≠ state.B.region)
    (hQW : Q.region ≠ W3.region)
    (hRW : state.T.region ≠ W3.region)
    (hRQ : state.T.region ≠ Q.region)
    (hRB : state.T.region ≠ state.B.region)
    (hcfg : DensePreinvConfigured this)
    (hrun : hgcdIterLoop this m Q W3 scratch state = .ok final) :
    ∃ heap1 lenQ lenR quotient remainder row23 row01,
      Generated.StrictDivrem.dense_upoly_zp__poly_divrem_ir this Q state.T
        state.A state.lenA state.B state.lenB W3 state.heap =
          .ok (heap1, lenQ, lenR) ∧
      RawDensePolyRep this heap1 Q lenQ quotient ∧
      RawDensePolyRep this heap1 state.B state.lenB divisor ∧
      RawDensePolyRep this heap1 state.T lenR remainder ∧
      dividend = quotient * divisor + remainder ∧
      normalize (EuclideanDomain.gcd dividend divisor) =
        normalize (EuclideanDomain.gcd divisor remainder) ∧
      dense_upoly_zp__mat_row_update_ir this state.matrix
        ⟨2, by omega⟩ ⟨3, by omega⟩ Q lenQ state.A state.lenT state.t
        scratch heap1 = .ok row23 ∧
      dense_upoly_zp__mat_row_update_ir this row23.matrix
        ⟨0, by omega⟩ ⟨1, by omega⟩ Q lenQ row23.T row23.lenT row23.t
        scratch row23.heap = .ok row01 ∧
      hgcdIterLoop this m Q W3 scratch {
        heap := row01.heap
        matrix := row01.matrix
        A := state.B
        lenA := state.lenB
        B := state.T
        lenB := lenR
        T := row01.T
        lenT := row01.lenT
        t := row01.t
        sgn := -state.sgn
      } = .ok final ∧
      lenR < state.lenB := by
  have hlenB : 0 < state.lenB := by omega
  rcases hgcdIterLoop_step_shape this m Q W3 scratch state final hguard hrun
      with ⟨heap1, lenQ, lenR, row23, row01, hdiv, hrow23, hrow01,
        htail, hlt⟩
  rcases polyDivrem_next_state this Q state.T state.A state.B state.lenA
      state.lenB W3 state.heap dividend divisor hlenB hARep hBRep hQ hR
      hW3 hqCapacity hRA hWA hWB hQB hQW hRW hRQ hRB hcfg with
    ⟨semanticHeap, semanticLenQ, semanticLenR, quotient, remainder,
      hsemantic, hQRep, hBRep1, hRRep, hdivision, hgcd, _, _, _, _⟩
  have heq : (semanticHeap, semanticLenQ, semanticLenR) =
      (heap1, lenQ, lenR) := Except.ok.inj (hsemantic.symm.trans hdiv)
  cases heq
  exact ⟨heap1, lenQ, lenR, quotient, remainder, row23, row01, hdiv,
    hQRep, hBRep1, hRRep, hdivision, hgcd, hrow23, hrow01, htail, hlt⟩

/-- The actual generated divrem call preserves all four live HGCD matrix
entries when its three write allocations are physically distinct from every
matrix allocation. -/
theorem polyDivrem_preserves_hgcdMatRawDenseRep (this : DenseUPolyZp)
    (M : HgcdMat) (Q R A B : RawPtr UInt64) (lenA lenB : Nat)
    (W3 : RawPtr Word3) (heap heap' : RawHeap) (lenQ lenR : Nat)
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat)) (hM : M.Valid)
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB)
    (hQ : heap.ValidU64Slice Q (lenA - (lenB - 1)))
    (hR : heap.ValidU64Slice R (Nat.min lenA (lenB - 1)))
    (hW3 : heap.ValidWord3Slice W3 lenA)
    (hQMatrix : ∀ i : Fin 4,
      Q.region ≠ (hgcdMatPtr M hM i).region)
    (hRMatrix : ∀ i : Fin 4,
      R.region ≠ (hgcdMatPtr M hM i).region)
    (hW3Matrix : ∀ i : Fin 4,
      W3.region ≠ (hgcdMatPtr M hM i).region)
    (hMatrix : HgcdMatRawDenseRep this heap M entries hM)
    (hrun : Generated.StrictDivrem.dense_upoly_zp__poly_divrem_ir this Q R
      A lenA B lenB W3 heap = .ok (heap', lenQ, lenR)) :
    HgcdMatRawDenseRep this heap' M entries hM := by
  intro i
  rcases polyDivrem_preserves_u64_region_ne this Q R A B
      (hgcdMatPtr M hM i) lenA lenB (hgcdMatLen M hM i) W3 heap heap'
      lenQ lenR hA hB hQ hR hW3 (hMatrix i).1 (hQMatrix i)
      (hRMatrix i) (hW3Matrix i) hrun with ⟨hlayout, hsame⟩
  exact rawDensePolyRep_of_same_prefix this heap heap'
    (hgcdMatPtr M hM i) (hgcdMatLen M hM i) (entries i) hlayout hsame
    (hMatrix i)

/-- The source's zero-quotient/zero-entry branch performs exactly the two
matrix-entry swaps and no heap access.  This exposes the real descriptor
state consumed by the next HGCD iteration. -/
theorem matRowUpdate_zero_exec (this : DenseUPolyZp) (M : HgcdMat)
    (i0 i1 : Fin 4) (Q : RawPtr UInt64) (lenQ : Nat)
    (T : RawPtr UInt64) (lenT : Nat) (t scratch : RawPtr UInt64)
    (heap : RawHeap) (hM : M.Valid)
    (hne : i0 ≠ i1)
    (hzero : lenQ = 0 ∨ hgcdMatLen M hM i0 = 0) :
    ∃ result,
      dense_upoly_zp__mat_row_update_ir this M i0 i1 Q lenQ T lenT t
          scratch heap = .ok result ∧
      result.heap = heap ∧ result.T = T ∧ result.lenT = lenT ∧
      result.t = t ∧
      ∃ hResult : result.matrix.Valid,
        hgcdMatPtr result.matrix hResult i0 = hgcdMatPtr M hM i1 ∧
        hgcdMatLen result.matrix hResult i0 = hgcdMatLen M hM i1 ∧
        hgcdMatPtr result.matrix hResult i1 = hgcdMatPtr M hM i0 ∧
        hgcdMatLen result.matrix hResult i1 = hgcdMatLen M hM i0 ∧
        ∀ j : Fin 4, j ≠ i0 → j ≠ i1 →
          hgcdMatPtr result.matrix hResult j = hgcdMatPtr M hM j ∧
          hgcdMatLen result.matrix hResult j = hgcdMatLen M hM j := by
  have hvalid : M.poly.size = 4 ∧ M.len.size = 4 := by
    simpa [HgcdMat.Valid] using hM
  have hneVal : i0.val ≠ i1.val := by
    intro heq
    exact hne (Fin.ext heq)
  let p0 := hgcdMatPtr M hM i0
  let p1 := hgcdMatPtr M hM i1
  let l0 := hgcdMatLen M hM i0
  let l1 := hgcdMatLen M hM i1
  let poly' := (M.poly.set i1.val p0 (by omega)).set i0.val p1 (by simp; omega)
  let len' := (M.len.set i1.val l0 (by omega)).set i0.val l1 (by simp; omega)
  let result := MatRowUpdateResult.mk heap
    ({ poly := poly', len := len' } : HgcdMat) T lenT t
  have hbranch : ¬(lenQ ≠ 0 ∧ l0 ≠ 0) := by
    intro h
    rcases hzero with hq | hentry
    · exact h.1 hq
    · exact h.2 (by simpa [l0, hgcdMatLen] using hentry)
  have hbranch' : ¬(lenQ ≠ 0 ∧ M.len[i0.val]'(by omega) ≠ 0) := by
    simpa [l0, hgcdMatLen] using hbranch
  have hrun : dense_upoly_zp__mat_row_update_ir this M i0 i1 Q lenQ T
      lenT t scratch heap = .ok result := by
    simp [dense_upoly_zp__mat_row_update_ir, hvalid, hbranch', result,
      poly', len', p0, p1, l0, l1, hgcdMatPtr, hgcdMatLen]
  have hResult : result.matrix.Valid := by
    simp [result, HgcdMat.Valid, poly', len', hvalid]
  refine ⟨result, hrun, rfl, rfl, rfl, rfl, hResult, ?_, ?_, ?_, ?_, ?_⟩
  · simp [result, poly', p1, hgcdMatPtr]
  · simp [result, len', l1, hgcdMatLen]
  · simp only [result, poly', hgcdMatPtr]
    rw [Array.getElem_set_ne
      (by simpa using (show i0.val < M.poly.size by
        rw [hvalid.1]; exact i0.isLt))
      (by simpa using (show i1.val < M.poly.size by
        rw [hvalid.1]; exact i1.isLt)) hneVal,
      Array.getElem_set_self]
    simp [p0, hgcdMatPtr]
  · simp only [result, len', hgcdMatLen]
    rw [Array.getElem_set_ne
      (by simpa using (show i0.val < M.len.size by
        rw [hvalid.2]; exact i0.isLt))
      (by simpa using (show i1.val < M.len.size by
        rw [hvalid.2]; exact i1.isLt)) hneVal,
      Array.getElem_set_self]
    simp [l0, hgcdMatLen]
  · intro j hji0 hji1
    have hi0j : i0.val ≠ j.val := by
      intro heq
      exact hji0 (Fin.ext heq.symm)
    have hi1j : i1.val ≠ j.val := by
      intro heq
      exact hji1 (Fin.ext heq.symm)
    have hi0PolySet : i0.val < (M.poly.set i1.val p0 (by omega)).size := by
      simpa [hvalid.1] using i0.isLt
    have hjPolySet : j.val < (M.poly.set i1.val p0 (by omega)).size := by
      simpa [hvalid.1] using j.isLt
    have hi1Poly : i1.val < M.poly.size := by simpa [hvalid.1] using i1.isLt
    have hjPoly : j.val < M.poly.size := by simpa [hvalid.1] using j.isLt
    have hi0LenSet : i0.val < (M.len.set i1.val l0 (by omega)).size := by
      simpa [hvalid.2] using i0.isLt
    have hjLenSet : j.val < (M.len.set i1.val l0 (by omega)).size := by
      simpa [hvalid.2] using j.isLt
    have hi1Len : i1.val < M.len.size := by simpa [hvalid.2] using i1.isLt
    have hjLen : j.val < M.len.size := by simpa [hvalid.2] using j.isLt
    constructor
    · simp only [result, poly', hgcdMatPtr]
      rw [Array.getElem_set_ne hi0PolySet hjPolySet hi0j,
        Array.getElem_set_ne hi1Poly hjPoly hi1j]
    · simp only [result, len', hgcdMatLen]
      rw [Array.getElem_set_ne hi0LenSet hjLenSet hi0j,
        Array.getElem_set_ne hi1Len hjLen hi1j]

/-- Complete semantic refinement of the source's zero row-update branch.
The physical descriptor swap represents `[entry1 + quotient*entry0,
entry0]`; the vanishing factor is derived from the actual zero-length raw
representation, not assumed at L2. -/
theorem matRowUpdate_zero_refines (this : DenseUPolyZp)
    (M : HgcdMat) (i0 i1 : Fin 4) (Q : RawPtr UInt64) (lenQ : Nat)
    (T : RawPtr UInt64) (lenT : Nat) (t scratch : RawPtr UInt64)
    (heap : RawHeap) (result : MatRowUpdateResult) (hM : M.Valid)
    (quotient entry0 entry1 : Polynomial (ZMod this._p.toNat))
    (hne : i0 ≠ i1)
    (hzero : lenQ = 0 ∨ hgcdMatLen M hM i0 = 0)
    (hQRep : RawDensePolyRep this heap Q lenQ quotient)
    (hEntry0Rep : RawDensePolyRep this heap (hgcdMatPtr M hM i0)
      (hgcdMatLen M hM i0) entry0)
    (hEntry1Rep : RawDensePolyRep this heap (hgcdMatPtr M hM i1)
      (hgcdMatLen M hM i1) entry1)
    (hrun : dense_upoly_zp__mat_row_update_ir this M i0 i1 Q lenQ T
      lenT t scratch heap = .ok result) :
    ∃ hResult : result.matrix.Valid,
      RawDensePolyRep this result.heap
        (hgcdMatPtr result.matrix hResult i0)
        (hgcdMatLen result.matrix hResult i0)
        (entry1 + quotient * entry0) ∧
      RawDensePolyRep this result.heap
        (hgcdMatPtr result.matrix hResult i1)
        (hgcdMatLen result.matrix hResult i1) entry0 := by
  rcases matRowUpdate_zero_exec this M i0 i1 Q lenQ T lenT t scratch
      heap hM hne hzero with
    ⟨actual, hactual, hheap, _, _, _, hActualM, hptr0, hlen0,
      hptr1, hlen1, _⟩
  have heq : actual = result := Except.ok.inj (hactual.symm.trans hrun)
  subst actual
  have hvanish : quotient = 0 ∨ entry0 = 0 := by
    rcases hzero with hq | he
    · left
      subst lenQ
      exact slicePolyRep_zero_length heap Q this._p.toNat quotient
        hQRep.2.2.1
    · right
      have hrep0 : SlicePolyRep heap (hgcdMatPtr M hM i0) 0
          this._p.toNat entry0 := by
        simpa [he] using hEntry0Rep.2.2.1
      exact slicePolyRep_zero_length heap (hgcdMatPtr M hM i0)
        this._p.toNat entry0 hrep0
  refine ⟨hActualM, ?_, ?_⟩
  · rcases hvanish with hq | he
    · simpa [hheap, hptr0, hlen0, hq] using hEntry1Rep
    · simpa [hheap, hptr0, hlen0, he] using hEntry1Rep
  · simpa [hheap, hptr1, hlen1] using hEntry0Rep

/-- Successful execution of the nonzero source branch can only arise from
the actual `_mul` followed by the actual `_poly_add`; the returned pointer
state is exactly the three source swaps. -/
theorem matRowUpdate_nonzero_success_shape (this : DenseUPolyZp)
    (M : HgcdMat) (i0 i1 : Fin 4) (Q : RawPtr UInt64) (lenQ : Nat)
    (T : RawPtr UInt64) (lenT : Nat) (t scratch : RawPtr UInt64)
    (heap : RawHeap) (result : MatRowUpdateResult) (hM : M.Valid)
    (hne : i0 ≠ i1) (hQ : lenQ ≠ 0)
    (hEntry : hgcdMatLen M hM i0 ≠ 0)
    (hrun : dense_upoly_zp__mat_row_update_ir this M i0 i1 Q lenQ T
      lenT t scratch heap = .ok result) :
    ∃ heap1 heap2 sumLen,
      Generated.StrictMul.dense_upoly_zp__mul_ir this T
        (if lenQ ≥ hgcdMatLen M hM i0 then Q else hgcdMatPtr M hM i0)
        (if lenQ ≥ hgcdMatLen M hM i0 then lenQ else hgcdMatLen M hM i0)
        (if lenQ ≥ hgcdMatLen M hM i0 then hgcdMatPtr M hM i0 else Q)
        (if lenQ ≥ hgcdMatLen M hM i0 then hgcdMatLen M hM i0 else lenQ)
        scratch heap = .ok heap1 ∧
      Generated.StrictPolyAddSub.dense_upoly_zp__poly_add_ir this t
        (hgcdMatPtr M hM i1) (hgcdMatLen M hM i1) T
        (lenQ + hgcdMatLen M hM i0 - 1) heap1 = .ok (heap2, sumLen) ∧
      result.heap = heap2 ∧ result.T = T ∧
      result.lenT = lenQ + hgcdMatLen M hM i0 - 1 ∧
      result.t = hgcdMatPtr M hM i1 ∧
      ∃ hResult : result.matrix.Valid,
        hgcdMatPtr result.matrix hResult i0 = t ∧
        hgcdMatLen result.matrix hResult i0 = sumLen ∧
        hgcdMatPtr result.matrix hResult i1 = hgcdMatPtr M hM i0 ∧
        hgcdMatLen result.matrix hResult i1 = hgcdMatLen M hM i0 ∧
        ∀ j : Fin 4, j ≠ i0 → j ≠ i1 →
          hgcdMatPtr result.matrix hResult j = hgcdMatPtr M hM j ∧
          hgcdMatLen result.matrix hResult j = hgcdMatLen M hM j := by
  have hvalid : M.poly.size = 4 ∧ M.len.size = 4 := by
    simpa [HgcdMat.Valid] using hM
  have hneVal : i0.val ≠ i1.val := by
    intro heq
    apply hne
    exact Fin.ext heq
  let p0 := hgcdMatPtr M hM i0
  let p1 := hgcdMatPtr M hM i1
  let l0 := hgcdMatLen M hM i0
  let l1 := hgcdMatLen M hM i1
  let left := if lenQ ≥ l0 then Q else p0
  let leftLen := if lenQ ≥ l0 then lenQ else l0
  let right := if lenQ ≥ l0 then p0 else Q
  let rightLen := if lenQ ≥ l0 then l0 else lenQ
  have hEntry' : M.len[i0.val]'(by omega) ≠ 0 := by
    simpa [hgcdMatLen] using hEntry
  cases hmul : Generated.StrictMul.dense_upoly_zp__mul_ir this T left
      leftLen right rightLen scratch heap with
  | error fault =>
      have hmul' : Generated.StrictMul.dense_upoly_zp__mul_ir this T
          (if lenQ ≥ M.len[i0.val]'(by omega) then Q else M.poly[i0.val]'(by omega))
          (if lenQ ≥ M.len[i0.val]'(by omega) then lenQ else M.len[i0.val]'(by omega))
          (if lenQ ≥ M.len[i0.val]'(by omega) then M.poly[i0.val]'(by omega) else Q)
          (if lenQ ≥ M.len[i0.val]'(by omega) then M.len[i0.val]'(by omega) else lenQ)
          scratch heap = .error fault := by
        simpa [left, leftLen, right, rightLen, p0, l0, hgcdMatPtr,
          hgcdMatLen] using hmul
      simp [dense_upoly_zp__mat_row_update_ir, hvalid, hQ, hEntry', hmul'] at hrun
  | ok heap1 =>
      have hmul' : Generated.StrictMul.dense_upoly_zp__mul_ir this T
          (if lenQ ≥ M.len[i0.val]'(by omega) then Q else M.poly[i0.val]'(by omega))
          (if lenQ ≥ M.len[i0.val]'(by omega) then lenQ else M.len[i0.val]'(by omega))
          (if lenQ ≥ M.len[i0.val]'(by omega) then M.poly[i0.val]'(by omega) else Q)
          (if lenQ ≥ M.len[i0.val]'(by omega) then M.len[i0.val]'(by omega) else lenQ)
          scratch heap = .ok heap1 := by
        simpa [left, leftLen, right, rightLen, p0, l0, hgcdMatPtr,
          hgcdMatLen] using hmul
      cases hadd : Generated.StrictPolyAddSub.dense_upoly_zp__poly_add_ir
          this t p1 l1 T (lenQ + l0 - 1) heap1 with
      | error fault =>
          have hadd' : Generated.StrictPolyAddSub.dense_upoly_zp__poly_add_ir
              this t (M.poly[i1.val]'(by omega)) (M.len[i1.val]'(by omega)) T
              (lenQ + M.len[i0.val]'(by omega) - 1) heap1 = .error fault := by
            simpa [p1, l0, l1, hgcdMatPtr, hgcdMatLen] using hadd
          simp [dense_upoly_zp__mat_row_update_ir, hvalid, hQ, hEntry',
            hmul', hadd'] at hrun
      | ok pair =>
          rcases pair with ⟨heap2, sumLen⟩
          have hadd' : Generated.StrictPolyAddSub.dense_upoly_zp__poly_add_ir
              this t (M.poly[i1.val]'(by omega)) (M.len[i1.val]'(by omega)) T
              (lenQ + M.len[i0.val]'(by omega) - 1) heap1 =
                .ok (heap2, sumLen) := by
            simpa [p1, l0, l1, hgcdMatPtr, hgcdMatLen] using hadd
          have heq : result = MatRowUpdateResult.mk heap2
              ({
                poly := (M.poly.set i1.val p0 (by omega)).set i0.val t
                  (by simp; omega)
                len := (M.len.set i1.val l0 (by omega)).set i0.val sumLen
                  (by simp; omega)
              } : HgcdMat) T (lenQ + l0 - 1) p1 := by
            have hrun' := hrun
            simp [dense_upoly_zp__mat_row_update_ir, hvalid, hQ, hEntry',
              hmul', hadd'] at hrun'
            simpa [p0, p1, l0, hgcdMatPtr, hgcdMatLen] using hrun'.symm
          subst result
          refine ⟨heap1, heap2, sumLen, ?_, ?_, rfl, rfl, ?_, ?_, ?_⟩
          · simpa [left, leftLen, right, rightLen, p0, l0,
              hgcdMatPtr, hgcdMatLen] using hmul
          · simpa [p1, l0, l1, hgcdMatPtr, hgcdMatLen] using hadd
          · simp [l0, hgcdMatLen]
          · simp [p1, hgcdMatPtr]
          · have hResult :
                ({
                  poly := (M.poly.set i1.val p0 (by omega)).set i0.val t
                    (by simp; omega)
                  len := (M.len.set i1.val l0 (by omega)).set i0.val sumLen
                    (by simp; omega)
                } : HgcdMat).Valid := by
              simp [HgcdMat.Valid, hvalid]
            refine ⟨hResult, ?_, ?_, ?_, ?_, ?_⟩
            · simp [hgcdMatPtr]
            · simp [hgcdMatLen]
            · simp only [hgcdMatPtr]
              rw [Array.getElem_set_ne
                (by simpa using (show i0.val < M.poly.size by
                  rw [hvalid.1]; exact i0.isLt))
                (by simpa using (show i1.val < M.poly.size by
                  rw [hvalid.1]; exact i1.isLt)) hneVal,
                Array.getElem_set_self]
              rfl
            · simp only [hgcdMatLen]
              rw [Array.getElem_set_ne
                (by simpa using (show i0.val < M.len.size by
                  rw [hvalid.2]; exact i0.isLt))
                (by simpa using (show i1.val < M.len.size by
                  rw [hvalid.2]; exact i1.isLt)) hneVal,
                Array.getElem_set_self]
              rfl
            · intro j hji0 hji1
              have hi0j : i0.val ≠ j.val := by
                intro heq
                exact hji0 (Fin.ext heq.symm)
              have hi1j : i1.val ≠ j.val := by
                intro heq
                exact hji1 (Fin.ext heq.symm)
              have hi0PolySet : i0.val <
                  (M.poly.set i1.val p0 (by omega)).size := by
                simpa [hvalid.1] using i0.isLt
              have hjPolySet : j.val <
                  (M.poly.set i1.val p0 (by omega)).size := by
                simpa [hvalid.1] using j.isLt
              have hi1Poly : i1.val < M.poly.size := by
                simpa [hvalid.1] using i1.isLt
              have hjPoly : j.val < M.poly.size := by
                simpa [hvalid.1] using j.isLt
              have hi0LenSet : i0.val <
                  (M.len.set i1.val l0 (by omega)).size := by
                simpa [hvalid.2] using i0.isLt
              have hjLenSet : j.val <
                  (M.len.set i1.val l0 (by omega)).size := by
                simpa [hvalid.2] using j.isLt
              have hi1Len : i1.val < M.len.size := by
                simpa [hvalid.2] using i1.isLt
              have hjLen : j.val < M.len.size := by
                simpa [hvalid.2] using j.isLt
              constructor
              · simp only [hgcdMatPtr]
                rw [Array.getElem_set_ne hi0PolySet hjPolySet hi0j,
                  Array.getElem_set_ne hi1Poly hjPoly hi1j]
              · simp only [hgcdMatLen]
                rw [Array.getElem_set_ne hi0LenSet hjLenSet hi0j,
                  Array.getElem_set_ne hi1Len hjLen hi1j]

/-- Branch-complete descriptor frame for `_mat_row_update`: indices other
than the selected pair retain exactly their source pointer and length. -/
theorem matRowUpdate_preserves_other_descriptor (this : DenseUPolyZp)
    (M : HgcdMat) (i0 i1 j : Fin 4) (Q : RawPtr UInt64) (lenQ : Nat)
    (T : RawPtr UInt64) (lenT : Nat) (t scratch : RawPtr UInt64)
    (heap : RawHeap) (result : MatRowUpdateResult) (hM : M.Valid)
    (hne : i0 ≠ i1) (hji0 : j ≠ i0) (hji1 : j ≠ i1)
    (hrun : dense_upoly_zp__mat_row_update_ir this M i0 i1 Q lenQ T
      lenT t scratch heap = .ok result) :
    ∃ hResult : result.matrix.Valid,
      hgcdMatPtr result.matrix hResult j = hgcdMatPtr M hM j ∧
      hgcdMatLen result.matrix hResult j = hgcdMatLen M hM j := by
  by_cases hzero : lenQ = 0 ∨ hgcdMatLen M hM i0 = 0
  · rcases matRowUpdate_zero_exec this M i0 i1 Q lenQ T lenT t scratch
      heap hM hne hzero with
    ⟨actual, hactual, _, _, _, _, hActualM, _, _, _, _, hother⟩
    have heq : actual = result := Except.ok.inj (hactual.symm.trans hrun)
    subst actual
    exact ⟨hActualM, hother j hji0 hji1⟩
  · rcases matRowUpdate_nonzero_success_shape this M i0 i1 Q lenQ T lenT
      t scratch heap result hM hne (by omega) (by omega) hrun with
    ⟨_, _, _, _, _, _, _, _, _, hResult, _, _, _, _, hother⟩
    exact ⟨hResult, hother j hji0 hji1⟩

/-- Descriptor frame specialized to the two source row updates in one HGCD
iteration: `(2,3)` leaves `(0,1)` unchanged, then `(0,1)` leaves `(2,3)`
unchanged. -/
theorem hgcdTwoRowUpdates_descriptor_frame (this : DenseUPolyZp)
    (M : HgcdMat) (Q : RawPtr UInt64) (lenQ : Nat)
    (T : RawPtr UInt64) (lenT : Nat) (t scratch : RawPtr UInt64)
    (heap : RawHeap) (row23 row01 : MatRowUpdateResult) (hM : M.Valid)
    (hrow23 : dense_upoly_zp__mat_row_update_ir this M
      ⟨2, by omega⟩ ⟨3, by omega⟩ Q lenQ T lenT t scratch heap = .ok row23)
    (hrow01 : dense_upoly_zp__mat_row_update_ir this row23.matrix
      ⟨0, by omega⟩ ⟨1, by omega⟩ Q lenQ row23.T row23.lenT row23.t
      scratch row23.heap = .ok row01) :
    ∃ h23 : row23.matrix.Valid, ∃ h01 : row01.matrix.Valid,
      (∀ j : Fin 2,
        hgcdMatPtr row23.matrix h23 ⟨j.val, by omega⟩ =
          hgcdMatPtr M hM ⟨j.val, by omega⟩ ∧
        hgcdMatLen row23.matrix h23 ⟨j.val, by omega⟩ =
          hgcdMatLen M hM ⟨j.val, by omega⟩) ∧
      (∀ j : Fin 2,
        hgcdMatPtr row01.matrix h01 ⟨j.val + 2, by omega⟩ =
          hgcdMatPtr row23.matrix h23 ⟨j.val + 2, by omega⟩ ∧
        hgcdMatLen row01.matrix h01 ⟨j.val + 2, by omega⟩ =
          hgcdMatLen row23.matrix h23 ⟨j.val + 2, by omega⟩) := by
  rcases matRowUpdate_preserves_other_descriptor this M
      (2 : Fin 4) (3 : Fin 4) (0 : Fin 4) Q lenQ T lenT t
      scratch heap row23 hM (by decide) (by decide) (by decide)
      hrow23 with
    ⟨h23, hzero⟩
  rcases matRowUpdate_preserves_other_descriptor this row23.matrix
      (0 : Fin 4) (1 : Fin 4) (2 : Fin 4) Q lenQ row23.T
      row23.lenT row23.t scratch row23.heap row01 h23 (by decide)
      (by decide) (by decide) hrow01 with ⟨h01, htwo⟩
  refine ⟨h23, h01, ?_, ?_⟩
  · intro j
    fin_cases j
    · exact hzero
    · rcases matRowUpdate_preserves_other_descriptor this M
        (2 : Fin 4) (3 : Fin 4) (1 : Fin 4) Q lenQ T lenT t
        scratch heap row23 hM (by decide) (by decide) (by decide)
        hrow23 with
        ⟨h23', hone⟩
      simpa only using hone
  · intro j
    fin_cases j
    · exact htwo
    · rcases matRowUpdate_preserves_other_descriptor this row23.matrix
        (0 : Fin 4) (1 : Fin 4) (3 : Fin 4) Q lenQ row23.T
        row23.lenT row23.t scratch row23.heap row01 h23 (by decide)
        (by decide) (by decide) hrow01 with ⟨h01', hthree⟩
      simpa only using hthree

/-- The multiplication call exposed by the real nonzero row-update branch
computes exactly the quotient times the old `i0` entry.  The conditional
operand order is the one used by C++ to satisfy `_mul`'s length contract;
commutativity is used only after refining that actual dispatcher call. -/
theorem matRowUpdate_mul_result (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (M : HgcdMat) (hM : M.Valid) (i0 : Fin 4)
    (Q T scratch : RawPtr UInt64) (lenQ : Nat) (heap heap1 : RawHeap)
    (quotient entry0 : Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (hQpos : 0 < lenQ) (hEntryPos : 0 < hgcdMatLen M hM i0)
    (hLenWord : max lenQ (hgcdMatLen M hM i0) < limbBase)
    (hT : heap.ValidU64Slice T
      (2 * max lenQ (hgcdMatLen M hM i0) - 1))
    (hScratch : heap.ValidU64Slice scratch
      (8 * max lenQ (hgcdMatLen M hM i0)))
    (hTQ : U64SlicesDisjoint T
      (2 * max lenQ (hgcdMatLen M hM i0) - 1) Q lenQ)
    (hTEntry : U64SlicesDisjoint T
      (2 * max lenQ (hgcdMatLen M hM i0) - 1)
      (hgcdMatPtr M hM i0) (hgcdMatLen M hM i0))
    (hTScratch : U64SlicesDisjoint T
      (2 * max lenQ (hgcdMatLen M hM i0) - 1) scratch
      (8 * max lenQ (hgcdMatLen M hM i0)))
    (hScratchQ : U64SlicesDisjoint scratch
      (8 * max lenQ (hgcdMatLen M hM i0)) Q lenQ)
    (hScratchEntry : U64SlicesDisjoint scratch
      (8 * max lenQ (hgcdMatLen M hM i0))
      (hgcdMatPtr M hM i0) (hgcdMatLen M hM i0))
    (hQRep : RawDensePolyRep this heap Q lenQ quotient)
    (hEntryRep : RawDensePolyRep this heap (hgcdMatPtr M hM i0)
      (hgcdMatLen M hM i0) entry0)
    (hmul : Generated.StrictMul.dense_upoly_zp__mul_ir this T
      (if lenQ ≥ hgcdMatLen M hM i0 then Q else hgcdMatPtr M hM i0)
      (if lenQ ≥ hgcdMatLen M hM i0 then lenQ else hgcdMatLen M hM i0)
      (if lenQ ≥ hgcdMatLen M hM i0 then hgcdMatPtr M hM i0 else Q)
      (if lenQ ≥ hgcdMatLen M hM i0 then hgcdMatLen M hM i0 else lenQ)
      scratch heap = .ok heap1) :
    RawHeap.SameLayout heap heap1 ∧
      RawDensePolyRep this heap1 T
        (lenQ + hgcdMatLen M hM i0 - 1) (quotient * entry0) := by
  by_cases horder : lenQ ≥ hgcdMatLen M hM i0
  · rcases mul_refines_rawDense this T Q lenQ (hgcdMatPtr M hM i0)
        (hgcdMatLen M hM i0) scratch heap quotient entry0 hcfg hp hQpos
        hEntryPos horder
        (by simpa [max_eq_left horder] using hLenWord)
        (by simpa [max_eq_left horder] using hT)
        (by simpa [max_eq_left horder] using hScratch)
        (by simpa [max_eq_left horder] using hTQ)
        (by simpa [max_eq_left horder] using hTEntry)
        (by simpa [max_eq_left horder] using hTScratch)
        (by simpa [max_eq_left horder] using hScratchQ)
        (by simpa [max_eq_left horder] using hScratchEntry)
        hQRep hEntryRep with ⟨heap', hrun, hlayout, hrep⟩
    have heq : heap' = heap1 := Except.ok.inj (hrun.symm.trans (by
      simpa [horder] using hmul))
    subst heap'
    exact ⟨hlayout, hrep⟩
  · have hle : lenQ ≤ hgcdMatLen M hM i0 := by omega
    rcases mul_refines_rawDense this T (hgcdMatPtr M hM i0)
        (hgcdMatLen M hM i0) Q lenQ scratch heap entry0 quotient hcfg hp
        hEntryPos hQpos hle
        (by simpa [max_eq_right hle] using hLenWord)
        (by simpa [max_eq_right hle] using hT)
        (by simpa [max_eq_right hle] using hScratch)
        (by simpa [max_eq_right hle] using hTEntry)
        (by simpa [max_eq_right hle] using hTQ)
        (by simpa [max_eq_right hle] using hTScratch)
        (by simpa [max_eq_right hle] using hScratchEntry)
        (by simpa [max_eq_right hle] using hScratchQ)
        hEntryRep hQRep with ⟨heap', hrun, hlayout, hrep⟩
    have heq : heap' = heap1 := Except.ok.inj (hrun.symm.trans (by
      simpa [horder] using hmul))
    subst heap'
    exact ⟨hlayout, by simpa [mul_comm, Nat.add_comm] using hrep⟩

/-- The same generated multiplication call leaves the old `i1` entry
byte-for-byte unchanged when that entry is outside both writable areas. -/
theorem matRowUpdate_mul_preserves_entry1 (this : DenseUPolyZp)
    (M : HgcdMat) (hM : M.Valid) (i0 i1 : Fin 4)
    (Q T scratch : RawPtr UInt64) (lenQ : Nat) (heap heap1 : RawHeap)
    (entry1 : Polynomial (ZMod this._p.toNat))
    (hQpos : 0 < lenQ) (hEntryPos : 0 < hgcdMatLen M hM i0)
    (hT : heap.ValidU64Slice T
      (2 * max lenQ (hgcdMatLen M hM i0) - 1))
    (hScratch : heap.ValidU64Slice scratch
      (8 * max lenQ (hgcdMatLen M hM i0)))
    (hQValid : heap.ValidU64Slice Q lenQ)
    (hEntryValid : heap.ValidU64Slice (hgcdMatPtr M hM i0)
      (hgcdMatLen M hM i0))
    (hScratchQ : U64SlicesDisjoint scratch
      (8 * max lenQ (hgcdMatLen M hM i0)) Q lenQ)
    (hScratchEntry : U64SlicesDisjoint scratch
      (8 * max lenQ (hgcdMatLen M hM i0))
      (hgcdMatPtr M hM i0) (hgcdMatLen M hM i0))
    (hTEntry1 : U64SlicesDisjoint T
      (2 * max lenQ (hgcdMatLen M hM i0) - 1)
      (hgcdMatPtr M hM i1) (hgcdMatLen M hM i1))
    (hScratchEntry1 : U64SlicesDisjoint scratch
      (8 * max lenQ (hgcdMatLen M hM i0))
      (hgcdMatPtr M hM i1) (hgcdMatLen M hM i1))
    (hEntry1Rep : RawDensePolyRep this heap (hgcdMatPtr M hM i1)
      (hgcdMatLen M hM i1) entry1)
    (hlayout : RawHeap.SameLayout heap heap1)
    (hmul : Generated.StrictMul.dense_upoly_zp__mul_ir this T
      (if lenQ ≥ hgcdMatLen M hM i0 then Q else hgcdMatPtr M hM i0)
      (if lenQ ≥ hgcdMatLen M hM i0 then lenQ else hgcdMatLen M hM i0)
      (if lenQ ≥ hgcdMatLen M hM i0 then hgcdMatPtr M hM i0 else Q)
      (if lenQ ≥ hgcdMatLen M hM i0 then hgcdMatLen M hM i0 else lenQ)
      scratch heap = .ok heap1) :
    RawDensePolyRep this heap1 (hgcdMatPtr M hM i1)
      (hgcdMatLen M hM i1) entry1 := by
  by_cases horder : lenQ ≥ hgcdMatLen M hM i0
  · have hsame := mul_preserves_prefix this T Q lenQ
      (hgcdMatPtr M hM i0) (hgcdMatLen M hM i0) scratch
      (hgcdMatPtr M hM i1) (hgcdMatLen M hM i1) heap heap1 hQpos
      hEntryPos horder
      (by simpa [max_eq_left horder] using hT) hQValid hEntryValid
      (by simpa [max_eq_left horder] using hScratch)
      (by simpa [max_eq_left horder] using hScratchEntry)
      (by simpa [max_eq_left horder] using hTEntry1)
      (by simpa [max_eq_left horder] using hScratchEntry1)
      (by simpa [horder] using hmul)
    exact rawDensePolyRep_of_same_prefix this heap heap1
      (hgcdMatPtr M hM i1) (hgcdMatLen M hM i1) entry1 hlayout hsame
      hEntry1Rep
  · have hle : lenQ ≤ hgcdMatLen M hM i0 := by omega
    have hsame := mul_preserves_prefix this T (hgcdMatPtr M hM i0)
      (hgcdMatLen M hM i0) Q lenQ scratch (hgcdMatPtr M hM i1)
      (hgcdMatLen M hM i1) heap heap1 hEntryPos hQpos hle
      (by simpa [max_eq_right hle] using hT) hEntryValid hQValid
      (by simpa [max_eq_right hle] using hScratch)
      (by simpa [max_eq_right hle] using hScratchQ)
      (by simpa [max_eq_right hle] using hTEntry1)
      (by simpa [max_eq_right hle] using hScratchEntry1)
      (by simpa [horder] using hmul)
    exact rawDensePolyRep_of_same_prefix this heap heap1
      (hgcdMatPtr M hM i1) (hgcdMatLen M hM i1) entry1 hlayout hsame
      hEntry1Rep

/-- Generic frame form of the generated multiplication inside a row update.
Any raw polynomial slice disjoint from both physical write areas (`T` and
`scratch`) is preserved. -/
theorem matRowUpdate_mul_preserves_guard (this : DenseUPolyZp)
    (M : HgcdMat) (hM : M.Valid) (i0 : Fin 4)
    (Q T scratch guard : RawPtr UInt64) (lenQ guardLen : Nat)
    (heap heap1 : RawHeap) (guardPoly : Polynomial (ZMod this._p.toNat))
    (hQpos : 0 < lenQ) (hEntryPos : 0 < hgcdMatLen M hM i0)
    (hT : heap.ValidU64Slice T
      (2 * max lenQ (hgcdMatLen M hM i0) - 1))
    (hScratch : heap.ValidU64Slice scratch
      (8 * max lenQ (hgcdMatLen M hM i0)))
    (hQValid : heap.ValidU64Slice Q lenQ)
    (hEntryValid : heap.ValidU64Slice (hgcdMatPtr M hM i0)
      (hgcdMatLen M hM i0))
    (hScratchQ : U64SlicesDisjoint scratch
      (8 * max lenQ (hgcdMatLen M hM i0)) Q lenQ)
    (hScratchEntry : U64SlicesDisjoint scratch
      (8 * max lenQ (hgcdMatLen M hM i0))
      (hgcdMatPtr M hM i0) (hgcdMatLen M hM i0))
    (hTGuard : U64SlicesDisjoint T
      (2 * max lenQ (hgcdMatLen M hM i0) - 1) guard guardLen)
    (hScratchGuard : U64SlicesDisjoint scratch
      (8 * max lenQ (hgcdMatLen M hM i0)) guard guardLen)
    (hGuardRep : RawDensePolyRep this heap guard guardLen guardPoly)
    (hlayout : RawHeap.SameLayout heap heap1)
    (hmul : Generated.StrictMul.dense_upoly_zp__mul_ir this T
      (if lenQ ≥ hgcdMatLen M hM i0 then Q else hgcdMatPtr M hM i0)
      (if lenQ ≥ hgcdMatLen M hM i0 then lenQ else hgcdMatLen M hM i0)
      (if lenQ ≥ hgcdMatLen M hM i0 then hgcdMatPtr M hM i0 else Q)
      (if lenQ ≥ hgcdMatLen M hM i0 then hgcdMatLen M hM i0 else lenQ)
      scratch heap = .ok heap1) :
    RawDensePolyRep this heap1 guard guardLen guardPoly := by
  by_cases horder : lenQ ≥ hgcdMatLen M hM i0
  · have hsame := mul_preserves_prefix this T Q lenQ
      (hgcdMatPtr M hM i0) (hgcdMatLen M hM i0) scratch guard guardLen
      heap heap1 hQpos hEntryPos horder
      (by simpa [max_eq_left horder] using hT) hQValid hEntryValid
      (by simpa [max_eq_left horder] using hScratch)
      (by simpa [max_eq_left horder] using hScratchEntry)
      (by simpa [max_eq_left horder] using hTGuard)
      (by simpa [max_eq_left horder] using hScratchGuard)
      (by simpa [horder] using hmul)
    exact rawDensePolyRep_of_same_prefix this heap heap1 guard guardLen
      guardPoly hlayout hsame hGuardRep
  · have hle : lenQ ≤ hgcdMatLen M hM i0 := by omega
    have hsame := mul_preserves_prefix this T (hgcdMatPtr M hM i0)
      (hgcdMatLen M hM i0) Q lenQ scratch guard guardLen heap heap1
      hEntryPos hQpos hle
      (by simpa [max_eq_right hle] using hT) hEntryValid hQValid
      (by simpa [max_eq_right hle] using hScratch)
      (by simpa [max_eq_right hle] using hScratchQ)
      (by simpa [max_eq_right hle] using hTGuard)
      (by simpa [max_eq_right hle] using hScratchGuard)
      (by simpa [horder] using hmul)
    exact rawDensePolyRep_of_same_prefix this heap heap1 guard guardLen
      guardPoly hlayout hsame hGuardRep

/-- The generated row addition writes only `t`; therefore the old `i0`
buffer installed as the new `i1` entry retains its normalized polynomial. -/
theorem matRowUpdate_add_preserves_entry0 (this : DenseUPolyZp)
    (t p1 T p0 : RawPtr UInt64) (l1 productLen l0 sumLen : Nat)
    (heap1 heap2 : RawHeap)
    (entry0 entry1 product : Polynomial (ZMod this._p.toNat))
    (hOutput : heap1.ValidU64Slice t (max l1 productLen))
    (hEntry1 : RawDensePolyRep this heap1 p1 l1 entry1)
    (hProduct : RawDensePolyRep this heap1 T productLen product)
    (hEntry0 : RawDensePolyRep this heap1 p0 l0 entry0)
    (htp0 : t.region ≠ p0.region)
    (hadd : Generated.StrictPolyAddSub.dense_upoly_zp__poly_add_ir this t
      p1 l1 T productLen heap1 = .ok (heap2, sumLen)) :
    RawDensePolyRep this heap2 p0 l0 entry0 := by
  rcases polyAdd_ok this t p1 l1 T productLen heap1 hOutput hEntry1.1
      hProduct.1 with ⟨heap', length, hrun, hlayout, _⟩
  have heq : (heap', length) = (heap2, sumLen) :=
    Except.ok.inj (hrun.symm.trans hadd)
  have hheap : heap' = heap2 := congrArg Prod.fst heq
  subst heap'
  have hsame := polyAdd_preserves_prefix_region_ne this t p1 l1 T p0
    productLen l0 heap1 heap2 sumLen hOutput hEntry1.1 hProduct.1 htp0
    hadd
  exact rawDensePolyRep_of_same_prefix this heap1 heap2 p0 l0 entry0
    hlayout hsame hEntry0

/-- Algebraic result of the actual add call exposed by the nonzero row-update
branch.  The product premise is supplied by strict `_mul`; this theorem then
binds the generated `_poly_add` result to the descriptor installed in M[i0]. -/
theorem matRowUpdate_nonzero_sum_rep (this : DenseUPolyZp)
    (M' : HgcdMat) (hM' : M'.Valid) (i0 : Fin 4)
    (t p1 T : RawPtr UInt64) (l1 productLen sumLen : Nat)
    (heap1 heap2 : RawHeap)
    (entry1 product : Polynomial (ZMod this._p.toNat))
    (hp : this._p ≠ 0)
    (hOutput : heap1.ValidU64Slice t (max l1 productLen))
    (hEntry1 : RawDensePolyRep this heap1 p1 l1 entry1)
    (hProduct : RawDensePolyRep this heap1 T productLen product)
    (hAliasEntry : ExactOrDisjoint t p1)
    (hAliasProduct : ExactOrDisjoint t T)
    (hadd : Generated.StrictPolyAddSub.dense_upoly_zp__poly_add_ir this t
      p1 l1 T productLen heap1 = .ok (heap2, sumLen))
    (hptr : hgcdMatPtr M' hM' i0 = t)
    (hlen : hgcdMatLen M' hM' i0 = sumLen) :
    RawDensePolyRep this heap2 (hgcdMatPtr M' hM' i0)
      (hgcdMatLen M' hM' i0) (entry1 + product) := by
  have hsum := polyAdd_refines this t p1 l1 T productLen heap1 heap2
    sumLen entry1 product hp hOutput hEntry1 hProduct hAliasEntry
    hAliasProduct hadd
  simpa [hptr, hlen] using hsum

/-- Complete semantic refinement of the successful nonzero C++ row update.
Both installed matrix entries are tied to the same generated `_mul` and
`_poly_add` executions exposed by `hrun`; no specification execution is
substituted for either call. -/
theorem matRowUpdate_nonzero_refines (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (M : HgcdMat) (i0 i1 : Fin 4) (Q : RawPtr UInt64) (lenQ : Nat)
    (T : RawPtr UInt64) (lenT : Nat) (t scratch : RawPtr UInt64)
    (heap : RawHeap) (result : MatRowUpdateResult) (hM : M.Valid)
    (quotient entry0 entry1 : Polynomial (ZMod this._p.toNat))
    (hne : i0 ≠ i1) (hQpos : 0 < lenQ)
    (hEntryPos : 0 < hgcdMatLen M hM i0)
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (hLenWord : max lenQ (hgcdMatLen M hM i0) < limbBase)
    (hT : heap.ValidU64Slice T
      (2 * max lenQ (hgcdMatLen M hM i0) - 1))
    (hScratch : heap.ValidU64Slice scratch
      (8 * max lenQ (hgcdMatLen M hM i0)))
    (hAddOutput : heap.ValidU64Slice t
      (max (hgcdMatLen M hM i1)
        (lenQ + hgcdMatLen M hM i0 - 1)))
    (hTQ : U64SlicesDisjoint T
      (2 * max lenQ (hgcdMatLen M hM i0) - 1) Q lenQ)
    (hTEntry0 : U64SlicesDisjoint T
      (2 * max lenQ (hgcdMatLen M hM i0) - 1)
      (hgcdMatPtr M hM i0) (hgcdMatLen M hM i0))
    (hTEntry1 : U64SlicesDisjoint T
      (2 * max lenQ (hgcdMatLen M hM i0) - 1)
      (hgcdMatPtr M hM i1) (hgcdMatLen M hM i1))
    (hTScratch : U64SlicesDisjoint T
      (2 * max lenQ (hgcdMatLen M hM i0) - 1) scratch
      (8 * max lenQ (hgcdMatLen M hM i0)))
    (hScratchQ : U64SlicesDisjoint scratch
      (8 * max lenQ (hgcdMatLen M hM i0)) Q lenQ)
    (hScratchEntry0 : U64SlicesDisjoint scratch
      (8 * max lenQ (hgcdMatLen M hM i0))
      (hgcdMatPtr M hM i0) (hgcdMatLen M hM i0))
    (hScratchEntry1 : U64SlicesDisjoint scratch
      (8 * max lenQ (hgcdMatLen M hM i0))
      (hgcdMatPtr M hM i1) (hgcdMatLen M hM i1))
    (hAliasEntry1 : ExactOrDisjoint t (hgcdMatPtr M hM i1))
    (hAliasProduct : ExactOrDisjoint t T)
    (htEntry0 : t.region ≠ (hgcdMatPtr M hM i0).region)
    (hQRep : RawDensePolyRep this heap Q lenQ quotient)
    (hEntry0Rep : RawDensePolyRep this heap (hgcdMatPtr M hM i0)
      (hgcdMatLen M hM i0) entry0)
    (hEntry1Rep : RawDensePolyRep this heap (hgcdMatPtr M hM i1)
      (hgcdMatLen M hM i1) entry1)
    (hrun : dense_upoly_zp__mat_row_update_ir this M i0 i1 Q lenQ T
      lenT t scratch heap = .ok result) :
    ∃ hResult : result.matrix.Valid,
      RawDensePolyRep this result.heap
        (hgcdMatPtr result.matrix hResult i0)
        (hgcdMatLen result.matrix hResult i0)
        (entry1 + quotient * entry0) ∧
      RawDensePolyRep this result.heap
        (hgcdMatPtr result.matrix hResult i1)
        (hgcdMatLen result.matrix hResult i1) entry0 := by
  rcases matRowUpdate_nonzero_success_shape this M i0 i1 Q lenQ T lenT t
      scratch heap result hM hne (by omega) (by omega) hrun with
    ⟨heap1, heap2, sumLen, hmul, hadd, hResultHeap, _, _, _,
      hResult, hptr0, hlen0, hptr1, hlen1, _⟩
  have hmulSemantic := matRowUpdate_mul_result this M hM i0 Q T scratch
    lenQ heap heap1 quotient entry0 hcfg hp hQpos hEntryPos hLenWord hT
    hScratch hTQ hTEntry0 hTScratch hScratchQ hScratchEntry0 hQRep
    hEntry0Rep hmul
  rcases hmulSemantic with ⟨hlayoutMul, hProduct1⟩
  have hEntry1_1 := matRowUpdate_mul_preserves_entry1 this M hM i0 i1 Q T
    scratch lenQ heap heap1 entry1 hQpos hEntryPos hT hScratch hQRep.1
    hEntry0Rep.1 hScratchQ hScratchEntry0 hTEntry1 hScratchEntry1
    hEntry1Rep hlayoutMul hmul
  have hEntry0_1 := matRowUpdate_mul_preserves_entry1 this M hM i0 i0 Q T
    scratch lenQ heap heap1 entry0 hQpos hEntryPos hT hScratch hQRep.1
    hEntry0Rep.1 hScratchQ hScratchEntry0 hTEntry0 hScratchEntry0
    hEntry0Rep hlayoutMul hmul
  have hAddOutput1 := (hlayoutMul t
    (max (hgcdMatLen M hM i1)
      (lenQ + hgcdMatLen M hM i0 - 1))).mp hAddOutput
  have hpWord : this._p ≠ 0 := by
    intro hzero
    have := congrArg UInt64.toNat hzero
    simp at this
    omega
  have hNew0 := matRowUpdate_nonzero_sum_rep this result.matrix hResult i0
    t (hgcdMatPtr M hM i1) T (hgcdMatLen M hM i1)
    (lenQ + hgcdMatLen M hM i0 - 1) sumLen heap1 heap2 entry1
    (quotient * entry0) hpWord hAddOutput1 hEntry1_1 hProduct1
    hAliasEntry1 hAliasProduct hadd hptr0 hlen0
  have hOld0_2 := matRowUpdate_add_preserves_entry0 this t
    (hgcdMatPtr M hM i1) T (hgcdMatPtr M hM i0)
    (hgcdMatLen M hM i1) (lenQ + hgcdMatLen M hM i0 - 1)
    (hgcdMatLen M hM i0) sumLen heap1 heap2 entry0 entry1
    (quotient * entry0) hAddOutput1 hEntry1_1 hProduct1 hEntry0_1
    htEntry0 hadd
  refine ⟨hResult, ?_, ?_⟩
  · simpa [hResultHeap] using hNew0
  · simpa [hResultHeap, hptr1, hlen1] using hOld0_2

/-- Branch-complete semantic refinement of the generated C++ row update.
The same conclusion is obtained from either the physical swap-only branch or
the real `_mul`/`_poly_add` branch. -/
theorem matRowUpdate_refines (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (M : HgcdMat) (i0 i1 : Fin 4) (Q : RawPtr UInt64) (lenQ : Nat)
    (T : RawPtr UInt64) (lenT : Nat) (t scratch : RawPtr UInt64)
    (heap : RawHeap) (result : MatRowUpdateResult) (hM : M.Valid)
    (quotient entry0 entry1 : Polynomial (ZMod this._p.toNat))
    (hne : i0 ≠ i1)
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (hLenWord : max lenQ (hgcdMatLen M hM i0) < limbBase)
    (hT : heap.ValidU64Slice T
      (2 * max lenQ (hgcdMatLen M hM i0) - 1))
    (hScratch : heap.ValidU64Slice scratch
      (8 * max lenQ (hgcdMatLen M hM i0)))
    (hAddOutput : heap.ValidU64Slice t
      (max (hgcdMatLen M hM i1)
        (lenQ + hgcdMatLen M hM i0 - 1)))
    (hTQ : U64SlicesDisjoint T
      (2 * max lenQ (hgcdMatLen M hM i0) - 1) Q lenQ)
    (hTEntry0 : U64SlicesDisjoint T
      (2 * max lenQ (hgcdMatLen M hM i0) - 1)
      (hgcdMatPtr M hM i0) (hgcdMatLen M hM i0))
    (hTEntry1 : U64SlicesDisjoint T
      (2 * max lenQ (hgcdMatLen M hM i0) - 1)
      (hgcdMatPtr M hM i1) (hgcdMatLen M hM i1))
    (hTScratch : U64SlicesDisjoint T
      (2 * max lenQ (hgcdMatLen M hM i0) - 1) scratch
      (8 * max lenQ (hgcdMatLen M hM i0)))
    (hScratchQ : U64SlicesDisjoint scratch
      (8 * max lenQ (hgcdMatLen M hM i0)) Q lenQ)
    (hScratchEntry0 : U64SlicesDisjoint scratch
      (8 * max lenQ (hgcdMatLen M hM i0))
      (hgcdMatPtr M hM i0) (hgcdMatLen M hM i0))
    (hScratchEntry1 : U64SlicesDisjoint scratch
      (8 * max lenQ (hgcdMatLen M hM i0))
      (hgcdMatPtr M hM i1) (hgcdMatLen M hM i1))
    (hAliasEntry1 : ExactOrDisjoint t (hgcdMatPtr M hM i1))
    (hAliasProduct : ExactOrDisjoint t T)
    (htEntry0 : t.region ≠ (hgcdMatPtr M hM i0).region)
    (hQRep : RawDensePolyRep this heap Q lenQ quotient)
    (hEntry0Rep : RawDensePolyRep this heap (hgcdMatPtr M hM i0)
      (hgcdMatLen M hM i0) entry0)
    (hEntry1Rep : RawDensePolyRep this heap (hgcdMatPtr M hM i1)
      (hgcdMatLen M hM i1) entry1)
    (hrun : dense_upoly_zp__mat_row_update_ir this M i0 i1 Q lenQ T
      lenT t scratch heap = .ok result) :
    ∃ hResult : result.matrix.Valid,
      RawDensePolyRep this result.heap
        (hgcdMatPtr result.matrix hResult i0)
        (hgcdMatLen result.matrix hResult i0)
        (entry1 + quotient * entry0) ∧
      RawDensePolyRep this result.heap
        (hgcdMatPtr result.matrix hResult i1)
        (hgcdMatLen result.matrix hResult i1) entry0 := by
  by_cases hzero : lenQ = 0 ∨ hgcdMatLen M hM i0 = 0
  · exact matRowUpdate_zero_refines this M i0 i1 Q lenQ T lenT t scratch
      heap result hM quotient entry0 entry1 hne hzero hQRep hEntry0Rep
      hEntry1Rep hrun
  · have hQpos : 0 < lenQ := Nat.pos_of_ne_zero (by
      intro heq
      exact hzero (Or.inl heq))
    have hEntryPos : 0 < hgcdMatLen M hM i0 := Nat.pos_of_ne_zero (by
      intro heq
      exact hzero (Or.inr heq))
    exact matRowUpdate_nonzero_refines this M i0 i1 Q lenQ T lenT t
      scratch heap result hM quotient entry0 entry1 hne hQpos hEntryPos
      hcfg hp hLenWord hT hScratch hAddOutput hTQ hTEntry0 hTEntry1
      hTScratch hScratchQ hScratchEntry0 hScratchEntry1 hAliasEntry1
      hAliasProduct htEntry0 hQRep hEntry0Rep hEntry1Rep hrun

/-- Branch-complete frame theorem for a generated row update.  A normalized
raw polynomial disjoint from `T`, `scratch`, and the add destination `t`
survives the exact `_mul`/`_poly_add` execution (or the swap-only branch). -/
theorem matRowUpdate_preserves_guard (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (M : HgcdMat) (i0 i1 : Fin 4) (Q : RawPtr UInt64) (lenQ : Nat)
    (T : RawPtr UInt64) (lenT : Nat) (t scratch guard : RawPtr UInt64)
    (guardLen : Nat) (heap : RawHeap) (result : MatRowUpdateResult)
    (hM : M.Valid)
    (quotient entry0 entry1 guardPoly : Polynomial (ZMod this._p.toNat))
    (hne : i0 ≠ i1)
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (hLenWord : max lenQ (hgcdMatLen M hM i0) < limbBase)
    (hT : heap.ValidU64Slice T
      (2 * max lenQ (hgcdMatLen M hM i0) - 1))
    (hScratch : heap.ValidU64Slice scratch
      (8 * max lenQ (hgcdMatLen M hM i0)))
    (hAddOutput : heap.ValidU64Slice t
      (max (hgcdMatLen M hM i1)
        (lenQ + hgcdMatLen M hM i0 - 1)))
    (hTQ : U64SlicesDisjoint T
      (2 * max lenQ (hgcdMatLen M hM i0) - 1) Q lenQ)
    (hTEntry0 : U64SlicesDisjoint T
      (2 * max lenQ (hgcdMatLen M hM i0) - 1)
      (hgcdMatPtr M hM i0) (hgcdMatLen M hM i0))
    (hTEntry1 : U64SlicesDisjoint T
      (2 * max lenQ (hgcdMatLen M hM i0) - 1)
      (hgcdMatPtr M hM i1) (hgcdMatLen M hM i1))
    (hTScratch : U64SlicesDisjoint T
      (2 * max lenQ (hgcdMatLen M hM i0) - 1) scratch
      (8 * max lenQ (hgcdMatLen M hM i0)))
    (hScratchQ : U64SlicesDisjoint scratch
      (8 * max lenQ (hgcdMatLen M hM i0)) Q lenQ)
    (hScratchEntry0 : U64SlicesDisjoint scratch
      (8 * max lenQ (hgcdMatLen M hM i0))
      (hgcdMatPtr M hM i0) (hgcdMatLen M hM i0))
    (hScratchEntry1 : U64SlicesDisjoint scratch
      (8 * max lenQ (hgcdMatLen M hM i0))
      (hgcdMatPtr M hM i1) (hgcdMatLen M hM i1))
    (hTGuard : U64SlicesDisjoint T
      (2 * max lenQ (hgcdMatLen M hM i0) - 1) guard guardLen)
    (hScratchGuard : U64SlicesDisjoint scratch
      (8 * max lenQ (hgcdMatLen M hM i0)) guard guardLen)
    (htGuard : t.region ≠ guard.region)
    (hQRep : RawDensePolyRep this heap Q lenQ quotient)
    (hEntry0Rep : RawDensePolyRep this heap (hgcdMatPtr M hM i0)
      (hgcdMatLen M hM i0) entry0)
    (hEntry1Rep : RawDensePolyRep this heap (hgcdMatPtr M hM i1)
      (hgcdMatLen M hM i1) entry1)
    (hGuardRep : RawDensePolyRep this heap guard guardLen guardPoly)
    (hrun : dense_upoly_zp__mat_row_update_ir this M i0 i1 Q lenQ T
      lenT t scratch heap = .ok result) :
    RawDensePolyRep this result.heap guard guardLen guardPoly := by
  by_cases hzero : lenQ = 0 ∨ hgcdMatLen M hM i0 = 0
  · rcases matRowUpdate_zero_exec this M i0 i1 Q lenQ T lenT t scratch
      heap hM hne hzero with ⟨actual, hactual, hheap, _⟩
    have heq : actual = result := Except.ok.inj (hactual.symm.trans hrun)
    subst actual
    simpa [hheap] using hGuardRep
  · have hQpos : 0 < lenQ := Nat.pos_of_ne_zero (by
      intro heq
      exact hzero (Or.inl heq))
    have hEntryPos : 0 < hgcdMatLen M hM i0 := Nat.pos_of_ne_zero (by
      intro heq
      exact hzero (Or.inr heq))
    rcases matRowUpdate_nonzero_success_shape this M i0 i1 Q lenQ T lenT
        t scratch heap result hM hne (by omega) (by omega) hrun with
      ⟨heap1, heap2, sumLen, hmul, hadd, hResultHeap, _, _, _, _, _, _, _,
        _, _⟩
    rcases matRowUpdate_mul_result this M hM i0 Q T scratch lenQ heap
        heap1 quotient entry0 hcfg hp hQpos hEntryPos hLenWord hT hScratch
        hTQ hTEntry0 hTScratch hScratchQ hScratchEntry0 hQRep hEntry0Rep
        hmul with ⟨hlayoutMul, hProduct1⟩
    have hGuard1 := matRowUpdate_mul_preserves_guard this M hM i0 Q T
      scratch guard lenQ guardLen heap heap1 guardPoly hQpos hEntryPos hT
      hScratch hQRep.1 hEntry0Rep.1 hScratchQ hScratchEntry0 hTGuard
      hScratchGuard hGuardRep hlayoutMul hmul
    have hEntry1_1 := matRowUpdate_mul_preserves_guard this M hM i0 Q T
      scratch (hgcdMatPtr M hM i1) lenQ (hgcdMatLen M hM i1) heap heap1
      entry1 hQpos hEntryPos hT hScratch hQRep.1 hEntry0Rep.1 hScratchQ
      hScratchEntry0 hTEntry1 hScratchEntry1
      hEntry1Rep hlayoutMul hmul
    have hAddOutput1 := (hlayoutMul t
      (max (hgcdMatLen M hM i1)
        (lenQ + hgcdMatLen M hM i0 - 1))).mp hAddOutput
    have hGuard2 := matRowUpdate_add_preserves_entry0 this t
      (hgcdMatPtr M hM i1) T guard (hgcdMatLen M hM i1)
      (lenQ + hgcdMatLen M hM i0 - 1) guardLen sumLen heap1 heap2
      guardPoly entry1 (quotient * entry0) hAddOutput1 hEntry1_1
      hProduct1 hGuard1 htGuard hadd
    simpa [hResultHeap] using hGuard2

/-- Physical L1 workspace invariant required by one generated
`_mat_row_update` call.  It records only allocation/aliasing facts about the
actual quotient, matrix entries, multiplication buffer, scratch buffer, and
addition destination; no L2 polynomial result is assumed. -/
structure MatRowUpdateWorkspace (M : HgcdMat) (i0 i1 : Fin 4)
    (Q : RawPtr UInt64) (lenQ : Nat) (T : RawPtr UInt64)
    (t scratch : RawPtr UInt64) (heap : RawHeap) (hM : M.Valid) : Prop where
  lenWord : max lenQ (hgcdMatLen M hM i0) < limbBase
  validT : heap.ValidU64Slice T
    (2 * max lenQ (hgcdMatLen M hM i0) - 1)
  validScratch : heap.ValidU64Slice scratch
    (8 * max lenQ (hgcdMatLen M hM i0))
  validAddOutput : heap.ValidU64Slice t
    (max (hgcdMatLen M hM i1)
      (lenQ + hgcdMatLen M hM i0 - 1))
  disjointTQ : U64SlicesDisjoint T
    (2 * max lenQ (hgcdMatLen M hM i0) - 1) Q lenQ
  disjointTMatrix : ∀ i : Fin 4, U64SlicesDisjoint T
    (2 * max lenQ (hgcdMatLen M hM i0) - 1)
    (hgcdMatPtr M hM i) (hgcdMatLen M hM i)
  disjointTScratch : U64SlicesDisjoint T
    (2 * max lenQ (hgcdMatLen M hM i0) - 1) scratch
    (8 * max lenQ (hgcdMatLen M hM i0))
  disjointScratchQ : U64SlicesDisjoint scratch
    (8 * max lenQ (hgcdMatLen M hM i0)) Q lenQ
  disjointScratchMatrix : ∀ i : Fin 4, U64SlicesDisjoint scratch
    (8 * max lenQ (hgcdMatLen M hM i0))
    (hgcdMatPtr M hM i) (hgcdMatLen M hM i)
  aliasEntry1 : ExactOrDisjoint t (hgcdMatPtr M hM i1)
  aliasProduct : ExactOrDisjoint t T
  tDisjointQ : t.region ≠ Q.region
  tDisjointOther : ∀ i : Fin 4, i ≠ i1 →
    t.region ≠ (hgcdMatPtr M hM i).region

/-- Extra physical separation facts for one live polynomial outside the
matrix.  These facts let the exact multiplication/addition execution frame
that polynomial without assigning it any specification-level behavior. -/
structure MatRowUpdateGuardWorkspace (M : HgcdMat) (i0 : Fin 4)
    (Q : RawPtr UInt64) (lenQ : Nat) (T t scratch guard : RawPtr UInt64)
    (guardLen : Nat) (hM : M.Valid) : Prop where
  disjointT : U64SlicesDisjoint T
    (2 * max lenQ (hgcdMatLen M hM i0) - 1) guard guardLen
  disjointScratch : U64SlicesDisjoint scratch
    (8 * max lenQ (hgcdMatLen M hM i0)) guard guardLen
  disjointAddOutput : t.region ≠ guard.region

/-- Apply the branch-complete semantic row-update theorem from the packaged
L1 workspace invariant. -/
theorem matRowUpdate_refines_of_workspace (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (M : HgcdMat) (i0 i1 : Fin 4) (Q : RawPtr UInt64) (lenQ : Nat)
    (T : RawPtr UInt64) (lenT : Nat) (t scratch : RawPtr UInt64)
    (heap : RawHeap) (result : MatRowUpdateResult) (hM : M.Valid)
    (quotient entry0 entry1 : Polynomial (ZMod this._p.toNat))
    (hne : i0 ≠ i1) (hcfg : DensePreinvConfigured this)
    (hp : 1 < this._p.toNat)
    (workspace : MatRowUpdateWorkspace M i0 i1 Q lenQ T t scratch heap hM)
    (hQRep : RawDensePolyRep this heap Q lenQ quotient)
    (hEntry0Rep : RawDensePolyRep this heap (hgcdMatPtr M hM i0)
      (hgcdMatLen M hM i0) entry0)
    (hEntry1Rep : RawDensePolyRep this heap (hgcdMatPtr M hM i1)
      (hgcdMatLen M hM i1) entry1)
    (hrun : dense_upoly_zp__mat_row_update_ir this M i0 i1 Q lenQ T
      lenT t scratch heap = .ok result) :
    ∃ hResult : result.matrix.Valid,
      RawDensePolyRep this result.heap
        (hgcdMatPtr result.matrix hResult i0)
        (hgcdMatLen result.matrix hResult i0)
        (entry1 + quotient * entry0) ∧
      RawDensePolyRep this result.heap
        (hgcdMatPtr result.matrix hResult i1)
        (hgcdMatLen result.matrix hResult i1) entry0 := by
  exact matRowUpdate_refines this M i0 i1 Q lenQ T lenT t scratch heap
    result hM quotient entry0 entry1 hne hcfg hp workspace.lenWord
    workspace.validT workspace.validScratch workspace.validAddOutput
    workspace.disjointTQ (workspace.disjointTMatrix i0)
    (workspace.disjointTMatrix i1) workspace.disjointTScratch
    workspace.disjointScratchQ (workspace.disjointScratchMatrix i0)
    (workspace.disjointScratchMatrix i1) workspace.aliasEntry1
    workspace.aliasProduct (workspace.tDisjointOther i0 hne)
    hQRep hEntry0Rep hEntry1Rep hrun

/-- Descriptor-length effect of the same real row update.  The active branch
uses the actual `_mul` heap transition and the actual `_poly_add`
normalization bound; the inactive branch is the source's pure swap. -/
theorem matRowUpdate_length_bound_of_workspace (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (M : HgcdMat) (i0 i1 : Fin 4) (Q : RawPtr UInt64) (lenQ : Nat)
    (T : RawPtr UInt64) (lenT : Nat) (t scratch : RawPtr UInt64)
    (heap : RawHeap) (result : MatRowUpdateResult) (hM : M.Valid)
    (quotient entry0 entry1 : Polynomial (ZMod this._p.toNat))
    (hne : i0 ≠ i1) (hcfg : DensePreinvConfigured this)
    (hp : 1 < this._p.toNat)
    (workspace : MatRowUpdateWorkspace M i0 i1 Q lenQ T t scratch heap hM)
    (hQRep : RawDensePolyRep this heap Q lenQ quotient)
    (hEntry0Rep : RawDensePolyRep this heap (hgcdMatPtr M hM i0)
      (hgcdMatLen M hM i0) entry0)
    (hEntry1Rep : RawDensePolyRep this heap (hgcdMatPtr M hM i1)
      (hgcdMatLen M hM i1) entry1)
    (hrun : dense_upoly_zp__mat_row_update_ir this M i0 i1 Q lenQ T
      lenT t scratch heap = .ok result) :
    ∃ hResult : result.matrix.Valid,
      hgcdMatLen result.matrix hResult i0 ≤
        max (hgcdMatLen M hM i1)
          (lenQ + hgcdMatLen M hM i0 - 1) ∧
      hgcdMatLen result.matrix hResult i1 = hgcdMatLen M hM i0 ∧
      ∀ j : Fin 4, j ≠ i0 → j ≠ i1 →
        hgcdMatLen result.matrix hResult j = hgcdMatLen M hM j := by
  by_cases hzero : lenQ = 0 ∨ hgcdMatLen M hM i0 = 0
  · rcases matRowUpdate_zero_exec this M i0 i1 Q lenQ T lenT t scratch
      heap hM hne hzero with
    ⟨actual, hactual, _, _, _, _, hActual, _, hlen0, _, hlen1, hother⟩
    have heq : actual = result := Except.ok.inj (hactual.symm.trans hrun)
    subst actual
    exact ⟨hActual, hlen0.le.trans (le_max_left _ _), hlen1,
      fun j hj0 hj1 => (hother j hj0 hj1).2⟩
  · have hQPos : 0 < lenQ := by omega
    have hEntryPos : 0 < hgcdMatLen M hM i0 := by omega
    rcases matRowUpdate_nonzero_success_shape this M i0 i1 Q lenQ T lenT
      t scratch heap result hM hne (by omega) (by omega) hrun with
      ⟨heap1, heap2, sumLen, hmul, hadd, _, _, _, _, hResult,
        _, hlen0, _, hlen1, hother⟩
    have hProduct := matRowUpdate_mul_result this M hM i0 Q T scratch lenQ
      heap heap1 quotient entry0 hcfg hp hQPos hEntryPos workspace.lenWord
      workspace.validT workspace.validScratch workspace.disjointTQ
      (workspace.disjointTMatrix i0) workspace.disjointTScratch
      workspace.disjointScratchQ (workspace.disjointScratchMatrix i0)
      hQRep hEntry0Rep hmul
    have hOutput1 := (hProduct.1 _ _).mp workspace.validAddOutput
    have hEntry1Valid1 :=
      (hProduct.1 (hgcdMatPtr M hM i1) (hgcdMatLen M hM i1)).mp
        hEntry1Rep.1
    rcases polyAdd_ok this t (hgcdMatPtr M hM i1)
        (hgcdMatLen M hM i1) T
        (lenQ + hgcdMatLen M hM i0 - 1) heap1 hOutput1
        hEntry1Valid1 hProduct.2.1 with
      ⟨heapAdd, lenAdd, haddAdd, _, hsum⟩
    have heqAdd : (heapAdd, lenAdd) = (heap2, sumLen) :=
      Except.ok.inj (haddAdd.symm.trans hadd)
    cases heqAdd
    exact ⟨hResult, hlen0.le.trans hsum, hlen1,
      fun j hj0 hj1 => (hother j hj0 hj1).2⟩

/-- Any non-target raw polynomial covered by the same physical workspace
survives the exact generated row-update execution. -/
theorem matRowUpdate_preserves_matrix_entry_of_workspace
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (M : HgcdMat) (i0 i1 guardIndex : Fin 4)
    (Q : RawPtr UInt64) (lenQ : Nat) (T : RawPtr UInt64) (lenT : Nat)
    (t scratch : RawPtr UInt64) (heap : RawHeap)
    (result : MatRowUpdateResult) (hM : M.Valid)
    (quotient entry0 entry1 guardPoly : Polynomial (ZMod this._p.toNat))
    (hne : i0 ≠ i1) (hguard : guardIndex ≠ i1)
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (workspace : MatRowUpdateWorkspace M i0 i1 Q lenQ T t scratch heap hM)
    (hQRep : RawDensePolyRep this heap Q lenQ quotient)
    (hEntry0Rep : RawDensePolyRep this heap (hgcdMatPtr M hM i0)
      (hgcdMatLen M hM i0) entry0)
    (hEntry1Rep : RawDensePolyRep this heap (hgcdMatPtr M hM i1)
      (hgcdMatLen M hM i1) entry1)
    (hGuardRep : RawDensePolyRep this heap (hgcdMatPtr M hM guardIndex)
      (hgcdMatLen M hM guardIndex) guardPoly)
    (hrun : dense_upoly_zp__mat_row_update_ir this M i0 i1 Q lenQ T
      lenT t scratch heap = .ok result) :
    RawDensePolyRep this result.heap (hgcdMatPtr M hM guardIndex)
      (hgcdMatLen M hM guardIndex) guardPoly := by
  exact matRowUpdate_preserves_guard this M i0 i1 Q lenQ T lenT t scratch
    (hgcdMatPtr M hM guardIndex) (hgcdMatLen M hM guardIndex) heap result
    hM quotient entry0 entry1 guardPoly hne hcfg hp workspace.lenWord
    workspace.validT workspace.validScratch workspace.validAddOutput
    workspace.disjointTQ (workspace.disjointTMatrix i0)
    (workspace.disjointTMatrix i1) workspace.disjointTScratch
    workspace.disjointScratchQ (workspace.disjointScratchMatrix i0)
    (workspace.disjointScratchMatrix i1)
    (workspace.disjointTMatrix guardIndex)
    (workspace.disjointScratchMatrix guardIndex)
    (workspace.tDisjointOther guardIndex hguard) hQRep hEntry0Rep
    hEntry1Rep hGuardRep hrun

/-- The quotient itself survives one generated row update, which is needed
because the source reuses the same quotient for the second matrix row. -/
theorem matRowUpdate_preserves_quotient_of_workspace
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (M : HgcdMat) (i0 i1 : Fin 4) (Q : RawPtr UInt64) (lenQ : Nat)
    (T : RawPtr UInt64) (lenT : Nat) (t scratch : RawPtr UInt64)
    (heap : RawHeap) (result : MatRowUpdateResult) (hM : M.Valid)
    (quotient entry0 entry1 : Polynomial (ZMod this._p.toNat))
    (hne : i0 ≠ i1) (hcfg : DensePreinvConfigured this)
    (hp : 1 < this._p.toNat)
    (workspace : MatRowUpdateWorkspace M i0 i1 Q lenQ T t scratch heap hM)
    (hQRep : RawDensePolyRep this heap Q lenQ quotient)
    (hEntry0Rep : RawDensePolyRep this heap (hgcdMatPtr M hM i0)
      (hgcdMatLen M hM i0) entry0)
    (hEntry1Rep : RawDensePolyRep this heap (hgcdMatPtr M hM i1)
      (hgcdMatLen M hM i1) entry1)
    (hrun : dense_upoly_zp__mat_row_update_ir this M i0 i1 Q lenQ T
      lenT t scratch heap = .ok result) :
    RawDensePolyRep this result.heap Q lenQ quotient := by
  exact matRowUpdate_preserves_guard this M i0 i1 Q lenQ T lenT t scratch
    Q lenQ heap result hM quotient entry0 entry1 quotient hne hcfg hp
    workspace.lenWord workspace.validT workspace.validScratch
    workspace.validAddOutput workspace.disjointTQ
    (workspace.disjointTMatrix i0) (workspace.disjointTMatrix i1)
    workspace.disjointTScratch workspace.disjointScratchQ
    (workspace.disjointScratchMatrix i0)
    (workspace.disjointScratchMatrix i1) workspace.disjointTQ
    workspace.disjointScratchQ workspace.tDisjointQ hQRep hEntry0Rep
    hEntry1Rep hQRep hrun

/-- Apply the general frame theorem using the packaged row workspace and the
three additional separation facts for an external live polynomial. -/
theorem matRowUpdate_preserves_guard_of_workspaces
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (M : HgcdMat) (i0 i1 : Fin 4) (Q : RawPtr UInt64) (lenQ : Nat)
    (T : RawPtr UInt64) (lenT : Nat) (t scratch guard : RawPtr UInt64)
    (guardLen : Nat) (heap : RawHeap) (result : MatRowUpdateResult)
    (hM : M.Valid)
    (quotient entry0 entry1 guardPoly : Polynomial (ZMod this._p.toNat))
    (hne : i0 ≠ i1) (hcfg : DensePreinvConfigured this)
    (hp : 1 < this._p.toNat)
    (workspace : MatRowUpdateWorkspace M i0 i1 Q lenQ T t scratch heap hM)
    (guardWorkspace : MatRowUpdateGuardWorkspace M i0 Q lenQ T t scratch
      guard guardLen hM)
    (hQRep : RawDensePolyRep this heap Q lenQ quotient)
    (hEntry0Rep : RawDensePolyRep this heap (hgcdMatPtr M hM i0)
      (hgcdMatLen M hM i0) entry0)
    (hEntry1Rep : RawDensePolyRep this heap (hgcdMatPtr M hM i1)
      (hgcdMatLen M hM i1) entry1)
    (hGuardRep : RawDensePolyRep this heap guard guardLen guardPoly)
    (hrun : dense_upoly_zp__mat_row_update_ir this M i0 i1 Q lenQ T
      lenT t scratch heap = .ok result) :
    RawDensePolyRep this result.heap guard guardLen guardPoly := by
  exact matRowUpdate_preserves_guard this M i0 i1 Q lenQ T lenT t scratch
    guard guardLen heap result hM quotient entry0 entry1 guardPoly hne hcfg
    hp workspace.lenWord workspace.validT workspace.validScratch
    workspace.validAddOutput workspace.disjointTQ
    (workspace.disjointTMatrix i0) (workspace.disjointTMatrix i1)
    workspace.disjointTScratch workspace.disjointScratchQ
    (workspace.disjointScratchMatrix i0)
    (workspace.disjointScratchMatrix i1) guardWorkspace.disjointT
    guardWorkspace.disjointScratch guardWorkspace.disjointAddOutput hQRep
    hEntry0Rep hEntry1Rep hGuardRep hrun

/-- The two actual generated `_mat_row_update` calls refine the complete
four-entry mathematical HGCD matrix step.  In particular, this theorem
executes and frames both physical calls; it does not replace them with the
L2 matrix formula. -/
theorem hgcdTwoRowUpdates_refine_matrix (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (M : HgcdMat) (Q : RawPtr UInt64) (lenQ : Nat)
    (T : RawPtr UInt64) (lenT : Nat) (t scratch : RawPtr UInt64)
    (heap : RawHeap) (row23 row01 : MatRowUpdateResult) (hM : M.Valid)
    (h23 : row23.matrix.Valid)
    (quotient : Polynomial (ZMod this._p.toNat))
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (workspace23 : MatRowUpdateWorkspace M (2 : Fin 4) (3 : Fin 4)
      Q lenQ T t scratch heap hM)
    (workspace01 : MatRowUpdateWorkspace row23.matrix (0 : Fin 4)
      (1 : Fin 4) Q lenQ row23.T row23.t scratch row23.heap h23)
    (hQRep : RawDensePolyRep this heap Q lenQ quotient)
    (hMatrix : HgcdMatRawDenseRep this heap M entries hM)
    (hrow23 : dense_upoly_zp__mat_row_update_ir this M
      (2 : Fin 4) (3 : Fin 4) Q lenQ T lenT t scratch heap = .ok row23)
    (hrow01 : dense_upoly_zp__mat_row_update_ir this row23.matrix
      (0 : Fin 4) (1 : Fin 4) Q lenQ row23.T row23.lenT row23.t
      scratch row23.heap = .ok row01) :
    ∃ h01 : row01.matrix.Valid,
      HgcdMatRawDenseRep this row01.heap row01.matrix
        (CLPoly.Impl.StrictHGCDRefinement.hgcdStepEntries quotient entries)
        h01 := by
  have hE0 := hMatrix (0 : Fin 4)
  have hE1 := hMatrix (1 : Fin 4)
  have hE2 := hMatrix (2 : Fin 4)
  have hE3 := hMatrix (3 : Fin 4)
  rcases matRowUpdate_refines_of_workspace this M (2 : Fin 4)
      (3 : Fin 4) Q lenQ T lenT t scratch heap row23 hM quotient
      (entries 2) (entries 3) (by decide) hcfg hp workspace23 hQRep hE2
      hE3 hrow23 with ⟨h23', hNew2', hNew3'⟩
  have hh23 : h23' = h23 := Subsingleton.elim _ _
  subst h23'
  have hQ23 := matRowUpdate_preserves_quotient_of_workspace this M
    (2 : Fin 4) (3 : Fin 4) Q lenQ T lenT t scratch heap row23 hM
    quotient (entries 2) (entries 3) (by decide) hcfg hp workspace23
    hQRep hE2 hE3 hrow23
  have hE0Old23 := matRowUpdate_preserves_matrix_entry_of_workspace this M
    (2 : Fin 4) (3 : Fin 4) (0 : Fin 4) Q lenQ T lenT t scratch heap
    row23 hM quotient (entries 2) (entries 3) (entries 0) (by decide)
    (by decide) hcfg hp workspace23 hQRep hE2 hE3 hE0 hrow23
  have hE1Old23 := matRowUpdate_preserves_matrix_entry_of_workspace this M
    (2 : Fin 4) (3 : Fin 4) (1 : Fin 4) Q lenQ T lenT t scratch heap
    row23 hM quotient (entries 2) (entries 3) (entries 1) (by decide)
    (by decide) hcfg hp workspace23 hQRep hE2 hE3 hE1 hrow23
  rcases hgcdTwoRowUpdates_descriptor_frame this M Q lenQ T lenT t scratch
      heap row23 row01 hM hrow23 hrow01 with
    ⟨h23Frame, h01Frame, hframe23, hframe01⟩
  have hh23Frame : h23Frame = h23 := Subsingleton.elim _ _
  subst h23Frame
  have hE0_23 : RawDensePolyRep this row23.heap
      (hgcdMatPtr row23.matrix h23 (0 : Fin 4))
      (hgcdMatLen row23.matrix h23 (0 : Fin 4)) (entries 0) := by
    have hptr : hgcdMatPtr row23.matrix h23 (0 : Fin 4) =
        hgcdMatPtr M hM (0 : Fin 4) := by
      simpa only using (hframe23 (0 : Fin 2)).1
    have hlen : hgcdMatLen row23.matrix h23 (0 : Fin 4) =
        hgcdMatLen M hM (0 : Fin 4) := by
      simpa only using (hframe23 (0 : Fin 2)).2
    rw [hptr, hlen]
    exact hE0Old23
  have hE1_23 : RawDensePolyRep this row23.heap
      (hgcdMatPtr row23.matrix h23 (1 : Fin 4))
      (hgcdMatLen row23.matrix h23 (1 : Fin 4)) (entries 1) := by
    have hptr : hgcdMatPtr row23.matrix h23 (1 : Fin 4) =
        hgcdMatPtr M hM (1 : Fin 4) := by
      simpa only using (hframe23 (1 : Fin 2)).1
    have hlen : hgcdMatLen row23.matrix h23 (1 : Fin 4) =
        hgcdMatLen M hM (1 : Fin 4) := by
      simpa only using (hframe23 (1 : Fin 2)).2
    rw [hptr, hlen]
    exact hE1Old23
  rcases matRowUpdate_refines_of_workspace this row23.matrix (0 : Fin 4)
      (1 : Fin 4) Q lenQ row23.T row23.lenT row23.t scratch row23.heap
      row01 h23 quotient (entries 0) (entries 1) (by decide) hcfg hp
      workspace01 hQ23 hE0_23 hE1_23 hrow01 with
    ⟨h01, hNew0, hNew1⟩
  have hNew2 := matRowUpdate_preserves_matrix_entry_of_workspace this
    row23.matrix (0 : Fin 4) (1 : Fin 4) (2 : Fin 4) Q lenQ row23.T
    row23.lenT row23.t scratch row23.heap row01 h23 quotient
    (entries 0) (entries 1) (entries 3 + quotient * entries 2)
    (by decide) (by decide) hcfg hp workspace01 hQ23 hE0_23 hE1_23
    hNew2' hrow01
  have hNew3 := matRowUpdate_preserves_matrix_entry_of_workspace this
    row23.matrix (0 : Fin 4) (1 : Fin 4) (3 : Fin 4) Q lenQ row23.T
    row23.lenT row23.t scratch row23.heap row01 h23 quotient
    (entries 0) (entries 1) (entries 2) (by decide) (by decide) hcfg hp
    workspace01 hQ23 hE0_23 hE1_23 hNew3' hrow01
  have hh01Frame : h01Frame = h01 := Subsingleton.elim _ _
  subst h01Frame
  refine ⟨h01, ?_⟩
  intro i
  fin_cases i
  · simpa [CLPoly.Impl.StrictHGCDRefinement.hgcdStepEntries] using hNew0
  · simpa [CLPoly.Impl.StrictHGCDRefinement.hgcdStepEntries] using hNew1
  · have hptr : hgcdMatPtr row01.matrix h01 (2 : Fin 4) =
        hgcdMatPtr row23.matrix h23 (2 : Fin 4) := by
      simpa only using (hframe01 (0 : Fin 2)).1
    have hlen : hgcdMatLen row01.matrix h01 (2 : Fin 4) =
        hgcdMatLen row23.matrix h23 (2 : Fin 4) := by
      simpa only using (hframe01 (0 : Fin 2)).2
    change RawDensePolyRep this row01.heap
      (hgcdMatPtr row01.matrix h01 (2 : Fin 4))
      (hgcdMatLen row01.matrix h01 (2 : Fin 4))
      (entries 3 + quotient * entries 2)
    rw [hptr, hlen]
    exact hNew2
  · have hptr : hgcdMatPtr row01.matrix h01 (3 : Fin 4) =
        hgcdMatPtr row23.matrix h23 (3 : Fin 4) := by
      simpa only using (hframe01 (1 : Fin 2)).1
    have hlen : hgcdMatLen row01.matrix h01 (3 : Fin 4) =
        hgcdMatLen row23.matrix h23 (3 : Fin 4) := by
      simpa only using (hframe01 (1 : Fin 2)).2
    change RawDensePolyRep this row01.heap
      (hgcdMatPtr row01.matrix h01 (3 : Fin 4))
      (hgcdMatLen row01.matrix h01 (3 : Fin 4)) (entries 2)
    rw [hptr, hlen]
    exact hNew3

/-- An arbitrary external live raw polynomial survives both actual generated
row updates.  The proof follows the first call, reconstructs the unchanged
`(0,1)` descriptors and representations, and then follows the second call. -/
theorem hgcdTwoRowUpdates_preserve_guard (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (M : HgcdMat) (Q : RawPtr UInt64) (lenQ : Nat)
    (T : RawPtr UInt64) (lenT : Nat) (t scratch guard : RawPtr UInt64)
    (guardLen : Nat) (heap : RawHeap) (row23 row01 : MatRowUpdateResult)
    (hM : M.Valid) (h23 : row23.matrix.Valid)
    (quotient guardPoly : Polynomial (ZMod this._p.toNat))
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (workspace23 : MatRowUpdateWorkspace M (2 : Fin 4) (3 : Fin 4)
      Q lenQ T t scratch heap hM)
    (guard23 : MatRowUpdateGuardWorkspace M (2 : Fin 4) Q lenQ T t
      scratch guard guardLen hM)
    (workspace01 : MatRowUpdateWorkspace row23.matrix (0 : Fin 4)
      (1 : Fin 4) Q lenQ row23.T row23.t scratch row23.heap h23)
    (guard01 : MatRowUpdateGuardWorkspace row23.matrix (0 : Fin 4) Q
      lenQ row23.T row23.t scratch guard guardLen h23)
    (hQRep : RawDensePolyRep this heap Q lenQ quotient)
    (hMatrix : HgcdMatRawDenseRep this heap M entries hM)
    (hGuardRep : RawDensePolyRep this heap guard guardLen guardPoly)
    (hrow23 : dense_upoly_zp__mat_row_update_ir this M
      (2 : Fin 4) (3 : Fin 4) Q lenQ T lenT t scratch heap = .ok row23)
    (hrow01 : dense_upoly_zp__mat_row_update_ir this row23.matrix
      (0 : Fin 4) (1 : Fin 4) Q lenQ row23.T row23.lenT row23.t
      scratch row23.heap = .ok row01) :
    RawDensePolyRep this row01.heap guard guardLen guardPoly := by
  have hE0 := hMatrix (0 : Fin 4)
  have hE1 := hMatrix (1 : Fin 4)
  have hE2 := hMatrix (2 : Fin 4)
  have hE3 := hMatrix (3 : Fin 4)
  have hGuard23 := matRowUpdate_preserves_guard_of_workspaces this M
    (2 : Fin 4) (3 : Fin 4) Q lenQ T lenT t scratch guard guardLen heap
    row23 hM quotient (entries 2) (entries 3) guardPoly (by decide) hcfg hp
    workspace23 guard23 hQRep hE2 hE3 hGuardRep hrow23
  have hQ23 := matRowUpdate_preserves_quotient_of_workspace this M
    (2 : Fin 4) (3 : Fin 4) Q lenQ T lenT t scratch heap row23 hM
    quotient (entries 2) (entries 3) (by decide) hcfg hp workspace23
    hQRep hE2 hE3 hrow23
  have hE0Old23 := matRowUpdate_preserves_matrix_entry_of_workspace this M
    (2 : Fin 4) (3 : Fin 4) (0 : Fin 4) Q lenQ T lenT t scratch heap
    row23 hM quotient (entries 2) (entries 3) (entries 0) (by decide)
    (by decide) hcfg hp workspace23 hQRep hE2 hE3 hE0 hrow23
  have hE1Old23 := matRowUpdate_preserves_matrix_entry_of_workspace this M
    (2 : Fin 4) (3 : Fin 4) (1 : Fin 4) Q lenQ T lenT t scratch heap
    row23 hM quotient (entries 2) (entries 3) (entries 1) (by decide)
    (by decide) hcfg hp workspace23 hQRep hE2 hE3 hE1 hrow23
  rcases hgcdTwoRowUpdates_descriptor_frame this M Q lenQ T lenT t scratch
      heap row23 row01 hM hrow23 hrow01 with
    ⟨h23Frame, _, hframe23, _⟩
  have hh23Frame : h23Frame = h23 := Subsingleton.elim _ _
  subst h23Frame
  have hE0_23 : RawDensePolyRep this row23.heap
      (hgcdMatPtr row23.matrix h23 (0 : Fin 4))
      (hgcdMatLen row23.matrix h23 (0 : Fin 4)) (entries 0) := by
    have hptr : hgcdMatPtr row23.matrix h23 (0 : Fin 4) =
        hgcdMatPtr M hM (0 : Fin 4) := by
      simpa only using (hframe23 (0 : Fin 2)).1
    have hlen : hgcdMatLen row23.matrix h23 (0 : Fin 4) =
        hgcdMatLen M hM (0 : Fin 4) := by
      simpa only using (hframe23 (0 : Fin 2)).2
    rw [hptr, hlen]
    exact hE0Old23
  have hE1_23 : RawDensePolyRep this row23.heap
      (hgcdMatPtr row23.matrix h23 (1 : Fin 4))
      (hgcdMatLen row23.matrix h23 (1 : Fin 4)) (entries 1) := by
    have hptr : hgcdMatPtr row23.matrix h23 (1 : Fin 4) =
        hgcdMatPtr M hM (1 : Fin 4) := by
      simpa only using (hframe23 (1 : Fin 2)).1
    have hlen : hgcdMatLen row23.matrix h23 (1 : Fin 4) =
        hgcdMatLen M hM (1 : Fin 4) := by
      simpa only using (hframe23 (1 : Fin 2)).2
    rw [hptr, hlen]
    exact hE1Old23
  exact matRowUpdate_preserves_guard_of_workspaces this row23.matrix
    (0 : Fin 4) (1 : Fin 4) Q lenQ row23.T row23.lenT row23.t scratch
    guard guardLen row23.heap row01 h23 quotient (entries 0) (entries 1)
    guardPoly (by decide) hcfg hp workspace01 guard01 hQ23 hE0_23 hE1_23
    hGuard23 hrow01

/-- Descriptor-length effect of the two source-ordered matrix row updates.
This theorem follows both real generated calls; in particular the second
bound uses the quotient and the untouched `(0,1)` entries framed through the
first heap transition. -/
theorem hgcdTwoRowUpdates_length_bounds (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (M : HgcdMat) (Q : RawPtr UInt64) (lenQ : Nat)
    (T : RawPtr UInt64) (lenT : Nat) (t scratch : RawPtr UInt64)
    (heap : RawHeap) (row23 row01 : MatRowUpdateResult)
    (hM : M.Valid) (h23 : row23.matrix.Valid)
    (quotient : Polynomial (ZMod this._p.toNat))
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (workspace23 : MatRowUpdateWorkspace M (2 : Fin 4) (3 : Fin 4)
      Q lenQ T t scratch heap hM)
    (workspace01 : MatRowUpdateWorkspace row23.matrix (0 : Fin 4)
      (1 : Fin 4) Q lenQ row23.T row23.t scratch row23.heap h23)
    (hQRep : RawDensePolyRep this heap Q lenQ quotient)
    (hMatrix : HgcdMatRawDenseRep this heap M entries hM)
    (hrow23 : dense_upoly_zp__mat_row_update_ir this M
      (2 : Fin 4) (3 : Fin 4) Q lenQ T lenT t scratch heap = .ok row23)
    (hrow01 : dense_upoly_zp__mat_row_update_ir this row23.matrix
      (0 : Fin 4) (1 : Fin 4) Q lenQ row23.T row23.lenT row23.t
      scratch row23.heap = .ok row01) :
    ∃ h01 : row01.matrix.Valid,
      hgcdMatLen row01.matrix h01 (0 : Fin 4) ≤
        max (hgcdMatLen M hM (1 : Fin 4))
          (lenQ + hgcdMatLen M hM (0 : Fin 4) - 1) ∧
      hgcdMatLen row01.matrix h01 (1 : Fin 4) =
        hgcdMatLen M hM (0 : Fin 4) ∧
      hgcdMatLen row01.matrix h01 (2 : Fin 4) ≤
        max (hgcdMatLen M hM (3 : Fin 4))
          (lenQ + hgcdMatLen M hM (2 : Fin 4) - 1) ∧
      hgcdMatLen row01.matrix h01 (3 : Fin 4) =
        hgcdMatLen M hM (2 : Fin 4) := by
  have hE0 := hMatrix (0 : Fin 4)
  have hE1 := hMatrix (1 : Fin 4)
  have hE2 := hMatrix (2 : Fin 4)
  have hE3 := hMatrix (3 : Fin 4)
  rcases matRowUpdate_length_bound_of_workspace this M (2 : Fin 4)
      (3 : Fin 4) Q lenQ T lenT t scratch heap row23 hM quotient
      (entries 2) (entries 3) (by decide) hcfg hp workspace23 hQRep hE2 hE3
      hrow23 with ⟨h23', hlen2, hlen3, hother23⟩
  have hh23 : h23' = h23 := Subsingleton.elim _ _
  subst h23'
  have hQ23 := matRowUpdate_preserves_quotient_of_workspace this M
    (2 : Fin 4) (3 : Fin 4) Q lenQ T lenT t scratch heap row23 hM
    quotient (entries 2) (entries 3) (by decide) hcfg hp workspace23
    hQRep hE2 hE3 hrow23
  have hE0_23 := matRowUpdate_preserves_matrix_entry_of_workspace this M
    (2 : Fin 4) (3 : Fin 4) (0 : Fin 4) Q lenQ T lenT t scratch heap
    row23 hM quotient (entries 2) (entries 3) (entries 0) (by decide)
    (by decide) hcfg hp workspace23 hQRep hE2 hE3 hE0 hrow23
  have hE1_23 := matRowUpdate_preserves_matrix_entry_of_workspace this M
    (2 : Fin 4) (3 : Fin 4) (1 : Fin 4) Q lenQ T lenT t scratch heap
    row23 hM quotient (entries 2) (entries 3) (entries 1) (by decide)
    (by decide) hcfg hp workspace23 hQRep hE2 hE3 hE1 hrow23
  rcases hgcdTwoRowUpdates_descriptor_frame this M Q lenQ T lenT t scratch
      heap row23 row01 hM hrow23 hrow01 with
    ⟨h23Frame, _, hframe23, _⟩
  have hh23Frame : h23Frame = h23 := Subsingleton.elim _ _
  subst h23Frame
  have hE0_23' : RawDensePolyRep this row23.heap
      (hgcdMatPtr row23.matrix h23 (0 : Fin 4))
      (hgcdMatLen row23.matrix h23 (0 : Fin 4)) (entries 0) := by
    have hptr : hgcdMatPtr row23.matrix h23 (0 : Fin 4) =
        hgcdMatPtr M hM (0 : Fin 4) := by
      simpa only using (hframe23 (0 : Fin 2)).1
    have hlen : hgcdMatLen row23.matrix h23 (0 : Fin 4) =
        hgcdMatLen M hM (0 : Fin 4) := by
      simpa only using (hframe23 (0 : Fin 2)).2
    rw [hptr, hlen]
    exact hE0_23
  have hE1_23' : RawDensePolyRep this row23.heap
      (hgcdMatPtr row23.matrix h23 (1 : Fin 4))
      (hgcdMatLen row23.matrix h23 (1 : Fin 4)) (entries 1) := by
    have hptr : hgcdMatPtr row23.matrix h23 (1 : Fin 4) =
        hgcdMatPtr M hM (1 : Fin 4) := by
      simpa only using (hframe23 (1 : Fin 2)).1
    have hlen : hgcdMatLen row23.matrix h23 (1 : Fin 4) =
        hgcdMatLen M hM (1 : Fin 4) := by
      simpa only using (hframe23 (1 : Fin 2)).2
    rw [hptr, hlen]
    exact hE1_23
  rcases matRowUpdate_length_bound_of_workspace this row23.matrix
      (0 : Fin 4) (1 : Fin 4) Q lenQ row23.T row23.lenT row23.t scratch
      row23.heap row01 h23 quotient (entries 0) (entries 1) (by decide)
      hcfg hp workspace01 hQ23 hE0_23' hE1_23' hrow01 with
    ⟨h01, hlen0, hlen1, hother01⟩
  have h23len0 : hgcdMatLen row23.matrix h23 (0 : Fin 4) =
      hgcdMatLen M hM (0 : Fin 4) := hother23 (0 : Fin 4) (by decide)
        (by decide)
  have h23len1 : hgcdMatLen row23.matrix h23 (1 : Fin 4) =
      hgcdMatLen M hM (1 : Fin 4) := hother23 (1 : Fin 4) (by decide)
        (by decide)
  have h01len2 := hother01 (2 : Fin 4) (by decide) (by decide)
  have h01len3 := hother01 (3 : Fin 4) (by decide) (by decide)
  have hfinal0 : hgcdMatLen row01.matrix h01 (0 : Fin 4) ≤
      max (hgcdMatLen M hM (1 : Fin 4))
        (lenQ + hgcdMatLen M hM (0 : Fin 4) - 1) := by
    simpa only [h23len0, h23len1] using hlen0
  have hfinal1 : hgcdMatLen row01.matrix h01 (1 : Fin 4) =
      hgcdMatLen M hM (0 : Fin 4) := by
    simpa only [h23len0] using hlen1
  have hfinal2 : hgcdMatLen row01.matrix h01 (2 : Fin 4) ≤
      max (hgcdMatLen M hM (3 : Fin 4))
        (lenQ + hgcdMatLen M hM (2 : Fin 4) - 1) := h01len2 ▸ hlen2
  have hfinal3 : hgcdMatLen row01.matrix h01 (3 : Fin 4) =
      hgcdMatLen M hM (2 : Fin 4) := h01len3.trans hlen3
  exact ⟨h01, hfinal0, hfinal1, hfinal2, hfinal3⟩

/-- One complete nonterminal HGCD Euclidean iteration, from the actual
generated divrem call through both actual row updates, preserves the raw
state, matrix transform, signed determinant, normalized gcd, and the strict
well-founded measure. -/
theorem hgcdIterationCalls_refine (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (M : HgcdMat) (Q : RawPtr UInt64) (W3 : RawPtr Word3)
    (scratch A B R : RawPtr UInt64) (lenA lenB : Nat)
    (T : RawPtr UInt64) (lenT : Nat) (t : RawPtr UInt64)
    (heap heap1 : RawHeap) (lenQ lenR : Nat)
    (row23 row01 : MatRowUpdateResult) (hM : M.Valid)
    (h23 : row23.matrix.Valid)
    (left right dividend divisor : Polynomial (ZMod this._p.toNat))
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (sgn : Int) (hcfg : DensePreinvConfigured this)
    (hp : 1 < this._p.toNat) (hlenB : 0 < lenB)
    (hARep : RawDensePolyRep this heap A lenA dividend)
    (hBRep : RawDensePolyRep this heap B lenB divisor)
    (hMatrix : HgcdMatRawDenseRep this heap M entries hM)
    (htransform : CLPoly.Impl.StrictHGCDRefinement.HgcdTransform left right
      dividend divisor (entries 0) (entries 1) (entries 2) (entries 3))
    (hdet : CLPoly.Impl.StrictHGCDRefinement.HgcdSignedDet sgn
      (entries 0) (entries 1) (entries 2) (entries 3))
    (hQ : heap.ValidU64Slice Q (lenA - (lenB - 1)))
    (hR : heap.ValidU64Slice R (Nat.min lenA (lenB - 1)))
    (hW3 : heap.ValidWord3Slice W3 lenA)
    (hqCapacity : lenA - (lenB - 1) < limbBase)
    (hRA : R.region ≠ A.region) (hWA : W3.region ≠ A.region)
    (hWB : W3.region ≠ B.region) (hQB : Q.region ≠ B.region)
    (hQW : Q.region ≠ W3.region) (hRW : R.region ≠ W3.region)
    (hRQ : R.region ≠ Q.region) (hRB : R.region ≠ B.region)
    (hQMatrix : ∀ i : Fin 4, Q.region ≠ (hgcdMatPtr M hM i).region)
    (hRMatrix : ∀ i : Fin 4, R.region ≠ (hgcdMatPtr M hM i).region)
    (hW3Matrix : ∀ i : Fin 4, W3.region ≠ (hgcdMatPtr M hM i).region)
    (workspace23 : MatRowUpdateWorkspace M (2 : Fin 4) (3 : Fin 4)
      Q lenQ T t scratch heap1 hM)
    (workspace01 : MatRowUpdateWorkspace row23.matrix (0 : Fin 4)
      (1 : Fin 4) Q lenQ row23.T row23.t scratch row23.heap h23)
    (divisorGuard23 : MatRowUpdateGuardWorkspace M (2 : Fin 4) Q lenQ
      T t scratch B lenB hM)
    (divisorGuard01 : MatRowUpdateGuardWorkspace row23.matrix (0 : Fin 4)
      Q lenQ row23.T row23.t scratch B lenB h23)
    (remainderGuard23 : MatRowUpdateGuardWorkspace M (2 : Fin 4) Q lenQ
      T t scratch R lenR hM)
    (remainderGuard01 : MatRowUpdateGuardWorkspace row23.matrix (0 : Fin 4)
      Q lenQ row23.T row23.t scratch R lenR h23)
    (hdiv : Generated.StrictDivrem.dense_upoly_zp__poly_divrem_ir this Q R
      A lenA B lenB W3 heap = .ok (heap1, lenQ, lenR))
    (hrow23 : dense_upoly_zp__mat_row_update_ir this M (2 : Fin 4)
      (3 : Fin 4) Q lenQ T lenT t scratch heap1 = .ok row23)
    (hrow01 : dense_upoly_zp__mat_row_update_ir this row23.matrix
      (0 : Fin 4) (1 : Fin 4) Q lenQ row23.T row23.lenT row23.t scratch
      row23.heap = .ok row01) :
    ∃ quotient remainder h01,
      RawDensePolyRep this heap1 Q lenQ quotient ∧
      HgcdMatRawDenseRep this heap1 M entries hM ∧
      HgcdMatRawDenseRep this row01.heap row01.matrix
        (CLPoly.Impl.StrictHGCDRefinement.hgcdStepEntries quotient entries)
        h01 ∧
      RawDensePolyRep this row01.heap B lenB divisor ∧
      RawDensePolyRep this row01.heap R lenR remainder ∧
      CLPoly.Impl.StrictHGCDRefinement.HgcdTransform left right divisor
        remainder
        (CLPoly.Impl.StrictHGCDRefinement.hgcdStepEntries quotient entries 0)
        (CLPoly.Impl.StrictHGCDRefinement.hgcdStepEntries quotient entries 1)
        (CLPoly.Impl.StrictHGCDRefinement.hgcdStepEntries quotient entries 2)
        (CLPoly.Impl.StrictHGCDRefinement.hgcdStepEntries quotient entries 3) ∧
      CLPoly.Impl.StrictHGCDRefinement.HgcdSignedDet (-sgn)
        (CLPoly.Impl.StrictHGCDRefinement.hgcdStepEntries quotient entries 0)
        (CLPoly.Impl.StrictHGCDRefinement.hgcdStepEntries quotient entries 1)
        (CLPoly.Impl.StrictHGCDRefinement.hgcdStepEntries quotient entries 2)
        (CLPoly.Impl.StrictHGCDRefinement.hgcdStepEntries quotient entries 3) ∧
      normalize (EuclideanDomain.gcd dividend divisor) =
        normalize (EuclideanDomain.gcd divisor remainder) ∧
      lenQ ≤ lenA - (lenB - 1) ∧ lenR < lenB := by
  rcases polyDivrem_next_state this Q R A B lenA lenB W3 heap dividend
      divisor hlenB hARep hBRep hQ hR hW3 hqCapacity hRA hWA hWB hQB hQW
      hRW hRQ hRB hcfg with
    ⟨semanticHeap, semanticLenQ, semanticLenR, quotient, remainder,
      hsemantic, hQRep, hBRep1, hRRep1, hdivision, hgcd, _, hlenQ, _, hlt⟩
  have heq : (semanticHeap, semanticLenQ, semanticLenR) =
      (heap1, lenQ, lenR) := Except.ok.inj (hsemantic.symm.trans hdiv)
  cases heq
  have hMatrix1 := polyDivrem_preserves_hgcdMatRawDenseRep this M Q R A B
    lenA lenB W3 heap heap1 lenQ lenR entries hM hARep.1 hBRep.1 hQ hR
    hW3 hQMatrix hRMatrix hW3Matrix hMatrix hdiv
  rcases hgcdTwoRowUpdates_refine_matrix this M Q lenQ T lenT t scratch
      heap1 row23 row01 hM h23 quotient entries hcfg hp workspace23
      workspace01 hQRep hMatrix1 hrow23 hrow01 with ⟨h01, hMatrix01⟩
  have hDivisor01 := hgcdTwoRowUpdates_preserve_guard this M Q lenQ T
    lenT t scratch B lenB heap1 row23 row01 hM h23 quotient divisor entries
    hcfg hp workspace23 divisorGuard23 workspace01 divisorGuard01 hQRep
    hMatrix1 hBRep1 hrow23 hrow01
  have hRemainder01 := hgcdTwoRowUpdates_preserve_guard this M Q lenQ T
    lenT t scratch R lenR heap1 row23 row01 hM h23 quotient remainder
    entries hcfg hp workspace23 remainderGuard23 workspace01
    remainderGuard01 hQRep hMatrix1 hRRep1 hrow23 hrow01
  have htransform' :=
    CLPoly.Impl.StrictHGCDRefinement.hgcdStepEntries_preserves_transform
      left right dividend divisor remainder quotient entries htransform
      hdivision
  have hdet' :=
    CLPoly.Impl.StrictHGCDRefinement.hgcdStepEntries_preserves_signedDet
      sgn quotient entries hdet
  exact ⟨quotient, remainder, h01, hQRep, hMatrix1, hMatrix01, hDivisor01,
    hRemainder01,
    htransform', hdet', hgcd, hlenQ, hlt⟩

/-- Semantic invariant carried by the well-founded generated HGCD loop.
The original pair is related to the current raw pair by the actual matrix,
whose determinant sign is synchronized with the source `sgn` field. -/
structure HgcdIterRawInvariant (this : DenseUPolyZp)
    (left right currentA currentB : Polynomial (ZMod this._p.toNat))
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (state : HgcdIterState) (hM : state.matrix.Valid) : Prop where
  matrixRep : HgcdMatRawDenseRep this state.heap state.matrix entries hM
  aRep : RawDensePolyRep this state.heap state.A state.lenA currentA
  bRep : RawDensePolyRep this state.heap state.B state.lenB currentB
  transform : CLPoly.Impl.StrictHGCDRefinement.HgcdTransform left right
    currentA currentB (entries 0) (entries 1) (entries 2) (entries 3)
  signedDet : CLPoly.Impl.StrictHGCDRefinement.HgcdSignedDet state.sgn
    (entries 0) (entries 1) (entries 2) (entries 3)

/-- Descriptor-length invariant needed by the enclosing recursive HGCD.
Each matrix coefficient is paired with the current Euclidean operand whose
degree it complements. -/
structure HgcdMatrixLengthInvariant (inputLength : Nat)
    (state : HgcdIterState) (hM : state.matrix.Valid) : Prop where
  row0A : hgcdMatLen state.matrix hM (0 : Fin 4) + state.lenA ≤
    inputLength + 1
  row1B : hgcdMatLen state.matrix hM (1 : Fin 4) + state.lenB ≤
    inputLength + 1
  row2A : hgcdMatLen state.matrix hM (2 : Fin 4) + state.lenA ≤
    inputLength + 1
  row3B : hgcdMatLen state.matrix hM (3 : Fin 4) + state.lenB ≤
    inputLength + 1
  row1A : hgcdMatLen state.matrix hM (1 : Fin 4) + state.lenA ≤
    inputLength + 1
  row3A : hgcdMatLen state.matrix hM (3 : Fin 4) + state.lenA ≤
    inputLength + 1

/-- Uniform coefficient bound valid while HGCD processes only operands above
the source half-length threshold.  It is the degree-separation fact used by
the final low/high reconstruction. -/
def HgcdMatrixCoefficientBound (inputLength : Nat)
    (state : HgcdIterState) (hM : state.matrix.Valid) : Prop :=
  ∀ i : Fin 4,
    hgcdMatLen state.matrix hM i ≤ inputLength - inputLength / 2

/-- Unified induction result for the complete recursive HGCD call.  It
contains only semantics derived from the returned raw heap plus the exact
length facts required by the enclosing reconstruction and recursive calls. -/
structure HgcdRecursiveRawInvariant (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (left right finalA finalB : Polynomial (ZMod this._p.toNat))
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (computeM : Bool) (outA outB : RawPtr UInt64) (inputLength : Nat)
    (result : HgcdRecursiveResult) : Prop where
  aRep : RawDensePolyRep this result.heap outA result.lenA finalA
  bRep : RawDensePolyRep this result.heap outB result.lenB finalB
  matrixSemantics : computeM = true →
    HgcdMatRawDenseRep this result.heap result.matrix entries result.valid ∧
    CLPoly.Impl.StrictHGCDRefinement.HgcdTransform left right finalA finalB
      (entries 0) (entries 1) (entries 2) (entries 3) ∧
    CLPoly.Impl.StrictHGCDRefinement.HgcdSignedDet result.sgn
      (entries 0) (entries 1) (entries 2) (entries 3)
  gcdPreserved : normalize (EuclideanDomain.gcd left right) =
    normalize (EuclideanDomain.gcd finalA finalB)
  stopped : result.lenB < inputLength / 2 + 1
  lengths : computeM = true → HgcdRecursiveLengthInvariant inputLength result

/-- Package the exact early-return execution as the common recursive result.
The algebraic premises are precisely the facts established by the preceding
recursive call and paired reconstruction; this theorem only transports them
through the generated output and optional matrix copies. -/
theorem hgcdRecursiveEarlyReturn_rawInvariant (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (M R : HgcdMat) (hM : M.Valid) (hR : R.Valid) (computeM : Bool)
    (A B a2 b2 : RawPtr UInt64) (lenA2 lenB2 inputLength : Nat)
    (sgn : Int) (left right finalA finalB : Polynomial (ZMod this._p.toNat))
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (heap : RawHeap) (result : HgcdRecursiveEarlyResult)
    (hwork : HgcdEarlyReturnRefineWorkspace heap M R hM hR
      A B a2 b2 lenA2 lenB2)
    (hARep : RawDensePolyRep this heap a2 lenA2 finalA)
    (hBRep : RawDensePolyRep this heap b2 lenB2 finalB)
    (hMatrix : HgcdMatRawDenseRep this heap R entries hR)
    (hTransform : CLPoly.Impl.StrictHGCDRefinement.HgcdTransform left right
      finalA finalB (entries 0) (entries 1) (entries 2) (entries 3))
    (hDet : CLPoly.Impl.StrictHGCDRefinement.HgcdSignedDet sgn
      (entries 0) (entries 1) (entries 2) (entries 3))
    (hGcd : normalize (EuclideanDomain.gcd left right) =
      normalize (EuclideanDomain.gcd finalA finalB))
    (hstop : lenB2 < inputLength / 2 + 1)
    (hLength : HgcdRecursiveLengthInvariant inputLength
      ⟨heap, R, hR, lenA2, lenB2, sgn⟩)
    (hrun : hgcdRecursiveEarlyReturn M R hM hR computeM A B a2 b2
      lenA2 lenB2 sgn heap = .ok result) :
    ∃ hResult : result.matrix.Valid,
      HgcdRecursiveRawInvariant this left right finalA finalB entries computeM
        A B inputLength (result.toResult hResult) := by
  rcases hgcdRecursiveEarlyReturn_refines this M R hM hR computeM A B a2 b2
      lenA2 lenB2 sgn finalA finalB entries heap hwork hARep hBRep hMatrix with
    ⟨actual, hactual, _, hlenA, hlenB, hsgn, hAResult, hBResult,
      hActualM, hmatrixTrue, _⟩
  have heq : actual = result := Except.ok.inj (hactual.symm.trans hrun)
  subst actual
  refine ⟨hActualM, {
    aRep := by simpa [HgcdRecursiveEarlyResult.toResult, hlenA] using hAResult
    bRep := by simpa [HgcdRecursiveEarlyResult.toResult, hlenB] using hBResult
    matrixSemantics := ?_
    gcdPreserved := hGcd
    stopped := by
      simpa [HgcdRecursiveEarlyResult.toResult, hlenB] using hstop
    lengths := ?_ }⟩
  · intro hcompute
    have hcopied := hmatrixTrue hcompute
    exact ⟨by
        simpa [HgcdRecursiveEarlyResult.toResult] using hcopied.2,
      hTransform,
      by simpa [HgcdRecursiveEarlyResult.toResult, hsgn] using hDet⟩
  · intro hcompute
    have hcopied := hmatrixTrue hcompute
    constructor
    · simpa [HgcdRecursiveEarlyResult.toResult, hgcdMatLenRaw,
        hlenA, hcopied.1] using hLength.row0A
    · simpa [HgcdRecursiveEarlyResult.toResult, hgcdMatLenRaw,
        hlenB, hcopied.1] using hLength.row1B
    · simpa [HgcdRecursiveEarlyResult.toResult, hgcdMatLenRaw,
        hlenA, hcopied.1] using hLength.row2A
    · simpa [HgcdRecursiveEarlyResult.toResult, hgcdMatLenRaw,
        hlenB, hcopied.1] using hLength.row3B
    · simpa [HgcdRecursiveEarlyResult.toResult, hgcdMatLenRaw,
        hlenA, hcopied.1] using hLength.row1A
    · simpa [HgcdRecursiveEarlyResult.toResult, hgcdMatLenRaw,
        hlenA, hcopied.1] using hLength.row3A
    · simpa [HgcdRecursiveEarlyResult.toResult, hlenA] using
        hLength.inputBound
    · simpa [HgcdRecursiveEarlyResult.toResult, hlenB] using
        hLength.stopped
    · simpa [HgcdRecursiveEarlyResult.toResult, hlenA] using
        hLength.positive
    · simpa [HgcdRecursiveEarlyResult.toResult, hlenA] using
        hLength.aboveHalf
    · intro i
      simpa [HgcdRecursiveEarlyResult.toResult, hgcdMatLenRaw,
        hcopied.1] using hLength.coeffBound i

/-- Arithmetic closure for one Euclidean matrix step.  The quotient bound
comes from the real generated divrem call and the four descriptor bounds
come from the two real row updates. -/
theorem hgcdMatrixLengthInvariant_step
    (inputLength lenA lenB lenQ lenR : Nat)
    (M nextM : HgcdMat) (hM : M.Valid) (hNextM : nextM.Valid)
    (hinvariant :
      hgcdMatLen M hM (0 : Fin 4) + lenA ≤ inputLength + 1 ∧
      hgcdMatLen M hM (1 : Fin 4) + lenB ≤ inputLength + 1 ∧
      hgcdMatLen M hM (2 : Fin 4) + lenA ≤ inputLength + 1 ∧
      hgcdMatLen M hM (3 : Fin 4) + lenB ≤ inputLength + 1)
    (horder : lenB ≤ lenA)
    (hlenQ : lenQ ≤ lenA - (lenB - 1)) (hlenR : lenR < lenB)
    (hlen0 : hgcdMatLen nextM hNextM (0 : Fin 4) ≤
      max (hgcdMatLen M hM (1 : Fin 4))
        (lenQ + hgcdMatLen M hM (0 : Fin 4) - 1))
    (hlen1 : hgcdMatLen nextM hNextM (1 : Fin 4) =
      hgcdMatLen M hM (0 : Fin 4))
    (hlen2 : hgcdMatLen nextM hNextM (2 : Fin 4) ≤
      max (hgcdMatLen M hM (3 : Fin 4))
        (lenQ + hgcdMatLen M hM (2 : Fin 4) - 1))
    (hlen3 : hgcdMatLen nextM hNextM (3 : Fin 4) =
      hgcdMatLen M hM (2 : Fin 4)) :
    hgcdMatLen nextM hNextM (0 : Fin 4) + lenB ≤ inputLength + 1 ∧
    hgcdMatLen nextM hNextM (1 : Fin 4) + lenR ≤ inputLength + 1 ∧
    hgcdMatLen nextM hNextM (2 : Fin 4) + lenB ≤ inputLength + 1 ∧
    hgcdMatLen nextM hNextM (3 : Fin 4) + lenR ≤ inputLength + 1 := by
  rcases hinvariant with ⟨h0A, h1B, h2A, h3B⟩
  have hq0 : lenQ + hgcdMatLen M hM (0 : Fin 4) - 1 + lenB ≤
      inputLength + 1 := by omega
  have hq2 : lenQ + hgcdMatLen M hM (2 : Fin 4) - 1 + lenB ≤
      inputLength + 1 := by omega
  have hmax0 : max (hgcdMatLen M hM (1 : Fin 4))
        (lenQ + hgcdMatLen M hM (0 : Fin 4) - 1) + lenB ≤
      inputLength + 1 := by
    rcases le_total (hgcdMatLen M hM (1 : Fin 4))
        (lenQ + hgcdMatLen M hM (0 : Fin 4) - 1) with hle | hle
    · rw [max_eq_right hle]
      exact hq0
    · rw [max_eq_left hle]
      omega
  have hmax2 : max (hgcdMatLen M hM (3 : Fin 4))
        (lenQ + hgcdMatLen M hM (2 : Fin 4) - 1) + lenB ≤
      inputLength + 1 := by
    rcases le_total (hgcdMatLen M hM (3 : Fin 4))
        (lenQ + hgcdMatLen M hM (2 : Fin 4) - 1) with hle | hle
    · rw [max_eq_right hle]
      exact hq2
    · rw [max_eq_left hle]
      omega
  constructor
  · omega
  constructor
  · omega
  constructor <;> omega

/-- One concrete divrem/row-update step preserves the uniform coefficient
bound while the source loop guard keeps the divisor above half length. -/
theorem hgcdMatrixCoefficientBound_step
    (inputLength lenA lenB lenQ : Nat)
    (M nextM : HgcdMat) (hM : M.Valid) (hNextM : nextM.Valid)
    (hinvariant :
      hgcdMatLen M hM (0 : Fin 4) + lenA ≤ inputLength + 1 ∧
      hgcdMatLen M hM (2 : Fin 4) + lenA ≤ inputLength + 1)
    (hcoeff : ∀ i : Fin 4,
      hgcdMatLen M hM i ≤ inputLength - inputLength / 2)
    (hguard : inputLength / 2 + 1 ≤ lenB)
    (horder : lenB ≤ lenA)
    (hlenQ : lenQ ≤ lenA - (lenB - 1))
    (hlen0 : hgcdMatLen nextM hNextM (0 : Fin 4) ≤
      max (hgcdMatLen M hM (1 : Fin 4))
        (lenQ + hgcdMatLen M hM (0 : Fin 4) - 1))
    (hlen1 : hgcdMatLen nextM hNextM (1 : Fin 4) =
      hgcdMatLen M hM (0 : Fin 4))
    (hlen2 : hgcdMatLen nextM hNextM (2 : Fin 4) ≤
      max (hgcdMatLen M hM (3 : Fin 4))
        (lenQ + hgcdMatLen M hM (2 : Fin 4) - 1))
    (hlen3 : hgcdMatLen nextM hNextM (3 : Fin 4) =
      hgcdMatLen M hM (2 : Fin 4)) :
    ∀ i : Fin 4,
      hgcdMatLen nextM hNextM i ≤ inputLength - inputLength / 2 := by
  have hq0 : lenQ + hgcdMatLen M hM (0 : Fin 4) - 1 ≤
      inputLength - inputLength / 2 := by
    have := hinvariant.1
    omega
  have hq2 : lenQ + hgcdMatLen M hM (2 : Fin 4) - 1 ≤
      inputLength - inputLength / 2 := by
    have := hinvariant.2
    omega
  have h0 := hcoeff (0 : Fin 4)
  have h1 := hcoeff (1 : Fin 4)
  have h2 := hcoeff (2 : Fin 4)
  have h3 := hcoeff (3 : Fin 4)
  have hn0 : hgcdMatLen nextM hNextM (0 : Fin 4) ≤
      inputLength - inputLength / 2 := by
    exact hlen0.trans (max_le h1 hq0)
  have hn2 : hgcdMatLen nextM hNextM (2 : Fin 4) ≤
      inputLength - inputLength / 2 := by
    exact hlen2.trans (max_le h3 hq2)
  intro i
  fin_cases i
  · exact hn0
  · simpa [hlen1] using h0
  · exact hn2
  · simpa [hlen3] using h2

/-- The odd matrix entries also complement the current leading operand.
After a source Euclidean step they are the previous even entries, while the
new A descriptor is the previous B descriptor. -/
theorem hgcdMatrixOddALengthInvariant_step
    (inputLength lenA lenB : Nat)
    (M nextM : HgcdMat) (hM : M.Valid) (hNextM : nextM.Valid)
    (hrow0A : hgcdMatLen M hM (0 : Fin 4) + lenA ≤ inputLength + 1)
    (hrow2A : hgcdMatLen M hM (2 : Fin 4) + lenA ≤ inputLength + 1)
    (horder : lenB ≤ lenA)
    (hlen1 : hgcdMatLen nextM hNextM (1 : Fin 4) =
      hgcdMatLen M hM (0 : Fin 4))
    (hlen3 : hgcdMatLen nextM hNextM (3 : Fin 4) =
      hgcdMatLen M hM (2 : Fin 4)) :
    hgcdMatLen nextM hNextM (1 : Fin 4) + lenB ≤ inputLength + 1 ∧
      hgcdMatLen nextM hNextM (3 : Fin 4) + lenB ≤ inputLength + 1 := by
  constructor <;> omega

/-- The real `_mat_one` descriptor and the two source-ordered input copies
establish the matrix-length invariant. -/
theorem hgcdIterInit_matrixLengthInvariant
    (M : HgcdMat) (A B T t : RawPtr UInt64) (lenT : Nat)
    (a : RawPtr UInt64) (lenA : Nat) (b : RawPtr UInt64) (lenB : Nat)
    (heap : RawHeap) (initial : HgcdIterState) (hM : initial.matrix.Valid)
    (horder : lenB ≤ lenA)
    (hrun : hgcdIterInit M A B T t lenT a lenA b lenB heap = .ok initial) :
    HgcdMatrixLengthInvariant lenA initial hM := by
  have hlens := hgcdIterInit_lengths M A B T t lenT a lenA b lenB heap
    initial hrun
  rcases hlens with ⟨hA, hB, hMatrix⟩
  constructor <;>
    simp [hgcdMatLen, hA, hB, hMatrix] <;> omega

/-- The identity matrix created by the exact iterator initialization already
satisfies the uniform half-length coefficient bound. -/
theorem hgcdIterInit_matrixCoefficientBound
    (M : HgcdMat) (A B T t : RawPtr UInt64) (lenT : Nat)
    (a : RawPtr UInt64) (lenA : Nat) (b : RawPtr UInt64) (lenB : Nat)
    (heap : RawHeap) (initial : HgcdIterState) (hM : initial.matrix.Valid)
    (hlenAPos : 0 < lenA)
    (hrun : hgcdIterInit M A B T t lenT a lenA b lenB heap = .ok initial) :
    HgcdMatrixCoefficientBound lenA initial hM := by
  have hlens := hgcdIterInit_lengths M A B T t lenT a lenA b lenB heap
    initial hrun
  rcases hlens with ⟨_, _, hMatrix⟩
  intro i
  fin_cases i <;> simp [HgcdMatrixCoefficientBound, hgcdMatLen, hMatrix] <;>
    omega

/-- Purely physical obligations for one concrete successful nonterminal
iteration.  Every field concerns validity, capacity, or allocation
separation of the exact buffers used by the generated calls. -/
structure HgcdIterationWorkspace (this : DenseUPolyZp)
    (Q : RawPtr UInt64) (W3 : RawPtr Word3) (scratch : RawPtr UInt64)
    (state : HgcdIterState) (heap1 : RawHeap) (lenQ lenR : Nat)
    (row23 row01 : MatRowUpdateResult) (hM : state.matrix.Valid) : Prop where
  validQ : state.heap.ValidU64Slice Q
    (state.lenA - (state.lenB - 1))
  validR : state.heap.ValidU64Slice state.T
    (Nat.min state.lenA (state.lenB - 1))
  validW3 : state.heap.ValidWord3Slice W3 state.lenA
  quotientCapacity : state.lenA - (state.lenB - 1) < limbBase
  rA : state.T.region ≠ state.A.region
  wA : W3.region ≠ state.A.region
  wB : W3.region ≠ state.B.region
  qB : Q.region ≠ state.B.region
  qW : Q.region ≠ W3.region
  rW : state.T.region ≠ W3.region
  rQ : state.T.region ≠ Q.region
  rB : state.T.region ≠ state.B.region
  qMatrix : ∀ i : Fin 4,
    Q.region ≠ (hgcdMatPtr state.matrix hM i).region
  rMatrix : ∀ i : Fin 4,
    state.T.region ≠ (hgcdMatPtr state.matrix hM i).region
  wMatrix : ∀ i : Fin 4,
    W3.region ≠ (hgcdMatPtr state.matrix hM i).region
  matrix23Valid : row23.matrix.Valid
  row23Workspace : MatRowUpdateWorkspace state.matrix (2 : Fin 4)
    (3 : Fin 4) Q lenQ state.A state.t scratch heap1 hM
  row01Workspace : MatRowUpdateWorkspace row23.matrix (0 : Fin 4)
    (1 : Fin 4) Q lenQ row23.T row23.t scratch row23.heap matrix23Valid
  divisorGuard23 : MatRowUpdateGuardWorkspace state.matrix (2 : Fin 4)
    Q lenQ state.A state.t scratch state.B state.lenB hM
  divisorGuard01 : MatRowUpdateGuardWorkspace row23.matrix (0 : Fin 4)
    Q lenQ row23.T row23.t scratch state.B state.lenB matrix23Valid
  remainderGuard23 : MatRowUpdateGuardWorkspace state.matrix (2 : Fin 4)
    Q lenQ state.A state.t scratch state.T lenR hM
  remainderGuard01 : MatRowUpdateGuardWorkspace row23.matrix (0 : Fin 4)
    Q lenQ row23.T row23.t scratch state.T lenR matrix23Valid

/-- Separation-logic style precondition for the whole loop: for every
concrete successful source step, the caller can discharge the purely
physical workspace obligations above.  It contains no polynomial or L2
result. -/
def HgcdLoopWorkspaceProvider (this : DenseUPolyZp) (m : Nat)
    (Q : RawPtr UInt64) (W3 : RawPtr Word3)
    (scratch : RawPtr UInt64) : Prop :=
  ∀ (state : HgcdIterState) (hM : state.matrix.Valid)
    (heap1 : RawHeap) (lenQ lenR : Nat)
    (row23 row01 : MatRowUpdateResult),
    state.lenB ≥ m + 1 →
    Generated.StrictDivrem.dense_upoly_zp__poly_divrem_ir this Q state.T
      state.A state.lenA state.B state.lenB W3 state.heap =
        .ok (heap1, lenQ, lenR) →
    dense_upoly_zp__mat_row_update_ir this state.matrix (2 : Fin 4)
      (3 : Fin 4) Q lenQ state.A state.lenT state.t scratch heap1 =
        .ok row23 →
    dense_upoly_zp__mat_row_update_ir this row23.matrix (0 : Fin 4)
      (1 : Fin 4) Q lenQ row23.T row23.lenT row23.t scratch row23.heap =
        .ok row01 →
    HgcdIterationWorkspace this Q W3 scratch state heap1 lenQ lenR row23
      row01 hM

/-- End-to-end semantic refinement of the actual well-founded generated
`hgcdIterLoop`.  Recursion is on the source measure `state.lenB`; each
recursive call is justified by the exact divrem result's `lenR < lenB`; no
bounded execution counter or alternate L2 execution is introduced. -/
theorem hgcdIterLoop_refines (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (m : Nat) (Q : RawPtr UInt64) (W3 : RawPtr Word3)
    (scratch : RawPtr UInt64)
    (left right : Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (physical : HgcdLoopWorkspaceProvider this m Q W3 scratch) :
    ∀ (state final : HgcdIterState)
      (currentA currentB : Polynomial (ZMod this._p.toNat))
      (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
      (hM : state.matrix.Valid),
      HgcdIterRawInvariant this left right currentA currentB entries state hM →
      hgcdIterLoop this m Q W3 scratch state = .ok final →
      ∃ finalA finalB finalEntries hFinalM,
        HgcdIterRawInvariant this left right finalA finalB finalEntries final
          hFinalM ∧
        normalize (EuclideanDomain.gcd currentA currentB) =
          normalize (EuclideanDomain.gcd finalA finalB) ∧
        final.lenB < m + 1
  | state, final, currentA, currentB, entries, hM, hinvariant, hrun => by
    by_cases hguard : state.lenB ≥ m + 1
    · rcases hgcdIterLoop_step_shape this m Q W3 scratch state final hguard
        hrun with
      ⟨heap1, lenQ, lenR, row23, row01, hdiv, hrow23, hrow01, htail, hlt⟩
      have hworkspace := physical state hM heap1 lenQ lenR row23 row01
        hguard hdiv hrow23 hrow01
      rcases hgcdIterationCalls_refine this state.matrix Q W3 scratch state.A
          state.B state.T state.lenA state.lenB state.A state.lenT state.t
          state.heap heap1 lenQ lenR row23 row01 hM
          hworkspace.matrix23Valid left right currentA currentB entries
          state.sgn hcfg hp (by omega) hinvariant.aRep hinvariant.bRep
          hinvariant.matrixRep hinvariant.transform hinvariant.signedDet
          hworkspace.validQ hworkspace.validR hworkspace.validW3
          hworkspace.quotientCapacity hworkspace.rA hworkspace.wA
          hworkspace.wB hworkspace.qB hworkspace.qW hworkspace.rW
          hworkspace.rQ hworkspace.rB hworkspace.qMatrix
          hworkspace.rMatrix hworkspace.wMatrix hworkspace.row23Workspace
          hworkspace.row01Workspace hworkspace.divisorGuard23
          hworkspace.divisorGuard01 hworkspace.remainderGuard23
          hworkspace.remainderGuard01 hdiv hrow23 hrow01 with
        ⟨quotient, remainder, h01, _, _, hMatrix01, hDivisor01, hRemainder01,
          htransform, hdet, hgcdStep, _, hlt'⟩
      let next : HgcdIterState := {
        heap := row01.heap
        matrix := row01.matrix
        A := state.B
        lenA := state.lenB
        B := state.T
        lenB := lenR
        T := row01.T
        lenT := row01.lenT
        t := row01.t
        sgn := -state.sgn }
      have hnextInvariant : HgcdIterRawInvariant this left right currentB
          remainder
          (CLPoly.Impl.StrictHGCDRefinement.hgcdStepEntries quotient entries)
          next h01 := by
        exact ⟨hMatrix01, hDivisor01, hRemainder01, htransform, hdet⟩
      have htail' : hgcdIterLoop this m Q W3 scratch next = .ok final := by
        simpa [next] using htail
      rcases hgcdIterLoop_refines this m Q W3 scratch left right hcfg hp
          physical next final currentB remainder
          (CLPoly.Impl.StrictHGCDRefinement.hgcdStepEntries quotient entries)
          h01 hnextInvariant htail' with
        ⟨finalA, finalB, finalEntries, hFinalM, hfinalInvariant,
          hgcdRest, hstop⟩
      exact ⟨finalA, finalB, finalEntries, hFinalM, hfinalInvariant,
        hgcdStep.trans hgcdRest, hstop⟩
    · have hstop : state.lenB < m + 1 := by omega
      have hsame := hgcdIterLoop_stop this m Q W3 scratch state hstop
      have hfinal : state = final := Except.ok.inj (hsame.symm.trans hrun)
      subst final
      exact ⟨currentA, currentB, entries, hM, hinvariant, rfl, hstop⟩
termination_by state => state.lenB
decreasing_by exact hlt'

/-- The same well-founded source loop carries the descriptor-length
invariant needed by the enclosing recursive call.  This is deliberately a
separate theorem from semantic refinement so its result cannot be supplied
without following the generated heap execution. -/
theorem hgcdIterLoop_preserves_matrixLength (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (inputLength m : Nat) (Q : RawPtr UInt64) (W3 : RawPtr Word3)
    (scratch : RawPtr UInt64)
    (left right : Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (hm : m = inputLength / 2)
    (physical : HgcdLoopWorkspaceProvider this m Q W3 scratch) :
    ∀ (state final : HgcdIterState)
      (currentA currentB : Polynomial (ZMod this._p.toNat))
      (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
      (hM : state.matrix.Valid),
      HgcdIterRawInvariant this left right currentA currentB entries state hM →
      HgcdMatrixLengthInvariant inputLength state hM →
      HgcdMatrixCoefficientBound inputLength state hM →
      state.lenB ≤ state.lenA →
      state.lenA ≤ inputLength →
      0 < state.lenA →
      hgcdIterLoop this m Q W3 scratch state = .ok final →
      ∃ hFinalM : final.matrix.Valid,
        HgcdMatrixLengthInvariant inputLength final hFinalM ∧
        final.lenB ≤ final.lenA ∧ final.lenA ≤ inputLength ∧ 0 < final.lenA ∧
        HgcdMatrixCoefficientBound inputLength final hFinalM
  | state, final, currentA, currentB, entries, hM, hraw, hlength, hcoeff,
      horder, hinputBound, hpositive, hrun => by
    by_cases hguard : state.lenB ≥ m + 1
    · rcases hgcdIterLoop_step_shape this m Q W3 scratch state final hguard
        hrun with
      ⟨heap1, lenQ, lenR, row23, row01, hdiv, hrow23, hrow01, htail, hlt⟩
      have hworkspace := physical state hM heap1 lenQ lenR row23 row01
        hguard hdiv hrow23 hrow01
      rcases hgcdIterationCalls_refine this state.matrix Q W3 scratch state.A
          state.B state.T state.lenA state.lenB state.A state.lenT state.t
          state.heap heap1 lenQ lenR row23 row01 hM
          hworkspace.matrix23Valid left right currentA currentB entries
          state.sgn hcfg hp (by omega) hraw.aRep hraw.bRep hraw.matrixRep
          hraw.transform hraw.signedDet hworkspace.validQ hworkspace.validR
          hworkspace.validW3 hworkspace.quotientCapacity hworkspace.rA
          hworkspace.wA hworkspace.wB hworkspace.qB hworkspace.qW
          hworkspace.rW hworkspace.rQ hworkspace.rB hworkspace.qMatrix
          hworkspace.rMatrix hworkspace.wMatrix hworkspace.row23Workspace
          hworkspace.row01Workspace hworkspace.divisorGuard23
          hworkspace.divisorGuard01 hworkspace.remainderGuard23
          hworkspace.remainderGuard01 hdiv hrow23 hrow01 with
        ⟨quotient, remainder, h01, hQRep, hMatrix1, hMatrix01, hDivisor01,
          hRemainder01, htransform, hdet, _, hlenQ, hlenR⟩
      rcases hgcdTwoRowUpdates_length_bounds this state.matrix Q lenQ
          state.A state.lenT state.t scratch heap1 row23 row01 hM
          hworkspace.matrix23Valid quotient entries hcfg hp
          hworkspace.row23Workspace hworkspace.row01Workspace hQRep
          hMatrix1 hrow23 hrow01 with
        ⟨h01', hlen0, hlen1, hlen2, hlen3⟩
      have hh01 : h01' = h01 := Subsingleton.elim _ _
      subst h01'
      let next : HgcdIterState := {
        heap := row01.heap
        matrix := row01.matrix
        A := state.B
        lenA := state.lenB
        B := state.T
        lenB := lenR
        T := row01.T
        lenT := row01.lenT
        t := row01.t
        sgn := -state.sgn }
      have hnextRaw : HgcdIterRawInvariant this left right currentB remainder
          (CLPoly.Impl.StrictHGCDRefinement.hgcdStepEntries quotient entries)
          next h01 := ⟨hMatrix01, hDivisor01, hRemainder01, htransform, hdet⟩
      have hnextBounds := hgcdMatrixLengthInvariant_step inputLength
        state.lenA state.lenB lenQ lenR state.matrix row01.matrix hM h01
        ⟨hlength.row0A, hlength.row1B, hlength.row2A, hlength.row3B⟩
        horder hlenQ hlenR hlen0 hlen1 hlen2 hlen3
      have hnextOdd := hgcdMatrixOddALengthInvariant_step inputLength
        state.lenA state.lenB state.matrix row01.matrix hM h01
        hlength.row0A hlength.row2A horder hlen1 hlen3
      have hnextLength : HgcdMatrixLengthInvariant inputLength next h01 := by
        exact ⟨hnextBounds.1, hnextBounds.2.1, hnextBounds.2.2.1,
          hnextBounds.2.2.2, hnextOdd.1, hnextOdd.2⟩
      have hnextCoeff : HgcdMatrixCoefficientBound inputLength next h01 := by
        exact hgcdMatrixCoefficientBound_step inputLength state.lenA
          state.lenB lenQ state.matrix row01.matrix hM h01
          ⟨hlength.row0A, hlength.row2A⟩ hcoeff (by omega) horder hlenQ
          hlen0 hlen1 hlen2 hlen3
      have htail' : hgcdIterLoop this m Q W3 scratch next = .ok final := by
        simpa [next] using htail
      exact hgcdIterLoop_preserves_matrixLength this inputLength m Q W3
        scratch left right hcfg hp hm physical next final currentB remainder
        (CLPoly.Impl.StrictHGCDRefinement.hgcdStepEntries quotient entries)
        h01 hnextRaw hnextLength hnextCoeff
        (by simpa [next] using hlenR.le)
        (by simpa [next] using horder.trans hinputBound)
        (by simpa [next] using (show 0 < state.lenB by omega)) htail'
    · have hstop : state.lenB < m + 1 := by omega
      have hsame := hgcdIterLoop_stop this m Q W3 scratch state hstop
      have hfinal : state = final := Except.ok.inj (hsame.symm.trans hrun)
      subst final
      exact ⟨hM, hlength, horder, hinputBound, hpositive, hcoeff⟩
termination_by state => state.lenB
decreasing_by exact hlt

/-- The real loop cannot stop with its leading operand below the cutoff.
On every taken iteration the next A descriptor is exactly the old B
descriptor that satisfied the source guard `m + 1 ≤ lenB`. -/
theorem hgcdIterLoop_preserves_above_half (this : DenseUPolyZp)
    (m : Nat) (Q : RawPtr UInt64) (W3 : RawPtr Word3)
    (scratch : RawPtr UInt64) :
    ∀ (state final : HgcdIterState),
      m < state.lenA →
      hgcdIterLoop this m Q W3 scratch state = .ok final →
      m < final.lenA
  | state, final, habove, hrun => by
    by_cases hguard : state.lenB ≥ m + 1
    · rcases hgcdIterLoop_step_shape this m Q W3 scratch state final hguard
        hrun with
      ⟨heap1, lenQ, lenR, row23, row01, hdiv, hrow23, hrow01, htail, hlt⟩
      let next : HgcdIterState := {
        heap := row01.heap
        matrix := row01.matrix
        A := state.B
        lenA := state.lenB
        B := state.T
        lenB := lenR
        T := row01.T
        lenT := row01.lenT
        t := row01.t
        sgn := -state.sgn }
      have htail' : hgcdIterLoop this m Q W3 scratch next = .ok final := by
        simpa [next] using htail
      exact hgcdIterLoop_preserves_above_half this m Q W3 scratch next final
        (by simpa [next] using hguard) htail'
    · have hstop : state.lenB < m + 1 := by omega
      have hsame := hgcdIterLoop_stop this m Q W3 scratch state hstop
      have hfinal : state = final := Except.ok.inj (hsame.symm.trans hrun)
      subst final
      exact habove
termination_by state => state.lenB
decreasing_by exact hlt

/-- End-to-end refinement of generated C++ `_hgcd_iter`: the exact identity
initialization and ordered copies feed the exact well-founded loop theorem. -/
theorem hgcdIter_refines (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (M : HgcdMat) (A B T t : RawPtr UInt64) (lenT : Nat)
    (a : RawPtr UInt64) (lenA : Nat) (b : RawPtr UInt64) (lenB : Nat)
    (Q : RawPtr UInt64) (W3 : RawPtr Word3) (scratch : RawPtr UInt64)
    (heap : RawHeap) (final : HgcdIterState)
    (left right : Polynomial (ZMod this._p.toNat)) (hM : M.Valid)
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (horder : lenB ≤ lenA)
    (hlenAPos : 0 < lenA)
    (physical : HgcdLoopWorkspaceProvider this (lenA / 2) Q W3 scratch)
    (h0 : heap.ValidU64Slice (hgcdMatPtr M hM (0 : Fin 4)) 1)
    (h3 : heap.ValidU64Slice (hgcdMatPtr M hM (3 : Fin 4)) 1)
    (h03 : U64SlicesDisjoint (hgcdMatPtr M hM (0 : Fin 4)) 1
      (hgcdMatPtr M hM (3 : Fin 4)) 1)
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB)
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
    (hrun : dense_upoly_zp__hgcd_iter_ir this M A B T t lenT a lenA b
      lenB Q W3 scratch heap = .ok final) :
    ∃ finalA finalB finalEntries hFinalM,
      HgcdIterRawInvariant this left right finalA finalB finalEntries final
        hFinalM ∧
      normalize (EuclideanDomain.gcd left right) =
        normalize (EuclideanDomain.gcd finalA finalB) ∧
      final.lenB < lenA / 2 + 1 ∧
      HgcdMatrixLengthInvariant lenA final hFinalM ∧
      final.lenB ≤ final.lenA ∧ final.lenA ≤ lenA ∧ 0 < final.lenA ∧
      lenA / 2 < final.lenA ∧
      HgcdMatrixCoefficientBound lenA final hFinalM := by
  rcases hgcdIterInit_refines this M A B T t lenT a lenA b lenB heap left
      right hM hp h0 h3 h03 hA hB hAa hBb hAb hBA h0a h3a h0b h3b
      hAMatrix hBMatrix hMatrixValid hLeft hRight with
    ⟨initial, hinit, _, _, _, _, _, _, _, _, hInitialM,
      hInitialMatrix, hInitialA, hInitialB, hInitialTransform,
      hInitialDet⟩
  have hloop : hgcdIterLoop this (lenA / 2) Q W3 scratch initial =
      .ok final := by
    simpa [dense_upoly_zp__hgcd_iter_ir, hinit] using hrun
  have hInitialInvariant : HgcdIterRawInvariant this left right left right
      (identityEntries this._p.toNat) initial hInitialM :=
    ⟨hInitialMatrix, hInitialA, hInitialB, hInitialTransform, hInitialDet⟩
  have hInitialLength := hgcdIterInit_matrixLengthInvariant M A B T t lenT
    a lenA b lenB heap initial hInitialM horder hinit
  have hInitialCoeff := hgcdIterInit_matrixCoefficientBound M A B T t lenT
    a lenA b lenB heap initial hInitialM hlenAPos hinit
  have hInitialOrder : initial.lenB ≤ initial.lenA := by
    have hlens := hgcdIterInit_lengths M A B T t lenT a lenA b lenB heap
      initial hinit
    omega
  have hFinalAbove := hgcdIterLoop_preserves_above_half this (lenA / 2) Q
    W3 scratch initial final (by
      have hlens := hgcdIterInit_lengths M A B T t lenT a lenA b lenB heap
        initial hinit
      omega) hloop
  rcases hgcdIterLoop_refines this (lenA / 2) Q W3 scratch left right hcfg hp
      physical initial final left right (identityEntries this._p.toNat)
      hInitialM hInitialInvariant hloop with
    ⟨finalA, finalB, finalEntries, hFinalM, hFinalRaw, hGcd, hStop⟩
  rcases hgcdIterLoop_preserves_matrixLength this lenA (lenA / 2) Q W3
      scratch left right hcfg hp rfl physical initial final left right
      (identityEntries this._p.toNat) hInitialM hInitialInvariant
      hInitialLength hInitialCoeff hInitialOrder (by
        have hlens := hgcdIterInit_lengths M A B T t lenT a lenA b lenB heap
          initial hinit
        omega) (by
        have hlens := hgcdIterInit_lengths M A B T t lenT a lenA b lenB heap
          initial hinit
        omega) hloop with
    ⟨hFinalM', hFinalLength, hFinalOrder, hFinalInputBound,
      hFinalPositive, hFinalCoeff⟩
  have hhFinal : hFinalM' = hFinalM := Subsingleton.elim _ _
  subst hFinalM'
  exact ⟨finalA, finalB, finalEntries, hFinalM, hFinalRaw, hGcd, hStop,
    hFinalLength, hFinalOrder, hFinalInputBound, hFinalPositive, hFinalAbove,
    hFinalCoeff⟩

/-- Purely physical obligations needed after one concrete iterator result to
stabilize its matrix and normalize its two output pointers.  No polynomial
or expected L2 result occurs in this contract. -/
structure HgcdRecursiveIterFinalizeWorkspace
    (original : HgcdMat) (hOriginal : original.Valid)
    (a3 b3 stage : RawPtr UInt64)
    (iter : HgcdIterState) (hIter : iter.matrix.Valid)
    (stable : HgcdMatRestoreResult) (hStable : stable.matrix.Valid) : Prop where
  stabilize : HgcdMatStabilizeWorkspace iter.heap original iter.matrix
    hOriginal hIter stage
  stageA : stage.region ≠ iter.A.region
  stageB : stage.region ≠ iter.B.region
  originalA : ∀ j : Fin 4,
    (hgcdMatPtr original hOriginal j).region ≠ iter.A.region
  originalB : ∀ j : Fin 4,
    (hgcdMatPtr original hOriginal j).region ≠ iter.B.region
  validA3 : iter.heap.ValidU64Slice a3 iter.lenA
  validB3 : iter.heap.ValidU64Slice b3 iter.lenB
  pAEq : (iter.A == a3) = true → iter.A = a3
  pBEq : (iter.B == b3) = true → iter.B = b3
  b3PB : U64SlicesDisjoint b3 iter.lenB iter.B iter.lenB
  b3PA : U64SlicesDisjoint b3 iter.lenB iter.A iter.lenA
  a3PA : U64SlicesDisjoint a3 iter.lenA iter.A iter.lenA
  a3PB : U64SlicesDisjoint a3 iter.lenA iter.B iter.lenB
  a3B3 : U64SlicesDisjoint a3 iter.lenA b3 iter.lenB
  a3Matrix : ∀ j : Fin 4, U64SlicesDisjoint a3 iter.lenA
    (hgcdMatPtr stable.matrix hStable j) (hgcdMatLen stable.matrix hStable j)
  b3Matrix : ∀ j : Fin 4, U64SlicesDisjoint b3 iter.lenB
    (hgcdMatPtr stable.matrix hStable j) (hgcdMatLen stable.matrix hStable j)

/-- Separation-logic provider for the exact successful iterator and
stabilization executions exposed by `hgcdRecursiveIterBranch_exec`. -/
def HgcdRecursiveIterFinalizeWorkspaceProvider (this : DenseUPolyZp)
    (original : HgcdMat) (hOriginal : original.Valid)
    (a3 b3 inputA inputB : RawPtr UInt64) (lenInputA lenInputB : Nat)
    (Q : RawPtr UInt64) (W3 : RawPtr Word3)
    (T0 T1 scratch stage : RawPtr UInt64) (heap : RawHeap) : Prop :=
  ∀ (iter : HgcdIterState) (hIter : iter.matrix.Valid)
    (stable : HgcdMatRestoreResult) (hStable : stable.matrix.Valid),
    dense_upoly_zp__hgcd_iter_ir this original a3 b3 T0 T1 0 inputA
        lenInputA inputB lenInputB Q W3 scratch heap = .ok iter →
    hgcdMatStabilize original iter.matrix hOriginal hIter stage iter.heap =
        .ok stable →
    HgcdRecursiveIterFinalizeWorkspace original hOriginal a3 b3 stage iter
      hIter stable hStable

/-- End-to-end refinement of an actual recursive-HGCD iterator arm.  The
proof consumes exactly the generated iterator, stabilization, and
alias-sensitive store executions and returns their raw polynomial meaning,
matrix transform, signed determinant, GCD invariant, and stopping bound. -/
theorem hgcdRecursiveIterBranch_refines (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (original : HgcdMat) (hOriginal : original.Valid)
    (a3 b3 inputA inputB : RawPtr UInt64) (lenInputA lenInputB : Nat)
    (Q : RawPtr UInt64) (W3 : RawPtr Word3)
    (T0 T1 scratch stage : RawPtr UInt64) (heap : RawHeap)
    (result : HgcdRecursiveIterBranchResult)
    (left right : Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (hInputOrder : lenInputB ≤ lenInputA)
    (hInputAPos : 0 < lenInputA)
    (loopPhysical : HgcdLoopWorkspaceProvider this (lenInputA / 2) Q W3
      scratch)
    (finalizePhysical : HgcdRecursiveIterFinalizeWorkspaceProvider this
      original hOriginal a3 b3 inputA inputB lenInputA lenInputB Q W3 T0 T1
      scratch stage heap)
    (h0 : heap.ValidU64Slice
      (hgcdMatPtr original hOriginal (0 : Fin 4)) 1)
    (h3 : heap.ValidU64Slice
      (hgcdMatPtr original hOriginal (3 : Fin 4)) 1)
    (h03 : U64SlicesDisjoint
      (hgcdMatPtr original hOriginal (0 : Fin 4)) 1
      (hgcdMatPtr original hOriginal (3 : Fin 4)) 1)
    (hA3 : heap.ValidU64Slice a3 lenInputA)
    (hB3 : heap.ValidU64Slice b3 lenInputB)
    (hAInput : U64SlicesDisjoint a3 lenInputA inputA lenInputA)
    (hBInput : U64SlicesDisjoint b3 lenInputB inputB lenInputB)
    (hAInputB : U64SlicesDisjoint a3 lenInputA inputB lenInputB)
    (hB3A3 : U64SlicesDisjoint b3 lenInputB a3 lenInputA)
    (h0A : U64SlicesDisjoint
      (hgcdMatPtr original hOriginal (0 : Fin 4)) 1 inputA lenInputA)
    (h3A : U64SlicesDisjoint
      (hgcdMatPtr original hOriginal (3 : Fin 4)) 1 inputA lenInputA)
    (h0B : U64SlicesDisjoint
      (hgcdMatPtr original hOriginal (0 : Fin 4)) 1 inputB lenInputB)
    (h3B : U64SlicesDisjoint
      (hgcdMatPtr original hOriginal (3 : Fin 4)) 1 inputB lenInputB)
    (hA3Matrix : ∀ i : Fin 4, U64SlicesDisjoint a3 lenInputA
      (hgcdMatPtr original hOriginal i) (identityEntryLen i))
    (hB3Matrix : ∀ i : Fin 4, U64SlicesDisjoint b3 lenInputB
      (hgcdMatPtr original hOriginal i) (identityEntryLen i))
    (hMatrixValid : ∀ i : Fin 4, heap.ValidU64Slice
      (hgcdMatPtr original hOriginal i) (identityEntryLen i))
    (hLeft : RawDensePolyRep this heap inputA lenInputA left)
    (hRight : RawDensePolyRep this heap inputB lenInputB right)
    (hrun : hgcdRecursiveIterBranch this original hOriginal a3 b3 inputA
      inputB lenInputA lenInputB Q W3 T0 T1 scratch stage heap = .ok result) :
    ∃ finalA finalB finalEntries hResultM,
      HgcdRecursiveRawInvariant this left right finalA finalB finalEntries
        true a3 b3 lenInputA (result.toResult hResultM) := by
  rcases hgcdRecursiveIterBranch_exec this original hOriginal a3 b3 inputA
      inputB lenInputA lenInputB Q W3 T0 T1 scratch stage heap result hrun with
    ⟨iter, hIter, stable, hiter, hstable, hstore, hResultMatrix,
      hResultLenA, hResultLenB, hResultSgn⟩
  rcases hgcdIter_refines this original a3 b3 T0 T1 0 inputA lenInputA
      inputB lenInputB Q W3 scratch heap iter left right hOriginal hcfg hp
      hInputOrder hInputAPos loopPhysical h0 h3 h03 hA3 hB3 hAInput hBInput hAInputB hB3A3 h0A
      h3A h0B h3B hA3Matrix hB3Matrix hMatrixValid hLeft hRight hiter with
    ⟨finalA, finalB, finalEntries, hFinalIter,
      hInvariant, hGcd, hStop, hMatrixLength, hFinalOrder,
      hFinalInputBound, hFinalPositive, hFinalAbove, _hFinalCoeff⟩
  have hStableDescriptors := hgcdMatStabilize_preserves_descriptors original
    iter.matrix hOriginal hIter stage iter.heap stable hstable
  have hStableValid : stable.matrix.Valid := hStableDescriptors.1
  have hFinalize := finalizePhysical iter hIter stable hStableValid hiter hstable
  rcases hgcdMatStabilize_refines this original iter.matrix hOriginal hIter
      stage finalEntries iter.heap hFinalize.stabilize
      hInvariant.matrixRep with
    ⟨semanticStable, hSemanticStable, hSemanticValid, hStableMatrixRep⟩
  have hStableEq : semanticStable = stable :=
    Except.ok.inj (hSemanticStable.symm.trans hstable)
  subst semanticStable
  have hAFrame := hgcdMatStabilize_preserves_rawDenseRep this original
    iter.matrix hOriginal hIter stage finalEntries iter.heap
    hFinalize.stabilize hInvariant.matrixRep iter.A iter.lenA finalA
    hFinalize.stageA hFinalize.originalA hInvariant.aRep stable hstable
  have hBFrame := hgcdMatStabilize_preserves_rawDenseRep this original
    iter.matrix hOriginal hIter stage finalEntries iter.heap
    hFinalize.stabilize hInvariant.matrixRep iter.B iter.lenB finalB
    hFinalize.stageB hFinalize.originalB hInvariant.bRep stable hstable
  have hA3Stable : stable.heap.ValidU64Slice a3 iter.lenA :=
    (hAFrame.1 a3 iter.lenA).mp hFinalize.validA3
  have hB3Stable : stable.heap.ValidU64Slice b3 iter.lenB :=
    (hBFrame.1 b3 iter.lenB).mp hFinalize.validB3
  rcases hgcdRecursiveStoreIterOutputs_refines this a3 b3 iter.A iter.B
      iter.lenA iter.lenB finalA finalB stable.heap hFinalize.pAEq
      hFinalize.pBEq hA3Stable hB3Stable hAFrame.2 hBFrame.2
      hFinalize.b3PB hFinalize.b3PA hFinalize.a3PA hFinalize.a3PB
      hFinalize.a3B3 with
    ⟨semanticHeap, hSemanticStore, hFinalARep, hFinalBRep⟩
  have hHeapEq : semanticHeap = result.heap :=
    Except.ok.inj (hSemanticStore.symm.trans hstore)
  subst semanticHeap
  have hMatrixAfter := hgcdRecursiveStoreIterOutputs_preserves_matrix this
    a3 b3 iter.A iter.B iter.lenA iter.lenB stable.heap result.heap
    stable.matrix finalEntries hSemanticValid hA3Stable hB3Stable
    hAFrame.2.1 hBFrame.2.1 hFinalize.a3Matrix hFinalize.b3Matrix
    hStableMatrixRep hstore
  have hResultValid : result.matrix.Valid := by
    rw [hResultMatrix]
    exact hSemanticValid
  have hResultLength0 : hgcdMatLen result.matrix hResultValid (0 : Fin 4) +
      result.lenA ≤ lenInputA + 1 := by
    simp only [hgcdMatLen, hResultMatrix, hResultLenA,
      hStableDescriptors.2.2]
    exact hMatrixLength.row0A
  have hResultLength1 : hgcdMatLen result.matrix hResultValid (1 : Fin 4) +
      result.lenB ≤ lenInputA + 1 := by
    simp only [hgcdMatLen, hResultMatrix, hResultLenB,
      hStableDescriptors.2.2]
    exact hMatrixLength.row1B
  have hResultLength2 : hgcdMatLen result.matrix hResultValid (2 : Fin 4) +
      result.lenA ≤ lenInputA + 1 := by
    simp only [hgcdMatLen, hResultMatrix, hResultLenA,
      hStableDescriptors.2.2]
    exact hMatrixLength.row2A
  have hResultLength3 : hgcdMatLen result.matrix hResultValid (3 : Fin 4) +
      result.lenB ≤ lenInputA + 1 := by
    simp only [hgcdMatLen, hResultMatrix, hResultLenB,
      hStableDescriptors.2.2]
    exact hMatrixLength.row3B
  have hMatrixResult : HgcdMatRawDenseRep this result.heap result.matrix
      finalEntries hResultValid := by
    simpa only [hResultMatrix] using hMatrixAfter.2
  have hAResult : RawDensePolyRep this result.heap a3 result.lenA finalA := by
    simpa only [hResultLenA] using hFinalARep
  have hBResult : RawDensePolyRep this result.heap b3 result.lenB finalB := by
    simpa only [hResultLenB] using hFinalBRep
  have hResultLengthInvariant : HgcdRecursiveLengthInvariant lenInputA
      (result.toResult hResultValid) := by
    exact {
      row0A := by simpa [HgcdRecursiveIterBranchResult.toResult,
        hgcdMatLenRaw, hgcdMatLen] using hResultLength0
      row1B := by simpa [HgcdRecursiveIterBranchResult.toResult,
        hgcdMatLenRaw, hgcdMatLen] using hResultLength1
      row2A := by simpa [HgcdRecursiveIterBranchResult.toResult,
        hgcdMatLenRaw, hgcdMatLen] using hResultLength2
      row3B := by simpa [HgcdRecursiveIterBranchResult.toResult,
        hgcdMatLenRaw, hgcdMatLen] using hResultLength3
      row1A := by
        simpa [HgcdRecursiveIterBranchResult.toResult, hgcdMatLenRaw,
          hgcdMatLen, hResultMatrix, hStableDescriptors.2.2,
          hResultLenA] using hMatrixLength.row1A
      row3A := by
        simpa [HgcdRecursiveIterBranchResult.toResult, hgcdMatLenRaw,
          hgcdMatLen, hResultMatrix, hStableDescriptors.2.2,
          hResultLenA] using hMatrixLength.row3A
      inputBound := by simpa [HgcdRecursiveIterBranchResult.toResult,
        hResultLenA] using hFinalInputBound
      stopped := by simpa [HgcdRecursiveIterBranchResult.toResult,
        hResultLenB] using hStop
      positive := by simpa [HgcdRecursiveIterBranchResult.toResult,
        hResultLenA] using hFinalPositive
      aboveHalf := by simpa [HgcdRecursiveIterBranchResult.toResult,
        hResultLenA] using hFinalAbove
      coeffBound := by
        intro i
        simpa [HgcdRecursiveIterBranchResult.toResult, hgcdMatLenRaw,
          hgcdMatLen, hResultMatrix, hStableDescriptors.2.2] using
          _hFinalCoeff i }
  refine ⟨finalA, finalB, finalEntries, hResultValid, {
    aRep := by simpa [HgcdRecursiveIterBranchResult.toResult] using hAResult
    bRep := by simpa [HgcdRecursiveIterBranchResult.toResult] using hBResult
    matrixSemantics := by
      intro _
      exact ⟨by simpa [HgcdRecursiveIterBranchResult.toResult] using
          hMatrixResult,
        by simpa [HgcdRecursiveIterBranchResult.toResult, hResultLenA,
          hResultLenB] using hInvariant.transform,
        by simpa [HgcdRecursiveIterBranchResult.toResult, hResultSgn] using
          hInvariant.signedDet⟩
    gcdPreserved := hGcd
    stopped := by
      simpa [HgcdRecursiveIterBranchResult.toResult, hResultLenB] using hStop
    lengths := fun _ => hResultLengthInvariant }⟩

/-- Physical and representation facts needed when a recursive dispatch takes
its generated iterator arm.  The large arm does not inspect these fields;
they are retained so one parent workspace can serve either source branch. -/
structure HgcdRecursiveDispatchIterWorkspace (this : DenseUPolyZp)
    (original : HgcdMat) (hOriginal : original.Valid)
    (a3 b3 inputA inputB : RawPtr UInt64) (lenInputA lenInputB : Nat)
    (Q : RawPtr UInt64) (W3 : RawPtr Word3)
    (T0 T1 scratch stage : RawPtr UInt64) (heap : RawHeap)
    (left right : Polynomial (ZMod this._p.toNat)) : Prop where
  inputAPos : 0 < lenInputA
  loopPhysical : HgcdLoopWorkspaceProvider this (lenInputA / 2) Q W3 scratch
  finalizePhysical : HgcdRecursiveIterFinalizeWorkspaceProvider this original
    hOriginal a3 b3 inputA inputB lenInputA lenInputB Q W3 T0 T1 scratch
    stage heap
  valid0 : heap.ValidU64Slice
    (hgcdMatPtr original hOriginal (0 : Fin 4)) 1
  valid3 : heap.ValidU64Slice
    (hgcdMatPtr original hOriginal (3 : Fin 4)) 1
  disjoint03 : U64SlicesDisjoint
    (hgcdMatPtr original hOriginal (0 : Fin 4)) 1
    (hgcdMatPtr original hOriginal (3 : Fin 4)) 1
  validA3 : heap.ValidU64Slice a3 lenInputA
  validB3 : heap.ValidU64Slice b3 lenInputB
  a3InputA : U64SlicesDisjoint a3 lenInputA inputA lenInputA
  b3InputB : U64SlicesDisjoint b3 lenInputB inputB lenInputB
  a3InputB : U64SlicesDisjoint a3 lenInputA inputB lenInputB
  b3A3 : U64SlicesDisjoint b3 lenInputB a3 lenInputA
  row0InputA : U64SlicesDisjoint
    (hgcdMatPtr original hOriginal (0 : Fin 4)) 1 inputA lenInputA
  row3InputA : U64SlicesDisjoint
    (hgcdMatPtr original hOriginal (3 : Fin 4)) 1 inputA lenInputA
  row0InputB : U64SlicesDisjoint
    (hgcdMatPtr original hOriginal (0 : Fin 4)) 1 inputB lenInputB
  row3InputB : U64SlicesDisjoint
    (hgcdMatPtr original hOriginal (3 : Fin 4)) 1 inputB lenInputB
  a3Matrix : ∀ i : Fin 4, U64SlicesDisjoint a3 lenInputA
    (hgcdMatPtr original hOriginal i) (identityEntryLen i)
  b3Matrix : ∀ i : Fin 4, U64SlicesDisjoint b3 lenInputB
    (hgcdMatPtr original hOriginal i) (identityEntryLen i)
  matrixValid : ∀ i : Fin 4, heap.ValidU64Slice
    (hgcdMatPtr original hOriginal i) (identityEntryLen i)
  leftRep : RawDensePolyRep this heap inputA lenInputA left
  rightRep : RawDensePolyRep this heap inputB lenInputB right

/-- Concrete frame retained across the first cutoff dispatch.  These are the
two source low prefixes consumed by the immediately following reconstruction;
no semantic result is stored in this frame. -/
structure HgcdRecursiveFirstDispatchFrame
    (lowA lowB : RawPtr UInt64) (lenLowA lenLowB : Nat)
    (before : RawHeap) (result : HgcdRecursiveResult) : Prop where
  layout : RawHeap.SameLayout before result.heap
  lowAPrefix : SameU64Prefix before result.heap lowA lenLowA
  lowBPrefix : SameU64Prefix before result.heap lowB lenLowB

/-- Physical frame provider for the actual first dispatch execution. -/
def HgcdRecursiveFirstDispatchFrameProvider
    (this : DenseUPolyZp) (bound : Nat)
    (recurse : HgcdRecursiveCallBelow bound)
    (matrix : HgcdMat) (hMatrix : matrix.Valid)
    (a3 b3 inputA inputB : RawPtr UInt64) (lenInputA lenInputB : Nat)
    (Q : RawPtr UInt64) (W3 : RawPtr Word3)
    (T0 T1 scratch stage WNext : RawPtr UInt64) (heap : RawHeap)
    (horder : lenInputB < lenInputA) (hdecrease : lenInputA < bound)
    (lowA lowB : RawPtr UInt64) (lenLowA lenLowB : Nat) : Prop :=
  ∀ result,
    hgcdRecursiveDispatchBelow this bound recurse matrix hMatrix a3 b3 inputA
      inputB lenInputA lenInputB Q W3 T0 T1 scratch stage WNext heap horder
      hdecrease = .ok result →
    HgcdRecursiveFirstDispatchFrame lowA lowB lenLowA lenLowB heap result

/-- Eliminate the existential order/decrease witnesses of an actual dispatch
result and recover its physical first-call frame. -/
theorem hgcdRecursiveFirstDispatchResult_frame
    (this : DenseUPolyZp) (bound : Nat)
    (recurse : HgcdRecursiveCallBelow bound)
    (matrix : HgcdMat) (hMatrix : matrix.Valid)
    (a3 b3 inputA inputB : RawPtr UInt64) (lenInputA lenInputB : Nat)
    (Q : RawPtr UInt64) (W3 : RawPtr Word3)
    (T0 T1 scratch stage WNext : RawPtr UInt64) (heap : RawHeap)
    (lowA lowB : RawPtr UInt64) (lenLowA lenLowB : Nat)
    (frame : ∀ (horder : lenInputB < lenInputA)
      (hdecrease : lenInputA < bound),
      HgcdRecursiveFirstDispatchFrameProvider this bound recurse matrix
        hMatrix a3 b3 inputA inputB lenInputA lenInputB Q W3 T0 T1 scratch
        stage WNext heap horder hdecrease lowA lowB lenLowA lenLowB)
    (result : HgcdRecursiveResult)
    (hactual : ∃ (horder : lenInputB < lenInputA)
      (hdecrease : lenInputA < bound),
      hgcdRecursiveDispatchBelow this bound recurse matrix hMatrix a3 b3
        inputA inputB lenInputA lenInputB Q W3 T0 T1 scratch stage WNext heap
        horder hdecrease = .ok result) :
    HgcdRecursiveFirstDispatchFrame lowA lowB lenLowA lenLowB heap result := by
  rcases hactual with ⟨horder, hdecrease, hrun⟩
  exact frame horder hdecrease result hrun

/-- The first dispatch frame transports exactly the two low polynomial
prefixes needed by the generated paired reconstruction. -/
theorem hgcdRecursiveFirstDispatch_preserves_low_reps
    (this : DenseUPolyZp) (bound : Nat)
    (recurse : HgcdRecursiveCallBelow bound)
    (matrix : HgcdMat) (hMatrix : matrix.Valid)
    (a3 b3 inputA inputB : RawPtr UInt64) (lenInputA lenInputB : Nat)
    (Q : RawPtr UInt64) (W3 : RawPtr Word3)
    (T0 T1 scratch stage WNext : RawPtr UInt64) (heap : RawHeap)
    (horder : lenInputB < lenInputA) (hdecrease : lenInputA < bound)
    (lowA lowB : RawPtr UInt64) (lenLowA lenLowB : Nat)
    (polyLowA polyLowB : Polynomial (ZMod this._p.toNat))
    (frame : HgcdRecursiveFirstDispatchFrameProvider this bound recurse matrix
      hMatrix a3 b3 inputA inputB lenInputA lenInputB Q W3 T0 T1 scratch
      stage WNext heap horder hdecrease lowA lowB lenLowA lenLowB)
    (hLowA : RawDensePolyRep this heap lowA lenLowA polyLowA)
    (hLowB : RawDensePolyRep this heap lowB lenLowB polyLowB)
    (result : HgcdRecursiveResult)
    (hrun : hgcdRecursiveDispatchBelow this bound recurse matrix hMatrix a3
      b3 inputA inputB lenInputA lenInputB Q W3 T0 T1 scratch stage WNext
      heap horder hdecrease = .ok result) :
    RawDensePolyRep this result.heap lowA lenLowA polyLowA ∧
      RawDensePolyRep this result.heap lowB lenLowB polyLowB := by
  have hframe := frame result hrun
  exact ⟨rawDensePolyRep_of_same_prefix this heap result.heap lowA lenLowA
      polyLowA hframe.layout hframe.lowAPrefix hLowA,
    rawDensePolyRep_of_same_prefix this heap result.heap lowB lenLowB
      polyLowB hframe.layout hframe.lowBPrefix hLowB⟩

/-- Semantic induction hypothesis at one strictly smaller recursive call.
It is phrased over the actual callback execution, not a preselected result. -/
def HgcdRecursiveCallbackRefinesAt (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)] (bound : Nat)
    (recurse : HgcdRecursiveCallBelow bound)
    (matrix : HgcdMat) (hMatrix : matrix.Valid)
    (a3 b3 inputA inputB : RawPtr UInt64) (lenInputA lenInputB : Nat)
    (WNext scratch : RawPtr UInt64) (heap : RawHeap)
    (horder : lenInputB < lenInputA) (hdecrease : lenInputA < bound)
    (left right : Polynomial (ZMod this._p.toNat)) : Prop :=
  ∀ result,
    recurse matrix hMatrix true a3 b3 inputA inputB lenInputA lenInputB
        WNext scratch heap horder hdecrease = .ok result →
    ∃ finalA finalB entries,
      HgcdRecursiveRawInvariant this left right finalA finalB entries true
        a3 b3 lenInputA result

/-- Refinement of the exact cutoff dispatch used at both recursive call
sites.  The small branch is discharged by its generated iterator execution;
the large branch is exactly the well-founded induction hypothesis. -/
theorem hgcdRecursiveDispatchBelow_rawInvariant (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (bound : Nat) (recurse : HgcdRecursiveCallBelow bound)
    (matrix : HgcdMat) (hMatrix : matrix.Valid)
    (a3 b3 inputA inputB : RawPtr UInt64) (lenInputA lenInputB : Nat)
    (Q : RawPtr UInt64) (W3 : RawPtr Word3)
    (T0 T1 scratch stage WNext : RawPtr UInt64) (heap : RawHeap)
    (left right : Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (horder : lenInputB < lenInputA) (hdecrease : lenInputA < bound)
    (iterWorkspace : HgcdRecursiveDispatchIterWorkspace this matrix hMatrix
      a3 b3 inputA inputB lenInputA lenInputB Q W3 T0 T1 scratch stage heap
      left right)
    (recursiveRefines : HgcdRecursiveCallbackRefinesAt this bound recurse
      matrix hMatrix a3 b3 inputA inputB lenInputA lenInputB WNext scratch
      heap horder hdecrease left right)
    (result : HgcdRecursiveResult)
    (hrun : hgcdRecursiveDispatchBelow this bound recurse matrix hMatrix a3
      b3 inputA inputB lenInputA lenInputB Q W3 T0 T1 scratch stage WNext
      heap horder hdecrease = .ok result) :
    ∃ finalA finalB entries,
      HgcdRecursiveRawInvariant this left right finalA finalB entries true
        a3 b3 lenInputA result := by
  rw [hgcdRecursiveDispatchBelow] at hrun
  split at hrun
  next hsmall =>
    split at hrun
    next fault hiter => simp at hrun
    next iter hiter =>
      rcases hgcdRecursiveIterBranch_refines this matrix hMatrix a3 b3 inputA
          inputB lenInputA lenInputB Q W3 T0 T1 scratch stage heap iter left
          right hcfg hp (Nat.le_of_lt horder) iterWorkspace.inputAPos
          iterWorkspace.loopPhysical iterWorkspace.finalizePhysical
          iterWorkspace.valid0 iterWorkspace.valid3 iterWorkspace.disjoint03
          iterWorkspace.validA3 iterWorkspace.validB3 iterWorkspace.a3InputA
          iterWorkspace.b3InputB iterWorkspace.a3InputB iterWorkspace.b3A3
          iterWorkspace.row0InputA iterWorkspace.row3InputA
          iterWorkspace.row0InputB iterWorkspace.row3InputB
          iterWorkspace.a3Matrix iterWorkspace.b3Matrix
          iterWorkspace.matrixValid iterWorkspace.leftRep
          iterWorkspace.rightRep hiter with
        ⟨finalA, finalB, entries, hIterValid, hInvariant⟩
      have heq : iter.toResult hIterValid = result := by
        apply HgcdRecursiveResult.ext_value
        simpa only [HgcdRecursiveResult.value,
          HgcdRecursiveIterBranchResult.toResult] using
          congrArg HgcdRecursiveResult.value (Except.ok.inj hrun)
      subst result
      exact ⟨finalA, finalB, entries, hInvariant⟩
  next hlarge =>
    exact recursiveRefines result hrun

/-- Refinement of a dispatch result whose order/decrease witnesses are packed
existentially by the actual execution relation. -/
theorem hgcdRecursiveDispatchResult_rawInvariant (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (bound : Nat) (recurse : HgcdRecursiveCallBelow bound)
    (matrix : HgcdMat) (hMatrix : matrix.Valid)
    (a3 b3 inputA inputB : RawPtr UInt64) (lenInputA lenInputB : Nat)
    (Q : RawPtr UInt64) (W3 : RawPtr Word3)
    (T0 T1 scratch stage WNext : RawPtr UInt64) (heap : RawHeap)
    (left right : Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (iterWorkspace : HgcdRecursiveDispatchIterWorkspace this matrix hMatrix
      a3 b3 inputA inputB lenInputA lenInputB Q W3 T0 T1 scratch stage heap
      left right)
    (recursiveRefines : ∀ (horder : lenInputB < lenInputA)
      (hdecrease : lenInputA < bound),
      HgcdRecursiveCallbackRefinesAt this bound recurse matrix hMatrix a3 b3
        inputA inputB lenInputA lenInputB WNext scratch heap horder hdecrease
        left right)
    (result : HgcdRecursiveResult)
    (hactual : ∃ (horder : lenInputB < lenInputA)
      (hdecrease : lenInputA < bound),
      hgcdRecursiveDispatchBelow this bound recurse matrix hMatrix a3 b3
        inputA inputB lenInputA lenInputB Q W3 T0 T1 scratch stage WNext heap
        horder hdecrease = .ok result) :
    ∃ finalA finalB entries,
      HgcdRecursiveRawInvariant this left right finalA finalB entries true
        a3 b3 lenInputA result := by
  rcases hactual with ⟨horder, hdecrease, hrun⟩
  exact hgcdRecursiveDispatchBelow_rawInvariant this bound recurse matrix
    hMatrix a3 b3 inputA inputB lenInputA lenInputB Q W3 T0 T1 scratch stage
    WNext heap left right hcfg hp horder hdecrease iterWorkspace
    (recursiveRefines horder hdecrease) result hrun

/-- The exact first-dispatch length provider consumed by the well-founded
body, derived from that same dispatch execution. -/
theorem hgcdRecursiveDispatchBelow_lengthInvariant (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (bound : Nat) (recurse : HgcdRecursiveCallBelow bound)
    (matrix : HgcdMat) (hMatrix : matrix.Valid)
    (a3 b3 inputA inputB : RawPtr UInt64) (lenInputA lenInputB : Nat)
    (Q : RawPtr UInt64) (W3 : RawPtr Word3)
    (T0 T1 scratch stage WNext : RawPtr UInt64) (heap : RawHeap)
    (left right : Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (iterWorkspace : HgcdRecursiveDispatchIterWorkspace this matrix hMatrix
      a3 b3 inputA inputB lenInputA lenInputB Q W3 T0 T1 scratch stage heap
      left right)
    (recursiveRefines : ∀ (horder : lenInputB < lenInputA)
      (hdecrease : lenInputA < bound),
      HgcdRecursiveCallbackRefinesAt this bound recurse matrix hMatrix a3 b3
        inputA inputB lenInputA lenInputB WNext scratch heap horder hdecrease
        left right)
    (result : HgcdRecursiveResult)
    (horder : lenInputB < lenInputA) (hdecrease : lenInputA < bound)
    (hrun : hgcdRecursiveDispatchBelow this bound recurse matrix hMatrix a3
      b3 inputA inputB lenInputA lenInputB Q W3 T0 T1 scratch stage WNext
      heap horder hdecrease = .ok result) :
    HgcdRecursiveLengthInvariant lenInputA result := by
  rcases hgcdRecursiveDispatchBelow_rawInvariant this bound recurse matrix
      hMatrix a3 b3 inputA inputB lenInputA lenInputB Q W3 T0 T1 scratch
      stage WNext heap left right hcfg hp horder hdecrease iterWorkspace
      (recursiveRefines horder hdecrease) result hrun with
    ⟨finalA, finalB, entries, hInvariant⟩
  exact hInvariant.lengths rfl

/-- Frame facts for the exact second cutoff dispatch.  These facts mention
only source buffers that the second call must not overwrite: the first
matrix, the middle quotient, and the two low prefixes later consumed by
`hgcdRecursiveFinish`. -/
structure HgcdRecursiveSecondDispatchFrame
    (preserved : HgcdMat) (hPreserved : preserved.Valid)
    (quotient lowA lowB : RawPtr UInt64)
    (lenQuotient lenLowA lenLowB : Nat)
    (before : RawHeap) (result : HgcdRecursiveResult) : Prop where
  layout : RawHeap.SameLayout before result.heap
  matrixPrefix : ∀ i : Fin 4, SameU64Prefix before result.heap
    (hgcdMatPtr preserved hPreserved i) (hgcdMatLen preserved hPreserved i)
  quotientPrefix : SameU64Prefix before result.heap quotient lenQuotient
  lowAPrefix : SameU64Prefix before result.heap lowA lenLowA
  lowBPrefix : SameU64Prefix before result.heap lowB lenLowB

def HgcdRecursiveSecondDispatchFrameProvider
    (this : DenseUPolyZp) (bound : Nat)
    (recurse : HgcdRecursiveCallBelow bound)
    (matrix : HgcdMat) (hMatrix : matrix.Valid)
    (a3 b3 inputA inputB : RawPtr UInt64) (lenInputA lenInputB : Nat)
    (Q : RawPtr UInt64) (W3 : RawPtr Word3)
    (T0 T1 scratch stage WNext : RawPtr UInt64) (heap : RawHeap)
    (horder : lenInputB < lenInputA) (hdecrease : lenInputA < bound)
    (preserved : HgcdMat) (hPreserved : preserved.Valid)
    (quotient lowA lowB : RawPtr UInt64)
    (lenQuotient lenLowA lenLowB : Nat) : Prop :=
  ∀ result,
    hgcdRecursiveDispatchBelow this bound recurse matrix hMatrix a3 b3 inputA
      inputB lenInputA lenInputB Q W3 T0 T1 scratch stage WNext heap horder
      hdecrease = .ok result →
    HgcdRecursiveSecondDispatchFrame preserved hPreserved quotient lowA lowB
      lenQuotient lenLowA lenLowB heap result

/-- The actual second cutoff dispatch supplies its recursive semantic
invariant while its physical frame transports every value needed by the
following finish block to the returned heap. -/
theorem hgcdRecursiveSecondDispatch_refines (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (bound : Nat) (recurse : HgcdRecursiveCallBelow bound)
    (matrix : HgcdMat) (hMatrix : matrix.Valid)
    (a3 b3 inputA inputB : RawPtr UInt64) (lenInputA lenInputB : Nat)
    (Q : RawPtr UInt64) (W3 : RawPtr Word3)
    (T0 T1 scratch stage WNext : RawPtr UInt64) (heap : RawHeap)
    (left right : Polynomial (ZMod this._p.toNat))
    (preserved : HgcdMat) (hPreserved : preserved.Valid)
    (quotient lowA lowB : RawPtr UInt64)
    (lenQuotient lenLowA lenLowB : Nat)
    (preservedEntries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (quotientPoly lowAPoly lowBPoly : Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (horder : lenInputB < lenInputA) (hdecrease : lenInputA < bound)
    (iterWorkspace : HgcdRecursiveDispatchIterWorkspace this matrix hMatrix
      a3 b3 inputA inputB lenInputA lenInputB Q W3 T0 T1 scratch stage heap
      left right)
    (recursiveRefines : HgcdRecursiveCallbackRefinesAt this bound recurse
      matrix hMatrix a3 b3 inputA inputB lenInputA lenInputB WNext scratch
      heap horder hdecrease left right)
    (frame : HgcdRecursiveSecondDispatchFrameProvider this bound recurse matrix
      hMatrix a3 b3 inputA inputB lenInputA lenInputB Q W3 T0 T1 scratch stage
      WNext heap horder hdecrease preserved hPreserved quotient lowA lowB
      lenQuotient lenLowA lenLowB)
    (hPreservedRep : HgcdMatRawDenseRep this heap preserved preservedEntries
      hPreserved)
    (hQuotient : RawDensePolyRep this heap quotient lenQuotient quotientPoly)
    (hLowA : RawCanonicalPolySlice this heap lowA lenLowA lowAPoly)
    (hLowB : RawCanonicalPolySlice this heap lowB lenLowB lowBPoly)
    (result : HgcdRecursiveResult)
    (hrun : hgcdRecursiveDispatchBelow this bound recurse matrix hMatrix a3
      b3 inputA inputB lenInputA lenInputB Q W3 T0 T1 scratch stage WNext
      heap horder hdecrease = .ok result) :
    ∃ finalA finalB entries,
      HgcdRecursiveRawInvariant this left right finalA finalB entries true
          a3 b3 lenInputA result ∧
        HgcdMatRawDenseRep this result.heap preserved preservedEntries
          hPreserved ∧
        RawDensePolyRep this result.heap quotient lenQuotient quotientPoly ∧
        RawCanonicalPolySlice this result.heap lowA lenLowA lowAPoly ∧
        RawCanonicalPolySlice this result.heap lowB lenLowB lowBPoly := by
  rcases hgcdRecursiveDispatchBelow_rawInvariant this bound recurse matrix
      hMatrix a3 b3 inputA inputB lenInputA lenInputB Q W3 T0 T1 scratch
      stage WNext heap left right hcfg hp horder hdecrease iterWorkspace
      recursiveRefines result hrun with ⟨finalA, finalB, entries, hResult⟩
  have hFrame := frame result hrun
  have hMatrixAfter : HgcdMatRawDenseRep this result.heap preserved
      preservedEntries hPreserved := by
    intro i
    exact rawDensePolyRep_of_same_prefix this heap result.heap
      (hgcdMatPtr preserved hPreserved i) (hgcdMatLen preserved hPreserved i)
      (preservedEntries i) hFrame.layout (hFrame.matrixPrefix i)
      (hPreservedRep i)
  have hQuotientAfter := rawDensePolyRep_of_same_prefix this heap result.heap
    quotient lenQuotient quotientPoly hFrame.layout hFrame.quotientPrefix
    hQuotient
  have hLowAAfter := rawCanonicalPolySlice_of_same_prefix this heap result.heap
    lowA lenLowA lowAPoly hFrame.layout hFrame.lowAPrefix hLowA
  have hLowBAfter := rawCanonicalPolySlice_of_same_prefix this heap result.heap
    lowB lenLowB lowBPoly hFrame.layout hFrame.lowBPrefix hLowB
  exact ⟨finalA, finalB, entries, hResult, hMatrixAfter, hQuotientAfter,
    hLowAAfter, hLowBAfter⟩

/-- Physical obligations for the two concrete guarded products and the
source-selected tail of one `_mat_mul_entry`. -/
structure HgcdMatMulEntryWorkspace
    (heap : RawHeap) (C P Q R S T scratch : RawPtr UInt64)
    (lenP lenQ lenR lenS : Nat)
    (productPQ productRS : HgcdMulTermResult) : Prop where
  first : HgcdMulTermWorkspace heap C P lenP Q lenQ scratch
  second : HgcdMulTermWorkspace productPQ.heap T R lenR S lenS scratch
  firstDstR : U64SlicesDisjoint C (hgcdMulCapacity lenP lenQ) R lenR
  firstScratchR : U64SlicesDisjoint scratch (8 * max lenP lenQ) R lenR
  firstDstS : U64SlicesDisjoint C (hgcdMulCapacity lenP lenQ) S lenS
  firstScratchS : U64SlicesDisjoint scratch (8 * max lenP lenQ) S lenS
  secondDstC : U64SlicesDisjoint T (hgcdMulCapacity lenR lenS)
    C productPQ.length
  secondScratchC : U64SlicesDisjoint scratch (8 * max lenR lenS)
    C productPQ.length
  finalCValid : productRS.heap.ValidU64Slice C
    (max productPQ.length productRS.length)
  addAliasT : ExactOrDisjoint C T
  copyCT : U64SlicesDisjoint C productRS.length T productRS.length

def HgcdMatMulEntryWorkspaceProvider
    (this : DenseUPolyZp) (heap : RawHeap)
    (C P Q R S T scratch : RawPtr UInt64)
    (lenP lenQ lenR lenS : Nat) : Prop :=
  ∀ productPQ productRS,
    hgcdRecursiveMulTerm this C P lenP Q lenQ scratch heap = .ok productPQ →
    hgcdRecursiveMulTerm this T R lenR S lenS scratch productPQ.heap =
      .ok productRS →
    HgcdMatMulEntryWorkspace heap C P Q R S T scratch lenP lenQ lenR lenS
      productPQ productRS

/-- Complete raw semantic refinement of the actual `_mat_mul_entry` control
flow.  All four source tails produce `P*Q + R*S`; zero products are derived
from their real zero-length raw representations. -/
theorem hgcdMatMulEntry_refines (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (C P Q R S T scratch : RawPtr UInt64)
    (lenP lenQ lenR lenS : Nat) (heap : RawHeap)
    (result : HgcdMatMulEntryResult)
    (polyP polyQ polyR polyS : Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (physical : HgcdMatMulEntryWorkspaceProvider this heap C P Q R S T
      scratch lenP lenQ lenR lenS)
    (hP : RawDensePolyRep this heap P lenP polyP)
    (hQ : RawDensePolyRep this heap Q lenQ polyQ)
    (hR : RawDensePolyRep this heap R lenR polyR)
    (hS : RawDensePolyRep this heap S lenS polyS)
    (hrun : hgcdMatMulEntry this C P Q R S T scratch lenP lenQ lenR lenS
      heap = .ok result) :
    RawDensePolyRep this result.heap C result.length
      (polyP * polyQ + polyR * polyS) := by
  rcases hgcdMatMulEntry_exec this C P Q R S T scratch lenP lenQ lenR lenS
      heap result hrun with ⟨productPQ, productRS, hPQ, hRS, htail⟩
  have hwork := physical productPQ productRS hPQ hRS
  rcases hgcdRecursiveMulTerm_refines this C P lenP Q lenQ scratch heap
      polyP polyQ hcfg hp hwork.first hP hQ with
    ⟨semanticPQ, hSemanticPQ, hLayoutPQ, hPQRep⟩
  have hPQEq : semanticPQ = productPQ :=
    Except.ok.inj (hSemanticPQ.symm.trans hPQ)
  subst semanticPQ
  have hRPrefix := hgcdRecursiveMulTerm_preserves_guard this C P lenP Q lenQ
    scratch R lenR heap productPQ hwork.first hP.1 hQ.1 hwork.firstDstR
    hwork.firstScratchR hPQ
  have hSPrefix := hgcdRecursiveMulTerm_preserves_guard this C P lenP Q lenQ
    scratch S lenS heap productPQ hwork.first hP.1 hQ.1 hwork.firstDstS
    hwork.firstScratchS hPQ
  have hR1 := rawDensePolyRep_of_same_prefix this heap productPQ.heap R lenR
    polyR hLayoutPQ hRPrefix hR
  have hS1 := rawDensePolyRep_of_same_prefix this heap productPQ.heap S lenS
    polyS hLayoutPQ hSPrefix hS
  rcases hgcdRecursiveMulTerm_refines this T R lenR S lenS scratch
      productPQ.heap polyR polyS hcfg hp hwork.second hR1 hS1 with
    ⟨semanticRS, hSemanticRS, hLayoutRS, hRSRep⟩
  have hRSEq : semanticRS = productRS :=
    Except.ok.inj (hSemanticRS.symm.trans hRS)
  subst semanticRS
  have hCPrefix := hgcdRecursiveMulTerm_preserves_guard this T R lenR S lenS
    scratch C productPQ.length productPQ.heap productRS hwork.second hR1.1
    hS1.1 hwork.secondDstC hwork.secondScratchC hRS
  have hPQRep2 := rawDensePolyRep_of_same_prefix this productPQ.heap
    productRS.heap C productPQ.length (polyP * polyQ) hLayoutRS hCPrefix hPQRep
  split at htail
  next hboth =>
    rcases htail with ⟨length, hadd, hlength⟩
    subst length
    have hpWord : this._p ≠ 0 := by
      intro hzero
      have hzeroNat := congrArg UInt64.toNat hzero
      simp at hzeroNat
      omega
    exact polyAdd_refines this C C productPQ.length T productRS.length
      productRS.heap result.heap result.length (polyP * polyQ)
      (polyR * polyS) hpWord hwork.finalCValid hPQRep2 hRSRep
      (Or.inl rfl) hwork.addAliasT hadd
  next hnotBoth =>
    split at htail
    next hPQPos =>
      rcases htail with ⟨hheap, hlength⟩
      have hRSLength : productRS.length = 0 := by
        simp at hnotBoth
        omega
      have hzeroRS : polyR * polyS = 0 :=
        slicePolyRep_zero_length productRS.heap T this._p.toNat
          (polyR * polyS) (by simpa [hRSLength] using hRSRep.2.2.1)
      simpa [hheap, hlength, hzeroRS] using hPQRep2
    next hPQZero =>
      have hPQLength : productPQ.length = 0 := by omega
      split at htail
      next hRSPos =>
        rcases htail with ⟨hcopy, hlength⟩
        have hzeroPQ : polyP * polyQ = 0 :=
          slicePolyRep_zero_length productRS.heap C this._p.toNat
            (polyP * polyQ) (by simpa [hPQLength] using hPQRep2.2.2.1)
        have hCValid : productRS.heap.ValidU64Slice C productRS.length :=
          productRS.heap.validU64Slice_mono C
            (max productPQ.length productRS.length) productRS.length
            hwork.finalCValid (by omega)
        rcases copyU64_refines_rawDense this productRS.heap C T
            productRS.length (polyR * polyS) hCValid hwork.copyCT hRSRep with
          ⟨heap', hcopy', _, hrep⟩
        have heq : heap' = result.heap :=
          Except.ok.inj (hcopy'.symm.trans hcopy)
        subst heap'
        simpa [hlength, hzeroPQ] using hrep
      next hRSZero =>
        rcases htail with ⟨hheap, hlength⟩
        have hRSLength : productRS.length = 0 := by omega
        have hzeroPQ : polyP * polyQ = 0 :=
          slicePolyRep_zero_length productRS.heap C this._p.toNat
            (polyP * polyQ) (by simpa [hPQLength] using hPQRep2.2.2.1)
        have hzeroRS : polyR * polyS = 0 :=
          slicePolyRep_zero_length productRS.heap T this._p.toNat
            (polyR * polyS) (by simpa [hRSLength] using hRSRep.2.2.1)
        simpa [hheap, hlength, hzeroPQ, hzeroRS] using
          rawDensePolyRep_zero_length this productRS.heap C
            (productRS.heap.validU64Slice_mono C
              (max productPQ.length productRS.length) 0 hwork.finalCValid
              (by omega))

/-- Purely spatial obligations saying that one concrete `_mat_mul_entry`
does not overwrite an additional live polynomial.  In particular, this
contract contains no polynomial value and no expected matrix product. -/
structure HgcdMatMulEntryGuardWorkspace
    (C P Q R S T scratch guard : RawPtr UInt64)
    (lenP lenQ lenR lenS guardLen : Nat) : Prop where
  firstDst : U64SlicesDisjoint C (hgcdMulCapacity lenP lenQ) guard guardLen
  firstScratch : U64SlicesDisjoint scratch (8 * max lenP lenQ) guard guardLen
  secondDst : U64SlicesDisjoint T (hgcdMulCapacity lenR lenS) guard guardLen
  secondScratch : U64SlicesDisjoint scratch (8 * max lenR lenS) guard guardLen
  addDstRegion : C.region ≠ guard.region

/-- Frame rule for the complete, source-selected `_mat_mul_entry` execution.
The proof follows both guarded products and then the actual add/copy/skip
tail.  It is the bridge needed to retain both input matrices and previously
computed output entries while `_mat_mul` advances. -/
theorem hgcdMatMulEntry_preserves_rawDenseRep (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (C P Q R S T scratch guard : RawPtr UInt64)
    (lenP lenQ lenR lenS guardLen : Nat) (heap : RawHeap)
    (result : HgcdMatMulEntryResult)
    (polyP polyQ polyR polyS guardPoly : Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (physical : HgcdMatMulEntryWorkspaceProvider this heap C P Q R S T
      scratch lenP lenQ lenR lenS)
    (frame : HgcdMatMulEntryGuardWorkspace C P Q R S T scratch guard
      lenP lenQ lenR lenS guardLen)
    (hP : RawDensePolyRep this heap P lenP polyP)
    (hQ : RawDensePolyRep this heap Q lenQ polyQ)
    (hR : RawDensePolyRep this heap R lenR polyR)
    (hS : RawDensePolyRep this heap S lenS polyS)
    (hGuard : RawDensePolyRep this heap guard guardLen guardPoly)
    (hrun : hgcdMatMulEntry this C P Q R S T scratch lenP lenQ lenR lenS
      heap = .ok result) :
    RawHeap.SameLayout heap result.heap ∧
      RawDensePolyRep this result.heap guard guardLen guardPoly := by
  rcases hgcdMatMulEntry_exec this C P Q R S T scratch lenP lenQ lenR lenS
      heap result hrun with ⟨productPQ, productRS, hPQ, hRS, htail⟩
  have hwork := physical productPQ productRS hPQ hRS
  rcases hgcdRecursiveMulTerm_refines this C P lenP Q lenQ scratch heap
      polyP polyQ hcfg hp hwork.first hP hQ with
    ⟨semanticPQ, hSemanticPQ, hLayoutPQ, hPQRep⟩
  have hPQEq : semanticPQ = productPQ :=
    Except.ok.inj (hSemanticPQ.symm.trans hPQ)
  subst semanticPQ
  have hGuardPrefix1 := hgcdRecursiveMulTerm_preserves_guard this C P lenP
    Q lenQ scratch guard guardLen heap productPQ hwork.first hP.1 hQ.1
    frame.firstDst frame.firstScratch hPQ
  have hGuard1 := rawDensePolyRep_of_same_prefix this heap productPQ.heap
    guard guardLen guardPoly hLayoutPQ hGuardPrefix1 hGuard
  have hRPrefix := hgcdRecursiveMulTerm_preserves_guard this C P lenP Q lenQ
    scratch R lenR heap productPQ hwork.first hP.1 hQ.1 hwork.firstDstR
    hwork.firstScratchR hPQ
  have hSPrefix := hgcdRecursiveMulTerm_preserves_guard this C P lenP Q lenQ
    scratch S lenS heap productPQ hwork.first hP.1 hQ.1 hwork.firstDstS
    hwork.firstScratchS hPQ
  have hR1 := rawDensePolyRep_of_same_prefix this heap productPQ.heap R lenR
    polyR hLayoutPQ hRPrefix hR
  have hS1 := rawDensePolyRep_of_same_prefix this heap productPQ.heap S lenS
    polyS hLayoutPQ hSPrefix hS
  rcases hgcdRecursiveMulTerm_refines this T R lenR S lenS scratch
      productPQ.heap polyR polyS hcfg hp hwork.second hR1 hS1 with
    ⟨semanticRS, hSemanticRS, hLayoutRS, hRSRep⟩
  have hRSEq : semanticRS = productRS :=
    Except.ok.inj (hSemanticRS.symm.trans hRS)
  subst semanticRS
  have hGuardPrefix2 := hgcdRecursiveMulTerm_preserves_guard this T R lenR
    S lenS scratch guard guardLen productPQ.heap productRS hwork.second
    hR1.1 hS1.1 frame.secondDst frame.secondScratch hRS
  have hGuard2 := rawDensePolyRep_of_same_prefix this productPQ.heap
    productRS.heap guard guardLen guardPoly hLayoutRS hGuardPrefix2 hGuard1
  have hCPrefix := hgcdRecursiveMulTerm_preserves_guard this T R lenR S lenS
    scratch C productPQ.length productPQ.heap productRS hwork.second hR1.1
    hS1.1 hwork.secondDstC hwork.secondScratchC hRS
  have hPQRep2 := rawDensePolyRep_of_same_prefix this productPQ.heap
    productRS.heap C productPQ.length (polyP * polyQ) hLayoutRS hCPrefix
    hPQRep
  have hLayout12 : RawHeap.SameLayout heap productRS.heap :=
    fun ptr count => (hLayoutPQ ptr count).trans (hLayoutRS ptr count)
  split at htail
  next hboth =>
    rcases htail with ⟨length, hadd, hlength⟩
    subst length
    rcases polyAdd_ok this C C productPQ.length T productRS.length
        productRS.heap hwork.finalCValid hPQRep2.1 hRSRep.1 with
      ⟨heap', length, hadd', hLayoutTail, _⟩
    have heq : (heap', length) = (result.heap, result.length) :=
      Except.ok.inj (hadd'.symm.trans hadd)
    have hHeapEq : heap' = result.heap := congrArg Prod.fst heq
    subst heap'
    have hPrefixTail := polyAdd_preserves_prefix_region_ne this C C
      productPQ.length T guard productRS.length guardLen productRS.heap
      result.heap result.length hwork.finalCValid hPQRep2.1 hRSRep.1
      frame.addDstRegion hadd
    exact ⟨fun ptr count =>
        (hLayout12 ptr count).trans (hLayoutTail ptr count),
      rawDensePolyRep_of_same_prefix this productRS.heap result.heap guard
        guardLen guardPoly hLayoutTail hPrefixTail hGuard2⟩
  next hnotBoth =>
    split at htail
    next hPQPos =>
      rcases htail with ⟨hheap, _⟩
      simpa [hheap] using And.intro hLayout12 hGuard2
    next hPQZero =>
      split at htail
      next hRSPos =>
        rcases htail with ⟨hcopy, _⟩
        have hCValid : productRS.heap.ValidU64Slice C productRS.length :=
          productRS.heap.validU64Slice_mono C
            (max productPQ.length productRS.length) productRS.length
            hwork.finalCValid (by omega)
        have hCopyFrame := copyU64_preserves_rawDenseRep this productRS.heap
          result.heap C T productRS.length guard guardLen guardPoly hCValid
          hRSRep.1
          (by intro _ _ _ _; exact Or.inl frame.addDstRegion) hcopy hGuard2
        exact ⟨fun ptr count =>
            (hLayout12 ptr count).trans (hCopyFrame.1 ptr count),
          hCopyFrame.2⟩
      next hRSZero =>
        rcases htail with ⟨hheap, _⟩
        simpa [hheap] using And.intro hLayout12 hGuard2

/-- The L2 entry selected by source iteration `i` of `_mat_mul`. -/
noncomputable def hgcdMatProductEntry {p : Nat}
    (left right : Fin 4 → Polynomial (ZMod p)) (i : Fin 4) :
    Polynomial (ZMod p) :=
  left ⟨2 * (i.val / 2), by omega⟩ * right ⟨i.val % 2, by omega⟩ +
    left ⟨2 * (i.val / 2) + 1, by omega⟩ *
      right ⟨2 + i.val % 2, by omega⟩

/-- Physical obligations for one actual source iteration of `_mat_mul`.
The frame fields mention only addresses and lengths; the expected L2 matrix
product is deliberately absent. -/
structure HgcdMatMulStepWorkspace (this : DenseUPolyZp)
    (A B C : HgcdMat) (hA : A.Valid) (hB : B.Valid) (hC : C.Valid)
    (T scratch : RawPtr UInt64) (i : Nat) (hi : i < 4)
    (heap : RawHeap) : Prop where
  entry : HgcdMatMulEntryWorkspaceProvider this heap
    (hgcdMatPtr C hC ⟨i, hi⟩)
    (hgcdMatPtr A hA ⟨2 * (i / 2), by omega⟩)
    (hgcdMatPtr B hB ⟨i % 2, by omega⟩)
    (hgcdMatPtr A hA ⟨2 * (i / 2) + 1, by omega⟩)
    (hgcdMatPtr B hB ⟨2 + i % 2, by omega⟩) T scratch
    (hgcdMatLen A hA ⟨2 * (i / 2), by omega⟩)
    (hgcdMatLen B hB ⟨i % 2, by omega⟩)
    (hgcdMatLen A hA ⟨2 * (i / 2) + 1, by omega⟩)
    (hgcdMatLen B hB ⟨2 + i % 2, by omega⟩)
  frameA : ∀ j : Fin 4, HgcdMatMulEntryGuardWorkspace
    (hgcdMatPtr C hC ⟨i, hi⟩)
    (hgcdMatPtr A hA ⟨2 * (i / 2), by omega⟩)
    (hgcdMatPtr B hB ⟨i % 2, by omega⟩)
    (hgcdMatPtr A hA ⟨2 * (i / 2) + 1, by omega⟩)
    (hgcdMatPtr B hB ⟨2 + i % 2, by omega⟩) T scratch
    (hgcdMatPtr A hA j)
    (hgcdMatLen A hA ⟨2 * (i / 2), by omega⟩)
    (hgcdMatLen B hB ⟨i % 2, by omega⟩)
    (hgcdMatLen A hA ⟨2 * (i / 2) + 1, by omega⟩)
    (hgcdMatLen B hB ⟨2 + i % 2, by omega⟩)
    (hgcdMatLen A hA j)
  frameB : ∀ j : Fin 4, HgcdMatMulEntryGuardWorkspace
    (hgcdMatPtr C hC ⟨i, hi⟩)
    (hgcdMatPtr A hA ⟨2 * (i / 2), by omega⟩)
    (hgcdMatPtr B hB ⟨i % 2, by omega⟩)
    (hgcdMatPtr A hA ⟨2 * (i / 2) + 1, by omega⟩)
    (hgcdMatPtr B hB ⟨2 + i % 2, by omega⟩) T scratch
    (hgcdMatPtr B hB j)
    (hgcdMatLen A hA ⟨2 * (i / 2), by omega⟩)
    (hgcdMatLen B hB ⟨i % 2, by omega⟩)
    (hgcdMatLen A hA ⟨2 * (i / 2) + 1, by omega⟩)
    (hgcdMatLen B hB ⟨2 + i % 2, by omega⟩)
    (hgcdMatLen B hB j)
  frameDone : ∀ j : Fin 4, j.val < i → HgcdMatMulEntryGuardWorkspace
    (hgcdMatPtr C hC ⟨i, hi⟩)
    (hgcdMatPtr A hA ⟨2 * (i / 2), by omega⟩)
    (hgcdMatPtr B hB ⟨i % 2, by omega⟩)
    (hgcdMatPtr A hA ⟨2 * (i / 2) + 1, by omega⟩)
    (hgcdMatPtr B hB ⟨2 + i % 2, by omega⟩) T scratch
    (hgcdMatPtr C hC j)
    (hgcdMatLen A hA ⟨2 * (i / 2), by omega⟩)
    (hgcdMatLen B hB ⟨i % 2, by omega⟩)
    (hgcdMatLen A hA ⟨2 * (i / 2) + 1, by omega⟩)
    (hgcdMatLen B hB ⟨2 + i % 2, by omega⟩)
    (hgcdMatLen C hC j)

/-- A provider for every concrete state reached by the four-entry loop. -/
def HgcdMatMulLoopWorkspaceProvider (this : DenseUPolyZp)
    (A B : HgcdMat) (hA : A.Valid) (hB : B.Valid)
    (T scratch : RawPtr UInt64) : Prop :=
  ∀ (C : HgcdMat) (hC : C.Valid) (i : Nat) (hi : i < 4)
    (heap : RawHeap),
    HgcdMatMulStepWorkspace this A B C hA hB hC T scratch i hi heap

/-- Semantic refinement of the actual four-step `_mat_mul` loop.  The
induction follows `hgcdMatMulLoop`; each new result comes from executing
`hgcdMatMulEntry`, while frame proofs retain both inputs and earlier outputs. -/
theorem hgcdMatMulLoop_refines (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (A B : HgcdMat) (hA : A.Valid) (hB : B.Valid)
    (T scratch : RawPtr UInt64) (C : HgcdMat) (hC : C.Valid)
    (i : Nat) (heap : RawHeap) (result : HgcdMatMulResult)
    (left right : Fin 4 → Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (physical : HgcdMatMulLoopWorkspaceProvider this A B hA hB T scratch)
    (hLeft : HgcdMatRawDenseRep this heap A left hA)
    (hRight : HgcdMatRawDenseRep this heap B right hB)
    (hDone : ∀ j : Fin 4, j.val < i → RawDensePolyRep this heap
      (hgcdMatPtr C hC j) (hgcdMatLen C hC j)
      (hgcdMatProductEntry left right j))
    (hrun : hgcdMatMulLoop this A B hA hB T scratch C hC i heap =
      .ok result) :
    HgcdMatRawDenseRep this result.heap A left hA ∧
      HgcdMatRawDenseRep this result.heap B right hB ∧
      ∃ hResult : result.matrix.Valid,
        HgcdMatRawDenseRep this result.heap result.matrix
          (hgcdMatProductEntry left right) hResult := by
  rw [hgcdMatMulLoop] at hrun
  split at hrun
  next hi =>
    dsimp only at hrun
    split at hrun
    next fault hentry => simp at hrun
    next entry hentry =>
      let out : Fin 4 := ⟨i, hi⟩
      let rowBase : Fin 4 := ⟨2 * (i / 2), by omega⟩
      let rowNext : Fin 4 := ⟨2 * (i / 2) + 1, by omega⟩
      let col : Fin 4 := ⟨i % 2, by omega⟩
      let lowerCol : Fin 4 := ⟨2 + i % 2, by omega⟩
      have hstep := physical C hC i hi heap
      have hCurrent := hgcdMatMulEntry_refines this
        (hgcdMatPtr C hC out) (hgcdMatPtr A hA rowBase)
        (hgcdMatPtr B hB col) (hgcdMatPtr A hA rowNext)
        (hgcdMatPtr B hB lowerCol) T scratch
        (hgcdMatLen A hA rowBase) (hgcdMatLen B hB col)
        (hgcdMatLen A hA rowNext) (hgcdMatLen B hB lowerCol) heap entry
        (left rowBase) (right col) (left rowNext) (right lowerCol)
        hcfg hp hstep.entry (hLeft rowBase) (hRight col) (hLeft rowNext)
        (hRight lowerCol) hentry
      have hLeft1 : HgcdMatRawDenseRep this entry.heap A left hA := by
        intro j
        exact (hgcdMatMulEntry_preserves_rawDenseRep this
          (hgcdMatPtr C hC out) (hgcdMatPtr A hA rowBase)
          (hgcdMatPtr B hB col) (hgcdMatPtr A hA rowNext)
          (hgcdMatPtr B hB lowerCol) T scratch (hgcdMatPtr A hA j)
          (hgcdMatLen A hA rowBase) (hgcdMatLen B hB col)
          (hgcdMatLen A hA rowNext) (hgcdMatLen B hB lowerCol)
          (hgcdMatLen A hA j) heap entry (left rowBase) (right col)
          (left rowNext) (right lowerCol) (left j) hcfg hp hstep.entry
          (hstep.frameA j) (hLeft rowBase) (hRight col) (hLeft rowNext)
          (hRight lowerCol) (hLeft j) hentry).2
      have hRight1 : HgcdMatRawDenseRep this entry.heap B right hB := by
        intro j
        exact (hgcdMatMulEntry_preserves_rawDenseRep this
          (hgcdMatPtr C hC out) (hgcdMatPtr A hA rowBase)
          (hgcdMatPtr B hB col) (hgcdMatPtr A hA rowNext)
          (hgcdMatPtr B hB lowerCol) T scratch (hgcdMatPtr B hB j)
          (hgcdMatLen A hA rowBase) (hgcdMatLen B hB col)
          (hgcdMatLen A hA rowNext) (hgcdMatLen B hB lowerCol)
          (hgcdMatLen B hB j) heap entry (left rowBase) (right col)
          (left rowNext) (right lowerCol) (right j) hcfg hp hstep.entry
          (hstep.frameB j) (hLeft rowBase) (hRight col) (hLeft rowNext)
          (hRight lowerCol) (hRight j) hentry).2
      let nextLen := C.len.set i entry.length
        (by rw [hC.2]; exact hi)
      let next : HgcdMat := { C with len := nextLen }
      have hNext : next.Valid := by
        exact ⟨hC.1, by simp [next, nextLen, hC.2]⟩
      have hDone1 : ∀ j : Fin 4, j.val < i + 1 →
          RawDensePolyRep this entry.heap (hgcdMatPtr next hNext j)
            (hgcdMatLen next hNext j) (hgcdMatProductEntry left right j) := by
        intro j hj
        by_cases hji : j = out
        · subst j
          simpa [out, rowBase, rowNext, col, lowerCol, next, nextLen,
            hgcdMatProductEntry, hgcdMatPtr, hgcdMatLen] using hCurrent
        · have hjne : j.val ≠ i := by
            intro hval
            exact hji (Fin.ext hval)
          have hjlt : j.val < i := by omega
          have hpreserved := (hgcdMatMulEntry_preserves_rawDenseRep this
            (hgcdMatPtr C hC out) (hgcdMatPtr A hA rowBase)
            (hgcdMatPtr B hB col) (hgcdMatPtr A hA rowNext)
            (hgcdMatPtr B hB lowerCol) T scratch (hgcdMatPtr C hC j)
            (hgcdMatLen A hA rowBase) (hgcdMatLen B hB col)
            (hgcdMatLen A hA rowNext) (hgcdMatLen B hB lowerCol)
            (hgcdMatLen C hC j) heap entry (left rowBase) (right col)
            (left rowNext) (right lowerCol) (hgcdMatProductEntry left right j)
            hcfg hp hstep.entry (hstep.frameDone j hjlt) (hLeft rowBase)
            (hRight col) (hLeft rowNext) (hRight lowerCol) (hDone j hjlt)
            hentry).2
          simp only [next, nextLen, hgcdMatPtr, hgcdMatLen]
          rw [Array.getElem_set_ne
            (by simpa [hC.2] using hi)
            (by simpa [hC.2] using j.isLt) (Ne.symm hjne)]
          exact hpreserved
      exact hgcdMatMulLoop_refines this A B hA hB T scratch next hNext
        (i + 1) entry.heap result left right hcfg hp physical hLeft1 hRight1
        hDone1 hrun
  next hi =>
    have heq : result = HgcdMatMulResult.mk heap C :=
      (Except.ok.inj hrun).symm
    subst result
    exact ⟨hLeft, hRight, hC, fun j => hDone j (by omega)⟩
termination_by 4 - i
decreasing_by omega

/-- Entry-point theorem for the complete generated `_mat_mul`. -/
theorem hgcdMatMul_refines (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (C A B : HgcdMat) (hC : C.Valid) (hA : A.Valid) (hB : B.Valid)
    (T scratch : RawPtr UInt64) (heap : RawHeap) (result : HgcdMatMulResult)
    (left right : Fin 4 → Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (physical : HgcdMatMulLoopWorkspaceProvider this A B hA hB T scratch)
    (hLeft : HgcdMatRawDenseRep this heap A left hA)
    (hRight : HgcdMatRawDenseRep this heap B right hB)
    (hrun : hgcdMatMul this C A B hC hA hB T scratch heap = .ok result) :
    HgcdMatRawDenseRep this result.heap A left hA ∧
      HgcdMatRawDenseRep this result.heap B right hB ∧
      ∃ hResult : result.matrix.Valid,
        HgcdMatRawDenseRep this result.heap result.matrix
          (hgcdMatProductEntry left right) hResult := by
  exact hgcdMatMulLoop_refines this A B hA hB T scratch C hC 0 heap result
    left right hcfg hp physical hLeft hRight (by intro j hj; omega) hrun

/-- Every descriptor produced by the real four-call `_mat_mul` is bounded by
the two source product capacities for that entry.  Quantifying over `i`
simultaneously covers all four physical output descriptors. -/
theorem hgcdMatProductEntry_length_le (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (heap : RawHeap) (A B C : HgcdMat)
    (hA : A.Valid) (hB : B.Valid) (hC : C.Valid)
    (left right : Fin 4 → Polynomial (ZMod this._p.toNat))
    (hLeft : HgcdMatRawDenseRep this heap A left hA)
    (hRight : HgcdMatRawDenseRep this heap B right hB)
    (hProduct : HgcdMatRawDenseRep this heap C
      (hgcdMatProductEntry left right) hC) (i : Fin 4) :
    hgcdMatLen C hC i ≤
      max
        (hgcdMatLen A hA ⟨2 * (i.val / 2), by omega⟩ +
          hgcdMatLen B hB ⟨i.val % 2, by omega⟩ - 1)
        (hgcdMatLen A hA ⟨2 * (i.val / 2) + 1, by omega⟩ +
          hgcdMatLen B hB ⟨2 + i.val % 2, by omega⟩ - 1) := by
  simpa only [hgcdMatProductEntry] using
    rawDensePolyRep_sum_products_length_le this heap
      (hgcdMatPtr A hA ⟨2 * (i.val / 2), by omega⟩)
      (hgcdMatPtr B hB ⟨i.val % 2, by omega⟩)
      (hgcdMatPtr A hA ⟨2 * (i.val / 2) + 1, by omega⟩)
      (hgcdMatPtr B hB ⟨2 + i.val % 2, by omega⟩)
      (hgcdMatPtr C hC i)
      (hgcdMatLen A hA ⟨2 * (i.val / 2), by omega⟩)
      (hgcdMatLen B hB ⟨i.val % 2, by omega⟩)
      (hgcdMatLen A hA ⟨2 * (i.val / 2) + 1, by omega⟩)
      (hgcdMatLen B hB ⟨2 + i.val % 2, by omega⟩)
      (hgcdMatLen C hC i)
      (left ⟨2 * (i.val / 2), by omega⟩)
      (right ⟨i.val % 2, by omega⟩)
      (left ⟨2 * (i.val / 2) + 1, by omega⟩)
      (right ⟨2 + i.val % 2, by omega⟩)
      (hLeft ⟨2 * (i.val / 2), by omega⟩)
      (hRight ⟨i.val % 2, by omega⟩)
      (hLeft ⟨2 * (i.val / 2) + 1, by omega⟩)
      (hRight ⟨2 + i.val % 2, by omega⟩)
      (hProduct i)

/-- L2 effect of one quotient-matrix column update. -/
noncomputable def hgcdMatQuotientUpdateEntries {p : Nat}
    (entries : Fin 4 → Polynomial (ZMod p))
    (quotient : Polynomial (ZMod p)) (top bottom : Fin 4) :
    Fin 4 → Polynomial (ZMod p) :=
  Function.update entries top (entries top + quotient * entries bottom)

/-- The quotient-entry lowering has exactly the same physical multiplication
and addition requirements as a row update with `bottom` as the multiplied
entry and the existing `top` buffer as the in-place sum destination. -/
abbrev HgcdMatQuotientEntryWorkspace (S : HgcdMat) (hS : S.Valid)
    (top bottom : Fin 4) (q : RawPtr UInt64) (lenQ : Nat)
    (T scratch : RawPtr UInt64) (heap : RawHeap) : Prop :=
  MatRowUpdateWorkspace S bottom top q lenQ T (hgcdMatPtr S hS top)
    scratch heap hS

/-- Complete raw refinement of one actual guarded quotient-column update.
The inactive branch derives the vanishing product from a real zero-length
representation.  The active branch consumes the generated ordered `_mul`
and in-place `_poly_add`, then frames every non-target matrix entry. -/
theorem hgcdMatQuotientEntry_refines (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (S : HgcdMat) (hS : S.Valid) (top bottom : Fin 4)
    (q : RawPtr UInt64) (lenQ : Nat) (T scratch : RawPtr UInt64)
    (heap : RawHeap) (result : HgcdMatQuotientEntryResult)
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (quotient : Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (workspace : HgcdMatQuotientEntryWorkspace S hS top bottom q lenQ T
      scratch heap)
    (hQ : RawDensePolyRep this heap q lenQ quotient)
    (hMatrix : HgcdMatRawDenseRep this heap S entries hS)
    (hrun : hgcdMatQuotientEntry this S hS top bottom q lenQ T scratch
      heap = .ok result) :
    HgcdMatRawDenseRep this result.heap result.matrix
      (hgcdMatQuotientUpdateEntries entries quotient top bottom)
      result.valid ∧
    RawDensePolyRep this result.heap q lenQ quotient := by
  simp only [hgcdMatQuotientEntry] at hrun
  split at hrun
  next hactive =>
    have hparts : decide (lenQ > 0) = true ∧
        decide (hgcdMatLenRaw S hS bottom > 0) = true := by
      simpa using hactive
    have hQPos : 0 < lenQ := by
      simpa only [decide_eq_true_eq] using hparts.1
    have hBottomPos : 0 < hgcdMatLen S hS bottom := by
      simpa only [decide_eq_true_eq, hgcdMatLenRaw, hgcdMatLen] using
        hparts.2
    generalize hmul : Generated.StrictMul.dense_upoly_zp__mul_ir this T
        (if lenQ ≥ hgcdMatLenRaw S hS bottom then q
          else hgcdMatPtrRaw S hS bottom)
        (if lenQ ≥ hgcdMatLenRaw S hS bottom then lenQ
          else hgcdMatLenRaw S hS bottom)
        (if lenQ ≥ hgcdMatLenRaw S hS bottom then
          hgcdMatPtrRaw S hS bottom else q)
        (if lenQ ≥ hgcdMatLenRaw S hS bottom then
          hgcdMatLenRaw S hS bottom else lenQ)
        scratch heap = mulRun at hrun
    cases mulRun with
    | error fault => simp at hrun
    | ok heap1 =>
      simp only [hmul] at hrun
      generalize hadd : Generated.StrictPolyAddSub.dense_upoly_zp__poly_add_ir this
          (hgcdMatPtrRaw S hS top) (hgcdMatPtrRaw S hS top)
          (hgcdMatLenRaw S hS top) T
          (lenQ + hgcdMatLenRaw S hS bottom - 1) heap1 = addRun at hrun
      cases addRun with
      | error fault => simp [hadd] at hrun
      | ok pair =>
        rcases pair with ⟨heap2, sumLen⟩
        simp only [hadd] at hrun
        have heq := Except.ok.inj hrun
        cases heq
        have hProduct := matRowUpdate_mul_result this S hS bottom q T scratch
          lenQ heap heap1 quotient (entries bottom) hcfg hp hQPos hBottomPos
          workspace.lenWord workspace.validT workspace.validScratch
          workspace.disjointTQ (workspace.disjointTMatrix bottom)
          workspace.disjointTScratch workspace.disjointScratchQ
          (workspace.disjointScratchMatrix bottom) hQ (hMatrix bottom)
          (by simpa [hgcdMatPtrRaw, hgcdMatLenRaw, hgcdMatPtr,
              hgcdMatLen] using hmul)
        have hMatrix1 : HgcdMatRawDenseRep this heap1 S entries hS := by
          intro j
          exact matRowUpdate_mul_preserves_guard this S hS bottom q T scratch
            (hgcdMatPtr S hS j) lenQ (hgcdMatLen S hS j) heap heap1
            (entries j) hQPos hBottomPos workspace.validT
            workspace.validScratch hQ.1 (hMatrix bottom).1
            workspace.disjointScratchQ
            (workspace.disjointScratchMatrix bottom)
            (workspace.disjointTMatrix j)
            (workspace.disjointScratchMatrix j) (hMatrix j) hProduct.1
            (by simpa [hgcdMatPtrRaw, hgcdMatLenRaw, hgcdMatPtr,
                hgcdMatLen] using hmul)
        have hQ1 := matRowUpdate_mul_preserves_guard this S hS bottom q T
          scratch q lenQ lenQ heap heap1 quotient hQPos hBottomPos
          workspace.validT workspace.validScratch hQ.1 (hMatrix bottom).1
          workspace.disjointScratchQ
          (workspace.disjointScratchMatrix bottom) workspace.disjointTQ
          workspace.disjointScratchQ hQ hProduct.1
          (by simpa [hgcdMatPtrRaw, hgcdMatLenRaw, hgcdMatPtr,
              hgcdMatLen] using hmul)
        have hOutput1 := (hProduct.1 _ _).mp workspace.validAddOutput
        have hpWord : this._p ≠ 0 := by
          intro hzero
          have hzeroNat := congrArg UInt64.toNat hzero
          simp at hzeroNat
          omega
        have hSum : RawDensePolyRep this heap2 (hgcdMatPtr S hS top)
            sumLen (entries top + quotient * entries bottom) :=
          polyAdd_refines this (hgcdMatPtr S hS top)
            (hgcdMatPtr S hS top) (hgcdMatLen S hS top) T
            (lenQ + hgcdMatLen S hS bottom - 1) heap1 heap2 sumLen
            (entries top) (quotient * entries bottom) hpWord hOutput1
            (hMatrix1 top) hProduct.2 workspace.aliasEntry1
            workspace.aliasProduct
            (by simpa [hgcdMatPtrRaw, hgcdMatLenRaw, hgcdMatPtr,
                hgcdMatLen] using hadd)
        rcases polyAdd_ok this (hgcdMatPtr S hS top)
            (hgcdMatPtr S hS top) (hgcdMatLen S hS top) T
            (lenQ + hgcdMatLen S hS bottom - 1) heap1 hOutput1
            (hMatrix1 top).1 hProduct.2.1 with
          ⟨heapAdd, lenAdd, haddAdd, hLayoutAdd, _⟩
        have heqAdd : (heapAdd, lenAdd) = (heap2, sumLen) :=
          Except.ok.inj (haddAdd.symm.trans (by
            simpa [hgcdMatPtrRaw, hgcdMatLenRaw, hgcdMatPtr,
              hgcdMatLen] using hadd))
        have hHeapAdd : heapAdd = heap2 := congrArg Prod.fst heqAdd
        have hLayoutAdd2 : RawHeap.SameLayout heap1 heap2 := by
          intro ptr count
          rw [← hHeapAdd]
          exact hLayoutAdd ptr count
        have hQPrefix := polyAdd_preserves_prefix_region_ne this
          (hgcdMatPtr S hS top) (hgcdMatPtr S hS top)
          (hgcdMatLen S hS top) T q
          (lenQ + hgcdMatLen S hS bottom - 1) lenQ heap1 heap2 sumLen
          hOutput1 (hMatrix1 top).1 hProduct.2.1 workspace.tDisjointQ
          (by simpa [hgcdMatPtrRaw, hgcdMatLenRaw, hgcdMatPtr,
              hgcdMatLen] using hadd)
        have hQ2 := rawDensePolyRep_of_same_prefix this heap1 heap2 q lenQ
          quotient hLayoutAdd2 hQPrefix hQ1
        constructor
        · intro j
          by_cases hj : j = top
          · subst j
            simpa [hgcdMatQuotientUpdateEntries, hgcdMatPtr, hgcdMatLen,
              hS.2] using hSum
          · have hPrefix := polyAdd_preserves_prefix_region_ne this
              (hgcdMatPtr S hS top) (hgcdMatPtr S hS top)
              (hgcdMatLen S hS top) T (hgcdMatPtr S hS j)
              (lenQ + hgcdMatLen S hS bottom - 1) (hgcdMatLen S hS j)
              heap1 heap2 sumLen hOutput1 (hMatrix1 top).1 hProduct.2.1
              (workspace.tDisjointOther j hj)
              (by simpa [hgcdMatPtrRaw, hgcdMatLenRaw, hgcdMatPtr,
                  hgcdMatLen] using hadd)
            have hFrame := rawDensePolyRep_of_same_prefix this heap1 heap2
              (hgcdMatPtr S hS j) (hgcdMatLen S hS j) (entries j)
              hLayoutAdd2 hPrefix (hMatrix1 j)
            simp only [hgcdMatPtr, hgcdMatLen]
            rw [Array.getElem_set_ne
              (by simpa [hS.2] using top.isLt)
              (by simpa [hS.2] using j.isLt) (by
                intro hval
                exact hj (Fin.ext hval.symm))]
            simpa [hgcdMatQuotientUpdateEntries, hj] using hFrame
        · simpa using hQ2
  next hinactive =>
    have heq := Except.ok.inj hrun
    subst result
    have hzero : quotient = 0 ∨ entries bottom = 0 := by
      have hnot : ¬(0 < lenQ ∧ 0 < hgcdMatLen S hS bottom) := by
        simpa only [Bool.and_eq_true, decide_eq_true_eq, hgcdMatLenRaw,
          hgcdMatLen] using hinactive
      rcases not_and_or.mp hnot with hq | hb
      · left
        have hlen : lenQ = 0 := by omega
        exact slicePolyRep_zero_length heap q this._p.toNat quotient
          (by simpa [hlen] using hQ.2.2.1)
      · right
        have hlen : hgcdMatLen S hS bottom = 0 := by omega
        exact slicePolyRep_zero_length heap (hgcdMatPtr S hS bottom)
          this._p.toNat (entries bottom)
          (by simpa [hlen] using (hMatrix bottom).2.2.1)
    constructor
    · intro j
      by_cases hj : j = top
      · subst j
        rcases hzero with hq | hb
        · simpa [hgcdMatQuotientUpdateEntries, hq] using hMatrix top
        · simpa [hgcdMatQuotientUpdateEntries, hb] using hMatrix top
      · simpa [hgcdMatQuotientUpdateEntries, hj] using hMatrix j
    · exact hQ

/-- L2 descriptor permutation induced by the four source row swaps. -/
noncomputable def hgcdMatSwapEntries {p : Nat}
    (entries : Fin 4 → Polynomial (ZMod p)) :
    Fin 4 → Polynomial (ZMod p) :=
  fun j => entries ⟨(j.val + 2) % 4, by omega⟩

/-- The descriptor-only row swaps preserve the same four raw buffers and
only permute their L2 indexing. -/
theorem hgcdMatSwapRows_refines (this : DenseUPolyZp)
    (S : HgcdMat) (hS : S.Valid)
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat)) (heap : RawHeap)
    (hMatrix : HgcdMatRawDenseRep this heap S entries hS) :
    HgcdMatRawDenseRep this heap (hgcdMatSwapRows S hS)
      (hgcdMatSwapEntries entries) (hgcdMatSwapRows_valid S hS) := by
  intro j
  fin_cases j <;>
    simpa [hgcdMatSwapRows, hgcdMatSwapEntries, hgcdMatPtr, hgcdMatLen,
      hgcdMatPtrRaw, hgcdMatLenRaw] using
      hMatrix (⟨_, by omega⟩ : Fin 4)

/-- Final L2 entries after both source-order quotient column updates. -/
noncomputable def hgcdMatApplyQuotientEntries {p : Nat}
    (entries : Fin 4 → Polynomial (ZMod p))
    (quotient : Polynomial (ZMod p)) : Fin 4 → Polynomial (ZMod p) :=
  let swapped := hgcdMatSwapEntries entries
  let first := hgcdMatQuotientUpdateEntries swapped quotient
    (0 : Fin 4) (2 : Fin 4)
  hgcdMatQuotientUpdateEntries first quotient (1 : Fin 4) (3 : Fin 4)

/-- The quotient update represents the source matrix
`[[q,1],[1,0]] * S`, whose determinant is the negation of `det S`. -/
theorem hgcdMatApplyQuotientEntries_det
    {p : Nat} (entries : Fin 4 → Polynomial (ZMod p))
    (quotient : Polynomial (ZMod p)) :
    hgcdMatApplyQuotientEntries entries quotient 0 *
        hgcdMatApplyQuotientEntries entries quotient 3 -
      hgcdMatApplyQuotientEntries entries quotient 1 *
        hgcdMatApplyQuotientEntries entries quotient 2 =
      -(entries 0 * entries 3 - entries 1 * entries 2) := by
  simp [hgcdMatApplyQuotientEntries, hgcdMatSwapEntries,
    hgcdMatQuotientUpdateEntries, Function.update]
  ring

/-- Determinant of the exact flat-matrix product used by `_mat_mul`. -/
theorem hgcdMatProductEntry_det
    {p : Nat} (left right : Fin 4 → Polynomial (ZMod p)) :
    hgcdMatProductEntry left right 0 * hgcdMatProductEntry left right 3 -
      hgcdMatProductEntry left right 1 * hgcdMatProductEntry left right 2 =
    (left 0 * left 3 - left 1 * left 2) *
      (right 0 * right 3 - right 1 * right 2) := by
  simp [hgcdMatProductEntry]
  ring

/-- Algebraic composition of the two real recursive matrices around the
middle divrem step, in the exact orientation used by the generated code. -/
theorem hgcdRecursiveCombined_preserves_transform
    {p : Nat}
    (left right currentA currentB remainder quotient finalA finalB :
      Polynomial (ZMod p))
    (first second : Fin 4 → Polynomial (ZMod p))
    (hFirst : CLPoly.Impl.StrictHGCDRefinement.HgcdTransform left right
      currentA currentB (first 0) (first 1) (first 2) (first 3))
    (hDivision : currentA = quotient * currentB + remainder)
    (hSecond : CLPoly.Impl.StrictHGCDRefinement.HgcdTransform currentB
      remainder finalA finalB (second 0) (second 1) (second 2) (second 3)) :
    CLPoly.Impl.StrictHGCDRefinement.HgcdTransform left right finalA finalB
      (hgcdMatProductEntry first
        (hgcdMatApplyQuotientEntries second quotient) 0)
      (hgcdMatProductEntry first
        (hgcdMatApplyQuotientEntries second quotient) 1)
      (hgcdMatProductEntry first
        (hgcdMatApplyQuotientEntries second quotient) 2)
      (hgcdMatProductEntry first
        (hgcdMatApplyQuotientEntries second quotient) 3) := by
  rcases hFirst with ⟨hLeft, hRight⟩
  rcases hSecond with ⟨hCurrentB, hRemainder⟩
  constructor
  · rw [hLeft, hDivision, hCurrentB, hRemainder]
    simp [hgcdMatProductEntry, hgcdMatApplyQuotientEntries,
      hgcdMatSwapEntries, hgcdMatQuotientUpdateEntries, Function.update]
    ring
  · rw [hRight, hDivision, hCurrentB, hRemainder]
    simp [hgcdMatProductEntry, hgcdMatApplyQuotientEntries,
      hgcdMatSwapEntries, hgcdMatQuotientUpdateEntries, Function.update]
    ring

/-- The source return sign `-(sgnR*sgnS)` equals the determinant sign of
the exact product `R * ([[q,1],[1,0]] * S)`. -/
theorem hgcdRecursiveCombined_preserves_signedDet
    {p : Nat} (sgnR sgnS : Int) (quotient : Polynomial (ZMod p))
    (first second : Fin 4 → Polynomial (ZMod p))
    (hFirst : CLPoly.Impl.StrictHGCDRefinement.HgcdSignedDet sgnR
      (first 0) (first 1) (first 2) (first 3))
    (hSecond : CLPoly.Impl.StrictHGCDRefinement.HgcdSignedDet sgnS
      (second 0) (second 1) (second 2) (second 3)) :
    CLPoly.Impl.StrictHGCDRefinement.HgcdSignedDet (-(sgnR * sgnS))
      (hgcdMatProductEntry first
        (hgcdMatApplyQuotientEntries second quotient) 0)
      (hgcdMatProductEntry first
        (hgcdMatApplyQuotientEntries second quotient) 1)
      (hgcdMatProductEntry first
        (hgcdMatApplyQuotientEntries second quotient) 2)
      (hgcdMatProductEntry first
        (hgcdMatApplyQuotientEntries second quotient) 3) := by
  have hProduct := hgcdMatProductEntry_det first
    (hgcdMatApplyQuotientEntries second quotient)
  have hStep := hgcdMatApplyQuotientEntries_det second quotient
  rcases hFirst with ⟨hsgnR, hdetR⟩ | ⟨hsgnR, hdetR⟩
  · rcases hSecond with ⟨hsgnS, hdetS⟩ | ⟨hsgnS, hdetS⟩
    · subst sgnR
      subst sgnS
      right
      constructor
      · norm_num
      · rw [hProduct, hStep, hdetR, hdetS]
        ring
    · subst sgnR
      subst sgnS
      left
      constructor
      · norm_num
      · rw [hProduct, hStep, hdetR, hdetS]
        ring
  · rcases hSecond with ⟨hsgnS, hdetS⟩ | ⟨hsgnS, hdetS⟩
    · subst sgnR
      subst sgnS
      left
      constructor
      · norm_num
      · rw [hProduct, hStep, hdetR, hdetS]
        ring
    · subst sgnR
      subst sgnS
      right
      constructor
      · norm_num
      · rw [hProduct, hStep, hdetR, hdetS]
        ring

/-- Physical provider for the two concrete column executions.  It is
quantified over the actual first result and contains no L2 polynomial. -/
def HgcdMatApplyQuotientWorkspaceProvider (this : DenseUPolyZp)
    (S : HgcdMat) (hS : S.Valid) (q : RawPtr UInt64) (lenQ : Nat)
    (T scratch : RawPtr UInt64) (heap : RawHeap) : Prop :=
  ∀ first,
    hgcdMatQuotientEntry this (hgcdMatSwapRows S hS)
      (hgcdMatSwapRows_valid S hS) (0 : Fin 4) (2 : Fin 4) q lenQ T
      scratch heap = .ok first →
    HgcdMatQuotientEntryWorkspace (hgcdMatSwapRows S hS)
        (hgcdMatSwapRows_valid S hS) (0 : Fin 4) (2 : Fin 4) q lenQ T
        scratch heap ∧
      HgcdMatQuotientEntryWorkspace first.matrix first.valid
        (1 : Fin 4) (3 : Fin 4) q lenQ T scratch first.heap

/-- End-to-end refinement of the exact source block
`S := [[q,1],[1,0]] * S`: descriptor swaps and both real guarded
multiplication/addition executions are all consumed in source order. -/
theorem hgcdMatApplyQuotient_refines (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (S : HgcdMat) (hS : S.Valid) (q : RawPtr UInt64) (lenQ : Nat)
    (T scratch : RawPtr UInt64) (heap : RawHeap)
    (result : HgcdMatQuotientResult)
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (quotient : Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (physical : HgcdMatApplyQuotientWorkspaceProvider this S hS q lenQ T
      scratch heap)
    (hQ : RawDensePolyRep this heap q lenQ quotient)
    (hMatrix : HgcdMatRawDenseRep this heap S entries hS)
    (hrun : hgcdMatApplyQuotient this S hS q lenQ T scratch heap =
      .ok result) :
    HgcdMatRawDenseRep this result.heap result.matrix
      (hgcdMatApplyQuotientEntries entries quotient) result.valid ∧
    RawDensePolyRep this result.heap q lenQ quotient := by
  rcases hgcdMatApplyQuotient_exec this S hS q lenQ T scratch heap result
      hrun with ⟨first, hfirst, hsecond⟩
  have hwork := physical first hfirst
  have hSwapped := hgcdMatSwapRows_refines this S hS entries heap hMatrix
  have hFirst := hgcdMatQuotientEntry_refines this (hgcdMatSwapRows S hS)
    (hgcdMatSwapRows_valid S hS) (0 : Fin 4) (2 : Fin 4) q lenQ T
    scratch heap first (hgcdMatSwapEntries entries) quotient hcfg hp hwork.1
    hQ hSwapped hfirst
  have hSecond := hgcdMatQuotientEntry_refines this first.matrix first.valid
    (1 : Fin 4) (3 : Fin 4) q lenQ T scratch first.heap
    (HgcdMatQuotientEntryResult.mk result.heap result.matrix result.valid)
    (hgcdMatQuotientUpdateEntries (hgcdMatSwapEntries entries) quotient
      (0 : Fin 4) (2 : Fin 4)) quotient hcfg hp hwork.2 hFirst.2 hFirst.1
    hsecond
  simpa [hgcdMatApplyQuotientEntries] using hSecond

/-- Exact descriptor bounds induced by the source row swap followed by its
two quotient column updates. -/
structure HgcdMatApplyQuotientLengthBounds
    (S result : HgcdMat) (hS : S.Valid) (hResult : result.Valid)
    (lenQ : Nat) : Prop where
  row0 : hgcdMatLen result hResult (0 : Fin 4) ≤
    max (hgcdMatLen S hS (2 : Fin 4))
      (lenQ + hgcdMatLen S hS (0 : Fin 4) - 1)
  row1 : hgcdMatLen result hResult (1 : Fin 4) ≤
    max (hgcdMatLen S hS (3 : Fin 4))
      (lenQ + hgcdMatLen S hS (1 : Fin 4) - 1)
  row2 : hgcdMatLen result hResult (2 : Fin 4) =
    hgcdMatLen S hS (0 : Fin 4)
  row3 : hgcdMatLen result hResult (3 : Fin 4) =
    hgcdMatLen S hS (1 : Fin 4)

/-- The raw semantic result of the actual quotient block determines all four
source descriptor bounds, including exact preservation of the lower row. -/
theorem hgcdMatApplyQuotientEntries_length_bounds (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (inputHeap resultHeap quotientHeap : RawHeap)
    (S result : HgcdMat) (hS : S.Valid) (hResult : result.Valid)
    (q : RawPtr UInt64) (lenQ : Nat)
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (quotient : Polynomial (ZMod this._p.toNat))
    (hMatrix : HgcdMatRawDenseRep this inputHeap S entries hS)
    (hQ : RawDensePolyRep this quotientHeap q lenQ quotient)
    (hOutput : HgcdMatRawDenseRep this resultHeap result
      (hgcdMatApplyQuotientEntries entries quotient) hResult) :
    HgcdMatApplyQuotientLengthBounds S result hS hResult lenQ := by
  constructor
  · simpa [hgcdMatApplyQuotientEntries, hgcdMatSwapEntries,
      hgcdMatQuotientUpdateEntries, Function.update] using
      rawDensePolyRep_add_mul_length_le this inputHeap quotientHeap inputHeap
        resultHeap (hgcdMatPtr S hS (2 : Fin 4)) q
        (hgcdMatPtr S hS (0 : Fin 4))
        (hgcdMatPtr result hResult (0 : Fin 4))
        (hgcdMatLen S hS (2 : Fin 4)) lenQ
        (hgcdMatLen S hS (0 : Fin 4))
        (hgcdMatLen result hResult (0 : Fin 4))
        (entries 2) quotient (entries 0) (hMatrix 2) hQ (hMatrix 0)
        (hOutput 0)
  · simpa [hgcdMatApplyQuotientEntries, hgcdMatSwapEntries,
      hgcdMatQuotientUpdateEntries, Function.update] using
      rawDensePolyRep_add_mul_length_le this inputHeap quotientHeap inputHeap
        resultHeap (hgcdMatPtr S hS (3 : Fin 4)) q
        (hgcdMatPtr S hS (1 : Fin 4))
        (hgcdMatPtr result hResult (1 : Fin 4))
        (hgcdMatLen S hS (3 : Fin 4)) lenQ
        (hgcdMatLen S hS (1 : Fin 4))
        (hgcdMatLen result hResult (1 : Fin 4))
        (entries 3) quotient (entries 1) (hMatrix 3) hQ (hMatrix 1)
        (hOutput 1)
  · exact (rawDensePolyRep_length_eq this inputHeap resultHeap
      (hgcdMatPtr S hS (0 : Fin 4))
      (hgcdMatPtr result hResult (2 : Fin 4))
      (hgcdMatLen S hS (0 : Fin 4))
      (hgcdMatLen result hResult (2 : Fin 4)) (entries 0) (hMatrix 0)
      (by simpa [hgcdMatApplyQuotientEntries, hgcdMatSwapEntries,
        hgcdMatQuotientUpdateEntries, Function.update] using hOutput 2)).symm
  · exact (rawDensePolyRep_length_eq this inputHeap resultHeap
      (hgcdMatPtr S hS (1 : Fin 4))
      (hgcdMatPtr result hResult (3 : Fin 4))
      (hgcdMatLen S hS (1 : Fin 4))
      (hgcdMatLen result hResult (3 : Fin 4)) (entries 1) (hMatrix 1)
      (by simpa [hgcdMatApplyQuotientEntries, hgcdMatSwapEntries,
        hgcdMatQuotientUpdateEntries, Function.update] using hOutput 3)).symm

/-- L2 low parts computed by the two sign-selected reconstruction blocks. -/
noncomputable def hgcdReconstructedLowB {p : Nat}
    (entries : Fin 4 → Polynomial (ZMod p))
    (lowA lowB : Polynomial (ZMod p)) (sgn : Int) :=
  if sgn < 0 then entries 2 * lowA - entries 0 * lowB
  else entries 0 * lowB - entries 2 * lowA

noncomputable def hgcdReconstructedLowA {p : Nat}
    (entries : Fin 4 → Polynomial (ZMod p))
    (lowA lowB : Polynomial (ZMod p)) (sgn : Int) :=
  if sgn < 0 then entries 1 * lowB - entries 3 * lowA
  else entries 3 * lowA - entries 1 * lowB

/-- The source's sign-selected low reconstruction is the inverse of the
actual HGCD matrix.  Adding the shifted high outputs therefore extends the
high-half transform to the complete input polynomials. -/
theorem hgcdReconstruction_preserves_transform
    {R : Type*} [CommRing R]
    (lowA lowB highLeft highRight highA highB shiftTerm : R)
    (sgn : Int) (entries : Fin 4 → R)
    (hHigh : CLPoly.Impl.StrictHGCDRefinement.HgcdTransform
      highLeft highRight highA highB
      (entries 0) (entries 1) (entries 2) (entries 3))
    (hDet : CLPoly.Impl.StrictHGCDRefinement.HgcdSignedDet sgn
      (entries 0) (entries 1) (entries 2) (entries 3)) :
    CLPoly.Impl.StrictHGCDRefinement.HgcdTransform
      (lowA + shiftTerm * highLeft) (lowB + shiftTerm * highRight)
      ((if sgn < 0 then entries 1 * lowB - entries 3 * lowA
        else entries 3 * lowA - entries 1 * lowB) + shiftTerm * highA)
      ((if sgn < 0 then entries 2 * lowA - entries 0 * lowB
        else entries 0 * lowB - entries 2 * lowA) + shiftTerm * highB)
      (entries 0) (entries 1) (entries 2) (entries 3) := by
  rcases hHigh with ⟨hHighLeft, hHighRight⟩
  rcases hDet with ⟨hsgn, hdet⟩ | ⟨hsgn, hdet⟩
  · subst sgn
    constructor
    · rw [hHighLeft]
      simp
      linear_combination -hdet * lowA
    · rw [hHighRight]
      simp
      linear_combination -hdet * lowB
  · subst sgn
    constructor
    · rw [hHighLeft]
      simp
      linear_combination hdet * lowA
    · rw [hHighRight]
      simp
      linear_combination hdet * lowB

/-- Polynomial specialization matching exactly the values returned by
`hgcdRecursiveReconstructPair_refines`. -/
theorem hgcdReconstructedPair_preserves_transform
    {p : Nat} (lowA lowB highLeft highRight highA highB : Polynomial (ZMod p))
    (shift : Nat) (sgn : Int)
    (entries : Fin 4 → Polynomial (ZMod p))
    (hHigh : CLPoly.Impl.StrictHGCDRefinement.HgcdTransform
      highLeft highRight highA highB
      (entries 0) (entries 1) (entries 2) (entries 3))
    (hDet : CLPoly.Impl.StrictHGCDRefinement.HgcdSignedDet sgn
      (entries 0) (entries 1) (entries 2) (entries 3)) :
    CLPoly.Impl.StrictHGCDRefinement.HgcdTransform
      (lowA + Polynomial.X ^ shift * highLeft)
      (lowB + Polynomial.X ^ shift * highRight)
      (hgcdReconstructedLowA entries lowA lowB sgn +
        Polynomial.X ^ shift * highA)
      (hgcdReconstructedLowB entries lowA lowB sgn +
        Polynomial.X ^ shift * highB)
      (entries 0) (entries 1) (entries 2) (entries 3) := by
  simpa [hgcdReconstructedLowA, hgcdReconstructedLowB] using
    hgcdReconstruction_preserves_transform lowA lowB highLeft highRight
      highA highB (Polynomial.X ^ shift) sgn entries hHigh hDet

/-- Physical obligations for the exact four-call reconstruction pair.  All
frame fields mention only heap cells; no expected L2 polynomial occurs in
this contract. -/
structure HgcdRecursiveReconstructPairWorkspace (this : DenseUPolyZp)
    (A B T0 lowA lowB highA highB scratch : RawPtr UInt64)
    (lenLowA lenLowB lenHighA lenHighB shift : Nat)
    (M : HgcdMat) (hM : M.Valid) (sgn : Int) (heap : RawHeap)
    (heap1 : RawHeap) (lowLenB : Nat) (liftedB : HgcdLiftHighResult)
    (heap3 : RawHeap) (lowLenA : Nat) (liftedA : HgcdLiftHighResult) : Prop where
  reconstructB : HgcdReconstructWorkspace heap B T0
    (hgcdMatPtr M hM (2 : Fin 4)) (hgcdMatLen M hM (2 : Fin 4)) lowA lenLowA
    (hgcdMatPtr M hM (0 : Fin 4)) (hgcdMatLen M hM (0 : Fin 4)) lowB lenLowB
    scratch
  highBLayout : RawHeap.SameLayout heap heap1
  highBPrefix : SameU64Prefix heap heap1 highB lenHighB
  liftB : HgcdLiftHighWorkspace heap1 B highB lowLenB shift lenHighB
  aInputLayout : RawHeap.SameLayout heap liftedB.heap
  aMatrixPrefix : ∀ i : Fin 4, SameU64Prefix heap liftedB.heap
    (hgcdMatPtr M hM i) (hgcdMatLen M hM i)
  aLowAPrefix : SameU64Prefix heap liftedB.heap lowA lenLowA
  aLowBPrefix : SameU64Prefix heap liftedB.heap lowB lenLowB
  reconstructA : HgcdReconstructWorkspace liftedB.heap A T0
    (hgcdMatPtr M hM (3 : Fin 4)) (hgcdMatLen M hM (3 : Fin 4)) lowA lenLowA
    (hgcdMatPtr M hM (1 : Fin 4)) (hgcdMatLen M hM (1 : Fin 4)) lowB lenLowB
    scratch
  highALayout : RawHeap.SameLayout heap heap3
  highAPrefix : SameU64Prefix heap heap3 highA lenHighA
  liftA : HgcdLiftHighWorkspace heap3 A highA lowLenA shift lenHighA
  finalBLayout : RawHeap.SameLayout liftedB.heap liftedA.heap
  finalBPrefix : SameU64Prefix liftedB.heap liftedA.heap B liftedB.length
  finalMatrixLayout : RawHeap.SameLayout heap liftedA.heap
  finalMatrixPrefix : ∀ i : Fin 4, SameU64Prefix heap liftedA.heap
    (hgcdMatPtr M hM i) (hgcdMatLen M hM i)

def HgcdRecursiveReconstructPairWorkspaceProvider (this : DenseUPolyZp)
    (A B T0 lowA lowB highA highB scratch : RawPtr UInt64)
    (lenLowA lenLowB lenHighA lenHighB shift : Nat)
    (M : HgcdMat) (hM : M.Valid) (sgn : Int) (heap : RawHeap) : Prop :=
  ∀ heap1 lowLenB liftedB heap3 lowLenA liftedA,
    hgcdRecursiveReconstructB this B T0 (hgcdMatPtr M hM (2 : Fin 4))
        (hgcdMatPtr M hM (0 : Fin 4)) lowA lowB scratch
        (hgcdMatLen M hM (2 : Fin 4)) (hgcdMatLen M hM (0 : Fin 4))
        lenLowA lenLowB sgn heap = .ok (heap1, lowLenB) →
    hgcdRecursiveLiftHigh this B highB lowLenB shift lenHighB heap1 =
        .ok liftedB →
    hgcdRecursiveReconstructA this A T0 (hgcdMatPtr M hM (3 : Fin 4))
        (hgcdMatPtr M hM (1 : Fin 4)) lowA lowB scratch
        (hgcdMatLen M hM (3 : Fin 4)) (hgcdMatLen M hM (1 : Fin 4))
        lenLowA lenLowB sgn liftedB.heap = .ok (heap3, lowLenA) →
    hgcdRecursiveLiftHigh this A highA lowLenA shift lenHighA heap3 =
        .ok liftedA →
    HgcdRecursiveReconstructPairWorkspace this A B T0 lowA lowB highA highB
      scratch lenLowA lenLowB lenHighA lenHighB shift M hM sgn heap heap1
      lowLenB liftedB heap3 lowLenA liftedA

/-- The four physical reconstruction calls do not overwrite the matrix
returned by the first recursive child. -/
theorem hgcdRecursiveReconstructPair_preserves_matrix (this : DenseUPolyZp)
    (A B T0 lowA lowB highA highB scratch : RawPtr UInt64)
    (lenLowA lenLowB lenHighA lenHighB shift : Nat)
    (M : HgcdMat) (hM : M.Valid) (sgn : Int) (heap : RawHeap)
    (result : HgcdRecursiveReconstructPairResult)
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (physical : HgcdRecursiveReconstructPairWorkspaceProvider this A B T0
      lowA lowB highA highB scratch lenLowA lenLowB lenHighA lenHighB shift
      M hM sgn heap)
    (hMatrix : HgcdMatRawDenseRep this heap M entries hM)
    (hrun : hgcdRecursiveReconstructPair this A B T0 lowA lowB highA highB
      scratch lenLowA lenLowB lenHighA lenHighB shift M hM sgn heap =
        .ok result) :
    HgcdMatRawDenseRep this result.heap M entries hM := by
  rcases hgcdRecursiveReconstructPair_exec this A B T0 lowA lowB highA
      highB scratch lenLowA lenLowB lenHighA lenHighB shift M hM sgn heap
      result hrun with
    ⟨heap1, lowLenB, liftedB, heap3, lowLenA, liftedA,
      hBRun, hLiftBRun, hARun, hLiftARun, hHeap, _, _⟩
  have hwork := physical heap1 lowLenB liftedB heap3 lowLenA liftedA hBRun
    hLiftBRun hARun hLiftARun
  intro i
  rw [hHeap]
  exact rawDensePolyRep_of_same_prefix this heap liftedA.heap
    (hgcdMatPtr M hM i) (hgcdMatLen M hM i) (entries i)
    hwork.finalMatrixLayout (hwork.finalMatrixPrefix i) (hMatrix i)

/-- Complete raw semantic refinement of the actual paired reconstruction.
The returned polynomials arise from the four generated calls in source
order, including both final normalizations. -/
theorem hgcdRecursiveReconstructPair_refines (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (A B T0 lowA lowB highA highB scratch : RawPtr UInt64)
    (lenLowA lenLowB lenHighA lenHighB shift : Nat)
    (M : HgcdMat) (hM : M.Valid) (sgn : Int) (heap : RawHeap)
    (result : HgcdRecursiveReconstructPairResult)
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (polyLowA polyLowB polyHighA polyHighB :
      Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (physical : HgcdRecursiveReconstructPairWorkspaceProvider this A B T0
      lowA lowB highA highB scratch lenLowA lenLowB lenHighA lenHighB shift
      M hM sgn heap)
    (hMatrix : HgcdMatRawDenseRep this heap M entries hM)
    (hLowA : RawCanonicalPolySlice this heap lowA lenLowA polyLowA)
    (hLowB : RawCanonicalPolySlice this heap lowB lenLowB polyLowB)
    (hHighA : RawDensePolyRep this heap highA lenHighA polyHighA)
    (hHighB : RawDensePolyRep this heap highB lenHighB polyHighB)
    (hrun : hgcdRecursiveReconstructPair this A B T0 lowA lowB highA highB
      scratch lenLowA lenLowB lenHighA lenHighB shift M hM sgn heap =
      .ok result) :
    RawDensePolyRep this result.heap A result.lenA
        (hgcdReconstructedLowA entries polyLowA polyLowB sgn +
          Polynomial.X ^ shift * polyHighA) ∧
      RawDensePolyRep this result.heap B result.lenB
        (hgcdReconstructedLowB entries polyLowA polyLowB sgn +
          Polynomial.X ^ shift * polyHighB) ∧
      result.lenA ≤ max (shift + lenHighA)
        (max
          (hgcdMatLen M hM (3 : Fin 4) + lenLowA - 1)
          (hgcdMatLen M hM (1 : Fin 4) + lenLowB - 1)) ∧
      result.lenB ≤ max (shift + lenHighB)
        (max
          (hgcdMatLen M hM (2 : Fin 4) + lenLowA - 1)
          (hgcdMatLen M hM (0 : Fin 4) + lenLowB - 1)) ∧
      (max
          (hgcdMatLen M hM (3 : Fin 4) + lenLowA - 1)
          (hgcdMatLen M hM (1 : Fin 4) + lenLowB - 1) <
            shift + lenHighA →
        0 < lenHighA → result.lenA = shift + lenHighA) ∧
      (max
          (hgcdMatLen M hM (2 : Fin 4) + lenLowA - 1)
          (hgcdMatLen M hM (0 : Fin 4) + lenLowB - 1) <
            shift + lenHighB →
        0 < lenHighB → result.lenB = shift + lenHighB) := by
  rcases hgcdRecursiveReconstructPair_exec this A B T0 lowA lowB highA
      highB scratch lenLowA lenLowB lenHighA lenHighB shift M hM sgn heap
      result hrun with
    ⟨heap1, lowLenB, liftedB, heap3, lowLenA, liftedA,
      hBRun, hLiftBRun, hARun, hLiftARun, hHeap, hLenA, hLenB⟩
  have hwork := physical heap1 lowLenB liftedB heap3 lowLenA liftedA hBRun
    hLiftBRun hARun hLiftARun
  rcases hgcdRecursiveReconstructB_refines this B T0
      (hgcdMatPtr M hM (2 : Fin 4)) (hgcdMatPtr M hM (0 : Fin 4)) lowA
      lowB scratch (hgcdMatLen M hM (2 : Fin 4))
      (hgcdMatLen M hM (0 : Fin 4)) lenLowA lenLowB sgn heap
      (entries 2) (entries 0) polyLowA polyLowB hcfg hp
      hwork.reconstructB (hMatrix 2) (hMatrix 0) hLowA hLowB with
    ⟨heapB, lenB0, hBRun', _, hB0, hLowLenB⟩
  have hEqB : (heapB, lenB0) = (heap1, lowLenB) :=
    Except.ok.inj (hBRun'.symm.trans hBRun)
  cases hEqB
  have hHighB1 := rawDensePolyRep_of_same_prefix this heap heap1 highB
    lenHighB polyHighB hwork.highBLayout hwork.highBPrefix hHighB
  have hpWord : this._p ≠ 0 := by
    intro hzero
    have := congrArg UInt64.toNat hzero
    simp at this
    omega
  rcases hgcdRecursiveLiftHigh_terminates this B highB lowLenB shift lenHighB
      heap1 (hgcdReconstructedLowB entries polyLowA polyLowB sgn) hpWord
      hwork.liftB (by simpa [hgcdReconstructedLowB] using hB0) with
    ⟨boundedB, hBoundedBRun, _, hBoundB⟩
  have hEqBoundedB : boundedB = liftedB :=
    Except.ok.inj (hBoundedBRun.symm.trans hLiftBRun)
  subst boundedB
  rcases hgcdRecursiveLiftHigh_refines this B highB lowLenB shift lenHighB
      heap1 (hgcdReconstructedLowB entries polyLowA polyLowB sgn) polyHighB
      hpWord hwork.liftB (by simpa [hgcdReconstructedLowB] using hB0)
      hHighB1 with ⟨liftedB', hLiftBRun', _, hBFinal⟩
  have hEqLiftB : liftedB' = liftedB :=
    Except.ok.inj (hLiftBRun'.symm.trans hLiftBRun)
  subst liftedB'
  have matrixAtB (i : Fin 4) : RawDensePolyRep this liftedB.heap
      (hgcdMatPtr M hM i) (hgcdMatLen M hM i) (entries i) :=
    rawDensePolyRep_of_same_prefix this heap liftedB.heap
      (hgcdMatPtr M hM i) (hgcdMatLen M hM i) (entries i)
      hwork.aInputLayout (hwork.aMatrixPrefix i) (hMatrix i)
  have lowAAtB := rawCanonicalPolySlice_of_same_prefix this heap liftedB.heap
    lowA lenLowA polyLowA hwork.aInputLayout hwork.aLowAPrefix hLowA
  have lowBAtB := rawCanonicalPolySlice_of_same_prefix this heap liftedB.heap
    lowB lenLowB polyLowB hwork.aInputLayout hwork.aLowBPrefix hLowB
  rcases hgcdRecursiveReconstructA_refines this A T0
      (hgcdMatPtr M hM (3 : Fin 4)) (hgcdMatPtr M hM (1 : Fin 4)) lowA
      lowB scratch (hgcdMatLen M hM (3 : Fin 4))
      (hgcdMatLen M hM (1 : Fin 4)) lenLowA lenLowB sgn liftedB.heap
      (entries 3) (entries 1) polyLowA polyLowB hcfg hp hwork.reconstructA
      (matrixAtB 3) (matrixAtB 1) lowAAtB lowBAtB with
    ⟨heapA, lenA0, hARun', _, hA0, hLowLenA⟩
  have hEqA : (heapA, lenA0) = (heap3, lowLenA) :=
    Except.ok.inj (hARun'.symm.trans hARun)
  cases hEqA
  have hHighA3 := rawDensePolyRep_of_same_prefix this heap heap3 highA
    lenHighA polyHighA hwork.highALayout hwork.highAPrefix hHighA
  rcases hgcdRecursiveLiftHigh_terminates this A highA lowLenA shift lenHighA
      heap3 (hgcdReconstructedLowA entries polyLowA polyLowB sgn) hpWord
      hwork.liftA (by simpa [hgcdReconstructedLowA] using hA0) with
    ⟨boundedA, hBoundedARun, _, hBoundA⟩
  have hEqBoundedA : boundedA = liftedA :=
    Except.ok.inj (hBoundedARun.symm.trans hLiftARun)
  subst boundedA
  rcases hgcdRecursiveLiftHigh_refines this A highA lowLenA shift lenHighA
      heap3 (hgcdReconstructedLowA entries polyLowA polyLowB sgn) polyHighA
      hpWord hwork.liftA (by simpa [hgcdReconstructedLowA] using hA0)
      hHighA3 with ⟨liftedA', hLiftARun', _, hAFinal⟩
  have hEqLiftA : liftedA' = liftedA :=
    Except.ok.inj (hLiftARun'.symm.trans hLiftARun)
  subst liftedA'
  have hBAtFinal := rawDensePolyRep_of_same_prefix this liftedB.heap
    liftedA.heap B liftedB.length
    (hgcdReconstructedLowB entries polyLowA polyLowB sgn +
      Polynomial.X ^ shift * polyHighB) hwork.finalBLayout
      hwork.finalBPrefix hBFinal
  have hExactA : max
        (hgcdMatLen M hM (3 : Fin 4) + lenLowA - 1)
        (hgcdMatLen M hM (1 : Fin 4) + lenLowB - 1) <
          shift + lenHighA →
      0 < lenHighA → liftedA.length = shift + lenHighA := by
    intro hlow hhigh
    exact rawDensePolyRep_add_shift_length_eq_of_lt this heap3 heap3 liftedA.heap
      A highA A lowLenA lenHighA liftedA.length shift
      (hgcdReconstructedLowA entries polyLowA polyLowB sgn) polyHighA hA0
      hHighA3 hAFinal (hLowLenA.trans_lt hlow) hhigh
  have hExactB : max
        (hgcdMatLen M hM (2 : Fin 4) + lenLowA - 1)
        (hgcdMatLen M hM (0 : Fin 4) + lenLowB - 1) <
          shift + lenHighB →
      0 < lenHighB → liftedB.length = shift + lenHighB := by
    intro hlow hhigh
    exact rawDensePolyRep_add_shift_length_eq_of_lt this heap1 heap1 liftedB.heap
      B highB B lowLenB lenHighB liftedB.length shift
      (hgcdReconstructedLowB entries polyLowA polyLowB sgn) polyHighB hB0
      hHighB1 hBFinal (hLowLenB.trans_lt hlow) hhigh
  rw [hHeap, hLenA, hLenB]
  refine ⟨hAFinal, hBAtFinal, ?_, ?_, hExactA, hExactB⟩
  · exact hBoundA.trans (max_le_max (Nat.le_refl _) hLowLenA)
  · exact hBoundB.trans (max_le_max (Nat.le_refl _) hLowLenB)

/-- Semantic composition of a real first recursive result with the four-call
paired reconstruction.  The full-input transform and GCD theorem are derived
from the returned matrix, its signed determinant, and the actual low/high
decompositions of the two source operands. -/
theorem hgcdRecursiveReconstructPair_preserves_input (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (A B T0 lowA lowB highA highB scratch : RawPtr UInt64)
    (lenLowA lenLowB shift inputLength : Nat)
    (first : HgcdRecursiveResult)
    (result : HgcdRecursiveReconstructPairResult)
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (fullA fullB polyLowA polyLowB polyInputHighA polyInputHighB
      polyOutputHighA polyOutputHighB : Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (physical : HgcdRecursiveReconstructPairWorkspaceProvider this A B T0
      lowA lowB highA highB scratch lenLowA lenLowB first.lenA first.lenB
      shift first.matrix first.valid first.sgn first.heap)
    (hFirst : HgcdRecursiveRawInvariant this polyInputHighA polyInputHighB
      polyOutputHighA polyOutputHighB entries true highA highB inputLength
      first)
    (hLowA : RawCanonicalPolySlice this first.heap lowA lenLowA polyLowA)
    (hLowB : RawCanonicalPolySlice this first.heap lowB lenLowB polyLowB)
    (hFullA : fullA = polyLowA + Polynomial.X ^ shift * polyInputHighA)
    (hFullB : fullB = polyLowB + Polynomial.X ^ shift * polyInputHighB)
    (hrun : hgcdRecursiveReconstructPair this A B T0 lowA lowB highA highB
      scratch lenLowA lenLowB first.lenA first.lenB shift first.matrix
      first.valid first.sgn first.heap = .ok result) :
    RawDensePolyRep this result.heap A result.lenA
        (hgcdReconstructedLowA entries polyLowA polyLowB first.sgn +
          Polynomial.X ^ shift * polyOutputHighA) ∧
      RawDensePolyRep this result.heap B result.lenB
        (hgcdReconstructedLowB entries polyLowA polyLowB first.sgn +
          Polynomial.X ^ shift * polyOutputHighB) ∧
      CLPoly.Impl.StrictHGCDRefinement.HgcdTransform fullA fullB
        (hgcdReconstructedLowA entries polyLowA polyLowB first.sgn +
          Polynomial.X ^ shift * polyOutputHighA)
        (hgcdReconstructedLowB entries polyLowA polyLowB first.sgn +
          Polynomial.X ^ shift * polyOutputHighB)
        (entries 0) (entries 1) (entries 2) (entries 3) ∧
      CLPoly.Impl.StrictHGCDRefinement.HgcdSignedDet first.sgn
        (entries 0) (entries 1) (entries 2) (entries 3) ∧
      normalize (EuclideanDomain.gcd fullA fullB) =
        normalize (EuclideanDomain.gcd
          (hgcdReconstructedLowA entries polyLowA polyLowB first.sgn +
            Polynomial.X ^ shift * polyOutputHighA)
          (hgcdReconstructedLowB entries polyLowA polyLowB first.sgn +
            Polynomial.X ^ shift * polyOutputHighB)) ∧
      result.lenB ≤ max (shift + first.lenB)
        (max
          (hgcdMatLen first.matrix first.valid (2 : Fin 4) + lenLowA - 1)
          (hgcdMatLen first.matrix first.valid (0 : Fin 4) + lenLowB - 1)) := by
  have hMatrixSemantics := hFirst.matrixSemantics rfl
  rcases hgcdRecursiveReconstructPair_refines this A B T0 lowA lowB highA
      highB scratch lenLowA lenLowB first.lenA first.lenB shift first.matrix
      first.valid first.sgn first.heap result entries polyLowA polyLowB
      polyOutputHighA polyOutputHighB hcfg hp physical hMatrixSemantics.1
      hLowA hLowB hFirst.aRep hFirst.bRep hrun with
    ⟨hAResult, hBResult, _, hLength, _, _⟩
  let finalA := hgcdReconstructedLowA entries polyLowA polyLowB first.sgn +
    Polynomial.X ^ shift * polyOutputHighA
  let finalB := hgcdReconstructedLowB entries polyLowA polyLowB first.sgn +
    Polynomial.X ^ shift * polyOutputHighB
  have hTransform : CLPoly.Impl.StrictHGCDRefinement.HgcdTransform fullA fullB
      finalA finalB (entries 0) (entries 1) (entries 2) (entries 3) := by
    rw [hFullA, hFullB]
    exact hgcdReconstructedPair_preserves_transform polyLowA polyLowB
      polyInputHighA polyInputHighB polyOutputHighA polyOutputHighB shift
      first.sgn entries hMatrixSemantics.2.1 hMatrixSemantics.2.2
  have hGcd :=
    CLPoly.Impl.StrictHGCDRefinement.normalize_gcd_eq_of_hgcd_signed_transform
      first.sgn fullA fullB finalA finalB (entries 0) (entries 1)
      (entries 2) (entries 3) hTransform hMatrixSemantics.2.2
  simpa [finalA, finalB] using
    And.intro hAResult (And.intro hBResult
      (And.intro hTransform (And.intro hMatrixSemantics.2.2
        (And.intro hGcd hLength))))

/-- The physical paired-reconstruction bound closes against an enclosing
input length once its shifted high part and both real low products fit that
input.  This is the arithmetic form consumed by the well-founded recursive
call, not a runtime decrease test. -/
theorem hgcdRecursiveReconstructPair_lenB_le_input
    (resultLen inputLength shift lenHighB lenR2 lenLowA lenR0 lenLowB : Nat)
    (hresult : resultLen ≤ max (shift + lenHighB)
      (max (lenR2 + lenLowA - 1) (lenR0 + lenLowB - 1)))
    (hhigh : shift + lenHighB ≤ inputLength)
    (hsecond : lenR2 + lenLowA - 1 ≤ inputLength)
    (hzero : lenR0 + lenLowB - 1 ≤ inputLength) :
    resultLen ≤ inputLength := by
  exact hresult.trans (max_le hhigh (max_le hsecond hzero))

/-- Close the first reconstruction bound from the returned first-HGCD
matrix invariant.  Here `highLength = inputLength - inputLength / 2` and
the low slices have exactly the lengths selected by the source. -/
theorem hgcdRecursiveFirstReconstruct_lenB_lt_input
    (resultLen inputLength inputLengthB returnedLenA returnedLenB lenR2 lenR0 : Nat)
    (hresult : resultLen ≤
      max (inputLength / 2 + returnedLenB)
        (max
          (lenR2 + Nat.min inputLength (inputLength / 2) - 1)
          (lenR0 + Nat.min inputLengthB (inputLength / 2) - 1)))
    (hreturnedStop : returnedLenB <
      (inputLength - inputLength / 2) / 2 + 1)
    (hreturnedPos : 0 < returnedLenA)
    (hinputOrder : inputLengthB < inputLength)
    (hrow2 : lenR2 + returnedLenA ≤
      (inputLength - inputLength / 2) + 1)
    (hrow0 : lenR0 + returnedLenA ≤
      (inputLength - inputLength / 2) + 1) :
    resultLen < inputLength := by
  have hsplit : inputLength / 2 +
      (inputLength - inputLength / 2) = inputLength := by omega
  have hhigh : inputLength / 2 + returnedLenB < inputLength := by omega
  have hsecond : lenR2 + Nat.min inputLength (inputLength / 2) - 1 <
      inputLength := by
    have hmin : Nat.min inputLength (inputLength / 2) ≤
        inputLength / 2 := Nat.min_le_right _ _
    omega
  have hzero : lenR0 + Nat.min inputLengthB (inputLength / 2) - 1 <
      inputLength := by
    have hmin : Nat.min inputLengthB (inputLength / 2) ≤
        inputLength / 2 := Nat.min_le_right _ _
    omega
  exact hresult.trans_lt (max_lt hhigh (max_lt hsecond hzero))

/-- Concrete discharge of the proof-only first-reconstruction bound used by
the strictly decreasing recursive body.  Every premise is either a raw
representation, a physical workspace fact, or the length invariant returned
by the real first HGCD call. -/
theorem hgcdRecursiveFirstReconstruct_bound_of_invariant
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (A B T0 a b highA highB scratch : RawPtr UInt64)
    (lenA lenB : Nat) (first : HgcdRecursiveResult)
    (result : HgcdRecursiveReconstructPairResult)
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (polyLowA polyLowB polyHighA polyHighB :
      Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (hinputOrder : lenB < lenA)
    (hinvariant : HgcdRecursiveLengthInvariant (lenA - lenA / 2) first)
    (physical : HgcdRecursiveReconstructPairWorkspaceProvider this A B T0
      a b highA highB scratch (Nat.min lenA (lenA / 2))
      (Nat.min lenB (lenA / 2)) first.lenA first.lenB (lenA / 2)
      first.matrix first.valid first.sgn first.heap)
    (hMatrix : HgcdMatRawDenseRep this first.heap first.matrix entries
      first.valid)
    (hLowA : RawDensePolyRep this first.heap a
      (Nat.min lenA (lenA / 2)) polyLowA)
    (hLowB : RawDensePolyRep this first.heap b
      (Nat.min lenB (lenA / 2)) polyLowB)
    (hHighA : RawDensePolyRep this first.heap highA first.lenA polyHighA)
    (hHighB : RawDensePolyRep this first.heap highB first.lenB polyHighB)
    (hrun : hgcdRecursiveReconstructPair this A B T0 a b highA highB scratch
      (Nat.min lenA (lenA / 2)) (Nat.min lenB (lenA / 2)) first.lenA
      first.lenB (lenA / 2) first.matrix first.valid first.sgn first.heap =
        .ok result) :
    result.lenB < lenA := by
  have hrefines := hgcdRecursiveReconstructPair_refines this A B T0 a b
    highA highB scratch (Nat.min lenA (lenA / 2))
    (Nat.min lenB (lenA / 2)) first.lenA first.lenB (lenA / 2)
    first.matrix first.valid first.sgn first.heap result entries polyLowA
    polyLowB polyHighA polyHighB hcfg hp physical hMatrix
    hLowA.toCanonicalSlice hLowB.toCanonicalSlice
    hHighA hHighB hrun
  exact hgcdRecursiveFirstReconstruct_lenB_lt_input result.lenB lenA lenB
    first.lenA first.lenB
    (hgcdMatLen first.matrix first.valid (2 : Fin 4))
    (hgcdMatLen first.matrix first.valid (0 : Fin 4)) hrefines.2.2.2.1
    hinvariant.stopped hinvariant.positive hinputOrder
    (by simpa [hgcdMatLen, hgcdMatLenRaw] using hinvariant.row2A)
    (by simpa [hgcdMatLen, hgcdMatLenRaw] using hinvariant.row0A)

/-- The actual first reconstruction retains Euclidean operand order.  Its A
length is the shifted first-child A length, while every concrete B branch is
bounded by the child stop and row-A descriptor invariants. -/
theorem hgcdRecursiveFirstReconstruct_order
    (inputLength inputLengthB highLength shift firstLenA firstLenB
      r0 r2 resultLenA resultLenB : Nat)
    (hhigh : highLength = inputLength - shift)
    (hfirstAbove : highLength / 2 < firstLenA)
    (hfirstStop : firstLenB < highLength / 2 + 1)
    (hrow0A : r0 + firstLenA ≤ highLength + 1)
    (hrow2A : r2 + firstLenA ≤ highLength + 1)
    (hresultA : resultLenA = shift + firstLenA)
    (hresultB : resultLenB ≤ max (shift + firstLenB)
      (max
        (r2 + Nat.min inputLength shift - 1)
        (r0 + Nat.min inputLengthB shift - 1))) :
    resultLenB ≤ resultLenA := by
  have hminA : Nat.min inputLength shift ≤ shift := Nat.min_le_right _ _
  have hminB : Nat.min inputLengthB shift ≤ shift := Nat.min_le_right _ _
  apply hresultB.trans
  rw [hresultA]
  apply max_le
  · omega
  · apply max_le <;> omega

/-- Arithmetic closure for the second paired reconstruction.  The source
chooses `k` so that `k + lenC0` is exactly the reconstructed divisor length;
the sharp product `- 1` bounds then fit both final operands inside the outer
input length. -/
theorem hgcdRecursiveFinalReconstruct_lengths_le_input
    (outerLength reconstructedLenB k lenC0 lenD : Nat)
    (secondLenA secondLenB s0 s1 s2 s3 resultLenA resultLenB : Nat)
    (hsplit : k + lenC0 = reconstructedLenB)
    (hreconstructed : reconstructedLenB ≤ outerLength)
    (hsecondBound : secondLenA ≤ lenC0)
    (hsecondStop : secondLenB < lenC0 / 2 + 1)
    (hrow0 : s0 + secondLenA ≤ lenC0 + 1)
    (hrow1 : s1 + secondLenB ≤ lenC0 + 1)
    (hrow2 : s2 + secondLenA ≤ lenC0 + 1)
    (hrow3 : s3 + secondLenB ≤ lenC0 + 1)
    (hresultA : resultLenA ≤ max (k + secondLenA)
      (max
        (s3 + Nat.min reconstructedLenB k - 1)
        (s1 + Nat.min lenD k - 1)))
    (hresultB : resultLenB ≤ max (k + secondLenB)
      (max
        (s2 + Nat.min reconstructedLenB k - 1)
        (s0 + Nat.min lenD k - 1))) :
    resultLenA ≤ outerLength ∧ resultLenB ≤ outerLength := by
  have hminReconstructed : Nat.min reconstructedLenB k ≤ k :=
    Nat.min_le_right _ _
  have hminD : Nat.min lenD k ≤ k := Nat.min_le_right _ _
  constructor
  · apply hresultA.trans
    apply max_le
    · omega
    · apply max_le <;> omega
  · apply hresultB.trans
    apply max_le
    · omega
    · apply max_le <;> omega

/-- The second HGCD result's half-length coefficient bound makes its shifted
leading A term strictly dominate both real low reconstruction products.  The
normalized final A descriptor is therefore exactly `shift + secondLenA`. -/
theorem hgcdRecursiveFinalReconstruct_lenA_eq
    (inputLength shift sourceLenA sourceLenB secondLenA s1 s3 resultLenA : Nat)
    (hsecondAbove : inputLength / 2 < secondLenA)
    (hcoeff1 : s1 ≤ inputLength - inputLength / 2)
    (hcoeff3 : s3 ≤ inputLength - inputLength / 2)
    (hsourceA : sourceLenA ≤ shift) (hsourceB : sourceLenB ≤ shift)
    (hexact : max
        (s3 + sourceLenA - 1)
        (s1 + sourceLenB - 1) < shift + secondLenA →
      0 < secondLenA → resultLenA = shift + secondLenA) :
    resultLenA = shift + secondLenA ∧ 0 < resultLenA := by
  have hlow3 : s3 + sourceLenA - 1 <
      shift + secondLenA := by omega
  have hlow1 : s1 + sourceLenB - 1 <
      shift + secondLenA := by omega
  have heq := hexact (max_lt hlow3 hlow1) (by omega)
  exact ⟨heq, by omega⟩

/-- Concrete final-A leading proof for the real four-call paired
reconstruction.  All polynomial and descriptor facts come from the returned
second HGCD invariant and the successful raw execution. -/
theorem hgcdRecursiveFinalReconstruct_lenA_eq_of_invariant
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (A B T0 lowA lowB highA highB scratch : RawPtr UInt64)
    (lenLowA lenLowB shift inputLength : Nat)
    (second : HgcdRecursiveResult)
    (result : HgcdRecursiveReconstructPairResult)
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (polyLowA polyLowB polyHighA polyHighB :
      Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (hinvariant : HgcdRecursiveLengthInvariant inputLength second)
    (hlenLowA : lenLowA ≤ shift) (hlenLowB : lenLowB ≤ shift)
    (physical : HgcdRecursiveReconstructPairWorkspaceProvider this A B T0
      lowA lowB highA highB scratch lenLowA lenLowB second.lenA second.lenB
      shift second.matrix second.valid second.sgn second.heap)
    (hMatrix : HgcdMatRawDenseRep this second.heap second.matrix entries
      second.valid)
    (hLowA : RawCanonicalPolySlice this second.heap lowA lenLowA polyLowA)
    (hLowB : RawCanonicalPolySlice this second.heap lowB lenLowB polyLowB)
    (hHighA : RawDensePolyRep this second.heap highA second.lenA polyHighA)
    (hHighB : RawDensePolyRep this second.heap highB second.lenB polyHighB)
    (hrun : hgcdRecursiveReconstructPair this A B T0 lowA lowB highA highB
      scratch lenLowA lenLowB second.lenA second.lenB shift second.matrix
      second.valid second.sgn second.heap = .ok result) :
    result.lenA = shift + second.lenA ∧ 0 < result.lenA := by
  have hrefines := hgcdRecursiveReconstructPair_refines this A B T0 lowA
    lowB highA highB scratch lenLowA lenLowB second.lenA second.lenB shift
    second.matrix second.valid second.sgn second.heap result entries polyLowA
    polyLowB polyHighA polyHighB hcfg hp physical hMatrix
    hLowA hLowB hHighA hHighB hrun
  exact hgcdRecursiveFinalReconstruct_lenA_eq inputLength shift lenLowA
    lenLowB second.lenA
    (hgcdMatLen second.matrix second.valid (1 : Fin 4))
    (hgcdMatLen second.matrix second.valid (3 : Fin 4)) result.lenA
    hinvariant.aboveHalf
    (by simpa [hgcdMatLen, hgcdMatLenRaw] using
      hinvariant.coeffBound (1 : Fin 4))
    (by simpa [hgcdMatLen, hgcdMatLenRaw] using
      hinvariant.coeffBound (3 : Fin 4))
    hlenLowA hlenLowB
    hrefines.2.2.2.2.1

/-- Execution-level order bridge for the first reconstructed pair.  Both the
exact leading-A length and the B upper bound are extracted from the same
successful four-call raw reconstruction, so the second recursive call gets
its operand order without a specification-side assumption. -/
theorem hgcdRecursiveFirstReconstruct_order_of_invariant
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (A B T0 a b highA highB scratch : RawPtr UInt64)
    (lenA lenB : Nat) (first : HgcdRecursiveResult)
    (result : HgcdRecursiveReconstructPairResult)
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (polyLowA polyLowB polyHighA polyHighB :
      Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (hinvariant : HgcdRecursiveLengthInvariant (lenA - lenA / 2) first)
    (physical : HgcdRecursiveReconstructPairWorkspaceProvider this A B T0
      a b highA highB scratch (Nat.min lenA (lenA / 2))
      (Nat.min lenB (lenA / 2)) first.lenA first.lenB (lenA / 2)
      first.matrix first.valid first.sgn first.heap)
    (hMatrix : HgcdMatRawDenseRep this first.heap first.matrix entries
      first.valid)
    (hLowA : RawDensePolyRep this first.heap a
      (Nat.min lenA (lenA / 2)) polyLowA)
    (hLowB : RawDensePolyRep this first.heap b
      (Nat.min lenB (lenA / 2)) polyLowB)
    (hHighA : RawDensePolyRep this first.heap highA first.lenA polyHighA)
    (hHighB : RawDensePolyRep this first.heap highB first.lenB polyHighB)
    (hrun : hgcdRecursiveReconstructPair this A B T0 a b highA highB scratch
      (Nat.min lenA (lenA / 2)) (Nat.min lenB (lenA / 2)) first.lenA
      first.lenB (lenA / 2) first.matrix first.valid first.sgn first.heap =
        .ok result) :
    result.lenA = lenA / 2 + first.lenA ∧
      0 < result.lenA ∧ result.lenB ≤ result.lenA := by
  have hlenLowA : Nat.min lenA (lenA / 2) ≤ lenA / 2 :=
    Nat.min_le_right _ _
  have hlenLowB : Nat.min lenB (lenA / 2) ≤ lenA / 2 :=
    Nat.min_le_right _ _
  have hleading := hgcdRecursiveFinalReconstruct_lenA_eq_of_invariant this A
    B T0 a b highA highB scratch (Nat.min lenA (lenA / 2))
    (Nat.min lenB (lenA / 2)) (lenA / 2) (lenA - lenA / 2) first result
    entries polyLowA polyLowB polyHighA polyHighB hcfg hp hinvariant
    hlenLowA hlenLowB physical hMatrix hLowA.toCanonicalSlice
    hLowB.toCanonicalSlice hHighA hHighB hrun
  have hrefines := hgcdRecursiveReconstructPair_refines this A B T0 a b
    highA highB scratch (Nat.min lenA (lenA / 2))
    (Nat.min lenB (lenA / 2)) first.lenA first.lenB (lenA / 2)
    first.matrix first.valid first.sgn first.heap result entries polyLowA
    polyLowB polyHighA polyHighB hcfg hp physical hMatrix
    hLowA.toCanonicalSlice hLowB.toCanonicalSlice
    hHighA hHighB hrun
  refine ⟨hleading.1, hleading.2, ?_⟩
  exact hgcdRecursiveFirstReconstruct_order lenA lenB
    (lenA - lenA / 2) (lenA / 2) first.lenA first.lenB
    (hgcdMatLen first.matrix first.valid (0 : Fin 4))
    (hgcdMatLen first.matrix first.valid (2 : Fin 4)) result.lenA
    result.lenB rfl hinvariant.aboveHalf hinvariant.stopped
    (by simpa [hgcdMatLen, hgcdMatLenRaw] using hinvariant.row0A)
    (by simpa [hgcdMatLen, hgcdMatLenRaw] using hinvariant.row2A)
    hleading.1 hrefines.2.2.2.1

/-- Complete proof-erased provider value consumed by
`hgcdRecursiveBodyBelow`.  Its four fields all come from one successful raw
paired reconstruction plus the real first-child length invariant. -/
theorem hgcdRecursiveFirstReconstruct_invariant_of_execution
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (A B T0 a b highA highB scratch : RawPtr UInt64)
    (lenA lenB : Nat) (first : HgcdRecursiveResult)
    (result : HgcdRecursiveReconstructPairResult)
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (polyLowA polyLowB polyHighA polyHighB :
      Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (hinputOrder : lenB < lenA)
    (hinvariant : HgcdRecursiveLengthInvariant (lenA - lenA / 2) first)
    (physical : HgcdRecursiveReconstructPairWorkspaceProvider this A B T0
      a b highA highB scratch (Nat.min lenA (lenA / 2))
      (Nat.min lenB (lenA / 2)) first.lenA first.lenB (lenA / 2)
      first.matrix first.valid first.sgn first.heap)
    (hMatrix : HgcdMatRawDenseRep this first.heap first.matrix entries
      first.valid)
    (hLowA : RawDensePolyRep this first.heap a
      (Nat.min lenA (lenA / 2)) polyLowA)
    (hLowB : RawDensePolyRep this first.heap b
      (Nat.min lenB (lenA / 2)) polyLowB)
    (hHighA : RawDensePolyRep this first.heap highA first.lenA polyHighA)
    (hHighB : RawDensePolyRep this first.heap highB first.lenB polyHighB)
    (hrun : hgcdRecursiveReconstructPair this A B T0 a b highA highB scratch
      (Nat.min lenA (lenA / 2)) (Nat.min lenB (lenA / 2)) first.lenA
      first.lenB (lenA / 2) first.matrix first.valid first.sgn first.heap =
        .ok result) :
    HgcdFirstReconstructionInvariant lenA first result := by
  have horder := hgcdRecursiveFirstReconstruct_order_of_invariant this A B
    T0 a b highA highB scratch lenA lenB first result entries polyLowA
    polyLowB polyHighA polyHighB hcfg hp hinvariant physical hMatrix hLowA
    hLowB hHighA hHighB hrun
  exact {
    leadingA := horder.1
    positiveA := horder.2.1
    ordered := horder.2.2
    decreases := hgcdRecursiveFirstReconstruct_bound_of_invariant this A B
      T0 a b highA highB scratch lenA lenB first result entries polyLowA
      polyLowB polyHighA polyHighB hcfg hp hinputOrder hinvariant physical
      hMatrix hLowA hLowB hHighA hHighB hrun }

/-- Physical workspace for the reconstruction that follows an actual first
dispatch result. -/
def HgcdRecursiveFirstReconstructWorkspaceProvider (this : DenseUPolyZp)
    (a b W scratch : RawPtr UInt64) (lenA lenB : Nat)
    (actualFirst : HgcdRecursiveResult → Prop) : Prop :=
  let ws := hgcdRecursiveWorkspace W lenA
  ∀ first, actualFirst first →
    HgcdRecursiveReconstructPairWorkspaceProvider this ws.a2 ws.b2 ws.T0 a b
      ws.a3 ws.b3 scratch (Nat.min lenA (lenA / 2))
      (Nat.min lenB (lenA / 2)) first.lenA first.lenB (lenA / 2)
      first.matrix first.valid first.sgn first.heap

/-- Construct the first-reconstruction bound from the actual child execution,
its induction invariant, and physical frame/workspace facts. -/
theorem hgcdFirstReconstructionBoundProvider_of_actual_dispatch
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (a b W scratch : RawPtr UInt64) (lenA lenB : Nat) (heap : RawHeap)
    (actualFirst : HgcdRecursiveResult → Prop)
    (lowPolyA lowPolyB : Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (horder : lenB < lenA)
    (hLowA : RawDensePolyRep this heap a (Nat.min lenA (lenA / 2)) lowPolyA)
    (hLowB : RawDensePolyRep this heap b (Nat.min lenB (lenA / 2)) lowPolyB)
    (frame : ∀ first, actualFirst first →
      HgcdRecursiveFirstDispatchFrame a b (Nat.min lenA (lenA / 2))
        (Nat.min lenB (lenA / 2)) heap first)
    (reconstructWorkspace :
      HgcdRecursiveFirstReconstructWorkspaceProvider this a b W scratch lenA
        lenB actualFirst)
    (firstInvariant : ∀ first, actualFirst first →
      ∃ inputHighA inputHighB outputHighA outputHighB entries,
        let ws := hgcdRecursiveWorkspace W lenA
        HgcdRecursiveRawInvariant this inputHighA inputHighB outputHighA
          outputHighB entries true ws.a3 ws.b3 (lenA - lenA / 2) first) :
    HgcdFirstReconstructionBoundProvider this a b W scratch lenA lenB
      actualFirst := by
  intro first reconstructed hactual hlength hreconstruct
  let ws := hgcdRecursiveWorkspace W lenA
  rcases firstInvariant first hactual with
    ⟨inputHighA, inputHighB, outputHighA, outputHighB, entries, hFirst⟩
  have hframe := frame first hactual
  have hLowAAfter : RawDensePolyRep this first.heap a
      (Nat.min lenA (lenA / 2)) lowPolyA :=
    rawDensePolyRep_of_same_prefix this heap first.heap a
      (Nat.min lenA (lenA / 2)) lowPolyA hframe.layout hframe.lowAPrefix hLowA
  have hLowBAfter : RawDensePolyRep this first.heap b
      (Nat.min lenB (lenA / 2)) lowPolyB :=
    rawDensePolyRep_of_same_prefix this heap first.heap b
      (Nat.min lenB (lenA / 2)) lowPolyB hframe.layout hframe.lowBPrefix hLowB
  exact hgcdRecursiveFirstReconstruct_invariant_of_execution this ws.a2 ws.b2
    ws.T0 a b ws.a3 ws.b3 scratch lenA lenB first reconstructed entries
    lowPolyA lowPolyB outputHighA outputHighB hcfg hp horder hlength
    (reconstructWorkspace first hactual) (hFirst.matrixSemantics rfl).1
    hLowAAfter hLowBAfter hFirst.aRep hFirst.bRep hreconstruct

/-- Assemble the body reconstruction provider directly from the real first
cutoff dispatch, its recursive induction hypothesis, and physical frames. -/
theorem hgcdFirstReconstructionBoundProvider_of_dispatch
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (bound : Nat) (recurse : HgcdRecursiveCallBelow bound)
    (a b W scratch : RawPtr UInt64) (lenA lenB : Nat) (heap : RawHeap)
    (inputHighA inputHighB lowPolyA lowPolyB :
      Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (horder : lenB < lenA)
    (iterWorkspace :
      let ws := hgcdRecursiveWorkspace W lenA
      let high := hgcdRecursiveHighInput a b lenA lenB
      HgcdRecursiveDispatchIterWorkspace this ws.R
        (hgcdRecursiveWorkspace_R_valid W lenA) ws.a3 ws.b3 high.a0 high.b0
        high.lenA0 high.lenB0 ws.q ws.W3 ws.T0 ws.T1 scratch ws.a2 heap
        inputHighA inputHighB)
    (recursiveRefines :
      let ws := hgcdRecursiveWorkspace W lenA
      let high := hgcdRecursiveHighInput a b lenA lenB
      ∀ (hchildOrder : high.lenB0 < high.lenA0)
        (hchildDecrease : high.lenA0 < bound),
      HgcdRecursiveCallbackRefinesAt this bound recurse ws.R
        (hgcdRecursiveWorkspace_R_valid W lenA) ws.a3 ws.b3 high.a0 high.b0
        high.lenA0 high.lenB0 ws.next scratch heap hchildOrder hchildDecrease
        inputHighA inputHighB)
    (frame :
      let ws := hgcdRecursiveWorkspace W lenA
      let high := hgcdRecursiveHighInput a b lenA lenB
      ∀ (hchildOrder : high.lenB0 < high.lenA0)
        (hchildDecrease : high.lenA0 < bound),
      HgcdRecursiveFirstDispatchFrameProvider this bound recurse ws.R
        (hgcdRecursiveWorkspace_R_valid W lenA) ws.a3 ws.b3 high.a0 high.b0
        high.lenA0 high.lenB0 ws.q ws.W3 ws.T0 ws.T1 scratch ws.a2 ws.next
        heap hchildOrder hchildDecrease a b (Nat.min lenA (lenA / 2))
        (Nat.min lenB (lenA / 2)))
    (reconstructWorkspace : HgcdRecursiveFirstReconstructWorkspaceProvider
      this a b W scratch lenA lenB
        (HgcdFirstDispatchResult this bound recurse a b W scratch lenA lenB
          heap))
    (hLowA : RawDensePolyRep this heap a (Nat.min lenA (lenA / 2)) lowPolyA)
    (hLowB : RawDensePolyRep this heap b (Nat.min lenB (lenA / 2)) lowPolyB) :
    HgcdFirstReconstructionBoundProvider this a b W scratch lenA lenB
      (HgcdFirstDispatchResult this bound recurse a b W scratch lenA lenB
        heap) := by
  let ws := hgcdRecursiveWorkspace W lenA
  let high := hgcdRecursiveHighInput a b lenA lenB
  apply hgcdFirstReconstructionBoundProvider_of_actual_dispatch this a b W
    scratch lenA lenB heap
      (HgcdFirstDispatchResult this bound recurse a b W scratch lenA lenB heap)
      lowPolyA lowPolyB hcfg hp horder hLowA hLowB
  · intro first hactual
    apply hgcdRecursiveFirstDispatchResult_frame this bound recurse ws.R
      (hgcdRecursiveWorkspace_R_valid W lenA) ws.a3 ws.b3 high.a0 high.b0
      high.lenA0 high.lenB0 ws.q ws.W3 ws.T0 ws.T1 scratch ws.a2 ws.next
      heap a b (Nat.min lenA (lenA / 2)) (Nat.min lenB (lenA / 2)) frame
      first
    simpa [HgcdFirstDispatchResult, ws, high] using hactual
  · exact reconstructWorkspace
  · intro first hactual
    rcases hgcdRecursiveDispatchResult_rawInvariant this bound recurse ws.R
        (hgcdRecursiveWorkspace_R_valid W lenA) ws.a3 ws.b3 high.a0 high.b0
        high.lenA0 high.lenB0 ws.q ws.W3 ws.T0 ws.T1 scratch ws.a2 ws.next
        heap inputHighA inputHighB hcfg hp iterWorkspace recursiveRefines first
        (by simpa [HgcdFirstDispatchResult, ws, high] using hactual) with
      ⟨outputHighA, outputHighB, entries, hInvariant⟩
    exact ⟨inputHighA, inputHighB, outputHighA, outputHighB, entries, by
      simpa [high, hgcdRecursiveHighInput] using hInvariant⟩

/-- Admissible physical/representation workspace for the first child call of
one recursive body.  Recursive semantics are deliberately not stored here;
they are supplied only by the well-founded induction hypothesis. -/
structure HgcdRecursiveFirstCallWorkspace (this : DenseUPolyZp)
    (bound : Nat) (recurse : HgcdRecursiveCallBelow bound)
    (a b W scratch : RawPtr UInt64) (lenA lenB : Nat) (heap : RawHeap)
    (inputHighA inputHighB lowPolyA lowPolyB :
      Polynomial (ZMod this._p.toNat)) : Prop where
  iter :
    let ws := hgcdRecursiveWorkspace W lenA
    let high := hgcdRecursiveHighInput a b lenA lenB
    HgcdRecursiveDispatchIterWorkspace this ws.R
      (hgcdRecursiveWorkspace_R_valid W lenA) ws.a3 ws.b3 high.a0 high.b0
      high.lenA0 high.lenB0 ws.q ws.W3 ws.T0 ws.T1 scratch ws.a2 heap
      inputHighA inputHighB
  frame :
    let ws := hgcdRecursiveWorkspace W lenA
    let high := hgcdRecursiveHighInput a b lenA lenB
    ∀ (hchildOrder : high.lenB0 < high.lenA0)
      (hchildDecrease : high.lenA0 < bound),
    HgcdRecursiveFirstDispatchFrameProvider this bound recurse ws.R
      (hgcdRecursiveWorkspace_R_valid W lenA) ws.a3 ws.b3 high.a0 high.b0
      high.lenA0 high.lenB0 ws.q ws.W3 ws.T0 ws.T1 scratch ws.a2 ws.next
      heap hchildOrder hchildDecrease a b (Nat.min lenA (lenA / 2))
      (Nat.min lenB (lenA / 2))
  reconstruct : HgcdRecursiveFirstReconstructWorkspaceProvider this a b W
    scratch lenA lenB
      (HgcdFirstDispatchResult this bound recurse a b W scratch lenA lenB heap)
  lowA : RawDensePolyRep this heap a (Nat.min lenA (lenA / 2)) lowPolyA
  lowB : RawDensePolyRep this heap b (Nat.min lenB (lenA / 2)) lowPolyB

/-- First-call data available at one well-founded body invocation.  The
semantic field is precisely the smaller-call induction hypothesis. -/
structure HgcdRecursiveFirstCallAdmissible (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (bound : Nat) (recurse : HgcdRecursiveCallBelow bound)
    (a b W scratch : RawPtr UInt64) (lenA lenB : Nat) (heap : RawHeap)
    (inputHighA inputHighB lowPolyA lowPolyB :
      Polynomial (ZMod this._p.toNat)) : Prop where
  workspace : HgcdRecursiveFirstCallWorkspace this bound recurse a b W
    scratch lenA lenB heap inputHighA inputHighB lowPolyA lowPolyB
  recursiveRefines :
    let ws := hgcdRecursiveWorkspace W lenA
    let high := hgcdRecursiveHighInput a b lenA lenB
    ∀ (hchildOrder : high.lenB0 < high.lenA0)
      (hchildDecrease : high.lenA0 < bound),
    HgcdRecursiveCallbackRefinesAt this bound recurse ws.R
      (hgcdRecursiveWorkspace_R_valid W lenA) ws.a3 ws.b3 high.a0 high.b0
      high.lenA0 high.lenB0 ws.next scratch heap hchildOrder hchildDecrease
      inputHighA inputHighB

/-- Attach only the well-founded child semantic theorem to an already
established physical first-call workspace. -/
theorem HgcdRecursiveFirstCallWorkspace.admissible (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (bound : Nat) (recurse : HgcdRecursiveCallBelow bound)
    (a b W scratch : RawPtr UInt64) (lenA lenB : Nat) (heap : RawHeap)
    (inputHighA inputHighB lowPolyA lowPolyB :
      Polynomial (ZMod this._p.toNat))
    (workspace : HgcdRecursiveFirstCallWorkspace this bound recurse a b W
      scratch lenA lenB heap inputHighA inputHighB lowPolyA lowPolyB)
    (recursiveRefines :
      let ws := hgcdRecursiveWorkspace W lenA
      let high := hgcdRecursiveHighInput a b lenA lenB
      ∀ (hchildOrder : high.lenB0 < high.lenA0)
        (hchildDecrease : high.lenA0 < bound),
      HgcdRecursiveCallbackRefinesAt this bound recurse ws.R
        (hgcdRecursiveWorkspace_R_valid W lenA) ws.a3 ws.b3 high.a0 high.b0
        high.lenA0 high.lenB0 ws.next scratch heap hchildOrder hchildDecrease
        inputHighA inputHighB) :
    HgcdRecursiveFirstCallAdmissible this bound recurse a b W scratch lenA
      lenB heap inputHighA inputHighB lowPolyA lowPolyB :=
  ⟨workspace, recursiveRefines⟩

/-- A first-call admissible workspace plus the child induction hypothesis
supplies exactly both proof arguments consumed by `hgcdRecursiveBodyBelow`. -/
theorem hgcdRecursiveFirstCall_providers (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (bound : Nat) (recurse : HgcdRecursiveCallBelow bound)
    (a b W scratch : RawPtr UInt64) (lenA lenB : Nat) (heap : RawHeap)
    (inputHighA inputHighB lowPolyA lowPolyB :
      Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (horder : lenB < lenA)
    (workspace : HgcdRecursiveFirstCallWorkspace this bound recurse a b W
      scratch lenA lenB heap inputHighA inputHighB lowPolyA lowPolyB)
    (recursiveRefines :
      let ws := hgcdRecursiveWorkspace W lenA
      let high := hgcdRecursiveHighInput a b lenA lenB
      ∀ (hchildOrder : high.lenB0 < high.lenA0)
        (hchildDecrease : high.lenA0 < bound),
      HgcdRecursiveCallbackRefinesAt this bound recurse ws.R
        (hgcdRecursiveWorkspace_R_valid W lenA) ws.a3 ws.b3 high.a0 high.b0
        high.lenA0 high.lenB0 ws.next scratch heap hchildOrder hchildDecrease
        inputHighA inputHighB) :
    (∀ first,
      let ws := hgcdRecursiveWorkspace W lenA
      let high := hgcdRecursiveHighInput a b lenA lenB
      ∀ (hchildOrder : high.lenB0 < high.lenA0)
        (hchildDecrease : high.lenA0 < bound),
      hgcdRecursiveDispatchBelow this bound recurse ws.R
        (hgcdRecursiveWorkspace_R_valid W lenA) ws.a3 ws.b3 high.a0 high.b0
        high.lenA0 high.lenB0 ws.q ws.W3 ws.T0 ws.T1 scratch ws.a2 ws.next
        heap hchildOrder hchildDecrease = .ok first →
      HgcdRecursiveLengthInvariant high.lenA0 first) ∧
    HgcdFirstReconstructionBoundProvider this a b W scratch lenA lenB
      (HgcdFirstDispatchResult this bound recurse a b W scratch lenA lenB
        heap) := by
  let ws := hgcdRecursiveWorkspace W lenA
  let high := hgcdRecursiveHighInput a b lenA lenB
  constructor
  · intro first
    dsimp only
    intro hchildOrder hchildDecrease hrun
    exact hgcdRecursiveDispatchBelow_lengthInvariant this bound recurse ws.R
      (hgcdRecursiveWorkspace_R_valid W lenA) ws.a3 ws.b3 high.a0 high.b0
      high.lenA0 high.lenB0 ws.q ws.W3 ws.T0 ws.T1 scratch ws.a2 ws.next
      heap inputHighA inputHighB hcfg hp workspace.iter recursiveRefines first
      hchildOrder hchildDecrease hrun
  · exact hgcdFirstReconstructionBoundProvider_of_dispatch this bound
      recurse a b W scratch lenA lenB heap inputHighA inputHighB lowPolyA
      lowPolyB hcfg hp horder workspace.iter recursiveRefines workspace.frame
      workspace.reconstruct workspace.lowA workspace.lowB

/-- The actual first dispatch simultaneously supplies its recursive semantic
invariant and transports both source low prefixes into the returned heap. -/
theorem hgcdRecursiveFirstCall_refines (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (bound : Nat) (recurse : HgcdRecursiveCallBelow bound)
    (a b W scratch : RawPtr UInt64) (lenA lenB : Nat) (heap : RawHeap)
    (inputHighA inputHighB lowPolyA lowPolyB :
      Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (admissible : HgcdRecursiveFirstCallAdmissible this bound recurse a b W
      scratch lenA lenB heap inputHighA inputHighB lowPolyA lowPolyB)
    (first : HgcdRecursiveResult)
    (hchildOrder :
      (hgcdRecursiveHighInput a b lenA lenB).lenB0 <
        (hgcdRecursiveHighInput a b lenA lenB).lenA0)
    (hchildDecrease :
      (hgcdRecursiveHighInput a b lenA lenB).lenA0 < bound)
    (hrun :
      let ws := hgcdRecursiveWorkspace W lenA
      let high := hgcdRecursiveHighInput a b lenA lenB
      hgcdRecursiveDispatchBelow this bound recurse ws.R
        (hgcdRecursiveWorkspace_R_valid W lenA) ws.a3 ws.b3 high.a0 high.b0
        high.lenA0 high.lenB0 ws.q ws.W3 ws.T0 ws.T1 scratch ws.a2 ws.next
        heap hchildOrder hchildDecrease = .ok first) :
    ∃ outputHighA outputHighB entries,
      let ws := hgcdRecursiveWorkspace W lenA
      let high := hgcdRecursiveHighInput a b lenA lenB
      HgcdRecursiveRawInvariant this inputHighA inputHighB outputHighA
          outputHighB entries true ws.a3 ws.b3 high.lenA0 first ∧
        RawDensePolyRep this first.heap a (Nat.min lenA (lenA / 2)) lowPolyA ∧
        RawDensePolyRep this first.heap b (Nat.min lenB (lenA / 2)) lowPolyB := by
  let ws := hgcdRecursiveWorkspace W lenA
  let high := hgcdRecursiveHighInput a b lenA lenB
  rcases hgcdRecursiveDispatchBelow_rawInvariant this bound recurse ws.R
      (hgcdRecursiveWorkspace_R_valid W lenA) ws.a3 ws.b3 high.a0 high.b0
      high.lenA0 high.lenB0 ws.q ws.W3 ws.T0 ws.T1 scratch ws.a2 ws.next
      heap inputHighA inputHighB hcfg hp hchildOrder hchildDecrease
      admissible.workspace.iter
      (admissible.recursiveRefines hchildOrder hchildDecrease) first hrun with
    ⟨outputHighA, outputHighB, entries, hInvariant⟩
  refine ⟨outputHighA, outputHighB, entries, ?_⟩
  dsimp only
  refine ⟨hInvariant, ?_⟩
  have hframe : HgcdRecursiveFirstDispatchFrame a b
      (Nat.min lenA (lenA / 2)) (Nat.min lenB (lenA / 2)) heap first := by
    exact admissible.workspace.frame hchildOrder hchildDecrease first hrun
  have hLowAAfter : RawDensePolyRep this first.heap a
      (Nat.min lenA (lenA / 2)) lowPolyA :=
    rawDensePolyRep_of_same_prefix this heap first.heap a
    (Nat.min lenA (lenA / 2)) lowPolyA hframe.layout hframe.lowAPrefix
    admissible.workspace.lowA
  have hLowBAfter : RawDensePolyRep this first.heap b
      (Nat.min lenB (lenA / 2)) lowPolyB :=
    rawDensePolyRep_of_same_prefix this heap first.heap b
    (Nat.min lenB (lenA / 2)) lowPolyB hframe.layout hframe.lowBPrefix
    admissible.workspace.lowB
  exact ⟨hLowAAfter, hLowBAfter⟩

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

/-- Erasing the strict-decrease arguments from one cutoff dispatch gives the
same generated execution as the source-shaped callback dispatch, provided
the two callbacks agree on that actual recursive call. -/
theorem hgcdRecursiveDispatchBelow_eq_dispatch (this : DenseUPolyZp)
    (bound : Nat) (below : HgcdRecursiveCallBelow bound)
    (plain : HgcdRecursiveCall)
    (matrix : HgcdMat) (hMatrix : matrix.Valid)
    (a3 b3 inputA inputB : RawPtr UInt64) (lenInputA lenInputB : Nat)
    (Q : RawPtr UInt64) (W3 : RawPtr Word3)
    (T0 T1 scratch stage WNext : RawPtr UInt64) (heap : RawHeap)
    (horder : lenInputB < lenInputA) (hdecrease : lenInputA < bound)
    (hagrees : below matrix hMatrix true a3 b3 inputA inputB lenInputA
      lenInputB WNext scratch heap horder hdecrease =
        plain matrix hMatrix true a3 b3 inputA inputB lenInputA lenInputB
          WNext scratch heap) :
    hgcdRecursiveDispatchBelow this bound below matrix hMatrix a3 b3 inputA
        inputB lenInputA lenInputB Q W3 T0 T1 scratch stage WNext heap horder
        hdecrease =
      hgcdRecursiveDispatch this plain matrix hMatrix a3 b3 inputA inputB
        lenInputA lenInputB Q W3 T0 T1 scratch stage WNext heap := by
  simp only [hgcdRecursiveDispatchBelow, hgcdRecursiveDispatch]
  split
  · rfl
  · exact hagrees

/-- Proof-erased observation of a recursive execution.  Faults are retained
exactly, and every computational field of a successful result is retained. -/
def hgcdRecursiveExecValue (run : RawExec HgcdRecursiveResult) :
    RawExec HgcdRecursiveValue :=
  run.map HgcdRecursiveResult.value

/-- Equality of proof-erased recursive executions is enough for equality of
the original executions because `HgcdRecursiveResult.valid` is proof-only. -/
theorem hgcdRecursiveExec_ext_value
    (left right : RawExec HgcdRecursiveResult)
    (hvalue : hgcdRecursiveExecValue left = hgcdRecursiveExecValue right) :
    left = right := by
  cases left with
  | error leftFault =>
      cases right with
      | error rightFault =>
          simp only [hgcdRecursiveExecValue, Except.map] at hvalue
          exact congrArg Except.error (Except.error.inj hvalue)
      | ok rightResult =>
          simp only [hgcdRecursiveExecValue, Except.map] at hvalue
          contradiction
  | ok leftResult =>
      cases right with
      | error rightFault =>
          simp only [hgcdRecursiveExecValue, Except.map] at hvalue
          contradiction
      | ok rightResult =>
          simp only [hgcdRecursiveExecValue, Except.map, Except.ok.injEq]
            at hvalue
          exact congrArg Except.ok
            (HgcdRecursiveResult.ext_value leftResult rightResult hvalue)

/-- The proof-indexed well-founded body is the exact source-shaped body after
erasing only its strict-decrease proofs. -/
theorem hgcdRecursiveBodyBelow_eq_body (this : DenseUPolyZp)
    (bound : Nat) (below : HgcdRecursiveCallBelow bound)
    (plain : HgcdRecursiveCall)
    (M : HgcdMat) (hM : M.Valid) (computeM : Bool)
    (A B a b : RawPtr UInt64) (lenA lenB : Nat)
    (W scratch : RawPtr UInt64) (heap : RawHeap)
    (hbound : lenA = bound) (horder : lenB < lenA)
    (firstLength : ¬ lenB < lenA / 2 + 1 → ∀ first,
      let ws := hgcdRecursiveWorkspace W lenA
      let high := hgcdRecursiveHighInput a b lenA lenB
      ∀ (hchildOrder : high.lenB0 < high.lenA0)
        (hchildDecrease : high.lenA0 < bound),
      hgcdRecursiveDispatchBelow this bound below ws.R
        (hgcdRecursiveWorkspace_R_valid W lenA) ws.a3 ws.b3 high.a0 high.b0
        high.lenA0 high.lenB0 ws.q ws.W3 ws.T0 ws.T1 scratch ws.a2 ws.next
        heap hchildOrder hchildDecrease = .ok first →
      HgcdRecursiveLengthInvariant high.lenA0 first)
    (reconstructionBound : ¬ lenB < lenA / 2 + 1 →
      HgcdFirstReconstructionBoundProvider this a b W
      scratch lenA lenB (HgcdFirstDispatchResult this bound below a b W
        scratch lenA lenB heap))
    (hagrees : ∀ (matrix : HgcdMat) (hMatrix : matrix.Valid)
      (a3 b3 inputA inputB : RawPtr UInt64) (lenInputA lenInputB : Nat)
      (WNext childScratch : RawPtr UInt64) (childHeap : RawHeap)
      (hchildOrder : lenInputB < lenInputA)
      (hchildDecrease : lenInputA < bound),
      below matrix hMatrix true a3 b3 inputA inputB lenInputA lenInputB
          WNext childScratch childHeap hchildOrder hchildDecrease =
        plain matrix hMatrix true a3 b3 inputA inputB lenInputA lenInputB
          WNext childScratch childHeap) :
    hgcdRecursiveBodyBelow this bound below M hM computeM A B a b lenA lenB
        W scratch heap hbound horder firstLength reconstructionBound =
      hgcdRecursiveBody this plain M hM computeM A B a b lenA lenB W scratch
        heap := by
  apply hgcdRecursiveExec_ext_value
  rw [hgcdRecursiveBodyBelow, hgcdRecursiveBody]
  split
  · rfl
  · simp only [hgcdRecursiveExecValue, hgcdRecursiveDispatchBelow,
      hgcdRecursiveDispatch, hagrees]
    split <;> simp_all
    split <;> simp_all
    split <;> simp_all
    split <;> simp_all

/-- Invoke the genuine well-founded body from admissible first-call data.
Neither of the body's semantic proof providers is accepted from the caller. -/
def hgcdRecursiveBodyAdmissible (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (bound : Nat) (recurse : HgcdRecursiveCallBelow bound)
    (M : HgcdMat) (hM : M.Valid) (computeM : Bool)
    (A B a b : RawPtr UInt64) (lenA lenB : Nat)
    (W scratch : RawPtr UInt64) (heap : RawHeap)
    (inputHighA inputHighB lowPolyA lowPolyB :
      Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (hbound : lenA = bound) (horder : lenB < lenA)
    (admissible : HgcdRecursiveFirstCallAdmissible this bound recurse a b W
      scratch lenA lenB heap inputHighA inputHighB lowPolyA lowPolyB) :
    RawExec HgcdRecursiveResult :=
  let providers := hgcdRecursiveFirstCall_providers this bound recurse a b W
    scratch lenA lenB heap inputHighA inputHighB lowPolyA lowPolyB hcfg hp
    horder admissible.workspace admissible.recursiveRefines
  hgcdRecursiveBodyBelow this bound recurse M hM computeM A B a b lenA lenB
    W scratch heap hbound horder (fun _ => providers.1)
      (fun _ => providers.2)

/-- The admissible wrapper erases to the exact generated source body whenever
the well-founded child callback erases to the plain recursive callback. -/
theorem hgcdRecursiveBodyAdmissible_eq_body (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (bound : Nat) (below : HgcdRecursiveCallBelow bound)
    (plain : HgcdRecursiveCall)
    (M : HgcdMat) (hM : M.Valid) (computeM : Bool)
    (A B a b : RawPtr UInt64) (lenA lenB : Nat)
    (W scratch : RawPtr UInt64) (heap : RawHeap)
    (inputHighA inputHighB lowPolyA lowPolyB :
      Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (hbound : lenA = bound) (horder : lenB < lenA)
    (admissible : HgcdRecursiveFirstCallAdmissible this bound below a b W
      scratch lenA lenB heap inputHighA inputHighB lowPolyA lowPolyB)
    (hagrees : ∀ (matrix : HgcdMat) (hMatrix : matrix.Valid)
      (a3 b3 inputA inputB : RawPtr UInt64) (lenInputA lenInputB : Nat)
      (WNext childScratch : RawPtr UInt64) (childHeap : RawHeap)
      (hchildOrder : lenInputB < lenInputA)
      (hchildDecrease : lenInputA < bound),
      below matrix hMatrix true a3 b3 inputA inputB lenInputA lenInputB
          WNext childScratch childHeap hchildOrder hchildDecrease =
        plain matrix hMatrix true a3 b3 inputA inputB lenInputA lenInputB
          WNext childScratch childHeap) :
    hgcdRecursiveBodyAdmissible this bound below M hM computeM A B a b lenA
        lenB W scratch heap inputHighA inputHighB lowPolyA lowPolyB hcfg hp
        hbound horder admissible =
      hgcdRecursiveBody this plain M hM computeM A B a b lenA lenB W scratch
        heap := by
  unfold hgcdRecursiveBodyAdmissible
  exact hgcdRecursiveBodyBelow_eq_body this bound below plain M hM computeM A
    B a b lenA lenB W scratch heap hbound horder _ _ hagrees

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

/-- Branch-local admissibility changes no generated computation.  Erasing the
strict-decrease proofs yields exactly the source-shaped recursive body on both
the base and non-base paths. -/
theorem hgcdRecursiveBodyBranchAdmissible_eq_body (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (bound : Nat) (below : HgcdRecursiveCallBelow bound)
    (plain : HgcdRecursiveCall)
    (M : HgcdMat) (hM : M.Valid) (computeM : Bool)
    (A B a b : RawPtr UInt64) (lenA lenB : Nat)
    (W scratch : RawPtr UInt64) (heap : RawHeap)
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (hbound : lenA = bound) (horder : lenB < lenA)
    (nonBase : ¬ lenB < lenA / 2 + 1 →
      HgcdRecursiveNonBasePackage this bound below a b W scratch lenA lenB
        heap)
    (recursiveRefines : ∀ hbase,
      HgcdRecursiveNonBaseCallbackRefines this bound below a b W scratch lenA
        lenB heap (nonBase hbase))
    (hagrees : ∀ (matrix : HgcdMat) (hMatrix : matrix.Valid)
      (a3 b3 inputA inputB : RawPtr UInt64) (lenInputA lenInputB : Nat)
      (WNext childScratch : RawPtr UInt64) (childHeap : RawHeap)
      (hchildOrder : lenInputB < lenInputA)
      (hchildDecrease : lenInputA < bound),
      below matrix hMatrix true a3 b3 inputA inputB lenInputA lenInputB
          WNext childScratch childHeap hchildOrder hchildDecrease =
        plain matrix hMatrix true a3 b3 inputA inputB lenInputA lenInputB
          WNext childScratch childHeap) :
    hgcdRecursiveBodyBranchAdmissible this bound below M hM computeM A B a b
        lenA lenB W scratch heap hcfg hp hbound horder nonBase
        recursiveRefines =
      hgcdRecursiveBody this plain M hM computeM A B a b lenA lenB W scratch
        heap := by
  unfold hgcdRecursiveBodyBranchAdmissible
  split
  next hbase =>
    exact hgcdRecursiveBodyBelow_eq_body this bound below plain M hM computeM
      A B a b lenA lenB W scratch heap hbound horder
      (fun hnonbase => False.elim (hnonbase hbase))
      (fun hnonbase => False.elim (hnonbase hbase)) hagrees
  next hbase =>
    let package := nonBase hbase
    let providers := hgcdRecursiveFirstCall_providers this bound below a b W
      scratch lenA lenB heap package.inputHighA package.inputHighB
      package.lowPolyA package.lowPolyB hcfg hp horder
      package.workspace (recursiveRefines hbase)
    exact hgcdRecursiveBodyBelow_eq_body this bound below plain M hM computeM
      A B a b lenA lenB W scratch heap hbound horder
      (fun _ => providers.1) (fun _ => providers.2) hagrees

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

end CLPoly.Impl.StrictHGCDRawRefinement
