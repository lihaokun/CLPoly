import CLPoly.Generated.StrictMul
import CLPoly.Impl.StrictDivremRefinement
import CLPoly.Impl.StrictEuclidRefinement
import CLPoly.Impl.StrictPolyAddSubRefinement

set_option autoImplicit false

namespace CLPoly.Impl.StrictMulRefinement

open Generated.StrictMul
open Generated.StrictPolyAddSub
open CLPoly.Impl.StrictWordArithmetic
open CLPoly.Impl.StrictDivremRefinement
open CLPoly.Impl.RawPolynomialRep
open CLPoly.Impl.StrictEuclidRefinement
open CLPoly.Impl.StrictPolyAddSubRefinement

theorem rawPtr_add_add (ptr : RawPtr UInt64) (first second : Nat) :
    (ptr.add first).add second = ptr.add (first + second) := by
  cases ptr
  have hwidth : RawLimbWidth.width UInt64 = 1 := rfl
  simp [RawPtr.add, hwidth]
  omega

/-- Address-level non-aliasing for two UInt64 slices.  Unlike comparing
allocation regions, this also represents disjoint subarrays of one C++
scratch allocation. -/
def U64SlicesDisjoint (left : RawPtr UInt64) (leftLength : Nat)
    (right : RawPtr UInt64) (rightLength : Nat) : Prop :=
  ∀ i, i < leftLength → ∀ j, j < rightLength →
    left.region ≠ right.region ∨
      left.limbOffset + i ≠ right.limbOffset + j

theorem u64SlicesDisjoint_symm {left right : RawPtr UInt64}
    {leftLength rightLength : Nat}
    (h : U64SlicesDisjoint left leftLength right rightLength) :
    U64SlicesDisjoint right rightLength left leftLength := by
  intro j hj i hi
  rcases h i hi j hj with hregion | hoffset
  · exact Or.inl (Ne.symm hregion)
  · exact Or.inr (Ne.symm hoffset)

theorem u64SlicesDisjoint_mono {left right : RawPtr UInt64}
    {leftLength rightLength smallLeft smallRight : Nat}
    (h : U64SlicesDisjoint left leftLength right rightLength)
    (hleft : smallLeft ≤ leftLength) (hright : smallRight ≤ rightLength) :
    U64SlicesDisjoint left smallLeft right smallRight := by
  intro i hi j hj
  exact h i (by omega) j (by omega)

theorem u64SlicesDisjoint_of_region_ne {left right : RawPtr UInt64}
    {leftLength rightLength : Nat} (hregion : left.region ≠ right.region) :
    U64SlicesDisjoint left leftLength right rightLength := by
  intro _ _ _ _
  exact Or.inl hregion

theorem u64SlicesDisjoint_adjacent (base : RawPtr UInt64)
    (leftLength rightLength : Nat) :
    U64SlicesDisjoint base leftLength (base.add leftLength) rightLength := by
  intro i hi j hj
  right
  have hwidth : RawLimbWidth.width UInt64 = 1 := rfl
  simp [RawPtr.add, hwidth]
  omega

theorem u64SlicesDisjoint_add_of_le (base : RawPtr UInt64)
    (leftStart leftLength rightStart rightLength : Nat)
    (hle : leftStart + leftLength ≤ rightStart) :
    U64SlicesDisjoint (base.add leftStart) leftLength
      (base.add rightStart) rightLength := by
  intro i hi j hj
  right
  have hwidth : RawLimbWidth.width UInt64 = 1 := rfl
  simp [RawPtr.add, hwidth]
  omega

/-- Pairwise non-aliasing of the five consecutive regions used by one
Karatsuba frame, expressed at their canonical offsets in the shared scratch
allocation. -/
theorem karScratchRanges_pairwise (scratch : RawPtr UInt64)
    (m h recLength : Nat) :
    let t1 := scratch
    let t2 := scratch.add h
    let sP0 := scratch.add (2 * h)
    let sP1 := scratch.add (2 * h + (2 * m - 1))
    let recScratch := scratch.add
      (2 * h + (2 * m - 1) + (2 * h - 1))
    U64SlicesDisjoint t1 h t2 h ∧
      U64SlicesDisjoint t1 h sP0 (2 * m - 1) ∧
      U64SlicesDisjoint t1 h sP1 (2 * h - 1) ∧
      U64SlicesDisjoint t2 h sP0 (2 * m - 1) ∧
      U64SlicesDisjoint t2 h sP1 (2 * h - 1) ∧
      U64SlicesDisjoint sP0 (2 * m - 1) sP1 (2 * h - 1) ∧
      U64SlicesDisjoint sP0 (2 * m - 1) recScratch recLength ∧
      U64SlicesDisjoint sP1 (2 * h - 1) recScratch recLength := by
  dsimp
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · exact u64SlicesDisjoint_add_of_le scratch 0 h h h (by omega)
  · exact u64SlicesDisjoint_add_of_le scratch 0 h (2 * h)
      (2 * m - 1) (by omega)
  · exact u64SlicesDisjoint_add_of_le scratch 0 h
      (2 * h + (2 * m - 1)) (2 * h - 1) (by omega)
  · exact u64SlicesDisjoint_add_of_le scratch h h (2 * h)
      (2 * m - 1) (by omega)
  · exact u64SlicesDisjoint_add_of_le scratch h h
      (2 * h + (2 * m - 1)) (2 * h - 1) (by omega)
  · exact u64SlicesDisjoint_add_of_le scratch (2 * h) (2 * m - 1)
      (2 * h + (2 * m - 1)) (2 * h - 1) (by omega)
  · exact u64SlicesDisjoint_add_of_le scratch (2 * h) (2 * m - 1)
      (2 * h + (2 * m - 1) + (2 * h - 1))
      recLength (by omega)
  · exact u64SlicesDisjoint_add_of_le scratch
      (2 * h + (2 * m - 1)) (2 * h - 1)
      (2 * h + (2 * m - 1) + (2 * h - 1))
      recLength (by omega)

/-- The same pairwise partition, normalized to the nested pointer additions
that occur literally in generated `_kar_mul`. -/
theorem karScratchNested_pairwise (scratch : RawPtr UInt64)
    (m h recLength : Nat) :
    let t1 := scratch
    let t2 := t1.add h
    let sP0 := t2.add h
    let sP1 := sP0.add (2 * m - 1)
    let recScratch := sP1.add (2 * h - 1)
    U64SlicesDisjoint t1 h t2 h ∧
      U64SlicesDisjoint t1 h sP0 (2 * m - 1) ∧
      U64SlicesDisjoint t1 h sP1 (2 * h - 1) ∧
      U64SlicesDisjoint t2 h sP0 (2 * m - 1) ∧
      U64SlicesDisjoint t2 h sP1 (2 * h - 1) ∧
      U64SlicesDisjoint sP0 (2 * m - 1) sP1 (2 * h - 1) ∧
      U64SlicesDisjoint sP0 (2 * m - 1) recScratch recLength ∧
      U64SlicesDisjoint sP1 (2 * h - 1) recScratch recLength := by
  simpa only [rawPtr_add_add, two_mul, Nat.add_assoc] using
    karScratchRanges_pairwise scratch m h recLength

theorem u64SlicesDisjoint_add_left {base guard : RawPtr UInt64}
    {length guardLength start count : Nat}
    (h : U64SlicesDisjoint base length guard guardLength)
    (hrange : start + count ≤ length) :
    U64SlicesDisjoint (base.add start) count guard guardLength := by
  intro i hi j hj
  have hbase := h (start + i) (by omega) j hj
  rcases hbase with hregion | hoffset
  · exact Or.inl (by simpa [RawPtr.add] using hregion)
  · right
    have hwidth : RawLimbWidth.width UInt64 = 1 := rfl
    simp [RawPtr.add, hwidth] at hoffset ⊢
    omega

theorem sameU64Prefix_trans {heap1 heap2 heap3 : RawHeap}
    {ptr : RawPtr UInt64} {length : Nat}
    (h12 : SameU64Prefix heap1 heap2 ptr length)
    (h23 : SameU64Prefix heap2 heap3 ptr length) :
    SameU64Prefix heap1 heap3 ptr length := by
  intro i value hi hread
  exact h23 i value hi (h12 i value hi hread)

theorem writeU64_preserves_prefix (heap heap' : RawHeap)
    (dst guard : RawPtr UInt64) (dstLength guardLength writeIndex : Nat)
    (value : UInt64)
    (hdisjoint : U64SlicesDisjoint dst dstLength guard guardLength)
    (hwriteIndex : writeIndex < dstLength)
    (hwrite : heap.writeU64 dst writeIndex value = .ok heap') :
    SameU64Prefix heap heap' guard guardLength := by
  intro readIndex old hreadIndex hread
  exact RawHeap.readU64_writeU64_ne heap heap' dst guard writeIndex
    readIndex value old hwrite hread
    (hdisjoint writeIndex hwriteIndex readIndex hreadIndex)

theorem copyU64_preserves_prefix (heap heap' : RawHeap)
    (dst src guard : RawPtr UInt64) (count guardLength : Nat)
    (hDst : heap.ValidU64Slice dst count)
    (hSrc : heap.ValidU64Slice src count)
    (hdisjoint : U64SlicesDisjoint dst count guard guardLength)
    (hcopy : heap.copyU64 dst src count = .ok heap') :
    SameU64Prefix heap heap' guard guardLength := by
  intro readIndex old hreadIndex hread
  exact copyU64_preserves_read heap heap' dst src guard count readIndex old
    hDst hSrc hread
    (by
      intro writeIndex hwriteIndex
      exact hdisjoint writeIndex hwriteIndex readIndex hreadIndex)
    hcopy

/-- Content semantics for the actual recursive copy under address-level
non-overlap.  Unlike the older region-level lemma, this also covers adjacent
sub-slices of one Karatsuba scratch allocation. -/
theorem copyU64_refines_disjoint (heap : RawHeap)
    (dst src : RawPtr UInt64) (count : Nat)
    (hDst : heap.ValidU64Slice dst count)
    (hSrc : heap.ValidU64Slice src count)
    (hdisjoint : U64SlicesDisjoint dst count src count) :
    ∃ heap', heap.copyU64 dst src count = .ok heap' ∧
      RawHeap.SameLayout heap heap' ∧
      CopyU64Contents heap heap' dst src count := by
  cases count with
  | zero =>
    exact ⟨heap, rfl, fun _ _ => Iff.rfl, by intro k _ hk; omega⟩
  | succ n =>
    simp only [RawHeap.copyU64]
    rcases heap.readU64_of_valid src (n + 1) 0 hSrc (by omega) with
      ⟨value, hread⟩
    simp only [hread]
    rcases heap.writeU64_of_valid dst (n + 1) 0 value hDst (by omega) with
      ⟨heap1, hwrite⟩
    simp only [hwrite]
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 dst 0 value hwrite
    have hDstTail0 := heap.validU64Slice_add dst (n + 1) 1 n hDst (by omega)
    have hSrcTail0 := heap.validU64Slice_add src (n + 1) 1 n hSrc (by omega)
    have hDstTail1 := (hlayout1 (dst.add 1) n).mp hDstTail0
    have hSrcTail1 := (hlayout1 (src.add 1) n).mp hSrcTail0
    have htailDisjoint : U64SlicesDisjoint (dst.add 1) n (src.add 1) n := by
      intro i hi j hj
      simpa [RawPtr.add, Nat.add_assoc, Nat.add_comm, Nat.add_left_comm] using
        hdisjoint (i + 1) (by omega) (j + 1) (by omega)
    rcases copyU64_refines_disjoint heap1 (dst.add 1) (src.add 1) n
      hDstTail1 hSrcTail1 htailDisjoint with
      ⟨heap2, hcopy, hlayout2, hcontents⟩
    simp only [hcopy]
    refine ⟨heap2, rfl, fun ptr length =>
      (hlayout1 ptr length).trans (hlayout2 ptr length), ?_⟩
    intro k old hk hold
    cases k with
    | zero =>
      have hnow := RawHeap.readU64_writeU64_same heap heap1 dst 0 value hwrite
      have hpres := copyU64_preserves_read heap1 heap2 (dst.add 1)
        (src.add 1) dst n 0 value hDstTail1 hSrcTail1 hnow (by
          intro j hj
          right
          dsimp [RawPtr.add]
          change dst.limbOffset + 1 + j ≠ dst.limbOffset + 0
          omega) hcopy
      have holdEq : old = value := Except.ok.inj (hold.symm.trans hread)
      simpa [holdEq] using hpres
    | succ k =>
      have hk' : k < n := by omega
      have hold1 := RawHeap.readU64_writeU64_ne heap heap1 dst src
        0 (k + 1) value old hwrite hold
        (hdisjoint 0 (by omega) (k + 1) (by omega))
      have htail := hcontents k old hk'
      rw [RawHeap.readU64_add, RawHeap.readU64_add] at htail
      have hout := htail (by simpa [Nat.add_comm] using hold1)
      simpa [Nat.add_comm] using hout
termination_by count

theorem copyU64_slicePolyRep_disjoint (heap : RawHeap)
    (dst src : RawPtr UInt64) (count p : Nat)
    (poly : Polynomial (ZMod p))
    (hDst : heap.ValidU64Slice dst count)
    (hSrc : heap.ValidU64Slice src count)
    (hdisjoint : U64SlicesDisjoint dst count src count)
    (hrep : SlicePolyRep heap src count p poly) :
    ∃ heap', heap.copyU64 dst src count = .ok heap' ∧
      RawHeap.SameLayout heap heap' ∧
      SlicePolyRep heap' dst count p poly := by
  rcases hrep with ⟨coeffs, hslice, hsize, hpoly⟩
  rcases copyU64_refines_disjoint heap dst src count hDst hSrc hdisjoint with
    ⟨heap', hcopy, hlayout, hcontents⟩
  have hDst' := (hlayout dst count).mp hDst
  rcases readU64s_ok heap' dst count hDst' with
    ⟨other, hother, hotherSize⟩
  have heq : other = coeffs := by
    apply Array.ext (hotherSize.trans hsize.symm)
    intro i hiOther hiCoeffs
    have hi : i < count := by simpa [hsize] using hiCoeffs
    have hsrc := readU64s_get heap src count coeffs hslice hsize i hi
    have hdst := hcontents i coeffs[i] hi hsrc
    have hdstOther := readU64s_get heap' dst count other hother
      hotherSize i hi
    exact Except.ok.inj (hdstOther.symm.trans hdst)
  subst other
  exact ⟨heap', hcopy, hlayout, coeffs, hother, hsize, hpoly⟩

theorem copyU64_refines_slice_canonical (heap : RawHeap)
    (dst src : RawPtr UInt64) (count p : Nat)
    (poly : Polynomial (ZMod p)) (modulus : UInt64)
    (hDst : heap.ValidU64Slice dst count)
    (hSrc : heap.ValidU64Slice src count)
    (hdisjoint : U64SlicesDisjoint dst count src count)
    (hrep : SlicePolyRep heap src count p poly)
    (hcanonical : CanonicalU64Prefix heap src count modulus) :
    ∃ heap', heap.copyU64 dst src count = .ok heap' ∧
      RawHeap.SameLayout heap heap' ∧
      SlicePolyRep heap' dst count p poly ∧
      CanonicalU64Prefix heap' dst count modulus := by
  rcases copyU64_refines_disjoint heap dst src count hDst hSrc hdisjoint with
    ⟨heap', hcopy, hlayout, hcontents⟩
  rcases copyU64_slicePolyRep_disjoint heap dst src count p poly hDst hSrc
      hdisjoint
      hrep with ⟨repHeap, hcopyRep, _, hrep'⟩
  have heq : repHeap = heap' := Except.ok.inj (hcopyRep.symm.trans hcopy)
  subst repHeap
  refine ⟨heap', hcopy, hlayout, hrep', ?_⟩
  intro i value hi hread
  rcases heap.readU64_of_valid src count i hSrc hi with ⟨source, hsource⟩
  have hcopied := hcontents i source hi hsource
  have hvalue : value = source := Except.ok.inj (hread.symm.trans hcopied)
  subst value
  exact hcanonical i source hi hsource

theorem slicePolyRep_prefix_exists (heap : RawHeap) (ptr : RawPtr UInt64)
    (length prefixLength p : Nat) (poly : Polynomial (ZMod p))
    (hvalid : heap.ValidU64Slice ptr length)
    (hle : prefixLength ≤ length)
    (hrep : SlicePolyRep heap ptr length p poly) :
    ∃ prefixPoly : Polynomial (ZMod p),
      SlicePolyRep heap ptr prefixLength p prefixPoly ∧
      ∀ i, i < prefixLength → prefixPoly.coeff i = poly.coeff i := by
  have hvalidPrefix := heap.validU64Slice_mono ptr length prefixLength hvalid hle
  rcases slicePolyRep_exists_unique heap ptr prefixLength p hvalidPrefix with
    ⟨prefixPoly, hprefix, _⟩
  refine ⟨prefixPoly, hprefix, ?_⟩
  intro i hi
  rcases slicePolyRep_coeff heap ptr prefixLength p prefixPoly hprefix i hi with
    ⟨prefixValue, hreadPrefix, hcoeffPrefix⟩
  rcases slicePolyRep_coeff heap ptr length p poly hrep i (by omega) with
    ⟨fullValue, hreadFull, hcoeffFull⟩
  have heq : prefixValue = fullValue :=
    Except.ok.inj (hreadPrefix.symm.trans hreadFull)
  rw [hcoeffPrefix, hcoeffFull, heq]

theorem slicePolyRep_split_exists (heap : RawHeap) (ptr : RawPtr UInt64)
    (lowLength highLength p : Nat) (poly : Polynomial (ZMod p))
    (hvalid : heap.ValidU64Slice ptr (lowLength + highLength))
    (hrep : SlicePolyRep heap ptr (lowLength + highLength) p poly) :
    ∃ low high : Polynomial (ZMod p),
      SlicePolyRep heap ptr lowLength p low ∧
      SlicePolyRep heap (ptr.add lowLength) highLength p high ∧
      poly = low + Polynomial.X ^ lowLength * high := by
  have hvalidLow := heap.validU64Slice_mono ptr (lowLength + highLength)
    lowLength hvalid (by omega)
  have hvalidHigh := heap.validU64Slice_add ptr (lowLength + highLength)
    lowLength highLength hvalid (by omega)
  rcases slicePolyRep_exists_unique heap ptr lowLength p hvalidLow with
    ⟨low, hrepLow, _⟩
  rcases slicePolyRep_exists_unique heap (ptr.add lowLength) highLength p
      hvalidHigh with ⟨high, hrepHigh, _⟩
  refine ⟨low, high, hrepLow, hrepHigh, ?_⟩
  ext degree
  by_cases hdlow : degree < lowLength
  · rcases slicePolyRep_coeff heap ptr (lowLength + highLength) p poly hrep
      degree (by omega) with ⟨fullValue, hfull, hcoeffFull⟩
    rcases slicePolyRep_coeff heap ptr lowLength p low hrepLow degree hdlow with
      ⟨lowValue, hlow, hcoeffLow⟩
    have hvalue : fullValue = lowValue := Except.ok.inj (hfull.symm.trans hlow)
    rw [Polynomial.coeff_add, Polynomial.coeff_X_pow_mul', if_neg (by omega),
      hcoeffFull, hcoeffLow, hvalue, add_zero]
  · by_cases hdfull : degree < lowLength + highLength
    · let highDegree := degree - lowLength
      have hdHigh : highDegree < highLength := by
        dsimp [highDegree]
        omega
      have hdegree : lowLength + highDegree = degree := by
        dsimp [highDegree]
        omega
      rcases slicePolyRep_coeff heap ptr (lowLength + highLength) p poly hrep
        degree hdfull with ⟨fullValue, hfull, hcoeffFull⟩
      rcases slicePolyRep_coeff heap (ptr.add lowLength) highLength p high
        hrepHigh highDegree hdHigh with ⟨highValue, hhigh, hcoeffHigh⟩
      have hhighBase : heap.readU64 ptr degree = .ok highValue := by
        rw [← hdegree, ← RawHeap.readU64_add]
        exact hhigh
      have hvalue : fullValue = highValue :=
        Except.ok.inj (hfull.symm.trans hhighBase)
      rw [Polynomial.coeff_add,
        slicePolyRep_coeff_zero_of_length_le heap ptr lowLength p low hrepLow
          degree (by omega),
        Polynomial.coeff_X_pow_mul', if_pos (by omega), hcoeffFull,
        hcoeffHigh, hvalue, zero_add]
    · rw [slicePolyRep_coeff_zero_of_length_le heap ptr
          (lowLength + highLength) p poly hrep degree (by omega),
        Polynomial.coeff_add,
        slicePolyRep_coeff_zero_of_length_le heap ptr lowLength p low hrepLow
          degree (by omega),
        Polynomial.coeff_X_pow_mul', if_pos (by omega),
        slicePolyRep_coeff_zero_of_length_le heap (ptr.add lowLength)
          highLength p high hrepHigh (degree - lowLength) (by omega),
        zero_add]

theorem slicePolyRep_join (heap : RawHeap) (ptr : RawPtr UInt64)
    (lowLength highLength p : Nat)
    (low high : Polynomial (ZMod p))
    (hvalid : heap.ValidU64Slice ptr (lowLength + highLength))
    (hrepLow : SlicePolyRep heap ptr lowLength p low)
    (hrepHigh : SlicePolyRep heap (ptr.add lowLength) highLength p high) :
    SlicePolyRep heap ptr (lowLength + highLength) p
      (low + Polynomial.X ^ lowLength * high) := by
  rcases slicePolyRep_exists_unique heap ptr (lowLength + highLength) p hvalid with
    ⟨whole, hrepWhole, huniqueWhole⟩
  rcases slicePolyRep_split_exists heap ptr lowLength highLength p whole hvalid
      hrepWhole with ⟨actualLow, actualHigh, hrepActualLow,
        hrepActualHigh, hsplit⟩
  have hvalidLow := heap.validU64Slice_mono ptr (lowLength + highLength)
    lowLength hvalid (by omega)
  have hvalidHigh := heap.validU64Slice_add ptr (lowLength + highLength)
    lowLength highLength hvalid (by omega)
  rcases slicePolyRep_exists_unique heap ptr lowLength p hvalidLow with
    ⟨canonicalLow, _, huniqueLow⟩
  rcases slicePolyRep_exists_unique heap (ptr.add lowLength) highLength p
      hvalidHigh with ⟨canonicalHigh, _, huniqueHigh⟩
  have hlow : low = actualLow :=
    (huniqueLow low hrepLow).trans (huniqueLow actualLow hrepActualLow).symm
  have hhigh : high = actualHigh :=
    (huniqueHigh high hrepHigh).trans
      (huniqueHigh actualHigh hrepActualHigh).symm
  have hwhole : whole = low + Polynomial.X ^ lowLength * high := by
    simpa [hlow, hhigh] using hsplit
  rw [← hwhole]
  exact hrepWhole

theorem canonicalU64Prefix_join (heap : RawHeap) (ptr : RawPtr UInt64)
    (lowLength highLength : Nat) (modulus : UInt64)
    (hvalid : heap.ValidU64Slice ptr (lowLength + highLength))
    (hcanonicalLow : CanonicalU64Prefix heap ptr lowLength modulus)
    (hcanonicalHigh : CanonicalU64Prefix heap (ptr.add lowLength)
      highLength modulus) :
    CanonicalU64Prefix heap ptr (lowLength + highLength) modulus := by
  intro i value hi hread
  by_cases hilow : i < lowLength
  · exact hcanonicalLow i value hilow hread
  · let j := i - lowLength
    have hj : j < highLength := by
      dsimp [j]
      omega
    have hij : lowLength + j = i := by
      dsimp [j]
      omega
    apply hcanonicalHigh j value hj
    rw [RawHeap.readU64_add, hij]
    exact hread

theorem writeZero_extends_slice (heap heap' : RawHeap)
    (ptr : RawPtr UInt64) (length p : Nat) (modulus : UInt64)
    (poly : Polynomial (ZMod p))
    (hmodulus : modulus ≠ 0)
    (hvalid : heap.ValidU64Slice ptr (length + 1))
    (hrep : SlicePolyRep heap ptr length p poly)
    (hcanonical : CanonicalU64Prefix heap ptr length modulus)
    (hwrite : heap.writeU64 ptr length 0 = .ok heap') :
    RawHeap.SameLayout heap heap' ∧
      SlicePolyRep heap' ptr (length + 1) p poly ∧
      CanonicalU64Prefix heap' ptr (length + 1) modulus := by
  have hlayout := RawHeap.writeU64_sameLayout heap heap' ptr length 0 hwrite
  have hvalidPrefix := heap.validU64Slice_mono ptr (length + 1) length hvalid
    (by omega)
  have hvalid' := (hlayout ptr (length + 1)).mp hvalid
  have hvalidPrefix' := (hlayout ptr length).mp hvalidPrefix
  have hsame : SameU64Prefix heap heap' ptr length := by
    intro i old hi hread
    exact RawHeap.readU64_writeU64_ne heap heap' ptr ptr length i 0 old
      hwrite hread (Or.inr (by omega))
  have hrepPrefix' := slicePolyRep_of_same_prefix heap heap' ptr length p poly
    hvalidPrefix hvalidPrefix' hsame hrep
  have hzero := RawHeap.readU64_writeU64_same heap heap' ptr length 0 hwrite
  rcases slicePolyRep_extend_exists heap' ptr length p 0 poly hvalid'
      hrepPrefix' hzero with ⟨full, hrepFull, hfull⟩
  have hfullEq : full = poly := by simpa using hfull
  rw [hfullEq] at hrepFull
  refine ⟨hlayout, hrepFull, ?_⟩
  intro i value hi hread
  by_cases hilength : i < length
  · rcases heap.readU64_of_valid ptr length i hvalidPrefix hilength with
      ⟨old, hold⟩
    have hpreserved := hsame i old hilength hold
    have heq : value = old := Except.ok.inj (hread.symm.trans hpreserved)
    subst value
    exact hcanonical i old hilength hold
  · have hiEq : i = length := by omega
    subst i
    have heq : value = 0 := Except.ok.inj (hread.symm.trans hzero)
    subst value
    have hpositive : 0 < modulus.toNat := by
      apply Nat.pos_of_ne_zero
      intro hzeroNat
      apply hmodulus
      exact UInt64.toNat_inj.mp (by simpa using hzeroNat)
    simpa using hpositive

theorem karCopyZero_refines_base (heap heap8 heap9 : RawHeap)
    (C P0 : RawPtr UInt64) (m highLength p : Nat)
    (lowPoly highPoly : Polynomial (ZMod p)) (modulus : UInt64)
    (hm : 0 < m) (hmodulus : modulus ≠ 0)
    (hC : heap.ValidU64Slice C (2 * m + highLength))
    (hP0 : heap.ValidU64Slice P0 (2 * m - 1))
    (hCP0 : U64SlicesDisjoint C (2 * m - 1) P0 (2 * m - 1))
    (hRepP0 : SlicePolyRep heap P0 (2 * m - 1) p lowPoly)
    (hCanonicalP0 : CanonicalU64Prefix heap P0 (2 * m - 1) modulus)
    (hRepHigh : SlicePolyRep heap (C.add (2 * m)) highLength p highPoly)
    (hCanonicalHigh : CanonicalU64Prefix heap (C.add (2 * m))
      highLength modulus)
    (hcopy : heap.copyU64 C P0 (2 * m - 1) = .ok heap8)
    (hzero : heap8.writeU64 C (2 * m - 1) 0 = .ok heap9) :
    RawHeap.SameLayout heap heap9 ∧
      SlicePolyRep heap9 C (2 * m + highLength) p
        (lowPoly + Polynomial.X ^ (2 * m) * highPoly) ∧
      CanonicalU64Prefix heap9 C (2 * m + highLength) modulus := by
  have hCPrefix := heap.validU64Slice_mono C (2 * m + highLength)
    (2 * m - 1) hC (by omega)
  have hCThroughZero := heap.validU64Slice_mono C (2 * m + highLength)
    ((2 * m - 1) + 1) hC (by omega)
  have hCHigh := heap.validU64Slice_add C (2 * m + highLength)
    (2 * m) highLength hC (by omega)
  rcases copyU64_refines_slice_canonical heap C P0 (2 * m - 1) p lowPoly
      modulus hCPrefix hP0 hCP0 hRepP0 hCanonicalP0 with
    ⟨copyHeap, hcopy', hlayout8, hRepLow8, hCanonicalLow8⟩
  have hcopyHeap : copyHeap = heap8 := Except.ok.inj (hcopy'.symm.trans hcopy)
  subst copyHeap
  have hLowHighDisjoint : U64SlicesDisjoint C (2 * m - 1)
      (C.add (2 * m)) highLength := by
    have h := u64SlicesDisjoint_add_of_le C 0 (2 * m - 1) (2 * m)
      highLength (by omega)
    simpa [RawPtr.add] using h
  have hsameHigh8 := copyU64_preserves_prefix heap heap8 C P0
    (C.add (2 * m)) (2 * m - 1) highLength hCPrefix hP0
    hLowHighDisjoint hcopy
  have hCHigh8 := (hlayout8 (C.add (2 * m)) highLength).mp hCHigh
  have hRepHigh8 := slicePolyRep_of_same_prefix heap heap8 (C.add (2 * m))
    highLength p highPoly hCHigh hCHigh8 hsameHigh8 hRepHigh
  have hCanonicalHigh8 : CanonicalU64Prefix heap8 (C.add (2 * m))
      highLength modulus := by
    intro k value hk hread8
    rcases heap.readU64_of_valid (C.add (2 * m)) highLength k hCHigh hk with
      ⟨old, hread⟩
    have hpreserved := hsameHigh8 k old hk hread
    have hvalue : value = old := Except.ok.inj (hread8.symm.trans hpreserved)
    subst value
    exact hCanonicalHigh k old hk hread
  have hCThroughZero8 := (hlayout8 C ((2 * m - 1) + 1)).mp hCThroughZero
  rcases writeZero_extends_slice heap8 heap9 C (2 * m - 1) p modulus
      lowPoly hmodulus hCThroughZero8 hRepLow8 hCanonicalLow8 hzero with
    ⟨hlayout9, hRepLow9, hCanonicalLow9⟩
  have hPrefixHighDisjoint : U64SlicesDisjoint C ((2 * m - 1) + 1)
      (C.add (2 * m)) highLength := by
    have h := u64SlicesDisjoint_add_of_le C 0 ((2 * m - 1) + 1)
      (2 * m) highLength (by omega)
    simpa [RawPtr.add] using h
  have hsameHigh9 := writeU64_preserves_prefix heap8 heap9 C
    (C.add (2 * m)) ((2 * m - 1) + 1) highLength (2 * m - 1) 0
    hPrefixHighDisjoint (by omega) hzero
  have hCHigh9 := (hlayout9 (C.add (2 * m)) highLength).mp hCHigh8
  have hRepHigh9 := slicePolyRep_of_same_prefix heap8 heap9 (C.add (2 * m))
    highLength p highPoly hCHigh8 hCHigh9 hsameHigh9 hRepHigh8
  have hCanonicalHigh9 : CanonicalU64Prefix heap9 (C.add (2 * m))
      highLength modulus := by
    intro k value hk hread9
    rcases heap8.readU64_of_valid (C.add (2 * m)) highLength k hCHigh8 hk with
      ⟨old, hread⟩
    have hpreserved := hsameHigh9 k old hk hread
    have hvalue : value = old := Except.ok.inj (hread9.symm.trans hpreserved)
    subst value
    exact hCanonicalHigh8 k old hk hread
  have hC9 := (hlayout9 C (2 * m + highLength)).mp
    ((hlayout8 C (2 * m + highLength)).mp hC)
  have hRepLow9' : SlicePolyRep heap9 C (2 * m) p lowPoly := by
    simpa [Nat.sub_add_cancel (by omega : 1 ≤ 2 * m)] using hRepLow9
  have hCanonicalLow9' : CanonicalU64Prefix heap9 C (2 * m) modulus := by
    simpa [Nat.sub_add_cancel (by omega : 1 ≤ 2 * m)] using hCanonicalLow9
  refine ⟨fun ptr length => (hlayout8 ptr length).trans (hlayout9 ptr length),
    slicePolyRep_join heap9 C (2 * m) highLength p lowPoly highPoly hC9
      hRepLow9' hRepHigh9,
    canonicalU64Prefix_join heap9 C (2 * m) highLength modulus hC9
      hCanonicalLow9' hCanonicalHigh9⟩

/-- Exact number of UInt64 scratch cells reachable by the generated
Karatsuba recursion.  The three child products share `recScratch`, so the
recursive contribution is a maximum rather than a sum. -/
def karScratchNeed (n : Nat) : Nat :=
  if n < 16 then
    0
  else
    let m := n / 2
    let h := n - m
    2 * h + (2 * m - 1) + (2 * h - 1) +
      max (karScratchNeed m) (karScratchNeed h)
termination_by n
decreasing_by
  all_goals
    have hn : 0 < n := by omega
    have hm : 0 < n / 2 := Nat.div_pos (by omega) (by omega)
    omega

theorem karScratchNeed_base (n : Nat) (hn : n < 16) :
    karScratchNeed n = 0 := by
  rw [karScratchNeed, if_pos hn]

theorem karScratchNeed_step (n : Nat) (hn : ¬n < 16) :
    karScratchNeed n =
      let m := n / 2
      let h := n - m
      2 * h + (2 * m - 1) + (2 * h - 1) +
        max (karScratchNeed m) (karScratchNeed h) := by
  rw [karScratchNeed, if_neg hn]

theorem kar_split_shape (n : Nat) :
    let m := n / 2
    let h := n - m
    h = m ∨ h = m + 1 := by
  omega

theorem kar_split_children_lt (n : Nat) (hn : 16 ≤ n) :
    n / 2 < n ∧ n - n / 2 < n := by
  have hm : 0 < n / 2 := Nat.div_pos (by omega) (by omega)
  constructor
  · exact Nat.div_lt_self (by omega) (by omega)
  · omega

theorem karScratchNeed_child_le_rec (n : Nat) (hn : 16 ≤ n) :
    let m := n / 2
    let h := n - m
    karScratchNeed m ≤ max (karScratchNeed m) (karScratchNeed h) ∧
      karScratchNeed h ≤ max (karScratchNeed m) (karScratchNeed h) := by
  exact ⟨Nat.le_max_left _ _, Nat.le_max_right _ _⟩

theorem karScratchNeed_current_le (n : Nat) (hn : 16 ≤ n) :
    let m := n / 2
    let h := n - m
    2 * h + (2 * m - 1) + (2 * h - 1) ≤ karScratchNeed n := by
  rw [karScratchNeed_step n (by omega)]
  dsimp
  omega

/-- The single scratch allocation accepted by `_kar_mul` contains every
concrete segment carved out by the generated C++ pointer arithmetic. -/
theorem karScratchSlices (heap : RawHeap) (scratch : RawPtr UInt64)
    (n : Nat) (hn : 16 ≤ n)
    (hScratch : heap.ValidU64Slice scratch (karScratchNeed n)) :
    let m := n / 2
    let h := n - m
    let t1 := scratch
    let t2 := t1.add h
    let sP0 := t2.add h
    let sP1 := sP0.add (2 * m - 1)
    let recScratch := sP1.add (2 * h - 1)
    heap.ValidU64Slice t1 h ∧
      heap.ValidU64Slice t2 h ∧
      heap.ValidU64Slice sP0 (2 * m - 1) ∧
      heap.ValidU64Slice sP1 (2 * h - 1) ∧
      heap.ValidU64Slice recScratch
        (max (karScratchNeed m) (karScratchNeed h)) := by
  let m := n / 2
  let h := n - m
  have hneed : karScratchNeed n =
      2 * h + (2 * m - 1) + (2 * h - 1) +
        max (karScratchNeed m) (karScratchNeed h) := by
    simpa [m, h] using karScratchNeed_step n (by omega)
  have ht1 := heap.validU64Slice_mono scratch (karScratchNeed n) h
    hScratch (by rw [hneed]; omega)
  have ht2 := heap.validU64Slice_add scratch (karScratchNeed n) h h
    hScratch (by rw [hneed]; omega)
  have hsP0 := heap.validU64Slice_add scratch (karScratchNeed n) (2 * h)
    (2 * m - 1) hScratch (by rw [hneed]; omega)
  have hsP1 := heap.validU64Slice_add scratch (karScratchNeed n)
    (2 * h + (2 * m - 1)) (2 * h - 1) hScratch
    (by rw [hneed]; omega)
  have hrec := heap.validU64Slice_add scratch (karScratchNeed n)
    (2 * h + (2 * m - 1) + (2 * h - 1))
    (max (karScratchNeed m) (karScratchNeed h)) hScratch
    (by rw [hneed])
  dsimp [m, h] at ht1 ht2 hsP0 hsP1 hrec ⊢
  refine ⟨ht1, ht2, ?_, ?_, ?_⟩
  · have hp : (scratch.add (n - n / 2)).add (n - n / 2) =
        scratch.add (2 * (n - n / 2)) := by
      cases scratch
      have hwidth : RawLimbWidth.width UInt64 = 1 := rfl
      simp [RawPtr.add, hwidth]
      omega
    rw [hp]
    exact hsP0
  · have hp : ((scratch.add (n - n / 2)).add (n - n / 2)).add
          (2 * (n / 2) - 1) =
        scratch.add (2 * (n - n / 2) + (2 * (n / 2) - 1)) := by
      cases scratch
      have hwidth : RawLimbWidth.width UInt64 = 1 := rfl
      simp [RawPtr.add, hwidth]
      omega
    rw [hp]
    exact hsP1
  · have hp : (((scratch.add (n - n / 2)).add (n - n / 2)).add
          (2 * (n / 2) - 1)).add (2 * (n - n / 2) - 1) =
        scratch.add (2 * (n - n / 2) + (2 * (n / 2) - 1) +
          (2 * (n - n / 2) - 1)) := by
      cases scratch
      have hwidth : RawLimbWidth.width UInt64 = 1 := rfl
      simp [RawPtr.add, hwidth]
      omega
    rw [hp]
    exact hrec

theorem karScratchSlices_disjoint_guard (scratch guard : RawPtr UInt64)
    (n guardLength : Nat) (hn : 16 ≤ n)
    (hdisjoint : U64SlicesDisjoint scratch (karScratchNeed n)
      guard guardLength) :
    let m := n / 2
    let h := n - m
    let t1 := scratch
    let t2 := t1.add h
    let sP0 := t2.add h
    let sP1 := sP0.add (2 * m - 1)
    let recScratch := sP1.add (2 * h - 1)
    U64SlicesDisjoint t1 h guard guardLength ∧
      U64SlicesDisjoint t2 h guard guardLength ∧
      U64SlicesDisjoint sP0 (2 * m - 1) guard guardLength ∧
      U64SlicesDisjoint sP1 (2 * h - 1) guard guardLength ∧
      U64SlicesDisjoint recScratch
        (max (karScratchNeed m) (karScratchNeed h)) guard guardLength := by
  let m := n / 2
  let h := n - m
  have hneed : karScratchNeed n =
      2 * h + (2 * m - 1) + (2 * h - 1) +
        max (karScratchNeed m) (karScratchNeed h) := by
    simpa [m, h] using karScratchNeed_step n (by omega)
  have ht1 := u64SlicesDisjoint_mono hdisjoint (smallLeft := h)
    (by rw [hneed]; omega) (Nat.le_refl guardLength)
  have ht2 := u64SlicesDisjoint_add_left hdisjoint (start := h) (count := h)
    (by rw [hneed]; omega)
  have hp0 := u64SlicesDisjoint_add_left hdisjoint (start := 2 * h)
    (count := 2 * m - 1) (by rw [hneed]; omega)
  have hp1 := u64SlicesDisjoint_add_left hdisjoint
    (start := 2 * h + (2 * m - 1)) (count := 2 * h - 1)
    (by rw [hneed]; omega)
  have hrec := u64SlicesDisjoint_add_left hdisjoint
    (start := 2 * h + (2 * m - 1) + (2 * h - 1))
    (count := max (karScratchNeed m) (karScratchNeed h))
    (by rw [hneed])
  dsimp [m, h] at ht1 ht2 hp0 hp1 hrec ⊢
  refine ⟨ht1, ht2, ?_, ?_, ?_⟩
  · rw [rawPtr_add_add, show (n - n / 2) + (n - n / 2) =
        2 * (n - n / 2) by omega]
    exact hp0
  · rw [rawPtr_add_add, rawPtr_add_add,
      show (n - n / 2) + ((n - n / 2) + (2 * (n / 2) - 1)) =
        2 * (n - n / 2) + (2 * (n / 2) - 1) by omega]
    exact hp1
  · rw [rawPtr_add_add, rawPtr_add_add, rawPtr_add_add,
      show (n - n / 2) + ((n - n / 2) +
          ((2 * (n / 2) - 1) + (2 * (n - n / 2) - 1))) =
        2 * (n - n / 2) + (2 * (n / 2) - 1) +
          (2 * (n - n / 2) - 1) by omega]
    exact hrec

theorem karScratchNeed_le_seven (n : Nat) : karScratchNeed n ≤ 7 * n := by
  rw [karScratchNeed]
  split
  next hbase => omega
  next hrec =>
    let m := n / 2
    let h := n - m
    have hn16 : 16 ≤ n := by omega
    have hmLt : m < n := (kar_split_children_lt n hn16).1
    have hhLt : h < n := by
      dsimp [h, m]
      exact (kar_split_children_lt n hn16).2
    have ihm := karScratchNeed_le_seven m
    have ihh := karScratchNeed_le_seven h
    have hmLeH : m ≤ h := by
      dsimp [m, h]
      omega
    have hmax : max (karScratchNeed m) (karScratchNeed h) ≤ 7 * h := by
      apply max_le
      · exact le_trans ihm (by omega)
      · exact ihh
    have hshape : h = m ∨ h = m + 1 := by
      simpa [m, h] using kar_split_shape n
    rcases hshape with heven | hodd
    · dsimp [m, h] at hmax ⊢
      omega
    · dsimp [m, h] at hmax ⊢
      omega
termination_by n
decreasing_by
  · exact hmLt
  · exact hhLt

theorem karAddHalvesLoop_ok (this : DenseUPolyZp)
    (A B t1 t2 : RawPtr UInt64) (m i : Nat) (heap : RawHeap)
    (hA : heap.ValidU64Slice A (2 * m))
    (hB : heap.ValidU64Slice B (2 * m))
    (hT1 : heap.ValidU64Slice t1 m)
    (hT2 : heap.ValidU64Slice t2 m) :
    ∃ heap', karAddHalvesLoop this A B t1 t2 m i heap = .ok heap' ∧
      RawHeap.SameLayout heap heap' := by
  unfold karAddHalvesLoop
  split
  next hi =>
    rcases heap.readU64_of_valid A (2 * m) i hA (by omega) with ⟨alo, halo⟩
    simp only [halo]
    rcases heap.readU64_of_valid A (2 * m) (m + i) hA (by omega) with
      ⟨ahi, hahi⟩
    simp only [hahi]
    rcases heap.readU64_of_valid B (2 * m) i hB (by omega) with ⟨blo, hblo⟩
    simp only [hblo]
    rcases heap.readU64_of_valid B (2 * m) (m + i) hB (by omega) with
      ⟨bhi, hbhi⟩
    simp only [hbhi]
    let av := dense_upoly_zp_nmod_add_ir this alo ahi
    rcases heap.writeU64_of_valid t1 m i av hT1 hi with ⟨heap1, hw1⟩
    simp only [av, hw1]
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 t1 i av hw1
    let bv := dense_upoly_zp_nmod_add_ir this blo bhi
    rcases heap1.writeU64_of_valid t2 m i bv ((hlayout1 t2 m).mp hT2) hi with
      ⟨heap2, hw2⟩
    simp only [bv, hw2]
    have hlayout2 := RawHeap.writeU64_sameLayout heap1 heap2 t2 i bv hw2
    rcases karAddHalvesLoop_ok this A B t1 t2 m (i + 1) heap2
      ((hlayout2 A (2 * m)).mp ((hlayout1 A (2 * m)).mp hA))
      ((hlayout2 B (2 * m)).mp ((hlayout1 B (2 * m)).mp hB))
      ((hlayout2 t1 m).mp ((hlayout1 t1 m).mp hT1))
      ((hlayout2 t2 m).mp ((hlayout1 t2 m).mp hT2)) with
      ⟨heap3, hrun, hlayout3⟩
    rw [hrun]
    exact ⟨heap3, rfl, fun ptr length =>
      (hlayout1 ptr length).trans
        ((hlayout2 ptr length).trans (hlayout3 ptr length))⟩
  next => exact ⟨heap, rfl, fun _ _ => Iff.rfl⟩
termination_by m - i
decreasing_by omega

theorem karSubLoop_ok (this : DenseUPolyZp)
    (dst sub : RawPtr UInt64) (count i : Nat) (heap : RawHeap)
    (hDst : heap.ValidU64Slice dst count)
    (hSub : heap.ValidU64Slice sub count) :
    ∃ heap', karSubLoop this dst sub count i heap = .ok heap' ∧
      RawHeap.SameLayout heap heap' := by
  unfold karSubLoop
  split
  next hi =>
    rcases heap.readU64_of_valid dst count i hDst hi with ⟨a, ha⟩
    simp only [ha]
    rcases heap.readU64_of_valid sub count i hSub hi with ⟨b, hb⟩
    simp only [hb]
    let value := dense_upoly_zp_nmod_sub_ir this a b
    rcases heap.writeU64_of_valid dst count i value hDst hi with ⟨heap1, hw⟩
    simp only [value, hw]
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 dst i value hw
    rcases karSubLoop_ok this dst sub count (i + 1) heap1
      ((hlayout1 dst count).mp hDst) ((hlayout1 sub count).mp hSub) with
      ⟨heap2, hrun, hlayout2⟩
    rw [hrun]
    exact ⟨heap2, rfl, fun ptr length =>
      (hlayout1 ptr length).trans (hlayout2 ptr length)⟩
  next => exact ⟨heap, rfl, fun _ _ => Iff.rfl⟩
termination_by count - i
decreasing_by omega

theorem karAssembleLoop_ok (this : DenseUPolyZp)
    (C sP1 : RawPtr UInt64) (m count i : Nat) (heap : RawHeap)
    (hC : heap.ValidU64Slice C (m + count))
    (hP1 : heap.ValidU64Slice sP1 count) :
    ∃ heap', karAssembleLoop this C sP1 m count i heap = .ok heap' ∧
      RawHeap.SameLayout heap heap' := by
  unfold karAssembleLoop
  split
  next hi =>
    rcases heap.readU64_of_valid C (m + count) (m + i) hC (by omega) with
      ⟨base, hbase⟩
    simp only [hbase]
    rcases heap.readU64_of_valid sP1 count i hP1 hi with ⟨cross, hcross⟩
    simp only [hcross]
    let value := dense_upoly_zp_nmod_add_ir this base cross
    rcases heap.writeU64_of_valid C (m + count) (m + i) value hC (by omega) with
      ⟨heap1, hw⟩
    simp only [value, hw]
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 C (m + i) value hw
    rcases karAssembleLoop_ok this C sP1 m count (i + 1) heap1
      ((hlayout1 C (m + count)).mp hC) ((hlayout1 sP1 count).mp hP1) with
      ⟨heap2, hrun, hlayout2⟩
    rw [hrun]
    exact ⟨heap2, rfl, fun ptr length =>
      (hlayout1 ptr length).trans (hlayout2 ptr length)⟩
  next => exact ⟨heap, rfl, fun _ _ => Iff.rfl⟩
termination_by count - i
decreasing_by omega

theorem karOddTail_ok (A B t1 t2 : RawPtr UInt64) (m h : Nat)
    (heap : RawHeap)
    (hA : heap.ValidU64Slice A (m + h))
    (hB : heap.ValidU64Slice B (m + h))
    (hT1 : heap.ValidU64Slice t1 h)
    (hT2 : heap.ValidU64Slice t2 h) :
    ∃ heap', karOddTail A B t1 t2 m h heap = .ok heap' ∧
      RawHeap.SameLayout heap heap' := by
  unfold karOddTail
  split
  next hodd =>
    rcases heap.readU64_of_valid A (m + h) (m + m) hA (by omega) with
      ⟨aTail, ha⟩
    simp only [ha]
    rcases heap.readU64_of_valid B (m + h) (m + m) hB (by omega) with
      ⟨bTail, hb⟩
    simp only [hb]
    rcases heap.writeU64_of_valid t1 h m aTail hT1 hodd with ⟨heap1, hw1⟩
    simp only [hw1]
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 t1 m aTail hw1
    rcases heap1.writeU64_of_valid t2 h m bTail ((hlayout1 t2 h).mp hT2)
      hodd with ⟨heap2, hw2⟩
    simp only [hw2]
    have hlayout2 := RawHeap.writeU64_sameLayout heap1 heap2 t2 m bTail hw2
    exact ⟨heap2, rfl, fun ptr length =>
      (hlayout1 ptr length).trans (hlayout2 ptr length)⟩
  next => exact ⟨heap, rfl, fun _ _ => Iff.rfl⟩

theorem karOddTail_preserves_outside (A B t1 t2 guard : RawPtr UInt64)
    (m h guardLen : Nat) (heap heap' : RawHeap)
    (hA : heap.ValidU64Slice A (m + h))
    (hB : heap.ValidU64Slice B (m + h))
    (hT1 : heap.ValidU64Slice t1 h)
    (hT2 : heap.ValidU64Slice t2 h)
    (hT1Guard : U64SlicesDisjoint t1 h guard guardLen)
    (hT2Guard : U64SlicesDisjoint t2 h guard guardLen)
    (hrun : karOddTail A B t1 t2 m h heap = .ok heap') :
    SameU64Prefix heap heap' guard guardLen := by
  intro i old hi hread
  unfold karOddTail at hrun
  split at hrun
  next hodd =>
    rcases heap.readU64_of_valid A (m + h) (m + m) hA (by omega) with
      ⟨aTail, ha⟩
    simp only [ha] at hrun
    rcases heap.readU64_of_valid B (m + h) (m + m) hB (by omega) with
      ⟨bTail, hb⟩
    simp only [hb] at hrun
    rcases heap.writeU64_of_valid t1 h m aTail hT1 hodd with ⟨heap1, hw1⟩
    simp only [hw1] at hrun
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 t1 m aTail hw1
    rcases heap1.writeU64_of_valid t2 h m bTail ((hlayout1 t2 h).mp hT2)
      hodd with ⟨heap2, hw2⟩
    simp only [hw2] at hrun
    have heq : heap' = heap2 := Except.ok.inj hrun.symm
    subst heap'
    have hread1 := RawHeap.readU64_writeU64_ne heap heap1 t1 guard m i
      aTail old hw1 hread (hT1Guard m hodd i hi)
    exact RawHeap.readU64_writeU64_ne heap1 heap2 t2 guard m i bTail old
      hw2 hread1 (hT2Guard m hodd i hi)
  next hnot =>
    have heq : heap' = heap := Except.ok.inj hrun.symm
    simpa [heq] using hread

theorem karPrepareHalves_ok (this : DenseUPolyZp)
    (A B t1 t2 : RawPtr UInt64) (m h : Nat) (heap : RawHeap)
    (hmh : m ≤ h)
    (hA : heap.ValidU64Slice A (m + h))
    (hB : heap.ValidU64Slice B (m + h))
    (hT1 : heap.ValidU64Slice t1 h)
    (hT2 : heap.ValidU64Slice t2 h) :
    ∃ heap', karPrepareHalves this A B t1 t2 m h heap = .ok heap' ∧
      RawHeap.SameLayout heap heap' := by
  have hA2m := heap.validU64Slice_mono A (m + h) (2 * m) hA (by omega)
  have hB2m := heap.validU64Slice_mono B (m + h) (2 * m) hB (by omega)
  have hT1m := heap.validU64Slice_mono t1 h m hT1 hmh
  have hT2m := heap.validU64Slice_mono t2 h m hT2 hmh
  rcases karAddHalvesLoop_ok this A B t1 t2 m 0 heap hA2m hB2m hT1m
      hT2m with ⟨heap1, hadd, hlayout1⟩
  have hA1 := (hlayout1 A (m + h)).mp hA
  have hB1 := (hlayout1 B (m + h)).mp hB
  have hT11 := (hlayout1 t1 h).mp hT1
  have hT21 := (hlayout1 t2 h).mp hT2
  rcases karOddTail_ok A B t1 t2 m h heap1 hA1 hB1 hT11 hT21 with
    ⟨heap2, htail, hlayout2⟩
  refine ⟨heap2, ?_, fun ptr length =>
    (hlayout1 ptr length).trans (hlayout2 ptr length)⟩
  simp [karPrepareHalves, hadd, htail]

theorem karOddTail_values (A B t1 t2 : RawPtr UInt64) (m h : Nat)
    (heap heap' : RawHeap)
    (hA : heap.ValidU64Slice A (m + h))
    (hB : heap.ValidU64Slice B (m + h))
    (hT1 : heap.ValidU64Slice t1 h)
    (hT2 : heap.ValidU64Slice t2 h)
    (hT1T2 : U64SlicesDisjoint t1 h t2 h)
    (hodd : h > m)
    (hrun : karOddTail A B t1 t2 m h heap = .ok heap') :
    ∃ aTail bTail,
      heap.readU64 A (m + m) = .ok aTail ∧
      heap.readU64 B (m + m) = .ok bTail ∧
      heap'.readU64 t1 m = .ok aTail ∧
      heap'.readU64 t2 m = .ok bTail := by
  unfold karOddTail at hrun
  simp only [hodd, ↓reduceIte] at hrun
  rcases heap.readU64_of_valid A (m + h) (m + m) hA (by omega) with
    ⟨aTail, ha⟩
  simp only [ha] at hrun
  rcases heap.readU64_of_valid B (m + h) (m + m) hB (by omega) with
    ⟨bTail, hb⟩
  simp only [hb] at hrun
  rcases heap.writeU64_of_valid t1 h m aTail hT1 hodd with ⟨heap1, hw1⟩
  simp only [hw1] at hrun
  have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 t1 m aTail hw1
  rcases heap1.writeU64_of_valid t2 h m bTail ((hlayout1 t2 h).mp hT2)
    hodd with ⟨heap2, hw2⟩
  simp only [hw2] at hrun
  have heq : heap' = heap2 := Except.ok.inj hrun.symm
  subst heap'
  have hread1 := RawHeap.readU64_writeU64_same heap heap1 t1 m aTail hw1
  have hread1' := RawHeap.readU64_writeU64_ne heap1 heap2 t2 t1 m m
    bTail aTail hw2 hread1
      (u64SlicesDisjoint_symm hT1T2 m hodd m hodd)
  have hread2 := RawHeap.readU64_writeU64_same heap1 heap2 t2 m bTail hw2
  exact ⟨aTail, bTail, ha, hb, hread1', hread2⟩

theorem karOddTail_preserves_own_prefixes
    (A B t1 t2 : RawPtr UInt64) (m h : Nat) (heap heap' : RawHeap)
    (hA : heap.ValidU64Slice A (m + h))
    (hB : heap.ValidU64Slice B (m + h))
    (hT1 : heap.ValidU64Slice t1 h)
    (hT2 : heap.ValidU64Slice t2 h)
    (hT1T2 : U64SlicesDisjoint t1 h t2 h)
    (hrun : karOddTail A B t1 t2 m h heap = .ok heap') :
    SameU64Prefix heap heap' t1 m ∧ SameU64Prefix heap heap' t2 m := by
  constructor <;> intro i old hi hread
  · unfold karOddTail at hrun
    split at hrun
    next hodd =>
      rcases heap.readU64_of_valid A (m + h) (m + m) hA (by omega) with
        ⟨aTail, ha⟩
      simp only [ha] at hrun
      rcases heap.readU64_of_valid B (m + h) (m + m) hB (by omega) with
        ⟨bTail, hb⟩
      simp only [hb] at hrun
      rcases heap.writeU64_of_valid t1 h m aTail hT1 hodd with ⟨heap1, hw1⟩
      simp only [hw1] at hrun
      have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 t1 m aTail hw1
      rcases heap1.writeU64_of_valid t2 h m bTail ((hlayout1 t2 h).mp hT2)
        hodd with ⟨heap2, hw2⟩
      simp only [hw2] at hrun
      have heq : heap' = heap2 := Except.ok.inj hrun.symm
      subst heap'
      have hread1 := RawHeap.readU64_writeU64_ne heap heap1 t1 t1 m i
        aTail old hw1 hread (Or.inr (by omega))
      exact RawHeap.readU64_writeU64_ne heap1 heap2 t2 t1 m i bTail old
        hw2 hread1 (u64SlicesDisjoint_symm hT1T2 m hodd i (by omega))
    next hnot =>
      have heq : heap' = heap := Except.ok.inj hrun.symm
      simpa [heq] using hread
  · unfold karOddTail at hrun
    split at hrun
    next hodd =>
      rcases heap.readU64_of_valid A (m + h) (m + m) hA (by omega) with
        ⟨aTail, ha⟩
      simp only [ha] at hrun
      rcases heap.readU64_of_valid B (m + h) (m + m) hB (by omega) with
        ⟨bTail, hb⟩
      simp only [hb] at hrun
      rcases heap.writeU64_of_valid t1 h m aTail hT1 hodd with ⟨heap1, hw1⟩
      simp only [hw1] at hrun
      have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 t1 m aTail hw1
      rcases heap1.writeU64_of_valid t2 h m bTail ((hlayout1 t2 h).mp hT2)
        hodd with ⟨heap2, hw2⟩
      simp only [hw2] at hrun
      have heq : heap' = heap2 := Except.ok.inj hrun.symm
      subst heap'
      have hread1 := RawHeap.readU64_writeU64_ne heap heap1 t1 t2 m i
        aTail old hw1 hread (hT1T2 m hodd i (by omega))
      exact RawHeap.readU64_writeU64_ne heap1 heap2 t2 t2 m i bTail old
        hw2 hread1 (Or.inr (by omega))
    next hnot =>
      have heq : heap' = heap := Except.ok.inj hrun.symm
      simpa [heq] using hread

theorem karOddTail_coeffs (this : DenseUPolyZp)
    (A B t1 t2 : RawPtr UInt64) (m h : Nat) (heap heap' : RawHeap)
    (left right : Polynomial (ZMod this._p.toNat))
    (hA : heap.ValidU64Slice A (m + h))
    (hB : heap.ValidU64Slice B (m + h))
    (hT1 : heap.ValidU64Slice t1 h)
    (hT2 : heap.ValidU64Slice t2 h)
    (hT1T2 : U64SlicesDisjoint t1 h t2 h)
    (hCanonicalA : CanonicalU64Prefix heap A (m + h) this._p)
    (hCanonicalB : CanonicalU64Prefix heap B (m + h) this._p)
    (hRepA : SlicePolyRep heap A (m + h) this._p.toNat left)
    (hRepB : SlicePolyRep heap B (m + h) this._p.toNat right)
    (hodd : h > m)
    (hrun : karOddTail A B t1 t2 m h heap = .ok heap') :
    ∃ aTail bTail,
      heap'.readU64 t1 m = .ok aTail ∧
      heap'.readU64 t2 m = .ok bTail ∧
      (aTail.toNat : ZMod this._p.toNat) = left.coeff (m + m) ∧
      (bTail.toNat : ZMod this._p.toNat) = right.coeff (m + m) ∧
      aTail.toNat < this._p.toNat ∧ bTail.toNat < this._p.toNat := by
  rcases karOddTail_values A B t1 t2 m h heap heap' hA hB hT1 hT2
      hT1T2 hodd hrun with ⟨aTail, bTail, ha, hb, ht1, ht2⟩
  rcases slicePolyRep_coeff heap A (m + h) this._p.toNat left hRepA
      (m + m) (by omega) with ⟨a', ha', hcoeffA⟩
  have haEq : a' = aTail := Except.ok.inj (ha'.symm.trans ha)
  subst a'
  rcases slicePolyRep_coeff heap B (m + h) this._p.toNat right hRepB
      (m + m) (by omega) with ⟨b', hb', hcoeffB⟩
  have hbEq : b' = bTail := Except.ok.inj (hb'.symm.trans hb)
  subst b'
  exact ⟨aTail, bTail, ht1, ht2, hcoeffA.symm, hcoeffB.symm,
    hCanonicalA (m + m) aTail (by omega) ha,
    hCanonicalB (m + m) bTail (by omega) hb⟩

theorem karAddHalvesLoop_preserves_outside (this : DenseUPolyZp)
    (A B t1 t2 guard : RawPtr UInt64) (m i readIndex : Nat)
    (heap heap' : RawHeap) (old : UInt64)
    (hA : heap.ValidU64Slice A (2 * m))
    (hB : heap.ValidU64Slice B (2 * m))
    (hT1 : heap.ValidU64Slice t1 m)
    (hT2 : heap.ValidU64Slice t2 m)
    (hread : heap.readU64 guard readIndex = .ok old)
    (hout1 : ∀ j, i ≤ j → j < m → t1.region ≠ guard.region ∨
      t1.limbOffset + j ≠ guard.limbOffset + readIndex)
    (hout2 : ∀ j, i ≤ j → j < m → t2.region ≠ guard.region ∨
      t2.limbOffset + j ≠ guard.limbOffset + readIndex)
    (hrun : karAddHalvesLoop this A B t1 t2 m i heap = .ok heap') :
    heap'.readU64 guard readIndex = .ok old := by
  unfold karAddHalvesLoop at hrun
  split at hrun
  next hi =>
    rcases heap.readU64_of_valid A (2 * m) i hA (by omega) with ⟨alo, halo⟩
    simp only [halo] at hrun
    rcases heap.readU64_of_valid A (2 * m) (m + i) hA (by omega) with
      ⟨ahi, hahi⟩
    simp only [hahi] at hrun
    rcases heap.readU64_of_valid B (2 * m) i hB (by omega) with ⟨blo, hblo⟩
    simp only [hblo] at hrun
    rcases heap.readU64_of_valid B (2 * m) (m + i) hB (by omega) with
      ⟨bhi, hbhi⟩
    simp only [hbhi] at hrun
    let av := dense_upoly_zp_nmod_add_ir this alo ahi
    rcases heap.writeU64_of_valid t1 m i av hT1 hi with ⟨heap1, hw1⟩
    simp only [av, hw1] at hrun
    have hread1 := RawHeap.readU64_writeU64_ne heap heap1 t1 guard i
      readIndex av old hw1 hread (hout1 i (Nat.le_refl _) hi)
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 t1 i av hw1
    let bv := dense_upoly_zp_nmod_add_ir this blo bhi
    rcases heap1.writeU64_of_valid t2 m i bv ((hlayout1 t2 m).mp hT2) hi with
      ⟨heap2, hw2⟩
    simp only [bv, hw2] at hrun
    have hread2 := RawHeap.readU64_writeU64_ne heap1 heap2 t2 guard i
      readIndex bv old hw2 hread1 (hout2 i (Nat.le_refl _) hi)
    have hlayout2 := RawHeap.writeU64_sameLayout heap1 heap2 t2 i bv hw2
    apply karAddHalvesLoop_preserves_outside this A B t1 t2 guard m (i + 1)
      readIndex heap2 heap' old
      ((hlayout2 A (2 * m)).mp ((hlayout1 A (2 * m)).mp hA))
      ((hlayout2 B (2 * m)).mp ((hlayout1 B (2 * m)).mp hB))
      ((hlayout2 t1 m).mp ((hlayout1 t1 m).mp hT1))
      ((hlayout2 t2 m).mp ((hlayout1 t2 m).mp hT2)) hread2
    · intro j hij hjm
      exact hout1 j (by omega) hjm
    · intro j hij hjm
      exact hout2 j (by omega) hjm
    · exact hrun
  next hnot =>
    have heq : heap' = heap := Except.ok.inj hrun.symm
    simpa [heq] using hread
termination_by m - i
decreasing_by omega

theorem karAddHalvesLoop_current_values (this : DenseUPolyZp)
    (A B t1 t2 : RawPtr UInt64) (m i : Nat) (heap heap' : RawHeap)
    (hA : heap.ValidU64Slice A (2 * m))
    (hB : heap.ValidU64Slice B (2 * m))
    (hT1 : heap.ValidU64Slice t1 m)
    (hT2 : heap.ValidU64Slice t2 m)
    (hT1T2 : U64SlicesDisjoint t1 m t2 m)
    (hi : i < m)
    (hrun : karAddHalvesLoop this A B t1 t2 m i heap = .ok heap') :
    ∃ alo ahi blo bhi,
      heap.readU64 A i = .ok alo ∧
      heap.readU64 A (m + i) = .ok ahi ∧
      heap.readU64 B i = .ok blo ∧
      heap.readU64 B (m + i) = .ok bhi ∧
      heap'.readU64 t1 i =
        .ok (dense_upoly_zp_nmod_add_ir this alo ahi) ∧
      heap'.readU64 t2 i =
        .ok (dense_upoly_zp_nmod_add_ir this blo bhi) := by
  unfold karAddHalvesLoop at hrun
  simp only [hi, ↓reduceDIte] at hrun
  rcases heap.readU64_of_valid A (2 * m) i hA (by omega) with ⟨alo, halo⟩
  simp only [halo] at hrun
  rcases heap.readU64_of_valid A (2 * m) (m + i) hA (by omega) with
    ⟨ahi, hahi⟩
  simp only [hahi] at hrun
  rcases heap.readU64_of_valid B (2 * m) i hB (by omega) with ⟨blo, hblo⟩
  simp only [hblo] at hrun
  rcases heap.readU64_of_valid B (2 * m) (m + i) hB (by omega) with
    ⟨bhi, hbhi⟩
  simp only [hbhi] at hrun
  let av := dense_upoly_zp_nmod_add_ir this alo ahi
  rcases heap.writeU64_of_valid t1 m i av hT1 hi with ⟨heap1, hw1⟩
  simp only [av, hw1] at hrun
  have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 t1 i av hw1
  let bv := dense_upoly_zp_nmod_add_ir this blo bhi
  rcases heap1.writeU64_of_valid t2 m i bv ((hlayout1 t2 m).mp hT2) hi with
    ⟨heap2, hw2⟩
  simp only [bv, hw2] at hrun
  have hlayout2 := RawHeap.writeU64_sameLayout heap1 heap2 t2 i bv hw2
  have hav1 := RawHeap.readU64_writeU64_same heap heap1 t1 i av hw1
  have hav2 := RawHeap.readU64_writeU64_ne heap1 heap2 t2 t1 i i bv av
    hw2 hav1 (u64SlicesDisjoint_symm hT1T2 i hi i hi)
  have hbv2 := RawHeap.readU64_writeU64_same heap1 heap2 t2 i bv hw2
  have havFinal := karAddHalvesLoop_preserves_outside this A B t1 t2 t1 m
    (i + 1) i heap2 heap' av
    ((hlayout2 A (2 * m)).mp ((hlayout1 A (2 * m)).mp hA))
    ((hlayout2 B (2 * m)).mp ((hlayout1 B (2 * m)).mp hB))
    ((hlayout2 t1 m).mp ((hlayout1 t1 m).mp hT1))
    ((hlayout2 t2 m).mp ((hlayout1 t2 m).mp hT2)) hav2
    (by intro j _ _; exact Or.inr (by omega))
    (by intro j _ hj; exact u64SlicesDisjoint_symm hT1T2 j hj i hi) hrun
  have hbvFinal := karAddHalvesLoop_preserves_outside this A B t1 t2 t2 m
    (i + 1) i heap2 heap' bv
    ((hlayout2 A (2 * m)).mp ((hlayout1 A (2 * m)).mp hA))
    ((hlayout2 B (2 * m)).mp ((hlayout1 B (2 * m)).mp hB))
    ((hlayout2 t1 m).mp ((hlayout1 t1 m).mp hT1))
    ((hlayout2 t2 m).mp ((hlayout1 t2 m).mp hT2)) hbv2
    (by intro j _ hj; exact hT1T2 j hj i hi)
    (by intro j _ _; exact Or.inr (by omega)) hrun
  exact ⟨alo, ahi, blo, bhi, halo, hahi, hblo, hbhi,
    by simpa [av] using havFinal, by simpa [bv] using hbvFinal⟩

theorem karAddHalvesLoop_current_coeffs (this : DenseUPolyZp)
    (A B t1 t2 : RawPtr UInt64) (m i : Nat) (heap heap' : RawHeap)
    (left right : Polynomial (ZMod this._p.toNat))
    (hp : 1 < this._p.toNat)
    (hA : heap.ValidU64Slice A (2 * m))
    (hB : heap.ValidU64Slice B (2 * m))
    (hT1 : heap.ValidU64Slice t1 m)
    (hT2 : heap.ValidU64Slice t2 m)
    (hT1T2 : U64SlicesDisjoint t1 m t2 m)
    (hCanonicalA : CanonicalU64Prefix heap A (2 * m) this._p)
    (hCanonicalB : CanonicalU64Prefix heap B (2 * m) this._p)
    (hRepA : SlicePolyRep heap A (2 * m) this._p.toNat left)
    (hRepB : SlicePolyRep heap B (2 * m) this._p.toNat right)
    (hi : i < m)
    (hrun : karAddHalvesLoop this A B t1 t2 m i heap = .ok heap') :
    ∃ value1 value2,
      heap'.readU64 t1 i = .ok value1 ∧
      heap'.readU64 t2 i = .ok value2 ∧
      (value1.toNat : ZMod this._p.toNat) =
        left.coeff i + left.coeff (m + i) ∧
      (value2.toNat : ZMod this._p.toNat) =
        right.coeff i + right.coeff (m + i) ∧
      value1.toNat < this._p.toNat ∧ value2.toNat < this._p.toNat := by
  rcases karAddHalvesLoop_current_values this A B t1 t2 m i heap heap'
      hA hB hT1 hT2 hT1T2 hi hrun with
    ⟨alo, ahi, blo, bhi, halo, hahi, hblo, hbhi, ht1, ht2⟩
  have haloLt : alo < this._p := by
    simpa [UInt64.lt_iff_toNat_lt] using hCanonicalA i alo (by omega) halo
  have hahiLt : ahi < this._p := by
    simpa [UInt64.lt_iff_toNat_lt] using hCanonicalA (m + i) ahi (by omega) hahi
  have hbloLt : blo < this._p := by
    simpa [UInt64.lt_iff_toNat_lt] using hCanonicalB i blo (by omega) hblo
  have hbhiLt : bhi < this._p := by
    simpa [UInt64.lt_iff_toNat_lt] using hCanonicalB (m + i) bhi (by omega) hbhi
  have hpWord : this._p ≠ 0 := by
    intro hzero
    have hzeroNat := congrArg UInt64.toNat hzero
    simp at hzeroNat
    omega
  rcases slicePolyRep_coeff heap A (2 * m) this._p.toNat left hRepA i
      (by omega) with ⟨alo', halo', hcoeffAlo⟩
  have haloEq : alo' = alo := Except.ok.inj (halo'.symm.trans halo)
  subst alo'
  rcases slicePolyRep_coeff heap A (2 * m) this._p.toNat left hRepA (m + i)
      (by omega) with ⟨ahi', hahi', hcoeffAhi⟩
  have hahiEq : ahi' = ahi := Except.ok.inj (hahi'.symm.trans hahi)
  subst ahi'
  rcases slicePolyRep_coeff heap B (2 * m) this._p.toNat right hRepB i
      (by omega) with ⟨blo', hblo', hcoeffBlo⟩
  have hbloEq : blo' = blo := Except.ok.inj (hblo'.symm.trans hblo)
  subst blo'
  rcases slicePolyRep_coeff heap B (2 * m) this._p.toNat right hRepB (m + i)
      (by omega) with ⟨bhi', hbhi', hcoeffBhi⟩
  have hbhiEq : bhi' = bhi := Except.ok.inj (hbhi'.symm.trans hbhi)
  subst bhi'
  let value1 := dense_upoly_zp_nmod_add_ir this alo ahi
  let value2 := dense_upoly_zp_nmod_add_ir this blo bhi
  refine ⟨value1, value2, by simpa [value1] using ht1,
    by simpa [value2] using ht2, ?_, ?_, ?_, ?_⟩
  · simpa [value1, hcoeffAlo, hcoeffAhi] using
      nmodAdd_cast this alo ahi hpWord haloLt hahiLt
  · simpa [value2, hcoeffBlo, hcoeffBhi] using
      nmodAdd_cast this blo bhi hpWord hbloLt hbhiLt
  · simpa [value1, UInt64.lt_iff_toNat_lt] using
      nmodAdd_lt this alo ahi hpWord haloLt hahiLt
  · simpa [value2, UInt64.lt_iff_toNat_lt] using
      nmodAdd_lt this blo bhi hpWord hbloLt hbhiLt

theorem karAddHalvesLoop_coeffs (this : DenseUPolyZp)
    (A B t1 t2 : RawPtr UInt64) (m i : Nat) (heap heap' : RawHeap)
    (left right : Polynomial (ZMod this._p.toNat))
    (hp : 1 < this._p.toNat)
    (hA : heap.ValidU64Slice A (2 * m))
    (hB : heap.ValidU64Slice B (2 * m))
    (hT1 : heap.ValidU64Slice t1 m)
    (hT2 : heap.ValidU64Slice t2 m)
    (hT1T2 : U64SlicesDisjoint t1 m t2 m)
    (hT1A : U64SlicesDisjoint t1 m A (2 * m))
    (hT1B : U64SlicesDisjoint t1 m B (2 * m))
    (hT2A : U64SlicesDisjoint t2 m A (2 * m))
    (hT2B : U64SlicesDisjoint t2 m B (2 * m))
    (hCanonicalA : CanonicalU64Prefix heap A (2 * m) this._p)
    (hCanonicalB : CanonicalU64Prefix heap B (2 * m) this._p)
    (hRepA : SlicePolyRep heap A (2 * m) this._p.toNat left)
    (hRepB : SlicePolyRep heap B (2 * m) this._p.toNat right)
    (hrun : karAddHalvesLoop this A B t1 t2 m i heap = .ok heap') :
    ∀ k, i ≤ k → k < m → ∃ value1 value2,
      heap'.readU64 t1 k = .ok value1 ∧
      heap'.readU64 t2 k = .ok value2 ∧
      (value1.toNat : ZMod this._p.toNat) =
        left.coeff k + left.coeff (m + k) ∧
      (value2.toNat : ZMod this._p.toNat) =
        right.coeff k + right.coeff (m + k) ∧
      value1.toNat < this._p.toNat ∧ value2.toNat < this._p.toNat := by
  intro k hik hkm
  by_cases heq : k = i
  · subst k
    exact karAddHalvesLoop_current_coeffs this A B t1 t2 m i heap heap'
      left right hp hA hB hT1 hT2 hT1T2 hCanonicalA hCanonicalB hRepA
      hRepB hkm hrun
  · have hik' : i + 1 ≤ k := by omega
    unfold karAddHalvesLoop at hrun
    have hi : i < m := by omega
    simp only [hi, ↓reduceDIte] at hrun
    rcases heap.readU64_of_valid A (2 * m) i hA (by omega) with ⟨alo, halo⟩
    simp only [halo] at hrun
    rcases heap.readU64_of_valid A (2 * m) (m + i) hA (by omega) with
      ⟨ahi, hahi⟩
    simp only [hahi] at hrun
    rcases heap.readU64_of_valid B (2 * m) i hB (by omega) with ⟨blo, hblo⟩
    simp only [hblo] at hrun
    rcases heap.readU64_of_valid B (2 * m) (m + i) hB (by omega) with
      ⟨bhi, hbhi⟩
    simp only [hbhi] at hrun
    let av := dense_upoly_zp_nmod_add_ir this alo ahi
    rcases heap.writeU64_of_valid t1 m i av hT1 hi with ⟨heap1, hw1⟩
    simp only [av, hw1] at hrun
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 t1 i av hw1
    let bv := dense_upoly_zp_nmod_add_ir this blo bhi
    rcases heap1.writeU64_of_valid t2 m i bv ((hlayout1 t2 m).mp hT2) hi with
      ⟨heap2, hw2⟩
    simp only [bv, hw2] at hrun
    have hlayout2 := RawHeap.writeU64_sameLayout heap1 heap2 t2 i bv hw2
    have hA1 := (hlayout1 A (2 * m)).mp hA
    have hB1 := (hlayout1 B (2 * m)).mp hB
    have hA2 := (hlayout2 A (2 * m)).mp hA1
    have hB2 := (hlayout2 B (2 * m)).mp hB1
    have hT12 := (hlayout2 t1 m).mp ((hlayout1 t1 m).mp hT1)
    have hT22 := (hlayout2 t2 m).mp ((hlayout1 t2 m).mp hT2)
    have hsameA : SameU64Prefix heap heap2 A (2 * m) := by
      intro j old hj hread
      have hread1 := RawHeap.readU64_writeU64_ne heap heap1 t1 A i j av old
        hw1 hread (hT1A i hi j hj)
      exact RawHeap.readU64_writeU64_ne heap1 heap2 t2 A i j bv old
        hw2 hread1 (hT2A i hi j hj)
    have hsameB : SameU64Prefix heap heap2 B (2 * m) := by
      intro j old hj hread
      have hread1 := RawHeap.readU64_writeU64_ne heap heap1 t1 B i j av old
        hw1 hread (hT1B i hi j hj)
      exact RawHeap.readU64_writeU64_ne heap1 heap2 t2 B i j bv old
        hw2 hread1 (hT2B i hi j hj)
    have hCanonicalA2 : CanonicalU64Prefix heap2 A (2 * m) this._p := by
      intro j value hj hread2
      rcases heap.readU64_of_valid A (2 * m) j hA hj with ⟨old, hread⟩
      have hpreserved := hsameA j old hj hread
      have heq : value = old := Except.ok.inj (hread2.symm.trans hpreserved)
      subst value
      exact hCanonicalA j old hj hread
    have hCanonicalB2 : CanonicalU64Prefix heap2 B (2 * m) this._p := by
      intro j value hj hread2
      rcases heap.readU64_of_valid B (2 * m) j hB hj with ⟨old, hread⟩
      have hpreserved := hsameB j old hj hread
      have heq : value = old := Except.ok.inj (hread2.symm.trans hpreserved)
      subst value
      exact hCanonicalB j old hj hread
    have hRepA2 := slicePolyRep_of_same_prefix heap heap2 A (2 * m)
      this._p.toNat left hA hA2 hsameA hRepA
    have hRepB2 := slicePolyRep_of_same_prefix heap heap2 B (2 * m)
      this._p.toNat right hB hB2 hsameB hRepB
    exact karAddHalvesLoop_coeffs this A B t1 t2 m (i + 1) heap2 heap'
      left right hp hA2 hB2 hT12 hT22 hT1T2 hT1A hT1B hT2A hT2B
      hCanonicalA2 hCanonicalB2 hRepA2 hRepB2 hrun k hik' hkm
termination_by m - i
decreasing_by omega

theorem karSubLoop_eq_subCommonLoop (this : DenseUPolyZp)
    (dst sub : RawPtr UInt64) (count i : Nat) (heap : RawHeap) :
    karSubLoop this dst sub count i heap =
      subCommonLoop this dst dst sub count i heap := by
  unfold karSubLoop subCommonLoop
  split
  next hi =>
    rcases hdst : heap.readU64 dst i with fault | a
    · simp [hdst]
    · rcases hsub : heap.readU64 sub i with fault | b
      · simp [hdst, hsub]
      · rcases hwrite : heap.writeU64 dst i
          (dense_upoly_zp_nmod_sub_ir this a b) with fault | heap1
        · simp [hdst, hsub, hwrite]
        · simp only [hdst, hsub, hwrite]
          exact karSubLoop_eq_subCommonLoop this dst sub count (i + 1) heap1
  next hdone => simp [hdone]
termination_by count - i
decreasing_by omega

theorem karSubLoop_value (this : DenseUPolyZp)
    (dst sub : RawPtr UInt64) (count i k : Nat) (heap heap' : RawHeap)
    (a b : UInt64)
    (hDst : heap.ValidU64Slice dst count)
    (hSub : heap.ValidU64Slice sub count)
    (hDstSub : U64SlicesDisjoint dst count sub count)
    (hik : i ≤ k) (hkc : k < count)
    (ha : heap.readU64 dst k = .ok a)
    (hb : heap.readU64 sub k = .ok b)
    (hrun : karSubLoop this dst sub count i heap = .ok heap') :
    heap'.readU64 dst k =
      .ok (dense_upoly_zp_nmod_sub_ir this a b) := by
  unfold karSubLoop at hrun
  split at hrun
  next hi =>
    rcases heap.readU64_of_valid dst count i hDst hi with ⟨ai, hai⟩
    simp only [hai] at hrun
    rcases heap.readU64_of_valid sub count i hSub hi with ⟨bi, hbi⟩
    simp only [hbi] at hrun
    let value := dense_upoly_zp_nmod_sub_ir this ai bi
    rcases heap.writeU64_of_valid dst count i value hDst hi with ⟨heap1, hw⟩
    simp only [value, hw] at hrun
    have hlayout := RawHeap.writeU64_sameLayout heap heap1 dst i value hw
    by_cases heq : k = i
    · subst k
      have haiEq : ai = a := Except.ok.inj (hai.symm.trans ha)
      have hbiEq : bi = b := Except.ok.inj (hbi.symm.trans hb)
      subst ai
      subst bi
      have hnow := RawHeap.readU64_writeU64_same heap heap1 dst i value hw
      apply subCommonLoop_preserves_outside this dst dst sub dst count
        (i + 1) i heap1 heap' value ((hlayout dst count).mp hDst)
        ((hlayout dst count).mp hDst) ((hlayout sub count).mp hSub) hnow
      · intro j _ _
        exact Or.inr (by omega)
      · rw [← karSubLoop_eq_subCommonLoop]
        simpa [value] using hrun
    · have hik' : i + 1 ≤ k := by omega
      have ha1 := RawHeap.readU64_writeU64_ne heap heap1 dst dst i k value a
        hw ha (Or.inr (by omega))
      have hb1 := RawHeap.readU64_writeU64_ne heap heap1 dst sub i k value b
        hw hb (hDstSub i hi k hkc)
      exact karSubLoop_value this dst sub count (i + 1) k heap1 heap'
        a b ((hlayout dst count).mp hDst) ((hlayout sub count).mp hSub)
        hDstSub hik' hkc ha1 hb1 (by simpa [value] using hrun)
  next hdone => omega
termination_by count - i
decreasing_by omega

theorem karSubLoop_preserves_outside (this : DenseUPolyZp)
    (dst sub guard : RawPtr UInt64) (count i readIndex : Nat)
    (heap heap' : RawHeap) (old : UInt64)
    (hDst : heap.ValidU64Slice dst count)
    (hSub : heap.ValidU64Slice sub count)
    (hread : heap.readU64 guard readIndex = .ok old)
    (hout : ∀ j, i ≤ j → j < count → dst.region ≠ guard.region ∨
      dst.limbOffset + j ≠ guard.limbOffset + readIndex)
    (hrun : karSubLoop this dst sub count i heap = .ok heap') :
    heap'.readU64 guard readIndex = .ok old := by
  unfold karSubLoop at hrun
  split at hrun
  next hi =>
    rcases heap.readU64_of_valid dst count i hDst hi with ⟨a, ha⟩
    simp only [ha] at hrun
    rcases heap.readU64_of_valid sub count i hSub hi with ⟨b, hb⟩
    simp only [hb] at hrun
    let value := dense_upoly_zp_nmod_sub_ir this a b
    rcases heap.writeU64_of_valid dst count i value hDst hi with ⟨heap1, hw⟩
    simp only [value, hw] at hrun
    have hread1 := RawHeap.readU64_writeU64_ne heap heap1 dst guard i
      readIndex value old hw hread (hout i (Nat.le_refl _) hi)
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 dst i value hw
    apply karSubLoop_preserves_outside this dst sub guard count (i + 1)
      readIndex heap1 heap' old ((hlayout1 dst count).mp hDst)
      ((hlayout1 sub count).mp hSub) hread1
    · intro j hij hjc
      exact hout j (by omega) hjc
    · exact hrun
  next hnot =>
    have heq : heap' = heap := Except.ok.inj hrun.symm
    simpa [heq] using hread
termination_by count - i
decreasing_by omega

theorem karSubLoop_preserves_prefix (this : DenseUPolyZp)
    (dst sub guard : RawPtr UInt64) (count guardLength : Nat)
    (heap heap' : RawHeap)
    (hDst : heap.ValidU64Slice dst count)
    (hSub : heap.ValidU64Slice sub count)
    (hDstGuard : U64SlicesDisjoint dst count guard guardLength)
    (hrun : karSubLoop this dst sub count 0 heap = .ok heap') :
    SameU64Prefix heap heap' guard guardLength := by
  intro readIndex old hreadIndex hread
  exact karSubLoop_preserves_outside this dst sub guard count 0 readIndex
    heap heap' old hDst hSub hread
    (by
      intro writeIndex _ hwriteIndex
      exact hDstGuard writeIndex hwriteIndex readIndex hreadIndex)
    hrun

theorem karSubLoop_refines_slice (this : DenseUPolyZp)
    (dst sub : RawPtr UInt64) (count : Nat) (heap heap' : RawHeap)
    (left right : Polynomial (ZMod this._p.toNat))
    (hp : this._p ≠ 0)
    (hDst : heap.ValidU64Slice dst count)
    (hSub : heap.ValidU64Slice sub count)
    (hDstSub : U64SlicesDisjoint dst count sub count)
    (hCanonicalDst : CanonicalU64Prefix heap dst count this._p)
    (hCanonicalSub : CanonicalU64Prefix heap sub count this._p)
    (hRepDst : SlicePolyRep heap dst count this._p.toNat left)
    (hRepSub : SlicePolyRep heap sub count this._p.toNat right)
    (hrun : karSubLoop this dst sub count 0 heap = .ok heap') :
    RawHeap.SameLayout heap heap' ∧
      SlicePolyRep heap' dst count this._p.toNat (left - right) ∧
      CanonicalU64Prefix heap' dst count this._p := by
  rcases karSubLoop_ok this dst sub count 0 heap hDst hSub with
    ⟨okHeap, hok, hlayout⟩
  have heq : okHeap = heap' := Except.ok.inj (hok.symm.trans hrun)
  subst okHeap
  have hDst' := (hlayout dst count).mp hDst
  rcases slicePolyRep_exists_unique heap' dst count this._p.toNat hDst' with
    ⟨output, hRepOutput, _⟩
  have houtput : output = left - right := by
    ext degree
    by_cases hd : degree < count
    · rcases slicePolyRep_coeff heap dst count this._p.toNat left hRepDst
        degree hd with ⟨a, ha, hcoeffA⟩
      rcases slicePolyRep_coeff heap sub count this._p.toNat right hRepSub
        degree hd with ⟨b, hb, hcoeffB⟩
      rcases slicePolyRep_coeff heap' dst count this._p.toNat output
        hRepOutput degree hd with ⟨c, hc, hcoeffC⟩
      have hvalue := karSubLoop_value this dst sub count 0 degree heap heap'
        a b hDst hSub hDstSub (by omega) hd ha hb hrun
      have hcEq : c = dense_upoly_zp_nmod_sub_ir this a b :=
        Except.ok.inj (hc.symm.trans hvalue)
      have haLt : a < this._p := by
        simpa [UInt64.lt_iff_toNat_lt] using
          hCanonicalDst degree a hd ha
      have hbLt : b < this._p := by
        simpa [UInt64.lt_iff_toNat_lt] using
          hCanonicalSub degree b hd hb
      rw [Polynomial.coeff_sub, hcoeffC, hcEq,
        nmodSub_cast this a b hp haLt hbLt, hcoeffA, hcoeffB]
    · rw [slicePolyRep_coeff_zero_of_length_le heap' dst count
          this._p.toNat output hRepOutput degree (by omega),
        Polynomial.coeff_sub,
        slicePolyRep_coeff_zero_of_length_le heap dst count
          this._p.toNat left hRepDst degree (by omega),
        slicePolyRep_coeff_zero_of_length_le heap sub count
          this._p.toNat right hRepSub degree (by omega), sub_zero]
  subst output
  refine ⟨hlayout, hRepOutput, ?_⟩
  intro k value hk hread
  rcases slicePolyRep_coeff heap dst count this._p.toNat left hRepDst k hk with
    ⟨a, ha, _⟩
  rcases slicePolyRep_coeff heap sub count this._p.toNat right hRepSub k hk with
    ⟨b, hb, _⟩
  have hvalue := karSubLoop_value this dst sub count 0 k heap heap'
    a b hDst hSub hDstSub (by omega) hk ha hb hrun
  have heqValue : value = dense_upoly_zp_nmod_sub_ir this a b :=
    Except.ok.inj (hread.symm.trans hvalue)
  subst value
  exact nmodSub_lt this a b hp
    (by simpa [UInt64.lt_iff_toNat_lt] using hCanonicalDst k a hk ha)
    (by simpa [UInt64.lt_iff_toNat_lt] using hCanonicalSub k b hk hb)

theorem karAssembleLoop_preserves_outside (this : DenseUPolyZp)
    (C sP1 guard : RawPtr UInt64) (m count i readIndex : Nat)
    (heap heap' : RawHeap) (old : UInt64)
    (hC : heap.ValidU64Slice C (m + count))
    (hP1 : heap.ValidU64Slice sP1 count)
    (hread : heap.readU64 guard readIndex = .ok old)
    (hout : ∀ j, i ≤ j → j < count → C.region ≠ guard.region ∨
      C.limbOffset + (m + j) ≠ guard.limbOffset + readIndex)
    (hrun : karAssembleLoop this C sP1 m count i heap = .ok heap') :
    heap'.readU64 guard readIndex = .ok old := by
  unfold karAssembleLoop at hrun
  split at hrun
  next hi =>
    rcases heap.readU64_of_valid C (m + count) (m + i) hC (by omega) with
      ⟨base, hbase⟩
    simp only [hbase] at hrun
    rcases heap.readU64_of_valid sP1 count i hP1 hi with ⟨cross, hcross⟩
    simp only [hcross] at hrun
    let value := dense_upoly_zp_nmod_add_ir this base cross
    rcases heap.writeU64_of_valid C (m + count) (m + i) value hC (by omega) with
      ⟨heap1, hw⟩
    simp only [value, hw] at hrun
    have hread1 := RawHeap.readU64_writeU64_ne heap heap1 C guard (m + i)
      readIndex value old hw hread (hout i (Nat.le_refl _) hi)
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 C (m + i) value hw
    apply karAssembleLoop_preserves_outside this C sP1 guard m count (i + 1)
      readIndex heap1 heap' old ((hlayout1 C (m + count)).mp hC)
      ((hlayout1 sP1 count).mp hP1) hread1
    · intro j hij hjc
      exact hout j (by omega) hjc
    · exact hrun
  next hnot =>
    have heq : heap' = heap := Except.ok.inj hrun.symm
    simpa [heq] using hread
termination_by count - i
decreasing_by omega

theorem karAssembleLoop_preserves_prefix (this : DenseUPolyZp)
    (C sP1 guard : RawPtr UInt64) (m count guardLength : Nat)
    (heap heap' : RawHeap)
    (hC : heap.ValidU64Slice C (m + count))
    (hP1 : heap.ValidU64Slice sP1 count)
    (hCGuard : U64SlicesDisjoint C (m + count) guard guardLength)
    (hrun : karAssembleLoop this C sP1 m count 0 heap = .ok heap') :
    SameU64Prefix heap heap' guard guardLength := by
  intro readIndex old hreadIndex hread
  exact karAssembleLoop_preserves_outside this C sP1 guard m count 0
    readIndex heap heap' old hC hP1 hread
    (by
      intro writeIndex _ hwriteIndex
      exact hCGuard (m + writeIndex) (by omega) readIndex hreadIndex)
    hrun

theorem karAssembleLoop_value (this : DenseUPolyZp)
    (C sP1 : RawPtr UInt64) (m count i k : Nat)
    (heap heap' : RawHeap) (base cross : UInt64)
    (hC : heap.ValidU64Slice C (m + count))
    (hP1 : heap.ValidU64Slice sP1 count)
    (hCP1 : U64SlicesDisjoint C (m + count) sP1 count)
    (hik : i ≤ k) (hkc : k < count)
    (hbase : heap.readU64 C (m + k) = .ok base)
    (hcross : heap.readU64 sP1 k = .ok cross)
    (hrun : karAssembleLoop this C sP1 m count i heap = .ok heap') :
    heap'.readU64 C (m + k) =
      .ok (dense_upoly_zp_nmod_add_ir this base cross) := by
  unfold karAssembleLoop at hrun
  split at hrun
  next hi =>
    rcases heap.readU64_of_valid C (m + count) (m + i) hC (by omega) with
      ⟨baseI, hbaseI⟩
    simp only [hbaseI] at hrun
    rcases heap.readU64_of_valid sP1 count i hP1 hi with
      ⟨crossI, hcrossI⟩
    simp only [hcrossI] at hrun
    let value := dense_upoly_zp_nmod_add_ir this baseI crossI
    rcases heap.writeU64_of_valid C (m + count) (m + i) value hC
      (by omega) with ⟨heap1, hw⟩
    simp only [value, hw] at hrun
    have hlayout := RawHeap.writeU64_sameLayout heap heap1 C (m + i) value hw
    by_cases heq : k = i
    · subst k
      have hbaseEq : baseI = base := Except.ok.inj (hbaseI.symm.trans hbase)
      have hcrossEq : crossI = cross := Except.ok.inj
        (hcrossI.symm.trans hcross)
      subst baseI
      subst crossI
      have hnow := RawHeap.readU64_writeU64_same heap heap1 C (m + i)
        value hw
      apply karAssembleLoop_preserves_outside this C sP1 C m count (i + 1)
        (m + i) heap1 heap' value ((hlayout C (m + count)).mp hC)
        ((hlayout sP1 count).mp hP1) hnow
      · intro j _ _
        exact Or.inr (by omega)
      · simpa [value] using hrun
    · have hik' : i + 1 ≤ k := by omega
      have hbase1 := RawHeap.readU64_writeU64_ne heap heap1 C C (m + i)
        (m + k) value base hw hbase (Or.inr (by omega))
      have hcross1 := RawHeap.readU64_writeU64_ne heap heap1 C sP1
        (m + i) k value cross hw hcross (hCP1 (m + i) (by omega) k hkc)
      exact karAssembleLoop_value this C sP1 m count (i + 1) k heap1
        heap' base cross ((hlayout C (m + count)).mp hC)
        ((hlayout sP1 count).mp hP1) hCP1 hik' hkc hbase1 hcross1
        (by simpa [value] using hrun)
  next hdone => omega
termination_by count - i
decreasing_by omega

theorem karAssembleLoop_refines_slice (this : DenseUPolyZp)
    (C sP1 : RawPtr UInt64) (m count : Nat) (heap heap' : RawHeap)
    (basePoly crossPoly : Polynomial (ZMod this._p.toNat))
    (hp : this._p ≠ 0)
    (hC : heap.ValidU64Slice C (m + count))
    (hP1 : heap.ValidU64Slice sP1 count)
    (hCP1 : U64SlicesDisjoint C (m + count) sP1 count)
    (hCanonicalC : CanonicalU64Prefix heap C (m + count) this._p)
    (hCanonicalP1 : CanonicalU64Prefix heap sP1 count this._p)
    (hRepC : SlicePolyRep heap C (m + count) this._p.toNat basePoly)
    (hRepP1 : SlicePolyRep heap sP1 count this._p.toNat crossPoly)
    (hrun : karAssembleLoop this C sP1 m count 0 heap = .ok heap') :
    RawHeap.SameLayout heap heap' ∧
      SlicePolyRep heap' C (m + count) this._p.toNat
        (basePoly + Polynomial.X ^ m * crossPoly) ∧
      CanonicalU64Prefix heap' C (m + count) this._p := by
  rcases karAssembleLoop_ok this C sP1 m count 0 heap hC hP1 with
    ⟨okHeap, hok, hlayout⟩
  have heq : okHeap = heap' := Except.ok.inj (hok.symm.trans hrun)
  subst okHeap
  have hC' := (hlayout C (m + count)).mp hC
  rcases slicePolyRep_exists_unique heap' C (m + count) this._p.toNat hC' with
    ⟨output, hRepOutput, _⟩
  have houtput : output = basePoly + Polynomial.X ^ m * crossPoly := by
    ext degree
    by_cases hd : degree < m + count
    · rcases slicePolyRep_coeff heap C (m + count) this._p.toNat basePoly
        hRepC degree hd with ⟨base, hbase, hcoeffBase⟩
      rcases slicePolyRep_coeff heap' C (m + count) this._p.toNat output
        hRepOutput degree hd with ⟨result, hresult, hcoeffResult⟩
      by_cases hdm : degree < m
      · have hpreserved := karAssembleLoop_preserves_outside this C sP1 C
          m count 0 degree heap heap' base hC hP1 hbase
          (by
            intro j _ _
            exact Or.inr (by omega)) hrun
        have hresultEq : result = base :=
          Except.ok.inj (hresult.symm.trans hpreserved)
        rw [Polynomial.coeff_add, Polynomial.coeff_X_pow_mul', if_neg (by omega),
          hcoeffResult, hresultEq, hcoeffBase, add_zero]
      · let k := degree - m
        have hk : k < count := by
          dsimp [k]
          omega
        have hdegree : m + k = degree := by
          dsimp [k]
          omega
        rcases slicePolyRep_coeff heap sP1 count this._p.toNat crossPoly
          hRepP1 k hk with ⟨cross, hcross, hcoeffCross⟩
        have hvalue := karAssembleLoop_value this C sP1 m count 0 k heap
          heap' base cross hC hP1 hCP1 (by omega) hk
          (by simpa [hdegree] using hbase) hcross hrun
        have hresultEq : result = dense_upoly_zp_nmod_add_ir this base cross :=
          Except.ok.inj (hresult.symm.trans (by simpa [hdegree] using hvalue))
        have hbaseLt : base < this._p := by
          simpa [UInt64.lt_iff_toNat_lt] using
            hCanonicalC degree base hd hbase
        have hcrossLt : cross < this._p := by
          simpa [UInt64.lt_iff_toNat_lt] using
            hCanonicalP1 k cross hk hcross
        rw [Polynomial.coeff_add, Polynomial.coeff_X_pow_mul', if_pos (by omega),
          hcoeffResult, hresultEq, nmodAdd_cast this base cross hp hbaseLt
            hcrossLt,
          hcoeffBase]
        simpa [k, hdegree] using hcoeffCross.symm
    · rw [slicePolyRep_coeff_zero_of_length_le heap' C (m + count)
          this._p.toNat output hRepOutput degree (by omega),
        Polynomial.coeff_add,
        slicePolyRep_coeff_zero_of_length_le heap C (m + count)
          this._p.toNat basePoly hRepC degree (by omega),
        Polynomial.coeff_X_pow_mul', if_pos (by omega),
        slicePolyRep_coeff_zero_of_length_le heap sP1 count
          this._p.toNat crossPoly hRepP1 (degree - m) (by omega), add_zero]
  subst output
  refine ⟨hlayout, hRepOutput, ?_⟩
  intro degree value hd hread
  rcases slicePolyRep_coeff heap C (m + count) this._p.toNat basePoly
    hRepC degree hd with ⟨base, hbase, _⟩
  by_cases hdm : degree < m
  · have hpreserved := karAssembleLoop_preserves_outside this C sP1 C
      m count 0 degree heap heap' base hC hP1 hbase
      (by
        intro j _ _
        exact Or.inr (by omega)) hrun
    have hvalueEq : value = base := Except.ok.inj (hread.symm.trans hpreserved)
    subst value
    exact hCanonicalC degree base hd hbase
  · let k := degree - m
    have hk : k < count := by
      dsimp [k]
      omega
    have hdegree : m + k = degree := by
      dsimp [k]
      omega
    rcases slicePolyRep_coeff heap sP1 count this._p.toNat crossPoly
      hRepP1 k hk with ⟨cross, hcross, _⟩
    have hresult := karAssembleLoop_value this C sP1 m count 0 k heap
      heap' base cross hC hP1 hCP1 (by omega) hk
      (by simpa [hdegree] using hbase) hcross hrun
    have hvalueEq : value = dense_upoly_zp_nmod_add_ir this base cross :=
      Except.ok.inj (hread.symm.trans (by simpa [hdegree] using hresult))
    subst value
    exact nmodAdd_lt this base cross hp
      (by simpa [UInt64.lt_iff_toNat_lt] using hCanonicalC degree base hd hbase)
      (by simpa [UInt64.lt_iff_toNat_lt] using hCanonicalP1 k cross hk hcross)

theorem karPrepareHalves_preserves_outside (this : DenseUPolyZp)
    (A B t1 t2 guard : RawPtr UInt64) (m h guardLen : Nat)
    (heap heap' : RawHeap)
    (hmh : m ≤ h)
    (hA : heap.ValidU64Slice A (m + h))
    (hB : heap.ValidU64Slice B (m + h))
    (hT1 : heap.ValidU64Slice t1 h)
    (hT2 : heap.ValidU64Slice t2 h)
    (hT1Guard : U64SlicesDisjoint t1 h guard guardLen)
    (hT2Guard : U64SlicesDisjoint t2 h guard guardLen)
    (hrun : karPrepareHalves this A B t1 t2 m h heap = .ok heap') :
    SameU64Prefix heap heap' guard guardLen := by
  have hA2m := heap.validU64Slice_mono A (m + h) (2 * m) hA (by omega)
  have hB2m := heap.validU64Slice_mono B (m + h) (2 * m) hB (by omega)
  have hT1m := heap.validU64Slice_mono t1 h m hT1 hmh
  have hT2m := heap.validU64Slice_mono t2 h m hT2 hmh
  rcases karAddHalvesLoop_ok this A B t1 t2 m 0 heap hA2m hB2m hT1m
      hT2m with ⟨heap1, hadd, hlayout1⟩
  have htail : karOddTail A B t1 t2 m h heap1 = .ok heap' := by
    simpa [karPrepareHalves, hadd] using hrun
  have hsameAdd : SameU64Prefix heap heap1 guard guardLen := by
    intro i old hi hread
    apply karAddHalvesLoop_preserves_outside this A B t1 t2 guard m 0 i
      heap heap1 old hA2m hB2m hT1m hT2m hread
    · intro _ _ _
      exact hT1Guard _ (by omega) i hi
    · intro _ _ _
      exact hT2Guard _ (by omega) i hi
    · exact hadd
  have hsameTail := karOddTail_preserves_outside A B t1 t2 guard m h
    guardLen heap1 heap'
    ((hlayout1 A (m + h)).mp hA) ((hlayout1 B (m + h)).mp hB)
    ((hlayout1 t1 h).mp hT1) ((hlayout1 t2 h).mp hT2)
    hT1Guard hT2Guard htail
  intro i old hi hread
  exact hsameTail i old hi (hsameAdd i old hi hread)

/-- Mathematical value of the exact raw cells visited by the C++ dot loop.
It shares the loop's reads and failure behavior, but performs unbounded natural
addition so that machine accumulation can be related to it explicitly. -/
def classicalDotNat (heap : RawHeap) (A B : RawPtr UInt64)
    (k stop j : Nat) : RawExec Nat :=
  if h : j ≤ stop then
    match heap.readU64 A j with
    | .error fault => .error fault
    | .ok a =>
      match heap.readU64 B (k - j) with
      | .error fault => .error fault
      | .ok b =>
        match classicalDotNat heap A B k stop (j + 1) with
        | .error fault => .error fault
        | .ok tail => .ok (a.toNat * b.toNat + tail)
  else
    .ok 0
termination_by stop + 1 - j
decreasing_by omega

def classicalDotPoly {p : Nat} (left right : Polynomial (ZMod p))
    (k stop j : Nat) : ZMod p :=
  if h : j ≤ stop then
    left.coeff j * right.coeff (k - j) +
      classicalDotPoly left right k stop (j + 1)
  else
    0
termination_by stop + 1 - j
decreasing_by omega

theorem classicalDotPoly_eq_sum_Icc {p : Nat}
    (left right : Polynomial (ZMod p)) (k stop j : Nat) :
    classicalDotPoly left right k stop j =
      ∑ t ∈ Finset.Icc j stop, left.coeff t * right.coeff (k - t) := by
  rw [classicalDotPoly]
  split
  next hle =>
    rw [← Finset.insert_Icc_succ_left_eq_Icc hle]
    simp
    exact classicalDotPoly_eq_sum_Icc left right k stop (j + 1)
  next hnot =>
    symm
    apply Finset.sum_eq_zero
    intro t ht
    simp only [Finset.mem_Icc] at ht
    omega
termination_by stop + 1 - j
decreasing_by omega

theorem classicalDotPoly_source_eq_coeff {p : Nat}
    (heap : RawHeap) (A B : RawPtr UInt64) (lenA lenB k : Nat)
    (left right : Polynomial (ZMod p))
    (hLenA : 0 < lenA)
    (hRepA : SlicePolyRep heap A lenA p left)
    (hRepB : SlicePolyRep heap B lenB p right) :
    classicalDotPoly left right k
      (if k < lenA then k else lenA - 1)
      (if k ≥ lenB then k - lenB + 1 else 0) =
        (left * right).coeff k := by
  rw [classicalDotPoly_eq_sum_Icc]
  rw [Polynomial.coeff_mul,
    Finset.Nat.sum_antidiagonal_eq_sum_range_succ_mk]
  by_cases hkA : k < lenA <;> by_cases hkB : k ≥ lenB
  all_goals simp only [hkA, hkB, if_true, if_false]
  all_goals
  apply Finset.sum_subset
  · intro t ht
    simp only [Finset.mem_Icc] at ht
    simp only [Finset.mem_range]
    omega
  · intro t htRange htIcc
    simp only [Finset.mem_range] at htRange
    simp only [Finset.mem_Icc, not_and_or, not_le] at htIcc
    rcases htIcc with hBelow | hAbove
    · by_cases hkt : lenB ≤ k
      · have hzero := slicePolyRep_coeff_zero_of_length_le heap B lenB p
          right hRepB (k - t) (by omega)
        rw [hzero, mul_zero]
      · omega
    · by_cases hkt : k < lenA
      · omega
      · have hzero := slicePolyRep_coeff_zero_of_length_le heap A lenA p
          left hRepA t (by omega)
        rw [hzero, zero_mul]

theorem classicalDotNat_cast_eq_poly (heap : RawHeap)
    (A B : RawPtr UInt64) (lenA lenB p k stop j sum : Nat)
    (left right : Polynomial (ZMod p))
    (hRepA : SlicePolyRep heap A lenA p left)
    (hRepB : SlicePolyRep heap B lenB p right)
    (hAIndex : ∀ t, j ≤ t → t ≤ stop → t < lenA)
    (hBIndex : ∀ t, j ≤ t → t ≤ stop → k - t < lenB)
    (hrun : classicalDotNat heap A B k stop j = .ok sum) :
    (sum : ZMod p) = classicalDotPoly left right k stop j := by
  unfold classicalDotNat at hrun
  split at hrun
  next hle =>
    rw [classicalDotPoly, dif_pos hle]
    have hjA := hAIndex j (Nat.le_refl _) hle
    have hjB := hBIndex j (Nat.le_refl _) hle
    rcases slicePolyRep_coeff heap A lenA p left hRepA j hjA with
      ⟨a, ha, hcoeffA⟩
    simp only [ha] at hrun
    rcases slicePolyRep_coeff heap B lenB p right hRepB (k - j) hjB with
      ⟨b, hb, hcoeffB⟩
    simp only [hb] at hrun
    cases ht : classicalDotNat heap A B k stop (j + 1) with
    | error fault => simp [ht] at hrun
    | ok tail =>
      simp only [ht] at hrun
      have hsum : sum = a.toNat * b.toNat + tail :=
        Except.ok.inj hrun.symm
      have htail := classicalDotNat_cast_eq_poly heap A B lenA lenB p k
        stop (j + 1) tail left right hRepA hRepB
        (by intro t hjt hts; exact hAIndex t (by omega) hts)
        (by intro t hjt hts; exact hBIndex t (by omega) hts) ht
      rw [hsum, Nat.cast_add, Nat.cast_mul, hcoeffA, hcoeffB, htail]
  next hnot =>
    rw [classicalDotPoly, dif_neg hnot]
    have hzero : sum = 0 := Except.ok.inj hrun.symm
    subst sum
    simp
termination_by stop + 1 - j
decreasing_by omega

theorem classicalDotLoop_modEq (heap : RawHeap) (A B : RawPtr UInt64)
    (k stop j : Nat) (acc result : Word3) (sum : Nat)
    (hrun : classicalDotLoop heap A B k stop j acc = .ok result)
    (hsum : classicalDotNat heap A B k stop j = .ok sum) :
    Nat.ModEq (limbBase ^ 3) (word3Value result)
      (word3Value acc + sum) := by
  unfold classicalDotLoop at hrun
  unfold classicalDotNat at hsum
  split at hrun
  next hle =>
    simp only [hle, ↓reduceDIte] at hsum
    cases ha : heap.readU64 A j with
    | error fault => simp [ha] at hrun
    | ok a =>
      simp only [ha] at hrun hsum
      cases hb : heap.readU64 B (k - j) with
      | error fault => simp [hb] at hrun
      | ok b =>
        simp only [hb] at hrun hsum
        cases ht : classicalDotNat heap A B k stop (j + 1) with
        | error fault => simp [ht] at hsum
        | ok tail =>
          simp only [ht] at hsum
          have hsumEq : sum = a.toNat * b.toNat + tail :=
            Except.ok.inj hsum.symm
          let product := Generated.StrictGCD.dense_upoly_zp__umul128_ir 0 0 a b
          let acc' := Generated.StrictGCD.dense_upoly_zp__add_carry3_ir
            acc product.fst product.snd
          have hrec := classicalDotLoop_modEq heap A B k stop (j + 1)
            acc' result tail hrun ht
          have hstep : Nat.ModEq (limbBase ^ 3) (word3Value acc')
              (word3Value acc + a.toNat * b.toNat) := by
            simpa [product, acc'] using addMulWord3_modEq acc a b
          have htotal := hrec.trans (hstep.add_right tail)
          simpa [hsumEq, Nat.add_assoc] using htotal
  next hnot =>
    simp only [hnot, ↓reduceDIte] at hsum
    have hresult : result = acc := Except.ok.inj hrun.symm
    have hzero : sum = 0 := Except.ok.inj hsum.symm
    subst result
    subst sum
    simpa using (Nat.ModEq.refl (word3Value acc) :
      Nat.ModEq (limbBase ^ 3) (word3Value acc) (word3Value acc))
termination_by stop + 1 - j
decreasing_by omega

theorem classicalDotNat_ok (heap : RawHeap) (A B : RawPtr UInt64)
    (lenA lenB k stop j : Nat)
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB)
    (hAIndex : ∀ t, j ≤ t → t ≤ stop → t < lenA)
    (hBIndex : ∀ t, j ≤ t → t ≤ stop → k - t < lenB) :
    ∃ sum, classicalDotNat heap A B k stop j = .ok sum := by
  unfold classicalDotNat
  split
  next hle =>
    rcases heap.readU64_of_valid A lenA j hA
      (hAIndex j (Nat.le_refl _) hle) with ⟨a, ha⟩
    simp only [ha]
    rcases heap.readU64_of_valid B lenB (k - j) hB
      (hBIndex j (Nat.le_refl _) hle) with ⟨b, hb⟩
    simp only [hb]
    rcases classicalDotNat_ok heap A B lenA lenB k stop (j + 1)
      hA hB
      (by intro t hjt hts; exact hAIndex t (by omega) hts)
      (by intro t hjt hts; exact hBIndex t (by omega) hts) with
      ⟨tail, htail⟩
    rw [htail]
    exact ⟨a.toNat * b.toNat + tail, rfl⟩
  next => exact ⟨0, rfl⟩
termination_by stop + 1 - j
decreasing_by omega

theorem classicalDotNat_bound (heap : RawHeap) (A B : RawPtr UInt64)
    (lenA lenB k stop j : Nat) (p : UInt64) (sum : Nat)
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB)
    (hCanonicalA : CanonicalU64Prefix heap A lenA p)
    (hCanonicalB : CanonicalU64Prefix heap B lenB p)
    (hAIndex : ∀ t, j ≤ t → t ≤ stop → t < lenA)
    (hBIndex : ∀ t, j ≤ t → t ≤ stop → k - t < lenB)
    (hrun : classicalDotNat heap A B k stop j = .ok sum) :
    sum ≤ (stop + 1 - j) * (p.toNat - 1) ^ 2 := by
  unfold classicalDotNat at hrun
  split at hrun
  next hle =>
    have hjA := hAIndex j (Nat.le_refl _) hle
    have hjB := hBIndex j (Nat.le_refl _) hle
    rcases heap.readU64_of_valid A lenA j hA hjA with ⟨a, ha⟩
    simp only [ha] at hrun
    rcases heap.readU64_of_valid B lenB (k - j) hB hjB with ⟨b, hb⟩
    simp only [hb] at hrun
    cases ht : classicalDotNat heap A B k stop (j + 1) with
    | error fault => simp [ht] at hrun
    | ok tail =>
      simp only [ht] at hrun
      have hsum : sum = a.toNat * b.toNat + tail :=
        Except.ok.inj hrun.symm
      have haLe : a.toNat ≤ p.toNat - 1 := by
        have := hCanonicalA j a hjA ha
        omega
      have hbLe : b.toNat ≤ p.toNat - 1 := by
        have := hCanonicalB (k - j) b hjB hb
        omega
      have hprod : a.toNat * b.toNat ≤ (p.toNat - 1) ^ 2 := by
        rw [pow_two]
        exact Nat.mul_le_mul haLe hbLe
      have htail := classicalDotNat_bound heap A B lenA lenB k stop (j + 1)
        p tail hA hB hCanonicalA hCanonicalB
        (by intro t hjt hts; exact hAIndex t (by omega) hts)
        (by intro t hjt hts; exact hBIndex t (by omega) hts) ht
      rw [hsum]
      calc
        a.toNat * b.toNat + tail ≤
            (p.toNat - 1) ^ 2 +
              (stop + 1 - (j + 1)) * (p.toNat - 1) ^ 2 :=
          Nat.add_le_add hprod htail
        _ = (stop + 1 - j) * (p.toNat - 1) ^ 2 := by
          have : stop + 1 - j = (stop + 1 - (j + 1)) + 1 := by omega
          rw [this]
          ring
  next hnot =>
    have hzero : sum = 0 := Except.ok.inj hrun.symm
    subst sum
    exact Nat.zero_le _
termination_by stop + 1 - j
decreasing_by omega

theorem classicalDotLoop_ok (heap : RawHeap) (A B : RawPtr UInt64)
    (lenA lenB k stop j : Nat) (acc : Word3)
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB)
    (hAIndex : ∀ t, j ≤ t → t ≤ stop → t < lenA)
    (hBIndex : ∀ t, j ≤ t → t ≤ stop → k - t < lenB) :
    ∃ result, classicalDotLoop heap A B k stop j acc = .ok result := by
  unfold classicalDotLoop
  split
  next hle =>
    have hjA := hAIndex j (Nat.le_refl _) hle
    have hjB := hBIndex j (Nat.le_refl _) hle
    rcases heap.readU64_of_valid A lenA j hA hjA with ⟨a, ha⟩
    simp only [ha]
    rcases heap.readU64_of_valid B lenB (k - j) hB hjB with ⟨b, hb⟩
    simp only [hb]
    let product := Generated.StrictGCD.dense_upoly_zp__umul128_ir 0 0 a b
    let acc' := Generated.StrictGCD.dense_upoly_zp__add_carry3_ir
      acc product.fst product.snd
    apply classicalDotLoop_ok heap A B lenA lenB k stop (j + 1) acc'
      hA hB
    · intro t hjt hts
      exact hAIndex t (by omega) hts
    · intro t hjt hts
      exact hBIndex t (by omega) hts
  next => exact ⟨acc, rfl⟩
termination_by stop + 1 - j
decreasing_by omega

theorem classicalDotLoop_raw_sum (heap : RawHeap) (A B : RawPtr UInt64)
    (lenA lenB k stop j : Nat) (acc : Word3)
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB)
    (hAIndex : ∀ t, j ≤ t → t ≤ stop → t < lenA)
    (hBIndex : ∀ t, j ≤ t → t ≤ stop → k - t < lenB) :
    ∃ result sum,
      classicalDotLoop heap A B k stop j acc = .ok result ∧
      classicalDotNat heap A B k stop j = .ok sum ∧
      Nat.ModEq (limbBase ^ 3) (word3Value result)
        (word3Value acc + sum) := by
  rcases classicalDotLoop_ok heap A B lenA lenB k stop j acc hA hB
    hAIndex hBIndex with ⟨result, hrun⟩
  rcases classicalDotNat_ok heap A B lenA lenB k stop j hA hB
    hAIndex hBIndex with ⟨sum, hsum⟩
  exact ⟨result, sum, hrun, hsum,
    classicalDotLoop_modEq heap A B k stop j acc result sum hrun hsum⟩

theorem classicalDotLoop_exact_zero (heap : RawHeap) (A B : RawPtr UInt64)
    (lenA lenB k stop j : Nat) (p : UInt64) (result : Word3) (sum : Nat)
    (hp : 1 < p.toNat)
    (hcount : stop + 1 - j < limbBase)
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB)
    (hCanonicalA : CanonicalU64Prefix heap A lenA p)
    (hCanonicalB : CanonicalU64Prefix heap B lenB p)
    (hAIndex : ∀ t, j ≤ t → t ≤ stop → t < lenA)
    (hBIndex : ∀ t, j ≤ t → t ≤ stop → k - t < lenB)
    (hrun : classicalDotLoop heap A B k stop j
      { lo := 0, mid := 0, hi := 0 } = .ok result)
    (hsum : classicalDotNat heap A B k stop j = .ok sum) :
    word3Value result = sum := by
  have hbound := classicalDotNat_bound heap A B lenA lenB k stop j p sum
    hA hB hCanonicalA hCanonicalB hAIndex hBIndex hsum
  have hlazy := lazyAccumulation_word3_budget p (stop + 1 - j) 0
    hp hcount (by omega)
  have hpB : p.toNat < limbBase := by
    simpa [limbBase] using UInt64.toNat_lt p
  have hsumLt : sum < limbBase ^ 3 := by
    calc
      sum ≤ (stop + 1 - j) * (p.toNat - 1) ^ 2 := hbound
      _ < p.toNat * limbBase ^ 2 := by simpa using hlazy
      _ < limbBase ^ 3 := by
        have hpowPos : 0 < limbBase ^ 2 :=
          pow_pos (by norm_num [limbBase]) 2
        have hmul : p.toNat * limbBase ^ 2 < limbBase * limbBase ^ 2 :=
          Nat.mul_lt_mul_of_pos_right hpB hpowPos
        simpa [pow_succ, Nat.mul_comm, Nat.mul_left_comm,
          Nat.mul_assoc] using hmul
  have hmod := classicalDotLoop_modEq heap A B k stop j
    { lo := 0, mid := 0, hi := 0 } result sum hrun hsum
  have hmod' : Nat.ModEq (limbBase ^ 3) (word3Value result) sum := by
    simpa [word3Value] using hmod
  exact hmod'.eq_of_lt_of_lt (word3Value_lt result) hsumLt

theorem classicalDotReduced_toNat (this : DenseUPolyZp)
    (heap : RawHeap) (A B : RawPtr UInt64)
    (lenA lenB k stop j : Nat) (result : Word3) (sum : Nat)
    (hcfg : DensePreinvConfigured this)
    (hp : 1 < this._p.toNat)
    (hcount : stop + 1 - j < limbBase)
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB)
    (hCanonicalA : CanonicalU64Prefix heap A lenA this._p)
    (hCanonicalB : CanonicalU64Prefix heap B lenB this._p)
    (hAIndex : ∀ t, j ≤ t → t ≤ stop → t < lenA)
    (hBIndex : ∀ t, j ≤ t → t ≤ stop → k - t < lenB)
    (hrun : classicalDotLoop heap A B k stop j
      { lo := 0, mid := 0, hi := 0 } = .ok result)
    (hsum : classicalDotNat heap A B k stop j = .ok sum) :
    (Generated.StrictGCD.dense_upoly_zp__lll_mod_preinv_ir
      result.hi result.mid result.lo this._p this._ninv this._norm).toNat =
      sum % this._p.toNat := by
  have hexact := classicalDotLoop_exact_zero heap A B lenA lenB k stop j
    this._p result sum hp hcount hA hB hCanonicalA hCanonicalB
    hAIndex hBIndex hrun hsum
  have hsumBound := classicalDotNat_bound heap A B lenA lenB k stop j
    this._p sum hA hB hCanonicalA hCanonicalB hAIndex hBIndex hsum
  have hlazy := lazyAccumulation_word3_budget this._p (stop + 1 - j) 0
    hp hcount (by omega)
  have hvalue : word3Value result < this._p.toNat * limbBase ^ 2 := by
    rw [hexact]
    exact lt_of_le_of_lt hsumBound (by simpa using hlazy)
  have hhi : result.hi.toNat < this._p.toNat :=
    word3_hi_lt_of_value_lt result this._p hvalue
  rw [lll_mod_preinv_ir_correct_of_configured this result.hi result.mid
    result.lo hcfg hhi, hexact]

theorem classicalDotReduced_cast (this : DenseUPolyZp)
    (heap : RawHeap) (A B : RawPtr UInt64)
    (lenA lenB k stop j : Nat) (result : Word3) (sum : Nat)
    (hcfg : DensePreinvConfigured this)
    (hp : 1 < this._p.toNat)
    (hcount : stop + 1 - j < limbBase)
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB)
    (hCanonicalA : CanonicalU64Prefix heap A lenA this._p)
    (hCanonicalB : CanonicalU64Prefix heap B lenB this._p)
    (hAIndex : ∀ t, j ≤ t → t ≤ stop → t < lenA)
    (hBIndex : ∀ t, j ≤ t → t ≤ stop → k - t < lenB)
    (hrun : classicalDotLoop heap A B k stop j
      { lo := 0, mid := 0, hi := 0 } = .ok result)
    (hsum : classicalDotNat heap A B k stop j = .ok sum) :
    ((Generated.StrictGCD.dense_upoly_zp__lll_mod_preinv_ir
      result.hi result.mid result.lo this._p this._ninv this._norm).toNat :
        ZMod this._p.toNat) = (sum : ZMod this._p.toNat) := by
  rw [classicalDotReduced_toNat this heap A B lenA lenB k stop j result sum
    hcfg hp hcount hA hB hCanonicalA hCanonicalB hAIndex hBIndex hrun hsum,
    ZMod.natCast_mod]

theorem classical_index_bounds (lenA lenB k j : Nat)
    (hApos : 0 < lenA) (hBpos : 0 < lenB)
    (hk : k < lenA + lenB - 1)
    (hjMin : (if k ≥ lenB then k - lenB + 1 else 0) ≤ j)
    (hjMax : j ≤ (if k < lenA then k else lenA - 1)) :
    j < lenA ∧ k - j < lenB := by
  split at hjMax
  next hkA =>
    constructor
    · omega
    · split at hjMin <;> omega
  next hkNotA =>
    constructor
    · omega
    · split at hjMin <;> omega

theorem classicalReduced_source_eq_coeff (this : DenseUPolyZp)
    (heap : RawHeap) (A B : RawPtr UInt64) (lenA lenB lenC k : Nat)
    (left right : Polynomial (ZMod this._p.toNat)) (acc : Word3)
    (hcfg : DensePreinvConfigured this)
    (hp : 1 < this._p.toNat)
    (hApos : 0 < lenA) (hBpos : 0 < lenB)
    (hLenAWord : lenA < limbBase)
    (hlenC : lenC = lenA + lenB - 1) (hk : k < lenC)
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB)
    (hCanonicalA : CanonicalU64Prefix heap A lenA this._p)
    (hCanonicalB : CanonicalU64Prefix heap B lenB this._p)
    (hRepA : SlicePolyRep heap A lenA this._p.toNat left)
    (hRepB : SlicePolyRep heap B lenB this._p.toNat right)
    (hdot : classicalDotLoop heap A B k
      (if k < lenA then k else lenA - 1)
      (if k ≥ lenB then k - lenB + 1 else 0)
      { lo := 0, mid := 0, hi := 0 } = .ok acc) :
    ((Generated.StrictGCD.dense_upoly_zp__lll_mod_preinv_ir
      acc.hi acc.mid acc.lo this._p this._ninv this._norm).toNat :
        ZMod this._p.toNat) = (left * right).coeff k ∧
      (Generated.StrictGCD.dense_upoly_zp__lll_mod_preinv_ir
        acc.hi acc.mid acc.lo this._p this._ninv this._norm).toNat <
          this._p.toNat := by
  let jMin := if k ≥ lenB then k - lenB + 1 else 0
  let jMax := if k < lenA then k else lenA - 1
  have hkSource : k < lenA + lenB - 1 := by omega
  have hrange : jMin ≤ jMax := by
    dsimp [jMin, jMax]
    split <;> split <;> omega
  have hAIndex : ∀ t, jMin ≤ t → t ≤ jMax → t < lenA := by
    intro t hjt hts
    exact (classical_index_bounds lenA lenB k t hApos hBpos hkSource
      hjt hts).1
  have hBIndex : ∀ t, jMin ≤ t → t ≤ jMax → k - t < lenB := by
    intro t hjt hts
    exact (classical_index_bounds lenA lenB k t hApos hBpos hkSource
      hjt hts).2
  rcases classicalDotNat_ok heap A B lenA lenB k jMax jMin hA hB
    hAIndex hBIndex with ⟨sum, hsum⟩
  have hcount : jMax + 1 - jMin < limbBase := by
    have hjMaxA := hAIndex jMax hrange (Nat.le_refl _)
    omega
  constructor
  · calc
      ((Generated.StrictGCD.dense_upoly_zp__lll_mod_preinv_ir
        acc.hi acc.mid acc.lo this._p this._ninv this._norm).toNat :
          ZMod this._p.toNat) = (sum : ZMod this._p.toNat) :=
        classicalDotReduced_cast this heap A B lenA lenB k jMax jMin acc sum
          hcfg hp hcount hA hB hCanonicalA hCanonicalB hAIndex hBIndex
          (by simpa [jMin, jMax] using hdot) hsum
      _ = classicalDotPoly left right k jMax jMin :=
        classicalDotNat_cast_eq_poly heap A B lenA lenB this._p.toNat k
          jMax jMin sum left right hRepA hRepB hAIndex hBIndex hsum
      _ = (left * right).coeff k := by
        simpa [jMin, jMax] using classicalDotPoly_source_eq_coeff heap A B
          lenA lenB k left right hApos hRepA hRepB
  · rw [classicalDotReduced_toNat this heap A B lenA lenB k jMax jMin acc sum
      hcfg hp hcount hA hB hCanonicalA hCanonicalB hAIndex hBIndex
      (by simpa [jMin, jMax] using hdot) hsum]
    exact Nat.mod_lt sum (by omega)

/-- The generated schoolbook outer loop writes only `C[k..lenC)`.  This
memory fact keeps both source buffers and already-produced coefficients tied
to their original raw cells throughout the actual C++ loop. -/
theorem classicalOuterLoop_preserves_outside (this : DenseUPolyZp)
    (C A B guard : RawPtr UInt64) (lenA lenB lenC k readIndex : Nat)
    (heap heap' : RawHeap) (old : UInt64)
    (hC : heap.ValidU64Slice C lenC)
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB)
    (hread : heap.readU64 guard readIndex = .ok old)
    (houtside : ∀ i, k ≤ i → i < lenC →
      C.region ≠ guard.region ∨
        C.limbOffset + i ≠ guard.limbOffset + readIndex)
    (hrun : classicalOuterLoop this C A B lenA lenB lenC k heap =
      .ok heap') :
    heap'.readU64 guard readIndex = .ok old := by
  unfold classicalOuterLoop at hrun
  split at hrun
  next hk =>
    let jMin := if k ≥ lenB then k - lenB + 1 else 0
    let jMax := if k < lenA then k else lenA - 1
    cases hdot : classicalDotLoop heap A B k jMax jMin
        { lo := 0, mid := 0, hi := 0 } with
    | error fault => simp [jMin, jMax, hdot] at hrun
    | ok acc =>
      simp only [jMin, jMax, hdot] at hrun
      let value := Generated.StrictGCD.dense_upoly_zp__lll_mod_preinv_ir
        acc.hi acc.mid acc.lo this._p this._ninv this._norm
      rcases heap.writeU64_of_valid C lenC k value hC hk with ⟨heap1, hw⟩
      simp only [value, hw] at hrun
      have hread1 := RawHeap.readU64_writeU64_ne heap heap1 C guard
        k readIndex value old hw hread (houtside k (Nat.le_refl _) hk)
      have hlayout := RawHeap.writeU64_sameLayout heap heap1 C k value hw
      apply classicalOuterLoop_preserves_outside this C A B guard lenA lenB
        lenC (k + 1) readIndex heap1 heap' old
        ((hlayout C lenC).mp hC) ((hlayout A lenA).mp hA)
        ((hlayout B lenB).mp hB) hread1
      · intro i hik hil
        exact houtside i (by omega) hil
      · exact hrun
  next hnot =>
    have heq : heap' = heap := Except.ok.inj hrun.symm
    simpa [heq] using hread
termination_by lenC - k
decreasing_by omega

theorem classicalOuterLoop_same_prefix_region_ne (this : DenseUPolyZp)
    (C A B guard : RawPtr UInt64) (lenA lenB lenC k guardLen : Nat)
    (heap heap' : RawHeap)
    (hC : heap.ValidU64Slice C lenC)
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB)
    (hregions : C.region ≠ guard.region)
    (hrun : classicalOuterLoop this C A B lenA lenB lenC k heap =
      .ok heap') :
    SameU64Prefix heap heap' guard guardLen := by
  intro i old hi hread
  apply classicalOuterLoop_preserves_outside this C A B guard lenA lenB
    lenC k i heap heap' old hC hA hB hread
  · intro _ _ _
    exact Or.inl hregions
  · exact hrun

theorem canonicalU64Prefix_of_same_prefix (before after : RawHeap)
    (ptr : RawPtr UInt64) (length : Nat) (p : UInt64)
    (hvalid : before.ValidU64Slice ptr length)
    (hsame : SameU64Prefix before after ptr length)
    (hcanonical : CanonicalU64Prefix before ptr length p) :
    CanonicalU64Prefix after ptr length p := by
  intro k value hk hreadAfter
  rcases before.readU64_of_valid ptr length k
      hvalid hk with ⟨old, hreadBefore⟩
  have hpreserved := hsame k old hk hreadBefore
  have hvalue : value = old := Except.ok.inj (hreadAfter.symm.trans hpreserved)
  subst value
  exact hcanonical k old hk hreadBefore

/-- Coefficient-level output invariant for the prefix already produced by
the C++ schoolbook outer loop. -/
def ClassicalCoeffPrefix {p : Nat} (heap : RawHeap)
    (C : RawPtr UInt64) (upto : Nat) (poly : Polynomial (ZMod p)) : Prop :=
  ∀ i, i < upto → ∃ value : UInt64,
    heap.readU64 C i = .ok value ∧
      (value.toNat : ZMod p) = poly.coeff i ∧ value.toNat < p

theorem slicePolyRep_of_classicalCoeffPrefix {p : Nat}
    (heap : RawHeap) (C : RawPtr UInt64) (length : Nat)
    (poly : Polynomial (ZMod p))
    (hvalid : heap.ValidU64Slice C length)
    (hprefix : ClassicalCoeffPrefix heap C length poly)
    (hzero : ∀ i, length ≤ i → poly.coeff i = 0) :
    SlicePolyRep heap C length p poly := by
  rcases slicePolyRep_exists_unique heap C length p hvalid with
    ⟨observed, hObserved, _⟩
  have heq : observed = poly := by
    ext i
    by_cases hi : i < length
    · rcases hprefix i hi with ⟨value, hread, hcoeff, _⟩
      rcases slicePolyRep_coeff heap C length p observed hObserved i hi with
        ⟨observedValue, hObservedRead, hObservedCoeff⟩
      have hvalue : observedValue = value :=
        Except.ok.inj (hObservedRead.symm.trans hread)
      rw [hObservedCoeff, hvalue, hcoeff]
    · rw [slicePolyRep_coeff_zero_of_length_le heap C length p observed
          hObserved i (by omega),
        hzero i (by omega)]
  rw [heq] at hObserved
  exact hObserved

theorem canonicalU64Prefix_of_classicalCoeffPrefix
    (heap : RawHeap) (C : RawPtr UInt64) (length : Nat)
    (modulus : UInt64) (poly : Polynomial (ZMod modulus.toNat))
    (hprefix : ClassicalCoeffPrefix heap C length poly) :
    CanonicalU64Prefix heap C length modulus := by
  intro i value hi hread
  rcases hprefix i hi with ⟨stored, hstored, _, hstoredLt⟩
  have heq : value = stored := Except.ok.inj (hread.symm.trans hstored)
  subst value
  exact hstoredLt

noncomputable def karHalfSumPoly {p : Nat}
    (poly : Polynomial (ZMod p)) (m : Nat) : Polynomial (ZMod p) :=
  ∑ i : Fin m,
    Polynomial.monomial i.val (poly.coeff i.val + poly.coeff (m + i.val))

theorem coeff_karHalfSumPoly {p : Nat}
    (poly : Polynomial (ZMod p)) (m degree : Nat) :
    (karHalfSumPoly poly m).coeff degree =
      if degree < m then poly.coeff degree + poly.coeff (m + degree) else 0 := by
  classical
  unfold karHalfSumPoly
  rw [Polynomial.finset_sum_coeff]
  by_cases hdegree : degree < m
  · rw [if_pos hdegree, Finset.sum_eq_single ⟨degree, hdegree⟩]
    · simp
    · intro index _ hne
      have hval : index.val ≠ degree := by
        intro heq
        apply hne
        exact Fin.ext heq
      simp [Polynomial.coeff_monomial, hval]
    · simp
  · rw [if_neg hdegree]
    apply Finset.sum_eq_zero
    intro index _
    have hval : index.val ≠ degree := by
      intro heq
      apply hdegree
      simpa [heq] using index.isLt
    simp [Polynomial.coeff_monomial, hval]

theorem karHalfSumPoly_congr_of_prefix {p : Nat}
    (prefixPoly fullPoly : Polynomial (ZMod p)) (m : Nat)
    (hcoeff : ∀ i, i < 2 * m → prefixPoly.coeff i = fullPoly.coeff i) :
    karHalfSumPoly prefixPoly m = karHalfSumPoly fullPoly m := by
  ext degree
  rw [coeff_karHalfSumPoly, coeff_karHalfSumPoly]
  by_cases hdegree : degree < m
  · rw [if_pos hdegree, if_pos hdegree,
      hcoeff degree (by omega), hcoeff (m + degree) (by omega)]
  · rw [if_neg hdegree, if_neg hdegree]

theorem karAddHalvesLoop_refines_slices (this : DenseUPolyZp)
    (A B t1 t2 : RawPtr UInt64) (m : Nat) (heap heap' : RawHeap)
    (left right : Polynomial (ZMod this._p.toNat))
    (hp : 1 < this._p.toNat)
    (hA : heap.ValidU64Slice A (2 * m))
    (hB : heap.ValidU64Slice B (2 * m))
    (hT1 : heap.ValidU64Slice t1 m)
    (hT2 : heap.ValidU64Slice t2 m)
    (hT1T2 : U64SlicesDisjoint t1 m t2 m)
    (hT1A : U64SlicesDisjoint t1 m A (2 * m))
    (hT1B : U64SlicesDisjoint t1 m B (2 * m))
    (hT2A : U64SlicesDisjoint t2 m A (2 * m))
    (hT2B : U64SlicesDisjoint t2 m B (2 * m))
    (hCanonicalA : CanonicalU64Prefix heap A (2 * m) this._p)
    (hCanonicalB : CanonicalU64Prefix heap B (2 * m) this._p)
    (hRepA : SlicePolyRep heap A (2 * m) this._p.toNat left)
    (hRepB : SlicePolyRep heap B (2 * m) this._p.toNat right)
    (hrun : karAddHalvesLoop this A B t1 t2 m 0 heap = .ok heap') :
    SlicePolyRep heap' t1 m this._p.toNat (karHalfSumPoly left m) ∧
      SlicePolyRep heap' t2 m this._p.toNat (karHalfSumPoly right m) ∧
      CanonicalU64Prefix heap' t1 m this._p ∧
      CanonicalU64Prefix heap' t2 m this._p := by
  have hvalues := karAddHalvesLoop_coeffs this A B t1 t2 m 0 heap heap'
    left right hp hA hB hT1 hT2 hT1T2 hT1A hT1B hT2A hT2B
    hCanonicalA hCanonicalB hRepA hRepB hrun
  rcases karAddHalvesLoop_ok this A B t1 t2 m 0 heap hA hB hT1 hT2 with
    ⟨okHeap, hok, hlayout⟩
  have hokHeap : okHeap = heap' := Except.ok.inj (hok.symm.trans hrun)
  subst okHeap
  have hT1' := (hlayout t1 m).mp hT1
  have hT2' := (hlayout t2 m).mp hT2
  have hprefix1 : ClassicalCoeffPrefix heap' t1 m (karHalfSumPoly left m) := by
    intro k hk
    rcases hvalues k (Nat.zero_le _) hk with
      ⟨value1, value2, hread1, _, hcoeff1, _, hlt1, _⟩
    refine ⟨value1, hread1, ?_, hlt1⟩
    simpa [coeff_karHalfSumPoly, hk] using hcoeff1
  have hprefix2 : ClassicalCoeffPrefix heap' t2 m (karHalfSumPoly right m) := by
    intro k hk
    rcases hvalues k (Nat.zero_le _) hk with
      ⟨value1, value2, _, hread2, _, hcoeff2, _, hlt2⟩
    refine ⟨value2, hread2, ?_, hlt2⟩
    simpa [coeff_karHalfSumPoly, hk] using hcoeff2
  refine ⟨slicePolyRep_of_classicalCoeffPrefix heap' t1 m
      (karHalfSumPoly left m) hT1' hprefix1 ?_,
    slicePolyRep_of_classicalCoeffPrefix heap' t2 m
      (karHalfSumPoly right m) hT2' hprefix2 ?_,
    canonicalU64Prefix_of_classicalCoeffPrefix heap' t1 m this._p
      (karHalfSumPoly left m) hprefix1,
    canonicalU64Prefix_of_classicalCoeffPrefix heap' t2 m this._p
      (karHalfSumPoly right m) hprefix2⟩
  · intro degree hdegree
    rw [coeff_karHalfSumPoly, if_neg (by omega)]
  · intro degree hdegree
    rw [coeff_karHalfSumPoly, if_neg (by omega)]

theorem karOddTail_refines_slices_odd (this : DenseUPolyZp)
    (A B t1 t2 : RawPtr UInt64) (m h : Nat) (heap heap' : RawHeap)
    (left right : Polynomial (ZMod this._p.toNat))
    (hsucc : h = m + 1)
    (hA : heap.ValidU64Slice A (m + h))
    (hB : heap.ValidU64Slice B (m + h))
    (hT1 : heap.ValidU64Slice t1 h)
    (hT2 : heap.ValidU64Slice t2 h)
    (hT1T2 : U64SlicesDisjoint t1 h t2 h)
    (hCanonicalA : CanonicalU64Prefix heap A (m + h) this._p)
    (hCanonicalB : CanonicalU64Prefix heap B (m + h) this._p)
    (hCanonicalT1 : CanonicalU64Prefix heap t1 m this._p)
    (hCanonicalT2 : CanonicalU64Prefix heap t2 m this._p)
    (hRepA : SlicePolyRep heap A (m + h) this._p.toNat left)
    (hRepB : SlicePolyRep heap B (m + h) this._p.toNat right)
    (hRepT1 : SlicePolyRep heap t1 m this._p.toNat (karHalfSumPoly left m))
    (hRepT2 : SlicePolyRep heap t2 m this._p.toNat (karHalfSumPoly right m))
    (hrun : karOddTail A B t1 t2 m h heap = .ok heap') :
    SlicePolyRep heap' t1 h this._p.toNat
        (karHalfSumPoly left m +
          Polynomial.monomial m (left.coeff (m + m))) ∧
      SlicePolyRep heap' t2 h this._p.toNat
        (karHalfSumPoly right m +
          Polynomial.monomial m (right.coeff (m + m))) ∧
      CanonicalU64Prefix heap' t1 h this._p ∧
      CanonicalU64Prefix heap' t2 h this._p := by
  have hodd : h > m := by omega
  rcases karOddTail_coeffs this A B t1 t2 m h heap heap' left right
      hA hB hT1 hT2 hT1T2 hCanonicalA hCanonicalB hRepA hRepB hodd hrun with
    ⟨aTail, bTail, ht1, ht2, hcoeffA, hcoeffB, haLt, hbLt⟩
  rcases karOddTail_preserves_own_prefixes A B t1 t2 m h heap heap'
      hA hB hT1 hT2 hT1T2 hrun with ⟨hsameT1, hsameT2⟩
  rcases karOddTail_ok A B t1 t2 m h heap hA hB hT1 hT2 with
    ⟨okHeap, hok, hlayout⟩
  have hokHeap : okHeap = heap' := Except.ok.inj (hok.symm.trans hrun)
  subst okHeap
  have hT1m := heap.validU64Slice_mono t1 h m hT1 (by omega)
  have hT2m := heap.validU64Slice_mono t2 h m hT2 (by omega)
  have hT1' := (hlayout t1 h).mp hT1
  have hT2' := (hlayout t2 h).mp hT2
  have hT1m' := (hlayout t1 m).mp hT1m
  have hT2m' := (hlayout t2 m).mp hT2m
  have hRepT1' := slicePolyRep_of_same_prefix heap heap' t1 m
    this._p.toNat (karHalfSumPoly left m) hT1m hT1m' hsameT1 hRepT1
  have hRepT2' := slicePolyRep_of_same_prefix heap heap' t2 m
    this._p.toNat (karHalfSumPoly right m) hT2m hT2m' hsameT2 hRepT2
  rcases slicePolyRep_extend_exists heap' t1 m this._p.toNat aTail
      (karHalfSumPoly left m) (by simpa [hsucc] using hT1') hRepT1' ht1 with
    ⟨full1, hfull1, heq1⟩
  rcases slicePolyRep_extend_exists heap' t2 m this._p.toNat bTail
      (karHalfSumPoly right m) (by simpa [hsucc] using hT2') hRepT2' ht2 with
    ⟨full2, hfull2, heq2⟩
  have htarget1 : full1 = karHalfSumPoly left m +
      Polynomial.monomial m (left.coeff (m + m)) := by
    simpa [hcoeffA] using heq1
  have htarget2 : full2 = karHalfSumPoly right m +
      Polynomial.monomial m (right.coeff (m + m)) := by
    simpa [hcoeffB] using heq2
  rw [htarget1] at hfull1
  rw [htarget2] at hfull2
  refine ⟨by simpa [hsucc] using hfull1, by simpa [hsucc] using hfull2, ?_, ?_⟩
  · intro k value hk hread
    by_cases hkm : k < m
    · rcases heap.readU64_of_valid t1 m k hT1m hkm with ⟨old, hold⟩
      have hpreserved := hsameT1 k old hkm hold
      have heq : value = old := Except.ok.inj (hread.symm.trans hpreserved)
      subst value
      exact hCanonicalT1 k old hkm hold
    · have hkEq : k = m := by omega
      subst k
      have heq : value = aTail := Except.ok.inj (hread.symm.trans ht1)
      subst value
      exact haLt
  · intro k value hk hread
    by_cases hkm : k < m
    · rcases heap.readU64_of_valid t2 m k hT2m hkm with ⟨old, hold⟩
      have hpreserved := hsameT2 k old hkm hold
      have heq : value = old := Except.ok.inj (hread.symm.trans hpreserved)
      subst value
      exact hCanonicalT2 k old hkm hold
    · have hkEq : k = m := by omega
      subst k
      have heq : value = bTail := Except.ok.inj (hread.symm.trans ht2)
      subst value
      exact hbLt

noncomputable def karPreparedPoly {p : Nat}
    (poly : Polynomial (ZMod p)) (m h : Nat) : Polynomial (ZMod p) :=
  if h > m then
    karHalfSumPoly poly m + Polynomial.monomial m (poly.coeff (m + m))
  else
    karHalfSumPoly poly m

theorem karPreparedPoly_eq_low_add_high {p : Nat}
    (heap : RawHeap) (ptr : RawPtr UInt64) (m h : Nat)
    (poly low high : Polynomial (ZMod p))
    (hshape : h = m ∨ h = m + 1)
    (hRepPoly : SlicePolyRep heap ptr (m + h) p poly)
    (hRepLow : SlicePolyRep heap ptr m p low)
    (hRepHigh : SlicePolyRep heap (ptr.add m) h p high)
    (hsplit : poly = low + Polynomial.X ^ m * high) :
    karPreparedPoly poly m h = low + high := by
  ext degree
  rcases hshape with heven | hodd
  · subst h
    have hnot : ¬m > m := by omega
    rw [karPreparedPoly, if_neg hnot, coeff_karHalfSumPoly]
    by_cases hd : degree < m
    · rw [if_pos hd, Polynomial.coeff_add]
      have hpolyLow : poly.coeff degree = low.coeff degree := by
        rw [hsplit, Polynomial.coeff_add, Polynomial.coeff_X_pow_mul',
          if_neg (by omega), add_zero]
      have hpolyHigh : poly.coeff (m + degree) = high.coeff degree := by
        rw [hsplit, Polynomial.coeff_add,
          slicePolyRep_coeff_zero_of_length_le heap ptr m p low hRepLow
            (m + degree) (by omega),
          Polynomial.coeff_X_pow_mul', if_pos (by omega), zero_add]
        simp

      rw [hpolyLow, hpolyHigh]
    · rw [if_neg hd, Polynomial.coeff_add,
        slicePolyRep_coeff_zero_of_length_le heap ptr m p low hRepLow degree
          (by omega),
        slicePolyRep_coeff_zero_of_length_le heap (ptr.add m) m p high
          hRepHigh degree (by omega), zero_add]
  · subst h
    have hgt : m + 1 > m := by omega
    rw [karPreparedPoly, if_pos hgt, Polynomial.coeff_add,
      Polynomial.coeff_add, coeff_karHalfSumPoly]
    by_cases hd : degree < m
    · rw [if_pos hd, Polynomial.coeff_monomial, if_neg (by omega)]
      have hpolyLow : poly.coeff degree = low.coeff degree := by
        rw [hsplit, Polynomial.coeff_add, Polynomial.coeff_X_pow_mul',
          if_neg (by omega), add_zero]
      have hpolyHigh : poly.coeff (m + degree) = high.coeff degree := by
        rw [hsplit, Polynomial.coeff_add,
          slicePolyRep_coeff_zero_of_length_le heap ptr m p low hRepLow
            (m + degree) (by omega),
          Polynomial.coeff_X_pow_mul', if_pos (by omega), zero_add]
        simp

      rw [hpolyLow, hpolyHigh, add_zero]
    · by_cases hdm : degree = m
      · subst degree
        rw [if_neg (by omega), Polynomial.coeff_monomial, if_pos rfl,
          slicePolyRep_coeff_zero_of_length_le heap ptr m p low hRepLow m
            (by omega), zero_add]
        have hpolyHigh : poly.coeff (m + m) = high.coeff m := by
          rw [hsplit, Polynomial.coeff_add,
            slicePolyRep_coeff_zero_of_length_le heap ptr m p low hRepLow
              (m + m) (by omega),
            Polynomial.coeff_X_pow_mul', if_pos (by omega), zero_add]
          simp
        simpa using hpolyHigh
      · have hdhigh : m + 1 ≤ degree := by omega
        rw [if_neg hd, Polynomial.coeff_monomial, if_neg (by omega), add_zero,
          slicePolyRep_coeff_zero_of_length_le heap ptr m p low hRepLow degree
            (by omega),
          slicePolyRep_coeff_zero_of_length_le heap (ptr.add m) (m + 1) p high
            hRepHigh degree hdhigh]
        simp

theorem karatsuba_polynomial_identity {p : Nat}
    (left right leftLow leftHigh rightLow rightHigh : Polynomial (ZMod p))
    (m : Nat)
    (hleft : left = leftLow + Polynomial.X ^ m * leftHigh)
    (hright : right = rightLow + Polynomial.X ^ m * rightHigh) :
    leftLow * rightLow +
        Polynomial.X ^ m *
          ((leftLow + leftHigh) * (rightLow + rightHigh) -
            leftLow * rightLow - leftHigh * rightHigh) +
        Polynomial.X ^ (2 * m) * (leftHigh * rightHigh) =
      left * right := by
  rw [hleft, hright]
  have hpow : Polynomial.X ^ (2 * m) =
      (Polynomial.X ^ m : Polynomial (ZMod p)) * Polynomial.X ^ m := by
    rw [two_mul, pow_add]
  rw [hpow]
  ring

theorem karPrepareHalves_refines (this : DenseUPolyZp)
    (A B t1 t2 : RawPtr UInt64) (m h : Nat) (heap heap' : RawHeap)
    (left right : Polynomial (ZMod this._p.toNat))
    (hp : 1 < this._p.toNat)
    (hshape : h = m ∨ h = m + 1)
    (hA : heap.ValidU64Slice A (m + h))
    (hB : heap.ValidU64Slice B (m + h))
    (hT1 : heap.ValidU64Slice t1 h)
    (hT2 : heap.ValidU64Slice t2 h)
    (hT1T2 : U64SlicesDisjoint t1 h t2 h)
    (hT1A : U64SlicesDisjoint t1 h A (m + h))
    (hT1B : U64SlicesDisjoint t1 h B (m + h))
    (hT2A : U64SlicesDisjoint t2 h A (m + h))
    (hT2B : U64SlicesDisjoint t2 h B (m + h))
    (hCanonicalA : CanonicalU64Prefix heap A (m + h) this._p)
    (hCanonicalB : CanonicalU64Prefix heap B (m + h) this._p)
    (hRepA : SlicePolyRep heap A (m + h) this._p.toNat left)
    (hRepB : SlicePolyRep heap B (m + h) this._p.toNat right)
    (hrun : karPrepareHalves this A B t1 t2 m h heap = .ok heap') :
    SlicePolyRep heap' t1 h this._p.toNat (karPreparedPoly left m h) ∧
      SlicePolyRep heap' t2 h this._p.toNat (karPreparedPoly right m h) ∧
      CanonicalU64Prefix heap' t1 h this._p ∧
      CanonicalU64Prefix heap' t2 h this._p := by
  have hmh : m ≤ h := by rcases hshape with h | h <;> omega
  have h2m : 2 * m ≤ m + h := by omega
  have hA2m := heap.validU64Slice_mono A (m + h) (2 * m) hA h2m
  have hB2m := heap.validU64Slice_mono B (m + h) (2 * m) hB h2m
  have hT1m := heap.validU64Slice_mono t1 h m hT1 hmh
  have hT2m := heap.validU64Slice_mono t2 h m hT2 hmh
  have hT1T2m := u64SlicesDisjoint_mono hT1T2 hmh hmh
  have hT1A2m := u64SlicesDisjoint_mono hT1A hmh h2m
  have hT1B2m := u64SlicesDisjoint_mono hT1B hmh h2m
  have hT2A2m := u64SlicesDisjoint_mono hT2A hmh h2m
  have hT2B2m := u64SlicesDisjoint_mono hT2B hmh h2m
  rcases slicePolyRep_prefix_exists heap A (m + h) (2 * m)
      this._p.toNat left hA h2m hRepA with
    ⟨leftPrefix, hRepAPrefix, hleftPrefix⟩
  rcases slicePolyRep_prefix_exists heap B (m + h) (2 * m)
      this._p.toNat right hB h2m hRepB with
    ⟨rightPrefix, hRepBPrefix, hrightPrefix⟩
  rcases karAddHalvesLoop_ok this A B t1 t2 m 0 heap hA2m hB2m hT1m
      hT2m with ⟨heap1, hadd, hlayout1⟩
  have htail : karOddTail A B t1 t2 m h heap1 = .ok heap' := by
    simpa [karPrepareHalves, hadd] using hrun
  have hCanonicalAPrefix : CanonicalU64Prefix heap A (2 * m) this._p := by
    intro i value hi hread
    exact hCanonicalA i value (by omega) hread
  have hCanonicalBPrefix : CanonicalU64Prefix heap B (2 * m) this._p := by
    intro i value hi hread
    exact hCanonicalB i value (by omega) hread
  rcases karAddHalvesLoop_refines_slices this A B t1 t2 m heap heap1
      leftPrefix rightPrefix hp hA2m hB2m hT1m hT2m hT1T2m
      hT1A2m hT1B2m hT2A2m hT2B2m hCanonicalAPrefix hCanonicalBPrefix
      hRepAPrefix hRepBPrefix hadd with
    ⟨hRepT1Prefix, hRepT2Prefix, hCanonicalT1, hCanonicalT2⟩
  have hsumLeft : karHalfSumPoly leftPrefix m = karHalfSumPoly left m :=
    karHalfSumPoly_congr_of_prefix leftPrefix left m hleftPrefix
  have hsumRight : karHalfSumPoly rightPrefix m = karHalfSumPoly right m :=
    karHalfSumPoly_congr_of_prefix rightPrefix right m hrightPrefix
  rw [hsumLeft] at hRepT1Prefix
  rw [hsumRight] at hRepT2Prefix
  have hsameA : SameU64Prefix heap heap1 A (m + h) := by
    intro i old hi hread
    apply karAddHalvesLoop_preserves_outside this A B t1 t2 A m 0 i
      heap heap1 old hA2m hB2m hT1m hT2m hread
    · intro j _ hj
      exact hT1A j (by omega) i hi
    · intro j _ hj
      exact hT2A j (by omega) i hi
    · exact hadd
  have hsameB : SameU64Prefix heap heap1 B (m + h) := by
    intro i old hi hread
    apply karAddHalvesLoop_preserves_outside this A B t1 t2 B m 0 i
      heap heap1 old hA2m hB2m hT1m hT2m hread
    · intro j _ hj
      exact hT1B j (by omega) i hi
    · intro j _ hj
      exact hT2B j (by omega) i hi
    · exact hadd
  have hA1 := (hlayout1 A (m + h)).mp hA
  have hB1 := (hlayout1 B (m + h)).mp hB
  have hT11 := (hlayout1 t1 h).mp hT1
  have hT21 := (hlayout1 t2 h).mp hT2
  have hRepA1 := slicePolyRep_of_same_prefix heap heap1 A (m + h)
    this._p.toNat left hA hA1 hsameA hRepA
  have hRepB1 := slicePolyRep_of_same_prefix heap heap1 B (m + h)
    this._p.toNat right hB hB1 hsameB hRepB
  have hCanonicalA1 : CanonicalU64Prefix heap1 A (m + h) this._p := by
    intro i value hi hread1
    rcases heap.readU64_of_valid A (m + h) i hA hi with ⟨old, hread⟩
    have hpreserved := hsameA i old hi hread
    have heq : value = old := Except.ok.inj (hread1.symm.trans hpreserved)
    subst value
    exact hCanonicalA i old hi hread
  have hCanonicalB1 : CanonicalU64Prefix heap1 B (m + h) this._p := by
    intro i value hi hread1
    rcases heap.readU64_of_valid B (m + h) i hB hi with ⟨old, hread⟩
    have hpreserved := hsameB i old hi hread
    have heq : value = old := Except.ok.inj (hread1.symm.trans hpreserved)
    subst value
    exact hCanonicalB i old hi hread
  rcases hshape with heven | hodd
  · have hnot : ¬h > m := by omega
    unfold karOddTail at htail
    simp only [hnot, ↓reduceIte] at htail
    have heq : heap' = heap1 := Except.ok.inj htail.symm
    subst heap'
    exact ⟨by simpa [karPreparedPoly, hnot, heven] using hRepT1Prefix,
      by simpa [karPreparedPoly, hnot, heven] using hRepT2Prefix,
      by simpa [heven] using hCanonicalT1,
      by simpa [heven] using hCanonicalT2⟩
  · have hoddGt : h > m := by omega
    rcases karOddTail_refines_slices_odd this A B t1 t2 m h heap1 heap'
      left right hodd hA1 hB1 hT11 hT21 hT1T2 hCanonicalA1
      hCanonicalB1 hCanonicalT1 hCanonicalT2 hRepA1 hRepB1
      hRepT1Prefix hRepT2Prefix htail with
      ⟨hT1Final, hT2Final, hCanonicalT1Final, hCanonicalT2Final⟩
    exact ⟨by simpa [karPreparedPoly, hoddGt] using hT1Final,
      by simpa [karPreparedPoly, hoddGt] using hT2Final,
      hCanonicalT1Final, hCanonicalT2Final⟩

/-- Full interface needed by the first recursive Karatsuba child: the actual
half-preparation execution produces the two sum slices while preserving both
input slices in the resulting heap. -/
theorem karPrepareHalves_refines_and_preserves_inputs (this : DenseUPolyZp)
    (A B t1 t2 : RawPtr UInt64) (m h : Nat) (heap heap' : RawHeap)
    (left right : Polynomial (ZMod this._p.toNat))
    (hp : 1 < this._p.toNat)
    (hshape : h = m ∨ h = m + 1)
    (hA : heap.ValidU64Slice A (m + h))
    (hB : heap.ValidU64Slice B (m + h))
    (hT1 : heap.ValidU64Slice t1 h)
    (hT2 : heap.ValidU64Slice t2 h)
    (hT1T2 : U64SlicesDisjoint t1 h t2 h)
    (hT1A : U64SlicesDisjoint t1 h A (m + h))
    (hT1B : U64SlicesDisjoint t1 h B (m + h))
    (hT2A : U64SlicesDisjoint t2 h A (m + h))
    (hT2B : U64SlicesDisjoint t2 h B (m + h))
    (hCanonicalA : CanonicalU64Prefix heap A (m + h) this._p)
    (hCanonicalB : CanonicalU64Prefix heap B (m + h) this._p)
    (hRepA : SlicePolyRep heap A (m + h) this._p.toNat left)
    (hRepB : SlicePolyRep heap B (m + h) this._p.toNat right)
    (hrun : karPrepareHalves this A B t1 t2 m h heap = .ok heap') :
    RawHeap.SameLayout heap heap' ∧
      SlicePolyRep heap' A (m + h) this._p.toNat left ∧
      SlicePolyRep heap' B (m + h) this._p.toNat right ∧
      CanonicalU64Prefix heap' A (m + h) this._p ∧
      CanonicalU64Prefix heap' B (m + h) this._p ∧
      SlicePolyRep heap' t1 h this._p.toNat
        (karPreparedPoly left m h) ∧
      SlicePolyRep heap' t2 h this._p.toNat
        (karPreparedPoly right m h) ∧
      CanonicalU64Prefix heap' t1 h this._p ∧
      CanonicalU64Prefix heap' t2 h this._p := by
  have hmh : m ≤ h := by rcases hshape with h | h <;> omega
  rcases karPrepareHalves_ok this A B t1 t2 m h heap hmh hA hB hT1 hT2 with
    ⟨okHeap, hok, hlayout⟩
  have heq : okHeap = heap' := Except.ok.inj (hok.symm.trans hrun)
  subst okHeap
  rcases karPrepareHalves_refines this A B t1 t2 m h heap heap' left right
      hp hshape hA hB hT1 hT2 hT1T2 hT1A hT1B hT2A hT2B
      hCanonicalA hCanonicalB hRepA hRepB hrun with
    ⟨hRepT1, hRepT2, hCanonicalT1, hCanonicalT2⟩
  have hsameA := karPrepareHalves_preserves_outside this A B t1 t2 A m h
    (m + h) heap heap' hmh hA hB hT1 hT2 hT1A hT2A hrun
  have hsameB := karPrepareHalves_preserves_outside this A B t1 t2 B m h
    (m + h) heap heap' hmh hA hB hT1 hT2 hT1B hT2B hrun
  have hA' := (hlayout A (m + h)).mp hA
  have hB' := (hlayout B (m + h)).mp hB
  have hRepA' := slicePolyRep_of_same_prefix heap heap' A (m + h)
    this._p.toNat left hA hA' hsameA hRepA
  have hRepB' := slicePolyRep_of_same_prefix heap heap' B (m + h)
    this._p.toNat right hB hB' hsameB hRepB
  have hCanonicalA' : CanonicalU64Prefix heap' A (m + h) this._p := by
    intro k value hk hread'
    rcases heap.readU64_of_valid A (m + h) k hA hk with ⟨old, hread⟩
    have hvalue : value = old :=
      Except.ok.inj (hread'.symm.trans (hsameA k old hk hread))
    subst value
    exact hCanonicalA k old hk hread
  have hCanonicalB' : CanonicalU64Prefix heap' B (m + h) this._p := by
    intro k value hk hread'
    rcases heap.readU64_of_valid B (m + h) k hB hk with ⟨old, hread⟩
    have hvalue : value = old :=
      Except.ok.inj (hread'.symm.trans (hsameB k old hk hread))
    subst value
    exact hCanonicalB k old hk hread
  exact ⟨hlayout, hRepA', hRepB', hCanonicalA', hCanonicalB', hRepT1,
    hRepT2, hCanonicalT1, hCanonicalT2⟩

theorem normaliseU64_eq_length_of_classicalCoeffPrefix {p : Nat}
    (heap : RawHeap) (C : RawPtr UInt64) (length : Nat)
    (poly : Polynomial (ZMod p))
    (hpos : 0 < length)
    (hprefix : ClassicalCoeffPrefix heap C length poly)
    (hlast : poly.coeff (length - 1) ≠ 0) :
    heap.normaliseU64 C length = .ok length := by
  cases length with
  | zero => omega
  | succ n =>
    rcases hprefix n (by omega) with ⟨value, hread, hcoeff, _⟩
    have hvalue : value ≠ 0 := by
      intro hzero
      subst value
      apply hlast
      simpa using hcoeff.symm
    simp [RawHeap.normaliseU64, hread, hvalue]

theorem mul_coeff_zero_of_slice_lengths {p : Nat}
    (heap : RawHeap) (A B : RawPtr UInt64) (lenA lenB degree : Nat)
    (left right : Polynomial (ZMod p))
    (hApos : 0 < lenA) (hBpos : 0 < lenB)
    (hRepA : SlicePolyRep heap A lenA p left)
    (hRepB : SlicePolyRep heap B lenB p right)
    (hdegree : lenA + lenB - 1 ≤ degree) :
    (left * right).coeff degree = 0 := by
  rw [Polynomial.coeff_mul]
  apply Finset.sum_eq_zero
  intro pair hpair
  simp only [Finset.mem_antidiagonal] at hpair
  by_cases hleft : lenA ≤ pair.1
  · rw [slicePolyRep_coeff_zero_of_length_le heap A lenA p left hRepA
      pair.1 hleft, zero_mul]
  · have hright : lenB ≤ pair.2 := by omega
    rw [slicePolyRep_coeff_zero_of_length_le heap B lenB p right hRepB
      pair.2 hright, mul_zero]

theorem mul_last_coeff_ne_zero_of_rawDense (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (heap : RawHeap) (A B : RawPtr UInt64) (lenA lenB : Nat)
    (left right : Polynomial (ZMod this._p.toNat))
    (hApos : 0 < lenA) (hBpos : 0 < lenB)
    (hLeft : RawDensePolyRep this heap A lenA left)
    (hRight : RawDensePolyRep this heap B lenB right) :
    (left * right).coeff (lenA + lenB - 1 - 1) ≠ 0 := by
  rcases hLeft with ⟨hA, hCanonicalA, hRepA, hNormA⟩
  rcases hRight with ⟨hB, hCanonicalB, hRepB, hNormB⟩
  have hleftLast := normaliseU64_poly_last_coeff_ne_zero heap A lenA
    this._p.toNat lenA left hA hRepA hCanonicalA hNormA
    (Nat.ne_of_gt hApos)
  have hrightLast := normaliseU64_poly_last_coeff_ne_zero heap B lenB
    this._p.toNat lenB right hB hRepB hCanonicalB hNormB
    (Nat.ne_of_gt hBpos)
  have hleftDegree := normaliseU64_poly_natDegree_eq heap A lenA
    this._p.toNat lenA left hA hRepA hCanonicalA hNormA
    (Nat.ne_of_gt hApos)
  have hrightDegree := normaliseU64_poly_natDegree_eq heap B lenB
    this._p.toNat lenB right hB hRepB hCanonicalB hNormB
    (Nat.ne_of_gt hBpos)
  have hcoeff := Polynomial.coeff_mul_add_eq_of_natDegree_le
    (show left.natDegree ≤ lenA - 1 by omega)
    (show right.natDegree ≤ lenB - 1 by omega)
  rw [show lenA + lenB - 1 - 1 = (lenA - 1) + (lenB - 1) by omega,
    hcoeff]
  exact mul_ne_zero hleftLast hrightLast

theorem classicalCoeffPrefix_succ_of_write {p : Nat}
    (before after : RawHeap) (C : RawPtr UInt64) (upto : Nat)
    (poly : Polynomial (ZMod p)) (value : UInt64)
    (hprefix : ClassicalCoeffPrefix before C upto poly)
    (hwrite : before.writeU64 C upto value = .ok after)
    (hvalue : (value.toNat : ZMod p) = poly.coeff upto)
    (hvalueLt : value.toNat < p) :
    ClassicalCoeffPrefix after C (upto + 1) poly := by
  intro i hi
  by_cases heq : i = upto
  · subst i
    exact ⟨value, RawHeap.readU64_writeU64_same before after C upto value
      hwrite, hvalue, hvalueLt⟩
  · have hiOld : i < upto := by omega
    rcases hprefix i hiOld with ⟨old, hread, hold, holdLt⟩
    refine ⟨old, ?_, hold, holdLt⟩
    exact RawHeap.readU64_writeU64_ne before after C C upto i value old
      hwrite hread (Or.inr (by omega))

theorem classicalOuterLoop_preserves_coeff_prefix {p : Nat}
    (this : DenseUPolyZp) (C A B : RawPtr UInt64)
    (lenA lenB lenC k : Nat) (heap heap' : RawHeap)
    (poly : Polynomial (ZMod p))
    (hC : heap.ValidU64Slice C lenC)
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB)
    (hprefix : ClassicalCoeffPrefix heap C k poly)
    (hrun : classicalOuterLoop this C A B lenA lenB lenC k heap =
      .ok heap') :
    ClassicalCoeffPrefix heap' C k poly := by
  intro i hi
  rcases hprefix i hi with ⟨old, hread, hold, holdLt⟩
  refine ⟨old, ?_, hold, holdLt⟩
  apply classicalOuterLoop_preserves_outside this C A B C lenA lenB lenC
    k i heap heap' old hC hA hB hread
  · intro target hktarget _
    exact Or.inr (by omega)
  · exact hrun

theorem classicalOuterLoop_refines_coeff_prefix (this : DenseUPolyZp)
    (C A B : RawPtr UInt64) (lenA lenB lenC k : Nat) (heap : RawHeap)
    (left right : Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this)
    (hp : 1 < this._p.toNat)
    (hApos : 0 < lenA) (hBpos : 0 < lenB)
    (hLenAWord : lenA < limbBase)
    (hlenC : lenC = lenA + lenB - 1)
    (hC : heap.ValidU64Slice C lenC)
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB)
    (hCA : U64SlicesDisjoint C lenC A lenA)
    (hCB : U64SlicesDisjoint C lenC B lenB)
    (hCanonicalA : CanonicalU64Prefix heap A lenA this._p)
    (hCanonicalB : CanonicalU64Prefix heap B lenB this._p)
    (hRepA : SlicePolyRep heap A lenA this._p.toNat left)
    (hRepB : SlicePolyRep heap B lenB this._p.toNat right)
    (hprefix : ClassicalCoeffPrefix heap C k (left * right)) :
    ∃ heap', classicalOuterLoop this C A B lenA lenB lenC k heap =
        .ok heap' ∧ RawHeap.SameLayout heap heap' ∧
      ClassicalCoeffPrefix heap' C lenC (left * right) := by
  unfold classicalOuterLoop
  split
  next hk =>
    let jMin := if k ≥ lenB then k - lenB + 1 else 0
    let jMax := if k < lenA then k else lenA - 1
    rcases classicalDotLoop_ok heap A B lenA lenB k jMax jMin
      { lo := 0, mid := 0, hi := 0 } hA hB
      (by
        intro t hjt hts
        exact (classical_index_bounds lenA lenB k t hApos hBpos
          (by omega) hjt hts).1)
      (by
        intro t hjt hts
        exact (classical_index_bounds lenA lenB k t hApos hBpos
          (by omega) hjt hts).2) with ⟨acc, hdot⟩
    simp only [jMin, jMax, hdot]
    let value := Generated.StrictGCD.dense_upoly_zp__lll_mod_preinv_ir
      acc.hi acc.mid acc.lo this._p this._ninv this._norm
    have hvalueResult : (value.toNat : ZMod this._p.toNat) =
          (left * right).coeff k ∧ value.toNat < this._p.toNat := by
      dsimp [value]
      exact classicalReduced_source_eq_coeff this heap A B lenA lenB lenC k
        left right acc hcfg hp hApos hBpos hLenAWord hlenC hk hA hB
        hCanonicalA hCanonicalB hRepA hRepB (by simpa [jMin, jMax] using hdot)
    have hvalue := hvalueResult.1
    have hvalueLt := hvalueResult.2
    rcases heap.writeU64_of_valid C lenC k value hC hk with ⟨heap1, hw⟩
    simp only [value, hw]
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 C k value hw
    have hC1 := (hlayout1 C lenC).mp hC
    have hA1 := (hlayout1 A lenA).mp hA
    have hB1 := (hlayout1 B lenB).mp hB
    have hsameA : SameU64Prefix heap heap1 A lenA := by
      intro i old hi hread
      exact RawHeap.readU64_writeU64_ne heap heap1 C A k i value old hw
        hread (hCA k hk i hi)
    have hsameB : SameU64Prefix heap heap1 B lenB := by
      intro i old hi hread
      exact RawHeap.readU64_writeU64_ne heap heap1 C B k i value old hw
        hread (hCB k hk i hi)
    have hCanonicalA1 := canonicalU64Prefix_of_same_prefix heap heap1 A lenA
      this._p hA hsameA hCanonicalA
    have hCanonicalB1 := canonicalU64Prefix_of_same_prefix heap heap1 B lenB
      this._p hB hsameB hCanonicalB
    have hRepA1 := slicePolyRep_of_same_prefix heap heap1 A lenA
      this._p.toNat left hA hA1 hsameA hRepA
    have hRepB1 := slicePolyRep_of_same_prefix heap heap1 B lenB
      this._p.toNat right hB hB1 hsameB hRepB
    have hprefix1 := classicalCoeffPrefix_succ_of_write heap heap1 C k
      (left * right) value hprefix hw hvalue hvalueLt
    rcases classicalOuterLoop_refines_coeff_prefix this C A B lenA lenB lenC
      (k + 1) heap1 left right hcfg hp hApos hBpos hLenAWord hlenC hC1
      hA1 hB1 hCA hCB hCanonicalA1 hCanonicalB1 hRepA1 hRepB1 hprefix1 with
      ⟨heap2, hrun, hlayout2, hfull⟩
    rw [hrun]
    exact ⟨heap2, rfl, fun ptr length =>
      (hlayout1 ptr length).trans (hlayout2 ptr length), hfull⟩
  next hnot =>
    refine ⟨heap, rfl, fun _ _ => Iff.rfl, ?_⟩
    intro i hi
    exact hprefix i (by omega)
termination_by lenC - k
decreasing_by omega

theorem classicalOuterLoop_ok (this : DenseUPolyZp)
    (C A B : RawPtr UInt64) (lenA lenB lenC k : Nat) (heap : RawHeap)
    (hApos : 0 < lenA) (hBpos : 0 < lenB)
    (hlenC : lenC = lenA + lenB - 1)
    (hC : heap.ValidU64Slice C lenC)
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB) :
    ∃ heap', classicalOuterLoop this C A B lenA lenB lenC k heap =
        .ok heap' ∧ RawHeap.SameLayout heap heap' := by
  unfold classicalOuterLoop
  split
  next hk =>
    let jMin := if k ≥ lenB then k - lenB + 1 else 0
    let jMax := if k < lenA then k else lenA - 1
    have hrange : jMin ≤ jMax := by
      dsimp [jMin, jMax]
      split <;> split <;> omega
    rcases classicalDotLoop_ok heap A B lenA lenB k jMax jMin
      { lo := 0, mid := 0, hi := 0 } hA hB
      (by
        intro t hjt hts
        exact (classical_index_bounds lenA lenB k t hApos hBpos
          (by omega) hjt hts).1)
      (by
        intro t hjt hts
        exact (classical_index_bounds lenA lenB k t hApos hBpos
          (by omega) hjt hts).2) with ⟨acc, hdot⟩
    simp only [jMin, jMax, hdot]
    let value := Generated.StrictGCD.dense_upoly_zp__lll_mod_preinv_ir
      acc.hi acc.mid acc.lo this._p this._ninv this._norm
    rcases heap.writeU64_of_valid C lenC k value hC hk with ⟨heap1, hw⟩
    simp only [value, hw]
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 C k value hw
    rcases classicalOuterLoop_ok this C A B lenA lenB lenC (k + 1) heap1
      hApos hBpos hlenC ((hlayout1 C lenC).mp hC)
      ((hlayout1 A lenA).mp hA) ((hlayout1 B lenB).mp hB) with
      ⟨heap2, hrun, hlayout2⟩
    rw [hrun]
    exact ⟨heap2, rfl, fun ptr length =>
      (hlayout1 ptr length).trans (hlayout2 ptr length)⟩
  next => exact ⟨heap, rfl, fun _ _ => Iff.rfl⟩
termination_by lenC - k
decreasing_by omega

theorem classicalMul_ok (this : DenseUPolyZp)
    (C A : RawPtr UInt64) (lenA : Nat) (B : RawPtr UInt64) (lenB : Nat)
    (heap : RawHeap) (hApos : 0 < lenA) (hBpos : 0 < lenB)
    (hC : heap.ValidU64Slice C (lenA + lenB - 1))
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB) :
    ∃ heap', dense_upoly_zp__classical_mul_ir this C A lenA B lenB heap =
        .ok heap' ∧ RawHeap.SameLayout heap heap' := by
  rcases classicalOuterLoop_ok this C A B lenA lenB
    (lenA + lenB - 1) 0 heap hApos hBpos rfl hC hA hB with
    ⟨heap', hrun, hlayout⟩
  refine ⟨heap', ?_, hlayout⟩
  simp [dense_upoly_zp__classical_mul_ir, Nat.ne_of_gt hApos,
    Nat.ne_of_gt hBpos, hrun]

theorem classicalMul_preserves_outside (this : DenseUPolyZp)
    (C A : RawPtr UInt64) (lenA : Nat) (B : RawPtr UInt64) (lenB : Nat)
    (guard : RawPtr UInt64) (guardLen : Nat) (heap heap' : RawHeap)
    (hApos : 0 < lenA) (hBpos : 0 < lenB)
    (hC : heap.ValidU64Slice C (lenA + lenB - 1))
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB)
    (hCGuard : U64SlicesDisjoint C (lenA + lenB - 1) guard guardLen)
    (hrun : dense_upoly_zp__classical_mul_ir this C A lenA B lenB heap =
      .ok heap') :
    SameU64Prefix heap heap' guard guardLen := by
  have houter : classicalOuterLoop this C A B lenA lenB
      (lenA + lenB - 1) 0 heap = .ok heap' := by
    simpa [dense_upoly_zp__classical_mul_ir, Nat.ne_of_gt hApos,
      Nat.ne_of_gt hBpos] using hrun
  intro readIndex old hreadIndex hread
  exact classicalOuterLoop_preserves_outside this C A B guard lenA lenB
    (lenA + lenB - 1) 0 readIndex heap heap' old hC hA hB hread
    (by
      intro writeIndex _ hwriteIndex
      exact hCGuard writeIndex hwriteIndex readIndex hreadIndex)
    houter

/-- The generated Karatsuba routine terminates on its exact raw slices.
This follows the three C++ recursive calls and uses `n` itself as the
well-founded measure with no auxiliary recursion counter. -/
theorem karMul_ok (this : DenseUPolyZp)
    (C A B : RawPtr UInt64) (n : Nat) (scratch : RawPtr UInt64)
    (heap : RawHeap) (hn : 0 < n)
    (hC : heap.ValidU64Slice C (2 * n - 1))
    (hA : heap.ValidU64Slice A n)
    (hB : heap.ValidU64Slice B n)
    (hScratch : heap.ValidU64Slice scratch (karScratchNeed n)) :
    ∃ heap', dense_upoly_zp__kar_mul_ir this C A B n scratch heap =
        .ok heap' ∧ RawHeap.SameLayout heap heap' ∧
      ∀ guard guardLength,
        U64SlicesDisjoint C (2 * n - 1) guard guardLength →
        U64SlicesDisjoint scratch (karScratchNeed n) guard guardLength →
        SameU64Prefix heap heap' guard guardLength := by
  unfold dense_upoly_zp__kar_mul_ir
  split
  next hbase =>
    rcases classicalMul_ok this C A n B n heap hn hn
      (by simpa [two_mul] using hC) hA hB with ⟨heap', hrun, hlayout⟩
    refine ⟨heap', hrun, hlayout, ?_⟩
    intro guard guardLength hCGuard _
    exact classicalMul_preserves_outside this C A n B n guard guardLength
      heap heap' hn hn (by simpa [two_mul] using hC) hA hB
      (by simpa [two_mul] using hCGuard) hrun
  next hrecCase =>
    let m := n / 2
    let h := n - m
    let t1 := scratch
    let t2 := t1.add h
    let sP0 := t2.add h
    let sP1 := sP0.add (2 * m - 1)
    let recScratch := sP1.add (2 * h - 1)
    have hn16 : 16 ≤ n := by omega
    have hmPos : 0 < m := by
      dsimp [m]
      exact Nat.div_pos (by omega) (by omega)
    have hhPos : 0 < h := by
      dsimp [h, m]
      omega
    have hmLt : m < n := by
      simpa [m] using (kar_split_children_lt n hn16).1
    have hhLt : h < n := by
      simpa [h, m] using (kar_split_children_lt n hn16).2
    have hmh : m ≤ h := by
      rcases kar_split_shape n with heven | hodd
      · simpa [m, h] using le_of_eq heven.symm
      · dsimp [m, h]
        omega
    rcases karScratchSlices heap scratch n hn16 hScratch with
      ⟨hT1, hT2, hP0, hP1, hRec⟩
    have hAm := heap.validU64Slice_mono A n m hA (by omega)
    have hBm := heap.validU64Slice_mono B n m hB (by omega)
    have hAh := heap.validU64Slice_add A n m h hA (by
      dsimp [m, h]
      omega)
    have hBh := heap.validU64Slice_add B n m h hB (by
      dsimp [m, h]
      omega)
    have hCHigh := heap.validU64Slice_add C (2 * n - 1) (2 * m)
      (2 * h - 1) hC (by
        dsimp [m, h]
        omega)
    have hRecM := heap.validU64Slice_mono recScratch
      (max (karScratchNeed m) (karScratchNeed h)) (karScratchNeed m)
      hRec (Nat.le_max_left _ _)
    have hRecH := heap.validU64Slice_mono recScratch
      (max (karScratchNeed m) (karScratchNeed h)) (karScratchNeed h)
      hRec (Nat.le_max_right _ _)
    have hmn : m + h = n := by
      dsimp [m, h]
      omega
    have hAfull : heap.ValidU64Slice A (m + h) := by
      rw [hmn]
      exact hA
    have hBfull : heap.ValidU64Slice B (m + h) := by
      rw [hmn]
      exact hB
    rcases karPrepareHalves_ok this A B t1 t2 m h heap hmh
      hAfull hBfull hT1 hT2 with
      ⟨heap2, hprep, hlay2⟩
    rcases karMul_ok this sP0 A B m recScratch heap2 hmPos
      ((hlay2 sP0 (2 * m - 1)).mp hP0)
      ((hlay2 A m).mp hAm) ((hlay2 B m).mp hBm)
      ((hlay2 recScratch (karScratchNeed m)).mp hRecM) with
      ⟨heap3, hp0, hlay3, hframe3⟩
    rcases karMul_ok this sP1 t1 t2 h recScratch heap3 hhPos
      ((hlay3 sP1 (2 * h - 1)).mp ((hlay2 sP1 (2 * h - 1)).mp hP1))
      ((hlay3 t1 h).mp ((hlay2 t1 h).mp hT1))
      ((hlay3 t2 h).mp ((hlay2 t2 h).mp hT2))
      ((hlay3 recScratch (karScratchNeed h)).mp
        ((hlay2 recScratch (karScratchNeed h)).mp hRecH)) with
      ⟨heap4, hp1, hlay4, hframe4⟩
    rcases karMul_ok this (C.add (2 * m)) (A.add m) (B.add m) h
      recScratch heap4 hhPos
      ((hlay4 (C.add (2 * m)) (2 * h - 1)).mp
        ((hlay3 (C.add (2 * m)) (2 * h - 1)).mp
          ((hlay2 (C.add (2 * m)) (2 * h - 1)).mp hCHigh)))
      ((hlay4 (A.add m) h).mp ((hlay3 (A.add m) h).mp
        ((hlay2 (A.add m) h).mp hAh)))
      ((hlay4 (B.add m) h).mp ((hlay3 (B.add m) h).mp
        ((hlay2 (B.add m) h).mp hBh)))
      ((hlay4 recScratch (karScratchNeed h)).mp
        ((hlay3 recScratch (karScratchNeed h)).mp
          ((hlay2 recScratch (karScratchNeed h)).mp hRecH))) with
      ⟨heap5, phigh, hlay5, hframe5⟩
    have hP1short := heap5.validU64Slice_mono sP1 (2 * h - 1)
      (2 * m - 1)
      ((hlay5 sP1 (2 * h - 1)).mp ((hlay4 sP1 (2 * h - 1)).mp
        ((hlay3 sP1 (2 * h - 1)).mp ((hlay2 sP1 (2 * h - 1)).mp hP1))))
      (by omega)
    have hP05 := (hlay5 sP0 (2 * m - 1)).mp
      ((hlay4 sP0 (2 * m - 1)).mp ((hlay3 sP0 (2 * m - 1)).mp
        ((hlay2 sP0 (2 * m - 1)).mp hP0)))
    rcases karSubLoop_ok this sP1 sP0 (2 * m - 1) 0 heap5 hP1short hP05 with
      ⟨heap6, hsub0, hlay6⟩
    have hP16 := (hlay6 sP1 (2 * h - 1)).mp
      ((hlay5 sP1 (2 * h - 1)).mp ((hlay4 sP1 (2 * h - 1)).mp
        ((hlay3 sP1 (2 * h - 1)).mp ((hlay2 sP1 (2 * h - 1)).mp hP1))))
    have hCHigh6 := (hlay6 (C.add (2 * m)) (2 * h - 1)).mp
      ((hlay5 (C.add (2 * m)) (2 * h - 1)).mp
        ((hlay4 (C.add (2 * m)) (2 * h - 1)).mp
          ((hlay3 (C.add (2 * m)) (2 * h - 1)).mp
            ((hlay2 (C.add (2 * m)) (2 * h - 1)).mp hCHigh))))
    rcases karSubLoop_ok this sP1 (C.add (2 * m)) (2 * h - 1) 0 heap6
      hP16 hCHigh6 with ⟨heap7, hsub1, hlay7⟩
    have hC7 := (hlay7 C (2 * n - 1)).mp ((hlay6 C (2 * n - 1)).mp
      ((hlay5 C (2 * n - 1)).mp ((hlay4 C (2 * n - 1)).mp
        ((hlay3 C (2 * n - 1)).mp ((hlay2 C (2 * n - 1)).mp hC)))))
    have hCprefix := heap7.validU64Slice_mono C (2 * n - 1)
      (2 * m - 1) hC7 (by omega)
    have hP07 := (hlay7 sP0 (2 * m - 1)).mp ((hlay6 sP0 (2 * m - 1)).mp
      ((hlay5 sP0 (2 * m - 1)).mp ((hlay4 sP0 (2 * m - 1)).mp
        ((hlay3 sP0 (2 * m - 1)).mp ((hlay2 sP0 (2 * m - 1)).mp hP0)))))
    rcases copyU64_ok heap7 C sP0 (2 * m - 1) hCprefix hP07 with
      ⟨heap8, hcopy, hlay8⟩
    rcases heap8.writeU64_of_valid C (2 * n - 1) (2 * m - 1) 0
      ((hlay8 C (2 * n - 1)).mp hC7) (by omega) with
      ⟨heap9, hzero⟩
    have hlay9 := RawHeap.writeU64_sameLayout heap8 heap9 C (2 * m - 1) 0
      hzero
    have hCassemble := heap9.validU64Slice_mono C (2 * n - 1)
      (m + (2 * h - 1)) ((hlay9 C (2 * n - 1)).mp
        ((hlay8 C (2 * n - 1)).mp hC7)) (by
          dsimp [m, h]
          omega)
    have hP19 := (hlay9 sP1 (2 * h - 1)).mp
      ((hlay8 sP1 (2 * h - 1)).mp ((hlay7 sP1 (2 * h - 1)).mp hP16))
    rcases karAssembleLoop_ok this C sP1 m (2 * h - 1) 0 heap9
      hCassemble hP19 with ⟨heap10, hassemble, hlay10⟩
    refine ⟨heap10, ?_, ?_, ?_⟩
    · simp [m, h, t1, t2, sP0, sP1, recScratch, hprep, hp0, hp1,
        phigh, hsub0, hsub1, hcopy, hzero, hassemble]
    · intro ptr length
      exact (hlay2 ptr length).trans ((hlay3 ptr length).trans
        ((hlay4 ptr length).trans ((hlay5 ptr length).trans
          ((hlay6 ptr length).trans ((hlay7 ptr length).trans
            ((hlay8 ptr length).trans ((hlay9 ptr length).trans
              (hlay10 ptr length))))))))
    · intro guard guardLength hCGuard hScratchGuard
      rcases karScratchSlices_disjoint_guard scratch guard n guardLength hn16
          hScratchGuard with ⟨hT1Guard, hT2Guard, hP0Guard, hP1Guard,
            hRecGuard⟩
      have hRecMGuard := u64SlicesDisjoint_mono hRecGuard
        (Nat.le_max_left _ _) (Nat.le_refl guardLength)
      have hRecHGuard := u64SlicesDisjoint_mono hRecGuard
        (Nat.le_max_right _ _) (Nat.le_refl guardLength)
      have hCHighGuard := u64SlicesDisjoint_add_left hCGuard
        (start := 2 * m) (count := 2 * h - 1) (by
          dsimp [m, h]
          omega)
      have hPrepFrame := karPrepareHalves_preserves_outside this A B t1 t2
        guard m h guardLength heap heap2 hmh hAfull hBfull hT1 hT2
        hT1Guard hT2Guard hprep
      have hP0Frame := hframe3 guard guardLength hP0Guard hRecMGuard
      have hP1Frame := hframe4 guard guardLength hP1Guard hRecHGuard
      have hHighFrame := hframe5 guard guardLength hCHighGuard hRecHGuard
      have hP1ShortGuard := u64SlicesDisjoint_mono hP1Guard
        (smallLeft := 2 * m - 1) (by omega)
        (Nat.le_refl guardLength)
      have hSub0Frame := karSubLoop_preserves_prefix this sP1 sP0 guard
        (2 * m - 1) guardLength heap5 heap6 hP1short hP05
        hP1ShortGuard hsub0
      have hSub1Frame := karSubLoop_preserves_prefix this sP1
        (C.add (2 * m)) guard (2 * h - 1) guardLength heap6 heap7
        hP16 hCHigh6 hP1Guard hsub1
      have hCPrefixGuard := u64SlicesDisjoint_mono hCGuard
        (smallLeft := 2 * m - 1) (by omega)
        (Nat.le_refl guardLength)
      have hCopyFrame := copyU64_preserves_prefix heap7 heap8 C sP0 guard
        (2 * m - 1) guardLength hCprefix hP07 hCPrefixGuard hcopy
      have hZeroFrame := writeU64_preserves_prefix heap8 heap9 C guard
        (2 * n - 1) guardLength (2 * m - 1) 0 hCGuard (by omega) hzero
      have hCAssembleGuard := u64SlicesDisjoint_mono hCGuard
        (smallLeft := m + (2 * h - 1)) (by
          dsimp [m, h]
          omega) (Nat.le_refl guardLength)
      have hAssembleFrame := karAssembleLoop_preserves_prefix this C sP1
        guard m (2 * h - 1) guardLength heap9 heap10 hCassemble hP19
        hCAssembleGuard hassemble
      exact sameU64Prefix_trans hPrepFrame
        (sameU64Prefix_trans hP0Frame
          (sameU64Prefix_trans hP1Frame
            (sameU64Prefix_trans hHighFrame
              (sameU64Prefix_trans hSub0Frame
                (sameU64Prefix_trans hSub1Frame
                  (sameU64Prefix_trans hCopyFrame
                    (sameU64Prefix_trans hZeroFrame hAssembleFrame)))))))
termination_by n
decreasing_by
  · exact hmLt
  · exact hhLt
  · exact hhLt

theorem karMul_preserves_prefix (this : DenseUPolyZp)
    (C A B : RawPtr UInt64) (n : Nat) (scratch guard : RawPtr UInt64)
    (guardLength : Nat) (heap heap' : RawHeap) (hn : 0 < n)
    (hC : heap.ValidU64Slice C (2 * n - 1))
    (hA : heap.ValidU64Slice A n)
    (hB : heap.ValidU64Slice B n)
    (hScratch : heap.ValidU64Slice scratch (karScratchNeed n))
    (hCGuard : U64SlicesDisjoint C (2 * n - 1) guard guardLength)
    (hScratchGuard : U64SlicesDisjoint scratch (karScratchNeed n)
      guard guardLength)
    (hrun : dense_upoly_zp__kar_mul_ir this C A B n scratch heap = .ok heap') :
    SameU64Prefix heap heap' guard guardLength := by
  rcases karMul_ok this C A B n scratch heap hn hC hA hB hScratch with
    ⟨okHeap, hok, _, hframe⟩
  have heq : okHeap = heap' := Except.ok.inj (hok.symm.trans hrun)
  subst okHeap
  exact hframe guard guardLength hCGuard hScratchGuard

theorem classicalMul_refines_slice (this : DenseUPolyZp)
    (C A : RawPtr UInt64) (lenA : Nat) (B : RawPtr UInt64) (lenB : Nat)
    (heap : RawHeap) (left right : Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this)
    (hp : 1 < this._p.toNat)
    (hApos : 0 < lenA) (hBpos : 0 < lenB)
    (hLenAWord : lenA < limbBase)
    (hC : heap.ValidU64Slice C (lenA + lenB - 1))
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB)
    (hCA : U64SlicesDisjoint C (lenA + lenB - 1) A lenA)
    (hCB : U64SlicesDisjoint C (lenA + lenB - 1) B lenB)
    (hCanonicalA : CanonicalU64Prefix heap A lenA this._p)
    (hCanonicalB : CanonicalU64Prefix heap B lenB this._p)
    (hRepA : SlicePolyRep heap A lenA this._p.toNat left)
    (hRepB : SlicePolyRep heap B lenB this._p.toNat right) :
    ∃ heap', dense_upoly_zp__classical_mul_ir this C A lenA B lenB heap =
        .ok heap' ∧ RawHeap.SameLayout heap heap' ∧
      SlicePolyRep heap' C (lenA + lenB - 1) this._p.toNat
        (left * right) ∧
      CanonicalU64Prefix heap' C (lenA + lenB - 1) this._p := by
  have hempty : ClassicalCoeffPrefix heap C 0 (left * right) := by
    intro _ hi
    omega
  rcases classicalOuterLoop_refines_coeff_prefix this C A B lenA lenB
      (lenA + lenB - 1) 0 heap left right hcfg hp hApos hBpos hLenAWord
      rfl hC hA hB hCA hCB hCanonicalA hCanonicalB hRepA hRepB hempty with
    ⟨heap', hrun, hlayout, hprefix⟩
  have hvalid' := (hlayout C (lenA + lenB - 1)).mp hC
  have hcanonical' := canonicalU64Prefix_of_classicalCoeffPrefix heap' C
    (lenA + lenB - 1) this._p (left * right) hprefix
  have hslice := slicePolyRep_of_classicalCoeffPrefix heap' C
    (lenA + lenB - 1) (left * right) hvalid' hprefix
    (by
      intro degree hdegree
      exact mul_coeff_zero_of_slice_lengths heap A B lenA lenB degree
        left right hApos hBpos hRepA hRepB hdegree)
  refine ⟨heap', ?_, hlayout, hslice, hcanonical'⟩
  simp [dense_upoly_zp__classical_mul_ir, Nat.ne_of_gt hApos,
    Nat.ne_of_gt hBpos, hrun]

theorem classicalMul_refines (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (C A : RawPtr UInt64) (lenA : Nat) (B : RawPtr UInt64) (lenB : Nat)
    (heap : RawHeap) (left right : Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this)
    (hp : 1 < this._p.toNat)
    (hApos : 0 < lenA) (hBpos : 0 < lenB)
    (hLenAWord : lenA < limbBase)
    (hC : heap.ValidU64Slice C (lenA + lenB - 1))
    (hCA : C.region ≠ A.region) (hCB : C.region ≠ B.region)
    (hLeft : RawDensePolyRep this heap A lenA left)
    (hRight : RawDensePolyRep this heap B lenB right) :
    ∃ heap', dense_upoly_zp__classical_mul_ir this C A lenA B lenB heap =
        .ok heap' ∧ RawHeap.SameLayout heap heap' ∧
      RawDensePolyRep this heap' C (lenA + lenB - 1) (left * right) := by
  have hempty : ClassicalCoeffPrefix heap C 0 (left * right) := by
    intro _ hi
    omega
  rcases classicalOuterLoop_refines_coeff_prefix this C A B lenA lenB
      (lenA + lenB - 1) 0 heap left right hcfg hp hApos hBpos hLenAWord
      rfl hC hLeft.1 hRight.1
      (u64SlicesDisjoint_of_region_ne hCA)
      (u64SlicesDisjoint_of_region_ne hCB) hLeft.2.1 hRight.2.1
      hLeft.2.2.1 hRight.2.2.1 hempty with
    ⟨heap', hrun, hlayout, hprefix⟩
  have hvalid' := (hlayout C (lenA + lenB - 1)).mp hC
  have hcanonical' := canonicalU64Prefix_of_classicalCoeffPrefix heap' C
    (lenA + lenB - 1) this._p (left * right) hprefix
  have hslice' := slicePolyRep_of_classicalCoeffPrefix heap' C
    (lenA + lenB - 1) (left * right) hvalid' hprefix
    (by
      intro degree hdegree
      exact mul_coeff_zero_of_slice_lengths heap A B lenA lenB degree
        left right hApos hBpos hLeft.2.2.1 hRight.2.2.1 hdegree)
  have hlast := mul_last_coeff_ne_zero_of_rawDense this heap A B lenA lenB
    left right hApos hBpos hLeft hRight
  have hnorm' := normaliseU64_eq_length_of_classicalCoeffPrefix heap' C
    (lenA + lenB - 1) (left * right) (by omega) hprefix hlast
  refine ⟨heap', ?_, hlayout, hvalid', hcanonical', hslice', hnorm'⟩
  simp [dense_upoly_zp__classical_mul_ir, Nat.ne_of_gt hApos,
    Nat.ne_of_gt hBpos, hrun]

end CLPoly.Impl.StrictMulRefinement
