import CLPoly.Impl.RawPolynomialRep

namespace CLPoly.Impl.StrictDivremRefinement

open Generated.StrictDivrem
open CLPoly.Impl.RawPolynomialRep

/-- A valid raw coefficient slice has a safe observation of exactly the C++
declared length. -/
theorem readU64s_ok (heap : RawHeap) (ptr : RawPtr UInt64) (count : Nat)
    (hvalid : heap.ValidU64Slice ptr count) :
    ∃ coeffs, heap.readU64s ptr count = .ok coeffs ∧ coeffs.size = count := by
  cases count with
  | zero => exact ⟨#[], rfl, rfl⟩
  | succ n =>
    simp only [RawHeap.readU64s]
    rcases heap.readU64_of_valid ptr (n + 1) 0 hvalid (by omega) with
      ⟨head, hhead⟩
    simp only [hhead]
    have htail := heap.validU64Slice_add ptr (n + 1) 1 n hvalid (by omega)
    rcases readU64s_ok heap (RawPtr.add ptr 1) n htail with
      ⟨tail, hreadTail, hsizeTail⟩
    simp only [hreadTail]
    refine ⟨#[head] ++ tail, rfl, ?_⟩
    simp [hsizeTail, Nat.add_comm]

theorem readU64s_get (heap : RawHeap) (ptr : RawPtr UInt64) (count : Nat)
    (coeffs : Array UInt64) (hread : heap.readU64s ptr count = .ok coeffs)
    (hsize : coeffs.size = count) (index : Nat) (hindex : index < count) :
    heap.readU64 ptr index = .ok coeffs[index] := by
  cases count with
  | zero => omega
  | succ n =>
    simp only [RawHeap.readU64s] at hread
    cases hhead : heap.readU64 ptr 0 with
    | error fault => simp [hhead] at hread
    | ok head =>
      cases htail : heap.readU64s (RawPtr.add ptr 1) n with
      | error fault => simp [hhead, htail] at hread
      | ok tail =>
        simp only [hhead, htail] at hread
        have hcoeffs : coeffs = #[head] ++ tail := (Except.ok.inj hread).symm
        subst coeffs
        have htailSizeEq : 1 + tail.size = n + 1 := by simpa using hsize
        have htailSize : tail.size = n := by omega
        cases index with
        | zero => simpa using hhead
        | succ j =>
          have hj : j < n := by omega
          have ih := readU64s_get heap (RawPtr.add ptr 1) n tail htail
            htailSize j hj
          rw [RawHeap.readU64_add] at ih
          simpa [Nat.add_comm] using ih

theorem sliceRep_exists_unique (heap : RawHeap) (ptr : RawPtr UInt64)
    (length : Nat) (hvalid : heap.ValidU64Slice ptr length) :
    ∃ coeffs : Array UInt64,
      (heap.SliceRep ptr length coeffs ∧ coeffs.size = length) ∧
      ∀ other : Array UInt64,
        heap.SliceRep ptr length other ∧ other.size = length → other = coeffs := by
  rcases readU64s_ok heap ptr length hvalid with ⟨coeffs, hread, hsize⟩
  refine ⟨coeffs, ⟨hread, hsize⟩, ?_⟩
  intro other hother
  exact (Except.ok.inj (hread.symm.trans hother.1)).symm

theorem slicePolyRep_exists_unique (heap : RawHeap) (ptr : RawPtr UInt64)
    (length p : Nat) (hvalid : heap.ValidU64Slice ptr length) :
    ∃ poly : Polynomial (ZMod p),
      SlicePolyRep heap ptr length p poly ∧
      ∀ other : Polynomial (ZMod p),
        SlicePolyRep heap ptr length p other → other = poly := by
  rcases sliceRep_exists_unique heap ptr length hvalid with
    ⟨coeffs, hcoeffs, hunique⟩
  let poly := coeffArrayPoly p coeffs
  refine ⟨poly, ⟨coeffs, hcoeffs.1, hcoeffs.2, rfl⟩, ?_⟩
  intro other hother
  rcases hother with ⟨otherCoeffs, hreadOther, hsizeOther, hpolyOther⟩
  have heq : otherCoeffs = coeffs :=
    hunique otherCoeffs ⟨hreadOther, hsizeOther⟩
  subst otherCoeffs
  exact hpolyOther

theorem slicePolyRep_coeff (heap : RawHeap) (ptr : RawPtr UInt64)
    (length p : Nat) (poly : Polynomial (ZMod p))
    (hrep : SlicePolyRep heap ptr length p poly)
    (degree : Nat) (hdegree : degree < length) :
    ∃ value : UInt64,
      heap.readU64 ptr degree = .ok value ∧
      poly.coeff degree = (value.toNat : ZMod p) := by
  rcases hrep with ⟨coeffs, hread, hsize, hpoly⟩
  have hraw := readU64s_get heap ptr length coeffs hread hsize degree hdegree
  refine ⟨coeffs[degree], hraw, ?_⟩
  rw [hpoly, coeff_coeffArrayPoly, dif_pos]

/-- `_poly_normalise` only reads within its declared prefix.  Structural
recursion on `len` proves both successful execution and the returned prefix
bound; structural recursion produces no alternate execution result. -/
theorem normaliseU64_ok (heap : RawHeap) (ptr : RawPtr UInt64) (len : Nat)
    (hvalid : heap.ValidU64Slice ptr len) :
    ∃ result, heap.normaliseU64 ptr len = .ok result ∧ result ≤ len := by
  cases len with
  | zero => exact ⟨0, rfl, Nat.le_refl 0⟩
  | succ n =>
    simp only [RawHeap.normaliseU64]
    rcases heap.readU64_of_valid ptr (n + 1) n hvalid (by omega) with
      ⟨value, hread⟩
    simp only [hread]
    split
    next hzero =>
      have hprefix := heap.validU64Slice_mono ptr (n + 1) n hvalid (by omega)
      rcases normaliseU64_ok heap ptr n hprefix with ⟨result, hresult, hle⟩
      exact ⟨result, hresult, by omega⟩
    next hnonzero => exact ⟨n + 1, rfl, Nat.le_refl _⟩

/-- Limb `memcpy` succeeds for valid source and destination slices.  The
returned heap has exactly the same allocation layout as the input heap, so
all caller slice invariants can be transported through the copy. -/
theorem copyU64_ok (heap : RawHeap) (dst src : RawPtr UInt64) (count : Nat)
    (hDst : heap.ValidU64Slice dst count)
    (hSrc : heap.ValidU64Slice src count) :
    ∃ heap', heap.copyU64 dst src count = .ok heap' ∧
      RawHeap.SameLayout heap heap' := by
  cases count with
  | zero =>
    refine ⟨heap, rfl, ?_⟩
    intro ptr length
    exact Iff.rfl
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
    have hDstTail1 : heap1.ValidU64Slice (RawPtr.add dst 1) n :=
      (hlayout1 (RawPtr.add dst 1) n).mp hDstTail0
    have hSrcTail1 : heap1.ValidU64Slice (RawPtr.add src 1) n :=
      (hlayout1 (RawPtr.add src 1) n).mp hSrcTail0
    rcases copyU64_ok heap1 (RawPtr.add dst 1) (RawPtr.add src 1) n
      hDstTail1 hSrcTail1 with ⟨heap2, hcopy, hlayout2⟩
    simp only [hcopy]
    refine ⟨heap2, rfl, ?_⟩
    intro ptr length
    exact (hlayout1 ptr length).trans (hlayout2 ptr length)

/-- Natural-language proof outline:

At an iteration with `i < lenA`, validity of `A[0..lenA)` gives a successful
read of `A[i]`; validity of the `lenA`-element W3 slice gives a successful
three-limb write at `W3[i]`.  `writeWord3_preserves_valid` shows that this
write changes no allocation sizes, so both slice invariants hold for the
recursive heap.  The recursive measure `lenA - (i+1)` is smaller than
`lenA-i`.  When `i ≥ lenA`, the source loop exits without a heap access. -/
theorem initW3Loop_ok (heap : RawHeap) (A : RawPtr UInt64)
    (W3 : RawPtr Word3) (lenA i : Nat)
    (hA : heap.ValidU64Slice A lenA)
    (hW3 : heap.ValidWord3Slice W3 lenA) (hi : i ≤ lenA) :
    ∃ heap', initW3Loop heap A W3 lenA i = .ok heap' ∧
      heap'.ValidU64Slice A lenA ∧ heap'.ValidWord3Slice W3 lenA ∧
      RawHeap.SameLayout heap heap' := by
  rw [initW3Loop]
  split
  next hlt =>
    rcases heap.readU64_of_valid A lenA i hA hlt with ⟨value, hread⟩
    simp only [hread]
    rcases heap.writeWord3_of_valid W3 lenA i
      { lo := value, mid := 0, hi := 0 } hW3 hlt with ⟨heap1, hwrite⟩
    simp only [hwrite]
    have hA1 : heap1.ValidU64Slice A lenA :=
      (RawHeap.writeWord3_preserves_valid heap heap1 W3 i
        { lo := value, mid := 0, hi := 0 } hwrite A lenA).mp hA
    have hW31 : heap1.ValidWord3Slice W3 lenA :=
      (RawHeap.writeWord3_preserves_valid heap heap1 W3 i
        { lo := value, mid := 0, hi := 0 } hwrite
        (RawPtr.reinterpret W3) (3 * lenA)).mp hW3
    have hlayout1 := RawHeap.writeWord3_sameLayout heap heap1 W3 i
      { lo := value, mid := 0, hi := 0 } hwrite
    rcases initW3Loop_ok heap1 A W3 lenA (i + 1) hA1 hW31 (by omega) with
      ⟨heap2, hloop, hA2, hW32, hlayout2⟩
    refine ⟨heap2, hloop, hA2, hW32, ?_⟩
    intro ptr length
    exact (hlayout1 ptr length).trans (hlayout2 ptr length)
  next hnot =>
    exact ⟨heap, rfl, hA, hW3, fun _ _ => Iff.rfl⟩
termination_by lenA - i
decreasing_by omega

/-- The already initialized W3 prefix is the exact zero-extended image of
the corresponding source coefficient prefix. -/
def InitW3Prefix (heap : RawHeap) (A : RawPtr UInt64)
    (W3 : RawPtr Word3) (upto : Nat) : Prop :=
  ∀ j, j < upto → ∃ value : UInt64,
    heap.readU64 A j = .ok value ∧
    heap.readWord3 W3 j = .ok { lo := value, mid := 0, hi := 0 }

/-- Content-level refinement of the generated initialization loop.  With
the C++ non-aliasing allocation precondition, every output W3 cell is exactly
the corresponding A limb zero-extended to three limbs. -/
theorem initW3Loop_refines (heap : RawHeap) (A : RawPtr UInt64)
    (W3 : RawPtr Word3) (lenA i : Nat)
    (hA : heap.ValidU64Slice A lenA)
    (hW3 : heap.ValidWord3Slice W3 lenA) (hi : i ≤ lenA)
    (hregions : W3.region ≠ A.region)
    (hprefix : InitW3Prefix heap A W3 i) :
    ∃ heap', initW3Loop heap A W3 lenA i = .ok heap' ∧
      heap'.ValidU64Slice A lenA ∧ heap'.ValidWord3Slice W3 lenA ∧
      RawHeap.SameLayout heap heap' ∧ InitW3Prefix heap' A W3 lenA := by
  rw [initW3Loop]
  split
  next hlt =>
    rcases heap.readU64_of_valid A lenA i hA hlt with ⟨value, hread⟩
    simp only [hread]
    let word : Word3 := { lo := value, mid := 0, hi := 0 }
    rcases heap.writeWord3_of_valid W3 lenA i word hW3 hlt with
      ⟨heap1, hwrite⟩
    dsimp [word] at hwrite ⊢
    simp only [hwrite]
    have hA1 : heap1.ValidU64Slice A lenA :=
      (RawHeap.writeWord3_preserves_valid heap heap1 W3 i word
        hwrite A lenA).mp hA
    have hW31 : heap1.ValidWord3Slice W3 lenA :=
      (RawHeap.writeWord3_preserves_valid heap heap1 W3 i word
        hwrite (RawPtr.reinterpret W3) (3 * lenA)).mp hW3
    have hlayout1 := RawHeap.writeWord3_sameLayout heap heap1 W3 i word hwrite
    have hprefix1 : InitW3Prefix heap1 A W3 (i + 1) := by
      intro j hj
      by_cases hji : j = i
      · subst j
        have hreadA := RawHeap.readU64_writeWord3_region_ne heap heap1
          W3 A i i word value hwrite hread hregions
        have hreadW := RawHeap.readWord3_writeWord3_same heap heap1
          W3 i word hwrite
        exact ⟨value, hreadA, hreadW⟩
      · have hjlt : j < i := by omega
        rcases hprefix j hjlt with ⟨old, hreadA, hreadW⟩
        have hreadA1 := RawHeap.readU64_writeWord3_region_ne heap heap1
          W3 A i j word old hwrite hreadA hregions
        have hreadW1 := RawHeap.readWord3_writeWord3_ne heap heap1
          W3 i j word { lo := old, mid := 0, hi := 0 }
          hwrite hreadW (Ne.symm hji)
        exact ⟨old, hreadA1, hreadW1⟩
    rcases initW3Loop_refines heap1 A W3 lenA (i + 1) hA1 hW31
      (by omega) hregions hprefix1 with
      ⟨heap2, hloop, hA2, hW32, hlayout2, hfull⟩
    refine ⟨heap2, hloop, hA2, hW32, ?_, hfull⟩
    intro ptr length
    exact (hlayout1 ptr length).trans (hlayout2 ptr length)
  next hnot =>
    have hieq : i = lenA := by omega
    subst i
    exact ⟨heap, rfl, hA, hW3, fun _ _ => Iff.rfl, hprefix⟩
termination_by lenA - i
decreasing_by omega

/-- Natural-language proof outline:

For `j ≤ d`, the divisor invariant makes `B[j]` readable.  The bound
`i+d < lenW3` implies `i+j < lenW3`, so the accumulator is readable and
writable.  The multiplication and carry routines are total word operations.
The W3 write preserves both B and W3 allocation invariants; recurse at `j+1`
with the strictly smaller measure `d+1-j`.  At `j>d`, the source loop exits. -/
theorem addMulLoop_ok (heap : RawHeap) (B : RawPtr UInt64)
    (W3 : RawPtr Word3) (lenW3 i d j : Nat) (c : UInt64)
    (other : RawPtr UInt64) (otherLen : Nat)
    (hB : heap.ValidU64Slice B (d + 1))
    (hW3 : heap.ValidWord3Slice W3 lenW3)
    (hOther : heap.ValidU64Slice other otherLen)
    (htop : i + d < lenW3) (hj : j ≤ d + 1) :
    ∃ heap', addMulLoop heap B W3 i d j c = .ok heap' ∧
      heap'.ValidU64Slice B (d + 1) ∧
      heap'.ValidWord3Slice W3 lenW3 ∧
      heap'.ValidU64Slice other otherLen ∧
      RawHeap.SameLayout heap heap' := by
  rw [addMulLoop]
  split
  next hle =>
    have hjB : j < d + 1 := by omega
    have hijW : i + j < lenW3 := by omega
    rcases heap.readU64_of_valid B (d + 1) j hB hjB with ⟨bj, hreadB⟩
    simp only [hreadB]
    rcases heap.readWord3_of_valid W3 lenW3 (i + j) hW3 hijW with
      ⟨accum, hreadW⟩
    simp only [hreadW]
    let product := Generated.StrictGCD.dense_upoly_zp__umul128_ir 0 0 c bj
    let accum' := Generated.StrictGCD.dense_upoly_zp__add_carry3_ir
      accum product.fst product.snd
    rcases heap.writeWord3_of_valid W3 lenW3 (i + j) accum' hW3 hijW with
      ⟨heap1, hwrite⟩
    dsimp [product, accum'] at hwrite ⊢
    simp only [hwrite]
    have hB1 : heap1.ValidU64Slice B (d + 1) :=
      (RawHeap.writeWord3_preserves_valid heap heap1 W3 (i + j) accum'
        hwrite B (d + 1)).mp hB
    have hW31 : heap1.ValidWord3Slice W3 lenW3 :=
      (RawHeap.writeWord3_preserves_valid heap heap1 W3 (i + j) accum'
        hwrite (RawPtr.reinterpret W3) (3 * lenW3)).mp hW3
    have hOther1 : heap1.ValidU64Slice other otherLen :=
      (RawHeap.writeWord3_preserves_valid heap heap1 W3 (i + j) accum'
        hwrite other otherLen).mp hOther
    have hlayout1 := RawHeap.writeWord3_sameLayout heap heap1 W3 (i + j)
      accum' hwrite
    rcases addMulLoop_ok heap1 B W3 lenW3 i d (j + 1) c other otherLen
      hB1 hW31 hOther1 htop (by omega) with
      ⟨heap2, hloop, hB2, hW32, hOther2, hlayout2⟩
    refine ⟨heap2, hloop, hB2, hW32, hOther2, ?_⟩
    intro ptr length
    exact (hlayout1 ptr length).trans (hlayout2 ptr length)
  next hnot =>
    exact ⟨heap, rfl, hB, hW3, hOther, fun _ _ => Iff.rfl⟩
termination_by d + 1 - j
decreasing_by omega

/-- The descending quotient loop cannot fault when Q, B and W3 have the
capacities documented by the C++ raw API.  Its induction variable is exactly
the source `ii`: the successor case processes coefficient `i = ii-1`, then
recurses on the predecessor. -/
theorem quotientLoop_ok (this : DenseUPolyZp) (Q B : RawPtr UInt64)
    (W3 : RawPtr Word3) (qLen d lenW3 : Nat) (invLc : UInt64)
    (heap : RawHeap) (ii : Nat)
    (hQ : heap.ValidU64Slice Q qLen)
    (hB : heap.ValidU64Slice B (d + 1))
    (hW3 : heap.ValidWord3Slice W3 lenW3)
    (hii : ii ≤ qLen) (hspan : qLen + d ≤ lenW3) :
    ∃ heap', quotientLoop this Q B W3 d invLc heap ii = .ok heap' ∧
      heap'.ValidU64Slice Q qLen ∧
      heap'.ValidU64Slice B (d + 1) ∧
      heap'.ValidWord3Slice W3 lenW3 ∧
      RawHeap.SameLayout heap heap' := by
  cases ii with
  | zero => exact ⟨heap, rfl, hQ, hB, hW3, fun _ _ => Iff.rfl⟩
  | succ i =>
    have hiQ : i < qLen := by omega
    have hiW : i + d < lenW3 := by omega
    simp only [quotientLoop]
    rcases heap.readWord3_of_valid W3 lenW3 (i + d) hW3 hiW with
      ⟨accum, hread⟩
    simp only [hread]
    let r := Generated.StrictGCD.dense_upoly_zp__lll_mod_preinv_ir
      accum.hi accum.mid accum.lo this._p this._ninv this._norm
    let qi := Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this r invLc
    rcases heap.writeU64_of_valid Q qLen i qi hQ hiQ with ⟨heap1, hwrite⟩
    dsimp [r, qi] at hwrite ⊢
    simp only [hwrite]
    have hQ1 : heap1.ValidU64Slice Q qLen :=
      (RawHeap.writeU64_preserves_valid heap heap1 Q i _ hwrite Q qLen).mp hQ
    have hB1 : heap1.ValidU64Slice B (d + 1) :=
      (RawHeap.writeU64_preserves_valid heap heap1 Q i _ hwrite B (d + 1)).mp hB
    have hW31 : heap1.ValidWord3Slice W3 lenW3 :=
      (RawHeap.writeU64_preserves_valid heap heap1 Q i _ hwrite
        (RawPtr.reinterpret W3) (3 * lenW3)).mp hW3
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 Q i _ hwrite
    split
    next hnonzero =>
      rcases addMulLoop_ok heap1 B W3 lenW3 i d 0 (this._p -
          Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this
            (Generated.StrictGCD.dense_upoly_zp__lll_mod_preinv_ir
              accum.hi accum.mid accum.lo this._p this._ninv this._norm) invLc)
          Q qLen hB1 hW31 hQ1 hiW (by omega) with
        ⟨heap2, hadd, hB2, hW32, hQ2, hlayout2⟩
      simp only [hadd]
      rcases quotientLoop_ok this Q B W3 qLen d lenW3 invLc heap2 i
        hQ2 hB2 hW32 (by omega) hspan with
        ⟨heap3, hloop, hQ3, hB3, hW33, hlayout3⟩
      refine ⟨heap3, hloop, hQ3, hB3, hW33, ?_⟩
      intro ptr length
      exact (hlayout1 ptr length).trans
        ((hlayout2 ptr length).trans (hlayout3 ptr length))
    next hzero =>
      rcases quotientLoop_ok this Q B W3 qLen d lenW3 invLc heap1 i
        hQ1 hB1 hW31 (by omega) hspan with
        ⟨heap2, hloop, hQ2, hB2, hW32, hlayout2⟩
      refine ⟨heap2, hloop, hQ2, hB2, hW32, ?_⟩
      intro ptr length
      exact (hlayout1 ptr length).trans (hlayout2 ptr length)

/-- The final C++ remainder loop reads exactly W3[0..d) and writes R[0..d).
Both layouts are preserved at each iteration, and `d-i` decreases. -/
theorem remainderLoop_ok (this : DenseUPolyZp) (R : RawPtr UInt64)
    (W3 : RawPtr Word3) (d lenW3 i : Nat) (heap : RawHeap)
    (hR : heap.ValidU64Slice R d)
    (hW3 : heap.ValidWord3Slice W3 lenW3)
    (hdW : d ≤ lenW3) (hi : i ≤ d) :
    ∃ heap', remainderLoop this R W3 d i heap = .ok heap' ∧
      heap'.ValidU64Slice R d ∧ heap'.ValidWord3Slice W3 lenW3 ∧
      RawHeap.SameLayout heap heap' := by
  rw [remainderLoop]
  split
  next hlt =>
    have hiW : i < lenW3 := by omega
    rcases heap.readWord3_of_valid W3 lenW3 i hW3 hiW with ⟨accum, hread⟩
    simp only [hread]
    let value := Generated.StrictGCD.dense_upoly_zp__lll_mod_preinv_ir
      accum.hi accum.mid accum.lo this._p this._ninv this._norm
    rcases heap.writeU64_of_valid R d i value hR hlt with ⟨heap1, hwrite⟩
    dsimp [value] at hwrite ⊢
    simp only [hwrite]
    have hR1 : heap1.ValidU64Slice R d :=
      (RawHeap.writeU64_preserves_valid heap heap1 R i _ hwrite R d).mp hR
    have hW31 : heap1.ValidWord3Slice W3 lenW3 :=
      (RawHeap.writeU64_preserves_valid heap heap1 R i _ hwrite
        (RawPtr.reinterpret W3) (3 * lenW3)).mp hW3
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 R i _ hwrite
    rcases remainderLoop_ok this R W3 d lenW3 (i + 1) heap1
      hR1 hW31 hdW (by omega) with
      ⟨heap2, hloop, hR2, hW32, hlayout2⟩
    refine ⟨heap2, hloop, hR2, hW32, ?_⟩
    intro ptr length
    exact (hlayout1 ptr length).trans (hlayout2 ptr length)
  next hnot => exact ⟨heap, rfl, hR, hW3, fun _ _ => Iff.rfl⟩
termination_by d - i
decreasing_by omega

/-- Under exactly the capacities documented on the C++ raw API,
`_poly_divrem` cannot take `RawFault`.  This theorem is only the termination
and memory-safety bridge; the quotient/remainder algebraic invariant is the
next refinement obligation. -/
theorem polyDivrem_ok (this : DenseUPolyZp) (Q R A B : RawPtr UInt64)
    (lenA lenB : Nat) (W3 : RawPtr Word3) (heap : RawHeap)
    (hlenB : 0 < lenB)
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB)
    (hQ : heap.ValidU64Slice Q (lenA - (lenB - 1)))
    (hR : heap.ValidU64Slice R (Nat.min lenA (lenB - 1)))
    (hW3 : heap.ValidWord3Slice W3 lenA) :
    ∃ heap' lenQ lenR,
      dense_upoly_zp__poly_divrem_ir this Q R A lenA B lenB W3 heap =
        .ok (heap', lenQ, lenR) ∧
      lenQ ≤ lenA - (lenB - 1) ∧
      lenR ≤ Nat.min lenA (lenB - 1) := by
  cases lenB with
  | zero => omega
  | succ d =>
    simp only [dense_upoly_zp__poly_divrem_ir]
    split
    next hshort =>
      have hlenAd : lenA ≤ d := by omega
      have hRfull : heap.ValidU64Slice R lenA := by
        simpa [Nat.min_eq_left hlenAd] using hR
      rcases copyU64_ok heap R A lenA hRfull hA with
        ⟨heap1, hcopy, hlayout⟩
      simp only [hcopy]
      exact ⟨heap1, 0, lenA, rfl, by omega, by simpa [Nat.min_eq_left hlenAd]⟩
    next hlong =>
      have hdA : d < lenA := by omega
      have hdleA : d ≤ lenA := Nat.le_of_lt hdA
      have hRfull : heap.ValidU64Slice R d := by
        simpa [Nat.min_eq_right hdleA] using hR
      rcases heap.readU64_of_valid B (d + 1) d hB (by omega) with
        ⟨lead, hlead⟩
      simp only [hlead]
      let invLc := Generated.StrictGCD.dense_upoly_zp_nmod_inv_ir this lead
      rcases initW3Loop_ok heap A W3 lenA 0 hA hW3 (by omega) with
        ⟨heap1, hinit, hA1, hW31, hlayout1⟩
      simp only [hinit]
      have hQ1 : heap1.ValidU64Slice Q (lenA - d) :=
        (hlayout1 Q (lenA - d)).mp (by simpa using hQ)
      have hB1 : heap1.ValidU64Slice B (d + 1) :=
        (hlayout1 B (d + 1)).mp hB
      have hR1 : heap1.ValidU64Slice R d := (hlayout1 R d).mp hRfull
      rcases quotientLoop_ok this Q B W3 (lenA - d) d lenA invLc heap1
        (lenA - d) hQ1 hB1 hW31 (Nat.le_refl _) (by omega) with
        ⟨heap2, hquot, hQ2, hB2, hW32, hlayout2⟩
      dsimp [invLc] at hquot ⊢
      simp only [hquot]
      have hR2 : heap2.ValidU64Slice R d := (hlayout2 R d).mp hR1
      rcases remainderLoop_ok this R W3 d lenA 0 heap2 hR2 hW32
        (by omega) (by omega) with
        ⟨heap3, hrem, hR3, hW33, hlayout3⟩
      simp only [hrem]
      have hQ3 : heap3.ValidU64Slice Q (lenA - d) :=
        (hlayout3 Q (lenA - d)).mp ((hlayout2 Q (lenA - d)).mp hQ1)
      rcases normaliseU64_ok heap3 Q (lenA - d) hQ3 with
        ⟨lenQ, hnormQ, hlenQ⟩
      simp only [hnormQ]
      rcases normaliseU64_ok heap3 R d hR3 with ⟨lenR, hnormR, hlenR⟩
      simp only [hnormR]
      refine ⟨heap3, lenQ, lenR, rfl, ?_, ?_⟩
      · simpa using hlenQ
      · simpa [Nat.min_eq_right hdleA] using hlenR

end CLPoly.Impl.StrictDivremRefinement
