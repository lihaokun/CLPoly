import CLPoly.Impl.StrictDivremRefinement

set_option autoImplicit false

namespace CLPoly.Impl.StrictEuclidRefinement

open Generated.StrictDivrem
open Generated.StrictEuclidGCD
open CLPoly.Impl.RawPolynomialRep
open CLPoly.Impl.StrictDivremRefinement
open CLPoly.Impl.StrictWordArithmetic

/-- The raw representation invariant of a normalized C++ `dense_upoly_zp`
coefficient buffer.  Every field is an observable L1 property. -/
def RawDensePolyRep (this : DenseUPolyZp) (heap : RawHeap)
    (ptr : RawPtr UInt64) (length : Nat)
    (poly : Polynomial (ZMod this._p.toNat)) : Prop :=
  heap.ValidU64Slice ptr length ∧
    CanonicalU64Prefix heap ptr length this._p ∧
    SlicePolyRep heap ptr length this._p.toNat poly ∧
    heap.normaliseU64 ptr length = .ok length

/-- Observable representation of a fixed-length canonical coefficient slice.
Unlike `RawDensePolyRep`, trailing zero limbs are permitted. -/
def RawCanonicalPolySlice (this : DenseUPolyZp) (heap : RawHeap)
    (ptr : RawPtr UInt64) (length : Nat)
    (poly : Polynomial (ZMod this._p.toNat)) : Prop :=
  heap.ValidU64Slice ptr length ∧
    CanonicalU64Prefix heap ptr length this._p ∧
    SlicePolyRep heap ptr length this._p.toNat poly

/-- A normalized raw polynomial survives a heap transition that preserves
its allocation and every coefficient in its declared prefix. -/
theorem rawDensePolyRep_of_same_prefix (this : DenseUPolyZp)
    (before after : RawHeap) (ptr : RawPtr UInt64) (length : Nat)
    (poly : Polynomial (ZMod this._p.toNat))
    (hlayout : RawHeap.SameLayout before after)
    (hsame : SameU64Prefix before after ptr length)
    (hrep : RawDensePolyRep this before ptr length poly) :
    RawDensePolyRep this after ptr length poly := by
  have hvalidAfter := (hlayout ptr length).mp hrep.1
  refine ⟨hvalidAfter, ?_,
    slicePolyRep_of_same_prefix before after ptr length this._p.toNat poly
      hrep.1 hvalidAfter hsame hrep.2.2.1, ?_⟩
  · intro k value hk hreadAfter
    rcases before.readU64_of_valid ptr length k hrep.1 hk with
      ⟨old, hreadBefore⟩
    have hreadOld := hsame k old hk hreadBefore
    have hvalue : value = old := Except.ok.inj (hreadAfter.symm.trans hreadOld)
    subst value
    exact hrep.2.1 k old hk hreadBefore
  · have hnorm := normaliseU64_eq_of_prefix_map before after ptr ptr length
      hrep.1 hsame
    exact hnorm.symm.trans hrep.2.2.2

/-- The source `memcpy` transports the complete normalized polynomial
representation between distinct C++ allocations. -/
theorem copyU64_refines_rawDense_of_region_ne (this : DenseUPolyZp)
    (heap : RawHeap) (dst src : RawPtr UInt64) (length : Nat)
    (poly : Polynomial (ZMod this._p.toNat))
    (hDst : heap.ValidU64Slice dst length)
    (hregion : dst.region ≠ src.region)
    (hrep : RawDensePolyRep this heap src length poly) :
    ∃ heap', heap.copyU64 dst src length = .ok heap' ∧
      RawHeap.SameLayout heap heap' ∧
      RawDensePolyRep this heap' dst length poly := by
  rcases copyU64_refines heap dst src length hDst hrep.1 hregion with
    ⟨heap', hcopy, hlayout, hcontents⟩
  rcases copyU64_slicePolyRep heap dst src length this._p.toNat poly hDst
      hrep.1 hregion hrep.2.2.1 with
    ⟨repHeap, hcopyRep, _, hslice⟩
  have heq : repHeap = heap' := Except.ok.inj (hcopyRep.symm.trans hcopy)
  subst repHeap
  have hvalidDst := (hlayout dst length).mp hDst
  have hcanonical : CanonicalU64Prefix heap' dst length this._p := by
    intro k value hk hreadAfter
    rcases heap.readU64_of_valid src length k hrep.1 hk with
      ⟨old, hreadBefore⟩
    have hreadOld := hcontents k old hk hreadBefore
    have hvalue : value = old := Except.ok.inj (hreadAfter.symm.trans hreadOld)
    subst value
    exact hrep.2.1 k old hk hreadBefore
  have hnorm := normaliseU64_eq_of_prefix_map heap heap' src dst length
    hrep.1 hcontents
  exact ⟨heap', hcopy, hlayout, hvalidDst, hcanonical, hslice,
    hnorm.symm.trans hrep.2.2.2⟩

/-- A represented polynomial in an allocation distinct from the destination
is a frame of the actual source `memcpy`. -/
theorem copyU64_preserves_rawDense_of_region_ne (this : DenseUPolyZp)
    (heap heap' : RawHeap) (dst src : RawPtr UInt64) (count : Nat)
    (ptr : RawPtr UInt64) (length : Nat)
    (poly : Polynomial (ZMod this._p.toNat))
    (hDst : heap.ValidU64Slice dst count)
    (hSrc : heap.ValidU64Slice src count)
    (hregion : dst.region ≠ ptr.region)
    (hcopy : heap.copyU64 dst src count = .ok heap')
    (hrep : RawDensePolyRep this heap ptr length poly) :
    RawHeap.SameLayout heap heap' ∧
      RawDensePolyRep this heap' ptr length poly := by
  rcases copyU64_ok heap dst src count hDst hSrc with
    ⟨copyHeap, hcopy', hlayout⟩
  have heq : copyHeap = heap' := Except.ok.inj (hcopy'.symm.trans hcopy)
  subst copyHeap
  refine ⟨hlayout, rawDensePolyRep_of_same_prefix this heap heap' ptr length
    poly hlayout ?_ hrep⟩
  intro k value hk hread
  exact copyU64_preserves_read heap heap' dst src ptr count k value hDst hSrc
    hread (by intro _ _; exact Or.inl hregion) hcopy

theorem RawDensePolyRep.toCanonicalSlice (this : DenseUPolyZp)
    (heap : RawHeap) (ptr : RawPtr UInt64) (length : Nat)
    (poly : Polynomial (ZMod this._p.toNat))
    (hrep : RawDensePolyRep this heap ptr length poly) :
    RawCanonicalPolySlice this heap ptr length poly :=
  ⟨hrep.1, hrep.2.1, hrep.2.2.1⟩

/-- Pairwise non-aliasing of the five physical work areas used by the source
Euclid loop. -/
structure EuclidRegions (Q : RawPtr UInt64) (W3 : RawPtr Word3)
    (A B R : RawPtr UInt64) : Prop where
  q_w3 : Q.region ≠ W3.region
  q_a : Q.region ≠ A.region
  q_b : Q.region ≠ B.region
  q_r : Q.region ≠ R.region
  w3_a : W3.region ≠ A.region
  w3_b : W3.region ≠ B.region
  w3_r : W3.region ≠ R.region
  a_b : A.region ≠ B.region
  a_r : A.region ≠ R.region
  b_r : B.region ≠ R.region

/-- Fixed physical capacities allocated for the whole source Euclid loop.
Logical lengths shrink, while these allocations and their disjointness remain
unchanged under `(A,B,R) → (B,R,A)`. -/
structure EuclidWorkspace (heap : RawHeap) (Q : RawPtr UInt64)
    (W3 : RawPtr Word3) (A B R : RawPtr UInt64) (capacity : Nat) : Prop where
  validQ : heap.ValidU64Slice Q capacity
  validW3 : heap.ValidWord3Slice W3 capacity
  validA : heap.ValidU64Slice A capacity
  validB : heap.ValidU64Slice B capacity
  validR : heap.ValidU64Slice R capacity
  regions : EuclidRegions Q W3 A B R

theorem EuclidRegions.rotate {Q : RawPtr UInt64} {W3 : RawPtr Word3}
    {A B R : RawPtr UInt64} (h : EuclidRegions Q W3 A B R) :
    EuclidRegions Q W3 B R A := by
  exact {
    q_w3 := h.q_w3
    q_a := h.q_b
    q_b := h.q_r
    q_r := h.q_a
    w3_a := h.w3_b
    w3_b := h.w3_r
    w3_r := h.w3_a
    a_b := h.b_r
    a_r := h.a_b.symm
    b_r := h.a_r.symm
  }

theorem EuclidWorkspace.rotate_of_sameLayout
    {before after : RawHeap} {Q : RawPtr UInt64} {W3 : RawPtr Word3}
    {A B R : RawPtr UInt64} {capacity : Nat}
    (h : EuclidWorkspace before Q W3 A B R capacity)
    (hlayout : RawHeap.SameLayout before after) :
    EuclidWorkspace after Q W3 B R A capacity := by
  exact {
    validQ := (hlayout Q capacity).mp h.validQ
    validW3 := (hlayout (RawPtr.reinterpret W3) (3 * capacity)).mp h.validW3
    validA := (hlayout B capacity).mp h.validB
    validB := (hlayout R capacity).mp h.validR
    validR := (hlayout A capacity).mp h.validA
    regions := h.regions.rotate
  }

theorem EuclidWorkspace.of_sameLayout
    {before after : RawHeap} {Q : RawPtr UInt64} {W3 : RawPtr Word3}
    {A B R : RawPtr UInt64} {capacity : Nat}
    (h : EuclidWorkspace before Q W3 A B R capacity)
    (hlayout : RawHeap.SameLayout before after) :
    EuclidWorkspace after Q W3 A B R capacity := by
  exact {
    validQ := (hlayout Q capacity).mp h.validQ
    validW3 := (hlayout (RawPtr.reinterpret W3) (3 * capacity)).mp h.validW3
    validA := (hlayout A capacity).mp h.validA
    validB := (hlayout B capacity).mp h.validB
    validR := (hlayout R capacity).mp h.validR
    regions := h.regions
  }

/-- One source Euclid rotation preserves the normalized mathematical gcd.

The premise is precisely the algebraic identity established for the output of
the generated C++ `_poly_divrem`; it does not evaluate an L2 remainder or call
an L2 gcd implementation. -/
theorem normalize_gcd_eq_of_division_identity
    {R : Type*} [EuclideanDomain R] [NormalizationMonoid R] [DecidableEq R]
    (a b q r : R) (hdivision : a = q * b + r) :
    normalize (EuclideanDomain.gcd a b) =
      normalize (EuclideanDomain.gcd b r) := by
  rw [normalize_eq_normalize_iff_associated]
  apply associated_of_dvd_dvd
  · apply EuclideanDomain.dvd_gcd
    · exact EuclideanDomain.gcd_dvd_right a b
    · have ha := EuclideanDomain.gcd_dvd_left a b
      have hb := EuclideanDomain.gcd_dvd_right a b
      have hr : r = a - q * b := by rw [hdivision]; ring
      rw [hr]
      exact dvd_sub ha (dvd_mul_of_dvd_right hb q)
  · apply EuclideanDomain.dvd_gcd
    · rw [hdivision]
      exact dvd_add (dvd_mul_of_dvd_right
        (EuclideanDomain.gcd_dvd_left b r) q)
        (EuclideanDomain.gcd_dvd_right b r)
    · exact EuclideanDomain.gcd_dvd_left b r

/-- A successful, generated C++ division step preserves the normalized L2
gcd across the exact buffer rotation performed by the source Euclid loop. -/
theorem polyDivrem_preserves_normalized_gcd (this : DenseUPolyZp)
    [hprime : Fact (Nat.Prime this._p.toNat)]
    (Q R A B : RawPtr UInt64) (lenA lenB : Nat)
    (W3 : RawPtr Word3) (heap : RawHeap)
    (dividend divisor : Polynomial (ZMod this._p.toNat))
    (hlenB : 0 < lenB)
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB)
    (hQ : heap.ValidU64Slice Q (lenA - (lenB - 1)))
    (hR : heap.ValidU64Slice R (Nat.min lenA (lenB - 1)))
    (hW3 : heap.ValidWord3Slice W3 lenA)
    (hcanonicalA : CanonicalU64Prefix heap A lenA this._p)
    (hcanonicalB : CanonicalU64Prefix heap B lenB this._p)
    (hdividend : SlicePolyRep heap A lenA this._p.toNat dividend)
    (hdivisor : SlicePolyRep heap B lenB this._p.toNat divisor)
    (hnormA : heap.normaliseU64 A lenA = .ok lenA)
    (hnormB : heap.normaliseU64 B lenB = .ok lenB)
    (hqCapacity : lenA - (lenB - 1) < limbBase)
    (hRA : R.region ≠ A.region)
    (hWA : W3.region ≠ A.region) (hWB : W3.region ≠ B.region)
    (hQB : Q.region ≠ B.region) (hQW : Q.region ≠ W3.region)
    (hRW : R.region ≠ W3.region) (hRQ : R.region ≠ Q.region)
    (hRB : R.region ≠ B.region)
    (hcfg : DensePreinvConfigured this) :
    ∃ heap' lenQ lenR quotient remainder,
      dense_upoly_zp__poly_divrem_ir this Q R A lenA B lenB W3 heap =
        .ok (heap', lenQ, lenR) ∧
      RawDensePolyRep this heap' Q lenQ quotient ∧
      SlicePolyRep heap' R lenR this._p.toNat remainder ∧
      CanonicalU64Prefix heap' R lenR this._p ∧
      heap'.normaliseU64 R lenR = .ok lenR ∧
      RawHeap.SameLayout heap heap' ∧
      SameU64Prefix heap heap' B lenB ∧
      dividend = quotient * divisor + remainder ∧
      normalize (EuclideanDomain.gcd dividend divisor) =
        normalize (EuclideanDomain.gcd divisor remainder) ∧
      (remainder = 0 ∨ remainder.natDegree < divisor.natDegree) ∧
      lenQ ≤ lenA - (lenB - 1) ∧
      lenR ≤ Nat.min lenA (lenB - 1) ∧ lenR < lenB := by
  rcases polyDivrem_refines this Q R A B lenA lenB W3 heap dividend divisor
      hlenB hA hB hQ hR hW3 hcanonicalA hcanonicalB hdividend hdivisor
      hnormA hnormB hqCapacity hRA hWA hWB hQB hQW hRW hRQ hRB hcfg
      hprime.out with
    ⟨heap', lenQ, lenR, quotient, remainder, hrun, hquotient, hcanonicalQ,
      hnormQ, hremainder, hcanonicalR, hnormR, hlayout, hsameB, hdivision,
      hdegree, hlenQ, hlenRCapacity, hlenR⟩
  exact ⟨heap', lenQ, lenR, quotient, remainder, hrun,
    ⟨(hlayout Q lenQ).mp
        (heap.validU64Slice_mono Q (lenA - (lenB - 1)) lenQ hQ hlenQ),
      hcanonicalQ, hquotient, hnormQ⟩,
    hremainder, hcanonicalR, hnormR, hlayout, hsameB, hdivision,
    normalize_gcd_eq_of_division_identity dividend divisor quotient remainder
      hdivision,
    hdegree, hlenQ, hlenRCapacity, hlenR⟩

/-- Closure theorem used by well-founded Euclid recursion: a real generated
division call turns the physical remainder allocation into a normalized raw
dense polynomial representing the next mathematical divisor. -/
theorem polyDivrem_next_state (this : DenseUPolyZp)
    [hprime : Fact (Nat.Prime this._p.toNat)]
    (Q R A B : RawPtr UInt64) (lenA lenB : Nat)
    (W3 : RawPtr Word3) (heap : RawHeap)
    (dividend divisor : Polynomial (ZMod this._p.toNat))
    (hlenB : 0 < lenB)
    (hARep : RawDensePolyRep this heap A lenA dividend)
    (hBRep : RawDensePolyRep this heap B lenB divisor)
    (hQ : heap.ValidU64Slice Q (lenA - (lenB - 1)))
    (hR : heap.ValidU64Slice R (Nat.min lenA (lenB - 1)))
    (hW3 : heap.ValidWord3Slice W3 lenA)
    (hqCapacity : lenA - (lenB - 1) < limbBase)
    (hRA : R.region ≠ A.region)
    (hWA : W3.region ≠ A.region) (hWB : W3.region ≠ B.region)
    (hQB : Q.region ≠ B.region) (hQW : Q.region ≠ W3.region)
    (hRW : R.region ≠ W3.region) (hRQ : R.region ≠ Q.region)
    (hRB : R.region ≠ B.region)
    (hcfg : DensePreinvConfigured this) :
    ∃ heap' lenQ lenR quotient remainder,
      dense_upoly_zp__poly_divrem_ir this Q R A lenA B lenB W3 heap =
        .ok (heap', lenQ, lenR) ∧
      RawDensePolyRep this heap' Q lenQ quotient ∧
      RawDensePolyRep this heap' B lenB divisor ∧
      RawDensePolyRep this heap' R lenR remainder ∧
      dividend = quotient * divisor + remainder ∧
      normalize (EuclideanDomain.gcd dividend divisor) =
        normalize (EuclideanDomain.gcd divisor remainder) ∧
      RawHeap.SameLayout heap heap' ∧
      lenQ ≤ lenA - (lenB - 1) ∧
      lenR ≤ Nat.min lenA (lenB - 1) ∧ lenR < lenB := by
  rcases hARep with ⟨hA, hcanonicalA, hdividend, hnormA⟩
  rcases hBRep with ⟨hB, hcanonicalB, hdivisor, hnormB⟩
  rcases polyDivrem_preserves_normalized_gcd this Q R A B lenA lenB W3
      heap dividend divisor hlenB hA hB hQ hR hW3 hcanonicalA hcanonicalB
      hdividend hdivisor hnormA hnormB hqCapacity hRA hWA hWB hQB hQW hRW
      hRQ hRB hcfg with
    ⟨heap', lenQ, lenR, quotient, remainder, hrun, hquotient, hremainder,
      hcanonicalR, hnormR, hlayout, hsameB, hdivision, hgcd, _, hlenQ,
      hlenRCapacity, hlenR⟩
  have hBResult : heap'.ValidU64Slice B lenB := (hlayout B lenB).mp hB
  have hcanonicalBResult : CanonicalU64Prefix heap' B lenB this._p := by
    intro k value hk hread'
    rcases heap.readU64_of_valid B lenB k hB hk with ⟨old, hread⟩
    have hvalue : value = old :=
      Except.ok.inj (hread'.symm.trans (hsameB k old hk hread))
    subst value
    exact hcanonicalB k old hk hread
  have hdivisorResult := slicePolyRep_of_same_prefix heap heap' B lenB
    this._p.toNat divisor hB hBResult hsameB hdivisor
  have hnormBEq := normaliseU64_eq_of_prefix_map heap heap' B B lenB hB
    hsameB
  have hnormBResult : heap'.normaliseU64 B lenB = .ok lenB := by
    rw [← hnormBEq]
    exact hnormB
  have hRFull : heap'.ValidU64Slice R (Nat.min lenA (lenB - 1)) :=
    (hlayout R (Nat.min lenA (lenB - 1))).mp hR
  have hRResult : heap'.ValidU64Slice R lenR :=
    heap'.validU64Slice_mono R (Nat.min lenA (lenB - 1)) lenR hRFull
      hlenRCapacity
  exact ⟨heap', lenQ, lenR, quotient, remainder, hrun, hquotient,
    ⟨hBResult, hcanonicalBResult, hdivisorResult, hnormBResult⟩,
    ⟨hRResult, hcanonicalR, hremainder, hnormR⟩, hdivision, hgcd, hlayout,
    hlenQ, hlenRCapacity, hlenR⟩

/-- End-to-end semantic refinement of the actual well-founded raw Euclid
loop.  The proof follows every generated `_poly_divrem` call and every source
buffer rotation; no L2 loop, assumed remainder, alternate path, or bounded execution
counter is used. -/
theorem strictEuclidLoop_refines (this : DenseUPolyZp)
    [hprime : Fact (Nat.Prime this._p.toNat)]
    (Q : RawPtr UInt64) (W3 : RawPtr Word3) (capacity : Nat) :
    (heap : RawHeap) → (A : RawPtr UInt64) → (lenA : Nat) →
      (B : RawPtr UInt64) → (lenB : Nat) → (R : RawPtr UInt64) →
      (dividend divisor : Polynomial (ZMod this._p.toNat)) →
      EuclidWorkspace heap Q W3 A B R capacity →
      lenA ≤ capacity → lenB ≤ capacity → capacity < limbBase →
      RawDensePolyRep this heap A lenA dividend →
      RawDensePolyRep this heap B lenB divisor →
      DensePreinvConfigured this →
      ∃ heap' out outLen result,
        strictEuclidLoop this Q W3 heap A lenA B lenB R =
          .ok (heap', out, outLen) ∧
        RawDensePolyRep this heap' out outLen result ∧
        normalize (EuclideanDomain.gcd dividend divisor) = normalize result ∧
        RawHeap.SameLayout heap heap'
  | heap, A, lenA, B, 0, R, dividend, divisor, hworkspace, hlenA, _, _,
      hARep, hBRep, _ => by
    have hdivisorZero : divisor = 0 :=
      slicePolyRep_zero_length heap B this._p.toNat divisor hBRep.2.2.1
    subst divisor
    refine ⟨heap, A, lenA, dividend, ?_, hARep, ?_, fun _ _ => Iff.rfl⟩
    · simp [strictEuclidLoop, Generated.StrictEuclidGCD.euclidLoop]
    · rw [EuclideanDomain.gcd_zero_right]
  | heap, A, lenA, B, lenB + 1, R, dividend, divisor, hworkspace, hlenA,
      hlenB, hcapacity, hARep, hBRep, hcfg => by
    have hQ : heap.ValidU64Slice Q (lenA - ((lenB + 1) - 1)) :=
      heap.validU64Slice_mono Q capacity _ hworkspace.validQ (by omega)
    have hR : heap.ValidU64Slice R (Nat.min lenA ((lenB + 1) - 1)) :=
      heap.validU64Slice_mono R capacity _ hworkspace.validR
        ((Nat.min_le_left lenA ((lenB + 1) - 1)).trans hlenA)
    have hW3 : heap.ValidWord3Slice W3 lenA :=
      heap.validU64Slice_mono (RawPtr.reinterpret W3) (3 * capacity)
        (3 * lenA) hworkspace.validW3 (by omega)
    have hqCapacity : lenA - ((lenB + 1) - 1) < limbBase := by omega
    rcases polyDivrem_next_state this Q R A B lenA (lenB + 1) W3 heap
        dividend divisor (by omega) hARep hBRep hQ hR hW3 hqCapacity
        hworkspace.regions.a_r.symm hworkspace.regions.w3_a
        hworkspace.regions.w3_b hworkspace.regions.q_b
        hworkspace.regions.q_w3 hworkspace.regions.w3_r.symm
        hworkspace.regions.q_r.symm hworkspace.regions.b_r.symm hcfg with
      ⟨heap1, lenQ, lenR, quotient, remainder, hrun, _, hBRep1, hRRep1,
        _, hgcdStep, hlayout1, _, _, hlenR⟩
    have hworkspace1 := hworkspace.rotate_of_sameLayout hlayout1
    rcases strictEuclidLoop_refines this Q W3 capacity heap1 B (lenB + 1) R
        lenR A divisor remainder hworkspace1 hlenB (Nat.le_trans
          (Nat.le_of_lt hlenR) hlenB) hcapacity hBRep1 hRRep1 hcfg with
      ⟨heap2, out, outLen, result, hloop, hresult, hgcdRest, hlayout2⟩
    refine ⟨heap2, out, outLen, result, ?_, hresult,
      hgcdStep.trans hgcdRest, ?_⟩
    · unfold strictEuclidLoop
      rw [Generated.StrictEuclidGCD.euclidLoop, hrun]
      simpa [strictEuclidLoop] using hloop
    · intro ptr length
      exact (hlayout1 ptr length).trans (hlayout2 ptr length)
termination_by heap A lenA B lenB R dividend divisor => lenB
decreasing_by exact hlenR

/-- The final live buffer of the generated Euclid rotation is one of its
three physical polynomial work regions, and its logical length never exceeds
their common capacity.  This follows the source pointer rotation itself. -/
theorem euclidLoop_output_mem_and_length (this : DenseUPolyZp)
    (Q : RawPtr UInt64) (W3 : RawPtr Word3)
    (hdec : DivremLengthDecreases this Q W3) (capacity : Nat) :
    (heap : RawHeap) → (A : RawPtr UInt64) → (lenA : Nat) →
      (B : RawPtr UInt64) → (lenB : Nat) → (R : RawPtr UInt64) →
      (heap' : RawHeap) → (out : RawPtr UInt64) → (outLen : Nat) →
      lenA ≤ capacity → lenB ≤ capacity →
      euclidLoop this Q W3 hdec heap A lenA B lenB R =
        .ok (heap', out, outLen) →
      (out = A ∨ out = B ∨ out = R) ∧ outLen ≤ capacity
  | heap, A, lenA, B, 0, R, heap', out, outLen, hlenA, _, hrun => by
      simp only [euclidLoop] at hrun
      have heq : (heap', out, outLen) = (heap, A, lenA) :=
        Except.ok.inj hrun.symm
      cases heq
      exact ⟨Or.inl rfl, hlenA⟩
  | heap, A, lenA, B, lenB + 1, R, heap', out, outLen, hlenA, hlenB,
      hrun => by
      rw [euclidLoop] at hrun
      generalize hcall : dense_upoly_zp__poly_divrem_ir this Q R A lenA B
        (lenB + 1) W3 heap = call at hrun
      cases call with
      | error fault => simp at hrun
      | ok result =>
          rcases result with ⟨heap1, lenQ, lenR⟩
          simp only at hrun
          have hlenR : lenR < lenB + 1 :=
            hdec heap A B R lenA (lenB + 1) heap1 lenQ lenR (by omega) hcall
          rcases euclidLoop_output_mem_and_length this Q W3 hdec capacity
              heap1 B (lenB + 1) R lenR A heap' out outLen hlenB
              (by omega) hrun with ⟨hloc, houtLen⟩
          exact ⟨by rcases hloc with hB | hR | hA
                    · exact Or.inr (Or.inl hB)
                    · exact Or.inr (Or.inr hR)
                    · exact Or.inl hA,
            houtLen⟩
termination_by heap A lenA B lenB R heap' out outLen => lenB
decreasing_by exact hlenR

/-- Physical allocation contract for the complete C++ `_gcd_euclid`
helper, including its two local input copies and final copy to `G`. -/
structure GcdEuclidRawWorkspace (heap : RawHeap)
    (G A B aBuf bBuf Q R : RawPtr UInt64) (W3 : RawPtr Word3)
    (capacity : Nat) : Prop where
  validG : heap.ValidU64Slice G capacity
  euclid : EuclidWorkspace heap Q W3 aBuf bBuf R capacity
  aBuf_A : aBuf.region ≠ A.region
  aBuf_B : aBuf.region ≠ B.region
  bBuf_B : bBuf.region ≠ B.region
  bBuf_aBuf : bBuf.region ≠ aBuf.region
  G_aBuf : G.region ≠ aBuf.region
  G_bBuf : G.region ≠ bBuf.region
  G_R : G.region ≠ R.region

/-- Proof-instantiated execution of the exact generated raw helper. -/
def strictGcdEuclidRaw (this : DenseUPolyZp)
    (G A B aBuf bBuf Q R : RawPtr UInt64) (W3 : RawPtr Word3)
    (lenA lenB : Nat) (heap : RawHeap) : RawExec GcdEuclidRawResult :=
  dense_upoly_zp__gcd_euclid_raw_ir this G A B aBuf bBuf Q R W3 lenA lenB
    (euclidDivremLengthDecreases this Q W3) heap

/-- End-to-end L1-to-L2 refinement of the actual C++ `_gcd_euclid` raw
execution.  Every copy and every division step is executed in `RawHeap`; the
L2 gcd appears only in the postcondition. -/
theorem strictGcdEuclidRaw_refines (this : DenseUPolyZp)
    [hprime : Fact (Nat.Prime this._p.toNat)]
    (G A B aBuf bBuf Q R : RawPtr UInt64) (W3 : RawPtr Word3)
    (lenA lenB capacity : Nat) (heap : RawHeap)
    (left right : Polynomial (ZMod this._p.toNat))
    (hworkspace : GcdEuclidRawWorkspace heap G A B aBuf bBuf Q R W3 capacity)
    (hlenA : lenA ≤ capacity) (hlenB : lenB ≤ capacity)
    (hcapacity : capacity < limbBase)
    (hleft : RawDensePolyRep this heap A lenA left)
    (hright : RawDensePolyRep this heap B lenB right)
    (hcfg : DensePreinvConfigured this) :
    ∃ heap' lenG result,
      strictGcdEuclidRaw this G A B aBuf bBuf Q R W3 lenA lenB heap =
        .ok ⟨heap', lenG⟩ ∧
      RawDensePolyRep this heap' G lenG result ∧
      normalize (EuclideanDomain.gcd left right) = normalize result ∧
      RawHeap.SameLayout heap heap' := by
  have haValid := heap.validU64Slice_mono aBuf capacity lenA
    hworkspace.euclid.validA hlenA
  rcases copyU64_refines_rawDense_of_region_ne this heap aBuf A lenA left
      haValid hworkspace.aBuf_A hleft with
    ⟨heap1, hcopyA, hlayout1, hleft1⟩
  have hright1 := (copyU64_preserves_rawDense_of_region_ne this heap heap1
    aBuf A lenA B lenB right haValid hleft.1 hworkspace.aBuf_B hcopyA
    hright).2
  have hworkspace1 := hworkspace.euclid.of_sameLayout hlayout1
  have hbValid := heap1.validU64Slice_mono bBuf capacity lenB
    hworkspace1.validB hlenB
  rcases copyU64_refines_rawDense_of_region_ne this heap1 bBuf B lenB right
      hbValid hworkspace.bBuf_B hright1 with
    ⟨heap2, hcopyB, hlayout2, hright2⟩
  have hleft2 := (copyU64_preserves_rawDense_of_region_ne this heap1 heap2
    bBuf B lenB aBuf lenA left hbValid hright1.1 hworkspace.bBuf_aBuf
    hcopyB hleft1).2
  have hworkspace2 := hworkspace1.of_sameLayout hlayout2
  rcases strictEuclidLoop_refines this Q W3 capacity heap2 aBuf lenA bBuf
      lenB R left right hworkspace2 hlenA hlenB hcapacity hleft2 hright2
      hcfg with
    ⟨heap3, out, outLen, result, hloop, hresult, hgcd, hlayout3⟩
  have hout := euclidLoop_output_mem_and_length this Q W3
    (euclidDivremLengthDecreases this Q W3) capacity heap2 aBuf lenA bBuf
    lenB R heap3 out outLen hlenA hlenB hloop
  have hGValid0 := heap.validU64Slice_mono G capacity outLen
    hworkspace.validG hout.2
  have hGValid1 := (hlayout1 G outLen).mp hGValid0
  have hGValid2 := (hlayout2 G outLen).mp hGValid1
  have hGValid3 := (hlayout3 G outLen).mp hGValid2
  have hGOut : G.region ≠ out.region := by
    rcases hout.1 with rfl | rfl | rfl
    · exact hworkspace.G_aBuf
    · exact hworkspace.G_bBuf
    · exact hworkspace.G_R
  rcases copyU64_refines_rawDense_of_region_ne this heap3 G out outLen result
      hGValid3 hGOut hresult with
    ⟨heap4, hcopyG, hlayout4, hresultG⟩
  refine ⟨heap4, outLen, result, ?_, hresultG, hgcd, ?_⟩
  · unfold strictGcdEuclidRaw dense_upoly_zp__gcd_euclid_raw_ir
    simp only [hcopyA]
    simp only [hcopyB]
    change euclidLoop this Q W3 (euclidDivremLengthDecreases this Q W3)
      heap2 aBuf lenA bBuf lenB R = .ok (heap3, out, outLen) at hloop
    simp only [hloop, hcopyG]
  · intro ptr length
    exact (hlayout1 ptr length).trans ((hlayout2 ptr length).trans
      ((hlayout3 ptr length).trans (hlayout4 ptr length)))

end CLPoly.Impl.StrictEuclidRefinement
