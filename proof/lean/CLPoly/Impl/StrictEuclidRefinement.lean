import CLPoly.Impl.StrictDivremRefinement

set_option autoImplicit false

namespace CLPoly.Impl.StrictEuclidRefinement

open Generated.StrictDivrem
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
      hcanonicalR, hnormR, hlayout, hsameB, _, hgcd, _, hlenQ,
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
    ⟨hRResult, hcanonicalR, hremainder, hnormR⟩, hgcd, hlayout,
    hlenQ, hlenRCapacity, hlenR⟩

/-- End-to-end semantic refinement of the actual well-founded raw Euclid
loop.  The proof follows every generated `_poly_divrem` call and every source
buffer rotation; no L2 loop, remainder oracle, fallback, or bounded execution
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
        hgcdStep, hlayout1, _, _, hlenR⟩
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

end CLPoly.Impl.StrictEuclidRefinement
