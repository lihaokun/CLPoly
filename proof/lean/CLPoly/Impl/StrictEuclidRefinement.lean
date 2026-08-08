import CLPoly.Impl.StrictDivremRefinement

set_option autoImplicit false

namespace CLPoly.Impl.StrictEuclidRefinement

open Generated.StrictDivrem
open CLPoly.Impl.RawPolynomialRep
open CLPoly.Impl.StrictDivremRefinement
open CLPoly.Impl.StrictWordArithmetic

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
    (hcfg : DensePreinvConfigured this) :
    ∃ heap' lenQ lenR quotient remainder,
      dense_upoly_zp__poly_divrem_ir this Q R A lenA B lenB W3 heap =
        .ok (heap', lenQ, lenR) ∧
      SlicePolyRep heap' Q lenQ this._p.toNat quotient ∧
      SlicePolyRep heap' R lenR this._p.toNat remainder ∧
      dividend = quotient * divisor + remainder ∧
      normalize (EuclideanDomain.gcd dividend divisor) =
        normalize (EuclideanDomain.gcd divisor remainder) ∧
      (remainder = 0 ∨ remainder.natDegree < divisor.natDegree) ∧
      lenQ ≤ lenA - (lenB - 1) ∧ lenR < lenB := by
  rcases polyDivrem_refines this Q R A B lenA lenB W3 heap dividend divisor
      hlenB hA hB hQ hR hW3 hcanonicalA hcanonicalB hdividend hdivisor
      hnormA hnormB hqCapacity hRA hWA hWB hQB hQW hRW hRQ hcfg
      hprime.out with
    ⟨heap', lenQ, lenR, quotient, remainder, hrun, hquotient, hremainder,
      hdivision, hdegree, hlenQ, hlenR⟩
  exact ⟨heap', lenQ, lenR, quotient, remainder, hrun, hquotient,
    hremainder, hdivision,
    normalize_gcd_eq_of_division_identity dividend divisor quotient remainder
      hdivision,
    hdegree, hlenQ, hlenR⟩

end CLPoly.Impl.StrictEuclidRefinement
