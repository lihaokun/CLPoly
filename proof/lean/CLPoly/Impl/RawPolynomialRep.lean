import CLPoly.Generated.StrictDivrem
import CLPoly.Spec

namespace CLPoly.Impl.RawPolynomialRep

open Polynomial

/-- Interpret a C++ little-endian coefficient buffer as the same L2 object
used by SQF, DDF and EDF: `Polynomial (ZMod p)`. -/
noncomputable def coeffArrayPoly (p : Nat) (coeffs : Array UInt64) :
    Polynomial (ZMod p) :=
  ∑ i : Fin coeffs.size,
    Polynomial.monomial i.val (coeffs[i].toNat : ZMod p)

/-- A raw slice represents an L2 polynomial iff its complete, failure-aware
observation has the declared C++ length and coefficient interpretation. -/
def SlicePolyRep (heap : RawHeap) (ptr : RawPtr UInt64) (length p : Nat)
    (poly : Polynomial (ZMod p)) : Prop :=
  ∃ coeffs : Array UInt64,
    heap.SliceRep ptr length coeffs ∧ coeffs.size = length ∧
    poly = coeffArrayPoly p coeffs

end CLPoly.Impl.RawPolynomialRep
