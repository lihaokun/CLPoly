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

/-- Coefficient-level characterization of the raw-array interpretation. -/
theorem coeff_coeffArrayPoly (p : Nat) (coeffs : Array UInt64) (degree : Nat) :
    (coeffArrayPoly p coeffs).coeff degree =
      if h : degree < coeffs.size then
        (coeffs[degree].toNat : ZMod p)
      else 0 := by
  classical
  unfold coeffArrayPoly
  rw [Polynomial.finset_sum_coeff]
  by_cases hdegree : degree < coeffs.size
  · rw [dif_pos hdegree, Finset.sum_eq_single ⟨degree, hdegree⟩]
    · simp
    · intro index _ hne
      have hval : index.val ≠ degree := by
        intro heq
        apply hne
        exact Fin.ext heq
      simp [Polynomial.coeff_monomial, hval]
    · simp
  · rw [dif_neg hdegree]
    apply Finset.sum_eq_zero
    intro index _
    have hval : index.val ≠ degree := by
      intro heq
      apply hdegree
      simpa [heq] using index.isLt
    simp [Polynomial.coeff_monomial, hval]

/-- A raw slice represents an L2 polynomial iff its complete, failure-aware
observation has the declared C++ length and coefficient interpretation. -/
def SlicePolyRep (heap : RawHeap) (ptr : RawPtr UInt64) (length p : Nat)
    (poly : Polynomial (ZMod p)) : Prop :=
  ∃ coeffs : Array UInt64,
    heap.SliceRep ptr length coeffs ∧ coeffs.size = length ∧
    poly = coeffArrayPoly p coeffs

end CLPoly.Impl.RawPolynomialRep
