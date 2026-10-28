import CLPoly.Refinement.VanHoeijSeparationInput
import CLPoly.Refinement.VanHoeijWidth
import CLPoly.Refinement.VanHoeijColumnCoverage

set_option autoImplicit false

namespace Refinement.StrictRecombine

/-- The same CLD execution bounds every physical head before the source's
Int32 width scan, so its result equals the natural full coefficient width. -/
theorem CldSeparationData.truncated_width_toInt
    {f : SparsePolyZZ} {active original : Array SparsePolyZZ} {modulus : Nat}
    (data : CldSeparationData f active original modulus)
    (truncated : Array SparsePolyZZ) (divisor : ZZ)
    (htruncate : Generated.StrictVanHoeij.truncateCld original divisor = .ok truncated)
    (hfits : (SparsePolyZZ.toPoly f).natDegree < 2 ^ 31) :
    (Generated.StrictVanHoeijParameters.cldWidth truncated).toInt =
      (Generated.StrictRecombine.cldSpiralWidth truncated : Int) := by
  have hw := truncateCld_width_le original truncated divisor htruncate
    (fun i hi => (data.values i hi).canonical) _ data.degrees
  apply (cldWidth_eq_natural truncated ?_).1
  intro i hi he
  have hhead := cldSpiralWidth_covers truncated i hi
  rw [getElem!_pos truncated i hi, he] at hhead
  simp only [Bool.false_eq_true, ↓reduceIte] at hhead
  omega

end Refinement.StrictRecombine
