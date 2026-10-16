import CLPoly.Refinement.VanHoeijAssembly

set_option autoImplicit false

namespace Refinement.StrictRecombine

/-- Actual same-class members have the same membership in a preserved
support. The nonempty-row premise is derived from the physical member. -/
theorem lll_class_support_membership
    (matrix : Generated.StrictRecombine.LLLMatrix)
    (hinput : ConcreteLLLInputValid matrix) (bound : ZZ) (factorCount : Nat)
    (hfactorCount : factorCount ≤ matrix.size)
    (output : Generated.StrictRecombine.LLLMatrix ×
      Generated.StrictRecombine.LLLMatrix × Array Nat × Array QQ)
    (hrun : Generated.StrictVanHoeij.lllReduceWithNorms concreteLLLExecution
      matrix hinput bound = .ok output)
    (indices : Array Nat) (coefficients : Fin matrix.size → ZZ)
    (hcoordinates : ∀ index, index.val < factorCount →
      coefficients index = if index.val ∈ indices.toList then 1 else 0)
    (hshort : dotProduct
      (fun column => ((Matrix.vecMul coefficients
        (basisPrefixMatrix matrix matrix.size matrix.size) column : ZZ) : QQ))
      (fun column => ((Matrix.vecMul coefficients
        (basisPrefixMatrix matrix matrix.size matrix.size) column : ZZ) : QQ)) ≤ (bound : QQ))
    (candidates : Array (Array Int32))
    (hextract : Generated.StrictRecombine.extractCandidates
      (Generated.StrictVanHoeij.shortSpanFromNorms output.2.2.2 bound)
      output.2.1 factorCount = .ok candidates)
    (classId : Nat) (leftValue rightValue : Int32)
    (hleft : CandidateMember candidates classId leftValue)
    (hright : CandidateMember candidates classId rightValue) :
    ∃ (left right : Nat), left < factorCount ∧ right < factorCount ∧
      leftValue = left.toUInt32.toInt32 ∧ rightValue = right.toUInt32.toInt32 ∧
      (left ∈ indices.toList ↔ right ∈ indices.toList) := by
  have hnonempty : (Generated.StrictVanHoeij.shortSpanFromNorms output.2.2.2 bound).isEmpty = false := by
    by_contra hn
    have he : (Generated.StrictVanHoeij.shortSpanFromNorms output.2.2.2 bound).isEmpty = true := by
      cases h : (Generated.StrictVanHoeij.shortSpanFromNorms output.2.2.2 bound).isEmpty <;> simp_all
    simp only [Generated.StrictRecombine.extractCandidates, he, ↓reduceIte] at hextract
    have hc := Except.ok.inj hextract
    have hb := hleft.1
    rw [← hc] at hb
    simp at hb
  obtain ⟨left, right, hl, hr, hlv, hrv, heq⟩ :=
    lllReduceWithNorms_class_members_coordinates matrix hinput bound factorCount
      hfactorCount output hrun coefficients hshort candidates hnonempty hextract
      classId leftValue rightValue hleft hright
  rw [hcoordinates _ hl, hcoordinates _ hr] at heq
  refine ⟨left, right, hl, hr, hlv, hrv, ?_⟩
  by_cases hlm : left ∈ indices.toList <;>
    by_cases hrm : right ∈ indices.toList <;> simp_all

end Refinement.StrictRecombine
