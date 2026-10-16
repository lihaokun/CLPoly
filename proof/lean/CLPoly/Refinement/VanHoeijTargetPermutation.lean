import CLPoly.Refinement.VanHoeijFactorGroupBijection
import CLPoly.Refinement.VanHoeijPartitionValidationComplete
import Mathlib.Algebra.BigOperators.Fin

set_option autoImplicit false

namespace Refinement.StrictRecombine

/-- Reindexing a mathematical factor family along the actual group bijection
constructs all validation targets and preserves their product. -/
theorem exactPartitionTargets_of_bijective_assignment
    (source : Polynomial Int) (factors : Array (Polynomial Int))
    (hproduct : factors.toList.prod = source)
    (hprimitive : ∀ i : Fin factors.size, factors[i.val].IsPrimitive)
    (hpositive : ∀ i : Fin factors.size, (0 : Int) < factors[i.val].leadingCoeff)
    (hdegree : ∀ i : Fin factors.size, 0 < factors[i.val].natDegree)
    (modulus : Nat) (lifted : Array SparsePolyZZ) (groups : Array (Array Int32))
    (assigned : Fin factors.size → Fin groups.size) (hbijective : Function.Bijective assigned)
    (hsupport : ∀ i, factors[i.val].map (Int.castRingHom (ZMod modulus)) =
      Polynomial.C (factors[i.val].leadingCoeff : ZMod modulus) *
        SelectedProductMod modulus groups[(assigned i).val]! lifted 0) :
    ∃ targets, ExactPartitionTargets source modulus lifted groups targets := by
  classical
  let e : Fin factors.size ≃ Fin groups.size := Equiv.ofBijective assigned hbijective
  let targets := Array.ofFn fun i : Fin groups.size => factors[(e.symm i).val]
  have hsize : targets.size = groups.size := Array.size_ofFn
  have hvalue : ∀ i (hi : i < groups.size), targets[i]! = factors[(e.symm ⟨i, hi⟩).val] := by
    intro i hi
    simp [targets, getElem!_pos, hi]
  have hmember : ∀ factor ∈ targets.toList, ∃ i : Fin factors.size, factor = factors[i.val] := by
    intro factor hf
    have h : factor ∈ List.ofFn (fun i : Fin groups.size => factors[(e.symm i).val]) := by
      simpa only [targets, Array.toList_ofFn] using hf
    obtain ⟨i, hi⟩ := List.mem_ofFn.mp h
    exact ⟨e.symm i, hi.symm⟩
  refine ⟨targets, ⟨hsize, ?_, ?_, ?_, ?_, ?_⟩⟩
  · have hlist : List.ofFn (fun i : Fin factors.size => factors[i.val]) = factors.toList := by
      simpa only [Array.toList_ofFn] using
        congrArg Array.toList (Array.ofFn_getElem (xs := factors))
    calc
      targets.toList.prod = ∏ i : Fin groups.size, factors[(e.symm i).val] := by
        rw [show targets.toList = List.ofFn (fun i : Fin groups.size => factors[(e.symm i).val])
          from Array.toList_ofFn, List.prod_ofFn]
      _ = ∏ i : Fin factors.size, factors[i.val] := e.symm.prod_comp _
      _ = factors.toList.prod := by rw [← List.prod_ofFn, hlist]
      _ = source := hproduct
  · intro factor hf
    obtain ⟨i, rfl⟩ := hmember factor hf
    exact hprimitive i
  · intro factor hf
    obtain ⟨i, rfl⟩ := hmember factor hf
    exact hpositive i
  · intro factor hf
    obtain ⟨i, rfl⟩ := hmember factor hf
    exact hdegree i
  · intro i hi
    rw [hvalue i hi]
    have he : assigned (e.symm ⟨i, hi⟩) = ⟨i, hi⟩ := e.apply_symm_apply _
    simpa only [he] using hsupport (e.symm ⟨i, hi⟩)

end Refinement.StrictRecombine
