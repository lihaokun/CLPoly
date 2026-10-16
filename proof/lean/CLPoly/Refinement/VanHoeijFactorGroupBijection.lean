import CLPoly.Refinement.VanHoeijFactorFamilyCoverage

set_option autoImplicit false

namespace Refinement.StrictRecombine

/-- The individual exact groups cover all actual output-class indices, with
the original factor-array indices retained as witnesses. -/
theorem factorGroupAssignment_surjective
    (source : Polynomial Int) (factors : Array (Polynomial Int))
    (hproduct : factors.toList.prod = source)
    (prime modulus : Nat) [Fact (Nat.Prime prime)] (hdivides : prime ∣ modulus)
    (lifted : Array SparsePolyZZ)
    (transform : Generated.StrictRecombine.LLLMatrix) (rows : Array Nat)
    (groups : Array (Array Int32)) (hfits : lifted.size ≤ 2 ^ 31)
    (hrun : Generated.StrictRecombine.extractCandidates rows transform lifted.size = .ok groups)
    (hprime : ∀ i : Fin lifted.size, Prime (StrictHensel.toPolyMod prime lifted[i.val]))
    (hseparate : ∀ i j : Fin lifted.size,
      Associated (StrictHensel.toPolyMod prime lifted[i.val])
        (StrictHensel.toPolyMod prime lifted[j.val]) → i = j)
    (hatom : ∀ i : Fin lifted.size, StrictHensel.toPolyMod prime lifted[i.val] ∣
      source.map (Int.castRingHom (ZMod prime)))
    (hleading : ∀ i : Fin factors.size, (factors[i.val].leadingCoeff : ZMod prime) ≠ 0)
    (assigned : Fin factors.size → Fin groups.size)
    (indices : Fin factors.size → Array Nat)
    (hlegal : ∀ i, LegalCombination lifted.size (indices i).size (indices i))
    (hsupport : ∀ i, factors[i.val].map (Int.castRingHom (ZMod modulus)) =
      Polynomial.C (factors[i.val].leadingCoeff : ZMod modulus) *
        ((indices i).toList.map fun j => StrictHensel.toPolyMod modulus lifted[j]!).prod)
    (hmembers : ∀ i index, index < lifted.size →
      (CandidateMember groups (assigned i).val index.toUInt32.toInt32 ↔ index ∈ (indices i).toList)) :
    Function.Surjective assigned := by
  intro classId
  obtain ⟨sourceIndex, hi, hm⟩ := extractCandidates_every_group_has_source
    transform rows lifted.size groups hrun classId.val classId.isLt
  let i : Fin lifted.size := ⟨sourceIndex, hi⟩
  let atoms := fun j : Fin lifted.size => StrictHensel.toPolyMod prime lifted[j.val]
  have hd : atoms i ∣ (factors.toList.map (Polynomial.map (Int.castRingHom (ZMod prime)))).prod := by
    rw [← Polynomial.map_list_prod, hproduct]
    exact hatom i
  obtain ⟨reduced, hmem, hdiv⟩ := prime_dvd_list_product_witness (atoms i) (hprime i) _ hd
  obtain ⟨factor, hf, he⟩ := List.mem_map.mp hmem
  obtain ⟨position, hp, hv⟩ := List.mem_iff_getElem.mp hf
  have hp' : position < factors.size := by simpa using hp
  have hv' : factors[position] = factor := by simpa only [Array.getElem_toList hp'] using hv
  let j : Fin factors.size := ⟨position, hp'⟩
  have hdiv' : atoms i ∣ factors[j.val].map (Int.castRingHom (ZMod prime)) := by
    simpa only [j, hv', he] using hdiv
  have ha := divisor_support_project_associated prime modulus hdivides lifted (indices j)
    (hlegal j) factors[j.val] (hleading j) (hsupport j)
  have hin : i ∈ legalPhysicalSupport lifted.size (indices j) (hlegal j) :=
    (prime_atom_dvd_support_iff atoms hprime hseparate _ i).mp (hdiv'.trans ha.dvd)
  have hin' : sourceIndex ∈ (indices j).toList :=
    (legalPhysicalSupport_mem lifted.size (indices j) (hlegal j) i).mp hin
  have hother := (hmembers j sourceIndex hi).mpr hin'
  refine ⟨j, Fin.ext ?_⟩
  exact extractCandidates_class_unique transform rows lifted.size groups hfits hrun
    sourceIndex (assigned j).val classId.val hi hother hm

/-- Distinct mathematical factors cannot be assigned the same exact group
once the explicit separation threshold is reached. -/
theorem factorGroupAssignment_injective
    (source : Polynomial Int) (factors : Array (Polynomial Int))
    (hproduct : factors.toList.prod = source) (hnodup : factors.toList.Nodup)
    (hf : source ≠ 0) (hn : 0 < source.natDegree)
    (hprimitive : ∀ i : Fin factors.size, factors[i.val].IsPrimitive)
    (hpositive : ∀ i : Fin factors.size, (0 : Int) < factors[i.val].leadingCoeff)
    (modulus : Nat) (hm : 0 < modulus) (lifted : Array SparsePolyZZ)
    (groups : Array (Array Int32)) (divisor : Int) (hd : 0 ≤ divisor)
    (hprecision : separationResultantThreshold source lifted.size divisor < (modulus : Int))
    (assigned : Fin factors.size → Fin groups.size)
    (hsupport : ∀ i, factors[i.val].map (Int.castRingHom (ZMod modulus)) =
      Polynomial.C (factors[i.val].leadingCoeff : ZMod modulus) *
        SelectedProductMod modulus groups[(assigned i).val]! lifted 0) :
    Function.Injective assigned := by
  intro i j he
  have hdi : factors[i.val] ∣ source := by
    rw [← hproduct]; exact List.dvd_prod (Array.getElem_mem_toList i.isLt)
  have hdj : factors[j.val] ∣ source := by
    rw [← hproduct]; exact List.dvd_prod (Array.getElem_mem_toList j.isLt)
  have heq : factors[i.val] = factors[j.val] :=
    positive_primitive_divisor_unique_of_support source _ _ lifted.size divisor hf hn hd hdi hdj
      (hprimitive i) (hprimitive j) (hpositive i) (hpositive j) modulus hm hprecision
      (SelectedProductMod modulus groups[(assigned i).val]! lifted 0) (hsupport i)
      (by simpa only [he] using hsupport j)
  apply Fin.ext
  exact hnodup.getElem_inj_iff.mp (by simpa only [Array.getElem_toList] using heq)

end Refinement.StrictRecombine
