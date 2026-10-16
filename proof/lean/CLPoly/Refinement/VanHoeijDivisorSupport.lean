import CLPoly.Refinement.VanHoeijNoSplit

set_option autoImplicit false

namespace Refinement.StrictRecombine

def legalPhysicalSupport (count : Nat) (indices : Array Nat)
    (hlegal : LegalCombination count indices.size indices) : List (Fin count) :=
  List.ofFn fun i : Fin indices.size => ⟨indices[i.val],
    by simpa only [getElem!_pos indices i.val i.isLt] using hlegal.2.2 i.val i.isLt⟩

theorem legalPhysicalSupport_values (count : Nat) (indices : Array Nat)
    (hlegal : LegalCombination count indices.size indices) :
    (legalPhysicalSupport count indices hlegal).map Fin.val = indices.toList := by
  unfold legalPhysicalSupport
  rw [List.map_ofFn]
  simpa only [Array.toList_ofFn] using
    congrArg Array.toList (Array.ofFn_getElem (xs := indices))

theorem legalPhysicalSupport_mem (count : Nat) (indices : Array Nat)
    (hlegal : LegalCombination count indices.size indices) (i : Fin count) :
    i ∈ legalPhysicalSupport count indices hlegal ↔ i.val ∈ indices.toList := by
  rw [← legalPhysicalSupport_values count indices hlegal, List.mem_map]
  constructor
  · intro hi
    exact ⟨i, hi, rfl⟩
  · rintro ⟨j, hj, he⟩
    have hij : j = i := Fin.ext he
    simpa only [hij] using hj

theorem legalPhysicalSupport_product
    (prime : Nat) (lifted : Array SparsePolyZZ) (indices : Array Nat)
    (hlegal : LegalCombination lifted.size indices.size indices) :
    ((legalPhysicalSupport lifted.size indices hlegal).map
      (fun i => StrictHensel.toPolyMod prime lifted[i.val])).prod =
      (indices.toList.map fun i => StrictHensel.toPolyMod prime lifted[i]!).prod := by
  rw [← legalPhysicalSupport_values lifted.size indices hlegal, List.map_map]
  congr 1
  apply List.map_congr_left
  intro i hi
  simp only [Function.comp_apply, getElem!_pos lifted i.val i.isLt]

theorem divisor_support_project_associated
    (prime modulus : Nat) [Fact (Nat.Prime prime)] (hd : prime ∣ modulus)
    (lifted : Array SparsePolyZZ) (indices : Array Nat)
    (hlegal : LegalCombination lifted.size indices.size indices)
    (factor : Polynomial Int)
    (hleading : (factor.leadingCoeff : ZMod prime) ≠ 0)
    (hsupport : factor.map (Int.castRingHom (ZMod modulus)) =
      Polynomial.C (factor.leadingCoeff : ZMod modulus) *
        (indices.toList.map fun i => StrictHensel.toPolyMod modulus lifted[i]!).prod) :
    Associated (factor.map (Int.castRingHom (ZMod prime)))
      (((legalPhysicalSupport lifted.size indices hlegal).map
        (fun i => StrictHensel.toPolyMod prime lifted[i.val])).prod) := by
  let product := (indices.toList.map fun i => SparsePolyZZ.toPoly lifted[i]!).prod
  have he : factor.map (Int.castRingHom (ZMod modulus)) =
      (Polynomial.C factor.leadingCoeff * product).map (Int.castRingHom (ZMod modulus)) := by
    simpa only [product, Polynomial.map_mul, Polynomial.map_C, Polynomial.map_list_prod,
      List.map_map, Function.comp_def, StrictHensel.toPolyMod] using hsupport
  have hp := polynomialMap_eq_of_modulus_dvd prime modulus hd factor
    (Polynomial.C factor.leadingCoeff * product) he
  have hp' : factor.map (Int.castRingHom (ZMod prime)) =
      Polynomial.C (factor.leadingCoeff : ZMod prime) *
        (indices.toList.map fun i => StrictHensel.toPolyMod prime lifted[i]!).prod := by
    simpa only [product, Polynomial.map_mul, Polynomial.map_C, Polynomial.map_list_prod,
      List.map_map, Function.comp_def, StrictHensel.toPolyMod] using hp
  rw [legalPhysicalSupport_product, hp']
  exact (associated_isUnit_mul_left_iff
    (Polynomial.isUnit_C.mpr (isUnit_iff_ne_zero.mpr hleading))).mpr (Associated.refl _)

end Refinement.StrictRecombine
