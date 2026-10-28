import CLPoly.Refinement.VanHoeijSupportUnique
import CLPoly.Refinement.VanHoeijDivisorSupport

set_option autoImplicit false

namespace Refinement.StrictRecombine

theorem prime_dvd_list_product_witness
    {R : Type*} [CommMonoidWithZero R] (atom : R) (hprime : Prime atom)
    (values : List R) (hd : atom ∣ values.prod) :
    ∃ value ∈ values, atom ∣ value := by
  induction values with
  | nil => exact (hprime.not_unit (isUnit_of_dvd_one (by simpa using hd))).elim
  | cons head tail ih =>
    rw [List.prod_cons] at hd
    rcases hprime.dvd_or_dvd hd with hh | ht
    · exact ⟨head, List.mem_cons_self, hh⟩
    · obtain ⟨value, hv, hdv⟩ := ih ht
      exact ⟨value, List.mem_cons_of_mem head hv, hdv⟩

/-- A mathematical factor family and its individual exact groups cover every
actual returned class. No factor correspondence for an arbitrary output class
is assumed, and no validation-success premise is used. -/
theorem extractCandidates_every_group_has_factor
    (source : Polynomial Int) (factors : List (Polynomial Int))
    (hproduct : factors.prod = source)
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
    (hleading : ∀ factor ∈ factors, (factor.leadingCoeff : ZMod prime) ≠ 0)
    (hexact : ∀ factor ∈ factors, ∃ (indices : Array Nat) (classId : Nat),
      LegalCombination lifted.size indices.size indices ∧ classId < groups.size ∧
      factor.map (Int.castRingHom (ZMod modulus)) =
        Polynomial.C (factor.leadingCoeff : ZMod modulus) *
          (indices.toList.map fun i => StrictHensel.toPolyMod modulus lifted[i]!).prod ∧
      factor.map (Int.castRingHom (ZMod modulus)) =
        Polynomial.C (factor.leadingCoeff : ZMod modulus) *
          SelectedProductMod modulus groups[classId]! lifted 0 ∧
      ∀ index, index < lifted.size →
        (CandidateMember groups classId index.toUInt32.toInt32 ↔ index ∈ indices.toList)) :
    ∀ classId, classId < groups.size → ∃ factor ∈ factors,
      factor.map (Int.castRingHom (ZMod modulus)) =
        Polynomial.C (factor.leadingCoeff : ZMod modulus) *
          SelectedProductMod modulus groups[classId]! lifted 0 := by
  intro classId hc
  obtain ⟨sourceIndex, hi, hm⟩ := extractCandidates_every_group_has_source
    transform rows lifted.size groups hrun classId hc
  let i : Fin lifted.size := ⟨sourceIndex, hi⟩
  let atoms := fun j : Fin lifted.size => StrictHensel.toPolyMod prime lifted[j.val]
  have hd : atoms i ∣ (factors.map (Polynomial.map (Int.castRingHom (ZMod prime)))).prod := by
    rw [← Polynomial.map_list_prod, hproduct]
    exact hatom i
  obtain ⟨reduced, hmem, hdiv⟩ := prime_dvd_list_product_witness (atoms i) (hprime i) _ hd
  obtain ⟨factor, hf, he⟩ := List.mem_map.mp hmem
  rw [← he] at hdiv
  obtain ⟨indices, other, hlegal, ho, hs, hg, hmembers⟩ := hexact factor hf
  have ha := divisor_support_project_associated prime modulus hdivides lifted indices hlegal
    factor (hleading factor hf) hs
  have hin : i ∈ legalPhysicalSupport lifted.size indices hlegal :=
    (prime_atom_dvd_support_iff atoms hprime hseparate _ i).mp (hdiv.trans ha.dvd)
  have hin' : sourceIndex ∈ indices.toList :=
    (legalPhysicalSupport_mem lifted.size indices hlegal i).mp hin
  have hother := (hmembers sourceIndex hi).mpr hin'
  have heq := extractCandidates_class_unique transform rows lifted.size groups hfits hrun
    sourceIndex other classId hi hother hm
  exact ⟨factor, hf, by simpa only [heq] using hg⟩

end Refinement.StrictRecombine
