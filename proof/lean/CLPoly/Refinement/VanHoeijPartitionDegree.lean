import CLPoly.Refinement.VanHoeijPartitionProductExecution
import CLPoly.Refinement.VanHoeijDivisorSupport

set_option autoImplicit false

namespace Refinement.StrictRecombine

theorem monic_list_product_natDegree {R : Type*} [CommRing R]
    (factors : List (Polynomial R)) (hmonic : ∀ factor ∈ factors, factor.Monic) :
    factors.prod.natDegree = (factors.map Polynomial.natDegree).sum := by
  have h := Polynomial.natDegree_multiset_prod_of_monic (↑factors : Multiset (Polynomial R))
    (by simpa using hmonic)
  simpa using h

theorem selectedProductMod_natDegree
    (prime : Nat) [Fact (Nat.Prime prime)] (group : Array Int32) (lifted : Array SparsePolyZZ)
    (hvalid : CandidateIndicesValid group (Array.replicate lifted.size false))
    (hmonic : ∀ index, index < lifted.size → (SparsePolyZZ.toPoly lifted[index]!).Monic) :
    (SelectedProductMod prime group lifted 0).natDegree =
      (group.toList.map fun source => (SparsePolyZZ.toPoly lifted[source.toInt64.toNat]!).natDegree).sum := by
  have hm : ∀ source ∈ group.toList, (SparsePolyZZ.toPoly lifted[source.toInt64.toNat]!).Monic := by
    intro source hs
    obtain ⟨position, hp, hv⟩ := List.mem_iff_getElem.mp hs
    have hp' : position < group.size := by simpa using hp
    have hb : source.toInt64.toNat < lifted.size := by
      have h := (hvalid position hp').2
      simpa only [Array.size_replicate, ← Array.getElem_toList hp', hv] using h
    exact hmonic _ hb
  unfold SelectedProductMod
  rw [List.drop_zero, monic_list_product_natDegree _ (by
    intro factor hf
    obtain ⟨source, hs, rfl⟩ := List.mem_map.mp hf
    exact (hm source hs).map (Int.castRingHom (ZMod prime)))]
  rw [List.map_map]
  congr 1
  apply List.map_congr_left
  intro source hs
  exact (hm source hs).natDegree_map (Int.castRingHom (ZMod prime))

theorem selected_support_project
    (prime modulus : Nat) (hd : prime ∣ modulus)
    (group : Array Int32) (lifted : Array SparsePolyZZ) (factor : Polynomial Int)
    (hsupport : factor.map (Int.castRingHom (ZMod modulus)) =
      Polynomial.C (factor.leadingCoeff : ZMod modulus) * SelectedProductMod modulus group lifted 0) :
    factor.map (Int.castRingHom (ZMod prime)) =
      Polynomial.C (factor.leadingCoeff : ZMod prime) * SelectedProductMod prime group lifted 0 := by
  let product := (group.toList.map fun source => SparsePolyZZ.toPoly lifted[source.toInt64.toNat]!).prod
  have he : factor.map (Int.castRingHom (ZMod modulus)) =
      (Polynomial.C factor.leadingCoeff * product).map (Int.castRingHom (ZMod modulus)) := by
    simpa only [product, SelectedProductMod, List.drop_zero, Polynomial.map_mul, Polynomial.map_C,
      Polynomial.map_list_prod, List.map_map, Function.comp_def, StrictHensel.toPolyMod] using hsupport
  have hp := polynomialMap_eq_of_modulus_dvd prime modulus hd factor
    (Polynomial.C factor.leadingCoeff * product) he
  simpa only [product, SelectedProductMod, List.drop_zero, Polynomial.map_mul, Polynomial.map_C,
    Polynomial.map_list_prod, List.map_map, Function.comp_def, StrictHensel.toPolyMod] using hp

theorem selected_support_degree_sum
    (prime : Nat) [Fact (Nat.Prime prime)] (group : Array Int32) (lifted : Array SparsePolyZZ)
    (factor : Polynomial Int) (hleading : (factor.leadingCoeff : ZMod prime) ≠ 0)
    (hvalid : CandidateIndicesValid group (Array.replicate lifted.size false))
    (hmonic : ∀ index, index < lifted.size → (SparsePolyZZ.toPoly lifted[index]!).Monic)
    (hsupport : factor.map (Int.castRingHom (ZMod prime)) =
      Polynomial.C (factor.leadingCoeff : ZMod prime) * SelectedProductMod prime group lifted 0) :
    (group.toList.map fun source => (SparsePolyZZ.toPoly lifted[source.toInt64.toNat]!).natDegree).sum =
      factor.natDegree := by
  have hd := congrArg Polynomial.natDegree hsupport
  rw [Polynomial.natDegree_map_of_leadingCoeff_ne_zero _ hleading,
    Polynomial.natDegree_C_mul hleading, selectedProductMod_natDegree prime group lifted hvalid hmonic] at hd
  exact hd.symm

/-- The product loop computes the true factor degree, including the Int64
positivity check, once its physical group is known to be an exact support. -/
theorem partitionProductLoop_succeeds_with_factor_degree
    (prime modulus : Nat) [Fact (Nat.Prime prime)] (hm : 0 < modulus) (hd : prime ∣ modulus)
    (group : Array Int32) (lifted : Array SparsePolyZZ) (initial : SparsePolyZZ)
    (factor : Polynomial Int) (hleading : (factor.leadingCoeff : ZMod prime) ≠ 0)
    (hvalid : CandidateIndicesValid group (Array.replicate lifted.size false))
    (hcanonical : ∀ index, index < lifted.size → StrictPolynomialMod.SparsePolyZZCanonical lifted[index]!)
    (hnonempty : ∀ index, index < lifted.size → 0 < lifted[index]!.size)
    (hmonic : ∀ index, index < lifted.size → (SparsePolyZZ.toPoly lifted[index]!).Monic)
    (hsupport : factor.map (Int.castRingHom (ZMod modulus)) =
      Polynomial.C (factor.leadingCoeff : ZMod modulus) * SelectedProductMod modulus group lifted 0)
    (hdegree : 0 < factor.natDegree) (hfits : factor.natDegree < 2 ^ 63) :
    ∃ product,
      Generated.StrictVanHoeij.partitionProductLoop group lifted modulus 0 initial 0 =
        .ok (product, Int64.ofNat factor.natDegree) ∧
      0 < Int64.ofNat factor.natDegree := by
  obtain ⟨product, degree, hrun⟩ := partitionProductLoop_succeeds group lifted modulus 0 initial 0
    (by exact_mod_cast Nat.ne_of_gt hm) hvalid
  have hcounter := partitionProductLoop_degree_ofNat group lifted modulus initial product degree
    hvalid hcanonical hnonempty hrun
  rw [selected_support_degree_sum prime group lifted factor hleading hvalid hmonic
    (selected_support_project prime modulus hd group lifted factor hsupport)] at hcounter
  subst degree
  refine ⟨product, hrun, ?_⟩
  rw [Int64.lt_iff_toInt_lt, Int64.toInt_ofNat_of_lt hfits]
  change (0 : Int) < (factor.natDegree : Int)
  exact_mod_cast hdegree

end Refinement.StrictRecombine
