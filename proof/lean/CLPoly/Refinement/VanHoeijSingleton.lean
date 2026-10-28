import CLPoly.Refinement.VanHoeijIrreducibility

set_option autoImplicit false

namespace Refinement.StrictRecombine

theorem primitive_irreducible_of_mod_irreducible
    (prime : Nat) [Fact (Nat.Prime prime)] (factor : Polynomial Int)
    (hp : factor.IsPrimitive) (hl : (factor.leadingCoeff : ZMod prime) ≠ 0)
    (hirr : Irreducible (factor.map (Int.castRingHom (ZMod prime)))) :
    Irreducible factor := by
  refine ⟨?_, ?_⟩
  · intro hu
    exact hirr.not_isUnit (hu.map (Polynomial.mapRingHom (Int.castRingHom (ZMod prime))))
  · intro left right heq
    by_cases hlu : IsUnit left
    · exact Or.inl hlu
    by_cases hru : IsUnit right
    · exact Or.inr hru
    obtain ⟨_, _, _, _, hlm, hrm⟩ := StrictFactorZZ.primitive_factorization_maps_nonunits
      prime factor left right hp hl heq hlu hru
    have hm : factor.map (Int.castRingHom (ZMod prime)) =
        left.map (Int.castRingHom (ZMod prime)) * right.map (Int.castRingHom (ZMod prime)) := by
      rw [heq, Polynomial.map_mul]
    exact False.elim ((hirr.isUnit_or_isUnit hm).elim hlm hrm)

theorem singleton_selected_product_irreducible
    (prime : Nat) [Fact (Nat.Prime prime)] (lifted : Array SparsePolyZZ)
    (group : Array Int32) (hsingle : group.size = 1)
    (hvalid : CandidateIndicesValid group (Array.replicate lifted.size false))
    (hirr : ∀ i (hi : i < lifted.size), Irreducible (StrictHensel.toPolyMod prime lifted[i])) :
    Irreducible (SelectedProductMod prime group lifted 0) := by
  have hn : 0 < group.size := by omega
  have hlist : group.toList = [group[0]] := by
    simpa only [Array.getElem_toList] using
      List.eq_getElem_of_length_eq_one group.toList (by simpa using hsingle)
  have hi : group[0].toInt64.toNat < lifted.size := by
    simpa only [Array.size_replicate] using (hvalid 0 hn).2
  simpa only [SelectedProductMod, List.drop_zero, hlist, List.map_cons, List.map_nil,
    List.prod_cons, List.prod_nil, mul_one, getElem!_pos lifted group[0].toInt64.toNat hi] using
    hirr group[0].toInt64.toNat hi

theorem validatePartition_singletons_irreducible
    (f : SparsePolyZZ) (lifted : Array SparsePolyZZ) (modulus prime : Nat)
    [Fact (Nat.Prime prime)] (groups : Array (Array Int32))
    (before after : Generated.StrictVanHoeij.VanHoeijResult)
    (hm : 0 < modulus) (hdivides : prime ∣ modulus)
    (hvalid : ∀ i (hi : i < groups.size),
      CandidateIndicesValid groups[i] (Array.replicate lifted.size false))
    (hsingle : ∀ group ∈ groups.toList, group.size = 1)
    (hirr : ∀ i (hi : i < lifted.size), Irreducible (StrictHensel.toPolyMod prime lifted[i]))
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hdegree : (SparsePolyZZ.toPoly f).natDegree < 2 ^ 63)
    (hprimitive : (SparsePolyZZ.toPoly f).IsPrimitive)
    (hleading : ((SparsePolyZZ.toPoly f).leadingCoeff : ZMod prime) ≠ 0)
    (hrun : Generated.StrictVanHoeij.validatePartition f lifted modulus groups before =
      .ok (true, after)) :
    ∀ factor ∈ after.factors.toList, Irreducible (SparsePolyZZ.toPoly factor) := by
  have hp := validatePartition_product f lifted modulus groups before after hcanonical hdegree hrun
  have hprim := validatePartition_factor_primitive f lifted modulus groups before after
    hcanonical hdegree hprimitive hrun
  have hs := validatePartition_support f lifted modulus prime groups before after hm hdivides
    hvalid hirr hcanonical hleading hrun
  intro factor hf
  have hd : SparsePolyZZ.toPoly factor ∣ SparsePolyZZ.toPoly f := by
    rw [← hp]
    exact List.dvd_prod (List.mem_map.mpr ⟨factor, hf, rfl⟩)
  obtain ⟨group, hg, ha⟩ := hs factor hf
  obtain ⟨i, hi, heq⟩ := List.mem_iff_getElem.mp hg
  have hi' : i < groups.size := by simpa using hi
  have heq' : groups[i] = group := by simpa only [Array.getElem_toList] using heq
  have hv : CandidateIndicesValid group (Array.replicate lifted.size false) := by
    rw [← heq']
    exact hvalid i hi'
  apply primitive_irreducible_of_mod_irreducible prime _ (hprim factor hf)
    (divisor_leading_mod_nonzero prime _ _ hleading hd)
  exact ha.irreducible_iff.mpr
    (singleton_selected_product_irreducible prime lifted group (hsingle group hg) hv hirr)

end Refinement.StrictRecombine
