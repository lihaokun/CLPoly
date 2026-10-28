import CLPoly.Refinement.VanHoeijHenselAtoms

set_option autoImplicit false

namespace Refinement.StrictRecombine

theorem support_nodup_of_squarefree_product
    {R I : Type*} [CommMonoid R] (atoms : I → R)
    (hnonunit : ∀ i, ¬ IsUnit (atoms i))
    (indices : List I) (hsquarefree : Squarefree (indices.map atoms).prod) :
    indices.Nodup := by
  induction indices with
  | nil => exact List.nodup_nil
  | cons i tail ih =>
    change Squarefree (atoms i * (tail.map atoms).prod) at hsquarefree
    refine List.nodup_cons.mpr ⟨?_, ih hsquarefree.of_mul_right⟩
    intro hi
    have hd : atoms i ∣ (tail.map atoms).prod :=
      List.dvd_prod (List.mem_map.mpr ⟨i, hi, rfl⟩)
    exact hnonunit i (hsquarefree (atoms i) (mul_dvd_mul (dvd_refl _) hd))

theorem validatePartition_factor_mod_squarefree
    (f : SparsePolyZZ) (lifted : Array SparsePolyZZ) (modulus : ZZ)
    (prime : Nat) (groups : Array (Array Int32))
    (before after : Generated.StrictVanHoeij.VanHoeijResult)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hdegree : (SparsePolyZZ.toPoly f).natDegree < 2 ^ 63)
    (hsquarefree : Squarefree (StrictHensel.toPolyMod prime f))
    (hrun : Generated.StrictVanHoeij.validatePartition f lifted modulus groups before =
      .ok (true, after)) :
    ∀ factor ∈ after.factors.toList, Squarefree (StrictHensel.toPolyMod prime factor) := by
  have hp := validatePartition_product f lifted modulus groups before after hcanonical hdegree hrun
  intro factor hf
  have hd : SparsePolyZZ.toPoly factor ∣ SparsePolyZZ.toPoly f := by
    rw [← hp]
    exact List.dvd_prod (List.mem_map.mpr ⟨factor, hf, rfl⟩)
  have hm : StrictHensel.toPolyMod prime factor ∣ StrictHensel.toPolyMod prime f := by
    obtain ⟨q, hq⟩ := hd
    refine ⟨Polynomial.map (Int.castRingHom (ZMod prime)) q, ?_⟩
    change Polynomial.map (Int.castRingHom (ZMod prime)) (SparsePolyZZ.toPoly f) = _
    rw [hq, Polynomial.map_mul]
    rfl
  exact hsquarefree.squarefree_of_dvd hm

theorem associated_support_nodup
    {R I : Type*} [CommMonoid R] (atoms : I → R)
    (hnonunit : ∀ i, ¬ IsUnit (atoms i)) (indices : List I) (factor : R)
    (hsquarefree : Squarefree factor)
    (hsupport : Associated factor (indices.map atoms).prod) : indices.Nodup :=
  support_nodup_of_squarefree_product atoms hnonunit indices
    (hsupport.squarefree_iff.mp hsquarefree)

end Refinement.StrictRecombine
