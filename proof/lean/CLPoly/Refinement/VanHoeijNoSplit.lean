import CLPoly.Refinement.VanHoeijPhysicalGrouping

set_option autoImplicit false

namespace Refinement.StrictRecombine

theorem uniform_support_factorization_unit
    {R I : Type*} [CommMonoidWithZero R] [IsCancelMulZero R] [Nontrivial R]
    (atoms : I → R) (hprime : ∀ i, Prime (atoms i))
    (hseparate : ∀ i j, Associated (atoms i) (atoms j) → i = j)
    (small block : List I) (left right complement : R)
    (hsquarefree : Squarefree right)
    (hleft : Associated left (small.map atoms).prod)
    (hright : Associated right (block.map atoms).prod)
    (hfactor : right = left * complement) (hleftNonunit : ¬ IsUnit left)
    (huniform : ∀ i ∈ block, ∀ j ∈ block, (i ∈ small ↔ j ∈ small)) :
    IsUnit complement := by
  have hd : left ∣ right := ⟨complement, hfactor⟩
  have hsmall := associated_support_nodup atoms (fun i => (hprime i).not_unit)
    small left (hsquarefree.squarefree_of_dvd hd) hleft
  have hblock := associated_support_nodup atoms (fun i => (hprime i).not_unit)
    block right hsquarefree hright
  have hne : small ≠ [] := by
    intro he
    apply hleftNonunit
    exact associated_one_iff_isUnit.mp (by simpa only [he, List.map_nil, List.prod_nil] using hleft)
  have ha := uniform_support_divisor_associated atoms hprime hseparate small block
    left right hsmall hblock hne hleft hright hd huniform
  have hleftNe : left ≠ 0 := by
    intro hz
    exact hsquarefree.ne_zero (by simpa only [hz, zero_mul] using hfactor)
  have hc : Associated (left * 1) (left * complement) := by
    simpa only [mul_one, ← hfactor] using ha
  exact associated_one_iff_isUnit.mp
    (Associated.of_mul_left hc (Associated.refl left) hleftNe).symm

theorem irreducible_of_uniform_divisor_support
    {I : Type*} (prime : Nat) [Fact (Nat.Prime prime)]
    (atoms : I → Polynomial (ZMod prime)) (hprime : ∀ i, Prime (atoms i))
    (hseparate : ∀ i j, Associated (atoms i) (atoms j) → i = j)
    (factor : Polynomial Int) (block : List I)
    (hprimitive : factor.IsPrimitive) (hnonunit : ¬ IsUnit factor)
    (hleading : (factor.leadingCoeff : ZMod prime) ≠ 0)
    (hsquarefree : Squarefree (factor.map (Int.castRingHom (ZMod prime))))
    (hblock : Associated (factor.map (Int.castRingHom (ZMod prime))) (block.map atoms).prod)
    (hdivisor : ∀ divisor : Polynomial Int, divisor ∣ factor → ¬ IsUnit divisor →
      ∃ small : List I,
        Associated (divisor.map (Int.castRingHom (ZMod prime))) (small.map atoms).prod ∧
        ∀ i ∈ block, ∀ j ∈ block, (i ∈ small ↔ j ∈ small)) :
    Irreducible factor := by
  refine ⟨hnonunit, ?_⟩
  intro left right hfactor
  by_cases hl : IsUnit left
  · exact Or.inl hl
  by_cases hr : IsUnit right
  · exact Or.inr hr
  exfalso
  obtain ⟨_, _, _, _, hlmap, hrmap⟩ :=
    StrictFactorZZ.primitive_factorization_maps_nonunits prime factor left right
      hprimitive hleading hfactor hl hr
  obtain ⟨small, hs, hu⟩ := hdivisor left ⟨right, hfactor⟩ hl
  apply hrmap
  apply uniform_support_factorization_unit atoms hprime hseparate small block
    (left.map (Int.castRingHom (ZMod prime)))
    (factor.map (Int.castRingHom (ZMod prime)))
    (right.map (Int.castRingHom (ZMod prime))) hsquarefree hs hblock
    (by rw [hfactor, Polynomial.map_mul]) hlmap hu

end Refinement.StrictRecombine
