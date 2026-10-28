import CLPoly.Refinement.VanHoeijCandidateValidation

set_option autoImplicit false

namespace Refinement.StrictRecombine

/-- Finite physical supports are recognised by divisibility by their distinct
prime atoms. This does not choose an integer irreducible factorization. -/
theorem prime_atom_dvd_support_iff
    {R I : Type*} [CommMonoidWithZero R] [IsCancelMulZero R]
    (atoms : I → R) (hprime : ∀ i, Prime (atoms i))
    (hseparate : ∀ i j, Associated (atoms i) (atoms j) → i = j)
    (indices : List I) (i : I) :
    atoms i ∣ (indices.map atoms).prod ↔ i ∈ indices := by
  induction indices with
  | nil => simp only [List.map_nil, List.prod_nil, List.not_mem_nil, iff_false]
           exact fun h => (hprime i).not_unit (isUnit_of_dvd_one h)
  | cons j tail ih =>
    simp only [List.map_cons, List.prod_cons, List.mem_cons]
    constructor
    · intro hd
      rcases (hprime i).dvd_or_dvd hd with hhead | htail
      · exact Or.inl (hseparate i j ((hprime i).associated_of_dvd (hprime j) hhead))
      · exact Or.inr (ih.mp htail)
    · rintro (rfl | htail)
      · exact dvd_mul_right _ _
      · exact dvd_mul_of_dvd_right (ih.mpr htail) _

theorem associated_support_subset_of_dvd
    {R I : Type*} [CommMonoidWithZero R] [IsCancelMulZero R]
    (atoms : I → R) (hprime : ∀ i, Prime (atoms i))
    (hseparate : ∀ i j, Associated (atoms i) (atoms j) → i = j)
    (small large : List I) (left right : R)
    (hleft : Associated left (small.map atoms).prod)
    (hright : Associated right (large.map atoms).prod)
    (hdivides : left ∣ right) : small ⊆ large := by
  intro i hi
  apply (prime_atom_dvd_support_iff atoms hprime hseparate large i).mp
  have hs := (prime_atom_dvd_support_iff atoms hprime hseparate small i).mpr hi
  exact hs.trans (hleft.symm.dvd.trans (hdivides.trans hright.dvd))

theorem uniform_support_divisor_associated
    {R I : Type*} [CommMonoidWithZero R] [IsCancelMulZero R]
    (atoms : I → R) (hprime : ∀ i, Prime (atoms i))
    (hseparate : ∀ i j, Associated (atoms i) (atoms j) → i = j)
    (small block : List I) (left right : R)
    (hsmall : small.Nodup) (hblock : block.Nodup) (hne : small ≠ [])
    (hleft : Associated left (small.map atoms).prod)
    (hright : Associated right (block.map atoms).prod)
    (hdivides : left ∣ right)
    (huniform : ∀ i ∈ block, ∀ j ∈ block, (i ∈ small ↔ j ∈ small)) :
    Associated left right := by
  have hsubset := associated_support_subset_of_dvd atoms hprime hseparate
    small block left right hleft hright hdivides
  have hexists : ∃ first, first ∈ small := List.exists_mem_of_ne_nil small hne
  obtain ⟨first, hfirst⟩ := hexists
  have hperm : small.Perm block := (List.perm_ext_iff_of_nodup hsmall hblock).mpr (by
    intro i
    exact ⟨fun hi => hsubset hi,
      fun hi => (huniform first (hsubset hfirst) i hi).mp hfirst⟩)
  exact hleft.trans ((Associated.of_eq (hperm.map atoms).prod_eq).trans hright.symm)

end Refinement.StrictRecombine
