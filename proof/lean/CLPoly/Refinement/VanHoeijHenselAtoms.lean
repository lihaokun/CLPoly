import CLPoly.Refinement.VanHoeijSupportSeparation

set_option autoImplicit false

namespace Refinement.StrictRecombine
open CLPoly.Math

theorem hensel_physical_atoms_prime_separated
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {aTarget : Int32} {output : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hfactors : ∀ factor ∈ selection.factors.toList,
      SparsePolyZp.Canonical selection.prime.toNat factor)
    (hleadingSemantic : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) =
        (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f
      selection.factors selection.prime aTarget output) :
    (∀ i : Fin output.1.size,
      Prime (StrictHensel.toPolyMod selection.prime.toNat output.1[i.val])) ∧
    (∀ i j : Fin output.1.size,
      Associated (StrictHensel.toPolyMod selection.prime.toNat output.1[i.val])
        (StrictHensel.toPolyMod selection.prime.toNat output.1[j.val]) → i = j) := by
  have hirr := StrictFactorZZ.selectionHenselFactors_mod_irreducible hcount hfactors
    hleadingSemantic hselection hentry
  have hcoprime := StrictFactorZZ.henselFactors_mod_pairwise_coprime hcount
    (fun adjusted hadjust => StrictFactorZZ.selectionAdjustedFactors_pairwise_isCoprime
      hfactors hleadingSemantic hselection hadjust) hentry
  constructor
  · intro i
    exact (hirr i.val i.isLt).prime
  · intro i j ha
    apply Fin.ext
    by_contra hne
    rcases lt_or_gt_of_ne hne with hij | hji
    · exact (hirr i.val i.isLt).not_isUnit
        ((hcoprime i.val j.val i.isLt j.isLt hij).isUnit_of_dvd ha.dvd)
    · exact (hirr j.val j.isLt).not_isUnit
        ((hcoprime j.val i.val j.isLt i.isLt hji).isUnit_of_dvd ha.symm.dvd)

end Refinement.StrictRecombine
