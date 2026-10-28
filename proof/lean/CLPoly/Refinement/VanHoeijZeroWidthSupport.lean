import CLPoly.Refinement.VanHoeijZeroWidth
import CLPoly.Refinement.VanHoeijRetainedNonempty

set_option autoImplicit false

namespace Refinement.StrictRecombine

/-- The same Hensel execution supplies a singleton full-modulus support for
every genuine irreducible divisor when all truncated CLDs vanish. -/
theorem hensel_zero_width_divisor_singleton_support
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {aTarget : Int32} {lifted : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hfactors : StrictSelectPrime.SelectionPhysical selection)
    (hleading : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) = (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f selection.factors selection.prime aTarget lifted)
    (modulus : Nat) (hmodulus : lifted.2 = (modulus : Int))
    (original truncated : Array SparsePolyZZ)
    (data : CldSeparationData f lifted.1 original modulus)
    (divisor : ZZ)
    (htruncate : Generated.StrictVanHoeij.truncateCld original divisor = .ok truncated)
    (hwidth : Generated.StrictRecombine.cldSpiralWidth truncated = 0)
    (factor complement : Polynomial Int) (hfactor : SparsePolyZZ.toPoly f = factor * complement)
    (hprimitive : factor.IsPrimitive)
    (hirreducible : Irreducible (factor.map (Int.castRingHom QQ)))
    (hfactorDegree : 0 < factor.natDegree)
    (hresultant : separationResultantThreshold (SparsePolyZZ.toPoly f) lifted.1.size divisor < (modulus : Int))
    (hcoordinates : 2 * separationCoordinateBound lifted.1.size (SparsePolyZZ.toPoly f).natDegree < (modulus : Int)) :
    ∃ index : Fin lifted.1.size,
      factor.map (Int.castRingHom (ZMod modulus)) =
        Polynomial.C (factor.leadingCoeff : ZMod modulus) *
          StrictHensel.toPolyMod modulus lifted.1[index.val]! := by
  obtain ⟨indices, exponent, hlegal, he, hp, hs⟩ :=
    StrictFactorZZ.hensel_divisor_physical_support hcount hfactors hleading hselection hentry
      factor complement hfactor
  have hpower : modulus = selection.prime.toNat ^ exponent := by
    exact_mod_cast hmodulus.symm.trans hp
  have hsupport : factor.map (Int.castRingHom (ZMod modulus)) =
      Polynomial.C (factor.leadingCoeff : ZMod modulus) *
        (indices.toList.map fun i => StrictHensel.toPolyMod modulus lifted.1[i]!).prod := by
    rw [hpower]
    exact hs
  have hprimeDvd : selection.prime.toNat ∣ modulus := by
    rw [hpower]
    exact dvd_pow_self _ (Nat.ne_of_gt he)
  have hfactorLeading : (factor.leadingCoeff : ZMod selection.prime.toNat) ≠ 0 :=
    divisor_leading_mod_nonzero selection.prime.toNat (SparsePolyZZ.toPoly f) factor
      hselection.goodPrime.lc_nonzero ⟨complement, hfactor⟩
  have hnonempty := divisor_modular_support_nonempty selection.prime.toNat modulus hprimeDvd
    lifted.1 indices hlegal factor hfactorDegree hfactorLeading hsupport
  have hsizePositive : 0 < indices.size := by
    by_contra h
    have hempty : indices = #[] := Array.size_eq_zero_iff.mp (by omega)
    exact hnonempty (by simp [hempty])
  have hsize : indices.size = 1 := by
    by_contra hne
    have hsecond : 1 < indices.size := by omega
    have hl := hlegal.2.2 0 hsizePositive
    have hr := hlegal.2.2 1 hsecond
    have hl' : indices[0]! < original.size := by rw [data.size_eq]; exact hl
    have hr' : indices[1]! < original.size := by rw [data.size_eq]; exact hr
    have hlmem : indices[0]! ∈ indices.toList := by
      rw [getElem!_pos indices 0 hsizePositive]
      exact Array.getElem_mem_toList hsizePositive
    have hrmem : indices[1]! ∈ indices.toList := by
      rw [getElem!_pos indices 1 hsecond]
      exact Array.getElem_mem_toList hsecond
    have heq := data.zero_width_local_indices_equal truncated divisor htruncate hwidth
      factor ⟨complement, hfactor⟩ hprimitive hirreducible (Nat.ne_of_gt hfactorDegree)
      ⟨indices[0]!, hl'⟩ ⟨indices[1]!, hr'⟩
      (local_factor_dvd_of_support lifted.1 modulus factor indices.toList _ hsupport _ hlmem)
      (local_factor_dvd_of_support lifted.1 modulus factor indices.toList _ hsupport _ hrmem)
      hresultant hcoordinates
    have hv : indices[0]! = indices[1]! := congrArg Fin.val heq
    have hgap := hlegal.2.1 0 1 hsizePositive hsecond (by decide)
    omega
  obtain ⟨index, hsingle⟩ := List.length_eq_one_iff.mp
    (show indices.toList.length = 1 by simpa using hsize)
  have hindex : index < lifted.1.size := by
    have heq : indices[0]! = index := by
      have hlist := congrArg (fun l : List Nat => l[0]!) hsingle
      simpa only [getElem!_pos indices.toList 0 (by simpa using hsizePositive),
        getElem!_pos indices 0 hsizePositive, Array.getElem_toList,
        List.getElem!_cons_zero] using hlist
    rw [← heq]
    exact hlegal.2.2 0 hsizePositive
  refine ⟨⟨index, hindex⟩, ?_⟩
  simpa only [hsingle, List.map_cons, List.map_nil, List.prod_cons, List.prod_nil,
    mul_one] using hsupport

end Refinement.StrictRecombine
