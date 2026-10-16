import CLPoly.Refinement.VanHoeijDivisorGrouping
import CLPoly.Refinement.VanHoeijValidationNonunit

set_option autoImplicit false

namespace Refinement.StrictRecombine
open CLPoly.Math

theorem hensel_lattice_validation_irreducible
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
      selection.factors selection.prime aTarget output)
    (hindexFits : output.1.size ≤ 2 ^ 31)
    (original truncated : Array SparsePolyZZ)
    (modulus : Nat) (hmod : output.2 = (modulus : Int))
    (divisor height scale : ZZ) (target columns : Nat)
    (initial matrix : Generated.StrictRecombine.LLLMatrix)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hnonempty : 0 < f.size)
    (hfits : (SparsePolyZZ.toPoly f).natDegree < 2 ^ 31)
    (hliftedDegree : ∀ index (hi : index < output.1.size),
      output.1[index][0]!.1.deg < 2 ^ 63)
    (hcldRun : Generated.StrictRecombine.cldPolys ⟨termination⟩ f output.1 modulus = .ok original)
    (htruncate : Generated.StrictVanHoeij.truncateCld original divisor = .ok truncated)
    (hscale : scale ≠ 0)
    (hinit : Generated.StrictRecombine.makeInitialMatrix truncated.size scale = .ok initial)
    (hbuild : Generated.StrictVanHoeij.buildCldMatrix initial truncated 0 target
      ((modulus : Int) / divisor) = .ok (matrix, columns))
    (hm : 1 < modulus) (hd : 0 < divisor) (hh : height < (modulus : Int))
    (hheightRun : Generated.StrictVanHoeij.cldBound f = .ok height)
    (bound : ZZ)
    (hbound : bound = (output.1.size : ZZ) * scale * scale + (columns : ZZ) *
      ((height + divisor - 1) / divisor + 2 * (output.1.size : Int)) *
      ((height + divisor - 1) / divisor + 2 * (output.1.size : Int)))
    (hinput : ConcreteLLLInputValid matrix)
    (reduced : Generated.StrictRecombine.LLLMatrix ×
      Generated.StrictRecombine.LLLMatrix × Array Nat × Array QQ)
    (hlll : Generated.StrictVanHoeij.lllReduceWithNorms concreteLLLExecution
      matrix hinput bound = .ok reduced)
    (candidates : Array (Array Int32))
    (hextract : Generated.StrictRecombine.extractCandidates
      (Generated.StrictVanHoeij.shortSpanFromNorms reduced.2.2.2 bound)
      reduced.2.1 output.1.size = .ok candidates)
    (hprimitive : (SparsePolyZZ.toPoly f).IsPrimitive)
    (before after : Generated.StrictVanHoeij.VanHoeijResult)
    (hvalidate : Generated.StrictVanHoeij.validatePartition f output.1 modulus
      candidates before = .ok (true, after)) :
    factorArrayProduct after.factors = SparsePolyZZ.toPoly f ∧
      ∀ factor ∈ after.factors.toList, Irreducible (SparsePolyZZ.toPoly factor) := by
  have hdegree63 : (SparsePolyZZ.toPoly f).natDegree < 2 ^ 63 := by omega
  have hp := validatePartition_product f output.1 modulus candidates before after
    hcanonical hdegree63 hvalidate
  refine ⟨hp, ?_⟩
  have hrows : (Generated.StrictVanHoeij.shortSpanFromNorms reduced.2.2.2 bound).isEmpty = false := by
    cases he : (Generated.StrictVanHoeij.shortSpanFromNorms reduced.2.2.2 bound).isEmpty with
    | false => rfl
    | true =>
      have hgroups : candidates = #[] := by
        simpa [Generated.StrictRecombine.extractCandidates, he] using hextract.symm
      simp [hgroups, Generated.StrictVanHoeij.validatePartition] at hvalidate
  have hv := extractCandidates_indices_valid
    (Generated.StrictVanHoeij.shortSpanFromNorms reduced.2.2.2 bound) reduced.2.1
    output.1.size candidates hindexFits hrows hextract
  have hirr := StrictFactorZZ.selectionHenselFactors_mod_irreducible hcount hfactors
    hleadingSemantic hselection hentry
  obtain ⟨exponent, he, hpower⟩ := hentry.outputModulus_eq_prime_pow
  have hmodulus : modulus = selection.prime.toNat ^ exponent := by
    exact_mod_cast hmod.symm.trans hpower
  have hdivides : selection.prime.toNat ∣ modulus := by
    rw [hmodulus]
    exact dvd_pow_self _ (Nat.ne_of_gt he)
  have hnonunit := validatePartition_nonunit f output.1 modulus selection.prime.toNat
    candidates before after (by omega) hdivides hv hirr hcanonical
    hselection.goodPrime.lc_nonzero hvalidate
  have hprim := validatePartition_factor_primitive f output.1 modulus candidates before after
    hcanonical hdegree63 hprimitive hvalidate
  have hsq := validatePartition_factor_mod_squarefree f output.1 modulus selection.prime.toNat
    candidates before after hcanonical hdegree63 hselection.goodPrime.sqfree hvalidate
  have hsupport := hensel_validated_factor_physical_support hcount hfactors hleadingSemantic
    hselection hentry modulus hmod candidates
    (Generated.StrictVanHoeij.shortSpanFromNorms reduced.2.2.2 bound) reduced.2.1 before after
    hindexFits hcanonical hdegree63 hextract hvalidate
  obtain ⟨hprimes, hseparate⟩ := hensel_physical_atoms_prime_separated hcount hfactors
    hleadingSemantic hselection hentry
  intro factor hf
  have hfdvd : SparsePolyZZ.toPoly factor ∣ SparsePolyZZ.toPoly f := by
    rw [← hp]
    exact List.dvd_prod (List.mem_map.mpr ⟨factor, hf, rfl⟩)
  have hleading := divisor_leading_mod_nonzero selection.prime.toNat
    (SparsePolyZZ.toPoly f) (SparsePolyZZ.toPoly factor) hselection.goodPrime.lc_nonzero hfdvd
  obtain ⟨group, hg, hvalid, _, hs⟩ := hsupport factor hf
  apply irreducible_of_uniform_divisor_support selection.prime.toNat
    (fun i : Fin output.1.size => StrictHensel.toPolyMod selection.prime.toNat output.1[i.val])
    hprimes hseparate (SparsePolyZZ.toPoly factor)
    (candidatePhysicalSupport output.1.size group hvalid)
    (hprim factor hf) (hnonunit factor hf) hleading (hsq factor hf) hs
  intro g hgdvd hgnonunit
  have hgf : g ∣ SparsePolyZZ.toPoly f := hgdvd.trans hfdvd
  have hgprimitive := Polynomial.isPrimitive_of_dvd hprimitive hgf
  have hgdegree := primitive_nonunit_natDegree_pos g hgprimitive hgnonunit
  obtain ⟨complement, hgproduct⟩ := hgf
  obtain ⟨small, hsmall, huniform⟩ := hensel_execution_physical_divisor_support
    hcount hfactors hleadingSemantic hselection hentry hindexFits original truncated modulus hmod
    divisor height scale target columns initial matrix hcanonical hnonempty hfits
    hliftedDegree hcldRun htruncate hscale hinit hbuild g complement hgproduct
    hm hd hh hgdegree hheightRun bound hbound hinput reduced hlll candidates hextract
  obtain ⟨classId, hc, heq⟩ := List.mem_iff_getElem.mp hg
  have hc' : classId < candidates.size := by simpa using hc
  have heq' : candidates[classId] = group := by simpa only [Array.getElem_toList] using heq
  subst group
  exact ⟨small, hsmall, huniform classId hc' hvalid⟩

end Refinement.StrictRecombine
