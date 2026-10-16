import CLPoly.Refinement.VanHoeijDivisorSupport

set_option autoImplicit false

namespace Refinement.StrictRecombine
open CLPoly.Math

theorem hensel_execution_physical_divisor_support
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
    (factor complement : Polynomial Int)
    (hfactor : SparsePolyZZ.toPoly f = factor * complement)
    (hm : 1 < modulus) (hd : 0 < divisor) (hh : height < (modulus : Int))
    (hdegree : 0 < factor.natDegree)
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
      reduced.2.1 output.1.size = .ok candidates) :
    ∃ small : List (Fin output.1.size),
      Associated (factor.map (Int.castRingHom (ZMod selection.prime.toNat)))
        ((small.map fun i => StrictHensel.toPolyMod selection.prime.toNat output.1[i.val]).prod) ∧
      ∀ classId (hc : classId < candidates.size)
        (hv : CandidateIndicesValid candidates[classId] (Array.replicate output.1.size false)),
        ∀ left ∈ candidatePhysicalSupport output.1.size candidates[classId] hv,
          ∀ right ∈ candidatePhysicalSupport output.1.size candidates[classId] hv,
            (left ∈ small ↔ right ∈ small) := by
  obtain ⟨indices, hlegal, hs, hu⟩ := hensel_execution_group_support
    hcount hfactors hleadingSemantic hselection hentry original truncated modulus hmod
    divisor height scale target columns initial matrix hcanonical hnonempty hfits
    hliftedDegree hcldRun htruncate hscale hinit hbuild factor complement hfactor
    hm hd hh hdegree hheightRun bound hbound hinput reduced hlll candidates hextract
  obtain ⟨exponent, he, hpower⟩ := hentry.outputModulus_eq_prime_pow
  have hmodulus : modulus = selection.prime.toNat ^ exponent := by
    exact_mod_cast hmod.symm.trans hpower
  have hdivides : selection.prime.toNat ∣ modulus := by
    rw [hmodulus]
    exact dvd_pow_self _ (Nat.ne_of_gt he)
  have hmul : (factor.leadingCoeff : ZMod selection.prime.toNat) *
      (complement.leadingCoeff : ZMod selection.prime.toNat) ≠ 0 := by
    simpa only [hfactor, Polynomial.leadingCoeff_mul, Int.cast_mul] using
      hselection.goodPrime.lc_nonzero
  refine ⟨legalPhysicalSupport output.1.size indices hlegal,
    divisor_support_project_associated selection.prime.toNat modulus hdivides
      output.1 indices hlegal factor (left_ne_zero_of_mul hmul) hs, ?_⟩
  intro classId hc hv left hl right hr
  rw [legalPhysicalSupport_mem, legalPhysicalSupport_mem]
  exact candidatePhysicalSupport_uniform output.1.size candidates classId hc hv
    hindexFits indices (hu classId) left hl right hr

end Refinement.StrictRecombine
