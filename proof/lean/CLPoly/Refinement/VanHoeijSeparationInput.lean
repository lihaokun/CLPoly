import CLPoly.Refinement.VanHoeijHenselLocal
import CLPoly.Refinement.VanHoeijCldDegree
import CLPoly.Refinement.VanHoeijExecution

set_option autoImplicit false

namespace Refinement.StrictRecombine

/-- Erased facts about the same physical Hensel and CLD arrays, constructed
below from their actual executions. It contains no factorization-output
irreducibility or assumed separation property. -/
structure CldSeparationData (f : SparsePolyZZ) (active cld : Array SparsePolyZZ)
    (modulus : Nat) : Prop where
  size_eq : cld.size = active.size
  modulus_gt_one : 1 < modulus
  source_degree_pos : 0 < (SparsePolyZZ.toPoly f).natDegree
  values : ∀ i, i < cld.size → ModularCldValue f active[i]! cld[i]! modulus
  degrees : ∀ i, i < cld.size →
    (SparsePolyZZ.toPoly cld[i]!).natDegree < (SparsePolyZZ.toPoly f).natDegree
  monic : ∀ i, i < cld.size → (StrictHensel.toPolyMod modulus active[i]!).Monic
  local_degree_pos : ∀ i, i < cld.size → 0 < (StrictHensel.toPolyMod modulus active[i]!).degree
  divides : ∀ i, i < cld.size →
    StrictHensel.toPolyMod modulus active[i]! ∣ StrictHensel.toPolyMod modulus f
  pairwise : ∀ i j, i < cld.size → j < cld.size → i ≠ j →
    IsCoprime (StrictHensel.toPolyMod modulus active[i]!) (StrictHensel.toPolyMod modulus active[j]!)
  derivative_coprime : ∀ i, i < cld.size →
    IsCoprime (StrictHensel.toPolyMod modulus active[i]!)
      (Polynomial.derivative (StrictHensel.toPolyMod modulus f))

theorem hensel_cld_execution_separation_data
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {aTarget : Int32} {output : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hfactors : StrictSelectPrime.SelectionPhysical selection)
    (hleading : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) = (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f selection.factors selection.prime aTarget output)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hnonempty : 0 < f.size) (hpositive : 0 < (SparsePolyZZ.toPoly f).natDegree)
    (hfits : (SparsePolyZZ.toPoly f).natDegree < 2 ^ 31)
    (modulus : Nat) (hmodulus : output.2 = (modulus : Int))
    (cld : Array SparsePolyZZ)
    (hrun : Generated.StrictRecombine.cldPolys ⟨termination⟩ f output.1 modulus = .ok cld) :
    CldSeparationData f output.1 cld modulus := by
  obtain ⟨exponent, he, hp⟩ := hentry.outputModulus_eq_prime_pow
  have hpower : modulus = selection.prime.toNat ^ exponent := by
    exact_mod_cast hmodulus.symm.trans hp
  have hm : 1 < modulus := by
    rw [hpower]
    exact Nat.one_lt_pow (Nat.ne_of_gt he) hselection.goodPrime.prime.one_lt
  have hfBound : StrictHensel.DegreesBound f := by
    intro term ht
    have hb := StrictHensel.canonical_headDominates f hnonempty hcanonical term ht
    rw [getElem!_pos f 0 hnonempty, ← sparsePolyZZ_natDegree_eq_head f hcanonical hnonempty] at hb
    omega
  have hdegrees := hensel_output_head_degree_fits hcount hfactors hleading hselection hentry hfits
  obtain ⟨hsize, hc⟩ := hensel_cld_run_valid hentry modulus hm hfBound hdegrees cld hrun
  have hmonic : ∀ i, i < cld.size → (StrictHensel.toPolyMod modulus output.1[i]!).Monic := by
    intro i hi
    have hio : i < output.1.size := by rwa [← hsize]
    simpa only [getElem!_pos output.1 i hio] using hentry.outputToPolyModMonic modulus i hio
  have hvalues : ∀ i, i < cld.size → ModularCldValue f output.1[i]! cld[i]! modulus := by
    intro i hi
    exact hc i (by rwa [← hsize])
  refine ⟨hsize, hm, hpositive, hvalues, ?_, hmonic, ?_, ?_, ?_, ?_⟩
  · intro i hi
    exact (hvalues i hi).natDegree_lt hm (hmonic i hi) hpositive
  · intro i hi
    exact hensel_output_modulus_degree_pos hcount hfactors hleading hselection hentry
      modulus hm i (by rwa [← hsize])
  · intro i hi
    exact hensel_output_modulus_factor_dvd hcount hfactors hleading hselection hentry
      modulus hmodulus i (by rwa [← hsize])
  · intro i j hi hj hne
    exact hensel_output_modulus_pairwise_coprime hcount hfactors hleading hselection hentry
      modulus hmodulus i j (by rwa [← hsize]) (by rwa [← hsize]) hne
  · intro i hi
    exact hensel_output_modulus_derivative_coprime hcount hfactors hleading hselection hentry
      modulus hmodulus i (by rwa [← hsize])

/-- All upstream local and degree conditions for the complete residual are
now consequences of the execution-derived data and actual truncation. -/
theorem CldSeparationData.residual_properties
    {f : SparsePolyZZ} {active original : Array SparsePolyZZ} {modulus : Nat}
    (data : CldSeparationData f active original modulus)
    (truncated : Array SparsePolyZZ) (divisor : ZZ)
    (htruncate : Generated.StrictVanHoeij.truncateCld original divisor = .ok truncated)
    (size : Nat) (coefficients : Nat → Int) :
    let width := Generated.StrictRecombine.cldSpiralWidth truncated
    let H := rowResidualPolynomial original size width coefficients (modulus : Int)
    width ≤ (SparsePolyZZ.toPoly f).natDegree ∧
      H.natDegree < (SparsePolyZZ.toPoly f).natDegree ∧
      ∀ selected : Fin size, selected.val < original.size →
        (StrictHensel.toPolyMod modulus active[selected.val]! ∣
          Polynomial.map (Int.castRingHom (ZMod modulus)) H -
            Polynomial.C ((coefficients selected.val : Int) : ZMod modulus) *
              Polynomial.derivative (StrictHensel.toPolyMod modulus f)) ∧
        IsCoprime (StrictHensel.toPolyMod modulus active[selected.val]!)
          (Polynomial.derivative (StrictHensel.toPolyMod modulus f)) := by
  dsimp only
  have hw := truncateCld_width_le original truncated divisor htruncate
    (fun i hi => (data.values i hi).canonical) _ data.degrees
  refine ⟨hw, ?_, ?_⟩
  · exact rowResidualPolynomial_natDegree_lt original size _ _ coefficients _
      data.source_degree_pos hw (fun i hi => (data.values i hi).canonical) data.degrees
  · intro selected hs
    refine ⟨rowResidualPolynomial_local_dvd f active original size _ modulus coefficients
      data.values selected hs data.monic data.divides ?_, data.derivative_coprime selected.val hs⟩
    intro i hi hne
    exact data.pairwise selected.val i hs hi (Ne.symm hne)

end Refinement.StrictRecombine
