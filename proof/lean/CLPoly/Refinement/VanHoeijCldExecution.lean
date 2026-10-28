import CLPoly.Refinement.VanHoeijDivmodComplete
import CLPoly.Refinement.VanHoeijHenselLocal

set_option autoImplicit false

namespace Refinement.StrictRecombine

open _root_.Generated.StrictRecombine

theorem cldPolysLoop_succeeds
    (f : SparsePolyZZ) (active : Array SparsePolyZZ) (m : Nat) (hm : 1 < m)
    (hf : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hfBound : StrictHensel.DegreesBound f)
    (hcanonical : ∀ i (hi : i < active.size),
      StrictPolynomialMod.SparsePolyZZCanonical active[i])
    (hnonempty : ∀ i (hi : i < active.size), 0 < active[i].size)
    (hone : ∀ i (hi : i < active.size), active[i][0]!.2 = 1)
    (hbound : ∀ i (hi : i < active.size), StrictHensel.DegreesBound active[i])
    (hdivides : ∀ i (hi : i < active.size),
      StrictHensel.toPolyMod m active[i] ∣ StrictHensel.toPolyMod m f)
    (index : Nat) (result : Array SparsePolyZZ) :
    ∃ output, cldPolysLoop ⟨StrictHensel.concreteDivmodTermination⟩
      f active (m : Int) index result = .ok output := by
  have hmInt : (0 : Int) < m := by exact_mod_cast (by omega : 0 < m)
  induction hmeasure : active.size - index using Nat.strong_induction_on
      generalizing index result with
  | h measure ih =>
    rw [cldPolysLoop]
    split
    next hi =>
      obtain ⟨q, hq⟩ := StrictHensel.monic_divmod_complete_of_dvd f active[index] m hm
        hf (hcanonical index hi) (hnonempty index hi) (hone index hi)
        hfBound (hbound index hi) (hdivides index hi)
      rw [hq]
      simp only [Array.isEmpty_empty, ↓reduceDIte, derivativeZZRaw]
      obtain ⟨derivative, hd⟩ := modCoeffLoop_complete
        (derivativeZZLoop active[index] 0 #[]) (m : Int) 0 #[] (ne_of_gt hmInt)
      rw [hd]
      dsimp only
      obtain ⟨product, hp⟩ := multiplyNormalizeModRaw_complete q derivative (m : Int) (ne_of_gt hmInt)
      rw [hp]
      dsimp only
      obtain ⟨cld, hc⟩ := symmetricModRaw_complete product (m : Int) hmInt
      rw [hc]
      dsimp only
      exact ih _ (by omega) (index + 1) (result.push cld) rfl
    next hi => exact ⟨result, rfl⟩

theorem hensel_cldPolys_succeeds
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {aTarget : Int32} {output : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hfactors : StrictSelectPrime.SelectionPhysical selection)
    (hleading : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) = (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect StrictHensel.concreteDivmodTermination
      f selection.factors selection.prime aTarget output)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hnonempty : 0 < f.size)
    (hfits : (SparsePolyZZ.toPoly f).natDegree < 2 ^ 31)
    (modulus : Nat) (hmodulus : output.2 = (modulus : Int)) :
    ∃ cld, cldPolys ⟨StrictHensel.concreteDivmodTermination⟩
      f output.1 (modulus : Int) = .ok cld := by
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
  have hphysical : ∀ i (hi : i < output.1.size),
      0 < output.1[i].size ∧ output.1[i][0]!.2 = 1 := by
    intro i hi
    obtain ⟨head, tail, hlist⟩ := hentry.outputOneHead i hi
    have hlength := congrArg List.length hlist
    have hsEq : output.1[i].size = tail.length + 1 := by simpa using hlength
    have hs : 0 < output.1[i].size := by omega
    refine ⟨hs, ?_⟩
    have hget := Array.getElem_toList hs
    rw [getElem!_pos output.1[i] 0 hs]
    simpa [hlist] using (congrArg Prod.snd hget).symm
  have hdegrees := hensel_output_head_degree_fits hcount hfactors hleading hselection hentry hfits
  apply cldPolysLoop_succeeds f output.1 modulus hm hcanonical hfBound
    hentry.outputCanonical (fun i hi => (hphysical i hi).1) (fun i hi => (hphysical i hi).2)
  · intro i hi term ht
    have hb := StrictHensel.canonical_headDominates _ (hphysical i hi).1
      (hentry.outputCanonical i hi) term ht
    have hd := hdegrees i hi
    omega
  · intro i hi
    simpa only [getElem!_pos output.1 i hi] using
      hensel_output_modulus_factor_dvd hcount hfactors hleading hselection hentry
        modulus hmodulus i hi

end Refinement.StrictRecombine
