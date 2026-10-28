import CLPoly.Refinement.VanHoeijControllerTotal
import CLPoly.Refinement.VanHoeijCldExecution
import CLPoly.Refinement.VanHoeijTruncateExecution
import CLPoly.Refinement.VanHoeijEntrySufficient

set_option autoImplicit false

namespace Refinement.StrictRecombine
open _root_.Generated.StrictVanHoeij _root_.Generated.StrictVanHoeijEntry

theorem cldBound_succeeds (f : SparsePolyZZ)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hnonempty : 0 < f.size) (hdegree : 0 < (SparsePolyZZ.toPoly f).natDegree)
    (hfits : (SparsePolyZZ.toPoly f).natDegree < 2 ^ 31) :
    ∃ height, cldBound f = .ok height := by
  have hcast : (get_deg f).toInt32 = Int32.ofNat (SparsePolyZZ.toPoly f).natDegree := by
    rw [canonical_get_deg_eq_ofNat f hcanonical hnonempty, Int64.toInt32_ofNat']
  have hp : 0 < (get_deg f).toInt32 := by
    rw [hcast, Int32.lt_iff_toInt_lt, Int32.toInt_zero, Int32.toInt_ofNat_of_lt hfits]
    exact_mod_cast hdegree
  unfold cldBound
  rw [if_pos hp]
  exact ⟨_, rfl⟩

theorem hensel_tryRecombine_succeeds
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {aTarget : Int32} {lifted : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hfactors : StrictSelectPrime.SelectionPhysical selection)
    (hleading : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) = (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect StrictHensel.concreteDivmodTermination
      f selection.factors selection.prime aTarget lifted)
    (hcountFits : lifted.1.size < 2 ^ 31)
    (modulus : Nat) (hmodulus : lifted.2 = (modulus : Int))
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hnonempty : 0 < f.size) (hdegree : 0 < (SparsePolyZZ.toPoly f).natDegree)
    (hfits : (SparsePolyZZ.toPoly f).natDegree < 2 ^ 31) :
    ∃ after, tryRecombine (concreteVanHoeijEntryExecution StrictHensel.concreteDivmodTermination)
      f lifted.1 modulus hcountFits = .ok after := by
  have hs := (StrictFactorZZ.selectionHenselFactors_pointwise_associated hcount hfactors
    hleading hselection hentry).1
  have hcountValue : (Int32.ofNat lifted.1.size).toInt = (lifted.1.size : Int) :=
    Int32.toInt_ofNat_of_lt hcountFits
  have hcountNat : (Int32.ofNat lifted.1.size).toInt.toNat = lifted.1.size := by
    rw [hcountValue]; exact Int.toNat_natCast _
  obtain ⟨exponent, he, hp⟩ := hentry.outputModulus_eq_prime_pow
  have hpower : modulus = selection.prime.toNat ^ exponent := by
    exact_mod_cast hmodulus.symm.trans hp
  have hm : 1 < modulus := by
    rw [hpower]
    exact Nat.one_lt_pow (Nat.ne_of_gt he) hselection.goodPrime.prime.one_lt
  have hguard : 2 ≤ Int32.ofNat lifted.1.size ∧ (1 : Int) < modulus := by
    constructor
    · rw [Int32.le_iff_toInt_le, hcountValue]
      change (2 : Int) ≤ (lifted.1.size : Int)
      exact_mod_cast (by omega : 2 ≤ lifted.1.size)
    · exact_mod_cast hm
  obtain ⟨height, hh⟩ := cldBound_succeeds f hcanonical hnonempty hdegree hfits
  obtain ⟨singletons, hsingleton, _, _⟩ := singletonGroups_succeeds (Int32.ofNat lifted.1.size) (by
    rw [Int32.le_iff_toInt_le, Int32.toInt_zero, hcountValue]; positivity)
  have hvalid := singletonGroups_indices_valid (Int32.ofNat lifted.1.size) singletons hsingleton
  rw [hcountNat] at hvalid
  obtain ⟨checked, hv⟩ := validatePartition_succeeds f lifted.1 modulus (by omega)
    singletons hvalid hcanonical {}
  unfold tryRecombine
  dsimp only
  rw [if_pos hguard]
  simp only [hh, hsingleton, hv, bind, Except.bind]
  split
  · exact ⟨checked.2, rfl⟩
  · split
    next hp => exact ⟨checked.2, rfl⟩
    next parameters hparams =>
      obtain ⟨original, hc⟩ := hensel_cldPolys_succeeds hcount hfactors hleading hselection hentry
        hcanonical hnonempty hfits modulus hmodulus
      have hd := (selectParameters_spec (Int32.ofNat lifted.1.size) height modulus parameters hparams).2.1
      obtain ⟨truncated, ht⟩ := truncateCld_succeeds original parameters.divisor hd
      split
      next fault hbad =>
        have he := hc.symm.trans hbad
        contradiction
      next original' hc' =>
        have he := Except.ok.inj (hc'.symm.trans hc)
        subst original'
        split
        next fault hbad => rw [ht] at hbad; contradiction
        next truncated' ht' =>
          have he := Except.ok.inj (ht'.symm.trans ht)
          subst truncated'
          apply columnLoop_succeeds f lifted.1 truncated modulus (by omega)
            parameters.modulus 32 parameters.residual (Int32.ofNat lifted.1.size)
            _ _ _ (by decide) _ hcountNat (by omega) (by omega) hcanonical

theorem tryRecombine_above_input_threshold_exists_irreducible
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {aTarget : Int32} {lifted : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hfactors : StrictSelectPrime.SelectionPhysical selection)
    (hleading : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) = (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect StrictHensel.concreteDivmodTermination
      f selection.factors selection.prime aTarget lifted)
    (hcountFits : lifted.1.size < 2 ^ 31)
    (modulus : Nat) (hmodulus : lifted.2 = (modulus : Int))
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hnonempty : 0 < f.size) (hdegree : 0 < (SparsePolyZZ.toPoly f).natDegree)
    (hfits : (SparsePolyZZ.toPoly f).natDegree < 2 ^ 31)
    (hprimitive : (SparsePolyZZ.toPoly f).IsPrimitive)
    (hpositive : (0 : Int) < (SparsePolyZZ.toPoly f).leadingCoeff)
    (height : Int)
    (hheight : height = ((SparsePolyZZ.toPoly f).natDegree *
      2 ^ ((SparsePolyZZ.toPoly f).natDegree - 1) *
      (Finset.range ((SparsePolyZZ.toPoly f).natDegree + 1)).sum
        (fun k => ((SparsePolyZZ.toPoly f).coeff k).natAbs) : Nat))
    (hprecision : vanHoeijInputThreshold (SparsePolyZZ.toPoly f)
      (Int32.ofNat lifted.1.size) height < (modulus : Int)) :
    ∃ after, tryRecombine (concreteVanHoeijEntryExecution StrictHensel.concreteDivmodTermination)
      f lifted.1 modulus hcountFits = .ok after ∧ after.complete = true ∧
      factorArrayProduct after.factors = SparsePolyZZ.toPoly f ∧
      ∀ factor ∈ after.factors.toList, Irreducible (SparsePolyZZ.toPoly factor) := by
  obtain ⟨after, ha⟩ := hensel_tryRecombine_succeeds hcount hfactors hleading hselection hentry
    hcountFits modulus hmodulus hcanonical hnonempty hdegree hfits
  exact ⟨after, ha, tryRecombine_above_input_threshold_irreducible hcount hfactors hleading
    hselection hentry hcountFits modulus hmodulus hcanonical hnonempty hdegree hfits hprimitive
    hpositive height hheight hprecision after ha⟩

end Refinement.StrictRecombine
