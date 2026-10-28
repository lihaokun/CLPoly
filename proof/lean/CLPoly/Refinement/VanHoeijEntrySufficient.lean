import CLPoly.Refinement.VanHoeijEntryWidth
import CLPoly.Refinement.VanHoeijEntryExecution
import CLPoly.Refinement.VanHoeijColumnLoopSufficient
import CLPoly.Refinement.VanHoeijZeroWidthValidation
import CLPoly.Refinement.VanHoeijHeight
import CLPoly.Refinement.VanHoeijEntryIrreducibility

set_option autoImplicit false

namespace Refinement.StrictRecombine
open Generated.StrictRecombine Generated.StrictVanHoeij
open Generated.StrictVanHoeijParameters Generated.StrictVanHoeijEntry

/-- The actual entry cannot normally return incomplete above the input-fixed
threshold. All CLD, truncation, width and column-loop facts are obtained from
that entry's execution, including the zero-width branch. -/
theorem tryRecombine_above_input_threshold_complete
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
      (Int32.ofNat lifted.1.size) height < (modulus : Int))
    (after : VanHoeijResult)
    (hrun : tryRecombine (concreteVanHoeijEntryExecution termination)
      f lifted.1 modulus hcountFits = .ok after) :
    after.complete = true := by
  have hcountValue : (Int32.ofNat lifted.1.size).toInt = (lifted.1.size : Int) :=
    Int32.toInt_ofNat_of_lt hcountFits
  have hcountNat : (Int32.ofNat lifted.1.size).toInt.toNat = lifted.1.size := by
    rw [hcountValue]; exact Int.toNat_natCast _
  have hcountNonneg : 0 ≤ (Int32.ofNat lifted.1.size).toInt := by rw [hcountValue]; positivity
  have hheightNonneg : 0 ≤ height := by rw [hheight]; positivity
  obtain ⟨selected, hselected, hresultant, hcoordinates⟩ := selected_parameters_above_input_threshold
    (SparsePolyZZ.toPoly f) (Int32.ofNat lifted.1.size) height (modulus : Int)
    hcountNonneg hheightNonneg hprecision
  rw [hcountNat] at hresultant hcoordinates
  have hsourceFits : (SparsePolyZZ.toPoly f).natDegree < 2 ^ 63 := by omega
  unfold tryRecombine at hrun
  dsimp only at hrun
  split at hrun
  next hguard =>
    cases hh : cldBound f with
    | error fault => simp only [hh, bind, Except.bind] at hrun; contradiction
    | ok actualHeight =>
      simp only [hh, bind, Except.bind] at hrun
      have heHeight : actualHeight = height :=
        (cldBound_polynomial f actualHeight hcanonical hnonempty hfits hh).trans hheight.symm
      subst actualHeight
      cases hs : singletonGroups (Int32.ofNat lifted.1.size) with
      | error fault => simp only [hs, bind, Except.bind] at hrun; contradiction
      | ok groups =>
        simp only [hs, bind, Except.bind] at hrun
        cases hv : validatePartition f lifted.1 modulus groups {} with
        | error fault => simp only [hv, bind, Except.bind] at hrun; contradiction
        | ok checked =>
          rcases checked with ⟨accepted, before⟩
          simp only [hv, bind, Except.bind] at hrun
          cases accepted with
          | true =>
            simp only [↓reduceIte] at hrun
            have he : before = after := Except.ok.inj hrun
            subst before
            exact (validatePartition_success_frame f lifted.1 modulus groups {} after hv).1
          | false =>
            simp only [Bool.false_eq_true, ↓reduceIte] at hrun
            split at hrun
            next hp => simp only [hselected] at hp; contradiction
            next parameters hp =>
              have he : parameters = selected := Option.some.inj (hp.symm.trans hselected)
              subst parameters
              split at hrun
              next fault hc => contradiction
              next original hc =>
                split at hrun
                next fault ht => contradiction
                next truncated ht =>
                  have data := hensel_cld_execution_separation_data hcount hfactors hleading
                    hselection hentry hcanonical hnonempty hdegree hfits modulus hmodulus original hc
                  have hw := data.truncated_width_toInt truncated selected.divisor ht hfits
                  by_cases hz : cldSpiralWidth truncated = 0
                  · obtain ⟨result, hsuccess⟩ := hensel_zero_width_singleton_validation_complete
                      hcount hfactors hleading hselection hentry modulus hmodulus original truncated
                      data selected.divisor ht hz (Int32.ofNat lifted.1.size) groups hcountNat hs
                      hcanonical hprimitive hpositive hsourceFits hresultant hcoordinates {}
                    have he := (Prod.mk.inj (Except.ok.inj (hv.symm.trans hsuccess))).1
                    cases he
                  · have htarget : 0 < min (2 : Int32) (cldWidth truncated) := by
                      rw [Int32.lt_iff_toInt_lt, Int32.toInt_zero, int32_min_toInt, hw]
                      change (0 : Int) < min 2 (cldSpiralWidth truncated : Int)
                      omega
                    have hwidth : (cldWidth truncated).toInt.toNat = cldSpiralWidth truncated := by
                      rw [hw]; exact Int.toNat_natCast _
                    exact columnLoop_above_input_threshold_complete hcount hfactors hleading
                      hselection hentry hcountFits.le modulus hmodulus original truncated data
                      height (Int32.ofNat lifted.1.size) (cldWidth truncated) selected hcountValue
                      hselected ht hwidth _ _ hheight hcanonical hprimitive hpositive hsourceFits
                      hprecision _ _ htarget before after hrun
  next hguard => contradiction

/-- At sufficient precision, every normal return of the actual entry is a
complete irreducible factorization. Completeness is proved, not assumed. -/
theorem tryRecombine_above_input_threshold_irreducible
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
      (Int32.ofNat lifted.1.size) height < (modulus : Int))
    (after : VanHoeijResult)
    (hrun : tryRecombine (concreteVanHoeijEntryExecution termination)
      f lifted.1 modulus hcountFits = .ok after) :
    after.complete = true ∧ factorArrayProduct after.factors = SparsePolyZZ.toPoly f ∧
      ∀ factor ∈ after.factors.toList, Irreducible (SparsePolyZZ.toPoly factor) := by
  have hc := tryRecombine_above_input_threshold_complete hcount hfactors hleading hselection
    hentry hcountFits modulus hmodulus hcanonical hnonempty hdegree hfits hprimitive
    hpositive height hheight hprecision after hrun
  exact ⟨hc, hensel_tryRecombine_complete_irreducible hcount hfactors hleading hselection
    hentry modulus hmodulus hcanonical hnonempty hfits hcountFits hprimitive after hc hrun⟩

end Refinement.StrictRecombine
