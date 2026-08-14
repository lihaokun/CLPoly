/-
  No strict recombination refinement theorem is exported yet.

  The former theorem equated the generated C++ result with an arbitrary UFD
  existence witness and contained a placeholder proof.  Both the witness-based
  implementation and the claimed refinement have been removed.  A replacement
  must prove the generated candidate/extraction loops directly.
-/
import CLPoly.Algorithm.Recombine
import CLPoly.Generated.StrictRecombine
import CLPoly.Refinement.Basic
import CLPoly.Refinement.Hensel
import Batteries.Data.Array.Lemmas

set_option autoImplicit false

open Polynomial
open CLPoly.Math

namespace Refinement.StrictRecombine

noncomputable def factorArrayToL2 (factors : Array SparsePolyZZ) :
    List (Polynomial Int) :=
  factors.toList.map SparsePolyZZ.toPoly

def ValidActiveIndices (active : Array Int32)
    (lifted : Array SparsePolyZZ) : Prop :=
  ∀ index (hindex : index < active.size),
    0 ≤ active[index]'hindex ∧
      (active[index]'hindex).toInt64.toNat < lifted.size

private def gatherSuffix (active : Array Int32)
    (lifted : Array SparsePolyZZ) (index : Nat) : Array SparsePolyZZ :=
  ((active.toList.drop index).map
    (fun sourceIndex => lifted[sourceIndex.toInt64.toNat]!)).toArray

theorem gatherActiveLoop_refines (active : Array Int32)
    (lifted : Array SparsePolyZZ) (index : Nat)
    (result : Array SparsePolyZZ)
    (hvalid : ValidActiveIndices active lifted) :
    Generated.StrictRecombine.gatherActiveLoop active lifted index result =
      .ok (result ++ gatherSuffix active lifted index) := by
  induction hmeasure : active.size - index using Nat.strong_induction_on
      generalizing index result with
  | h measure ih =>
      rw [Generated.StrictRecombine.gatherActiveLoop]
      split
      next hindex =>
        have hentry := hvalid index hindex
        rw [dif_pos hentry.1, dif_pos hentry.2]
        rw [ih (active.size - (index + 1)) (by omega)
          (index + 1) (result.push lifted[active[index].toInt64.toNat]) rfl]
        have hsuffix : active.toList.drop index = active[index] ::
            active.toList.drop (index + 1) := by
          simpa using List.drop_eq_getElem_cons
            (l := active.toList) (i := index) (by simpa using hindex)
        simp [gatherSuffix, hsuffix, Array.push, hentry.2]
      next hindex =>
        have hle : active.size ≤ index := Nat.le_of_not_gt hindex
        simp [gatherSuffix, List.drop_eq_nil_iff.mpr hle]

theorem gatherActive_refines (active : Array Int32)
    (lifted : Array SparsePolyZZ) (hvalid : ValidActiveIndices active lifted) :
    ∃ output,
      Generated.StrictRecombine.gatherActive active lifted = .ok output ∧
      factorArrayToL2 output = active.toList.map
        (fun sourceIndex => SparsePolyZZ.toPoly
          lifted[sourceIndex.toInt64.toNat]!) := by
  refine ⟨((active.toList.map
      (fun sourceIndex => lifted[sourceIndex.toInt64.toNat]!)).toArray), ?_, ?_⟩
  · simpa [Generated.StrictRecombine.gatherActive, gatherSuffix] using
      gatherActiveLoop_refines active lifted 0 #[] hvalid
  · simp [factorArrayToL2]

theorem gatherActive_size_of_success (active : Array Int32)
    (lifted output : Array SparsePolyZZ)
    (hrun : Generated.StrictRecombine.gatherActive active lifted = .ok output) :
    output.size = active.size := by
  unfold Generated.StrictRecombine.gatherActive at hrun
  have loopSize : ∀ index result,
      Generated.StrictRecombine.gatherActiveLoop active lifted index result =
        .ok output → output.size = result.size + (active.size - index) := by
    intro index result hloop
    induction hmeasure : active.size - index using Nat.strong_induction_on
        generalizing index result with
    | h measure ih =>
        rw [Generated.StrictRecombine.gatherActiveLoop] at hloop
        split at hloop
        next hindex =>
          dsimp at hloop
          split at hloop
          next hnonnegative =>
            split at hloop
            next hsource =>
              have htail := ih (active.size - (index + 1)) (by omega)
                (index + 1) (result.push lifted[active[index].toInt64.toNat])
                hloop rfl
              simp only [Array.size_push] at htail
              omega
            next hsource => simp at hloop
          next hnonnegative => simp at hloop
        next hindex =>
          have hle : active.size ≤ index := Nat.le_of_not_gt hindex
          have hout := Except.ok.inj hloop
          subst output
          have : measure = 0 := by omega
          simp [this]
  simpa using (loopSize 0 #[] hrun)

theorem appendFallbackLoop_refines (fallback : Array SparsePolyZZ)
    (index : Nat) (result : Array SparsePolyZZ) :
    Generated.StrictRecombine.appendFallbackLoop fallback index result =
      .ok (result ++ (fallback.toList.drop index).toArray) := by
  induction hmeasure : fallback.size - index using Nat.strong_induction_on
      generalizing index result with
  | h measure ih =>
      rw [Generated.StrictRecombine.appendFallbackLoop]
      split
      next hindex =>
        rw [ih (fallback.size - (index + 1)) (by omega)
          (index + 1) (result.push fallback[index]) rfl]
        have hsuffix : fallback.toList.drop index = fallback[index] ::
            fallback.toList.drop (index + 1) := by
          simpa using List.drop_eq_getElem_cons
            (l := fallback.toList) (i := index) (by simpa using hindex)
        simp [hsuffix, Array.push]
      next hindex =>
        have hle : fallback.size ≤ index := Nat.le_of_not_gt hindex
        simp [List.drop_eq_nil_iff.mpr hle]

theorem appendFallback_refines (fallback result : Array SparsePolyZZ) :
    ∃ output,
      Generated.StrictRecombine.appendFallback fallback result = .ok output ∧
      factorArrayToL2 output =
        factorArrayToL2 result ++ factorArrayToL2 fallback := by
  refine ⟨result ++ fallback, ?_, ?_⟩
  · simpa [Generated.StrictRecombine.appendFallback] using
      appendFallbackLoop_refines fallback 0 result
  · simp [factorArrayToL2]

/-- Pure list meaning of the source reverse-erasure loop. -/
def removeConsumedL2 (active : Array Int32) (consumed : Array Bool) :
    List Int32 :=
  (active.toList.zip consumed.toList).filterMap fun item =>
    if item.2 then none else some item.1

private theorem removeConsumedLoop_refines_prefix
    (consumed : Array Bool) (remaining : Nat) (active : Array Int32)
    (hremaining : remaining ≤ consumed.size)
    (hactiveRemaining : remaining ≤ active.size) :
    ∃ output,
      Generated.StrictRecombine.removeConsumedLoop consumed remaining active =
        .ok output ∧
      output.size ≤ active.size := by
  induction remaining generalizing active with
  | zero =>
      refine ⟨active, ?_, Nat.le_refl _⟩
      rw [Generated.StrictRecombine.removeConsumedLoop]
  | succ remaining ih =>
      rw [Generated.StrictRecombine.removeConsumedLoop]
      have hindex : remaining < consumed.size := by omega
      rw [dif_pos hindex]
      split
      next hmarked =>
        have hactive : remaining < active.size := by omega
        rw [dif_pos hactive]
        have heraseSize : (active.eraseIdxIfInBounds remaining).size =
            active.size - 1 := by
          simp [hactive]
        have hremaining' : remaining ≤
            (active.eraseIdxIfInBounds remaining).size := by
          rw [heraseSize]
          omega
        rcases ih (active.eraseIdxIfInBounds remaining) (by omega) hremaining' with
          ⟨output, hrun, hsize⟩
        exact ⟨output, hrun, hsize.trans (by rw [heraseSize]; omega)⟩
      next hkept =>
        rcases ih active (by omega) (by omega) with ⟨output, hrun, hsize⟩
        exact ⟨output, hrun, hsize⟩

private theorem removeConsumedLoop_strict_of_marked
    (consumed : Array Bool) (remaining : Nat) (active output : Array Int32)
    (hremaining : remaining ≤ consumed.size)
    (hactiveRemaining : remaining ≤ active.size)
    (hmarked : ∃ index, ∃ hindex : index < remaining,
      consumed[index] = true)
    (hrun : Generated.StrictRecombine.removeConsumedLoop consumed remaining active =
      .ok output) :
    output.size < active.size := by
  induction remaining generalizing active output with
  | zero =>
      rcases hmarked with ⟨index, hindex, _⟩
      omega
  | succ remaining ih =>
      rw [Generated.StrictRecombine.removeConsumedLoop] at hrun
      have hindexBound : remaining < consumed.size := by omega
      rw [dif_pos hindexBound] at hrun
      by_cases hlast : consumed[remaining] = true
      · rw [if_pos hlast] at hrun
        have hactive : remaining < active.size := by omega
        rw [dif_pos hactive] at hrun
        have htailLe : output.size ≤
            (active.eraseIdxIfInBounds remaining).size := by
          rcases removeConsumedLoop_refines_prefix consumed remaining
              (active.eraseIdxIfInBounds remaining) (by omega) (by
                simp [hactive]; omega) with ⟨tail, htailRun, htailSize⟩
          rw [hrun] at htailRun
          exact Except.ok.inj htailRun ▸ htailSize
        simpa [hactive] using lt_of_le_of_lt htailLe (by
          simp [hactive]
          omega)
      · rw [if_neg hlast] at hrun
        have hmarkedPrefix : ∃ index, ∃ hindex : index < remaining,
            consumed[index] = true := by
          rcases hmarked with ⟨index, hindex, hvalue⟩
          refine ⟨index, ?_, hvalue⟩
          have hne : index ≠ remaining := by
            intro heq
            subst index
            exact hlast hvalue
          omega
        exact ih active output (by omega) (by omega) hmarkedPrefix hrun

theorem removeConsumed_succeeds (active : Array Int32)
    (consumed : Array Bool) (hsizes : consumed.size = active.size) :
    ∃ output,
      Generated.StrictRecombine.removeConsumed active consumed = .ok output ∧
      output.size ≤ active.size := by
  unfold Generated.StrictRecombine.removeConsumed
  rw [dif_pos hsizes]
  exact removeConsumedLoop_refines_prefix consumed active.size active
    (by omega) (Nat.le_refl _)

theorem removeConsumed_strict_of_marked (active : Array Int32)
    (consumed : Array Bool) (hsizes : consumed.size = active.size)
    (hmarked : ∃ index, ∃ hindex : index < consumed.size,
      consumed[index] = true)
    (output : Array Int32)
    (hrun : Generated.StrictRecombine.removeConsumed active consumed =
      .ok output) :
    output.size < active.size := by
  unfold Generated.StrictRecombine.removeConsumed at hrun
  rw [dif_pos hsizes] at hrun
  apply removeConsumedLoop_strict_of_marked consumed active.size active output
    (by omega) (Nat.le_refl _) (by simpa [hsizes] using hmarked) hrun

def CandidateIndicesValid (candidate : Array Int32)
    (consumed : Array Bool) : Prop :=
  ∀ index (hindex : index < candidate.size),
    0 ≤ candidate[index] ∧
      candidate[index].toInt64.toNat < consumed.size

theorem candidateAvailableLoop_succeeds (candidate : Array Int32)
    (consumed : Array Bool) (index : Nat)
    (hvalid : CandidateIndicesValid candidate consumed) :
    ∃ available,
      Generated.StrictRecombine.candidateAvailableLoop candidate consumed index =
        .ok available := by
  induction hmeasure : candidate.size - index using Nat.strong_induction_on
      generalizing index with
  | h measure ih =>
      rw [Generated.StrictRecombine.candidateAvailableLoop]
      split
      next hindex =>
        have hentry := hvalid index hindex
        rw [dif_pos hentry.1, dif_pos hentry.2]
        split
        next hconsumed => exact ⟨false, rfl⟩
        next hfree =>
          exact ih (candidate.size - (index + 1)) (by omega)
            (index + 1) rfl
      next hindex => exact ⟨true, rfl⟩

theorem markConsumedLoop_succeeds_size (candidate : Array Int32)
    (consumed : Array Bool) (index : Nat)
    (hvalid : CandidateIndicesValid candidate consumed) :
    ∃ output,
      Generated.StrictRecombine.markConsumedLoop candidate index consumed =
        .ok output ∧ output.size = consumed.size := by
  induction hmeasure : candidate.size - index using Nat.strong_induction_on
      generalizing index consumed with
  | h measure ih =>
      rw [Generated.StrictRecombine.markConsumedLoop]
      split
      next hindex =>
        have hentry := hvalid index hindex
        let activeNat := candidate[index].toInt64.toNat
        rw [dif_pos hentry.1, dif_pos hentry.2]
        have hvalidSet : CandidateIndicesValid candidate
            (consumed.set activeNat true) := by
          intro next hnext
          have hv := hvalid next hnext
          exact ⟨hv.1, by simpa using hv.2⟩
        rcases ih (candidate.size - (index + 1)) (by omega)
            (consumed.set activeNat true) (index + 1) hvalidSet rfl with
          ⟨output, hrun, hsize⟩
        refine ⟨output, hrun, ?_⟩
        simpa using hsize
      next hindex => exact ⟨consumed, rfl, rfl⟩

theorem markConsumedLoop_size_of_success (candidate : Array Int32)
    (index : Nat) (consumed output : Array Bool)
    (hrun : Generated.StrictRecombine.markConsumedLoop candidate index consumed =
      .ok output) :
    output.size = consumed.size := by
  induction hmeasure : candidate.size - index using Nat.strong_induction_on
      generalizing index consumed output with
  | h measure ih =>
      rw [Generated.StrictRecombine.markConsumedLoop] at hrun
      split at hrun
      next hindex =>
        dsimp at hrun
        split at hrun
        next hnonnegative =>
          split at hrun
          next hactive =>
            have htail := ih (candidate.size - (index + 1)) (by omega)
              (index + 1)
              (consumed.set candidate[index].toInt64.toNat true) output
              hrun rfl
            simpa using htail
          next hactive => contradiction
        next hnonnegative => contradiction
      next hindex => exact Except.ok.inj hrun ▸ rfl

theorem validateCandidatesLoop_consumed_size
    (ops : Generated.StrictRecombine.CandidateValidationRawOps)
    (candidates : Array (Array Int32)) (candidateIndex : Nat)
    (activeLifted : Array SparsePolyZZ) (modulus : ZZ)
    (fStar fStar' : SparsePolyZZ) (result result' : Array SparsePolyZZ)
    (consumed consumed' : Array Bool) (remaining : Nat)
    (hrun : Generated.StrictRecombine.validateCandidatesLoop ops candidates
      candidateIndex activeLifted modulus fStar result consumed remaining =
        .ok (fStar', result', consumed')) :
    consumed'.size = consumed.size := by
  induction hmeasure : candidates.size - candidateIndex using Nat.strong_induction_on
      generalizing candidateIndex fStar result consumed remaining fStar' result' consumed' with
  | h measure ih =>
      rw [Generated.StrictRecombine.validateCandidatesLoop] at hrun
      split at hrun
      next hcandidates =>
        dsimp at hrun
        split at hrun
        next hempty =>
          exact ih (candidates.size - (candidateIndex + 1)) (by omega)
            (candidateIndex := candidateIndex + 1) (fStar := fStar)
            (result := result) (consumed := consumed) (remaining := remaining)
            (fStar' := fStar') (result' := result') (consumed' := consumed') hrun rfl
        next hempty =>
          split at hrun
          next htrivial =>
            exact ih (candidates.size - (candidateIndex + 1)) (by omega)
              (candidateIndex := candidateIndex + 1) (fStar := fStar)
              (result := result) (consumed := consumed) (remaining := remaining)
              (fStar' := fStar') (result' := result') (consumed' := consumed') hrun rfl
          next hnontrivial =>
            cases havailable : Generated.StrictRecombine.candidateAvailable
                candidates[candidateIndex] consumed with
            | error fault => simp [havailable] at hrun
            | ok available =>
              cases available with
              | false =>
                simp only [havailable] at hrun
                exact ih (candidates.size - (candidateIndex + 1)) (by omega)
                  (candidateIndex := candidateIndex + 1) (fStar := fStar)
                  (result := result) (consumed := consumed) (remaining := remaining)
                  (fStar' := fStar') (result' := result')
                  (consumed' := consumed') hrun rfl
              | true =>
                simp only [havailable] at hrun
                split at hrun
                next hfstar =>
                  cases hproduct : Generated.StrictRecombine.trialProductLoop
                      ops.product candidates[candidateIndex] activeLifted modulus 0
                      #[(⟨0⟩, fStar[0].2)] with
                  | error fault => simp [hproduct] at hrun
                  | ok product =>
                    simp only [hproduct] at hrun
                    cases hsymmetric : ops.symmetricMod product modulus with
                    | error fault => simp [hsymmetric] at hrun
                    | ok symmetric =>
                      simp only [hsymmetric] at hrun
                      cases hprimitive : ops.primitive symmetric with
                      | error fault => simp [hprimitive] at hrun
                      | ok primitiveResult =>
                        rcases primitiveResult with ⟨content, factor⟩
                        simp only [hprimitive] at hrun
                        cases hdivmod : ops.divmod fStar factor with
                        | error fault => simp [hdivmod] at hrun
                        | ok divResult =>
                          rcases divResult with ⟨quotient, remainder⟩
                          simp only [hdivmod] at hrun
                          by_cases hremainder : remainder.isEmpty = true
                          · simp only [hremainder, if_true] at hrun
                            cases hquotientPrimitive : ops.primitive quotient with
                            | error fault => simp [hquotientPrimitive] at hrun
                            | ok quotientResult =>
                              rcases quotientResult with ⟨quotientContent,
                                quotientPrimitive⟩
                              simp only [hquotientPrimitive] at hrun
                              cases hmark : Generated.StrictRecombine.markConsumedLoop
                                  candidates[candidateIndex] 0 consumed with
                              | error fault => simp [hmark] at hrun
                              | ok consumedNext =>
                                simp only [hmark] at hrun
                                have htail := ih
                                  (candidates.size - (candidateIndex + 1))
                                  (by omega) (candidateIndex := candidateIndex + 1)
                                  (fStar := quotientPrimitive)
                                  (result := result.push factor)
                                  (consumed := consumedNext)
                                  (remaining := remaining - candidates[candidateIndex].size)
                                  (fStar' := fStar') (result' := result')
                                  (consumed' := consumed') hrun rfl
                                exact htail.trans
                                  (markConsumedLoop_size_of_success
                                    candidates[candidateIndex] 0 consumed
                                    consumedNext hmark)
                          · simp only [hremainder, if_false] at hrun
                            exact ih (candidates.size - (candidateIndex + 1))
                              (by omega) (candidateIndex := candidateIndex + 1)
                              (fStar := fStar) (result := result)
                              (consumed := consumed) (remaining := remaining)
                              (fStar' := fStar') (result' := result')
                              (consumed' := consumed') hrun rfl
                next hfstar => contradiction
      next hcandidates =>
        exact (congrArg (fun output => output.2.2.size)
          (Except.ok.inj hrun)).symm

theorem validateCandidates_consumed_size
    (ops : Generated.StrictRecombine.CandidateValidationRawOps)
    (fStar : SparsePolyZZ) (activeLifted : Array SparsePolyZZ) (modulus : ZZ)
    (candidates : Array (Array Int32)) (result : Array SparsePolyZZ)
    (fStar' : SparsePolyZZ) (result' : Array SparsePolyZZ)
    (consumed : Array Bool)
    (hrun : Generated.StrictRecombine.validateCandidates ops fStar activeLifted
      modulus candidates result = .ok (fStar', result', consumed)) :
    consumed.size = activeLifted.size := by
  unfold Generated.StrictRecombine.validateCandidates at hrun
  exact (validateCandidatesLoop_consumed_size ops candidates 0 activeLifted
    modulus fStar fStar' result result'
    (Array.replicate activeLifted.size false) consumed activeLifted.size hrun).trans
      (by simp)

/-- The concrete gather and validation loops themselves prove that every
successful extraction removes at least one active entry.  No external length
or semantic termination oracle is required. -/
def removalTermination (ops : Generated.StrictRecombine.VanHoeijRawOps) :
    Generated.StrictRecombine.VanHoeijTermination ops := {
  extraction_decreases := by
    intro lifted modulus state activeLifted candidates fStar' result' consumed
      active' hgather hvalidate hfound hremove
    exact removeConsumed_strict_of_marked state.active consumed
      ((validateCandidates_consumed_size ops.validation state.fStar activeLifted
        modulus candidates state.result fStar' result' consumed hvalidate).trans
        (gatherActive_size_of_success state.active lifted activeLifted hgather))
      hfound active' hremove }

/-- Algebraic execution contract for one actual multiply/normalize/mod step.
It relates concrete successful output modulo `m`; it does not supply or choose
a factorization. -/
def TrialProductStepCorrect
    (ops : Generated.StrictRecombine.TrialProductRawOps) : Prop :=
  ∀ left right (modulus : Nat) output, 0 < modulus →
    ops.multiplyNormalizeMod left right (modulus : ZZ) = .ok output →
    Refinement.StrictHensel.toPolyMod modulus output =
      Refinement.StrictHensel.toPolyMod modulus left *
        Refinement.StrictHensel.toPolyMod modulus right

noncomputable def SelectedProductMod (modulus : Nat) (candidate : Array Int32)
    (activeLifted : Array SparsePolyZZ) (index : Nat) :
    Polynomial (ZMod modulus) :=
  ((candidate.toList.drop index).map fun activeIndex =>
    Refinement.StrictHensel.toPolyMod modulus
      activeLifted[activeIndex.toInt64.toNat]!).prod

theorem trialProductLoop_refines
    (ops : Generated.StrictRecombine.TrialProductRawOps)
    (hstep : TrialProductStepCorrect ops)
    (candidate : Array Int32) (activeLifted : Array SparsePolyZZ)
    (modulus : Nat) (hmodulus : 0 < modulus) (index : Nat)
    (product output : SparsePolyZZ)
    (hvalid : CandidateIndicesValid candidate
      (Array.replicate activeLifted.size false))
    (hrun : Generated.StrictRecombine.trialProductLoop ops candidate
      activeLifted (modulus : ZZ) index product = .ok output) :
    Refinement.StrictHensel.toPolyMod modulus output =
      Refinement.StrictHensel.toPolyMod modulus product *
        SelectedProductMod modulus candidate activeLifted index := by
  induction hmeasure : candidate.size - index using Nat.strong_induction_on
      generalizing index product output with
  | h measure ih =>
      rw [Generated.StrictRecombine.trialProductLoop] at hrun
      split at hrun
      next hindex =>
        have hentry := hvalid index hindex
        have hactive : candidate[index].toInt64.toNat < activeLifted.size := by
          simpa using hentry.2
        rw [dif_pos hentry.1, dif_pos hactive] at hrun
        split at hrun
        next fault hcall => contradiction
        next product' hcall =>
          let selected := activeLifted[candidate[index].toInt64.toNat]
          have htail := ih (candidate.size - (index + 1)) (by omega)
            (index + 1) product' output hrun rfl
          have hsuffix : candidate.toList.drop index = candidate[index] ::
              candidate.toList.drop (index + 1) := by
            simpa using List.drop_eq_getElem_cons
              (l := candidate.toList) (i := index) (by simpa using hindex)
          have hstep' := hstep product selected modulus product'
            hmodulus hcall
          rw [htail, hstep']
          have hselected : selected =
              activeLifted[candidate[index].toInt64.toNat]! := by
            simp [selected, hactive]
          rw [hselected]
          simp [SelectedProductMod, hsuffix, mul_assoc]
      next hindex =>
        have hle : candidate.size ≤ index := Nat.le_of_not_gt hindex
        have hout := Except.ok.inj hrun
        subst output
        simp [SelectedProductMod, List.drop_eq_nil_iff.mpr hle]

private noncomputable def intTermsToPoly (terms : List (UMonomial × Int)) :
    Polynomial Int :=
  (terms.map fun term => Polynomial.monomial term.1.deg term.2).sum

theorem multiplyRowLoop_toPoly (left : UMonomial × Int)
    (right : SparsePolyZZ) (rightIndex : Nat) (terms : SparsePolyZZ) :
    intTermsToPoly (Generated.StrictRecombine.multiplyRowLoop
      left right rightIndex terms).toList =
      intTermsToPoly terms.toList +
        Polynomial.monomial left.1.deg left.2 *
          intTermsToPoly (right.toList.drop rightIndex) := by
  induction hmeasure : right.size - rightIndex using Nat.strong_induction_on
      generalizing rightIndex terms with
  | h measure ih =>
      rw [Generated.StrictRecombine.multiplyRowLoop]
      split
      next hright =>
        rw [ih (right.size - (rightIndex + 1)) (by omega)
          (rightIndex + 1) (terms.push
            (⟨left.1.deg + right[rightIndex].1.deg⟩,
              left.2 * right[rightIndex].2)) rfl]
        have hsuffix : right.toList.drop rightIndex = right[rightIndex] ::
            right.toList.drop (rightIndex + 1) := by
          simpa using List.drop_eq_getElem_cons
            (l := right.toList) (i := rightIndex) (by simpa using hright)
        simp [intTermsToPoly, hsuffix, add_comm, add_left_comm,
          mul_add]
        exact (Polynomial.monomial_mul_monomial _ _ _ _).symm
      next hright =>
        have hle : right.size ≤ rightIndex := Nat.le_of_not_gt hright
        simp [intTermsToPoly, List.drop_eq_nil_iff.mpr hle]

theorem multiplyTermsLoop_toPoly (left right : SparsePolyZZ)
    (leftIndex : Nat) (terms : SparsePolyZZ) :
    intTermsToPoly (Generated.StrictRecombine.multiplyTermsLoop
      left right leftIndex terms).toList =
      intTermsToPoly terms.toList +
        intTermsToPoly (left.toList.drop leftIndex) *
          intTermsToPoly right.toList := by
  induction hmeasure : left.size - leftIndex using Nat.strong_induction_on
      generalizing leftIndex terms with
  | h measure ih =>
      rw [Generated.StrictRecombine.multiplyTermsLoop]
      split
      next hleft =>
        rw [ih (left.size - (leftIndex + 1)) (by omega)
          (leftIndex + 1)
          (Generated.StrictRecombine.multiplyRowLoop left[leftIndex]
            right 0 terms) rfl]
        rw [multiplyRowLoop_toPoly]
        have hsuffix : left.toList.drop leftIndex = left[leftIndex] ::
            left.toList.drop (leftIndex + 1) := by
          simpa using List.drop_eq_getElem_cons
            (l := left.toList) (i := leftIndex) (by simpa using hleft)
        simp [intTermsToPoly, hsuffix]
        ring
      next hleft =>
        have hle : left.size ≤ leftIndex := Nat.le_of_not_gt hleft
        simp [intTermsToPoly, List.drop_eq_nil_iff.mpr hle]

private theorem intTermsToPoly_modify_add (terms : List (UMonomial × Int))
    (index : Nat) (term : UMonomial × Int) (hindex : index < terms.length)
    (hdegree : terms[index].1.deg = term.1.deg) :
    intTermsToPoly (terms.modify index
      (fun existing => (existing.1, existing.2 + term.2))) =
      intTermsToPoly terms + Polynomial.monomial term.1.deg term.2 := by
  induction terms generalizing index with
  | nil => simp at hindex
  | cons head tail ih =>
      cases index with
      | zero =>
          simp [intTermsToPoly] at hdegree ⊢
          rw [hdegree]
          abel
      | succ index =>
          have htail : index < tail.length := by simpa using hindex
          have hdegree' : tail[index].1.deg = term.1.deg := by
            simpa using hdegree
          simp only [List.modify_succ_cons, intTermsToPoly, List.map_cons,
            List.sum_cons]
          change Polynomial.monomial head.1.deg head.2 +
              intTermsToPoly (tail.modify index
                (fun existing => (existing.1, existing.2 + term.2))) = _
          rw [ih index htail hdegree']
          change _ + (intTermsToPoly tail + _) =
            _ + intTermsToPoly tail + _
          abel

private theorem groupTermsStep_toPoly (acc : SparsePolyZZ)
    (term : UMonomial × Int) :
    intTermsToPoly
      ((match acc.findIdx? (fun t : UMonomial × Int =>
          t.fst.deg = term.fst.deg) with
        | some index => acc.modify index
            (fun existing => (existing.1, existing.2 + term.2))
        | none => acc.push term).toList) =
      intTermsToPoly acc.toList + Polynomial.monomial term.1.deg term.2 := by
  split
  next index hfind =>
    obtain ⟨hindex, hdegree, _⟩ :=
      Array.findIdx?_eq_some_iff_getElem.mp hfind
    simpa using intTermsToPoly_modify_add acc.toList index term
      (by simpa using hindex) (by simpa using hdegree)
  next hfind =>
    simp [intTermsToPoly]

private theorem groupTerms_toPoly_aux (source : List (UMonomial × Int))
    (acc : SparsePolyZZ) :
    intTermsToPoly
      (source.foldl (fun acc term =>
        match acc.findIdx? (fun t => t.fst.deg = term.fst.deg) with
        | some index => acc.modify index
            (fun existing => (existing.1, existing.2 + term.2))
        | none => acc.push term) acc).toList =
      intTermsToPoly acc.toList + intTermsToPoly source := by
  induction source generalizing acc with
  | nil => simp [intTermsToPoly]
  | cons head tail ih =>
      simp only [List.foldl_cons]
      rw [ih]
      rw [groupTermsStep_toPoly]
      simp [intTermsToPoly]
      abel

private theorem groupTerms_toPoly (terms : SparsePolyZZ) :
    intTermsToPoly
      (terms.foldl (fun acc term =>
        match acc.findIdx? (fun t => t.fst.deg = term.fst.deg) with
        | some index => acc.modify index
            (fun existing => (existing.1, existing.2 + term.2))
        | none => acc.push term) #[]).toList =
      intTermsToPoly terms.toList := by
  rw [← Array.foldl_toList]
  simpa [intTermsToPoly] using groupTerms_toPoly_aux terms.toList #[]

private theorem filterNonzero_toPoly (terms : SparsePolyZZ) :
    intTermsToPoly (terms.filter (fun term => term.2 ≠ 0)).toList =
      intTermsToPoly terms.toList := by
  rw [Array.toList_filter]
  induction terms.toList with
  | nil => simp [intTermsToPoly]
  | cons head tail ih =>
      rcases head with ⟨monomial, coefficient⟩
      by_cases hzero : coefficient = 0
      · simpa [hzero, intTermsToPoly] using ih
      · simpa [hzero, intTermsToPoly] using ih

private theorem mergeSort_toPoly (terms : List (UMonomial × Int)) :
    intTermsToPoly
        (terms.mergeSort (fun a b => a.fst.deg > b.fst.deg)) =
      intTermsToPoly terms := by
  exact ((List.mergeSort_perm terms _).map
    (fun term => Polynomial.monomial term.1.deg term.2)).sum_eq

theorem normalization_toPoly (terms : SparsePolyZZ) :
    SparsePolyZZ.toPoly (SparsePolyZZ.normalization terms) =
      SparsePolyZZ.toPoly terms := by
  unfold SparsePolyZZ.normalization SparsePolyZZ.toPoly
  rw [List.toList_toArray]
  change intTermsToPoly
      (((terms.foldl (fun acc term =>
          match acc.findIdx? (fun t => t.fst.deg = term.fst.deg) with
          | some index => acc.modify index
              (fun existing => (existing.1, existing.2 + term.2))
          | none => acc.push term) #[]).filter
        (fun term => term.2 ≠ 0)).toList.mergeSort
          (fun a b => a.fst.deg > b.fst.deg)) = intTermsToPoly terms.toList
  rw [mergeSort_toPoly, filterNonzero_toPoly, groupTerms_toPoly]

theorem multiplyNormalizeRaw_toPoly (left right output : SparsePolyZZ)
    (hrun : Generated.StrictRecombine.multiplyNormalizeRaw left right =
      .ok output) :
    SparsePolyZZ.toPoly output =
      SparsePolyZZ.toPoly left * SparsePolyZZ.toPoly right := by
  unfold Generated.StrictRecombine.multiplyNormalizeRaw at hrun
  have houtput := Except.ok.inj hrun
  subst output
  rw [normalization_toPoly]
  change intTermsToPoly
      (Generated.StrictRecombine.multiplyTermsLoop left right 0 #[]).toList = _
  rw [multiplyTermsLoop_toPoly]
  simp [intTermsToPoly, SparsePolyZZ.toPoly]

end Refinement.StrictRecombine
