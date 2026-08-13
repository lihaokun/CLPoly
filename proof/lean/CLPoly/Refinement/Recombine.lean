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

end Refinement.StrictRecombine
