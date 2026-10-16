import CLPoly.Refinement.VanHoeijEntryExecution
import CLPoly.Refinement.VanHoeijValidation

set_option autoImplicit false

namespace Refinement.StrictRecombine
open Generated.StrictRecombine Generated.StrictVanHoeij
open Generated.StrictVanHoeijParameters Generated.StrictVanHoeijController
open Generated.StrictVanHoeijEntry

theorem tryRecombine_complete_trace (execution : EntryExecution) (f : SparsePolyZZ)
    (lifted : Array SparsePolyZZ) (modulus : ZZ) (countFits : lifted.size < 2 ^ 31)
    (after : VanHoeijResult) (hcomplete : after.complete = true)
    (hrun : tryRecombine execution f lifted modulus countFits = .ok after) :
    (2 ≤ Int32.ofNat lifted.size ∧ 1 < modulus) ∧
    ∃ height singletons,
      cldBound f = .ok height ∧
      singletonGroups (Int32.ofNat lifted.size) = .ok singletons ∧
      (validatePartition f lifted modulus singletons {} = .ok (true, after) ∨
      ∃ parameters original truncated before,
        selectParameters (Int32.ofNat lifted.size) height modulus = some parameters ∧
        cldPolys execution.cld f lifted modulus = .ok original ∧
        truncateCld original parameters.divisor = .ok truncated ∧
        before.complete = false ∧
        ∃ countSize : (Int32.ofNat lifted.size).toInt.toNat = truncated.size,
        ∃ modulusNonzero : parameters.modulus ≠ 0,
        columnLoop execution.controller f lifted truncated modulus parameters.modulus 32
          parameters.residual (Int32.ofNat lifted.size) (cldWidth truncated)
          (min (2 : Int32) (cldWidth truncated)) countSize (by decide) modulusNonzero
          (execution.initialBound (cldWidth truncated)) before = .ok after) := by
  unfold tryRecombine at hrun
  dsimp only at hrun
  split at hrun
  next hguard =>
    refine ⟨hguard, ?_⟩
    cases hh : cldBound f with
    | error fault => simp only [hh, bind, Except.bind] at hrun; contradiction
    | ok height =>
      simp only [hh, bind, Except.bind] at hrun
      cases hs : singletonGroups (Int32.ofNat lifted.size) with
      | error fault => simp only [hs, bind, Except.bind] at hrun; contradiction
      | ok singletons =>
        simp only [hs, bind, Except.bind] at hrun
        refine ⟨height, singletons, rfl, rfl, ?_⟩
        cases hv : validatePartition f lifted modulus singletons {} with
        | error fault => simp only [hv, bind, Except.bind] at hrun; contradiction
        | ok checked =>
          rcases checked with ⟨accepted, before⟩
          simp only [hv, bind, Except.bind] at hrun
          cases accepted with
          | true =>
            simp only [↓reduceIte] at hrun
            have he : before = after := Except.ok.inj hrun
            subst before
            exact Or.inl rfl
          | false =>
            have hbefore : before.complete = false := by
              have he := validatePartition_failure_preserves f lifted modulus singletons {} before hv
              have hc := congrArg VanHoeijResult.complete he
              exact hc
            simp only [Bool.false_eq_true, ↓reduceIte] at hrun
            split at hrun
            next hp =>
              have he : before = after := Except.ok.inj hrun
              subst before
              simp [hbefore] at hcomplete
            next parameters hp =>
              split at hrun
              next fault hc => contradiction
              next original hc =>
                split at hrun
                next fault ht => contradiction
                next truncated ht =>
                  exact Or.inr ⟨parameters, original, truncated, before, hp, hc, ht,
                    hbefore, _, _, hrun⟩
  next hguard => contradiction

end Refinement.StrictRecombine
