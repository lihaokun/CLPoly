import CLPoly.Refinement.VanHoeijControllerExecution

set_option autoImplicit false

namespace Refinement.StrictRecombine
open Generated.StrictRecombine Generated.StrictVanHoeij
open Generated.StrictVanHoeijController

theorem columnRound_trace (execution : ControllerExecution)
    (f : SparsePolyZZ) (lifted cld : Array SparsePolyZZ)
    (originalModulus modulus scale residual : ZZ) (count target : Int32)
    (countSize : count.toInt.toNat = cld.size)
    (scaleNonzero : scale ≠ 0) (modulusNonzero : modulus ≠ 0)
    (before after : VanHoeijResult) (accepted : Bool)
    (hrun : columnRound execution f lifted cld originalModulus modulus scale residual
      count target countSize scaleNonzero modulusNonzero before = .ok (accepted, after)) :
    ∃ initial built,
      makeInitialMatrix count.toInt.toNat scale = .ok initial ∧
      buildCldMatrix initial cld 0 target.toInt.toNat modulus = .ok built ∧
      ∃ ready : execution.lll.inputValid built.1,
      let added := Int32.ofNat built.2
      let bound := count.toInt * scale * scale + added.toInt * residual * residual
      ∃ reduced candidates,
        lllReduceWithNorms execution.lll built.1 ready bound = .ok reduced ∧
        extractCandidates (shortSpanFromNorms reduced.2.2.2 bound)
          reduced.2.1 count.toInt.toNat = .ok candidates ∧
        validatePartition f lifted originalModulus candidates
          { before with cldColumns := added, retainedRows := Int32.ofNat (shortSpanFromNorms reduced.2.2.2 bound).size } =
          .ok (accepted, after) := by
  unfold columnRound at hrun
  split at hrun
  next fault hi => contradiction
  next initial hi =>
    split at hrun
    next fault hb => contradiction
    next built hb =>
      let ready := execution.buildValid cld target.toInt.toNat scale modulus initial built
        scaleNonzero modulusNonzero (countSize ▸ hi) hb
      let bound := count.toInt * scale * scale +
        (Int32.ofNat built.2).toInt * residual * residual
      cases hl : lllReduceWithNorms execution.lll built.1 ready bound with
      | error fault =>
        dsimp only [bound] at hl
        simp only [hl, bind, Except.bind] at hrun
        contradiction
      | ok reduced =>
        dsimp only [bound] at hl
        cases hc : extractCandidates (shortSpanFromNorms reduced.2.2.2 bound)
            reduced.2.1 count.toInt.toNat with
        | error fault =>
          dsimp only [bound] at hc
          simp only [hl, hc, bind, Except.bind] at hrun
          contradiction
        | ok candidates =>
          refine ⟨initial, built, hi, hb, ready, reduced, candidates, hl, hc, ?_⟩
          dsimp only [bound] at hc
          simpa only [hl, hc, bind, Except.bind] using hrun

end Refinement.StrictRecombine
