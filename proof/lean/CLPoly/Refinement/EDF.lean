/-
  Strict EDF refinement is intentionally not exported yet.

  The former contents of this module selected an L2 existence witness when a
  bounded generated run failed or when its result did not satisfy
  `EDFCorrect`.  That construction was a certified specification fallback, not
  a refinement of the C++ `__edf_Zp` execution, and has therefore been removed.

  A future theorem in this module must start from a cpp2lean-generated EDF trace
  (including RNG state transitions), prove every successful split and recursive
  output invariant, and expose failure/nontermination explicitly.  It must not
  manufacture an L2 result on an unproved execution branch.
-/
import CLPoly.Algorithm.EDF
import CLPoly.Generated.StrictEDF
import CLPoly.Refinement.Basic

set_option autoImplicit false

open Polynomial
open CLPoly.Math

namespace Refinement

namespace StrictEDF

/-- Decode the concrete C++ EDF accumulator without manufacturing or replacing
any factor. -/
noncomputable def edfResultToL2 (p : Nat) (result : Array SparsePolyZp) :
    List (Polynomial (ZMod p)) :=
  result.toList.map (SparsePolyZp.toPoly p)

@[simp] theorem edfResultToL2_empty (p : Nat) :
    edfResultToL2 p (#[] : Array SparsePolyZp) = [] := by
  simp [edfResultToL2]

@[simp] theorem edfResultToL2_push (p : Nat) (result : Array SparsePolyZp)
    (factor : SparsePolyZp) :
    edfResultToL2 p (result.push factor) =
      edfResultToL2 p result ++ [SparsePolyZp.toPoly p factor] := by
  simp [edfResultToL2]

private theorem certifyRawExec_ok {α : Type} (run : RawExec α) (output : α)
    (hrun : run = .ok output) :
    ∃ certified, Generated.StrictEDF.certifyRawExec run = .ok certified ∧
      certified.val = output := by
  cases run with
  | error fault => simp at hrun
  | ok value =>
      cases hrun
      exact ⟨⟨output, rfl⟩, rfl, rfl⟩

/-- Exact execution of the generated C++ base branch.  The theorem refers to
the generated state machine itself and to the concrete `makeMonic` run; it is
not an L2 execution fallback. -/
theorem rawState_base_run
    (ops : Generated.StrictEDF.EDFRawOps)
    (law : Generated.StrictEDF.EDFRetryLaw ops)
    (state : Generated.StrictEDF.EDFState ops)
    (hdegree : ((get_deg state.f).toUInt64 == state.d) = true)
    (hmonic : ops.makeMonic state.f = .ok state.f) :
    Generated.StrictEDF.__edf_Zp_raw_ir_state ops law state =
      .ok (state.result.push state.f, state.rng) := by
  rw [Generated.StrictEDF.__edf_Zp_raw_ir_state.eq_1]
  simp only [hdegree, ↓reduceIte]
  rcases certifyRawExec_ok _ _ hmonic with ⟨certified, hcertified, hvalue⟩
  rw [hcertified]
  simp [hvalue]

/-- L2 semantic fact used for the same C++ base branch: under the real EDF
preconditions, the factor appended by that branch is irreducible of degree
`d`. -/
theorem base_factor_correct
    {p : Nat} [Fact (Nat.Prime p)]
    (f : Polynomial (ZMod p)) (d : Nat)
    (hmonic : f.Monic) (hd : 0 < d) (hdegree : f.natDegree = d)
    (hfactors : ∀ q, Irreducible q → q ∣ f → q.natDegree = d) :
    EDFCorrect f d [f] := by
  have hirreducible := edf_base_irred f d hmonic hd hdegree hfactors
  exact ⟨by simp, by simp [hirreducible, hmonic, hdegree]⟩

end StrictEDF

-- The public L1→L2 EDF theorem remains deliberately absent until the retry
-- and recursive split simulations below this base layer are closed.

end Refinement
