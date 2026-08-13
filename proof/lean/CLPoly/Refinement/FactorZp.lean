/-
  Genuine refinement of the generated C++ `__factor_Zp` entry.

  This file reasons about the strict, source-shaped L1 loops in
  `Generated.StrictFactorZp`.  Component executions are instantiated with the
  already verified strict SQF, DDF and EDF entries; no L2 algorithm is used as
  an executable fallback.
-/
import CLPoly.Generated.StrictFactorZp
import CLPoly.Refinement.Generated
import CLPoly.Pipeline.FactorZp

set_option autoImplicit false

open Polynomial
open CLPoly.Math

namespace Refinement.StrictFactorZp

/-- Decode the concrete output array of C++ `__factor_Zp`. -/
noncomputable def factorResultToL2 (p : Nat)
    (result : Array (SparsePolyZp × UInt64)) :
    List (Polynomial (ZMod p) × Nat) :=
  result.toList.map fun item =>
    (SparsePolyZp.toPoly p item.1, item.2.toNat)

/-- Decode an EDF factor array after the C++ innermost loop has attached one
SQF multiplicity to every element. -/
noncomputable def attachMultiplicityToL2 (p : Nat)
    (factors : Array SparsePolyZp) (multiplicity : UInt64) :
    List (Polynomial (ZMod p) × Nat) :=
  factors.toList.map fun factor =>
    (SparsePolyZp.toPoly p factor, multiplicity.toNat)

theorem factorResultToL2_push (p : Nat)
    (result : Array (SparsePolyZp × UInt64))
    (factor : SparsePolyZp) (multiplicity : UInt64) :
    factorResultToL2 p (result.push (factor, multiplicity)) =
      factorResultToL2 p result ++
        [(SparsePolyZp.toPoly p factor, multiplicity.toNat)] := by
  simp [factorResultToL2]

/-- Exact semantic account of the innermost generated range-for.  Its output
is the old accumulator followed by precisely the unvisited suffix of the EDF
array, with the source SQF multiplicity attached. -/
theorem factorLoop0_toL2 (p : Nat) (factors : Array SparsePolyZp)
    (multiplicity : UInt64) (index : Nat)
    (result : Array (SparsePolyZp × UInt64)) :
    factorResultToL2 p
        (Generated.StrictFactorZp._loop___factor_Zp_0_raw_ir
          factors multiplicity index result) =
      factorResultToL2 p result ++
        (factors.toList.drop index).map fun factor =>
          (SparsePolyZp.toPoly p factor, multiplicity.toNat) := by
  induction hmeasure : factors.size - index using Nat.strong_induction_on
    generalizing index result with
  | h measure ih =>
      rw [Generated.StrictFactorZp._loop___factor_Zp_0_raw_ir]
      split
      next hindex =>
        rw [ih (factors.size - (index + 1)) (by omega)
          (index + 1) (result.push (factors[index], multiplicity)) rfl]
        rw [factorResultToL2_push, List.append_assoc]
        congr 1
        have hsuffix : factors.toList.drop index =
            factors[index] :: factors.toList.drop (index + 1) := by
          simpa using (List.drop_eq_getElem_cons
            (l := factors.toList) (i := index) (by simpa using hindex))
        rw [hsuffix]
        rfl
      next hindex =>
        have hle : factors.size ≤ index := Nat.le_of_not_gt hindex
        have hdrop : factors.toList.drop index = [] := by
          apply List.drop_eq_nil_iff.mpr
          simpa using hle
        simp [hdrop]

theorem factorLoop0_toL2_zero (p : Nat) (factors : Array SparsePolyZp)
    (multiplicity : UInt64)
    (result : Array (SparsePolyZp × UInt64)) :
    factorResultToL2 p
        (Generated.StrictFactorZp._loop___factor_Zp_0_raw_ir
          factors multiplicity 0 result) =
      factorResultToL2 p result ++
        attachMultiplicityToL2 p factors multiplicity := by
  simpa [attachMultiplicityToL2] using
    factorLoop0_toL2 p factors multiplicity 0 result

/-- Structural execution theorem for the generated DDF-component loop.  The
only way to advance is to supply the successful equation of the concrete EDF
call at the current array position.  `Preserved` may express the accumulated
L2 product/irreducibility invariant, but it cannot influence execution. -/
theorem factorLoop1_preserves {State : Type}
    (ops : Generated.StrictFactorZp.FactorZpRawOps State)
    (components : Array (SparsePolyZp × UInt64))
    (Preserved : Array (SparsePolyZp × UInt64) → Prop)
    (hstep : ∀ (index : Nat) (hindex : index < components.size)
        (rng : State) (result : Array (SparsePolyZp × UInt64)),
      Preserved result →
      ∃ factors rng',
        ops.edf #[] components[index].1 components[index].2 rng =
          .ok (factors, rng') ∧
        Preserved
          (Generated.StrictFactorZp._loop___factor_Zp_0_raw_ir
            factors components[index].2 0 result))
    (index : Nat) (rng : State)
    (result : Array (SparsePolyZp × UInt64)) (hresult : Preserved result) :
    ∃ output rng',
      Generated.StrictFactorZp._loop___factor_Zp_1_raw_ir
          ops components index rng result = .ok (output, rng') ∧
      Preserved output := by
  induction hmeasure : components.size - index using Nat.strong_induction_on
    generalizing index rng result with
  | h measure ih =>
      rw [Generated.StrictFactorZp._loop___factor_Zp_1_raw_ir]
      split
      next hindex =>
        rcases hstep index hindex rng result hresult with
          ⟨factors, rng', hedf, hnext⟩
        simp only [hedf]
        exact ih (components.size - (index + 1)) (by omega)
          (index + 1) rng'
          (Generated.StrictFactorZp._loop___factor_Zp_0_raw_ir
            factors components[index].2 0 result) hnext rfl
      next hindex =>
        exact ⟨result, rng, rfl, hresult⟩

/-- Structural execution theorem for the generated SQF-component loop.  Each
step must exhibit the actual DDF run and the actual nested EDF-loop run. -/
theorem factorLoop2_preserves {State : Type}
    (ops : Generated.StrictFactorZp.FactorZpRawOps State)
    (squarefreeFactors : Array (SparsePolyZp × UInt64))
    (Preserved : Array (SparsePolyZp × UInt64) → Prop)
    (hstep : ∀ (index : Nat) (hindex : index < squarefreeFactors.size)
        (rng : State) (result : Array (SparsePolyZp × UInt64)),
      Preserved result →
      ∃ components output rng',
        ops.ddf squarefreeFactors[index].1 = .ok components ∧
        Generated.StrictFactorZp._loop___factor_Zp_1_raw_ir
          ops components 0 rng result = .ok (output, rng') ∧
        Preserved output)
    (index : Nat) (rng : State)
    (result : Array (SparsePolyZp × UInt64)) (hresult : Preserved result) :
    ∃ output,
      Generated.StrictFactorZp._loop___factor_Zp_2_raw_ir
          ops squarefreeFactors index rng result = .ok output ∧
      Preserved output := by
  induction hmeasure : squarefreeFactors.size - index using
    Nat.strong_induction_on generalizing index rng result with
  | h measure ih =>
      rw [Generated.StrictFactorZp._loop___factor_Zp_2_raw_ir]
      split
      next hindex =>
        rcases hstep index hindex rng result hresult with
          ⟨components, output, rng', hddf, hedf, houtput⟩
        simp only [hddf, hedf]
        exact ih (squarefreeFactors.size - (index + 1)) (by omega)
          (index + 1) rng' output houtput rfl
      next hindex =>
        exact ⟨result, rfl, hresult⟩

end Refinement.StrictFactorZp
