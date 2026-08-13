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

private theorem monicListProd {p : Nat} [Fact (Nat.Prime p)]
    (factors : List (Polynomial (ZMod p)))
    (hmonic : ∀ factor ∈ factors, factor.Monic) : factors.prod.Monic := by
  induction factors with
  | nil => exact monic_one
  | cons head tail ih =>
      rw [List.prod_cons]
      exact (hmonic head (by simp)).mul
        (ih fun factor hfactor => hmonic factor (by simp [hfactor]))

private theorem listProdPow {p : Nat} [Fact (Nat.Prime p)]
    (factors : List (Polynomial (ZMod p))) (multiplicity : Nat) :
    factors.prod ^ multiplicity =
      (factors.map fun factor => factor ^ multiplicity).prod := by
  induction factors with
  | nil => simp
  | cons head tail ih => simp [mul_pow, ih]

/-- Representation facts established for every concrete SQF output before it
is passed to the source DDF call. -/
structure DDFReady (this : DenseUPolyZp) (factor : SparsePolyZp) : Prop where
  primeMatches : factor[0]!.2.prime = this._p
  canonical : SparsePolyZp.Canonical this._p.toNat factor
  degreeBound : ∀ term ∈ factor.toList, term.1.deg < 2 ^ 62
  monic : (SparsePolyZp.toPoly this._p.toNat factor).Monic
  squarefree : Squarefree (SparsePolyZp.toPoly this._p.toNat factor)

/-- Erased-proof adapter for the C++ DDF call site.  On the verified pipeline
path this reduces to the exact strict DDF entry; the negative branch only
models an assertion failure for an invalid internal state. -/
noncomputable def strictDDFCall (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (providers : StrictDDF.DDFRawProviders this) (f : SparsePolyZp) :
    RawExec (Array (SparsePolyZp × UInt64)) := by
  classical
  exact if hready : DDFReady this f then
    Generated.StrictDDF.__ddf_Zp_raw_ir
      (StrictDDF.strictDDFRawOps this providers
        (SparsePolyZp.toPoly this._p.toNat f)) f
      (fun _ => StrictDDF.DDFLoopInvariant.initial this f
        f[0]!.2.prime hready.primeMatches hready.canonical
        hready.degreeBound hready.monic hready.squarefree)
  else
    .error .assertionFailure

/-- Erased-proof adapter for the C++ EDF call site.  The executable state and
RNG transition are precisely those of the strict generated EDF entry. -/
noncomputable def strictEDFCall {State : Type}
    (engine : Generated.StrictEDF.RandomEngine State)
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (providers : StrictDDF.DDFRawProviders this)
    (termination : Generated.StrictEDF.EDFTermination
      (StrictEDF.strictEDFRawOps engine this providers))
    (result : Array SparsePolyZp) (f : SparsePolyZp) (d : UInt64)
    (rng : State) : RawExec (Array SparsePolyZp × State) := by
  classical
  exact if hinvariant : StrictEDF.EDFEntryInvariant this f d then
    Generated.StrictEDF.__edf_Zp_raw_ir
      (StrictEDF.strictEDFRawOps engine this providers)
      (StrictEDF.strictEDFSplitLaw engine this providers) termination
      result f d rng hinvariant
  else
    .error .assertionFailure

theorem strictDDFCall_eq (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (providers : StrictDDF.DDFRawProviders this) (f : SparsePolyZp)
    (hready : DDFReady this f) :
    strictDDFCall this providers f =
      Generated.StrictDDF.__ddf_Zp_raw_ir
        (StrictDDF.strictDDFRawOps this providers
          (SparsePolyZp.toPoly this._p.toNat f)) f
        (fun _ => StrictDDF.DDFLoopInvariant.initial this f
          f[0]!.2.prime hready.primeMatches hready.canonical
          hready.degreeBound hready.monic hready.squarefree) := by
  classical
  simp [strictDDFCall, hready]

theorem strictEDFCall_eq {State : Type}
    (engine : Generated.StrictEDF.RandomEngine State)
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (providers : StrictDDF.DDFRawProviders this)
    (termination : Generated.StrictEDF.EDFTermination
      (StrictEDF.strictEDFRawOps engine this providers))
    (result : Array SparsePolyZp) (f : SparsePolyZp) (d : UInt64)
    (rng : State) (hinvariant : StrictEDF.EDFEntryInvariant this f d) :
    strictEDFCall engine this providers termination result f d rng =
      Generated.StrictEDF.__edf_Zp_raw_ir
        (StrictEDF.strictEDFRawOps engine this providers)
        (StrictEDF.strictEDFSplitLaw engine this providers) termination
        result f d rng hinvariant := by
  classical
  simp [strictEDFCall, hinvariant]

/-- The proof-erased DDF call executes successfully on every SQF-produced
component and equips each concrete output element for the following EDF call. -/
theorem strictDDFCall_refines
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (providers : StrictDDF.DDFRawProviders this) (f : SparsePolyZp)
    (hready : DDFReady this f) :
    ∃ output,
      strictDDFCall this providers f = .ok output ∧
      StrictDDF.ddfResultToL2 this._p.toNat output =
        ddf (SparsePolyZp.toPoly this._p.toNat f) ∧
      ∀ item ∈ output.toList,
        StrictEDF.EDFEntryInvariant this item.1 item.2 := by
  rcases Refinement.__ddf_Zp_raw_ir_refines_ddf this providers f
      hready.primeMatches hready.canonical hready.degreeBound hready.monic
      hready.squarefree with
    ⟨output, hrun, hsemantic, hcanonical, hdegreePositive,
      hindexPositive⟩
  refine ⟨output, ?_, hsemantic, ?_⟩
  · rw [strictDDFCall_eq this providers f hready]
    exact hrun
  · intro item hitem
    have hdecoded :
        (SparsePolyZp.toPoly this._p.toNat item.1, item.2.toNat) ∈
          ddf (SparsePolyZp.toPoly this._p.toNat f) := by
      rw [← hsemantic]
      exact List.mem_map.mpr ⟨item, hitem, rfl⟩
    have hcorrect := ddf_correct
      (SparsePolyZp.toPoly this._p.toNat f) hready.monic hready.squarefree
    have hcomponentDvd := hcorrect.1 _ hdecoded
    have hcomponentMonic := hcorrect.2.2.2.2 _ hdecoded
    have hcomponentCanonical := hcanonical item hitem
    have hcomponentPositive := hdegreePositive item hitem
    have hcomponentNonempty :=
      StrictSquarefreeZp.sparsePolyZp_size_pos_of_toPoly_ne_zero
        this._p.toNat item.1 hcomponentMonic.ne_zero
    have hheadMem : item.1[0] ∈ item.1.toList := by
      simpa using Array.getElem_mem item.1 0 hcomponentNonempty
    have hhead : item.1[0]! = item.1[0] := by
      rw [getElem!_def, getElem?_def]
      simp [hcomponentNonempty]
    have hprimeNat : item.1[0]!.2.prime.toNat = this._p.toNat := by
      rw [hhead]
      exact (hcomponentCanonical.1 item.1[0] hheadMem).1
    have hfDegree := StrictDDF.canonical_natDegree_lt_of_terms_lt
      this._p.toNat f hready.canonical hready.monic.ne_zero (2 ^ 62)
      hready.degreeBound
    have hcomponentDegree :
        (SparsePolyZp.toPoly this._p.toNat item.1).natDegree < 2 ^ 62 :=
      lt_of_le_of_lt
        (Polynomial.natDegree_le_of_dvd hcomponentDvd hready.monic.ne_zero)
        hfDegree
    refine ⟨hcomponentCanonical, ?_, ?_, hcomponentMonic,
      hcomponentPositive, hindexPositive item hitem,
      hready.squarefree.squarefree_of_dvd hcomponentDvd, ?_⟩
    · intro _
      exact UInt64.toNat_inj.mp hprimeNat
    · intro term hterm
      exact lt_of_le_of_lt
        (StrictDDF.canonical_term_degree_le_natDegree this._p.toNat item.1
          hcomponentCanonical term hterm) hcomponentDegree
    · intro q hq hqDvd
      exact hcorrect.2.1 _ hdecoded q hq hqDvd

/-- The proof-erased EDF call is exactly the strict generated execution and
returns a concrete `EDFCorrect` factor list. -/
theorem strictEDFCall_refines {State : Type}
    (engine : Generated.StrictEDF.RandomEngine State)
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (providers : StrictDDF.DDFRawProviders this)
    (termination : Generated.StrictEDF.EDFTermination
      (StrictEDF.strictEDFRawOps engine this providers))
    (f : SparsePolyZp) (d : UInt64) (rng : State)
    (hinvariant : StrictEDF.EDFEntryInvariant this f d) :
    ∃ (output : Array SparsePolyZp) (rng' : State)
        (factors : List (Polynomial (ZMod this._p.toNat))),
      strictEDFCall engine this providers termination #[] f d rng =
        .ok (output, rng') ∧
      StrictEDF.edfResultToL2 this._p.toNat output = factors ∧
      EDFCorrect (SparsePolyZp.toPoly this._p.toNat f) d.toNat factors := by
  rcases Refinement.__edf_Zp_raw_ir_refines_edf engine this providers
      termination #[] f d rng hinvariant with
    ⟨output, rng', factors, hrun, hdecode, hcorrect⟩
  refine ⟨output, rng', factors, ?_, ?_, hcorrect⟩
  · rw [strictEDFCall_eq engine this providers termination #[] f d rng
      hinvariant]
    exact hrun
  · simpa using hdecode

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
    (multiplicity : UInt64)
    (Preserved : Array (SparsePolyZp × UInt64) → Prop)
    (hstep : ∀ (index : Nat) (hindex : index < components.size)
        (rng : State) (result : Array (SparsePolyZp × UInt64)),
      Preserved result →
      ∃ factors rng',
        ops.edf #[] components[index].1 components[index].2 rng =
          .ok (factors, rng') ∧
        Preserved
          (Generated.StrictFactorZp._loop___factor_Zp_0_raw_ir
            factors multiplicity 0 result))
    (index : Nat) (rng : State)
    (result : Array (SparsePolyZp × UInt64)) (hresult : Preserved result) :
    ∃ output rng',
      Generated.StrictFactorZp._loop___factor_Zp_1_raw_ir
          ops components multiplicity index rng result = .ok (output, rng') ∧
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
            factors multiplicity 0 result) hnext rfl
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
          ops components squarefreeFactors[index].2 0 rng result =
            .ok (output, rng') ∧
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

/-- Product of the concrete polynomial components in an unvisited DDF suffix. -/
noncomputable def componentSuffixProduct (p : Nat)
    (components : Array (SparsePolyZp × UInt64)) (index : Nat) :
    Polynomial (ZMod p) :=
  ((components.toList.drop index).map fun item =>
    SparsePolyZp.toPoly p item.1).prod

/-- Genuine refinement of the generated middle loop.  Each EDF factor list is
obtained from the actual strict call at the current DDF array element. -/
theorem factorLoop1_refines {State : Type}
    (engine : Generated.StrictEDF.RandomEngine State)
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (providers : StrictDDF.DDFRawProviders this)
    (termination : Generated.StrictEDF.EDFTermination
      (StrictEDF.strictEDFRawOps engine this providers))
    (components : Array (SparsePolyZp × UInt64))
    (hready : ∀ item ∈ components.toList,
      StrictEDF.EDFEntryInvariant this item.1 item.2)
    (multiplicity : UInt64) (index : Nat) (rng : State)
    (result : Array (SparsePolyZp × UInt64)) :
    let ops : Generated.StrictFactorZp.FactorZpRawOps State := {
      makeMonic := fun _ => .error .assertionFailure
      squarefree := fun _ => .error .assertionFailure
      ddf := strictDDFCall this providers
      edf := strictEDFCall engine this providers termination
      sortByDegree := fun _ => .error .assertionFailure }
    ∃ (output : Array (SparsePolyZp × UInt64)) (rng' : State)
        (factors : List (Polynomial (ZMod this._p.toNat))),
      Generated.StrictFactorZp._loop___factor_Zp_1_raw_ir
          ops components multiplicity index rng result = .ok (output, rng') ∧
      factorResultToL2 this._p.toNat output =
        factorResultToL2 this._p.toNat result ++
          (factors.map fun q => (q, multiplicity.toNat)) ∧
      componentSuffixProduct this._p.toNat components index = factors.prod ∧
      ∀ q ∈ factors, Irreducible q ∧ Monic q := by
  dsimp only
  induction hmeasure : components.size - index using Nat.strong_induction_on
    generalizing index rng result with
  | h measure ih =>
      rw [Generated.StrictFactorZp._loop___factor_Zp_1_raw_ir]
      split
      next hindex =>
        let component := components[index]
        have hcomponentMem : component ∈ components.toList := by
          exact List.getElem_mem (by simpa using hindex)
        have hinvariant := hready component hcomponentMem
        rcases strictEDFCall_refines engine this providers termination
            component.1 component.2 rng hinvariant with
          ⟨edfOutput, rngNext, edfFactors, hedfRun, hedfDecode,
            hedfCorrect⟩
        change ∃ (output : Array (SparsePolyZp × UInt64)) (rng' : State)
            (factors : List (Polynomial (ZMod this._p.toNat))),
          (match strictEDFCall engine this providers termination #[]
              component.1 component.2 rng with
            | Except.error fault => Except.error fault
            | Except.ok (factors, rng') =>
              Generated.StrictFactorZp._loop___factor_Zp_1_raw_ir
                { makeMonic := fun _ => .error .assertionFailure
                  squarefree := fun _ => .error .assertionFailure
                  ddf := strictDDFCall this providers
                  edf := strictEDFCall engine this providers termination
                  sortByDegree := fun _ => .error .assertionFailure }
                components multiplicity (index + 1) rng'
                (Generated.StrictFactorZp._loop___factor_Zp_0_raw_ir
                  factors multiplicity 0 result)) = Except.ok (output, rng') ∧
          factorResultToL2 this._p.toNat output =
            factorResultToL2 this._p.toNat result ++
              (factors.map fun q => (q, multiplicity.toNat)) ∧
          componentSuffixProduct this._p.toNat components index =
            factors.prod ∧
          ∀ q ∈ factors, Irreducible q ∧ Monic q
        rw [hedfRun]
        rcases ih (components.size - (index + 1)) (by omega)
            (index + 1) rngNext
            (Generated.StrictFactorZp._loop___factor_Zp_0_raw_ir
              edfOutput multiplicity 0 result) rfl with
          ⟨output, rng', tailFactors, htailRun, htailDecode,
            htailProduct, htailQuality⟩
        refine ⟨output, rng', edfFactors ++ tailFactors, htailRun, ?_, ?_, ?_⟩
        · rw [htailDecode, factorLoop0_toL2_zero]
          simp only [attachMultiplicityToL2, StrictEDF.edfResultToL2,
            List.map_append]
          have hmapped := congrArg
            (List.map fun q : Polynomial (ZMod this._p.toNat) =>
              (q, multiplicity.toNat)) hedfDecode
          simpa [StrictEDF.edfResultToL2, List.map_map,
            Function.comp_def, List.append_assoc] using hmapped
        · have hsuffix : components.toList.drop index =
              component :: components.toList.drop (index + 1) := by
            simpa [component] using (List.drop_eq_getElem_cons
              (l := components.toList) (i := index) (by simpa using hindex))
          rw [componentSuffixProduct, hsuffix, List.map_cons, List.prod_cons,
            List.prod_append, ← htailProduct]
          have hleftMonic := hinvariant.monic
          have hrightMonic := monicListProd edfFactors
            (fun q hq => (hedfCorrect.2 q hq).2.1)
          have hcomponentEq : SparsePolyZp.toPoly this._p.toNat component.1 =
              edfFactors.prod :=
            _root_.eq_of_associated_monic _ _ hleftMonic hrightMonic
              hedfCorrect.1
          rw [hcomponentEq]
          rfl
        · intro q hq
          rcases List.mem_append.mp hq with hq | hq
          · exact ⟨(hedfCorrect.2 q hq).1, (hedfCorrect.2 q hq).2.1⟩
          · exact htailQuality q hq
      next hindex =>
        have hle : components.size ≤ index := Nat.le_of_not_gt hindex
        have hdrop : components.toList.drop index = [] := by
          exact List.drop_eq_nil_iff.mpr (by simpa using hle)
        refine ⟨result, rng, [], rfl, by simp, ?_, by simp⟩
        simp [componentSuffixProduct, hdrop]

/-- A concrete SQF component together with exactly the representation facts
needed by the strict DDF entry. -/
structure SQFOutputReady (this : DenseUPolyZp)
    (item : SparsePolyZp × UInt64) : Prop extends DDFReady this item.1

/-- Product represented by the unvisited suffix of the concrete SQF array. -/
noncomputable def squarefreeSuffixProduct (p : Nat)
    (squarefreeFactors : Array (SparsePolyZp × UInt64)) (index : Nat) :
    Polynomial (ZMod p) :=
  ((squarefreeFactors.toList.drop index).map fun item =>
    SparsePolyZp.toPoly p item.1 ^ item.2.toNat).prod

/-- Genuine refinement of the generated outer SQF loop.  The returned list is
the exact flattening of strict DDF/EDF executions for the concrete unvisited
SQF suffix, with each original SQF multiplicity preserved. -/
theorem factorLoop2_refines {State : Type}
    (engine : Generated.StrictEDF.RandomEngine State)
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (providers : StrictDDF.DDFRawProviders this)
    (termination : Generated.StrictEDF.EDFTermination
      (StrictEDF.strictEDFRawOps engine this providers))
    (squarefreeFactors : Array (SparsePolyZp × UInt64))
    (hready : ∀ item ∈ squarefreeFactors.toList,
      SQFOutputReady this item)
    (index : Nat) (rng : State)
    (result : Array (SparsePolyZp × UInt64)) :
    let ops : Generated.StrictFactorZp.FactorZpRawOps State := {
      makeMonic := fun _ => .error .assertionFailure
      squarefree := fun _ => .error .assertionFailure
      ddf := strictDDFCall this providers
      edf := strictEDFCall engine this providers termination
      sortByDegree := fun _ => .error .assertionFailure }
    ∃ (output : Array (SparsePolyZp × UInt64))
        (factors : List (Polynomial (ZMod this._p.toNat) × Nat)),
      Generated.StrictFactorZp._loop___factor_Zp_2_raw_ir
          ops squarefreeFactors index rng result = .ok output ∧
      factorResultToL2 this._p.toNat output =
        factorResultToL2 this._p.toNat result ++ factors ∧
      squarefreeSuffixProduct this._p.toNat squarefreeFactors index =
        (factors.map fun item => item.1 ^ item.2).prod ∧
      ∀ item ∈ factors, Irreducible item.1 ∧ Monic item.1 := by
  dsimp only
  induction hmeasure : squarefreeFactors.size - index using
    Nat.strong_induction_on generalizing index rng result with
  | h measure ih =>
      rw [Generated.StrictFactorZp._loop___factor_Zp_2_raw_ir]
      split
      next hindex =>
        let sqfItem := squarefreeFactors[index]
        have hsqfMem : sqfItem ∈ squarefreeFactors.toList := by
          exact List.getElem_mem (by simpa using hindex)
        have hsqfReady := hready sqfItem hsqfMem
        rcases strictDDFCall_refines this providers sqfItem.1 hsqfReady.toDDFReady with
          ⟨components, hddfRun, hddfDecode, hedfReady⟩
        change ∃ output factors,
          (match strictDDFCall this providers sqfItem.1 with
            | Except.error fault => Except.error fault
            | Except.ok components =>
              match Generated.StrictFactorZp._loop___factor_Zp_1_raw_ir
                  { makeMonic := fun _ => .error .assertionFailure
                    squarefree := fun _ => .error .assertionFailure
                    ddf := strictDDFCall this providers
                    edf := strictEDFCall engine this providers termination
                    sortByDegree := fun _ => .error .assertionFailure }
                  components sqfItem.2 0 rng result with
              | Except.error fault => Except.error fault
              | Except.ok (result', rng') =>
                Generated.StrictFactorZp._loop___factor_Zp_2_raw_ir
                  { makeMonic := fun _ => .error .assertionFailure
                    squarefree := fun _ => .error .assertionFailure
                    ddf := strictDDFCall this providers
                    edf := strictEDFCall engine this providers termination
                    sortByDegree := fun _ => .error .assertionFailure }
                  squarefreeFactors (index + 1) rng' result') =
              Except.ok output ∧
          factorResultToL2 this._p.toNat output =
            factorResultToL2 this._p.toNat result ++ factors ∧
          squarefreeSuffixProduct this._p.toNat squarefreeFactors index =
            (factors.map fun item => item.1 ^ item.2).prod ∧
          ∀ item ∈ factors, Irreducible item.1 ∧ Monic item.1
        rw [hddfRun]
        rcases factorLoop1_refines engine this providers termination components
            hedfReady sqfItem.2 0 rng result with
          ⟨middleOutput, rngNext, componentFactors, hmiddleRun,
            hmiddleDecode, hmiddleProduct, hmiddleQuality⟩
        simp only
        rw [hmiddleRun]
        rcases ih (squarefreeFactors.size - (index + 1)) (by omega)
            (index + 1) rngNext middleOutput rfl with
          ⟨output, tailFactors, htailRun, htailDecode,
            htailProduct, htailQuality⟩
        let currentFactors := componentFactors.map fun q =>
          (q, sqfItem.2.toNat)
        refine ⟨output, currentFactors ++ tailFactors, htailRun, ?_, ?_, ?_⟩
        · rw [htailDecode, hmiddleDecode, List.append_assoc]
        · have hsuffix : squarefreeFactors.toList.drop index =
              sqfItem :: squarefreeFactors.toList.drop (index + 1) := by
            simpa [sqfItem] using (List.drop_eq_getElem_cons
              (l := squarefreeFactors.toList) (i := index)
              (by simpa using hindex))
          rw [squarefreeSuffixProduct, hsuffix, List.map_cons, List.prod_cons,
            List.map_append, List.prod_append, ← htailProduct]
          have hcomponentProduct :
              SparsePolyZp.toPoly this._p.toNat sqfItem.1 =
                componentFactors.prod := by
            have hddfCorrect := ddf_correct
              (SparsePolyZp.toPoly this._p.toNat sqfItem.1)
              hsqfReady.monic hsqfReady.squarefree
            rw [← hddfDecode] at hddfCorrect
            have hddfAssociated := hddfCorrect.2.2.2.1
            have hddfMonic := monicListProd
              ((StrictDDF.ddfResultToL2 this._p.toNat components).map
                Prod.fst) (by
                  intro polynomial hpolynomial
                  rcases List.mem_map.mp hpolynomial with ⟨item, hitem, rfl⟩
                  exact hddfCorrect.2.2.2.2 item hitem)
            have hddfEq := _root_.eq_of_associated_monic _ _
              hsqfReady.monic hddfMonic hddfAssociated
            have hdecodedProduct :
                ((StrictDDF.ddfResultToL2 this._p.toNat components).map
                    Prod.fst).prod = componentFactors.prod := by
              simpa [componentSuffixProduct, StrictDDF.ddfResultToL2,
                List.map_map, Function.comp_def] using hmiddleProduct
            exact hddfEq.trans hdecodedProduct
          simp only [currentFactors, List.map_map, Function.comp_def]
          have hpow := listProdPow componentFactors sqfItem.2.toNat
          rw [← hpow, ← hcomponentProduct]
          rfl
        · intro item hitem
          rcases List.mem_append.mp hitem with hitem | hitem
          · rcases List.mem_map.mp hitem with ⟨q, hq, rfl⟩
            exact hmiddleQuality q hq
          · exact htailQuality item hitem
      next hindex =>
        have hle : squarefreeFactors.size ≤ index := Nat.le_of_not_gt hindex
        have hdrop : squarefreeFactors.toList.drop index = [] := by
          exact List.drop_eq_nil_iff.mpr (by simpa using hle)
        refine ⟨result, [], rfl, by simp, ?_, by simp⟩
        simp [squarefreeSuffixProduct, hdrop]

end Refinement.StrictFactorZp
