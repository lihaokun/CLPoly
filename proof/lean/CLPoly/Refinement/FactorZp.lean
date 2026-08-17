/-
  Genuine refinement of the generated C++ `__factor_Zp` entry.

  This file reasons about the strict, source-shaped L1 loops in
  `Generated.StrictFactorZp`.  Component executions are instantiated with the
  already verified strict SQF, DDF and EDF entries; no L2 algorithm is used as
  an executable fallback.
-/
import CLPoly.Generated.StrictFactorZp
import CLPoly.Refinement.SquarefreeZpEntry
import CLPoly.Refinement.DDF
import CLPoly.Refinement.EDF
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

/-- A concrete SQF component together with exactly the representation facts
needed by the strict DDF entry. -/
structure SQFOutputReady (this : DenseUPolyZp)
    (item : SparsePolyZp × UInt64) : Prop extends DDFReady this item.1 where
  multiplicityPositive : 1 ≤ item.2.toNat

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
  rcases StrictDDF.strictDDFEntryIR_refines_ddf this providers f
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
      EDFCorrect (SparsePolyZp.toPoly this._p.toNat f) d.toNat factors ∧
      ∀ factor ∈ output.toList,
        SparsePolyZp.Canonical this._p.toNat factor := by
  rcases StrictEDF.strictEDFEntryIR_refines_edf engine this providers
      termination #[] f d rng hinvariant (by simp) with
    ⟨output, rng', factors, hrun, hdecode, hcorrect, hcanonical⟩
  refine ⟨output, rng', factors, ?_, ?_, hcorrect, hcanonical⟩
  · rw [strictEDFCall_eq engine this providers termination #[] f d rng
      hinvariant]
    exact hrun
  · simpa using hdecode

/-- General (not already-monic) correctness of the exact C++ make-monic
entry used by `__factor_Zp`. -/
theorem upolyMakeMonicIR_refines
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (f : SparsePolyZp)
    (hcanonical : SparsePolyZp.Canonical this._p.toNat f)
    (hnonempty : 0 < f.size) :
    ∃ lc monic,
      StrictSquarefreeZp.upolyMakeMonicIR this f = .ok (lc, monic) ∧
      SparsePolyZp.Canonical this._p.toNat monic ∧
      (SparsePolyZp.toPoly this._p.toNat monic).Monic ∧
      SparsePolyZp.toPoly this._p.toNat f =
        Polynomial.C (Zp.toZMod this._p.toNat lc) *
          SparsePolyZp.toPoly this._p.toNat monic ∧
      monic.size = f.size ∧
      (SparsePolyZp.toPoly this._p.toNat monic).natDegree =
        (SparsePolyZp.toPoly this._p.toNat f).natDegree := by
  let lc := f[0].2
  have hlcMem : lc ∈ f.toList.map Prod.snd := by
    apply List.mem_map.mpr
    exact ⟨f[0], by simpa using Array.getElem_mem f 0 hnonempty, rfl⟩
  have hheadMem : f[0] ∈ f.toList := by
    simpa using Array.getElem_mem f 0 hnonempty
  have hlcReduced : Zp.Reduced this._p.toNat lc :=
    hcanonical.1 f[0] hheadMem
  have hlcNonzero : lc.val ≠ 0 := hcanonical.2.2 f[0] hheadMem
  have hlcPos : 0 < lc.val.toNat := Nat.pos_of_ne_zero (by
    intro hzero
    exact hlcNonzero (UInt64.toNat_inj.mp (by simpa using hzero)))
  have hpWord : this._p.toNat < UInt64.size := UInt64.toNat_lt_size this._p
  by_cases hone : lc.val = 1
  · have hrun : StrictSquarefreeZp.upolyMakeMonicIR this f = .ok (lc, f) := by
      simp [StrictSquarefreeZp.upolyMakeMonicIR, hnonempty, lc, hone]
    have hleadField : Zp.toZMod this._p.toNat lc = 1 := by
      unfold Zp.toZMod
      rw [hone]
      simp
    have hdegree := StrictSquarefreeZp.sparsePolyZp_toPoly_degree_eq_head
      this._p.toNat f hcanonical hnonempty
    have hnatDegree : (SparsePolyZp.toPoly this._p.toNat f).natDegree =
        f[0].1.deg := Polynomial.natDegree_eq_of_degree_eq_some hdegree
    have hleading : (SparsePolyZp.toPoly this._p.toNat f).leadingCoeff = 1 := by
      rw [Polynomial.leadingCoeff, hnatDegree]
      unfold SparsePolyZp.toPoly
      rw [StrictSquarefreeZp.listSum_coeff_of_mem_chain this._p.toNat
        f.toList f[0] hcanonical.2.1 hheadMem]
      exact hleadField
    refine ⟨lc, f, hrun, hcanonical, hleading, ?_, rfl, rfl⟩
    simp [hleadField]
  · let inverse := CLPoly.Impl.StrictPolynomialGCDRefinement.generatedZpInvIR
        this lc
    have hinverse := CLPoly.Impl.StrictWordArithmetic.dense_upoly_zp_nmod_inv_ir_correct
      this lc.val (Fact.out : Nat.Prime this._p.toNat) hlcPos hlcReduced.2
    change inverse.val.toNat < this._p.toNat ∧
      (inverse.val.toNat : ZMod this._p.toNat) *
        (lc.val.toNat : ZMod this._p.toNat) = 1 at hinverse
    have hinverseReduced : Zp.Reduced this._p.toNat inverse := by
      exact ⟨hlcReduced.1, hinverse.1⟩
    have hinverseNonzero : inverse.val ≠ 0 := by
      intro hzero
      have hfieldZero : (inverse.val.toNat : ZMod this._p.toNat) = 0 := by
        rw [hzero]
        simp
      have hbad := hinverse.2
      rw [hfieldZero, zero_mul] at hbad
      exact zero_ne_one hbad
    let monic := CLPoly.Impl.StrictPolynomialGCDRefinement.sparseMonicLoop
      0 f inverse
    have hrun : StrictSquarefreeZp.upolyMakeMonicIR this f =
        .ok (lc, monic) := by
      simp [StrictSquarefreeZp.upolyMakeMonicIR, hnonempty, lc, hone,
        inverse, monic]
    have hcanonicalMonic : SparsePolyZp.Canonical this._p.toNat monic :=
      CLPoly.Impl.StrictPolynomialGCDRefinement.sparseMonicLoop_canonical
        this._p.toNat f inverse hcanonical hinverseReduced hinverseNonzero
        hpWord
    have hsemantic : SparsePolyZp.toPoly this._p.toNat monic =
        Polynomial.C (Zp.toZMod this._p.toNat inverse) *
          SparsePolyZp.toPoly this._p.toNat f :=
      CLPoly.Impl.StrictPolynomialGCDRefinement.sparseMonicLoop_toPoly
        this._p.toNat f inverse hcanonical hinverseReduced.1
        (Fact.out : Nat.Prime this._p.toNat).pos hpWord
    have hlcInverse : Zp.toZMod this._p.toNat inverse *
        Zp.toZMod this._p.toNat lc = 1 := by
      unfold Zp.toZMod
      simpa [inverse] using hinverse.2
    have hfDegree := StrictSquarefreeZp.sparsePolyZp_toPoly_degree_eq_head
      this._p.toNat f hcanonical hnonempty
    have hmonicDegree : (SparsePolyZp.toPoly this._p.toNat monic).degree =
        f[0].1.deg := by
      rw [hsemantic, Polynomial.degree_mul]
      · rw [Polynomial.degree_C]
        · simp [hfDegree]
        · exact Zp.toZMod_ne_zero_of_val_ne_zero this._p.toNat inverse
            hinverseReduced hinverseNonzero
    have hmonicNatDegree :
        (SparsePolyZp.toPoly this._p.toNat monic).natDegree = f[0].1.deg :=
      Polynomial.natDegree_eq_of_degree_eq_some hmonicDegree
    have hmonic : (SparsePolyZp.toPoly this._p.toNat monic).Monic := by
      rw [Polynomial.Monic.def, Polynomial.leadingCoeff,
        hmonicNatDegree, hsemantic]
      rw [Polynomial.coeff_C_mul]
      unfold SparsePolyZp.toPoly
      rw [StrictSquarefreeZp.listSum_coeff_of_mem_chain this._p.toNat
        f.toList f[0] hcanonical.2.1 hheadMem]
      exact hlcInverse
    have hreconstruct : SparsePolyZp.toPoly this._p.toNat f =
        Polynomial.C (Zp.toZMod this._p.toNat lc) *
          SparsePolyZp.toPoly this._p.toNat monic := by
      rw [hsemantic, ← mul_assoc, ← Polynomial.C_mul]
      have hcomm : Zp.toZMod this._p.toNat lc *
          Zp.toZMod this._p.toNat inverse = 1 := by
        simpa [mul_comm] using hlcInverse
      rw [hcomm, map_one, one_mul]
    have hsize : monic.size = f.size := by
      exact CLPoly.Impl.StrictPolynomialGCDRefinement.sparseMonicLoop_size
        0 f inverse
    have hfNatDegree :
        (SparsePolyZp.toPoly this._p.toNat f).natDegree = f[0].1.deg :=
      Polynomial.natDegree_eq_of_degree_eq_some hfDegree
    exact ⟨lc, monic, hrun, hcanonicalMonic, hmonic, hreconstruct, hsize,
      hmonicNatDegree.trans hfNatDegree.symm⟩

/-- Entry facts for the strict SQF call made by the top-level factorization. -/
structure SQFReady (this : DenseUPolyZp) (source : SparsePolyZp) : Prop where
  canonical : SparsePolyZp.Canonical this._p.toNat source
  monic : (SparsePolyZp.toPoly this._p.toNat source).Monic
  nonempty : 0 < source.size
  degreePositive : 0 < (SparsePolyZp.toPoly this._p.toNat source).natDegree
  denseBound : CLPoly.Impl.StrictPolynomialGCDRefinement.sparseDenseLength
    source ≤ 2 ^ 63
  sourceDegreeBound :
    (SparsePolyZp.toPoly this._p.toNat source).natDegree < 2 ^ 62

noncomputable def strictSQFCall
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (physical : StrictSquarefreeZp.YunRawGCDWorkspaceProvider this hcfg)
    (source : SparsePolyZp) : RawExec (Array (SparsePolyZp × UInt64)) := by
  classical
  exact if hready : SQFReady this source then
    Generated.StrictSquarefreeZp.__squarefree_Zp_raw_ir
      (StrictSquarefreeEntry.strictSQFRawOps this hcfg physical) source
      (fun _ => ⟨hready.canonical, hready.monic, hready.nonempty,
        hready.degreePositive, hready.denseBound⟩)
  else .error .assertionFailure

theorem strictSQFCall_refines
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (physical : StrictSquarefreeZp.YunRawGCDWorkspaceProvider this hcfg)
    (source : SparsePolyZp) (hready : SQFReady this source) :
    ∃ output,
      strictSQFCall this hcfg physical source = .ok output ∧
      toPolyList output this._p.toNat =
        sqfZp (SparsePolyZp.toPoly this._p.toNat source) ∧
      ∀ item ∈ output.toList, SQFOutputReady this item := by
  rcases StrictSquarefreeEntry.__squarefree_Zp_raw_ir_refines_sqfZp
      this hcfg physical
      source hready.canonical hready.monic hready.nonempty
      hready.degreePositive hready.denseBound with
    ⟨output, hrun, hdecode, hcanonical⟩
  have hsquarefreeCorrect := sqf_correct
    (SparsePolyZp.toPoly this._p.toNat source) hready.monic.ne_zero
  refine ⟨output, ?_, hdecode, ?_⟩
  · classical
    simp only [strictSQFCall, dif_pos hready]
    exact hrun
  · intro item hitem
    have hdecoded :
        (SparsePolyZp.toPoly this._p.toNat item.1, item.2.toNat) ∈
          sqfZp (SparsePolyZp.toPoly this._p.toNat source) := by
      rw [← hdecode]
      rw [toPolyList, Array.toList_map]
      exact List.mem_map.mpr ⟨item, hitem, by cases item; rfl⟩
    have hcanonicalItem := hcanonical item hitem
    let decodedItem : Polynomial (ZMod this._p.toNat) × Nat :=
      (SparsePolyZp.toPoly this._p.toNat item.1, item.2.toNat)
    have hdecoded' : decodedItem ∈
        sqfZp (SparsePolyZp.toPoly this._p.toNat source) := hdecoded
    have hmonicItem := (hsquarefreeCorrect.2.1 decodedItem hdecoded').2
    have hsquarefreeItem := (hsquarefreeCorrect.2.1 decodedItem hdecoded').1
    have hmultPositive := hsquarefreeCorrect.2.2.1 decodedItem hdecoded'
    have hitemDvd : SparsePolyZp.toPoly this._p.toNat item.1 ∣
        SparsePolyZp.toPoly this._p.toNat source := by
      have hpowMem : (SparsePolyZp.toPoly this._p.toNat item.1) ^
            item.2.toNat ∈
          (sqfZp (SparsePolyZp.toPoly this._p.toNat source)).map
            (fun factor => factor.1 ^ factor.2) :=
        List.mem_map.mpr ⟨_, hdecoded, rfl⟩
      exact (dvd_pow_self _ (Nat.one_le_iff_ne_zero.mp hmultPositive)).trans
        ((List.dvd_prod hpowMem).trans hsquarefreeCorrect.1.symm.dvd)
    have hdegreeItem :
        (SparsePolyZp.toPoly this._p.toNat item.1).natDegree < 2 ^ 62 :=
      lt_of_le_of_lt
        (Polynomial.natDegree_le_of_dvd hitemDvd hready.monic.ne_zero)
        hready.sourceDegreeBound
    have hnonemptyItem :=
      StrictSquarefreeZp.sparsePolyZp_size_pos_of_toPoly_ne_zero
        this._p.toNat item.1 hmonicItem.ne_zero
    have hheadMem : item.1[0] ∈ item.1.toList := by
      simpa using Array.getElem_mem item.1 0 hnonemptyItem
    have hhead : item.1[0]! = item.1[0] := by
      rw [getElem!_def, getElem?_def]
      simp [hnonemptyItem]
    refine ⟨⟨?_, hcanonicalItem, ?_, hmonicItem, hsquarefreeItem⟩,
      hmultPositive⟩
    · rw [hhead]
      exact UInt64.toNat_inj.mp (hcanonicalItem.1 item.1[0] hheadMem).1
    · intro term hterm
      exact lt_of_le_of_lt
        (StrictDDF.canonical_term_degree_le_natDegree this._p.toNat item.1
          hcanonicalItem term hterm) hdegreeItem

/-- Decode the concrete output array of C++ `__factor_Zp`. -/
noncomputable def factorResultToL2 (p : Nat)
    (result : Array (SparsePolyZp × UInt64)) :
    List (Polynomial (ZMod p) × Nat) :=
  result.toList.map fun item =>
    (SparsePolyZp.toPoly p item.1, item.2.toNat)

/-- Physical contract of the C++ `std::sort` library call used by
`__factor_Zp`.  It may choose any concrete sorting implementation, but its
returned buffer must be a permutation of the input buffer.  No factorization
result or polynomial property is supplied by this boundary. -/
structure SortByDegreeProvider where
  run : Array (SparsePolyZp × UInt64) →
    RawExec (Array (SparsePolyZp × UInt64))
  permutation : ∀ input, ∃ output,
    run input = .ok output ∧ output.toList.Perm input.toList

theorem factorResultToL2_perm {p : Nat}
    {left right : Array (SparsePolyZp × UInt64)}
    (hperm : left.toList.Perm right.toList) :
    (factorResultToL2 p left).Perm (factorResultToL2 p right) := by
  exact hperm.map (fun item =>
    (SparsePolyZp.toPoly p item.1, item.2.toNat))

theorem factorZpCorrect_perm {p : Nat} [Fact (Nat.Prime p)]
    {f : Polynomial (ZMod p)} {lc : ZMod p}
    {left right : List (Polynomial (ZMod p) × Nat)}
    (hcorrect : FactorZpCorrect f lc left) (hperm : right.Perm left) :
    FactorZpCorrect f lc right := by
  refine ⟨?_, ?_⟩
  · rw [hcorrect.1]
    rw [(hperm.map fun item => item.1 ^ item.2).prod_eq]
  · intro item hitem
    exact hcorrect.2 item (hperm.mem_iff.mp hitem)

theorem sparseDenseLength_eq_natDegree_add_one
    (p : Nat) [Fact (Nat.Prime p)] (poly : SparsePolyZp)
    (hcanonical : SparsePolyZp.Canonical p poly)
    (hnonempty : 0 < poly.size) :
    CLPoly.Impl.StrictPolynomialGCDRefinement.sparseDenseLength poly =
      (SparsePolyZp.toPoly p poly).natDegree + 1 := by
  rw [CLPoly.Impl.StrictPolynomialGCDRefinement.sparseDenseLength,
    dif_pos hnonempty]
  have hdegree := StrictSquarefreeZp.sparsePolyZp_toPoly_degree_eq_head
    p poly hcanonical hnonempty
  rw [Polynomial.natDegree_eq_of_degree_eq_some hdegree]

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

theorem factorLoop0_preserves_canonical
    (p : Nat) (factors : Array SparsePolyZp) (multiplicity : UInt64)
    (index : Nat) (result : Array (SparsePolyZp × UInt64))
    (hfactors : ∀ factor ∈ factors.toList, SparsePolyZp.Canonical p factor)
    (hresult : ∀ item ∈ result.toList, SparsePolyZp.Canonical p item.1) :
    ∀ item ∈ (Generated.StrictFactorZp._loop___factor_Zp_0_raw_ir
      factors multiplicity index result).toList,
      SparsePolyZp.Canonical p item.1 := by
  induction hmeasure : factors.size - index using Nat.strong_induction_on
      generalizing index result with
  | h measure ih =>
      rw [Generated.StrictFactorZp._loop___factor_Zp_0_raw_ir]
      split
      next hindex =>
        apply ih (factors.size - (index + 1)) (by omega)
          (index + 1) (result.push (factors[index], multiplicity))
        · intro item hitem
          rw [Array.toList_push, List.mem_append] at hitem
          rcases hitem with hitem | hitem
          · exact hresult item hitem
          · have heq : item = (factors[index], multiplicity) := by
              simpa using hitem
            subst item
            exact hfactors factors[index] (by
              simpa using Array.getElem_mem factors index hindex)
        · rfl
      next _ => exact hresult

/-- Genuine refinement of the generated middle loop.  Each EDF factor list is
obtained from the actual strict call at the current DDF array element. -/
theorem factorLoop1_refines {State : Type}
    (engine : Generated.StrictEDF.RandomEngine State)
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (providers : StrictDDF.DDFRawProviders this)
    (termination : Generated.StrictEDF.EDFTermination
      (StrictEDF.strictEDFRawOps engine this providers))
    (makeMonic : SparsePolyZp → RawExec (Zp × SparsePolyZp))
    (squarefree : SparsePolyZp → RawExec (Array (SparsePolyZp × UInt64)))
    (sortByDegree : Array (SparsePolyZp × UInt64) →
      RawExec (Array (SparsePolyZp × UInt64)))
    (components : Array (SparsePolyZp × UInt64))
    (hready : ∀ item ∈ components.toList,
      StrictEDF.EDFEntryInvariant this item.1 item.2)
    (multiplicity : UInt64) (index : Nat) (rng : State)
    (result : Array (SparsePolyZp × UInt64))
    (hresultCanonical : ∀ item ∈ result.toList,
      SparsePolyZp.Canonical this._p.toNat item.1) :
    let ops : Generated.StrictFactorZp.FactorZpRawOps State := {
      makeMonic := makeMonic
      squarefree := squarefree
      ddf := strictDDFCall this providers
      edf := strictEDFCall engine this providers termination
      sortByDegree := sortByDegree }
    ∃ (output : Array (SparsePolyZp × UInt64)) (rng' : State)
        (factors : List (Polynomial (ZMod this._p.toNat))),
      Generated.StrictFactorZp._loop___factor_Zp_1_raw_ir
          ops components multiplicity index rng result = .ok (output, rng') ∧
      factorResultToL2 this._p.toNat output =
        factorResultToL2 this._p.toNat result ++
          (factors.map fun q => (q, multiplicity.toNat)) ∧
      componentSuffixProduct this._p.toNat components index = factors.prod ∧
      (∀ q ∈ factors, Irreducible q ∧ Monic q) ∧
      ∀ item ∈ output.toList,
        SparsePolyZp.Canonical this._p.toNat item.1 := by
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
            hedfCorrect, hedfCanonical⟩
        change ∃ (output : Array (SparsePolyZp × UInt64)) (rng' : State)
            (factors : List (Polynomial (ZMod this._p.toNat))),
          (match strictEDFCall engine this providers termination #[]
              component.1 component.2 rng with
            | Except.error fault => Except.error fault
            | Except.ok (factors, rng') =>
              Generated.StrictFactorZp._loop___factor_Zp_1_raw_ir
                { makeMonic := makeMonic
                  squarefree := squarefree
                  ddf := strictDDFCall this providers
                  edf := strictEDFCall engine this providers termination
                  sortByDegree := sortByDegree }
                components multiplicity (index + 1) rng'
                (Generated.StrictFactorZp._loop___factor_Zp_0_raw_ir
                  factors multiplicity 0 result)) = Except.ok (output, rng') ∧
          factorResultToL2 this._p.toNat output =
            factorResultToL2 this._p.toNat result ++
              (factors.map fun q => (q, multiplicity.toNat)) ∧
          componentSuffixProduct this._p.toNat components index =
            factors.prod ∧
          (∀ q ∈ factors, Irreducible q ∧ Monic q) ∧
          ∀ item ∈ output.toList,
            SparsePolyZp.Canonical this._p.toNat item.1
        rw [hedfRun]
        rcases ih (components.size - (index + 1)) (by omega)
            (index + 1) rngNext
            (Generated.StrictFactorZp._loop___factor_Zp_0_raw_ir
              edfOutput multiplicity 0 result)
            (factorLoop0_preserves_canonical this._p.toNat edfOutput
              multiplicity 0 result hedfCanonical hresultCanonical) rfl with
          ⟨output, rng', tailFactors, htailRun, htailDecode,
            htailProduct, htailQuality, houtputCanonical⟩
        refine ⟨output, rng', edfFactors ++ tailFactors, htailRun, ?_, ?_, ?_,
          houtputCanonical⟩
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
        refine ⟨result, rng, [], rfl, by simp, ?_, by simp, hresultCanonical⟩
        simp [componentSuffixProduct, hdrop]

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
    (makeMonic : SparsePolyZp → RawExec (Zp × SparsePolyZp))
    (squarefree : SparsePolyZp → RawExec (Array (SparsePolyZp × UInt64)))
    (sortByDegree : Array (SparsePolyZp × UInt64) →
      RawExec (Array (SparsePolyZp × UInt64)))
    (squarefreeFactors : Array (SparsePolyZp × UInt64))
    (hready : ∀ item ∈ squarefreeFactors.toList,
      SQFOutputReady this item)
    (index : Nat) (rng : State)
    (result : Array (SparsePolyZp × UInt64))
    (hresultCanonical : ∀ item ∈ result.toList,
      SparsePolyZp.Canonical this._p.toNat item.1) :
    let ops : Generated.StrictFactorZp.FactorZpRawOps State := {
      makeMonic := makeMonic
      squarefree := squarefree
      ddf := strictDDFCall this providers
      edf := strictEDFCall engine this providers termination
      sortByDegree := sortByDegree }
    ∃ (output : Array (SparsePolyZp × UInt64))
        (factors : List (Polynomial (ZMod this._p.toNat) × Nat)),
      Generated.StrictFactorZp._loop___factor_Zp_2_raw_ir
          ops squarefreeFactors index rng result = .ok output ∧
      factorResultToL2 this._p.toNat output =
        factorResultToL2 this._p.toNat result ++ factors ∧
      squarefreeSuffixProduct this._p.toNat squarefreeFactors index =
        (factors.map fun item => item.1 ^ item.2).prod ∧
      (∀ item ∈ factors,
        Irreducible item.1 ∧ Monic item.1 ∧ 1 ≤ item.2) ∧
      ∀ item ∈ output.toList,
        SparsePolyZp.Canonical this._p.toNat item.1 := by
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
                  { makeMonic := makeMonic
                    squarefree := squarefree
                    ddf := strictDDFCall this providers
                    edf := strictEDFCall engine this providers termination
                    sortByDegree := sortByDegree }
                  components sqfItem.2 0 rng result with
              | Except.error fault => Except.error fault
              | Except.ok (result', rng') =>
                Generated.StrictFactorZp._loop___factor_Zp_2_raw_ir
                  { makeMonic := makeMonic
                    squarefree := squarefree
                    ddf := strictDDFCall this providers
                    edf := strictEDFCall engine this providers termination
                    sortByDegree := sortByDegree }
                  squarefreeFactors (index + 1) rng' result') =
              Except.ok output ∧
          factorResultToL2 this._p.toNat output =
            factorResultToL2 this._p.toNat result ++ factors ∧
          squarefreeSuffixProduct this._p.toNat squarefreeFactors index =
            (factors.map fun item => item.1 ^ item.2).prod ∧
          (∀ item ∈ factors,
            Irreducible item.1 ∧ Monic item.1 ∧ 1 ≤ item.2) ∧
          ∀ item ∈ output.toList,
            SparsePolyZp.Canonical this._p.toNat item.1
        rw [hddfRun]
        rcases factorLoop1_refines engine this providers termination makeMonic
            squarefree sortByDegree components hedfReady sqfItem.2 0 rng result
            hresultCanonical with
          ⟨middleOutput, rngNext, componentFactors, hmiddleRun,
            hmiddleDecode, hmiddleProduct, hmiddleQuality, hmiddleCanonical⟩
        simp only
        rw [hmiddleRun]
        rcases ih (squarefreeFactors.size - (index + 1)) (by omega)
            (index + 1) rngNext middleOutput hmiddleCanonical rfl with
          ⟨output, tailFactors, htailRun, htailDecode,
            htailProduct, htailQuality, houtputCanonical⟩
        let currentFactors := componentFactors.map fun q =>
          (q, sqfItem.2.toNat)
        refine ⟨output, currentFactors ++ tailFactors, htailRun, ?_, ?_, ?_,
          houtputCanonical⟩
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
            exact ⟨(hmiddleQuality q hq).1, (hmiddleQuality q hq).2,
              hsqfReady.multiplicityPositive⟩
          · exact htailQuality item hitem
      next hindex =>
        have hle : squarefreeFactors.size ≤ index := Nat.le_of_not_gt hindex
        have hdrop : squarefreeFactors.toList.drop index = [] := by
          exact List.drop_eq_nil_iff.mpr (by simpa using hle)
        refine ⟨result, [], rfl, by simp, ?_, by simp, hresultCanonical⟩
        simp [squarefreeSuffixProduct, hdrop]

/-- End-to-end refinement of the exact generated C++ `__factor_Zp` control
flow.  Every executable component is a strict L1 entry; the only library
boundary is `std::sort`, whose contract exposes only buffer permutation. -/
theorem __factor_Zp_raw_ir_refines_FactorZpCorrect
    {State : Type} (engine : Generated.StrictEDF.RandomEngine State)
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (sqfPhysical : StrictSquarefreeZp.YunRawGCDWorkspaceProvider this hcfg)
    (providers : StrictDDF.DDFRawProviders this)
    (termination : Generated.StrictEDF.EDFTermination
      (StrictEDF.strictEDFRawOps engine this providers))
    (sort : SortByDegreeProvider) (initialRng : State)
    (f : SparsePolyZp)
    (hcanonical : SparsePolyZp.Canonical this._p.toNat f)
    (hnonempty : 0 < f.size)
    (hdegreePositive :
      0 < (SparsePolyZp.toPoly this._p.toNat f).natDegree)
    (hdegreeBound :
      (SparsePolyZp.toPoly this._p.toNat f).natDegree < 2 ^ 62) :
    let ops : Generated.StrictFactorZp.FactorZpRawOps State := {
      makeMonic := StrictSquarefreeZp.upolyMakeMonicIR this
      squarefree := strictSQFCall this hcfg sqfPhysical
      ddf := strictDDFCall this providers
      edf := strictEDFCall engine this providers termination
      sortByDegree := sort.run }
    ∃ lc output,
      Generated.StrictFactorZp.__factor_Zp_raw_ir ops initialRng f =
        .ok (lc, output) ∧
      FactorZpCorrect (SparsePolyZp.toPoly this._p.toNat f)
        (Zp.toZMod this._p.toNat lc)
        (factorResultToL2 this._p.toNat output) ∧
      (∀ item ∈ output.toList,
        SparsePolyZp.Canonical this._p.toNat item.1) := by
  dsimp only
  have hempty : f.isEmpty = false := by
    simpa [Array.isEmpty] using (Nat.ne_of_gt hnonempty)
  have hgetDeg : get_deg f > 0 :=
    (StrictDDF.strict_get_deg_pos_iff_natDegree_pos this._p.toNat f
      hcanonical (lt_trans hdegreeBound (by omega))).2 hdegreePositive
  rcases upolyMakeMonicIR_refines this f hcanonical hnonempty with
    ⟨lc, monic, hmonicRun, hmonicCanonical, hmonic, hreconstruct,
      hmonicSize, hmonicDegree⟩
  have hmonicNonempty : 0 < monic.size := by omega
  have hmonicDenseBound :
      CLPoly.Impl.StrictPolynomialGCDRefinement.sparseDenseLength monic ≤
        2 ^ 63 := by
    rw [sparseDenseLength_eq_natDegree_add_one this._p.toNat monic
      hmonicCanonical hmonicNonempty, hmonicDegree]
    omega
  have hsqfReady : SQFReady this monic :=
    ⟨hmonicCanonical, hmonic, hmonicNonempty,
      hmonicDegree ▸ hdegreePositive, hmonicDenseBound,
      hmonicDegree ▸ hdegreeBound⟩
  rcases strictSQFCall_refines this hcfg sqfPhysical monic hsqfReady with
    ⟨sqfOutput, hsqfRun, hsqfDecode, hsqfOutputReady⟩
  rcases factorLoop2_refines engine this providers termination
      (StrictSquarefreeZp.upolyMakeMonicIR this)
      (strictSQFCall this hcfg sqfPhysical) sort.run sqfOutput
      hsqfOutputReady 0 initialRng #[] (by simp) with
    ⟨unsorted, factors, hloops, hdecode, hproduct, hquality,
      hunsortedCanonical⟩
  rcases sort.permutation unsorted with ⟨sorted, hsort, hpermRaw⟩
  have hunsortedDecode :
      factorResultToL2 this._p.toNat unsorted = factors := by
    simpa using hdecode
  have hpermDecoded :
      (factorResultToL2 this._p.toNat sorted).Perm factors := by
    have hp := factorResultToL2_perm (p := this._p.toNat) hpermRaw
    rw [hunsortedDecode] at hp
    exact hp
  have hsortedCanonical :
      ∀ item ∈ sorted.toList,
        SparsePolyZp.Canonical this._p.toNat item.1 := by
    intro item hitem
    exact hunsortedCanonical item (hpermRaw.mem_iff.mp hitem)
  have hcorrectUnsorted :
      FactorZpCorrect (SparsePolyZp.toPoly this._p.toNat f)
        (Zp.toZMod this._p.toNat lc) factors := by
    refine ⟨?_, ?_⟩
    · rw [hreconstruct]
      congr 1
      have hsqfCorrect := sqf_correct
        (SparsePolyZp.toPoly this._p.toNat monic) hmonic.ne_zero
      let powered : Polynomial (ZMod this._p.toNat) × Nat →
          Polynomial (ZMod this._p.toNat) := fun item => item.1 ^ item.2
      have hsqfProductMonic :
          ((sqfZp (SparsePolyZp.toPoly this._p.toNat monic)).map
            powered).prod.Monic := by
        exact monicListProd _ (by
          intro q hq
          rcases List.mem_map.mp hq with
            ⟨item : Polynomial (ZMod this._p.toNat) × Nat, hitem, rfl⟩
          exact (hsqfCorrect.2.1 item hitem).2.pow item.2)
      have hrawProduct :
          ((sqfZp (SparsePolyZp.toPoly this._p.toNat monic)).map
            (fun item => item.1 ^ item.2)).prod =
          (factors.map fun item => item.1 ^ item.2).prod := by
        have hraw := hproduct
        rw [squarefreeSuffixProduct, List.drop_zero] at hraw
        have hdecodedMap :
            sqfOutput.toList.map (fun item =>
              SparsePolyZp.toPoly this._p.toNat item.1 ^ item.2.toNat) =
            (toPolyList sqfOutput this._p.toNat).map
              (fun item => item.1 ^ item.2) := by
          simp [toPolyList, Array.toList_map, List.map_map,
            Function.comp_def]
        rw [hdecodedMap, hsqfDecode] at hraw
        exact hraw
      have hsqfEq := _root_.eq_of_associated_monic _ _ hmonic
        hsqfProductMonic hsqfCorrect.1
      exact hsqfEq.trans hrawProduct
    · exact hquality
  refine ⟨lc, sorted, ?_,
    factorZpCorrect_perm hcorrectUnsorted hpermDecoded, hsortedCanonical⟩
  simp only [Generated.StrictFactorZp.__factor_Zp_raw_ir, hempty,
    Bool.false_eq_true, ↓reduceIte]
  have hnotConstant : ¬get_deg f ≤ 0 := by
    intro hle
    have hltInt : (0 : Int) < (get_deg f).toInt := by
      simpa only [Int64.lt_iff_toInt_lt] using hgetDeg
    have hleInt : (get_deg f).toInt ≤ (0 : Int) := by
      simpa only [Int64.le_iff_toInt_le] using hle
    omega
  rw [if_neg hnotConstant, hmonicRun]
  simp only
  rw [hsqfRun]
  simp only
  rw [hloops]
  simp only
  rw [hsort]

end Refinement.StrictFactorZp
