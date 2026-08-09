/-
  Strict squarefree refinement boundary.

  The earlier proof operated on legacy sparse-array entries whose division and
  GCD primitives could reduce definitionally to hand-written L2 algorithms.
  It therefore did not prove the raw C++ dense-polynomial execution and is not
  exported.

  Squarefree refinement will resume after strict RawHeap division and Euclidean
  GCD are connected to their L2 polynomial specifications.  All recursive
  source loops must remain well-founded; no bounded execution wrapper is part
  of this boundary.
-/
import CLPoly.Algorithm.SquarefreeZp
import CLPoly.Impl.StrictPolynomialGCDRefinement
import CLPoly.Impl.StrictPolyAddSubRefinement
import CLPoly.Refinement.Basic

set_option autoImplicit false

namespace Refinement

namespace StrictSquarefreeZp

open CLPoly.Impl.StrictPolynomialGCDRefinement
open CLPoly.Math

/-- Exact checked entry for the source `__upoly_make_monic`. -/
def upolyMakeMonicIR (this : DenseUPolyZp) (f : SparsePolyZp) :
    RawExec (Zp × SparsePolyZp) :=
  if hnonempty : 0 < f.size then
    let lc := f[0].2
    if lc.val == 1 then
      .ok (lc, f)
    else
      let lcInv := generatedZpInvIR this lc
      .ok (lc, sparseMonicLoop 0 f lcInv)
  else
    .error .assertionFailure

/-- A concrete monic sparse input takes the actual early-return comparison in
`__upoly_make_monic`; no mutation loop is erased from the definition. -/
theorem upolyMakeMonicIR_eq_of_monic (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (f : SparsePolyZp)
    (hcanonical : SparsePolyZp.Canonical this._p.toNat f)
    (hnonempty : 0 < f.size)
    (hmonic : (SparsePolyZp.toPoly this._p.toNat f).Monic) :
    upolyMakeMonicIR this f = .ok (f[0].2, f) := by
  have hlead := sparse_leading_val_eq_one_of_monic this._p.toNat f hcanonical
    hnonempty hmonic
  simp [upolyMakeMonicIR, hnonempty, hlead]

/-- One observable iteration of the current C++ `derivative` template.  The
integer degree is first constructed as a `Zp` residue, then the generated
normalized modular multiplier is executed. -/
def derivativeValueIR (this : DenseUPolyZp) (term : UMonomial × Zp) : UInt64 :=
  let degreeResidue := term.1.deg.toUInt64 % this._p
  Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this term.2.val degreeResidue

def derivativeTermIR (this : DenseUPolyZp) (term : UMonomial × Zp) :
    Option (UMonomial × Zp) :=
  if term.1.deg = 0 then
    none
  else
    let value := derivativeValueIR this term
    if value = 0 then none
    else some (⟨term.1.deg - 1⟩, ⟨value, term.2.prime⟩)

def derivativeLoopIR (this : DenseUPolyZp) (index : Nat)
    (output source : SparsePolyZp) : SparsePolyZp :=
  if h : index < source.size then
    let next := match derivativeTermIR this source[index] with
      | none => output
      | some term => output.push term
    derivativeLoopIR this (index + 1) next source
  else
    output
termination_by source.size - index
decreasing_by omega

def derivativeIR (this : DenseUPolyZp) (source : SparsePolyZp) :
    SparsePolyZp := derivativeLoopIR this 0 #[] source

theorem derivativeLoopIR_toList (this : DenseUPolyZp) (index : Nat)
    (output source : SparsePolyZp) :
    (derivativeLoopIR this index output source).toList =
      output.toList ++ (source.toList.drop index).filterMap
        (derivativeTermIR this) := by
  rw [derivativeLoopIR]
  split
  next hmore =>
    rw [derivativeLoopIR_toList]
    have hdrop := List.drop_eq_getElem_cons
      (l := source.toList) (i := index) (by simpa using hmore)
    rw [hdrop, List.filterMap_cons, Array.getElem_toList hmore]
    cases hterm : derivativeTermIR this source[index] with
    | none => simp
    | some term => simp [List.append_assoc]
  next hdone =>
    have hindex : source.toList.length ≤ index := by simpa using hdone
    rw [List.drop_eq_nil_iff.mpr hindex]
    simp
termination_by source.size - index
decreasing_by omega

theorem derivativeIR_toList (this : DenseUPolyZp) (source : SparsePolyZp) :
    (derivativeIR this source).toList =
      source.toList.filterMap (derivativeTermIR this) := by
  simpa [derivativeIR] using derivativeLoopIR_toList this 0 #[] source

theorem derivativeValueIR_toZMod (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)] (term : UMonomial × Zp)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hreduced : Zp.Reduced this._p.toNat term.2)
    (hdegree : term.1.deg < UInt64.size) :
    ((derivativeValueIR this term).toNat : ZMod this._p.toNat) =
      Zp.toZMod this._p.toNat term.2 * term.1.deg := by
  let degreeWord := term.1.deg.toUInt64
  let degreeResidue := degreeWord % this._p
  have hdegreeWord : degreeWord.toNat = term.1.deg := by
    change (OfNat.ofNat term.1.deg : UInt64).toNat = term.1.deg
    rw [UInt64.toNat_ofNat, Nat.mod_eq_of_lt hdegree]
  have hdegreeResidue : degreeResidue.toNat =
      term.1.deg % this._p.toNat := by
    simp [degreeResidue, UInt64.toNat_mod, hdegreeWord]
  have hvalue :=
    CLPoly.Impl.StrictWordArithmetic.nmod_mul_ir_correct_of_configured this
      term.2.val degreeResidue hcfg hreduced.2
  change (Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this term.2.val
    degreeResidue).toNat = _ at hvalue
  change ((Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this term.2.val
    degreeResidue).toNat : ZMod this._p.toNat) = _
  rw [hvalue]
  unfold Zp.toZMod
  rw [hdegreeResidue, ZMod.natCast_mod]
  simp [Nat.cast_mul, ZMod.natCast_mod]

theorem listSum_derivativeTermIR (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this) :
    ∀ terms : List (UMonomial × Zp),
      (∀ term ∈ terms, Zp.Reduced this._p.toNat term.2) →
      (∀ term ∈ terms, term.1.deg < UInt64.size) →
      listSum this._p.toNat
          (terms.filterMap (derivativeTermIR this)) =
        Polynomial.derivative (listSum this._p.toNat terms) := by
  intro terms hreduced hdegree
  induction terms with
  | nil => simp [listSum]
  | cons term rest ih =>
      have htermReduced := hreduced term List.mem_cons_self
      have hrestReduced : ∀ item ∈ rest,
          Zp.Reduced this._p.toNat item.2 := by
        intro item hitem
        exact hreduced item (List.mem_cons_of_mem term hitem)
      have htermDegree := hdegree term List.mem_cons_self
      have hrestDegree : ∀ item ∈ rest, item.1.deg < UInt64.size := by
        intro item hitem
        exact hdegree item (List.mem_cons_of_mem term hitem)
      have ih' := ih hrestReduced hrestDegree
      by_cases hdegreeZero : term.1.deg = 0
      · have htermIR : derivativeTermIR this term = none := by
          simp [derivativeTermIR, hdegreeZero]
        rw [List.filterMap_cons, htermIR, listSum_cons,
          Polynomial.derivative_add, Polynomial.derivative_monomial, ih']
        simp [hdegreeZero]
      · let value := derivativeValueIR this term
        have hvalueCast := derivativeValueIR_toZMod this term hcfg
          htermReduced htermDegree
        change (value.toNat : ZMod this._p.toNat) =
          Zp.toZMod this._p.toNat term.2 * term.1.deg at hvalueCast
        by_cases hvalueZero : value = 0
        · have hcoefficientZero :
              Zp.toZMod this._p.toNat term.2 * term.1.deg = 0 := by
            rw [← hvalueCast, hvalueZero]
            simp
          have htermIR : derivativeTermIR this term = none := by
            simp [derivativeTermIR, hdegreeZero, value, hvalueZero]
          rw [List.filterMap_cons, htermIR, listSum_cons,
            Polynomial.derivative_add, Polynomial.derivative_monomial, ih',
            hcoefficientZero]
          simp
        · have htermIR : derivativeTermIR this term =
              some (⟨term.1.deg - 1⟩,
                ⟨value, term.2.prime⟩) := by
            simp [derivativeTermIR, hdegreeZero, value, hvalueZero]
          rw [List.filterMap_cons, htermIR, listSum_cons, ih',
            listSum_cons, Polynomial.derivative_add,
            Polynomial.derivative_monomial]
          change Polynomial.monomial (term.1.deg - 1)
              (value.toNat : ZMod this._p.toNat) +
              Polynomial.derivative (listSum this._p.toNat rest) = _
          rw [hvalueCast]

theorem derivativeIR_toPoly (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (source : SparsePolyZp)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hcanonical : SparsePolyZp.Canonical this._p.toNat source)
    (hdegree : ∀ term ∈ source.toList, term.1.deg < UInt64.size) :
    SparsePolyZp.toPoly this._p.toNat (derivativeIR this source) =
      Polynomial.derivative (SparsePolyZp.toPoly this._p.toNat source) := by
  unfold SparsePolyZp.toPoly
  rw [derivativeIR_toList]
  exact listSum_derivativeTermIR this hcfg source.toList hcanonical.1 hdegree

theorem derivativeTermIR_some_degree (this : DenseUPolyZp)
    (term output : UMonomial × Zp)
    (hrun : derivativeTermIR this term = some output) :
    output.1.deg = term.1.deg - 1 ∧ 0 < term.1.deg := by
  unfold derivativeTermIR at hrun
  split at hrun
  next hzero => simp at hrun
  next hpositive =>
    dsimp only at hrun
    split at hrun
    next hvalueZero => simp at hrun
    next hvalueNonzero =>
      simp only [Option.some.injEq] at hrun
      subst output
      exact ⟨rfl, Nat.pos_of_ne_zero hpositive⟩

theorem derivativeTermIR_some_properties (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (term output : UMonomial × Zp)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hreduced : Zp.Reduced this._p.toNat term.2)
    (hrun : derivativeTermIR this term = some output) :
    output.1.deg = term.1.deg - 1 ∧
      output.2.prime.toNat = this._p.toNat ∧
      output.2.val ≠ 0 ∧ output.2.val.toNat < this._p.toNat := by
  unfold derivativeTermIR at hrun
  split at hrun
  next hzero => simp at hrun
  next hpositive =>
    dsimp only at hrun
    split at hrun
    next hvalueZero => simp at hrun
    next hvalueNonzero =>
      have hvalueNat :=
        CLPoly.Impl.StrictWordArithmetic.nmod_mul_ir_correct_of_configured
          this term.2.val (term.1.deg.toUInt64 % this._p) hcfg hreduced.2
      have hvalueLt : (derivativeValueIR this term).toNat < this._p.toNat := by
        rw [derivativeValueIR, hvalueNat]
        exact Nat.mod_lt _ (Fact.out : Nat.Prime this._p.toNat).pos
      simp only [Option.some.injEq] at hrun
      subst output
      exact ⟨rfl, hreduced.1, hvalueNonzero, hvalueLt⟩

theorem derivativeIR_canonical (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (source : SparsePolyZp)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hcanonical : SparsePolyZp.Canonical this._p.toNat source) :
    SparsePolyZp.Canonical this._p.toNat (derivativeIR this source) := by
  rw [SparsePolyZp.Canonical, SparsePolyZp.WellFormed_arr,
    derivativeIR_toList]
  refine ⟨?_, ?_, ?_⟩
  · intro output houtput
    rw [List.mem_filterMap] at houtput
    rcases houtput with ⟨term, hterm, hrun⟩
    have hproperties := derivativeTermIR_some_properties this term output hcfg
      (hcanonical.1 term hterm) hrun
    exact ⟨hproperties.2.1, hproperties.2.2.2⟩
  · apply List.isChain_iff_pairwise.mpr
    apply List.Pairwise.filterMap
      (R := fun a b : UMonomial × Zp => a.1.deg > b.1.deg)
      (S := fun a b : UMonomial × Zp => a.1.deg > b.1.deg)
      (derivativeTermIR this)
    · intro left right horder leftOut hleft rightOut hright
      have hleftProperties := derivativeTermIR_some_degree this left leftOut
        hleft
      have hrightProperties := derivativeTermIR_some_degree this right
        rightOut hright
      rw [hleftProperties.1, hrightProperties.1]
      omega
    · exact List.isChain_iff_pairwise.mp hcanonical.2.1
  · intro output houtput
    rw [List.mem_filterMap] at houtput
    rcases houtput with ⟨term, hterm, hrun⟩
    exact (derivativeTermIR_some_properties this term output hcfg
      (hcanonical.1 term hterm) hrun).2.2.1

theorem derivativeIR_eq_empty_iff (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (source : SparsePolyZp)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hcanonical : SparsePolyZp.Canonical this._p.toNat source)
    (hdegree : ∀ term ∈ source.toList, term.1.deg < UInt64.size) :
    derivativeIR this source = #[] ↔
      Polynomial.derivative (SparsePolyZp.toPoly this._p.toNat source) = 0 := by
  constructor
  · intro hempty
    rw [← derivativeIR_toPoly this source hcfg hcanonical hdegree, hempty]
    simp
  · intro hzero
    apply SparsePolyZp.toPoly_inj_canonical this._p.toNat
    · exact derivativeIR_canonical this source hcfg hcanonical
    · refine ⟨?_, ?_, ?_⟩
      · intro term hterm
        simp at hterm
      · simp
      · simp
    · rw [derivativeIR_toPoly this source hcfg hcanonical hdegree, hzero]
      simp

theorem derivativeIR_isEmpty_iff (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (source : SparsePolyZp)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hcanonical : SparsePolyZp.Canonical this._p.toNat source)
    (hdegree : ∀ term ∈ source.toList, term.1.deg < UInt64.size) :
    (derivativeIR this source).isEmpty = true ↔
      Polynomial.derivative (SparsePolyZp.toPoly this._p.toNat source) = 0 := by
  rw [Array.isEmpty_iff]
  exact derivativeIR_eq_empty_iff this source hcfg hcanonical hdegree

/-- Concrete quotient cell produced by the C++ `v2_.size()==1` branch of
`pair_vec_div`. -/
def pairVecDivSingleTermIR (this : DenseUPolyZp)
    (divisor term : UMonomial × Zp) : Option (UMonomial × Zp) :=
  if divisor.1.deg ≤ term.1.deg then
    let inverse := Generated.StrictGCD.dense_upoly_zp_nmod_inv_ir this
      divisor.2.val
    let value := Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this
      term.2.val inverse
    some (⟨term.1.deg - divisor.1.deg⟩, ⟨value, term.2.prime⟩)
  else
    none

def pairVecDivSingleLoopIR (this : DenseUPolyZp) (index : Nat)
    (output dividend : SparsePolyZp) (divisor : UMonomial × Zp) :
    SparsePolyZp :=
  if h : index < dividend.size then
    let next := match pairVecDivSingleTermIR this divisor dividend[index] with
      | none => output
      | some term => output.push term
    pairVecDivSingleLoopIR this (index + 1) next dividend divisor
  else
    output
termination_by dividend.size - index
decreasing_by omega

theorem pairVecDivSingleLoopIR_toList (this : DenseUPolyZp) (index : Nat)
    (output dividend : SparsePolyZp) (divisor : UMonomial × Zp) :
    (pairVecDivSingleLoopIR this index output dividend divisor).toList =
      output.toList ++ (dividend.toList.drop index).filterMap
        (pairVecDivSingleTermIR this divisor) := by
  rw [pairVecDivSingleLoopIR]
  split
  next hmore =>
    rw [pairVecDivSingleLoopIR_toList]
    have hdrop := List.drop_eq_getElem_cons
      (l := dividend.toList) (i := index) (by simpa using hmore)
    rw [hdrop, List.filterMap_cons, Array.getElem_toList hmore]
    cases hterm : pairVecDivSingleTermIR this divisor dividend[index] with
    | none => simp
    | some term => simp [List.append_assoc]
  next hdone =>
    have hindex : dividend.toList.length ≤ index := by simpa using hdone
    rw [List.drop_eq_nil_iff.mpr hindex]
    simp
termination_by dividend.size - index
decreasing_by omega

theorem pairVecDivSingleTermIR_refines (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (divisor term output : UMonomial × Zp)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hdivisorReduced : Zp.Reduced this._p.toNat divisor.2)
    (htermReduced : Zp.Reduced this._p.toNat term.2)
    (hdivisorNonzero : divisor.2.val ≠ 0)
    (htermNonzero : term.2.val ≠ 0)
    (hrun : pairVecDivSingleTermIR this divisor term = some output) :
    output.1.deg = term.1.deg - divisor.1.deg ∧
      Zp.Reduced this._p.toNat output.2 ∧ output.2.val ≠ 0 ∧
      Zp.toZMod this._p.toNat output.2 =
        Zp.toZMod this._p.toNat term.2 /
          Zp.toZMod this._p.toNat divisor.2 := by
  unfold pairVecDivSingleTermIR at hrun
  split at hrun
  next hdegree =>
    dsimp only at hrun
    simp only [Option.some.injEq] at hrun
    subst output
    let inverse := Generated.StrictGCD.dense_upoly_zp_nmod_inv_ir this
      divisor.2.val
    let value := Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this
      term.2.val inverse
    have hdivisorPos : 0 < divisor.2.val.toNat := by
      exact Nat.pos_of_ne_zero (fun hzero => hdivisorNonzero
        (UInt64.toNat_inj.mp (by simpa using hzero)))
    have hquotient :=
      CLPoly.Impl.StrictWordArithmetic.quotientCoeff_eliminates_lead this
        term.2.val divisor.2.val hcfg
        (Fact.out : Nat.Prime this._p.toNat) htermReduced.2 hdivisorPos
        hdivisorReduced.2
    change value.toNat < this._p.toNat ∧
      (value.toNat * divisor.2.val.toNat) % this._p.toNat =
        term.2.val.toNat at hquotient
    have hfieldMul : (value.toNat : ZMod this._p.toNat) *
        (divisor.2.val.toNat : ZMod this._p.toNat) =
          (term.2.val.toNat : ZMod this._p.toNat) := by
      rw [← Nat.cast_mul, ← ZMod.natCast_mod]
      exact congrArg (fun n : Nat => (n : ZMod this._p.toNat)) hquotient.2
    have hdivisorFieldNonzero : Zp.toZMod this._p.toNat divisor.2 ≠ 0 :=
      Zp.toZMod_ne_zero_of_val_ne_zero this._p.toNat divisor.2
        hdivisorReduced hdivisorNonzero
    have hvalueField : (value.toNat : ZMod this._p.toNat) =
        Zp.toZMod this._p.toNat term.2 /
          Zp.toZMod this._p.toNat divisor.2 := by
      apply (eq_div_iff hdivisorFieldNonzero).2
      unfold Zp.toZMod
      exact hfieldMul
    have hvalueNonzero : value ≠ 0 := by
      intro hzero
      have htermFieldNonzero := Zp.toZMod_ne_zero_of_val_ne_zero
        this._p.toNat term.2 htermReduced htermNonzero
      have hquotientFieldNonzero := div_ne_zero htermFieldNonzero
        hdivisorFieldNonzero
      apply hquotientFieldNonzero
      rw [← hvalueField, hzero]
      simp
    refine ⟨rfl, ⟨htermReduced.1, hquotient.1⟩, ?_, ?_⟩
    · exact hvalueNonzero
    · exact hvalueField
  next hdegree => simp at hrun

theorem pairVecDivSingleTermIR_some_degree (this : DenseUPolyZp)
    (divisor term output : UMonomial × Zp)
    (hrun : pairVecDivSingleTermIR this divisor term = some output) :
    output.1.deg = term.1.deg - divisor.1.deg ∧
      divisor.1.deg ≤ term.1.deg := by
  unfold pairVecDivSingleTermIR at hrun
  split at hrun
  next hdegree =>
    dsimp only at hrun
    simp only [Option.some.injEq] at hrun
    subst output
    exact ⟨rfl, hdegree⟩
  next hdegree => simp at hrun

theorem pairVecDivSingleLoopIR_zero_canonical (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (dividend : SparsePolyZp) (divisor : UMonomial × Zp)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hdividendCanonical : SparsePolyZp.Canonical this._p.toNat dividend)
    (hdivisorReduced : Zp.Reduced this._p.toNat divisor.2)
    (hdivisorNonzero : divisor.2.val ≠ 0) :
    SparsePolyZp.Canonical this._p.toNat
      (pairVecDivSingleLoopIR this 0 #[] dividend divisor) := by
  rw [SparsePolyZp.Canonical, SparsePolyZp.WellFormed_arr,
    pairVecDivSingleLoopIR_toList]
  simp only [List.drop_zero, List.nil_append]
  refine ⟨?_, ?_, ?_⟩
  · intro output houtput
    rw [List.mem_filterMap] at houtput
    rcases houtput with ⟨term, hterm, hrun⟩
    exact (pairVecDivSingleTermIR_refines this divisor term output hcfg
      hdivisorReduced (hdividendCanonical.1 term hterm) hdivisorNonzero
      (hdividendCanonical.2.2 term hterm) hrun).2.1
  · rw [List.isChain_iff_pairwise]
    apply List.Pairwise.filterMap
      (R := fun a b : UMonomial × Zp => a.1.deg > b.1.deg)
      (S := fun a b : UMonomial × Zp => a.1.deg > b.1.deg)
      (pairVecDivSingleTermIR this divisor)
    · intro left right horder leftOut hleft rightOut hright
      have hleftDegree := pairVecDivSingleTermIR_some_degree this divisor
        left leftOut hleft
      have hrightDegree := pairVecDivSingleTermIR_some_degree this divisor
        right rightOut hright
      rw [hleftDegree.1, hrightDegree.1]
      omega
    · exact List.isChain_iff_pairwise.mp hdividendCanonical.2.1
  · intro output houtput
    rw [List.mem_filterMap] at houtput
    rcases houtput with ⟨term, hterm, hrun⟩
    exact (pairVecDivSingleTermIR_refines this divisor term output hcfg
      hdivisorReduced (hdividendCanonical.1 term hterm) hdivisorNonzero
      (hdividendCanonical.2.2 term hterm) hrun).2.2.1

theorem pairVecDivSingleTermIR_monomial_mul (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (divisor term output : UMonomial × Zp)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hdivisorReduced : Zp.Reduced this._p.toNat divisor.2)
    (htermReduced : Zp.Reduced this._p.toNat term.2)
    (hdivisorNonzero : divisor.2.val ≠ 0)
    (htermNonzero : term.2.val ≠ 0)
    (hrun : pairVecDivSingleTermIR this divisor term = some output) :
    Polynomial.monomial output.1.deg
          (Zp.toZMod this._p.toNat output.2) *
        Polynomial.monomial divisor.1.deg
          (Zp.toZMod this._p.toNat divisor.2) =
      Polynomial.monomial term.1.deg
        (Zp.toZMod this._p.toNat term.2) := by
  have hrefines := pairVecDivSingleTermIR_refines this divisor term output
    hcfg hdivisorReduced htermReduced hdivisorNonzero htermNonzero hrun
  have hdivisorFieldNonzero : Zp.toZMod this._p.toNat divisor.2 ≠ 0 :=
    Zp.toZMod_ne_zero_of_val_ne_zero this._p.toNat divisor.2
      hdivisorReduced hdivisorNonzero
  rw [Polynomial.monomial_mul_monomial, hrefines.1,
    Nat.sub_add_cancel (pairVecDivSingleTermIR_some_degree this divisor term
      output hrun).2, hrefines.2.2.2,
    div_mul_cancel₀ _ hdivisorFieldNonzero]

theorem listSum_pairVecDivSingleTermIR_mul (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (divisor : UMonomial × Zp)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hdivisorReduced : Zp.Reduced this._p.toNat divisor.2)
    (hdivisorNonzero : divisor.2.val ≠ 0) :
    ∀ terms : List (UMonomial × Zp),
      (∀ term ∈ terms, Zp.Reduced this._p.toNat term.2) →
      (∀ term ∈ terms, term.2.val ≠ 0) →
      (∀ term ∈ terms, divisor.1.deg ≤ term.1.deg) →
      listSum this._p.toNat
          (terms.filterMap (pairVecDivSingleTermIR this divisor)) *
        Polynomial.monomial divisor.1.deg
          (Zp.toZMod this._p.toNat divisor.2) =
      listSum this._p.toNat terms := by
  intro terms hreduced hnonzero hdegree
  induction terms with
  | nil => simp [listSum]
  | cons term rest ih =>
      have htermReduced := hreduced term List.mem_cons_self
      have htermNonzero := hnonzero term List.mem_cons_self
      have htermDegree := hdegree term List.mem_cons_self
      have hrestReduced : ∀ item ∈ rest,
          Zp.Reduced this._p.toNat item.2 := by
        intro item hitem
        exact hreduced item (List.mem_cons_of_mem term hitem)
      have hrestNonzero : ∀ item ∈ rest, item.2.val ≠ 0 := by
        intro item hitem
        exact hnonzero item (List.mem_cons_of_mem term hitem)
      have hrestDegree : ∀ item ∈ rest, divisor.1.deg ≤ item.1.deg := by
        intro item hitem
        exact hdegree item (List.mem_cons_of_mem term hitem)
      let inverse := Generated.StrictGCD.dense_upoly_zp_nmod_inv_ir this
        divisor.2.val
      let value := Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this
        term.2.val inverse
      let output : UMonomial × Zp :=
        (⟨term.1.deg - divisor.1.deg⟩, ⟨value, term.2.prime⟩)
      have hrun : pairVecDivSingleTermIR this divisor term = some output := by
        simp [pairVecDivSingleTermIR, htermDegree, output, value, inverse]
      rw [List.filterMap_cons, hrun, listSum_cons, add_mul,
        pairVecDivSingleTermIR_monomial_mul this divisor term output hcfg
          hdivisorReduced htermReduced hdivisorNonzero htermNonzero hrun,
        ih hrestReduced hrestNonzero hrestDegree, listSum_cons]

def pairVecDivSingleBranchIR (this : DenseUPolyZp)
    (dividend divisor : SparsePolyZp) : RawExec SparsePolyZp :=
  if hone : divisor.size = 1 then
    .ok (pairVecDivSingleLoopIR this 0 #[] dividend divisor[0])
  else
    .error .assertionFailure

theorem pairVecDivSingleBranchIR_refines (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (dividend divisor : SparsePolyZp)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hdividendCanonical : SparsePolyZp.Canonical this._p.toNat dividend)
    (hdivisorCanonical : SparsePolyZp.Canonical this._p.toNat divisor)
    (hdivisorSize : divisor.size = 1)
    (hdegree : ∀ term ∈ dividend.toList,
      divisor[0].1.deg ≤ term.1.deg) :
    ∃ quotient,
      pairVecDivSingleBranchIR this dividend divisor = .ok quotient ∧
      SparsePolyZp.Canonical this._p.toNat quotient ∧
      SparsePolyZp.toPoly this._p.toNat quotient *
          SparsePolyZp.toPoly this._p.toNat divisor =
        SparsePolyZp.toPoly this._p.toNat dividend := by
  rcases Array.size_eq_one_iff.mp hdivisorSize with ⟨divisorTerm, rfl⟩
  have hdivisorMem : divisorTerm ∈ (#[divisorTerm] : SparsePolyZp).toList := by
    simp
  have hdivisorReduced := hdivisorCanonical.1 divisorTerm hdivisorMem
  have hdivisorNonzero := hdivisorCanonical.2.2 divisorTerm hdivisorMem
  let quotient := pairVecDivSingleLoopIR this 0 #[] dividend divisorTerm
  refine ⟨quotient, by simp [pairVecDivSingleBranchIR, quotient], ?_, ?_⟩
  · exact pairVecDivSingleLoopIR_zero_canonical this dividend divisorTerm hcfg
      hdividendCanonical hdivisorReduced hdivisorNonzero
  · unfold SparsePolyZp.toPoly quotient
    rw [pairVecDivSingleLoopIR_toList]
    simp only [List.drop_zero, List.nil_append]
    have hdivisorPoly : listSum this._p.toNat [divisorTerm] =
        Polynomial.monomial divisorTerm.1.deg
          (Zp.toZMod this._p.toNat divisorTerm.2) := by
      simp [listSum]
    rw [hdivisorPoly]
    rw [listSum_pairVecDivSingleTermIR_mul this divisorTerm hcfg
      hdivisorReduced hdivisorNonzero dividend.toList hdividendCanonical.1
      hdividendCanonical.2.2 (by simpa using hdegree)]

theorem eq_divByMonic_of_mul_eq {p : Nat}
    [Fact (Nat.Prime p)]
    (quotient dividend divisor : Polynomial (ZMod p))
    (hdivisorMonic : divisor.Monic) (hdivisorDvd : divisor ∣ dividend)
    (hmul : quotient * divisor = dividend) :
    quotient = dividend /ₘ divisor := by
  have hmodZero : dividend %ₘ divisor = 0 :=
    (Polynomial.modByMonic_eq_zero_iff_dvd hdivisorMonic).2 hdivisorDvd
  have hstandard := Polynomial.modByMonic_add_div dividend divisor
  rw [hmodZero, zero_add] at hstandard
  apply mul_left_cancel₀ hdivisorMonic.ne_zero
  rw [mul_comm divisor quotient, hmul]
  exact hstandard.symm

theorem pairVecDivSingleBranchIR_refines_divByMonic
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (dividend divisor : SparsePolyZp)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hdividendCanonical : SparsePolyZp.Canonical this._p.toNat dividend)
    (hdivisorCanonical : SparsePolyZp.Canonical this._p.toNat divisor)
    (hdivisorSize : divisor.size = 1)
    (hdegree : ∀ term ∈ dividend.toList,
      divisor[0].1.deg ≤ term.1.deg)
    (hdivisorMonic : (SparsePolyZp.toPoly this._p.toNat divisor).Monic)
    (hdivisorDvd : SparsePolyZp.toPoly this._p.toNat divisor ∣
      SparsePolyZp.toPoly this._p.toNat dividend) :
    ∃ quotient,
      pairVecDivSingleBranchIR this dividend divisor = .ok quotient ∧
      SparsePolyZp.Canonical this._p.toNat quotient ∧
      SparsePolyZp.toPoly this._p.toNat quotient =
        SparsePolyZp.toPoly this._p.toNat dividend /ₘ
          SparsePolyZp.toPoly this._p.toNat divisor := by
  rcases pairVecDivSingleBranchIR_refines this dividend divisor hcfg
      hdividendCanonical hdivisorCanonical hdivisorSize hdegree with
    ⟨quotient, hrun, hcanonical, hmul⟩
  exact ⟨quotient, hrun, hcanonical,
    eq_divByMonic_of_mul_eq _ _ _ hdivisorMonic hdivisorDvd hmul⟩

/-! ### General `pair_vec_div` priority-heap path

The source allocates one `VHC` node for every non-leading divisor term.  A
node is not assigned a monomial until `reset_h` activates it, so reading that
field earlier would be undefined C++.  The safe execution records this fact
with `Option`; an active node always has `some mono`.  Array indices replace
the three source pointer fields without changing their observable meaning.
-/

/-- Safe representation of the source
`VHC<umonomial, size_t, divisor_iterator>` allocation. -/
structure PairVecDivVHCNode where
  mono : Option UMonomial
  quotientIndex : Nat
  divisorIndex : Nat
  next : Option Nat
  deriving Repr

/-- Exact initial value of node `i`: `v1_ptr = 0`, while `v2_ptr` points to
the `(i + 1)`th divisor cell.  `mono` and `next` are not observable before the
node is activated by the `reset_h` loop. -/
def pairVecDivVHCInitialNode (i : Nat) : PairVecDivVHCNode where
  mono := none
  quotientIndex := 0
  divisorIndex := i + 1
  next := none

/-- The source node-allocation loop, with `divisor.size - 1` nodes. -/
def pairVecDivVHCInitLoop (i count : Nat)
    (nodes : Array PairVecDivVHCNode) : Array PairVecDivVHCNode :=
  if h : i < count then
    pairVecDivVHCInitLoop (i + 1) count
      (nodes.push (pairVecDivVHCInitialNode i))
  else
    nodes
termination_by count - i
decreasing_by omega

def pairVecDivVHCInit (divisor : SparsePolyZp) :
    Array PairVecDivVHCNode :=
  pairVecDivVHCInitLoop 0 (divisor.size - 1) #[]

theorem pairVecDivVHCInitLoop_toList (i count : Nat)
    (nodes : Array PairVecDivVHCNode) :
    (pairVecDivVHCInitLoop i count nodes).toList =
      nodes.toList ++
        (List.range' i (count - i)).map pairVecDivVHCInitialNode := by
  rw [pairVecDivVHCInitLoop]
  split
  next hmore =>
    rw [pairVecDivVHCInitLoop_toList]
    rw [Array.toList_push, List.append_assoc]
    congr 1
    rw [show count - i = (count - (i + 1)) + 1 by omega,
      List.range'_succ]
    simp
  next hdone =>
    have hzero : count - i = 0 := Nat.sub_eq_zero_of_le (by omega)
    simp [hzero]
termination_by count - i
decreasing_by omega

theorem pairVecDivVHCInit_toList (divisor : SparsePolyZp) :
    (pairVecDivVHCInit divisor).toList =
      (List.range (divisor.size - 1)).map pairVecDivVHCInitialNode := by
  simp [pairVecDivVHCInit, pairVecDivVHCInitLoop_toList, List.range_eq_range']

theorem pairVecDivVHCInit_size (divisor : SparsePolyZp) :
    (pairVecDivVHCInit divisor).size = divisor.size - 1 := by
  simpa using congrArg List.length (pairVecDivVHCInit_toList divisor)

/-- One execution of the source `reset_h` activation body.  The checked
indices are precisely the three pointer/array dereferences in that body. -/
def pairVecDivVHCActivate (nodeIndex : Nat)
    (nodes : Array PairVecDivVHCNode) (quotient divisor : SparsePolyZp) :
    RawExec (Array PairVecDivVHCNode) :=
  if hn : nodeIndex < nodes.size then
    let node := nodes[nodeIndex]
    if hq : node.quotientIndex < quotient.size then
      if hd : node.divisorIndex < divisor.size then
        let updated := { node with
          mono := some ⟨quotient[node.quotientIndex].1.deg +
            divisor[node.divisorIndex].1.deg⟩
          next := none }
        .ok (nodes.set nodeIndex updated)
      else
        .error .assertionFailure
    else
      .error .assertionFailure
  else
    .error .assertionFailure

theorem pairVecDivVHCActivate_size (nodeIndex : Nat)
    (nodes nodes' : Array PairVecDivVHCNode)
    (quotient divisor : SparsePolyZp)
    (hrun : pairVecDivVHCActivate nodeIndex nodes quotient divisor =
      .ok nodes') :
    nodes'.size = nodes.size := by
  unfold pairVecDivVHCActivate at hrun
  split at hrun <;> try contradiction
  next hn =>
    dsimp only at hrun
    split at hrun <;> try contradiction
    next hq =>
      split at hrun <;> try contradiction
      next hd =>
        simp only [Except.ok.injEq] at hrun
        subst nodes'
        simp

theorem pairVecDivVHCActivate_get (nodeIndex : Nat)
    (nodes nodes' : Array PairVecDivVHCNode)
    (quotient divisor : SparsePolyZp)
    (hrun : pairVecDivVHCActivate nodeIndex nodes quotient divisor =
      .ok nodes') :
    ∃ hn : nodeIndex < nodes.size,
      ∃ hq : nodes[nodeIndex].quotientIndex < quotient.size,
      ∃ hd : nodes[nodeIndex].divisorIndex < divisor.size,
        nodes'[nodeIndex]? = some { nodes[nodeIndex] with
            mono := some ⟨quotient[nodes[nodeIndex].quotientIndex].1.deg +
              divisor[nodes[nodeIndex].divisorIndex].1.deg⟩
            next := none } := by
  unfold pairVecDivVHCActivate at hrun
  split at hrun
  next hn =>
    dsimp only at hrun
    split at hrun
    next hq =>
      split at hrun
      next hd =>
        simp only [Except.ok.injEq] at hrun
        subst nodes'
        refine ⟨hn, hq, hd, ?_⟩
        simp
      next hd => contradiction
    next hq => contradiction
  next hn => contradiction

/-- Checked write of the source `VHC::next` pointer field. -/
def pairVecDivVHCSetNext (nodeIndex : Nat) (next : Option Nat)
    (nodes : Array PairVecDivVHCNode) :
    RawExec (Array PairVecDivVHCNode) :=
  if hn : nodeIndex < nodes.size then
    .ok (nodes.set nodeIndex { nodes[nodeIndex] with next := next })
  else
    .error .assertionFailure

/-- Checked active monomial read.  A `none` monomial is exactly an attempted
read of a source node before its first `reset_h` activation. -/
def pairVecDivVHCMono (nodeIndex : Nat)
    (nodes : Array PairVecDivVHCNode) : RawExec UMonomial :=
  if hn : nodeIndex < nodes.size then
    match nodes[nodeIndex].mono with
    | some mono => .ok mono
    | none => .error .assertionFailure
  else
    .error .assertionFailure

/-- Source parent-slot calculation `(i - 1) >> 1`. -/
def pairVecDivVHCParent (i : Nat) : Nat := (i - 1) / 2

/-- The pointer-copying `for` body shared by the two non-bucket insertion
branches.  `stop` is the source slot `0` or `i1`; the newly appended final
slot makes every write valid. -/
def pairVecDivVHCBubble (i stop newNode : Nat) (heap : Array Nat) :
    RawExec (Array Nat) :=
  if hi : i < heap.size then
    if hstop : stop ≤ i then
      if heq : i = stop then
        .ok (heap.set i newNode)
      else
        let parent := pairVecDivVHCParent i
        if hp : parent < heap.size then
          pairVecDivVHCBubble parent stop newNode
            (heap.set i heap[parent])
        else
          .error .assertionFailure
    else
      .error .assertionFailure
  else
    .error .assertionFailure
termination_by i
decreasing_by
  have hipos : 0 < i := by omega
  unfold pairVecDivVHCParent
  have hhalf : (i - 1) / 2 ≤ i - 1 := Nat.div_le_self _ _
  omega

theorem pairVecDivVHCBubble_size (i stop newNode : Nat)
    (heap heap' : Array Nat)
    (hrun : pairVecDivVHCBubble i stop newNode heap = .ok heap') :
    heap'.size = heap.size := by
  rw [pairVecDivVHCBubble] at hrun
  split at hrun <;> try contradiction
  next hi =>
    split at hrun <;> try contradiction
    next hstop =>
      split at hrun
      next heq =>
        simp only [Except.ok.injEq] at hrun
        subst heap'
        simp
      next heq =>
        dsimp only at hrun
        split at hrun <;> try contradiction
        next hp =>
          have hrec := pairVecDivVHCBubble_size
            (pairVecDivVHCParent i) stop newNode
            (heap.set i heap[pairVecDivVHCParent i]) heap' hrun
          simpa using hrec
termination_by i
decreasing_by
  have hipos : 0 < i := by omega
  unfold pairVecDivVHCParent
  have hhalf : (i - 1) / 2 ≤ i - 1 := Nat.div_le_self _ _
  omega

/-- The source search
`for (i1 = (heap_size-1)>>1; comp(new, heap[i1]); i1=(i1-1)>>1)`.
The root comparison performed by the caller makes reaching `i1 = 0` with a
strict comparison an invariant failure rather than unsigned wraparound. -/
def pairVecDivVHCFindAnchor (newDegree i : Nat) (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode) : RawExec Nat :=
  if hi : i < heap.size then
    match pairVecDivVHCMono heap[i] nodes with
    | .error fault => .error fault
    | .ok mono =>
        if newDegree > mono.deg then
          if hzero : i = 0 then
            .error .assertionFailure
          else
            pairVecDivVHCFindAnchor newDegree (pairVecDivVHCParent i)
              heap nodes
        else
          .ok i
  else
    .error .assertionFailure
termination_by i
decreasing_by
  have hipos : 0 < i := Nat.pos_of_ne_zero hzero
  unfold pairVecDivVHCParent
  have hhalf : (i - 1) / 2 ≤ i - 1 := Nat.div_le_self _ _
  omega

/-- Internal insertion shift.  Unlike root insertion, the source stops when
the *parent* is `anchor`, then writes the new node into the current child. -/
def pairVecDivVHCBubbleBelow (i anchor newNode : Nat)
    (heap : Array Nat) : RawExec (Array Nat) :=
  if hi : i < heap.size then
    if hpos : 0 < i then
      let parent := pairVecDivVHCParent i
      if heq : parent = anchor then
        .ok (heap.set i newNode)
      else if hp : parent < heap.size then
        pairVecDivVHCBubbleBelow parent anchor newNode
          (heap.set i heap[parent])
      else
        .error .assertionFailure
    else
      .error .assertionFailure
  else
    .error .assertionFailure
termination_by i
decreasing_by
  unfold pairVecDivVHCParent
  have hhalf : (i - 1) / 2 ≤ i - 1 := Nat.div_le_self _ _
  omega

/-- Exact checked execution of current C++ `VHC_insert`.  Heap cells contain
indices into the separately allocated node array; equal-degree entries are
linked through `next`, exactly as in the source. -/
def pairVecDivVHCInsert (newNode : Nat) (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode) :
    RawExec (Array Nat × Array PairVecDivVHCNode) := do
  let newMono ← pairVecDivVHCMono newNode nodes
  if hempty : heap.size = 0 then
    let nodes' ← pairVecDivVHCSetNext newNode none nodes
    .ok (#[newNode], nodes')
  else
    let rootNode := heap[0]'(Nat.pos_of_ne_zero hempty)
    let rootMono ← pairVecDivVHCMono rootNode nodes
    if hequal : newMono.deg = rootMono.deg then
      let nodes' ← pairVecDivVHCSetNext newNode (some rootNode) nodes
      .ok (heap.set 0 newNode, nodes')
    else if hgreater : newMono.deg > rootMono.deg then
      let nodes' ← pairVecDivVHCSetNext newNode none nodes
      let heap' ← pairVecDivVHCBubble heap.size 0 newNode
        (heap.push newNode)
      .ok (heap', nodes')
    else
      let firstAnchor := pairVecDivVHCParent heap.size
      let anchor ← pairVecDivVHCFindAnchor newMono.deg firstAnchor heap nodes
      if ha : anchor < heap.size then
        let anchorNode := heap[anchor]
        let anchorMono ← pairVecDivVHCMono anchorNode nodes
        if hequalAnchor : newMono.deg = anchorMono.deg then
          let nodes' ← pairVecDivVHCSetNext newNode (some anchorNode) nodes
          .ok (heap.set anchor newNode, nodes')
        else
          let nodes' ← pairVecDivVHCSetNext newNode none nodes
          let heap' ← pairVecDivVHCBubbleBelow heap.size anchor newNode
            (heap.push newNode)
          .ok (heap', nodes')
      else
        .error .assertionFailure

/-- Downward pointer-copy loop of `VHC_extract`.  `limit` is the decremented
source `heap_size` (`s`), and `lastNode = heap[s]` remains readable as the
sentinel until the loop finishes. -/
def pairVecDivVHCSiftDown (i child limit lastNode : Nat)
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode) :
    RawExec (Array Nat) :=
  if hi : i < heap.size then
    if hlimit : limit < heap.size then
      if hchild : child < limit then
        have hright : child + 1 < heap.size := by omega
        let leftNode := heap[child]
        let rightNode := heap[child + 1]
        match pairVecDivVHCMono leftNode nodes,
            pairVecDivVHCMono rightNode nodes,
            pairVecDivVHCMono lastNode nodes with
        | .ok leftMono, .ok rightMono, .ok lastMono =>
            let selected := if leftMono.deg > rightMono.deg then
              child else child + 1
            have hselected : selected < heap.size := by
              dsimp only [selected]
              split <;> omega
            let selectedNode := heap[selected]
            let selectedMono := if leftMono.deg > rightMono.deg then
              leftMono else rightMono
            if hgreater : selectedMono.deg > lastMono.deg then
              pairVecDivVHCSiftDown selected (selected * 2 + 1) limit
                lastNode (heap.set i selectedNode) nodes
            else
              .ok (heap.set i lastNode)
        | .error fault, _, _ => .error fault
        | _, .error fault, _ => .error fault
        | _, _, .error fault => .error fault
      else
        .ok (heap.set i lastNode)
    else
      .error .assertionFailure
  else
    .error .assertionFailure
termination_by limit - child
decreasing_by
  have hselected : child ≤ selected := by
    simp only [selected]
    split <;> omega
  omega

/-- Exact checked execution of current C++ `VHC_extract`.  The returned array
is the active prefix after `--heap_size`; storage beyond that logical size is
unobservable and is removed with `pop`. -/
def pairVecDivVHCExtract (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode) : RawExec (Array Nat) :=
  if hnonempty : 0 < heap.size then
    let limit := heap.size - 1
    let lastNode := heap[limit]
    match pairVecDivVHCSiftDown 0 1 limit lastNode heap nodes with
    | .ok shifted => .ok shifted.pop
    | .error fault => .error fault
  else
    .error .assertionFailure

theorem pairVecDivVHCSiftDown_size (i child limit lastNode : Nat)
    (heap heap' : Array Nat) (nodes : Array PairVecDivVHCNode)
    (hrun : pairVecDivVHCSiftDown i child limit lastNode heap nodes =
      .ok heap') :
    heap'.size = heap.size := by
  rw [pairVecDivVHCSiftDown] at hrun
  split at hrun <;> try contradiction
  next hi =>
    split at hrun <;> try contradiction
    next hlimit =>
      split at hrun
      next hchild =>
        dsimp only at hrun
        cases hleft : pairVecDivVHCMono heap[child] nodes with
        | error fault => simp [hleft] at hrun
        | ok leftMono =>
          cases hright : pairVecDivVHCMono heap[child + 1] nodes with
          | error fault => simp [hleft, hright] at hrun
          | ok rightMono =>
            cases hlast : pairVecDivVHCMono lastNode nodes with
            | error fault => simp [hleft, hright, hlast] at hrun
            | ok lastMono =>
              simp only [hleft, hright, hlast] at hrun
              by_cases hselected : leftMono.deg > rightMono.deg
              · simp only [hselected, ↓reduceIte] at hrun
                split at hrun
                next hgreater =>
                  have hrec := pairVecDivVHCSiftDown_size
                    child (child * 2 + 1) limit lastNode
                    (heap.set i heap[child]) heap' nodes hrun
                  simpa using hrec
                next hgreater =>
                  simp only [Except.ok.injEq] at hrun
                  subst heap'
                  simp
              · simp only [hselected, ↓reduceIte] at hrun
                split at hrun
                next hgreater =>
                  have hrec := pairVecDivVHCSiftDown_size
                    (child + 1) ((child + 1) * 2 + 1) limit lastNode
                    (heap.set i heap[child + 1]) heap' nodes hrun
                  simpa using hrec
                next hgreater =>
                  simp only [Except.ok.injEq] at hrun
                  subst heap'
                  simp
      next hchild =>
        simp only [Except.ok.injEq] at hrun
        subst heap'
        simp
termination_by limit - child
decreasing_by
  all_goals omega

theorem pairVecDivVHCExtract_size (heap heap' : Array Nat)
    (nodes : Array PairVecDivVHCNode)
    (hrun : pairVecDivVHCExtract heap nodes = .ok heap') :
    heap'.size + 1 = heap.size := by
  unfold pairVecDivVHCExtract at hrun
  split at hrun <;> try contradiction
  next hnonempty =>
    dsimp only at hrun
    split at hrun <;> try contradiction
    next shifted hsift =>
      simp only [Except.ok.injEq] at hrun
      subst heap'
      have hsize := pairVecDivVHCSiftDown_size 0 1 (heap.size - 1)
        heap[heap.size - 1] heap shifted nodes hsift
      rw [Array.size_pop, hsize]
      omega

/-- Proof-carrying safe boundary for `VHC_extract`; it executes the same raw
definition and packages its established logical-size decrement. -/
def pairVecDivVHCExtractChecked (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode) :
    RawExec { heap' : Array Nat // heap'.size + 1 = heap.size } :=
  match hrun : pairVecDivVHCExtract heap nodes with
  | .error fault => .error fault
  | .ok heap' => .ok ⟨heap', pairVecDivVHCExtract_size heap heap' nodes hrun⟩

/-- Every active heap node denotes the concrete product cell addressed by its
two source pointer fields. -/
def PairVecDivVHCNodeDenotes (quotient divisor : SparsePolyZp)
    (node : PairVecDivVHCNode) : Prop :=
  ∃ quotientTerm divisorTerm,
    quotient[node.quotientIndex]? = some quotientTerm ∧
    divisor[node.divisorIndex]? = some divisorTerm ∧
    node.mono = some ⟨quotientTerm.1.deg + divisorTerm.1.deg⟩

/-- All active heap slots are valid initialized node pointers. -/
def PairVecDivVHCHeapPointersValid (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode) : Prop :=
  ∀ slot : Nat, slot < heap.size →
    ∃ nodeIndex node mono,
      heap[slot]? = some nodeIndex ∧ nodes[nodeIndex]? = some node ∧
      node.mono = some mono

/-- Binary max-heap order used by the source comparator. -/
def PairVecDivVHCHeapOrdered (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode) : Prop :=
  ∀ child parent : Nat, child < heap.size →
    parent = pairVecDivVHCParent child → 0 < child →
    ∀ childNode parentNode childMono parentMono,
        heap[child]? = some childNode →
        heap[parent]? = some parentNode →
        (nodes[childNode]?.map PairVecDivVHCNode.mono).join =
          some childMono →
        (nodes[parentNode]?.map PairVecDivVHCNode.mono).join =
          some parentMono →
        childMono.deg ≤ parentMono.deg

/-- A `next` edge is a valid equal-monomial bucket link. -/
def PairVecDivVHCNextValid (nodes : Array PairVecDivVHCNode) : Prop :=
  ∀ (i : Nat) (node : PairVecDivVHCNode) (next : Nat),
    nodes[i]? = some node → node.next = some next →
    ∃ nextNode, nodes[next]? = some nextNode ∧ node.mono = nextNode.mono

/-- Main representation invariant for the source general-division heap. -/
def PairVecDivVHCInvariant (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode)
    (quotient divisor : SparsePolyZp) : Prop :=
  nodes.size = divisor.size - 1 ∧
    PairVecDivVHCHeapPointersValid heap nodes ∧
    PairVecDivVHCHeapOrdered heap nodes ∧
    PairVecDivVHCNextValid nodes ∧
    (∀ (i : Nat) (node : PairVecDivVHCNode),
      nodes[i]? = some node → node.mono ≠ none →
      PairVecDivVHCNodeDenotes quotient divisor node)

/-- Generated coefficient execution of the source
`submul(k, new_v[q].second, divisor[d].second)`. -/
def pairVecDivSubmulIR (this : DenseUPolyZp) (k a b : UInt64) : UInt64 :=
  let product := Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this a b
  Generated.StrictPolyAddSub.dense_upoly_zp_nmod_sub_ir this k product

theorem pairVecDivSubmulIR_toZMod (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)] (k a b : UInt64)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hk : k.toNat < this._p.toNat) (ha : a.toNat < this._p.toNat) :
    ((pairVecDivSubmulIR this k a b).toNat : ZMod this._p.toNat) =
      (k.toNat : ZMod this._p.toNat) -
        (a.toNat : ZMod this._p.toNat) *
          (b.toNat : ZMod this._p.toNat) := by
  let product := Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this a b
  have hpPos : 0 < this._p.toNat := (Fact.out : Nat.Prime this._p.toNat).pos
  have hpWord : this._p ≠ 0 := by
    intro hp
    have : this._p.toNat = 0 := congrArg UInt64.toNat hp
    omega
  have hproductNat :=
    CLPoly.Impl.StrictWordArithmetic.nmod_mul_ir_correct_of_configured
      this a b hcfg ha
  change product.toNat = (a.toNat * b.toNat) % this._p.toNat at hproductNat
  have hproductLtNat : product.toNat < this._p.toNat := by
    rw [hproductNat]
    exact Nat.mod_lt _ hpPos
  have hkWord : k < this._p := by
    simpa [UInt64.lt_iff_toNat_lt] using hk
  have hproductWord : product < this._p := by
    simpa [UInt64.lt_iff_toNat_lt] using hproductLtNat
  have hsub := CLPoly.Impl.StrictPolyAddSubRefinement.nmodSub_cast
    this k product hpWord hkWord hproductWord
  change ((pairVecDivSubmulIR this k a b).toNat : ZMod this._p.toNat) = _
  rw [show pairVecDivSubmulIR this k a b =
      Generated.StrictPolyAddSub.dense_upoly_zp_nmod_sub_ir this k product
    by rfl, hsub]
  congr 1
  rw [hproductNat, ZMod.natCast_mod]
  push_cast
  rfl

/-- One pass through a non-null node in the source inner `while (heap[0])`
bucket chain.  The heap root itself is represented by the returned `next`:
the C++ code temporarily allows it to become null before `VHC_extract`.
Nodes that advance are appended to `lin`; exhausted nodes increment
`reset_h`. -/
def pairVecDivVHCConsumeNode (this : DenseUPolyZp) (nodeIndex : Nat)
    (k : UInt64) (nodes : Array PairVecDivVHCNode) (lin : Array Nat)
    (resetH : Nat) (quotient divisor : SparsePolyZp) :
    RawExec (UInt64 × Array PairVecDivVHCNode × Array Nat × Nat × Option Nat) :=
  if hn : nodeIndex < nodes.size then
    let node := nodes[nodeIndex]
    if hq : node.quotientIndex < quotient.size then
      if hd : node.divisorIndex < divisor.size then
        let k' := pairVecDivSubmulIR this k
          quotient[node.quotientIndex].2.val divisor[node.divisorIndex].2.val
        let quotientIndex' := node.quotientIndex + 1
        if hadvance : quotientIndex' < quotient.size then
          let node' := { node with
            quotientIndex := quotientIndex'
            mono := some ⟨quotient[quotientIndex'].1.deg +
              divisor[node.divisorIndex].1.deg⟩ }
          .ok (k', nodes.set nodeIndex node', lin.push nodeIndex,
            resetH, node.next)
        else if hexhausted : quotientIndex' = quotient.size then
          .ok (k', nodes, lin, resetH + 1, node.next)
        else
          .error .assertionFailure
      else
        .error .assertionFailure
    else
      .error .assertionFailure
  else
    .error .assertionFailure

theorem pairVecDivVHCConsumeNode_progress (this : DenseUPolyZp)
    (nodeIndex : Nat) (k k' : UInt64)
    (nodes nodes' : Array PairVecDivVHCNode) (lin lin' : Array Nat)
    (resetH resetH' : Nat) (next : Option Nat)
    (quotient divisor : SparsePolyZp)
    (hrun : pairVecDivVHCConsumeNode this nodeIndex k nodes lin resetH
      quotient divisor = .ok (k', nodes', lin', resetH', next)) :
    (lin'.size = lin.size + 1 ∧ resetH' = resetH) ∨
      (lin'.size = lin.size ∧ resetH' = resetH + 1) := by
  unfold pairVecDivVHCConsumeNode at hrun
  split at hrun <;> try contradiction
  next hn =>
    dsimp only at hrun
    split at hrun <;> try contradiction
    next hq =>
      split at hrun <;> try contradiction
      next hd =>
        split at hrun
        next hadvance =>
          simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
          rcases hrun with ⟨rfl, rfl, rfl, rfl, rfl⟩
          left
          simp
        next hadvance =>
          split at hrun <;> try contradiction
          next hexhausted =>
            simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
            rcases hrun with ⟨rfl, rfl, rfl, rfl, rfl⟩
            right
            simp

/-- Result of consuming one complete equal-monomial `next` bucket. -/
structure PairVecDivVHCBucketResult where
  coefficient : UInt64
  nodes : Array PairVecDivVHCNode
  lin : Array Nat
  resetH : Nat

/-- Exact `while (heap[0] != nullptr)` execution.  `unvisited` is the safe
ownership view of the allocated node block.  Removing the current node is the
well-founded measure and detects a malformed cyclic `next` chain; it is not a
step bound and does not truncate any source execution satisfying the heap
invariant. -/
def pairVecDivVHCConsumeChain (this : DenseUPolyZp)
    (current : Option Nat) (unvisited : Finset Nat) (k : UInt64)
    (nodes : Array PairVecDivVHCNode) (lin : Array Nat) (resetH : Nat)
    (quotient divisor : SparsePolyZp) : RawExec PairVecDivVHCBucketResult :=
  match current with
  | none => .ok (PairVecDivVHCBucketResult.mk k nodes lin resetH)
  | some nodeIndex =>
      if hmem : nodeIndex ∈ unvisited then
        match pairVecDivVHCConsumeNode this nodeIndex k nodes lin resetH
            quotient divisor with
        | .error fault => .error fault
        | .ok (k', nodes', lin', resetH', next) =>
            pairVecDivVHCConsumeChain this next (unvisited.erase nodeIndex)
              k' nodes' lin' resetH' quotient divisor
      else
        .error .assertionFailure
termination_by unvisited.card
decreasing_by
  exact Finset.card_erase_lt_of_mem hmem

def pairVecDivVHCConsumeRootBucket (this : DenseUPolyZp)
    (heap : Array Nat) (k : UInt64) (nodes : Array PairVecDivVHCNode)
    (lin : Array Nat) (resetH : Nat) (quotient divisor : SparsePolyZp) :
    RawExec PairVecDivVHCBucketResult :=
  if hheap : 0 < heap.size then
    pairVecDivVHCConsumeChain this (some heap[0])
      (Finset.range nodes.size) k nodes lin resetH quotient divisor
  else
    .error .assertionFailure

theorem pairVecDivVHCConsumeChain_none (this : DenseUPolyZp)
    (unvisited : Finset Nat) (k : UInt64)
    (nodes : Array PairVecDivVHCNode) (lin : Array Nat) (resetH : Nat)
    (quotient divisor : SparsePolyZp) :
    pairVecDivVHCConsumeChain this none unvisited k nodes lin resetH
      quotient divisor =
        .ok (PairVecDivVHCBucketResult.mk k nodes lin resetH) := by
  rw [pairVecDivVHCConsumeChain]

/-- Result after the source has consumed every heap bucket whose root
monomial equals the current outer-loop monomial. -/
structure PairVecDivVHCEqualDegreeResult where
  heap : Array Nat
  coefficient : UInt64
  nodes : Array PairVecDivVHCNode
  lin : Array Nat
  resetH : Nat

/-- Exact source loop
`while (heap_size > 0 && heap[0]->mono == m)`: consume the root `next` chain,
then execute `VHC_extract`, and repeat only while the new root has the same
degree. -/
def pairVecDivVHCConsumeEqualDegree (this : DenseUPolyZp) (degree : Nat)
    (heap : Array Nat) (k : UInt64) (nodes : Array PairVecDivVHCNode)
    (lin : Array Nat) (resetH : Nat) (quotient divisor : SparsePolyZp) :
    RawExec PairVecDivVHCEqualDegreeResult :=
  if hheap : 0 < heap.size then
    match pairVecDivVHCMono heap[0] nodes with
    | .error fault => .error fault
    | .ok rootMono =>
        if hequal : rootMono.deg = degree then
          match pairVecDivVHCConsumeRootBucket this heap k nodes lin resetH
              quotient divisor with
          | .error fault => .error fault
          | .ok bucket =>
              match pairVecDivVHCExtractChecked heap bucket.nodes with
              | .error fault => .error fault
              | .ok extracted =>
                  pairVecDivVHCConsumeEqualDegree this degree extracted.1
                    bucket.coefficient bucket.nodes bucket.lin bucket.resetH
                    quotient divisor
        else
          .ok (PairVecDivVHCEqualDegreeResult.mk heap k nodes lin resetH)
  else
    .ok (PairVecDivVHCEqualDegreeResult.mk heap k nodes lin resetH)
termination_by heap.size
decreasing_by
  have hsize := extracted.2
  omega

theorem pairVecDivVHCConsumeEqualDegree_empty (this : DenseUPolyZp)
    (degree : Nat) (k : UInt64) (nodes : Array PairVecDivVHCNode)
    (lin : Array Nat) (resetH : Nat) (quotient divisor : SparsePolyZp) :
    pairVecDivVHCConsumeEqualDegree this degree #[] k nodes lin resetH
      quotient divisor =
        .ok (PairVecDivVHCEqualDegreeResult.mk #[] k nodes lin resetH) := by
  rw [pairVecDivVHCConsumeEqualDegree]
  simp

structure PairVecDivVHCHeapState where
  heap : Array Nat
  nodes : Array PairVecDivVHCNode

/-- Exact source loop `while (lin_size > 0) VHC_insert(...,
lin[--lin_size], ...)`. -/
def pairVecDivVHCReinsertLin (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode) (lin : Array Nat) :
    RawExec PairVecDivVHCHeapState :=
  if hlin : 0 < lin.size then
    let nodeIndex := lin[lin.size - 1]
    match pairVecDivVHCInsert nodeIndex heap nodes with
    | .error fault => .error fault
    | .ok (heap', nodes') =>
        pairVecDivVHCReinsertLin heap' nodes' lin.pop
  else
    .ok (PairVecDivVHCHeapState.mk heap nodes)
termination_by lin.size
decreasing_by
  simp only [Array.size_pop]
  omega

theorem pairVecDivVHCReinsertLin_empty (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode) :
    pairVecDivVHCReinsertLin heap nodes #[] =
      .ok (PairVecDivVHCHeapState.mk heap nodes) := by
  rw [pairVecDivVHCReinsertLin]
  simp

/-- Exact `while (reset_h > 0)` activation loop after emitting a quotient
cell: decrement `reset_h`, compute that node's concrete product monomial, and
insert it with the real `VHC_insert`. -/
def pairVecDivVHCActivateReset (resetH : Nat) (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode) (quotient divisor : SparsePolyZp) :
    RawExec PairVecDivVHCHeapState :=
  if hreset : 0 < resetH then
    let nodeIndex := resetH - 1
    match pairVecDivVHCActivate nodeIndex nodes quotient divisor with
    | .error fault => .error fault
    | .ok nodes' =>
        match pairVecDivVHCInsert nodeIndex heap nodes' with
        | .error fault => .error fault
        | .ok (heap', nodes'') =>
            pairVecDivVHCActivateReset nodeIndex heap' nodes'' quotient divisor
  else
    .ok (PairVecDivVHCHeapState.mk heap nodes)
termination_by resetH
decreasing_by omega

theorem pairVecDivVHCActivateReset_zero (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode) (quotient divisor : SparsePolyZp) :
    pairVecDivVHCActivateReset 0 heap nodes quotient divisor =
      .ok (PairVecDivVHCHeapState.mk heap nodes) := by
  rw [pairVecDivVHCActivateReset]
  simp

structure PairVecDivVHCFrontier where
  degree : Nat
  coefficient : UInt64
  dividendIndex : Nat

/-- Source choice between `v1_ptr` and the heap root.  For the univariate
comparator, `!comp(heapRoot, dividendTerm)` is exactly
`heapRoot.degree ≤ dividendTerm.degree`; ties therefore consume the dividend
cell first, as in C++. -/
def pairVecDivVHCSelectFrontier (dividendIndex : Nat)
    (dividend : SparsePolyZp) (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode) : RawExec PairVecDivVHCFrontier :=
  if hdividend : dividendIndex < dividend.size then
    if hheap : 0 < heap.size then
      match pairVecDivVHCMono heap[0] nodes with
      | .error fault => .error fault
      | .ok rootMono =>
          if rootMono.deg ≤ dividend[dividendIndex].1.deg then
            .ok (PairVecDivVHCFrontier.mk dividend[dividendIndex].1.deg
              dividend[dividendIndex].2.val (dividendIndex + 1))
          else
            .ok (PairVecDivVHCFrontier.mk rootMono.deg 0 dividendIndex)
    else
      .ok (PairVecDivVHCFrontier.mk dividend[dividendIndex].1.deg
        dividend[dividendIndex].2.val (dividendIndex + 1))
  else if hheap : 0 < heap.size then
    match pairVecDivVHCMono heap[0] nodes with
    | .error fault => .error fault
    | .ok rootMono =>
        .ok (PairVecDivVHCFrontier.mk rootMono.deg 0 dividendIndex)
  else
    .error .assertionFailure

theorem pairVecDivVHCSelectFrontier_index (dividendIndex : Nat)
    (dividend : SparsePolyZp) (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode) (frontier : PairVecDivVHCFrontier)
    (hrun : pairVecDivVHCSelectFrontier dividendIndex dividend heap nodes =
      .ok frontier) :
    frontier.dividendIndex = dividendIndex ∨
      frontier.dividendIndex = dividendIndex + 1 := by
  unfold pairVecDivVHCSelectFrontier at hrun
  split at hrun
  next hdividend =>
    split at hrun
    next hheap =>
      split at hrun <;> try contradiction
      next rootMono hmono =>
        split at hrun <;>
          simp only [Except.ok.injEq] at hrun <;>
          subst frontier <;> simp
    next hheap =>
      simp only [Except.ok.injEq] at hrun
      subst frontier
      simp
  next hdividend =>
    split at hrun <;> try contradiction
    next hheap =>
      split at hrun <;> try contradiction
      next rootMono hmono =>
        simp only [Except.ok.injEq] at hrun
        subst frontier
        simp

structure PairVecDivVHCIterationResult where
  dividendIndex : Nat
  heap : Array Nat
  nodes : Array PairVecDivVHCNode
  quotient : SparsePolyZp
  resetH : Nat

/-- One complete body of the general `pair_vec_div` outer loop, including
frontier selection, equal-degree bucket subtraction, quotient emission,
`reset_h` activation, and the final reverse `lin` reinsertion. -/
def pairVecDivVHCOuterIteration (this : DenseUPolyZp)
    (dividendIndex : Nat) (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode) (quotient dividend divisor : SparsePolyZp)
    (resetH : Nat) : RawExec PairVecDivVHCIterationResult := do
  if hdivisor : 0 < divisor.size then
    let frontier ← pairVecDivVHCSelectFrontier dividendIndex dividend heap nodes
    let consumed ← pairVecDivVHCConsumeEqualDegree this frontier.degree heap
      frontier.coefficient nodes #[] resetH quotient divisor
    let emit : RawExec (SparsePolyZp × PairVecDivVHCHeapState × Nat) :=
      if hcoefficient : consumed.coefficient ≠ 0 then
        if hdegree : divisor[0].1.deg ≤ frontier.degree then
          let inverse := Generated.StrictGCD.dense_upoly_zp_nmod_inv_ir this
            divisor[0].2.val
          let value := Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this
            consumed.coefficient inverse
          if hvalue : value ≠ 0 then
            let quotient' := quotient.push
              (⟨frontier.degree - divisor[0].1.deg⟩, ⟨value, this._p⟩)
            match pairVecDivVHCActivateReset consumed.resetH consumed.heap
                consumed.nodes quotient' divisor with
            | .error fault => .error fault
            | .ok activated => .ok (quotient', activated, 0)
          else
            .ok (quotient,
              PairVecDivVHCHeapState.mk consumed.heap consumed.nodes,
              consumed.resetH)
        else
          .ok (quotient,
            PairVecDivVHCHeapState.mk consumed.heap consumed.nodes,
            consumed.resetH)
      else
        .ok (quotient,
          PairVecDivVHCHeapState.mk consumed.heap consumed.nodes,
          consumed.resetH)
    let (quotient', activated, resetH') ← emit
    let reinserted ← pairVecDivVHCReinsertLin activated.heap activated.nodes
      consumed.lin
    .ok (PairVecDivVHCIterationResult.mk frontier.dividendIndex
      reinserted.heap reinserted.nodes quotient' resetH')
  else
    .error .assertionFailure

/-- Complete general-path outer `while`.  `degreeLimit` is a proof-relevant
strict upper bound on the next source frontier, not a step counter: after one
body the recursive limit becomes the degree that was actually selected.
Failure of the decrease check exposes a broken heap/canonical invariant rather
than silently truncating execution. -/
def pairVecDivVHCOuterLoop (this : DenseUPolyZp) (degreeLimit : Nat)
    (dividendIndex : Nat) (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode) (quotient dividend divisor : SparsePolyZp)
    (resetH : Nat) : RawExec SparsePolyZp :=
  if hdone : dividend.size ≤ dividendIndex ∧ heap.size = 0 then
    .ok quotient
  else
    match pairVecDivVHCSelectFrontier dividendIndex dividend heap nodes with
    | .error fault => .error fault
    | .ok frontier =>
        if hdecrease : frontier.degree < degreeLimit then
          match pairVecDivVHCOuterIteration this dividendIndex heap nodes
              quotient dividend divisor resetH with
          | .error fault => .error fault
          | .ok next =>
              pairVecDivVHCOuterLoop this frontier.degree next.dividendIndex
                next.heap next.nodes next.quotient dividend divisor next.resetH
        else
          .error .assertionFailure
termination_by degreeLimit
decreasing_by exact hdecrease

theorem pairVecDivVHCOuterLoop_done (this : DenseUPolyZp)
    (degreeLimit dividendIndex : Nat) (nodes : Array PairVecDivVHCNode)
    (quotient dividend divisor : SparsePolyZp) (resetH : Nat)
    (hdone : dividend.size ≤ dividendIndex) :
    pairVecDivVHCOuterLoop this degreeLimit dividendIndex #[] nodes quotient
      dividend divisor resetH = .ok quotient := by
  rw [pairVecDivVHCOuterLoop]
  simp [hdone]

/-- Checked entry to the current C++ general priority-heap branch.  The source
has already ruled out empty/single-term divisors and empty dividends. -/
def pairVecDivGeneralBranchIR (this : DenseUPolyZp)
    (dividend divisor : SparsePolyZp) : RawExec SparsePolyZp :=
  if hdividend : 0 < dividend.size then
    if hdivisor : 1 < divisor.size then
      let nodes := pairVecDivVHCInit divisor
      pairVecDivVHCOuterLoop this (dividend[0].1.deg + 1) 0 #[] nodes #[]
        dividend divisor (divisor.size - 1)
    else
      .error .assertionFailure
  else
    .error .assertionFailure

/-- Non-aliasing source entry used by SQF.  Its branch order is the current
C++ `pair_vec_div`: reject zero divisor, clear the fresh output and return for
zero dividend, execute the size-one loop, otherwise pass the compile-time
false univariate compression hook and enter the priority heap. -/
def pairVecDivIR (this : DenseUPolyZp)
    (dividend divisor : SparsePolyZp) : RawExec SparsePolyZp :=
  if hdivisorEmpty : divisor.size = 0 then
    .error .assertionFailure
  else if hdividendEmpty : dividend.size = 0 then
    .ok #[]
  else if hsingle : divisor.size = 1 then
    pairVecDivSingleBranchIR this dividend divisor
  else
    pairVecDivGeneralBranchIR this dividend divisor

theorem pairVecDivIR_empty_dividend (this : DenseUPolyZp)
    (dividend divisor : SparsePolyZp)
    (hdividend : dividend.size = 0) (hdivisor : divisor.size ≠ 0) :
    pairVecDivIR this dividend divisor = .ok #[] := by
  simp [pairVecDivIR, hdivisor, hdividend]

theorem pairVecDivIR_single (this : DenseUPolyZp)
    (dividend divisor : SparsePolyZp)
    (hdividend : dividend.size ≠ 0) (hdivisor : divisor.size = 1) :
    pairVecDivIR this dividend divisor =
      pairVecDivSingleBranchIR this dividend divisor := by
  simp [pairVecDivIR, hdividend, hdivisor]

theorem pairVecDivIR_general (this : DenseUPolyZp)
    (dividend divisor : SparsePolyZp)
    (hdividend : dividend.size ≠ 0) (hdivisor : 1 < divisor.size) :
    pairVecDivIR this dividend divisor =
      pairVecDivGeneralBranchIR this dividend divisor := by
  simp [pairVecDivIR, hdividend, Nat.ne_of_gt (lt_trans Nat.zero_lt_one
    hdivisor), Nat.ne_of_gt hdivisor]

/-- Exact source-order range-for loop of `__extract_pth_root`. -/
def pthRootTerm (prime : UInt64) (term : UMonomial × Zp) : UMonomial × Zp :=
  (UMonomial.mk ((term.1.deg.toUInt64 / prime).toInt64), term.2)

theorem pthRootTerm_degree (prime : UInt64) (term : UMonomial × Zp)
    (hdegree : term.1.deg < 2 ^ 63) :
    (pthRootTerm prime term).1.deg = term.1.deg / prime.toNat := by
  have hdegreeWord : term.1.deg.toUInt64.toNat = term.1.deg := by
    change (OfNat.ofNat term.1.deg : UInt64).toNat = term.1.deg
    rw [UInt64.toNat_ofNat, Nat.mod_eq_of_lt]
    omega
  have hquotient : (term.1.deg.toUInt64 / prime).toNat =
      term.1.deg / prime.toNat := by
    rw [UInt64.toNat_div, hdegreeWord]
  have hquotientBound : (term.1.deg.toUInt64 / prime).toNat < 2 ^ 63 := by
    rw [hquotient]
    exact Nat.lt_of_le_of_lt (Nat.div_le_self _ _) hdegree
  change (term.1.deg.toUInt64 / prime).toInt64.toNatClampNeg = _
  rw [UInt64_toInt64_toNatClampNeg_eq_toNat_of_lt hquotientBound,
    hquotient]

def extractPthRootLoop (index : Nat) (output source : SparsePolyZp)
    (prime : UInt64) : SparsePolyZp :=
  if h : index < source.size then
    let term := source[index]
    let next := output.push (pthRootTerm prime term)
    extractPthRootLoop (index + 1) next source prime
  else
    output
termination_by source.size - index
decreasing_by omega

theorem extractPthRootLoop_toList (index : Nat)
    (output source : SparsePolyZp) (prime : UInt64) :
    (extractPthRootLoop index output source prime).toList =
      output.toList ++ (source.toList.drop index).map (pthRootTerm prime) := by
  rw [extractPthRootLoop]
  split
  next hmore =>
    rw [extractPthRootLoop_toList]
    have hdrop := List.drop_eq_getElem_cons
      (l := source.toList) (i := index) (by simpa using hmore)
    rw [hdrop, List.map_cons, Array.getElem_toList hmore]
    simp [List.append_assoc]
  next hdone =>
    have hindex : source.toList.length ≤ index := by simpa using hdone
    rw [List.drop_eq_nil_iff.mpr hindex]
    simp
termination_by source.size - index
decreasing_by omega

theorem extractPthRootLoop_zero_canonical (p : Nat) (prime : UInt64)
    (source : SparsePolyZp) (hprime : prime.toNat = p) (hp : 0 < p)
    (hcanonical : SparsePolyZp.Canonical p source)
    (hdegree : ∀ term ∈ source.toList, term.1.deg < 2 ^ 63)
    (hdivisible : ∀ term ∈ source.toList, p ∣ term.1.deg) :
    SparsePolyZp.Canonical p (extractPthRootLoop 0 #[] source prime) := by
  rw [SparsePolyZp.Canonical, SparsePolyZp.WellFormed_arr,
    extractPthRootLoop_toList]
  simp only [List.drop_zero, List.nil_append]
  refine ⟨?_, ?_, ?_⟩
  · intro output houtput
    rw [List.mem_map] at houtput
    rcases houtput with ⟨term, hterm, rfl⟩
    exact hcanonical.1 term hterm
  · rw [List.isChain_iff_pairwise, List.pairwise_map]
    have hpairwise := List.isChain_iff_pairwise.mp hcanonical.2.1
    apply hpairwise.imp_of_mem
    intro left right hleft hright horder
    rw [pthRootTerm_degree prime left (hdegree left hleft),
      pthRootTerm_degree prime right (hdegree right hright), hprime]
    rcases hdivisible left hleft with ⟨leftQ, hleftQ⟩
    rcases hdivisible right hright with ⟨rightQ, hrightQ⟩
    have hmulOrder : p * rightQ < p * leftQ := by
      simpa [hleftQ, hrightQ] using horder
    have hquotientOrder : leftQ > rightQ :=
      (Nat.mul_lt_mul_left hp).mp hmulOrder
    rw [hleftQ, hrightQ]
    simpa [Nat.mul_div_cancel_left _ hp] using hquotientOrder
  · intro output houtput
    rw [List.mem_map] at houtput
    rcases houtput with ⟨term, hterm, rfl⟩
    exact hcanonical.2.2 term hterm

theorem contract_monomial_of_dvd {R : Type} [CommSemiring R]
    (p degree : Nat) (coefficient : R) (hp : 0 < p) (hdiv : p ∣ degree) :
    Polynomial.contract p (Polynomial.monomial degree coefficient) =
      Polynomial.monomial (degree / p) coefficient := by
  ext exponent
  rw [Polynomial.coeff_contract (Nat.ne_of_gt hp)]
  simp only [Polynomial.coeff_monomial]
  rcases hdiv with ⟨quotient, rfl⟩
  rw [Nat.mul_div_cancel_left quotient hp]
  by_cases heq : p * quotient = exponent * p
  · have hquotient : quotient = exponent :=
      Nat.mul_left_cancel hp (heq.trans (Nat.mul_comm exponent p))
    subst quotient
    simp [Nat.mul_comm]
  · have hquotient : quotient ≠ exponent := by
      intro h
      apply heq
      simp [h, Nat.mul_comm]
    simp [heq, hquotient]

theorem listSum_pthRootTerm (p : Nat) (prime : UInt64)
    (hprime : prime.toNat = p) (hp : 0 < p) :
    ∀ terms : List (UMonomial × Zp),
      (∀ term ∈ terms, term.1.deg < 2 ^ 63) →
      (∀ term ∈ terms, p ∣ term.1.deg) →
      listSum p (terms.map (pthRootTerm prime)) =
        Polynomial.contract p (listSum p terms) := by
  intro terms hdegree hdivisible
  induction terms with
  | nil => simp [listSum, Polynomial.contract]
  | cons term rest ih =>
      have htermDegree := hdegree term List.mem_cons_self
      have hrestDegree : ∀ item ∈ rest, item.1.deg < 2 ^ 63 := by
        intro item hitem
        exact hdegree item (List.mem_cons_of_mem term hitem)
      have htermDiv := hdivisible term List.mem_cons_self
      have hrestDiv : ∀ item ∈ rest, p ∣ item.1.deg := by
        intro item hitem
        exact hdivisible item (List.mem_cons_of_mem term hitem)
      rw [List.map_cons, listSum_cons, listSum_cons,
        pthRootTerm_degree prime term htermDegree,
        hprime,
        Polynomial.contract_add (Nat.ne_of_gt hp), ih hrestDegree hrestDiv,
        contract_monomial_of_dvd p term.1.deg
          (Zp.toZMod p term.2) hp htermDiv]
      rfl

def extractPthRootIR (f : SparsePolyZp) : RawExec SparsePolyZp :=
  if hnonempty : 0 < f.size then
    .ok (extractPthRootLoop 0 #[] f f[0].2.prime)
  else
    .error .assertionFailure

theorem extractPthRootIR_refines (p : Nat) [Fact (Nat.Prime p)]
    (source : SparsePolyZp) (hnonempty : 0 < source.size)
    (hcanonical : SparsePolyZp.Canonical p source)
    (hdegree : ∀ term ∈ source.toList, term.1.deg < 2 ^ 63)
    (hdivisible : ∀ term ∈ source.toList, p ∣ term.1.deg) :
    ∃ root,
      extractPthRootIR source = .ok root ∧
      SparsePolyZp.Canonical p root ∧
      SparsePolyZp.toPoly p root =
        Polynomial.contract p (SparsePolyZp.toPoly p source) := by
  have hheadMem : source[0] ∈ source.toList := by
    have hget := Array.getElem_toList hnonempty
    rw [← hget]
    exact List.getElem_mem (by simpa using hnonempty)
  have hprime : source[0].2.prime.toNat = p :=
    (hcanonical.1 source[0] hheadMem).1
  let root := extractPthRootLoop 0 #[] source source[0].2.prime
  refine ⟨root, by simp [extractPthRootIR, hnonempty, root], ?_, ?_⟩
  · exact extractPthRootLoop_zero_canonical p source[0].2.prime source
      hprime (Fact.out : Nat.Prime p).pos hcanonical hdegree hdivisible
  · unfold SparsePolyZp.toPoly root
    rw [extractPthRootLoop_toList]
    simp only [List.drop_zero, List.nil_append]
    exact listSum_pthRootTerm p source[0].2.prime hprime
      (Fact.out : Nat.Prime p).pos source.toList hdegree hdivisible

theorem listSum_coeff_of_mem_chain (p : Nat)
    (terms : List (UMonomial × Zp)) (term : UMonomial × Zp)
    (hchain : List.IsChain
      (fun left right : UMonomial × Zp => left.1.deg > right.1.deg) terms)
    (hterm : term ∈ terms) :
    (listSum p terms).coeff term.1.deg = Zp.toZMod p term.2 := by
  induction terms with
  | nil => simp at hterm
  | cons head rest ih =>
      rcases List.mem_cons.mp hterm with heq | hrest
      · subst head
        exact listSum_coeff_at_head p term rest hchain
      · rw [listSum_cons, Polynomial.coeff_add,
          Polynomial.coeff_monomial]
        have hheadGreater : head.1.deg > term.1.deg :=
          chain_gt_all_after_head head rest hchain term hrest
        rw [if_neg (by omega), zero_add]
        exact ih hchain.tail hrest

theorem canonical_degrees_dvd_of_derivative_eq_zero (p : Nat)
    [Fact (Nat.Prime p)] (source : SparsePolyZp)
    (hcanonical : SparsePolyZp.Canonical p source)
    (hderivative : Polynomial.derivative (SparsePolyZp.toPoly p source) = 0) :
    ∀ term ∈ source.toList, p ∣ term.1.deg := by
  intro term hterm
  by_cases hdegreeZero : term.1.deg = 0
  · simp [hdegreeZero]
  have hcoefficient := listSum_coeff_of_mem_chain p source.toList term
    hcanonical.2.1 hterm
  have hcoefficientNonzero :
      (SparsePolyZp.toPoly p source).coeff term.1.deg ≠ 0 := by
    unfold SparsePolyZp.toPoly
    rw [hcoefficient]
    exact Zp.toZMod_ne_zero_of_val_ne_zero p term.2
      (hcanonical.1 term hterm) (hcanonical.2.2 term hterm)
  have hderivativeCoefficient := congrArg
    (fun poly : Polynomial (ZMod p) => poly.coeff (term.1.deg - 1))
    hderivative
  change (Polynomial.derivative (SparsePolyZp.toPoly p source)).coeff
      (term.1.deg - 1) = 0 at hderivativeCoefficient
  rw [Polynomial.coeff_derivative] at hderivativeCoefficient
  have hindex : term.1.deg - 1 + 1 = term.1.deg := by omega
  rw [hindex] at hderivativeCoefficient
  have hfactor : ((term.1.deg - 1 : Nat) : ZMod p) + 1 =
      (term.1.deg : ZMod p) := by
    calc
      ((term.1.deg - 1 : Nat) : ZMod p) + 1 =
          ((term.1.deg - 1 : Nat) : ZMod p) + ((1 : Nat) : ZMod p) := by
        norm_num
      _ = ((term.1.deg - 1 + 1 : Nat) : ZMod p) := by push_cast; rfl
      _ = (term.1.deg : ZMod p) := congrArg (fun n : Nat => (n : ZMod p)) hindex
  rw [hfactor] at hderivativeCoefficient
  have hdegreeCast : (term.1.deg : ZMod p) = 0 := by
    rcases mul_eq_zero.mp hderivativeCoefficient with hbad | hzero
    · exact False.elim (hcoefficientNonzero hbad)
    · exact hzero
  exact (ZMod.natCast_eq_zero_iff term.1.deg p).mp hdegreeCast

theorem extractPthRootIR_refines_of_derivative_zero (p : Nat)
    [Fact (Nat.Prime p)] (source : SparsePolyZp)
    (hnonempty : 0 < source.size)
    (hcanonical : SparsePolyZp.Canonical p source)
    (hdegree : ∀ term ∈ source.toList, term.1.deg < 2 ^ 63)
    (hderivative : Polynomial.derivative (SparsePolyZp.toPoly p source) = 0) :
    ∃ root,
      extractPthRootIR source = .ok root ∧
      SparsePolyZp.Canonical p root ∧
      SparsePolyZp.toPoly p root =
        Polynomial.contract p (SparsePolyZp.toPoly p source) := by
  apply extractPthRootIR_refines p source hnonempty hcanonical hdegree
  exact canonical_degrees_dvd_of_derivative_eq_zero p source hcanonical
    hderivative

theorem extractPthRootLoop_size (index : Nat) (output source : SparsePolyZp)
    (prime : UInt64) :
    (extractPthRootLoop index output source prime).size =
      output.size + (source.size - index) := by
  rw [extractPthRootLoop]
  split
  next hmore =>
    rw [extractPthRootLoop_size]
    simp
    omega
  next hdone =>
    simp
    omega
termination_by source.size - index
decreasing_by omega

/-- Exact result-copy loop used after recursive SQF contraction. -/
def scaleMultiplicityLoop (index : Nat)
    (source output : Array (SparsePolyZp × UInt64)) (prime : UInt64) :
    Array (SparsePolyZp × UInt64) :=
  if h : index < source.size then
    let item := source[index]
    scaleMultiplicityLoop (index + 1) source
      (output.push (item.1, item.2 * prime)) prime
  else
    output
termination_by source.size - index
decreasing_by omega

theorem scaleMultiplicityLoop_toList (index : Nat)
    (source output : Array (SparsePolyZp × UInt64)) (prime : UInt64) :
    (scaleMultiplicityLoop index source output prime).toList =
      output.toList ++ (source.toList.drop index).map
        (fun item => (item.1, item.2 * prime)) := by
  rw [scaleMultiplicityLoop]
  split
  next hmore =>
    rw [scaleMultiplicityLoop_toList]
    have hdrop := List.drop_eq_getElem_cons
      (l := source.toList) (i := index) (by simpa using hmore)
    rw [hdrop, List.map_cons, Array.getElem_toList hmore]
    simp [List.append_assoc]
  next hdone =>
    have hindex : source.toList.length ≤ index := by simpa using hdone
    rw [List.drop_eq_nil_iff.mpr hindex]
    simp
termination_by source.size - index
decreasing_by omega

theorem scaleMultiplicityLoop_size (index : Nat)
    (source output : Array (SparsePolyZp × UInt64)) (prime : UInt64) :
    (scaleMultiplicityLoop index source output prime).size =
      output.size + (source.size - index) := by
  rw [scaleMultiplicityLoop]
  split
  next hmore =>
    rw [scaleMultiplicityLoop_size]
    simp
    omega
  next hdone =>
    simp
    omega
termination_by source.size - index
decreasing_by omega

theorem scaleMultiplicityLoop_toPolyList (p : Nat)
    (source : Array (SparsePolyZp × UInt64)) (prime : UInt64)
    (hprime : prime.toNat = p)
    (hnowrap : ∀ item ∈ source.toList,
      item.2.toNat * prime.toNat < UInt64.size) :
    toPolyList (scaleMultiplicityLoop 0 source #[] prime) p =
      (toPolyList source p).map (fun item => (item.1, item.2 * p)) := by
  unfold toPolyList
  rw [Array.toList_map, scaleMultiplicityLoop_toList]
  simp only [List.drop_zero, List.nil_append, List.map_map]
  rw [Array.toList_map, List.map_map]
  apply List.map_congr_left
  intro item hitem
  simp only [Function.comp_apply]
  congr 2
  rw [UInt64.toNat_mul, Nat.mod_eq_of_lt (hnowrap item hitem), hprime]

end StrictSquarefreeZp
end Refinement
