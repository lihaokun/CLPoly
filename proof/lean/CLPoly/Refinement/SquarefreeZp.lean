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

theorem pairVecDivVHCInit_get (divisor : SparsePolyZp) (i : Nat)
    (hi : i < divisor.size - 1) :
    (pairVecDivVHCInit divisor)[i]? = some (pairVecDivVHCInitialNode i) := by
  have harr : pairVecDivVHCInit divisor =
      (Array.range (divisor.size - 1)).map pairVecDivVHCInitialNode := by
    apply Array.toList_inj.mp
    simp [pairVecDivVHCInit_toList]
  rw [harr]
  simp [hi]

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

theorem pairVecDivVHCSetNext_nodes_size
    (nodeIndex : Nat) (next : Option Nat)
    (nodes nodes' : Array PairVecDivVHCNode)
    (hrun : pairVecDivVHCSetNext nodeIndex next nodes = .ok nodes') :
    nodes'.size = nodes.size := by
  unfold pairVecDivVHCSetNext at hrun
  split at hrun <;> try contradiction
  next hn =>
    simp only [Except.ok.injEq] at hrun
    subst nodes'
    simp

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

theorem pairVecDivVHCMono_eq_ok_iff (nodeIndex : Nat)
    (nodes : Array PairVecDivVHCNode) (mono : UMonomial) :
    pairVecDivVHCMono nodeIndex nodes = .ok mono ↔
      ∃ node, nodes[nodeIndex]? = some node ∧ node.mono = some mono := by
  unfold pairVecDivVHCMono
  constructor
  · intro hrun
    split at hrun <;> try contradiction
    next hn =>
      split at hrun <;> try contradiction
      next storedMono hstoredMono =>
        simp only [Except.ok.injEq] at hrun
        subst storedMono
        exact ⟨nodes[nodeIndex], Array.getElem?_eq_getElem hn, hstoredMono⟩
  · rintro ⟨node, hget, hmono⟩
    have hn : nodeIndex < nodes.size := by
      by_contra hnot
      rw [Array.getElem?_eq_none (by omega)] at hget
      contradiction
    simp only [hn, ↓reduceDIte]
    rw [Array.getElem?_eq_getElem hn] at hget
    simp only [Option.some.injEq] at hget
    subst node
    rw [hmono]

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

theorem pairVecDivVHCBubble_stop_get
    (i stop newNode : Nat) (heap heap' : Array Nat)
    (hrun : pairVecDivVHCBubble i stop newNode heap = .ok heap') :
    heap'[stop]? = some newNode := by
  rw [pairVecDivVHCBubble] at hrun
  split at hrun <;> try contradiction
  next hi =>
    split at hrun <;> try contradiction
    next hstop =>
      split at hrun
      next heq =>
        subst i
        simp only [Except.ok.injEq] at hrun
        subst heap'
        exact Array.getElem?_set_self hi
      next heq =>
        dsimp only at hrun
        split at hrun <;> try contradiction
        next hp =>
          exact pairVecDivVHCBubble_stop_get (pairVecDivVHCParent i) stop
            newNode (heap.set i heap[pairVecDivVHCParent i]) heap' hrun
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

/-- Exact comparison trace of the generated upward anchor search.  A climb
records that the visited key is strictly below the new degree; the stop records
the first ancestor whose key is at least the new degree. -/
inductive PairVecDivVHCFindAnchorTrace (newDegree : Nat)
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode) : Nat → Nat → Prop
  | stop (i head : Nat) (mono : UMonomial)
      (hhead : heap[i]? = some head)
      (hmono : pairVecDivVHCMono head nodes = .ok mono)
      (hle : newDegree ≤ mono.deg) :
      PairVecDivVHCFindAnchorTrace newDegree heap nodes i i
  | climb (i head anchor : Nat) (mono : UMonomial)
      (hhead : heap[i]? = some head)
      (hmono : pairVecDivVHCMono head nodes = .ok mono)
      (hlt : mono.deg < newDegree) (hpos : 0 < i)
      (htail : PairVecDivVHCFindAnchorTrace newDegree heap nodes
        (pairVecDivVHCParent i) anchor) :
      PairVecDivVHCFindAnchorTrace newDegree heap nodes i anchor

theorem pairVecDivVHCFindAnchor_trace
    (newDegree i anchor : Nat) (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode)
    (hrun : pairVecDivVHCFindAnchor newDegree i heap nodes = .ok anchor) :
    PairVecDivVHCFindAnchorTrace newDegree heap nodes i anchor := by
  unfold pairVecDivVHCFindAnchor at hrun
  split at hrun <;> try contradiction
  next hi =>
    cases hread : pairVecDivVHCMono heap[i] nodes with
    | error fault => simp [hread] at hrun
    | ok mono =>
        simp only [hread] at hrun
        split at hrun
        next hgreater =>
          split at hrun <;> try contradiction
          next hzero =>
            apply PairVecDivVHCFindAnchorTrace.climb i heap[i] anchor mono
              (Array.getElem?_eq_getElem hi) hread (by omega)
              (Nat.pos_of_ne_zero hzero)
            exact pairVecDivVHCFindAnchor_trace newDegree
              (pairVecDivVHCParent i) anchor heap nodes hrun
        next hgreater =>
          simp only [Except.ok.injEq] at hrun
          subst anchor
          exact PairVecDivVHCFindAnchorTrace.stop i heap[i] mono
            (Array.getElem?_eq_getElem hi) hread (by omega)
termination_by i
decreasing_by
  unfold pairVecDivVHCParent
  have hhalf : (i - 1) / 2 ≤ i - 1 := Nat.div_le_self _ _
  omega

theorem PairVecDivVHCFindAnchorTrace.anchor_read
    (newDegree : Nat) (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode) (i anchor : Nat)
    (htrace : PairVecDivVHCFindAnchorTrace newDegree heap nodes i anchor) :
    ∃ head mono, heap[anchor]? = some head ∧
      pairVecDivVHCMono head nodes = .ok mono ∧ newDegree ≤ mono.deg := by
  induction htrace with
  | stop i head mono hhead hmono hle => exact ⟨head, mono, hhead, hmono, hle⟩
  | climb i head anchor mono hhead hmono hlt hpos htail ih => exact ih

theorem PairVecDivVHCFindAnchorTrace.anchor_le_start
    (newDegree : Nat) (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode) (i anchor : Nat)
    (htrace : PairVecDivVHCFindAnchorTrace newDegree heap nodes i anchor) :
    anchor ≤ i := by
  induction htrace with
  | stop => exact Nat.le_refl _
  | climb i head anchor mono hhead hmono hlt hpos htail ih =>
      have hparentLt : pairVecDivVHCParent i < i := by
        unfold pairVecDivVHCParent
        have hhalf : (i - 1) / 2 ≤ i - 1 := Nat.div_le_self _ _
        omega
      exact Nat.le_trans ih (Nat.le_of_lt hparentLt)

/-- A write strictly below the remaining ancestor-search path cannot change
that path's recorded comparisons. -/
theorem PairVecDivVHCFindAnchorTrace.set_above
    (newDegree : Nat) (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode) (i anchor target replacement : Nat)
    (htarget : target < heap.size) (hiTarget : i < target)
    (htrace : PairVecDivVHCFindAnchorTrace newDegree heap nodes i anchor) :
    PairVecDivVHCFindAnchorTrace newDegree (heap.set target replacement)
      nodes i anchor := by
  induction htrace with
  | stop i head mono hhead hmono hle =>
      apply PairVecDivVHCFindAnchorTrace.stop i head mono
      · rw [Array.getElem?_set_ne htarget (by omega)]
        exact hhead
      · exact hmono
      · exact hle
  | climb i head anchor mono hhead hmono hlt hpos htail ih =>
      apply PairVecDivVHCFindAnchorTrace.climb i head anchor mono
      · rw [Array.getElem?_set_ne htarget (by omega)]
        exact hhead
      · exact hmono
      · exact hlt
      · exact hpos
      · apply ih
        have hparentLt : pairVecDivVHCParent i < i := by
          unfold pairVecDivVHCParent
          have hhalf : (i - 1) / 2 ≤ i - 1 := Nat.div_le_self _ _
          omega
        exact Nat.lt_trans hparentLt hiTarget

theorem PairVecDivVHCFindAnchorTrace.push
    (newDegree : Nat) (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode) (i anchor replacement : Nat)
    (hi : i < heap.size)
    (htrace : PairVecDivVHCFindAnchorTrace newDegree heap nodes i anchor) :
    PairVecDivVHCFindAnchorTrace newDegree (heap.push replacement)
      nodes i anchor := by
  induction htrace with
  | stop i head mono hhead hmono hle =>
      apply PairVecDivVHCFindAnchorTrace.stop i head mono
      · simpa [Array.getElem?_push, Nat.ne_of_lt hi] using hhead
      · exact hmono
      · exact hle
  | climb i head anchor mono hhead hmono hlt hpos htail ih =>
      apply PairVecDivVHCFindAnchorTrace.climb i head anchor mono
      · simpa [Array.getElem?_push, Nat.ne_of_lt hi] using hhead
      · exact hmono
      · exact hlt
      · exact hpos
      · apply ih
        have hparentLt : pairVecDivVHCParent i < i := by
          unfold pairVecDivVHCParent
          have hhalf : (i - 1) / 2 ≤ i - 1 := Nat.div_le_self _ _
          omega
        exact Nat.lt_trans hparentLt hi

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
    RawExec (Array Nat × Array PairVecDivVHCNode) :=
  match pairVecDivVHCMono newNode nodes with
  | .error fault => .error fault
  | .ok newMono =>
      if hempty : heap.size = 0 then
        match pairVecDivVHCSetNext newNode none nodes with
        | .error fault => .error fault
        | .ok nodes' => .ok (#[newNode], nodes')
      else
        let rootNode := heap[0]'(Nat.pos_of_ne_zero hempty)
        match pairVecDivVHCMono rootNode nodes with
        | .error fault => .error fault
        | .ok rootMono =>
            if hequal : newMono.deg = rootMono.deg then
              match pairVecDivVHCSetNext newNode (some rootNode) nodes with
              | .error fault => .error fault
              | .ok nodes' => .ok (heap.set 0 newNode, nodes')
            else if hgreater : newMono.deg > rootMono.deg then
              match pairVecDivVHCSetNext newNode none nodes with
              | .error fault => .error fault
              | .ok nodes' =>
                  match pairVecDivVHCBubble heap.size 0 newNode
                      (heap.push newNode) with
                  | .error fault => .error fault
                  | .ok heap' => .ok (heap', nodes')
            else
              let firstAnchor := pairVecDivVHCParent heap.size
              match pairVecDivVHCFindAnchor newMono.deg firstAnchor heap nodes with
              | .error fault => .error fault
              | .ok anchor =>
                  if ha : anchor < heap.size then
                    let anchorNode := heap[anchor]
                    match pairVecDivVHCMono anchorNode nodes with
                    | .error fault => .error fault
                    | .ok anchorMono =>
                        if hequalAnchor : newMono.deg = anchorMono.deg then
                          match pairVecDivVHCSetNext newNode (some anchorNode)
                              nodes with
                          | .error fault => .error fault
                          | .ok nodes' =>
                              .ok (heap.set anchor newNode, nodes')
                        else
                          match pairVecDivVHCSetNext newNode none nodes with
                          | .error fault => .error fault
                          | .ok nodes' =>
                              match pairVecDivVHCBubbleBelow heap.size anchor
                                  newNode (heap.push newNode) with
                              | .error fault => .error fault
                              | .ok heap' => .ok (heap', nodes')
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

/-- Sift-down only writes its current path.  Every slot strictly before the
current node is unchanged by the complete recursive execution. -/
theorem pairVecDivVHCSiftDown_get_before
    (i child limit lastNode : Nat) (heap heap' : Array Nat)
    (nodes : Array PairVecDivVHCNode) (hipath : i < child)
    (hrun : pairVecDivVHCSiftDown i child limit lastNode heap nodes =
      .ok heap') :
    ∀ j < i, heap'[j]? = heap[j]? := by
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
                  intro j hj
                  have hrec := pairVecDivVHCSiftDown_get_before child
                    (child * 2 + 1) limit lastNode
                    (heap.set i heap[child]) heap' nodes (by omega) hrun j
                    (by omega)
                  rw [hrec, Array.getElem?_set_ne hi (by omega)]
                next hgreater =>
                  simp only [Except.ok.injEq] at hrun
                  subst heap'
                  intro j hj
                  rw [Array.getElem?_set_ne hi (by omega)]
              · simp only [hselected, ↓reduceIte] at hrun
                split at hrun
                next hgreater =>
                  intro j hj
                  have hrec := pairVecDivVHCSiftDown_get_before (child + 1)
                    ((child + 1) * 2 + 1) limit lastNode
                    (heap.set i heap[child + 1]) heap' nodes (by omega) hrun j
                    (by omega)
                  rw [hrec, Array.getElem?_set_ne hi (by omega)]
                next hgreater =>
                  simp only [Except.ok.injEq] at hrun
                  subst heap'
                  intro j hj
                  rw [Array.getElem?_set_ne hi (by omega)]
      next hchild =>
        simp only [Except.ok.injEq] at hrun
        subst heap'
        intro j hj
        rw [Array.getElem?_set_ne hi (by omega)]
termination_by limit - child
decreasing_by
  all_goals omega

/-- On a heap with both root children present, the first source sift step
writes the maximum of the two children and the saved last sentinel at slot
zero; recursive descent cannot overwrite that slot. -/
theorem pairVecDivVHCSiftDown_root_dominates_candidates
    (heap shifted : Array Nat) (nodes : Array PairVecDivVHCNode)
    (leftMono rightMono lastMono : UMonomial)
    (hsize : 2 < heap.size)
    (hleft : pairVecDivVHCMono heap[1] nodes = .ok leftMono)
    (hright : pairVecDivVHCMono heap[2] nodes = .ok rightMono)
    (hlast : pairVecDivVHCMono heap[heap.size - 1] nodes = .ok lastMono)
    (hrun : pairVecDivVHCSiftDown 0 1 (heap.size - 1)
      heap[heap.size - 1] heap nodes = .ok shifted) :
    ∃ rootHead rootMono,
      shifted[0]? = some rootHead ∧
        pairVecDivVHCMono rootHead nodes = .ok rootMono ∧
        leftMono.deg ≤ rootMono.deg ∧ rightMono.deg ≤ rootMono.deg ∧
        lastMono.deg ≤ rootMono.deg := by
  have hzero : 0 < heap.size := by omega
  have hlimit : heap.size - 1 < heap.size := by omega
  have hchild : 1 < heap.size - 1 := by omega
  rw [pairVecDivVHCSiftDown] at hrun
  simp only [hzero, hlimit, hchild, ↓reduceDIte] at hrun
  simp only [hleft, hright, hlast] at hrun
  by_cases hselected : leftMono.deg > rightMono.deg
  · simp only [hselected, ↓reduceIte] at hrun
    by_cases hgreater : leftMono.deg > lastMono.deg
    · simp only [hgreater, ↓reduceDIte] at hrun
      have hroot := pairVecDivVHCSiftDown_get_before 1 3
        (heap.size - 1) heap[heap.size - 1]
        (heap.set 0 heap[1]) shifted nodes (by omega) hrun 0 (by omega)
      rw [Array.getElem?_set_self hzero] at hroot
      exact ⟨heap[1], leftMono, hroot, hleft, Nat.le_refl _,
        Nat.le_of_lt hselected, Nat.le_of_lt hgreater⟩
    · simp only [hgreater, ↓reduceDIte, Except.ok.injEq] at hrun
      subst shifted
      refine ⟨heap[heap.size - 1], lastMono,
        Array.getElem?_set_self hzero, hlast, ?_, ?_, Nat.le_refl _⟩
      · exact Nat.le_of_not_gt hgreater
      · exact Nat.le_trans (Nat.le_of_lt hselected)
          (Nat.le_of_not_gt hgreater)
  · simp only [hselected, ↓reduceIte] at hrun
    by_cases hgreater : rightMono.deg > lastMono.deg
    · simp only [hgreater, ↓reduceDIte] at hrun
      have hroot := pairVecDivVHCSiftDown_get_before 2 5
        (heap.size - 1) heap[heap.size - 1]
        (heap.set 0 heap[2]) shifted nodes (by omega) hrun 0 (by omega)
      rw [Array.getElem?_set_self hzero] at hroot
      exact ⟨heap[2], rightMono, hroot, hright, Nat.le_of_not_gt hselected,
        Nat.le_refl _, Nat.le_of_lt hgreater⟩
    · simp only [hgreater, ↓reduceDIte, Except.ok.injEq] at hrun
      subst shifted
      refine ⟨heap[heap.size - 1], lastMono,
        Array.getElem?_set_self hzero, hlast, ?_, ?_, Nat.le_refl _⟩
      · exact Nat.le_trans (Nat.le_of_not_gt hselected)
          (Nat.le_of_not_gt hgreater)
      · exact Nat.le_of_not_gt hgreater

/-- The root-candidate argument holds at every recursive sift hole.  The
value finally stored at `i` dominates both children inspected at that level
and the saved last sentinel. -/
theorem pairVecDivVHCSiftDown_hole_dominates_candidates
    (i child limit lastNode : Nat) (heap shifted : Array Nat)
    (nodes : Array PairVecDivVHCNode)
    (leftMono rightMono lastMono : UMonomial)
    (hipath : i < child) (hchildLimit : child < limit)
    (hlimitBound : limit < heap.size)
    (hleft : pairVecDivVHCMono heap[child] nodes = .ok leftMono)
    (hright : pairVecDivVHCMono heap[child + 1] nodes = .ok rightMono)
    (hlast : pairVecDivVHCMono lastNode nodes = .ok lastMono)
    (hrun : pairVecDivVHCSiftDown i child limit lastNode heap nodes =
      .ok shifted) :
    ∃ holeHead holeMono,
      shifted[i]? = some holeHead ∧
        pairVecDivVHCMono holeHead nodes = .ok holeMono ∧
        leftMono.deg ≤ holeMono.deg ∧ rightMono.deg ≤ holeMono.deg ∧
        lastMono.deg ≤ holeMono.deg := by
  rw [pairVecDivVHCSiftDown] at hrun
  split at hrun <;> try contradiction
  next hi =>
      simp only [hleft, hright, hlast] at hrun
      by_cases hselected : leftMono.deg > rightMono.deg
      · simp only [hselected, ↓reduceIte] at hrun
        by_cases hgreater : leftMono.deg > lastMono.deg
        · simp only [hgreater, ↓reduceDIte] at hrun
          have hroot := pairVecDivVHCSiftDown_get_before child
            (child * 2 + 1) limit lastNode (heap.set i heap[child]) shifted
            nodes (by omega) hrun i hipath
          rw [Array.getElem?_set_self hi] at hroot
          exact ⟨heap[child], leftMono, hroot, hleft, Nat.le_refl _,
            Nat.le_of_lt hselected, Nat.le_of_lt hgreater⟩
        · simp only [hgreater, ↓reduceDIte, Except.ok.injEq] at hrun
          subst shifted
          exact ⟨lastNode, lastMono, Array.getElem?_set_self hi, hlast,
            Nat.le_of_not_gt hgreater,
            Nat.le_trans (Nat.le_of_lt hselected)
              (Nat.le_of_not_gt hgreater), Nat.le_refl _⟩
      · simp only [hselected, ↓reduceIte] at hrun
        by_cases hgreater : rightMono.deg > lastMono.deg
        · simp only [hgreater, ↓reduceDIte] at hrun
          have hroot := pairVecDivVHCSiftDown_get_before (child + 1)
            ((child + 1) * 2 + 1) limit lastNode
            (heap.set i heap[child + 1]) shifted nodes (by omega) hrun i
            (by omega)
          rw [Array.getElem?_set_self hi] at hroot
          exact ⟨heap[child + 1], rightMono, hroot, hright,
            Nat.le_of_not_gt hselected, Nat.le_refl _, Nat.le_of_lt hgreater⟩
        · simp only [hgreater, ↓reduceDIte, Except.ok.injEq] at hrun
          subst shifted
          exact ⟨lastNode, lastMono, Array.getElem?_set_self hi, hlast,
            Nat.le_trans (Nat.le_of_not_gt hselected)
              (Nat.le_of_not_gt hgreater), Nat.le_of_not_gt hgreater,
            Nat.le_refl _⟩

/-- A value occurs in an array only at one distinguished slot.  This is the
local form of heap-head uniqueness needed to show that extracting slot zero
cannot leave the consumed bucket head elsewhere in the active heap. -/
def PairVecDivVHCOnlyAt (heap : Array Nat) (nodeIndex slot : Nat) : Prop :=
  ∀ (i : Nat), heap[i]? = some nodeIndex → i = slot

/-- A node index is absent from every observable array slot. -/
def PairVecDivVHCExcludes (heap : Array Nat) (nodeIndex : Nat) : Prop :=
  ∀ (i : Nat), heap[i]? ≠ some nodeIndex

/-- Every active value in `target` originates at an active slot of `source`.
This intentionally records provenance rather than merely a range bound. -/
def PairVecDivVHCValuesFrom (target source : Array Nat) : Prop :=
  ∀ (i value : Nat), target[i]? = some value →
    ∃ j : Nat, source[j]? = some value

/-- Values below `limit` are unique, except that the current sift hole may
temporarily duplicate the child that was copied into it. -/
def PairVecDivVHCPrefixUniqueExcept (heap : Array Nat)
    (limit hole : Nat) : Prop :=
  ∀ (left right value : Nat), left < limit → right < limit →
    heap[left]? = some value → heap[right]? = some value →
    left = right ∨ left = hole ∨ right = hole

def PairVecDivVHCPrefixUnique (heap : Array Nat) (limit : Nat) : Prop :=
  ∀ (left right value : Nat), left < limit → right < limit →
    heap[left]? = some value → heap[right]? = some value → left = right

/-- Within the active prefix, the saved last node may occur only at the
current hole.  Initially it lies just outside the prefix. -/
def PairVecDivVHCPrefixOnlyAt (heap : Array Nat) (limit saved hole : Nat) :
    Prop :=
  ∀ (slot : Nat), slot < limit → heap[slot]? = some saved → slot = hole

theorem pairVecDivVHCSet_child_prefixUniqueExcept
    (heap : Array Nat) (limit hole child : Nat)
    (hhole : hole < heap.size) (hchild : child < heap.size)
    (hchildLimit : child < limit) (hne : hole ≠ child)
    (hunique : PairVecDivVHCPrefixUniqueExcept heap limit hole) :
    PairVecDivVHCPrefixUniqueExcept (heap.set hole heap[child]) limit child := by
  intro left right value hleft hright hleftValue hrightValue
  by_cases hleftHole : hole = left
  · subst left
    by_cases hrightHole : hole = right
    · exact Or.inl hrightHole
    · rw [Array.getElem?_set_ne hhole hrightHole] at hrightValue
      have hcopied : heap[child]? = some value := by
        have heq : heap[child] = value := by
          simpa [Array.getElem?_set, hhole] using hleftValue
        rw [Array.getElem?_eq_getElem hchild, heq]
      rcases hunique child right value hchildLimit hright hcopied
        hrightValue with heq | heq | heq
      · exact Or.inr (Or.inr heq.symm)
      · exact False.elim (hne heq.symm)
      · exact False.elim (hrightHole heq.symm)
  · rw [Array.getElem?_set_ne hhole hleftHole] at hleftValue
    by_cases hrightHole : hole = right
    · subst right
      have hcopied : heap[child]? = some value := by
        have heq : heap[child] = value := by
          simpa [Array.getElem?_set, hhole] using hrightValue
        rw [Array.getElem?_eq_getElem hchild, heq]
      rcases hunique left child value hleft hchildLimit hleftValue
        hcopied with heq | heq | heq
      · exact Or.inr (Or.inl heq)
      · exact False.elim (hleftHole heq.symm)
      · exact False.elim (hne heq.symm)
    · rw [Array.getElem?_set_ne hhole hrightHole] at hrightValue
      rcases hunique left right value hleft hright hleftValue
        hrightValue with heq | heq | heq
      · exact Or.inl heq
      · exact False.elim (hleftHole heq.symm)
      · exact False.elim (hrightHole heq.symm)

theorem pairVecDivVHCSet_child_prefixOnlyAt
    (heap : Array Nat) (limit saved hole child : Nat)
    (hhole : hole < heap.size) (hchild : child < heap.size)
    (hchildLimit : child < limit)
    (hne : hole ≠ child)
    (hsaved : PairVecDivVHCPrefixOnlyAt heap limit saved hole) :
    PairVecDivVHCPrefixOnlyAt (heap.set hole heap[child]) limit saved child := by
  intro slot hslot hvalue
  by_cases hslotHole : hole = slot
  · subst slot
    have hchildSaved : heap[child]? = some saved := by
      have heq : heap[child] = saved := by
        simpa [Array.getElem?_set, hhole] using hvalue
      rw [Array.getElem?_eq_getElem hchild, heq]
    exact False.elim (hne (hsaved child hchildLimit hchildSaved).symm)
  · rw [Array.getElem?_set_ne hhole hslotHole] at hvalue
    exact False.elim (hslotHole (hsaved slot hslot hvalue).symm)

theorem pairVecDivVHCSet_saved_prefixUnique
    (heap : Array Nat) (limit saved hole : Nat)
    (hhole : hole < heap.size) (hholeLimit : hole < limit)
    (hunique : PairVecDivVHCPrefixUniqueExcept heap limit hole)
    (hsaved : PairVecDivVHCPrefixOnlyAt heap limit saved hole) :
    PairVecDivVHCPrefixUnique (heap.set hole saved) limit := by
  intro left right value hleft hright hleftValue hrightValue
  by_cases hleftHole : hole = left
  · subst left
    have hvalueSaved : saved = value := by
      simpa [Array.getElem?_set, hhole] using hleftValue
    subst value
    by_cases hrightHole : hole = right
    · exact hrightHole
    · rw [Array.getElem?_set_ne hhole hrightHole] at hrightValue
      exact (hsaved right hright hrightValue).symm
  · rw [Array.getElem?_set_ne hhole hleftHole] at hleftValue
    by_cases hrightHole : hole = right
    · subst right
      have hvalueSaved : saved = value := by
        simpa [Array.getElem?_set, hhole] using hrightValue
      subst value
      exact hsaved left hleft hleftValue
    · rw [Array.getElem?_set_ne hhole hrightHole] at hrightValue
      rcases hunique left right value hleft hright hleftValue hrightValue with
        heq | heq | heq
      · exact heq
      · exact False.elim (hleftHole heq.symm)
      · exact False.elim (hrightHole heq.symm)

theorem pairVecDivVHCValuesFrom_refl (heap : Array Nat) :
    PairVecDivVHCValuesFrom heap heap := by
  intro i value hi
  exact ⟨i, hi⟩

theorem pairVecDivVHCSet_valuesFrom (target source : Array Nat)
    (slot replacement : Nat) (hslot : slot < target.size)
    (hfrom : PairVecDivVHCValuesFrom target source)
    (hreplacement : ∃ j : Nat, source[j]? = some replacement) :
    PairVecDivVHCValuesFrom (target.set slot replacement) source := by
  intro i value hi
  by_cases hislot : slot = i
  · subst i
    have hvalue : replacement = value := by
      simpa [Array.getElem?_set, hslot] using hi
    subst value
    exact hreplacement
  · rw [Array.getElem?_set_ne hslot hislot] at hi
    exact hfrom i value hi

theorem pairVecDivVHCBubble_valuesFrom
    (i stop newNode : Nat) (heap heap' source : Array Nat)
    (hfrom : PairVecDivVHCValuesFrom heap source)
    (hnew : ∃ slot : Nat, source[slot]? = some newNode)
    (hrun : pairVecDivVHCBubble i stop newNode heap = .ok heap') :
    PairVecDivVHCValuesFrom heap' source := by
  rw [pairVecDivVHCBubble] at hrun
  split at hrun <;> try contradiction
  next hi =>
    split at hrun <;> try contradiction
    next hstop =>
      split at hrun
      next heq =>
        simp only [Except.ok.injEq] at hrun
        subst heap'
        exact pairVecDivVHCSet_valuesFrom heap source i newNode hi hfrom hnew
      next heq =>
        dsimp only at hrun
        split at hrun <;> try contradiction
        next hp =>
          have hparent := hfrom (pairVecDivVHCParent i)
            heap[pairVecDivVHCParent i] (by
              rw [Array.getElem?_eq_getElem hp])
          exact pairVecDivVHCBubble_valuesFrom (pairVecDivVHCParent i) stop
            newNode (heap.set i heap[pairVecDivVHCParent i]) heap' source
            (pairVecDivVHCSet_valuesFrom heap source i
              heap[pairVecDivVHCParent i] hi hfrom hparent) hnew hrun
termination_by i
decreasing_by
  have hipos : 0 < i := by omega
  unfold pairVecDivVHCParent
  have hhalf : (i - 1) / 2 ≤ i - 1 := Nat.div_le_self _ _
  omega

theorem pairVecDivVHCBubbleBelow_valuesFrom
    (i anchor newNode : Nat) (heap heap' source : Array Nat)
    (hfrom : PairVecDivVHCValuesFrom heap source)
    (hnew : ∃ slot : Nat, source[slot]? = some newNode)
    (hrun : pairVecDivVHCBubbleBelow i anchor newNode heap = .ok heap') :
    PairVecDivVHCValuesFrom heap' source := by
  rw [pairVecDivVHCBubbleBelow] at hrun
  split at hrun <;> try contradiction
  next hi =>
    split at hrun <;> try contradiction
    next hpos =>
      dsimp only at hrun
      split at hrun
      next heq =>
        simp only [Except.ok.injEq] at hrun
        subst heap'
        exact pairVecDivVHCSet_valuesFrom heap source i newNode hi hfrom hnew
      next heq =>
        split at hrun <;> try contradiction
        next hp =>
          have hparent := hfrom (pairVecDivVHCParent i)
            heap[pairVecDivVHCParent i] (by
              rw [Array.getElem?_eq_getElem hp])
          exact pairVecDivVHCBubbleBelow_valuesFrom (pairVecDivVHCParent i)
            anchor newNode (heap.set i heap[pairVecDivVHCParent i]) heap'
            source (pairVecDivVHCSet_valuesFrom heap source i
              heap[pairVecDivVHCParent i] hi hfrom hparent) hnew hrun
termination_by i
decreasing_by
  unfold pairVecDivVHCParent
  have hhalf : (i - 1) / 2 ≤ i - 1 := Nat.div_le_self _ _
  omega

theorem pairVecDivVHCBubbleBelow_size (i anchor newNode : Nat)
    (heap heap' : Array Nat)
    (hrun : pairVecDivVHCBubbleBelow i anchor newNode heap = .ok heap') :
    heap'.size = heap.size := by
  rw [pairVecDivVHCBubbleBelow] at hrun
  split at hrun <;> try contradiction
  next hi =>
    split at hrun <;> try contradiction
    next hpos =>
      dsimp only at hrun
      split at hrun
      next heq =>
        simp only [Except.ok.injEq] at hrun
        subst heap'
        simp
      next heq =>
        split at hrun <;> try contradiction
        next hp =>
          have hrec := pairVecDivVHCBubbleBelow_size
            (pairVecDivVHCParent i) anchor newNode
            (heap.set i heap[pairVecDivVHCParent i]) heap' hrun
          simpa using hrec
termination_by i
decreasing_by
  unfold pairVecDivVHCParent
  have hhalf : (i - 1) / 2 ≤ i - 1 := Nat.div_le_self _ _
  omega

theorem pairVecDivVHCBubble_prefixUnique
    (i stop newNode limit : Nat) (heap heap' : Array Nat)
    (hlimit : limit ≤ heap.size) (hiLimit : i < limit)
    (hunique : PairVecDivVHCPrefixUniqueExcept heap limit i)
    (hnewOnly : PairVecDivVHCPrefixOnlyAt heap limit newNode i)
    (hrun : pairVecDivVHCBubble i stop newNode heap = .ok heap') :
    PairVecDivVHCPrefixUnique heap' limit := by
  rw [pairVecDivVHCBubble] at hrun
  split at hrun <;> try contradiction
  next hi =>
    split at hrun <;> try contradiction
    next hstop =>
      split at hrun
      next heq =>
        simp only [Except.ok.injEq] at hrun
        subst heap'
        exact pairVecDivVHCSet_saved_prefixUnique heap limit newNode i hi
          hiLimit hunique hnewOnly
      next heq =>
        dsimp only at hrun
        split at hrun <;> try contradiction
        next hp =>
          have hipos : 0 < i := by omega
          have hparentLt : pairVecDivVHCParent i < i := by
            unfold pairVecDivVHCParent
            have hhalf : (i - 1) / 2 ≤ i - 1 := Nat.div_le_self _ _
            omega
          have hparentLimit : pairVecDivVHCParent i < limit :=
            Nat.lt_trans hparentLt hiLimit
          have hunique' := pairVecDivVHCSet_child_prefixUniqueExcept heap
            limit i (pairVecDivVHCParent i) hi hp hparentLimit
            (Nat.ne_of_gt hparentLt) hunique
          have hnewOnly' := pairVecDivVHCSet_child_prefixOnlyAt heap limit
            newNode i (pairVecDivVHCParent i) hi hp hparentLimit
            (Nat.ne_of_gt hparentLt) hnewOnly
          exact pairVecDivVHCBubble_prefixUnique (pairVecDivVHCParent i) stop
            newNode limit (heap.set i heap[pairVecDivVHCParent i]) heap'
            (by simpa using hlimit) hparentLimit hunique' hnewOnly' hrun
termination_by i
decreasing_by
  have hipos : 0 < i := by omega
  unfold pairVecDivVHCParent
  have hhalf : (i - 1) / 2 ≤ i - 1 := Nat.div_le_self _ _
  omega

theorem pairVecDivVHCBubbleBelow_prefixUnique
    (i anchor newNode limit : Nat) (heap heap' : Array Nat)
    (hlimit : limit ≤ heap.size) (hiLimit : i < limit)
    (hunique : PairVecDivVHCPrefixUniqueExcept heap limit i)
    (hnewOnly : PairVecDivVHCPrefixOnlyAt heap limit newNode i)
    (hrun : pairVecDivVHCBubbleBelow i anchor newNode heap = .ok heap') :
    PairVecDivVHCPrefixUnique heap' limit := by
  rw [pairVecDivVHCBubbleBelow] at hrun
  split at hrun <;> try contradiction
  next hi =>
    split at hrun <;> try contradiction
    next hpos =>
      dsimp only at hrun
      have hparentLt : pairVecDivVHCParent i < i := by
        unfold pairVecDivVHCParent
        have hhalf : (i - 1) / 2 ≤ i - 1 := Nat.div_le_self _ _
        omega
      have hparentLimit : pairVecDivVHCParent i < limit :=
        Nat.lt_trans hparentLt hiLimit
      split at hrun
      next heq =>
        simp only [Except.ok.injEq] at hrun
        subst heap'
        exact pairVecDivVHCSet_saved_prefixUnique heap limit newNode i hi
          hiLimit hunique hnewOnly
      next heq =>
        split at hrun <;> try contradiction
        next hp =>
          have hunique' := pairVecDivVHCSet_child_prefixUniqueExcept heap
            limit i (pairVecDivVHCParent i) hi hp hparentLimit
            (Nat.ne_of_gt hparentLt) hunique
          have hnewOnly' := pairVecDivVHCSet_child_prefixOnlyAt heap limit
            newNode i (pairVecDivVHCParent i) hi hp hparentLimit
            (Nat.ne_of_gt hparentLt) hnewOnly
          exact pairVecDivVHCBubbleBelow_prefixUnique
            (pairVecDivVHCParent i) anchor newNode limit
            (heap.set i heap[pairVecDivVHCParent i]) heap'
            (by simpa using hlimit) hparentLimit hunique' hnewOnly' hrun
termination_by i
decreasing_by
  unfold pairVecDivVHCParent
  have hhalf : (i - 1) / 2 ≤ i - 1 := Nat.div_le_self _ _
  omega

theorem pairVecDivVHCPush_prefixUniqueExcept
    (heap : Array Nat) (newNode : Nat)
    (hunique : ∀ (left right head : Nat), heap[left]? = some head →
      heap[right]? = some head → left = right)
    (hfresh : ∀ (slot : Nat), heap[slot]? ≠ some newNode) :
    PairVecDivVHCPrefixUniqueExcept (heap.push newNode)
      (heap.push newNode).size heap.size := by
  intro left right head hleft hright hleftGet hrightGet
  rw [Array.getElem?_push] at hleftGet hrightGet
  by_cases hleftLast : left = heap.size
  · subst left
    simp only [ite_true, Option.some.injEq] at hleftGet
    subst head
    by_cases hrightLast : right = heap.size
    · exact Or.inl hrightLast.symm
    · simp only [hrightLast, ↓reduceIte] at hrightGet
      exact False.elim (hfresh right hrightGet)
  · simp only [hleftLast, ↓reduceIte] at hleftGet
    by_cases hrightLast : right = heap.size
    · exact Or.inr (Or.inr hrightLast)
    · simp only [hrightLast, ↓reduceIte] at hrightGet
      exact Or.inl (hunique left right head hleftGet hrightGet)

theorem pairVecDivVHCPush_prefixOnlyAt (heap : Array Nat) (newNode : Nat)
    (hfresh : ∀ (slot : Nat), heap[slot]? ≠ some newNode) :
    PairVecDivVHCPrefixOnlyAt (heap.push newNode) (heap.push newNode).size
      newNode heap.size := by
  intro slot hslot hget
  rw [Array.getElem?_push] at hget
  by_cases hlast : slot = heap.size
  · exact hlast
  · simp only [hlast, ↓reduceIte] at hget
    exact False.elim (hfresh slot hget)

theorem pairVecDivVHCPush_preserves_unique
    (heap : Array Nat) (newNode : Nat)
    (hunique : ∀ (left right head : Nat), heap[left]? = some head →
      heap[right]? = some head → left = right)
    (hfresh : ∀ (slot : Nat), heap[slot]? ≠ some newNode) :
    ∀ (left right head : Nat), (heap.push newNode)[left]? = some head →
      (heap.push newNode)[right]? = some head → left = right := by
  intro left right head hleft hright
  rw [Array.getElem?_push] at hleft hright
  by_cases hleftLast : left = heap.size
  · subst left
    simp only [ite_true, Option.some.injEq] at hleft
    subst head
    by_cases hrightLast : right = heap.size
    · exact hrightLast.symm
    · simp only [hrightLast, ↓reduceIte] at hright
      exact False.elim (hfresh right hright)
  · simp only [hleftLast, ↓reduceIte] at hleft
    by_cases hrightLast : right = heap.size
    · subst right
      simp only [ite_true, Option.some.injEq] at hright
      subst head
      exact False.elim (hfresh left hleft)
    · simp only [hrightLast, ↓reduceIte] at hright
      exact hunique left right head hleft hright

theorem pairVecDivVHCSet_fresh_preserves_unique
    (heap : Array Nat) (slot newNode : Nat) (hslot : slot < heap.size)
    (hunique : ∀ (left right head : Nat), heap[left]? = some head →
      heap[right]? = some head → left = right)
    (hfresh : ∀ (i : Nat), heap[i]? ≠ some newNode) :
    ∀ (left right head : Nat), (heap.set slot newNode)[left]? = some head →
      (heap.set slot newNode)[right]? = some head → left = right := by
  intro left right head hleft hright
  by_cases hleftSlot : slot = left
  · subst left
    rw [Array.getElem?_set_self hslot] at hleft
    simp only [Option.some.injEq] at hleft
    subst head
    by_cases hrightSlot : slot = right
    · exact hrightSlot
    · rw [Array.getElem?_set_ne hslot hrightSlot] at hright
      exact False.elim (hfresh right hright)
  · rw [Array.getElem?_set_ne hslot hleftSlot] at hleft
    by_cases hrightSlot : slot = right
    · subst right
      rw [Array.getElem?_set_self hslot] at hright
      simp only [Option.some.injEq] at hright
      subst head
      exact False.elim (hfresh left hleft)
    · rw [Array.getElem?_set_ne hslot hrightSlot] at hright
      exact hunique left right head hleft hright

theorem pairVecDivVHCBubble_push_preserves_unique
    (stop newNode : Nat) (heap heap' : Array Nat)
    (hunique : ∀ (left right head : Nat), heap[left]? = some head →
      heap[right]? = some head → left = right)
    (hfresh : ∀ (slot : Nat), heap[slot]? ≠ some newNode)
    (hrun : pairVecDivVHCBubble heap.size stop newNode (heap.push newNode) =
      .ok heap') :
    ∀ (left right head : Nat), heap'[left]? = some head →
      heap'[right]? = some head → left = right := by
  have hprefix := pairVecDivVHCBubble_prefixUnique heap.size stop newNode
    (heap.push newNode).size (heap.push newNode) heap' (by omega) (by simp)
    (pairVecDivVHCPush_prefixUniqueExcept heap newNode hunique hfresh)
    (pairVecDivVHCPush_prefixOnlyAt heap newNode hfresh) hrun
  have hsize := pairVecDivVHCBubble_size heap.size stop newNode
    (heap.push newNode) heap' hrun
  intro left right head hleft hright
  have hleftBound : left < heap'.size := by
    by_contra hnot
    rw [Array.getElem?_eq_none (by omega)] at hleft
    contradiction
  have hrightBound : right < heap'.size := by
    by_contra hnot
    rw [Array.getElem?_eq_none (by omega)] at hright
    contradiction
  exact hprefix left right head (by simpa [hsize] using hleftBound)
    (by simpa [hsize] using hrightBound) hleft hright

theorem pairVecDivVHCBubbleBelow_push_preserves_unique
    (anchor newNode : Nat) (heap heap' : Array Nat)
    (hunique : ∀ (left right head : Nat), heap[left]? = some head →
      heap[right]? = some head → left = right)
    (hfresh : ∀ (slot : Nat), heap[slot]? ≠ some newNode)
    (hrun : pairVecDivVHCBubbleBelow heap.size anchor newNode
      (heap.push newNode) = .ok heap') :
    ∀ (left right head : Nat), heap'[left]? = some head →
      heap'[right]? = some head → left = right := by
  have hprefix := pairVecDivVHCBubbleBelow_prefixUnique heap.size anchor newNode
    (heap.push newNode).size (heap.push newNode) heap' (by omega) (by simp)
    (pairVecDivVHCPush_prefixUniqueExcept heap newNode hunique hfresh)
    (pairVecDivVHCPush_prefixOnlyAt heap newNode hfresh) hrun
  have hsize := pairVecDivVHCBubbleBelow_size heap.size anchor newNode
    (heap.push newNode) heap' hrun
  intro left right head hleft hright
  have hleftBound : left < heap'.size := by
    by_contra hnot
    rw [Array.getElem?_eq_none (by omega)] at hleft
    contradiction
  have hrightBound : right < heap'.size := by
    by_contra hnot
    rw [Array.getElem?_eq_none (by omega)] at hright
    contradiction
  exact hprefix left right head (by simpa [hsize] using hleftBound)
    (by simpa [hsize] using hrightBound) hleft hright

theorem pairVecDivVHCSet_excludes_of_onlyAt
    (heap : Array Nat) (nodeIndex slot replacement : Nat)
    (hslot : slot < heap.size)
    (honly : PairVecDivVHCOnlyAt heap nodeIndex slot)
    (hne : replacement ≠ nodeIndex) :
    PairVecDivVHCExcludes (heap.set slot replacement) nodeIndex := by
  intro i hi
  by_cases hislot : slot = i
  · subst i
    have hvalue : replacement = nodeIndex := by
      simpa [Array.getElem?_set, hslot] using hi
    exact hne hvalue
  · rw [Array.getElem?_set_ne hslot hislot] at hi
    exact hislot (honly i hi).symm

theorem pairVecDivVHCExcludes_onlyAt (heap : Array Nat) (nodeIndex slot : Nat)
    (hexcludes : PairVecDivVHCExcludes heap nodeIndex) :
    PairVecDivVHCOnlyAt heap nodeIndex slot := by
  intro i hi
  exact False.elim (hexcludes i hi)

theorem pairVecDivVHCSiftDown_valuesFrom
    (i child limit lastNode : Nat) (heap heap' source : Array Nat)
    (nodes : Array PairVecDivVHCNode)
    (hfrom : PairVecDivVHCValuesFrom heap source)
    (hlast : ∃ j : Nat, source[j]? = some lastNode)
    (hrun : pairVecDivVHCSiftDown i child limit lastNode heap nodes =
      .ok heap') :
    PairVecDivVHCValuesFrom heap' source := by
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
            cases hlastMono : pairVecDivVHCMono lastNode nodes with
            | error fault => simp [hleft, hright, hlastMono] at hrun
            | ok lastMono =>
              simp only [hleft, hright, hlastMono] at hrun
              have hleftFrom := hfrom child heap[child] (by
                rw [Array.getElem?_eq_getElem (by omega)])
              have hrightFrom := hfrom (child + 1) heap[child + 1] (by
                rw [Array.getElem?_eq_getElem (by omega)])
              by_cases hselected : leftMono.deg > rightMono.deg
              · simp only [hselected, ↓reduceIte] at hrun
                split at hrun
                next hgreater =>
                  exact pairVecDivVHCSiftDown_valuesFrom child
                    (child * 2 + 1) limit lastNode (heap.set i heap[child])
                    heap' source nodes
                    (pairVecDivVHCSet_valuesFrom heap source i heap[child] hi
                      hfrom hleftFrom)
                    hlast hrun
                next hgreater =>
                  simp only [Except.ok.injEq] at hrun
                  subst heap'
                  exact pairVecDivVHCSet_valuesFrom heap source i lastNode hi
                    hfrom hlast
              · simp only [hselected, ↓reduceIte] at hrun
                split at hrun
                next hgreater =>
                  exact pairVecDivVHCSiftDown_valuesFrom (child + 1)
                    ((child + 1) * 2 + 1) limit lastNode
                    (heap.set i heap[child + 1]) heap' source nodes
                    (pairVecDivVHCSet_valuesFrom heap source i heap[child + 1]
                      hi hfrom hrightFrom)
                    hlast hrun
                next hgreater =>
                  simp only [Except.ok.injEq] at hrun
                  subst heap'
                  exact pairVecDivVHCSet_valuesFrom heap source i lastNode hi
                    hfrom hlast
      next hchild =>
        simp only [Except.ok.injEq] at hrun
        subst heap'
        exact pairVecDivVHCSet_valuesFrom heap source i lastNode hi hfrom hlast
termination_by limit - child
decreasing_by
  all_goals omega

theorem pairVecDivVHCSiftDown_prefixUnique
    (i child limit lastNode : Nat) (heap heap' : Array Nat)
    (nodes : Array PairVecDivVHCNode)
    (hlimit : limit ≤ heap.size) (hiLimit : i < limit)
    (hichild : i < child)
    (hunique : PairVecDivVHCPrefixUniqueExcept heap limit i)
    (hsaved : PairVecDivVHCPrefixOnlyAt heap limit lastNode i)
    (hlastSlot : heap[limit]? = some lastNode)
    (hrun : pairVecDivVHCSiftDown i child limit lastNode heap nodes =
      .ok heap') :
    PairVecDivVHCPrefixUnique heap' limit := by
  rw [pairVecDivVHCSiftDown] at hrun
  split at hrun <;> try contradiction
  next hi =>
    split at hrun <;> try contradiction
    next hlimitGuard =>
      split at hrun
      next hchildGuard =>
        dsimp only at hrun
        cases hleft : pairVecDivVHCMono heap[child] nodes with
        | error fault => simp [hleft] at hrun
        | ok leftMono =>
          cases hright : pairVecDivVHCMono heap[child + 1] nodes with
          | error fault => simp [hleft, hright] at hrun
          | ok rightMono =>
            cases hlastMono : pairVecDivVHCMono lastNode nodes with
            | error fault => simp [hleft, hright, hlastMono] at hrun
            | ok lastMono =>
              simp only [hleft, hright, hlastMono] at hrun
              have hleftBound : child < heap.size := by omega
              have hrightBound : child + 1 < heap.size := by omega
              have hleftLimit : child < limit := hchildGuard
              have hiLeft : i ≠ child := Nat.ne_of_lt hichild
              have hleftUnique := pairVecDivVHCSet_child_prefixUniqueExcept heap
                limit i child hi hleftBound hleftLimit hiLeft hunique
              have hleftSaved := pairVecDivVHCSet_child_prefixOnlyAt heap limit
                lastNode i child hi hleftBound hleftLimit hiLeft hsaved
              have hleftLast : (heap.set i heap[child])[limit]? =
                  some lastNode := by
                rw [Array.getElem?_set_ne hi (by omega)]
                exact hlastSlot
              by_cases hselected : leftMono.deg > rightMono.deg
              · simp only [hselected, ↓reduceIte] at hrun
                split at hrun
                next hgreater =>
                  exact pairVecDivVHCSiftDown_prefixUnique child
                    (child * 2 + 1) limit lastNode (heap.set i heap[child])
                    heap' nodes (by simpa using hlimit) hleftLimit (by omega)
                    hleftUnique hleftSaved hleftLast hrun
                next hgreater =>
                  simp only [Except.ok.injEq] at hrun
                  subst heap'
                  exact pairVecDivVHCSet_saved_prefixUnique heap limit lastNode i
                    hi hiLimit hunique hsaved
              · simp only [hselected, ↓reduceIte] at hrun
                split at hrun
                next hgreater =>
                  have hrightLimit : child + 1 < limit := by
                    by_contra hnot
                    have heq : child + 1 = limit := by omega
                    have hrightNode : heap[child + 1] = lastNode := by
                      have hlastGet : heap[limit] = lastNode := by
                        rw [Array.getElem?_eq_getElem hlimitGuard] at hlastSlot
                        exact Option.some.inj hlastSlot
                      simpa [heq] using hlastGet
                    have hmonoEq : rightMono = lastMono := by
                      rw [hrightNode] at hright
                      rw [hright] at hlastMono
                      exact Except.ok.inj hlastMono
                    subst rightMono
                    omega
                  have hiRight : i ≠ child + 1 := by omega
                  have hrightUnique :=
                    pairVecDivVHCSet_child_prefixUniqueExcept heap limit i
                      (child + 1) hi hrightBound hrightLimit hiRight hunique
                  have hrightSaved := pairVecDivVHCSet_child_prefixOnlyAt heap
                    limit lastNode i (child + 1) hi hrightBound hrightLimit
                    hiRight hsaved
                  have hrightLast : (heap.set i heap[child + 1])[limit]? =
                      some lastNode := by
                    rw [Array.getElem?_set_ne hi (by omega)]
                    exact hlastSlot
                  exact pairVecDivVHCSiftDown_prefixUnique (child + 1)
                    ((child + 1) * 2 + 1) limit lastNode
                    (heap.set i heap[child + 1]) heap' nodes
                    (by simpa using hlimit) hrightLimit (by omega)
                    hrightUnique hrightSaved hrightLast hrun
                next hgreater =>
                  simp only [Except.ok.injEq] at hrun
                  subst heap'
                  exact pairVecDivVHCSet_saved_prefixUnique heap limit lastNode i
                    hi hiLimit hunique hsaved
      next hchildGuard =>
        simp only [Except.ok.injEq] at hrun
        subst heap'
        exact pairVecDivVHCSet_saved_prefixUnique heap limit lastNode i hi
          hiLimit hunique hsaved
termination_by limit - child
decreasing_by
  all_goals omega

/-- `sift_down` removes the distinguished value at its current hole.  The
proof follows the real recursive writes: a child is copied into the hole and
the saved last node is eventually installed, neither of which can be the
removed root under heap-head uniqueness. -/
theorem pairVecDivVHCSiftDown_excludes_of_onlyAt
    (i child limit lastNode root : Nat) (heap heap' : Array Nat)
    (nodes : Array PairVecDivVHCNode)
    (hichild : i < child)
    (honly : PairVecDivVHCOnlyAt heap root i)
    (hlast : lastNode ≠ root)
    (hrun : pairVecDivVHCSiftDown i child limit lastNode heap nodes =
      .ok heap') :
    PairVecDivVHCExcludes heap' root := by
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
            cases hlastMono : pairVecDivVHCMono lastNode nodes with
            | error fault => simp [hleft, hright, hlastMono] at hrun
            | ok lastMono =>
              simp only [hleft, hright, hlastMono] at hrun
              have hchildNe : heap[child] ≠ root := by
                intro heq
                have hget : heap[child]? = some root := by
                  rw [Array.getElem?_eq_getElem (by omega), heq]
                exact (Nat.ne_of_lt hichild) (honly child hget).symm
              have hchildSet := pairVecDivVHCSet_excludes_of_onlyAt heap root i
                heap[child] hi honly hchildNe
              by_cases hselected : leftMono.deg > rightMono.deg
              · simp only [hselected, ↓reduceIte] at hrun
                split at hrun
                next hgreater =>
                  exact pairVecDivVHCSiftDown_excludes_of_onlyAt child
                    (child * 2 + 1) limit lastNode root
                    (heap.set i heap[child]) heap' nodes (by omega)
                    (pairVecDivVHCExcludes_onlyAt _ _ _ hchildSet) hlast hrun
                next hgreater =>
                  simp only [Except.ok.injEq] at hrun
                  subst heap'
                  exact pairVecDivVHCSet_excludes_of_onlyAt heap root i
                    lastNode hi honly hlast
              · simp only [hselected, ↓reduceIte] at hrun
                have hrightNe : heap[child + 1] ≠ root := by
                  intro heq
                  have hget : heap[child + 1]? = some root := by
                    rw [Array.getElem?_eq_getElem (by omega), heq]
                  exact (by omega : child + 1 ≠ i) (honly (child + 1) hget)
                have hrightSet := pairVecDivVHCSet_excludes_of_onlyAt heap root
                  i heap[child + 1] hi honly hrightNe
                split at hrun
                next hgreater =>
                  exact pairVecDivVHCSiftDown_excludes_of_onlyAt (child + 1)
                    ((child + 1) * 2 + 1) limit lastNode root
                    (heap.set i heap[child + 1]) heap' nodes (by omega)
                    (pairVecDivVHCExcludes_onlyAt _ _ _ hrightSet) hlast hrun
                next hgreater =>
                  simp only [Except.ok.injEq] at hrun
                  subst heap'
                  exact pairVecDivVHCSet_excludes_of_onlyAt heap root i
                    lastNode hi honly hlast
      next hchild =>
        simp only [Except.ok.injEq] at hrun
        subst heap'
        exact pairVecDivVHCSet_excludes_of_onlyAt heap root i lastNode hi honly
          hlast
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

theorem pairVecDivVHCExtract_two_root (heap heap' : Array Nat)
    (nodes : Array PairVecDivVHCNode) (hsize : heap.size = 2)
    (hrun : pairVecDivVHCExtract heap nodes = .ok heap') :
    heap'[0]? = some heap[1] := by
  have hnonempty : 0 < heap.size := by omega
  have hlimit : heap.size - 1 < heap.size := by omega
  have hchild : ¬1 < heap.size - 1 := by omega
  unfold pairVecDivVHCExtract at hrun
  simp only [hnonempty, ↓reduceDIte] at hrun
  rw [pairVecDivVHCSiftDown] at hrun
  simp only [hnonempty, hlimit, hchild, ↓reduceDIte,
    Except.ok.injEq] at hrun
  subst heap'
  rw [Array.getElem?_pop]
  simp only [Array.size_set, hsize, Nat.reduceSub, Nat.reduceLT,
    ↓reduceDIte, Array.getElem?_set_self hnonempty]
  simp

theorem pairVecDivVHCExtract_root_dominates_candidates
    (heap heap' : Array Nat) (nodes : Array PairVecDivVHCNode)
    (leftMono rightMono lastMono : UMonomial) (hsize : 2 < heap.size)
    (hleft : pairVecDivVHCMono heap[1] nodes = .ok leftMono)
    (hright : pairVecDivVHCMono heap[2] nodes = .ok rightMono)
    (hlast : pairVecDivVHCMono heap[heap.size - 1] nodes = .ok lastMono)
    (hrun : pairVecDivVHCExtract heap nodes = .ok heap') :
    ∃ rootHead rootMono,
      heap'[0]? = some rootHead ∧
        pairVecDivVHCMono rootHead nodes = .ok rootMono ∧
        leftMono.deg ≤ rootMono.deg ∧ rightMono.deg ≤ rootMono.deg ∧
        lastMono.deg ≤ rootMono.deg := by
  have hnonempty : 0 < heap.size := by omega
  unfold pairVecDivVHCExtract at hrun
  simp only [hnonempty, ↓reduceDIte] at hrun
  cases hsift : pairVecDivVHCSiftDown 0 1 (heap.size - 1)
      heap[heap.size - 1] heap nodes with
  | error fault => simp [hsift] at hrun
  | ok shifted =>
      rw [hsift] at hrun
      simp only [Except.ok.injEq] at hrun
      subst heap'
      rcases pairVecDivVHCSiftDown_root_dominates_candidates heap shifted
          nodes leftMono rightMono lastMono hsize hleft hright hlast hsift with
        ⟨rootHead, rootMono, hroot, hrootMono, hleftLe, hrightLe, hlastLe⟩
      refine ⟨rootHead, rootMono, ?_, hrootMono, hleftLe, hrightLe, hlastLe⟩
      rw [Array.getElem?_pop]
      have hshiftSize := pairVecDivVHCSiftDown_size 0 1 (heap.size - 1)
        heap[heap.size - 1] heap shifted nodes hsift
      have hpop : 0 < heap.size - 1 := by omega
      simp only [hshiftSize, hpop, ↓reduceDIte]
      exact hroot

theorem pairVecDivVHCExtract_valuesFrom (heap heap' : Array Nat)
    (nodes : Array PairVecDivVHCNode)
    (hrun : pairVecDivVHCExtract heap nodes = .ok heap') :
    PairVecDivVHCValuesFrom heap' heap := by
  unfold pairVecDivVHCExtract at hrun
  split at hrun <;> try contradiction
  next hnonempty =>
    dsimp only at hrun
    cases hsift : pairVecDivVHCSiftDown 0 1 (heap.size - 1)
        heap[heap.size - 1] heap nodes with
    | error fault => simp [hsift] at hrun
    | ok shifted =>
        rw [hsift] at hrun
        simp only [Except.ok.injEq] at hrun
        subst heap'
        have hlastBound : heap.size - 1 < heap.size := by omega
        have hshifted := pairVecDivVHCSiftDown_valuesFrom 0 1
          (heap.size - 1) heap[heap.size - 1] heap shifted heap nodes
          (pairVecDivVHCValuesFrom_refl heap)
          ⟨heap.size - 1, by
            rw [Array.getElem?_eq_getElem hlastBound]⟩ hsift
        intro i value hi
        rw [Array.getElem?_pop] at hi
        split at hi
        · exact hshifted i value hi
        · contradiction

/-- Extracting a unique heap root removes that node index from the whole
active heap, not merely from slot zero. -/
theorem pairVecDivVHCExtract_excludes_root (heap heap' : Array Nat)
    (nodes : Array PairVecDivVHCNode)
    (hnonempty : 0 < heap.size)
    (hunique : PairVecDivVHCOnlyAt heap heap[0] 0)
    (hrun : pairVecDivVHCExtract heap nodes = .ok heap') :
    PairVecDivVHCExcludes heap' heap[0] := by
  have hsize := pairVecDivVHCExtract_size heap heap' nodes hrun
  by_cases hone : heap.size = 1
  · intro i
    have hempty : heap'.size = 0 := by omega
    rw [Array.getElem?_eq_none (by omega)]
    simp
  · unfold pairVecDivVHCExtract at hrun
    simp only [hnonempty, ↓reduceDIte] at hrun
    cases hsift : pairVecDivVHCSiftDown 0 1 (heap.size - 1)
        heap[heap.size - 1] heap nodes with
    | error fault => simp [hsift] at hrun
    | ok shifted =>
        rw [hsift] at hrun
        simp only [Except.ok.injEq] at hrun
        subst heap'
        have hlastNe : heap[heap.size - 1] ≠ heap[0] := by
          intro heq
          have hlastBound : heap.size - 1 < heap.size := by omega
          have hget : heap[heap.size - 1]? = some heap[0] := by
            rw [Array.getElem?_eq_getElem hlastBound, heq]
          have := hunique (heap.size - 1) hget
          omega
        have hexcludes := pairVecDivVHCSiftDown_excludes_of_onlyAt 0 1
          (heap.size - 1) heap[heap.size - 1] heap[0] heap shifted nodes
          (by omega) hunique hlastNe hsift
        intro i hi
        rw [Array.getElem?_pop] at hi
        split at hi
        · exact hexcludes i hi
        · contradiction

/-- The pointer-copy sift followed by dropping the saved sentinel preserves
uniqueness of every surviving heap head. -/
theorem pairVecDivVHCExtract_preserves_unique (heap heap' : Array Nat)
    (nodes : Array PairVecDivVHCNode)
    (hunique : ∀ (left right head : Nat), heap[left]? = some head →
      heap[right]? = some head → left = right)
    (hrun : pairVecDivVHCExtract heap nodes = .ok heap') :
    ∀ (left right head : Nat), heap'[left]? = some head →
      heap'[right]? = some head → left = right := by
  unfold pairVecDivVHCExtract at hrun
  split at hrun <;> try contradiction
  next hnonempty =>
    dsimp only at hrun
    cases hsift : pairVecDivVHCSiftDown 0 1 (heap.size - 1)
        heap[heap.size - 1] heap nodes with
    | error fault => simp [hsift] at hrun
    | ok shifted =>
        rw [hsift] at hrun
        simp only [Except.ok.injEq] at hrun
        subst heap'
        have hshiftSize := pairVecDivVHCSiftDown_size 0 1 (heap.size - 1)
          heap[heap.size - 1] heap shifted nodes hsift
        by_cases hone : heap.size = 1
        · intro left right head hleft hright
          rw [Array.getElem?_pop] at hleft
          split at hleft <;> try contradiction
          next hleftBound =>
            rw [hshiftSize, hone] at hleftBound
            omega
        · have hlimitBound : heap.size - 1 < heap.size := by omega
          have hprefixExcept : PairVecDivVHCPrefixUniqueExcept heap
              (heap.size - 1) 0 := by
            intro left right head hleft hright hleftGet hrightGet
            exact Or.inl (hunique left right head hleftGet hrightGet)
          have hsaved : PairVecDivVHCPrefixOnlyAt heap (heap.size - 1)
              heap[heap.size - 1] 0 := by
            intro slot hslot hget
            have hlastGet : heap[heap.size - 1]? =
                some heap[heap.size - 1] := by
              rw [Array.getElem?_eq_getElem hlimitBound]
            have heq := hunique slot (heap.size - 1) heap[heap.size - 1]
              hget hlastGet
            omega
          have hprefix := pairVecDivVHCSiftDown_prefixUnique 0 1
            (heap.size - 1) heap[heap.size - 1] heap shifted nodes (by omega)
            (by omega) (by omega) hprefixExcept hsaved
            (by rw [Array.getElem?_eq_getElem hlimitBound]) hsift
          intro left right head hleft hright
          rw [Array.getElem?_pop] at hleft hright
          split at hleft <;> try contradiction
          next hleftBound =>
            split at hright <;> try contradiction
            next hrightBound =>
              exact hprefix left right head
                (by simpa [hshiftSize] using hleftBound)
                (by simpa [hshiftSize] using hrightBound) hleft hright

theorem pairVecDivVHCHeapUnique_toList_nodup (heap : Array Nat)
    (hunique : ∀ (left right head : Nat), heap[left]? = some head →
      heap[right]? = some head → left = right) :
    heap.toList.Nodup := by
  rw [List.nodup_iff_injective_getElem]
  intro left right heq
  have hleft : left.val < heap.size := by simpa using left.isLt
  have hright : right.val < heap.size := by simpa using right.isLt
  have hvalueEq : heap[left.val] = heap[right.val] := by
    rw [← Array.getElem_toList hleft, ← Array.getElem_toList hright]
    exact heq
  have hslots := hunique left.val right.val heap[left.val]
    (Array.getElem?_eq_getElem hleft) (by
      rw [Array.getElem?_eq_getElem hright, ← hvalueEq])
  exact Fin.ext hslots

/-- A same-size heap rewrite that only copies source values and preserves slot
uniqueness is a permutation of the concrete heads.  Hence every source head
has a concrete target slot; this supplies the reverse direction absent from
`PairVecDivVHCValuesFrom`. -/
theorem pairVecDivVHCValuesFrom_preserves_every_head
    (target source : Array Nat)
    (hfrom : PairVecDivVHCValuesFrom target source)
    (hsize : target.size = source.size)
    (hsourceUnique : ∀ (left right head : Nat),
      source[left]? = some head → source[right]? = some head → left = right)
    (htargetUnique : ∀ (left right head : Nat),
      target[left]? = some head → target[right]? = some head → left = right) :
    ∀ (slot head : Nat), source[slot]? = some head →
      ∃ targetSlot : Nat, target[targetSlot]? = some head := by
  have hsourceNodup :=
    pairVecDivVHCHeapUnique_toList_nodup source hsourceUnique
  have htargetNodup :=
    pairVecDivVHCHeapUnique_toList_nodup target htargetUnique
  have hsubset : target.toList.toFinset ⊆ source.toList.toFinset := by
    intro head hhead
    simp only [List.mem_toFinset] at hhead ⊢
    rcases List.getElem_of_mem hhead with ⟨slot, hslotList, hslotValue⟩
    have hslot : slot < target.size := by simpa using hslotList
    have hget : target[slot]? = some head := by
      rw [Array.getElem?_eq_getElem hslot, ← Array.getElem_toList hslot]
      exact congrArg some hslotValue
    rcases hfrom slot head hget with ⟨sourceSlot, hsourceGet⟩
    have hsourceBound : sourceSlot < source.size := by
      by_contra hnot
      rw [Array.getElem?_eq_none (by omega)] at hsourceGet
      contradiction
    rw [Array.getElem?_eq_getElem hsourceBound] at hsourceGet
    have hsourceEq := Option.some.inj hsourceGet
    rw [← hsourceEq]
    exact Array.getElem_mem_toList hsourceBound
  have htargetCard : target.toList.toFinset.card = target.size := by
    rw [List.toFinset_card_of_nodup htargetNodup]
    rfl
  have hsourceCard : source.toList.toFinset.card = source.size := by
    rw [List.toFinset_card_of_nodup hsourceNodup]
    rfl
  have hsets : target.toList.toFinset = source.toList.toFinset := by
    apply Finset.eq_of_subset_of_card_le hsubset
    rw [htargetCard, hsourceCard, hsize]
  intro slot head hget
  have hslot : slot < source.size := by
    by_contra hnot
    rw [Array.getElem?_eq_none (by omega)] at hget
    contradiction
  have hheadSource : head ∈ source.toList.toFinset := by
    simp only [List.mem_toFinset]
    rw [Array.getElem?_eq_getElem hslot] at hget
    have hheadEq := Option.some.inj hget
    rw [← hheadEq]
    exact Array.getElem_mem_toList hslot
  have hheadTarget : head ∈ target.toList.toFinset := by
    rw [hsets]
    exact hheadSource
  simp only [List.mem_toFinset] at hheadTarget
  rcases List.getElem_of_mem hheadTarget with
    ⟨targetSlot, htargetBound, htargetValue⟩
  refine ⟨targetSlot, ?_⟩
  have htargetSlot : targetSlot < target.size := by simpa using htargetBound
  rw [Array.getElem?_eq_getElem htargetSlot,
    ← Array.getElem_toList htargetSlot]
  exact congrArg some htargetValue

theorem pairVecDivVHCBubble_push_preserves_every_head
    (stop newNode : Nat) (heap heap' : Array Nat)
    (hunique : ∀ (left right head : Nat), heap[left]? = some head →
      heap[right]? = some head → left = right)
    (hfresh : ∀ (slot : Nat), heap[slot]? ≠ some newNode)
    (hrun : pairVecDivVHCBubble heap.size stop newNode (heap.push newNode) =
      .ok heap') :
    ∀ (slot head : Nat), (heap.push newNode)[slot]? = some head →
      ∃ targetSlot : Nat, heap'[targetSlot]? = some head := by
  have hfrom := pairVecDivVHCBubble_valuesFrom heap.size stop newNode
    (heap.push newNode) heap' (heap.push newNode)
    (pairVecDivVHCValuesFrom_refl _) ⟨heap.size, by simp⟩ hrun
  exact pairVecDivVHCValuesFrom_preserves_every_head heap'
    (heap.push newNode) hfrom
    (pairVecDivVHCBubble_size heap.size stop newNode (heap.push newNode)
      heap' hrun)
    (pairVecDivVHCPush_preserves_unique heap newNode hunique hfresh)
    (pairVecDivVHCBubble_push_preserves_unique stop newNode heap heap'
      hunique hfresh hrun)

theorem pairVecDivVHCBubbleBelow_push_preserves_every_head
    (anchor newNode : Nat) (heap heap' : Array Nat)
    (hunique : ∀ (left right head : Nat), heap[left]? = some head →
      heap[right]? = some head → left = right)
    (hfresh : ∀ (slot : Nat), heap[slot]? ≠ some newNode)
    (hrun : pairVecDivVHCBubbleBelow heap.size anchor newNode
      (heap.push newNode) = .ok heap') :
    ∀ (slot head : Nat), (heap.push newNode)[slot]? = some head →
      ∃ targetSlot : Nat, heap'[targetSlot]? = some head := by
  have hfrom := pairVecDivVHCBubbleBelow_valuesFrom heap.size anchor newNode
    (heap.push newNode) heap' (heap.push newNode)
    (pairVecDivVHCValuesFrom_refl _) ⟨heap.size, by simp⟩ hrun
  exact pairVecDivVHCValuesFrom_preserves_every_head heap'
    (heap.push newNode) hfrom
    (pairVecDivVHCBubbleBelow_size heap.size anchor newNode
      (heap.push newNode) heap' hrun)
    (pairVecDivVHCPush_preserves_unique heap newNode hunique hfresh)
    (pairVecDivVHCBubbleBelow_push_preserves_unique anchor newNode heap heap'
      hunique hfresh hrun)

/-- The checked extract removes exactly the unique old root: every other old
heap head occurs at some concrete slot of the returned heap.  This is the
reverse direction missing from `pairVecDivVHCExtract_valuesFrom`. -/
theorem pairVecDivVHCExtract_preserves_nonroot_head
    (heap heap' : Array Nat) (nodes : Array PairVecDivVHCNode)
    (hnonempty : 0 < heap.size)
    (hunique : ∀ (left right head : Nat), heap[left]? = some head →
      heap[right]? = some head → left = right)
    (hrun : pairVecDivVHCExtract heap nodes = .ok heap') :
    ∀ (slot head : Nat), heap[slot]? = some head → head ≠ heap[0] →
      ∃ newSlot : Nat, heap'[newSlot]? = some head := by
  have hsize := pairVecDivVHCExtract_size heap heap' nodes hrun
  have hfrom := pairVecDivVHCExtract_valuesFrom heap heap' nodes hrun
  have hrootOnly : PairVecDivVHCOnlyAt heap heap[0] 0 := by
    intro slot hget
    exact hunique slot 0 heap[0] hget
      (Array.getElem?_eq_getElem hnonempty)
  have hexcludes := pairVecDivVHCExtract_excludes_root heap heap' nodes
    hnonempty hrootOnly hrun
  have hunique' := pairVecDivVHCExtract_preserves_unique heap heap' nodes
    hunique hrun
  have hnodup := pairVecDivVHCHeapUnique_toList_nodup heap hunique
  have hnodup' := pairVecDivVHCHeapUnique_toList_nodup heap' hunique'
  have hrootMem : heap[0] ∈ heap.toList.toFinset := by
    simp only [List.mem_toFinset]
    exact Array.getElem_mem_toList hnonempty
  have hsubset : heap'.toList.toFinset ⊆
      heap.toList.toFinset.erase heap[0] := by
    intro head hhead
    simp only [List.mem_toFinset] at hhead
    rcases List.getElem_of_mem hhead with ⟨slot, hslotList, hslotValue⟩
    have hslot : slot < heap'.size := by simpa using hslotList
    have hget : heap'[slot]? = some head := by
      rw [Array.getElem?_eq_getElem hslot]
      rw [← Array.getElem_toList hslot]
      exact congrArg some hslotValue
    rcases hfrom slot head hget with ⟨oldSlot, holdGet⟩
    apply Finset.mem_erase.mpr
    refine ⟨?_, ?_⟩
    · intro heq
      apply hexcludes slot
      rw [← heq]
      exact hget
    · simp only [List.mem_toFinset]
      have holdBound : oldSlot < heap.size := by
        by_contra hnot
        rw [Array.getElem?_eq_none (by omega)] at holdGet
        contradiction
      rw [Array.getElem?_eq_getElem holdBound] at holdGet
      have holdEq := Option.some.inj holdGet
      rw [← holdEq]
      exact Array.getElem_mem_toList holdBound
  have htargetCard : heap'.toList.toFinset.card = heap'.size := by
    rw [List.toFinset_card_of_nodup hnodup']
    rfl
  have hsourceCard : (heap.toList.toFinset.erase heap[0]).card =
      heap.size - 1 := by
    rw [Finset.card_erase_of_mem hrootMem,
      List.toFinset_card_of_nodup hnodup]
    rfl
  have hsets : heap'.toList.toFinset =
      heap.toList.toFinset.erase heap[0] := by
    apply Finset.eq_of_subset_of_card_le hsubset
    rw [hsourceCard, htargetCard]
    omega
  intro slot head hget hne
  have hslot : slot < heap.size := by
    by_contra hnot
    rw [Array.getElem?_eq_none (by omega)] at hget
    contradiction
  have hheadOld : head ∈ heap.toList.toFinset := by
    simp only [List.mem_toFinset]
    rw [Array.getElem?_eq_getElem hslot] at hget
    have hheadEq := Option.some.inj hget
    rw [← hheadEq]
    exact Array.getElem_mem_toList hslot
  have hheadNew : head ∈ heap'.toList.toFinset := by
    rw [hsets]
    exact Finset.mem_erase.mpr ⟨hne, hheadOld⟩
  simp only [List.mem_toFinset] at hheadNew
  rcases List.getElem_of_mem hheadNew with ⟨newSlot, hnewBound, hnewValue⟩
  refine ⟨newSlot, ?_⟩
  have hnew : newSlot < heap'.size := by simpa using hnewBound
  rw [Array.getElem?_eq_getElem hnew]
  rw [← Array.getElem_toList hnew]
  exact congrArg some hnewValue

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

theorem PairVecDivVHCNodeDenotes.quotientIndex_lt
    (quotient divisor : SparsePolyZp) (node : PairVecDivVHCNode)
    (hdenotes : PairVecDivVHCNodeDenotes quotient divisor node) :
    node.quotientIndex < quotient.size := by
  rcases hdenotes with ⟨quotientTerm, _, hquotient, _, _⟩
  by_contra hnot
  rw [Array.getElem?_eq_none (by omega)] at hquotient
  contradiction

theorem PairVecDivVHCNodeDenotes.product_degree_eq
    (quotient divisor : SparsePolyZp) (node : PairVecDivVHCNode)
    (degree : Nat) (hdenotes : PairVecDivVHCNodeDenotes quotient divisor node)
    (hmono : node.mono = some ⟨degree⟩)
    (quotientTerm divisorTerm : UMonomial × Zp)
    (hquotient : quotient[node.quotientIndex]? = some quotientTerm)
    (hdivisor : divisor[node.divisorIndex]? = some divisorTerm) :
    quotientTerm.1.deg + divisorTerm.1.deg = degree := by
  rcases hdenotes with
    ⟨storedQuotient, storedDivisor, hstoredQuotient, hstoredDivisor,
      hstoredMono⟩
  rw [hquotient] at hstoredQuotient
  rw [hdivisor] at hstoredDivisor
  simp only [Option.some.injEq] at hstoredQuotient hstoredDivisor
  subst storedQuotient
  subst storedDivisor
  rw [hmono] at hstoredMono
  exact congrArg UMonomial.deg (Option.some.inj hstoredMono).symm

/-- For every divisor-tail cursor, all quotient cells strictly before its
current pointer have product degree at least the current outer-loop bound.
This is the semantic statement that the cursor prefix has already been
processed; unlike a fuel counter, the bound is the actual frontier degree
used by the well-founded source loop. -/
def PairVecDivVHCCursorPrefixAbove (degreeLimit : Nat)
    (nodes : Array PairVecDivVHCNode)
    (quotient divisor : SparsePolyZp) : Prop :=
  ∀ (i : Nat) (node : PairVecDivVHCNode), nodes[i]? = some node →
    ∀ (q : Nat) (quotientTerm divisorTerm : UMonomial × Zp),
      q < node.quotientIndex →
      quotient[q]? = some quotientTerm →
      divisor[node.divisorIndex]? = some divisorTerm →
      degreeLimit ≤ quotientTerm.1.deg + divisorTerm.1.deg

/-- Every concrete cursor is either within the current quotient or exactly at
its one-past-the-end reset position. -/
def PairVecDivVHCCursorIndicesBounded (quotientSize : Nat)
    (nodes : Array PairVecDivVHCNode) : Prop :=
  ∀ (i : Nat) (node : PairVecDivVHCNode), nodes[i]? = some node →
    node.quotientIndex ≤ quotientSize

theorem PairVecDivVHCCursorPrefixAbove.mono
    (largeLimit smallLimit : Nat) (nodes : Array PairVecDivVHCNode)
    (quotient divisor : SparsePolyZp)
    (hprefix : PairVecDivVHCCursorPrefixAbove largeLimit nodes quotient divisor)
    (hlimits : smallLimit ≤ largeLimit) :
    PairVecDivVHCCursorPrefixAbove smallLimit nodes quotient divisor := by
  intro i node hget q quotientTerm divisorTerm hq hquotient hdivisor
  exact Nat.le_trans hlimits
    (hprefix i node hget q quotientTerm divisorTerm hq hquotient hdivisor)

theorem PairVecDivVHCCursorPrefixAbove.earlier_product_degree_gt
    (degreeLimit frontierDegree : Nat)
    (nodes : Array PairVecDivVHCNode) (quotient divisor : SparsePolyZp)
    (hprefix : PairVecDivVHCCursorPrefixAbove degreeLimit nodes quotient divisor)
    (hfrontier : frontierDegree < degreeLimit)
    (i : Nat) (node : PairVecDivVHCNode) (hget : nodes[i]? = some node)
    (q : Nat) (quotientTerm divisorTerm : UMonomial × Zp)
    (hearlier : q < node.quotientIndex)
    (hquotient : quotient[q]? = some quotientTerm)
    (hdivisor : divisor[node.divisorIndex]? = some divisorTerm) :
    frontierDegree < quotientTerm.1.deg + divisorTerm.1.deg := by
  have habove := hprefix i node hget q quotientTerm divisorTerm hearlier
    hquotient hdivisor
  omega

theorem PairVecDivVHCCursorPrefixAbove.push
    (degreeLimit : Nat) (nodes : Array PairVecDivVHCNode)
    (quotient divisor : SparsePolyZp) (term : UMonomial × Zp)
    (hprefix : PairVecDivVHCCursorPrefixAbove degreeLimit nodes quotient divisor)
    (hbounded : PairVecDivVHCCursorIndicesBounded quotient.size nodes) :
    PairVecDivVHCCursorPrefixAbove degreeLimit nodes (quotient.push term)
      divisor := by
  intro i node hget q quotientTerm divisorTerm hq hquotient hdivisor
  rw [Array.getElem?_push] at hquotient
  have hne : q ≠ quotient.size := by
    have hcursor := hbounded i node hget
    omega
  simp only [hne, ↓reduceIte] at hquotient
  exact hprefix i node hget q quotientTerm divisorTerm hq hquotient hdivisor

theorem PairVecDivVHCCursorPrefixAbove.set_fields
    (degreeLimit nodeIndex : Nat) (nodes : Array PairVecDivVHCNode)
    (quotient divisor : SparsePolyZp) (node updated : PairVecDivVHCNode)
    (hn : nodeIndex < nodes.size) (hget : nodes[nodeIndex]? = some node)
    (hprefix : PairVecDivVHCCursorPrefixAbove degreeLimit nodes quotient divisor)
    (hquotientIndex : updated.quotientIndex = node.quotientIndex)
    (hdivisorIndex : updated.divisorIndex = node.divisorIndex) :
    PairVecDivVHCCursorPrefixAbove degreeLimit
      (nodes.set nodeIndex updated) quotient divisor := by
  intro i current hcurrentGet q quotientTerm divisorTerm hq hquotient hdivisor
  by_cases heq : nodeIndex = i
  · subst i
    rw [Array.getElem?_set_self hn] at hcurrentGet
    simp only [Option.some.injEq] at hcurrentGet
    subst current
    rw [hquotientIndex] at hq
    rw [hdivisorIndex] at hdivisor
    exact hprefix nodeIndex node hget q quotientTerm divisorTerm hq hquotient
      hdivisor
  · rw [Array.getElem?_set_ne hn heq] at hcurrentGet
    exact hprefix i current hcurrentGet q quotientTerm divisorTerm hq hquotient
      hdivisor

theorem pairVecDivVHCInit_cursorPrefixAbove
    (degreeLimit : Nat) (quotient divisor : SparsePolyZp) :
    PairVecDivVHCCursorPrefixAbove degreeLimit (pairVecDivVHCInit divisor)
      quotient divisor := by
  intro i node hget q quotientTerm divisorTerm hq hquotient hdivisor
  have hi : i < divisor.size - 1 := by
    by_contra hnot
    rw [Array.getElem?_eq_none (by
      rw [pairVecDivVHCInit_size]
      omega)] at hget
    contradiction
  rw [pairVecDivVHCInit_get divisor i hi] at hget
  simp only [Option.some.injEq] at hget
  subst node
  simp [pairVecDivVHCInitialNode] at hq

theorem pairVecDivVHCInit_cursorIndicesBounded (divisor : SparsePolyZp) :
    PairVecDivVHCCursorIndicesBounded 0 (pairVecDivVHCInit divisor) := by
  intro i node hget
  have hi : i < divisor.size - 1 := by
    by_contra hnot
    rw [Array.getElem?_eq_none (by
      rw [pairVecDivVHCInit_size]
      omega)] at hget
    contradiction
  rw [pairVecDivVHCInit_get divisor i hi] at hget
  simp only [Option.some.injEq] at hget
  subst node
  simp [pairVecDivVHCInitialNode]

theorem PairVecDivVHCCursorPrefixAbove.set_advance
    (degreeLimit nodeIndex : Nat) (nodes : Array PairVecDivVHCNode)
    (quotient divisor : SparsePolyZp) (node updated : PairVecDivVHCNode)
    (hn : nodeIndex < nodes.size) (hget : nodes[nodeIndex]? = some node)
    (hprefix : PairVecDivVHCCursorPrefixAbove degreeLimit nodes quotient divisor)
    (hquotientIndex : updated.quotientIndex = node.quotientIndex + 1)
    (hdivisorIndex : updated.divisorIndex = node.divisorIndex)
    (hcurrent : ∀ (quotientTerm divisorTerm : UMonomial × Zp),
      quotient[node.quotientIndex]? = some quotientTerm →
      divisor[node.divisorIndex]? = some divisorTerm →
      degreeLimit ≤ quotientTerm.1.deg + divisorTerm.1.deg) :
    PairVecDivVHCCursorPrefixAbove degreeLimit
      (nodes.set nodeIndex updated) quotient divisor := by
  intro i current hcurrentGet q quotientTerm divisorTerm hq hquotient hdivisor
  by_cases heq : nodeIndex = i
  · subst i
    rw [Array.getElem?_set_self hn] at hcurrentGet
    simp only [Option.some.injEq] at hcurrentGet
    subst current
    rw [hquotientIndex] at hq
    by_cases hbefore : q < node.quotientIndex
    · exact hprefix nodeIndex node hget q quotientTerm divisorTerm hbefore
        hquotient (by simpa [hdivisorIndex] using hdivisor)
    · have hqeq : q = node.quotientIndex := by omega
      subst q
      exact hcurrent quotientTerm divisorTerm hquotient
        (by simpa [hdivisorIndex] using hdivisor)
  · rw [Array.getElem?_set_ne hn heq] at hcurrentGet
    exact hprefix i current hcurrentGet q quotientTerm divisorTerm hq
      hquotient hdivisor

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

/-- Degree-level form of heap order on the observable prefix.  It names the
actual heap heads and successful mono reads, which makes concrete `Array.set`
effects during sift-down tractable. -/
def PairVecDivVHCHeapDegreesOrderedUpTo (limit : Nat) (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode) : Prop :=
  ∀ (child : Nat), child < limit → 0 < child →
    ∀ childHead parentHead childMono parentMono,
      heap[child]? = some childHead →
      heap[pairVecDivVHCParent child]? = some parentHead →
      pairVecDivVHCMono childHead nodes = .ok childMono →
      pairVecDivVHCMono parentHead nodes = .ok parentMono →
      childMono.deg ≤ parentMono.deg

theorem PairVecDivVHCHeapOrdered.degreesUpTo
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode)
    (limit : Nat) (hlimit : limit ≤ heap.size)
    (hordered : PairVecDivVHCHeapOrdered heap nodes) :
    PairVecDivVHCHeapDegreesOrderedUpTo limit heap nodes := by
  intro child hchild hpos childHead parentHead childMono parentMono
    hchildGet hparentGet hchildMono hparentMono
  have hchildSlot : child < heap.size := Nat.lt_of_lt_of_le hchild hlimit
  have hparentSlot : pairVecDivVHCParent child < heap.size :=
    Nat.lt_trans (by
      unfold pairVecDivVHCParent
      have hhalf : (child - 1) / 2 ≤ child - 1 := Nat.div_le_self _ _
      omega) hchildSlot
  have hchildEq : heap[child] = childHead := by
    rw [Array.getElem?_eq_getElem hchildSlot] at hchildGet
    exact Option.some.inj hchildGet
  have hparentEq : heap[pairVecDivVHCParent child] = parentHead := by
    rw [Array.getElem?_eq_getElem hparentSlot] at hparentGet
    exact Option.some.inj hparentGet
  rw [← hchildEq] at hchildMono
  rw [← hparentEq] at hparentMono
  rcases (pairVecDivVHCMono_eq_ok_iff childHead nodes childMono).mp
      (by simpa [hchildEq] using hchildMono) with
    ⟨childNode, hchildNode, hchildActive⟩
  rcases (pairVecDivVHCMono_eq_ok_iff parentHead nodes parentMono).mp
      (by simpa [hparentEq] using hparentMono) with
    ⟨parentNode, hparentNode, hparentActive⟩
  have hchildMap :
      (nodes[childHead]?.map PairVecDivVHCNode.mono).join =
        some childMono := by
    rw [hchildNode]
    simp [hchildActive]
  have hparentMap :
      (nodes[parentHead]?.map PairVecDivVHCNode.mono).join =
        some parentMono := by
    rw [hparentNode]
    simp [hparentActive]
  exact hordered child (pairVecDivVHCParent child) hchildSlot rfl hpos
    childHead parentHead childMono parentMono hchildGet hparentGet hchildMap
    hparentMap

theorem PairVecDivVHCHeapDegreesOrderedUpTo.toHeapOrdered
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode)
    (hordered : PairVecDivVHCHeapDegreesOrderedUpTo heap.size heap nodes) :
    PairVecDivVHCHeapOrdered heap nodes := by
  intro child parent hchild hparent hpos childHead parentHead childMono
    parentMono hchildGet hparentGet hchildMap hparentMap
  subst parent
  have hchildMono : pairVecDivVHCMono childHead nodes = .ok childMono := by
    apply (pairVecDivVHCMono_eq_ok_iff childHead nodes childMono).mpr
    cases hnode : nodes[childHead]? with
    | none => simp [hnode] at hchildMap
    | some node =>
        cases hactive : node.mono with
        | none => simp [hnode, hactive] at hchildMap
        | some mono =>
            simp only [hnode, hactive, Option.map_some, Option.join_some,
              Option.some.injEq] at hchildMap
            subst mono
            exact ⟨node, rfl, hactive⟩
  have hparentMono : pairVecDivVHCMono parentHead nodes = .ok parentMono := by
    apply (pairVecDivVHCMono_eq_ok_iff parentHead nodes parentMono).mpr
    cases hnode : nodes[parentHead]? with
    | none => simp [hnode] at hparentMap
    | some node =>
        cases hactive : node.mono with
        | none => simp [hnode, hactive] at hparentMap
        | some mono =>
            simp only [hnode, hactive, Option.map_some, Option.join_some,
              Option.some.injEq] at hparentMap
            subst mono
            exact ⟨node, rfl, hactive⟩
  exact hordered child hchild hpos childHead parentHead childMono parentMono
    hchildGet hparentGet hchildMono hparentMono

theorem PairVecDivVHCHeapDegreesOrderedUpTo.pop
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode)
    (hordered : PairVecDivVHCHeapDegreesOrderedUpTo (heap.size - 1)
      heap nodes) :
    PairVecDivVHCHeapDegreesOrderedUpTo heap.pop.size heap.pop nodes := by
  intro child hchild hpos childHead parentHead childMono parentMono
    hchildGet hparentGet hchildMono hparentMono
  have hchildBound : child < heap.size - 1 := by
    simpa [Array.size_pop] using hchild
  have hparentLt : pairVecDivVHCParent child < child := by
    unfold pairVecDivVHCParent
    have hhalf : (child - 1) / 2 ≤ child - 1 := Nat.div_le_self _ _
    omega
  have hparentBound : pairVecDivVHCParent child < heap.size - 1 :=
    Nat.lt_trans hparentLt hchildBound
  rw [Array.getElem?_pop] at hchildGet hparentGet
  simp only [hchildBound, hparentBound, ↓reduceDIte] at hchildGet hparentGet
  exact hordered child hchildBound hpos childHead parentHead childMono
    parentMono hchildGet hparentGet hchildMono hparentMono

theorem PairVecDivVHCHeapDegreesOrderedUpTo.set_of_affected
    (limit i newHead : Nat) (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode) (hi : i < heap.size)
    (hordered : PairVecDivVHCHeapDegreesOrderedUpTo limit heap nodes)
    (haffected : ∀ (child : Nat), child < limit → 0 < child →
      (child = i ∨ pairVecDivVHCParent child = i) →
      ∀ childHead parentHead childMono parentMono,
        (heap.set i newHead)[child]? = some childHead →
        (heap.set i newHead)[pairVecDivVHCParent child]? = some parentHead →
        pairVecDivVHCMono childHead nodes = .ok childMono →
        pairVecDivVHCMono parentHead nodes = .ok parentMono →
        childMono.deg ≤ parentMono.deg) :
    PairVecDivVHCHeapDegreesOrderedUpTo limit (heap.set i newHead) nodes := by
  intro child hchild hpos childHead parentHead childMono parentMono
    hchildGet hparentGet hchildMono hparentMono
  by_cases hchanged : child = i ∨ pairVecDivVHCParent child = i
  · exact haffected child hchild hpos hchanged childHead parentHead childMono
      parentMono hchildGet hparentGet hchildMono hparentMono
  · have hchildNe : i ≠ child := by
      intro heq
      exact hchanged (Or.inl heq.symm)
    have hparentNe : i ≠ pairVecDivVHCParent child := by
      intro heq
      exact hchanged (Or.inr heq.symm)
    rw [Array.getElem?_set_ne hi hchildNe] at hchildGet
    rw [Array.getElem?_set_ne hi hparentNe] at hparentGet
    exact hordered child hchild hpos childHead parentHead childMono parentMono
      hchildGet hparentGet hchildMono hparentMono

theorem PairVecDivVHCHeapDegreesOrderedUpTo.set_leaf
    (limit i newHead : Nat) (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode) (newMono : UMonomial)
    (hi : i < heap.size)
    (hordered : PairVecDivVHCHeapDegreesOrderedUpTo limit heap nodes)
    (hnewMono : pairVecDivVHCMono newHead nodes = .ok newMono)
    (hnoChildren : ∀ child < limit, 0 < child →
      pairVecDivVHCParent child ≠ i)
    (hup : ∀ parentHead parentMono,
      heap[pairVecDivVHCParent i]? = some parentHead →
      pairVecDivVHCMono parentHead nodes = .ok parentMono →
      newMono.deg ≤ parentMono.deg) :
    PairVecDivVHCHeapDegreesOrderedUpTo limit (heap.set i newHead) nodes := by
  apply hordered.set_of_affected limit i newHead heap nodes hi
  intro child hchild hpos hchanged childHead parentHead childMono parentMono
    hchildGet hparentGet hchildMono hparentMono
  rcases hchanged with rfl | hparentEq
  · rw [Array.getElem?_set_self hi] at hchildGet
    have hchildHeadEq : childHead = newHead := (Option.some.inj hchildGet).symm
    subst childHead
    rw [hnewMono] at hchildMono
    have hchildMonoEq : childMono = newMono := (Except.ok.inj hchildMono).symm
    subst childMono
    have hparentLt : pairVecDivVHCParent child < child := by
      unfold pairVecDivVHCParent
      have hhalf : (child - 1) / 2 ≤ child - 1 := Nat.div_le_self _ _
      omega
    rw [Array.getElem?_set_ne hi (by omega)] at hparentGet
    exact hup parentHead parentMono hparentGet hparentMono
  · exact False.elim (hnoChildren child hchild hpos hparentEq)

theorem PairVecDivVHCHeapDegreesOrderedUpTo.set_parent
    (limit i newHead : Nat) (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode) (newMono : UMonomial)
    (hi : i < heap.size)
    (hordered : PairVecDivVHCHeapDegreesOrderedUpTo limit heap nodes)
    (hnewMono : pairVecDivVHCMono newHead nodes = .ok newMono)
    (hup : ∀ parentHead parentMono,
      heap[pairVecDivVHCParent i]? = some parentHead →
      pairVecDivVHCMono parentHead nodes = .ok parentMono →
      newMono.deg ≤ parentMono.deg)
    (hdown : ∀ child < limit, 0 < child → pairVecDivVHCParent child = i →
      ∀ childHead childMono,
        heap[child]? = some childHead →
        pairVecDivVHCMono childHead nodes = .ok childMono →
        childMono.deg ≤ newMono.deg) :
    PairVecDivVHCHeapDegreesOrderedUpTo limit (heap.set i newHead) nodes := by
  apply hordered.set_of_affected limit i newHead heap nodes hi
  intro child hchild hpos hchanged childHead parentHead childMono parentMono
    hchildGet hparentGet hchildMono hparentMono
  rcases hchanged with rfl | hparentEq
  · rw [Array.getElem?_set_self hi] at hchildGet
    have hchildHeadEq : childHead = newHead := (Option.some.inj hchildGet).symm
    subst childHead
    rw [hnewMono] at hchildMono
    have hchildMonoEq : childMono = newMono := (Except.ok.inj hchildMono).symm
    subst childMono
    have hparentLt : pairVecDivVHCParent child < child := by
      unfold pairVecDivVHCParent
      have hhalf : (child - 1) / 2 ≤ child - 1 := Nat.div_le_self _ _
      omega
    rw [Array.getElem?_set_ne hi (by omega)] at hparentGet
    exact hup parentHead parentMono hparentGet hparentMono
  · rw [hparentEq, Array.getElem?_set_self hi] at hparentGet
    have hparentHeadEq : parentHead = newHead :=
      (Option.some.inj hparentGet).symm
    subst parentHead
    rw [hnewMono] at hparentMono
    have hparentMonoEq : parentMono = newMono :=
      (Except.ok.inj hparentMono).symm
    subst parentMono
    have hchildNe : i ≠ child := by
      have hparentLt : pairVecDivVHCParent child < child := by
        unfold pairVecDivVHCParent
        have hhalf : (child - 1) / 2 ≤ child - 1 := Nat.div_le_self _ _
        omega
      omega
    rw [Array.getElem?_set_ne hi hchildNe] at hchildGet
    exact hdown child hchild hpos hparentEq childHead childMono hchildGet hchildMono

theorem pairVecDivVHCParent_eq_iff_children (i child : Nat)
    (hchild : 0 < child) :
    pairVecDivVHCParent child = i ↔
      child = i * 2 + 1 ∨ child = i * 2 + 2 := by
  constructor
  · intro hparent
    have hdiv : (child - 1) / 2 = i := by
      simpa [pairVecDivVHCParent] using hparent
    have hlower : i * 2 ≤ child - 1 := by
      apply (Nat.le_div_iff_mul_le (by omega : 0 < 2)).mp
      rw [hdiv]
    have hupper : child - 1 < (i + 1) * 2 := by
      apply (Nat.div_lt_iff_lt_mul (by omega : 0 < 2)).mp
      rw [hdiv]
      omega
    omega
  · intro hchildren
    rcases hchildren with rfl | rfl <;>
      unfold pairVecDivVHCParent <;> omega

theorem PairVecDivVHCHeapDegreesOrderedUpTo.copy_selected_child
    (limit i left selected holeHead leftHead rightHead selectedHead : Nat)
    (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode)
    (holeMono leftMono rightMono selectedMono : UMonomial)
    (hi : i < heap.size) (hileft : left = i * 2 + 1)
    (hrightBound : left + 1 ≤ limit) (hselectedLimit : selected < limit)
    (hordered : PairVecDivVHCHeapDegreesOrderedUpTo limit heap nodes)
    (hholeGet : heap[i]? = some holeHead)
    (hholeMono : pairVecDivVHCMono holeHead nodes = .ok holeMono)
    (hleftGet : heap[left]? = some leftHead)
    (hleftMono : pairVecDivVHCMono leftHead nodes = .ok leftMono)
    (hrightGet : heap[left + 1]? = some rightHead)
    (hrightMono : pairVecDivVHCMono rightHead nodes = .ok rightMono)
    (hselected : selected = left ∨ selected = left + 1)
    (hselectedGet : heap[selected]? = some selectedHead)
    (hselectedMono : pairVecDivVHCMono selectedHead nodes = .ok selectedMono)
    (hleftLe : leftMono.deg ≤ selectedMono.deg)
    (hrightLe : rightMono.deg ≤ selectedMono.deg) :
    PairVecDivVHCHeapDegreesOrderedUpTo limit
      (heap.set i selectedHead) nodes := by
  have hiLimit : i < limit := by omega
  have hselectedPos : 0 < selected := by rcases hselected with rfl | rfl <;> omega
  have hselectedParent : pairVecDivVHCParent selected = i := by
    apply (pairVecDivVHCParent_eq_iff_children i selected hselectedPos).mpr
    rcases hselected with rfl | rfl
    · exact Or.inl hileft
    · exact Or.inr (by omega)
  have hselectedLeHole : selectedMono.deg ≤ holeMono.deg :=
    hordered selected hselectedLimit hselectedPos selectedHead holeHead
      selectedMono holeMono hselectedGet (by simpa [hselectedParent] using hholeGet)
      hselectedMono hholeMono
  apply hordered.set_parent limit i selectedHead heap nodes selectedMono hi
    hselectedMono
  · intro parentHead parentMono hparentGet hparentMono
    by_cases hizero : i = 0
    · have hparentZero : pairVecDivVHCParent i = i := by
        simp [hizero, pairVecDivVHCParent]
      rw [hparentZero] at hparentGet
      rw [hholeGet] at hparentGet
      have hheadEq : parentHead = holeHead := (Option.some.inj hparentGet).symm
      subst parentHead
      rw [hholeMono] at hparentMono
      have hmonoEq : parentMono = holeMono :=
        (Except.ok.inj hparentMono).symm
      subst parentMono
      exact hselectedLeHole
    · have hipos : 0 < i := Nat.pos_of_ne_zero hizero
      exact Nat.le_trans hselectedLeHole
        (hordered i hiLimit hipos holeHead parentHead holeMono parentMono
          hholeGet hparentGet hholeMono hparentMono)
  · intro child hchild hpos hparent childHead childMono hchildGet hchildMono
    rcases (pairVecDivVHCParent_eq_iff_children i child hpos).mp hparent with
      hchildEq | hchildEq
    · have heq : child = left := by omega
      rw [← heq] at hleftGet
      rw [hleftGet] at hchildGet
      have hheadEq : childHead = leftHead :=
        (Option.some.inj hchildGet).symm
      subst childHead
      rw [hleftMono] at hchildMono
      have hmonoEq : childMono = leftMono :=
        (Except.ok.inj hchildMono).symm
      subst childMono
      exact hleftLe
    · have heq : child = left + 1 := by omega
      rw [← heq] at hrightGet
      rw [hrightGet] at hchildGet
      have hheadEq : childHead = rightHead :=
        (Option.some.inj hchildGet).symm
      subst childHead
      rw [hrightMono] at hchildMono
      have hmonoEq : childMono = rightMono :=
        (Except.ok.inj hchildMono).symm
      subst childMono
      exact hrightLe

/-- The generated well-founded sift loop preserves degree heap order on the
active prefix.  `lastNode` is the source sentinel at `limit`; `hup` records
the only order fact needed when that sentinel is finally written into the
current hole. -/
theorem pairVecDivVHCSiftDown_degreesOrderedUpTo
    (i child limit lastNode : Nat) (heap shifted : Array Nat)
    (nodes : Array PairVecDivVHCNode) (holeMono lastMono : UMonomial)
    (hchildEq : child = i * 2 + 1)
    (hi : i < heap.size) (hlimit : limit < heap.size)
    (hlastSlot : heap[limit]? = some lastNode)
    (hholeMono : pairVecDivVHCMono heap[i] nodes = .ok holeMono)
    (hlastMono : pairVecDivVHCMono lastNode nodes = .ok lastMono)
    (hordered : PairVecDivVHCHeapDegreesOrderedUpTo limit heap nodes)
    (hup : ∀ parentHead parentMono,
      heap[pairVecDivVHCParent i]? = some parentHead →
      pairVecDivVHCMono parentHead nodes = .ok parentMono →
      lastMono.deg ≤ parentMono.deg)
    (hrun : pairVecDivVHCSiftDown i child limit lastNode heap nodes =
      .ok shifted) :
    PairVecDivVHCHeapDegreesOrderedUpTo limit shifted nodes := by
  rw [pairVecDivVHCSiftDown] at hrun
  simp only [hi, hlimit, hlastMono, true_and] at hrun
  by_cases hchild : child < limit
  · simp only [hchild, ↓reduceDIte] at hrun
    have hrightSize : child + 1 < heap.size := by omega
    cases hleftMono : pairVecDivVHCMono heap[child] nodes with
    | error fault => simp [hleftMono] at hrun
    | ok leftMono =>
      cases hrightMono : pairVecDivVHCMono heap[child + 1] nodes with
      | error fault => simp [hleftMono, hrightMono] at hrun
      | ok rightMono =>
        simp only [hleftMono, hrightMono] at hrun
        by_cases hselected : leftMono.deg > rightMono.deg
        · simp only [hselected, ↓reduceIte] at hrun
          by_cases hgreater : leftMono.deg > lastMono.deg
          · simp only [hgreater, ↓reduceDIte] at hrun
            have hcopy := hordered.copy_selected_child limit i child child
              heap[i] heap[child] heap[child + 1] heap[child] heap nodes
              holeMono leftMono rightMono leftMono hi hchildEq (by omega) hchild
              (Array.getElem?_eq_getElem hi)
              hholeMono
              (Array.getElem?_eq_getElem (by omega)) hleftMono
              (Array.getElem?_eq_getElem hrightSize) hrightMono
              (Or.inl rfl) (Array.getElem?_eq_getElem (by omega)) hleftMono
              (Nat.le_refl _) (Nat.le_of_lt hselected)
            exact pairVecDivVHCSiftDown_degreesOrderedUpTo child
              (child * 2 + 1) limit lastNode (heap.set i heap[child]) shifted
              nodes leftMono lastMono rfl
              (by simpa using (show child < heap.size by omega))
              (by simpa using hlimit) (by
                rw [Array.getElem?_set_ne hi (by omega)]
                exact hlastSlot)
              (by simpa [Array.getElem_set, hi] using hleftMono)
              hlastMono hcopy (by
                  intro parentHead parentMono hparentGet hparentMono
                  have hp : pairVecDivVHCParent child = i :=
                    (pairVecDivVHCParent_eq_iff_children i child (by omega)).2
                      (Or.inl hchildEq)
                  rw [hp, Array.getElem?_set_self hi] at hparentGet
                  have hhead : parentHead = heap[child] :=
                    (Option.some.inj hparentGet).symm
                  subst parentHead
                  rw [hleftMono] at hparentMono
                  have hm : parentMono = leftMono :=
                    (Except.ok.inj hparentMono).symm
                  subst parentMono
                  exact Nat.le_of_lt hgreater) hrun
          · simp only [hgreater, ↓reduceDIte, Except.ok.injEq] at hrun
            subst shifted
            apply hordered.set_parent limit i lastNode heap nodes lastMono hi
              hlastMono hup
            intro active hactive hpos hparent activeHead activeMono
              hactiveGet hactiveMono
            rcases (pairVecDivVHCParent_eq_iff_children i active hpos).1
              hparent with hleft | hright
            · have : active = child := by omega
              rw [this] at hactiveGet
              rw [Array.getElem?_eq_getElem (by omega)] at hactiveGet
              have hh : activeHead = heap[child] :=
                (Option.some.inj hactiveGet).symm
              subst activeHead
              rw [hleftMono] at hactiveMono
              have hm : activeMono = leftMono :=
                (Except.ok.inj hactiveMono).symm
              subst activeMono
              exact Nat.le_trans (Nat.le_refl _)
                (Nat.le_of_not_gt hgreater)
            · have : active = child + 1 := by omega
              rw [this] at hactiveGet
              rw [Array.getElem?_eq_getElem hrightSize] at hactiveGet
              have hh : activeHead = heap[child + 1] :=
                (Option.some.inj hactiveGet).symm
              subst activeHead
              rw [hrightMono] at hactiveMono
              have hm : activeMono = rightMono :=
                (Except.ok.inj hactiveMono).symm
              subst activeMono
              exact Nat.le_trans (Nat.le_of_lt hselected)
                (Nat.le_of_not_gt hgreater)
        · simp only [hselected, ↓reduceIte] at hrun
          by_cases hgreater : rightMono.deg > lastMono.deg
          · simp only [hgreater, ↓reduceDIte] at hrun
            have hselectedLimit : child + 1 < limit := by
              by_contra hnot
              have heq : child + 1 = limit := by omega
              have hrightHead : heap[child + 1]? = some lastNode := by
                simpa [heq] using hlastSlot
              have harray : heap[child + 1]? = some heap[child + 1] :=
                Array.getElem?_eq_getElem hrightSize
              rw [hrightHead] at harray
              have hheadEq : heap[child + 1] = lastNode :=
                (Option.some.inj harray).symm
              rw [hheadEq, hlastMono] at hrightMono
              have hm : rightMono = lastMono :=
                (Except.ok.inj hrightMono).symm
              subst rightMono
              omega
            have hcopy := hordered.copy_selected_child limit i child
              (child + 1) heap[i] heap[child] heap[child + 1]
              heap[child + 1] heap nodes holeMono leftMono rightMono rightMono
              hi hchildEq (by omega) hselectedLimit
              (Array.getElem?_eq_getElem hi) hholeMono
              (Array.getElem?_eq_getElem (by omega)) hleftMono
              (Array.getElem?_eq_getElem hrightSize) hrightMono
              (Or.inr rfl) (Array.getElem?_eq_getElem hrightSize) hrightMono
              (Nat.le_of_not_gt hselected) (Nat.le_refl _)
            exact pairVecDivVHCSiftDown_degreesOrderedUpTo (child + 1)
              ((child + 1) * 2 + 1) limit lastNode
              (heap.set i heap[child + 1]) shifted nodes rightMono lastMono rfl
              (by simpa using (show child + 1 < heap.size by omega))
              (by simpa using hlimit) (by
                rw [Array.getElem?_set_ne hi (by omega)]
                exact hlastSlot)
              (by simpa [Array.getElem_set, hi] using hrightMono)
              hlastMono hcopy (by
                  intro parentHead parentMono hparentGet hparentMono
                  have hp : pairVecDivVHCParent (child + 1) = i :=
                    (pairVecDivVHCParent_eq_iff_children i (child + 1)
                      (by omega)).2 (Or.inr (by omega))
                  rw [hp, Array.getElem?_set_self hi] at hparentGet
                  have hhead : parentHead = heap[child + 1] :=
                    (Option.some.inj hparentGet).symm
                  subst parentHead
                  rw [hrightMono] at hparentMono
                  have hm : parentMono = rightMono :=
                    (Except.ok.inj hparentMono).symm
                  subst parentMono
                  exact Nat.le_of_lt hgreater) hrun
          · simp only [hgreater, ↓reduceDIte, Except.ok.injEq] at hrun
            subst shifted
            apply hordered.set_parent limit i lastNode heap nodes lastMono hi
              hlastMono hup
            intro active hactive hpos hparent activeHead activeMono
              hactiveGet hactiveMono
            rcases (pairVecDivVHCParent_eq_iff_children i active hpos).1
              hparent with hleft | hright
            · have heq : active = child := by omega
              rw [heq] at hactiveGet
              rw [Array.getElem?_eq_getElem (by omega)] at hactiveGet
              have hh : activeHead = heap[child] :=
                (Option.some.inj hactiveGet).symm
              subst activeHead
              rw [hleftMono] at hactiveMono
              have hm : activeMono = leftMono :=
                (Except.ok.inj hactiveMono).symm
              subst activeMono
              exact Nat.le_trans (Nat.le_of_not_gt hselected)
                (Nat.le_of_not_gt hgreater)
            · have heq : active = child + 1 := by omega
              rw [heq] at hactiveGet
              rw [Array.getElem?_eq_getElem hrightSize] at hactiveGet
              have hh : activeHead = heap[child + 1] :=
                (Option.some.inj hactiveGet).symm
              subst activeHead
              rw [hrightMono] at hactiveMono
              have hm : activeMono = rightMono :=
                (Except.ok.inj hactiveMono).symm
              subst activeMono
              exact Nat.le_of_not_gt hgreater
  · simp only [hchild, ↓reduceDIte, Except.ok.injEq] at hrun
    subst shifted
    apply hordered.set_leaf limit i lastNode heap nodes lastMono hi hlastMono
    · intro active hactive hpos hparent
      rcases (pairVecDivVHCParent_eq_iff_children i active hpos).1 hparent with
        hleft | hright <;> omega
    · exact hup
termination_by limit - child
decreasing_by
  all_goals omega

/-- Degree-prefix heap order plus pointer validity makes the root dominate
every concrete slot. -/
theorem PairVecDivVHCHeapDegreesOrderedUpTo.slot_le_root
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode)
    (hvalid : PairVecDivVHCHeapPointersValid heap nodes)
    (hordered : PairVecDivVHCHeapDegreesOrderedUpTo heap.size heap nodes)
    (slot : Nat) (hslot : slot < heap.size)
    (mono rootMono : UMonomial)
    (hmono : pairVecDivVHCMono heap[slot] nodes = .ok mono)
    (hroot : pairVecDivVHCMono heap[0] nodes = .ok rootMono) :
    mono.deg ≤ rootMono.deg := by
  induction slot using Nat.strong_induction_on generalizing mono with
  | h slot ih =>
      by_cases hzero : slot = 0
      · subst slot
        rw [hroot] at hmono
        have hm : mono = rootMono := (Except.ok.inj hmono).symm
        subst mono
        exact Nat.le_refl _
      · have hpos : 0 < slot := Nat.pos_of_ne_zero hzero
        let parent := pairVecDivVHCParent slot
        have hparentLt : parent < slot := by
          unfold parent pairVecDivVHCParent
          have hhalf : (slot - 1) / 2 ≤ slot - 1 := Nat.div_le_self _ _
          omega
        have hparentSlot : parent < heap.size :=
          Nat.lt_trans hparentLt hslot
        rcases hvalid parent hparentSlot with
          ⟨parentHead, parentNode, parentMono, hparentGet, hparentNode,
            hparentActive⟩
        have hparentHead : heap[parent] = parentHead := by
          rw [Array.getElem?_eq_getElem hparentSlot] at hparentGet
          exact Option.some.inj hparentGet
        have hparentMono : pairVecDivVHCMono heap[parent] nodes =
            .ok parentMono := by
          apply (pairVecDivVHCMono_eq_ok_iff heap[parent] nodes parentMono).2
          rw [hparentHead]
          exact ⟨parentNode, hparentNode, hparentActive⟩
        have hstep := hordered slot hslot hpos heap[slot] heap[parent] mono
          parentMono (Array.getElem?_eq_getElem hslot)
          (Array.getElem?_eq_getElem hparentSlot) hmono hparentMono
        exact Nat.le_trans hstep
          (ih parent hparentLt hparentSlot parentMono hparentMono hroot)

/-- Root entry for sift after the old root bucket has been consumed.  It only
assumes order on edges strictly below the root; the first concrete root write
re-establishes the two missing root edges from the source comparisons. -/
theorem pairVecDivVHCSiftDown_root_of_nonroot_order
    (limit lastNode : Nat) (heap shifted : Array Nat)
    (nodes : Array PairVecDivVHCNode) (lastMono : UMonomial)
    (hroot : 0 < heap.size) (hlimit : limit < heap.size)
    (hlastSlot : heap[limit]? = some lastNode)
    (hlastMono : pairVecDivVHCMono lastNode nodes = .ok lastMono)
    (haway : ∀ (child : Nat), child < limit → 0 < child →
      0 < pairVecDivVHCParent child →
      ∀ childHead parentHead childMono parentMono,
        heap[child]? = some childHead →
        heap[pairVecDivVHCParent child]? = some parentHead →
        pairVecDivVHCMono childHead nodes = .ok childMono →
        pairVecDivVHCMono parentHead nodes = .ok parentMono →
        childMono.deg ≤ parentMono.deg)
    (hrun : pairVecDivVHCSiftDown 0 1 limit lastNode heap nodes =
      .ok shifted) :
    PairVecDivVHCHeapDegreesOrderedUpTo limit shifted nodes := by
  have rootSetOrdered (newHead : Nat) (newMono : UMonomial)
      (hnewMono : pairVecDivVHCMono newHead nodes = .ok newMono)
      (hchildren : ∀ child < limit, 0 < child →
        pairVecDivVHCParent child = 0 →
        ∀ childHead childMono, heap[child]? = some childHead →
          pairVecDivVHCMono childHead nodes = .ok childMono →
          childMono.deg ≤ newMono.deg) :
      PairVecDivVHCHeapDegreesOrderedUpTo limit
        (heap.set 0 newHead) nodes := by
    intro child hchild hpos childHead parentHead childMono parentMono
      hchildGet hparentGet hchildMono hparentMono
    by_cases hpzero : pairVecDivVHCParent child = 0
    · rw [hpzero, Array.getElem?_set_self hroot] at hparentGet
      have hparentHead : parentHead = newHead :=
        (Option.some.inj hparentGet).symm
      subst parentHead
      rw [hnewMono] at hparentMono
      have hm : parentMono = newMono :=
        (Except.ok.inj hparentMono).symm
      subst parentMono
      rw [Array.getElem?_set_ne hroot (by omega)] at hchildGet
      exact hchildren child hchild hpos hpzero childHead childMono hchildGet
        hchildMono
    · have hparentPos : 0 < pairVecDivVHCParent child :=
        Nat.pos_of_ne_zero hpzero
      rw [Array.getElem?_set_ne hroot (by omega)] at hchildGet
      rw [Array.getElem?_set_ne hroot (fun heq => hpzero heq.symm)] at hparentGet
      exact haway child hchild hpos hparentPos childHead parentHead childMono
        parentMono hchildGet hparentGet hchildMono hparentMono
  rw [pairVecDivVHCSiftDown] at hrun
  simp only [hroot, hlimit, hlastMono] at hrun
  by_cases hchild : 1 < limit
  · simp only [hchild, ↓reduceDIte] at hrun
    have hleftSize : 1 < heap.size := by omega
    have hrightSize : 2 < heap.size := by omega
    cases hleftMono : pairVecDivVHCMono heap[1] nodes with
    | error fault => simp [hleftMono] at hrun
    | ok leftMono =>
      cases hrightMono : pairVecDivVHCMono heap[2] nodes with
      | error fault => simp [hleftMono, hrightMono] at hrun
      | ok rightMono =>
        simp only [hleftMono, hrightMono] at hrun
        by_cases hselected : leftMono.deg > rightMono.deg
        · simp only [hselected, ↓reduceIte] at hrun
          by_cases hgreater : leftMono.deg > lastMono.deg
          · simp only [hgreater, ↓reduceDIte] at hrun
            have hordered' := rootSetOrdered heap[1] leftMono hleftMono (by
              intro child hchild' hpos hp childHead childMono hget hmono
              rcases (pairVecDivVHCParent_eq_iff_children 0 child hpos).1 hp
                with hc | hc
              · have : child = 1 := by omega
                rw [this, Array.getElem?_eq_getElem hleftSize] at hget
                have hh : childHead = heap[1] :=
                  (Option.some.inj hget).symm
                subst childHead
                rw [hleftMono] at hmono
                have hm : childMono = leftMono :=
                  (Except.ok.inj hmono).symm
                subst childMono
                exact Nat.le_refl _
              · have : child = 2 := by omega
                rw [this, Array.getElem?_eq_getElem hrightSize] at hget
                have hh : childHead = heap[2] :=
                  (Option.some.inj hget).symm
                subst childHead
                rw [hrightMono] at hmono
                have hm : childMono = rightMono :=
                  (Except.ok.inj hmono).symm
                subst childMono
                exact Nat.le_of_lt hselected)
            exact pairVecDivVHCSiftDown_degreesOrderedUpTo 1 3 limit
              lastNode (heap.set 0 heap[1]) shifted nodes leftMono lastMono
              rfl (by simpa) (by simpa using hlimit) (by
                rw [Array.getElem?_set_ne hroot (by omega)]
                exact hlastSlot)
              (by simpa [Array.getElem_set, hroot] using hleftMono)
              hlastMono hordered' (by
                intro parentHead parentMono hget hmono
                have hp : pairVecDivVHCParent 1 = 0 := by
                  simp [pairVecDivVHCParent]
                rw [hp, Array.getElem?_set_self hroot] at hget
                have hh : parentHead = heap[1] :=
                  (Option.some.inj hget).symm
                subst parentHead
                rw [hleftMono] at hmono
                have hm : parentMono = leftMono :=
                  (Except.ok.inj hmono).symm
                subst parentMono
                exact Nat.le_of_lt hgreater) hrun
          · simp only [hgreater, ↓reduceDIte, Except.ok.injEq] at hrun
            subst shifted
            exact rootSetOrdered lastNode lastMono hlastMono (by
              intro child hchild' hpos hp childHead childMono hget hmono
              rcases (pairVecDivVHCParent_eq_iff_children 0 child hpos).1 hp
                with hc | hc
              · have : child = 1 := by omega
                rw [this, Array.getElem?_eq_getElem hleftSize] at hget
                have hh : childHead = heap[1] :=
                  (Option.some.inj hget).symm
                subst childHead
                rw [hleftMono] at hmono
                have hm : childMono = leftMono :=
                  (Except.ok.inj hmono).symm
                subst childMono
                exact Nat.le_of_not_gt hgreater
              · have : child = 2 := by omega
                rw [this, Array.getElem?_eq_getElem hrightSize] at hget
                have hh : childHead = heap[2] :=
                  (Option.some.inj hget).symm
                subst childHead
                rw [hrightMono] at hmono
                have hm : childMono = rightMono :=
                  (Except.ok.inj hmono).symm
                subst childMono
                exact Nat.le_trans (Nat.le_of_lt hselected)
                  (Nat.le_of_not_gt hgreater))
        · simp only [hselected, ↓reduceIte] at hrun
          by_cases hgreater : rightMono.deg > lastMono.deg
          · simp only [hgreater, ↓reduceDIte] at hrun
            have hselectedLimit : 2 < limit := by
              by_contra hnot
              have heq : 2 = limit := by omega
              have hsaved : heap[2]? = some lastNode := by
                simpa [heq] using hlastSlot
              rw [Array.getElem?_eq_getElem hrightSize] at hsaved
              have hh : heap[2] = lastNode := Option.some.inj hsaved
              rw [hh, hlastMono] at hrightMono
              have hm : rightMono = lastMono :=
                (Except.ok.inj hrightMono).symm
              subst rightMono
              omega
            have hordered' := rootSetOrdered heap[2] rightMono hrightMono (by
              intro child hchild' hpos hp childHead childMono hget hmono
              rcases (pairVecDivVHCParent_eq_iff_children 0 child hpos).1 hp
                with hc | hc
              · have : child = 1 := by omega
                rw [this, Array.getElem?_eq_getElem hleftSize] at hget
                have hh : childHead = heap[1] :=
                  (Option.some.inj hget).symm
                subst childHead
                rw [hleftMono] at hmono
                have hm : childMono = leftMono :=
                  (Except.ok.inj hmono).symm
                subst childMono
                exact Nat.le_of_not_gt hselected
              · have : child = 2 := by omega
                rw [this, Array.getElem?_eq_getElem hrightSize] at hget
                have hh : childHead = heap[2] :=
                  (Option.some.inj hget).symm
                subst childHead
                rw [hrightMono] at hmono
                have hm : childMono = rightMono :=
                  (Except.ok.inj hmono).symm
                subst childMono
                exact Nat.le_refl _)
            exact pairVecDivVHCSiftDown_degreesOrderedUpTo 2 5 limit
              lastNode (heap.set 0 heap[2]) shifted nodes rightMono lastMono
              rfl (by simpa) (by simpa using hlimit) (by
                rw [Array.getElem?_set_ne hroot (by omega)]
                exact hlastSlot)
              (by simpa [Array.getElem_set, hroot] using hrightMono)
              hlastMono hordered' (by
                intro parentHead parentMono hget hmono
                have hp : pairVecDivVHCParent 2 = 0 := by
                  simp [pairVecDivVHCParent]
                rw [hp, Array.getElem?_set_self hroot] at hget
                have hh : parentHead = heap[2] :=
                  (Option.some.inj hget).symm
                subst parentHead
                rw [hrightMono] at hmono
                have hm : parentMono = rightMono :=
                  (Except.ok.inj hmono).symm
                subst parentMono
                exact Nat.le_of_lt hgreater) hrun
          · simp only [hgreater, ↓reduceDIte, Except.ok.injEq] at hrun
            subst shifted
            exact rootSetOrdered lastNode lastMono hlastMono (by
              intro child hchild' hpos hp childHead childMono hget hmono
              rcases (pairVecDivVHCParent_eq_iff_children 0 child hpos).1 hp
                with hc | hc
              · have : child = 1 := by omega
                rw [this, Array.getElem?_eq_getElem hleftSize] at hget
                have hh : childHead = heap[1] :=
                  (Option.some.inj hget).symm
                subst childHead
                rw [hleftMono] at hmono
                have hm : childMono = leftMono :=
                  (Except.ok.inj hmono).symm
                subst childMono
                exact Nat.le_trans (Nat.le_of_not_gt hselected)
                  (Nat.le_of_not_gt hgreater)
              · have : child = 2 := by omega
                rw [this, Array.getElem?_eq_getElem hrightSize] at hget
                have hh : childHead = heap[2] :=
                  (Option.some.inj hget).symm
                subst childHead
                rw [hrightMono] at hmono
                have hm : childMono = rightMono :=
                  (Except.ok.inj hmono).symm
                subst childMono
                exact Nat.le_of_not_gt hgreater)
  · simp only [hchild, ↓reduceDIte, Except.ok.injEq] at hrun
    subst shifted
    exact rootSetOrdered lastNode lastMono hlastMono (by
      intro child hchild' hpos hp
      omega)

/-- A successful generated extract preserves the complete source max-heap
order after the concrete `pop`. -/
theorem pairVecDivVHCExtract_preserves_heapOrdered
    (heap heap' : Array Nat) (nodes : Array PairVecDivVHCNode)
    (hvalid : PairVecDivVHCHeapPointersValid heap nodes)
    (hordered : PairVecDivVHCHeapOrdered heap nodes)
    (hrun : pairVecDivVHCExtract heap nodes = .ok heap') :
    PairVecDivVHCHeapOrdered heap' nodes := by
  have hsize := pairVecDivVHCExtract_size heap heap' nodes hrun
  by_cases hone : heap.size = 1
  · apply PairVecDivVHCHeapDegreesOrderedUpTo.toHeapOrdered
    intro child hchild
    have : heap'.size = 0 := by omega
    omega
  · have hnonempty : 0 < heap.size := by
      by_contra hempty
      have : heap.size = 0 := by omega
      unfold pairVecDivVHCExtract at hrun
      simp [this] at hrun
    have hlimit : heap.size - 1 < heap.size := by omega
    rcases hvalid 0 hnonempty with
      ⟨rootHead, rootNode, rootMono, hrootGet, hrootNode, hrootActive⟩
    have hrootHead : heap[0] = rootHead := by
      rw [Array.getElem?_eq_getElem hnonempty] at hrootGet
      exact Option.some.inj hrootGet
    have hrootMono : pairVecDivVHCMono heap[0] nodes = .ok rootMono := by
      apply (pairVecDivVHCMono_eq_ok_iff heap[0] nodes rootMono).2
      rw [hrootHead]
      exact ⟨rootNode, hrootNode, hrootActive⟩
    rcases hvalid (heap.size - 1) hlimit with
      ⟨lastHead, lastNodeData, lastMono, hlastGet, hlastNode,
        hlastActive⟩
    have hlastHead : heap[heap.size - 1] = lastHead := by
      rw [Array.getElem?_eq_getElem hlimit] at hlastGet
      exact Option.some.inj hlastGet
    have hlastMono : pairVecDivVHCMono heap[heap.size - 1] nodes =
        .ok lastMono := by
      apply (pairVecDivVHCMono_eq_ok_iff heap[heap.size - 1] nodes lastMono).2
      rw [hlastHead]
      exact ⟨lastNodeData, hlastNode, hlastActive⟩
    unfold pairVecDivVHCExtract at hrun
    simp only [hnonempty, ↓reduceDIte] at hrun
    cases hsift : pairVecDivVHCSiftDown 0 1 (heap.size - 1)
        heap[heap.size - 1] heap nodes with
    | error fault => simp [hsift] at hrun
    | ok shifted =>
        rw [hsift] at hrun
        simp only [Except.ok.injEq] at hrun
        subst heap'
        have hsiftOrdered := pairVecDivVHCSiftDown_degreesOrderedUpTo 0 1
          (heap.size - 1) heap[heap.size - 1] heap shifted nodes rootMono
          lastMono rfl hnonempty hlimit
          (Array.getElem?_eq_getElem hlimit) hrootMono hlastMono
          (hordered.degreesUpTo heap nodes (heap.size - 1) (by omega)) (by
            intro parentHead parentMono hparentGet hparentMono
            have hp : pairVecDivVHCParent 0 = 0 := by
              simp [pairVecDivVHCParent]
            rw [hp] at hparentGet
            rw [Array.getElem?_eq_getElem hnonempty] at hparentGet
            have hhead : parentHead = heap[0] :=
              (Option.some.inj hparentGet).symm
            subst parentHead
            rw [hrootMono] at hparentMono
            have hm : parentMono = rootMono :=
              (Except.ok.inj hparentMono).symm
            subst parentMono
            exact (hordered.degreesUpTo heap nodes heap.size (Nat.le_refl _)).slot_le_root
              heap nodes hvalid (heap.size - 1) hlimit lastMono rootMono
              hlastMono hrootMono) hsift
        have hshiftSize := pairVecDivVHCSiftDown_size 0 1 (heap.size - 1)
          heap[heap.size - 1] heap shifted nodes hsift
        have hsiftOrdered' : PairVecDivVHCHeapDegreesOrderedUpTo
            (shifted.size - 1) shifted nodes := by
          simpa [hshiftSize] using hsiftOrdered
        exact hsiftOrdered'.pop shifted nodes |>.toHeapOrdered

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

theorem pairVecDivVHCParent_lt (i : Nat) (hi : 0 < i) :
    pairVecDivVHCParent i < i := by
  unfold pairVecDivVHCParent
  have hhalf : (i - 1) / 2 ≤ i - 1 := Nat.div_le_self _ _
  omega

/-- Copying a heap parent into its child slot, as one upward-bubble step does,
preserves max-heap order.  The overwritten head already dominated its children,
and its parent dominates it, so the copied parent also dominates those children. -/
theorem pairVecDivVHCSet_child_to_parent_preserves_heapOrdered
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode)
    (slot parentHead : Nat)
    (hslot : slot < heap.size) (hpos : 0 < slot)
    (hparentGet : heap[pairVecDivVHCParent slot]? = some parentHead)
    (hvalid : PairVecDivVHCHeapPointersValid heap nodes)
    (hordered : PairVecDivVHCHeapOrdered heap nodes) :
    PairVecDivVHCHeapOrdered (heap.set slot parentHead) nodes := by
  rcases hvalid (pairVecDivVHCParent slot)
      (Nat.lt_trans (pairVecDivVHCParent_lt slot hpos) hslot) with
    ⟨actualParentHead, parentNode, parentMono, hactualParentGet,
      hparentNode, hparentActive⟩
  have hparentHeadEq : actualParentHead = parentHead := by
    rw [hparentGet] at hactualParentGet
    exact (Option.some.inj hactualParentGet).symm
  subst actualParentHead
  have hparentMono : pairVecDivVHCMono parentHead nodes = .ok parentMono := by
    exact (pairVecDivVHCMono_eq_ok_iff parentHead nodes parentMono).mpr
      ⟨parentNode, hparentNode, hparentActive⟩
  have hdegrees := hordered.degreesUpTo heap nodes heap.size (Nat.le_refl _)
  apply PairVecDivVHCHeapDegreesOrderedUpTo.toHeapOrdered
  have hsetOrdered : PairVecDivVHCHeapDegreesOrderedUpTo heap.size
      (heap.set slot parentHead) nodes := by
    apply hdegrees.set_parent heap.size slot parentHead heap nodes parentMono
      hslot hparentMono
    · intro grandparentHead grandparentMono hgrandparentGet hgrandparentMono
      rw [hparentGet] at hgrandparentGet
      have hheadEq : grandparentHead = parentHead :=
        (Option.some.inj hgrandparentGet).symm
      subst grandparentHead
      rw [hparentMono] at hgrandparentMono
      exact Nat.le_of_eq
        (congrArg UMonomial.deg (Except.ok.inj hgrandparentMono))
    · intro child hchild hchildPos hchildParent childHead childMono
        hchildGet hchildMono
      rcases hvalid slot hslot with
        ⟨oldHead, oldNode, oldMono, holdHeadGet, holdNode, holdActive⟩
      have holdMono : pairVecDivVHCMono oldHead nodes = .ok oldMono :=
        (pairVecDivVHCMono_eq_ok_iff oldHead nodes oldMono).mpr
          ⟨oldNode, holdNode, holdActive⟩
      have hchildLeOld : childMono.deg ≤ oldMono.deg :=
        hdegrees child hchild hchildPos childHead oldHead childMono oldMono
          hchildGet (by simpa [hchildParent] using holdHeadGet) hchildMono holdMono
      have holdLeParent : oldMono.deg ≤ parentMono.deg :=
        hdegrees slot hslot hpos oldHead parentHead oldMono parentMono
          holdHeadGet hparentGet holdMono hparentMono
      exact Nat.le_trans hchildLeOld holdLeParent
  simpa only [Array.size_set] using hsetOrdered

/-- Replacing a heap head by another active head of the same degree preserves
all parent-child comparisons.  This is the heap side of equal-degree bucketing. -/
theorem pairVecDivVHCSet_sameDegree_preserves_heapOrdered
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode)
    (slot newHead oldHead : Nat) (newMono oldMono : UMonomial)
    (hslot : slot < heap.size) (holdGet : heap[slot]? = some oldHead)
    (hnew : pairVecDivVHCMono newHead nodes = .ok newMono)
    (hold : pairVecDivVHCMono oldHead nodes = .ok oldMono)
    (hdegree : newMono.deg = oldMono.deg)
    (hordered : PairVecDivVHCHeapOrdered heap nodes) :
    PairVecDivVHCHeapOrdered (heap.set slot newHead) nodes := by
  have hdegrees := hordered.degreesUpTo heap nodes heap.size (Nat.le_refl _)
  apply PairVecDivVHCHeapDegreesOrderedUpTo.toHeapOrdered
  have hsetOrdered : PairVecDivVHCHeapDegreesOrderedUpTo heap.size
      (heap.set slot newHead) nodes := by
    apply hdegrees.set_parent heap.size slot newHead heap nodes newMono hslot hnew
    · intro parentHead parentMono hparentGet hparentMono
      by_cases hzero : slot = 0
      · subst slot
        simp only [pairVecDivVHCParent] at hparentGet
        rw [holdGet] at hparentGet
        have hparentHeadEq : parentHead = oldHead :=
          (Option.some.inj hparentGet).symm
        subst parentHead
        rw [hold] at hparentMono
        have hparentMonoEq : parentMono = oldMono :=
          (Except.ok.inj hparentMono).symm
        subst parentMono
        exact Nat.le_of_eq hdegree
      · have hpos : 0 < slot := Nat.pos_of_ne_zero hzero
        have holdLeParent := hdegrees slot hslot hpos oldHead parentHead
          oldMono parentMono holdGet hparentGet hold hparentMono
        exact hdegree.le.trans holdLeParent
    · intro child hchild hchildPos hchildParent childHead childMono
        hchildGet hchildMono
      have hchildLeOld := hdegrees child hchild hchildPos childHead oldHead
        childMono oldMono hchildGet (by simpa [hchildParent] using holdGet)
        hchildMono hold
      exact hchildLeOld.trans (Nat.le_of_eq hdegree.symm)
  simpa only [Array.size_set] using hsetOrdered

/-- A pointwise degree bound on every active head currently stored in a heap. -/
def PairVecDivVHCHeapBoundedBy (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode) (bound : UMonomial) : Prop :=
  ∀ (slot head : Nat) (mono : UMonomial), heap[slot]? = some head →
    pairVecDivVHCMono head nodes = .ok mono → mono.deg ≤ bound.deg

theorem PairVecDivVHCHeapPointersValid.set_from_slot
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode)
    (target source : Nat) (htarget : target < heap.size)
    (hsource : source < heap.size)
    (hvalid : PairVecDivVHCHeapPointersValid heap nodes) :
    PairVecDivVHCHeapPointersValid (heap.set target heap[source]) nodes := by
  intro slot hslot
  by_cases heq : target = slot
  · subst slot
    rw [Array.getElem?_set_self htarget]
    simpa [Array.getElem?_eq_getElem hsource] using hvalid source hsource
  · rw [Array.getElem?_set_ne htarget heq]
    exact hvalid slot (by simpa only [Array.size_set] using hslot)

theorem PairVecDivVHCHeapBoundedBy.set_from_slot
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode) (bound : UMonomial)
    (target source : Nat) (htarget : target < heap.size)
    (hsource : source < heap.size)
    (hbound : PairVecDivVHCHeapBoundedBy heap nodes bound) :
    PairVecDivVHCHeapBoundedBy (heap.set target heap[source]) nodes bound := by
  intro slot head mono hget hmono
  by_cases heq : target = slot
  · subst slot
    rw [Array.getElem?_set_self htarget] at hget
    exact hbound source head mono
      (by simpa [Array.getElem?_eq_getElem hsource] using hget) hmono
  · rw [Array.getElem?_set_ne htarget heq] at hget
    exact hbound slot head mono hget hmono

theorem PairVecDivVHCHeapPointersValid.push_from_slot
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode)
    (source : Nat) (hsource : source < heap.size)
    (hvalid : PairVecDivVHCHeapPointersValid heap nodes) :
    PairVecDivVHCHeapPointersValid (heap.push heap[source]) nodes := by
  intro slot hslot
  by_cases heq : slot = heap.size
  · subst slot
    simp only [Array.getElem?_push, if_pos rfl]
    simpa [Array.getElem?_eq_getElem hsource] using hvalid source hsource
  · have hold : slot < heap.size := by
      simp only [Array.size_push] at hslot
      omega
    simp only [Array.getElem?_push, if_neg heq]
    exact hvalid slot hold

theorem PairVecDivVHCHeapBoundedBy.push_from_slot
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode) (bound : UMonomial)
    (source : Nat) (hsource : source < heap.size)
    (hbound : PairVecDivVHCHeapBoundedBy heap nodes bound) :
    PairVecDivVHCHeapBoundedBy (heap.push heap[source]) nodes bound := by
  intro slot head mono hget hmono
  by_cases heq : slot = heap.size
  · subst slot
    simp only [Array.getElem?_push, if_pos rfl] at hget
    exact hbound source head mono
      (by simpa [Array.getElem?_eq_getElem hsource] using hget) hmono
  · have hold : slot < heap.size := by
      have hslot : slot < (heap.push heap[source]).size :=
        Array.getElem?_eq_some_iff.mp hget |>.1
      simp only [Array.size_push] at hslot
      omega
    simp only [Array.getElem?_push, if_neg heq] at hget
    exact hbound slot head mono hget hmono

/-- Appending a duplicate of the future last slot's parent preserves heap
order.  This is the ordered state after the first pointer copy of root bubble. -/
theorem pairVecDivVHCPush_parent_preserves_heapOrdered
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode) (parentHead : Nat)
    (hordered : PairVecDivVHCHeapOrdered heap nodes)
    (hparent : heap[pairVecDivVHCParent heap.size]? = some parentHead) :
    PairVecDivVHCHeapOrdered (heap.push parentHead) nodes := by
  apply PairVecDivVHCHeapDegreesOrderedUpTo.toHeapOrdered
  intro child hchild hpos childHead targetParentHead childMono parentMono
    hchildGet htargetParentGet hchildMono hparentMono
  by_cases hlast : child = heap.size
  · subst child
    simp only [Array.getElem?_push, if_pos rfl] at hchildGet
    have hchildHeadEq : childHead = parentHead :=
      (Option.some.inj hchildGet).symm
    subst childHead
    have hparentSlot : pairVecDivVHCParent heap.size < heap.size :=
      pairVecDivVHCParent_lt heap.size hpos
    simp only [Array.getElem?_push,
      if_neg (Nat.ne_of_lt hparentSlot)] at htargetParentGet
    rw [hparent] at htargetParentGet
    have htargetHeadEq : targetParentHead = parentHead :=
      (Option.some.inj htargetParentGet).symm
    subst targetParentHead
    rw [hchildMono] at hparentMono
    exact Nat.le_of_eq
      (congrArg UMonomial.deg (Except.ok.inj hparentMono))
  · have hchildOld : child < heap.size := by
      simp only [Array.size_push] at hchild
      omega
    have hparentOld : pairVecDivVHCParent child < heap.size :=
      Nat.lt_trans (pairVecDivVHCParent_lt child hpos) hchildOld
    simp only [Array.getElem?_push, if_neg hlast] at hchildGet
    simp only [Array.getElem?_push,
      if_neg (Nat.ne_of_lt hparentOld)] at htargetParentGet
    exact (hordered.degreesUpTo heap nodes heap.size (Nat.le_refl _)) child
      hchildOld hpos childHead targetParentHead childMono parentMono
      hchildGet htargetParentGet hchildMono hparentMono

theorem pairVecDivVHCPush_le_parent_preserves_heapOrdered
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode)
    (newHead expectedParentHead : Nat)
    (newMono parentMono : UMonomial) (hheap : 0 < heap.size)
    (hordered : PairVecDivVHCHeapOrdered heap nodes)
    (hnew : pairVecDivVHCMono newHead nodes = .ok newMono)
    (hparentGet : heap[pairVecDivVHCParent heap.size]? = some expectedParentHead)
    (hparent : pairVecDivVHCMono expectedParentHead nodes = .ok parentMono)
    (hle : newMono.deg ≤ parentMono.deg) :
    PairVecDivVHCHeapOrdered (heap.push newHead) nodes := by
  apply PairVecDivVHCHeapDegreesOrderedUpTo.toHeapOrdered
  intro child hchild hpos childHead actualParentHead childMono targetParentMono
    hchildGet hactualParentGet hchildMono htargetParentMono
  by_cases hlast : child = heap.size
  · subst child
    simp only [Array.getElem?_push, if_pos rfl] at hchildGet
    have hchildHeadEq : childHead = newHead :=
      (Option.some.inj hchildGet).symm
    subst childHead
    rw [hnew] at hchildMono
    have hchildMonoEq : childMono = newMono :=
      (Except.ok.inj hchildMono).symm
    subst childMono
    have hparentSlot : pairVecDivVHCParent heap.size < heap.size :=
      pairVecDivVHCParent_lt heap.size hheap
    simp only [Array.getElem?_push,
      if_neg (Nat.ne_of_lt hparentSlot)] at hactualParentGet
    rw [hparentGet] at hactualParentGet
    have hparentHeadEq : actualParentHead = expectedParentHead :=
      (Option.some.inj hactualParentGet).symm
    subst actualParentHead
    rw [hparent] at htargetParentMono
    have hparentMonoEq : targetParentMono = parentMono :=
      (Except.ok.inj htargetParentMono).symm
    subst targetParentMono
    exact hle
  · have hchildOld : child < heap.size := by
      simp only [Array.size_push] at hchild
      omega
    have hparentOld : pairVecDivVHCParent child < heap.size :=
      Nat.lt_trans (pairVecDivVHCParent_lt child hpos) hchildOld
    simp only [Array.getElem?_push, if_neg hlast] at hchildGet
    simp only [Array.getElem?_push,
      if_neg (Nat.ne_of_lt hparentOld)] at hactualParentGet
    exact (hordered.degreesUpTo heap nodes heap.size (Nat.le_refl _)) child
      hchildOld hpos childHead actualParentHead childMono targetParentMono
      hchildGet hactualParentGet hchildMono htargetParentMono

/-- Root-directed pointer-copying bubble preserves heap order once its current
hole state is ordered and every stored key is bounded by the inserted key. -/
theorem pairVecDivVHCBubble_to_root_preserves_heapOrdered
    (i newNode : Nat) (heap heap' : Array Nat)
    (nodes : Array PairVecDivVHCNode) (newMono : UMonomial)
    (hvalid : PairVecDivVHCHeapPointersValid heap nodes)
    (hordered : PairVecDivVHCHeapOrdered heap nodes)
    (hbound : PairVecDivVHCHeapBoundedBy heap nodes newMono)
    (hnew : pairVecDivVHCMono newNode nodes = .ok newMono)
    (hrun : pairVecDivVHCBubble i 0 newNode heap = .ok heap') :
    PairVecDivVHCHeapOrdered heap' nodes := by
  rw [pairVecDivVHCBubble] at hrun
  split at hrun <;> try contradiction
  next hi =>
    split at hrun <;> try contradiction
    next hstop =>
      split at hrun
      next heq =>
        simp only [Except.ok.injEq] at hrun
        subst i
        subst heap'
        apply PairVecDivVHCHeapDegreesOrderedUpTo.toHeapOrdered
        intro child hchild hpos childHead parentHead childMono parentMono
          hchildGet hparentGet hchildMono hparentMono
        have hchildOld : child < heap.size := by
          simpa only [Array.size_set] using hchild
        have hrootNeChild : 0 ≠ child := Nat.ne_of_lt hpos
        rw [Array.getElem?_set_ne hi hrootNeChild] at hchildGet
        by_cases hparentRoot : pairVecDivVHCParent child = 0
        · rw [hparentRoot, Array.getElem?_set_self hi] at hparentGet
          have hparentHeadEq : parentHead = newNode :=
            (Option.some.inj hparentGet).symm
          subst parentHead
          rw [hnew] at hparentMono
          have hparentMonoEq : parentMono = newMono :=
            (Except.ok.inj hparentMono).symm
          subst parentMono
          exact hbound child childHead childMono hchildGet hchildMono
        · have hrootNeParent : 0 ≠ pairVecDivVHCParent child :=
            Ne.symm hparentRoot
          rw [Array.getElem?_set_ne hi hrootNeParent] at hparentGet
          exact (hordered.degreesUpTo heap nodes heap.size (Nat.le_refl _))
            child hchildOld hpos childHead parentHead childMono parentMono
            hchildGet hparentGet hchildMono hparentMono
      next heq =>
        dsimp only at hrun
        split at hrun <;> try contradiction
        next hp =>
          have hpos : 0 < i := by omega
          let parent := pairVecDivVHCParent i
          have hnextValid : PairVecDivVHCHeapPointersValid
              (heap.set i heap[parent]) nodes := by
            exact hvalid.set_from_slot heap nodes i parent hi (by simpa [parent])
          have hnextOrdered : PairVecDivVHCHeapOrdered
              (heap.set i heap[parent]) nodes := by
            apply pairVecDivVHCSet_child_to_parent_preserves_heapOrdered
              heap nodes i heap[parent] hi hpos
            · simpa [parent, Array.getElem?_eq_getElem hp]
            · exact hvalid
            · exact hordered
          have hnextBound : PairVecDivVHCHeapBoundedBy
              (heap.set i heap[parent]) nodes newMono := by
            exact hbound.set_from_slot heap nodes newMono i parent hi
              (by simpa [parent])
          exact pairVecDivVHCBubble_to_root_preserves_heapOrdered parent newNode
            (heap.set i heap[parent]) heap' nodes newMono hnextValid
            hnextOrdered hnextBound hnew (by simpa [parent] using hrun)
termination_by i
decreasing_by
  exact pairVecDivVHCParent_lt i (by omega)

/-- `BubbleBelow` preserves max-heap order when driven by the exact comparison
trace produced by `FindAnchor`.  The climb comparison bounds the overwritten
slot's children below the inserted key, while the trace's stop comparison
bounds the inserted key below its anchor parent. -/
theorem pairVecDivVHCBubbleBelow_trace_preserves_heapOrdered
    (i anchor newNode : Nat) (heap heap' : Array Nat)
    (nodes : Array PairVecDivVHCNode) (newMono : UMonomial)
    (hvalid : PairVecDivVHCHeapPointersValid heap nodes)
    (hordered : PairVecDivVHCHeapOrdered heap nodes)
    (hnew : pairVecDivVHCMono newNode nodes = .ok newMono)
    (htrace : PairVecDivVHCFindAnchorTrace newMono.deg heap nodes i anchor)
    (hne : i ≠ anchor)
    (hrun : pairVecDivVHCBubbleBelow i anchor newNode heap = .ok heap') :
    PairVecDivVHCHeapOrdered heap' nodes := by
  cases htrace with
  | stop i head mono hhead hmono hle => exact False.elim (hne rfl)
  | climb i head anchor mono hhead hmono hlt hpos htail =>
      rw [pairVecDivVHCBubbleBelow] at hrun
      split at hrun <;> try contradiction
      next hi =>
          dsimp only at hrun
          split at hrun
          next hparentAnchor =>
            simp only [Except.ok.injEq] at hrun
            subst heap'
            rcases htail.anchor_read newMono.deg heap nodes
                (pairVecDivVHCParent i) anchor with
              ⟨anchorHead, anchorMono, hanchorHead, hanchorMono, hnewLeAnchor⟩
            have hdegrees := hordered.degreesUpTo heap nodes heap.size
              (Nat.le_refl _)
            apply PairVecDivVHCHeapDegreesOrderedUpTo.toHeapOrdered
            have hsetOrdered : PairVecDivVHCHeapDegreesOrderedUpTo heap.size
                (heap.set i newNode) nodes := by
              apply hdegrees.set_parent heap.size i newNode heap nodes newMono
                hi hnew
              · intro parentHead parentMono hparentGet hparentMono
                rw [hparentAnchor, hanchorHead] at hparentGet
                have hparentHeadEq : parentHead = anchorHead :=
                  (Option.some.inj hparentGet).symm
                subst parentHead
                rw [hanchorMono] at hparentMono
                have hparentMonoEq : parentMono = anchorMono :=
                  (Except.ok.inj hparentMono).symm
                subst parentMono
                exact hnewLeAnchor
              · intro child hchild hchildPos hchildParent childHead childMono
                  hchildGet hchildMono
                have hchildLeOld : childMono.deg ≤ mono.deg :=
                  hdegrees child hchild hchildPos childHead head childMono mono
                    hchildGet (by simpa [hchildParent] using hhead)
                    hchildMono hmono
                exact Nat.le_trans hchildLeOld (Nat.le_of_lt hlt)
            simpa only [Array.size_set] using hsetOrdered
          next hparentAnchor =>
            split at hrun <;> try contradiction
            next hp =>
              let parent := pairVecDivVHCParent i
              have hparentLt : parent < i := pairVecDivVHCParent_lt i hpos
              have htail' : PairVecDivVHCFindAnchorTrace newMono.deg
                  (heap.set i heap[parent]) nodes parent anchor := by
                exact htail.set_above newMono.deg heap nodes parent anchor i
                  heap[parent] hi hparentLt
              have hnextValid : PairVecDivVHCHeapPointersValid
                  (heap.set i heap[parent]) nodes :=
                hvalid.set_from_slot heap nodes i parent hi (by simpa [parent])
              have hnextOrdered : PairVecDivVHCHeapOrdered
                  (heap.set i heap[parent]) nodes := by
                apply pairVecDivVHCSet_child_to_parent_preserves_heapOrdered
                  heap nodes i heap[parent] hi hpos
                · simpa [parent, Array.getElem?_eq_getElem hp]
                · exact hvalid
                · exact hordered
              exact pairVecDivVHCBubbleBelow_trace_preserves_heapOrdered parent
                anchor newNode (heap.set i heap[parent]) heap' nodes newMono
                hnextValid hnextOrdered hnew htail'
                (by simpa [parent] using hparentAnchor) (by simpa [parent] using hrun)
termination_by i
decreasing_by
  apply pairVecDivVHCParent_lt
  omega

/-- Complete non-root bubble from the appended last slot, driven by the trace
of the preceding generated `FindAnchor` call. -/
theorem pairVecDivVHCBubbleBelow_push_trace_preserves_heapOrdered
    (anchor newNode : Nat) (heap heap' : Array Nat)
    (nodes : Array PairVecDivVHCNode) (newMono : UMonomial)
    (hheap : 0 < heap.size)
    (hvalid : PairVecDivVHCHeapPointersValid heap nodes)
    (hordered : PairVecDivVHCHeapOrdered heap nodes)
    (hnew : pairVecDivVHCMono newNode nodes = .ok newMono)
    (htrace : PairVecDivVHCFindAnchorTrace newMono.deg heap nodes
      (pairVecDivVHCParent heap.size) anchor)
    (hrun : pairVecDivVHCBubbleBelow heap.size anchor newNode
      (heap.push newNode) = .ok heap') :
    PairVecDivVHCHeapOrdered heap' nodes := by
  let first := pairVecDivVHCParent heap.size
  have hfirst : first < heap.size := pairVecDivVHCParent_lt heap.size hheap
  rw [pairVecDivVHCBubbleBelow] at hrun
  split at hrun <;> try contradiction
  next hi =>
    dsimp only at hrun
    split at hrun
    next hfirstAnchor =>
      simp only [Array.set_push, lt_self_iff_false, ↓reduceDIte,
        Except.ok.injEq] at hrun
      subst heap'
      rcases htrace.anchor_read newMono.deg heap nodes first anchor with
        ⟨anchorHead, anchorMono, hanchorGet, hanchorMono, hle⟩
      have hparentGet : heap[pairVecDivVHCParent heap.size]? =
          some anchorHead := by
        simpa [first, hfirstAnchor] using hanchorGet
      exact pairVecDivVHCPush_le_parent_preserves_heapOrdered heap nodes
        newNode anchorHead newMono anchorMono hheap hordered hnew hparentGet
        hanchorMono hle
    next hfirstAnchor =>
      split at hrun <;> try contradiction
      next hp =>
        rw [Array.getElem_push_lt hfirst] at hrun
        have hrun' : pairVecDivVHCBubbleBelow first anchor newNode
            (heap.push heap[first]) = .ok heap' := by
          simpa [first, Array.set_push, hfirst] using hrun
        have htracePush : PairVecDivVHCFindAnchorTrace newMono.deg
            (heap.push heap[first]) nodes first anchor :=
          htrace.push newMono.deg heap nodes first anchor heap[first] hfirst
        exact pairVecDivVHCBubbleBelow_trace_preserves_heapOrdered first anchor
          newNode (heap.push heap[first]) heap' nodes newMono
          (hvalid.push_from_slot heap nodes first hfirst)
          (pairVecDivVHCPush_parent_preserves_heapOrdered heap nodes heap[first]
            hordered (by simpa [first, Array.getElem?_eq_getElem hfirst]))
          hnew htracePush (by simpa [first] using hfirstAnchor) hrun'

theorem pairVecDivVHCHeapOrdered_slot_le_root
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode)
    (hvalid : PairVecDivVHCHeapPointersValid heap nodes)
    (hordered : PairVecDivVHCHeapOrdered heap nodes)
    (slot : Nat) (hslot : slot < heap.size)
    (mono rootMono : UMonomial)
    (hmono : pairVecDivVHCMono heap[slot] nodes = .ok mono)
    (hroot : pairVecDivVHCMono heap[0] nodes = .ok rootMono) :
    mono.deg ≤ rootMono.deg := by
  induction slot using Nat.strong_induction_on generalizing mono with
  | h slot ih =>
      by_cases hzero : slot = 0
      · subst slot
        rw [hroot] at hmono
        simp only [Except.ok.injEq] at hmono
        subst mono
        exact Nat.le_refl _
      · have hpos : 0 < slot := Nat.pos_of_ne_zero hzero
        let parent := pairVecDivVHCParent slot
        have hparentLt : parent < slot := pairVecDivVHCParent_lt slot hpos
        have hparentSlot : parent < heap.size := Nat.lt_trans hparentLt hslot
        rcases hvalid parent hparentSlot with
          ⟨parentIndex, parentNode, parentMono, hheapParent, hparentNode,
            hparentActive⟩
        have hparentIndexEq : heap[parent] = parentIndex := by
          rw [Array.getElem?_eq_getElem hparentSlot] at hheapParent
          exact Option.some.inj hheapParent
        have hparentRun : pairVecDivVHCMono heap[parent] nodes =
            .ok parentMono := by
          apply (pairVecDivVHCMono_eq_ok_iff heap[parent] nodes
            parentMono).mpr
          rw [hparentIndexEq]
          exact ⟨parentNode, hparentNode, hparentActive⟩
        rcases (pairVecDivVHCMono_eq_ok_iff heap[slot] nodes mono).mp hmono with
          ⟨childNode, hchildNode, hchildActive⟩
        have hchildMap :
            (nodes[heap[slot]]?.map PairVecDivVHCNode.mono).join =
              some mono := by
          rw [hchildNode]
          simp [hchildActive]
        have hparentMap :
            (nodes[heap[parent]]?.map PairVecDivVHCNode.mono).join =
              some parentMono := by
          rw [hparentIndexEq, hparentNode]
          simp [hparentActive]
        have hstep := hordered slot parent hslot rfl hpos heap[slot]
          heap[parent] mono parentMono
          (Array.getElem?_eq_getElem hslot)
          (Array.getElem?_eq_getElem hparentSlot) hchildMap hparentMap
        exact Nat.le_trans hstep
          (ih parent hparentLt hparentSlot parentMono hparentRun hroot)

theorem pairVecDivVHCHeapOrdered_child_le_parent
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode)
    (hordered : PairVecDivVHCHeapOrdered heap nodes)
    (child : Nat) (hchild : child < heap.size) (hpos : 0 < child)
    (hparent : pairVecDivVHCParent child < heap.size)
    (childMono parentMono : UMonomial)
    (hchildMono : pairVecDivVHCMono heap[child] nodes = .ok childMono)
    (hparentMono : pairVecDivVHCMono heap[pairVecDivVHCParent child] nodes =
      .ok parentMono) :
    childMono.deg ≤ parentMono.deg := by
  rcases (pairVecDivVHCMono_eq_ok_iff heap[child] nodes childMono).mp
      hchildMono with ⟨childNode, hchildNode, hchildActive⟩
  rcases (pairVecDivVHCMono_eq_ok_iff heap[pairVecDivVHCParent child] nodes
      parentMono).mp hparentMono with
    ⟨parentNode, hparentNode, hparentActive⟩
  have hchildMap :
      (nodes[heap[child]]?.map PairVecDivVHCNode.mono).join =
        some childMono := by
    rw [hchildNode]
    simp [hchildActive]
  have hparentMap :
      (nodes[heap[pairVecDivVHCParent child]]?.map
        PairVecDivVHCNode.mono).join = some parentMono := by
    rw [hparentNode]
    simp [hparentActive]
  exact hordered child (pairVecDivVHCParent child) hchild rfl hpos
    heap[child] heap[pairVecDivVHCParent child] childMono parentMono
    (Array.getElem?_eq_getElem hchild) (Array.getElem?_eq_getElem hparent)
    hchildMap hparentMap

/-- Every non-root heap slot lies below one of the root's two direct children,
and its active degree is bounded by that child's degree. -/
theorem pairVecDivVHCHeapOrdered_slot_le_root_child
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode)
    (hvalid : PairVecDivVHCHeapPointersValid heap nodes)
    (hordered : PairVecDivVHCHeapOrdered heap nodes)
    (slot : Nat) (hslot : slot < heap.size) (hpos : 0 < slot)
    (mono : UMonomial)
    (hmono : pairVecDivVHCMono heap[slot] nodes = .ok mono) :
    ∃ child head childMono,
      (child = 1 ∨ child = 2) ∧ heap[child]? = some head ∧
        pairVecDivVHCMono head nodes = .ok childMono ∧
        mono.deg ≤ childMono.deg := by
  induction slot using Nat.strong_induction_on generalizing mono with
  | h slot ih =>
      let parent := pairVecDivVHCParent slot
      have hparentLt : parent < slot := pairVecDivVHCParent_lt slot hpos
      have hparentSlot : parent < heap.size := Nat.lt_trans hparentLt hslot
      rcases hvalid parent hparentSlot with
        ⟨parentIndex, parentNode, parentMono, hheapParent, hparentNode,
          hparentActive⟩
      have hparentIndexEq : heap[parent] = parentIndex := by
        rw [Array.getElem?_eq_getElem hparentSlot] at hheapParent
        exact Option.some.inj hheapParent
      have hparentRun : pairVecDivVHCMono heap[parent] nodes =
          .ok parentMono := by
        apply (pairVecDivVHCMono_eq_ok_iff heap[parent] nodes parentMono).mpr
        rw [hparentIndexEq]
        exact ⟨parentNode, hparentNode, hparentActive⟩
      have hstep : mono.deg ≤ parentMono.deg :=
        pairVecDivVHCHeapOrdered_child_le_parent heap nodes hordered slot
          hslot hpos hparentSlot mono parentMono hmono
          (by simpa [parent] using hparentRun)
      by_cases hrootParent : parent = 0
      · have hchild : slot = 1 ∨ slot = 2 := by
          have hrootParent' : (slot - 1) / 2 = 0 := by
            simpa [parent, pairVecDivVHCParent] using hrootParent
          have hsmall : slot - 1 < 2 := by
            apply (Nat.div_eq_zero_iff_lt (by omega : 0 < 2)).mp
            exact hrootParent'
          omega
        exact ⟨slot, heap[slot], mono, hchild,
          Array.getElem?_eq_getElem hslot, hmono, Nat.le_refl _⟩
      · have hparentPos : 0 < parent := Nat.pos_of_ne_zero hrootParent
        rcases ih parent hparentLt hparentSlot hparentPos parentMono hparentRun
          with ⟨child, head, childMono, hchild, hchildGet, hchildRun,
            hparentLe⟩
        exact ⟨child, head, childMono, hchild, hchildGet, hchildRun,
          Nat.le_trans hstep hparentLe⟩

theorem pairVecDivVHCBubble_new_root_bounds_all
    (newNode : Nat) (heap heap' : Array Nat)
    (nodes : Array PairVecDivVHCNode) (newMono rootMono : UMonomial)
    (hheap : 0 < heap.size)
    (hvalid : PairVecDivVHCHeapPointersValid heap nodes)
    (hordered : PairVecDivVHCHeapOrdered heap nodes)
    (hnew : pairVecDivVHCMono newNode nodes = .ok newMono)
    (hroot : pairVecDivVHCMono heap[0] nodes = .ok rootMono)
    (hrootLe : rootMono.deg ≤ newMono.deg)
    (hrun : pairVecDivVHCBubble heap.size 0 newNode
      (heap.push newNode) = .ok heap') :
    ∀ (slot head : Nat) (mono : UMonomial),
      heap'[slot]? = some head → pairVecDivVHCMono head nodes = .ok mono →
        mono.deg ≤ newMono.deg := by
  have hfrom := pairVecDivVHCBubble_valuesFrom heap.size 0 newNode
    (heap.push newNode) heap' (heap.push newNode)
    (pairVecDivVHCValuesFrom_refl _) ⟨heap.size, by simp⟩ hrun
  intro slot head mono hget hmono
  rcases hfrom slot head hget with ⟨sourceSlot, hsource⟩
  have hsourceBound : sourceSlot < (heap.push newNode).size := by
    by_contra hnot
    rw [Array.getElem?_eq_none (by omega)] at hsource
    contradiction
  by_cases hold : sourceSlot < heap.size
  · have hsourceOld : heap[sourceSlot]? = some head := by
      have hne : sourceSlot ≠ heap.size := by omega
      simp only [Array.getElem?_push, hne, ↓reduceIte] at hsource
      exact hsource
    have hheadEq : heap[sourceSlot] = head := by
      rw [Array.getElem?_eq_getElem hold] at hsourceOld
      exact Option.some.inj hsourceOld
    have hslotMono : pairVecDivVHCMono heap[sourceSlot] nodes = .ok mono := by
      simpa [hheadEq] using hmono
    have hleRoot := pairVecDivVHCHeapOrdered_slot_le_root heap nodes hvalid
      hordered sourceSlot hold mono rootMono hslotMono hroot
    exact Nat.le_trans hleRoot hrootLe
  · have heq : sourceSlot = heap.size := by
      simp only [Array.size_push] at hsourceBound
      omega
    subst sourceSlot
    simp only [Array.getElem?_push, ↓reduceIte, if_pos rfl] at hsource
    have hheadEq : head = newNode := (Option.some.inj hsource).symm
    subst head
    rw [hnew] at hmono
    exact Nat.le_of_eq (congrArg UMonomial.deg
      (Except.ok.inj hmono).symm)

/-- The complete generated greater-root bubble, including its initially
unordered appended key, returns a max-heap. -/
theorem pairVecDivVHCBubble_new_root_preserves_heapOrdered
    (newNode : Nat) (heap heap' : Array Nat)
    (nodes : Array PairVecDivVHCNode) (newMono rootMono : UMonomial)
    (hheap : 0 < heap.size)
    (hvalid : PairVecDivVHCHeapPointersValid heap nodes)
    (hordered : PairVecDivVHCHeapOrdered heap nodes)
    (hnew : pairVecDivVHCMono newNode nodes = .ok newMono)
    (hroot : pairVecDivVHCMono heap[0] nodes = .ok rootMono)
    (hrootLe : rootMono.deg ≤ newMono.deg)
    (hrun : pairVecDivVHCBubble heap.size 0 newNode
      (heap.push newNode) = .ok heap') :
    PairVecDivVHCHeapOrdered heap' nodes := by
  have holdBound : PairVecDivVHCHeapBoundedBy heap nodes newMono := by
    intro slot head mono hget hmono
    have hslot : slot < heap.size := Array.getElem?_eq_some_iff.mp hget |>.1
    have hheadEq : heap[slot] = head := by
      rw [Array.getElem?_eq_getElem hslot] at hget
      exact Option.some.inj hget
    have hslotMono : pairVecDivVHCMono heap[slot] nodes = .ok mono := by
      simpa [hheadEq] using hmono
    exact Nat.le_trans
      (pairVecDivVHCHeapOrdered_slot_le_root heap nodes hvalid hordered slot
        hslot mono rootMono hslotMono hroot) hrootLe
  let parent := pairVecDivVHCParent heap.size
  have hparentOld : parent < heap.size := by
    exact pairVecDivVHCParent_lt heap.size hheap
  have hparentGet : heap[parent]? = some heap[parent] :=
    Array.getElem?_eq_getElem hparentOld
  rw [pairVecDivVHCBubble] at hrun
  split at hrun <;> try contradiction
  next hi =>
    split at hrun <;> try contradiction
    next hstop =>
      split at hrun
      next heq => omega
      next heq =>
        dsimp only at hrun
        split at hrun <;> try contradiction
        next hp =>
          rw [Array.getElem_push_lt hparentOld] at hrun
          have hrun' : pairVecDivVHCBubble parent 0 newNode
              (heap.push heap[parent]) = .ok heap' := by
            simpa [parent, Array.set_push, hparentOld] using hrun
          apply pairVecDivVHCBubble_to_root_preserves_heapOrdered parent
            newNode (heap.push heap[parent]) heap' nodes newMono
          · exact hvalid.push_from_slot heap nodes parent hparentOld
          · exact pairVecDivVHCPush_parent_preserves_heapOrdered heap nodes
              heap[parent] hordered (by simpa [parent] using hparentGet)
          · exact holdBound.push_from_slot heap nodes newMono parent hparentOld
          · exact hnew
          · exact hrun'

theorem pairVecDivVHCExtract_root_dominates
    (heap heap' : Array Nat) (sourceNodes nodes : Array PairVecDivVHCNode)
    (hvalid : PairVecDivVHCHeapPointersValid heap sourceNodes)
    (hordered : PairVecDivVHCHeapOrdered heap sourceNodes)
    (hunique : ∀ (left right head : Nat), heap[left]? = some head →
      heap[right]? = some head → left = right)
    (hsame : ∀ (slot head : Nat), 0 < slot →
      heap[slot]? = some head → ∀ mono,
        pairVecDivVHCMono head sourceNodes = .ok mono ↔
          pairVecDivVHCMono head nodes = .ok mono)
    (hresult : 0 < heap'.size)
    (hrun : pairVecDivVHCExtract heap nodes = .ok heap')
    (slot : Nat) (hslot : slot < heap'.size) (mono rootMono : UMonomial)
    (hmono : pairVecDivVHCMono heap'[slot] nodes = .ok mono)
    (hroot : pairVecDivVHCMono heap'[0] nodes = .ok rootMono) :
    mono.deg ≤ rootMono.deg := by
  have hsourceSize := pairVecDivVHCExtract_size heap heap' nodes hrun
  have hsourceNonempty : 0 < heap.size := by omega
  have hrootOnly : PairVecDivVHCOnlyAt heap heap[0] 0 := by
    intro i hget
    exact hunique i 0 heap[0] hget
      (Array.getElem?_eq_getElem hsourceNonempty)
  by_cases htwo : heap.size = 2
  · have htargetZero : slot = 0 := by omega
    subst slot
    rw [hroot] at hmono
    exact Nat.le_of_eq (congrArg UMonomial.deg
      (Except.ok.inj hmono).symm)
  · have hmany : 2 < heap.size := by omega
    rcases hvalid 1 (by omega) with
      ⟨leftHead, leftNode, leftMono, hleftHeap, hleftNode, hleftActive⟩
    rcases hvalid 2 (by omega) with
      ⟨rightHead, rightNode, rightMono, hrightHeap, hrightNode,
        hrightActive⟩
    rcases hvalid (heap.size - 1) (by omega) with
      ⟨lastHead, lastNode, lastMono, hlastHeap, hlastNode, hlastActive⟩
    have hleftHead : heap[1] = leftHead := by
      rw [Array.getElem?_eq_getElem (by omega)] at hleftHeap
      exact Option.some.inj hleftHeap
    have hrightHead : heap[2] = rightHead := by
      rw [Array.getElem?_eq_getElem (by omega)] at hrightHeap
      exact Option.some.inj hrightHeap
    have hlastHead : heap[heap.size - 1] = lastHead := by
      rw [Array.getElem?_eq_getElem (by omega)] at hlastHeap
      exact Option.some.inj hlastHeap
    have hleftSource : pairVecDivVHCMono heap[1] sourceNodes =
        .ok leftMono := by
      apply (pairVecDivVHCMono_eq_ok_iff heap[1] sourceNodes leftMono).mpr
      rw [hleftHead]
      exact ⟨leftNode, hleftNode, hleftActive⟩
    have hrightSource : pairVecDivVHCMono heap[2] sourceNodes =
        .ok rightMono := by
      apply (pairVecDivVHCMono_eq_ok_iff heap[2] sourceNodes rightMono).mpr
      rw [hrightHead]
      exact ⟨rightNode, hrightNode, hrightActive⟩
    have hlastSource : pairVecDivVHCMono heap[heap.size - 1] sourceNodes =
        .ok lastMono := by
      apply (pairVecDivVHCMono_eq_ok_iff heap[heap.size - 1] sourceNodes
        lastMono).mpr
      rw [hlastHead]
      exact ⟨lastNode, hlastNode, hlastActive⟩
    have hleftRun := (hsame 1 heap[1] (by omega)
      (Array.getElem?_eq_getElem (by omega)) leftMono).mp hleftSource
    have hrightRun := (hsame 2 heap[2] (by omega)
      (Array.getElem?_eq_getElem (by omega)) rightMono).mp hrightSource
    have hlastRun := (hsame (heap.size - 1) heap[heap.size - 1] (by omega)
      (Array.getElem?_eq_getElem (by omega)) lastMono).mp hlastSource
    rcases pairVecDivVHCExtract_root_dominates_candidates heap heap' nodes
        leftMono rightMono lastMono hmany hleftRun hrightRun hlastRun hrun with
      ⟨newRootHead, newRootMono, hnewRoot, hnewRootRun, hleftLe,
        hrightLe, hlastLe⟩
    have hnewRootEq : heap'[0] = newRootHead := by
      rw [Array.getElem?_eq_getElem hresult] at hnewRoot
      exact Option.some.inj hnewRoot
    have hrootMonoEq : rootMono = newRootMono := by
      rw [hnewRootEq, hnewRootRun] at hroot
      exact (Except.ok.inj hroot).symm
    subst newRootMono
    have htargetGet : heap'[slot]? = some heap'[slot] :=
      Array.getElem?_eq_getElem hslot
    rcases pairVecDivVHCExtract_valuesFrom heap heap' nodes hrun slot
        heap'[slot] htargetGet with ⟨oldSlot, holdGet⟩
    have holdSlot : oldSlot < heap.size := by
      by_contra hnot
      rw [Array.getElem?_eq_none (by omega)] at holdGet
      contradiction
    have holdHead : heap[oldSlot] = heap'[slot] := by
      rw [Array.getElem?_eq_getElem holdSlot] at holdGet
      exact Option.some.inj holdGet
    have holdPos : 0 < oldSlot := by
      by_contra hnot
      have holdZero : oldSlot = 0 := by omega
      subst oldSlot
      have hrootHead : heap'[slot] = heap[0] := by
        rw [Array.getElem?_eq_getElem hsourceNonempty] at holdGet
        exact (Option.some.inj holdGet).symm
      have hexcludes := pairVecDivVHCExtract_excludes_root heap heap' nodes
        hsourceNonempty hrootOnly hrun
      exact hexcludes slot (by simpa [hrootHead] using htargetGet)
    have holdRunCurrent : pairVecDivVHCMono heap[oldSlot] nodes =
        .ok mono := by
      rw [holdHead]
      exact hmono
    have holdRun := (hsame oldSlot heap[oldSlot] holdPos
      (Array.getElem?_eq_getElem holdSlot) mono).mpr holdRunCurrent
    rcases pairVecDivVHCHeapOrdered_slot_le_root_child heap sourceNodes hvalid
        hordered oldSlot holdSlot holdPos mono holdRun with
      ⟨child, childHead, childMono, hchild, hchildGet, hchildRun, hmonoLe⟩
    rcases hchild with rfl | rfl
    · rw [Array.getElem?_eq_getElem (by omega)] at hchildGet
      have hheadEq := Option.some.inj hchildGet
      rw [← hheadEq, hleftSource] at hchildRun
      have hmonoEq := Except.ok.inj hchildRun
      subst childMono
      exact Nat.le_trans hmonoLe hleftLe
    · rw [Array.getElem?_eq_getElem (by omega)] at hchildGet
      have hheadEq := Option.some.inj hchildGet
      rw [← hheadEq, hrightSource] at hchildRun
      have hmonoEq := Except.ok.inj hchildRun
      subst childMono
      exact Nat.le_trans hmonoLe hrightLe

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
          if horder : nodeIndex = resetH then
            /- The C++ preincrement leaves `v1_ptr` one-past-end before
            placing this node in the `reset_h` prefix.  The checked equality
            exposes the row-order invariant on which the source's bare
            `++reset_h` relies.  Its stale `mono`/`next` fields are
            unobservable until activation rewrites both. -/
            let node' := { node with
              quotientIndex := quotientIndex'
              mono := none
              next := none }
            .ok (k', nodes.set nodeIndex node', lin, resetH + 1, node.next)
          else
            .error .assertionFailure
        else
          .error .assertionFailure
      else
        .error .assertionFailure
    else
      .error .assertionFailure
  else
    .error .assertionFailure

theorem pairVecDivVHCConsumeNode_preserves_cursorPrefixAbove
    (this : DenseUPolyZp) (degree nodeIndex : Nat) (k k' : UInt64)
    (nodes nodes' : Array PairVecDivVHCNode) (lin lin' : Array Nat)
    (resetH resetH' : Nat) (next : Option Nat)
    (quotient divisor : SparsePolyZp) (node : PairVecDivVHCNode)
    (hprefix : PairVecDivVHCCursorPrefixAbove degree nodes quotient divisor)
    (hget : nodes[nodeIndex]? = some node)
    (hdenotes : PairVecDivVHCNodeDenotes quotient divisor node)
    (hmono : node.mono = some ⟨degree⟩)
    (hrun : pairVecDivVHCConsumeNode this nodeIndex k nodes lin resetH
      quotient divisor = .ok (k', nodes', lin', resetH', next)) :
    PairVecDivVHCCursorPrefixAbove degree nodes' quotient divisor := by
  unfold pairVecDivVHCConsumeNode at hrun
  split at hrun <;> try contradiction
  next hn =>
    rw [Array.getElem?_eq_getElem hn] at hget
    have hnodeEq : nodes[nodeIndex] = node := Option.some.inj hget
    dsimp only at hrun
    split at hrun <;> try contradiction
    next hq =>
      split at hrun <;> try contradiction
      next hd =>
        split at hrun
        next hadvance =>
          simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
          rcases hrun with ⟨rfl, rfl, rfl, rfl, rfl⟩
          apply PairVecDivVHCCursorPrefixAbove.set_advance degree nodeIndex
            nodes quotient divisor node _ hn
            (by simpa [hnodeEq] using Array.getElem?_eq_getElem hn) hprefix
            (by simp [hnodeEq]) (by simp [hnodeEq])
          intro quotientTerm divisorTerm hquotient hdivisor
          rw [hdenotes.product_degree_eq quotient divisor node
            degree hmono quotientTerm divisorTerm hquotient hdivisor]
        next hadvance =>
          split at hrun <;> try contradiction
          next hexhausted =>
            split at hrun <;> try contradiction
            next horder =>
              simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
              rcases hrun with ⟨rfl, rfl, rfl, rfl, rfl⟩
              apply PairVecDivVHCCursorPrefixAbove.set_advance degree nodeIndex
                nodes quotient divisor node _ hn
                (by simpa [hnodeEq] using Array.getElem?_eq_getElem hn) hprefix
                (by simp [hnodeEq]) (by simp [hnodeEq])
              intro quotientTerm divisorTerm hquotient hdivisor
              rw [hdenotes.product_degree_eq quotient divisor node
                degree hmono quotientTerm divisorTerm hquotient hdivisor]

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
            split at hrun <;> try contradiction
            next horder =>
              simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
              rcases hrun with ⟨rfl, rfl, rfl, rfl, rfl⟩
              right
              simp

/-- A consumed node is never lost from the source cursor partition: an
advanced cursor is appended to `lin`, while an exhausted cursor is exactly
the old `reset_h` slot and therefore enters the enlarged reset prefix.  Both
destination regions grow monotonically. -/
theorem pairVecDivVHCConsumeNode_location_progress
    (this : DenseUPolyZp) (nodeIndex : Nat) (k k' : UInt64)
    (nodes nodes' : Array PairVecDivVHCNode) (lin lin' : Array Nat)
    (resetH resetH' : Nat) (next : Option Nat)
    (quotient divisor : SparsePolyZp)
    (hrun : pairVecDivVHCConsumeNode this nodeIndex k nodes lin resetH
      quotient divisor = .ok (k', nodes', lin', resetH', next)) :
    lin.toList.toFinset ⊆ lin'.toList.toFinset ∧
      resetH ≤ resetH' ∧
      (nodeIndex ∈ lin'.toList.toFinset ∨ nodeIndex < resetH') := by
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
          refine ⟨?_, Nat.le_refl _, Or.inl ?_⟩
          · intro i hi
            simp only [List.mem_toFinset, Array.toList_push,
              List.mem_append, List.mem_singleton] at hi ⊢
            exact Or.inl hi
          · simp
        next hadvance =>
          split at hrun <;> try contradiction
          next hexhausted =>
            split at hrun <;> try contradiction
            next horder =>
              simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
              rcases hrun with ⟨rfl, rfl, rfl, rfl, rfl⟩
              refine ⟨Finset.Subset.rfl, by omega, Or.inr ?_⟩
              omega

theorem pairVecDivVHCConsumeNode_next_of_success
    (this : DenseUPolyZp) (nodeIndex : Nat) (k k' : UInt64)
    (nodes nodes' : Array PairVecDivVHCNode) (lin lin' : Array Nat)
    (resetH resetH' : Nat) (next : Option Nat)
    (quotient divisor : SparsePolyZp) (hn : nodeIndex < nodes.size)
    (hrun : pairVecDivVHCConsumeNode this nodeIndex k nodes lin resetH
      quotient divisor = .ok (k', nodes', lin', resetH', next)) :
    next = nodes[nodeIndex].next := by
  unfold pairVecDivVHCConsumeNode at hrun
  simp only [hn, ↓reduceDIte] at hrun
  split at hrun <;> try contradiction
  next hq =>
    split at hrun <;> try contradiction
    next hd =>
      split at hrun
      next hadvance =>
        simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
        exact hrun.2.2.2.2.symm
      next hadvance =>
        split at hrun <;> try contradiction
        next hexhausted =>
          split at hrun <;> try contradiction
          next horder =>
            simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
            exact hrun.2.2.2.2.symm

theorem pairVecDivVHCConsumeNode_nodes_size
    (this : DenseUPolyZp) (nodeIndex : Nat) (k k' : UInt64)
    (nodes nodes' : Array PairVecDivVHCNode) (lin lin' : Array Nat)
    (resetH resetH' : Nat) (next : Option Nat)
    (quotient divisor : SparsePolyZp)
    (hrun : pairVecDivVHCConsumeNode this nodeIndex k nodes lin resetH
      quotient divisor = .ok (k', nodes', lin', resetH', next)) :
    nodes'.size = nodes.size := by
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
          simp
        next hadvance =>
          split at hrun <;> try contradiction
          next hexhausted =>
            split at hrun <;> try contradiction
            next horder =>
              simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
              rcases hrun with ⟨rfl, rfl, rfl, rfl, rfl⟩
              simp

theorem pairVecDivVHCConsumeNode_coefficient_reduced
    (this : DenseUPolyZp) (nodeIndex : Nat) (k k' : UInt64)
    (nodes nodes' : Array PairVecDivVHCNode) (lin lin' : Array Nat)
    (resetH resetH' : Nat) (next : Option Nat)
    (quotient divisor : SparsePolyZp) (hp : this._p ≠ 0)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hcanonical : SparsePolyZp.Canonical this._p.toNat quotient)
    (hk : k.toNat < this._p.toNat)
    (hrun : pairVecDivVHCConsumeNode this nodeIndex k nodes lin resetH
      quotient divisor = .ok (k', nodes', lin', resetH', next)) :
    k'.toNat < this._p.toNat := by
  unfold pairVecDivVHCConsumeNode at hrun
  split at hrun <;> try contradiction
  next hn =>
    dsimp only at hrun
    split at hrun <;> try contradiction
    next hq =>
      split at hrun <;> try contradiction
      next hd =>
        let product := Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this
          quotient[nodes[nodeIndex].quotientIndex].2.val
          divisor[nodes[nodeIndex].divisorIndex].2.val
        have hquotientMem : quotient[nodes[nodeIndex].quotientIndex] ∈
            quotient.toList := Array.getElem_mem_toList hq
        have ha := (hcanonical.1 quotient[nodes[nodeIndex].quotientIndex]
          hquotientMem).2
        have hproductNat :=
          CLPoly.Impl.StrictWordArithmetic.nmod_mul_ir_correct_of_configured
            this quotient[nodes[nodeIndex].quotientIndex].2.val
              divisor[nodes[nodeIndex].divisorIndex].2.val hcfg ha
        change product.toNat = _ at hproductNat
        have hpNat : 0 < this._p.toNat := by
          by_contra hzero
          have : this._p.toNat = 0 := by omega
          exact hp (UInt64.toNat_inj.mp (by simpa using this))
        have hproduct : product < this._p := by
          rw [UInt64.lt_iff_toNat_lt, hproductNat]
          exact Nat.mod_lt _ hpNat
        have hkWord : k < this._p := by
          simpa [UInt64.lt_iff_toNat_lt] using hk
        split at hrun
        next hadvance =>
          simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
          rcases hrun with ⟨rfl, rfl, rfl, rfl, rfl⟩
          have hlt := CLPoly.Impl.StrictPolyAddSubRefinement.nmodSub_lt this k
            product hp hkWord hproduct
          simpa [pairVecDivSubmulIR, UInt64.lt_iff_toNat_lt] using hlt
        next hadvance =>
          split at hrun <;> try contradiction
          next hexhausted =>
            split at hrun <;> try contradiction
            next horder =>
              simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
              rcases hrun with ⟨rfl, rfl, rfl, rfl, rfl⟩
              have hlt :=
                CLPoly.Impl.StrictPolyAddSubRefinement.nmodSub_lt this k
                  product hp hkWord hproduct
              simpa [pairVecDivSubmulIR, UInt64.lt_iff_toNat_lt] using hlt

theorem pairVecDivVHCConsumeNode_coefficient_toZMod
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (nodeIndex : Nat) (k k' : UInt64)
    (nodes nodes' : Array PairVecDivVHCNode) (lin lin' : Array Nat)
    (resetH resetH' : Nat) (next : Option Nat)
    (quotient divisor : SparsePolyZp)
    (hn : nodeIndex < nodes.size)
    (hq : nodes[nodeIndex].quotientIndex < quotient.size)
    (hd : nodes[nodeIndex].divisorIndex < divisor.size)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hcanonical : SparsePolyZp.Canonical this._p.toNat quotient)
    (hk : k.toNat < this._p.toNat)
    (hrun : pairVecDivVHCConsumeNode this nodeIndex k nodes lin resetH
      quotient divisor = .ok (k', nodes', lin', resetH', next)) :
    (k'.toNat : ZMod this._p.toNat) =
      (k.toNat : ZMod this._p.toNat) -
        Zp.toZMod this._p.toNat
            quotient[nodes[nodeIndex].quotientIndex].2 *
          Zp.toZMod this._p.toNat
            divisor[nodes[nodeIndex].divisorIndex].2 := by
  unfold pairVecDivVHCConsumeNode at hrun
  simp only [hn, ↓reduceDIte] at hrun
  simp only [hq, hd, ↓reduceDIte] at hrun
  have hquotientMem : quotient[nodes[nodeIndex].quotientIndex] ∈
      quotient.toList := Array.getElem_mem_toList hq
  have ha := (hcanonical.1 quotient[nodes[nodeIndex].quotientIndex]
    hquotientMem).2
  have hsubmul := pairVecDivSubmulIR_toZMod this k
    quotient[nodes[nodeIndex].quotientIndex].2.val
    divisor[nodes[nodeIndex].divisorIndex].2.val hcfg hk ha
  split at hrun
  next hadvance =>
    simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
    rcases hrun with ⟨rfl, rfl, rfl, rfl, rfl⟩
    simpa [Zp.toZMod] using hsubmul
  next hadvance =>
    split at hrun <;> try contradiction
    next hexhausted =>
      split at hrun <;> try contradiction
      next horder =>
        simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
        rcases hrun with ⟨rfl, rfl, rfl, rfl, rfl⟩
        simpa [Zp.toZMod] using hsubmul

theorem pairVecDivVHCConsumeNode_indices
    (this : DenseUPolyZp) (nodeIndex : Nat) (k k' : UInt64)
    (nodes nodes' : Array PairVecDivVHCNode) (lin lin' : Array Nat)
    (resetH resetH' : Nat) (next : Option Nat)
    (quotient divisor : SparsePolyZp)
    (hrun : pairVecDivVHCConsumeNode this nodeIndex k nodes lin resetH
      quotient divisor = .ok (k', nodes', lin', resetH', next)) :
    ∃ hn : nodeIndex < nodes.size,
      nodes[nodeIndex].quotientIndex < quotient.size ∧
        nodes[nodeIndex].divisorIndex < divisor.size := by
  unfold pairVecDivVHCConsumeNode at hrun
  split at hrun <;> try contradiction
  next hn =>
    dsimp only at hrun
    split at hrun <;> try contradiction
    next hq =>
      split at hrun <;> try contradiction
      next hd => exact ⟨hn, hq, hd⟩

theorem pairVecDivVHCConsumeNode_coefficient_semantics
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (nodeIndex : Nat) (k k' : UInt64)
    (nodes nodes' : Array PairVecDivVHCNode) (lin lin' : Array Nat)
    (resetH resetH' : Nat) (next : Option Nat)
    (quotient divisor : SparsePolyZp)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hcanonical : SparsePolyZp.Canonical this._p.toNat quotient)
    (hk : k.toNat < this._p.toNat)
    (hrun : pairVecDivVHCConsumeNode this nodeIndex k nodes lin resetH
      quotient divisor = .ok (k', nodes', lin', resetH', next)) :
    ∃ hn : nodeIndex < nodes.size,
      ∃ hq : nodes[nodeIndex].quotientIndex < quotient.size,
      ∃ hd : nodes[nodeIndex].divisorIndex < divisor.size,
        (k'.toNat : ZMod this._p.toNat) =
          (k.toNat : ZMod this._p.toNat) -
            Zp.toZMod this._p.toNat
                quotient[nodes[nodeIndex].quotientIndex].2 *
              Zp.toZMod this._p.toNat
                divisor[nodes[nodeIndex].divisorIndex].2 := by
  rcases pairVecDivVHCConsumeNode_indices this nodeIndex k k' nodes nodes' lin
      lin' resetH resetH' next quotient divisor hrun with ⟨hn, hq, hd⟩
  exact ⟨hn, hq, hd,
    pairVecDivVHCConsumeNode_coefficient_toZMod this nodeIndex k k' nodes
      nodes' lin lin' resetH resetH' next quotient divisor hn hq hd hcfg
      hcanonical hk hrun⟩

theorem pairVecDivVHCConsumeNode_get_ne (this : DenseUPolyZp)
    (nodeIndex : Nat) (k k' : UInt64)
    (nodes nodes' : Array PairVecDivVHCNode) (lin lin' : Array Nat)
    (resetH resetH' : Nat) (next : Option Nat)
    (quotient divisor : SparsePolyZp)
    (hrun : pairVecDivVHCConsumeNode this nodeIndex k nodes lin resetH
      quotient divisor = .ok (k', nodes', lin', resetH', next))
    (i : Nat) (hne : nodeIndex ≠ i) :
    nodes'[i]? = nodes[i]? := by
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
          exact Array.getElem?_set_ne hn hne
        next hadvance =>
          split at hrun <;> try contradiction
          next hexhausted =>
            split at hrun <;> try contradiction
            next horder =>
              simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
              rcases hrun with ⟨rfl, rfl, rfl, rfl, rfl⟩
              exact Array.getElem?_set_ne hn hne

/-- Result of consuming one complete equal-monomial `next` bucket. -/
structure PairVecDivVHCBucketResult where
  coefficient : UInt64
  nodes : Array PairVecDivVHCNode
  lin : Array Nat
  resetH : Nat
  /-- Proof-visible remainder of the finite node ownership set.  The C++
  result does not observe this field; it records exactly which `next`-chain
  nodes the safe execution consumed. -/
  unvisited : Finset Nat

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
  | none => .ok (PairVecDivVHCBucketResult.mk k nodes lin resetH unvisited)
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

/-- Execution-indexed record of the concrete coefficient pairs subtracted by
one real `next`-chain traversal.  Each `step` stores the values read before
`ConsumeNode` advances or retires that node. -/
inductive PairVecDivVHCConsumeTrace (this : DenseUPolyZp)
    (quotient divisor : SparsePolyZp) :
    Option Nat → Finset Nat → UInt64 → Array PairVecDivVHCNode →
      Array Nat → Nat → PairVecDivVHCBucketResult →
      List (UInt64 × UInt64) → Prop
  | done (unvisited : Finset Nat) (k : UInt64)
      (nodes : Array PairVecDivVHCNode) (lin : Array Nat) (resetH : Nat) :
      PairVecDivVHCConsumeTrace this quotient divisor none unvisited k nodes
        lin resetH
        (PairVecDivVHCBucketResult.mk k nodes lin resetH unvisited) []
  | step (nodeIndex : Nat) (unvisited : Finset Nat) (k k' : UInt64)
      (nodes nodes' : Array PairVecDivVHCNode) (lin lin' : Array Nat)
      (resetH resetH' : Nat) (next : Option Nat)
      (result : PairVecDivVHCBucketResult)
      (products : List (UInt64 × UInt64))
      (hmem : nodeIndex ∈ unvisited) (hn : nodeIndex < nodes.size)
      (hq : nodes[nodeIndex].quotientIndex < quotient.size)
      (hd : nodes[nodeIndex].divisorIndex < divisor.size)
      (hconsume : pairVecDivVHCConsumeNode this nodeIndex k nodes lin resetH
        quotient divisor = .ok (k', nodes', lin', resetH', next))
      (htail : PairVecDivVHCConsumeTrace this quotient divisor next
        (unvisited.erase nodeIndex) k' nodes' lin' resetH' result products) :
      PairVecDivVHCConsumeTrace this quotient divisor (some nodeIndex)
        unvisited k nodes lin resetH result
        ((quotient[nodes[nodeIndex].quotientIndex].2.val,
          divisor[nodes[nodeIndex].divisorIndex].2.val) :: products)

theorem pairVecDivVHCConsumeChain_has_trace
    (this : DenseUPolyZp) (current : Option Nat) (unvisited : Finset Nat)
    (k : UInt64) (nodes : Array PairVecDivVHCNode) (lin : Array Nat)
    (resetH : Nat) (quotient divisor : SparsePolyZp)
    (result : PairVecDivVHCBucketResult)
    (hrun : pairVecDivVHCConsumeChain this current unvisited k nodes lin
      resetH quotient divisor = .ok result) :
    ∃ products : List (UInt64 × UInt64),
      PairVecDivVHCConsumeTrace this quotient divisor current unvisited k nodes
        lin resetH result products := by
  cases current with
  | none =>
      rw [pairVecDivVHCConsumeChain] at hrun
      simp only [Except.ok.injEq] at hrun
      subst result
      exact ⟨[], PairVecDivVHCConsumeTrace.done unvisited k nodes lin resetH⟩
  | some nodeIndex =>
      rw [pairVecDivVHCConsumeChain] at hrun
      split at hrun <;> try contradiction
      next hmem =>
        cases hconsume : pairVecDivVHCConsumeNode this nodeIndex k nodes lin
            resetH quotient divisor with
        | error fault => simp [hconsume] at hrun
        | ok step =>
            rcases step with ⟨k', nodes', lin', resetH', next⟩
            rw [hconsume] at hrun
            rcases pairVecDivVHCConsumeNode_indices this nodeIndex k k' nodes
                nodes' lin lin' resetH resetH' next quotient divisor hconsume
              with ⟨hn, hq, hd⟩
            rcases pairVecDivVHCConsumeChain_has_trace this next
                (unvisited.erase nodeIndex) k' nodes' lin' resetH' quotient
                divisor result hrun with ⟨products, htrace⟩
            exact ⟨_, PairVecDivVHCConsumeTrace.step nodeIndex unvisited k k'
              nodes nodes' lin lin' resetH resetH' next result products hmem hn
              hq hd hconsume htrace⟩
termination_by unvisited.card
decreasing_by
  exact Finset.card_erase_lt_of_mem (by assumption)

def pairVecDivVHCProductsValue (p : Nat) :
    List (UInt64 × UInt64) → ZMod p
  | [] => 0
  | product :: products =>
      (product.1.toNat : ZMod p) * (product.2.toNat : ZMod p) +
        pairVecDivVHCProductsValue p products

theorem pairVecDivVHCProductsValue_append (p : Nat)
    (left right : List (UInt64 × UInt64)) :
    pairVecDivVHCProductsValue p (left ++ right) =
      pairVecDivVHCProductsValue p left +
        pairVecDivVHCProductsValue p right := by
  induction left with
  | nil => simp [pairVecDivVHCProductsValue]
  | cons product products ih =>
      simp only [List.cons_append, pairVecDivVHCProductsValue, ih]
      ring

theorem PairVecDivVHCConsumeTrace.coefficient_semantics
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (quotient divisor : SparsePolyZp) (current : Option Nat)
    (unvisited : Finset Nat) (k : UInt64)
    (nodes : Array PairVecDivVHCNode) (lin : Array Nat) (resetH : Nat)
    (result : PairVecDivVHCBucketResult)
    (products : List (UInt64 × UInt64))
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hcanonical : SparsePolyZp.Canonical this._p.toNat quotient)
    (hk : k.toNat < this._p.toNat)
    (htrace : PairVecDivVHCConsumeTrace this quotient divisor current unvisited
      k nodes lin resetH result products) :
    (result.coefficient.toNat : ZMod this._p.toNat) =
      (k.toNat : ZMod this._p.toNat) -
        pairVecDivVHCProductsValue this._p.toNat products := by
  induction htrace with
  | done => simp [pairVecDivVHCProductsValue]
  | @step nodeIndex unvisited k k' nodes nodes' lin lin' resetH resetH' next
      result products hmem hn hq hd hconsume htail ih =>
      have hk' := pairVecDivVHCConsumeNode_coefficient_reduced this nodeIndex
        k k' nodes nodes' lin lin' resetH resetH' next quotient divisor
        (by
          intro hp
          have hzero : this._p.toNat = 0 := congrArg UInt64.toNat hp
          exact (Fact.out : Nat.Prime this._p.toNat).ne_zero hzero)
        hcfg hcanonical hk hconsume
      have hnode := pairVecDivVHCConsumeNode_coefficient_toZMod this nodeIndex
        k k' nodes nodes' lin lin' resetH resetH' next quotient divisor hn hq
        hd hcfg hcanonical hk hconsume
      calc
        (result.coefficient.toNat : ZMod this._p.toNat) =
            (k'.toNat : ZMod this._p.toNat) -
              pairVecDivVHCProductsValue this._p.toNat products := ih hk'
        _ = (k.toNat : ZMod this._p.toNat) -
              pairVecDivVHCProductsValue this._p.toNat
                ((quotient[nodes[nodeIndex].quotientIndex].2.val,
                  divisor[nodes[nodeIndex].divisorIndex].2.val) :: products) := by
            rw [hnode]
            simp only [pairVecDivVHCProductsValue, Zp.toZMod]
            ring

theorem pairVecDivVHCConsumeChain_coefficient_semantics
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (current : Option Nat) (unvisited : Finset Nat) (k : UInt64)
    (nodes : Array PairVecDivVHCNode) (lin : Array Nat) (resetH : Nat)
    (quotient divisor : SparsePolyZp) (result : PairVecDivVHCBucketResult)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hcanonical : SparsePolyZp.Canonical this._p.toNat quotient)
    (hk : k.toNat < this._p.toNat)
    (hrun : pairVecDivVHCConsumeChain this current unvisited k nodes lin
      resetH quotient divisor = .ok result) :
    ∃ products : List (UInt64 × UInt64),
      PairVecDivVHCConsumeTrace this quotient divisor current unvisited k nodes
          lin resetH result products ∧
        (result.coefficient.toNat : ZMod this._p.toNat) =
          (k.toNat : ZMod this._p.toNat) -
            pairVecDivVHCProductsValue this._p.toNat products := by
  rcases pairVecDivVHCConsumeChain_has_trace this current unvisited k nodes lin
      resetH quotient divisor result hrun with ⟨products, htrace⟩
  exact ⟨products, htrace,
    htrace.coefficient_semantics this quotient divisor current unvisited k
      nodes lin resetH result products hcfg hcanonical hk⟩

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
        .ok (PairVecDivVHCBucketResult.mk k nodes lin resetH unvisited) := by
  rw [pairVecDivVHCConsumeChain]

theorem pairVecDivVHCConsumeChain_unvisited (this : DenseUPolyZp)
    (current : Option Nat) (unvisited : Finset Nat) (k : UInt64)
    (nodes : Array PairVecDivVHCNode) (lin : Array Nat) (resetH : Nat)
    (quotient divisor : SparsePolyZp) (result : PairVecDivVHCBucketResult)
    (hrun : pairVecDivVHCConsumeChain this current unvisited k nodes lin
      resetH quotient divisor = .ok result) :
    result.unvisited ⊆ unvisited ∧
      ∀ i ∈ result.unvisited, result.nodes[i]? = nodes[i]? := by
  cases current with
  | none =>
      rw [pairVecDivVHCConsumeChain] at hrun
      simp only [Except.ok.injEq] at hrun
      subst result
      exact ⟨Finset.Subset.rfl, by simp⟩
  | some nodeIndex =>
      rw [pairVecDivVHCConsumeChain] at hrun
      split at hrun <;> try contradiction
      next hmem =>
        cases hconsume : pairVecDivVHCConsumeNode this nodeIndex k nodes lin
            resetH quotient divisor with
        | error fault => simp [hconsume] at hrun
        | ok step =>
            rcases step with ⟨k', nodes', lin', resetH', next⟩
            rw [hconsume] at hrun
            have ih := pairVecDivVHCConsumeChain_unvisited this next
              (unvisited.erase nodeIndex) k' nodes' lin' resetH' quotient
              divisor result hrun
            refine ⟨?_, ?_⟩
            · intro i hi
              exact Finset.mem_of_mem_erase (ih.1 hi)
            · intro i hi
              have hiErase := ih.1 hi
              have hine : nodeIndex ≠ i := by
                exact (Finset.mem_erase.mp hiErase).1.symm
              rw [ih.2 i hi]
              exact pairVecDivVHCConsumeNode_get_ne this nodeIndex k k' nodes
                nodes' lin lin' resetH resetH' next quotient divisor hconsume
                i hine
termination_by unvisited.card
decreasing_by
  exact Finset.card_erase_lt_of_mem (by assumption)

theorem pairVecDivVHCConsumeChain_location_monotone
    (this : DenseUPolyZp) (current : Option Nat)
    (unvisited : Finset Nat) (k : UInt64)
    (nodes : Array PairVecDivVHCNode) (lin : Array Nat) (resetH : Nat)
    (quotient divisor : SparsePolyZp) (result : PairVecDivVHCBucketResult)
    (hrun : pairVecDivVHCConsumeChain this current unvisited k nodes lin
      resetH quotient divisor = .ok result) :
    lin.toList.toFinset ⊆ result.lin.toList.toFinset ∧
      resetH ≤ result.resetH := by
  cases current with
  | none =>
      rw [pairVecDivVHCConsumeChain] at hrun
      simp only [Except.ok.injEq] at hrun
      subst result
      exact ⟨Finset.Subset.rfl, Nat.le_refl _⟩
  | some nodeIndex =>
      rw [pairVecDivVHCConsumeChain] at hrun
      split at hrun <;> try contradiction
      next hmem =>
        cases hconsume : pairVecDivVHCConsumeNode this nodeIndex k nodes lin
            resetH quotient divisor with
        | error fault => simp [hconsume] at hrun
        | ok step =>
            rcases step with ⟨k', nodes', lin', resetH', next⟩
            rw [hconsume] at hrun
            have hstep := pairVecDivVHCConsumeNode_location_progress this
              nodeIndex k k' nodes nodes' lin lin' resetH resetH' next quotient
              divisor hconsume
            have htail := pairVecDivVHCConsumeChain_location_monotone this next
              (unvisited.erase nodeIndex) k' nodes' lin' resetH' quotient
              divisor result hrun
            exact ⟨Finset.Subset.trans hstep.1 htail.1,
              Nat.le_trans hstep.2.1 htail.2⟩
termination_by unvisited.card
decreasing_by
  exact Finset.card_erase_lt_of_mem (by assumption)

theorem pairVecDivVHCConsumeChain_nodes_size
    (this : DenseUPolyZp) (current : Option Nat)
    (unvisited : Finset Nat) (k : UInt64)
    (nodes : Array PairVecDivVHCNode) (lin : Array Nat) (resetH : Nat)
    (quotient divisor : SparsePolyZp) (result : PairVecDivVHCBucketResult)
    (hrun : pairVecDivVHCConsumeChain this current unvisited k nodes lin
      resetH quotient divisor = .ok result) :
    result.nodes.size = nodes.size := by
  cases current with
  | none =>
      rw [pairVecDivVHCConsumeChain] at hrun
      simp only [Except.ok.injEq] at hrun
      subst result
      rfl
  | some nodeIndex =>
      rw [pairVecDivVHCConsumeChain] at hrun
      split at hrun <;> try contradiction
      next hmem =>
        cases hconsume : pairVecDivVHCConsumeNode this nodeIndex k nodes lin
            resetH quotient divisor with
        | error fault => simp [hconsume] at hrun
        | ok step =>
            rcases step with ⟨k', nodes', lin', resetH', next⟩
            rw [hconsume] at hrun
            rw [pairVecDivVHCConsumeChain_nodes_size this next
              (unvisited.erase nodeIndex) k' nodes' lin' resetH' quotient
              divisor result hrun]
            exact pairVecDivVHCConsumeNode_nodes_size this nodeIndex k k' nodes
              nodes' lin lin' resetH resetH' next quotient divisor hconsume
termination_by unvisited.card
decreasing_by
  exact Finset.card_erase_lt_of_mem (by assumption)

theorem pairVecDivVHCConsumeChain_coefficient_reduced
    (this : DenseUPolyZp) (current : Option Nat) (unvisited : Finset Nat)
    (k : UInt64) (nodes : Array PairVecDivVHCNode) (lin : Array Nat)
    (resetH : Nat) (quotient divisor : SparsePolyZp)
    (result : PairVecDivVHCBucketResult) (hp : this._p ≠ 0)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hcanonical : SparsePolyZp.Canonical this._p.toNat quotient)
    (hk : k.toNat < this._p.toNat)
    (hrun : pairVecDivVHCConsumeChain this current unvisited k nodes lin
      resetH quotient divisor = .ok result) :
    result.coefficient.toNat < this._p.toNat := by
  cases current with
  | none =>
      rw [pairVecDivVHCConsumeChain] at hrun
      simp only [Except.ok.injEq] at hrun
      subst result
      exact hk
  | some nodeIndex =>
      rw [pairVecDivVHCConsumeChain] at hrun
      split at hrun <;> try contradiction
      next hmem =>
        cases hconsume : pairVecDivVHCConsumeNode this nodeIndex k nodes lin
            resetH quotient divisor with
        | error fault => simp [hconsume] at hrun
        | ok step =>
            rcases step with ⟨k', nodes', lin', resetH', next⟩
            rw [hconsume] at hrun
            have hk' := pairVecDivVHCConsumeNode_coefficient_reduced this
              nodeIndex k k' nodes nodes' lin lin' resetH resetH' next quotient
              divisor hp hcfg hcanonical hk hconsume
            exact pairVecDivVHCConsumeChain_coefficient_reduced this next
              (unvisited.erase nodeIndex) k' nodes' lin' resetH' quotient
              divisor result hp hcfg hcanonical hk' hrun
termination_by unvisited.card
decreasing_by
  exact Finset.card_erase_lt_of_mem (by assumption)

theorem pairVecDivVHCConsumeRootBucket_coefficient_reduced
    (this : DenseUPolyZp) (heap : Array Nat) (k : UInt64)
    (nodes : Array PairVecDivVHCNode) (lin : Array Nat) (resetH : Nat)
    (quotient divisor : SparsePolyZp) (result : PairVecDivVHCBucketResult)
    (hp : this._p ≠ 0)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hcanonical : SparsePolyZp.Canonical this._p.toNat quotient)
    (hk : k.toNat < this._p.toNat)
    (hrun : pairVecDivVHCConsumeRootBucket this heap k nodes lin resetH
      quotient divisor = .ok result) :
    result.coefficient.toNat < this._p.toNat := by
  unfold pairVecDivVHCConsumeRootBucket at hrun
  split at hrun <;> try contradiction
  next hheap =>
    exact pairVecDivVHCConsumeChain_coefficient_reduced this (some heap[0])
      (Finset.range nodes.size) k nodes lin resetH quotient divisor result hp
      hcfg hcanonical hk hrun

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

theorem pairVecDivVHCConsumeEqualDegree_coefficient_reduced
    (this : DenseUPolyZp) (degree : Nat) (heap : Array Nat) (k : UInt64)
    (nodes : Array PairVecDivVHCNode) (lin : Array Nat) (resetH : Nat)
    (quotient divisor : SparsePolyZp)
    (result : PairVecDivVHCEqualDegreeResult) (hp : this._p ≠ 0)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hcanonical : SparsePolyZp.Canonical this._p.toNat quotient)
    (hk : k.toNat < this._p.toNat)
    (hrun : pairVecDivVHCConsumeEqualDegree this degree heap k nodes lin resetH
      quotient divisor = .ok result) :
    result.coefficient.toNat < this._p.toNat := by
  induction hsize : heap.size using Nat.strong_induction_on generalizing heap k
      nodes lin resetH result with
  | h size ih =>
      rw [pairVecDivVHCConsumeEqualDegree] at hrun
      by_cases hheap : 0 < heap.size
      · simp only [hheap, ↓reduceDIte] at hrun
        cases hmono : pairVecDivVHCMono heap[0] nodes with
        | error fault => simp [hmono] at hrun
        | ok rootMono =>
            simp only [hmono] at hrun
            by_cases hequal : rootMono.deg = degree
            · simp only [hequal, ↓reduceDIte] at hrun
              cases hconsume : pairVecDivVHCConsumeRootBucket this heap k nodes
                  lin resetH quotient divisor with
              | error fault => simp [hconsume] at hrun
              | ok bucket =>
                  simp only [dif_pos trivial, hconsume] at hrun
                  cases hchecked : pairVecDivVHCExtractChecked heap
                      bucket.nodes with
                  | error fault => simp [hchecked] at hrun
                  | ok extracted =>
                      rw [hchecked] at hrun
                      have hk' :=
                        pairVecDivVHCConsumeRootBucket_coefficient_reduced this
                          heap k nodes lin resetH quotient divisor bucket hp
                          hcfg hcanonical hk hconsume
                      have hsmaller : extracted.1.size < size := by
                        rw [← hsize]
                        omega
                      exact ih extracted.1.size hsmaller extracted.1
                        bucket.coefficient bucket.nodes bucket.lin bucket.resetH
                        result hk' hrun rfl
            · simp only [hequal, ↓reduceDIte, Except.ok.injEq] at hrun
              subst result
              exact hk
      · simp only [hheap, ↓reduceDIte, Except.ok.injEq] at hrun
        subst result
        exact hk

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

inductive PairVecDivVHCFrontierSource (dividendIndex : Nat)
    (dividend : SparsePolyZp) (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode) : PairVecDivVHCFrontier → Prop
  | dividend (hindex : dividendIndex < dividend.size) :
      PairVecDivVHCFrontierSource dividendIndex dividend heap nodes
        (PairVecDivVHCFrontier.mk dividend[dividendIndex].1.deg
          dividend[dividendIndex].2.val (dividendIndex + 1))
  | heap (hheap : 0 < heap.size) (rootMono : UMonomial)
      (hmono : pairVecDivVHCMono heap[0] nodes = .ok rootMono)
      (hdominates : ∀ hindex : dividendIndex < dividend.size,
        dividend[dividendIndex].1.deg < rootMono.deg) :
      PairVecDivVHCFrontierSource dividendIndex dividend heap nodes
        (PairVecDivVHCFrontier.mk rootMono.deg 0 dividendIndex)

/-- A successful selector result is exactly one of the two concrete source
reads.  In particular, a zero coefficient is tied to an actual heap-root
monomial; it is not a specification-level default for an arbitrary degree. -/
theorem pairVecDivVHCSelectFrontier_has_source (dividendIndex : Nat)
    (dividend : SparsePolyZp) (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode) (frontier : PairVecDivVHCFrontier)
    (hrun : pairVecDivVHCSelectFrontier dividendIndex dividend heap nodes =
      .ok frontier) :
    PairVecDivVHCFrontierSource dividendIndex dividend heap nodes frontier := by
  unfold pairVecDivVHCSelectFrontier at hrun
  split at hrun
  next hdividend =>
    split at hrun
    next hheap =>
      split at hrun <;> try contradiction
      next rootMono hmono =>
        split at hrun
        next hdegree =>
          simp only [Except.ok.injEq] at hrun
          rw [← hrun]
          exact PairVecDivVHCFrontierSource.dividend hdividend
        next hdegree =>
          simp only [Except.ok.injEq] at hrun
          rw [← hrun]
          exact PairVecDivVHCFrontierSource.heap hheap rootMono hmono
            (fun _ => Nat.lt_of_not_ge hdegree)
    next hheap =>
      simp only [Except.ok.injEq] at hrun
      rw [← hrun]
      exact PairVecDivVHCFrontierSource.dividend hdividend
  next hdividend =>
    split at hrun <;> try contradiction
    next hheap =>
      split at hrun <;> try contradiction
      next rootMono hmono =>
        simp only [Except.ok.injEq] at hrun
        rw [← hrun]
        exact PairVecDivVHCFrontierSource.heap hheap rootMono hmono
          (fun hindex => by omega)

theorem pairVecDivVHCSelectFrontier_root_degree_le
    (dividendIndex : Nat) (dividend : SparsePolyZp)
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode)
    (frontier : PairVecDivVHCFrontier) (rootMono : UMonomial)
    (hheap : 0 < heap.size)
    (hroot : pairVecDivVHCMono heap[0] nodes = .ok rootMono)
    (hrun : pairVecDivVHCSelectFrontier dividendIndex dividend heap nodes =
      .ok frontier) :
    rootMono.deg ≤ frontier.degree := by
  unfold pairVecDivVHCSelectFrontier at hrun
  by_cases hdividend : dividendIndex < dividend.size
  · rw [dif_pos hdividend, dif_pos hheap, hroot] at hrun
    by_cases hdegree : rootMono.deg ≤ dividend[dividendIndex].1.deg
    · simp only [Except.bind, hdegree, if_pos] at hrun
      rw [← Except.ok.inj hrun]
      exact hdegree
    · simp only [Except.bind, hdegree, if_neg] at hrun
      rw [← Except.ok.inj hrun]
  · rw [dif_neg hdividend, dif_pos hheap, hroot] at hrun
    simp only [Except.bind] at hrun
    rw [← Except.ok.inj hrun]

theorem pairVecDivVHCSelectFrontier_heap_slot_degree_le
    (dividendIndex : Nat) (dividend : SparsePolyZp)
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode)
    (frontier : PairVecDivVHCFrontier)
    (hvalid : PairVecDivVHCHeapPointersValid heap nodes)
    (hordered : PairVecDivVHCHeapOrdered heap nodes)
    (slot : Nat) (hslot : slot < heap.size) (mono : UMonomial)
    (hmono : pairVecDivVHCMono heap[slot] nodes = .ok mono)
    (hrun : pairVecDivVHCSelectFrontier dividendIndex dividend heap nodes =
      .ok frontier) :
    mono.deg ≤ frontier.degree := by
  have hheap : 0 < heap.size := Nat.lt_of_le_of_lt (Nat.zero_le slot) hslot
  rcases hvalid 0 hheap with
    ⟨rootIndex, rootNode, rootMono, hrootHeap, hrootNode, hrootActive⟩
  have hrootIndexEq : heap[0] = rootIndex := by
    rw [Array.getElem?_eq_getElem hheap] at hrootHeap
    exact Option.some.inj hrootHeap
  have hroot : pairVecDivVHCMono heap[0] nodes = .ok rootMono := by
    apply (pairVecDivVHCMono_eq_ok_iff heap[0] nodes rootMono).mpr
    rw [hrootIndexEq]
    exact ⟨rootNode, hrootNode, hrootActive⟩
  exact Nat.le_trans
    (pairVecDivVHCHeapOrdered_slot_le_root heap nodes hvalid hordered slot
      hslot mono rootMono hmono hroot)
    (pairVecDivVHCSelectFrontier_root_degree_le dividendIndex dividend heap
      nodes frontier rootMono hheap hroot hrun)

theorem pairVecDivVHCSelectFrontier_coefficient_reduced
    (p dividendIndex : Nat) (dividend : SparsePolyZp) (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode) (frontier : PairVecDivVHCFrontier)
    (hp : 0 < p) (hcanonical : SparsePolyZp.Canonical p dividend)
    (hrun : pairVecDivVHCSelectFrontier dividendIndex dividend heap nodes =
      .ok frontier) :
    frontier.coefficient.toNat < p := by
  unfold pairVecDivVHCSelectFrontier at hrun
  split at hrun
  next hdividend =>
    have htermMem : dividend[dividendIndex] ∈ dividend.toList :=
      Array.getElem_mem_toList hdividend
    have hreduced := (hcanonical.1 dividend[dividendIndex] htermMem).2
    split at hrun
    next hheap =>
      split at hrun <;> try contradiction
      next rootMono hmono =>
        split at hrun
        next hdegree =>
          simp only [Except.ok.injEq] at hrun
          subst frontier
          exact hreduced
        next hdegree =>
          simp only [Except.ok.injEq] at hrun
          subst frontier
          simpa using hp
    next hheap =>
      simp only [Except.ok.injEq] at hrun
      subst frontier
      exact hreduced
  next hdividend =>
    split at hrun <;> try contradiction
    next hheap =>
      split at hrun <;> try contradiction
      next rootMono hmono =>
        simp only [Except.ok.injEq] at hrun
        subst frontier
        simpa using hp

/-- Every source candidate still reachable by the outer loop lies strictly
below `degreeLimit`.  This is the termination component of the full heap
representation invariant. -/
def PairVecDivVHCFrontierBelow (degreeLimit dividendIndex : Nat)
    (dividend : SparsePolyZp) (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode) : Prop :=
  (∀ (i : Nat) (term : UMonomial × Zp), dividendIndex ≤ i →
      dividend[i]? = some term → term.1.deg < degreeLimit) ∧
    (∀ (slot nodeIndex : Nat) (node : PairVecDivVHCNode)
        (mono : UMonomial),
      heap[slot]? = some nodeIndex → nodes[nodeIndex]? = some node →
      node.mono = some mono → mono.deg < degreeLimit)

def PairVecDivVHCRemainingDividendBelow (degreeLimit dividendIndex : Nat)
    (dividend : SparsePolyZp) : Prop :=
  ∀ (i : Nat) (term : UMonomial × Zp), dividendIndex ≤ i →
    dividend[i]? = some term → term.1.deg < degreeLimit

/-- Every dividend cell before the source pointer belongs to an already
processed degree, hence is at or above the next strict frontier bound. -/
def PairVecDivVHCConsumedDividendAbove (degreeLimit dividendIndex : Nat)
    (dividend : SparsePolyZp) : Prop :=
  ∀ (i : Nat) (term : UMonomial × Zp), i < dividendIndex →
    dividend[i]? = some term → degreeLimit ≤ term.1.deg

theorem pairVecDivVHCConsumedDividendAbove_zero (degreeLimit : Nat)
    (dividend : SparsePolyZp) :
    PairVecDivVHCConsumedDividendAbove degreeLimit 0 dividend := by
  intro i term hi
  omega

theorem pairVecDivVHCSelectFrontier_preserves_consumed_above
    (degreeLimit dividendIndex : Nat) (dividend : SparsePolyZp)
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode)
    (frontier : PairVecDivVHCFrontier)
    (hconsumed : PairVecDivVHCConsumedDividendAbove degreeLimit
      dividendIndex dividend)
    (hdecrease : frontier.degree < degreeLimit)
    (hrun : pairVecDivVHCSelectFrontier dividendIndex dividend heap nodes =
      .ok frontier) :
    PairVecDivVHCConsumedDividendAbove frontier.degree
      frontier.dividendIndex dividend := by
  have hsource := pairVecDivVHCSelectFrontier_has_source dividendIndex
    dividend heap nodes frontier hrun
  cases hsource with
  | dividend hindex =>
      intro i term hi hget
      change i < dividendIndex + 1 at hi
      by_cases heq : i = dividendIndex
      · subst i
        rw [Array.getElem?_eq_getElem hindex] at hget
        simp only [Option.some.injEq] at hget
        subst term
        exact Nat.le_refl _
      · have hold : i < dividendIndex := by omega
        exact Nat.le_trans (Nat.le_of_lt hdecrease)
          (hconsumed i term hold hget)
  | heap hheap rootMono hmono hdominates =>
      intro i term hi hget
      change i < dividendIndex at hi
      exact Nat.le_trans (Nat.le_of_lt hdecrease)
        (hconsumed i term hi hget)

theorem pairVecDivVHCSelectFrontier_preserves_remaining_below
    (degreeLimit dividendIndex : Nat) (dividend : SparsePolyZp)
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode)
    (frontier : PairVecDivVHCFrontier)
    (hremaining : PairVecDivVHCRemainingDividendBelow degreeLimit
      dividendIndex dividend)
    (hrun : pairVecDivVHCSelectFrontier dividendIndex dividend heap nodes =
      .ok frontier) :
    PairVecDivVHCRemainingDividendBelow degreeLimit frontier.dividendIndex
      dividend := by
  have hindex := pairVecDivVHCSelectFrontier_index dividendIndex dividend heap
    nodes frontier hrun
  have hle : dividendIndex ≤ frontier.dividendIndex := by
    rcases hindex with hsame | hnext
    · omega
    · omega
  intro i term hi hterm
  exact hremaining i term (Nat.le_trans hle hi) hterm

theorem pairVecDivVHCSelectFrontier_degree_lt (degreeLimit dividendIndex : Nat)
    (dividend : SparsePolyZp) (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode) (frontier : PairVecDivVHCFrontier)
    (hbelow : PairVecDivVHCFrontierBelow degreeLimit dividendIndex
      dividend heap nodes)
    (hrun : pairVecDivVHCSelectFrontier dividendIndex dividend heap nodes =
      .ok frontier) :
    frontier.degree < degreeLimit := by
  unfold pairVecDivVHCSelectFrontier at hrun
  split at hrun
  next hdividend =>
    have hdividendGet : dividend[dividendIndex]? =
        some dividend[dividendIndex] := Array.getElem?_eq_getElem hdividend
    have hdividendDegree := hbelow.1 dividendIndex dividend[dividendIndex]
      (Nat.le_refl _) hdividendGet
    split at hrun
    next hheap =>
      have hheapGet : heap[0]? = some heap[0] :=
        Array.getElem?_eq_getElem hheap
      split at hrun <;> try contradiction
      next rootMono hmono =>
        have hnode : ∃ node, nodes[heap[0]]? = some node ∧
            node.mono = some rootMono := by
          unfold pairVecDivVHCMono at hmono
          split at hmono <;> try contradiction
          next hn =>
            split at hmono <;> try contradiction
            next mono hnodeMono =>
              simp only [Except.ok.injEq] at hmono
              subst mono
              exact ⟨nodes[heap[0]], Array.getElem?_eq_getElem hn, hnodeMono⟩
        rcases hnode with ⟨node, hnodeGet, hnodeMono⟩
        have hrootDegree := hbelow.2 0 heap[0] node rootMono hheapGet
          hnodeGet hnodeMono
        split at hrun <;>
          simp only [Except.ok.injEq] at hrun <;>
          subst frontier
        · exact hdividendDegree
        · exact hrootDegree
    next hheap =>
      simp only [Except.ok.injEq] at hrun
      subst frontier
      exact hdividendDegree
  next hdividend =>
    split at hrun <;> try contradiction
    next hheap =>
      have hheapGet : heap[0]? = some heap[0] :=
        Array.getElem?_eq_getElem hheap
      split at hrun <;> try contradiction
      next rootMono hmono =>
        unfold pairVecDivVHCMono at hmono
        split at hmono <;> try contradiction
        next hn =>
          split at hmono <;> try contradiction
          next mono hnodeMono =>
            simp only [Except.ok.injEq] at hmono
            subst mono
            simp only [Except.ok.injEq] at hrun
            subst frontier
            exact hbelow.2 0 heap[0] nodes[heap[0]] rootMono hheapGet
              (Array.getElem?_eq_getElem hn) hnodeMono

theorem pairVecDivVHCInitialFrontierBelow (dividend : SparsePolyZp)
    (nodes : Array PairVecDivVHCNode) (hnonempty : 0 < dividend.size)
    (hdegree : ∀ (i : Nat) (term : UMonomial × Zp),
      dividend[i]? = some term → term.1.deg ≤ dividend[0].1.deg) :
    PairVecDivVHCFrontierBelow (dividend[0].1.deg + 1) 0 dividend #[]
      nodes := by
  refine ⟨?_, ?_⟩
  · intro i term hi hget
    have := hdegree i term hget
    omega
  · intro slot nodeIndex node mono hheap
    simp at hheap

theorem canonical_degree_le_head (p : Nat) (dividend : SparsePolyZp)
    (hcanonical : SparsePolyZp.Canonical p dividend)
    (hnonempty : 0 < dividend.size) :
    ∀ (i : Nat) (term : UMonomial × Zp), dividend[i]? = some term →
      term.1.deg ≤ dividend[0].1.deg := by
  intro i term hget
  have hi : i < dividend.size := by
    by_contra hnot
    have hnone : dividend[i]? = none :=
      Array.getElem?_eq_none (by omega)
    rw [hnone] at hget
    contradiction
  rw [Array.getElem?_eq_getElem hi] at hget
  simp only [Option.some.injEq] at hget
  subst term
  by_cases hizero : i = 0
  · subst i
    exact Nat.le_refl _
  · have hpairwise := List.isChain_iff_pairwise.mp hcanonical.2.1
    have horder := List.pairwise_iff_getElem.mp hpairwise 0 i
      (by simpa using hnonempty) (by simpa using hi) (by omega)
    rw [Array.getElem_toList hnonempty, Array.getElem_toList hi] at horder
    exact Nat.le_of_lt horder

theorem canonical_degree_lt_of_index_lt (p : Nat)
    (dividend : SparsePolyZp)
    (hcanonical : SparsePolyZp.Canonical p dividend)
    (i j : Nat) (hi : i < dividend.size) (hj : j < dividend.size)
    (hij : i < j) :
    dividend[j].1.deg < dividend[i].1.deg := by
  have hpairwise := List.isChain_iff_pairwise.mp hcanonical.2.1
  have horder := List.pairwise_iff_getElem.mp hpairwise i j
    (by simpa using hi) (by simpa using hj) hij
  rw [Array.getElem_toList hi, Array.getElem_toList hj] at horder
  exact horder

theorem PairVecDivVHCNodeDenotes.later_product_degree_lt
    (p : Nat) (quotient divisor : SparsePolyZp)
    (node : PairVecDivVHCNode) (degree q : Nat)
    (quotientTerm divisorTerm : UMonomial × Zp)
    (hcanonical : SparsePolyZp.Canonical p quotient)
    (hdenotes : PairVecDivVHCNodeDenotes quotient divisor node)
    (hmono : node.mono = some ⟨degree⟩)
    (hlater : node.quotientIndex < q)
    (hquotient : quotient[q]? = some quotientTerm)
    (hdivisor : divisor[node.divisorIndex]? = some divisorTerm) :
    quotientTerm.1.deg + divisorTerm.1.deg < degree := by
  rcases hdenotes with
    ⟨currentQuotient, currentDivisor, hcurrentQuotient, hcurrentDivisor,
      hstoredMono⟩
  have hcurrentIndex : node.quotientIndex < quotient.size := by
    by_contra hnot
    rw [Array.getElem?_eq_none (by omega)] at hcurrentQuotient
    contradiction
  have hq : q < quotient.size := by
    by_contra hnot
    rw [Array.getElem?_eq_none (by omega)] at hquotient
    contradiction
  have hdegreeLt := canonical_degree_lt_of_index_lt p quotient hcanonical
    node.quotientIndex q hcurrentIndex hq hlater
  rw [Array.getElem?_eq_getElem hcurrentIndex] at hcurrentQuotient
  rw [Array.getElem?_eq_getElem hq] at hquotient
  simp only [Option.some.injEq] at hcurrentQuotient hquotient
  have hcurrentDegree : currentQuotient.1.deg + currentDivisor.1.deg =
      degree := by
    rw [hmono] at hstoredMono
    exact congrArg UMonomial.deg (Option.some.inj hstoredMono).symm
  rw [hdivisor] at hcurrentDivisor
  simp only [Option.some.injEq] at hcurrentDivisor
  subst currentDivisor
  subst currentQuotient
  subst quotientTerm
  omega

theorem pairVecDivVHCProductAtFrontier_eq_cursor
    (p degreeLimit frontierDegree : Nat)
    (nodes : Array PairVecDivVHCNode)
    (quotient divisor : SparsePolyZp) (i q : Nat)
    (node : PairVecDivVHCNode) (quotientTerm divisorTerm : UMonomial × Zp)
    (hcanonical : SparsePolyZp.Canonical p quotient)
    (hprefix : PairVecDivVHCCursorPrefixAbove degreeLimit nodes quotient divisor)
    (hfrontier : frontierDegree < degreeLimit)
    (hget : nodes[i]? = some node)
    (hdenotes : PairVecDivVHCNodeDenotes quotient divisor node)
    (hmono : node.mono = some ⟨frontierDegree⟩)
    (hquotient : quotient[q]? = some quotientTerm)
    (hdivisor : divisor[node.divisorIndex]? = some divisorTerm)
    (hdegree : quotientTerm.1.deg + divisorTerm.1.deg = frontierDegree) :
    q = node.quotientIndex := by
  rcases Nat.lt_trichotomy q node.quotientIndex with hearlier | heq | hlater
  · have habove := hprefix.earlier_product_degree_gt degreeLimit
      frontierDegree nodes quotient divisor hfrontier i node hget q
      quotientTerm divisorTerm hearlier hquotient hdivisor
    omega
  · exact heq
  · have hbelow := hdenotes.later_product_degree_lt p quotient divisor node
      frontierDegree q quotientTerm divisorTerm hcanonical hmono hlater
      hquotient hdivisor
    omega

theorem pairVecDivVHCProductAtFrontier_eq_cursor_of_mono_le
    (p degreeLimit frontierDegree : Nat)
    (nodes : Array PairVecDivVHCNode)
    (quotient divisor : SparsePolyZp) (i q : Nat)
    (node : PairVecDivVHCNode) (mono : UMonomial)
    (quotientTerm divisorTerm : UMonomial × Zp)
    (hcanonical : SparsePolyZp.Canonical p quotient)
    (hprefix : PairVecDivVHCCursorPrefixAbove degreeLimit nodes quotient divisor)
    (hfrontier : frontierDegree < degreeLimit)
    (hget : nodes[i]? = some node)
    (hdenotes : PairVecDivVHCNodeDenotes quotient divisor node)
    (hmono : node.mono = some mono)
    (hmonoLe : mono.deg ≤ frontierDegree)
    (hquotient : quotient[q]? = some quotientTerm)
    (hdivisor : divisor[node.divisorIndex]? = some divisorTerm)
    (hdegree : quotientTerm.1.deg + divisorTerm.1.deg = frontierDegree) :
    q = node.quotientIndex := by
  rcases Nat.lt_trichotomy q node.quotientIndex with hearlier | heq | hlater
  · have habove := hprefix.earlier_product_degree_gt degreeLimit
      frontierDegree nodes quotient divisor hfrontier i node hget q
      quotientTerm divisorTerm hearlier hquotient hdivisor
    omega
  · exact heq
  · have hbelow := hdenotes.later_product_degree_lt p quotient divisor node
      mono.deg q quotientTerm divisorTerm hcanonical hmono hlater hquotient
      hdivisor
    omega

theorem pairVecDivVHCNode_advanced_degree_lt (p : Nat)
    (quotient divisor : SparsePolyZp)
    (node : PairVecDivVHCNode) (currentMono : UMonomial)
    (hcanonical : SparsePolyZp.Canonical p quotient)
    (hdenotes : PairVecDivVHCNodeDenotes quotient divisor node)
    (hactive : node.mono = some currentMono)
    (hd : node.divisorIndex < divisor.size)
    (hadvance : node.quotientIndex + 1 < quotient.size) :
    quotient[node.quotientIndex + 1].1.deg +
        divisor[node.divisorIndex].1.deg < currentMono.deg := by
  rcases hdenotes with ⟨quotientTerm, divisorTerm, hquotient,
    hdivisor, hmono⟩
  have hq : node.quotientIndex < quotient.size := by
    by_contra hnot
    rw [Array.getElem?_eq_none (by omega)] at hquotient
    contradiction
  rw [Array.getElem?_eq_getElem hq] at hquotient
  rw [Array.getElem?_eq_getElem hd] at hdivisor
  simp only [Option.some.injEq] at hquotient hdivisor
  subst quotientTerm
  subst divisorTerm
  rw [hactive] at hmono
  simp only [Option.some.injEq] at hmono
  have hcurrentDegree : currentMono.deg =
      quotient[node.quotientIndex].1.deg +
        divisor[node.divisorIndex].1.deg := by
    simpa using congrArg UMonomial.deg hmono
  have hdegree := canonical_degree_lt_of_index_lt p quotient hcanonical
    node.quotientIndex (node.quotientIndex + 1) hq hadvance (by omega)
  omega

theorem canonical_remaining_below_advanced (p : Nat)
    (dividend : SparsePolyZp)
    (hcanonical : SparsePolyZp.Canonical p dividend)
    (dividendIndex : Nat) (hindex : dividendIndex < dividend.size) :
    ∀ (i : Nat) (term : UMonomial × Zp), dividendIndex + 1 ≤ i →
      dividend[i]? = some term →
      term.1.deg < dividend[dividendIndex].1.deg := by
  intro i term hafter hget
  have hi : i < dividend.size := by
    by_contra hnot
    have hnone : dividend[i]? = none := Array.getElem?_eq_none (by omega)
    rw [hnone] at hget
    contradiction
  rw [Array.getElem?_eq_getElem hi] at hget
  simp only [Option.some.injEq] at hget
  subst term
  exact canonical_degree_lt_of_index_lt p dividend hcanonical
    dividendIndex i hindex hi (by omega)

theorem canonical_remaining_below_larger_frontier (p : Nat)
    (dividend : SparsePolyZp)
    (hcanonical : SparsePolyZp.Canonical p dividend)
    (dividendIndex frontierDegree : Nat)
    (hindex : dividendIndex < dividend.size)
    (hfrontier : dividend[dividendIndex].1.deg < frontierDegree) :
    ∀ (i : Nat) (term : UMonomial × Zp), dividendIndex ≤ i →
      dividend[i]? = some term → term.1.deg < frontierDegree := by
  intro i term hafter hget
  by_cases heq : i = dividendIndex
  · subst i
    rw [Array.getElem?_eq_getElem hindex] at hget
    simp only [Option.some.injEq] at hget
    subst term
    exact hfrontier
  · have hstrict := canonical_remaining_below_advanced p dividend hcanonical
      dividendIndex hindex i term (by omega) hget
    exact lt_trans hstrict hfrontier

theorem pairVecDivVHCSelectFrontier_remaining_dividend_below (p : Nat)
    (dividend : SparsePolyZp) (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode)
    (hcanonical : SparsePolyZp.Canonical p dividend)
    (dividendIndex : Nat) (frontier : PairVecDivVHCFrontier)
    (hselect : pairVecDivVHCSelectFrontier dividendIndex dividend heap nodes =
      .ok frontier) :
    ∀ (i : Nat) (term : UMonomial × Zp), frontier.dividendIndex ≤ i →
      dividend[i]? = some term → term.1.deg < frontier.degree := by
  unfold pairVecDivVHCSelectFrontier at hselect
  split at hselect
  next hdividend =>
    split at hselect
    next hheap =>
      split at hselect <;> try contradiction
      next rootMono hroot =>
        split at hselect
        next htakeDividend =>
          simp only [Except.ok.injEq] at hselect
          subst frontier
          exact canonical_remaining_below_advanced p dividend hcanonical
            dividendIndex hdividend
        next htakeHeap =>
          simp only [Except.ok.injEq] at hselect
          subst frontier
          exact canonical_remaining_below_larger_frontier p dividend hcanonical
            dividendIndex rootMono.deg hdividend (by omega)
    next hheap =>
      simp only [Except.ok.injEq] at hselect
      subst frontier
      exact canonical_remaining_below_advanced p dividend hcanonical
        dividendIndex hdividend
  next hdividend =>
    split at hselect <;> try contradiction
    next hheap =>
      split at hselect <;> try contradiction
      next rootMono hroot =>
        simp only [Except.ok.injEq] at hselect
        subst frontier
        intro i term hafter hget
        have hstart : dividend.size ≤ dividendIndex := by omega
        have hsize : dividend.size ≤ i := le_trans hstart hafter
        have hnone : dividend[i]? = none := Array.getElem?_eq_none hsize
        rw [hnone] at hget
        contradiction

/-- Strong heap-side termination invariant: every initialized node in the
allocated block is below the current frontier, independently of its temporary
location in `heap`, `next`, `lin`, or the retired suffix. -/
def PairVecDivVHCAllActiveNodesBelow (degreeLimit : Nat)
    (nodes : Array PairVecDivVHCNode) : Prop :=
  ∀ (i : Nat) (node : PairVecDivVHCNode) (mono : UMonomial),
    nodes[i]? = some node → node.mono = some mono → mono.deg < degreeLimit

/-- The source `reset_h` prefix consists exactly of nodes whose incremented
quotient pointer is one-past the old quotient and whose product fields are
therefore not observable until activation after the next append. -/
def PairVecDivVHCResetReady (resetH quotientSize : Nat)
    (nodes : Array PairVecDivVHCNode) : Prop :=
  resetH ≤ nodes.size ∧
    ∀ (i : Nat), i < resetH →
      ∃ node, nodes[i]? = some node ∧
        node.quotientIndex = quotientSize ∧
        node.divisorIndex = i + 1 ∧ node.mono = none

theorem PairVecDivVHCResetReady.shrink_set
    (resetH resetH' quotientSize nodeIndex : Nat)
    (nodes : Array PairVecDivVHCNode) (updated : PairVecDivVHCNode)
    (hn : nodeIndex < nodes.size)
    (hready : PairVecDivVHCResetReady resetH quotientSize nodes)
    (hshrink : resetH' ≤ resetH) (houtside : resetH' ≤ nodeIndex) :
    PairVecDivVHCResetReady resetH' quotientSize
      (nodes.set nodeIndex updated) := by
  refine ⟨by simpa using Nat.le_trans hshrink hready.1, ?_⟩
  intro i hi
  rcases hready.2 i (Nat.lt_of_lt_of_le hi hshrink) with
    ⟨node, hget, hqIndex, hdIndex, hmono⟩
  refine ⟨node, ?_, hqIndex, hdIndex, hmono⟩
  rw [Array.getElem?_set_ne hn (by omega)]
  exact hget

/-- Allocation identity of every source node: array slot `i` permanently
points at divisor tail cell `i + 1`.  No source transition mutates `v2_ptr`. -/
def PairVecDivVHCNodeDivisorIndicesFixed
    (nodes : Array PairVecDivVHCNode) : Prop :=
  ∀ (i : Nat) (node : PairVecDivVHCNode),
    nodes[i]? = some node → node.divisorIndex = i + 1

theorem PairVecDivVHCNodeDivisorIndicesFixed.node_for_tail
    (nodes : Array PairVecDivVHCNode) (divisorSize d : Nat)
    (hsize : nodes.size = divisorSize - 1)
    (hfixed : PairVecDivVHCNodeDivisorIndicesFixed nodes)
    (hdpos : 0 < d) (hd : d < divisorSize) :
    ∃ node, nodes[d - 1]? = some node ∧ node.divisorIndex = d := by
  have hi : d - 1 < nodes.size := by
    rw [hsize]
    omega
  let node := nodes[d - 1]
  have hget : nodes[d - 1]? = some node := Array.getElem?_eq_getElem hi
  refine ⟨node, hget, ?_⟩
  rw [hfixed (d - 1) node hget]
  omega

theorem PairVecDivVHCNodeDivisorIndicesFixed.index_eq_of_divisorIndex
    (nodes : Array PairVecDivVHCNode) (i d : Nat)
    (node : PairVecDivVHCNode)
    (hfixed : PairVecDivVHCNodeDivisorIndicesFixed nodes)
    (hget : nodes[i]? = some node) (hd : node.divisorIndex = d) :
    i = d - 1 := by
  have hidentity := hfixed i node hget
  omega

theorem PairVecDivVHCNodeDivisorIndicesFixed.unique_for_tail
    (nodes : Array PairVecDivVHCNode) (left right d : Nat)
    (leftNode rightNode : PairVecDivVHCNode)
    (hfixed : PairVecDivVHCNodeDivisorIndicesFixed nodes)
    (hleft : nodes[left]? = some leftNode)
    (hright : nodes[right]? = some rightNode)
    (hleftD : leftNode.divisorIndex = d)
    (hrightD : rightNode.divisorIndex = d) :
    left = right := by
  rw [hfixed left leftNode hleft] at hleftD
  rw [hfixed right rightNode hright] at hrightD
  omega

theorem pairVecDivVHCResetRow_product_degree_gt_frontier
    (degreeLimit frontierDegree resetH d q : Nat)
    (nodes : Array PairVecDivVHCNode)
    (quotient divisor : SparsePolyZp)
    (quotientTerm divisorTerm : UMonomial × Zp)
    (hsize : nodes.size = divisor.size - 1)
    (hfixed : PairVecDivVHCNodeDivisorIndicesFixed nodes)
    (hready : PairVecDivVHCResetReady resetH quotient.size nodes)
    (hprefix : PairVecDivVHCCursorPrefixAbove degreeLimit nodes quotient divisor)
    (hfrontier : frontierDegree < degreeLimit)
    (hdpos : 0 < d) (hd : d < divisor.size)
    (hreset : d - 1 < resetH)
    (hquotient : quotient[q]? = some quotientTerm)
    (hdivisor : divisor[d]? = some divisorTerm) :
    frontierDegree < quotientTerm.1.deg + divisorTerm.1.deg := by
  rcases hfixed.node_for_tail nodes divisor.size d hsize hdpos hd with
    ⟨node, hnode, hnodeD⟩
  rcases hready.2 (d - 1) hreset with
    ⟨readyNode, hreadyNode, hcursor, hreadyD, hmono⟩
  rw [hnode] at hreadyNode
  simp only [Option.some.injEq] at hreadyNode
  subst readyNode
  have hq : q < quotient.size := by
    by_contra hnot
    rw [Array.getElem?_eq_none (by omega)] at hquotient
    contradiction
  have hearlier : q < node.quotientIndex := by omega
  apply hprefix.earlier_product_degree_gt degreeLimit frontierDegree nodes
    quotient divisor hfrontier (d - 1) node hnode q quotientTerm divisorTerm
    hearlier hquotient
  simpa [hnodeD] using hdivisor

theorem PairVecDivVHCNodeDivisorIndicesFixed.set
    (nodes : Array PairVecDivVHCNode) (nodeIndex : Nat)
    (updated : PairVecDivVHCNode) (hn : nodeIndex < nodes.size)
    (hfixed : PairVecDivVHCNodeDivisorIndicesFixed nodes)
    (hsame : updated.divisorIndex = nodes[nodeIndex].divisorIndex) :
    PairVecDivVHCNodeDivisorIndicesFixed (nodes.set nodeIndex updated) := by
  intro i node hget
  by_cases heq : nodeIndex = i
  · subst i
    rw [Array.getElem?_set_self hn] at hget
    simp only [Option.some.injEq] at hget
    subst node
    rw [hsame]
    exact hfixed nodeIndex nodes[nodeIndex]
      (Array.getElem?_eq_getElem hn)
  · rw [Array.getElem?_set_ne hn heq] at hget
    exact hfixed i node hget

theorem pairVecDivVHCInit_resetReady (divisor : SparsePolyZp) :
    PairVecDivVHCResetReady (divisor.size - 1) 0
      (pairVecDivVHCInit divisor) := by
  refine ⟨?_, ?_⟩
  · rw [pairVecDivVHCInit_size]
  · intro i hi
    refine ⟨pairVecDivVHCInitialNode i,
      pairVecDivVHCInit_get divisor i hi, ?_⟩
    simp [pairVecDivVHCInitialNode]

theorem pairVecDivVHCInit_divisorIndicesFixed (divisor : SparsePolyZp) :
    PairVecDivVHCNodeDivisorIndicesFixed (pairVecDivVHCInit divisor) := by
  intro i node hget
  have hi : i < divisor.size - 1 := by
    by_contra hnot
    rw [Array.getElem?_eq_none (by
      rw [pairVecDivVHCInit_size]
      omega)] at hget
    contradiction
  rw [pairVecDivVHCInit_get divisor i hi] at hget
  simp only [Option.some.injEq] at hget
  subst node
  rfl

/-- Precisely the initialized `next` chain owned by the inner source loop.
Unlike a blanket node-block premise, this permits the genuine uninitialized
`reset_h` suffix.  Erasing the current node is both ownership transfer and the
well-founded recursion measure. -/
def PairVecDivVHCChainValid (current : Option Nat) (unvisited : Finset Nat)
    (nodes : Array PairVecDivVHCNode) : Prop :=
  match current with
  | none => True
  | some nodeIndex =>
      if hmem : nodeIndex ∈ unvisited then
        ∃ node mono, nodes[nodeIndex]? = some node ∧
          node.mono = some mono ∧
          PairVecDivVHCChainValid node.next (unvisited.erase nodeIndex) nodes
      else
        False
termination_by unvisited.card
decreasing_by
  exact Finset.card_erase_lt_of_mem hmem

def PairVecDivVHCChainAtDegree (current : Option Nat)
    (unvisited : Finset Nat) (nodes : Array PairVecDivVHCNode)
    (degree : Nat) : Prop :=
  match current with
  | none => True
  | some nodeIndex =>
      if hmem : nodeIndex ∈ unvisited then
        ∃ node mono, nodes[nodeIndex]? = some node ∧
          node.mono = some mono ∧ mono.deg = degree ∧
          PairVecDivVHCChainAtDegree node.next (unvisited.erase nodeIndex)
            nodes degree
      else
        False
termination_by unvisited.card
decreasing_by
  exact Finset.card_erase_lt_of_mem hmem

theorem pairVecDivVHCChainAtDegree_congr_on
    (current : Option Nat) (unvisited : Finset Nat)
    (nodes nodes' : Array PairVecDivVHCNode) (degree : Nat)
    (hdegree : PairVecDivVHCChainAtDegree current unvisited nodes degree)
    (hsame : ∀ i ∈ unvisited, nodes'[i]? = nodes[i]?) :
    PairVecDivVHCChainAtDegree current unvisited nodes' degree := by
  cases current with
  | none => simp [PairVecDivVHCChainAtDegree]
  | some nodeIndex =>
      rw [PairVecDivVHCChainAtDegree] at hdegree ⊢
      split at hdegree <;> try contradiction
      next hmem =>
        simp only [hmem, ↓reduceDIte]
        rcases hdegree with ⟨node, mono, hget, hmono, hmonoDegree, htail⟩
        refine ⟨node, mono, ?_, hmono, hmonoDegree, ?_⟩
        · rw [hsame nodeIndex hmem]
          exact hget
        · exact pairVecDivVHCChainAtDegree_congr_on node.next
            (unvisited.erase nodeIndex) nodes nodes' degree htail (by
              intro i hi
              exact hsame i (Finset.mem_of_mem_erase hi))
termination_by unvisited.card
decreasing_by
  exact Finset.card_erase_lt_of_mem (by assumption)

theorem pairVecDivVHCChainAtDegree_mono
    (current : Option Nat) (small large : Finset Nat)
    (nodes : Array PairVecDivVHCNode) (degree : Nat)
    (hdegree : PairVecDivVHCChainAtDegree current small nodes degree)
    (hsubset : small ⊆ large) :
    PairVecDivVHCChainAtDegree current large nodes degree := by
  cases current with
  | none => simp [PairVecDivVHCChainAtDegree]
  | some nodeIndex =>
      rw [PairVecDivVHCChainAtDegree] at hdegree ⊢
      split at hdegree <;> try contradiction
      next hmem =>
        have hmemLarge := hsubset hmem
        simp only [hmemLarge, ↓reduceDIte]
        rcases hdegree with
          ⟨node, mono, hget, hmono, hmonoDegree, htail⟩
        exact ⟨node, mono, hget, hmono, hmonoDegree,
          pairVecDivVHCChainAtDegree_mono node.next (small.erase nodeIndex)
            (large.erase nodeIndex) nodes degree htail (by
              intro i hi
              exact Finset.mem_erase.mpr ⟨(Finset.mem_erase.mp hi).1,
                hsubset (Finset.mem_of_mem_erase hi)⟩)⟩
termination_by small.card
decreasing_by
  exact Finset.card_erase_lt_of_mem (by assumption)

theorem pairVecDivVHCChainValid_atDegree_of_nextValid
    (nodeIndex : Nat) (unvisited : Finset Nat)
    (nodes : Array PairVecDivVHCNode) (degree : Nat)
    (hvalid : PairVecDivVHCChainValid (some nodeIndex) unvisited nodes)
    (hnextValid : PairVecDivVHCNextValid nodes)
    (hhead : ∃ node mono, nodes[nodeIndex]? = some node ∧
      node.mono = some mono ∧ mono.deg = degree) :
    PairVecDivVHCChainAtDegree (some nodeIndex) unvisited nodes degree := by
  rw [PairVecDivVHCChainValid] at hvalid
  rw [PairVecDivVHCChainAtDegree]
  split at hvalid <;> try contradiction
  next hmem =>
    simp only [hmem, ↓reduceDIte]
    rcases hvalid with ⟨node, mono, hget, hmono, htail⟩
    rcases hhead with ⟨headNode, headMono, hheadGet, hheadMono,
      hheadDegree⟩
    rw [hget] at hheadGet
    simp only [Option.some.injEq] at hheadGet
    subst headNode
    rw [hmono] at hheadMono
    simp only [Option.some.injEq] at hheadMono
    subst headMono
    refine ⟨node, mono, hget, hmono, hheadDegree, ?_⟩
    cases hnext : node.next with
    | none => simp [PairVecDivVHCChainAtDegree]
    | some nextIndex =>
        rcases hnextValid nodeIndex node nextIndex hget hnext with
          ⟨nextNode, hnextGet, hmonoEq⟩
        apply pairVecDivVHCChainValid_atDegree_of_nextValid nextIndex
          (unvisited.erase nodeIndex) nodes degree (by simpa [hnext] using htail)
          hnextValid
        refine ⟨nextNode, mono, hnextGet, ?_, hheadDegree⟩
        rw [hmono] at hmonoEq
        exact hmonoEq.symm
termination_by unvisited.card
decreasing_by
  exact Finset.card_erase_lt_of_mem (by assumption)

/-- Exact ownership variant of `PairVecDivVHCChainValid`: after following the
whole `next` chain, no owner node remains.  This makes disjoint heap buckets
stable under destructive consumption. -/
def PairVecDivVHCChainOwns (current : Option Nat) (owner : Finset Nat)
    (nodes : Array PairVecDivVHCNode) : Prop :=
  match current with
  | none => owner = ∅
  | some nodeIndex =>
      if hmem : nodeIndex ∈ owner then
        ∃ node mono, nodes[nodeIndex]? = some node ∧
          node.mono = some mono ∧
          PairVecDivVHCChainOwns node.next (owner.erase nodeIndex) nodes
      else
        False
termination_by owner.card
decreasing_by
  exact Finset.card_erase_lt_of_mem hmem

/-- Exact chain ownership is functional: for fixed concrete node storage and
head, following `next` determines the owner set uniquely.  This lets later
insertion proofs reconcile a freshly constructed ownership witness with the
owner sets used by total node coverage. -/
theorem pairVecDivVHCChainOwns_unique
    (current : Option Nat) (left right : Finset Nat)
    (nodes : Array PairVecDivVHCNode)
    (hleft : PairVecDivVHCChainOwns current left nodes)
    (hright : PairVecDivVHCChainOwns current right nodes) :
    left = right := by
  cases current with
  | none =>
      simp only [PairVecDivVHCChainOwns] at hleft hright
      rw [hleft, hright]
  | some nodeIndex =>
      rw [PairVecDivVHCChainOwns] at hleft hright
      split at hleft <;> try contradiction
      next hleftMem =>
        split at hright <;> try contradiction
        next hrightMem =>
          rcases hleft with
            ⟨leftNode, leftMono, hleftGet, hleftMono, hleftTail⟩
          rcases hright with
            ⟨rightNode, rightMono, hrightGet, hrightMono, hrightTail⟩
          rw [hleftGet] at hrightGet
          simp only [Option.some.injEq] at hrightGet
          subst rightNode
          have htailEq := pairVecDivVHCChainOwns_unique leftNode.next
            (left.erase nodeIndex) (right.erase nodeIndex) nodes hleftTail
            hrightTail
          rw [← Finset.insert_erase hleftMem,
            ← Finset.insert_erase hrightMem, htailEq]
termination_by left.card
decreasing_by
  exact Finset.card_erase_lt_of_mem (by assumption)

theorem pairVecDivVHCChainValid_set_of_not_mem
    (current : Option Nat) (unvisited : Finset Nat)
    (nodes : Array PairVecDivVHCNode) (nodeIndex : Nat)
    (updated : PairVecDivVHCNode) (hn : nodeIndex < nodes.size)
    (hnotmem : nodeIndex ∉ unvisited)
    (hvalid : PairVecDivVHCChainValid current unvisited nodes) :
    PairVecDivVHCChainValid current unvisited
      (nodes.set nodeIndex updated) := by
  cases current with
  | none => simp [PairVecDivVHCChainValid]
  | some currentIndex =>
      rw [PairVecDivVHCChainValid] at hvalid ⊢
      split at hvalid <;> try contradiction
      next hmem =>
        simp only [hmem, ↓reduceDIte]
        rcases hvalid with ⟨node, mono, hget, hmono, htail⟩
        have hne : nodeIndex ≠ currentIndex := by
          intro heq
          subst currentIndex
          exact hnotmem hmem
        refine ⟨node, mono, ?_, hmono, ?_⟩
        · rwa [Array.getElem?_set_ne hn hne]
        · exact pairVecDivVHCChainValid_set_of_not_mem node.next
            (unvisited.erase currentIndex) nodes nodeIndex updated hn
            (by simp [hnotmem]) htail
termination_by unvisited.card
decreasing_by
  exact Finset.card_erase_lt_of_mem (by assumption)

theorem pairVecDivVHCChainValid_congr_on
    (current : Option Nat) (unvisited : Finset Nat)
    (nodes nodes' : Array PairVecDivVHCNode)
    (hvalid : PairVecDivVHCChainValid current unvisited nodes)
    (hsame : ∀ i ∈ unvisited, nodes'[i]? = nodes[i]?) :
    PairVecDivVHCChainValid current unvisited nodes' := by
  cases current with
  | none => simp [PairVecDivVHCChainValid]
  | some nodeIndex =>
      rw [PairVecDivVHCChainValid] at hvalid ⊢
      split at hvalid <;> try contradiction
      next hmem =>
        simp only [hmem, ↓reduceDIte]
        rcases hvalid with ⟨node, mono, hget, hmono, htail⟩
        refine ⟨node, mono, ?_, hmono, ?_⟩
        · rw [hsame nodeIndex hmem]
          exact hget
        · exact pairVecDivVHCChainValid_congr_on node.next
            (unvisited.erase nodeIndex) nodes nodes' htail (by
              intro i hi
              exact hsame i (Finset.mem_of_mem_erase hi))
termination_by unvisited.card
decreasing_by
  exact Finset.card_erase_lt_of_mem (by assumption)

theorem pairVecDivVHCChainValid_mono
    (current : Option Nat) (small large : Finset Nat)
    (nodes : Array PairVecDivVHCNode)
    (hvalid : PairVecDivVHCChainValid current small nodes)
    (hsubset : small ⊆ large) :
    PairVecDivVHCChainValid current large nodes := by
  cases current with
  | none => simp [PairVecDivVHCChainValid]
  | some nodeIndex =>
      rw [PairVecDivVHCChainValid] at hvalid ⊢
      split at hvalid <;> try contradiction
      next hmem =>
        have hmemLarge := hsubset hmem
        simp only [hmemLarge, ↓reduceDIte]
        rcases hvalid with ⟨node, mono, hget, hmono, htail⟩
        exact ⟨node, mono, hget, hmono,
          pairVecDivVHCChainValid_mono node.next (small.erase nodeIndex)
            (large.erase nodeIndex) nodes htail (by
              intro i hi
              exact Finset.mem_erase.mpr ⟨(Finset.mem_erase.mp hi).1,
                hsubset (Finset.mem_of_mem_erase hi)⟩)⟩
termination_by small.card
decreasing_by
  exact Finset.card_erase_lt_of_mem (by assumption)

theorem pairVecDivVHCChainOwns_congr_on
    (current : Option Nat) (owner : Finset Nat)
    (nodes nodes' : Array PairVecDivVHCNode)
    (howns : PairVecDivVHCChainOwns current owner nodes)
    (hsame : ∀ i ∈ owner, nodes'[i]? = nodes[i]?) :
    PairVecDivVHCChainOwns current owner nodes' := by
  cases current with
  | none => simpa [PairVecDivVHCChainOwns] using howns
  | some nodeIndex =>
      rw [PairVecDivVHCChainOwns] at howns ⊢
      split at howns <;> try contradiction
      next hmem =>
        simp only [hmem, ↓reduceDIte]
        rcases howns with ⟨node, mono, hget, hmono, htail⟩
        refine ⟨node, mono, ?_, hmono, ?_⟩
        · rw [hsame nodeIndex hmem]
          exact hget
        · exact pairVecDivVHCChainOwns_congr_on node.next
            (owner.erase nodeIndex) nodes nodes' htail (by
              intro i hi
              exact hsame i (Finset.mem_of_mem_erase hi))
termination_by owner.card
decreasing_by
  exact Finset.card_erase_lt_of_mem (by assumption)

theorem pairVecDivVHCConsumeChain_owner_reclassified
    (this : DenseUPolyZp) (current : Option Nat)
    (owner unvisited : Finset Nat) (k : UInt64)
    (nodes : Array PairVecDivVHCNode) (lin : Array Nat) (resetH : Nat)
    (quotient divisor : SparsePolyZp) (result : PairVecDivVHCBucketResult)
    (howns : PairVecDivVHCChainOwns current owner nodes)
    (hownerSubset : owner ⊆ unvisited)
    (hrun : pairVecDivVHCConsumeChain this current unvisited k nodes lin
      resetH quotient divisor = .ok result) :
    owner ⊆ result.lin.toList.toFinset ∪ Finset.range result.resetH := by
  cases current with
  | none =>
      rw [PairVecDivVHCChainOwns] at howns
      subst owner
      simp
  | some nodeIndex =>
      rw [PairVecDivVHCChainOwns] at howns
      split at howns <;> try contradiction
      next hownerMem =>
        rcases howns with ⟨node, mono, hget, hmono, htailOwns⟩
        have hmem : nodeIndex ∈ unvisited := hownerSubset hownerMem
        rw [pairVecDivVHCConsumeChain] at hrun
        simp only [hmem, ↓reduceDIte] at hrun
        cases hconsume : pairVecDivVHCConsumeNode this nodeIndex k nodes lin
            resetH quotient divisor with
        | error fault => simp [hconsume] at hrun
        | ok step =>
            rcases step with ⟨k', nodes', lin', resetH', next⟩
            rw [hconsume] at hrun
            have hn : nodeIndex < nodes.size := by
              by_contra hnot
              rw [Array.getElem?_eq_none (by omega)] at hget
              contradiction
            have hnext := pairVecDivVHCConsumeNode_next_of_success this nodeIndex k k'
              nodes nodes' lin lin' resetH resetH' next quotient divisor hn
              hconsume
            have hnodeEq : nodes[nodeIndex] = node := by
              rw [Array.getElem?_eq_getElem hn] at hget
              exact Option.some.inj hget
            rw [hnodeEq] at hnext
            subst next
            have htailOwns' := pairVecDivVHCChainOwns_congr_on node.next
              (owner.erase nodeIndex) nodes nodes' htailOwns (by
                intro i hi
                exact pairVecDivVHCConsumeNode_get_ne this nodeIndex k k'
                  nodes nodes' lin lin' resetH resetH' node.next quotient
                  divisor hconsume i (Finset.mem_erase.mp hi).1.symm)
            have hownerErase : owner.erase nodeIndex ⊆
                unvisited.erase nodeIndex := by
              intro i hi
              exact Finset.mem_erase.mpr ⟨(Finset.mem_erase.mp hi).1,
                hownerSubset (Finset.mem_of_mem_erase hi)⟩
            have ih := pairVecDivVHCConsumeChain_owner_reclassified this
              node.next (owner.erase nodeIndex) (unvisited.erase nodeIndex) k'
              nodes' lin' resetH' quotient divisor result htailOwns'
              hownerErase hrun
            have hstep := pairVecDivVHCConsumeNode_location_progress this
              nodeIndex k k' nodes nodes' lin lin' resetH resetH' node.next
              quotient divisor hconsume
            have htailMonotone := pairVecDivVHCConsumeChain_location_monotone
              this node.next (unvisited.erase nodeIndex) k' nodes' lin'
              resetH' quotient divisor result hrun
            intro i hi
            by_cases heq : i = nodeIndex
            · subst i
              rcases hstep.2.2 with hlin | hreset
              · exact Finset.mem_union_left _ (htailMonotone.1 hlin)
              · exact Finset.mem_union_right _ (by
                  rw [Finset.mem_range]
                  exact Nat.lt_of_lt_of_le hreset htailMonotone.2)
            · exact ih (Finset.mem_erase.mpr ⟨heq, hi⟩)
termination_by owner.card
decreasing_by
  exact Finset.card_erase_lt_of_mem (by assumption)

theorem pairVecDivVHCConsumeChain_preserves_disjoint_owner
    (this : DenseUPolyZp) (current : Option Nat)
    (owner protectedSet unvisited : Finset Nat) (k : UInt64)
    (nodes : Array PairVecDivVHCNode) (lin : Array Nat) (resetH : Nat)
    (quotient divisor : SparsePolyZp) (result : PairVecDivVHCBucketResult)
    (howns : PairVecDivVHCChainOwns current owner nodes)
    (hownerSubset : owner ⊆ unvisited)
    (hprotectedSubset : protectedSet ⊆ unvisited)
    (hdisjoint : Disjoint protectedSet owner)
    (hrun : pairVecDivVHCConsumeChain this current unvisited k nodes lin
      resetH quotient divisor = .ok result) :
    protectedSet ⊆ result.unvisited ∧
      ∀ i ∈ protectedSet, result.nodes[i]? = nodes[i]? := by
  cases current with
  | none =>
      rw [pairVecDivVHCConsumeChain] at hrun
      simp only [Except.ok.injEq] at hrun
      subst result
      exact ⟨hprotectedSubset, by simp⟩
  | some nodeIndex =>
      rw [PairVecDivVHCChainOwns] at howns
      split at howns <;> try contradiction
      next hownerMem =>
        rcases howns with ⟨node, mono, hget, hmono, htailOwns⟩
        have hmem : nodeIndex ∈ unvisited := hownerSubset hownerMem
        rw [pairVecDivVHCConsumeChain] at hrun
        simp only [hmem, ↓reduceDIte] at hrun
        cases hconsume : pairVecDivVHCConsumeNode this nodeIndex k nodes lin
            resetH quotient divisor with
        | error fault => simp [hconsume] at hrun
        | ok step =>
            rcases step with ⟨k', nodes', lin', resetH', next⟩
            rw [hconsume] at hrun
            have hnext : next = node.next := by
              unfold pairVecDivVHCConsumeNode at hconsume
              have hn : nodeIndex < nodes.size := by
                by_contra hnot
                rw [Array.getElem?_eq_none (by omega)] at hget
                contradiction
              simp only [hn, ↓reduceDIte] at hconsume
              rw [Array.getElem?_eq_getElem hn] at hget
              simp only [Option.some.injEq] at hget
              subst node
              split at hconsume <;> try contradiction
              next hq =>
                split at hconsume <;> try contradiction
                next hd =>
                  split at hconsume
                  next hadvance =>
                    simp only [Except.ok.injEq, Prod.mk.injEq] at hconsume
                    exact hconsume.2.2.2.2.symm
                  next hadvance =>
                    split at hconsume <;> try contradiction
                    next hexhausted =>
                      split at hconsume <;> try contradiction
                      next horder =>
                        simp only [Except.ok.injEq, Prod.mk.injEq] at hconsume
                        exact hconsume.2.2.2.2.symm
            subst next
            have htailOwns' := pairVecDivVHCChainOwns_congr_on node.next
              (owner.erase nodeIndex) nodes nodes' htailOwns (by
                intro i hi
                exact (pairVecDivVHCConsumeNode_get_ne this nodeIndex k k'
                  nodes nodes' lin lin' resetH resetH' node.next quotient
                  divisor hconsume i (Finset.mem_erase.mp hi).1.symm))
            have hnodeNotProtected : nodeIndex ∉ protectedSet := by
              intro hp
              exact Finset.disjoint_left.mp hdisjoint hp hownerMem
            have hprotectedErase : protectedSet ⊆ unvisited.erase nodeIndex := by
              intro i hi
              exact Finset.mem_erase.mpr ⟨by
                intro heq
                subst i
                exact hnodeNotProtected hi, hprotectedSubset hi⟩
            have hownerErase : owner.erase nodeIndex ⊆
                unvisited.erase nodeIndex := by
              intro i hi
              exact Finset.mem_erase.mpr ⟨(Finset.mem_erase.mp hi).1,
                hownerSubset (Finset.mem_of_mem_erase hi)⟩
            have hdisjointErase : Disjoint protectedSet
                (owner.erase nodeIndex) := by
              exact Finset.disjoint_left.mpr (by
                intro i hp ho
                exact Finset.disjoint_left.mp hdisjoint hp
                  (Finset.mem_of_mem_erase ho))
            have ih := pairVecDivVHCConsumeChain_preserves_disjoint_owner this
              node.next (owner.erase nodeIndex) protectedSet
              (unvisited.erase nodeIndex) k' nodes' lin' resetH' quotient
              divisor result htailOwns' hownerErase hprotectedErase
              hdisjointErase hrun
            refine ⟨ih.1, ?_⟩
            intro i hi
            rw [ih.2 i hi]
            exact pairVecDivVHCConsumeNode_get_ne this nodeIndex k k' nodes
              nodes' lin lin' resetH resetH' node.next quotient divisor
              hconsume i (by
                intro heq
                subst i
                exact hnodeNotProtected hi)
termination_by owner.card
decreasing_by
  exact Finset.card_erase_lt_of_mem (by assumption)

theorem pairVecDivVHCChainOwns_valid
    (current : Option Nat) (owner : Finset Nat)
    (nodes : Array PairVecDivVHCNode)
    (howns : PairVecDivVHCChainOwns current owner nodes) :
    PairVecDivVHCChainValid current owner nodes := by
  cases current with
  | none => simp [PairVecDivVHCChainValid]
  | some nodeIndex =>
      rw [PairVecDivVHCChainOwns] at howns
      rw [PairVecDivVHCChainValid]
      split at howns <;> try contradiction
      next hmem =>
        simp only [hmem, ↓reduceDIte]
        rcases howns with ⟨node, mono, hget, hmono, htail⟩
        exact ⟨node, mono, hget, hmono,
          pairVecDivVHCChainOwns_valid node.next (owner.erase nodeIndex)
            nodes htail⟩
termination_by owner.card
decreasing_by
  exact Finset.card_erase_lt_of_mem (by assumption)

theorem pairVecDivVHCChainOwns_subset_range
    (current : Option Nat) (owner : Finset Nat)
    (nodes : Array PairVecDivVHCNode)
    (howns : PairVecDivVHCChainOwns current owner nodes) :
    owner ⊆ Finset.range nodes.size := by
  cases current with
  | none =>
      simp [PairVecDivVHCChainOwns] at howns
      subst owner
      simp
  | some nodeIndex =>
      rw [PairVecDivVHCChainOwns] at howns
      split at howns <;> try contradiction
      next hmem =>
        rcases howns with ⟨node, mono, hget, hmono, htail⟩
        have hn : nodeIndex < nodes.size := by
          by_contra hnot
          rw [Array.getElem?_eq_none (by omega)] at hget
          contradiction
        have htailSubset := pairVecDivVHCChainOwns_subset_range node.next
          (owner.erase nodeIndex) nodes htail
        intro i hi
        by_cases heq : i = nodeIndex
        · subst i
          simpa using hn
        · exact htailSubset (Finset.mem_erase.mpr ⟨heq, hi⟩)
termination_by owner.card
decreasing_by
  exact Finset.card_erase_lt_of_mem (by assumption)

theorem pairVecDivVHCChainOwns_mem_active
    (current : Option Nat) (owner : Finset Nat)
    (nodes : Array PairVecDivVHCNode)
    (howns : PairVecDivVHCChainOwns current owner nodes)
    (i : Nat) (hi : i ∈ owner) :
    ∃ node mono, nodes[i]? = some node ∧ node.mono = some mono := by
  cases current with
  | none =>
      simp [PairVecDivVHCChainOwns] at howns
      subst owner
      simp at hi
  | some nodeIndex =>
      rw [PairVecDivVHCChainOwns] at howns
      split at howns <;> try contradiction
      next hmem =>
        rcases howns with ⟨node, mono, hget, hmono, htail⟩
        by_cases heq : i = nodeIndex
        · subst i
          exact ⟨node, mono, hget, hmono⟩
        · exact pairVecDivVHCChainOwns_mem_active node.next
            (owner.erase nodeIndex) nodes htail i
            (Finset.mem_erase.mpr ⟨heq, hi⟩)
termination_by owner.card
decreasing_by
  exact Finset.card_erase_lt_of_mem (by assumption)

theorem pairVecDivVHCChainOwns_mem_degree
    (current : Option Nat) (owner : Finset Nat)
    (nodes : Array PairVecDivVHCNode) (degree i : Nat)
    (howns : PairVecDivVHCChainOwns current owner nodes)
    (hdegree : PairVecDivVHCChainAtDegree current owner nodes degree)
    (hi : i ∈ owner) :
    ∃ node mono, nodes[i]? = some node ∧ node.mono = some mono ∧
      mono.deg = degree := by
  cases current with
  | none =>
      simp [PairVecDivVHCChainOwns] at howns
      subst owner
      simp at hi
  | some nodeIndex =>
      rw [PairVecDivVHCChainOwns] at howns
      rw [PairVecDivVHCChainAtDegree] at hdegree
      split at howns <;> try contradiction
      next hmem =>
        simp only [hmem, ↓reduceDIte] at hdegree
        rcases howns with ⟨node, mono, hget, hmono, htailOwns⟩
        rcases hdegree with
          ⟨degreeNode, degreeMono, hdegreeGet, hdegreeMono,
            hmonoDegree, htailDegree⟩
        rw [hget] at hdegreeGet
        simp only [Option.some.injEq] at hdegreeGet
        subst degreeNode
        rw [hmono] at hdegreeMono
        simp only [Option.some.injEq] at hdegreeMono
        subst degreeMono
        by_cases heq : i = nodeIndex
        · subst i
          exact ⟨node, mono, hget, hmono, hmonoDegree⟩
        · exact pairVecDivVHCChainOwns_mem_degree node.next
            (owner.erase nodeIndex) nodes degree i htailOwns htailDegree
            (Finset.mem_erase.mpr ⟨heq, hi⟩)
termination_by owner.card
decreasing_by
  exact Finset.card_erase_lt_of_mem (by assumption)

theorem pairVecDivVHCChainOwns_head_mem
    (nodeIndex : Nat) (owner : Finset Nat)
    (nodes : Array PairVecDivVHCNode)
    (howns : PairVecDivVHCChainOwns (some nodeIndex) owner nodes) :
    nodeIndex ∈ owner := by
  rw [PairVecDivVHCChainOwns] at howns
  split at howns
  next hmem => exact hmem
  next hmem => contradiction

/-- Linking a fresh node in front of an existing exact chain transfers
ownership to the inserted head.  This is the ownership effect of the real
`VHC::next` field write used by heap insertion. -/
theorem pairVecDivVHCSetNext_chainOwns_insert
    (nodeIndex : Nat) (next : Option Nat) (nodes nodes' : Array PairVecDivVHCNode)
    (node : PairVecDivVHCNode) (mono : UMonomial) (tailOwner : Finset Nat)
    (hnode : nodes[nodeIndex]? = some node) (hmono : node.mono = some mono)
    (hfresh : nodeIndex ∉ tailOwner)
    (htail : PairVecDivVHCChainOwns next tailOwner nodes)
    (hrun : pairVecDivVHCSetNext nodeIndex next nodes = .ok nodes') :
    PairVecDivVHCChainOwns (some nodeIndex) (insert nodeIndex tailOwner)
      nodes' := by
  unfold pairVecDivVHCSetNext at hrun
  split at hrun <;> try contradiction
  next hn =>
    simp only [Except.ok.injEq] at hrun
    subst nodes'
    rw [Array.getElem?_eq_getElem hn] at hnode
    simp only [Option.some.injEq] at hnode
    subst node
    rw [PairVecDivVHCChainOwns]
    simp only [Finset.mem_insert, true_or, ↓reduceDIte]
    refine ⟨{ nodes[nodeIndex] with next := next }, mono, ?_, hmono, ?_⟩
    · rw [Array.getElem?_set_self hn]
    · have htail' := pairVecDivVHCChainOwns_congr_on next tailOwner nodes
        (nodes.set nodeIndex { nodes[nodeIndex] with next := next }) htail (by
          intro i hi
          rw [Array.getElem?_set_ne hn]
          intro heq
          subst i
          exact hfresh hi)
      simpa [hfresh] using htail'

theorem pairVecDivVHCSetNext_chainAtDegree_insert
    (nodeIndex : Nat) (next : Option Nat) (nodes nodes' : Array PairVecDivVHCNode)
    (node : PairVecDivVHCNode) (mono : UMonomial) (tailOwner : Finset Nat)
    (degree : Nat) (hnode : nodes[nodeIndex]? = some node)
    (hmono : node.mono = some mono) (hmonoDegree : mono.deg = degree)
    (hfresh : nodeIndex ∉ tailOwner)
    (htail : PairVecDivVHCChainAtDegree next tailOwner nodes degree)
    (hrun : pairVecDivVHCSetNext nodeIndex next nodes = .ok nodes') :
    PairVecDivVHCChainAtDegree (some nodeIndex) (insert nodeIndex tailOwner)
      nodes' degree := by
  unfold pairVecDivVHCSetNext at hrun
  split at hrun <;> try contradiction
  next hn =>
    simp only [Except.ok.injEq] at hrun
    subst nodes'
    rw [Array.getElem?_eq_getElem hn] at hnode
    simp only [Option.some.injEq] at hnode
    subst node
    rw [PairVecDivVHCChainAtDegree]
    simp only [Finset.mem_insert, true_or, ↓reduceDIte]
    refine ⟨{ nodes[nodeIndex] with next := next }, mono, ?_, hmono,
      hmonoDegree, ?_⟩
    · rw [Array.getElem?_set_self hn]
    · have htail' := pairVecDivVHCChainAtDegree_congr_on next tailOwner nodes
        (nodes.set nodeIndex { nodes[nodeIndex] with next := next }) degree htail
        (by
          intro i hi
          rw [Array.getElem?_set_ne hn]
          intro heq
          subst i
          exact hfresh hi)
      simpa [hfresh] using htail'

theorem pairVecDivVHCSetNext_preserves_disjoint_chain
    (nodeIndex : Nat) (next : Option Nat) (nodes nodes' : Array PairVecDivVHCNode)
    (current : Option Nat) (owner : Finset Nat)
    (hfresh : nodeIndex ∉ owner)
    (howns : PairVecDivVHCChainOwns current owner nodes)
    (hrun : pairVecDivVHCSetNext nodeIndex next nodes = .ok nodes') :
    PairVecDivVHCChainOwns current owner nodes' := by
  unfold pairVecDivVHCSetNext at hrun
  split at hrun <;> try contradiction
  next hn =>
    simp only [Except.ok.injEq] at hrun
    subst nodes'
    exact pairVecDivVHCChainOwns_congr_on current owner nodes
      (nodes.set nodeIndex { nodes[nodeIndex] with next := next }) howns (by
        intro i hi
        rw [Array.getElem?_set_ne hn]
        intro heq
        subst i
        exact hfresh hi)

/-- Exact ownership of every active heap bucket, keyed by its head node.
The source heap contains each head once, and different heads own disjoint
`next` chains. -/
def PairVecDivVHCHeapChainOwnership (heap : Array Nat)
    (owners : Nat → Finset Nat) (nodes : Array PairVecDivVHCNode) : Prop :=
  (∀ (slot head : Nat), heap[slot]? = some head →
      PairVecDivVHCChainOwns (some head) (owners head) nodes) ∧
    (∀ (left right head : Nat), heap[left]? = some head →
      heap[right]? = some head → left = right) ∧
    (∀ (left right leftHead rightHead : Nat), heap[left]? = some leftHead →
      heap[right]? = some rightHead → leftHead ≠ rightHead →
      Disjoint (owners leftHead) (owners rightHead))

/-- The finite union of node owners named by the concrete active heap heads. -/
def PairVecDivVHCHeapOwnedNodes (heap : Array Nat)
    (owners : Nat → Finset Nat) : Finset Nat :=
  heap.toList.toFinset.biUnion owners

def PairVecDivVHCNodeAtDegree (degree i : Nat)
    (nodes : Array PairVecDivVHCNode) : Prop :=
  match pairVecDivVHCMono i nodes with
  | .ok mono => mono.deg = degree
  | .error _ => False

instance (degree i : Nat) (nodes : Array PairVecDivVHCNode) :
    Decidable (PairVecDivVHCNodeAtDegree degree i nodes) := by
  unfold PairVecDivVHCNodeAtDegree
  split <;> infer_instance

theorem pairVecDivVHCNodeAtDegree_iff (degree i : Nat)
    (nodes : Array PairVecDivVHCNode) :
    PairVecDivVHCNodeAtDegree degree i nodes ↔
      ∃ mono, pairVecDivVHCMono i nodes = .ok mono ∧ mono.deg = degree := by
  unfold PairVecDivVHCNodeAtDegree
  cases hmono : pairVecDivVHCMono i nodes with
  | error fault => simp [hmono]
  | ok mono => simp [hmono]

def PairVecDivVHCHeapOwnedNodesAtDegree (degree : Nat)
    (heap : Array Nat) (owners : Nat → Finset Nat)
    (nodes : Array PairVecDivVHCNode) : Finset Nat :=
  (PairVecDivVHCHeapOwnedNodes heap owners).filter fun i =>
    PairVecDivVHCNodeAtDegree degree i nodes

/-- A real extract removes exactly the old root owner from the finite union of
heap-owned nodes. -/
theorem pairVecDivVHCExtract_heapOwnedNodes
    (heap heap' : Array Nat) (sourceNodes nodes : Array PairVecDivVHCNode)
    (owners : Nat → Finset Nat) (hheap : 0 < heap.size)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners sourceNodes)
    (hrun : pairVecDivVHCExtract heap nodes = .ok heap') :
    PairVecDivVHCHeapOwnedNodes heap' owners =
      PairVecDivVHCHeapOwnedNodes heap owners \ owners heap[0] := by
  apply Finset.ext
  intro i
  simp only [PairVecDivVHCHeapOwnedNodes, Finset.mem_sdiff,
    Finset.mem_biUnion, List.mem_toFinset]
  constructor
  · rintro ⟨head, hheadMem, hi⟩
    rcases List.getElem_of_mem hheadMem with ⟨slot, hslot, hvalue⟩
    have hslotBound : slot < heap'.size := by simpa using hslot
    have harrayValue : heap'[slot] = head := by
      rw [← Array.getElem_toList hslotBound]
      exact hvalue
    have hheadGet : heap'[slot]? = some head := by
      rw [Array.getElem?_eq_getElem hslotBound, harrayValue]
    rcases pairVecDivVHCExtract_valuesFrom heap heap' nodes hrun slot head
        hheadGet with ⟨oldSlot, holdGet⟩
    have holdBound : oldSlot < heap.size := by
      by_contra hn
      rw [Array.getElem?_eq_none (by omega)] at holdGet
      contradiction
    have hrootNe : head ≠ heap[0] := by
      intro heq
      have hrootOnly : PairVecDivVHCOnlyAt heap heap[0] 0 := by
        intro current hget
        exact hownership.2.1 current 0 heap[0] hget
          (Array.getElem?_eq_getElem hheap)
      have hgone := pairVecDivVHCExtract_excludes_root heap heap' nodes hheap
        hrootOnly hrun slot
      apply hgone
      rw [hheadGet, heq]
    refine ⟨⟨head, ?_, hi⟩, ?_⟩
    · have hvalueEq : heap[oldSlot] = head := by
        rw [Array.getElem?_eq_getElem holdBound] at holdGet
        exact Option.some.inj holdGet
      rw [← hvalueEq]
      exact Array.getElem_mem_toList holdBound
    · have hdisjoint : Disjoint (owners head) (owners heap[0]) :=
        hownership.2.2 oldSlot 0 head heap[0] holdGet
          (Array.getElem?_eq_getElem hheap) hrootNe
      exact fun hrootMem => Finset.disjoint_left.mp hdisjoint hi hrootMem
  · rintro ⟨⟨head, hheadMem, hi⟩, hnotRoot⟩
    rcases List.getElem_of_mem hheadMem with ⟨slot, hslot, hvalue⟩
    have hslotBound : slot < heap.size := by simpa using hslot
    have harrayValue : heap[slot] = head := by
      rw [← Array.getElem_toList hslotBound]
      exact hvalue
    have hheadGet : heap[slot]? = some head := by
      rw [Array.getElem?_eq_getElem hslotBound, harrayValue]
    have hrootNe : head ≠ heap[0] := by
      intro heq
      apply hnotRoot
      rw [← heq]
      exact hi
    rcases pairVecDivVHCExtract_preserves_nonroot_head heap heap' nodes hheap
        hownership.2.1 hrun slot head hheadGet hrootNe with
      ⟨newSlot, hnewGet⟩
    have hnewBound : newSlot < heap'.size := by
      by_contra hn
      rw [Array.getElem?_eq_none (by omega)] at hnewGet
      contradiction
    refine ⟨head, ?_, hi⟩
    have hvalueEq : heap'[newSlot] = head := by
      rw [Array.getElem?_eq_getElem hnewBound] at hnewGet
      exact Option.some.inj hnewGet
    rw [← hvalueEq]
    exact Array.getElem_mem_toList hnewBound

theorem PairVecDivVHCHeapChainOwnership.heapPointersValid
    (heap : Array Nat) (owners : Nat → Finset Nat)
    (nodes : Array PairVecDivVHCNode)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes) :
    PairVecDivVHCHeapPointersValid heap nodes := by
  intro slot hslot
  let head := heap[slot]
  have hheap : heap[slot]? = some head := Array.getElem?_eq_getElem hslot
  have howns := hownership.1 slot head hheap
  have hheadMem := pairVecDivVHCChainOwns_head_mem head (owners head) nodes
    howns
  rcases pairVecDivVHCChainOwns_mem_active (some head) (owners head) nodes
      howns head hheadMem with ⟨node, mono, hnode, hmono⟩
  exact ⟨head, node, mono, hheap, hnode, hmono⟩

theorem pairVecDivVHCHeapChainOwnership_owner_eq_at
    (heap : Array Nat) (leftOwners rightOwners : Nat → Finset Nat)
    (nodes : Array PairVecDivVHCNode)
    (hleft : PairVecDivVHCHeapChainOwnership heap leftOwners nodes)
    (hright : PairVecDivVHCHeapChainOwnership heap rightOwners nodes)
    (slot head : Nat) (hget : heap[slot]? = some head) :
    leftOwners head = rightOwners head :=
  pairVecDivVHCChainOwns_unique (some head) (leftOwners head)
    (rightOwners head) nodes (hleft.1 slot head hget)
    (hright.1 slot head hget)

/-- Total location coverage of the allocated divisor-tail cursor block.
Every node is either in the exhausted `reset_h` prefix, on the temporary
reverse-reinsertion stack, or owned by one concrete heap bucket.  Uniqueness
and separation of these regions are supplied by `ResetReady`, `LinReady`, and
`HeapChainsOwnedAway`; this predicate records the missing totality direction. -/
def PairVecDivVHCNodesCovered (heap : Array Nat)
    (owners : Nat → Finset Nat) (lin : Array Nat) (resetH : Nat)
    (nodes : Array PairVecDivVHCNode) : Prop :=
  ∀ i : Nat, i < nodes.size →
    i < resetH ∨ i ∈ lin.toList.toFinset ∨
      ∃ (slot head : Nat), heap[slot]? = some head ∧ i ∈ owners head

theorem PairVecDivVHCNodesCovered.congr_owners
    (heap : Array Nat) (leftOwners rightOwners : Nat → Finset Nat)
    (lin : Array Nat) (resetH : Nat) (nodes : Array PairVecDivVHCNode)
    (hleftOwnership :
      PairVecDivVHCHeapChainOwnership heap leftOwners nodes)
    (hrightOwnership :
      PairVecDivVHCHeapChainOwnership heap rightOwners nodes)
    (hcovered : PairVecDivVHCNodesCovered heap leftOwners lin resetH nodes) :
    PairVecDivVHCNodesCovered heap rightOwners lin resetH nodes := by
  intro i hi
  rcases hcovered i hi with hreset | hlin | ⟨slot, head, hget, hmem⟩
  · exact Or.inl hreset
  · exact Or.inr (Or.inl hlin)
  · refine Or.inr (Or.inr ⟨slot, head, hget, ?_⟩)
    rw [← pairVecDivVHCHeapChainOwnership_owner_eq_at heap leftOwners
      rightOwners nodes hleftOwnership hrightOwnership slot head hget]
    exact hmem

/-- A single ownership witness must justify both the concrete heap chains and
total coverage of the allocated node block.  Keeping the existential outside
the conjunction is essential: heap insertion can replace a bucket head and
therefore changes the owner map, so two unrelated existential witnesses would
not be strong enough to transport coverage through the real execution. -/
def PairVecDivVHCStateCovered (heap : Array Nat) (nodes : Array PairVecDivVHCNode)
    (lin : Array Nat) (resetH : Nat) : Prop :=
  ∃ owners : Nat → Finset Nat,
    PairVecDivVHCHeapChainOwnership heap owners nodes ∧
      PairVecDivVHCNodesCovered heap owners lin resetH nodes

theorem PairVecDivVHCStateCovered.covered_with
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode)
    (lin : Array Nat) (resetH : Nat) (owners : Nat → Finset Nat)
    (hstate : PairVecDivVHCStateCovered heap nodes lin resetH)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes) :
    PairVecDivVHCNodesCovered heap owners lin resetH nodes := by
  rcases hstate with ⟨stateOwners, hstateOwnership, hcovered⟩
  exact hcovered.congr_owners heap stateOwners owners lin resetH nodes
    hstateOwnership hownership

theorem pairVecDivVHCInit_nodesCovered (divisor : SparsePolyZp)
    (owners : Nat → Finset Nat) :
    PairVecDivVHCNodesCovered #[] owners #[] (divisor.size - 1)
      (pairVecDivVHCInit divisor) := by
  intro i hi
  left
  rw [pairVecDivVHCInit_size] at hi
  exact hi

theorem pairVecDivVHCInit_stateCovered (divisor : SparsePolyZp) :
    PairVecDivVHCStateCovered #[] (pairVecDivVHCInit divisor) #[]
      (divisor.size - 1) := by
  refine ⟨fun _ => ∅, ?_, pairVecDivVHCInit_nodesCovered divisor _⟩
  constructor
  · simp
  constructor <;> simp

theorem PairVecDivVHCNodesCovered.owned_of_active_not_lin
    (heap : Array Nat) (owners : Nat → Finset Nat)
    (lin : Array Nat) (resetH quotientSize : Nat)
    (nodes : Array PairVecDivVHCNode) (nodeIndex : Nat)
    (node : PairVecDivVHCNode)
    (hcovered : PairVecDivVHCNodesCovered heap owners lin resetH nodes)
    (hready : PairVecDivVHCResetReady resetH quotientSize nodes)
    (hget : nodes[nodeIndex]? = some node) (hactive : node.mono ≠ none)
    (hnotLin : nodeIndex ∉ lin.toList.toFinset) :
    ∃ (slot head : Nat), heap[slot]? = some head ∧
      nodeIndex ∈ owners head := by
  have hn : nodeIndex < nodes.size := by
    by_contra hnot
    rw [Array.getElem?_eq_none (by omega)] at hget
    contradiction
  rcases hcovered nodeIndex hn with hreset | hlin | howned
  · rcases hready.2 nodeIndex hreset with
      ⟨readyNode, hreadyGet, hquotientIndex, hdivisorIndex, hmono⟩
    rw [hget] at hreadyGet
    simp only [Option.some.injEq] at hreadyGet
    subst readyNode
    exact (hactive hmono).elim
  · exact (hnotLin hlin).elim
  · exact howned

/-- Every currently heap-owned bucket is homogeneous at the degree reported
by its actual head read.  This deliberately excludes nodes already moved to
`lin`, whose stale `next` fields are not observable by the source heap. -/
def PairVecDivVHCHeapChainsHomogeneous (heap : Array Nat)
    (owners : Nat → Finset Nat) (nodes : Array PairVecDivVHCNode) : Prop :=
  ∀ (slot head : Nat) (mono : UMonomial), heap[slot]? = some head →
    pairVecDivVHCMono head nodes = .ok mono →
      PairVecDivVHCChainAtDegree (some head) (owners head) nodes mono.deg

theorem PairVecDivVHCHeapChainsHomogeneous.of_valuesFrom
    (source target : Array Nat) (owners : Nat → Finset Nat)
    (nodes : Array PairVecDivVHCNode)
    (hhomogeneous : PairVecDivVHCHeapChainsHomogeneous source owners nodes)
    (hfrom : PairVecDivVHCValuesFrom target source) :
    PairVecDivVHCHeapChainsHomogeneous target owners nodes := by
  intro slot head mono hget hmono
  rcases hfrom slot head hget with ⟨sourceSlot, hsource⟩
  exact hhomogeneous sourceSlot head mono hsource hmono

theorem PairVecDivVHCHeapChainsHomogeneous.congr_owners
    (heap : Array Nat) (leftOwners rightOwners : Nat → Finset Nat)
    (nodes : Array PairVecDivVHCNode)
    (hleftOwnership : PairVecDivVHCHeapChainOwnership heap leftOwners nodes)
    (hrightOwnership : PairVecDivVHCHeapChainOwnership heap rightOwners nodes)
    (hhomogeneous : PairVecDivVHCHeapChainsHomogeneous heap leftOwners nodes) :
    PairVecDivVHCHeapChainsHomogeneous heap rightOwners nodes := by
  intro slot head mono hget hmono
  rw [← pairVecDivVHCHeapChainOwnership_owner_eq_at heap leftOwners
    rightOwners nodes hleftOwnership hrightOwnership slot head hget]
  exact hhomogeneous slot head mono hget hmono

theorem pairVecDivVHCRootOwner_subset_ownedNodesAtDegree
    (degree : Nat) (heap : Array Nat) (owners : Nat → Finset Nat)
    (nodes : Array PairVecDivVHCNode) (rootMono : UMonomial)
    (hheap : 0 < heap.size)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hhomogeneous : PairVecDivVHCHeapChainsHomogeneous heap owners nodes)
    (hrootMono : pairVecDivVHCMono heap[0] nodes = .ok rootMono)
    (hdegree : rootMono.deg = degree) :
    owners heap[0] ⊆
      PairVecDivVHCHeapOwnedNodesAtDegree degree heap owners nodes := by
  intro i hi
  have hrootGet : heap[0]? = some heap[0] :=
    Array.getElem?_eq_getElem hheap
  have howns := hownership.1 0 heap[0] hrootGet
  have hchain := hhomogeneous 0 heap[0] rootMono hrootGet hrootMono
  rcases pairVecDivVHCChainOwns_mem_degree (some heap[0]) (owners heap[0])
      nodes rootMono.deg i howns hchain hi with
    ⟨node, mono, hnode, hmono, hmonoDegree⟩
  simp only [PairVecDivVHCHeapOwnedNodesAtDegree, Finset.mem_filter,
    pairVecDivVHCNodeAtDegree_iff]
  refine ⟨?_, mono, ?_, ?_⟩
  · simp only [PairVecDivVHCHeapOwnedNodes, Finset.mem_biUnion,
      List.mem_toFinset]
    exact ⟨heap[0], Array.getElem_mem_toList hheap, hi⟩
  · exact (pairVecDivVHCMono_eq_ok_iff i nodes mono).2
      ⟨node, hnode, hmono⟩
  · omega

theorem pairVecDivVHCHeapOwnedNodesAtDegree_eq_empty_of_root_ne
    (degree : Nat) (heap : Array Nat) (owners : Nat → Finset Nat)
    (nodes : Array PairVecDivVHCNode) (rootMono : UMonomial)
    (hheap : 0 < heap.size)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hhomogeneous : PairVecDivVHCHeapChainsHomogeneous heap owners nodes)
    (hordered : PairVecDivVHCHeapOrdered heap nodes)
    (hrootMono : pairVecDivVHCMono heap[0] nodes = .ok rootMono)
    (hmax : ∀ (slot head : Nat) (mono : UMonomial),
      heap[slot]? = some head → pairVecDivVHCMono head nodes = .ok mono →
        mono.deg ≤ degree)
    (hne : rootMono.deg ≠ degree) :
    PairVecDivVHCHeapOwnedNodesAtDegree degree heap owners nodes = ∅ := by
  apply Finset.ext
  intro i
  simp
  intro hi
  simp only [PairVecDivVHCHeapOwnedNodesAtDegree, Finset.mem_filter,
    PairVecDivVHCHeapOwnedNodes, Finset.mem_biUnion, List.mem_toFinset,
    pairVecDivVHCNodeAtDegree_iff] at hi
  rcases hi with ⟨⟨head, hheadMem, hiOwner⟩, nodeMono, hnodeMono,
    hnodeDegree⟩
  rcases List.getElem_of_mem hheadMem with ⟨slot, hslot, hvalue⟩
  have hslotBound : slot < heap.size := by simpa using hslot
  have hheadEq : heap[slot] = head := by
    rw [← Array.getElem_toList hslotBound]
    exact hvalue
  have hheadGet : heap[slot]? = some head := by
    rw [Array.getElem?_eq_getElem hslotBound, hheadEq]
  rcases hownership.heapPointersValid heap owners nodes slot hslotBound with
    ⟨validHead, headNode, headMono, hvalidGet, hheadNode, hheadActive⟩
  rw [hheadGet] at hvalidGet
  have hvalidHead : validHead = head := (Option.some.inj hvalidGet).symm
  subst validHead
  have hheadMono : pairVecDivVHCMono head nodes = .ok headMono :=
    (pairVecDivVHCMono_eq_ok_iff head nodes headMono).2
      ⟨headNode, hheadNode, hheadActive⟩
  have howns := hownership.1 slot head hheadGet
  have hchain := hhomogeneous slot head headMono hheadGet hheadMono
  rcases pairVecDivVHCChainOwns_mem_degree (some head) (owners head) nodes
      headMono.deg i howns hchain hiOwner with
    ⟨ownedNode, ownedMono, hownedNode, hownedMono, hownedDegree⟩
  have hmonoEq : ownedMono = nodeMono := by
    have hownedRun := (pairVecDivVHCMono_eq_ok_iff i nodes ownedMono).2
      ⟨ownedNode, hownedNode, hownedMono⟩
    rw [hnodeMono] at hownedRun
    exact (Except.ok.inj hownedRun).symm
  have hheadDegree : headMono.deg = degree := by
    rw [hmonoEq, hnodeDegree] at hownedDegree
    exact hownedDegree.symm
  have hslotRun : pairVecDivVHCMono heap[slot] nodes = .ok headMono := by
    simpa [hheadEq] using hheadMono
  have hleRoot := pairVecDivVHCHeapOrdered_slot_le_root heap nodes
    (hownership.heapPointersValid heap owners nodes) hordered slot hslotBound
    headMono rootMono hslotRun hrootMono
  have hrootLe := hmax 0 heap[0] rootMono
    (Array.getElem?_eq_getElem hheap) hrootMono
  exact hne (by omega)

theorem pairVecDivVHCOwnedNode_degree_le_frontier
    (dividendIndex : Nat) (dividend : SparsePolyZp)
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode)
    (owners : Nat → Finset Nat) (frontier : PairVecDivVHCFrontier)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hhomogeneous : PairVecDivVHCHeapChainsHomogeneous heap owners nodes)
    (hordered : PairVecDivVHCHeapOrdered heap nodes)
    (slot head i : Nat) (node : PairVecDivVHCNode)
    (hheap : heap[slot]? = some head) (hmem : i ∈ owners head)
    (hget : nodes[i]? = some node)
    (hselect : pairVecDivVHCSelectFrontier dividendIndex dividend heap nodes =
      .ok frontier) :
    ∃ mono, node.mono = some mono ∧ mono.deg ≤ frontier.degree := by
  have hvalid := hownership.heapPointersValid heap owners nodes
  have hslot : slot < heap.size := by
    by_contra hnot
    rw [Array.getElem?_eq_none (by omega)] at hheap
    contradiction
  rcases hvalid slot hslot with
    ⟨validHead, headNode, headMono, hvalidHeap, hheadNode, hheadActive⟩
  rw [hheap] at hvalidHeap
  simp only [Option.some.injEq] at hvalidHeap
  subst validHead
  have hheadRun : pairVecDivVHCMono head nodes = .ok headMono := by
    apply (pairVecDivVHCMono_eq_ok_iff head nodes headMono).mpr
    exact ⟨headNode, hheadNode, hheadActive⟩
  have hchainDegree := hhomogeneous slot head headMono hheap hheadRun
  have howns := hownership.1 slot head hheap
  rcases pairVecDivVHCChainOwns_mem_degree (some head) (owners head) nodes
      headMono.deg i howns hchainDegree hmem with
    ⟨ownedNode, mono, hownedGet, hmono, hmonoDegree⟩
  rw [hget] at hownedGet
  simp only [Option.some.injEq] at hownedGet
  subst ownedNode
  refine ⟨mono, hmono, ?_⟩
  rw [hmonoDegree]
  have hslotRun : pairVecDivVHCMono heap[slot] nodes = .ok headMono := by
    have hheadEq : heap[slot] = head := by
      rw [Array.getElem?_eq_getElem hslot] at hheap
      exact Option.some.inj hheap
    rw [hheadEq]
    exact hheadRun
  exact pairVecDivVHCSelectFrontier_heap_slot_degree_le dividendIndex dividend
    heap nodes frontier hvalid hordered slot hslot headMono hslotRun hselect

theorem pairVecDivVHCOwnedRow_productAtFrontier_eq_cursor
    (p degreeLimit dividendIndex : Nat) (dividend : SparsePolyZp)
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode)
    (owners : Nat → Finset Nat) (frontier : PairVecDivVHCFrontier)
    (quotient divisor : SparsePolyZp) (slot head i q : Nat)
    (node : PairVecDivVHCNode)
    (quotientTerm divisorTerm : UMonomial × Zp)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hhomogeneous : PairVecDivVHCHeapChainsHomogeneous heap owners nodes)
    (hordered : PairVecDivVHCHeapOrdered heap nodes)
    (hcanonical : SparsePolyZp.Canonical p quotient)
    (hprefix : PairVecDivVHCCursorPrefixAbove degreeLimit nodes quotient divisor)
    (hfrontier : frontier.degree < degreeLimit)
    (hheap : heap[slot]? = some head) (hmem : i ∈ owners head)
    (hget : nodes[i]? = some node)
    (hdenotes : PairVecDivVHCNodeDenotes quotient divisor node)
    (hquotient : quotient[q]? = some quotientTerm)
    (hdivisor : divisor[node.divisorIndex]? = some divisorTerm)
    (hdegree : quotientTerm.1.deg + divisorTerm.1.deg = frontier.degree)
    (hselect : pairVecDivVHCSelectFrontier dividendIndex dividend heap nodes =
      .ok frontier) :
    q = node.quotientIndex := by
  rcases pairVecDivVHCOwnedNode_degree_le_frontier dividendIndex dividend heap
      nodes owners frontier hownership hhomogeneous hordered slot head i
      node hheap hmem hget hselect with ⟨mono, hmono, hmonoLe⟩
  exact pairVecDivVHCProductAtFrontier_eq_cursor_of_mono_le p degreeLimit
    frontier.degree nodes quotient divisor i q node mono quotientTerm
    divisorTerm hcanonical hprefix hfrontier hget hdenotes hmono hmonoLe
    hquotient hdivisor hdegree

/-- At an outer-loop state with no deferred `lin` nodes, every concrete
divisor-tail product at the selected frontier belongs to a live heap bucket
and is exactly at that row's current quotient cursor. -/
theorem pairVecDivVHCFrontierProduct_owned_cursor
    (p degreeLimit dividendIndex resetH d q : Nat)
    (dividend quotient divisor : SparsePolyZp)
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode)
    (owners : Nat → Finset Nat) (frontier : PairVecDivVHCFrontier)
    (quotientTerm divisorTerm : UMonomial × Zp)
    (hsize : nodes.size = divisor.size - 1)
    (hfixed : PairVecDivVHCNodeDivisorIndicesFixed nodes)
    (hstate : PairVecDivVHCStateCovered heap nodes #[] resetH)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hhomogeneous : PairVecDivVHCHeapChainsHomogeneous heap owners nodes)
    (hresetReady : PairVecDivVHCResetReady resetH quotient.size nodes)
    (hordered : PairVecDivVHCHeapOrdered heap nodes)
    (hdenotes : ∀ (i : Nat) (node : PairVecDivVHCNode),
      nodes[i]? = some node → node.mono ≠ none →
        PairVecDivVHCNodeDenotes quotient divisor node)
    (hcanonical : SparsePolyZp.Canonical p quotient)
    (hprefix : PairVecDivVHCCursorPrefixAbove degreeLimit nodes quotient divisor)
    (hfrontier : frontier.degree < degreeLimit)
    (hdpos : 0 < d) (hd : d < divisor.size)
    (hquotient : quotient[q]? = some quotientTerm)
    (hdivisor : divisor[d]? = some divisorTerm)
    (hdegree : quotientTerm.1.deg + divisorTerm.1.deg = frontier.degree)
    (hselect : pairVecDivVHCSelectFrontier dividendIndex dividend heap nodes =
      .ok frontier) :
    ∃ (slot head : Nat) (node : PairVecDivVHCNode),
      heap[slot]? = some head ∧ nodes[d - 1]? = some node ∧
        d - 1 ∈ owners head ∧ node.divisorIndex = d ∧
        q = node.quotientIndex := by
  rcases hfixed.node_for_tail nodes divisor.size d hsize hdpos hd with
    ⟨node, hnode, hnodeD⟩
  have hi : d - 1 < nodes.size := by
    rw [hsize]
    omega
  have hcovered := hstate.covered_with heap nodes #[] resetH owners
    hownership
  rcases hcovered (d - 1) hi with hreset | hlin | howned
  · have habove := pairVecDivVHCResetRow_product_degree_gt_frontier
      degreeLimit frontier.degree resetH d q nodes quotient divisor
      quotientTerm divisorTerm hsize hfixed hresetReady hprefix hfrontier
      hdpos hd hreset hquotient hdivisor
    omega
  · simp at hlin
  · rcases howned with ⟨slot, head, hheap, hmem⟩
    rcases pairVecDivVHCChainOwns_mem_active (some head) (owners head) nodes
        (hownership.1 slot head hheap) (d - 1) hmem with
      ⟨activeNode, mono, hactiveGet, hactiveMono⟩
    rw [hnode] at hactiveGet
    simp only [Option.some.injEq] at hactiveGet
    subst activeNode
    have hnodeDenotes := hdenotes (d - 1) node hnode (by
      simp [hactiveMono])
    have hcursor := pairVecDivVHCOwnedRow_productAtFrontier_eq_cursor p
      degreeLimit dividendIndex dividend heap nodes owners frontier quotient
      divisor slot head (d - 1) q node quotientTerm divisorTerm hownership
      hhomogeneous hordered hcanonical hprefix hfrontier hheap hmem
      hnode hnodeDenotes hquotient (by simpa [hnodeD] using hdivisor) hdegree
      hselect
    exact ⟨slot, head, node, hheap, hnode, hmem, hnodeD, hcursor⟩

/-- Every concrete quotient/divisor-tail pair at the selected degree is not
merely represented by some cursor: its unique row node belongs to the exact
finite set consumed by the equal-degree heap loop.  Keeping the membership
statement at the concrete node index preserves multiplicity when distinct
rows happen to carry equal coefficient values. -/
theorem pairVecDivVHCFrontierProduct_mem_ownedNodesAtDegree
    (p degreeLimit dividendIndex resetH d q : Nat)
    (dividend quotient divisor : SparsePolyZp)
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode)
    (owners : Nat → Finset Nat) (frontier : PairVecDivVHCFrontier)
    (quotientTerm divisorTerm : UMonomial × Zp)
    (hsize : nodes.size = divisor.size - 1)
    (hfixed : PairVecDivVHCNodeDivisorIndicesFixed nodes)
    (hstate : PairVecDivVHCStateCovered heap nodes #[] resetH)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hhomogeneous : PairVecDivVHCHeapChainsHomogeneous heap owners nodes)
    (hresetReady : PairVecDivVHCResetReady resetH quotient.size nodes)
    (hordered : PairVecDivVHCHeapOrdered heap nodes)
    (hdenotes : ∀ (i : Nat) (node : PairVecDivVHCNode),
      nodes[i]? = some node → node.mono ≠ none →
        PairVecDivVHCNodeDenotes quotient divisor node)
    (hcanonical : SparsePolyZp.Canonical p quotient)
    (hprefix : PairVecDivVHCCursorPrefixAbove degreeLimit nodes quotient divisor)
    (hfrontier : frontier.degree < degreeLimit)
    (hdpos : 0 < d) (hd : d < divisor.size)
    (hquotient : quotient[q]? = some quotientTerm)
    (hdivisor : divisor[d]? = some divisorTerm)
    (hdegree : quotientTerm.1.deg + divisorTerm.1.deg = frontier.degree)
    (hselect : pairVecDivVHCSelectFrontier dividendIndex dividend heap nodes =
      .ok frontier) :
    d - 1 ∈ PairVecDivVHCHeapOwnedNodesAtDegree frontier.degree heap owners
      nodes := by
  rcases pairVecDivVHCFrontierProduct_owned_cursor p degreeLimit dividendIndex
      resetH d q dividend quotient divisor heap nodes owners frontier
      quotientTerm divisorTerm hsize hfixed hstate hownership hhomogeneous
      hresetReady hordered hdenotes hcanonical hprefix hfrontier hdpos hd
      hquotient hdivisor hdegree hselect with
    ⟨slot, head, node, hheap, hnode, hmem, hnodeD, hcursor⟩
  have howns := hownership.1 slot head hheap
  rcases pairVecDivVHCChainOwns_mem_active (some head) (owners head) nodes
      howns (d - 1) hmem with
    ⟨activeNode, mono, hactiveGet, hactiveMono⟩
  rw [hnode] at hactiveGet
  simp only [Option.some.injEq] at hactiveGet
  subst activeNode
  rcases hdenotes (d - 1) node hnode (by simp [hactiveMono]) with
    ⟨storedQuotient, storedDivisor, hstoredQuotient, hstoredDivisor,
      hstoredMono⟩
  rw [← hcursor, hquotient] at hstoredQuotient
  rw [hnodeD, hdivisor] at hstoredDivisor
  simp only [Option.some.injEq] at hstoredQuotient hstoredDivisor
  subst storedQuotient
  subst storedDivisor
  rw [hactiveMono] at hstoredMono
  have hmonoDegree := congrArg UMonomial.deg (Option.some.inj hstoredMono)
  simp only [PairVecDivVHCHeapOwnedNodesAtDegree, Finset.mem_filter,
    PairVecDivVHCHeapOwnedNodes, Finset.mem_biUnion, List.mem_toFinset,
    pairVecDivVHCNodeAtDegree_iff]
  refine ⟨⟨head, ?_, hmem⟩, mono, ?_, ?_⟩
  · have hslot : slot < heap.size := by
      by_contra hnot
      rw [Array.getElem?_eq_none (by omega)] at hheap
      contradiction
    have hheadEq : heap[slot] = head := by
      rw [Array.getElem?_eq_getElem hslot] at hheap
      exact Option.some.inj hheap
    rw [← hheadEq]
    exact Array.getElem_mem_toList hslot
  · exact (pairVecDivVHCMono_eq_ok_iff (d - 1) nodes mono).2
      ⟨node, hnode, hactiveMono⟩
  · exact hmonoDegree.trans hdegree

theorem pairVecDivVHCHeapChainsHomogeneous_of_nextValid
    (heap : Array Nat) (owners : Nat → Finset Nat)
    (nodes : Array PairVecDivVHCNode)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hnextValid : PairVecDivVHCNextValid nodes) :
    PairVecDivVHCHeapChainsHomogeneous heap owners nodes := by
  intro slot head mono hheap hmono
  have howns := hownership.1 slot head hheap
  have hvalid := pairVecDivVHCChainOwns_valid (some head) (owners head) nodes
    howns
  unfold pairVecDivVHCMono at hmono
  split at hmono <;> try contradiction
  next hn =>
    split at hmono <;> try contradiction
    next storedMono hstoredMono =>
      simp only [Except.ok.injEq] at hmono
      subst storedMono
      exact pairVecDivVHCChainValid_atDegree_of_nextValid head (owners head)
        nodes mono.deg hvalid hnextValid
          ⟨nodes[head], mono, Array.getElem?_eq_getElem hn, hstoredMono, rfl⟩

theorem pairVecDivVHCFrontierBelow_of_remaining_owned
    (degreeLimit dividendIndex : Nat) (dividend : SparsePolyZp)
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode)
    (owners : Nat → Finset Nat)
    (hremaining : PairVecDivVHCRemainingDividendBelow degreeLimit
      dividendIndex dividend)
    (hbelow : PairVecDivVHCAllActiveNodesBelow degreeLimit nodes)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes) :
    PairVecDivVHCFrontierBelow degreeLimit dividendIndex dividend heap nodes := by
  refine ⟨hremaining, ?_⟩
  intro slot nodeIndex node mono hheap hnode hmono
  have howns := hownership.1 slot nodeIndex hheap
  have hheadMem := pairVecDivVHCChainOwns_head_mem nodeIndex
    (owners nodeIndex) nodes howns
  rcases pairVecDivVHCChainOwns_mem_active (some nodeIndex)
      (owners nodeIndex) nodes howns nodeIndex hheadMem with
    ⟨ownedNode, ownedMono, hownedNode, hownedMono⟩
  rw [hnode] at hownedNode
  simp only [Option.some.injEq] at hownedNode
  subst ownedNode
  rw [hmono] at hownedMono
  simp only [Option.some.injEq] at hownedMono
  subst ownedMono
  exact hbelow nodeIndex node mono hnode hmono

/-- Topology-independent ownership invariant used across insertion: the owner
map may change when a newly activated head is linked in front of an existing
equal-monomial bucket. -/
def PairVecDivVHCHeapChainsOwned (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode) : Prop :=
  ∃ owners : Nat → Finset Nat,
    PairVecDivVHCHeapChainOwnership heap owners nodes

/-- Heap ownership together with a protected set (used for the source `lin`
stack) that is absent from every active bucket and from every heap head. -/
def PairVecDivVHCHeapChainsOwnedAway (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode) (protectedSet : Finset Nat) : Prop :=
  ∃ owners : Nat → Finset Nat,
    PairVecDivVHCHeapChainOwnership heap owners nodes ∧
      ∀ (slot head : Nat), heap[slot]? = some head →
        Disjoint protectedSet (owners head) ∧ head ∉ protectedSet

theorem pairVecDivVHCHeapChainsOwned_empty (nodes : Array PairVecDivVHCNode) :
    PairVecDivVHCHeapChainsOwned #[] nodes := by
  refine ⟨fun _ => ∅, ?_⟩
  refine ⟨?_, ?_, ?_⟩
  · intro slot head hget
    simp at hget
  · intro left right head hleft
    simp at hleft
  · intro left right leftHead rightHead hleft
    simp at hleft

theorem pairVecDivVHCHeapChainsOwnedAway_empty
    (nodes : Array PairVecDivVHCNode) (protectedSet : Finset Nat) :
    PairVecDivVHCHeapChainsOwnedAway #[] nodes protectedSet := by
  rcases pairVecDivVHCHeapChainsOwned_empty nodes with ⟨owners, hownership⟩
  refine ⟨owners, hownership, ?_⟩
  intro slot head hget
  simp at hget

theorem PairVecDivVHCHeapChainsOwned.away_empty
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode)
    (howned : PairVecDivVHCHeapChainsOwned heap nodes) :
    PairVecDivVHCHeapChainsOwnedAway heap nodes ∅ := by
  rcases howned with ⟨owners, hownership⟩
  refine ⟨owners, hownership, ?_⟩
  intro slot head hget
  simp

theorem pairVecDivVHCHeapChainsOwnedAway.owned
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode)
    (protectedSet : Finset Nat)
    (haway : PairVecDivVHCHeapChainsOwnedAway heap nodes protectedSet) :
    PairVecDivVHCHeapChainsOwned heap nodes := by
  rcases haway with ⟨owners, hownership, hprotected⟩
  exact ⟨owners, hownership⟩

theorem pairVecDivVHCHeapChainsOwnedAway.mono
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode)
    (small large : Finset Nat)
    (haway : PairVecDivVHCHeapChainsOwnedAway heap nodes large)
    (hsubset : small ⊆ large) :
    PairVecDivVHCHeapChainsOwnedAway heap nodes small := by
  rcases haway with ⟨owners, hownership, hprotected⟩
  refine ⟨owners, hownership, ?_⟩
  intro slot head hget
  have hold := hprotected slot head hget
  exact ⟨Finset.disjoint_left.mpr (by
    intro i hiSmall hiOwner
    exact Finset.disjoint_left.mp hold.1 (hsubset hiSmall) hiOwner),
    fun hiSmall => hold.2 (hsubset hiSmall)⟩

/-- Ownership is invariant under a head-array reordering once provenance and
target uniqueness have been proved for the concrete pointer-copy routine. -/
theorem pairVecDivVHCHeapChainOwnership_of_valuesFrom
    (source target : Array Nat) (owners : Nat → Finset Nat)
    (nodes : Array PairVecDivVHCNode)
    (hownership : PairVecDivVHCHeapChainOwnership source owners nodes)
    (hfrom : PairVecDivVHCValuesFrom target source)
    (hunique : ∀ (left right head : Nat), target[left]? = some head →
      target[right]? = some head → left = right) :
    PairVecDivVHCHeapChainOwnership target owners nodes := by
  refine ⟨?_, hunique, ?_⟩
  · intro slot head hget
    obtain ⟨sourceSlot, hsource⟩ := hfrom slot head hget
    exact hownership.1 sourceSlot head hsource
  · intro left right leftHead rightHead hleft hright hne
    obtain ⟨sourceLeft, hsourceLeft⟩ := hfrom left leftHead hleft
    obtain ⟨sourceRight, hsourceRight⟩ := hfrom right rightHead hright
    exact hownership.2.2 sourceLeft sourceRight leftHead rightHead hsourceLeft
      hsourceRight hne

theorem pairVecDivVHCHeapChainsOwnedAway_of_valuesFrom
    (source target : Array Nat) (owners : Nat → Finset Nat)
    (nodes : Array PairVecDivVHCNode) (protectedSet : Finset Nat)
    (hownership : PairVecDivVHCHeapChainOwnership target owners nodes)
    (hsourceAway : ∀ (slot head : Nat), source[slot]? = some head →
      Disjoint protectedSet (owners head) ∧ head ∉ protectedSet)
    (hfrom : PairVecDivVHCValuesFrom target source) :
    PairVecDivVHCHeapChainsOwnedAway target nodes protectedSet := by
  refine ⟨owners, hownership, ?_⟩
  intro slot head hget
  obtain ⟨sourceSlot, hsource⟩ := hfrom slot head hget
  exact hsourceAway sourceSlot head hsource

theorem pairVecDivVHCHeapChainOwnership_push_fresh
    (heap : Array Nat) (owners : Nat → Finset Nat)
    (newNode : Nat) (nodes nodes' : Array PairVecDivVHCNode)
    (node : PairVecDivVHCNode) (mono : UMonomial)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hnode : nodes[newNode]? = some node) (hmono : node.mono = some mono)
    (hfreshHead : ∀ (slot : Nat), heap[slot]? ≠ some newNode)
    (hfreshOwners : ∀ (slot head : Nat), heap[slot]? = some head →
      newNode ∉ owners head)
    (hset : pairVecDivVHCSetNext newNode none nodes = .ok nodes') :
    PairVecDivVHCHeapChainOwnership (heap.push newNode)
      (fun head => if head = newNode then {newNode} else owners head) nodes' := by
  let owners' := fun head => if head = newNode then {newNode} else owners head
  have hnewOwns : PairVecDivVHCChainOwns (some newNode) {newNode} nodes' := by
    simpa using pairVecDivVHCSetNext_chainOwns_insert newNode none nodes nodes'
      node mono ∅ hnode hmono (by simp) (by simp [PairVecDivVHCChainOwns]) hset
  refine ⟨?_, pairVecDivVHCPush_preserves_unique heap newNode
    hownership.2.1 hfreshHead, ?_⟩
  · intro slot head hget
    rw [Array.getElem?_push] at hget
    by_cases hlast : slot = heap.size
    · subst slot
      simp only [ite_true, Option.some.injEq] at hget
      subst head
      simpa [owners'] using hnewOwns
    · simp only [hlast, ↓reduceIte] at hget
      have hne : head ≠ newNode := by
        intro heq
        subst head
        exact hfreshHead slot hget
      have hold := hownership.1 slot head hget
      have hold' := pairVecDivVHCSetNext_preserves_disjoint_chain newNode none
        nodes nodes' (some head) (owners head) (hfreshOwners slot head hget)
        hold hset
      simpa [owners', hne] using hold'
  · intro left right leftHead rightHead hleft hright hne
    rw [Array.getElem?_push] at hleft hright
    by_cases hleftNew : leftHead = newNode
    · subst leftHead
      have hrightNew : rightHead ≠ newNode := by
        exact fun heq => hne heq.symm
      have hrightOld : right ≠ heap.size := by
        intro heq
        subst right
        simp only [ite_true, Option.some.injEq] at hright
        exact hrightNew hright.symm
      simp only [hrightOld, ↓reduceIte] at hright
      simp only [owners', if_pos rfl, hrightNew, ↓reduceIte]
      exact Finset.disjoint_left.mpr (by
        intro x hx hxo
        simp only [Finset.mem_singleton] at hx
        subst x
        exact hfreshOwners right rightHead hright hxo)
    · by_cases hrightNew : rightHead = newNode
      · subst rightHead
        have hleftOld : left ≠ heap.size := by
          intro heq
          subst left
          simp only [ite_true, Option.some.injEq] at hleft
          exact hleftNew hleft.symm
        simp only [hleftOld, ↓reduceIte] at hleft
        simp only [owners', hleftNew, ↓reduceIte, if_pos rfl]
        exact Finset.disjoint_left.mpr (by
          intro x hxo hx
          simp only [Finset.mem_singleton] at hx
          subst x
          exact hfreshOwners left leftHead hleft hxo)
      · have hleftOld : left ≠ heap.size := by
          intro heq
          subst left
          simp only [ite_true, Option.some.injEq] at hleft
          exact hleftNew hleft.symm
        have hrightOld : right ≠ heap.size := by
          intro heq
          subst right
          simp only [ite_true, Option.some.injEq] at hright
          exact hrightNew hright.symm
        simp only [hleftOld, ↓reduceIte] at hleft
        simp only [hrightOld, ↓reduceIte] at hright
        simpa [owners', hleftNew, hrightNew] using
          hownership.2.2 left right leftHead rightHead hleft hright hne

theorem pairVecDivVHCHeapChainsOwnedAway_push_fresh
    (heap : Array Nat) (owners : Nat → Finset Nat)
    (newNode : Nat) (nodes nodes' : Array PairVecDivVHCNode)
    (node : PairVecDivVHCNode) (mono : UMonomial)
    (protectedSet remaining : Finset Nat)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (haway : ∀ (slot head : Nat), heap[slot]? = some head →
      Disjoint protectedSet (owners head) ∧ head ∉ protectedSet)
    (hremaining : remaining ⊆ protectedSet) (hnewRemaining : newNode ∉ remaining)
    (hnode : nodes[newNode]? = some node) (hmono : node.mono = some mono)
    (hfreshHead : ∀ (slot : Nat), heap[slot]? ≠ some newNode)
    (hfreshOwners : ∀ (slot head : Nat), heap[slot]? = some head →
      newNode ∉ owners head)
    (hset : pairVecDivVHCSetNext newNode none nodes = .ok nodes') :
    PairVecDivVHCHeapChainsOwnedAway (heap.push newNode) nodes' remaining := by
  let owners' := fun head => if head = newNode then {newNode} else owners head
  have hownership' := pairVecDivVHCHeapChainOwnership_push_fresh heap owners
    newNode nodes nodes' node mono hownership hnode hmono hfreshHead
    hfreshOwners hset
  refine ⟨owners', hownership', ?_⟩
  intro slot head hget
  rw [Array.getElem?_push] at hget
  by_cases hlast : slot = heap.size
  · subst slot
    simp only [ite_true, Option.some.injEq] at hget
    subst head
    simp only [owners', if_pos rfl]
    exact ⟨Finset.disjoint_left.mpr (by
      intro i hi hsingle
      simp only [Finset.mem_singleton] at hsingle
      subst i
      exact hnewRemaining hi), hnewRemaining⟩
  · simp only [hlast, ↓reduceIte] at hget
    have hne : head ≠ newNode := by
      intro heq
      subst head
      exact hfreshHead slot hget
    simp only [owners', hne, ↓reduceIte]
    have hold := haway slot head hget
    exact ⟨Finset.disjoint_left.mpr (by
      intro i hi howner
      exact Finset.disjoint_left.mp hold.1 (hremaining hi) howner),
      fun hi => hold.2 (hremaining hi)⟩

theorem pairVecDivVHCPush_fresh_separated
    (heap : Array Nat) (owners : Nat → Finset Nat) (newNode : Nat)
    (protectedSet remaining : Finset Nat)
    (haway : ∀ (slot head : Nat), heap[slot]? = some head →
      Disjoint protectedSet (owners head) ∧ head ∉ protectedSet)
    (hremaining : remaining ⊆ protectedSet) (hnewRemaining : newNode ∉ remaining)
    (hfreshHead : ∀ (slot : Nat), heap[slot]? ≠ some newNode) :
    ∀ (slot head : Nat), (heap.push newNode)[slot]? = some head →
      Disjoint remaining
        ((fun h => if h = newNode then {newNode} else owners h) head) ∧
        head ∉ remaining := by
  intro slot head hget
  rw [Array.getElem?_push] at hget
  by_cases hlast : slot = heap.size
  · subst slot
    simp only [ite_true, Option.some.injEq] at hget
    subst head
    exact ⟨Finset.disjoint_left.mpr (by
      intro i hi hsingle
      have hiNew : i = newNode := by simpa using hsingle
      subst i
      exact hnewRemaining hi), hnewRemaining⟩
  · simp only [hlast, ↓reduceIte] at hget
    have hne : head ≠ newNode := by
      intro heq
      subst head
      exact hfreshHead slot hget
    simp only [hne, ↓reduceIte]
    have hold := haway slot head hget
    exact ⟨Finset.disjoint_left.mpr (by
      intro i hi howner
      exact Finset.disjoint_left.mp hold.1 (hremaining hi) howner),
      fun hi => hold.2 (hremaining hi)⟩

theorem pairVecDivVHCHeapChainOwnership_bubble_fresh
    (heap heap' : Array Nat) (owners : Nat → Finset Nat)
    (stop newNode : Nat) (nodes nodes' : Array PairVecDivVHCNode)
    (node : PairVecDivVHCNode) (mono : UMonomial)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hnode : nodes[newNode]? = some node) (hmono : node.mono = some mono)
    (hfreshHead : ∀ (slot : Nat), heap[slot]? ≠ some newNode)
    (hfreshOwners : ∀ (slot head : Nat), heap[slot]? = some head →
      newNode ∉ owners head)
    (hset : pairVecDivVHCSetNext newNode none nodes = .ok nodes')
    (hbubble : pairVecDivVHCBubble heap.size stop newNode
      (heap.push newNode) = .ok heap') :
    PairVecDivVHCHeapChainOwnership heap'
      (fun head => if head = newNode then {newNode} else owners head) nodes' := by
  let owners' := fun head => if head = newNode then {newNode} else owners head
  have hpush := pairVecDivVHCHeapChainOwnership_push_fresh heap owners newNode
    nodes nodes' node mono hownership hnode hmono hfreshHead hfreshOwners hset
  have hfrom := pairVecDivVHCBubble_valuesFrom heap.size stop newNode
    (heap.push newNode) heap' (heap.push newNode)
    (pairVecDivVHCValuesFrom_refl _) ⟨heap.size, by simp⟩ hbubble
  exact pairVecDivVHCHeapChainOwnership_of_valuesFrom (heap.push newNode) heap'
    owners' nodes' hpush hfrom
    (pairVecDivVHCBubble_push_preserves_unique stop newNode heap heap'
      hownership.2.1 hfreshHead hbubble)

theorem pairVecDivVHCHeapChainOwnership_bubbleBelow_fresh
    (heap heap' : Array Nat) (owners : Nat → Finset Nat)
    (anchor newNode : Nat) (nodes nodes' : Array PairVecDivVHCNode)
    (node : PairVecDivVHCNode) (mono : UMonomial)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hnode : nodes[newNode]? = some node) (hmono : node.mono = some mono)
    (hfreshHead : ∀ (slot : Nat), heap[slot]? ≠ some newNode)
    (hfreshOwners : ∀ (slot head : Nat), heap[slot]? = some head →
      newNode ∉ owners head)
    (hset : pairVecDivVHCSetNext newNode none nodes = .ok nodes')
    (hbubble : pairVecDivVHCBubbleBelow heap.size anchor newNode
      (heap.push newNode) = .ok heap') :
    PairVecDivVHCHeapChainOwnership heap'
      (fun head => if head = newNode then {newNode} else owners head) nodes' := by
  let owners' := fun head => if head = newNode then {newNode} else owners head
  have hpush := pairVecDivVHCHeapChainOwnership_push_fresh heap owners newNode
    nodes nodes' node mono hownership hnode hmono hfreshHead hfreshOwners hset
  have hfrom := pairVecDivVHCBubbleBelow_valuesFrom heap.size anchor newNode
    (heap.push newNode) heap' (heap.push newNode)
    (pairVecDivVHCValuesFrom_refl _) ⟨heap.size, by simp⟩ hbubble
  exact pairVecDivVHCHeapChainOwnership_of_valuesFrom (heap.push newNode) heap'
    owners' nodes' hpush hfrom
    (pairVecDivVHCBubbleBelow_push_preserves_unique anchor newNode heap heap'
      hownership.2.1 hfreshHead hbubble)

theorem pairVecDivVHCHeapChainsOwnedAway_bubble_fresh
    (heap heap' : Array Nat) (owners : Nat → Finset Nat)
    (stop newNode : Nat) (nodes nodes' : Array PairVecDivVHCNode)
    (node : PairVecDivVHCNode) (mono : UMonomial)
    (protectedSet remaining : Finset Nat)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (haway : ∀ (slot head : Nat), heap[slot]? = some head →
      Disjoint protectedSet (owners head) ∧ head ∉ protectedSet)
    (hremaining : remaining ⊆ protectedSet) (hnewRemaining : newNode ∉ remaining)
    (hnode : nodes[newNode]? = some node) (hmono : node.mono = some mono)
    (hfreshHead : ∀ (slot : Nat), heap[slot]? ≠ some newNode)
    (hfreshOwners : ∀ (slot head : Nat), heap[slot]? = some head →
      newNode ∉ owners head)
    (hset : pairVecDivVHCSetNext newNode none nodes = .ok nodes')
    (hbubble : pairVecDivVHCBubble heap.size stop newNode
      (heap.push newNode) = .ok heap') :
    PairVecDivVHCHeapChainsOwnedAway heap' nodes' remaining := by
  let owners' := fun head => if head = newNode then {newNode} else owners head
  have hpushSeparated := pairVecDivVHCPush_fresh_separated heap owners newNode
    protectedSet remaining haway hremaining hnewRemaining hfreshHead
  have htargetOwnership := pairVecDivVHCHeapChainOwnership_bubble_fresh heap
    heap' owners stop newNode nodes nodes' node mono hownership hnode hmono
    hfreshHead hfreshOwners hset hbubble
  have hfrom := pairVecDivVHCBubble_valuesFrom heap.size stop newNode
    (heap.push newNode) heap' (heap.push newNode)
    (pairVecDivVHCValuesFrom_refl _) ⟨heap.size, by simp⟩ hbubble
  exact pairVecDivVHCHeapChainsOwnedAway_of_valuesFrom (heap.push newNode)
    heap' owners' nodes' remaining htargetOwnership hpushSeparated hfrom

theorem pairVecDivVHCHeapChainsOwnedAway_bubbleBelow_fresh
    (heap heap' : Array Nat) (owners : Nat → Finset Nat)
    (anchor newNode : Nat) (nodes nodes' : Array PairVecDivVHCNode)
    (node : PairVecDivVHCNode) (mono : UMonomial)
    (protectedSet remaining : Finset Nat)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (haway : ∀ (slot head : Nat), heap[slot]? = some head →
      Disjoint protectedSet (owners head) ∧ head ∉ protectedSet)
    (hremaining : remaining ⊆ protectedSet) (hnewRemaining : newNode ∉ remaining)
    (hnode : nodes[newNode]? = some node) (hmono : node.mono = some mono)
    (hfreshHead : ∀ (slot : Nat), heap[slot]? ≠ some newNode)
    (hfreshOwners : ∀ (slot head : Nat), heap[slot]? = some head →
      newNode ∉ owners head)
    (hset : pairVecDivVHCSetNext newNode none nodes = .ok nodes')
    (hbubble : pairVecDivVHCBubbleBelow heap.size anchor newNode
      (heap.push newNode) = .ok heap') :
    PairVecDivVHCHeapChainsOwnedAway heap' nodes' remaining := by
  let owners' := fun head => if head = newNode then {newNode} else owners head
  have hpushSeparated := pairVecDivVHCPush_fresh_separated heap owners newNode
    protectedSet remaining haway hremaining hnewRemaining hfreshHead
  have htargetOwnership :=
    pairVecDivVHCHeapChainOwnership_bubbleBelow_fresh heap heap' owners anchor
      newNode nodes nodes' node mono hownership hnode hmono hfreshHead
      hfreshOwners hset hbubble
  have hfrom := pairVecDivVHCBubbleBelow_valuesFrom heap.size anchor newNode
    (heap.push newNode) heap' (heap.push newNode)
    (pairVecDivVHCValuesFrom_refl _) ⟨heap.size, by simp⟩ hbubble
  exact pairVecDivVHCHeapChainsOwnedAway_of_valuesFrom (heap.push newNode)
    heap' owners' nodes' remaining htargetOwnership hpushSeparated hfrom

theorem pairVecDivVHCHeapChainOwnership_merge_fresh
    (heap : Array Nat) (owners : Nat → Finset Nat) (slot newNode : Nat)
    (nodes nodes' : Array PairVecDivVHCNode) (node : PairVecDivVHCNode)
    (mono : UMonomial) (hslot : slot < heap.size)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hnode : nodes[newNode]? = some node) (hmono : node.mono = some mono)
    (hfreshHead : ∀ (i : Nat), heap[i]? ≠ some newNode)
    (hfreshOwners : ∀ (i head : Nat), heap[i]? = some head →
      newNode ∉ owners head)
    (hset : pairVecDivVHCSetNext newNode (some heap[slot]) nodes = .ok nodes') :
    PairVecDivVHCHeapChainOwnership (heap.set slot newNode)
      (fun head => if head = newNode then insert newNode (owners heap[slot])
        else owners head) nodes' := by
  let owners' := fun head => if head = newNode then
    insert newNode (owners heap[slot]) else owners head
  have holdGet : heap[slot]? = some heap[slot] := by
    rw [Array.getElem?_eq_getElem hslot]
  have holdOwns := hownership.1 slot heap[slot] holdGet
  have hfreshOld := hfreshOwners slot heap[slot] holdGet
  have hnewOwns := pairVecDivVHCSetNext_chainOwns_insert newNode
    (some heap[slot]) nodes nodes' node mono (owners heap[slot]) hnode hmono
    hfreshOld holdOwns hset
  refine ⟨?_, pairVecDivVHCSet_fresh_preserves_unique heap slot newNode hslot
    hownership.2.1 hfreshHead, ?_⟩
  · intro targetSlot head hget
    by_cases hat : slot = targetSlot
    · subst targetSlot
      rw [Array.getElem?_set_self hslot] at hget
      simp only [Option.some.injEq] at hget
      subst head
      simpa [owners'] using hnewOwns
    · rw [Array.getElem?_set_ne hslot hat] at hget
      have hne : head ≠ newNode := by
        intro heq
        subst head
        exact hfreshHead targetSlot hget
      have hold := hownership.1 targetSlot head hget
      have hold' := pairVecDivVHCSetNext_preserves_disjoint_chain newNode
        (some heap[slot]) nodes nodes' (some head) (owners head)
        (hfreshOwners targetSlot head hget) hold hset
      simpa [owners', hne] using hold'
  · intro left right leftHead rightHead hleft hright hne
    by_cases hleftNew : leftHead = newNode
    · subst leftHead
      have hrightNew : rightHead ≠ newNode := fun heq => hne heq.symm
      have hrightNotSlot : slot ≠ right := by
        intro heq
        subst right
        rw [Array.getElem?_set_self hslot] at hright
        exact hrightNew (Option.some.inj hright).symm
      rw [Array.getElem?_set_ne hslot hrightNotSlot] at hright
      have hrightOldNe : rightHead ≠ heap[slot] := by
        intro heq
        subst rightHead
        have := hownership.2.1 right slot heap[slot] hright holdGet
        exact hrightNotSlot this.symm
      simp only [owners', if_pos rfl, hrightNew, ↓reduceIte]
      exact Finset.disjoint_left.mpr (by
        intro x hx hxo
        rw [Finset.mem_insert] at hx
        rcases hx with rfl | hx
        · exact hfreshOwners right rightHead hright hxo
        · exact Finset.disjoint_left.mp
            (hownership.2.2 slot right heap[slot] rightHead holdGet hright
              hrightOldNe.symm) hx hxo)
    · by_cases hrightNew : rightHead = newNode
      · subst rightHead
        have hleftNotSlot : slot ≠ left := by
          intro heq
          subst left
          rw [Array.getElem?_set_self hslot] at hleft
          exact hleftNew (Option.some.inj hleft).symm
        rw [Array.getElem?_set_ne hslot hleftNotSlot] at hleft
        have hleftOldNe : leftHead ≠ heap[slot] := by
          intro heq
          subst leftHead
          have := hownership.2.1 left slot heap[slot] hleft holdGet
          exact hleftNotSlot this.symm
        simp only [owners', hleftNew, ↓reduceIte, if_pos rfl]
        exact Finset.disjoint_left.mpr (by
          intro x hxo hx
          rw [Finset.mem_insert] at hx
          rcases hx with rfl | hx
          · exact hfreshOwners left leftHead hleft hxo
          · exact Finset.disjoint_left.mp
              (hownership.2.2 left slot leftHead heap[slot] hleft holdGet
                hleftOldNe) hxo hx)
      · have hleftNotSlot : slot ≠ left := by
          intro heq
          subst left
          rw [Array.getElem?_set_self hslot] at hleft
          exact hleftNew (Option.some.inj hleft).symm
        have hrightNotSlot : slot ≠ right := by
          intro heq
          subst right
          rw [Array.getElem?_set_self hslot] at hright
          exact hrightNew (Option.some.inj hright).symm
        rw [Array.getElem?_set_ne hslot hleftNotSlot] at hleft
        rw [Array.getElem?_set_ne hslot hrightNotSlot] at hright
        simpa [owners', hleftNew, hrightNew] using
          hownership.2.2 left right leftHead rightHead hleft hright hne

theorem pairVecDivVHCSet_fresh_separated
    (heap : Array Nat) (owners : Nat → Finset Nat) (slot newNode : Nat)
    (protectedSet remaining : Finset Nat) (hslot : slot < heap.size)
    (haway : ∀ (i head : Nat), heap[i]? = some head →
      Disjoint protectedSet (owners head) ∧ head ∉ protectedSet)
    (hremaining : remaining ⊆ protectedSet)
    (hnewRemaining : newNode ∉ remaining)
    (hfreshHead : ∀ (i : Nat), heap[i]? ≠ some newNode) :
    ∀ (targetSlot head : Nat), (heap.set slot newNode)[targetSlot]? = some head →
      Disjoint remaining
        ((fun h => if h = newNode then insert newNode (owners heap[slot])
          else owners h) head) ∧
        head ∉ remaining := by
  have holdGet : heap[slot]? = some heap[slot] := by
    rw [Array.getElem?_eq_getElem hslot]
  intro targetSlot head hget
  by_cases hat : slot = targetSlot
  · subst targetSlot
    rw [Array.getElem?_set_self hslot] at hget
    simp only [Option.some.injEq] at hget
    subst head
    refine ⟨Finset.disjoint_left.mpr ?_, hnewRemaining⟩
    intro i hi howner
    simp at howner
    rcases howner with rfl | howner
    · exact hnewRemaining hi
    · have hold := haway slot heap[slot] holdGet
      exact Finset.disjoint_left.mp hold.1 (hremaining hi) howner
  · rw [Array.getElem?_set_ne hslot hat] at hget
    have hne : head ≠ newNode := by
      intro heq
      subst head
      exact hfreshHead targetSlot hget
    simp only [hne, ↓reduceIte]
    have hold := haway targetSlot head hget
    exact ⟨Finset.disjoint_left.mpr (by
      intro i hi howner
      exact Finset.disjoint_left.mp hold.1 (hremaining hi) howner),
      fun hi => hold.2 (hremaining hi)⟩

theorem pairVecDivVHCHeapChainsOwnedAway_merge_fresh
    (heap : Array Nat) (owners : Nat → Finset Nat) (slot newNode : Nat)
    (nodes nodes' : Array PairVecDivVHCNode) (node : PairVecDivVHCNode)
    (mono : UMonomial) (protectedSet remaining : Finset Nat)
    (hslot : slot < heap.size)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (haway : ∀ (i head : Nat), heap[i]? = some head →
      Disjoint protectedSet (owners head) ∧ head ∉ protectedSet)
    (hremaining : remaining ⊆ protectedSet)
    (hnewRemaining : newNode ∉ remaining)
    (hnode : nodes[newNode]? = some node) (hmono : node.mono = some mono)
    (hfreshHead : ∀ (i : Nat), heap[i]? ≠ some newNode)
    (hfreshOwners : ∀ (i head : Nat), heap[i]? = some head →
      newNode ∉ owners head)
    (hset : pairVecDivVHCSetNext newNode (some heap[slot]) nodes = .ok nodes') :
    PairVecDivVHCHeapChainsOwnedAway (heap.set slot newNode) nodes' remaining := by
  let owners' := fun head => if head = newNode then
    insert newNode (owners heap[slot]) else owners head
  have hownership' := pairVecDivVHCHeapChainOwnership_merge_fresh heap owners
    slot newNode nodes nodes' node mono hslot hownership hnode hmono
    hfreshHead hfreshOwners hset
  exact ⟨owners', hownership', pairVecDivVHCSet_fresh_separated heap owners
    slot newNode protectedSet remaining hslot haway hremaining hnewRemaining
    hfreshHead⟩

theorem pairVecDivVHCInsert_preserves_heapChainOwnership_of_fresh
    (newNode : Nat) (heap heap' : Array Nat) (nodes nodes' : Array PairVecDivVHCNode)
    (owners : Nat → Finset Nat) (node : PairVecDivVHCNode) (mono : UMonomial)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hnode : nodes[newNode]? = some node) (hmono : node.mono = some mono)
    (hfreshHead : ∀ (slot : Nat), heap[slot]? ≠ some newNode)
    (hfreshOwners : ∀ (slot head : Nat), heap[slot]? = some head →
      newNode ∉ owners head)
    (hrun : pairVecDivVHCInsert newNode heap nodes = .ok (heap', nodes')) :
    PairVecDivVHCHeapChainsOwned heap' nodes' := by
  cases hnew : pairVecDivVHCMono newNode nodes with
  | error fault => simp [pairVecDivVHCInsert, hnew] at hrun
  | ok newMono =>
      by_cases hempty : heap.size = 0
      · have heq : heap = #[] := Array.eq_empty_of_size_eq_zero hempty
        subst heap
        unfold pairVecDivVHCInsert at hrun
        simp only [hnew, Array.size_empty, ↓reduceDIte] at hrun
        cases hset : pairVecDivVHCSetNext newNode none nodes with
        | error fault => simp [hset] at hrun
        | ok updated =>
            rw [hset] at hrun
            simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
            rcases hrun with ⟨rfl, rfl⟩
            exact ⟨_, pairVecDivVHCHeapChainOwnership_push_fresh #[] owners
              newNode nodes updated node mono hownership hnode hmono
              (by simp) (by simp) hset⟩
      · have hheap : 0 < heap.size := Nat.pos_of_ne_zero hempty
        cases hroot : pairVecDivVHCMono heap[0] nodes with
        | error fault =>
            simp [pairVecDivVHCInsert, hnew, hempty, hroot] at hrun
        | ok rootMono =>
            by_cases hequal : newMono.deg = rootMono.deg
            · unfold pairVecDivVHCInsert at hrun
              simp only [hnew, hempty, ↓reduceDIte, hroot, hequal] at hrun
              cases hset : pairVecDivVHCSetNext newNode (some heap[0]) nodes with
              | error fault => simp [hset] at hrun
              | ok updated =>
                  rw [hset] at hrun
                  simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
                  rcases hrun with ⟨rfl, rfl⟩
                  exact ⟨_, pairVecDivVHCHeapChainOwnership_merge_fresh heap
                    owners 0 newNode nodes updated node mono hheap hownership
                    hnode hmono hfreshHead hfreshOwners hset⟩
            · by_cases hgreater : newMono.deg > rootMono.deg
              · unfold pairVecDivVHCInsert at hrun
                simp only [hnew, hempty, ↓reduceDIte, hroot, hequal,
                  hgreater] at hrun
                cases hset : pairVecDivVHCSetNext newNode none nodes with
                | error fault => simp [hset] at hrun
                | ok updated =>
                    rw [hset] at hrun
                    cases hbubble : pairVecDivVHCBubble heap.size 0 newNode
                        (heap.push newNode) with
                    | error fault => simp [hbubble] at hrun
                    | ok shifted =>
                        rw [hbubble] at hrun
                        simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
                        rcases hrun with ⟨rfl, rfl⟩
                        exact ⟨_,
                          pairVecDivVHCHeapChainOwnership_bubble_fresh heap
                            shifted owners 0 newNode nodes updated node mono
                            hownership hnode hmono hfreshHead hfreshOwners hset
                            hbubble⟩
              · cases hanchor : pairVecDivVHCFindAnchor newMono.deg
                    (pairVecDivVHCParent heap.size) heap nodes with
                | error fault =>
                    simp [pairVecDivVHCInsert, hnew, hempty, hroot, hequal,
                      hgreater, hanchor] at hrun
                | ok anchor =>
                    by_cases ha : anchor < heap.size
                    · cases hanchorMono : pairVecDivVHCMono heap[anchor] nodes with
                      | error fault =>
                          simp [pairVecDivVHCInsert, hnew, hempty, hroot,
                            hequal, hgreater, hanchor, ha, hanchorMono] at hrun
                      | ok anchorMono =>
                          by_cases hequalAnchor :
                              newMono.deg = anchorMono.deg
                          · unfold pairVecDivVHCInsert at hrun
                            simp only [hnew, hempty, ↓reduceDIte, hroot,
                              hequal, hgreater, hanchor, ha, hanchorMono] at hrun
                            simp only [hequalAnchor] at hrun
                            cases hset : pairVecDivVHCSetNext newNode
                                (some heap[anchor]) nodes with
                            | error fault => simp [hset] at hrun
                            | ok updated =>
                                rw [hset] at hrun
                                simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
                                rcases hrun with ⟨rfl, rfl⟩
                                exact ⟨_,
                                  pairVecDivVHCHeapChainOwnership_merge_fresh
                                    heap owners anchor newNode nodes nodes' node
                                    mono ha hownership hnode hmono hfreshHead
                                    hfreshOwners hset⟩
                          · unfold pairVecDivVHCInsert at hrun
                            simp only [hnew, hempty, ↓reduceDIte, hroot,
                              hequal, hgreater, hanchor, ha, hanchorMono] at hrun
                            simp only [hequalAnchor] at hrun
                            cases hset : pairVecDivVHCSetNext newNode none nodes with
                            | error fault => simp [hset] at hrun
                            | ok updated =>
                                rw [hset] at hrun
                                cases hbubble : pairVecDivVHCBubbleBelow
                                    heap.size anchor newNode (heap.push newNode) with
                                | error fault => simp [hbubble] at hrun
                                | ok shifted =>
                                    rw [hbubble] at hrun
                                    simp only [Except.ok.injEq,
                                      Prod.mk.injEq] at hrun
                                    rcases hrun with ⟨rfl, rfl⟩
                                    exact ⟨_,
                                      pairVecDivVHCHeapChainOwnership_bubbleBelow_fresh
                                        heap heap' owners anchor newNode nodes
                                        nodes' node mono hownership hnode hmono
                                        hfreshHead hfreshOwners hset hbubble⟩
                    · simp [pairVecDivVHCInsert, hnew, hempty, hroot, hequal,
                        hgreater, hanchor, ha] at hrun

theorem pairVecDivVHCInsert_preserves_away_of_fresh
    (newNode : Nat) (heap heap' : Array Nat)
    (nodes nodes' : Array PairVecDivVHCNode)
    (protectedSet remaining : Finset Nat)
    (owners : Nat → Finset Nat) (node : PairVecDivVHCNode)
    (mono : UMonomial)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (haway : ∀ (slot head : Nat), heap[slot]? = some head →
      Disjoint protectedSet (owners head) ∧ head ∉ protectedSet)
    (hremaining : remaining ⊆ protectedSet)
    (hnewRemaining : newNode ∉ remaining)
    (hnode : nodes[newNode]? = some node) (hmono : node.mono = some mono)
    (hfreshHead : ∀ (slot : Nat), heap[slot]? ≠ some newNode)
    (hfreshOwners : ∀ (slot head : Nat), heap[slot]? = some head →
      newNode ∉ owners head)
    (hrun : pairVecDivVHCInsert newNode heap nodes = .ok (heap', nodes')) :
    PairVecDivVHCHeapChainsOwnedAway heap' nodes' remaining := by
  cases hnew : pairVecDivVHCMono newNode nodes with
  | error fault => simp [pairVecDivVHCInsert, hnew] at hrun
  | ok newMono =>
      by_cases hempty : heap.size = 0
      · have heq : heap = #[] := Array.eq_empty_of_size_eq_zero hempty
        subst heap
        unfold pairVecDivVHCInsert at hrun
        simp only [hnew, Array.size_empty, ↓reduceDIte] at hrun
        cases hset : pairVecDivVHCSetNext newNode none nodes with
        | error fault => simp [hset] at hrun
        | ok updated =>
            rw [hset] at hrun
            simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
            rcases hrun with ⟨rfl, rfl⟩
            exact pairVecDivVHCHeapChainsOwnedAway_push_fresh #[] owners
              newNode nodes updated node mono protectedSet
              remaining hownership (by simp)
              hremaining hnewRemaining hnode hmono (by simp) (by simp) hset
      · have hheap : 0 < heap.size := Nat.pos_of_ne_zero hempty
        cases hroot : pairVecDivVHCMono heap[0] nodes with
        | error fault =>
            simp [pairVecDivVHCInsert, hnew, hempty, hroot] at hrun
        | ok rootMono =>
            by_cases hequal : newMono.deg = rootMono.deg
            · unfold pairVecDivVHCInsert at hrun
              simp only [hnew, hempty, ↓reduceDIte, hroot, hequal] at hrun
              cases hset : pairVecDivVHCSetNext newNode (some heap[0]) nodes with
              | error fault => simp [hset] at hrun
              | ok updated =>
                  rw [hset] at hrun
                  simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
                  rcases hrun with ⟨rfl, rfl⟩
                  exact pairVecDivVHCHeapChainsOwnedAway_merge_fresh heap
                    owners 0 newNode nodes updated node mono protectedSet
                    remaining hheap hownership haway
                    hremaining hnewRemaining hnode hmono hfreshHead
                    hfreshOwners hset
            · by_cases hgreater : newMono.deg > rootMono.deg
              · unfold pairVecDivVHCInsert at hrun
                simp only [hnew, hempty, ↓reduceDIte, hroot, hequal,
                  hgreater] at hrun
                cases hset : pairVecDivVHCSetNext newNode none nodes with
                | error fault => simp [hset] at hrun
                | ok updated =>
                    rw [hset] at hrun
                    cases hbubble : pairVecDivVHCBubble heap.size 0 newNode
                        (heap.push newNode) with
                    | error fault => simp [hbubble] at hrun
                    | ok shifted =>
                        rw [hbubble] at hrun
                        simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
                        rcases hrun with ⟨rfl, rfl⟩
                        exact pairVecDivVHCHeapChainsOwnedAway_bubble_fresh heap
                          shifted owners 0 newNode nodes updated node mono
                          protectedSet remaining hownership
                          haway hremaining hnewRemaining hnode hmono
                          hfreshHead hfreshOwners hset hbubble
              · cases hanchor : pairVecDivVHCFindAnchor newMono.deg
                    (pairVecDivVHCParent heap.size) heap nodes with
                | error fault =>
                    simp [pairVecDivVHCInsert, hnew, hempty, hroot, hequal,
                      hgreater, hanchor] at hrun
                | ok anchor =>
                    by_cases ha : anchor < heap.size
                    · cases hanchorMono : pairVecDivVHCMono heap[anchor] nodes with
                      | error fault =>
                          simp [pairVecDivVHCInsert, hnew, hempty, hroot,
                            hequal, hgreater, hanchor, ha, hanchorMono] at hrun
                      | ok anchorMono =>
                          by_cases hequalAnchor :
                              newMono.deg = anchorMono.deg
                          · unfold pairVecDivVHCInsert at hrun
                            simp only [hnew, hempty, ↓reduceDIte, hroot,
                              hequal, hgreater, hanchor, ha, hanchorMono] at hrun
                            simp only [hequalAnchor] at hrun
                            cases hset : pairVecDivVHCSetNext newNode
                                (some heap[anchor]) nodes with
                            | error fault => simp [hset] at hrun
                            | ok updated =>
                                rw [hset] at hrun
                                simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
                                rcases hrun with ⟨rfl, rfl⟩
                                exact pairVecDivVHCHeapChainsOwnedAway_merge_fresh
                                  heap owners anchor newNode nodes nodes' node
                                  mono protectedSet remaining
                                  ha hownership haway hremaining hnewRemaining
                                  hnode hmono hfreshHead hfreshOwners hset
                          · unfold pairVecDivVHCInsert at hrun
                            simp only [hnew, hempty, ↓reduceDIte, hroot,
                              hequal, hgreater, hanchor, ha, hanchorMono] at hrun
                            simp only [hequalAnchor] at hrun
                            cases hset : pairVecDivVHCSetNext newNode none nodes with
                            | error fault => simp [hset] at hrun
                            | ok updated =>
                                rw [hset] at hrun
                                cases hbubble : pairVecDivVHCBubbleBelow
                                    heap.size anchor newNode (heap.push newNode) with
                                | error fault => simp [hbubble] at hrun
                                | ok shifted =>
                                    rw [hbubble] at hrun
                                    simp only [Except.ok.injEq,
                                      Prod.mk.injEq] at hrun
                                    rcases hrun with ⟨rfl, rfl⟩
                                    exact
                                      pairVecDivVHCHeapChainsOwnedAway_bubbleBelow_fresh
                                        heap heap' owners anchor newNode nodes
                                        nodes' node mono protectedSet
                                        remaining hownership
                                        haway hremaining hnewRemaining hnode hmono
                                        hfreshHead hfreshOwners hset hbubble
                    · simp [pairVecDivVHCInsert, hnew, hempty, hroot, hequal,
                        hgreater, hanchor, ha] at hrun

theorem pairVecDivVHCInsert_preserves_away_of_protected
    (newNode : Nat) (heap heap' : Array Nat)
    (nodes nodes' : Array PairVecDivVHCNode)
    (protectedSet remaining : Finset Nat)
    (owners : Nat → Finset Nat) (node : PairVecDivVHCNode)
    (mono : UMonomial)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (haway : ∀ (slot head : Nat), heap[slot]? = some head →
      Disjoint protectedSet (owners head) ∧ head ∉ protectedSet)
    (hprotected : newNode ∈ protectedSet)
    (hremaining : remaining ⊆ protectedSet)
    (hnewRemaining : newNode ∉ remaining)
    (hnode : nodes[newNode]? = some node) (hmono : node.mono = some mono)
    (hrun : pairVecDivVHCInsert newNode heap nodes = .ok (heap', nodes')) :
    PairVecDivVHCHeapChainsOwnedAway heap' nodes' remaining := by
  have hfreshHead : ∀ (slot : Nat), heap[slot]? ≠ some newNode := by
    intro slot hget
    exact (haway slot newNode hget).2 hprotected
  have hfreshOwners : ∀ (slot head : Nat), heap[slot]? = some head →
      newNode ∉ owners head := by
    intro slot head hget hmem
    exact Finset.disjoint_left.mp (haway slot head hget).1 hprotected hmem
  exact pairVecDivVHCInsert_preserves_away_of_fresh newNode heap heap' nodes
    nodes' protectedSet remaining owners node mono hownership haway hremaining
    hnewRemaining hnode hmono hfreshHead hfreshOwners hrun

theorem pairVecDivVHCInsert_preserves_owned_of_away
    (newNode : Nat) (heap heap' : Array Nat)
    (nodes nodes' : Array PairVecDivVHCNode) (protectedSet : Finset Nat)
    (node : PairVecDivVHCNode) (mono : UMonomial)
    (haway : PairVecDivVHCHeapChainsOwnedAway heap nodes protectedSet)
    (hprotected : newNode ∈ protectedSet)
    (hnode : nodes[newNode]? = some node) (hmono : node.mono = some mono)
    (hrun : pairVecDivVHCInsert newNode heap nodes = .ok (heap', nodes')) :
    PairVecDivVHCHeapChainsOwned heap' nodes' := by
  rcases haway with ⟨owners, hownership, hseparated⟩
  have hfreshHead : ∀ (slot : Nat), heap[slot]? ≠ some newNode := by
    intro slot hget
    exact (hseparated slot newNode hget).2 hprotected
  have hfreshOwners : ∀ (slot head : Nat), heap[slot]? = some head →
      newNode ∉ owners head := by
    intro slot head hget hmem
    exact Finset.disjoint_left.mp (hseparated slot head hget).1 hprotected hmem
  exact pairVecDivVHCInsert_preserves_heapChainOwnership_of_fresh newNode heap
    heap' nodes nodes' owners node mono hownership hnode hmono hfreshHead
    hfreshOwners hrun

theorem pairVecDivVHCHeapChainOwnership_root_onlyAt
    (heap : Array Nat) (owners : Nat → Finset Nat)
    (nodes : Array PairVecDivVHCNode)
    (howns : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hnonempty : 0 < heap.size) :
    PairVecDivVHCOnlyAt heap heap[0] 0 := by
  intro slot hget
  have hslot : slot < heap.size := by
    by_contra hnot
    rw [Array.getElem?_eq_none (by omega)] at hget
    contradiction
  exact howns.2.1 slot 0 heap[0] hget
    (by rw [Array.getElem?_eq_getElem hnonempty])

theorem pairVecDivVHCHeapChainOwnership_root_owns
    (heap : Array Nat) (owners : Nat → Finset Nat)
    (nodes : Array PairVecDivVHCNode)
    (howns : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hnonempty : 0 < heap.size) :
    PairVecDivVHCChainOwns (some heap[0]) (owners heap[0]) nodes :=
  howns.1 0 heap[0] (by rw [Array.getElem?_eq_getElem hnonempty])

theorem pairVecDivVHCHeapChainOwnership_root_disjoint
    (heap : Array Nat) (owners : Nat → Finset Nat)
    (nodes : Array PairVecDivVHCNode)
    (howns : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hnonempty : 0 < heap.size) (slot : Nat) (hslot : slot < heap.size)
    (hne : heap[slot] ≠ heap[0]) :
    Disjoint (owners heap[slot]) (owners heap[0]) :=
  howns.2.2 slot 0 heap[slot] heap[0]
    (by rw [Array.getElem?_eq_getElem hslot])
    (by rw [Array.getElem?_eq_getElem hnonempty]) hne

theorem pairVecDivVHCHeapChainOwnership_fresh_of_mono_none
    (heap : Array Nat) (owners : Nat → Finset Nat)
    (nodes : Array PairVecDivVHCNode) (nodeIndex : Nat)
    (node : PairVecDivVHCNode)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hnode : nodes[nodeIndex]? = some node) (hmono : node.mono = none) :
    (∀ (slot : Nat), heap[slot]? ≠ some nodeIndex) ∧
      (∀ (slot head : Nat), heap[slot]? = some head →
        nodeIndex ∉ owners head) := by
  have hnotActive : ¬ ∃ activeNode mono,
      nodes[nodeIndex]? = some activeNode ∧ activeNode.mono = some mono := by
    rintro ⟨activeNode, mono, hactiveGet, hactiveMono⟩
    rw [hnode] at hactiveGet
    simp only [Option.some.injEq] at hactiveGet
    subst activeNode
    rw [hmono] at hactiveMono
    contradiction
  refine ⟨?_, ?_⟩
  · intro slot hget
    have hchain := hownership.1 slot nodeIndex hget
    have hmem := pairVecDivVHCChainOwns_head_mem nodeIndex
      (owners nodeIndex) nodes hchain
    exact hnotActive (pairVecDivVHCChainOwns_mem_active (some nodeIndex)
      (owners nodeIndex) nodes hchain nodeIndex hmem)
  · intro slot head hget hmem
    have hchain := hownership.1 slot head hget
    exact hnotActive (pairVecDivVHCChainOwns_mem_active (some head)
      (owners head) nodes hchain nodeIndex hmem)

theorem pairVecDivVHCActivate_get_ne
    (nodeIndex : Nat) (nodes nodes' : Array PairVecDivVHCNode)
    (quotient divisor : SparsePolyZp)
    (hrun : pairVecDivVHCActivate nodeIndex nodes quotient divisor = .ok nodes')
    (i : Nat) (hne : nodeIndex ≠ i) :
    nodes'[i]? = nodes[i]? := by
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
        exact Array.getElem?_set_ne hn hne

theorem pairVecDivVHCActivate_nodes_size
    (nodeIndex : Nat) (nodes nodes' : Array PairVecDivVHCNode)
    (quotient divisor : SparsePolyZp)
    (hrun : pairVecDivVHCActivate nodeIndex nodes quotient divisor = .ok nodes') :
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

theorem pairVecDivVHCActivate_preserves_mono_read_ne
    (nodeIndex : Nat) (nodes nodes' : Array PairVecDivVHCNode)
    (quotient divisor : SparsePolyZp)
    (hrun : pairVecDivVHCActivate nodeIndex nodes quotient divisor = .ok nodes')
    (i : Nat) (hne : nodeIndex ≠ i) :
    pairVecDivVHCMono i nodes' = pairVecDivVHCMono i nodes := by
  have hsize := pairVecDivVHCActivate_nodes_size nodeIndex nodes nodes'
    quotient divisor hrun
  unfold pairVecDivVHCMono
  by_cases hi : i < nodes.size
  · have hi' : i < nodes'.size := by omega
    simp only [hi, hi', ↓reduceDIte]
    have hget := pairVecDivVHCActivate_get_ne nodeIndex nodes nodes' quotient
      divisor hrun i hne
    rw [Array.getElem?_eq_getElem hi, Array.getElem?_eq_getElem hi'] at hget
    rw [Option.some.inj hget]
  · have hi' : ¬ i < nodes'.size := by omega
    simp only [hi, hi', ↓reduceDIte]

theorem pairVecDivVHCActivate_preserves_heapOrdered_of_freshHead
    (nodeIndex : Nat) (heap : Array Nat)
    (nodes nodes' : Array PairVecDivVHCNode)
    (quotient divisor : SparsePolyZp)
    (hfresh : ∀ slot : Nat, heap[slot]? ≠ some nodeIndex)
    (hordered : PairVecDivVHCHeapOrdered heap nodes)
    (hrun : pairVecDivVHCActivate nodeIndex nodes quotient divisor = .ok nodes') :
    PairVecDivVHCHeapOrdered heap nodes' := by
  apply PairVecDivVHCHeapDegreesOrderedUpTo.toHeapOrdered
  intro child hchild hpos childHead parentHead childMono parentMono
    hchildGet hparentGet hchildMono hparentMono
  have hchildNe : nodeIndex ≠ childHead := by
    intro heq
    subst childHead
    exact hfresh child hchildGet
  have hparentNe : nodeIndex ≠ parentHead := by
    intro heq
    subst parentHead
    exact hfresh (pairVecDivVHCParent child) hparentGet
  rw [pairVecDivVHCActivate_preserves_mono_read_ne nodeIndex nodes nodes'
    quotient divisor hrun childHead hchildNe] at hchildMono
  rw [pairVecDivVHCActivate_preserves_mono_read_ne nodeIndex nodes nodes'
    quotient divisor hrun parentHead hparentNe] at hparentMono
  exact hordered.degreesUpTo heap nodes heap.size (Nat.le_refl _) child hchild
    hpos childHead parentHead childMono parentMono hchildGet hparentGet
    hchildMono hparentMono

theorem pairVecDivVHCActivate_preserves_heapChainsHomogeneous_of_fresh
    (nodeIndex : Nat) (heap : Array Nat) (owners : Nat → Finset Nat)
    (nodes nodes' : Array PairVecDivVHCNode)
    (quotient divisor : SparsePolyZp)
    (hhomogeneous : PairVecDivVHCHeapChainsHomogeneous heap owners nodes)
    (hfreshHead : ∀ slot : Nat, heap[slot]? ≠ some nodeIndex)
    (hfreshOwners : ∀ (slot head : Nat), heap[slot]? = some head →
      nodeIndex ∉ owners head)
    (hrun : pairVecDivVHCActivate nodeIndex nodes quotient divisor = .ok nodes') :
    PairVecDivVHCHeapChainsHomogeneous heap owners nodes' := by
  intro slot head mono hget hmono
  have hheadNe : nodeIndex ≠ head := by
    intro heq
    subst head
    exact hfreshHead slot hget
  rw [pairVecDivVHCActivate_preserves_mono_read_ne nodeIndex nodes nodes'
    quotient divisor hrun head hheadNe] at hmono
  have hold := hhomogeneous slot head mono hget hmono
  exact pairVecDivVHCChainAtDegree_congr_on (some head) (owners head) nodes
    nodes' mono.deg hold (by
      intro i hi
      exact pairVecDivVHCActivate_get_ne nodeIndex nodes nodes' quotient divisor
        hrun i (by
          intro heq
          subst i
          exact hfreshOwners slot head hget hi))

theorem pairVecDivVHCActivate_preserves_cursorPrefixAbove
    (degreeLimit nodeIndex : Nat) (nodes nodes' : Array PairVecDivVHCNode)
    (quotient divisor : SparsePolyZp)
    (hprefix : PairVecDivVHCCursorPrefixAbove degreeLimit nodes quotient divisor)
    (hrun : pairVecDivVHCActivate nodeIndex nodes quotient divisor =
      .ok nodes') :
    PairVecDivVHCCursorPrefixAbove degreeLimit nodes' quotient divisor := by
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
        exact hprefix.set_fields degreeLimit nodeIndex nodes quotient divisor
          nodes[nodeIndex] _ hn (Array.getElem?_eq_getElem hn) rfl rfl

theorem pairVecDivVHCActivateInsert_preserves_heapChainsOwned
    (nodeIndex : Nat) (heap heap' : Array Nat)
    (nodes activated inserted : Array PairVecDivVHCNode)
    (quotient divisor : SparsePolyZp) (oldNode : PairVecDivVHCNode)
    (howned : PairVecDivVHCHeapChainsOwned heap nodes)
    (hold : nodes[nodeIndex]? = some oldNode) (hinactive : oldNode.mono = none)
    (hactivate : pairVecDivVHCActivate nodeIndex nodes quotient divisor =
      .ok activated)
    (hinsert : pairVecDivVHCInsert nodeIndex heap activated =
      .ok (heap', inserted)) :
    PairVecDivVHCHeapChainsOwned heap' inserted := by
  rcases howned with ⟨owners, hownership⟩
  have hfresh := pairVecDivVHCHeapChainOwnership_fresh_of_mono_none heap owners
    nodes nodeIndex oldNode hownership hold hinactive
  have hownershipActivated :
      PairVecDivVHCHeapChainOwnership heap owners activated := by
    refine ⟨?_, hownership.2.1, hownership.2.2⟩
    intro slot head hget
    have hchain := hownership.1 slot head hget
    exact pairVecDivVHCChainOwns_congr_on (some head) (owners head) nodes
      activated hchain (by
        intro i hi
        exact pairVecDivVHCActivate_get_ne nodeIndex nodes activated quotient
          divisor hactivate i (by
            intro heq
            subst i
            exact hfresh.2 slot head hget hi))
  rcases pairVecDivVHCActivate_get nodeIndex nodes activated quotient divisor
    hactivate with ⟨hn, hq, hd, hnewGet⟩
  have holdEq : nodes[nodeIndex] = oldNode := by
    rw [Array.getElem?_eq_getElem hn] at hold
    exact Option.some.inj hold
  have hqOld : oldNode.quotientIndex < quotient.size := by
    simpa [holdEq] using hq
  have hdOld : oldNode.divisorIndex < divisor.size := by
    simpa [holdEq] using hd
  let newNode : PairVecDivVHCNode := { oldNode with
    mono := some ⟨(quotient[oldNode.quotientIndex]'hqOld).1.deg +
      (divisor[oldNode.divisorIndex]'hdOld).1.deg⟩
    next := none }
  have hnewGet' : activated[nodeIndex]? = some newNode := by
    simpa [newNode, holdEq] using hnewGet
  have hnewMono : newNode.mono = some
      ⟨(quotient[oldNode.quotientIndex]'hqOld).1.deg +
        (divisor[oldNode.divisorIndex]'hdOld).1.deg⟩ := by
    rfl
  exact pairVecDivVHCInsert_preserves_heapChainOwnership_of_fresh nodeIndex
    heap heap' activated inserted owners newNode _ hownershipActivated hnewGet'
    hnewMono hfresh.1 hfresh.2 hinsert

theorem pairVecDivVHCActivateInsert_preserves_away
    (nodeIndex : Nat) (heap heap' : Array Nat)
    (nodes activated inserted : Array PairVecDivVHCNode)
    (quotient divisor : SparsePolyZp) (oldNode : PairVecDivVHCNode)
    (protectedSet : Finset Nat)
    (haway : PairVecDivVHCHeapChainsOwnedAway heap nodes protectedSet)
    (hprotectedActive : ∀ i ∈ protectedSet,
      ∃ node mono, nodes[i]? = some node ∧ node.mono = some mono)
    (hold : nodes[nodeIndex]? = some oldNode) (hinactive : oldNode.mono = none)
    (hactivate : pairVecDivVHCActivate nodeIndex nodes quotient divisor =
      .ok activated)
    (hinsert : pairVecDivVHCInsert nodeIndex heap activated =
      .ok (heap', inserted)) :
    PairVecDivVHCHeapChainsOwnedAway heap' inserted protectedSet := by
  rcases haway with ⟨owners, hownership, hseparated⟩
  have hfresh := pairVecDivVHCHeapChainOwnership_fresh_of_mono_none heap owners
    nodes nodeIndex oldNode hownership hold hinactive
  have hnotProtected : nodeIndex ∉ protectedSet := by
    intro hmem
    rcases hprotectedActive nodeIndex hmem with
      ⟨activeNode, mono, hactiveGet, hactiveMono⟩
    rw [hold] at hactiveGet
    simp only [Option.some.injEq] at hactiveGet
    subst activeNode
    rw [hinactive] at hactiveMono
    contradiction
  have hownershipActivated :
      PairVecDivVHCHeapChainOwnership heap owners activated := by
    refine ⟨?_, hownership.2.1, hownership.2.2⟩
    intro slot head hget
    have hchain := hownership.1 slot head hget
    exact pairVecDivVHCChainOwns_congr_on (some head) (owners head) nodes
      activated hchain (by
        intro i hi
        exact pairVecDivVHCActivate_get_ne nodeIndex nodes activated quotient
          divisor hactivate i (by
            intro heq
            subst i
            exact hfresh.2 slot head hget hi))
  rcases pairVecDivVHCActivate_get nodeIndex nodes activated quotient divisor
    hactivate with ⟨hn, hq, hd, hnewGet⟩
  have holdEq : nodes[nodeIndex] = oldNode := by
    rw [Array.getElem?_eq_getElem hn] at hold
    exact Option.some.inj hold
  have hqOld : oldNode.quotientIndex < quotient.size := by
    simpa [holdEq] using hq
  have hdOld : oldNode.divisorIndex < divisor.size := by
    simpa [holdEq] using hd
  let newNode : PairVecDivVHCNode := { oldNode with
    mono := some ⟨(quotient[oldNode.quotientIndex]'hqOld).1.deg +
      (divisor[oldNode.divisorIndex]'hdOld).1.deg⟩
    next := none }
  have hnewGet' : activated[nodeIndex]? = some newNode := by
    simpa [newNode, holdEq] using hnewGet
  have hnewMono : newNode.mono = some
      ⟨(quotient[oldNode.quotientIndex]'hqOld).1.deg +
        (divisor[oldNode.divisorIndex]'hdOld).1.deg⟩ := by
    rfl
  exact pairVecDivVHCInsert_preserves_away_of_fresh nodeIndex heap heap'
    activated inserted protectedSet protectedSet owners newNode _
    hownershipActivated hseparated Finset.Subset.rfl hnotProtected hnewGet'
    hnewMono hfresh.1 hfresh.2 hinsert

theorem pairVecDivVHCAllActiveNodesBelow.heapBelow (degreeLimit : Nat)
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode)
    (hbelow : PairVecDivVHCAllActiveNodesBelow degreeLimit nodes) :
    ∀ (slot nodeIndex : Nat) (node : PairVecDivVHCNode) (mono : UMonomial),
      heap[slot]? = some nodeIndex → nodes[nodeIndex]? = some node →
      node.mono = some mono → mono.deg < degreeLimit := by
  intro slot nodeIndex node mono hheap hnode hmono
  exact hbelow nodeIndex node mono hnode hmono

theorem pairVecDivVHCSetNext_preserves_allActiveNodesBelow
    (degreeLimit nodeIndex : Nat) (next : Option Nat)
    (nodes nodes' : Array PairVecDivVHCNode)
    (hbelow : PairVecDivVHCAllActiveNodesBelow degreeLimit nodes)
    (hrun : pairVecDivVHCSetNext nodeIndex next nodes = .ok nodes') :
    PairVecDivVHCAllActiveNodesBelow degreeLimit nodes' := by
  unfold pairVecDivVHCSetNext at hrun
  split at hrun <;> try contradiction
  next hn =>
    simp only [Except.ok.injEq] at hrun
    subst nodes'
    intro i node mono hget hmono
    by_cases heq : nodeIndex = i
    · subst i
      rw [Array.getElem?_set_self hn] at hget
      simp only [Option.some.injEq] at hget
      subst node
      exact hbelow nodeIndex nodes[nodeIndex] mono
        (Array.getElem?_eq_getElem hn) hmono
    · rw [Array.getElem?_set_ne hn heq] at hget
      exact hbelow i node mono hget hmono

theorem pairVecDivVHCSetNext_preserves_divisorIndicesFixed
    (nodeIndex : Nat) (next : Option Nat)
    (nodes nodes' : Array PairVecDivVHCNode)
    (hfixed : PairVecDivVHCNodeDivisorIndicesFixed nodes)
    (hrun : pairVecDivVHCSetNext nodeIndex next nodes = .ok nodes') :
    PairVecDivVHCNodeDivisorIndicesFixed nodes' := by
  unfold pairVecDivVHCSetNext at hrun
  split at hrun <;> try contradiction
  next hn =>
    simp only [Except.ok.injEq] at hrun
    subst nodes'
    exact PairVecDivVHCNodeDivisorIndicesFixed.set nodes nodeIndex _ hn hfixed
      rfl

theorem pairVecDivVHCSetNext_preserves_denotes
    (nodeIndex : Nat) (next : Option Nat)
    (nodes nodes' : Array PairVecDivVHCNode)
    (quotient divisor : SparsePolyZp)
    (hdenotes : ∀ (i : Nat) (node : PairVecDivVHCNode),
      nodes[i]? = some node → node.mono ≠ none →
        PairVecDivVHCNodeDenotes quotient divisor node)
    (hrun : pairVecDivVHCSetNext nodeIndex next nodes = .ok nodes') :
    ∀ (i : Nat) (node : PairVecDivVHCNode),
      nodes'[i]? = some node → node.mono ≠ none →
        PairVecDivVHCNodeDenotes quotient divisor node := by
  unfold pairVecDivVHCSetNext at hrun
  split at hrun <;> try contradiction
  next hn =>
    simp only [Except.ok.injEq] at hrun
    subst nodes'
    intro i node hget hactive
    by_cases heq : nodeIndex = i
    · subst i
      rw [Array.getElem?_set_self hn] at hget
      simp only [Option.some.injEq] at hget
      subst node
      have hold := hdenotes nodeIndex nodes[nodeIndex]
        (Array.getElem?_eq_getElem hn) (by simpa using hactive)
      simpa [PairVecDivVHCNodeDenotes] using hold
    · rw [Array.getElem?_set_ne hn heq] at hget
      exact hdenotes i node hget hactive

theorem pairVecDivVHCInsert_nodes_result
    (newNode : Nat) (heap heap' : Array Nat)
    (nodes nodes' : Array PairVecDivVHCNode)
    (hrun : pairVecDivVHCInsert newNode heap nodes = .ok (heap', nodes')) :
    ∃ next, pairVecDivVHCSetNext newNode next nodes = .ok nodes' := by
  cases hnew : pairVecDivVHCMono newNode nodes with
  | error fault => simp [pairVecDivVHCInsert, hnew] at hrun
  | ok newMono =>
      by_cases hempty : heap.size = 0
      · unfold pairVecDivVHCInsert at hrun
        simp only [hnew, hempty, ↓reduceDIte] at hrun
        cases hset : pairVecDivVHCSetNext newNode none nodes with
        | error fault => simp [hset] at hrun
        | ok updated =>
            rw [hset] at hrun
            simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
            rcases hrun with ⟨rfl, rfl⟩
            exact ⟨none, hset⟩
      · have hheap : 0 < heap.size := Nat.pos_of_ne_zero hempty
        cases hroot : pairVecDivVHCMono heap[0] nodes with
        | error fault =>
            simp [pairVecDivVHCInsert, hnew, hempty, hroot] at hrun
        | ok rootMono =>
            by_cases hequal : newMono.deg = rootMono.deg
            · unfold pairVecDivVHCInsert at hrun
              simp only [hnew, hempty, ↓reduceDIte, hroot, hequal] at hrun
              cases hset : pairVecDivVHCSetNext newNode (some heap[0]) nodes with
              | error fault => simp [hset] at hrun
              | ok updated =>
                  rw [hset] at hrun
                  simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
                  rcases hrun with ⟨rfl, rfl⟩
                  exact ⟨some heap[0], hset⟩
            · by_cases hgreater : newMono.deg > rootMono.deg
              · unfold pairVecDivVHCInsert at hrun
                simp only [hnew, hempty, ↓reduceDIte, hroot, hequal,
                  hgreater] at hrun
                cases hset : pairVecDivVHCSetNext newNode none nodes with
                | error fault => simp [hset] at hrun
                | ok updated =>
                    rw [hset] at hrun
                    cases hbubble : pairVecDivVHCBubble heap.size 0 newNode
                        (heap.push newNode) with
                    | error fault => simp [hbubble] at hrun
                    | ok shifted =>
                        rw [hbubble] at hrun
                        simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
                        rcases hrun with ⟨rfl, rfl⟩
                        exact ⟨none, hset⟩
              · cases hanchor : pairVecDivVHCFindAnchor newMono.deg
                    (pairVecDivVHCParent heap.size) heap nodes with
                | error fault =>
                    simp [pairVecDivVHCInsert, hnew, hempty, hroot, hequal,
                      hgreater, hanchor] at hrun
                | ok anchor =>
                    by_cases ha : anchor < heap.size
                    · cases hanchorMono : pairVecDivVHCMono heap[anchor] nodes with
                      | error fault =>
                          simp [pairVecDivVHCInsert, hnew, hempty, hroot,
                            hequal, hgreater, hanchor, ha, hanchorMono] at hrun
                      | ok anchorMono =>
                          by_cases hequalAnchor :
                              newMono.deg = anchorMono.deg
                          · unfold pairVecDivVHCInsert at hrun
                            simp only [hnew, hempty, ↓reduceDIte, hroot,
                              hequal, hgreater, hanchor, ha, hanchorMono] at hrun
                            simp only [hequalAnchor] at hrun
                            cases hset : pairVecDivVHCSetNext newNode
                                (some heap[anchor]) nodes with
                            | error fault => simp [hset] at hrun
                            | ok updated =>
                                rw [hset] at hrun
                                simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
                                rcases hrun with ⟨rfl, rfl⟩
                                exact ⟨some heap[anchor], hset⟩
                          · unfold pairVecDivVHCInsert at hrun
                            simp only [hnew, hempty, ↓reduceDIte, hroot,
                              hequal, hgreater, hanchor, ha, hanchorMono] at hrun
                            simp only [hequalAnchor] at hrun
                            cases hset : pairVecDivVHCSetNext newNode none nodes with
                            | error fault => simp [hset] at hrun
                            | ok updated =>
                                rw [hset] at hrun
                                cases hbubble : pairVecDivVHCBubbleBelow
                                    heap.size anchor newNode
                                    (heap.push newNode) with
                                | error fault => simp [hbubble] at hrun
                                | ok shifted =>
                                    rw [hbubble] at hrun
                                    simp only [Except.ok.injEq,
                                      Prod.mk.injEq] at hrun
                                    rcases hrun with ⟨rfl, rfl⟩
                                    exact ⟨none, hset⟩
                    · simp [pairVecDivVHCInsert, hnew, hempty, hroot, hequal,
                        hgreater, hanchor, ha] at hrun

theorem pairVecDivVHCInsert_root_of_greater
    (newNode : Nat) (heap heap' : Array Nat)
    (nodes nodes' : Array PairVecDivVHCNode)
    (newMono rootMono : UMonomial)
    (hheap : 0 < heap.size)
    (hnew : pairVecDivVHCMono newNode nodes = .ok newMono)
    (hroot : pairVecDivVHCMono heap[0] nodes = .ok rootMono)
    (hgreater : rootMono.deg < newMono.deg)
    (hrun : pairVecDivVHCInsert newNode heap nodes = .ok (heap', nodes')) :
    heap'[0]? = some newNode := by
  have hempty : heap.size ≠ 0 := Nat.ne_of_gt hheap
  have hequal : newMono.deg ≠ rootMono.deg := by omega
  unfold pairVecDivVHCInsert at hrun
  simp only [hnew, hempty, ↓reduceDIte, hroot, hequal, hgreater] at hrun
  cases hset : pairVecDivVHCSetNext newNode none nodes with
  | error fault => simp [hset] at hrun
  | ok updated =>
      rw [hset] at hrun
      cases hbubble : pairVecDivVHCBubble heap.size 0 newNode
          (heap.push newNode) with
      | error fault => simp [hbubble] at hrun
      | ok shifted =>
          rw [hbubble] at hrun
          simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
          rcases hrun with ⟨hheapResult, hnodesResult⟩
          subst heap'
          subst nodes'
          exact pairVecDivVHCBubble_stop_get heap.size 0 newNode
            (heap.push newNode) shifted hbubble

theorem pairVecDivVHCInsert_new_root_bounds_all
    (newNode : Nat) (heap heap' : Array Nat)
    (nodes nodes' : Array PairVecDivVHCNode)
    (newMono rootMono : UMonomial)
    (hheap : 0 < heap.size)
    (hvalid : PairVecDivVHCHeapPointersValid heap nodes)
    (hordered : PairVecDivVHCHeapOrdered heap nodes)
    (hnew : pairVecDivVHCMono newNode nodes = .ok newMono)
    (hroot : pairVecDivVHCMono heap[0] nodes = .ok rootMono)
    (hgreater : rootMono.deg < newMono.deg)
    (hrun : pairVecDivVHCInsert newNode heap nodes = .ok (heap', nodes')) :
    ∀ (slot head : Nat) (mono : UMonomial),
      heap'[slot]? = some head → pairVecDivVHCMono head nodes = .ok mono →
        mono.deg ≤ newMono.deg := by
  have hempty : heap.size ≠ 0 := Nat.ne_of_gt hheap
  have hequal : newMono.deg ≠ rootMono.deg := by omega
  unfold pairVecDivVHCInsert at hrun
  simp only [hnew, hempty, ↓reduceDIte, hroot, hequal, hgreater] at hrun
  cases hset : pairVecDivVHCSetNext newNode none nodes with
  | error fault => simp [hset] at hrun
  | ok updated =>
      rw [hset] at hrun
      cases hbubble : pairVecDivVHCBubble heap.size 0 newNode
          (heap.push newNode) with
      | error fault => simp [hbubble] at hrun
      | ok shifted =>
          rw [hbubble] at hrun
          simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
          rcases hrun with ⟨rfl, rfl⟩
          exact pairVecDivVHCBubble_new_root_bounds_all newNode heap shifted
            nodes newMono rootMono hheap hvalid hordered hnew hroot
            (Nat.le_of_lt hgreater) hbubble

theorem pairVecDivVHCSetNext_get_ne
    (nodeIndex : Nat) (next : Option Nat)
    (nodes nodes' : Array PairVecDivVHCNode)
    (hrun : pairVecDivVHCSetNext nodeIndex next nodes = .ok nodes')
    (i : Nat) (hne : nodeIndex ≠ i) :
    nodes'[i]? = nodes[i]? := by
  unfold pairVecDivVHCSetNext at hrun
  split at hrun <;> try contradiction
  next hn =>
    simp only [Except.ok.injEq] at hrun
    subst nodes'
    rw [Array.getElem?_set_ne hn hne]

theorem pairVecDivVHCSetNext_preserves_mono_read
    (nodeIndex : Nat) (next : Option Nat)
    (nodes nodes' : Array PairVecDivVHCNode)
    (hrun : pairVecDivVHCSetNext nodeIndex next nodes = .ok nodes')
    (i : Nat) :
    pairVecDivVHCMono i nodes' = pairVecDivVHCMono i nodes := by
  unfold pairVecDivVHCSetNext at hrun
  split at hrun <;> try contradiction
  next hn =>
    simp only [Except.ok.injEq] at hrun
    subst nodes'
    unfold pairVecDivVHCMono
    by_cases hi : i < nodes.size
    · simp only [Array.size_set, hi, ↓reduceDIte]
      by_cases heq : nodeIndex = i
      · subst i
        rw [Array.getElem_set_self]
      · rw [Array.getElem_set_ne hn hi heq]
    · have hiSet : ¬ i < (nodes.set nodeIndex
          { nodes[nodeIndex] with next := next }).size := by
          simpa only [Array.size_set] using hi
      simp only [hi, hiSet, ↓reduceDIte]

theorem pairVecDivVHCHeapChainsHomogeneous_push_fresh
    (heap : Array Nat) (owners : Nat → Finset Nat) (newNode : Nat)
    (nodes nodes' : Array PairVecDivVHCNode) (node : PairVecDivVHCNode)
    (mono : UMonomial)
    (hhomogeneous : PairVecDivVHCHeapChainsHomogeneous heap owners nodes)
    (hnode : nodes[newNode]? = some node) (hmono : node.mono = some mono)
    (hfreshHead : ∀ slot : Nat, heap[slot]? ≠ some newNode)
    (hfreshOwners : ∀ (slot head : Nat), heap[slot]? = some head →
      newNode ∉ owners head)
    (hset : pairVecDivVHCSetNext newNode none nodes = .ok nodes') :
    PairVecDivVHCHeapChainsHomogeneous (heap.push newNode)
      (fun head => if head = newNode then {newNode} else owners head) nodes' := by
  intro slot head headMono hget hheadMono
  rw [Array.getElem?_push] at hget
  by_cases hlast : slot = heap.size
  · subst slot
    simp only [ite_true, Option.some.injEq] at hget
    subst head
    have hnewRun : pairVecDivVHCMono newNode nodes = .ok mono :=
      (pairVecDivVHCMono_eq_ok_iff newNode nodes mono).2 ⟨node, hnode, hmono⟩
    rw [pairVecDivVHCSetNext_preserves_mono_read newNode none nodes nodes' hset
      newNode, hnewRun] at hheadMono
    have hdegree : mono.deg = headMono.deg :=
      congrArg UMonomial.deg (Except.ok.inj hheadMono)
    simpa only [hdegree] using
      pairVecDivVHCSetNext_chainAtDegree_insert newNode none nodes
      nodes' node mono ∅ mono.deg hnode hmono rfl (by simp)
      (by simp [PairVecDivVHCChainAtDegree]) hset
  · simp only [hlast, ↓reduceIte] at hget
    have hheadNe : head ≠ newNode := by
      intro heq
      subst head
      exact hfreshHead slot hget
    simp only [hheadNe, ↓reduceIte]
    rw [pairVecDivVHCSetNext_preserves_mono_read newNode none nodes nodes' hset
      head] at hheadMono
    have hold := hhomogeneous slot head headMono hget hheadMono
    exact pairVecDivVHCChainAtDegree_congr_on (some head) (owners head) nodes
      nodes' headMono.deg hold (by
        intro i hi
        exact pairVecDivVHCSetNext_get_ne newNode none nodes nodes' hset i
          (by
            intro heq
            subst i
            exact hfreshOwners slot head hget hi))

theorem pairVecDivVHCHeapChainsHomogeneous_merge_fresh
    (heap : Array Nat) (owners : Nat → Finset Nat) (slot newNode : Nat)
    (nodes nodes' : Array PairVecDivVHCNode) (node : PairVecDivVHCNode)
    (newMono oldMono : UMonomial) (hslot : slot < heap.size)
    (hhomogeneous : PairVecDivVHCHeapChainsHomogeneous heap owners nodes)
    (hnode : nodes[newNode]? = some node) (hmono : node.mono = some newMono)
    (holdMono : pairVecDivVHCMono heap[slot] nodes = .ok oldMono)
    (hdegree : newMono.deg = oldMono.deg)
    (hfreshHead : ∀ i : Nat, heap[i]? ≠ some newNode)
    (hfreshOwners : ∀ (i head : Nat), heap[i]? = some head →
      newNode ∉ owners head)
    (hset : pairVecDivVHCSetNext newNode (some heap[slot]) nodes = .ok nodes') :
    PairVecDivVHCHeapChainsHomogeneous (heap.set slot newNode)
      (fun head => if head = newNode then insert newNode (owners heap[slot])
        else owners head) nodes' := by
  intro targetSlot head headMono hget hheadMono
  by_cases hat : slot = targetSlot
  · subst targetSlot
    rw [Array.getElem?_set_self hslot] at hget
    simp only [Option.some.injEq] at hget
    subst head
    have hnewRun : pairVecDivVHCMono newNode nodes = .ok newMono :=
      (pairVecDivVHCMono_eq_ok_iff newNode nodes newMono).2
        ⟨node, hnode, hmono⟩
    rw [pairVecDivVHCSetNext_preserves_mono_read newNode
      (some heap[slot]) nodes nodes' hset newNode, hnewRun] at hheadMono
    have hheadDegree : newMono.deg = headMono.deg :=
      congrArg UMonomial.deg (Except.ok.inj hheadMono)
    have htail := hhomogeneous slot heap[slot] oldMono
      (Array.getElem?_eq_getElem hslot) holdMono
    have hnewDegree := pairVecDivVHCSetNext_chainAtDegree_insert newNode
      (some heap[slot]) nodes nodes' node newMono (owners heap[slot])
      oldMono.deg hnode hmono hdegree (hfreshOwners slot heap[slot]
        (Array.getElem?_eq_getElem hslot)) htail hset
    have holdHeadDegree : oldMono.deg = headMono.deg :=
      hdegree.symm.trans hheadDegree
    simpa only [holdHeadDegree] using hnewDegree
  · rw [Array.getElem?_set_ne hslot hat] at hget
    have hheadNe : head ≠ newNode := by
      intro heq
      subst head
      exact hfreshHead targetSlot hget
    simp only [hheadNe, ↓reduceIte]
    rw [pairVecDivVHCSetNext_preserves_mono_read newNode
      (some heap[slot]) nodes nodes' hset head] at hheadMono
    have hold := hhomogeneous targetSlot head headMono hget hheadMono
    exact pairVecDivVHCChainAtDegree_congr_on (some head) (owners head) nodes
      nodes' headMono.deg hold (by
        intro i hi
        exact pairVecDivVHCSetNext_get_ne newNode (some heap[slot]) nodes nodes'
          hset i (by
            intro heq
            subst i
            exact hfreshOwners targetSlot head hget hi))

theorem pairVecDivVHCInsert_preserves_heapChainsHomogeneous_of_fresh
    (newNode : Nat) (heap heap' : Array Nat)
    (nodes nodes' : Array PairVecDivVHCNode)
    (owners : Nat → Finset Nat) (node : PairVecDivVHCNode) (mono : UMonomial)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hhomogeneous : PairVecDivVHCHeapChainsHomogeneous heap owners nodes)
    (hnode : nodes[newNode]? = some node) (hmono : node.mono = some mono)
    (hfreshHead : ∀ slot : Nat, heap[slot]? ≠ some newNode)
    (hfreshOwners : ∀ (slot head : Nat), heap[slot]? = some head →
      newNode ∉ owners head)
    (hrun : pairVecDivVHCInsert newNode heap nodes = .ok (heap', nodes')) :
    ∃ owners', PairVecDivVHCHeapChainOwnership heap' owners' nodes' ∧
      PairVecDivVHCHeapChainsHomogeneous heap' owners' nodes' := by
  cases hnew : pairVecDivVHCMono newNode nodes with
  | error fault => simp [pairVecDivVHCInsert, hnew] at hrun
  | ok newMono =>
      have hmonoRun : pairVecDivVHCMono newNode nodes = .ok mono :=
        (pairVecDivVHCMono_eq_ok_iff newNode nodes mono).2 ⟨node, hnode, hmono⟩
      rw [hmonoRun] at hnew
      have hnewEq : newMono = mono := (Except.ok.inj hnew).symm
      subst newMono
      by_cases hempty : heap.size = 0
      · have heq : heap = #[] := Array.eq_empty_of_size_eq_zero hempty
        subst heap
        unfold pairVecDivVHCInsert at hrun
        simp only [hmonoRun, Array.size_empty, ↓reduceDIte] at hrun
        cases hset : pairVecDivVHCSetNext newNode none nodes with
        | error fault => simp [hset] at hrun
        | ok updated =>
            rw [hset] at hrun
            simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
            rcases hrun with ⟨rfl, rfl⟩
            let owners' := fun head => if head = newNode then {newNode}
              else owners head
            exact ⟨owners',
              pairVecDivVHCHeapChainOwnership_push_fresh #[] owners newNode
                nodes updated node mono hownership hnode hmono (by simp)
                (by simp) hset,
              pairVecDivVHCHeapChainsHomogeneous_push_fresh #[] owners newNode
                nodes updated node mono hhomogeneous hnode hmono (by simp)
                (by simp) hset⟩
      · have hheap : 0 < heap.size := Nat.pos_of_ne_zero hempty
        cases hroot : pairVecDivVHCMono heap[0] nodes with
        | error fault =>
            simp [pairVecDivVHCInsert, hmonoRun, hempty, hroot] at hrun
        | ok rootMono =>
            by_cases hequal : mono.deg = rootMono.deg
            · unfold pairVecDivVHCInsert at hrun
              simp only [hmonoRun, hempty, ↓reduceDIte, hroot, hequal] at hrun
              cases hset : pairVecDivVHCSetNext newNode (some heap[0]) nodes with
              | error fault => simp [hset] at hrun
              | ok updated =>
                  rw [hset] at hrun
                  simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
                  rcases hrun with ⟨rfl, rfl⟩
                  let owners' := fun head => if head = newNode then
                    insert newNode (owners heap[0]) else owners head
                  exact ⟨owners',
                    pairVecDivVHCHeapChainOwnership_merge_fresh heap owners 0
                      newNode nodes updated node mono hheap hownership hnode
                      hmono hfreshHead hfreshOwners hset,
                    pairVecDivVHCHeapChainsHomogeneous_merge_fresh heap owners 0
                      newNode nodes updated node mono rootMono hheap
                      hhomogeneous hnode hmono hroot hequal hfreshHead
                      hfreshOwners hset⟩
            · by_cases hgreater : mono.deg > rootMono.deg
              · unfold pairVecDivVHCInsert at hrun
                simp only [hmonoRun, hempty, ↓reduceDIte, hroot, hequal,
                  hgreater] at hrun
                cases hset : pairVecDivVHCSetNext newNode none nodes with
                | error fault => simp [hset] at hrun
                | ok updated =>
                    rw [hset] at hrun
                    cases hbubble : pairVecDivVHCBubble heap.size 0 newNode
                        (heap.push newNode) with
                    | error fault => simp [hbubble] at hrun
                    | ok shifted =>
                        rw [hbubble] at hrun
                        let owners' := fun head => if head = newNode then
                          {newNode} else owners head
                        have hpushHomogeneous :=
                          pairVecDivVHCHeapChainsHomogeneous_push_fresh heap
                            owners newNode nodes updated node mono hhomogeneous
                            hnode hmono hfreshHead hfreshOwners hset
                        have hfrom := pairVecDivVHCBubble_valuesFrom heap.size 0
                          newNode (heap.push newNode) shifted (heap.push newNode)
                          (pairVecDivVHCValuesFrom_refl _) ⟨heap.size, by simp⟩
                          hbubble
                        have hshiftHomogeneous := hpushHomogeneous.of_valuesFrom
                          (heap.push newNode) shifted owners' updated hfrom
                        have hshiftOwnership :=
                          pairVecDivVHCHeapChainOwnership_bubble_fresh heap
                            shifted owners 0 newNode nodes updated node mono
                            hownership hnode hmono hfreshHead hfreshOwners hset
                            hbubble
                        simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
                        rcases hrun with ⟨rfl, rfl⟩
                        exact ⟨owners', hshiftOwnership, hshiftHomogeneous⟩
              · cases hanchor : pairVecDivVHCFindAnchor mono.deg
                    (pairVecDivVHCParent heap.size) heap nodes with
                | error fault =>
                    simp [pairVecDivVHCInsert, hmonoRun, hempty, hroot, hequal,
                      hgreater, hanchor] at hrun
                | ok anchor =>
                    by_cases ha : anchor < heap.size
                    · cases hanchorMono : pairVecDivVHCMono heap[anchor] nodes with
                      | error fault =>
                          simp [pairVecDivVHCInsert, hmonoRun, hempty, hroot,
                            hequal, hgreater, hanchor, ha, hanchorMono] at hrun
                      | ok anchorMono =>
                          by_cases hequalAnchor : mono.deg = anchorMono.deg
                          · unfold pairVecDivVHCInsert at hrun
                            simp only [hmonoRun, hempty, ↓reduceDIte, hroot,
                              hequal, hgreater, hanchor, ha, hanchorMono] at hrun
                            simp only [hequalAnchor] at hrun
                            cases hset : pairVecDivVHCSetNext newNode
                                (some heap[anchor]) nodes with
                            | error fault => simp [hset] at hrun
                            | ok updated =>
                                rw [hset] at hrun
                                let owners' := fun head => if head = newNode then
                                  insert newNode (owners heap[anchor])
                                  else owners head
                                have hresult :
                                    PairVecDivVHCHeapChainOwnership
                                        (heap.set anchor newNode) owners' updated ∧
                                      PairVecDivVHCHeapChainsHomogeneous
                                        (heap.set anchor newNode) owners' updated := ⟨
                                  pairVecDivVHCHeapChainOwnership_merge_fresh
                                    heap owners anchor newNode nodes updated node
                                    mono ha hownership hnode hmono hfreshHead
                                    hfreshOwners hset,
                                  pairVecDivVHCHeapChainsHomogeneous_merge_fresh
                                    heap owners anchor newNode nodes updated node
                                    mono anchorMono ha hhomogeneous hnode hmono
                                    hanchorMono hequalAnchor hfreshHead
                                    hfreshOwners hset⟩
                                simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
                                rcases hrun with ⟨rfl, rfl⟩
                                exact ⟨owners', hresult⟩
                          · unfold pairVecDivVHCInsert at hrun
                            simp only [hmonoRun, hempty, ↓reduceDIte, hroot,
                              hequal, hgreater, hanchor, ha, hanchorMono] at hrun
                            simp only [hequalAnchor] at hrun
                            cases hset : pairVecDivVHCSetNext newNode none nodes with
                            | error fault => simp [hset] at hrun
                            | ok updated =>
                                rw [hset] at hrun
                                cases hbubble : pairVecDivVHCBubbleBelow
                                    heap.size anchor newNode
                                    (heap.push newNode) with
                                | error fault => simp [hbubble] at hrun
                                | ok shifted =>
                                    rw [hbubble] at hrun
                                    let owners' := fun head =>
                                      if head = newNode then {newNode}
                                      else owners head
                                    have hpushHomogeneous :=
                                      pairVecDivVHCHeapChainsHomogeneous_push_fresh
                                        heap owners newNode nodes updated node mono
                                        hhomogeneous hnode hmono hfreshHead
                                        hfreshOwners hset
                                    have hfrom :=
                                      pairVecDivVHCBubbleBelow_valuesFrom
                                        heap.size anchor newNode
                                        (heap.push newNode) shifted
                                        (heap.push newNode)
                                        (pairVecDivVHCValuesFrom_refl _)
                                        ⟨heap.size, by simp⟩ hbubble
                                    have hshiftHomogeneous :=
                                      hpushHomogeneous.of_valuesFrom
                                        (heap.push newNode) shifted owners'
                                        updated hfrom
                                    have hshiftOwnership :=
                                      pairVecDivVHCHeapChainOwnership_bubbleBelow_fresh
                                        heap shifted owners anchor newNode nodes
                                        updated node mono hownership hnode hmono
                                        hfreshHead hfreshOwners hset hbubble
                                    simp only [Except.ok.injEq,
                                      Prod.mk.injEq] at hrun
                                    rcases hrun with ⟨rfl, rfl⟩
                                    exact ⟨owners', hshiftOwnership,
                                      hshiftHomogeneous⟩
                    · simp [pairVecDivVHCInsert, hmonoRun, hempty, hroot,
                        hequal, hgreater, hanchor, ha] at hrun

theorem pairVecDivVHCSetNext_preserves_heapOrdered
    (nodeIndex : Nat) (next : Option Nat) (heap : Array Nat)
    (nodes nodes' : Array PairVecDivVHCNode)
    (hordered : PairVecDivVHCHeapOrdered heap nodes)
    (hrun : pairVecDivVHCSetNext nodeIndex next nodes = .ok nodes') :
    PairVecDivVHCHeapOrdered heap nodes' := by
  apply PairVecDivVHCHeapDegreesOrderedUpTo.toHeapOrdered
  intro child hchild hpos childHead parentHead childMono parentMono
      hchildGet hparentGet hchildMono hparentMono
  rw [pairVecDivVHCSetNext_preserves_mono_read nodeIndex next nodes nodes'
    hrun childHead] at hchildMono
  rw [pairVecDivVHCSetNext_preserves_mono_read nodeIndex next nodes nodes'
    hrun parentHead] at hparentMono
  exact hordered.degreesUpTo heap nodes heap.size (Nat.le_refl _) child hchild
    hpos childHead parentHead childMono parentMono hchildGet hparentGet
    hchildMono hparentMono

/-- The node-array side effect of insertion changes only `next`; consequently
any heap ordering established for the returned heap against the pre-write
monomial array transports to the actual returned node array. -/
theorem pairVecDivVHCInsert_nodes_preserve_heapOrdered
    (newNode : Nat) (heap heap' : Array Nat)
    (nodes nodes' : Array PairVecDivVHCNode)
    (hordered : PairVecDivVHCHeapOrdered heap' nodes)
    (hrun : pairVecDivVHCInsert newNode heap nodes = .ok (heap', nodes')) :
    PairVecDivVHCHeapOrdered heap' nodes' := by
  rcases pairVecDivVHCInsert_nodes_result newNode heap heap' nodes nodes' hrun
    with ⟨next, hset⟩
  exact pairVecDivVHCSetNext_preserves_heapOrdered newNode next heap' nodes
    nodes' hordered hset

/-- The generated greater-root insertion branch preserves max-heap order all
the way through its `next := none` node write and root bubble. -/
theorem pairVecDivVHCInsert_newRoot_preserves_heapOrdered
    (newNode : Nat) (heap heap' : Array Nat)
    (nodes nodes' : Array PairVecDivVHCNode)
    (newMono rootMono : UMonomial)
    (hheap : 0 < heap.size)
    (hvalid : PairVecDivVHCHeapPointersValid heap nodes)
    (hordered : PairVecDivVHCHeapOrdered heap nodes)
    (hnew : pairVecDivVHCMono newNode nodes = .ok newMono)
    (hroot : pairVecDivVHCMono heap[0] nodes = .ok rootMono)
    (hgreater : rootMono.deg < newMono.deg)
    (hrun : pairVecDivVHCInsert newNode heap nodes = .ok (heap', nodes')) :
    PairVecDivVHCHeapOrdered heap' nodes' := by
  have hempty : heap.size ≠ 0 := Nat.ne_of_gt hheap
  have hequal : newMono.deg ≠ rootMono.deg := by omega
  unfold pairVecDivVHCInsert at hrun
  simp only [hnew, hempty, ↓reduceDIte, hroot, hequal, hgreater] at hrun
  cases hset : pairVecDivVHCSetNext newNode none nodes with
  | error fault => simp [hset] at hrun
  | ok updated =>
      rw [hset] at hrun
      cases hbubble : pairVecDivVHCBubble heap.size 0 newNode
          (heap.push newNode) with
      | error fault => simp [hbubble] at hrun
      | ok shifted =>
          rw [hbubble] at hrun
          simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
          rcases hrun with ⟨rfl, rfl⟩
          have hshifted : PairVecDivVHCHeapOrdered shifted nodes :=
            pairVecDivVHCBubble_new_root_preserves_heapOrdered newNode heap
              shifted nodes newMono rootMono hheap hvalid hordered hnew hroot
              (Nat.le_of_lt hgreater) hbubble
          exact pairVecDivVHCSetNext_preserves_heapOrdered newNode none shifted
            nodes updated hshifted hset

/-- The generated unequal-anchor insertion branch preserves max-heap order,
using the exact comparison trace of its preceding `FindAnchor` execution. -/
theorem pairVecDivVHCInsert_bubbleBelow_preserves_heapOrdered
    (newNode anchor : Nat) (heap heap' : Array Nat)
    (nodes nodes' : Array PairVecDivVHCNode)
    (newMono rootMono anchorMono : UMonomial)
    (hheap : 0 < heap.size)
    (hvalid : PairVecDivVHCHeapPointersValid heap nodes)
    (hordered : PairVecDivVHCHeapOrdered heap nodes)
    (hnew : pairVecDivVHCMono newNode nodes = .ok newMono)
    (hroot : pairVecDivVHCMono heap[0] nodes = .ok rootMono)
    (hnequal : newMono.deg ≠ rootMono.deg)
    (hgreater : ¬ newMono.deg > rootMono.deg)
    (hanchor : pairVecDivVHCFindAnchor newMono.deg
      (pairVecDivVHCParent heap.size) heap nodes = .ok anchor)
    (ha : anchor < heap.size)
    (hanchorMono : pairVecDivVHCMono heap[anchor] nodes = .ok anchorMono)
    (hnequalAnchor : newMono.deg ≠ anchorMono.deg)
    (hrun : pairVecDivVHCInsert newNode heap nodes = .ok (heap', nodes')) :
    PairVecDivVHCHeapOrdered heap' nodes' := by
  unfold pairVecDivVHCInsert at hrun
  simp only [hnew, Nat.ne_of_gt hheap, ↓reduceDIte, hroot, hnequal,
    hgreater, hanchor, ha, hanchorMono] at hrun
  simp only [hnequalAnchor] at hrun
  cases hset : pairVecDivVHCSetNext newNode none nodes with
  | error fault => simp [hset] at hrun
  | ok updated =>
      rw [hset] at hrun
      cases hbubble : pairVecDivVHCBubbleBelow heap.size anchor newNode
          (heap.push newNode) with
      | error fault => simp [hbubble] at hrun
      | ok shifted =>
          rw [hbubble] at hrun
          have htrace := pairVecDivVHCFindAnchor_trace newMono.deg
            (pairVecDivVHCParent heap.size) anchor heap nodes hanchor
          have hshifted : PairVecDivVHCHeapOrdered shifted nodes :=
            pairVecDivVHCBubbleBelow_push_trace_preserves_heapOrdered anchor
              newNode heap shifted nodes newMono hheap hvalid hordered hnew
              htrace hbubble
          have hactual := pairVecDivVHCSetNext_preserves_heapOrdered newNode
            none shifted nodes updated hshifted hset
          simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
          rcases hrun with ⟨rfl, rfl⟩
          exact hactual

theorem pairVecDivVHCInsert_equalRoot_preserves_heapOrdered
    (newNode : Nat) (heap heap' : Array Nat)
    (nodes nodes' : Array PairVecDivVHCNode)
    (newMono rootMono : UMonomial)
    (hheap : 0 < heap.size)
    (hordered : PairVecDivVHCHeapOrdered heap nodes)
    (hnew : pairVecDivVHCMono newNode nodes = .ok newMono)
    (hroot : pairVecDivVHCMono heap[0] nodes = .ok rootMono)
    (hequal : newMono.deg = rootMono.deg)
    (hrun : pairVecDivVHCInsert newNode heap nodes = .ok (heap', nodes')) :
    PairVecDivVHCHeapOrdered heap' nodes' := by
  unfold pairVecDivVHCInsert at hrun
  simp only [hnew, Nat.ne_of_gt hheap, ↓reduceDIte, hroot, hequal] at hrun
  cases hset : pairVecDivVHCSetNext newNode (some heap[0]) nodes with
  | error fault => simp [hset] at hrun
  | ok updated =>
      rw [hset] at hrun
      have hheapOrdered : PairVecDivVHCHeapOrdered
          (heap.set 0 newNode) nodes :=
        pairVecDivVHCSet_sameDegree_preserves_heapOrdered heap nodes 0 newNode
          heap[0] newMono rootMono hheap
          (Array.getElem?_eq_getElem hheap) hnew hroot hequal hordered
      have hactual := pairVecDivVHCSetNext_preserves_heapOrdered newNode
        (some heap[0]) (heap.set 0 newNode) nodes updated hheapOrdered hset
      simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
      rcases hrun with ⟨rfl, rfl⟩
      exact hactual

theorem pairVecDivVHCInsert_equalAnchor_preserves_heapOrdered
    (newNode anchor : Nat) (heap heap' : Array Nat)
    (nodes nodes' : Array PairVecDivVHCNode)
    (newMono rootMono anchorMono : UMonomial)
    (hheap : 0 < heap.size)
    (hordered : PairVecDivVHCHeapOrdered heap nodes)
    (hnew : pairVecDivVHCMono newNode nodes = .ok newMono)
    (hroot : pairVecDivVHCMono heap[0] nodes = .ok rootMono)
    (hnequal : newMono.deg ≠ rootMono.deg)
    (hgreater : ¬ newMono.deg > rootMono.deg)
    (hanchor : pairVecDivVHCFindAnchor newMono.deg
      (pairVecDivVHCParent heap.size) heap nodes = .ok anchor)
    (ha : anchor < heap.size)
    (hanchorMono : pairVecDivVHCMono heap[anchor] nodes = .ok anchorMono)
    (hequalAnchor : newMono.deg = anchorMono.deg)
    (hrun : pairVecDivVHCInsert newNode heap nodes = .ok (heap', nodes')) :
    PairVecDivVHCHeapOrdered heap' nodes' := by
  unfold pairVecDivVHCInsert at hrun
  simp only [hnew, Nat.ne_of_gt hheap, ↓reduceDIte, hroot, hnequal,
    hgreater, hanchor, ha, hanchorMono] at hrun
  simp only [hequalAnchor] at hrun
  cases hset : pairVecDivVHCSetNext newNode (some heap[anchor]) nodes with
  | error fault => simp [hset] at hrun
  | ok updated =>
      rw [hset] at hrun
      have hheapOrdered : PairVecDivVHCHeapOrdered
          (heap.set anchor newNode) nodes :=
        pairVecDivVHCSet_sameDegree_preserves_heapOrdered heap nodes anchor
          newNode heap[anchor] newMono anchorMono ha
          (Array.getElem?_eq_getElem ha) hnew hanchorMono hequalAnchor hordered
      have hactual := pairVecDivVHCSetNext_preserves_heapOrdered newNode
        (some heap[anchor]) (heap.set anchor newNode) nodes updated hheapOrdered
        hset
      simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
      rcases hrun with ⟨rfl, rfl⟩
      exact hactual

theorem pairVecDivVHCSetNext_preserves_cursorPrefixAbove
    (degreeLimit nodeIndex : Nat) (next : Option Nat)
    (nodes nodes' : Array PairVecDivVHCNode)
    (quotient divisor : SparsePolyZp)
    (hprefix : PairVecDivVHCCursorPrefixAbove degreeLimit nodes quotient divisor)
    (hrun : pairVecDivVHCSetNext nodeIndex next nodes = .ok nodes') :
    PairVecDivVHCCursorPrefixAbove degreeLimit nodes' quotient divisor := by
  unfold pairVecDivVHCSetNext at hrun
  split at hrun <;> try contradiction
  next hn =>
    simp only [Except.ok.injEq] at hrun
    subst nodes'
    exact hprefix.set_fields degreeLimit nodeIndex nodes quotient divisor
      nodes[nodeIndex] _ hn (Array.getElem?_eq_getElem hn) rfl rfl

theorem pairVecDivVHCInsert_get_ne
    (newNode : Nat) (heap heap' : Array Nat)
    (nodes nodes' : Array PairVecDivVHCNode)
    (hrun : pairVecDivVHCInsert newNode heap nodes = .ok (heap', nodes'))
    (i : Nat) (hne : newNode ≠ i) :
    nodes'[i]? = nodes[i]? := by
  rcases pairVecDivVHCInsert_nodes_result newNode heap heap' nodes nodes' hrun
    with ⟨next, hset⟩
  exact pairVecDivVHCSetNext_get_ne newNode next nodes nodes' hset i hne

theorem pairVecDivVHCInsert_preserves_cursorPrefixAbove
    (degreeLimit newNode : Nat) (heap heap' : Array Nat)
    (nodes nodes' : Array PairVecDivVHCNode)
    (quotient divisor : SparsePolyZp)
    (hprefix : PairVecDivVHCCursorPrefixAbove degreeLimit nodes quotient divisor)
    (hrun : pairVecDivVHCInsert newNode heap nodes = .ok (heap', nodes')) :
    PairVecDivVHCCursorPrefixAbove degreeLimit nodes' quotient divisor := by
  rcases pairVecDivVHCInsert_nodes_result newNode heap heap' nodes nodes' hrun
    with ⟨next, hset⟩
  exact pairVecDivVHCSetNext_preserves_cursorPrefixAbove degreeLimit newNode
    next nodes nodes' quotient divisor hprefix hset

theorem pairVecDivVHCInsert_preserves_divisorIndicesFixed
    (newNode : Nat) (heap heap' : Array Nat)
    (nodes nodes' : Array PairVecDivVHCNode)
    (hfixed : PairVecDivVHCNodeDivisorIndicesFixed nodes)
    (hrun : pairVecDivVHCInsert newNode heap nodes = .ok (heap', nodes')) :
    PairVecDivVHCNodeDivisorIndicesFixed nodes' := by
  rcases pairVecDivVHCInsert_nodes_result newNode heap heap' nodes nodes' hrun
    with ⟨next, hset⟩
  exact pairVecDivVHCSetNext_preserves_divisorIndicesFixed newNode next nodes
    nodes' hfixed hset

theorem pairVecDivVHCInsert_preserves_denotes
    (newNode : Nat) (heap heap' : Array Nat)
    (nodes nodes' : Array PairVecDivVHCNode)
    (quotient divisor : SparsePolyZp)
    (hdenotes : ∀ (i : Nat) (node : PairVecDivVHCNode),
      nodes[i]? = some node → node.mono ≠ none →
        PairVecDivVHCNodeDenotes quotient divisor node)
    (hrun : pairVecDivVHCInsert newNode heap nodes = .ok (heap', nodes')) :
    ∀ (i : Nat) (node : PairVecDivVHCNode),
      nodes'[i]? = some node → node.mono ≠ none →
        PairVecDivVHCNodeDenotes quotient divisor node := by
  rcases pairVecDivVHCInsert_nodes_result newNode heap heap' nodes nodes' hrun
    with ⟨next, hset⟩
  exact pairVecDivVHCSetNext_preserves_denotes newNode next nodes nodes'
    quotient divisor hdenotes hset

theorem pairVecDivVHCSetNext_preserves_resetReady
    (resetH quotientSize nodeIndex : Nat) (next : Option Nat)
    (nodes nodes' : Array PairVecDivVHCNode)
    (hready : PairVecDivVHCResetReady resetH quotientSize nodes)
    (houtside : resetH ≤ nodeIndex)
    (hrun : pairVecDivVHCSetNext nodeIndex next nodes = .ok nodes') :
    PairVecDivVHCResetReady resetH quotientSize nodes' := by
  unfold pairVecDivVHCSetNext at hrun
  split at hrun <;> try contradiction
  next hn =>
    simp only [Except.ok.injEq] at hrun
    subst nodes'
    exact PairVecDivVHCResetReady.shrink_set resetH resetH quotientSize
      nodeIndex nodes _ hn hready (Nat.le_refl _) houtside

theorem pairVecDivVHCInsert_preserves_resetReady
    (resetH quotientSize newNode : Nat) (heap heap' : Array Nat)
    (nodes nodes' : Array PairVecDivVHCNode)
    (hready : PairVecDivVHCResetReady resetH quotientSize nodes)
    (houtside : resetH ≤ newNode)
    (hrun : pairVecDivVHCInsert newNode heap nodes = .ok (heap', nodes')) :
    PairVecDivVHCResetReady resetH quotientSize nodes' := by
  rcases pairVecDivVHCInsert_nodes_result newNode heap heap' nodes nodes' hrun
    with ⟨next, hset⟩
  exact pairVecDivVHCSetNext_preserves_resetReady resetH quotientSize newNode
    next nodes nodes' hready houtside hset

theorem pairVecDivVHCInsert_empty_preserves_allActiveNodesBelow
    (degreeLimit newNode : Nat) (heap heap' : Array Nat)
    (nodes nodes' : Array PairVecDivVHCNode)
    (hempty : heap.size = 0)
    (hbelow : PairVecDivVHCAllActiveNodesBelow degreeLimit nodes)
    (hrun : pairVecDivVHCInsert newNode heap nodes = .ok (heap', nodes')) :
    PairVecDivVHCAllActiveNodesBelow degreeLimit nodes' := by
  unfold pairVecDivVHCInsert at hrun
  cases hmono : pairVecDivVHCMono newNode nodes with
  | error fault => simp [hmono] at hrun
  | ok mono =>
      simp only [hmono, hempty, ↓reduceDIte] at hrun
      cases hset : pairVecDivVHCSetNext newNode none nodes with
      | error fault => simp [hset] at hrun
      | ok updated =>
          rw [hset] at hrun
          simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
          rcases hrun with ⟨rfl, rfl⟩
          exact pairVecDivVHCSetNext_preserves_allActiveNodesBelow
            degreeLimit newNode none nodes updated hbelow hset

theorem pairVecDivVHCInsert_equalRoot_preserves_allActiveNodesBelow
    (degreeLimit newNode : Nat) (heap heap' : Array Nat)
    (nodes nodes' : Array PairVecDivVHCNode) (newMono rootMono : UMonomial)
    (hheap : 0 < heap.size)
    (hnew : pairVecDivVHCMono newNode nodes = .ok newMono)
    (hroot : pairVecDivVHCMono heap[0] nodes = .ok rootMono)
    (hequal : newMono.deg = rootMono.deg)
    (hbelow : PairVecDivVHCAllActiveNodesBelow degreeLimit nodes)
    (hrun : pairVecDivVHCInsert newNode heap nodes = .ok (heap', nodes')) :
    PairVecDivVHCAllActiveNodesBelow degreeLimit nodes' := by
  unfold pairVecDivVHCInsert at hrun
  simp only [hnew, Nat.ne_of_gt hheap, ↓reduceDIte, hroot, hequal] at hrun
  cases hset : pairVecDivVHCSetNext newNode (some heap[0]) nodes with
  | error fault => simp [hset] at hrun
  | ok updated =>
      rw [hset] at hrun
      simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
      rcases hrun with ⟨rfl, rfl⟩
      exact pairVecDivVHCSetNext_preserves_allActiveNodesBelow
        degreeLimit newNode (some heap[0]) nodes updated hbelow hset

theorem pairVecDivVHCInsert_newRoot_preserves_allActiveNodesBelow
    (degreeLimit newNode : Nat) (heap heap' : Array Nat)
    (nodes nodes' : Array PairVecDivVHCNode) (newMono rootMono : UMonomial)
    (hheap : 0 < heap.size)
    (hnew : pairVecDivVHCMono newNode nodes = .ok newMono)
    (hroot : pairVecDivVHCMono heap[0] nodes = .ok rootMono)
    (hnequal : newMono.deg ≠ rootMono.deg)
    (hgreater : newMono.deg > rootMono.deg)
    (hbelow : PairVecDivVHCAllActiveNodesBelow degreeLimit nodes)
    (hrun : pairVecDivVHCInsert newNode heap nodes = .ok (heap', nodes')) :
    PairVecDivVHCAllActiveNodesBelow degreeLimit nodes' := by
  unfold pairVecDivVHCInsert at hrun
  simp only [hnew, Nat.ne_of_gt hheap, ↓reduceDIte, hroot, hnequal,
    hgreater] at hrun
  cases hset : pairVecDivVHCSetNext newNode none nodes with
  | error fault => simp [hset] at hrun
  | ok updated =>
      rw [hset] at hrun
      cases hbubble : pairVecDivVHCBubble heap.size 0 newNode
          (heap.push newNode) with
      | error fault => simp [hbubble] at hrun
      | ok shifted =>
          rw [hbubble] at hrun
          simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
          rcases hrun with ⟨rfl, rfl⟩
          exact pairVecDivVHCSetNext_preserves_allActiveNodesBelow
            degreeLimit newNode none nodes updated hbelow hset

theorem pairVecDivVHCInsert_equalAnchor_preserves_allActiveNodesBelow
    (degreeLimit newNode anchor : Nat) (heap heap' : Array Nat)
    (nodes nodes' : Array PairVecDivVHCNode)
    (newMono rootMono anchorMono : UMonomial)
    (hheap : 0 < heap.size)
    (hnew : pairVecDivVHCMono newNode nodes = .ok newMono)
    (hroot : pairVecDivVHCMono heap[0] nodes = .ok rootMono)
    (hnequal : newMono.deg ≠ rootMono.deg)
    (hgreater : ¬ newMono.deg > rootMono.deg)
    (hanchor : pairVecDivVHCFindAnchor newMono.deg
      (pairVecDivVHCParent heap.size) heap nodes = .ok anchor)
    (ha : anchor < heap.size)
    (hanchorMono : pairVecDivVHCMono heap[anchor] nodes = .ok anchorMono)
    (hequalAnchor : newMono.deg = anchorMono.deg)
    (hbelow : PairVecDivVHCAllActiveNodesBelow degreeLimit nodes)
    (hrun : pairVecDivVHCInsert newNode heap nodes = .ok (heap', nodes')) :
    PairVecDivVHCAllActiveNodesBelow degreeLimit nodes' := by
  unfold pairVecDivVHCInsert at hrun
  simp only [hnew, Nat.ne_of_gt hheap, ↓reduceDIte, hroot, hnequal,
    hgreater, hanchor, ha, hanchorMono] at hrun
  simp only [hequalAnchor] at hrun
  cases hset : pairVecDivVHCSetNext newNode (some heap[anchor]) nodes with
  | error fault => simp [hset] at hrun
  | ok updated =>
      rw [hset] at hrun
      simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
      rcases hrun with ⟨rfl, rfl⟩
      exact pairVecDivVHCSetNext_preserves_allActiveNodesBelow
        degreeLimit newNode (some heap[anchor]) nodes nodes' hbelow hset

theorem pairVecDivVHCInsert_bubbleBelow_preserves_allActiveNodesBelow
    (degreeLimit newNode anchor : Nat) (heap heap' : Array Nat)
    (nodes nodes' : Array PairVecDivVHCNode)
    (newMono rootMono anchorMono : UMonomial)
    (hheap : 0 < heap.size)
    (hnew : pairVecDivVHCMono newNode nodes = .ok newMono)
    (hroot : pairVecDivVHCMono heap[0] nodes = .ok rootMono)
    (hnequal : newMono.deg ≠ rootMono.deg)
    (hgreater : ¬ newMono.deg > rootMono.deg)
    (hanchor : pairVecDivVHCFindAnchor newMono.deg
      (pairVecDivVHCParent heap.size) heap nodes = .ok anchor)
    (ha : anchor < heap.size)
    (hanchorMono : pairVecDivVHCMono heap[anchor] nodes = .ok anchorMono)
    (hnequalAnchor : newMono.deg ≠ anchorMono.deg)
    (hbelow : PairVecDivVHCAllActiveNodesBelow degreeLimit nodes)
    (hrun : pairVecDivVHCInsert newNode heap nodes = .ok (heap', nodes')) :
    PairVecDivVHCAllActiveNodesBelow degreeLimit nodes' := by
  unfold pairVecDivVHCInsert at hrun
  simp only [hnew, Nat.ne_of_gt hheap, ↓reduceDIte, hroot, hnequal,
    hgreater, hanchor, ha, hanchorMono] at hrun
  simp only [hnequalAnchor] at hrun
  cases hset : pairVecDivVHCSetNext newNode none nodes with
  | error fault => simp [hset] at hrun
  | ok updated =>
      rw [hset] at hrun
      cases hbubble : pairVecDivVHCBubbleBelow heap.size anchor newNode
          (heap.push newNode) with
      | error fault => simp [hbubble] at hrun
      | ok shifted =>
          rw [hbubble] at hrun
          simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
          rcases hrun with ⟨rfl, rfl⟩
          exact pairVecDivVHCSetNext_preserves_allActiveNodesBelow
            degreeLimit newNode none nodes nodes' hbelow hset

theorem pairVecDivVHCInsert_preserves_allActiveNodesBelow
    (degreeLimit newNode : Nat) (heap heap' : Array Nat)
    (nodes nodes' : Array PairVecDivVHCNode)
    (hbelow : PairVecDivVHCAllActiveNodesBelow degreeLimit nodes)
    (hrun : pairVecDivVHCInsert newNode heap nodes = .ok (heap', nodes')) :
    PairVecDivVHCAllActiveNodesBelow degreeLimit nodes' := by
  cases hnew : pairVecDivVHCMono newNode nodes with
  | error fault => simp [pairVecDivVHCInsert, hnew] at hrun
  | ok newMono =>
      by_cases hempty : heap.size = 0
      · exact pairVecDivVHCInsert_empty_preserves_allActiveNodesBelow
          degreeLimit newNode heap heap' nodes nodes' hempty hbelow hrun
      · have hheap : 0 < heap.size := Nat.pos_of_ne_zero hempty
        cases hroot : pairVecDivVHCMono heap[0] nodes with
        | error fault =>
            simp [pairVecDivVHCInsert, hnew, hempty, hroot] at hrun
        | ok rootMono =>
            by_cases hequal : newMono.deg = rootMono.deg
            · exact pairVecDivVHCInsert_equalRoot_preserves_allActiveNodesBelow
                degreeLimit newNode heap heap' nodes nodes' newMono rootMono
                hheap hnew hroot hequal hbelow hrun
            · by_cases hgreater : newMono.deg > rootMono.deg
              · exact pairVecDivVHCInsert_newRoot_preserves_allActiveNodesBelow
                  degreeLimit newNode heap heap' nodes nodes' newMono rootMono
                  hheap hnew hroot hequal hgreater hbelow hrun
              · cases hanchor : pairVecDivVHCFindAnchor newMono.deg
                    (pairVecDivVHCParent heap.size) heap nodes with
                | error fault =>
                    simp [pairVecDivVHCInsert, hnew, hempty, hroot, hequal,
                      hgreater, hanchor] at hrun
                | ok anchor =>
                    by_cases ha : anchor < heap.size
                    · cases hanchorMono : pairVecDivVHCMono heap[anchor] nodes with
                      | error fault =>
                          simp [pairVecDivVHCInsert, hnew, hempty, hroot,
                            hequal, hgreater, hanchor, ha, hanchorMono] at hrun
                      | ok anchorMono =>
                          by_cases hequalAnchor :
                              newMono.deg = anchorMono.deg
                          · exact
                              pairVecDivVHCInsert_equalAnchor_preserves_allActiveNodesBelow
                                degreeLimit newNode anchor heap heap' nodes
                                nodes' newMono rootMono anchorMono hheap hnew
                                hroot hequal hgreater hanchor ha hanchorMono
                                hequalAnchor hbelow hrun
                          · exact
                              pairVecDivVHCInsert_bubbleBelow_preserves_allActiveNodesBelow
                                degreeLimit newNode anchor heap heap' nodes
                                nodes' newMono rootMono anchorMono hheap hnew
                                hroot hequal hgreater hanchor ha hanchorMono
                                hequalAnchor hbelow hrun
                    · simp [pairVecDivVHCInsert, hnew, hempty, hroot, hequal,
                        hgreater, hanchor, ha] at hrun

/-- Complete heap-order preservation for every successful generated insertion
branch. -/
theorem pairVecDivVHCInsert_preserves_heapOrdered
    (newNode : Nat) (heap heap' : Array Nat)
    (nodes nodes' : Array PairVecDivVHCNode)
    (hvalid : PairVecDivVHCHeapPointersValid heap nodes)
    (hordered : PairVecDivVHCHeapOrdered heap nodes)
    (hrun : pairVecDivVHCInsert newNode heap nodes = .ok (heap', nodes')) :
    PairVecDivVHCHeapOrdered heap' nodes' := by
  cases hnew : pairVecDivVHCMono newNode nodes with
  | error fault => simp [pairVecDivVHCInsert, hnew] at hrun
  | ok newMono =>
      by_cases hempty : heap.size = 0
      · unfold pairVecDivVHCInsert at hrun
        simp only [hnew, hempty, ↓reduceDIte] at hrun
        cases hset : pairVecDivVHCSetNext newNode none nodes with
        | error fault => simp [hset] at hrun
        | ok updated =>
            rw [hset] at hrun
            simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
            rcases hrun with ⟨rfl, rfl⟩
            intro child parent hchild hparent hpos
            simp at hchild
            omega
      · have hheap : 0 < heap.size := Nat.pos_of_ne_zero hempty
        cases hroot : pairVecDivVHCMono heap[0] nodes with
        | error fault =>
            simp [pairVecDivVHCInsert, hnew, hempty, hroot] at hrun
        | ok rootMono =>
            by_cases hequal : newMono.deg = rootMono.deg
            · exact pairVecDivVHCInsert_equalRoot_preserves_heapOrdered
                newNode heap heap' nodes nodes' newMono rootMono hheap hordered
                hnew hroot hequal hrun
            · by_cases hgreater : newMono.deg > rootMono.deg
              · exact pairVecDivVHCInsert_newRoot_preserves_heapOrdered
                  newNode heap heap' nodes nodes' newMono rootMono hheap hvalid
                  hordered hnew hroot hgreater hrun
              · cases hanchor : pairVecDivVHCFindAnchor newMono.deg
                    (pairVecDivVHCParent heap.size) heap nodes with
                | error fault =>
                    simp [pairVecDivVHCInsert, hnew, hempty, hroot, hequal,
                      hgreater, hanchor] at hrun
                | ok anchor =>
                    by_cases ha : anchor < heap.size
                    · cases hanchorMono : pairVecDivVHCMono heap[anchor] nodes with
                      | error fault =>
                          simp [pairVecDivVHCInsert, hnew, hempty, hroot,
                            hequal, hgreater, hanchor, ha, hanchorMono] at hrun
                      | ok anchorMono =>
                          by_cases hequalAnchor :
                              newMono.deg = anchorMono.deg
                          · exact
                              pairVecDivVHCInsert_equalAnchor_preserves_heapOrdered
                                newNode anchor heap heap' nodes nodes' newMono
                                rootMono anchorMono hheap hordered hnew hroot
                                hequal hgreater hanchor ha hanchorMono
                                hequalAnchor hrun
                          · exact
                              pairVecDivVHCInsert_bubbleBelow_preserves_heapOrdered
                                newNode anchor heap heap' nodes nodes' newMono
                                rootMono anchorMono hheap hvalid hordered hnew
                                hroot hequal hgreater hanchor ha hanchorMono
                                hequalAnchor hrun
                    · simp [pairVecDivVHCInsert, hnew, hempty, hroot, hequal,
                        hgreater, hanchor, ha] at hrun

theorem pairVecDivVHCActivateInsert_preserves_heapOrdered
    (nodeIndex : Nat) (heap heap' : Array Nat)
    (nodes activated inserted : Array PairVecDivVHCNode)
    (quotient divisor : SparsePolyZp) (oldNode : PairVecDivVHCNode)
    (howned : PairVecDivVHCHeapChainsOwned heap nodes)
    (hordered : PairVecDivVHCHeapOrdered heap nodes)
    (hold : nodes[nodeIndex]? = some oldNode) (hinactive : oldNode.mono = none)
    (hactivate : pairVecDivVHCActivate nodeIndex nodes quotient divisor =
      .ok activated)
    (hinsert : pairVecDivVHCInsert nodeIndex heap activated =
      .ok (heap', inserted)) :
    PairVecDivVHCHeapOrdered heap' inserted := by
  rcases howned with ⟨owners, hownership⟩
  have hfresh := pairVecDivVHCHeapChainOwnership_fresh_of_mono_none heap owners
    nodes nodeIndex oldNode hownership hold hinactive
  have hactivatedOrdered :=
    pairVecDivVHCActivate_preserves_heapOrdered_of_freshHead nodeIndex heap
      nodes activated quotient divisor hfresh.1 hordered hactivate
  have hownershipActivated :
      PairVecDivVHCHeapChainOwnership heap owners activated := by
    refine ⟨?_, hownership.2.1, hownership.2.2⟩
    intro slot head hget
    have hchain := hownership.1 slot head hget
    exact pairVecDivVHCChainOwns_congr_on (some head) (owners head) nodes
      activated hchain (by
        intro i hi
        exact pairVecDivVHCActivate_get_ne nodeIndex nodes activated quotient
          divisor hactivate i (by
            intro heq
            subst i
            exact hfresh.2 slot head hget hi))
  exact pairVecDivVHCInsert_preserves_heapOrdered nodeIndex heap heap' activated
    inserted (hownershipActivated.heapPointersValid heap owners activated)
    hactivatedOrdered hinsert

theorem pairVecDivVHCActivateInsert_preserves_heapChainsHomogeneous
    (nodeIndex : Nat) (heap heap' : Array Nat)
    (nodes activated inserted : Array PairVecDivVHCNode)
    (quotient divisor : SparsePolyZp) (oldNode : PairVecDivVHCNode)
    (owners : Nat → Finset Nat)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hhomogeneous : PairVecDivVHCHeapChainsHomogeneous heap owners nodes)
    (hold : nodes[nodeIndex]? = some oldNode) (hinactive : oldNode.mono = none)
    (hactivate : pairVecDivVHCActivate nodeIndex nodes quotient divisor =
      .ok activated)
    (hinsert : pairVecDivVHCInsert nodeIndex heap activated =
      .ok (heap', inserted)) :
    ∃ owners', PairVecDivVHCHeapChainOwnership heap' owners' inserted ∧
      PairVecDivVHCHeapChainsHomogeneous heap' owners' inserted := by
  have hfresh := pairVecDivVHCHeapChainOwnership_fresh_of_mono_none heap owners
    nodes nodeIndex oldNode hownership hold hinactive
  have hownershipActivated :
      PairVecDivVHCHeapChainOwnership heap owners activated := by
    refine ⟨?_, hownership.2.1, hownership.2.2⟩
    intro slot head hget
    have hchain := hownership.1 slot head hget
    exact pairVecDivVHCChainOwns_congr_on (some head) (owners head) nodes
      activated hchain (by
        intro i hi
        exact pairVecDivVHCActivate_get_ne nodeIndex nodes activated quotient
          divisor hactivate i (by
            intro heq
            subst i
            exact hfresh.2 slot head hget hi))
  have hhomogeneousActivated :=
    pairVecDivVHCActivate_preserves_heapChainsHomogeneous_of_fresh nodeIndex
      heap owners nodes activated quotient divisor hhomogeneous hfresh.1
      hfresh.2 hactivate
  rcases pairVecDivVHCActivate_get nodeIndex nodes activated quotient divisor
      hactivate with ⟨hn, hq, hd, hnewGet⟩
  have holdEq : nodes[nodeIndex] = oldNode := by
    rw [Array.getElem?_eq_getElem hn] at hold
    exact Option.some.inj hold
  have hqOld : oldNode.quotientIndex < quotient.size := by
    simpa [holdEq] using hq
  have hdOld : oldNode.divisorIndex < divisor.size := by
    simpa [holdEq] using hd
  let newRecord : PairVecDivVHCNode := { oldNode with
    mono := some ⟨(quotient[oldNode.quotientIndex]'hqOld).1.deg +
      (divisor[oldNode.divisorIndex]'hdOld).1.deg⟩
    next := none }
  have hnewGet' : activated[nodeIndex]? = some newRecord := by
    simpa [newRecord, holdEq] using hnewGet
  have hnewMono : newRecord.mono = some
      ⟨(quotient[oldNode.quotientIndex]'hqOld).1.deg +
        (divisor[oldNode.divisorIndex]'hdOld).1.deg⟩ := rfl
  exact pairVecDivVHCInsert_preserves_heapChainsHomogeneous_of_fresh nodeIndex
    heap heap' activated inserted owners newRecord _ hownershipActivated
    hhomogeneousActivated hnewGet' hnewMono hfresh.1 hfresh.2 hinsert

/-- Every deferred source node occurs exactly once in `lin` and is already
active.  This is the proof state required by the literal reverse reinsertion
loop; it contains no execution budget or specification-level result. -/
def PairVecDivVHCLinReady (lin : Array Nat)
    (nodes : Array PairVecDivVHCNode) : Prop :=
  lin.toList.Nodup ∧
    ∀ nodeIndex ∈ lin.toList,
      ∃ node mono, nodes[nodeIndex]? = some node ∧ node.mono = some mono

theorem pairVecDivVHCCursorIndicesBounded_of_state
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode)
    (lin : Array Nat) (resetH : Nat) (quotient divisor : SparsePolyZp)
    (hstate : PairVecDivVHCStateCovered heap nodes lin resetH)
    (hlinReady : PairVecDivVHCLinReady lin nodes)
    (hresetReady : PairVecDivVHCResetReady resetH quotient.size nodes)
    (hdenotes : ∀ (i : Nat) (node : PairVecDivVHCNode),
      nodes[i]? = some node → node.mono ≠ none →
        PairVecDivVHCNodeDenotes quotient divisor node) :
    PairVecDivVHCCursorIndicesBounded quotient.size nodes := by
  rcases hstate with ⟨owners, hownership, hcovered⟩
  intro i node hget
  have hi : i < nodes.size := by
    by_contra hnot
    rw [Array.getElem?_eq_none (by omega)] at hget
    contradiction
  rcases hcovered i hi with hreset | hlin | ⟨slot, head, hheap, hmem⟩
  · rcases hresetReady.2 i hreset with
      ⟨readyNode, hreadyGet, hquotientIndex, _, _⟩
    rw [hget] at hreadyGet
    simp only [Option.some.injEq] at hreadyGet
    subst readyNode
    exact hquotientIndex.le
  · rcases hlinReady.2 i (by simpa using hlin) with
      ⟨activeNode, mono, hactiveGet, hactiveMono⟩
    rw [hget] at hactiveGet
    simp only [Option.some.injEq] at hactiveGet
    subst activeNode
    exact (hdenotes i node hget (by simp [hactiveMono])).quotientIndex_lt.le
  · have howns := hownership.1 slot head hheap
    rcases pairVecDivVHCChainOwns_mem_active (some head) (owners head) nodes
        howns i hmem with ⟨activeNode, mono, hactiveGet, hactiveMono⟩
    rw [hget] at hactiveGet
    simp only [Option.some.injEq] at hactiveGet
    subst activeNode
    exact (hdenotes i node hget (by simp [hactiveMono])).quotientIndex_lt.le

theorem PairVecDivVHCLinReady.set_outside
    (lin : Array Nat) (nodes : Array PairVecDivVHCNode)
    (nodeIndex : Nat) (updated : PairVecDivVHCNode)
    (hn : nodeIndex < nodes.size)
    (houtside : nodeIndex ∉ lin.toList)
    (hready : PairVecDivVHCLinReady lin nodes) :
    PairVecDivVHCLinReady lin (nodes.set nodeIndex updated) := by
  refine ⟨hready.1, ?_⟩
  intro i hmem
  rcases hready.2 i hmem with ⟨node, mono, hnode, hmono⟩
  refine ⟨node, mono, ?_, hmono⟩
  rw [Array.getElem?_set_ne hn (by
    intro heq
    subst i
    exact houtside hmem)]
  exact hnode

theorem PairVecDivVHCLinReady.push_set
    (lin : Array Nat) (nodes : Array PairVecDivVHCNode)
    (nodeIndex : Nat) (updated : PairVecDivVHCNode) (mono : UMonomial)
    (hn : nodeIndex < nodes.size)
    (houtside : nodeIndex ∉ lin.toList)
    (hready : PairVecDivVHCLinReady lin nodes)
    (hmono : updated.mono = some mono) :
    PairVecDivVHCLinReady (lin.push nodeIndex)
      (nodes.set nodeIndex updated) := by
  refine ⟨?_, ?_⟩
  · simp only [Array.toList_push, List.nodup_append, hready.1,
      List.nodup_singleton, true_and, List.mem_singleton]
    intro a ha b hb
    subst b
    intro heq
    subst a
    exact houtside ha
  · intro i hmem
    simp only [Array.toList_push, List.mem_append, List.mem_singleton] at hmem
    rcases hmem with hmem | rfl
    · rcases hready.2 i hmem with ⟨node, oldMono, hnode, holdMono⟩
      refine ⟨node, oldMono, ?_, holdMono⟩
      rw [Array.getElem?_set_ne hn (by
        intro heq
        subst i
        exact houtside hmem)]
      exact hnode
    · exact ⟨updated, mono, Array.getElem?_set_self hn, hmono⟩

theorem pairVecDivVHCConsumeNode_preserves_linReady
    (this : DenseUPolyZp) (nodeIndex : Nat) (currentMono : UMonomial)
    (k k' : UInt64) (nodes nodes' : Array PairVecDivVHCNode)
    (lin lin' : Array Nat) (resetH resetH' : Nat) (next : Option Nat)
    (quotient divisor : SparsePolyZp)
    (hn : nodeIndex < nodes.size)
    (hactive : nodes[nodeIndex].mono = some currentMono)
    (houtside : nodeIndex ∉ lin.toList)
    (hready : PairVecDivVHCLinReady lin nodes)
    (hrun : pairVecDivVHCConsumeNode this nodeIndex k nodes lin resetH
      quotient divisor = .ok (k', nodes', lin', resetH', next)) :
    PairVecDivVHCLinReady lin' nodes' ∧
      lin'.toList.toFinset ⊆ insert nodeIndex lin.toList.toFinset := by
  unfold pairVecDivVHCConsumeNode at hrun
  simp only [hn, ↓reduceDIte] at hrun
  split at hrun <;> try contradiction
  next hq =>
    split at hrun <;> try contradiction
    next hd =>
      split at hrun
      next hadvance =>
        simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
        rcases hrun with ⟨rfl, rfl, rfl, rfl, rfl⟩
        refine ⟨PairVecDivVHCLinReady.push_set lin nodes nodeIndex _ _ hn
          houtside hready rfl, ?_⟩
        intro i hmem
        simp only [Array.toList_push, List.mem_toFinset, List.mem_append,
          List.mem_singleton, Finset.mem_insert] at hmem ⊢
        exact hmem.symm
      next hadvance =>
        split at hrun <;> try contradiction
        next hexhausted =>
          split at hrun <;> try contradiction
          next horder =>
            simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
            rcases hrun with ⟨rfl, rfl, rfl, rfl, rfl⟩
            refine ⟨PairVecDivVHCLinReady.set_outside lin nodes nodeIndex _
              hn houtside hready, ?_⟩
            intro i hmem
            exact Finset.mem_insert_of_mem hmem

theorem pairVecDivVHCConsumeNode_next
    (this : DenseUPolyZp) (nodeIndex : Nat) (k k' : UInt64)
    (nodes nodes' : Array PairVecDivVHCNode) (lin lin' : Array Nat)
    (resetH resetH' : Nat) (next : Option Nat)
    (quotient divisor : SparsePolyZp) (hn : nodeIndex < nodes.size)
    (hrun : pairVecDivVHCConsumeNode this nodeIndex k nodes lin resetH
      quotient divisor = .ok (k', nodes', lin', resetH', next)) :
    next = nodes[nodeIndex].next := by
  unfold pairVecDivVHCConsumeNode at hrun
  simp only [hn, ↓reduceDIte] at hrun
  split at hrun <;> try contradiction
  next hq =>
    split at hrun <;> try contradiction
    next hd =>
      split at hrun
      next hadvance =>
        simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
        exact hrun.2.2.2.2.symm
      next hadvance =>
        split at hrun <;> try contradiction
        next hexhausted =>
          split at hrun <;> try contradiction
          next horder =>
            simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
            exact hrun.2.2.2.2.symm

theorem pairVecDivVHCConsumeNode_preserves_chainAtDegree_tail
    (this : DenseUPolyZp) (nodeIndex : Nat) (unvisited : Finset Nat)
    (degree : Nat) (k k' : UInt64)
    (nodes nodes' : Array PairVecDivVHCNode) (lin lin' : Array Nat)
    (resetH resetH' : Nat) (next : Option Nat)
    (quotient divisor : SparsePolyZp) (hn : nodeIndex < nodes.size)
    (htail : PairVecDivVHCChainAtDegree nodes[nodeIndex].next
      (unvisited.erase nodeIndex) nodes degree)
    (hrun : pairVecDivVHCConsumeNode this nodeIndex k nodes lin resetH
      quotient divisor = .ok (k', nodes', lin', resetH', next)) :
    PairVecDivVHCChainAtDegree next (unvisited.erase nodeIndex) nodes'
      degree := by
  have hnext := pairVecDivVHCConsumeNode_next this nodeIndex k k' nodes nodes'
    lin lin' resetH resetH' next quotient divisor hn hrun
  rw [hnext]
  exact pairVecDivVHCChainAtDegree_congr_on nodes[nodeIndex].next
    (unvisited.erase nodeIndex) nodes nodes' degree htail (by
      intro i hi
      exact pairVecDivVHCConsumeNode_get_ne this nodeIndex k k' nodes nodes'
        lin lin' resetH resetH' next quotient divisor hrun i
          (Finset.mem_erase.mp hi).1.symm)

def PairVecDivVHCStoredProductAtDegree (degree : Nat)
    (quotient divisor : SparsePolyZp) (product : UInt64 × UInt64) : Prop :=
  ∃ quotientTerm ∈ quotient.toList, ∃ divisorTerm ∈ divisor.toList,
    quotientTerm.2.val = product.1 ∧ divisorTerm.2.val = product.2 ∧
      quotientTerm.1.deg + divisorTerm.1.deg = degree

/-- The concrete coefficient product currently named by one cursor node.
Invalid storage contributes zero; refinement invariants rule that case out for
owned active nodes. -/
def pairVecDivVHCNodeProductValue (p : Nat)
    (nodes : Array PairVecDivVHCNode) (quotient divisor : SparsePolyZp)
    (nodeIndex : Nat) : ZMod p :=
  match nodes[nodeIndex]? with
  | none => 0
  | some node =>
      match quotient[node.quotientIndex]?, divisor[node.divisorIndex]? with
      | some quotientTerm, some divisorTerm =>
          (quotientTerm.2.val.toNat : ZMod p) *
            (divisorTerm.2.val.toNat : ZMod p)
      | _, _ => 0

def PairVecDivVHCPairAtDegree (degree : Nat)
    (quotient divisor : SparsePolyZp) (pair : Nat × Nat) : Prop :=
  match quotient[pair.1]?, divisor[pair.2]? with
  | some quotientTerm, some divisorTerm =>
      quotientTerm.1.deg + divisorTerm.1.deg = degree
  | _, _ => False

instance (degree : Nat) (quotient divisor : SparsePolyZp) (pair : Nat × Nat) :
    Decidable (PairVecDivVHCPairAtDegree degree quotient divisor pair) := by
  unfold PairVecDivVHCPairAtDegree
  split <;> infer_instance

/-- Concrete indexed quotient/divisor-tail pairs whose monomial degrees hit
the current source frontier.  Indices, rather than term values, are retained
so equal coefficients in different rows remain distinct contributions. -/
def PairVecDivVHCTargetPairsAtDegree (degree : Nat)
    (quotient divisor : SparsePolyZp) : Finset (Nat × Nat) :=
  ((Finset.range quotient.size) ×ˢ (Finset.Ico 1 divisor.size)).filter
    fun pair => PairVecDivVHCPairAtDegree degree quotient divisor pair

/-- Coefficient contribution of one concrete indexed source pair. -/
def pairVecDivVHCIndexedPairProductValue (p : Nat)
    (quotient divisor : SparsePolyZp) (pair : Nat × Nat) : ZMod p :=
  match quotient[pair.1]?, divisor[pair.2]? with
  | some quotientTerm, some divisorTerm =>
      (quotientTerm.2.val.toNat : ZMod p) *
        (divisorTerm.2.val.toNat : ZMod p)
  | _, _ => 0

/-- Source-index pair stored in a cursor node; invalid indices receive the
irrelevant default pair and are excluded by ownership invariants. -/
def pairVecDivVHCNodeSourcePair (nodes : Array PairVecDivVHCNode)
    (nodeIndex : Nat) : Nat × Nat :=
  match nodes[nodeIndex]? with
  | some node => (node.quotientIndex, node.divisorIndex)
  | none => (0, 0)

theorem pairVecDivVHCTargetPairsAtDegree_mem_iff
    (degree q d : Nat) (quotient divisor : SparsePolyZp)
    (quotientTerm divisorTerm : UMonomial × Zp)
    (hquotient : quotient[q]? = some quotientTerm)
    (hdivisor : divisor[d]? = some divisorTerm) :
    (q, d) ∈ PairVecDivVHCTargetPairsAtDegree degree quotient divisor ↔
      0 < d ∧
        quotientTerm.1.deg + divisorTerm.1.deg = degree := by
  have hq : q < quotient.size := by
    by_contra hnot
    rw [Array.getElem?_eq_none (by omega)] at hquotient
    contradiction
  have hd : d < divisor.size := by
    by_contra hnot
    rw [Array.getElem?_eq_none (by omega)] at hdivisor
    contradiction
  rw [Array.getElem?_eq_getElem hq] at hquotient
  rw [Array.getElem?_eq_getElem hd] at hdivisor
  simp only [Option.some.injEq] at hquotient hdivisor
  subst quotientTerm
  subst divisorTerm
  simp [PairVecDivVHCTargetPairsAtDegree, PairVecDivVHCPairAtDegree,
    hq, hd] <;> omega

theorem pairVecDivVHCIndexedPairProductValue_of_get
    (p q d : Nat) (quotient divisor : SparsePolyZp)
    (quotientTerm divisorTerm : UMonomial × Zp)
    (hquotient : quotient[q]? = some quotientTerm)
    (hdivisor : divisor[d]? = some divisorTerm) :
    pairVecDivVHCIndexedPairProductValue p quotient divisor (q, d) =
      (quotientTerm.2.val.toNat : ZMod p) *
        (divisorTerm.2.val.toNat : ZMod p) := by
  simp [pairVecDivVHCIndexedPairProductValue, hquotient, hdivisor]

/-- Conversely, every heap-owned node at `degree` stores a genuine indexed
quotient/divisor-tail pair at that degree.  The fixed allocation identity
also recovers the node index from the divisor index, giving the inverse map
needed for a multiplicity-preserving finite-sum bijection. -/
theorem pairVecDivVHCOwnedNode_sourcePairAtDegree
    (p degree i : Nat) (quotient divisor : SparsePolyZp)
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode)
    (owners : Nat → Finset Nat)
    (hfixed : PairVecDivVHCNodeDivisorIndicesFixed nodes)
    (hdenotes : ∀ (j : Nat) (node : PairVecDivVHCNode),
      nodes[j]? = some node → node.mono ≠ none →
        PairVecDivVHCNodeDenotes quotient divisor node)
    (hi : i ∈ PairVecDivVHCHeapOwnedNodesAtDegree degree heap owners nodes) :
    pairVecDivVHCNodeSourcePair nodes i ∈
        PairVecDivVHCTargetPairsAtDegree degree quotient divisor ∧
      (pairVecDivVHCNodeSourcePair nodes i).2 - 1 = i ∧
      pairVecDivVHCIndexedPairProductValue p quotient divisor
          (pairVecDivVHCNodeSourcePair nodes i) =
        pairVecDivVHCNodeProductValue p nodes quotient divisor i := by
  simp only [PairVecDivVHCHeapOwnedNodesAtDegree, Finset.mem_filter,
    PairVecDivVHCHeapOwnedNodes, Finset.mem_biUnion, List.mem_toFinset,
    pairVecDivVHCNodeAtDegree_iff] at hi
  rcases hi with ⟨⟨_, _, _⟩, mono, hmonoRun, hmonoDegree⟩
  rcases (pairVecDivVHCMono_eq_ok_iff i nodes mono).1 hmonoRun with
    ⟨node, hnode, hnodeMono⟩
  rcases hdenotes i node hnode (by simp [hnodeMono]) with
    ⟨quotientTerm, divisorTerm, hquotient, hdivisor, hstoredMono⟩
  have hdivisorIndex := hfixed i node hnode
  have hdpos : 0 < node.divisorIndex := by omega
  have hsourceDegree :
      quotientTerm.1.deg + divisorTerm.1.deg = degree := by
    rw [hnodeMono] at hstoredMono
    have hdegrees := congrArg UMonomial.deg (Option.some.inj hstoredMono)
    exact hdegrees.symm.trans hmonoDegree
  have hpairMem :=
    (pairVecDivVHCTargetPairsAtDegree_mem_iff degree node.quotientIndex
      node.divisorIndex quotient divisor quotientTerm divisorTerm hquotient
      hdivisor).2 ⟨hdpos, hsourceDegree⟩
  have hsourcePair : pairVecDivVHCNodeSourcePair nodes i =
      (node.quotientIndex, node.divisorIndex) := by
    simp [pairVecDivVHCNodeSourcePair, hnode]
  rw [hsourcePair]
  refine ⟨hpairMem, by omega, ?_⟩
  rw [pairVecDivVHCIndexedPairProductValue_of_get p node.quotientIndex
    node.divisorIndex quotient divisor quotientTerm divisorTerm hquotient
    hdivisor]
  simp [pairVecDivVHCNodeProductValue, hnode, hquotient, hdivisor]

/-- The reverse half of the indexed bijection: a target pair maps to its
fixed divisor-row node, and reading that node recovers the original pair and
the identical coefficient contribution. -/
theorem pairVecDivVHCTargetPair_ownedNode
    (p degreeLimit dividendIndex resetH : Nat)
    (dividend quotient divisor : SparsePolyZp)
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode)
    (owners : Nat → Finset Nat) (frontier : PairVecDivVHCFrontier)
    (hsize : nodes.size = divisor.size - 1)
    (hfixed : PairVecDivVHCNodeDivisorIndicesFixed nodes)
    (hstate : PairVecDivVHCStateCovered heap nodes #[] resetH)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hhomogeneous : PairVecDivVHCHeapChainsHomogeneous heap owners nodes)
    (hresetReady : PairVecDivVHCResetReady resetH quotient.size nodes)
    (hordered : PairVecDivVHCHeapOrdered heap nodes)
    (hdenotes : ∀ (i : Nat) (node : PairVecDivVHCNode),
      nodes[i]? = some node → node.mono ≠ none →
        PairVecDivVHCNodeDenotes quotient divisor node)
    (hcanonical : SparsePolyZp.Canonical p quotient)
    (hprefix : PairVecDivVHCCursorPrefixAbove degreeLimit nodes quotient divisor)
    (hfrontier : frontier.degree < degreeLimit)
    (hselect : pairVecDivVHCSelectFrontier dividendIndex dividend heap nodes =
      .ok frontier)
    (pair : Nat × Nat)
    (hpair : pair ∈ PairVecDivVHCTargetPairsAtDegree frontier.degree
      quotient divisor) :
    pair.2 - 1 ∈ PairVecDivVHCHeapOwnedNodesAtDegree frontier.degree heap
        owners nodes ∧
      pairVecDivVHCNodeSourcePair nodes (pair.2 - 1) = pair ∧
      pairVecDivVHCNodeProductValue p nodes quotient divisor (pair.2 - 1) =
        pairVecDivVHCIndexedPairProductValue p quotient divisor pair := by
  rcases pair with ⟨q, d⟩
  rw [PairVecDivVHCTargetPairsAtDegree] at hpair
  rcases Finset.mem_filter.mp hpair with ⟨hindices, hdegree⟩
  rcases Finset.mem_product.mp hindices with ⟨hq, hd⟩
  rw [Finset.mem_range] at hq
  rw [Finset.mem_Ico] at hd
  rcases hd with ⟨hdpos, hd⟩
  have hquotient : quotient[q]? = some quotient[q] :=
    Array.getElem?_eq_getElem hq
  have hdivisor : divisor[d]? = some divisor[d] :=
    Array.getElem?_eq_getElem hd
  simp [PairVecDivVHCPairAtDegree, hquotient, hdivisor] at hdegree
  rcases pairVecDivVHCFrontierProduct_owned_cursor p degreeLimit dividendIndex
      resetH d q dividend quotient divisor heap nodes owners frontier
      quotient[q] divisor[d] hsize hfixed hstate hownership hhomogeneous
      hresetReady hordered hdenotes hcanonical hprefix hfrontier hdpos hd
      hquotient hdivisor hdegree hselect with
    ⟨slot, head, node, hheap, hnode, hmem, hnodeD, hcursor⟩
  have howned := pairVecDivVHCFrontierProduct_mem_ownedNodesAtDegree p
    degreeLimit dividendIndex resetH d q dividend quotient divisor heap nodes
    owners frontier quotient[q] divisor[d] hsize hfixed hstate hownership
    hhomogeneous hresetReady hordered hdenotes hcanonical hprefix hfrontier
    hdpos hd hquotient hdivisor hdegree hselect
  refine ⟨howned, ?_, ?_⟩
  · simp [pairVecDivVHCNodeSourcePair, hnode, ← hcursor, hnodeD]
  · simp [pairVecDivVHCNodeProductValue,
      pairVecDivVHCIndexedPairProductValue, hnode, ← hcursor, hnodeD,
      hquotient, hdivisor]

/-- Exact change of variables from generated-C++ heap owner nodes to the L2
indexed quotient/divisor-tail pairs.  Both inverse laws are explicit, so this
is a true bijective sum rewrite rather than a set-level coverage argument. -/
theorem pairVecDivVHCHeapOwnerSum_eq_targetPairSum
    (p degreeLimit dividendIndex resetH : Nat)
    (dividend quotient divisor : SparsePolyZp)
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode)
    (owners : Nat → Finset Nat) (frontier : PairVecDivVHCFrontier)
    (hsize : nodes.size = divisor.size - 1)
    (hfixed : PairVecDivVHCNodeDivisorIndicesFixed nodes)
    (hstate : PairVecDivVHCStateCovered heap nodes #[] resetH)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hhomogeneous : PairVecDivVHCHeapChainsHomogeneous heap owners nodes)
    (hresetReady : PairVecDivVHCResetReady resetH quotient.size nodes)
    (hordered : PairVecDivVHCHeapOrdered heap nodes)
    (hdenotes : ∀ (i : Nat) (node : PairVecDivVHCNode),
      nodes[i]? = some node → node.mono ≠ none →
        PairVecDivVHCNodeDenotes quotient divisor node)
    (hcanonical : SparsePolyZp.Canonical p quotient)
    (hprefix : PairVecDivVHCCursorPrefixAbove degreeLimit nodes quotient divisor)
    (hfrontier : frontier.degree < degreeLimit)
    (hselect : pairVecDivVHCSelectFrontier dividendIndex dividend heap nodes =
      .ok frontier) :
    (∑ i ∈ PairVecDivVHCHeapOwnedNodesAtDegree frontier.degree heap owners
        nodes, pairVecDivVHCNodeProductValue p nodes quotient divisor i) =
      ∑ pair ∈ PairVecDivVHCTargetPairsAtDegree frontier.degree quotient
        divisor, pairVecDivVHCIndexedPairProductValue p quotient divisor pair := by
  symm
  refine Finset.sum_bij'
    (s := PairVecDivVHCTargetPairsAtDegree frontier.degree quotient divisor)
    (t := PairVecDivVHCHeapOwnedNodesAtDegree frontier.degree heap owners nodes)
    (f := pairVecDivVHCIndexedPairProductValue p quotient divisor)
    (g := pairVecDivVHCNodeProductValue p nodes quotient divisor)
    (fun pair _ => pair.2 - 1)
    (fun i _ => pairVecDivVHCNodeSourcePair nodes i)
    ?_ ?_ ?_ ?_ ?_
  · intro pair hpair
    exact (pairVecDivVHCTargetPair_ownedNode p degreeLimit dividendIndex resetH
      dividend quotient divisor heap nodes owners frontier hsize hfixed hstate
      hownership hhomogeneous hresetReady hordered hdenotes hcanonical hprefix
      hfrontier hselect pair hpair).1
  · intro i hi
    exact (pairVecDivVHCOwnedNode_sourcePairAtDegree p frontier.degree i
      quotient divisor heap nodes owners hfixed hdenotes hi).1
  · intro pair hpair
    exact (pairVecDivVHCTargetPair_ownedNode p degreeLimit dividendIndex resetH
      dividend quotient divisor heap nodes owners frontier hsize hfixed hstate
      hownership hhomogeneous hresetReady hordered hdenotes hcanonical hprefix
      hfrontier hselect pair hpair).2.1
  · intro i hi
    exact (pairVecDivVHCOwnedNode_sourcePairAtDegree p frontier.degree i
      quotient divisor heap nodes owners hfixed hdenotes hi).2.1
  · intro pair hpair
    exact (pairVecDivVHCTargetPair_ownedNode p degreeLimit dividendIndex resetH
      dividend quotient divisor heap nodes owners frontier hsize hfixed hstate
      hownership hhomogeneous hresetReady hordered hdenotes hcanonical hprefix
      hfrontier hselect pair hpair).2.2.symm

def pairVecDivVHCListProductCoeffValue (p degree : Nat) :
    List (UMonomial × Zp) → List (UMonomial × Zp) → ZMod p
  | [], _ => 0
  | quotientTerm :: quotientTerms, divisorTerms =>
      divisorTerms.foldr (fun divisorTerm sum =>
        if quotientTerm.1.deg + divisorTerm.1.deg = degree then
          (quotientTerm.2.val.toNat : ZMod p) *
              (divisorTerm.2.val.toNat : ZMod p) + sum
        else sum) 0 +
        pairVecDivVHCListProductCoeffValue p degree quotientTerms divisorTerms

theorem pairVecDivVHCListProductCoeffValue_row_eq_finSum
    (p degree : Nat) (quotientTerm : UMonomial × Zp)
    (divisorTerms : List (UMonomial × Zp)) :
    divisorTerms.foldr (fun divisorTerm sum =>
        if quotientTerm.1.deg + divisorTerm.1.deg = degree then
          (quotientTerm.2.val.toNat : ZMod p) *
              (divisorTerm.2.val.toNat : ZMod p) + sum
        else sum) 0 =
      ∑ d : Fin divisorTerms.length,
        if quotientTerm.1.deg + divisorTerms[d].1.deg = degree then
          (quotientTerm.2.val.toNat : ZMod p) *
            (divisorTerms[d].2.val.toNat : ZMod p)
        else 0 := by
  rw [← List.sum_ofFn]
  change _ = (List.ofFn ((fun divisorTerm =>
    if quotientTerm.1.deg + divisorTerm.1.deg = degree then
      (quotientTerm.2.val.toNat : ZMod p) *
        (divisorTerm.2.val.toNat : ZMod p)
    else 0) ∘ divisorTerms.get)).sum
  rw [← List.map_ofFn, List.ofFn_get]
  induction divisorTerms with
  | nil => simp
  | cons divisorTerm divisorTerms ih =>
      by_cases hdegree :
          quotientTerm.1.deg + divisorTerm.1.deg = degree
      · simp [hdegree, ih]
      · simp [hdegree, ih]

theorem pairVecDivVHCListProductCoeffValue_eq_finSum
    (p degree : Nat) (quotientTerms divisorTerms : List (UMonomial × Zp)) :
    pairVecDivVHCListProductCoeffValue p degree quotientTerms divisorTerms =
      ∑ q : Fin quotientTerms.length, ∑ d : Fin divisorTerms.length,
        if quotientTerms[q].1.deg + divisorTerms[d].1.deg = degree then
          (quotientTerms[q].2.val.toNat : ZMod p) *
            (divisorTerms[d].2.val.toNat : ZMod p)
        else 0 := by
  rw [← List.sum_ofFn]
  change _ = (List.ofFn ((fun quotientTerm =>
    ∑ d : Fin divisorTerms.length,
      if quotientTerm.1.deg + divisorTerms[d].1.deg = degree then
        (quotientTerm.2.val.toNat : ZMod p) *
          (divisorTerms[d].2.val.toNat : ZMod p)
      else 0) ∘ quotientTerms.get)).sum
  rw [← List.map_ofFn, List.ofFn_get]
  induction quotientTerms with
  | nil => simp [pairVecDivVHCListProductCoeffValue]
  | cons quotientTerm quotientTerms ih =>
      simp [pairVecDivVHCListProductCoeffValue, ih,
        pairVecDivVHCListProductCoeffValue_row_eq_finSum]

theorem pairVecDivVHCTargetPairSum_eq_listProductCoeffValue
    (p degree : Nat) (quotient divisor : SparsePolyZp) :
    (∑ pair ∈ PairVecDivVHCTargetPairsAtDegree degree quotient divisor,
        pairVecDivVHCIndexedPairProductValue p quotient divisor pair) =
      pairVecDivVHCListProductCoeffValue p degree quotient.toList
        divisor.toList.tail := by
  rw [pairVecDivVHCListProductCoeffValue_eq_finSum]
  rw [PairVecDivVHCTargetPairsAtDegree, Finset.sum_filter]
  rw [Finset.sum_product (Finset.range quotient.size)
    (Finset.Ico 1 divisor.size) (fun pair =>
      if PairVecDivVHCPairAtDegree degree quotient divisor pair then
        pairVecDivVHCIndexedPairProductValue p quotient divisor pair
      else 0)]
  rw [Finset.sum_fin_eq_sum_range]
  apply Finset.sum_congr (by simp)
  intro q hqMem
  have hq : q < quotient.size := Finset.mem_range.mp hqMem
  simp only [Array.length_toList, hq, dif_pos]
  have hfinRows :
      (∑ d : Fin divisor.toList.tail.length,
        if quotient[q].1.deg +
            divisor.toList.tail[d].1.deg = degree then
          (quotient[q].2.val.toNat : ZMod p) *
            (divisor.toList.tail[d].2.val.toNat : ZMod p)
        else 0) =
      ∑ d : Fin divisor.toList.tail.length,
        if quotient[q].1.deg +
            (divisor[d.val + 1]'(by
              have hdlt := d.isLt
              simp only [List.length_tail, Array.length_toList] at hdlt
              omega)).1.deg = degree then
          (quotient[q].2.val.toNat : ZMod p) *
            ((divisor[d.val + 1]'(by
              have hdlt := d.isLt
              simp only [List.length_tail, Array.length_toList] at hdlt
              omega)).2.val.toNat : ZMod p)
        else 0 := by
    apply Finset.sum_congr rfl
    intro d _
    have hd : d.val + 1 < divisor.size := by
      have hdlt := d.isLt
      simp only [List.length_tail, Array.length_toList] at hdlt
      omega
    have htailGet : divisor.toList.tail[d] = divisor[d.val + 1] := by
      change divisor.toList.tail[d.val] = divisor[d.val + 1]
      rw [List.getElem_tail d.isLt, Array.getElem_toList hd]
    rw [htailGet]
  let qFin : Fin quotient.toList.length := ⟨q, by simpa using hq⟩
  have hquotientGet : quotient.toList[qFin] = quotient[q] := by
    change quotient.toList[qFin.val] = quotient[q]
    exact Array.getElem_toList hq
  rw [hquotientGet]
  rw [hfinRows]
  refine Finset.sum_bij'
    (s := Finset.Ico 1 divisor.size)
    (t := (Finset.univ : Finset (Fin divisor.toList.tail.length)))
    (f := fun d =>
      if PairVecDivVHCPairAtDegree degree quotient divisor (q, d) then
        pairVecDivVHCIndexedPairProductValue p quotient divisor (q, d)
      else 0)
    (g := fun d =>
      if quotient[q].1.deg +
          (divisor[d.val + 1]'(by
            have hdlt := d.isLt
            simp only [List.length_tail, Array.length_toList] at hdlt
            omega)).1.deg = degree then
        (quotient[q].2.val.toNat : ZMod p) *
          ((divisor[d.val + 1]'(by
            have hdlt := d.isLt
            simp only [List.length_tail, Array.length_toList] at hdlt
            omega)).2.val.toNat : ZMod p)
      else 0)
    (fun d hd => ⟨d - 1, by
      simp only [Finset.mem_Ico] at hd
      simp only [List.length_tail, Array.length_toList]
      omega⟩)
    (fun d _ => d.val + 1)
    ?_ ?_ ?_ ?_ ?_
  · intro d hd
    simp
  · intro d hd
    simp only [Finset.mem_univ, Finset.mem_Ico]
    have hlength : divisor.toList.tail.length = divisor.size - 1 := by simp
    have hdBound := d.isLt
    omega
  · intro d hd
    simp only [Finset.mem_Ico] at hd
    exact Nat.sub_add_cancel hd.1
  · intro d hd
    apply Fin.ext
    simp
  · intro d hd
    simp only [Finset.mem_Ico] at hd
    have hdpos : 0 < d := by omega
    have hdbound : d < divisor.size := by omega
    have hquotient := Array.getElem?_eq_getElem hq
    have hdivisor := Array.getElem?_eq_getElem hdbound
    have hsucc : d - 1 + 1 = d := by omega
    have hsourceBound : d - 1 + 1 < divisor.size := by omega
    have hsourceGet : divisor[d - 1 + 1]'hsourceBound = divisor[d] :=
      getElem_congr rfl hsucc _
    dsimp
    rw [hsourceGet]
    simp [PairVecDivVHCPairAtDegree,
      pairVecDivVHCIndexedPairProductValue, hquotient, hdivisor]

theorem pairVecDivVHCListProductCoeffValue_row (p degree : Nat)
    (quotientMono : UMonomial) (quotientCoefficient : Zp)
    (divisorTerms : List (UMonomial × Zp)) :
    divisorTerms.foldr (fun divisorTerm sum =>
        if quotientMono.deg + divisorTerm.1.deg = degree then
          (quotientCoefficient.val.toNat : ZMod p) *
              (divisorTerm.2.val.toNat : ZMod p) + sum
        else sum) 0 =
      (Polynomial.monomial quotientMono.deg
          (quotientCoefficient.val.toNat : ZMod p) *
        listSum p divisorTerms).coeff degree := by
  induction divisorTerms with
  | nil => simp
  | cons divisorTerm divisorTerms ih =>
      rcases divisorTerm with ⟨divisorMono, divisorCoefficient⟩
      simp only [List.foldr_cons, listSum_cons, mul_add,
        Polynomial.coeff_add, ih]
      by_cases hdegree : quotientMono.deg + divisorMono.deg = degree
      · simp only [hdegree, if_pos]
        subst degree
        simp [Zp.toZMod, Polynomial.monomial_mul_monomial]
      · simp only [hdegree, if_neg]
        rw [Polynomial.monomial_mul_monomial, Polynomial.coeff_monomial]
        simp [hdegree]

theorem pairVecDivVHCListProductCoeffValue_eq_coeff (p degree : Nat)
    (quotientTerms divisorTerms : List (UMonomial × Zp)) :
    pairVecDivVHCListProductCoeffValue p degree quotientTerms divisorTerms =
      (listSum p quotientTerms * listSum p divisorTerms).coeff degree := by
  induction quotientTerms with
  | nil => simp [pairVecDivVHCListProductCoeffValue]
  | cons quotientTerm quotientTerms ih =>
      rcases quotientTerm with ⟨quotientMono, quotientCoefficient⟩
      simp only [pairVecDivVHCListProductCoeffValue, listSum_cons, add_mul,
        Polynomial.coeff_add, ih]
      have hrow := pairVecDivVHCListProductCoeffValue_row p degree quotientMono
        quotientCoefficient divisorTerms
      simpa [Zp.toZMod] using congrArg
        (fun value => value +
          (listSum p quotientTerms * listSum p divisorTerms).coeff degree) hrow

theorem pairVecDivVHCListProductCoeffValue_toPoly (p degree : Nat)
    (quotient divisor : SparsePolyZp) :
    pairVecDivVHCListProductCoeffValue p degree quotient.toList divisor.toList =
      (SparsePolyZp.toPoly p quotient *
        SparsePolyZp.toPoly p divisor).coeff degree := by
  exact pairVecDivVHCListProductCoeffValue_eq_coeff p degree quotient.toList
    divisor.toList

theorem pairVecDivVHCListProductCoeffValue_divisorTail (p degree : Nat)
    (quotient divisor : SparsePolyZp) :
    pairVecDivVHCListProductCoeffValue p degree quotient.toList
        divisor.toList.tail =
      (SparsePolyZp.toPoly p quotient *
        listSum p divisor.toList.tail).coeff degree := by
  exact pairVecDivVHCListProductCoeffValue_eq_coeff p degree quotient.toList
    divisor.toList.tail

theorem pairVecDivVHCTargetPairSum_eq_productCoeffTail
    (p degree : Nat) (quotient divisor : SparsePolyZp) :
    (∑ pair ∈ PairVecDivVHCTargetPairsAtDegree degree quotient divisor,
        pairVecDivVHCIndexedPairProductValue p quotient divisor pair) =
      (SparsePolyZp.toPoly p quotient *
        listSum p divisor.toList.tail).coeff degree := by
  rw [pairVecDivVHCTargetPairSum_eq_listProductCoeffValue,
    pairVecDivVHCListProductCoeffValue_divisorTail]

/-- End-to-end coefficient bridge for the heap frontier: the exact sum of
products consumed from generated-C++ owner chains is the corresponding L2
coefficient of `quotient * divisor.tail`. -/
theorem pairVecDivVHCHeapOwnerSum_eq_productCoeffTail
    (p degreeLimit dividendIndex resetH : Nat)
    (dividend quotient divisor : SparsePolyZp)
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode)
    (owners : Nat → Finset Nat) (frontier : PairVecDivVHCFrontier)
    (hsize : nodes.size = divisor.size - 1)
    (hfixed : PairVecDivVHCNodeDivisorIndicesFixed nodes)
    (hstate : PairVecDivVHCStateCovered heap nodes #[] resetH)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hhomogeneous : PairVecDivVHCHeapChainsHomogeneous heap owners nodes)
    (hresetReady : PairVecDivVHCResetReady resetH quotient.size nodes)
    (hordered : PairVecDivVHCHeapOrdered heap nodes)
    (hdenotes : ∀ (i : Nat) (node : PairVecDivVHCNode),
      nodes[i]? = some node → node.mono ≠ none →
        PairVecDivVHCNodeDenotes quotient divisor node)
    (hcanonical : SparsePolyZp.Canonical p quotient)
    (hprefix : PairVecDivVHCCursorPrefixAbove degreeLimit nodes quotient divisor)
    (hfrontier : frontier.degree < degreeLimit)
    (hselect : pairVecDivVHCSelectFrontier dividendIndex dividend heap nodes =
      .ok frontier) :
    (∑ i ∈ PairVecDivVHCHeapOwnedNodesAtDegree frontier.degree heap owners
        nodes, pairVecDivVHCNodeProductValue p nodes quotient divisor i) =
      (SparsePolyZp.toPoly p quotient *
        listSum p divisor.toList.tail).coeff frontier.degree := by
  rw [pairVecDivVHCHeapOwnerSum_eq_targetPairSum p degreeLimit dividendIndex
    resetH dividend quotient divisor heap nodes owners frontier hsize hfixed
    hstate hownership hhomogeneous hresetReady hordered hdenotes hcanonical
    hprefix hfrontier hselect]
  exact pairVecDivVHCTargetPairSum_eq_productCoeffTail p frontier.degree
    quotient divisor

/-- No quotient/divisor-tail product can lie strictly between the selected
frontier and the previous outer-loop bound.  This is the semantic "no skipped
degree" fact behind the sparse heap traversal: reset rows point one past the
quotient and therefore place every existing quotient cell in the processed
cursor prefix, while every heap-owned row is at its cursor (at most the
selected frontier), before it (at least the old bound), or after it (strictly
below the cursor by canonical quotient order). -/
theorem pairVecDivVHCTargetPairsAtDegree_eq_empty_of_gap
    (p degreeLimit targetDegree dividendIndex resetH : Nat)
    (dividend quotient divisor : SparsePolyZp)
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode)
    (owners : Nat → Finset Nat) (frontier : PairVecDivVHCFrontier)
    (hsize : nodes.size = divisor.size - 1)
    (hfixed : PairVecDivVHCNodeDivisorIndicesFixed nodes)
    (hstate : PairVecDivVHCStateCovered heap nodes #[] resetH)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hhomogeneous : PairVecDivVHCHeapChainsHomogeneous heap owners nodes)
    (hresetReady : PairVecDivVHCResetReady resetH quotient.size nodes)
    (hordered : PairVecDivVHCHeapOrdered heap nodes)
    (hdenotes : ∀ (i : Nat) (node : PairVecDivVHCNode),
      nodes[i]? = some node → node.mono ≠ none →
        PairVecDivVHCNodeDenotes quotient divisor node)
    (hcanonical : SparsePolyZp.Canonical p quotient)
    (hprefix : PairVecDivVHCCursorPrefixAbove degreeLimit nodes quotient divisor)
    (hfrontier : frontier.degree < targetDegree)
    (htarget : targetDegree < degreeLimit)
    (hselect : pairVecDivVHCSelectFrontier dividendIndex dividend heap nodes =
      .ok frontier) :
    PairVecDivVHCTargetPairsAtDegree targetDegree quotient divisor = ∅ := by
  rw [← Finset.not_nonempty_iff_eq_empty]
  rintro ⟨pair, hpair⟩
  rcases pair with ⟨q, d⟩
  rw [PairVecDivVHCTargetPairsAtDegree] at hpair
  rcases Finset.mem_filter.mp hpair with ⟨hindices, hdegree⟩
  rcases Finset.mem_product.mp hindices with ⟨hq, hd⟩
  rw [Finset.mem_range] at hq
  rw [Finset.mem_Ico] at hd
  rcases hd with ⟨hdpos, hd⟩
  have hquotient : quotient[q]? = some quotient[q] :=
    Array.getElem?_eq_getElem hq
  have hdivisor : divisor[d]? = some divisor[d] :=
    Array.getElem?_eq_getElem hd
  simp [PairVecDivVHCPairAtDegree, hquotient, hdivisor] at hdegree
  rcases hfixed.node_for_tail nodes divisor.size d hsize hdpos hd with
    ⟨node, hnode, hnodeD⟩
  have hnodeIndex : d - 1 < nodes.size := by
    rw [hsize]
    omega
  have hcovered := hstate.covered_with heap nodes #[] resetH owners
    hownership
  rcases hcovered (d - 1) hnodeIndex with hreset | hlin | howned
  · rcases hresetReady.2 (d - 1) hreset with
      ⟨readyNode, hreadyNode, hcursor, hreadyD, hmono⟩
    rw [hnode] at hreadyNode
    simp only [Option.some.injEq] at hreadyNode
    subst readyNode
    have hearlier : q < node.quotientIndex := by omega
    have habove := hprefix (d - 1) node hnode q quotient[q] divisor[d]
      hearlier hquotient (by simpa [hnodeD] using hdivisor)
    omega
  · simp at hlin
  · rcases howned with ⟨slot, head, hheap, hmem⟩
    have howns := hownership.1 slot head hheap
    rcases pairVecDivVHCChainOwns_mem_active (some head) (owners head) nodes
        howns (d - 1) hmem with
      ⟨activeNode, mono, hactiveGet, hactiveMono⟩
    rw [hnode] at hactiveGet
    simp only [Option.some.injEq] at hactiveGet
    subst activeNode
    have hnodeDenotes := hdenotes (d - 1) node hnode (by
      simp [hactiveMono])
    rcases pairVecDivVHCOwnedNode_degree_le_frontier dividendIndex dividend
        heap nodes owners frontier hownership hhomogeneous hordered slot head
        (d - 1) node hheap hmem hnode hselect with
      ⟨storedMono, hstoredMono, hstoredLe⟩
    have hmonoEq : storedMono = mono := by
      rw [hactiveMono] at hstoredMono
      exact (Option.some.inj hstoredMono).symm
    subst storedMono
    rcases Nat.lt_trichotomy q node.quotientIndex with hearlier | heq | hlater
    · have habove := hprefix (d - 1) node hnode q quotient[q] divisor[d]
        hearlier hquotient (by simpa [hnodeD] using hdivisor)
      omega
    · subst q
      have hcursorDegree := hnodeDenotes.product_degree_eq quotient divisor node
        mono.deg (by simpa using hactiveMono) quotient[node.quotientIndex]
        divisor[d] hquotient (by simpa [hnodeD] using hdivisor)
      omega
    · have hbelow := hnodeDenotes.later_product_degree_lt p quotient divisor
        node mono.deg q quotient[q] divisor[d] hcanonical
        (by simpa using hactiveMono) hlater hquotient
        (by simpa [hnodeD] using hdivisor)
      omega

theorem pairVecDivVHCTail_product_coeff_eq_zero_of_gap
    (p degreeLimit targetDegree dividendIndex resetH : Nat)
    (dividend quotient divisor : SparsePolyZp)
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode)
    (owners : Nat → Finset Nat) (frontier : PairVecDivVHCFrontier)
    (hsize : nodes.size = divisor.size - 1)
    (hfixed : PairVecDivVHCNodeDivisorIndicesFixed nodes)
    (hstate : PairVecDivVHCStateCovered heap nodes #[] resetH)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hhomogeneous : PairVecDivVHCHeapChainsHomogeneous heap owners nodes)
    (hresetReady : PairVecDivVHCResetReady resetH quotient.size nodes)
    (hordered : PairVecDivVHCHeapOrdered heap nodes)
    (hdenotes : ∀ (i : Nat) (node : PairVecDivVHCNode),
      nodes[i]? = some node → node.mono ≠ none →
        PairVecDivVHCNodeDenotes quotient divisor node)
    (hcanonical : SparsePolyZp.Canonical p quotient)
    (hprefix : PairVecDivVHCCursorPrefixAbove degreeLimit nodes quotient divisor)
    (hfrontier : frontier.degree < targetDegree)
    (htarget : targetDegree < degreeLimit)
    (hselect : pairVecDivVHCSelectFrontier dividendIndex dividend heap nodes =
      .ok frontier) :
    (SparsePolyZp.toPoly p quotient *
      listSum p divisor.toList.tail).coeff targetDegree = 0 := by
  rw [← pairVecDivVHCTargetPairSum_eq_productCoeffTail p targetDegree
    quotient divisor,
    pairVecDivVHCTargetPairsAtDegree_eq_empty_of_gap p degreeLimit targetDegree
      dividendIndex resetH dividend quotient divisor heap nodes owners frontier
      hsize hfixed hstate hownership hhomogeneous hresetReady hordered hdenotes
      hcanonical hprefix hfrontier htarget hselect]
  simp

theorem pairVecDivVHCDividend_coeff_eq_zero_of_gap
    (p degreeLimit targetDegree dividendIndex : Nat)
    (dividend : SparsePolyZp) (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode) (frontier : PairVecDivVHCFrontier)
    (hcanonical : SparsePolyZp.Canonical p dividend)
    (hconsumed : PairVecDivVHCConsumedDividendAbove degreeLimit dividendIndex
      dividend)
    (hfrontier : frontier.degree < targetDegree)
    (htarget : targetDegree < degreeLimit)
    (hselect : pairVecDivVHCSelectFrontier dividendIndex dividend heap nodes =
      .ok frontier) :
    (SparsePolyZp.toPoly p dividend).coeff targetDegree = 0 := by
  unfold SparsePolyZp.toPoly
  suffices habsent : ∀ term ∈ dividend.toList,
      term.1.deg ≠ targetDegree by
    have go : ∀ terms : List (UMonomial × Zp),
        (∀ term ∈ terms, term.1.deg ≠ targetDegree) →
        (listSum p terms).coeff targetDegree = 0 := by
      intro terms hterms
      induction terms with
      | nil => simp [listSum]
      | cons term rest ih =>
          rw [listSum_cons, Polynomial.coeff_add,
            Polynomial.coeff_monomial,
            if_neg (hterms term List.mem_cons_self),
            ih (by
              intro item hitem
              exact hterms item (List.mem_cons_of_mem term hitem))]
          simp
    exact go dividend.toList habsent
  intro term hterm
  rcases List.getElem_of_mem hterm with ⟨i, hiList, htermEq⟩
  have hi : i < dividend.size := by simpa using hiList
  have htermEq' : term = dividend[i] := by
    rw [← Array.getElem_toList hi]
    exact htermEq.symm
  subst term
  change dividend[i].1.deg ≠ targetDegree
  have htermGet : dividend[i]? = some dividend[i] :=
    Array.getElem?_eq_getElem hi
  have hsource := pairVecDivVHCSelectFrontier_has_source dividendIndex
    dividend heap nodes frontier hselect
  cases hsource with
  | dividend hindex =>
      dsimp only at hfrontier
      by_cases hbefore : i < dividendIndex
      · have habove := hconsumed i dividend[i] hbefore htermGet
        omega
      · by_cases heq : i = dividendIndex
        · subst i
          omega
        · have hbelow := canonical_degree_lt_of_index_lt p dividend hcanonical
            dividendIndex i hindex hi (by omega)
          omega
  | heap hheap rootMono hmono hdominates =>
      dsimp only at hfrontier
      by_cases hbefore : i < dividendIndex
      · have habove := hconsumed i dividend[i] hbefore htermGet
        omega
      · have hindex : dividendIndex < dividend.size := by
          by_contra hnot
          omega
        have hheadBelow := hdominates hindex
        by_cases heq : i = dividendIndex
        · subst i
          omega
        · have hbelow := canonical_degree_lt_of_index_lt p dividend hcanonical
            dividendIndex i hindex hi (by omega)
          omega

theorem pairVecDivVHCConsumeChain_preserves_linReady
    (this : DenseUPolyZp) (current : Option Nat)
    (owner unvisited : Finset Nat) (k : UInt64)
    (nodes : Array PairVecDivVHCNode) (lin : Array Nat) (resetH : Nat)
    (quotient divisor : SparsePolyZp) (result : PairVecDivVHCBucketResult)
    (howns : PairVecDivVHCChainOwns current owner nodes)
    (hownerSubset : owner ⊆ unvisited)
    (hdisjoint : Disjoint lin.toList.toFinset owner)
    (hready : PairVecDivVHCLinReady lin nodes)
    (hrun : pairVecDivVHCConsumeChain this current unvisited k nodes lin
      resetH quotient divisor = .ok result) :
    PairVecDivVHCLinReady result.lin result.nodes ∧
      result.lin.toList.toFinset ⊆ lin.toList.toFinset ∪ owner := by
  cases current with
  | none =>
      rw [PairVecDivVHCChainOwns] at howns
      subst owner
      rw [pairVecDivVHCConsumeChain] at hrun
      simp only [Except.ok.injEq] at hrun
      subst result
      exact ⟨hready, by simp⟩
  | some nodeIndex =>
      rw [PairVecDivVHCChainOwns] at howns
      split at howns <;> try contradiction
      next hownerMem =>
        rcases howns with ⟨node, mono, hget, hmono, htailOwns⟩
        have hmem : nodeIndex ∈ unvisited := hownerSubset hownerMem
        have hn : nodeIndex < nodes.size := by
          by_contra hnot
          rw [Array.getElem?_eq_none (by omega)] at hget
          contradiction
        rw [Array.getElem?_eq_getElem hn] at hget
        simp only [Option.some.injEq] at hget
        subst node
        have houtside : nodeIndex ∉ lin.toList := by
          intro hlinMem
          exact Finset.disjoint_left.mp hdisjoint
            (by simpa only [List.mem_toFinset] using hlinMem) hownerMem
        rw [pairVecDivVHCConsumeChain] at hrun
        simp only [hmem, ↓reduceDIte] at hrun
        cases hconsume : pairVecDivVHCConsumeNode this nodeIndex k nodes lin
            resetH quotient divisor with
        | error fault => simp [hconsume] at hrun
        | ok step =>
            rcases step with ⟨k', nodes', lin', resetH', next⟩
            rw [hconsume] at hrun
            have hnext := pairVecDivVHCConsumeNode_next this nodeIndex k k'
              nodes nodes' lin lin' resetH resetH' next quotient divisor hn
              hconsume
            subst next
            have hlinStep := pairVecDivVHCConsumeNode_preserves_linReady this
              nodeIndex mono k k' nodes nodes' lin lin' resetH resetH'
              nodes[nodeIndex].next quotient divisor hn hmono houtside hready
              hconsume
            have htailOwns' := pairVecDivVHCChainOwns_congr_on
              nodes[nodeIndex].next (owner.erase nodeIndex) nodes nodes'
              htailOwns (by
                intro i hi
                exact pairVecDivVHCConsumeNode_get_ne this nodeIndex k k'
                  nodes nodes' lin lin' resetH resetH' nodes[nodeIndex].next
                  quotient divisor hconsume i
                  (Finset.mem_erase.mp hi).1.symm)
            have hownerErase : owner.erase nodeIndex ⊆
                unvisited.erase nodeIndex := by
              intro i hi
              exact Finset.mem_erase.mpr ⟨(Finset.mem_erase.mp hi).1,
                hownerSubset (Finset.mem_of_mem_erase hi)⟩
            have hdisjoint' : Disjoint lin'.toList.toFinset
                (owner.erase nodeIndex) := by
              exact Finset.disjoint_left.mpr (by
                intro i hiLin hiOwner
                have hiStep := hlinStep.2 hiLin
                rw [Finset.mem_insert] at hiStep
                rcases hiStep with rfl | hiOld
                · exact (Finset.mem_erase.mp hiOwner).1 rfl
                · exact Finset.disjoint_left.mp hdisjoint hiOld
                    (Finset.mem_of_mem_erase hiOwner))
            have ih := pairVecDivVHCConsumeChain_preserves_linReady this
              nodes[nodeIndex].next (owner.erase nodeIndex)
              (unvisited.erase nodeIndex) k' nodes' lin' resetH' quotient
              divisor result htailOwns' hownerErase hdisjoint' hlinStep.1 hrun
            refine ⟨ih.1, ?_⟩
            intro i hi
            have hiResult := ih.2 hi
            rw [Finset.mem_union] at hiResult ⊢
            rcases hiResult with hiLin | hiOwner
            · have hiStep := hlinStep.2 hiLin
              rw [Finset.mem_insert] at hiStep
              rcases hiStep with rfl | hiOld
              · exact Or.inr hownerMem
              · exact Or.inl hiOld
            · exact Or.inr (Finset.mem_of_mem_erase hiOwner)
termination_by owner.card
decreasing_by
  exact Finset.card_erase_lt_of_mem (by assumption)

theorem pairVecDivVHCConsumeRootBucket_preserves_linReady
    (this : DenseUPolyZp) (heap : Array Nat) (k : UInt64)
    (nodes : Array PairVecDivVHCNode) (lin : Array Nat) (resetH : Nat)
    (quotient divisor : SparsePolyZp) (result : PairVecDivVHCBucketResult)
    (rootOwner : Finset Nat) (hheap : 0 < heap.size)
    (hrootOwns : PairVecDivVHCChainOwns (some heap[0]) rootOwner nodes)
    (hdisjoint : Disjoint lin.toList.toFinset rootOwner)
    (hready : PairVecDivVHCLinReady lin nodes)
    (hrun : pairVecDivVHCConsumeRootBucket this heap k nodes lin resetH
      quotient divisor = .ok result) :
    PairVecDivVHCLinReady result.lin result.nodes ∧
      result.lin.toList.toFinset ⊆ lin.toList.toFinset ∪ rootOwner := by
  unfold pairVecDivVHCConsumeRootBucket at hrun
  simp only [hheap, ↓reduceDIte] at hrun
  exact pairVecDivVHCConsumeChain_preserves_linReady this (some heap[0])
    rootOwner (Finset.range nodes.size) k nodes lin resetH quotient divisor
    result hrootOwns
    (pairVecDivVHCChainOwns_subset_range _ _ _ hrootOwns) hdisjoint hready
    hrun

theorem pairVecDivVHCLast_mem_toFinset (lin : Array Nat)
    (hlin : 0 < lin.size) :
    lin[lin.size - 1] ∈ lin.toList.toFinset := by
  simp only [List.mem_toFinset]
  exact Array.getElem_mem_toList (by omega)

theorem pairVecDivVHCPop_toFinset_subset (lin : Array Nat) :
    lin.pop.toList.toFinset ⊆ lin.toList.toFinset := by
  intro nodeIndex hmem
  simp only [List.mem_toFinset] at hmem ⊢
  rw [Array.toList_pop] at hmem
  exact List.mem_of_mem_dropLast hmem

theorem pairVecDivVHCLast_not_mem_pop (lin : Array Nat)
    (hlin : 0 < lin.size) (hnodup : lin.toList.Nodup) :
    lin[lin.size - 1] ∉ lin.pop.toList.toFinset := by
  simp only [List.mem_toFinset, Array.toList_pop]
  intro hmem
  have hnonempty : lin.toList ≠ [] :=
    List.ne_nil_of_length_pos (by
      rw [Array.length_toList]
      exact hlin)
  have hdecomp := List.dropLast_append_getLast hnonempty
  have hnodupAppend :
      (lin.toList.dropLast ++ [lin.toList.getLast hnonempty]).Nodup := by
    rw [hdecomp]
    exact hnodup
  have hlastEq : lin.toList.getLast hnonempty = lin[lin.size - 1] := by
    rw [List.getLast_eq_getElem]
    simpa only [Array.length_toList] using
      (Array.getElem_toList (xs := lin) (i := lin.size - 1) (by omega))
  have hne := (List.nodup_append.mp hnodupAppend).2.2
    lin[lin.size - 1] hmem (lin.toList.getLast hnonempty) (by simp)
  exact hne hlastEq.symm

theorem pairVecDivVHCMem_pop_of_mem_ne_last (lin : Array Nat)
    (hlin : 0 < lin.size) (nodeIndex : Nat)
    (hmem : nodeIndex ∈ lin.toList.toFinset)
    (hne : nodeIndex ≠ lin[lin.size - 1]) :
    nodeIndex ∈ lin.pop.toList.toFinset := by
  simp only [List.mem_toFinset] at hmem ⊢
  rw [Array.toList_pop]
  have hnonempty : lin.toList ≠ [] :=
    List.ne_nil_of_length_pos (by
      rw [Array.length_toList]
      exact hlin)
  have hdecomp := List.dropLast_append_getLast hnonempty
  rw [← hdecomp] at hmem
  rcases List.mem_append.mp hmem with hprefix | hlast
  · exact hprefix
  · simp only [List.mem_singleton] at hlast
    have hlastEq : lin.toList.getLast hnonempty =
        lin[lin.size - 1] := by
      rw [List.getLast_eq_getElem]
      simpa only [Array.length_toList] using
        (Array.getElem_toList (xs := lin) (i := lin.size - 1) (by omega))
    exact False.elim (hne (hlast.trans hlastEq))

/-- Coverage effect of the equal-degree insertion branch.  The concrete heap
slot head is replaced by `newNode`, whose written `next` points to the old
head; consequently its exact owner is the inserted node together with the
old bucket owner.  The last `lin` member is thereby moved into that bucket. -/
theorem pairVecDivVHCSetMerge_preserves_stateCovered_pop
    (newNode slot : Nat) (heap : Array Nat)
    (owners : Nat → Finset Nat) (nodes nodes' : Array PairVecDivVHCNode)
    (lin : Array Nat) (resetH : Nat) (node : PairVecDivVHCNode)
    (mono : UMonomial) (hlin : 0 < lin.size)
    (hnew : newNode = lin[lin.size - 1]) (hslot : slot < heap.size)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hcovered : PairVecDivVHCNodesCovered heap owners lin resetH nodes)
    (hnode : nodes[newNode]? = some node) (hmono : node.mono = some mono)
    (hfreshHead : ∀ i : Nat, heap[i]? ≠ some newNode)
    (hfreshOwners : ∀ (i head : Nat), heap[i]? = some head →
      newNode ∉ owners head)
    (hset : pairVecDivVHCSetNext newNode (some heap[slot]) nodes =
      .ok nodes') :
    PairVecDivVHCStateCovered (heap.set slot newNode) nodes' lin.pop resetH := by
  let owners' := fun head => if head = newNode then
    insert newNode (owners heap[slot]) else owners head
  have hownership' := pairVecDivVHCHeapChainOwnership_merge_fresh heap owners
    slot newNode nodes nodes' node mono hslot hownership hnode hmono
    hfreshHead hfreshOwners hset
  refine ⟨owners', hownership', ?_⟩
  intro i hi
  have hsize := pairVecDivVHCSetNext_nodes_size newNode (some heap[slot])
    nodes nodes' hset
  have hiOld : i < nodes.size := by omega
  rcases hcovered i hiOld with hreset | hlinMem |
      ⟨oldSlot, head, hget, hmem⟩
  · exact Or.inl hreset
  · by_cases heq : i = newNode
    · subst i
      refine Or.inr (Or.inr ⟨slot, newNode, ?_, ?_⟩)
      · rw [Array.getElem?_set_self hslot]
      · simp [owners']
    · exact Or.inr (Or.inl
        (pairVecDivVHCMem_pop_of_mem_ne_last lin hlin i hlinMem (by
          rw [← hnew]
          exact heq)))
  · by_cases hhead : head = heap[slot]
    · refine Or.inr (Or.inr ⟨slot, newNode, ?_, ?_⟩)
      · rw [Array.getElem?_set_self hslot]
      · simp only [owners', if_pos rfl, Finset.mem_insert]
        exact Or.inr (by simpa [hhead] using hmem)
    · have hslotNe : slot ≠ oldSlot := by
        intro heq
        subst oldSlot
        rw [Array.getElem?_eq_getElem hslot] at hget
        exact hhead (Option.some.inj hget).symm
      have hheadNew : head ≠ newNode := by
        intro heq
        subst head
        exact hfreshHead oldSlot hget
      refine Or.inr (Or.inr ⟨oldSlot, head, ?_, ?_⟩)
      · rw [Array.getElem?_set_ne hslot hslotNe]
        exact hget
      · simp only [owners', hheadNew, ↓reduceIte]
        exact hmem

/-- Coverage effect shared by the fresh-head insertion branches (empty heap,
root bubble, and bubble-below).  The heap rewrite must provide the concrete
reverse-survival witness for every head of `heap.push newNode`. -/
theorem pairVecDivVHCFreshHead_preserves_stateCovered_pop
    (newNode : Nat) (heap heap' : Array Nat)
    (owners : Nat → Finset Nat) (nodes nodes' : Array PairVecDivVHCNode)
    (lin : Array Nat) (resetH : Nat) (hlin : 0 < lin.size)
    (hnew : newNode = lin[lin.size - 1])
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hcovered : PairVecDivVHCNodesCovered heap owners lin resetH nodes)
    (hfreshHead : ∀ i : Nat, heap[i]? ≠ some newNode)
    (hownership' : PairVecDivVHCHeapChainOwnership heap'
      (fun head => if head = newNode then {newNode} else owners head) nodes')
    (hsize : nodes'.size = nodes.size)
    (hsurvives : ∀ (slot head : Nat),
      (heap.push newNode)[slot]? = some head →
        ∃ targetSlot : Nat, heap'[targetSlot]? = some head) :
    PairVecDivVHCStateCovered heap' nodes' lin.pop resetH := by
  let owners' := fun head => if head = newNode then {newNode} else owners head
  refine ⟨owners', hownership', ?_⟩
  intro i hi
  have hiOld : i < nodes.size := by omega
  rcases hcovered i hiOld with hreset | hlinMem |
      ⟨oldSlot, head, hget, hmem⟩
  · exact Or.inl hreset
  · by_cases heq : i = newNode
    · subst i
      have hpushGet : (heap.push newNode)[heap.size]? = some newNode := by simp
      rcases hsurvives heap.size newNode hpushGet with
        ⟨targetSlot, htargetGet⟩
      refine Or.inr (Or.inr ⟨targetSlot, newNode, htargetGet, ?_⟩)
      simp [owners']
    · exact Or.inr (Or.inl
        (pairVecDivVHCMem_pop_of_mem_ne_last lin hlin i hlinMem (by
          rw [← hnew]
          exact heq)))
  · have holdBound : oldSlot < heap.size := by
        by_contra hnot
        rw [Array.getElem?_eq_none (by omega)] at hget
        contradiction
    have hpushGet : (heap.push newNode)[oldSlot]? = some head := by
      rw [Array.getElem?_push]
      simp only [Nat.ne_of_lt holdBound, ↓reduceIte]
      exact hget
    rcases hsurvives oldSlot head hpushGet with
      ⟨targetSlot, htargetGet⟩
    have hheadNew : head ≠ newNode := by
      intro heq
      subst head
      exact hfreshHead oldSlot hget
    refine Or.inr (Or.inr ⟨targetSlot, head, htargetGet, ?_⟩)
    simp only [owners', hheadNew, ↓reduceIte]
    exact hmem

theorem pairVecDivVHCInsert_preserves_stateCovered_pop
    (heap heap' : Array Nat) (nodes nodes' : Array PairVecDivVHCNode)
    (lin : Array Nat) (resetH : Nat) (hlin : 0 < lin.size)
    (haway : PairVecDivVHCHeapChainsOwnedAway heap nodes
      lin.toList.toFinset)
    (hready : PairVecDivVHCLinReady lin nodes)
    (hcovered : PairVecDivVHCStateCovered heap nodes lin resetH)
    (hrun : pairVecDivVHCInsert lin[lin.size - 1] heap nodes =
      .ok (heap', nodes')) :
    PairVecDivVHCStateCovered heap' nodes' lin.pop resetH := by
  let newNode := lin[lin.size - 1]
  change pairVecDivVHCInsert newNode heap nodes = .ok (heap', nodes') at hrun
  unfold pairVecDivVHCInsert at hrun
  have hlastMem := pairVecDivVHCLast_mem_toFinset lin hlin
  rcases hready.2 newNode (by
      simpa only [newNode, List.mem_toFinset] using hlastMem) with
    ⟨node, mono, hnode, hmono⟩
  rcases haway with ⟨owners, hownership, hseparated⟩
  have hnodesCovered := hcovered.covered_with heap nodes lin resetH owners
    hownership
  have hfreshHead : ∀ slot : Nat, heap[slot]? ≠ some newNode := by
    intro slot hget
    exact (hseparated slot newNode hget).2 hlastMem
  have hfreshOwners : ∀ (slot head : Nat), heap[slot]? = some head →
      newNode ∉ owners head := by
    intro slot head hget hmem
    exact Finset.disjoint_left.mp (hseparated slot head hget).1 hlastMem hmem
  cases hnewMono : pairVecDivVHCMono newNode nodes with
  | error fault => simp [hnewMono] at hrun
  | ok newMono =>
      rw [hnewMono] at hrun
      by_cases hempty : heap.size = 0
      · have heq : heap = #[] := Array.eq_empty_of_size_eq_zero hempty
        subst heap
        simp only [Array.size_empty, ↓reduceDIte] at hrun
        cases hset : pairVecDivVHCSetNext newNode none nodes with
        | error fault => simp [hset] at hrun
        | ok updated =>
            rw [hset] at hrun
            simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
            rcases hrun with ⟨rfl, rfl⟩
            have hownership' :=
              pairVecDivVHCHeapChainOwnership_push_fresh #[] owners newNode
                nodes updated node mono hownership hnode hmono (by simp)
                (by simp) hset
            exact pairVecDivVHCFreshHead_preserves_stateCovered_pop newNode
              #[] #[newNode] owners nodes updated lin resetH hlin rfl
              hownership hnodesCovered (by simp) hownership'
              (pairVecDivVHCSetNext_nodes_size newNode none nodes updated hset)
              (by intro slot head hget; refine ⟨slot, ?_⟩; simpa using hget)
      · have hheap : 0 < heap.size := Nat.pos_of_ne_zero hempty
        cases hroot : pairVecDivVHCMono heap[0] nodes with
        | error fault =>
            simp [hempty, hroot] at hrun
        | ok rootMono =>
            by_cases hequal : newMono.deg = rootMono.deg
            · simp only [hempty, ↓reduceDIte, hroot, hequal] at hrun
              cases hset : pairVecDivVHCSetNext newNode (some heap[0]) nodes with
              | error fault => simp [hset] at hrun
              | ok updated =>
                  rw [hset] at hrun
                  simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
                  rcases hrun with ⟨rfl, rfl⟩
                  exact pairVecDivVHCSetMerge_preserves_stateCovered_pop
                    newNode 0 heap owners nodes updated lin resetH node mono
                    hlin rfl hheap hownership hnodesCovered hnode hmono
                    hfreshHead hfreshOwners hset
            · by_cases hgreater : newMono.deg > rootMono.deg
              · simp only [hempty, ↓reduceDIte, hroot, hequal,
                  hgreater] at hrun
                cases hset : pairVecDivVHCSetNext newNode none nodes with
                | error fault => simp [hset] at hrun
                | ok updated =>
                    rw [hset] at hrun
                    cases hbubble : pairVecDivVHCBubble heap.size 0 newNode
                        (heap.push newNode) with
                    | error fault => simp [hbubble] at hrun
                    | ok shifted =>
                        rw [hbubble] at hrun
                        simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
                        rcases hrun with ⟨rfl, rfl⟩
                        have hownership' :=
                          pairVecDivVHCHeapChainOwnership_bubble_fresh heap
                            shifted owners 0 newNode nodes updated node mono
                            hownership hnode hmono hfreshHead hfreshOwners hset
                            hbubble
                        exact pairVecDivVHCFreshHead_preserves_stateCovered_pop
                          newNode heap shifted owners nodes updated lin resetH
                          hlin rfl hownership hnodesCovered hfreshHead
                          hownership'
                          (pairVecDivVHCSetNext_nodes_size newNode none nodes
                            updated hset)
                          (pairVecDivVHCBubble_push_preserves_every_head 0
                            newNode heap shifted hownership.2.1 hfreshHead
                            hbubble)
              · simp only [hempty, ↓reduceDIte, hroot, hequal,
                  hgreater] at hrun
                cases hanchor : pairVecDivVHCFindAnchor newMono.deg
                    (pairVecDivVHCParent heap.size) heap nodes with
                | error fault =>
                    rw [hanchor] at hrun
                    contradiction
                | ok anchor =>
                    rw [hanchor] at hrun
                    by_cases ha : anchor < heap.size
                    · simp only [ha, ↓reduceDIte] at hrun
                      cases hanchorMono : pairVecDivVHCMono heap[anchor]
                          nodes with
                      | error fault =>
                          rw [hanchorMono] at hrun
                          contradiction
                      | ok anchorMono =>
                          rw [hanchorMono] at hrun
                          by_cases hequalAnchor :
                              newMono.deg = anchorMono.deg
                          · simp only [hequalAnchor, ↓reduceDIte] at hrun
                            cases hset : pairVecDivVHCSetNext newNode
                                (some heap[anchor]) nodes with
                            | error fault => simp [hset] at hrun
                            | ok updated =>
                                rw [hset] at hrun
                                simp only [Except.ok.injEq,
                                  Prod.mk.injEq] at hrun
                                rcases hrun with ⟨rfl, rfl⟩
                                exact
                                  pairVecDivVHCSetMerge_preserves_stateCovered_pop
                                    newNode anchor heap owners nodes updated lin
                                    resetH node mono hlin rfl ha hownership
                                    hnodesCovered hnode hmono hfreshHead
                                    hfreshOwners hset
                          · simp only [hequalAnchor, ↓reduceDIte] at hrun
                            cases hset : pairVecDivVHCSetNext newNode none nodes with
                            | error fault => simp [hset] at hrun
                            | ok updated =>
                                rw [hset] at hrun
                                cases hbubble : pairVecDivVHCBubbleBelow
                                    heap.size anchor newNode
                                    (heap.push newNode) with
                                | error fault => simp [hbubble] at hrun
                                | ok shifted =>
                                    rw [hbubble] at hrun
                                    simp only [Except.ok.injEq,
                                      Prod.mk.injEq] at hrun
                                    rcases hrun with ⟨rfl, rfl⟩
                                    have hownership' :=
                                      pairVecDivVHCHeapChainOwnership_bubbleBelow_fresh
                                        heap shifted owners anchor newNode nodes
                                        updated node mono hownership hnode hmono
                                        hfreshHead hfreshOwners hset hbubble
                                    exact
                                      pairVecDivVHCFreshHead_preserves_stateCovered_pop
                                        newNode heap shifted owners nodes updated
                                        lin resetH hlin rfl hownership
                                        hnodesCovered hfreshHead hownership'
                                        (pairVecDivVHCSetNext_nodes_size newNode
                                          none nodes updated hset)
                                        (pairVecDivVHCBubbleBelow_push_preserves_every_head
                                          anchor newNode heap shifted
                                          hownership.2.1 hfreshHead hbubble)
                    · simp only [ha, ↓reduceDIte] at hrun
                      contradiction

/-- One concrete reset activation followed by its concrete heap insertion
moves exactly `resetH - 1` out of the inactive prefix and into the heap while
preserving the pre-existing `lin` stack.  The proof reuses the complete
`VHC_insert` branch theorem by representing the intermediate location as
`lin.push nodeIndex`. -/
theorem pairVecDivVHCActivateInsert_preserves_stateCovered
    (resetH quotientSize : Nat) (heap heap' : Array Nat)
    (nodes activated inserted : Array PairVecDivVHCNode)
    (lin : Array Nat) (quotient divisor : SparsePolyZp)
    (hreset : 0 < resetH)
    (haway : PairVecDivVHCHeapChainsOwnedAway heap nodes
      lin.toList.toFinset)
    (hlinReady : PairVecDivVHCLinReady lin nodes)
    (hresetReady : PairVecDivVHCResetReady resetH quotientSize nodes)
    (hcovered : PairVecDivVHCStateCovered heap nodes lin resetH)
    (hactivate : pairVecDivVHCActivate (resetH - 1) nodes quotient divisor =
      .ok activated)
    (hinsert : pairVecDivVHCInsert (resetH - 1) heap activated =
      .ok (heap', inserted)) :
    PairVecDivVHCStateCovered heap' inserted lin (resetH - 1) := by
  let nodeIndex := resetH - 1
  have hindex : nodeIndex < resetH := by omega
  rcases hresetReady.2 nodeIndex hindex with
    ⟨oldNode, holdGet, hqIndex, hdIndex, hinactive⟩
  have hnotLin : nodeIndex ∉ lin.toList := by
    intro hmem
    rcases hlinReady.2 nodeIndex hmem with
      ⟨activeNode, activeMono, hactiveGet, hactiveMono⟩
    rw [holdGet] at hactiveGet
    simp only [Option.some.injEq] at hactiveGet
    subst activeNode
    rw [hinactive] at hactiveMono
    contradiction
  rcases haway with ⟨owners, hownership, hseparated⟩
  have hfresh := pairVecDivVHCHeapChainOwnership_fresh_of_mono_none heap
    owners nodes nodeIndex oldNode hownership holdGet hinactive
  have hownershipActivated :
      PairVecDivVHCHeapChainOwnership heap owners activated := by
    refine ⟨?_, hownership.2.1, hownership.2.2⟩
    intro slot head hget
    exact pairVecDivVHCChainOwns_congr_on (some head) (owners head) nodes
      activated (hownership.1 slot head hget) (by
        intro i hi
        exact pairVecDivVHCActivate_get_ne nodeIndex nodes activated quotient
          divisor hactivate i (by
            intro heq
            subst i
            exact hfresh.2 slot head hget hi))
  rcases pairVecDivVHCActivate_get nodeIndex nodes activated quotient divisor
      hactivate with ⟨hn, hq, hd, hnewGet⟩
  let activeMono : UMonomial :=
    ⟨quotient[nodes[nodeIndex].quotientIndex].1.deg +
      divisor[nodes[nodeIndex].divisorIndex].1.deg⟩
  let newNode : PairVecDivVHCNode := { nodes[nodeIndex] with
    mono := some activeMono
    next := none }
  have hnewGet' : activated[nodeIndex]? = some newNode := by
    simpa [newNode, activeMono] using hnewGet
  have hnewMono : newNode.mono = some activeMono := by rfl
  have hlinReadyPush :
      PairVecDivVHCLinReady (lin.push nodeIndex) activated := by
    refine ⟨?_, ?_⟩
    · simp only [Array.toList_push, List.nodup_append, hlinReady.1,
        List.nodup_singleton, true_and, List.mem_singleton]
      intro i hi j hj
      subst j
      intro heq
      subst i
      exact hnotLin hi
    · intro i hi
      simp only [Array.toList_push, List.mem_append,
        List.mem_singleton] at hi
      rcases hi with hi | rfl
      · rcases hlinReady.2 i hi with
          ⟨oldActive, oldMono, hget, hmono⟩
        refine ⟨oldActive, oldMono, ?_, hmono⟩
        rw [pairVecDivVHCActivate_get_ne nodeIndex nodes activated quotient
          divisor hactivate i (by
            intro heq
            subst i
            exact hnotLin hi)]
        exact hget
      · exact ⟨newNode, activeMono, hnewGet', hnewMono⟩
  have hprotected : (lin.push nodeIndex).toList.toFinset =
      insert nodeIndex lin.toList.toFinset := by
    ext i
    simp [Array.toList_push, or_comm]
  have hawayPush : PairVecDivVHCHeapChainsOwnedAway heap activated
      (lin.push nodeIndex).toList.toFinset := by
    refine ⟨owners, hownershipActivated, ?_⟩
    intro slot head hget
    have hold := hseparated slot head hget
    rw [hprotected]
    constructor
    · exact Finset.disjoint_left.mpr (by
        intro i hi howner
        simp only [Finset.mem_insert] at hi
        rcases hi with rfl | hi
        · exact hfresh.2 slot head hget howner
        · exact Finset.disjoint_left.mp hold.1 hi howner)
    · intro hi
      simp only [Finset.mem_insert] at hi
      rcases hi with heq | hi
      · subst head
        exact hfresh.1 slot hget
      · exact hold.2 hi
  have hnodesCovered := hcovered.covered_with heap nodes lin resetH owners
    hownership
  have hstatePush : PairVecDivVHCStateCovered heap activated
      (lin.push nodeIndex) nodeIndex := by
    refine ⟨owners, hownershipActivated, ?_⟩
    intro i hi
    have hsize := pairVecDivVHCActivate_nodes_size nodeIndex nodes activated
      quotient divisor hactivate
    have hiOld : i < nodes.size := by omega
    rcases hnodesCovered i hiOld with hprefix | hlin | hheap
    · by_cases heq : i = nodeIndex
      · subst i
        exact Or.inr (Or.inl (by simp [Array.toList_push]))
      · exact Or.inl (by omega)
    · exact Or.inr (Or.inl (by
        simp only [List.mem_toFinset, Array.toList_push,
          List.mem_append, List.mem_singleton]
        exact Or.inl (by simpa only [List.mem_toFinset] using hlin)))
    · exact Or.inr (Or.inr hheap)
  have hstep := pairVecDivVHCInsert_preserves_stateCovered_pop heap heap'
    activated inserted (lin.push nodeIndex) nodeIndex (by simp) hawayPush
    hlinReadyPush hstatePush (by simpa [nodeIndex] using hinsert)
  simpa [Array.pop_push] using hstep

theorem PairVecDivVHCLinReady.pop_after_insert
    (lin : Array Nat) (heap heap' : Array Nat)
    (nodes nodes' : Array PairVecDivVHCNode)
    (hlin : 0 < lin.size) (hready : PairVecDivVHCLinReady lin nodes)
    (hinsert : pairVecDivVHCInsert lin[lin.size - 1] heap nodes =
      .ok (heap', nodes')) :
    PairVecDivVHCLinReady lin.pop nodes' := by
  refine ⟨?_, ?_⟩
  · rw [Array.toList_pop]
    have hnonempty : lin.toList ≠ [] :=
      List.ne_nil_of_length_pos (by
        rw [Array.length_toList]
        exact hlin)
    have happ :
        (lin.toList.dropLast ++ [lin.toList.getLast hnonempty]).Nodup := by
      rw [List.dropLast_append_getLast hnonempty]
      exact hready.1
    exact (List.nodup_append.mp happ).1
  · intro nodeIndex hmem
    have hmemOld : nodeIndex ∈ lin.toList := by
      rw [Array.toList_pop] at hmem
      exact List.mem_of_mem_dropLast hmem
    rcases hready.2 nodeIndex hmemOld with ⟨node, mono, hnode, hmono⟩
    have hne : lin[lin.size - 1] ≠ nodeIndex := by
      intro heq
      subst nodeIndex
      exact pairVecDivVHCLast_not_mem_pop lin hlin hready.1 (by
        simpa only [List.mem_toFinset] using hmem)
    refine ⟨node, mono, ?_, hmono⟩
    rw [pairVecDivVHCInsert_get_ne lin[lin.size - 1] heap heap' nodes nodes'
      hinsert nodeIndex hne]
    exact hnode

/-- The literal C++ reverse-`lin` loop restores exact ownership of every heap
chain.  Termination is Lean's well-founded recursion on `lin.size`; the
statement and proof use only the concrete generated execution state. -/
theorem pairVecDivVHCReinsertLin_preserves_heapChainsOwned
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode) (lin : Array Nat)
    (state : PairVecDivVHCHeapState)
    (haway : PairVecDivVHCHeapChainsOwnedAway heap nodes
      lin.toList.toFinset)
    (hready : PairVecDivVHCLinReady lin nodes)
    (hrun : pairVecDivVHCReinsertLin heap nodes lin = .ok state) :
    PairVecDivVHCHeapChainsOwned state.heap state.nodes := by
  rw [pairVecDivVHCReinsertLin] at hrun
  split at hrun
  next hlin =>
    dsimp only at hrun
    cases hinsert : pairVecDivVHCInsert lin[lin.size - 1] heap nodes with
    | error fault => simp [hinsert] at hrun
    | ok inserted =>
        rcases inserted with ⟨heap', nodes'⟩
        rw [hinsert] at hrun
        have hlastMem := pairVecDivVHCLast_mem_toFinset lin hlin
        rcases hready.2 lin[lin.size - 1]
            (by simpa only [List.mem_toFinset] using hlastMem) with
          ⟨node, mono, hnode, hmono⟩
        rcases haway with ⟨owners, hownership, hseparated⟩
        have haway' := pairVecDivVHCInsert_preserves_away_of_protected
          lin[lin.size - 1] heap heap' nodes nodes' lin.toList.toFinset
          lin.pop.toList.toFinset owners node mono hownership hseparated
          hlastMem (pairVecDivVHCPop_toFinset_subset lin)
          (pairVecDivVHCLast_not_mem_pop lin hlin hready.1) hnode hmono hinsert
        have hready' := hready.pop_after_insert lin heap heap' nodes nodes'
          hlin hinsert
        exact pairVecDivVHCReinsertLin_preserves_heapChainsOwned heap' nodes'
          lin.pop state haway' hready' hrun
  next hlin =>
    simp only [Except.ok.injEq] at hrun
    subst state
    rcases haway with ⟨owners, hownership, hseparated⟩
    exact ⟨owners, hownership⟩
termination_by lin.size
decreasing_by
  simp only [Array.size_pop]
  omega

/-- The literal reverse-`lin` reinsertion loop preserves max-heap order at
every generated insertion step. -/
theorem pairVecDivVHCReinsertLin_preserves_heapOrdered
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode) (lin : Array Nat)
    (state : PairVecDivVHCHeapState)
    (haway : PairVecDivVHCHeapChainsOwnedAway heap nodes
      lin.toList.toFinset)
    (hready : PairVecDivVHCLinReady lin nodes)
    (hordered : PairVecDivVHCHeapOrdered heap nodes)
    (hrun : pairVecDivVHCReinsertLin heap nodes lin = .ok state) :
    PairVecDivVHCHeapOrdered state.heap state.nodes := by
  rw [pairVecDivVHCReinsertLin] at hrun
  split at hrun
  next hlin =>
    dsimp only at hrun
    cases hinsert : pairVecDivVHCInsert lin[lin.size - 1] heap nodes with
    | error fault => simp [hinsert] at hrun
    | ok inserted =>
        rcases inserted with ⟨heap', nodes'⟩
        rw [hinsert] at hrun
        have hlastMem := pairVecDivVHCLast_mem_toFinset lin hlin
        rcases hready.2 lin[lin.size - 1]
            (by simpa only [List.mem_toFinset] using hlastMem) with
          ⟨node, mono, hnode, hmono⟩
        rcases haway with ⟨owners, hownership, hseparated⟩
        have hordered' := pairVecDivVHCInsert_preserves_heapOrdered
          lin[lin.size - 1] heap heap' nodes nodes'
          (hownership.heapPointersValid heap owners nodes) hordered hinsert
        have haway' := pairVecDivVHCInsert_preserves_away_of_protected
          lin[lin.size - 1] heap heap' nodes nodes' lin.toList.toFinset
          lin.pop.toList.toFinset owners node mono hownership hseparated
          hlastMem (pairVecDivVHCPop_toFinset_subset lin)
          (pairVecDivVHCLast_not_mem_pop lin hlin hready.1) hnode hmono hinsert
        have hready' := hready.pop_after_insert lin heap heap' nodes nodes'
          hlin hinsert
        exact pairVecDivVHCReinsertLin_preserves_heapOrdered heap' nodes'
          lin.pop state haway' hready' hordered' hrun
  next hlin =>
    simp only [Except.ok.injEq] at hrun
    subst state
    exact hordered
termination_by lin.size
decreasing_by
  simp only [Array.size_pop]
  omega

theorem pairVecDivVHCReinsertLin_preserves_heapChainsHomogeneous
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode) (lin : Array Nat)
    (state : PairVecDivVHCHeapState) (owners : Nat → Finset Nat)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hseparated : ∀ (slot head : Nat), heap[slot]? = some head →
      Disjoint lin.toList.toFinset (owners head) ∧
        head ∉ lin.toList.toFinset)
    (hhomogeneous : PairVecDivVHCHeapChainsHomogeneous heap owners nodes)
    (hready : PairVecDivVHCLinReady lin nodes)
    (hrun : pairVecDivVHCReinsertLin heap nodes lin = .ok state) :
    ∃ owners', PairVecDivVHCHeapChainOwnership state.heap owners' state.nodes ∧
      PairVecDivVHCHeapChainsHomogeneous state.heap owners' state.nodes := by
  rw [pairVecDivVHCReinsertLin] at hrun
  split at hrun
  next hlin =>
    dsimp only at hrun
    cases hinsert : pairVecDivVHCInsert lin[lin.size - 1] heap nodes with
    | error fault => simp [hinsert] at hrun
    | ok inserted =>
        rcases inserted with ⟨heap', nodes'⟩
        rw [hinsert] at hrun
        have hlastMem := pairVecDivVHCLast_mem_toFinset lin hlin
        rcases hready.2 lin[lin.size - 1]
            (by simpa only [List.mem_toFinset] using hlastMem) with
          ⟨node, mono, hnode, hmono⟩
        have hfreshHead : ∀ slot : Nat,
            heap[slot]? ≠ some lin[lin.size - 1] := by
          intro slot hget
          exact (hseparated slot lin[lin.size - 1] hget).2 hlastMem
        have hfreshOwners : ∀ (slot head : Nat), heap[slot]? = some head →
            lin[lin.size - 1] ∉ owners head := by
          intro slot head hget hmem
          exact Finset.disjoint_left.mp (hseparated slot head hget).1
            hlastMem hmem
        rcases pairVecDivVHCInsert_preserves_heapChainsHomogeneous_of_fresh
            lin[lin.size - 1] heap heap' nodes nodes' owners node mono
            hownership hhomogeneous hnode hmono hfreshHead hfreshOwners hinsert
          with ⟨homOwners, hhomOwnership, hhomogeneous'⟩
        have haway' := pairVecDivVHCInsert_preserves_away_of_protected
          lin[lin.size - 1] heap heap' nodes nodes' lin.toList.toFinset
          lin.pop.toList.toFinset owners node mono hownership hseparated
          hlastMem (pairVecDivVHCPop_toFinset_subset lin)
          (pairVecDivVHCLast_not_mem_pop lin hlin hready.1) hnode hmono hinsert
        rcases haway' with ⟨nextOwners, hnextOwnership, hnextSeparated⟩
        have hnextHomogeneous := hhomogeneous'.congr_owners heap' homOwners
          nextOwners nodes' hhomOwnership hnextOwnership
        have hready' := hready.pop_after_insert lin heap heap' nodes nodes'
          hlin hinsert
        exact pairVecDivVHCReinsertLin_preserves_heapChainsHomogeneous heap'
          nodes' lin.pop state nextOwners hnextOwnership hnextSeparated
          hnextHomogeneous hready' hrun
  next hlin =>
    simp only [Except.ok.injEq] at hrun
    subst state
    exact ⟨owners, hownership, hhomogeneous⟩
termination_by lin.size
decreasing_by
  simp only [Array.size_pop]
  omega

theorem pairVecDivVHCReinsertLin_preserves_cursorPrefixAbove
    (degreeLimit : Nat) (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode) (lin : Array Nat)
    (quotient divisor : SparsePolyZp) (state : PairVecDivVHCHeapState)
    (hprefix : PairVecDivVHCCursorPrefixAbove degreeLimit nodes quotient divisor)
    (hrun : pairVecDivVHCReinsertLin heap nodes lin = .ok state) :
    PairVecDivVHCCursorPrefixAbove degreeLimit state.nodes quotient divisor := by
  rw [pairVecDivVHCReinsertLin] at hrun
  split at hrun
  next hlin =>
    dsimp only at hrun
    cases hinsert : pairVecDivVHCInsert lin[lin.size - 1] heap nodes with
    | error fault => simp [hinsert] at hrun
    | ok inserted =>
        rcases inserted with ⟨heap', nodes'⟩
        rw [hinsert] at hrun
        have hprefix' := pairVecDivVHCInsert_preserves_cursorPrefixAbove
          degreeLimit lin[lin.size - 1] heap heap' nodes nodes' quotient
          divisor hprefix hinsert
        exact pairVecDivVHCReinsertLin_preserves_cursorPrefixAbove degreeLimit
          heap' nodes' lin.pop quotient divisor state hprefix' hrun
  next hlin =>
    simp only [Except.ok.injEq] at hrun
    subst state
    exact hprefix
termination_by lin.size
decreasing_by
  simp only [Array.size_pop]
  omega

/-- The literal reverse-`lin` loop preserves total node coverage together with
the exact heap ownership witness it constructs.  Recursion is on the concrete
shrinking `lin.size`; each step uses the full branch analysis of the executed
`VHC_insert`. -/
theorem pairVecDivVHCReinsertLin_preserves_stateCovered
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode) (lin : Array Nat)
    (resetH : Nat) (state : PairVecDivVHCHeapState)
    (haway : PairVecDivVHCHeapChainsOwnedAway heap nodes
      lin.toList.toFinset)
    (hready : PairVecDivVHCLinReady lin nodes)
    (hcovered : PairVecDivVHCStateCovered heap nodes lin resetH)
    (hrun : pairVecDivVHCReinsertLin heap nodes lin = .ok state) :
    PairVecDivVHCStateCovered state.heap state.nodes #[] resetH := by
  rw [pairVecDivVHCReinsertLin] at hrun
  split at hrun
  next hlin =>
    dsimp only at hrun
    cases hinsert : pairVecDivVHCInsert lin[lin.size - 1] heap nodes with
    | error fault => simp [hinsert] at hrun
    | ok inserted =>
        rcases inserted with ⟨heap', nodes'⟩
        rw [hinsert] at hrun
        have hlastMem := pairVecDivVHCLast_mem_toFinset lin hlin
        rcases hready.2 lin[lin.size - 1]
            (by simpa only [List.mem_toFinset] using hlastMem) with
          ⟨node, mono, hnode, hmono⟩
        rcases haway with ⟨owners, hownership, hseparated⟩
        have haway' := pairVecDivVHCInsert_preserves_away_of_protected
          lin[lin.size - 1] heap heap' nodes nodes' lin.toList.toFinset
          lin.pop.toList.toFinset owners node mono hownership hseparated
          hlastMem (pairVecDivVHCPop_toFinset_subset lin)
          (pairVecDivVHCLast_not_mem_pop lin hlin hready.1) hnode hmono hinsert
        have hready' := hready.pop_after_insert lin heap heap' nodes nodes'
          hlin hinsert
        have hcovered' := pairVecDivVHCInsert_preserves_stateCovered_pop
          heap heap' nodes nodes' lin resetH hlin
          ⟨owners, hownership, hseparated⟩ hready hcovered hinsert
        exact pairVecDivVHCReinsertLin_preserves_stateCovered heap' nodes'
          lin.pop resetH state haway' hready' hcovered' hrun
  next hlin =>
    simp only [Except.ok.injEq] at hrun
    subst state
    have hempty : lin = #[] := Array.eq_empty_of_size_eq_zero (by omega)
    subst lin
    exact hcovered
termination_by lin.size
decreasing_by
  simp only [Array.size_pop]
  omega

theorem pairVecDivVHCReinsertLin_preserves_allActiveNodesBelow
    (degreeLimit : Nat) (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode) (lin : Array Nat)
    (state : PairVecDivVHCHeapState)
    (hbelow : PairVecDivVHCAllActiveNodesBelow degreeLimit nodes)
    (hrun : pairVecDivVHCReinsertLin heap nodes lin = .ok state) :
    PairVecDivVHCAllActiveNodesBelow degreeLimit state.nodes := by
  rw [pairVecDivVHCReinsertLin] at hrun
  split at hrun
  next hlin =>
    dsimp only at hrun
    cases hinsert : pairVecDivVHCInsert lin[lin.size - 1] heap nodes with
    | error fault => simp [hinsert] at hrun
    | ok inserted =>
        rcases inserted with ⟨heap', nodes'⟩
        rw [hinsert] at hrun
        have hbelow' := pairVecDivVHCInsert_preserves_allActiveNodesBelow
          degreeLimit lin[lin.size - 1] heap heap' nodes nodes' hbelow hinsert
        exact pairVecDivVHCReinsertLin_preserves_allActiveNodesBelow
          degreeLimit heap' nodes' lin.pop state hbelow' hrun
  next hlin =>
    simp only [Except.ok.injEq] at hrun
    subst state
    exact hbelow
termination_by lin.size
decreasing_by
  simp only [Array.size_pop]
  omega

theorem pairVecDivVHCReinsertLin_preserves_divisorIndicesFixed
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode) (lin : Array Nat)
    (state : PairVecDivVHCHeapState)
    (hfixed : PairVecDivVHCNodeDivisorIndicesFixed nodes)
    (hrun : pairVecDivVHCReinsertLin heap nodes lin = .ok state) :
    PairVecDivVHCNodeDivisorIndicesFixed state.nodes := by
  rw [pairVecDivVHCReinsertLin] at hrun
  split at hrun
  next hlin =>
    dsimp only at hrun
    cases hinsert : pairVecDivVHCInsert lin[lin.size - 1] heap nodes with
    | error fault => simp [hinsert] at hrun
    | ok inserted =>
        rcases inserted with ⟨heap', nodes'⟩
        rw [hinsert] at hrun
        have hfixed' := pairVecDivVHCInsert_preserves_divisorIndicesFixed
          lin[lin.size - 1] heap heap' nodes nodes' hfixed hinsert
        exact pairVecDivVHCReinsertLin_preserves_divisorIndicesFixed heap'
          nodes' lin.pop state hfixed' hrun
  next hlin =>
    simp only [Except.ok.injEq] at hrun
    subst state
    exact hfixed
termination_by lin.size
decreasing_by
  simp only [Array.size_pop]
  omega

theorem pairVecDivVHCReinsertLin_preserves_node_invariants
    (degreeLimit resetH quotientSize : Nat) (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode) (lin : Array Nat)
    (quotient divisor : SparsePolyZp) (state : PairVecDivVHCHeapState)
    (hbelow : PairVecDivVHCAllActiveNodesBelow degreeLimit nodes)
    (hdenotes : ∀ (i : Nat) (node : PairVecDivVHCNode),
      nodes[i]? = some node → node.mono ≠ none →
        PairVecDivVHCNodeDenotes quotient divisor node)
    (hfixed : PairVecDivVHCNodeDivisorIndicesFixed nodes)
    (hready : PairVecDivVHCResetReady resetH quotientSize nodes)
    (hlinReady : PairVecDivVHCLinReady lin nodes)
    (hrun : pairVecDivVHCReinsertLin heap nodes lin = .ok state) :
    PairVecDivVHCAllActiveNodesBelow degreeLimit state.nodes ∧
      (∀ (i : Nat) (node : PairVecDivVHCNode),
        state.nodes[i]? = some node → node.mono ≠ none →
          PairVecDivVHCNodeDenotes quotient divisor node) ∧
      PairVecDivVHCNodeDivisorIndicesFixed state.nodes ∧
      PairVecDivVHCResetReady resetH quotientSize state.nodes := by
  rw [pairVecDivVHCReinsertLin] at hrun
  split at hrun
  next hlin =>
    dsimp only at hrun
    cases hinsert : pairVecDivVHCInsert lin[lin.size - 1] heap nodes with
    | error fault => simp [hinsert] at hrun
    | ok inserted =>
        rcases inserted with ⟨heap', nodes'⟩
        rw [hinsert] at hrun
        have hbelow' := pairVecDivVHCInsert_preserves_allActiveNodesBelow
          degreeLimit lin[lin.size - 1] heap heap' nodes nodes' hbelow hinsert
        have hdenotes' := pairVecDivVHCInsert_preserves_denotes
          lin[lin.size - 1] heap heap' nodes nodes' quotient divisor hdenotes
          hinsert
        have hfixed' := pairVecDivVHCInsert_preserves_divisorIndicesFixed
          lin[lin.size - 1] heap heap' nodes nodes' hfixed hinsert
        have hlastMem := pairVecDivVHCLast_mem_toFinset lin hlin
        have houtside : resetH ≤ lin[lin.size - 1] := by
          by_contra hnot
          have hlt : lin[lin.size - 1] < resetH := by omega
          rcases hready.2 lin[lin.size - 1] hlt with
            ⟨readyNode, hreadyNode, hq, hd, hnone⟩
          rcases hlinReady.2 lin[lin.size - 1] (by
            simpa only [List.mem_toFinset] using hlastMem) with
            ⟨activeNode, mono, hactiveNode, hsome⟩
          rw [hreadyNode] at hactiveNode
          simp only [Option.some.injEq] at hactiveNode
          subst activeNode
          rw [hnone] at hsome
          contradiction
        have hready' := pairVecDivVHCInsert_preserves_resetReady
          resetH quotientSize lin[lin.size - 1] heap heap' nodes nodes' hready
          houtside hinsert
        have hlinReady' := hlinReady.pop_after_insert lin heap heap' nodes
          nodes' hlin hinsert
        exact pairVecDivVHCReinsertLin_preserves_node_invariants degreeLimit
          resetH quotientSize heap' nodes' lin.pop quotient divisor state
          hbelow' hdenotes' hfixed' hready' hlinReady' hrun
  next hlin =>
    simp only [Except.ok.injEq] at hrun
    subst state
    exact ⟨hbelow, hdenotes, hfixed, hready⟩
termination_by lin.size
decreasing_by
  simp only [Array.size_pop]
  omega

theorem pairVecDivVHCConsumeNode_preserves_allActiveNodesBelow
    (this : DenseUPolyZp) (p degreeLimit nodeIndex : Nat)
    (currentMono : UMonomial) (k k' : UInt64)
    (nodes nodes' : Array PairVecDivVHCNode) (lin lin' : Array Nat)
    (resetH resetH' : Nat) (next : Option Nat)
    (quotient divisor : SparsePolyZp)
    (hn : nodeIndex < nodes.size)
    (hactive : nodes[nodeIndex].mono = some currentMono)
    (hcanonical : SparsePolyZp.Canonical p quotient)
    (hbelow : PairVecDivVHCAllActiveNodesBelow degreeLimit nodes)
    (hdenotes : ∀ (i : Nat) (node : PairVecDivVHCNode),
      nodes[i]? = some node → node.mono ≠ none →
        PairVecDivVHCNodeDenotes quotient divisor node)
    (hrun : pairVecDivVHCConsumeNode this nodeIndex k nodes lin resetH
      quotient divisor = .ok (k', nodes', lin', resetH', next)) :
    PairVecDivVHCAllActiveNodesBelow degreeLimit nodes' := by
  unfold pairVecDivVHCConsumeNode at hrun
  simp only [hn, ↓reduceDIte] at hrun
  split at hrun <;> try contradiction
  next hq =>
    split at hrun <;> try contradiction
    next hd =>
      split at hrun
      next hadvance =>
        simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
        rcases hrun with ⟨rfl, rfl, rfl, rfl, rfl⟩
        have hnodeGet : nodes[nodeIndex]? = some nodes[nodeIndex] :=
          Array.getElem?_eq_getElem hn
        have hnodeDenotes := hdenotes nodeIndex nodes[nodeIndex] hnodeGet
          (by rw [hactive]; simp)
        have hadvanced := pairVecDivVHCNode_advanced_degree_lt p quotient
          divisor nodes[nodeIndex] currentMono hcanonical hnodeDenotes hactive
          hd hadvance
        have hcurrentBelow := hbelow nodeIndex nodes[nodeIndex] currentMono
          hnodeGet hactive
        intro i node mono hget hmono
        by_cases heq : nodeIndex = i
        · subst i
          rw [Array.getElem?_set_self hn] at hget
          simp only [Option.some.injEq] at hget
          subst node
          simp only at hmono
          have hmonoDegree := congrArg UMonomial.deg
            (Option.some.inj hmono)
          dsimp only at hmonoDegree
          omega
        · rw [Array.getElem?_set_ne hn heq] at hget
          exact hbelow i node mono hget hmono
      next hadvance =>
        split at hrun <;> try contradiction
        next hexhausted =>
          split at hrun <;> try contradiction
          next horder =>
            simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
            rcases hrun with ⟨rfl, rfl, rfl, rfl, rfl⟩
            intro i node mono hget hmono
            by_cases heq : nodeIndex = i
            · subst i
              rw [Array.getElem?_set_self hn] at hget
              simp only [Option.some.injEq] at hget
              subst node
              contradiction
            · rw [Array.getElem?_set_ne hn heq] at hget
              exact hbelow i node mono hget hmono

theorem pairVecDivVHCConsumeNode_preserves_denotes
    (this : DenseUPolyZp) (nodeIndex : Nat) (k k' : UInt64)
    (nodes nodes' : Array PairVecDivVHCNode) (lin lin' : Array Nat)
    (resetH resetH' : Nat) (next : Option Nat)
    (quotient divisor : SparsePolyZp)
    (hdenotes : ∀ (i : Nat) (node : PairVecDivVHCNode),
      nodes[i]? = some node → node.mono ≠ none →
        PairVecDivVHCNodeDenotes quotient divisor node)
    (hrun : pairVecDivVHCConsumeNode this nodeIndex k nodes lin resetH
      quotient divisor = .ok (k', nodes', lin', resetH', next)) :
    ∀ (i : Nat) (node : PairVecDivVHCNode),
      nodes'[i]? = some node → node.mono ≠ none →
        PairVecDivVHCNodeDenotes quotient divisor node := by
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
          intro i node hget hactive
          by_cases heq : nodeIndex = i
          · subst i
            rw [Array.getElem?_set_self hn] at hget
            simp only [Option.some.injEq] at hget
            subst node
            refine ⟨quotient[nodes[nodeIndex].quotientIndex + 1],
              divisor[nodes[nodeIndex].divisorIndex], ?_, ?_, rfl⟩
            · exact Array.getElem?_eq_getElem hadvance
            · exact Array.getElem?_eq_getElem hd
          · rw [Array.getElem?_set_ne hn heq] at hget
            exact hdenotes i node hget hactive
        next hadvance =>
          split at hrun <;> try contradiction
          next hexhausted =>
            split at hrun <;> try contradiction
            next horder =>
              simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
              rcases hrun with ⟨rfl, rfl, rfl, rfl, rfl⟩
              intro i node hget hactive
              by_cases heq : nodeIndex = i
              · subst i
                rw [Array.getElem?_set_self hn] at hget
                simp only [Option.some.injEq] at hget
                subst node
                contradiction
              · rw [Array.getElem?_set_ne hn heq] at hget
                exact hdenotes i node hget hactive

theorem PairVecDivVHCConsumeTrace.products_at_degree
    (this : DenseUPolyZp) (quotient divisor : SparsePolyZp)
    (current : Option Nat) (unvisited : Finset Nat) (degree : Nat)
    (k : UInt64) (nodes : Array PairVecDivVHCNode) (lin : Array Nat)
    (resetH : Nat) (result : PairVecDivVHCBucketResult)
    (products : List (UInt64 × UInt64))
    (hdegree : PairVecDivVHCChainAtDegree current unvisited nodes degree)
    (hdenotes : ∀ (i : Nat) (node : PairVecDivVHCNode),
      nodes[i]? = some node → node.mono ≠ none →
        PairVecDivVHCNodeDenotes quotient divisor node)
    (htrace : PairVecDivVHCConsumeTrace this quotient divisor current unvisited
      k nodes lin resetH result products) :
    ∀ product ∈ products,
      PairVecDivVHCStoredProductAtDegree degree quotient divisor product := by
  induction htrace with
  | done => simp
  | @step nodeIndex unvisited k k' nodes nodes' lin lin' resetH resetH' next
      result products hmem hn hq hd hconsume htail ih =>
      rw [PairVecDivVHCChainAtDegree] at hdegree
      simp only [hmem, ↓reduceDIte] at hdegree
      rcases hdegree with
        ⟨node, mono, hnodeGet, hnodeMono, hmonoDegree, hchainTail⟩
      have hnodeEq : nodes[nodeIndex] = node := by
        rw [Array.getElem?_eq_getElem hn] at hnodeGet
        exact Option.some.inj hnodeGet
      subst node
      have hnodeDenotes := hdenotes nodeIndex nodes[nodeIndex]
        (Array.getElem?_eq_getElem hn) (by rw [hnodeMono]; simp)
      rcases hnodeDenotes with
        ⟨quotientTerm, divisorTerm, hquotientGet, hdivisorGet, hdenotesMono⟩
      have hquotientEq : quotient[nodes[nodeIndex].quotientIndex] =
          quotientTerm := by
        rw [Array.getElem?_eq_getElem hq] at hquotientGet
        exact Option.some.inj hquotientGet
      have hdivisorEq : divisor[nodes[nodeIndex].divisorIndex] =
          divisorTerm := by
        rw [Array.getElem?_eq_getElem hd] at hdivisorGet
        exact Option.some.inj hdivisorGet
      have htailDegree :=
        pairVecDivVHCConsumeNode_preserves_chainAtDegree_tail this nodeIndex
          unvisited degree k k' nodes nodes' lin lin' resetH resetH' next
          quotient divisor hn hchainTail hconsume
      have hdenotes' := pairVecDivVHCConsumeNode_preserves_denotes this
        nodeIndex k k' nodes nodes' lin lin' resetH resetH' next quotient
        divisor hdenotes hconsume
      intro product hproduct
      simp only [List.mem_cons] at hproduct
      rcases hproduct with rfl | hproduct
      · refine ⟨quotientTerm, ?_, divisorTerm, ?_, ?_, ?_, ?_⟩
        · rw [← hquotientEq]
          exact Array.getElem_mem_toList hq
        · rw [← hdivisorEq]
          exact Array.getElem_mem_toList hd
        · exact congrArg (fun term => term.2.val) hquotientEq.symm
        · exact congrArg (fun term => term.2.val) hdivisorEq.symm
        · rw [hnodeMono] at hdenotesMono
          have hmonoEq := Option.some.inj hdenotesMono
          have hdegreeEq := congrArg UMonomial.deg hmonoEq
          change quotientTerm.1.deg + divisorTerm.1.deg = degree
          simpa using hdegreeEq.symm.trans hmonoDegree
      · exact ih htailDegree hdenotes' product hproduct

theorem pairVecDivVHCConsumeRootBucket_coefficient_semantics_at_degree
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (degree : Nat) (heap : Array Nat) (k : UInt64)
    (nodes : Array PairVecDivVHCNode) (lin : Array Nat) (resetH : Nat)
    (quotient divisor : SparsePolyZp) (result : PairVecDivVHCBucketResult)
    (rootMono : UMonomial) (owner : Finset Nat)
    (hheap : 0 < heap.size)
    (hmono : pairVecDivVHCMono heap[0] nodes = .ok rootMono)
    (hdegree : rootMono.deg = degree)
    (howns : PairVecDivVHCChainOwns (some heap[0]) owner nodes)
    (hnextValid : PairVecDivVHCNextValid nodes)
    (hdenotes : ∀ (i : Nat) (node : PairVecDivVHCNode),
      nodes[i]? = some node → node.mono ≠ none →
        PairVecDivVHCNodeDenotes quotient divisor node)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hcanonical : SparsePolyZp.Canonical this._p.toNat quotient)
    (hk : k.toNat < this._p.toNat)
    (hrun : pairVecDivVHCConsumeRootBucket this heap k nodes lin resetH
      quotient divisor = .ok result) :
    ∃ products : List (UInt64 × UInt64),
      PairVecDivVHCConsumeTrace this quotient divisor (some heap[0])
          (Finset.range nodes.size) k nodes lin resetH result products ∧
        (result.coefficient.toNat : ZMod this._p.toNat) =
          (k.toNat : ZMod this._p.toNat) -
            pairVecDivVHCProductsValue this._p.toNat products ∧
        ∀ product ∈ products,
          PairVecDivVHCStoredProductAtDegree degree quotient divisor
            product := by
  have hownerSubset : owner ⊆ Finset.range nodes.size := by
    exact pairVecDivVHCChainOwns_subset_range (some heap[0]) owner nodes howns
  have hvalidOwner := pairVecDivVHCChainOwns_valid (some heap[0]) owner nodes
    howns
  have hvalid := pairVecDivVHCChainValid_mono (some heap[0]) owner
    (Finset.range nodes.size) nodes hvalidOwner hownerSubset
  unfold pairVecDivVHCMono at hmono
  split at hmono <;> try contradiction
  next hn =>
    split at hmono <;> try contradiction
    next mono hnodeMono =>
      simp only [Except.ok.injEq] at hmono
      subst mono
      have hchainDegree :=
        pairVecDivVHCChainValid_atDegree_of_nextValid heap[0]
          (Finset.range nodes.size) nodes degree hvalid hnextValid
          ⟨nodes[heap[0]], rootMono, Array.getElem?_eq_getElem hn,
            hnodeMono, hdegree⟩
      unfold pairVecDivVHCConsumeRootBucket at hrun
      simp only [hheap, ↓reduceDIte] at hrun
      rcases pairVecDivVHCConsumeChain_coefficient_semantics this
          (some heap[0]) (Finset.range nodes.size) k nodes lin resetH quotient
          divisor result hcfg hcanonical hk hrun with
        ⟨products, htrace, hcoefficient⟩
      exact ⟨products, htrace, hcoefficient,
        htrace.products_at_degree this quotient divisor (some heap[0])
          (Finset.range nodes.size) degree k nodes lin resetH result products
          hchainDegree hdenotes⟩

theorem pairVecDivVHCConsumeRootBucket_coefficient_semantics_of_chainDegree
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (degree : Nat) (heap : Array Nat) (k : UInt64)
    (nodes : Array PairVecDivVHCNode) (lin : Array Nat) (resetH : Nat)
    (quotient divisor : SparsePolyZp) (result : PairVecDivVHCBucketResult)
    (owner : Finset Nat) (hheap : 0 < heap.size)
    (howns : PairVecDivVHCChainOwns (some heap[0]) owner nodes)
    (hchainDegree : PairVecDivVHCChainAtDegree (some heap[0]) owner nodes
      degree)
    (hdenotes : ∀ (i : Nat) (node : PairVecDivVHCNode),
      nodes[i]? = some node → node.mono ≠ none →
        PairVecDivVHCNodeDenotes quotient divisor node)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hcanonical : SparsePolyZp.Canonical this._p.toNat quotient)
    (hk : k.toNat < this._p.toNat)
    (hrun : pairVecDivVHCConsumeRootBucket this heap k nodes lin resetH
      quotient divisor = .ok result) :
    ∃ products : List (UInt64 × UInt64),
      PairVecDivVHCConsumeTrace this quotient divisor (some heap[0])
          (Finset.range nodes.size) k nodes lin resetH result products ∧
        (result.coefficient.toNat : ZMod this._p.toNat) =
          (k.toNat : ZMod this._p.toNat) -
            pairVecDivVHCProductsValue this._p.toNat products ∧
        ∀ product ∈ products,
          PairVecDivVHCStoredProductAtDegree degree quotient divisor
            product := by
  have hdegreeRange := pairVecDivVHCChainAtDegree_mono (some heap[0]) owner
    (Finset.range nodes.size) nodes degree hchainDegree
      (pairVecDivVHCChainOwns_subset_range _ _ _ howns)
  unfold pairVecDivVHCConsumeRootBucket at hrun
  simp only [hheap, ↓reduceDIte] at hrun
  rcases pairVecDivVHCConsumeChain_coefficient_semantics this
      (some heap[0]) (Finset.range nodes.size) k nodes lin resetH quotient
      divisor result hcfg hcanonical hk hrun with
    ⟨products, htrace, hcoefficient⟩
  exact ⟨products, htrace, hcoefficient,
    htrace.products_at_degree this quotient divisor (some heap[0])
      (Finset.range nodes.size) degree k nodes lin resetH result products
      hdegreeRange hdenotes⟩

theorem pairVecDivVHCConsumeChain_preserves_denotes
    (this : DenseUPolyZp) (current : Option Nat)
    (unvisited : Finset Nat) (k : UInt64)
    (nodes : Array PairVecDivVHCNode) (lin : Array Nat) (resetH : Nat)
    (quotient divisor : SparsePolyZp) (result : PairVecDivVHCBucketResult)
    (hdenotes : ∀ (i : Nat) (node : PairVecDivVHCNode),
      nodes[i]? = some node → node.mono ≠ none →
        PairVecDivVHCNodeDenotes quotient divisor node)
    (hrun : pairVecDivVHCConsumeChain this current unvisited k nodes lin
      resetH quotient divisor = .ok result) :
    ∀ (i : Nat) (node : PairVecDivVHCNode),
      result.nodes[i]? = some node → node.mono ≠ none →
        PairVecDivVHCNodeDenotes quotient divisor node := by
  cases current with
  | none =>
      rw [pairVecDivVHCConsumeChain] at hrun
      simp only [Except.ok.injEq] at hrun
      subst result
      exact hdenotes
  | some nodeIndex =>
      rw [pairVecDivVHCConsumeChain] at hrun
      split at hrun <;> try contradiction
      next hmem =>
        cases hconsume : pairVecDivVHCConsumeNode this nodeIndex k nodes lin
            resetH quotient divisor with
        | error fault => simp [hconsume] at hrun
        | ok step =>
            rcases step with ⟨k', nodes', lin', resetH', next⟩
            rw [hconsume] at hrun
            have hdenotes' := pairVecDivVHCConsumeNode_preserves_denotes this
              nodeIndex k k' nodes nodes' lin lin' resetH resetH' next quotient
              divisor hdenotes hconsume
            exact pairVecDivVHCConsumeChain_preserves_denotes this next
              (unvisited.erase nodeIndex) k' nodes' lin' resetH' quotient
              divisor result hdenotes' hrun
termination_by unvisited.card
decreasing_by
  exact Finset.card_erase_lt_of_mem (by assumption)

theorem PairVecDivVHCConsumeTrace.products_cover_owner
    (this : DenseUPolyZp) (quotient divisor : SparsePolyZp)
    (current : Option Nat) (unvisited owner : Finset Nat) (k : UInt64)
    (nodes : Array PairVecDivVHCNode) (lin : Array Nat) (resetH : Nat)
    (result : PairVecDivVHCBucketResult)
    (products : List (UInt64 × UInt64))
    (htrace : PairVecDivVHCConsumeTrace this quotient divisor current unvisited
      k nodes lin resetH result products)
    (howns : PairVecDivVHCChainOwns current owner nodes)
    (hdenotes : ∀ (i : Nat) (node : PairVecDivVHCNode),
      nodes[i]? = some node → node.mono ≠ none →
        PairVecDivVHCNodeDenotes quotient divisor node) :
    ∀ i ∈ owner, ∃ node quotientTerm divisorTerm,
      nodes[i]? = some node ∧
        quotient[node.quotientIndex]? = some quotientTerm ∧
        divisor[node.divisorIndex]? = some divisorTerm ∧
        (quotientTerm.2.val, divisorTerm.2.val) ∈ products := by
  induction htrace generalizing owner with
  | done unvisited k nodes lin resetH =>
      simp only [PairVecDivVHCChainOwns] at howns
      subst owner
      simp
  | @step nodeIndex unvisited k k' nodes nodes' lin lin' resetH resetH' next
      result products hmem hn hq hd hconsume htail ih =>
      rw [PairVecDivVHCChainOwns] at howns
      split at howns <;> try contradiction
      next hownerMem =>
        rcases howns with ⟨ownerNode, ownerMono, hownerGet, hownerMono,
          htailOwns⟩
        have hownerEq : ownerNode = nodes[nodeIndex] := by
          rw [Array.getElem?_eq_getElem hn] at hownerGet
          exact (Option.some.inj hownerGet).symm
        subst ownerNode
        intro i hi
        by_cases heq : i = nodeIndex
        · subst i
          refine ⟨nodes[nodeIndex], quotient[nodes[nodeIndex].quotientIndex],
            divisor[nodes[nodeIndex].divisorIndex],
            Array.getElem?_eq_getElem hn, Array.getElem?_eq_getElem hq,
            Array.getElem?_eq_getElem hd, ?_⟩
          exact List.mem_cons_self
        · have hnext := pairVecDivVHCConsumeNode_next this nodeIndex k k'
            nodes nodes' lin lin' resetH resetH' next quotient divisor hn
            hconsume
          have htailOwns' : PairVecDivVHCChainOwns next
              (owner.erase nodeIndex) nodes' := by
            rw [hnext]
            exact pairVecDivVHCChainOwns_congr_on nodes[nodeIndex].next
              (owner.erase nodeIndex) nodes nodes' htailOwns (by
                intro j hj
                exact pairVecDivVHCConsumeNode_get_ne this nodeIndex k k'
                  nodes nodes' lin lin' resetH resetH' next quotient divisor
                  hconsume j (Finset.mem_erase.mp hj).1.symm)
          have hdenotes' := pairVecDivVHCConsumeNode_preserves_denotes this
            nodeIndex k k' nodes nodes' lin lin' resetH resetH' next quotient
            divisor hdenotes hconsume
          rcases ih (owner.erase nodeIndex) htailOwns' hdenotes' i
              (Finset.mem_erase.mpr ⟨heq, hi⟩) with
            ⟨node, quotientTerm, divisorTerm, hnode, hquotient, hdivisor,
              hproduct⟩
          have hsame := pairVecDivVHCConsumeNode_get_ne this nodeIndex k k'
            nodes nodes' lin lin' resetH resetH' next quotient divisor hconsume
            i (fun h => heq h.symm)
          rw [hsame] at hnode
          exact ⟨node, quotientTerm, divisorTerm, hnode, hquotient, hdivisor,
            List.mem_cons_of_mem _ hproduct⟩

/-- A consume trace counts every owned row exactly once, including duplicate
coefficient values.  This is the multiplicity-strengthened form of owner
coverage needed for polynomial coefficients. -/
theorem PairVecDivVHCConsumeTrace.productsValue_eq_owner_sum
    (this : DenseUPolyZp) (p : Nat)
    (quotient divisor : SparsePolyZp) (current : Option Nat)
    (unvisited owner : Finset Nat) (k : UInt64)
    (nodes : Array PairVecDivVHCNode) (lin : Array Nat) (resetH : Nat)
    (result : PairVecDivVHCBucketResult)
    (products : List (UInt64 × UInt64))
    (howns : PairVecDivVHCChainOwns current owner nodes)
    (htrace : PairVecDivVHCConsumeTrace this quotient divisor current unvisited
      k nodes lin resetH result products) :
    pairVecDivVHCProductsValue p products =
      ∑ i ∈ owner, pairVecDivVHCNodeProductValue p nodes quotient divisor i := by
  induction htrace generalizing owner with
  | done =>
      simp only [PairVecDivVHCChainOwns] at howns
      subst owner
      simp [pairVecDivVHCProductsValue]
  | @step nodeIndex unvisited k k' nodes nodes' lin lin' resetH resetH' next
      result products hmem hn hq hd hconsume htail ih =>
      rw [PairVecDivVHCChainOwns] at howns
      split at howns <;> try contradiction
      next hownerMem =>
        rcases howns with ⟨node, mono, hnodeGet, hnodeMono, htailOwns⟩
        have hnodeEq : nodes[nodeIndex] = node := by
          rw [Array.getElem?_eq_getElem hn] at hnodeGet
          exact Option.some.inj hnodeGet
        subst node
        have hnext := pairVecDivVHCConsumeNode_next this nodeIndex k k' nodes
          nodes' lin lin' resetH resetH' next quotient divisor hn hconsume
        have htailOwns' : PairVecDivVHCChainOwns next
            (owner.erase nodeIndex) nodes' := by
          rw [hnext]
          exact pairVecDivVHCChainOwns_congr_on nodes[nodeIndex].next
            (owner.erase nodeIndex) nodes nodes' htailOwns (by
              intro i hi
              exact pairVecDivVHCConsumeNode_get_ne this nodeIndex k k' nodes
                nodes' lin lin' resetH resetH' next quotient divisor hconsume
                i (Finset.mem_erase.mp hi).1.symm)
        have htailValue := ih (owner.erase nodeIndex) htailOwns'
        have htailSame :
            (∑ i ∈ owner.erase nodeIndex,
                pairVecDivVHCNodeProductValue p nodes' quotient divisor i) =
              ∑ i ∈ owner.erase nodeIndex,
                pairVecDivVHCNodeProductValue p nodes quotient divisor i := by
          apply Finset.sum_congr rfl
          intro i hi
          have hne : nodeIndex ≠ i := (Finset.mem_erase.mp hi).1.symm
          unfold pairVecDivVHCNodeProductValue
          rw [pairVecDivVHCConsumeNode_get_ne this nodeIndex k k' nodes nodes'
            lin lin' resetH resetH' next quotient divisor hconsume i hne]
        have hownerEq : owner = insert nodeIndex (owner.erase nodeIndex) := by
          exact (Finset.insert_erase hownerMem).symm
        rw [hownerEq, Finset.sum_insert (by simp)]
        simp only [pairVecDivVHCProductsValue]
        rw [htailValue, htailSame]
        unfold pairVecDivVHCNodeProductValue
        rw [Array.getElem?_eq_getElem hn]
        simp [Array.getElem?_eq_getElem hq,
          Array.getElem?_eq_getElem hd]

theorem pairVecDivVHCConsumeRootBucket_products_complete
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (degree : Nat) (heap : Array Nat) (k : UInt64)
    (nodes : Array PairVecDivVHCNode) (lin : Array Nat) (resetH : Nat)
    (quotient divisor : SparsePolyZp) (result : PairVecDivVHCBucketResult)
    (owner : Finset Nat) (hheap : 0 < heap.size)
    (howns : PairVecDivVHCChainOwns (some heap[0]) owner nodes)
    (hchainDegree : PairVecDivVHCChainAtDegree (some heap[0]) owner nodes
      degree)
    (hdenotes : ∀ (i : Nat) (node : PairVecDivVHCNode),
      nodes[i]? = some node → node.mono ≠ none →
        PairVecDivVHCNodeDenotes quotient divisor node)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hcanonical : SparsePolyZp.Canonical this._p.toNat quotient)
    (hk : k.toNat < this._p.toNat)
    (hrun : pairVecDivVHCConsumeRootBucket this heap k nodes lin resetH
      quotient divisor = .ok result) :
    ∃ products : List (UInt64 × UInt64),
      PairVecDivVHCConsumeTrace this quotient divisor (some heap[0])
          (Finset.range nodes.size) k nodes lin resetH result products ∧
        (result.coefficient.toNat : ZMod this._p.toNat) =
          (k.toNat : ZMod this._p.toNat) -
            pairVecDivVHCProductsValue this._p.toNat products ∧
        (∀ product ∈ products,
          PairVecDivVHCStoredProductAtDegree degree quotient divisor product) ∧
        ∀ i ∈ owner, ∃ node quotientTerm divisorTerm,
          nodes[i]? = some node ∧
            quotient[node.quotientIndex]? = some quotientTerm ∧
            divisor[node.divisorIndex]? = some divisorTerm ∧
            (quotientTerm.2.val, divisorTerm.2.val) ∈ products := by
  rcases pairVecDivVHCConsumeRootBucket_coefficient_semantics_of_chainDegree
      this degree heap k nodes lin resetH quotient divisor result owner hheap
      howns hchainDegree hdenotes hcfg hcanonical hk hrun with
    ⟨products, htrace, hcoefficient, hsound⟩
  refine ⟨products, htrace, hcoefficient, hsound, ?_⟩
  exact htrace.products_cover_owner this quotient divisor (some heap[0])
    (Finset.range nodes.size) owner k nodes lin resetH result products howns
    hdenotes

/-- Every row in every heap-owned chain at `degree` contributes its concrete
coefficient pair to `products`. -/
def PairVecDivVHCProductsCoverDegreeOwners (degree : Nat)
    (heap : Array Nat) (owners : Nat → Finset Nat)
    (nodes : Array PairVecDivVHCNode)
    (quotient divisor : SparsePolyZp)
    (products : List (UInt64 × UInt64)) : Prop :=
  ∀ (slot head : Nat), heap[slot]? = some head →
    PairVecDivVHCChainAtDegree (some head) (owners head) nodes degree →
    ∀ i ∈ owners head, ∃ node quotientTerm divisorTerm,
      nodes[i]? = some node ∧
        quotient[node.quotientIndex]? = some quotientTerm ∧
        divisor[node.divisorIndex]? = some divisorTerm ∧
        (quotientTerm.2.val, divisorTerm.2.val) ∈ products

theorem pairVecDivVHCConsumeChain_preserves_cursorPrefixAbove
    (this : DenseUPolyZp) (degree : Nat) (current : Option Nat)
    (unvisited : Finset Nat) (k : UInt64)
    (nodes : Array PairVecDivVHCNode) (lin : Array Nat) (resetH : Nat)
    (quotient divisor : SparsePolyZp) (result : PairVecDivVHCBucketResult)
    (hdegree : PairVecDivVHCChainAtDegree current unvisited nodes degree)
    (hprefix : PairVecDivVHCCursorPrefixAbove degree nodes quotient divisor)
    (hdenotes : ∀ (i : Nat) (node : PairVecDivVHCNode),
      nodes[i]? = some node → node.mono ≠ none →
        PairVecDivVHCNodeDenotes quotient divisor node)
    (hrun : pairVecDivVHCConsumeChain this current unvisited k nodes lin
      resetH quotient divisor = .ok result) :
    PairVecDivVHCCursorPrefixAbove degree result.nodes quotient divisor := by
  cases current with
  | none =>
      rw [pairVecDivVHCConsumeChain] at hrun
      simp only [Except.ok.injEq] at hrun
      subst result
      exact hprefix
  | some nodeIndex =>
      rw [PairVecDivVHCChainAtDegree] at hdegree
      split at hdegree <;> try contradiction
      next hmem =>
        rcases hdegree with
          ⟨node, mono, hget, hmono, hmonoDegree, htailDegree⟩
        rw [pairVecDivVHCConsumeChain] at hrun
        simp only [hmem, ↓reduceDIte] at hrun
        cases hconsume : pairVecDivVHCConsumeNode this nodeIndex k nodes lin
            resetH quotient divisor with
        | error fault => simp [hconsume] at hrun
        | ok step =>
            rcases step with ⟨k', nodes', lin', resetH', next⟩
            rw [hconsume] at hrun
            have hnodeDenotes := hdenotes nodeIndex node hget (by
              rw [hmono]
              simp)
            have hmonoEq : mono = ⟨degree⟩ := by
              cases mono with
              | mk monoDegree =>
                  simp only [UMonomial.deg] at hmonoDegree
                  subst monoDegree
                  rfl
            have hprefix' :=
              pairVecDivVHCConsumeNode_preserves_cursorPrefixAbove this degree
                nodeIndex k k' nodes nodes' lin lin' resetH resetH' next
                quotient divisor node hprefix hget hnodeDenotes (by
                  simpa [hmonoEq] using hmono) hconsume
            have hn : nodeIndex < nodes.size := by
              by_contra hnot
              rw [Array.getElem?_eq_none (by omega)] at hget
              contradiction
            have hnodeEq : nodes[nodeIndex] = node := by
              rw [Array.getElem?_eq_getElem hn] at hget
              exact Option.some.inj hget
            have hnext := pairVecDivVHCConsumeNode_next this nodeIndex k k'
              nodes nodes' lin lin' resetH resetH' next quotient divisor
              hn hconsume
            have htailDegree' : PairVecDivVHCChainAtDegree next
                (unvisited.erase nodeIndex) nodes' degree := by
              rw [hnext, hnodeEq]
              exact pairVecDivVHCChainAtDegree_congr_on node.next
                (unvisited.erase nodeIndex) nodes nodes' degree htailDegree (by
                  intro i hi
                  exact pairVecDivVHCConsumeNode_get_ne this nodeIndex k k'
                    nodes nodes' lin lin' resetH resetH' next quotient divisor
                    hconsume i (Finset.mem_erase.mp hi).1.symm)
            have hdenotes' := pairVecDivVHCConsumeNode_preserves_denotes this
              nodeIndex k k' nodes nodes' lin lin' resetH resetH' next quotient
              divisor hdenotes hconsume
            exact pairVecDivVHCConsumeChain_preserves_cursorPrefixAbove this
              degree next (unvisited.erase nodeIndex) k' nodes' lin' resetH'
              quotient divisor result htailDegree' hprefix' hdenotes' hrun
termination_by unvisited.card
decreasing_by
  exact Finset.card_erase_lt_of_mem (by assumption)

theorem pairVecDivVHCConsumeRootBucket_preserves_cursorPrefixAbove
    (this : DenseUPolyZp) (degree : Nat) (heap : Array Nat) (k : UInt64)
    (nodes : Array PairVecDivVHCNode) (lin : Array Nat) (resetH : Nat)
    (quotient divisor : SparsePolyZp) (result : PairVecDivVHCBucketResult)
    (hheap : 0 < heap.size)
    (hdegree : PairVecDivVHCChainAtDegree (some heap[0])
      (Finset.range nodes.size) nodes degree)
    (hprefix : PairVecDivVHCCursorPrefixAbove degree nodes quotient divisor)
    (hdenotes : ∀ (i : Nat) (node : PairVecDivVHCNode),
      nodes[i]? = some node → node.mono ≠ none →
        PairVecDivVHCNodeDenotes quotient divisor node)
    (hrun : pairVecDivVHCConsumeRootBucket this heap k nodes lin resetH
      quotient divisor = .ok result) :
    PairVecDivVHCCursorPrefixAbove degree result.nodes quotient divisor := by
  unfold pairVecDivVHCConsumeRootBucket at hrun
  simp only [hheap, ↓reduceDIte] at hrun
  exact pairVecDivVHCConsumeChain_preserves_cursorPrefixAbove this degree
    (some heap[0]) (Finset.range nodes.size) k nodes lin resetH quotient
    divisor result hdegree hprefix hdenotes hrun

theorem pairVecDivVHCConsumeRootBucket_preserves_denotes
    (this : DenseUPolyZp) (heap : Array Nat) (k : UInt64)
    (nodes : Array PairVecDivVHCNode) (lin : Array Nat) (resetH : Nat)
    (quotient divisor : SparsePolyZp) (result : PairVecDivVHCBucketResult)
    (hdenotes : ∀ (i : Nat) (node : PairVecDivVHCNode),
      nodes[i]? = some node → node.mono ≠ none →
        PairVecDivVHCNodeDenotes quotient divisor node)
    (hrun : pairVecDivVHCConsumeRootBucket this heap k nodes lin resetH
      quotient divisor = .ok result) :
    ∀ (i : Nat) (node : PairVecDivVHCNode),
      result.nodes[i]? = some node → node.mono ≠ none →
        PairVecDivVHCNodeDenotes quotient divisor node := by
  unfold pairVecDivVHCConsumeRootBucket at hrun
  split at hrun <;> try contradiction
  next hheap =>
    exact pairVecDivVHCConsumeChain_preserves_denotes this (some heap[0])
      (Finset.range nodes.size) k nodes lin resetH quotient divisor result
      hdenotes hrun

theorem pairVecDivVHCConsumeNode_preserves_divisorIndicesFixed
    (this : DenseUPolyZp) (nodeIndex : Nat) (k k' : UInt64)
    (nodes nodes' : Array PairVecDivVHCNode) (lin lin' : Array Nat)
    (resetH resetH' : Nat) (next : Option Nat)
    (quotient divisor : SparsePolyZp)
    (hfixed : PairVecDivVHCNodeDivisorIndicesFixed nodes)
    (hrun : pairVecDivVHCConsumeNode this nodeIndex k nodes lin resetH
      quotient divisor = .ok (k', nodes', lin', resetH', next)) :
    PairVecDivVHCNodeDivisorIndicesFixed nodes' := by
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
          exact PairVecDivVHCNodeDivisorIndicesFixed.set nodes nodeIndex _ hn
            hfixed rfl
        next hadvance =>
          split at hrun <;> try contradiction
          next hexhausted =>
            split at hrun <;> try contradiction
            next horder =>
              simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
              rcases hrun with ⟨rfl, rfl, rfl, rfl, rfl⟩
              exact PairVecDivVHCNodeDivisorIndicesFixed.set nodes nodeIndex _ hn
                hfixed rfl

theorem pairVecDivVHCConsumeNode_preserves_chain_tail
    (this : DenseUPolyZp) (nodeIndex : Nat) (unvisited : Finset Nat)
    (k k' : UInt64) (nodes nodes' : Array PairVecDivVHCNode)
    (lin lin' : Array Nat) (resetH resetH' : Nat) (next : Option Nat)
    (quotient divisor : SparsePolyZp)
    (hn : nodeIndex < nodes.size)
    (htail : PairVecDivVHCChainValid nodes[nodeIndex].next
      (unvisited.erase nodeIndex) nodes)
    (hrun : pairVecDivVHCConsumeNode this nodeIndex k nodes lin resetH
      quotient divisor = .ok (k', nodes', lin', resetH', next)) :
    PairVecDivVHCChainValid next (unvisited.erase nodeIndex) nodes' := by
  unfold pairVecDivVHCConsumeNode at hrun
  simp only [hn, ↓reduceDIte] at hrun
  split at hrun <;> try contradiction
  next hq =>
    split at hrun <;> try contradiction
    next hd =>
      split at hrun
      next hadvance =>
        simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
        rcases hrun with ⟨rfl, rfl, rfl, rfl, rfl⟩
        exact pairVecDivVHCChainValid_set_of_not_mem nodes[nodeIndex].next
          (unvisited.erase nodeIndex) nodes nodeIndex
          { nodes[nodeIndex] with
            quotientIndex := nodes[nodeIndex].quotientIndex + 1
            mono := some ⟨quotient[nodes[nodeIndex].quotientIndex + 1].1.deg +
              divisor[nodes[nodeIndex].divisorIndex].1.deg⟩ }
          hn (by simp) htail
      next hadvance =>
        split at hrun <;> try contradiction
        next hexhausted =>
          split at hrun <;> try contradiction
          next horder =>
            simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
            rcases hrun with ⟨rfl, rfl, rfl, rfl, rfl⟩
            exact pairVecDivVHCChainValid_set_of_not_mem nodes[nodeIndex].next
              (unvisited.erase nodeIndex) nodes nodeIndex
              { nodes[nodeIndex] with
                quotientIndex := nodes[nodeIndex].quotientIndex + 1
                mono := none
                next := none }
              hn (by simp) htail

theorem pairVecDivVHCConsumeNode_exhausted_preserves_resetReady
    (this : DenseUPolyZp) (nodeIndex : Nat) (k k' : UInt64)
    (nodes nodes' : Array PairVecDivVHCNode) (lin lin' : Array Nat)
    (resetH resetH' : Nat) (next : Option Nat)
    (quotient divisor : SparsePolyZp)
    (hn : nodeIndex < nodes.size)
    (hdivisorIndex : nodes[nodeIndex].divisorIndex = nodeIndex + 1)
    (hready : PairVecDivVHCResetReady resetH quotient.size nodes)
    (hlin : lin'.size = lin.size)
    (hrun : pairVecDivVHCConsumeNode this nodeIndex k nodes lin resetH
      quotient divisor = .ok (k', nodes', lin', resetH', next)) :
    PairVecDivVHCResetReady resetH' quotient.size nodes' := by
  unfold pairVecDivVHCConsumeNode at hrun
  simp only [hn, ↓reduceDIte] at hrun
  split at hrun <;> try contradiction
  next hq =>
    split at hrun <;> try contradiction
    next hd =>
      split at hrun
      next hadvance =>
        simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
        rcases hrun with ⟨rfl, rfl, rfl, rfl, rfl⟩
        simp at hlin
      next hadvance =>
        split at hrun <;> try contradiction
        next hexhausted =>
          split at hrun <;> try contradiction
          next horder =>
            simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
            rcases hrun with ⟨rfl, rfl, rfl, rfl, rfl⟩
            subst nodeIndex
            refine ⟨?_, ?_⟩
            · simpa using hn
            · intro i hi
              by_cases hold : i < resetH
              · rcases hready.2 i hold with ⟨node, hget, hqIndex,
                  hdIndex, hmono⟩
                refine ⟨node, ?_, hqIndex, hdIndex, hmono⟩
                rw [Array.getElem?_set_ne hn (by omega)]
                exact hget
              · have heq : i = resetH := by omega
                subst i
                refine ⟨{ nodes[resetH] with
                    quotientIndex := nodes[resetH].quotientIndex + 1
                    mono := none
                    next := none }, ?_, ?_, ?_, rfl⟩
                · exact Array.getElem?_set_self hn
                · exact hexhausted
                · exact hdivisorIndex

theorem pairVecDivVHCConsumeNode_preserves_resetReady
    (this : DenseUPolyZp) (nodeIndex : Nat) (currentMono : UMonomial)
    (k k' : UInt64) (nodes nodes' : Array PairVecDivVHCNode)
    (lin lin' : Array Nat) (resetH resetH' : Nat) (next : Option Nat)
    (quotient divisor : SparsePolyZp)
    (hn : nodeIndex < nodes.size)
    (hactive : nodes[nodeIndex].mono = some currentMono)
    (hready : PairVecDivVHCResetReady resetH quotient.size nodes)
    (hfixed : PairVecDivVHCNodeDivisorIndicesFixed nodes)
    (hrun : pairVecDivVHCConsumeNode this nodeIndex k nodes lin resetH
      quotient divisor = .ok (k', nodes', lin', resetH', next)) :
    PairVecDivVHCResetReady resetH' quotient.size nodes' := by
  have hnotPrefix : ¬ nodeIndex < resetH := by
    intro hlt
    rcases hready.2 nodeIndex hlt with ⟨node, hget, hqIndex,
      hdIndex, hmono⟩
    rw [Array.getElem?_eq_getElem hn] at hget
    simp only [Option.some.injEq] at hget
    subst node
    rw [hactive] at hmono
    contradiction
  unfold pairVecDivVHCConsumeNode at hrun
  simp only [hn, ↓reduceDIte] at hrun
  split at hrun <;> try contradiction
  next hq =>
    split at hrun <;> try contradiction
    next hd =>
      split at hrun
      next hadvance =>
        simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
        rcases hrun with ⟨rfl, rfl, rfl, rfl, rfl⟩
        refine ⟨by simpa using hready.1, ?_⟩
        intro i hi
        rcases hready.2 i hi with ⟨node, hget, hqIndex,
          hdIndex, hmono⟩
        refine ⟨node, ?_, hqIndex, hdIndex, hmono⟩
        rw [Array.getElem?_set_ne hn (by omega)]
        exact hget
      next hadvance =>
        split at hrun <;> try contradiction
        next hexhausted =>
          split at hrun <;> try contradiction
          next horder =>
            simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
            rcases hrun with ⟨rfl, rfl, rfl, rfl, rfl⟩
            have hdivisorIndex := hfixed nodeIndex nodes[nodeIndex]
              (Array.getElem?_eq_getElem hn)
            subst nodeIndex
            refine ⟨by simpa using hn, ?_⟩
            intro i hi
            by_cases hold : i < resetH
            · rcases hready.2 i hold with ⟨node, hget, hqIndex,
                hdIndex, hmono⟩
              refine ⟨node, ?_, hqIndex, hdIndex, hmono⟩
              rw [Array.getElem?_set_ne hn (by omega)]
              exact hget
            · have heq : i = resetH := by omega
              subst i
              refine ⟨{ nodes[resetH] with
                  quotientIndex := nodes[resetH].quotientIndex + 1
                  mono := none
                  next := none }, Array.getElem?_set_self hn,
                hexhausted, hdivisorIndex, rfl⟩

theorem pairVecDivVHCConsumeChain_preserves_node_invariants
    (this : DenseUPolyZp) (p degreeLimit : Nat)
    (current : Option Nat) (unvisited : Finset Nat) (k : UInt64)
    (nodes : Array PairVecDivVHCNode) (lin : Array Nat) (resetH : Nat)
    (quotient divisor : SparsePolyZp) (result : PairVecDivVHCBucketResult)
    (hcanonical : SparsePolyZp.Canonical p quotient)
    (hbelow : PairVecDivVHCAllActiveNodesBelow degreeLimit nodes)
    (hdenotes : ∀ (i : Nat) (node : PairVecDivVHCNode),
      nodes[i]? = some node → node.mono ≠ none →
        PairVecDivVHCNodeDenotes quotient divisor node)
    (hfixed : PairVecDivVHCNodeDivisorIndicesFixed nodes)
    (hready : PairVecDivVHCResetReady resetH quotient.size nodes)
    (hchain : PairVecDivVHCChainValid current unvisited nodes)
    (hrun : pairVecDivVHCConsumeChain this current unvisited k nodes lin
      resetH quotient divisor = .ok result) :
    PairVecDivVHCAllActiveNodesBelow degreeLimit result.nodes ∧
      (∀ (i : Nat) (node : PairVecDivVHCNode),
        result.nodes[i]? = some node → node.mono ≠ none →
          PairVecDivVHCNodeDenotes quotient divisor node) ∧
      PairVecDivVHCNodeDivisorIndicesFixed result.nodes ∧
      PairVecDivVHCResetReady result.resetH quotient.size result.nodes := by
  cases current with
  | none =>
      rw [pairVecDivVHCConsumeChain] at hrun
      simp only [Except.ok.injEq] at hrun
      subst result
      exact ⟨hbelow, hdenotes, hfixed, hready⟩
  | some nodeIndex =>
      rw [pairVecDivVHCConsumeChain] at hrun
      split at hrun <;> try contradiction
      next hmem =>
        rw [PairVecDivVHCChainValid] at hchain
        simp only [hmem, ↓reduceDIte] at hchain
        rcases hchain with ⟨node, mono, hget, hactive, htail⟩
        have hn : nodeIndex < nodes.size := by
          by_contra hnot
          rw [Array.getElem?_eq_none (by omega)] at hget
          contradiction
        rw [Array.getElem?_eq_getElem hn] at hget
        simp only [Option.some.injEq] at hget
        subst node
        cases hconsume : pairVecDivVHCConsumeNode this nodeIndex k nodes lin
            resetH quotient divisor with
        | error fault => simp [hconsume] at hrun
        | ok step =>
            rcases step with ⟨k', nodes', lin', resetH', next⟩
            rw [hconsume] at hrun
            have hbelow' :=
              pairVecDivVHCConsumeNode_preserves_allActiveNodesBelow this p
                degreeLimit nodeIndex mono k k' nodes nodes' lin lin' resetH
                resetH' next quotient divisor hn hactive hcanonical hbelow
                hdenotes hconsume
            have hdenotes' := pairVecDivVHCConsumeNode_preserves_denotes this
              nodeIndex k k' nodes nodes' lin lin' resetH resetH' next
              quotient divisor hdenotes hconsume
            have hfixed' :=
              pairVecDivVHCConsumeNode_preserves_divisorIndicesFixed this
                nodeIndex k k' nodes nodes' lin lin' resetH resetH' next
                quotient divisor hfixed hconsume
            have hready' := pairVecDivVHCConsumeNode_preserves_resetReady this
              nodeIndex mono k k' nodes nodes' lin lin' resetH resetH' next
              quotient divisor hn hactive hready hfixed hconsume
            have htail' := pairVecDivVHCConsumeNode_preserves_chain_tail this
              nodeIndex unvisited k k' nodes nodes' lin lin' resetH resetH'
              next quotient divisor hn htail hconsume
            exact pairVecDivVHCConsumeChain_preserves_node_invariants this p
              degreeLimit next (unvisited.erase nodeIndex) k' nodes' lin'
              resetH' quotient divisor result hcanonical hbelow' hdenotes'
              hfixed' hready' htail' hrun
termination_by unvisited.card
decreasing_by
  exact Finset.card_erase_lt_of_mem (by assumption)

theorem pairVecDivVHCConsumeRootBucket_preserves_node_invariants
    (this : DenseUPolyZp) (p degreeLimit : Nat)
    (heap : Array Nat) (k : UInt64) (nodes : Array PairVecDivVHCNode)
    (lin : Array Nat) (resetH : Nat) (quotient divisor : SparsePolyZp)
    (result : PairVecDivVHCBucketResult)
    (hheap : 0 < heap.size)
    (hcanonical : SparsePolyZp.Canonical p quotient)
    (hbelow : PairVecDivVHCAllActiveNodesBelow degreeLimit nodes)
    (hdenotes : ∀ (i : Nat) (node : PairVecDivVHCNode),
      nodes[i]? = some node → node.mono ≠ none →
        PairVecDivVHCNodeDenotes quotient divisor node)
    (hfixed : PairVecDivVHCNodeDivisorIndicesFixed nodes)
    (hready : PairVecDivVHCResetReady resetH quotient.size nodes)
    (hchain : PairVecDivVHCChainValid (some heap[0])
      (Finset.range nodes.size) nodes)
    (hrun : pairVecDivVHCConsumeRootBucket this heap k nodes lin resetH
      quotient divisor = .ok result) :
    PairVecDivVHCAllActiveNodesBelow degreeLimit result.nodes ∧
      (∀ (i : Nat) (node : PairVecDivVHCNode),
        result.nodes[i]? = some node → node.mono ≠ none →
          PairVecDivVHCNodeDenotes quotient divisor node) ∧
      PairVecDivVHCNodeDivisorIndicesFixed result.nodes ∧
      PairVecDivVHCResetReady result.resetH quotient.size result.nodes := by
  unfold pairVecDivVHCConsumeRootBucket at hrun
  simp only [hheap, ↓reduceDIte] at hrun
  exact pairVecDivVHCConsumeChain_preserves_node_invariants this p degreeLimit
    (some heap[0]) (Finset.range nodes.size) k nodes lin resetH quotient
    divisor result hcanonical hbelow hdenotes hfixed hready hchain hrun

theorem pairVecDivVHCConsumeRootBucket_unvisited (this : DenseUPolyZp)
    (heap : Array Nat) (k : UInt64) (nodes : Array PairVecDivVHCNode)
    (lin : Array Nat) (resetH : Nat) (quotient divisor : SparsePolyZp)
    (result : PairVecDivVHCBucketResult)
    (hheap : 0 < heap.size)
    (hrun : pairVecDivVHCConsumeRootBucket this heap k nodes lin resetH
      quotient divisor = .ok result) :
    result.unvisited ⊆ Finset.range nodes.size ∧
      ∀ i ∈ result.unvisited, result.nodes[i]? = nodes[i]? := by
  unfold pairVecDivVHCConsumeRootBucket at hrun
  simp only [hheap, ↓reduceDIte] at hrun
  exact pairVecDivVHCConsumeChain_unvisited this (some heap[0])
    (Finset.range nodes.size) k nodes lin resetH quotient divisor result hrun

theorem pairVecDivVHCConsumeRootBucket_preserves_disjoint_chain
    (this : DenseUPolyZp) (heap : Array Nat) (k : UInt64)
    (nodes : Array PairVecDivVHCNode) (lin : Array Nat) (resetH : Nat)
    (quotient divisor : SparsePolyZp) (result : PairVecDivVHCBucketResult)
    (rootOwner otherOwner : Finset Nat) (otherCurrent : Option Nat)
    (hheap : 0 < heap.size)
    (hrootOwns : PairVecDivVHCChainOwns (some heap[0]) rootOwner nodes)
    (hotherOwns : PairVecDivVHCChainOwns otherCurrent otherOwner nodes)
    (hdisjoint : Disjoint otherOwner rootOwner)
    (hrun : pairVecDivVHCConsumeRootBucket this heap k nodes lin resetH
      quotient divisor = .ok result) :
    otherOwner ⊆ result.unvisited ∧
      PairVecDivVHCChainOwns otherCurrent otherOwner result.nodes := by
  unfold pairVecDivVHCConsumeRootBucket at hrun
  simp only [hheap, ↓reduceDIte] at hrun
  have hpreserved := pairVecDivVHCConsumeChain_preserves_disjoint_owner this
    (some heap[0]) rootOwner otherOwner (Finset.range nodes.size) k nodes lin
    resetH quotient divisor result hrootOwns
    (pairVecDivVHCChainOwns_subset_range _ _ _ hrootOwns)
    (pairVecDivVHCChainOwns_subset_range _ _ _ hotherOwns) hdisjoint hrun
  exact ⟨hpreserved.1, pairVecDivVHCChainOwns_congr_on otherCurrent
    otherOwner nodes result.nodes hotherOwns hpreserved.2⟩

theorem pairVecDivVHCConsumeRootBucket_preserves_disjoint_chainAtDegree
    (this : DenseUPolyZp) (heap : Array Nat) (k : UInt64)
    (nodes : Array PairVecDivVHCNode) (lin : Array Nat) (resetH : Nat)
    (quotient divisor : SparsePolyZp) (result : PairVecDivVHCBucketResult)
    (rootOwner otherOwner : Finset Nat) (otherCurrent : Option Nat)
    (degree : Nat) (hheap : 0 < heap.size)
    (hrootOwns : PairVecDivVHCChainOwns (some heap[0]) rootOwner nodes)
    (hotherOwns : PairVecDivVHCChainOwns otherCurrent otherOwner nodes)
    (hotherDegree : PairVecDivVHCChainAtDegree otherCurrent otherOwner nodes
      degree)
    (hdisjoint : Disjoint otherOwner rootOwner)
    (hrun : pairVecDivVHCConsumeRootBucket this heap k nodes lin resetH
      quotient divisor = .ok result) :
    (∀ i ∈ otherOwner, result.nodes[i]? = nodes[i]?) ∧
      PairVecDivVHCChainAtDegree otherCurrent otherOwner result.nodes
        degree := by
  unfold pairVecDivVHCConsumeRootBucket at hrun
  simp only [hheap, ↓reduceDIte] at hrun
  have hpreserved := pairVecDivVHCConsumeChain_preserves_disjoint_owner this
    (some heap[0]) rootOwner otherOwner (Finset.range nodes.size) k nodes lin
    resetH quotient divisor result hrootOwns
    (pairVecDivVHCChainOwns_subset_range _ _ _ hrootOwns)
    (pairVecDivVHCChainOwns_subset_range _ _ _ hotherOwns) hdisjoint hrun
  exact ⟨hpreserved.2, pairVecDivVHCChainAtDegree_congr_on otherCurrent
    otherOwner nodes result.nodes degree hotherDegree hpreserved.2⟩

/-- Consuming the real root bucket preserves every other heap bucket's exact
chain ownership.  This is the heap-wide lifting of the local destructive
update isolation theorem. -/
theorem pairVecDivVHCConsumeRootBucket_preserves_other_heap_chains
    (this : DenseUPolyZp) (heap : Array Nat) (k : UInt64)
    (nodes : Array PairVecDivVHCNode) (lin : Array Nat) (resetH : Nat)
    (quotient divisor : SparsePolyZp) (result : PairVecDivVHCBucketResult)
    (owners : Nat → Finset Nat)
    (hheap : 0 < heap.size)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hrun : pairVecDivVHCConsumeRootBucket this heap k nodes lin resetH
      quotient divisor = .ok result) :
    ∀ (slot head : Nat), heap[slot]? = some head → head ≠ heap[0] →
      PairVecDivVHCChainOwns (some head) (owners head) result.nodes := by
  intro slot head hget hne
  have hslot : slot < heap.size := by
    by_contra hnot
    rw [Array.getElem?_eq_none (by omega)] at hget
    contradiction
  have hrootOwns := pairVecDivVHCHeapChainOwnership_root_owns heap owners nodes
    hownership hheap
  have hotherOwns := hownership.1 slot head hget
  have hhead : heap[slot] = head := by
    rw [Array.getElem?_eq_getElem hslot] at hget
    exact Option.some.inj hget
  have hdisjoint := pairVecDivVHCHeapChainOwnership_root_disjoint heap owners
    nodes hownership hheap slot hslot (by simpa [hhead] using hne)
  exact (pairVecDivVHCConsumeRootBucket_preserves_disjoint_chain this heap k
    nodes lin resetH quotient divisor result (owners heap[0])
    (owners head) (some head) hheap hrootOwns hotherOwns (by simpa [hhead]
    using hdisjoint)
    hrun).2

theorem pairVecDivVHCConsumeRootBucket_nonroot_mono_iff
    (this : DenseUPolyZp) (heap : Array Nat) (k : UInt64)
    (nodes : Array PairVecDivVHCNode) (lin : Array Nat) (resetH : Nat)
    (quotient divisor : SparsePolyZp) (result : PairVecDivVHCBucketResult)
    (owners : Nat → Finset Nat) (hheap : 0 < heap.size)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hrun : pairVecDivVHCConsumeRootBucket this heap k nodes lin resetH
      quotient divisor = .ok result) :
    ∀ (slot head : Nat), 0 < slot → heap[slot]? = some head →
      ∀ mono, pairVecDivVHCMono head nodes = .ok mono ↔
        pairVecDivVHCMono head result.nodes = .ok mono := by
  intro slot head hslotPos hget mono
  have hslot : slot < heap.size := by
    by_contra hnot
    rw [Array.getElem?_eq_none (by omega)] at hget
    contradiction
  have hrootGet : heap[0]? = some heap[0] :=
    Array.getElem?_eq_getElem hheap
  have hne : head ≠ heap[0] := by
    intro heq
    subst head
    exact (by omega : slot ≠ 0) (hownership.2.1 slot 0 heap[0] hget hrootGet)
  have hrootOwns := pairVecDivVHCHeapChainOwnership_root_owns heap owners nodes
    hownership hheap
  have hotherOwns := hownership.1 slot head hget
  have hdisjoint : Disjoint (owners head) (owners heap[0]) :=
    hownership.2.2 slot 0 head heap[0] hget hrootGet hne
  have hpreserved := pairVecDivVHCConsumeRootBucket_preserves_disjoint_chain
    this heap k nodes lin resetH quotient divisor result (owners heap[0])
    (owners head) (some head) hheap hrootOwns hotherOwns hdisjoint hrun
  have hunvisited := pairVecDivVHCConsumeRootBucket_unvisited this heap k nodes
    lin resetH quotient divisor result hheap hrun
  have hheadMem := pairVecDivVHCChainOwns_head_mem head (owners head) nodes
    hotherOwns
  have hsame : result.nodes[head]? = nodes[head]? :=
    hunvisited.2 head (hpreserved.1 hheadMem)
  constructor
  · intro hmono
    rcases (pairVecDivVHCMono_eq_ok_iff head nodes mono).mp hmono with
      ⟨node, hnode, hactive⟩
    apply (pairVecDivVHCMono_eq_ok_iff head result.nodes mono).mpr
    refine ⟨node, ?_, hactive⟩
    rw [hsame]
    exact hnode
  · intro hmono
    rcases (pairVecDivVHCMono_eq_ok_iff head result.nodes mono).mp hmono with
      ⟨node, hnode, hactive⟩
    apply (pairVecDivVHCMono_eq_ok_iff head nodes mono).mpr
    refine ⟨node, ?_, hactive⟩
    rw [← hsame]
    exact hnode

/-- Consuming the root bucket leaves every old heap edge whose parent is not
the root ordered in the updated node array.  The root itself is intentionally
excluded: its chain has just been consumed and may no longer be active. -/
theorem pairVecDivVHCConsumeRootBucket_preserves_nonroot_heap_order
    (this : DenseUPolyZp) (heap : Array Nat) (k : UInt64)
    (nodes : Array PairVecDivVHCNode) (lin : Array Nat) (resetH : Nat)
    (quotient divisor : SparsePolyZp) (result : PairVecDivVHCBucketResult)
    (owners : Nat → Finset Nat) (hheap : 0 < heap.size)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hordered : PairVecDivVHCHeapOrdered heap nodes)
    (hrun : pairVecDivVHCConsumeRootBucket this heap k nodes lin resetH
      quotient divisor = .ok result) :
    ∀ (child : Nat), child < heap.size → 0 < child →
      0 < pairVecDivVHCParent child →
      ∀ childHead parentHead childMono parentMono,
        heap[child]? = some childHead →
        heap[pairVecDivVHCParent child]? = some parentHead →
        pairVecDivVHCMono childHead result.nodes = .ok childMono →
        pairVecDivVHCMono parentHead result.nodes = .ok parentMono →
        childMono.deg ≤ parentMono.deg := by
  intro child hchild hchildPos hparentPos childHead parentHead childMono
    parentMono hchildGet hparentGet hchildMono hparentMono
  have hchildSource :=
    (pairVecDivVHCConsumeRootBucket_nonroot_mono_iff this heap k nodes lin
      resetH quotient divisor result owners hheap hownership hrun child
      childHead hchildPos hchildGet childMono).mpr hchildMono
  have hparentSlot : pairVecDivVHCParent child < heap.size := by
    have hparentLt : pairVecDivVHCParent child < child :=
      pairVecDivVHCParent_lt child hchildPos
    exact Nat.lt_trans hparentLt hchild
  have hparentSource :=
    (pairVecDivVHCConsumeRootBucket_nonroot_mono_iff this heap k nodes lin
      resetH quotient divisor result owners hheap hownership hrun
      (pairVecDivVHCParent child) parentHead hparentPos hparentGet
      parentMono).mpr hparentMono
  exact hordered.degreesUpTo heap nodes heap.size (Nat.le_refl _) child hchild
    hchildPos childHead parentHead childMono parentMono hchildGet hparentGet
    hchildSource hparentSource

theorem pairVecDivVHCConsumeRootBucket_owned_nonroot_get
    (this : DenseUPolyZp) (heap : Array Nat) (k : UInt64)
    (nodes : Array PairVecDivVHCNode) (lin : Array Nat) (resetH : Nat)
    (quotient divisor : SparsePolyZp) (result : PairVecDivVHCBucketResult)
    (owners : Nat → Finset Nat) (hheap : 0 < heap.size)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hrun : pairVecDivVHCConsumeRootBucket this heap k nodes lin resetH
      quotient divisor = .ok result) :
    ∀ i ∈ PairVecDivVHCHeapOwnedNodes heap owners \ owners heap[0],
      result.nodes[i]? = nodes[i]? := by
  intro i hi
  rcases Finset.mem_sdiff.mp hi with ⟨hiOwned, hiRoot⟩
  simp only [PairVecDivVHCHeapOwnedNodes, Finset.mem_biUnion,
    List.mem_toFinset] at hiOwned
  rcases hiOwned with ⟨head, hheadMem, hiHead⟩
  rcases List.getElem_of_mem hheadMem with ⟨slot, hslot, hvalue⟩
  have hslotBound : slot < heap.size := by simpa using hslot
  have harrayValue : heap[slot] = head := by
    rw [← Array.getElem_toList hslotBound]
    exact hvalue
  have hheadGet : heap[slot]? = some head := by
    rw [Array.getElem?_eq_getElem hslotBound, harrayValue]
  have hheadNe : head ≠ heap[0] := by
    intro heq
    apply hiRoot
    rw [← heq]
    exact hiHead
  have hrootOwns := pairVecDivVHCHeapChainOwnership_root_owns heap owners nodes
    hownership hheap
  have hotherOwns := hownership.1 slot head hheadGet
  have hdisjoint : Disjoint (owners head) (owners heap[0]) :=
    hownership.2.2 slot 0 head heap[0] hheadGet
      (Array.getElem?_eq_getElem hheap) hheadNe
  have hpreserved := pairVecDivVHCConsumeRootBucket_preserves_disjoint_chain
    this heap k nodes lin resetH quotient divisor result (owners heap[0])
    (owners head) (some head) hheap hrootOwns hotherOwns hdisjoint hrun
  have hunvisited := pairVecDivVHCConsumeRootBucket_unvisited this heap k nodes
    lin resetH quotient divisor result hheap hrun
  exact hunvisited.2 i (hpreserved.1 hiHead)

theorem pairVecDivVHCConsumeRootExtract_ownedNodesAtDegree
    (this : DenseUPolyZp) (degree : Nat) (heap heap' : Array Nat)
    (k : UInt64) (nodes : Array PairVecDivVHCNode) (lin : Array Nat)
    (resetH : Nat) (quotient divisor : SparsePolyZp)
    (bucket : PairVecDivVHCBucketResult) (owners : Nat → Finset Nat)
    (hheap : 0 < heap.size)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hconsume : pairVecDivVHCConsumeRootBucket this heap k nodes lin resetH
      quotient divisor = .ok bucket)
    (hextract : pairVecDivVHCExtract heap bucket.nodes = .ok heap') :
    PairVecDivVHCHeapOwnedNodesAtDegree degree heap' owners bucket.nodes =
      PairVecDivVHCHeapOwnedNodesAtDegree degree heap owners nodes \
        owners heap[0] := by
  have howned := pairVecDivVHCExtract_heapOwnedNodes heap heap' nodes
    bucket.nodes owners hheap hownership hextract
  have hsame := pairVecDivVHCConsumeRootBucket_owned_nonroot_get this heap k
    nodes lin resetH quotient divisor bucket owners hheap hownership hconsume
  apply Finset.ext
  intro i
  unfold PairVecDivVHCHeapOwnedNodesAtDegree
  rw [howned]
  simp only [Finset.mem_filter, Finset.mem_sdiff,
    pairVecDivVHCNodeAtDegree_iff]
  constructor
  · rintro ⟨⟨hiOwned, hiRoot⟩, mono, hmono, hdegree⟩
    refine ⟨⟨hiOwned, ?_⟩, hiRoot⟩
    rcases (pairVecDivVHCMono_eq_ok_iff i bucket.nodes mono).1 hmono with
      ⟨node, hnode, hactive⟩
    rw [hsame i (Finset.mem_sdiff.mpr ⟨hiOwned, hiRoot⟩)] at hnode
    exact ⟨mono, (pairVecDivVHCMono_eq_ok_iff i nodes mono).2
      ⟨node, hnode, hactive⟩, hdegree⟩
  · rintro ⟨⟨hiOwned, mono, hmono, hdegree⟩, hiRoot⟩
    refine ⟨⟨hiOwned, hiRoot⟩, mono, ?_, hdegree⟩
    rcases (pairVecDivVHCMono_eq_ok_iff i nodes mono).1 hmono with
      ⟨node, hnode, hactive⟩
    rw [← hsame i (Finset.mem_sdiff.mpr ⟨hiOwned, hiRoot⟩)] at hnode
    exact (pairVecDivVHCMono_eq_ok_iff i bucket.nodes mono).2
      ⟨node, hnode, hactive⟩

/-- Consuming the root chain and then executing the generated extract restores
a complete heap order in the updated node array, despite the old root having
become inactive. -/
theorem pairVecDivVHCConsumeRootExtract_preserves_heapOrdered
    (this : DenseUPolyZp) (heap heap' : Array Nat) (k : UInt64)
    (nodes : Array PairVecDivVHCNode) (lin : Array Nat) (resetH : Nat)
    (quotient divisor : SparsePolyZp) (bucket : PairVecDivVHCBucketResult)
    (owners : Nat → Finset Nat) (hheap : 0 < heap.size)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hordered : PairVecDivVHCHeapOrdered heap nodes)
    (hconsume : pairVecDivVHCConsumeRootBucket this heap k nodes lin resetH
      quotient divisor = .ok bucket)
    (hextract : pairVecDivVHCExtract heap bucket.nodes = .ok heap') :
    PairVecDivVHCHeapOrdered heap' bucket.nodes := by
  have hsize := pairVecDivVHCExtract_size heap heap' bucket.nodes hextract
  by_cases hone : heap.size = 1
  · apply PairVecDivVHCHeapDegreesOrderedUpTo.toHeapOrdered
    intro child hchild'
    have : heap'.size = 0 := by omega
    omega
  · have hlimit : heap.size - 1 < heap.size := by omega
    have hlimitPos : 0 < heap.size - 1 := by omega
    have hvalid := hownership.heapPointersValid heap owners nodes
    rcases hvalid (heap.size - 1) hlimit with
      ⟨lastHead, lastNode, lastMono, hlastGet, hlastNode, hlastActive⟩
    have hlastHead : heap[heap.size - 1] = lastHead := by
      rw [Array.getElem?_eq_getElem hlimit] at hlastGet
      exact Option.some.inj hlastGet
    have hlastSource : pairVecDivVHCMono heap[heap.size - 1] nodes =
        .ok lastMono := by
      apply (pairVecDivVHCMono_eq_ok_iff heap[heap.size - 1] nodes lastMono).2
      rw [hlastHead]
      exact ⟨lastNode, hlastNode, hlastActive⟩
    have hlastRun : pairVecDivVHCMono heap[heap.size - 1] bucket.nodes =
        .ok lastMono :=
      (pairVecDivVHCConsumeRootBucket_nonroot_mono_iff this heap k nodes lin
        resetH quotient divisor bucket owners hheap hownership hconsume
        (heap.size - 1) heap[heap.size - 1] hlimitPos
        (Array.getElem?_eq_getElem hlimit) lastMono).1 hlastSource
    unfold pairVecDivVHCExtract at hextract
    simp only [hheap, ↓reduceDIte] at hextract
    cases hsift : pairVecDivVHCSiftDown 0 1 (heap.size - 1)
        heap[heap.size - 1] heap bucket.nodes with
    | error fault => simp [hsift] at hextract
    | ok shifted =>
        rw [hsift] at hextract
        simp only [Except.ok.injEq] at hextract
        subst heap'
        have haway :=
          pairVecDivVHCConsumeRootBucket_preserves_nonroot_heap_order this
            heap k nodes lin resetH quotient divisor bucket owners hheap
            hownership hordered hconsume
        have hsiftOrdered := pairVecDivVHCSiftDown_root_of_nonroot_order
          (heap.size - 1) heap[heap.size - 1] heap shifted bucket.nodes
          lastMono hheap hlimit (Array.getElem?_eq_getElem hlimit) hlastRun
          (by
            intro child hchild hpos hparentPos
            exact haway child (by omega) hpos hparentPos) hsift
        have hshiftSize := pairVecDivVHCSiftDown_size 0 1 (heap.size - 1)
          heap[heap.size - 1] heap shifted bucket.nodes hsift
        have hsiftOrdered' : PairVecDivVHCHeapDegreesOrderedUpTo
            (shifted.size - 1) shifted bucket.nodes := by
          simpa [hshiftSize] using hsiftOrdered
        exact hsiftOrdered'.pop shifted bucket.nodes |>.toHeapOrdered

theorem pairVecDivVHCConsumeRootExtract_root_dominates
    (this : DenseUPolyZp) (heap heap' : Array Nat) (k : UInt64)
    (nodes : Array PairVecDivVHCNode) (lin : Array Nat) (resetH : Nat)
    (quotient divisor : SparsePolyZp) (bucket : PairVecDivVHCBucketResult)
    (owners : Nat → Finset Nat) (hheap : 0 < heap.size)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hordered : PairVecDivVHCHeapOrdered heap nodes)
    (hconsume : pairVecDivVHCConsumeRootBucket this heap k nodes lin resetH
      quotient divisor = .ok bucket)
    (hextract : pairVecDivVHCExtract heap bucket.nodes = .ok heap')
    (hresult : 0 < heap'.size) (slot : Nat) (hslot : slot < heap'.size)
    (mono rootMono : UMonomial)
    (hmono : pairVecDivVHCMono heap'[slot] bucket.nodes = .ok mono)
    (hroot : pairVecDivVHCMono heap'[0] bucket.nodes = .ok rootMono) :
    mono.deg ≤ rootMono.deg := by
  apply pairVecDivVHCExtract_root_dominates heap heap' nodes bucket.nodes
    (hownership.heapPointersValid heap owners nodes) hordered hownership.2.1
    (pairVecDivVHCConsumeRootBucket_nonroot_mono_iff this heap k nodes lin
      resetH quotient divisor bucket owners hheap hownership hconsume)
    hresult hextract slot hslot mono rootMono hmono hroot

/-- After consuming and extracting the root bucket, the surviving heap again
has exact, unique, pairwise-disjoint chain ownership. -/
theorem pairVecDivVHCConsumeRootExtract_preserves_heapChainOwnership
    (this : DenseUPolyZp) (heap heap' : Array Nat) (k : UInt64)
    (nodes : Array PairVecDivVHCNode) (lin : Array Nat) (resetH : Nat)
    (quotient divisor : SparsePolyZp) (bucket : PairVecDivVHCBucketResult)
    (owners : Nat → Finset Nat)
    (hheap : 0 < heap.size)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hconsume : pairVecDivVHCConsumeRootBucket this heap k nodes lin resetH
      quotient divisor = .ok bucket)
    (hextract : pairVecDivVHCExtract heap bucket.nodes = .ok heap') :
    PairVecDivVHCHeapChainOwnership heap' owners bucket.nodes := by
  have hfrom := pairVecDivVHCExtract_valuesFrom heap heap' bucket.nodes hextract
  have hrootOnly := pairVecDivVHCHeapChainOwnership_root_onlyAt heap owners nodes
    hownership hheap
  have hrootGone := pairVecDivVHCExtract_excludes_root heap heap' bucket.nodes
    hheap hrootOnly hextract
  have hsurvives :=
    pairVecDivVHCConsumeRootBucket_preserves_other_heap_chains this heap k nodes
      lin resetH quotient divisor bucket owners hheap hownership hconsume
  refine ⟨?_, ?_, ?_⟩
  · intro slot head hget
    obtain ⟨oldSlot, holdGet⟩ := hfrom slot head hget
    have hne : head ≠ heap[0] := by
      intro heq
      subst head
      exact hrootGone slot hget
    exact hsurvives oldSlot head holdGet hne
  · exact pairVecDivVHCExtract_preserves_unique heap heap' bucket.nodes
      hownership.2.1 hextract
  · intro left right leftHead rightHead hleft hright hne
    obtain ⟨oldLeft, holdLeft⟩ := hfrom left leftHead hleft
    obtain ⟨oldRight, holdRight⟩ := hfrom right rightHead hright
    exact hownership.2.2 oldLeft oldRight leftHead rightHead holdLeft holdRight hne

theorem pairVecDivVHCConsumeRootExtract_preserves_nodesCovered
    (this : DenseUPolyZp) (heap heap' : Array Nat) (k : UInt64)
    (nodes : Array PairVecDivVHCNode) (lin : Array Nat) (resetH : Nat)
    (quotient divisor : SparsePolyZp) (bucket : PairVecDivVHCBucketResult)
    (owners : Nat → Finset Nat) (hheap : 0 < heap.size)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hcovered : PairVecDivVHCNodesCovered heap owners lin resetH nodes)
    (hconsume : pairVecDivVHCConsumeRootBucket this heap k nodes lin resetH
      quotient divisor = .ok bucket)
    (hextract : pairVecDivVHCExtract heap bucket.nodes = .ok heap') :
    PairVecDivVHCNodesCovered heap' owners bucket.lin bucket.resetH
      bucket.nodes := by
  have hrootOwns := pairVecDivVHCHeapChainOwnership_root_owns heap owners nodes
    hownership hheap
  have hownerSubset := pairVecDivVHCChainOwns_subset_range (some heap[0])
    (owners heap[0]) nodes hrootOwns
  have hchainRun : pairVecDivVHCConsumeChain this (some heap[0])
      (Finset.range nodes.size) k nodes lin resetH quotient divisor =
        .ok bucket := by
    unfold pairVecDivVHCConsumeRootBucket at hconsume
    simpa only [hheap, ↓reduceDIte] using hconsume
  have hreclassified := pairVecDivVHCConsumeChain_owner_reclassified this
    (some heap[0]) (owners heap[0]) (Finset.range nodes.size) k nodes lin
    resetH quotient divisor bucket hrootOwns hownerSubset hchainRun
  have hmonotone := pairVecDivVHCConsumeChain_location_monotone this
    (some heap[0]) (Finset.range nodes.size) k nodes lin resetH quotient
    divisor bucket hchainRun
  have hnodesSize := pairVecDivVHCConsumeChain_nodes_size this (some heap[0])
    (Finset.range nodes.size) k nodes lin resetH quotient divisor bucket
    hchainRun
  have hsurvives := pairVecDivVHCExtract_preserves_nonroot_head heap heap'
    bucket.nodes hheap hownership.2.1 hextract
  intro i hi
  have hiOld : i < nodes.size := by rw [hnodesSize] at hi; exact hi
  rcases hcovered i hiOld with hreset | hlin | howned
  · exact Or.inl (Nat.lt_of_lt_of_le hreset hmonotone.2)
  · exact Or.inr (Or.inl (hmonotone.1 hlin))
  · rcases howned with ⟨slot, head, hhead, hiOwner⟩
    by_cases heq : head = heap[0]
    · subst head
      have hclassified := hreclassified hiOwner
      rw [Finset.mem_union] at hclassified
      rcases hclassified with hiLin | hiReset
      · exact Or.inr (Or.inl hiLin)
      · exact Or.inl (by simpa only [Finset.mem_range] using hiReset)
    · rcases hsurvives slot head hhead heq with ⟨newSlot, hnewHead⟩
      exact Or.inr (Or.inr ⟨newSlot, head, hnewHead, hiOwner⟩)

theorem pairVecDivVHCConsumeRootExtract_preserves_heapChainsHomogeneous
    (this : DenseUPolyZp) (heap heap' : Array Nat) (k : UInt64)
    (nodes : Array PairVecDivVHCNode) (lin : Array Nat) (resetH : Nat)
    (quotient divisor : SparsePolyZp) (bucket : PairVecDivVHCBucketResult)
    (owners : Nat → Finset Nat)
    (hheap : 0 < heap.size)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hhomogeneous : PairVecDivVHCHeapChainsHomogeneous heap owners nodes)
    (hconsume : pairVecDivVHCConsumeRootBucket this heap k nodes lin resetH
      quotient divisor = .ok bucket)
    (hextract : pairVecDivVHCExtract heap bucket.nodes = .ok heap') :
    PairVecDivVHCHeapChainsHomogeneous heap' owners bucket.nodes := by
  have hfrom := pairVecDivVHCExtract_valuesFrom heap heap' bucket.nodes hextract
  have hrootGet : heap[0]? = some heap[0] :=
    Array.getElem?_eq_getElem hheap
  have hrootOwns := pairVecDivVHCHeapChainOwnership_root_owns heap owners nodes
    hownership hheap
  have hrootOnly := pairVecDivVHCHeapChainOwnership_root_onlyAt heap owners nodes
    hownership hheap
  have hrootGone := pairVecDivVHCExtract_excludes_root heap heap' bucket.nodes
    hheap hrootOnly hextract
  intro slot head mono hget hmono
  obtain ⟨oldSlot, holdGet⟩ := hfrom slot head hget
  have hne : head ≠ heap[0] := by
    intro heq
    subst head
    exact hrootGone slot hget
  have hotherOwns := hownership.1 oldSlot head holdGet
  have hdisjoint : Disjoint (owners head) (owners heap[0]) :=
    hownership.2.2 oldSlot 0 head heap[0] holdGet hrootGet hne
  have hpreserved :=
    pairVecDivVHCConsumeRootBucket_preserves_disjoint_chain this heap k nodes
      lin resetH quotient divisor bucket (owners heap[0]) (owners head)
      (some head) hheap hrootOwns hotherOwns hdisjoint hconsume
  have hunvisited := pairVecDivVHCConsumeRootBucket_unvisited this heap k nodes
    lin resetH quotient divisor bucket hheap hconsume
  have heqOn : ∀ i ∈ owners head, bucket.nodes[i]? = nodes[i]? := by
    intro i hi
    exact hunvisited.2 i (hpreserved.1 hi)
  rcases pairVecDivVHCMono_eq_ok_iff head bucket.nodes mono |>.mp hmono with
    ⟨headNode, hheadNode, hheadMono⟩
  have hheadMem := pairVecDivVHCChainOwns_head_mem head (owners head) nodes
    hotherOwns
  have hheadNodeOld : nodes[head]? = some headNode := by
    rw [← heqOn head hheadMem]
    exact hheadNode
  have hmonoOld : pairVecDivVHCMono head nodes = .ok mono :=
    pairVecDivVHCMono_eq_ok_iff head nodes mono |>.mpr
      ⟨headNode, hheadNodeOld, hheadMono⟩
  have holdDegree := hhomogeneous oldSlot head mono holdGet hmonoOld
  exact pairVecDivVHCChainAtDegree_congr_on (some head) (owners head) nodes
    bucket.nodes mono.deg holdDegree heqOn

theorem pairVecDivVHCConsumeRootExtract_preserves_away_linReady
    (this : DenseUPolyZp) (heap heap' : Array Nat) (k : UInt64)
    (nodes : Array PairVecDivVHCNode) (lin : Array Nat) (resetH : Nat)
    (quotient divisor : SparsePolyZp) (bucket : PairVecDivVHCBucketResult)
    (owners : Nat → Finset Nat) (hheap : 0 < heap.size)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hseparated : ∀ (slot head : Nat), heap[slot]? = some head →
      Disjoint lin.toList.toFinset (owners head) ∧
        head ∉ lin.toList.toFinset)
    (hlinReady : PairVecDivVHCLinReady lin nodes)
    (hconsume : pairVecDivVHCConsumeRootBucket this heap k nodes lin resetH
      quotient divisor = .ok bucket)
    (hextract : pairVecDivVHCExtract heap bucket.nodes = .ok heap') :
    PairVecDivVHCHeapChainsOwnedAway heap' bucket.nodes
        bucket.lin.toList.toFinset ∧
      PairVecDivVHCLinReady bucket.lin bucket.nodes := by
  have hrootGet : heap[0]? = some heap[0] := by
    rw [Array.getElem?_eq_getElem hheap]
  have hrootOwns := pairVecDivVHCHeapChainOwnership_root_owns heap owners nodes
    hownership hheap
  have hlin := pairVecDivVHCConsumeRootBucket_preserves_linReady this heap k
    nodes lin resetH quotient divisor bucket (owners heap[0]) hheap hrootOwns
    (hseparated 0 heap[0] hrootGet).1 hlinReady hconsume
  have hownership' :=
    pairVecDivVHCConsumeRootExtract_preserves_heapChainOwnership this heap
      heap' k nodes lin resetH quotient divisor bucket owners hheap hownership
      hconsume hextract
  have hfrom := pairVecDivVHCExtract_valuesFrom heap heap' bucket.nodes hextract
  have hrootOnly := pairVecDivVHCHeapChainOwnership_root_onlyAt heap owners nodes
    hownership hheap
  have hrootGone := pairVecDivVHCExtract_excludes_root heap heap' bucket.nodes
    hheap hrootOnly hextract
  refine ⟨⟨owners, hownership', ?_⟩, hlin.1⟩
  intro slot head hget
  obtain ⟨oldSlot, holdGet⟩ := hfrom slot head hget
  have hne : head ≠ heap[0] := by
    intro heq
    subst head
    exact hrootGone slot hget
  have holdSeparated := hseparated oldSlot head holdGet
  have hrootDisjoint : Disjoint (owners heap[0]) (owners head) :=
    hownership.2.2 0 oldSlot heap[0] head hrootGet holdGet hne.symm
  have hheadOwns := hownership.1 oldSlot head holdGet
  have hheadMem := pairVecDivVHCChainOwns_head_mem head (owners head) nodes
    hheadOwns
  refine ⟨Finset.disjoint_left.mpr ?_, ?_⟩
  · intro i hiLin hiOwner
    have hiSource := hlin.2 hiLin
    rw [Finset.mem_union] at hiSource
    rcases hiSource with hiOld | hiRoot
    · exact Finset.disjoint_left.mp holdSeparated.1 hiOld hiOwner
    · exact Finset.disjoint_left.mp hrootDisjoint hiRoot hiOwner
  · intro hheadLin
    have hheadSource := hlin.2 hheadLin
    rw [Finset.mem_union] at hheadSource
    rcases hheadSource with hheadOld | hheadRoot
    · exact holdSeparated.2 hheadOld
    · exact Finset.disjoint_left.mp hrootDisjoint hheadRoot hheadMem

theorem pairVecDivVHCExtractChecked_raw
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode)
    (extracted : { heap' : Array Nat // heap'.size + 1 = heap.size })
    (hrun : pairVecDivVHCExtractChecked heap nodes = .ok extracted) :
    pairVecDivVHCExtract heap nodes = .ok extracted.1 := by
  unfold pairVecDivVHCExtractChecked at hrun
  split at hrun <;> try contradiction
  next heap' hraw =>
    have hval : heap' = extracted.1 := by
      exact congrArg Subtype.val (Except.ok.inj hrun)
    simpa [hval] using hraw

/-- The full equal-degree source loop preserves exact heap-chain ownership.
The induction measure is the real active heap size, strictly decreased by
the checked execution of `VHC_extract`. -/
theorem pairVecDivVHCConsumeEqualDegree_preserves_heapChainOwnership
    (this : DenseUPolyZp) (degree : Nat) (heap : Array Nat) (k : UInt64)
    (nodes : Array PairVecDivVHCNode) (lin : Array Nat) (resetH : Nat)
    (quotient divisor : SparsePolyZp)
    (result : PairVecDivVHCEqualDegreeResult) (owners : Nat → Finset Nat)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hrun : pairVecDivVHCConsumeEqualDegree this degree heap k nodes lin resetH
      quotient divisor = .ok result) :
    PairVecDivVHCHeapChainOwnership result.heap owners result.nodes := by
  induction hsize : heap.size using Nat.strong_induction_on generalizing heap k
      nodes lin resetH result with
  | h size ih =>
      rw [pairVecDivVHCConsumeEqualDegree] at hrun
      by_cases hheap : 0 < heap.size
      · simp only [hheap, ↓reduceDIte] at hrun
        cases hmono : pairVecDivVHCMono heap[0] nodes with
        | error fault => simp [hmono] at hrun
        | ok rootMono =>
            simp only [hmono] at hrun
            by_cases hequal : rootMono.deg = degree
            · simp only [hequal, ↓reduceDIte] at hrun
              cases hconsume : pairVecDivVHCConsumeRootBucket this heap k nodes
                  lin resetH quotient divisor with
              | error fault => simp [hconsume] at hrun
              | ok bucket =>
                  simp only [dif_pos trivial, hconsume] at hrun
                  cases hchecked : pairVecDivVHCExtractChecked heap bucket.nodes with
                  | error fault => simp [hchecked] at hrun
                  | ok extracted =>
                      rw [hchecked] at hrun
                      have hraw := pairVecDivVHCExtractChecked_raw heap
                        bucket.nodes extracted hchecked
                      have hownership' :=
                        pairVecDivVHCConsumeRootExtract_preserves_heapChainOwnership
                          this heap extracted.1 k nodes lin resetH quotient
                          divisor bucket owners hheap hownership hconsume hraw
                      have hsmaller : extracted.1.size < size := by
                        rw [← hsize]
                        omega
                      exact ih extracted.1.size hsmaller extracted.1
                        bucket.coefficient bucket.nodes bucket.lin bucket.resetH
                        result hownership' hrun rfl
            · simp only [hequal, ↓reduceDIte, Except.ok.injEq] at hrun
              subst result
              exact hownership
      · simp only [hheap, ↓reduceDIte, Except.ok.injEq] at hrun
        subst result
        exact hownership

/-- The complete well-founded equal-degree loop carries both exact chain
ownership and the max-heap order needed by every subsequent extraction. -/
theorem pairVecDivVHCConsumeEqualDegree_preserves_heapOwnershipOrdered
    (this : DenseUPolyZp) (degree : Nat) (heap : Array Nat) (k : UInt64)
    (nodes : Array PairVecDivVHCNode) (lin : Array Nat) (resetH : Nat)
    (quotient divisor : SparsePolyZp)
    (result : PairVecDivVHCEqualDegreeResult) (owners : Nat → Finset Nat)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hordered : PairVecDivVHCHeapOrdered heap nodes)
    (hrun : pairVecDivVHCConsumeEqualDegree this degree heap k nodes lin resetH
      quotient divisor = .ok result) :
    PairVecDivVHCHeapChainOwnership result.heap owners result.nodes ∧
      PairVecDivVHCHeapOrdered result.heap result.nodes := by
  induction hsize : heap.size using Nat.strong_induction_on generalizing heap k
      nodes lin resetH result with
  | h size ih =>
      rw [pairVecDivVHCConsumeEqualDegree] at hrun
      by_cases hheap : 0 < heap.size
      · simp only [hheap, ↓reduceDIte] at hrun
        cases hmono : pairVecDivVHCMono heap[0] nodes with
        | error fault => simp [hmono] at hrun
        | ok rootMono =>
            simp only [hmono] at hrun
            by_cases hequal : rootMono.deg = degree
            · simp only [hequal, ↓reduceDIte] at hrun
              cases hconsume : pairVecDivVHCConsumeRootBucket this heap k nodes
                  lin resetH quotient divisor with
              | error fault => simp [hconsume] at hrun
              | ok bucket =>
                  simp only [dif_pos trivial, hconsume] at hrun
                  cases hchecked : pairVecDivVHCExtractChecked heap bucket.nodes with
                  | error fault => simp [hchecked] at hrun
                  | ok extracted =>
                      rw [hchecked] at hrun
                      have hraw := pairVecDivVHCExtractChecked_raw heap
                        bucket.nodes extracted hchecked
                      have hownership' :=
                        pairVecDivVHCConsumeRootExtract_preserves_heapChainOwnership
                          this heap extracted.1 k nodes lin resetH quotient
                          divisor bucket owners hheap hownership hconsume hraw
                      have hordered' :=
                        pairVecDivVHCConsumeRootExtract_preserves_heapOrdered
                          this heap extracted.1 k nodes lin resetH quotient
                          divisor bucket owners hheap hownership hordered
                          hconsume hraw
                      have hsmaller : extracted.1.size < size := by
                        rw [← hsize]
                        omega
                      exact ih extracted.1.size hsmaller extracted.1
                        bucket.coefficient bucket.nodes bucket.lin bucket.resetH
                        result hownership' hordered' hrun rfl
            · simp only [hequal, ↓reduceDIte, Except.ok.injEq] at hrun
              subst result
              exact ⟨hownership, hordered⟩
      · simp only [hheap, ↓reduceDIte, Except.ok.injEq] at hrun
        subst result
        exact ⟨hownership, hordered⟩

theorem pairVecDivVHCConsumeEqualDegree_preserves_heapChainsHomogeneous
    (this : DenseUPolyZp) (degree : Nat) (heap : Array Nat) (k : UInt64)
    (nodes : Array PairVecDivVHCNode) (lin : Array Nat) (resetH : Nat)
    (quotient divisor : SparsePolyZp)
    (result : PairVecDivVHCEqualDegreeResult) (owners : Nat → Finset Nat)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hhomogeneous : PairVecDivVHCHeapChainsHomogeneous heap owners nodes)
    (hrun : pairVecDivVHCConsumeEqualDegree this degree heap k nodes lin resetH
      quotient divisor = .ok result) :
    PairVecDivVHCHeapChainOwnership result.heap owners result.nodes ∧
      PairVecDivVHCHeapChainsHomogeneous result.heap owners result.nodes := by
  induction hsize : heap.size using Nat.strong_induction_on generalizing heap k
      nodes lin resetH result with
  | h size ih =>
      rw [pairVecDivVHCConsumeEqualDegree] at hrun
      by_cases hheap : 0 < heap.size
      · simp only [hheap, ↓reduceDIte] at hrun
        cases hmono : pairVecDivVHCMono heap[0] nodes with
        | error fault => simp [hmono] at hrun
        | ok rootMono =>
            simp only [hmono] at hrun
            by_cases hequal : rootMono.deg = degree
            · simp only [hequal, ↓reduceDIte] at hrun
              cases hconsume : pairVecDivVHCConsumeRootBucket this heap k nodes
                  lin resetH quotient divisor with
              | error fault => simp [hconsume] at hrun
              | ok bucket =>
                  simp only [dif_pos trivial, hconsume] at hrun
                  cases hchecked : pairVecDivVHCExtractChecked heap bucket.nodes with
                  | error fault => simp [hchecked] at hrun
                  | ok extracted =>
                      rw [hchecked] at hrun
                      have hraw := pairVecDivVHCExtractChecked_raw heap
                        bucket.nodes extracted hchecked
                      have hownership' :=
                        pairVecDivVHCConsumeRootExtract_preserves_heapChainOwnership
                          this heap extracted.1 k nodes lin resetH quotient
                          divisor bucket owners hheap hownership hconsume hraw
                      have hhomogeneous' :=
                        pairVecDivVHCConsumeRootExtract_preserves_heapChainsHomogeneous
                          this heap extracted.1 k nodes lin resetH quotient
                          divisor bucket owners hheap hownership hhomogeneous
                          hconsume hraw
                      have hsmaller : extracted.1.size < size := by
                        rw [← hsize]
                        omega
                      exact ih extracted.1.size hsmaller extracted.1
                        bucket.coefficient bucket.nodes bucket.lin bucket.resetH
                        result hownership' hhomogeneous' hrun rfl
            · simp only [hequal, ↓reduceDIte, Except.ok.injEq] at hrun
              subst result
              exact ⟨hownership, hhomogeneous⟩
      · simp only [hheap, ↓reduceDIte, Except.ok.injEq] at hrun
        subst result
        exact ⟨hownership, hhomogeneous⟩

theorem pairVecDivVHCConsumeEqualDegree_preserves_cursorPrefixAbove
    (this : DenseUPolyZp) (degree : Nat) (heap : Array Nat) (k : UInt64)
    (nodes : Array PairVecDivVHCNode) (lin : Array Nat) (resetH : Nat)
    (quotient divisor : SparsePolyZp)
    (result : PairVecDivVHCEqualDegreeResult) (owners : Nat → Finset Nat)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hhomogeneous : PairVecDivVHCHeapChainsHomogeneous heap owners nodes)
    (hprefix : PairVecDivVHCCursorPrefixAbove degree nodes quotient divisor)
    (hdenotes : ∀ (i : Nat) (node : PairVecDivVHCNode),
      nodes[i]? = some node → node.mono ≠ none →
        PairVecDivVHCNodeDenotes quotient divisor node)
    (hrun : pairVecDivVHCConsumeEqualDegree this degree heap k nodes lin resetH
      quotient divisor = .ok result) :
    PairVecDivVHCCursorPrefixAbove degree result.nodes quotient divisor := by
  induction hsize : heap.size using Nat.strong_induction_on generalizing heap k
      nodes lin resetH result with
  | h size ih =>
      rw [pairVecDivVHCConsumeEqualDegree] at hrun
      by_cases hheap : 0 < heap.size
      · simp only [hheap, ↓reduceDIte] at hrun
        cases hmono : pairVecDivVHCMono heap[0] nodes with
        | error fault => simp [hmono] at hrun
        | ok rootMono =>
            simp only [hmono] at hrun
            by_cases hequal : rootMono.deg = degree
            · simp only [hequal, ↓reduceDIte] at hrun
              cases hconsume : pairVecDivVHCConsumeRootBucket this heap k nodes
                  lin resetH quotient divisor with
              | error fault => simp [hconsume] at hrun
              | ok bucket =>
                  simp only [dif_pos trivial, hconsume] at hrun
                  cases hchecked : pairVecDivVHCExtractChecked heap bucket.nodes with
                  | error fault => simp [hchecked] at hrun
                  | ok extracted =>
                      rw [hchecked] at hrun
                      have hrootOwns :=
                        pairVecDivVHCHeapChainOwnership_root_owns heap owners
                          nodes hownership hheap
                      have hrootDegree :=
                        hhomogeneous 0 heap[0] rootMono
                          (by rw [Array.getElem?_eq_getElem hheap]) hmono
                      have hdegree : PairVecDivVHCChainAtDegree
                          (some heap[0]) (Finset.range nodes.size) nodes degree := by
                        rw [← hequal]
                        exact pairVecDivVHCChainAtDegree_mono
                          (some heap[0]) (owners heap[0])
                          (Finset.range nodes.size) nodes rootMono.deg hrootDegree
                          (pairVecDivVHCChainOwns_subset_range
                            (some heap[0]) (owners heap[0]) nodes hrootOwns)
                      have hprefix' :=
                        pairVecDivVHCConsumeRootBucket_preserves_cursorPrefixAbove
                          this degree heap k nodes lin resetH quotient divisor
                          bucket hheap hdegree hprefix hdenotes hconsume
                      have hdenotes' :=
                        pairVecDivVHCConsumeRootBucket_preserves_denotes this heap
                          k nodes lin resetH quotient divisor bucket hdenotes
                          hconsume
                      have hraw := pairVecDivVHCExtractChecked_raw heap
                        bucket.nodes extracted hchecked
                      have hownership' :=
                        pairVecDivVHCConsumeRootExtract_preserves_heapChainOwnership
                          this heap extracted.1 k nodes lin resetH quotient
                          divisor bucket owners hheap hownership hconsume hraw
                      have hhomogeneous' :=
                        pairVecDivVHCConsumeRootExtract_preserves_heapChainsHomogeneous
                          this heap extracted.1 k nodes lin resetH quotient
                          divisor bucket owners hheap hownership hhomogeneous
                          hconsume hraw
                      have hsmaller : extracted.1.size < size := by
                        rw [← hsize]
                        omega
                      exact ih extracted.1.size hsmaller extracted.1
                        bucket.coefficient bucket.nodes bucket.lin bucket.resetH
                        result hownership' hhomogeneous' hprefix' hdenotes' hrun rfl
            · simp only [hequal, ↓reduceDIte, Except.ok.injEq] at hrun
              subst result
              exact hprefix
      · simp only [hheap, ↓reduceDIte, Except.ok.injEq] at hrun
        subst result
        exact hprefix

theorem pairVecDivVHCConsumeEqualDegree_preserves_nodesCovered
    (this : DenseUPolyZp) (degree : Nat) (heap : Array Nat) (k : UInt64)
    (nodes : Array PairVecDivVHCNode) (lin : Array Nat) (resetH : Nat)
    (quotient divisor : SparsePolyZp)
    (result : PairVecDivVHCEqualDegreeResult) (owners : Nat → Finset Nat)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hcovered : PairVecDivVHCNodesCovered heap owners lin resetH nodes)
    (hrun : pairVecDivVHCConsumeEqualDegree this degree heap k nodes lin resetH
      quotient divisor = .ok result) :
    PairVecDivVHCHeapChainOwnership result.heap owners result.nodes ∧
      PairVecDivVHCNodesCovered result.heap owners result.lin result.resetH
        result.nodes := by
  induction hsize : heap.size using Nat.strong_induction_on generalizing heap k
      nodes lin resetH result with
  | h size ih =>
      rw [pairVecDivVHCConsumeEqualDegree] at hrun
      by_cases hheap : 0 < heap.size
      · simp only [hheap, ↓reduceDIte] at hrun
        cases hmono : pairVecDivVHCMono heap[0] nodes with
        | error fault => simp [hmono] at hrun
        | ok rootMono =>
            simp only [hmono] at hrun
            by_cases hequal : rootMono.deg = degree
            · simp only [hequal, ↓reduceDIte] at hrun
              cases hconsume : pairVecDivVHCConsumeRootBucket this heap k nodes
                  lin resetH quotient divisor with
              | error fault => simp [hconsume] at hrun
              | ok bucket =>
                  simp only [dif_pos trivial, hconsume] at hrun
                  cases hchecked : pairVecDivVHCExtractChecked heap bucket.nodes with
                  | error fault => simp [hchecked] at hrun
                  | ok extracted =>
                      rw [hchecked] at hrun
                      have hraw := pairVecDivVHCExtractChecked_raw heap
                        bucket.nodes extracted hchecked
                      have hownership' :=
                        pairVecDivVHCConsumeRootExtract_preserves_heapChainOwnership
                          this heap extracted.1 k nodes lin resetH quotient
                          divisor bucket owners hheap hownership hconsume hraw
                      have hcovered' :=
                        pairVecDivVHCConsumeRootExtract_preserves_nodesCovered
                          this heap extracted.1 k nodes lin resetH quotient
                          divisor bucket owners hheap hownership hcovered
                          hconsume hraw
                      have hsmaller : extracted.1.size < size := by
                        rw [← hsize]
                        omega
                      exact ih extracted.1.size hsmaller extracted.1
                        bucket.coefficient bucket.nodes bucket.lin bucket.resetH
                        result hownership' hcovered' hrun rfl
            · simp only [hequal, ↓reduceDIte, Except.ok.injEq] at hrun
              subst result
              exact ⟨hownership, hcovered⟩
      · simp only [hheap, ↓reduceDIte, Except.ok.injEq] at hrun
        subst result
        exact ⟨hownership, hcovered⟩

theorem pairVecDivVHCConsumeEqualDegree_preserves_stateCovered
    (this : DenseUPolyZp) (degree : Nat) (heap : Array Nat) (k : UInt64)
    (nodes : Array PairVecDivVHCNode) (lin : Array Nat) (resetH : Nat)
    (quotient divisor : SparsePolyZp)
    (result : PairVecDivVHCEqualDegreeResult)
    (hcovered : PairVecDivVHCStateCovered heap nodes lin resetH)
    (hrun : pairVecDivVHCConsumeEqualDegree this degree heap k nodes lin resetH
      quotient divisor = .ok result) :
    PairVecDivVHCStateCovered result.heap result.nodes result.lin
      result.resetH := by
  rcases hcovered with ⟨owners, hownership, hnodesCovered⟩
  rcases pairVecDivVHCConsumeEqualDegree_preserves_nodesCovered this degree
      heap k nodes lin resetH quotient divisor result owners hownership
      hnodesCovered hrun with ⟨hownership', hnodesCovered'⟩
  exact ⟨owners, hownership', hnodesCovered'⟩

theorem pairVecDivVHCConsumeEqualDegree_coefficient_semantics
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (degree : Nat) (heap : Array Nat) (k : UInt64)
    (nodes : Array PairVecDivVHCNode) (lin : Array Nat) (resetH : Nat)
    (quotient divisor : SparsePolyZp)
    (result : PairVecDivVHCEqualDegreeResult) (owners : Nat → Finset Nat)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hhomogeneous : PairVecDivVHCHeapChainsHomogeneous heap owners nodes)
    (hdenotes : ∀ (i : Nat) (node : PairVecDivVHCNode),
      nodes[i]? = some node → node.mono ≠ none →
        PairVecDivVHCNodeDenotes quotient divisor node)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hcanonical : SparsePolyZp.Canonical this._p.toNat quotient)
    (hk : k.toNat < this._p.toNat)
    (hrun : pairVecDivVHCConsumeEqualDegree this degree heap k nodes lin resetH
      quotient divisor = .ok result) :
    ∃ products : List (UInt64 × UInt64),
      (result.coefficient.toNat : ZMod this._p.toNat) =
          (k.toNat : ZMod this._p.toNat) -
            pairVecDivVHCProductsValue this._p.toNat products ∧
        ∀ product ∈ products,
          PairVecDivVHCStoredProductAtDegree degree quotient divisor
            product := by
  induction hsize : heap.size using Nat.strong_induction_on generalizing heap k
      nodes lin resetH result with
  | h size ih =>
      rw [pairVecDivVHCConsumeEqualDegree] at hrun
      by_cases hheap : 0 < heap.size
      · simp only [hheap, ↓reduceDIte] at hrun
        cases hmono : pairVecDivVHCMono heap[0] nodes with
        | error fault => simp [hmono] at hrun
        | ok rootMono =>
            simp only [hmono] at hrun
            by_cases hequal : rootMono.deg = degree
            · simp only [hequal, ↓reduceDIte] at hrun
              cases hconsume : pairVecDivVHCConsumeRootBucket this heap k nodes
                  lin resetH quotient divisor with
              | error fault => simp [hconsume] at hrun
              | ok bucket =>
                  simp only [dif_pos trivial, hconsume] at hrun
                  cases hchecked : pairVecDivVHCExtractChecked heap bucket.nodes with
                  | error fault => simp [hchecked] at hrun
                  | ok extracted =>
                      rw [hchecked] at hrun
                      have hrootGet : heap[0]? = some heap[0] :=
                        Array.getElem?_eq_getElem hheap
                      have hrootOwns :=
                        pairVecDivVHCHeapChainOwnership_root_owns heap owners
                          nodes hownership hheap
                      have hrootDegree := hhomogeneous 0 heap[0] rootMono
                        hrootGet hmono
                      rw [hequal] at hrootDegree
                      rcases
                          pairVecDivVHCConsumeRootBucket_coefficient_semantics_of_chainDegree
                            this degree heap k nodes lin resetH quotient divisor
                            bucket (owners heap[0]) hheap hrootOwns hrootDegree
                            hdenotes hcfg hcanonical hk hconsume with
                        ⟨rootProducts, hrootTrace, hrootCoefficient,
                          hrootProductsDegree⟩
                      have hraw := pairVecDivVHCExtractChecked_raw heap
                        bucket.nodes extracted hchecked
                      have hownership' :=
                        pairVecDivVHCConsumeRootExtract_preserves_heapChainOwnership
                          this heap extracted.1 k nodes lin resetH quotient
                          divisor bucket owners hheap hownership hconsume hraw
                      have hhomogeneous' :=
                        pairVecDivVHCConsumeRootExtract_preserves_heapChainsHomogeneous
                          this heap extracted.1 k nodes lin resetH quotient
                          divisor bucket owners hheap hownership hhomogeneous
                          hconsume hraw
                      have hdenotes' :=
                        pairVecDivVHCConsumeRootBucket_preserves_denotes this
                          heap k nodes lin resetH quotient divisor bucket
                          hdenotes hconsume
                      have hp : this._p ≠ 0 := by
                        intro hp
                        have hzero : this._p.toNat = 0 := congrArg UInt64.toNat hp
                        exact (Fact.out : Nat.Prime this._p.toNat).ne_zero hzero
                      have hk' :=
                        pairVecDivVHCConsumeRootBucket_coefficient_reduced this
                          heap k nodes lin resetH quotient divisor bucket hp
                          hcfg hcanonical hk hconsume
                      have hsmaller : extracted.1.size < size := by
                        rw [← hsize]
                        omega
                      rcases ih extracted.1.size hsmaller extracted.1
                          bucket.coefficient bucket.nodes bucket.lin
                          bucket.resetH result hownership' hhomogeneous'
                          hdenotes' hk' hrun rfl with
                        ⟨tailProducts, htailCoefficient,
                          htailProductsDegree⟩
                      refine ⟨rootProducts ++ tailProducts, ?_, ?_⟩
                      · rw [pairVecDivVHCProductsValue_append,
                          htailCoefficient, hrootCoefficient]
                        ring
                      · intro product hproduct
                        rw [List.mem_append] at hproduct
                        exact hproduct.elim
                          (hrootProductsDegree product)
                          (htailProductsDegree product)
            · simp only [hequal, ↓reduceDIte, Except.ok.injEq] at hrun
              subst result
              refine ⟨[], ?_, by simp⟩
              simp [pairVecDivVHCProductsValue]
      · simp only [hheap, ↓reduceDIte, Except.ok.injEq] at hrun
        subst result
        refine ⟨[], ?_, by simp⟩
        simp [pairVecDivVHCProductsValue]

/-- Full equal-degree coefficient semantics: the existential trace is both
sound and complete for every initially heap-owned chain at the target degree. -/
theorem pairVecDivVHCConsumeEqualDegree_products_complete
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (degree : Nat) (heap : Array Nat) (k : UInt64)
    (nodes : Array PairVecDivVHCNode) (lin : Array Nat) (resetH : Nat)
    (quotient divisor : SparsePolyZp)
    (result : PairVecDivVHCEqualDegreeResult) (owners : Nat → Finset Nat)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hhomogeneous : PairVecDivVHCHeapChainsHomogeneous heap owners nodes)
    (hordered : PairVecDivVHCHeapOrdered heap nodes)
    (hmax : ∀ (slot head : Nat) (mono : UMonomial), heap[slot]? = some head →
      pairVecDivVHCMono head nodes = .ok mono → mono.deg ≤ degree)
    (hdenotes : ∀ (i : Nat) (node : PairVecDivVHCNode),
      nodes[i]? = some node → node.mono ≠ none →
        PairVecDivVHCNodeDenotes quotient divisor node)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hcanonical : SparsePolyZp.Canonical this._p.toNat quotient)
    (hk : k.toNat < this._p.toNat)
    (hrun : pairVecDivVHCConsumeEqualDegree this degree heap k nodes lin resetH
      quotient divisor = .ok result) :
    ∃ products : List (UInt64 × UInt64),
      (result.coefficient.toNat : ZMod this._p.toNat) =
          (k.toNat : ZMod this._p.toNat) -
            pairVecDivVHCProductsValue this._p.toNat products ∧
      (∀ product ∈ products,
        PairVecDivVHCStoredProductAtDegree degree quotient divisor product) ∧
      PairVecDivVHCProductsCoverDegreeOwners degree heap owners nodes quotient
        divisor products ∧
      pairVecDivVHCProductsValue this._p.toNat products =
        ∑ i ∈ PairVecDivVHCHeapOwnedNodesAtDegree degree heap owners nodes,
          pairVecDivVHCNodeProductValue this._p.toNat nodes quotient divisor i := by
  induction hsize : heap.size using Nat.strong_induction_on generalizing heap k
      nodes lin resetH result with
  | h size ih =>
      rw [pairVecDivVHCConsumeEqualDegree] at hrun
      by_cases hheap : 0 < heap.size
      · simp only [hheap, ↓reduceDIte] at hrun
        cases hmono : pairVecDivVHCMono heap[0] nodes with
        | error fault => simp [hmono] at hrun
        | ok rootMono =>
            simp only [hmono] at hrun
            by_cases hequal : rootMono.deg = degree
            · simp only [hequal, ↓reduceDIte] at hrun
              cases hconsume : pairVecDivVHCConsumeRootBucket this heap k nodes
                  lin resetH quotient divisor with
              | error fault => simp [hconsume] at hrun
              | ok bucket =>
                  simp only [dif_pos trivial, hconsume] at hrun
                  cases hchecked : pairVecDivVHCExtractChecked heap bucket.nodes with
                  | error fault => simp [hchecked] at hrun
                  | ok extracted =>
                      rw [hchecked] at hrun
                      have hrootOwns :=
                        pairVecDivVHCHeapChainOwnership_root_owns heap owners
                          nodes hownership hheap
                      have hrootDegree := hhomogeneous 0 heap[0] rootMono
                        (Array.getElem?_eq_getElem hheap) hmono
                      rw [hequal] at hrootDegree
                      rcases pairVecDivVHCConsumeRootBucket_products_complete
                          this degree heap k nodes lin resetH quotient divisor
                          bucket (owners heap[0]) hheap hrootOwns hrootDegree
                          hdenotes hcfg hcanonical hk hconsume with
                        ⟨rootProducts, hrootTrace, hrootCoefficient, hrootSound,
                          hrootCover⟩
                      have hraw := pairVecDivVHCExtractChecked_raw heap
                        bucket.nodes extracted hchecked
                      have hownership' :=
                        pairVecDivVHCConsumeRootExtract_preserves_heapChainOwnership
                          this heap extracted.1 k nodes lin resetH quotient
                          divisor bucket owners hheap hownership hconsume hraw
                      have hhomogeneous' :=
                        pairVecDivVHCConsumeRootExtract_preserves_heapChainsHomogeneous
                          this heap extracted.1 k nodes lin resetH quotient
                          divisor bucket owners hheap hownership hhomogeneous
                          hconsume hraw
                      have hordered' :=
                        pairVecDivVHCConsumeRootExtract_preserves_heapOrdered
                          this heap extracted.1 k nodes lin resetH quotient
                          divisor bucket owners hheap hownership hordered
                          hconsume hraw
                      have hdenotes' :=
                        pairVecDivVHCConsumeRootBucket_preserves_denotes this
                          heap k nodes lin resetH quotient divisor bucket
                          hdenotes hconsume
                      have hp : this._p ≠ 0 := by
                        intro hp
                        have hz : this._p.toNat = 0 := congrArg UInt64.toNat hp
                        exact (Fact.out : Nat.Prime this._p.toNat).ne_zero hz
                      have hk' :=
                        pairVecDivVHCConsumeRootBucket_coefficient_reduced this
                          heap k nodes lin resetH quotient divisor bucket hp
                          hcfg hcanonical hk hconsume
                      have hsurvives :=
                        pairVecDivVHCExtract_preserves_nonroot_head heap
                          extracted.1 bucket.nodes hheap hownership.2.1 hraw
                      have hsameOther :=
                        pairVecDivVHCConsumeRootBucket_preserves_disjoint_chainAtDegree
                      have hmax' : ∀ (slot head : Nat) (mono : UMonomial),
                          extracted.1[slot]? = some head →
                          pairVecDivVHCMono head bucket.nodes = .ok mono →
                          mono.deg ≤ degree := by
                        intro slot head mono hget hrunMono
                        rcases pairVecDivVHCExtract_valuesFrom heap extracted.1
                            bucket.nodes hraw slot head hget with
                          ⟨oldSlot, hold⟩
                        have hslotPos : 0 < oldSlot := by
                          by_contra hz
                          have : oldSlot = 0 := by omega
                          subst oldSlot
                          have hgone := pairVecDivVHCExtract_excludes_root heap
                            extracted.1 bucket.nodes hheap
                            (pairVecDivVHCHeapChainOwnership_root_onlyAt heap
                              owners nodes hownership hheap) hraw slot
                          have hheadRoot : head = heap[0] := by
                            rw [Array.getElem?_eq_getElem hheap] at hold
                            exact (Option.some.inj hold).symm
                          exact hgone (by simpa [hheadRoot] using hget)
                        have hsource :=
                          (pairVecDivVHCConsumeRootBucket_nonroot_mono_iff this
                            heap k nodes lin resetH quotient divisor bucket
                            owners hheap hownership hconsume oldSlot head
                            hslotPos hold mono).mpr hrunMono
                        exact hmax oldSlot head mono hold hsource
                      have hsmaller : extracted.1.size < size := by
                        rw [← hsize]
                        omega
                      rcases ih extracted.1.size hsmaller extracted.1
                          bucket.coefficient bucket.nodes bucket.lin
                          bucket.resetH result hownership' hhomogeneous'
                          hordered' hmax' hdenotes' hk' hrun rfl with
                        ⟨tailProducts, htailCoefficient, htailSound,
                          htailCover, htailValue⟩
                      have hrootValue :=
                        hrootTrace.productsValue_eq_owner_sum this
                          this._p.toNat quotient divisor (some heap[0])
                          (Finset.range nodes.size) (owners heap[0]) k nodes lin
                          resetH bucket rootProducts hrootOwns
                      have htargetStep :=
                        pairVecDivVHCConsumeRootExtract_ownedNodesAtDegree this
                          degree heap extracted.1 k nodes lin resetH quotient
                          divisor bucket owners hheap hownership hconsume hraw
                      have hsame :=
                        pairVecDivVHCConsumeRootBucket_owned_nonroot_get this
                          heap k nodes lin resetH quotient divisor bucket owners
                          hheap hownership hconsume
                      have htailSame :
                          (∑ i ∈ PairVecDivVHCHeapOwnedNodesAtDegree degree
                              heap owners nodes \ owners heap[0],
                            pairVecDivVHCNodeProductValue this._p.toNat
                              bucket.nodes quotient divisor i) =
                          ∑ i ∈ PairVecDivVHCHeapOwnedNodesAtDegree degree
                              heap owners nodes \ owners heap[0],
                            pairVecDivVHCNodeProductValue this._p.toNat nodes
                              quotient divisor i := by
                        apply Finset.sum_congr rfl
                        intro i hi
                        unfold pairVecDivVHCNodeProductValue
                        rw [hsame i (Finset.mem_sdiff.mpr ⟨
                          (Finset.mem_filter.mp (Finset.mem_sdiff.mp hi).1).1,
                          (Finset.mem_sdiff.mp hi).2⟩)]
                      have hrootSubset :=
                        pairVecDivVHCRootOwner_subset_ownedNodesAtDegree degree
                          heap owners nodes rootMono hheap hownership
                          hhomogeneous hmono hequal
                      have hsplit := Finset.sum_sdiff
                        (f := pairVecDivVHCNodeProductValue this._p.toNat nodes
                          quotient divisor) hrootSubset
                      refine ⟨rootProducts ++ tailProducts, ?_, ?_, ?_, ?_⟩
                      · rw [pairVecDivVHCProductsValue_append,
                          htailCoefficient, hrootCoefficient]
                        ring
                      · intro product hproduct
                        rw [List.mem_append] at hproduct
                        exact hproduct.elim (hrootSound product)
                          (htailSound product)
                      · intro slot head hhead hdegree i hi
                        by_cases hrootHead : head = heap[0]
                        · subst head
                          rcases hrootCover i hi with
                            ⟨node, quotientTerm, divisorTerm, hnode,
                              hquotient, hdivisor, hproduct⟩
                          exact ⟨node, quotientTerm, divisorTerm, hnode,
                            hquotient, hdivisor,
                            List.mem_append_left _ hproduct⟩
                        · rcases hsurvives slot head hhead hrootHead with
                            ⟨newSlot, hnewHead⟩
                          have hrootOwns' := hrootOwns
                          have hotherOwns := hownership.1 slot head hhead
                          have hdisjoint : Disjoint (owners head)
                              (owners heap[0]) := by
                            have hslot : slot < heap.size := by
                              by_contra hn
                              rw [Array.getElem?_eq_none (by omega)] at hhead
                              contradiction
                            have hheadEq : heap[slot] = head := by
                              rw [Array.getElem?_eq_getElem hslot] at hhead
                              exact Option.some.inj hhead
                            exact hownership.2.2 slot 0 head heap[0] hhead
                              (Array.getElem?_eq_getElem hheap) hrootHead
                          have hpreserved := hsameOther this heap k nodes lin
                            resetH quotient divisor bucket (owners heap[0])
                            (owners head) (some head) degree hheap hrootOwns'
                            hotherOwns hdegree hdisjoint hconsume
                          rcases htailCover newSlot head hnewHead hpreserved.2
                              i hi with
                            ⟨node, quotientTerm, divisorTerm, hnode,
                              hquotient, hdivisor, hproduct⟩
                          rw [hpreserved.1 i hi] at hnode
                          exact ⟨node, quotientTerm, divisorTerm, hnode,
                            hquotient, hdivisor,
                            List.mem_append_right _ hproduct⟩
                      · rw [pairVecDivVHCProductsValue_append, hrootValue,
                          htailValue, htargetStep, htailSame]
                        rw [add_comm]
                        exact hsplit
            · simp only [hequal, ↓reduceDIte, Except.ok.injEq] at hrun
              subst result
              refine ⟨[], ?_, by simp, ?_, ?_⟩
              · simp [pairVecDivVHCProductsValue]
              · intro slot head hhead hdegree
                rw [PairVecDivVHCChainAtDegree] at hdegree
                split at hdegree <;> try contradiction
                next hmem =>
                  rcases hdegree with
                    ⟨headNode, headMono, hheadNode, hheadMono,
                      hheadDegree, htail⟩
                  have hheadRun : pairVecDivVHCMono head nodes =
                      .ok headMono :=
                    (pairVecDivVHCMono_eq_ok_iff head nodes headMono).2
                      ⟨headNode, hheadNode, hheadMono⟩
                  have hslot : slot < heap.size := by
                    by_contra hn
                    rw [Array.getElem?_eq_none (by omega)] at hhead
                    contradiction
                  have hheadEq : heap[slot] = head := by
                    rw [Array.getElem?_eq_getElem hslot] at hhead
                    exact Option.some.inj hhead
                  have hslotRun : pairVecDivVHCMono heap[slot] nodes =
                      .ok headMono := by simpa [hheadEq] using hheadRun
                  have hleRoot := pairVecDivVHCHeapOrdered_slot_le_root heap
                    nodes (hownership.heapPointersValid heap owners nodes)
                    hordered slot hslot headMono rootMono hslotRun hmono
                  have hrootLe := hmax 0 heap[0] rootMono
                    (Array.getElem?_eq_getElem hheap) hmono
                  have : rootMono.deg = degree := by omega
                  exact (hequal this).elim
              · have hempty :=
                  pairVecDivVHCHeapOwnedNodesAtDegree_eq_empty_of_root_ne
                    degree heap owners nodes rootMono hheap hownership
                    hhomogeneous hordered hmono hmax hequal
                simp [pairVecDivVHCProductsValue, hempty]
      · simp only [hheap, ↓reduceDIte, Except.ok.injEq] at hrun
        subst result
        refine ⟨[], ?_, by simp, ?_, ?_⟩
        · simp [pairVecDivVHCProductsValue]
        · intro slot head hget
          rw [Array.getElem?_eq_none (by omega)] at hget
          contradiction
        · have hempty : heap.size = 0 := by omega
          have heq : heap = #[] := Array.eq_empty_of_size_eq_zero hempty
          subst heap
          simp [pairVecDivVHCProductsValue,
            PairVecDivVHCHeapOwnedNodesAtDegree,
            PairVecDivVHCHeapOwnedNodes]

theorem pairVecDivVHCConsumeEqualDegree_preserves_heapChainsOwned
    (this : DenseUPolyZp) (degree : Nat) (heap : Array Nat) (k : UInt64)
    (nodes : Array PairVecDivVHCNode) (lin : Array Nat) (resetH : Nat)
    (quotient divisor : SparsePolyZp)
    (result : PairVecDivVHCEqualDegreeResult)
    (howned : PairVecDivVHCHeapChainsOwned heap nodes)
    (hrun : pairVecDivVHCConsumeEqualDegree this degree heap k nodes lin resetH
      quotient divisor = .ok result) :
    PairVecDivVHCHeapChainsOwned result.heap result.nodes := by
  rcases howned with ⟨owners, hownership⟩
  exact ⟨owners,
    pairVecDivVHCConsumeEqualDegree_preserves_heapChainOwnership this degree
      heap k nodes lin resetH quotient divisor result owners hownership hrun⟩

theorem pairVecDivVHCConsumeEqualDegree_preserves_away_linReady
    (this : DenseUPolyZp) (degree : Nat) (heap : Array Nat) (k : UInt64)
    (nodes : Array PairVecDivVHCNode) (lin : Array Nat) (resetH : Nat)
    (quotient divisor : SparsePolyZp)
    (result : PairVecDivVHCEqualDegreeResult)
    (haway : PairVecDivVHCHeapChainsOwnedAway heap nodes
      lin.toList.toFinset)
    (hlinReady : PairVecDivVHCLinReady lin nodes)
    (hrun : pairVecDivVHCConsumeEqualDegree this degree heap k nodes lin resetH
      quotient divisor = .ok result) :
    PairVecDivVHCHeapChainsOwnedAway result.heap result.nodes
        result.lin.toList.toFinset ∧
      PairVecDivVHCLinReady result.lin result.nodes := by
  induction hsize : heap.size using Nat.strong_induction_on generalizing heap k
      nodes lin resetH result with
  | h size ih =>
      rw [pairVecDivVHCConsumeEqualDegree] at hrun
      by_cases hheap : 0 < heap.size
      · simp only [hheap, ↓reduceDIte] at hrun
        cases hmono : pairVecDivVHCMono heap[0] nodes with
        | error fault => simp [hmono] at hrun
        | ok rootMono =>
            simp only [hmono] at hrun
            by_cases hequal : rootMono.deg = degree
            · simp only [hequal, ↓reduceDIte] at hrun
              cases hconsume : pairVecDivVHCConsumeRootBucket this heap k nodes
                  lin resetH quotient divisor with
              | error fault => simp [hconsume] at hrun
              | ok bucket =>
                  simp only [dif_pos trivial, hconsume] at hrun
                  cases hchecked : pairVecDivVHCExtractChecked heap bucket.nodes with
                  | error fault => simp [hchecked] at hrun
                  | ok extracted =>
                      rw [hchecked] at hrun
                      have hraw := pairVecDivVHCExtractChecked_raw heap
                        bucket.nodes extracted hchecked
                      rcases haway with ⟨owners, hownership, hseparated⟩
                      have hstep :=
                        pairVecDivVHCConsumeRootExtract_preserves_away_linReady
                          this heap extracted.1 k nodes lin resetH quotient
                          divisor bucket owners hheap hownership hseparated
                          hlinReady hconsume hraw
                      have hsmaller : extracted.1.size < size := by
                        rw [← hsize]
                        omega
                      exact ih extracted.1.size hsmaller extracted.1
                        bucket.coefficient bucket.nodes bucket.lin bucket.resetH
                        result hstep.1 hstep.2 hrun rfl
            · simp only [hequal, ↓reduceDIte, Except.ok.injEq] at hrun
              subst result
              exact ⟨haway, hlinReady⟩
      · simp only [hheap, ↓reduceDIte, Except.ok.injEq] at hrun
        subst result
        exact ⟨haway, hlinReady⟩

theorem pairVecDivVHCConsumeEqualDegree_preserves_node_invariants
    (this : DenseUPolyZp) (p degreeLimit degree : Nat)
    (heap : Array Nat) (k : UInt64) (nodes : Array PairVecDivVHCNode)
    (lin : Array Nat) (resetH : Nat) (quotient divisor : SparsePolyZp)
    (result : PairVecDivVHCEqualDegreeResult) (owners : Nat → Finset Nat)
    (hcanonical : SparsePolyZp.Canonical p quotient)
    (hbelow : PairVecDivVHCAllActiveNodesBelow degreeLimit nodes)
    (hdenotes : ∀ (i : Nat) (node : PairVecDivVHCNode),
      nodes[i]? = some node → node.mono ≠ none →
        PairVecDivVHCNodeDenotes quotient divisor node)
    (hfixed : PairVecDivVHCNodeDivisorIndicesFixed nodes)
    (hready : PairVecDivVHCResetReady resetH quotient.size nodes)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hrun : pairVecDivVHCConsumeEqualDegree this degree heap k nodes lin resetH
      quotient divisor = .ok result) :
    PairVecDivVHCAllActiveNodesBelow degreeLimit result.nodes ∧
      (∀ (i : Nat) (node : PairVecDivVHCNode),
        result.nodes[i]? = some node → node.mono ≠ none →
          PairVecDivVHCNodeDenotes quotient divisor node) ∧
      PairVecDivVHCNodeDivisorIndicesFixed result.nodes ∧
      PairVecDivVHCResetReady result.resetH quotient.size result.nodes := by
  induction hsize : heap.size using Nat.strong_induction_on generalizing heap k
      nodes lin resetH result with
  | h size ih =>
      rw [pairVecDivVHCConsumeEqualDegree] at hrun
      by_cases hheap : 0 < heap.size
      · simp only [hheap, ↓reduceDIte] at hrun
        cases hmono : pairVecDivVHCMono heap[0] nodes with
        | error fault => simp [hmono] at hrun
        | ok rootMono =>
            simp only [hmono] at hrun
            by_cases hequal : rootMono.deg = degree
            · simp only [hequal, ↓reduceDIte] at hrun
              cases hconsume : pairVecDivVHCConsumeRootBucket this heap k nodes
                  lin resetH quotient divisor with
              | error fault => simp [hconsume] at hrun
              | ok bucket =>
                  simp only [dif_pos trivial, hconsume] at hrun
                  cases hchecked : pairVecDivVHCExtractChecked heap bucket.nodes with
                  | error fault => simp [hchecked] at hrun
                  | ok extracted =>
                      rw [hchecked] at hrun
                      have hrootOwns :=
                        pairVecDivVHCHeapChainOwnership_root_owns heap owners
                          nodes hownership hheap
                      have hrootValidOwner := pairVecDivVHCChainOwns_valid
                        (some heap[0]) (owners heap[0]) nodes hrootOwns
                      have hrootValid := pairVecDivVHCChainValid_mono
                        (some heap[0]) (owners heap[0])
                        (Finset.range nodes.size) nodes hrootValidOwner
                        (pairVecDivVHCChainOwns_subset_range _ _ _ hrootOwns)
                      have hinvariants :=
                        pairVecDivVHCConsumeRootBucket_preserves_node_invariants
                          this p degreeLimit heap k nodes lin resetH quotient
                          divisor bucket hheap hcanonical hbelow hdenotes hfixed
                          hready hrootValid hconsume
                      have hraw := pairVecDivVHCExtractChecked_raw heap
                        bucket.nodes extracted hchecked
                      have hownership' :=
                        pairVecDivVHCConsumeRootExtract_preserves_heapChainOwnership
                          this heap extracted.1 k nodes lin resetH quotient
                          divisor bucket owners hheap hownership hconsume hraw
                      have hsmaller : extracted.1.size < size := by
                        rw [← hsize]
                        omega
                      exact ih extracted.1.size hsmaller extracted.1
                        bucket.coefficient bucket.nodes bucket.lin bucket.resetH
                        result hinvariants.1 hinvariants.2.1 hinvariants.2.2.1
                        hinvariants.2.2.2 hownership' hrun rfl
            · simp only [hequal, ↓reduceDIte, Except.ok.injEq] at hrun
              subst result
              exact ⟨hbelow, hdenotes, hfixed, hready⟩
      · simp only [hheap, ↓reduceDIte, Except.ok.injEq] at hrun
        subst result
        exact ⟨hbelow, hdenotes, hfixed, hready⟩

theorem pairVecDivVHCActivatedTail_degree_lt_frontier (p frontierDegree : Nat)
    (quotient divisor : SparsePolyZp) (node : PairVecDivVHCNode)
    (hdivisorCanonical : SparsePolyZp.Canonical p divisor)
    (hdivisorNonempty : 0 < divisor.size)
    (hq : node.quotientIndex < quotient.size)
    (hd : node.divisorIndex < divisor.size)
    (htail : 0 < node.divisorIndex)
    (hlead : divisor[0].1.deg ≤ frontierDegree)
    (hquotientDegree : quotient[node.quotientIndex].1.deg =
      frontierDegree - divisor[0].1.deg) :
    quotient[node.quotientIndex].1.deg +
        divisor[node.divisorIndex].1.deg < frontierDegree := by
  have htailDegree := canonical_degree_lt_of_index_lt p divisor
    hdivisorCanonical 0 node.divisorIndex hdivisorNonempty hd htail
  omega

theorem pairVecDivVHCActivate_preserves_allActiveNodesBelow
    (degreeLimit nodeIndex : Nat) (nodes nodes' : Array PairVecDivVHCNode)
    (quotient divisor : SparsePolyZp)
    (hbelow : PairVecDivVHCAllActiveNodesBelow degreeLimit nodes)
    (hnewBelow : ∀ (node : PairVecDivVHCNode),
      nodes[nodeIndex]? = some node →
      ∀ (hq : node.quotientIndex < quotient.size)
        (hd : node.divisorIndex < divisor.size),
        quotient[node.quotientIndex].1.deg +
          divisor[node.divisorIndex].1.deg < degreeLimit)
    (hrun : pairVecDivVHCActivate nodeIndex nodes quotient divisor =
      .ok nodes') :
    PairVecDivVHCAllActiveNodesBelow degreeLimit nodes' := by
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
        intro i node mono hget hmono
        by_cases heq : nodeIndex = i
        · subst i
          rw [Array.getElem?_set_self hn] at hget
          simp only [Option.some.injEq] at hget
          subst node
          have hdegree := hnewBelow nodes[nodeIndex]
            (Array.getElem?_eq_getElem hn) hq hd
          have hmonoDegree := congrArg UMonomial.deg (Option.some.inj hmono)
          dsimp only at hmonoDegree
          omega
        · rw [Array.getElem?_set_ne hn heq] at hget
          exact hbelow i node mono hget hmono

theorem pairVecDivVHCActivate_preserves_denotes
    (nodeIndex : Nat) (nodes nodes' : Array PairVecDivVHCNode)
    (quotient divisor : SparsePolyZp)
    (hdenotes : ∀ (i : Nat) (node : PairVecDivVHCNode),
      nodes[i]? = some node → node.mono ≠ none →
        PairVecDivVHCNodeDenotes quotient divisor node)
    (hrun : pairVecDivVHCActivate nodeIndex nodes quotient divisor =
      .ok nodes') :
    ∀ (i : Nat) (node : PairVecDivVHCNode),
      nodes'[i]? = some node → node.mono ≠ none →
        PairVecDivVHCNodeDenotes quotient divisor node := by
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
        intro i node hget hactive
        by_cases heq : nodeIndex = i
        · subst i
          rw [Array.getElem?_set_self hn] at hget
          simp only [Option.some.injEq] at hget
          subst node
          exact ⟨quotient[nodes[nodeIndex].quotientIndex],
            divisor[nodes[nodeIndex].divisorIndex],
            Array.getElem?_eq_getElem hq, Array.getElem?_eq_getElem hd, rfl⟩
        · rw [Array.getElem?_set_ne hn heq] at hget
          exact hdenotes i node hget hactive

theorem pairVecDivVHCActivate_preserves_divisorIndicesFixed
    (nodeIndex : Nat) (nodes nodes' : Array PairVecDivVHCNode)
    (quotient divisor : SparsePolyZp)
    (hfixed : PairVecDivVHCNodeDivisorIndicesFixed nodes)
    (hrun : pairVecDivVHCActivate nodeIndex nodes quotient divisor =
      .ok nodes') :
    PairVecDivVHCNodeDivisorIndicesFixed nodes' := by
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
        exact PairVecDivVHCNodeDivisorIndicesFixed.set nodes nodeIndex _ hn
          hfixed rfl

theorem pairVecDivVHCActivate_shrinks_resetReady
    (resetH resetH' quotientSize nodeIndex : Nat)
    (nodes nodes' : Array PairVecDivVHCNode)
    (quotient divisor : SparsePolyZp)
    (hready : PairVecDivVHCResetReady resetH quotientSize nodes)
    (hshrink : resetH' ≤ resetH) (houtside : resetH' ≤ nodeIndex)
    (hrun : pairVecDivVHCActivate nodeIndex nodes quotient divisor =
      .ok nodes') :
    PairVecDivVHCResetReady resetH' quotientSize nodes' := by
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
        exact PairVecDivVHCResetReady.shrink_set resetH resetH' quotientSize
          nodeIndex nodes _ hn hready hshrink houtside

theorem pairVecDivVHCActivateReset_preserves_divisorIndicesFixed
    (resetH : Nat) (heap : Array Nat) (nodes : Array PairVecDivVHCNode)
    (quotient divisor : SparsePolyZp) (state : PairVecDivVHCHeapState)
    (hfixed : PairVecDivVHCNodeDivisorIndicesFixed nodes)
    (hrun : pairVecDivVHCActivateReset resetH heap nodes quotient divisor =
      .ok state) :
    PairVecDivVHCNodeDivisorIndicesFixed state.nodes := by
  rw [pairVecDivVHCActivateReset] at hrun
  split at hrun
  next hreset =>
    dsimp only at hrun
    cases hactivate : pairVecDivVHCActivate (resetH - 1) nodes quotient
        divisor with
    | error fault => simp [hactivate] at hrun
    | ok nodes' =>
        simp only [hactivate] at hrun
        cases hinsert : pairVecDivVHCInsert (resetH - 1) heap nodes' with
        | error fault => simp [hinsert] at hrun
        | ok inserted =>
            rcases inserted with ⟨heap', nodes''⟩
            simp only [hinsert] at hrun
            have hfixed' :=
              pairVecDivVHCActivate_preserves_divisorIndicesFixed
                (resetH - 1) nodes nodes' quotient divisor hfixed hactivate
            have hfixed'' := pairVecDivVHCInsert_preserves_divisorIndicesFixed
              (resetH - 1) heap heap' nodes' nodes'' hfixed' hinsert
            exact pairVecDivVHCActivateReset_preserves_divisorIndicesFixed
              (resetH - 1) heap' nodes'' quotient divisor state hfixed'' hrun
  next hreset =>
    simp only [Except.ok.injEq] at hrun
    subst state
    exact hfixed
termination_by resetH
decreasing_by omega

theorem pairVecDivVHCActivateReset_preserves_cursorPrefixAbove
    (degreeLimit resetH : Nat) (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode) (quotient divisor : SparsePolyZp)
    (state : PairVecDivVHCHeapState)
    (hprefix : PairVecDivVHCCursorPrefixAbove degreeLimit nodes quotient divisor)
    (hrun : pairVecDivVHCActivateReset resetH heap nodes quotient divisor =
      .ok state) :
    PairVecDivVHCCursorPrefixAbove degreeLimit state.nodes quotient divisor := by
  rw [pairVecDivVHCActivateReset] at hrun
  split at hrun
  next hreset =>
    dsimp only at hrun
    cases hactivate : pairVecDivVHCActivate (resetH - 1) nodes quotient
        divisor with
    | error fault => simp [hactivate] at hrun
    | ok nodes' =>
        simp only [hactivate] at hrun
        cases hinsert : pairVecDivVHCInsert (resetH - 1) heap nodes' with
        | error fault => simp [hinsert] at hrun
        | ok inserted =>
            rcases inserted with ⟨heap', nodes''⟩
            simp only [hinsert] at hrun
            have hprefix' :=
              pairVecDivVHCActivate_preserves_cursorPrefixAbove degreeLimit
                (resetH - 1) nodes nodes' quotient divisor hprefix hactivate
            have hprefix'' :=
              pairVecDivVHCInsert_preserves_cursorPrefixAbove degreeLimit
                (resetH - 1) heap heap' nodes' nodes'' quotient divisor
                hprefix' hinsert
            exact pairVecDivVHCActivateReset_preserves_cursorPrefixAbove
              degreeLimit (resetH - 1) heap' nodes'' quotient divisor state
              hprefix'' hrun
  next hreset =>
    simp only [Except.ok.injEq] at hrun
    subst state
    exact hprefix
termination_by resetH
decreasing_by omega

theorem pairVecDivVHCActivateReset_clears_resetReady
    (resetH quotientSize : Nat) (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode) (quotient divisor : SparsePolyZp)
    (state : PairVecDivVHCHeapState)
    (hready : PairVecDivVHCResetReady resetH quotientSize nodes)
    (hrun : pairVecDivVHCActivateReset resetH heap nodes quotient divisor =
      .ok state) :
    PairVecDivVHCResetReady 0 quotientSize state.nodes := by
  rw [pairVecDivVHCActivateReset] at hrun
  split at hrun
  next hreset =>
    dsimp only at hrun
    cases hactivate : pairVecDivVHCActivate (resetH - 1) nodes quotient
        divisor with
    | error fault => simp [hactivate] at hrun
    | ok nodes' =>
        simp only [hactivate] at hrun
        cases hinsert : pairVecDivVHCInsert (resetH - 1) heap nodes' with
        | error fault => simp [hinsert] at hrun
        | ok inserted =>
            rcases inserted with ⟨heap', nodes''⟩
            simp only [hinsert] at hrun
            have hready' := pairVecDivVHCActivate_shrinks_resetReady resetH
              (resetH - 1) quotientSize (resetH - 1) nodes nodes' quotient
              divisor hready (by omega) (Nat.le_refl _) hactivate
            have hready'' := pairVecDivVHCInsert_preserves_resetReady
              (resetH - 1) quotientSize (resetH - 1) heap heap' nodes'
              nodes'' hready' (Nat.le_refl _) hinsert
            exact pairVecDivVHCActivateReset_clears_resetReady
              (resetH - 1) quotientSize heap' nodes'' quotient divisor state
              hready'' hrun
  next hreset =>
    simp only [Except.ok.injEq] at hrun
    subst state
    simpa [PairVecDivVHCResetReady]
termination_by resetH
decreasing_by omega

theorem pairVecDivVHCActivateReset_preserves_heapChainsOwned
    (resetH quotientSize : Nat) (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode) (quotient divisor : SparsePolyZp)
    (state : PairVecDivVHCHeapState)
    (howned : PairVecDivVHCHeapChainsOwned heap nodes)
    (hready : PairVecDivVHCResetReady resetH quotientSize nodes)
    (hrun : pairVecDivVHCActivateReset resetH heap nodes quotient divisor =
      .ok state) :
    PairVecDivVHCHeapChainsOwned state.heap state.nodes := by
  rw [pairVecDivVHCActivateReset] at hrun
  split at hrun
  next hreset =>
    dsimp only at hrun
    have hindex : resetH - 1 < resetH := by omega
    rcases hready.2 (resetH - 1) hindex with
      ⟨oldNode, holdGet, hqIndex, hdIndex, hinactive⟩
    cases hactivate : pairVecDivVHCActivate (resetH - 1) nodes quotient
        divisor with
    | error fault => simp [hactivate] at hrun
    | ok activated =>
        simp only [hactivate] at hrun
        cases hinsert : pairVecDivVHCInsert (resetH - 1) heap activated with
        | error fault => simp [hinsert] at hrun
        | ok inserted =>
            rcases inserted with ⟨heap', nodes'⟩
            simp only [hinsert] at hrun
            have howned' := pairVecDivVHCActivateInsert_preserves_heapChainsOwned
              (resetH - 1) heap heap' nodes activated nodes' quotient divisor
              oldNode howned holdGet hinactive hactivate hinsert
            have hready' := pairVecDivVHCActivate_shrinks_resetReady resetH
              (resetH - 1) quotientSize (resetH - 1) nodes activated quotient
              divisor hready (by omega) (Nat.le_refl _) hactivate
            have hready'' := pairVecDivVHCInsert_preserves_resetReady
              (resetH - 1) quotientSize (resetH - 1) heap heap' activated
              nodes' hready' (Nat.le_refl _) hinsert
            exact pairVecDivVHCActivateReset_preserves_heapChainsOwned
              (resetH - 1) quotientSize heap' nodes' quotient divisor state
              howned' hready'' hrun
  next hreset =>
    simp only [Except.ok.injEq] at hrun
    subst state
    exact howned
termination_by resetH
decreasing_by omega

theorem pairVecDivVHCActivateReset_preserves_heapChainsHomogeneous
    (resetH quotientSize : Nat) (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode) (quotient divisor : SparsePolyZp)
    (state : PairVecDivVHCHeapState) (owners : Nat → Finset Nat)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hhomogeneous : PairVecDivVHCHeapChainsHomogeneous heap owners nodes)
    (hready : PairVecDivVHCResetReady resetH quotientSize nodes)
    (hrun : pairVecDivVHCActivateReset resetH heap nodes quotient divisor =
      .ok state) :
    ∃ owners', PairVecDivVHCHeapChainOwnership state.heap owners' state.nodes ∧
      PairVecDivVHCHeapChainsHomogeneous state.heap owners' state.nodes := by
  rw [pairVecDivVHCActivateReset] at hrun
  split at hrun
  next hreset =>
    dsimp only at hrun
    have hindex : resetH - 1 < resetH := by omega
    rcases hready.2 (resetH - 1) hindex with
      ⟨oldNode, holdGet, hqIndex, hdIndex, hinactive⟩
    cases hactivate : pairVecDivVHCActivate (resetH - 1) nodes quotient
        divisor with
    | error fault => simp [hactivate] at hrun
    | ok activated =>
        simp only [hactivate] at hrun
        cases hinsert : pairVecDivVHCInsert (resetH - 1) heap activated with
        | error fault => simp [hinsert] at hrun
        | ok inserted =>
            rcases inserted with ⟨heap', nodes'⟩
            simp only [hinsert] at hrun
            rcases pairVecDivVHCActivateInsert_preserves_heapChainsHomogeneous
                (resetH - 1) heap heap' nodes activated nodes' quotient divisor
                oldNode owners hownership hhomogeneous holdGet hinactive
                hactivate hinsert with
              ⟨owners', hownership', hhomogeneous'⟩
            have hready' := pairVecDivVHCActivate_shrinks_resetReady resetH
              (resetH - 1) quotientSize (resetH - 1) nodes activated quotient
              divisor hready (by omega) (Nat.le_refl _) hactivate
            have hready'' := pairVecDivVHCInsert_preserves_resetReady
              (resetH - 1) quotientSize (resetH - 1) heap heap' activated
              nodes' hready' (Nat.le_refl _) hinsert
            exact pairVecDivVHCActivateReset_preserves_heapChainsHomogeneous
              (resetH - 1) quotientSize heap' nodes' quotient divisor state
              owners' hownership' hhomogeneous' hready'' hrun
  next hreset =>
    simp only [Except.ok.injEq] at hrun
    subst state
    exact ⟨owners, hownership, hhomogeneous⟩
termination_by resetH
decreasing_by omega

theorem pairVecDivVHCActivateReset_preserves_heapOrdered
    (resetH quotientSize : Nat) (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode) (quotient divisor : SparsePolyZp)
    (state : PairVecDivVHCHeapState)
    (howned : PairVecDivVHCHeapChainsOwned heap nodes)
    (hready : PairVecDivVHCResetReady resetH quotientSize nodes)
    (hordered : PairVecDivVHCHeapOrdered heap nodes)
    (hrun : pairVecDivVHCActivateReset resetH heap nodes quotient divisor =
      .ok state) :
    PairVecDivVHCHeapOrdered state.heap state.nodes := by
  rw [pairVecDivVHCActivateReset] at hrun
  split at hrun
  next hreset =>
    dsimp only at hrun
    have hindex : resetH - 1 < resetH := by omega
    rcases hready.2 (resetH - 1) hindex with
      ⟨oldNode, holdGet, hqIndex, hdIndex, hinactive⟩
    cases hactivate : pairVecDivVHCActivate (resetH - 1) nodes quotient
        divisor with
    | error fault => simp [hactivate] at hrun
    | ok activated =>
        simp only [hactivate] at hrun
        cases hinsert : pairVecDivVHCInsert (resetH - 1) heap activated with
        | error fault => simp [hinsert] at hrun
        | ok inserted =>
            rcases inserted with ⟨heap', nodes'⟩
            simp only [hinsert] at hrun
            have hordered' := pairVecDivVHCActivateInsert_preserves_heapOrdered
              (resetH - 1) heap heap' nodes activated nodes' quotient divisor
              oldNode howned hordered holdGet hinactive hactivate hinsert
            have howned' := pairVecDivVHCActivateInsert_preserves_heapChainsOwned
              (resetH - 1) heap heap' nodes activated nodes' quotient divisor
              oldNode howned holdGet hinactive hactivate hinsert
            have hready' := pairVecDivVHCActivate_shrinks_resetReady resetH
              (resetH - 1) quotientSize (resetH - 1) nodes activated quotient
              divisor hready (by omega) (Nat.le_refl _) hactivate
            have hready'' := pairVecDivVHCInsert_preserves_resetReady
              (resetH - 1) quotientSize (resetH - 1) heap heap' activated
              nodes' hready' (Nat.le_refl _) hinsert
            exact pairVecDivVHCActivateReset_preserves_heapOrdered
              (resetH - 1) quotientSize heap' nodes' quotient divisor state
              howned' hready'' hordered' hrun
  next hreset =>
    simp only [Except.ok.injEq] at hrun
    subst state
    exact hordered
termination_by resetH
decreasing_by omega

theorem pairVecDivVHCActivateReset_preserves_away_linReady
    (resetH quotientSize : Nat) (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode) (lin : Array Nat)
    (quotient divisor : SparsePolyZp) (state : PairVecDivVHCHeapState)
    (haway : PairVecDivVHCHeapChainsOwnedAway heap nodes
      lin.toList.toFinset)
    (hlinReady : PairVecDivVHCLinReady lin nodes)
    (hready : PairVecDivVHCResetReady resetH quotientSize nodes)
    (hrun : pairVecDivVHCActivateReset resetH heap nodes quotient divisor =
      .ok state) :
    PairVecDivVHCHeapChainsOwnedAway state.heap state.nodes
        lin.toList.toFinset ∧
      PairVecDivVHCLinReady lin state.nodes := by
  rw [pairVecDivVHCActivateReset] at hrun
  split at hrun
  next hreset =>
    dsimp only at hrun
    have hindex : resetH - 1 < resetH := by omega
    rcases hready.2 (resetH - 1) hindex with
      ⟨oldNode, holdGet, hqIndex, hdIndex, hinactive⟩
    have hnotLin : resetH - 1 ∉ lin.toList := by
      intro hmem
      rcases hlinReady.2 (resetH - 1) hmem with
        ⟨activeNode, mono, hactiveGet, hactiveMono⟩
      rw [holdGet] at hactiveGet
      simp only [Option.some.injEq] at hactiveGet
      subst activeNode
      rw [hinactive] at hactiveMono
      contradiction
    cases hactivate : pairVecDivVHCActivate (resetH - 1) nodes quotient
        divisor with
    | error fault => simp [hactivate] at hrun
    | ok activated =>
        simp only [hactivate] at hrun
        cases hinsert : pairVecDivVHCInsert (resetH - 1) heap activated with
        | error fault => simp [hinsert] at hrun
        | ok inserted =>
            rcases inserted with ⟨heap', nodes'⟩
            simp only [hinsert] at hrun
            have haway' := pairVecDivVHCActivateInsert_preserves_away
              (resetH - 1) heap heap' nodes activated nodes' quotient divisor
              oldNode lin.toList.toFinset haway (by
                intro i hi
                exact hlinReady.2 i (by
                  simpa only [List.mem_toFinset] using hi)) holdGet hinactive
              hactivate hinsert
            have hlinReady' : PairVecDivVHCLinReady lin nodes' := by
              refine ⟨hlinReady.1, ?_⟩
              intro i hi
              rcases hlinReady.2 i hi with ⟨node, mono, hnode, hmono⟩
              have hne : resetH - 1 ≠ i := by
                intro heq
                subst i
                exact hnotLin hi
              refine ⟨node, mono, ?_, hmono⟩
              rw [pairVecDivVHCInsert_get_ne (resetH - 1) heap heap'
                activated nodes' hinsert i hne]
              rw [pairVecDivVHCActivate_get_ne (resetH - 1) nodes activated
                quotient divisor hactivate i hne]
              exact hnode
            have hready' := pairVecDivVHCActivate_shrinks_resetReady resetH
              (resetH - 1) quotientSize (resetH - 1) nodes activated quotient
              divisor hready (by omega) (Nat.le_refl _) hactivate
            have hready'' := pairVecDivVHCInsert_preserves_resetReady
              (resetH - 1) quotientSize (resetH - 1) heap heap' activated
              nodes' hready' (Nat.le_refl _) hinsert
            exact pairVecDivVHCActivateReset_preserves_away_linReady
              (resetH - 1) quotientSize heap' nodes' lin quotient divisor state
              haway' hlinReady' hready'' hrun
  next hreset =>
    simp only [Except.ok.injEq] at hrun
    subst state
    exact ⟨haway, hlinReady⟩
termination_by resetH
decreasing_by omega

theorem pairVecDivVHCActivateReset_preserves_stateCovered
    (resetH quotientSize : Nat) (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode) (lin : Array Nat)
    (quotient divisor : SparsePolyZp) (state : PairVecDivVHCHeapState)
    (haway : PairVecDivVHCHeapChainsOwnedAway heap nodes
      lin.toList.toFinset)
    (hlinReady : PairVecDivVHCLinReady lin nodes)
    (hready : PairVecDivVHCResetReady resetH quotientSize nodes)
    (hcovered : PairVecDivVHCStateCovered heap nodes lin resetH)
    (hrun : pairVecDivVHCActivateReset resetH heap nodes quotient divisor =
      .ok state) :
    PairVecDivVHCStateCovered state.heap state.nodes lin 0 := by
  rw [pairVecDivVHCActivateReset] at hrun
  split at hrun
  next hreset =>
    dsimp only at hrun
    have hindex : resetH - 1 < resetH := by omega
    rcases hready.2 (resetH - 1) hindex with
      ⟨oldNode, holdGet, hqIndex, hdIndex, hinactive⟩
    have hnotLin : resetH - 1 ∉ lin.toList := by
      intro hmem
      rcases hlinReady.2 (resetH - 1) hmem with
        ⟨activeNode, mono, hactiveGet, hactiveMono⟩
      rw [holdGet] at hactiveGet
      simp only [Option.some.injEq] at hactiveGet
      subst activeNode
      rw [hinactive] at hactiveMono
      contradiction
    cases hactivate : pairVecDivVHCActivate (resetH - 1) nodes quotient
        divisor with
    | error fault => simp [hactivate] at hrun
    | ok activated =>
        simp only [hactivate] at hrun
        cases hinsert : pairVecDivVHCInsert (resetH - 1) heap activated with
        | error fault => simp [hinsert] at hrun
        | ok inserted =>
            rcases inserted with ⟨heap', nodes'⟩
            simp only [hinsert] at hrun
            have haway' := pairVecDivVHCActivateInsert_preserves_away
              (resetH - 1) heap heap' nodes activated nodes' quotient divisor
              oldNode lin.toList.toFinset haway (by
                intro i hi
                exact hlinReady.2 i (by
                  simpa only [List.mem_toFinset] using hi)) holdGet hinactive
              hactivate hinsert
            have hlinReady' : PairVecDivVHCLinReady lin nodes' := by
              refine ⟨hlinReady.1, ?_⟩
              intro i hi
              rcases hlinReady.2 i hi with ⟨node, mono, hnode, hmono⟩
              have hne : resetH - 1 ≠ i := by
                intro heq
                subst i
                exact hnotLin hi
              refine ⟨node, mono, ?_, hmono⟩
              rw [pairVecDivVHCInsert_get_ne (resetH - 1) heap heap'
                activated nodes' hinsert i hne]
              rw [pairVecDivVHCActivate_get_ne (resetH - 1) nodes activated
                quotient divisor hactivate i hne]
              exact hnode
            have hready' := pairVecDivVHCActivate_shrinks_resetReady resetH
              (resetH - 1) quotientSize (resetH - 1) nodes activated quotient
              divisor hready (by omega) (Nat.le_refl _) hactivate
            have hready'' := pairVecDivVHCInsert_preserves_resetReady
              (resetH - 1) quotientSize (resetH - 1) heap heap' activated
              nodes' hready' (Nat.le_refl _) hinsert
            have hcovered' :=
              pairVecDivVHCActivateInsert_preserves_stateCovered resetH
                quotientSize heap heap' nodes activated nodes' lin quotient
                divisor hreset haway hlinReady hready hcovered hactivate hinsert
            exact pairVecDivVHCActivateReset_preserves_stateCovered
              (resetH - 1) quotientSize heap' nodes' lin quotient divisor state
              haway' hlinReady' hready'' hcovered' hrun
  next hreset =>
    simp only [Except.ok.injEq] at hrun
    subst state
    have hzero : resetH = 0 := by omega
    simpa [hzero] using hcovered
termination_by resetH
decreasing_by omega

theorem pairVecDivVHCActivateReset_preserves_node_invariants
    (p frontierDegree degreeLimit resetH quotientSize : Nat)
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode)
    (quotient divisor : SparsePolyZp) (state : PairVecDivVHCHeapState)
    (hdivisorCanonical : SparsePolyZp.Canonical p divisor)
    (hdivisorNonempty : 0 < divisor.size)
    (hlead : divisor[0].1.deg ≤ frontierDegree)
    (hfrontierBelow : frontierDegree ≤ degreeLimit)
    (hquotientIndex : quotientSize < quotient.size)
    (hquotientDegree : quotient[quotientSize].1.deg =
      frontierDegree - divisor[0].1.deg)
    (hbelow : PairVecDivVHCAllActiveNodesBelow degreeLimit nodes)
    (hdenotes : ∀ (i : Nat) (node : PairVecDivVHCNode),
      nodes[i]? = some node → node.mono ≠ none →
        PairVecDivVHCNodeDenotes quotient divisor node)
    (hfixed : PairVecDivVHCNodeDivisorIndicesFixed nodes)
    (hready : PairVecDivVHCResetReady resetH quotientSize nodes)
    (hrun : pairVecDivVHCActivateReset resetH heap nodes quotient divisor =
      .ok state) :
    PairVecDivVHCAllActiveNodesBelow degreeLimit state.nodes ∧
      (∀ (i : Nat) (node : PairVecDivVHCNode),
        state.nodes[i]? = some node → node.mono ≠ none →
          PairVecDivVHCNodeDenotes quotient divisor node) ∧
      PairVecDivVHCNodeDivisorIndicesFixed state.nodes ∧
      PairVecDivVHCResetReady 0 quotientSize state.nodes := by
  rw [pairVecDivVHCActivateReset] at hrun
  split at hrun
  next hreset =>
    let nodeIndex := resetH - 1
    have hnodeIndexLt : nodeIndex < resetH := by
      dsimp only [nodeIndex]
      omega
    rcases hready.2 nodeIndex hnodeIndexLt with
      ⟨readyNode, hreadyGet, hreadyQuotient, hreadyDivisor, hreadyMono⟩
    dsimp only at hrun
    cases hactivate : pairVecDivVHCActivate nodeIndex nodes quotient divisor with
    | error fault => simp [nodeIndex, hactivate] at hrun
    | ok nodes' =>
        simp only [nodeIndex, hactivate] at hrun
        cases hinsert : pairVecDivVHCInsert nodeIndex heap nodes' with
        | error fault => simp [nodeIndex, hinsert] at hrun
        | ok inserted =>
            rcases inserted with ⟨heap', nodes''⟩
            simp only [nodeIndex, hinsert] at hrun
            have hactivationBelow : ∀ (node : PairVecDivVHCNode),
                nodes[nodeIndex]? = some node →
                ∀ (hq : node.quotientIndex < quotient.size)
                  (hd : node.divisorIndex < divisor.size),
                  quotient[node.quotientIndex].1.deg +
                    divisor[node.divisorIndex].1.deg < degreeLimit := by
              intro node hget hq hd
              rw [hreadyGet] at hget
              simp only [Option.some.injEq] at hget
              subst node
              have htail : 0 < readyNode.divisorIndex := by
                rw [hreadyDivisor]
                omega
              have hdegree := pairVecDivVHCActivatedTail_degree_lt_frontier p
                frontierDegree quotient divisor readyNode hdivisorCanonical
                hdivisorNonempty hq hd htail hlead (by
                  simpa only [hreadyQuotient] using hquotientDegree)
              exact Nat.lt_of_lt_of_le hdegree hfrontierBelow
            have hbelow' :=
              pairVecDivVHCActivate_preserves_allActiveNodesBelow degreeLimit
                nodeIndex nodes nodes' quotient divisor hbelow
                hactivationBelow hactivate
            have hdenotes' := pairVecDivVHCActivate_preserves_denotes nodeIndex
              nodes nodes' quotient divisor hdenotes hactivate
            have hfixed' :=
              pairVecDivVHCActivate_preserves_divisorIndicesFixed nodeIndex
                nodes nodes' quotient divisor hfixed hactivate
            have hready' := pairVecDivVHCActivate_shrinks_resetReady resetH
              nodeIndex quotientSize nodeIndex nodes nodes' quotient divisor
              hready (by dsimp [nodeIndex]; omega) (Nat.le_refl _) hactivate
            have hbelow'' := pairVecDivVHCInsert_preserves_allActiveNodesBelow
              degreeLimit nodeIndex heap heap' nodes' nodes'' hbelow' hinsert
            have hdenotes'' := pairVecDivVHCInsert_preserves_denotes nodeIndex
              heap heap' nodes' nodes'' quotient divisor hdenotes' hinsert
            have hfixed'' := pairVecDivVHCInsert_preserves_divisorIndicesFixed
              nodeIndex heap heap' nodes' nodes'' hfixed' hinsert
            have hready'' := pairVecDivVHCInsert_preserves_resetReady nodeIndex
              quotientSize nodeIndex heap heap' nodes' nodes'' hready'
              (Nat.le_refl _) hinsert
            exact pairVecDivVHCActivateReset_preserves_node_invariants p
              frontierDegree degreeLimit nodeIndex quotientSize heap' nodes''
              quotient divisor state hdivisorCanonical hdivisorNonempty hlead
              hfrontierBelow hquotientIndex hquotientDegree hbelow''
              hdenotes'' hfixed'' hready'' hrun
  next hreset =>
    simp only [Except.ok.injEq] at hrun
    subst state
    have hzero : resetH = 0 := by omega
    subst resetH
    exact ⟨hbelow, hdenotes, hfixed, hready⟩
termination_by resetH
decreasing_by
  omega

theorem pairVecDivVHCCanonicalInitialFrontierBelow (p : Nat)
    (dividend : SparsePolyZp) (nodes : Array PairVecDivVHCNode)
    (hcanonical : SparsePolyZp.Canonical p dividend)
    (hnonempty : 0 < dividend.size) :
    PairVecDivVHCFrontierBelow (dividend[0].1.deg + 1) 0 dividend #[]
      nodes :=
  pairVecDivVHCInitialFrontierBelow dividend nodes hnonempty
    (canonical_degree_le_head p dividend hcanonical hnonempty)

theorem pairVecDivVHCCanonicalInitialRemainingBelow (p : Nat)
    (dividend : SparsePolyZp) (hcanonical : SparsePolyZp.Canonical p dividend)
    (hnonempty : 0 < dividend.size) :
    PairVecDivVHCRemainingDividendBelow (dividend[0].1.deg + 1) 0
      dividend := by
  exact (pairVecDivVHCCanonicalInitialFrontierBelow p dividend #[] hcanonical
    hnonempty).1

/-- The selector's well-founded guard is a theorem consequence of the
frontier-bound invariant; it is not an independently assumed execution
budget. -/
theorem pairVecDivVHCDecreaseGuard (degreeLimit dividendIndex : Nat)
    (dividend : SparsePolyZp) (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode) (frontier : PairVecDivVHCFrontier)
    (hbelow : PairVecDivVHCFrontierBelow degreeLimit dividendIndex
      dividend heap nodes)
    (hselect : pairVecDivVHCSelectFrontier dividendIndex dividend heap nodes =
      .ok frontier) :
    frontier.degree < degreeLimit :=
  pairVecDivVHCSelectFrontier_degree_lt degreeLimit dividendIndex dividend
    heap nodes frontier hbelow hselect

structure PairVecDivVHCIterationResult where
  dividendIndex : Nat
  heap : Array Nat
  nodes : Array PairVecDivVHCNode
  quotient : SparsePolyZp
  resetH : Nat

/-- Exact quotient-emission sub-block of one outer iteration. -/
def pairVecDivVHCEmit (this : DenseUPolyZp)
    (frontier : PairVecDivVHCFrontier)
    (consumed : PairVecDivVHCEqualDegreeResult)
    (quotient divisor : SparsePolyZp) (hdivisor : 0 < divisor.size) :
    RawExec (SparsePolyZp × PairVecDivVHCHeapState × Nat) :=
  if consumed.coefficient ≠ 0 then
    if divisor[0].1.deg ≤ frontier.degree then
      let inverse := Generated.StrictGCD.dense_upoly_zp_nmod_inv_ir this
        divisor[0].2.val
      let value := Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this
        consumed.coefficient inverse
      if value ≠ 0 then
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

def PairVecDivVHCQuotientAbove (frontierDegree leadDegree : Nat)
    (quotient : SparsePolyZp) : Prop :=
  leadDegree ≤ frontierDegree →
    ∀ term ∈ quotient.toList, frontierDegree - leadDegree < term.1.deg

theorem PairVecDivVHCQuotientAbove.empty (frontierDegree leadDegree : Nat) :
    PairVecDivVHCQuotientAbove frontierDegree leadDegree #[] := by
  intro hlead term hterm
  simp at hterm

/-- Every emitted quotient term still present at an outer-loop state names a
lead-product degree already processed at or above `degreeLimit`. -/
def PairVecDivVHCQuotientLeadAbove (degreeLimit leadDegree : Nat)
    (quotient : SparsePolyZp) : Prop :=
  ∀ term ∈ quotient.toList, degreeLimit ≤ term.1.deg + leadDegree

theorem PairVecDivVHCQuotientLeadAbove.empty
    (degreeLimit leadDegree : Nat) :
    PairVecDivVHCQuotientLeadAbove degreeLimit leadDegree #[] := by
  intro term hterm
  simp at hterm

theorem pairVecDivVHCEmit_quotient_eq_or_push
    (this : DenseUPolyZp) (frontier : PairVecDivVHCFrontier)
    (consumed : PairVecDivVHCEqualDegreeResult)
    (quotient divisor quotient' : SparsePolyZp)
    (activated : PairVecDivVHCHeapState) (resetH' : Nat)
    (hdivisor : 0 < divisor.size)
    (hrun : pairVecDivVHCEmit this frontier consumed quotient divisor hdivisor =
      .ok (quotient', activated, resetH')) :
    quotient' = quotient ∨
      divisor[0].1.deg ≤ frontier.degree ∧
        ∃ value : UInt64, value ≠ 0 ∧
          quotient' = quotient.push
            (⟨frontier.degree - divisor[0].1.deg⟩, ⟨value, this._p⟩) := by
  unfold pairVecDivVHCEmit at hrun
  by_cases hcoefficient : consumed.coefficient ≠ 0
  · rw [if_pos hcoefficient] at hrun
    by_cases hdegree : divisor[0].1.deg ≤ frontier.degree
    · rw [if_pos hdegree] at hrun
      let inverse := Generated.StrictGCD.dense_upoly_zp_nmod_inv_ir this
        divisor[0].2.val
      let value := Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this
        consumed.coefficient inverse
      by_cases hvalue : value ≠ 0
      · rw [if_pos hvalue] at hrun
        let emitted := quotient.push
          (⟨frontier.degree - divisor[0].1.deg⟩, ⟨value, this._p⟩)
        cases hactivate : pairVecDivVHCActivateReset consumed.resetH
            consumed.heap consumed.nodes emitted divisor with
        | error fault =>
            dsimp only [emitted] at hactivate
            dsimp only at hrun
            rw [hactivate] at hrun
            contradiction
        | ok state =>
            dsimp only [emitted] at hactivate
            dsimp only at hrun
            rw [hactivate] at hrun
            simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
            rcases hrun with ⟨rfl, rfl, rfl⟩
            exact Or.inr ⟨hdegree, value, hvalue, rfl⟩
      · rw [if_neg hvalue] at hrun
        simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
        exact Or.inl hrun.1.symm
    · rw [if_neg hdegree] at hrun
      simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
      exact Or.inl hrun.1.symm
  · rw [if_neg hcoefficient] at hrun
    simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
    exact Or.inl hrun.1.symm

theorem pairVecDivVHCEmit_preserves_quotientLeadAbove
    (this : DenseUPolyZp) (degreeLimit : Nat)
    (frontier : PairVecDivVHCFrontier)
    (consumed : PairVecDivVHCEqualDegreeResult)
    (quotient divisor quotient' : SparsePolyZp)
    (activated : PairVecDivVHCHeapState) (resetH' : Nat)
    (hdivisor : 0 < divisor.size)
    (hprocessed : PairVecDivVHCQuotientLeadAbove degreeLimit
      divisor[0].1.deg quotient)
    (hdecrease : frontier.degree < degreeLimit)
    (hrun : pairVecDivVHCEmit this frontier consumed quotient divisor hdivisor =
      .ok (quotient', activated, resetH')) :
    PairVecDivVHCQuotientLeadAbove frontier.degree divisor[0].1.deg
      quotient' := by
  rcases pairVecDivVHCEmit_quotient_eq_or_push this frontier consumed quotient
      divisor quotient' activated resetH' hdivisor hrun with
    hsame | ⟨hlead, value, hvalue, hpush⟩
  · subst quotient'
    intro term hterm
    exact Nat.le_trans (Nat.le_of_lt hdecrease)
      (hprocessed term hterm)
  · subst quotient'
    intro term hterm
    simp only [Array.toList_push, List.mem_append, List.mem_singleton] at hterm
    rcases hterm with hold | hnew
    · exact Nat.le_trans (Nat.le_of_lt hdecrease)
        (hprocessed term hold)
    · subst term
      dsimp only
      omega

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
    let (quotient', activated, resetH') ←
      pairVecDivVHCEmit this frontier consumed quotient divisor hdivisor
    let reinserted ← pairVecDivVHCReinsertLin activated.heap activated.nodes
      consumed.lin
    .ok (PairVecDivVHCIterationResult.mk frontier.dividendIndex
      reinserted.heap reinserted.nodes quotient' resetH')
  else
    .error .assertionFailure

theorem pairVecDivVHCOuterIteration_dividendIndex
    (this : DenseUPolyZp) (dividendIndex : Nat) (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode) (quotient dividend divisor : SparsePolyZp)
    (resetH : Nat) (frontier : PairVecDivVHCFrontier)
    (result : PairVecDivVHCIterationResult) (hdivisor : 0 < divisor.size)
    (hselect : pairVecDivVHCSelectFrontier dividendIndex dividend heap nodes =
      .ok frontier)
    (hrun : pairVecDivVHCOuterIteration this dividendIndex heap nodes quotient
      dividend divisor resetH = .ok result) :
    result.dividendIndex = frontier.dividendIndex := by
  unfold pairVecDivVHCOuterIteration at hrun
  simp only [hdivisor, ↓reduceDIte, hselect, Bind.bind, Except.bind] at hrun
  generalize hconsume : pairVecDivVHCConsumeEqualDegree this frontier.degree
      heap frontier.coefficient nodes #[] resetH quotient divisor =
        consumeExec at hrun
  cases consumeExec with
  | error fault => contradiction
  | ok consumed =>
      simp only [Bind.bind, Except.bind] at hrun
      generalize hemit : pairVecDivVHCEmit this frontier consumed quotient
          divisor hdivisor = emitExec at hrun
      cases emitExec with
      | error fault => contradiction
      | ok emitted =>
          rcases emitted with ⟨quotient', activated, resetH'⟩
          simp only [Bind.bind, Except.bind] at hrun
          generalize hreinsert : pairVecDivVHCReinsertLin activated.heap
              activated.nodes consumed.lin = reinsertExec at hrun
          cases reinsertExec with
          | error fault => contradiction
          | ok reinserted =>
              simp only [Bind.bind, Except.bind, Except.ok.injEq] at hrun
              rw [← hrun]

theorem pairVecDivVHCOuterIteration_components
    (this : DenseUPolyZp) (dividendIndex : Nat) (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode) (quotient dividend divisor : SparsePolyZp)
    (resetH : Nat) (frontier : PairVecDivVHCFrontier)
    (result : PairVecDivVHCIterationResult) (hdivisor : 0 < divisor.size)
    (hselect : pairVecDivVHCSelectFrontier dividendIndex dividend heap nodes =
      .ok frontier)
    (hrun : pairVecDivVHCOuterIteration this dividendIndex heap nodes quotient
      dividend divisor resetH = .ok result) :
    ∃ consumed quotient' activated resetH' reinserted,
      pairVecDivVHCConsumeEqualDegree this frontier.degree heap
          frontier.coefficient nodes #[] resetH quotient divisor =
        .ok consumed ∧
      pairVecDivVHCEmit this frontier consumed quotient divisor hdivisor =
        .ok (quotient', activated, resetH') ∧
      pairVecDivVHCReinsertLin activated.heap activated.nodes consumed.lin =
        .ok reinserted ∧
      result = PairVecDivVHCIterationResult.mk frontier.dividendIndex
        reinserted.heap reinserted.nodes quotient' resetH' := by
  unfold pairVecDivVHCOuterIteration at hrun
  simp only [hdivisor, ↓reduceDIte, hselect, Bind.bind, Except.bind] at hrun
  cases hconsume : pairVecDivVHCConsumeEqualDegree this frontier.degree heap
      frontier.coefficient nodes #[] resetH quotient divisor with
  | error fault => simp [hconsume] at hrun
  | ok consumed =>
      simp only [hconsume, Bind.bind, Except.bind] at hrun
      cases hemit : pairVecDivVHCEmit this frontier consumed quotient divisor
          hdivisor with
      | error fault => simp [hemit] at hrun
      | ok emitted =>
          rcases emitted with ⟨quotient', activated, resetH'⟩
          simp only [hemit, Bind.bind, Except.bind] at hrun
          cases hreinsert : pairVecDivVHCReinsertLin activated.heap
              activated.nodes consumed.lin with
          | error fault => simp [hreinsert] at hrun
          | ok reinserted =>
              simp only [hreinsert, Bind.bind, Except.bind,
                Except.ok.injEq] at hrun
              exact ⟨consumed, quotient', activated, resetH', reinserted,
                rfl, hemit, hreinsert, hrun.symm⟩

theorem pairVecDivVHCOuterIteration_preserves_consumed_above
    (this : DenseUPolyZp) (degreeLimit dividendIndex : Nat)
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode)
    (quotient dividend divisor : SparsePolyZp) (resetH : Nat)
    (frontier : PairVecDivVHCFrontier)
    (result : PairVecDivVHCIterationResult) (hdivisor : 0 < divisor.size)
    (hconsumed : PairVecDivVHCConsumedDividendAbove degreeLimit dividendIndex
      dividend)
    (hdecrease : frontier.degree < degreeLimit)
    (hselect : pairVecDivVHCSelectFrontier dividendIndex dividend heap nodes =
      .ok frontier)
    (hrun : pairVecDivVHCOuterIteration this dividendIndex heap nodes quotient
      dividend divisor resetH = .ok result) :
    PairVecDivVHCConsumedDividendAbove frontier.degree result.dividendIndex
      dividend := by
  have hnext := pairVecDivVHCSelectFrontier_preserves_consumed_above
    degreeLimit dividendIndex dividend heap nodes frontier hconsumed hdecrease
    hselect
  rw [pairVecDivVHCOuterIteration_dividendIndex this dividendIndex heap nodes
    quotient dividend divisor resetH frontier result hdivisor hselect hrun]
  exact hnext

theorem pairVecDivVHCEmit_preserves_away_linReady
    (this : DenseUPolyZp) (frontier : PairVecDivVHCFrontier)
    (consumed : PairVecDivVHCEqualDegreeResult)
    (quotient divisor quotient' : SparsePolyZp)
    (activated : PairVecDivVHCHeapState) (resetH' : Nat)
    (haway : PairVecDivVHCHeapChainsOwnedAway consumed.heap consumed.nodes
      consumed.lin.toList.toFinset)
    (hlinReady : PairVecDivVHCLinReady consumed.lin consumed.nodes)
    (hready : PairVecDivVHCResetReady consumed.resetH quotient.size
      consumed.nodes)
    (hdivisor : 0 < divisor.size)
    (hrun : pairVecDivVHCEmit this frontier consumed quotient divisor hdivisor =
      .ok (quotient', activated, resetH')) :
    PairVecDivVHCHeapChainsOwnedAway activated.heap activated.nodes
        consumed.lin.toList.toFinset ∧
      PairVecDivVHCLinReady consumed.lin activated.nodes := by
  unfold pairVecDivVHCEmit at hrun
  by_cases hcoefficient : consumed.coefficient ≠ 0
  · rw [if_pos hcoefficient] at hrun
    by_cases hdegree : divisor[0].1.deg ≤ frontier.degree
    · rw [if_pos hdegree] at hrun
      let inverse := Generated.StrictGCD.dense_upoly_zp_nmod_inv_ir this
        divisor[0].2.val
      let value := Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this
        consumed.coefficient inverse
      by_cases hvalue : value ≠ 0
      · rw [if_pos hvalue] at hrun
        let emitted := quotient.push
          (⟨frontier.degree - divisor[0].1.deg⟩, ⟨value, this._p⟩)
        cases hactivate : pairVecDivVHCActivateReset consumed.resetH
            consumed.heap consumed.nodes emitted divisor with
        | error fault =>
            dsimp only [emitted] at hactivate
            dsimp only at hrun
            rw [hactivate] at hrun
            contradiction
        | ok state =>
            dsimp only [emitted] at hactivate
            dsimp only at hrun
            rw [hactivate] at hrun
            simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
            rcases hrun with ⟨rfl, rfl, rfl⟩
            exact pairVecDivVHCActivateReset_preserves_away_linReady
              consumed.resetH quotient.size consumed.heap consumed.nodes
              consumed.lin emitted divisor state haway hlinReady hready
              hactivate
      · rw [if_neg hvalue] at hrun
        simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
        rcases hrun with ⟨rfl, rfl, rfl⟩
        exact ⟨haway, hlinReady⟩
    · rw [if_neg hdegree] at hrun
      simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
      rcases hrun with ⟨rfl, rfl, rfl⟩
      exact ⟨haway, hlinReady⟩
  · rw [if_neg hcoefficient] at hrun
    simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
    rcases hrun with ⟨rfl, rfl, rfl⟩
    exact ⟨haway, hlinReady⟩

theorem pairVecDivVHCEmit_preserves_heapOrdered
    (this : DenseUPolyZp) (frontier : PairVecDivVHCFrontier)
    (consumed : PairVecDivVHCEqualDegreeResult)
    (quotient divisor quotient' : SparsePolyZp)
    (activated : PairVecDivVHCHeapState) (resetH' : Nat)
    (howned : PairVecDivVHCHeapChainsOwned consumed.heap consumed.nodes)
    (hready : PairVecDivVHCResetReady consumed.resetH quotient.size
      consumed.nodes)
    (hordered : PairVecDivVHCHeapOrdered consumed.heap consumed.nodes)
    (hdivisor : 0 < divisor.size)
    (hrun : pairVecDivVHCEmit this frontier consumed quotient divisor hdivisor =
      .ok (quotient', activated, resetH')) :
    PairVecDivVHCHeapOrdered activated.heap activated.nodes := by
  unfold pairVecDivVHCEmit at hrun
  by_cases hcoefficient : consumed.coefficient ≠ 0
  · rw [if_pos hcoefficient] at hrun
    by_cases hdegree : divisor[0].1.deg ≤ frontier.degree
    · rw [if_pos hdegree] at hrun
      let inverse := Generated.StrictGCD.dense_upoly_zp_nmod_inv_ir this
        divisor[0].2.val
      let value := Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this
        consumed.coefficient inverse
      by_cases hvalue : value ≠ 0
      · rw [if_pos hvalue] at hrun
        let emitted := quotient.push
          (⟨frontier.degree - divisor[0].1.deg⟩, ⟨value, this._p⟩)
        cases hactivate : pairVecDivVHCActivateReset consumed.resetH
            consumed.heap consumed.nodes emitted divisor with
        | error fault =>
            dsimp only [emitted] at hactivate
            dsimp only at hrun
            rw [hactivate] at hrun
            contradiction
        | ok state =>
            dsimp only [emitted] at hactivate
            dsimp only at hrun
            rw [hactivate] at hrun
            have hstateOrdered :=
              pairVecDivVHCActivateReset_preserves_heapOrdered consumed.resetH
                quotient.size consumed.heap consumed.nodes emitted divisor state
                howned hready hordered hactivate
            simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
            rcases hrun with ⟨rfl, rfl, rfl⟩
            exact hstateOrdered
      · rw [if_neg hvalue] at hrun
        simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
        rcases hrun with ⟨rfl, rfl, rfl⟩
        exact hordered
    · rw [if_neg hdegree] at hrun
      simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
      rcases hrun with ⟨rfl, rfl, rfl⟩
      exact hordered
  · rw [if_neg hcoefficient] at hrun
    simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
    rcases hrun with ⟨rfl, rfl, rfl⟩
    exact hordered

theorem pairVecDivVHCEmit_preserves_heapChainsHomogeneous
    (this : DenseUPolyZp) (frontier : PairVecDivVHCFrontier)
    (consumed : PairVecDivVHCEqualDegreeResult)
    (quotient divisor quotient' : SparsePolyZp)
    (activated : PairVecDivVHCHeapState) (resetH' : Nat)
    (owners : Nat → Finset Nat)
    (hownership : PairVecDivVHCHeapChainOwnership consumed.heap owners
      consumed.nodes)
    (hhomogeneous : PairVecDivVHCHeapChainsHomogeneous consumed.heap owners
      consumed.nodes)
    (hready : PairVecDivVHCResetReady consumed.resetH quotient.size
      consumed.nodes)
    (hdivisor : 0 < divisor.size)
    (hrun : pairVecDivVHCEmit this frontier consumed quotient divisor hdivisor =
      .ok (quotient', activated, resetH')) :
    ∃ owners', PairVecDivVHCHeapChainOwnership activated.heap owners'
        activated.nodes ∧
      PairVecDivVHCHeapChainsHomogeneous activated.heap owners'
        activated.nodes := by
  unfold pairVecDivVHCEmit at hrun
  by_cases hcoefficient : consumed.coefficient ≠ 0
  · rw [if_pos hcoefficient] at hrun
    by_cases hdegree : divisor[0].1.deg ≤ frontier.degree
    · rw [if_pos hdegree] at hrun
      let inverse := Generated.StrictGCD.dense_upoly_zp_nmod_inv_ir this
        divisor[0].2.val
      let value := Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this
        consumed.coefficient inverse
      by_cases hvalue : value ≠ 0
      · rw [if_pos hvalue] at hrun
        let emitted := quotient.push
          (⟨frontier.degree - divisor[0].1.deg⟩, ⟨value, this._p⟩)
        cases hactivate : pairVecDivVHCActivateReset consumed.resetH
            consumed.heap consumed.nodes emitted divisor with
        | error fault =>
            dsimp only [emitted] at hactivate
            dsimp only at hrun
            rw [hactivate] at hrun
            contradiction
        | ok state =>
            dsimp only [emitted] at hactivate
            dsimp only at hrun
            rw [hactivate] at hrun
            have hstate :=
              pairVecDivVHCActivateReset_preserves_heapChainsHomogeneous
                consumed.resetH quotient.size consumed.heap consumed.nodes
                emitted divisor state owners hownership hhomogeneous hready
                hactivate
            simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
            rcases hrun with ⟨rfl, rfl, rfl⟩
            exact hstate
      · rw [if_neg hvalue] at hrun
        simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
        rcases hrun with ⟨rfl, rfl, rfl⟩
        exact ⟨owners, hownership, hhomogeneous⟩
    · rw [if_neg hdegree] at hrun
      simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
      rcases hrun with ⟨rfl, rfl, rfl⟩
      exact ⟨owners, hownership, hhomogeneous⟩
  · rw [if_neg hcoefficient] at hrun
    simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
    rcases hrun with ⟨rfl, rfl, rfl⟩
    exact ⟨owners, hownership, hhomogeneous⟩

/-- One complete generated outer iteration preserves max-heap order through
consume/extract, optional emission activation, and reverse reinsertion. -/
theorem pairVecDivVHCOuterIteration_preserves_heapOrdered
    (this : DenseUPolyZp) (p degreeLimit dividendIndex : Nat)
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode)
    (quotient dividend divisor : SparsePolyZp) (resetH : Nat)
    (frontier : PairVecDivVHCFrontier)
    (result : PairVecDivVHCIterationResult) (owners : Nat → Finset Nat)
    (hdivisor : 0 < divisor.size)
    (hcanonical : SparsePolyZp.Canonical p quotient)
    (hbelow : PairVecDivVHCAllActiveNodesBelow degreeLimit nodes)
    (hdenotes : ∀ (i : Nat) (node : PairVecDivVHCNode),
      nodes[i]? = some node → node.mono ≠ none →
        PairVecDivVHCNodeDenotes quotient divisor node)
    (hfixed : PairVecDivVHCNodeDivisorIndicesFixed nodes)
    (hready : PairVecDivVHCResetReady resetH quotient.size nodes)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hordered : PairVecDivVHCHeapOrdered heap nodes)
    (hselect : pairVecDivVHCSelectFrontier dividendIndex dividend heap nodes =
      .ok frontier)
    (hrun : pairVecDivVHCOuterIteration this dividendIndex heap nodes quotient
      dividend divisor resetH = .ok result) :
    PairVecDivVHCHeapOrdered result.heap result.nodes := by
  rcases pairVecDivVHCOuterIteration_components this dividendIndex heap nodes
      quotient dividend divisor resetH frontier result hdivisor hselect hrun with
    ⟨consumed, quotient', activated, resetH', reinserted, hconsume, hemit,
      hreinsert, hresult⟩
  have hconsumedOrder :=
    pairVecDivVHCConsumeEqualDegree_preserves_heapOwnershipOrdered this
      frontier.degree heap frontier.coefficient nodes #[] resetH quotient
      divisor consumed owners hownership hordered hconsume
  have hconsumedNodes :=
    pairVecDivVHCConsumeEqualDegree_preserves_node_invariants this p degreeLimit
      frontier.degree heap frontier.coefficient nodes #[] resetH quotient divisor
      consumed owners hcanonical hbelow hdenotes hfixed hready hownership
      hconsume
  have haway0 : PairVecDivVHCHeapChainsOwnedAway heap nodes
      (#[] : Array Nat).toList.toFinset := by
    simpa using (show PairVecDivVHCHeapChainsOwned heap nodes from
      ⟨owners, hownership⟩).away_empty heap nodes
  have hlin0 : PairVecDivVHCLinReady (#[] : Array Nat) nodes := by
    simp [PairVecDivVHCLinReady]
  have hconsumedAway :=
    pairVecDivVHCConsumeEqualDegree_preserves_away_linReady this
      frontier.degree heap frontier.coefficient nodes #[] resetH quotient divisor
      consumed haway0 hlin0 hconsume
  have hemittedAway := pairVecDivVHCEmit_preserves_away_linReady this frontier
    consumed quotient divisor quotient' activated resetH' hconsumedAway.1
    hconsumedAway.2 hconsumedNodes.2.2.2 hdivisor hemit
  have hemittedOrder := pairVecDivVHCEmit_preserves_heapOrdered this frontier
    consumed quotient divisor quotient' activated resetH'
    ⟨owners, hconsumedOrder.1⟩ hconsumedNodes.2.2.2 hconsumedOrder.2
    hdivisor hemit
  have hreinsertedOrder := pairVecDivVHCReinsertLin_preserves_heapOrdered
    activated.heap activated.nodes consumed.lin reinserted hemittedAway.1
    hemittedAway.2 hemittedOrder hreinsert
  rw [hresult]
  exact hreinsertedOrder

theorem pairVecDivVHCOuterIteration_preserves_heapChainsHomogeneous
    (this : DenseUPolyZp) (p degreeLimit dividendIndex : Nat)
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode)
    (quotient dividend divisor : SparsePolyZp) (resetH : Nat)
    (frontier : PairVecDivVHCFrontier)
    (result : PairVecDivVHCIterationResult) (owners : Nat → Finset Nat)
    (hdivisor : 0 < divisor.size)
    (hcanonical : SparsePolyZp.Canonical p quotient)
    (hbelow : PairVecDivVHCAllActiveNodesBelow degreeLimit nodes)
    (hdenotes : ∀ (i : Nat) (node : PairVecDivVHCNode),
      nodes[i]? = some node → node.mono ≠ none →
        PairVecDivVHCNodeDenotes quotient divisor node)
    (hfixed : PairVecDivVHCNodeDivisorIndicesFixed nodes)
    (hready : PairVecDivVHCResetReady resetH quotient.size nodes)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hhomogeneous : PairVecDivVHCHeapChainsHomogeneous heap owners nodes)
    (hselect : pairVecDivVHCSelectFrontier dividendIndex dividend heap nodes =
      .ok frontier)
    (hrun : pairVecDivVHCOuterIteration this dividendIndex heap nodes quotient
      dividend divisor resetH = .ok result) :
    ∃ owners', PairVecDivVHCHeapChainOwnership result.heap owners'
        result.nodes ∧
      PairVecDivVHCHeapChainsHomogeneous result.heap owners' result.nodes := by
  rcases pairVecDivVHCOuterIteration_components this dividendIndex heap nodes
      quotient dividend divisor resetH frontier result hdivisor hselect hrun with
    ⟨consumed, quotient', activated, resetH', reinserted, hconsume, hemit,
      hreinsert, hresult⟩
  have hconsumedHomogeneous :=
    pairVecDivVHCConsumeEqualDegree_preserves_heapChainsHomogeneous this
      frontier.degree heap frontier.coefficient nodes #[] resetH quotient
      divisor consumed owners hownership hhomogeneous hconsume
  have hconsumedNodes :=
    pairVecDivVHCConsumeEqualDegree_preserves_node_invariants this p degreeLimit
      frontier.degree heap frontier.coefficient nodes #[] resetH quotient divisor
      consumed owners hcanonical hbelow hdenotes hfixed hready hownership
      hconsume
  have haway0 : PairVecDivVHCHeapChainsOwnedAway heap nodes
      (#[] : Array Nat).toList.toFinset := by
    simpa using (show PairVecDivVHCHeapChainsOwned heap nodes from
      ⟨owners, hownership⟩).away_empty heap nodes
  have hlin0 : PairVecDivVHCLinReady (#[] : Array Nat) nodes := by
    simp [PairVecDivVHCLinReady]
  rcases pairVecDivVHCConsumeEqualDegree_preserves_away_linReady this
      frontier.degree heap frontier.coefficient nodes #[] resetH quotient divisor
      consumed haway0 hlin0 hconsume with ⟨hconsumedAway, hconsumedLin⟩
  rcases hconsumedAway with
    ⟨consumedOwners, hconsumedOwnership, hconsumedSeparated⟩
  have hconsumedHomogeneous' := hconsumedHomogeneous.2.congr_owners
    consumed.heap owners consumedOwners consumed.nodes hconsumedHomogeneous.1
    hconsumedOwnership
  rcases pairVecDivVHCEmit_preserves_heapChainsHomogeneous this frontier
      consumed quotient divisor quotient' activated resetH' consumedOwners
      hconsumedOwnership hconsumedHomogeneous' hconsumedNodes.2.2.2 hdivisor
      hemit with ⟨emitOwners, hemitOwnership, hemitHomogeneous⟩
  rcases pairVecDivVHCEmit_preserves_away_linReady this frontier consumed
      quotient divisor quotient' activated resetH'
      ⟨consumedOwners, hconsumedOwnership, hconsumedSeparated⟩ hconsumedLin
      hconsumedNodes.2.2.2 hdivisor hemit with
    ⟨hemitAway, hemitLin⟩
  rcases hemitAway with ⟨reinsertOwners, hreinsertOwnership,
    hreinsertSeparated⟩
  have hemitHomogeneous' := hemitHomogeneous.congr_owners activated.heap
    emitOwners reinsertOwners activated.nodes hemitOwnership hreinsertOwnership
  have hfinal := pairVecDivVHCReinsertLin_preserves_heapChainsHomogeneous
    activated.heap activated.nodes consumed.lin reinserted reinsertOwners
    hreinsertOwnership hreinsertSeparated hemitHomogeneous' hemitLin hreinsert
  rw [hresult]
  exact hfinal

theorem pairVecDivVHCEmit_preserves_stateCovered
    (this : DenseUPolyZp) (frontier : PairVecDivVHCFrontier)
    (consumed : PairVecDivVHCEqualDegreeResult)
    (quotient divisor quotient' : SparsePolyZp)
    (activated : PairVecDivVHCHeapState) (resetH' : Nat)
    (haway : PairVecDivVHCHeapChainsOwnedAway consumed.heap consumed.nodes
      consumed.lin.toList.toFinset)
    (hlinReady : PairVecDivVHCLinReady consumed.lin consumed.nodes)
    (hready : PairVecDivVHCResetReady consumed.resetH quotient.size
      consumed.nodes)
    (hcovered : PairVecDivVHCStateCovered consumed.heap consumed.nodes
      consumed.lin consumed.resetH)
    (hdivisor : 0 < divisor.size)
    (hrun : pairVecDivVHCEmit this frontier consumed quotient divisor hdivisor =
      .ok (quotient', activated, resetH')) :
    PairVecDivVHCStateCovered activated.heap activated.nodes consumed.lin
      resetH' := by
  unfold pairVecDivVHCEmit at hrun
  by_cases hcoefficient : consumed.coefficient ≠ 0
  · rw [if_pos hcoefficient] at hrun
    by_cases hdegree : divisor[0].1.deg ≤ frontier.degree
    · rw [if_pos hdegree] at hrun
      let inverse := Generated.StrictGCD.dense_upoly_zp_nmod_inv_ir this
        divisor[0].2.val
      let value := Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this
        consumed.coefficient inverse
      by_cases hvalue : value ≠ 0
      · rw [if_pos hvalue] at hrun
        let emitted := quotient.push
          (⟨frontier.degree - divisor[0].1.deg⟩, ⟨value, this._p⟩)
        cases hactivate : pairVecDivVHCActivateReset consumed.resetH
            consumed.heap consumed.nodes emitted divisor with
        | error fault =>
            dsimp only [emitted] at hactivate
            dsimp only at hrun
            rw [hactivate] at hrun
            contradiction
        | ok state =>
            dsimp only [emitted] at hactivate
            dsimp only at hrun
            rw [hactivate] at hrun
            simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
            rcases hrun with ⟨rfl, rfl, rfl⟩
            exact pairVecDivVHCActivateReset_preserves_stateCovered
              consumed.resetH quotient.size consumed.heap consumed.nodes
              consumed.lin emitted divisor state haway hlinReady hready
              hcovered hactivate
      · rw [if_neg hvalue] at hrun
        simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
        rcases hrun with ⟨rfl, rfl, rfl⟩
        exact hcovered
    · rw [if_neg hdegree] at hrun
      simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
      rcases hrun with ⟨rfl, rfl, rfl⟩
      exact hcovered
  · rw [if_neg hcoefficient] at hrun
    simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
    rcases hrun with ⟨rfl, rfl, rfl⟩
    exact hcovered

theorem pairVecDivVHCEmit_preserves_cursorPrefixAbove
    (this : DenseUPolyZp) (frontier : PairVecDivVHCFrontier)
    (consumed : PairVecDivVHCEqualDegreeResult)
    (quotient divisor quotient' : SparsePolyZp)
    (activated : PairVecDivVHCHeapState) (resetH' : Nat)
    (hprefix : PairVecDivVHCCursorPrefixAbove frontier.degree consumed.nodes
      quotient divisor)
    (hbounded : PairVecDivVHCCursorIndicesBounded quotient.size
      consumed.nodes)
    (hdivisor : 0 < divisor.size)
    (hrun : pairVecDivVHCEmit this frontier consumed quotient divisor hdivisor =
      .ok (quotient', activated, resetH')) :
    PairVecDivVHCCursorPrefixAbove frontier.degree activated.nodes quotient'
      divisor := by
  unfold pairVecDivVHCEmit at hrun
  by_cases hcoefficient : consumed.coefficient ≠ 0
  · rw [if_pos hcoefficient] at hrun
    by_cases hdegree : divisor[0].1.deg ≤ frontier.degree
    · rw [if_pos hdegree] at hrun
      let inverse := Generated.StrictGCD.dense_upoly_zp_nmod_inv_ir this
        divisor[0].2.val
      let value := Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this
        consumed.coefficient inverse
      by_cases hvalue : value ≠ 0
      · rw [if_pos hvalue] at hrun
        let emitted := quotient.push
          (⟨frontier.degree - divisor[0].1.deg⟩, ⟨value, this._p⟩)
        cases hactivate : pairVecDivVHCActivateReset consumed.resetH
            consumed.heap consumed.nodes emitted divisor with
        | error fault =>
            dsimp only [emitted] at hactivate
            dsimp only at hrun
            rw [hactivate] at hrun
            contradiction
        | ok state =>
            dsimp only [emitted] at hactivate
            dsimp only at hrun
            rw [hactivate] at hrun
            simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
            rcases hrun with ⟨rfl, rfl, rfl⟩
            have hprefixEmitted := hprefix.push frontier.degree consumed.nodes
              quotient divisor
              (⟨frontier.degree - divisor[0].1.deg⟩, ⟨value, this._p⟩)
              hbounded
            exact pairVecDivVHCActivateReset_preserves_cursorPrefixAbove
              frontier.degree consumed.resetH consumed.heap consumed.nodes
              emitted divisor state hprefixEmitted hactivate
      · rw [if_neg hvalue] at hrun
        simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
        rcases hrun with ⟨rfl, rfl, rfl⟩
        exact hprefix
    · rw [if_neg hdegree] at hrun
      simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
      rcases hrun with ⟨rfl, rfl, rfl⟩
      exact hprefix
  · rw [if_neg hcoefficient] at hrun
    simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
    rcases hrun with ⟨rfl, rfl, rfl⟩
    exact hprefix

theorem pairVecDivVHCOuterIteration_preserves_stateCovered
    (this : DenseUPolyZp) (p degreeLimit dividendIndex : Nat)
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode)
    (quotient dividend divisor : SparsePolyZp) (resetH : Nat)
    (frontier : PairVecDivVHCFrontier)
    (result : PairVecDivVHCIterationResult) (hdivisor : 0 < divisor.size)
    (hselect : pairVecDivVHCSelectFrontier dividendIndex dividend heap nodes =
      .ok frontier)
    (hcanonical : SparsePolyZp.Canonical p quotient)
    (hbelow : PairVecDivVHCAllActiveNodesBelow degreeLimit nodes)
    (hdenotes : ∀ (i : Nat) (node : PairVecDivVHCNode),
      nodes[i]? = some node → node.mono ≠ none →
        PairVecDivVHCNodeDenotes quotient divisor node)
    (hfixed : PairVecDivVHCNodeDivisorIndicesFixed nodes)
    (hready : PairVecDivVHCResetReady resetH quotient.size nodes)
    (hcovered : PairVecDivVHCStateCovered heap nodes #[] resetH)
    (hrun : pairVecDivVHCOuterIteration this dividendIndex heap nodes quotient
      dividend divisor resetH = .ok result) :
    PairVecDivVHCStateCovered result.heap result.nodes #[] result.resetH := by
  rcases hcovered with ⟨owners, hownership, hnodesCovered⟩
  rcases pairVecDivVHCOuterIteration_components this dividendIndex heap nodes
      quotient dividend divisor resetH frontier result hdivisor hselect hrun with
    ⟨consumed, quotient', activated, resetH', reinserted, hconsume, hemit,
      hreinsert, hresult⟩
  have hstate0 : PairVecDivVHCStateCovered heap nodes #[] resetH :=
    ⟨owners, hownership, hnodesCovered⟩
  have hconsumedState :=
    pairVecDivVHCConsumeEqualDegree_preserves_stateCovered this
      frontier.degree heap frontier.coefficient nodes #[] resetH quotient
      divisor consumed hstate0 hconsume
  have howned : PairVecDivVHCHeapChainsOwned heap nodes :=
    ⟨owners, hownership⟩
  have haway0 : PairVecDivVHCHeapChainsOwnedAway heap nodes
      (#[] : Array Nat).toList.toFinset := by
    simpa using howned.away_empty
  have hlin0 : PairVecDivVHCLinReady (#[] : Array Nat) nodes := by
    simp [PairVecDivVHCLinReady]
  have hconsumedAway :=
    pairVecDivVHCConsumeEqualDegree_preserves_away_linReady this
      frontier.degree heap frontier.coefficient nodes #[] resetH quotient
      divisor consumed haway0 hlin0 hconsume
  have hconsumedInvariants :=
    pairVecDivVHCConsumeEqualDegree_preserves_node_invariants this p
      degreeLimit frontier.degree heap frontier.coefficient nodes #[] resetH
      quotient divisor consumed owners hcanonical hbelow hdenotes hfixed
      hready hownership hconsume
  have hemittedState := pairVecDivVHCEmit_preserves_stateCovered this frontier
    consumed quotient divisor quotient' activated resetH' hconsumedAway.1
    hconsumedAway.2 hconsumedInvariants.2.2.2 hconsumedState hdivisor hemit
  have hemittedAway := pairVecDivVHCEmit_preserves_away_linReady this frontier
    consumed quotient divisor quotient' activated resetH' hconsumedAway.1
    hconsumedAway.2 hconsumedInvariants.2.2.2 hdivisor hemit
  have hreinserted := pairVecDivVHCReinsertLin_preserves_stateCovered
    activated.heap activated.nodes consumed.lin resetH' reinserted
    hemittedAway.1 hemittedAway.2 hemittedState hreinsert
  subst result
  exact hreinserted

theorem pairVecDivVHCEmit_preserves_node_invariants
    (this : DenseUPolyZp) (p degreeLimit : Nat)
    (frontier : PairVecDivVHCFrontier)
    (consumed : PairVecDivVHCEqualDegreeResult)
    (quotient divisor quotient' : SparsePolyZp)
    (activated : PairVecDivVHCHeapState) (resetH' : Nat)
    (hdivisorCanonical : SparsePolyZp.Canonical p divisor)
    (hdivisor : 0 < divisor.size)
    (hfrontierBelow : frontier.degree ≤ degreeLimit)
    (hbelow : PairVecDivVHCAllActiveNodesBelow degreeLimit consumed.nodes)
    (hdenotes : ∀ (i : Nat) (node : PairVecDivVHCNode),
      consumed.nodes[i]? = some node → node.mono ≠ none →
        PairVecDivVHCNodeDenotes quotient divisor node)
    (hfixed : PairVecDivVHCNodeDivisorIndicesFixed consumed.nodes)
    (hready : PairVecDivVHCResetReady consumed.resetH quotient.size
      consumed.nodes)
    (hrun : pairVecDivVHCEmit this frontier consumed quotient divisor hdivisor =
      .ok (quotient', activated, resetH')) :
    PairVecDivVHCAllActiveNodesBelow degreeLimit activated.nodes ∧
      (∀ (i : Nat) (node : PairVecDivVHCNode),
        activated.nodes[i]? = some node → node.mono ≠ none →
          PairVecDivVHCNodeDenotes quotient' divisor node) ∧
      PairVecDivVHCNodeDivisorIndicesFixed activated.nodes ∧
      PairVecDivVHCResetReady resetH' quotient'.size activated.nodes := by
  unfold pairVecDivVHCEmit at hrun
  by_cases hcoefficient : consumed.coefficient ≠ 0
  · rw [if_pos hcoefficient] at hrun
    by_cases hdegree : divisor[0].1.deg ≤ frontier.degree
    · rw [if_pos hdegree] at hrun
      let inverse := Generated.StrictGCD.dense_upoly_zp_nmod_inv_ir this
        divisor[0].2.val
      let value := Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this
        consumed.coefficient inverse
      by_cases hvalue : value ≠ 0
      · rw [if_pos hvalue] at hrun
        let emitted := quotient.push
          (⟨frontier.degree - divisor[0].1.deg⟩, ⟨value, this._p⟩)
        cases hactivate : pairVecDivVHCActivateReset consumed.resetH
            consumed.heap consumed.nodes emitted divisor with
        | error fault =>
            dsimp only [emitted] at hactivate
            dsimp only at hrun
            rw [hactivate] at hrun
            contradiction
        | ok state =>
            dsimp only [emitted] at hactivate
            dsimp only at hrun
            rw [hactivate] at hrun
            simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
            rcases hrun with ⟨rfl, rfl, rfl⟩
            have hdenotesEmitted : ∀ (i : Nat)
                (node : PairVecDivVHCNode),
                consumed.nodes[i]? = some node → node.mono ≠ none →
                  PairVecDivVHCNodeDenotes emitted divisor node := by
              intro i node hnode hmono
              rcases hdenotes i node hnode hmono with
                ⟨quotientTerm, divisorTerm, hquotientTerm, hdivisorTerm,
                  hnodeMono⟩
              refine ⟨quotientTerm, divisorTerm, ?_, hdivisorTerm,
                hnodeMono⟩
              rw [Array.getElem?_push]
              have hne : node.quotientIndex ≠ quotient.size := by
                intro heq
                rw [heq, Array.getElem?_eq_none (Nat.le_refl _)] at hquotientTerm
                contradiction
              simp only [hne, ↓reduceIte, hquotientTerm]
            have hinvariants :=
              pairVecDivVHCActivateReset_preserves_node_invariants p
                frontier.degree degreeLimit consumed.resetH quotient.size
                consumed.heap consumed.nodes emitted divisor state
                hdivisorCanonical hdivisor hdegree hfrontierBelow (by simp
                  [emitted]) (by simp [emitted]) hbelow hdenotesEmitted hfixed
                hready hactivate
            refine ⟨hinvariants.1, hinvariants.2.1, hinvariants.2.2.1, ?_⟩
            exact ⟨Nat.zero_le _, fun i hi => (Nat.not_lt_zero i hi).elim⟩
      · rw [if_neg hvalue] at hrun
        simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
        rcases hrun with ⟨rfl, rfl, rfl⟩
        exact ⟨hbelow, hdenotes, hfixed, hready⟩
    · rw [if_neg hdegree] at hrun
      simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
      rcases hrun with ⟨rfl, rfl, rfl⟩
      exact ⟨hbelow, hdenotes, hfixed, hready⟩
  · rw [if_neg hcoefficient] at hrun
    simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
    rcases hrun with ⟨rfl, rfl, rfl⟩
    exact ⟨hbelow, hdenotes, hfixed, hready⟩

theorem pairVecDivVHCOuterIteration_preserves_cursorPrefixAbove
    (this : DenseUPolyZp) (p degreeLimit dividendIndex : Nat)
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode)
    (quotient dividend divisor : SparsePolyZp) (resetH : Nat)
    (frontier : PairVecDivVHCFrontier)
    (result : PairVecDivVHCIterationResult) (owners : Nat → Finset Nat)
    (hdivisor : 0 < divisor.size)
    (hselect : pairVecDivVHCSelectFrontier dividendIndex dividend heap nodes =
      .ok frontier)
    (hdecrease : frontier.degree < degreeLimit)
    (hcanonical : SparsePolyZp.Canonical p quotient)
    (hbelow : PairVecDivVHCAllActiveNodesBelow degreeLimit nodes)
    (hdenotes : ∀ (i : Nat) (node : PairVecDivVHCNode),
      nodes[i]? = some node → node.mono ≠ none →
        PairVecDivVHCNodeDenotes quotient divisor node)
    (hfixed : PairVecDivVHCNodeDivisorIndicesFixed nodes)
    (hready : PairVecDivVHCResetReady resetH quotient.size nodes)
    (hstate : PairVecDivVHCStateCovered heap nodes #[] resetH)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hhomogeneous : PairVecDivVHCHeapChainsHomogeneous heap owners nodes)
    (hprefix : PairVecDivVHCCursorPrefixAbove degreeLimit nodes quotient divisor)
    (hrun : pairVecDivVHCOuterIteration this dividendIndex heap nodes quotient
      dividend divisor resetH = .ok result) :
    PairVecDivVHCCursorPrefixAbove frontier.degree result.nodes
      result.quotient divisor := by
  rcases pairVecDivVHCOuterIteration_components this dividendIndex heap nodes
      quotient dividend divisor resetH frontier result hdivisor hselect hrun with
    ⟨consumed, quotient', activated, resetH', reinserted, hconsume, hemit,
      hreinsert, hresult⟩
  have hprefixCurrent := hprefix.mono degreeLimit frontier.degree nodes
    quotient divisor (Nat.le_of_lt hdecrease)
  have hconsumedPrefix :=
    pairVecDivVHCConsumeEqualDegree_preserves_cursorPrefixAbove this
      frontier.degree heap frontier.coefficient nodes #[] resetH quotient
      divisor consumed owners hownership hhomogeneous hprefixCurrent hdenotes
      hconsume
  have hconsumedState :=
    pairVecDivVHCConsumeEqualDegree_preserves_stateCovered this
      frontier.degree heap frontier.coefficient nodes #[] resetH quotient
      divisor consumed hstate hconsume
  have howned : PairVecDivVHCHeapChainsOwned heap nodes :=
    ⟨owners, hownership⟩
  have haway0 : PairVecDivVHCHeapChainsOwnedAway heap nodes
      (#[] : Array Nat).toList.toFinset := by
    simpa using howned.away_empty
  have hlin0 : PairVecDivVHCLinReady (#[] : Array Nat) nodes := by
    simp [PairVecDivVHCLinReady]
  have hconsumedAway :=
    pairVecDivVHCConsumeEqualDegree_preserves_away_linReady this
      frontier.degree heap frontier.coefficient nodes #[] resetH quotient
      divisor consumed haway0 hlin0 hconsume
  have hconsumedInvariants :=
    pairVecDivVHCConsumeEqualDegree_preserves_node_invariants this p
      degreeLimit frontier.degree heap frontier.coefficient nodes #[] resetH
      quotient divisor consumed owners hcanonical hbelow hdenotes hfixed
      hready hownership hconsume
  have hbounded := pairVecDivVHCCursorIndicesBounded_of_state consumed.heap
    consumed.nodes consumed.lin consumed.resetH quotient divisor hconsumedState
    hconsumedAway.2 hconsumedInvariants.2.2.2 hconsumedInvariants.2.1
  have hemittedPrefix := pairVecDivVHCEmit_preserves_cursorPrefixAbove this
    frontier consumed quotient divisor quotient' activated resetH'
    hconsumedPrefix hbounded hdivisor hemit
  have hreinsertedPrefix :=
    pairVecDivVHCReinsertLin_preserves_cursorPrefixAbove frontier.degree
      activated.heap activated.nodes consumed.lin quotient' divisor reinserted
      hemittedPrefix hreinsert
  subst result
  exact hreinsertedPrefix

theorem pairVecDivVHCEmit_preserves_canonical
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (frontier : PairVecDivVHCFrontier)
    (consumed : PairVecDivVHCEqualDegreeResult)
    (quotient divisor quotient' : SparsePolyZp)
    (activated : PairVecDivVHCHeapState) (resetH' : Nat)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hcanonical : SparsePolyZp.Canonical this._p.toNat quotient)
    (hcoefficientReduced : consumed.coefficient.toNat < this._p.toNat)
    (hdivisor : 0 < divisor.size)
    (hdegrees : PairVecDivVHCQuotientAbove frontier.degree divisor[0].1.deg
      quotient)
    (hrun : pairVecDivVHCEmit this frontier consumed quotient divisor hdivisor =
      .ok (quotient', activated, resetH')) :
    SparsePolyZp.Canonical this._p.toNat quotient' := by
  unfold pairVecDivVHCEmit at hrun
  by_cases hcoefficient : consumed.coefficient ≠ 0
  · rw [if_pos hcoefficient] at hrun
    by_cases hdegree : divisor[0].1.deg ≤ frontier.degree
    · rw [if_pos hdegree] at hrun
      let inverse := Generated.StrictGCD.dense_upoly_zp_nmod_inv_ir this
        divisor[0].2.val
      let value := Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this
        consumed.coefficient inverse
      have hvalueNat :=
        CLPoly.Impl.StrictWordArithmetic.nmod_mul_ir_correct_of_configured this
          consumed.coefficient inverse hcfg hcoefficientReduced
      change value.toNat = _ at hvalueNat
      have hvalueReduced : value.toNat < this._p.toNat := by
        rw [hvalueNat]
        exact Nat.mod_lt _ (Fact.out : Nat.Prime this._p.toNat).pos
      by_cases hvalue : value ≠ 0
      · rw [if_pos hvalue] at hrun
        let emitted := quotient.push
          (⟨frontier.degree - divisor[0].1.deg⟩, ⟨value, this._p⟩)
        cases hactivate : pairVecDivVHCActivateReset consumed.resetH
            consumed.heap consumed.nodes emitted divisor with
        | error fault =>
            dsimp only [emitted] at hactivate
            dsimp only at hrun
            rw [hactivate] at hrun
            contradiction
        | ok state =>
            dsimp only [emitted] at hactivate
            dsimp only at hrun
            rw [hactivate] at hrun
            simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
            rcases hrun with ⟨rfl, rfl, rfl⟩
            exact
              CLPoly.Impl.StrictPolynomialGCDRefinement.canonical_push_lower
                this._p quotient (frontier.degree - divisor[0].1.deg) value
                hcanonical (hdegrees hdegree) hvalueReduced hvalue
      · rw [if_neg hvalue] at hrun
        simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
        rcases hrun with ⟨rfl, rfl, rfl⟩
        exact hcanonical
    · rw [if_neg hdegree] at hrun
      simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
      rcases hrun with ⟨rfl, rfl, rfl⟩
      exact hcanonical
  · rw [if_neg hcoefficient] at hrun
    simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
    rcases hrun with ⟨rfl, rfl, rfl⟩
    exact hcanonical

/-- Exact polynomial delta of the generated emit branch whenever the divisor
lead monomial divides the selected frontier monomial.  The coefficient is
computed by the same generated modular inverse/multiply operations used by
the source; `pairVecDivSingleTermIR_refines` supplies their field semantics. -/
theorem pairVecDivVHCEmit_toPoly_of_lead_le
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (frontier : PairVecDivVHCFrontier)
    (consumed : PairVecDivVHCEqualDegreeResult)
    (quotient divisor quotient' : SparsePolyZp)
    (activated : PairVecDivVHCHeapState) (resetH' : Nat)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hdivisorCanonical : SparsePolyZp.Canonical this._p.toNat divisor)
    (hcoefficientReduced : consumed.coefficient.toNat < this._p.toNat)
    (hdivisor : 0 < divisor.size)
    (hlead : divisor[0].1.deg ≤ frontier.degree)
    (hrun : pairVecDivVHCEmit this frontier consumed quotient divisor hdivisor =
      .ok (quotient', activated, resetH')) :
    SparsePolyZp.toPoly this._p.toNat quotient' =
      SparsePolyZp.toPoly this._p.toNat quotient +
        Polynomial.monomial (frontier.degree - divisor[0].1.deg)
          (Zp.toZMod this._p.toNat ⟨consumed.coefficient, this._p⟩ /
            Zp.toZMod this._p.toNat divisor[0].2) := by
  have hdivisorMem : divisor[0] ∈ divisor.toList :=
    Array.getElem_mem_toList hdivisor
  have hdivisorReduced := hdivisorCanonical.1 divisor[0] hdivisorMem
  have hdivisorNonzero := hdivisorCanonical.2.2 divisor[0] hdivisorMem
  unfold pairVecDivVHCEmit at hrun
  by_cases hcoefficient : consumed.coefficient ≠ 0
  · rw [if_pos hcoefficient, if_pos hlead] at hrun
    let sourceTerm : UMonomial × Zp :=
      (⟨frontier.degree⟩, ⟨consumed.coefficient, this._p⟩)
    let inverse := Generated.StrictGCD.dense_upoly_zp_nmod_inv_ir this
      divisor[0].2.val
    let value := Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this
      consumed.coefficient inverse
    let emittedTerm : UMonomial × Zp :=
      (⟨frontier.degree - divisor[0].1.deg⟩, ⟨value, this._p⟩)
    have hsourceReduced : Zp.Reduced this._p.toNat sourceTerm.2 := by
      exact ⟨rfl, hcoefficientReduced⟩
    have hsourceNonzero : sourceTerm.2.val ≠ 0 := hcoefficient
    have hsingle : pairVecDivSingleTermIR this divisor[0] sourceTerm =
        some emittedTerm := by
      simp [pairVecDivSingleTermIR, sourceTerm, emittedTerm, inverse, value,
        hlead]
    have hrefines := pairVecDivSingleTermIR_refines this divisor[0] sourceTerm
      emittedTerm hcfg hdivisorReduced hsourceReduced hdivisorNonzero
      hsourceNonzero hsingle
    have hvalue : value ≠ 0 := by
      exact hrefines.2.2.1
    have hvalueField : (value.toNat : ZMod this._p.toNat) =
        Zp.toZMod this._p.toNat ⟨consumed.coefficient, this._p⟩ /
          Zp.toZMod this._p.toNat divisor[0].2 := by
      simpa [Zp.toZMod, emittedTerm, sourceTerm] using hrefines.2.2.2
    rw [if_pos hvalue] at hrun
    let emitted := quotient.push emittedTerm
    cases hactivate : pairVecDivVHCActivateReset consumed.resetH consumed.heap
        consumed.nodes emitted divisor with
    | error fault =>
        dsimp only [emitted, emittedTerm] at hactivate
        dsimp only at hrun
        rw [hactivate] at hrun
        contradiction
    | ok state =>
        dsimp only [emitted, emittedTerm] at hactivate
        dsimp only at hrun
        rw [hactivate] at hrun
        simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
        rcases hrun with ⟨rfl, rfl, rfl⟩
        rw [CLPoly.Impl.StrictPolynomialGCDRefinement.sparseToPoly_push_raw]
        rw [hvalueField]
  · rw [if_neg hcoefficient] at hrun
    simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
    rcases hrun with ⟨rfl, rfl, rfl⟩
    have hzero : consumed.coefficient = 0 := by
      exact Decidable.not_not.mp hcoefficient
    rw [hzero]
    simp [Zp.toZMod]

theorem sparsePolyZp_toPoly_eq_head_add_tail (p : Nat)
    (divisor : SparsePolyZp) (hdivisor : 0 < divisor.size) :
    SparsePolyZp.toPoly p divisor =
      Polynomial.monomial divisor[0].1.deg
          (Zp.toZMod p divisor[0].2) +
        listSum p divisor.toList.tail := by
  have hnonempty : divisor.toList ≠ [] := by
    intro hempty
    have : divisor.toList.length = 0 := by rw [hempty]; rfl
    simp only [Array.length_toList] at this
    omega
  have hlist : divisor.toList = divisor[0] :: divisor.toList.tail := by
    have hdrop := List.drop_eq_getElem_cons (l := divisor.toList) (i := 0)
      (by simpa using hdivisor)
    simpa [Array.getElem_toList hdivisor] using hdrop
  unfold SparsePolyZp.toPoly
  rw [hlist, listSum_cons]
  rfl

theorem pairVecDivVHCQuotient_mul_lead_coeff_eq_zero
    (p frontierDegree leadDegree : Nat) (quotient : SparsePolyZp)
    (leadCoefficient : Zp)
    (habove : PairVecDivVHCQuotientAbove frontierDegree leadDegree quotient)
    (hlead : leadDegree ≤ frontierDegree) :
    (SparsePolyZp.toPoly p quotient *
      Polynomial.monomial leadDegree (Zp.toZMod p leadCoefficient)).coeff
        frontierDegree = 0 := by
  unfold SparsePolyZp.toPoly
  rw [show Polynomial.monomial leadDegree (Zp.toZMod p leadCoefficient) =
      listSum p [(⟨leadDegree⟩, leadCoefficient)] by
    simp [listSum, Zp.toZMod]]
  rw [← pairVecDivVHCListProductCoeffValue_eq_coeff p frontierDegree
    quotient.toList [(⟨leadDegree⟩, leadCoefficient)]]
  have hdegrees := habove hlead
  have go : ∀ terms : List (UMonomial × Zp),
      (∀ term ∈ terms, frontierDegree - leadDegree < term.1.deg) →
      pairVecDivVHCListProductCoeffValue p frontierDegree terms
        [(⟨leadDegree⟩, leadCoefficient)] = 0 := by
    intro terms hterms
    induction terms with
    | nil => simp [pairVecDivVHCListProductCoeffValue]
    | cons term terms ih =>
        have htermAbove := hterms term List.mem_cons_self
        have hrestAbove : ∀ item ∈ terms,
            frontierDegree - leadDegree < item.1.deg := by
          intro item hitem
          exact hterms item (List.mem_cons_of_mem term hitem)
        have hdegree : term.1.deg + leadDegree ≠ frontierDegree := by omega
        simp [pairVecDivVHCListProductCoeffValue, hdegree, ih hrestAbove]
  exact go quotient.toList hdegrees

theorem pairVecDivVHCQuotient_mul_lead_coeff_eq_zero_below_processed
    (p degreeLimit targetDegree leadDegree : Nat)
    (quotient : SparsePolyZp) (leadCoefficient : Zp)
    (hprocessed : PairVecDivVHCQuotientLeadAbove degreeLimit leadDegree
      quotient)
    (htarget : targetDegree < degreeLimit) :
    (SparsePolyZp.toPoly p quotient *
      Polynomial.monomial leadDegree (Zp.toZMod p leadCoefficient)).coeff
        targetDegree = 0 := by
  unfold SparsePolyZp.toPoly
  rw [show Polynomial.monomial leadDegree (Zp.toZMod p leadCoefficient) =
      listSum p [(⟨leadDegree⟩, leadCoefficient)] by
    simp [listSum, Zp.toZMod]]
  rw [← pairVecDivVHCListProductCoeffValue_eq_coeff p targetDegree
    quotient.toList [(⟨leadDegree⟩, leadCoefficient)]]
  have go : ∀ terms : List (UMonomial × Zp),
      (∀ term ∈ terms, degreeLimit ≤ term.1.deg + leadDegree) →
      pairVecDivVHCListProductCoeffValue p targetDegree terms
        [(⟨leadDegree⟩, leadCoefficient)] = 0 := by
    intro terms hterms
    induction terms with
    | nil => simp [pairVecDivVHCListProductCoeffValue]
    | cons term terms ih =>
        have htermAbove := hterms term List.mem_cons_self
        have hrestAbove : ∀ item ∈ terms,
            degreeLimit ≤ item.1.deg + leadDegree := by
          intro item hitem
          exact hterms item (List.mem_cons_of_mem term hitem)
        have hdegree : term.1.deg + leadDegree ≠ targetDegree := by omega
        simp [pairVecDivVHCListProductCoeffValue, hdegree, ih hrestAbove]
  exact go quotient.toList hprocessed

theorem pairVecDivVHCProduct_coeff_eq_zero_of_gap
    (p degreeLimit targetDegree dividendIndex resetH : Nat)
    (dividend quotient divisor : SparsePolyZp)
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode)
    (owners : Nat → Finset Nat) (frontier : PairVecDivVHCFrontier)
    (hdivisor : 0 < divisor.size)
    (hsize : nodes.size = divisor.size - 1)
    (hfixed : PairVecDivVHCNodeDivisorIndicesFixed nodes)
    (hstate : PairVecDivVHCStateCovered heap nodes #[] resetH)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hhomogeneous : PairVecDivVHCHeapChainsHomogeneous heap owners nodes)
    (hresetReady : PairVecDivVHCResetReady resetH quotient.size nodes)
    (hordered : PairVecDivVHCHeapOrdered heap nodes)
    (hdenotes : ∀ (i : Nat) (node : PairVecDivVHCNode),
      nodes[i]? = some node → node.mono ≠ none →
        PairVecDivVHCNodeDenotes quotient divisor node)
    (hcanonical : SparsePolyZp.Canonical p quotient)
    (hprefix : PairVecDivVHCCursorPrefixAbove degreeLimit nodes quotient divisor)
    (hprocessed : PairVecDivVHCQuotientLeadAbove degreeLimit
      divisor[0].1.deg quotient)
    (hfrontier : frontier.degree < targetDegree)
    (htarget : targetDegree < degreeLimit)
    (hselect : pairVecDivVHCSelectFrontier dividendIndex dividend heap nodes =
      .ok frontier) :
    (SparsePolyZp.toPoly p quotient *
      SparsePolyZp.toPoly p divisor).coeff targetDegree = 0 := by
  rw [sparsePolyZp_toPoly_eq_head_add_tail p divisor hdivisor, mul_add,
    Polynomial.coeff_add,
    pairVecDivVHCQuotient_mul_lead_coeff_eq_zero_below_processed p
      degreeLimit targetDegree divisor[0].1.deg quotient divisor[0].2
      hprocessed htarget,
    pairVecDivVHCTail_product_coeff_eq_zero_of_gap p degreeLimit targetDegree
      dividendIndex resetH dividend quotient divisor heap nodes owners frontier
      hsize hfixed hstate hownership hhomogeneous hresetReady hordered hdenotes
      hcanonical hprefix hfrontier htarget hselect]
  simp

/-- When the generated outer loop has exhausted both source queues, every
dividend coefficient below the current well-founded bound is zero.  This is
read directly from the consumed-prefix invariant and the concrete source
index; no divisibility specification is used. -/
theorem pairVecDivVHCDividend_coeff_eq_zero_of_done
    (p degreeLimit targetDegree dividendIndex : Nat)
    (dividend : SparsePolyZp)
    (hcanonical : SparsePolyZp.Canonical p dividend)
    (hconsumed : PairVecDivVHCConsumedDividendAbove degreeLimit dividendIndex
      dividend)
    (hdone : dividend.size ≤ dividendIndex)
    (htarget : targetDegree < degreeLimit) :
    (SparsePolyZp.toPoly p dividend).coeff targetDegree = 0 := by
  unfold SparsePolyZp.toPoly
  suffices habsent : ∀ term ∈ dividend.toList,
      term.1.deg ≠ targetDegree by
    have go : ∀ terms : List (UMonomial × Zp),
        (∀ term ∈ terms, term.1.deg ≠ targetDegree) →
        (listSum p terms).coeff targetDegree = 0 := by
      intro terms hterms
      induction terms with
      | nil => simp [listSum]
      | cons term rest ih =>
          rw [listSum_cons, Polynomial.coeff_add,
            Polynomial.coeff_monomial,
            if_neg (hterms term List.mem_cons_self),
            ih (by
              intro item hitem
              exact hterms item (List.mem_cons_of_mem term hitem))]
          simp
    exact go dividend.toList habsent
  intro term hterm
  rcases List.getElem_of_mem hterm with ⟨i, hiList, htermEq⟩
  have hi : i < dividend.size := by simpa using hiList
  have htermEq' : term = dividend[i] := by
    rw [← Array.getElem_toList hi]
    exact htermEq.symm
  have habove := hconsumed i dividend[i] (by omega)
    (Array.getElem?_eq_getElem hi)
  rw [htermEq']
  omega

/-- At a concrete terminal heap state every divisor-tail product has already
crossed the cursor.  Consequently no such product can occur below the current
outer-loop bound. -/
theorem pairVecDivVHCTail_product_coeff_eq_zero_of_done
    (p degreeLimit targetDegree resetH : Nat)
    (quotient divisor : SparsePolyZp) (nodes : Array PairVecDivVHCNode)
    (owners : Nat → Finset Nat)
    (hsize : nodes.size = divisor.size - 1)
    (hfixed : PairVecDivVHCNodeDivisorIndicesFixed nodes)
    (hstate : PairVecDivVHCStateCovered #[] nodes #[] resetH)
    (hownership : PairVecDivVHCHeapChainOwnership #[] owners nodes)
    (hresetReady : PairVecDivVHCResetReady resetH quotient.size nodes)
    (hprefix : PairVecDivVHCCursorPrefixAbove degreeLimit nodes quotient divisor)
    (htarget : targetDegree < degreeLimit) :
    (SparsePolyZp.toPoly p quotient *
      listSum p divisor.toList.tail).coeff targetDegree = 0 := by
  rw [← pairVecDivVHCTargetPairSum_eq_productCoeffTail p targetDegree
    quotient divisor]
  suffices hempty : PairVecDivVHCTargetPairsAtDegree targetDegree quotient
      divisor = ∅ by simp [hempty]
  rw [← Finset.not_nonempty_iff_eq_empty]
  rintro ⟨⟨q, d⟩, hpair⟩
  rw [PairVecDivVHCTargetPairsAtDegree] at hpair
  rcases Finset.mem_filter.mp hpair with ⟨hindices, hdegree⟩
  rcases Finset.mem_product.mp hindices with ⟨hq, hd⟩
  rw [Finset.mem_range] at hq
  rw [Finset.mem_Ico] at hd
  rcases hd with ⟨hdpos, hd⟩
  have hquotient : quotient[q]? = some quotient[q] :=
    Array.getElem?_eq_getElem hq
  have hdivisor : divisor[d]? = some divisor[d] :=
    Array.getElem?_eq_getElem hd
  simp [PairVecDivVHCPairAtDegree, hquotient, hdivisor] at hdegree
  rcases hfixed.node_for_tail nodes divisor.size d hsize hdpos hd with
    ⟨node, hnode, hnodeD⟩
  have hnodeIndex : d - 1 < nodes.size := by rw [hsize]; omega
  have hcovered := hstate.covered_with #[] nodes #[] resetH owners hownership
  rcases hcovered (d - 1) hnodeIndex with hreset | hlin | howned
  · rcases hresetReady.2 (d - 1) hreset with
      ⟨readyNode, hreadyNode, hcursor, hreadyD, hmono⟩
    rw [hnode] at hreadyNode
    simp only [Option.some.injEq] at hreadyNode
    subst readyNode
    have habove := hprefix (d - 1) node hnode q quotient[q] divisor[d]
      (by omega) hquotient (by simpa [hnodeD] using hdivisor)
    omega
  · simp at hlin
  · rcases howned with ⟨slot, head, hheap, _⟩
    simp at hheap

/-- The complete generated product has no coefficient below the current bound
once both concrete outer-loop queues are exhausted. -/
theorem pairVecDivVHCProduct_coeff_eq_zero_of_done
    (p degreeLimit targetDegree resetH : Nat)
    (quotient divisor : SparsePolyZp) (nodes : Array PairVecDivVHCNode)
    (owners : Nat → Finset Nat) (hdivisor : 0 < divisor.size)
    (hsize : nodes.size = divisor.size - 1)
    (hfixed : PairVecDivVHCNodeDivisorIndicesFixed nodes)
    (hstate : PairVecDivVHCStateCovered #[] nodes #[] resetH)
    (hownership : PairVecDivVHCHeapChainOwnership #[] owners nodes)
    (hresetReady : PairVecDivVHCResetReady resetH quotient.size nodes)
    (hprefix : PairVecDivVHCCursorPrefixAbove degreeLimit nodes quotient divisor)
    (hprocessed : PairVecDivVHCQuotientLeadAbove degreeLimit
      divisor[0].1.deg quotient)
    (htarget : targetDegree < degreeLimit) :
    (SparsePolyZp.toPoly p quotient * SparsePolyZp.toPoly p divisor).coeff
      targetDegree = 0 := by
  rw [sparsePolyZp_toPoly_eq_head_add_tail p divisor hdivisor, mul_add,
    Polynomial.coeff_add,
    pairVecDivVHCQuotient_mul_lead_coeff_eq_zero_below_processed p
      degreeLimit targetDegree divisor[0].1.deg quotient divisor[0].2
      hprocessed htarget,
    pairVecDivVHCTail_product_coeff_eq_zero_of_done p degreeLimit targetDegree
      resetH quotient divisor nodes owners hsize hfixed hstate hownership
      hresetReady hprefix htarget]
  simp

theorem pairVecDivVHCEmitted_monomial_mul_tail_coeff_eq_zero
    (p frontierDegree : Nat) (divisor : SparsePolyZp)
    (coefficient : ZMod p) (hdivisor : 0 < divisor.size)
    (hcanonical : SparsePolyZp.Canonical p divisor)
    (hlead : divisor[0].1.deg ≤ frontierDegree) :
    (Polynomial.monomial (frontierDegree - divisor[0].1.deg)
        coefficient *
      listSum p divisor.toList.tail).coeff frontierDegree = 0 := by
  have hlist : divisor.toList = divisor[0] :: divisor.toList.tail := by
    have hnonempty : divisor.toList ≠ [] := by
      intro hempty
      have : divisor.toList.length = 0 := by rw [hempty]; rfl
      simp only [Array.length_toList] at this
      omega
    have hdrop := List.drop_eq_getElem_cons (l := divisor.toList) (i := 0)
      (by simpa using hdivisor)
    simpa [Array.getElem_toList hdivisor] using hdrop
  have hchain : List.IsChain
      (fun left right : UMonomial × Zp => left.1.deg > right.1.deg)
      (divisor[0] :: divisor.toList.tail) := by
    rw [← hlist]
    exact hcanonical.2.1
  have htailDegree : ∀ term ∈ divisor.toList.tail,
      term.1.deg < divisor[0].1.deg := by
    intro term hterm
    exact chain_gt_all_after_head divisor[0] divisor.toList.tail hchain term
      hterm
  have go : ∀ terms : List (UMonomial × Zp),
      (∀ term ∈ terms, term.1.deg < divisor[0].1.deg) →
      (Polynomial.monomial (frontierDegree - divisor[0].1.deg) coefficient *
        listSum p terms).coeff frontierDegree = 0 := by
    intro terms hterms
    induction terms with
    | nil => simp [listSum]
    | cons term terms ih =>
        have htermBelow := hterms term List.mem_cons_self
        have hrestBelow : ∀ item ∈ terms,
            item.1.deg < divisor[0].1.deg := by
          intro item hitem
          exact hterms item (List.mem_cons_of_mem term hitem)
        have hdegree : frontierDegree - divisor[0].1.deg + term.1.deg ≠
            frontierDegree := by omega
        rw [listSum_cons, mul_add, Polynomial.coeff_add,
          Polynomial.monomial_mul_monomial, Polynomial.coeff_monomial,
          if_neg hdegree, ih hrestBelow]
        simp
  exact go divisor.toList.tail htailDegree

theorem pairVecDivVHCEmitted_monomial_mul_divisor_coeff_eq_zero_above
    (p frontierDegree targetDegree : Nat) (divisor : SparsePolyZp)
    (coefficient : ZMod p) (hdivisor : 0 < divisor.size)
    (hcanonical : SparsePolyZp.Canonical p divisor)
    (hlead : divisor[0].1.deg ≤ frontierDegree)
    (habove : frontierDegree < targetDegree) :
    (Polynomial.monomial (frontierDegree - divisor[0].1.deg) coefficient *
      SparsePolyZp.toPoly p divisor).coeff targetDegree = 0 := by
  have htermDegree : ∀ term ∈ divisor.toList,
      term.1.deg ≤ divisor[0].1.deg := by
    intro term hterm
    rcases List.getElem_of_mem hterm with ⟨i, hi, htermEq⟩
    have hiArray : i < divisor.size := by simpa using hi
    have hget : divisor[i]? = some term := by
      rw [Array.getElem?_eq_getElem hiArray]
      rw [← Array.getElem_toList hiArray]
      exact congrArg some htermEq
    exact canonical_degree_le_head p divisor hcanonical hdivisor i term hget
  unfold SparsePolyZp.toPoly
  have go : ∀ terms : List (UMonomial × Zp),
      (∀ term ∈ terms, term.1.deg ≤ divisor[0].1.deg) →
      (Polynomial.monomial (frontierDegree - divisor[0].1.deg) coefficient *
        listSum p terms).coeff targetDegree = 0 := by
    intro terms hterms
    induction terms with
    | nil => simp [listSum]
    | cons term terms ih =>
        have htermLe := hterms term List.mem_cons_self
        have hrest : ∀ item ∈ terms,
            item.1.deg ≤ divisor[0].1.deg := by
          intro item hitem
          exact hterms item (List.mem_cons_of_mem term hitem)
        have hdegree : frontierDegree - divisor[0].1.deg + term.1.deg ≠
            targetDegree := by omega
        rw [listSum_cons, mul_add, Polynomial.coeff_add,
          Polynomial.monomial_mul_monomial, Polynomial.coeff_monomial,
          if_neg hdegree, ih hrest]
        simp
  exact go divisor.toList htermDegree

theorem pairVecDivVHCEmit_product_coeff_above
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (frontier : PairVecDivVHCFrontier)
    (consumed : PairVecDivVHCEqualDegreeResult)
    (quotient divisor quotient' : SparsePolyZp)
    (activated : PairVecDivVHCHeapState) (resetH' targetDegree : Nat)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hdivisorCanonical : SparsePolyZp.Canonical this._p.toNat divisor)
    (hcoefficientReduced : consumed.coefficient.toNat < this._p.toNat)
    (hdivisor : 0 < divisor.size)
    (hlead : divisor[0].1.deg ≤ frontier.degree)
    (habove : frontier.degree < targetDegree)
    (hrun : pairVecDivVHCEmit this frontier consumed quotient divisor hdivisor =
      .ok (quotient', activated, resetH')) :
    (SparsePolyZp.toPoly this._p.toNat quotient' *
        SparsePolyZp.toPoly this._p.toNat divisor).coeff targetDegree =
      (SparsePolyZp.toPoly this._p.toNat quotient *
        SparsePolyZp.toPoly this._p.toNat divisor).coeff targetDegree := by
  rw [pairVecDivVHCEmit_toPoly_of_lead_le this frontier consumed quotient
    divisor quotient' activated resetH' hcfg hdivisorCanonical
    hcoefficientReduced hdivisor hlead hrun, add_mul,
    Polynomial.coeff_add,
    pairVecDivVHCEmitted_monomial_mul_divisor_coeff_eq_zero_above
      this._p.toNat frontier.degree targetDegree divisor
      (Zp.toZMod this._p.toNat ⟨consumed.coefficient, this._p⟩ /
        Zp.toZMod this._p.toNat divisor[0].2) hdivisor hdivisorCanonical
      hlead habove, add_zero]

theorem pairVecDivVHCEmit_preserves_quotientAbove_of_lt
    (this : DenseUPolyZp) (nextDegree : Nat)
    (frontier : PairVecDivVHCFrontier)
    (consumed : PairVecDivVHCEqualDegreeResult)
    (quotient divisor quotient' : SparsePolyZp)
    (activated : PairVecDivVHCHeapState) (resetH' : Nat)
    (hdivisor : 0 < divisor.size)
    (habove : PairVecDivVHCQuotientAbove frontier.degree divisor[0].1.deg
      quotient)
    (hnext : nextDegree < frontier.degree)
    (hrun : pairVecDivVHCEmit this frontier consumed quotient divisor hdivisor =
      .ok (quotient', activated, resetH')) :
    PairVecDivVHCQuotientAbove nextDegree divisor[0].1.deg quotient' := by
  have hold : PairVecDivVHCQuotientAbove nextDegree divisor[0].1.deg
      quotient := by
    intro hleadNext term hterm
    have hleadFrontier : divisor[0].1.deg ≤ frontier.degree :=
      Nat.le_trans hleadNext (Nat.le_of_lt hnext)
    have htermDegree := habove hleadFrontier term hterm
    omega
  unfold pairVecDivVHCEmit at hrun
  by_cases hcoefficient : consumed.coefficient ≠ 0
  · rw [if_pos hcoefficient] at hrun
    by_cases hdegree : divisor[0].1.deg ≤ frontier.degree
    · rw [if_pos hdegree] at hrun
      let inverse := Generated.StrictGCD.dense_upoly_zp_nmod_inv_ir this
        divisor[0].2.val
      let value := Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this
        consumed.coefficient inverse
      by_cases hvalue : value ≠ 0
      · rw [if_pos hvalue] at hrun
        let emitted := quotient.push
          (⟨frontier.degree - divisor[0].1.deg⟩, ⟨value, this._p⟩)
        cases hactivate : pairVecDivVHCActivateReset consumed.resetH
            consumed.heap consumed.nodes emitted divisor with
        | error fault =>
            dsimp only [emitted] at hactivate
            dsimp only at hrun
            rw [hactivate] at hrun
            contradiction
        | ok state =>
            dsimp only [emitted] at hactivate
            dsimp only at hrun
            rw [hactivate] at hrun
            simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
            rcases hrun with ⟨rfl, rfl, rfl⟩
            intro hleadNext term hterm
            simp only [Array.toList_push, List.mem_append,
              List.mem_singleton] at hterm
            rcases hterm with hterm | rfl
            · exact hold hleadNext term hterm
            · change nextDegree - divisor[0].1.deg <
                frontier.degree - divisor[0].1.deg
              omega
      · rw [if_neg hvalue] at hrun
        simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
        rcases hrun with ⟨rfl, rfl, rfl⟩
        exact hold
    · rw [if_neg hdegree] at hrun
      simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
      rcases hrun with ⟨rfl, rfl, rfl⟩
      exact hold
  · rw [if_neg hcoefficient] at hrun
    simp only [Except.ok.injEq, Prod.mk.injEq] at hrun
    rcases hrun with ⟨rfl, rfl, rfl⟩
    exact hold

theorem pairVecDivVHCOuterIteration_preserves_heapChainsOwned
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (degreeLimit dividendIndex : Nat)
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode)
    (quotient dividend divisor : SparsePolyZp) (resetH : Nat)
    (result : PairVecDivVHCIterationResult) (owners : Nat → Finset Nat)
    (hdivisor : 0 < divisor.size)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hcanonical : SparsePolyZp.Canonical this._p.toNat quotient)
    (hdividendCanonical : SparsePolyZp.Canonical this._p.toNat dividend)
    (hdivisorCanonical : SparsePolyZp.Canonical this._p.toNat divisor)
    (hquotientReady : ∀ frontier : PairVecDivVHCFrontier,
      pairVecDivVHCSelectFrontier dividendIndex dividend heap nodes =
        .ok frontier →
      PairVecDivVHCQuotientAbove frontier.degree divisor[0].1.deg quotient)
    (hremaining : PairVecDivVHCRemainingDividendBelow degreeLimit
      dividendIndex dividend)
    (hbelow : PairVecDivVHCAllActiveNodesBelow degreeLimit nodes)
    (hdenotes : ∀ (i : Nat) (node : PairVecDivVHCNode),
      nodes[i]? = some node → node.mono ≠ none →
        PairVecDivVHCNodeDenotes quotient divisor node)
    (hfixed : PairVecDivVHCNodeDivisorIndicesFixed nodes)
    (hready : PairVecDivVHCResetReady resetH quotient.size nodes)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hrun : pairVecDivVHCOuterIteration this dividendIndex heap nodes quotient
      dividend divisor resetH = .ok result) :
    SparsePolyZp.Canonical this._p.toNat result.quotient ∧
      PairVecDivVHCHeapChainsOwned result.heap result.nodes ∧
      PairVecDivVHCRemainingDividendBelow degreeLimit result.dividendIndex
        dividend ∧
      PairVecDivVHCAllActiveNodesBelow degreeLimit result.nodes ∧
      (∀ (i : Nat) (node : PairVecDivVHCNode),
        result.nodes[i]? = some node → node.mono ≠ none →
          PairVecDivVHCNodeDenotes result.quotient divisor node) ∧
      PairVecDivVHCNodeDivisorIndicesFixed result.nodes ∧
      PairVecDivVHCResetReady result.resetH result.quotient.size result.nodes := by
  unfold pairVecDivVHCOuterIteration at hrun
  simp only [hdivisor, ↓reduceDIte] at hrun
  generalize hselect :
      pairVecDivVHCSelectFrontier dividendIndex dividend heap nodes = selected
      at hrun
  cases selected with
  | error fault =>
      simp only [Bind.bind, Except.bind] at hrun
      contradiction
  | ok frontier =>
      simp only [Bind.bind, Except.bind] at hrun
      have hfrontierBelow := pairVecDivVHCFrontierBelow_of_remaining_owned
        degreeLimit dividendIndex dividend heap nodes owners hremaining hbelow
        hownership
      have hselectedBelow := pairVecDivVHCSelectFrontier_degree_lt degreeLimit
        dividendIndex dividend heap nodes frontier hfrontierBelow hselect
      have hremaining' :=
        pairVecDivVHCSelectFrontier_preserves_remaining_below degreeLimit
          dividendIndex dividend heap nodes frontier hremaining hselect
      have hfrontierCoefficient :=
        pairVecDivVHCSelectFrontier_coefficient_reduced this._p.toNat
          dividendIndex dividend heap nodes frontier
          (Fact.out : Nat.Prime this._p.toNat).pos hdividendCanonical hselect
      have hquotientAbove := hquotientReady frontier hselect
      generalize hconsume : pairVecDivVHCConsumeEqualDegree this
          frontier.degree heap frontier.coefficient nodes #[] resetH quotient
          divisor = consumedExec at hrun
      cases consumedExec with
      | error fault =>
          simp only [Bind.bind, Except.bind] at hrun
          contradiction
      | ok consumed =>
          simp only [Bind.bind, Except.bind] at hrun
          have howned : PairVecDivVHCHeapChainsOwned heap nodes :=
            ⟨owners, hownership⟩
          have haway0 : PairVecDivVHCHeapChainsOwnedAway heap nodes
              (#[] : Array Nat).toList.toFinset := by
            simpa using howned.away_empty
          have hlin0 : PairVecDivVHCLinReady (#[] : Array Nat) nodes := by
            simp [PairVecDivVHCLinReady]
          have hconsumed :=
            pairVecDivVHCConsumeEqualDegree_preserves_away_linReady this
              frontier.degree heap frontier.coefficient nodes #[] resetH
              quotient divisor consumed haway0 hlin0 hconsume
          have hnodeInvariants :=
            pairVecDivVHCConsumeEqualDegree_preserves_node_invariants this
              this._p.toNat degreeLimit frontier.degree heap
              frontier.coefficient nodes #[]
              resetH quotient divisor consumed owners hcanonical hbelow
              hdenotes hfixed hready hownership hconsume
          generalize hemit : pairVecDivVHCEmit this frontier consumed quotient
              divisor hdivisor = emitExec at hrun
          cases emitExec with
          | error fault =>
              simp only [Bind.bind, Except.bind] at hrun
              contradiction
          | ok emitted =>
              simp only [Bind.bind, Except.bind] at hrun
              rcases emitted with ⟨quotient', activated, resetH'⟩
              have hemitted :=
                pairVecDivVHCEmit_preserves_away_linReady this frontier consumed
                  quotient divisor quotient' activated resetH' hconsumed.1
                  hconsumed.2 hnodeInvariants.2.2.2 hdivisor hemit
              have hemittedNodes :=
                pairVecDivVHCEmit_preserves_node_invariants this this._p.toNat
                  degreeLimit frontier consumed quotient divisor quotient'
                  activated resetH' hdivisorCanonical hdivisor
                  (Nat.le_of_lt hselectedBelow) hnodeInvariants.1
                  hnodeInvariants.2.1 hnodeInvariants.2.2.1
                  hnodeInvariants.2.2.2 hemit
              have hconsumedCoefficient :=
                pairVecDivVHCConsumeEqualDegree_coefficient_reduced this
                  frontier.degree heap frontier.coefficient nodes #[] resetH
                  quotient divisor consumed
                  (by
                    intro hp
                    have hzero : this._p.toNat = 0 :=
                      congrArg UInt64.toNat hp
                    exact (Fact.out : Nat.Prime this._p.toNat).ne_zero hzero)
                  hcfg hcanonical hfrontierCoefficient hconsume
              have hemittedCanonical := pairVecDivVHCEmit_preserves_canonical
                this frontier consumed quotient divisor quotient' activated
                resetH' hcfg hcanonical hconsumedCoefficient hdivisor
                hquotientAbove hemit
              generalize hreinsert : pairVecDivVHCReinsertLin activated.heap
                  activated.nodes consumed.lin = reinsertExec at hrun
              cases reinsertExec with
              | error fault =>
                  simp only [Bind.bind, Except.bind] at hrun
                  contradiction
              | ok reinserted =>
                  simp only [Bind.bind, Except.bind] at hrun
                  rw [← Except.ok.inj hrun]
                  have howned' :=
                    pairVecDivVHCReinsertLin_preserves_heapChainsOwned
                      activated.heap activated.nodes consumed.lin reinserted
                      hemitted.1 hemitted.2 hreinsert
                  have hinvariants' :=
                    pairVecDivVHCReinsertLin_preserves_node_invariants
                      degreeLimit resetH' quotient'.size activated.heap
                      activated.nodes consumed.lin quotient' divisor reinserted
                      hemittedNodes.1 hemittedNodes.2.1
                      hemittedNodes.2.2.1 hemittedNodes.2.2.2 hemitted.2
                      hreinsert
                  exact ⟨hemittedCanonical, howned', hremaining',
                    hinvariants'⟩

theorem pairVecDivVHCOuterIteration_preserves_quotientAbove_of_lt
    (this : DenseUPolyZp) (dividendIndex nextDegree : Nat)
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode)
    (quotient dividend divisor : SparsePolyZp) (resetH : Nat)
    (frontier : PairVecDivVHCFrontier)
    (result : PairVecDivVHCIterationResult) (hdivisor : 0 < divisor.size)
    (habove : PairVecDivVHCQuotientAbove frontier.degree divisor[0].1.deg
      quotient)
    (hnext : nextDegree < frontier.degree)
    (hselect : pairVecDivVHCSelectFrontier dividendIndex dividend heap nodes =
      .ok frontier)
    (hrun : pairVecDivVHCOuterIteration this dividendIndex heap nodes quotient
      dividend divisor resetH = .ok result) :
    PairVecDivVHCQuotientAbove nextDegree divisor[0].1.deg result.quotient := by
  unfold pairVecDivVHCOuterIteration at hrun
  simp only [hdivisor, ↓reduceDIte, hselect, Bind.bind, Except.bind] at hrun
  generalize hconsume : pairVecDivVHCConsumeEqualDegree this frontier.degree
      heap frontier.coefficient nodes #[] resetH quotient divisor =
        consumedExec at hrun
  cases consumedExec with
  | error fault => contradiction
  | ok consumed =>
      simp only [Bind.bind, Except.bind] at hrun
      generalize hemit : pairVecDivVHCEmit this frontier consumed quotient
          divisor hdivisor = emitExec at hrun
      cases emitExec with
      | error fault => contradiction
      | ok emitted =>
          simp only [Bind.bind, Except.bind] at hrun
          rcases emitted with ⟨quotient', activated, resetH'⟩
          generalize hreinsert : pairVecDivVHCReinsertLin activated.heap
              activated.nodes consumed.lin = reinsertExec at hrun
          cases reinsertExec with
          | error fault => contradiction
          | ok reinserted =>
              simp only [Bind.bind, Except.bind] at hrun
              rw [← Except.ok.inj hrun]
              exact pairVecDivVHCEmit_preserves_quotientAbove_of_lt this
                nextDegree frontier consumed quotient divisor quotient'
                activated resetH' hdivisor habove hnext hemit

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

theorem pairVecDivVHCOuterLoop_step_of_success (this : DenseUPolyZp)
    (degreeLimit dividendIndex : Nat) (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode) (quotient dividend divisor : SparsePolyZp)
    (resetH : Nat) (frontier : PairVecDivVHCFrontier)
    (output : SparsePolyZp)
    (hnotDone : ¬ (dividend.size ≤ dividendIndex ∧ heap.size = 0))
    (hselect : pairVecDivVHCSelectFrontier dividendIndex dividend heap nodes =
      .ok frontier)
    (hrun : pairVecDivVHCOuterLoop this degreeLimit dividendIndex heap nodes
      quotient dividend divisor resetH = .ok output) :
    frontier.degree < degreeLimit ∧
      ∃ next : PairVecDivVHCIterationResult,
        pairVecDivVHCOuterIteration this dividendIndex heap nodes quotient
            dividend divisor resetH = .ok next ∧
          pairVecDivVHCOuterLoop this frontier.degree next.dividendIndex
            next.heap next.nodes next.quotient dividend divisor next.resetH =
              .ok output := by
  rw [pairVecDivVHCOuterLoop] at hrun
  simp only [hnotDone, ↓reduceDIte, hselect] at hrun
  by_cases hdecrease : frontier.degree < degreeLimit
  · simp only [hdecrease, ↓reduceDIte] at hrun
    cases hiteration : pairVecDivVHCOuterIteration this dividendIndex heap
        nodes quotient dividend divisor resetH with
    | error fault => simp [hiteration] at hrun
    | ok next =>
        have hrecursive : pairVecDivVHCOuterLoop this frontier.degree
            next.dividendIndex next.heap next.nodes next.quotient dividend
              divisor next.resetH = .ok output := by
          simpa only [hiteration] using hrun
        exact ⟨hdecrease, next, rfl, hrecursive⟩
  · simp [hdecrease] at hrun

theorem pairVecDivVHCOuterLoop_selected_degree_lt_of_success
    (this : DenseUPolyZp) (degreeLimit dividendIndex : Nat)
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode)
    (quotient dividend divisor output : SparsePolyZp) (resetH : Nat)
    (frontier : PairVecDivVHCFrontier)
    (hselect : pairVecDivVHCSelectFrontier dividendIndex dividend heap nodes =
      .ok frontier)
    (hrun : pairVecDivVHCOuterLoop this degreeLimit dividendIndex heap nodes
      quotient dividend divisor resetH = .ok output) :
    frontier.degree < degreeLimit := by
  by_cases hdone : dividend.size ≤ dividendIndex ∧ heap.size = 0
  · rcases hdone with ⟨hdividend, hheap⟩
    unfold pairVecDivVHCSelectFrontier at hselect
    have hdividendDone : ¬ dividendIndex < dividend.size := by omega
    simp [hdividendDone, hheap] at hselect
  · exact (pairVecDivVHCOuterLoop_step_of_success this degreeLimit
      dividendIndex heap nodes quotient dividend divisor resetH frontier output
      hdone hselect hrun).1

theorem pairVecDivVHCOuterLoop_preserves_canonical_of_success
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (globalLimit degreeLimit dividendIndex : Nat)
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode)
    (quotient dividend divisor output : SparsePolyZp) (resetH : Nat)
    (owners : Nat → Finset Nat) (hdivisor : 0 < divisor.size)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hcanonical : SparsePolyZp.Canonical this._p.toNat quotient)
    (hdividendCanonical : SparsePolyZp.Canonical this._p.toNat dividend)
    (hdivisorCanonical : SparsePolyZp.Canonical this._p.toNat divisor)
    (hquotientReady : ∀ frontier : PairVecDivVHCFrontier,
      pairVecDivVHCSelectFrontier dividendIndex dividend heap nodes =
        .ok frontier →
      PairVecDivVHCQuotientAbove frontier.degree divisor[0].1.deg quotient)
    (hremaining : PairVecDivVHCRemainingDividendBelow globalLimit
      dividendIndex dividend)
    (hbelow : PairVecDivVHCAllActiveNodesBelow globalLimit nodes)
    (hdenotes : ∀ (i : Nat) (node : PairVecDivVHCNode),
      nodes[i]? = some node → node.mono ≠ none →
        PairVecDivVHCNodeDenotes quotient divisor node)
    (hfixed : PairVecDivVHCNodeDivisorIndicesFixed nodes)
    (hready : PairVecDivVHCResetReady resetH quotient.size nodes)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hcovered : PairVecDivVHCStateCovered heap nodes #[] resetH)
    (hrun : pairVecDivVHCOuterLoop this degreeLimit dividendIndex heap nodes
      quotient dividend divisor resetH = .ok output) :
    SparsePolyZp.Canonical this._p.toNat output := by
  induction degreeLimit using Nat.strong_induction_on generalizing dividendIndex
      heap nodes quotient resetH owners with
  | h degreeLimit ih =>
      by_cases hdone : dividend.size ≤ dividendIndex ∧ heap.size = 0
      · rw [pairVecDivVHCOuterLoop] at hrun
        rw [dif_pos hdone] at hrun
        rw [← Except.ok.inj hrun]
        exact hcanonical
      · cases hselect : pairVecDivVHCSelectFrontier dividendIndex dividend
            heap nodes with
        | error fault =>
            rw [pairVecDivVHCOuterLoop] at hrun
            rw [dif_neg hdone, hselect] at hrun
            contradiction
        | ok frontier =>
            rcases pairVecDivVHCOuterLoop_step_of_success this degreeLimit
                dividendIndex heap nodes quotient dividend divisor resetH
                frontier output hdone hselect hrun with
              ⟨hdecrease, next, hiteration, hrecursive⟩
            have hiterationInv :=
              pairVecDivVHCOuterIteration_preserves_heapChainsOwned this
                globalLimit dividendIndex heap nodes quotient dividend divisor
                resetH next owners hdivisor hcfg hcanonical
                hdividendCanonical hdivisorCanonical hquotientReady hremaining
                hbelow hdenotes hfixed hready hownership hiteration
            have hnextCovered :=
              pairVecDivVHCOuterIteration_preserves_stateCovered this
                this._p.toNat globalLimit dividendIndex heap nodes quotient
                dividend divisor resetH frontier next hdivisor hselect
                hcanonical hbelow hdenotes hfixed hready hcovered hiteration
            rcases hiterationInv with
              ⟨hnextCanonical, ⟨nextOwners, hnextOwnership⟩,
                hnextRemaining, hnextBelow, hnextDenotes, hnextFixed,
                hnextReady⟩
            have hnextQuotientReady : ∀ nextFrontier :
                PairVecDivVHCFrontier,
                pairVecDivVHCSelectFrontier next.dividendIndex dividend
                    next.heap next.nodes = .ok nextFrontier →
                  PairVecDivVHCQuotientAbove nextFrontier.degree
                    divisor[0].1.deg next.quotient := by
              intro nextFrontier hnextSelect
              have hnextDecrease :=
                pairVecDivVHCOuterLoop_selected_degree_lt_of_success this
                  frontier.degree next.dividendIndex next.heap next.nodes
                  next.quotient dividend divisor output next.resetH nextFrontier
                  hnextSelect hrecursive
              exact
                pairVecDivVHCOuterIteration_preserves_quotientAbove_of_lt this
                  dividendIndex nextFrontier.degree heap nodes quotient dividend
                  divisor resetH frontier next hdivisor
                  (hquotientReady frontier hselect) hnextDecrease hselect
                  hiteration
            exact ih frontier.degree hdecrease next.dividendIndex next.heap
              next.nodes next.quotient next.resetH nextOwners hnextCanonical
              hnextQuotientReady hnextRemaining hnextBelow hnextDenotes
              hnextFixed hnextReady hnextOwnership hnextCovered hrecursive

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

theorem listSum_coeff_zero_of_degree_absent (p degree : Nat)
    (terms : List (UMonomial × Zp))
    (habsent : ∀ term ∈ terms, term.1.deg ≠ degree) :
    (listSum p terms).coeff degree = 0 := by
  induction terms with
  | nil => simp [listSum]
  | cons term rest ih =>
      rw [listSum_cons, Polynomial.coeff_add, Polynomial.coeff_monomial,
        if_neg (habsent term List.mem_cons_self),
        ih (by
          intro item hitem
          exact habsent item (List.mem_cons_of_mem term hitem))]
      simp

/-- The selector's machine-word coefficient is the actual dividend
polynomial coefficient at the selected degree.  In the heap-source case the
proof excludes that degree from both sides of the concrete dividend pointer:
the processed prefix is above the strict loop bound, while canonical ordering
and the selector comparison put the unprocessed suffix below the heap root. -/
theorem pairVecDivVHCSelectFrontier_coefficient_toPoly
    (p degreeLimit dividendIndex : Nat) (dividend : SparsePolyZp)
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode)
    (frontier : PairVecDivVHCFrontier)
    (hcanonical : SparsePolyZp.Canonical p dividend)
    (hconsumed : PairVecDivVHCConsumedDividendAbove degreeLimit dividendIndex
      dividend)
    (hdecrease : frontier.degree < degreeLimit)
    (hrun : pairVecDivVHCSelectFrontier dividendIndex dividend heap nodes =
      .ok frontier) :
    (frontier.coefficient.toNat : ZMod p) =
      (SparsePolyZp.toPoly p dividend).coeff frontier.degree := by
  have hsource := pairVecDivVHCSelectFrontier_has_source dividendIndex
    dividend heap nodes frontier hrun
  cases hsource with
  | dividend hindex =>
      have htermMem : dividend[dividendIndex] ∈ dividend.toList :=
        Array.getElem_mem_toList hindex
      have hcoefficient := listSum_coeff_of_mem_chain p dividend.toList
        dividend[dividendIndex] hcanonical.2.1 htermMem
      unfold SparsePolyZp.toPoly
      rw [hcoefficient]
      rfl
  | heap hheap rootMono hmono hdominates =>
      change rootMono.deg < degreeLimit at hdecrease
      change ((0 : UInt64).toNat : ZMod p) =
        (SparsePolyZp.toPoly p dividend).coeff rootMono.deg
      rw [show (0 : UInt64).toNat = 0 by rfl, Nat.cast_zero]
      symm
      unfold SparsePolyZp.toPoly
      apply listSum_coeff_zero_of_degree_absent
      intro term hterm
      rcases List.getElem_of_mem hterm with ⟨i, hiList, hgetList⟩
      have hi : i < dividend.size := by simpa using hiList
      have htermEq : term = dividend[i] := by
        rw [← Array.getElem_toList hi]
        exact hgetList.symm
      clear hgetList hterm
      subst term
      have hget : dividend[i]? = some dividend[i] :=
        Array.getElem?_eq_getElem hi
      by_cases hprefix : i < dividendIndex
      · have habove := hconsumed i dividend[i] hprefix hget
        omega
      · have hsuffix : dividendIndex ≤ i := by omega
        have hindex : dividendIndex < dividend.size :=
          Nat.lt_of_le_of_lt hsuffix hi
        have hheadBelow := hdominates hindex
        by_cases heq : i = dividendIndex
        · subst i
          omega
        · have htail := canonical_degree_lt_of_index_lt p dividend
            hcanonical dividendIndex i hindex hi (by omega)
          omega

theorem pairVecDivVHCOuterIteration_residual_coefficient
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (degreeLimit dividendIndex : Nat) (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode)
    (quotient dividend divisor : SparsePolyZp) (resetH : Nat)
    (frontier : PairVecDivVHCFrontier)
    (result : PairVecDivVHCIterationResult) (owners : Nat → Finset Nat)
    (hdivisor : 0 < divisor.size)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hhomogeneous : PairVecDivVHCHeapChainsHomogeneous heap owners nodes)
    (hdenotes : ∀ (i : Nat) (node : PairVecDivVHCNode),
      nodes[i]? = some node → node.mono ≠ none →
        PairVecDivVHCNodeDenotes quotient divisor node)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hquotientCanonical : SparsePolyZp.Canonical this._p.toNat quotient)
    (hdividendCanonical : SparsePolyZp.Canonical this._p.toNat dividend)
    (hconsumed : PairVecDivVHCConsumedDividendAbove degreeLimit dividendIndex
      dividend)
    (hdecrease : frontier.degree < degreeLimit)
    (hselect : pairVecDivVHCSelectFrontier dividendIndex dividend heap nodes =
      .ok frontier)
    (hrun : pairVecDivVHCOuterIteration this dividendIndex heap nodes quotient
      dividend divisor resetH = .ok result) :
    ∃ consumed products,
      pairVecDivVHCConsumeEqualDegree this frontier.degree heap
          frontier.coefficient nodes #[] resetH quotient divisor =
        .ok consumed ∧
      (consumed.coefficient.toNat : ZMod this._p.toNat) =
        (SparsePolyZp.toPoly this._p.toNat dividend).coeff frontier.degree -
          pairVecDivVHCProductsValue this._p.toNat products ∧
      ∀ product ∈ products,
        PairVecDivVHCStoredProductAtDegree frontier.degree quotient divisor
          product := by
  rcases pairVecDivVHCOuterIteration_components this dividendIndex heap nodes
      quotient dividend divisor resetH frontier result hdivisor hselect hrun with
    ⟨consumed, quotient', activated, resetH', reinserted, hconsume, hemit,
      hreinsert, hresult⟩
  have hk := pairVecDivVHCSelectFrontier_coefficient_reduced
    this._p.toNat dividendIndex dividend heap nodes frontier
      (Fact.out : Nat.Prime this._p.toNat).pos hdividendCanonical hselect
  rcases pairVecDivVHCConsumeEqualDegree_coefficient_semantics this
      frontier.degree heap frontier.coefficient nodes #[] resetH quotient
      divisor consumed owners hownership hhomogeneous hdenotes hcfg
      hquotientCanonical hk hconsume with
    ⟨products, hcoefficient, hproducts⟩
  have hfrontier := pairVecDivVHCSelectFrontier_coefficient_toPoly
    this._p.toNat degreeLimit dividendIndex dividend heap nodes frontier
      hdividendCanonical hconsumed hdecrease hselect
  rw [hfrontier] at hcoefficient
  exact ⟨consumed, products, hconsume, hcoefficient, hproducts⟩

/-- Polynomial form of the real outer-iteration subtraction.  Unlike the
legacy product-list statement, this theorem uses exact owner multiplicities
and the cursor/index bijection to identify the generated heap consumption
with the L2 coefficient of `quotient * divisor.tail`. -/
theorem pairVecDivVHCOuterIteration_residual_coefficient_toPoly
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (degreeLimit dividendIndex : Nat) (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode)
    (quotient dividend divisor : SparsePolyZp) (resetH : Nat)
    (frontier : PairVecDivVHCFrontier)
    (result : PairVecDivVHCIterationResult) (owners : Nat → Finset Nat)
    (hdivisor : 0 < divisor.size)
    (hsize : nodes.size = divisor.size - 1)
    (hstate : PairVecDivVHCStateCovered heap nodes #[] resetH)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hhomogeneous : PairVecDivVHCHeapChainsHomogeneous heap owners nodes)
    (hordered : PairVecDivVHCHeapOrdered heap nodes)
    (hdenotes : ∀ (i : Nat) (node : PairVecDivVHCNode),
      nodes[i]? = some node → node.mono ≠ none →
        PairVecDivVHCNodeDenotes quotient divisor node)
    (hfixed : PairVecDivVHCNodeDivisorIndicesFixed nodes)
    (hresetReady : PairVecDivVHCResetReady resetH quotient.size nodes)
    (hprefix : PairVecDivVHCCursorPrefixAbove degreeLimit nodes quotient divisor)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hquotientCanonical : SparsePolyZp.Canonical this._p.toNat quotient)
    (hdividendCanonical : SparsePolyZp.Canonical this._p.toNat dividend)
    (hconsumed : PairVecDivVHCConsumedDividendAbove degreeLimit dividendIndex
      dividend)
    (hdecrease : frontier.degree < degreeLimit)
    (hselect : pairVecDivVHCSelectFrontier dividendIndex dividend heap nodes =
      .ok frontier)
    (hrun : pairVecDivVHCOuterIteration this dividendIndex heap nodes quotient
      dividend divisor resetH = .ok result) :
    ∃ consumed,
      pairVecDivVHCConsumeEqualDegree this frontier.degree heap
          frontier.coefficient nodes #[] resetH quotient divisor =
        .ok consumed ∧
      (consumed.coefficient.toNat : ZMod this._p.toNat) =
        (SparsePolyZp.toPoly this._p.toNat dividend).coeff frontier.degree -
          (SparsePolyZp.toPoly this._p.toNat quotient *
            listSum this._p.toNat divisor.toList.tail).coeff
              frontier.degree := by
  rcases pairVecDivVHCOuterIteration_components this dividendIndex heap nodes
      quotient dividend divisor resetH frontier result hdivisor hselect hrun with
    ⟨consumed, quotient', activated, resetH', reinserted, hconsume, hemit,
      hreinsert, hresult⟩
  have hk := pairVecDivVHCSelectFrontier_coefficient_reduced
    this._p.toNat dividendIndex dividend heap nodes frontier
      (Fact.out : Nat.Prime this._p.toNat).pos hdividendCanonical hselect
  have hmax : ∀ (slot head : Nat) (mono : UMonomial),
      heap[slot]? = some head → pairVecDivVHCMono head nodes = .ok mono →
        mono.deg ≤ frontier.degree := by
    intro slot head mono hheap hmono
    have hslot : slot < heap.size := by
      by_contra hnot
      rw [Array.getElem?_eq_none (by omega)] at hheap
      contradiction
    have hheadEq : heap[slot] = head := by
      rw [Array.getElem?_eq_getElem hslot] at hheap
      exact Option.some.inj hheap
    have hslotMono : pairVecDivVHCMono heap[slot] nodes = .ok mono := by
      simpa [hheadEq] using hmono
    exact pairVecDivVHCSelectFrontier_heap_slot_degree_le dividendIndex
      dividend heap nodes frontier (hownership.heapPointersValid heap owners
        nodes) hordered slot hslot mono hslotMono hselect
  rcases pairVecDivVHCConsumeEqualDegree_products_complete this
      frontier.degree heap frontier.coefficient nodes #[] resetH quotient
      divisor consumed owners hownership hhomogeneous hordered hmax hdenotes
      hcfg hquotientCanonical hk hconsume with
    ⟨products, hcoefficient, hsound, hcover, hproductsValue⟩
  have hfrontierCoefficient := pairVecDivVHCSelectFrontier_coefficient_toPoly
    this._p.toNat degreeLimit dividendIndex dividend heap nodes frontier
      hdividendCanonical hconsumed hdecrease hselect
  have hownerCoefficient := pairVecDivVHCHeapOwnerSum_eq_productCoeffTail
    this._p.toNat degreeLimit dividendIndex resetH dividend quotient divisor
    heap nodes owners frontier hsize hfixed hstate hownership hhomogeneous
    hresetReady hordered hdenotes hquotientCanonical hprefix hdecrease hselect
  refine ⟨consumed, hconsume, ?_⟩
  rw [hcoefficient, hproductsValue, hfrontierCoefficient, hownerCoefficient]

/-- One successful generated outer iteration makes the complete
`quotient * divisor` coefficient agree with the dividend at the selected
frontier.  This combines the old quotient/tail heap consumption with the
newly emitted quotient/lead contribution; both cross terms are proved zero
from the quotient-above and canonical divisor-order invariants. -/
theorem pairVecDivVHCOuterIteration_product_coefficient
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (degreeLimit dividendIndex : Nat) (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode)
    (quotient dividend divisor : SparsePolyZp) (resetH : Nat)
    (frontier : PairVecDivVHCFrontier)
    (result : PairVecDivVHCIterationResult) (owners : Nat → Finset Nat)
    (hdivisor : 0 < divisor.size)
    (hsize : nodes.size = divisor.size - 1)
    (hstate : PairVecDivVHCStateCovered heap nodes #[] resetH)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hhomogeneous : PairVecDivVHCHeapChainsHomogeneous heap owners nodes)
    (hordered : PairVecDivVHCHeapOrdered heap nodes)
    (hdenotes : ∀ (i : Nat) (node : PairVecDivVHCNode),
      nodes[i]? = some node → node.mono ≠ none →
        PairVecDivVHCNodeDenotes quotient divisor node)
    (hfixed : PairVecDivVHCNodeDivisorIndicesFixed nodes)
    (hresetReady : PairVecDivVHCResetReady resetH quotient.size nodes)
    (hprefix : PairVecDivVHCCursorPrefixAbove degreeLimit nodes quotient divisor)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hquotientCanonical : SparsePolyZp.Canonical this._p.toNat quotient)
    (hdividendCanonical : SparsePolyZp.Canonical this._p.toNat dividend)
    (hdivisorCanonical : SparsePolyZp.Canonical this._p.toNat divisor)
    (habove : PairVecDivVHCQuotientAbove frontier.degree divisor[0].1.deg
      quotient)
    (hconsumed : PairVecDivVHCConsumedDividendAbove degreeLimit dividendIndex
      dividend)
    (hdecrease : frontier.degree < degreeLimit)
    (hlead : divisor[0].1.deg ≤ frontier.degree)
    (hselect : pairVecDivVHCSelectFrontier dividendIndex dividend heap nodes =
      .ok frontier)
    (hrun : pairVecDivVHCOuterIteration this dividendIndex heap nodes quotient
      dividend divisor resetH = .ok result) :
    (SparsePolyZp.toPoly this._p.toNat result.quotient *
      SparsePolyZp.toPoly this._p.toNat divisor).coeff frontier.degree =
        (SparsePolyZp.toPoly this._p.toNat dividend).coeff frontier.degree := by
  rcases pairVecDivVHCOuterIteration_components this dividendIndex heap nodes
      quotient dividend divisor resetH frontier result hdivisor hselect hrun with
    ⟨consumed, quotient', activated, resetH', reinserted, hconsume, hemit,
      hreinsert, hresult⟩
  rcases pairVecDivVHCOuterIteration_residual_coefficient_toPoly this
      degreeLimit dividendIndex heap nodes quotient dividend divisor resetH
      frontier result owners hdivisor hsize hstate hownership hhomogeneous
      hordered hdenotes hfixed hresetReady hprefix hcfg hquotientCanonical
      hdividendCanonical hconsumed hdecrease hselect hrun with
    ⟨residualConsumed, hresidualConsume, hresidual⟩
  have hconsumedEq : residualConsumed = consumed := by
    rw [hconsume] at hresidualConsume
    exact (Except.ok.inj hresidualConsume).symm
  subst residualConsumed
  have hk := pairVecDivVHCSelectFrontier_coefficient_reduced
    this._p.toNat dividendIndex dividend heap nodes frontier
      (Fact.out : Nat.Prime this._p.toNat).pos hdividendCanonical hselect
  have hp : this._p ≠ 0 := by
    intro hp
    have hz : this._p.toNat = 0 := congrArg UInt64.toNat hp
    exact (Fact.out : Nat.Prime this._p.toNat).ne_zero hz
  have hcoefficientReduced :=
    pairVecDivVHCConsumeEqualDegree_coefficient_reduced this frontier.degree
      heap frontier.coefficient nodes #[] resetH quotient divisor consumed hp
      hcfg hquotientCanonical hk hconsume
  have hemitPoly := pairVecDivVHCEmit_toPoly_of_lead_le this frontier consumed
    quotient divisor quotient' activated resetH' hcfg hdivisorCanonical
    hcoefficientReduced hdivisor hlead hemit
  have hdivisorPoly := sparsePolyZp_toPoly_eq_head_add_tail this._p.toNat
    divisor hdivisor
  have holdLead := pairVecDivVHCQuotient_mul_lead_coeff_eq_zero this._p.toNat
    frontier.degree divisor[0].1.deg quotient divisor[0].2 habove hlead
  let emittedCoefficient : Zp := ⟨consumed.coefficient, this._p⟩
  have hnewTail := pairVecDivVHCEmitted_monomial_mul_tail_coeff_eq_zero
    this._p.toNat frontier.degree divisor
    (Zp.toZMod this._p.toNat emittedCoefficient /
      Zp.toZMod this._p.toNat divisor[0].2) hdivisor hdivisorCanonical hlead
  have hleadFieldNonzero : Zp.toZMod this._p.toNat divisor[0].2 ≠ 0 :=
    Zp.toZMod_ne_zero_of_val_ne_zero this._p.toNat divisor[0].2
      (hdivisorCanonical.1 divisor[0]
        (Array.getElem_mem_toList hdivisor))
      (hdivisorCanonical.2.2 divisor[0]
        (Array.getElem_mem_toList hdivisor))
  have hnewLead :
      (Polynomial.monomial (frontier.degree - divisor[0].1.deg)
          (Zp.toZMod this._p.toNat emittedCoefficient /
            Zp.toZMod this._p.toNat divisor[0].2) *
        Polynomial.monomial divisor[0].1.deg
          (Zp.toZMod this._p.toNat divisor[0].2)).coeff frontier.degree =
        Zp.toZMod this._p.toNat emittedCoefficient := by
    rw [Polynomial.monomial_mul_monomial, Nat.sub_add_cancel hlead,
      div_mul_cancel₀ _ hleadFieldNonzero, Polynomial.coeff_monomial]
    simp
  rw [hresult]
  simp only [PairVecDivVHCIterationResult.quotient]
  rw [hemitPoly, hdivisorPoly, add_mul, mul_add, mul_add,
    Polynomial.coeff_add, Polynomial.coeff_add, Polynomial.coeff_add,
    holdLead, hnewTail, hnewLead]
  change 0 + (SparsePolyZp.toPoly this._p.toNat quotient *
      listSum this._p.toNat divisor.toList.tail).coeff frontier.degree +
      (Zp.toZMod this._p.toNat emittedCoefficient + 0) = _
  simp only [zero_add, add_zero]
  change _ + (consumed.coefficient.toNat : ZMod this._p.toNat) = _
  rw [hresidual]
  ring

theorem pairVecDivVHCOuterIteration_product_coeff_above
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (dividendIndex targetDegree : Nat) (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode)
    (quotient dividend divisor : SparsePolyZp) (resetH : Nat)
    (frontier : PairVecDivVHCFrontier)
    (result : PairVecDivVHCIterationResult)
    (hdivisor : 0 < divisor.size)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hquotientCanonical : SparsePolyZp.Canonical this._p.toNat quotient)
    (hdividendCanonical : SparsePolyZp.Canonical this._p.toNat dividend)
    (hdivisorCanonical : SparsePolyZp.Canonical this._p.toNat divisor)
    (hlead : divisor[0].1.deg ≤ frontier.degree)
    (habove : frontier.degree < targetDegree)
    (hselect : pairVecDivVHCSelectFrontier dividendIndex dividend heap nodes =
      .ok frontier)
    (hrun : pairVecDivVHCOuterIteration this dividendIndex heap nodes quotient
      dividend divisor resetH = .ok result) :
    (SparsePolyZp.toPoly this._p.toNat result.quotient *
        SparsePolyZp.toPoly this._p.toNat divisor).coeff targetDegree =
      (SparsePolyZp.toPoly this._p.toNat quotient *
        SparsePolyZp.toPoly this._p.toNat divisor).coeff targetDegree := by
  rcases pairVecDivVHCOuterIteration_components this dividendIndex heap nodes
      quotient dividend divisor resetH frontier result hdivisor hselect hrun with
    ⟨consumed, quotient', activated, resetH', reinserted, hconsume, hemit,
      hreinsert, hresult⟩
  have hk := pairVecDivVHCSelectFrontier_coefficient_reduced
    this._p.toNat dividendIndex dividend heap nodes frontier
      (Fact.out : Nat.Prime this._p.toNat).pos hdividendCanonical hselect
  have hp : this._p ≠ 0 := by
    intro hp
    have hz : this._p.toNat = 0 := congrArg UInt64.toNat hp
    exact (Fact.out : Nat.Prime this._p.toNat).ne_zero hz
  have hcoefficientReduced :=
    pairVecDivVHCConsumeEqualDegree_coefficient_reduced this frontier.degree
      heap frontier.coefficient nodes #[] resetH quotient divisor consumed hp
      hcfg hquotientCanonical hk hconsume
  have hemitAbove := pairVecDivVHCEmit_product_coeff_above this frontier
    consumed quotient divisor quotient' activated resetH' targetDegree hcfg
    hdivisorCanonical hcoefficientReduced hdivisor hlead habove hemit
  rw [hresult]
  exact hemitAbove

def PairVecDivVHCProductAgreesAbove (p degreeLimit : Nat)
    (quotient dividend divisor : SparsePolyZp) : Prop :=
  ∀ degree, degreeLimit ≤ degree →
    (SparsePolyZp.toPoly p quotient * SparsePolyZp.toPoly p divisor).coeff
        degree =
      (SparsePolyZp.toPoly p dividend).coeff degree

/-- One real outer iteration extends coefficient agreement from the previous
strict bound down through the selected frontier.  Degrees above the old bound
use the induction hypothesis, the open sparse gap uses the concrete
dividend/heap exclusion lemmas, and the selected degree uses the full
generated consume/emit coefficient theorem. -/
theorem pairVecDivVHCOuterIteration_extends_productAgreement
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (degreeLimit dividendIndex : Nat) (heap : Array Nat)
    (nodes : Array PairVecDivVHCNode)
    (quotient dividend divisor : SparsePolyZp) (resetH : Nat)
    (frontier : PairVecDivVHCFrontier)
    (result : PairVecDivVHCIterationResult) (owners : Nat → Finset Nat)
    (hdivisor : 0 < divisor.size)
    (hsize : nodes.size = divisor.size - 1)
    (hstate : PairVecDivVHCStateCovered heap nodes #[] resetH)
    (hownership : PairVecDivVHCHeapChainOwnership heap owners nodes)
    (hhomogeneous : PairVecDivVHCHeapChainsHomogeneous heap owners nodes)
    (hordered : PairVecDivVHCHeapOrdered heap nodes)
    (hdenotes : ∀ (i : Nat) (node : PairVecDivVHCNode),
      nodes[i]? = some node → node.mono ≠ none →
        PairVecDivVHCNodeDenotes quotient divisor node)
    (hfixed : PairVecDivVHCNodeDivisorIndicesFixed nodes)
    (hresetReady : PairVecDivVHCResetReady resetH quotient.size nodes)
    (hprefix : PairVecDivVHCCursorPrefixAbove degreeLimit nodes quotient divisor)
    (hprocessed : PairVecDivVHCQuotientLeadAbove degreeLimit
      divisor[0].1.deg quotient)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hquotientCanonical : SparsePolyZp.Canonical this._p.toNat quotient)
    (hdividendCanonical : SparsePolyZp.Canonical this._p.toNat dividend)
    (hdivisorCanonical : SparsePolyZp.Canonical this._p.toNat divisor)
    (habove : PairVecDivVHCQuotientAbove frontier.degree divisor[0].1.deg
      quotient)
    (hconsumed : PairVecDivVHCConsumedDividendAbove degreeLimit dividendIndex
      dividend)
    (hdecrease : frontier.degree < degreeLimit)
    (hlead : divisor[0].1.deg ≤ frontier.degree)
    (hagrees : PairVecDivVHCProductAgreesAbove this._p.toNat degreeLimit
      quotient dividend divisor)
    (hselect : pairVecDivVHCSelectFrontier dividendIndex dividend heap nodes =
      .ok frontier)
    (hrun : pairVecDivVHCOuterIteration this dividendIndex heap nodes quotient
      dividend divisor resetH = .ok result) :
    PairVecDivVHCProductAgreesAbove this._p.toNat frontier.degree
      result.quotient dividend divisor := by
  intro degree hdegree
  by_cases heq : degree = frontier.degree
  · subst degree
    exact pairVecDivVHCOuterIteration_product_coefficient this degreeLimit
      dividendIndex heap nodes quotient dividend divisor resetH frontier result
      owners hdivisor hsize hstate hownership hhomogeneous hordered hdenotes
      hfixed hresetReady hprefix hcfg hquotientCanonical hdividendCanonical
      hdivisorCanonical habove hconsumed hdecrease hlead hselect hrun
  · have hstrict : frontier.degree < degree := by omega
    rw [pairVecDivVHCOuterIteration_product_coeff_above this dividendIndex degree
      heap nodes quotient dividend divisor resetH frontier result hdivisor hcfg
      hquotientCanonical hdividendCanonical hdivisorCanonical hlead hstrict
      hselect hrun]
    by_cases hold : degreeLimit ≤ degree
    · exact hagrees degree hold
    · have hgapProduct := pairVecDivVHCProduct_coeff_eq_zero_of_gap
        this._p.toNat degreeLimit degree dividendIndex resetH dividend quotient
        divisor heap nodes owners frontier hdivisor hsize hfixed hstate
        hownership hhomogeneous hresetReady hordered hdenotes
        hquotientCanonical hprefix hprocessed hstrict (by omega) hselect
      have hgapDividend := pairVecDivVHCDividend_coeff_eq_zero_of_gap
        this._p.toNat degreeLimit degree dividendIndex dividend heap nodes
        frontier hdividendCanonical hconsumed hstrict (by omega) hselect
      rw [hgapProduct, hgapDividend]

theorem pairVecDivVHCOuterIteration_preserves_quotientLeadAbove
    (this : DenseUPolyZp) (degreeLimit dividendIndex : Nat)
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode)
    (quotient dividend divisor : SparsePolyZp) (resetH : Nat)
    (frontier : PairVecDivVHCFrontier)
    (result : PairVecDivVHCIterationResult)
    (hdivisor : 0 < divisor.size)
    (hprocessed : PairVecDivVHCQuotientLeadAbove degreeLimit
      divisor[0].1.deg quotient)
    (hdecrease : frontier.degree < degreeLimit)
    (hselect : pairVecDivVHCSelectFrontier dividendIndex dividend heap nodes =
      .ok frontier)
    (hrun : pairVecDivVHCOuterIteration this dividendIndex heap nodes quotient
      dividend divisor resetH = .ok result) :
    PairVecDivVHCQuotientLeadAbove frontier.degree divisor[0].1.deg
      result.quotient := by
  rcases pairVecDivVHCOuterIteration_components this dividendIndex heap nodes
      quotient dividend divisor resetH frontier result hdivisor hselect hrun with
    ⟨consumed, quotient', activated, resetH', reinserted, hconsume, hemit,
      hreinsert, hresult⟩
  have hemitted := pairVecDivVHCEmit_preserves_quotientLeadAbove this
    degreeLimit frontier consumed quotient divisor quotient' activated resetH'
    hdivisor hprocessed hdecrease hemit
  rw [hresult]
  exact hemitted

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
