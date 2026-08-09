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

theorem pairVecDivVHCInsert_get_ne
    (newNode : Nat) (heap heap' : Array Nat)
    (nodes nodes' : Array PairVecDivVHCNode)
    (hrun : pairVecDivVHCInsert newNode heap nodes = .ok (heap', nodes'))
    (i : Nat) (hne : newNode ≠ i) :
    nodes'[i]? = nodes[i]? := by
  rcases pairVecDivVHCInsert_nodes_result newNode heap heap' nodes nodes' hrun
    with ⟨next, hset⟩
  exact pairVecDivVHCSetNext_get_ne newNode next nodes nodes' hset i hne

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

/-- Every deferred source node occurs exactly once in `lin` and is already
active.  This is the proof state required by the literal reverse reinsertion
loop; it contains no execution budget or specification-level result. -/
def PairVecDivVHCLinReady (lin : Array Nat)
    (nodes : Array PairVecDivVHCNode) : Prop :=
  lin.toList.Nodup ∧
    ∀ nodeIndex ∈ lin.toList,
      ∃ node mono, nodes[nodeIndex]? = some node ∧ node.mono = some mono

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

theorem pairVecDivVHCOuterIteration_preserves_heapChainsOwned
    (this : DenseUPolyZp) (p degreeLimit dividendIndex : Nat)
    (heap : Array Nat) (nodes : Array PairVecDivVHCNode)
    (quotient dividend divisor : SparsePolyZp) (resetH : Nat)
    (result : PairVecDivVHCIterationResult) (owners : Nat → Finset Nat)
    (hdivisor : 0 < divisor.size)
    (hcanonical : SparsePolyZp.Canonical p quotient)
    (hdivisorCanonical : SparsePolyZp.Canonical p divisor)
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
            pairVecDivVHCConsumeEqualDegree_preserves_node_invariants this p
              degreeLimit frontier.degree heap frontier.coefficient nodes #[]
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
                pairVecDivVHCEmit_preserves_node_invariants this p degreeLimit
                  frontier consumed quotient divisor quotient' activated
                  resetH' hdivisorCanonical hdivisor
                  (Nat.le_of_lt hselectedBelow) hnodeInvariants.1
                  hnodeInvariants.2.1 hnodeInvariants.2.2.1
                  hnodeInvariants.2.2.2 hemit
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
                  exact ⟨howned', hremaining', hinvariants'⟩

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
