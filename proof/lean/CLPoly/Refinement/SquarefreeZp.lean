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

set_option autoImplicit false

namespace Refinement

namespace StrictSquarefreeZp

open CLPoly.Impl.StrictPolynomialGCDRefinement
open CLPoly.Math

/-- Exact checked entry for the source `__upoly_make_monic`. -/
def upolyMakeMonicIR (f : SparsePolyZp) : RawExec (Zp × SparsePolyZp) :=
  if hnonempty : 0 < f.size then
    let lc := f[0].2
    if lc.val == 1 then
      .ok (lc, f)
    else
      let lcInv := Zp.inv lc
      .ok (lc, sparseMonicLoop 0 f lcInv)
  else
    .error .assertionFailure

/-- A concrete monic sparse input takes the actual early-return comparison in
`__upoly_make_monic`; no mutation loop is erased from the definition. -/
theorem upolyMakeMonicIR_eq_of_monic (p : Nat) [Fact (Nat.Prime p)]
    (f : SparsePolyZp) (hcanonical : SparsePolyZp.Canonical p f)
    (hnonempty : 0 < f.size) (hmonic : (SparsePolyZp.toPoly p f).Monic) :
    upolyMakeMonicIR f = .ok (f[0].2, f) := by
  have hlead := sparse_leading_val_eq_one_of_monic p f hcanonical hnonempty
    hmonic
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

/-- Exact source-order range-for loop of `__extract_pth_root`. -/
def pthRootTerm (prime : UInt64) (term : UMonomial × Zp) : UMonomial × Zp :=
  (UMonomial.mk ((term.1.deg.toUInt64 / prime).toInt64), term.2)

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

def extractPthRootIR (f : SparsePolyZp) : RawExec SparsePolyZp :=
  if hnonempty : 0 < f.size then
    .ok (extractPthRootLoop 0 #[] f f[0].2.prime)
  else
    .error .assertionFailure

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

end StrictSquarefreeZp
end Refinement
