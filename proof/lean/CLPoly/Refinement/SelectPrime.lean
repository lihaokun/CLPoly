/- Genuine refinement of the generated `__select_prime` candidate pipeline. -/
import CLPoly.Generated.StrictSelectPrime
import CLPoly.Generated.StrictPrimeEnumeration
import CLPoly.Refinement.FactorZp
import CLPoly.Refinement.PolynomialMod
import Mathlib.FieldTheory.Separable
import Mathlib.Tactic.NormNum.Prime

set_option autoImplicit false
set_option maxHeartbeats 0

open Polynomial
open CLPoly.Math

namespace Refinement.StrictSelectPrime

/-- Actual prime-enumeration half of the C++ `__select_prime` operation
bundle.  The candidate callback remains explicit because it also carries the
RNG and dense arithmetic workspaces for the current prime. -/
def selectPrimeRawOps {State : Type}
    (tryCandidate : SparsePolyZZ → Int64 → ZZ → UInt64 → State →
      RawExec (Generated.StrictSelectPrime.CandidateResult State)) :
    Generated.StrictSelectPrime.SelectPrimeRawOps State := {
  nextPrime := Generated.StrictPrimeEnumeration.nextPrimeRaw
  tryCandidate := tryCandidate }

def selectPrimeTermination {State : Type}
    (tryCandidate : SparsePolyZZ → Int64 → ZZ → UInt64 → State →
      RawExec (Generated.StrictSelectPrime.CandidateResult State)) :
    Generated.StrictSelectPrime.PrimeEnumerationTermination
      (selectPrimeRawOps tryCandidate) := {
  rank := Generated.StrictPrimeEnumeration.rank
  next_decreases := by
    intro useLargePrime p p' hrun
    exact Generated.StrictPrimeEnumeration.nextPrimeRaw_decreases
      useLargePrime p p' hrun }

/-- Raw storage of the C++ dense arithmetic object, indexed by its modulus.
The reconstructed object's `_p` field is definitionally `p`; refinement
proofs therefore never transport polynomials between propositionally equal
`ZMod` types. -/
structure CandidateArithmetic (p : UInt64) where
  coeffs : Array UInt64
  ninv : UInt64
  norm : UInt32

def CandidateArithmetic.dense {p : UInt64}
    (arithmetic : CandidateArithmetic p) : DenseUPolyZp :=
  { _coeffs := arithmetic.coeffs
    _p := p
    _ninv := arithmetic.ninv
    _norm := arithmetic.norm }

/-- Physical arithmetic state for one machine prime returned by the source
prime iterator.  It contains only primality/configuration/workspace data; no
field states that the input polynomial is squarefree or supplies factors. -/
structure CandidatePhysical (p : UInt64) where
  prime : Nat.Prime p.toNat
  arithmetic : CandidateArithmetic p
  configured : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured
    arithmetic.dense
  providers : @StrictDDF.DDFRawProviders arithmetic.dense ⟨prime⟩

abbrev CandidatePhysical.dense {p : UInt64}
    (physical : CandidatePhysical p) : DenseUPolyZp :=
  physical.arithmetic.dense

@[simp] theorem CandidatePhysical.dense_prime {p : UInt64}
    (physical : CandidatePhysical p) : physical.dense._p = p := rfl

/-- Physical candidate service for a word already proved prime by the concrete
iterator.  `atPrime` is data-independent: it may allocate arithmetic buffers
but receives no polynomial and therefore cannot encode a factorization. -/
structure CandidatePhysicalProvider where
  atPrime : (p : UInt64) → Nat.Prime p.toNat → CandidatePhysical p

/-- Concrete per-prime workspaces and well-founded EDF certificate used by the
runtime callback.  They are indexed only by the machine prime, never by an L2
factorization witness. -/
structure CandidateRuntimeProvider {State : Type}
    (engine : Generated.StrictEDF.RandomEngine State) where
  physical : (p : UInt64) → Nat.Prime p.toNat → CandidatePhysical p
  termination : ∀ p hp,
    let candidate := physical p hp
    letI : Fact (Nat.Prime candidate.dense._p.toNat) := ⟨candidate.prime⟩
    Generated.StrictEDF.EDFTermination
      (StrictEDF.strictEDFRawOps engine candidate.dense candidate.providers)

instance candidateFact {p : UInt64} (physical : CandidatePhysical p) :
    Fact (Nat.Prime physical.dense._p.toNat) := by
  exact ⟨physical.prime⟩

/-- Actual strict candidate operation bundle at one runtime prime. -/
noncomputable def strictCandidateRawOps {State : Type}
    (engine : Generated.StrictEDF.RandomEngine State)
    {p : UInt64} (physical : CandidatePhysical p)
    (edfTermination : Generated.StrictEDF.EDFTermination
      (StrictEDF.strictEDFRawOps engine physical.dense physical.providers)) :
    Generated.StrictSelectPrime.CandidateRawOps State := {
  lcMod := fun coefficient modulus =>
    .ok (ZZ.fdiv_r 0 coefficient modulus.toNat)
  polynomialMod := Generated.StrictPolynomialMod.polynomial_mod_raw_ir
  derivative := fun source => .ok
    (StrictSquarefreeZp.derivativeIR physical.dense source)
  gcd := fun left right =>
    StrictDDF.strictGCDIR physical.dense left right
      (physical.providers.gcd left right)
  makeMonic := StrictSquarefreeZp.upolyMakeMonicIR physical.dense
  ddf := StrictFactorZp.strictDDFCall physical.dense physical.providers
  edf := StrictFactorZp.strictEDFCall engine physical.dense
    physical.providers edfTermination }

/-- Actual runtime callback: reject a non-prime machine word and otherwise
execute the complete strict candidate body at that very word. -/
noncomputable def concreteTryCandidate {State : Type}
    (engine : Generated.StrictEDF.RandomEngine State)
    (provider : CandidateRuntimeProvider engine) :
    SparsePolyZZ → Int64 → ZZ → UInt64 → State →
      RawExec (Generated.StrictSelectPrime.CandidateResult State) :=
  fun f degF lcF p rng =>
    if hp : Nat.Prime p.toNat then
      let candidate := provider.physical p hp
      letI : Fact (Nat.Prime candidate.dense._p.toNat) := ⟨candidate.prime⟩
      Generated.StrictSelectPrime.tryCandidateRaw
        (strictCandidateRawOps engine candidate (provider.termination p hp))
        f degF lcF p rng
    else .error .assertionFailure

noncomputable def factorArrayToL2 (p : Nat) (factors : Array SparsePolyZp) :
    List (Polynomial (ZMod p)) :=
  factors.toList.map (SparsePolyZp.toPoly p)

theorem appendFactorsLoop_refines (p : Nat) (source : Array SparsePolyZp)
    (index : Nat) (result : Array SparsePolyZp) :
    factorArrayToL2 p
        (Generated.StrictSelectPrime.appendFactorsLoop source index result) =
      factorArrayToL2 p result ++
        (source.toList.drop index).map (SparsePolyZp.toPoly p) := by
  induction hmeasure : source.size - index using Nat.strong_induction_on
      generalizing index result with
  | h measure ih =>
      rw [Generated.StrictSelectPrime.appendFactorsLoop]
      split
      next hindex =>
        rw [ih (source.size - (index + 1)) (by omega)
          (index + 1) (result.push source[index]) rfl]
        have hsuffix : source.toList.drop index =
            source[index] :: source.toList.drop (index + 1) := by
          simpa using List.drop_eq_getElem_cons
            (l := source.toList) (i := index) (by simpa using hindex)
        simp [factorArrayToL2, hsuffix, List.append_assoc]
      next hindex =>
        have hdrop : source.toList.drop index = [] :=
          List.drop_eq_nil_iff.mpr (Nat.le_of_not_gt hindex)
        simp [factorArrayToL2, hdrop]

/-- The exact nested DDF-component/EDF-factor range-for returns the flattening
of the actual strict EDF calls, including the concrete final RNG state. -/
theorem factorComponentsLoop_refines {State : Type}
    (engine : Generated.StrictEDF.RandomEngine State)
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (providers : StrictDDF.DDFRawProviders this)
    (termination : Generated.StrictEDF.EDFTermination
      (StrictEDF.strictEDFRawOps engine this providers))
    (ops : Generated.StrictSelectPrime.CandidateRawOps State)
    (hedf : ops.edf = StrictFactorZp.strictEDFCall engine this providers termination)
    (components : Array (SparsePolyZp × UInt64))
    (hready : ∀ item ∈ components.toList,
      StrictEDF.EDFEntryInvariant this item.1 item.2)
    (index : Nat) (rng : State) (result : Array SparsePolyZp) :
    ∃ output rng' factors,
      Generated.StrictSelectPrime.factorComponentsLoop ops components index
          rng result = .ok (output, rng') ∧
      factorArrayToL2 this._p.toNat output =
        factorArrayToL2 this._p.toNat result ++ factors ∧
      StrictFactorZp.componentSuffixProduct this._p.toNat components index =
        factors.prod ∧
      ∀ q ∈ factors, Irreducible q ∧ Monic q := by
  induction hmeasure : components.size - index using Nat.strong_induction_on
      generalizing index rng result with
  | h measure ih =>
      rw [Generated.StrictSelectPrime.factorComponentsLoop]
      split
      next hindex =>
        let component := components[index]
        have hmem : component ∈ components.toList :=
          List.getElem_mem (by simpa using hindex)
        have hinvariant := hready component hmem
        rcases StrictFactorZp.strictEDFCall_refines engine this providers
            termination component.1 component.2 rng hinvariant with
          ⟨edfOutput, rngNext, edfFactors, hedfRun, hedfDecode, hedfCorrect⟩
        simp only
        rw [hedf]
        rw [hedfRun]
        rcases ih (components.size - (index + 1)) (by omega)
            (index + 1) rngNext
            (Generated.StrictSelectPrime.appendFactorsLoop edfOutput 0 result)
            rfl with ⟨output, rng', tail, htailRun, htailDecode,
              htailProduct, htailQuality⟩
        refine ⟨output, rng', edfFactors ++ tail, htailRun, ?_, ?_, ?_⟩
        · rw [htailDecode, appendFactorsLoop_refines]
          simpa [factorArrayToL2, hedfDecode, List.append_assoc]
        · have hsuffix : components.toList.drop index =
              component :: components.toList.drop (index + 1) := by
            simpa [component] using List.drop_eq_getElem_cons
              (l := components.toList) (i := index) (by simpa using hindex)
          rw [StrictFactorZp.componentSuffixProduct, hsuffix, List.map_cons,
            List.prod_cons, List.prod_append, ← htailProduct]
          have hcomponentEq : SparsePolyZp.toPoly this._p.toNat component.1 =
              edfFactors.prod :=
            _root_.eq_of_associated_monic _ _ hinvariant.monic
              (monic_list_prod edfFactors
                (fun q hq => (hedfCorrect.2 q hq).2.1)) hedfCorrect.1
          rw [hcomponentEq]
          simpa [StrictFactorZp.componentSuffixProduct] using
            congrArg (edfFactors.prod * ·) htailProduct
        · intro q hq
          rcases List.mem_append.mp hq with hq | hq
          · exact ⟨(hedfCorrect.2 q hq).1, (hedfCorrect.2 q hq).2.1⟩
          · exact htailQuality q hq
      next hindex =>
        have hdrop : components.toList.drop index = [] :=
          List.drop_eq_nil_iff.mpr (Nat.le_of_not_gt hindex)
        exact ⟨result, rng, [], rfl, by simp, by
          simp [StrictFactorZp.componentSuffixProduct, hdrop], by simp⟩

/-- A normalized GCD of degree zero is a unit, hence the source polynomial
and its derivative are coprime and the polynomial is squarefree. -/
theorem squarefree_of_normalized_gcd_degree_zero
    {p : Nat} [Fact (Nat.Prime p)] (f derivative gcd : Polynomial (ZMod p))
    (hderivative : derivative = Polynomial.derivative f)
    (hgcd : gcd = normalize (EuclideanDomain.gcd f derivative))
    (hgcdMonic : gcd.Monic) (hdegree : gcd.natDegree = 0) : Squarefree f := by
  have hgcdOne : gcd = 1 := by
    exact eq_one_of_monic_natDegree_zero hgcdMonic hdegree
  have hgcdUnit : IsUnit (EuclideanDomain.gcd f derivative) := by
    rw [← normalize_eq_one]
    rw [← hgcd, hgcdOne]
  have hcoprime : IsCoprime f derivative :=
    (EuclideanDomain.gcd_isUnit_iff).mp hgcdUnit
  exact (Polynomial.Separable.squarefree
    ((Polynomial.separable_def f).2 (hderivative ▸ hcoprime)))

private theorem length_le_sum_map {α : Type} (degree : α → Nat)
    (items : List α) (hpositive : ∀ item ∈ items, 1 ≤ degree item) :
    items.length ≤ (items.map degree).sum := by
  induction items with
  | nil => simp
  | cons item tail ih =>
      have hhead := hpositive item (by simp)
      have htail : ∀ x ∈ tail, 1 ≤ degree x := by
        intro x hx
        exact hpositive x (by simp [hx])
      have hi := ih htail
      simp only [List.length_cons, List.map_cons, List.sum_cons]
      omega

private theorem monic_list_prod_generic {R : Type} [Semiring R]
    (items : List (Polynomial R)) (hmonic : ∀ q ∈ items, q.Monic) :
    items.prod.Monic := by
  induction items with
  | nil => exact Polynomial.monic_one
  | cons q tail ih =>
      rw [List.prod_cons]
      exact (hmonic q (by simp)).mul (ih (by
        intro r hr
        exact hmonic r (by simp [hr])))

private theorem natDegree_list_prod_of_monic {R : Type} [Semiring R]
    [Nontrivial R] [NoZeroDivisors R]
    (items : List (Polynomial R)) (hmonic : ∀ q ∈ items, q.Monic) :
    items.prod.natDegree = (items.map Polynomial.natDegree).sum := by
  induction items with
  | nil => simp
  | cons q tail ih =>
      have hq : q.Monic := hmonic q (by simp)
      have ht : ∀ r ∈ tail, r.Monic := by
        intro r hr
        exact hmonic r (by simp [hr])
      have htprod : tail.prod.Monic := monic_list_prod_generic tail ht
      rw [List.prod_cons, Polynomial.natDegree_mul hq.ne_zero
        htprod.ne_zero, List.map_cons, List.sum_cons, ih ht]

/-- Mathematical payload proved for an actual successful candidate. -/
structure CandidateCorrect (f : Polynomial Int) (p : Nat)
    (factors : List (Polynomial (ZMod p))) : Prop where
  goodPrime : GoodPrime f p
  productAssociated : Associated
    (Polynomial.map (Int.castRingHom (ZMod p)) f) factors.prod
  quality : ∀ q ∈ factors, Irreducible q ∧ Monic q
  sizeFits : factors.length < 18446744073709551615

theorem CandidateCorrect.factors_nonempty {f : Polynomial Int} {p : Nat}
    {factors : List (Polynomial (ZMod p))}
    (correct : CandidateCorrect f p factors) (hdegree : 2 ≤ f.natDegree) :
    factors ≠ [] := by
  letI : Fact (Nat.Prime p) := ⟨correct.goodPrime.prime⟩
  intro hempty
  subst factors
  have hmapDegree :
      (Polynomial.map (Int.castRingHom (ZMod p)) f).natDegree = f.natDegree := by
    exact Polynomial.natDegree_map_of_leadingCoeff_ne_zero
      (Int.castRingHom (ZMod p)) correct.goodPrime.lc_nonzero
  have hunit : IsUnit (Polynomial.map (Int.castRingHom (ZMod p)) f) :=
    correct.productAssociated.isUnit_iff.mpr (by simp)
  have hdegreeZero :
      (Polynomial.map (Int.castRingHom (ZMod p)) f).natDegree = 0 := by
    exact Polynomial.natDegree_eq_zero_of_isUnit hunit
  omega

/-- L2 meaning of a concrete C++ prime-selection result. -/
def SelectionCorrect (f : Polynomial Int) (result : PrimeSelectionResult) : Prop :=
  CandidateCorrect f result.prime.toNat
    (factorArrayToL2 result.prime.toNat result.factors)

/-- The mutable `best` fields are either still at their exact C++ initial
state, or contain a genuinely refined candidate whose stored count agrees
with the concrete array. -/
def BestInvariant {State : Type} (f : Polynomial Int)
    (state : Generated.StrictSelectPrime.LoopState State) : Prop :=
  (state.tried = 0 ∧ state.bestCount = 18446744073709551615) ∨
  (0 < state.tried ∧ SelectionCorrect f state.best ∧
    state.bestCount = state.best.factors.size)

/-- Execution-only contract used by the outer generated loop.  It speaks
solely about successful runs of the concrete callback; the strict candidate
pipeline theorem below is what discharges it. -/
def CandidateExecutionCorrect {State : Type}
    (tryCandidate : SparsePolyZZ → Int64 → ZZ → UInt64 → State →
      RawExec (Generated.StrictSelectPrime.CandidateResult State))
    (f : SparsePolyZZ) (degF : Int64) (lcF : ZZ) : Prop :=
  ∀ p rng fp factors rng', Nat.Prime p.toNat →
    tryCandidate f degF lcF p rng = .ok (.factored fp factors rng') →
    CandidateCorrect (SparsePolyZZ.toPoly f) p.toNat
      (factorArrayToL2 p.toNat factors)

theorem selectPrimeLoop_refines {State : Type}
    (tryCandidate : SparsePolyZZ → Int64 → ZZ → UInt64 → State →
      RawExec (Generated.StrictSelectPrime.CandidateResult State))
    (f : SparsePolyZZ) (degF : Int64) (lcF : ZZ)
    (useLargePrime : Bool) (maxTries : Nat)
    (hmaxTries : 0 < maxTries)
    (hcandidate : CandidateExecutionCorrect tryCandidate f degF lcF)
    (hdegree : 2 ≤ (SparsePolyZZ.toPoly f).natDegree)
    (state : Generated.StrictSelectPrime.LoopState State)
    (hprime : Nat.Prime state.p.toNat)
    (hinvariant : BestInvariant (SparsePolyZZ.toPoly f) state)
    (result : PrimeSelectionResult)
    (hrun : Generated.StrictSelectPrime.selectPrimeLoop
      (selectPrimeRawOps tryCandidate) (selectPrimeTermination tryCandidate)
      f degF lcF useLargePrime maxTries state = .ok result) :
    SelectionCorrect (SparsePolyZZ.toPoly f) result := by
  let motive : Generated.StrictSelectPrime.LoopState State → Prop := fun state =>
    Nat.Prime state.p.toNat → BestInvariant (SparsePolyZZ.toPoly f) state →
      ∀ result, Generated.StrictSelectPrime.selectPrimeLoop
        (selectPrimeRawOps tryCandidate) (selectPrimeTermination tryCandidate)
        f degF lcF useLargePrime maxTries state = .ok result →
        SelectionCorrect (SparsePolyZZ.toPoly f) result
  have hall : ∀ state, motive state := by
    apply Generated.StrictSelectPrime.selectPrimeLoop.induct
      (ops := selectPrimeRawOps tryCandidate)
      (termination := selectPrimeTermination tryCandidate)
      (f := f) (degF := degF) (lcF := lcF)
      (useLargePrime := useLargePrime) (maxTries := maxTries)
      (motive := motive)
    · intro state hdone _ hinvariant result hrun
      rw [Generated.StrictSelectPrime.selectPrimeLoop.eq_1, dif_pos hdone] at hrun
      have hresult := Except.ok.inj hrun
      clear hrun
      cases hresult
      rcases hinvariant with hinitial | hbest
      · omega
      · simpa [SelectionCorrect] using hbest.2.1
    · intro state hnotdone fault hcandError _ _ result hrun
      rw [Generated.StrictSelectPrime.selectPrimeLoop.eq_1,
        dif_neg hnotdone, hcandError] at hrun
      contradiction
    · intro state hnotdone rng' hcand fault hnextError _ _ result hrun
      rw [Generated.StrictSelectPrime.selectPrimeLoop.eq_1,
        dif_neg hnotdone, hcand, hnextError] at hrun
      contradiction
    · intro state hnotdone rng' hcand p' hnext ih _ hinvariant result hrun
      rw [Generated.StrictSelectPrime.selectPrimeLoop.eq_1,
        dif_neg hnotdone, hcand, hnext] at hrun
      exact ih (Generated.StrictPrimeEnumeration.nextPrimeRaw_prime
        useLargePrime state.p p' hnext) (by simpa [BestInvariant] using hinvariant)
        result hrun
    · intro state hnotdone fp factors rng' hcand hsmall hprime _ result hrun
      have hcorrect := hcandidate state.p state.rng fp factors rng' hprime hcand
      have hnonempty := hcorrect.factors_nonempty hdegree
      have harrayNonempty : factors.isEmpty = false := by
        apply Bool.eq_false_of_not_eq_true
        intro hempty
        apply hnonempty
        simpa [factorArrayToL2, Array.isEmpty] using hempty
      simp [Generated.StrictSelectPrime.selectPrimeLoop.eq_1,
        hnotdone, hcand, hsmall, harrayNonempty] at hrun
      subst result
      simpa [SelectionCorrect] using hcorrect
    · intro state hnotdone fp factors rng' hcand hlarge fault hnextError
        _ _ result hrun
      rw [Generated.StrictSelectPrime.selectPrimeLoop.eq_1,
        dif_neg hnotdone] at hrun
      simp only [hcand] at hrun
      rw [if_neg hlarge] at hrun
      split at hrun
      next fault' hcase =>
        rw [hnextError] at hcase
        cases hrun
      next p'' hcase =>
        rw [hnextError] at hcase
        contradiction
    · intro state hnotdone fp factors rng' hcand hlarge best bestCount p' hnext ih
        hprime hinvariant result hrun
      have hcorrect := hcandidate state.p state.rng fp factors rng' hprime hcand
      rw [Generated.StrictSelectPrime.selectPrimeLoop.eq_1,
        dif_neg hnotdone] at hrun
      simp only [hcand] at hrun
      rw [if_neg hlarge] at hrun
      split at hrun
      next fault hcase =>
        rw [hnext] at hcase
        contradiction
      next p'' hcase =>
        have hp' : p'' = p' := Except.ok.inj (hcase.symm.trans hnext)
        subst p''
        apply ih (Generated.StrictPrimeEnumeration.nextPrimeRaw_prime
          useLargePrime state.p p' hnext)
        · dsimp [best, bestCount]
          unfold BestInvariant
          right
          dsimp only
          refine ⟨by omega, ?_, ?_⟩
          · dsimp [best]
            split
            next hbetter => simpa [SelectionCorrect] using hcorrect
            next hnotBetter =>
              rcases hinvariant with hinitial | hbest
              · exact False.elim (hnotBetter (by
                  rw [hinitial.2]
                  simpa [factorArrayToL2] using hcorrect.sizeFits))
              · exact hbest.2.1
          · dsimp [best, bestCount]
            split
            next hbetter => simp [Nat.min_eq_right (Nat.le_of_lt hbetter)]
            next hnotBetter =>
              rcases hinvariant with hinitial | hbest
              · exact False.elim (hnotBetter (by
                  rw [hinitial.2]
                  simpa [factorArrayToL2] using hcorrect.sizeFits))
              · rw [hbest.2.2]
                exact Nat.min_eq_left (Nat.le_of_not_gt
                  (hbest.2.2 ▸ hnotBetter))
        · exact hrun
  exact hall state hprime hinvariant result hrun

private theorem canonical_prime_word_eq
    (p : UInt64) (f : SparsePolyZp)
    (hcanonical : SparsePolyZp.Canonical p.toNat f)
    (hnonempty : 0 < f.size) : f[0]!.2.prime = p := by
  have hmem : f[0] ∈ f.toList := by
    simpa using Array.getElem_mem f 0 hnonempty
  have hnat := (hcanonical.1 f[0] hmem).1
  rw [getElem!_pos f 0 hnonempty]
  exact UInt64.toNat_inj.mp hnat

private theorem canonical_degree_bound
    (p : Nat) [Fact (Nat.Prime p)] (f : SparsePolyZp)
    (hcanonical : SparsePolyZp.Canonical p f)
    (hbound : (SparsePolyZp.toPoly p f).natDegree < 2 ^ 62) :
    ∀ term ∈ f.toList, term.1.deg < 2 ^ 62 := by
  intro term hterm
  exact lt_of_le_of_lt
    (StrictDDF.canonical_term_degree_le_natDegree p f hcanonical term hterm)
    hbound

/-- A successful execution of the fully expanded candidate body establishes
the good-prime and irreducible-factor facts used by prime selection. -/
theorem tryCandidateRaw_factored_refines {State : Type}
    (engine : Generated.StrictEDF.RandomEngine State)
    {p : UInt64} (physical : CandidatePhysical p)
    (edfTermination : Generated.StrictEDF.EDFTermination
      (StrictEDF.strictEDFRawOps engine physical.dense physical.providers))
    (f : SparsePolyZZ) (degF : Int64) (lcF : ZZ) (rng rng' : State)
    (hcanonicalZZ : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hdegreeBound : (SparsePolyZZ.toPoly f).natDegree < 2 ^ 62)
    (hlcSemantic : (lcF : ZMod p.toNat) =
      ((SparsePolyZZ.toPoly f).leadingCoeff : ZMod p.toNat))
    (fp : SparsePolyZp) (factors : Array SparsePolyZp)
    (hrun : Generated.StrictSelectPrime.tryCandidateRaw
      (strictCandidateRawOps engine physical edfTermination)
      f degF lcF p rng = .ok (.factored fp factors rng')) :
    CandidateCorrect (SparsePolyZZ.toPoly f) p.toNat
      (factorArrayToL2 p.toNat factors) := by
  let rawMod := polynomial_mod f p
  letI : Fact (Nat.Prime p.toNat) := ⟨physical.prime⟩
  have hp : 0 < p.toNat := physical.prime.pos
  have hmodRun : Generated.StrictPolynomialMod.polynomial_mod_raw_ir f p =
      .ok rawMod := by
    simp [Generated.StrictPolynomialMod.polynomial_mod_raw_ir,
      rawMod, StrictPolynomialMod.polynomialModLoop_eq_model]
  have hrawCanonical : SparsePolyZp.Canonical p.toNat rawMod :=
    StrictPolynomialMod.polynomial_mod_canonical f p hp hcanonicalZZ
  have hrawSemantic : SparsePolyZp.toPoly p.toNat rawMod =
      Polynomial.map (Int.castRingHom (ZMod p.toNat))
        (SparsePolyZZ.toPoly f) :=
    StrictPolynomialMod.polynomial_mod_toPoly f p hp
  unfold Generated.StrictSelectPrime.tryCandidateRaw at hrun
  simp only [strictCandidateRawOps] at hrun
  split at hrun
  · simp_all
  next hlcNonzero =>
    rw [hmodRun] at hrun
    simp only at hrun
    split at hrun
    · simp_all
    next hvalidRaw =>
      let derivative := StrictSquarefreeZp.derivativeIR physical.dense rawMod
      split at hrun
      · simp_all
      next hderivativeNonempty =>
        let workspace := physical.providers.gcd rawMod derivative
        have hrawNonempty : 0 < rawMod.size := by
          apply Nat.pos_of_ne_zero
          have hguard : (rawMod.isEmpty || get_deg rawMod != degF) = false := by
            simpa [Bool.not_eq_true] using hvalidRaw
          have hempty : rawMod.isEmpty = false := by
            exact (Bool.or_eq_false_iff.mp hguard).1
          simpa [Array.isEmpty] using hempty
        have hderivativeCanonical : SparsePolyZp.Canonical p.toNat derivative := by
          exact StrictSquarefreeZp.derivativeIR_canonical physical.dense rawMod
            physical.configured hrawCanonical
        have hrawCanonicalDense :
            SparsePolyZp.Canonical physical.dense._p.toNat rawMod := by
          exact hrawCanonical
        have hderivativeCanonicalDense :
            SparsePolyZp.Canonical physical.dense._p.toNat derivative := by
          exact hderivativeCanonical
        have hderivativeNonzero : SparsePolyZp.toPoly p.toNat derivative ≠ 0 := by
          intro hzero
          have hempty := StrictSquarefreeZp.sparsePolyZp_size_pos_of_toPoly_ne_zero
            p.toNat derivative
          have : derivative.isEmpty = false :=
            Bool.eq_false_of_not_eq_true hderivativeNonempty
          have hpos : 0 < derivative.size := by
            apply Nat.pos_of_ne_zero
            simpa [Array.isEmpty] using this
          exact (Nat.ne_of_gt hpos) (by
            have := hderivativeCanonical
            simpa [hzero, SparsePolyZp.toPoly] using
              (show derivative.size = 0 from by
                by_contra hne
                have hd := StrictSquarefreeZp.sparsePolyZp_toPoly_degree_eq_head
                  p.toNat derivative hderivativeCanonical (Nat.pos_of_ne_zero hne)
                rw [hzero] at hd
                simp at hd))
        rcases StrictDDF.strictGCDIR_refines physical.dense rawMod derivative
            workspace hrawCanonicalDense hderivativeCanonicalDense
            (by
              simpa using (show SparsePolyZp.toPoly p.toNat rawMod ≠ 0 from by
                intro hzero
                have hd := StrictSquarefreeZp.sparsePolyZp_toPoly_degree_eq_head
                  p.toNat rawMod hrawCanonical hrawNonempty
                rw [hzero] at hd
                simp at hd))
            hderivativeNonzero with
          ⟨gcd, hgcdRun, hgcdCanonical, hgcdSemantic⟩
        rw [hgcdRun] at hrun
        simp only at hrun
        split at hrun
        · simp_all
        next hgcdDegree =>
          rcases StrictFactorZp.upolyMakeMonicIR_refines physical.dense rawMod
              hrawCanonical hrawNonempty with
            ⟨lc, monic, hmonicRun, hmonicCanonical, hmonic,
              hreconstruct, hmonicSize, hmonicDegree⟩
          rw [hmonicRun] at hrun
          simp only at hrun
          have hrawDegree : (SparsePolyZp.toPoly p.toNat rawMod).natDegree <
              2 ^ 62 := by
            rw [hrawSemantic]
            exact lt_of_le_of_lt Polynomial.natDegree_map_le hdegreeBound
          have hrawNonzero : SparsePolyZp.toPoly p.toNat rawMod ≠ 0 := by
            intro hzero
            have hd := StrictSquarefreeZp.sparsePolyZp_toPoly_degree_eq_head
              p.toNat rawMod hrawCanonical hrawNonempty
            rw [hzero] at hd
            simp at hd
          have hsquarefreeRaw : Squarefree
              (SparsePolyZp.toPoly p.toNat rawMod) := by
            have hgcdNonzero : EuclideanDomain.gcd
                (SparsePolyZp.toPoly p.toNat rawMod)
                (SparsePolyZp.toPoly p.toNat derivative) ≠ 0 := by
              intro hzero
              exact hrawNonzero (EuclideanDomain.gcd_eq_zero_iff.mp hzero).1
            have hgcdMonic : (SparsePolyZp.toPoly p.toNat gcd).Monic := by
              rw [show SparsePolyZp.toPoly p.toNat gcd =
                  normalize (EuclideanDomain.gcd
                    (SparsePolyZp.toPoly p.toNat rawMod)
                    (SparsePolyZp.toPoly p.toNat derivative)) by
                simpa using hgcdSemantic]
              exact Polynomial.monic_normalize hgcdNonzero
            have hgcdNatDegree :
                (SparsePolyZp.toPoly p.toNat gcd).natDegree = 0 := by
              have hgcdCanonicalP : SparsePolyZp.Canonical p.toNat gcd := by
                simpa using hgcdCanonical
              have hgcdDvd : SparsePolyZp.toPoly p.toNat gcd ∣
                  SparsePolyZp.toPoly p.toNat rawMod := by
                rw [show SparsePolyZp.toPoly p.toNat gcd =
                    normalize (EuclideanDomain.gcd
                      (SparsePolyZp.toPoly p.toNat rawMod)
                      (SparsePolyZp.toPoly p.toNat derivative)) by
                  simpa using hgcdSemantic, normalize_dvd_iff]
                exact EuclideanDomain.gcd_dvd_left _ _
              have hgcdDegreeLe := Polynomial.natDegree_le_of_dvd hgcdDvd
                hrawNonzero
              have hdegree63 := StrictDDF.canonical_natDegree_lt_of_terms_lt
                p.toNat gcd hgcdCanonicalP hgcdMonic.ne_zero (2 ^ 62)
                (canonical_degree_bound p.toNat gcd hgcdCanonicalP
                  (lt_of_le_of_lt hgcdDegreeLe hrawDegree))
              have hnonempty := StrictSquarefreeZp.sparsePolyZp_size_pos_of_toPoly_ne_zero
                p.toNat gcd hgcdMonic.ne_zero
              have hiff := StrictDDF.strict_get_deg_pos_iff_natDegree_pos
                p.toNat gcd hgcdCanonicalP (lt_trans hdegree63 (by omega))
              by_contra hne
              have hpos : 0 < (SparsePolyZp.toPoly p.toNat gcd).natDegree :=
                Nat.pos_of_ne_zero hne
              exact hgcdDegree (hiff.mpr hpos)
            exact squarefree_of_normalized_gcd_degree_zero
              (SparsePolyZp.toPoly p.toNat rawMod)
              (SparsePolyZp.toPoly p.toNat derivative)
              (SparsePolyZp.toPoly p.toNat gcd)
              (StrictSquarefreeZp.derivativeIR_toPoly physical.dense
                rawMod physical.configured hrawCanonical
                (fun term hterm => lt_trans
                  (canonical_degree_bound p.toNat rawMod hrawCanonical hrawDegree
                    term hterm) (by norm_num [UInt64.size])))
              (by simpa using hgcdSemantic)
              hgcdMonic hgcdNatDegree
          have hsquarefreeMonic : Squarefree
              (SparsePolyZp.toPoly physical.dense._p.toNat monic) := by
            have hsquarefreeRawDense : Squarefree
                (SparsePolyZp.toPoly physical.dense._p.toNat rawMod) := by
              exact hsquarefreeRaw
            exact hsquarefreeRawDense.squarefree_of_dvd
              ⟨Polynomial.C (Zp.toZMod physical.dense._p.toNat lc), by
                rw [mul_comm]
                exact hreconstruct⟩
          have hrawDegreeDense :
              (SparsePolyZp.toPoly physical.dense._p.toNat rawMod).natDegree <
                2 ^ 62 := by
            exact hrawDegree
          let ready : StrictFactorZp.DDFReady physical.dense monic :=
            ⟨canonical_prime_word_eq physical.dense._p monic hmonicCanonical
                (by simpa [hmonicSize] using hrawNonempty),
              hmonicCanonical,
              canonical_degree_bound physical.dense._p.toNat monic hmonicCanonical
                (hmonicDegree.trans_lt hrawDegreeDense),
              hmonic, hsquarefreeMonic⟩
          rcases StrictFactorZp.strictDDFCall_refines physical.dense
              physical.providers monic ready with
            ⟨components, hddfRun, hddfDecode, hedfReady⟩
          rw [hddfRun] at hrun
          simp only at hrun
          rcases factorComponentsLoop_refines engine physical.dense
              physical.providers edfTermination
              (strictCandidateRawOps engine physical edfTermination) rfl
              components hedfReady 0 rng #[] with
            ⟨output, rngOut, decoded, hloops, hdecode, hcomponentProduct,
              hquality⟩
          have hloops' := hloops
          simp only [strictCandidateRawOps] at hloops'
          rw [hloops'] at hrun
          simp only at hrun
          have hout : monic = fp ∧ output = factors ∧ rngOut = rng' := by
            exact Generated.StrictSelectPrime.CandidateResult.factored.inj
              (Except.ok.inj hrun)
          rcases hout with ⟨rfl, rfl, rfl⟩
          have hdecoded : factorArrayToL2 p.toNat output = decoded := by
            simpa using hdecode
          have hddfCorrect := ddf_correct
            (SparsePolyZp.toPoly p.toNat monic) hmonic hsquarefreeMonic
          have hddfDecodeP : StrictDDF.ddfResultToL2 p.toNat components =
              ddf (SparsePolyZp.toPoly p.toNat monic) := by
            simpa using hddfDecode
          rw [← hddfDecodeP] at hddfCorrect
          have hproductMonic := hmonic
          have hproductDecoded : Associated
              (SparsePolyZp.toPoly p.toNat monic) decoded.prod := by
            have hddfProduct := hddfCorrect.2.2.2.1
            have hcomponentProductP :
                StrictFactorZp.componentSuffixProduct p.toNat components 0 =
                  decoded.prod := by simpa using hcomponentProduct
            have hcomponentAll :
                StrictFactorZp.componentSuffixProduct p.toNat components 0 =
                  ((StrictDDF.ddfResultToL2 p.toNat components).map Prod.fst).prod := by
              simp [StrictFactorZp.componentSuffixProduct,
                StrictDDF.ddfResultToL2, Function.comp_def]
            exact Associated.trans hddfProduct (by
              rw [← hcomponentAll, hcomponentProductP])
          refine ⟨?_, ?_, ?_, ?_⟩
          · have hlcModNonzero : (lcF : ZMod p.toNat) ≠ 0 := by
              intro hlcZero
              apply hlcNonzero
              have hdiv : (p.toNat : Int) ∣ lcF :=
                (ZMod.intCast_zmod_eq_zero_iff_dvd lcF p.toNat).mp hlcZero
              have hfmod : Int.fmod lcF p.toNat = 0 := by
                rw [Int.fmod_eq_emod]
                simp [Int.emod_eq_zero_of_dvd hdiv]
              simp [ZZ.fdiv_r, hfmod, ZZ.toBool]
            exact ⟨physical.prime, hlcSemantic ▸ hlcModNonzero, by
              rw [← hrawSemantic]
              exact hsquarefreeRaw⟩
          · rw [← hrawSemantic]
            have hreconstructP : SparsePolyZp.toPoly p.toNat rawMod =
                Polynomial.C (Zp.toZMod p.toNat lc) *
                  SparsePolyZp.toPoly p.toNat monic := by
              simpa using hreconstruct
            have hrawMonicAssociated : Associated
                (SparsePolyZp.toPoly p.toNat rawMod)
                (SparsePolyZp.toPoly p.toNat monic) := by
              have hlcField : Zp.toZMod p.toNat lc ≠ 0 := by
                intro hlcZero
                simp [hlcZero] at hreconstructP
                exact hrawNonzero hreconstructP
              rw [hreconstructP, mul_comm]
              exact associated_mul_unit_left
                (SparsePolyZp.toPoly p.toNat monic)
                (Polynomial.C (Zp.toZMod p.toNat lc)) (by
                  simpa [Polynomial.isUnit_C] using hlcField)
            exact Associated.trans hrawMonicAssociated (by
              rw [hdecoded]
              exact hproductDecoded)
          · simpa [hdecoded] using hquality
          · have hfactorMonic : ∀ q ∈ decoded, q.Monic :=
              fun q hq => (hquality q hq).2
            have hfactorDegree : ∀ q ∈ decoded, 1 ≤ q.natDegree := by
              intro q hq
              have hirr := (hquality q hq).1
              by_contra hzero
              have hdegree : q.natDegree = 0 := by omega
              have hqOne : q = 1 := eq_one_of_monic_natDegree_zero
                (hfactorMonic q hq) hdegree
              exact hirr.not_isUnit (hqOne ▸ isUnit_one)
            have hprodMonic := monic_list_prod decoded hfactorMonic
            have hproductEq : SparsePolyZp.toPoly p.toNat monic = decoded.prod :=
              eq_of_associated_monic _ _ hmonic hprodMonic hproductDecoded
            have hdegreeSum : decoded.prod.natDegree =
                (decoded.map Polynomial.natDegree).sum := by
              exact natDegree_list_prod_of_monic decoded hfactorMonic
            have hlengthLe : decoded.length ≤
                (decoded.map Polynomial.natDegree).sum := by
              exact length_le_sum_map Polynomial.natDegree decoded hfactorDegree
            rw [hdecoded]
            exact lt_of_le_of_lt hlengthLe (by
              rw [← hdegreeSum, ← hproductEq]
              have hmonicDegreeP :
                  (SparsePolyZp.toPoly p.toNat monic).natDegree =
                    (SparsePolyZp.toPoly p.toNat rawMod).natDegree := by
                simpa using hmonicDegree
              exact lt_trans (hmonicDegreeP.trans_lt hrawDegree) (by norm_num))

/-- The concrete runtime callback satisfies the outer loop's execution
contract by running the already-refined strict candidate body at the prime
proved by the generated iterator. -/
theorem concreteTryCandidate_correct {State : Type}
    (engine : Generated.StrictEDF.RandomEngine State)
    (provider : CandidateRuntimeProvider engine)
    (f : SparsePolyZZ) (degF : Int64) (lcF : ZZ)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hdegreeBound : (SparsePolyZZ.toPoly f).natDegree < 2 ^ 62)
    (hlcSemantic : ∀ p : UInt64, Nat.Prime p.toNat →
      (lcF : ZMod p.toNat) =
        ((SparsePolyZZ.toPoly f).leadingCoeff : ZMod p.toNat)) :
    CandidateExecutionCorrect (concreteTryCandidate engine provider)
      f degF lcF := by
  intro p rng fp factors rng' hp hrun
  unfold concreteTryCandidate at hrun
  simp only [dif_pos hp] at hrun
  let candidate := provider.physical p hp
  letI : Fact (Nat.Prime candidate.dense._p.toNat) := ⟨candidate.prime⟩
  exact tryCandidateRaw_factored_refines engine candidate
    (provider.termination p hp) f degF lcF rng rng' hcanonical
    hdegreeBound (hlcSemantic p hp) fp factors hrun

/-- Genuine L1→L2 refinement of the generated original C++
`__select_prime` entry, instantiated with the concrete modular
reduction/GCD/DDF/EDF callback. -/
theorem __select_prime_raw_ir_refines {State : Type}
    (engine : Generated.StrictEDF.RandomEngine State)
    (provider : CandidateRuntimeProvider engine) (initialRng : State)
    (useLargePrime : Bool) (f : SparsePolyZZ)
    (hinitialPrimeCorrect : Nat.Prime
      (if useLargePrime then
        ((18446744073709551615 : UInt64) - 58).toNat
      else (2 : UInt64).toNat))
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hnonempty : 0 < f.size)
    (hdegree : 2 ≤ (SparsePolyZZ.toPoly f).natDegree)
    (hdegreeBound : (SparsePolyZZ.toPoly f).natDegree < 2 ^ 62)
    (hlcSemantic : ∀ p : UInt64, Nat.Prime p.toNat →
      ((SparsePolyZZ.front! f).2 : ZMod p.toNat) =
        ((SparsePolyZZ.toPoly f).leadingCoeff : ZMod p.toNat))
    (result : PrimeSelectionResult)
    (hrun : Generated.StrictSelectPrime.__select_prime_raw_ir
      (selectPrimeRawOps (concreteTryCandidate engine provider))
      (selectPrimeTermination (concreteTryCandidate engine provider))
      initialRng useLargePrime f = .ok result) :
    SelectionCorrect (SparsePolyZZ.toPoly f) result := by
  unfold Generated.StrictSelectPrime.__select_prime_raw_ir at hrun
  split at hrun
  next hguard => contradiction
  next hguard =>
    let initialPrime : UInt64 :=
      if useLargePrime then (18446744073709551615 : UInt64) - 58 else 2
    have hinitialPrime : Nat.Prime initialPrime.toNat := by
      cases useLargePrime <;> simpa [initialPrime] using hinitialPrimeCorrect
    apply selectPrimeLoop_refines
      (tryCandidate := concreteTryCandidate engine provider)
      f (get_deg f) (SparsePolyZZ.front! f).2 useLargePrime 3 (by omega)
      (concreteTryCandidate_correct engine provider f (get_deg f)
        (SparsePolyZZ.front! f).2 hcanonical hdegreeBound hlcSemantic)
      hdegree
      { tried := 0, p := initialPrime, rng := initialRng,
        bestCount := 18446744073709551615, best := default }
      hinitialPrime (by left; exact ⟨rfl, rfl⟩) result
    simpa [initialPrime] using hrun

end Refinement.StrictSelectPrime
