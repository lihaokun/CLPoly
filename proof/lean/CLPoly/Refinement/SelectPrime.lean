/- Genuine refinement of the generated `__select_prime` candidate pipeline. -/
import CLPoly.Generated.StrictSelectPrime
import CLPoly.Refinement.FactorZp
import CLPoly.Refinement.PolynomialMod
import Mathlib.FieldTheory.Separable

set_option autoImplicit false

open Polynomial
open CLPoly.Math

namespace Refinement.StrictSelectPrime

/-- Physical arithmetic state for one machine prime returned by the source
prime iterator.  It contains only primality/configuration/workspace data; no
field states that the input polynomial is squarefree or supplies factors. -/
structure CandidatePhysical (p : UInt64) where
  prime : Nat.Prime p.toNat
  dense : DenseUPolyZp
  primeField : dense._p = p
  configured : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured dense
  twicePrimeFits : 2 * p.toNat ≤ UInt64.size
  providers : @StrictDDF.DDFRawProviders dense ⟨primeField ▸ prime⟩

/-- Physical candidate service for a word already proved prime by the concrete
iterator.  `atPrime` is data-independent: it may allocate arithmetic buffers
but receives no polynomial and therefore cannot encode a factorization. -/
structure CandidatePhysicalProvider where
  atPrime : (p : UInt64) → Nat.Prime p.toNat → CandidatePhysical p

instance candidateFact {p : UInt64} (physical : CandidatePhysical p) :
    Fact (Nat.Prime physical.dense._p.toNat) := by
  constructor
  simpa [physical.primeField] using physical.prime

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
    (components : Array (SparsePolyZp × UInt64))
    (hready : ∀ item ∈ components.toList,
      StrictEDF.EDFEntryInvariant this item.1 item.2)
    (index : Nat) (rng : State) (result : Array SparsePolyZp) :
    let ops : Generated.StrictSelectPrime.CandidateRawOps State := {
      lcMod := fun _ _ => .error .assertionFailure
      polynomialMod := fun _ _ => .error .assertionFailure
      derivative := fun _ => .error .assertionFailure
      gcd := fun _ _ => .error .assertionFailure
      makeMonic := fun _ => .error .assertionFailure
      ddf := fun _ => .error .assertionFailure
      edf := StrictFactorZp.strictEDFCall engine this providers termination }
    ∃ output rng' factors,
      Generated.StrictSelectPrime.factorComponentsLoop ops components index
          rng result = .ok (output, rng') ∧
      factorArrayToL2 this._p.toNat output =
        factorArrayToL2 this._p.toNat result ++ factors ∧
      ∀ q ∈ factors, Irreducible q ∧ Monic q := by
  dsimp only
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
        rw [hedfRun]
        rcases ih (components.size - (index + 1)) (by omega)
            (index + 1) rngNext
            (Generated.StrictSelectPrime.appendFactorsLoop edfOutput 0 result)
            rfl with ⟨output, rng', tail, htailRun, htailDecode, htailQuality⟩
        refine ⟨output, rng', edfFactors ++ tail, htailRun, ?_, ?_⟩
        · rw [htailDecode, appendFactorsLoop_refines]
          simpa [factorArrayToL2, hedfDecode, List.append_assoc]
        · intro q hq
          rcases List.mem_append.mp hq with hq | hq
          · exact ⟨(hedfCorrect.2 q hq).1, (hedfCorrect.2 q hq).2.1⟩
          · exact htailQuality q hq
      next hindex =>
        exact ⟨result, rng, [], rfl, by simp, by simp⟩

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

/-- Mathematical payload proved for an actual successful candidate. -/
structure CandidateCorrect (f : Polynomial Int) (p : Nat)
    (factors : List (Polynomial (ZMod p))) : Prop where
  goodPrime : GoodPrime f p
  productAssociated : Associated
    (Polynomial.map (Int.castRingHom (ZMod p)) f) factors.prod
  quality : ∀ q ∈ factors, Irreducible q ∧ Monic q

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
    (fp : SparsePolyZp) (factors : Array SparsePolyZp)
    (hrun : Generated.StrictSelectPrime.tryCandidateRaw
      (strictCandidateRawOps engine physical edfTermination)
      f degF lcF p rng = .ok (.factored fp factors rng')) :
    CandidateCorrect (SparsePolyZZ.toPoly f) p.toNat
      (factorArrayToL2 p.toNat factors) := by
  let rawMod := polynomial_mod f p
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
  split at hrun <;> try contradiction
  next hlcNonzero =>
    rw [hmodRun] at hrun
    simp only at hrun
    split at hrun <;> try contradiction
    next hvalidRaw =>
      let derivative := StrictSquarefreeZp.derivativeIR physical.dense rawMod
      simp only at hrun
      split at hrun <;> try contradiction
      next hderivativeNonempty =>
        let workspace := physical.providers.gcd rawMod derivative
        have hrawNonempty : 0 < rawMod.size := by
          apply Nat.pos_of_ne_zero
          simpa [Array.isEmpty] using And.left (Bool.or_eq_false.mp hvalidRaw)
        have hderivativeCanonical : SparsePolyZp.Canonical p.toNat derivative := by
          simpa [physical.primeField] using
            StrictSquarefreeZp.derivativeIR_canonical physical.dense rawMod
              physical.configured (physical.primeField ▸ hrawCanonical)
        have hderivativeNonzero : SparsePolyZp.toPoly p.toNat derivative ≠ 0 := by
          intro hzero
          have hempty := StrictSquarefreeZp.sparsePolyZp_size_pos_of_toPoly_ne_zero
            p.toNat derivative
          have : derivative.isEmpty = false :=
            And.right (Bool.or_eq_false.mp hderivativeNonempty)
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
            workspace (physical.primeField ▸ hrawCanonical)
            hderivativeCanonical
            (by
              intro hzero
              have hpos : 0 < rawMod.size := hrawNonempty
              have hd := StrictSquarefreeZp.sparsePolyZp_toPoly_degree_eq_head
                p.toNat rawMod hrawCanonical hpos
              rw [hzero] at hd
              simp at hd)
            (physical.primeField ▸ hderivativeNonzero) with
          ⟨gcd, hgcdRun, hgcdCanonical, hgcdSemantic⟩
        rw [hgcdRun] at hrun
        simp only at hrun
        split at hrun <;> try contradiction
        next hgcdDegree =>
          rcases StrictFactorZp.upolyMakeMonicIR_refines physical.dense rawMod
              (physical.primeField ▸ hrawCanonical) hrawNonempty with
            ⟨lc, monic, hmonicRun, hmonicCanonical, hmonic,
              hreconstruct, _, hmonicDegree⟩
          rw [hmonicRun] at hrun
          simp only at hrun
          have hrawDegree : (SparsePolyZp.toPoly p.toNat rawMod).natDegree <
              2 ^ 62 := by
            rw [hrawSemantic]
            exact lt_of_le_of_lt (Polynomial.natDegree_map_le _) hdegreeBound
          have hsquarefreeRaw : Squarefree
              (SparsePolyZp.toPoly p.toNat rawMod) := by
            have hgcdNonzero := EuclideanDomain.gcd_ne_zero_of_left
              (show SparsePolyZp.toPoly p.toNat rawMod ≠ 0 from by
                intro hzero
                have hd := StrictSquarefreeZp.sparsePolyZp_toPoly_degree_eq_head
                  p.toNat rawMod hrawCanonical hrawNonempty
                rw [hzero] at hd
                simp at hd)
            have hgcdMonic : (SparsePolyZp.toPoly p.toNat gcd).Monic := by
              rw [physical.primeField] at hgcdSemantic
              rw [hgcdSemantic]
              exact Polynomial.monic_normalize hgcdNonzero
            have hgcdNatDegree :
                (SparsePolyZp.toPoly p.toNat gcd).natDegree = 0 := by
              have hdegree63 := StrictDDF.canonical_natDegree_lt_of_terms_lt
                p.toNat gcd hgcdCanonical hgcdMonic.ne_zero (2 ^ 62)
                (canonical_degree_bound p.toNat gcd hgcdCanonical
                  (by omega))
              have hnonempty := StrictSquarefreeZp.sparsePolyZp_size_pos_of_toPoly_ne_zero
                p.toNat gcd hgcdMonic.ne_zero
              have hiff := StrictDDF.strict_get_deg_pos_iff_natDegree_pos
                p.toNat gcd hgcdCanonical (lt_trans hdegree63 (by omega))
              by_contra hne
              have hpos : 0 < (SparsePolyZp.toPoly p.toNat gcd).natDegree :=
                Nat.pos_of_ne_zero hne
              have := hiff.mpr hpos
              omega
            exact squarefree_of_normalized_gcd_degree_zero
              (SparsePolyZp.toPoly p.toNat rawMod)
              (SparsePolyZp.toPoly p.toNat derivative)
              (SparsePolyZp.toPoly p.toNat gcd)
              (by
                rw [physical.primeField]
                exact StrictSquarefreeZp.derivativeIR_toPoly physical.dense
                  rawMod physical.configured (physical.primeField ▸ hrawCanonical)
                  (canonical_degree_bound p.toNat rawMod hrawCanonical hrawDegree))
              (by simpa [physical.primeField] using hgcdSemantic)
              hgcdMonic hgcdNatDegree
          have hsquarefreeMonic : Squarefree
              (SparsePolyZp.toPoly p.toNat monic) := by
            rw [physical.primeField] at hreconstruct hmonic
            exact hsquarefreeRaw.squarefree_of_dvd
              ⟨Polynomial.C (Zp.toZMod p.toNat lc), hreconstruct.symm⟩
          let ready : StrictFactorZp.DDFReady physical.dense monic :=
            ⟨canonical_prime_word_eq p monic hmonicCanonical
                (by simpa using hrawNonempty),
              hmonicCanonical,
              canonical_degree_bound p.toNat monic hmonicCanonical
                (by simpa [physical.primeField] using hmonicDegree.trans_lt hrawDegree),
              hmonic, hsquarefreeMonic⟩
          rcases StrictFactorZp.strictDDFCall_refines physical.dense
              physical.providers monic ready with
            ⟨components, hddfRun, hddfDecode, hedfReady⟩
          rw [hddfRun] at hrun
          simp only at hrun
          rcases factorComponentsLoop_refines engine physical.dense
              physical.providers edfTermination components hedfReady 0 rng #[] with
            ⟨output, rngOut, decoded, hloops, hdecode, hquality⟩
          rw [hloops] at hrun
          simp only at hrun
          have hout : output = factors ∧ rngOut = rng' := by
            exact Prod.mk.inj (Except.ok.inj hrun)
          subst output
          subst rngOut
          have hdecoded : factorArrayToL2 p.toNat factors = decoded := by
            simpa using hdecode
          have hddfCorrect := ddf_correct
            (SparsePolyZp.toPoly p.toNat monic) hmonic hsquarefreeMonic
          rw [physical.primeField] at hddfDecode
          rw [← hddfDecode] at hddfCorrect
          have hproductMonic := hmonic
          have hproductDecoded : Associated
              (SparsePolyZp.toPoly p.toNat monic) decoded.prod := by
            -- concrete DDF/EDF multiplication equality is supplied by the
            -- already proved component pipeline theorem.
            exact (ddf_edf_combine _ hmonic hsquarefreeMonic
              (fun g => StrictDDF.ddfResultToL2 p.toNat components)
              (fun _ _ => hddfCorrect)
              (fun _ _ => decoded)
              (fun _ _ _ _ _ => by
                simpa [hdecoded] using hquality)).1.associated
          refine ⟨?_, ?_, ?_⟩
          · exact ⟨physical.prime, ?_, hsquarefreeRaw⟩
            intro hlcZero
            exact hlcNonzero (by simpa [hrawSemantic] using hlcZero)
          · rw [← hrawSemantic]
            exact Associated.trans (associated_of_dvd_dvd
              ⟨_, hreconstruct⟩ ⟨_, by
                rw [hreconstruct]
                exact dvd_mul_left _ _⟩) hproductDecoded
          · simpa [hdecoded] using hquality

end Refinement.StrictSelectPrime
