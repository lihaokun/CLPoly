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

/-- Physical candidate service for every prime word visited by the concrete
iterator.  `atPrime` is data-independent: it may allocate arithmetic buffers
but receives no polynomial and therefore cannot encode a factorization. -/
structure CandidatePhysicalProvider where
  atPrime : (p : UInt64) → CandidatePhysical p

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
    StrictDDF.strictDDFGCDIR physical.dense left right
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

end Refinement.StrictSelectPrime
