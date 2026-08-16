/-
  Strict end-to-end refinement assembly for the original C++
  `__factor_squarefree_primitive_ZZ` entry.

  This module is intentionally downstream of both SelectPrime and Hensel.  It
  contains cross-stage composition facts; neither component module imports the
  other or assumes a semantic result produced by a later stage.
-/
import CLPoly.Generated.StrictFactorZZ
import CLPoly.Refinement.Hensel
import CLPoly.Refinement.Recombine
import CLPoly.Refinement.SelectPrime

set_option autoImplicit false

namespace Refinement

namespace StrictFactorZZ

open CLPoly.Math

private theorem polynomial_eq_C_leadingCoeff_mul_of_associated_monic
    {K : Type*} [Field K] {f g : Polynomial K}
    (hf : f ≠ 0) (hg : g.Monic) (hassociated : Associated f g) :
    f = Polynomial.C f.leadingCoeff * g := by
  rcases hassociated with ⟨unit, hunit⟩
  rcases Polynomial.isUnit_iff.mp unit.isUnit with
    ⟨coefficient, _hcoefficientUnit, hcoefficient⟩
  have hproduct : f * Polynomial.C coefficient = g := by
    rw [hcoefficient]
    exact hunit
  have hleading := congrArg Polynomial.leadingCoeff hproduct
  rw [Polynomial.leadingCoeff_mul, Polynomial.leadingCoeff_C,
    hg.leadingCoeff] at hleading
  symm
  rw [← hproduct]
  calc
    Polynomial.C f.leadingCoeff *
          (f * Polynomial.C coefficient) =
        f * Polynomial.C (f.leadingCoeff * coefficient) := by
      rw [Polynomial.C_mul]
      ring
    _ = f := by rw [hleading]; simp

/-- Resolve the unit ambiguity of a modular subproduct using the leading
coefficient that the concrete Zassenhaus trial places in front of it.

If `source = divisor * quotient`, the monic modular representative associated
with `divisor` becomes exactly `C quotient.leadingCoeff * divisor` after it is
scaled by `source.leadingCoeff`.  This is the normalization needed by the
generated trial product; keeping the quotient's leading coefficient here is
essential for non-monic integer factors. -/
theorem leading_scaled_monic_associated_divisor
    (p : Nat) [Fact (Nat.Prime p)]
    (source divisor quotient : Polynomial Int)
    (candidate : Polynomial (ZMod p))
    (hfactor : source = divisor * quotient)
    (hdivisorModNonzero : Polynomial.map
      (Int.castRingHom (ZMod p)) divisor ≠ 0)
    (hdivisorLeading :
      (Polynomial.map (Int.castRingHom (ZMod p)) divisor).leadingCoeff =
        (divisor.leadingCoeff : ZMod p))
    (hcandidateMonic : candidate.Monic)
    (hassociated : Associated
      (Polynomial.map (Int.castRingHom (ZMod p)) divisor) candidate) :
    Polynomial.C (source.leadingCoeff : ZMod p) * candidate =
      Polynomial.map (Int.castRingHom (ZMod p))
        (Polynomial.C quotient.leadingCoeff * divisor) := by
  have hexact := polynomial_eq_C_leadingCoeff_mul_of_associated_monic
    hdivisorModNonzero hcandidateMonic hassociated
  rw [hdivisorLeading] at hexact
  have hsourceLeading :
      source.leadingCoeff = divisor.leadingCoeff * quotient.leadingCoeff := by
    rw [hfactor, Polynomial.leadingCoeff_mul]
  rw [hsourceLeading, Int.cast_mul, Polynomial.C_mul,
    Polynomial.map_mul, Polynomial.map_C]
  rw [hexact]
  change Polynomial.C (divisor.leadingCoeff : ZMod p) *
      Polynomial.C (quotient.leadingCoeff : ZMod p) * candidate =
    Polynomial.C (quotient.leadingCoeff : ZMod p) *
      (Polynomial.C (divisor.leadingCoeff : ZMod p) * candidate)
  ring

private theorem monic_list_product {K : Type*} [Semiring K]
    (factors : List (Polynomial K))
    (hmonic : ∀ factor ∈ factors, factor.Monic) : factors.prod.Monic := by
  induction factors with
  | nil => exact Polynomial.monic_one
  | cons factor tail ih =>
      exact (hmonic factor (by simp)).mul
        (ih fun candidate hcandidate => hmonic candidate (by simp [hcandidate]))

/-- Every divisor of a concrete finite product of irreducibles is associated
to the product of an actual sublist.  The proof recursively follows the
source-order list: if the head atom divides the divisor it is cancelled from
both sides; otherwise irreducibility makes it coprime to the divisor and it
is cancelled from the ambient product. -/
theorem divisor_associated_sublist_product
    {R : Type*} [CommRing R] [IsDomain R] [IsBezout R]
    (atoms : List R) (hirreducible : ∀ atom ∈ atoms, Irreducible atom)
    (divisor : R) (hdivisor : divisor ∣ atoms.prod) :
    ∃ chosen : List R, List.Sublist chosen atoms ∧
      Associated divisor chosen.prod := by
  induction atoms generalizing divisor with
  | nil =>
      refine ⟨[], List.Sublist.refl [], associated_one_iff_isUnit.mpr ?_⟩
      simpa using isUnit_of_dvd_one hdivisor
  | cons atom atoms ih =>
      have hatom := hirreducible atom (by simp)
      have htail : ∀ candidate ∈ atoms, Irreducible candidate := by
        intro candidate hcandidate
        exact hirreducible candidate (by simp [hcandidate])
      by_cases hheadDivisor : atom ∣ divisor
      · rcases hheadDivisor with ⟨quotient, hquotient⟩
        subst divisor
        have hquotientDivides : quotient ∣ atoms.prod := by
          exact (mul_dvd_mul_iff_left hatom.ne_zero).mp (by
            simpa [List.prod_cons] using hdivisor)
        rcases ih htail quotient hquotientDivides with
          ⟨chosen, hchosen, hassociated⟩
        refine ⟨atom :: chosen, hchosen.cons_cons atom, ?_⟩
        simpa [List.prod_cons] using
          Associated.mul_mul (Associated.refl atom) hassociated
      · have hcoprime : IsCoprime divisor atom :=
          (hatom.coprime_iff_not_dvd.mpr hheadDivisor).symm
        have hdivisorTail : divisor ∣ atoms.prod := by
          exact hcoprime.dvd_of_dvd_mul_left (by
            simpa [List.prod_cons] using hdivisor)
        rcases ih htail divisor hdivisorTail with
          ⟨chosen, hchosen, hassociated⟩
        exact ⟨chosen, hchosen.cons atom, hassociated⟩

/-- A primitive integer polynomial is irreducible when its reduction at a
prime is irreducible and the leading coefficient survives reduction.  This is
the non-monic form needed for the actual primitive factors emitted by C++
recombination.  The leading-coefficient hypothesis prevents a factor from
becoming a constant merely through degree loss modulo `p`. -/
theorem primitive_irreducible_of_irreducible_mod
    (p : Nat) [Fact (Nat.Prime p)] (g : Polynomial Int)
    (hprimitive : g.IsPrimitive)
    (hleading : (g.leadingCoeff : ZMod p) ≠ 0)
    (hmodIrreducible :
      Irreducible (g.map (Int.castRingHom (ZMod p)))) :
    Irreducible g := by
  refine ⟨fun hunit => hmodIrreducible.not_isUnit
      (hunit.map (Polynomial.mapRingHom (Int.castRingHom (ZMod p)))), ?_⟩
  intro left right hfactor
  have hleadingProduct :
      (left.leadingCoeff : ZMod p) * (right.leadingCoeff : ZMod p) =
        (g.leadingCoeff : ZMod p) := by
    rw [← Int.cast_mul]
    congr 1
    rw [← Polynomial.leadingCoeff_mul, hfactor]
  have hleftLeading : (left.leadingCoeff : ZMod p) ≠ 0 := by
    intro hzero
    rw [hzero, zero_mul] at hleadingProduct
    exact hleading hleadingProduct.symm
  have hrightLeading : (right.leadingCoeff : ZMod p) ≠ 0 := by
    intro hzero
    rw [hzero, mul_zero] at hleadingProduct
    exact hleading hleadingProduct.symm
  have hmapped :
      g.map (Int.castRingHom (ZMod p)) =
        left.map (Int.castRingHom (ZMod p)) *
          right.map (Int.castRingHom (ZMod p)) := by
    rw [← Polynomial.map_mul, hfactor]
  rcases hmodIrreducible.isUnit_or_isUnit hmapped with hleft | hright
  · have hleftPrimitive : left.IsPrimitive :=
      Polynomial.isPrimitive_of_dvd hprimitive
        (Dvd.intro right hfactor.symm)
    have hdegreeMapped :
        (left.map (Int.castRingHom (ZMod p))).degree = left.degree :=
      Polynomial.degree_map_eq_of_leadingCoeff_ne_zero _ hleftLeading
    have hdegreeZero : left.degree = 0 := by
      rw [← hdegreeMapped]
      exact Polynomial.isUnit_iff_degree_eq_zero.mp hleft
    rw [Polynomial.eq_C_of_degree_eq_zero hdegreeZero] at hleftPrimitive ⊢
    exact Or.inl (Polynomial.isUnit_C.mpr
      (Polynomial.isPrimitive_iff_isUnit_of_C_dvd.mp hleftPrimitive
        (left.coeff 0) dvd_rfl))
  · have hrightPrimitive : right.IsPrimitive :=
      Polynomial.isPrimitive_of_dvd hprimitive
        (Dvd.intro_left left hfactor.symm)
    have hdegreeMapped :
        (right.map (Int.castRingHom (ZMod p))).degree = right.degree :=
      Polynomial.degree_map_eq_of_leadingCoeff_ne_zero _ hrightLeading
    have hdegreeZero : right.degree = 0 := by
      rw [← hdegreeMapped]
      exact Polynomial.isUnit_iff_degree_eq_zero.mp hright
    rw [Polynomial.eq_C_of_degree_eq_zero hdegreeZero] at hrightPrimitive ⊢
    exact Or.inr (Polynomial.isUnit_C.mpr
      (Polynomial.isPrimitive_iff_isUnit_of_C_dvd.mp hrightPrimitive
        (right.coeff 0) dvd_rfl))

/-- The sparse factor array returned by a genuinely refined successful prime
selection consists pointwise of the irreducible factors recorded by that
selection's concrete candidate certificate. -/
theorem selectionFactors_irreducible
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    (hselection : StrictSelectPrime.SelectionCorrect
      (SparsePolyZZ.toPoly f) selection) :
    ∀ index (hindex : index < selection.factors.size),
      Irreducible (SparsePolyZp.toPoly selection.prime.toNat
        selection.factors[index]) := by
  intro index hindex
  have hmember : SparsePolyZp.toPoly selection.prime.toNat
      selection.factors[index] ∈
        StrictSelectPrime.factorArrayToL2 selection.prime.toNat
          selection.factors := by
    simp only [StrictSelectPrime.factorArrayToL2, List.mem_map]
    exact ⟨selection.factors[index], Array.getElem_mem_toList hindex, rfl⟩
  exact (hselection.quality _ hmember).1

/-- The exact first-factor adjustment performed by Hensel preserves the
pointwise irreducibility supplied by the actual SelectPrime result. -/
theorem selectionFactors_adjusted_irreducible
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {adjusted : Array SparsePolyZp}
    (hp2 : selection.prime.toNat * selection.prime.toNat ≤ UInt64.size)
    (hfactors : ∀ factor ∈ selection.factors.toList,
      SparsePolyZp.Canonical selection.prime.toNat factor)
    (hleadingSemantic : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) =
        (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect
      (SparsePolyZZ.toPoly f) selection)
    (hadjust : StrictHensel.HenselAdjustFirstFactorCorrect f
      selection.factors selection.prime adjusted) :
    ∀ index (hindex : index < adjusted.size),
      Irreducible (SparsePolyZp.toPoly selection.prime.toNat
        adjusted[index]) := by
  have hcandidate : StrictSelectPrime.CandidateCorrect
      (SparsePolyZZ.toPoly f) selection.prime.toNat
      (StrictSelectPrime.factorArrayToL2 selection.prime.toNat
        selection.factors) := hselection
  letI : Fact (Nat.Prime selection.prime.toNat) :=
    ⟨hcandidate.goodPrime.prime⟩
  have hleadingNonzero : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) ≠ 0 := by
    intro leading hleading
    rw [hleadingSemantic leading hleading]
    exact hcandidate.goodPrime.lc_nonzero
  have hrel := hadjust.unitRel hp2 hfactors hleadingNonzero
  intro index hindex
  exact hrel.irreducible (selectionFactors_irreducible hselection) index hindex

/-- The actual first-factor adjustment turns the monic finite-field factors
returned by the selected `__factor_Zp` execution into an *exact* factorization
of the integer source reduced modulo the selected prime.  This recovers the
leading unit from the concrete C++ write instead of weakening the Hensel input
to associatedness. -/
theorem selectionAdjustedFactors_product_eq_source
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {adjusted : Array SparsePolyZp}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hp2 : selection.prime.toNat * selection.prime.toNat ≤ UInt64.size)
    (hfactors : ∀ factor ∈ selection.factors.toList,
      SparsePolyZp.Canonical selection.prime.toNat factor)
    (hleadingSemantic : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) =
        (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect
      (SparsePolyZZ.toPoly f) selection)
    (hadjust : StrictHensel.HenselAdjustFirstFactorCorrect f
      selection.factors selection.prime adjusted) :
    (adjusted.toList.map
        (SparsePolyZp.toPoly selection.prime.toNat)).prod =
      Polynomial.map (Int.castRingHom (ZMod selection.prime.toNat))
        (SparsePolyZZ.toPoly f) := by
  let sourceMod := Polynomial.map
    (Int.castRingHom (ZMod selection.prime.toNat)) (SparsePolyZZ.toPoly f)
  let selected := StrictSelectPrime.factorArrayToL2
    selection.prime.toNat selection.factors
  have hselectedMonic : selected.prod.Monic := by
    apply monic_list_product
    intro factor hfactor
    exact (hselection.quality factor hfactor).2
  have hsourceLeading : sourceMod.leadingCoeff =
      ((SparsePolyZZ.toPoly f).leadingCoeff :
        ZMod selection.prime.toNat) := by
    exact Polynomial.leadingCoeff_map_of_leadingCoeff_ne_zero _
      hselection.goodPrime.lc_nonzero
  have hsourceNonzero : sourceMod ≠ 0 := by
    intro hzero
    have := congrArg Polynomial.leadingCoeff hzero
    rw [hsourceLeading] at this
    exact hselection.goodPrime.lc_nonzero (by simpa using this)
  have hselectedExact :
      sourceMod = Polynomial.C sourceMod.leadingCoeff * selected.prod :=
    polynomial_eq_C_leadingCoeff_mul_of_associated_monic hsourceNonzero
      hselectedMonic hselection.productAssociated
  rcases hadjust.product_eq hp2 hfactors with
    ⟨leading, hleading, hadjustedProduct⟩
  rw [hadjustedProduct]
  change Polynomial.C (leading.2 : ZMod selection.prime.toNat) *
      selected.prod = sourceMod
  rw [hleadingSemantic leading hleading, ← hsourceLeading]
  exact hselectedExact.symm

/-- The exact selected-prime factorization now discharges the formerly
explicit source-product premise of the real Hensel loop.  Consequently the
pre-normalized factors extracted from the generated C++ tree multiply to the
integer source at the exact modulus returned by that loop. -/
theorem selectionHenselFactors_preNormalization_product
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {aTarget : Int32} {output : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hp2 : selection.prime.toNat * selection.prime.toNat ≤ UInt64.size)
    (hfactors : ∀ factor ∈ selection.factors.toList,
      SparsePolyZp.Canonical selection.prime.toNat factor)
    (hleadingSemantic : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) =
        (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect
      (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f
      selection.factors selection.prime aTarget output) :
    ∃ adjusted, ∃ extracted, ∃ outputM : Nat,
      StrictHensel.HenselAdjustFirstFactorCorrect f selection.factors
        selection.prime adjusted ∧
      StrictHensel.HenselNormalizeCorrect extracted (outputM : Int) output.1 ∧
      output.2 = (outputM : Int) ∧
      (extracted.toList.map (StrictHensel.toPolyMod outputM)).prod =
        StrictHensel.toPolyMod outputM f := by
  apply hentry.preNormalizationProduct hcount
  intro adjusted hadjust
  have hadjustSize : adjusted.size = selection.factors.size := by
    cases hadjust
    simp [Array.set!]
  rw [← hadjustSize]
  rw [StrictHensel.henselFactorRangeList_full]
  have hexact := selectionAdjustedFactors_product_eq_source hp2 hfactors
    hleadingSemantic hselection hadjust
  simpa [StrictHensel.toPolyMod] using hexact

/-- Product invariant of the public normalized Hensel output.  The only
remaining discrepancy from the source is the concrete unit computed by the
generated normalization branch; both the modulus and the output array are the
ones returned by the actual C++-shaped entry execution. -/
theorem selectionHenselFactors_normalized_product_eq_unit_mul_source
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {aTarget : Int32} {output : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hp2 : selection.prime.toNat * selection.prime.toNat ≤ UInt64.size)
    (hfactors : ∀ factor ∈ selection.factors.toList,
      SparsePolyZp.Canonical selection.prime.toNat factor)
    (hleadingSemantic : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) =
        (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect
      (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f
      selection.factors selection.prime aTarget output) :
    ∃ outputM : Nat, ∃ scale : ZMod outputM,
      output.2 = (outputM : Int) ∧ IsUnit scale ∧
      (output.1.toList.map (StrictHensel.toPolyMod outputM)).prod =
        Polynomial.C scale * StrictHensel.toPolyMod outputM f := by
  rcases selectionHenselFactors_preNormalization_product hcount hp2 hfactors
      hleadingSemantic hselection hentry with
    ⟨_adjusted, extracted, outputM, _hadjust, hnormalize, houtputM,
      hpreProduct⟩
  rcases hnormalize.product_eq_unit_mul with
    ⟨scale, hscaleUnit, hnormalizedProduct⟩
  refine ⟨outputM, scale, houtputM, hscaleUnit, ?_⟩
  rw [hnormalizedProduct, hpreProduct]

/-- Prime-power specialization of the normalized full-product invariant.  The
exponent comes from the actual well-founded Hensel execution and the scale is
the concrete normalization unit at that same returned modulus. -/
theorem selectionHenselFactors_primePower_product_eq_unit_mul_source
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {aTarget : Int32} {output : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hp2 : selection.prime.toNat * selection.prime.toNat ≤ UInt64.size)
    (hfactors : ∀ factor ∈ selection.factors.toList,
      SparsePolyZp.Canonical selection.prime.toNat factor)
    (hleadingSemantic : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) =
        (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect
      (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f
      selection.factors selection.prime aTarget output) :
    ∃ exponent : Nat, 0 < exponent ∧
      ∃ scale : ZMod (selection.prime.toNat ^ exponent),
        output.2 = ((selection.prime.toNat ^ exponent : Nat) : Int) ∧
        IsUnit scale ∧
        (output.1.toList.map
          (StrictHensel.toPolyMod
            (selection.prime.toNat ^ exponent))).prod =
          Polynomial.C scale *
            StrictHensel.toPolyMod
              (selection.prime.toNat ^ exponent) f := by
  rcases selectionHenselFactors_normalized_product_eq_unit_mul_source
      hcount hp2 hfactors hleadingSemantic hselection hentry with
    ⟨outputM, scale, houtputM, hscaleUnit, hproduct⟩
  rcases hentry.outputModulus_eq_prime_pow with
    ⟨exponent, hexponent, hprimePower⟩
  have hmodulus : outputM = selection.prime.toNat ^ exponent := by
    exact_mod_cast houtputM.symm.trans hprimePower
  subst outputM
  exact ⟨exponent, hexponent, scale, hprimePower, hscaleUnit, hproduct⟩

/-- Reduce the concrete prime-power normalization unit and full-product
identity back to the selected prime.  Both unit witnesses are images of the
same physical normalization scalar. -/
theorem selectionHenselFactors_prime_product_eq_unit_mul_source
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {aTarget : Int32} {output : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hp2 : selection.prime.toNat * selection.prime.toNat ≤ UInt64.size)
    (hfactors : ∀ factor ∈ selection.factors.toList,
      SparsePolyZp.Canonical selection.prime.toNat factor)
    (hleadingSemantic : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) =
        (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect
      (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f
      selection.factors selection.prime aTarget output) :
    ∃ exponent : Nat, 0 < exponent ∧
      ∃ scale : ZMod (selection.prime.toNat ^ exponent),
      ∃ scaleAtPrime : ZMod selection.prime.toNat,
        output.2 = ((selection.prime.toNat ^ exponent : Nat) : Int) ∧
        IsUnit scale ∧ IsUnit scaleAtPrime ∧
        scaleAtPrime = (scale.val : ZMod selection.prime.toNat) ∧
        (output.1.toList.map
          (StrictHensel.toPolyMod
            (selection.prime.toNat ^ exponent))).prod =
          Polynomial.C scale * StrictHensel.toPolyMod
            (selection.prime.toNat ^ exponent) f ∧
        (output.1.toList.map
          (StrictHensel.toPolyMod selection.prime.toNat)).prod =
          Polynomial.C scaleAtPrime *
            StrictHensel.toPolyMod selection.prime.toNat f := by
  rcases selectionHenselFactors_primePower_product_eq_unit_mul_source
      hcount hp2 hfactors hleadingSemantic hselection hentry with
    ⟨exponent, hexponent, scale, houtput, hscaleUnit, hlargeProduct⟩
  let prime := selection.prime.toNat
  have hprime : Nat.Prime prime := Fact.out
  have hpowerNe : prime ^ exponent ≠ 0 :=
    pow_ne_zero exponent hprime.ne_zero
  letI : NeZero (prime ^ exponent) := ⟨hpowerNe⟩
  let scaleAtPrime : ZMod prime := scale.val
  have hdivides : prime ∣ prime ^ exponent :=
    dvd_pow_self prime (Nat.ne_of_gt hexponent)
  have hscaleAtPrimeUnit : IsUnit scaleAtPrime := by
    have hmapped := hscaleUnit.map
      (ZMod.castHom hdivides (ZMod prime))
    rw [ZMod.castHom_apply, ZMod.cast_eq_val] at hmapped
    exact hmapped
  let scaleInt : Int := Int.ofNat scale.val
  let fullInteger : Polynomial Int :=
    (output.1.toList.map SparsePolyZZ.toPoly).prod
  let scaledSource : Polynomial Int :=
    Polynomial.C scaleInt * SparsePolyZZ.toPoly f
  have hscaleLarge : (scaleInt : ZMod (prime ^ exponent)) = scale := by
    simp [scaleInt, ZMod.natCast_zmod_val]
  have hmapFullLarge : Polynomial.map
        (Int.castRingHom (ZMod (prime ^ exponent))) fullInteger =
      (output.1.toList.map
        (StrictHensel.toPolyMod (prime ^ exponent))).prod := by
    simp only [fullInteger, Polynomial.map_list_prod]
    apply congrArg List.prod
    rw [List.map_map]
    apply List.map_congr_left
    intro factor hfactor
    rfl
  have hmapScaledLarge : Polynomial.map
        (Int.castRingHom (ZMod (prime ^ exponent))) scaledSource =
      Polynomial.C scale * StrictHensel.toPolyMod
        (prime ^ exponent) f := by
    simp only [scaledSource, Polynomial.map_mul, Polynomial.map_C]
    change Polynomial.C (scaleInt : ZMod (prime ^ exponent)) *
        StrictHensel.toPolyMod (prime ^ exponent) f = _
    rw [hscaleLarge]
  have hlargeInteger : Polynomial.map
        (Int.castRingHom (ZMod (prime ^ exponent)))
        fullInteger =
      Polynomial.map (Int.castRingHom (ZMod (prime ^ exponent)))
        scaledSource :=
    hmapFullLarge.trans (hlargeProduct.trans hmapScaledLarge.symm)
  have hprimeInteger := StrictRecombine.polynomialMap_eq_of_modulus_dvd
    prime (prime ^ exponent) hdivides
    fullInteger scaledSource hlargeInteger
  have hmapFullPrime : Polynomial.map
        (Int.castRingHom (ZMod prime)) fullInteger =
      (output.1.toList.map (StrictHensel.toPolyMod prime)).prod := by
    simp only [fullInteger, Polynomial.map_list_prod]
    apply congrArg List.prod
    rw [List.map_map]
    apply List.map_congr_left
    intro factor hfactor
    rfl
  have hscalePrime : (scaleInt : ZMod prime) = scaleAtPrime := by
    simp [scaleInt, scaleAtPrime]
  have hmapScaledPrime : Polynomial.map
        (Int.castRingHom (ZMod prime)) scaledSource =
      Polynomial.C scaleAtPrime * StrictHensel.toPolyMod prime f := by
    simp only [scaledSource, Polynomial.map_mul, Polynomial.map_C]
    change Polynomial.C (scaleInt : ZMod prime) *
        StrictHensel.toPolyMod prime f = _
    rw [hscalePrime]
  have hprimeProduct :
      (output.1.toList.map (StrictHensel.toPolyMod prime)).prod =
        Polynomial.C scaleAtPrime * StrictHensel.toPolyMod prime f := by
    rw [← hmapFullPrime, hprimeInteger, hmapScaledPrime]
  exact ⟨exponent, hexponent, scale, scaleAtPrime, houtput, hscaleUnit,
    hscaleAtPrimeUnit, rfl, hlargeProduct, hprimeProduct⟩

/-- The exact selected product and the exact physical complement computed by
`removeCombination` are coprime at the selected prime after the source-leading
scaling used by Zassenhaus.  Their product differs from the squarefree source
only by the two concrete units supplied by prime selection and Hensel
normalization. -/
theorem henselCandidate_physicalComplement_coprime
    (p : Nat) [Fact (Nat.Prime p)]
    (source : Polynomial Int) (active complement : Array SparsePolyZZ)
    (candidate : Array Nat) (scale : ZMod p)
    (hlegal : StrictRecombine.LegalCombination active.size candidate.size
      candidate)
    (hremove : Generated.StrictRecombine.removeCombination candidate active =
      .ok complement)
    (hfull : (active.toList.map (StrictHensel.toPolyMod p)).prod =
      Polynomial.C scale *
        Polynomial.map (Int.castRingHom (ZMod p)) source)
    (hscaleUnit : IsUnit scale)
    (hleadingNonzero : (source.leadingCoeff : ZMod p) ≠ 0)
    (hsquarefree : Squarefree
      (Polynomial.map (Int.castRingHom (ZMod p)) source)) :
    IsCoprime
      (complement.toList.map (StrictHensel.toPolyMod p)).prod
      (Polynomial.C (source.leadingCoeff : ZMod p) *
        ((StrictRecombine.selectSourceIndices active.toList candidate.toList).map
          (StrictHensel.toPolyMod p)).prod) := by
  let selected :=
    ((StrictRecombine.selectSourceIndices active.toList candidate.toList).map
      (StrictHensel.toPolyMod p)).prod
  let remainder :=
    (complement.toList.map (StrictHensel.toPolyMod p)).prod
  let sourceMod := Polynomial.map (Int.castRingHom (ZMod p)) source
  have hpartition := StrictRecombine.removeCombination_toPolyMod_product_partition
    p candidate active complement hlegal hremove
  have hproduct : remainder *
        (Polynomial.C (source.leadingCoeff : ZMod p) * selected) =
      Polynomial.C ((source.leadingCoeff : ZMod p) * scale) * sourceMod := by
    dsimp [selected, remainder, sourceMod]
    calc
      _ = Polynomial.C (source.leadingCoeff : ZMod p) *
          (((StrictRecombine.selectSourceIndices active.toList
            candidate.toList).map (StrictHensel.toPolyMod p)).prod *
            (complement.toList.map (StrictHensel.toPolyMod p)).prod) := by
              ring
      _ = Polynomial.C (source.leadingCoeff : ZMod p) *
          (active.toList.map (StrictHensel.toPolyMod p)).prod := by
            rw [hpartition]
      _ = Polynomial.C (source.leadingCoeff : ZMod p) *
          (Polynomial.C scale *
            Polynomial.map (Int.castRingHom (ZMod p)) source) := by
              rw [hfull]
      _ = _ := by rw [Polynomial.C_mul]; ring
  have hleadingUnit : IsUnit (source.leadingCoeff : ZMod p) :=
    isUnit_iff_ne_zero.mpr hleadingNonzero
  have hscalarUnit : IsUnit
      (Polynomial.C ((source.leadingCoeff : ZMod p) * scale)) :=
    Polynomial.isUnit_C.mpr (hleadingUnit.mul hscaleUnit)
  have hassociated : Associated
      (remainder * (Polynomial.C (source.leadingCoeff : ZMod p) * selected))
      sourceMod := by
    rw [hproduct]
    exact (associated_isUnit_mul_left_iff hscalarUnit).mpr
      (Associated.refl sourceMod)
  have hproductSquarefree : Squarefree
      (remainder *
        (Polynomial.C (source.leadingCoeff : ZMod p) * selected)) :=
    hassociated.squarefree_iff.mpr hsquarefree
  exact (IsRelPrime.of_squarefree_mul hproductSquarefree).isCoprime

private theorem origins_preserve_irreducible
    {p : Nat} {inputs : List SparsePolyZp} {outputs : List SparsePolyZZ}
    (horigins : List.Forall₂
      (fun input output => StrictHensel.toPolyMod p output =
        SparsePolyZp.toPoly p input) inputs outputs)
    (hirreducible : ∀ input ∈ inputs,
      Irreducible (SparsePolyZp.toPoly p input)) :
    ∀ output ∈ outputs, Irreducible (StrictHensel.toPolyMod p output) := by
  induction horigins with
  | nil => simp
  | @cons input output inputs outputs horigin horigins ih =>
      intro candidate hcandidate
      rcases List.mem_cons.mp hcandidate with rfl | htail
      · rw [horigin]
        exact hirreducible input (List.mem_cons_self)
      · exact ih
          (fun factor hfactor =>
            hirreducible factor (List.mem_cons_of_mem input hfactor))
          candidate htail

/-- The factors returned by the actual generated Hensel entry remain
irreducible after reduction modulo the selected prime.  The proof follows the
concrete adjustment, tree leaves, lift extraction, and final normalization in
source order. -/
theorem selectionHenselFactors_mod_irreducible
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {aTarget : Int32} {output : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hp2 : selection.prime.toNat * selection.prime.toNat ≤ UInt64.size)
    (hfactors : ∀ factor ∈ selection.factors.toList,
      SparsePolyZp.Canonical selection.prime.toNat factor)
    (hleadingSemantic : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) =
        (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect
      (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f
      selection.factors selection.prime aTarget output) :
    ∀ index (hindex : index < output.1.size),
      Irreducible (StrictHensel.toPolyMod selection.prime.toNat
        output.1[index]) := by
  rcases hentry.preNormalizationOrigins hcount with
    ⟨adjusted, extracted, outputM, hadjust, hnormalize, houtputM,
      horigins, hnormalizeRel⟩
  have hadjusted := selectionFactors_adjusted_irreducible
    hp2 hfactors hleadingSemantic hselection hadjust
  have hadjustSize : adjusted.size = selection.factors.size := by
    cases hadjust with
    | adjusted leading first adjusted hsource hfirst hadjustedEq =>
        have hzero : 0 < selection.factors.size := by omega
        simp [Array.set!, hzero]
  have hfull : StrictHensel.henselFactorRangeList adjusted
      selection.factors.size 0 = adjusted.toList := by
    rw [← hadjustSize]
    exact StrictHensel.henselFactorRangeList_full adjusted
  rw [hfull] at horigins
  have hadjustedList : ∀ factor ∈ adjusted.toList,
      Irreducible (SparsePolyZp.toPoly selection.prime.toNat factor) := by
    intro factor hfactor
    obtain ⟨index, hindex, hfactorEq⟩ := List.mem_iff_getElem.mp hfactor
    subst factor
    simpa using hadjusted index (by simpa using hindex)
  have hextractedList := origins_preserve_irreducible horigins hadjustedList
  have hextracted : ∀ index (hindex : index < extracted.size),
      Irreducible (StrictHensel.toPolyMod selection.prime.toNat
        extracted[index]) := by
    intro index hindex
    exact hextractedList extracted[index]
      (Array.getElem_mem_toList hindex)
  exact hnormalizeRel.irreducible hextracted

/-- The actual Hensel output array is pointwise associated, in the same
source order, with the factor array returned by the concrete selected-prime
execution.  The association is assembled from the generated first-factor
adjustment, exact lifted-leaf origins, and generated final normalization. -/
theorem selectionHenselFactors_pointwise_associated
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {aTarget : Int32} {output : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hp2 : selection.prime.toNat * selection.prime.toNat ≤ UInt64.size)
    (hfactors : ∀ factor ∈ selection.factors.toList,
      SparsePolyZp.Canonical selection.prime.toNat factor)
    (hleadingSemantic : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) =
        (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect
      (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f
      selection.factors selection.prime aTarget output) :
    selection.factors.size = output.1.size ∧
      ∀ index (hinput : index < selection.factors.size)
        (houtput : index < output.1.size),
        Associated
          (SparsePolyZp.toPoly selection.prime.toNat
            selection.factors[index])
          (StrictHensel.toPolyMod selection.prime.toNat output.1[index]) := by
  have hleadingNonzero : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) ≠ 0 := by
    intro leading hleading
    rw [hleadingSemantic leading hleading]
    exact hselection.goodPrime.lc_nonzero
  rcases hentry.preNormalizationOrigins hcount with
    ⟨adjusted, extracted, outputM, hadjust, hnormalize, houtputM,
      horigins, hnormalizeRel⟩
  have hadjustRel := hadjust.unitRel hp2 hfactors hleadingNonzero
  have hadjustSize : adjusted.size = selection.factors.size := hadjustRel.1.symm
  have hfull : StrictHensel.henselFactorRangeList adjusted
      selection.factors.size 0 = adjusted.toList := by
    rw [← hadjustSize]
    exact StrictHensel.henselFactorRangeList_full adjusted
  rw [hfull] at horigins
  have horiginSize : adjusted.size = extracted.size := by
    simpa using horigins.length_eq
  have hsize : selection.factors.size = output.1.size := by
    rw [hadjustRel.1, horiginSize, hnormalizeRel.1]
  refine ⟨hsize, ?_⟩
  intro index hinput houtput
  have hadjusted : index < adjusted.size := by
    rw [← hadjustRel.1]
    exact hinput
  have hextracted : index < extracted.size := by omega
  rcases hadjustRel.2 index hinput hadjusted with
    ⟨adjustScale, hadjustScale, hadjustEq⟩
  have horigin := horigins.get (by simpa using hadjusted)
    (by simpa using hextracted)
  have horiginEq : StrictHensel.toPolyMod selection.prime.toNat
      extracted[index] = SparsePolyZp.toPoly selection.prime.toNat
        adjusted[index] := by
    simpa [Array.getElem_toList] using horigin
  rcases hnormalizeRel.2 index hextracted houtput with
    ⟨normalizeScale, hnormalizeScale, hnormalizeEq⟩
  have hadjustAssociated : Associated
      (SparsePolyZp.toPoly selection.prime.toNat selection.factors[index])
      (SparsePolyZp.toPoly selection.prime.toNat adjusted[index]) := by
    rw [hadjustEq]
    exact associated_unit_mul_right _ _
      (Polynomial.isUnit_C.mpr hadjustScale)
  have horiginAssociated : Associated
      (SparsePolyZp.toPoly selection.prime.toNat adjusted[index])
      (StrictHensel.toPolyMod selection.prime.toNat extracted[index]) :=
    Associated.of_eq horiginEq.symm
  have hnormalizeAssociated : Associated
      (StrictHensel.toPolyMod selection.prime.toNat extracted[index])
      (StrictHensel.toPolyMod selection.prime.toNat output.1[index]) := by
    rw [hnormalizeEq]
    exact associated_unit_mul_right _ _
      (Polynomial.isUnit_C.mpr hnormalizeScale)
  exact hadjustAssociated.trans
    (horiginAssociated.trans hnormalizeAssociated)

private theorem associated_prod_of_forall₂ {R : Type*} [CommMonoid R]
    {left right : List R} (hrel : List.Forall₂ Associated left right) :
    Associated left.prod right.prod := by
  induction hrel with
  | nil => exact Associated.refl 1
  | cons hhead htail ih =>
      simpa only [List.prod_cons] using hhead.mul_mul ih

/-- Product form of `selectionHenselFactors_pointwise_associated`: the
modular product of the concrete selected-prime array is associated with the
modular product of the concrete normalized Hensel output array. -/
theorem selectionHenselFactors_product_associated
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {aTarget : Int32} {output : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hp2 : selection.prime.toNat * selection.prime.toNat ≤ UInt64.size)
    (hfactors : ∀ factor ∈ selection.factors.toList,
      SparsePolyZp.Canonical selection.prime.toNat factor)
    (hleadingSemantic : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) =
        (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect
      (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f
      selection.factors selection.prime aTarget output) :
    Associated
      (StrictSelectPrime.factorArrayToL2 selection.prime.toNat
        selection.factors).prod
      ((output.1.toList.map
        (StrictHensel.toPolyMod selection.prime.toNat)).prod) := by
  have hpointwise := selectionHenselFactors_pointwise_associated hcount hp2
    hfactors hleadingSemantic hselection hentry
  apply associated_prod_of_forall₂
  rw [List.forall₂_iff_get]
  constructor
  · simpa [StrictSelectPrime.factorArrayToL2] using hpointwise.1
  · intro index hleft hright
    have hinput : index < selection.factors.size := by
      simpa [StrictSelectPrime.factorArrayToL2] using hleft
    have houtput : index < output.1.size := by simpa using hright
    simpa [StrictSelectPrime.factorArrayToL2, Array.getElem_toList] using
      hpointwise.2 index hinput houtput

/-- Every genuine integer divisor reduces, at the actually selected prime,
to the product of an occurrence-sensitive sublist of the actual normalized
Hensel output array, up to a modular unit.  The sublist comes from recursive
irreducible divisibility over that concrete finite list; no semantic factor
oracle chooses it. -/
theorem integer_divisor_mod_associated_hensel_sublist
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {aTarget : Int32} {output : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hp2 : selection.prime.toNat * selection.prime.toNat ≤ UInt64.size)
    (hfactors : ∀ factor ∈ selection.factors.toList,
      SparsePolyZp.Canonical selection.prime.toNat factor)
    (hleadingSemantic : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) =
        (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect
      (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f
      selection.factors selection.prime aTarget output)
    (g : Polynomial Int) (hg : g ∣ SparsePolyZZ.toPoly f) :
    ∃ chosen : List SparsePolyZZ,
      chosen.Sublist output.1.toList ∧
      Associated
        (Polynomial.map
          (Int.castRingHom (ZMod selection.prime.toNat)) g)
        ((chosen.map
          (StrictHensel.toPolyMod selection.prime.toNat)).prod) := by
  let atoms := output.1.toList.map
    (StrictHensel.toPolyMod selection.prime.toNat)
  have hirreducible : ∀ atom ∈ atoms, Irreducible atom := by
    intro atom hatom
    rcases List.mem_map.mp hatom with ⟨lifted, hlifted, rfl⟩
    rcases List.mem_iff_getElem.mp hlifted with ⟨index, hindex, rfl⟩
    have hindexArray : index < output.1.size := by simpa using hindex
    simpa [Array.getElem_toList] using
      selectionHenselFactors_mod_irreducible hcount hp2 hfactors
        hleadingSemantic hselection hentry index hindexArray
  have hmapDvd : Polynomial.map
      (Int.castRingHom (ZMod selection.prime.toNat)) g ∣
      Polynomial.map (Int.castRingHom (ZMod selection.prime.toNat))
        (SparsePolyZZ.toPoly f) := by
    rcases hg with ⟨quotient, hfactor⟩
    refine ⟨Polynomial.map
      (Int.castRingHom (ZMod selection.prime.toNat)) quotient, ?_⟩
    rw [hfactor, Polynomial.map_mul]
  have hselectedDvd : Polynomial.map
      (Int.castRingHom (ZMod selection.prime.toNat)) g ∣
      (StrictSelectPrime.factorArrayToL2 selection.prime.toNat
        selection.factors).prod :=
    dvd_trans hmapDvd hselection.productAssociated.dvd
  have hproductAssociated := selectionHenselFactors_product_associated hcount
    hp2 hfactors hleadingSemantic hselection hentry
  have hatomsDvd : Polynomial.map
      (Int.castRingHom (ZMod selection.prime.toNat)) g ∣ atoms.prod :=
    dvd_trans hselectedDvd hproductAssociated.dvd
  rcases divisor_associated_sublist_product atoms hirreducible _ hatomsDvd with
    ⟨chosenMod, hchosenMod, hassociated⟩
  rcases List.sublist_map_iff.mp hchosenMod with
    ⟨chosen, hchosen, rfl⟩
  exact ⟨chosen, hchosen, hassociated⟩

/-- Array-index form consumed by the generated Zassenhaus scanner.  Every
genuine divisor supplies a legal strictly increasing candidate over the
actual Hensel output array, and the selected occurrence-sensitive product is
associated with that divisor modulo the selected prime. -/
theorem integer_divisor_mod_has_legal_hensel_candidate
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {aTarget : Int32} {output : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hp2 : selection.prime.toNat * selection.prime.toNat ≤ UInt64.size)
    (hfactors : ∀ factor ∈ selection.factors.toList,
      SparsePolyZp.Canonical selection.prime.toNat factor)
    (hleadingSemantic : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) =
        (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect
      (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f
      selection.factors selection.prime aTarget output)
    (g : Polynomial Int) (hg : g ∣ SparsePolyZZ.toPoly f) :
    ∃ indices : Array Nat,
      StrictRecombine.LegalCombination output.1.size indices.size indices ∧
      Associated
        (Polynomial.map
          (Int.castRingHom (ZMod selection.prime.toNat)) g)
        (((StrictRecombine.selectSourceIndices output.1.toList indices.toList).map
          (StrictHensel.toPolyMod selection.prime.toNat)).prod) := by
  rcases integer_divisor_mod_associated_hensel_sublist hcount hp2 hfactors
      hleadingSemantic hselection hentry g hg with
    ⟨chosen, hchosen, hassociated⟩
  rcases StrictRecombine.sublist_exists_legal_combination hchosen with
    ⟨indices, hlegal, hselected⟩
  have hlegal' : StrictRecombine.LegalCombination output.1.size
      indices.size indices := by
    simpa [hlegal.1] using hlegal
  refine ⟨indices, hlegal', ?_⟩
  rw [hselected]
  exact hassociated

/-- Every occurrence-sensitive subproduct selected from the actual normalized
Hensel output is monic.  The proof reads the literal coefficient-one head and
canonical sparse representation carried by `HenselLiftEntryCorrect` for each
selected array cell. -/
theorem henselSelectedProduct_monic
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {factors : Array SparsePolyZp} {p : UInt64}
    [Fact (Nat.Prime p.toNat)]
    {aTarget : Int32} {output : Array SparsePolyZZ × ZZ}
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f factors p
      aTarget output)
    (modulus : Nat) (indices : Array Nat)
    (hbound : ∀ position (hposition : position < indices.size),
      indices[position] < output.1.size) :
    ((StrictRecombine.selectSourceIndices output.1.toList indices.toList).map
      (StrictHensel.toPolyMod modulus)).prod.Monic := by
  apply monic_list_product
  intro candidate hcandidate
  rcases List.mem_map.mp hcandidate with ⟨lifted, hlifted, rfl⟩
  unfold StrictRecombine.selectSourceIndices at hlifted
  rcases List.mem_map.mp hlifted with ⟨index, hindex, rfl⟩
  rcases List.mem_iff_getElem.mp hindex with
    ⟨position, hposition, hindexEq⟩
  have hpositionArray : position < indices.size := by simpa using hposition
  have hactive := hbound position hpositionArray
  have hselected : indices[position] = index := by
    rw [← Array.getElem_toList hpositionArray]
    exact hindexEq
  rw [← hselected, getElem!_pos output.1.toList indices[position]
    (by simpa using hactive), Array.getElem_toList hactive]
  exact hentry.outputToPolyModMonic modulus indices[position] hactive

/-- Integer-polynomial version used by the concrete scalar-pruning loops. -/
theorem henselSelectedIntegerProduct_monic
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {factors : Array SparsePolyZp} {p : UInt64}
    [Fact (Nat.Prime p.toNat)]
    {aTarget : Int32} {output : Array SparsePolyZZ × ZZ}
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f factors p
      aTarget output)
    (indices : Array Nat)
    (hbound : ∀ position (hposition : position < indices.size),
      indices[position] < output.1.size) :
    ((StrictRecombine.selectSourceIndices output.1.toList indices.toList).map
      SparsePolyZZ.toPoly).prod.Monic := by
  apply monic_list_product
  intro candidate hcandidate
  rcases List.mem_map.mp hcandidate with ⟨lifted, hlifted, rfl⟩
  unfold StrictRecombine.selectSourceIndices at hlifted
  rcases List.mem_map.mp hlifted with ⟨index, hindex, rfl⟩
  rcases List.mem_iff_getElem.mp hindex with
    ⟨position, hposition, hindexEq⟩
  have hpositionArray : position < indices.size := by simpa using hposition
  have hactive := hbound position hpositionArray
  have hselected : indices[position] = index := by
    rw [← Array.getElem_toList hpositionArray]
    exact hindexEq
  rw [← hselected, getElem!_pos output.1.toList indices[position]
    (by simpa using hactive), Array.getElem_toList hactive]
  exact hentry.outputToPolyMonic indices[position] hactive

/-- At the selected prime, the exact integer trial polynomial built from a
legal physical Hensel candidate is the quotient-leading-scaled genuine
divisor.  This is the base congruence supplied to `hensel_unique`. -/
theorem henselCandidate_scaled_eq_divisor_mod_prime
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {factors : Array SparsePolyZp} {p : UInt64}
    [Fact (Nat.Prime p.toNat)]
    {aTarget : Int32} {output : Array SparsePolyZZ × ZZ}
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f factors p
      aTarget output)
    (divisor quotient : Polynomial Int)
    (hfactor : SparsePolyZZ.toPoly f = divisor * quotient)
    (hdivisorModNonzero : Polynomial.map
      (Int.castRingHom (ZMod p.toNat)) divisor ≠ 0)
    (hdivisorLeading :
      (Polynomial.map
        (Int.castRingHom (ZMod p.toNat)) divisor).leadingCoeff =
          (divisor.leadingCoeff : ZMod p.toNat))
    (candidate : Array Nat)
    (hlegal : StrictRecombine.LegalCombination output.1.size candidate.size
      candidate)
    (hassociated : Associated
      (Polynomial.map (Int.castRingHom (ZMod p.toNat)) divisor)
      (((StrictRecombine.selectSourceIndices output.1.toList
        candidate.toList).map
          (StrictHensel.toPolyMod p.toNat)).prod)) :
    Polynomial.map (Int.castRingHom (ZMod p.toNat))
        (Polynomial.C (SparsePolyZZ.toPoly f).leadingCoeff *
          ((StrictRecombine.selectSourceIndices output.1.toList
            candidate.toList).map SparsePolyZZ.toPoly).prod) =
      Polynomial.map (Int.castRingHom (ZMod p.toNat))
        (Polynomial.C quotient.leadingCoeff * divisor) := by
  have hbound : ∀ position (hposition : position < candidate.size),
      candidate[position] < output.1.size := by
    intro position hposition
    simpa [getElem!_pos candidate position hposition] using
      hlegal.2.2 position hposition
  have hmonic := henselSelectedProduct_monic hentry p.toNat candidate hbound
  have hscaled := leading_scaled_monic_associated_divisor p.toNat
    (SparsePolyZZ.toPoly f) divisor quotient
    (((StrictRecombine.selectSourceIndices output.1.toList candidate.toList).map
      (StrictHensel.toPolyMod p.toNat)).prod)
    hfactor hdivisorModNonzero hdivisorLeading hmonic hassociated
  simpa [Polynomial.map_mul, Polynomial.map_C,
    Polynomial.map_list_prod, List.map_map, StrictHensel.toPolyMod] using
      hscaled

/-- The two scaled factors compared by Hensel uniqueness have literally equal
integer leading coefficients, and that coefficient is prime to the selected
prime.  The first equality uses physical candidate monicity; the second is
exactly the `GoodPrime` leading-coefficient guard. -/
theorem henselCandidate_scaled_leadingCoeff
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {aTarget : Int32} {output : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hselection : StrictSelectPrime.SelectionCorrect
      (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f
      selection.factors selection.prime aTarget output)
    (divisor quotient : Polynomial Int)
    (hfactor : SparsePolyZZ.toPoly f = divisor * quotient)
    (candidate : Array Nat)
    (hlegal : StrictRecombine.LegalCombination output.1.size candidate.size
      candidate) :
    let selectedInteger :=
      ((StrictRecombine.selectSourceIndices output.1.toList candidate.toList).map
        SparsePolyZZ.toPoly).prod
    let liftedScaled :=
      Polynomial.C (SparsePolyZZ.toPoly f).leadingCoeff * selectedInteger
    let divisorScaled := Polynomial.C quotient.leadingCoeff * divisor
    liftedScaled.leadingCoeff = divisorScaled.leadingCoeff ∧
      ¬((selection.prime.toNat : Int) ∣ liftedScaled.leadingCoeff) := by
  dsimp
  have hbound : ∀ position (hposition : position < candidate.size),
      candidate[position] < output.1.size := by
    intro position hposition
    simpa [getElem!_pos candidate position hposition] using
      hlegal.2.2 position hposition
  have hmonic := henselSelectedIntegerProduct_monic hentry candidate hbound
  have hsourceLeading : (SparsePolyZZ.toPoly f).leadingCoeff =
      divisor.leadingCoeff * quotient.leadingCoeff := by
    rw [hfactor, Polynomial.leadingCoeff_mul]
  have hlifted :
      (Polynomial.C (SparsePolyZZ.toPoly f).leadingCoeff *
        ((StrictRecombine.selectSourceIndices output.1.toList
          candidate.toList).map SparsePolyZZ.toPoly).prod).leadingCoeff =
        (SparsePolyZZ.toPoly f).leadingCoeff := by
    rw [Polynomial.leadingCoeff_mul, Polynomial.leadingCoeff_C,
      hmonic.leadingCoeff, mul_one]
  have hdivisor :
      (Polynomial.C quotient.leadingCoeff * divisor).leadingCoeff =
        (SparsePolyZZ.toPoly f).leadingCoeff := by
    rw [Polynomial.leadingCoeff_mul, Polynomial.leadingCoeff_C,
      hsourceLeading]
    ring
  refine ⟨hlifted.trans hdivisor.symm, ?_⟩
  rw [hlifted]
  intro hdivides
  apply hselection.goodPrime.lc_nonzero
  rw [ZMod.intCast_zmod_eq_zero_iff_dvd]
  exact hdivides

/-- Hensel uniqueness for the concrete occurrence-sensitive candidate.  The
first factorization is built from the literal selected product and the exact
array returned by generated reverse erasure; the second is built from the
genuine integer divisor and quotient.  All normalization scalars are explicit
representatives of the unit produced by the physical Hensel entry. -/
theorem henselCandidate_scaled_eq_divisor_mod_primePower
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {aTarget : Int32} {output : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hp2 : selection.prime.toNat * selection.prime.toNat ≤ UInt64.size)
    (hfactors : ∀ factor ∈ selection.factors.toList,
      SparsePolyZp.Canonical selection.prime.toNat factor)
    (hleadingSemantic : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) =
        (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect
      (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f
      selection.factors selection.prime aTarget output)
    (divisor quotient : Polynomial Int)
    (hfactor : SparsePolyZZ.toPoly f = divisor * quotient)
    (hdivisorModNonzero : Polynomial.map
      (Int.castRingHom (ZMod selection.prime.toNat)) divisor ≠ 0)
    (hdivisorLeading :
      (Polynomial.map
        (Int.castRingHom (ZMod selection.prime.toNat)) divisor).leadingCoeff =
          (divisor.leadingCoeff : ZMod selection.prime.toNat))
    (candidate : Array Nat)
    (hlegal : StrictRecombine.LegalCombination output.1.size candidate.size
      candidate)
    (hassociated : Associated
      (Polynomial.map
        (Int.castRingHom (ZMod selection.prime.toNat)) divisor)
      (((StrictRecombine.selectSourceIndices output.1.toList
        candidate.toList).map
          (StrictHensel.toPolyMod selection.prime.toNat)).prod)) :
    ∃ exponent : Nat, 0 < exponent ∧
      output.2 = ((selection.prime.toNat ^ exponent : Nat) : Int) ∧
      Polynomial.map
        (Int.castRingHom
          (ZMod (selection.prime.toNat ^ exponent)))
        (Polynomial.C (SparsePolyZZ.toPoly f).leadingCoeff *
          ((StrictRecombine.selectSourceIndices output.1.toList
            candidate.toList).map SparsePolyZZ.toPoly).prod) =
      Polynomial.map
        (Int.castRingHom
          (ZMod (selection.prime.toNat ^ exponent)))
        (Polynomial.C quotient.leadingCoeff * divisor) := by
  rcases selectionHenselFactors_prime_product_eq_unit_mul_source
      hcount hp2 hfactors hleadingSemantic hselection hentry with
    ⟨exponent, hexponent, scale, scaleAtPrime, houtput, hscaleUnit,
      hscaleAtPrimeUnit, hscaleAtPrime, hfullLarge, hfullPrime⟩
  rcases StrictRecombine.removeCombination_succeeds candidate output.1 hlegal with
    ⟨complement, hremove⟩
  let prime := selection.prime.toNat
  let modulus := prime ^ exponent
  let source := SparsePolyZZ.toPoly f
  let selectedInteger :=
    ((StrictRecombine.selectSourceIndices output.1.toList candidate.toList).map
      SparsePolyZZ.toPoly).prod
  let complementInteger :=
    (complement.toList.map SparsePolyZZ.toPoly).prod
  let liftedScaled := Polynomial.C source.leadingCoeff * selectedInteger
  let divisorScaled := Polynomial.C quotient.leadingCoeff * divisor
  let scaleInt : Int := Int.ofNat scale.val
  let complementScaled :=
    Polynomial.C (scaleInt * divisor.leadingCoeff) * quotient
  let commonSource :=
    Polynomial.C scaleInt * Polynomial.C source.leadingCoeff * source
  have hprime : Nat.Prime prime := Fact.out
  have hmodulusNe : modulus ≠ 0 :=
    pow_ne_zero exponent hprime.ne_zero
  letI : NeZero modulus := ⟨hmodulusNe⟩
  have hdivides : prime ∣ modulus :=
    dvd_pow_self prime (Nat.ne_of_gt hexponent)
  have hscaleLarge : (scaleInt : ZMod modulus) = scale := by
    simp [scaleInt, modulus, prime]
  have hsourceLeading : source.leadingCoeff =
      divisor.leadingCoeff * quotient.leadingCoeff := by
    dsimp [source]
    rw [hfactor, Polynomial.leadingCoeff_mul]
  have hpartitionLarge :=
    StrictRecombine.removeCombination_toPolyMod_product_partition
      modulus candidate output.1 complement hlegal hremove
  have hselectedMapLarge : Polynomial.map
        (Int.castRingHom (ZMod modulus)) selectedInteger =
      ((StrictRecombine.selectSourceIndices output.1.toList
        candidate.toList).map (StrictHensel.toPolyMod modulus)).prod := by
    dsimp [selectedInteger]
    rw [Polynomial.map_list_prod, List.map_map]
    apply congrArg List.prod
    apply List.map_congr_left
    intro factor hfactorMem
    rfl
  have hcomplementMapLarge : Polynomial.map
        (Int.castRingHom (ZMod modulus)) complementInteger =
      (complement.toList.map (StrictHensel.toPolyMod modulus)).prod := by
    dsimp [complementInteger]
    rw [Polynomial.map_list_prod, List.map_map]
    apply congrArg List.prod
    apply List.map_congr_left
    intro factor hfactorMem
    rfl
  have hsourceMapLarge : Polynomial.map
      (Int.castRingHom (ZMod modulus)) source =
        StrictHensel.toPolyMod modulus f := rfl
  have hproductLifted : Polynomial.map
        (Int.castRingHom (ZMod modulus))
        (complementInteger * liftedScaled) =
      Polynomial.map (Int.castRingHom (ZMod modulus)) commonSource := by
    simp only [Polynomial.map_mul, liftedScaled, commonSource,
      Polynomial.map_C]
    rw [hselectedMapLarge, hcomplementMapLarge]
    change _ * (Polynomial.C (source.leadingCoeff : ZMod modulus) * _) = _
    calc
      _ = Polynomial.C (source.leadingCoeff : ZMod modulus) *
          (((StrictRecombine.selectSourceIndices output.1.toList
            candidate.toList).map (StrictHensel.toPolyMod modulus)).prod *
            (complement.toList.map
              (StrictHensel.toPolyMod modulus)).prod) := by ring
      _ = Polynomial.C (source.leadingCoeff : ZMod modulus) *
          (output.1.toList.map
            (StrictHensel.toPolyMod modulus)).prod := by
              rw [hpartitionLarge]
      _ = Polynomial.C (source.leadingCoeff : ZMod modulus) *
          (Polynomial.C scale * StrictHensel.toPolyMod modulus f) := by
            simpa [modulus, prime] using congrArg
              (fun value => Polynomial.C (source.leadingCoeff : ZMod modulus) * value)
              hfullLarge
      _ = Polynomial.C (scaleInt : ZMod modulus) *
          Polynomial.C (source.leadingCoeff : ZMod modulus) *
          Polynomial.map (Int.castRingHom (ZMod modulus)) source := by
            rw [hsourceMapLarge, hscaleLarge]
            ring
  have hfactorMapLarge := congrArg
    (Polynomial.map (Int.castRingHom (ZMod modulus))) hfactor
  have hproductDivisor : Polynomial.map
        (Int.castRingHom (ZMod modulus))
        (complementScaled * divisorScaled) =
      Polynomial.map (Int.castRingHom (ZMod modulus)) commonSource := by
    simp only [complementScaled, divisorScaled, commonSource,
      Polynomial.map_mul, Polynomial.map_C]
    rw [Polynomial.map_mul] at hfactorMapLarge
    change (Polynomial.C ((scaleInt * divisor.leadingCoeff : Int) : ZMod modulus) *
        Polynomial.map (Int.castRingHom (ZMod modulus)) quotient) *
      (Polynomial.C (quotient.leadingCoeff : ZMod modulus) *
        Polynomial.map (Int.castRingHom (ZMod modulus)) divisor) = _
    rw [Int.cast_mul, Polynomial.C_mul]
    rw [hsourceLeading]
    simp only [map_mul]
    change _ = _ * Polynomial.map (Int.castRingHom (ZMod modulus)) source
    have hfactorMapLarge' :
        Polynomial.map (Int.castRingHom (ZMod modulus)) source =
          Polynomial.map (Int.castRingHom (ZMod modulus)) divisor *
            Polynomial.map (Int.castRingHom (ZMod modulus)) quotient := by
      simpa [source] using hfactorMapLarge
    rw [hfactorMapLarge']
    have hcast (value : Int) :
        (Int.castRingHom (ZMod modulus)) value =
          (value : ZMod modulus) := rfl
    rw [hcast scaleInt, hcast divisor.leadingCoeff,
      hcast quotient.leadingCoeff]
    ring
  have hbaseB := henselCandidate_scaled_eq_divisor_mod_prime hentry
    divisor quotient hfactor hdivisorModNonzero hdivisorLeading candidate
    hlegal hassociated
  have hcop := henselCandidate_physicalComplement_coprime prime source
    output.1 complement candidate scaleAtPrime hlegal hremove
    (by simpa [prime, source] using hfullPrime) hscaleAtPrimeUnit
    (by simpa [prime, source] using hselection.goodPrime.lc_nonzero)
    (by simpa [prime, source] using hselection.goodPrime.sqfree)
  have hleading := henselCandidate_scaled_leadingCoeff hselection hentry
    divisor quotient hfactor candidate hlegal
  dsimp only at hleading
  have hproductLiftedPrime := StrictRecombine.polynomialMap_eq_of_modulus_dvd
    prime modulus hdivides (complementInteger * liftedScaled) commonSource
    hproductLifted
  have hproductDivisorPrime := StrictRecombine.polynomialMap_eq_of_modulus_dvd
    prime modulus hdivides (complementScaled * divisorScaled) commonSource
    hproductDivisor
  have hselectedMapPrime : Polynomial.map
        (Int.castRingHom (ZMod prime)) selectedInteger =
      ((StrictRecombine.selectSourceIndices output.1.toList
        candidate.toList).map (StrictHensel.toPolyMod prime)).prod := by
    dsimp [selectedInteger]
    rw [Polynomial.map_list_prod, List.map_map]
    apply congrArg List.prod
    apply List.map_congr_left
    intro factor hfactorMem
    rfl
  have hcomplementMapPrime : Polynomial.map
        (Int.castRingHom (ZMod prime)) complementInteger =
      (complement.toList.map (StrictHensel.toPolyMod prime)).prod := by
    dsimp [complementInteger]
    rw [Polynomial.map_list_prod, List.map_map]
    apply congrArg List.prod
    apply List.map_congr_left
    intro factor hfactorMem
    rfl
  have hbaseB' : Polynomial.map (Int.castRingHom (ZMod prime)) liftedScaled =
      Polynomial.map (Int.castRingHom (ZMod prime)) divisorScaled := by
    simpa [prime, source, liftedScaled, divisorScaled, selectedInteger] using hbaseB
  have hbaseBNonzero :
      Polynomial.map (Int.castRingHom (ZMod prime)) divisorScaled ≠ 0 := by
    rw [← hbaseB']
    simp only [liftedScaled, Polynomial.map_mul, Polynomial.map_C]
    rw [hselectedMapPrime]
    exact mul_ne_zero
      (Polynomial.C_ne_zero.mpr (by
        simpa [prime, source] using hselection.goodPrime.lc_nonzero))
      (henselSelectedProduct_monic hentry prime candidate (by
        intro position hposition
        simpa [getElem!_pos candidate position hposition] using
          hlegal.2.2 position hposition)).ne_zero
  have hbaseA :
      Polynomial.map (Int.castRingHom (ZMod prime)) complementInteger =
        Polynomial.map (Int.castRingHom (ZMod prime)) complementScaled := by
    have heq :
        Polynomial.map (Int.castRingHom (ZMod prime)) complementInteger *
            Polynomial.map (Int.castRingHom (ZMod prime)) liftedScaled =
          Polynomial.map (Int.castRingHom (ZMod prime)) complementScaled *
            Polynomial.map (Int.castRingHom (ZMod prime)) divisorScaled := by
      calc
        _ = Polynomial.map (Int.castRingHom (ZMod prime))
              (complementInteger * liftedScaled) := by
                exact (Polynomial.map_mul (Int.castRingHom (ZMod prime))).symm
        _ = Polynomial.map (Int.castRingHom (ZMod prime)) commonSource :=
              hproductLiftedPrime
        _ = Polynomial.map (Int.castRingHom (ZMod prime))
              (complementScaled * divisorScaled) := hproductDivisorPrime.symm
        _ = _ := by rw [Polynomial.map_mul]
    rw [hbaseB'] at heq
    exact mul_right_cancel₀ hbaseBNonzero heq
  have hcop' : IsCoprime
      (Polynomial.map (Int.castRingHom (ZMod prime)) complementInteger)
      (Polynomial.map (Int.castRingHom (ZMod prime)) liftedScaled) := by
    simpa [prime, source, liftedScaled, hselectedMapPrime,
      hcomplementMapPrime] using hcop
  have hunique := hensel_unique prime hprime exponent hexponent commonSource
    complementInteger liftedScaled complementScaled divisorScaled
    hproductLifted hproductDivisor hbaseA hbaseB' hcop'
    hleading.1 hleading.2
  exact ⟨exponent, hexponent, houtput, hunique.2⟩

/-- The first scalar-pruning branch of the literal `zassenhausAttempt` cannot
reject a bounded candidate drawn from the actual normalized Hensel output.
Every selected physical head is one, so the generated accumulator is exactly
the source leading coefficient; the generated Mignotte precision then makes
`ZZ.symmetricMod` recover that coefficient literally. -/
theorem zassenhausLeadingPrune_accepts_hensel_candidate
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {factors : Array SparsePolyZp} {p : UInt64}
    [Fact (Nat.Prime p.toNat)]
    {output : Array SparsePolyZZ × ZZ}
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f factors p 0
      output)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hnonempty : 0 < f.size) (hdegree : f[0].1.deg < 2 ^ 63)
    (leading : UMonomial × ZZ) (hleading : f[0]? = some leading)
    (indices : Array Nat)
    (hbound : ∀ position (hposition : position < indices.size),
      indices[position] < output.1.size) :
    ∃ leadingProduct,
      Generated.StrictRecombine.selectedLeadingProductLoop indices output.1
        0 leading.2 = .ok leadingProduct ∧
      ¬(ZZ.symmetricMod leadingProduct output.2 ≠ 0 ∧
        ZZ.fdiv_r 0 (leading.2 * leading.2)
          (ZZ.symmetricMod leadingProduct output.2) ≠ 0) := by
  have hselectedMonic := henselSelectedIntegerProduct_monic hentry indices
    hbound
  have hselectedCanonical : ∀ position
      (hposition : position < indices.size),
      StrictPolynomialMod.SparsePolyZZCanonical
        output.1[indices[position]]! := by
    intro position hposition
    have hactive := hbound position hposition
    simpa [getElem!_pos output.1 indices[position] hactive] using
      hentry.outputCanonical indices[position] hactive
  have hselectedNonempty : ∀ position
      (hposition : position < indices.size),
      output.1[indices[position]]!.isEmpty = false := by
    intro position hposition
    have hactive := hbound position hposition
    rcases hentry.outputOneHead indices[position] hactive with
      ⟨head, tail, hlist⟩
    have hsize : 0 < output.1[indices[position]].size := by
      have hlength := congrArg List.length hlist
      simp at hlength
      omega
    simpa [getElem!_pos output.1 indices[position] hactive,
      Array.isEmpty, Nat.ne_of_gt hsize]
  have hleadingValues :=
    StrictRecombine.selectedLeadingValues_prod_eq_leadingCoeff_of_canonical
      indices output.1 hbound hselectedCanonical hselectedNonempty
  have hvaluesOne :
      (StrictRecombine.selectedLeadingValues indices output.1 0).prod = 1 := by
    rw [hleadingValues, hselectedMonic.leadingCoeff]
  have hloop := StrictRecombine.selectedLeadingProductLoop_succeeds indices
    output.1 0 leading.2 hbound hselectedNonempty
  rw [hvaluesOne, mul_one] at hloop
  refine ⟨leading.2, hloop, ?_⟩
  have hsourceNe : SparsePolyZZ.toPoly f ≠ 0 := by
    intro hzero
    have hhead := StrictRecombine.sparsePolyZZ_leadingCoeff_eq_head f
      hcanonical hnonempty
    rw [hzero] at hhead
    exact hcanonical.2 f[0] (Array.getElem_mem_toList hnonempty)
      (by simpa using hhead.symm)
  have hprecision := StrictRecombine.hensel_output_modulus_bounds_scaled_divisor
    hentry rfl hcanonical hnonempty hdegree leading hleading
      (1 : Polynomial Int) hsourceNe (one_dvd _)
  let modulus := output.2.toNat
  have hmodulusCast : (modulus : Int) = output.2 := by
    exact Int.toNat_of_nonneg hprecision.1.le
  have hmodulus : 0 < modulus := by
    exact Int.pos_iff_toNat_pos.mp hprecision.1
  have hleadingBound : leading.2.natAbs * 2 < modulus := by
    have hzeroBound := hprecision.2 0
    simpa [modulus, hmodulusCast] using hzeroBound
  have hrecovered : ZZ.symmetricMod leading.2 output.2 = leading.2 := by
    rw [← hmodulusCast]
    exact StrictRecombine.symmetricMod_eq_of_strict_bound leading.2 modulus
      hmodulus hleadingBound
  rw [hrecovered]
  exact StrictRecombine.zassenhaus_prune_condition_false_of_dvd
    (leading.2 * leading.2) leading.2 (dvd_mul_right _ _)

/-- The literal constant-coefficient pruning branch accepts the genuine
divisor candidate.  Hensel uniqueness identifies the generated accumulator
with the quotient-leading-scaled divisor at the full returned modulus; the
generated Mignotte precision then recovers its constant coefficient exactly. -/
theorem zassenhausConstantPrune_accepts_hensel_divisor_candidate
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {output : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hp2 : selection.prime.toNat * selection.prime.toNat ≤ UInt64.size)
    (hfactors : ∀ factor ∈ selection.factors.toList,
      SparsePolyZp.Canonical selection.prime.toNat factor)
    (hleadingSemantic : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) =
        (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect
      (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f
      selection.factors selection.prime 0 output)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hnonempty : 0 < f.size) (hdegree : f[0].1.deg < 2 ^ 63)
    (leading : UMonomial × ZZ) (hleading : f[0]? = some leading)
    (divisor quotient : Polynomial Int)
    (hfactor : SparsePolyZZ.toPoly f = divisor * quotient)
    (hdivisorModNonzero : Polynomial.map
      (Int.castRingHom (ZMod selection.prime.toNat)) divisor ≠ 0)
    (hdivisorLeading :
      (Polynomial.map
        (Int.castRingHom (ZMod selection.prime.toNat)) divisor).leadingCoeff =
          (divisor.leadingCoeff : ZMod selection.prime.toNat))
    (candidate : Array Nat)
    (hlegal : StrictRecombine.LegalCombination output.1.size candidate.size
      candidate)
    (hassociated : Associated
      (Polynomial.map
        (Int.castRingHom (ZMod selection.prime.toNat)) divisor)
      (((StrictRecombine.selectSourceIndices output.1.toList
        candidate.toList).map
          (StrictHensel.toPolyMod selection.prime.toNat)).prod)) :
    ∃ constantProduct,
      Generated.StrictRecombine.selectedConstantProductLoop candidate output.1
        0 leading.2 = .ok constantProduct ∧
      ZZ.symmetricMod constantProduct output.2 =
        (Polynomial.C quotient.leadingCoeff * divisor).coeff 0 ∧
      ¬(ZZ.symmetricMod constantProduct output.2 ≠ 0 ∧
        ZZ.fdiv_r 0 (leading.2 *
          Generated.StrictRecombine.constantTerm f)
          (ZZ.symmetricMod constantProduct output.2) ≠ 0) := by
  have hbound : ∀ position (hposition : position < candidate.size),
      candidate[position] < output.1.size := by
    intro position hposition
    simpa [getElem!_pos candidate position hposition] using
      hlegal.2.2 position hposition
  have hselectedCanonical : ∀ position
      (hposition : position < candidate.size),
      StrictPolynomialMod.SparsePolyZZCanonical
        output.1[candidate[position]]! := by
    intro position hposition
    have hactive := hbound position hposition
    simpa [getElem!_pos output.1 candidate[position] hactive] using
      hentry.outputCanonical candidate[position] hactive
  have hconstantValues :=
    StrictRecombine.selectedConstantValues_prod_eq_coeff_zero_of_canonical
      candidate output.1 hbound hselectedCanonical
  let selectedInteger :=
    ((StrictRecombine.selectSourceIndices output.1.toList candidate.toList).map
      SparsePolyZZ.toPoly).prod
  have hloop := StrictRecombine.selectedConstantProductLoop_succeeds candidate
    output.1 0 leading.2 hbound
  rw [hconstantValues] at hloop
  change Generated.StrictRecombine.selectedConstantProductLoop candidate
      output.1 0 leading.2 =
    .ok (leading.2 * selectedInteger.coeff 0) at hloop
  rcases henselCandidate_scaled_eq_divisor_mod_primePower hcount hp2 hfactors
      hleadingSemantic hselection hentry divisor quotient hfactor
      hdivisorModNonzero hdivisorLeading candidate hlegal hassociated with
    ⟨exponent, hexponent, houtput, hunique⟩
  let modulus := selection.prime.toNat ^ exponent
  let target := Polynomial.C quotient.leadingCoeff * divisor
  have hmodulus : 0 < modulus := by
    exact pow_pos (Fact.out : Nat.Prime selection.prime.toNat).pos exponent
  have hsourceLeading : (SparsePolyZZ.toPoly f).leadingCoeff = leading.2 := by
    have hfront : f[0] = leading := by
      rw [Array.getElem?_eq_getElem hnonempty] at hleading
      exact Option.some.inj hleading
    rw [StrictRecombine.sparsePolyZZ_leadingCoeff_eq_head f hcanonical
      hnonempty, hfront]
  have hcoefficientCongruence :
      ((leading.2 * selectedInteger.coeff 0 : Int) : ZMod modulus) =
        (target.coeff 0 : ZMod modulus) := by
    have hcoeff := congrArg
      (fun polynomial : Polynomial (ZMod modulus) => polynomial.coeff 0)
      hunique
    simpa [selectedInteger, target, Polynomial.map_mul, Polynomial.map_C,
      Polynomial.coeff_map, hsourceLeading] using hcoeff
  have hsourceNe : SparsePolyZZ.toPoly f ≠ 0 := by
    intro hzero
    apply hselection.goodPrime.lc_nonzero
    rw [hzero]
    simp
  have hprecision := StrictRecombine.hensel_output_modulus_bounds_scaled_divisor
    hentry rfl hcanonical hnonempty hdegree leading hleading divisor
      hsourceNe (hfactor ▸ dvd_mul_right divisor quotient)
  have hdivisorNe : divisor ≠ 0 := by
    intro hzero
    apply hdivisorModNonzero
    simp [hzero]
  have hdivisorLeadingNe : divisor.leadingCoeff ≠ 0 := by
    exact Polynomial.leadingCoeff_ne_zero.mpr hdivisorNe
  have hsourceLeadingFactor : leading.2 =
      divisor.leadingCoeff * quotient.leadingCoeff := by
    rw [← hsourceLeading, hfactor, Polynomial.leadingCoeff_mul]
  have htargetBound : (target.coeff 0).natAbs * 2 < modulus := by
    have hlargeInt := hprecision.2 0
    rw [houtput] at hlargeInt
    have hlarge : (leading.2 * divisor.coeff 0).natAbs * 2 < modulus := by
      exact_mod_cast hlargeInt
    have hdivisorAbs : 0 < divisor.leadingCoeff.natAbs :=
      Int.natAbs_pos.mpr hdivisorLeadingNe
    have hle : (quotient.leadingCoeff * divisor.coeff 0).natAbs * 2 ≤
        (leading.2 * divisor.coeff 0).natAbs * 2 := by
      rw [hsourceLeadingFactor, Int.natAbs_mul, Int.natAbs_mul,
        Int.natAbs_mul]
      have hbase := Nat.le_mul_of_pos_left
        (quotient.leadingCoeff.natAbs * (divisor.coeff 0).natAbs)
        hdivisorAbs
      have hscaled := Nat.mul_le_mul_right 2 hbase
      simpa [Nat.mul_assoc, Nat.mul_left_comm, Nat.mul_comm] using hscaled
    have htargetCoeff : target.coeff 0 =
        quotient.leadingCoeff * divisor.coeff 0 := by
      simp [target]
    rw [htargetCoeff]
    exact lt_of_le_of_lt hle hlarge
  have hrecovered := StrictRecombine.symmetricMod_eq_of_congruent_strict_bound
    (leading.2 * selectedInteger.coeff 0) (target.coeff 0) modulus hmodulus
    hcoefficientCongruence htargetBound
  have hmodulusInt : (modulus : Int) = output.2 := by
    simpa [modulus] using houtput.symm
  rw [hmodulusInt] at hrecovered
  refine ⟨leading.2 * selectedInteger.coeff 0, hloop, hrecovered, ?_⟩
  rw [hrecovered]
  apply StrictRecombine.zassenhaus_prune_condition_false_of_dvd
  have hsourceConstant : (SparsePolyZZ.toPoly f).coeff 0 =
      divisor.coeff 0 * quotient.coeff 0 := by
    rw [hfactor]
    simp
  have hfConstant : Generated.StrictRecombine.constantTerm f =
      (SparsePolyZZ.toPoly f).coeff 0 :=
    StrictRecombine.sparsePolyZZ_constantTerm_eq_coeff_zero f hcanonical
  rw [hfConstant, hsourceConstant, hsourceLeadingFactor]
  refine ⟨divisor.leadingCoeff * quotient.coeff 0, ?_⟩
  simp [target]
  ring

/-- After both literal scalar prunes accept a genuine Hensel candidate, the
generated checked-index conversion, modular trial product, symmetric recovery,
and primitive-part calls all execute.  The physical symmetric array denotes
exactly the quotient-leading-scaled genuine divisor. -/
theorem zassenhausCandidate_executes_through_primitive
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {output : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hp2 : selection.prime.toNat * selection.prime.toNat ≤ UInt64.size)
    (hfactors : ∀ factor ∈ selection.factors.toList,
      SparsePolyZp.Canonical selection.prime.toNat factor)
    (hleadingSemantic : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) =
        (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect
      (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f
      selection.factors selection.prime 0 output)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hnonempty : 0 < f.size) (hdegree : f[0].1.deg < 2 ^ 63)
    (leading : UMonomial × ZZ) (hleading : f[0]? = some leading)
    (divisor quotient : Polynomial Int)
    (hfactor : SparsePolyZZ.toPoly f = divisor * quotient)
    (hdivisorModNonzero : Polynomial.map
      (Int.castRingHom (ZMod selection.prime.toNat)) divisor ≠ 0)
    (hdivisorLeading :
      (Polynomial.map
        (Int.castRingHom (ZMod selection.prime.toNat)) divisor).leadingCoeff =
          (divisor.leadingCoeff : ZMod selection.prime.toNat))
    (candidate : Array Nat)
    (hlegal : StrictRecombine.LegalCombination output.1.size candidate.size
      candidate)
    (hactiveFits : output.1.size ≤ 2 ^ 31)
    (hassociated : Associated
      (Polynomial.map
        (Int.castRingHom (ZMod selection.prime.toNat)) divisor)
      (((StrictRecombine.selectSourceIndices output.1.toList
        candidate.toList).map
          (StrictHensel.toPolyMod selection.prime.toNat)).prod)) :
    ∃ candidate32 product symmetric content recoveredFactor,
      Generated.StrictRecombine.combinationToInt32 candidate = .ok candidate32 ∧
      Generated.StrictRecombine.trialProductLoop ⟨()⟩ candidate32 output.1
        output.2 0 #[(⟨0⟩, leading.2)] = .ok product ∧
      Generated.StrictRecombine.symmetricModRaw product output.2 = .ok symmetric ∧
      SparsePolyZZ.toPoly symmetric =
        Polynomial.C quotient.leadingCoeff * divisor ∧
      Generated.StrictRecombine.primitiveRaw symmetric =
        .ok (content, recoveredFactor) := by
  have hbound : ∀ position (hposition : position < candidate.size),
      candidate[position] < output.1.size := by
    intro position hposition
    simpa [getElem!_pos candidate position hposition] using
      hlegal.2.2 position hposition
  have hfits : ∀ position (hposition : position < candidate.size),
      candidate[position] < 2 ^ 31 := by
    intro position hposition
    exact lt_of_lt_of_le (hbound position hposition) hactiveFits
  rcases StrictRecombine.combinationToInt32_toList candidate hfits with
    ⟨candidate32, hconvert, _⟩
  have hvalid := StrictRecombine.combinationToInt32_candidate_valid candidate
    output.1.size candidate32 hbound hactiveFits hconvert
  have houtputPositive :=
    (StrictRecombine.hensel_output_modulus_bounds_scaled_divisor hentry rfl
      hcanonical hnonempty hdegree leading hleading divisor
      (by
        intro hzero
        apply hselection.goodPrime.lc_nonzero
        rw [hzero]
        simp)
      (hfactor ▸ dvd_mul_right divisor quotient)).1
  rcases StrictRecombine.trialProductLoop_complete ⟨()⟩ candidate32 output.1
      output.2 0 #[(⟨0⟩, leading.2)] houtputPositive.ne' hvalid with
    ⟨product, hproduct⟩
  rcases henselCandidate_scaled_eq_divisor_mod_primePower hcount hp2 hfactors
      hleadingSemantic hselection hentry divisor quotient hfactor
      hdivisorModNonzero hdivisorLeading candidate hlegal hassociated with
    ⟨exponent, hexponent, houtput, hunique⟩
  let modulus := selection.prime.toNat ^ exponent
  let target := Polynomial.C quotient.leadingCoeff * divisor
  have hmodulus : 0 < modulus :=
    pow_pos (Fact.out : Nat.Prime selection.prime.toNat).pos exponent
  have houtputCast : (modulus : Int) = output.2 := by
    simpa [modulus] using houtput.symm
  have htrial := StrictRecombine.trialProductLoop_source_indices_refines
    modulus hmodulus candidate output.1 candidate32
    #[(⟨0⟩, leading.2)] product hbound hactiveFits hconvert
    (by simpa [houtputCast] using hproduct)
  have hsourceLeading : (SparsePolyZZ.toPoly f).leadingCoeff = leading.2 := by
    have hfront : f[0] = leading := by
      rw [Array.getElem?_eq_getElem hnonempty] at hleading
      exact Option.some.inj hleading
    rw [StrictRecombine.sparsePolyZZ_leadingCoeff_eq_head f hcanonical
      hnonempty, hfront]
  have hcongruent : Polynomial.map (Int.castRingHom (ZMod modulus))
      (SparsePolyZZ.toPoly product) =
      Polynomial.map (Int.castRingHom (ZMod modulus)) target := by
    change StrictHensel.toPolyMod modulus product = _
    rw [htrial]
    have hinitialMod : StrictHensel.toPolyMod modulus
        #[(⟨0⟩, leading.2)] = Polynomial.C (leading.2 : ZMod modulus) := by
      simp [StrictHensel.toPolyMod, SparsePolyZZ.toPoly]
    rw [hinitialMod]
    simpa [modulus, StrictHensel.toPolyMod, target, Polynomial.map_mul,
      Polynomial.map_C, Polynomial.map_list_prod, hsourceLeading] using hunique
  have hinitialCanonical :
      StrictPolynomialMod.SparsePolyZZCanonical #[(⟨0⟩, leading.2)] := by
    constructor
    · simp
    · intro term hterm
      simp at hterm
      subst term
      have hfront : f[0] = leading := by
        rw [Array.getElem?_eq_getElem hnonempty] at hleading
        exact Option.some.inj hleading
      rw [← hfront]
      exact hcanonical.2 f[0] (Array.getElem_mem_toList hnonempty)
  have hproductCanonical := StrictRecombine.trialProductLoop_canonical ⟨()⟩
    candidate32 output.1 output.2 0 #[(⟨0⟩, leading.2)] product
    hinitialCanonical hproduct
  rcases StrictRecombine.symmetricModRaw_complete product output.2
      houtputPositive with ⟨symmetric, hsymmetric⟩
  have hsourceNe : SparsePolyZZ.toPoly f ≠ 0 := by
    intro hzero
    apply hselection.goodPrime.lc_nonzero
    rw [hzero]
    simp
  have hprecision := StrictRecombine.hensel_output_modulus_bounds_scaled_divisor
    hentry rfl hcanonical hnonempty hdegree leading hleading divisor hsourceNe
      (hfactor ▸ dvd_mul_right divisor quotient)
  have hdivisorNe : divisor ≠ 0 := by
    intro hzero
    apply hdivisorModNonzero
    simp [hzero]
  have hsourceLeadingFactor : leading.2 =
      divisor.leadingCoeff * quotient.leadingCoeff := by
    rw [← hsourceLeading, hfactor, Polynomial.leadingCoeff_mul]
  have htargetBound : ∀ degree, (target.coeff degree).natAbs * 2 < modulus := by
    intro degree
    have hlargeInt := hprecision.2 degree
    rw [← houtputCast] at hlargeInt
    have hlarge : (leading.2 * divisor.coeff degree).natAbs * 2 < modulus := by
      exact_mod_cast hlargeInt
    have hdivisorLeadingAbs : 0 < divisor.leadingCoeff.natAbs :=
      Int.natAbs_pos.mpr (Polynomial.leadingCoeff_ne_zero.mpr hdivisorNe)
    have hle : (quotient.leadingCoeff * divisor.coeff degree).natAbs * 2 ≤
        (leading.2 * divisor.coeff degree).natAbs * 2 := by
      rw [hsourceLeadingFactor, Int.natAbs_mul, Int.natAbs_mul,
        Int.natAbs_mul]
      have hbase := Nat.le_mul_of_pos_left
        (quotient.leadingCoeff.natAbs * (divisor.coeff degree).natAbs)
        hdivisorLeadingAbs
      have hscaled := Nat.mul_le_mul_right 2 hbase
      simpa [Nat.mul_assoc, Nat.mul_left_comm, Nat.mul_comm] using hscaled
    have htargetCoeff : target.coeff degree =
        quotient.leadingCoeff * divisor.coeff degree := by simp [target]
    rw [htargetCoeff]
    exact lt_of_le_of_lt hle hlarge
  have hsymmetricPoly :=
    StrictRecombine.symmetricModRaw_recovers_strictly_bounded_target product
      symmetric target modulus hmodulus hproductCanonical
      (by simpa [houtputCast] using hsymmetric) hcongruent htargetBound
  have hsymmetricCanonical := StrictRecombine.symmetricModRaw_canonical product
    symmetric modulus hmodulus hproductCanonical
    (by simpa [houtputCast] using hsymmetric)
  rcases StrictRecombine.primitiveRaw_complete symmetric hsymmetricCanonical with
    ⟨content, recoveredFactor, hprimitive⟩
  exact ⟨candidate32, product, symmetric, content, recoveredFactor, hconvert,
    hproduct, hsymmetric, hsymmetricPoly, hprimitive⟩

/-- If the actual generated fixed-size Zassenhaus scan exhausts, then the
occurrence-sensitive candidate supplied by any genuine integer divisor was
not omitted: that exact index array was executed and rejected.  This is the
bridge from divisor existence to the concrete combination enumerator; the
next completeness step rules out the rejection using symmetric recovery and
exact division. -/
theorem integer_divisor_candidate_rejected_of_scan_exhausted
    {termination : Generated.StrictHensel.DivmodTermination}
    {f fStar : SparsePolyZZ} {selection : PrimeSelectionResult}
    {aTarget : Int32} {output : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hp2 : selection.prime.toNat * selection.prime.toNat ≤ UInt64.size)
    (hfactors : ∀ factor ∈ selection.factors.toList,
      SparsePolyZp.Canonical selection.prime.toNat factor)
    (hleadingSemantic : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) =
        (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect
      (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f
      selection.factors selection.prime aTarget output)
    (g : Polynomial Int) (hg : g ∣ SparsePolyZZ.toPoly f)
    (hrun : StrictRecombine.FixedSizeScanExhausted fStar output.1 output.2
      (integer_divisor_mod_has_legal_hensel_candidate hcount hp2 hfactors
        hleadingSemantic hselection hentry g hg).choose.size) :
    ∃ indices : Array Nat,
      StrictRecombine.LegalCombination output.1.size indices.size indices ∧
      Associated
        (Polynomial.map
          (Int.castRingHom (ZMod selection.prime.toNat)) g)
        (((StrictRecombine.selectSourceIndices output.1.toList indices.toList).map
          (StrictHensel.toPolyMod selection.prime.toNat)).prod) ∧
      Generated.StrictRecombine.zassenhausAttempt fStar output.1 output.2
        indices = .ok .rejected := by
  let witness := integer_divisor_mod_has_legal_hensel_candidate hcount hp2
    hfactors hleadingSemantic hselection hentry g hg
  let indices := witness.choose
  have hspec := witness.choose_spec
  have hrejected := hrun.rejects indices hspec.1
  exact ⟨indices, hspec.1, hspec.2, hrejected⟩

/-- Pairwise coprimality supplied for the concrete adjusted input array is
transported through the actual Hensel leaf origins, lift extraction, and
final source normalization to the returned lifted-factor array. -/
theorem henselFactors_mod_pairwise_coprime
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {factors : Array SparsePolyZp} {p : UInt64}
    {aTarget : Int32} {output : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime p.toNat)]
    (hcount : 2 ≤ factors.size)
    (hpairwise : ∀ adjusted,
      StrictHensel.HenselAdjustFirstFactorCorrect f factors p adjusted →
      ∀ i j (hi : i < adjusted.size) (hj : j < adjusted.size), i < j →
        IsCoprime
          (SparsePolyZp.toPoly p.toNat adjusted[i])
          (SparsePolyZp.toPoly p.toNat adjusted[j]))
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f factors p
      aTarget output) :
    ∀ i j (hi : i < output.1.size) (hj : j < output.1.size), i < j →
      IsCoprime (StrictHensel.toPolyMod p.toNat output.1[i])
        (StrictHensel.toPolyMod p.toNat output.1[j]) := by
  rcases hentry.preNormalizationOrigins hcount with
    ⟨adjusted, extracted, outputM, hadjust, hnormalize, houtputM,
      horigins, hnormalizeRel⟩
  have hadjustSize : adjusted.size = factors.size := by
    cases hadjust with
    | adjusted leading first adjusted hsource hfirst hadjustedEq =>
        have hzero : 0 < factors.size := by omega
        simp [Array.set!, hzero]
  have hfull : StrictHensel.henselFactorRangeList adjusted factors.size 0 =
      adjusted.toList := by
    rw [← hadjustSize]
    exact StrictHensel.henselFactorRangeList_full adjusted
  rw [hfull] at horigins
  have hlength : adjusted.size = extracted.size := by
    simpa using horigins.length_eq
  intro i j hi hj hij
  have hiExtracted : i < extracted.size := by rw [hnormalizeRel.1]; exact hi
  have hjExtracted : j < extracted.size := by rw [hnormalizeRel.1]; exact hj
  have hiAdjusted : i < adjusted.size := by omega
  have hjAdjusted : j < adjusted.size := by omega
  have hiOrigin := horigins.get (by simpa using hiAdjusted)
    (by simpa using hiExtracted)
  have hjOrigin := horigins.get (by simpa using hjAdjusted)
    (by simpa using hjExtracted)
  have hiOrigin' : StrictHensel.toPolyMod p.toNat extracted[i] =
      SparsePolyZp.toPoly p.toNat adjusted[i] := by
    simpa using hiOrigin
  have hjOrigin' : StrictHensel.toPolyMod p.toNat extracted[j] =
      SparsePolyZp.toPoly p.toNat adjusted[j] := by
    simpa using hjOrigin
  have hadjustedCoprime := hpairwise adjusted hadjust i j hiAdjusted
    hjAdjusted hij
  have hextractedCoprime : IsCoprime
      (StrictHensel.toPolyMod p.toNat extracted[i]!)
      (StrictHensel.toPolyMod p.toNat extracted[j]!) := by
    rw [getElem!_pos extracted i hiExtracted,
      getElem!_pos extracted j hjExtracted]
    rw [hiOrigin', hjOrigin']
    simpa [getElem!_pos adjusted i hiAdjusted,
      getElem!_pos adjusted j hjAdjusted] using hadjustedCoprime
  simpa [getElem!_pos output.1 i hi, getElem!_pos output.1 j hj] using
    hnormalizeRel.isCoprime hi hj hextractedCoprime

end StrictFactorZZ

end Refinement
