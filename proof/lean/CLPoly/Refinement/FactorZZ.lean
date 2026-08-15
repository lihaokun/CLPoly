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
