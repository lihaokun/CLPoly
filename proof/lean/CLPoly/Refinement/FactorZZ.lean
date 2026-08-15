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

end StrictFactorZZ

end Refinement
