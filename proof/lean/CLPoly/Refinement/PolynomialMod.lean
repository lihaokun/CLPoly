/- Genuine refinement of the generated C++ `polynomial_mod` range-for. -/
import CLPoly.Generated.StrictPolynomialMod
import CLPoly.Math.Univariate
import Mathlib.Algebra.Polynomial.Eval.Defs

set_option autoImplicit false

open Polynomial
open CLPoly.Math

namespace Refinement.StrictPolynomialMod

/-- Representation invariant of the C++ sparse integer polynomial input. -/
def SparsePolyZZCanonical (f : SparsePolyZZ) : Prop :=
  List.IsChain (fun a b : UMonomial × Int => a.1.deg > b.1.deg) f.toList ∧
  ∀ term ∈ f.toList, term.2 ≠ 0

private def reduceTerm (p : UInt64) (term : UMonomial × Int) :
    Option (UMonomial × Zp) :=
  let coefficient := Zp.ofInt term.2 p
  if coefficient.val = 0 then none else some (term.1, coefficient)

theorem polynomialModLoop_eq_filterMap (f : SparsePolyZZ) (p : UInt64)
    (index : Nat) (result : SparsePolyZp) :
    Generated.StrictPolynomialMod.polynomialModLoop f p index result =
      result ++ ((f.toList.drop index).filterMap (reduceTerm p)).toArray := by
  induction hmeasure : f.size - index using Nat.strong_induction_on
      generalizing index result with
  | h measure ih =>
      rw [Generated.StrictPolynomialMod.polynomialModLoop]
      split
      next hindex =>
        let term := f[index]
        let coefficient := Zp.ofInt term.2 p
        let result' := if coefficient.val = 0 then result
          else result.push (term.1, coefficient)
        rw [ih (f.size - (index + 1)) (by omega) (index + 1) result' rfl]
        have hsuffix : f.toList.drop index =
            term :: f.toList.drop (index + 1) := by
          simpa [term] using List.drop_eq_getElem_cons
            (l := f.toList) (i := index) (by simpa using hindex)
        rw [hsuffix, List.filterMap_cons]
        by_cases hz : coefficient.val = 0
        · simp [reduceTerm, coefficient, result', hz]
        · simp [reduceTerm, coefficient, result', hz]
      next hindex =>
        have hle : f.size ≤ index := Nat.le_of_not_gt hindex
        simp [List.drop_eq_nil_iff.mpr hle]

theorem polynomialModLoop_eq_model (f : SparsePolyZZ) (p : UInt64) :
    Generated.StrictPolynomialMod.polynomialModLoop f p 0 #[] =
      polynomial_mod f p := by
  rw [polynomialModLoop_eq_filterMap]
  unfold polynomial_mod reduceTerm
  apply Array.toList_inj.mp
  simp [Array.toList_filterMap]

private theorem zmod_intCast_ofInt (coefficient : Int) (p : UInt64)
    (hp : 0 < p.toNat) :
    Zp.toZMod p.toNat (Zp.ofInt coefficient p) =
      (coefficient : ZMod p.toNat) := by
  unfold Zp.toZMod Zp.ofInt
  have hpInt : (0 : Int) < p.toNat := by exact_mod_cast hp
  have hnonneg : 0 ≤ coefficient.emod p.toNat := Int.emod_nonneg _ (by omega)
  have hltInt : coefficient.emod p.toNat < p.toNat := Int.emod_lt_of_pos _ hpInt
  have hltNatP : (coefficient.emod p.toNat).toNat < p.toNat := by
    exact (Int.toNat_lt hnonneg).2 hltInt
  have hltNat : (coefficient.emod p.toNat).toNat < UInt64.size :=
    lt_trans hltNatP (UInt64.toNat_lt_size p)
  have hremNonneg : ¬coefficient % (p.toNat : Int) < 0 :=
    not_lt_of_ge hnonneg
  simp only [hremNonneg, ↓reduceIte]
  change (((coefficient.emod p.toNat).toNat.toUInt64).toNat :
      ZMod p.toNat) = (coefficient : ZMod p.toNat)
  have hword : ((coefficient.emod p.toNat).toNat.toUInt64).toNat =
      (coefficient.emod p.toNat).toNat := by
    exact UInt64.toNat_ofNat_of_lt hltNat
  rw [hword]
  have htoInt : ((coefficient.emod p.toNat).toNat : Int) =
      coefficient.emod p.toNat := Int.toNat_of_nonneg hnonneg
  calc
    ((coefficient.emod p.toNat).toNat : ZMod p.toNat) =
        (((coefficient.emod p.toNat).toNat : Int) : ZMod p.toNat) := by
          exact_mod_cast rfl
    _ = (coefficient.emod p.toNat : Int) := by rw [htoInt]
    _ = (coefficient : ZMod p.toNat) := ZMod.intCast_mod coefficient p.toNat

theorem polynomial_mod_toPoly (f : SparsePolyZZ) (p : UInt64)
    (hp : 0 < p.toNat) :
    SparsePolyZp.toPoly p.toNat (polynomial_mod f p) =
      Polynomial.map (Int.castRingHom (ZMod p.toNat))
        (SparsePolyZZ.toPoly f) := by
  unfold polynomial_mod SparsePolyZZ.toPoly SparsePolyZp.toPoly
  rw [Array.toList_filterMap]
  induction f.toList with
  | nil => simp [listSum]
  | cons term rest ih =>
      rcases term with ⟨monomial, coefficient⟩
      simp only [List.filterMap_cons, List.map_cons, List.sum_cons,
        Polynomial.map_add, Polynomial.map_monomial, listSum_cons]
      let zp := Zp.ofInt coefficient p
      by_cases hz : zp.val = 0
      · have hcastZero : (coefficient : ZMod p.toNat) = 0 := by
          rw [← zmod_intCast_ofInt coefficient p hp]
          simp [zp, hz, Zp.toZMod]
        simp [zp, hz, hcastZero, ih]
      · simp [zp, hz, zmod_intCast_ofInt coefficient p hp, ih]

theorem polynomial_mod_canonical (f : SparsePolyZZ) (p : UInt64)
    (hp : 0 < p.toNat) (hcanonical : SparsePolyZZCanonical f) :
    SparsePolyZp.Canonical p.toNat (polynomial_mod f p) := by
  unfold polynomial_mod SparsePolyZp.Canonical SparsePolyZp.WellFormed_arr
  rw [Array.toList_filterMap]
  refine ⟨?_, ?_, ?_⟩
  · intro output houtput
    rw [List.mem_filterMap] at houtput
    rcases houtput with ⟨term, hterm, hreduce⟩
    dsimp only at hreduce
    let coefficient := Zp.ofInt term.2 p
    split at hreduce
    next hz => contradiction
    next hz =>
      simp only [Option.some.injEq] at hreduce
      subst output
      refine ⟨rfl, ?_⟩
      simp only [Prod.snd]
      unfold Zp.ofInt
      have hpInt : (0 : Int) < p.toNat := by exact_mod_cast hp
      have hnonneg : 0 ≤ term.2.emod p.toNat :=
        Int.emod_nonneg _ (by omega)
      have hlt : term.2.emod p.toNat < p.toNat :=
        Int.emod_lt_of_pos _ hpInt
      have hremNonneg : ¬term.2 % (p.toNat : Int) < 0 :=
        not_lt_of_ge hnonneg
      simp only [hremNonneg, ↓reduceIte]
      have hltNatP : (term.2.emod p.toNat).toNat < p.toNat :=
        (Int.toNat_lt hnonneg).2 hlt
      have hltWord : (term.2.emod p.toNat).toNat < UInt64.size :=
        lt_trans hltNatP (UInt64.toNat_lt_size p)
      have hword : ((term.2.emod p.toNat).toNat.toUInt64).toNat =
          (term.2.emod p.toNat).toNat :=
        UInt64.toNat_ofNat_of_lt hltWord
      change ((term.2 % (p.toNat : Int)).toNat.toUInt64).toNat < p.toNat
      change ((term.2.emod p.toNat).toNat.toUInt64).toNat < p.toNat
      rw [hword]
      exact hltNatP
  · apply List.isChain_iff_pairwise.mpr
    let degreeGreater : (UMonomial × Int) → (UMonomial × Int) → Prop :=
      fun a b => a.1.deg > b.1.deg
    letI : Trans degreeGreater degreeGreater degreeGreater :=
      ⟨by
        intro a b c hab hbc
        dsimp [degreeGreater] at hab hbc ⊢
        omega⟩
    apply List.Pairwise.filterMap (R := fun a b : UMonomial × Int =>
        a.1.deg > b.1.deg)
      (S := fun a b : UMonomial × Zp => a.1.deg > b.1.deg)
      (fun term =>
        let coefficient := Zp.ofInt term.2 p
        if coefficient.val = 0 then none else some (term.1, coefficient))
      (fun a b hab outputA ha outputB hb => by
        dsimp only at ha hb
        split at ha <;> try contradiction
        split at hb <;> try contradiction
        simp only [Option.some.injEq] at ha hb
        subst outputA
        subst outputB
        exact hab)
    exact List.isChain_iff_pairwise.mp hcanonical.1
  · intro output houtput
    rw [List.mem_filterMap] at houtput
    rcases houtput with ⟨term, _, hreduce⟩
    dsimp only at hreduce
    split at hreduce
    next hz => contradiction
    next hz =>
      simp only [Option.some.injEq] at hreduce
      subst output
      exact hz

theorem polynomial_mod_raw_ir_refines (f : SparsePolyZZ) (p : UInt64)
    (hp : 0 < p.toNat) :
    ∃ result,
      Generated.StrictPolynomialMod.polynomial_mod_raw_ir f p = .ok result ∧
      SparsePolyZp.toPoly p.toNat result =
      Polynomial.map (Int.castRingHom (ZMod p.toNat))
          (SparsePolyZZ.toPoly f) := by
  refine ⟨polynomial_mod f p, ?_, polynomial_mod_toPoly f p hp⟩
  simp [Generated.StrictPolynomialMod.polynomial_mod_raw_ir,
    polynomialModLoop_eq_model]

end Refinement.StrictPolynomialMod
