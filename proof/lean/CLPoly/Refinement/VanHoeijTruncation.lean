import CLPoly.Refinement.VanHoeijLattice

set_option autoImplicit false

namespace Refinement.StrictRecombine

/-- Exact contents of the physical coefficient loop, before normalization.
The statement deliberately does not commute flooring with sums of duplicate
degrees; canonical input is needed for that later semantic bridge. -/
theorem truncateTermsLoop_spec (input output : SparsePolyZZ) (divisor : ZZ) (index : Nat)
    (hrun : Generated.StrictVanHoeij.truncateTermsLoop input divisor index = .ok output) :
    output.size = input.size ∧
    ∀ position, position < input.size → output[position]! =
      if position < index then input[position]!
      else (input[position]!.1, ZZ.fdiv_q 0 input[position]!.2 divisor) := by
  induction hmeasure : input.size - index using Nat.strong_induction_on
      generalizing index input with
  | h measure ih =>
    rw [Generated.StrictVanHoeij.truncateTermsLoop] at hrun
    split at hrun
    next hi =>
      split at hrun
      next hz => contradiction
      next hnz =>
        let updated := input.set index (input[index].1, ZZ.fdiv_q 0 input[index].2 divisor)
        have htail := ih (input.size - (index + 1)) (by omega)
          updated (index + 1) hrun (by simp [updated])
        refine ⟨by simpa [updated] using htail.1, ?_⟩
        intro position hp
        rw [htail.2 position (by simpa [updated] using hp)]
        by_cases heq : position = index
        · subst position
          simp [updated, hi, getElem!_pos]
        · have hentry : updated[position]! = input[position]! := by
            simp [updated, getElem!_pos, hp, Ne.symm heq]
          rw [hentry]
          have hiff : position < index + 1 ↔ position < index := by omega
          simp only [hiff]
    next hi =>
      cases Except.ok.inj hrun
      refine ⟨rfl, ?_⟩
      intro position hp
      simp [show position < index by omega]

theorem truncateTermsLoop_zero_entry (input output : SparsePolyZZ) (divisor : ZZ)
    (hrun : Generated.StrictVanHoeij.truncateTermsLoop input divisor 0 = .ok output)
    (position : Nat) (hp : position < input.size) :
    output[position]! = (input[position]!.1, ZZ.fdiv_q 0 input[position]!.2 divisor) := by
  simpa using (truncateTermsLoop_spec input output divisor 0 hrun).2 position hp

theorem truncateTermsLoop_zero_eq_map (input output : SparsePolyZZ) (divisor : ZZ)
    (hrun : Generated.StrictVanHoeij.truncateTermsLoop input divisor 0 = .ok output) :
    output = input.map (fun term => (term.1, ZZ.fdiv_q 0 term.2 divisor)) := by
  have hsize := (truncateTermsLoop_spec input output divisor 0 hrun).1
  apply Array.ext (by simpa using hsize)
  intro position hp hq
  have hi : position < input.size := by omega
  have he := truncateTermsLoop_zero_entry input output divisor hrun position hi
  simpa only [getElem!_pos output position hp, getElem!_pos input position hi,
    Array.getElem_map] using he

/-- Positive-divisor floor in the C++/Lean model is the same integer
division used in the already proved truncation-error theorem. -/
theorem truncateTermsLoop_zero_floor (input output : SparsePolyZZ) (divisor : ZZ)
    (hd : 0 < divisor)
    (hrun : Generated.StrictVanHoeij.truncateTermsLoop input divisor 0 = .ok output) :
    output = input.map (fun term => (term.1, term.2 / divisor)) := by
  rw [truncateTermsLoop_zero_eq_map input output divisor hrun]
  simp only [ZZ.fdiv_q, Int.fdiv_eq_ediv_of_nonneg _ (le_of_lt hd)]

theorem truncateCldLoop_spec (input output : Array SparsePolyZZ) (divisor : ZZ) (index : Nat)
    (hrun : Generated.StrictVanHoeij.truncateCldLoop input divisor index = .ok output) :
    output.size = input.size ∧
    ∀ position, position < input.size → output[position]! =
      if position < index then input[position]!
      else SparsePolyZZ.normalization
        (input[position]!.map fun term => (term.1, ZZ.fdiv_q 0 term.2 divisor)) := by
  induction hmeasure : input.size - index using Nat.strong_induction_on
      generalizing input index with
  | h measure ih =>
    rw [Generated.StrictVanHoeij.truncateCldLoop] at hrun
    split at hrun
    next hi =>
      split at hrun
      next fault hterms => contradiction
      next poly hterms =>
        have hpoly := truncateTermsLoop_zero_eq_map input[index] poly divisor hterms
        let updated := input.set index (SparsePolyZZ.normalization poly)
        have htail := ih (input.size - (index + 1)) (by omega)
          updated (index + 1) hrun (by simp [updated])
        refine ⟨by simpa [updated] using htail.1, ?_⟩
        intro position hp
        rw [htail.2 position (by simpa [updated] using hp)]
        by_cases heq : position = index
        · subst position
          simp [updated, hpoly, hi, getElem!_pos]
        · have hentry : updated[position]! = input[position]! := by
            simp [updated, getElem!_pos, hp, Ne.symm heq]
          rw [hentry]
          have hiff : position < index + 1 ↔ position < index := by omega
          simp only [hiff]
    next hi =>
      cases Except.ok.inj hrun
      refine ⟨rfl, ?_⟩
      intro position hp
      simp [show position < index by omega]

theorem truncateCld_spec (input output : Array SparsePolyZZ) (divisor : ZZ)
    (hrun : Generated.StrictVanHoeij.truncateCld input divisor = .ok output) :
    0 < divisor ∧ output.size = input.size ∧
    ∀ position, position < input.size → output[position]! =
      SparsePolyZZ.normalization
        (input[position]!.map fun term => (term.1, term.2 / divisor)) := by
  unfold Generated.StrictVanHoeij.truncateCld at hrun
  split at hrun
  next hd =>
    obtain ⟨hs, he⟩ := truncateCldLoop_spec input output divisor 0 hrun
    refine ⟨hd, hs, ?_⟩
    intro position hp
    simpa [ZZ.fdiv_q, Int.fdiv_eq_ediv_of_nonneg _ (le_of_lt hd)] using he position hp
  next hd => contradiction

theorem truncateCld_toPoly (input output : Array SparsePolyZZ) (divisor : ZZ)
    (hrun : Generated.StrictVanHoeij.truncateCld input divisor = .ok output)
    (position : Nat) (hp : position < input.size) :
    SparsePolyZZ.toPoly output[position]! =
      SparsePolyZZ.toPoly (input[position]!.map fun term => (term.1, term.2 / divisor)) := by
  rw [(truncateCld_spec input output divisor hrun).2.2 position hp, normalization_toPoly]

private theorem descending_head_gt (head : UMonomial × Int)
    (tail : List (UMonomial × Int))
    (hchain : List.IsChain (fun a b : UMonomial × Int => a.1.deg > b.1.deg) (head :: tail)) :
    ∀ term ∈ tail, term.1.deg < head.1.deg := by
  induction tail generalizing head with
  | nil => simp
  | cons next rest ih =>
    rw [List.isChain_cons_cons] at hchain
    intro term hterm
    rcases List.mem_cons.mp hterm with rfl | hrest
    · exact hchain.1
    · exact lt_trans (ih next hchain.2 term hrest) hchain.1

private theorem terms_coeff_zero (terms : List (UMonomial × Int))
    (degree : Nat) (coefficient : Int → Int)
    (habsent : ∀ term ∈ terms, term.1.deg ≠ degree) :
    ((terms.map fun term => Polynomial.monomial term.1.deg (coefficient term.2)).sum).coeff degree = 0 := by
  induction terms with
  | nil => simp
  | cons head tail ih =>
    have hh := habsent head (by simp)
    have ht : ∀ term ∈ tail, term.1.deg ≠ degree :=
      fun term hterm => habsent term (List.mem_cons_of_mem head hterm)
    simp [Polynomial.coeff_monomial, hh, ih ht]

theorem descending_terms_map_coeff (terms : List (UMonomial × Int))
    (hchain : List.IsChain (fun a b : UMonomial × Int => a.1.deg > b.1.deg) terms)
    (coefficient : Int → Int) (hzero : coefficient 0 = 0) (degree : Nat) :
    ((terms.map fun term => Polynomial.monomial term.1.deg (coefficient term.2)).sum).coeff degree =
      coefficient (((terms.map fun term => Polynomial.monomial term.1.deg term.2).sum).coeff degree) := by
  induction terms with
  | nil => simp [hzero]
  | cons head tail ih =>
    have ht := List.IsChain.tail hchain
    by_cases heq : head.1.deg = degree
    · have habsent : ∀ term ∈ tail, term.1.deg ≠ degree := by
        intro term hterm
        have hlt := descending_head_gt head tail hchain term hterm
        omega
      have hbefore := terms_coeff_zero tail degree id habsent
      have hafter := terms_coeff_zero tail degree coefficient habsent
      simp only [id_eq] at hbefore
      simp [Polynomial.coeff_monomial, heq, hbefore, hafter]
    · simpa [Polynomial.coeff_monomial, heq] using ih ht

theorem canonical_floor_map_coeff (input : SparsePolyZZ)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical input)
    (divisor : ZZ) (degree : Nat) :
    (SparsePolyZZ.toPoly (input.map fun term => (term.1, term.2 / divisor))).coeff degree =
      (SparsePolyZZ.toPoly input).coeff degree / divisor := by
  simpa [SparsePolyZZ.toPoly, List.map_map, Function.comp_def] using
    descending_terms_map_coeff input.toList hcanonical.1 (fun value => value / divisor)
      (by simp) degree

/-- Exact coefficient bridge used by the canonical lattice short-vector
theorem. Canonical input prevents duplicate degrees from invalidating floor. -/
theorem truncateCld_sparseCoeff (input output : Array SparsePolyZZ) (divisor : ZZ)
    (hrun : Generated.StrictVanHoeij.truncateCld input divisor = .ok output)
    (position : Nat) (hp : position < input.size)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical input[position]!) (degree : Nat) :
    Generated.StrictRecombine.sparseCoeff output[position]! degree =
      Generated.StrictRecombine.sparseCoeff input[position]! degree / divisor := by
  have hout := (truncateCld_spec input output divisor hrun).2.2 position hp
  have houtCanonical : StrictPolynomialMod.SparsePolyZZCanonical output[position]! := by
    rw [hout]
    exact normalization_canonical _
  rw [sparseCoeff_eq_coeff _ houtCanonical, truncateCld_toPoly input output divisor hrun position hp,
    canonical_floor_map_coeff _ hcanonical, sparseCoeff_eq_coeff _ hcanonical]

end Refinement.StrictRecombine
