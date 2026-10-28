import CLPoly.Refinement.VanHoeijResidualBounds
import CLPoly.Refinement.VanHoeijColumnBounds

set_option autoImplicit false

namespace Refinement.StrictRecombine

/-- The physical alternating low/high coefficient scan has this explicit
inverse on all degrees below its width. -/
def cldSpiralIndex (width degree : Nat) : Nat :=
  if 2 * degree < width then 2 * degree else 2 * (width - 1 - degree) + 1

theorem cldSpiralDegree_lt (width index : Nat) (hi : index < width) :
    cldSpiralDegree width index < width := by
  unfold cldSpiralDegree
  split <;> omega

theorem cldSpiralIndex_lt (width degree : Nat) (hd : degree < width) :
    cldSpiralIndex width degree < width := by
  unfold cldSpiralIndex
  split <;> omega

theorem cldSpiralDegree_index (width degree : Nat) (hd : degree < width) :
    cldSpiralDegree width (cldSpiralIndex width degree) = degree := by
  unfold cldSpiralIndex
  split <;> unfold cldSpiralDegree <;> split <;> omega

/-- A successful canonical rebuild cannot silently omit a requested column
when the requested prefix fits the actual width. -/
theorem buildCldMatrixLoop_added_eq_target
    (cld : Array SparsePolyZZ) (target width : Nat) (modulus : ZZ)
    (added : Nat) (matrix : Generated.StrictVanHoeij.LLLMatrix)
    (output : Generated.StrictVanHoeij.LLLMatrix × Nat)
    (htarget : target ≤ width) (hadded : added ≤ target)
    (hrun : Generated.StrictVanHoeij.buildCldMatrixLoop cld 0 target width
      modulus added matrix = .ok output) : output.2 = target := by
  induction hmeasure : target - added using Nat.strong_induction_on generalizing added matrix with
  | h measure ih =>
    rw [Generated.StrictVanHoeij.buildCldMatrixLoop] at hrun
    split at hrun
    next ha =>
      dsimp only at hrun
      split at hrun
      next hw =>
        split at hrun
        next fault he => contradiction
        next nextMatrix he =>
          exact ih (target - (added + 1)) (by omega) (added + 1) nextMatrix
            (by omega) hrun rfl
      next hw => omega
    next ha =>
      cases Except.ok.inj hrun
      omega

theorem buildCldMatrix_full_width
    (cld : Array SparsePolyZZ) (modulus : ZZ)
    (matrix : Generated.StrictVanHoeij.LLLMatrix)
    (output : Generated.StrictVanHoeij.LLLMatrix × Nat)
    (hrun : Generated.StrictVanHoeij.buildCldMatrix matrix cld 0
      (Generated.StrictRecombine.cldSpiralWidth cld) modulus = .ok output) :
    output.2 = Generated.StrictRecombine.cldSpiralWidth cld := by
  exact buildCldMatrixLoop_added_eq_target cld _ _ modulus 0 matrix output
    le_rfl (Nat.zero_le _) hrun

/-- Bounds for the actual maximum-head-degree scan, including the lower
bound for every unprocessed nonempty physical polynomial. -/
theorem cldSpiralWidthLoop_bounds (cld : Array SparsePolyZZ)
    (index maximum limit : Nat) (hmaximum : maximum ≤ limit)
    (hterms : ∀ i, i < cld.size → index ≤ i →
      (if cld[i]!.isEmpty then 0 else (cld[i]!)[0]!.1.deg + 1) ≤ limit) :
    maximum ≤ Generated.StrictRecombine.cldSpiralWidthLoop cld index maximum ∧
      Generated.StrictRecombine.cldSpiralWidthLoop cld index maximum ≤ limit ∧
      ∀ i, i < cld.size → index ≤ i →
        (if cld[i]!.isEmpty then 0 else (cld[i]!)[0]!.1.deg + 1) ≤
          Generated.StrictRecombine.cldSpiralWidthLoop cld index maximum := by
  induction hmeasure : cld.size - index using Nat.strong_induction_on generalizing index maximum with
  | h measure ih =>
    rw [Generated.StrictRecombine.cldSpiralWidthLoop]
    split
    next hi =>
      let termWidth := if cld[index]!.isEmpty then 0 else (cld[index]!)[0]!.1.deg + 1
      have hw : termWidth ≤ limit := hterms index hi le_rfl
      have hnext := ih (cld.size - (index + 1)) (by omega) (index + 1)
        (max maximum termWidth) (max_le hmaximum hw)
        (fun i hi hafter => hterms i hi (by omega)) rfl
      simp only [← getElem!_pos cld index hi]
      change maximum ≤ Generated.StrictRecombine.cldSpiralWidthLoop cld (index + 1) (max maximum termWidth) ∧ _
      refine ⟨(le_max_left _ _).trans hnext.1, hnext.2.1, ?_⟩
      intro i hib hafter
      by_cases heq : i = index
      · subst i
        exact (le_max_right _ _).trans hnext.1
      · exact hnext.2.2 i hib (by omega)
    next hi =>
      exact ⟨le_rfl, hmaximum, by intro i hib hafter; omega⟩

theorem cldSpiralWidth_bounds (cld : Array SparsePolyZZ) (limit : Nat)
    (hterms : ∀ i, i < cld.size →
      (if cld[i]!.isEmpty then 0 else (cld[i]!)[0]!.1.deg + 1) ≤ limit) :
    Generated.StrictRecombine.cldSpiralWidth cld ≤ limit ∧
      ∀ i, i < cld.size →
        (if cld[i]!.isEmpty then 0 else (cld[i]!)[0]!.1.deg + 1) ≤
          Generated.StrictRecombine.cldSpiralWidth cld := by
  have hb := cldSpiralWidthLoop_bounds cld 0 0 limit (Nat.zero_le _)
    (fun i hi _ => hterms i hi)
  exact ⟨hb.2.1, fun i hi => hb.2.2 i hi (Nat.zero_le _)⟩

theorem cldSpiralWidth_covers (cld : Array SparsePolyZZ) (i : Nat) (hi : i < cld.size) :
    (if cld[i]!.isEmpty then 0 else (cld[i]!)[0]!.1.deg + 1) ≤
      Generated.StrictRecombine.cldSpiralWidth cld := by
  let termWidth : Fin cld.size → Nat := fun j =>
    if cld[j.val]!.isEmpty then 0 else (cld[j.val]!)[0]!.1.deg + 1
  have hterms : ∀ j, j < cld.size →
      (if cld[j]!.isEmpty then 0 else (cld[j]!)[0]!.1.deg + 1) ≤ ∑ k, termWidth k := by
    intro j hj
    exact Finset.single_le_sum (fun k _ => Nat.zero_le (termWidth k))
      (Finset.mem_univ (⟨j, hj⟩ : Fin cld.size))
  exact (cldSpiralWidth_bounds cld _ hterms).2 i hi

/-- Every coefficient omitted above the actual scanned width is zero,
including the empty-polynomial case. -/
theorem cldSpiralWidth_sparseCoeff_zero (cld : Array SparsePolyZZ)
    (i : Nat) (hi : i < cld.size)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical cld[i]!)
    (degree : Nat) (hdegree : Generated.StrictRecombine.cldSpiralWidth cld ≤ degree) :
    Generated.StrictRecombine.sparseCoeff cld[i]! degree = 0 := by
  rw [sparseCoeff_eq_coeff _ hcanonical]
  by_cases he : cld[i]!.isEmpty = true
  · have hz := Array.isEmpty_iff.mp he
    simp [hz, SparsePolyZZ.toPoly]
  · have hn : 0 < cld[i]!.size := by
      have : cld[i]!.size ≠ 0 := by simpa [Array.isEmpty] using he
      omega
    have hb := cldSpiralWidth_covers cld i hi
    simp [he] at hb
    have hdeg := sparsePolyZZ_natDegree_eq_head cld[i]! hcanonical hn
    rw [← getElem!_pos cld[i]! 0 hn] at hdeg
    exact Polynomial.coeff_eq_zero_of_natDegree_lt (by omega)

/-- An omitted floor coefficient is an ordinary nonnegative remainder.
Negative input coefficients therefore cannot be silently rounded to zero. -/
theorem floor_zero_abs_lt (value divisor : Int) (hd : 0 < divisor)
    (hz : value / divisor = 0) : |value| < divisor := by
  have he := Int.mul_ediv_add_emod value divisor
  rw [hz, mul_zero, zero_add] at he
  rw [← he, abs_of_nonneg (Int.emod_nonneg value (ne_of_gt hd))]
  exact Int.emod_lt_of_pos value hd

/-- The actual floor-and-normalize path introduces no new degrees, so the
actual width scan respects any strict degree bound on the original CLDs. -/
theorem truncateCld_width_le (original truncated : Array SparsePolyZZ) (divisor : ZZ)
    (hrun : Generated.StrictVanHoeij.truncateCld original divisor = .ok truncated)
    (hcanonical : ∀ i, i < original.size → StrictPolynomialMod.SparsePolyZZCanonical original[i]!)
    (limit : Nat) (hdegree : ∀ i, i < original.size →
      (SparsePolyZZ.toPoly original[i]!).natDegree < limit) :
    Generated.StrictRecombine.cldSpiralWidth truncated ≤ limit := by
  have ht := truncateCld_spec original truncated divisor hrun
  apply (cldSpiralWidth_bounds truncated limit ?_).1
  intro i hi
  have hio : i < original.size := by rw [← ht.2.1]; exact hi
  have htc : StrictPolynomialMod.SparsePolyZZCanonical truncated[i]! := by
    rw [ht.2.2 i hio]
    exact normalization_canonical _
  have hdeg : (SparsePolyZZ.toPoly truncated[i]!).natDegree ≤
      (SparsePolyZZ.toPoly original[i]!).natDegree := by
    apply Polynomial.natDegree_le_iff_coeff_eq_zero.mpr
    intro d hd
    rw [← sparseCoeff_eq_coeff _ htc,
      truncateCld_sparseCoeff original truncated divisor hrun i hio (hcanonical i hio),
      sparseCoeff_eq_coeff _ (hcanonical i hio), Polynomial.coeff_eq_zero_of_natDegree_lt hd]
    simp
  split
  next he => exact Nat.zero_le _
  next he =>
    have hn : 0 < truncated[i]!.size := by
      have : truncated[i]!.size ≠ 0 := by simpa [Array.isEmpty] using he
      omega
    have hh := sparsePolyZZ_natDegree_eq_head truncated[i]! htc hn
    rw [← getElem!_pos truncated[i]! 0 hn] at hh
    have hb := hdegree i hio
    omega

end Refinement.StrictRecombine
