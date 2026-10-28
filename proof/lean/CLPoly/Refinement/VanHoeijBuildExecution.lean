import CLPoly.Refinement.VanHoeijLattice
import CLPoly.Refinement.VanHoeijInitialExecution

set_option autoImplicit false

namespace Refinement.StrictRecombine
open Generated.StrictRecombine Generated.StrictVanHoeij

theorem fillCldColumnLoop_succeeds (cld : Array SparsePolyZZ) (degree column index : Nat)
    (matrix : Generated.StrictRecombine.LLLMatrix) (hsize : cld.size ≤ matrix.size)
    (hcolumns : ∀ row (hr : row < matrix.size), column < matrix[row].size) :
    ∃ output, fillCldColumnLoop cld degree column index matrix = .ok output := by
  induction hmeasure : cld.size - index using Nat.strong_induction_on generalizing index matrix with
  | h measure ih =>
    rw [fillCldColumnLoop]
    by_cases hi : index < cld.size
    · have hm : index < matrix.size := lt_of_lt_of_le hi hsize
      have hc := hcolumns index hm
      rw [dif_pos hi, dif_pos hm, dif_pos hc]
      apply ih _ (by omega) (index + 1) _ (by simpa using hsize) ?_ rfl
      intro row hr
      have hr' : row < matrix.size := by simpa using hr
      rw [Array.getElem_set]
      split
      · simpa only [Array.size_set] using hc
      · exact hcolumns row hr'
    · exact ⟨matrix, by rw [dif_neg hi]⟩

theorem appendCldColumn_succeeds (matrix : Generated.StrictRecombine.LLLMatrix)
    (cld : Array SparsePolyZZ) (columns degree : Nat) (modulus : ZZ)
    (hsize : matrix.size = cld.size + columns)
    (hrows : ∀ row, row < matrix.size → matrix[row]!.size = matrix.size) :
    ∃ output, Generated.StrictVanHoeij.appendCldColumn matrix cld columns degree modulus = .ok output := by
  obtain ⟨extended, hext⟩ : ∃ output, appendZeroColumn matrix = .ok output :=
    ⟨_, appendZeroColumn_eq matrix⟩
  have heq := Except.ok.inj (hext.symm.trans (appendZeroColumn_eq matrix))
  have hes : extended.size = matrix.size := by rw [heq]; exact appendZeroSuffix_size matrix
  have hec : ∀ row (hr : row < extended.size), cld.size + columns < extended[row].size := by
    intro row hr
    have hr' : row < matrix.size := by omega
    have he := appendZeroSuffix_row matrix row hr'
    rw [← heq] at he
    rw [getElem!_pos extended row hr] at he
    rw [he, Array.size_push, hrows row hr', hsize]
    omega
  obtain ⟨filled, hf⟩ := fillCldColumnLoop_succeeds cld degree (cld.size + columns) 0 extended
    (by omega) hec
  have hc : cld.size + columns < (zeroMatrixRow (cld.size + columns + 1)).size := by
    rw [zeroMatrixRow_size]; omega
  refine ⟨filled.push ((zeroMatrixRow (cld.size + columns + 1)).set (cld.size + columns) modulus), ?_⟩
  unfold Generated.StrictVanHoeij.appendCldColumn
  rw [hext]
  simp only [bind, Except.bind]
  rw [hf]
  simp only [bind, Except.bind]
  rw [dif_pos hc]

theorem buildCldMatrixLoop_succeeds (cld : Array SparsePolyZZ)
    (width current target added : Nat) (scale modulus : ZZ)
    (matrix : Generated.StrictRecombine.LLLMatrix)
    (hcanonical : CanonicalCldBasis cld width scale modulus (current + added) matrix) :
    ∃ output, Generated.StrictVanHoeij.buildCldMatrixLoop cld current target width modulus
      added matrix = .ok output := by
  induction hmeasure : target - added using Nat.strong_induction_on generalizing added matrix with
  | h measure ih =>
    rw [Generated.StrictVanHoeij.buildCldMatrixLoop]
    by_cases ht : added < target
    · rw [if_pos ht]
      by_cases hw : current + added < width
      · rw [if_pos hw]
        obtain ⟨next, hn⟩ := appendCldColumn_succeeds matrix cld (current + added)
          (cldSpiralDegree width (current + added)) modulus hcanonical.size_eq hcanonical.row_size
        have hnext := appendCldColumn_canonical cld width (current + added) scale modulus
          matrix next hcanonical hn
        obtain ⟨output, ho⟩ := ih _ (by omega) (added + 1) next
          (by simpa only [Nat.add_assoc] using hnext) rfl
        refine ⟨output, ?_⟩
        change (match Generated.StrictVanHoeij.appendCldColumn matrix cld (current + added)
          (cldSpiralDegree width (current + added)) modulus with
          | .error fault => Except.error fault
          | .ok next => Generated.StrictVanHoeij.buildCldMatrixLoop cld current target width modulus
            (added + 1) next) = _
        rw [hn]
        exact ho
      · exact ⟨(matrix, added), by rw [if_neg hw]⟩
    · exact ⟨(matrix, added), by rw [if_neg ht]⟩

theorem buildCldMatrix_succeeds (cld : Array SparsePolyZZ) (current target : Nat)
    (scale modulus : ZZ) (matrix : Generated.StrictRecombine.LLLMatrix)
    (hcanonical : CanonicalCldBasis cld (cldSpiralWidth cld) scale modulus current matrix) :
    ∃ output, Generated.StrictVanHoeij.buildCldMatrix matrix cld current target modulus = .ok output := by
  exact buildCldMatrixLoop_succeeds cld (cldSpiralWidth cld) current target 0 scale modulus matrix
    (by simpa using hcanonical)

end Refinement.StrictRecombine
