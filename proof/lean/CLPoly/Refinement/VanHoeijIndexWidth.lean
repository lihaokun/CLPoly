import CLPoly.Refinement.VanHoeijLattice
import CLPoly.Refinement.VanHoeijColumnBounds

set_option autoImplicit false

namespace Refinement.StrictRecombine

/-- Matrix coordinates can exceed signed 32-bit even though both source
factor counts and CLD column counts fit. They fit the 64-bit size_t target. -/
theorem lattice_dimension_fits_word (count columns : Nat)
    (hcount : count < 2 ^ 31) (hcolumns : columns < 2 ^ 31) :
    count + columns < 2 ^ 32 ∧ count + columns < 2 ^ 64 := by omega

theorem lattice_dimension_word_add (count columns : Nat)
    (hcount : count < 2 ^ 31) (hcolumns : columns < 2 ^ 31) :
    (UInt64.ofNat count + UInt64.ofNat columns).toNat = count + columns := by
  rw [← UInt64.ofNat_add]
  exact UInt64.toNat_ofNat_of_lt (lattice_dimension_fits_word count columns hcount hcolumns).2

theorem lattice_dimension_usize_add (count columns : Nat)
    (hcount : count < 2 ^ 31) (hcolumns : columns < 2 ^ 31) :
    (USize.ofNat count + USize.ofNat columns).toNat = count + columns := by
  rw [← USize.ofNat_add]
  exact USize.toNat_ofNat_of_lt_32 (lattice_dimension_fits_word count columns hcount hcolumns).1

/-- Apply the width bound to the actual initialized and extended matrix;
there is no separate abstract matrix or dimension assumption. -/
theorem initialized_buildCldMatrix_size_fits_word
    (cld : Array SparsePolyZZ) (target : Nat) (scale modulus : ZZ)
    (matrix : Generated.StrictRecombine.LLLMatrix)
    (output : Generated.StrictRecombine.LLLMatrix × Nat)
    (hcount : cld.size < 2 ^ 31) (htarget : target < 2 ^ 31)
    (hscale : scale ≠ 0)
    (hinit : Generated.StrictRecombine.makeInitialMatrix cld.size scale = .ok matrix)
    (hbuild : Generated.StrictVanHoeij.buildCldMatrix matrix cld 0 target modulus = .ok output) :
    output.1.size < 2 ^ 32 ∧ output.1.size < 2 ^ 64 ∧
      ∀ row (hr : row < output.1.size), output.1[row].size < 2 ^ 64 := by
  have hc := buildCldMatrix_canonical cld 0 target scale modulus matrix output
    (makeInitialMatrix_canonical cld _ scale modulus matrix hscale hinit) hbuild
  have ha := buildCldMatrix_added_le cld 0 target modulus matrix output hbuild
  have hb := lattice_dimension_fits_word cld.size output.2 hcount (by omega)
  have hs : output.1.size = cld.size + output.2 := by simpa using hc.size_eq
  refine ⟨by omega, by omega, ?_⟩
  intro row hr
  have hh := hc.row_size row hr
  rw [getElem!_pos output.1 row hr] at hh
  omega

/-- Exact successor for each bounded container loop. -/
theorem lattice_index_word_increment (index dimension : Nat)
    (hindex : index < dimension) (hdimension : dimension < 2 ^ 64) :
    (UInt64.ofNat index + 1).toNat = index + 1 := by
  change (UInt64.ofNat index + UInt64.ofNat 1).toNat = index + 1
  rw [← UInt64.ofNat_add]
  exact UInt64.toNat_ofNat_of_lt (by change index + 1 < 2 ^ 64; omega)

/-- The unsigned countdown visits exactly k-2,...,0 without underflow. -/
theorem lattice_index_word_predecessor (remaining : Nat)
    (hpositive : 0 < remaining) (hfits : remaining < 2 ^ 64) :
    (UInt64.ofNat remaining - 1).toNat = remaining - 1 := by
  rw [UInt64.toNat_sub_of_le]
  · rw [UInt64.toNat_ofNat_of_lt' hfits]
    rfl
  · change (1 : Nat) ≤ (UInt64.ofNat remaining).toNat
    rw [UInt64.toNat_ofNat_of_lt' hfits]
    omega

end Refinement.StrictRecombine
