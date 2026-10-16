import CLPoly.Generated.StrictVanHoeij

set_option autoImplicit false

namespace Refinement.StrictRecombine
open Generated.StrictVanHoeij

theorem buildCldMatrixLoop_added_le (cld : Array SparsePolyZZ)
    (current target width : Nat) (modulus : ZZ) (added : Nat)
    (matrix : Generated.StrictVanHoeij.LLLMatrix) (output : Generated.StrictVanHoeij.LLLMatrix × Nat)
    (hadded : added ≤ target)
    (hrun : buildCldMatrixLoop cld current target width modulus added matrix = .ok output) :
    output.2 ≤ target := by
  induction hmeasure : target - added using Nat.strong_induction_on generalizing added matrix with
  | h measure ih =>
    rw [buildCldMatrixLoop] at hrun
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
      next hw =>
        have he := Except.ok.inj hrun
        rw [← he]
        exact hadded
    next ha =>
      have he := Except.ok.inj hrun
      rw [← he]
      exact hadded

theorem buildCldMatrix_added_le (cld : Array SparsePolyZZ)
    (current target : Nat) (modulus : ZZ)
    (matrix : Generated.StrictVanHoeij.LLLMatrix) (output : Generated.StrictVanHoeij.LLLMatrix × Nat)
    (hrun : buildCldMatrix matrix cld current target modulus = .ok output) :
    output.2 ≤ target :=
  buildCldMatrixLoop_added_le cld current target _ modulus 0 matrix output (by omega) hrun

theorem buildCldMatrix_added_int32 (cld : Array SparsePolyZZ)
    (current : Nat) (target : Int32) (modulus : ZZ)
    (matrix : Generated.StrictVanHoeij.LLLMatrix) (output : Generated.StrictVanHoeij.LLLMatrix × Nat)
    (hrun : buildCldMatrix matrix cld current target.toInt.toNat modulus = .ok output) :
    (Int32.ofNat output.2).toInt = (output.2 : Int) := by
  have hb := buildCldMatrix_added_le cld current target.toInt.toNat modulus matrix output hrun
  have ht := Int32.toInt_lt target
  have hf : (output.2 : Int) < 2 ^ 31 := by omega
  exact Int32.toInt_ofNat_of_lt (by omega)

end Refinement.StrictRecombine
