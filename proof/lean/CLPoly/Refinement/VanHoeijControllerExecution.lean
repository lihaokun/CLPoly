import CLPoly.Generated.StrictVanHoeijController
import CLPoly.Refinement.VanHoeijControllerArithmetic
import CLPoly.Refinement.VanHoeijLattice

set_option autoImplicit false

namespace Refinement.StrictRecombine

theorem column_target_lt (target width : Int32) (hne : target ≠ width)
    (hbound : target.toInt ≤ width.toInt) : target < width := by
  apply Int32.lt_iff_toInt_lt.mpr
  have hneq : target.toInt ≠ width.toInt := by
    intro heq
    apply hne
    exact Int32.toInt_inj.mp heq
  omega

/-- All erased controller obligations discharged by existing concrete proofs. -/
def concreteVanHoeijControllerExecution : Generated.StrictVanHoeijController.ControllerExecution where
  lll := concreteLLLExecution
  buildValid := by
    intro cld target scale modulus initial built hs hm hi hb
    exact initialized_buildCldMatrix_input_valid cld target scale modulus initial built hs hm hi hb
  targetBound := by
    intro target width ht hne hb
    exact (nextCldTarget_progress target width ht (column_target_lt target width hne hb)).2
  targetDecreases := by
    intro target width ht hne hb
    exact nextCldTarget_decreases_remaining target width ht
      (column_target_lt target width hne hb)

end Refinement.StrictRecombine
