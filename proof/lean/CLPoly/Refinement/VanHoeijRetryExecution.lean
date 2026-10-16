import CLPoly.Generated.StrictVanHoeijRetry
import CLPoly.Refinement.VanHoeijEntryExecution
import CLPoly.Refinement.VanHoeijHenselCount

set_option autoImplicit false

namespace Refinement.StrictRecombine

/-- All operations execute the existing strict arithmetic implementations.
The lifted-count premise follows from the same successful Hensel execution. -/
def concreteVanHoeijRetryExecution (this : DenseUPolyZp)
    (termination : Generated.StrictHensel.DivmodTermination)
    (mulProvider : StrictDDF.RawMulWorkspaceProvider this) :
    Generated.StrictVanHoeijRetry.RetryExecution where
  step := StrictHensel.strictHenselRawOps termination
  tree := StrictHensel.strictHenselTreeBuildRawOps this mulProvider
  entry := concreteVanHoeijEntryExecution termination
  liftedCountFits := by
    intro f factors p target hp output hrun
    exact (StrictHensel.hensel_entry_success_count
      (StrictHensel.strictHenselRawOps termination)
      (StrictHensel.strictHenselTreeBuildRawOps this mulProvider)
      f factors p target hp output hrun).2
  precisionNonnegative := by
    intro f factors p target hp output hrun
    rcases StrictHensel.__hensel_lift_upoly_raw_ir_exponent_of_success
      (StrictHensel.strictHenselRawOps termination)
      (StrictHensel.strictHenselTreeBuildRawOps this mulProvider)
      f factors p target hp output hrun with hz | hp
    · subst target; exact Int32.le_refl _
    · exact Int32.le_iff_toInt_le.mpr (le_of_lt (Int32.lt_iff_toInt_lt.mp hp))

end Refinement.StrictRecombine
