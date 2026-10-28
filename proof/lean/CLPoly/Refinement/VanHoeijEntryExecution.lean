import CLPoly.Generated.StrictVanHoeijEntry
import CLPoly.Refinement.VanHoeijControllerExecution
import CLPoly.Refinement.VanHoeijParameters
import CLPoly.Refinement.VanHoeijWidth
import CLPoly.Refinement.VanHoeijTruncation

set_option autoImplicit false

namespace Refinement.StrictRecombine

/-- Concrete proofs discharge every entry obligation; none can choose data. -/
def concreteVanHoeijEntryExecution (termination : Generated.StrictHensel.DivmodTermination) :
    Generated.StrictVanHoeijEntry.EntryExecution where
  controller := concreteVanHoeijControllerExecution
  cld := ⟨termination⟩
  cldSize := by
    intro f lifted modulus output hrun
    exact cldPolys_size_of_success ⟨termination⟩ f lifted output modulus hrun
  truncateSize := by
    intro input divisor output hrun
    exact (truncateCld_spec input output divisor hrun).2.1
  parameterModulus := by
    intro count height modulus parameters hrun
    have hh := (selectParameters_spec count height modulus parameters hrun).2.2.2.2.2.1
    dsimp only [ZZ] at *
    omega
  initialBound := initialCldTarget_bound

end Refinement.StrictRecombine
