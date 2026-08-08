/-
  Strict integer-arithmetic refinement boundary.

  The former theorem in this module unfolded `Generated.__symmetric_mod_ir`
  from the legacy aggregate corpus.  That corpus is not an admissible source
  for the strict C++ L1→Lean L2 dependency closure: unrelated entries may
  contain partial or specification-dispatched definitions.  Scalar integer
  refinement must instead start from a separately generated, reproducible,
  placeholder-free C++ root, just like `Generated.StrictGCD`.
-/
import CLPoly.Refinement.Basic

set_option autoImplicit false

namespace Refinement

-- No integer L1→L2 theorem is exported until its exact C++ root is admitted
-- to a strict generated module and proved from that generated execution.

end Refinement
