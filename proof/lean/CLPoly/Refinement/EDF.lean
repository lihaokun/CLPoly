/-
  Strict EDF refinement is intentionally not exported yet.

  The former contents of this module selected an L2 existence witness when a
  bounded generated run failed or when its result did not satisfy
  `EDFCorrect`.  That construction was a certified specification fallback, not
  a refinement of the C++ `__edf_Zp` execution, and has therefore been removed.

  A future theorem in this module must start from a cpp2lean-generated EDF trace
  (including RNG state transitions), prove every successful split and recursive
  output invariant, and expose failure/nontermination explicitly.  It must not
  manufacture an L2 result on an unproved execution branch.
-/
import CLPoly.Algorithm.EDF
import CLPoly.Generated.Corpus
import CLPoly.Refinement.Basic

set_option autoImplicit false

namespace Refinement

-- No L1→L2 EDF theorem is exported until the generated execution proof closes.

end Refinement
