/-
  Strict Z/p arithmetic refinement boundary.

  The legacy Corpus polynomial operations are expressed through typeclass
  dispatch to hand-written `SparsePolyZp` operations.  Equality obtained by
  unfolding that dispatch is not evidence that the corresponding C++ loops
  execute correctly, so this module exports no L1→L2 polynomial theorem.

  Concrete scalar and dense-polynomial proofs are developed from
  `Generated.StrictGCD` and `Generated.StrictDivrem`, then lifted through
  `RawHeap` observations in `CLPoly.Impl`.
-/
import CLPoly.Generated.StrictGCD
import CLPoly.Generated.StrictDivrem
import CLPoly.Impl.RawPolynomialRep

set_option autoImplicit false

namespace Refinement

-- No legacy dispatch-based polynomial refinement theorem is exported.

end Refinement
