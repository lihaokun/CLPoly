/-
  Strict Hensel refinement boundary.

  This module deliberately exports no candidate from the legacy generated
  corpus.  The former safe wrapper selected the L2 `hensel_lift`
  implementation whenever the candidate failed a post-hoc `HenselCorrect`
  check.  That was not an L1→L2 refinement and has been removed.
-/
import CLPoly.Algorithm.Hensel
import CLPoly.Refinement.Basic

set_option autoImplicit false

namespace Refinement

-- No Hensel L1→L2 theorem or legacy candidate is exported until a strict
-- cpp2lean-generated entry and its direct execution proof are available.

end Refinement
