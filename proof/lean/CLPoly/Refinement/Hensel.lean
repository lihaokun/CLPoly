/-
  Strict Hensel refinement boundary.

  This module deliberately exports only the observation map for the generated
  C++ candidate.  The former safe wrapper selected the L2 `hensel_lift`
  implementation whenever the candidate failed a post-hoc `HenselCorrect`
  check.  That was not an L1→L2 refinement and has been removed.
-/
import CLPoly.Generated.Corpus
import CLPoly.Algorithm.Hensel
import CLPoly.Refinement.Basic

set_option autoImplicit false

namespace Refinement

/-- Map the integer representatives returned by generated Hensel lifting into
the ring in which the L2 result is observed.  This is only a representation
map; it makes no correctness claim about the generated candidate. -/
noncomputable def henselCandidateToPk (p k : Nat)
    (candidate : Array SparsePolyZZ × Int) :
    List (Polynomial (ZMod (p ^ k))) :=
  (candidate.1.map fun factor =>
    Polynomial.map (Int.castRingHom (ZMod (p ^ k)))
      (SparsePolyZZ.toPoly factor)).toList

/-- Observable output of the raw cpp2lean Hensel entry for an explicit target
exponent.  Correctness must be proved from the generated tree/lift execution;
there is no L2 fallback. -/
noncomputable def henselGeneratedCandidate (p k : Nat) (f : SparsePolyZZ)
    (factors : Array SparsePolyZp) : List (Polynomial (ZMod (p ^ k))) :=
  henselCandidateToPk p k
    (Generated.__hensel_lift_upoly_ir f factors (UInt64.ofNat p) (k : Int32))

end Refinement
