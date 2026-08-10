/-
Execution facts for the generated C++ `__upoly_random` loop.
-/
import CLPoly.Generated.StrictEDF

set_option autoImplicit false

namespace Refinement.StrictEDF

/-- The exact generated random-polynomial entry is total.  Its output and
advanced RNG state are the values computed by the well-founded source loop. -/
theorem __upoly_random_raw_ir_terminates
    (maxDegree : Int64) (p : UInt64) (rng : Rng) :
    ∃ output rngNext,
      Generated.StrictEDF.__upoly_random_raw_ir maxDegree p rng =
        .ok (output, rngNext) := by
  unfold Generated.StrictEDF.__upoly_random_raw_ir
  split
  · exact ⟨_, _, rfl⟩
  · exact ⟨#[], rng, rfl⟩

end Refinement.StrictEDF
