# Strict Hensel quadratic-precision loop

## Scope

This stage closes only the source `while (m <= target)` segment inside
`__hensel_lift`.  It is intentionally not published as the final
`__hensel_lift` contract: target computation, leading-coefficient adjustment,
tree construction, factor extraction, and final monic normalization remain
separate source segments that must still be connected.

## Exact L1 execution

`Generated.StrictHensel.__hensel_lift_loop_raw_ir` executes the already strict
`__hensel_lift_recursive_raw_ir` traversal and only then updates `m` to
`m * m`.  The loop uses the nonnegative `Nat` view of the source `ZZ` modulus
and target.  Its `m >= 2` proof is erased termination evidence, not a fuel
counter.

Termination is well founded on `target + 1 - m`.  On every continuing
iteration, `m >= 2` implies `m < m * m`, so the measure decreases strictly.

## Refinement evidence

- `HenselLiftLoopPrefix` records only exact generated traversal equations.
- The execution and refinement invariants quantify over every state reachable
  by such a prefix, so they cannot select an expected node array.
- `HenselLiftLoopCorrect` records the source guard, exact traversal run,
  `HenselLiftRecursiveCorrect` proof, squaring update, and remaining loop trace.
- `__hensel_lift_loop_raw_ir_terminates` is the raw-to-safe bridge.
- `__hensel_lift_loop_raw_ir_refines` is the strict semantic refinement theorem
  for this source loop segment.

There is no fuel, specification oracle, L2 execution fallback, or partial
definition in this closure.
