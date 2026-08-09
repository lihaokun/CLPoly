# HGCD physical continuation packages

## Date

2026-08-09

## What changed

- Added separate physical packages for the successful early and non-early
  continuations of the generated HGCD recursive body.
- Recorded the actual first dispatch, reconstruction, middle division,
  second dispatch, and finish executions in the appropriate package.
- Added wrapper theorems that attach the first and second strictly-smaller
  recursive semantic hypotheses to those physical packages and derive the
  common raw HGCD invariant.
- Added a dependent physical sum over the early and non-early continuations,
  with a branch-selected second-child obligation and one uniform raw-invariant
  theorem.
- Added a uniform theorem for the complete branch-admissible body.  It selects
  the base workspace or the successful non-base continuation from the actual
  generated guard and result.
- Added the proof-erasing adapter from an ordinary recursive callback to a
  strictly-smaller callback, the genuine generated recursive-equation
  predicate, and execution equalities back to both the source body and any
  callback satisfying that equation.
- Added a reusable physical invocation node and a theorem transferring its
  branch-local single-step invariant to an actual source-recursive call.
- Added the hereditary smaller-call workspace provider and the complete
  well-founded source-call refinement theorem.
- Strengthened the pinned-source gate so deletion of any physical
  continuation, source-call equality, invocation workspace, or final
  well-founded theorem fails CI.
- Removed stale translator documentation that described the GCD route as a
  mathematical-spec/L1-placeholder path; it now points exclusively at the
  strict raw-heap GCD/HGCD closure.

## Why

The top-level well-founded proof must case-split on the generated execution
and recover branch-local memory evidence without accepting the desired
polynomial semantics as an oracle.  These packages provide that physical
boundary while keeping recursive semantics exclusively in the induction
hypotheses.

## Key decisions

- The packages contain no `HgcdRecursiveCallbackRefinesAt` field.
- The early package has only the first-child semantic hypothesis at its
  wrapper theorem.
- The non-early package requires exactly two semantic hypotheses, each
  indexed by the actual smaller input length used by its generated call.
- The dependent sum reduces the second-child obligation to `True` on the
  actual early branch; no unused recursive premise is demanded there.
- The non-base package must prove that its reconstructed mathematical inputs
  are exactly the caller's represented polynomials; an arbitrary preselected
  specification value cannot be substituted.
- The decrease witnesses passed through the adapter are erased (`rfl` at an
  actual call), so this bridge adds no executable counter or alternate path.
- Invocation nodes quantify over the first-child induction theorem instead of
  storing one.  Consequently, recursive semantics still have to be supplied
  by the forthcoming strong induction, while all memory evidence remains
  reusable and non-semantic.
- `hgcdRecursiveCall_rawInvariant_wf` recursively discharges both generated
  child calls.  Its termination measure is the actual `lenA`; the two
  decrease obligations are the source-derived first-child bound and
  `middle.lenC0 < lenA`.  No fuel or executable termination guard exists.
- No fuel parameter or L2 fallback was introduced.

## Problems and resolution

The first-child decrease proof also needs the non-base guard.  The guard is
therefore an index of both continuation packages rather than an untracked
assumption reconstructed later.

## Files

- `proof/lean/CLPoly/Impl/StrictHGCDRawRefinement.lean`
- `proof/cpp2lean_v2/tests/check_strict_hgcd_source.py`
- `proof/cpp2lean_v2/tests/build_strict_gcd.py`
- `proof/cpp2lean_v2/class_map.py`
- `docs/devlog/2026-08-09-hgcd-physical-continuation-packages.md`
