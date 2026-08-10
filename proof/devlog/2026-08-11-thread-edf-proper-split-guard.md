# Thread the concrete EDF proper-split guard

The generated `SplitState` already stores the successful C++ guard
`0 < get_deg g && get_deg g < get_deg f`, but the old `EDFRawOps.splitStep`
interface discarded it before asking for recursive invariant preservation and
measure decrease.

The generator now passes `split.proper` directly to `splitStep`. This is the
source execution fact needed to prove both recursive degrees decrease; it
prevents an implementation from hiding that obligation inside an unrelated
provider assumption.

Verification:

```text
python3 ../cpp2lean_v2/tests/build_strict_edf.py --check
lake build CLPoly.Refinement.EDF
```
