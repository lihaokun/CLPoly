import B2B.StrictRuntime

namespace B2B.TestStrictRuntime

unsafe def draws (count : Nat) (upper : UInt64)
    (state : StrictRuntime.MT19937State) :
    List UInt64 × StrictRuntime.MT19937State :=
  match count with
  | 0 => ([], state)
  | count + 1 =>
    let (word, next) := StrictRuntime.uniformBelow state upper
    let (tail, final) := draws count upper next
    (word :: tail, final)

def expected : List (UInt64 × List UInt64) := [
  (2, [0, 1, 0, 0, 0, 1, 0, 0, 0, 1]),
  (3, [2, 0, 2, 2, 0, 0, 2, 1, 2, 2]),
  (5, [3, 4, 2, 4, 4, 1, 2, 2, 2, 4]),
  (7, [6, 3, 4, 6, 2, 4, 4, 6, 1, 2]),
  (11, [6, 3, 10, 7, 4, 6, 9, 2, 6, 10]),
  (4294967311, [112825614, 2208691414, 2100207746, 220921140,
    100854103, 696226747]),
  (9223372036854775783, [6909045637428952499, 8314211556539077902,
    4279532810384561223, 1819927849474927636, 2878035897379592313,
    2877591057541362902]),
  (18446744073709551557, [6909045637428952499, 17537583593393853710,
    13502904847239337031, 11043299886329703444, 2878035897379592313,
    2877591057541362902])
]

unsafe def testCases : List (UInt64 × List UInt64) → IO Bool
  | [] => pure true
  | (upper, wanted) :: rest => do
    let actual := (draws wanted.length upper (StrictRuntime.mtSeed 42)).1
    if actual = wanted then
      IO.println s!"PASS: mt19937/uniform [0,{upper})"
      testCases rest
    else
      IO.eprintln s!"FAIL: mt19937/uniform [0,{upper}): {actual}"
      pure false

end B2B.TestStrictRuntime

unsafe def main : IO UInt32 := do
  if ← B2B.TestStrictRuntime.testCases B2B.TestStrictRuntime.expected then
    return 0
  else
    return 1
