-- Auto-generated strict raw C++ control flow for
-- `__lll_factorize` and `__factor_squarefree_primitive_ZZ`.
import CLPoly.Model
import CLPoly.Generated.StrictHensel

set_option autoImplicit false

namespace Generated.StrictFactorZZ

/-- Concrete effectful callees of the two C++ entries.  These fields expose
only raw executions and cannot manufacture an L2 factorization witness. -/
structure FactorZZRawOps where
  selectPrime : SparsePolyZZ → Bool → RawExec PrimeSelectionResult
  henselLift : SparsePolyZZ → Array SparsePolyZp → UInt64 → Int32 →
    RawExec (Array SparsePolyZZ × ZZ)
  vanHoeijRecombine : SparsePolyZZ → Array SparsePolyZZ → ZZ →
    RawExec (Array SparsePolyZZ)
  zassenhausRecombine : SparsePolyZZ → Array SparsePolyZZ → ZZ →
    RawExec (Array SparsePolyZZ)

/-- Exact safety-net decision shared with the C++ top-level controller. -/
def __needs_zassenhaus_safety_net_ir (resultCount modularFactorCount : Nat)
    (atFullPrecision : Bool) : Bool :=
  atFullPrecision && resultCount < modularFactorCount

/-- Exact source power loop from `__heuristic_starting_precision`.  The C++
caller supplies a prime, hence `p ≥ 2`; making that source precondition
explicit gives the genuine decreasing distance `target + 1 - pa`. -/
def heuristicPrecisionLoop (p target : Nat) (hp : 2 ≤ p) :
    (pa : Nat) → 0 < pa → Nat → Nat
  | pa, hpa, exponent =>
      if hcontinue : pa ≤ target then
        heuristicPrecisionLoop p target hp (pa * p) (Nat.mul_pos hpa (by omega))
          (exponent + 1)
      else exponent
termination_by pa _ _ => target + 1 - pa
decreasing_by
  have hdouble : pa + pa ≤ pa * p := by
    calc
      pa + pa = 2 * pa := by omega
      _ ≤ p * pa := Nat.mul_le_mul_right pa hp
      _ = pa * p := Nat.mul_comm _ _
  have hgrow : pa < pa * p := by omega
  omega

/-- Strict generated lowering of C++ `__heuristic_starting_precision`.  It
executes the floating-point heuristic and the concrete Mignotte-bound helper,
then computes the full-precision exponent with the well-founded source loop
above. -/
def __heuristic_starting_precision_raw_ir (f : SparsePolyZZ) (r : Int32)
    (p : UInt64) : RawExec (Int32 × Int32) :=
  if hp : 2 ≤ p.toNat then
    match f[0]? with
    | none => .error .assertionFailure
    | some leading =>
      match Generated.StrictHensel.__mignotte_bound_upoly_raw_ir f with
      | .error fault => .error fault
      | .ok bound =>
        let logp := Float.log (Nat.toFloat p.toNat)
        let minB : Int32 := (ZZ.sizeinbase (p.toNat : Int) 2).toUInt32.toInt32
        let degree : Int32 := f.size.toUInt32.toInt32 - 1
        let heuristicFloat := Float.ceil
          (((2.5 * Int.toFloat r + Int.toFloat minB) * Float.log 2.0 / logp) +
            (Float.log (Int.toFloat (degree + 1)) / (2.0 * logp)))
        let heuristic : Int32 := max 1 (Float.toInt32 heuristicFloat)
        let leadingAbs := if leading.2 < 0 then -leading.2 else leading.2
        let target := (2 * leadingAbs * bound).natAbs
        let fullNat := heuristicPrecisionLoop p.toNat target hp 1 (by omega) 0
        if hfit : fullNat ≤ UInt32.size - 1 then
          let full := fullNat.toUInt32.toInt32
          .ok (min full heuristic, full)
        else .error .arithmeticOverflow
  else .error .assertionFailure

/-- Exact source-shaped lowering of C++ `__lll_factorize`: heuristic
precision, first Hensel/recombine pass, and the conditional full-precision
second pass. -/
def __lll_factorize_raw_ir (ops : FactorZZRawOps) (f : SparsePolyZZ)
    (factors : Array SparsePolyZp) (p : UInt64) :
    RawExec (Array SparsePolyZZ) :=
  let r : Int32 := factors.size.toUInt32.toInt32
  match __heuristic_starting_precision_raw_ir f r p with
  | .error fault => .error fault
  | .ok (aH, aMig) =>
    match ops.henselLift f factors p aH with
    | .error fault => .error fault
    | .ok (liftedH, mH) =>
      match ops.vanHoeijRecombine f liftedH mH with
      | .error fault => .error fault
      | .ok result =>
        if result.size.toUInt32.toInt32 < r && aH < aMig then
          match ops.henselLift f factors p 0 with
          | .error fault => .error fault
          | .ok (liftedMig, mMig) =>
            match ops.vanHoeijRecombine f liftedMig mMig with
            | .error fault => .error fault
            | .ok resultMig =>
              if __needs_zassenhaus_safety_net_ir resultMig.size factors.size
                  true then
                ops.zassenhausRecombine f liftedMig mMig
              else .ok resultMig
        else
          if __needs_zassenhaus_safety_net_ir result.size factors.size
              (aMig ≤ aH) then
            ops.zassenhausRecombine f liftedH mH
          else .ok result

/-- Exact source-shaped lowering of C++
`__factor_squarefree_primitive_ZZ`: select a prime, return the input on the
irreducible/single-factor branch, otherwise execute `__lll_factorize`. -/
def __factor_squarefree_primitive_ZZ_raw_ir (ops : FactorZZRawOps)
    (useLargePrime : Bool) (f : SparsePolyZZ) : RawExec (Array SparsePolyZZ) :=
  if f.isEmpty || get_deg f < 2 then
    .error .assertionFailure
  else
    match ops.selectPrime f useLargePrime with
    | .error fault => .error fault
    | .ok selection =>
      if selection.irreducible || selection.factors.size ≤ 1 then
        .ok #[f]
      else
        __lll_factorize_raw_ir ops f selection.factors selection.prime

end Generated.StrictFactorZZ
