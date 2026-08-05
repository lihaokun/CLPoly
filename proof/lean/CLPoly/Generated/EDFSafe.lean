/- Total, fuel-bounded counterpart of the generated EDF control flow. -/
import CLPoly.Generated.Corpus

set_option autoImplicit false

namespace Generated

/-- Total counterpart of `__upoly_mod_ir`. -/
def upolyModSafe (f g : SparsePolyZp) : SparsePolyZp :=
  (SparsePolyZp.divmod f g).2

/-- Generate coefficients for degrees `degree, degree - 1, ...`. -/
def upolyRandomLoopSafe : Nat → Int64 → SparsePolyZp → UInt64 → Rng →
    SparsePolyZp × Rng
  | 0, _, result, _, rng => (result, rng)
  | count + 1, degree, result, p, rng =>
      -- `Rng.next` uses an exclusive upper bound; C++ samples `[0, p - 1]`.
      let next := Rng.next_advance rng p
      let result' := if next.1 != 0 then
        result.push (UMonomial.mk degree, Zp.ofInt next.1.toInt p)
      else result
      upolyRandomLoopSafe count (degree - 1) result' p next.2

/-- Total counterpart of `__upoly_random_ir` for a nonnegative maximum degree. -/
def upolyRandomSafe (maxDeg : Int64) (p : UInt64) (rng : Rng) :
    SparsePolyZp × Rng :=
  if maxDeg > 0 then
    upolyRandomLoopSafe maxDeg.toNat (maxDeg - 1) #[] p rng
  else
    (#[], rng)

/-- Characteristic-two trace loop, expressed by the exact number of remaining
iterations rather than a wrapping `UInt64` counter. -/
def edfTraceLoopSafe : Nat → SparsePolyZp → SparsePolyZp → SparsePolyZp →
    SparsePolyZp
  | 0, g, _, _ => g
  | n + 1, g, f, r =>
      edfTraceLoopSafe n (upolyModSafe (g * g + r) f) f r

/-- Total canonical counterpart of `__upoly_subtract_one_ir`. -/
def upolySubtractOneSafe (h : SparsePolyZp) (p : UInt64) : SparsePolyZp :=
  h - #[(UMonomial.mk (0 : Int32), Zp.ofInt (1 : Int) p)]

/-- One generated Cantor–Zassenhaus attempt. `none` means that the random
candidate was empty or produced a trivial gcd. -/
def edfSplitAttemptSafe (f : SparsePolyZp) (d : UInt64) (rng : Rng) :
    Option SparsePolyZp × Rng :=
  let p := (SparsePolyZp.front! f).snd.prime
  let sample := upolyRandomSafe (get_deg f) p rng
  let r := sample.1
  let rng' := sample.2
  if r.isEmpty then (none, rng')
  else
    let g := if p == 2 then
      let g0 := upolyModSafe r f
      polynomial_GCD (edfTraceLoopSafe (d.toNat - 1) g0 f r) f
    else
      let exponent : Int := ((p.toNat : Int) ^ d.toNat - 1) / 2
      let pow := upolyPowmodSafe upolyModSafe r exponent f
      polynomial_GCD (upolySubtractOneSafe pow p) f
    if get_deg g > 0 && get_deg g < get_deg f then (some g, rng') else (none, rng')

/-- Fuel-bounded totalization of `__edf_Zp_ir`. Every retry and every recursive
split consumes one unit of fuel. -/
def edfZpTrySafe : Nat → Array SparsePolyZp → SparsePolyZp → UInt64 → Rng →
    Option (Array SparsePolyZp × Rng)
  | fuel, result, f, d, rng =>
      if (get_deg f).toUInt64 == d then
        some (result.push (Generated.__upoly_make_monic_ir_def f).2, rng)
      else if get_deg f <= 0 then
        some (result, rng)
      else match fuel with
        | 0 => none
        | remaining + 1 =>
          let attempt := edfSplitAttemptSafe f d rng
          match attempt.1 with
          | none => edfZpTrySafe remaining result f d attempt.2
          | some g =>
            let quotient := pair_vec_div #[] f g (SparsePolyZp.comp f)
            let h := SparsePolyZp.normalization quotient
            let gMonic := (Generated.__upoly_make_monic_ir_def g).2
            let hMonic := (Generated.__upoly_make_monic_ir_def h).2
            match edfZpTrySafe remaining result gMonic d attempt.2 with
            | none => none
            | some (result', rng'') =>
              edfZpTrySafe remaining result' hMonic d rng''

end Generated
