-- Auto-generated strict C++ control flow for `polynomial_mod`.
import CLPoly.Model

set_option autoImplicit false

namespace Generated.StrictPolynomialMod

/-- Source range-for: reduce one integer coefficient and omit zero residues. -/
def polynomialModLoop (f : SparsePolyZZ) (p : UInt64) (index : Nat)
    (result : SparsePolyZp) : SparsePolyZp :=
  if h : index < f.size then
    let term := f[index]
    let coefficient := Zp.ofInt term.2 p
    let result' := if coefficient.val = 0 then result
      else result.push (term.1, coefficient)
    polynomialModLoop f p (index + 1) result'
  else result
termination_by f.size - index
decreasing_by omega

/-- Strict generated C++ entry. -/
def polynomial_mod_raw_ir (f : SparsePolyZZ) (p : UInt64) :
    RawExec SparsePolyZp :=
  .ok (polynomialModLoop f p 0 #[])

end Generated.StrictPolynomialMod
