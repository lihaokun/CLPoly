-- Auto-generated strict refinement contract for `__ddf_Zp_raw_ir`.
import CLPoly.Refinement.DDF

set_option autoImplicit false

namespace Refinement

open CLPoly.Math

/-- The cpp2lean-generated C++ DDF entry terminates and decodes to L2 `ddf`.
The double underscores are retained verbatim from the original C++ symbol. -/
theorem __ddf_Zp_raw_ir_refines_ddf
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (providers : StrictDDF.DDFRawProviders this) (f : SparsePolyZp)
    (hfPrime : f[0]!.2.prime = this._p)
    (hfCanonical : SparsePolyZp.Canonical this._p.toNat f)
    (hfDegree : ∀ term ∈ f.toList, term.1.deg < 2 ^ 62)
    (hfMonic : (SparsePolyZp.toPoly this._p.toNat f).Monic)
    (hfSquarefree : Squarefree (SparsePolyZp.toPoly this._p.toNat f)) :
    ∃ output,
      Generated.StrictDDF.__ddf_Zp_raw_ir
          (StrictDDF.strictDDFRawOps this providers
            (SparsePolyZp.toPoly this._p.toNat f)) f
          (fun _ => StrictDDF.DDFLoopInvariant.initial this f
            f[0]!.2.prime hfPrime
            hfCanonical hfDegree hfMonic hfSquarefree) = .ok output ∧
      StrictDDF.ddfResultToL2 this._p.toNat output =
        ddf (SparsePolyZp.toPoly this._p.toNat f) := by
  exact StrictDDF.strictDDFEntryIR_refines_ddf this providers f hfPrime
    hfCanonical hfDegree hfMonic hfSquarefree

end Refinement
