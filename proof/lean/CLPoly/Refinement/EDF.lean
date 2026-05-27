/-
  CLPoly/Refinement/EDF.lean — EDF 精化：__edf_Zp_ir（C++） ≃ edf（L2）

  将 `__edf_Zp_ir` 的输出映射到 L2 `edf` 输出。
-/
import CLPoly.Generated.Corpus
import CLPoly.Algorithm.EDF
import CLPoly.Refinement.Basic
import CLPoly.Model
import CLPoly.Math.Univariate

set_option autoImplicit false

open Refinement
open CLPoly.Math

theorem __edf_Zp_ir_refines (p : ℕ) [hp : Fact (Nat.Prime p)]
    (f : SparsePolyZp) (d : UInt64) (rng : Rng)
    (hwf_f : SparsePolyZp.WellFormed p f)
    (hred_f : SparsePolyZp.AllReduced p f.toList)
    (hp_size : 2 * p ≤ UInt64.size)
    : ((Generated.__edf_Zp_ir #[] f d rng).1.map (SparsePolyZp.toPoly p)).toList =
      edf (SparsePolyZp.toPoly p f) d.toNat [] := by
  sorry
