/-
  CLPoly/Refinement/Recombine.lean — Recombine 精化：__factor_recombine_upoly_ir（C++） ≃ recombine_correct（L2）
-/
import CLPoly.Generated.Corpus
import CLPoly.Algorithm.Recombine
import CLPoly.Refinement.Basic
import CLPoly.Model
import CLPoly.Math.Univariate

set_option autoImplicit false

open Refinement
open CLPoly.Math

noncomputable def recombine_l2 (f : Polynomial ℤ) (hf : f ≠ 0) : List (Polynomial ℤ) :=
  (recombine_correct f hf).choose

theorem __factor_recombine_upoly_ir_refines (f : Polynomial ℤ) (hf : f ≠ 0)
    (f_sparse : SparsePolyZZ) (lifted : Array SparsePolyZZ) (m : ZZ)
    (h_f_sparse : SparsePolyZZ.toPoly f_sparse = f)
    : ((Generated.__factor_recombine_upoly_ir f_sparse lifted m).map SparsePolyZZ.toPoly).toList =
      recombine_l2 f hf := by
  sorry
