/-
  Refinement/Basic.lean — L1→L2 精化证明的基础 bridge 引理

  本文件是手写可信基，不被 Pass 9 覆盖。

  精化定理通常形如：

    theorem <l1_ir>_refines (f : SparsePolyZp) (hwf : WellFormed p f) ... :
      (l1_ir f).map (SparsePolyZp.toPoly p ∘ Prod.fst, ...) = l2 (SparsePolyZp.toPoly p f) := by
      <direct generated-execution proof>
-/

import CLPoly.Model
import CLPoly.Generated.Corpus
import CLPoly.Math.Univariate

set_option autoImplicit false

open Polynomial
open CLPoly.Math

namespace Refinement

noncomputable section BridgeLemmas

/-- Array (SparsePolyZp × UInt64) → List (Polynomial (ZMod p) × ℕ) 的映射。
    L1 返回 Array (SparsePolyZp × UInt64)，L2 返回 List (Polynomial (ZMod p) × ℕ)。 -/
noncomputable def toPolyList (arr : Array (SparsePolyZp × UInt64)) (p : ℕ) :
    List (Polynomial (ZMod p) × ℕ) :=
  (arr.map (fun (g, e) => (SparsePolyZp.toPoly p g, e.toNat))).toList

/-- toPolyList 在空数组上是空列表 -/
@[simp] theorem toPolyList_empty (p : ℕ) :
    toPolyList (#[] : Array (SparsePolyZp × UInt64)) p = [] := by
  simp [toPolyList]

end BridgeLemmas

end Refinement
