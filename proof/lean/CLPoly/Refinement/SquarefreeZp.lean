/-
  Strict squarefree refinement boundary.

  The earlier proof operated on legacy sparse-array entries whose division and
  GCD primitives could reduce definitionally to hand-written L2 algorithms.
  It therefore did not prove the raw C++ dense-polynomial execution and is not
  exported.

  Squarefree refinement will resume after strict RawHeap division and Euclidean
  GCD are connected to their L2 polynomial specifications.  All recursive
  source loops must remain well-founded; no bounded execution wrapper is part
  of this boundary.
-/
import CLPoly.Algorithm.SquarefreeZp
import CLPoly.Impl.StrictPolynomialGCDRefinement

set_option autoImplicit false

namespace Refinement

namespace StrictSquarefreeZp

open CLPoly.Impl.StrictPolynomialGCDRefinement
open CLPoly.Math

/-- Exact checked entry for the source `__upoly_make_monic`. -/
def upolyMakeMonicIR (f : SparsePolyZp) : RawExec (Zp × SparsePolyZp) :=
  if hnonempty : 0 < f.size then
    let lc := f[0].2
    if lc.val == 1 then
      .ok (lc, f)
    else
      let lcInv := Zp.inv lc
      .ok (lc, sparseMonicLoop 0 f lcInv)
  else
    .error .assertionFailure

/-- A concrete monic sparse input takes the actual early-return comparison in
`__upoly_make_monic`; no mutation loop is erased from the definition. -/
theorem upolyMakeMonicIR_eq_of_monic (p : Nat) [Fact (Nat.Prime p)]
    (f : SparsePolyZp) (hcanonical : SparsePolyZp.Canonical p f)
    (hnonempty : 0 < f.size) (hmonic : (SparsePolyZp.toPoly p f).Monic) :
    upolyMakeMonicIR f = .ok (f[0].2, f) := by
  have hlead := sparse_leading_val_eq_one_of_monic p f hcanonical hnonempty
    hmonic
  simp [upolyMakeMonicIR, hnonempty, hlead]

/-- Exact source-order range-for loop of `__extract_pth_root`. -/
def pthRootTerm (prime : UInt64) (term : UMonomial × Zp) : UMonomial × Zp :=
  (UMonomial.mk ((term.1.deg.toUInt64 / prime).toInt64), term.2)

def extractPthRootLoop (index : Nat) (output source : SparsePolyZp)
    (prime : UInt64) : SparsePolyZp :=
  if h : index < source.size then
    let term := source[index]
    let next := output.push (pthRootTerm prime term)
    extractPthRootLoop (index + 1) next source prime
  else
    output
termination_by source.size - index
decreasing_by omega

theorem extractPthRootLoop_toList (index : Nat)
    (output source : SparsePolyZp) (prime : UInt64) :
    (extractPthRootLoop index output source prime).toList =
      output.toList ++ (source.toList.drop index).map (pthRootTerm prime) := by
  rw [extractPthRootLoop]
  split
  next hmore =>
    rw [extractPthRootLoop_toList]
    have hdrop := List.drop_eq_getElem_cons
      (l := source.toList) (i := index) (by simpa using hmore)
    rw [hdrop, List.map_cons, Array.getElem_toList hmore]
    simp [List.append_assoc]
  next hdone =>
    have hindex : source.toList.length ≤ index := by simpa using hdone
    rw [List.drop_eq_nil_iff.mpr hindex]
    simp
termination_by source.size - index
decreasing_by omega

def extractPthRootIR (f : SparsePolyZp) : RawExec SparsePolyZp :=
  if hnonempty : 0 < f.size then
    .ok (extractPthRootLoop 0 #[] f f[0].2.prime)
  else
    .error .assertionFailure

theorem extractPthRootLoop_size (index : Nat) (output source : SparsePolyZp)
    (prime : UInt64) :
    (extractPthRootLoop index output source prime).size =
      output.size + (source.size - index) := by
  rw [extractPthRootLoop]
  split
  next hmore =>
    rw [extractPthRootLoop_size]
    simp
    omega
  next hdone =>
    simp
    omega
termination_by source.size - index
decreasing_by omega

/-- Exact result-copy loop used after recursive SQF contraction. -/
def scaleMultiplicityLoop (index : Nat)
    (source output : Array (SparsePolyZp × UInt64)) (prime : UInt64) :
    Array (SparsePolyZp × UInt64) :=
  if h : index < source.size then
    let item := source[index]
    scaleMultiplicityLoop (index + 1) source
      (output.push (item.1, item.2 * prime)) prime
  else
    output
termination_by source.size - index
decreasing_by omega

theorem scaleMultiplicityLoop_toList (index : Nat)
    (source output : Array (SparsePolyZp × UInt64)) (prime : UInt64) :
    (scaleMultiplicityLoop index source output prime).toList =
      output.toList ++ (source.toList.drop index).map
        (fun item => (item.1, item.2 * prime)) := by
  rw [scaleMultiplicityLoop]
  split
  next hmore =>
    rw [scaleMultiplicityLoop_toList]
    have hdrop := List.drop_eq_getElem_cons
      (l := source.toList) (i := index) (by simpa using hmore)
    rw [hdrop, List.map_cons, Array.getElem_toList hmore]
    simp [List.append_assoc]
  next hdone =>
    have hindex : source.toList.length ≤ index := by simpa using hdone
    rw [List.drop_eq_nil_iff.mpr hindex]
    simp
termination_by source.size - index
decreasing_by omega

theorem scaleMultiplicityLoop_size (index : Nat)
    (source output : Array (SparsePolyZp × UInt64)) (prime : UInt64) :
    (scaleMultiplicityLoop index source output prime).size =
      output.size + (source.size - index) := by
  rw [scaleMultiplicityLoop]
  split
  next hmore =>
    rw [scaleMultiplicityLoop_size]
    simp
    omega
  next hdone =>
    simp
    omega
termination_by source.size - index
decreasing_by omega

end StrictSquarefreeZp
end Refinement
