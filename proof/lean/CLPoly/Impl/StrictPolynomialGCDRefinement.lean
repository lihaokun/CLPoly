import CLPoly.Impl.StrictGCDHGCDRefinement
import CLPoly.Math.Univariate

set_option autoImplicit false

namespace CLPoly.Impl.StrictPolynomialGCDRefinement

open CLPoly.Impl.StrictEuclidRefinement
open CLPoly.Impl.StrictDivremRefinement
open CLPoly.Impl.RawPolynomialRep
open CLPoly.Math

/-- Exact raw lowering of the reverse coefficient scan in
`dense_upoly_zp::to_upoly`. -/
def denseToSparseLoop (this : DenseUPolyZp) (heap : RawHeap)
    (ptr : RawPtr UInt64) : Nat → SparsePolyZp → RawExec SparsePolyZp
  | 0, output => .ok output
  | remaining + 1, output =>
      match heap.readU64 ptr remaining with
      | .error fault => .error fault
      | .ok value =>
          let output' := if value = 0 then output else
            output.push (⟨remaining⟩, ⟨value, this._p⟩)
          denseToSparseLoop this heap ptr remaining output'

/-- Exact raw object conversion corresponding to `to_upoly()`. -/
def dense_upoly_zp_to_upoly_raw_ir (this : DenseUPolyZp) (heap : RawHeap)
    (ptr : RawPtr UInt64) (length : Nat) : RawExec SparsePolyZp :=
  denseToSparseLoop this heap ptr length #[]

/-- Every read performed by the exact dense-to-sparse scan is in bounds when
the source coefficient vector is valid. -/
theorem denseToSparseLoop_succeeds (this : DenseUPolyZp) (heap : RawHeap)
    (ptr : RawPtr UInt64) (length remaining : Nat) (output : SparsePolyZp)
    (hvalid : heap.ValidU64Slice ptr length) (hremaining : remaining ≤ length) :
    ∃ result, denseToSparseLoop this heap ptr remaining output = .ok result := by
  induction remaining generalizing output with
  | zero => exact ⟨output, rfl⟩
  | succ remaining ih =>
      rcases heap.readU64_of_valid ptr length remaining hvalid (by omega) with
        ⟨value, hread⟩
      by_cases hzero : value = 0
      · rcases ih output (by omega) with ⟨result, hrun⟩
        exact ⟨result, by
          simp [denseToSparseLoop, hread, hzero, hrun]⟩
      · rcases ih (output.push (⟨remaining⟩, ⟨value, this._p⟩)) (by omega) with
          ⟨result, hrun⟩
        exact ⟨result, by
          simp [denseToSparseLoop, hread, hzero, hrun]⟩

/-- The exact C++ dense-to-sparse conversion cannot fault on a valid dense
coefficient allocation. -/
theorem dense_upoly_zp_to_upoly_raw_ir_succeeds (this : DenseUPolyZp)
    (heap : RawHeap) (ptr : RawPtr UInt64) (length : Nat)
    (hvalid : heap.ValidU64Slice ptr length) :
    ∃ result, dense_upoly_zp_to_upoly_raw_ir this heap ptr length = .ok result :=
  denseToSparseLoop_succeeds this heap ptr length length #[] hvalid le_rfl

theorem sparseToPoly_push_raw (p : UInt64) (output : SparsePolyZp)
    (degree : Nat) (value : UInt64) :
    SparsePolyZp.toPoly p.toNat
        (output.push (⟨degree⟩, ⟨value, p⟩)) =
      SparsePolyZp.toPoly p.toNat output +
        Polynomial.monomial degree (value.toNat : ZMod p.toNat) := by
  simp [SparsePolyZp.toPoly, Array.toList_push, listSum_append, listSum,
    Zp.toZMod, add_comm]

/-- Semantic invariant of the reverse scan: scanning a represented low
prefix appends exactly that prefix polynomial after the already accumulated
higher sparse terms. -/
theorem denseToSparseLoop_refines (this : DenseUPolyZp) (heap : RawHeap)
    (ptr : RawPtr UInt64) (remaining : Nat) (output : SparsePolyZp)
    (poly : Polynomial (ZMod this._p.toNat))
    (hvalid : heap.ValidU64Slice ptr remaining)
    (hrep : SlicePolyRep heap ptr remaining this._p.toNat poly) :
    ∃ result,
      denseToSparseLoop this heap ptr remaining output = .ok result ∧
      SparsePolyZp.toPoly this._p.toNat result =
        SparsePolyZp.toPoly this._p.toNat output + poly := by
  induction remaining generalizing output poly with
  | zero =>
      have hpoly : poly = 0 :=
        slicePolyRep_zero_length heap ptr this._p.toNat poly hrep
      subst poly
      exact ⟨output, rfl, by simp⟩
  | succ remaining ih =>
      have hvalidPrefix : heap.ValidU64Slice ptr remaining :=
        heap.validU64Slice_mono ptr (remaining + 1) remaining hvalid (by omega)
      rcases slicePolyRep_exists_unique heap ptr remaining this._p.toNat
          hvalidPrefix with ⟨prefixPoly, hprefix, _⟩
      rcases heap.readU64_of_valid ptr (remaining + 1) remaining hvalid
          (by omega) with ⟨value, hread⟩
      have hpolyStep := slicePolyRep_succ_eq_add_monomial heap ptr remaining
        this._p.toNat value prefixPoly poly hprefix hrep hread
      by_cases hzero : value = 0
      · subst value
        rcases ih output prefixPoly hvalidPrefix hprefix with
          ⟨result, hrun, hsem⟩
        refine ⟨result, ?_, ?_⟩
        · simp [denseToSparseLoop, hread, hrun]
        · rw [hsem, hpolyStep]
          simp
      · let output' := output.push (⟨remaining⟩, ⟨value, this._p⟩)
        rcases ih output' prefixPoly hvalidPrefix hprefix with
          ⟨result, hrun, hsem⟩
        refine ⟨result, ?_, ?_⟩
        · simp [denseToSparseLoop, hread, hzero, output', hrun]
        · rw [hsem, sparseToPoly_push_raw, hpolyStep]
          ac_rfl

/-- The exact dense-to-sparse execution denotes the same polynomial as its
normalized raw input buffer. -/
theorem dense_upoly_zp_to_upoly_raw_ir_refines (this : DenseUPolyZp)
    (heap : RawHeap) (ptr : RawPtr UInt64) (length : Nat)
    (poly : Polynomial (ZMod this._p.toNat))
    (hrep : RawDensePolyRep this heap ptr length poly) :
    ∃ result,
      dense_upoly_zp_to_upoly_raw_ir this heap ptr length = .ok result ∧
      SparsePolyZp.toPoly this._p.toNat result = poly := by
  rcases denseToSparseLoop_refines this heap ptr length #[] poly hrep.1
      hrep.2.2.1 with ⟨result, hrun, hsem⟩
  exact ⟨result, hrun, by simpa using hsem⟩

/-- Appending the next lower nonzero raw coefficient preserves the canonical
descending sparse representation. -/
theorem canonical_push_lower (p : UInt64) (output : SparsePolyZp)
    (degree : Nat) (value : UInt64)
    (hcanonical : SparsePolyZp.Canonical p.toNat output)
    (hdegrees : ∀ term ∈ output.toList, degree < term.1.deg)
    (hvalue : value.toNat < p.toNat) (hnonzero : value ≠ 0) :
    SparsePolyZp.Canonical p.toNat
      (output.push (⟨degree⟩, ⟨value, p⟩)) := by
  rcases hcanonical with ⟨hreduced, hchain, hnonzeroOutput⟩
  constructor
  · intro term hterm
    simp only [Array.toList_push, List.mem_append, List.mem_singleton] at hterm
    rcases hterm with hterm | rfl
    · exact hreduced term hterm
    · exact ⟨rfl, hvalue⟩
  constructor
  · rw [Array.toList_push, List.isChain_append]
    refine ⟨hchain, by simp, ?_⟩
    intro left hleft right hright
    have hrightEq : right = (⟨degree⟩, ⟨value, p⟩) := by
      symm
      simpa using hright
    subst right
    exact hdegrees left (List.mem_of_mem_getLast? hleft)
  · intro term hterm
    simp only [Array.toList_push, List.mem_append, List.mem_singleton] at hterm
    rcases hterm with hterm | rfl
    · exact hnonzeroOutput term hterm
    · exact hnonzero

/-- Shape invariant of the reverse scan.  Existing output terms are exactly
the already-scanned higher coefficients, hence all lie above `remaining`. -/
theorem denseToSparseLoop_canonical (this : DenseUPolyZp) (heap : RawHeap)
    (ptr : RawPtr UInt64) (length remaining : Nat) (output : SparsePolyZp)
    (hvalid : heap.ValidU64Slice ptr length) (hremaining : remaining ≤ length)
    (hrawCanonical : CanonicalU64Prefix heap ptr length this._p)
    (houtput : SparsePolyZp.Canonical this._p.toNat output)
    (hdegrees : ∀ term ∈ output.toList, remaining ≤ term.1.deg) :
    ∃ result,
      denseToSparseLoop this heap ptr remaining output = .ok result ∧
      SparsePolyZp.Canonical this._p.toNat result := by
  induction remaining generalizing output with
  | zero => exact ⟨output, rfl, houtput⟩
  | succ remaining ih =>
      rcases heap.readU64_of_valid ptr length remaining hvalid (by omega) with
        ⟨value, hread⟩
      have hvalue : value.toNat < this._p.toNat :=
        hrawCanonical remaining value (by omega) hread
      by_cases hzero : value = 0
      · rcases ih output (by omega) houtput
          (fun term hterm => by
            exact Nat.le_trans (Nat.le_succ remaining) (hdegrees term hterm)) with
          ⟨result, hrun, hcanonical⟩
        exact ⟨result, by
          simp [denseToSparseLoop, hread, hzero, hrun], hcanonical⟩
      · let output' := output.push (⟨remaining⟩, ⟨value, this._p⟩)
        have houtput' : SparsePolyZp.Canonical this._p.toNat output' :=
          canonical_push_lower this._p output remaining value houtput
            (fun term hterm => Nat.lt_of_succ_le (hdegrees term hterm))
            hvalue hzero
        have hdegrees' : ∀ term ∈ output'.toList,
            remaining ≤ term.1.deg := by
          intro term hterm
          simp only [output', Array.toList_push, List.mem_append,
            List.mem_singleton] at hterm
          rcases hterm with hterm | rfl
          · exact Nat.le_trans (Nat.le_succ remaining) (hdegrees term hterm)
          · exact le_rfl
        rcases ih output' (by omega) houtput' hdegrees' with
          ⟨result, hrun, hcanonical⟩
        exact ⟨result, by
          simp [denseToSparseLoop, hread, hzero, output', hrun], hcanonical⟩

/-- The exact dense-to-sparse scan produces canonical descending sparse
terms when the raw dense coefficients are canonical residues. -/
theorem dense_upoly_zp_to_upoly_raw_ir_canonical (this : DenseUPolyZp)
    (heap : RawHeap) (ptr : RawPtr UInt64) (length : Nat)
    (hvalid : heap.ValidU64Slice ptr length)
    (hcanonical : CanonicalU64Prefix heap ptr length this._p) :
    ∃ result,
      dense_upoly_zp_to_upoly_raw_ir this heap ptr length = .ok result ∧
      SparsePolyZp.Canonical this._p.toNat result := by
  apply denseToSparseLoop_canonical this heap ptr length length #[] hvalid
    le_rfl hcanonical
  · exact ⟨by
      intro term hterm
      simp at hterm,
      by simp,
      by simp⟩
  · intro term hterm
    simp at hterm

/-- A canonical sparse C++ polynomial and the normalized raw dense buffer
created from it denote exactly the same mathematical polynomial.  This is a
representation relation only; it performs no GCD computation. -/
structure SparseRawDenseRep (this : DenseUPolyZp)
    (heap : RawHeap) (ptr : RawPtr UInt64) (length : Nat)
    (sparse : SparsePolyZp) : Prop where
  canonical : SparsePolyZp.Canonical this._p.toNat sparse
  dense : RawDensePolyRep this heap ptr length
    (SparsePolyZp.toPoly this._p.toNat sparse)

/-- Exact postcondition of the C++ dense-to-sparse `to_upoly` conversion. -/
structure RawDenseSparseResult (this : DenseUPolyZp)
    (heap : RawHeap) (ptr : RawPtr UInt64) (length : Nat)
    (sparse : SparsePolyZp) : Prop where
  canonical : SparsePolyZp.Canonical this._p.toNat sparse
  dense : RawDensePolyRep this heap ptr length
    (SparsePolyZp.toPoly this._p.toNat sparse)

theorem SparseRawDenseRep.raw (this : DenseUPolyZp)
    (heap : RawHeap) (ptr : RawPtr UInt64) (length : Nat)
    (sparse : SparsePolyZp)
    (hrep : SparseRawDenseRep this heap ptr length sparse) :
    RawDensePolyRep this heap ptr length
      (SparsePolyZp.toPoly this._p.toNat sparse) :=
  hrep.dense

theorem RawDenseSparseResult.toPoly_unique (this : DenseUPolyZp)
    (heap : RawHeap) (ptr : RawPtr UInt64) (length : Nat)
    (left right : SparsePolyZp)
    (hleft : RawDenseSparseResult this heap ptr length left)
    (hright : RawDenseSparseResult this heap ptr length right) :
    left = right := by
  apply SparsePolyZp.toPoly_inj_canonical this._p.toNat left right
      hleft.canonical hright.canonical
  rcases slicePolyRep_exists_unique heap ptr length this._p.toNat
      hleft.dense.1 with ⟨poly, _, hunique⟩
  exact (hunique _ hleft.dense.2.2.1).trans
    (hunique _ hright.dense.2.2.1).symm

/-- The actual `to_upoly` result satisfies the complete canonical
dense-to-sparse output relation. -/
theorem dense_upoly_zp_to_upoly_raw_ir_result (this : DenseUPolyZp)
    (heap : RawHeap) (ptr : RawPtr UInt64) (length : Nat)
    (poly : Polynomial (ZMod this._p.toNat))
    (hrep : RawDensePolyRep this heap ptr length poly) :
    ∃ result,
      dense_upoly_zp_to_upoly_raw_ir this heap ptr length = .ok result ∧
      RawDenseSparseResult this heap ptr length result := by
  rcases dense_upoly_zp_to_upoly_raw_ir_refines this heap ptr length poly hrep
      with ⟨semanticResult, hsemanticRun, hsemantic⟩
  rcases dense_upoly_zp_to_upoly_raw_ir_canonical this heap ptr length hrep.1
      hrep.2.1 with ⟨canonicalResult, hcanonicalRun, hcanonical⟩
  have heq : canonicalResult = semanticResult :=
    Except.ok.inj (hcanonicalRun.symm.trans hsemanticRun)
  subst canonicalResult
  exact ⟨semanticResult, hsemanticRun, hcanonical, by
    simpa [hsemantic] using hrep⟩

end CLPoly.Impl.StrictPolynomialGCDRefinement
