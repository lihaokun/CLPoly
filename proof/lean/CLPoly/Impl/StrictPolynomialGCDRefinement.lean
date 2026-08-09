import CLPoly.Impl.StrictGCDHGCDRefinement
import CLPoly.Math.Univariate

set_option autoImplicit false

namespace CLPoly.Impl.StrictPolynomialGCDRefinement

open CLPoly.Impl.StrictEuclidRefinement
open CLPoly.Impl.StrictDivremRefinement
open CLPoly.Impl.RawPolynomialRep
open CLPoly.Math
open Generated.StrictMul

/-- One concrete raw coefficient write into a zero coefficient adds exactly
the corresponding L2 monomial. -/
theorem slicePolyRep_write_add_monomial (heap heap' : RawHeap)
    (ptr : RawPtr UInt64) (length p degree : Nat) (value : UInt64)
    (poly : Polynomial (ZMod p))
    (hvalid : heap.ValidU64Slice ptr length)
    (hdegree : degree < length)
    (hrep : SlicePolyRep heap ptr length p poly)
    (hzero : poly.coeff degree = 0)
    (hwrite : heap.writeU64 ptr degree value = .ok heap') :
    SlicePolyRep heap' ptr length p
      (poly + Polynomial.monomial degree (value.toNat : ZMod p)) := by
  have hvalid' : heap'.ValidU64Slice ptr length :=
    (RawHeap.writeU64_preserves_valid heap heap' ptr degree value hwrite
      ptr length).mp hvalid
  rcases slicePolyRep_exists_unique heap' ptr length p hvalid' with
    ⟨resultPoly, hresultRep, hunique⟩
  have heq : poly + Polynomial.monomial degree
      (value.toNat : ZMod p) = resultPoly := by
    ext observedDegree
    by_cases hobserved : observedDegree < length
    · rcases slicePolyRep_coeff heap ptr length p poly hrep observedDegree
          hobserved with ⟨old, hreadOld, hcoeffOld⟩
      rcases slicePolyRep_coeff heap' ptr length p resultPoly hresultRep
          observedDegree hobserved with ⟨observed, hreadObserved, hcoeffObserved⟩
      by_cases heqDegree : observedDegree = degree
      · subst observedDegree
        have hreadNow := RawHeap.readU64_writeU64_same heap heap' ptr degree
          value hwrite
        have hvalue : observed = value :=
          Except.ok.inj (hreadObserved.symm.trans hreadNow)
        subst observed
        rw [Polynomial.coeff_add, hzero, Polynomial.coeff_monomial]
        simp [hcoeffObserved]
      · have hreadStill := RawHeap.readU64_writeU64_ne heap heap' ptr ptr
          degree observedDegree value old hwrite hreadOld
          (Or.inr (by omega))
        have hvalue : observed = old :=
          Except.ok.inj (hreadObserved.symm.trans hreadStill)
        subst observed
        rw [Polynomial.coeff_add, hcoeffOld, hcoeffObserved,
          Polynomial.coeff_monomial]
        simp [show degree ≠ observedDegree by omega]
    · rw [Polynomial.coeff_add,
        slicePolyRep_coeff_zero_of_length_le heap ptr length p poly hrep
          observedDegree (by omega),
        slicePolyRep_coeff_zero_of_length_le heap' ptr length p resultPoly
          hresultRep observedDegree (by omega)]
      simp [Polynomial.coeff_monomial,
        show degree ≠ observedDegree by omega]
  rw [heq]
  exact hresultRep

theorem canonicalU64Prefix_write (heap heap' : RawHeap)
    (ptr : RawPtr UInt64) (length degree : Nat) (value modulus : UInt64)
    (hvalid : heap.ValidU64Slice ptr length)
    (hcanonical : CanonicalU64Prefix heap ptr length modulus)
    (hvalue : value.toNat < modulus.toNat)
    (hwrite : heap.writeU64 ptr degree value = .ok heap') :
    CanonicalU64Prefix heap' ptr length modulus := by
  intro index observed hindex hreadObserved
  by_cases heq : index = degree
  · subst index
    have hreadNow := RawHeap.readU64_writeU64_same heap heap' ptr degree
      value hwrite
    have hobserved : observed = value :=
      Except.ok.inj (hreadObserved.symm.trans hreadNow)
    simpa [hobserved] using hvalue
  · rcases heap.readU64_of_valid ptr length index hvalid hindex with
      ⟨old, hreadOld⟩
    have hreadStill := RawHeap.readU64_writeU64_ne heap heap' ptr ptr degree
      index value old hwrite hreadOld (Or.inr (by omega))
    have hobserved : observed = old :=
      Except.ok.inj (hreadObserved.symm.trans hreadStill)
    subst observed
    exact hcanonical index old hindex hreadOld

/-- L2 polynomial denoted by the sparse terms already consumed by the source
forward iterator. -/
noncomputable def sparsePrefixPoly (p : Nat) (sparse : SparsePolyZp)
    (index : Nat) : Polynomial (ZMod p) :=
  SparsePolyZp.toPoly p (sparse.extract 0 index)

/-- Extending the consumed sparse prefix by its next actual array element
adds exactly that element's monomial. -/
theorem sparsePrefixPoly_succ (p : Nat) (sparse : SparsePolyZp)
    (index : Nat) (hindex : index < sparse.size) :
    sparsePrefixPoly p sparse (index + 1) =
      sparsePrefixPoly p sparse index +
        Polynomial.monomial sparse[index].1.deg
          (Zp.toZMod p sparse[index].2) := by
  simp [sparsePrefixPoly, SparsePolyZp.toPoly, Array.toList_extract,
    List.extract, List.take_add_one, hindex, listSum_append, listSum]

theorem sparsePrefixPoly_size (p : Nat) (sparse : SparsePolyZp) :
    sparsePrefixPoly p sparse sparse.size = SparsePolyZp.toPoly p sparse := by
  simp [sparsePrefixPoly]

theorem listSum_coeff_zero_of_all_gt (p degree : Nat)
    (terms : List (UMonomial × Zp))
    (hgreater : ∀ term ∈ terms, degree < term.1.deg) :
    (listSum p terms).coeff degree = 0 := by
  induction terms with
  | nil => simp [listSum]
  | cons term rest ih =>
      have hterm := hgreater term List.mem_cons_self
      have hrest : ∀ item ∈ rest, degree < item.1.deg := by
        intro item hitem
        exact hgreater item (List.mem_cons_of_mem term hitem)
      rw [listSum_cons, Polynomial.coeff_add, Polynomial.coeff_monomial,
        ih hrest]
      simp [show term.1.deg ≠ degree by omega]

/-- Strict descending canonical degrees establish the fresh-cell premise used
by the raw sparse iterator. -/
theorem sparsePrefixPoly_coeff_current_eq_zero (p : Nat)
    (sparse : SparsePolyZp) (hcanonical : SparsePolyZp.Canonical p sparse)
    (index : Nat) (hindex : index < sparse.size) :
    (sparsePrefixPoly p sparse index).coeff sparse[index].1.deg = 0 := by
  unfold sparsePrefixPoly SparsePolyZp.toPoly
  apply listSum_coeff_zero_of_all_gt
  intro term hterm
  have htermTake : term ∈ sparse.toList.take index := by
    simpa [Array.toList_extract, List.extract] using hterm
  have hcurrentDrop : sparse[index] ∈ sparse.toList.drop index := by
    have hlistIndex : sparse.toList[index] = sparse[index] :=
      Array.getElem_toList hindex
    rw [← hlistIndex]
    apply List.mem_drop_iff_getElem.mpr
    exact ⟨0, by simpa, rfl⟩
  have hpairwise : List.Pairwise
      (fun a b : UMonomial × Zp => a.1.deg > b.1.deg) sparse.toList :=
    List.isChain_iff_pairwise.mp hcanonical.2.1
  exact hpairwise.rel_of_mem_take_of_mem_drop htermTake hcurrentDrop

/-- Exact raw lowering of the sparse-term write loop in the
`dense_upoly_zp(upolynomial,Zp)` constructor.  The preceding vector resize is
represented by `mulZeroPadLoop` below, which performs the same forward zero
writes on the allocated coefficient buffer. -/
def sparseToDenseWriteLoop (ptr : RawPtr UInt64) (length : Nat)
    (sparse : SparsePolyZp) (index : Nat) (heap : RawHeap) : RawExec RawHeap :=
  if h : index < sparse.size then
    let term := sparse[index]
    match heap.writeU64 ptr term.1.deg term.2.val with
    | .error fault => .error fault
    | .ok heap' =>
        sparseToDenseWriteLoop ptr length sparse (index + 1) heap'
  else
    .ok heap
termination_by sparse.size - index
decreasing_by omega

/-- Exact raw constructor path: zero-initialize the resized dense vector,
then write every stored sparse coefficient at its source degree. -/
def sparse_upoly_zp_to_dense_raw_ir (ptr : RawPtr UInt64) (length : Nat)
    (sparse : SparsePolyZp) (heap : RawHeap) : RawExec RawHeap :=
  match mulZeroPadLoop ptr 0 length 0 heap with
  | .error fault => .error fault
  | .ok heap' => sparseToDenseWriteLoop ptr length sparse 0 heap'

/-- The actual sparse-term loop cannot fault when the target allocation is
valid and every source degree lies inside the constructor length. -/
theorem sparseToDenseWriteLoop_succeeds (ptr : RawPtr UInt64) (length : Nat)
    (sparse : SparsePolyZp) (index : Nat) (heap : RawHeap)
    (hvalid : heap.ValidU64Slice ptr length)
    (hdegree : ∀ i (hi : i < sparse.size), sparse[i].1.deg < length) :
    ∃ heap', sparseToDenseWriteLoop ptr length sparse index heap = .ok heap' := by
  rw [sparseToDenseWriteLoop]
  split
  next hmore =>
    let term := sparse[index]
    have htermDegree : term.1.deg < length := hdegree index hmore
    rcases heap.writeU64_of_valid ptr length term.1.deg term.2.val hvalid
        htermDegree with ⟨heap1, hwrite⟩
    simp only [term, hwrite]
    have hvalid1 : heap1.ValidU64Slice ptr length :=
      (RawHeap.writeU64_preserves_valid heap heap1 ptr term.1.deg term.2.val
        hwrite ptr length).mp hvalid
    exact sparseToDenseWriteLoop_succeeds ptr length sparse (index + 1)
      heap1 hvalid1 hdegree
  next hdone => exact ⟨heap, rfl⟩
termination_by sparse.size - index
decreasing_by omega

/-- Sparse coefficient writes change values only; all raw allocations retain
their exact layout throughout the generated constructor loop. -/
theorem sparseToDenseWriteLoop_sameLayout (ptr : RawPtr UInt64)
    (length : Nat) (sparse : SparsePolyZp) (index : Nat)
    (heap heap' : RawHeap) (hvalid : heap.ValidU64Slice ptr length)
    (hdegree : ∀ i (hi : i < sparse.size), sparse[i].1.deg < length)
    (hrun : sparseToDenseWriteLoop ptr length sparse index heap = .ok heap') :
    RawHeap.SameLayout heap heap' := by
  rw [sparseToDenseWriteLoop] at hrun
  split at hrun
  next hmore =>
    let term := sparse[index]
    have htermDegree : term.1.deg < length := hdegree index hmore
    rcases heap.writeU64_of_valid ptr length term.1.deg term.2.val hvalid
        htermDegree with ⟨heap1, hwrite⟩
    simp only [term, hwrite] at hrun
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 ptr
      term.1.deg term.2.val hwrite
    have hvalid1 := (hlayout1 ptr length).mp hvalid
    have hlayout2 := sparseToDenseWriteLoop_sameLayout ptr length sparse
      (index + 1) heap1 heap' hvalid1 hdegree hrun
    exact fun other count =>
      (hlayout1 other count).trans (hlayout2 other count)
  next hdone =>
    have heq : heap' = heap := Except.ok.inj hrun.symm
    subst heap'
    exact fun _ _ => Iff.rfl
termination_by sparse.size - index
decreasing_by omega

/-- Semantic invariant of the actual sparse write loop.  The explicit
fresh-degree premise records that the canonical source iterator never writes
the same monomial degree twice. -/
theorem sparseToDenseWriteLoop_refines (ptr : RawPtr UInt64) (length p : Nat)
    (sparse : SparsePolyZp) (index : Nat) (heap : RawHeap)
    (hindex : index ≤ sparse.size)
    (hvalid : heap.ValidU64Slice ptr length)
    (hdegree : ∀ i (hi : i < sparse.size), sparse[i].1.deg < length)
    (hfresh : ∀ i (hi : i < sparse.size),
      (sparsePrefixPoly p sparse i).coeff sparse[i].1.deg = 0)
    (hrep : SlicePolyRep heap ptr length p
      (sparsePrefixPoly p sparse index)) :
    ∃ heap', sparseToDenseWriteLoop ptr length sparse index heap = .ok heap' ∧
      RawHeap.SameLayout heap heap' ∧
      SlicePolyRep heap' ptr length p (SparsePolyZp.toPoly p sparse) := by
  rw [sparseToDenseWriteLoop]
  split
  next hmore =>
    let term := sparse[index]
    have htermDegree : term.1.deg < length := hdegree index hmore
    rcases heap.writeU64_of_valid ptr length term.1.deg term.2.val hvalid
        htermDegree with ⟨heap1, hwrite⟩
    simp only [term, hwrite]
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 ptr
      term.1.deg term.2.val hwrite
    have hvalid1 := (hlayout1 ptr length).mp hvalid
    have hrep1 : SlicePolyRep heap1 ptr length p
        (sparsePrefixPoly p sparse (index + 1)) := by
      rw [sparsePrefixPoly_succ p sparse index hmore]
      simpa [term, Zp.toZMod] using
        slicePolyRep_write_add_monomial heap heap1 ptr length p
          term.1.deg term.2.val (sparsePrefixPoly p sparse index) hvalid
          htermDegree hrep (hfresh index hmore) hwrite
    rcases sparseToDenseWriteLoop_refines ptr length p sparse (index + 1)
        heap1 (by omega) hvalid1 hdegree hfresh hrep1 with
      ⟨heap2, hrun, hlayout2, hresult⟩
    exact ⟨heap2, hrun, fun other count =>
      (hlayout1 other count).trans (hlayout2 other count), hresult⟩
  next hdone =>
    have heq : index = sparse.size := by omega
    subst index
    exact ⟨heap, rfl, fun _ _ => Iff.rfl, by
      simpa [sparsePrefixPoly_size] using hrep⟩
termination_by sparse.size - index
decreasing_by omega

theorem sparseToDenseWriteLoop_canonical (ptr : RawPtr UInt64)
    (length : Nat) (modulus : UInt64) (sparse : SparsePolyZp)
    (index : Nat) (heap heap' : RawHeap)
    (hvalid : heap.ValidU64Slice ptr length)
    (hdegree : ∀ i (hi : i < sparse.size), sparse[i].1.deg < length)
    (hreduced : ∀ i (hi : i < sparse.size),
      sparse[i].2.val.toNat < modulus.toNat)
    (hcanonical : CanonicalU64Prefix heap ptr length modulus)
    (hrun : sparseToDenseWriteLoop ptr length sparse index heap = .ok heap') :
    CanonicalU64Prefix heap' ptr length modulus := by
  rw [sparseToDenseWriteLoop] at hrun
  split at hrun
  next hmore =>
    let term := sparse[index]
    have htermDegree : term.1.deg < length := hdegree index hmore
    rcases heap.writeU64_of_valid ptr length term.1.deg term.2.val hvalid
        htermDegree with ⟨heap1, hwrite⟩
    simp only [term, hwrite] at hrun
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 ptr
      term.1.deg term.2.val hwrite
    have hvalid1 := (hlayout1 ptr length).mp hvalid
    have hcanonical1 := canonicalU64Prefix_write heap heap1 ptr length
      term.1.deg term.2.val modulus hvalid hcanonical
      (hreduced index hmore) hwrite
    exact sparseToDenseWriteLoop_canonical ptr length modulus sparse
      (index + 1) heap1 heap' hvalid1 hdegree hreduced hcanonical1 hrun
  next hdone =>
    have heq : heap' = heap := Except.ok.inj hrun.symm
    simpa [heq] using hcanonical
termination_by sparse.size - index
decreasing_by omega

/-- Both physical phases of the sparse-to-dense constructor terminate without
a raw access fault under their exact allocation and degree bounds. -/
theorem sparse_upoly_zp_to_dense_raw_ir_succeeds
    (ptr : RawPtr UInt64) (length : Nat) (sparse : SparsePolyZp)
    (heap : RawHeap) (hvalid : heap.ValidU64Slice ptr length)
    (hdegree : ∀ i (hi : i < sparse.size), sparse[i].1.deg < length) :
    ∃ heap', sparse_upoly_zp_to_dense_raw_ir ptr length sparse heap =
      .ok heap' := by
  rcases CLPoly.Impl.StrictMulRefinement.mulZeroPadLoop_refines ptr 0 length 0
      1 heap 0 1 (by omega) (by decide) (by simpa using hvalid)
      (CLPoly.Impl.StrictHGCDRawRefinement.slicePolyRep_zero_length_any
        heap ptr 1)
      (by
        intro i value hi
        omega) with ⟨heap1, hzero, hlayout, _, _⟩
  have hvalid1 := (hlayout ptr length).mp hvalid
  rcases sparseToDenseWriteLoop_succeeds ptr length sparse 0 heap1 hvalid1
      hdegree with ⟨heap2, hwrites⟩
  exact ⟨heap2, by
    simp [sparse_upoly_zp_to_dense_raw_ir, hzero, hwrites]⟩

/-- The concrete constructor result keeps the target dense allocation valid;
this is the raw→safe boundary needed by the following representation proof. -/
theorem sparse_upoly_zp_to_dense_raw_ir_valid
    (ptr : RawPtr UInt64) (length : Nat) (sparse : SparsePolyZp)
    (heap heap' : RawHeap) (hvalid : heap.ValidU64Slice ptr length)
    (hdegree : ∀ i (hi : i < sparse.size), sparse[i].1.deg < length)
    (hrun : sparse_upoly_zp_to_dense_raw_ir ptr length sparse heap =
      .ok heap') :
    heap'.ValidU64Slice ptr length := by
  rcases CLPoly.Impl.StrictMulRefinement.mulZeroPadLoop_refines ptr 0 length 0
      1 heap 0 1 (by omega) (by decide) (by simpa using hvalid)
      (CLPoly.Impl.StrictHGCDRawRefinement.slicePolyRep_zero_length_any
        heap ptr 1)
      (by
        intro i value hi
        omega) with ⟨heap1, hzero, hlayoutZero, _, _⟩
  have hvalid1 := (hlayoutZero ptr length).mp hvalid
  have hwrites : sparseToDenseWriteLoop ptr length sparse 0 heap1 =
      .ok heap' := by
    simpa [sparse_upoly_zp_to_dense_raw_ir, hzero] using hrun
  have hlayoutWrites := sparseToDenseWriteLoop_sameLayout ptr length sparse
    0 heap1 heap' hvalid1 hdegree hwrites
  exact (hlayoutWrites ptr length).mp hvalid1

/-- End-to-end semantic result of the concrete sparse-to-dense writes. -/
theorem sparse_upoly_zp_to_dense_raw_ir_refines
    (ptr : RawPtr UInt64) (length p : Nat) (sparse : SparsePolyZp)
    (heap : RawHeap) (hvalid : heap.ValidU64Slice ptr length)
    (hcanonical : SparsePolyZp.Canonical p sparse)
    (hdegree : ∀ i (hi : i < sparse.size), sparse[i].1.deg < length) :
    ∃ heap', sparse_upoly_zp_to_dense_raw_ir ptr length sparse heap =
        .ok heap' ∧
      RawHeap.SameLayout heap heap' ∧
      SlicePolyRep heap' ptr length p (SparsePolyZp.toPoly p sparse) := by
  rcases CLPoly.Impl.StrictMulRefinement.mulZeroPadLoop_refines ptr 0 length 0
      p heap 0 1 (by omega) (by decide) (by simpa using hvalid)
      (CLPoly.Impl.StrictHGCDRawRefinement.slicePolyRep_zero_length_any
        heap ptr p)
      (by
        intro i value hi
        omega) with ⟨heap1, hzero, hlayoutZero, hzeroRep, _⟩
  have hvalid1 := (hlayoutZero ptr length).mp hvalid
  have hprefixZero : sparsePrefixPoly p sparse 0 = 0 := by
    simp [sparsePrefixPoly, SparsePolyZp.toPoly]
  have hprefixRep : SlicePolyRep heap1 ptr length p
      (sparsePrefixPoly p sparse 0) := by
    simpa [hprefixZero] using hzeroRep
  rcases sparseToDenseWriteLoop_refines ptr length p sparse 0 heap1
      (Nat.zero_le _) hvalid1 hdegree
      (sparsePrefixPoly_coeff_current_eq_zero p sparse hcanonical)
      hprefixRep with ⟨heap2, hwrites, hlayoutWrites, hresult⟩
  exact ⟨heap2, by
    simp [sparse_upoly_zp_to_dense_raw_ir, hzero, hwrites],
    fun other count =>
      (hlayoutZero other count).trans (hlayoutWrites other count), hresult⟩

/-- Every raw value produced by the constructor is either its explicit zero
initialization or a reduced residue copied from a canonical sparse term. -/
theorem sparse_upoly_zp_to_dense_raw_ir_canonical
    (this : DenseUPolyZp) (ptr : RawPtr UInt64) (length : Nat)
    (sparse : SparsePolyZp) (heap heap' : RawHeap)
    (hmodulus : this._p ≠ 0)
    (hvalid : heap.ValidU64Slice ptr length)
    (hcanonical : SparsePolyZp.Canonical this._p.toNat sparse)
    (hdegree : ∀ i (hi : i < sparse.size), sparse[i].1.deg < length)
    (hrun : sparse_upoly_zp_to_dense_raw_ir ptr length sparse heap =
      .ok heap') :
    CanonicalU64Prefix heap' ptr length this._p := by
  rcases CLPoly.Impl.StrictMulRefinement.mulZeroPadLoop_refines ptr 0 length 0
      this._p.toNat heap 0 this._p (by omega) hmodulus
      (by simpa using hvalid)
      (CLPoly.Impl.StrictHGCDRawRefinement.slicePolyRep_zero_length_any
        heap ptr this._p.toNat)
      (by
        intro i value hi
        omega) with ⟨heap1, hzero, hlayoutZero, _, hzeroCanonical⟩
  have hvalid1 := (hlayoutZero ptr length).mp hvalid
  have hwrites : sparseToDenseWriteLoop ptr length sparse 0 heap1 =
      .ok heap' := by
    simpa [sparse_upoly_zp_to_dense_raw_ir, hzero] using hrun
  apply sparseToDenseWriteLoop_canonical ptr length this._p sparse 0 heap1
    heap' hvalid1 hdegree
  · intro i hi
    have hmember : sparse[i] ∈ sparse.toList := by simp
    exact (hcanonical.1 sparse[i] hmember).2
  · simpa using hzeroCanonical
  · exact hwrites

/-- Before the constructor's final normalization check, the returned physical
slice already has the exact polynomial and canonical residue semantics. -/
theorem sparse_upoly_zp_to_dense_raw_ir_canonicalSlice
    (this : DenseUPolyZp) (ptr : RawPtr UInt64) (length : Nat)
    (sparse : SparsePolyZp) (heap : RawHeap)
    (hmodulus : this._p ≠ 0)
    (hvalid : heap.ValidU64Slice ptr length)
    (hcanonical : SparsePolyZp.Canonical this._p.toNat sparse)
    (hdegree : ∀ i (hi : i < sparse.size), sparse[i].1.deg < length) :
    ∃ heap', sparse_upoly_zp_to_dense_raw_ir ptr length sparse heap =
        .ok heap' ∧
      RawHeap.SameLayout heap heap' ∧
      RawCanonicalPolySlice this heap' ptr length
        (SparsePolyZp.toPoly this._p.toNat sparse) := by
  rcases sparse_upoly_zp_to_dense_raw_ir_refines ptr length this._p.toNat
      sparse heap hvalid hcanonical hdegree with
    ⟨heap', hrun, hlayout, hslice⟩
  exact ⟨heap', hrun, hlayout,
    sparse_upoly_zp_to_dense_raw_ir_valid ptr length sparse heap heap'
      hvalid hdegree hrun,
    sparse_upoly_zp_to_dense_raw_ir_canonical this ptr length sparse heap
      heap' hmodulus hvalid hcanonical hdegree hrun,
    hslice⟩

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
