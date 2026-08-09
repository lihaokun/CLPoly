import CLPoly.Impl.StrictGCDHGCDRefinement
import CLPoly.Math.Univariate

set_option autoImplicit false

namespace CLPoly.Impl.StrictPolynomialGCDRefinement

open CLPoly.Impl.StrictEuclidRefinement
open CLPoly.Impl.StrictDivremRefinement
open CLPoly.Impl.StrictGCDHGCDRefinement
open CLPoly.Impl.StrictWordArithmetic
open CLPoly.Impl.RawPolynomialRep
open CLPoly.Math
open Generated.StrictMul

/-- Observable result shared by the three object-level dense GCD branches. -/
structure DenseGcdRawResult where
  heap : RawHeap
  lenG : Nat

def denseGcdCutoff : Nat := 340

/-- Exact object-level dense GCD dispatch after the source degree-order swap
has selected `large` and `small`. -/
def dense_upoly_zp_gcd_ordered_raw_ir (this : DenseUPolyZp)
    (M : HgcdMat) (hM : M.Valid)
    (G large small aBuf bBuf J Q R : RawPtr UInt64) (W3 : RawPtr Word3)
    (W scratch : RawPtr UInt64)
    (euclidQ euclidR : RawPtr UInt64) (euclidW3 : RawPtr Word3)
    (lenLarge lenSmall : Nat)
    (loopDecrease : Generated.StrictGCDHGCD.HgcdGcdLoopLengthDecreases
      this M hM W scratch)
    (heap : RawHeap) : RawExec DenseGcdRawResult :=
  if lenSmall = 0 then
    match heap.copyU64 G large lenLarge with
    | .error fault => .error fault
    | .ok heap' => .ok ⟨heap', lenLarge⟩
  else if lenSmall < denseGcdCutoff then
    match CLPoly.Impl.StrictEuclidRefinement.strictGcdEuclidRaw this G large
        small aBuf bBuf Q R W3 lenLarge lenSmall heap with
    | .error fault => .error fault
    | .ok result => .ok ⟨result.heap, result.lenG⟩
  else
    match Generated.StrictGCDHGCD.dense_upoly_zp__gcd_hgcd_raw_ir this M hM
        G large small J Q R W3 W scratch aBuf bBuf euclidQ euclidR euclidW3
        lenLarge lenSmall
        (CLPoly.Impl.StrictDivremRefinement.euclidDivremLengthDecreases this
          euclidQ euclidW3)
        loopDecrease heap with
    | .error fault => .error fault
    | .ok result => .ok ⟨result.heap, result.lenG⟩

/-- Exact source swap followed by the ordered dense GCD dispatcher. -/
def dense_upoly_zp_gcd_raw_ir (this : DenseUPolyZp)
    (M : HgcdMat) (hM : M.Valid)
    (G A B aBuf bBuf J Q R : RawPtr UInt64) (W3 : RawPtr Word3)
    (W scratch : RawPtr UInt64)
    (euclidQ euclidR : RawPtr UInt64) (euclidW3 : RawPtr Word3)
    (lenA lenB : Nat)
    (loopDecrease : Generated.StrictGCDHGCD.HgcdGcdLoopLengthDecreases
      this M hM W scratch)
    (heap : RawHeap) : RawExec DenseGcdRawResult :=
  if lenA < lenB then
    dense_upoly_zp_gcd_ordered_raw_ir this M hM G B A aBuf bBuf J Q R W3 W
      scratch euclidQ euclidR euclidW3 lenB lenA loopDecrease heap
  else
    dense_upoly_zp_gcd_ordered_raw_ir this M hM G A B aBuf bBuf J Q R W3 W
      scratch euclidQ euclidR euclidW3 lenA lenB loopDecrease heap

/-- Refinement of the exact already-ordered object dispatcher. -/
theorem dense_upoly_zp_gcd_ordered_raw_ir_refines (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (M : HgcdMat) (hM : M.Valid)
    (G large small aBuf bBuf J Q R : RawPtr UInt64) (W3 : RawPtr Word3)
    (W scratch : RawPtr UInt64)
    (euclidQ euclidR : RawPtr UInt64) (euclidW3 : RawPtr Word3)
    (lenLarge lenSmall capacity : Nat)
    (loopDecrease : Generated.StrictGCDHGCD.HgcdGcdLoopLengthDecreases
      this M hM W scratch)
    (initialDivrem : GcdHgcdDynamicDivremWorkspaceProvider this large small
      Q R W3)
    (initialHgcd : GcdHgcdInitialHgcdWorkspaceProvider this hcfg hp M hM G J
      small R W scratch lenSmall)
    (loopDivrem : GcdHgcdDynamicDivremWorkspaceProvider this G J Q R W3)
    (loopEuclid : GcdHgcdDynamicEuclidWorkspaceProvider this G J R aBuf bBuf
      euclidQ euclidR euclidW3)
    (loopHgcd : GcdHgcdDynamicHgcdWorkspaceProvider this hcfg hp M hM G J R W
      scratch)
    (heap : RawHeap) (left right : Polynomial (ZMod this._p.toNat))
    (hleft : RawDensePolyRep this heap large lenLarge left)
    (hright : RawDensePolyRep this heap small lenSmall right)
    (hGlarge : heap.ValidU64Slice G lenLarge)
    (hGsmall : heap.ValidU64Slice G lenSmall)
    (hGlargeRegion : G.region ≠ large.region)
    (hGsmallRegion : G.region ≠ small.region)
    (heuclid : GcdEuclidRawWorkspace heap G large small aBuf bBuf Q R W3
      capacity)
    (hlenLarge : lenLarge ≤ capacity) (hlenSmall : lenSmall ≤ capacity)
    (hcapacity : capacity < limbBase) :
    ∃ out result,
      dense_upoly_zp_gcd_ordered_raw_ir this M hM G large small aBuf bBuf J
          Q R W3 W scratch euclidQ euclidR euclidW3 lenLarge lenSmall
          loopDecrease heap = .ok out ∧
      RawDensePolyRep this out.heap G out.lenG result ∧
      normalize (EuclideanDomain.gcd left right) = normalize result := by
  by_cases hzero : lenSmall = 0
  · have hrightZero : right = 0 :=
      slicePolyRep_zero_length heap small this._p.toNat right (by
        simpa [hzero] using hright.2.2.1)
    subst right
    rcases copyU64_refines_rawDense_of_region_ne this heap G large lenLarge
        left hGlarge hGlargeRegion hleft with
      ⟨heap', hcopy, _, hresult⟩
    exact ⟨⟨heap', lenLarge⟩, left, by
      simp [dense_upoly_zp_gcd_ordered_raw_ir, hzero, hcopy],
      hresult, by simp⟩
  · by_cases hsmall : lenSmall < denseGcdCutoff
    · rcases strictGcdEuclidRaw_refines this G large small aBuf bBuf Q R
          W3 lenLarge lenSmall capacity heap left right heuclid hlenLarge
          hlenSmall hcapacity hleft hright hcfg with
        ⟨heap', lenG, result, hrun, hresult, hgcd, _⟩
      exact ⟨⟨heap', lenG⟩, result, by
        simp [dense_upoly_zp_gcd_ordered_raw_ir, hzero, hsmall, hrun],
        hresult, hgcd⟩
    · rcases dense_upoly_zp__gcd_hgcd_raw_ir_succeeds this hcfg hp M hM
          G large small J Q R W3 W scratch aBuf bBuf euclidQ euclidR
          euclidW3 lenLarge lenSmall
          (CLPoly.Impl.StrictDivremRefinement.euclidDivremLengthDecreases this
            euclidQ euclidW3)
          loopDecrease initialDivrem initialHgcd loopDivrem loopEuclid loopHgcd
          heap left right (Nat.pos_of_ne_zero hzero) hleft hright hGsmall
          hGsmallRegion with ⟨raw, result, hrun, hresult, hgcd⟩
      exact ⟨⟨raw.heap, raw.lenG⟩, result, by
        simp [dense_upoly_zp_gcd_ordered_raw_ir, hzero, hsmall, hrun],
        hresult, hgcd⟩

/-- All physical obligations for one concrete ordering of the object-level
dense GCD call.  No result or polynomial operation is stored here. -/
structure DenseGcdOrderedReady (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (M : HgcdMat) (hM : M.Valid)
    (G large small aBuf bBuf J Q R : RawPtr UInt64) (W3 : RawPtr Word3)
    (W scratch : RawPtr UInt64)
    (euclidQ euclidR : RawPtr UInt64) (euclidW3 : RawPtr Word3)
    (lenLarge lenSmall : Nat)
    (heap : RawHeap) (left right : Polynomial (ZMod this._p.toNat)) where
  capacity : Nat
  initialDivrem : GcdHgcdDynamicDivremWorkspaceProvider this large small Q R W3
  initialHgcd : GcdHgcdInitialHgcdWorkspaceProvider this hcfg hp M hM G J
    small R W scratch lenSmall
  loopDivrem : GcdHgcdDynamicDivremWorkspaceProvider this G J Q R W3
  loopEuclid : GcdHgcdDynamicEuclidWorkspaceProvider this G J R aBuf bBuf
    euclidQ euclidR euclidW3
  loopHgcd : GcdHgcdDynamicHgcdWorkspaceProvider this hcfg hp M hM G J R W
    scratch
  leftRep : RawDensePolyRep this heap large lenLarge left
  rightRep : RawDensePolyRep this heap small lenSmall right
  validGLarge : heap.ValidU64Slice G lenLarge
  validGSmall : heap.ValidU64Slice G lenSmall
  gLargeRegion : G.region ≠ large.region
  gSmallRegion : G.region ≠ small.region
  euclid : GcdEuclidRawWorkspace heap G large small aBuf bBuf Q R W3 capacity
  lenLargeCapacity : lenLarge ≤ capacity
  lenSmallCapacity : lenSmall ≤ capacity
  capacityBound : capacity < limbBase

theorem dense_upoly_zp_gcd_ordered_raw_ir_refines_of_ready
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (M : HgcdMat) (hM : M.Valid)
    (G large small aBuf bBuf J Q R : RawPtr UInt64) (W3 : RawPtr Word3)
    (W scratch : RawPtr UInt64)
    (euclidQ euclidR : RawPtr UInt64) (euclidW3 : RawPtr Word3)
    (lenLarge lenSmall : Nat)
    (loopDecrease : Generated.StrictGCDHGCD.HgcdGcdLoopLengthDecreases
      this M hM W scratch)
    (heap : RawHeap) (left right : Polynomial (ZMod this._p.toNat))
    (ready : DenseGcdOrderedReady this hcfg hp M hM G large small aBuf bBuf J
      Q R W3 W scratch euclidQ euclidR euclidW3 lenLarge lenSmall heap left
      right) :
    ∃ out result,
      dense_upoly_zp_gcd_ordered_raw_ir this M hM G large small aBuf bBuf J
          Q R W3 W scratch euclidQ euclidR euclidW3 lenLarge lenSmall
          loopDecrease heap = .ok out ∧
      RawDensePolyRep this out.heap G out.lenG result ∧
      normalize (EuclideanDomain.gcd left right) = normalize result := by
  exact dense_upoly_zp_gcd_ordered_raw_ir_refines this hcfg hp M hM G large
    small aBuf bBuf J Q R W3 W scratch euclidQ euclidR euclidW3 lenLarge
    lenSmall ready.capacity loopDecrease ready.initialDivrem ready.initialHgcd
    ready.loopDivrem ready.loopEuclid ready.loopHgcd heap left right
    ready.leftRep ready.rightRep ready.validGLarge ready.validGSmall
    ready.gLargeRegion ready.gSmallRegion ready.euclid
    ready.lenLargeCapacity ready.lenSmallCapacity ready.capacityBound

/-- End-to-end refinement of the source degree comparison and swap followed by
the exact zero/Euclid/HGCD object dispatcher. -/
theorem dense_upoly_zp_gcd_raw_ir_refines (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (M : HgcdMat) (hM : M.Valid)
    (G A B aBuf bBuf J Q R : RawPtr UInt64) (W3 : RawPtr Word3)
    (W scratch : RawPtr UInt64)
    (euclidQ euclidR : RawPtr UInt64) (euclidW3 : RawPtr Word3)
    (lenA lenB : Nat)
    (loopDecrease : Generated.StrictGCDHGCD.HgcdGcdLoopLengthDecreases
      this M hM W scratch)
    (heap : RawHeap) (left right : Polynomial (ZMod this._p.toNat))
    (readyAB : DenseGcdOrderedReady this hcfg hp M hM G A B aBuf bBuf J Q R
      W3 W scratch euclidQ euclidR euclidW3 lenA lenB heap left right)
    (readyBA : DenseGcdOrderedReady this hcfg hp M hM G B A aBuf bBuf J Q R
      W3 W scratch euclidQ euclidR euclidW3 lenB lenA heap right left) :
    ∃ out result,
      dense_upoly_zp_gcd_raw_ir this M hM G A B aBuf bBuf J Q R W3 W scratch
          euclidQ euclidR euclidW3 lenA lenB loopDecrease heap = .ok out ∧
      RawDensePolyRep this out.heap G out.lenG result ∧
      normalize (EuclideanDomain.gcd left right) = normalize result := by
  by_cases hswap : lenA < lenB
  · rcases dense_upoly_zp_gcd_ordered_raw_ir_refines_of_ready this hcfg hp
        M hM G B A aBuf bBuf J Q R W3 W scratch euclidQ euclidR euclidW3
        lenB lenA loopDecrease heap right left readyBA with
      ⟨out, result, hrun, hresult, hgcd⟩
    exact ⟨out, result, by
      simpa [dense_upoly_zp_gcd_raw_ir, hswap] using hrun,
      hresult, by
        have hcomm : normalize (EuclideanDomain.gcd left right) =
            normalize (EuclideanDomain.gcd right left) :=
          normalize_eq_normalize
            (EuclideanDomain.dvd_gcd
              (EuclideanDomain.gcd_dvd_right left right)
              (EuclideanDomain.gcd_dvd_left left right))
            (EuclideanDomain.dvd_gcd
              (EuclideanDomain.gcd_dvd_right right left)
              (EuclideanDomain.gcd_dvd_left right left))
        rw [hcomm]
        exact hgcd⟩
  · rcases dense_upoly_zp_gcd_ordered_raw_ir_refines_of_ready this hcfg hp
        M hM G A B aBuf bBuf J Q R W3 W scratch euclidQ euclidR euclidW3
        lenA lenB loopDecrease heap left right readyAB with
      ⟨out, result, hrun, hresult, hgcd⟩
    exact ⟨out, result, by
      simpa [dense_upoly_zp_gcd_raw_ir, hswap] using hrun,
      hresult, hgcd⟩

structure DenseScalarRawResult where
  heap : RawHeap
  length : Nat

/-- Exact forward coefficient loop in `dense_upoly_zp::scalar_mul`. -/
def denseScalarMulLoop (this : DenseUPolyZp) (ptr : RawPtr UInt64)
    (length : Nat) (scalar : UInt64) (index : Nat)
    (heap : RawHeap) : RawExec RawHeap :=
  if index < length then
    match heap.readU64 ptr index with
    | .error fault => .error fault
    | .ok coeff =>
      match heap.writeU64 ptr index
          (Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this coeff scalar) with
      | .error fault => .error fault
      | .ok heap' => denseScalarMulLoop this ptr length scalar (index + 1) heap'
  else
    .ok heap
termination_by length - index
decreasing_by omega

/-- Exact zero/one/general dispatch of the source scalar multiplication. -/
def dense_upoly_zp_scalar_mul_raw_ir (this : DenseUPolyZp)
    (ptr : RawPtr UInt64) (length : Nat) (scalar : UInt64)
    (heap : RawHeap) : RawExec DenseScalarRawResult :=
  if scalar = 0 then
    .ok ⟨heap, 0⟩
  else if scalar = 1 then
    .ok ⟨heap, length⟩
  else
    match denseScalarMulLoop this ptr length scalar 0 heap with
    | .error fault => .error fault
    | .ok heap' => .ok ⟨heap', length⟩

theorem denseScalarMulLoop_succeeds (this : DenseUPolyZp)
    (ptr : RawPtr UInt64) (length : Nat) (scalar : UInt64)
    (index : Nat) (heap : RawHeap)
    (hvalid : heap.ValidU64Slice ptr length) :
    ∃ heap', denseScalarMulLoop this ptr length scalar index heap =
      .ok heap' := by
  rw [denseScalarMulLoop]
  split
  next hmore =>
    rcases heap.readU64_of_valid ptr length index hvalid hmore with
      ⟨coeff, hread⟩
    simp only [hread]
    let product := Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this coeff
      scalar
    rcases heap.writeU64_of_valid ptr length index product hvalid hmore with
      ⟨heap1, hwrite⟩
    simp only [product, hwrite]
    have hvalid1 : heap1.ValidU64Slice ptr length :=
      (RawHeap.writeU64_preserves_valid heap heap1 ptr index product hwrite
        ptr length).mp hvalid
    exact denseScalarMulLoop_succeeds this ptr length scalar (index + 1)
      heap1 hvalid1
  next hdone => exact ⟨heap, rfl⟩
termination_by length - index
decreasing_by omega

theorem denseScalarMulLoop_preserves_read_before (this : DenseUPolyZp)
    (ptr : RawPtr UInt64) (length : Nat) (scalar : UInt64)
    (index guard : Nat) (heap heap' : RawHeap) (value : UInt64)
    (hvalid : heap.ValidU64Slice ptr length) (hguard : guard < index)
    (hread : heap.readU64 ptr guard = .ok value)
    (hrun : denseScalarMulLoop this ptr length scalar index heap = .ok heap') :
    heap'.readU64 ptr guard = .ok value := by
  rw [denseScalarMulLoop] at hrun
  split at hrun
  next hmore =>
    rcases heap.readU64_of_valid ptr length index hvalid hmore with
      ⟨coeff, hcoeff⟩
    simp only [hcoeff] at hrun
    let product := Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this coeff
      scalar
    rcases heap.writeU64_of_valid ptr length index product hvalid hmore with
      ⟨heap1, hwrite⟩
    simp only [product, hwrite] at hrun
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 ptr index product
      hwrite
    have hvalid1 := (hlayout1 ptr length).mp hvalid
    have hread1 := RawHeap.readU64_writeU64_ne heap heap1 ptr ptr index guard
      product value hwrite hread (Or.inr (by omega))
    exact denseScalarMulLoop_preserves_read_before this ptr length scalar
      (index + 1) guard heap1 heap' value hvalid1 (by omega) hread1 hrun
  next hdone =>
    have heq : heap' = heap := Except.ok.inj hrun.symm
    simpa [heq] using hread
termination_by length - index
decreasing_by omega

/-- Every coefficient read before the concrete scalar loop is transformed by
the exact generated modular multiplication at the same raw cell. -/
theorem denseScalarMulLoop_reads (this : DenseUPolyZp)
    (ptr : RawPtr UInt64) (length : Nat) (scalar : UInt64)
    (index : Nat) (heap heap' : RawHeap)
    (hvalid : heap.ValidU64Slice ptr length)
    (hrun : denseScalarMulLoop this ptr length scalar index heap = .ok heap') :
    ∀ k coeff, index ≤ k → k < length →
      heap.readU64 ptr k = .ok coeff →
      heap'.readU64 ptr k = .ok
        (Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this coeff scalar) := by
  intro k original hindex hk hreadOriginal
  rw [denseScalarMulLoop] at hrun
  split at hrun
  next hmore =>
    rcases heap.readU64_of_valid ptr length index hvalid hmore with
      ⟨coeff, hcoeff⟩
    simp only [hcoeff] at hrun
    let product := Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this coeff
      scalar
    rcases heap.writeU64_of_valid ptr length index product hvalid hmore with
      ⟨heap1, hwrite⟩
    simp only [product, hwrite] at hrun
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 ptr index product
      hwrite
    have hvalid1 := (hlayout1 ptr length).mp hvalid
    by_cases heq : k = index
    · subst k
      have hcoeffEq : original = coeff :=
        Except.ok.inj (hreadOriginal.symm.trans hcoeff)
      subst original
      have hreadNow := RawHeap.readU64_writeU64_same heap heap1 ptr index
        product hwrite
      exact denseScalarMulLoop_preserves_read_before this ptr length scalar
        (index + 1) index heap1 heap' product hvalid1 (by omega) hreadNow hrun
    · have hread1 := RawHeap.readU64_writeU64_ne heap heap1 ptr ptr index k
          product original hwrite hreadOriginal (Or.inr (by omega))
      exact denseScalarMulLoop_reads this ptr length scalar (index + 1)
        heap1 heap' hvalid1 hrun k original (by omega) hk hread1
  next hdone => omega
termination_by length - index
decreasing_by omega

theorem denseScalarMulLoop_sameLayout (this : DenseUPolyZp)
    (ptr : RawPtr UInt64) (length : Nat) (scalar : UInt64)
    (index : Nat) (heap heap' : RawHeap)
    (hvalid : heap.ValidU64Slice ptr length)
    (hrun : denseScalarMulLoop this ptr length scalar index heap = .ok heap') :
    RawHeap.SameLayout heap heap' := by
  rw [denseScalarMulLoop] at hrun
  split at hrun
  next hmore =>
    rcases heap.readU64_of_valid ptr length index hvalid hmore with
      ⟨coeff, hcoeff⟩
    simp only [hcoeff] at hrun
    let product := Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this coeff
      scalar
    rcases heap.writeU64_of_valid ptr length index product hvalid hmore with
      ⟨heap1, hwrite⟩
    simp only [product, hwrite] at hrun
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 ptr index product
      hwrite
    have hvalid1 := (hlayout1 ptr length).mp hvalid
    have hlayout2 := denseScalarMulLoop_sameLayout this ptr length scalar
      (index + 1) heap1 heap' hvalid1 hrun
    exact fun other count =>
      (hlayout1 other count).trans (hlayout2 other count)
  next hdone =>
    have heq : heap' = heap := Except.ok.inj hrun.symm
    subst heap'
    exact fun _ _ => Iff.rfl
termination_by length - index
decreasing_by omega

/-- Semantic and residue invariant of the concrete general scalar loop. -/
theorem denseScalarMulLoop_refines (this : DenseUPolyZp)
    (ptr : RawPtr UInt64) (length : Nat) (scalar : UInt64)
    (heap heap' : RawHeap) (poly : Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 0 < this._p.toNat)
    (hvalid : heap.ValidU64Slice ptr length)
    (hcanonical : CanonicalU64Prefix heap ptr length this._p)
    (hrep : SlicePolyRep heap ptr length this._p.toNat poly)
    (hrun : denseScalarMulLoop this ptr length scalar 0 heap = .ok heap') :
    RawHeap.SameLayout heap heap' ∧
      RawCanonicalPolySlice this heap' ptr length
        (Polynomial.C (scalar.toNat : ZMod this._p.toNat) * poly) := by
  have hlayout := denseScalarMulLoop_sameLayout this ptr length scalar 0 heap
    heap' hvalid hrun
  have hvalid' := (hlayout ptr length).mp hvalid
  rcases slicePolyRep_exists_unique heap' ptr length this._p.toNat hvalid' with
    ⟨resultPoly, hresultRep, hunique⟩
  have hpolyEq : Polynomial.C (scalar.toNat : ZMod this._p.toNat) * poly =
      resultPoly := by
    ext degree
    by_cases hdegree : degree < length
    · rcases slicePolyRep_coeff heap ptr length this._p.toNat poly hrep
          degree hdegree with ⟨old, hreadOld, hcoeffOld⟩
      let product := Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this old
        scalar
      have hreadProduct := denseScalarMulLoop_reads this ptr length scalar 0
        heap heap' hvalid hrun degree old (Nat.zero_le _) hdegree hreadOld
      rcases slicePolyRep_coeff heap' ptr length this._p.toNat resultPoly
          hresultRep degree hdegree with
        ⟨observed, hreadObserved, hcoeffObserved⟩
      have hobserved : observed = product :=
        Except.ok.inj (hreadObserved.symm.trans hreadProduct)
      subst observed
      have hproduct := nmod_mul_ir_correct_of_configured this old scalar hcfg
        (hcanonical degree old hdegree hreadOld)
      rw [Polynomial.coeff_C_mul, hcoeffOld, hcoeffObserved]
      change (scalar.toNat : ZMod this._p.toNat) * old.toNat =
        (product.toNat : ZMod this._p.toNat)
      rw [hproduct]
      simp [ZMod.natCast_mod, mul_comm]
    · rw [slicePolyRep_coeff_zero_of_length_le heap' ptr length
          this._p.toNat resultPoly hresultRep degree (by omega)]
      have hzero := slicePolyRep_coeff_zero_of_length_le heap ptr length
        this._p.toNat poly hrep degree (by omega)
      simp [hzero]
  have hscaledRep : SlicePolyRep heap' ptr length this._p.toNat
      (Polynomial.C (scalar.toNat : ZMod this._p.toNat) * poly) := by
    rw [hpolyEq]
    exact hresultRep
  have hcanonical' : CanonicalU64Prefix heap' ptr length this._p := by
    intro degree observed hdegree hreadObserved
    rcases heap.readU64_of_valid ptr length degree hvalid hdegree with
      ⟨old, hreadOld⟩
    let product := Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this old
      scalar
    have hreadProduct := denseScalarMulLoop_reads this ptr length scalar 0
      heap heap' hvalid hrun degree old (Nat.zero_le _) hdegree hreadOld
    have hobserved : observed = product :=
      Except.ok.inj (hreadObserved.symm.trans hreadProduct)
    subst observed
    rw [nmod_mul_ir_correct_of_configured this old scalar hcfg
      (hcanonical degree old hdegree hreadOld)]
    exact Nat.mod_lt _ hp
  exact ⟨hlayout, hvalid', hcanonical', hscaledRep⟩

theorem dense_upoly_zp_scalar_mul_raw_ir_succeeds (this : DenseUPolyZp)
    (ptr : RawPtr UInt64) (length : Nat) (scalar : UInt64)
    (heap : RawHeap) (hvalid : heap.ValidU64Slice ptr length) :
    ∃ result, dense_upoly_zp_scalar_mul_raw_ir this ptr length scalar heap =
      .ok result := by
  by_cases hzero : scalar = 0
  · exact ⟨⟨heap, 0⟩, by
      simp [dense_upoly_zp_scalar_mul_raw_ir, hzero]⟩
  · by_cases hone : scalar = 1
    · exact ⟨⟨heap, length⟩, by
        simp [dense_upoly_zp_scalar_mul_raw_ir, hone]⟩
    · rcases denseScalarMulLoop_succeeds this ptr length scalar 0 heap
          hvalid with ⟨heap', hrun⟩
      exact ⟨⟨heap', length⟩, by
        simp [dense_upoly_zp_scalar_mul_raw_ir, hzero, hone, hrun]⟩

/-- The exact three-way source scalar dispatcher preserves the normalized raw
representation and denotes multiplication by the supplied canonical residue. -/
theorem dense_upoly_zp_scalar_mul_raw_ir_refines (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (ptr : RawPtr UInt64) (length : Nat) (scalar : UInt64)
    (heap : RawHeap) (poly : Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 0 < this._p.toNat)
    (hscalar : scalar.toNat < this._p.toNat)
    (hrep : RawDensePolyRep this heap ptr length poly) :
    ∃ result,
      dense_upoly_zp_scalar_mul_raw_ir this ptr length scalar heap =
        .ok result ∧
      RawDensePolyRep this result.heap ptr result.length
        (Polynomial.C (scalar.toNat : ZMod this._p.toNat) * poly) := by
  by_cases hzero : scalar = 0
  · subst scalar
    have hvalid0 : heap.ValidU64Slice ptr 0 :=
      heap.validU64Slice_mono ptr length 0 hrep.1 (Nat.zero_le _)
    have hzeroRep :=
      CLPoly.Impl.StrictHGCDRawRefinement.rawDensePolyRep_zero_length
        this heap ptr hvalid0
    exact ⟨⟨heap, 0⟩, by
      simp [dense_upoly_zp_scalar_mul_raw_ir],
      by simpa using hzeroRep⟩
  · by_cases hone : scalar = 1
    · subst scalar
      exact ⟨⟨heap, length⟩, by
        simp [dense_upoly_zp_scalar_mul_raw_ir], by simpa using hrep⟩
    · rcases denseScalarMulLoop_succeeds this ptr length scalar 0 heap hrep.1
          with ⟨heap', hrun⟩
      have hscaled := denseScalarMulLoop_refines this ptr length scalar heap
        heap' poly hcfg hp hrep.1 hrep.2.1 hrep.2.2.1 hrun
      have hnormalise : heap'.normaliseU64 ptr length = .ok length := by
        by_cases hlength : length = 0
        · subst length
          simp [denseScalarMulLoop] at hrun
          subst heap'
          exact hrep.2.2.2
        · have hlast := normaliseU64_poly_last_coeff_ne_zero heap ptr length
            this._p.toNat length poly hrep.1 hrep.2.2.1 hrep.2.1
            hrep.2.2.2 hlength
          have hscalarPos : 0 < scalar.toNat := by
            apply Nat.pos_of_ne_zero
            intro hscalarZero
            apply hzero
            apply UInt64.toNat_inj.mp
            simpa using hscalarZero
          have hscalarCast : (scalar.toNat : ZMod this._p.toNat) ≠ 0 := by
            intro hcast
            have hdvd : this._p.toNat ∣ scalar.toNat :=
              (ZMod.natCast_eq_zero_iff scalar.toNat this._p.toNat).mp hcast
            exact (Nat.not_dvd_of_pos_of_lt hscalarPos hscalar) hdvd
          have hlastScaled :
              (Polynomial.C (scalar.toNat : ZMod this._p.toNat) * poly).coeff
                  (length - 1) ≠ 0 := by
            rw [Polynomial.coeff_C_mul]
            exact mul_ne_zero hscalarCast hlast
          rcases slicePolyRep_coeff heap' ptr length this._p.toNat
              (Polynomial.C (scalar.toNat : ZMod this._p.toNat) * poly)
              hscaled.2.2.2 (length - 1) (by omega) with
            ⟨value, hread, hcoeff⟩
          have hvalue : value ≠ 0 := by
            intro hvalueZero
            subst value
            apply hlastScaled
            simpa using hcoeff.symm
          cases length with
          | zero => contradiction
          | succ n =>
            have hread' : heap'.readU64 ptr n = .ok value := by
              simpa using hread
            simp [RawHeap.normaliseU64, hread', hvalue]
      exact ⟨⟨heap', length⟩, by
        simp [dense_upoly_zp_scalar_mul_raw_ir, hzero, hone, hrun],
        hscaled.2.1, hscaled.2.2.1, hscaled.2.2.2, hnormalise⟩

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

/-- Dense vector length selected by the source sparse constructor. -/
def sparseDenseLength (sparse : SparsePolyZp) : Nat :=
  if h : 0 < sparse.size then sparse[0].1.deg + 1 else 0

/-- Object-level constructor entry with its source-derived vector length. -/
def sparse_upoly_zp_dense_constructor_raw_ir (ptr : RawPtr UInt64)
    (sparse : SparsePolyZp) (heap : RawHeap) : RawExec RawHeap :=
  sparse_upoly_zp_to_dense_raw_ir ptr (sparseDenseLength sparse) sparse heap

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

/-- A selected raw cell is framed by the remainder of the sparse iterator
when no remaining source term targets that degree. -/
theorem sparseToDenseWriteLoop_preserves_read (ptr : RawPtr UInt64)
    (length : Nat) (sparse : SparsePolyZp) (index guard : Nat)
    (heap heap' : RawHeap) (value : UInt64)
    (hvalid : heap.ValidU64Slice ptr length)
    (hdegree : ∀ i (hi : i < sparse.size), sparse[i].1.deg < length)
    (hne : ∀ i (hi : i < sparse.size), index ≤ i →
      sparse[i].1.deg ≠ guard)
    (hread : heap.readU64 ptr guard = .ok value)
    (hrun : sparseToDenseWriteLoop ptr length sparse index heap = .ok heap') :
    heap'.readU64 ptr guard = .ok value := by
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
    have hread1 := RawHeap.readU64_writeU64_ne heap heap1 ptr ptr
      term.1.deg guard term.2.val value hwrite hread (Or.inr (by
        have htarget := hne index hmore le_rfl
        intro habsolute
        apply htarget
        exact Nat.add_left_cancel habsolute))
    apply sparseToDenseWriteLoop_preserves_read ptr length sparse
      (index + 1) guard heap1 heap' value hvalid1 hdegree
    · intro i hi hle
      exact hne i hi (by omega)
    · exact hread1
    · exact hrun
  next hdone =>
    have heq : heap' = heap := Except.ok.inj hrun.symm
    simpa [heq] using hread
termination_by sparse.size - index
decreasing_by omega

/-- The first sparse term's concrete coefficient write remains observable at
the end of the iterator because canonical later degrees are strictly lower. -/
theorem sparseToDenseWriteLoop_leading_read (ptr : RawPtr UInt64)
    (length p : Nat) (sparse : SparsePolyZp) (heap heap' : RawHeap)
    (hnonempty : 0 < sparse.size)
    (hvalid : heap.ValidU64Slice ptr length)
    (hcanonical : SparsePolyZp.Canonical p sparse)
    (hdegree : ∀ i (hi : i < sparse.size), sparse[i].1.deg < length)
    (hrun : sparseToDenseWriteLoop ptr length sparse 0 heap = .ok heap') :
    heap'.readU64 ptr sparse[0].1.deg = .ok sparse[0].2.val := by
  rw [sparseToDenseWriteLoop] at hrun
  split at hrun
  next hmore =>
    let leading := sparse[0]
    have hleadingDegree : leading.1.deg < length := hdegree 0 hmore
    rcases heap.writeU64_of_valid ptr length leading.1.deg leading.2.val
        hvalid hleadingDegree with ⟨heap1, hwrite⟩
    simp only [leading, hwrite] at hrun
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 ptr
      leading.1.deg leading.2.val hwrite
    have hvalid1 := (hlayout1 ptr length).mp hvalid
    have hread1 := RawHeap.readU64_writeU64_same heap heap1 ptr
      leading.1.deg leading.2.val hwrite
    apply sparseToDenseWriteLoop_preserves_read ptr length sparse 1
      leading.1.deg heap1 heap' leading.2.val hvalid1 hdegree
    · intro i hi hpositive
      have hpairwise : List.Pairwise
          (fun a b : UMonomial × Zp => a.1.deg > b.1.deg) sparse.toList :=
        List.isChain_iff_pairwise.mp hcanonical.2.1
      have hrel := hpairwise.rel_get_of_lt
        (a := ⟨0, by simpa using hnonempty⟩)
        (b := ⟨i, by simpa using hi⟩) (by simpa using hpositive)
      have hzeroGet : sparse.toList[0] = sparse[0] :=
        Array.getElem_toList hmore
      have hiGet : sparse.toList[i] = sparse[i] :=
        Array.getElem_toList hi
      simp only [List.get_eq_getElem] at hrel
      rw [hzeroGet, hiGet] at hrel
      exact Nat.ne_of_lt hrel
    · exact hread1
    · exact hrun
  next hdone => omega

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

/-- For a nonempty canonical sparse input, the constructor length is one past
the first degree, so the real leading write makes raw normalization return the
entire dense vector. -/
theorem sparse_upoly_zp_to_dense_raw_ir_rawDense_nonempty
    (this : DenseUPolyZp) (ptr : RawPtr UInt64) (length : Nat)
    (sparse : SparsePolyZp) (heap : RawHeap)
    (hmodulus : this._p ≠ 0)
    (hvalid : heap.ValidU64Slice ptr length)
    (hcanonical : SparsePolyZp.Canonical this._p.toNat sparse)
    (hnonempty : 0 < sparse.size)
    (hlength : length = sparse[0].1.deg + 1) :
    ∃ heap', sparse_upoly_zp_to_dense_raw_ir ptr length sparse heap =
        .ok heap' ∧
      RawHeap.SameLayout heap heap' ∧
      RawDensePolyRep this heap' ptr length
        (SparsePolyZp.toPoly this._p.toNat sparse) := by
  have hdegree : ∀ i (hi : i < sparse.size),
      sparse[i].1.deg < length := by
    intro i hi
    by_cases hizero : i = 0
    · subst i
      omega
    · have hpairwise : List.Pairwise
          (fun a b : UMonomial × Zp => a.1.deg > b.1.deg) sparse.toList :=
        List.isChain_iff_pairwise.mp hcanonical.2.1
      have hrel := hpairwise.rel_get_of_lt
        (a := ⟨0, by simpa using hnonempty⟩)
        (b := ⟨i, by simpa using hi⟩) (by
          simp
          omega)
      simp only [List.get_eq_getElem] at hrel
      have hzeroGet : sparse.toList[0] = sparse[0] :=
        Array.getElem_toList hnonempty
      have hiGet : sparse.toList[i] = sparse[i] := Array.getElem_toList hi
      rw [hzeroGet, hiGet] at hrel
      omega
  rcases sparse_upoly_zp_to_dense_raw_ir_canonicalSlice this ptr length
      sparse heap hmodulus hvalid hcanonical hdegree with
    ⟨heap', hrun, hlayout, hslice⟩
  rcases CLPoly.Impl.StrictMulRefinement.mulZeroPadLoop_refines ptr 0 length 0
      this._p.toNat heap 0 this._p (by omega) hmodulus
      (by simpa using hvalid)
      (CLPoly.Impl.StrictHGCDRawRefinement.slicePolyRep_zero_length_any
        heap ptr this._p.toNat)
      (by
        intro i value hi
        omega) with ⟨heap1, hzero, hlayoutZero, _, _⟩
  have hvalid1 := (hlayoutZero ptr length).mp hvalid
  have hwrites : sparseToDenseWriteLoop ptr length sparse 0 heap1 =
      .ok heap' := by
    simpa [sparse_upoly_zp_to_dense_raw_ir, hzero] using hrun
  have hleadRead := sparseToDenseWriteLoop_leading_read ptr length
    this._p.toNat sparse heap1 heap' hnonempty hvalid1 hcanonical hdegree
    hwrites
  have hleadNonzero : sparse[0].2.val ≠ 0 := by
    apply hcanonical.2.2 sparse[0]
    simp
  have hnormalise : heap'.normaliseU64 ptr length = .ok length := by
    rw [hlength]
    simp [RawHeap.normaliseU64, hleadRead, hleadNonzero]
  exact ⟨heap', hrun, hlayout, hslice.1, hslice.2.1, hslice.2.2,
    hnormalise⟩

theorem sparse_upoly_zp_to_dense_raw_ir_rawDense_empty
    (this : DenseUPolyZp) (ptr : RawPtr UInt64) (heap : RawHeap)
    (hvalid : heap.ValidU64Slice ptr 0) :
    sparse_upoly_zp_to_dense_raw_ir ptr 0 #[] heap = .ok heap ∧
      RawDensePolyRep this heap ptr 0
        (SparsePolyZp.toPoly this._p.toNat #[]) := by
  constructor
  · simp [sparse_upoly_zp_to_dense_raw_ir, mulZeroPadLoop,
      sparseToDenseWriteLoop]
  · refine ⟨hvalid, ?_, ?_, rfl⟩
    · intro i value hi
      omega
    · simpa using
        (CLPoly.Impl.StrictHGCDRawRefinement.slicePolyRep_zero_length_any
          heap ptr this._p.toNat)

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

/-- The complete source-derived sparse constructor produces the strict input
relation consumed by the raw dense GCD path. -/
theorem sparse_upoly_zp_dense_constructor_raw_ir_result
    (this : DenseUPolyZp) (ptr : RawPtr UInt64) (sparse : SparsePolyZp)
    (heap : RawHeap) (hmodulus : this._p ≠ 0)
    (hvalid : heap.ValidU64Slice ptr (sparseDenseLength sparse))
    (hcanonical : SparsePolyZp.Canonical this._p.toNat sparse) :
    ∃ heap', sparse_upoly_zp_dense_constructor_raw_ir ptr sparse heap =
        .ok heap' ∧
      RawHeap.SameLayout heap heap' ∧
      SparseRawDenseRep this heap' ptr (sparseDenseLength sparse) sparse := by
  by_cases hnonempty : 0 < sparse.size
  · have hlength : sparseDenseLength sparse = sparse[0].1.deg + 1 := by
      simp [sparseDenseLength, hnonempty]
    rcases sparse_upoly_zp_to_dense_raw_ir_rawDense_nonempty this ptr
        (sparseDenseLength sparse) sparse heap hmodulus hvalid hcanonical
        hnonempty hlength with ⟨heap', hrun, hlayout, hdense⟩
    exact ⟨heap', by
      simpa [sparse_upoly_zp_dense_constructor_raw_ir] using hrun,
      hlayout, hcanonical, hdense⟩
  · have hsize : sparse.size = 0 := Nat.eq_zero_of_not_pos hnonempty
    have hsparse : sparse = #[] := by
      apply Array.ext
      · simp [hsize]
      · intro i hiSparse hiEmpty
        omega
    subst sparse
    rcases sparse_upoly_zp_to_dense_raw_ir_rawDense_empty this ptr heap
        (by simpa [sparseDenseLength] using hvalid) with ⟨hrun, hdense⟩
    exact ⟨heap, by
      simpa [sparse_upoly_zp_dense_constructor_raw_ir, sparseDenseLength]
        using hrun,
      fun _ _ => Iff.rfl, hcanonical, by
        simpa [sparseDenseLength] using hdense⟩

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

/-- A raw normalized GCD result is nonzero whenever either source input is
nonzero.  This discharges the concrete `lead()` access in
`__polynomial_GCD`; it is not an additional algorithm. -/
theorem rawDense_gcd_result_nonzero (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (heap : RawHeap) (ptr : RawPtr UInt64) (length : Nat)
    (left right result : Polynomial (ZMod this._p.toNat))
    (hleft : left ≠ 0)
    (hresult : RawDensePolyRep this heap ptr length result)
    (hgcd : normalize (EuclideanDomain.gcd left right) = normalize result) :
    result ≠ 0 ∧ length ≠ 0 := by
  have hresultNonzero : result ≠ 0 := by
    intro hresultZero
    have hgcdZero : EuclideanDomain.gcd left right = 0 := by
      apply normalize_eq_zero.mp
      rw [hgcd, hresultZero]
      simp
    have hdivides : (0 : Polynomial (ZMod this._p.toNat)) ∣ left := by
      simpa [hgcdZero] using EuclideanDomain.gcd_dvd_left left right
    exact hleft (zero_dvd_iff.mp hdivides)
  exact ⟨hresultNonzero, fun hlength =>
    hresultNonzero
      ((CLPoly.Impl.StrictHGCDRawRefinement.rawDensePolyRep_length_zero_iff
        this heap ptr length result hresult).mp hlength)⟩

/-- Observable success payload of the source `__polynomial_GCD` wrapper.
`output = none` is exactly the early degree-bound return, which leaves the
caller-owned sparse output unchanged. -/
structure InternalPolynomialGCDRawResult where
  heap : RawHeap
  denseLength : Nat
  returnCode : Int64
  output : Option SparsePolyZp

/-- The source tail after dense GCD has produced its normalized result. -/
def internalPolynomialGCDFinishRaw (this : DenseUPolyZp)
    (resultPtr : RawPtr UInt64) (gcdResult : DenseGcdRawResult)
    (lcGcd : UInt64) (degreeBound : Int64) :
    RawExec InternalPolynomialGCDRawResult :=
  match gcdResult.heap.readU64 resultPtr (gcdResult.lenG - 1) with
  | .error fault => .error fault
  | .ok lead =>
    let inverse := Generated.StrictGCD.dense_upoly_zp_nmod_inv_ir this lead
    let scale := Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this inverse
      lcGcd
    match dense_upoly_zp_scalar_mul_raw_ir this resultPtr gcdResult.lenG scale
        gcdResult.heap with
    | .error fault => .error fault
    | .ok scaled =>
      let degree : Int64 := (scaled.length - 1).toUInt64.toInt64
      if degree > degreeBound then
        .ok ⟨scaled.heap, scaled.length, (-1 : Int64), none⟩
      else
        match dense_upoly_zp_to_upoly_raw_ir this scaled.heap resultPtr
            scaled.length with
        | .error fault => .error fault
        | .ok sparse =>
          .ok ⟨scaled.heap, scaled.length, degree, some sparse⟩

theorem dense_upoly_zp_scalar_mul_raw_ir_length_of_ne_zero
    (this : DenseUPolyZp) (ptr : RawPtr UInt64) (length : Nat)
    (scalar : UInt64) (heap : RawHeap) (result : DenseScalarRawResult)
    (hscalar : scalar ≠ 0)
    (hrun : dense_upoly_zp_scalar_mul_raw_ir this ptr length scalar heap =
      .ok result) :
    result.length = length := by
  by_cases hone : scalar = 1
  · subst scalar
    simp [dense_upoly_zp_scalar_mul_raw_ir] at hrun
    rw [← hrun]
  · simp [dense_upoly_zp_scalar_mul_raw_ir, hscalar, hone] at hrun
    split at hrun
    · contradiction
    · have heq : (⟨_, length⟩ : DenseScalarRawResult) = result :=
        Except.ok.inj hrun
      rw [← heq]

/-- With the public wrapper's target leading coefficient `1`, every concrete
operation in the post-GCD tail succeeds and the degree-accepted branch emits
the actual reverse-scan result. -/
theorem internalPolynomialGCDFinishRaw_monic_succeeds
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (resultPtr : RawPtr UInt64) (gcdResult : DenseGcdRawResult)
    (degreeBound : Int64) (result : Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (hresult : RawDensePolyRep this gcdResult.heap resultPtr gcdResult.lenG
      result)
    (hresultNonzero : result ≠ 0)
    (hbound : ¬degreeBound < Int64.ofNat (gcdResult.lenG - 1)) :
    ∃ out sparse,
      internalPolynomialGCDFinishRaw this resultPtr gcdResult 1 degreeBound =
        .ok out ∧
      out.output = some sparse ∧
      out.denseLength = gcdResult.lenG ∧
      RawDenseSparseResult this out.heap resultPtr out.denseLength sparse ∧
      SparsePolyZp.toPoly this._p.toNat sparse = normalize result := by
  have hlength : gcdResult.lenG ≠ 0 := fun hzero =>
    hresultNonzero
      ((CLPoly.Impl.StrictHGCDRawRefinement.rawDensePolyRep_length_zero_iff
        this gcdResult.heap resultPtr gcdResult.lenG result hresult).mp hzero)
  rcases gcdResult.heap.readU64_of_valid resultPtr gcdResult.lenG
      (gcdResult.lenG - 1) hresult.1 (by omega) with ⟨lead, hreadLead⟩
  have hleadCoeff := normaliseU64_poly_last_coeff_ne_zero gcdResult.heap
    resultPtr gcdResult.lenG this._p.toNat gcdResult.lenG result hresult.1
    hresult.2.2.1 hresult.2.1 hresult.2.2.2 hlength
  rcases slicePolyRep_coeff gcdResult.heap resultPtr gcdResult.lenG
      this._p.toNat result hresult.2.2.1 (gcdResult.lenG - 1) (by omega) with
    ⟨observed, hreadObserved, hcoeffObserved⟩
  have hobserved : observed = lead :=
    Except.ok.inj (hreadObserved.symm.trans hreadLead)
  subst observed
  have hleadPos : 0 < lead.toNat := by
    apply Nat.pos_of_ne_zero
    intro hzeroNat
    apply hleadCoeff
    rw [hcoeffObserved]
    have hleadZero : lead = 0 := UInt64.toNat_inj.mp (by simpa using hzeroNat)
    simp [hleadZero]
  have hleadLt : lead.toNat < this._p.toNat :=
    hresult.2.1 (gcdResult.lenG - 1) lead (by omega) hreadLead
  let inverse := Generated.StrictGCD.dense_upoly_zp_nmod_inv_ir this lead
  have hinverse := dense_upoly_zp_nmod_inv_ir_correct this lead
    (Fact.out : Nat.Prime this._p.toNat) hleadPos hleadLt
  change inverse.toNat < this._p.toNat ∧
    (inverse.toNat : ZMod this._p.toNat) *
      (lead.toNat : ZMod this._p.toNat) = 1 at hinverse
  let scale := Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this inverse 1
  have hscaleNat := nmod_mul_ir_correct_of_configured this inverse 1 hcfg
    hinverse.1
  have hscaleLt : scale.toNat < this._p.toNat := by
    rw [hscaleNat]
    exact Nat.mod_lt _ (by omega)
  have hscaleCast : (scale.toNat : ZMod this._p.toNat) =
      (inverse.toNat : ZMod this._p.toNat) := by
    rw [hscaleNat]
    simp [ZMod.natCast_mod]
  have hinverseNonzero : (inverse.toNat : ZMod this._p.toNat) ≠ 0 := by
    intro hzero
    have hbad := hinverse.2
    rw [hzero, zero_mul] at hbad
    exact zero_ne_one hbad
  have hscaleCastNonzero : (scale.toNat : ZMod this._p.toNat) ≠ 0 := by
    intro hzero
    exact hinverseNonzero (hscaleCast.symm.trans hzero)
  have hscaleNonzero : scale ≠ 0 := by
    intro hzero
    apply hscaleCastNonzero
    simp [hzero]
  have hnatDegree :=
    CLPoly.Impl.StrictHGCDRawRefinement.rawDensePolyRep_natDegree_add_one
      this gcdResult.heap resultPtr gcdResult.lenG result hresult
      (Nat.pos_of_ne_zero hlength)
  have hleading : result.leadingCoeff =
      (lead.toNat : ZMod this._p.toNat) := by
    rw [← Polynomial.coeff_natDegree]
    rw [show result.natDegree = gcdResult.lenG - 1 by omega,
      hcoeffObserved]
  have hinverseEq : (inverse.toNat : ZMod this._p.toNat) =
      (lead.toNat : ZMod this._p.toNat)⁻¹ :=
    eq_inv_of_mul_eq_one_left hinverse.2
  have hscaledNormalize :
      Polynomial.C (scale.toNat : ZMod this._p.toNat) * result =
        normalize result := by
    rw [normalize_apply, Polynomial.coe_normUnit,
      CommGroupWithZero.coe_normUnit _
        (Polynomial.leadingCoeff_ne_zero.mpr hresultNonzero), hleading]
    rw [hscaleCast, hinverseEq]
    ac_rfl
  rcases dense_upoly_zp_scalar_mul_raw_ir_refines this resultPtr
      gcdResult.lenG scale gcdResult.heap result hcfg (by omega)
      hscaleLt hresult with ⟨scaled, hrunScaled, hscaledRep⟩
  have hscaledLength := dense_upoly_zp_scalar_mul_raw_ir_length_of_ne_zero
    this resultPtr gcdResult.lenG scale gcdResult.heap scaled hscaleNonzero
    hrunScaled
  have hscaledRep' : RawDensePolyRep this scaled.heap resultPtr
      scaled.length
      (Polynomial.C (scale.toNat : ZMod this._p.toNat) * result) := hscaledRep
  rcases dense_upoly_zp_to_upoly_raw_ir_result this scaled.heap resultPtr
      scaled.length
      (Polynomial.C (scale.toNat : ZMod this._p.toNat) * result) hscaledRep'
      with ⟨sparse, hrunSparse, hsparse⟩
  let out : InternalPolynomialGCDRawResult :=
    ⟨scaled.heap, scaled.length,
      (scaled.length - 1).toUInt64.toInt64, some sparse⟩
  have hrunSparse' : dense_upoly_zp_to_upoly_raw_ir this scaled.heap
      resultPtr gcdResult.lenG = .ok sparse := by
    simpa [hscaledLength] using hrunSparse
  have houtputSemantic : SparsePolyZp.toPoly this._p.toNat sparse =
      Polynomial.C (scale.toNat : ZMod this._p.toNat) * result := by
    rcases slicePolyRep_exists_unique scaled.heap resultPtr scaled.length
        this._p.toNat hsparse.dense.1 with ⟨observed, _, hunique⟩
    exact (hunique _ hsparse.dense.2.2.1).trans
      (hunique _ hscaledRep'.2.2.1).symm
  refine ⟨out, sparse, ?_, rfl, hscaledLength, hsparse,
    houtputSemantic.trans hscaledNormalize⟩
  simp [internalPolynomialGCDFinishRaw, hreadLead, inverse, scale,
    hrunScaled, hscaledLength, hbound, hrunSparse', out]

/-- Exact post-constructor control flow of `__polynomial_GCD`: dense GCD,
leading-cell read, generated inverse and multiply, concrete scalar loop,
signed degree guard, then the real reverse dense scan. -/
def internal_polynomial_GCD_postconstruct_raw_ir (this : DenseUPolyZp)
    (M : HgcdMat) (hM : M.Valid)
    (resultPtr leftPtr rightPtr aBuf bBuf J Q R : RawPtr UInt64)
    (W3 : RawPtr Word3) (W scratch : RawPtr UInt64)
    (euclidQ euclidR : RawPtr UInt64) (euclidW3 : RawPtr Word3)
    (lenLeft lenRight : Nat) (lcGcd : UInt64) (degreeBound : Int64)
    (loopDecrease : Generated.StrictGCDHGCD.HgcdGcdLoopLengthDecreases
      this M hM W scratch)
    (heap : RawHeap) : RawExec InternalPolynomialGCDRawResult :=
  match dense_upoly_zp_gcd_raw_ir this M hM resultPtr leftPtr rightPtr aBuf
      bBuf J Q R W3 W scratch euclidQ euclidR euclidW3 lenLeft lenRight
      loopDecrease heap with
  | .error fault => .error fault
  | .ok gcdResult => internalPolynomialGCDFinishRaw this resultPtr gcdResult
      lcGcd degreeBound

/-- Composition bridge from one witnessed real dense-GCD execution through
the complete public-monic tail. -/
theorem internal_polynomial_GCD_postconstruct_raw_ir_monic_refines
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (M : HgcdMat) (hM : M.Valid)
    (resultPtr leftPtr rightPtr aBuf bBuf J Q R : RawPtr UInt64)
    (W3 : RawPtr Word3) (W scratch : RawPtr UInt64)
    (euclidQ euclidR : RawPtr UInt64) (euclidW3 : RawPtr Word3)
    (lenLeft lenRight : Nat) (degreeBound : Int64)
    (loopDecrease : Generated.StrictGCDHGCD.HgcdGcdLoopLengthDecreases
      this M hM W scratch)
    (heap : RawHeap) (gcdResult : DenseGcdRawResult)
    (left right result : Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (hrunGcd : dense_upoly_zp_gcd_raw_ir this M hM resultPtr leftPtr rightPtr
      aBuf bBuf J Q R W3 W scratch euclidQ euclidR euclidW3 lenLeft lenRight
      loopDecrease heap = .ok gcdResult)
    (hresult : RawDensePolyRep this gcdResult.heap resultPtr gcdResult.lenG
      result)
    (hgcd : normalize (EuclideanDomain.gcd left right) = normalize result)
    (hleftNonzero : left ≠ 0)
    (hbound : ¬degreeBound < Int64.ofNat (gcdResult.lenG - 1)) :
    ∃ out sparse,
      internal_polynomial_GCD_postconstruct_raw_ir this M hM resultPtr
          leftPtr rightPtr aBuf bBuf J Q R W3 W scratch euclidQ euclidR
          euclidW3 lenLeft lenRight 1 degreeBound loopDecrease heap = .ok out ∧
      out.output = some sparse ∧
      RawDenseSparseResult this out.heap resultPtr out.denseLength sparse ∧
      SparsePolyZp.toPoly this._p.toNat sparse =
        normalize (EuclideanDomain.gcd left right) := by
  have hnonzero := rawDense_gcd_result_nonzero this gcdResult.heap resultPtr
    gcdResult.lenG left right result hleftNonzero hresult hgcd
  rcases internalPolynomialGCDFinishRaw_monic_succeeds this resultPtr
      gcdResult degreeBound result hcfg hp hresult hnonzero.1 hbound with
    ⟨out, sparse, hfinish, houtput, _, hsparse, hsemantic⟩
  exact ⟨out, sparse, by
    simp [internal_polynomial_GCD_postconstruct_raw_ir, hrunGcd, hfinish],
    houtput, hsparse, hsemantic.trans hgcd.symm⟩

end CLPoly.Impl.StrictPolynomialGCDRefinement
