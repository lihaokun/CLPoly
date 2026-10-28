import CLPoly.Impl.StrictMulRefinement

set_option autoImplicit false
set_option maxRecDepth 4000

namespace CLPoly.Impl.StrictMulRefinement

open Generated.StrictMul
open CLPoly.Impl.StrictWordArithmetic
open CLPoly.Impl.StrictDivremRefinement
open CLPoly.Impl.RawPolynomialRep

theorem karMul_refines_slice (this : DenseUPolyZp)
    (C A B : RawPtr UInt64) (n : Nat) (scratch : RawPtr UInt64)
    (heap : RawHeap)
    (left right : Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this)
    (hp : 1 < this._p.toNat)
    (hn : 0 < n) (hNWord : n < limbBase)
    (hC : heap.ValidU64Slice C (2 * n - 1))
    (hA : heap.ValidU64Slice A n)
    (hB : heap.ValidU64Slice B n)
    (hScratch : heap.ValidU64Slice scratch (karScratchNeed n))
    (hCA : U64SlicesDisjoint C (2 * n - 1) A n)
    (hCB : U64SlicesDisjoint C (2 * n - 1) B n)
    (hCScratch : U64SlicesDisjoint C (2 * n - 1) scratch
      (karScratchNeed n))
    (hScratchA : U64SlicesDisjoint scratch (karScratchNeed n) A n)
    (hScratchB : U64SlicesDisjoint scratch (karScratchNeed n) B n)
    (hCanonicalA : CanonicalU64Prefix heap A n this._p)
    (hCanonicalB : CanonicalU64Prefix heap B n this._p)
    (hRepA : SlicePolyRep heap A n this._p.toNat left)
    (hRepB : SlicePolyRep heap B n this._p.toNat right) :
    ∃ heap', dense_upoly_zp__kar_mul_ir this C A B n scratch heap =
        .ok heap' ∧ RawHeap.SameLayout heap heap' ∧
      SlicePolyRep heap' C (2 * n - 1) this._p.toNat (left * right) ∧
      CanonicalU64Prefix heap' C (2 * n - 1) this._p := by
  rw [dense_upoly_zp__kar_mul_ir]
  split
  next hbase =>
    rcases classicalMul_refines_slice this C A n B n heap left right hcfg hp
      hn hn hNWord (by simpa [two_mul] using hC) hA hB
      (by simpa [two_mul] using hCA) (by simpa [two_mul] using hCB)
      hCanonicalA hCanonicalB hRepA hRepB with
      ⟨heap', hrun, hlayout, hrep, hcanonical⟩
    exact ⟨heap', hrun, hlayout, by simpa [two_mul] using hrep,
      by simpa [two_mul] using hcanonical⟩
  next hrec =>
    let m := n / 2
    let h := n - m
    let t1 := scratch
    let t2 := t1.add h
    let sP0 := t2.add h
    let sP1 := sP0.add (2 * m - 1)
    let recScratch := sP1.add (2 * h - 1)
    have hn16 : 16 ≤ n := by omega
    have hmPos : 0 < m := by
      dsimp [m]
      exact Nat.div_pos (by omega) (by omega)
    have hhPos : 0 < h := by dsimp [h, m]; omega
    have hmLt : m < n := by
      simpa [m] using (kar_split_children_lt n hn16).1
    have hhLt : h < n := by
      simpa [h, m] using (kar_split_children_lt n hn16).2
    have hshape : h = m ∨ h = m + 1 := by
      simpa [m, h] using kar_split_shape n
    have hmn : m + h = n := by dsimp [m, h]; omega
    have hAfull : heap.ValidU64Slice A (m + h) := by simpa [hmn] using hA
    have hBfull : heap.ValidU64Slice B (m + h) := by simpa [hmn] using hB
    have hCanonicalAfull : CanonicalU64Prefix heap A (m + h) this._p := by
      simpa [hmn] using hCanonicalA
    have hCanonicalBfull : CanonicalU64Prefix heap B (m + h) this._p := by
      simpa [hmn] using hCanonicalB
    have hRepAfull : SlicePolyRep heap A (m + h) this._p.toNat left := by
      simpa [hmn] using hRepA
    have hRepBfull : SlicePolyRep heap B (m + h) this._p.toNat right := by
      simpa [hmn] using hRepB
    rcases splitSlicePolyRepCanonical heap A m h this._p.toNat left this._p
      hAfull hRepAfull hCanonicalAfull with
      ⟨leftLow, leftHigh, hRepALow, hRepAHigh, hCanonicalALow,
        hCanonicalAHigh, hleftSplit⟩
    rcases splitSlicePolyRepCanonical heap B m h this._p.toNat right this._p
      hBfull hRepBfull hCanonicalBfull with
      ⟨rightLow, rightHigh, hRepBLow, hRepBHigh, hCanonicalBLow,
        hCanonicalBHigh, hrightSplit⟩
    rcases karScratchSlices heap scratch n hn16 hScratch with
      ⟨hT1, hT2, hP0, hP1, hRec⟩
    rcases karScratchNested_pairwise scratch m h
      (max (karScratchNeed m) (karScratchNeed h)) with
      ⟨hT1T2, hT1P0, hT1P1, hT2P0, hT2P1, hP0P1,
        hT1Rec, hT2Rec, hP0Rec, hP1Rec⟩
    have hT1A : U64SlicesDisjoint t1 h A (m + h) := by
      dsimp [t1]
      apply u64SlicesDisjoint_mono hScratchA (smallLeft := h)
        (smallRight := m + h)
      · rw [karScratchNeed_step n (by omega)]
        dsimp [m, h]
        omega
      · omega
    have hT1B : U64SlicesDisjoint t1 h B (m + h) := by
      dsimp [t1]
      apply u64SlicesDisjoint_mono hScratchB (smallLeft := h)
        (smallRight := m + h)
      · rw [karScratchNeed_step n (by omega)]
        dsimp [m, h]
        omega
      · omega
    have hT2A := u64SlicesDisjoint_add_left hScratchA
      (start := h) (count := h) (by
        rw [karScratchNeed_step n (by omega)]
        dsimp [m, h]
        omega)
    have hT2Afull : U64SlicesDisjoint t2 h A (m + h) := by
      exact u64SlicesDisjoint_mono hT2A (Nat.le_refl h) (by omega)
    have hT2B := u64SlicesDisjoint_add_left hScratchB
      (start := h) (count := h) (by
        rw [karScratchNeed_step n (by omega)]
        dsimp [m, h]
        omega)
    have hT2Bfull : U64SlicesDisjoint t2 h B (m + h) := by
      exact u64SlicesDisjoint_mono hT2B (Nat.le_refl h) (by omega)
    rcases karPrepareHalves_ok this A B t1 t2 m h heap (by omega)
      hAfull hBfull hT1 hT2 with ⟨heap2, hprep, hlay2⟩
    rcases karPrepareHalves_refines_and_preserves_inputs this A B t1 t2 m h
      heap heap2 left right hp hshape hAfull hBfull hT1 hT2 hT1T2
      hT1A hT1B hT2Afull hT2Bfull hCanonicalAfull hCanonicalBfull
      hRepAfull hRepBfull hprep with
      ⟨_, hsameA2, hsameB2, hRepA2, hRepB2, hCanonicalA2, hCanonicalB2,
        hRepT12, hRepT22, hCanonicalT12, hCanonicalT22⟩
    have hAm := heap.validU64Slice_mono A (m + h) m hAfull (by omega)
    have hBm := heap.validU64Slice_mono B (m + h) m hBfull (by omega)
    have hAh := heap.validU64Slice_add A (m + h) m h hAfull (by omega)
    have hBh := heap.validU64Slice_add B (m + h) m h hBfull (by omega)
    have hAm2 := (hlay2 A m).mp hAm
    have hBm2 := (hlay2 B m).mp hBm
    have hAh2 := (hlay2 (A.add m) h).mp hAh
    have hBh2 := (hlay2 (B.add m) h).mp hBh
    have hsameALow2 := sameU64Prefix_subslice hsameA2
      (start := 0) (count := m) (by omega)
    have hsameBLow2 := sameU64Prefix_subslice hsameB2
      (start := 0) (count := m) (by omega)
    have hsameAHigh2 := sameU64Prefix_subslice hsameA2
      (start := m) (count := h) (by omega)
    have hsameBHigh2 := sameU64Prefix_subslice hsameB2
      (start := m) (count := h) (by omega)
    have hRepALow2 : SlicePolyRep heap2 A m this._p.toNat leftLow := by
      simpa using slicePolyRep_of_same_prefix heap heap2 (A.add 0) m
        this._p.toNat leftLow (by simpa using hAm) (by simpa using hAm2)
        hsameALow2 (by simpa using hRepALow)
    have hRepBLow2 : SlicePolyRep heap2 B m this._p.toNat rightLow := by
      simpa using slicePolyRep_of_same_prefix heap heap2 (B.add 0) m
        this._p.toNat rightLow (by simpa using hBm) (by simpa using hBm2)
        hsameBLow2 (by simpa using hRepBLow)
    have hRepAHigh2 := slicePolyRep_of_same_prefix heap heap2 (A.add m) h
      this._p.toNat leftHigh hAh hAh2 hsameAHigh2 hRepAHigh
    have hRepBHigh2 := slicePolyRep_of_same_prefix heap heap2 (B.add m) h
      this._p.toNat rightHigh hBh hBh2 hsameBHigh2 hRepBHigh
    have hCanonicalALow2 : CanonicalU64Prefix heap2 A m this._p := by
      simpa using canonicalU64Prefix_of_same_prefix heap heap2 (A.add 0) m
        this._p (by simpa using hAm) hsameALow2 (by simpa using hCanonicalALow)
    have hCanonicalBLow2 : CanonicalU64Prefix heap2 B m this._p := by
      simpa using canonicalU64Prefix_of_same_prefix heap heap2 (B.add 0) m
        this._p (by simpa using hBm) hsameBLow2 (by simpa using hCanonicalBLow)
    have hCanonicalAHigh2 := canonicalU64Prefix_of_same_prefix heap heap2
      (A.add m) h this._p hAh hsameAHigh2 hCanonicalAHigh
    have hCanonicalBHigh2 := canonicalU64Prefix_of_same_prefix heap heap2
      (B.add m) h this._p hBh hsameBHigh2 hCanonicalBHigh
    have hRecM := heap.validU64Slice_mono recScratch
      (max (karScratchNeed m) (karScratchNeed h)) (karScratchNeed m)
      hRec (Nat.le_max_left _ _)
    have hP0A : U64SlicesDisjoint sP0 (2 * m - 1) A m := by
      have hx := u64SlicesDisjoint_subslices hScratchA
        (leftStart := 2 * h) (leftCount := 2 * m - 1)
        (rightStart := 0) (rightCount := m) (by
          rw [karScratchNeed_step n (by omega)]
          dsimp [m, h]
          omega) (by omega)
      simpa [sP0, t2, t1, rawPtr_add_add, two_mul] using hx
    have hP0B : U64SlicesDisjoint sP0 (2 * m - 1) B m := by
      have hx := u64SlicesDisjoint_subslices hScratchB
        (leftStart := 2 * h) (leftCount := 2 * m - 1)
        (rightStart := 0) (rightCount := m) (by
          rw [karScratchNeed_step n (by omega)]
          dsimp [m, h]
          omega) (by omega)
      simpa [sP0, t2, t1, rawPtr_add_add, two_mul] using hx
    have hP0RecM := u64SlicesDisjoint_mono hP0Rec
      (Nat.le_refl (2 * m - 1)) (Nat.le_max_left _ _)
    rcases karScratchSlices_disjoint_guard scratch A n m hn16
      (u64SlicesDisjoint_mono hScratchA (Nat.le_refl _) (by omega)) with
      ⟨_, _, _, _, hRecA⟩
    rcases karScratchSlices_disjoint_guard scratch B n m hn16
      (u64SlicesDisjoint_mono hScratchB (Nat.le_refl _) (by omega)) with
      ⟨_, _, _, _, hRecB⟩
    have hRecMA := u64SlicesDisjoint_mono hRecA
      (Nat.le_max_left _ _) (Nat.le_refl m)
    have hRecMB := u64SlicesDisjoint_mono hRecB
      (Nat.le_max_left _ _) (Nat.le_refl m)
    rcases karMul_refines_slice this sP0 A B m recScratch heap2
      leftLow rightLow hcfg hp hmPos (by omega)
      ((hlay2 sP0 (2 * m - 1)).mp hP0) hAm2 hBm2
      ((hlay2 recScratch (karScratchNeed m)).mp hRecM)
      hP0A hP0B hP0RecM hRecMA hRecMB hCanonicalALow2
      hCanonicalBLow2 hRepALow2 hRepBLow2 with
      ⟨heap3, hp0, hlay3, hRepP03, hCanonicalP03⟩
    have hRecMT1 : U64SlicesDisjoint recScratch (karScratchNeed m) t1 h :=
      u64SlicesDisjoint_symm (u64SlicesDisjoint_mono hT1Rec
        (Nat.le_refl h) (Nat.le_max_left _ _))
    have hRecMT2 : U64SlicesDisjoint recScratch (karScratchNeed m) t2 h :=
      u64SlicesDisjoint_symm (u64SlicesDisjoint_mono hT2Rec
        (Nat.le_refl h) (Nat.le_max_left _ _))
    rcases karMul_preserves_slice_canonical this sP0 A B m recScratch t1 h
      this._p.toNat heap2 heap3 (karPreparedPoly left m h) this._p hmPos
      ((hlay2 sP0 (2 * m - 1)).mp hP0) hAm2 hBm2
      ((hlay2 recScratch (karScratchNeed m)).mp hRecM)
      ((hlay2 t1 h).mp hT1) (u64SlicesDisjoint_symm hT1P0) hRecMT1
      hRepT12 hCanonicalT12 hp0 with
      ⟨_, hRepT13, hCanonicalT13⟩
    rcases karMul_preserves_slice_canonical this sP0 A B m recScratch t2 h
      this._p.toNat heap2 heap3 (karPreparedPoly right m h) this._p hmPos
      ((hlay2 sP0 (2 * m - 1)).mp hP0) hAm2 hBm2
      ((hlay2 recScratch (karScratchNeed m)).mp hRecM)
      ((hlay2 t2 h).mp hT2) (u64SlicesDisjoint_symm hT2P0) hRecMT2
      hRepT22 hCanonicalT22 hp0 with
      ⟨_, hRepT23, hCanonicalT23⟩
    have hScratchAHigh := u64SlicesDisjoint_add_right hScratchA
      (start := m) (count := h) (by omega)
    have hScratchBHigh := u64SlicesDisjoint_add_right hScratchB
      (start := m) (count := h) (by omega)
    rcases karScratchSlices_disjoint_guard scratch (A.add m) n h hn16
      hScratchAHigh with ⟨_, _, hP0AHigh, _, hRecAHigh⟩
    rcases karScratchSlices_disjoint_guard scratch (B.add m) n h hn16
      hScratchBHigh with ⟨_, _, hP0BHigh, _, hRecBHigh⟩
    have hRecMAHigh := u64SlicesDisjoint_mono hRecAHigh
      (Nat.le_max_left _ _) (Nat.le_refl h)
    have hRecMBHigh := u64SlicesDisjoint_mono hRecBHigh
      (Nat.le_max_left _ _) (Nat.le_refl h)
    rcases karMul_preserves_slice_canonical this sP0 A B m recScratch
      (A.add m) h this._p.toNat heap2 heap3 leftHigh this._p hmPos
      ((hlay2 sP0 (2 * m - 1)).mp hP0) hAm2 hBm2
      ((hlay2 recScratch (karScratchNeed m)).mp hRecM) hAh2
      hP0AHigh hRecMAHigh hRepAHigh2 hCanonicalAHigh2 hp0 with
      ⟨_, hRepAHigh3, hCanonicalAHigh3⟩
    rcases karMul_preserves_slice_canonical this sP0 A B m recScratch
      (B.add m) h this._p.toNat heap2 heap3 rightHigh this._p hmPos
      ((hlay2 sP0 (2 * m - 1)).mp hP0) hAm2 hBm2
      ((hlay2 recScratch (karScratchNeed m)).mp hRecM) hBh2
      hP0BHigh hRecMBHigh hRepBHigh2 hCanonicalBHigh2 hp0 with
      ⟨_, hRepBHigh3, hCanonicalBHigh3⟩
    have hRecH := heap.validU64Slice_mono recScratch
      (max (karScratchNeed m) (karScratchNeed h)) (karScratchNeed h)
      hRec (Nat.le_max_right _ _)
    have hRecH3 := (hlay3 recScratch (karScratchNeed h)).mp
      ((hlay2 recScratch (karScratchNeed h)).mp hRecH)
    have hP13 := (hlay3 sP1 (2 * h - 1)).mp
      ((hlay2 sP1 (2 * h - 1)).mp hP1)
    have hT13 := (hlay3 t1 h).mp ((hlay2 t1 h).mp hT1)
    have hT23 := (hlay3 t2 h).mp ((hlay2 t2 h).mp hT2)
    have hP1RecH := u64SlicesDisjoint_mono hP1Rec
      (Nat.le_refl (2 * h - 1)) (Nat.le_max_right _ _)
    have hRecHT1 : U64SlicesDisjoint recScratch (karScratchNeed h) t1 h :=
      u64SlicesDisjoint_symm (u64SlicesDisjoint_mono hT1Rec
        (Nat.le_refl h) (Nat.le_max_right _ _))
    have hRecHT2 : U64SlicesDisjoint recScratch (karScratchNeed h) t2 h :=
      u64SlicesDisjoint_symm (u64SlicesDisjoint_mono hT2Rec
        (Nat.le_refl h) (Nat.le_max_right _ _))
    rcases karMul_refines_slice this sP1 t1 t2 h recScratch heap3
      (karPreparedPoly left m h) (karPreparedPoly right m h) hcfg hp hhPos
      (by omega) hP13 hT13 hT23 hRecH3
      (u64SlicesDisjoint_symm hT1P1)
      (u64SlicesDisjoint_symm hT2P1) hP1RecH hRecHT1 hRecHT2
      hCanonicalT13 hCanonicalT23 hRepT13 hRepT23 with
      ⟨heap4, hp1, hlay4, hRepP14, hCanonicalP14⟩
    have hRecHP0 : U64SlicesDisjoint recScratch (karScratchNeed h)
        sP0 (2 * m - 1) :=
      u64SlicesDisjoint_symm (u64SlicesDisjoint_mono hP0Rec
        (Nat.le_refl (2 * m - 1)) (Nat.le_max_right _ _))
    have hP03 := (hlay3 sP0 (2 * m - 1)).mp
      ((hlay2 sP0 (2 * m - 1)).mp hP0)
    rcases karMul_preserves_slice_canonical this sP1 t1 t2 h recScratch
      sP0 (2 * m - 1) this._p.toNat heap3 heap4
      (leftLow * rightLow) this._p hhPos hP13 hT13 hT23 hRecH3 hP03
      (u64SlicesDisjoint_symm hP0P1) hRecHP0 hRepP03 hCanonicalP03 hp1 with
      ⟨_, hRepP04, hCanonicalP04⟩
    have hAh3 := (hlay3 (A.add m) h).mp hAh2
    have hBh3 := (hlay3 (B.add m) h).mp hBh2
    rcases karScratchSlices_disjoint_guard scratch (A.add m) n h hn16
      hScratchAHigh with ⟨_, _, _, hP1AHigh, hRecAHigh'⟩
    rcases karScratchSlices_disjoint_guard scratch (B.add m) n h hn16
      hScratchBHigh with ⟨_, _, _, hP1BHigh, hRecBHigh'⟩
    have hRecHAHigh := u64SlicesDisjoint_mono hRecAHigh'
      (Nat.le_max_right _ _) (Nat.le_refl h)
    have hRecHBHigh := u64SlicesDisjoint_mono hRecBHigh'
      (Nat.le_max_right _ _) (Nat.le_refl h)
    rcases karMul_preserves_slice_canonical this sP1 t1 t2 h recScratch
      (A.add m) h this._p.toNat heap3 heap4 leftHigh this._p hhPos
      hP13 hT13 hT23 hRecH3 hAh3 hP1AHigh hRecHAHigh
      hRepAHigh3 hCanonicalAHigh3 hp1 with
      ⟨_, hRepAHigh4, hCanonicalAHigh4⟩
    rcases karMul_preserves_slice_canonical this sP1 t1 t2 h recScratch
      (B.add m) h this._p.toNat heap3 heap4 rightHigh this._p hhPos
      hP13 hT13 hT23 hRecH3 hBh3 hP1BHigh hRecHBHigh
      hRepBHigh3 hCanonicalBHigh3 hp1 with
      ⟨_, hRepBHigh4, hCanonicalBHigh4⟩
    have hCHigh := heap.validU64Slice_add C (2 * n - 1) (2 * m)
      (2 * h - 1) hC (by dsimp [m, h]; omega)
    have hCHigh4 := (hlay4 (C.add (2 * m)) (2 * h - 1)).mp
      ((hlay3 (C.add (2 * m)) (2 * h - 1)).mp
        ((hlay2 (C.add (2 * m)) (2 * h - 1)).mp hCHigh))
    have hAh4 := (hlay4 (A.add m) h).mp hAh3
    have hBh4 := (hlay4 (B.add m) h).mp hBh3
    have hRecH4 := (hlay4 recScratch (karScratchNeed h)).mp hRecH3
    have hHighA := u64SlicesDisjoint_subslices hCA
      (leftStart := 2 * m) (leftCount := 2 * h - 1)
      (rightStart := m) (rightCount := h) (by dsimp [m, h]; omega)
      (by omega)
    have hHighB := u64SlicesDisjoint_subslices hCB
      (leftStart := 2 * m) (leftCount := 2 * h - 1)
      (rightStart := m) (rightCount := h) (by dsimp [m, h]; omega)
      (by omega)
    have hScratchHigh := u64SlicesDisjoint_add_right
      (u64SlicesDisjoint_symm hCScratch)
      (start := 2 * m) (count := 2 * h - 1) (by dsimp [m, h]; omega)
    rcases karScratchSlices_disjoint_guard scratch (C.add (2 * m)) n
      (2 * h - 1) hn16 hScratchHigh with
      ⟨_, _, _, _, hRecHigh⟩
    have hHighRecH : U64SlicesDisjoint (C.add (2 * m)) (2 * h - 1)
        recScratch (karScratchNeed h) :=
      u64SlicesDisjoint_symm (u64SlicesDisjoint_mono hRecHigh
        (Nat.le_max_right _ _) (Nat.le_refl (2 * h - 1)))
    rcases karMul_refines_slice this (C.add (2 * m)) (A.add m)
      (B.add m) h recScratch heap4 leftHigh rightHigh hcfg hp hhPos
      (by omega) hCHigh4 hAh4 hBh4 hRecH4 hHighA hHighB hHighRecH
      hRecHAHigh hRecHBHigh hCanonicalAHigh4 hCanonicalBHigh4
      hRepAHigh4 hRepBHigh4 with
      ⟨heap5, phigh, hlay5, hRepHigh5, hCanonicalHigh5⟩
    have hHighP0 : U64SlicesDisjoint (C.add (2 * m)) (2 * h - 1)
        sP0 (2 * m - 1) := by
      have hx := u64SlicesDisjoint_subslices hCScratch
        (leftStart := 2 * m) (leftCount := 2 * h - 1)
        (rightStart := 2 * h) (rightCount := 2 * m - 1)
        (by dsimp [m, h]; omega) (by
          rw [karScratchNeed_step n (by omega)]
          dsimp [m, h]
          omega)
      simpa [sP0, t2, t1, rawPtr_add_add, two_mul] using hx
    have hHighP1 : U64SlicesDisjoint (C.add (2 * m)) (2 * h - 1)
        sP1 (2 * h - 1) := by
      have hx := u64SlicesDisjoint_subslices hCScratch
        (leftStart := 2 * m) (leftCount := 2 * h - 1)
        (rightStart := 2 * h + (2 * m - 1))
        (rightCount := 2 * h - 1) (by dsimp [m, h]; omega) (by
          rw [karScratchNeed_step n (by omega)]
          dsimp [m, h]
          omega)
      simpa [sP1, sP0, t2, t1, rawPtr_add_add, two_mul,
        Nat.add_assoc] using hx
    have hRecHP0' : U64SlicesDisjoint recScratch (karScratchNeed h)
        sP0 (2 * m - 1) :=
      u64SlicesDisjoint_symm (u64SlicesDisjoint_mono hP0Rec
        (Nat.le_refl (2 * m - 1)) (Nat.le_max_right _ _))
    have hRecHP1' : U64SlicesDisjoint recScratch (karScratchNeed h)
        sP1 (2 * h - 1) :=
      u64SlicesDisjoint_symm (u64SlicesDisjoint_mono hP1Rec
        (Nat.le_refl (2 * h - 1)) (Nat.le_max_right _ _))
    have hP04 := (hlay4 sP0 (2 * m - 1)).mp hP03
    rcases karMul_preserves_slice_canonical this (C.add (2 * m))
      (A.add m) (B.add m) h recScratch sP0 (2 * m - 1) this._p.toNat
      heap4 heap5 (leftLow * rightLow) this._p hhPos hCHigh4 hAh4 hBh4
      hRecH4 hP04 hHighP0 hRecHP0' hRepP04 hCanonicalP04 phigh with
      ⟨_, hRepP05, hCanonicalP05⟩
    have hP14 := (hlay4 sP1 (2 * h - 1)).mp hP13
    rcases karMul_preserves_slice_canonical this (C.add (2 * m))
      (A.add m) (B.add m) h recScratch sP1 (2 * h - 1) this._p.toNat
      heap4 heap5
      (karPreparedPoly left m h * karPreparedPoly right m h) this._p
      hhPos hCHigh4 hAh4 hBh4 hRecH4 hP14 hHighP1 hRecHP1'
      hRepP14 hCanonicalP14 phigh with
      ⟨_, hRepP15, hCanonicalP15⟩
    let tailLength := (2 * h - 1) - (2 * m - 1)
    have hP1Length : (2 * m - 1) + tailLength = 2 * h - 1 := by
      dsimp [tailLength]
      rcases hshape with heq | heq <;> omega
    have hP15 := (hlay5 sP1 (2 * h - 1)).mp hP14
    have hP05 := (hlay5 sP0 (2 * m - 1)).mp hP04
    have hP1Prefix5 := heap5.validU64Slice_mono sP1 (2 * h - 1)
      (2 * m - 1) hP15 (by omega)
    rcases karSubLoop_ok this sP1 sP0 (2 * m - 1) 0 heap5
      hP1Prefix5 hP05 with ⟨heap6, hsub0, hlay6⟩
    have hP1P0short : U64SlicesDisjoint sP1 (2 * m - 1)
        sP0 (2 * m - 1) :=
      u64SlicesDisjoint_mono (u64SlicesDisjoint_symm hP0P1)
        (by omega) (Nat.le_refl _)
    have hpNonzero : this._p ≠ 0 := by
      intro heq
      have hz : this._p.toNat = 0 := congrArg UInt64.toNat heq
      omega
    rcases karSubLoop_refines_full_of_prefix this sP1 sP0 (2 * m - 1)
      tailLength heap5 heap6
      (karPreparedPoly left m h * karPreparedPoly right m h)
      (leftLow * rightLow) hpNonzero (by rw [hP1Length]; exact hP15)
      hP05 hP1P0short (by rw [hP1Length]; exact hCanonicalP15)
      hCanonicalP05 (by rw [hP1Length]; exact hRepP15) hRepP05 hsub0 with
      ⟨_, hRepP1Sub06, hCanonicalP1Sub06⟩
    have hP0Frame6 := karSubLoop_preserves_prefix this sP1 sP0 sP0
      (2 * m - 1) (2 * m - 1) heap5 heap6 hP1Prefix5 hP05
      hP1P0short hsub0
    have hP06 := (hlay6 sP0 (2 * m - 1)).mp hP05
    have hRepP06 := slicePolyRep_of_same_prefix heap5 heap6 sP0
      (2 * m - 1) this._p.toNat (leftLow * rightLow) hP05 hP06
      hP0Frame6 hRepP05
    have hCanonicalP06 := canonicalU64Prefix_of_same_prefix heap5 heap6 sP0
      (2 * m - 1) this._p hP05 hP0Frame6 hCanonicalP05
    have hHigh5 := (hlay5 (C.add (2 * m)) (2 * h - 1)).mp hCHigh4
    have hP1High := u64SlicesDisjoint_symm hHighP1
    have hHighFrame6 := karSubLoop_preserves_prefix this sP1 sP0
      (C.add (2 * m)) (2 * m - 1) (2 * h - 1) heap5 heap6
      hP1Prefix5 hP05
      (u64SlicesDisjoint_mono hP1High (by omega) (Nat.le_refl _)) hsub0
    have hHigh6 := (hlay6 (C.add (2 * m)) (2 * h - 1)).mp hHigh5
    have hRepHigh6 := slicePolyRep_of_same_prefix heap5 heap6
      (C.add (2 * m)) (2 * h - 1) this._p.toNat
      (leftHigh * rightHigh) hHigh5 hHigh6 hHighFrame6 hRepHigh5
    have hCanonicalHigh6 := canonicalU64Prefix_of_same_prefix heap5 heap6
      (C.add (2 * m)) (2 * h - 1) this._p hHigh5 hHighFrame6
      hCanonicalHigh5
    have hP16 := (hlay6 sP1 (2 * h - 1)).mp hP15
    rcases karSubLoop_ok this sP1 (C.add (2 * m)) (2 * h - 1) 0 heap6
      hP16 hHigh6 with ⟨heap7, hsub1, hlay7⟩
    rcases karSubLoop_refines_slice this sP1 (C.add (2 * m))
      (2 * h - 1) heap6 heap7
      (karPreparedPoly left m h * karPreparedPoly right m h -
        leftLow * rightLow) (leftHigh * rightHigh) hpNonzero hP16
      hHigh6 (u64SlicesDisjoint_symm hHighP1)
      (by rw [← hP1Length]; exact hCanonicalP1Sub06)
      hCanonicalHigh6 (by rw [← hP1Length]; exact hRepP1Sub06)
      hRepHigh6 hsub1 with
      ⟨_, hRepCross7, hCanonicalCross7⟩
    have hP1P0Full := u64SlicesDisjoint_symm hP0P1
    have hP0Frame7 := karSubLoop_preserves_prefix this sP1
      (C.add (2 * m)) sP0 (2 * h - 1) (2 * m - 1) heap6 heap7
      hP16 hHigh6 hP1P0Full hsub1
    have hP07 := (hlay7 sP0 (2 * m - 1)).mp hP06
    have hRepP07 := slicePolyRep_of_same_prefix heap6 heap7 sP0
      (2 * m - 1) this._p.toNat (leftLow * rightLow) hP06 hP07
      hP0Frame7 hRepP06
    have hCanonicalP07 := canonicalU64Prefix_of_same_prefix heap6 heap7 sP0
      (2 * m - 1) this._p hP06 hP0Frame7 hCanonicalP06
    have hHighFrame7 := karSubLoop_preserves_prefix this sP1
      (C.add (2 * m)) (C.add (2 * m)) (2 * h - 1) (2 * h - 1)
      heap6 heap7 hP16 hHigh6 (u64SlicesDisjoint_symm hHighP1) hsub1
    have hHigh7 := (hlay7 (C.add (2 * m)) (2 * h - 1)).mp hHigh6
    have hRepHigh7 := slicePolyRep_of_same_prefix heap6 heap7
      (C.add (2 * m)) (2 * h - 1) this._p.toNat
      (leftHigh * rightHigh) hHigh6 hHigh7 hHighFrame7 hRepHigh6
    have hCanonicalHigh7 := canonicalU64Prefix_of_same_prefix heap6 heap7
      (C.add (2 * m)) (2 * h - 1) this._p hHigh6 hHighFrame7
      hCanonicalHigh6
    have hC7 := (hlay7 C (2 * n - 1)).mp ((hlay6 C (2 * n - 1)).mp
      ((hlay5 C (2 * n - 1)).mp ((hlay4 C (2 * n - 1)).mp
        ((hlay3 C (2 * n - 1)).mp ((hlay2 C (2 * n - 1)).mp hC)))))
    have hOutputLength : 2 * m + (2 * h - 1) = 2 * n - 1 := by
      omega
    have hCprefix7 := heap7.validU64Slice_mono C (2 * n - 1)
      (2 * m - 1) hC7 (by omega)
    have hCP0 : U64SlicesDisjoint C (2 * m - 1) sP0 (2 * m - 1) := by
      have hx := u64SlicesDisjoint_subslices hCScratch
        (leftStart := 0) (leftCount := 2 * m - 1)
        (rightStart := 2 * h) (rightCount := 2 * m - 1)
        (by omega) (by
          rw [karScratchNeed_step n (by omega)]
          dsimp [m, h]
          omega)
      simpa [sP0, t2, t1, rawPtr_add_add, two_mul] using hx
    rcases copyU64_ok heap7 C sP0 (2 * m - 1) hCprefix7 hP07 with
      ⟨heap8, hcopy, hlay8⟩
    rcases heap8.writeU64_of_valid C (2 * n - 1) (2 * m - 1) 0
      ((hlay8 C (2 * n - 1)).mp hC7) (by omega) with ⟨heap9, hzero⟩
    have hlay9 := RawHeap.writeU64_sameLayout heap8 heap9 C (2 * m - 1) 0
      hzero
    have hCAsBase : heap7.ValidU64Slice C (2 * m + (2 * h - 1)) := by
      rw [hOutputLength]
      exact hC7
    rcases karCopyZero_refines_base heap7 heap8 heap9 C sP0 m
      (2 * h - 1) this._p.toNat (leftLow * rightLow)
      (leftHigh * rightHigh) this._p hmPos hpNonzero hCAsBase hP07 hCP0
      hRepP07 hCanonicalP07 hRepHigh7 hCanonicalHigh7 hcopy hzero with
      ⟨_, hRepBase9, hCanonicalBase9⟩
    have hCP1Full : U64SlicesDisjoint C (2 * n - 1) sP1 (2 * h - 1) := by
      have hx := u64SlicesDisjoint_add_right hCScratch
        (start := 2 * h + (2 * m - 1)) (count := 2 * h - 1) (by
          rw [karScratchNeed_step n (by omega)]
          dsimp [m, h]
          omega)
      simpa [sP1, sP0, t2, t1, rawPtr_add_add, two_mul,
        Nat.add_assoc] using hx
    have hP17 := (hlay7 sP1 (2 * h - 1)).mp hP16
    have hCopyCross := copyU64_preserves_prefix heap7 heap8 C sP0 sP1
      (2 * m - 1) (2 * h - 1) hCprefix7 hP07
      (u64SlicesDisjoint_mono hCP1Full (by omega) (Nat.le_refl _)) hcopy
    have hP18 := (hlay8 sP1 (2 * h - 1)).mp hP17
    have hRepCross8 := slicePolyRep_of_same_prefix heap7 heap8 sP1
      (2 * h - 1) this._p.toNat
      (karPreparedPoly left m h * karPreparedPoly right m h -
        leftLow * rightLow - leftHigh * rightHigh)
      hP17 hP18 hCopyCross hRepCross7
    have hCanonicalCross8 := canonicalU64Prefix_of_same_prefix heap7 heap8
      sP1 (2 * h - 1) this._p hP17 hCopyCross hCanonicalCross7
    have hZeroCross := writeU64_preserves_prefix heap8 heap9 C sP1
      (2 * n - 1) (2 * h - 1) (2 * m - 1) 0 hCP1Full
      (by omega) hzero
    have hP19 := (hlay9 sP1 (2 * h - 1)).mp hP18
    have hRepCross9 := slicePolyRep_of_same_prefix heap8 heap9 sP1
      (2 * h - 1) this._p.toNat
      (karPreparedPoly left m h * karPreparedPoly right m h -
        leftLow * rightLow - leftHigh * rightHigh)
      hP18 hP19 hZeroCross hRepCross8
    have hCanonicalCross9 := canonicalU64Prefix_of_same_prefix heap8 heap9
      sP1 (2 * h - 1) this._p hP18 hZeroCross hCanonicalCross8
    have hC9 := (hlay9 C (2 * n - 1)).mp ((hlay8 C (2 * n - 1)).mp hC7)
    have hCAssemble := heap9.validU64Slice_mono C (2 * n - 1)
      (m + (2 * h - 1)) hC9 (by dsimp [m, h]; omega)
    rcases karAssembleLoop_ok this C sP1 m (2 * h - 1) 0 heap9
      hCAssemble hP19 with ⟨heap10, hassemble, hlay10⟩
    have hFullLength : (m + (2 * h - 1)) + m = 2 * n - 1 := by
      dsimp [m, h]
      omega
    have hCP1Prefix := u64SlicesDisjoint_mono hCP1Full
      (smallLeft := m + (2 * h - 1)) (by dsimp [m, h]; omega)
      (Nat.le_refl (2 * h - 1))
    have hCanonicalBase9' : CanonicalU64Prefix heap9 C (2 * n - 1)
        this._p := by
      rw [← hOutputLength]
      exact hCanonicalBase9
    have hRepBase9' : SlicePolyRep heap9 C (2 * n - 1) this._p.toNat
        (leftLow * rightLow + Polynomial.X ^ (2 * m) *
          (leftHigh * rightHigh)) := by
      rw [← hOutputLength]
      exact hRepBase9
    rcases karAssembleLoop_refines_full this C sP1 m (2 * h - 1) m
      heap9 heap10
      (leftLow * rightLow + Polynomial.X ^ (2 * m) *
        (leftHigh * rightHigh))
      (karPreparedPoly left m h * karPreparedPoly right m h -
        leftLow * rightLow - leftHigh * rightHigh)
      hpNonzero (by rw [hFullLength]; exact hC9) hP19 hCP1Prefix
      (by rw [hFullLength]; exact hCanonicalBase9') hCanonicalCross9
      (by rw [hFullLength]; exact hRepBase9') hRepCross9 hassemble with
      ⟨_, hRepFinal, hCanonicalFinal⟩
    have hPreparedLeft := karPreparedPoly_eq_low_add_high heap A m h left
      leftLow leftHigh hshape hRepAfull hRepALow hRepAHigh hleftSplit
    have hPreparedRight := karPreparedPoly_eq_low_add_high heap B m h right
      rightLow rightHigh hshape hRepBfull hRepBLow hRepBHigh hrightSplit
    have hKarIdentity := karatsuba_polynomial_identity left right leftLow
      leftHigh rightLow rightHigh m hleftSplit hrightSplit
    have hPolyFinal :
        leftLow * rightLow + Polynomial.X ^ (2 * m) *
            (leftHigh * rightHigh) +
          Polynomial.X ^ m *
            (karPreparedPoly left m h * karPreparedPoly right m h -
              leftLow * rightLow - leftHigh * rightHigh) =
        left * right := by
      rw [hPreparedLeft, hPreparedRight]
      calc
        leftLow * rightLow + Polynomial.X ^ (2 * m) *
              (leftHigh * rightHigh) +
            Polynomial.X ^ m *
              ((leftLow + leftHigh) * (rightLow + rightHigh) -
                leftLow * rightLow - leftHigh * rightHigh) =
            leftLow * rightLow + Polynomial.X ^ m *
              ((leftLow + leftHigh) * (rightLow + rightHigh) -
                leftLow * rightLow - leftHigh * rightHigh) +
              Polynomial.X ^ (2 * m) * (leftHigh * rightHigh) := by ring
        _ = left * right := hKarIdentity
    rw [hPolyFinal] at hRepFinal
    refine ⟨heap10, ?_, ?_, ?_, ?_⟩
    · simp [m, h, t1, t2, sP0, sP1, recScratch, hprep, hp0, hp1,
        phigh, hsub0, hsub1, hcopy, hzero, hassemble]
    · intro ptr length
      exact (hlay2 ptr length).trans ((hlay3 ptr length).trans
        ((hlay4 ptr length).trans ((hlay5 ptr length).trans
          ((hlay6 ptr length).trans ((hlay7 ptr length).trans
            ((hlay8 ptr length).trans ((hlay9 ptr length).trans
              (hlay10 ptr length))))))))
    · rw [← hFullLength]
      exact hRepFinal
    · rw [← hFullLength]
      exact hCanonicalFinal
termination_by n
decreasing_by
  · exact hmLt
  · exact hhLt
  · exact hhLt

end CLPoly.Impl.StrictMulRefinement
