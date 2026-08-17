import CLPoly.Impl.StrictKarMulRefinement
import CLPoly.Impl.StrictEuclidRefinement

set_option autoImplicit false

namespace CLPoly.Impl.StrictMulRefinement

open Generated.StrictMul
open CLPoly.Impl.StrictWordArithmetic
open CLPoly.Impl.StrictDivremRefinement
open CLPoly.Impl.StrictEuclidRefinement
open CLPoly.Impl.RawPolynomialRep

theorem canonicalU64Prefix_mono (heap : RawHeap) (ptr : RawPtr UInt64)
    (large small : Nat) (modulus : UInt64) (hle : small ≤ large)
    (h : CanonicalU64Prefix heap ptr large modulus) :
    CanonicalU64Prefix heap ptr small modulus := by
  intro i value hi hread
  exact h i value (by omega) hread

/-- A full raw representation may be shortened without changing its
polynomial when the omitted coefficients are mathematically zero. -/
theorem slicePolyRep_prefix_of_coeff_zero (heap : RawHeap)
    (ptr : RawPtr UInt64) (length prefixLength p : Nat)
    (poly : Polynomial (ZMod p))
    (hvalid : heap.ValidU64Slice ptr length) (hle : prefixLength ≤ length)
    (hrep : SlicePolyRep heap ptr length p poly)
    (hzero : ∀ degree, prefixLength ≤ degree → poly.coeff degree = 0) :
    SlicePolyRep heap ptr prefixLength p poly := by
  rcases slicePolyRep_prefix_exists heap ptr length prefixLength p poly hvalid
      hle hrep with ⟨prefixPoly, hprefix, hcoeff⟩
  have heq : prefixPoly = poly := by
    ext degree
    by_cases hd : degree < prefixLength
    · exact hcoeff degree hd
    · rw [slicePolyRep_coeff_zero_of_length_le heap ptr prefixLength p
          prefixPoly hprefix degree (by omega), hzero degree (by omega)]
  simpa [heq] using hprefix

/-- The actual zero-padding loop extends a represented prefix without
changing its L2 polynomial.  Recursion follows the generated loop index. -/
theorem mulZeroPadLoop_refines (bPad : RawPtr UInt64)
    (start count i p : Nat) (heap : RawHeap)
    (poly : Polynomial (ZMod p)) (modulus : UInt64)
    (hi : i ≤ count) (hmodulus : modulus ≠ 0)
    (hPad : heap.ValidU64Slice bPad (start + count))
    (hRep : SlicePolyRep heap bPad (start + i) p poly)
    (hCanonical : CanonicalU64Prefix heap bPad (start + i) modulus) :
    ∃ heap', mulZeroPadLoop bPad start count i heap = .ok heap' ∧
      RawHeap.SameLayout heap heap' ∧
      SlicePolyRep heap' bPad (start + count) p poly ∧
      CanonicalU64Prefix heap' bPad (start + count) modulus := by
  rw [mulZeroPadLoop]
  split
  next hmore =>
    have hThroughNext := heap.validU64Slice_mono bPad (start + count)
      ((start + i) + 1) hPad (by omega)
    rcases heap.writeU64_of_valid bPad (start + count) (start + i) 0 hPad
      (by omega) with ⟨heap1, hwrite⟩
    simp only [hwrite]
    rcases writeZero_extends_slice heap heap1 bPad (start + i) p modulus
      poly hmodulus hThroughNext hRep hCanonical hwrite with
      ⟨hlayout1, hRep1, hCanonical1⟩
    have hPad1 := (hlayout1 bPad (start + count)).mp hPad
    rcases mulZeroPadLoop_refines bPad start count (i + 1) p heap1 poly
      modulus (by omega) hmodulus hPad1 (by simpa [Nat.add_assoc] using hRep1)
      (by simpa [Nat.add_assoc] using hCanonical1) with
      ⟨heap2, hrun, hlayout2, hRep2, hCanonical2⟩
    rw [hrun]
    exact ⟨heap2, rfl, fun ptr length =>
      (hlayout1 ptr length).trans (hlayout2 ptr length), hRep2, hCanonical2⟩
  next hdone =>
    have hieq : i = count := by omega
    subst i
    exact ⟨heap, rfl, fun _ _ => Iff.rfl, hRep, hCanonical⟩
termination_by count - i
decreasing_by omega

theorem mulZeroPadLoop_preserves_prefix (bPad guard : RawPtr UInt64)
    (start count i guardLength : Nat) (heap heap' : RawHeap)
    (hPad : heap.ValidU64Slice bPad (start + count))
    (hdisjoint : U64SlicesDisjoint bPad (start + count) guard guardLength)
    (hrun : mulZeroPadLoop bPad start count i heap = .ok heap') :
    SameU64Prefix heap heap' guard guardLength := by
  rw [mulZeroPadLoop] at hrun
  split at hrun
  next hmore =>
    rcases heap.writeU64_of_valid bPad (start + count) (start + i) 0 hPad
      (by omega) with ⟨heap1, hwrite⟩
    simp only [hwrite] at hrun
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 bPad
      (start + i) 0 hwrite
    have hPad1 := (hlayout1 bPad (start + count)).mp hPad
    exact sameU64Prefix_trans
      (writeU64_preserves_prefix heap heap1 bPad guard (start + count)
        guardLength (start + i) 0 hdisjoint (by omega) hwrite)
      (mulZeroPadLoop_preserves_prefix bPad guard start count (i + 1)
        guardLength heap1 heap' hPad1 hdisjoint hrun)
  next hdone =>
    simp only [Except.ok.injEq] at hrun
    subst heap'
    intro _ _ _ hread
    exact hread
termination_by count - i
decreasing_by omega

/-- The zero-fill loop starts at `start`, so every earlier cell of the same
allocation is preserved.  This sharper frame is needed when a C++ buffer is
extended in place. -/
theorem mulZeroPadLoop_preserves_before_start (bPad : RawPtr UInt64)
    (start count i prefixLength : Nat) (heap heap' : RawHeap)
    (hprefix : prefixLength ≤ start)
    (hPad : heap.ValidU64Slice bPad (start + count))
    (hrun : mulZeroPadLoop bPad start count i heap = .ok heap') :
    SameU64Prefix heap heap' bPad prefixLength := by
  rw [mulZeroPadLoop] at hrun
  split at hrun
  next hmore =>
    rcases heap.writeU64_of_valid bPad (start + count) (start + i) 0 hPad
      (by omega) with ⟨heap1, hwrite⟩
    simp only [hwrite] at hrun
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 bPad
      (start + i) 0 hwrite
    have hPad1 := (hlayout1 bPad (start + count)).mp hPad
    exact sameU64Prefix_trans
      (by
        intro k value hk hread
        exact RawHeap.readU64_writeU64_ne heap heap1 bPad bPad
          (start + i) k 0 value hwrite hread (by
            right
            omega))
      (mulZeroPadLoop_preserves_before_start bPad start count (i + 1)
        prefixLength heap1 heap' hprefix hPad1 hrun)
  next hdone =>
    have heq : heap' = heap := Except.ok.inj hrun.symm
    subst heap'
    exact fun _ _ _ hread => hread
termination_by count - i
decreasing_by omega

theorem mulZeroPadLoop_sameLayout (bPad : RawPtr UInt64)
    (start count i : Nat) (heap heap' : RawHeap)
    (hPad : heap.ValidU64Slice bPad (start + count))
    (hrun : mulZeroPadLoop bPad start count i heap = .ok heap') :
    RawHeap.SameLayout heap heap' := by
  rw [mulZeroPadLoop] at hrun
  split at hrun
  next hmore =>
    rcases heap.writeU64_of_valid bPad (start + count) (start + i) 0 hPad
      (by omega) with ⟨heap1, hwrite⟩
    simp only [hwrite] at hrun
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 bPad
      (start + i) 0 hwrite
    have hPad1 := (hlayout1 bPad (start + count)).mp hPad
    have hlayout2 := mulZeroPadLoop_sameLayout bPad start count (i + 1)
      heap1 heap' hPad1 hrun
    exact fun ptr length =>
      (hlayout1 ptr length).trans (hlayout2 ptr length)
  next hdone =>
    simp only [Except.ok.injEq] at hrun
    subst heap'
    exact fun _ _ => Iff.rfl
termination_by count - i
decreasing_by omega

/-- The complete generated `_mul` dispatcher only writes C and scratch. -/
theorem mul_preserves_prefix (this : DenseUPolyZp)
    (C A : RawPtr UInt64) (lenA : Nat) (B : RawPtr UInt64) (lenB : Nat)
    (scratch guard : RawPtr UInt64) (guardLen : Nat)
    (heap heap' : RawHeap)
    (hApos : 0 < lenA) (hBpos : 0 < lenB) (hBA : lenB ≤ lenA)
    (hC : heap.ValidU64Slice C (2 * lenA - 1))
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB)
    (hScratch : heap.ValidU64Slice scratch (8 * lenA))
    (hScratchB : U64SlicesDisjoint scratch (8 * lenA) B lenB)
    (hCGuard : U64SlicesDisjoint C (2 * lenA - 1) guard guardLen)
    (hScratchGuard : U64SlicesDisjoint scratch (8 * lenA) guard guardLen)
    (hrun : dense_upoly_zp__mul_ir this C A lenA B lenB scratch heap =
      .ok heap') :
    SameU64Prefix heap heap' guard guardLen := by
  have hassert : lenB ≤ lenA ∧ 0 < lenB := ⟨hBA, hBpos⟩
  by_cases hschool : lenB < 16
  · have hrunSchool : dense_upoly_zp__classical_mul_ir this C A lenA B lenB
        heap = .ok heap' := by
      simpa [dense_upoly_zp__mul_ir, hassert, hschool] using hrun
    exact classicalMul_preserves_outside this C A lenA B lenB guard guardLen
      heap heap' hApos hBpos
      (heap.validU64Slice_mono C (2 * lenA - 1) (lenA + lenB - 1) hC
        (by omega)) hA hB
      (u64SlicesDisjoint_mono hCGuard (by omega) (by omega)) hrunSchool
  · let bPad := scratch
    let karScratch := scratch.add lenA
    cases hcopy : heap.copyU64 bPad B lenB with
    | error fault =>
        simp [dense_upoly_zp__mul_ir, hassert, hschool, bPad, karScratch,
          hcopy] at hrun
    | ok heap1 =>
        cases hzero : mulZeroPadLoop bPad lenB (lenA - lenB) 0 heap1 with
        | error fault =>
            simp [dense_upoly_zp__mul_ir, hassert, hschool, bPad, karScratch,
              hcopy, hzero] at hrun
        | ok heap2 =>
            have hkar : dense_upoly_zp__kar_mul_ir this C A bPad lenA
                karScratch heap2 = .ok heap' := by
              simpa [dense_upoly_zp__mul_ir, hassert, hschool, bPad,
                karScratch, hcopy, hzero] using hrun
            have hBPad : heap.ValidU64Slice bPad lenA := by
              dsimp [bPad]
              exact heap.validU64Slice_mono scratch (8 * lenA) lenA hScratch
                (by omega)
            have hBPadPrefix := heap.validU64Slice_mono bPad lenA lenB hBPad hBA
            have hPadB : U64SlicesDisjoint bPad lenB B lenB := by
              dsimp [bPad]
              exact u64SlicesDisjoint_mono hScratchB (by omega) (by omega)
            rcases copyU64_refines_disjoint heap bPad B lenB hBPadPrefix hB
              hPadB with ⟨copyHeap, hcopy', hlayout1, _⟩
            have heq1 : copyHeap = heap1 := Except.ok.inj (hcopy'.symm.trans hcopy)
            subst copyHeap
            have hBPad1 := (hlayout1 bPad lenA).mp hBPad
            have hlayout2 := mulZeroPadLoop_sameLayout bPad lenB
              (lenA - lenB) 0 heap1 heap2
              (by simpa [Nat.add_sub_of_le hBA] using hBPad1) hzero
            have hcopyGuard := copyU64_preserves_prefix heap heap1 bPad B guard
              lenB guardLen hBPadPrefix hB
              (by
                dsimp [bPad]
                exact u64SlicesDisjoint_mono hScratchGuard (by omega) (by omega))
              hcopy
            have hzeroGuard := mulZeroPadLoop_preserves_prefix bPad guard lenB
              (lenA - lenB) 0 guardLen heap1 heap2
              (by simpa [Nat.add_sub_of_le hBA] using hBPad1)
              (by
                dsimp [bPad]
                simpa [Nat.add_sub_of_le hBA] using
                  u64SlicesDisjoint_mono hScratchGuard
                    (smallLeft := lenA) (smallRight := guardLen)
                    (by omega) (by omega)) hzero
            have hNeed := karScratchNeed_le_seven lenA
            have hKarScratch : heap.ValidU64Slice karScratch
                (karScratchNeed lenA) := by
              dsimp [karScratch]
              exact heap.validU64Slice_add scratch (8 * lenA) lenA
                (karScratchNeed lenA) hScratch (by omega)
            have hkarGuard := karMul_preserves_prefix this C A bPad lenA
              karScratch guard guardLen heap2 heap' hApos
              ((hlayout2 C (2 * lenA - 1)).mp
                ((hlayout1 C (2 * lenA - 1)).mp hC))
              ((hlayout2 A lenA).mp ((hlayout1 A lenA).mp hA))
              ((hlayout2 bPad lenA).mp hBPad1)
              ((hlayout2 karScratch (karScratchNeed lenA)).mp
                ((hlayout1 karScratch (karScratchNeed lenA)).mp hKarScratch))
              hCGuard
              (by
                dsimp [karScratch]
                exact u64SlicesDisjoint_add_left hScratchGuard (by
                  calc
                    lenA + karScratchNeed lenA ≤ lenA + 7 * lenA :=
                      Nat.add_le_add_left hNeed lenA
                    _ = 8 * lenA := by omega)) hkar
            exact sameU64Prefix_trans hcopyGuard
              (sameU64Prefix_trans hzeroGuard hkarGuard)

/-- Semantic refinement of the actual generated `_mul` dispatcher.  The
Karatsuba branch copies and pads B in the C++ scratch allocation before
calling the already-refined well-founded `_kar_mul`; the larger symmetric
Karatsuba result is then shortened only after proving its omitted product
coefficients are zero. -/
theorem mul_refines_slice (this : DenseUPolyZp)
    (C A : RawPtr UInt64) (lenA : Nat) (B : RawPtr UInt64) (lenB : Nat)
    (scratch : RawPtr UInt64) (heap : RawHeap)
    (left right : Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (hApos : 0 < lenA) (hBpos : 0 < lenB) (hBA : lenB ≤ lenA)
    (hLenAWord : lenA < limbBase)
    (hC : heap.ValidU64Slice C (2 * lenA - 1))
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB)
    (hScratch : heap.ValidU64Slice scratch (8 * lenA))
    (hCA : U64SlicesDisjoint C (2 * lenA - 1) A lenA)
    (hCB : U64SlicesDisjoint C (2 * lenA - 1) B lenB)
    (hCScratch : U64SlicesDisjoint C (2 * lenA - 1) scratch (8 * lenA))
    (hScratchA : U64SlicesDisjoint scratch (8 * lenA) A lenA)
    (hScratchB : U64SlicesDisjoint scratch (8 * lenA) B lenB)
    (hCanonicalA : CanonicalU64Prefix heap A lenA this._p)
    (hCanonicalB : CanonicalU64Prefix heap B lenB this._p)
    (hRepA : SlicePolyRep heap A lenA this._p.toNat left)
    (hRepB : SlicePolyRep heap B lenB this._p.toNat right) :
    ∃ heap', dense_upoly_zp__mul_ir this C A lenA B lenB scratch heap =
        .ok heap' ∧ RawHeap.SameLayout heap heap' ∧
      SlicePolyRep heap' C (lenA + lenB - 1) this._p.toNat
        (left * right) ∧
      CanonicalU64Prefix heap' C (lenA + lenB - 1) this._p := by
  rw [dense_upoly_zp__mul_ir]
  rw [if_neg (by
    intro hnot
    exact hnot ⟨hBA, hBpos⟩)]
  by_cases hschool : lenB < 16
  · rw [if_pos hschool]
    have hCshort := heap.validU64Slice_mono C (2 * lenA - 1)
      (lenA + lenB - 1) hC (by omega)
    have hCAshort := u64SlicesDisjoint_mono hCA
      (smallLeft := lenA + lenB - 1) (smallRight := lenA)
      (by omega) (by omega)
    have hCBshort := u64SlicesDisjoint_mono hCB
      (smallLeft := lenA + lenB - 1) (smallRight := lenB)
      (by omega) (by omega)
    exact classicalMul_refines_slice this C A lenA B lenB heap left right
      hcfg hp hApos hBpos hLenAWord hCshort hA hB hCAshort hCBshort
      hCanonicalA hCanonicalB hRepA hRepB
  · rw [if_neg hschool]
    let bPad := scratch
    let karScratch := scratch.add lenA
    have hNeed : karScratchNeed lenA ≤ 7 * lenA :=
      karScratchNeed_le_seven lenA
    have hBPad : heap.ValidU64Slice bPad lenA := by
      dsimp [bPad]
      exact heap.validU64Slice_mono scratch (8 * lenA) lenA hScratch (by omega)
    have hKarScratch : heap.ValidU64Slice karScratch
        (karScratchNeed lenA) := by
      dsimp [karScratch]
      exact heap.validU64Slice_add scratch (8 * lenA) lenA
        (karScratchNeed lenA) hScratch (by omega)
    have hBPadB : U64SlicesDisjoint bPad lenB B lenB := by
      dsimp [bPad]
      exact u64SlicesDisjoint_mono hScratchB (by omega) (by omega)
    rcases copyU64_refines_slice_canonical heap bPad B lenB
      this._p.toNat right this._p
      (heap.validU64Slice_mono bPad lenA lenB hBPad hBA) hB hBPadB
      hRepB hCanonicalB with
      ⟨heap1, hcopy, hlayout1, hRepPad1, hCanonicalPad1⟩
    have hBPad1 := (hlayout1 bPad lenA).mp hBPad
    have hPadA : U64SlicesDisjoint bPad lenA A lenA := by
      dsimp [bPad]
      exact u64SlicesDisjoint_mono hScratchA (by omega) (by omega)
    have hcopyA := copyU64_preserves_prefix heap heap1 bPad B A lenB lenA
      (heap.validU64Slice_mono bPad lenA lenB hBPad hBA) hB
      (u64SlicesDisjoint_mono hPadA hBA (by omega)) hcopy
    have hA1 := (hlayout1 A lenA).mp hA
    have hRepA1 := slicePolyRep_of_same_prefix heap heap1 A lenA
      this._p.toNat left hA hA1 hcopyA hRepA
    have hCanonicalA1 := canonicalU64Prefix_of_same_prefix heap heap1 A lenA
      this._p hA hcopyA hCanonicalA
    have hpWord : this._p ≠ 0 := by
      intro hzero
      have hzeroNat := congrArg UInt64.toNat hzero
      simp at hzeroNat
      omega
    rcases mulZeroPadLoop_refines bPad lenB (lenA - lenB) 0
      this._p.toNat heap1 right this._p (by omega) hpWord
      (by simpa [Nat.add_sub_of_le hBA] using hBPad1)
      (by simpa using hRepPad1) (by simpa using hCanonicalPad1) with
      ⟨heap2, hzero, hlayout2, hRepPad2, hCanonicalPad2⟩
    have hzeroA := mulZeroPadLoop_preserves_prefix bPad A lenB
      (lenA - lenB) 0 lenA heap1 heap2
      (by simpa [Nat.add_sub_of_le hBA] using hBPad1)
      (by simpa [Nat.add_sub_of_le hBA] using hPadA) hzero
    have hA2 := (hlayout2 A lenA).mp hA1
    have hRepA2 := slicePolyRep_of_same_prefix heap1 heap2 A lenA
      this._p.toNat left hA1 hA2 hzeroA hRepA1
    have hCanonicalA2 := canonicalU64Prefix_of_same_prefix heap1 heap2 A
      lenA this._p hA1 hzeroA hCanonicalA1
    have hC2 := (hlayout2 C (2 * lenA - 1)).mp
      ((hlayout1 C (2 * lenA - 1)).mp hC)
    have hKarScratch2 := (hlayout2 karScratch (karScratchNeed lenA)).mp
      ((hlayout1 karScratch (karScratchNeed lenA)).mp hKarScratch)
    have hCKar : U64SlicesDisjoint C (2 * lenA - 1) karScratch
        (karScratchNeed lenA) := by
      dsimp [karScratch]
      exact u64SlicesDisjoint_add_right hCScratch (by omega)
    have hKarA : U64SlicesDisjoint karScratch (karScratchNeed lenA) A lenA := by
      dsimp [karScratch]
      exact u64SlicesDisjoint_add_left
        (start := lenA) (count := karScratchNeed lenA) hScratchA (by
          calc
            lenA + karScratchNeed lenA ≤ lenA + 7 * lenA :=
              Nat.add_le_add_left hNeed lenA
            _ = 8 * lenA := by omega)
    have hKarPad : U64SlicesDisjoint karScratch (karScratchNeed lenA)
        bPad lenA := by
      dsimp [karScratch, bPad]
      exact u64SlicesDisjoint_symm
        (u64SlicesDisjoint_adjacent scratch lenA (karScratchNeed lenA))
    have hCPad : U64SlicesDisjoint C (2 * lenA - 1) bPad lenA := by
      dsimp [bPad]
      exact u64SlicesDisjoint_mono hCScratch (by omega) (by omega)
    rcases karMul_refines_slice this C A bPad lenA karScratch heap2 left
      right hcfg hp hApos hLenAWord hC2 hA2
      ((hlayout2 bPad lenA).mp hBPad1) hKarScratch2 hCA hCPad hCKar
      hKarA hKarPad hCanonicalA2
      (by simpa [Nat.add_sub_of_le hBA] using hCanonicalPad2)
      hRepA2 (by simpa [Nat.add_sub_of_le hBA] using hRepPad2) with
      ⟨heap3, hmul, hlayout3, hRepFull, hCanonicalFull⟩
    have hlogical : lenA + lenB - 1 ≤ 2 * lenA - 1 := by omega
    have hC3 := (hlayout3 C (2 * lenA - 1)).mp hC2
    have hRepShort := slicePolyRep_prefix_of_coeff_zero heap3 C
      (2 * lenA - 1) (lenA + lenB - 1) this._p.toNat (left * right)
      hC3 hlogical hRepFull (by
        intro degree hdegree
        exact mul_coeff_zero_of_slice_lengths heap A B lenA lenB degree
          left right hApos hBpos hRepA hRepB hdegree)
    refine ⟨heap3, ?_, fun ptr length =>
      (hlayout1 ptr length).trans ((hlayout2 ptr length).trans
        (hlayout3 ptr length)), hRepShort,
      canonicalU64Prefix_mono heap3 C (2 * lenA - 1)
        (lenA + lenB - 1) this._p hlogical hCanonicalFull⟩
    simpa [bPad, karScratch, hcopy, hzero] using hmul

/-- Normalized-object form of the real `_mul` dispatcher, needed by HGCD
matrix arithmetic.  Nonzero leading output follows from the two normalized
input buffers and the actual final-cell read. -/
theorem mul_refines_rawDense (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (C A : RawPtr UInt64) (lenA : Nat) (B : RawPtr UInt64) (lenB : Nat)
    (scratch : RawPtr UInt64) (heap : RawHeap)
    (left right : Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (hApos : 0 < lenA) (hBpos : 0 < lenB) (hBA : lenB ≤ lenA)
    (hLenAWord : lenA < limbBase)
    (hC : heap.ValidU64Slice C (2 * lenA - 1))
    (hScratch : heap.ValidU64Slice scratch (8 * lenA))
    (hCA : U64SlicesDisjoint C (2 * lenA - 1) A lenA)
    (hCB : U64SlicesDisjoint C (2 * lenA - 1) B lenB)
    (hCScratch : U64SlicesDisjoint C (2 * lenA - 1) scratch (8 * lenA))
    (hScratchA : U64SlicesDisjoint scratch (8 * lenA) A lenA)
    (hScratchB : U64SlicesDisjoint scratch (8 * lenA) B lenB)
    (hLeft : RawDensePolyRep this heap A lenA left)
    (hRight : RawDensePolyRep this heap B lenB right) :
    ∃ heap', dense_upoly_zp__mul_ir this C A lenA B lenB scratch heap =
        .ok heap' ∧ RawHeap.SameLayout heap heap' ∧
      RawDensePolyRep this heap' C (lenA + lenB - 1) (left * right) := by
  rcases mul_refines_slice this C A lenA B lenB scratch heap left right hcfg
      hp hApos hBpos hBA hLenAWord hC hLeft.1 hRight.1 hScratch hCA hCB
      hCScratch hScratchA hScratchB hLeft.2.1 hRight.2.1 hLeft.2.2.1
      hRight.2.2.1 with
    ⟨heap', hrun, hlayout, hrep, hcanonical⟩
  let outLen := lenA + lenB - 1
  have houtPos : 0 < outLen := by dsimp [outLen]; omega
  have hlogical : outLen ≤ 2 * lenA - 1 := by dsimp [outLen]; omega
  have hvalid := (hlayout C outLen).mp
    (heap.validU64Slice_mono C (2 * lenA - 1) outLen hC hlogical)
  have hlast := mul_last_coeff_ne_zero_of_rawDense this heap A B lenA lenB
    left right hApos hBpos hLeft hRight
  rcases slicePolyRep_coeff heap' C outLen this._p.toNat (left * right) hrep
      (outLen - 1) (by omega) with ⟨value, hread, hcoeff⟩
  have hvalue : value ≠ 0 := by
    intro hz
    subst value
    apply hlast
    simpa [outLen] using hcoeff
  have hnorm : heap'.normaliseU64 C outLen = .ok outLen := by
    rw [show outLen = (outLen - 1) + 1 by omega]
    simp [RawHeap.normaliseU64, hread, hvalue]
  exact ⟨heap', hrun, hlayout, hvalid, hcanonical, hrep, hnorm⟩

end CLPoly.Impl.StrictMulRefinement
