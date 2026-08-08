import CLPoly.Generated.StrictMul
import CLPoly.Impl.StrictDivremRefinement

set_option autoImplicit false

namespace CLPoly.Impl.StrictMulRefinement

open Generated.StrictMul
open CLPoly.Impl.StrictWordArithmetic
open CLPoly.Impl.StrictDivremRefinement
open CLPoly.Impl.RawPolynomialRep

/-- Mathematical value of the exact raw cells visited by the C++ dot loop.
It shares the loop's reads and failure behavior, but performs unbounded natural
addition so that machine accumulation can be related to it explicitly. -/
def classicalDotNat (heap : RawHeap) (A B : RawPtr UInt64)
    (k stop j : Nat) : RawExec Nat :=
  if h : j ≤ stop then
    match heap.readU64 A j with
    | .error fault => .error fault
    | .ok a =>
      match heap.readU64 B (k - j) with
      | .error fault => .error fault
      | .ok b =>
        match classicalDotNat heap A B k stop (j + 1) with
        | .error fault => .error fault
        | .ok tail => .ok (a.toNat * b.toNat + tail)
  else
    .ok 0
termination_by stop + 1 - j
decreasing_by omega

def classicalDotPoly {p : Nat} (left right : Polynomial (ZMod p))
    (k stop j : Nat) : ZMod p :=
  if h : j ≤ stop then
    left.coeff j * right.coeff (k - j) +
      classicalDotPoly left right k stop (j + 1)
  else
    0
termination_by stop + 1 - j
decreasing_by omega

theorem classicalDotPoly_eq_sum_Icc {p : Nat}
    (left right : Polynomial (ZMod p)) (k stop j : Nat) :
    classicalDotPoly left right k stop j =
      ∑ t ∈ Finset.Icc j stop, left.coeff t * right.coeff (k - t) := by
  rw [classicalDotPoly]
  split
  next hle =>
    rw [← Finset.insert_Icc_succ_left_eq_Icc hle]
    simp
    exact classicalDotPoly_eq_sum_Icc left right k stop (j + 1)
  next hnot =>
    symm
    apply Finset.sum_eq_zero
    intro t ht
    simp only [Finset.mem_Icc] at ht
    omega
termination_by stop + 1 - j
decreasing_by omega

theorem classicalDotPoly_source_eq_coeff {p : Nat}
    (heap : RawHeap) (A B : RawPtr UInt64) (lenA lenB k : Nat)
    (left right : Polynomial (ZMod p))
    (hLenA : 0 < lenA)
    (hRepA : SlicePolyRep heap A lenA p left)
    (hRepB : SlicePolyRep heap B lenB p right) :
    classicalDotPoly left right k
      (if k < lenA then k else lenA - 1)
      (if k ≥ lenB then k - lenB + 1 else 0) =
        (left * right).coeff k := by
  rw [classicalDotPoly_eq_sum_Icc]
  rw [Polynomial.coeff_mul,
    Finset.Nat.sum_antidiagonal_eq_sum_range_succ_mk]
  by_cases hkA : k < lenA <;> by_cases hkB : k ≥ lenB
  all_goals simp only [hkA, hkB, if_true, if_false]
  all_goals
  apply Finset.sum_subset
  · intro t ht
    simp only [Finset.mem_Icc] at ht
    simp only [Finset.mem_range]
    omega
  · intro t htRange htIcc
    simp only [Finset.mem_range] at htRange
    simp only [Finset.mem_Icc, not_and_or, not_le] at htIcc
    rcases htIcc with hBelow | hAbove
    · by_cases hkt : lenB ≤ k
      · have hzero := slicePolyRep_coeff_zero_of_length_le heap B lenB p
          right hRepB (k - t) (by omega)
        rw [hzero, mul_zero]
      · omega
    · by_cases hkt : k < lenA
      · omega
      · have hzero := slicePolyRep_coeff_zero_of_length_le heap A lenA p
          left hRepA t (by omega)
        rw [hzero, zero_mul]

theorem classicalDotNat_cast_eq_poly (heap : RawHeap)
    (A B : RawPtr UInt64) (lenA lenB p k stop j sum : Nat)
    (left right : Polynomial (ZMod p))
    (hRepA : SlicePolyRep heap A lenA p left)
    (hRepB : SlicePolyRep heap B lenB p right)
    (hAIndex : ∀ t, j ≤ t → t ≤ stop → t < lenA)
    (hBIndex : ∀ t, j ≤ t → t ≤ stop → k - t < lenB)
    (hrun : classicalDotNat heap A B k stop j = .ok sum) :
    (sum : ZMod p) = classicalDotPoly left right k stop j := by
  unfold classicalDotNat at hrun
  split at hrun
  next hle =>
    rw [classicalDotPoly, dif_pos hle]
    have hjA := hAIndex j (Nat.le_refl _) hle
    have hjB := hBIndex j (Nat.le_refl _) hle
    rcases slicePolyRep_coeff heap A lenA p left hRepA j hjA with
      ⟨a, ha, hcoeffA⟩
    simp only [ha] at hrun
    rcases slicePolyRep_coeff heap B lenB p right hRepB (k - j) hjB with
      ⟨b, hb, hcoeffB⟩
    simp only [hb] at hrun
    cases ht : classicalDotNat heap A B k stop (j + 1) with
    | error fault => simp [ht] at hrun
    | ok tail =>
      simp only [ht] at hrun
      have hsum : sum = a.toNat * b.toNat + tail :=
        Except.ok.inj hrun.symm
      have htail := classicalDotNat_cast_eq_poly heap A B lenA lenB p k
        stop (j + 1) tail left right hRepA hRepB
        (by intro t hjt hts; exact hAIndex t (by omega) hts)
        (by intro t hjt hts; exact hBIndex t (by omega) hts) ht
      rw [hsum, Nat.cast_add, Nat.cast_mul, hcoeffA, hcoeffB, htail]
  next hnot =>
    rw [classicalDotPoly, dif_neg hnot]
    have hzero : sum = 0 := Except.ok.inj hrun.symm
    subst sum
    simp
termination_by stop + 1 - j
decreasing_by omega

theorem classicalDotLoop_modEq (heap : RawHeap) (A B : RawPtr UInt64)
    (k stop j : Nat) (acc result : Word3) (sum : Nat)
    (hrun : classicalDotLoop heap A B k stop j acc = .ok result)
    (hsum : classicalDotNat heap A B k stop j = .ok sum) :
    Nat.ModEq (limbBase ^ 3) (word3Value result)
      (word3Value acc + sum) := by
  unfold classicalDotLoop at hrun
  unfold classicalDotNat at hsum
  split at hrun
  next hle =>
    simp only [hle, ↓reduceDIte] at hsum
    cases ha : heap.readU64 A j with
    | error fault => simp [ha] at hrun
    | ok a =>
      simp only [ha] at hrun hsum
      cases hb : heap.readU64 B (k - j) with
      | error fault => simp [hb] at hrun
      | ok b =>
        simp only [hb] at hrun hsum
        cases ht : classicalDotNat heap A B k stop (j + 1) with
        | error fault => simp [ht] at hsum
        | ok tail =>
          simp only [ht] at hsum
          have hsumEq : sum = a.toNat * b.toNat + tail :=
            Except.ok.inj hsum.symm
          let product := Generated.StrictGCD.dense_upoly_zp__umul128_ir 0 0 a b
          let acc' := Generated.StrictGCD.dense_upoly_zp__add_carry3_ir
            acc product.fst product.snd
          have hrec := classicalDotLoop_modEq heap A B k stop (j + 1)
            acc' result tail hrun ht
          have hstep : Nat.ModEq (limbBase ^ 3) (word3Value acc')
              (word3Value acc + a.toNat * b.toNat) := by
            simpa [product, acc'] using addMulWord3_modEq acc a b
          have htotal := hrec.trans (hstep.add_right tail)
          simpa [hsumEq, Nat.add_assoc] using htotal
  next hnot =>
    simp only [hnot, ↓reduceDIte] at hsum
    have hresult : result = acc := Except.ok.inj hrun.symm
    have hzero : sum = 0 := Except.ok.inj hsum.symm
    subst result
    subst sum
    simpa using (Nat.ModEq.refl (word3Value acc) :
      Nat.ModEq (limbBase ^ 3) (word3Value acc) (word3Value acc))
termination_by stop + 1 - j
decreasing_by omega

theorem classicalDotNat_ok (heap : RawHeap) (A B : RawPtr UInt64)
    (lenA lenB k stop j : Nat)
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB)
    (hAIndex : ∀ t, j ≤ t → t ≤ stop → t < lenA)
    (hBIndex : ∀ t, j ≤ t → t ≤ stop → k - t < lenB) :
    ∃ sum, classicalDotNat heap A B k stop j = .ok sum := by
  unfold classicalDotNat
  split
  next hle =>
    rcases heap.readU64_of_valid A lenA j hA
      (hAIndex j (Nat.le_refl _) hle) with ⟨a, ha⟩
    simp only [ha]
    rcases heap.readU64_of_valid B lenB (k - j) hB
      (hBIndex j (Nat.le_refl _) hle) with ⟨b, hb⟩
    simp only [hb]
    rcases classicalDotNat_ok heap A B lenA lenB k stop (j + 1)
      hA hB
      (by intro t hjt hts; exact hAIndex t (by omega) hts)
      (by intro t hjt hts; exact hBIndex t (by omega) hts) with
      ⟨tail, htail⟩
    rw [htail]
    exact ⟨a.toNat * b.toNat + tail, rfl⟩
  next => exact ⟨0, rfl⟩
termination_by stop + 1 - j
decreasing_by omega

theorem classicalDotNat_bound (heap : RawHeap) (A B : RawPtr UInt64)
    (lenA lenB k stop j : Nat) (p : UInt64) (sum : Nat)
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB)
    (hCanonicalA : CanonicalU64Prefix heap A lenA p)
    (hCanonicalB : CanonicalU64Prefix heap B lenB p)
    (hAIndex : ∀ t, j ≤ t → t ≤ stop → t < lenA)
    (hBIndex : ∀ t, j ≤ t → t ≤ stop → k - t < lenB)
    (hrun : classicalDotNat heap A B k stop j = .ok sum) :
    sum ≤ (stop + 1 - j) * (p.toNat - 1) ^ 2 := by
  unfold classicalDotNat at hrun
  split at hrun
  next hle =>
    have hjA := hAIndex j (Nat.le_refl _) hle
    have hjB := hBIndex j (Nat.le_refl _) hle
    rcases heap.readU64_of_valid A lenA j hA hjA with ⟨a, ha⟩
    simp only [ha] at hrun
    rcases heap.readU64_of_valid B lenB (k - j) hB hjB with ⟨b, hb⟩
    simp only [hb] at hrun
    cases ht : classicalDotNat heap A B k stop (j + 1) with
    | error fault => simp [ht] at hrun
    | ok tail =>
      simp only [ht] at hrun
      have hsum : sum = a.toNat * b.toNat + tail :=
        Except.ok.inj hrun.symm
      have haLe : a.toNat ≤ p.toNat - 1 := by
        have := hCanonicalA j a hjA ha
        omega
      have hbLe : b.toNat ≤ p.toNat - 1 := by
        have := hCanonicalB (k - j) b hjB hb
        omega
      have hprod : a.toNat * b.toNat ≤ (p.toNat - 1) ^ 2 := by
        rw [pow_two]
        exact Nat.mul_le_mul haLe hbLe
      have htail := classicalDotNat_bound heap A B lenA lenB k stop (j + 1)
        p tail hA hB hCanonicalA hCanonicalB
        (by intro t hjt hts; exact hAIndex t (by omega) hts)
        (by intro t hjt hts; exact hBIndex t (by omega) hts) ht
      rw [hsum]
      calc
        a.toNat * b.toNat + tail ≤
            (p.toNat - 1) ^ 2 +
              (stop + 1 - (j + 1)) * (p.toNat - 1) ^ 2 :=
          Nat.add_le_add hprod htail
        _ = (stop + 1 - j) * (p.toNat - 1) ^ 2 := by
          have : stop + 1 - j = (stop + 1 - (j + 1)) + 1 := by omega
          rw [this]
          ring
  next hnot =>
    have hzero : sum = 0 := Except.ok.inj hrun.symm
    subst sum
    exact Nat.zero_le _
termination_by stop + 1 - j
decreasing_by omega

theorem classicalDotLoop_ok (heap : RawHeap) (A B : RawPtr UInt64)
    (lenA lenB k stop j : Nat) (acc : Word3)
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB)
    (hAIndex : ∀ t, j ≤ t → t ≤ stop → t < lenA)
    (hBIndex : ∀ t, j ≤ t → t ≤ stop → k - t < lenB) :
    ∃ result, classicalDotLoop heap A B k stop j acc = .ok result := by
  unfold classicalDotLoop
  split
  next hle =>
    have hjA := hAIndex j (Nat.le_refl _) hle
    have hjB := hBIndex j (Nat.le_refl _) hle
    rcases heap.readU64_of_valid A lenA j hA hjA with ⟨a, ha⟩
    simp only [ha]
    rcases heap.readU64_of_valid B lenB (k - j) hB hjB with ⟨b, hb⟩
    simp only [hb]
    let product := Generated.StrictGCD.dense_upoly_zp__umul128_ir 0 0 a b
    let acc' := Generated.StrictGCD.dense_upoly_zp__add_carry3_ir
      acc product.fst product.snd
    apply classicalDotLoop_ok heap A B lenA lenB k stop (j + 1) acc'
      hA hB
    · intro t hjt hts
      exact hAIndex t (by omega) hts
    · intro t hjt hts
      exact hBIndex t (by omega) hts
  next => exact ⟨acc, rfl⟩
termination_by stop + 1 - j
decreasing_by omega

theorem classicalDotLoop_raw_sum (heap : RawHeap) (A B : RawPtr UInt64)
    (lenA lenB k stop j : Nat) (acc : Word3)
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB)
    (hAIndex : ∀ t, j ≤ t → t ≤ stop → t < lenA)
    (hBIndex : ∀ t, j ≤ t → t ≤ stop → k - t < lenB) :
    ∃ result sum,
      classicalDotLoop heap A B k stop j acc = .ok result ∧
      classicalDotNat heap A B k stop j = .ok sum ∧
      Nat.ModEq (limbBase ^ 3) (word3Value result)
        (word3Value acc + sum) := by
  rcases classicalDotLoop_ok heap A B lenA lenB k stop j acc hA hB
    hAIndex hBIndex with ⟨result, hrun⟩
  rcases classicalDotNat_ok heap A B lenA lenB k stop j hA hB
    hAIndex hBIndex with ⟨sum, hsum⟩
  exact ⟨result, sum, hrun, hsum,
    classicalDotLoop_modEq heap A B k stop j acc result sum hrun hsum⟩

theorem classicalDotLoop_exact_zero (heap : RawHeap) (A B : RawPtr UInt64)
    (lenA lenB k stop j : Nat) (p : UInt64) (result : Word3) (sum : Nat)
    (hp : 1 < p.toNat)
    (hcount : stop + 1 - j < limbBase)
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB)
    (hCanonicalA : CanonicalU64Prefix heap A lenA p)
    (hCanonicalB : CanonicalU64Prefix heap B lenB p)
    (hAIndex : ∀ t, j ≤ t → t ≤ stop → t < lenA)
    (hBIndex : ∀ t, j ≤ t → t ≤ stop → k - t < lenB)
    (hrun : classicalDotLoop heap A B k stop j
      { lo := 0, mid := 0, hi := 0 } = .ok result)
    (hsum : classicalDotNat heap A B k stop j = .ok sum) :
    word3Value result = sum := by
  have hbound := classicalDotNat_bound heap A B lenA lenB k stop j p sum
    hA hB hCanonicalA hCanonicalB hAIndex hBIndex hsum
  have hlazy := lazyAccumulation_word3_budget p (stop + 1 - j) 0
    hp hcount (by omega)
  have hpB : p.toNat < limbBase := by
    simpa [limbBase] using UInt64.toNat_lt p
  have hsumLt : sum < limbBase ^ 3 := by
    calc
      sum ≤ (stop + 1 - j) * (p.toNat - 1) ^ 2 := hbound
      _ < p.toNat * limbBase ^ 2 := by simpa using hlazy
      _ < limbBase ^ 3 := by
        have hpowPos : 0 < limbBase ^ 2 :=
          pow_pos (by norm_num [limbBase]) 2
        have hmul : p.toNat * limbBase ^ 2 < limbBase * limbBase ^ 2 :=
          Nat.mul_lt_mul_of_pos_right hpB hpowPos
        simpa [pow_succ, Nat.mul_comm, Nat.mul_left_comm,
          Nat.mul_assoc] using hmul
  have hmod := classicalDotLoop_modEq heap A B k stop j
    { lo := 0, mid := 0, hi := 0 } result sum hrun hsum
  have hmod' : Nat.ModEq (limbBase ^ 3) (word3Value result) sum := by
    simpa [word3Value] using hmod
  exact hmod'.eq_of_lt_of_lt (word3Value_lt result) hsumLt

theorem classicalDotReduced_toNat (this : DenseUPolyZp)
    (heap : RawHeap) (A B : RawPtr UInt64)
    (lenA lenB k stop j : Nat) (result : Word3) (sum : Nat)
    (hcfg : DensePreinvConfigured this)
    (hp : 1 < this._p.toNat)
    (hcount : stop + 1 - j < limbBase)
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB)
    (hCanonicalA : CanonicalU64Prefix heap A lenA this._p)
    (hCanonicalB : CanonicalU64Prefix heap B lenB this._p)
    (hAIndex : ∀ t, j ≤ t → t ≤ stop → t < lenA)
    (hBIndex : ∀ t, j ≤ t → t ≤ stop → k - t < lenB)
    (hrun : classicalDotLoop heap A B k stop j
      { lo := 0, mid := 0, hi := 0 } = .ok result)
    (hsum : classicalDotNat heap A B k stop j = .ok sum) :
    (Generated.StrictGCD.dense_upoly_zp__lll_mod_preinv_ir
      result.hi result.mid result.lo this._p this._ninv this._norm).toNat =
      sum % this._p.toNat := by
  have hexact := classicalDotLoop_exact_zero heap A B lenA lenB k stop j
    this._p result sum hp hcount hA hB hCanonicalA hCanonicalB
    hAIndex hBIndex hrun hsum
  have hsumBound := classicalDotNat_bound heap A B lenA lenB k stop j
    this._p sum hA hB hCanonicalA hCanonicalB hAIndex hBIndex hsum
  have hlazy := lazyAccumulation_word3_budget this._p (stop + 1 - j) 0
    hp hcount (by omega)
  have hvalue : word3Value result < this._p.toNat * limbBase ^ 2 := by
    rw [hexact]
    exact lt_of_le_of_lt hsumBound (by simpa using hlazy)
  have hhi : result.hi.toNat < this._p.toNat :=
    word3_hi_lt_of_value_lt result this._p hvalue
  rw [lll_mod_preinv_ir_correct_of_configured this result.hi result.mid
    result.lo hcfg hhi, hexact]

theorem classicalDotReduced_cast (this : DenseUPolyZp)
    (heap : RawHeap) (A B : RawPtr UInt64)
    (lenA lenB k stop j : Nat) (result : Word3) (sum : Nat)
    (hcfg : DensePreinvConfigured this)
    (hp : 1 < this._p.toNat)
    (hcount : stop + 1 - j < limbBase)
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB)
    (hCanonicalA : CanonicalU64Prefix heap A lenA this._p)
    (hCanonicalB : CanonicalU64Prefix heap B lenB this._p)
    (hAIndex : ∀ t, j ≤ t → t ≤ stop → t < lenA)
    (hBIndex : ∀ t, j ≤ t → t ≤ stop → k - t < lenB)
    (hrun : classicalDotLoop heap A B k stop j
      { lo := 0, mid := 0, hi := 0 } = .ok result)
    (hsum : classicalDotNat heap A B k stop j = .ok sum) :
    ((Generated.StrictGCD.dense_upoly_zp__lll_mod_preinv_ir
      result.hi result.mid result.lo this._p this._ninv this._norm).toNat :
        ZMod this._p.toNat) = (sum : ZMod this._p.toNat) := by
  rw [classicalDotReduced_toNat this heap A B lenA lenB k stop j result sum
    hcfg hp hcount hA hB hCanonicalA hCanonicalB hAIndex hBIndex hrun hsum,
    ZMod.natCast_mod]

theorem classical_index_bounds (lenA lenB k j : Nat)
    (hApos : 0 < lenA) (hBpos : 0 < lenB)
    (hk : k < lenA + lenB - 1)
    (hjMin : (if k ≥ lenB then k - lenB + 1 else 0) ≤ j)
    (hjMax : j ≤ (if k < lenA then k else lenA - 1)) :
    j < lenA ∧ k - j < lenB := by
  split at hjMax
  next hkA =>
    constructor
    · omega
    · split at hjMin <;> omega
  next hkNotA =>
    constructor
    · omega
    · split at hjMin <;> omega

theorem classicalReduced_source_eq_coeff (this : DenseUPolyZp)
    (heap : RawHeap) (A B : RawPtr UInt64) (lenA lenB lenC k : Nat)
    (left right : Polynomial (ZMod this._p.toNat)) (acc : Word3)
    (hcfg : DensePreinvConfigured this)
    (hp : 1 < this._p.toNat)
    (hApos : 0 < lenA) (hBpos : 0 < lenB)
    (hLenAWord : lenA < limbBase)
    (hlenC : lenC = lenA + lenB - 1) (hk : k < lenC)
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB)
    (hCanonicalA : CanonicalU64Prefix heap A lenA this._p)
    (hCanonicalB : CanonicalU64Prefix heap B lenB this._p)
    (hRepA : SlicePolyRep heap A lenA this._p.toNat left)
    (hRepB : SlicePolyRep heap B lenB this._p.toNat right)
    (hdot : classicalDotLoop heap A B k
      (if k < lenA then k else lenA - 1)
      (if k ≥ lenB then k - lenB + 1 else 0)
      { lo := 0, mid := 0, hi := 0 } = .ok acc) :
    ((Generated.StrictGCD.dense_upoly_zp__lll_mod_preinv_ir
      acc.hi acc.mid acc.lo this._p this._ninv this._norm).toNat :
        ZMod this._p.toNat) = (left * right).coeff k := by
  let jMin := if k ≥ lenB then k - lenB + 1 else 0
  let jMax := if k < lenA then k else lenA - 1
  have hkSource : k < lenA + lenB - 1 := by omega
  have hrange : jMin ≤ jMax := by
    dsimp [jMin, jMax]
    split <;> split <;> omega
  have hAIndex : ∀ t, jMin ≤ t → t ≤ jMax → t < lenA := by
    intro t hjt hts
    exact (classical_index_bounds lenA lenB k t hApos hBpos hkSource
      hjt hts).1
  have hBIndex : ∀ t, jMin ≤ t → t ≤ jMax → k - t < lenB := by
    intro t hjt hts
    exact (classical_index_bounds lenA lenB k t hApos hBpos hkSource
      hjt hts).2
  rcases classicalDotNat_ok heap A B lenA lenB k jMax jMin hA hB
    hAIndex hBIndex with ⟨sum, hsum⟩
  have hcount : jMax + 1 - jMin < limbBase := by
    have hjMaxA := hAIndex jMax hrange (Nat.le_refl _)
    omega
  calc
    ((Generated.StrictGCD.dense_upoly_zp__lll_mod_preinv_ir
      acc.hi acc.mid acc.lo this._p this._ninv this._norm).toNat :
        ZMod this._p.toNat) = (sum : ZMod this._p.toNat) :=
      classicalDotReduced_cast this heap A B lenA lenB k jMax jMin acc sum
        hcfg hp hcount hA hB hCanonicalA hCanonicalB hAIndex hBIndex
        (by simpa [jMin, jMax] using hdot) hsum
    _ = classicalDotPoly left right k jMax jMin :=
      classicalDotNat_cast_eq_poly heap A B lenA lenB this._p.toNat k
        jMax jMin sum left right hRepA hRepB hAIndex hBIndex hsum
    _ = (left * right).coeff k := by
      simpa [jMin, jMax] using classicalDotPoly_source_eq_coeff heap A B
        lenA lenB k left right hApos hRepA hRepB

/-- The generated schoolbook outer loop writes only `C[k..lenC)`.  This
memory fact keeps both source buffers and already-produced coefficients tied
to their original raw cells throughout the actual C++ loop. -/
theorem classicalOuterLoop_preserves_outside (this : DenseUPolyZp)
    (C A B guard : RawPtr UInt64) (lenA lenB lenC k readIndex : Nat)
    (heap heap' : RawHeap) (old : UInt64)
    (hC : heap.ValidU64Slice C lenC)
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB)
    (hread : heap.readU64 guard readIndex = .ok old)
    (houtside : ∀ i, k ≤ i → i < lenC →
      C.region ≠ guard.region ∨
        C.limbOffset + i ≠ guard.limbOffset + readIndex)
    (hrun : classicalOuterLoop this C A B lenA lenB lenC k heap =
      .ok heap') :
    heap'.readU64 guard readIndex = .ok old := by
  unfold classicalOuterLoop at hrun
  split at hrun
  next hk =>
    let jMin := if k ≥ lenB then k - lenB + 1 else 0
    let jMax := if k < lenA then k else lenA - 1
    cases hdot : classicalDotLoop heap A B k jMax jMin
        { lo := 0, mid := 0, hi := 0 } with
    | error fault => simp [jMin, jMax, hdot] at hrun
    | ok acc =>
      simp only [jMin, jMax, hdot] at hrun
      let value := Generated.StrictGCD.dense_upoly_zp__lll_mod_preinv_ir
        acc.hi acc.mid acc.lo this._p this._ninv this._norm
      rcases heap.writeU64_of_valid C lenC k value hC hk with ⟨heap1, hw⟩
      simp only [value, hw] at hrun
      have hread1 := RawHeap.readU64_writeU64_ne heap heap1 C guard
        k readIndex value old hw hread (houtside k (Nat.le_refl _) hk)
      have hlayout := RawHeap.writeU64_sameLayout heap heap1 C k value hw
      apply classicalOuterLoop_preserves_outside this C A B guard lenA lenB
        lenC (k + 1) readIndex heap1 heap' old
        ((hlayout C lenC).mp hC) ((hlayout A lenA).mp hA)
        ((hlayout B lenB).mp hB) hread1
      · intro i hik hil
        exact houtside i (by omega) hil
      · exact hrun
  next hnot =>
    have heq : heap' = heap := Except.ok.inj hrun.symm
    simpa [heq] using hread
termination_by lenC - k
decreasing_by omega

theorem classicalOuterLoop_same_prefix_region_ne (this : DenseUPolyZp)
    (C A B guard : RawPtr UInt64) (lenA lenB lenC k guardLen : Nat)
    (heap heap' : RawHeap)
    (hC : heap.ValidU64Slice C lenC)
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB)
    (hregions : C.region ≠ guard.region)
    (hrun : classicalOuterLoop this C A B lenA lenB lenC k heap =
      .ok heap') :
    SameU64Prefix heap heap' guard guardLen := by
  intro i old hi hread
  apply classicalOuterLoop_preserves_outside this C A B guard lenA lenB
    lenC k i heap heap' old hC hA hB hread
  · intro _ _ _
    exact Or.inl hregions
  · exact hrun

theorem canonicalU64Prefix_of_same_prefix (before after : RawHeap)
    (ptr : RawPtr UInt64) (length : Nat) (p : UInt64)
    (hvalid : before.ValidU64Slice ptr length)
    (hsame : SameU64Prefix before after ptr length)
    (hcanonical : CanonicalU64Prefix before ptr length p) :
    CanonicalU64Prefix after ptr length p := by
  intro k value hk hreadAfter
  rcases before.readU64_of_valid ptr length k
      hvalid hk with ⟨old, hreadBefore⟩
  have hpreserved := hsame k old hk hreadBefore
  have hvalue : value = old := Except.ok.inj (hreadAfter.symm.trans hpreserved)
  subst value
  exact hcanonical k old hk hreadBefore

/-- Coefficient-level output invariant for the prefix already produced by
the C++ schoolbook outer loop. -/
def ClassicalCoeffPrefix {p : Nat} (heap : RawHeap)
    (C : RawPtr UInt64) (upto : Nat) (poly : Polynomial (ZMod p)) : Prop :=
  ∀ i, i < upto → ∃ value : UInt64,
    heap.readU64 C i = .ok value ∧
      (value.toNat : ZMod p) = poly.coeff i

theorem classicalCoeffPrefix_succ_of_write {p : Nat}
    (before after : RawHeap) (C : RawPtr UInt64) (upto : Nat)
    (poly : Polynomial (ZMod p)) (value : UInt64)
    (hprefix : ClassicalCoeffPrefix before C upto poly)
    (hwrite : before.writeU64 C upto value = .ok after)
    (hvalue : (value.toNat : ZMod p) = poly.coeff upto) :
    ClassicalCoeffPrefix after C (upto + 1) poly := by
  intro i hi
  by_cases heq : i = upto
  · subst i
    exact ⟨value, RawHeap.readU64_writeU64_same before after C upto value
      hwrite, hvalue⟩
  · have hiOld : i < upto := by omega
    rcases hprefix i hiOld with ⟨old, hread, hold⟩
    refine ⟨old, ?_, hold⟩
    exact RawHeap.readU64_writeU64_ne before after C C upto i value old
      hwrite hread (Or.inr (by omega))

theorem classicalOuterLoop_preserves_coeff_prefix {p : Nat}
    (this : DenseUPolyZp) (C A B : RawPtr UInt64)
    (lenA lenB lenC k : Nat) (heap heap' : RawHeap)
    (poly : Polynomial (ZMod p))
    (hC : heap.ValidU64Slice C lenC)
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB)
    (hprefix : ClassicalCoeffPrefix heap C k poly)
    (hrun : classicalOuterLoop this C A B lenA lenB lenC k heap =
      .ok heap') :
    ClassicalCoeffPrefix heap' C k poly := by
  intro i hi
  rcases hprefix i hi with ⟨old, hread, hold⟩
  refine ⟨old, ?_, hold⟩
  apply classicalOuterLoop_preserves_outside this C A B C lenA lenB lenC
    k i heap heap' old hC hA hB hread
  · intro target hktarget _
    exact Or.inr (by omega)
  · exact hrun

theorem classicalOuterLoop_ok (this : DenseUPolyZp)
    (C A B : RawPtr UInt64) (lenA lenB lenC k : Nat) (heap : RawHeap)
    (hApos : 0 < lenA) (hBpos : 0 < lenB)
    (hlenC : lenC = lenA + lenB - 1)
    (hC : heap.ValidU64Slice C lenC)
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB) :
    ∃ heap', classicalOuterLoop this C A B lenA lenB lenC k heap =
        .ok heap' ∧ RawHeap.SameLayout heap heap' := by
  unfold classicalOuterLoop
  split
  next hk =>
    let jMin := if k ≥ lenB then k - lenB + 1 else 0
    let jMax := if k < lenA then k else lenA - 1
    have hrange : jMin ≤ jMax := by
      dsimp [jMin, jMax]
      split <;> split <;> omega
    rcases classicalDotLoop_ok heap A B lenA lenB k jMax jMin
      { lo := 0, mid := 0, hi := 0 } hA hB
      (by
        intro t hjt hts
        exact (classical_index_bounds lenA lenB k t hApos hBpos
          (by omega) hjt hts).1)
      (by
        intro t hjt hts
        exact (classical_index_bounds lenA lenB k t hApos hBpos
          (by omega) hjt hts).2) with ⟨acc, hdot⟩
    simp only [jMin, jMax, hdot]
    let value := Generated.StrictGCD.dense_upoly_zp__lll_mod_preinv_ir
      acc.hi acc.mid acc.lo this._p this._ninv this._norm
    rcases heap.writeU64_of_valid C lenC k value hC hk with ⟨heap1, hw⟩
    simp only [value, hw]
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 C k value hw
    rcases classicalOuterLoop_ok this C A B lenA lenB lenC (k + 1) heap1
      hApos hBpos hlenC ((hlayout1 C lenC).mp hC)
      ((hlayout1 A lenA).mp hA) ((hlayout1 B lenB).mp hB) with
      ⟨heap2, hrun, hlayout2⟩
    rw [hrun]
    exact ⟨heap2, rfl, fun ptr length =>
      (hlayout1 ptr length).trans (hlayout2 ptr length)⟩
  next => exact ⟨heap, rfl, fun _ _ => Iff.rfl⟩
termination_by lenC - k
decreasing_by omega

theorem classicalMul_ok (this : DenseUPolyZp)
    (C A : RawPtr UInt64) (lenA : Nat) (B : RawPtr UInt64) (lenB : Nat)
    (heap : RawHeap) (hApos : 0 < lenA) (hBpos : 0 < lenB)
    (hC : heap.ValidU64Slice C (lenA + lenB - 1))
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB) :
    ∃ heap', dense_upoly_zp__classical_mul_ir this C A lenA B lenB heap =
        .ok heap' ∧ RawHeap.SameLayout heap heap' := by
  rcases classicalOuterLoop_ok this C A B lenA lenB
    (lenA + lenB - 1) 0 heap hApos hBpos rfl hC hA hB with
    ⟨heap', hrun, hlayout⟩
  refine ⟨heap', ?_, hlayout⟩
  simp [dense_upoly_zp__classical_mul_ir, Nat.ne_of_gt hApos,
    Nat.ne_of_gt hBpos, hrun]

end CLPoly.Impl.StrictMulRefinement
