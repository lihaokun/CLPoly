/-
Execution facts for the generated C++ `__upoly_random` loop.
-/
import CLPoly.Generated.StrictEDF
import CLPoly.Math.Univariate

set_option autoImplicit false

namespace Refinement.StrictEDF

open CLPoly.Math

private theorem zpOfRandomReduced (p coefficient : UInt64)
    (hcoefficient : coefficient < p) :
    Zp.Reduced p.toNat (Zp.ofInt coefficient.toInt p) := by
  constructor
  · simp [Zp.ofInt]
  · simp only [Zp.ofInt, UInt64.toInt]
    rw [Int.emod_eq_of_lt]
    · change coefficient.toNat.toUInt64.toNat < p.toNat
      simpa [Nat.mod_eq_of_lt coefficient.toNat_lt_size] using hcoefficient
    · omega
    · exact_mod_cast hcoefficient

private theorem rngNextLt (p rng : UInt64) (hp : 0 < p) :
    Rng.next rng p < p := by
  unfold Rng.next
  have hp0 : p ≠ 0 := by
    intro h
    subst p
    simp at hp
  have hb : (p == 0) = false := beq_eq_false_iff_ne.mpr hp0
  rw [hb]
  simp only [Bool.false_eq_true, ↓reduceIte]
  exact UInt64.mod_lt _ hp

private theorem canonicalPushDescending (p degree : Nat)
    (result : SparsePolyZp) (coefficient : Zp)
    (hresult : SparsePolyZp.Canonical p result)
    (hdegrees : ∀ term ∈ result.toList, degree < term.1.deg)
    (hreduced : Zp.Reduced p coefficient)
    (hnonzero : coefficient.val ≠ 0) :
    SparsePolyZp.Canonical p
      (result.push (UMonomial.mk degree, coefficient)) := by
  rcases hresult with ⟨hwellFormed, hchain, hnonzeroResult⟩
  refine ⟨?_, ?_, ?_⟩
  · rw [SparsePolyZp.WellFormed_arr, Array.toList_push]
    intro term hterm
    rcases List.mem_append.mp hterm with hterm | hterm
    · exact hwellFormed term hterm
    · simp only [List.mem_singleton] at hterm
      subst term
      exact hreduced
  · rw [Array.toList_push]
    apply List.isChain_iff_pairwise.mpr
    rw [List.pairwise_append]
    refine ⟨List.isChain_iff_pairwise.mp hchain, by simp, ?_⟩
    intro earlier hearlier final hfinal
    simp only [List.mem_singleton] at hfinal
    subst final
    exact hdegrees earlier hearlier
  · rw [Array.toList_push]
    intro term hterm
    rcases List.mem_append.mp hterm with hterm | hterm
    · exact hnonzeroResult term hterm
    · simp only [List.mem_singleton] at hterm
      subst term
      exact hnonzero

private theorem randomLoopCanonical (remaining degree : Nat)
    (result : SparsePolyZp) (p : UInt64) (rng : Rng)
    (hp : 0 < p)
    (hposition : remaining = 0 ∨ degree + 1 = remaining)
    (hcanonical : SparsePolyZp.Canonical p.toNat result)
    (hdegrees : ∀ term ∈ result.toList, degree < term.1.deg) :
    SparsePolyZp.Canonical p.toNat
      (Generated.StrictEDF._loop___upoly_random_0_raw_ir
        remaining degree result p rng).1 := by
  induction remaining generalizing degree result rng with
  | zero =>
      simpa [Generated.StrictEDF._loop___upoly_random_0_raw_ir] using hcanonical
  | succ remaining ih =>
      have hdegree : degree = remaining := by omega
      subst degree
      rw [Generated.StrictEDF._loop___upoly_random_0_raw_ir]
      simp only [Nat.succ_ne_zero, ↓reduceDIte]
      by_cases hzero : Rng.next rng p = 0
      · have hb : (Rng.next rng p != 0) = false :=
          bne_eq_false_iff_eq.mpr hzero
        simp only [Rng.next_advance, hb, Bool.false_eq_true, ↓reduceIte,
          Nat.succ_sub_one]
        apply ih (degree := remaining - 1) (result := result)
          (rng := Rng.step rng)
        · omega
        · exact hcanonical
        · intro term hterm
          exact lt_of_le_of_lt (Nat.sub_le _ _) (hdegrees term hterm)
      · let coefficient := Zp.ofInt (Rng.next rng p).toInt p
        have hcoefficientReduced : Zp.Reduced p.toNat coefficient :=
          zpOfRandomReduced p (Rng.next rng p) (rngNextLt p rng hp)
        have hcoefficientNonzero : coefficient.val ≠ 0 := by
          intro h
          have hvalue : coefficient.val.toNat = (Rng.next rng p).toNat := by
            unfold coefficient
            simp only [Zp.ofInt, UInt64.toInt]
            rw [Int.emod_eq_of_lt]
            · change (Rng.next rng p).toNat.toUInt64.toNat =
                (Rng.next rng p).toNat
              simp
            · omega
            · exact_mod_cast rngNextLt p rng hp
          have : (Rng.next rng p).toNat = 0 := by
            rw [← hvalue, h]
            rfl
          exact hzero (UInt64.toNat_inj.mp this)
        have hnextCanonical : SparsePolyZp.Canonical p.toNat
            (result.push (UMonomial.mk remaining, coefficient)) :=
          canonicalPushDescending p.toNat remaining result coefficient
            hcanonical hdegrees hcoefficientReduced hcoefficientNonzero
        have hb : (Rng.next rng p != 0) = true :=
          Bool.eq_true_of_not_eq_false (fun hb =>
            hzero (bne_eq_false_iff_eq.mp hb))
        simp only [Rng.next_advance, hb, ↓reduceIte, Nat.succ_sub_one]
        change SparsePolyZp.Canonical p.toNat
          (Generated.StrictEDF._loop___upoly_random_0_raw_ir remaining
            (remaining - 1)
            (result.push (UMonomial.mk remaining, coefficient)) p
            (Rng.step rng)).1
        by_cases hremaining : remaining = 0
        · subst remaining
          simpa [Generated.StrictEDF._loop___upoly_random_0_raw_ir] using
            hnextCanonical
        apply ih (degree := remaining - 1)
          (result := result.push (UMonomial.mk remaining, coefficient))
          (rng := Rng.step rng)
        · omega
        · exact hnextCanonical
        · intro term hterm
          rw [Array.toList_push] at hterm
          rcases List.mem_append.mp hterm with hterm | hterm
          · exact lt_of_le_of_lt (Nat.sub_le _ _) (hdegrees term hterm)
          · simp only [List.mem_singleton] at hterm
            subst term
            simp
            exact Nat.pos_of_ne_zero hremaining

private theorem randomLoopDegreeBound (upper remaining degree : Nat)
    (result : SparsePolyZp) (p : UInt64) (rng : Rng)
    (hposition : remaining = 0 ∨ degree + 1 = remaining)
    (hcursor : remaining = 0 ∨ degree < upper)
    (hresult : ∀ term ∈ result.toList, term.1.deg < upper) :
    ∀ term ∈
      (Generated.StrictEDF._loop___upoly_random_0_raw_ir
        remaining degree result p rng).1.toList,
      term.1.deg < upper := by
  induction remaining generalizing degree result rng with
  | zero =>
      simpa [Generated.StrictEDF._loop___upoly_random_0_raw_ir] using hresult
  | succ remaining ih =>
      have hdegree : degree = remaining := by omega
      have hdegreeUpper : degree < upper := by omega
      subst degree
      rw [Generated.StrictEDF._loop___upoly_random_0_raw_ir]
      simp only [Nat.succ_ne_zero, ↓reduceDIte, Rng.next_advance,
        Nat.succ_sub_one]
      by_cases hnonzero : (Rng.next rng p != 0) = true
      · simp only [hnonzero, ↓reduceIte]
        apply ih (degree := remaining - 1)
          (result := result.push (UMonomial.mk remaining,
            Zp.ofInt (Rng.next rng p).toInt p))
          (rng := Rng.step rng)
        · omega
        · omega
        · intro term hterm
          rw [Array.toList_push] at hterm
          rcases List.mem_append.mp hterm with hterm | hterm
          · exact hresult term hterm
          · simp only [List.mem_singleton] at hterm
            subst term
            exact hdegreeUpper
      · have hzero : (Rng.next rng p != 0) = false :=
          Bool.eq_false_of_not_eq_true hnonzero
        simp only [hzero, Bool.false_eq_true, ↓reduceIte]
        apply ih (degree := remaining - 1) (result := result)
          (rng := Rng.step rng)
        · omega
        · omega
        · exact hresult

/-- The exact generated random-polynomial entry is total.  Its output and
advanced RNG state are the values computed by the well-founded source loop. -/
theorem __upoly_random_raw_ir_terminates
    (maxDegree : Int64) (p : UInt64) (rng : Rng) :
    ∃ output rngNext,
      Generated.StrictEDF.__upoly_random_raw_ir maxDegree p rng =
        .ok (output, rngNext) := by
  unfold Generated.StrictEDF.__upoly_random_raw_ir
  split
  · exact ⟨_, _, rfl⟩
  · exact ⟨#[], rng, rfl⟩

/-- The generated random entry returns a canonical sparse polynomial.  This
is proved from its actual coefficient draws and descending push order; no
normalization or L2 result substitution is inserted after execution. -/
theorem __upoly_random_raw_ir_canonical
    (maxDegree : Int64) (p : UInt64) (rng : Rng) (hp : 0 < p) :
    ∃ output rngNext,
      Generated.StrictEDF.__upoly_random_raw_ir maxDegree p rng =
        .ok (output, rngNext) ∧
      SparsePolyZp.Canonical p.toNat output ∧
      (∀ term ∈ output.toList,
        term.1.deg < maxDegree.toNatClampNeg) := by
  unfold Generated.StrictEDF.__upoly_random_raw_ir
  split
  · let count := maxDegree.toNatClampNeg
    let run := Generated.StrictEDF._loop___upoly_random_0_raw_ir
      count (count - 1) #[] p rng
    refine ⟨run.1, run.2, rfl, ?_, ?_⟩
    · apply randomLoopCanonical count (count - 1) #[] p rng hp
      · omega
      · simp [SparsePolyZp.Canonical, SparsePolyZp.WellFormed_arr,
          SparsePolyZp.AllReduced]
      · simp
    · apply randomLoopDegreeBound count count (count - 1) #[] p rng
      · omega
      · omega
      · simp
  · refine ⟨#[], rng, rfl, ?_, ?_⟩
    · simp [SparsePolyZp.Canonical, SparsePolyZp.WellFormed_arr,
        SparsePolyZp.AllReduced]
    · simp

end Refinement.StrictEDF
