import CLPoly.Refinement.VanHoeijShortVector

set_option autoImplicit false

namespace Refinement.StrictRecombine

theorem normL1Loop_sum (f : SparsePolyZZ) (index : Nat) (sum : ZZ) :
    Generated.StrictVanHoeij.normL1Loop f index sum =
      sum + ∑ i ∈ Finset.Ico index f.size, |f[i]!.2| := by
  induction hmeasure : f.size - index using Nat.strong_induction_on
      generalizing index sum with
  | h measure ih =>
    rw [Generated.StrictVanHoeij.normL1Loop]
    split
    next hi =>
      rw [ih (f.size - (index + 1)) (by omega) (index + 1) _ rfl]
      rw [Finset.sum_eq_sum_Ico_succ_bot hi]
      have habs : (if f[index].2 < 0 then -f[index].2 else f[index].2) =
          |f[index].2| := by split <;> simp_all [abs_of_neg, abs_of_nonneg]
      rw [habs]
      simp only [getElem!_pos f index hi]
      ring
    next hi =>
      simp [Finset.Ico_eq_empty_of_le (by omega : f.size ≤ index)]

theorem normL1_sum (f : SparsePolyZZ) :
    Generated.StrictVanHoeij.normL1 f = ∑ i ∈ Finset.range f.size, |f[i]!.2| := by
  simpa [Generated.StrictVanHoeij.normL1] using normL1Loop_sum f 0 0

/-- Literal successful height computation, including the narrowed machine
degree. Relating this degree to natDegree requires the caller's range bound. -/
theorem cldBound_output (f : SparsePolyZZ) (height : ZZ)
    (hrun : Generated.StrictVanHoeij.cldBound f = .ok height) :
    0 < (get_deg f).toInt32 ∧
      height = (get_deg f).toInt32.toInt *
        ((1 : ZZ) <<< ((get_deg f).toInt32 - 1).toInt.toNat) *
        ∑ i ∈ Finset.range f.size, |f[i]!.2| := by
  unfold Generated.StrictVanHoeij.cldBound at hrun
  dsimp only at hrun
  split at hrun
  next hn =>
    refine ⟨hn, ?_⟩
    have he := Except.ok.inj hrun
    rw [← he, normL1_sum]
  next hn => contradiction

theorem normL1_polynomial (f : SparsePolyZZ)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f) :
    Generated.StrictVanHoeij.normL1 f =
      ∑ degree ∈ Finset.range ((SparsePolyZZ.toPoly f).natDegree + 1),
        |(SparsePolyZZ.toPoly f).coeff degree| := by
  let p : Polynomial Int :=
    (f.toList.map fun term => Polynomial.monomial term.1.deg |term.2|).sum
  have hcoeff : ∀ degree, p.coeff degree = |(SparsePolyZZ.toPoly f).coeff degree| := by
    intro degree
    exact descending_terms_map_coeff f.toList hcanonical.1 (fun c => |c|)
      (by simp) degree
  have hdegree : p.natDegree ≤ (SparsePolyZZ.toPoly f).natDegree := by
    apply Polynomial.natDegree_le_iff_coeff_eq_zero.mpr
    intro degree hd
    rw [hcoeff, Polynomial.coeff_eq_zero_of_natDegree_lt hd, abs_zero]
  have hevalTerms : ∀ terms : List (UMonomial × Int),
      ((terms.map fun term => Polynomial.monomial term.1.deg |term.2|).sum).eval 1 =
        (terms.map fun term => |term.2|).sum := by
    intro terms
    induction terms with
    | nil => simp
    | cons head tail ih => simp [Polynomial.eval_add, Polynomial.eval_monomial, ih]
  have hlist : f.toList = List.ofFn (fun i : Fin f.size => f[i.val]) := by
    apply List.ext_getElem <;> simp
  have heval : p.eval 1 = Generated.StrictVanHoeij.normL1 f := by
    rw [normL1_sum]
    dsimp only [p]
    rw [hevalTerms, hlist, List.map_ofFn, List.sum_ofFn]
    rw [← Fin.sum_univ_eq_sum_range]
    apply Finset.sum_congr rfl
    intro i hi
    rw [getElem!_pos f i.val i.isLt]
    rfl
  rw [← heval, Polynomial.eval_eq_sum_range' (by omega :
    p.natDegree < (SparsePolyZZ.toPoly f).natDegree + 1)]
  simp only [hcoeff, one_pow, mul_one]

theorem cldBound_polynomial (f : SparsePolyZZ) (height : ZZ)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hnonempty : 0 < f.size)
    (hfits : (SparsePolyZZ.toPoly f).natDegree < 2 ^ 31)
    (hrun : Generated.StrictVanHoeij.cldBound f = .ok height) :
    height = ((SparsePolyZZ.toPoly f).natDegree *
      2 ^ ((SparsePolyZZ.toPoly f).natDegree - 1) *
      (Finset.range ((SparsePolyZZ.toPoly f).natDegree + 1)).sum
        (fun degree => ((SparsePolyZZ.toPoly f).coeff degree).natAbs) : Nat) := by
  let n := (SparsePolyZZ.toPoly f).natDegree
  have hcast : (get_deg f).toInt32 = Int32.ofNat n := by
    have hempty : f.isEmpty = false := by simp [Array.isEmpty, Nat.ne_of_gt hnonempty]
    unfold get_deg
    rw [hempty]
    simp only [Bool.false_eq_true, ↓reduceIte]
    change (Int64.ofNat f[0]!.1.deg).toInt32 = Int32.ofNat n
    rw [Int64.toInt32_ofNat', getElem!_pos f 0 hnonempty]
    rw [← sparsePolyZZ_natDegree_eq_head f hcanonical hnonempty]
  obtain ⟨hpositive, houtput⟩ := cldBound_output f height hrun
  rw [hcast] at hpositive houtput
  have hnpos : 0 < n := by
    rw [Int32.lt_iff_toInt_lt, Int32.toInt_ofNat_of_lt hfits] at hpositive
    simpa using hpositive
  have hsub : (Int32.ofNat n - 1).toInt.toNat = n - 1 := by
    change (Int32.ofNat n - Int32.ofNat 1).toInt.toNat = n - 1
    rw [← Int32.ofNat_sub n 1 (by omega)]
    rw [Int32.toInt_ofNat_of_lt (by omega : n - 1 < 2 ^ 31)]
    exact Int.toNat_natCast _
  rw [Int32.toInt_ofNat_of_lt hfits, hsub, Int.shiftLeft_eq, one_mul] at houtput
  have hl1 := normL1_polynomial f hcanonical
  rw [normL1_sum] at hl1
  rw [hl1] at houtput
  rw [houtput]
  simp only [n, Nat.cast_mul, Nat.cast_pow, Nat.cast_ofNat, Nat.cast_sum, Int.natCast_natAbs]

end Refinement.StrictRecombine
