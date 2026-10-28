import CLPoly.Refinement.VanHoeijHeuristicLoop
import CLPoly.Refinement.Recombine

set_option autoImplicit false

namespace Refinement.StrictRecombine
open _root_.Generated.StrictVanHoeijHeuristic

/-- The recovery target computed from the input, before any factor is returned. -/
noncomputable def heuristicRecoveryTarget (f : SparsePolyZZ) : Int :=
  2 * |(SparsePolyZZ.toPoly f).leadingCoeff| *
    (((SparsePolyZZ.toPoly f).natDegree.choose
      ((SparsePolyZZ.toPoly f).natDegree / 2) : Int) *
      Generated.StrictHensel.__isqrt_ceil_raw_ir
        (Generated.StrictHensel.__upoly_norm_l2_sq_upoly_raw_ir f))

/-- Construct the actual prefix run and its positive start. The only precision
premise is representability of the input recovery target; there is no assumed
heuristic execution, floating-point estimate accuracy, or factorization oracle. -/
theorem startingPrecision_succeeds_positive (f : SparsePolyZZ) (r : Int32)
    (p : UInt64) (hp : 2 ≤ p.toNat)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hnonempty : 0 < f.size) (hdegree : (SparsePolyZZ.toPoly f).natDegree < 2 ^ 31)
    (hpositive : (0 : Int) < (SparsePolyZZ.toPoly f).leadingCoeff)
    (hfeasible : heuristicRecoveryTarget f < (p.toNat : Int) ^ 2147483647) :
    ∃ precision full, startingPrecision f r p = .ok (precision, full) ∧
      0 < precision := by
  have hhead : f[0]? = some f[0] := Array.getElem?_eq_getElem hnonempty
  have hd := sparsePolyZZ_natDegree_eq_head f hcanonical hnonempty
  have hl := sparsePolyZZ_leadingCoeff_eq_head f hcanonical hnonempty
  have hheadDegree : f[0].1.deg < 2 ^ 63 := by omega
  have hheadPositive : (0 : Int) < f[0].2 := by rw [← hl]; exact hpositive
  have hf : SparsePolyZZ.toPoly f ≠ 0 := by
    intro he
    simp only [he, Polynomial.leadingCoeff_zero] at hpositive
    omega
  obtain ⟨bound, hb, hbNonnegative, hbCoeff⟩ :=
    mignotteBoundRaw_bounds_divisor f hcanonical hnonempty hheadDegree 1 hf (one_dvd _)
  have hbOne : (1 : Int) ≤ bound := by simpa using hbCoeff 0
  have hbFormula := mignotteBoundRaw_eq_choose_isqrt f f[0] hhead hheadDegree
  have heBound := Except.ok.inj (hb.symm.trans hbFormula)
  have ht : heuristicRecoveryTarget f = (2 : Int) * f[0].2 * bound := by
    unfold heuristicRecoveryTarget
    rw [hd, hl, abs_of_pos hheadPositive, heBound]
  have htPositive : (1 : Int) ≤ (2 : Int) * f[0].2 * bound := by
    simp only [ZZ] at hbOne
    nlinarith
  obtain ⟨full, hfull⟩ := fullPrecision_succeeds p hp ((2 : Int) * f[0].2 * bound)
    (by rwa [← ht])
  have hfullPositive := fullPrecision_positive p hp ((2 : Int) * f[0].2 * bound)
    htPositive full hfull
  unfold startingPrecision
  rw [dif_pos hp, hb, hhead]
  simp only [not_lt_of_ge (le_of_lt hheadPositive), ↓reduceIte, hfull]
  refine ⟨_, full, rfl, ?_⟩
  apply Int32.lt_iff_toInt_lt.mpr
  rw [Int32.toInt_zero, int32_min_toInt]
  have hfpos := Int32.lt_iff_toInt_lt.mp hfullPositive
  change (0 : Int) < full.toInt at hfpos
  apply lt_min hfpos
  exact Int32.lt_iff_toInt_lt.mp (clampPrecision_positive _)

end Refinement.StrictRecombine
