import CLPoly.Refinement.VanHoeijResidualPolynomial

set_option autoImplicit false

namespace Refinement.StrictRecombine

/-- Monic logarithmic derivatives have degree below the source bound,
also over rings with zero divisors. -/
theorem monic_cld_natDegree_lt {R : Type*} [CommRing R] [Nontrivial R]
    (f q : Polynomial R) (hq : q.Monic) (limit : Nat)
    (hlimit : 0 < limit) (hf : f.natDegree ≤ limit) :
    ((f /ₘ q) * Polynomial.derivative q).natDegree < limit := by
  by_cases hq0 : q.natDegree = 0
  · rw [Polynomial.eq_C_of_natDegree_eq_zero hq0, Polynomial.derivative_C, mul_zero,
      Polynomial.natDegree_zero]
    exact hlimit
  by_cases hqf : q.natDegree ≤ f.natDegree
  · have hm : ((f /ₘ q) * Polynomial.derivative q).natDegree ≤
        (f /ₘ q).natDegree + (Polynomial.derivative q).natDegree := Polynomial.natDegree_mul_le
    rw [Polynomial.natDegree_divByMonic f hq] at hm
    have hd := Polynomial.natDegree_derivative_le q
    omega
  · have hz : f /ₘ q = 0 := (Polynomial.divByMonic_eq_zero_iff hq).mpr
      (Polynomial.degree_lt_degree (by omega))
    simpa only [hz, zero_mul, Polynomial.natDegree_zero] using hlimit

/-- The actual symmetric coefficient bound rules out hidden high-degree
integer coefficients whose residues vanish. -/
theorem ModularCldValue.natDegree_lt
    {f divisor output : SparsePolyZZ} {modulus : Nat}
    (hvalid : ModularCldValue f divisor output modulus) (hm : 1 < modulus)
    (hmonic : (StrictHensel.toPolyMod modulus divisor).Monic)
    (hdegree : 0 < (SparsePolyZZ.toPoly f).natDegree) :
    (SparsePolyZZ.toPoly output).natDegree < (SparsePolyZZ.toPoly f).natDegree := by
  haveI : Fact (1 < modulus) := ⟨hm⟩
  have hcld := monic_cld_natDegree_lt (StrictHensel.toPolyMod modulus f)
    (StrictHensel.toPolyMod modulus divisor) hmonic (SparsePolyZZ.toPoly f).natDegree
    hdegree Polynomial.natDegree_map_le
  have hb : (SparsePolyZZ.toPoly output).natDegree ≤ (SparsePolyZZ.toPoly f).natDegree - 1 := by
    apply Polynomial.natDegree_le_iff_coeff_eq_zero.mpr
    intro d hd
    have hres := congrArg (fun p : Polynomial (ZMod modulus) => p.coeff d) hvalid.semantic
    dsimp only at hres
    have hz := Polynomial.coeff_eq_zero_of_natDegree_lt
      (lt_of_lt_of_le hcld (by omega : (SparsePolyZZ.toPoly f).natDegree ≤ d))
    rw [hz] at hres
    rw [StrictHensel.toPolyMod, Polynomial.coeff_map] at hres
    change (((SparsePolyZZ.toPoly output).coeff d : Int) : ZMod modulus) = 0 at hres
    have hdiv := (ZMod.intCast_zmod_eq_zero_iff_dvd _ modulus).mp hres
    have hbound := hvalid.bounded d
    by_contra hne
    have hle := Int.natAbs_le_of_dvd_ne_zero hdiv hne
    have hle' : modulus ≤ ((SparsePolyZZ.toPoly output).coeff d).natAbs := by
      simpa using hle
    omega
  omega

/-- Actual CLD certificates and actual truncation now suffice to bound
the actual scanned width by the input degree. -/
theorem modularCld_truncate_width_le
    (f : SparsePolyZZ) (active original truncated : Array SparsePolyZZ)
    (modulus : Nat) (divisor : ZZ) (hm : 1 < modulus)
    (hdegree : 0 < (SparsePolyZZ.toPoly f).natDegree)
    (hcld : ∀ i, i < original.size → ModularCldValue f active[i]! original[i]! modulus)
    (hmonic : ∀ i, i < original.size → (StrictHensel.toPolyMod modulus active[i]!).Monic)
    (htruncate : Generated.StrictVanHoeij.truncateCld original divisor = .ok truncated) :
    Generated.StrictRecombine.cldSpiralWidth truncated ≤ (SparsePolyZZ.toPoly f).natDegree := by
  exact truncateCld_width_le original truncated divisor htruncate
    (fun i hi => (hcld i hi).canonical) _
    (fun i hi => (hcld i hi).natDegree_lt hm (hmonic i hi) hdegree)

theorem modularCld_residualPolynomial_natDegree_lt
    (f : SparsePolyZZ) (active original truncated : Array SparsePolyZZ)
    (modulus size : Nat) (divisor : ZZ) (coefficients : Nat → Int) (hm : 1 < modulus)
    (hdegree : 0 < (SparsePolyZZ.toPoly f).natDegree)
    (hcld : ∀ i, i < original.size → ModularCldValue f active[i]! original[i]! modulus)
    (hmonic : ∀ i, i < original.size → (StrictHensel.toPolyMod modulus active[i]!).Monic)
    (htruncate : Generated.StrictVanHoeij.truncateCld original divisor = .ok truncated) :
    (rowResidualPolynomial original size (Generated.StrictRecombine.cldSpiralWidth truncated)
      coefficients (modulus : Int)).natDegree < (SparsePolyZZ.toPoly f).natDegree := by
  exact rowResidualPolynomial_natDegree_lt original size _ _ coefficients _ hdegree
    (modularCld_truncate_width_le f active original truncated modulus divisor hm hdegree hcld hmonic htruncate)
    (fun i hi => (hcld i hi).canonical)
    (fun i hi => (hcld i hi).natDegree_lt hm (hmonic i hi) hdegree)

end Refinement.StrictRecombine
