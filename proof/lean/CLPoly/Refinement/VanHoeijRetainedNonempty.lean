import CLPoly.Refinement.VanHoeijSeparationInput
import CLPoly.Refinement.VanHoeijShortVector
import CLPoly.Refinement.VanHoeijDivisorSupport
import CLPoly.Refinement.VanHoeijValidationNonunit

set_option autoImplicit false

namespace Refinement.StrictRecombine

/-- A nonzero short initial lattice vector prevents the actual suffix scan
from removing every row. -/
theorem lllReduceWithNorms_retained_nonempty_of_short
    (matrix : Generated.StrictRecombine.LLLMatrix)
    (hinput : ConcreteLLLInputValid matrix) (bound : ZZ)
    (output : Generated.StrictRecombine.LLLMatrix ×
      Generated.StrictRecombine.LLLMatrix × Array Nat × Array QQ)
    (hrun : Generated.StrictVanHoeij.lllReduceWithNorms concreteLLLExecution matrix hinput bound = .ok output)
    (coefficients : Fin matrix.size → ZZ) (hne : coefficients ≠ 0)
    (hshort : dotProduct
      (fun column => ((Matrix.vecMul coefficients
        (basisPrefixMatrix matrix matrix.size matrix.size) column : ZZ) : QQ))
      (fun column => ((Matrix.vecMul coefficients
        (basisPrefixMatrix matrix matrix.size matrix.size) column : ZZ) : QQ)) ≤ (bound : QQ)) :
    (Generated.StrictVanHoeij.shortSpanFromNorms output.2.2.2 bound).isEmpty = false := by
  obtain ⟨state, hv, _, _, hn, hrel, hu⟩ :=
    lllReduceWithNorms_invariants matrix hinput bound output hrun
  rw [← hn]
  let reduced := integerReducedCoordinates matrix state.transform coefficients
  have hne' : reduced ≠ 0 := integerReducedCoordinates_ne_zero matrix state.matrix
    state.transform hinput hrel hu coefficients hne
  have hquadratic : dotProduct (fun i => (reduced i : QQ))
      (Matrix.mulVec (gramPrefixMatrixQQ state.matrix matrix.size)
        (fun i => (reduced i : QQ))) ≤ (bound : QQ) := by
    rw [gramPrefix_quadratic_eq_vecMul_norm, hrel.1]
    have hreconstruct := integerReducedCoordinates_reconstruct matrix state.matrix
      state.transform hrel hu coefficients
    change Matrix.vecMul reduced _ = _ at hreconstruct
    rw [hreconstruct]
    exact hshort
  obtain ⟨last, hlast, hafter⟩ := exists_last_nonzero reduced hne'
  have hlastBound := gramSchmidt_norm_le_of_last_nonzero state hv matrix.size
    (by rw [hrel.1]) reduced last hlast hafter
  cases he : (Generated.StrictVanHoeij.shortSpanFromNorms state.norms bound).isEmpty with
  | false => rfl
  | true =>
    have hempty := Array.isEmpty_iff.mp he
    have hlarge := (shortSpanFromNorms_spec state.norms bound).2.2 last.val
      (by rw [hv.norms_size, hrel.1]; exact last.isLt)
      (by rw [hempty]; exact Nat.zero_le _)
    exact ((not_le_of_gt hlarge) (hlastBound.trans hquadratic)).elim

/-- A positive-degree divisor with nonvanishing leading coefficient cannot
have an empty physical modular support. -/
theorem divisor_modular_support_nonempty
    (prime modulus : Nat) [Fact (Nat.Prime prime)] (hd : prime ∣ modulus)
    (lifted : Array SparsePolyZZ) (indices : Array Nat)
    (hlegal : LegalCombination lifted.size indices.size indices)
    (factor : Polynomial Int) (hdegree : 0 < factor.natDegree)
    (hleading : (factor.leadingCoeff : ZMod prime) ≠ 0)
    (hsupport : factor.map (Int.castRingHom (ZMod modulus)) =
      Polynomial.C (factor.leadingCoeff : ZMod modulus) *
        (indices.toList.map fun i => StrictHensel.toPolyMod modulus lifted[i]!).prod) :
    indices.toList ≠ [] := by
  have ha := divisor_support_project_associated prime modulus hd lifted indices
    hlegal factor hleading hsupport
  rw [legalPhysicalSupport_product] at ha
  intro he
  simp only [he, List.map_nil, List.prod_nil] at ha
  have hu := ha.isUnit_iff.mpr isUnit_one
  have hz := Polynomial.natDegree_eq_zero_of_isUnit hu
  rw [Polynomial.natDegree_map_of_leadingCoeff_ne_zero _ hleading] at hz
  omega

/-- The whole source polynomial supplies a nonzero preserved vector. Thus
actual Hensel/CLD/LLL execution derives nonemptiness without an irreducible
output assumption or a caller-provided retained-space certificate. -/
theorem hensel_lll_retained_nonempty
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {aTarget : Int32} {lifted : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hfactors : StrictSelectPrime.SelectionPhysical selection)
    (hleading : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) = (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f selection.factors selection.prime aTarget lifted)
    (modulus : Nat) (hmodulus : lifted.2 = (modulus : Int))
    (original truncated : Array SparsePolyZZ)
    (data : CldSeparationData f lifted.1 original modulus)
    (scale divisor height bound : ZZ) (width columns : Nat)
    (matrix : Generated.StrictRecombine.LLLMatrix)
    (htruncate : Generated.StrictVanHoeij.truncateCld original divisor = .ok truncated)
    (hcanonical : CanonicalCldBasis truncated width scale ((modulus : Int) / divisor) columns matrix)
    (hheight : height = ((SparsePolyZZ.toPoly f).natDegree *
      2 ^ ((SparsePolyZZ.toPoly f).natDegree - 1) *
      (Finset.range ((SparsePolyZZ.toPoly f).natDegree + 1)).sum
        (fun k => ((SparsePolyZZ.toPoly f).coeff k).natAbs) : Nat))
    (hh : height < (modulus : Int))
    (hcutoff : (lifted.1.size : ZZ) * scale * scale + (columns : ZZ) *
      ((height + divisor - 1) / divisor + 2 * (lifted.1.size : Int)) *
      ((height + divisor - 1) / divisor + 2 * (lifted.1.size : Int)) ≤ bound)
    (hinput : ConcreteLLLInputValid matrix)
    (output : Generated.StrictRecombine.LLLMatrix ×
      Generated.StrictRecombine.LLLMatrix × Array Nat × Array QQ)
    (hrun : Generated.StrictVanHoeij.lllReduceWithNorms concreteLLLExecution matrix hinput bound = .ok output) :
    (Generated.StrictVanHoeij.shortSpanFromNorms output.2.2.2 bound).isEmpty = false := by
  obtain ⟨indices, coefficients, hlegal, hsupport, hcoordinates, hshort⟩ :=
    hensel_factor_truncated_short hcount hfactors hleading hselection hentry original truncated
      modulus hmodulus divisor height scale width columns matrix data.size_eq
      (fun i hi => data.values i (by rw [data.size_eq]; exact hi)) htruncate hcanonical
      (SparsePolyZZ.toPoly f) 1 (by simp) (by have h := data.modulus_gt_one; omega)
      (truncateCld_spec original truncated divisor htruncate).1 hh data.source_degree_pos hheight
  obtain ⟨exponent, he, hp⟩ := hentry.outputModulus_eq_prime_pow
  have hp' : modulus = selection.prime.toNat ^ exponent := by
    exact_mod_cast hmodulus.symm.trans hp
  have hd : selection.prime.toNat ∣ modulus := by
    rw [hp']; exact dvd_pow_self _ (Nat.ne_of_gt he)
  have hnonempty := divisor_modular_support_nonempty selection.prime.toNat modulus hd
    lifted.1 indices hlegal (SparsePolyZZ.toPoly f) data.source_degree_pos
    hselection.goodPrime.lc_nonzero hsupport
  obtain ⟨first, hfirst⟩ := List.exists_mem_of_ne_nil indices.toList hnonempty
  obtain ⟨position, hposition, hvalue⟩ := List.mem_iff_getElem.mp hfirst
  have hpSize : position < indices.size := by simpa using hposition
  have hfirstBound : first < lifted.1.size := by
    have hb := hlegal.2.2 position hpSize
    simpa only [getElem!_pos indices position hpSize, ← Array.getElem_toList hpSize, hvalue] using hb
  have htruncSize := (truncateCld_spec original truncated divisor htruncate).2.1
  have hfactorCount : lifted.1.size ≤ matrix.size := by
    rw [hcanonical.size_eq, htruncSize, data.size_eq]
    omega
  have hne : coefficients ≠ 0 := by
    intro hz
    have hc := hcoordinates ⟨first, lt_of_lt_of_le hfirstBound hfactorCount⟩ hfirstBound
    simp [hz, hfirst] at hc
  apply lllReduceWithNorms_retained_nonempty_of_short matrix hinput bound output hrun coefficients hne
  exact hshort.trans (by exact_mod_cast hcutoff)

end Refinement.StrictRecombine
