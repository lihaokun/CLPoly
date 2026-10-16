import CLPoly.Refinement.VanHoeijHensel

set_option autoImplicit false

namespace Refinement.StrictRecombine

/-- Preserve physical occurrences when embedding a bounded, duplicate-free
index list into the finite row type of the lattice. -/
theorem physical_support_sum (indices : List Nat) (size : Nat)
    (hbound : ∀ index ∈ indices, index < size) (hnodup : indices.Nodup)
    (value : Nat → Int) :
    (∑ row ∈ (Finset.univ.filter fun row : Fin size => row.val ∈ indices), value row.val) =
      (indices.map value).sum := by
  classical
  calc
    _ = ∑ index ∈ indices.toFinset, value index := by
      apply Finset.sum_bij (fun row _ => row.val)
      · intro row hr
        simpa using (Finset.mem_filter.mp hr).2
      · intro left hl right hr heq
        exact Fin.ext heq
      · intro index hi
        have hm : index ∈ indices := by simpa using hi
        exact ⟨⟨index, hbound index hm⟩, by simp [hm], rfl⟩
      · intro row hr
        rfl
    _ = _ := List.sum_toFinset _ hnodup

theorem legal_support_lattice_data
    (indices : Array Nat) (upper size : Nat)
    (hlegal : LegalCombination upper indices.size indices) (hsize : upper ≤ size)
    (original : Array SparsePolyZZ) (actual : Polynomial Int) (modulus : Int)
    (hcongr : ∀ degree, ∃ carry : Int,
      (indices.toList.map fun row =>
        Generated.StrictRecombine.sparseCoeff original[row]! degree).sum =
          actual.coeff degree + carry * modulus) :
    let support := Finset.univ.filter fun row : Fin size => row.val ∈ indices.toList
    support.card ≤ upper ∧
      (∀ row ∈ support, row.val < upper) ∧
      ∀ degree, ∃ carry : Int,
        (∑ row ∈ support, Generated.StrictRecombine.sparseCoeff original[row.val]! degree) =
          actual.coeff degree + carry * modulus := by
  classical
  dsimp only
  have hbound : ∀ index ∈ indices.toList, index < upper := by
    intro index hi
    obtain ⟨position, hp, heq⟩ := List.mem_iff_getElem.mp hi
    have hp' : position < indices.size := by simpa using hp
    have hb := hlegal.2.2 position hp'
    simpa only [getElem!_pos indices position hp', ← Array.getElem_toList hp', heq] using hb
  have hboundSize : ∀ index ∈ indices.toList, index < size :=
    fun index hi => lt_of_lt_of_le (hbound index hi) hsize
  have hn := StrictFactorZZ.legalCombination_indices_nodup hlegal
  have hc := physical_support_sum indices.toList size hboundSize hn (fun _ => (1 : Int))
  have hcard : (Finset.univ.filter fun row : Fin size => row.val ∈ indices.toList).card =
      indices.size := by
    simpa using hc
  refine ⟨by rw [hcard]; exact hlegal.count_le, ?_, ?_⟩
  · intro row hr
    exact hbound row.val (Finset.mem_filter.mp hr).2
  · intro degree
    rw [physical_support_sum indices.toList size hboundSize hn
      (fun row => Generated.StrictRecombine.sparseCoeff original[row]! degree)]
    exact hcongr degree

/-- Real CLD outputs and the executed truncation loop supply the raw bounds
and floor correspondence needed by the lattice theorem. -/
theorem physical_cld_truncated_short
    (f : SparsePolyZZ) (active original truncated : Array SparsePolyZZ)
    (modulus : Nat) (divisor height scale : ZZ) (width columns : Nat)
    (matrix : Generated.StrictRecombine.LLLMatrix)
    (hsize : original.size = active.size)
    (hcld : ∀ index, index < active.size →
      ModularCldValue f active[index]! original[index]! modulus)
    (htruncate : Generated.StrictVanHoeij.truncateCld original divisor = .ok truncated)
    (hbasis : CanonicalCldBasis truncated width scale
      ((modulus : Int) / divisor) columns matrix)
    (support : Finset (Fin matrix.size))
    (hcard : support.card ≤ active.size)
    (hsupport : ∀ index ∈ support, index.val < active.size)
    (actual : Polynomial Int)
    (hcongr : ∀ degree, ∃ carry : Int,
      (∑ row ∈ support, Generated.StrictRecombine.sparseCoeff original[row.val]! degree) =
        actual.coeff degree + carry * (modulus : Int))
    (hm : 0 < modulus) (hd : 0 < divisor) (hh : height < (modulus : Int))
    (hactual : ∀ degree, |actual.coeff degree| ≤ height) :
    ∃ coefficients : Fin matrix.size → ZZ,
      (∀ index, index.val < active.size →
        coefficients index = if index ∈ support then 1 else 0) ∧
      dotProduct
        (fun index => ((Matrix.vecMul coefficients
          (basisPrefixMatrix matrix matrix.size matrix.size) index : ZZ) : QQ))
        (fun index => ((Matrix.vecMul coefficients
          (basisPrefixMatrix matrix matrix.size matrix.size) index : ZZ) : QQ)) ≤
        ((active.size : ZZ) * scale * scale + (columns : ZZ) *
          ((height + divisor - 1) / divisor + 2 * (active.size : Int)) *
          ((height + divisor - 1) / divisor + 2 * (active.size : Int)) : ZZ) := by
  classical
  have htruncSize : truncated.size = active.size :=
    (truncateCld_spec original truncated divisor htruncate).2.1.trans hsize
  let degree := fun column : Fin matrix.size =>
    cldSpiralDegree width (column.val - truncated.size)
  let carry := fun column : Fin matrix.size => Classical.choose (hcongr (degree column))
  let coefficients := fun index : Fin matrix.size =>
    if index.val < truncated.size then (if index ∈ support then (1 : ZZ) else 0)
    else -carry index
  refine ⟨coefficients, ?_, ?_⟩
  · intro index hi
    simp only [coefficients, htruncSize, if_pos hi]
  · have hbound := hbasis.truncated_indicator_short (modulus : Int) divisor height support
      (fun row column => Generated.StrictRecombine.sparseCoeff original[row.val]! (degree column))
      (fun column => actual.coeff (degree column)) carry
      (by simpa only [htruncSize] using hcard)
      (by simpa only [htruncSize] using hsupport)
      (by exact_mod_cast hm) hd hh
      (fun column _ => hactual (degree column))
      (fun column _ row hr => (hcld row.val (hsupport row hr)).sparseCoeff_abs_le (degree column))
      (fun column _ => Classical.choose_spec (hcongr (degree column)))
      (fun column _ row hr => truncateCld_sparseCoeff original truncated divisor htruncate
        row.val (by rw [hsize]; exact hsupport row hr)
        (hcld row.val (hsupport row hr)).canonical (degree column))
    simpa only [coefficients, htruncSize] using hbound

open CLPoly.Math

/-- One Hensel witness supplies both the true modular support and the short
vector. Keeping the support equation in the result is essential for the
later candidate irreducibility argument. -/
theorem hensel_factor_truncated_short
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {aTarget : Int32} {output : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hfactors : ∀ factor ∈ selection.factors.toList,
      SparsePolyZp.Canonical selection.prime.toNat factor)
    (hleadingSemantic : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) =
        (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f
      selection.factors selection.prime aTarget output)
    (original truncated : Array SparsePolyZZ)
    (modulus : Nat) (hmod : output.2 = (modulus : Int))
    (divisor height scale : ZZ) (width columns : Nat)
    (matrix : Generated.StrictRecombine.LLLMatrix)
    (hsize : original.size = output.1.size)
    (hcld : ∀ index, index < output.1.size →
      ModularCldValue f output.1[index]! original[index]! modulus)
    (htruncate : Generated.StrictVanHoeij.truncateCld original divisor = .ok truncated)
    (hbasis : CanonicalCldBasis truncated width scale
      ((modulus : Int) / divisor) columns matrix)
    (factor complement : Polynomial Int)
    (hfactor : SparsePolyZZ.toPoly f = factor * complement)
    (hm : 0 < modulus) (hd : 0 < divisor) (hh : height < (modulus : Int))
    (hdegree : 0 < factor.natDegree)
    (hheight : height = ((SparsePolyZZ.toPoly f).natDegree *
      2 ^ ((SparsePolyZZ.toPoly f).natDegree - 1) *
      (Finset.range ((SparsePolyZZ.toPoly f).natDegree + 1)).sum
        (fun k => ((SparsePolyZZ.toPoly f).coeff k).natAbs) : Nat)) :
    ∃ (indices : Array Nat) (coefficients : Fin matrix.size → ZZ),
      LegalCombination output.1.size indices.size indices ∧
      Polynomial.map (Int.castRingHom (ZMod modulus)) factor =
        Polynomial.C (factor.leadingCoeff : ZMod modulus) *
          (indices.toList.map fun index => StrictHensel.toPolyMod modulus output.1[index]!).prod ∧
      (∀ index, index.val < output.1.size →
        coefficients index = if index.val ∈ indices.toList then 1 else 0) ∧
      dotProduct
        (fun index => ((Matrix.vecMul coefficients
          (basisPrefixMatrix matrix matrix.size matrix.size) index : ZZ) : QQ))
        (fun index => ((Matrix.vecMul coefficients
          (basisPrefixMatrix matrix matrix.size matrix.size) index : ZZ) : QQ)) ≤
        ((output.1.size : ZZ) * scale * scale + (columns : ZZ) *
          ((height + divisor - 1) / divisor + 2 * (output.1.size : Int)) *
          ((height + divisor - 1) / divisor + 2 * (output.1.size : Int)) : ZZ) := by
  classical
  have hfNe : SparsePolyZZ.toPoly f ≠ 0 := by
    intro hz
    have hl := hselection.goodPrime.lc_nonzero
    simp [hz] at hl
  have hmulNe : factor * complement ≠ 0 := by rwa [← hfactor]
  have hactual : ∀ degree,
      |(complement * Polynomial.derivative factor).coeff degree| ≤ height := by
    intro degree
    have hb := CLPoly.VanHoeij.integer_cld_coeff_bound factor complement
      (left_ne_zero_of_mul hmulNe) (right_ne_zero_of_mul hmulNe) hdegree degree
    have hsource : complement * factor = SparsePolyZZ.toPoly f := by
      rw [mul_comm, ← hfactor]
    rw [hsource] at hb
    rw [hheight, ← Int.natCast_natAbs]
    exact_mod_cast hb
  obtain ⟨indices, exponent, hlegal, he, hm', hsupport⟩ :=
    StrictFactorZZ.hensel_divisor_physical_support hcount hfactors hleadingSemantic
      hselection hentry factor complement hfactor
  have hpow : modulus = selection.prime.toNat ^ exponent := by
    exact_mod_cast hmod.symm.trans hm'
  rw [← hpow] at hsupport
  have hindices : ∀ index ∈ indices.toList, index < output.1.size := by
    intro index hi
    obtain ⟨position, hp, heq⟩ := List.mem_iff_getElem.mp hi
    have hp' : position < indices.size := by simpa using hp
    have hb := hlegal.2.2 position hp'
    simpa only [getElem!_pos indices position hp', ← Array.getElem_toList hp', heq] using hb
  have hmonic : ∀ index, index < output.1.size →
      (StrictHensel.toPolyMod modulus output.1[index]!).Monic := by
    intro index hi
    rw [getElem!_pos output.1 index hi]
    exact hentry.outputToPolyModMonic modulus index hi
  have hcongr := fun degree => modularCld_factor_coeff_carry
    f output.1 original modulus hcld hmonic indices.toList hindices factor complement _
      (by simpa only [mul_comm] using hfactor) hsupport degree
  have htruncSize : truncated.size = output.1.size :=
    (truncateCld_spec original truncated divisor htruncate).2.1.trans hsize
  have hrows : output.1.size ≤ matrix.size := by
    rw [hbasis.size_eq, htruncSize]
    omega
  obtain ⟨hcard, hbounded, hsum⟩ := legal_support_lattice_data indices output.1.size
    matrix.size hlegal hrows original (complement * Polynomial.derivative factor) modulus hcongr
  obtain ⟨coefficients, hcoordinates, hnorm⟩ := physical_cld_truncated_short
    f output.1 original truncated modulus divisor height scale width columns matrix
    hsize hcld htruncate hbasis
    (Finset.univ.filter fun row : Fin matrix.size => row.val ∈ indices.toList)
    hcard hbounded (complement * Polynomial.derivative factor) hsum hm hd hh hactual
  refine ⟨indices, coefficients, hlegal, hsupport, ?_, hnorm⟩
  intro index hi
  simpa only [Finset.mem_filter, Finset.mem_univ, true_and] using hcoordinates index hi

end Refinement.StrictRecombine
