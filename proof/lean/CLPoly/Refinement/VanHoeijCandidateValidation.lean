import CLPoly.Refinement.VanHoeijValidationProvenance

set_option autoImplicit false

namespace Refinement.StrictRecombine
open CLPoly.Math

theorem extractCandidates_indices_valid
    (rows : Array Nat) (transform : Generated.StrictRecombine.LLLMatrix)
    (count : Nat) (groups : Array (Array Int32))
    (hfits : count ≤ 2 ^ 31) (hnonempty : rows.isEmpty = false)
    (hextract : Generated.StrictRecombine.extractCandidates rows transform count = .ok groups) :
    ∀ i (hi : i < groups.size),
      CandidateIndicesValid groups[i] (Array.replicate count false) := by
  intro i hi position hp
  have hmember : CandidateMember groups i groups[i][position] := by
    refine ⟨hi, position, ?_, ?_⟩
    · simpa only [getElem!_pos groups i hi] using hp
    · simp only [getElem!_pos groups i hi, getElem!_pos groups[i] position hp]
  obtain ⟨source, hsource, heq⟩ := extractCandidates_members_sound rows transform count
    groups hnonempty hextract i groups[i][position] hmember
  have hroundtrip := nat_toUInt32_toInt32_nonnegative_and_toNat source
    (lt_of_lt_of_le hsource hfits)
  rw [heq]
  refine ⟨hroundtrip.1, ?_⟩
  simpa only [hroundtrip.2, Array.size_replicate] using hsource

theorem extractCandidates_validatePartition_support
    (f : SparsePolyZZ) (lifted : Array SparsePolyZZ) (modulus prime : Nat)
    [Fact (Nat.Prime prime)] (groups : Array (Array Int32))
    (rows : Array Nat) (transform : Generated.StrictRecombine.LLLMatrix)
    (before after : Generated.StrictVanHoeij.VanHoeijResult)
    (hm : 0 < modulus) (hdivides : prime ∣ modulus)
    (hfits : lifted.size ≤ 2 ^ 31)
    (hirreducible : ∀ i (hi : i < lifted.size),
      Irreducible (StrictHensel.toPolyMod prime lifted[i]))
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hleading : ((SparsePolyZZ.toPoly f).leadingCoeff : ZMod prime) ≠ 0)
    (hextract : Generated.StrictRecombine.extractCandidates rows transform lifted.size =
      .ok groups)
    (hrun : Generated.StrictVanHoeij.validatePartition f lifted modulus groups before =
      .ok (true, after)) :
    ∀ factor ∈ after.factors.toList, HasPartitionSupport prime lifted groups factor := by
  have hnonempty : rows.isEmpty = false := by
    cases he : rows.isEmpty with
    | false => rfl
    | true =>
      have hgroups : groups = #[] := by
        simpa [Generated.StrictRecombine.extractCandidates, he] using hextract.symm
      subst groups
      simp [Generated.StrictVanHoeij.validatePartition] at hrun
  exact validatePartition_support f lifted modulus prime groups before after hm hdivides
    (extractCandidates_indices_valid rows transform lifted.size groups hfits hnonempty hextract)
    hirreducible hcanonical hleading hrun

theorem hensel_extract_validatePartition_support
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
    (modulus : Nat) (hmod : output.2 = (modulus : Int))
    (groups : Array (Array Int32)) (rows : Array Nat)
    (transform : Generated.StrictRecombine.LLLMatrix)
    (before after : Generated.StrictVanHoeij.VanHoeijResult)
    (hfits : output.1.size ≤ 2 ^ 31)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hextract : Generated.StrictRecombine.extractCandidates rows transform output.1.size =
      .ok groups)
    (hrun : Generated.StrictVanHoeij.validatePartition f output.1 modulus groups before =
      .ok (true, after)) :
    ∀ factor ∈ after.factors.toList,
      HasPartitionSupport selection.prime.toNat output.1 groups factor := by
  obtain ⟨exponent, he, hpower⟩ := hentry.outputModulus_eq_prime_pow
  have hmodulus : modulus = selection.prime.toNat ^ exponent := by
    exact_mod_cast hmod.symm.trans hpower
  have hprime := (Fact.out : Nat.Prime selection.prime.toNat).pos
  have hm : 0 < modulus := by rw [hmodulus]; exact pow_pos hprime _
  have hd : selection.prime.toNat ∣ modulus := by
    rw [hmodulus]
    exact dvd_pow_self _ (Nat.ne_of_gt he)
  exact extractCandidates_validatePartition_support f output.1 modulus selection.prime.toNat
    groups rows transform before after hm hd hfits
    (StrictFactorZZ.selectionHenselFactors_mod_irreducible hcount hfactors hleadingSemantic
      hselection hentry)
    hcanonical hselection.goodPrime.lc_nonzero hextract hrun

end Refinement.StrictRecombine
