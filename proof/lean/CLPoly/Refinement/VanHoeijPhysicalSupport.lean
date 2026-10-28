import CLPoly.Refinement.VanHoeijSupportMultiplicity

set_option autoImplicit false

namespace Refinement.StrictRecombine
open CLPoly.Math

def candidatePhysicalSupport (count : Nat) (group : Array Int32)
    (hvalid : CandidateIndicesValid group (Array.replicate count false)) : List (Fin count) :=
  List.ofFn fun i : Fin group.size => ⟨group[i.val].toInt64.toNat,
    by simpa only [Array.size_replicate] using (hvalid i.val i.isLt).2⟩

theorem candidatePhysicalSupport_product
    (prime : Nat) (lifted : Array SparsePolyZZ) (group : Array Int32)
    (hvalid : CandidateIndicesValid group (Array.replicate lifted.size false)) :
    ((candidatePhysicalSupport lifted.size group hvalid).map
      (fun i => StrictHensel.toPolyMod prime lifted[i.val])).prod =
      SelectedProductMod prime group lifted 0 := by
  have hlist : group.toList = List.ofFn (fun i : Fin group.size => group[i.val]) := by
    simpa only [Array.toList_ofFn] using
      (congrArg Array.toList (Array.ofFn_getElem (xs := group))).symm
  unfold SelectedProductMod candidatePhysicalSupport
  rw [List.drop_zero, hlist, List.map_ofFn, List.map_ofFn]
  congr 1
  apply congrArg List.ofFn
  funext i
  have hi : group[i.val].toInt64.toNat < lifted.size := by
    simpa only [Array.size_replicate] using (hvalid i.val i.isLt).2
  change StrictHensel.toPolyMod prime lifted[group[i.val].toInt64.toNat] =
    StrictHensel.toPolyMod prime lifted[group[i.val].toInt64.toNat]!
  rw [getElem!_pos lifted group[i.val].toInt64.toNat hi]

theorem candidatePhysicalSupport_nodup
    (prime : Nat) [Fact (Nat.Prime prime)]
    (lifted : Array SparsePolyZZ) (group : Array Int32)
    (hvalid : CandidateIndicesValid group (Array.replicate lifted.size false))
    (hirr : ∀ i (hi : i < lifted.size),
      Irreducible (StrictHensel.toPolyMod prime lifted[i]))
    (factor : SparsePolyZZ) (hsquarefree : Squarefree (StrictHensel.toPolyMod prime factor))
    (hsupport : Associated (StrictHensel.toPolyMod prime factor)
      (SelectedProductMod prime group lifted 0)) :
    (candidatePhysicalSupport lifted.size group hvalid).Nodup := by
  apply associated_support_nodup
    (fun i : Fin lifted.size => StrictHensel.toPolyMod prime lifted[i.val])
    (fun i => (hirr i.val i.isLt).not_isUnit) _ _ hsquarefree
  rw [candidatePhysicalSupport_product]
  exact hsupport

theorem hensel_validated_factor_physical_support
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
    (hdegree : (SparsePolyZZ.toPoly f).natDegree < 2 ^ 63)
    (hextract : Generated.StrictRecombine.extractCandidates rows transform output.1.size =
      .ok groups)
    (hrun : Generated.StrictVanHoeij.validatePartition f output.1 modulus groups before =
      .ok (true, after)) :
    ∀ factor ∈ after.factors.toList, ∃ group ∈ groups.toList,
      ∃ hvalid : CandidateIndicesValid group (Array.replicate output.1.size false),
        (candidatePhysicalSupport output.1.size group hvalid).Nodup ∧
        Associated (StrictHensel.toPolyMod selection.prime.toNat factor)
          (((candidatePhysicalSupport output.1.size group hvalid).map
            (fun i => StrictHensel.toPolyMod selection.prime.toNat output.1[i.val])).prod) := by
  have hs := hensel_extract_validatePartition_support hcount hfactors hleadingSemantic
    hselection hentry modulus hmod groups rows transform before after hfits hcanonical
    hextract hrun
  have hsq := validatePartition_factor_mod_squarefree f output.1 modulus selection.prime.toNat
    groups before after hcanonical hdegree hselection.goodPrime.sqfree hrun
  have hirr := StrictFactorZZ.selectionHenselFactors_mod_irreducible hcount hfactors
    hleadingSemantic hselection hentry
  intro factor hf
  obtain ⟨group, hg, ha⟩ := hs factor hf
  have hn : rows.isEmpty = false := by
    cases he : rows.isEmpty with
    | false => rfl
    | true =>
      have hgroups : groups = #[] := by
        simpa [Generated.StrictRecombine.extractCandidates, he] using hextract.symm
      simp [hgroups] at hg
  obtain ⟨i, hi, heq⟩ := List.mem_iff_getElem.mp hg
  have hi' : i < groups.size := by simpa using hi
  have heq' : groups[i] = group := by simpa only [Array.getElem_toList] using heq
  have hv : CandidateIndicesValid group (Array.replicate output.1.size false) := by
    rw [← heq']
    exact extractCandidates_indices_valid rows transform output.1.size groups hfits hn hextract i hi'
  refine ⟨group, List.mem_iff_getElem.mpr ⟨i, hi, heq⟩, hv,
    candidatePhysicalSupport_nodup selection.prime.toNat output.1 group hv hirr
      factor (hsq factor hf) ha, ?_⟩
  rw [candidatePhysicalSupport_product]
  exact ha

end Refinement.StrictRecombine
