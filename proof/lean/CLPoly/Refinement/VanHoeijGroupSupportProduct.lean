import CLPoly.Refinement.VanHoeijGroupMultiplicity

set_option autoImplicit false

namespace Refinement.StrictRecombine

theorem candidateMember_iff_mem_row
    (groups : Array (Array Int32)) (classId : Nat) (hc : classId < groups.size) (value : Int32) :
    CandidateMember groups classId value ↔ value ∈ groups[classId]!.toList := by
  constructor
  · rintro ⟨_, position, hp, hv⟩
    apply List.mem_iff_getElem.mpr
    refine ⟨position, by simpa using hp, ?_⟩
    simpa only [getElem!_pos groups[classId]! position hp, Array.getElem_toList] using hv
  · intro hm
    obtain ⟨position, hp, hv⟩ := List.mem_iff_getElem.mp hm
    have hp' : position < groups[classId]!.size := by simpa using hp
    refine ⟨hc, position, hp', ?_⟩
    simpa only [getElem!_pos groups[classId]! position hp', Array.getElem_toList] using hv

theorem extractCandidates_decoded_member_iff
    (transform : Generated.StrictRecombine.LLLMatrix) (rows : Array Nat)
    (count : Nat) (groups : Array (Array Int32)) (hfits : count ≤ 2 ^ 31)
    (hrun : Generated.StrictRecombine.extractCandidates rows transform count = .ok groups)
    (classId : Nat) (hc : classId < groups.size) (index : Nat) :
    index ∈ groups[classId]!.toList.map (fun value => value.toInt64.toNat) ↔
      index < count ∧ CandidateMember groups classId index.toUInt32.toInt32 := by
  have hrows : rows.isEmpty = false := by
    cases he : rows.isEmpty with
    | false => rfl
    | true =>
      have hg : groups = #[] := by
        simpa [Generated.StrictRecombine.extractCandidates, he] using hrun.symm
      simp [hg] at hc
  constructor
  · intro hi
    obtain ⟨value, hv, he⟩ := List.mem_map.mp hi
    have hm := (candidateMember_iff_mem_row groups classId hc value).mpr hv
    obtain ⟨source, hs, hsource⟩ := extractCandidates_members_sound rows transform count groups
      hrows hrun classId value hm
    rw [hsource, (nat_toUInt32_toInt32_nonnegative_and_toNat source (lt_of_lt_of_le hs hfits)).2] at he
    subst index
    exact ⟨hs, by simpa only [hsource] using hm⟩
  · rintro ⟨hi, hm⟩
    exact List.mem_map.mpr ⟨index.toUInt32.toInt32,
      (candidateMember_iff_mem_row groups classId hc _).mp hm,
      (nat_toUInt32_toInt32_nonnegative_and_toNat index (lt_of_lt_of_le hi hfits)).2⟩

/-- Exact membership plus multiplicity from the executed collector yields
permutation equality with the same legal Hensel support. -/
theorem extractCandidates_exact_support_perm
    (transform : Generated.StrictRecombine.LLLMatrix) (rows : Array Nat)
    (count : Nat) (groups : Array (Array Int32)) (hfits : count ≤ 2 ^ 31)
    (hrun : Generated.StrictRecombine.extractCandidates rows transform count = .ok groups)
    (classId : Nat) (hc : classId < groups.size)
    (indices : Array Nat) (hlegal : LegalCombination count indices.size indices)
    (hmembers : ∀ index, index < count →
      (CandidateMember groups classId index.toUInt32.toInt32 ↔ index ∈ indices.toList)) :
    (groups[classId]!.toList.map (fun value => value.toInt64.toNat)).Perm indices.toList := by
  apply (List.perm_ext_iff_of_nodup
    (extractCandidates_decoded_row_nodup transform rows count groups hfits hrun classId hc)
    (StrictFactorZZ.legalCombination_indices_nodup hlegal)).mpr
  intro index
  rw [extractCandidates_decoded_member_iff transform rows count groups hfits hrun classId hc index]
  constructor
  · rintro ⟨hi, hm⟩
    exact (hmembers index hi).mp hm
  · intro hm
    obtain ⟨position, hp, he⟩ := List.mem_iff_getElem.mp hm
    have hp' : position < indices.size := by simpa using hp
    have hi : index < count := by
      have h := hlegal.2.2 position hp'
      simpa only [getElem!_pos indices position hp', ← Array.getElem_toList hp', he] using h
    exact ⟨hi, (hmembers index hi).mpr hm⟩

theorem extractCandidates_exact_support_product
    (transform : Generated.StrictRecombine.LLLMatrix) (rows : Array Nat)
    (lifted : Array SparsePolyZZ) (groups : Array (Array Int32)) (hfits : lifted.size ≤ 2 ^ 31)
    (hrun : Generated.StrictRecombine.extractCandidates rows transform lifted.size = .ok groups)
    (classId : Nat) (hc : classId < groups.size)
    (indices : Array Nat) (hlegal : LegalCombination lifted.size indices.size indices)
    (hmembers : ∀ index, index < lifted.size →
      (CandidateMember groups classId index.toUInt32.toInt32 ↔ index ∈ indices.toList))
    (modulus : Nat) :
    SelectedProductMod modulus groups[classId]! lifted 0 =
      (indices.toList.map fun index => StrictHensel.toPolyMod modulus lifted[index]!).prod := by
  have hp := extractCandidates_exact_support_perm transform rows lifted.size groups hfits
    hrun classId hc indices hlegal hmembers
  have hprod := (hp.map (fun index => StrictHensel.toPolyMod modulus lifted[index]!)).prod_eq
  simpa only [SelectedProductMod, List.drop_zero, List.map_map, Function.comp_def] using hprod

end Refinement.StrictRecombine
