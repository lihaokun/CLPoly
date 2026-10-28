import CLPoly.Refinement.VanHoeijValidationSupport

set_option autoImplicit false

namespace Refinement.StrictRecombine

/-- Valid physical indices and a nonzero modulus make the actual product
and degree loop return normally. -/
theorem partitionProductLoop_succeeds
    (group : Array Int32) (lifted : Array SparsePolyZZ) (modulus : ZZ)
    (index : Nat) (product : SparsePolyZZ) (degree : Int64)
    (hm : modulus ≠ 0)
    (hvalid : CandidateIndicesValid group (Array.replicate lifted.size false)) :
    ∃ output finalDegree, Generated.StrictVanHoeij.partitionProductLoop group lifted modulus
      index product degree = .ok (output, finalDegree) := by
  induction hmeasure : group.size - index using Nat.strong_induction_on
      generalizing index product degree with
  | h measure ih =>
    rw [Generated.StrictVanHoeij.partitionProductLoop]
    split
    next hi =>
      have hv := hvalid index hi
      have hn : ¬ group[index] < 0 := by
        have h := hv.1
        rw [Int32.le_iff_toInt_le] at h
        rw [Int32.lt_iff_toInt_lt]
        omega
      have hp : group[index].toInt64.toNat < lifted.size := by simpa using hv.2
      dsimp only
      rw [if_neg hn, dif_pos hp]
      obtain ⟨next, hnext⟩ := multiplyNormalizeModRaw_complete product
        lifted[group[index].toInt64.toNat] modulus hm
      simp only [hnext, bind, Except.bind]
      exact ih (group.size - (index + 1)) (by omega) (index + 1) next
        (degree + get_deg lifted[group[index].toInt64.toNat]) rfl
    next hi => exact ⟨product, degree, rfl⟩

/-- Literal machine addition is retained in this exact execution equation. -/
theorem partitionProductLoop_degree_fold
    (group : Array Int32) (lifted : Array SparsePolyZZ) (modulus : ZZ)
    (index : Nat) (product output : SparsePolyZZ) (degree finalDegree : Int64)
    (hrun : Generated.StrictVanHoeij.partitionProductLoop group lifted modulus index
      product degree = .ok (output, finalDegree)) :
    finalDegree = (group.toList.drop index).foldl
      (fun acc source => acc + get_deg lifted[source.toInt64.toNat]!) degree := by
  induction hmeasure : group.size - index using Nat.strong_induction_on
      generalizing index product degree with
  | h measure ih =>
    rw [Generated.StrictVanHoeij.partitionProductLoop] at hrun
    split at hrun
    next hi =>
      dsimp only at hrun
      split at hrun
      next hn => simp at hrun
      next hn =>
        split at hrun
        next hp =>
          cases hm : Generated.StrictRecombine.multiplyNormalizeModRaw product
              lifted[group[index].toInt64.toNat] modulus with
          | error fault => simp [hm, bind, Except.bind] at hrun
          | ok next =>
            simp only [hm, bind, Except.bind] at hrun
            have ht := ih (group.size - (index + 1)) (by omega) (index + 1) next
              (degree + get_deg lifted[group[index].toInt64.toNat]) hrun rfl
            have hdrop : group.toList.drop index = group[index] :: group.toList.drop (index + 1) := by
              simpa using List.drop_eq_getElem_cons (l := group.toList) (i := index) (by simpa using hi)
            simpa only [hdrop, List.foldl_cons,
              getElem!_pos lifted group[index].toInt64.toNat hp] using ht
        next hp => simp at hrun
    next hi =>
      have he := congrArg Prod.snd (Except.ok.inj hrun)
      have hdrop : group.toList.drop index = [] :=
        List.drop_eq_nil_iff.mpr (by simpa using Nat.le_of_not_gt hi)
      simpa only [hdrop, List.foldl_nil] using he.symm

theorem canonical_get_deg_eq_ofNat (f : SparsePolyZZ)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f) (hne : 0 < f.size) :
    get_deg f = Int64.ofNat (SparsePolyZZ.toPoly f).natDegree := by
  have hempty : f.isEmpty = false := by simp [Array.isEmpty, Nat.ne_of_gt hne]
  unfold get_deg
  rw [hempty]
  simp only [Bool.false_eq_true, ↓reduceIte, getElem!_pos f 0 hne]
  rw [sparsePolyZZ_natDegree_eq_head f hcanonical hne]
  rfl

theorem int64_degree_fold_ofNat {α : Type*} (entries : List α) (degree : α → Nat) (initial : Nat) :
    entries.foldl (fun acc entry => acc + Int64.ofNat (degree entry)) (Int64.ofNat initial) =
      Int64.ofNat (initial + (entries.map degree).sum) := by
  induction entries generalizing initial with
  | nil => simp
  | cons head tail ih =>
    simp only [List.foldl_cons, ← Int64.ofNat_add]
    rw [ih]
    simp only [List.map_cons, List.sum_cons, Nat.add_assoc]

/-- The degree counter equals the natural degree sum converted once to
Int64. Order comparisons require a separate no-overflow bound. -/
theorem partitionProductLoop_degree_ofNat
    (group : Array Int32) (lifted : Array SparsePolyZZ) (modulus : ZZ)
    (product output : SparsePolyZZ) (finalDegree : Int64)
    (hvalid : CandidateIndicesValid group (Array.replicate lifted.size false))
    (hcanonical : ∀ index, index < lifted.size → StrictPolynomialMod.SparsePolyZZCanonical lifted[index]!)
    (hnonempty : ∀ index, index < lifted.size → 0 < lifted[index]!.size)
    (hrun : Generated.StrictVanHoeij.partitionProductLoop group lifted modulus 0
      product 0 = .ok (output, finalDegree)) :
    finalDegree = Int64.ofNat
      ((group.toList.map fun source => (SparsePolyZZ.toPoly lifted[source.toInt64.toNat]!).natDegree).sum) := by
  rw [partitionProductLoop_degree_fold group lifted modulus 0 product output 0 finalDegree hrun]
  simp only [List.drop_zero]
  have he : ∀ source ∈ group.toList,
      get_deg lifted[source.toInt64.toNat]! =
        Int64.ofNat (SparsePolyZZ.toPoly lifted[source.toInt64.toNat]!).natDegree := by
    intro source hs
    obtain ⟨position, hp, hv⟩ := List.mem_iff_getElem.mp hs
    have hp' : position < group.size := by simpa using hp
    have hb : source.toInt64.toNat < lifted.size := by
      have h := (hvalid position hp').2
      simpa only [Array.size_replicate, ← Array.getElem_toList hp', hv] using h
    exact canonical_get_deg_eq_ofNat _ (hcanonical _ hb) (hnonempty _ hb)
  have hf : ∀ (entries : List Int32), entries ⊆ group.toList → ∀ initial : Int64,
      entries.foldl (fun acc source => acc + get_deg lifted[source.toInt64.toNat]!) initial =
      entries.foldl (fun acc source => acc +
        Int64.ofNat (SparsePolyZZ.toPoly lifted[source.toInt64.toNat]!).natDegree) initial := by
    intro entries
    induction entries with
    | nil => intros; rfl
    | cons head tail ih =>
      intro hsub initial
      simp only [List.foldl_cons, he head (hsub List.mem_cons_self)]
      exact ih (fun x hx => hsub (List.mem_cons_of_mem head hx)) _
  rw [hf group.toList (List.Subset.refl _)]
  simpa only [Nat.zero_add] using int64_degree_fold_ofNat group.toList
    (fun source => (SparsePolyZZ.toPoly lifted[source.toInt64.toNat]!).natDegree) 0

end Refinement.StrictRecombine
