import CLPoly.Refinement.VanHoeijControllerArithmetic

set_option autoImplicit false

namespace Refinement.StrictRecombine

theorem singletonGroupsLoop_size (count index : Int32)
    (groups output : Array (Array Int32))
    (hrun : Generated.StrictVanHoeij.singletonGroupsLoop count index groups = .ok output) :
    output.size = groups.size := by
  induction hmeasure : (count.toInt - index.toInt).toNat using Nat.strong_induction_on
      generalizing index groups output with
  | h measure ih =>
    rw [Generated.StrictVanHoeij.singletonGroupsLoop] at hrun
    split at hrun
    next hlt =>
      split at hrun
      next hn =>
        split at hrun
        next hi =>
          have hdecrease : (count.toInt - (index + 1).toInt).toNat < measure := by
            rw [Generated.StrictVanHoeij.int32_increment_toInt_of_lt index count hlt]
            have hl := Int32.lt_iff_toInt_lt.mp hlt
            omega
          have hs := ih _ hdecrease (index + 1) _ output hrun rfl
          simpa only [Array.size_set] using hs
        next hout => simp at hrun
      next hnegative => simp at hrun
    next hdone =>
      exact congrArg Array.size (Except.ok.inj hrun.symm)

theorem singletonGroups_size (count : Int32) (output : Array (Array Int32))
    (hrun : Generated.StrictVanHoeij.singletonGroups count = .ok output) :
    output.size = count.toInt.toNat := by
  unfold Generated.StrictVanHoeij.singletonGroups at hrun
  split at hrun
  · simpa only [Array.size_replicate] using singletonGroupsLoop_size count 0 _ output hrun
  · simp at hrun

theorem singletonGroupsLoop_entries (count index : Int32)
    (groups output : Array (Array Int32))
    (hnonnegative : 0 ≤ index) (hle : index.toInt ≤ count.toInt)
    (hsize : groups.size = count.toInt.toNat)
    (hprefix : ∀ p (hp : p < groups.size), groups[p] =
      if p < index.toInt.toNat then #[Int32.ofNat p] else #[])
    (hrun : Generated.StrictVanHoeij.singletonGroupsLoop count index groups = .ok output) :
    ∀ p (hp : p < output.size), output[p] = #[Int32.ofNat p] := by
  induction hmeasure : (count.toInt - index.toInt).toNat using Nat.strong_induction_on
      generalizing index groups output with
  | h measure ih =>
    rw [Generated.StrictVanHoeij.singletonGroupsLoop] at hrun
    split at hrun
    next hlt =>
      split at hrun
      next hi =>
        have hn := hnonnegative
        have hl := Int32.lt_iff_toInt_lt.mp hlt
        have hz := Int32.le_iff_toInt_le.mp hn
        change (0 : Int) ≤ index.toInt at hz
        have hincr := Generated.StrictVanHoeij.int32_increment_toInt_of_lt index count hlt
        have hnat : index.toInt64.toNat = index.toInt.toNat := by
          change index.toInt64.toNatClampNeg = index.toInt.toNat
          simp only [Int32.toNatClampNeg_toInt64, Int32.toNat_toInt]
        have hnextNat : (index + 1).toInt.toNat = index.toInt.toNat + 1 := by
          rw [hincr]
          omega
        have hround : Int32.ofNat index.toInt.toNat = index :=
          Int32.ofNat_toNatClampNeg index hn
        have hnextNonnegative : 0 ≤ index + 1 := by
          rw [Int32.le_iff_toInt_le, hincr]
          change (0 : Int) ≤ index.toInt + 1
          omega
        have hnextLe : (index + 1).toInt ≤ count.toInt := by rw [hincr]; omega
        have hnextPrefix : ∀ p (hp : p < (groups.set index.toInt64.toNat
            (groups[index.toInt64.toNat].push index)).size),
            (groups.set index.toInt64.toNat (groups[index.toInt64.toNat].push index))[p] =
              if p < (index + 1).toInt.toNat then #[Int32.ofNat p] else #[] := by
          intro p hp
          have hp' : p < groups.size := by simpa only [Array.size_set] using hp
          rw [Array.getElem_set]
          by_cases heq : index.toInt64.toNat = p
          · subst p
            have hold : groups[index.toInt64.toNat] = #[] := by
              simpa only [hnat, Nat.lt_irrefl, ↓reduceIte] using hprefix _ hi
            rw [if_pos rfl, hold]
            simp only [Array.push_empty, hnextNat, hnat, Nat.lt_succ_self, ↓reduceIte, hround]
          · rw [if_neg heq, hprefix p hp']
            have hequiv : p < index.toInt.toNat ↔ p < (index + 1).toInt.toNat := by
              rw [hnextNat]
              rw [hnat] at heq
              omega
            simp only [hequiv]
        have hdecrease : (count.toInt - (index + 1).toInt).toNat < measure := by
          rw [hincr]
          omega
        exact ih _ hdecrease (index + 1) _ output hnextNonnegative hnextLe
          (by simpa only [Array.size_set] using hsize) hnextPrefix hrun rfl
      next hout => simp at hrun
    next hdone =>
      have heq := Except.ok.inj hrun
      subst output
      intro p hp
      rw [hprefix p hp, if_pos]
      have hz := Int32.le_iff_toInt_le.mp hnonnegative
      rw [Int32.lt_iff_toInt_lt] at hdone
      omega

theorem singletonGroups_entries (count : Int32) (output : Array (Array Int32))
    (hrun : Generated.StrictVanHoeij.singletonGroups count = .ok output) :
    ∀ p (hp : p < output.size), output[p] = #[Int32.ofNat p] := by
  unfold Generated.StrictVanHoeij.singletonGroups at hrun
  split at hrun
  next hn =>
    apply singletonGroupsLoop_entries count 0 _ output (by decide)
      (Int32.le_iff_toInt_le.mp hn) (by simp) ?_ hrun
    intro p hp
    simp only [Array.getElem_replicate, Int32.toInt_zero, Int.toNat_zero, Nat.not_lt_zero,
      ↓reduceIte]
  next hnegative => simp at hrun

theorem singletonGroupsLoop_succeeds (count index : Int32)
    (groups : Array (Array Int32))
    (hnonnegative : 0 ≤ index) (hle : index.toInt ≤ count.toInt)
    (hsize : groups.size = count.toInt.toNat) :
    ∃ output, Generated.StrictVanHoeij.singletonGroupsLoop count index groups = .ok output := by
  induction hmeasure : (count.toInt - index.toInt).toNat using Nat.strong_induction_on
      generalizing index groups with
  | h measure ih =>
    by_cases hlt : index < count
    · have hl := Int32.lt_iff_toInt_lt.mp hlt
      have hz := Int32.le_iff_toInt_le.mp hnonnegative
      change (0 : Int) ≤ index.toInt at hz
      have hnat : index.toInt64.toNat = index.toInt.toNat := by
        change index.toInt64.toNatClampNeg = index.toInt.toNat
        simp only [Int32.toNatClampNeg_toInt64, Int32.toNat_toInt]
      have hi : index.toInt64.toNat < groups.size := by rw [hnat, hsize]; omega
      have hincr := Generated.StrictVanHoeij.int32_increment_toInt_of_lt index count hlt
      have hnextNonnegative : 0 ≤ index + 1 := by
        rw [Int32.le_iff_toInt_le, hincr]
        change (0 : Int) ≤ index.toInt + 1
        omega
      have hnextLe : (index + 1).toInt ≤ count.toInt := by rw [hincr]; omega
      have hdecrease : (count.toInt - (index + 1).toInt).toNat < measure := by
        rw [hincr]
        omega
      obtain ⟨output, hout⟩ := ih _ hdecrease (index + 1)
        (groups.set index.toInt64.toNat (groups[index.toInt64.toNat].push index))
        hnextNonnegative hnextLe (by simpa only [Array.size_set] using hsize) rfl
      refine ⟨output, ?_⟩
      rw [Generated.StrictVanHoeij.singletonGroupsLoop, dif_pos hlt,
        if_pos hnonnegative, dif_pos hi]
      exact hout
    · exact ⟨groups, by rw [Generated.StrictVanHoeij.singletonGroupsLoop, dif_neg hlt]⟩

theorem singletonGroups_succeeds (count : Int32) (hnonnegative : 0 ≤ count) :
    ∃ output, Generated.StrictVanHoeij.singletonGroups count = .ok output ∧
      output.size = count.toInt.toNat ∧
      ∀ p (hp : p < output.size), output[p] = #[Int32.ofNat p] := by
  obtain ⟨output, hout⟩ := singletonGroupsLoop_succeeds count 0
    (Array.replicate count.toInt.toNat #[]) (by decide)
    (Int32.le_iff_toInt_le.mp hnonnegative) (by simp)
  have hrun : Generated.StrictVanHoeij.singletonGroups count = .ok output := by
    rw [Generated.StrictVanHoeij.singletonGroups, if_pos hnonnegative]
    exact hout
  exact ⟨output, hrun, singletonGroups_size count output hrun,
    singletonGroups_entries count output hrun⟩

end Refinement.StrictRecombine
