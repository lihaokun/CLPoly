import CLPoly.Refinement.VanHoeijHeight

set_option autoImplicit false

namespace Refinement.StrictRecombine

theorem invert_one_success (modulus : Nat) (hm : 1 < modulus) :
    (ZZ.invert 0 1 (modulus : Int)).1 = true := by
  have hmInt : (1 : Int) < modulus := by exact_mod_cast hm
  have hone : (1 : Int) % (modulus : Int) = 1 := Int.emod_eq_of_lt (by omega) hmInt
  simp only [ZZ.invert, ZZ.invertImpl, not_le.mpr hmInt, ↓reduceIte,
    Int.add_emod_right, hone, Int.natAbs_one, Int.natAbs_natCast]
  have heg : Nat.extGcd 1 modulus = (1, 1, 0) := by
    have hz : Nat.extGcd 1 0 = (1, 1, 0) := by rw [Nat.extGcd]; simp
    have hstep : Nat.extGcd modulus 1 = (1, 0, 1) := by
      rw [Nat.extGcd]
      simp only [Nat.one_ne_zero, ↓reduceDIte, Nat.mod_one, Nat.div_one, hz]
      norm_num
    rw [Nat.extGcd]
    simp [Nat.ne_of_gt (by omega : 0 < modulus), Nat.mod_eq_of_lt hm,
      Nat.div_eq_of_lt hm, hstep]
  simp [heg]

theorem hensel_cld_run_valid
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {factors : Array SparsePolyZp} {prime : UInt64}
    {aTarget : Int32} {output : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime prime.toNat)]
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f factors prime aTarget output)
    (modulus : Nat) (hm : 1 < modulus)
    (hfBound : StrictHensel.DegreesBound f)
    (hdegree : ∀ index (hi : index < output.1.size),
      output.1[index][0]!.1.deg < 2 ^ 63)
    (cld : Array SparsePolyZZ)
    (hrun : Generated.StrictRecombine.cldPolys ⟨termination⟩ f output.1 modulus = .ok cld) :
    cld.size = output.1.size ∧
      ∀ index, index < output.1.size →
        ModularCldValue f output.1[index]! cld[index]! modulus := by
  have hphysical : ∀ index (hi : index < output.1.size),
      0 < output.1[index].size ∧ output.1[index][0]!.2 = 1 := by
    intro index hi
    obtain ⟨head, tail, hlist⟩ := hentry.outputOneHead index hi
    have hlength := congrArg List.length hlist
    have hsEq : output.1[index].size = tail.length + 1 := by simpa using hlength
    have hs : 0 < output.1[index].size := by omega
    refine ⟨hs, ?_⟩
    have hget := Array.getElem_toList hs
    rw [getElem!_pos output.1[index] 0 hs]
    simpa [hlist] using (congrArg Prod.snd hget).symm
  apply cldPolys_modular_valid ⟨termination⟩ f output.1 cld modulus (by omega) hfBound
    (fun index hi => (hphysical index hi).1) hdegree
  · intro index hi
    exact StrictHensel.canonical_headDominates _ (hphysical index hi).1
      (hentry.outputCanonical index hi)
  · intro index hi
    rw [(hphysical index hi).2]
    exact invert_one_success modulus hm
  · intro index hi
    exact hentry.outputToPolyModMonic modulus index hi
  · exact hrun

end Refinement.StrictRecombine
