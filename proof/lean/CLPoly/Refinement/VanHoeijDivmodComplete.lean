import CLPoly.Refinement.VanHoeijDivmodReduced

set_option autoImplicit false

namespace Refinement.StrictHensel
open _root_.Generated.StrictHensel

theorem canonical_nonempty_mod_degree (poly : SparsePolyZZ) (m : Nat)
    (hc : StrictPolynomialMod.SparsePolyZZCanonical poly) (hn : 0 < poly.size)
    (hcast : (poly[0].2 : ZMod m) ≠ 0) :
    toPolyMod m poly ≠ 0 ∧ (toPolyMod m poly).natDegree = poly[0].1.deg := by
  have hlc : ((SparsePolyZZ.toPoly poly).leadingCoeff : ZMod m) ≠ 0 := by
    rw [henselSparsePolyZZ_leadingCoeff_eq_head poly hc hn]
    exact hcast
  have hleading : (toPolyMod m poly).leadingCoeff = (poly[0].2 : ZMod m) := by
    unfold toPolyMod
    rw [Polynomial.leadingCoeff_map_of_leadingCoeff_ne_zero _ hlc,
      henselSparsePolyZZ_leadingCoeff_eq_head poly hc hn]
    rfl
  constructor
  · intro hz
    apply hcast
    rw [hz, Polynomial.leadingCoeff_zero] at hleading
    exact hleading.symm
  · unfold toPolyMod
    rw [Polynomial.natDegree_map_of_leadingCoeff_ne_zero _ hlc,
      StrictRecombine.sparsePolyZZ_natDegree_eq_head poly hc hn]

theorem reduced_canonical_nonempty_mod_degree (poly : SparsePolyZZ) (m : Nat)
    (hc : StrictPolynomialMod.SparsePolyZZCanonical poly) (hn : 0 < poly.size)
    (hr : CoefficientsReduced m poly) :
    toPolyMod m poly ≠ 0 ∧ (toPolyMod m poly).natDegree = poly[0].1.deg := by
  apply canonical_nonempty_mod_degree poly m hc hn
  have hb := hr poly[0] (Array.getElem_mem_toList hn)
  have hne := hc.2 poly[0] (Array.getElem_mem_toList hn)
  intro hz
  have hd := (ZMod.intCast_zmod_eq_zero_iff_dvd poly[0].2 m).mp hz
  have hm := Int.emod_eq_zero_of_dvd hd
  have hs := Int.emod_eq_of_lt hb.1 hb.2
  exact hne (hs.symm.trans hm)

theorem monic_divmod_remainder_empty_of_dvd
    (termination : DivmodTermination) (f g q r : SparsePolyZZ) (m : Nat) (hm : 1 < m)
    (hf : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hg : StrictPolynomialMod.SparsePolyZZCanonical g)
    (hsize : 0 < g.size) (hone : g[0]!.2 = 1)
    (hfBound : DegreesBound f) (hgBound : DegreesBound g)
    (hdivides : toPolyMod m g ∣ toPolyMod m f)
    (hrun : __upoly_divmod_mod_raw_ir termination f g (m : Int) = .ok (q, r)) :
    r = #[] := by
  letI : Fact (1 < m) := ⟨hm⟩
  letI : NeZero m := ⟨by omega⟩
  have hhead := canonical_headDominates g hsize hg
  have hdegree : g[0]!.1.deg < 2 ^ 63 := by
    apply hgBound g[0]!
    rw [getElem!_pos g 0 hsize]
    exact Array.getElem_mem_toList hsize
  have hinvert : (ZZ.invert 0 g[0]!.2 (m : Int)).1 = true := by
    rw [hone]; exact StrictRecombine.invert_one_success m hm
  have hsemantic := __upoly_divmod_mod_raw_ir_refines_of_run termination f g q r m
    hsize hdegree hhead hfBound hinvert hrun
  have hrc := __upoly_divmod_mod_raw_ir_remainder_canonical_of_run termination f g q r (m : Int) hf hg hrun
  have hrb := __upoly_divmod_mod_raw_ir_remainder_degreesBound_of_run termination f g q r m
    hsize hdegree hhead hfBound hrun
  have hbelow := __upoly_divmod_mod_raw_ir_remainder_terms_below_of_run termination f g q r
    (m : Int) hrun hsize hdegree hrc hrb
  have hred := divmod_raw_remainder_coefficientsReduced termination f g q r m (by omega) hrun
  by_contra hnotEmpty
  have hrsize : 0 < r.size := by
    by_contra h
    exact hnotEmpty (Array.size_eq_zero_iff.mp (by omega))
  obtain ⟨hrne, hrdegree⟩ := reduced_canonical_nonempty_mod_degree r m hrc hrsize hred
  have hgcast : (g[0].2 : ZMod m) ≠ 0 := by
    have h : g[0].2 = 1 := by simpa only [getElem!_pos g 0 hsize] using hone
    rw [h]; simp
  have hgdegree := (canonical_nonempty_mod_degree g m hg hsize hgcast).2
  have hrem : toPolyMod m r = toPolyMod m f - toPolyMod m q * toPolyMod m g := by
    rw [← hsemantic]; ring
  have hdvd : toPolyMod m g ∣ toPolyMod m r := by
    rw [hrem]
    exact dvd_sub hdivides (dvd_mul_left _ _)
  have hgmonic : (toPolyMod m g).Monic := by
    unfold Polynomial.Monic toPolyMod
    rw [Polynomial.leadingCoeff_map_of_leadingCoeff_ne_zero _ (by
      rw [henselSparsePolyZZ_leadingCoeff_eq_head g hg hsize]
      exact hgcast), henselSparsePolyZZ_leadingCoeff_eq_head g hg hsize]
    have h : g[0].2 = 1 := by simpa only [getElem!_pos g 0 hsize] using hone
    simp [h]
  obtain ⟨s, hs⟩ := hdvd
  have hsne : s ≠ 0 := by
    intro hz
    apply hrne
    simpa [hz] using hs
  have hle : (toPolyMod m g).natDegree ≤ (toPolyMod m r).natDegree := by
    rw [hs, hgmonic.natDegree_mul' hsne]
    omega
  rw [hgdegree, hrdegree] at hle
  have hlt := hbelow r[0] (Array.getElem_mem_toList hrsize)
  rw [getElem!_pos g 0 hsize] at hlt
  omega

/-- Concrete modular division by a physical monic divisor of the source
executes normally and produces an actually empty remainder array. -/
theorem monic_divmod_complete_of_dvd (f g : SparsePolyZZ) (m : Nat) (hm : 1 < m)
    (hf : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hg : StrictPolynomialMod.SparsePolyZZCanonical g)
    (hsize : 0 < g.size) (hone : g[0]!.2 = 1)
    (hfBound : DegreesBound f) (hgBound : DegreesBound g)
    (hdivides : toPolyMod m g ∣ toPolyMod m f) :
    ∃ q, __upoly_divmod_mod_raw_ir concreteDivmodTermination f g (m : Int) = .ok (q, #[]) := by
  have hvalid : concreteDivmodTermination.inputValid g (modCoeffOutput f (m : Int)) :=
    (concreteDivmodInputValid_iff g _).mpr ⟨hsize, hg, hgBound,
      modCoeffOutput_canonical f (m : Int) hf, modCoeffOutput_degreesBound f (m : Int) hfBound⟩
  have hinvert : (ZZ.invert 0 g[0]!.2 (m : Int)).1 = true := by
    rw [hone]; exact StrictRecombine.invert_one_success m hm
  obtain ⟨q, r, hrun⟩ := __upoly_divmod_mod_raw_ir_terminates concreteDivmodTermination f g (m : Int)
    (by simp [Array.isEmpty, Nat.ne_of_gt hsize]) hvalid hinvert
  have he := monic_divmod_remainder_empty_of_dvd concreteDivmodTermination f g q r m hm
    hf hg hsize hone hfBound hgBound hdivides hrun
  exact ⟨q, by simpa only [he] using hrun⟩

end Refinement.StrictHensel
