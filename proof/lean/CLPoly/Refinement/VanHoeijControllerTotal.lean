import CLPoly.Refinement.VanHoeijBuildExecution
import CLPoly.Refinement.VanHoeijLLLExecutionComplete
import CLPoly.Refinement.VanHoeijGroupingTotal
import CLPoly.Refinement.VanHoeijCandidateValidation
import CLPoly.Refinement.VanHoeijValidationExecution
import CLPoly.Refinement.VanHoeijControllerExecution

set_option autoImplicit false

namespace Refinement.StrictRecombine
open _root_.Generated.StrictVanHoeijController

theorem columnRound_succeeds
    (f : SparsePolyZZ) (lifted cld : Array SparsePolyZZ)
    (originalModulus : Nat) (hm : 0 < originalModulus)
    (modulus scale residual : ZZ) (count target : Int32)
    (countSize : count.toInt.toNat = cld.size)
    (scaleNonzero : scale ≠ 0) (modulusNonzero : modulus ≠ 0)
    (hcount : count.toInt.toNat = lifted.size)
    (hnonempty : 0 < lifted.size) (hfits : lifted.size ≤ 2 ^ 31)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (answer : Generated.StrictVanHoeij.VanHoeijResult) :
    ∃ output, columnRound concreteVanHoeijControllerExecution f lifted cld
      (originalModulus : Int) modulus scale residual count target countSize
      scaleNonzero modulusNonzero answer = .ok output := by
  obtain ⟨initial, hi⟩ := makeInitialMatrix_succeeds count.toInt.toNat scale
  have hi' : Generated.StrictRecombine.makeInitialMatrix cld.size scale = .ok initial := by
    simpa only [countSize] using hi
  have hic := makeInitialMatrix_canonical cld (Generated.StrictRecombine.cldSpiralWidth cld)
    scale modulus initial scaleNonzero hi'
  obtain ⟨built, hb⟩ := buildCldMatrix_succeeds cld 0 target.toInt.toNat scale modulus initial hic
  have hbc := buildCldMatrix_canonical cld 0 target.toInt.toNat scale modulus initial built hic hb
  have hbn : 0 < built.1.size := by
    have hs := hbc.size_eq
    omega
  let ready := concreteVanHoeijControllerExecution.buildValid cld target.toInt.toNat scale modulus
    initial built scaleNonzero modulusNonzero hi' hb
  let bound := count.toInt * scale * scale + (Int32.ofNat built.2).toInt * residual * residual
  obtain ⟨reduced, hl⟩ := lllReduceWithNorms_succeeds built.1 ready hbn bound
  have hle : count.toInt.toNat ≤ built.1.size := by
    have hs := hbc.size_eq
    omega
  have hcolumns := lllReduceWithNorms_candidate_columns_valid built.1 ready bound
    count.toInt.toNat hle reduced hl
  obtain ⟨groups, hg⟩ := extractCandidates_succeeds reduced.2.1
    (Generated.StrictVanHoeij.shortSpanFromNorms reduced.2.2.2 bound) count.toInt.toNat hcolumns
  have hvalid : ∀ i (hi : i < groups.size),
      CandidateIndicesValid groups[i] (Array.replicate lifted.size false) := by
    by_cases he : (Generated.StrictVanHoeij.shortSpanFromNorms reduced.2.2.2 bound).isEmpty = true
    · have hz : groups = #[] := by
        simpa only [Generated.StrictRecombine.extractCandidates, he, ↓reduceIte, Except.ok.injEq] using hg.symm
      subst groups
      intro i hi
      simp at hi
    · have hn : (Generated.StrictVanHoeij.shortSpanFromNorms reduced.2.2.2 bound).isEmpty = false :=
        Bool.eq_false_iff.mpr he
      simpa only [hcount] using extractCandidates_indices_valid _ _ count.toInt.toNat groups
        (by omega) hn hg
  obtain ⟨result, hv⟩ := validatePartition_succeeds f lifted originalModulus hm groups hvalid
    hcanonical { answer with cldColumns := Int32.ofNat built.2, retainedRows := Int32.ofNat (Generated.StrictVanHoeij.shortSpanFromNorms reduced.2.2.2 bound).size }
  refine ⟨result, ?_⟩
  unfold columnRound
  split
  next fault hbad => rw [hi] at hbad; contradiction
  next initial' hi'' =>
    have he := Except.ok.inj (hi''.symm.trans hi)
    subst initial'
    split
    next fault hbad => rw [hb] at hbad; contradiction
    next built' hb' =>
      have he := Except.ok.inj (hb'.symm.trans hb)
      subst built'
      change (Generated.StrictVanHoeij.lllReduceWithNorms concreteLLLExecution built.1 ready bound >>= _) = _
      simp only [hl, bind, Except.bind]
      dsimp only [bound] at hg hv
      rw [hg]
      exact hv

theorem columnLoop_succeeds
    (f : SparsePolyZZ) (lifted cld : Array SparsePolyZZ)
    (originalModulus : Nat) (hm : 0 < originalModulus)
    (modulus scale residual : ZZ) (count width target : Int32)
    (countSize : count.toInt.toNat = cld.size)
    (scaleNonzero : scale ≠ 0) (modulusNonzero : modulus ≠ 0)
    (hcount : count.toInt.toNat = lifted.size)
    (hnonempty : 0 < lifted.size) (hfits : lifted.size ≤ 2 ^ 31)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (targetBound : target.toInt ≤ width.toInt)
    (answer : Generated.StrictVanHoeij.VanHoeijResult) :
    ∃ output, columnLoop concreteVanHoeijControllerExecution f lifted cld
      (originalModulus : Int) modulus scale residual count width target countSize
      scaleNonzero modulusNonzero targetBound answer = .ok output := by
  induction hmeasure : (width.toInt - target.toInt).toNat using Nat.strong_induction_on
      generalizing target answer with
  | h measure ih =>
    rw [columnLoop]
    split
    next ht =>
      obtain ⟨result, hr⟩ := columnRound_succeeds f lifted cld originalModulus hm modulus scale residual
        count target countSize scaleNonzero modulusNonzero hcount hnonempty hfits hcanonical answer
      rw [hr]
      dsimp only
      split
      · exact ⟨result.2, rfl⟩
      · split
        · exact ⟨result.2, rfl⟩
        next he =>
          apply ih _ ?_ _ _ result.2 rfl
          have hd := concreteVanHoeijControllerExecution.targetDecreases target width ht he targetBound
          omega
    next ht => exact ⟨answer, rfl⟩

end Refinement.StrictRecombine
