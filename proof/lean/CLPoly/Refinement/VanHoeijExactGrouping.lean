import CLPoly.Refinement.VanHoeijGroupingTotal
import CLPoly.Refinement.VanHoeijRetainedNonempty

set_option autoImplicit false

namespace Refinement.StrictRecombine

/-- The actual collector's membership test exactly recovers its class slot.
The machine-index bound is essential for injectivity of the emitted values. -/
theorem collectCandidateClasses_member_iff_class
    (classes : Array (Option Nat)) (classCount : Nat) (groups : Array (Array Int32))
    (hfits : classes.size ≤ 2 ^ 31)
    (hrun : Generated.StrictRecombine.collectCandidateClasses classes 0
      (Array.replicate classCount #[]) = .ok groups)
    (source classId : Nat) (hs : source < classes.size) :
    CandidateMember groups classId source.toUInt32.toInt32 ↔ classes[source]! = some classId := by
  constructor
  · intro hm
    obtain ⟨other, ho, hc, hv⟩ := collectCandidateClasses_from_empty_member_source
      classes classCount groups hrun classId source.toUInt32.toInt32 hm
    have he := congrArg (fun value : Int32 => value.toInt64.toNat) hv
    dsimp only at he
    rw [(nat_toUInt32_toInt32_nonnegative_and_toNat source (lt_of_lt_of_le hs hfits)).2,
      (nat_toUInt32_toInt32_nonnegative_and_toNat other (lt_of_lt_of_le ho hfits)).2] at he
    simpa only [← he] using hc
  · intro hc
    exact collectCandidateClasses_includes_source classes 0 (Array.replicate classCount #[])
      groups hrun source classId (Nat.zero_le _) hs hc

/-- A returned group containing a fixed source column consists of exactly
the columns accepted by its actual equality comparator. -/
theorem extractCandidates_class_of_column
    (transform : Generated.StrictRecombine.LLLMatrix) (rows : Array Nat)
    (count : Nat) (groups : Array (Array Int32)) (hfits : count ≤ 2 ^ 31)
    (hvalid : CandidateColumnsValid transform rows count) (hnonempty : rows.isEmpty = false)
    (hrun : Generated.StrictRecombine.extractCandidates rows transform count = .ok groups)
    (source : Nat) (hs : source < count) :
    ∃ classId, ∀ other, other < count →
      (CandidateMember groups classId other.toUInt32.toInt32 ↔
        Generated.StrictRecombine.candidateColumnsEqual transform rows source other 0 = .ok true) := by
  simp only [Generated.StrictRecombine.extractCandidates, hnonempty,
    Bool.false_eq_true, ↓reduceIte] at hrun
  cases hp : Generated.StrictRecombine.partitionCandidateColumns transform rows count
      0 0 (Array.replicate count none) with
  | error fault => rw [hp] at hrun; contradiction
  | ok result =>
    obtain ⟨classes, classCount⟩ := result
    rw [hp] at hrun
    obtain ⟨hsize, hassigned⟩ := partitionCandidateColumns_from_empty_all_assigned transform
      rows count classes classCount hp
    have hclosed := partitionCandidateColumns_preserves_closed transform rows count 0 0
      (Array.replicate count none) classes classCount hvalid (by simp)
      (by intro index hi; omega)
      (by intro i j hi hj _; simp [getElem!_pos, hi, hj]) hp
    have hpairwise := partitionCandidateColumns_from_empty_pairwise transform rows count
      classes classCount hvalid hp
    cases hc : classes[source]! with
    | none => exact (hassigned source hs hc).elim
    | some classId =>
      refine ⟨classId, ?_⟩
      intro other ho
      rw [collectCandidateClasses_member_iff_class classes classCount groups
        (by omega) hrun other classId (by omega)]
      constructor
      · intro hother
        exact hpairwise source other classId (by omega) (by omega) hc hother
      · intro heq
        exact (hclosed source other hs ho heq).symm.trans hc

/-- A preserved nonempty indicator whose support is constant on retained
rows is recovered as one exact physical output group. -/
theorem lllReduceWithNorms_exact_support_group
    (matrix : Generated.StrictRecombine.LLLMatrix) (hinput : ConcreteLLLInputValid matrix)
    (bound : ZZ) (output : Generated.StrictRecombine.LLLMatrix ×
      Generated.StrictRecombine.LLLMatrix × Array Nat × Array QQ)
    (hrun : Generated.StrictVanHoeij.lllReduceWithNorms concreteLLLExecution matrix hinput bound = .ok output)
    (count : Nat) (hcount : count ≤ matrix.size) (hfits : count ≤ 2 ^ 31)
    (indices : Array Nat) (hlegal : LegalCombination count indices.size indices)
    (hnonempty : indices.toList ≠ [])
    (coefficients : Fin matrix.size → ZZ)
    (hcoordinates : ∀ index, index.val < count →
      coefficients index = if index.val ∈ indices.toList then 1 else 0)
    (hshort : dotProduct
      (fun column => ((Matrix.vecMul coefficients
        (basisPrefixMatrix matrix matrix.size matrix.size) column : ZZ) : QQ))
      (fun column => ((Matrix.vecMul coefficients
        (basisPrefixMatrix matrix matrix.size matrix.size) column : ZZ) : QQ)) ≤ (bound : QQ))
    (hconstant : ∀ row : Fin matrix.size,
      row.val < Generated.StrictVanHoeij.retainedPrefixCount output.2.2.2 (bound : QQ) output.2.2.2.size →
      ∀ left right : Fin count, left.val ∈ indices.toList → right.val ∈ indices.toList →
        (output.2.1[row.val]!)[left.val]! = (output.2.1[row.val]!)[right.val]!) :
    ∃ groups classId,
      Generated.StrictRecombine.extractCandidates
        (Generated.StrictVanHoeij.shortSpanFromNorms output.2.2.2 bound) output.2.1 count = .ok groups ∧
      classId < groups.size ∧
      ∀ index, index < count →
        (CandidateMember groups classId index.toUInt32.toInt32 ↔ index ∈ indices.toList) := by
  obtain ⟨first, hf⟩ := List.exists_mem_of_ne_nil indices.toList hnonempty
  obtain ⟨position, hp, hv⟩ := List.mem_iff_getElem.mp hf
  have hp' : position < indices.size := by simpa using hp
  have hfirst : first < count := by
    have h := hlegal.2.2 position hp'
    simpa only [getElem!_pos indices position hp', ← Array.getElem_toList hp', hv] using h
  let firstFin : Fin matrix.size := ⟨first, lt_of_lt_of_le hfirst hcount⟩
  have hne : coefficients ≠ 0 := by
    intro hz
    have h := hcoordinates firstFin hfirst
    simp [firstFin, hz, hf] at h
  have hrows := lllReduceWithNorms_retained_nonempty_of_short matrix hinput bound output hrun
    coefficients hne hshort
  have hvalid := lllReduceWithNorms_candidate_columns_valid matrix hinput bound count hcount output hrun
  obtain ⟨groups, hg⟩ := extractCandidates_succeeds output.2.1
    (Generated.StrictVanHoeij.shortSpanFromNorms output.2.2.2 bound) count hvalid
  obtain ⟨classId, hclass⟩ := extractCandidates_class_of_column output.2.1
    (Generated.StrictVanHoeij.shortSpanFromNorms output.2.2.2 bound) count groups hfits
    hvalid hrows hg first hfirst
  have hmember := (hclass first hfirst).mpr
    (candidateColumnsEqual_refl _ _ count first hvalid hfirst)
  refine ⟨groups, classId, hg, hmember.1, ?_⟩
  intro index hi
  rw [hclass index hi]
  constructor
  · intro hcompare
    have hequal : coefficients firstFin = coefficients ⟨index, lt_of_lt_of_le hi hcount⟩ := by
      apply lllReduceWithNorms_equal_original_coordinates matrix hinput bound output hrun
        coefficients hshort firstFin ⟨index, lt_of_lt_of_le hi hcount⟩
      intro row hrow
      obtain ⟨p, hp, he⟩ := Array.mem_iff_getElem.mp hrow
      have hs := candidateColumnsEqual_true_sound output.2.1
        (Generated.StrictVanHoeij.shortSpanFromNorms output.2.2.2 bound) first index 0 hcompare
        p (Nat.zero_le _) hp
      have hr : (Generated.StrictVanHoeij.shortSpanFromNorms output.2.2.2 bound)[p]! = row := by
        simpa only [getElem!_pos
          (Generated.StrictVanHoeij.shortSpanFromNorms output.2.2.2 bound) p hp] using he
      dsimp only at hs
      rw [hr] at hs
      exact hs.2.2.2
    rw [hcoordinates firstFin hfirst, hcoordinates _ hi] at hequal
    dsimp only [firstFin] at hequal
    by_contra hn
    simp [hf, hn] at hequal
  · intro hm
    apply lllReduceWithNorms_columns_equal_of_retained matrix hinput bound output hrun
      firstFin ⟨index, lt_of_lt_of_le hi hcount⟩
    intro row hkeep
    exact hconstant row hkeep ⟨first, hfirst⟩ ⟨index, hi⟩ hf hm

end Refinement.StrictRecombine
