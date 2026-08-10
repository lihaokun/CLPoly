/-
Direct structural refinement lemmas for the generated C++
`__upoly_subtract_one` range loop.
-/
import CLPoly.Generated.StrictEDF
import CLPoly.Refinement.DDFSubtractX

set_option autoImplicit false
set_option maxHeartbeats 0

namespace Refinement.StrictEDF

private lemma drop_eq_getElem_cons (input : SparsePolyZp) (i : Nat)
    (hi : i < input.size) :
    input.toList.drop i = input[i] :: input.toList.drop (i + 1) := by
  have hlen : i < input.toList.length := by simpa using hi
  have hdrop := drop_eq_get_cons_general input.toList i hlen
  simpa using hdrop

/-- Once the constant term has been processed, strict ordering means the
remaining suffix has no second constant term.  The generated loop therefore
copies that suffix verbatim. -/
theorem subtractOneLoop_after_found
    (p : UInt64) : ∀ (i : Nat) (result input : SparsePolyZp),
    (∀ term ∈ input.toList.drop i,
      term.1.deg.toUInt64.toInt64 ≠ (0 : Int64)) →
    let out := Generated.StrictEDF._loop___upoly_subtract_one_0_raw_ir
      i true result input p
    out.1 = true ∧
      out.2.toList = result.toList ++ input.toList.drop i := by
  intro i result input hno
  induction hn : input.size - i using Nat.strongRecOn generalizing i result with
  | ind n ih =>
      rw [Generated.StrictEDF._loop___upoly_subtract_one_0_raw_ir.eq_1]
      by_cases hi : i < input.size
      · simp only [hi, ↓reduceDIte]
        let term := input[i]
        have hdrop : input.toList.drop i =
            term :: input.toList.drop (i + 1) := by
          simpa [term] using drop_eq_getElem_cons input i hi
        have hterm : term.1.deg.toUInt64.toInt64 ≠ (0 : Int64) :=
          hno term (by rw [hdrop]; simp)
        have hdegree : term.1.deg ≠ 0 := by
          intro hzero
          rw [hzero] at hterm
          simp at hterm
        have hguard :
            (input[i].1.deg == Int64.toNatClampNeg (0 : Int64)) = false := by
          simp [term, hdegree]
        have htail : ∀ x ∈ input.toList.drop (i + 1),
            x.1.deg.toUInt64.toInt64 ≠ (0 : Int64) := by
          intro x hx
          exact hno x (by rw [hdrop]; exact List.mem_cons_of_mem _ hx)
        have hdec : input.size - (i + 1) < n := by omega
        have hrec := ih (input.size - (i + 1)) hdec (i + 1)
          (result.push term) htail rfl
        simp only [hguard, Bool.false_eq_true, ↓reduceIte]
        refine ⟨hrec.1, ?_⟩
        rw [hrec.2, Array.toList_push, hdrop]
        simp [List.append_assoc]
      · have hdrop : input.toList.drop i = [] := by
          apply List.drop_eq_nil_of_le
          simpa using Nat.le_of_not_gt hi
        simp [hi, hdrop]

end Refinement.StrictEDF
