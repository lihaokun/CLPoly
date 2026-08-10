/-
Direct structural refinement lemmas for the generated C++
`__upoly_subtract_one` range loop.
-/
import CLPoly.Generated.StrictEDF
import CLPoly.Refinement.DDFSubtractX

set_option autoImplicit false
set_option maxHeartbeats 0

namespace Refinement.StrictEDF

open Polynomial
open CLPoly.Math

private lemma drop_eq_getElem_cons (input : SparsePolyZp) (i : Nat)
    (hi : i < input.size) :
    input.toList.drop i = input[i] :: input.toList.drop (i + 1) := by
  have hlen : i < input.toList.length := by simpa using hi
  have hdrop := drop_eq_get_cons_general input.toList i hlen
  simpa using hdrop

private lemma pairwise_drop_tail (input : SparsePolyZp) (i : Nat)
    (hi : i < input.size)
    (hs : (input.toList.drop i).Pairwise
      (fun a b => b.1.deg < a.1.deg)) :
    (input.toList.drop (i + 1)).Pairwise
      (fun a b => b.1.deg < a.1.deg) := by
  rw [drop_eq_getElem_cons input i hi, List.pairwise_cons] at hs
  exact hs.2

private lemma no_constant_in_tail_of_constant
    (input : SparsePolyZp) (i : Nat) (hi : i < input.size)
    (hs : (input.toList.drop i).Pairwise
      (fun a b => b.1.deg < a.1.deg))
    (hzero : input[i].1.deg = 0) :
    ∀ x ∈ input.toList.drop (i + 1),
      x.1.deg.toUInt64.toInt64 ≠ (0 : Int64) := by
  rw [drop_eq_getElem_cons input i hi, List.pairwise_cons] at hs
  intro x hx
  have hlt := hs.1 x hx
  rw [hzero] at hlt
  omega

/-- List trace of the generated branch order.  This is not an L2 polynomial
operation: it mirrors the concrete coefficient tests and array actions. -/
private def subtractOneTerms (p : UInt64) :
    List (UMonomial × Zp) → List (UMonomial × Zp)
  | [] => [(UMonomial.mk 0, Zp.ofInt (p - 1).toInt p)]
  | term :: rest =>
      if term.1.deg = 0 then
        let coefficient := term.2 - Generated.StrictEDF.__make_zp_ir 1 p
        if coefficient.val = 0 then rest
        else (UMonomial.mk 0, coefficient) :: rest
      else
        term :: subtractOneTerms p rest

/-- Polynomial meaning of the structural trace.  The induction follows the
same coefficient tests as the generated loop. -/
theorem subtractOneTerms_toPoly (q : UInt64) [Fact (Nat.Prime q.toNat)]
    (h2q : 2 * q.toNat ≤ UInt64.size) :
    ∀ terms : List (UMonomial × Zp),
      SparsePolyZp.AllReduced q.toNat terms →
      listSum q.toNat (subtractOneTerms q terms) =
        listSum q.toNat terms - 1 := by
  intro terms
  induction terms with
  | nil =>
      intro _
      simp [subtractOneTerms, listSum, strict_minus_one_toZMod q rfl,
        sub_eq_add_neg]
  | cons term rest ih =>
      intro hreduced
      have hterm := hreduced term List.mem_cons_self
      have hrest : SparsePolyZp.AllReduced q.toNat rest :=
        fun x hx => hreduced x (List.mem_cons_of_mem _ hx)
      by_cases hzero : term.1.deg = 0
      · let coefficient :=
          term.2 - Generated.StrictEDF.__make_zp_ir 1 q
        have hcoefficient : Zp.toZMod q.toNat coefficient =
            Zp.toZMod q.toNat term.2 - 1 := by
          simpa [coefficient, Generated.StrictEDF.__make_zp_ir,
            Generated.StrictDDF.__make_zp_ir] using
            strict_sub_one_toZMod h2q term.2 q hterm.1 hterm.2 rfl
        by_cases hc : coefficient.val = 0
        · have hczero : Zp.toZMod q.toNat coefficient = 0 := by
            simp [Zp.toZMod, hc]
          have htermOne : Zp.toZMod q.toNat term.2 - 1 = 0 :=
            hcoefficient.symm.trans hczero
          rw [subtractOneTerms, if_pos hzero, if_pos hc]
          simp only [listSum]
          rw [hzero]
          simp only [Polynomial.monomial_zero_left]
          have ht : Zp.toZMod q.toNat term.2 = 1 := sub_eq_zero.mp htermOne
          rw [ht]
          simp
        · rw [subtractOneTerms, if_pos hzero, if_neg hc]
          simp only [listSum]
          rw [hzero, hcoefficient]
          simp only [Polynomial.monomial_zero_left]
          rw [map_sub, map_one]
          ring
      · rw [subtractOneTerms, if_neg hzero, listSum, ih hrest, listSum]
        ring

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

/-- If the concrete input suffix contains no constant term, the generated loop
never takes its update branch: it preserves `found = false` and copies the
suffix verbatim.  This is the exact execution fact used by the source
epilogue which appends the representation of `-1`. -/
theorem subtractOneLoop_no_constant
    (p : UInt64) : ∀ (i : Nat) (result input : SparsePolyZp),
    (∀ term ∈ input.toList.drop i,
      term.1.deg.toUInt64.toInt64 ≠ (0 : Int64)) →
    let out := Generated.StrictEDF._loop___upoly_subtract_one_0_raw_ir
      i false result input p
    out.1 = false ∧
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

/-- Starting from `found = false`, the exact generated loop and its source
epilogue produce the structural trace above for every strictly ordered suffix. -/
theorem subtractOneLoop_trace
    (p : UInt64) : ∀ (i : Nat) (result input : SparsePolyZp),
    (input.toList.drop i).Pairwise
      (fun a b => b.1.deg < a.1.deg) →
    let out := Generated.StrictEDF._loop___upoly_subtract_one_0_raw_ir
      i false result input p
    (if !out.1 then
      (out.2.push (UMonomial.mk 0, Zp.ofInt (p - 1).toInt p)).toList
    else out.2.toList) =
      result.toList ++ subtractOneTerms p (input.toList.drop i) := by
  intro i result input hs
  induction hn : input.size - i using Nat.strongRecOn generalizing i result with
  | ind n ih =>
      rw [Generated.StrictEDF._loop___upoly_subtract_one_0_raw_ir.eq_1]
      by_cases hi : i < input.size
      · simp only [hi, ↓reduceDIte]
        let term := input[i]
        have hdrop : input.toList.drop i =
            term :: input.toList.drop (i + 1) := by
          simpa [term] using drop_eq_getElem_cons input i hi
        have hsTail := pairwise_drop_tail input i hi hs
        have hdec : input.size - (i + 1) < n := by omega
        by_cases hzero : term.1.deg = 0
        · have hguard :
              (input[i].1.deg == Int64.toNatClampNeg (0 : Int64)) = true := by
            simp [term, hzero]
          let coefficient := term.2 - Generated.StrictEDF.__make_zp_ir 1 p
          let nextResult := if coefficient.val = 0 then result
            else result.push (UMonomial.mk 0, coefficient)
          have hnoTail := no_constant_in_tail_of_constant input i hi hs (by
            simpa [term] using hzero)
          have hpost := subtractOneLoop_after_found p (i + 1)
            nextResult input hnoTail
          simp only [hguard, ↓reduceIte]
          by_cases hcoefficient : coefficient.val = 0
          · have hbool :
                ((input[i].2 - Generated.StrictEDF.__make_zp_ir 1 p).val != 0) =
                  false := by
              simp [coefficient, term, hcoefficient]
            simp only [hbool, Bool.false_eq_true, ↓reduceIte]
            have hnext : nextResult = result := by
              simp [nextResult, hcoefficient]
            rw [hnext] at hpost
            rw [hpost.1]
            simp only [Bool.not_true, ↓reduceIte]
            rw [hpost.2, hdrop, subtractOneTerms, if_pos hzero]
            simp [coefficient, hcoefficient]
          · have hbool :
                ((input[i].2 - Generated.StrictEDF.__make_zp_ir 1 p).val != 0) =
                  true := by
              simp [coefficient, term, hcoefficient]
            simp only [hbool, ↓reduceIte]
            have hnext : nextResult =
                result.push (UMonomial.mk 0, coefficient) := by
              simp [nextResult, hcoefficient]
            rw [hnext] at hpost
            rw [hpost.1]
            simp only [Bool.not_true, ↓reduceIte]
            rw [hpost.2, hdrop, subtractOneTerms, if_pos hzero]
            simp [coefficient, hcoefficient, Array.toList_push,
              List.append_assoc]
        · have hguard :
              (input[i].1.deg == Int64.toNatClampNeg (0 : Int64)) = false := by
            simp [term, hzero]
          have hrec := ih (input.size - (i + 1)) hdec (i + 1)
            (result.push term) hsTail rfl
          simp only [hguard, Bool.false_eq_true, ↓reduceIte]
          rw [hrec, hdrop, subtractOneTerms, if_neg hzero,
            Array.toList_push]
          simp [List.append_assoc]
      · have hdrop : input.toList.drop i = [] := by
          apply List.drop_eq_nil_of_le
          simpa using Nat.le_of_not_gt hi
        simp [hi, hdrop, subtractOneTerms, Array.toList_push]

/-- Direct raw-entry equation for the C++ branch where the input has no
constant term.  The result is obtained from the actual generated loop and its
actual normalization epilogue. -/
theorem subtractOneRaw_no_constant (h : SparsePolyZp) (p : UInt64)
    (hno : ∀ term ∈ h.toList,
      term.1.deg.toUInt64.toInt64 ≠ (0 : Int64)) :
    Generated.StrictEDF.__upoly_subtract_one_raw_ir h p =
      .ok (SparsePolyZp.normalization
        (h.push (UMonomial.mk 0, Zp.ofInt (p - 1).toInt p))) := by
  let out := Generated.StrictEDF._loop___upoly_subtract_one_0_raw_ir
    0 false #[] h p
  have hloop := subtractOneLoop_no_constant p 0 #[] h (by simpa using hno)
  change out.1 = false ∧ out.2.toList = h.toList at hloop
  have hout : out.2 = h := Array.toList_inj.mp hloop.2
  simp [Generated.StrictEDF.__upoly_subtract_one_raw_ir, out,
    hloop.1, hout]

/-- Semantic refinement of the same concrete no-constant execution branch. -/
theorem subtractOneRaw_no_constant_refines
    (h : SparsePolyZp) (p : UInt64) [Fact (Nat.Prime p.toNat)]
    (hno : ∀ term ∈ h.toList,
      term.1.deg.toUInt64.toInt64 ≠ (0 : Int64)) :
    ∃ output,
      Generated.StrictEDF.__upoly_subtract_one_raw_ir h p = .ok output ∧
      SparsePolyZp.toPoly p.toNat output =
        SparsePolyZp.toPoly p.toNat h - 1 := by
  let output := SparsePolyZp.normalization
    (h.push (UMonomial.mk 0, Zp.ofInt (p - 1).toInt p))
  refine ⟨output, subtractOneRaw_no_constant h p hno, ?_⟩
  simp only [output, normalization_toPoly, toPoly_push]
  rw [strict_minus_one_toZMod p rfl]
  simp [sub_eq_add_neg]

/-- The exact generated C++ entry always terminates on a canonical sparse
input and decodes to subtraction by one.  Both source epilogue branches are
covered: normalization is used only when the generated loop found no constant
term, exactly as in `__upoly_subtract_one_raw_ir`. -/
theorem __upoly_subtract_one_raw_ir_refines
    (h : SparsePolyZp) (p : UInt64) [Fact (Nat.Prime p.toNat)]
    (h2p : 2 * p.toNat ≤ UInt64.size)
    (hcanonical : SparsePolyZp.Canonical p.toNat h) :
    ∃ output,
      Generated.StrictEDF.__upoly_subtract_one_raw_ir h p = .ok output ∧
      SparsePolyZp.toPoly p.toNat output =
        SparsePolyZp.toPoly p.toNat h - 1 := by
  let out := Generated.StrictEDF._loop___upoly_subtract_one_0_raw_ir
    0 false #[] h p
  have htrace := subtractOneLoop_trace p 0 #[] h
    (List.isChain_iff_pairwise.mp hcanonical.2.1)
  change (if !out.1 then
      (out.2.push (UMonomial.mk 0,
        Zp.ofInt (p - 1).toInt p)).toList
    else out.2.toList) = subtractOneTerms p h.toList at htrace
  by_cases hfound : out.1 = true
  · refine ⟨out.2, ?_, ?_⟩
    · simp [Generated.StrictEDF.__upoly_subtract_one_raw_ir, out,
        hfound]
    · simp [hfound] at htrace
      change listSum p.toNat out.2.toList =
        listSum p.toNat h.toList - 1
      rw [htrace]
      exact subtractOneTerms_toPoly p h2p h.toList hcanonical.1
  · have hnotFound : out.1 = false := Bool.eq_false_of_not_eq_true hfound
    let output := SparsePolyZp.normalization
      (out.2.push (UMonomial.mk 0, Zp.ofInt (p - 1).toInt p))
    refine ⟨output, ?_, ?_⟩
    · simp [Generated.StrictEDF.__upoly_subtract_one_raw_ir, out,
        hnotFound, output]
    · rw [show SparsePolyZp.toPoly p.toNat output =
          SparsePolyZp.toPoly p.toNat
            (out.2.push (UMonomial.mk 0,
              Zp.ofInt (p - 1).toInt p)) by
          exact normalization_toPoly _]
      simp [hnotFound] at htrace
      change listSum p.toNat
          (out.2.push (UMonomial.mk 0,
            Zp.ofInt (p - 1).toInt p)).toList =
        listSum p.toNat h.toList - 1
      rw [Array.toList_push]
      rw [htrace]
      exact subtractOneTerms_toPoly p h2p h.toList hcanonical.1

end Refinement.StrictEDF
