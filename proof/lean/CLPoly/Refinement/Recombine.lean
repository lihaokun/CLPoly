/-
  No strict recombination refinement theorem is exported yet.

  The former theorem equated the generated C++ result with an arbitrary UFD
  existence witness and contained a placeholder proof.  Both the witness-based
  implementation and the claimed refinement have been removed.  A replacement
  must prove the generated candidate/extraction loops directly.
-/
import CLPoly.Algorithm.Recombine
import CLPoly.Generated.StrictRecombine
import CLPoly.Refinement.Basic
import CLPoly.Refinement.Hensel
import Batteries.Data.Array.Lemmas
import Mathlib.Analysis.Matrix.PosDef
import Mathlib.Data.Rat.Star
import Mathlib.Data.Real.StarOrdered
import Mathlib.LinearAlgebra.LinearIndependent.BaseChange
import Mathlib.LinearAlgebra.Matrix.Block

set_option autoImplicit false

open Polynomial
open CLPoly.Math

namespace Refinement.StrictRecombine

/-- Project an integer-polynomial congruence from a larger modulus to any
divisor modulus.  This is the transport used to compare the actual `p^k`
Hensel execution with the selected-prime factorization at `p`. -/
theorem polynomialMap_eq_of_modulus_dvd (small large : Nat)
    (hdivides : small ∣ large) (left right : Polynomial Int)
    (heq : Polynomial.map (Int.castRingHom (ZMod large)) left =
      Polynomial.map (Int.castRingHom (ZMod large)) right) :
    Polynomial.map (Int.castRingHom (ZMod small)) left =
      Polynomial.map (Int.castRingHom (ZMod small)) right := by
  have projected : Polynomial.map (ZMod.castHom hdivides (ZMod small))
        (Polynomial.map (Int.castRingHom (ZMod large)) left) =
      Polynomial.map (ZMod.castHom hdivides (ZMod small))
        (Polynomial.map (Int.castRingHom (ZMod large)) right) :=
    congrArg _ heq
  simp only [Polynomial.map_map] at projected
  convert projected using 1 <;> congr 1 <;> ext value <;> simp

private def appendZeroSuffix
    (matrix : Generated.StrictRecombine.LLLMatrix) (index : Nat) :
    Generated.StrictRecombine.LLLMatrix :=
  ((matrix.toList.drop index).map fun row => row.push 0).toArray

theorem zeroMatrixRowLoop_size (columns index : Nat) (row : Array ZZ)
    (hindex : index ≤ columns) :
    (Generated.StrictRecombine.zeroMatrixRowLoop columns index row).size =
      row.size + (columns - index) := by
  induction hremaining : columns - index using Nat.strong_induction_on
      generalizing index row with
  | h remaining ih =>
      rw [Generated.StrictRecombine.zeroMatrixRowLoop]
      split
      next hmore =>
        rw [ih (columns - (index + 1)) (by omega) (index + 1)
          (row.push 0) (by omega) rfl]
        simp only [Array.size_push]
        omega
      next hdone => omega

theorem zeroMatrixRow_size (columns : Nat) :
    (Generated.StrictRecombine.zeroMatrixRow columns).size = columns := by
  unfold Generated.StrictRecombine.zeroMatrixRow
  simpa using zeroMatrixRowLoop_size columns 0 (#[] : Array ZZ) (by omega)

theorem zeroMatrixRowLoop_get_of_lt (columns index : Nat) (row : Array ZZ)
    (position : Nat) (hposition : position < row.size) :
    (Generated.StrictRecombine.zeroMatrixRowLoop columns index row)[position]! =
      row[position] := by
  induction hremaining : columns - index using Nat.strong_induction_on
      generalizing index row with
  | h remaining ih =>
      rw [Generated.StrictRecombine.zeroMatrixRowLoop]
      split
      next hmore =>
        have hdecrease : columns - (index + 1) < remaining := by omega
        rw [ih (columns - (index + 1)) hdecrease (index + 1)
          (row.push 0) (by simp only [Array.size_push]; omega) rfl]
        simpa only [Array.getElem_push_lt hposition]
      next hdone =>
        simpa only [getElem!_pos row position hposition]

theorem zeroMatrixRowLoop_suffix_get (columns index : Nat) (row : Array ZZ)
    (hindex : index ≤ columns) (offset : Nat)
    (hoffset : offset < columns - index) :
    (Generated.StrictRecombine.zeroMatrixRowLoop columns index row)[row.size + offset]! =
      0 := by
  induction offset generalizing index row with
  | zero =>
      rw [Generated.StrictRecombine.zeroMatrixRowLoop]
      split
      next hmore =>
        simp only [Nat.add_zero]
        rw [zeroMatrixRowLoop_get_of_lt columns (index + 1) (row.push 0)
          row.size (by simp)]
        simp
      next hdone => omega
  | succ offset ih =>
      rw [Generated.StrictRecombine.zeroMatrixRowLoop]
      split
      next hmore =>
        have htail := ih (index + 1) (row.push 0) (by omega) (by omega)
        simpa only [Array.size_push, Nat.add_assoc, Nat.add_left_comm,
          Nat.add_comm] using htail
      next hdone => omega

theorem zeroMatrixRow_get (columns position : Nat)
    (hposition : position < columns) :
    (Generated.StrictRecombine.zeroMatrixRow columns)[position]! = 0 := by
  unfold Generated.StrictRecombine.zeroMatrixRow
  simpa using zeroMatrixRowLoop_suffix_get columns 0 (#[] : Array ZZ)
    (by omega) position (by omega)

theorem zeroMatrixLoop_size (rows columns index : Nat)
    (matrix : Generated.StrictRecombine.LLLMatrix) (hindex : index ≤ rows) :
    (Generated.StrictRecombine.zeroMatrixLoop rows columns index matrix).size =
      matrix.size + (rows - index) := by
  induction hremaining : rows - index using Nat.strong_induction_on
      generalizing index matrix with
  | h remaining ih =>
      rw [Generated.StrictRecombine.zeroMatrixLoop]
      split
      next hmore =>
        rw [ih (rows - (index + 1)) (by omega) (index + 1)
          (matrix.push (Generated.StrictRecombine.zeroMatrixRow columns))
          (by omega) rfl]
        simp only [Array.size_push]
        omega
      next hdone => omega

theorem zeroMatrix_size (rows columns : Nat) :
    (Generated.StrictRecombine.zeroMatrix rows columns).size = rows := by
  unfold Generated.StrictRecombine.zeroMatrix
  simpa using zeroMatrixLoop_size rows columns 0
    (#[] : Generated.StrictRecombine.LLLMatrix) (by omega)

theorem zeroMatrixLoop_get_of_lt (rows columns index : Nat)
    (matrix : Generated.StrictRecombine.LLLMatrix) (position : Nat)
    (hposition : position < matrix.size) :
    (Generated.StrictRecombine.zeroMatrixLoop rows columns index matrix)[position]! =
      matrix[position] := by
  induction hremaining : rows - index using Nat.strong_induction_on
      generalizing index matrix with
  | h remaining ih =>
      rw [Generated.StrictRecombine.zeroMatrixLoop]
      split
      next hmore =>
        have hdecrease : rows - (index + 1) < remaining := by omega
        rw [ih (rows - (index + 1)) hdecrease (index + 1)
          (matrix.push (Generated.StrictRecombine.zeroMatrixRow columns))
          (by simp only [Array.size_push]; omega) rfl]
        simpa only [Array.getElem_push_lt hposition]
      next hdone =>
        simpa only [getElem!_pos matrix position hposition]

theorem zeroMatrixLoop_suffix_get (rows columns index : Nat)
    (matrix : Generated.StrictRecombine.LLLMatrix) (hindex : index ≤ rows)
    (offset : Nat) (hoffset : offset < rows - index) :
    (Generated.StrictRecombine.zeroMatrixLoop rows columns index matrix)[matrix.size + offset]! =
      Generated.StrictRecombine.zeroMatrixRow columns := by
  induction offset generalizing index matrix with
  | zero =>
      rw [Generated.StrictRecombine.zeroMatrixLoop]
      split
      next hmore =>
        simp only [Nat.add_zero]
        rw [zeroMatrixLoop_get_of_lt rows columns (index + 1)
          (matrix.push (Generated.StrictRecombine.zeroMatrixRow columns))
          matrix.size (by simp)]
        simp
      next hdone => omega
  | succ offset ih =>
      rw [Generated.StrictRecombine.zeroMatrixLoop]
      split
      next hmore =>
        have htail := ih (index + 1)
          (matrix.push (Generated.StrictRecombine.zeroMatrixRow columns))
          (by omega) (by omega)
        simpa only [Array.size_push, Nat.add_assoc, Nat.add_left_comm,
          Nat.add_comm] using htail
      next hdone => omega

theorem zeroMatrix_get (rows columns row : Nat) (hrow : row < rows) :
    (Generated.StrictRecombine.zeroMatrix rows columns)[row]! =
      Generated.StrictRecombine.zeroMatrixRow columns := by
  unfold Generated.StrictRecombine.zeroMatrix
  simpa using zeroMatrixLoop_suffix_get rows columns 0
    (#[] : Generated.StrictRecombine.LLLMatrix) (by omega) row (by omega)

theorem zeroMatrix_entry (rows columns row column : Nat)
    (hrow : row < rows) (hcolumn : column < columns) :
    ((Generated.StrictRecombine.zeroMatrix rows columns)[row]!)[column]! = 0 := by
  rw [zeroMatrix_get rows columns row hrow,
    zeroMatrixRow_get columns column hcolumn]

theorem appendZeroColumnLoop_eq
    (matrix : Generated.StrictRecombine.LLLMatrix) (index : Nat)
    (result : Generated.StrictRecombine.LLLMatrix) :
    Generated.StrictRecombine.appendZeroColumnLoop matrix index result =
      .ok (result ++ appendZeroSuffix matrix index) := by
  induction hmeasure : matrix.size - index using Nat.strong_induction_on
      generalizing index result with
  | h measure ih =>
      rw [Generated.StrictRecombine.appendZeroColumnLoop]
      split
      next hindex =>
        rw [ih (matrix.size - (index + 1)) (by omega)
          (index + 1) (result.push (matrix[index].push 0)) rfl]
        have hsuffix : matrix.toList.drop index = matrix[index] ::
            matrix.toList.drop (index + 1) := by
          simpa using List.drop_eq_getElem_cons
            (l := matrix.toList) (i := index) (by simpa using hindex)
        simp [appendZeroSuffix, hsuffix, Array.push]
      next hindex =>
        have hle : matrix.size ≤ index := Nat.le_of_not_gt hindex
        simp [appendZeroSuffix, List.drop_eq_nil_iff.mpr hle]

theorem appendZeroColumn_eq
    (matrix : Generated.StrictRecombine.LLLMatrix) :
    Generated.StrictRecombine.appendZeroColumn matrix =
      .ok (appendZeroSuffix matrix 0) := by
  simpa [Generated.StrictRecombine.appendZeroColumn] using
    appendZeroColumnLoop_eq matrix 0 #[]

theorem appendZeroSuffix_size
    (matrix : Generated.StrictRecombine.LLLMatrix) :
    (appendZeroSuffix matrix 0).size = matrix.size := by
  simp [appendZeroSuffix]

theorem appendZeroSuffix_row
    (matrix : Generated.StrictRecombine.LLLMatrix) (row : Nat)
    (hrow : row < matrix.size) :
    (appendZeroSuffix matrix 0)[row]! = matrix[row]!.push 0 := by
  simp [appendZeroSuffix, hrow]

theorem arrayPush_getElem!_lt {α : Type*} [Inhabited α]
    (input : Array α) (value : α) (index : Nat) (hindex : index < input.size) :
    (input.push value)[index]! = input[index]! := by
  rw [getElem!_pos _ index (by simp only [Array.size_push]; omega),
    Array.getElem_push_lt hindex, ← getElem!_pos _ index hindex]

theorem arrayPush_getElem!_last {α : Type*} [Inhabited α]
    (input : Array α) (value : α) :
    (input.push value)[input.size]! = value := by
  rw [getElem!_pos _ input.size (by simp), Array.getElem_push_eq]

theorem fillCldDataRowLoop_size
    (cld : Array SparsePolyZZ) (degree index : Nat) (row output : Array ZZ)
    (hrun : Generated.StrictRecombine.fillCldDataRowLoop cld degree index row =
      .ok output) :
    output.size = row.size := by
  induction hmeasure : cld.size - index using Nat.strong_induction_on
      generalizing index row with
  | h measure ih =>
      rw [Generated.StrictRecombine.fillCldDataRowLoop] at hrun
      split at hrun
      next hindex =>
        split at hrun
        next hrow =>
          have htail := ih (cld.size - (index + 1)) (by omega)
            (index + 1)
            (row.set index
              (Generated.StrictRecombine.sparseCoeff cld[index] degree))
            hrun rfl
          simpa using htail
        next hrow => contradiction
      next hindex =>
        exact congrArg Array.size (Except.ok.inj hrun).symm

theorem appendCldColumn_shape
    (matrix output : Generated.StrictRecombine.LLLMatrix)
    (cld : Array SparsePolyZZ) (existingColumns spiralDegree : Nat)
    (hdimension : matrix.size = cld.size + existingColumns)
    (hrun : Generated.StrictRecombine.appendCldColumn matrix cld
      existingColumns spiralDegree = .ok output) :
    ∃ finalRow : Array ZZ,
      finalRow.size = matrix.size + 1 ∧
      finalRow[matrix.size]! = 1 ∧
      output = (appendZeroSuffix matrix 0).push finalRow := by
  unfold Generated.StrictRecombine.appendCldColumn at hrun
  rw [appendZeroColumn_eq] at hrun
  simp only [Except.bind] at hrun
  cases hfill : Generated.StrictRecombine.fillCldDataRowLoop cld spiralDegree
      0 (Generated.StrictRecombine.zeroMatrixRow
        (cld.size + existingColumns + 1)) with
  | error fault =>
      rw [hfill] at hrun
      change (Except.error fault : RawExec
        Generated.StrictRecombine.LLLMatrix) = .ok output at hrun
      contradiction
  | ok row =>
      rw [hfill] at hrun
      change (if hidentity : cld.size + existingColumns < row.size then
          (Except.ok ((appendZeroSuffix matrix 0).push
            (row.set (cld.size + existingColumns) 1 hidentity)) :
              RawExec Generated.StrictRecombine.LLLMatrix)
        else Except.error (RawFault.outOfBounds
          (cld.size + existingColumns) row.size)) =
          Except.ok output at hrun
      split at hrun
      next hidentity =>
        have hout := Except.ok.inj hrun
        have hrowSize := fillCldDataRowLoop_size cld spiralDegree 0
          (Generated.StrictRecombine.zeroMatrixRow
            (cld.size + existingColumns + 1)) row hfill
        rw [zeroMatrixRow_size, ← hdimension] at hrowSize
        let finalRow := row.set (cld.size + existingColumns) 1 hidentity
        refine ⟨finalRow, by simpa [finalRow] using hrowSize,
          ?_, hout.symm⟩
        simp [finalRow, hdimension, hrowSize]
      next hidentity => contradiction

theorem zeroQQRowLoop_size (columns index : Nat) (row : Array QQ)
    (hindex : index ≤ columns) :
    (Generated.StrictRecombine.zeroQQRowLoop columns index row).size =
      row.size + (columns - index) := by
  induction hremaining : columns - index using Nat.strong_induction_on
      generalizing index row with
  | h remaining ih =>
      rw [Generated.StrictRecombine.zeroQQRowLoop]
      split
      next hmore =>
        have hdecrease : columns - (index + 1) < remaining := by omega
        rw [ih (columns - (index + 1)) hdecrease (index + 1) (row.push 0)
          (by omega) rfl]
        simp only [Array.size_push]
        omega
      next hdone => omega

theorem zeroQQRow_size (columns : Nat) :
    (Generated.StrictRecombine.zeroQQRow columns).size = columns := by
  unfold Generated.StrictRecombine.zeroQQRow
  simpa using zeroQQRowLoop_size columns 0 (#[] : Array QQ) (by omega)

theorem zeroQQRowLoop_get_of_lt (columns index : Nat) (row : Array QQ)
    (position : Nat) (hposition : position < row.size) :
    (Generated.StrictRecombine.zeroQQRowLoop columns index row)[position]! =
      row[position] := by
  induction hremaining : columns - index using Nat.strong_induction_on
      generalizing index row with
  | h remaining ih =>
      rw [Generated.StrictRecombine.zeroQQRowLoop]
      split
      next hmore =>
        have hdecrease : columns - (index + 1) < remaining := by omega
        rw [ih (columns - (index + 1)) hdecrease (index + 1) (row.push 0)
          (by simp only [Array.size_push]; omega) rfl]
        simpa only [Array.getElem_push_lt hposition]
      next hdone =>
        simpa only [getElem!_pos row position hposition]

theorem zeroQQRowLoop_suffix_get (columns index : Nat) (row : Array QQ)
    (hindex : index ≤ columns) (offset : Nat)
    (hoffset : offset < columns - index) :
    (Generated.StrictRecombine.zeroQQRowLoop columns index row)[row.size + offset]! =
      0 := by
  induction offset generalizing index row with
  | zero =>
      rw [Generated.StrictRecombine.zeroQQRowLoop]
      split
      next hmore =>
        have hposition : row.size < (row.push (0 : QQ)).size := by simp
        simp only [Nat.add_zero]
        rw [zeroQQRowLoop_get_of_lt columns (index + 1) (row.push 0)
          row.size hposition]
        simp
      next hdone => omega
  | succ offset ih =>
      rw [Generated.StrictRecombine.zeroQQRowLoop]
      split
      next hmore =>
        have hoffset' : offset < columns - (index + 1) := by omega
        have htail := ih (index + 1) (row.push 0) (by omega) hoffset'
        simpa only [Array.size_push, Nat.add_assoc, Nat.add_left_comm,
          Nat.add_comm] using htail
      next hdone => omega

theorem zeroQQRow_get (columns position : Nat) (hposition : position < columns) :
    (Generated.StrictRecombine.zeroQQRow columns)[position]! = 0 := by
  unfold Generated.StrictRecombine.zeroQQRow
  simpa using zeroQQRowLoop_suffix_get columns 0 (#[] : Array QQ)
    (by omega) position (by omega)

theorem zeroQQMatrixLoop_size (rows columns index : Nat)
    (matrix : Generated.StrictRecombine.QQMatrix) (hindex : index ≤ rows) :
    (Generated.StrictRecombine.zeroQQMatrixLoop rows columns index matrix).size =
      matrix.size + (rows - index) := by
  induction hremaining : rows - index using Nat.strong_induction_on
      generalizing index matrix with
  | h remaining ih =>
      rw [Generated.StrictRecombine.zeroQQMatrixLoop]
      split
      next hmore =>
        have hdecrease : rows - (index + 1) < remaining := by omega
        rw [ih (rows - (index + 1)) hdecrease (index + 1)
          (matrix.push (Generated.StrictRecombine.zeroQQRow columns))
          (by omega) rfl]
        simp only [Array.size_push]
        omega
      next hdone => omega

theorem zeroQQMatrix_size (rows columns : Nat) :
    (Generated.StrictRecombine.zeroQQMatrix rows columns).size = rows := by
  unfold Generated.StrictRecombine.zeroQQMatrix
  simpa using zeroQQMatrixLoop_size rows columns 0
    (#[] : Generated.StrictRecombine.QQMatrix) (by omega)

theorem zeroQQMatrixLoop_get_of_lt (rows columns index : Nat)
    (matrix : Generated.StrictRecombine.QQMatrix) (position : Nat)
    (hposition : position < matrix.size) :
    (Generated.StrictRecombine.zeroQQMatrixLoop rows columns index matrix)[position]! =
      matrix[position] := by
  induction hremaining : rows - index using Nat.strong_induction_on
      generalizing index matrix with
  | h remaining ih =>
      rw [Generated.StrictRecombine.zeroQQMatrixLoop]
      split
      next hmore =>
        have hdecrease : rows - (index + 1) < remaining := by omega
        rw [ih (rows - (index + 1)) hdecrease (index + 1)
          (matrix.push (Generated.StrictRecombine.zeroQQRow columns))
          (by simp only [Array.size_push]; omega) rfl]
        simpa only [Array.getElem_push_lt hposition]
      next hdone =>
        simpa only [getElem!_pos matrix position hposition]

theorem zeroQQMatrixLoop_suffix_get (rows columns index : Nat)
    (matrix : Generated.StrictRecombine.QQMatrix) (hindex : index ≤ rows)
    (offset : Nat) (hoffset : offset < rows - index) :
    (Generated.StrictRecombine.zeroQQMatrixLoop rows columns index matrix)[matrix.size + offset]! =
      Generated.StrictRecombine.zeroQQRow columns := by
  induction offset generalizing index matrix with
  | zero =>
      rw [Generated.StrictRecombine.zeroQQMatrixLoop]
      split
      next hmore =>
        have hposition : matrix.size <
            (matrix.push (Generated.StrictRecombine.zeroQQRow columns)).size := by simp
        simp only [Nat.add_zero]
        rw [zeroQQMatrixLoop_get_of_lt rows columns (index + 1)
          (matrix.push (Generated.StrictRecombine.zeroQQRow columns))
          matrix.size hposition]
        simp
      next hdone => omega
  | succ offset ih =>
      rw [Generated.StrictRecombine.zeroQQMatrixLoop]
      split
      next hmore =>
        have hoffset' : offset < rows - (index + 1) := by omega
        have htail := ih (index + 1)
          (matrix.push (Generated.StrictRecombine.zeroQQRow columns))
          (by omega) hoffset'
        simpa only [Array.size_push, Nat.add_assoc, Nat.add_left_comm,
          Nat.add_comm] using htail
      next hdone => omega

theorem zeroQQMatrix_get (rows columns row : Nat) (hrow : row < rows) :
    (Generated.StrictRecombine.zeroQQMatrix rows columns)[row]! =
      Generated.StrictRecombine.zeroQQRow columns := by
  unfold Generated.StrictRecombine.zeroQQMatrix
  simpa using zeroQQMatrixLoop_suffix_get rows columns 0
    (#[] : Generated.StrictRecombine.QQMatrix) (by omega) row (by omega)

theorem zeroQQMatrix_entry (rows columns row column : Nat)
    (hrow : row < rows) (hcolumn : column < columns) :
    ((Generated.StrictRecombine.zeroQQMatrix rows columns)[row]!)[column]! = 0 := by
  rw [zeroQQMatrix_get rows columns row hrow]
  exact zeroQQRow_get columns column hcolumn

theorem zeroQQMatrixLoop_rows (rows columns index : Nat)
    (matrix : Generated.StrictRecombine.QQMatrix) (hindex : index ≤ rows)
    (hmatrixRows : ∀ row (hrow : row < matrix.size),
      matrix[row].size = columns) :
    ∀ row
      (hrow : row <
        (Generated.StrictRecombine.zeroQQMatrixLoop rows columns index matrix).size),
      (Generated.StrictRecombine.zeroQQMatrixLoop rows columns index matrix)[row].size =
        columns := by
  induction hremaining : rows - index using Nat.strong_induction_on
      generalizing index matrix with
  | h remaining ih =>
      rw [Generated.StrictRecombine.zeroQQMatrixLoop]
      split
      next hmore =>
        have hdecrease : rows - (index + 1) < remaining := by omega
        have hpushRows : ∀ row
            (hrow : row <
              (matrix.push (Generated.StrictRecombine.zeroQQRow columns)).size),
            (matrix.push
              (Generated.StrictRecombine.zeroQQRow columns))[row].size = columns := by
          intro row hrow
          by_cases hlast : row = matrix.size
          · subst row
            simp [zeroQQRow_size]
          · have hrowOld : row < matrix.size := by
              have hrow' : row < matrix.size + 1 := by
                simpa only [Array.size_push] using hrow
              omega
            simpa only [Array.getElem_push_lt hrowOld] using
              hmatrixRows row hrowOld
        exact ih (rows - (index + 1)) hdecrease (index + 1)
          (matrix.push (Generated.StrictRecombine.zeroQQRow columns))
          (by omega) hpushRows rfl
      next hdone => exact hmatrixRows

theorem zeroQQMatrix_rows (rows columns : Nat) (row : Nat)
    (hrow : row < rows) :
    (Generated.StrictRecombine.zeroQQMatrix rows columns)[row]!.size = columns := by
  have hshape := zeroQQMatrixLoop_rows rows columns 0
    (#[] : Generated.StrictRecombine.QQMatrix) (by omega)
    (by intro row hrow; simp at hrow)
  have hbound : row <
      (Generated.StrictRecombine.zeroQQMatrix rows columns).size := by
    simpa only [zeroQQMatrix_size] using hrow
  rw [getElem!_pos (Generated.StrictRecombine.zeroQQMatrix rows columns)
    row hbound]
  exact hshape row hbound

theorem dotRowsLoop_eq_Ico_sum (left right : Array ZZ)
    (index : Nat) (acc output : ZZ)
    (hrun : Generated.StrictRecombine.dotRowsLoop left right index acc =
      .ok output) :
    output = acc + ∑ k ∈ Finset.Ico index left.size, left[k]! * right[k]! := by
  induction hremaining : left.size - index using Nat.strong_induction_on
      generalizing index acc output with
  | h remaining ih =>
      rw [Generated.StrictRecombine.dotRowsLoop] at hrun
      split at hrun
      next hindex =>
        split at hrun
        next hrightIndex =>
          have hdecrease : left.size - (index + 1) < remaining := by omega
          have htail := ih (left.size - (index + 1)) hdecrease (index + 1)
            (acc + left[index] * right[index]) output hrun rfl
          rw [htail]
          rw [Finset.sum_eq_sum_Ico_succ_bot (by omega : index < left.size)]
          simp only [getElem!_pos left index hindex,
            getElem!_pos right index hrightIndex]
          ring
        next hrightIndex => contradiction
      next hindex =>
        have hempty : Finset.Ico index left.size = ∅ := by
          exact Finset.Ico_eq_empty hindex
        have hout := Except.ok.inj hrun
        subst output
        simp [hempty]

theorem dotRows_eq_fin_sum (left right : Array ZZ)
    (hright : left.size ≤ right.size) :
    Generated.StrictRecombine.dotRows left right =
      .ok (∑ k : Fin left.size, left[k.val] * right[k.val]) := by
  unfold Generated.StrictRecombine.dotRows
  have hsucceeds : ∃ output,
      Generated.StrictRecombine.dotRowsLoop left right 0 0 = .ok output := by
    have hloop : ∀ index acc, index ≤ left.size →
        ∃ output, Generated.StrictRecombine.dotRowsLoop left right index acc =
          .ok output := by
      intro index acc hindex
      induction hremaining : left.size - index using Nat.strong_induction_on
          generalizing index acc with
      | h remaining ih =>
          rw [Generated.StrictRecombine.dotRowsLoop]
          split
          next hmore =>
            split
            next hrightIndex =>
              have hdecrease : left.size - (index + 1) < remaining := by omega
              exact ih (left.size - (index + 1)) hdecrease (index + 1)
                (acc + left[index] * right[index]) (by omega) rfl
            next hrightIndex => omega
          next hmore => exact ⟨acc, rfl⟩
    exact hloop 0 0 (by omega)
  obtain ⟨output, houtput⟩ := hsucceeds
  have hrun := dotRowsLoop_eq_Ico_sum left right 0 0 output houtput
  rw [houtput, Except.ok.injEq]
  rw [hrun]
  congr 1
  simp only [zero_add]
  rw [show Finset.Ico 0 left.size = Finset.range left.size by
    ext value
    simp]
  rw [← Fin.sum_univ_eq_sum_range]
  apply Finset.sum_congr rfl
  intro index _
  simp [getElem!_pos left index.val index.isLt,
    getElem!_pos right index.val (lt_of_lt_of_le index.isLt hright)]

theorem gramNumeratorLoop_eq_Ico_sum
    (mu : Generated.StrictRecombine.QQMatrix) (norms : Array QQ)
    (i j l : Nat) (numerator output : QQ)
    (hrun : Generated.StrictRecombine.gramNumeratorLoop mu norms i j l
      numerator = .ok output) :
    output = numerator - ∑ k ∈ Finset.Ico l j,
      (mu[i]!)[k]! * (mu[j]!)[k]! * norms[k]! := by
  induction hremaining : j - l using Nat.strong_induction_on
      generalizing l numerator output with
  | h remaining ih =>
      rw [Generated.StrictRecombine.gramNumeratorLoop] at hrun
      split at hrun
      next hl =>
        split at hrun
        next hi =>
          split at hrun
          next hj =>
            split at hrun
            next hil =>
              split at hrun
              next hjl =>
                split at hrun
                next hn =>
                  have hdecrease : j - (l + 1) < remaining := by omega
                  have htail := ih (j - (l + 1)) hdecrease (l + 1)
                    (numerator - mu[i][l] * mu[j][l] * norms[l]) output
                    hrun rfl
                  rw [htail]
                  rw [Finset.sum_eq_sum_Ico_succ_bot hl]
                  simp only [getElem!_pos mu i hi, getElem!_pos mu j hj,
                    getElem!_pos mu[i] l hil, getElem!_pos mu[j] l hjl,
                    getElem!_pos norms l hn]
                  ring
                next hn => contradiction
              next hjl => contradiction
            next hil => contradiction
          next hj => contradiction
        next hi => contradiction
      next hl =>
        have hempty : Finset.Ico l j = ∅ := Finset.Ico_eq_empty hl
        have hout := Except.ok.inj hrun
        subst output
        simp [hempty]

theorem gramNumeratorLoop_succeeds
    (mu : Generated.StrictRecombine.QQMatrix) (norms : Array QQ)
    (i j l : Nat) (numerator : QQ)
    (hi : i < mu.size) (hj : j < mu.size)
    (hmuRows : ∀ row (hrow : row < mu.size), mu[row].size = mu.size)
    (hnorms : mu.size ≤ norms.size) (hl : l ≤ j) :
    ∃ output, Generated.StrictRecombine.gramNumeratorLoop mu norms i j l
      numerator = .ok output := by
  induction hremaining : j - l using Nat.strong_induction_on
      generalizing l numerator with
  | h remaining ih =>
      rw [Generated.StrictRecombine.gramNumeratorLoop]
      split
      next hmore =>
        have hil : l < mu[i].size := by rw [hmuRows i hi]; omega
        have hjl : l < mu[j].size := by rw [hmuRows j hj]; omega
        have hn : l < norms.size := lt_of_lt_of_le (by omega) hnorms
        rw [dif_pos hil, dif_pos hjl, dif_pos hn]
        have hdecrease : j - (l + 1) < remaining := by omega
        exact ih (j - (l + 1)) hdecrease (l + 1)
          (numerator - mu[i][l] * mu[j][l] * norms[l]) (by omega) rfl
      next hmore => exact ⟨numerator, rfl⟩

theorem gramNumeratorLoop_exact
    (mu : Generated.StrictRecombine.QQMatrix) (norms : Array QQ)
    (i j : Nat) (numerator : QQ)
    (hi : i < mu.size) (hj : j < mu.size)
    (hmuRows : ∀ row (hrow : row < mu.size), mu[row].size = mu.size)
    (hnorms : mu.size ≤ norms.size) :
    Generated.StrictRecombine.gramNumeratorLoop mu norms i j 0 numerator =
      .ok (numerator - ∑ k : Fin j,
        (mu[i][k.val]'(by rw [hmuRows i hi]; exact lt_trans k.isLt hj)) *
        (mu[j][k.val]'(by rw [hmuRows j hj]; exact lt_trans k.isLt hj)) *
        norms[k.val]'(lt_of_lt_of_le (lt_trans k.isLt hj) hnorms)) := by
  obtain ⟨output, houtput⟩ := gramNumeratorLoop_succeeds mu norms i j 0
    numerator hi hj hmuRows hnorms (by omega)
  have hexact := gramNumeratorLoop_eq_Ico_sum mu norms i j 0 numerator output
    houtput
  rw [houtput, Except.ok.injEq]
  rw [hexact]
  congr 1
  rw [show Finset.Ico 0 j = Finset.range j by ext value; simp]
  rw [← Fin.sum_univ_eq_sum_range]
  apply Finset.sum_congr rfl
  intro index _
  simp [getElem!_pos mu i hi, getElem!_pos mu j hj,
    getElem!_pos mu[i] index.val (by rw [hmuRows i hi]; exact lt_trans index.isLt hj),
    getElem!_pos mu[j] index.val (by rw [hmuRows j hj]; exact lt_trans index.isLt hj),
    getElem!_pos norms index.val
      (lt_of_lt_of_le (lt_trans index.isLt hj) hnorms)]

theorem gramNormLoop_eq_Ico_sum
    (mu : Generated.StrictRecombine.QQMatrix) (norms : Array QQ)
    (i j : Nat) (norm output : QQ)
    (hrun : Generated.StrictRecombine.gramNormLoop mu norms i j norm =
      .ok output) :
    output = norm - ∑ k ∈ Finset.Ico j i,
      (mu[i]!)[k]! * (mu[i]!)[k]! * norms[k]! := by
  induction hremaining : i - j using Nat.strong_induction_on
      generalizing j norm output with
  | h remaining ih =>
      rw [Generated.StrictRecombine.gramNormLoop] at hrun
      split at hrun
      next hj =>
        split at hrun
        next hi =>
          split at hrun
          next hij =>
            split at hrun
            next hn =>
              have hdecrease : i - (j + 1) < remaining := by omega
              have htail := ih (i - (j + 1)) hdecrease (j + 1)
                (norm - mu[i][j] * mu[i][j] * norms[j]) output hrun rfl
              rw [htail]
              rw [Finset.sum_eq_sum_Ico_succ_bot hj]
              simp only [getElem!_pos mu i hi, getElem!_pos mu[i] j hij,
                getElem!_pos norms j hn]
              ring
            next hn => contradiction
          next hij => contradiction
        next hi => contradiction
      next hj =>
        have hempty : Finset.Ico j i = ∅ := Finset.Ico_eq_empty hj
        have hout := Except.ok.inj hrun
        subst output
        simp [hempty]

theorem gramNormLoop_succeeds
    (mu : Generated.StrictRecombine.QQMatrix) (norms : Array QQ)
    (i j : Nat) (norm : QQ)
    (hi : i < mu.size)
    (hmuRows : ∀ row (hrow : row < mu.size), mu[row].size = mu.size)
    (hnorms : mu.size ≤ norms.size) (hj : j ≤ i) :
    ∃ output, Generated.StrictRecombine.gramNormLoop mu norms i j norm =
      .ok output := by
  induction hremaining : i - j using Nat.strong_induction_on
      generalizing j norm with
  | h remaining ih =>
      rw [Generated.StrictRecombine.gramNormLoop]
      split
      next hmore =>
        have hij : j < mu[i].size := by rw [hmuRows i hi]; omega
        have hn : j < norms.size := lt_of_lt_of_le (by omega) hnorms
        rw [dif_pos hij, dif_pos hn]
        have hdecrease : i - (j + 1) < remaining := by omega
        exact ih (i - (j + 1)) hdecrease (j + 1)
          (norm - mu[i][j] * mu[i][j] * norms[j]) (by omega) rfl
      next hmore => exact ⟨norm, rfl⟩

theorem gramNormLoop_exact
    (mu : Generated.StrictRecombine.QQMatrix) (norms : Array QQ)
    (i : Nat) (norm : QQ)
    (hi : i < mu.size)
    (hmuRows : ∀ row (hrow : row < mu.size), mu[row].size = mu.size)
    (hnorms : mu.size ≤ norms.size) :
    Generated.StrictRecombine.gramNormLoop mu norms i 0 norm =
      .ok (norm - ∑ k : Fin i,
        (mu[i][k.val]'(by rw [hmuRows i hi]; exact lt_trans k.isLt hi)) ^ 2 *
        norms[k.val]'(lt_of_lt_of_le (lt_trans k.isLt hi) hnorms)) := by
  obtain ⟨output, houtput⟩ := gramNormLoop_succeeds mu norms i 0 norm hi
    hmuRows hnorms (by omega)
  have hexact := gramNormLoop_eq_Ico_sum mu norms i 0 norm output houtput
  rw [houtput, Except.ok.injEq, hexact]
  congr 1
  rw [show Finset.Ico 0 i = Finset.range i by ext value; simp]
  rw [← Fin.sum_univ_eq_sum_range]
  apply Finset.sum_congr rfl
  intro index _
  have hindex : index.val < mu[i].size := by
    rw [hmuRows i hi]
    exact lt_trans index.isLt hi
  simp [pow_two, getElem!_pos mu i hi, getElem!_pos mu[i] index.val hindex,
    getElem!_pos norms index.val
      (lt_of_lt_of_le (lt_trans index.isLt hi) hnorms)]

/-- Array-shape invariant carried by the generated Gram--Schmidt
initialization loops. -/
structure GramStorageShape (mu : Generated.StrictRecombine.QQMatrix)
    (norms : Array QQ) (size : Nat) : Prop where
  mu_size : mu.size = size
  mu_rows : ∀ row (hrow : row < mu.size), mu[row].size = size
  norms_size : norms.size = size

theorem GramStorageShape.setMu
    {mu : Generated.StrictRecombine.QQMatrix} {norms : Array QQ} {size : Nat}
    (hshape : GramStorageShape mu norms size)
    (row column : Nat) (value : QQ)
    (hrow : row < mu.size) (hcolumn : column < mu[row].size) :
    GramStorageShape (mu.set row (mu[row].set column value)) norms size := by
  refine ⟨by simp [hshape.mu_size], ?_, hshape.norms_size⟩
  intro other hother
  by_cases heq : other = row
  · subst other
    simp [hshape.mu_rows row hrow]
  · have hotherOld : other < mu.size := by simpa using hother
    rw [Array.getElem_set_ne hrow hotherOld (Ne.symm heq)]
    exact hshape.mu_rows other hotherOld

theorem GramStorageShape.setNorm
    {mu : Generated.StrictRecombine.QQMatrix} {norms : Array QQ} {size : Nat}
    (hshape : GramStorageShape mu norms size)
    (index : Nat) (value : QQ) (hindex : index < norms.size) :
    GramStorageShape mu (norms.set index value) size := by
  refine ⟨hshape.mu_size, hshape.mu_rows, ?_⟩
  simp [hshape.norms_size]

/-- The exact scalar written by one source `j < i` Gram--Schmidt iteration,
expressed only through the cells read by the generated loops. -/
noncomputable def sourceGramCoefficient
    (matrix : Generated.StrictRecombine.LLLMatrix)
    (mu : Generated.StrictRecombine.QQMatrix) (norms : Array QQ)
    (i j rowSize : Nat) : QQ :=
  let dot : QQ := ((∑ k : Fin rowSize,
    (matrix[i]!)[k.val]! * (matrix[j]!)[k.val]! : ZZ) : QQ)
  let numerator := dot - ∑ k : Fin j,
    (mu[i]!)[k.val]! * (mu[j]!)[k.val]! * norms[k.val]!
  if norms[j]! = 0 then 0 else numerator / norms[j]!

/-- Exact rational dot product of two rows as read by the generated source.
The left row fixes the loop bound, matching `dotRows`. -/
noncomputable def sourceRowDot
    (matrix : Generated.StrictRecombine.LLLMatrix) (i j rowSize : Nat) : QQ :=
  ((∑ k : Fin rowSize,
    (matrix[i]!)[k.val]! * (matrix[j]!)[k.val]! : ZZ) : QQ)

/-- The recurrence already established for the written prefix of source row
`i`.  Column `j` states the exact lower-triangular LDL equation for the
`(i,j)` Gram entry; it talks only about cells in the generated arrays. -/
def GramMuPrefixCorrect
    (matrix : Generated.StrictRecombine.LLLMatrix)
    (mu : Generated.StrictRecombine.QQMatrix) (norms : Array QQ)
    (i columns : Nat) : Prop :=
  ∀ j, j < columns →
    sourceRowDot matrix i j matrix[i]!.size =
      (∑ k : Fin j,
        (mu[i]!)[k.val]! * (mu[j]!)[k.val]! * norms[k.val]!) +
        (mu[i]!)[j]! * norms[j]!

theorem gramMuPrefixCorrect_zero
    (matrix : Generated.StrictRecombine.LLLMatrix)
    (mu : Generated.StrictRecombine.QQMatrix) (norms : Array QQ) (i : Nat) :
    GramMuPrefixCorrect matrix mu norms i 0 := by
  intro j hj
  omega

theorem sourceGramCoefficient_closes_column
    (matrix : Generated.StrictRecombine.LLLMatrix)
    (mu : Generated.StrictRecombine.QQMatrix) (norms : Array QQ)
    (i j rowSize : Nat)
    (hnorm : norms[j]! ≠ 0) :
    sourceRowDot matrix i j rowSize =
      (∑ k : Fin j,
        (mu[i]!)[k.val]! * (mu[j]!)[k.val]! * norms[k.val]!) +
        sourceGramCoefficient matrix mu norms i j rowSize * norms[j]! := by
  unfold sourceRowDot sourceGramCoefficient
  rw [if_neg hnorm]
  field_simp [hnorm]
  ring

theorem gramMuPrefixCorrect_set_next
    (matrix : Generated.StrictRecombine.LLLMatrix)
    (mu : Generated.StrictRecombine.QQMatrix) (norms : Array QQ)
    (i j : Nat) (coefficient : QQ)
    (hshape : GramStorageShape mu norms matrix.size)
    (hi : i < matrix.size) (hj : j < i)
    (hiMu : i < mu.size) (hjMu : j < mu.size)
    (hjRow : j < mu[i].size)
    (hprefix : GramMuPrefixCorrect matrix mu norms i j)
    (hcolumn : sourceRowDot matrix i j matrix[i]!.size =
      (∑ k : Fin j,
        (mu[i]!)[k.val]! * (mu[j]!)[k.val]! * norms[k.val]!) +
        coefficient * norms[j]!) :
    let mu' := mu.set i (mu[i].set j coefficient hjRow) hiMu
    GramMuPrefixCorrect matrix mu' norms i (j + 1) := by
  dsimp only
  let changedRow := mu[i].set j coefficient hjRow
  let mu' := mu.set i changedRow hiMu
  have hmuI : mu'[i]! = changedRow := by
    rw [getElem!_pos mu' i (by simp [mu', hiMu])]
    simp [mu']
  have hmuOther : ∀ row, row < mu.size → row ≠ i → mu'[row]! = mu[row]! := by
    intro row hrow hne
    rw [getElem!_pos mu' row (by simpa [mu'] using hrow)]
    rw [Array.getElem_set_ne hiMu hrow (Ne.symm hne)]
    exact (getElem!_pos mu row hrow).symm
  have hchangedPrefix : ∀ column, column < j →
      changedRow[column]! = (mu[i]!)[column]! := by
    intro column hcolumn
    have hcolumnRow : column < mu[i].size := lt_trans hcolumn hjRow
    rw [getElem!_pos changedRow column (by simpa [changedRow] using hcolumnRow)]
    rw [Array.getElem_set_ne hjRow hcolumnRow (by omega)]
    rw [getElem!_pos mu i hiMu]
    exact (getElem!_pos mu[i] column hcolumnRow).symm
  intro column hcolumnBound
  by_cases hcurrent : column = j
  · subst column
    rw [hcolumn]
    congr 1
    · apply Fintype.sum_congr
      intro k
      have hkI : k.val < mu[i].size := lt_trans k.isLt hjRow
      have hkJ : k.val < mu[j].size := by
        rw [hshape.mu_rows j hjMu]
        omega
      rw [hmuI]
      rw [hmuOther j hjMu (by omega)]
      rw [hchangedPrefix k.val k.isLt]
    · rw [hmuI]
      simp [changedRow, getElem!_pos, hjRow]
  · have hprevious : column < j := by omega
    rw [hprefix column hprevious]
    have hcolumnMu : column < mu.size := lt_trans hprevious hjMu
    have hcolumnRow : column < mu[i].size := lt_trans hprevious hjRow
    have hcolumnNorm : column < norms.size := by
      rw [hshape.norms_size]
      omega
    congr 1
    · apply Fintype.sum_congr
      intro k
      have hkI : k.val < mu[i].size := lt_trans k.isLt hcolumnRow
      have hkColumn : k.val < mu[column].size := by
        rw [hshape.mu_rows column hcolumnMu]
        omega
      rw [hmuI]
      rw [hmuOther column hcolumnMu (by omega)]
      rw [hchangedPrefix k.val (lt_trans k.isLt hprevious)]
    · rw [hmuI]
      rw [hchangedPrefix column hprevious]

theorem gramMuPrefixCorrect_set_norm_at_end
    (matrix : Generated.StrictRecombine.LLLMatrix)
    (mu : Generated.StrictRecombine.QQMatrix) (norms : Array QQ)
    (i columns : Nat) (value : QQ)
    (hcolumns : columns ≤ i) (hiNorm : i < norms.size)
    (hprefix : GramMuPrefixCorrect matrix mu norms i columns) :
    GramMuPrefixCorrect matrix mu (norms.set i value hiNorm) i columns := by
  have hframe : ∀ index, index < i →
      (norms.set i value hiNorm)[index]! = norms[index]! := by
    intro index hindex
    rw [getElem!_pos (norms.set i value hiNorm) index (by simpa using
      (lt_trans hindex hiNorm))]
    rw [Array.getElem_set_ne hiNorm (lt_trans hindex hiNorm) (by omega)]
    exact (getElem!_pos norms index (lt_trans hindex hiNorm)).symm
  intro column hcolumn
  rw [hprefix column hcolumn]
  congr 1
  · apply Fintype.sum_congr
    intro index
    rw [hframe index.val (by omega)]
  · rw [hframe column (lt_of_lt_of_le hcolumn hcolumns)]

theorem gramMuRowLoop_succeeds
    (matrix : Generated.StrictRecombine.LLLMatrix) (i j size : Nat)
    (mu : Generated.StrictRecombine.QQMatrix) (norms : Array QQ)
    (hmatrixSize : matrix.size = size)
    (hmatrixRows : ∀ row (hrow : row < matrix.size),
      matrix[row].size = size)
    (hshape : GramStorageShape mu norms size)
    (hi : i < size) (hj : j ≤ i) :
    ∃ output,
      Generated.StrictRecombine.gramMuRowLoop matrix i j mu norms = .ok output ∧
      GramStorageShape output.1 output.2 size ∧ output.2 = norms ∧
      (∀ row (hrow : row < mu.size), row ≠ i →
        output.1[row]! = mu[row]) ∧
      (∀ column, column < j →
        (output.1[i]!)[column]! = (mu[i]!)[column]!) := by
  induction hremaining : i - j using Nat.strong_induction_on
      generalizing j mu norms with
  | h remaining ih =>
      rw [Generated.StrictRecombine.gramMuRowLoop]
      split
      next hmore =>
        have hiM : i < matrix.size := by omega
        have hjM : j < matrix.size := by omega
        rw [dif_pos hiM, dif_pos hjM]
        have hrowI : matrix[i].size = size := hmatrixRows i hiM
        have hrowJ : matrix[j].size = size := hmatrixRows j hjM
        have hdotSize : matrix[i].size ≤ matrix[j].size := by
          rw [hrowI, hrowJ]
        rw [dotRows_eq_fin_sum matrix[i] matrix[j] hdotSize]
        have hiMu : i < mu.size := by rw [hshape.mu_size]; exact hi
        have hjMu : j < mu.size := by rw [hshape.mu_size]; omega
        have hmuRows : ∀ row (hrow : row < mu.size),
            mu[row].size = mu.size := by
          intro row hrow
          rw [hshape.mu_rows row hrow, hshape.mu_size]
        have hnorms : mu.size ≤ norms.size := by
          rw [hshape.mu_size, hshape.norms_size]
        obtain ⟨numerator, hnumerator⟩ := gramNumeratorLoop_succeeds mu norms
          i j 0 ((∑ k : Fin matrix[i].size,
            matrix[i][k.val] * matrix[j][k.val] : ZZ) : QQ)
          hiMu hjMu hmuRows hnorms (by omega)
        simp only [hnumerator]
        have hjMuRow : j < mu[i].size := by
          rw [hshape.mu_rows i hiMu]
          omega
        have hjNorm : j < norms.size := by omega
        rw [dif_pos hiMu, dif_pos hjMuRow, dif_pos hjNorm]
        let coefficient := if norms[j] = 0 then 0 else numerator / norms[j]
        let mu' := mu.set i (mu[i].set j coefficient)
        have hshape' : GramStorageShape mu' norms size :=
          hshape.setMu i j coefficient hiMu hjMuRow
        have hdecrease : i - (j + 1) < remaining := by omega
        obtain ⟨output, hrun, houtShape, houtNorms, houtRows, houtPrefix⟩ :=
          ih (i - (j + 1)) hdecrease (j + 1) mu' norms hshape'
            (by omega) rfl
        refine ⟨output, hrun, houtShape, houtNorms, ?_, ?_⟩
        · intro row hrow hrowI
          have hrowMu' : row < mu'.size := by
            simpa [mu'] using hrow
          rw [houtRows row hrowMu' hrowI]
          simp only [mu']
          rw [Array.getElem_set_ne hiMu hrow (Ne.symm hrowI)]
        · intro column hcolumn
          have hcolumnMu : column < mu[i].size := by omega
          rw [houtPrefix column (by omega)]
          simp only [mu']
          rw [getElem!_pos (mu.set i (mu[i].set j coefficient)) i (by simpa)]
          rw [Array.getElem_set_self]
          rw [getElem!_pos (mu[i].set j coefficient) column (by simpa)]
          rw [Array.getElem_set_ne hjMuRow hcolumnMu (by omega)]
          rw [getElem!_pos mu i hiMu, getElem!_pos mu[i] column hcolumnMu]
      next hmore =>
        refine ⟨(mu, norms), rfl, hshape, rfl, ?_, ?_⟩
        · intro row hrow hrowI
          exact getElem!_pos mu row hrow
        · intro column hcolumn
          rfl

theorem gramMuRowLoop_written_coefficient
    (matrix : Generated.StrictRecombine.LLLMatrix) (i j size : Nat)
    (mu : Generated.StrictRecombine.QQMatrix) (norms : Array QQ)
    (output : Generated.StrictRecombine.QQMatrix × Array QQ)
    (hmatrixSize : matrix.size = size)
    (hmatrixRows : ∀ row (hrow : row < matrix.size),
      matrix[row].size = size)
    (hshape : GramStorageShape mu norms size)
    (hi : i < size) (hj : j < i)
    (hrun : Generated.StrictRecombine.gramMuRowLoop matrix i j mu norms =
      .ok output) :
    (output.1[i]!)[j]! =
      sourceGramCoefficient matrix mu norms i j matrix[i].size := by
  rw [Generated.StrictRecombine.gramMuRowLoop] at hrun
  rw [dif_pos hj] at hrun
  have hiM : i < matrix.size := by omega
  have hjM : j < matrix.size := by omega
  rw [dif_pos hiM, dif_pos hjM] at hrun
  have hrowI : matrix[i].size = size := hmatrixRows i hiM
  have hrowJ : matrix[j].size = size := hmatrixRows j hjM
  have hdotSize : matrix[i].size ≤ matrix[j].size := by rw [hrowI, hrowJ]
  rw [dotRows_eq_fin_sum matrix[i] matrix[j] hdotSize] at hrun
  simp only at hrun
  have hiMu : i < mu.size := by rw [hshape.mu_size]; exact hi
  have hjMu : j < mu.size := by rw [hshape.mu_size]; omega
  have hmuRows : ∀ row (hrow : row < mu.size),
      mu[row].size = mu.size := by
    intro row hrow
    rw [hshape.mu_rows row hrow, hshape.mu_size]
  have hnorms : mu.size ≤ norms.size := by
    rw [hshape.mu_size, hshape.norms_size]
  rw [gramNumeratorLoop_exact mu norms i j
    ((∑ k : Fin matrix[i].size,
      matrix[i][k.val] * matrix[j][k.val] : ZZ) : QQ)
    hiMu hjMu hmuRows hnorms] at hrun
  simp only at hrun
  have hjMuRow : j < mu[i].size := by
    rw [hshape.mu_rows i hiMu]
    omega
  have hjNorm : j < norms.size := by omega
  rw [dif_pos hiMu, dif_pos hjMuRow, dif_pos hjNorm] at hrun
  let numerator : QQ :=
    ((∑ k : Fin matrix[i].size,
      matrix[i][k.val] * matrix[j][k.val] : ZZ) : QQ) -
      ∑ k : Fin j,
        mu[i][k.val]'(by rw [hmuRows i hiMu]; exact lt_trans k.isLt hjMu) *
        mu[j][k.val]'(by rw [hmuRows j hjMu]; exact lt_trans k.isLt hjMu) *
        norms[k.val]'(lt_of_lt_of_le (lt_trans k.isLt hjMu) hnorms)
  let coefficient := if norms[j] = 0 then 0 else numerator / norms[j]
  let mu' := mu.set i (mu[i].set j coefficient)
  have hshape' : GramStorageShape mu' norms size :=
    hshape.setMu i j coefficient hiMu hjMuRow
  obtain ⟨tailOutput, htail, htailShape, htailNorms, htailRows,
      htailPrefix⟩ :=
    gramMuRowLoop_succeeds matrix i (j + 1) size mu' norms hmatrixSize
      hmatrixRows hshape' hi (by omega)
  rw [hrun] at htail
  have htailEq := Except.ok.inj htail
  subst tailOutput
  rw [htailPrefix j (by omega)]
  simp only [mu']
  rw [getElem!_pos (mu.set i (mu[i].set j coefficient)) i (by simpa)]
  rw [Array.getElem_set_self]
  rw [getElem!_pos (mu[i].set j coefficient) j (by simpa)]
  rw [Array.getElem_set_self]
  have hdotEq :
      ((∑ k : Fin matrix[i].size,
        matrix[i][k.val] * matrix[j][k.val] : ZZ) : QQ) =
      ((∑ k : Fin matrix[i].size,
        (matrix[i]!)[k.val]! * (matrix[j]!)[k.val]! : ZZ) : QQ) := by
    simp only [getElem!_pos matrix i hiM, getElem!_pos matrix j hjM]
    congr 1
    apply Fintype.sum_congr
    intro index
    simp [getElem!_pos matrix[i] index.val index.isLt,
      getElem!_pos matrix[j] index.val
        (lt_of_lt_of_le index.isLt hdotSize)]
  have hmuSumEq :
      (∑ k : Fin j,
        mu[i][k.val]'(by rw [hmuRows i hiMu]; exact lt_trans k.isLt hjMu) *
        mu[j][k.val]'(by rw [hmuRows j hjMu]; exact lt_trans k.isLt hjMu) *
        norms[k.val]'(lt_of_lt_of_le (lt_trans k.isLt hjMu) hnorms)) =
      ∑ k : Fin j,
        (mu[i]!)[k.val]! * (mu[j]!)[k.val]! * norms[k.val]! := by
    apply Fintype.sum_congr
    intro index
    have hki : index.val < mu[i].size := by
      rw [hmuRows i hiMu]
      exact lt_trans index.isLt hjMu
    have hkj : index.val < mu[j].size := by
      rw [hmuRows j hjMu]
      exact lt_trans index.isLt hjMu
    have hkn : index.val < norms.size :=
      lt_of_lt_of_le (lt_trans index.isLt hjMu) hnorms
    simp [getElem!_pos mu i hiMu, getElem!_pos mu j hjMu,
      getElem!_pos mu[i] index.val hki, getElem!_pos mu[j] index.val hkj,
      getElem!_pos norms index.val hkn]
  unfold sourceGramCoefficient
  simp only [coefficient, numerator]
  rw [getElem!_pos norms j hjNorm]
  rw [hdotEq, hmuSumEq]

theorem gramMuRowLoop_prefix_correct
    (matrix : Generated.StrictRecombine.LLLMatrix) (i j size : Nat)
    (mu : Generated.StrictRecombine.QQMatrix) (norms : Array QQ)
    (output : Generated.StrictRecombine.QQMatrix × Array QQ)
    (hmatrixSize : matrix.size = size)
    (hmatrixRows : ∀ row (hrow : row < matrix.size),
      matrix[row].size = size)
    (hshape : GramStorageShape mu norms size)
    (hi : i < size) (hj : j ≤ i)
    (hpositive : ∀ column, column < i → 0 < norms[column]!)
    (hprefix : GramMuPrefixCorrect matrix mu norms i j)
    (hrun : Generated.StrictRecombine.gramMuRowLoop matrix i j mu norms =
      .ok output) :
    GramMuPrefixCorrect matrix output.1 output.2 i i := by
  induction hremaining : i - j using Nat.strong_induction_on
      generalizing j mu norms with
  | h remaining ih =>
      rw [Generated.StrictRecombine.gramMuRowLoop] at hrun
      split at hrun
      next hmore =>
        have hiM : i < matrix.size := by omega
        have hjM : j < matrix.size := by omega
        rw [dif_pos hiM, dif_pos hjM] at hrun
        have hrowI : matrix[i].size = size := hmatrixRows i hiM
        have hrowJ : matrix[j].size = size := hmatrixRows j hjM
        have hdotSize : matrix[i].size ≤ matrix[j].size := by
          rw [hrowI, hrowJ]
        rw [dotRows_eq_fin_sum matrix[i] matrix[j] hdotSize] at hrun
        simp only at hrun
        have hiMu : i < mu.size := by rw [hshape.mu_size]; exact hi
        have hjMu : j < mu.size := by rw [hshape.mu_size]; omega
        have hmuRows : ∀ row (hrow : row < mu.size),
            mu[row].size = mu.size := by
          intro row hrow
          rw [hshape.mu_rows row hrow, hshape.mu_size]
        have hnorms : mu.size ≤ norms.size := by
          rw [hshape.mu_size, hshape.norms_size]
        rw [gramNumeratorLoop_exact mu norms i j
          ((∑ k : Fin matrix[i].size,
            matrix[i][k.val] * matrix[j][k.val] : ZZ) : QQ)
          hiMu hjMu hmuRows hnorms] at hrun
        simp only at hrun
        have hjMuRow : j < mu[i].size := by
          rw [hshape.mu_rows i hiMu]
          omega
        have hjNorm : j < norms.size := by omega
        rw [dif_pos hiMu, dif_pos hjMuRow, dif_pos hjNorm] at hrun
        let numerator : QQ :=
          ((∑ k : Fin matrix[i].size,
            matrix[i][k.val] * matrix[j][k.val] : ZZ) : QQ) -
            ∑ k : Fin j,
              mu[i][k.val]'(by rw [hmuRows i hiMu]; exact lt_trans k.isLt hjMu) *
              mu[j][k.val]'(by rw [hmuRows j hjMu]; exact lt_trans k.isLt hjMu) *
              norms[k.val]'(lt_of_lt_of_le (lt_trans k.isLt hjMu) hnorms)
        let coefficient := if norms[j] = 0 then 0 else numerator / norms[j]
        let mu' := mu.set i (mu[i].set j coefficient)
        have hshape' : GramStorageShape mu' norms size :=
          hshape.setMu i j coefficient hiMu hjMuRow
        have hdotEq :
            ((∑ k : Fin matrix[i].size,
              matrix[i][k.val] * matrix[j][k.val] : ZZ) : QQ) =
            ((∑ k : Fin matrix[i].size,
              (matrix[i]!)[k.val]! * (matrix[j]!)[k.val]! : ZZ) : QQ) := by
          simp only [getElem!_pos matrix i hiM, getElem!_pos matrix j hjM]
          congr 1
          apply Fintype.sum_congr
          intro index
          simp [getElem!_pos matrix[i] index.val index.isLt,
            getElem!_pos matrix[j] index.val
              (lt_of_lt_of_le index.isLt hdotSize)]
        have hmuSumEq :
            (∑ k : Fin j,
              mu[i][k.val]'(by rw [hmuRows i hiMu]; exact lt_trans k.isLt hjMu) *
              mu[j][k.val]'(by rw [hmuRows j hjMu]; exact lt_trans k.isLt hjMu) *
              norms[k.val]'(lt_of_lt_of_le (lt_trans k.isLt hjMu) hnorms)) =
            ∑ k : Fin j,
              (mu[i]!)[k.val]! * (mu[j]!)[k.val]! * norms[k.val]! := by
          apply Fintype.sum_congr
          intro index
          have hki : index.val < mu[i].size := by
            rw [hmuRows i hiMu]
            exact lt_trans index.isLt hjMu
          have hkj : index.val < mu[j].size := by
            rw [hmuRows j hjMu]
            exact lt_trans index.isLt hjMu
          have hkn : index.val < norms.size :=
            lt_of_lt_of_le (lt_trans index.isLt hjMu) hnorms
          simp [getElem!_pos mu i hiMu, getElem!_pos mu j hjMu,
            getElem!_pos mu[i] index.val hki,
            getElem!_pos mu[j] index.val hkj,
            getElem!_pos norms index.val hkn]
        have hcoefficient : coefficient =
            sourceGramCoefficient matrix mu norms i j matrix[i].size := by
          unfold coefficient sourceGramCoefficient numerator
          rw [getElem!_pos norms j hjNorm]
          rw [hdotEq, hmuSumEq]
        have hnormNe : norms[j]! ≠ 0 :=
          ne_of_gt (hpositive j (by omega))
        have hcolumn := sourceGramCoefficient_closes_column matrix mu norms
          i j matrix[i].size hnormNe
        rw [← hcoefficient] at hcolumn
        have hprefix' : GramMuPrefixCorrect matrix mu' norms i (j + 1) := by
          exact gramMuPrefixCorrect_set_next matrix mu norms i j coefficient
            (by simpa [hmatrixSize] using hshape) (by omega) (by omega)
            hiMu hjMu hjMuRow hprefix (by
              simpa [getElem!_pos matrix i hiM] using hcolumn)
        have hdecrease : i - (j + 1) < remaining := by omega
        exact ih (i - (j + 1)) hdecrease (j + 1) mu' norms hshape'
          (by omega) hpositive hprefix' hrun rfl
      next hmore =>
        have houtput := Except.ok.inj hrun
        subst output
        intro column hcolumn
        exact hprefix column (by omega)

/-- Exact diagonal equation produced by the generated `gramNormLoop` after
the μ row has been completed. -/
theorem gramNormLoop_closes_diagonal
    (matrix : Generated.StrictRecombine.LLLMatrix)
    (mu : Generated.StrictRecombine.QQMatrix) (norms : Array QQ)
    (i size : Nat) (norm : QQ)
    (hmatrixSize : matrix.size = size)
    (hmatrixRows : ∀ row (hrow : row < matrix.size),
      matrix[row].size = size)
    (hshape : GramStorageShape mu norms size) (hi : i < size)
    (hrun : Generated.StrictRecombine.gramNormLoop mu norms i 0
      (((∑ k : Fin matrix[i].size,
        matrix[i][k.val] * matrix[i][k.val] : ZZ) : QQ)) = .ok norm) :
    sourceRowDot matrix i i matrix[i]!.size =
      (∑ k : Fin i,
        (mu[i]!)[k.val]! * (mu[i]!)[k.val]! * norms[k.val]!) + norm := by
  have hiM : i < matrix.size := by omega
  have hiMu : i < mu.size := by rw [hshape.mu_size]; exact hi
  have hmuRows : ∀ row (hrow : row < mu.size),
      mu[row].size = mu.size := by
    intro row hrow
    rw [hshape.mu_rows row hrow, hshape.mu_size]
  have hnorms : mu.size ≤ norms.size := by
    rw [hshape.mu_size, hshape.norms_size]
  have hexact := gramNormLoop_exact mu norms i
    (((∑ k : Fin matrix[i].size,
      matrix[i][k.val] * matrix[i][k.val] : ZZ) : QQ))
    hiMu hmuRows hnorms
  rw [hexact] at hrun
  have hnorm := Except.ok.inj hrun
  subst norm
  have hdotEq :
      ((∑ k : Fin matrix[i].size,
        matrix[i][k.val] * matrix[i][k.val] : ZZ) : QQ) =
      sourceRowDot matrix i i matrix[i]!.size := by
    unfold sourceRowDot
    rw [getElem!_pos matrix i hiM]
    congr 1
    apply Fintype.sum_congr
    intro index
    simp [getElem!_pos matrix[i] index.val index.isLt]
  have hsumEq :
      (∑ k : Fin i,
        mu[i][k.val]'(by rw [hmuRows i hiMu]; exact lt_trans k.isLt hiMu) ^ 2 *
          norms[k.val]'(lt_of_lt_of_le (lt_trans k.isLt hiMu) hnorms)) =
      ∑ k : Fin i,
        (mu[i]!)[k.val]! * (mu[i]!)[k.val]! * norms[k.val]! := by
    apply Fintype.sum_congr
    intro index
    have hki : index.val < mu[i].size := by
      rw [hmuRows i hiMu]
      exact lt_trans index.isLt hiMu
    have hkn : index.val < norms.size :=
      lt_of_lt_of_le (lt_trans index.isLt hiMu) hnorms
    simp [pow_two, getElem!_pos mu i hiMu,
      getElem!_pos mu[i] index.val hki,
      getElem!_pos norms index.val hkn]
  rw [← hdotEq, ← hsumEq]
  ring

theorem initializeGramSchmidtLoop_succeeds
    (matrix : Generated.StrictRecombine.LLLMatrix) (i size : Nat)
    (mu : Generated.StrictRecombine.QQMatrix) (norms : Array QQ)
    (hmatrixSize : matrix.size = size)
    (hmatrixRows : ∀ row (hrow : row < matrix.size),
      matrix[row].size = size)
    (hshape : GramStorageShape mu norms size) (hi : i ≤ size) :
    ∃ output,
      Generated.StrictRecombine.initializeGramSchmidtLoop matrix i mu norms =
        .ok output ∧
      GramStorageShape output.1 output.2 size := by
  induction hremaining : size - i using Nat.strong_induction_on
      generalizing i mu norms with
  | h remaining ih =>
      rw [Generated.StrictRecombine.initializeGramSchmidtLoop]
      split
      next hiMatrix =>
        have hiSize : i < size := by omega
        obtain ⟨muOutput, hmuRun, hmuShape, hnormsUnchanged, hmuRowsFrame,
          hmuPrefixFrame⟩ :=
          gramMuRowLoop_succeeds matrix i 0 size mu norms hmatrixSize
            hmatrixRows hshape hiSize (by omega)
        rw [hmuRun]
        simp only
        subst hnormsUnchanged
        have hrowSize := hmatrixRows i hiMatrix
        rw [dotRows_eq_fin_sum matrix[i] matrix[i] (by omega)]
        simp only
        have hiMu : i < muOutput.1.size := by
          rw [hmuShape.mu_size]
          exact hiSize
        have hmuRows : ∀ row (hrow : row < muOutput.1.size),
            muOutput.1[row].size = muOutput.1.size := by
          intro row hrow
          rw [hmuShape.mu_rows row hrow, hmuShape.mu_size]
        have hnormsSize : muOutput.1.size ≤ muOutput.2.size := by
          rw [hmuShape.mu_size, hmuShape.norms_size]
        obtain ⟨normOutput, hnormRun⟩ := gramNormLoop_succeeds muOutput.1
          muOutput.2 i 0
          ((∑ k : Fin matrix[i].size,
            matrix[i][k.val] * matrix[i][k.val] : ZZ) : QQ)
          hiMu hmuRows hnormsSize (by omega)
        simp only [hnormRun]
        have hiNorm : i < muOutput.2.size := by
          rw [hmuShape.norms_size]
          exact hiSize
        rw [dif_pos hiNorm]
        let norms' := muOutput.2.set i normOutput
        have hshape' : GramStorageShape muOutput.1 norms' size :=
          hmuShape.setNorm i normOutput hiNorm
        have hdecrease : size - (i + 1) < remaining := by omega
        obtain ⟨output, hrun, houtShape⟩ := ih (size - (i + 1)) hdecrease
          (i + 1) muOutput.1 norms' hshape' (by omega) rfl
        exact ⟨output, hrun, houtShape⟩
      next hiMatrix =>
        exact ⟨(mu, norms), rfl, hshape⟩

private def positionalCode (base : Nat) : List Nat → Nat
  | [] => 0
  | digit :: rest => digit * base ^ rest.length + positionalCode base rest

private theorem positionalCode_lt_pow (base : Nat) (hbase : 0 < base)
    (digits : List Nat) (hdigits : ∀ digit ∈ digits, digit < base) :
    positionalCode base digits < base ^ digits.length := by
  induction digits with
  | nil => simp [positionalCode, hbase]
  | cons digit rest ih =>
      have hdigit := hdigits digit (by simp)
      have hrest : ∀ value ∈ rest, value < base := by
        intro value hvalue
        exact hdigits value (by simp [hvalue])
      have htail := ih hrest
      simp only [positionalCode, List.length_cons, pow_succ]
      have hpow : 0 < base ^ rest.length := Nat.pow_pos hbase
      nlinarith

private theorem positionalCode_append (base : Nat) (leading suffix : List Nat) :
    positionalCode base (leading ++ suffix) =
      positionalCode base leading * base ^ suffix.length +
        positionalCode base suffix := by
  induction leading with
  | nil => simp [positionalCode]
  | cons digit rest ih =>
      simp only [List.cons_append, positionalCode, List.length_append,
        List.length_cons]
      rw [ih, pow_add]
      ring

private theorem positionalCode_pivot_lt (base : Nat) (hbase : 0 < base)
    (leading leftSuffix rightSuffix : List Nat) (leftPivot rightPivot : Nat)
    (hlength : leftSuffix.length = rightSuffix.length)
    (hpivot : leftPivot < rightPivot)
    (hleftDigits : ∀ digit ∈ leftSuffix, digit < base)
    (hrightDigits : ∀ digit ∈ rightSuffix, digit < base) :
    positionalCode base (leading ++ leftPivot :: leftSuffix) <
      positionalCode base (leading ++ rightPivot :: rightSuffix) := by
  rw [positionalCode_append, positionalCode_append]
  simp only [positionalCode, List.length_cons, hlength]
  have hleftTail := positionalCode_lt_pow base hbase leftSuffix hleftDigits
  have hleftTail' : positionalCode base leftSuffix < base ^ rightSuffix.length := by
    simpa [hlength] using hleftTail
  apply Nat.add_lt_add_left
  calc
    leftPivot * base ^ rightSuffix.length + positionalCode base leftSuffix <
        leftPivot * base ^ rightSuffix.length + base ^ rightSuffix.length :=
      Nat.add_lt_add_left hleftTail' _
    _ = (leftPivot + 1) * base ^ rightSuffix.length := by ring
    _ ≤ rightPivot * base ^ rightSuffix.length :=
      Nat.mul_le_mul_right _ (by omega)
    _ ≤ rightPivot * base ^ rightSuffix.length + positionalCode base rightSuffix :=
      Nat.le_add_right _ _

private theorem array_toList_pivot (values : Array Nat) (pivot : Nat)
    (hpivot : pivot < values.size) :
    values.toList = values.toList.take pivot ++ values[pivot]! ::
      values.toList.drop (pivot + 1) := by
  calc
    values.toList = values.toList.take pivot ++ values.toList.drop pivot :=
      (List.take_append_drop pivot values.toList).symm
    _ = values.toList.take pivot ++ values.toList[pivot] ::
        values.toList.drop (pivot + 1) := by
      rw [List.drop_eq_getElem_cons (by simpa using hpivot)]
    _ = values.toList.take pivot ++ values[pivot]! ::
        values.toList.drop (pivot + 1) := by
      congr 2
      rw [getElem!_pos values pivot hpivot]
      exact Array.getElem_toList hpivot

private theorem positionalCode_array_pivot_lt (base : Nat) (hbase : 0 < base)
    (left right : Array Nat) (pivot : Nat)
    (hpivotLeft : pivot < left.size) (hpivotRight : pivot < right.size)
    (hsize : left.size = right.size)
    (hprefix : left.toList.take pivot = right.toList.take pivot)
    (hpivot : left[pivot]! < right[pivot]!)
    (hleftDigits : ∀ digit ∈ left.toList, digit < base)
    (hrightDigits : ∀ digit ∈ right.toList, digit < base) :
    positionalCode base left.toList < positionalCode base right.toList := by
  rw [array_toList_pivot left pivot hpivotLeft,
    array_toList_pivot right pivot hpivotRight, hprefix]
  apply positionalCode_pivot_lt base hbase
  · simp [hsize]
  · exact hpivot
  · intro digit hdigit
    exact hleftDigits digit (List.mem_of_mem_drop hdigit)
  · intro digit hdigit
    exact hrightDigits digit (List.mem_of_mem_drop hdigit)

theorem initialCombinationLoop_toList (count index : Nat)
    (result : Array Nat) :
    (Generated.StrictRecombine.initialCombinationLoop count index result).toList =
      result.toList ++ List.range' index (count - index) := by
  rw [Generated.StrictRecombine.initialCombinationLoop]
  split
  next hmore =>
    rw [initialCombinationLoop_toList]
    rw [Array.toList_push, List.append_assoc]
    congr 1
    rw [show count - index = (count - (index + 1)) + 1 by omega,
      List.range'_succ]
    simp
  next hdone =>
    have hzero : count - index = 0 := Nat.sub_eq_zero_of_le (by omega)
    simp [hzero]
termination_by count - index
decreasing_by omega

theorem initialCombination_toList (count : Nat) :
    (Generated.StrictRecombine.initialCombination count).toList =
      List.range count := by
  simp [Generated.StrictRecombine.initialCombination,
    initialCombinationLoop_toList, List.range_eq_range']

theorem initialCombination_size (count : Nat) :
    (Generated.StrictRecombine.initialCombination count).size = count := by
  simpa using congrArg List.length (initialCombination_toList count)

theorem resetCombinationSuffix_size (indices : Array Nat)
    (pivot offset : Nat) :
    (Generated.StrictRecombine.resetCombinationSuffix indices pivot offset).size =
      indices.size := by
  rw [Generated.StrictRecombine.resetCombinationSuffix]
  split
  next hposition =>
    rw [resetCombinationSuffix_size]
    simp
  next hposition => rfl
termination_by indices.size - (pivot + 1 + offset)
decreasing_by simp only [Array.size_set]; omega

theorem resetCombinationSuffix_getElem_le (indices : Array Nat)
    (pivot offset index : Nat) (hindex : index < indices.size)
    (hle : index ≤ pivot) :
    (Generated.StrictRecombine.resetCombinationSuffix indices pivot offset)[index]! =
      indices[index]! := by
  rw [Generated.StrictRecombine.resetCombinationSuffix]
  split
  next hposition =>
    rw [resetCombinationSuffix_getElem_le
      (indices.set (pivot + 1 + offset) (indices[pivot + offset] + 1)) pivot
      (offset + 1) index (by simpa using hindex) hle]
    have hne : index ≠ pivot + 1 + offset := by omega
    rw [getElem!_pos _ index (by simpa using hindex),
      getElem!_pos indices index hindex]
    rw [Array.getElem_set]
    rw [if_neg (by omega)]
  next hposition => rfl
termination_by indices.size - (pivot + 1 + offset)
decreasing_by simp only [Array.size_set]; omega

/-- Cells already lying before the next suffix-write position are preserved
by all later source writes. -/
theorem resetCombinationSuffix_getElem_before_position (indices : Array Nat)
    (pivot offset index : Nat) (hindex : index < indices.size)
    (hbefore : index < pivot + 1 + offset) :
    (Generated.StrictRecombine.resetCombinationSuffix
      indices pivot offset)[index]! = indices[index]! := by
  rw [Generated.StrictRecombine.resetCombinationSuffix]
  split
  next hposition =>
    rw [resetCombinationSuffix_getElem_before_position
      (indices.set (pivot + 1 + offset)
        (indices[pivot + offset]'(by omega) + 1))
      pivot (offset + 1) index (by simpa using hindex) (by omega)]
    rw [getElem!_pos _ index (by simpa using hindex), Array.getElem_set]
    rw [if_neg (by omega)]
    rw [getElem!_pos indices index hindex]
  next hposition => rfl
termination_by indices.size - (pivot + 1 + offset)
decreasing_by simp only [Array.size_set]; omega

/-- Every not-yet-reset suffix cell receives the exact smallest value allowed
after the pivot: the previous reset cell plus one.  This is the source-level
minimal-suffix fact needed to identify `nextCombination` as an immediate
lexicographic successor. -/
theorem resetCombinationSuffix_getElem_suffix (indices : Array Nat)
    (pivot offset index : Nat)
    (hpivotOffset : pivot + offset < indices.size)
    (hindex : index < indices.size)
    (hsuffix : pivot + 1 + offset ≤ index) :
    (Generated.StrictRecombine.resetCombinationSuffix
      indices pivot offset)[index]! =
      indices[pivot + offset]! + (index - (pivot + offset)) := by
  rw [Generated.StrictRecombine.resetCombinationSuffix]
  split
  next hposition =>
    let updated := indices.set (pivot + 1 + offset)
      (indices[pivot + offset]'(by omega) + 1)
    have hupdated :
        updated[pivot + offset + 1]! = indices[pivot + offset]! + 1 := by
      rw [getElem!_pos updated (pivot + offset + 1) (by
        dsimp only [updated]
        simp only [Array.size_set]
        omega)]
      simp only [updated, Array.getElem_set]
      rw [if_pos (by omega)]
      rw [getElem!_pos indices (pivot + offset) hpivotOffset]
    by_cases hcurrent : index = pivot + 1 + offset
    · subst index
      rw [resetCombinationSuffix_getElem_before_position updated pivot
        (offset + 1) (pivot + 1 + offset)
        (by
          dsimp only [updated]
          simpa only [Array.size_set] using hposition)
        (by omega)]
      rw [show pivot + 1 + offset = pivot + offset + 1 by omega, hupdated]
      omega
    · rw [resetCombinationSuffix_getElem_suffix updated pivot (offset + 1)
        index (by dsimp only [updated]; simp only [Array.size_set]; omega)
        (by simpa only [updated, Array.size_set] using hindex) (by omega)]
      rw [show pivot + (offset + 1) = pivot + offset + 1 by omega, hupdated]
      omega
  next hposition => omega
termination_by indices.size - (pivot + 1 + offset)
decreasing_by simp only [Array.size_set]; omega

theorem resetCombinationSuffix_getElem_after_pivot (indices : Array Nat)
    (pivot index : Nat) (hpivot : pivot < indices.size)
    (hindex : index < indices.size) (hafter : pivot < index) :
    (Generated.StrictRecombine.resetCombinationSuffix indices pivot 0)[index]! =
      indices[pivot]! + (index - pivot) := by
  exact resetCombinationSuffix_getElem_suffix indices pivot 0 index hpivot
    hindex (by omega)

private def ValidCombination (upper count : Nat) (indices : Array Nat) : Prop :=
  indices.size = count ∧
    ∀ index (hindex : index < indices.size),
      indices[index]! ≤ upper - count + index

private theorem resetCombinationSuffix_preserves_bounds (indices : Array Nat)
    (upper count pivot offset : Nat)
    (hbound : ∀ index (hindex : index < indices.size),
      indices[index]! ≤ upper - count + index) :
    ∀ index (hindex : index <
      (Generated.StrictRecombine.resetCombinationSuffix indices pivot offset).size),
      (Generated.StrictRecombine.resetCombinationSuffix indices pivot offset)[index]! ≤
        upper - count + index := by
  rw [Generated.StrictRecombine.resetCombinationSuffix]
  split
  next hposition =>
    apply resetCombinationSuffix_preserves_bounds
    intro index hindex
    rw [getElem!_pos _ index hindex, Array.getElem_set]
    split
    next heq =>
      subst index
      have hprevious : pivot + offset < indices.size := by omega
      have := hbound (pivot + offset) hprevious
      rw [getElem!_pos indices (pivot + offset) hprevious] at this
      omega
    next hne =>
      have hi : index < indices.size := by simpa using hindex
      have hb := hbound index hi
      rw [getElem!_pos indices index hi] at hb
      exact hb
  next hposition => exact hbound
termination_by indices.size - (pivot + 1 + offset)
decreasing_by simp only [Array.size_set]; omega

theorem nextCombination_size (indices : Array Nat) (upper : Nat) :
    (Generated.StrictRecombine.nextCombination indices upper).2.size =
      indices.size := by
  unfold Generated.StrictRecombine.nextCombination
  split
  next hfits =>
    split
    next hpivot => rfl
    next pivot hpivot =>
      split
      next hpivotBounds =>
        dsimp
        rw [resetCombinationSuffix_size]
        simp
      next hpivotBounds => rfl
  next hfits => rfl

theorem nextCombinationPivot_some_lt (indices : Array Nat)
    (upper inspected pivot : Nat)
    (hresult : Generated.StrictRecombine.nextCombinationPivot indices upper
      inspected = some pivot) :
    pivot < indices.size := by
  induction hmeasure : indices.size - inspected using Nat.strong_induction_on
      generalizing inspected with
  | h measure ih =>
      rw [Generated.StrictRecombine.nextCombinationPivot] at hresult
      split at hresult
      next hinspected =>
        dsimp at hresult
        split at hresult
        next hmaximal =>
          exact ih (indices.size - (inspected + 1)) (by omega)
            (inspected + 1) hresult rfl
        next havailable =>
          have hpivot := Option.some.inj hresult
          subst pivot
          omega
      next hinspected => contradiction

theorem nextCombinationPivot_some_suffix (indices : Array Nat)
    (upper inspected pivot : Nat)
    (hresult : Generated.StrictRecombine.nextCombinationPivot indices upper
    inspected = some pivot) :
    ∀ position (hpivot : pivot < position) (hposition : position < indices.size)
      (hunscanned : inspected ≤ indices.size - 1 - position),
      indices[position] = upper - indices.size + position := by
  induction hmeasure : indices.size - inspected using Nat.strong_induction_on
      generalizing inspected with
  | h measure ih =>
      rw [Generated.StrictRecombine.nextCombinationPivot] at hresult
      split at hresult
      next hinspected =>
        dsimp at hresult
        split at hresult
        next hmaximal =>
          intro position hpivot hposition hunscanned
          by_cases hcurrent : position = indices.size - 1 - inspected
          · subst position
            exact hmaximal
          · exact ih (indices.size - (inspected + 1)) (by omega)
              (inspected + 1) hresult rfl position hpivot hposition (by omega)
        next havailable =>
          intro position hpivot hposition hunscanned
          have hpivotValue := Option.some.inj hresult
          subst pivot
          omega
      next hinspected => contradiction

theorem nextCombinationPivot_some_suffix_zero (indices : Array Nat)
    (upper pivot : Nat)
    (hresult : Generated.StrictRecombine.nextCombinationPivot indices upper 0 =
      some pivot) :
    ∀ position (hpivot : pivot < position) (hposition : position < indices.size),
      indices[position] = upper - indices.size + position := by
  intro position hpivot hposition
  exact nextCombinationPivot_some_suffix indices upper 0 pivot hresult
    position hpivot hposition (by omega)

/-- If the generated right-to-left pivot search finds no incrementable
position, every fixed-size combination digit is at its unique maximal value.
This identifies the concrete terminal state of the C++ enumerator. -/
theorem nextCombinationPivot_none_all_maximal (indices : Array Nat)
    (upper inspected : Nat)
    (hresult : Generated.StrictRecombine.nextCombinationPivot indices upper
      inspected = none) :
    ∀ position (hposition : position < indices.size)
      (hunscanned : inspected ≤ indices.size - 1 - position),
      indices[position] = upper - indices.size + position := by
  induction hmeasure : indices.size - inspected using Nat.strong_induction_on
      generalizing inspected with
  | h measure ih =>
      rw [Generated.StrictRecombine.nextCombinationPivot] at hresult
      split at hresult
      next hinspected =>
        dsimp at hresult
        split at hresult
        next hmaximal =>
          intro position hposition hunscanned
          by_cases hcurrent : position = indices.size - 1 - inspected
          · subst position
            exact hmaximal
          · exact ih (indices.size - (inspected + 1)) (by omega)
              (inspected + 1) hresult rfl position hposition (by omega)
        next havailable => contradiction
      next hinspected =>
        intro position hposition hunscanned
        omega

theorem nextCombinationPivot_none_all_maximal_zero (indices : Array Nat)
    (upper : Nat)
    (hresult : Generated.StrictRecombine.nextCombinationPivot indices upper 0 =
      none) :
    ∀ position (hposition : position < indices.size),
      indices[position] = upper - indices.size + position := by
  intro position hposition
  exact nextCombinationPivot_none_all_maximal indices upper 0 hresult
    position hposition (by omega)

/-- A successful `false` return is exactly the last fixed-size combination;
the generated enumerator cannot stop at an interior array. -/
theorem nextCombination_false_is_final (indices next : Array Nat)
    (upper count : Nat) (hsize : indices.size = count)
    (hfits : count ≤ upper)
    (hrun : Generated.StrictRecombine.nextCombination indices upper =
      (false, next)) :
    next = indices ∧
      ∀ position (hposition : position < indices.size),
        indices[position] = upper - count + position := by
  unfold Generated.StrictRecombine.nextCombination at hrun
  split at hrun
  next harrayFits =>
    split at hrun
    next hpivotNone =>
      have hout := Prod.mk.inj hrun
      constructor
      · exact hout.2.symm
      · intro position hposition
        rw [← hsize]
        exact nextCombinationPivot_none_all_maximal_zero indices upper
          hpivotNone position hposition
    next pivot hpivotSome =>
      split at hrun
      next hpivotBounds => simp at hrun
      next hpivotBounds =>
        exact absurd
          (nextCombinationPivot_some_lt indices upper 0 pivot hpivotSome)
          hpivotBounds
  next harrayFits =>
    exact absurd (hsize.trans_le hfits) harrayFits

theorem nextCombinationPivot_some_ne (indices : Array Nat)
    (upper inspected pivot : Nat)
    (hresult : Generated.StrictRecombine.nextCombinationPivot indices upper
      inspected = some pivot) :
    indices[pivot]'(nextCombinationPivot_some_lt indices upper inspected pivot
      hresult) ≠ upper - indices.size + pivot := by
  induction hmeasure : indices.size - inspected using Nat.strong_induction_on
      generalizing inspected with
  | h measure ih =>
      rw [Generated.StrictRecombine.nextCombinationPivot] at hresult
      split at hresult
      next hinspected =>
        dsimp at hresult
        split at hresult
        next hmaximal =>
          exact ih (indices.size - (inspected + 1)) (by omega)
            (inspected + 1) hresult rfl
        next havailable =>
          have hpivotValue := Option.some.inj hresult
          subst pivot
          exact havailable
      next hinspected => contradiction

theorem nextCombination_true_pivot (indices next : Array Nat) (upper : Nat)
  (hrun : Generated.StrictRecombine.nextCombination indices upper =
      (true, next)) :
    ∃ pivot, ∃ hpivot : pivot < indices.size,
      next[pivot]! = indices[pivot]! + 1 ∧
        ∀ index, index < pivot → next[index]! = indices[index]! := by
  unfold Generated.StrictRecombine.nextCombination at hrun
  split at hrun
  next hfits =>
    split at hrun
    next hpivotNone => simp at hrun
    next pivot hpivotSome =>
      split at hrun
      next hpivotBounds =>
        have hpivotLt := nextCombinationPivot_some_lt indices upper 0 pivot
          hpivotSome
        have hout := Prod.mk.inj hrun
        have hnext := hout.2
        subst next
        refine ⟨pivot, hpivotLt, ?_, ?_⟩
        · rw [resetCombinationSuffix_getElem_le _ pivot 0 pivot
            (by simpa using hpivotLt) (Nat.le_refl _)]
          simp [getElem!_pos, hpivotBounds]
        · intro index hindex
          have hindexBounds : index < indices.size := by omega
          rw [resetCombinationSuffix_getElem_le _ pivot 0 index
            (by simpa using hindexBounds) (Nat.le_of_lt hindex)]
          rw [getElem!_pos _ index (by simpa using hindexBounds),
            getElem!_pos indices index hindexBounds]
          rw [Array.getElem_set]
          rw [if_neg (by omega)]
      next hpivotBounds => simp at hrun
  next hfits => simp at hrun

/-- The complete value-level shape of a successful source successor: prefix
cells are unchanged, the rightmost available pivot is incremented once, and
every suffix cell is reset to the unique consecutive minimum. -/
theorem nextCombination_true_minimal_suffix (indices next : Array Nat)
    (upper : Nat)
    (hrun : Generated.StrictRecombine.nextCombination indices upper =
      (true, next)) :
    ∃ pivot, ∃ hpivot : pivot < indices.size,
      indices[pivot]! ≠ upper - indices.size + pivot ∧
      next[pivot]! = indices[pivot]! + 1 ∧
      (∀ index, index < pivot → next[index]! = indices[index]!) ∧
      (∀ index (hindex : index < indices.size), pivot < index →
        indices[index] = upper - indices.size + index) ∧
      (∀ index (hindex : index < indices.size), pivot < index →
        next[index]! = indices[pivot]! + 1 + (index - pivot)) := by
  unfold Generated.StrictRecombine.nextCombination at hrun
  split at hrun
  next hfits =>
    split at hrun
    next hpivotNone => simp at hrun
    next pivot hpivotSome =>
      split at hrun
      next hpivotBounds =>
        have hpivotLt := nextCombinationPivot_some_lt indices upper 0 pivot
          hpivotSome
        have hout := Prod.mk.inj hrun
        have hnext := hout.2
        subst next
        refine ⟨pivot, hpivotLt, ?_, ?_, ?_, ?_, ?_⟩
        · rw [getElem!_pos indices pivot hpivotLt]
          exact nextCombinationPivot_some_ne indices upper 0 pivot hpivotSome
        · rw [resetCombinationSuffix_getElem_le _ pivot 0 pivot
            (by simpa using hpivotLt) (Nat.le_refl _)]
          simp [getElem!_pos, hpivotBounds]
        · intro index hindex
          have hindexBounds : index < indices.size := by omega
          rw [resetCombinationSuffix_getElem_le _ pivot 0 index
            (by simpa using hindexBounds) (Nat.le_of_lt hindex)]
          rw [getElem!_pos _ index (by simpa using hindexBounds),
            getElem!_pos indices index hindexBounds, Array.getElem_set]
          rw [if_neg (by omega)]
        · intro index hindex hafter
          exact nextCombinationPivot_some_suffix_zero indices upper pivot
            hpivotSome index hafter hindex
        · intro index hindex hafter
          rw [resetCombinationSuffix_getElem_after_pivot _ pivot index
            (by simpa using hpivotLt) (by simpa using hindex) hafter]
          rw [getElem!_pos _ pivot (by simpa using hpivotLt),
            Array.getElem_set, if_pos rfl,
            getElem!_pos indices pivot hpivotLt]
      next hpivotBounds => simp at hrun
  next hfits => simp at hrun

/-- The actual fixed-size subset representation consumed by recombination:
the array has the requested size, its entries are strictly increasing (the
`gap` form also records the accumulated distance), and every entry is below
`upper`. -/
def LegalCombination (upper count : Nat) (indices : Array Nat) : Prop :=
  indices.size = count ∧
  (∀ left right (hleft : left < indices.size) (hright : right < indices.size),
    left < right → indices[left]! + (right - left) ≤ indices[right]!) ∧
  ∀ index (hindex : index < indices.size), indices[index]! < upper

/-- Read the source occurrences named by a concrete list of natural-number
indices.  Bounds are kept as a separate proposition so this definition has
the same total lookup behaviour as the generated array code. -/
def selectSourceIndices {α : Type*} [Inhabited α]
    (source : List α) (indices : List Nat) : List α :=
  indices.map fun index => source[index]!

/-- A legal physical combination selects an occurrence-sensitive sublist of
the source.  Strictly increasing array indices form the finite order
embedding witnessing the sublist relation, so duplicate source values retain
their distinct positions. -/
theorem selectSourceIndices_sublist {α : Type*} [Inhabited α]
    (source : List α) (indices : Array Nat)
    (hlegal : LegalCombination source.length indices.size indices) :
    (selectSourceIndices source indices.toList).Sublist source := by
  rw [List.sublist_iff_exists_fin_orderEmbedding_get_eq]
  let embedding : Fin (selectSourceIndices source indices.toList).length ↪o
      Fin source.length :=
    OrderEmbedding.ofStrictMono
      (fun position =>
        ⟨indices[position.val]!, hlegal.2.2 position.val (by
          simpa [selectSourceIndices] using position.isLt)⟩)
      (by
        intro left right hlt
        have hleft : left.val < indices.size := by
          simpa [selectSourceIndices] using left.isLt
        have hright : right.val < indices.size := by
          simpa [selectSourceIndices] using right.isLt
        have hgap := hlegal.2.1 left.val right.val hleft hright hlt
        have hltVal : left.val < right.val := hlt
        exact Fin.mk_lt_mk.mpr (by omega))
  refine ⟨embedding, ?_⟩
  intro position
  have hposition : position.val < indices.size := by
    simpa [selectSourceIndices] using position.isLt
  simp only [selectSourceIndices, List.length_map, List.get_eq_getElem,
    List.getElem_map]
  rw [Array.getElem_toList hposition]
  have hembedding : (embedding position).val = indices[position.val] := by
    change indices[position.val]! = indices[position.val]
    rw [getElem!_pos indices position.val hposition]
  have hbound : indices[position.val] < source.length := by
    simpa [getElem!_pos indices position.val hposition] using
      hlegal.2.2 position.val hposition
  rw [getElem!_pos source indices[position.val]
    hbound]
  congr 1
  exact hembedding.symm

private theorem getElem!_map_of_lt {α β : Type*} [Inhabited α] [Inhabited β]
    (function : α → β) (values : List α) (position : Nat)
    (hposition : position < values.length) :
    (values.map function)[position]! = function values[position]! := by
  rw [getElem!_pos _ position (by simpa using hposition),
    getElem!_pos values position hposition]
  simp

/-- A source-order sublist has concrete, strictly increasing occurrence
indices.  The accumulated-gap form is exactly the legality invariant used by
the generated C++ combination enumerator; duplicates in `source` therefore
retain their occurrence identity. -/
theorem sublist_exists_legal_source_indices {α : Type*} [Inhabited α]
    {chosen source : List α} (hsublist : chosen.Sublist source) :
    ∃ indices : List Nat,
      indices.length = chosen.length ∧
      (∀ position (hposition : position < indices.length),
        position ≤ indices[position]!) ∧
      (∀ left right (hleft : left < indices.length)
          (hright : right < indices.length), left < right →
        indices[left]! + (right - left) ≤ indices[right]!) ∧
      (∀ position (hposition : position < indices.length),
        indices[position]! < source.length) ∧
      selectSourceIndices source indices = chosen := by
  induction hsublist with
  | slnil => exact ⟨[], rfl, by simp, by simp, by simp, rfl⟩
  | cons element hsublist ih =>
      rcases ih with ⟨indices, hlength, hposition, hgaps, hbounds, hselect⟩
      refine ⟨indices.map Nat.succ, by simp [hlength], ?_, ?_, ?_, ?_⟩
      · intro position hposition'
        simp only [List.length_map] at hposition'
        rw [getElem!_map_of_lt Nat.succ indices position hposition']
        exact Nat.le_trans (hposition position hposition') (Nat.le_succ _)
      · intro left right hleft hright hlt
        simp only [List.length_map] at hleft hright
        rw [getElem!_map_of_lt Nat.succ indices left hleft,
          getElem!_map_of_lt Nat.succ indices right hright]
        have := hgaps left right hleft hright hlt
        omega
      · intro position hposition'
        simp only [List.length_map] at hposition'
        rw [getElem!_map_of_lt Nat.succ indices position hposition']
        have := hbounds position hposition'
        simp only [List.length_cons]
        omega
      · simpa [selectSourceIndices] using hselect
  | cons₂ element hsublist ih =>
      rcases ih with ⟨indices, hlength, hposition, hgaps, hbounds, hselect⟩
      let shifted := indices.map Nat.succ
      refine ⟨0 :: shifted, by simp [shifted, hlength], ?_, ?_, ?_, ?_⟩
      · intro position hposition'
        cases position with
        | zero => simp
        | succ position =>
            simp only [List.length_cons, Nat.succ_lt_succ_iff] at hposition'
            have hposition'' : position < indices.length := by
              simpa [shifted] using hposition'
            simp only [List.getElem!_cons_succ, shifted]
            rw [getElem!_map_of_lt Nat.succ indices position hposition'']
            have := hposition position hposition''
            omega
      · intro left right hleft hright hlt
        cases left with
        | zero =>
            cases right with
            | zero => omega
            | succ right =>
                have hright' : right < indices.length := by
                  simpa [shifted] using hright
                simp only [List.getElem!_cons_zero, List.getElem!_cons_succ,
                  shifted]
                rw [getElem!_map_of_lt Nat.succ indices right hright']
                have := hposition right hright'
                omega
        | succ left =>
            cases right with
            | zero => omega
            | succ right =>
                simp only [List.length_cons, Nat.succ_lt_succ_iff] at hleft hright
                have hleft' : left < indices.length := by
                  simpa [shifted] using hleft
                have hright' : right < indices.length := by
                  simpa [shifted] using hright
                simp only [List.getElem!_cons_succ, shifted]
                rw [getElem!_map_of_lt Nat.succ indices left hleft',
                  getElem!_map_of_lt Nat.succ indices right hright']
                have := hgaps left right hleft' hright' (by omega)
                omega
      · intro position hposition'
        cases position with
        | zero => simp
        | succ position =>
            simp only [List.length_cons, Nat.succ_lt_succ_iff] at hposition'
            have hposition'' : position < indices.length := by
              simpa [shifted] using hposition'
            simp only [List.getElem!_cons_succ, shifted, List.length_cons]
            rw [getElem!_map_of_lt Nat.succ indices position hposition'']
            have := hbounds position hposition''
            omega
      · simpa [selectSourceIndices, shifted] using congrArg (List.cons element) hselect

/-- Array form consumed directly by `scanZassenhausCombinations`: every
source sublist supplies a legal candidate of exactly the same cardinality,
and selecting its array entries recovers that very sublist. -/
theorem sublist_exists_legal_combination {α : Type*} [Inhabited α]
    {chosen source : List α} (hsublist : chosen.Sublist source) :
    ∃ indices : Array Nat,
      LegalCombination source.length chosen.length indices ∧
      selectSourceIndices source indices.toList = chosen := by
  rcases sublist_exists_legal_source_indices hsublist with
    ⟨indices, hlength, _hposition, hgaps, hbounds, hselect⟩
  refine ⟨indices.toArray, ?_, by simpa using hselect⟩
  refine ⟨by simpa using hlength, ?_, ?_⟩
  · intro left right hleft hright hlt
    have hleft' : left < indices.length := by simpa using hleft
    have hright' : right < indices.length := by simpa using hright
    rw [getElem!_pos (indices.toArray) left hleft,
      getElem!_pos (indices.toArray) right hright]
    simpa [getElem!_pos indices left hleft',
      getElem!_pos indices right hright'] using
      hgaps left right hleft' hright' hlt
  · intro position hposition'
    have hposition : position < indices.length := by simpa using hposition'
    rw [getElem!_pos (indices.toArray) position hposition']
    simpa [getElem!_pos indices position hposition] using
      hbounds position hposition

/-- Exact successful output of the generated checked `Nat` to `Int32`
candidate lowering loop.  The theorem follows the actual source recursion
and records every emitted element in order. -/
theorem combinationToInt32Loop_toList (indices : Array Nat) (index : Nat)
    (result : Array Int32)
    (hfits : ∀ position (hposition : position < indices.size),
      indices[position] < 2 ^ 31) :
    ∃ output,
      Generated.StrictRecombine.combinationToInt32Loop indices index result =
        .ok output ∧
      output.toList = result.toList ++
        (indices.toList.drop index).map
          (fun value => value.toUInt32.toInt32) := by
  induction hmeasure : indices.size - index using Nat.strong_induction_on
      generalizing index result with
  | h measure ih =>
      rw [Generated.StrictRecombine.combinationToInt32Loop]
      split
      next hindex =>
        rw [dif_pos (hfits index hindex)]
        rcases ih (indices.size - (index + 1)) (by omega) (index + 1)
            (result.push indices[index].toUInt32.toInt32) rfl with
          ⟨output, hrun, houtput⟩
        refine ⟨output, hrun, ?_⟩
        rw [houtput, Array.toList_push, List.append_assoc]
        have hdrop : indices.toList.drop index = indices[index] ::
            indices.toList.drop (index + 1) := by
          simpa using List.drop_eq_getElem_cons
            (l := indices.toList) (i := index) (by simpa using hindex)
        rw [hdrop]
        rfl
      next hindex =>
        refine ⟨result, rfl, ?_⟩
        have hdrop : indices.toList.drop index = [] :=
          List.drop_eq_nil_iff.mpr (by simpa using Nat.not_lt.mp hindex)
        simp [hdrop]

/-- Exact top-level checked conversion used by `zassenhausAttempt`. -/
theorem combinationToInt32_toList (indices : Array Nat)
    (hfits : ∀ position (hposition : position < indices.size),
      indices[position] < 2 ^ 31) :
    ∃ output,
      Generated.StrictRecombine.combinationToInt32 indices = .ok output ∧
      output.toList = indices.toList.map
        (fun value => value.toUInt32.toInt32) := by
  simpa [Generated.StrictRecombine.combinationToInt32] using
    combinationToInt32Loop_toList indices 0 #[] hfits

theorem nat_toUInt32_toInt32_nonnegative_and_toNat (value : Nat)
    (hfits : value < 2 ^ 31) :
    0 ≤ value.toUInt32.toInt32 ∧
      value.toUInt32.toInt32.toInt64.toNat = value := by
  have htoInt : value.toUInt32.toInt32.toInt = (value : Int) := by
    change value.toInt32.toInt = (value : Int)
    exact Int32.toInt_ofNat_of_lt hfits
  have hnonnegative : 0 ≤ value.toUInt32.toInt32 := by
    rw [Int32.le_iff_toInt_le, Int32.toInt_zero, htoInt]
    omega
  refine ⟨hnonnegative, ?_⟩
  change value.toUInt32.toInt32.toInt64.toNatClampNeg = value
  rw [Int32.toNatClampNeg_toInt64]
  unfold Int32.toNatClampNeg
  rw [htoInt]
  simp

def selectedLeadingValues (candidate : Array Nat)
    (activeLifted : Array SparsePolyZZ) (index : Nat) : List ZZ :=
  (candidate.toList.drop index).map fun activeIndex =>
    (activeLifted[activeIndex]!)[0]!.2

def selectedConstantValues (candidate : Array Nat)
    (activeLifted : Array SparsePolyZZ) (index : Nat) : List ZZ :=
  (candidate.toList.drop index).map fun activeIndex =>
    Generated.StrictRecombine.constantTerm activeLifted[activeIndex]!

/-- Exact execution and value of the generated leading-coefficient pruning
loop for a bounded candidate of nonempty active factors. -/
theorem selectedLeadingProductLoop_succeeds
    (candidate : Array Nat) (activeLifted : Array SparsePolyZZ)
    (index : Nat) (acc : ZZ)
    (hbound : ∀ position (hposition : position < candidate.size),
      candidate[position] < activeLifted.size)
    (hnonempty : ∀ position (hposition : position < candidate.size),
      activeLifted[candidate[position]]!.isEmpty = false) :
    Generated.StrictRecombine.selectedLeadingProductLoop candidate
      activeLifted index acc =
        .ok (acc * (selectedLeadingValues candidate activeLifted index).prod) := by
  induction hmeasure : candidate.size - index using Nat.strong_induction_on
      generalizing index acc with
  | h measure ih =>
      rw [Generated.StrictRecombine.selectedLeadingProductLoop]
      split
      next hindex =>
        have hactive := hbound index hindex
        rw [dif_pos hactive]
        have hfactor : activeLifted[candidate[index]].isEmpty = false := by
          rw [← getElem!_pos activeLifted candidate[index] hactive]
          exact hnonempty index hindex
        have hfactorFalse : ¬ activeLifted[candidate[index]].isEmpty := by
          simp [hfactor]
        have hnonemptySize : 0 < activeLifted[candidate[index]].size := by
          have hsizeNe : activeLifted[candidate[index]].size ≠ 0 := by
            intro hsize
            apply hfactorFalse
            simp [Array.isEmpty, hsize]
          omega
        rw [dif_neg hfactorFalse]
        rw [ih (candidate.size - (index + 1)) (by omega) (index + 1)
          (acc * ((activeLifted[candidate[index]]'hactive)[0]'hnonemptySize).2) rfl]
        have hdrop : candidate.toList.drop index = candidate[index] ::
            candidate.toList.drop (index + 1) := by
          simpa using List.drop_eq_getElem_cons
            (l := candidate.toList) (i := index) (by simpa using hindex)
        simp only [selectedLeadingValues, hdrop, List.map_cons, List.prod_cons]
        rw [getElem!_pos activeLifted candidate[index] hactive,
          getElem!_pos activeLifted[candidate[index]] 0 hnonemptySize]
        ring
      next hindex =>
        have hdrop : candidate.toList.drop index = [] :=
          List.drop_eq_nil_iff.mpr (by simpa using Nat.le_of_not_gt hindex)
        simp [selectedLeadingValues, hdrop]

/-- Exact execution and value of the generated constant-coefficient pruning
loop for every bounded candidate. -/
theorem selectedConstantProductLoop_succeeds
    (candidate : Array Nat) (activeLifted : Array SparsePolyZZ)
    (index : Nat) (acc : ZZ)
    (hbound : ∀ position (hposition : position < candidate.size),
      candidate[position] < activeLifted.size) :
    Generated.StrictRecombine.selectedConstantProductLoop candidate
      activeLifted index acc =
        .ok (acc * (selectedConstantValues candidate activeLifted index).prod) := by
  induction hmeasure : candidate.size - index using Nat.strong_induction_on
      generalizing index acc with
  | h measure ih =>
      rw [Generated.StrictRecombine.selectedConstantProductLoop]
      split
      next hindex =>
        have hactive := hbound index hindex
        rw [dif_pos hactive]
        rw [ih (candidate.size - (index + 1)) (by omega) (index + 1)
          (acc * Generated.StrictRecombine.constantTerm
            (activeLifted[candidate[index]]'hactive)) rfl]
        have hdrop : candidate.toList.drop index = candidate[index] ::
            candidate.toList.drop (index + 1) := by
          simpa using List.drop_eq_getElem_cons
            (l := candidate.toList) (i := index) (by simpa using hindex)
        simp only [selectedConstantValues, hdrop, List.map_cons, List.prod_cons]
        rw [getElem!_pos activeLifted candidate[index] hactive]
        ring
      next hindex =>
        have hdrop : candidate.toList.drop index = [] :=
          List.drop_eq_nil_iff.mpr (by simpa using Nat.le_of_not_gt hindex)
        simp [selectedConstantValues, hdrop]

/-- The exact arithmetic test used by both generated Zassenhaus pruning
branches cannot reject an integer that divides its target. -/
theorem zassenhaus_prune_condition_false_of_dvd
    (target recovered : ZZ) (hdivides : recovered ∣ target) :
    ¬(recovered ≠ 0 ∧ ZZ.fdiv_r 0 target recovered ≠ 0) := by
  intro hreject
  exact hreject.2 (by
    unfold ZZ.fdiv_r
    exact Int.fmod_eq_zero_of_dvd hdivides)

/-- Divisibility of integer polynomials reaches both boundary
coefficients.  These are precisely the two coefficients inspected by the
generated pruning code. -/
theorem polynomial_divisor_boundary_coefficients
    {divisor dividend : Polynomial Int} (hdivides : divisor ∣ dividend) :
    divisor.leadingCoeff ∣ dividend.leadingCoeff ∧
      divisor.coeff 0 ∣ dividend.coeff 0 := by
  constructor
  · exact Polynomial.leadingCoeff_dvd_leadingCoeff hdivides
  · rcases hdivides with ⟨quotient, rfl⟩
    refine ⟨quotient.coeff 0, ?_⟩
    simp

theorem zassenhaus_leading_prune_accepts_associated_divisor
    {divisor dividend : Polynomial Int} (recovered : ZZ)
    (hdivides : divisor ∣ dividend)
    (hassociated : Associated recovered divisor.leadingCoeff) :
    ¬(recovered ≠ 0 ∧
      ZZ.fdiv_r 0 (dividend.leadingCoeff * dividend.leadingCoeff)
        recovered ≠ 0) := by
  apply zassenhaus_prune_condition_false_of_dvd
  have hboundary := (polynomial_divisor_boundary_coefficients hdivides).1
  exact dvd_mul_of_dvd_left (dvd_trans hassociated.dvd hboundary) _

theorem zassenhaus_constant_prune_accepts_associated_divisor
    {divisor dividend : Polynomial Int} (leading recovered : ZZ)
    (hdivides : divisor ∣ dividend)
    (hassociated : Associated recovered (divisor.coeff 0)) :
    ¬(recovered ≠ 0 ∧
      ZZ.fdiv_r 0 (leading * dividend.coeff 0) recovered ≠ 0) := by
  apply zassenhaus_prune_condition_false_of_dvd
  have hboundary := (polynomial_divisor_boundary_coefficients hdivides).2
  exact dvd_mul_of_dvd_right (dvd_trans hassociated.dvd hboundary) _

theorem leadingCoeff_list_prod (factors : List (Polynomial Int)) :
    factors.prod.leadingCoeff = (factors.map Polynomial.leadingCoeff).prod := by
  induction factors with
  | nil => simp
  | cons factor factors ih => simp [ih]

theorem coeff_zero_list_prod (factors : List (Polynomial Int)) :
    factors.prod.coeff 0 = (factors.map fun factor => factor.coeff 0).prod := by
  induction factors with
  | nil => simp
  | cons factor factors ih => simp [ih]

private theorem sparseZZTail_coeff_zero_above
    (limit : Nat) (terms : List (UMonomial × Int))
    (hbelow : ∀ term ∈ terms, term.1.deg < limit) :
    ((terms.map fun term =>
      Polynomial.monomial term.1.deg term.2).sum).coeff limit = 0 := by
  induction terms with
  | nil => simp
  | cons term terms ih =>
      have hterm := hbelow term (by simp)
      have htail : ∀ item ∈ terms, item.1.deg < limit := by
        intro item hitem
        exact hbelow item (by simp [hitem])
      simp [Polynomial.coeff_monomial, ne_of_lt hterm, ih htail]

private theorem sparseZZ_chain_head_gt_all
    (head : UMonomial × Int) (rest : List (UMonomial × Int))
    (hchain : List.IsChain
      (fun a b : UMonomial × Int => a.1.deg > b.1.deg) (head :: rest)) :
    ∀ item ∈ rest, item.1.deg < head.1.deg := by
  induction rest generalizing head with
  | nil => simp
  | cons next rest ih =>
      rw [List.isChain_cons_cons] at hchain
      intro item hitem
      rcases List.mem_cons.mp hitem with rfl | htail
      · exact hchain.1
      · exact Nat.lt_trans (ih next hchain.2 item htail) hchain.1

/-- A nonempty canonical sparse integer polynomial has mathematical degree
equal to the degree stored in its physical first cell. -/
theorem sparsePolyZZ_toPoly_degree_eq_head (poly : SparsePolyZZ)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical poly)
    (hnonempty : 0 < poly.size) :
    (SparsePolyZZ.toPoly poly).degree = poly[0]!.1.deg := by
  have hlistNonempty : poly.toList ≠ [] := by
    intro hempty
    have hlength := congrArg List.length hempty
    have hsizeZero : poly.size = 0 := by simpa using hlength
    omega
  obtain ⟨head, rest, hlist⟩ := List.exists_cons_of_ne_nil hlistNonempty
  have hheadEq : head = poly[0] := by
    have hget := Array.getElem_toList hnonempty
    simpa [hlist] using hget
  have hheadMem : head ∈ poly.toList := by simp [hlist]
  have hheadNonzero : head.2 ≠ 0 := hcanonical.2 head hheadMem
  have hchain : List.IsChain
      (fun a b : UMonomial × Int => a.1.deg > b.1.deg)
      (head :: rest) := by
    simpa [hlist] using hcanonical.1
  have hrestLt : ∀ item ∈ rest, item.1.deg < head.1.deg :=
    sparseZZ_chain_head_gt_all head rest hchain
  have hrestDegree :
      ((rest.map fun term =>
        Polynomial.monomial term.1.deg term.2).sum).degree < head.1.deg := by
    rw [Polynomial.degree_lt_iff_coeff_zero]
    intro degree hdegree
    apply sparseZZTail_coeff_zero_above
    intro item hitem
    exact Nat.lt_of_lt_of_le (hrestLt item hitem) hdegree
  have hheadDegree :
      (Polynomial.monomial head.1.deg head.2).degree = head.1.deg :=
    Polynomial.degree_monomial _ hheadNonzero
  unfold SparsePolyZZ.toPoly
  rw [hlist, List.map_cons, List.sum_cons,
    Polynomial.degree_add_eq_left_of_degree_lt (by
      simpa [hheadDegree] using hrestDegree), hheadDegree]
  simp [hheadEq, getElem!_pos poly 0 hnonempty]

/-- The first stored coefficient of a nonempty canonical sparse integer
polynomial is exactly its mathematical leading coefficient. -/
theorem sparsePolyZZ_leadingCoeff_eq_head (poly : SparsePolyZZ)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical poly)
    (hnonempty : 0 < poly.size) :
    (SparsePolyZZ.toPoly poly).leadingCoeff = poly[0].2 := by
  have hlistNonempty : poly.toList ≠ [] := by
    intro hempty
    have hlength := congrArg List.length hempty
    have hsizeZero : poly.size = 0 := by simpa using hlength
    omega
  obtain ⟨head, rest, hlist⟩ := List.exists_cons_of_ne_nil hlistNonempty
  have hheadEq : head = poly[0] := by
    have hget := Array.getElem_toList hnonempty
    simpa [hlist] using hget
  have hheadMem : head ∈ poly.toList := by simp [hlist]
  have hheadNonzero : head.2 ≠ 0 := hcanonical.2 head hheadMem
  have hchain : List.IsChain
      (fun a b : UMonomial × Int => a.1.deg > b.1.deg)
      (head :: rest) := by
    simpa [hlist] using hcanonical.1
  have hrestLt : ∀ item ∈ rest, item.1.deg < head.1.deg :=
    sparseZZ_chain_head_gt_all head rest hchain
  have hrestDegree :
      ((rest.map fun term =>
        Polynomial.monomial term.1.deg term.2).sum).degree < head.1.deg := by
    rw [Polynomial.degree_lt_iff_coeff_zero]
    intro degree hdegree
    apply sparseZZTail_coeff_zero_above
    intro item hitem
    exact Nat.lt_of_lt_of_le (hrestLt item hitem) hdegree
  have hheadDegree :
      (Polynomial.monomial head.1.deg head.2).degree = head.1.deg :=
    Polynomial.degree_monomial _ hheadNonzero
  have hdegree : (SparsePolyZZ.toPoly poly).degree = head.1.deg := by
    unfold SparsePolyZZ.toPoly
    rw [hlist, List.map_cons, List.sum_cons,
      Polynomial.degree_add_eq_left_of_degree_lt (by
        simpa [hheadDegree] using hrestDegree), hheadDegree]
  have hnatDegree : (SparsePolyZZ.toPoly poly).natDegree = head.1.deg :=
    Polynomial.natDegree_eq_of_degree_eq_some hdegree
  rw [Polynomial.leadingCoeff, hnatDegree]
  unfold SparsePolyZZ.toPoly
  rw [hlist, List.map_cons, List.sum_cons, Polynomial.coeff_add,
    Polynomial.coeff_monomial, if_pos rfl,
    sparseZZTail_coeff_zero_above head.1.deg rest hrestLt, add_zero,
    hheadEq]

private def sparseZZLastConstant : List (UMonomial × Int) → Int
  | [] => 0
  | term :: [] => if term.1.deg = 0 then term.2 else 0
  | _ :: next :: rest => sparseZZLastConstant (next :: rest)

private theorem sparseZZLastConstant_eq_last
    (head : UMonomial × Int) (rest : List (UMonomial × Int)) :
    sparseZZLastConstant (head :: rest) =
      let last := (head :: rest)[(head :: rest).length - 1]!
      if last.1.deg = 0 then last.2 else 0 := by
  induction rest generalizing head with
  | nil => simp [sparseZZLastConstant]
  | cons next rest ih =>
      rw [sparseZZLastConstant]
      simpa using ih next

private theorem sparseZZ_coeff_zero_eq_lastConstant
    (terms : List (UMonomial × Int))
    (hchain : List.IsChain
      (fun a b : UMonomial × Int => a.1.deg > b.1.deg) terms) :
    ((terms.map fun term =>
      Polynomial.monomial term.1.deg term.2).sum).coeff 0 =
        sparseZZLastConstant terms := by
  induction terms with
  | nil => simp [sparseZZLastConstant]
  | cons head terms ih =>
      cases terms with
      | nil => simp [sparseZZLastConstant, Polynomial.coeff_monomial]
      | cons next rest =>
          rw [List.isChain_cons_cons] at hchain
          have hheadDegree : head.1.deg ≠ 0 := by omega
          simp only [List.map_cons, List.sum_cons, Polynomial.coeff_add,
            Polynomial.coeff_monomial, if_neg hheadDegree, zero_add,
            sparseZZLastConstant]
          simpa [Polynomial.coeff_monomial] using ih hchain.2

/-- The generated sparse-array last-term lookup computes exactly the
mathematical constant coefficient on every canonical integer polynomial. -/
theorem sparsePolyZZ_constantTerm_eq_coeff_zero (poly : SparsePolyZZ)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical poly) :
    Generated.StrictRecombine.constantTerm poly =
      (SparsePolyZZ.toPoly poly).coeff 0 := by
  by_cases hempty : poly.isEmpty
  · have hpoly : poly = #[] := Array.isEmpty_iff.mp hempty
    subst poly
    simp [Generated.StrictRecombine.constantTerm, SparsePolyZZ.toPoly]
  · have hsize : 0 < poly.size := by
      have hsizeNe : poly.size ≠ 0 := by
        intro hzero
        apply hempty
        simp [Array.isEmpty, hzero]
      omega
    have hlastModel := sparseZZLastConstant_eq_last poly[0]!
      poly.toList.tail
    have hlistNonempty : poly.toList ≠ [] := by
      intro hnil
      have hlengthZero : poly.toList.length = 0 := by
        simpa using congrArg List.length hnil
      have hsizeZero : poly.size = 0 := by simpa using hlengthZero
      omega
    have hhead : poly.toList.head hlistNonempty = poly[0]! := by
      rw [List.head_eq_getElem]
      simpa [getElem!_pos poly 0 hsize] using Array.getElem_toList hsize
    have hfull : poly[0]! :: poly.toList.tail = poly.toList := by
      rw [← hhead]
      exact List.cons_head_tail hlistNonempty
    rw [hfull] at hlastModel
    have hcoeff := sparseZZ_coeff_zero_eq_lastConstant poly.toList hcanonical.1
    unfold Generated.StrictRecombine.constantTerm
    rw [dif_neg hempty]
    symm
    unfold SparsePolyZZ.toPoly
    rw [hcoeff, hlastModel]
    congr 1
    have hlast : poly.size - 1 < poly.size := by omega
    rw [getElem!_pos poly.toList (poly.toList.length - 1) (by simpa using hlast)]
    exact congrArg (fun term : UMonomial × Int =>
      if term.1.deg = 0 then term.2 else 0) (Array.getElem_toList hlast)

/-- Once the canonical sparse representation identifies each selected head
coefficient, the concrete leading-pruning product is the leading coefficient
of the exact source-sublist polynomial product. -/
theorem selectedLeadingValues_prod_eq_leadingCoeff
    (candidate : Array Nat) (activeLifted : Array SparsePolyZZ)
    (hbound : ∀ position (hposition : position < candidate.size),
      candidate[position] < activeLifted.size)
    (hheads : ∀ position (hposition : position < candidate.size),
      (SparsePolyZZ.toPoly activeLifted[candidate[position]]!).leadingCoeff =
        ((activeLifted[candidate[position]]!)[0]!).2) :
    (selectedLeadingValues candidate activeLifted 0).prod =
      ((selectSourceIndices activeLifted.toList candidate.toList).map
        SparsePolyZZ.toPoly).prod.leadingCoeff := by
  rw [leadingCoeff_list_prod]
  unfold selectedLeadingValues selectSourceIndices
  simp only [List.drop_zero, List.map_map]
  congr 1
  apply List.map_congr_left
  intro activeIndex hactiveIndex
  rcases List.mem_iff_getElem.mp hactiveIndex with
    ⟨position, hposition, rfl⟩
  have hpositionArray : position < candidate.size := by simpa using hposition
  have hactive := hbound position hpositionArray
  simp only [Function.comp_apply]
  simp only [Array.getElem_toList]
  rw [getElem!_pos activeLifted.toList candidate[position] (by simpa using hactive),
    Array.getElem_toList hactive]
  simpa [getElem!_pos activeLifted candidate[position] hactive] using
    (hheads position hpositionArray).symm

/-- Constant-coefficient analogue of
`selectedLeadingValues_prod_eq_leadingCoeff`. -/
theorem selectedConstantValues_prod_eq_coeff_zero
    (candidate : Array Nat) (activeLifted : Array SparsePolyZZ)
    (hbound : ∀ position (hposition : position < candidate.size),
      candidate[position] < activeLifted.size)
    (hconstants : ∀ position (hposition : position < candidate.size),
      (SparsePolyZZ.toPoly activeLifted[candidate[position]]!).coeff 0 =
        Generated.StrictRecombine.constantTerm
          activeLifted[candidate[position]]!) :
    (selectedConstantValues candidate activeLifted 0).prod =
      ((selectSourceIndices activeLifted.toList candidate.toList).map
        SparsePolyZZ.toPoly).prod.coeff 0 := by
  rw [coeff_zero_list_prod]
  unfold selectedConstantValues selectSourceIndices
  simp only [List.drop_zero, List.map_map]
  congr 1
  apply List.map_congr_left
  intro activeIndex hactiveIndex
  rcases List.mem_iff_getElem.mp hactiveIndex with
    ⟨position, hposition, rfl⟩
  have hpositionArray : position < candidate.size := by simpa using hposition
  have hactive := hbound position hpositionArray
  simp only [Function.comp_apply]
  simp only [Array.getElem_toList]
  rw [getElem!_pos activeLifted.toList candidate[position] (by simpa using hactive),
    Array.getElem_toList hactive]
  simpa [getElem!_pos activeLifted candidate[position] hactive] using
    (hconstants position hpositionArray).symm

theorem selectedLeadingValues_prod_eq_leadingCoeff_of_canonical
    (candidate : Array Nat) (activeLifted : Array SparsePolyZZ)
    (hbound : ∀ position (hposition : position < candidate.size),
      candidate[position] < activeLifted.size)
    (hcanonical : ∀ position (hposition : position < candidate.size),
      StrictPolynomialMod.SparsePolyZZCanonical
        activeLifted[candidate[position]]!)
    (hnonempty : ∀ position (hposition : position < candidate.size),
      activeLifted[candidate[position]]!.isEmpty = false) :
    (selectedLeadingValues candidate activeLifted 0).prod =
      ((selectSourceIndices activeLifted.toList candidate.toList).map
        SparsePolyZZ.toPoly).prod.leadingCoeff := by
  apply selectedLeadingValues_prod_eq_leadingCoeff candidate activeLifted hbound
  intro position hposition
  have hsize : 0 < activeLifted[candidate[position]]!.size := by
    have hsizeNe : activeLifted[candidate[position]]!.size ≠ 0 := by
      intro hzero
      have hempty : activeLifted[candidate[position]]!.isEmpty = true := by
        simp [Array.isEmpty, hzero]
      rw [hnonempty position hposition] at hempty
      contradiction
    omega
  simpa [getElem!_pos activeLifted[candidate[position]]! 0 hsize] using
    sparsePolyZZ_leadingCoeff_eq_head activeLifted[candidate[position]]!
      (hcanonical position hposition) hsize

theorem selectedConstantValues_prod_eq_coeff_zero_of_canonical
    (candidate : Array Nat) (activeLifted : Array SparsePolyZZ)
    (hbound : ∀ position (hposition : position < candidate.size),
      candidate[position] < activeLifted.size)
    (hcanonical : ∀ position (hposition : position < candidate.size),
      StrictPolynomialMod.SparsePolyZZCanonical
        activeLifted[candidate[position]]!) :
    (selectedConstantValues candidate activeLifted 0).prod =
      ((selectSourceIndices activeLifted.toList candidate.toList).map
        SparsePolyZZ.toPoly).prod.coeff 0 := by
  apply selectedConstantValues_prod_eq_coeff_zero candidate activeLifted hbound
  intro position hposition
  exact (sparsePolyZZ_constantTerm_eq_coeff_zero
    activeLifted[candidate[position]]! (hcanonical position hposition)).symm

/-- Lexicographic order on equal-size source arrays, stated at the concrete
first differing position.  This avoids assigning an order to arrays of
different lengths, which never occur in the C++ combination scan. -/
private def ArrayLexLT (left right : Array Nat) : Prop :=
  left.size = right.size ∧ ∃ pivot, pivot < left.size ∧
    (∀ index, index < pivot → left[index]! = right[index]!) ∧
    left[pivot]! < right[pivot]!

private theorem arrayLex_trichotomy (left right : Array Nat)
    (hsize : left.size = right.size) :
    left = right ∨ ArrayLexLT left right ∨ ArrayLexLT right left := by
  by_cases heq : left = right
  · exact Or.inl heq
  · have hexists : ∃ index, index < left.size ∧
        left[index]! ≠ right[index]! := by
      by_contra hall
      push Not at hall
      apply heq
      apply Array.ext hsize
      intro index hleft hright
      have heqBang := hall index hleft
      rw [getElem!_pos left index hleft,
        getElem!_pos right index (by omega)] at heqBang
      exact heqBang
    let pivot := Nat.find hexists
    have hpivotSpec := Nat.find_spec hexists
    have hprefix : ∀ index, index < pivot →
        left[index]! = right[index]! := by
      intro index hindex
      by_contra hne
      exact (Nat.find_min hexists hindex) ⟨by omega, hne⟩
    rcases lt_or_gt_of_ne hpivotSpec.2 with hpivotValue | hpivotValue
    · exact Or.inr (Or.inl
        ⟨hsize, pivot, hpivotSpec.1, hprefix, hpivotValue⟩)
    · exact Or.inr (Or.inr
        ⟨hsize.symm, pivot, by omega, fun index hindex =>
          (hprefix index hindex).symm, hpivotValue⟩)

private theorem legalCombination_positional_bound (upper count : Nat)
    (indices : Array Nat) (hlegal : LegalCombination upper count indices) :
    ∀ index (hindex : index < indices.size),
      indices[index]! ≤ upper - count + index := by
  intro index hindex
  rcases hlegal with ⟨hsize, hgap, hupper⟩
  have hlast : indices.size - 1 < indices.size := by omega
  have hlastUpper := hupper (indices.size - 1) hlast
  by_cases hindexLast : index = indices.size - 1
  · subst index
    omega
  · have htoLast := hgap index (indices.size - 1) hindex hlast (by omega)
    omega

private theorem legalCombination_valid (upper count : Nat)
    (indices : Array Nat) (hlegal : LegalCombination upper count indices) :
    ValidCombination upper count indices :=
  ⟨hlegal.1, legalCombination_positional_bound upper count indices hlegal⟩

private theorem legalCombination_lower_bound (upper count : Nat)
    (indices : Array Nat) (hlegal : LegalCombination upper count indices) :
    ∀ index (hindex : index < indices.size), index ≤ indices[index]! := by
  intro index hindex
  rcases hlegal with ⟨hsize, hgap, hupper⟩
  by_cases hzero : index = 0
  · subst index
    omega
  · have hzeroBounds : 0 < indices.size := by omega
    have hfromZero := hgap 0 index hzeroBounds hindex (by omega)
    omega

private theorem legalCombination_count_le (upper count : Nat)
    (indices : Array Nat) (hlegal : LegalCombination upper count indices) :
    count ≤ upper := by
  by_cases hzero : count = 0
  · omega
  · have hlast : count - 1 < indices.size := by rw [hlegal.1]; omega
    have hlower := legalCombination_lower_bound upper count indices hlegal
      (count - 1) hlast
    have hupper := hlegal.2.2 (count - 1) hlast
    omega

/-- Public cardinality bound carried by every concrete legal combination. -/
theorem LegalCombination.count_le {upper count : Nat} {indices : Array Nat}
    (hlegal : LegalCombination upper count indices) : count ≤ upper :=
  legalCombination_count_le upper count indices hlegal

private theorem initialCombination_lex_le (upper count : Nat)
    (target : Array Nat) (hlegal : LegalCombination upper count target) :
    Generated.StrictRecombine.initialCombination count = target ∨
      ArrayLexLT (Generated.StrictRecombine.initialCombination count) target := by
  have hsize : (Generated.StrictRecombine.initialCombination count).size =
      target.size := by rw [initialCombination_size, hlegal.1]
  rcases arrayLex_trichotomy
      (Generated.StrictRecombine.initialCombination count) target hsize with
    heq | hforward | hbackward
  · exact Or.inl heq
  · exact Or.inr hforward
  · exfalso
    rcases hbackward with ⟨_, pivot, hpivot, hprefix, hpivotValue⟩
    have hpivotTarget : pivot < target.size := hpivot
    have hpivotInitial :
        pivot < (Generated.StrictRecombine.initialCombination count).size := by
      omega
    have hlower := legalCombination_lower_bound upper count target hlegal
      pivot hpivotTarget
    have hpivotCount : pivot < count := by
      simpa [initialCombination_size] using hpivotInitial
    have hinitialValue :
        (Generated.StrictRecombine.initialCombination count)[pivot]! = pivot := by
      have harrayList :
          (Generated.StrictRecombine.initialCombination count)[pivot]! =
            (Generated.StrictRecombine.initialCombination count).toList[pivot]! := by
        rw [getElem!_pos _ pivot hpivotInitial,
          getElem!_pos _ pivot (by simpa using hpivotInitial)]
        exact (Array.getElem_toList hpivotInitial).symm
      rw [harrayList, initialCombination_toList]
      simp [getElem!_pos, hpivotCount]
    rw [hinitialValue] at hpivotValue
    omega

private theorem finalCombination_not_lt_legal (indices target : Array Nat)
    (upper count : Nat) (hlegal : LegalCombination upper count target)
    (hsize : indices.size = count)
    (hfinal : ∀ position (hposition : position < indices.size),
      indices[position] = upper - count + position) :
    ¬ ArrayLexLT indices target := by
  rintro ⟨hsizes, pivot, hpivot, hprefix, hpivotValue⟩
  have hpivotTarget : pivot < target.size := by omega
  have hbound := legalCombination_positional_bound upper count target hlegal
    pivot hpivotTarget
  rw [getElem!_pos indices pivot hpivot, hfinal pivot hpivot] at hpivotValue
  omega

/-- A successful generated successor is immediate among all legal
fixed-size combinations: no legal array lies strictly between the input and
output.  The proof uses the real rightmost-pivot state, including maximal old
suffix digits, and the real reset loop's exact minimal new suffix. -/
theorem nextCombination_true_no_legal_between (indices next middle : Array Nat)
    (upper count : Nat) (hlegal : LegalCombination upper count middle)
    (hrun : Generated.StrictRecombine.nextCombination indices upper =
      (true, next))
    (hleft : ArrayLexLT indices middle) (hright : ArrayLexLT middle next) :
    False := by
  obtain ⟨pivot, hpivot, hpivotNotMax, hpivotNext, hprefixNext, hmaxSuffix,
    hminSuffix⟩ :=
    nextCombination_true_minimal_suffix indices next upper hrun
  rcases hleft with ⟨hsizeMiddle, leftPivot, hleftPivot,
    hprefixMiddle, hleftValue⟩
  rcases hright with ⟨hsizeNext, rightPivot, hrightPivot,
    hmiddleNextPrefix, hrightValue⟩
  rcases hlegal with ⟨hmiddleCount, hmiddleGap, hmiddleUpper⟩
  have hpivotMiddle : pivot < middle.size := by omega
  have hpivotNextBounds : pivot < next.size := by omega
  by_cases hleftBefore : leftPivot < pivot
  · have hmiddleAtLeft : middle[leftPivot]! > next[leftPivot]! := by
      rw [hprefixNext leftPivot hleftBefore]
      exact hleftValue
    by_cases hrightBefore : rightPivot < leftPivot
    · have hindicesMiddleAtRight := hprefixMiddle rightPivot hrightBefore
      have hnextIndicesAtRight := hprefixNext rightPivot (by omega)
      omega
    · by_cases hsame : rightPivot = leftPivot
      · subst rightPivot
        omega
      · have hmiddleNextAtLeft := hmiddleNextPrefix leftPivot (by omega)
        omega
  · by_cases hpivotBefore : pivot < leftPivot
    · have hmaxCurrent := hmaxSuffix leftPivot (by omega) hpivotBefore
      have hlast : middle.size - 1 < middle.size := by omega
      have hlastUpper := hmiddleUpper (middle.size - 1) hlast
      rw [getElem!_pos indices leftPivot hleftPivot, hmaxCurrent] at hleftValue
      rw [hsizeMiddle, hmiddleCount] at hleftValue
      by_cases hleftLast : leftPivot = middle.size - 1
      · subst leftPivot
        omega
      · have hgapToLast := hmiddleGap leftPivot (middle.size - 1)
          (by omega) hlast (by omega)
        omega
    · have hleftEq : leftPivot = pivot := by omega
      subst leftPivot
      have hmiddlePivotLower : indices[pivot]! + 1 ≤ middle[pivot]! := by omega
      by_cases hrightBefore : rightPivot < pivot
      · have hindicesMiddleAtRight := hprefixMiddle rightPivot hrightBefore
        have hnextIndicesAtRight := hprefixNext rightPivot hrightBefore
        omega
      · by_cases hpivotRight : pivot < rightPivot
        · have hmiddleNextAtPivot := hmiddleNextPrefix pivot hpivotRight
          rw [hpivotNext] at hmiddleNextAtPivot
          have hmiddlePivotEq : middle[pivot]! = indices[pivot]! + 1 := by omega
          have hrightAfterPivot : pivot < rightPivot := hpivotRight
          have hmiddleLower := hmiddleGap pivot rightPivot hpivotMiddle
            (by omega) hrightAfterPivot
          have hnextExact := hminSuffix rightPivot (by omega) hrightAfterPivot
          rw [hmiddlePivotEq] at hmiddleLower
          rw [hnextExact] at hrightValue
          omega
        · have hrightEq : rightPivot = pivot := by omega
          subst rightPivot
          rw [hpivotNext] at hrightValue
          omega

private theorem initialCombination_legal (upper count : Nat)
    (hfits : count ≤ upper) :
    LegalCombination upper count
      (Generated.StrictRecombine.initialCombination count) := by
  refine ⟨initialCombination_size count, ?_, ?_⟩
  · intro left right hleft hright hlr
    have hleftCount : left < count := by
      simpa [initialCombination_size] using hleft
    have hrightCount : right < count := by
      simpa [initialCombination_size] using hright
    have hleftArray :
        (Generated.StrictRecombine.initialCombination count)[left]! =
          (Generated.StrictRecombine.initialCombination count).toList[left]! := by
      rw [getElem!_pos _ left hleft,
        getElem!_pos _ left (by simpa using hleft)]
      exact (Array.getElem_toList hleft).symm
    have hrightArray :
        (Generated.StrictRecombine.initialCombination count)[right]! =
          (Generated.StrictRecombine.initialCombination count).toList[right]! := by
      rw [getElem!_pos _ right hright,
        getElem!_pos _ right (by simpa using hright)]
      exact (Array.getElem_toList hright).symm
    rw [hleftArray, hrightArray, initialCombination_toList]
    simp [getElem!_pos, hleftCount, hrightCount]
    omega
  · intro index hindex
    have hindexCount : index < count := by
      simpa [initialCombination_size] using hindex
    have harray :
        (Generated.StrictRecombine.initialCombination count)[index]! =
          (Generated.StrictRecombine.initialCombination count).toList[index]! := by
      rw [getElem!_pos _ index hindex,
        getElem!_pos _ index (by simpa using hindex)]
      exact (Array.getElem_toList hindex).symm
    rw [harray, initialCombination_toList]
    simp [getElem!_pos, hindexCount]
    omega

private theorem nextCombination_preserves_legal (upper count : Nat)
    (indices next : Array Nat) (hlegal : LegalCombination upper count indices)
    (hrun : Generated.StrictRecombine.nextCombination indices upper =
      (true, next)) : LegalCombination upper count next := by
  obtain ⟨pivot, hpivot, hpivotNotMax, hpivotNext, hprefixNext, hmaxSuffix,
    hminSuffix⟩ :=
    nextCombination_true_minimal_suffix indices next upper hrun
  rcases hlegal with ⟨hsize, hgap, hupper⟩
  have hnextSizeRaw := nextCombination_size indices upper
  rw [hrun] at hnextSizeRaw
  have hnextSize : next.size = indices.size := by simpa using hnextSizeRaw
  have hfits : indices.size ≤ upper := by
    have hcopy := hrun
    unfold Generated.StrictRecombine.nextCombination at hcopy
    split at hcopy
    next hsizeFits => exact hsizeFits
    next hsizeFits => simp at hcopy
  have hpivotNextBounds : pivot < next.size := by omega
  have hpivotBound := legalCombination_positional_bound upper count indices
    ⟨hsize, hgap, hupper⟩ pivot hpivot
  refine ⟨hnextSize.trans hsize, ?_, ?_⟩
  · intro left right hleft hright hlr
    have hleftOld : left < indices.size := by omega
    have hrightOld : right < indices.size := by omega
    by_cases hrightPrefix : right < pivot
    · rw [hprefixNext left (by omega), hprefixNext right hrightPrefix]
      exact hgap left right hleftOld hrightOld hlr
    · by_cases hleftPrefix : left < pivot
      · rw [hprefixNext left hleftPrefix]
        by_cases hrightPivot : right = pivot
        · subst right
          rw [hpivotNext]
          have := hgap left pivot hleftOld hpivot hleftPrefix
          omega
        · rw [hminSuffix right hrightOld (by omega)]
          have := hgap left pivot hleftOld hpivot hleftPrefix
          omega
      · have hpivotLeLeft : pivot ≤ left := by omega
        by_cases hleftPivot : left = pivot
        · subst left
          rw [hpivotNext, hminSuffix right hrightOld hlr]
        · rw [hminSuffix left hleftOld (by omega),
            hminSuffix right hrightOld (by omega)]
          omega
  · intro index hindex
    have hindexOld : index < indices.size := by omega
    by_cases hprefix : index < pivot
    · rw [hprefixNext index hprefix]
      exact hupper index hindexOld
    · by_cases hpivotIndex : index = pivot
      · subst index
        rw [hpivotNext]
        rw [hsize] at hpivotNotMax
        have hpivotCount : pivot < count := by omega
        omega
      · rw [hminSuffix index hindexOld (by omega)]
        rw [hsize] at hpivotNotMax
        have hindexCount : index < count := by omega
        omega

private theorem nextCombination_valid (upper count : Nat)
    (indices next : Array Nat) (hvalid : ValidCombination upper count indices)
    (hrun : Generated.StrictRecombine.nextCombination indices upper =
      (true, next)) :
    ValidCombination upper count next := by
  rcases hvalid with ⟨hsize, hbound⟩
  unfold Generated.StrictRecombine.nextCombination at hrun
  split at hrun
  next hfits =>
    split at hrun
    next hpivotNone => simp at hrun
    next pivot hpivotSome =>
      split at hrun
      next hpivotBounds =>
        have hout := Prod.mk.inj hrun
        rw [← hout.2]
        constructor
        · simpa [resetCombinationSuffix_size, hsize]
        · apply resetCombinationSuffix_preserves_bounds
          intro index hindex
          rw [getElem!_pos _ index hindex, Array.getElem_set]
          split
          next heq =>
            subst index
            have hpivotLt := nextCombinationPivot_some_lt indices upper 0 pivot
              hpivotSome
            have hne := nextCombinationPivot_some_ne indices upper 0 pivot
              hpivotSome
            have hb := hbound pivot hpivotLt
            rw [getElem!_pos indices pivot hpivotLt] at hb
            rw [hsize] at hne
            omega
          next heq =>
            have hi : index < indices.size := by simpa using hindex
            have hb := hbound index hi
            rw [getElem!_pos indices index hi] at hb
            exact hb
      next hpivotBounds => simp at hrun
  next hfits => simp at hrun

private theorem validCombination_digits (upper count : Nat)
    (indices : Array Nat) (hvalid : ValidCombination upper count indices)
    (hfits : count ≤ upper) :
    ∀ digit ∈ indices.toList, digit < upper + 1 := by
  intro digit hdigit
  obtain ⟨index, hindex, hvalue⟩ := List.mem_iff_getElem.mp hdigit
  rcases hvalid with ⟨hsize, hbound⟩
  have hi : index < indices.size := by simpa using hindex
  have hb := hbound index hi
  rw [getElem!_pos indices index hi] at hb
  have harray := Array.getElem_toList (xs := indices) hi
  rw [hvalue] at harray
  rw [← harray] at hb
  omega

private theorem nextCombination_code_increases (upper count : Nat)
    (indices next : Array Nat) (hvalid : ValidCombination upper count indices)
    (hrun : Generated.StrictRecombine.nextCombination indices upper =
      (true, next)) :
    positionalCode (upper + 1) indices.toList <
      positionalCode (upper + 1) next.toList := by
  have hfits : count ≤ upper := by
    have hcopy := hrun
    unfold Generated.StrictRecombine.nextCombination at hcopy
    split at hcopy
    next hsizeFits => exact hvalid.1 ▸ hsizeFits
    next hsizeFits => simp at hcopy
  have hnextValid := nextCombination_valid upper count indices next hvalid hrun
  obtain ⟨pivot, hpivot, hpivotValue, hprefixValues⟩ :=
    nextCombination_true_pivot indices next upper hrun
  have hnextSize := nextCombination_size indices upper
  rw [hrun] at hnextSize
  have hnextSize' : next.size = indices.size := by simpa using hnextSize
  have hpivotNext : pivot < next.size := by omega
  apply positionalCode_array_pivot_lt (upper + 1) (by omega) indices next pivot
    hpivot hpivotNext hnextSize'.symm
  · apply List.ext_getElem
    · simp [List.length_take]
      omega
    · intro index hleft hright
      have hi : index < pivot := by
        simp [List.length_take] at hleft
        exact hleft.1
      rw [List.getElem_take, List.getElem_take,
        Array.getElem_toList (xs := indices) (by omega),
        Array.getElem_toList (xs := next) (by omega)]
      have heq := hprefixValues index hi
      rw [getElem!_pos indices index (by omega),
        getElem!_pos next index (by omega)] at heq
      exact heq.symm
  · omega
  · exact validCombination_digits upper count indices hvalid hfits
  · exact validCombination_digits upper count next hnextValid hfits

private theorem initialCombination_valid (upper count : Nat)
    (hfits : count ≤ upper) :
    ValidCombination upper count
      (Generated.StrictRecombine.initialCombination count) := by
  constructor
  · exact initialCombination_size count
  · intro index hindex
    have hindexCount : index < count := by
      simpa [initialCombination_size] using hindex
    have harrayList :
        (Generated.StrictRecombine.initialCombination count)[index]! =
          (Generated.StrictRecombine.initialCombination count).toList[index]! := by
      rw [getElem!_pos _ index hindex,
        getElem!_pos _ index (by simpa using hindex)]
      exact (Array.getElem_toList hindex).symm
    have htoList := congrArg (fun values : List Nat => values[index]!)
      (initialCombination_toList count)
    dsimp at htoList
    have hrange : (List.range count)[index]! = index := by
      rw [getElem!_pos _ index (by simpa using hindexCount)]
      simp
    rw [harrayList, htoList, hrange]
    omega

private def concreteCombinationTermination (upper count : Nat) :
    Generated.StrictRecombine.CombinationTermination upper count where
  valid := LegalCombination upper count
  valid_size := fun _ hvalid => hvalid.1
  rank := fun current =>
    (upper + 1) ^ count - positionalCode (upper + 1) current.toList
  next_valid := fun current next hvalid hrun =>
    nextCombination_preserves_legal upper count current next hvalid hrun
  next_decreases := by
    intro current next hvalid hrun
    have hfits : count ≤ upper := by
      have hcopy := hrun
      unfold Generated.StrictRecombine.nextCombination at hcopy
      split at hcopy
      next hsizeFits => exact hvalid.1 ▸ hsizeFits
      next hsizeFits => simp at hcopy
    have hvalidWeak := legalCombination_valid upper count current hvalid
    have hnextLegal := nextCombination_preserves_legal upper count current next
      hvalid hrun
    have hnextValid := legalCombination_valid upper count next hnextLegal
    have hcode := nextCombination_code_increases upper count current next
      hvalidWeak hrun
    apply Nat.sub_lt_sub_left
    · have hdigits := validCombination_digits upper count current hvalidWeak hfits
      simpa [hvalidWeak.1] using
        positionalCode_lt_pow (upper + 1) (by omega) current.toList hdigits
    · exact hcode

theorem removeCombinationLoop_size (candidate : Array Nat)
    (remaining : Nat) (active output : Array SparsePolyZZ)
    (hrun : Generated.StrictRecombine.removeCombinationLoop candidate remaining
      active = .ok output) :
    output.size + remaining = active.size := by
  induction remaining generalizing active output with
  | zero =>
      rw [Generated.StrictRecombine.removeCombinationLoop] at hrun
      have hout := Except.ok.inj hrun
      subst output
      simp
  | succ remaining ih =>
      rw [Generated.StrictRecombine.removeCombinationLoop] at hrun
      split at hrun
      next hcand =>
        dsimp at hrun
        split at hrun
        next hactive =>
          have htail := ih (active.eraseIdxIfInBounds candidate[remaining])
            output hrun
          simp only [Array.size_eraseIdxIfInBounds, if_pos hactive] at htail
          omega
        next hactive => contradiction
      next hcand => contradiction

theorem removeCombination_strict (candidate : Array Nat)
    (active output : Array SparsePolyZZ) (hnonempty : 0 < candidate.size)
    (hrun : Generated.StrictRecombine.removeCombination candidate active =
      .ok output) :
    output.size < active.size := by
  unfold Generated.StrictRecombine.removeCombination at hrun
  have hsize := removeCombinationLoop_size candidate candidate.size active
    output hrun
  omega

private theorem selectSourceIndices_eraseIdx_of_lt {α : Type*} [Inhabited α]
    (source : List α) (indices : List Nat) (erased : Nat)
    (hlt : ∀ index ∈ indices, index < erased) :
    selectSourceIndices (source.eraseIdx erased) indices =
      selectSourceIndices source indices := by
  unfold selectSourceIndices
  apply List.map_congr_left
  intro index hindex
  simp [List.getElem!_eq_getElem?_getD, List.getElem?_eraseIdx,
    hlt index hindex]

/-- Every legal prefix makes the literal reverse-erasure loop succeed.  The
returned array is obtained by executing the generated function; it is not an
independently chosen complement. -/
theorem removeCombinationLoop_succeeds (candidate : Array Nat)
    (remaining : Nat) (active : Array SparsePolyZZ)
    (hremaining : remaining ≤ candidate.size)
    (hgaps : ∀ left right (hleft : left < remaining)
      (hright : right < remaining), left < right →
        candidate[left] < candidate[right])
    (hbound : ∀ position (hposition : position < remaining),
      candidate[position] < active.size) :
    ∃ output, Generated.StrictRecombine.removeCombinationLoop candidate
      remaining active = .ok output := by
  induction remaining generalizing active with
  | zero =>
      refine ⟨active, ?_⟩
      simp [Generated.StrictRecombine.removeCombinationLoop]
  | succ remaining ih =>
      have hcand : remaining < candidate.size := by omega
      have hcurrent : candidate[remaining] < active.size :=
        hbound remaining (by omega)
      let erased := active.eraseIdxIfInBounds candidate[remaining]
      have hlower : ∀ position (hposition : position < remaining),
          candidate[position] < candidate[remaining] := by
        intro position hposition
        exact hgaps position remaining (by omega) (by omega) hposition
      have herasedSize : erased.size + 1 = active.size := by
        have : 0 < active.size := by omega
        simp only [erased, Array.size_eraseIdxIfInBounds, if_pos hcurrent]
        exact Nat.sub_add_cancel (by omega)
      have hbound' : ∀ position (hposition : position < remaining),
          candidate[position] < erased.size := by
        intro position hposition
        have := hlower position hposition
        omega
      have hgaps' : ∀ left right (hleft : left < remaining)
          (hright : right < remaining), left < right →
            candidate[left] < candidate[right] := by
        intro left right hleft hright hlt
        exact hgaps left right (by omega) (by omega) hlt
      rcases ih erased (by omega) hgaps' hbound' with ⟨output, houtput⟩
      refine ⟨output, ?_⟩
      rw [Generated.StrictRecombine.removeCombinationLoop]
      simp only [hcand, hcurrent, ↓reduceDIte, erased]
      exact houtput

/-- The literal reverse-erasure loop partitions the physical active array:
the occurrence-sensitive product named by the still-live candidate prefix,
times the product of the returned array, is exactly the input product.  The
proof follows every generated erase and uses strict candidate ordering to
show that deleting a later position leaves all earlier selected positions
unchanged. -/
theorem removeCombinationLoop_product_partition (candidate : Array Nat)
    (remaining : Nat) (active output : Array SparsePolyZZ)
    (hremaining : remaining ≤ candidate.size)
    (hgaps : ∀ left right (hleft : left < remaining)
      (hright : right < remaining), left < right →
        candidate[left] < candidate[right])
    (hbound : ∀ position (hposition : position < remaining),
      candidate[position] < active.size)
    (hrun : Generated.StrictRecombine.removeCombinationLoop candidate remaining
      active = .ok output) :
    ((selectSourceIndices active.toList
        (candidate.toList.take remaining)).map SparsePolyZZ.toPoly).prod *
      (output.toList.map SparsePolyZZ.toPoly).prod =
        (active.toList.map SparsePolyZZ.toPoly).prod := by
  induction remaining generalizing active output with
  | zero =>
      rw [Generated.StrictRecombine.removeCombinationLoop] at hrun
      have hout := Except.ok.inj hrun
      subst output
      simp [selectSourceIndices]
  | succ remaining ih =>
      rw [Generated.StrictRecombine.removeCombinationLoop] at hrun
      split at hrun
      next hcand =>
        dsimp at hrun
        split at hrun
        next hactive =>
          let erased := active.eraseIdxIfInBounds candidate[remaining]
          have hcurrent : candidate[remaining] < active.size :=
            hbound remaining (by omega)
          have hlower : ∀ position (hposition : position < remaining),
              candidate[position] < candidate[remaining] := by
            intro position hposition
            exact hgaps position remaining (by omega) (by omega) hposition
          have herasedSize : erased.size + 1 = active.size := by
            have : 0 < active.size := by omega
            simp only [erased, Array.size_eraseIdxIfInBounds, if_pos hcurrent]
            exact Nat.sub_add_cancel (by omega)
          have hbound' : ∀ position (hposition : position < remaining),
              candidate[position] < erased.size := by
            intro position hposition
            have := hlower position hposition
            omega
          have hgaps' : ∀ left right (hleft : left < remaining)
              (hright : right < remaining), left < right →
                candidate[left] < candidate[right] := by
            intro left right hleft hright hlt
            exact hgaps left right (by omega) (by omega) hlt
          have hih := ih erased output (by omega) hgaps' hbound' hrun
          have hselectErase :
              selectSourceIndices erased.toList
                  (candidate.toList.take remaining) =
                selectSourceIndices active.toList
                  (candidate.toList.take remaining) := by
            rw [Array.toList_eraseIdxIfInBounds]
            apply selectSourceIndices_eraseIdx_of_lt
            intro index hindex
            rcases List.mem_iff_getElem.mp hindex with
              ⟨position, hposition, hindexEq⟩
            have hpositionLt : position < remaining := by
              simpa only [List.length_take, Array.length_toList,
                Nat.min_eq_left (show remaining ≤ candidate.size by omega)]
                using hposition
            have hcandidatePos : candidate[position] = index := by
              rw [← Array.getElem_toList (by omega)]
              rw [List.getElem_take] at hindexEq
              exact hindexEq
            rw [← hcandidatePos]
            exact hlower position hpositionLt
          rw [hselectErase] at hih
          have htake : candidate.toList.take (remaining + 1) =
              candidate.toList.take remaining ++ [candidate[remaining]] := by
            rw [← List.take_append_getElem (l := candidate.toList)
              (i := remaining) (by simpa using hcand)]
            simp [Array.getElem_toList hcand]
          rw [htake]
          simp only [selectSourceIndices, List.map_append, List.map_singleton,
            List.prod_append, List.prod_singleton]
          have hselectedIndex : active.toList[candidate[remaining]]! =
              active[candidate[remaining]] := by
            rw [getElem!_pos _ _ (by simpa using hcurrent),
              Array.getElem_toList hcurrent]
          rw [hselectedIndex]
          have heraseProd : SparsePolyZZ.toPoly active[candidate[remaining]] *
                ((erased.toList.map SparsePolyZZ.toPoly).prod) =
              (active.toList.map SparsePolyZZ.toPoly).prod := by
            rw [Array.toList_eraseIdxIfInBounds]
            have hmapBound : candidate[remaining] <
                (active.toList.map SparsePolyZZ.toPoly).length := by
              simpa only [List.length_map, Array.length_toList] using hcurrent
            have hmapIndex :
                (active.toList.map SparsePolyZZ.toPoly)[candidate[remaining]]'hmapBound =
                  SparsePolyZZ.toPoly active[candidate[remaining]] := by
              simp [List.getElem_map, Array.getElem_toList hcurrent]
            have hprod := List.CommMonoid.mul_prod_eraseIdx
              (l := active.toList.map SparsePolyZZ.toPoly)
              (i := candidate[remaining]) (by
                simp only [List.length_map, Array.length_toList]
                exact hcurrent)
            rw [hmapIndex] at hprod
            simpa [List.eraseIdx_map] using hprod
          calc
            _ = SparsePolyZZ.toPoly active[candidate[remaining]] *
                (((selectSourceIndices active.toList
                    (candidate.toList.take remaining)).map
                      SparsePolyZZ.toPoly).prod *
                    (output.toList.map SparsePolyZZ.toPoly).prod) := by
                  ac_rfl
            _ = SparsePolyZZ.toPoly active[candidate[remaining]] *
                (erased.toList.map SparsePolyZZ.toPoly).prod := by rw [hih]
            _ = _ := heraseProd
        next hactive => contradiction
      next hcand => contradiction

/-- Full-array form for the generated successful-candidate removal.  This is
the concrete selected/complement factorization consumed by Hensel uniqueness;
the complement is the actual array returned by C++, not an existential list. -/
theorem removeCombination_product_partition (candidate : Array Nat)
    (active output : Array SparsePolyZZ)
    (hlegal : LegalCombination active.size candidate.size candidate)
    (hrun : Generated.StrictRecombine.removeCombination candidate active =
      .ok output) :
    ((selectSourceIndices active.toList candidate.toList).map
        SparsePolyZZ.toPoly).prod *
      (output.toList.map SparsePolyZZ.toPoly).prod =
        (active.toList.map SparsePolyZZ.toPoly).prod := by
  unfold Generated.StrictRecombine.removeCombination at hrun
  have hgaps : ∀ left right (hleft : left < candidate.size)
      (hright : right < candidate.size), left < right →
        candidate[left] < candidate[right] := by
    intro left right hleft hright hlt
    have hgap := hlegal.2.1 left right hleft hright hlt
    rw [getElem!_pos candidate left hleft,
      getElem!_pos candidate right hright] at hgap
    omega
  have hpartition := removeCombinationLoop_product_partition candidate
    candidate.size active output (Nat.le_refl _) hgaps
    (fun position hposition => by
      simpa [getElem!_pos candidate position hposition] using
        hlegal.2.2 position hposition) hrun
  rw [show candidate.size = candidate.toList.length by simp,
    List.take_length] at hpartition
  exact hpartition

/-- Modular image of the exact physical selected/complement partition. -/
theorem removeCombination_toPolyMod_product_partition (modulus : Nat)
    (candidate : Array Nat) (active output : Array SparsePolyZZ)
    (hlegal : LegalCombination active.size candidate.size candidate)
    (hrun : Generated.StrictRecombine.removeCombination candidate active =
      .ok output) :
    ((selectSourceIndices active.toList candidate.toList).map
        (Refinement.StrictHensel.toPolyMod modulus)).prod *
      (output.toList.map
        (Refinement.StrictHensel.toPolyMod modulus)).prod =
      (active.toList.map
        (Refinement.StrictHensel.toPolyMod modulus)).prod := by
  have hpartition := congrArg
    (Polynomial.map (Int.castRingHom (ZMod modulus)))
    (removeCombination_product_partition candidate active output hlegal hrun)
  simpa [Refinement.StrictHensel.toPolyMod, Polynomial.map_mul,
    Polynomial.map_list_prod, List.map_map] using hpartition

/-- Full legal-candidate entry form of `removeCombinationLoop_succeeds`. -/
theorem removeCombination_succeeds (candidate : Array Nat)
    (active : Array SparsePolyZZ)
    (hlegal : LegalCombination active.size candidate.size candidate) :
    ∃ output, Generated.StrictRecombine.removeCombination candidate active =
      .ok output := by
  have hgaps : ∀ left right (hleft : left < candidate.size)
      (hright : right < candidate.size), left < right →
        candidate[left] < candidate[right] := by
    intro left right hleft hright hlt
    have hgap := hlegal.2.1 left right hleft hright hlt
    rw [getElem!_pos candidate left hleft,
      getElem!_pos candidate right hright] at hgap
    omega
  have hbound : ∀ position (hposition : position < candidate.size),
      candidate[position] < active.size := by
    intro position hposition
    simpa [getElem!_pos candidate position hposition] using
      hlegal.2.2 position hposition
  simpa [Generated.StrictRecombine.removeCombination] using
    removeCombinationLoop_succeeds candidate candidate.size active
      (Nat.le_refl _) hgaps hbound

private theorem mem_toList_of_mem_eraseIdxIfInBounds_toList
    {α : Type*} (value : α) (active : Array α) (index : Nat)
    (hmember : value ∈ (active.eraseIdxIfInBounds index).toList) :
    value ∈ active.toList := by
  grind

/-- Reverse erasure in the generated successful-candidate path only removes
source occurrences; every surviving active factor is an actual member of the
input active array. -/
theorem removeCombinationLoop_member
    (candidate : Array Nat) (remaining : Nat)
    (active output : Array SparsePolyZZ)
    (hrun : Generated.StrictRecombine.removeCombinationLoop candidate
      remaining active = .ok output) :
    ∀ factor ∈ output.toList, factor ∈ active.toList := by
  induction remaining generalizing active output with
  | zero =>
      rw [Generated.StrictRecombine.removeCombinationLoop] at hrun
      have hout := Except.ok.inj hrun
      subst output
      exact fun factor hfactor => hfactor
  | succ remaining ih =>
      rw [Generated.StrictRecombine.removeCombinationLoop] at hrun
      split at hrun
      next hcand =>
        dsimp at hrun
        split at hrun
        next hactive =>
          intro factor hfactor
          exact mem_toList_of_mem_eraseIdxIfInBounds_toList factor active
            candidate[remaining] (ih _ _ hrun factor hfactor)
        next hactive => contradiction
      next hcand => contradiction

theorem removeCombination_member
    (candidate : Array Nat) (active output : Array SparsePolyZZ)
    (hrun : Generated.StrictRecombine.removeCombination candidate active =
      .ok output) :
    ∀ factor ∈ output.toList, factor ∈ active.toList := by
  exact removeCombinationLoop_member candidate candidate.size active output hrun

/-- The physical array returned by the generated reverse-erasure loop is an
occurrence-sensitive sublist of its input.  Unlike `removeCombination_member`,
this retains multiplicity and order, which are needed to re-encode the literal
complement as another legal scan candidate. -/
theorem removeCombinationLoop_sublist
    (candidate : Array Nat) (remaining : Nat)
    (active output : Array SparsePolyZZ)
    (hrun : Generated.StrictRecombine.removeCombinationLoop candidate
      remaining active = .ok output) :
    output.toList.Sublist active.toList := by
  induction remaining generalizing active output with
  | zero =>
      rw [Generated.StrictRecombine.removeCombinationLoop] at hrun
      have hout := Except.ok.inj hrun
      subst output
      exact List.Sublist.refl _
  | succ remaining ih =>
      rw [Generated.StrictRecombine.removeCombinationLoop] at hrun
      split at hrun
      next hcand =>
        dsimp at hrun
        split at hrun
        next hactive =>
          exact (ih _ _ hrun).trans (by
            rw [Array.toList_eraseIdxIfInBounds]
            exact List.eraseIdx_sublist _ _)
        next hactive => contradiction
      next hcand => contradiction

/-- Full generated removal preserves the exact occurrence order of the
surviving complement. -/
theorem removeCombination_sublist
    (candidate : Array Nat) (active output : Array SparsePolyZZ)
    (hrun : Generated.StrictRecombine.removeCombination candidate active =
      .ok output) :
    output.toList.Sublist active.toList := by
  exact removeCombinationLoop_sublist candidate candidate.size active output
    hrun

/-- Executing removal for a legal candidate produces a physical complement
which can itself be encoded as a legal candidate over the original array.
The returned size equation is the exact reverse-erasure accounting, so one of
the selected candidate and this complement has size at most half of `active`.
-/
theorem removeCombination_complement_candidate
    (candidate : Array Nat) (active : Array SparsePolyZZ)
    (hlegal : LegalCombination active.size candidate.size candidate) :
    ∃ output complement,
      Generated.StrictRecombine.removeCombination candidate active =
        .ok output ∧
      LegalCombination active.size output.size complement ∧
      selectSourceIndices active.toList complement.toList = output.toList ∧
      output.size + candidate.size = active.size := by
  rcases removeCombination_succeeds candidate active hlegal with
    ⟨output, houtput⟩
  have hsublist := removeCombination_sublist candidate active output houtput
  rcases sublist_exists_legal_combination hsublist with
    ⟨complement, hcomplement, hselected⟩
  have hsize : output.size + candidate.size = active.size := by
    unfold Generated.StrictRecombine.removeCombination at houtput
    simpa using removeCombinationLoop_size candidate candidate.size active
      output houtput
  exact ⟨output, complement, houtput, by simpa using hcomplement,
    hselected, hsize⟩

/-- Pointwise properties, in particular selected-prime irreducibility, survive
the actual generated reverse-erasure execution. -/
theorem removeCombination_preserves_pointwise
    (candidate : Array Nat) (active output : Array SparsePolyZZ)
    (property : SparsePolyZZ → Prop)
    (hproperty : ∀ index (hindex : index < active.size),
      property active[index])
    (hrun : Generated.StrictRecombine.removeCombination candidate active =
      .ok output) :
    ∀ index (hindex : index < output.size), property output[index] := by
  intro index hindex
  have hmember := removeCombination_member candidate active output hrun
    output[index] (Array.getElem_mem_toList hindex)
  rcases List.mem_iff_getElem.mp hmember with
    ⟨sourceIndex, hsourceIndex, hvalue⟩
  have hsourceArray : sourceIndex < active.size := by simpa using hsourceIndex
  rw [Array.getElem_toList hsourceArray] at hvalue
  rw [← hvalue]
  exact hproperty sourceIndex hsourceArray

/-- Representation and selected-prime facts carried by the physical active
array throughout the generated Zassenhaus outer loop.  This contains only
properties of actual array cells, not a semantic recombination result. -/
structure LiveActiveFactors (base : Nat) (active : Array SparsePolyZZ) : Prop where
  fitsInt32 : active.size ≤ 2 ^ 31
  canonical : ∀ index (hindex : index < active.size),
    StrictPolynomialMod.SparsePolyZZCanonical active[index]
  nonempty : ∀ index (hindex : index < active.size), 0 < active[index].size
  monic : ∀ index (hindex : index < active.size),
    (SparsePolyZZ.toPoly active[index]).Monic
  irreducible : ∀ index (hindex : index < active.size),
    Irreducible (Refinement.StrictHensel.toPolyMod base active[index])

/-- Every physical occurrence-sensitive selection from a live active array is
monic at any observation modulus. -/
theorem LiveActiveFactors.selectedToPolyModMonic
    {base modulus : Nat} {active : Array SparsePolyZZ}
    (state : LiveActiveFactors base active) (candidate : Array Nat)
    (hbound : ∀ position (hposition : position < candidate.size),
      candidate[position] < active.size) :
    ((selectSourceIndices active.toList candidate.toList).map
      (Refinement.StrictHensel.toPolyMod modulus)).prod.Monic := by
  let selected := (selectSourceIndices active.toList candidate.toList).map
    (Refinement.StrictHensel.toPolyMod modulus)
  have hall : ∀ mapped ∈ selected, mapped.Monic := by
    intro mapped hmapped
    rcases List.mem_map.mp hmapped with ⟨factor, hfactor, rfl⟩
    unfold selectSourceIndices at hfactor
    rcases List.mem_map.mp hfactor with ⟨index, hindex, rfl⟩
    rcases List.mem_iff_getElem.mp hindex with
      ⟨position, hposition, hindexEq⟩
    have hpositionArray : position < candidate.size := by simpa using hposition
    have hactive := hbound position hpositionArray
    have hselected : candidate[position] = index := by
      rw [← Array.getElem_toList hpositionArray]
      exact hindexEq
    rw [← hselected, getElem!_pos active.toList candidate[position]
      (by simpa using hactive), Array.getElem_toList hactive]
    simpa [Refinement.StrictHensel.toPolyMod] using
      (state.monic candidate[position] hactive).map
        (Int.castRingHom (ZMod modulus))
  change selected.prod.Monic
  have prodMonic : ∀ list : List (Polynomial (ZMod modulus)),
      (∀ polynomial ∈ list, polynomial.Monic) → list.prod.Monic := by
    intro list hlist
    induction list with
    | nil => simp
    | cons head tail ih =>
        rw [List.prod_cons]
        exact (hlist head (List.mem_cons_self)).mul
          (ih (fun polynomial hpolynomial =>
            hlist polynomial (List.mem_cons_of_mem head hpolynomial)))
  exact prodMonic selected hall

/-- Every physical occurrence-sensitive selection from a live active array is
monic as an integer polynomial. -/
theorem LiveActiveFactors.selectedToPolyMonic
    {base : Nat} {active : Array SparsePolyZZ}
    (state : LiveActiveFactors base active) (candidate : Array Nat)
    (hbound : ∀ position (hposition : position < candidate.size),
      candidate[position] < active.size) :
    ((selectSourceIndices active.toList candidate.toList).map
      SparsePolyZZ.toPoly).prod.Monic := by
  let selected := (selectSourceIndices active.toList candidate.toList).map
    SparsePolyZZ.toPoly
  have hall : ∀ mapped ∈ selected, mapped.Monic := by
    intro mapped hmapped
    rcases List.mem_map.mp hmapped with ⟨factor, hfactor, rfl⟩
    unfold selectSourceIndices at hfactor
    rcases List.mem_map.mp hfactor with ⟨index, hindex, rfl⟩
    rcases List.mem_iff_getElem.mp hindex with
      ⟨position, hposition, hindexEq⟩
    have hpositionArray : position < candidate.size := by simpa using hposition
    have hactive := hbound position hpositionArray
    have hselected : candidate[position] = index := by
      rw [← Array.getElem_toList hpositionArray]
      exact hindexEq
    rw [← hselected, getElem!_pos active.toList candidate[position]
      (by simpa using hactive), Array.getElem_toList hactive]
    exact state.monic candidate[position] hactive
  change selected.prod.Monic
  have prodMonic : ∀ list : List (Polynomial Int),
      (∀ polynomial ∈ list, polynomial.Monic) → list.prod.Monic := by
    intro list hlist
    induction list with
    | nil => simp
    | cons head tail ih =>
        rw [List.prod_cons]
        exact (hlist head (List.mem_cons_self)).mul
          (ih (fun polynomial hpolynomial =>
            hlist polynomial (List.mem_cons_of_mem head hpolynomial)))
  exact prodMonic selected hall

/-- The literal reverse erasure used after an extraction preserves every
field of `LiveActiveFactors`; the size bound follows from the exact generated
array-size equation. -/
theorem LiveActiveFactors.removeCombination
    {base : Nat} {active output : Array SparsePolyZZ}
    (state : LiveActiveFactors base active) (candidate : Array Nat)
    (hrun : Generated.StrictRecombine.removeCombination candidate active =
      .ok output) :
    LiveActiveFactors base output := by
  have hpointwise (property : SparsePolyZZ → Prop)
      (hproperty : ∀ index (hindex : index < active.size),
        property active[index]) :
      ∀ index (hindex : index < output.size), property output[index] :=
    removeCombination_preserves_pointwise candidate active output property
      hproperty hrun
  have hloop := hrun
  unfold Generated.StrictRecombine.removeCombination at hloop
  have hsize := removeCombinationLoop_size candidate candidate.size active
    output hloop
  have houtputLe : output.size ≤ active.size := by omega
  exact {
    fitsInt32 := houtputLe.trans state.fitsInt32
    canonical := hpointwise
      StrictPolynomialMod.SparsePolyZZCanonical state.canonical
    nonempty := hpointwise (fun factor => 0 < factor.size) state.nonempty
    monic := hpointwise (fun factor => (SparsePolyZZ.toPoly factor).Monic)
      state.monic
    irreducible := hpointwise
      (fun factor => Irreducible
        (Refinement.StrictHensel.toPolyMod base factor)) state.irreducible }

/-- The concrete termination package for the source Zassenhaus loops.  Its
combination rank is the complement of the base-`upper+1` positional code;
the outer rank is discharged by the actual successful subset removal. -/
def concreteZassenhausTermination :
    Generated.StrictRecombine.ZassenhausTermination where
  combinations := concreteCombinationTermination
  initial_valid := fun upper count hfits =>
    initialCombination_legal upper count hfits
  removal_decreases := fun active candidate output hnonempty hrun =>
    removeCombination_strict candidate active output hnonempty hrun

/-- Exhaustion of the actual generated fixed-size scan means that every
legal candidate at or after the supplied current candidate was executed and
rejected.  This follows the source recursion and its concrete legality/rank
certificate; it does not quantify over an abstract successful-factor oracle. -/
theorem scanZassenhausCombinations_exhausted_rejects_legal
    {upper count : Nat} (fStar : SparsePolyZZ)
    (activeLifted : Array SparsePolyZZ) (modulus : ZZ)
    (current target : Array Nat)
    (hcurrent : LegalCombination upper count current)
    (htarget : LegalCombination upper count target)
    (hle : current = target ∨ ArrayLexLT current target)
    (hrun : Generated.StrictRecombine.scanZassenhausCombinations
      (concreteCombinationTermination upper count) fStar activeLifted modulus
      current hcurrent = .ok .exhausted) :
    Generated.StrictRecombine.zassenhausAttempt fStar activeLifted modulus
      target = .ok .rejected := by
  let termination := concreteCombinationTermination upper count
  induction hmeasure : termination.rank current using Nat.strong_induction_on
      generalizing current with
  | h measure ih =>
      rw [Generated.StrictRecombine.scanZassenhausCombinations] at hrun
      cases hattempt : Generated.StrictRecombine.zassenhausAttempt fStar
          activeLifted modulus current with
      | error fault =>
          rw [hattempt] at hrun
          contradiction
      | ok attemptResult =>
          cases attemptResult with
          | extracted factor quotient =>
              simp [hattempt] at hrun
          | rejected =>
              simp only [hattempt] at hrun
              rcases hle with rfl | hbefore
              · exact hattempt
              · split at hrun
                next next hnext =>
                  have hfits := legalCombination_count_le upper count current
                    hcurrent
                  have hfinal := nextCombination_false_is_final current next
                    upper count hcurrent.1 hfits hnext
                  exact (finalCombination_not_lt_legal current target upper
                    count htarget hcurrent.1 hfinal.2 hbefore).elim
                next next hnext =>
                  have hnextLegal := nextCombination_preserves_legal upper
                    count current next hcurrent hnext
                  have hnextLe : next = target ∨ ArrayLexLT next target := by
                    rcases arrayLex_trichotomy next target
                        (hnextLegal.1.trans htarget.1.symm) with
                      heq | hforward | hbackward
                    · exact Or.inl heq
                    · exact Or.inr hforward
                    · exact (nextCombination_true_no_legal_between current
                        next target upper count htarget hnext hbefore
                        hbackward).elim
                  exact ih (termination.rank next) (by
                    rw [← hmeasure]
                    exact termination.next_decreases current next hcurrent
                      hnext) next hnextLegal hnextLe hrun rfl

/-- Public no-omission form for the generated scan from its real iota
initial state: an exhausted run has executed and rejected every legal
fixed-size subset candidate. -/
theorem scanZassenhausCombinations_exhausted_rejects_all
    {upper count : Nat} (fStar : SparsePolyZZ)
    (activeLifted : Array SparsePolyZZ) (modulus : ZZ)
    (hfits : count ≤ upper)
    (hrun : Generated.StrictRecombine.scanZassenhausCombinations
      (concreteCombinationTermination upper count) fStar activeLifted modulus
      (Generated.StrictRecombine.initialCombination count)
      (initialCombination_legal upper count hfits) = .ok .exhausted) :
    ∀ target, LegalCombination upper count target →
      Generated.StrictRecombine.zassenhausAttempt fStar activeLifted modulus
        target = .ok .rejected := by
  intro target htarget
  exact scanZassenhausCombinations_exhausted_rejects_legal fStar activeLifted
    modulus (Generated.StrictRecombine.initialCombination count) target
    (initialCombination_legal upper count hfits) htarget
    (initialCombination_lex_le upper count target htarget) hrun

/-- The literal execution proposition for an exhausted generated fixed-size
Zassenhaus scan from its source iota initial combination. -/
def FixedSizeScanExhausted (fStar : SparsePolyZZ)
    (activeLifted : Array SparsePolyZZ) (modulus : ZZ) (count : Nat) : Prop :=
  ∃ hfits : count ≤ activeLifted.size,
    Generated.StrictRecombine.scanZassenhausCombinations
      (concreteCombinationTermination activeLifted.size count)
      fStar activeLifted modulus
      (Generated.StrictRecombine.initialCombination count)
      (initialCombination_legal activeLifted.size count hfits) = .ok .exhausted

/-- Public wrapper exposing the no-omission consequence of the literal
fixed-size generated scan without exporting its private rank construction. -/
theorem FixedSizeScanExhausted.rejects
    {fStar : SparsePolyZZ} {activeLifted : Array SparsePolyZZ}
    {modulus : ZZ} {count : Nat}
    (hrun : FixedSizeScanExhausted fStar activeLifted modulus count) :
    ∀ target, LegalCombination activeLifted.size count target →
      Generated.StrictRecombine.zassenhausAttempt fStar activeLifted modulus
        target = .ok .rejected := by
  rcases hrun with ⟨hfits, hscan⟩
  exact scanZassenhausCombinations_exhausted_rejects_all fStar activeLifted
    modulus hfits hscan

/-- Literal execution history carried while the generated outer loop advances
its subset size.  Every strictly smaller positive size is backed by an actual
`.ok .exhausted` run of the corresponding generated fixed-size scan. -/
def SmallerZassenhausScansExhausted (fStar : SparsePolyZZ)
    (activeLifted : Array SparsePolyZZ) (modulus : ZZ)
    (subsetSize : Nat) : Prop :=
  ∀ count, 0 < count → count < subsetSize →
    FixedSizeScanExhausted fStar activeLifted modulus count

/-- At the source restart size one there are no smaller positive scans. -/
theorem SmallerZassenhausScansExhausted.one
    (fStar : SparsePolyZZ) (activeLifted : Array SparsePolyZZ)
    (modulus : ZZ) :
    SmallerZassenhausScansExhausted fStar activeLifted modulus 1 := by
  intro count hpositive hsmall
  omega

/-- Advancing the source subset size records precisely the fixed-size scan
which just returned `.exhausted`, while retaining all earlier run equations. -/
theorem SmallerZassenhausScansExhausted.succ
    {fStar : SparsePolyZZ} {activeLifted : Array SparsePolyZZ}
    {modulus : ZZ} {subsetSize : Nat}
    (hprevious : SmallerZassenhausScansExhausted fStar activeLifted modulus
      subsetSize)
    (hcurrent : FixedSizeScanExhausted fStar activeLifted modulus subsetSize) :
    SmallerZassenhausScansExhausted fStar activeLifted modulus
      (subsetSize + 1) := by
  intro count hpositive hsmall
  by_cases heq : count = subsetSize
  · simpa [heq] using hcurrent
  · exact hprevious count hpositive (by omega)

/-- Every legal candidate of a previously exhausted smaller size was
physically executed and rejected. -/
theorem SmallerZassenhausScansExhausted.rejects
    {fStar : SparsePolyZZ} {activeLifted : Array SparsePolyZZ}
    {modulus : ZZ} {subsetSize : Nat}
    (hhistory : SmallerZassenhausScansExhausted fStar activeLifted modulus
      subsetSize)
    (candidate : Array Nat)
    (hpositive : 0 < candidate.size)
    (hsmall : candidate.size < subsetSize)
    (hlegal : LegalCombination activeLifted.size candidate.size candidate) :
    Generated.StrictRecombine.zassenhausAttempt fStar activeLifted modulus
      candidate = .ok .rejected :=
  (hhistory candidate.size hpositive hsmall).rejects candidate hlegal

/-- Pure array value computed by the generated source-shaped μ prefix loop.
This is used only to state its exact execution theorem; it is itself strictly
well-founded on the same remaining-prefix measure. -/
def reduceMuPrefixArray (mu : Generated.StrictRecombine.QQMatrix)
    (k source : Nat) (q : ZZ) (limit index : Nat) :
    Generated.StrictRecombine.QQMatrix :=
  if hindex : index < limit then
    let value := (mu[k]!)[index]! - (q : QQ) * (mu[source]!)[index]!
    reduceMuPrefixArray
      (mu.setIfInBounds k (mu[k]!.setIfInBounds index value))
      k source q limit (index + 1)
  else mu
termination_by limit - index
decreasing_by omega

theorem reduceMuPrefixLoop_eq_array
    (mu : Generated.StrictRecombine.QQMatrix)
    (k source : Nat) (q : ZZ) (limit index : Nat)
    (hk : k < mu.size) (hsource : source < mu.size)
    (hlimitK : limit ≤ mu[k].size)
    (hlimitSource : limit ≤ mu[source].size)
    (hne : k ≠ source) :
    Generated.StrictRecombine.reduceMuPrefixLoop mu k source q limit index =
      .ok (reduceMuPrefixArray mu k source q limit index) := by
  induction hmeasure : limit - index using Nat.strong_induction_on
      generalizing mu index with
  | h measure ih =>
      rw [Generated.StrictRecombine.reduceMuPrefixLoop, reduceMuPrefixArray]
      by_cases hindex : index < limit
      · rw [dif_pos hindex, dif_pos hk, dif_pos hsource]
        have hindexK : index < mu[k].size := lt_of_lt_of_le hindex hlimitK
        have hindexSource : index < mu[source].size :=
          lt_of_lt_of_le hindex hlimitSource
        rw [dif_pos hindexK, dif_pos hindexSource]
        dsimp only
        rw [getElem!_pos mu k hk, getElem!_pos mu source hsource,
          getElem!_pos mu[k] index hindexK,
          getElem!_pos mu[source] index hindexSource]
        let value := mu[k][index] - (q : QQ) * mu[source][index]
        let nextMu := mu.set k (mu[k].set index value)
        have hnextSize : nextMu.size = mu.size := by simp [nextMu]
        have hkNext : k < nextMu.size := by simpa [hnextSize]
        have hsourceNext : source < nextMu.size := by simpa [hnextSize]
        have hkRowNext : nextMu[k].size = mu[k].size := by
          simp [nextMu, hk, value]
        have hsourceRowNext : nextMu[source].size = mu[source].size := by
          simp [nextMu, hk, hsource, hne, value]
        have htail := ih (limit - (index + 1)) (by omega)
          nextMu (index + 1) hkNext hsourceNext
          (by simpa [hkRowNext]) (by simpa [hsourceRowNext]) rfl
        simpa [dif_pos hindex, nextMu, value, Array.setIfInBounds,
          hk, hindexK] using htail
      · rw [dif_neg hindex, dif_neg hindex]

theorem reduceMuPrefixArray_size
    (mu : Generated.StrictRecombine.QQMatrix)
    (k source : Nat) (q : ZZ) (limit index : Nat) :
    (reduceMuPrefixArray mu k source q limit index).size = mu.size := by
  induction hmeasure : limit - index using Nat.strong_induction_on
      generalizing mu index with
  | h measure ih =>
      rw [reduceMuPrefixArray]
      split
      next hindex =>
        dsimp only
        rw [ih (limit - (index + 1)) (by omega) _ (index + 1) rfl]
        simp
      next hindex => rfl

theorem reduceMuPrefixArray_get_of_ne
    (mu : Generated.StrictRecombine.QQMatrix)
    (k source row : Nat) (q : ZZ) (limit index : Nat)
    (hrowBound : row < mu.size) (hrow : row ≠ k) :
    (reduceMuPrefixArray mu k source q limit index)[row]! = mu[row]! := by
  induction hmeasure : limit - index using Nat.strong_induction_on
      generalizing mu index with
  | h measure ih =>
      rw [reduceMuPrefixArray]
      split
      next hindex =>
        dsimp only
        rw [ih (limit - (index + 1)) (by omega) _ (index + 1)
          (by simpa) rfl]
        by_cases hk : k < mu.size
        · rw [getElem!_pos _ row (by simpa),
            getElem!_pos mu row hrowBound]
          simp only [Array.setIfInBounds, dif_pos hk]
          rw [Array.getElem_set]
          simp [Ne.symm hrow]
        · simp [Array.setIfInBounds, hk]
      next hindex => rfl

theorem reduceMuPrefixArray_row_size
    (mu : Generated.StrictRecombine.QQMatrix)
    (k source row : Nat) (q : ZZ) (limit index : Nat) :
    row < mu.size →
    (reduceMuPrefixArray mu k source q limit index)[row]!.size =
      mu[row]!.size := by
  intro hrowBound
  by_cases hrow : row = k
  · subst row
    induction hmeasure : limit - index using Nat.strong_induction_on
        generalizing mu index with
    | h measure ih =>
        rw [reduceMuPrefixArray]
        split
        next hindex =>
          dsimp only
          rw [ih (limit - (index + 1)) (by omega) _ (index + 1)
            (by simpa) rfl]
          by_cases hk : k < mu.size
          · by_cases hi : index < mu[k].size
            · simp [Array.setIfInBounds, hk, hi]
            · simp [Array.setIfInBounds, hk, hi]
          · simp [Array.setIfInBounds, hk]
        next hindex => rfl
  · rw [reduceMuPrefixArray_get_of_ne mu k source row q limit index
      hrowBound hrow]

theorem reduceMuPrefixArray_target_get
    (mu : Generated.StrictRecombine.QQMatrix)
    (k source : Nat) (q : ZZ) (limit index column : Nat)
    (hk : k < mu.size) (hsource : source < mu.size)
    (hlimitK : limit ≤ mu[k].size)
    (hlimitSource : limit ≤ mu[source].size)
    (hcolumn : column < mu[k].size) (hne : k ≠ source) :
    ((reduceMuPrefixArray mu k source q limit index)[k]!)[column]! =
      if index ≤ column ∧ column < limit then
        (mu[k]!)[column]! - (q : QQ) * (mu[source]!)[column]!
      else (mu[k]!)[column]! := by
  induction hmeasure : limit - index using Nat.strong_induction_on
      generalizing mu index with
  | h measure ih =>
      rw [reduceMuPrefixArray]
      by_cases hindex : index < limit
      · rw [dif_pos hindex]
        dsimp only
        have hindexK : index < mu[k].size := lt_of_lt_of_le hindex hlimitK
        have hindexSource : index < mu[source].size :=
          lt_of_lt_of_le hindex hlimitSource
        let value := mu[k][index] - (q : QQ) * mu[source][index]
        let nextMu := mu.set k (mu[k].set index value)
        have hnextSize : nextMu.size = mu.size := by simp [nextMu]
        have hkNext : k < nextMu.size := by simpa [hnextSize]
        have hsourceNext : source < nextMu.size := by simpa [hnextSize]
        have hkRowNext : nextMu[k].size = mu[k].size := by
          simp [nextMu, value]
        have hsourceRowNext : nextMu[source].size = mu[source].size := by
          simp [nextMu, hsource, hne, value]
        have htail := ih (limit - (index + 1)) (by omega)
          nextMu (index + 1) hkNext hsourceNext
          (by simpa [hkRowNext]) (by simpa [hsourceRowNext])
          (by simpa [hkRowNext]) rfl
        simp only [Array.setIfInBounds, dif_pos hk, getElem!_pos mu k hk,
          dif_pos hindexK, getElem!_pos mu[k] index hindexK,
          getElem!_pos mu source hsource,
          getElem!_pos mu[source] index hindexSource] at ⊢
        dsimp only [nextMu, value] at htail
        have houterK :
            (mu.set k (mu[k].set index value))[k]! =
              mu[k].set index value := by
          rw [getElem!_pos _ k (by simpa using hk), Array.getElem_set]
          simp
        have houterSource :
            (mu.set k (mu[k].set index value))[source]! = mu[source] := by
          rw [getElem!_pos _ source (by simp [hsource]), Array.getElem_set]
          simp only [if_neg hne]
        by_cases hcolumnIndex : column = index
        · subst column
          simp only [Nat.not_succ_le_self, false_and, if_false] at htail
          rw [htail]
          simp [nextMu, value, hk, hindexK, hindex, hindexSource]
        · by_cases hindexColumn : index < column
          · have hsuccColumn : index + 1 ≤ column := by omega
            rw [htail]
            simp only [hsuccColumn, true_and]
            rw [houterK, houterSource]
            have hsetGet :
                (mu[k].set index value)[column]! = mu[k][column]! := by
              rw [getElem!_pos _ column (by simpa), Array.getElem_set]
              rw [if_neg (Ne.symm hcolumnIndex)]
              exact (getElem!_pos mu[k] column hcolumn).symm
            rw [hsetGet]
            by_cases hcolumnLimit : column < limit
            · rw [if_pos hcolumnLimit,
                if_pos ⟨Nat.le_of_lt hindexColumn, hcolumnLimit⟩]
            · rw [if_neg hcolumnLimit,
                if_neg (fun condition => hcolumnLimit condition.2)]
          · have hcolumnIndexLt : column < index := by omega
            rw [htail]
            rw [houterK, houterSource]
            have hsetGet :
                (mu[k].set index value)[column]! = mu[k][column]! := by
              rw [getElem!_pos _ column (by simpa), Array.getElem_set]
              rw [if_neg (Ne.symm hcolumnIndex)]
              exact (getElem!_pos mu[k] column hcolumn).symm
            rw [hsetGet]
            rw [if_neg (by omega), if_neg (by omega)]
      · rw [dif_neg hindex]
        have hcond : ¬(index ≤ column ∧ column < limit) := by omega
        simp [hcond]

theorem reduceMuPrefixLoop_output_eq
    (mu output : Generated.StrictRecombine.QQMatrix)
    (k source : Nat) (q : ZZ) (limit index : Nat)
    (hk : k < mu.size) (hsource : source < mu.size)
    (hlimitK : limit ≤ mu[k].size)
    (hlimitSource : limit ≤ mu[source].size)
    (hne : k ≠ source)
    (hrun : Generated.StrictRecombine.reduceMuPrefixLoop
      mu k source q limit index = .ok output) :
    output = reduceMuPrefixArray mu k source q limit index := by
  have hexact := reduceMuPrefixLoop_eq_array mu k source q limit index
    hk hsource hlimitK hlimitSource hne
  rw [hrun] at hexact
  exact Except.ok.inj hexact

theorem reduceMuPrefixLoop_target_get
    (mu output : Generated.StrictRecombine.QQMatrix)
    (k source : Nat) (q : ZZ) (limit index column : Nat)
    (hk : k < mu.size) (hsource : source < mu.size)
    (hlimitK : limit ≤ mu[k].size)
    (hlimitSource : limit ≤ mu[source].size)
    (hcolumn : column < mu[k].size) (hne : k ≠ source)
    (hrun : Generated.StrictRecombine.reduceMuPrefixLoop
      mu k source q limit index = .ok output) :
    (output[k]!)[column]! =
      if index ≤ column ∧ column < limit then
        (mu[k]!)[column]! - (q : QQ) * (mu[source]!)[column]!
      else (mu[k]!)[column]! := by
  rw [reduceMuPrefixLoop_output_eq mu output k source q limit index
    hk hsource hlimitK hlimitSource hne hrun]
  exact reduceMuPrefixArray_target_get mu k source q limit index column
    hk hsource hlimitK hlimitSource hcolumn hne

/-- One exact row update used by the generated post-swap mu loop. -/
def updateMuAfterSwapRow (input : Array QQ) (k : Nat)
    (muOld muNew : QQ) : Array QQ :=
  let oldAtK := input[k]!
  let newAtK := input[k - 1]! - muOld * oldAtK
  let updated := input.setIfInBounds k newAtK
  updated.setIfInBounds (k - 1) (oldAtK + muNew * newAtK)

theorem updateMuAfterSwapRow_size (input : Array QQ) (k : Nat)
    (muOld muNew : QQ) :
    (updateMuAfterSwapRow input k muOld muNew).size = input.size := by
  unfold updateMuAfterSwapRow
  simp [Array.size_setIfInBounds]

theorem updateMuAfterSwapRow_get (input : Array QQ) (k column : Nat)
    (muOld muNew : QQ) (hk : k < input.size)
    (hpred : k - 1 < input.size) (hkPositive : 0 < k)
    (hcolumn : column < input.size) :
    (updateMuAfterSwapRow input k muOld muNew)[column]! =
      if column = k - 1 then
        input[k]! + muNew * (input[k - 1]! - muOld * input[k]!)
      else if column = k then input[k - 1]! - muOld * input[k]!
      else input[column]! := by
  unfold updateMuAfterSwapRow
  simp only [Array.setIfInBounds, dif_pos hk]
  rw [dif_pos (show k - 1 <
    (input.set k (input[k - 1]! - muOld * input[k]!) hk).size by
      simpa using hpred)]
  rw [getElem!_pos _ column (by simpa), Array.getElem_set]
  by_cases hcolumnPred : column = k - 1
  · subst column
    have hkPred : k ≠ k - 1 := by omega
    simp [hkPred, hkPred.symm]
  ·
    rw [getElem!_pos _ column (by simpa), Array.getElem_set]
    by_cases hcolumnK : column = k
    · subst column
      have hkPred : k ≠ k - 1 := by omega
      simp [hkPred, hkPred.symm]
    · have hpredColumn : k - 1 ≠ column := Ne.symm hcolumnPred
      have hkColumn : k ≠ column := Ne.symm hcolumnK
      simp [hcolumnPred, hcolumnK, hpredColumn, hkColumn]

/-- Pure array value computed by the generated post-swap mu correction loop. -/
def updateMuAfterSwapArray (mu : Generated.StrictRecombine.QQMatrix)
    (k : Nat) (muOld muNew : QQ) (row : Nat) :
    Generated.StrictRecombine.QQMatrix :=
  if hrow : row < mu.size then
    updateMuAfterSwapArray (mu.setIfInBounds row
      (updateMuAfterSwapRow mu[row]! k muOld muNew))
      k muOld muNew (row + 1)
  else mu
termination_by mu.size - row
decreasing_by
  simp [Array.size_setIfInBounds]
  omega

theorem updateMuAfterSwapLoop_eq_array
    (mu : Generated.StrictRecombine.QQMatrix)
    (k : Nat) (muOld muNew : QQ) (row : Nat)
    (hk : ∀ index (hindex : index < mu.size), k < mu[index].size)
    (hpred : ∀ index (hindex : index < mu.size),
      k - 1 < mu[index].size) :
    Generated.StrictRecombine.updateMuAfterSwapLoop
      mu k muOld muNew row =
        .ok (updateMuAfterSwapArray mu k muOld muNew row) := by
  induction hmeasure : mu.size - row using Nat.strong_induction_on
      generalizing mu row with
  | h measure ih =>
      rw [Generated.StrictRecombine.updateMuAfterSwapLoop,
        updateMuAfterSwapArray]
      by_cases hrow : row < mu.size
      · rw [dif_pos hrow, dif_pos (hk row hrow), dif_pos (hpred row hrow)]
        dsimp only
        rw [getElem!_pos mu row hrow]
        let oldAtK := (mu[row]'hrow)[k]'(hk row hrow)
        let newAtK := (mu[row]'hrow)[k - 1]'(hpred row hrow) -
          muOld * oldAtK
        let updated := ((mu[row]'hrow).set k newAtK (hk row hrow)).set (k - 1)
          (oldAtK + muNew * newAtK) (by
            rw [Array.size_set]
            exact hpred row hrow)
        let nextMu := mu.set row updated
        have hnextSize : nextMu.size = mu.size := by simp [nextMu]
        have hnextRows : ∀ index (hindex : index < nextMu.size),
            nextMu[index].size = mu[index].size := by
          intro index hindex
          by_cases hindexRow : index = row
          · subst index
            simp [nextMu, updated]
          · have hindexMu : index < mu.size := by simpa [hnextSize] using hindex
            simp [nextMu, hindexMu, Ne.symm hindexRow]
        have hkNext : ∀ index (hindex : index < nextMu.size),
            k < nextMu[index].size := by
          intro index hindex
          rw [hnextRows index hindex]
          exact hk index (by simpa [hnextSize] using hindex)
        have hpredNext : ∀ index (hindex : index < nextMu.size),
            k - 1 < nextMu[index].size := by
          intro index hindex
          rw [hnextRows index hindex]
          exact hpred index (by simpa [hnextSize] using hindex)
        have htail := ih (nextMu.size - (row + 1)) (by
          rw [hnextSize]
          omega) nextMu (row + 1) hkNext hpredNext rfl
        simpa [nextMu, updated, newAtK, oldAtK, updateMuAfterSwapRow,
          Array.setIfInBounds,
          hrow, hk row hrow, hpred row hrow] using htail
      · rw [dif_neg hrow, dif_neg hrow]

theorem updateMuAfterSwapArray_size
    (mu : Generated.StrictRecombine.QQMatrix)
    (k : Nat) (muOld muNew : QQ) (row : Nat) :
    (updateMuAfterSwapArray mu k muOld muNew row).size = mu.size := by
  induction hmeasure : mu.size - row using Nat.strong_induction_on
      generalizing mu row with
  | h measure ih =>
      rw [updateMuAfterSwapArray]
      split
      next hrow =>
        rw [ih _ (by simp; omega) _ (row + 1) rfl]
        simp
      next hrow => rfl

theorem updateMuAfterSwapArray_get_before
    (mu : Generated.StrictRecombine.QQMatrix)
    (k target : Nat) (muOld muNew : QQ) (row : Nat)
    (htarget : target < mu.size) (hbefore : target < row) :
    (updateMuAfterSwapArray mu k muOld muNew row)[target]! = mu[target]! := by
  induction hmeasure : mu.size - row using Nat.strong_induction_on
      generalizing mu row with
  | h measure ih =>
      rw [updateMuAfterSwapArray]
      split
      next hrow =>
        have hnextSize :
            (mu.setIfInBounds row
              (updateMuAfterSwapRow mu[row]! k muOld muNew)).size = mu.size := by
          simp
        rw [ih _ (by simp [hnextSize]; omega) _ (row + 1)
          (by simpa [hnextSize]) (by omega) rfl]
        simp only [Array.setIfInBounds, dif_pos hrow]
        rw [getElem!_pos _ target (by simpa), Array.getElem_set]
        rw [if_neg (by omega)]
        exact (getElem!_pos mu target htarget).symm
      next hrow => rfl

theorem updateMuAfterSwapArray_get_at
    (mu : Generated.StrictRecombine.QQMatrix)
    (k target : Nat) (muOld muNew : QQ) (row : Nat)
    (htarget : target < mu.size) (hrowTarget : row ≤ target) :
    (updateMuAfterSwapArray mu k muOld muNew row)[target]! =
      updateMuAfterSwapRow mu[target]! k muOld muNew := by
  induction hmeasure : mu.size - row using Nat.strong_induction_on
      generalizing mu row with
  | h measure ih =>
      rw [updateMuAfterSwapArray]
      have hrow : row < mu.size := lt_of_le_of_lt hrowTarget htarget
      rw [dif_pos hrow]
      let nextMu := mu.set row
        (updateMuAfterSwapRow mu[row] k muOld muNew)
      have hnextSize : nextMu.size = mu.size := by simp [nextMu]
      have htargetNext : target < nextMu.size := by simpa [hnextSize]
      simp only [Array.setIfInBounds, dif_pos hrow,
        getElem!_pos mu row hrow] at ⊢
      change (updateMuAfterSwapArray nextMu k muOld muNew
        (row + 1))[target]! = updateMuAfterSwapRow mu[target]! k muOld muNew
      by_cases htargetRow : target = row
      · subst target
        have htail := updateMuAfterSwapArray_get_before nextMu k row
          muOld muNew (row + 1) htargetNext (by omega)
        rw [htail]
        rw [getElem!_pos nextMu row htargetNext]
        simp [nextMu]
        rw [getElem!_pos mu row hrow]
      · have hrowTarget' : row + 1 ≤ target := by omega
        have htail := ih (nextMu.size - (row + 1)) (by
          rw [hnextSize]
          omega) nextMu (row + 1) htargetNext hrowTarget' rfl
        rw [htail]
        have hnextTarget : nextMu[target]! = mu[target]! := by
          rw [getElem!_pos nextMu target htargetNext]
          simp [nextMu, htarget, Ne.symm htargetRow]
        rw [hnextTarget]

theorem updateMuAfterSwapLoop_output_eq
    (mu output : Generated.StrictRecombine.QQMatrix)
    (k : Nat) (muOld muNew : QQ) (row : Nat)
    (hk : ∀ index (hindex : index < mu.size), k < mu[index].size)
    (hpred : ∀ index (hindex : index < mu.size),
      k - 1 < mu[index].size)
    (hrun : Generated.StrictRecombine.updateMuAfterSwapLoop
      mu k muOld muNew row = .ok output) :
    output = updateMuAfterSwapArray mu k muOld muNew row := by
  have hexact := updateMuAfterSwapLoop_eq_array mu k muOld muNew row hk hpred
  rw [hrun] at hexact
  exact Except.ok.inj hexact

theorem updateMuAfterSwapLoop_get_at
    (mu output : Generated.StrictRecombine.QQMatrix)
    (k target : Nat) (muOld muNew : QQ) (row : Nat)
    (hk : ∀ index (hindex : index < mu.size), k < mu[index].size)
    (hpred : ∀ index (hindex : index < mu.size),
      k - 1 < mu[index].size)
    (htarget : target < mu.size) (hrowTarget : row ≤ target)
    (hrun : Generated.StrictRecombine.updateMuAfterSwapLoop
      mu k muOld muNew row = .ok output) :
    output[target]! = updateMuAfterSwapRow mu[target]! k muOld muNew := by
  rw [updateMuAfterSwapLoop_output_eq mu output k muOld muNew row hk hpred hrun]
  exact updateMuAfterSwapArray_get_at mu k target muOld muNew row
    htarget hrowTarget

theorem swapQQRows_ok
    (matrix output : Generated.StrictRecombine.QQMatrix) (left right : Nat)
    (hrun : Generated.StrictRecombine.swapQQRows matrix left right =
      .ok output) :
    ∃ (hleft : left < matrix.size) (hright : right < matrix.size),
      output = (matrix.set left matrix[right]).set right matrix[left] (by
        simpa) := by
  unfold Generated.StrictRecombine.swapQQRows at hrun
  split at hrun
  next hleft =>
    split at hrun
    next hright => exact ⟨hleft, hright, (Except.ok.inj hrun).symm⟩
    next hright => contradiction
  next hleft => contradiction

theorem swapQQRows_size
    (matrix output : Generated.StrictRecombine.QQMatrix) (left right : Nat)
    (hrun : Generated.StrictRecombine.swapQQRows matrix left right =
      .ok output) :
    output.size = matrix.size := by
  rcases swapQQRows_ok matrix output left right hrun with
    ⟨hleft, hright, rfl⟩
  simp

theorem swapQQRows_get
    (matrix output : Generated.StrictRecombine.QQMatrix)
    (left right row : Nat)
    (hrun : Generated.StrictRecombine.swapQQRows matrix left right =
      .ok output) (hrow : row < matrix.size) :
    output[row]! =
      if right = row then matrix[left]!
      else if left = row then matrix[right]!
      else matrix[row]! := by
  rcases swapQQRows_ok matrix output left right hrun with
    ⟨hleft, hright, rfl⟩
  rw [getElem!_pos _ row (by simp; exact hrow)]
  simp only [Array.getElem_set]
  rw [getElem!_pos matrix left hleft, getElem!_pos matrix right hright,
    getElem!_pos matrix row hrow]

def swapRowsArray {element : Type}
    (matrix : Array (Array element)) (left right : Nat) :
    Array (Array element) :=
  matrix.setIfInBounds left matrix[right]! |>.setIfInBounds right matrix[left]!

theorem swapRowsArray_size {element : Type}
    (matrix : Array (Array element)) (left right : Nat) :
    (swapRowsArray matrix left right).size = matrix.size := by
  simp [swapRowsArray]

theorem swapRowsArray_get {element : Type} [Inhabited element]
    (matrix : Array (Array element)) (left right row : Nat)
    (hleft : left < matrix.size) (hright : right < matrix.size)
    (hrow : row < matrix.size) :
    (swapRowsArray matrix left right)[row]! =
      if right = row then matrix[left]!
      else if left = row then matrix[right]!
      else matrix[row]! := by
  unfold swapRowsArray
  simp only [Array.setIfInBounds, dif_pos hleft]
  rw [dif_pos (show right < (matrix.set left matrix[right]! hleft).size by
    simpa using hright)]
  rw [getElem!_pos _ row (by simp; exact hrow), Array.getElem_set]
  rw [getElem!_pos matrix left hleft, getElem!_pos matrix right hright,
    getElem!_pos matrix row hrow]
  by_cases hrightRow : right = row
  · simp [hrightRow]
  · by_cases hleftRow : left = row
    · simp [hrightRow, hleftRow]
    · simp only [hrightRow, hleftRow, if_false]
      rw [Array.getElem_set]
      simp [hleftRow]

theorem swapQQRows_output_eq
    (matrix output : Generated.StrictRecombine.QQMatrix) (left right : Nat)
    (hrun : Generated.StrictRecombine.swapQQRows matrix left right =
      .ok output) :
    output = swapRowsArray matrix left right := by
  rcases swapQQRows_ok matrix output left right hrun with
    ⟨hleft, hright, rfl⟩
  simp [swapRowsArray, Array.setIfInBounds, hleft, hright]

def lovaszSwapCorrectedMu (mu : Generated.StrictRecombine.QQMatrix)
    (k : Nat) (muNew : QQ) : Generated.StrictRecombine.QQMatrix :=
  let swapped := swapRowsArray mu k (k - 1)
  swapped.setIfInBounds k
    (swapped[k]!.setIfInBounds (k - 1) muNew)

def lovaszSwapMuResult (mu : Generated.StrictRecombine.QQMatrix)
    (k : Nat) (muOld muNew : QQ) : Generated.StrictRecombine.QQMatrix :=
  updateMuAfterSwapArray (lovaszSwapCorrectedMu mu k muNew)
    k muOld muNew (k + 1)

theorem lovaszSwapMuResult_of_generated
    (mu swappedMu finalMu : Generated.StrictRecombine.QQMatrix)
    (k : Nat) (muOld muNew : QQ)
    (hswap : Generated.StrictRecombine.swapQQRows mu k (k - 1) =
      .ok swappedMu)
    (hkSwapped : k < swappedMu.size)
    (hpredSwapped : k - 1 < swappedMu[k].size)
    (hrowsK : ∀ row (hrow : row <
      (swappedMu.set k
        (swappedMu[k].set (k - 1) muNew)).size),
      k < (swappedMu.set k
        (swappedMu[k].set (k - 1) muNew))[row].size)
    (hrowsPred : ∀ row (hrow : row <
      (swappedMu.set k
        (swappedMu[k].set (k - 1) muNew)).size),
      k - 1 < (swappedMu.set k
        (swappedMu[k].set (k - 1) muNew))[row].size)
    (hloop : Generated.StrictRecombine.updateMuAfterSwapLoop
      (swappedMu.set k (swappedMu[k].set (k - 1) muNew))
      k muOld muNew (k + 1) = .ok finalMu) :
    finalMu = lovaszSwapMuResult mu k muOld muNew := by
  have hswapExact := swapQQRows_output_eq mu swappedMu k (k - 1) hswap
  have hloopExact := updateMuAfterSwapLoop_output_eq
    (swappedMu.set k (swappedMu[k].set (k - 1) muNew)) finalMu
    k muOld muNew (k + 1) hrowsK hrowsPred hloop
  rw [hloopExact]
  unfold lovaszSwapMuResult lovaszSwapCorrectedMu
  rw [← hswapExact]
  simp [Array.setIfInBounds, hkSwapped, hpredSwapped]

theorem lovaszSwapCorrectedMu_size
    (mu : Generated.StrictRecombine.QQMatrix) (k : Nat) (muNew : QQ) :
    (lovaszSwapCorrectedMu mu k muNew).size = mu.size := by
  simp [lovaszSwapCorrectedMu, swapRowsArray_size]

theorem lovaszSwapCorrectedMu_get
    (mu : Generated.StrictRecombine.QQMatrix) (k row : Nat) (muNew : QQ)
    (hkPositive : 0 < k) (hk : k < mu.size) (hrow : row < mu.size) :
    (lovaszSwapCorrectedMu mu k muNew)[row]! =
      if row = k then
        mu[k - 1]!.setIfInBounds (k - 1) muNew
      else if row = k - 1 then mu[k]!
      else mu[row]! := by
  have hkPred : k - 1 < mu.size := by omega
  have hkSwapped : k < (swapRowsArray mu k (k - 1)).size := by
    simpa [swapRowsArray_size] using hk
  have hrowSwapped : row < (swapRowsArray mu k (k - 1)).size := by
    simpa [swapRowsArray_size] using hrow
  unfold lovaszSwapCorrectedMu
  simp only [Array.setIfInBounds, dif_pos hkSwapped]
  rw [getElem!_pos _ row (by simpa [swapRowsArray_size] using hrow),
    Array.getElem_set]
  by_cases hrowK : row = k
  · subst row
    simp only [if_pos]
    rw [swapRowsArray_get mu k (k - 1) k hk hkPred hk]
    have hkNePred : k - 1 ≠ k := by omega
    simp [hkNePred]
  · rw [if_neg hrowK]
    rw [← getElem!_pos (swapRowsArray mu k (k - 1)) row hrowSwapped]
    rw [swapRowsArray_get mu k (k - 1) row hk hkPred hrow]
    by_cases hrowPred : row = k - 1
    · subst row
      have hkNePred : k - 1 ≠ k := by omega
      have hkPredNe : k ≠ k - 1 := by omega
      simp [hkNePred, hkPredNe]
    · have hpredNeRow : k - 1 ≠ row := Ne.symm hrowPred
      have hkNeRow : k ≠ row := Ne.symm hrowK
      simp [hrowPred, hpredNeRow, hkNeRow]

theorem lovaszSwapMuResult_get
    (mu : Generated.StrictRecombine.QQMatrix) (k row : Nat)
    (muOld muNew : QQ) (hkPositive : 0 < k) (hk : k < mu.size)
    (hrow : row < mu.size) :
    (lovaszSwapMuResult mu k muOld muNew)[row]! =
      if row = k then
        mu[k - 1]!.setIfInBounds (k - 1) muNew
      else if row = k - 1 then mu[k]!
      else if k < row then
        updateMuAfterSwapRow mu[row]! k muOld muNew
      else mu[row]! := by
  have hcorrectedRow : row < (lovaszSwapCorrectedMu mu k muNew).size := by
    simpa [lovaszSwapCorrectedMu_size] using hrow
  unfold lovaszSwapMuResult
  by_cases hafter : k < row
  · rw [updateMuAfterSwapArray_get_at
      (lovaszSwapCorrectedMu mu k muNew) k row muOld muNew (k + 1)
      hcorrectedRow (by omega)]
    have hcorrected := lovaszSwapCorrectedMu_get mu k row muNew
      hkPositive hk hrow
    have hrowK : row ≠ k := by omega
    have hrowPred : row ≠ k - 1 := by omega
    simp only [hrowK, hrowPred, if_false] at hcorrected ⊢
    rw [hcorrected]
    simp [hafter]
  · rw [updateMuAfterSwapArray_get_before
      (lovaszSwapCorrectedMu mu k muNew) k row muOld muNew (k + 1)
      hcorrectedRow (by omega)]
    rw [lovaszSwapCorrectedMu_get mu k row muNew hkPositive hk hrow]
    simp [hafter]

theorem lovaszSwapMuResult_size
    (mu : Generated.StrictRecombine.QQMatrix) (k : Nat)
    (muOld muNew : QQ) :
    (lovaszSwapMuResult mu k muOld muNew).size = mu.size := by
  unfold lovaszSwapMuResult
  rw [updateMuAfterSwapArray_size, lovaszSwapCorrectedMu_size]

theorem lovaszSwapMuResult_entry
    (mu : Generated.StrictRecombine.QQMatrix) (k row column : Nat)
    (muOld muNew : QQ) (hkPositive : 0 < k) (hk : k < mu.size)
    (hrow : row < mu.size) (hcolumn : column < mu.size)
    (hrowsSquare : ∀ index (hindex : index < mu.size),
      mu[index]!.size = mu.size) :
    ((lovaszSwapMuResult mu k muOld muNew)[row]!)[column]! =
      if row = k then
        if column = k - 1 then muNew else (mu[k - 1]!)[column]!
      else if row = k - 1 then (mu[k]!)[column]!
      else if k < row then
        if column = k - 1 then
          (mu[row]!)[k]! + muNew *
            ((mu[row]!)[k - 1]! - muOld * (mu[row]!)[k]!)
        else if column = k then
          (mu[row]!)[k - 1]! - muOld * (mu[row]!)[k]!
        else (mu[row]!)[column]!
      else (mu[row]!)[column]! := by
  have hkPred : k - 1 < mu.size := by omega
  have hresult := lovaszSwapMuResult_get mu k row muOld muNew
    hkPositive hk hrow
  by_cases hrowK : row = k
  · subst row
    simp only [if_pos] at hresult ⊢
    rw [hresult]
    have hpredInRow : k - 1 < mu[k - 1]!.size := by
      rw [hrowsSquare (k - 1) hkPred]
      exact hkPred
    have hcolumnInRow : column < mu[k - 1]!.size := by
      rw [hrowsSquare (k - 1) hkPred]
      exact hcolumn
    simp only [Array.setIfInBounds, dif_pos hpredInRow]
    rw [getElem!_pos _ column (by simpa using hcolumnInRow),
      Array.getElem_set]
    by_cases hcolumnPred : column = k - 1
    · simp [hcolumnPred]
    · rw [if_neg (Ne.symm hcolumnPred)]
      rw [if_neg hcolumnPred]
      exact (getElem!_pos mu[k - 1]! column hcolumnInRow).symm
  · simp only [hrowK, if_false] at hresult ⊢
    by_cases hrowPred : row = k - 1
    · simp [hrowPred] at hresult ⊢
      rw [hresult]
    · simp only [hrowPred, if_false] at hresult ⊢
      by_cases hafter : k < row
      · rw [if_pos hafter] at hresult ⊢
        rw [hresult]
        have hkInRow : k < mu[row]!.size := by
          rw [hrowsSquare row hrow]
          exact hk
        have hpredInRow : k - 1 < mu[row]!.size := by
          rw [hrowsSquare row hrow]
          exact hkPred
        have hcolumnInRow : column < mu[row]!.size := by
          rw [hrowsSquare row hrow]
          exact hcolumn
        exact updateMuAfterSwapRow_get mu[row]! k column muOld muNew
          hkInRow hpredInRow hkPositive hcolumnInRow
      · rw [if_neg hafter] at hresult ⊢
        rw [hresult]

/-- Exact mu array produced by one generated size-reduction coefficient. -/
def sizeReduceMuResult (mu : Generated.StrictRecombine.QQMatrix)
    (k source : Nat) (q : ZZ) : Generated.StrictRecombine.QQMatrix :=
  if q = 0 then mu
  else
    let muPrefix := reduceMuPrefixArray mu k source q source 0
    muPrefix.setIfInBounds k
      (muPrefix[k]!.setIfInBounds source
        ((muPrefix[k]!)[source]! - (q : QQ)))

theorem sizeReduceMuResult_target_get
    (mu : Generated.StrictRecombine.QQMatrix)
    (k source column : Nat) (q : ZZ)
    (hk : k < mu.size) (hsource : source < mu.size)
    (hsourceKRow : source < mu[k].size)
    (hlimitSource : source ≤ mu[source].size)
    (hcolumn : column < mu[k].size) (hne : k ≠ source) :
    ((sizeReduceMuResult mu k source q)[k]!)[column]! =
      if column < source then
        (mu[k]!)[column]! - (q : QQ) * (mu[source]!)[column]!
      else if column = source then
        (mu[k]!)[column]! - (q : QQ)
      else (mu[k]!)[column]! := by
  by_cases hq : q = 0
  · subst q
    simp [sizeReduceMuResult]
  · let muPrefix := reduceMuPrefixArray mu k source q source 0
    have hkPrefix : k < muPrefix.size := by
      simpa [muPrefix, reduceMuPrefixArray_size] using hk
    have hsourcePrefix : source < muPrefix[k]!.size := by
      have hsize := reduceMuPrefixArray_row_size mu k source k q source 0 hk
      simpa [muPrefix, getElem!_pos mu k hk] using
        (show source < (reduceMuPrefixArray mu k source q source 0)[k]!.size by
          rw [hsize, getElem!_pos mu k hk]
          exact hsourceKRow)
    have hcolumnPrefix : column < muPrefix[k]!.size := by
      have hsize := reduceMuPrefixArray_row_size mu k source k q source 0 hk
      simpa [muPrefix, getElem!_pos mu k hk] using
        (show column < (reduceMuPrefixArray mu k source q source 0)[k]!.size by
          rw [hsize, getElem!_pos mu k hk]
          exact hcolumn)
    have hsourcePrefix' : source < muPrefix[k].size := by
      simpa [getElem!_pos muPrefix k hkPrefix] using hsourcePrefix
    have hcolumnPrefix' : column < muPrefix[k].size := by
      simpa [getElem!_pos muPrefix k hkPrefix] using hcolumnPrefix
    unfold sizeReduceMuResult
    rw [if_neg hq]
    change (((muPrefix.setIfInBounds k
      (muPrefix[k]!.setIfInBounds source
        ((muPrefix[k]!)[source]! - (q : QQ))))[k]!)[column]!) = _
    simp only [Array.setIfInBounds, dif_pos hkPrefix,
      getElem!_pos muPrefix k hkPrefix, dif_pos hsourcePrefix']
    rw [getElem!_pos _ k (by simpa), Array.getElem_set]
    simp only [if_pos]
    rw [getElem!_pos _ column (by simpa), Array.getElem_set]
    by_cases hcolumnSource : column = source
    · subst column
      rw [if_pos rfl]
      have hprefixGet := reduceMuPrefixArray_target_get mu k source q
        source 0 source hk hsource (Nat.le_of_lt hsourceKRow)
        hlimitSource hsourceKRow hne
      simp only [Nat.zero_le, true_and, Nat.lt_irrefl, if_false] at hprefixGet
      have hprefixGet' : muPrefix[k][source]! = (mu[k]!)[source]! := by
        simpa [muPrefix, getElem!_pos muPrefix k hkPrefix] using hprefixGet
      rw [hprefixGet']
      simp
    · rw [if_neg hcolumnSource]
      have hprefixGet := reduceMuPrefixArray_target_get mu k source q
        source 0 column hk hsource (Nat.le_of_lt hsourceKRow)
        hlimitSource hcolumn hne
      simp only [Nat.zero_le, true_and] at hprefixGet
      rw [if_neg (Ne.symm hcolumnSource)]
      have hprefixGet' : muPrefix[k][column]! =
          if column < source then
            (mu[k]!)[column]! - (q : QQ) * (mu[source]!)[column]!
          else (mu[k]!)[column]! := by
        simpa [muPrefix, getElem!_pos muPrefix k hkPrefix] using hprefixGet
      have hprefixGet'' : muPrefix[k][column] =
          if column < source then
            (mu[k]!)[column]! - (q : QQ) * (mu[source]!)[column]!
          else (mu[k]!)[column]! := by
        simpa only [getElem!_pos muPrefix[k] column hcolumnPrefix'] using
          hprefixGet'
      rw [hprefixGet'']

theorem sizeReduceMuResult_get_of_ne
    (mu : Generated.StrictRecombine.QQMatrix)
    (k source row : Nat) (q : ZZ)
    (hrow : row < mu.size) (hrowK : row ≠ k) :
    (sizeReduceMuResult mu k source q)[row]! = mu[row]! := by
  by_cases hq : q = 0
  · subst q
    simp [sizeReduceMuResult]
  · let muPrefix := reduceMuPrefixArray mu k source q source 0
    have hprefixSize : muPrefix.size = mu.size := by
      exact reduceMuPrefixArray_size mu k source q source 0
    have hrowPrefix : row < muPrefix.size := by simpa [hprefixSize]
    have hprefixRow : muPrefix[row]! = mu[row]! := by
      simpa [muPrefix] using reduceMuPrefixArray_get_of_ne
        mu k source row q source 0 hrow hrowK
    unfold sizeReduceMuResult
    rw [if_neg hq]
    change (muPrefix.setIfInBounds k
      (muPrefix[k]!.setIfInBounds source
        ((muPrefix[k]!)[source]! - (q : QQ))))[row]! = mu[row]!
    by_cases hkPrefix : k < muPrefix.size
    · simp only [Array.setIfInBounds, dif_pos hkPrefix]
      rw [getElem!_pos _ row (by simpa [hprefixSize]), Array.getElem_set]
      rw [if_neg (Ne.symm hrowK)]
      exact (getElem!_pos muPrefix row hrowPrefix).symm.trans hprefixRow
    · simp [Array.setIfInBounds, hkPrefix, hprefixRow]

theorem sizeReduceMuResult_size
    (mu : Generated.StrictRecombine.QQMatrix)
    (k source : Nat) (q : ZZ) :
    (sizeReduceMuResult mu k source q).size = mu.size := by
  by_cases hq : q = 0
  · subst q
    simp [sizeReduceMuResult]
  · unfold sizeReduceMuResult
    rw [if_neg hq]
    simp [reduceMuPrefixArray_size]

theorem sizeReduceMuResult_row_size
    (mu : Generated.StrictRecombine.QQMatrix)
    (k source row : Nat) (q : ZZ)
    (hrow : row < mu.size) :
    (sizeReduceMuResult mu k source q)[row]!.size = mu[row]!.size := by
  by_cases hrowK : row = k
  · subst row
    by_cases hq : q = 0
    · subst q
      simp [sizeReduceMuResult]
    · let muPrefix := reduceMuPrefixArray mu k source q source 0
      have hkPrefix : k < muPrefix.size := by
        simpa [muPrefix, reduceMuPrefixArray_size] using hrow
      have hprefixRowSize := reduceMuPrefixArray_row_size
        mu k source k q source 0 hrow
      unfold sizeReduceMuResult
      rw [if_neg hq]
      change (muPrefix.setIfInBounds k
        (muPrefix[k]!.setIfInBounds source
          ((muPrefix[k]!)[source]! - (q : QQ))))[k]!.size = mu[k]!.size
      simp only [Array.setIfInBounds, dif_pos hkPrefix]
      rw [getElem!_pos _ k (by simpa using hkPrefix), Array.getElem_set]
      simp only [if_pos]
      split
      next hsource =>
        simp only [Array.size_set]
        simpa [muPrefix] using hprefixRowSize
      next hsource => simpa [muPrefix] using hprefixRowSize
  · rw [sizeReduceMuResult_get_of_ne mu k source row q hrow hrowK]

theorem sizeReduceAt_preserves_norms_k
    (state output : Generated.StrictRecombine.LLLState) (j : Nat)
    (hrun : Generated.StrictRecombine.sizeReduceAt state j = .ok output) :
    output.norms = state.norms ∧ output.k = state.k := by
  unfold Generated.StrictRecombine.sizeReduceAt at hrun
  split at hrun
  next hk =>
    split at hrun
    next hj =>
      cases hround : Generated.StrictRecombine.roundQQ state.mu[state.k][j] with
      | error fault => simp [hround] at hrun
      | ok q =>
        simp only [hround] at hrun
        split at hrun
        next hzero =>
          have hout := Except.ok.inj hrun
          subst output
          exact ⟨rfl, rfl⟩
        next hnonzero =>
          cases hsubtract : Generated.StrictRecombine.subtractMatrixRows
              state.matrix state.transform state.k j q with
          | error fault => simp [hsubtract] at hrun
          | ok matrices =>
            rcases matrices with ⟨matrix', transform'⟩
            simp only [hsubtract] at hrun
            cases hmu : Generated.StrictRecombine.reduceMuPrefixLoop
                state.mu state.k j q j 0 with
            | error fault => simp [hmu] at hrun
            | ok mu' =>
              simp only [hmu] at hrun
              split at hrun
              next hk' =>
                split at hrun
                next hj' =>
                  have hout := Except.ok.inj hrun
                  subst output
                  exact ⟨rfl, rfl⟩
                next hj' => contradiction
              next hk' => contradiction
    next hj => contradiction
  next hk => contradiction

theorem extraSizeReduceLoop_preserves_norms_k
    (remaining : Nat) (state output : Generated.StrictRecombine.LLLState)
    (hrun : Generated.StrictRecombine.extraSizeReduceLoop remaining state =
      .ok output) :
    output.norms = state.norms ∧ output.k = state.k := by
  induction remaining generalizing state output with
  | zero =>
      simp [Generated.StrictRecombine.extraSizeReduceLoop] at hrun
      subst output
      exact ⟨rfl, rfl⟩
  | succ remaining ih =>
      rw [Generated.StrictRecombine.extraSizeReduceLoop] at hrun
      cases hstep : Generated.StrictRecombine.sizeReduceAt state remaining with
      | error fault => simp [hstep] at hrun
      | ok next =>
        simp only [hstep] at hrun
        have hhead := sizeReduceAt_preserves_norms_k state next remaining hstep
        have htail := ih next output hrun
        exact ⟨htail.1.trans hhead.1, htail.2.trans hhead.2⟩

theorem lllStep_advanced_control
    (state output : Generated.StrictRecombine.LLLState)
    (hrun : Generated.StrictRecombine.lllStep state =
      .ok (.advanced output)) :
    output.norms = state.norms ∧ output.k = state.k + 1 := by
  rw [Generated.StrictRecombine.lllStep] at hrun
  by_cases hkPositive : 0 < state.k
  · rw [dif_pos hkPositive] at hrun
    by_cases hkMatrix : state.k < state.matrix.size
    · rw [dif_pos hkMatrix] at hrun
      cases hreduce : Generated.StrictRecombine.sizeReduceAt state
          (state.k - 1) with
      | error fault => simp [hreduce] at hrun
      | ok reduced =>
        simp only [hreduce] at hrun
        have hreduceControl := sizeReduceAt_preserves_norms_k state reduced
          (state.k - 1) hreduce
        split at hrun
        next hkNorm =>
          split at hrun
          next hpredNorm =>
            split at hrun
            next hkMu =>
              split at hrun
              next hpredMu =>
                dsimp at hrun
                split at hrun
                next hlovasz =>
                  cases hextra : Generated.StrictRecombine.extraSizeReduceLoop
                      (reduced.k - 1) reduced with
                  | error fault => simp [hextra] at hrun
                  | ok fullyReduced =>
                    simp only [hextra] at hrun
                    have hextraControl :=
                      extraSizeReduceLoop_preserves_norms_k
                        (reduced.k - 1) reduced fullyReduced hextra
                    have hout := Except.ok.inj hrun
                    injection hout with hstate
                    subst output
                    constructor
                    · exact hextraControl.1.trans hreduceControl.1
                    · simp only
                      rw [hextraControl.2, hreduceControl.2]
                next hlovasz =>
                  repeat' first | split at hrun | simp_all
              next hpredMu => contradiction
            next hkMu => contradiction
          next hpredNorm => contradiction
        next hkNorm => contradiction
    · rw [dif_neg hkMatrix] at hrun
      simp [hkMatrix] at hrun
  · rw [dif_neg hkPositive] at hrun
    simp [hkPositive] at hrun

theorem lllStep_swapped_k
    (state output : Generated.StrictRecombine.LLLState)
    (hrun : Generated.StrictRecombine.lllStep state =
      .ok (.swapped output)) :
    output.k = Nat.max (state.k - 1) 1 := by
  rw [Generated.StrictRecombine.lllStep] at hrun
  by_cases hkPositive : 0 < state.k
  · rw [dif_pos hkPositive] at hrun
    by_cases hkMatrix : state.k < state.matrix.size
    · rw [dif_pos hkMatrix] at hrun
      cases hreduce : Generated.StrictRecombine.sizeReduceAt state
          (state.k - 1) with
      | error fault => simp [hreduce] at hrun
      | ok reduced =>
        simp only [hreduce] at hrun
        have hcontrol := sizeReduceAt_preserves_norms_k state reduced
          (state.k - 1) hreduce
        repeat' first | split at hrun | simp_all
        all_goals cases hrun; rfl
    · rw [dif_neg hkMatrix] at hrun
      contradiction
  · rw [dif_neg hkPositive] at hrun
    contradiction

theorem lllStep_swapped_norms
    (state output : Generated.StrictRecombine.LLLState)
    (hrun : Generated.StrictRecombine.lllStep state =
      .ok (.swapped output)) :
    ∃ reduced,
      Generated.StrictRecombine.sizeReduceAt state (state.k - 1) = .ok reduced ∧
      reduced.norms = state.norms ∧ reduced.k = state.k ∧
      output.norms = Generated.StrictRecombine.normsAfterLovaszSwap
        reduced.norms reduced.k
          ((reduced.mu[reduced.k]!)[reduced.k - 1]!) := by
  rw [Generated.StrictRecombine.lllStep] at hrun
  by_cases hkPositive : 0 < state.k
  · rw [dif_pos hkPositive] at hrun
    by_cases hkMatrix : state.k < state.matrix.size
    · rw [dif_pos hkMatrix] at hrun
      cases hreduce : Generated.StrictRecombine.sizeReduceAt state
          (state.k - 1) with
      | error fault => simp [hreduce] at hrun
      | ok reduced =>
        simp only [hreduce] at hrun
        have hcontrol := sizeReduceAt_preserves_norms_k state reduced
          (state.k - 1) hreduce
        repeat' first | split at hrun | simp_all
        all_goals cases hrun
        all_goals
          have hkMuState : state.k < reduced.mu.size := by
            rw [← hcontrol.2]
            assumption
          have hpredMuState : state.k - 1 < reduced.mu[state.k].size := by
            simpa only [hcontrol.2] using
              (‹reduced.k - 1 < reduced.mu[reduced.k].size›)
          change Generated.StrictRecombine.normsAfterLovaszSwap state.norms
              state.k ((reduced.mu[state.k]'hkMuState)[state.k - 1]'hpredMuState) =
            Generated.StrictRecombine.normsAfterLovaszSwap state.norms
              state.k ((reduced.mu[state.k]!)[state.k - 1]!)
          congr 1
          rw [getElem!_pos reduced.mu state.k hkMuState]
          rw [getElem!_pos (reduced.mu[state.k]) (state.k - 1) hpredMuState]
    · rw [dif_neg hkMatrix] at hrun
      simp [hkMatrix] at hrun
  · rw [dif_neg hkPositive] at hrun
    simp [hkPositive] at hrun

/-- The classical LLL potential, evaluated on the exact Gram–Schmidt norm
array carried by the generated C++ state. -/
noncomputable def arrayLLLPotential (norms : Array QQ) : QQ :=
  ∏ i : Fin norms.size, norms[i] ^ (norms.size - i.val)

theorem arrayLLLPotential_pos (norms : Array QQ)
    (hpositive : ∀ i (hi : i < norms.size), 0 < norms[i]'hi) :
    0 < arrayLLLPotential norms := by
  unfold arrayLLLPotential
  apply Finset.prod_pos
  intro i _
  exact pow_pos (hpositive i i.isLt) _

/-- View an operational array state as the finite-function state used by the
pure L2 LLL potential theorem.  Matrix and μ reads are total only because this
view is used for the termination algebra; the execution theorem separately
establishes every source-level bound before a read. -/
noncomputable def toPotentialState
    (state : Generated.StrictRecombine.LLLState) :
    LLLState state.norms.size state.matrix.size where
  basis := fun i j => (state.matrix[i.val]!)[j.val]!
  mu := fun i j => (state.mu[i.val]!)[j.val]!
  gs_norm_sq := fun i => state.norms[i]
  k := state.k

theorem toPotentialState_potential
    (state : Generated.StrictRecombine.LLLState) :
    lllPotential' (toPotentialState state) = arrayLLLPotential state.norms := by
  rfl

theorem normsAfterLovaszSwap_size (norms : Array QQ) (k : Nat) (mu : QQ) :
    (Generated.StrictRecombine.normsAfterLovaszSwap norms k mu).size =
      norms.size := by
  unfold Generated.StrictRecombine.normsAfterLovaszSwap
  dsimp
  split <;> simp

theorem normsAfterLovaszSwap_get (norms : Array QQ) (k : Nat) (mu : QQ)
    (hk : k < norms.size) (hpred : k - 1 < norms.size)
    (hnew : norms[k] + mu * mu * norms[k - 1] ≠ 0)
    (i : Nat) (hi : i < norms.size) :
    (Generated.StrictRecombine.normsAfterLovaszSwap norms k mu)[i]! =
      if k - 1 = i then norms[k] + mu * mu * norms[k - 1]
      else if k = i then
        norms[k] * norms[k - 1] /
          (norms[k] + mu * mu * norms[k - 1])
      else norms[i] := by
  unfold Generated.StrictRecombine.normsAfterLovaszSwap
  dsimp
  rw [getElem!_pos norms k hk, getElem!_pos norms (k - 1) hpred]
  rw [if_pos hnew]
  rw [getElem!_pos _ i (by simpa using hi)]
  rw [Array.getElem_setIfInBounds (by simpa using hi)]
  rw [Array.getElem_setIfInBounds hi]

theorem normsAfterLovaszSwap_positive (norms : Array QQ) (k : Nat) (mu : QQ)
    (hk : k < norms.size) (hpred : k - 1 < norms.size)
    (hpositiveK : 0 < k)
    (hpositive : ∀ i (hi : i < norms.size), 0 < norms[i]) :
    ∀ i (hi : i <
      (Generated.StrictRecombine.normsAfterLovaszSwap norms k mu).size),
      0 < (Generated.StrictRecombine.normsAfterLovaszSwap norms k mu)[i] := by
  intro i hi
  have hiOld : i < norms.size := by
    simpa [normsAfterLovaszSwap_size] using hi
  have hkPos := hpositive k hk
  have hpredPos := hpositive (k - 1) hpred
  have hnewPos : 0 < norms[k] + mu * mu * norms[k - 1] := by
    nlinarith [sq_nonneg mu]
  have hvalue :
      (Generated.StrictRecombine.normsAfterLovaszSwap norms k mu)[i] =
        (if k - 1 = i then norms[k] + mu * mu * norms[k - 1]
        else if k = i then
          norms[k] * norms[k - 1] /
            (norms[k] + mu * mu * norms[k - 1])
        else norms[i]) := by
    have hraw := normsAfterLovaszSwap_get norms k mu hk hpred
      (ne_of_gt hnewPos) i hiOld
    rw [getElem!_pos
      (Generated.StrictRecombine.normsAfterLovaszSwap norms k mu) i hi] at hraw
    exact hraw
  rw [hvalue]
  split
  next hpredI => exact hnewPos
  next hpredI =>
    split
    next hkI => positivity
    next hkI => exact hpositive i hiOld

theorem normsAfterLovaszSwap_pair_product (norms : Array QQ)
    (k : Nat) (mu : QQ)
    (hk : k < norms.size) (hpred : k - 1 < norms.size)
    (hpositiveK : 0 < k)
    (hpositive : ∀ i (hi : i < norms.size), 0 < norms[i]) :
    (Generated.StrictRecombine.normsAfterLovaszSwap norms k mu)[k - 1]! *
        (Generated.StrictRecombine.normsAfterLovaszSwap norms k mu)[k]! =
      norms[k - 1] * norms[k] := by
  have hnewPos : 0 < norms[k] + mu * mu * norms[k - 1] := by
    have hkPos := hpositive k hk
    have hpredPos := hpositive (k - 1) hpred
    nlinarith [sq_nonneg mu]
  rw [normsAfterLovaszSwap_get norms k mu hk hpred
      (ne_of_gt hnewPos) (k - 1) hpred,
    normsAfterLovaszSwap_get norms k mu hk hpred
      (ne_of_gt hnewPos) k hk]
  have hkNe : k ≠ k - 1 := by omega
  have hpredNe : k - 1 ≠ k := hkNe.symm
  simp only [if_pos, hpredNe, if_false, hkNe]
  have hdenom : norms[k] + norms[k - 1] * mu ^ 2 ≠ 0 := by
    nlinarith [sq_nonneg mu]
  apply (mul_left_cancel₀ hdenom)
  field_simp [hdenom]

/-- The exact `2 × 2` rational identity behind the generated adjacent
Lovász swap.  These are the four entries of `S · diag(a,b) · Sᵀ`, where
`S = [[mu, 1], [1 - muNew*mu, -muNew]]`. -/
theorem lovaszLocalLDLAlgebra (a b mu muNew delta : QQ)
    (hdelta : delta = b + mu * mu * a)
    (hdeltaNe : delta ≠ 0)
    (hmuNew : muNew = mu * a / delta) :
    mu * mu * a + b = delta ∧
      mu * a * (1 - muNew * mu) - b * muNew = 0 ∧
      (1 - muNew * mu) * a * mu - muNew * b = 0 ∧
      (1 - muNew * mu) * (1 - muNew * mu) * a +
          muNew * muNew * b = a * b / delta := by
  constructor
  · rw [hdelta]
    ring
  constructor
  · rw [hmuNew]
    field_simp [hdeltaNe]
    rw [hdelta]
    ring
  constructor
  · rw [hmuNew]
    field_simp [hdeltaNe]
    rw [hdelta]
    ring
  · rw [hmuNew]
    field_simp [hdeltaNe]
    rw [hdelta]
    ring

/-- Once the determinant potential is converted to a natural number, this is
the concrete scalar rank used by the generated while-loop: an advancing step
spends the low-order progress component, while a swap must lower the dominant
lattice component. -/
def lllLexRank (latticePotential dimension index : Nat) : Nat :=
  latticePotential * (dimension + 1) + (dimension - index)

theorem lllLexRank_advanced (potential dimension index : Nat)
    (hindex : index < dimension) :
    lllLexRank potential dimension (index + 1) <
      lllLexRank potential dimension index := by
  unfold lllLexRank
  omega

theorem lllLexRank_swap (oldPotential newPotential dimension oldIndex newIndex : Nat)
    (hpotential : newPotential < oldPotential) :
    lllLexRank newPotential dimension newIndex <
      lllLexRank oldPotential dimension oldIndex := by
  unfold lllLexRank
  calc
    newPotential * (dimension + 1) + (dimension - newIndex) ≤
        newPotential * (dimension + 1) + dimension :=
      Nat.add_le_add_left (Nat.sub_le dimension newIndex) _
    _ < newPotential * (dimension + 1) + (dimension + 1) := by omega
    _ = (newPotential + 1) * (dimension + 1) := by
      rw [Nat.add_mul, one_mul]
    _ ≤ oldPotential * (dimension + 1) :=
      Nat.mul_le_mul_right _ (Nat.succ_le_iff.mpr hpotential)
    _ ≤ oldPotential * (dimension + 1) + (dimension - oldIndex) :=
      Nat.le_add_right _ _

/-- Integer Gram matrix of the first `prefix` rows of the current C++ lattice
basis.  Total reads make the definition available on every raw state; the
operational validity invariant later states that the basis is square. -/
noncomputable def gramPrefixMatrix
    (matrix : Generated.StrictRecombine.LLLMatrix) (rowCount : Nat) :
    Matrix (Fin rowCount) (Fin rowCount) Int :=
  fun i j => ∑ c : Fin matrix.size,
    (matrix[i.val]!)[c.val]! * (matrix[j.val]!)[c.val]!

noncomputable def basisPrefixMatrix
    (matrix : Generated.StrictRecombine.LLLMatrix)
    (rowCount columnCount : Nat) : Matrix (Fin rowCount) (Fin columnCount) Int :=
  fun i j => (matrix[i.val]!)[j.val]!

/-- Appending a zero coordinate to every old row and then a final row with
last coordinate one preserves the determinant.  This is exactly the block
lower-triangular shape of one generated CLD lattice-column extension. -/
theorem det_append_unit_lower {R : Type*} [CommRing R] (n : Nat)
    (A : Matrix (Fin n) (Fin n) R) (c : Fin n → R)
    (E : Matrix (Fin (n + 1)) (Fin (n + 1)) R)
    (hold : ∀ i j : Fin n, E i.castSucc j.castSucc = A i j)
    (hzero : ∀ i : Fin n, E i.castSucc (Fin.last n) = 0)
    (hlower : ∀ j : Fin n, E (Fin.last n) j.castSucc = c j)
    (hone : E (Fin.last n) (Fin.last n) = 1) :
    Matrix.det E = Matrix.det A := by
  let block : Matrix (Fin n ⊕ Fin 1) (Fin n ⊕ Fin 1) R :=
    Matrix.fromBlocks A 0 (fun _ j => c j) 1
  have hextension :
      E = block.submatrix finSumFinEquiv.symm finSumFinEquiv.symm := by
    ext i j
    refine Fin.lastCases ?_ (fun i => ?_) i <;>
      refine Fin.lastCases ?_ (fun j => ?_) j
    · simp [block, hone]
    · simp [block, hlower]
    · simp [block, hzero]
    · simp [block, hold]
  rw [hextension, Matrix.det_submatrix_equiv_self]
  simp [block]

theorem gramPrefixMatrix_eq_mul_transpose
    (matrix : Generated.StrictRecombine.LLLMatrix) (rowCount : Nat) :
    gramPrefixMatrix matrix rowCount =
      basisPrefixMatrix matrix rowCount matrix.size *
        (basisPrefixMatrix matrix rowCount matrix.size).transpose := by
  rfl

/-- Each prefix Gram determinant is a nonnegative integer.  Their product is
the discrete version of the rational LLL potential
`∏ᵢ Bᵢ^(n-i)`: Gram–Schmidt identifies the determinant of prefix `r`
with `∏ i<r, Bᵢ`. -/
noncomputable def lllDeterminantPotential
    (matrix : Generated.StrictRecombine.LLLMatrix) : Nat :=
  ∏ i : Fin matrix.size, (Matrix.det (gramPrefixMatrix matrix (i.val + 1))).natAbs

noncomputable def concreteLLLRank
    (state : Generated.StrictRecombine.LLLState) : Nat :=
  lllLexRank (lllDeterminantPotential state.matrix) state.matrix.size state.k

noncomputable def prefixNormProduct (norms : Array QQ) (rowCount : Nat) : QQ :=
  ∏ i : Fin rowCount, norms[i.val]!

/-- Unit lower-triangular Gram–Schmidt coefficient matrix carried by the
generated C++ arrays.  Reads are total here; the execution invariant below
will separately establish the source bounds that make them meaningful. -/
noncomputable def gsLowerPrefix
    (state : Generated.StrictRecombine.LLLState) (rowCount : Nat) :
    Matrix (Fin rowCount) (Fin rowCount) QQ :=
  fun i j =>
    if j.val < i.val then (state.mu[i.val]!)[j.val]!
    else if i = j then 1 else 0

noncomputable def gsNormDiagonal
    (state : Generated.StrictRecombine.LLLState) (rowCount : Nat) :
    Matrix (Fin rowCount) (Fin rowCount) QQ :=
  Matrix.diagonal fun i => state.norms[i.val]!

/-- Column basis change used by an adjacent Lovász swap.  It is the identity
outside the two affected Gram--Schmidt columns. -/
noncomputable def lovaszLocalTransform {dimension : Nat}
    (previous current : Fin dimension) (mu muNew : QQ) :
    Matrix (Fin dimension) (Fin dimension) QQ :=
  fun row column =>
    if row = previous then
      if column = previous then mu else if column = current then 1 else 0
    else if row = current then
      if column = previous then 1 - muNew * mu
      else if column = current then -muNew else 0
    else if row = column then 1 else 0

theorem finSumIteTwo {dimension : Nat} (previous current : Fin dimension)
    (hne : previous ≠ current) (left right : Fin dimension → QQ) :
    (∑ index, if index = previous then left index
      else if index = current then right index else 0) =
      left previous + right current := by
  calc
    (∑ index, if index = previous then left index
        else if index = current then right index else 0) =
        ∑ index, ((if index = previous then left index else 0) +
          (if index = current then right index else 0)) := by
      apply Finset.sum_congr rfl
      intro index _
      by_cases hprevious : index = previous
      · subst index
        simp [hne]
      · simp [hprevious]
    _ = (∑ index, if index = previous then left index else 0) +
          ∑ index, if index = current then right index else 0 := by
      rw [Finset.sum_add_distrib]
    _ = left previous + right current := by simp

theorem finSumIteTwoRedundantLeft {dimension : Nat}
    (previous current : Fin dimension) (hne : previous ≠ current)
    (left right redundant : Fin dimension → QQ) :
    (∑ index, if index = previous then left index
      else if index = current then right index
      else if index = previous then redundant index else 0) =
      left previous + right current := by
  rw [← finSumIteTwo previous current hne left right]
  apply Finset.sum_congr rfl
  intro index _
  by_cases hprevious : index = previous <;> simp [hprevious]

theorem finSumIteTwoRedundantRight {dimension : Nat}
    (previous current : Fin dimension) (hne : previous ≠ current)
    (left right redundant : Fin dimension → QQ) :
    (∑ index, if index = previous then left index
      else if index = current then right index
      else if index = current then redundant index else 0) =
      left previous + right current := by
  rw [← finSumIteTwo previous current hne left right]
  apply Finset.sum_congr rfl
  intro index _
  by_cases hprevious : index = previous
  · simp [hprevious]
  · by_cases hcurrent : index = current <;> simp [hprevious, hcurrent]

theorem finSumIteThree {dimension : Nat}
    (first second third : Fin dimension)
    (hfirstSecond : first ≠ second) (hthirdFirst : third ≠ first)
    (hthirdSecond : third ≠ second)
    (one two three : Fin dimension → QQ) :
    (∑ index, if index = first then one index
      else if index = second then two index
      else if index = third then three index else 0) =
      one first + two second + three third := by
  calc
    _ = ∑ index, ((if index = first then one index else 0) +
        (if index = second then two index else 0) +
        (if index = third then three index else 0)) := by
      apply Finset.sum_congr rfl
      intro index _
      by_cases hfirst : index = first
      · subst index
        simp [hfirstSecond, hthirdFirst, Ne.symm hthirdFirst]
      · by_cases hsecond : index = second
        · subst index
          simp [hfirst, hthirdSecond, Ne.symm hthirdSecond]
        · simp [hfirst, hsecond]
    _ = _ := by simp [Finset.sum_add_distrib]

theorem mul_lovaszLocalTransform_apply {dimension : Nat}
    (matrix : Matrix (Fin dimension) (Fin dimension) QQ)
    (previous current : Fin dimension) (mu muNew : QQ)
    (hne : previous ≠ current) (row column : Fin dimension) :
    (matrix * lovaszLocalTransform previous current mu muNew) row column =
      if column = previous then
        matrix row previous * mu +
          matrix row current * (1 - muNew * mu)
      else if column = current then
        matrix row previous - matrix row current * muNew
      else matrix row column := by
  rw [Matrix.mul_apply]
  by_cases hcolumnPrevious : column = previous
  · subst column
    simp only [lovaszLocalTransform, if_pos, Matrix.transpose_apply]
    simp only [mul_ite, mul_zero]
    rw [finSumIteTwoRedundantLeft previous current hne]
  · by_cases hcolumnCurrent : column = current
    · subst column
      simp only [lovaszLocalTransform, if_pos, if_neg hcolumnPrevious,
        Matrix.transpose_apply]
      simp only [mul_ite, mul_zero]
      rw [finSumIteTwoRedundantRight previous current hne]
      ring
    · simp only [lovaszLocalTransform, if_neg hcolumnPrevious,
        if_neg hcolumnCurrent, Matrix.transpose_apply]
      simp only [mul_ite, mul_zero]
      rw [finSumIteThree previous current column hne
        hcolumnPrevious hcolumnCurrent]
      ring

noncomputable def rowReindexMatrix {dimension : Nat}
    (equivalence : Fin dimension ≃ Fin dimension)
    (matrix : Matrix (Fin dimension) (Fin dimension) QQ) :
    Matrix (Fin dimension) (Fin dimension) QQ :=
  fun row column => matrix (equivalence.symm row) column

theorem reindex_ldl_eq_rowReindex {dimension : Nat}
    (equivalence : Fin dimension ≃ Fin dimension)
    (lower diagonal : Matrix (Fin dimension) (Fin dimension) QQ) :
    Matrix.reindex equivalence equivalence
        (lower * diagonal * lower.transpose) =
      rowReindexMatrix equivalence lower * diagonal *
        (rowReindexMatrix equivalence lower).transpose := by
  funext row column
  simp only [rowReindexMatrix, Matrix.reindex_apply, Matrix.submatrix_apply,
    Equiv.symm_symm, Matrix.mul_apply, Matrix.transpose_apply]

theorem lovaszLocalTransform_diagonal
    {dimension : Nat} (previous current : Fin dimension)
    (diagonal : Fin dimension → QQ) (a b mu muNew delta : QQ)
    (hne : previous ≠ current) (ha : diagonal previous = a)
    (hb : diagonal current = b)
    (hdelta : delta = b + mu * mu * a) (hdeltaNe : delta ≠ 0)
    (hmuNew : muNew = mu * a / delta) :
    lovaszLocalTransform previous current mu muNew *
          Matrix.diagonal diagonal *
          (lovaszLocalTransform previous current mu muNew).transpose =
      Matrix.diagonal (fun index =>
        if index = previous then delta
        else if index = current then a * b / delta
        else diagonal index) := by
  have halgebra := lovaszLocalLDLAlgebra a b mu muNew delta
    hdelta hdeltaNe hmuNew
  funext row column
  simp only [Matrix.mul_apply, Matrix.mul_diagonal, Matrix.transpose_apply,
    Matrix.diagonal_apply]
  by_cases hrowPrevious : row = previous
  · subst row
    by_cases hcolumnPrevious : column = previous
    · subst column
      simp [lovaszLocalTransform, hne, hne.symm, ha, hb, halgebra.1,
        finSumIteTwo]
      convert halgebra.1 using 1 <;> ring
    · by_cases hcolumnCurrent : column = current
      · subst column
        simp [lovaszLocalTransform, hne, hne.symm, ha, hb, halgebra.2.1,
          finSumIteTwo]
        convert halgebra.2.1 using 1 <;> ring
      · simp [lovaszLocalTransform, hne, hne.symm, hcolumnPrevious,
          Ne.symm hcolumnPrevious, hcolumnCurrent, Ne.symm hcolumnCurrent,
          finSumIteTwo]
  · by_cases hrowCurrent : row = current
    · subst row
      by_cases hcolumnPrevious : column = previous
      · subst column
        simp [lovaszLocalTransform, hne, hne.symm, ha, hb, halgebra.2.2.1,
          finSumIteTwo]
        convert halgebra.2.2.1 using 1 <;> ring
      · by_cases hcolumnCurrent : column = current
        · subst column
          simp [lovaszLocalTransform, hne, hne.symm, ha, hb, halgebra.2.2.2,
            finSumIteTwo]
          convert halgebra.2.2.2 using 1 <;> ring
        · simp [lovaszLocalTransform, hne, hne.symm, hcolumnPrevious,
            Ne.symm hcolumnPrevious, hcolumnCurrent, Ne.symm hcolumnCurrent,
            finSumIteTwo]
    · by_cases hcolumnPrevious : column = previous
      · subst column
        simp [lovaszLocalTransform, hne, hne.symm, hrowPrevious,
          Ne.symm hrowPrevious, hrowCurrent, Ne.symm hrowCurrent,
          finSumIteTwo]
      · by_cases hcolumnCurrent : column = current
        · subst column
          simp [lovaszLocalTransform, hne, hne.symm, hrowPrevious,
            Ne.symm hrowPrevious, hrowCurrent, Ne.symm hrowCurrent,
            finSumIteTwo]
        · by_cases hrowColumn : row = column
          · subst column
            simp [lovaszLocalTransform, hne, hne.symm, hrowPrevious,
              Ne.symm hrowPrevious, hrowCurrent, Ne.symm hrowCurrent,
              hcolumnPrevious, Ne.symm hcolumnPrevious, hcolumnCurrent,
              Ne.symm hcolumnCurrent, finSumIteTwo]
          · simp [lovaszLocalTransform, hne, hne.symm, hrowPrevious,
              Ne.symm hrowPrevious, hrowCurrent, Ne.symm hrowCurrent,
              hcolumnPrevious, Ne.symm hcolumnPrevious, hcolumnCurrent,
              Ne.symm hcolumnCurrent, hrowColumn, Ne.symm hrowColumn,
              finSumIteTwo]

noncomputable def gramPrefixMatrixQQ
    (matrix : Generated.StrictRecombine.LLLMatrix) (rowCount : Nat) :
    Matrix (Fin rowCount) (Fin rowCount) QQ :=
  (gramPrefixMatrix matrix rowCount).map fun value : Int => (value : QQ)

/-- The first `rowCount` integer basis rows, cast to rationals but retaining
all source columns.  This is the rectangular matrix whose Gram matrix is the
prefix used by generated Gram--Schmidt. -/
noncomputable def basisPrefixMatrixQQ
    (matrix : Generated.StrictRecombine.LLLMatrix) (rowCount : Nat) :
    Matrix (Fin rowCount) (Fin matrix.size) QQ :=
  (basisPrefixMatrix matrix rowCount matrix.size).map
    fun value : Int => (value : QQ)

theorem gramPrefixMatrixQQ_eq_mul_transpose
    (matrix : Generated.StrictRecombine.LLLMatrix) (rowCount : Nat) :
    gramPrefixMatrixQQ matrix rowCount =
      basisPrefixMatrixQQ matrix rowCount *
        (basisPrefixMatrixQQ matrix rowCount).transpose := by
  funext row column
  simp [gramPrefixMatrixQQ, gramPrefixMatrix, basisPrefixMatrixQQ,
    basisPrefixMatrix, Matrix.mul_apply, Matrix.transpose_apply,
    Int.cast_sum]

/-- Exact, checkable `G = L D Lᵀ` meaning of the generated `matrix`, `mu`,
and `norms` arrays.  This is execution state, not an existence witness and not
a semantic result oracle. -/
def ConcreteGramSchmidt
    (state : Generated.StrictRecombine.LLLState) : Prop :=
  ∀ rowCount, rowCount ≤ state.matrix.size →
    gramPrefixMatrixQQ state.matrix rowCount =
      gsLowerPrefix state rowCount * gsNormDiagonal state rowCount *
        (gsLowerPrefix state rowCount).transpose

/-- Exact generated Gram--Schmidt semantics restricted to the rows already
processed by the source initialization loop. -/
def ConcreteGramSchmidtUpTo
    (state : Generated.StrictRecombine.LLLState) (processed : Nat) : Prop :=
  ∀ rowCount, rowCount ≤ processed → rowCount ≤ state.matrix.size →
    gramPrefixMatrixQQ state.matrix rowCount =
      gsLowerPrefix state rowCount * gsNormDiagonal state rowCount *
        (gsLowerPrefix state rowCount).transpose

/-- The induction invariant at the head of the generated outer initialization
loop.  All data are the actual mutable arrays returned by the preceding source
iterations. -/
structure ProcessedGramSchmidtValid
    (state : Generated.StrictRecombine.LLLState) (processed : Nat) : Prop where
  shape : GramStorageShape state.mu state.norms state.matrix.size
  processed_le : processed ≤ state.matrix.size
  gram_schmidt : ConcreteGramSchmidtUpTo state processed
  norms_positive : ∀ index, index < processed → 0 < state.norms[index]!

theorem concreteGramSchmidtUpTo_zero
    (state : Generated.StrictRecombine.LLLState) :
    ConcreteGramSchmidtUpTo state 0 := by
  intro rowCount hprocessed hmatrix
  have hzero : rowCount = 0 := by omega
  subst rowCount
  ext row
  exact Fin.elim0 row

theorem concreteGramSchmidtUpTo_one_initial
    (matrix : Generated.StrictRecombine.LLLMatrix) (hsize : 0 < matrix.size)
    (hrowZero : matrix[0].size = matrix.size) :
    let mu := Generated.StrictRecombine.zeroQQMatrix matrix.size matrix.size
    let dot : QQ :=
      ((∑ k : Fin matrix.size,
        (matrix[0]!)[k.val]! * (matrix[0]!)[k.val]! : ZZ) : QQ)
    let norms := (Array.replicate matrix.size (0 : QQ)).set 0 dot (by simp [hsize])
    ConcreteGramSchmidtUpTo
      (Generated.StrictRecombine.LLLState.mk matrix #[] mu norms 1) 1 := by
  dsimp only
  intro rowCount hprocessed hmatrix
  by_cases hzero : rowCount = 0
  · subst rowCount
    ext row
    exact Fin.elim0 row
  have hone : rowCount = 1 := by omega
  subst rowCount
  funext row column
  have hrowFin : row = (0 : Fin 1) := Subsingleton.elim _ _
  have hcolumnFin : column = (0 : Fin 1) := Subsingleton.elim _ _
  subst row
  subst column
  simp [gramPrefixMatrixQQ, gramPrefixMatrix,
    Int.cast_sum, Matrix.mul_apply, Matrix.transpose_apply,
    Fin.sum_univ_one, gsLowerPrefix, gsNormDiagonal, Matrix.diagonal_apply,
    zeroQQMatrix_entry, getElem!_pos matrix 0 hsize, hrowZero]
  let actualDot : QQ := ∑ k : Fin matrix.size,
    (matrix[0][k.val]'(lt_of_lt_of_eq k.isLt hrowZero.symm) : QQ) *
      (matrix[0][k.val]'(lt_of_lt_of_eq k.isLt hrowZero.symm) : QQ)
  let actualNorms := (Array.replicate matrix.size (0 : QQ)).set 0 actualDot
    (by simp [hsize])
  have hactual : 0 < actualNorms.size := by simp [actualNorms, hsize]
  have hget : actualNorms[0]! = actualDot := by
    rw [getElem!_pos actualNorms 0 hactual]
    simp [actualNorms]
  simpa [actualDot, actualNorms] using hget.symm

noncomputable def prefixProductPotential (norms : Array QQ) (dimension : Nat) : QQ :=
  ∏ i : Fin dimension, prefixNormProduct norms (i.val + 1)

noncomputable def weightedNormPotential (norms : Array QQ) (dimension : Nat) : QQ :=
  ∏ i : Fin dimension, norms[i.val]! ^ (dimension - i.val)

/-- Concrete semantic invariant carried by the generated LLL execution.  It
states square integer storage, positive Gram–Schmidt norms, and the exact
Gram determinant identity for every prefix.  No result or termination trace
is stored here. -/
structure ConcreteLLLValid
    (state : Generated.StrictRecombine.LLLState) : Prop where
  norms_size : state.norms.size = state.matrix.size
  rows_square : ∀ row (hrow : row < state.matrix.size),
    state.matrix[row].size = state.matrix.size
  norms_positive : ∀ index (hindex : index < state.norms.size),
    0 < state.norms[index]
  gram_prefix : ∀ rowCount, rowCount ≤ state.matrix.size →
    ((Matrix.det (gramPrefixMatrix state.matrix rowCount) : Int) : QQ) =
      prefixNormProduct state.norms rowCount

theorem gsLowerPrefix_blockTriangular
    (state : Generated.StrictRecombine.LLLState) (rowCount : Nat) :
    (gsLowerPrefix state rowCount).BlockTriangular OrderDual.toDual := by
  intro i j hij
  unfold gsLowerPrefix
  have hij' : i.val < j.val := hij
  have hne : i ≠ j := by
    intro heq
    subst j
    omega
  rw [if_neg (by omega), if_neg hne]

theorem gsLowerPrefix_det
    (state : Generated.StrictRecombine.LLLState) (rowCount : Nat) :
    Matrix.det (gsLowerPrefix state rowCount) = 1 := by
  rw [Matrix.det_of_lowerTriangular _
    (gsLowerPrefix_blockTriangular state rowCount)]
  apply Finset.prod_eq_one
  intro i _
  simp [gsLowerPrefix]

theorem gsLowerNormMul_apply
    (state : Generated.StrictRecombine.LLLState) (rowCount : Nat)
    (row column : Fin rowCount) :
    (gsLowerPrefix state rowCount * gsNormDiagonal state rowCount *
        (gsLowerPrefix state rowCount).transpose) row column =
      ∑ index : Fin rowCount,
        gsLowerPrefix state rowCount row index * state.norms[index.val]! *
          gsLowerPrefix state rowCount column index := by
  rw [Matrix.mul_apply]
  apply Fintype.sum_congr
  intro index
  rw [Matrix.transpose_apply]
  congr 1
  simp [Matrix.mul_apply, gsNormDiagonal, Matrix.diagonal_apply]

theorem gsLowerPrefix_eq_of_mu_prefix
    (matrix : Generated.StrictRecombine.LLLMatrix)
    (transform : Generated.StrictRecombine.LLLMatrix)
    (mu beforeMu : Generated.StrictRecombine.QQMatrix)
    (norms : Array QQ) (k processed rowCount : Nat)
    (hrowCount : rowCount ≤ processed)
    (hmu : ∀ row, row < processed → mu[row]! = beforeMu[row]!) :
    gsLowerPrefix
        (Generated.StrictRecombine.LLLState.mk matrix transform mu norms k)
        rowCount =
      gsLowerPrefix
        (Generated.StrictRecombine.LLLState.mk matrix transform beforeMu norms k)
        rowCount := by
  funext row column
  unfold gsLowerPrefix
  split
  next hbelow =>
    rw [hmu row.val (lt_of_lt_of_le row.isLt hrowCount)]
  next hbelow => rfl

theorem gsNormDiagonal_eq_of_norm_prefix
    (matrix transform : Generated.StrictRecombine.LLLMatrix)
    (mu : Generated.StrictRecombine.QQMatrix)
    (norms beforeNorms : Array QQ) (k processed rowCount : Nat)
    (hrowCount : rowCount ≤ processed)
    (hnorms : ∀ index, index < processed →
      norms[index]! = beforeNorms[index]!) :
    gsNormDiagonal
        (Generated.StrictRecombine.LLLState.mk matrix transform mu norms k)
        rowCount =
      gsNormDiagonal
        (Generated.StrictRecombine.LLLState.mk matrix transform mu beforeNorms k)
        rowCount := by
  unfold gsNormDiagonal
  congr 1
  funext index
  exact hnorms index.val (lt_of_lt_of_le index.isLt hrowCount)

theorem concreteGramSchmidtUpTo_of_prefix_frames
    (matrix transform : Generated.StrictRecombine.LLLMatrix)
    (mu beforeMu : Generated.StrictRecombine.QQMatrix)
    (norms beforeNorms : Array QQ) (k processed : Nat)
    (hbefore : ConcreteGramSchmidtUpTo
      (Generated.StrictRecombine.LLLState.mk matrix transform beforeMu
        beforeNorms k) processed)
    (hmu : ∀ row, row < processed → mu[row]! = beforeMu[row]!)
    (hnorms : ∀ index, index < processed →
      norms[index]! = beforeNorms[index]!) :
    ConcreteGramSchmidtUpTo
      (Generated.StrictRecombine.LLLState.mk matrix transform mu norms k)
      processed := by
  intro rowCount hprocessed hmatrix
  have hfactor := hbefore rowCount hprocessed hmatrix
  have hlower := gsLowerPrefix_eq_of_mu_prefix matrix transform mu beforeMu
    norms k processed rowCount hprocessed hmu
  have hdiagonal := gsNormDiagonal_eq_of_norm_prefix matrix transform mu
    norms beforeNorms k processed rowCount hprocessed hnorms
  rw [hlower, hdiagonal]
  simpa [gsLowerPrefix, gsNormDiagonal] using hfactor

theorem finSum_eq_prefix_add_at_of_zero_after
    {size : Nat} (values : Fin size → QQ) (index : Nat)
    (hindex : index < size)
    (hzero : ∀ position : Fin size, index < position.val →
      values position = 0) :
    (∑ position : Fin size, values position) =
      (∑ position : Fin index,
        values ⟨position.val, lt_trans position.isLt hindex⟩) +
        values ⟨index, hindex⟩ := by
  let natValues : Nat → QQ := fun position =>
    if hposition : position < size then values ⟨position, hposition⟩ else 0
  have hfull : (∑ position : Fin size, values position) =
      ∑ position ∈ Finset.range size, natValues position := by
    rw [← Fin.sum_univ_eq_sum_range natValues size]
    apply Fintype.sum_congr
    intro position
    simp [natValues]
  have hprefix : (∑ position : Fin index,
      values ⟨position.val, lt_trans position.isLt hindex⟩) =
      ∑ position ∈ Finset.range index, natValues position := by
    rw [← Fin.sum_univ_eq_sum_range natValues index]
    apply Fintype.sum_congr
    intro position
    simp only [natValues, dif_pos (lt_trans position.isLt hindex)]
  have hat : natValues index = values ⟨index, hindex⟩ := by
    simp [natValues, hindex]
  rw [hfull, hprefix, ← hat, ← Finset.sum_range_succ]
  symm
  apply Finset.sum_subset (Finset.range_mono (by omega))
  intro position hposition hnotPrefix
  have hpositionSize : position < size := Finset.mem_range.mp hposition
  have hafter : index < position := by
    have hnot : ¬position < index + 1 := by
      simpa only [Finset.mem_range] using hnotPrefix
    omega
  simp only [natValues, dif_pos hpositionSize]
  exact hzero ⟨position, hpositionSize⟩ hafter

theorem gsLowerNormSum_new_row
    (state : Generated.StrictRecombine.LLLState) (i j : Nat) (hj : j < i) :
    (∑ index : Fin (i + 1),
      gsLowerPrefix state (i + 1) ⟨i, by omega⟩ index *
        state.norms[index.val]! *
        gsLowerPrefix state (i + 1) ⟨j, by omega⟩ index) =
      (∑ index : Fin j,
        (state.mu[i]!)[index.val]! * (state.mu[j]!)[index.val]! *
          state.norms[index.val]!) +
        (state.mu[i]!)[j]! * state.norms[j]! := by
  let values : Fin (i + 1) → QQ := fun index =>
    gsLowerPrefix state (i + 1) ⟨i, by omega⟩ index *
      state.norms[index.val]! *
      gsLowerPrefix state (i + 1) ⟨j, by omega⟩ index
  rw [show (∑ index : Fin (i + 1),
      gsLowerPrefix state (i + 1) ⟨i, by omega⟩ index *
        state.norms[index.val]! *
        gsLowerPrefix state (i + 1) ⟨j, by omega⟩ index) =
      ∑ index, values index by rfl]
  rw [finSum_eq_prefix_add_at_of_zero_after values j (by omega)]
  · congr 1
    · apply Fintype.sum_congr
      intro index
      simp [values, gsLowerPrefix, hj, index.isLt,
        show index.val < i by omega]
      ring
    · simp [values, gsLowerPrefix, hj]
  · intro position hafter
    have hne : (⟨j, by omega⟩ : Fin (i + 1)) ≠ position := by
      intro heq
      have hvalues : j = position.val :=
        congrArg (fun value : Fin (i + 1) => value.val) heq
      omega
    simp [values, gsLowerPrefix, hafter, show ¬position.val < j by omega,
      hne]

theorem gsLowerNormSum_new_diagonal
    (state : Generated.StrictRecombine.LLLState) (i : Nat) :
    (∑ index : Fin (i + 1),
      gsLowerPrefix state (i + 1) ⟨i, by omega⟩ index *
        state.norms[index.val]! *
        gsLowerPrefix state (i + 1) ⟨i, by omega⟩ index) =
      (∑ index : Fin i,
        (state.mu[i]!)[index.val]! * (state.mu[i]!)[index.val]! *
          state.norms[index.val]!) + state.norms[i]! := by
  rw [Fin.sum_univ_castSucc]
  congr 1
  · apply Fintype.sum_congr
    intro index
    simp [gsLowerPrefix, index.isLt]
    ring
  · have hlast : (⟨i, by omega⟩ : Fin (i + 1)) = Fin.last i := Fin.ext rfl
    simp [gsLowerPrefix, hlast]

theorem gsLowerNormSum_old_entry
    (state : Generated.StrictRecombine.LLLState) (i : Nat)
    (row column : Fin i) :
    (∑ index : Fin (i + 1),
      gsLowerPrefix state (i + 1) ⟨row.val, by omega⟩ index *
        state.norms[index.val]! *
        gsLowerPrefix state (i + 1) ⟨column.val, by omega⟩ index) =
      ∑ index : Fin i,
        gsLowerPrefix state i row index * state.norms[index.val]! *
          gsLowerPrefix state i column index := by
  rw [Fin.sum_univ_castSucc]
  have hrowNe : (⟨row.val, by omega⟩ : Fin (i + 1)) ≠ Fin.last i := by
    intro heq
    have := congrArg Fin.val heq
    simp only [Fin.val_last] at this
    omega
  have hlastZero :
      gsLowerPrefix state (i + 1) ⟨row.val, by omega⟩ (Fin.last i) *
          state.norms[(Fin.last i).val]! *
          gsLowerPrefix state (i + 1) ⟨column.val, by omega⟩ (Fin.last i) = 0 := by
    have hleft :
        gsLowerPrefix state (i + 1) ⟨row.val, by omega⟩ (Fin.last i) = 0 := by
      unfold gsLowerPrefix
      have hnotBelow : ¬(Fin.last i).val < row.val := by
        simp only [Fin.val_last]
        omega
      rw [if_neg hnotBelow, if_neg hrowNe]
    simp [hleft]
  rw [hlastZero, add_zero]
  apply Fintype.sum_congr
  intro index
  have hrowEq :
      gsLowerPrefix state (i + 1) ⟨row.val, by omega⟩ index.castSucc =
        gsLowerPrefix state i row index := by
    unfold gsLowerPrefix
    simp only [Fin.val_castSucc]
    by_cases hbelow : index.val < row.val
    · simp [hbelow]
    · by_cases heq : row = index
      · subst row
        have hcastEq : (⟨index.val, by omega⟩ : Fin (i + 1)) =
            index.castSucc := Fin.ext rfl
        simp [hbelow, hcastEq]
      · have hbigNe : (⟨row.val, by omega⟩ : Fin (i + 1)) ≠
            index.castSucc := by
          intro hbig
          apply heq
          apply Fin.ext
          exact congrArg (fun value : Fin (i + 1) => value.val) hbig
        simp [hbelow, heq, hbigNe]
  have hcolumnEq :
      gsLowerPrefix state (i + 1) ⟨column.val, by omega⟩ index.castSucc =
        gsLowerPrefix state i column index := by
    unfold gsLowerPrefix
    simp only [Fin.val_castSucc]
    by_cases hbelow : index.val < column.val
    · simp [hbelow]
    · by_cases heq : column = index
      · subst column
        have hcastEq : (⟨index.val, by omega⟩ : Fin (i + 1)) =
            index.castSucc := Fin.ext rfl
        simp [hbelow, hcastEq]
      · have hbigNe : (⟨column.val, by omega⟩ : Fin (i + 1)) ≠
            index.castSucc := by
          intro hbig
          apply heq
          apply Fin.ext
          exact congrArg (fun value : Fin (i + 1) => value.val) hbig
        simp [hbelow, heq, hbigNe]
  rw [hrowEq, hcolumnEq]
  simp only [Fin.val_castSucc]

theorem gramPrefixMatrixQQ_apply_eq_sourceRowDot
    (matrix : Generated.StrictRecombine.LLLMatrix) (rowCount : Nat)
    (row column : Fin rowCount)
    (hrow : row.val < matrix.size)
    (hcolumn : column.val < matrix.size)
    (hrowSize : matrix[row.val].size = matrix.size)
    (hcolumnSize : matrix[column.val].size = matrix.size) :
    gramPrefixMatrixQQ matrix rowCount row column =
      sourceRowDot matrix row.val column.val matrix[row.val]!.size := by
  unfold gramPrefixMatrixQQ gramPrefixMatrix sourceRowDot
  simp only [Matrix.map_apply]
  rw [getElem!_pos matrix row.val hrow, getElem!_pos matrix column.val hcolumn,
    hrowSize]

theorem concreteGramSchmidtUpTo_extend_one
    (state : Generated.StrictRecombine.LLLState) (i : Nat)
    (hi : i < state.matrix.size)
    (hrowsSquare : ∀ row (hrow : row < state.matrix.size),
      state.matrix[row].size = state.matrix.size)
    (hprefix : ConcreteGramSchmidtUpTo state i)
    (hrow : GramMuPrefixCorrect state.matrix state.mu state.norms i i)
    (hdiagonal : sourceRowDot state.matrix i i state.matrix[i]!.size =
      (∑ index : Fin i,
        (state.mu[i]!)[index.val]! * (state.mu[i]!)[index.val]! *
          state.norms[index.val]!) + state.norms[i]!) :
    ConcreteGramSchmidtUpTo state (i + 1) := by
  intro rowCount hprocessed hmatrix
  by_cases hold : rowCount ≤ i
  · exact hprefix rowCount hold hmatrix
  have hrowCount : rowCount = i + 1 := by omega
  subst rowCount
  funext row column
  rw [gsLowerNormMul_apply]
  let last : Fin (i + 1) := ⟨i, by omega⟩
  have hlastVal : last.val = i := rfl
  have hnewRow : ∀ position : Fin (i + 1), position.val < i →
      gramPrefixMatrixQQ state.matrix (i + 1) last position =
        ∑ index : Fin (i + 1),
          gsLowerPrefix state (i + 1) last index *
            state.norms[index.val]! *
            gsLowerPrefix state (i + 1) position index := by
    intro position hposition
    calc
      gramPrefixMatrixQQ state.matrix (i + 1) last position =
          sourceRowDot state.matrix i position.val state.matrix[i]!.size := by
        exact gramPrefixMatrixQQ_apply_eq_sourceRowDot state.matrix (i + 1)
          last position hi (lt_trans hposition hi)
          (hrowsSquare i hi) (hrowsSquare position.val (lt_trans hposition hi))
      _ = (∑ index : Fin position.val,
            (state.mu[i]!)[index.val]! *
              (state.mu[position.val]!)[index.val]! *
              state.norms[index.val]!) +
            (state.mu[i]!)[position.val]! * state.norms[position.val]! :=
        hrow position.val hposition
      _ = ∑ index : Fin (i + 1),
          gsLowerPrefix state (i + 1) last index *
            state.norms[index.val]! *
            gsLowerPrefix state (i + 1) position index := by
        simpa [last] using (gsLowerNormSum_new_row state i position.val hposition).symm
  by_cases hrowLast : row.val = i
  · have hrowEq : row = last := by apply Fin.ext; exact hrowLast
    subst row
    by_cases hcolumnLast : column.val = i
    · have hcolumnEq : column = last := by apply Fin.ext; exact hcolumnLast
      subst column
      calc
        gramPrefixMatrixQQ state.matrix (i + 1) last last =
            sourceRowDot state.matrix i i state.matrix[i]!.size := by
          exact gramPrefixMatrixQQ_apply_eq_sourceRowDot state.matrix (i + 1)
            last last hi hi (hrowsSquare i hi) (hrowsSquare i hi)
        _ = (∑ index : Fin i,
              (state.mu[i]!)[index.val]! * (state.mu[i]!)[index.val]! *
                state.norms[index.val]!) + state.norms[i]! := hdiagonal
        _ = ∑ index : Fin (i + 1),
            gsLowerPrefix state (i + 1) last index *
              state.norms[index.val]! *
              gsLowerPrefix state (i + 1) last index := by
          simpa [last] using (gsLowerNormSum_new_diagonal state i).symm
    · exact hnewRow column (by omega)
  · have hrowBefore : row.val < i := by omega
    by_cases hcolumnLast : column.val = i
    · have hcolumnEq : column = last := by apply Fin.ext; exact hcolumnLast
      subst column
      calc
        gramPrefixMatrixQQ state.matrix (i + 1) row last =
            gramPrefixMatrixQQ state.matrix (i + 1) last row := by
          unfold gramPrefixMatrixQQ gramPrefixMatrix
          simp only [Matrix.map_apply]
          congr 1
          apply Fintype.sum_congr
          intro index
          ring
        _ = ∑ index : Fin (i + 1),
            gsLowerPrefix state (i + 1) last index *
              state.norms[index.val]! *
              gsLowerPrefix state (i + 1) row index :=
          hnewRow row hrowBefore
        _ = ∑ index : Fin (i + 1),
            gsLowerPrefix state (i + 1) row index *
              state.norms[index.val]! *
              gsLowerPrefix state (i + 1) last index := by
          apply Fintype.sum_congr
          intro index
          ring
    · have hcolumnBefore : column.val < i := by omega
      let oldRow : Fin i := ⟨row.val, hrowBefore⟩
      let oldColumn : Fin i := ⟨column.val, hcolumnBefore⟩
      have holdFactor := hprefix i le_rfl (by omega)
      have holdEntry := congrArg
        (fun matrix => matrix oldRow oldColumn) holdFactor
      dsimp only at holdEntry
      rw [gsLowerNormMul_apply] at holdEntry
      calc
        gramPrefixMatrixQQ state.matrix (i + 1) row column =
            gramPrefixMatrixQQ state.matrix i oldRow oldColumn := by
          rfl
        _ = ∑ index : Fin i,
            gsLowerPrefix state i oldRow index * state.norms[index.val]! *
              gsLowerPrefix state i oldColumn index := holdEntry
        _ = ∑ index : Fin (i + 1),
            gsLowerPrefix state (i + 1) row index *
              state.norms[index.val]! *
              gsLowerPrefix state (i + 1) column index := by
          simpa [oldRow, oldColumn] using
            (gsLowerNormSum_old_entry state i oldRow oldColumn).symm

theorem ConcreteGramSchmidtUpTo.gram_prefix
    {state : Generated.StrictRecombine.LLLState} {processed : Nat}
    (hgs : ConcreteGramSchmidtUpTo state processed)
    (rowCount : Nat) (hprocessed : rowCount ≤ processed)
    (hrowCount : rowCount ≤ state.matrix.size) :
    ((Matrix.det (gramPrefixMatrix state.matrix rowCount) : Int) : QQ) =
      prefixNormProduct state.norms rowCount := by
  have hfactor := hgs rowCount hprocessed hrowCount
  calc
    ((Matrix.det (gramPrefixMatrix state.matrix rowCount) : Int) : QQ) =
        Matrix.det (gramPrefixMatrixQQ state.matrix rowCount) := by
      exact Int.cast_det _
    _ = Matrix.det (gsLowerPrefix state rowCount *
          gsNormDiagonal state rowCount *
            (gsLowerPrefix state rowCount).transpose) := by rw [hfactor]
    _ = Matrix.det (gsLowerPrefix state rowCount) *
          Matrix.det (gsNormDiagonal state rowCount) *
            Matrix.det (gsLowerPrefix state rowCount) := by
      rw [Matrix.det_mul, Matrix.det_mul, Matrix.det_transpose]
    _ = prefixNormProduct state.norms rowCount := by
      rw [gsLowerPrefix_det]
      simp only [one_mul, mul_one]
      unfold gsNormDiagonal
      rw [Matrix.det_diagonal]
      unfold prefixNormProduct
      apply Finset.prod_congr rfl
      intro index _
      rfl

theorem ConcreteGramSchmidt.gram_prefix
    {state : Generated.StrictRecombine.LLLState}
    (hgs : ConcreteGramSchmidt state)
    (rowCount : Nat) (hrowCount : rowCount ≤ state.matrix.size) :
    ((Matrix.det (gramPrefixMatrix state.matrix rowCount) : Int) : QQ) =
      prefixNormProduct state.norms rowCount := by
  have hfactor := hgs rowCount hrowCount
  calc
    ((Matrix.det (gramPrefixMatrix state.matrix rowCount) : Int) : QQ) =
        Matrix.det (gramPrefixMatrixQQ state.matrix rowCount) := by
      exact Int.cast_det _
    _ = Matrix.det (gsLowerPrefix state rowCount *
          gsNormDiagonal state rowCount *
            (gsLowerPrefix state rowCount).transpose) := by rw [hfactor]
    _ = Matrix.det (gsLowerPrefix state rowCount) *
          Matrix.det (gsNormDiagonal state rowCount) *
            Matrix.det (gsLowerPrefix state rowCount) := by
      rw [Matrix.det_mul, Matrix.det_mul, Matrix.det_transpose]
    _ = prefixNormProduct state.norms rowCount := by
      rw [gsLowerPrefix_det]
      simp only [one_mul, mul_one]
      unfold gsNormDiagonal
      rw [Matrix.det_diagonal]
      unfold prefixNormProduct
      apply Finset.prod_congr rfl
      intro i _
      rfl

/-- Full invariant intended for the generated well-founded loop.  Unlike
`ConcreteLLLValid`, this stores the executable Gram–Schmidt relationship from
which the determinant field is derived. -/
structure ConcreteLLLExecutionValid
    (state : Generated.StrictRecombine.LLLState) : Prop where
  norms_size : state.norms.size = state.matrix.size
  mu_size : state.mu.size = state.matrix.size
  rows_square : ∀ row (hrow : row < state.matrix.size),
    state.matrix[row].size = state.matrix.size
  mu_rows_square : ∀ row (hrow : row < state.mu.size),
    state.mu[row].size = state.matrix.size
  norms_positive : ∀ index (hindex : index < state.norms.size),
    0 < state.norms[index]
  gram_schmidt : ConcreteGramSchmidt state

/-- The source precondition of `__lll_reduce`: a square integer basis of full
row rank.  Unlike the execution invariant, this mentions no preselected
Gram--Schmidt output and therefore can be reused when C++ reinitializes LLL on
the matrix retained from a previous van-Hoeij round. -/
structure ConcreteLLLInputValid
    (matrix : Generated.StrictRecombine.LLLMatrix) : Prop where
  rows_square : ∀ row (hrow : row < matrix.size),
    matrix[row].size = matrix.size
  determinant_ne : Matrix.det
    (basisPrefixMatrix matrix matrix.size matrix.size) ≠ 0

private structure InitialMatrixPrefix (scale : ZZ) (size index : Nat)
    (matrix : Generated.StrictRecombine.LLLMatrix) : Prop where
  matrix_size : matrix.size = size
  rows_square : ∀ row (hrow : row < size), matrix[row]!.size = size
  entry : ∀ row column (hrow : row < size) (hcolumn : column < size),
    (matrix[row]!)[column]! = if row = column ∧ row < index then scale else 0

private theorem zeroMatrix_initial_prefix (scale : ZZ) (size : Nat) :
    InitialMatrixPrefix scale size 0
      (Generated.StrictRecombine.zeroMatrix size size) := by
  refine ⟨zeroMatrix_size size size, ?_, ?_⟩
  · intro row hrow
    rw [zeroMatrix_get size size row hrow, zeroMatrixRow_size]
  · intro row column hrow hcolumn
    rw [zeroMatrix_entry size size row column hrow hcolumn]
    simp

private theorem InitialMatrixPrefix.set_next
    {scale : ZZ} {size index : Nat}
    {matrix : Generated.StrictRecombine.LLLMatrix}
    (hprefix : InitialMatrixPrefix scale size index matrix)
    (hindex : index < size)
    (hmatrixIndex : index < matrix.size)
    (hrowIndex : index < matrix[index].size) :
    InitialMatrixPrefix scale size (index + 1)
      (matrix.set index (matrix[index].set index scale hrowIndex)
        hmatrixIndex) := by
  refine ⟨by simp [hprefix.matrix_size], ?_, ?_⟩
  · intro row hrow
    have hrowMatrix : row < matrix.size := by
      rw [hprefix.matrix_size]
      exact hrow
    have hnewRow : row <
        (matrix.set index (matrix[index].set index scale hrowIndex)
          hmatrixIndex).size := by simp [hprefix.matrix_size, hrow]
    by_cases heq : index = row
    · subst row
      rw [getElem!_pos _ index hnewRow, Array.getElem_set_self,
        Array.size_set]
      simpa [getElem!_pos matrix index hmatrixIndex] using
        hprefix.rows_square index hindex
    · rw [getElem!_pos _ row hnewRow,
        Array.getElem_set_ne hmatrixIndex hrowMatrix heq]
      have hold := hprefix.rows_square row hrow
      rw [getElem!_pos matrix row hrowMatrix] at hold
      exact hold
  · intro row column hrow hcolumn
    have hrowMatrix : row < matrix.size := by
      rw [hprefix.matrix_size]
      exact hrow
    have hnewRow : row <
        (matrix.set index (matrix[index].set index scale hrowIndex)
          hmatrixIndex).size := by simp [hprefix.matrix_size, hrow]
    by_cases hrowEq : index = row
    · subst row
      have hcolumnOld : column < matrix[index].size := by
        have hsize := hprefix.rows_square index hindex
        rw [getElem!_pos matrix index hmatrixIndex] at hsize
        omega
      by_cases hcolumnEq : index = column
      · subst column
        rw [getElem!_pos
          (matrix.set index (matrix[index].set index scale hrowIndex)
            hmatrixIndex) index hnewRow]
        rw [Array.getElem_set_self]
        rw [getElem!_pos _ index (by simp [hrowIndex]),
          Array.getElem_set_self]
        simp [hindex]
      · rw [getElem!_pos _ index hnewRow, Array.getElem_set_self,
          getElem!_pos _ column (by simp [Array.size_set, hcolumnOld])]
        rw [Array.getElem_set_ne hrowIndex hcolumnOld hcolumnEq]
        have hold := hprefix.entry index column hindex hcolumn
        rw [getElem!_pos matrix index hmatrixIndex,
          getElem!_pos matrix[index] column hcolumnOld] at hold
        rw [hold]
        simp [hcolumnEq]
    · rw [getElem!_pos _ row hnewRow,
        Array.getElem_set_ne hmatrixIndex hrowMatrix hrowEq]
      have hcolumnOld : column < matrix[row].size := by
        have hsize := hprefix.rows_square row hrow
        rw [getElem!_pos matrix row hrowMatrix] at hsize
        omega
      rw [getElem!_pos matrix[row] column hcolumnOld]
      have hold := hprefix.entry row column hrow hcolumn
      rw [getElem!_pos matrix row hrowMatrix,
        getElem!_pos matrix[row] column hcolumnOld] at hold
      rw [hold]
      by_cases hdiag : row = column
      · subst column
        simp only [true_and]
        have : row ≠ index := by exact fun h => hrowEq h.symm
        have hLt : row < index + 1 ↔ row < index := by omega
        by_cases hri : row < index
        · have hrnext : row < index + 1 := hLt.mpr hri
          simp [hri, hrnext]
        · have hrnext : ¬ row < index + 1 := by rwa [hLt]
          simp [hri, hrnext]
      · simp [hdiag]

private theorem setInitialDiagonalLoop_prefix
    (scale : ZZ) (size index : Nat)
    (matrix output : Generated.StrictRecombine.LLLMatrix)
    (hprefix : InitialMatrixPrefix scale size index matrix)
    (hindexLe : index ≤ size)
    (hrun : Generated.StrictRecombine.setInitialDiagonalLoop scale size index
      matrix = .ok output) :
    InitialMatrixPrefix scale size size output := by
  induction hmeasure : size - index using Nat.strong_induction_on
      generalizing index matrix output with
  | h measure ih =>
      rw [Generated.StrictRecombine.setInitialDiagonalLoop] at hrun
      by_cases hindex : index < size
      · rw [dif_pos hindex] at hrun
        have hmatrixIndex : index < matrix.size := by
          rw [hprefix.matrix_size]
          exact hindex
        rw [dif_pos hmatrixIndex] at hrun
        have hrowIndex : index < matrix[index].size := by
          simpa [getElem!_pos matrix index hmatrixIndex] using
            (show index < matrix[index]!.size by
              rw [hprefix.rows_square index hindex]
              exact hindex)
        rw [dif_pos hrowIndex] at hrun
        exact ih (size - (index + 1)) (by omega) (index + 1)
          (matrix.set index (matrix[index].set index scale)) output
          (hprefix.set_next hindex hmatrixIndex hrowIndex) (by omega) hrun rfl
      · rw [dif_neg hindex] at hrun
        have hout := Except.ok.inj hrun
        subst output
        have hsize : index = size := by omega
        simpa [hsize] using hprefix

theorem makeInitialMatrix_prefix (size : Nat) (scale : ZZ)
    (matrix : Generated.StrictRecombine.LLLMatrix)
    (hrun : Generated.StrictRecombine.makeInitialMatrix size scale = .ok matrix) :
    InitialMatrixPrefix scale size size matrix := by
  unfold Generated.StrictRecombine.makeInitialMatrix at hrun
  exact setInitialDiagonalLoop_prefix scale size 0
    (Generated.StrictRecombine.zeroMatrix size size) matrix
    (zeroMatrix_initial_prefix scale size) (by omega) hrun

theorem makeInitialMatrix_input_valid (size : Nat) (scale : ZZ)
    (matrix : Generated.StrictRecombine.LLLMatrix) (hscale : scale ≠ 0)
    (hrun : Generated.StrictRecombine.makeInitialMatrix size scale = .ok matrix) :
    ConcreteLLLInputValid matrix ∧ matrix.size = size := by
  have hprefix := makeInitialMatrix_prefix size scale matrix hrun
  have hmatrixSize := hprefix.matrix_size
  refine ⟨⟨?_, ?_⟩, hmatrixSize⟩
  · intro row hrow
    have hrowSize := hprefix.rows_square row (by omega)
    rw [getElem!_pos matrix row hrow] at hrowSize
    exact hrowSize.trans hmatrixSize.symm
  · have hbasis : basisPrefixMatrix matrix matrix.size matrix.size =
        Matrix.diagonal (fun _ : Fin matrix.size => scale) := by
      apply Matrix.ext
      intro row column
      simp only [basisPrefixMatrix, Matrix.diagonal_apply]
      have hentry := hprefix.entry row.val column.val (by omega) (by omega)
      by_cases heq : row = column
      · subst column
        have hrowSize : row.val < size := by omega
        simpa [hrowSize] using hentry
      · have hval : row.val ≠ column.val := by
          exact fun h => heq (Fin.ext h)
        simpa [heq, hval] using hentry
    rw [hbasis, Matrix.det_diagonal]
    exact Finset.prod_ne_zero_iff.mpr (by
      intro index _
      exact hscale)

theorem resetVanHoeijLattice_input_valid (factorCount : Nat)
    (matrix : Generated.StrictRecombine.LLLMatrix) (bound : ZZ)
    (hrun : Generated.StrictRecombine.resetVanHoeijLattice factorCount =
      .ok (matrix, bound)) :
    ConcreteLLLInputValid matrix ∧ matrix.size = factorCount := by
  unfold Generated.StrictRecombine.resetVanHoeijLattice at hrun
  let exponent := Generated.StrictRecombine.vanHoeijExponent factorCount
  change (match Generated.StrictRecombine.makeInitialMatrix factorCount
      ((2 : ZZ) ^ exponent) with
    | Except.error fault => (Except.error fault :
        RawExec (Generated.StrictRecombine.LLLMatrix × ZZ))
    | Except.ok initial => Except.ok
        (initial, Generated.StrictRecombine.vanHoeijBound factorCount)) =
      Except.ok (matrix, bound) at hrun
  cases hmake : Generated.StrictRecombine.makeInitialMatrix factorCount
      ((2 : ZZ) ^ exponent) with
  | error fault =>
      rw [hmake] at hrun
      contradiction
  | ok initial =>
      rw [hmake] at hrun
      have hout := Except.ok.inj hrun
      injection hout with hmatrix hbound
      subst matrix
      exact makeInitialMatrix_input_valid factorCount
        ((2 : ZZ) ^ exponent)
        initial (Int.pow_ne_zero (by norm_num)) hmake

theorem appendCldColumn_input_valid
    (matrix output : Generated.StrictRecombine.LLLMatrix)
    (cld : Array SparsePolyZZ) (existingColumns spiralDegree : Nat)
    (hinput : ConcreteLLLInputValid matrix)
    (hdimension : matrix.size = cld.size + existingColumns)
    (hrun : Generated.StrictRecombine.appendCldColumn matrix cld
      existingColumns spiralDegree = .ok output) :
    ConcreteLLLInputValid output ∧ output.size = matrix.size + 1 := by
  rcases appendCldColumn_shape matrix output cld existingColumns spiralDegree
      hdimension hrun with ⟨finalRow, hfinalSize, hfinalLast, rfl⟩
  have houtputSize :
      ((appendZeroSuffix matrix 0).push finalRow).size = matrix.size + 1 := by
    rw [Array.size_push, appendZeroSuffix_size]
  refine ⟨?_, houtputSize⟩
  constructor
  · intro row hrow
    by_cases hold : row < matrix.size
    · have hrowSize : matrix[row]!.size = matrix.size := by
        rw [getElem!_pos _ row hold]
        exact hinput.rows_square row hold
      rw [← getElem!_pos _ row hrow,
        arrayPush_getElem!_lt _ _ row (by
          simpa [appendZeroSuffix_size] using hold),
        appendZeroSuffix_row matrix row hold, Array.size_push,
        hrowSize, houtputSize]
    · have hrowLast : row = matrix.size := by
        rw [houtputSize] at hrow
        omega
      subst row
      rw [← getElem!_pos _ matrix.size hrow,
        show matrix.size = (appendZeroSuffix matrix 0).size by
          rw [appendZeroSuffix_size], arrayPush_getElem!_last]
      exact hfinalSize.trans houtputSize.symm
  · rw [houtputSize]
    have hdet : Matrix.det
        (basisPrefixMatrix ((appendZeroSuffix matrix 0).push finalRow)
          (matrix.size + 1) (matrix.size + 1)) =
        Matrix.det (basisPrefixMatrix matrix matrix.size matrix.size) := by
      apply det_append_unit_lower matrix.size
        (basisPrefixMatrix matrix matrix.size matrix.size)
        (fun column => finalRow[column.val]!)
      · intro i j
        simp only [basisPrefixMatrix]
        have houter :
            ((appendZeroSuffix matrix 0).push finalRow)[i.val]! =
              matrix[i.val]!.push 0 := by
          rw [arrayPush_getElem!_lt _ _ i.val (by
            simpa [appendZeroSuffix_size] using i.isLt),
            appendZeroSuffix_row matrix i.val i.isLt]
        simp only [Fin.val_castSucc]
        rw [houter]
        have hrowSize : matrix[i.val]!.size = matrix.size := by
          rw [getElem!_pos _ i.val i.isLt]
          exact hinput.rows_square i.val i.isLt
        have hcolumn : j.val < matrix[i.val]!.size := by
          rw [hrowSize]
          exact j.isLt
        exact arrayPush_getElem!_lt _ _ j.val hcolumn
      · intro i
        simp only [basisPrefixMatrix]
        have houter :
            ((appendZeroSuffix matrix 0).push finalRow)[i.val]! =
              matrix[i.val]!.push 0 := by
          rw [arrayPush_getElem!_lt _ _ i.val (by
            simpa [appendZeroSuffix_size] using i.isLt),
            appendZeroSuffix_row matrix i.val i.isLt]
        simp only [Fin.val_castSucc, Fin.val_last]
        rw [houter]
        have hrowSize : matrix[i.val]!.size = matrix.size := by
          rw [getElem!_pos _ i.val i.isLt]
          exact hinput.rows_square i.val i.isLt
        have hlast := arrayPush_getElem!_last matrix[i.val]! (0 : ZZ)
        rw [hrowSize] at hlast
        exact hlast
      · intro j
        simp only [basisPrefixMatrix]
        have houter :
            ((appendZeroSuffix matrix 0).push finalRow)[matrix.size]! =
              finalRow := by
          rw [show matrix.size = (appendZeroSuffix matrix 0).size by
            rw [appendZeroSuffix_size], arrayPush_getElem!_last]
        simp only [Fin.val_last, Fin.val_castSucc]
        rw [houter]
      · simp only [basisPrefixMatrix]
        have houter :
            ((appendZeroSuffix matrix 0).push finalRow)[matrix.size]! =
              finalRow := by
          rw [show matrix.size = (appendZeroSuffix matrix 0).size by
            rw [appendZeroSuffix_size], arrayPush_getElem!_last]
        simp only [Fin.val_last]
        rw [houter]
        exact hfinalLast
    rw [hdet]
    exact hinput.determinant_ne

theorem buildCldMatrixLoop_input_valid
    (cld : Array SparsePolyZZ) (current target width added : Nat)
    (matrix output : Generated.StrictRecombine.LLLMatrix)
    (finalAdded : Nat) (hinput : ConcreteLLLInputValid matrix)
    (hdimension : matrix.size = cld.size + current + added)
    (hrun : Generated.StrictRecombine.buildCldMatrixLoop cld current target
      width added matrix = .ok (output, finalAdded)) :
    ConcreteLLLInputValid output ∧
      output.size = cld.size + current + finalAdded := by
  induction hmeasure : target - added using Nat.strong_induction_on
      generalizing added matrix output finalAdded with
  | h measure ih =>
      rw [Generated.StrictRecombine.buildCldMatrixLoop] at hrun
      by_cases htarget : added < target
      · rw [dif_pos htarget] at hrun
        dsimp only at hrun
        by_cases hwidth : current + added < width
        · rw [dif_pos hwidth] at hrun
          let degree := if (current + added) % 2 = 0 then
              (current + added) / 2
            else width - 1 - (current + added - 1) / 2
          cases happend : Generated.StrictRecombine.appendCldColumn matrix cld
              (current + added) degree with
          | error fault =>
              rw [happend] at hrun
              change (Except.error fault : RawExec
                (Generated.StrictRecombine.LLLMatrix × Nat)) =
                  .ok (output, finalAdded) at hrun
              contradiction
          | ok next =>
              rw [happend] at hrun
              change Generated.StrictRecombine.buildCldMatrixLoop cld current
                target width (added + 1) next = .ok (output, finalAdded) at hrun
              have hstep := appendCldColumn_input_valid matrix next cld
                (current + added) degree hinput (by omega) happend
              exact ih (target - (added + 1)) (by omega) (added + 1) next
                output finalAdded hstep.1 (by omega) hrun rfl
        · rw [dif_neg hwidth] at hrun
          have hout := Except.ok.inj hrun
          injection hout with hmatrix hadd
          subst output
          subst finalAdded
          exact ⟨hinput, hdimension⟩
      · rw [dif_neg htarget] at hrun
        have hout := Except.ok.inj hrun
        injection hout with hmatrix hadd
        subst output
        subst finalAdded
        exact ⟨hinput, hdimension⟩

theorem buildCldMatrix_input_valid
    (matrix output : Generated.StrictRecombine.LLLMatrix)
    (cld : Array SparsePolyZZ) (current target added : Nat)
    (hinput : ConcreteLLLInputValid matrix)
    (hdimension : matrix.size = cld.size + current)
    (hrun : Generated.StrictRecombine.buildCldMatrix matrix cld current target =
      .ok (output, added)) :
    ConcreteLLLInputValid output ∧
      output.size = cld.size + current + added := by
  unfold Generated.StrictRecombine.buildCldMatrix at hrun
  exact buildCldMatrixLoop_input_valid cld current target
    (Generated.StrictRecombine.cldSpiralWidth cld) 0 matrix output added
    hinput (by simpa using hdimension) hrun

set_option maxHeartbeats 800000 in
theorem ConcreteLLLInputValid.rational_prefix_rows_linearIndependent
    {matrix : Generated.StrictRecombine.LLLMatrix}
    (hinput : ConcreteLLLInputValid matrix) (rowCount : Nat)
    (hrowCount : rowCount ≤ matrix.size) :
    LinearIndependent QQ (basisPrefixMatrixQQ matrix rowCount).row := by
  have hrowsInt : LinearIndependent ZZ
      (basisPrefixMatrix matrix matrix.size matrix.size).row :=
    Matrix.linearIndependent_rows_of_det_ne_zero hinput.determinant_ne
  have hrowsQQ : LinearIndependent QQ
      (basisPrefixMatrixQQ matrix matrix.size).row := by
    change LinearIndependent QQ
      (fun row : Fin matrix.size => fun column : Fin matrix.size =>
        ((basisPrefixMatrix matrix matrix.size matrix.size row column : ZZ) : QQ))
    have hcast :=
      (linearIndependent_algebraMap_comp_iff (R := ZZ) (S := QQ)).2 hrowsInt
    simpa [Function.comp_def] using hcast
  let embed : Fin rowCount → Fin matrix.size :=
    fun row => ⟨row.val, lt_of_lt_of_le row.isLt hrowCount⟩
  have hembed : Function.Injective embed := by
    intro left right hequal
    have hvalues : left.val = right.val :=
      congrArg (fun value : Fin matrix.size => value.val) hequal
    exact Fin.ext hvalues
  have hprefix := hrowsQQ.comp embed hembed
  have hrowEq : (basisPrefixMatrixQQ matrix rowCount).row =
      fun row => (basisPrefixMatrixQQ matrix matrix.size).row (embed row) := by
    funext row column
    simp [basisPrefixMatrixQQ, basisPrefixMatrix, embed,
      getElem!_pos matrix row.val (lt_of_lt_of_le row.isLt hrowCount)]
  rw [hrowEq]
  exact hprefix

theorem ConcreteLLLInputValid.gramPrefixMatrixQQ_posDef
    {matrix : Generated.StrictRecombine.LLLMatrix}
    (hinput : ConcreteLLLInputValid matrix) (rowCount : Nat)
    (hrowCount : rowCount ≤ matrix.size) :
    Matrix.PosDef (gramPrefixMatrixQQ matrix rowCount) := by
  rw [gramPrefixMatrixQQ_eq_mul_transpose]
  apply Matrix.PosDef.mul_conjTranspose_self
  rw [Matrix.vecMul_injective_iff]
  exact hinput.rational_prefix_rows_linearIndependent rowCount hrowCount

set_option maxHeartbeats 800000 in
theorem ConcreteLLLInputValid.gramPrefixMatrixQQ_det_pos
    {matrix : Generated.StrictRecombine.LLLMatrix}
    (hinput : ConcreteLLLInputValid matrix) (rowCount : Nat)
    (hrowCount : rowCount ≤ matrix.size) :
    0 < Matrix.det (gramPrefixMatrixQQ matrix rowCount) := by
  let rationalBasis := basisPrefixMatrixQQ matrix rowCount
  let realBasis : Matrix (Fin rowCount) (Fin matrix.size) ℝ :=
    rationalBasis.map fun value : QQ => (value : ℝ)
  have hrowsQQ := hinput.rational_prefix_rows_linearIndependent
    rowCount hrowCount
  have hrowsReal : LinearIndependent ℝ realBasis.row := by
    change LinearIndependent ℝ
      (fun row column => ((rationalBasis row column : QQ) : ℝ))
    have hcast :=
      (linearIndependent_algebraMap_comp_iff (R := QQ) (S := ℝ)).2 hrowsQQ
    simpa [Function.comp_def, rationalBasis] using hcast
  have hrealPos : Matrix.PosDef (realBasis * realBasis.transpose) := by
    apply Matrix.PosDef.mul_conjTranspose_self
    rw [Matrix.vecMul_injective_iff]
    exact hrowsReal
  have hcastGram :
      (gramPrefixMatrixQQ matrix rowCount).map (fun value : QQ => (value : ℝ)) =
        realBasis * realBasis.transpose := by
    rw [gramPrefixMatrixQQ_eq_mul_transpose]
    funext row column
    simp [realBasis, rationalBasis, Matrix.mul_apply, Matrix.transpose_apply,
      Rat.cast_sum]
  apply (Rat.cast_pos (K := ℝ)).mp
  rw [Rat.cast_det]
  rw [hcastGram]
  exact hrealPos.det_pos

theorem ConcreteLLLInputValid.first_norm_positive
    {matrix : Generated.StrictRecombine.LLLMatrix}
    (hinput : ConcreteLLLInputValid matrix) (hsize : 0 < matrix.size) :
    0 < (((∑ k : Fin matrix.size,
      (matrix[0]!)[k.val]! * (matrix[0]!)[k.val]!) : ZZ) : QQ) := by
  have hrowsLI : LinearIndependent ZZ
      (fun row : Fin matrix.size =>
        basisPrefixMatrix matrix matrix.size matrix.size row) :=
    Matrix.linearIndependent_rows_of_det_ne_zero hinput.determinant_ne
  have hrowNe := LinearIndependent.ne_zero (⟨0, hsize⟩ : Fin matrix.size) hrowsLI
  have hentry : ∃ column : Fin matrix.size,
      (matrix[0]!)[column.val]! ≠ 0 := by
    by_contra hnone
    push_neg at hnone
    apply hrowNe
    funext column
    unfold basisPrefixMatrix
    simpa using hnone column
  have hsum : (0 : ZZ) < ∑ k : Fin matrix.size,
      (matrix[0]!)[k.val]! * (matrix[0]!)[k.val]! := by
    have hnonnegative : (0 : ZZ) ≤ ∑ k : Fin matrix.size,
        (matrix[0]!)[k.val]! * (matrix[0]!)[k.val]! := by
      apply Finset.sum_nonneg
      intro column hcolumn
      show (0 : ZZ) ≤ (matrix[0]!)[column.val]! * (matrix[0]!)[column.val]!
      exact mul_self_nonneg _
    have hnonzero : (∑ k : Fin matrix.size,
        (matrix[0]!)[k.val]! * (matrix[0]!)[k.val]! : ZZ) ≠ 0 := by
      intro hzero
      have hall := (Finset.sum_mul_self_eq_zero_iff Finset.univ
        (fun k : Fin matrix.size => (matrix[0]!)[k.val]!)).mp hzero
      obtain ⟨column, hcolumn⟩ := hentry
      exact hcolumn (hall column (Finset.mem_univ column))
    rcases lt_or_eq_of_le hnonnegative with hpositive | hzero
    · exact hpositive
    · exact (hnonzero hzero.symm).elim
  exact Int.cast_pos.mpr hsum

theorem initialProcessedGramSchmidtValid
    {matrix : Generated.StrictRecombine.LLLMatrix}
    (hinput : ConcreteLLLInputValid matrix) (hsize : 0 < matrix.size) :
    let mu := Generated.StrictRecombine.zeroQQMatrix matrix.size matrix.size
    let dot : QQ :=
      ((∑ k : Fin matrix.size,
        (matrix[0]!)[k.val]! * (matrix[0]!)[k.val]! : ZZ) : QQ)
    let norms := (Array.replicate matrix.size (0 : QQ)).set 0 dot (by simp [hsize])
    ProcessedGramSchmidtValid
      (Generated.StrictRecombine.LLLState.mk matrix #[] mu norms 1) 1 := by
  dsimp only
  have hrowZero := hinput.rows_square 0 hsize
  refine ⟨?_, hsize,
    concreteGramSchmidtUpTo_one_initial matrix hsize hrowZero, ?_⟩
  · refine ⟨zeroQQMatrix_size matrix.size matrix.size, ?_, by simp⟩
    intro row hrow
    have hrow' : row < matrix.size := by
      simpa only [zeroQQMatrix_size] using hrow
    have hzeroRow := zeroQQMatrix_rows matrix.size matrix.size row hrow'
    rw [getElem!_pos _ row (by simpa only [zeroQQMatrix_size] using hrow')] at hzeroRow
    exact hzeroRow
  · intro index hindex
    have hindexZero : index = 0 := by omega
    subst index
    rw [getElem!_pos _ 0 (by simp [hsize])]
    rw [Array.getElem_set]
    simp only [if_pos, hrowZero]
    exact hinput.first_norm_positive hsize

theorem ConcreteLLLExecutionValid.toInputValid
    {state : Generated.StrictRecombine.LLLState}
    (hvalid : ConcreteLLLExecutionValid state) :
    ConcreteLLLInputValid state.matrix := by
  refine ⟨hvalid.rows_square, ?_⟩
  intro hdet
  have hgramMatrix := gramPrefixMatrix_eq_mul_transpose state.matrix
    state.matrix.size
  have hgramDet : Matrix.det
      (gramPrefixMatrix state.matrix state.matrix.size) = 0 := by
    rw [hgramMatrix, Matrix.det_mul, Matrix.det_transpose, hdet, zero_mul]
  have hgramDetQQ :
      ((Matrix.det (gramPrefixMatrix state.matrix state.matrix.size) : Int) : QQ) =
        0 := by rw [hgramDet]; norm_num
  rw [hvalid.gram_schmidt.gram_prefix state.matrix.size le_rfl] at hgramDetQQ
  have hpositive : 0 < prefixNormProduct state.norms state.matrix.size := by
    unfold prefixNormProduct
    apply Finset.prod_pos
    intro index _
    have hindex : index.val < state.norms.size := by
      simpa [hvalid.norms_size] using index.isLt
    simpa [getElem!_pos state.norms index.val hindex] using
      hvalid.norms_positive index.val hindex
  linarith

theorem ConcreteLLLExecutionValid.toConcreteLLLValid
    {state : Generated.StrictRecombine.LLLState}
    (hvalid : ConcreteLLLExecutionValid state) : ConcreteLLLValid state where
  norms_size := hvalid.norms_size
  rows_square := hvalid.rows_square
  norms_positive := hvalid.norms_positive
  gram_prefix := hvalid.gram_schmidt.gram_prefix

theorem ConcreteLLLExecutionValid.withK
    {state : Generated.StrictRecombine.LLLState}
    (hvalid : ConcreteLLLExecutionValid state) (k : Nat) :
    ConcreteLLLExecutionValid { state with k := k } where
  norms_size := hvalid.norms_size
  mu_size := hvalid.mu_size
  rows_square := hvalid.rows_square
  mu_rows_square := hvalid.mu_rows_square
  norms_positive := hvalid.norms_positive
  gram_schmidt := by
    intro rowCount hrowCount
    simpa [ConcreteGramSchmidt, gsLowerPrefix, gsNormDiagonal] using
      hvalid.gram_schmidt rowCount hrowCount

theorem ConcreteLLLExecutionValid.withTransform
    {state : Generated.StrictRecombine.LLLState}
    (hvalid : ConcreteLLLExecutionValid state)
    (transform : Generated.StrictRecombine.LLLMatrix) :
    ConcreteLLLExecutionValid { state with transform := transform } where
  norms_size := hvalid.norms_size
  mu_size := hvalid.mu_size
  rows_square := hvalid.rows_square
  mu_rows_square := hvalid.mu_rows_square
  norms_positive := hvalid.norms_positive
  gram_schmidt := by
    intro rowCount hrowCount
    simpa [ConcreteGramSchmidt, gsLowerPrefix, gsNormDiagonal] using
      hvalid.gram_schmidt rowCount hrowCount

/-- The generated μ update is exactly the lower factor required after
swapping adjacent basis rows.  The left side swaps the two old rows; the
right side is the generated new lower factor followed by the concrete local
column basis change. -/
theorem gsLowerPrefix_lovaszSwapMuResult_mul
    (state : Generated.StrictRecombine.LLLState) (k rowCount : Nat)
    (muNew : QQ) (hvalid : ConcreteLLLExecutionValid state)
    (hkPositive : 0 < k) (hkCount : k < rowCount)
    (hrowCount : rowCount ≤ state.matrix.size) :
    (fun row column => gsLowerPrefix state rowCount
      ((Equiv.swap
        (⟨k - 1, by omega⟩ : Fin rowCount)
        (⟨k, hkCount⟩ : Fin rowCount)) row) column) =
      gsLowerPrefix
          { state with mu := (lovaszSwapMuResult state.mu k
              ((state.mu[k]!)[k - 1]!) muNew) }
          rowCount *
        lovaszLocalTransform
          (⟨k - 1, by omega⟩ : Fin rowCount)
          (⟨k, hkCount⟩ : Fin rowCount)
          ((state.mu[k]!)[k - 1]!) muNew := by
  have hpreviousCurrent :
      (⟨k - 1, by omega⟩ : Fin rowCount) ≠ ⟨k, hkCount⟩ := by
    intro heq
    injection heq
    omega
  have hkMu : k < state.mu.size := by
    rw [hvalid.mu_size]
    exact lt_of_lt_of_le hkCount hrowCount
  have hrowsSquare : ∀ index (hindex : index < state.mu.size),
      state.mu[index]!.size = state.mu.size := by
    intro index hindex
    rw [getElem!_pos state.mu index hindex]
    rw [hvalid.mu_rows_square index hindex, hvalid.mu_size]
  funext row column
  rw [mul_lovaszLocalTransform_apply _
    (⟨k - 1, by omega⟩ : Fin rowCount) (⟨k, hkCount⟩ : Fin rowCount)
    ((state.mu[k]!)[k - 1]!) muNew hpreviousCurrent row column]
  have hrowMu : row.val < state.mu.size := by
    rw [hvalid.mu_size]
    exact lt_of_lt_of_le row.isLt hrowCount
  have hcolumnMu : column.val < state.mu.size := by
    rw [hvalid.mu_size]
    exact lt_of_lt_of_le column.isLt hrowCount
  have hentry := lovaszSwapMuResult_entry state.mu k row.val column.val
    ((state.mu[k]!)[k - 1]!) muNew hkPositive hkMu hrowMu hcolumnMu
    hrowsSquare
  have hpredMu : k - 1 < state.mu.size := by omega
  have hentryPrevious := lovaszSwapMuResult_entry state.mu k row.val (k - 1)
    ((state.mu[k]!)[k - 1]!) muNew hkPositive hkMu hrowMu hpredMu
    hrowsSquare
  have hentryCurrent := lovaszSwapMuResult_entry state.mu k row.val k
    ((state.mu[k]!)[k - 1]!) muNew hkPositive hkMu hrowMu hkMu
    hrowsSquare
  unfold gsLowerPrefix
  rw [hentry]
  by_cases hrowPrevious : row = (⟨k - 1, by omega⟩ : Fin rowCount)
  · rw [hrowPrevious] at hentry ⊢
    rw [Equiv.swap_apply_left]
    by_cases hcolumnPrevious :
        column = (⟨k - 1, by omega⟩ : Fin rowCount)
    · rw [hcolumnPrevious] at hentry ⊢
      simp [show ¬ k < k - 1 by omega] at *
      simp_all [Fin.ext_iff] <;> ring
    · by_cases hcolumnCurrent : column = (⟨k, hkCount⟩ : Fin rowCount)
      · rw [hcolumnCurrent] at hentry ⊢
        simp [show ¬ k < k - 1 by omega] at *
        simp_all [Fin.ext_iff] <;> ring
      · simp_all [Fin.ext_iff]
        split_ifs <;> simp_all [Fin.ext_iff] <;> try omega <;> ring
  · by_cases hrowCurrent : row = (⟨k, hkCount⟩ : Fin rowCount)
    · rw [hrowCurrent] at hentry ⊢
      rw [Equiv.swap_apply_right]
      by_cases hcolumnPrevious :
          column = (⟨k - 1, by omega⟩ : Fin rowCount)
      · rw [hcolumnPrevious] at hentry ⊢
        simp [show ¬ k < k - 1 by omega] at *
        simp_all [Fin.ext_iff] <;> ring
      · by_cases hcolumnCurrent : column = (⟨k, hkCount⟩ : Fin rowCount)
        · rw [hcolumnCurrent] at hentry ⊢
          simp [show ¬ k < k - 1 by omega] at *
          simp_all [Fin.ext_iff] <;> ring
        · simp_all [Fin.ext_iff]
          split_ifs <;> simp_all [Fin.ext_iff] <;> try omega <;> ring
    · rw [Equiv.swap_apply_of_ne_of_ne hrowPrevious hrowCurrent]
      have hrowPreviousVal : row.val ≠ k - 1 := by
        intro heq
        apply hrowPrevious
        exact Fin.ext heq
      have hrowCurrentVal : row.val ≠ k := by
        intro heq
        apply hrowCurrent
        exact Fin.ext heq
      have hpredLtRow : k - 1 < row.val ↔ k < row.val := by omega
      have hkPrevious : k ≠ k - 1 := by omega
      by_cases hcolumnPrevious :
          column = (⟨k - 1, by omega⟩ : Fin rowCount)
      · rw [hcolumnPrevious] at hentry ⊢
        by_cases hrowAfter : k < row.val
        · simp [hrowAfter, hkPrevious, hrowPreviousVal, hrowCurrentVal,
            hpredLtRow] at hentry hentryPrevious hentryCurrent ⊢
          rw [hentryPrevious, hentryCurrent]
          ring
        · simp [hrowAfter, hkPrevious, hrowPreviousVal, hrowCurrentVal,
            hpredLtRow, hrowPrevious, hrowCurrent]
      · by_cases hcolumnCurrent : column = (⟨k, hkCount⟩ : Fin rowCount)
        · rw [hcolumnCurrent] at hentry ⊢
          by_cases hrowAfter : k < row.val
          · simp [hrowAfter, hkPrevious, hrowPreviousVal, hrowCurrentVal,
              hpredLtRow] at hentry hentryPrevious hentryCurrent ⊢
            rw [hentryPrevious, hentryCurrent]
            ring
          · simp [hrowAfter, hkPrevious, hrowPreviousVal, hrowCurrentVal,
              hpredLtRow, hrowPrevious, hrowCurrent]
        · by_cases hrowAfter : k < row.val
          · simp [hrowAfter, hkPrevious, hrowPreviousVal, hrowCurrentVal,
              hpredLtRow] at hentry hentryPrevious hentryCurrent ⊢
            split_ifs <;> simp_all [Fin.ext_iff] <;> try omega <;> ring
          · simp [hrowAfter, hkPrevious, hrowPreviousVal, hrowCurrentVal,
              hpredLtRow, hrowPrevious, hrowCurrent, hcolumnPrevious,
              hcolumnCurrent]

theorem gsNormDiagonal_normsAfterLovaszSwap
    (state : Generated.StrictRecombine.LLLState) (k rowCount : Nat)
    (muNew : QQ) (hvalid : ConcreteLLLExecutionValid state)
    (hkPositive : 0 < k) (hkCount : k < rowCount)
    (hrowCount : rowCount ≤ state.matrix.size)
    (hmuNew : muNew =
      ((state.mu[k]!)[k - 1]!) * state.norms[k - 1]! /
        (state.norms[k]! + ((state.mu[k]!)[k - 1]!) *
          ((state.mu[k]!)[k - 1]!) * state.norms[k - 1]!)) :
    lovaszLocalTransform
          (⟨k - 1, by omega⟩ : Fin rowCount)
          (⟨k, hkCount⟩ : Fin rowCount)
          ((state.mu[k]!)[k - 1]!) muNew *
        gsNormDiagonal state rowCount *
        (lovaszLocalTransform
          (⟨k - 1, by omega⟩ : Fin rowCount)
          (⟨k, hkCount⟩ : Fin rowCount)
          ((state.mu[k]!)[k - 1]!) muNew).transpose =
      gsNormDiagonal
        { state with norms := (Generated.StrictRecombine.normsAfterLovaszSwap
            state.norms k ((state.mu[k]!)[k - 1]!)) }
        rowCount := by
  have hkNorm : k < state.norms.size := by
    rw [hvalid.norms_size]
    exact lt_of_lt_of_le hkCount hrowCount
  have hpredNorm : k - 1 < state.norms.size := by omega
  have hdeltaPos : 0 < state.norms[k]! +
      ((state.mu[k]!)[k - 1]!) * ((state.mu[k]!)[k - 1]!) *
        state.norms[k - 1]! := by
    have hkPos := hvalid.norms_positive k hkNorm
    have hpPos := hvalid.norms_positive (k - 1) hpredNorm
    have hkPos' : 0 < state.norms[k]! := by
      simpa only [getElem!_pos state.norms k hkNorm] using hkPos
    have hpPos' : 0 < state.norms[k - 1]! := by
      simpa only [getElem!_pos state.norms (k - 1) hpredNorm] using hpPos
    nlinarith [sq_nonneg ((state.mu[k]!)[k - 1]!)]
  have hlocal := lovaszLocalTransform_diagonal
    (⟨k - 1, by omega⟩ : Fin rowCount)
    (⟨k, hkCount⟩ : Fin rowCount)
    (fun index : Fin rowCount => state.norms[index.val]!)
    state.norms[k - 1]! state.norms[k]!
    ((state.mu[k]!)[k - 1]!) muNew
    (state.norms[k]! + ((state.mu[k]!)[k - 1]!) *
      ((state.mu[k]!)[k - 1]!) * state.norms[k - 1]!)
    (by intro heq; injection heq; omega) rfl rfl rfl
    (ne_of_gt hdeltaPos) hmuNew
  unfold gsNormDiagonal
  rw [hlocal]
  congr 1
  funext index
  have hindexNorm : index.val < state.norms.size := by
    rw [hvalid.norms_size]
    exact lt_of_lt_of_le index.isLt hrowCount
  rw [normsAfterLovaszSwap_get state.norms k
    ((state.mu[k]!)[k - 1]!) hkNorm hpredNorm (by
      have := ne_of_gt hdeltaPos
      simpa only [getElem!_pos state.norms k hkNorm,
        getElem!_pos state.norms (k - 1) hpredNorm] using this)
    index.val hindexNorm]
  have hkNePrevious : k ≠ k - 1 := by omega
  by_cases hindexPrevious : index.val = k - 1
  · simp [hindexPrevious, hkNePrevious, hkNePrevious.symm, Fin.ext_iff,
      getElem!_pos state.norms k hkNorm,
      getElem!_pos state.norms (k - 1) hpredNorm]
  · by_cases hindexCurrent : index.val = k
    · simp [hindexPrevious, hindexCurrent, hkNePrevious,
        hkNePrevious.symm, Fin.ext_iff,
        getElem!_pos state.norms k hkNorm,
        getElem!_pos state.norms (k - 1) hpredNorm]
      ring
    · simp [hindexPrevious, Ne.symm hindexPrevious, hindexCurrent,
        Ne.symm hindexCurrent, hkNePrevious, hkNePrevious.symm, Fin.ext_iff,
        getElem!_pos state.norms index.val hindexNorm]

theorem sizeReduceAt_mu_eq
    (state output : Generated.StrictRecombine.LLLState) (source : Nat)
    (hvalid : ConcreteLLLExecutionValid state)
    (hsourceK : source < state.k)
    (hrun : Generated.StrictRecombine.sizeReduceAt state source = .ok output) :
    ∃ q : ZZ,
      Generated.StrictRecombine.roundQQ
        ((state.mu[state.k]!)[source]!) = .ok q ∧
      output.mu = sizeReduceMuResult state.mu state.k source q := by
  unfold Generated.StrictRecombine.sizeReduceAt at hrun
  split at hrun
  next hk =>
    split at hrun
    next hsource =>
      cases hround : Generated.StrictRecombine.roundQQ
          state.mu[state.k][source] with
      | error fault => simp [hround] at hrun
      | ok q =>
        simp only [hround] at hrun
        split at hrun
        next hzero =>
          have hout := Except.ok.inj hrun
          subst output
          refine ⟨q, ?_, by simp [sizeReduceMuResult, hzero]⟩
          simpa [getElem!_pos state.mu state.k hk,
            getElem!_pos state.mu[state.k] source hsource] using hround
        next hnonzero =>
          cases hsubtract : Generated.StrictRecombine.subtractMatrixRows
              state.matrix state.transform state.k source q with
          | error fault => simp [hsubtract] at hrun
          | ok matrices =>
            rcases matrices with ⟨matrix', transform'⟩
            simp only [hsubtract] at hrun
            cases hmu : Generated.StrictRecombine.reduceMuPrefixLoop
                state.mu state.k source q source 0 with
            | error fault => simp [hmu] at hrun
            | ok mu' =>
              simp only [hmu] at hrun
              split at hrun
              next hk' =>
                split at hrun
                next hsource' =>
                  have hout := Except.ok.inj hrun
                  subst output
                  have hsourceMu : source < state.mu.size := by omega
                  have hlimitK : source ≤ state.mu[state.k].size := by
                    rw [hvalid.mu_rows_square state.k hk]
                    rw [← hvalid.mu_size]
                    omega
                  have hlimitSource : source ≤ state.mu[source].size := by
                    rw [hvalid.mu_rows_square source hsourceMu]
                    rw [← hvalid.mu_size]
                    omega
                  have hmuExact := reduceMuPrefixLoop_output_eq
                    state.mu mu' state.k source q source 0 hk hsourceMu
                    hlimitK hlimitSource (by omega) hmu
                  subst mu'
                  refine ⟨q, ?_, ?_⟩
                  · simpa [getElem!_pos state.mu state.k hk,
                      getElem!_pos state.mu[state.k] source hsource] using hround
                  · simp only [sizeReduceMuResult, if_neg hnonzero]
                    simp [Array.setIfInBounds, hk', hsource']
                next hsource' => contradiction
              next hk' => contradiction
    next hsource => contradiction
  next hk => contradiction

theorem gsLowerPrefix_sizeReduceMuResult
    (state : Generated.StrictRecombine.LLLState)
    (source rowCount : Nat) (q : ZZ)
    (hvalid : ConcreteLLLExecutionValid state)
    (hsourceK : source < state.k)
    (hkCount : state.k < rowCount)
    (hrowCount : rowCount ≤ state.matrix.size) :
    gsLowerPrefix
        { state with mu := sizeReduceMuResult state.mu state.k source q }
        rowCount =
      Matrix.transvection
          (⟨state.k, hkCount⟩ : Fin rowCount)
          (⟨source, lt_trans hsourceK hkCount⟩ : Fin rowCount)
          (-(q : QQ)) * gsLowerPrefix state rowCount := by
  rw [Matrix.transvection, Matrix.add_mul, Matrix.one_mul]
  funext row column
  by_cases hrowK : row = (⟨state.k, hkCount⟩ : Fin rowCount)
  · subst row
    have hkMu : state.k < state.mu.size := by
      rw [hvalid.mu_size]
      exact lt_of_lt_of_le hkCount hrowCount
    have hsourceMu : source < state.mu.size := lt_trans hsourceK hkMu
    have hkRowSize := hvalid.mu_rows_square state.k hkMu
    have hsourceRowSize := hvalid.mu_rows_square source hsourceMu
    have hsourceKRow : source < state.mu[state.k].size := by
      rw [hkRowSize]
      exact lt_of_lt_of_le hsourceK
        (lt_of_lt_of_le hkCount hrowCount).le
    have hlimitSource : source ≤ state.mu[source].size := by
      rw [hsourceRowSize]
      exact (lt_of_lt_of_le hsourceK
        (lt_of_lt_of_le hkCount hrowCount).le).le
    have hcolumnKRow : column.val < state.mu[state.k].size := by
      rw [hkRowSize]
      exact lt_of_lt_of_le column.isLt hrowCount
    have hmuGet := sizeReduceMuResult_target_get state.mu state.k
      source column.val q hkMu hsourceMu hsourceKRow hlimitSource
      hcolumnKRow (by omega)
    unfold gsLowerPrefix
    rw [Matrix.add_apply]
    rw [Matrix.single_mul_apply_same]
    by_cases hcolumnSource : column.val < source
    · have hcolumnK : column.val < state.k := lt_trans hcolumnSource hsourceK
      rw [if_pos hcolumnK, hmuGet]
      simp [hcolumnSource, hcolumnK]
      ring
    · by_cases hcolumnSourceEq : column.val = source
      · have hcolumnK : column.val < state.k := by omega
        rw [if_pos hcolumnK, hmuGet]
        have hsourceFin : column =
            (⟨source, lt_trans hsourceK hkCount⟩ : Fin rowCount) :=
          Fin.ext hcolumnSourceEq
        simp [hcolumnSource, hcolumnSourceEq, hsourceFin, hsourceK]
        ring
      · by_cases hcolumnK : column.val < state.k
        · have hsourceColumn : source < column.val := by omega
          rw [if_pos hcolumnK, hmuGet]
          have hkFinNe :
              (⟨state.k, hkCount⟩ : Fin rowCount) ≠ column := by
            intro heq
            have heqVal := congrArg Fin.val heq
            simp only [Fin.val_mk] at heqVal
            omega
          have hsourceFinNe :
              (⟨source, lt_trans hsourceK hkCount⟩ : Fin rowCount) ≠
                column := by
            intro heq
            have heqVal := congrArg Fin.val heq
            simp only [Fin.val_mk] at heqVal
            omega
          simp [hcolumnSource, hcolumnSourceEq, hcolumnK,
            hkFinNe, hsourceFinNe]
        · by_cases hcolumnKEq : column.val = state.k
          · have hcolumnFin : column =
                (⟨state.k, hkCount⟩ : Fin rowCount) :=
              Fin.ext hcolumnKEq
            rw [hcolumnFin]
            have hnotKSource : ¬ state.k < source := by omega
            have hsourceFinNe :
                (⟨source, lt_trans hsourceK hkCount⟩ : Fin rowCount) ≠
                  ⟨state.k, hkCount⟩ := by
              intro heq
              have heqVal := congrArg Fin.val heq
              simp only [Fin.val_mk] at heqVal
              omega
            simp [hnotKSource, hsourceFinNe]
          · have hkColumn : state.k < column.val := by omega
            have hsourceColumn : source < column.val := lt_trans hsourceK hkColumn
            have hkFinNe :
                (⟨state.k, hkCount⟩ : Fin rowCount) ≠ column := by
              intro heq
              have heqVal := congrArg Fin.val heq
              simp only [Fin.val_mk] at heqVal
              omega
            have hsourceFinNe :
                (⟨source, lt_trans hsourceK hkCount⟩ : Fin rowCount) ≠
                  column := by
              intro heq
              have heqVal := congrArg Fin.val heq
              simp only [Fin.val_mk] at heqVal
              omega
            simp [hcolumnK, hcolumnKEq, Nat.not_lt.mpr hsourceColumn.le,
              hkColumn.ne, hkColumn.ne', hsourceColumn.ne, hsourceColumn.ne',
              hkFinNe, hsourceFinNe]
  · unfold gsLowerPrefix
    rw [Matrix.add_apply]
    have hsingle := Matrix.single_mul_apply_of_ne (-(q : QQ))
      (⟨state.k, hkCount⟩ : Fin rowCount)
      (⟨source, lt_trans hsourceK hkCount⟩ : Fin rowCount)
      row column hrowK (gsLowerPrefix state rowCount)
    unfold gsLowerPrefix at hsingle
    rw [hsingle, add_zero]
    by_cases hcolumnRow : column.val < row.val
    · rw [if_pos hcolumnRow, if_pos hcolumnRow]
      have hrowMu : row.val < state.mu.size := by
        rw [hvalid.mu_size]
        exact lt_of_lt_of_le row.isLt hrowCount
      have hrowNe : row.val ≠ state.k := by
        intro heq
        apply hrowK
        exact Fin.ext heq
      rw [sizeReduceMuResult_get_of_ne state.mu state.k source row.val q
        hrowMu hrowNe]
    · rw [if_neg hcolumnRow, if_neg hcolumnRow]

theorem gsLowerPrefix_sizeReduceMuResult_before
    (state : Generated.StrictRecombine.LLLState)
    (source rowCount : Nat) (q : ZZ)
    (hvalid : ConcreteLLLExecutionValid state)
    (hrowCount : rowCount ≤ state.k)
    (hkMatrix : state.k ≤ state.matrix.size) :
    gsLowerPrefix
        { state with mu := sizeReduceMuResult state.mu state.k source q }
        rowCount = gsLowerPrefix state rowCount := by
  funext row column
  unfold gsLowerPrefix
  by_cases hcolumnRow : column.val < row.val
  · rw [if_pos hcolumnRow, if_pos hcolumnRow]
    have hrowMatrix : row.val < state.matrix.size := by omega
    have hrowMu : row.val < state.mu.size := by
      rw [hvalid.mu_size]
      exact hrowMatrix
    have hrowNe : row.val ≠ state.k := by omega
    rw [sizeReduceMuResult_get_of_ne state.mu state.k source row.val q
      hrowMu hrowNe]
  · rw [if_neg hcolumnRow, if_neg hcolumnRow]

theorem ConcreteLLLValid.withK
    {state : Generated.StrictRecombine.LLLState}
    (hvalid : ConcreteLLLValid state) (k : Nat) :
    ConcreteLLLValid { state with k := k } where
  norms_size := hvalid.norms_size
  rows_square := hvalid.rows_square
  norms_positive := hvalid.norms_positive
  gram_prefix := hvalid.gram_prefix

theorem prefixNormProduct_pos
    (norms : Array QQ) (rowCount : Nat)
    (hrowCount : rowCount ≤ norms.size)
    (hpositive : ∀ index (hindex : index < norms.size), 0 < norms[index]) :
    0 < prefixNormProduct norms rowCount := by
  unfold prefixNormProduct
  apply Finset.prod_pos
  intro i _
  rw [getElem!_pos norms i.val (lt_of_lt_of_le i.isLt hrowCount)]
  exact hpositive i.val (lt_of_lt_of_le i.isLt hrowCount)

theorem prefixNormProduct_succ (norms : Array QQ) (rowCount : Nat) :
    prefixNormProduct norms (rowCount + 1) =
      prefixNormProduct norms rowCount * norms[rowCount]! := by
  unfold prefixNormProduct
  rw [Fin.prod_univ_castSucc]
  rfl

theorem generatedGramPivot_positive
    (state : Generated.StrictRecombine.LLLState) (i : Nat)
    (hinput : ConcreteLLLInputValid state.matrix)
    (hi : i < state.matrix.size)
    (hnormsSize : state.norms.size = state.matrix.size)
    (hgs : ConcreteGramSchmidtUpTo state (i + 1))
    (hprevious : ∀ index, index < i → 0 < state.norms[index]!) :
    0 < state.norms[i]! := by
  have hdetPositive := hinput.gramPrefixMatrixQQ_det_pos (i + 1) (by omega)
  have hproduct := hgs.gram_prefix (i + 1) le_rfl (by omega)
  have hproductEq : Matrix.det
      (gramPrefixMatrixQQ state.matrix (i + 1)) =
      prefixNormProduct state.norms (i + 1) := by
    calc
      Matrix.det (gramPrefixMatrixQQ state.matrix (i + 1)) =
          ((Matrix.det (gramPrefixMatrix state.matrix (i + 1)) : Int) : QQ) := by
        exact (Int.cast_det _).symm
      _ = prefixNormProduct state.norms (i + 1) := hproduct
  have hproductPositive : 0 < prefixNormProduct state.norms (i + 1) := by
    rw [← hproductEq]
    exact hdetPositive
  have hpreviousProduct : 0 < prefixNormProduct state.norms i := by
    unfold prefixNormProduct
    apply Finset.prod_pos
    intro index _
    exact hprevious index.val index.isLt
  rw [prefixNormProduct_succ] at hproductPositive
  nlinarith

theorem processedGramSchmidtValid_step
    (matrix : Generated.StrictRecombine.LLLMatrix)
    (mu : Generated.StrictRecombine.QQMatrix) (norms : Array QQ)
    (i : Nat) (muResult : Generated.StrictRecombine.QQMatrix × Array QQ)
    (norm : QQ)
    (hinput : ConcreteLLLInputValid matrix)
    (hvalid : ProcessedGramSchmidtValid
      (Generated.StrictRecombine.LLLState.mk matrix #[] mu norms 1) i)
    (hi : i < matrix.size)
    (hmuRun : Generated.StrictRecombine.gramMuRowLoop matrix i 0 mu norms =
      .ok muResult)
    (hiResultNorm : i < muResult.2.size)
    (hnormRun : Generated.StrictRecombine.gramNormLoop
      muResult.1 muResult.2 i 0
        (((∑ k : Fin matrix[i].size,
          matrix[i][k.val] * matrix[i][k.val] : ZZ) : QQ)) = .ok norm) :
    let norms' := muResult.2.set i norm hiResultNorm
    ProcessedGramSchmidtValid
      (Generated.StrictRecombine.LLLState.mk matrix #[] muResult.1 norms' 1)
      (i + 1) := by
  dsimp only
  have hmatrixRows := hinput.rows_square
  obtain ⟨known, hknown, hresultShape, hresultNorms, hresultRows,
      hresultPrefix⟩ :=
    gramMuRowLoop_succeeds matrix i 0 matrix.size mu norms rfl
      hmatrixRows hvalid.shape hi (by omega)
  rw [hmuRun] at hknown
  have hknownEq := Except.ok.inj hknown
  subst known
  have hiNorm : i < muResult.2.size := by
    exact hiResultNorm
  let norms' := muResult.2.set i norm hiNorm
  have hshape' : GramStorageShape muResult.1 norms' matrix.size :=
    hresultShape.setNorm i norm hiNorm
  have hrowCorrect : GramMuPrefixCorrect matrix muResult.1 muResult.2 i i := by
    exact gramMuRowLoop_prefix_correct matrix i 0 matrix.size mu norms
      muResult rfl hmatrixRows hvalid.shape hi (by omega)
      (fun column hcolumn => by
        exact hvalid.norms_positive column hcolumn)
      (gramMuPrefixCorrect_zero matrix mu norms i) hmuRun
  have hprefixAfterMu : ConcreteGramSchmidtUpTo
      (Generated.StrictRecombine.LLLState.mk matrix #[] muResult.1
        muResult.2 1) i := by
    apply concreteGramSchmidtUpTo_of_prefix_frames matrix #[] muResult.1 mu
      muResult.2 norms 1 i hvalid.gram_schmidt
    · intro row hrow
      have hrowMu : row < mu.size := by
        rw [hvalid.shape.mu_size]
        exact lt_trans hrow hi
      rw [hresultRows row hrowMu (by omega)]
      exact (getElem!_pos mu row hrowMu).symm
    · intro index hindex
      rw [hresultNorms]
  have hprefixFinal : ConcreteGramSchmidtUpTo
      (Generated.StrictRecombine.LLLState.mk matrix #[] muResult.1 norms' 1) i := by
    apply concreteGramSchmidtUpTo_of_prefix_frames matrix #[] muResult.1
      muResult.1 norms' muResult.2 1 i hprefixAfterMu
    · intro row hrow
      rfl
    · intro index hindex
      rw [getElem!_pos norms' index (by simpa [norms'] using
        (lt_trans hindex hiNorm))]
      rw [Array.getElem_set_ne hiNorm (lt_trans hindex hiNorm) (by omega)]
      exact (getElem!_pos muResult.2 index (lt_trans hindex hiNorm)).symm
  have hrowFinal : GramMuPrefixCorrect matrix muResult.1 norms' i i := by
    exact gramMuPrefixCorrect_set_norm_at_end matrix muResult.1 muResult.2
      i i norm le_rfl hiNorm hrowCorrect
  have hdiagonalOld := gramNormLoop_closes_diagonal matrix muResult.1
    muResult.2 i matrix.size norm rfl hmatrixRows hresultShape hi hnormRun
  have hdiagonalFinal : sourceRowDot matrix i i matrix[i]!.size =
      (∑ index : Fin i,
        (muResult.1[i]!)[index.val]! * (muResult.1[i]!)[index.val]! *
          norms'[index.val]!) + norms'[i]! := by
    rw [hdiagonalOld]
    congr 1
    · apply Fintype.sum_congr
      intro index
      have hindex : index.val < muResult.2.size := lt_trans index.isLt hiNorm
      rw [getElem!_pos norms' index.val (by simpa [norms'] using hindex)]
      rw [Array.getElem_set_ne hiNorm hindex (by omega)]
      rw [getElem!_pos muResult.2 index.val hindex]
    · rw [getElem!_pos norms' i (by simpa [norms'] using hiNorm)]
      simp [norms']
  have hnormsOld : ∀ index, index < i →
      norms'[index]! = norms[index]! := by
    intro index hindex
    have hindexResult : index < muResult.2.size := lt_trans hindex hiNorm
    rw [getElem!_pos norms' index (by simpa [norms'] using hindexResult)]
    rw [Array.getElem_set_ne hiNorm hindexResult (by omega)]
    have hindexNorms : index < norms.size := by
      rw [hvalid.shape.norms_size]
      exact lt_trans hindex hi
    calc
      muResult.2[index] = muResult.2[index]! :=
        (getElem!_pos muResult.2 index hindexResult).symm
      _ = norms[index]! := congrArg (fun values : Array QQ => values[index]!)
        hresultNorms
  have hgsFinal : ConcreteGramSchmidtUpTo
      (Generated.StrictRecombine.LLLState.mk matrix #[] muResult.1 norms' 1)
      (i + 1) := by
    exact concreteGramSchmidtUpTo_extend_one
      (Generated.StrictRecombine.LLLState.mk matrix #[] muResult.1 norms' 1)
      i hi hmatrixRows hprefixFinal hrowFinal hdiagonalFinal
  have hpivot : 0 < norms'[i]! := by
    apply generatedGramPivot_positive
      (Generated.StrictRecombine.LLLState.mk matrix #[] muResult.1 norms' 1)
      i hinput hi hshape'.norms_size hgsFinal
    intro index hindex
    rw [hnormsOld index hindex]
    exact hvalid.norms_positive index hindex
  refine ⟨hshape', Nat.succ_le_iff.mpr hi, hgsFinal, ?_⟩
  intro index hindex
  by_cases hcurrent : index = i
  · subst index
    exact hpivot
  · have hbefore : index < i := by omega
    rw [hnormsOld index hbefore]
    exact hvalid.norms_positive index hbefore

theorem initializeGramSchmidtLoop_processed_valid
    (matrix : Generated.StrictRecombine.LLLMatrix)
    (i : Nat) (mu : Generated.StrictRecombine.QQMatrix) (norms : Array QQ)
    (output : Generated.StrictRecombine.QQMatrix × Array QQ)
    (hinput : ConcreteLLLInputValid matrix)
    (hvalid : ProcessedGramSchmidtValid
      (Generated.StrictRecombine.LLLState.mk matrix #[] mu norms 1) i)
    (hrun : Generated.StrictRecombine.initializeGramSchmidtLoop matrix i
      mu norms = .ok output) :
    ProcessedGramSchmidtValid
      (Generated.StrictRecombine.LLLState.mk matrix #[] output.1 output.2 1)
      matrix.size := by
  induction hremaining : matrix.size - i using Nat.strong_induction_on
      generalizing i mu norms output with
  | h remaining ih =>
      rw [Generated.StrictRecombine.initializeGramSchmidtLoop] at hrun
      by_cases hi : i < matrix.size
      · rw [dif_pos hi] at hrun
        cases hmu : Generated.StrictRecombine.gramMuRowLoop matrix i 0 mu norms with
        | error fault => simp [hmu] at hrun
        | ok muResult =>
          simp only [hmu] at hrun
          have hrowSize := hinput.rows_square i hi
          rw [dotRows_eq_fin_sum matrix[i] matrix[i] (by omega)] at hrun
          simp only at hrun
          cases hnorm : Generated.StrictRecombine.gramNormLoop
              muResult.1 muResult.2 i 0
              (((∑ k : Fin matrix[i].size,
                matrix[i][k.val] * matrix[i][k.val] : ZZ) : QQ)) with
          | error fault =>
            rw [hnorm] at hrun
            simp only at hrun
            contradiction
          | ok norm =>
            rw [hnorm] at hrun
            simp only at hrun
            by_cases hiNorm : i < muResult.2.size
            · rw [dif_pos hiNorm] at hrun
              let norms' := muResult.2.set i norm hiNorm
              have hnext := processedGramSchmidtValid_step matrix mu norms i
                muResult norm hinput hvalid hi hmu hiNorm hnorm
              have hdecrease : matrix.size - (i + 1) < remaining := by omega
              exact ih (matrix.size - (i + 1)) hdecrease (i + 1)
                muResult.1 norms' output (by simpa [norms'] using hnext)
                hrun rfl
            · rw [dif_neg hiNorm] at hrun
              contradiction
      · rw [dif_neg hi] at hrun
        have houtput := Except.ok.inj hrun
        subst output
        have hiEq : i = matrix.size := by
          have hle : i ≤ matrix.size := by simpa using hvalid.processed_le
          omega
        simpa [hiEq] using hvalid

theorem weightedNormPotential_succ (norms : Array QQ) (dimension : Nat) :
    weightedNormPotential norms (dimension + 1) =
      weightedNormPotential norms dimension *
        prefixNormProduct norms (dimension + 1) := by
  unfold weightedNormPotential
  rw [Fin.prod_univ_castSucc]
  rw [prefixNormProduct_succ]
  simp only [Fin.val_castSucc, Fin.val_last]
  have hlast : dimension + 1 - dimension = 1 := by omega
  rw [hlast, pow_one]
  change (∏ i : Fin dimension,
      norms[i.val]! ^ (dimension + 1 - i.val)) * norms[dimension]! =
    (∏ i : Fin dimension, norms[i.val]! ^ (dimension - i.val)) *
      (prefixNormProduct norms dimension * norms[dimension]!)
  have hpowers : (∏ i : Fin dimension,
      norms[i.val]! ^ (dimension + 1 - i.val)) =
      (∏ i : Fin dimension,
        norms[i.val]! ^ (dimension - i.val) * norms[i.val]!) := by
    apply Finset.prod_congr rfl
    intro i _
    have hexponent : dimension + 1 - i.val = (dimension - i.val) + 1 := by
      omega
    rw [hexponent, pow_succ]
  rw [hpowers, Finset.prod_mul_distrib]
  unfold prefixNormProduct
  ring

theorem prefixProductPotential_eq_weightedNormPotential
    (norms : Array QQ) (dimension : Nat) :
    prefixProductPotential norms dimension =
      weightedNormPotential norms dimension := by
  induction dimension with
  | zero => rfl
  | succ dimension ih =>
    unfold prefixProductPotential
    rw [Fin.prod_univ_castSucc]
    change prefixProductPotential norms dimension *
        prefixNormProduct norms (dimension + 1) = _
    rw [ih, weightedNormPotential_succ]

theorem weightedNormPotential_size_eq_arrayLLLPotential (norms : Array QQ) :
    weightedNormPotential norms norms.size = arrayLLLPotential norms := by
  unfold weightedNormPotential arrayLLLPotential
  apply Finset.prod_congr rfl
  intro i _
  rw [getElem!_pos norms i.val i.isLt]
  congr

theorem prefixNormProduct_normsAfterLovaszSwap
    (norms : Array QQ) (k : Nat) (mu : QQ)
    (hk : k < norms.size) (hpred : k - 1 < norms.size)
    (hpositive : 0 < k)
    (hnew : norms[k] + mu * mu * norms[k - 1] ≠ 0) :
    prefixNormProduct
        (Generated.StrictRecombine.normsAfterLovaszSwap norms k mu) k =
      prefixNormProduct norms (k - 1) *
        (norms[k] + mu * mu * norms[k - 1]) := by
  have hkSplit : k - 1 + 1 = k := Nat.sub_add_cancel (by omega)
  calc
    prefixNormProduct
        (Generated.StrictRecombine.normsAfterLovaszSwap norms k mu) k =
        prefixNormProduct
          (Generated.StrictRecombine.normsAfterLovaszSwap norms k mu) (k - 1) *
          (Generated.StrictRecombine.normsAfterLovaszSwap norms k mu)[k - 1]! := by
      simpa only [hkSplit] using
        (prefixNormProduct_succ
          (Generated.StrictRecombine.normsAfterLovaszSwap norms k mu) (k - 1))
    _ = prefixNormProduct norms (k - 1) *
        (norms[k] + mu * mu * norms[k - 1]) := by
      congr 1
      · unfold prefixNormProduct
        apply Finset.prod_congr rfl
        intro i _
        rw [normsAfterLovaszSwap_get norms k mu hk hpred hnew
          i.val (lt_trans i.isLt hpred)]
        have hpredNe : k - 1 ≠ i.val := by omega
        have hkNe : k ≠ i.val := by omega
        simp [hpredNe, hkNe, getElem!_pos norms i.val (lt_trans i.isLt hpred)]
      · rw [normsAfterLovaszSwap_get norms k mu hk hpred hnew
          (k - 1) hpred]
        simp

theorem prefixNormProduct_normsAfterLovaszSwap_lt
    (norms : Array QQ) (k : Nat) (mu : QQ)
    (hk : k < norms.size) (hpred : k - 1 < norms.size)
    (hpositiveK : 0 < k)
    (hnormsPositive : ∀ index (hindex : index < norms.size), 0 < norms[index])
    (hfail : norms[k] < ((3 : QQ) / 4 - mu * mu) * norms[k - 1]) :
    prefixNormProduct
        (Generated.StrictRecombine.normsAfterLovaszSwap norms k mu) k <
      prefixNormProduct norms k := by
  have hkNormPos := hnormsPositive k hk
  have hpredNormPos := hnormsPositive (k - 1) hpred
  have hnewPos : 0 < norms[k] + mu * mu * norms[k - 1] := by
    nlinarith [sq_nonneg mu]
  have hnewLt : norms[k] + mu * mu * norms[k - 1] < norms[k - 1] := by
    nlinarith
  have hprefixPos : 0 < prefixNormProduct norms (k - 1) :=
    prefixNormProduct_pos norms (k - 1) (Nat.le_of_lt hpred) hnormsPositive
  have hkSplit : k - 1 + 1 = k := Nat.sub_add_cancel (by omega)
  rw [prefixNormProduct_normsAfterLovaszSwap norms k mu hk hpred hpositiveK
    (ne_of_gt hnewPos)]
  calc
    prefixNormProduct norms (k - 1) *
        (norms[k] + mu * mu * norms[k - 1]) <
        prefixNormProduct norms (k - 1) * norms[k - 1] :=
      mul_lt_mul_of_pos_left hnewLt hprefixPos
    _ = prefixNormProduct norms k := by
      simpa only [hkSplit, getElem!_pos norms (k - 1) hpred] using
        (prefixNormProduct_succ norms (k - 1)).symm

theorem lllStep_swapped_prefixNormProduct_lt
    (state output : Generated.StrictRecombine.LLLState)
    (hvalid : ConcreteLLLValid state)
    (hrun : Generated.StrictRecombine.lllStep state =
      .ok (.swapped output)) :
    prefixNormProduct output.norms state.k <
      prefixNormProduct state.norms state.k := by
  rw [Generated.StrictRecombine.lllStep] at hrun
  by_cases hkPositive : 0 < state.k
  · rw [dif_pos hkPositive] at hrun
    by_cases hkMatrix : state.k < state.matrix.size
    · rw [dif_pos hkMatrix] at hrun
      cases hreduce : Generated.StrictRecombine.sizeReduceAt state
          (state.k - 1) with
      | error fault => simp [hreduce] at hrun
      | ok reduced =>
        simp only [hreduce] at hrun
        have hcontrol := sizeReduceAt_preserves_norms_k state reduced
          (state.k - 1) hreduce
        repeat' first | split at hrun | simp_all
        all_goals cases hrun
        all_goals
          have hreducedPositive : ∀ index (hindex : index < reduced.norms.size),
              0 < reduced.norms[index] := by
            intro index hindex
            have hindexState : index < state.norms.size := by
              simpa only [hcontrol.1] using hindex
            have hpositiveState := hvalid.norms_positive index hindexState
            simpa only [← hcontrol.1] using hpositiveState
          have hkStateNorm : state.k < state.norms.size := by
            simpa only [← hcontrol.1, ← hcontrol.2] using
              (‹reduced.k < reduced.norms.size›)
          have hpredStateNorm : state.k - 1 < state.norms.size := by
            simpa only [← hcontrol.1, ← hcontrol.2] using
              (‹reduced.k - 1 < reduced.norms.size›)
          have hkStateMu : state.k < reduced.mu.size := by
            simpa only [← hcontrol.2] using (‹reduced.k < reduced.mu.size›)
          have hpredStateMu : state.k - 1 < reduced.mu[state.k].size := by
            simpa only [← hcontrol.2] using
              (‹reduced.k - 1 < reduced.mu[reduced.k].size›)
          have hfail : reduced.norms[reduced.k] <
              ((3 : QQ) / 4 -
                reduced.mu[reduced.k][reduced.k - 1] *
                  reduced.mu[reduced.k][reduced.k - 1]) *
                reduced.norms[reduced.k - 1] := by
            have hfailState :
                state.norms[state.k]'hkStateNorm <
                  ((3 : QQ) / 4 -
                    (reduced.mu[state.k]'hkStateMu)[state.k - 1]'hpredStateMu *
                      (reduced.mu[state.k]'hkStateMu)[state.k - 1]'hpredStateMu) *
                    state.norms[state.k - 1]'hpredStateNorm := by assumption
            simpa only [← hcontrol.1, ← hcontrol.2] using hfailState
          have hdecrease := prefixNormProduct_normsAfterLovaszSwap_lt
            reduced.norms reduced.k reduced.mu[reduced.k][reduced.k - 1]
            (by assumption) (by assumption) (by simpa [hcontrol.2] using hkPositive)
            hreducedPositive hfail
          simpa only [hcontrol.1, hcontrol.2] using hdecrease
    · rw [dif_neg hkMatrix] at hrun
      simp at hrun
  · rw [dif_neg hkPositive] at hrun
    simp at hrun

theorem lllStep_swapped_source_index_lt
    (state output : Generated.StrictRecombine.LLLState)
    (hrun : Generated.StrictRecombine.lllStep state =
      .ok (.swapped output)) :
    state.k < state.matrix.size := by
  rw [Generated.StrictRecombine.lllStep] at hrun
  by_cases hkPositive : 0 < state.k
  · rw [dif_pos hkPositive] at hrun
    by_cases hkMatrix : state.k < state.matrix.size
    · exact hkMatrix
    · rw [dif_neg hkMatrix] at hrun
      contradiction
  · rw [dif_neg hkPositive] at hrun
    contradiction

theorem lllStep_swapped_matrix_size
    (state output : Generated.StrictRecombine.LLLState)
    (hvalid : ConcreteLLLValid state)
    (houtputValid : ConcreteLLLValid output)
    (hrun : Generated.StrictRecombine.lllStep state =
      .ok (.swapped output)) :
    output.matrix.size = state.matrix.size := by
  rcases lllStep_swapped_norms state output hrun with
    ⟨reduced, hreduce, hnorms, hk, hout⟩
  calc
    output.matrix.size = output.norms.size := houtputValid.norms_size.symm
    _ = reduced.norms.size := by
      rw [hout, normsAfterLovaszSwap_size]
    _ = state.norms.size := congrArg Array.size hnorms
    _ = state.matrix.size := hvalid.norms_size

theorem lllStep_swapped_boundaryGramDet_lt
    (state output : Generated.StrictRecombine.LLLState)
    (hvalid : ConcreteLLLValid state)
    (houtputValid : ConcreteLLLValid output)
    (hrun : Generated.StrictRecombine.lllStep state =
      .ok (.swapped output)) :
    (Matrix.det (gramPrefixMatrix output.matrix state.k)).natAbs <
      (Matrix.det (gramPrefixMatrix state.matrix state.k)).natAbs := by
  have hkMatrix := lllStep_swapped_source_index_lt state output hrun
  have hsize := lllStep_swapped_matrix_size state output hvalid houtputValid hrun
  have hkOutput : state.k ≤ output.matrix.size := by omega
  have hprefix := lllStep_swapped_prefixNormProduct_lt state output hvalid hrun
  have houtGram := houtputValid.gram_prefix state.k hkOutput
  have hinGram := hvalid.gram_prefix state.k (Nat.le_of_lt hkMatrix)
  rw [← houtGram, ← hinGram] at hprefix
  have hdet : Matrix.det (gramPrefixMatrix output.matrix state.k) <
      Matrix.det (gramPrefixMatrix state.matrix state.k) := by
    exact_mod_cast hprefix
  have houtProductPos := prefixNormProduct_pos output.norms state.k
    (by rw [houtputValid.norms_size]; exact hkOutput)
    houtputValid.norms_positive
  have houtDetPos : 0 < Matrix.det
      (gramPrefixMatrix output.matrix state.k) := by
    exact_mod_cast (show (0 : QQ) <
      ((Matrix.det (gramPrefixMatrix output.matrix state.k) : Int) : QQ) by
        rw [houtGram]
        exact houtProductPos)
  exact Int.natAbs_lt_natAbs_of_nonneg_of_lt
    (Int.le_of_lt houtDetPos) hdet

theorem ConcreteLLLValid.gram_prefix_pos
    {state : Generated.StrictRecombine.LLLState}
    (hvalid : ConcreteLLLValid state) (rowCount : Nat)
    (hrowCount : rowCount ≤ state.matrix.size) :
    0 < Matrix.det (gramPrefixMatrix state.matrix rowCount) := by
  have hproduct := prefixNormProduct_pos state.norms rowCount
    (by rw [hvalid.norms_size]; exact hrowCount) hvalid.norms_positive
  have heq := hvalid.gram_prefix rowCount hrowCount
  exact_mod_cast (show (0 : QQ) <
    ((Matrix.det (gramPrefixMatrix state.matrix rowCount) : Int) : QQ) by
      rw [heq]
      exact hproduct)

theorem ConcreteLLLValid.gram_prefix_natAbs
    {state : Generated.StrictRecombine.LLLState}
    (hvalid : ConcreteLLLValid state) (rowCount : Nat)
    (hrowCount : rowCount ≤ state.matrix.size) :
    ((Matrix.det (gramPrefixMatrix state.matrix rowCount)).natAbs : Int) =
      Matrix.det (gramPrefixMatrix state.matrix rowCount) :=
  Int.natAbs_of_nonneg
    (Int.le_of_lt (hvalid.gram_prefix_pos rowCount hrowCount))

theorem ConcreteLLLValid.determinantPotential_cast
    {state : Generated.StrictRecombine.LLLState}
    (hvalid : ConcreteLLLValid state) :
    (lllDeterminantPotential state.matrix : QQ) =
      prefixProductPotential state.norms state.matrix.size := by
  unfold lllDeterminantPotential prefixProductPotential
  push_cast
  apply Finset.prod_congr rfl
  intro i _
  have hrowCount : i.val + 1 ≤ state.matrix.size := i.isLt
  rw [← hvalid.gram_prefix (i.val + 1) hrowCount]
  have hnatAbs := hvalid.gram_prefix_natAbs (i.val + 1) hrowCount
  have hcast := congrArg (fun value : Int => (value : QQ)) hnatAbs
  simpa using hcast

theorem ConcreteLLLValid.determinantPotential_cast_weighted
    {state : Generated.StrictRecombine.LLLState}
    (hvalid : ConcreteLLLValid state) :
    (lllDeterminantPotential state.matrix : QQ) =
      weightedNormPotential state.norms state.matrix.size := by
  rw [hvalid.determinantPotential_cast,
    prefixProductPotential_eq_weightedNormPotential]

theorem lllStep_advanced_determinantPotential_eq
    (state output : Generated.StrictRecombine.LLLState)
    (hvalid : ConcreteLLLValid state)
    (houtputValid : ConcreteLLLValid output)
    (hrun : Generated.StrictRecombine.lllStep state =
      .ok (.advanced output)) :
    lllDeterminantPotential output.matrix =
      lllDeterminantPotential state.matrix := by
  have hcontrol := lllStep_advanced_control state output hrun
  have hsize : output.matrix.size = state.matrix.size := by
    calc
      output.matrix.size = output.norms.size := houtputValid.norms_size.symm
      _ = state.norms.size := congrArg Array.size hcontrol.1
      _ = state.matrix.size := hvalid.norms_size
  have hcastOut := houtputValid.determinantPotential_cast_weighted
  have hcastIn := hvalid.determinantPotential_cast_weighted
  have hcast : (lllDeterminantPotential output.matrix : QQ) =
      (lllDeterminantPotential state.matrix : QQ) := by
    rw [hcastOut, hcastIn, hsize, hcontrol.1]
  exact_mod_cast hcast

theorem weightedNormPotential_normsAfterLovaszSwap_eq
    (state : Generated.StrictRecombine.LLLState) (k : Nat)
    (hkNorm : k < state.norms.size)
    (hpredNorm : k - 1 < state.norms.size)
    (hkMu : k < state.mu.size)
    (hpredMu : k - 1 < state.mu[k].size)
    (hkPositive : 0 < k)
    (hnew : state.norms[k] +
      state.mu[k][k - 1] * state.mu[k][k - 1] * state.norms[k - 1] ≠ 0) :
    weightedNormPotential
        (Generated.StrictRecombine.normsAfterLovaszSwap state.norms k
          state.mu[k][k - 1]) state.norms.size =
      lllPotential' (lllSwapStep' (toPotentialState state)
        (⟨k, hkNorm⟩ : Fin state.norms.size)
        (⟨k - 1, hpredNorm⟩ : Fin state.norms.size)) := by
  unfold weightedNormPotential lllPotential'
  apply Finset.prod_congr rfl
  intro i _
  congr 1
  rw [normsAfterLovaszSwap_get state.norms k state.mu[k][k - 1]
    hkNorm hpredNorm hnew i.val i.isLt]
  simp only [lllSwapStep', toPotentialState]
  rw [getElem!_pos state.mu k hkMu,
    getElem!_pos state.mu[k] (k - 1) hpredMu]
  by_cases hpredI : (⟨k - 1, hpredNorm⟩ : Fin state.norms.size) = i
  · subst i
    simp [pow_two]
  · by_cases hkI : (⟨k, hkNorm⟩ : Fin state.norms.size) = i
    · subst i
      have hkPred : k ≠ k - 1 := by omega
      have hpredK : k - 1 ≠ k := Ne.symm hkPred
      simp [hkPred, hpredK, pow_two]
      ring
    · have hpredVal : k - 1 ≠ i.val := by
        intro heq
        apply hpredI
        exact Fin.ext heq
      have hkVal : k ≠ i.val := by
        intro heq
        apply hkI
        exact Fin.ext heq
      have hiPred : i ≠ (⟨k - 1, hpredNorm⟩ : Fin state.norms.size) :=
        Ne.symm hpredI
      have hiK : i ≠ (⟨k, hkNorm⟩ : Fin state.norms.size) :=
        Ne.symm hkI
      simp [hpredI, hkI, hiPred, hiK, hpredVal, hkVal,
        getElem!_pos state.norms i.val i.isLt]

theorem lllStep_swapped_weightedNormPotential_lt
    (state output : Generated.StrictRecombine.LLLState)
    (hvalid : ConcreteLLLValid state)
    (hrun : Generated.StrictRecombine.lllStep state =
      .ok (.swapped output)) :
    weightedNormPotential output.norms output.norms.size <
      weightedNormPotential state.norms state.norms.size := by
  rw [Generated.StrictRecombine.lllStep] at hrun
  by_cases hkPositive : 0 < state.k
  · rw [dif_pos hkPositive] at hrun
    by_cases hkMatrix : state.k < state.matrix.size
    · rw [dif_pos hkMatrix] at hrun
      cases hreduce : Generated.StrictRecombine.sizeReduceAt state
          (state.k - 1) with
      | error fault => simp [hreduce] at hrun
      | ok reduced =>
        simp only [hreduce] at hrun
        have hcontrol := sizeReduceAt_preserves_norms_k state reduced
          (state.k - 1) hreduce
        repeat' first | split at hrun | simp_all
        all_goals cases hrun
        all_goals
          have hreducedPositive : ∀ index (hindex : index < reduced.norms.size),
              0 < reduced.norms[index] := by
            intro index hindex
            have hindexState : index < state.norms.size := by
              simpa only [hcontrol.1] using hindex
            have hpositiveState := hvalid.norms_positive index hindexState
            simpa only [← hcontrol.1] using hpositiveState
          have hfail : reduced.norms[reduced.k] <
              ((3 : QQ) / 4 -
                reduced.mu[reduced.k][reduced.k - 1] ^ 2) *
                reduced.norms[reduced.k - 1] := by
            have hkStateNorm : state.k < state.norms.size := by
              simpa only [← hcontrol.1, ← hcontrol.2] using
                (‹reduced.k < reduced.norms.size›)
            have hpredStateNorm : state.k - 1 < state.norms.size := by
              simpa only [← hcontrol.1, ← hcontrol.2] using
                (‹reduced.k - 1 < reduced.norms.size›)
            have hkStateMu : state.k < reduced.mu.size := by
              simpa only [← hcontrol.2] using (‹reduced.k < reduced.mu.size›)
            have hpredStateMu : state.k - 1 < reduced.mu[state.k].size := by
              simpa only [← hcontrol.2] using
                (‹reduced.k - 1 < reduced.mu[reduced.k].size›)
            have hfailState :
                state.norms[state.k]'hkStateNorm <
                  ((3 : QQ) / 4 -
                    (reduced.mu[state.k]'hkStateMu)[state.k - 1]'hpredStateMu *
                      (reduced.mu[state.k]'hkStateMu)[state.k - 1]'hpredStateMu) *
                    state.norms[state.k - 1]'hpredStateNorm := by assumption
            simpa only [pow_two, ← hcontrol.1, ← hcontrol.2] using hfailState
          let ki : Fin reduced.norms.size := ⟨reduced.k, by assumption⟩
          let kp : Fin reduced.norms.size := ⟨reduced.k - 1, by assumption⟩
          have hne : ki ≠ kp := by
            intro heq
            have := Fin.ext_iff.mp heq
            dsimp [ki, kp] at this
            omega
          have hadj : ki.val = kp.val + 1 := by
            dsimp [ki, kp]
            omega
          have hfailAbstract : (toPotentialState reduced).gs_norm_sq ki <
              ((3 : QQ) / 4 - (toPotentialState reduced).mu ki kp ^ 2) *
                (toPotentialState reduced).gs_norm_sq kp := by
            dsimp [toPotentialState, ki, kp]
            rw [getElem!_pos reduced.mu reduced.k (by assumption),
              getElem!_pos reduced.mu[reduced.k] (reduced.k - 1) (by assumption)]
            exact hfail
          have habstract := lll_swap_potential_decrease'
            (toPotentialState reduced) ki kp hfailAbstract hne hadj
            (hreducedPositive kp.val kp.isLt)
            (hreducedPositive ki.val ki.isLt)
            (fun i => hreducedPositive i.val i.isLt)
          have hnewPos : 0 < reduced.norms[reduced.k] +
              reduced.mu[reduced.k][reduced.k - 1] *
                reduced.mu[reduced.k][reduced.k - 1] *
                  reduced.norms[reduced.k - 1] := by
            have hkNormPos := hreducedPositive reduced.k (by assumption)
            have hpredNormPos := hreducedPositive (reduced.k - 1) (by assumption)
            nlinarith [sq_nonneg (reduced.mu[reduced.k][reduced.k - 1])]
          have hhelper := weightedNormPotential_normsAfterLovaszSwap_eq
            reduced reduced.k (by assumption) (by assumption)
            (by assumption) (by assumption) (by simpa [hcontrol.2] using hkPositive)
            (ne_of_gt hnewPos)
          have hweighted : weightedNormPotential
                (Generated.StrictRecombine.normsAfterLovaszSwap reduced.norms
                  reduced.k reduced.mu[reduced.k][reduced.k - 1])
                reduced.norms.size <
              weightedNormPotential reduced.norms reduced.norms.size := by
            calc
              _ = lllPotential' (lllSwapStep' (toPotentialState reduced) ki kp) :=
                hhelper
              _ < lllPotential' (toPotentialState reduced) := habstract
              _ = arrayLLLPotential reduced.norms :=
                toPotentialState_potential reduced
              _ = weightedNormPotential reduced.norms reduced.norms.size :=
                (weightedNormPotential_size_eq_arrayLLLPotential reduced.norms).symm
          simpa only [normsAfterLovaszSwap_size, hcontrol.1, hcontrol.2] using hweighted
    · rw [dif_neg hkMatrix] at hrun
      simp at hrun
  · rw [dif_neg hkPositive] at hrun
    simp at hrun

theorem lllStep_swapped_determinantPotential_lt
    (state output : Generated.StrictRecombine.LLLState)
    (hvalid : ConcreteLLLValid state)
    (houtputValid : ConcreteLLLValid output)
    (hrun : Generated.StrictRecombine.lllStep state =
      .ok (.swapped output)) :
    lllDeterminantPotential output.matrix <
      lllDeterminantPotential state.matrix := by
  have hweighted := lllStep_swapped_weightedNormPotential_lt
    state output hvalid hrun
  have hcastOut := houtputValid.determinantPotential_cast_weighted
  have hcastIn := hvalid.determinantPotential_cast_weighted
  have hcast : (lllDeterminantPotential output.matrix : QQ) <
      (lllDeterminantPotential state.matrix : QQ) := by
    rw [hcastOut, hcastIn]
    simpa only [houtputValid.norms_size, hvalid.norms_size] using hweighted
  exact_mod_cast hcast

theorem lllStep_advanced_source_index_lt
    (state output : Generated.StrictRecombine.LLLState)
    (hrun : Generated.StrictRecombine.lllStep state =
      .ok (.advanced output)) :
    state.k < state.matrix.size := by
  rw [Generated.StrictRecombine.lllStep] at hrun
  by_cases hkPositive : 0 < state.k
  · rw [dif_pos hkPositive] at hrun
    by_cases hkMatrix : state.k < state.matrix.size
    · exact hkMatrix
    · rw [dif_neg hkMatrix] at hrun
      contradiction
  · rw [dif_neg hkPositive] at hrun
    contradiction

theorem lllStep_advanced_concreteRank_lt
    (state output : Generated.StrictRecombine.LLLState)
    (hvalid : ConcreteLLLValid state)
    (houtputValid : ConcreteLLLValid output)
    (hrun : Generated.StrictRecombine.lllStep state =
      .ok (.advanced output)) :
    concreteLLLRank output < concreteLLLRank state := by
  have hcontrol := lllStep_advanced_control state output hrun
  have hpotential := lllStep_advanced_determinantPotential_eq
    state output hvalid houtputValid hrun
  have hsize : output.matrix.size = state.matrix.size := by
    calc
      output.matrix.size = output.norms.size := houtputValid.norms_size.symm
      _ = state.norms.size := congrArg Array.size hcontrol.1
      _ = state.matrix.size := hvalid.norms_size
  unfold concreteLLLRank
  rw [hpotential, hsize, hcontrol.2]
  exact lllLexRank_advanced _ _ _
    (lllStep_advanced_source_index_lt state output hrun)

theorem lllStep_swapped_concreteRank_lt
    (state output : Generated.StrictRecombine.LLLState)
    (hvalid : ConcreteLLLValid state)
    (houtputValid : ConcreteLLLValid output)
    (hrun : Generated.StrictRecombine.lllStep state =
      .ok (.swapped output)) :
    concreteLLLRank output < concreteLLLRank state := by
  have hpotential := lllStep_swapped_determinantPotential_lt
    state output hvalid houtputValid hrun
  have hsize := lllStep_swapped_matrix_size state output hvalid houtputValid hrun
  unfold concreteLLLRank
  rw [hsize]
  exact lllLexRank_swap _ _ _ _ _ hpotential

theorem lllStep_concreteRank_lt_of_valid
    (state : Generated.StrictRecombine.LLLState)
    (branch : Generated.StrictRecombine.LLLStepResult)
    (hvalid : ConcreteLLLValid state)
    (houtputValid : ConcreteLLLValid branch.state)
    (hrun : Generated.StrictRecombine.lllStep state = .ok branch) :
    concreteLLLRank branch.state < concreteLLLRank state := by
  cases branch with
  | advanced output =>
      exact lllStep_advanced_concreteRank_lt state output hvalid houtputValid hrun
  | swapped output =>
      exact lllStep_swapped_concreteRank_lt state output hvalid houtputValid hrun

def subtractRowSuffixList (target source : Array ZZ) (coefficient : ZZ)
    (size index : Nat) : List ZZ :=
  if hindex : index < size then
    (target[index]! - coefficient * source[index]!) ::
      subtractRowSuffixList target source coefficient size (index + 1)
  else []
termination_by size - index
decreasing_by omega

theorem subtractRowsLoop_toList
    (targetM sourceM targetU sourceU resultM resultU outputM outputU : Array ZZ)
    (coefficient : ZZ) (size index : Nat)
    (hrun : Generated.StrictRecombine.subtractRowsLoop
      targetM sourceM targetU sourceU coefficient size index resultM resultU =
        .ok (outputM, outputU)) :
    outputM.toList = resultM.toList ++
        subtractRowSuffixList targetM sourceM coefficient size index ∧
      outputU.toList = resultU.toList ++
        subtractRowSuffixList targetU sourceU coefficient size index := by
  induction hmeasure : size - index generalizing index resultM resultU with
  | zero =>
      rw [Generated.StrictRecombine.subtractRowsLoop] at hrun
      rw [subtractRowSuffixList]
      have hdone : ¬ index < size := by omega
      simp only [dif_neg hdone] at hrun ⊢
      have hout := Except.ok.inj hrun
      injection hout with hM hU
      subst outputM
      subst outputU
      simp [subtractRowSuffixList, hdone]
  | succ remaining ih =>
      rw [Generated.StrictRecombine.subtractRowsLoop] at hrun
      rw [subtractRowSuffixList]
      have hindex : index < size := by omega
      simp only [dif_pos hindex] at hrun ⊢
      split at hrun
      next htM =>
        split at hrun
        next hsM =>
          split at hrun
          next htU =>
            split at hrun
            next hsU =>
              have hnext : size - (index + 1) = remaining := by omega
              have htail := ih
                (resultM.push (targetM[index] - coefficient * sourceM[index]))
                (resultU.push (targetU[index] - coefficient * sourceU[index]))
                (index + 1) hrun hnext
              constructor
              · rw [htail.1]
                simp only [Array.toList_push, List.append_assoc,
                  List.singleton_append]
                rw [getElem!_pos targetM index htM,
                  getElem!_pos sourceM index hsM]
              · rw [htail.2]
                simp only [Array.toList_push, List.append_assoc,
                  List.singleton_append]
                conv_rhs => rw [subtractRowSuffixList, dif_pos hindex]
                rw [getElem!_pos targetU index htU,
                  getElem!_pos sourceU index hsU]
            next hsU => contradiction
          next htU => contradiction
        next hsM => contradiction
      next htM => contradiction

def subtractRowArray (target source : Array ZZ) (coefficient : ZZ)
    (size : Nat) : Array ZZ :=
  (subtractRowSuffixList target source coefficient size 0).toArray

theorem subtractRowSuffixList_length (target source : Array ZZ)
    (coefficient : ZZ) (size index : Nat) :
    (subtractRowSuffixList target source coefficient size index).length =
      size - index := by
  rw [subtractRowSuffixList]
  split
  next hindex =>
    simp only [List.length_cons]
    rw [subtractRowSuffixList_length]
    omega
  next hindex =>
    simp
    omega
termination_by size - index
decreasing_by omega

theorem subtractRowSuffixList_get (target source : Array ZZ)
    (coefficient : ZZ) (size index offset : Nat)
    (hoffset : offset < size - index) :
    (subtractRowSuffixList target source coefficient size index)[offset]! =
      target[index + offset]! - coefficient * source[index + offset]! := by
  rw [subtractRowSuffixList]
  have hindex : index < size := by omega
  simp only [dif_pos hindex]
  cases offset with
  | zero => simp
  | succ offset =>
      simp
      have hrec := subtractRowSuffixList_get target source coefficient size (index + 1)
        offset (by omega)
      simp at hrec
      convert hrec using 1
      rw [show index + (offset + 1) = index + 1 + offset by omega]
termination_by size - index
decreasing_by omega

theorem subtractRowArray_size (target source : Array ZZ) (coefficient : ZZ)
    (size : Nat) :
    (subtractRowArray target source coefficient size).size = size := by
  simp [subtractRowArray, subtractRowSuffixList_length]

theorem subtractRowArray_get (target source : Array ZZ) (coefficient : ZZ)
    (size index : Nat) (hindex : index < size) :
    (subtractRowArray target source coefficient size)[index]! =
      target[index]! - coefficient * source[index]! := by
  unfold subtractRowArray
  rw [getElem!_pos _ index (by
    simpa [subtractRowSuffixList_length] using hindex)]
  have hget := subtractRowSuffixList_get target source coefficient size 0 index
    (by omega)
  rw [getElem!_pos _ index (by
    simpa [subtractRowSuffixList_length] using hindex)] at hget
  simpa using hget

theorem subtractMatrixRows_ok
    (matrix transform outputMatrix outputTransform : Generated.StrictRecombine.LLLMatrix)
    (target source : Nat) (coefficient : ZZ)
    (hrun : Generated.StrictRecombine.subtractMatrixRows matrix transform
      target source coefficient = .ok (outputMatrix, outputTransform)) :
    ∃ (htM : target < matrix.size) (hsM : source < matrix.size)
      (htU : target < transform.size) (hsU : source < transform.size),
      outputMatrix = matrix.set target
          (subtractRowArray matrix[target] matrix[source] coefficient matrix.size) ∧
        outputTransform = transform.set target
          (subtractRowArray transform[target] transform[source] coefficient matrix.size) := by
  unfold Generated.StrictRecombine.subtractMatrixRows at hrun
  split at hrun
  next htM =>
    split at hrun
    next hsM =>
      split at hrun
      next htU =>
        split at hrun
        next hsU =>
          cases hloop : Generated.StrictRecombine.subtractRowsLoop
              matrix[target] matrix[source] transform[target] transform[source]
              coefficient matrix.size 0 #[] #[] with
          | error fault => simp [hloop] at hrun
          | ok rows =>
            rcases rows with ⟨rowM, rowU⟩
            simp only [hloop] at hrun
            have hrows := subtractRowsLoop_toList
              matrix[target] matrix[source] transform[target] transform[source]
              #[] #[] rowM rowU coefficient matrix.size 0 hloop
            have hrowM : rowM = subtractRowArray matrix[target] matrix[source]
                coefficient matrix.size := by
              apply Array.toList_inj.mp
              simpa [subtractRowArray] using hrows.1
            have hrowU : rowU = subtractRowArray transform[target] transform[source]
                coefficient matrix.size := by
              apply Array.toList_inj.mp
              simpa [subtractRowArray] using hrows.2
            have hout := Except.ok.inj hrun
            injection hout with hmatrix htransform
            exact ⟨htM, hsM, htU, hsU,
              hmatrix.symm.trans (congrArg (matrix.set target) hrowM),
              htransform.symm.trans (congrArg (transform.set target) hrowU)⟩
        next hsU => contradiction
      next htU => contradiction
    next hsM => contradiction
  next htM => contradiction

theorem subtractMatrixRows_matrix_size
    (matrix transform outputMatrix outputTransform : Generated.StrictRecombine.LLLMatrix)
    (target source : Nat) (coefficient : ZZ)
    (hrun : Generated.StrictRecombine.subtractMatrixRows matrix transform
      target source coefficient = .ok (outputMatrix, outputTransform)) :
    outputMatrix.size = matrix.size := by
  rcases subtractMatrixRows_ok matrix transform outputMatrix outputTransform
    target source coefficient hrun with ⟨htM, hsM, htU, hsU, rfl, hU⟩
  simp

theorem subtractMatrixRows_matrix_get
    (matrix transform outputMatrix outputTransform : Generated.StrictRecombine.LLLMatrix)
    (target source row : Nat) (coefficient : ZZ)
    (hrun : Generated.StrictRecombine.subtractMatrixRows matrix transform
      target source coefficient = .ok (outputMatrix, outputTransform))
    (hrow : row < matrix.size) :
    outputMatrix[row]! =
      if target = row then
        subtractRowArray matrix[target]! matrix[source]! coefficient matrix.size
      else matrix[row]! := by
  rcases subtractMatrixRows_ok matrix transform outputMatrix outputTransform
    target source coefficient hrun with ⟨htM, hsM, htU, hsU, rfl, hU⟩
  rw [getElem!_pos _ row (by simpa using hrow)]
  rw [Array.getElem_set]
  rw [getElem!_pos matrix target htM, getElem!_pos matrix source hsM,
    getElem!_pos matrix row hrow]

theorem subtractMatrixRows_target_get
    (matrix transform outputMatrix outputTransform : Generated.StrictRecombine.LLLMatrix)
    (target source column : Nat) (coefficient : ZZ)
    (hrun : Generated.StrictRecombine.subtractMatrixRows matrix transform
      target source coefficient = .ok (outputMatrix, outputTransform))
    (hcolumn : column < matrix.size) :
    (outputMatrix[target]!)[column]! =
      (matrix[target]!)[column]! - coefficient * (matrix[source]!)[column]! := by
  rcases subtractMatrixRows_ok matrix transform outputMatrix outputTransform
    target source coefficient hrun with ⟨htM, hsM, htU, hsU, rfl, hU⟩
  rw [getElem!_pos _ target (by simp; exact htM)]
  simp only [Array.getElem_set, if_pos]
  rw [getElem!_pos matrix target htM, getElem!_pos matrix source hsM]
  rw [subtractRowArray_get _ _ _ _ column hcolumn]

theorem basisPrefixMatrix_subtractMatrixRows
    (matrix transform outputMatrix outputTransform : Generated.StrictRecombine.LLLMatrix)
    (target source rowCount : Nat) (coefficient : ZZ)
    (hrun : Generated.StrictRecombine.subtractMatrixRows matrix transform
      target source coefficient = .ok (outputMatrix, outputTransform))
    (htarget : target < rowCount) (hsource : source < rowCount)
    (hrowCount : rowCount ≤ matrix.size) :
    basisPrefixMatrix outputMatrix rowCount matrix.size =
      Matrix.transvection (⟨target, htarget⟩ : Fin rowCount)
        (⟨source, hsource⟩ : Fin rowCount) (-coefficient) *
          basisPrefixMatrix matrix rowCount matrix.size := by
  rw [Matrix.transvection, Matrix.add_mul, Matrix.one_mul]
  funext i column
  by_cases hi : i = (⟨target, htarget⟩ : Fin rowCount)
  · subst i
    unfold basisPrefixMatrix
    rw [subtractMatrixRows_target_get matrix transform outputMatrix outputTransform
      target source column.val coefficient hrun column.isLt]
    simp
    ring
  · unfold basisPrefixMatrix
    have hiMatrix : i.val < matrix.size := lt_of_lt_of_le i.isLt hrowCount
    rw [subtractMatrixRows_matrix_get matrix transform outputMatrix outputTransform
      target source i.val coefficient hrun hiMatrix]
    have htargetI : target ≠ i.val := by
      intro heq
      apply hi
      exact Fin.ext heq.symm
    simp [htargetI, getElem!_pos matrix i.val hiMatrix]
    exact Matrix.single_mul_apply_of_ne (-coefficient)
      (⟨target, htarget⟩ : Fin rowCount) (⟨source, hsource⟩ : Fin rowCount)
      i column hi (basisPrefixMatrix matrix rowCount matrix.size)

theorem gramPrefixMatrix_subtractMatrixRows
    (matrix transform outputMatrix outputTransform : Generated.StrictRecombine.LLLMatrix)
    (target source rowCount : Nat) (coefficient : ZZ)
    (hrun : Generated.StrictRecombine.subtractMatrixRows matrix transform
      target source coefficient = .ok (outputMatrix, outputTransform))
    (htarget : target < rowCount) (hsource : source < rowCount)
    (hrowCount : rowCount ≤ matrix.size) :
    gramPrefixMatrix outputMatrix rowCount =
      Matrix.transvection (⟨target, htarget⟩ : Fin rowCount)
          (⟨source, hsource⟩ : Fin rowCount) (-coefficient) *
        gramPrefixMatrix matrix rowCount *
          (Matrix.transvection (⟨target, htarget⟩ : Fin rowCount)
            (⟨source, hsource⟩ : Fin rowCount) (-coefficient)).transpose := by
  have hsize := subtractMatrixRows_matrix_size matrix transform outputMatrix
    outputTransform target source coefficient hrun
  rw [gramPrefixMatrix_eq_mul_transpose,
    gramPrefixMatrix_eq_mul_transpose, hsize]
  rw [basisPrefixMatrix_subtractMatrixRows matrix transform outputMatrix
    outputTransform target source rowCount coefficient hrun htarget hsource
    hrowCount]
  rw [Matrix.transpose_mul]
  simp only [Matrix.mul_assoc]

theorem map_transvection_intCast
    {dimension : Nat} (target source : Fin dimension) (coefficient : ZZ) :
    (Matrix.transvection target source (-coefficient) :
        Matrix (Fin dimension) (Fin dimension) Int).map
          (fun value : Int => (value : QQ)) =
      Matrix.transvection target source (-(coefficient : QQ)) := by
  ext row column
  by_cases hrowColumn : row = column <;>
    simp [Matrix.transvection, Matrix.single, hrowColumn]

theorem map_mul_intCast_local
    {dimension : Nat}
    (left right : Matrix (Fin dimension) (Fin dimension) Int) :
    (left * right).map (fun value : Int => (value : QQ)) =
      left.map (fun value : Int => (value : QQ)) *
        right.map (fun value : Int => (value : QQ)) := by
  exact Matrix.map_mul (f := Int.castRingHom QQ)

theorem gramPrefixMatrixQQ_subtractMatrixRows
    (matrix transform outputMatrix outputTransform : Generated.StrictRecombine.LLLMatrix)
    (target source rowCount : Nat) (coefficient : ZZ)
    (hrun : Generated.StrictRecombine.subtractMatrixRows matrix transform
      target source coefficient = .ok (outputMatrix, outputTransform))
    (htarget : target < rowCount) (hsource : source < rowCount)
    (hrowCount : rowCount ≤ matrix.size) :
    gramPrefixMatrixQQ outputMatrix rowCount =
      Matrix.transvection (⟨target, htarget⟩ : Fin rowCount)
          (⟨source, hsource⟩ : Fin rowCount) (-(coefficient : QQ)) *
        gramPrefixMatrixQQ matrix rowCount *
          (Matrix.transvection (⟨target, htarget⟩ : Fin rowCount)
            (⟨source, hsource⟩ : Fin rowCount)
            (-(coefficient : QQ))).transpose := by
  have hgram := gramPrefixMatrix_subtractMatrixRows matrix transform
    outputMatrix outputTransform target source rowCount coefficient hrun
    htarget hsource hrowCount
  unfold gramPrefixMatrixQQ
  rw [hgram, map_mul_intCast_local, map_mul_intCast_local]
  rw [map_transvection_intCast]
  congr 2
  ext row column
  have hmap := map_transvection_intCast
    (⟨target, htarget⟩ : Fin rowCount)
    (⟨source, hsource⟩ : Fin rowCount) coefficient
  exact congrArg (fun matrix => matrix column row) hmap

theorem gramPrefixDet_subtractMatrixRows_preserved
    (matrix transform outputMatrix outputTransform : Generated.StrictRecombine.LLLMatrix)
    (target source rowCount : Nat) (coefficient : ZZ)
    (hrun : Generated.StrictRecombine.subtractMatrixRows matrix transform
      target source coefficient = .ok (outputMatrix, outputTransform))
    (htarget : target < rowCount) (hsource : source < rowCount)
    (hrowCount : rowCount ≤ matrix.size) (hne : target ≠ source) :
    Matrix.det (gramPrefixMatrix outputMatrix rowCount) =
      Matrix.det (gramPrefixMatrix matrix rowCount) := by
  have hsize := subtractMatrixRows_matrix_size matrix transform outputMatrix
    outputTransform target source coefficient hrun
  rw [gramPrefixMatrix_eq_mul_transpose,
    gramPrefixMatrix_eq_mul_transpose, hsize]
  rw [basisPrefixMatrix_subtractMatrixRows matrix transform outputMatrix
    outputTransform target source rowCount coefficient hrun htarget hsource hrowCount]
  let E : Matrix (Fin rowCount) (Fin rowCount) Int :=
    Matrix.transvection ⟨target, htarget⟩ ⟨source, hsource⟩ (-coefficient)
  let B := basisPrefixMatrix matrix rowCount matrix.size
  change Matrix.det ((E * B) * (E * B).transpose) =
    Matrix.det (B * B.transpose)
  rw [Matrix.transpose_mul]
  have hassoc : (E * B) * (B.transpose * E.transpose) =
      (E * (B * B.transpose)) * E.transpose := by
    simp only [Matrix.mul_assoc]
  rw [hassoc, Matrix.det_mul, Matrix.det_mul, Matrix.det_transpose]
  have hfinNe : (⟨target, htarget⟩ : Fin rowCount) ≠ ⟨source, hsource⟩ := by
    intro heq
    exact hne (Fin.ext_iff.mp heq)
  rw [show Matrix.det E = 1 by
    exact Matrix.det_transvection_of_ne
      (⟨target, htarget⟩ : Fin rowCount) (⟨source, hsource⟩ : Fin rowCount)
      hfinNe (-coefficient)]
  simp

theorem gramPrefixMatrix_subtractMatrixRows_preserved_before
    (matrix transform outputMatrix outputTransform : Generated.StrictRecombine.LLLMatrix)
    (target source rowCount : Nat) (coefficient : ZZ)
    (hrun : Generated.StrictRecombine.subtractMatrixRows matrix transform
      target source coefficient = .ok (outputMatrix, outputTransform))
    (hrowCount : rowCount ≤ matrix.size) (hbefore : rowCount ≤ target) :
    gramPrefixMatrix outputMatrix rowCount =
      gramPrefixMatrix matrix rowCount := by
  funext i j
  unfold gramPrefixMatrix
  have hsize := subtractMatrixRows_matrix_size matrix transform outputMatrix
    outputTransform target source coefficient hrun
  rw [hsize]
  apply Finset.sum_congr rfl
  intro column _
  have hiMatrix : i.val < matrix.size := lt_of_lt_of_le i.isLt hrowCount
  have hjMatrix : j.val < matrix.size := lt_of_lt_of_le j.isLt hrowCount
  rw [subtractMatrixRows_matrix_get matrix transform outputMatrix outputTransform
    target source i.val coefficient hrun hiMatrix]
  rw [subtractMatrixRows_matrix_get matrix transform outputMatrix outputTransform
    target source j.val coefficient hrun hjMatrix]
  have htargetI : target ≠ i.val := by omega
  have htargetJ : target ≠ j.val := by omega
  simp [htargetI, htargetJ, getElem!_pos matrix i.val hiMatrix,
    getElem!_pos matrix j.val hjMatrix]

theorem gramPrefixDet_subtractMatrixRows_preserved_before
    (matrix transform outputMatrix outputTransform : Generated.StrictRecombine.LLLMatrix)
    (target source rowCount : Nat) (coefficient : ZZ)
    (hrun : Generated.StrictRecombine.subtractMatrixRows matrix transform
      target source coefficient = .ok (outputMatrix, outputTransform))
    (hrowCount : rowCount ≤ matrix.size) (hbefore : rowCount ≤ target) :
    Matrix.det (gramPrefixMatrix outputMatrix rowCount) =
      Matrix.det (gramPrefixMatrix matrix rowCount) := by
  rw [gramPrefixMatrix_subtractMatrixRows_preserved_before matrix transform
    outputMatrix outputTransform target source rowCount coefficient hrun
    hrowCount hbefore]

theorem concreteGramSchmidt_subtractMatrixRows
    (state : Generated.StrictRecombine.LLLState)
    (outputMatrix outputTransform : Generated.StrictRecombine.LLLMatrix)
    (source : Nat) (coefficient : ZZ)
    (hvalid : ConcreteLLLExecutionValid state)
    (hsourceK : source < state.k)
    (hrun : Generated.StrictRecombine.subtractMatrixRows
      state.matrix state.transform state.k source coefficient =
        .ok (outputMatrix, outputTransform)) :
    ConcreteGramSchmidt
      { state with
        matrix := outputMatrix
        transform := outputTransform
        mu := sizeReduceMuResult state.mu state.k source coefficient } := by
  intro rowCount hrowCountOutput
  have hmatrixSize := subtractMatrixRows_matrix_size state.matrix
    state.transform outputMatrix outputTransform state.k source coefficient hrun
  have hrowCount : rowCount ≤ state.matrix.size := by
    simpa [hmatrixSize] using hrowCountOutput
  have hkMatrix : state.k < state.matrix.size := by
    rcases subtractMatrixRows_ok state.matrix state.transform outputMatrix
      outputTransform state.k source coefficient hrun with
      ⟨hk, hsource, htargetTransform, hsourceTransform, hmatrix, htransform⟩
    exact hk
  change gramPrefixMatrixQQ outputMatrix rowCount =
    gsLowerPrefix { state with
        mu := sizeReduceMuResult state.mu state.k source coefficient }
        rowCount *
      gsNormDiagonal state rowCount *
        (gsLowerPrefix { state with
          mu := sizeReduceMuResult state.mu state.k source coefficient }
          rowCount).transpose
  by_cases hbefore : rowCount ≤ state.k
  · have hgram := gramPrefixMatrix_subtractMatrixRows_preserved_before
      state.matrix state.transform outputMatrix outputTransform state.k source
      rowCount coefficient hrun hrowCount hbefore
    have hlower := gsLowerPrefix_sizeReduceMuResult_before state source
      rowCount coefficient hvalid hbefore hkMatrix.le
    have hfactor := hvalid.gram_schmidt rowCount hrowCount
    unfold gramPrefixMatrixQQ
    rw [hgram]
    change gramPrefixMatrixQQ state.matrix rowCount = _
    rw [hlower]
    exact hfactor
  · have hkCount : state.k < rowCount := by omega
    have hsourceCount : source < rowCount := lt_trans hsourceK hkCount
    have hgram := gramPrefixMatrixQQ_subtractMatrixRows state.matrix
      state.transform outputMatrix outputTransform state.k source rowCount
      coefficient hrun hkCount hsourceCount hrowCount
    have hlower := gsLowerPrefix_sizeReduceMuResult state source rowCount
      coefficient hvalid hsourceK hkCount hrowCount
    have hfactor := hvalid.gram_schmidt rowCount hrowCount
    rw [hgram, hlower, hfactor]
    rw [Matrix.transpose_mul]
    simp only [Matrix.mul_assoc]

theorem sizeReduceAt_preserves_valid
    (state output : Generated.StrictRecombine.LLLState) (j : Nat)
    (hvalid : ConcreteLLLValid state) (hj : j < state.k)
    (hrun : Generated.StrictRecombine.sizeReduceAt state j = .ok output) :
    ConcreteLLLValid output := by
  unfold Generated.StrictRecombine.sizeReduceAt at hrun
  split at hrun
  next hkMu =>
    split at hrun
    next hjMu =>
      cases hround : Generated.StrictRecombine.roundQQ state.mu[state.k][j] with
      | error fault => simp [hround] at hrun
      | ok q =>
        simp only [hround] at hrun
        split at hrun
        next hzero =>
          have hout := Except.ok.inj hrun
          subst output
          exact hvalid
        next hnonzero =>
          cases hsubtract : Generated.StrictRecombine.subtractMatrixRows
              state.matrix state.transform state.k j q with
          | error fault => simp [hsubtract] at hrun
          | ok matrices =>
            rcases matrices with ⟨matrix', transform'⟩
            simp only [hsubtract] at hrun
            cases hmu : Generated.StrictRecombine.reduceMuPrefixLoop
                state.mu state.k j q j 0 with
            | error fault => simp [hmu] at hrun
            | ok mu' =>
              simp only [hmu] at hrun
              split at hrun
              next hk' =>
                split at hrun
                next hj' =>
                  have hout := Except.ok.inj hrun
                  subst output
                  have hmatrixSize := subtractMatrixRows_matrix_size
                    state.matrix state.transform matrix' transform'
                    state.k j q hsubtract
                  rcases subtractMatrixRows_ok state.matrix state.transform
                    matrix' transform' state.k j q hsubtract with
                    ⟨htargetMatrix, hsourceMatrix, htargetTransform,
                      hsourceTransform, hmatrixExact, htransformExact⟩
                  refine {
                    norms_size := ?_
                    rows_square := ?_
                    norms_positive := hvalid.norms_positive
                    gram_prefix := ?_ }
                  · change state.norms.size = matrix'.size
                    rw [hmatrixSize]
                    exact hvalid.norms_size
                  · intro row hrow
                    change row < matrix'.size at hrow
                    change matrix'[row].size = matrix'.size
                    have hrowState : row < state.matrix.size := by
                      simpa [hmatrixSize] using hrow
                    by_cases htarget : state.k = row
                    · subst row
                      rw [← getElem!_pos matrix' state.k hrow]
                      rw [subtractMatrixRows_matrix_get state.matrix state.transform
                        matrix' transform' state.k j state.k q hsubtract
                        htargetMatrix]
                      simp only [if_pos]
                      rw [subtractRowArray_size]
                      exact hmatrixSize.symm
                    · rw [← getElem!_pos matrix' row hrow]
                      rw [subtractMatrixRows_matrix_get state.matrix state.transform
                        matrix' transform' state.k j row q hsubtract hrowState]
                      simp [htarget, getElem!_pos state.matrix row hrowState,
                        hvalid.rows_square row hrowState, hmatrixSize]
                  · intro rowCount hrowCount
                    change rowCount ≤ matrix'.size at hrowCount
                    change ((Matrix.det (gramPrefixMatrix matrix' rowCount) : Int) : QQ) =
                      prefixNormProduct state.norms rowCount
                    have hrowCountState : rowCount ≤ state.matrix.size := by
                      simpa [hmatrixSize] using hrowCount
                    have hdet : Matrix.det (gramPrefixMatrix matrix' rowCount) =
                        Matrix.det (gramPrefixMatrix state.matrix rowCount) := by
                      by_cases hbefore : rowCount ≤ state.k
                      · exact gramPrefixDet_subtractMatrixRows_preserved_before
                          state.matrix state.transform matrix' transform'
                          state.k j rowCount q hsubtract hrowCountState hbefore
                      · exact gramPrefixDet_subtractMatrixRows_preserved
                          state.matrix state.transform matrix' transform'
                          state.k j rowCount q hsubtract (by omega) (by omega)
                          hrowCountState (by omega)
                    rw [hdet]
                    exact hvalid.gram_prefix rowCount hrowCountState
                next hj' => contradiction
              next hk' => contradiction
    next hjMu => contradiction
  next hkMu => contradiction

theorem sizeReduceAt_preserves_execution_valid
    (state output : Generated.StrictRecombine.LLLState) (source : Nat)
    (hvalid : ConcreteLLLExecutionValid state) (hsourceK : source < state.k)
    (hrun : Generated.StrictRecombine.sizeReduceAt state source = .ok output) :
    ConcreteLLLExecutionValid output := by
  have hbasic := sizeReduceAt_preserves_valid state output source
    hvalid.toConcreteLLLValid hsourceK hrun
  rcases sizeReduceAt_mu_eq state output source hvalid hsourceK hrun with
    ⟨coefficient, hround, hmuOutput⟩
  refine {
    norms_size := hbasic.norms_size
    mu_size := ?_
    rows_square := hbasic.rows_square
    mu_rows_square := ?_
    norms_positive := hbasic.norms_positive
    gram_schmidt := ?_ }
  · rw [hmuOutput, sizeReduceMuResult_size]
    exact hvalid.mu_size.trans (by
      rw [← hbasic.norms_size, ← hvalid.norms_size]
      exact congrArg Array.size (sizeReduceAt_preserves_norms_k
        state output source hrun).1.symm)
  · intro row hrow
    have hrowState : row < state.mu.size := by
      rw [hmuOutput, sizeReduceMuResult_size] at hrow
      exact hrow
    rw [← getElem!_pos output.mu row hrow, hmuOutput]
    rw [sizeReduceMuResult_row_size state.mu state.k source row coefficient
      hrowState]
    rw [getElem!_pos state.mu row hrowState]
    rw [hvalid.mu_rows_square row hrowState]
    have hmatrixSize : output.matrix.size = state.matrix.size := by
      calc
        output.matrix.size = output.norms.size := hbasic.norms_size.symm
        _ = state.norms.size := congrArg Array.size
          (sizeReduceAt_preserves_norms_k state output source hrun).1
        _ = state.matrix.size := hvalid.norms_size
    exact hmatrixSize.symm
  · have hrun' := hrun
    unfold Generated.StrictRecombine.sizeReduceAt at hrun'
    split at hrun'
    next hk =>
      split at hrun'
      next hsource =>
        cases hroundCase : Generated.StrictRecombine.roundQQ
            state.mu[state.k][source] with
        | error fault => simp [hroundCase] at hrun'
        | ok q =>
          simp only [hroundCase] at hrun'
          split at hrun'
          next hzero =>
            have hout := Except.ok.inj hrun'
            subst output
            exact hvalid.gram_schmidt
          next hnonzero =>
            cases hsubtract : Generated.StrictRecombine.subtractMatrixRows
                state.matrix state.transform state.k source q with
            | error fault => simp [hsubtract] at hrun'
            | ok matrices =>
              rcases matrices with ⟨matrix', transform'⟩
              simp only [hsubtract] at hrun'
              cases hmu : Generated.StrictRecombine.reduceMuPrefixLoop
                  state.mu state.k source q source 0 with
              | error fault => simp [hmu] at hrun'
              | ok mu' =>
                simp only [hmu] at hrun'
                split at hrun'
                next hk' =>
                  split at hrun'
                  next hsource' =>
                    have hout := Except.ok.inj hrun'
                    subst output
                    have hsourceMu : source < state.mu.size := by omega
                    have hlimitK : source ≤ state.mu[state.k].size := by
                      rw [hvalid.mu_rows_square state.k hk]
                      rw [← hvalid.mu_size]
                      omega
                    have hlimitSource : source ≤ state.mu[source].size := by
                      rw [hvalid.mu_rows_square source hsourceMu]
                      rw [← hvalid.mu_size]
                      omega
                    have hmuExact := reduceMuPrefixLoop_output_eq state.mu mu'
                      state.k source q source 0 hk hsourceMu hlimitK
                      hlimitSource (by omega) hmu
                    subst mu'
                    change ConcreteGramSchmidt
                      { state with
                        matrix := matrix'
                        transform := transform'
                        mu := ((reduceMuPrefixArray state.mu state.k source q
                          source 0).set state.k
                          ((reduceMuPrefixArray state.mu state.k source q
                            source 0)[state.k].set source
                            ((reduceMuPrefixArray state.mu state.k source q
                              source 0)[state.k][source] - (q : QQ)))) }
                    have hmuFinal :
                        ((reduceMuPrefixArray state.mu state.k source q source 0).set
                          state.k
                          ((reduceMuPrefixArray state.mu state.k source q source 0)[state.k].set
                            source
                            ((reduceMuPrefixArray state.mu state.k source q source 0)[state.k][source] -
                              (q : QQ)))) =
                          sizeReduceMuResult state.mu state.k source q := by
                      simp [sizeReduceMuResult, hnonzero, Array.setIfInBounds,
                        hk', hsource']
                    rw [hmuFinal]
                    exact concreteGramSchmidt_subtractMatrixRows state matrix'
                      transform' source q hvalid hsourceK hsubtract
                  next hsource' => contradiction
                next hk' => contradiction
      next hsource => contradiction
    next hk => contradiction

theorem extraSizeReduceLoop_preserves_valid
    (remaining : Nat) (state output : Generated.StrictRecombine.LLLState)
    (hvalid : ConcreteLLLValid state) (hremaining : remaining < state.k)
    (hrun : Generated.StrictRecombine.extraSizeReduceLoop remaining state =
      .ok output) :
    ConcreteLLLValid output := by
  induction remaining generalizing state output with
  | zero =>
      simp [Generated.StrictRecombine.extraSizeReduceLoop] at hrun
      subst output
      exact hvalid
  | succ remaining ih =>
      rw [Generated.StrictRecombine.extraSizeReduceLoop] at hrun
      cases hstep : Generated.StrictRecombine.sizeReduceAt state remaining with
      | error fault => simp [hstep] at hrun
      | ok next =>
        simp only [hstep] at hrun
        have hnextValid := sizeReduceAt_preserves_valid state next remaining
          hvalid (by omega) hstep
        have hcontrol := sizeReduceAt_preserves_norms_k state next remaining hstep
        exact ih next output hnextValid (by rw [hcontrol.2]; omega)
          hrun

theorem extraSizeReduceLoop_preserves_execution_valid
    (remaining : Nat) (state output : Generated.StrictRecombine.LLLState)
    (hvalid : ConcreteLLLExecutionValid state) (hremaining : remaining < state.k)
    (hrun : Generated.StrictRecombine.extraSizeReduceLoop remaining state =
      .ok output) :
    ConcreteLLLExecutionValid output := by
  induction remaining generalizing state output with
  | zero =>
      simp [Generated.StrictRecombine.extraSizeReduceLoop] at hrun
      subst output
      exact hvalid
  | succ remaining ih =>
      rw [Generated.StrictRecombine.extraSizeReduceLoop] at hrun
      cases hstep : Generated.StrictRecombine.sizeReduceAt state remaining with
      | error fault => simp [hstep] at hrun
      | ok next =>
        simp only [hstep] at hrun
        have hnextValid := sizeReduceAt_preserves_execution_valid state next
          remaining hvalid (by omega) hstep
        have hcontrol := sizeReduceAt_preserves_norms_k state next remaining hstep
        exact ih next output hnextValid (by rw [hcontrol.2]; omega) hrun

theorem lllStep_advanced_preserves_valid
    (state output : Generated.StrictRecombine.LLLState)
    (hvalid : ConcreteLLLValid state)
    (hrun : Generated.StrictRecombine.lllStep state =
      .ok (.advanced output)) :
    ConcreteLLLValid output := by
  rw [Generated.StrictRecombine.lllStep] at hrun
  by_cases hkPositive : 0 < state.k
  · rw [dif_pos hkPositive] at hrun
    by_cases hkMatrix : state.k < state.matrix.size
    · rw [dif_pos hkMatrix] at hrun
      cases hreduce : Generated.StrictRecombine.sizeReduceAt state
          (state.k - 1) with
      | error fault => simp [hreduce] at hrun
      | ok reduced =>
        simp only [hreduce] at hrun
        have hreducedValid := sizeReduceAt_preserves_valid state reduced
          (state.k - 1) hvalid (by omega) hreduce
        have hcontrol := sizeReduceAt_preserves_norms_k state reduced
          (state.k - 1) hreduce
        split at hrun
        next hkNorm =>
          split at hrun
          next hpredNorm =>
            split at hrun
            next hkMu =>
              split at hrun
              next hpredMu =>
                dsimp at hrun
                split at hrun
                next hlovasz =>
                  cases hextra : Generated.StrictRecombine.extraSizeReduceLoop
                      (reduced.k - 1) reduced with
                  | error fault => simp [hextra] at hrun
                  | ok fullyReduced =>
                    simp only [hextra] at hrun
                    have hfullyValid := extraSizeReduceLoop_preserves_valid
                      (reduced.k - 1) reduced fullyReduced hreducedValid
                      (by omega) hextra
                    have hout := Except.ok.inj hrun
                    injection hout with hstate
                    subst output
                    exact hfullyValid.withK (fullyReduced.k + 1)
                next hlovasz =>
                  repeat' first | split at hrun | simp_all
              next hpredMu => contradiction
            next hkMu => contradiction
          next hpredNorm => contradiction
        next hkNorm => contradiction
    · rw [dif_neg hkMatrix] at hrun
      contradiction
  · rw [dif_neg hkPositive] at hrun
    contradiction

theorem lllStep_advanced_preserves_execution_valid
    (state output : Generated.StrictRecombine.LLLState)
    (hvalid : ConcreteLLLExecutionValid state)
    (hrun : Generated.StrictRecombine.lllStep state =
      .ok (.advanced output)) :
    ConcreteLLLExecutionValid output := by
  rw [Generated.StrictRecombine.lllStep] at hrun
  by_cases hkPositive : 0 < state.k
  · rw [dif_pos hkPositive] at hrun
    by_cases hkMatrix : state.k < state.matrix.size
    · rw [dif_pos hkMatrix] at hrun
      cases hreduce : Generated.StrictRecombine.sizeReduceAt state
          (state.k - 1) with
      | error fault => simp [hreduce] at hrun
      | ok reduced =>
        simp only [hreduce] at hrun
        have hreducedValid := sizeReduceAt_preserves_execution_valid state
          reduced (state.k - 1) hvalid (by omega) hreduce
        have hcontrol := sizeReduceAt_preserves_norms_k state reduced
          (state.k - 1) hreduce
        split at hrun
        next hkNorm =>
          split at hrun
          next hpredNorm =>
            split at hrun
            next hkMu =>
              split at hrun
              next hpredMu =>
                dsimp at hrun
                split at hrun
                next hlovasz =>
                  cases hextra : Generated.StrictRecombine.extraSizeReduceLoop
                      (reduced.k - 1) reduced with
                  | error fault => simp [hextra] at hrun
                  | ok fullyReduced =>
                    simp only [hextra] at hrun
                    have hfullyValid :=
                      extraSizeReduceLoop_preserves_execution_valid
                        (reduced.k - 1) reduced fullyReduced hreducedValid
                        (by omega) hextra
                    have hout := Except.ok.inj hrun
                    injection hout with hstate
                    subst output
                    exact hfullyValid.withK (fullyReduced.k + 1)
                next hlovasz =>
                  repeat' first | split at hrun | simp_all
              next hpredMu => contradiction
            next hkMu => contradiction
          next hpredNorm => contradiction
        next hkNorm => contradiction
    · rw [dif_neg hkMatrix] at hrun
      contradiction
  · rw [dif_neg hkPositive] at hrun
    contradiction

theorem swapMatrixRows_ok
    (matrix output : Generated.StrictRecombine.LLLMatrix) (left right : Nat)
    (hrun : Generated.StrictRecombine.swapMatrixRows matrix left right =
      .ok output) :
    ∃ (hleft : left < matrix.size) (hright : right < matrix.size),
      output = (matrix.set left matrix[right]).set right matrix[left] (by
        simpa) := by
  unfold Generated.StrictRecombine.swapMatrixRows at hrun
  split at hrun
  next hleft =>
    split at hrun
    next hright =>
      exact ⟨hleft, hright, (Except.ok.inj hrun).symm⟩
    next hright => contradiction
  next hleft => contradiction

theorem swapMatrixRows_size
    (matrix output : Generated.StrictRecombine.LLLMatrix) (left right : Nat)
    (hrun : Generated.StrictRecombine.swapMatrixRows matrix left right =
      .ok output) :
    output.size = matrix.size := by
  rcases swapMatrixRows_ok matrix output left right hrun with
    ⟨hleft, hright, rfl⟩
  simp

theorem swapMatrixRows_get
    (matrix output : Generated.StrictRecombine.LLLMatrix) (left right row : Nat)
    (hrun : Generated.StrictRecombine.swapMatrixRows matrix left right =
      .ok output)
    (hrow : row < matrix.size) :
    output[row]'(by simpa [swapMatrixRows_size matrix output left right hrun]) =
      if right = row then matrix[left]!
      else if left = row then matrix[right]!
      else matrix[row] := by
  rcases swapMatrixRows_ok matrix output left right hrun with
    ⟨hleft, hright, rfl⟩
  simp only [getElem!_pos matrix left hleft,
    getElem!_pos matrix right hright, Array.getElem_set]

theorem swapMatrixRows_get_fin
    (matrix output : Generated.StrictRecombine.LLLMatrix) (left right rowCount : Nat)
    (hrun : Generated.StrictRecombine.swapMatrixRows matrix left right =
      .ok output)
    (hleftCount : left < rowCount) (hrightCount : right < rowCount)
    (hrowCount : rowCount ≤ matrix.size) (row : Fin rowCount) :
    output[row.val]! = matrix[((Equiv.swap
      (⟨left, hleftCount⟩ : Fin rowCount)
      (⟨right, hrightCount⟩ : Fin rowCount)) row).val]! := by
  have hrow : row.val < matrix.size := lt_of_lt_of_le row.isLt hrowCount
  have houtSize := swapMatrixRows_size matrix output left right hrun
  rw [getElem!_pos output row.val (by simpa [houtSize])]
  rw [swapMatrixRows_get matrix output left right row hrun hrow]
  by_cases hrightRow : row = ⟨right, hrightCount⟩
  · subst row
    simp only [Equiv.swap_apply_right]
    rw [getElem!_pos matrix left (lt_of_lt_of_le hleftCount hrowCount)]
    simp
  · by_cases hleftRow : row = ⟨left, hleftCount⟩
    · subst row
      simp only [Equiv.swap_apply_left]
      rw [getElem!_pos matrix right (lt_of_lt_of_le hrightCount hrowCount)]
      have hrightLeft : right ≠ left := by
        intro heq
        apply hrightRow
        exact Fin.ext heq.symm
      simp [hrightLeft]
    · have hrightVal : right ≠ row.val := by
        intro heq
        apply hrightRow
        exact Fin.ext heq.symm
      have hleftVal : left ≠ row.val := by
        intro heq
        apply hleftRow
        exact Fin.ext heq.symm
      rw [Equiv.swap_apply_of_ne_of_ne hleftRow hrightRow]
      rw [getElem!_pos matrix row.val hrow]
      simp [hrightVal, hleftVal]

theorem gramPrefixMatrix_swap_reindex
    (matrix output : Generated.StrictRecombine.LLLMatrix) (left right rowCount : Nat)
    (hrun : Generated.StrictRecombine.swapMatrixRows matrix left right =
      .ok output)
    (hleftCount : left < rowCount) (hrightCount : right < rowCount)
    (hrowCount : rowCount ≤ matrix.size) :
    gramPrefixMatrix output rowCount = Matrix.reindex
      (Equiv.swap (⟨left, hleftCount⟩ : Fin rowCount)
        (⟨right, hrightCount⟩ : Fin rowCount))
      (Equiv.swap (⟨left, hleftCount⟩ : Fin rowCount)
        (⟨right, hrightCount⟩ : Fin rowCount))
      (gramPrefixMatrix matrix rowCount) := by
  funext i j
  unfold gramPrefixMatrix Matrix.reindex
  simp only [Equiv.coe_fn_mk]
  have hsize := swapMatrixRows_size matrix output left right hrun
  rw [hsize]
  apply Finset.sum_congr rfl
  intro c _
  rw [swapMatrixRows_get_fin matrix output left right rowCount hrun
    hleftCount hrightCount hrowCount i]
  rw [swapMatrixRows_get_fin matrix output left right rowCount hrun
    hleftCount hrightCount hrowCount j]
  simp

theorem concreteGramSchmidt_lovaszSwap_of_above
    (state : Generated.StrictRecombine.LLLState)
    (swappedMatrix : Generated.StrictRecombine.LLLMatrix)
    (k rowCount : Nat) (muNew : QQ)
    (hvalid : ConcreteLLLExecutionValid state)
    (hkPositive : 0 < k) (hkCount : k < rowCount)
    (hrowCount : rowCount ≤ state.matrix.size)
    (hswap : Generated.StrictRecombine.swapMatrixRows state.matrix k (k - 1) =
      .ok swappedMatrix)
    (hmuNew : muNew =
      ((state.mu[k]!)[k - 1]!) * state.norms[k - 1]! /
        (state.norms[k]! + ((state.mu[k]!)[k - 1]!) *
          ((state.mu[k]!)[k - 1]!) * state.norms[k - 1]!)) :
    gramPrefixMatrixQQ swappedMatrix rowCount =
      gsLowerPrefix
          { state with
            matrix := swappedMatrix
            mu := (lovaszSwapMuResult state.mu k
              ((state.mu[k]!)[k - 1]!) muNew)
            norms := (Generated.StrictRecombine.normsAfterLovaszSwap
              state.norms k ((state.mu[k]!)[k - 1]!)) }
          rowCount *
        gsNormDiagonal
          { state with
            matrix := swappedMatrix
            mu := (lovaszSwapMuResult state.mu k
              ((state.mu[k]!)[k - 1]!) muNew)
            norms := (Generated.StrictRecombine.normsAfterLovaszSwap
              state.norms k ((state.mu[k]!)[k - 1]!)) }
          rowCount *
        (gsLowerPrefix
          { state with
            matrix := swappedMatrix
            mu := (lovaszSwapMuResult state.mu k
              ((state.mu[k]!)[k - 1]!) muNew)
            norms := (Generated.StrictRecombine.normsAfterLovaszSwap
              state.norms k ((state.mu[k]!)[k - 1]!)) }
          rowCount).transpose := by
  let previous : Fin rowCount := ⟨k - 1, by omega⟩
  let current : Fin rowCount := ⟨k, hkCount⟩
  let permutation := Equiv.swap previous current
  have hgram := gramPrefixMatrix_swap_reindex state.matrix swappedMatrix
    k (k - 1) rowCount hswap hkCount (by omega) hrowCount
  have hold := hvalid.gram_schmidt rowCount hrowCount
  have hlower := gsLowerPrefix_lovaszSwapMuResult_mul state k rowCount
    muNew hvalid hkPositive hkCount hrowCount
  have hdiagonal := gsNormDiagonal_normsAfterLovaszSwap state k rowCount
    muNew hvalid hkPositive hkCount hrowCount hmuNew
  unfold gramPrefixMatrixQQ
  rw [hgram]
  rw [Equiv.swap_comm (⟨k, hkCount⟩ : Fin rowCount)
    (⟨k - 1, by omega⟩ : Fin rowCount)]
  change Matrix.reindex permutation permutation
      (gramPrefixMatrixQQ state.matrix rowCount) = _
  rw [hold]
  rw [reindex_ldl_eq_rowReindex permutation]
  change
    rowReindexMatrix permutation (gsLowerPrefix state rowCount) *
        gsNormDiagonal state rowCount *
      (rowReindexMatrix permutation
        (gsLowerPrefix state rowCount)).transpose = _
  unfold rowReindexMatrix
  dsimp [permutation, previous, current]
  rw [hlower]
  change
    (gsLowerPrefix
        { state with mu := (lovaszSwapMuResult state.mu k
          ((state.mu[k]!)[k - 1]!) muNew) } rowCount *
      lovaszLocalTransform (⟨k - 1, by omega⟩ : Fin rowCount)
        (⟨k, hkCount⟩ : Fin rowCount) ((state.mu[k]!)[k - 1]!) muNew) *
      gsNormDiagonal state rowCount *
      (gsLowerPrefix
        { state with mu := (lovaszSwapMuResult state.mu k
          ((state.mu[k]!)[k - 1]!) muNew) } rowCount *
      lovaszLocalTransform (⟨k - 1, by omega⟩ : Fin rowCount)
        (⟨k, hkCount⟩ : Fin rowCount) ((state.mu[k]!)[k - 1]!) muNew).transpose = _
  rw [Matrix.transpose_mul]
  calc
    _ = gsLowerPrefix
          { state with mu := (lovaszSwapMuResult state.mu k
            ((state.mu[k]!)[k - 1]!) muNew) } rowCount *
        (lovaszLocalTransform (⟨k - 1, by omega⟩ : Fin rowCount)
            (⟨k, hkCount⟩ : Fin rowCount) ((state.mu[k]!)[k - 1]!) muNew *
          gsNormDiagonal state rowCount *
          (lovaszLocalTransform (⟨k - 1, by omega⟩ : Fin rowCount)
            (⟨k, hkCount⟩ : Fin rowCount)
            ((state.mu[k]!)[k - 1]!) muNew).transpose) *
        (gsLowerPrefix
          { state with mu := (lovaszSwapMuResult state.mu k
            ((state.mu[k]!)[k - 1]!) muNew) } rowCount).transpose := by
          simp only [Matrix.mul_assoc]
    _ = gsLowerPrefix
          { state with mu := (lovaszSwapMuResult state.mu k
            ((state.mu[k]!)[k - 1]!) muNew) } rowCount *
        gsNormDiagonal
          { state with norms :=
            (Generated.StrictRecombine.normsAfterLovaszSwap state.norms k
              ((state.mu[k]!)[k - 1]!)) } rowCount *
        (gsLowerPrefix
          { state with mu := (lovaszSwapMuResult state.mu k
            ((state.mu[k]!)[k - 1]!) muNew) } rowCount).transpose := by
          rw [hdiagonal]
    _ = _ := by
      rfl

theorem gramPrefixDet_swap_preserved_of_both
    (matrix output : Generated.StrictRecombine.LLLMatrix) (left right rowCount : Nat)
    (hrun : Generated.StrictRecombine.swapMatrixRows matrix left right =
      .ok output)
    (hleftCount : left < rowCount) (hrightCount : right < rowCount)
    (hrowCount : rowCount ≤ matrix.size) :
    Matrix.det (gramPrefixMatrix output rowCount) =
      Matrix.det (gramPrefixMatrix matrix rowCount) := by
  rw [gramPrefixMatrix_swap_reindex matrix output left right rowCount hrun
    hleftCount hrightCount hrowCount]
  exact Matrix.det_reindex_self _ _

theorem swapMatrixRows_get_before
    (matrix output : Generated.StrictRecombine.LLLMatrix) (left right rowCount : Nat)
    (hrun : Generated.StrictRecombine.swapMatrixRows matrix left right =
      .ok output)
    (hrowCount : rowCount ≤ matrix.size)
    (hbeforeLeft : rowCount ≤ left) (hbeforeRight : rowCount ≤ right)
    (row : Fin rowCount) :
    output[row.val]! = matrix[row.val]! := by
  have hrow : row.val < matrix.size := lt_of_lt_of_le row.isLt hrowCount
  have houtSize := swapMatrixRows_size matrix output left right hrun
  rw [getElem!_pos output row.val (by simpa [houtSize])]
  rw [swapMatrixRows_get matrix output left right row hrun hrow]
  have hleft : left ≠ row.val := by omega
  have hright : right ≠ row.val := by omega
  simp [hleft, hright, getElem!_pos matrix row.val hrow]

theorem gramPrefixMatrix_swap_eq_of_before
    (matrix output : Generated.StrictRecombine.LLLMatrix)
    (left right rowCount : Nat)
    (hrun : Generated.StrictRecombine.swapMatrixRows matrix left right =
      .ok output)
    (hrowCount : rowCount ≤ matrix.size)
    (hbeforeLeft : rowCount ≤ left) (hbeforeRight : rowCount ≤ right) :
    gramPrefixMatrix output rowCount = gramPrefixMatrix matrix rowCount := by
  funext i j
  unfold gramPrefixMatrix
  have hsize := swapMatrixRows_size matrix output left right hrun
  rw [hsize]
  apply Finset.sum_congr rfl
  intro c _
  rw [swapMatrixRows_get_before matrix output left right rowCount hrun
    hrowCount hbeforeLeft hbeforeRight i]
  rw [swapMatrixRows_get_before matrix output left right rowCount hrun
    hrowCount hbeforeLeft hbeforeRight j]

theorem gsLowerPrefix_lovaszSwapMuResult_of_before
    (state : Generated.StrictRecombine.LLLState) (k rowCount : Nat)
    (muNew : QQ) (hvalid : ConcreteLLLExecutionValid state)
    (hkPositive : 0 < k) (hkMatrix : k < state.matrix.size)
    (hrowCount : rowCount ≤ k - 1) :
    gsLowerPrefix
        { state with mu := (lovaszSwapMuResult state.mu k
          ((state.mu[k]!)[k - 1]!) muNew) } rowCount =
      gsLowerPrefix state rowCount := by
  have hkMu : k < state.mu.size := by simpa [hvalid.mu_size] using hkMatrix
  have hrowsSquare : ∀ index (hindex : index < state.mu.size),
      state.mu[index]!.size = state.mu.size := by
    intro index hindex
    rw [getElem!_pos state.mu index hindex]
    rw [hvalid.mu_rows_square index (by simpa [hvalid.mu_size] using hindex)]
    exact hvalid.mu_size.symm
  funext row column
  unfold gsLowerPrefix
  by_cases hcolumnRow : column.val < row.val
  · simp only [hcolumnRow, if_pos]
    have hrowMu : row.val < state.mu.size := by omega
    have hcolumnMu : column.val < state.mu.size := by omega
    rw [lovaszSwapMuResult_entry state.mu k row.val column.val
      ((state.mu[k]!)[k - 1]!) muNew hkPositive hkMu hrowMu hcolumnMu
      hrowsSquare]
    have hrowBefore : row.val < k - 1 := lt_of_lt_of_le row.isLt hrowCount
    simp only [show row.val ≠ k by omega, show row.val ≠ k - 1 by omega,
      show ¬k < row.val by omega, if_false]
  · simp only [hcolumnRow, if_false]

theorem gsNormDiagonal_normsAfterLovaszSwap_of_before
    (state : Generated.StrictRecombine.LLLState) (k rowCount : Nat)
    (hvalid : ConcreteLLLExecutionValid state)
    (hkPositive : 0 < k) (hkMatrix : k < state.matrix.size)
    (hrowCount : rowCount ≤ k - 1) :
    gsNormDiagonal
        { state with norms := (Generated.StrictRecombine.normsAfterLovaszSwap
          state.norms k ((state.mu[k]!)[k - 1]!)) } rowCount =
      gsNormDiagonal state rowCount := by
  have hkNorm : k < state.norms.size := by simpa [hvalid.norms_size] using hkMatrix
  have hpredNorm : k - 1 < state.norms.size := by omega
  have hdeltaPos : 0 < state.norms[k]! +
      ((state.mu[k]!)[k - 1]!) * ((state.mu[k]!)[k - 1]!) *
        state.norms[k - 1]! := by
    have hkPos := hvalid.norms_positive k hkNorm
    have hpPos := hvalid.norms_positive (k - 1) hpredNorm
    have hkPos' : 0 < state.norms[k]! := by
      simpa only [getElem!_pos state.norms k hkNorm] using hkPos
    have hpPos' : 0 < state.norms[k - 1]! := by
      simpa only [getElem!_pos state.norms (k - 1) hpredNorm] using hpPos
    nlinarith [sq_nonneg ((state.mu[k]!)[k - 1]!)]
  have hdeltaNe : state.norms[k] +
      ((state.mu[k]!)[k - 1]!) * ((state.mu[k]!)[k - 1]!) *
        state.norms[k - 1] ≠ 0 := by
    simpa only [getElem!_pos state.norms k hkNorm,
      getElem!_pos state.norms (k - 1) hpredNorm] using ne_of_gt hdeltaPos
  unfold gsNormDiagonal
  congr 1
  funext index
  have hindexNorm : index.val < state.norms.size := by omega
  rw [normsAfterLovaszSwap_get state.norms k ((state.mu[k]!)[k - 1]!)
    hkNorm hpredNorm hdeltaNe index.val hindexNorm]
  have hindexBefore : index.val < k - 1 := lt_of_lt_of_le index.isLt hrowCount
  rw [if_neg (show k - 1 ≠ index.val by omega)]
  rw [if_neg (show k ≠ index.val by omega)]
  exact (getElem!_pos state.norms index.val hindexNorm).symm

theorem concreteGramSchmidt_lovaszSwap_of_before
    (state : Generated.StrictRecombine.LLLState)
    (swappedMatrix : Generated.StrictRecombine.LLLMatrix)
    (k rowCount : Nat) (muNew : QQ)
    (hvalid : ConcreteLLLExecutionValid state)
    (hkPositive : 0 < k) (hkMatrix : k < state.matrix.size)
    (hrowCount : rowCount ≤ k - 1)
    (hswap : Generated.StrictRecombine.swapMatrixRows state.matrix k (k - 1) =
      .ok swappedMatrix) :
    gramPrefixMatrixQQ swappedMatrix rowCount =
      gsLowerPrefix
          { state with
            matrix := swappedMatrix
            mu := lovaszSwapMuResult state.mu k
              ((state.mu[k]!)[k - 1]!) muNew
            norms := Generated.StrictRecombine.normsAfterLovaszSwap
              state.norms k ((state.mu[k]!)[k - 1]!) }
          rowCount *
        gsNormDiagonal
          { state with
            matrix := swappedMatrix
            mu := lovaszSwapMuResult state.mu k
              ((state.mu[k]!)[k - 1]!) muNew
            norms := Generated.StrictRecombine.normsAfterLovaszSwap
              state.norms k ((state.mu[k]!)[k - 1]!) }
          rowCount *
        (gsLowerPrefix
          { state with
            matrix := swappedMatrix
            mu := lovaszSwapMuResult state.mu k
              ((state.mu[k]!)[k - 1]!) muNew
            norms := Generated.StrictRecombine.normsAfterLovaszSwap
              state.norms k ((state.mu[k]!)[k - 1]!) }
          rowCount).transpose := by
  have hrowMatrix : rowCount ≤ state.matrix.size := by omega
  have hgram := gramPrefixMatrix_swap_eq_of_before state.matrix swappedMatrix
    k (k - 1) rowCount hswap hrowMatrix (by omega) hrowCount
  have hlower := gsLowerPrefix_lovaszSwapMuResult_of_before state k rowCount
    muNew hvalid hkPositive hkMatrix hrowCount
  have hdiagonal := gsNormDiagonal_normsAfterLovaszSwap_of_before state k
    rowCount hvalid hkPositive hkMatrix hrowCount
  have hlower' : gsLowerPrefix
        { state with
          matrix := swappedMatrix
          mu := lovaszSwapMuResult state.mu k
            ((state.mu[k]!)[k - 1]!) muNew
          norms := Generated.StrictRecombine.normsAfterLovaszSwap
            state.norms k ((state.mu[k]!)[k - 1]!) }
        rowCount = gsLowerPrefix state rowCount := by
    simpa only [gsLowerPrefix] using hlower
  have hdiagonal' : gsNormDiagonal
        { state with
          matrix := swappedMatrix
          mu := lovaszSwapMuResult state.mu k
            ((state.mu[k]!)[k - 1]!) muNew
          norms := Generated.StrictRecombine.normsAfterLovaszSwap
            state.norms k ((state.mu[k]!)[k - 1]!) }
        rowCount = gsNormDiagonal state rowCount := by
    simpa only [gsNormDiagonal] using hdiagonal
  change (gramPrefixMatrix swappedMatrix rowCount).map
      (fun value : Int => (value : QQ)) = _
  rw [hgram]
  change gramPrefixMatrixQQ state.matrix rowCount = _
  rw [hvalid.gram_schmidt rowCount hrowMatrix]
  rw [hlower', hdiagonal']

theorem concreteGramSchmidt_lovaszSwap_at_boundary
    (state : Generated.StrictRecombine.LLLState)
    (swappedMatrix : Generated.StrictRecombine.LLLMatrix)
    (k : Nat) (muNew : QQ)
    (hvalid : ConcreteLLLExecutionValid state)
    (hkPositive : 0 < k) (hkMatrix : k < state.matrix.size)
    (hswap : Generated.StrictRecombine.swapMatrixRows state.matrix k (k - 1) =
      .ok swappedMatrix)
    (hmuNew : muNew =
      ((state.mu[k]!)[k - 1]!) * state.norms[k - 1]! /
        (state.norms[k]! + ((state.mu[k]!)[k - 1]!) *
          ((state.mu[k]!)[k - 1]!) * state.norms[k - 1]!)) :
    gramPrefixMatrixQQ swappedMatrix k =
      gsLowerPrefix
          { state with
            matrix := swappedMatrix
            mu := lovaszSwapMuResult state.mu k
              ((state.mu[k]!)[k - 1]!) muNew
            norms := Generated.StrictRecombine.normsAfterLovaszSwap
              state.norms k ((state.mu[k]!)[k - 1]!) }
          k *
        gsNormDiagonal
          { state with
            matrix := swappedMatrix
            mu := lovaszSwapMuResult state.mu k
              ((state.mu[k]!)[k - 1]!) muNew
            norms := Generated.StrictRecombine.normsAfterLovaszSwap
              state.norms k ((state.mu[k]!)[k - 1]!) }
          k *
        (gsLowerPrefix
          { state with
            matrix := swappedMatrix
            mu := lovaszSwapMuResult state.mu k
              ((state.mu[k]!)[k - 1]!) muNew
            norms := Generated.StrictRecombine.normsAfterLovaszSwap
              state.norms k ((state.mu[k]!)[k - 1]!) }
          k).transpose := by
  have habove := concreteGramSchmidt_lovaszSwap_of_above state swappedMatrix
    k (k + 1) muNew hvalid hkPositive (by omega) (by omega) hswap hmuNew
  funext row column
  have hentry := congr_fun (congr_fun habove (Fin.castSucc row))
    (Fin.castSucc column)
  have hgramEntry : gramPrefixMatrixQQ swappedMatrix (k + 1)
      (Fin.castSucc row) (Fin.castSucc column) =
      gramPrefixMatrixQQ swappedMatrix k row column := by
    rfl
  unfold gsNormDiagonal at hentry ⊢
  rw [Matrix.mul_apply] at hentry ⊢
  simp_rw [Matrix.mul_diagonal, Matrix.transpose_apply] at hentry ⊢
  rw [Fin.sum_univ_castSucc] at hentry
  rw [hgramEntry] at hentry
  simp only [gsLowerPrefix, Fin.val_last, Fin.val_castSucc,
    show ¬k < row.val by omega, show ¬k < column.val by omega,
    row.castSucc_ne_last, column.castSucc_ne_last,
    if_false, mul_zero, add_zero, Fin.castSucc_inj] at hentry ⊢
  exact hentry

theorem concreteGramSchmidt_lovaszSwap
    (state : Generated.StrictRecombine.LLLState)
    (swappedMatrix : Generated.StrictRecombine.LLLMatrix)
    (k : Nat) (muNew : QQ)
    (hvalid : ConcreteLLLExecutionValid state)
    (hkPositive : 0 < k) (hkMatrix : k < state.matrix.size)
    (hswap : Generated.StrictRecombine.swapMatrixRows state.matrix k (k - 1) =
      .ok swappedMatrix)
    (hmuNew : muNew =
      ((state.mu[k]!)[k - 1]!) * state.norms[k - 1]! /
        (state.norms[k]! + ((state.mu[k]!)[k - 1]!) *
          ((state.mu[k]!)[k - 1]!) * state.norms[k - 1]!)) :
    ConcreteGramSchmidt
      { state with
        matrix := swappedMatrix
        mu := lovaszSwapMuResult state.mu k
          ((state.mu[k]!)[k - 1]!) muNew
        norms := Generated.StrictRecombine.normsAfterLovaszSwap
          state.norms k ((state.mu[k]!)[k - 1]!) } := by
  intro rowCount hrowCount
  have hsize := swapMatrixRows_size state.matrix swappedMatrix k (k - 1) hswap
  have hrowCountOld : rowCount ≤ state.matrix.size := by simpa [hsize] using hrowCount
  by_cases hbefore : rowCount ≤ k - 1
  · exact concreteGramSchmidt_lovaszSwap_of_before state swappedMatrix k
      rowCount muNew hvalid hkPositive hkMatrix hbefore hswap
  · by_cases habove : k < rowCount
    · exact concreteGramSchmidt_lovaszSwap_of_above state swappedMatrix k
        rowCount muNew hvalid hkPositive habove hrowCountOld hswap hmuNew
    · have hboundary : rowCount = k := by omega
      subst rowCount
      exact concreteGramSchmidt_lovaszSwap_at_boundary state swappedMatrix k
        muNew hvalid hkPositive hkMatrix hswap hmuNew

theorem lovaszSwapMuResult_row_size
    (mu : Generated.StrictRecombine.QQMatrix) (k row : Nat)
    (muOld muNew : QQ) (hkPositive : 0 < k) (hk : k < mu.size)
    (hrow : row < mu.size)
    (hrowsSquare : ∀ index (hindex : index < mu.size),
      mu[index]!.size = mu.size) :
    (lovaszSwapMuResult mu k muOld muNew)[row]!.size = mu.size := by
  rw [lovaszSwapMuResult_get mu k row muOld muNew hkPositive hk hrow]
  by_cases hrowK : row = k
  · subst row
    simp only [if_pos]
    simp [hrowsSquare (k - 1) (by omega)]
  · rw [if_neg hrowK]
    by_cases hrowPred : row = k - 1
    · rw [if_pos hrowPred]
      exact hrowsSquare k hk
    · rw [if_neg hrowPred]
      by_cases hafter : k < row
      · rw [if_pos hafter, updateMuAfterSwapRow_size]
        exact hrowsSquare row hrow
      · rw [if_neg hafter]
        exact hrowsSquare row hrow

theorem lovaszSwapCorrectedMu_row_size
    (mu : Generated.StrictRecombine.QQMatrix) (k row : Nat)
    (muNew : QQ) (hkPositive : 0 < k) (hk : k < mu.size)
    (hrow : row < mu.size)
    (hrowsSquare : ∀ index (hindex : index < mu.size),
      mu[index]!.size = mu.size) :
    (lovaszSwapCorrectedMu mu k muNew)[row]!.size = mu.size := by
  rw [lovaszSwapCorrectedMu_get mu k row muNew hkPositive hk hrow]
  by_cases hrowK : row = k
  · subst row
    simp only [if_pos]
    simp [hrowsSquare (k - 1) (by omega)]
  · rw [if_neg hrowK]
    by_cases hrowPred : row = k - 1
    · rw [if_pos hrowPred]
      exact hrowsSquare k hk
    · rw [if_neg hrowPred]
      exact hrowsSquare row hrow

theorem concreteLLLExecutionValid_lovaszSwap
    (state : Generated.StrictRecombine.LLLState)
    (swappedMatrix : Generated.StrictRecombine.LLLMatrix)
    (k : Nat) (muNew : QQ)
    (hvalid : ConcreteLLLExecutionValid state)
    (hkPositive : 0 < k) (hkMatrix : k < state.matrix.size)
    (hswap : Generated.StrictRecombine.swapMatrixRows state.matrix k (k - 1) =
      .ok swappedMatrix)
    (hmuNew : muNew =
      ((state.mu[k]!)[k - 1]!) * state.norms[k - 1]! /
        (state.norms[k]! + ((state.mu[k]!)[k - 1]!) *
          ((state.mu[k]!)[k - 1]!) * state.norms[k - 1]!)) :
    ConcreteLLLExecutionValid
      { state with
        matrix := swappedMatrix
        mu := lovaszSwapMuResult state.mu k
          ((state.mu[k]!)[k - 1]!) muNew
        norms := Generated.StrictRecombine.normsAfterLovaszSwap
          state.norms k ((state.mu[k]!)[k - 1]!) } := by
  have hsize := swapMatrixRows_size state.matrix swappedMatrix k (k - 1) hswap
  have hkNorm : k < state.norms.size := by simpa [hvalid.norms_size] using hkMatrix
  have hpredNorm : k - 1 < state.norms.size := by omega
  have hkMu : k < state.mu.size := by simpa [hvalid.mu_size] using hkMatrix
  have hmuRows : ∀ index (hindex : index < state.mu.size),
      state.mu[index]!.size = state.mu.size := by
    intro index hindex
    rw [getElem!_pos state.mu index hindex]
    rw [hvalid.mu_rows_square index (by simpa [hvalid.mu_size] using hindex)]
    exact hvalid.mu_size.symm
  constructor
  · rw [normsAfterLovaszSwap_size, hvalid.norms_size, hsize]
  · rw [lovaszSwapMuResult_size, hvalid.mu_size, hsize]
  · intro row hrow
    change row < swappedMatrix.size at hrow
    have hrowOld : row < state.matrix.size := by simpa [hsize] using hrow
    have hget := swapMatrixRows_get state.matrix swappedMatrix k (k - 1)
      row hswap hrowOld
    rw [hget, hsize]
    split
    next hright =>
      simpa only [getElem!_pos state.matrix k hkMatrix] using
        hvalid.rows_square k hkMatrix
    next hright =>
      split
      next hleft =>
        have hpredMatrix : k - 1 < state.matrix.size := by omega
        simpa only [getElem!_pos state.matrix (k - 1) hpredMatrix] using
          hvalid.rows_square (k - 1) hpredMatrix
      next hleft => exact hvalid.rows_square row hrowOld
  · intro row hrow
    change row < (lovaszSwapMuResult state.mu k
      ((state.mu[k]!)[k - 1]!) muNew).size at hrow
    have hrowOld : row < state.mu.size := by
      simpa [lovaszSwapMuResult_size] using hrow
    rw [show (lovaszSwapMuResult state.mu k
        ((state.mu[k]!)[k - 1]!) muNew)[row] =
      (lovaszSwapMuResult state.mu k
        ((state.mu[k]!)[k - 1]!) muNew)[row]! by
          rw [getElem!_pos _ row hrow]]
    rw [lovaszSwapMuResult_row_size state.mu k row
      ((state.mu[k]!)[k - 1]!) muNew hkPositive hkMu hrowOld hmuRows]
    exact hvalid.mu_size.trans hsize.symm
  · exact normsAfterLovaszSwap_positive state.norms k
      ((state.mu[k]!)[k - 1]!) hkNorm hpredNorm hkPositive
      hvalid.norms_positive
  · exact concreteGramSchmidt_lovaszSwap state swappedMatrix k muNew hvalid
      hkPositive hkMatrix hswap hmuNew

theorem lllStep_swapped_preserves_execution_valid
    (state output : Generated.StrictRecombine.LLLState)
    (hvalid : ConcreteLLLExecutionValid state)
    (hrun : Generated.StrictRecombine.lllStep state = .ok (.swapped output)) :
    ConcreteLLLExecutionValid output := by
  rw [Generated.StrictRecombine.lllStep] at hrun
  by_cases hkPositive : 0 < state.k
  · rw [dif_pos hkPositive] at hrun
    by_cases hkMatrix : state.k < state.matrix.size
    · rw [dif_pos hkMatrix] at hrun
      cases hreduce : Generated.StrictRecombine.sizeReduceAt state
          (state.k - 1) with
      | error fault => simp [hreduce] at hrun
      | ok reduced =>
        simp only [hreduce] at hrun
        have hreducedValid := sizeReduceAt_preserves_execution_valid state
          reduced (state.k - 1) hvalid (by omega) hreduce
        split at hrun
        next hkNorm =>
          split at hrun
          next hpredNorm =>
            split at hrun
            next hkMu =>
              split at hrun
              next hpredMu =>
                dsimp at hrun
                split at hrun
                next hlovasz =>
                  repeat' first | split at hrun | simp_all
                next hlovasz =>
                  have hnewPos : 0 < reduced.norms[reduced.k] +
                      reduced.mu[reduced.k][reduced.k - 1] *
                        reduced.mu[reduced.k][reduced.k - 1] *
                          reduced.norms[reduced.k - 1] := by
                    have hkPos := hreducedValid.norms_positive reduced.k hkNorm
                    have hpPos := hreducedValid.norms_positive
                      (reduced.k - 1) hpredNorm
                    nlinarith [sq_nonneg reduced.mu[reduced.k][reduced.k - 1]]
                  have hnewNe : reduced.norms[reduced.k] +
                      reduced.mu[reduced.k][reduced.k - 1] *
                        reduced.mu[reduced.k][reduced.k - 1] *
                          reduced.norms[reduced.k - 1] ≠ 0 := ne_of_gt hnewPos
                  simp only [if_pos hnewNe] at hrun
                  cases hswapMatrix : Generated.StrictRecombine.swapMatrixRows
                      reduced.matrix reduced.k (reduced.k - 1) with
                  | error fault => simp [hswapMatrix] at hrun
                  | ok matrix' =>
                    simp only [hswapMatrix] at hrun
                    cases hswapTransform : Generated.StrictRecombine.swapMatrixRows
                        reduced.transform reduced.k (reduced.k - 1) with
                    | error fault => simp [hswapTransform] at hrun
                    | ok transform' =>
                      simp only [hswapTransform] at hrun
                      cases hswapMu : Generated.StrictRecombine.swapQQRows
                          reduced.mu reduced.k (reduced.k - 1) with
                      | error fault => simp [hswapMu] at hrun
                      | ok swappedMu =>
                        simp only [hswapMu] at hrun
                        split at hrun
                        next hkSwapped =>
                          split at hrun
                          next hpredSwapped =>
                            let muOld := reduced.mu[reduced.k][reduced.k - 1]
                            let newNorm := reduced.norms[reduced.k] +
                              muOld * muOld * reduced.norms[reduced.k - 1]
                            let muNew := muOld * reduced.norms[reduced.k - 1] /
                              newNorm
                            let correctedMu := swappedMu.set reduced.k
                              (swappedMu[reduced.k].set (reduced.k - 1) muNew)
                            cases hupdate :
                                Generated.StrictRecombine.updateMuAfterSwapLoop
                                  correctedMu reduced.k muOld muNew
                                    (reduced.k + 1) with
                            | error fault => simp [correctedMu, muNew, newNorm,
                                muOld, hupdate] at hrun
                            | ok finalMu =>
                              simp only [correctedMu, muNew, newNorm, muOld,
                                hupdate] at hrun
                              have hkPositiveReduced : 0 < reduced.k := by
                                have hcontrol := sizeReduceAt_preserves_norms_k state
                                  reduced (state.k - 1) hreduce
                                rw [hcontrol.2]
                                exact hkPositive
                              have hmuRows : ∀ index
                                  (hindex : index < reduced.mu.size),
                                  reduced.mu[index]!.size = reduced.mu.size := by
                                intro index hindex
                                rw [getElem!_pos reduced.mu index hindex]
                                rw [hreducedValid.mu_rows_square index (by
                                  simpa [hreducedValid.mu_size] using hindex)]
                                exact hreducedValid.mu_size.symm
                              have hcorrected :
                                  swappedMu.set reduced.k
                                      (swappedMu[reduced.k].set
                                        (reduced.k - 1)
                                        (reduced.mu[reduced.k][reduced.k - 1] *
                                          reduced.norms[reduced.k - 1] /
                                          (reduced.norms[reduced.k] +
                                            reduced.mu[reduced.k][reduced.k - 1] *
                                            reduced.mu[reduced.k][reduced.k - 1] *
                                            reduced.norms[reduced.k - 1]))) =
                                    lovaszSwapCorrectedMu reduced.mu reduced.k
                                      (reduced.mu[reduced.k][reduced.k - 1] *
                                        reduced.norms[reduced.k - 1] /
                                        (reduced.norms[reduced.k] +
                                          reduced.mu[reduced.k][reduced.k - 1] *
                                          reduced.mu[reduced.k][reduced.k - 1] *
                                          reduced.norms[reduced.k - 1])) := by
                                unfold lovaszSwapCorrectedMu
                                rw [← swapQQRows_output_eq reduced.mu swappedMu
                                  reduced.k (reduced.k - 1) hswapMu]
                                simp [Array.setIfInBounds, hkSwapped, hpredSwapped]
                              have hcorrectedSize :
                                  (swappedMu.set reduced.k
                                    (swappedMu[reduced.k].set (reduced.k - 1)
                                      (reduced.mu[reduced.k][reduced.k - 1] *
                                        reduced.norms[reduced.k - 1] /
                                        (reduced.norms[reduced.k] +
                                          reduced.mu[reduced.k][reduced.k - 1] *
                                          reduced.mu[reduced.k][reduced.k - 1] *
                                          reduced.norms[reduced.k - 1])))).size =
                                      reduced.mu.size := by
                                rw [hcorrected, lovaszSwapCorrectedMu_size]
                              have hrowsK : ∀ row (hrow : row <
                                  (swappedMu.set reduced.k
                                    (swappedMu[reduced.k].set (reduced.k - 1)
                                      (reduced.mu[reduced.k][reduced.k - 1] *
                                        reduced.norms[reduced.k - 1] /
                                        (reduced.norms[reduced.k] +
                                          reduced.mu[reduced.k][reduced.k - 1] *
                                          reduced.mu[reduced.k][reduced.k - 1] *
                                          reduced.norms[reduced.k - 1])))).size),
                                  reduced.k <
                                    (swappedMu.set reduced.k
                                      (swappedMu[reduced.k].set (reduced.k - 1)
                                        (reduced.mu[reduced.k][reduced.k - 1] *
                                          reduced.norms[reduced.k - 1] /
                                          (reduced.norms[reduced.k] +
                                            reduced.mu[reduced.k][reduced.k - 1] *
                                            reduced.mu[reduced.k][reduced.k - 1] *
                                            reduced.norms[reduced.k - 1]))))[row].size := by
                                intro row hrow
                                have hrowOld : row < reduced.mu.size := by
                                  rw [← hcorrectedSize]
                                  exact hrow
                                rw [← getElem!_pos _ row hrow]
                                have hget := congrArg (fun matrix => matrix[row]!)
                                  hcorrected
                                have hsizeRow :
                                    (swappedMu.set reduced.k
                                      (swappedMu[reduced.k].set (reduced.k - 1)
                                        (reduced.mu[reduced.k][reduced.k - 1] *
                                          reduced.norms[reduced.k - 1] /
                                          (reduced.norms[reduced.k] +
                                            reduced.mu[reduced.k][reduced.k - 1] *
                                            reduced.mu[reduced.k][reduced.k - 1] *
                                            reduced.norms[reduced.k - 1]))))[row]!.size =
                                      reduced.mu.size := by
                                  calc
                                    _ = (lovaszSwapCorrectedMu reduced.mu
                                        reduced.k
                                        (reduced.mu[reduced.k][reduced.k - 1] *
                                          reduced.norms[reduced.k - 1] /
                                          (reduced.norms[reduced.k] +
                                            reduced.mu[reduced.k][reduced.k - 1] *
                                            reduced.mu[reduced.k][reduced.k - 1] *
                                            reduced.norms[reduced.k - 1])))[row]!.size :=
                                      congrArg Array.size hget
                                    _ = reduced.mu.size :=
                                      lovaszSwapCorrectedMu_row_size reduced.mu
                                        reduced.k row _ hkPositiveReduced hkMu
                                        hrowOld hmuRows
                                omega
                              have hrowsPred : ∀ row (hrow : row <
                                  (swappedMu.set reduced.k
                                    (swappedMu[reduced.k].set (reduced.k - 1)
                                      (reduced.mu[reduced.k][reduced.k - 1] *
                                        reduced.norms[reduced.k - 1] /
                                        (reduced.norms[reduced.k] +
                                          reduced.mu[reduced.k][reduced.k - 1] *
                                          reduced.mu[reduced.k][reduced.k - 1] *
                                          reduced.norms[reduced.k - 1])))).size),
                                  reduced.k - 1 <
                                    (swappedMu.set reduced.k
                                      (swappedMu[reduced.k].set (reduced.k - 1)
                                        (reduced.mu[reduced.k][reduced.k - 1] *
                                          reduced.norms[reduced.k - 1] /
                                          (reduced.norms[reduced.k] +
                                            reduced.mu[reduced.k][reduced.k - 1] *
                                            reduced.mu[reduced.k][reduced.k - 1] *
                                            reduced.norms[reduced.k - 1]))))[row].size := by
                                intro row hrow
                                exact lt_trans (by omega) (hrowsK row hrow)
                              have hfinalMu := lovaszSwapMuResult_of_generated
                                reduced.mu swappedMu finalMu reduced.k
                                reduced.mu[reduced.k][reduced.k - 1]
                                (reduced.mu[reduced.k][reduced.k - 1] *
                                  reduced.norms[reduced.k - 1] /
                                  (reduced.norms[reduced.k] +
                                    reduced.mu[reduced.k][reduced.k - 1] *
                                    reduced.mu[reduced.k][reduced.k - 1] *
                                    reduced.norms[reduced.k - 1]))
                                hswapMu hkSwapped hpredSwapped hrowsK hrowsPred
                                hupdate
                              have hout := Except.ok.inj hrun
                              injection hout with hstate
                              subst output
                              rw [hfinalMu]
                              have hswapValid :=
                                concreteLLLExecutionValid_lovaszSwap reduced
                                matrix' reduced.k
                                (reduced.mu[reduced.k][reduced.k - 1] *
                                  reduced.norms[reduced.k - 1] /
                                  (reduced.norms[reduced.k] +
                                    reduced.mu[reduced.k][reduced.k - 1] *
                                    reduced.mu[reduced.k][reduced.k - 1] *
                                    reduced.norms[reduced.k - 1]))
                                hreducedValid hkPositiveReduced
                                (by simpa [hreducedValid.norms_size] using hkNorm)
                                hswapMatrix (by
                                  simp only [getElem!_pos reduced.mu reduced.k hkMu,
                                    getElem!_pos reduced.mu[reduced.k]
                                      (reduced.k - 1) hpredMu,
                                    getElem!_pos reduced.norms reduced.k hkNorm,
                                    getElem!_pos reduced.norms (reduced.k - 1)
                                      hpredNorm])
                              have hswapValid' : ConcreteLLLExecutionValid
                                  { reduced with
                                    matrix := matrix'
                                    mu := lovaszSwapMuResult reduced.mu reduced.k
                                      reduced.mu[reduced.k][reduced.k - 1]
                                      (reduced.mu[reduced.k][reduced.k - 1] *
                                        reduced.norms[reduced.k - 1] /
                                        (reduced.norms[reduced.k] +
                                          reduced.mu[reduced.k][reduced.k - 1] *
                                          reduced.mu[reduced.k][reduced.k - 1] *
                                          reduced.norms[reduced.k - 1]))
                                    norms :=
                                      Generated.StrictRecombine.normsAfterLovaszSwap
                                        reduced.norms reduced.k
                                          reduced.mu[reduced.k][reduced.k - 1] } := by
                                simpa only [getElem!_pos reduced.mu reduced.k hkMu,
                                  getElem!_pos reduced.mu[reduced.k]
                                    (reduced.k - 1) hpredMu,
                                  getElem!_pos reduced.norms reduced.k hkNorm,
                                  getElem!_pos reduced.norms (reduced.k - 1)
                                    hpredNorm] using hswapValid
                              simpa only [getElem!_pos reduced.mu reduced.k hkMu,
                                getElem!_pos reduced.mu[reduced.k]
                                  (reduced.k - 1) hpredMu] using
                                (hswapValid'.withTransform transform').withK
                                  (Nat.max (reduced.k - 1) 1)
                          next hpredSwapped => contradiction
                        next hkSwapped => contradiction
              next hpredMu => contradiction
            next hkMu => contradiction
          next hpredNorm => contradiction
        next hkNorm => contradiction
    · rw [dif_neg hkMatrix] at hrun
      contradiction
  · rw [dif_neg hkPositive] at hrun
    contradiction

theorem lllStep_preserves_execution_valid
    (state : Generated.StrictRecombine.LLLState)
    (branch : Generated.StrictRecombine.LLLStepResult)
    (hvalid : ConcreteLLLExecutionValid state)
    (hrun : Generated.StrictRecombine.lllStep state = .ok branch) :
    ConcreteLLLExecutionValid branch.state := by
  cases branch with
  | advanced output =>
      exact lllStep_advanced_preserves_execution_valid state output hvalid hrun
  | swapped output =>
      exact lllStep_swapped_preserves_execution_valid state output hvalid hrun

/-- The genuine well-founded certificate for the generated C++ LLL loop.
Its validity predicate is the executable `G = L D Lᵀ` invariant, and its
rank is the concrete determinant/index lexicographic rank. -/
noncomputable def concreteLLLTermination :
    Generated.StrictRecombine.LLLTermination where
  valid := ConcreteLLLExecutionValid
  rank := concreteLLLRank
  step_valid := by
    intro current branch hvalid hrun
    exact lllStep_preserves_execution_valid current branch hvalid hrun
  step_decreases := by
    intro current branch hvalid hrun
    have houtput := lllStep_preserves_execution_valid current branch hvalid hrun
    exact lllStep_concreteRank_lt_of_valid current branch
      hvalid.toConcreteLLLValid houtput.toConcreteLLLValid hrun

theorem concreteLLLMainLoop_preserves_execution_valid
    (state output : Generated.StrictRecombine.LLLState)
    (hvalid : ConcreteLLLExecutionValid state)
    (hrun : Generated.StrictRecombine.lllMainLoop concreteLLLTermination
      state hvalid = .ok output) :
    ConcreteLLLExecutionValid output := by
  induction hmeasure : concreteLLLRank state using Nat.strong_induction_on
      generalizing state output with
  | h measure ih =>
      rw [Generated.StrictRecombine.lllMainLoop] at hrun
      split at hrun
      next hk =>
        split at hrun
        next hstep => contradiction
        next branch hstep =>
          have hnextValid := lllStep_preserves_execution_valid state branch
            hvalid hstep
          have hdecrease := lllStep_concreteRank_lt_of_valid state branch
            hvalid.toConcreteLLLValid hnextValid.toConcreteLLLValid hstep
          rw [hmeasure] at hdecrease
          exact ih (concreteLLLRank branch.state) hdecrease branch.state output
            hnextValid hrun rfl
      next hk =>
        have hout := Except.ok.inj hrun
        subst output
        exact hvalid

theorem lllStep_preserves_matrix_size
    (state : Generated.StrictRecombine.LLLState)
    (branch : Generated.StrictRecombine.LLLStepResult)
    (hvalid : ConcreteLLLExecutionValid state)
    (hrun : Generated.StrictRecombine.lllStep state = .ok branch) :
    branch.state.matrix.size = state.matrix.size := by
  have houtput := lllStep_preserves_execution_valid state branch hvalid hrun
  cases branch with
  | advanced output =>
      have hcontrol := lllStep_advanced_control state output hrun
      calc
        output.matrix.size = output.norms.size := houtput.norms_size.symm
        _ = state.norms.size := congrArg Array.size hcontrol.1
        _ = state.matrix.size := hvalid.norms_size
  | swapped output =>
      exact lllStep_swapped_matrix_size state output hvalid.toConcreteLLLValid
        houtput.toConcreteLLLValid hrun

theorem concreteLLLMainLoop_preserves_matrix_size
    (state output : Generated.StrictRecombine.LLLState)
    (hvalid : ConcreteLLLExecutionValid state)
    (hrun : Generated.StrictRecombine.lllMainLoop concreteLLLTermination
      state hvalid = .ok output) :
    output.matrix.size = state.matrix.size := by
  induction hmeasure : concreteLLLRank state using Nat.strong_induction_on
      generalizing state output with
  | h measure ih =>
      rw [Generated.StrictRecombine.lllMainLoop] at hrun
      split at hrun
      next hk =>
        split at hrun
        next hstep => contradiction
        next branch hstep =>
          have hnextValid := lllStep_preserves_execution_valid state branch
            hvalid hstep
          have hdecrease := lllStep_concreteRank_lt_of_valid state branch
            hvalid.toConcreteLLLValid hnextValid.toConcreteLLLValid hstep
          have htail := ih (concreteLLLRank branch.state) (by
              simpa [hmeasure] using hdecrease) branch.state output hnextValid
            hrun rfl
          exact htail.trans
            (lllStep_preserves_matrix_size state branch hvalid hstep)
      next hk =>
        have hout := Except.ok.inj hrun
        subst output
        rfl

theorem initializeLLL_concrete_valid
    (matrix : Generated.StrictRecombine.LLLMatrix)
    (mu : Generated.StrictRecombine.QQMatrix) (norms : Array QQ)
    (transform : Generated.StrictRecombine.LLLMatrix)
    (hinput : ConcreteLLLInputValid matrix)
    (hrun : Generated.StrictRecombine.initializeLLL matrix =
      .ok (mu, norms, transform)) :
    ConcreteLLLExecutionValid
      (Generated.StrictRecombine.LLLState.mk matrix transform mu norms 1) := by
  rw [Generated.StrictRecombine.initializeLLL] at hrun
  by_cases hsize : 0 < matrix.size
  · rw [dif_pos hsize] at hrun
    cases htransform : Generated.StrictRecombine.makeInitialMatrix matrix.size 1 with
    | error fault => simp [htransform] at hrun
    | ok initialTransform =>
      simp only [htransform] at hrun
      let initialMu := Generated.StrictRecombine.zeroQQMatrix
        matrix.size matrix.size
      let initialNorms := Array.replicate matrix.size (0 : QQ)
      have hrowZero := hinput.rows_square 0 hsize
      rw [dotRows_eq_fin_sum matrix[0] matrix[0] (by omega)] at hrun
      simp only at hrun
      let firstNorm : QQ :=
        ((∑ k : Fin matrix.size,
          (matrix[0]!)[k.val]! * (matrix[0]!)[k.val]! : ZZ) : QQ)
      have hfirstNorm :
          ((∑ k : Fin matrix[0].size,
            matrix[0][k.val] * matrix[0][k.val] : ZZ) : QQ) = firstNorm := by
        unfold firstNorm
        congr 1
        apply Fintype.sum_equiv (finCongr hrowZero)
        intro index
        simp [getElem!_pos matrix 0 hsize,
          getElem!_pos matrix[0] index.val index.isLt]
      rw [hfirstNorm] at hrun
      let seededNorms := initialNorms.set 0 firstNorm (by
        simp [initialNorms, hsize])
      cases hloop : Generated.StrictRecombine.initializeGramSchmidtLoop
          matrix 1 initialMu seededNorms with
      | error fault =>
        rw [hloop] at hrun
        contradiction
      | ok loopOutput =>
        rw [hloop] at hrun
        have hresult := Except.ok.inj hrun
        simp only at hresult
        rcases hresult with ⟨rfl, rfl, rfl⟩
        have hinitial : ProcessedGramSchmidtValid
            (Generated.StrictRecombine.LLLState.mk matrix #[] initialMu
              seededNorms 1) 1 := by
          simpa [initialMu, initialNorms, seededNorms, firstNorm,
            getElem!_pos matrix 0 hsize] using
            initialProcessedGramSchmidtValid hinput hsize
        have hfinal := initializeGramSchmidtLoop_processed_valid matrix 1
          initialMu seededNorms loopOutput hinput hinitial hloop
        refine ⟨hfinal.shape.norms_size, hfinal.shape.mu_size,
          hinput.rows_square, hfinal.shape.mu_rows, ?_, ?_⟩
        · intro index hindex
          have hpositive := hfinal.norms_positive index (by
            rw [hfinal.shape.norms_size] at hindex
            exact hindex)
          rw [getElem!_pos loopOutput.2 index hindex] at hpositive
          exact hpositive
        · intro rowCount hrowCount
          exact hfinal.gram_schmidt rowCount hrowCount hrowCount
  · rw [dif_neg hsize] at hrun
    contradiction

noncomputable def concreteLLLExecution :
    Generated.StrictRecombine.LLLExecution where
  inputValid := ConcreteLLLInputValid
  termination := concreteLLLTermination
  initialized_valid := by
    intro matrix mu norms transform hinput hrun
    exact initializeLLL_concrete_valid matrix mu norms transform hinput hrun
  output_input_valid := by
    intro initial output hvalid hrun
    exact (concreteLLLMainLoop_preserves_execution_valid initial output
      hvalid hrun).toInputValid
  output_size := by
    intro initial output hvalid hrun
    exact concreteLLLMainLoop_preserves_matrix_size initial output hvalid hrun

theorem gramPrefixDet_swap_preserved_of_before
    (matrix output : Generated.StrictRecombine.LLLMatrix) (left right rowCount : Nat)
    (hrun : Generated.StrictRecombine.swapMatrixRows matrix left right =
      .ok output)
    (hrowCount : rowCount ≤ matrix.size)
    (hbeforeLeft : rowCount ≤ left) (hbeforeRight : rowCount ≤ right) :
    Matrix.det (gramPrefixMatrix output rowCount) =
      Matrix.det (gramPrefixMatrix matrix rowCount) := by
  rw [gramPrefixMatrix_swap_eq_of_before matrix output left right rowCount
    hrun hrowCount hbeforeLeft hbeforeRight]

noncomputable def factorArrayToL2 (factors : Array SparsePolyZZ) :
    List (Polynomial Int) :=
  factors.toList.map SparsePolyZZ.toPoly

def ValidActiveIndices (active : Array Int32)
    (lifted : Array SparsePolyZZ) : Prop :=
  ∀ index (hindex : index < active.size),
    0 ≤ active[index]'hindex ∧
      (active[index]'hindex).toInt64.toNat < lifted.size

private def gatherSuffix (active : Array Int32)
    (lifted : Array SparsePolyZZ) (index : Nat) : Array SparsePolyZZ :=
  ((active.toList.drop index).map
    (fun sourceIndex => lifted[sourceIndex.toInt64.toNat]!)).toArray

theorem gatherActiveLoop_refines (active : Array Int32)
    (lifted : Array SparsePolyZZ) (index : Nat)
    (result : Array SparsePolyZZ)
    (hvalid : ValidActiveIndices active lifted) :
    Generated.StrictRecombine.gatherActiveLoop active lifted index result =
      .ok (result ++ gatherSuffix active lifted index) := by
  induction hmeasure : active.size - index using Nat.strong_induction_on
      generalizing index result with
  | h measure ih =>
      rw [Generated.StrictRecombine.gatherActiveLoop]
      split
      next hindex =>
        have hentry := hvalid index hindex
        rw [dif_pos hentry.1, dif_pos hentry.2]
        rw [ih (active.size - (index + 1)) (by omega)
          (index + 1) (result.push lifted[active[index].toInt64.toNat]) rfl]
        have hsuffix : active.toList.drop index = active[index] ::
            active.toList.drop (index + 1) := by
          simpa using List.drop_eq_getElem_cons
            (l := active.toList) (i := index) (by simpa using hindex)
        simp [gatherSuffix, hsuffix, Array.push, hentry.2]
      next hindex =>
        have hle : active.size ≤ index := Nat.le_of_not_gt hindex
        simp [gatherSuffix, List.drop_eq_nil_iff.mpr hle]

theorem gatherActive_refines (active : Array Int32)
    (lifted : Array SparsePolyZZ) (hvalid : ValidActiveIndices active lifted) :
    ∃ output,
      Generated.StrictRecombine.gatherActive active lifted = .ok output ∧
      factorArrayToL2 output = active.toList.map
        (fun sourceIndex => SparsePolyZZ.toPoly
          lifted[sourceIndex.toInt64.toNat]!) := by
  refine ⟨((active.toList.map
      (fun sourceIndex => lifted[sourceIndex.toInt64.toNat]!)).toArray), ?_, ?_⟩
  · simpa [Generated.StrictRecombine.gatherActive, gatherSuffix] using
      gatherActiveLoop_refines active lifted 0 #[] hvalid
  · simp [factorArrayToL2]

theorem gatherActive_size_of_success (active : Array Int32)
    (lifted output : Array SparsePolyZZ)
    (hrun : Generated.StrictRecombine.gatherActive active lifted = .ok output) :
    output.size = active.size := by
  unfold Generated.StrictRecombine.gatherActive at hrun
  have loopSize : ∀ index result,
      Generated.StrictRecombine.gatherActiveLoop active lifted index result =
        .ok output → output.size = result.size + (active.size - index) := by
    intro index result hloop
    induction hmeasure : active.size - index using Nat.strong_induction_on
        generalizing index result with
    | h measure ih =>
        rw [Generated.StrictRecombine.gatherActiveLoop] at hloop
        split at hloop
        next hindex =>
          dsimp at hloop
          split at hloop
          next hnonnegative =>
            split at hloop
            next hsource =>
              have htail := ih (active.size - (index + 1)) (by omega)
                (index + 1) (result.push lifted[active[index].toInt64.toNat])
                hloop rfl
              simp only [Array.size_push] at htail
              omega
            next hsource => simp at hloop
          next hnonnegative => simp at hloop
        next hindex =>
          have hle : active.size ≤ index := Nat.le_of_not_gt hindex
          have hout := Except.ok.inj hloop
          subst output
          have : measure = 0 := by omega
          simp [this]
  simpa using (loopSize 0 #[] hrun)

theorem cldPolysLoop_size_of_success
    (ops : Generated.StrictRecombine.CldRawOps)
    (fStar : SparsePolyZZ) (activeFactors : Array SparsePolyZZ)
    (modulus : ZZ) (index : Nat) (result output : Array SparsePolyZZ)
    (hrun : Generated.StrictRecombine.cldPolysLoop ops fStar activeFactors
      modulus index result = .ok output) :
    output.size = result.size + (activeFactors.size - index) := by
  induction hmeasure : activeFactors.size - index using Nat.strong_induction_on
      generalizing index result with
  | h measure ih =>
      rw [Generated.StrictRecombine.cldPolysLoop] at hrun
      split at hrun
      next hindex =>
        cases hdivmod : Generated.StrictHensel.__upoly_divmod_mod_raw_ir
            ops.divmodTermination fStar activeFactors[index] modulus with
        | error fault => simp [hdivmod] at hrun
        | ok division =>
            rcases division with ⟨quotient, remainder⟩
            simp only [hdivmod] at hrun
            split at hrun
            next hremainder =>
              cases hderivative : Generated.StrictRecombine.derivativeZZRaw
                  activeFactors[index] with
              | error fault => simp [hderivative] at hrun
              | ok derivativeRaw =>
                  simp only [hderivative] at hrun
                  cases hmod : Generated.StrictRecombine.modCoeffLoop
                      derivativeRaw modulus 0 #[] with
                  | error fault => simp [hmod] at hrun
                  | ok derivativeMod =>
                      simp only [hmod] at hrun
                      cases hproduct :
                          Generated.StrictRecombine.multiplyNormalizeModRaw
                            quotient derivativeMod modulus with
                      | error fault => simp [hproduct] at hrun
                      | ok product =>
                          simp only [hproduct] at hrun
                          cases hsymmetric :
                              Generated.StrictRecombine.symmetricModRaw product
                                modulus with
                          | error fault => simp [hsymmetric] at hrun
                          | ok cld =>
                              simp only [hsymmetric] at hrun
                              have htail := ih
                                (activeFactors.size - (index + 1)) (by omega)
                                (index + 1) (result.push cld) hrun rfl
                              simp only [Array.size_push] at htail
                              omega
            next hremainder => contradiction
      next hindex =>
        have hle : activeFactors.size ≤ index := Nat.le_of_not_gt hindex
        have hout := Except.ok.inj hrun
        subst output
        omega

theorem cldPolys_size_of_success
    (ops : Generated.StrictRecombine.CldRawOps)
    (fStar : SparsePolyZZ) (activeFactors output : Array SparsePolyZZ)
    (modulus : ZZ)
    (hrun : Generated.StrictRecombine.cldPolys ops fStar activeFactors
      modulus = .ok output) :
    output.size = activeFactors.size := by
  unfold Generated.StrictRecombine.cldPolys at hrun
  simpa using cldPolysLoop_size_of_success ops fStar activeFactors modulus 0
    #[] output hrun

/-- The generated class collector only updates existing result slots, so its
physical outer array size is invariant. -/
theorem collectCandidateClasses_size
    (classes : Array (Option Nat)) (column : Nat)
    (result output : Array (Array Int32))
    (hrun : Generated.StrictRecombine.collectCandidateClasses classes column
      result = .ok output) :
    output.size = result.size := by
  induction hmeasure : classes.size - column using Nat.strong_induction_on
      generalizing column result output with
  | h measure ih =>
      rw [Generated.StrictRecombine.collectCandidateClasses] at hrun
      split at hrun
      next hcolumn =>
        split at hrun
        next hnone => contradiction
        next classId hclassId =>
          split at hrun
          next hclass =>
            have htail := ih (classes.size - (column + 1)) (by omega)
              (column := column + 1)
              (result := result.set classId
                (result[classId].push column.toUInt32.toInt32))
              (output := output) hrun rfl
            simpa using htail
          next hclass => contradiction
      next hcolumn =>
        exact congrArg Array.size (Except.ok.inj hrun.symm)

/-- During the generated outer partition loop, at most one new class is
created per unprocessed physical column. -/
theorem partitionCandidateColumns_classCount_bound
    (transform : Generated.StrictRecombine.LLLMatrix)
    (shortRows : Array Nat) (factorCount column classCount : Nat)
    (classes outputClasses : Array (Option Nat)) (outputCount : Nat)
    (hcolumn : column ≤ factorCount)
    (hrun : Generated.StrictRecombine.partitionCandidateColumns transform
      shortRows factorCount column classCount classes =
        .ok (outputClasses, outputCount)) :
    outputCount ≤ classCount + (factorCount - column) := by
  induction hmeasure : factorCount - column using Nat.strong_induction_on
      generalizing column classCount classes outputClasses outputCount with
  | h measure ih =>
      rw [Generated.StrictRecombine.partitionCandidateColumns] at hrun
      split at hrun
      next hnext =>
        split at hrun
        next hclass =>
          split at hrun
          next existing =>
            have htail := ih (factorCount - (column + 1)) (by omega)
              (column := column + 1) (classCount := classCount)
              (classes := classes) (outputClasses := outputClasses)
              (outputCount := outputCount) (by omega) hrun rfl
            omega
          next hnone =>
            cases hassign : Generated.StrictRecombine.assignCandidateClass
                transform shortRows column factorCount classCount (column + 1)
                (classes.set column (some classCount)) with
            | error fault =>
              simp [hassign] at hrun
            | ok assigned =>
              simp only [hassign] at hrun
              have htail := ih (factorCount - (column + 1)) (by omega)
                (column := column + 1) (classCount := classCount + 1)
                (classes := assigned) (outputClasses := outputClasses)
                (outputCount := outputCount) (by omega) hrun rfl
              omega
        next hclass => contradiction
      next hnext =>
        have hout := Except.ok.inj hrun
        cases hout
        omega

/-- Every successful literal `__extract_candidates` execution returns no
more candidate classes than there are physical active-factor columns. -/
theorem extractCandidates_size_le
    (shortRows : Array Nat)
    (transform : Generated.StrictRecombine.LLLMatrix) (factorCount : Nat)
    (output : Array (Array Int32))
    (hrun : Generated.StrictRecombine.extractCandidates shortRows transform
      factorCount = .ok output) :
    output.size ≤ factorCount := by
  by_cases hempty : shortRows.isEmpty
  · simp [Generated.StrictRecombine.extractCandidates, hempty] at hrun
    subst output
    simp
  · cases hpartition : Generated.StrictRecombine.partitionCandidateColumns
        transform shortRows factorCount 0 0
          (Array.replicate factorCount none) with
    | error fault =>
      simp [Generated.StrictRecombine.extractCandidates, hempty,
        hpartition] at hrun
    | ok partition =>
      rcases partition with ⟨classes, classCount⟩
      simp only [Generated.StrictRecombine.extractCandidates, hempty,
        hpartition] at hrun
      have hcount := partitionCandidateColumns_classCount_bound transform
        shortRows factorCount 0 0 (Array.replicate factorCount none) classes
        classCount (by omega) hpartition
      have hsize := collectCandidateClasses_size classes 0
        (Array.replicate classCount #[]) output hrun
      simp at hsize
      simp at hcount
      exact hsize.le.trans hcount

theorem appendFallbackLoop_refines (fallback : Array SparsePolyZZ)
    (index : Nat) (result : Array SparsePolyZZ) :
    Generated.StrictRecombine.appendFallbackLoop fallback index result =
      .ok (result ++ (fallback.toList.drop index).toArray) := by
  induction hmeasure : fallback.size - index using Nat.strong_induction_on
      generalizing index result with
  | h measure ih =>
      rw [Generated.StrictRecombine.appendFallbackLoop]
      split
      next hindex =>
        rw [ih (fallback.size - (index + 1)) (by omega)
          (index + 1) (result.push fallback[index]) rfl]
        have hsuffix : fallback.toList.drop index = fallback[index] ::
            fallback.toList.drop (index + 1) := by
          simpa using List.drop_eq_getElem_cons
            (l := fallback.toList) (i := index) (by simpa using hindex)
        simp [hsuffix, Array.push]
      next hindex =>
        have hle : fallback.size ≤ index := Nat.le_of_not_gt hindex
        simp [List.drop_eq_nil_iff.mpr hle]

theorem appendFallback_refines (fallback result : Array SparsePolyZZ) :
    ∃ output,
      Generated.StrictRecombine.appendFallback fallback result = .ok output ∧
      factorArrayToL2 output =
        factorArrayToL2 result ++ factorArrayToL2 fallback := by
  refine ⟨result ++ fallback, ?_, ?_⟩
  · simpa [Generated.StrictRecombine.appendFallback] using
      appendFallbackLoop_refines fallback 0 result
  · simp [factorArrayToL2]

/-- Pure list meaning of the source reverse-erasure loop. -/
def removeConsumedL2 (active : Array Int32) (consumed : Array Bool) :
    List Int32 :=
  (active.toList.zip consumed.toList).filterMap fun item =>
    if item.2 then none else some item.1

private theorem removeConsumedLoop_refines_prefix
    (consumed : Array Bool) (remaining : Nat) (active : Array Int32)
    (hremaining : remaining ≤ consumed.size)
    (hactiveRemaining : remaining ≤ active.size) :
    ∃ output,
      Generated.StrictRecombine.removeConsumedLoop consumed remaining active =
        .ok output ∧
      output.size ≤ active.size := by
  induction remaining generalizing active with
  | zero =>
      refine ⟨active, ?_, Nat.le_refl _⟩
      rw [Generated.StrictRecombine.removeConsumedLoop]
  | succ remaining ih =>
      rw [Generated.StrictRecombine.removeConsumedLoop]
      have hindex : remaining < consumed.size := by omega
      rw [dif_pos hindex]
      split
      next hmarked =>
        have hactive : remaining < active.size := by omega
        rw [dif_pos hactive]
        have heraseSize : (active.eraseIdxIfInBounds remaining).size =
            active.size - 1 := by
          simp [hactive]
        have hremaining' : remaining ≤
            (active.eraseIdxIfInBounds remaining).size := by
          rw [heraseSize]
          omega
        rcases ih (active.eraseIdxIfInBounds remaining) (by omega) hremaining' with
          ⟨output, hrun, hsize⟩
        exact ⟨output, hrun, hsize.trans (by rw [heraseSize]; omega)⟩
      next hkept =>
        rcases ih active (by omega) (by omega) with ⟨output, hrun, hsize⟩
        exact ⟨output, hrun, hsize⟩

private theorem removeConsumedLoop_strict_of_marked
    (consumed : Array Bool) (remaining : Nat) (active output : Array Int32)
    (hremaining : remaining ≤ consumed.size)
    (hactiveRemaining : remaining ≤ active.size)
    (hmarked : ∃ index, ∃ hindex : index < remaining,
      consumed[index] = true)
    (hrun : Generated.StrictRecombine.removeConsumedLoop consumed remaining active =
      .ok output) :
    output.size < active.size := by
  induction remaining generalizing active output with
  | zero =>
      rcases hmarked with ⟨index, hindex, _⟩
      omega
  | succ remaining ih =>
      rw [Generated.StrictRecombine.removeConsumedLoop] at hrun
      have hindexBound : remaining < consumed.size := by omega
      rw [dif_pos hindexBound] at hrun
      by_cases hlast : consumed[remaining] = true
      · rw [if_pos hlast] at hrun
        have hactive : remaining < active.size := by omega
        rw [dif_pos hactive] at hrun
        have htailLe : output.size ≤
            (active.eraseIdxIfInBounds remaining).size := by
          rcases removeConsumedLoop_refines_prefix consumed remaining
              (active.eraseIdxIfInBounds remaining) (by omega) (by
                simp [hactive]; omega) with ⟨tail, htailRun, htailSize⟩
          rw [hrun] at htailRun
          exact Except.ok.inj htailRun ▸ htailSize
        simpa [hactive] using lt_of_le_of_lt htailLe (by
          simp [hactive]
          omega)
      · rw [if_neg hlast] at hrun
        have hmarkedPrefix : ∃ index, ∃ hindex : index < remaining,
            consumed[index] = true := by
          rcases hmarked with ⟨index, hindex, hvalue⟩
          refine ⟨index, ?_, hvalue⟩
          have hne : index ≠ remaining := by
            intro heq
            subst index
            exact hlast hvalue
          omega
        exact ih active output (by omega) (by omega) hmarkedPrefix hrun

theorem removeConsumed_succeeds (active : Array Int32)
    (consumed : Array Bool) (hsizes : consumed.size = active.size) :
    ∃ output,
      Generated.StrictRecombine.removeConsumed active consumed = .ok output ∧
      output.size ≤ active.size := by
  unfold Generated.StrictRecombine.removeConsumed
  rw [dif_pos hsizes]
  exact removeConsumedLoop_refines_prefix consumed active.size active
    (by omega) (Nat.le_refl _)

theorem removeConsumed_strict_of_marked (active : Array Int32)
    (consumed : Array Bool) (hsizes : consumed.size = active.size)
    (hmarked : ∃ index, ∃ hindex : index < consumed.size,
      consumed[index] = true)
    (output : Array Int32)
    (hrun : Generated.StrictRecombine.removeConsumed active consumed =
      .ok output) :
    output.size < active.size := by
  unfold Generated.StrictRecombine.removeConsumed at hrun
  rw [dif_pos hsizes] at hrun
  apply removeConsumedLoop_strict_of_marked consumed active.size active output
    (by omega) (Nat.le_refl _) (by simpa [hsizes] using hmarked) hrun

def CandidateIndicesValid (candidate : Array Int32)
    (consumed : Array Bool) : Prop :=
  ∀ index (hindex : index < candidate.size),
    0 ≤ candidate[index] ∧
      candidate[index].toInt64.toNat < consumed.size

/-- The checked conversion used by the generated Zassenhaus attempt produces
valid `Int32` lookup indices whenever the natural candidate is bounded by the
active array.  Both nonnegativity and round-trip equality are proved from the
actual `2^31` check. -/
theorem combinationToInt32_candidate_valid (indices : Array Nat)
    (activeSize : Nat) (output : Array Int32)
    (hbound : ∀ position (hposition : position < indices.size),
      indices[position] < activeSize)
    (hactiveFits : activeSize ≤ 2 ^ 31)
    (hrun : Generated.StrictRecombine.combinationToInt32 indices =
      .ok output) :
    CandidateIndicesValid output (Array.replicate activeSize false) := by
  have hfits : ∀ position (hposition : position < indices.size),
      indices[position] < 2 ^ 31 := by
    intro position hposition
    exact lt_of_lt_of_le (hbound position hposition) hactiveFits
  rcases combinationToInt32_toList indices hfits with
    ⟨expected, hexpected, hlist⟩
  have houtput : output = expected := by
    exact Except.ok.inj (hrun.symm.trans hexpected)
  subst expected
  have hsize : output.size = indices.size := by
    simpa using congrArg List.length hlist
  intro position hposition
  have hpositionIndices : position < indices.size := by omega
  have hentryDependent : output[position] =
      indices[position].toUInt32.toInt32 := by
    have hget := congrArg (fun values : List Int32 => values[position]!) hlist
    change output.toList[position]! =
      (indices.toList.map fun value => value.toUInt32.toInt32)[position]! at hget
    rw [getElem!_pos output.toList position (by simpa using hposition),
      Array.getElem_toList] at hget
    have hmap := getElem!_map_of_lt
      (fun value : Nat => value.toUInt32.toInt32) indices.toList position
      (by simpa using hpositionIndices)
    rw [hmap] at hget
    simpa [getElem!_pos indices.toList position (by simpa using hpositionIndices),
      Array.getElem_toList] using hget
  rw [hentryDependent]
  have hroundtrip := nat_toUInt32_toInt32_nonnegative_and_toNat
    indices[position] (hfits position hpositionIndices)
  refine ⟨hroundtrip.1, ?_⟩
  rw [hroundtrip.2]
  simpa using hbound position hpositionIndices

theorem legal_combination_toInt32_candidate_valid {upper count : Nat}
    (indices : Array Nat) (output : Array Int32)
    (hlegal : LegalCombination upper count indices)
    (hupperFits : upper ≤ 2 ^ 31)
    (hrun : Generated.StrictRecombine.combinationToInt32 indices =
      .ok output) :
    CandidateIndicesValid output (Array.replicate upper false) := by
  apply combinationToInt32_candidate_valid indices upper output
  · intro position hposition
    rw [← getElem!_pos indices position hposition]
    exact hlegal.2.2 position hposition
  · exact hupperFits
  · exact hrun

theorem candidateAvailableLoop_succeeds (candidate : Array Int32)
    (consumed : Array Bool) (index : Nat)
    (hvalid : CandidateIndicesValid candidate consumed) :
    ∃ available,
      Generated.StrictRecombine.candidateAvailableLoop candidate consumed index =
        .ok available := by
  induction hmeasure : candidate.size - index using Nat.strong_induction_on
      generalizing index with
  | h measure ih =>
      rw [Generated.StrictRecombine.candidateAvailableLoop]
      split
      next hindex =>
        have hentry := hvalid index hindex
        rw [dif_pos hentry.1, dif_pos hentry.2]
        split
        next hconsumed => exact ⟨false, rfl⟩
        next hfree =>
          exact ih (candidate.size - (index + 1)) (by omega)
            (index + 1) rfl
      next hindex => exact ⟨true, rfl⟩

theorem markConsumedLoop_succeeds_size (candidate : Array Int32)
    (consumed : Array Bool) (index : Nat)
    (hvalid : CandidateIndicesValid candidate consumed) :
    ∃ output,
      Generated.StrictRecombine.markConsumedLoop candidate index consumed =
        .ok output ∧ output.size = consumed.size := by
  induction hmeasure : candidate.size - index using Nat.strong_induction_on
      generalizing index consumed with
  | h measure ih =>
      rw [Generated.StrictRecombine.markConsumedLoop]
      split
      next hindex =>
        have hentry := hvalid index hindex
        let activeNat := candidate[index].toInt64.toNat
        rw [dif_pos hentry.1, dif_pos hentry.2]
        have hvalidSet : CandidateIndicesValid candidate
            (consumed.set activeNat true) := by
          intro next hnext
          have hv := hvalid next hnext
          exact ⟨hv.1, by simpa using hv.2⟩
        rcases ih (candidate.size - (index + 1)) (by omega)
            (consumed.set activeNat true) (index + 1) hvalidSet rfl with
          ⟨output, hrun, hsize⟩
        refine ⟨output, hrun, ?_⟩
        simpa using hsize
      next hindex => exact ⟨consumed, rfl, rfl⟩

theorem markConsumedLoop_size_of_success (candidate : Array Int32)
    (index : Nat) (consumed output : Array Bool)
    (hrun : Generated.StrictRecombine.markConsumedLoop candidate index consumed =
      .ok output) :
    output.size = consumed.size := by
  induction hmeasure : candidate.size - index using Nat.strong_induction_on
      generalizing index consumed output with
  | h measure ih =>
      rw [Generated.StrictRecombine.markConsumedLoop] at hrun
      split at hrun
      next hindex =>
        dsimp at hrun
        split at hrun
        next hnonnegative =>
          split at hrun
          next hactive =>
            have htail := ih (candidate.size - (index + 1)) (by omega)
              (index + 1)
              (consumed.set candidate[index].toInt64.toNat true) output
              hrun rfl
            simpa using htail
          next hactive => contradiction
        next hnonnegative => contradiction
      next hindex => exact Except.ok.inj hrun ▸ rfl

theorem validateCandidatesLoop_consumed_size
    (ops : Generated.StrictRecombine.CandidateValidationRawOps)
    (candidates : Array (Array Int32)) (candidateIndex : Nat)
    (activeLifted : Array SparsePolyZZ) (modulus : ZZ)
    (fStar fStar' : SparsePolyZZ) (result result' : Array SparsePolyZZ)
    (consumed consumed' : Array Bool) (remaining : Nat)
    (hrun : Generated.StrictRecombine.validateCandidatesLoop ops candidates
      candidateIndex activeLifted modulus fStar result consumed remaining =
        .ok (fStar', result', consumed')) :
    consumed'.size = consumed.size := by
  induction hmeasure : candidates.size - candidateIndex using Nat.strong_induction_on
      generalizing candidateIndex fStar result consumed remaining fStar' result' consumed' with
  | h measure ih =>
      rw [Generated.StrictRecombine.validateCandidatesLoop] at hrun
      split at hrun
      next hcandidates =>
        dsimp at hrun
        split at hrun
        next hempty =>
          exact ih (candidates.size - (candidateIndex + 1)) (by omega)
            (candidateIndex := candidateIndex + 1) (fStar := fStar)
            (result := result) (consumed := consumed) (remaining := remaining)
            (fStar' := fStar') (result' := result') (consumed' := consumed') hrun rfl
        next hempty =>
          split at hrun
          next htrivial =>
            exact ih (candidates.size - (candidateIndex + 1)) (by omega)
              (candidateIndex := candidateIndex + 1) (fStar := fStar)
              (result := result) (consumed := consumed) (remaining := remaining)
              (fStar' := fStar') (result' := result') (consumed' := consumed') hrun rfl
          next hnontrivial =>
            cases havailable : Generated.StrictRecombine.candidateAvailable
                candidates[candidateIndex] consumed with
            | error fault => simp [havailable] at hrun
            | ok available =>
              cases available with
              | false =>
                simp only [havailable] at hrun
                exact ih (candidates.size - (candidateIndex + 1)) (by omega)
                  (candidateIndex := candidateIndex + 1) (fStar := fStar)
                  (result := result) (consumed := consumed) (remaining := remaining)
                  (fStar' := fStar') (result' := result')
                  (consumed' := consumed') hrun rfl
              | true =>
                simp only [havailable] at hrun
                split at hrun
                next hfstar =>
                  cases hproduct : Generated.StrictRecombine.trialProductLoop
                      ops.product candidates[candidateIndex] activeLifted modulus 0
                      #[(⟨0⟩, fStar[0].2)] with
                  | error fault => simp [hproduct] at hrun
                  | ok product =>
                    simp only [hproduct] at hrun
                    cases hsymmetric : Generated.StrictRecombine.symmetricModRaw
                        product modulus with
                    | error fault => simp [hsymmetric] at hrun
                    | ok symmetric =>
                      simp only [hsymmetric] at hrun
                      cases hprimitive : Generated.StrictRecombine.primitiveRaw symmetric with
                      | error fault => simp [hprimitive] at hrun
                      | ok primitiveResult =>
                        rcases primitiveResult with ⟨content, factor⟩
                        simp only [hprimitive] at hrun
                        cases hdivmod : Generated.StrictRecombine.exactDivmodRaw
                            fStar factor with
                        | error fault => simp [hdivmod] at hrun
                        | ok divResult =>
                          rcases divResult with ⟨quotient, remainder⟩
                          simp only [hdivmod] at hrun
                          by_cases hremainder : remainder.isEmpty = true
                          · simp only [hremainder, if_true] at hrun
                            cases hquotientPrimitive :
                                Generated.StrictRecombine.primitiveRaw quotient with
                            | error fault => simp [hquotientPrimitive] at hrun
                            | ok quotientResult =>
                              rcases quotientResult with ⟨quotientContent,
                                quotientPrimitive⟩
                              simp only [hquotientPrimitive] at hrun
                              cases hmark : Generated.StrictRecombine.markConsumedLoop
                                  candidates[candidateIndex] 0 consumed with
                              | error fault => simp [hmark] at hrun
                              | ok consumedNext =>
                                simp only [hmark] at hrun
                                have htail := ih
                                  (candidates.size - (candidateIndex + 1))
                                  (by omega) (candidateIndex := candidateIndex + 1)
                                  (fStar := quotientPrimitive)
                                  (result := result.push factor)
                                  (consumed := consumedNext)
                                  (remaining := remaining - candidates[candidateIndex].size)
                                  (fStar' := fStar') (result' := result')
                                  (consumed' := consumed') hrun rfl
                                exact htail.trans
                                  (markConsumedLoop_size_of_success
                                    candidates[candidateIndex] 0 consumed
                                    consumedNext hmark)
                          · simp only [hremainder, if_false] at hrun
                            exact ih (candidates.size - (candidateIndex + 1))
                              (by omega) (candidateIndex := candidateIndex + 1)
                              (fStar := fStar) (result := result)
                              (consumed := consumed) (remaining := remaining)
                              (fStar' := fStar') (result' := result')
                              (consumed' := consumed') hrun rfl
                next hfstar => contradiction
      next hcandidates =>
        exact (congrArg (fun output => output.2.2.size)
          (Except.ok.inj hrun)).symm

/-- The literal validation loop appends at most one factor for each physical
candidate slot it scans.  This follows the actual success/rejection branches,
including exact division and consumed marking. -/
theorem validateCandidatesLoop_result_size_le
    (ops : Generated.StrictRecombine.CandidateValidationRawOps)
    (candidates : Array (Array Int32)) (candidateIndex : Nat)
    (activeLifted : Array SparsePolyZZ) (modulus : ZZ)
    (fStar fStar' : SparsePolyZZ) (result result' : Array SparsePolyZZ)
    (consumed consumed' : Array Bool) (remaining : Nat)
    (hrun : Generated.StrictRecombine.validateCandidatesLoop ops candidates
      candidateIndex activeLifted modulus fStar result consumed remaining =
        .ok (fStar', result', consumed')) :
    result'.size ≤ result.size + (candidates.size - candidateIndex) := by
  induction hmeasure : candidates.size - candidateIndex using Nat.strong_induction_on
      generalizing candidateIndex fStar result consumed remaining fStar' result' consumed' with
  | h measure ih =>
      rw [Generated.StrictRecombine.validateCandidatesLoop] at hrun
      split at hrun
      next hcandidates =>
        dsimp at hrun
        split at hrun
        next hempty =>
          have htail := ih (candidates.size - (candidateIndex + 1)) (by omega)
            (candidateIndex := candidateIndex + 1) (fStar := fStar)
            (result := result) (consumed := consumed) (remaining := remaining)
            (fStar' := fStar') (result' := result') (consumed' := consumed')
            hrun rfl
          omega
        next hempty =>
          split at hrun
          next htrivial =>
            have htail := ih (candidates.size - (candidateIndex + 1)) (by omega)
              (candidateIndex := candidateIndex + 1) (fStar := fStar)
              (result := result) (consumed := consumed) (remaining := remaining)
              (fStar' := fStar') (result' := result') (consumed' := consumed')
              hrun rfl
            omega
          next hnontrivial =>
            cases havailable : Generated.StrictRecombine.candidateAvailable
                candidates[candidateIndex] consumed with
            | error fault => simp [havailable] at hrun
            | ok available =>
              cases available with
              | false =>
                simp only [havailable] at hrun
                have htail := ih (candidates.size - (candidateIndex + 1))
                  (by omega) (candidateIndex := candidateIndex + 1)
                  (fStar := fStar) (result := result) (consumed := consumed)
                  (remaining := remaining) (fStar' := fStar')
                  (result' := result') (consumed' := consumed') hrun rfl
                omega
              | true =>
                simp only [havailable] at hrun
                split at hrun
                next hfstar =>
                  cases hproduct : Generated.StrictRecombine.trialProductLoop
                      ops.product candidates[candidateIndex] activeLifted modulus 0
                      #[(⟨0⟩, fStar[0].2)] with
                  | error fault => simp [hproduct] at hrun
                  | ok product =>
                    simp only [hproduct] at hrun
                    cases hsymmetric : Generated.StrictRecombine.symmetricModRaw
                        product modulus with
                    | error fault => simp [hsymmetric] at hrun
                    | ok symmetric =>
                      simp only [hsymmetric] at hrun
                      cases hprimitive : Generated.StrictRecombine.primitiveRaw
                          symmetric with
                      | error fault => simp [hprimitive] at hrun
                      | ok primitiveResult =>
                        rcases primitiveResult with ⟨content, factor⟩
                        simp only [hprimitive] at hrun
                        cases hdivmod : Generated.StrictRecombine.exactDivmodRaw
                            fStar factor with
                        | error fault => simp [hdivmod] at hrun
                        | ok divResult =>
                          rcases divResult with ⟨quotient, remainder⟩
                          simp only [hdivmod] at hrun
                          by_cases hremainder : remainder.isEmpty = true
                          · simp only [hremainder, if_true] at hrun
                            cases hquotientPrimitive :
                                Generated.StrictRecombine.primitiveRaw quotient with
                            | error fault => simp [hquotientPrimitive] at hrun
                            | ok quotientResult =>
                              rcases quotientResult with ⟨quotientContent,
                                quotientPrimitive⟩
                              simp only [hquotientPrimitive] at hrun
                              cases hmark : Generated.StrictRecombine.markConsumedLoop
                                  candidates[candidateIndex] 0 consumed with
                              | error fault => simp [hmark] at hrun
                              | ok consumedNext =>
                                simp only [hmark] at hrun
                                have htail := ih
                                  (candidates.size - (candidateIndex + 1))
                                  (by omega) (candidateIndex := candidateIndex + 1)
                                  (fStar := quotientPrimitive)
                                  (result := result.push factor)
                                  (consumed := consumedNext)
                                  (remaining := remaining -
                                    candidates[candidateIndex].size)
                                  (fStar' := fStar') (result' := result')
                                  (consumed' := consumed') hrun rfl
                                simp only [Array.size_push] at htail
                                omega
                          · simp only [hremainder, if_false] at hrun
                            have htail := ih
                              (candidates.size - (candidateIndex + 1))
                              (by omega) (candidateIndex := candidateIndex + 1)
                              (fStar := fStar) (result := result)
                              (consumed := consumed) (remaining := remaining)
                              (fStar' := fStar') (result' := result')
                              (consumed' := consumed') hrun rfl
                            omega
                next hfstar => contradiction
      next hcandidates =>
        have hout := Except.ok.inj hrun
        have hresult := congrArg (fun output => output.2.1.size) hout
        simp only at hresult
        omega

/-- A complete generated validation call returns no more new factors than
there are concrete candidate arrays. -/
theorem validateCandidates_result_size_le
    (ops : Generated.StrictRecombine.CandidateValidationRawOps)
    (fStar : SparsePolyZZ) (activeLifted : Array SparsePolyZZ) (modulus : ZZ)
    (candidates : Array (Array Int32)) (result : Array SparsePolyZZ)
    (fStar' : SparsePolyZZ) (result' : Array SparsePolyZZ)
    (consumed : Array Bool)
    (hrun : Generated.StrictRecombine.validateCandidates ops fStar activeLifted
      modulus candidates result = .ok (fStar', result', consumed)) :
    result'.size ≤ result.size + candidates.size := by
  unfold Generated.StrictRecombine.validateCandidates at hrun
  simpa using validateCandidatesLoop_result_size_le ops candidates 0
    activeLifted modulus fStar fStar' result result'
    (Array.replicate activeLifted.size false) consumed activeLifted.size hrun

/-- Candidate extraction followed by the literal validation loop can append
no more factors than the number of physical active-factor columns. -/
theorem extractAndValidate_result_size_le
    (ops : Generated.StrictRecombine.CandidateValidationRawOps)
    (shortRows : Array Nat)
    (transform : Generated.StrictRecombine.LLLMatrix) (factorCount : Nat)
    (candidates : Array (Array Int32))
    (fStar : SparsePolyZZ) (activeLifted : Array SparsePolyZZ) (modulus : ZZ)
    (result : Array SparsePolyZZ) (fStar' : SparsePolyZZ)
    (result' : Array SparsePolyZZ) (consumed : Array Bool)
    (hextract : Generated.StrictRecombine.extractCandidates shortRows transform
      factorCount = .ok candidates)
    (hvalidate : Generated.StrictRecombine.validateCandidates ops fStar
      activeLifted modulus candidates result =
        .ok (fStar', result', consumed)) :
    result'.size ≤ result.size + factorCount := by
  exact (validateCandidates_result_size_le ops fStar activeLifted modulus
    candidates result fStar' result' consumed hvalidate).trans
      (Nat.add_le_add_left
        (extractCandidates_size_le shortRows transform factorCount candidates
          hextract) result.size)

/-- Any successful literal long-division run on a nonempty dividend has
physically entered the positive-size divisor branch. -/
theorem exactDivmodRaw_divisor_nonempty_of_success
    (dividend divisor quotient remainder : SparsePolyZZ)
    (hdividend : 0 < dividend.size)
    (hrun : Generated.StrictRecombine.exactDivmodRaw dividend divisor =
      .ok (quotient, remainder)) :
    0 < divisor.size := by
  unfold Generated.StrictRecombine.exactDivmodRaw at hrun
  rw [Generated.StrictRecombine.exactDivmodLoop, dif_pos hdividend] at hrun
  split at hrun
  next hdivisor => exact hdivisor
  next hdivisor => contradiction

/-- Every factor physically present before validation or appended by a
successful exact-division branch is represented by a nonempty sparse array. -/
theorem validateCandidatesLoop_result_nonempty
    (ops : Generated.StrictRecombine.CandidateValidationRawOps)
    (candidates : Array (Array Int32)) (candidateIndex : Nat)
    (activeLifted : Array SparsePolyZZ) (modulus : ZZ)
    (fStar fStar' : SparsePolyZZ) (result result' : Array SparsePolyZZ)
    (consumed consumed' : Array Bool) (remaining : Nat)
    (hresult : ∀ factor ∈ result.toList, 0 < factor.size)
    (hrun : Generated.StrictRecombine.validateCandidatesLoop ops candidates
      candidateIndex activeLifted modulus fStar result consumed remaining =
        .ok (fStar', result', consumed')) :
    ∀ factor ∈ result'.toList, 0 < factor.size := by
  induction hmeasure : candidates.size - candidateIndex using Nat.strong_induction_on
      generalizing candidateIndex fStar result consumed remaining fStar' result' consumed' with
  | h measure ih =>
      rw [Generated.StrictRecombine.validateCandidatesLoop] at hrun
      split at hrun
      next hcandidates =>
        dsimp at hrun
        split at hrun
        next hempty =>
          exact ih (candidates.size - (candidateIndex + 1)) (by omega)
            (candidateIndex := candidateIndex + 1) (fStar := fStar)
            (result := result) (consumed := consumed) (remaining := remaining)
            (fStar' := fStar') (result' := result') (consumed' := consumed')
            hresult hrun rfl
        next hempty =>
          split at hrun
          next htrivial =>
            exact ih (candidates.size - (candidateIndex + 1)) (by omega)
              (candidateIndex := candidateIndex + 1) (fStar := fStar)
              (result := result) (consumed := consumed) (remaining := remaining)
              (fStar' := fStar') (result' := result') (consumed' := consumed')
              hresult hrun rfl
          next hnontrivial =>
            cases havailable : Generated.StrictRecombine.candidateAvailable
                candidates[candidateIndex] consumed with
            | error fault => simp [havailable] at hrun
            | ok available =>
              cases available with
              | false =>
                simp only [havailable] at hrun
                exact ih (candidates.size - (candidateIndex + 1)) (by omega)
                  (candidateIndex := candidateIndex + 1) (fStar := fStar)
                  (result := result) (consumed := consumed)
                  (remaining := remaining) (fStar' := fStar')
                  (result' := result') (consumed' := consumed') hresult hrun rfl
              | true =>
                simp only [havailable] at hrun
                split at hrun
                next hfstar =>
                  cases hproduct : Generated.StrictRecombine.trialProductLoop
                      ops.product candidates[candidateIndex] activeLifted modulus 0
                      #[(⟨0⟩, fStar[0].2)] with
                  | error fault => simp [hproduct] at hrun
                  | ok product =>
                    simp only [hproduct] at hrun
                    cases hsymmetric : Generated.StrictRecombine.symmetricModRaw
                        product modulus with
                    | error fault => simp [hsymmetric] at hrun
                    | ok symmetric =>
                      simp only [hsymmetric] at hrun
                      cases hprimitive : Generated.StrictRecombine.primitiveRaw
                          symmetric with
                      | error fault => simp [hprimitive] at hrun
                      | ok primitiveResult =>
                        rcases primitiveResult with ⟨content, factor⟩
                        simp only [hprimitive] at hrun
                        cases hdivmod : Generated.StrictRecombine.exactDivmodRaw
                            fStar factor with
                        | error fault => simp [hdivmod] at hrun
                        | ok divResult =>
                          rcases divResult with ⟨quotient, remainder⟩
                          simp only [hdivmod] at hrun
                          by_cases hremainder : remainder.isEmpty = true
                          · simp only [hremainder, if_true] at hrun
                            cases hquotientPrimitive :
                                Generated.StrictRecombine.primitiveRaw quotient with
                            | error fault => simp [hquotientPrimitive] at hrun
                            | ok quotientResult =>
                              rcases quotientResult with ⟨quotientContent,
                                quotientPrimitive⟩
                              simp only [hquotientPrimitive] at hrun
                              cases hmark : Generated.StrictRecombine.markConsumedLoop
                                  candidates[candidateIndex] 0 consumed with
                              | error fault => simp [hmark] at hrun
                              | ok consumedNext =>
                                simp only [hmark] at hrun
                                have hfactor : 0 < factor.size :=
                                  exactDivmodRaw_divisor_nonempty_of_success
                                    fStar factor quotient remainder hfstar hdivmod
                                have hresultPush : ∀ candidateFactor ∈
                                    (result.push factor).toList,
                                    0 < candidateFactor.size := by
                                  intro candidateFactor hcandidateFactor
                                  rw [Array.toList_push] at hcandidateFactor
                                  rcases List.mem_append.mp hcandidateFactor with
                                    hprefix | hlast
                                  · exact hresult candidateFactor hprefix
                                  · have hsame : candidateFactor = factor := by
                                      simpa using hlast
                                    subst candidateFactor
                                    exact hfactor
                                exact ih
                                  (candidates.size - (candidateIndex + 1))
                                  (by omega) (candidateIndex := candidateIndex + 1)
                                  (fStar := quotientPrimitive)
                                  (result := result.push factor)
                                  (consumed := consumedNext)
                                  (remaining := remaining -
                                    candidates[candidateIndex].size)
                                  (fStar' := fStar') (result' := result')
                                  (consumed' := consumed') hresultPush hrun rfl
                          · simp only [hremainder, if_false] at hrun
                            exact ih
                              (candidates.size - (candidateIndex + 1))
                              (by omega) (candidateIndex := candidateIndex + 1)
                              (fStar := fStar) (result := result)
                              (consumed := consumed) (remaining := remaining)
                              (fStar' := fStar') (result' := result')
                              (consumed' := consumed') hresult hrun rfl
                next hfstar => contradiction
      next hcandidates =>
        have hout := Except.ok.inj hrun
        cases hout
        exact hresult

theorem validateCandidates_result_nonempty
    (ops : Generated.StrictRecombine.CandidateValidationRawOps)
    (fStar : SparsePolyZZ) (activeLifted : Array SparsePolyZZ) (modulus : ZZ)
    (candidates : Array (Array Int32)) (result : Array SparsePolyZZ)
    (fStar' : SparsePolyZZ) (result' : Array SparsePolyZZ)
    (consumed : Array Bool)
    (hresult : ∀ factor ∈ result.toList, 0 < factor.size)
    (hrun : Generated.StrictRecombine.validateCandidates ops fStar activeLifted
      modulus candidates result = .ok (fStar', result', consumed)) :
    ∀ factor ∈ result'.toList, 0 < factor.size := by
  unfold Generated.StrictRecombine.validateCandidates at hrun
  exact validateCandidatesLoop_result_nonempty ops candidates 0 activeLifted
    modulus fStar fStar' result result'
    (Array.replicate activeLifted.size false) consumed activeLifted.size
    hresult hrun

theorem validateCandidates_consumed_size
    (ops : Generated.StrictRecombine.CandidateValidationRawOps)
    (fStar : SparsePolyZZ) (activeLifted : Array SparsePolyZZ) (modulus : ZZ)
    (candidates : Array (Array Int32)) (result : Array SparsePolyZZ)
    (fStar' : SparsePolyZZ) (result' : Array SparsePolyZZ)
    (consumed : Array Bool)
    (hrun : Generated.StrictRecombine.validateCandidates ops fStar activeLifted
      modulus candidates result = .ok (fStar', result', consumed)) :
    consumed.size = activeLifted.size := by
  unfold Generated.StrictRecombine.validateCandidates at hrun
  exact (validateCandidatesLoop_consumed_size ops candidates 0 activeLifted
    modulus fStar fStar' result result'
    (Array.replicate activeLifted.size false) consumed activeLifted.size hrun).trans
      (by simp)

/-- The concrete gather and validation loops themselves prove that every
successful extraction removes at least one active entry.  No external length
or semantic termination oracle is required. -/
def removalTermination (ops : Generated.StrictRecombine.VanHoeijRawOps) :
    Generated.StrictRecombine.VanHoeijTermination ops := {
  extraction_decreases := by
    intro lifted modulus state activeLifted candidates fStar' result' consumed
      active' hgather hvalidate hfound hremove
    exact removeConsumed_strict_of_marked state.active consumed
      ((validateCandidates_consumed_size ops.validation state.fStar activeLifted
        modulus candidates state.result fStar' result' consumed hvalidate).trans
        (gatherActive_size_of_success state.active lifted activeLifted hgather))
      hfound active' hremove }

noncomputable def SelectedProductMod (modulus : Nat) (candidate : Array Int32)
    (activeLifted : Array SparsePolyZZ) (index : Nat) :
    Polynomial (ZMod modulus) :=
  ((candidate.toList.drop index).map fun activeIndex =>
    Refinement.StrictHensel.toPolyMod modulus
      activeLifted[activeIndex.toInt64.toNat]!).prod

/-- The semantic product read by the generated `Int32` trial loop is exactly
the source-order sublist product named by the pre-conversion natural
candidate.  This is occurrence-sensitive and uses the checked conversion's
actual round trip for every source position. -/
theorem selectedProductMod_combinationToInt32
    (modulus : Nat) (indices : Array Nat)
    (activeLifted : Array SparsePolyZZ) (candidate : Array Int32)
    (hbound : ∀ position (hposition : position < indices.size),
      indices[position] < activeLifted.size)
    (hactiveFits : activeLifted.size ≤ 2 ^ 31)
    (hrun : Generated.StrictRecombine.combinationToInt32 indices =
      .ok candidate) :
    SelectedProductMod modulus candidate activeLifted 0 =
      ((selectSourceIndices activeLifted.toList indices.toList).map
        (Refinement.StrictHensel.toPolyMod modulus)).prod := by
  have hfits : ∀ position (hposition : position < indices.size),
      indices[position] < 2 ^ 31 := by
    intro position hposition
    exact lt_of_lt_of_le (hbound position hposition) hactiveFits
  rcases combinationToInt32_toList indices hfits with
    ⟨expected, hexpected, hlist⟩
  have hcand : candidate = expected :=
    Except.ok.inj (hrun.symm.trans hexpected)
  subst expected
  unfold SelectedProductMod selectSourceIndices
  rw [List.drop_zero, hlist]
  simp only [List.map_map]
  congr 1
  apply List.map_congr_left
  intro value hvalue
  rcases List.mem_iff_getElem.mp hvalue with ⟨position, hposition, rfl⟩
  have hpositionArray : position < indices.size := by simpa using hposition
  have hroundtrip := nat_toUInt32_toInt32_nonnegative_and_toNat
    indices.toList[position] (by
      simpa [Array.getElem_toList] using hfits position hpositionArray)
  simp only [Function.comp_apply]
  rw [hroundtrip.2]
  congr 1
  have hactive : indices.toList[position] < activeLifted.size := by
    simpa [Array.getElem_toList] using hbound position hpositionArray
  rw [getElem!_pos activeLifted indices.toList[position] hactive,
    getElem!_pos activeLifted.toList indices.toList[position] (by simpa using hactive)]
  exact Array.getElem_toList hactive

private noncomputable def intTermsToPoly (terms : List (UMonomial × Int)) :
    Polynomial Int :=
  (terms.map fun term => Polynomial.monomial term.1.deg term.2).sum

theorem multiplyRowLoop_toPoly (left : UMonomial × Int)
    (right : SparsePolyZZ) (rightIndex : Nat) (terms : SparsePolyZZ) :
    intTermsToPoly (Generated.StrictRecombine.multiplyRowLoop
      left right rightIndex terms).toList =
      intTermsToPoly terms.toList +
        Polynomial.monomial left.1.deg left.2 *
          intTermsToPoly (right.toList.drop rightIndex) := by
  induction hmeasure : right.size - rightIndex using Nat.strong_induction_on
      generalizing rightIndex terms with
  | h measure ih =>
      rw [Generated.StrictRecombine.multiplyRowLoop]
      split
      next hright =>
        rw [ih (right.size - (rightIndex + 1)) (by omega)
          (rightIndex + 1) (terms.push
            (⟨left.1.deg + right[rightIndex].1.deg⟩,
              left.2 * right[rightIndex].2)) rfl]
        have hsuffix : right.toList.drop rightIndex = right[rightIndex] ::
            right.toList.drop (rightIndex + 1) := by
          simpa using List.drop_eq_getElem_cons
            (l := right.toList) (i := rightIndex) (by simpa using hright)
        simp [intTermsToPoly, hsuffix, add_comm, add_left_comm,
          mul_add]
        exact (Polynomial.monomial_mul_monomial _ _ _ _).symm
      next hright =>
        have hle : right.size ≤ rightIndex := Nat.le_of_not_gt hright
        simp [intTermsToPoly, List.drop_eq_nil_iff.mpr hle]

theorem multiplyTermsLoop_toPoly (left right : SparsePolyZZ)
    (leftIndex : Nat) (terms : SparsePolyZZ) :
    intTermsToPoly (Generated.StrictRecombine.multiplyTermsLoop
      left right leftIndex terms).toList =
      intTermsToPoly terms.toList +
        intTermsToPoly (left.toList.drop leftIndex) *
          intTermsToPoly right.toList := by
  induction hmeasure : left.size - leftIndex using Nat.strong_induction_on
      generalizing leftIndex terms with
  | h measure ih =>
      rw [Generated.StrictRecombine.multiplyTermsLoop]
      split
      next hleft =>
        rw [ih (left.size - (leftIndex + 1)) (by omega)
          (leftIndex + 1)
          (Generated.StrictRecombine.multiplyRowLoop left[leftIndex]
            right 0 terms) rfl]
        rw [multiplyRowLoop_toPoly]
        have hsuffix : left.toList.drop leftIndex = left[leftIndex] ::
            left.toList.drop (leftIndex + 1) := by
          simpa using List.drop_eq_getElem_cons
            (l := left.toList) (i := leftIndex) (by simpa using hleft)
        simp [intTermsToPoly, hsuffix]
        ring
      next hleft =>
        have hle : left.size ≤ leftIndex := Nat.le_of_not_gt hleft
        simp [intTermsToPoly, List.drop_eq_nil_iff.mpr hle]

private theorem intTermsToPoly_modify_add (terms : List (UMonomial × Int))
    (index : Nat) (term : UMonomial × Int) (hindex : index < terms.length)
    (hdegree : terms[index].1.deg = term.1.deg) :
    intTermsToPoly (terms.modify index
      (fun existing => (existing.1, existing.2 + term.2))) =
      intTermsToPoly terms + Polynomial.monomial term.1.deg term.2 := by
  induction terms generalizing index with
  | nil => simp at hindex
  | cons head tail ih =>
      cases index with
      | zero =>
          simp [intTermsToPoly] at hdegree ⊢
          rw [hdegree]
          abel
      | succ index =>
          have htail : index < tail.length := by simpa using hindex
          have hdegree' : tail[index].1.deg = term.1.deg := by
            simpa using hdegree
          simp only [List.modify_succ_cons, intTermsToPoly, List.map_cons,
            List.sum_cons]
          change Polynomial.monomial head.1.deg head.2 +
              intTermsToPoly (tail.modify index
                (fun existing => (existing.1, existing.2 + term.2))) = _
          rw [ih index htail hdegree']
          change _ + (intTermsToPoly tail + _) =
            _ + intTermsToPoly tail + _
          abel

private theorem groupTermsStep_toPoly (acc : SparsePolyZZ)
    (term : UMonomial × Int) :
    intTermsToPoly
      ((match acc.findIdx? (fun t : UMonomial × Int =>
          t.fst.deg = term.fst.deg) with
        | some index => acc.modify index
            (fun existing => (existing.1, existing.2 + term.2))
        | none => acc.push term).toList) =
      intTermsToPoly acc.toList + Polynomial.monomial term.1.deg term.2 := by
  split
  next index hfind =>
    obtain ⟨hindex, hdegree, _⟩ :=
      Array.findIdx?_eq_some_iff_getElem.mp hfind
    simpa using intTermsToPoly_modify_add acc.toList index term
      (by simpa using hindex) (by simpa using hdegree)
  next hfind =>
    simp [intTermsToPoly]

private theorem groupTerms_toPoly_aux (source : List (UMonomial × Int))
    (acc : SparsePolyZZ) :
    intTermsToPoly
      (source.foldl (fun acc term =>
        match acc.findIdx? (fun t => t.fst.deg = term.fst.deg) with
        | some index => acc.modify index
            (fun existing => (existing.1, existing.2 + term.2))
        | none => acc.push term) acc).toList =
      intTermsToPoly acc.toList + intTermsToPoly source := by
  induction source generalizing acc with
  | nil => simp [intTermsToPoly]
  | cons head tail ih =>
      simp only [List.foldl_cons]
      rw [ih]
      rw [groupTermsStep_toPoly]
      simp [intTermsToPoly]
      abel

private theorem groupTerms_toPoly (terms : SparsePolyZZ) :
    intTermsToPoly
      (terms.foldl (fun acc term =>
        match acc.findIdx? (fun t => t.fst.deg = term.fst.deg) with
        | some index => acc.modify index
            (fun existing => (existing.1, existing.2 + term.2))
        | none => acc.push term) #[]).toList =
      intTermsToPoly terms.toList := by
  rw [← Array.foldl_toList]
  simpa [intTermsToPoly] using groupTerms_toPoly_aux terms.toList #[]

private def sparseZZDegrees (terms : SparsePolyZZ) : List Nat :=
  terms.toList.map (fun term => term.1.deg)

private theorem map_modify_degree_preserved
    (terms : List (UMonomial × Int)) (index : Nat)
    (update : UMonomial × Int → UMonomial × Int)
    (hupdate : ∀ term, (update term).1.deg = term.1.deg) :
    (terms.modify index update).map (fun term => term.1.deg) =
      terms.map (fun term => term.1.deg) := by
  induction terms generalizing index with
  | nil => simp
  | cons head tail ih =>
      cases index with
      | zero => simp [hupdate]
      | succ index => simp [ih index]

private theorem groupTermsStep_degrees_nodup (acc : SparsePolyZZ)
    (term : UMonomial × Int) (hacc : (sparseZZDegrees acc).Nodup) :
    (sparseZZDegrees
      (match acc.findIdx? (fun t => t.1.deg = term.1.deg) with
      | some index => acc.modify index
          (fun existing => (existing.1, existing.2 + term.2))
      | none => acc.push term)).Nodup := by
  split
  next index hfind =>
    unfold sparseZZDegrees at hacc ⊢
    rw [Array.toList_modify]
    rw [map_modify_degree_preserved acc.toList index _ (by simp)]
    exact hacc
  next hfind =>
    have hnone := Array.findIdx?_eq_none_iff.mp hfind
    unfold sparseZZDegrees at hacc ⊢
    simp only [Array.toList_push, List.map_append, List.map_cons, List.map_nil,
      List.nodup_append, List.nodup_singleton]
    refine ⟨hacc, trivial, ?_⟩
    intro degree hdegreeMem singleton hsingleton
    simp only [List.mem_singleton] at hsingleton
    subst singleton
    rw [List.mem_map] at hdegreeMem
    obtain ⟨existing, hexisting, hdegree⟩ := hdegreeMem
    have hnot := hnone existing (by simpa using hexisting)
    simpa [hdegree] using hnot

private theorem groupTermsFold_degrees_nodup
    (source : List (UMonomial × Int)) (acc : SparsePolyZZ)
    (hacc : (sparseZZDegrees acc).Nodup) :
    (sparseZZDegrees (source.foldl (fun acc term =>
      match acc.findIdx? (fun t => t.1.deg = term.1.deg) with
      | some index => acc.modify index
          (fun existing => (existing.1, existing.2 + term.2))
      | none => acc.push term) acc)).Nodup := by
  induction source generalizing acc with
  | nil => simpa
  | cons head tail ih =>
      simp only [List.foldl_cons]
      exact ih _ (groupTermsStep_degrees_nodup acc head hacc)

/-- The concrete grouping fold used by integer sparse normalization retains
exactly one cell for every encountered degree. -/
theorem groupedTerms_degrees_nodup (terms : SparsePolyZZ) :
    (sparseZZDegrees (terms.foldl (fun acc term =>
      match acc.findIdx? (fun t => t.1.deg = term.1.deg) with
      | some index => acc.modify index
          (fun existing => (existing.1, existing.2 + term.2))
      | none => acc.push term) #[])).Nodup := by
  rw [← Array.foldl_toList]
  apply groupTermsFold_degrees_nodup
  simp [sparseZZDegrees]

private def sparseZZDegreeGT (left right : UMonomial × Int) : Prop :=
  left.1.deg > right.1.deg

private theorem pairwise_merge_sparseZZDegreeGT
    (left right : List (UMonomial × Int))
    (hleft : left.Pairwise sparseZZDegreeGT)
    (hright : right.Pairwise sparseZZDegreeGT)
    (hcross : ∀ a ∈ left, ∀ b ∈ right, a.1.deg ≠ b.1.deg) :
    (List.merge left right
      (fun a b => a.1.deg > b.1.deg)).Pairwise sparseZZDegreeGT := by
  induction left generalizing right with
  | nil => simpa only [List.merge]
  | cons x xs ihLeft =>
      induction right with
      | nil => simpa only [List.merge]
      | cons y ys ihRight =>
          simp only [List.merge]
          split <;> rename_i hcompare
          · apply List.Pairwise.cons
            have hxy : x.1.deg > y.1.deg := of_decide_eq_true hcompare
            · intro z hz
              rw [List.mem_merge, List.mem_cons] at hz
              rcases hz with (hz | rfl | hz)
              · exact List.rel_of_pairwise_cons hleft hz
              · exact hxy
              · unfold sparseZZDegreeGT
                exact Nat.lt_trans
                  (List.rel_of_pairwise_cons hright hz) hxy
            · exact ihLeft _ hleft.tail hright (fun a ha b hb =>
                hcross a (List.mem_cons_of_mem x ha) b hb)
          · apply List.Pairwise.cons
            · intro z hz
              rw [List.mem_merge, List.mem_cons] at hz
              simp only [Bool.not_eq_true] at hcompare
              have hyx : y.1.deg > x.1.deg := by
                have hne := hcross x List.mem_cons_self y List.mem_cons_self
                have hnxy : ¬ x.1.deg > y.1.deg :=
                  of_decide_eq_false hcompare
                omega
              rcases hz with (⟨rfl | hz⟩ | hz)
              · exact hyx
              · unfold sparseZZDegreeGT
                exact Nat.lt_trans
                  (List.rel_of_pairwise_cons hleft hz) hyx
              · exact List.rel_of_pairwise_cons hright hz
            · exact ihRight hright.tail (fun a ha b hb =>
                hcross a ha b (List.mem_cons_of_mem y hb))

/-- Strict `deg >` is not a total comparator on duplicate degrees.  The
actual grouping pass removes precisely that obstruction, so the source's
strict-comparison merge sort is nevertheless a strict descending chain. -/
theorem pairwise_mergeSort_sparseZZDegreeGT
    (terms : List (UMonomial × Int))
    (hnodup : (terms.map fun term => term.1.deg).Nodup) :
    (terms.mergeSort
      (fun a b => a.1.deg > b.1.deg)).Pairwise sparseZZDegreeGT := by
  induction hlength : terms.length using Nat.strong_induction_on
      generalizing terms with
  | h length ih =>
    cases terms with
    | nil => simp
    | cons a tail =>
      cases tail with
      | nil => simp
      | cons b xs =>
      rw [List.mergeSort]
      let leftTerms :=
        (List.MergeSort.Internal.splitInTwo ⟨a :: b :: xs, rfl⟩).1.1
      let rightTerms :=
        (List.MergeSort.Internal.splitInTwo ⟨a :: b :: xs, rfl⟩).2.1
      have hleftLength : leftTerms.length < length := by
        rw [← hlength]
        simp [leftTerms, List.MergeSort.Internal.splitInTwo_fst]
        omega
      have hrightLength : rightTerms.length < length := by
        rw [← hlength]
        simp [rightTerms, List.MergeSort.Internal.splitInTwo_snd]
        omega
      have happend : leftTerms ++ rightTerms = a :: b :: xs :=
        List.MergeSort.Internal.splitInTwo_fst_append_splitInTwo_snd _
      have hdegreeAppend :
          ((leftTerms.map fun term => term.1.deg) ++
            (rightTerms.map fun term => term.1.deg)).Nodup := by
        rw [← List.map_append, happend]
        exact hnodup
      have hparts := List.nodup_append.mp hdegreeAppend
      apply pairwise_merge_sparseZZDegreeGT
      · exact ih leftTerms.length hleftLength leftTerms hparts.1 rfl
      · exact ih rightTerms.length hrightLength rightTerms hparts.2.1 rfl
      · intro x hx y hy hxy
        have hx' : x ∈ leftTerms := List.mem_mergeSort.mp hx
        have hy' : y ∈ rightTerms := List.mem_mergeSort.mp hy
        have hxdeg : x.1.deg ∈
            leftTerms.map (fun term => term.1.deg) :=
          List.mem_map.mpr ⟨x, hx', rfl⟩
        have hydeg : y.1.deg ∈
            rightTerms.map (fun term => term.1.deg) :=
          List.mem_map.mpr ⟨y, hy', rfl⟩
        have hxdegRight : x.1.deg ∈
            rightTerms.map (fun term => term.1.deg) := by
          rw [hxy]
          exact hydeg
        exact (hparts.2.2 x.1.deg hxdeg x.1.deg hxdegRight) rfl

private theorem filterNonzero_toPoly (terms : SparsePolyZZ) :
    intTermsToPoly (terms.filter (fun term => term.2 ≠ 0)).toList =
      intTermsToPoly terms.toList := by
  rw [Array.toList_filter]
  induction terms.toList with
  | nil => simp [intTermsToPoly]
  | cons head tail ih =>
      rcases head with ⟨monomial, coefficient⟩
      by_cases hzero : coefficient = 0
      · simpa [hzero, intTermsToPoly] using ih
      · simpa [hzero, intTermsToPoly] using ih

private theorem mergeSort_toPoly (terms : List (UMonomial × Int)) :
    intTermsToPoly
        (terms.mergeSort (fun a b => a.fst.deg > b.fst.deg)) =
      intTermsToPoly terms := by
  exact ((List.mergeSort_perm terms _).map
    (fun term => Polynomial.monomial term.1.deg term.2)).sum_eq

theorem normalization_toPoly (terms : SparsePolyZZ) :
    SparsePolyZZ.toPoly (SparsePolyZZ.normalization terms) =
      SparsePolyZZ.toPoly terms := by
  unfold SparsePolyZZ.normalization SparsePolyZZ.toPoly
  rw [List.toList_toArray]
  change intTermsToPoly
      (((terms.foldl (fun acc term =>
          match acc.findIdx? (fun t => t.fst.deg = term.fst.deg) with
          | some index => acc.modify index
              (fun existing => (existing.1, existing.2 + term.2))
          | none => acc.push term) #[]).filter
        (fun term => term.2 ≠ 0)).toList.mergeSort
          (fun a b => a.fst.deg > b.fst.deg)) = intTermsToPoly terms.toList
  rw [mergeSort_toPoly, filterNonzero_toPoly, groupTerms_toPoly]

/-- Every coefficient retained by the concrete sparse normalization is
nonzero.  This is the exact property needed by the checked long-division
loop; degree descent is supplied separately by its runtime rank check. -/
theorem normalization_coefficients_nonzero (terms : SparsePolyZZ) :
    ∀ term ∈ (SparsePolyZZ.normalization terms).toList, term.2 ≠ 0 := by
  intro term hterm
  unfold SparsePolyZZ.normalization at hterm
  simp only [List.toList_toArray] at hterm
  have hfiltered : term ∈
      ((terms.foldl (fun acc term =>
        match acc.findIdx? (fun t => t.fst.deg = term.fst.deg) with
        | some index => acc.modify index
            (fun existing => (existing.1, existing.2 + term.2))
        | none => acc.push term) #[]).filter
          (fun term => term.2 ≠ 0)).toList := by
    exact (List.Perm.mem_iff (List.mergeSort_perm _ _)).mp hterm
  rw [Array.toList_filter, List.mem_filter] at hfiltered
  simpa using hfiltered.2

/-- The exact group/filter/strict-merge-sort normalization used by C++ always
produces the canonical sparse integer representation.  In particular, this
does not replace the generated sort by a mathematical sorting oracle. -/
theorem normalization_canonical (terms : SparsePolyZZ) :
    StrictPolynomialMod.SparsePolyZZCanonical
      (SparsePolyZZ.normalization terms) := by
  let grouped : SparsePolyZZ := terms.foldl (fun acc term =>
    match acc.findIdx? (fun t => t.1.deg = term.1.deg) with
    | some index => acc.modify index
        (fun existing => (existing.1, existing.2 + term.2))
    | none => acc.push term) #[]
  let nonzero : SparsePolyZZ := grouped.filter (fun term => term.2 ≠ 0)
  have hgrouped : (sparseZZDegrees grouped).Nodup := by
    exact groupedTerms_degrees_nodup terms
  have hsub : List.Sublist nonzero.toList grouped.toList := by
    dsimp [nonzero]
    rw [Array.toList_filter]
    exact List.filter_sublist
  have hnonzeroDegrees :
      (nonzero.toList.map (fun term => term.1.deg)).Nodup := by
    exact hgrouped.sublist (hsub.map (fun term => term.1.deg))
  unfold StrictPolynomialMod.SparsePolyZZCanonical
  constructor
  · unfold SparsePolyZZ.normalization
    change List.IsChain (fun left right : UMonomial × Int =>
        left.1.deg > right.1.deg)
      (nonzero.toList.mergeSort (fun left right =>
        left.1.deg > right.1.deg))
    exact (pairwise_mergeSort_sparseZZDegreeGT nonzero.toList
      hnonzeroDegrees).isChain
  · exact normalization_coefficients_nonzero terms

theorem multiplyNormalizeRaw_toPoly (left right output : SparsePolyZZ)
    (hrun : Generated.StrictRecombine.multiplyNormalizeRaw left right =
      .ok output) :
    SparsePolyZZ.toPoly output =
      SparsePolyZZ.toPoly left * SparsePolyZZ.toPoly right := by
  unfold Generated.StrictRecombine.multiplyNormalizeRaw at hrun
  have houtput := Except.ok.inj hrun
  subst output
  rw [normalization_toPoly]
  change intTermsToPoly
      (Generated.StrictRecombine.multiplyTermsLoop left right 0 #[]).toList = _
  rw [multiplyTermsLoop_toPoly]
  simp [intTermsToPoly, SparsePolyZZ.toPoly]

theorem multiplyNormalizeRaw_canonical (left right output : SparsePolyZZ)
    (hrun : Generated.StrictRecombine.multiplyNormalizeRaw left right =
      .ok output) :
    StrictPolynomialMod.SparsePolyZZCanonical output := by
  unfold Generated.StrictRecombine.multiplyNormalizeRaw at hrun
  have hout := Except.ok.inj hrun
  subst output
  exact normalization_canonical _

private theorem emod_cast (modulus : Nat) (hmodulus : 0 < modulus)
    (coefficient : Int) :
    ((coefficient % (modulus : Int) : Int) : ZMod modulus) =
      (coefficient : ZMod modulus) := by
  rw [ZMod.intCast_eq_intCast_iff]
  refine Int.modEq_iff_dvd.mpr ?_
  use coefficient / (modulus : Int)
  have h := Int.mul_ediv_add_emod coefficient (modulus : Int)
  omega

private theorem emod_cast_of_dvd (base modulus : Nat) (hbase : 0 < base)
    (hdivides : base ∣ modulus) (coefficient : Int) :
    ((coefficient % (modulus : Int) : Int) : ZMod base) =
      (coefficient : ZMod base) := by
  rw [ZMod.intCast_eq_intCast_iff]
  rw [Int.modEq_iff_dvd]
  apply dvd_trans (b := (modulus : Int))
  · exact_mod_cast hdivides
  · use coefficient / (modulus : Int)
    have h := Int.mul_ediv_add_emod coefficient (modulus : Int)
    omega

/-- Exact coefficient cells emitted by the generated modular-reduction loop,
in source order. -/
theorem modCoeffLoop_toList (input : SparsePolyZZ) (modulus : ZZ)
    (index : Nat) (result output : SparsePolyZZ)
    (hrun : Generated.StrictRecombine.modCoeffLoop input modulus index result =
      .ok output) :
    output.toList = result.toList ++
      (input.toList.drop index).filterMap (fun term =>
        let coefficient := term.2 % modulus
        if coefficient = 0 then none else some (term.1, coefficient)) := by
  induction hmeasure : input.size - index using Nat.strong_induction_on
      generalizing index result output with
  | h measure ih =>
      rw [Generated.StrictRecombine.modCoeffLoop] at hrun
      split at hrun
      next hindex =>
        split at hrun
        next hmodulus =>
          let coefficient := input[index].2 % modulus
          by_cases hzero : coefficient = 0
          · change input[index].2 % modulus = 0 at hzero
            simp only [hzero, if_true] at hrun
            rw [ih (input.size - (index + 1)) (by omega)
              (index + 1) result output hrun rfl]
            have hsuffix : input.toList.drop index = input[index] ::
                input.toList.drop (index + 1) := by
              simpa using List.drop_eq_getElem_cons
                (l := input.toList) (i := index) (by simpa using hindex)
            rw [hsuffix, List.filterMap_cons]
            dsimp only
            rw [if_pos hzero]
          · change input[index].2 % modulus ≠ 0 at hzero
            simp only [hzero, if_false] at hrun
            rw [ih (input.size - (index + 1)) (by omega)
              (index + 1)
              (result.push (input[index].1, input[index].2 % modulus))
              output hrun rfl]
            have hsuffix : input.toList.drop index = input[index] ::
                input.toList.drop (index + 1) := by
              simpa using List.drop_eq_getElem_cons
                (l := input.toList) (i := index) (by simpa using hindex)
            rw [hsuffix, List.filterMap_cons]
            dsimp only
            rw [if_neg hzero, Array.toList_push, List.append_assoc]
            rfl
        next hmodulus => contradiction
      next hindex =>
        have hle : input.size ≤ index := Nat.le_of_not_gt hindex
        have hout := Except.ok.inj hrun
        subst output
        simp [List.drop_eq_nil_iff.mpr hle]

/-- Starting from the source's empty accumulator, modular coefficient
reduction preserves strict sparse order and omits every zero residue. -/
theorem modCoeffLoop_canonical (input output : SparsePolyZZ)
    (modulus : ZZ)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical input)
    (hrun : Generated.StrictRecombine.modCoeffLoop input modulus 0 #[] =
      .ok output) :
    StrictPolynomialMod.SparsePolyZZCanonical output := by
  let reduceTerm : UMonomial × Int → Option (UMonomial × Int) :=
    fun term =>
      let coefficient := term.2 % modulus
      if coefficient = 0 then none else some (term.1, coefficient)
  have hlist := modCoeffLoop_toList input modulus 0 #[] output hrun
  have houtput : output.toList = input.toList.filterMap reduceTerm := by
    simpa [reduceTerm] using hlist
  unfold StrictPolynomialMod.SparsePolyZZCanonical
  rw [houtput]
  constructor
  · letI : Trans
        (fun left right : UMonomial × Int => left.1.deg > right.1.deg)
        (fun left right : UMonomial × Int => left.1.deg > right.1.deg)
        (fun left right : UMonomial × Int => left.1.deg > right.1.deg) :=
      ⟨by intro left middle right hleft hright
          exact Nat.lt_trans hright hleft⟩
    exact (hcanonical.1.pairwise.filterMap reduceTerm (by
      intro left right hdegree left' hleft right' hright
      dsimp [reduceTerm] at hleft hright
      split at hleft <;> try contradiction
      split at hright <;> try contradiction
      injection hleft with hleft
      injection hright with hright
      subst left'
      subst right'
      exact hdegree)).isChain
  · intro term hterm
    rw [List.mem_filterMap] at hterm
    obtain ⟨source, hsource, hreduce⟩ := hterm
    dsimp [reduceTerm] at hreduce
    split at hreduce
    next hzero => contradiction
    next hzero =>
      injection hreduce with heq
      subst term
      exact hzero

theorem modCoeffLoop_toPolyMod (input : SparsePolyZZ)
    (modulus : Nat) (hmodulus : 0 < modulus) (index : Nat)
    (result output : SparsePolyZZ)
    (hrun : Generated.StrictRecombine.modCoeffLoop input (modulus : ZZ)
      index result = .ok output) :
    Refinement.StrictHensel.toPolyMod modulus output =
      Refinement.StrictHensel.toPolyMod modulus result +
        Refinement.StrictHensel.termsToPolyMod modulus
          (input.toList.drop index) := by
  induction hmeasure : input.size - index using Nat.strong_induction_on
      generalizing index result output with
  | h measure ih =>
      rw [Generated.StrictRecombine.modCoeffLoop] at hrun
      split at hrun
      next hindex =>
        split at hrun
        next hmodulus' =>
          let coefficient := input[index].2 % (modulus : Int)
          by_cases hzero : coefficient = 0
          · change input[index].2 % (modulus : Int) = 0 at hzero
            simp only [hzero, if_true] at hrun
            rw [ih (input.size - (index + 1)) (by omega)
              (index + 1) result output hrun rfl]
            have hsuffix : input.toList.drop index = input[index] ::
                input.toList.drop (index + 1) := by
              simpa using List.drop_eq_getElem_cons
                (l := input.toList) (i := index) (by simpa using hindex)
            have hcast := emod_cast modulus hmodulus input[index].2
            change (coefficient : ZMod modulus) =
              (input[index].2 : ZMod modulus) at hcast
            dsimp [coefficient] at hcast
            rw [hzero] at hcast
            rw [hsuffix, Refinement.StrictHensel.termsToPolyMod_cons]
            rw [← hcast]
            simp
          · change input[index].2 % (modulus : Int) ≠ 0 at hzero
            simp only [hzero, if_false] at hrun
            rw [ih (input.size - (index + 1)) (by omega)
              (index + 1) (result.push (input[index].1, coefficient))
              output hrun rfl]
            rw [Refinement.StrictHensel.toPolyMod_push]
            have hsuffix : input.toList.drop index = input[index] ::
                input.toList.drop (index + 1) := by
              simpa using List.drop_eq_getElem_cons
                (l := input.toList) (i := index) (by simpa using hindex)
            rw [emod_cast modulus hmodulus input[index].2]
            rw [hsuffix, Refinement.StrictHensel.termsToPolyMod_cons]
            abel
        next hmodulus' => contradiction
      next hindex =>
        have hle : input.size ≤ index := Nat.le_of_not_gt hindex
        have hout := Except.ok.inj hrun
        subst output
        simp [Refinement.StrictHensel.termsToPolyMod,
          List.drop_eq_nil_iff.mpr hle]

/-- The actual `% modulus` loop remains the identity after reduction at every
positive divisor `base` of that concrete modulus. -/
theorem modCoeffLoop_toPolyMod_of_dvd (input : SparsePolyZZ)
    (modulus base : Nat) (hmodulus : 0 < modulus) (hbase : 0 < base)
    (hdivides : base ∣ modulus)
    (index : Nat) (result output : SparsePolyZZ)
    (hrun : Generated.StrictRecombine.modCoeffLoop input (modulus : ZZ)
      index result = .ok output) :
    Refinement.StrictHensel.toPolyMod base output =
      Refinement.StrictHensel.toPolyMod base result +
        Refinement.StrictHensel.termsToPolyMod base
          (input.toList.drop index) := by
  induction hmeasure : input.size - index using Nat.strong_induction_on
      generalizing index result output with
  | h measure ih =>
      rw [Generated.StrictRecombine.modCoeffLoop] at hrun
      split at hrun
      next hindex =>
        have hmodulusInt : (modulus : Int) ≠ 0 := by
          exact_mod_cast (Nat.ne_of_gt hmodulus)
        rw [dif_pos hmodulusInt] at hrun
        let coefficient := input[index].2 % (modulus : Int)
        by_cases hzero : coefficient = 0
        · change input[index].2 % (modulus : Int) = 0 at hzero
          simp only [hzero, if_true] at hrun
          rw [ih (input.size - (index + 1)) (by omega)
            (index + 1) result output hrun rfl]
          have hsuffix : input.toList.drop index = input[index] ::
              input.toList.drop (index + 1) := by
            simpa using List.drop_eq_getElem_cons
              (l := input.toList) (i := index) (by simpa using hindex)
          have hcast := emod_cast_of_dvd base modulus hbase hdivides
            input[index].2
          change (coefficient : ZMod base) =
            (input[index].2 : ZMod base) at hcast
          dsimp [coefficient] at hcast
          rw [hzero] at hcast
          rw [hsuffix, Refinement.StrictHensel.termsToPolyMod_cons]
          rw [← hcast]
          simp
        · change input[index].2 % (modulus : Int) ≠ 0 at hzero
          simp only [hzero, if_false] at hrun
          rw [ih (input.size - (index + 1)) (by omega)
            (index + 1) (result.push (input[index].1, coefficient))
            output hrun rfl]
          rw [Refinement.StrictHensel.toPolyMod_push]
          have hsuffix : input.toList.drop index = input[index] ::
              input.toList.drop (index + 1) := by
            simpa using List.drop_eq_getElem_cons
              (l := input.toList) (i := index) (by simpa using hindex)
          rw [emod_cast_of_dvd base modulus hbase hdivides input[index].2]
          rw [hsuffix, Refinement.StrictHensel.termsToPolyMod_cons]
          abel
      next hindex =>
        have hle : input.size ≤ index := Nat.le_of_not_gt hindex
        have hout := Except.ok.inj hrun
        subst output
        simp [Refinement.StrictHensel.termsToPolyMod,
          List.drop_eq_nil_iff.mpr hle]

theorem multiplyNormalizeModRaw_correct
    (left right : SparsePolyZZ) (modulus : Nat) (output : SparsePolyZZ)
    (hmodulus : 0 < modulus)
    (hrun : Generated.StrictRecombine.multiplyNormalizeModRaw left right
      (modulus : ZZ) = .ok output) :
    Refinement.StrictHensel.toPolyMod modulus output =
      Refinement.StrictHensel.toPolyMod modulus left *
        Refinement.StrictHensel.toPolyMod modulus right := by
  unfold Generated.StrictRecombine.multiplyNormalizeModRaw at hrun
  split at hrun
  next fault hmultiply => contradiction
  next product hmultiply =>
    have hmod := modCoeffLoop_toPolyMod product modulus hmodulus 0 #[]
      output hrun
    have hmul := multiplyNormalizeRaw_toPoly left right product hmultiply
    have hmod' : Refinement.StrictHensel.toPolyMod modulus output =
        Refinement.StrictHensel.toPolyMod modulus product := by
      simpa [Refinement.StrictHensel.toPolyMod_eq_termsToPolyMod] using hmod
    rw [hmod']
    simpa [Refinement.StrictHensel.toPolyMod] using
      congrArg (Polynomial.map (Int.castRingHom (ZMod modulus))) hmul

theorem multiplyNormalizeModRaw_correct_of_dvd
    (left right : SparsePolyZZ) (modulus base : Nat)
    (output : SparsePolyZZ) (hmodulus : 0 < modulus) (hbase : 0 < base)
    (hdivides : base ∣ modulus)
    (hrun : Generated.StrictRecombine.multiplyNormalizeModRaw left right
      (modulus : ZZ) = .ok output) :
    Refinement.StrictHensel.toPolyMod base output =
      Refinement.StrictHensel.toPolyMod base left *
        Refinement.StrictHensel.toPolyMod base right := by
  unfold Generated.StrictRecombine.multiplyNormalizeModRaw at hrun
  split at hrun
  next fault hmultiply => contradiction
  next product hmultiply =>
    have hmod := modCoeffLoop_toPolyMod_of_dvd product modulus base hmodulus
      hbase hdivides 0 #[] output hrun
    have hmul := multiplyNormalizeRaw_toPoly left right product hmultiply
    have hmod' : Refinement.StrictHensel.toPolyMod base output =
        Refinement.StrictHensel.toPolyMod base product := by
      simpa [Refinement.StrictHensel.toPolyMod_eq_termsToPolyMod] using hmod
    rw [hmod']
    simpa [Refinement.StrictHensel.toPolyMod] using
      congrArg (Polynomial.map (Int.castRingHom (ZMod base))) hmul

theorem multiplyNormalizeModRaw_canonical
    (left right output : SparsePolyZZ) (modulus : ZZ)
    (hrun : Generated.StrictRecombine.multiplyNormalizeModRaw left right
      modulus = .ok output) :
    StrictPolynomialMod.SparsePolyZZCanonical output := by
  unfold Generated.StrictRecombine.multiplyNormalizeModRaw at hrun
  split at hrun
  next fault hmultiply => contradiction
  next product hmultiply =>
    exact modCoeffLoop_canonical product output modulus
      (multiplyNormalizeRaw_canonical left right product hmultiply) hrun

/-- The literal coefficient-reduction loop cannot fail at a nonzero modulus.
This is an execution theorem for the generated loop, not a semantic
replacement for its result. -/
theorem modCoeffLoop_complete (input : SparsePolyZZ) (modulus : ZZ)
    (index : Nat) (result : SparsePolyZZ) (hmodulus : modulus ≠ 0) :
    ∃ output, Generated.StrictRecombine.modCoeffLoop input modulus index result =
      .ok output := by
  induction hmeasure : input.size - index using Nat.strong_induction_on
      generalizing index result with
  | h measure ih =>
      rw [Generated.StrictRecombine.modCoeffLoop]
      split
      next hindex =>
        exact ih (input.size - (index + 1)) (by omega) (index + 1)
          (if input[index].2 % modulus = 0 then result
           else result.push (input[index].1, input[index].2 % modulus)) rfl
      next hindex => exact ⟨result, rfl⟩

/-- Multiplication followed by the generated modular normalization always
returns a physical sparse array when the modulus is nonzero. -/
theorem multiplyNormalizeModRaw_complete (left right : SparsePolyZZ)
    (modulus : ZZ) (hmodulus : modulus ≠ 0) :
    ∃ output,
      Generated.StrictRecombine.multiplyNormalizeModRaw left right modulus =
        .ok output := by
  unfold Generated.StrictRecombine.multiplyNormalizeModRaw
  simp only [Generated.StrictRecombine.multiplyNormalizeRaw]
  exact modCoeffLoop_complete
    (SparsePolyZZ.normalization
      (Generated.StrictRecombine.multiplyTermsLoop left right 0 #[]))
    modulus 0 #[] hmodulus

/-- Every valid candidate executes the complete generated trial-product loop.
The witness is the array computed by the loop itself; no product oracle is
present in `TrialProductRawOps`. -/
theorem trialProductLoop_complete
    (ops : Generated.StrictRecombine.TrialProductRawOps)
    (candidate : Array Int32) (activeLifted : Array SparsePolyZZ)
    (modulus : ZZ) (index : Nat) (product : SparsePolyZZ)
    (hmodulus : modulus ≠ 0)
    (hvalid : CandidateIndicesValid candidate
      (Array.replicate activeLifted.size false)) :
    ∃ output, Generated.StrictRecombine.trialProductLoop ops candidate
      activeLifted modulus index product = .ok output := by
  induction hmeasure : candidate.size - index using Nat.strong_induction_on
      generalizing index product with
  | h measure ih =>
      rw [Generated.StrictRecombine.trialProductLoop]
      split
      next hindex =>
        have hentry := hvalid index hindex
        have hactive : candidate[index].toInt64.toNat < activeLifted.size := by
          simpa using hentry.2
        rw [dif_pos hentry.1, dif_pos hactive]
        rcases multiplyNormalizeModRaw_complete product
            (activeLifted[candidate[index].toInt64.toNat]'hactive) modulus
            hmodulus with
          ⟨product', hproduct'⟩
        rw [hproduct']
        exact ih (candidate.size - (index + 1)) (by omega) (index + 1)
          product' rfl
      next hindex => exact ⟨product, rfl⟩

/-- The generated symmetric-coefficient traversal is total. -/
theorem symmetricModLoop_complete (input : SparsePolyZZ) (modulus : ZZ)
    (index : Nat) (result : SparsePolyZZ) :
    ∃ output, Generated.StrictRecombine.symmetricModLoop input modulus index
      result = .ok output := by
  induction hmeasure : input.size - index using Nat.strong_induction_on
      generalizing index result with
  | h measure ih =>
      rw [Generated.StrictRecombine.symmetricModLoop]
      split
      next hindex =>
        exact ih (input.size - (index + 1)) (by omega) (index + 1)
          (if ZZ.symmetricMod input[index].2 modulus = 0 then result
           else result.push
             (input[index].1, ZZ.symmetricMod input[index].2 modulus)) rfl
      next hindex => exact ⟨result, rfl⟩

/-- A positive modulus forces the literal generated symmetric-mod entry down
its successful branch. -/
theorem symmetricModRaw_complete (input : SparsePolyZZ) (modulus : ZZ)
    (hmodulus : 0 < modulus) :
    ∃ output,
      Generated.StrictRecombine.symmetricModRaw input modulus = .ok output := by
  unfold Generated.StrictRecombine.symmetricModRaw
  rw [dif_pos hmodulus]
  exact symmetricModLoop_complete input modulus 0 #[]

/-- Concrete candidate-validation dependencies.  Both generated operation
records are data-free, so this value cannot choose candidates, products, or
factorization witnesses: validation executes only the generated raw loops. -/
def concreteCandidateValidationRawOps :
    Generated.StrictRecombine.CandidateValidationRawOps where
  product := {}

/-- The sole concrete operation bundle for the generated van-Hoeij loop.
Every field is discharged by the generated raw execution proved above: CLD
uses the concrete well-founded modular division, LLL uses its determinant
rank, lattice extension/reset use the actual matrix constructors, candidate
validation executes the generated loops, and fallback uses the concrete
combination rank. -/
noncomputable def concreteVanHoeijRawOps :
    Generated.StrictRecombine.VanHoeijRawOps where
  cld := { divmodTermination := StrictHensel.concreteDivmodTermination }
  lll := concreteLLLExecution
  gather_size := by
    intro active lifted activeLifted hrun
    exact gatherActive_size_of_success active lifted activeLifted hrun
  cld_size := by
    intro fStar activeFactors modulus cldOutput hrun
    exact cldPolys_size_of_success
      { divmodTermination := StrictHensel.concreteDivmodTermination }
      fStar activeFactors cldOutput modulus hrun
  cld_extension_valid := by
    intro matrix cld current target matrix' added hinput hdimension hrun
    exact buildCldMatrix_input_valid matrix matrix' cld current target added
      hinput hdimension hrun
  reset_valid := by
    intro factorCount matrix bound hrun
    exact resetVanHoeijLattice_input_valid factorCount matrix bound hrun
  validation := concreteCandidateValidationRawOps
  zassenhausTermination := concreteZassenhausTermination

/-- Genuine active-set termination for the sole concrete van-Hoeij bundle. -/
noncomputable def concreteVanHoeijTermination :
    Generated.StrictRecombine.VanHoeijTermination concreteVanHoeijRawOps :=
  removalTermination concreteVanHoeijRawOps

theorem trialProductLoop_refines
    (ops : Generated.StrictRecombine.TrialProductRawOps)
    (candidate : Array Int32) (activeLifted : Array SparsePolyZZ)
    (modulus : Nat) (hmodulus : 0 < modulus) (index : Nat)
    (product output : SparsePolyZZ)
    (hvalid : CandidateIndicesValid candidate
      (Array.replicate activeLifted.size false))
    (hrun : Generated.StrictRecombine.trialProductLoop ops candidate
      activeLifted (modulus : ZZ) index product = .ok output) :
    Refinement.StrictHensel.toPolyMod modulus output =
      Refinement.StrictHensel.toPolyMod modulus product *
        SelectedProductMod modulus candidate activeLifted index := by
  induction hmeasure : candidate.size - index using Nat.strong_induction_on
      generalizing index product output with
  | h measure ih =>
      rw [Generated.StrictRecombine.trialProductLoop] at hrun
      split at hrun
      next hindex =>
        have hentry := hvalid index hindex
        have hactive : candidate[index].toInt64.toNat < activeLifted.size := by
          simpa using hentry.2
        rw [dif_pos hentry.1, dif_pos hactive] at hrun
        split at hrun
        next fault hcall => contradiction
        next product' hcall =>
          let selected := activeLifted[candidate[index].toInt64.toNat]
          have htail := ih (candidate.size - (index + 1)) (by omega)
            (index + 1) product' output hrun rfl
          have hsuffix : candidate.toList.drop index = candidate[index] ::
              candidate.toList.drop (index + 1) := by
            simpa using List.drop_eq_getElem_cons
              (l := candidate.toList) (i := index) (by simpa using hindex)
          have hstep' := multiplyNormalizeModRaw_correct product selected
            modulus product' hmodulus hcall
          rw [htail, hstep']
          have hselected : selected =
              activeLifted[candidate[index].toInt64.toNat]! := by
            simp [selected, hactive]
          rw [hselected]
          simp [SelectedProductMod, hsuffix, mul_assoc]

      next hindex =>
        have hle : candidate.size ≤ index := Nat.le_of_not_gt hindex
        have hout := Except.ok.inj hrun
        subst output
        simp [SelectedProductMod, List.drop_eq_nil_iff.mpr hle]

/-- The concrete trial product is unchanged when observed modulo any positive
divisor `base` of the actual Hensel modulus.  This follows the generated loop
and its actual `% modulus` multiplications rather than replacing the trial by
an abstract product. -/
theorem trialProductLoop_refines_of_dvd
    (ops : Generated.StrictRecombine.TrialProductRawOps)
    (candidate : Array Int32) (activeLifted : Array SparsePolyZZ)
    (modulus base : Nat) (hmodulus : 0 < modulus) (hbase : 0 < base)
    (hdivides : base ∣ modulus) (index : Nat)
    (product output : SparsePolyZZ)
    (hvalid : CandidateIndicesValid candidate
      (Array.replicate activeLifted.size false))
    (hrun : Generated.StrictRecombine.trialProductLoop ops candidate
      activeLifted (modulus : ZZ) index product = .ok output) :
    Refinement.StrictHensel.toPolyMod base output =
      Refinement.StrictHensel.toPolyMod base product *
        SelectedProductMod base candidate activeLifted index := by
  induction hmeasure : candidate.size - index using Nat.strong_induction_on
      generalizing index product output with
  | h measure ih =>
      rw [Generated.StrictRecombine.trialProductLoop] at hrun
      split at hrun
      next hindex =>
        have hentry := hvalid index hindex
        have hactive : candidate[index].toInt64.toNat < activeLifted.size := by
          simpa using hentry.2
        rw [dif_pos hentry.1, dif_pos hactive] at hrun
        split at hrun
        next fault hcall => contradiction
        next product' hcall =>
          let selected := activeLifted[candidate[index].toInt64.toNat]
          have htail := ih (candidate.size - (index + 1)) (by omega)
            (index + 1) product' output hrun rfl
          have hsuffix : candidate.toList.drop index = candidate[index] ::
              candidate.toList.drop (index + 1) := by
            simpa using List.drop_eq_getElem_cons
              (l := candidate.toList) (i := index) (by simpa using hindex)
          have hstep' := multiplyNormalizeModRaw_correct_of_dvd product
            selected modulus base product' hmodulus hbase hdivides hcall
          rw [htail, hstep']
          have hselected : selected =
              activeLifted[candidate[index].toInt64.toNat]! := by
            simp [selected, hactive]
          rw [hselected]
          simp [SelectedProductMod, hsuffix, mul_assoc]
      next hindex =>
        have hle : candidate.size ≤ index := Nat.le_of_not_gt hindex
        have hout := Except.ok.inj hrun
        subst output
        simp [SelectedProductMod, List.drop_eq_nil_iff.mpr hle]

/-- Canonicality follows the exact successful multiplication/reduction trace;
no validity oracle is needed because every checked source branch is retained. -/
theorem trialProductLoop_canonical
    (ops : Generated.StrictRecombine.TrialProductRawOps)
    (candidate : Array Int32) (activeLifted : Array SparsePolyZZ)
    (modulus : ZZ) (index : Nat) (product output : SparsePolyZZ)
    (hproduct : StrictPolynomialMod.SparsePolyZZCanonical product)
    (hrun : Generated.StrictRecombine.trialProductLoop ops candidate
      activeLifted modulus index product = .ok output) :
    StrictPolynomialMod.SparsePolyZZCanonical output := by
  induction hmeasure : candidate.size - index using Nat.strong_induction_on
      generalizing index product output with
  | h measure ih =>
      rw [Generated.StrictRecombine.trialProductLoop] at hrun
      split at hrun
      next hindex =>
        dsimp at hrun
        split at hrun
        next hnonnegative =>
          split at hrun
          next hactive =>
            split at hrun
            next fault hmultiply => contradiction
            next product' hmultiply =>
              exact ih (candidate.size - (index + 1)) (by omega)
                (index + 1) product' output
                (multiplyNormalizeModRaw_canonical product
                  activeLifted[candidate[index].toInt64.toNat]
                  product' modulus hmultiply) hrun rfl
          next hactive => contradiction
        next hnonnegative => contradiction
      next hindex =>
        have hout := Except.ok.inj hrun
        subst output
        exact hproduct

/-- End-to-end product statement for the candidate path used by
`zassenhausAttempt`: checked lowering followed by the generated multiplication
loop computes the product of the exact occurrence-sensitive source sublist. -/
theorem trialProductLoop_source_indices_refines
    (modulus : Nat) (hmodulus : 0 < modulus)
    (indices : Array Nat) (activeLifted : Array SparsePolyZZ)
    (candidate : Array Int32) (product output : SparsePolyZZ)
    (hbound : ∀ position (hposition : position < indices.size),
      indices[position] < activeLifted.size)
    (hactiveFits : activeLifted.size ≤ 2 ^ 31)
    (hconvert : Generated.StrictRecombine.combinationToInt32 indices =
      .ok candidate)
    (hrun : Generated.StrictRecombine.trialProductLoop ⟨()⟩ candidate
      activeLifted (modulus : ZZ) 0 product = .ok output) :
    Refinement.StrictHensel.toPolyMod modulus output =
      Refinement.StrictHensel.toPolyMod modulus product *
        ((selectSourceIndices activeLifted.toList indices.toList).map
          (Refinement.StrictHensel.toPolyMod modulus)).prod := by
  have hvalid := combinationToInt32_candidate_valid indices activeLifted.size
    candidate hbound hactiveFits hconvert
  have hrefines := trialProductLoop_refines ⟨()⟩ candidate activeLifted modulus
    hmodulus 0 product output hvalid hrun
  rw [selectedProductMod_combinationToInt32 modulus indices activeLifted
    candidate hbound hactiveFits hconvert] at hrefines
  exact hrefines

/-- Source-index form of `trialProductLoop_refines_of_dvd`, used at the exact
`combinationToInt32` boundary executed by `zassenhausAttempt`. -/
theorem trialProductLoop_source_indices_refines_of_dvd
    (modulus base : Nat) (hmodulus : 0 < modulus) (hbase : 0 < base)
    (hdivides : base ∣ modulus)
    (indices : Array Nat) (activeLifted : Array SparsePolyZZ)
    (candidate : Array Int32) (product output : SparsePolyZZ)
    (hbound : ∀ position (hposition : position < indices.size),
      indices[position] < activeLifted.size)
    (hactiveFits : activeLifted.size ≤ 2 ^ 31)
    (hconvert : Generated.StrictRecombine.combinationToInt32 indices =
      .ok candidate)
    (hrun : Generated.StrictRecombine.trialProductLoop ⟨()⟩ candidate
      activeLifted (modulus : ZZ) 0 product = .ok output) :
    Refinement.StrictHensel.toPolyMod base output =
      Refinement.StrictHensel.toPolyMod base product *
        ((selectSourceIndices activeLifted.toList indices.toList).map
          (Refinement.StrictHensel.toPolyMod base)).prod := by
  have hvalid := combinationToInt32_candidate_valid indices activeLifted.size
    candidate hbound hactiveFits hconvert
  have hrefines := trialProductLoop_refines_of_dvd ⟨()⟩ candidate
    activeLifted modulus base hmodulus hbase hdivides 0 product output hvalid
    hrun
  rw [selectedProductMod_combinationToInt32 base indices activeLifted
    candidate hbound hactiveFits hconvert] at hrefines
  exact hrefines

private theorem symmetricMod_cast (modulus : Nat) (hmodulus : 0 < modulus)
    (coefficient : Int) :
    (ZZ.symmetricMod coefficient (modulus : ZZ) : ZMod modulus) =
      (coefficient : ZMod modulus) := by
  unfold ZZ.symmetricMod
  have hfmod : ((Int.fmod coefficient (modulus : Int) : Int) : ZMod modulus) =
      (coefficient : ZMod modulus) := by
    rw [Int.fmod_eq_emod_of_nonneg coefficient (by omega)]
    rw [ZMod.intCast_eq_intCast_iff]
    refine Int.modEq_iff_dvd.mpr ?_
    use coefficient / (modulus : Int)
    have h := Int.mul_ediv_add_emod coefficient (modulus : Int)
    omega
  dsimp [ZZ.symmetricMod]
  split
  next hsmall => exact hfmod
  next hlarge =>
    rw [Int.cast_sub, hfmod]
    simp

private theorem symmetricMod_cast_of_dvd (base modulus : Nat)
    (hbase : 0 < base) (hmodulus : 0 < modulus) (hdivides : base ∣ modulus)
    (coefficient : Int) :
    (ZZ.symmetricMod coefficient (modulus : ZZ) : ZMod base) =
      (coefficient : ZMod base) := by
  unfold ZZ.symmetricMod
  have hfmod : ((Int.fmod coefficient (modulus : Int) : Int) : ZMod base) =
      (coefficient : ZMod base) := by
    rw [Int.fmod_eq_emod_of_nonneg coefficient (by omega)]
    exact emod_cast_of_dvd base modulus hbase hdivides coefficient
  dsimp [ZZ.symmetricMod]
  split
  next hsmall => exact hfmod
  next hlarge =>
    rw [Int.cast_sub, hfmod]
    have hmodulusZero : ((modulus : Int) : ZMod base) = 0 := by
      rw [← Int.cast_zero, ZMod.intCast_eq_intCast_iff,
        Int.modEq_iff_dvd]
      rw [zero_sub]
      apply dvd_neg.mpr
      exact_mod_cast hdivides
    rw [hmodulusZero, sub_zero]

/-- The exact representative selected by `ZZ.symmetricMod` lies in the closed
half-modulus interval.  The inequality is intentionally non-strict: for even
moduli the source keeps the positive midpoint. -/
private theorem symmetricMod_natAbs_mul_two_le (coefficient : Int)
    (modulus : Nat) (hmodulus : 0 < modulus) :
    (ZZ.symmetricMod coefficient (modulus : Int)).natAbs * 2 ≤ modulus := by
  unfold ZZ.symmetricMod
  rw [Int.fmod_eq_emod_of_nonneg coefficient (by omega)]
  let residue := coefficient % (modulus : Int)
  have hresidueNonnegative : 0 ≤ residue := by
    exact Int.emod_nonneg coefficient (by omega)
  have hresidueLt : residue < (modulus : Int) := by
    exact Int.emod_lt_of_pos coefficient (by omega)
  change (if residue * 2 ≤ (modulus : Int) then residue
    else residue - (modulus : Int)).natAbs * 2 ≤ modulus
  split
  next hhalf =>
    have habs : (residue.natAbs : Int) = residue :=
      Int.natAbs_of_nonneg hresidueNonnegative
    rw [← habs] at hhalf
    exact_mod_cast hhalf
  next hhalf =>
    have hdifferenceNonpositive : residue - (modulus : Int) ≤ 0 := by omega
    have hstrict : 2 * ((modulus : Int) - residue) < modulus := by omega
    have habs : ((residue - (modulus : Int)).natAbs : Int) =
        (modulus : Int) - residue := by
      calc
        ((residue - (modulus : Int)).natAbs : Int) =
            ((-(residue - (modulus : Int))).natAbs : Int) := by
              rw [Int.natAbs_neg]
        _ = -(residue - (modulus : Int)) :=
          Int.natAbs_of_nonneg (neg_nonneg.mpr hdifferenceNonpositive)
        _ = (modulus : Int) - residue := by ring
    rw [← habs] at hstrict
    have hcast : 2 * (residue - (modulus : Int)).natAbs ≤ modulus := by
      exact_mod_cast (le_of_lt hstrict)
    simpa [Nat.mul_comm] using hcast

/-- Exact coefficient cells emitted by the generated symmetric-mod loop. -/
theorem symmetricModLoop_toList (input : SparsePolyZZ) (modulus : ZZ)
    (index : Nat) (result output : SparsePolyZZ)
    (hrun : Generated.StrictRecombine.symmetricModLoop input modulus index
      result = .ok output) :
    output.toList = result.toList ++
      (input.toList.drop index).filterMap (fun term =>
        let coefficient := ZZ.symmetricMod term.2 modulus
        if coefficient = 0 then none else some (term.1, coefficient)) := by
  induction hmeasure : input.size - index using Nat.strong_induction_on
      generalizing index result output with
  | h measure ih =>
      rw [Generated.StrictRecombine.symmetricModLoop] at hrun
      split at hrun
      next hindex =>
        dsimp at hrun
        let coefficient := ZZ.symmetricMod input[index].2 modulus
        by_cases hzero : coefficient = 0
        · change ZZ.symmetricMod input[index].2 modulus = 0 at hzero
          simp only [hzero, if_true] at hrun
          rw [ih (input.size - (index + 1)) (by omega)
            (index + 1) result output hrun rfl]
          have hsuffix : input.toList.drop index = input[index] ::
              input.toList.drop (index + 1) := by
            simpa using List.drop_eq_getElem_cons
              (l := input.toList) (i := index) (by simpa using hindex)
          rw [hsuffix, List.filterMap_cons]
          dsimp only
          rw [if_pos hzero]
        · change ZZ.symmetricMod input[index].2 modulus ≠ 0 at hzero
          simp only [hzero, if_false] at hrun
          rw [ih (input.size - (index + 1)) (by omega)
            (index + 1)
            (result.push
              (input[index].1, ZZ.symmetricMod input[index].2 modulus))
            output hrun rfl]
          have hsuffix : input.toList.drop index = input[index] ::
              input.toList.drop (index + 1) := by
            simpa using List.drop_eq_getElem_cons
              (l := input.toList) (i := index) (by simpa using hindex)
          rw [hsuffix, List.filterMap_cons]
          dsimp only
          rw [if_neg hzero, Array.toList_push, List.append_assoc]
          rfl
      next hindex =>
        have hle : input.size ≤ index := Nat.le_of_not_gt hindex
        have hout := Except.ok.inj hrun
        subst output
        simp [List.drop_eq_nil_iff.mpr hle]

/-- Every coefficient physically emitted by the generated symmetric-mod loop
satisfies the source's closed half-modulus bound. -/
theorem symmetricModLoop_coefficients_bounded
    (input output : SparsePolyZZ) (modulus : Nat) (hmodulus : 0 < modulus)
    (hrun : Generated.StrictRecombine.symmetricModLoop input (modulus : ZZ)
      0 #[] = .ok output) :
    ∀ term ∈ output.toList, term.2.natAbs * 2 ≤ modulus := by
  have hlist := symmetricModLoop_toList input (modulus : ZZ) 0 #[] output hrun
  simp only [Array.toList_empty, List.nil_append, List.drop_zero] at hlist
  intro term hterm
  rw [hlist, List.mem_filterMap] at hterm
  rcases hterm with ⟨source, hsource, houtput⟩
  change (if ZZ.symmetricMod source.2 (modulus : ZZ) = 0 then none
    else some (source.1, ZZ.symmetricMod source.2 (modulus : ZZ))) =
      some term at houtput
  split at houtput
  next hzero => contradiction
  next hzero =>
    injection houtput with htermEq
    subst term
    exact symmetricMod_natAbs_mul_two_le source.2 modulus hmodulus

theorem symmetricModLoop_canonical (input output : SparsePolyZZ)
    (modulus : ZZ)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical input)
    (hrun : Generated.StrictRecombine.symmetricModLoop input modulus 0 #[] =
      .ok output) :
    StrictPolynomialMod.SparsePolyZZCanonical output := by
  let reduceTerm : UMonomial × Int → Option (UMonomial × Int) :=
    fun term =>
      let coefficient := ZZ.symmetricMod term.2 modulus
      if coefficient = 0 then none else some (term.1, coefficient)
  have hlist := symmetricModLoop_toList input modulus 0 #[] output hrun
  have houtput : output.toList = input.toList.filterMap reduceTerm := by
    simpa [reduceTerm] using hlist
  unfold StrictPolynomialMod.SparsePolyZZCanonical
  rw [houtput]
  constructor
  · letI : Trans
        (fun left right : UMonomial × Int => left.1.deg > right.1.deg)
        (fun left right : UMonomial × Int => left.1.deg > right.1.deg)
        (fun left right : UMonomial × Int => left.1.deg > right.1.deg) :=
      ⟨by intro left middle right hleft hright
          exact Nat.lt_trans hright hleft⟩
    exact (hcanonical.1.pairwise.filterMap reduceTerm (by
      intro left right hdegree left' hleft right' hright
      dsimp [reduceTerm] at hleft hright
      split at hleft <;> try contradiction
      split at hright <;> try contradiction
      injection hleft with hleft
      injection hright with hright
      subst left'
      subst right'
      exact hdegree)).isChain
  · intro term hterm
    rw [List.mem_filterMap] at hterm
    obtain ⟨source, hsource, hreduce⟩ := hterm
    dsimp [reduceTerm] at hreduce
    split at hreduce
    next hzero => contradiction
    next hzero =>
      injection hreduce with heq
      subst term
      exact hzero

theorem symmetricModLoop_toPolyMod (input : SparsePolyZZ)
    (modulus : Nat) (hmodulus : 0 < modulus) (index : Nat)
    (result output : SparsePolyZZ)
    (hrun : Generated.StrictRecombine.symmetricModLoop input (modulus : ZZ)
      index result = .ok output) :
    Refinement.StrictHensel.toPolyMod modulus output =
      Refinement.StrictHensel.toPolyMod modulus result +
        Refinement.StrictHensel.termsToPolyMod modulus
          (input.toList.drop index) := by
  induction hmeasure : input.size - index using Nat.strong_induction_on
      generalizing index result output with
  | h measure ih =>
      rw [Generated.StrictRecombine.symmetricModLoop] at hrun
      split at hrun
      next hindex =>
        dsimp at hrun
        let coefficient := ZZ.symmetricMod input[index].2 (modulus : ZZ)
        by_cases hzero : coefficient = 0
        · change ZZ.symmetricMod input[index].2 (modulus : ZZ) = 0 at hzero
          simp only [hzero, if_true] at hrun
          rw [ih (input.size - (index + 1)) (by omega)
            (index + 1) result output hrun rfl]
          have hsuffix : input.toList.drop index = input[index] ::
              input.toList.drop (index + 1) := by
            simpa using List.drop_eq_getElem_cons
              (l := input.toList) (i := index) (by simpa using hindex)
          have hcast := symmetricMod_cast modulus hmodulus input[index].2
          change (ZZ.symmetricMod input[index].2 (modulus : ZZ) : ZMod modulus) =
            (input[index].2 : ZMod modulus) at hcast
          rw [hzero] at hcast
          rw [hsuffix, Refinement.StrictHensel.termsToPolyMod_cons]
          rw [← hcast]
          simp
        · change ZZ.symmetricMod input[index].2 (modulus : ZZ) ≠ 0 at hzero
          simp only [hzero, if_false] at hrun
          rw [ih (input.size - (index + 1)) (by omega)
            (index + 1) (result.push (input[index].1,
              ZZ.symmetricMod input[index].2 (modulus : ZZ))) output
            hrun rfl]
          rw [Refinement.StrictHensel.toPolyMod_push]
          have hsuffix : input.toList.drop index = input[index] ::
              input.toList.drop (index + 1) := by
            simpa using List.drop_eq_getElem_cons
              (l := input.toList) (i := index) (by simpa using hindex)
          rw [symmetricMod_cast modulus hmodulus input[index].2]
          rw [hsuffix, Refinement.StrictHensel.termsToPolyMod_cons]
          abel
      next hindex =>
        have hle : input.size ≤ index := Nat.le_of_not_gt hindex
        have hout := Except.ok.inj hrun
        subst output
        simp [Refinement.StrictHensel.termsToPolyMod,
          List.drop_eq_nil_iff.mpr hle]

/-- The actual symmetric representative chosen modulo `modulus` has the same
polynomial value modulo every positive divisor `base`. -/
theorem symmetricModLoop_toPolyMod_of_dvd (input : SparsePolyZZ)
    (modulus base : Nat) (hmodulus : 0 < modulus) (hbase : 0 < base)
    (hdivides : base ∣ modulus) (index : Nat)
    (result output : SparsePolyZZ)
    (hrun : Generated.StrictRecombine.symmetricModLoop input
      (modulus : ZZ) index result = .ok output) :
    Refinement.StrictHensel.toPolyMod base output =
      Refinement.StrictHensel.toPolyMod base result +
        Refinement.StrictHensel.termsToPolyMod base
          (input.toList.drop index) := by
  induction hmeasure : input.size - index using Nat.strong_induction_on
      generalizing index result output with
  | h measure ih =>
      rw [Generated.StrictRecombine.symmetricModLoop] at hrun
      split at hrun
      next hindex =>
        dsimp at hrun
        let coefficient := ZZ.symmetricMod input[index].2 (modulus : ZZ)
        by_cases hzero : coefficient = 0
        · change ZZ.symmetricMod input[index].2 (modulus : ZZ) = 0 at hzero
          simp only [hzero, if_true] at hrun
          rw [ih (input.size - (index + 1)) (by omega)
            (index + 1) result output hrun rfl]
          have hsuffix : input.toList.drop index = input[index] ::
              input.toList.drop (index + 1) := by
            simpa using List.drop_eq_getElem_cons
              (l := input.toList) (i := index) (by simpa using hindex)
          have hcast := symmetricMod_cast_of_dvd base modulus hbase hmodulus
            hdivides input[index].2
          change (ZZ.symmetricMod input[index].2 (modulus : ZZ) : ZMod base) =
            (input[index].2 : ZMod base) at hcast
          rw [hzero] at hcast
          rw [hsuffix, Refinement.StrictHensel.termsToPolyMod_cons]
          rw [← hcast]
          simp
        · change ZZ.symmetricMod input[index].2 (modulus : ZZ) ≠ 0 at hzero
          simp only [hzero, if_false] at hrun
          rw [ih (input.size - (index + 1)) (by omega)
            (index + 1) (result.push (input[index].1,
              ZZ.symmetricMod input[index].2 (modulus : ZZ))) output
            hrun rfl]
          rw [Refinement.StrictHensel.toPolyMod_push]
          have hsuffix : input.toList.drop index = input[index] ::
              input.toList.drop (index + 1) := by
            simpa using List.drop_eq_getElem_cons
              (l := input.toList) (i := index) (by simpa using hindex)
          rw [symmetricMod_cast_of_dvd base modulus hbase hmodulus hdivides
            input[index].2]
          rw [hsuffix, Refinement.StrictHensel.termsToPolyMod_cons]
          abel
      next hindex =>
        have hle : input.size ≤ index := Nat.le_of_not_gt hindex
        have hout := Except.ok.inj hrun
        subst output
        simp [Refinement.StrictHensel.termsToPolyMod,
          List.drop_eq_nil_iff.mpr hle]

theorem symmetricModRaw_toPolyMod (input output : SparsePolyZZ)
    (modulus : Nat) (hmodulus : 0 < modulus)
    (hrun : Generated.StrictRecombine.symmetricModRaw input (modulus : ZZ) =
      .ok output) :
    Refinement.StrictHensel.toPolyMod modulus output =
      Refinement.StrictHensel.toPolyMod modulus input := by
  unfold Generated.StrictRecombine.symmetricModRaw at hrun
  rw [dif_pos (by exact_mod_cast hmodulus)] at hrun
  simpa [Refinement.StrictHensel.toPolyMod_eq_termsToPolyMod] using
    symmetricModLoop_toPolyMod input modulus hmodulus 0 #[] output hrun

theorem symmetricModRaw_coefficients_bounded
    (input output : SparsePolyZZ) (modulus : Nat) (hmodulus : 0 < modulus)
    (hrun : Generated.StrictRecombine.symmetricModRaw input (modulus : ZZ) =
      .ok output) :
    ∀ term ∈ output.toList, term.2.natAbs * 2 ≤ modulus := by
  unfold Generated.StrictRecombine.symmetricModRaw at hrun
  rw [dif_pos (by exact_mod_cast hmodulus)] at hrun
  exact symmetricModLoop_coefficients_bounded input output modulus hmodulus hrun

/-- Recovery with the exact boundary asymmetry produced by the C++ path:
the symmetric representative may attain `m/2`, while the true scaled factor
is strictly inside the half interval.  One strict side is sufficient for
uniqueness modulo `m`, including even moduli. -/
theorem symmetric_recovery_closed_left
    (representative target : Polynomial Int) (modulus : Nat)
    (hmodulus : 0 < modulus)
    (hmod : Polynomial.map (Int.castRingHom (ZMod modulus)) representative =
      Polynomial.map (Int.castRingHom (ZMod modulus)) target)
    (hrepresentative : ∀ degree,
      (representative.coeff degree).natAbs * 2 ≤ modulus)
    (htarget : ∀ degree,
      (target.coeff degree).natAbs * 2 < modulus) :
    representative = target := by
  ext degree
  let left := representative.coeff degree
  let right := target.coeff degree
  have hcoefficient : (left : ZMod modulus) = (right : ZMod modulus) := by
    have hcoefficient' := congrArg (fun polynomial => polynomial.coeff degree)
      hmod
    simpa only [Polynomial.coeff_map] using hcoefficient'
  have hdivides : (modulus : Int) ∣ left - right := by
    rw [← ZMod.intCast_zmod_eq_zero_iff_dvd]
    push_cast
    rw [sub_eq_zero]
    exact hcoefficient
  have habs : (left - right).natAbs < modulus := by
    calc
      (left - right).natAbs ≤ left.natAbs + right.natAbs :=
        Int.natAbs_sub_le left right
      _ < modulus := by
        have hleft := hrepresentative degree
        have hright := htarget degree
        omega
  obtain ⟨multiple, hmultiple⟩ := hdivides
  have habsProduct : (left - right).natAbs =
      modulus * multiple.natAbs := by
    rw [hmultiple, Int.natAbs_mul, Int.natAbs_natCast]
  rw [habsProduct] at habs
  have hmultipleZero : multiple.natAbs = 0 := by
    by_contra hnonzero
    exact Nat.not_lt.mpr
      (Nat.le_mul_of_pos_right modulus (Nat.pos_of_ne_zero hnonzero)) habs
  have hdifference : left - right = 0 := by
    rw [Int.natAbs_eq_zero.mp hmultipleZero, mul_zero] at hmultiple
    exact hmultiple
  exact Int.eq_of_sub_eq_zero hdifference

/-- Scalar specialization of the concrete C++ symmetric representative.  A
coefficient strictly inside the half-modulus interval is recovered literally,
including its sign; this follows from the same closed-left uniqueness argument
used by polynomial reconstruction. -/
theorem symmetricMod_eq_of_strict_bound
    (coefficient : Int) (modulus : Nat) (hmodulus : 0 < modulus)
    (hbound : coefficient.natAbs * 2 < modulus) :
    ZZ.symmetricMod coefficient (modulus : Int) = coefficient := by
  have hpolynomial := symmetric_recovery_closed_left
    (Polynomial.C (ZZ.symmetricMod coefficient (modulus : Int)))
    (Polynomial.C coefficient) modulus hmodulus
    (by
      ext degree
      by_cases hdegree : degree = 0
      · subst degree
        simpa using symmetricMod_cast modulus hmodulus coefficient
      · rw [Polynomial.map_C, Polynomial.map_C,
          Polynomial.coeff_C_ne_zero hdegree,
          Polynomial.coeff_C_ne_zero hdegree])
    (by
      intro degree
      by_cases hdegree : degree = 0
      · subst degree
        simpa using symmetricMod_natAbs_mul_two_le coefficient modulus hmodulus
      · rw [Polynomial.coeff_C_ne_zero hdegree]
        simp)
    (by
      intro degree
      by_cases hdegree : degree = 0
      · simpa [hdegree] using hbound
      · rw [Polynomial.coeff_C_ne_zero hdegree]
        simpa using hmodulus)
  have hcoefficient := congrArg (fun polynomial : Polynomial Int =>
    polynomial.coeff 0) hpolynomial
  simpa using hcoefficient

/-- Congruent-input form of scalar symmetric recovery.  The generated
representative is taken from `input`, while the strict half-modulus bound is
proved for the intended integer `target`. -/
theorem symmetricMod_eq_of_congruent_strict_bound
    (input target : Int) (modulus : Nat) (hmodulus : 0 < modulus)
    (hcongruent : (input : ZMod modulus) = (target : ZMod modulus))
    (hbound : target.natAbs * 2 < modulus) :
    ZZ.symmetricMod input (modulus : Int) = target := by
  have hpolynomial := symmetric_recovery_closed_left
    (Polynomial.C (ZZ.symmetricMod input (modulus : Int)))
    (Polynomial.C target) modulus hmodulus
    (by
      ext degree
      by_cases hdegree : degree = 0
      · subst degree
        simp only [Polynomial.map_C, Polynomial.coeff_C_zero]
        exact (symmetricMod_cast modulus hmodulus input).trans hcongruent
      · rw [Polynomial.map_C, Polynomial.map_C,
          Polynomial.coeff_C_ne_zero hdegree,
          Polynomial.coeff_C_ne_zero hdegree])
    (by
      intro degree
      by_cases hdegree : degree = 0
      · subst degree
        simpa using symmetricMod_natAbs_mul_two_le input modulus hmodulus
      · rw [Polynomial.coeff_C_ne_zero hdegree]
        simp)
    (by
      intro degree
      by_cases hdegree : degree = 0
      · simpa [hdegree] using hbound
      · rw [Polynomial.coeff_C_ne_zero hdegree]
        simpa using hmodulus)
  have hcoefficient := congrArg
    (fun polynomial : Polynomial Int => polynomial.coeff 0) hpolynomial
  simpa using hcoefficient

private theorem foldl_int_coefficient_squares
    (terms : List (UMonomial × Int)) (accumulator : Int) :
    terms.foldl (fun sum term => sum + term.2 * term.2) accumulator =
      accumulator + (terms.map fun term => term.2 * term.2).sum := by
  induction terms generalizing accumulator with
  | nil => simp
  | cons head tail ih =>
      rw [List.foldl_cons, ih]
      simp
      ring

/-- The generated range-for L2-square helper is exactly the source-order sum
of the squares of the physically stored integer coefficients. -/
theorem upolyNormL2SqRaw_eq_stored_sum (poly : SparsePolyZZ) :
    Generated.StrictHensel.__upoly_norm_l2_sq_upoly_raw_ir poly =
      (poly.toList.map fun term => term.2 * term.2).sum := by
  unfold Generated.StrictHensel.__upoly_norm_l2_sq_upoly_raw_ir
  rw [← Array.foldl_toList]
  simpa using foldl_int_coefficient_squares poly.toList 0

theorem upolyNormL2SqRaw_nonnegative (poly : SparsePolyZZ) :
    0 ≤ Generated.StrictHensel.__upoly_norm_l2_sq_upoly_raw_ir poly := by
  rw [upolyNormL2SqRaw_eq_stored_sum]
  apply List.sum_nonneg
  intro value hvalue
  rw [List.mem_map] at hvalue
  rcases hvalue with ⟨term, hterm, rfl⟩
  exact mul_self_nonneg term.2

/-- Exact source recurrence for the generated multiplicative binomial loop.
The checked integer division is justified at each iteration by the standard
binomial successor identity. -/
theorem binomialLoopRaw_eq_choose (n k i : Nat)
    (hi : i ≤ k) (hk : k ≤ n) :
    Generated.StrictHensel.__binomial_loop_raw_ir n k i (n.choose i : Int) =
      (n.choose k : Int) := by
  induction hmeasure : k - i using Nat.strong_induction_on generalizing i with
  | h measure ih =>
      rw [Generated.StrictHensel.__binomial_loop_raw_ir]
      split
      next hcontinue =>
        have hiN : i ≤ n := hi.trans hk
        have hchoose := Nat.choose_succ_right_eq n i
        have hstep : (n.choose i : Int) * ((n : Int) - i) =
            ((i : Int) + 1) * (n.choose (i + 1) : Int) := by
          have hcast : (n.choose i : Int) * ((n : Int) - i) =
              (n.choose (i + 1) : Int) * ((i : Int) + 1) := by
            exact_mod_cast hchoose.symm
          simpa [mul_comm] using hcast
        rw [hstep, Int.mul_ediv_cancel_left _ (by omega)]
        exact ih (k - (i + 1)) (by omega) (i + 1) (by omega) rfl
      next hstop =>
        have hieq : i = k := by omega
        subst i
        rfl

/-- The exact machine-integer call made by `__mignotte_bound` reaches the
well-founded natural-number binomial loop without signed wraparound,
truncation, or negative clamping. -/
theorem binomialRaw_degree_half_eq_choose (degree : Nat)
    (hdegree : degree < 2 ^ 63) :
    Generated.StrictHensel.__binomial_raw_ir
        (Int64.ofInt (degree : Int))
        (Int64.ofInt (degree : Int) / 2) =
      (degree.choose (degree / 2) : Int) := by
  have htwo : 2 < 2 ^ 63 := by norm_num
  have hhalf : degree / 2 < 2 ^ 63 :=
    (Nat.div_le_self degree 2).trans_lt hdegree
  have hdivision : Int64.ofNat degree / (2 : Int64) =
      Int64.ofNat (degree / 2) := by
    simpa using (Int64.ofNat_div hdegree htwo).symm
  rw [Int64.ofInt_eq_ofNat, hdivision]
  let n : Int64 := Int64.ofNat degree
  let k : Int64 := Int64.ofNat (degree / 2)
  have hnInt : n.toInt = degree := by
    exact Int64.toInt_ofNat_of_lt hdegree
  have hkInt : k.toInt = degree / 2 := by
    exact Int64.toInt_ofNat_of_lt hhalf
  have hnNonnegative : ¬(n < 0) := by
    rw [Int64.lt_iff_toInt_lt, hnInt, Int64.toInt_zero]
    omega
  have hkNonnegative : ¬(k < 0) := by
    rw [Int64.lt_iff_toInt_lt, hkInt, Int64.toInt_zero]
    omega
  have hkLeN : ¬(n < k) := by
    intro hlt
    rw [Int64.lt_iff_toInt_lt, hnInt, hkInt] at hlt
    omega
  have hnClamp : n.toNatClampNeg = degree := by
    exact Int64.toNatClampNeg_ofNat_of_lt hdegree
  have hkClamp : k.toNatClampNeg = degree / 2 := by
    exact Int64.toNatClampNeg_ofNat_of_lt hhalf
  change Generated.StrictHensel.__binomial_raw_ir n k = _
  unfold Generated.StrictHensel.__binomial_raw_ir
  rw [if_neg (by simp [hkNonnegative, hkLeN])]
  by_cases hkZero : degree / 2 = 0
  · have hkEqZero : k = 0 := by
      apply Int64.toInt.inj
      rw [hkInt, Int64.toInt_zero]
      exact_mod_cast hkZero
    rw [if_pos (by simp [hkEqZero])]
    simp [hkZero]
  · have hkNeZero : k ≠ 0 := by
      intro heq
      have := congrArg Int64.toInt heq
      rw [hkInt, Int64.toInt_zero] at this
      exact hkZero (by exact_mod_cast this)
    have hkLtDegree : degree / 2 < degree := by omega
    have hkNeN : k ≠ n := by
      intro heq
      have := congrArg Int64.toInt heq
      rw [hkInt, hnInt] at this
      omega
    rw [if_neg (by simp [hkNeZero, hkNeN])]
    have hkNatLe : degree / 2 ≤ degree := Nat.div_le_self degree 2
    have hsub : n - k = Int64.ofNat (degree - degree / 2) := by
      simpa [n, k] using (Int64.ofNat_sub degree (degree / 2) hkNatLe).symm
    have hnotSymmetric : ¬(n - k < k) := by
      rw [hsub]
      intro hlt
      rw [Int64.lt_iff_toInt_lt, hkInt,
        Int64.toInt_ofNat_of_lt (Nat.sub_lt_of_lt hdegree)] at hlt
      omega
    rw [if_neg hnotSymmetric]
    change Generated.StrictHensel.__binomial_loop_raw_ir
      n.toNatClampNeg k.toNatClampNeg 0 1 = _
    rw [hnClamp, hkClamp]
    simpa using
      (binomialLoopRaw_eq_choose degree (degree / 2) 0 (by omega) hkNatLe)

/-- The generated well-founded Newton loop is extensionally the kernel's
verified natural square-root iterator.  This theorem compares the two real
recurrences; it does not replace the generated execution. -/
theorem isqrtCeilLoopRaw_eq_sqrtIter (n current : Nat) :
    Generated.StrictHensel.__isqrt_ceil_loop_raw_ir n current =
      Nat.sqrt.iter n current := by
  induction current using Nat.strong_induction_on with
  | h current ih =>
      rw [Generated.StrictHensel.__isqrt_ceil_loop_raw_ir]
      unfold Nat.sqrt.iter
      by_cases hzero : current = 0
      · subst current
        simp
      · rw [dif_neg hzero]
        let next := (current + n / current) / 2
        by_cases hnext : next < current
        · simp only [next, dif_pos hnext]
          exact ih next hnext
        · simp only [next, dif_neg hnext]

/-- The concrete generated Newton loop returns the floor square root whenever
its actual starting estimate has the standard one-step upper bound. -/
theorem isqrtCeilLoopRaw_floor_bounds (n initial : Nat)
    (hinitial : n < (initial + 1) * (initial + 1)) :
    let root := Generated.StrictHensel.__isqrt_ceil_loop_raw_ir n initial
    root * root ≤ n ∧ n < (root + 1) * (root + 1) := by
  rw [isqrtCeilLoopRaw_eq_sqrtIter]
  exact ⟨Nat.sqrt.iter_sq_le n initial,
    Nat.sqrt.lt_iter_succ_sq n initial hinitial⟩

/-- The complete generated integer square-root helper is nonnegative and its
square bounds the actual integer input from above.  The proof executes its
bit-length initialisation, Newton recurrence, and final correction branch. -/
theorem isqrtCeilRaw_nonnegative_and_square_ge (value : Int) :
    0 ≤ Generated.StrictHensel.__isqrt_ceil_raw_ir value ∧
      value ≤ Generated.StrictHensel.__isqrt_ceil_raw_ir value *
        Generated.StrictHensel.__isqrt_ceil_raw_ir value := by
  rw [Generated.StrictHensel.__isqrt_ceil_raw_ir]
  by_cases hnonpositive : value ≤ 0
  · rw [dif_pos hnonpositive]
    simp [hnonpositive]
  · rw [dif_neg hnonpositive]
    have hpositive : 0 < value := lt_of_not_ge hnonpositive
    have hnonzero : value ≠ 0 := ne_of_gt hpositive
    let bits := ZZ.sizeinbase_nat value 2
    let initial : Nat := 2 ^ ((bits + 1) / 2)
    let root := Generated.StrictHensel.__isqrt_ceil_loop_raw_ir
      value.natAbs initial
    have hbits : value.natAbs < 2 ^ bits := by
      exact ZZ.lt_pow_sizeinbase_nat hnonzero (by decide)
    have hexponent : bits ≤ 2 * ((bits + 1) / 2) := by omega
    have hpow : 2 ^ bits ≤ initial * initial := by
      calc
        2 ^ bits ≤ 2 ^ (2 * ((bits + 1) / 2)) :=
          Nat.pow_le_pow_right (by decide) hexponent
        _ = initial * initial := by
          rw [show 2 * ((bits + 1) / 2) =
            (bits + 1) / 2 + (bits + 1) / 2 by omega, pow_add]
    have hinitial : value.natAbs < (initial + 1) * (initial + 1) := by
      have hinitialPositive : 0 < initial := by positivity
      exact lt_of_lt_of_le hbits (hpow.trans (by nlinarith))
    have hroot := isqrtCeilLoopRaw_floor_bounds value.natAbs initial hinitial
    change 0 ≤ (if (root : Int) * (root : Int) < value then
        (root : Int) + 1 else (root : Int)) ∧
      value ≤ (if (root : Int) * (root : Int) < value then
          (root : Int) + 1 else (root : Int)) *
        (if (root : Int) * (root : Int) < value then
          (root : Int) + 1 else (root : Int))
    have habs : (value.natAbs : Int) = value := by
      rw [Int.natAbs_of_nonneg hpositive.le]
    by_cases hbelow : (root : Int) * root < value
    · rw [if_pos hbelow]
      constructor
      · positivity
      · have hupper : value < ((root + 1 : Nat) : Int) * (root + 1) := by
          rw [← habs]
          exact_mod_cast hroot.2
        norm_num at hupper ⊢
        omega
    · rw [if_neg hbelow]
      constructor
      · positivity
      · exact le_of_not_gt hbelow

/-- The complete successful generated Mignotte helper exposes the exact
central binomial coefficient computed by its machine call, while retaining
the generated norm fold and generated Newton square root verbatim. -/
theorem mignotteBoundRaw_eq_choose_isqrt (f : SparsePolyZZ)
    (leading : UMonomial × Int) (hleading : f[0]? = some leading)
    (hdegree : leading.1.deg < 2 ^ 63) :
    Generated.StrictHensel.__mignotte_bound_upoly_raw_ir f = .ok
      ((leading.1.deg.choose (leading.1.deg / 2) : Int) *
        Generated.StrictHensel.__isqrt_ceil_raw_ir
          (Generated.StrictHensel.__upoly_norm_l2_sq_upoly_raw_ir f)) := by
  rw [Generated.StrictHensel.__mignotte_bound_upoly_raw_ir]
  simp only [hleading]
  rw [binomialRaw_degree_half_eq_choose leading.1.deg hdegree]

/-- The concrete generated square-root component used in the Mignotte bound
is a nonnegative upper square root of the concrete stored-coefficient fold. -/
theorem mignotteNormRaw_nonnegative_and_square_ge (f : SparsePolyZZ) :
    let norm := Generated.StrictHensel.__isqrt_ceil_raw_ir
      (Generated.StrictHensel.__upoly_norm_l2_sq_upoly_raw_ir f)
    0 ≤ norm ∧
      Generated.StrictHensel.__upoly_norm_l2_sq_upoly_raw_ir f ≤ norm * norm :=
  isqrtCeilRaw_nonnegative_and_square_ge
    (Generated.StrictHensel.__upoly_norm_l2_sq_upoly_raw_ir f)

private theorem intTermsToPoly_coeff_eq_zero_of_degrees_ne
    (terms : List (UMonomial × Int)) (degree : Nat)
    (hne : ∀ term ∈ terms, term.1.deg ≠ degree) :
    (intTermsToPoly terms).coeff degree = 0 := by
  induction terms with
  | nil => simp [intTermsToPoly]
  | cons head tail ih =>
      have hhead := hne head (by simp)
      have htail : ∀ term ∈ tail, term.1.deg ≠ degree := by
        intro term hterm
        exact hne term (by simp [hterm])
      change (Polynomial.monomial head.1.deg head.2).coeff degree +
        (intTermsToPoly tail).coeff degree = 0
      rw [Polynomial.coeff_monomial, if_neg hhead, ih htail]
      simp

/-- A strictly degree-ordered sparse list has exactly the same coefficient
square sum as its mathematical polynomial over any range containing every
stored degree.  This is the representation bridge needed by the concrete
Mignotte norm fold. -/
private theorem intTerms_square_sum_eq_range
    (terms : List (UMonomial × Int)) (limit : Nat)
    (hchain : List.IsChain
      (fun left right : UMonomial × Int => left.1.deg > right.1.deg) terms)
    (hdegrees : ∀ term ∈ terms, term.1.deg ≤ limit) :
    (terms.map fun term => term.2 * term.2).sum =
      (Finset.range (limit + 1)).sum fun degree =>
        (intTermsToPoly terms).coeff degree *
          (intTermsToPoly terms).coeff degree := by
  induction terms with
  | nil => simp [intTermsToPoly]
  | cons head tail ih =>
      have htailChain := List.IsChain.tail hchain
      have htailDegrees : ∀ term ∈ tail, term.1.deg ≤ limit := by
        intro term hterm
        exact hdegrees term (by simp [hterm])
      have hheadDegree : head.1.deg ≤ limit :=
        hdegrees head (by simp)
      have htailBelow : ∀ term ∈ tail, term.1.deg < head.1.deg :=
        sparseZZ_chain_head_gt_all head tail hchain
      have htailAtHead :
          (intTermsToPoly tail).coeff head.1.deg = 0 := by
        apply intTermsToPoly_coeff_eq_zero_of_degrees_ne
        intro term hterm
        exact Nat.ne_of_lt (htailBelow term hterm)
      rw [List.map_cons, List.sum_cons, ih htailChain htailDegrees]
      symm
      calc
        (Finset.range (limit + 1)).sum (fun degree =>
            (intTermsToPoly (head :: tail)).coeff degree *
              (intTermsToPoly (head :: tail)).coeff degree) =
            (Finset.range (limit + 1)).sum (fun degree =>
              (if degree = head.1.deg then head.2 * head.2 else 0) +
                (intTermsToPoly tail).coeff degree *
                  (intTermsToPoly tail).coeff degree) := by
              apply Finset.sum_congr rfl
              intro degree hdegree
              change ((Polynomial.monomial head.1.deg head.2 +
                intTermsToPoly tail).coeff degree) *
                  ((Polynomial.monomial head.1.deg head.2 +
                    intTermsToPoly tail).coeff degree) = _
              rw [Polynomial.coeff_add, Polynomial.coeff_monomial]
              by_cases heq : head.1.deg = degree
              · subst degree
                simp [htailAtHead]
              · simp [heq, Ne.symm heq]
        _ = (Finset.range (limit + 1)).sum (fun degree =>
              if degree = head.1.deg then head.2 * head.2 else 0) +
            (Finset.range (limit + 1)).sum (fun degree =>
              (intTermsToPoly tail).coeff degree *
                (intTermsToPoly tail).coeff degree) :=
              Finset.sum_add_distrib
        _ = head.2 * head.2 +
            (Finset.range (limit + 1)).sum (fun degree =>
              (intTermsToPoly tail).coeff degree *
                (intTermsToPoly tail).coeff degree) := by
              rw [Finset.sum_ite_eq']
              simp [hheadDegree]

theorem symmetricModRaw_toPolyMod_of_dvd (input output : SparsePolyZZ)
    (modulus base : Nat) (hmodulus : 0 < modulus) (hbase : 0 < base)
    (hdivides : base ∣ modulus)
    (hrun : Generated.StrictRecombine.symmetricModRaw input (modulus : ZZ) =
      .ok output) :
    Refinement.StrictHensel.toPolyMod base output =
      Refinement.StrictHensel.toPolyMod base input := by
  unfold Generated.StrictRecombine.symmetricModRaw at hrun
  rw [dif_pos (by exact_mod_cast hmodulus)] at hrun
  simpa [Refinement.StrictHensel.toPolyMod_eq_termsToPolyMod] using
    symmetricModLoop_toPolyMod_of_dvd input modulus base hmodulus hbase
      hdivides 0 #[] output hrun

theorem symmetricModRaw_canonical (input output : SparsePolyZZ)
    (modulus : Nat) (hmodulus : 0 < modulus)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical input)
    (hrun : Generated.StrictRecombine.symmetricModRaw input (modulus : ZZ) =
      .ok output) :
    StrictPolynomialMod.SparsePolyZZCanonical output := by
  unfold Generated.StrictRecombine.symmetricModRaw at hrun
  rw [dif_pos (by exact_mod_cast hmodulus)] at hrun
  exact symmetricModLoop_canonical input output (modulus : ZZ) hcanonical hrun

/-- The generated content loop is exactly the source-order gcd fold over the
remaining concrete coefficient cells.  This is the executable basis for
proving that `primitiveRaw` really returns a primitive polynomial; no
primitive-part specification is used to choose its result. -/
theorem contentLoop_eq_foldl (input : SparsePolyZZ) (index acc : Nat) :
    Generated.StrictRecombine.contentLoop input index acc =
      (input.toList.drop index).foldl
        (fun current term => Nat.gcd current term.2.natAbs) acc := by
  induction hmeasure : input.size - index using Nat.strong_induction_on
      generalizing index acc with
  | h measure ih =>
      rw [Generated.StrictRecombine.contentLoop]
      split
      next hindex =>
        have hsuffix : input.toList.drop index = input[index] ::
            input.toList.drop (index + 1) := by
          simpa using List.drop_eq_getElem_cons
            (l := input.toList) (i := index) (by simpa using hindex)
        rw [hsuffix, List.foldl_cons]
        exact ih (input.size - (index + 1)) (by omega)
          (index + 1) (Nat.gcd acc input[index].2.natAbs) rfl
      next hindex =>
        have hle : input.size ≤ index := Nat.le_of_not_gt hindex
        rw [List.drop_eq_nil_iff.mpr (by simpa using hle)]
        rfl

private theorem dvd_foldl_coefficient_gcd_iff
    (divisor accumulator : Nat) (terms : List (UMonomial × Int)) :
    divisor ∣ terms.foldl
        (fun current term => Nat.gcd current term.2.natAbs) accumulator ↔
      divisor ∣ accumulator ∧ ∀ term ∈ terms, divisor ∣ term.2.natAbs := by
  induction terms generalizing accumulator with
  | nil => simp
  | cons head tail ih =>
      rw [List.foldl_cons, ih]
      constructor
      · rintro ⟨hhead, htail⟩
        exact ⟨(Nat.dvd_gcd_iff.mp hhead).1, fun term hterm => by
          rcases List.mem_cons.mp hterm with rfl | htailMember
          · exact (Nat.dvd_gcd_iff.mp hhead).2
          · exact htail term htailMember⟩
      · rintro ⟨haccumulator, hall⟩
        exact ⟨Nat.dvd_gcd_iff.mpr
            ⟨haccumulator, hall head (List.mem_cons_self)⟩,
          fun term hterm => hall term (List.mem_cons_of_mem head hterm)⟩

/-- Divisibility characterization of the actual generated content loop.  It
states both halves of the gcd contract over the exact remaining source cells. -/
theorem contentLoop_dvd_iff (input : SparsePolyZZ) (index accumulator divisor : Nat) :
    divisor ∣ Generated.StrictRecombine.contentLoop input index accumulator ↔
      divisor ∣ accumulator ∧
        ∀ term ∈ input.toList.drop index, divisor ∣ term.2.natAbs := by
  rw [contentLoop_eq_foldl]
  exact dvd_foldl_coefficient_gcd_iff divisor accumulator
    (input.toList.drop index)

private theorem sparseZZ_sum_coeff_of_mem
    (terms : List (UMonomial × Int))
    (hchain : List.IsChain
      (fun left right : UMonomial × Int => left.1.deg > right.1.deg) terms)
    (term : UMonomial × Int) (hterm : term ∈ terms) :
    ((terms.map fun entry =>
      Polynomial.monomial entry.1.deg entry.2).sum).coeff term.1.deg =
        term.2 := by
  induction terms with
  | nil => simp at hterm
  | cons head tail ih =>
      have htailChain : List.IsChain
          (fun left right : UMonomial × Int =>
            left.1.deg > right.1.deg) tail :=
        List.IsChain.tail hchain
      have htailBelow : ∀ entry ∈ tail, entry.1.deg < head.1.deg :=
        sparseZZ_chain_head_gt_all head tail hchain
      rcases List.mem_cons.mp hterm with rfl | htailMember
      · simp only [List.map_cons, List.sum_cons, Polynomial.coeff_add,
          Polynomial.coeff_monomial, if_pos rfl]
        rw [sparseZZTail_coeff_zero_above term.1.deg tail (by
          simpa using htailBelow), add_zero]
        simp
      · have hdegrees : head.1.deg ≠ term.1.deg := by
          exact Nat.ne_of_gt (htailBelow term htailMember)
        simp only [List.map_cons, List.sum_cons, Polynomial.coeff_add,
          Polynomial.coeff_monomial, if_neg hdegrees, zero_add]
        exact ih htailChain htailMember

/-- Canonical sparse storage has no duplicate degree, so every concrete
stored coefficient is exactly the mathematical coefficient at that degree. -/
theorem sparsePolyZZ_toPoly_coeff_of_mem
    (input : SparsePolyZZ)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical input)
    (term : UMonomial × Int) (hterm : term ∈ input.toList) :
    (SparsePolyZZ.toPoly input).coeff term.1.deg = term.2 := by
  unfold SparsePolyZZ.toPoly
  exact sparseZZ_sum_coeff_of_mem input.toList hcanonical.1 term hterm

/-- The first physical sparse term carries the exact mathematical degree. -/
theorem sparsePolyZZ_natDegree_eq_head (poly : SparsePolyZZ)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical poly)
    (hnonempty : 0 < poly.size) :
    (SparsePolyZZ.toPoly poly).natDegree = poly[0].1.deg := by
  have hheadMem : poly[0] ∈ poly.toList := by
    exact Array.getElem_mem_toList hnonempty
  have hheadCoeff :
      (SparsePolyZZ.toPoly poly).coeff poly[0].1.deg = poly[0].2 :=
    sparsePolyZZ_toPoly_coeff_of_mem poly hcanonical poly[0] hheadMem
  have hheadNonzero : poly[0].2 ≠ 0 := hcanonical.2 poly[0] hheadMem
  apply Polynomial.natDegree_eq_of_le_of_coeff_ne_zero
  · rw [Polynomial.natDegree_le_iff_coeff_eq_zero]
    intro degree hdegree
    apply intTermsToPoly_coeff_eq_zero_of_degrees_ne
    intro term hterm heq
    have htermLe : term.1.deg ≤ poly[0].1.deg := by
      have hlistNonempty : poly.toList ≠ [] := by
        intro hempty
        have hlength := congrArg List.length hempty
        have hsizeZero : poly.size = 0 := by simpa using hlength
        omega
      have hhead : poly.toList.head hlistNonempty = poly[0] := by
        rw [List.head_eq_getElem]
        simpa using Array.getElem_toList hnonempty
      have hheadFirst : poly.toList = poly[0] :: poly.toList.tail := by
        rw [← hhead]
        exact (List.cons_head_tail hlistNonempty).symm
      have hchain : List.IsChain
          (fun left right : UMonomial × Int => left.1.deg > right.1.deg)
          (poly[0] :: poly.toList.tail) := by
        rw [← hheadFirst]
        exact hcanonical.1
      rw [hheadFirst] at hterm
      rcases List.mem_cons.mp hterm with heqHead | htail
      · simpa [heqHead]
      · exact (sparseZZ_chain_head_gt_all poly[0] poly.toList.tail hchain
          term htail).le
    omega
  · rw [hheadCoeff]
    exact hheadNonzero

/-- For a canonical sparse integer polynomial, the exact source-order norm
fold equals the mathematical L2 coefficient-square sum over its degree
range.  Strict stored degree order prevents duplicate coefficients. -/
theorem upolyNormL2SqRaw_eq_coeff_range_sum (poly : SparsePolyZZ)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical poly) :
    Generated.StrictHensel.__upoly_norm_l2_sq_upoly_raw_ir poly =
      (Finset.range ((SparsePolyZZ.toPoly poly).natDegree + 1)).sum
        (fun degree =>
          (SparsePolyZZ.toPoly poly).coeff degree *
            (SparsePolyZZ.toPoly poly).coeff degree) := by
  rw [upolyNormL2SqRaw_eq_stored_sum]
  have hdegrees : ∀ term ∈ poly.toList,
      term.1.deg ≤ (SparsePolyZZ.toPoly poly).natDegree := by
    intro term hterm
    apply Polynomial.le_natDegree_of_ne_zero
    rw [sparsePolyZZ_toPoly_coeff_of_mem poly hcanonical term hterm]
    exact hcanonical.2 term hterm
  simpa [SparsePolyZZ.toPoly, intTermsToPoly] using
    intTerms_square_sum_eq_range poly.toList
      (SparsePolyZZ.toPoly poly).natDegree hcanonical.1 hdegrees

/-- Every coefficient of a genuine integer divisor is bounded by the exact
central binomial and Newton norm computed by the generated C++ Mignotte
helper.  This converts the existing real L2 theorem back to an integer
inequality using the canonical sparse storage equality above. -/
theorem divisor_coeff_le_generated_mignotte_norm (f : SparsePolyZZ)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (g : Polynomial Int) (hf : SparsePolyZZ.toPoly f ≠ 0)
    (hg : g ∣ SparsePolyZZ.toPoly f) :
    let norm := Generated.StrictHensel.__isqrt_ceil_raw_ir
      (Generated.StrictHensel.__upoly_norm_l2_sq_upoly_raw_ir f)
    ∀ degree,
      ((g.coeff degree).natAbs : Int) ≤
        ((SparsePolyZZ.toPoly f).natDegree.choose
          ((SparsePolyZZ.toPoly f).natDegree / 2) : Int) * norm := by
  let source := SparsePolyZZ.toPoly f
  let squareSum := Generated.StrictHensel.__upoly_norm_l2_sq_upoly_raw_ir f
  let norm := Generated.StrictHensel.__isqrt_ceil_raw_ir squareSum
  have hnorm := mignotteNormRaw_nonnegative_and_square_ge f
  have hnorm' : 0 ≤ norm ∧ squareSum ≤ norm * norm := by
    simpa [squareSum, norm] using hnorm
  have hsquareSum : squareSum =
      (Finset.range (source.natDegree + 1)).sum
        (fun degree => source.coeff degree * source.coeff degree) := by
    simpa [source, squareSum] using
      upolyNormL2SqRaw_eq_coeff_range_sum f hcanonical
  have hrealSum :
      (Finset.range (source.natDegree + 1)).sum
          (fun degree => ((source.coeff degree).natAbs : Real) ^ 2) =
        (squareSum : Real) := by
    rw [hsquareSum]
    push_cast
    apply Finset.sum_congr rfl
    intro degree hdegree
    norm_num [sq]
  have hsqrt :
      Real.sqrt ((Finset.range (source.natDegree + 1)).sum
        (fun degree => ((source.coeff degree).natAbs : Real) ^ 2)) ≤
          (norm : Real) := by
    rw [hrealSum, Real.sqrt_le_left (by exact_mod_cast hnorm'.1)]
    have hpow : squareSum ≤ norm ^ 2 := by simpa [pow_two] using hnorm'.2
    exact_mod_cast hpow
  change ∀ degree : Nat, ((g.coeff degree).natAbs : Int) ≤
    ((source.natDegree.choose (source.natDegree / 2) : Nat) : Int) * norm
  intro degree
  have hbound := mignotte_bound_l2 source g (by simpa [source] using hf)
    (by simpa [source] using hg) degree
  have hreal : ((g.coeff degree).natAbs : Real) ≤
      ((source.natDegree.choose (source.natDegree / 2) : Nat) : Real) *
        (norm : Real) := hbound.trans (by gcongr)
  rw [Nat.cast_natAbs (α := Real), Int.cast_abs] at hreal
  apply (@Int.cast_le Real _ _ _).mp
  norm_num
  exact hreal

/-- The actual successful generated Mignotte call returns one concrete
nonnegative integer that bounds every coefficient of every genuine divisor
of the represented source polynomial. -/
theorem mignotteBoundRaw_bounds_divisor (f : SparsePolyZZ)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hnonempty : 0 < f.size) (hdegree : f[0].1.deg < 2 ^ 63)
    (g : Polynomial Int) (hf : SparsePolyZZ.toPoly f ≠ 0)
    (hg : g ∣ SparsePolyZZ.toPoly f) :
    ∃ bound,
      Generated.StrictHensel.__mignotte_bound_upoly_raw_ir f = .ok bound ∧
      0 ≤ bound ∧ ∀ degree, ((g.coeff degree).natAbs : Int) ≤ bound := by
  let norm := Generated.StrictHensel.__isqrt_ceil_raw_ir
    (Generated.StrictHensel.__upoly_norm_l2_sq_upoly_raw_ir f)
  let bound : Int :=
    (f[0].1.deg.choose (f[0].1.deg / 2) : Int) * norm
  have hleading : f[0]? = some f[0] := Array.getElem?_eq_getElem hnonempty
  have hrun := mignotteBoundRaw_eq_choose_isqrt f f[0] hleading hdegree
  have hnatDegree := sparsePolyZZ_natDegree_eq_head f hcanonical hnonempty
  have hnorm := mignotteNormRaw_nonnegative_and_square_ge f
  have hnormNonnegative : 0 ≤ norm := by simpa [norm] using hnorm.1
  have hboundNonnegative : 0 ≤ bound := by
    dsimp [bound]
    exact mul_nonneg (by positivity) hnormNonnegative
  refine ⟨bound, by simpa [bound, norm] using hrun,
    hboundNonnegative, ?_⟩
  intro degree
  have hcoefficient := divisor_coeff_le_generated_mignotte_norm
    f hcanonical g hf hg degree
  simpa [bound, norm, hnatDegree] using hcoefficient

/-- When the actual generated Hensel entry uses its default (`aTarget = 0`)
precision, the modulus returned by its well-founded lifting loop is already
strictly large enough to recover every coefficient of the leading-coefficient
scaled genuine divisor.  Both the coefficient bound and the final modulus are
values computed by the generated C++ lowering. -/
theorem hensel_output_modulus_bounds_scaled_divisor
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {factors : Array SparsePolyZp} {p : UInt64}
    [Fact (Nat.Prime p.toNat)]
    {aTarget : Int32} {output : Array SparsePolyZZ × ZZ}
    (hcorrect : Refinement.StrictHensel.HenselLiftEntryCorrect termination
      f factors p aTarget output)
    (hzero : aTarget = 0)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hnonempty : 0 < f.size) (hdegree : f[0].1.deg < 2 ^ 63)
    (leading : UMonomial × ZZ) (hleading : f[0]? = some leading)
    (g : Polynomial Int) (hf : SparsePolyZZ.toPoly f ≠ 0)
    (hg : g ∣ SparsePolyZZ.toPoly f) :
    0 < output.2 ∧ ∀ degree,
      ((leading.2 * g.coeff degree).natAbs : Int) * 2 < output.2 := by
  subst aTarget
  rcases mignotteBoundRaw_bounds_divisor f hcanonical hnonempty hdegree g hf hg
    with ⟨bound, hboundRun, hboundNonnegative, hbound⟩
  have htargetNonnegative : ∀ target,
      Refinement.StrictHensel.HenselLiftTargetCorrect f p 0 target →
        0 ≤ target := by
    intro target htarget
    cases htarget with
    | mignotte sourceLeading htargetZero hsourceLeading =>
        have hsame : sourceLeading = leading := by
          rw [hleading] at hsourceLeading
          exact Option.some.inj hsourceLeading.symm
        subst sourceLeading
        rw [Generated.StrictHensel.__mignotte_bound_upoly_raw_ir,
          hleading] at hboundRun
        have hboundValue := Except.ok.inj hboundRun
        rw [hboundValue]
        exact mul_nonneg
          (mul_nonneg (by norm_num) (Int.ofNat_zero_le leading.2.natAbs))
          hboundNonnegative
    | explicit hpositive =>
        exact False.elim ((by decide : ¬ (0 : Int32) > 0) hpositive)
  rcases hcorrect.outputModulus_gt_target htargetNonnegative with
    ⟨target, htarget, htargetLt⟩
  have htargetNonnegative' := htargetNonnegative target htarget
  constructor
  · exact lt_of_le_of_lt htargetNonnegative' htargetLt
  · intro degree
    cases htarget with
    | mignotte sourceLeading htargetZero hsourceLeading =>
        have hsame : sourceLeading = leading := by
          rw [hleading] at hsourceLeading
          exact Option.some.inj hsourceLeading.symm
        subst sourceLeading
        rw [Generated.StrictHensel.__mignotte_bound_upoly_raw_ir,
          hleading] at hboundRun
        have hboundValue := Except.ok.inj hboundRun
        rw [hboundValue] at htargetLt
        have hcoefficient := hbound degree
        rw [Int.natAbs_mul]
        norm_num only [Nat.cast_mul]
        have hscaled :
            (leading.2.natAbs : Int) * (g.coeff degree).natAbs ≤
              (leading.2.natAbs : Int) * bound :=
          mul_le_mul_of_nonneg_left hcoefficient
            (Int.ofNat_zero_le leading.2.natAbs)
        nlinarith
    | explicit hpositive =>
        exact False.elim ((by decide : ¬ (0 : Int32) > 0) hpositive)

/-- A bound on every physically stored coefficient of a canonical sparse
polynomial extends to every mathematical coefficient. -/
theorem sparsePolyZZ_coefficients_bounded_of_stored
    (poly : SparsePolyZZ)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical poly)
    (bound : Nat)
    (hstored : ∀ term ∈ poly.toList, term.2.natAbs * 2 ≤ bound) :
    ∀ degree, ((SparsePolyZZ.toPoly poly).coeff degree).natAbs * 2 ≤ bound := by
  intro degree
  by_cases hzero : (SparsePolyZZ.toPoly poly).coeff degree = 0
  · simp [hzero]
  · have hexists : ∃ term ∈ poly.toList, term.1.deg = degree := by
      by_contra hnone
      push Not at hnone
      apply hzero
      simpa [SparsePolyZZ.toPoly, intTermsToPoly] using
        intTermsToPoly_coeff_eq_zero_of_degrees_ne poly.toList degree hnone
    rcases hexists with ⟨term, hterm, hdegree⟩
    have hcoefficient := sparsePolyZZ_toPoly_coeff_of_mem poly hcanonical
      term hterm
    rw [hdegree] at hcoefficient
    rw [hcoefficient]
    exact hstored term hterm

/-- Polynomial-level closed half-modulus bound for the exact generated
`symmetricModRaw` result. -/
theorem symmetricModRaw_toPoly_coefficients_bounded
    (input output : SparsePolyZZ) (modulus : Nat) (hmodulus : 0 < modulus)
    (hinputCanonical : StrictPolynomialMod.SparsePolyZZCanonical input)
    (hrun : Generated.StrictRecombine.symmetricModRaw input (modulus : ZZ) =
      .ok output) :
    ∀ degree,
      ((SparsePolyZZ.toPoly output).coeff degree).natAbs * 2 ≤ modulus := by
  have houtputCanonical := symmetricModRaw_canonical input output modulus
    hmodulus hinputCanonical hrun
  exact sparsePolyZZ_coefficients_bounded_of_stored output houtputCanonical
    modulus (symmetricModRaw_coefficients_bounded input output modulus
      hmodulus hrun)

/-- The concrete generated `symmetricModRaw` execution recovers a unique
integer target whenever its input is congruent to that target and the target
lies strictly inside the half-modulus interval.  This is the executable
recovery lemma used after the generated Hensel precision proof. -/
theorem symmetricModRaw_recovers_strictly_bounded_target
    (input output : SparsePolyZZ) (target : Polynomial Int)
    (modulus : Nat) (hmodulus : 0 < modulus)
    (hinputCanonical : StrictPolynomialMod.SparsePolyZZCanonical input)
    (hrun : Generated.StrictRecombine.symmetricModRaw input (modulus : ZZ) =
      .ok output)
    (hcongruent : Polynomial.map (Int.castRingHom (ZMod modulus))
        (SparsePolyZZ.toPoly input) =
      Polynomial.map (Int.castRingHom (ZMod modulus)) target)
    (htarget : ∀ degree,
      (target.coeff degree).natAbs * 2 < modulus) :
    SparsePolyZZ.toPoly output = target := by
  apply symmetric_recovery_closed_left (SparsePolyZZ.toPoly output) target
    modulus hmodulus
  · calc
      Polynomial.map (Int.castRingHom (ZMod modulus))
          (SparsePolyZZ.toPoly output) =
          Polynomial.map (Int.castRingHom (ZMod modulus))
            (SparsePolyZZ.toPoly input) := by
              simpa [Refinement.StrictHensel.toPolyMod] using
                symmetricModRaw_toPolyMod input output modulus hmodulus hrun
      _ = Polynomial.map (Int.castRingHom (ZMod modulus)) target := hcongruent
  · exact symmetricModRaw_toPoly_coefficients_bounded input output modulus
      hmodulus hinputCanonical hrun
  · exact htarget

private theorem sparseZZ_sum_coeff_dvd
    (terms : List (UMonomial × Int)) (divisor : Int)
    (hall : ∀ term ∈ terms, divisor ∣ term.2) (degree : Nat) :
    divisor ∣ ((terms.map fun entry =>
      Polynomial.monomial entry.1.deg entry.2).sum).coeff degree := by
  induction terms with
  | nil => simp
  | cons head tail ih =>
      simp only [List.map_cons, List.sum_cons, Polynomial.coeff_add,
        Polynomial.coeff_monomial]
      apply dvd_add
      · split
        · exact hall head List.mem_cons_self
        · exact dvd_zero divisor
      · exact ih (fun term hterm =>
          hall term (List.mem_cons_of_mem head hterm))

private theorem sparsePolyZZ_toPoly_coeff_dvd
    (input : SparsePolyZZ) (divisor : Int)
    (hall : ∀ term ∈ input.toList, divisor ∣ term.2) (degree : Nat) :
    divisor ∣ (SparsePolyZZ.toPoly input).coeff degree := by
  unfold SparsePolyZZ.toPoly
  exact sparseZZ_sum_coeff_dvd input.toList divisor hall degree

/-- The nonnegative gcd returned by the exact generated content loop is the
normalized mathematical content of the canonical sparse polynomial. -/
theorem contentLoop_zero_eq_content
    (input : SparsePolyZZ)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical input) :
    (Generated.StrictRecombine.contentLoop input 0 0 : Int) =
      (SparsePolyZZ.toPoly input).content := by
  let divisor := Generated.StrictRecombine.contentLoop input 0 0
  have hloopSelf : divisor ∣ divisor := dvd_refl divisor
  have hdividesStoredNat : ∀ term ∈ input.toList,
      divisor ∣ term.2.natAbs := by
    have hcharacterization :=
      (contentLoop_dvd_iff input 0 0 divisor).mp hloopSelf
    simpa using hcharacterization.2
  have hdividesStoredInt : ∀ term ∈ input.toList,
      (divisor : Int) ∣ term.2 := by
    intro term hterm
    rw [← Int.natAbs_dvd_natAbs]
    simpa using hdividesStoredNat term hterm
  have hdivisorContent : (divisor : Int) ∣
      (SparsePolyZZ.toPoly input).content := by
    rw [Polynomial.content, Finset.dvd_gcd_iff]
    intro degree hdegree
    exact sparsePolyZZ_toPoly_coeff_dvd input (divisor : Int)
      hdividesStoredInt degree
  have hcontentDividesStored : ∀ term ∈ input.toList,
      (SparsePolyZZ.toPoly input).content.natAbs ∣ term.2.natAbs := by
    intro term hterm
    rw [Int.natAbs_dvd_natAbs]
    rw [← sparsePolyZZ_toPoly_coeff_of_mem input hcanonical term hterm]
    exact Polynomial.content_dvd_coeff term.1.deg
  have hcontentNatDivisor :
      (SparsePolyZZ.toPoly input).content.natAbs ∣ divisor := by
    apply (contentLoop_dvd_iff input 0 0
      (SparsePolyZZ.toPoly input).content.natAbs).mpr
    exact ⟨dvd_zero _, by simpa using hcontentDividesStored⟩
  have hcontentDivisor : (SparsePolyZZ.toPoly input).content ∣
      (divisor : Int) := by
    rw [← Int.natAbs_dvd_natAbs]
    simpa using hcontentNatDivisor
  exact dvd_antisymm_of_normalize_eq
    (Int.normalize_coe_nat divisor) Polynomial.normalize_content
      hdivisorContent hcontentDivisor

theorem primitiveDivideLoop_toPoly (input : SparsePolyZZ) (divisor : Int)
    (index : Nat) (result output : SparsePolyZZ)
    (hrun : Generated.StrictRecombine.primitiveDivideLoop input divisor index
      result = .ok output) :
    Polynomial.C divisor * SparsePolyZZ.toPoly output =
      Polynomial.C divisor * SparsePolyZZ.toPoly result +
        intTermsToPoly (input.toList.drop index) := by
  induction hmeasure : input.size - index using Nat.strong_induction_on
      generalizing index result output with
  | h measure ih =>
      rw [Generated.StrictRecombine.primitiveDivideLoop] at hrun
      split at hrun
      next hindex =>
        split at hrun
        next hdivisor =>
          split at hrun
          next hdivides =>
            rw [ih (input.size - (index + 1)) (by omega)
              (index + 1)
              (result.push (input[index].1, input[index].2 / divisor))
              output hrun rfl]
            rw [SparsePolyZZ.toPoly]
            simp only [Array.toList_push, List.map_append, List.map_singleton,
              List.sum_append, List.sum_singleton]
            have hsuffix : input.toList.drop index = input[index] ::
                input.toList.drop (index + 1) := by
              simpa using List.drop_eq_getElem_cons
                (l := input.toList) (i := index) (by simpa using hindex)
            rw [hsuffix]
            simp only [intTermsToPoly, List.map_cons, List.sum_cons]
            rw [mul_add, Polynomial.C_mul_monomial]
            have hcancel : divisor * (input[index].2 / divisor) =
                input[index].2 := by
              exact Int.mul_ediv_cancel_of_dvd hdivides
            rw [hcancel]
            change Polynomial.C divisor * SparsePolyZZ.toPoly result + _ + _ =
              Polynomial.C divisor * SparsePolyZZ.toPoly result + (_ + _)
            abel
          next hdivides => contradiction
        next hdivisor => contradiction
      next hindex =>
        have hle : input.size ≤ index := Nat.le_of_not_gt hindex
        have hout := Except.ok.inj hrun
        subst output
        simp [intTermsToPoly, List.drop_eq_nil_iff.mpr hle]

/-- Exact list produced by the source coefficient-wise primitive division
loop.  The successful trace itself certifies every checked divisibility. -/
theorem primitiveDivideLoop_toList (input : SparsePolyZZ) (divisor : Int)
    (index : Nat) (result output : SparsePolyZZ)
    (hrun : Generated.StrictRecombine.primitiveDivideLoop input divisor index
      result = .ok output) :
    output.toList = result.toList ++
      (input.toList.drop index).map
        (fun term => (term.1, term.2 / divisor)) := by
  induction hmeasure : input.size - index using Nat.strong_induction_on
      generalizing index result output with
  | h measure ih =>
      rw [Generated.StrictRecombine.primitiveDivideLoop] at hrun
      split at hrun
      next hindex =>
        split at hrun
        next hdivisor =>
          split at hrun
          next hdivides =>
            rw [ih (input.size - (index + 1)) (by omega)
              (index + 1)
              (result.push (input[index].1, input[index].2 / divisor))
              output hrun rfl]
            have hsuffix : input.toList.drop index = input[index] ::
                input.toList.drop (index + 1) := by
              simpa using List.drop_eq_getElem_cons
                (l := input.toList) (i := index) (by simpa using hindex)
            simp [hsuffix, List.append_assoc]
          next hdivides => contradiction
        next hdivisor => contradiction
      next hindex =>
        have hle : input.size ≤ index := Nat.le_of_not_gt hindex
        have hout := Except.ok.inj hrun
        subst output
        simp [List.drop_eq_nil_iff.mpr hle]

/-- Every coefficient copied by a successful primitive-division loop is an
exact quotient by a nonzero divisor. -/
theorem primitiveDivideLoop_constraints (input : SparsePolyZZ) (divisor : Int)
    (index : Nat) (result output : SparsePolyZZ)
    (hrun : Generated.StrictRecombine.primitiveDivideLoop input divisor index
      result = .ok output) :
    ∀ term ∈ input.toList.drop index, divisor ≠ 0 ∧ divisor ∣ term.2 := by
  induction hmeasure : input.size - index using Nat.strong_induction_on
      generalizing index result output with
  | h measure ih =>
      rw [Generated.StrictRecombine.primitiveDivideLoop] at hrun
      split at hrun
      next hindex =>
        split at hrun
        next hdivisor =>
          split at hrun
          next hdivides =>
            have htail := ih (input.size - (index + 1)) (by omega)
              (index + 1)
              (result.push (input[index].1, input[index].2 / divisor))
              output hrun rfl
            have hsuffix : input.toList.drop index = input[index] ::
                input.toList.drop (index + 1) := by
              simpa using List.drop_eq_getElem_cons
                (l := input.toList) (i := index) (by simpa using hindex)
            intro term hterm
            rw [hsuffix] at hterm
            rcases List.mem_cons.mp hterm with rfl | htailMember
            · exact ⟨hdivisor, hdivides⟩
            · exact htail term htailMember
          next hdivides => contradiction
        next hdivisor => contradiction
      next hindex =>
        intro term hterm
        exact False.elim (by
          rw [List.drop_eq_nil_iff.mpr (Nat.le_of_not_gt hindex)] at hterm
          simp at hterm)

/-- Coefficient-wise primitive division executes whenever its concrete
nonzero divisor divides every remaining stored coefficient. -/
theorem primitiveDivideLoop_complete (input : SparsePolyZZ) (divisor : Int)
    (index : Nat) (result : SparsePolyZZ) (hdivisor : divisor ≠ 0)
    (hdivides : ∀ term ∈ input.toList.drop index, divisor ∣ term.2) :
    ∃ output, Generated.StrictRecombine.primitiveDivideLoop input divisor
      index result = .ok output := by
  induction hmeasure : input.size - index using Nat.strong_induction_on
      generalizing index result with
  | h measure ih =>
      rw [Generated.StrictRecombine.primitiveDivideLoop]
      split
      next hindex =>
        have hsuffix : input.toList.drop index = input[index] ::
            input.toList.drop (index + 1) := by
          simpa using List.drop_eq_getElem_cons
            (l := input.toList) (i := index) (by simpa using hindex)
        have hhead : divisor ∣ input[index].2 := by
          exact hdivides input[index] (by simp [hsuffix])
        rw [dif_pos hhead]
        exact ih (input.size - (index + 1)) (by omega) (index + 1)
          (result.push (input[index].1, input[index].2 / divisor))
          (fun term hterm => hdivides term (by simp [hsuffix, hterm])) rfl
      next hindex => exact ⟨result, rfl⟩

/-- The exact generated primitive-part entry succeeds on every canonical
sparse polynomial.  In the nonempty branch its content is proved nonzero and
to divide every physical coefficient before the source loop is executed. -/
theorem primitiveRaw_complete (input : SparsePolyZZ)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical input) :
    ∃ content primitive,
      Generated.StrictRecombine.primitiveRaw input = .ok (content, primitive) := by
  by_cases hempty : input.isEmpty = true
  · exact ⟨1, input, by simp [Generated.StrictRecombine.primitiveRaw, hempty]⟩
  · have hnonempty : 0 < input.size := by
      have hsize : input.size ≠ 0 := by
        simpa [Array.isEmpty] using hempty
      omega
    let content : Int := Generated.StrictRecombine.contentLoop input 0 0
    let divisor : Int := if input[0].2 < 0 then -content else content
    have hpolyNe : SparsePolyZZ.toPoly input ≠ 0 := by
      intro hzero
      have hleading := sparsePolyZZ_leadingCoeff_eq_head input hcanonical
        hnonempty
      rw [hzero] at hleading
      exact (hcanonical.2 input[0]
        (Array.getElem_mem_toList hnonempty)) (by simpa using hleading.symm)
    have hcontentEq : content = (SparsePolyZZ.toPoly input).content := by
      exact contentLoop_zero_eq_content input hcanonical
    have hcontentNe : content ≠ 0 := by
      rw [hcontentEq]
      exact fun hzero => hpolyNe (Polynomial.content_eq_zero_iff.mp hzero)
    have hdivisorNe : divisor ≠ 0 := by
      simp only [divisor]
      split <;> simp_all
    have hdivides : ∀ term ∈ input.toList, divisor ∣ term.2 := by
      intro term hterm
      have hcoeff := Polynomial.content_dvd_coeff
        (p := SparsePolyZZ.toPoly input) term.1.deg
      rw [sparsePolyZZ_toPoly_coeff_of_mem input hcanonical term hterm] at hcoeff
      simp only [divisor]
      split
      · rw [hcontentEq]
        exact neg_dvd.mpr hcoeff
      · simpa [hcontentEq] using hcoeff
    rcases primitiveDivideLoop_complete input divisor 0 #[] hdivisorNe
        (by simpa using hdivides) with ⟨primitive, hprimitive⟩
    refine ⟨divisor, primitive, ?_⟩
    simp only [Generated.StrictRecombine.primitiveRaw, hempty, ↓reduceDIte]
    simpa [content, divisor, hprimitive]

theorem primitiveRaw_toPoly (input primitive : SparsePolyZZ) (content : Int)
    (hrun : Generated.StrictRecombine.primitiveRaw input =
      .ok (content, primitive)) :
    SparsePolyZZ.toPoly input =
      Polynomial.C content * SparsePolyZZ.toPoly primitive := by
  unfold Generated.StrictRecombine.primitiveRaw at hrun
  split at hrun
  next hempty =>
    have hinput : input = #[] := Array.isEmpty_iff.mp hempty
    subst input
    have hout := Except.ok.inj hrun
    cases hout
    simp [SparsePolyZZ.toPoly]
  next hempty =>
    dsimp at hrun
    split at hrun
    next fault hdivide => contradiction
    next primitive' hdivide =>
      have hout : (content, primitive) =
          (if (input[0]'(by
              have : input.size ≠ 0 := by simpa [Array.isEmpty] using hempty
              omega)).2 < 0 then
            -(Generated.StrictRecombine.contentLoop input 0 0 : Int)
          else (Generated.StrictRecombine.contentLoop input 0 0 : Int), primitive') :=
        (Except.ok.inj hrun).symm
      have hcontent := congrArg Prod.fst hout
      have hprimitive := congrArg Prod.snd hout
      have hsemantic := primitiveDivideLoop_toPoly input
        (if (input[0]'(by
            have : input.size ≠ 0 := by simpa [Array.isEmpty] using hempty
            omega)).2 < 0 then
          -(Generated.StrictRecombine.contentLoop input 0 0 : Int)
        else (Generated.StrictRecombine.contentLoop input 0 0 : Int)) 0 #[]
        primitive' hdivide
      have hcontent' : content =
          (if (input[0]'(by
              have : input.size ≠ 0 := by simpa [Array.isEmpty] using hempty
              omega)).2 < 0 then
            -(Generated.StrictRecombine.contentLoop input 0 0 : Int)
          else (Generated.StrictRecombine.contentLoop input 0 0 : Int)) := hcontent
      rw [← hcontent'] at hsemantic
      have hprimitive' : primitive = primitive' := hprimitive
      rw [← hprimitive'] at hsemantic
      simpa [SparsePolyZZ.toPoly, intTermsToPoly] using hsemantic.symm

/-- The actual primitive-part routine preserves the canonical sparse degree
order and cannot create a zero coefficient. -/
theorem primitiveRaw_canonical (input primitive : SparsePolyZZ) (content : Int)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical input)
    (hrun : Generated.StrictRecombine.primitiveRaw input =
      .ok (content, primitive)) :
    StrictPolynomialMod.SparsePolyZZCanonical primitive := by
  unfold Generated.StrictRecombine.primitiveRaw at hrun
  split at hrun
  next hempty =>
    have hout := Except.ok.inj hrun
    injection hout with hcontent hprimitive
    subst primitive
    exact hcanonical
  next hempty =>
    dsimp at hrun
    let divisor : Int := if (input[0]'(by
      have : input.size ≠ 0 := by simpa [Array.isEmpty] using hempty
      omega)).2 < 0 then
        -(Generated.StrictRecombine.contentLoop input 0 0 : Int)
      else Generated.StrictRecombine.contentLoop input 0 0
    cases hdivide : Generated.StrictRecombine.primitiveDivideLoop input
        divisor 0 #[] with
    | error fault => simp [divisor, hdivide] at hrun
    | ok primitive' =>
        simp only [divisor, hdivide] at hrun
        have hout := Except.ok.inj hrun
        injection hout with hcontent hprimitive
        subst primitive
        have hlist := primitiveDivideLoop_toList input divisor 0 #[]
          primitive' hdivide
        simp only [Array.toList_empty, List.nil_append, List.drop_zero] at hlist
        have hconstraints := primitiveDivideLoop_constraints input divisor 0 #[]
          primitive' hdivide
        simp only [List.drop_zero] at hconstraints
        unfold StrictPolynomialMod.SparsePolyZZCanonical
        rw [hlist]
        constructor
        · rw [List.isChain_map]
          simpa using hcanonical.1
        · intro output houtput
          rcases List.mem_map.mp houtput with ⟨term, hterm, rfl⟩
          have hinputNonzero := hcanonical.2 term hterm
          rcases hconstraints term hterm with ⟨hdivisor, hdivides⟩
          simp only [Prod.snd]
          intro hquotient
          have hcancel : divisor * (term.2 / divisor) = term.2 :=
            Int.mul_ediv_cancel_of_dvd hdivides
          rw [hquotient] at hcancel
          simp at hcancel
          exact hinputNonzero hcancel.symm

/-- On every nonempty canonical input, the exact generated content-gcd and
coefficient-division execution returns a mathematically primitive polynomial. -/
theorem primitiveRaw_isPrimitive
    (input primitive : SparsePolyZZ) (content : Int)
    (hnonempty : 0 < input.size)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical input)
    (hrun : Generated.StrictRecombine.primitiveRaw input =
      .ok (content, primitive)) :
    (SparsePolyZZ.toPoly primitive).IsPrimitive := by
  have hsemantic := primitiveRaw_toPoly input primitive content hrun
  have hloopContent := contentLoop_zero_eq_content input hcanonical
  have hinputNe : SparsePolyZZ.toPoly input ≠ 0 := by
    intro hzero
    have hleading := sparsePolyZZ_leadingCoeff_eq_head input hcanonical hnonempty
    rw [hzero] at hleading
    have hheadNonzero := hcanonical.2 input[0]
      (Array.getElem_mem_toList hnonempty)
    exact hheadNonzero (by simpa using hleading.symm)
  have hcontentNonzero : (SparsePolyZZ.toPoly input).content ≠ 0 :=
    fun hzero => hinputNe (Polynomial.content_eq_zero_iff.mp hzero)
  have hnormalizeContent : normalize content =
      (SparsePolyZZ.toPoly input).content := by
    rw [Generated.StrictRecombine.primitiveRaw] at hrun
    split at hrun
    next hempty =>
      have hsizeZero : input.size = 0 := by
        simpa [Array.isEmpty] using hempty
      omega
    next hempty =>
      dsimp at hrun
      split at hrun
      next fault hdivide => contradiction
      next primitive' hdivide =>
        have hout := Except.ok.inj hrun
        have hcontent := congrArg Prod.fst hout
        simp only [Prod.fst] at hcontent
        rw [← hcontent]
        split
        next hnegative =>
          rw [Int.normalize_of_nonpos (by
            exact neg_nonpos.mpr (Int.natCast_nonneg _))]
          rw [neg_neg, hloopContent]
        next hnonnegative =>
          rw [Int.normalize_coe_nat, hloopContent]
  have hcontentEquation := congrArg Polynomial.content hsemantic
  rw [Polynomial.content_C_mul, hnormalizeContent] at hcontentEquation
  apply Polynomial.isPrimitive_iff_content_eq_one.mpr
  apply mul_left_cancel₀ hcontentNonzero
  simpa using hcontentEquation.symm

theorem subtractScaledTermsLoop_toPoly (divisor : SparsePolyZZ)
    (scale : Int) (degreeShift index : Nat) (terms : SparsePolyZZ) :
    intTermsToPoly
        (Generated.StrictRecombine.subtractScaledTermsLoop divisor scale
          degreeShift index terms).toList =
      intTermsToPoly terms.toList -
        Polynomial.monomial degreeShift scale *
          intTermsToPoly (divisor.toList.drop index) := by
  induction hmeasure : divisor.size - index using Nat.strong_induction_on
      generalizing index terms with
  | h measure ih =>
      rw [Generated.StrictRecombine.subtractScaledTermsLoop]
      split
      next hindex =>
        rw [ih (divisor.size - (index + 1)) (by omega)
          (index + 1)
          (terms.push (⟨divisor[index].1.deg + degreeShift⟩,
            -(scale * divisor[index].2))) rfl]
        have hsuffix : divisor.toList.drop index = divisor[index] ::
            divisor.toList.drop (index + 1) := by
          simpa using List.drop_eq_getElem_cons
            (l := divisor.toList) (i := index) (by simpa using hindex)
        simp [intTermsToPoly, hsuffix]
        rw [mul_add, Polynomial.monomial_mul_monomial]
        rw [add_comm divisor[index].1.deg degreeShift]
        abel
      next hindex =>
        have hle : divisor.size ≤ index := Nat.le_of_not_gt hindex
        simp [intTermsToPoly, List.drop_eq_nil_iff.mpr hle]

theorem subtractScaledNormalize_toPoly (remainder divisor : SparsePolyZZ)
    (scale : Int) (degreeShift : Nat) :
    SparsePolyZZ.toPoly
        (Generated.StrictRecombine.subtractScaledNormalize remainder divisor
          scale degreeShift) =
      SparsePolyZZ.toPoly remainder -
        Polynomial.monomial degreeShift scale * SparsePolyZZ.toPoly divisor := by
  unfold Generated.StrictRecombine.subtractScaledNormalize
  rw [normalization_toPoly]
  change intTermsToPoly
      (Generated.StrictRecombine.subtractScaledTermsLoop divisor scale
        degreeShift 0 remainder).toList = _
  rw [subtractScaledTermsLoop_toPoly]
  simp [intTermsToPoly, SparsePolyZZ.toPoly]

/-- For canonical nonempty sparse inputs, genuine polynomial divisibility
forces the exact C++ head-coefficient divisibility test to succeed. -/
theorem canonical_head_coefficient_dvd_of_poly_dvd
    (divisor remainder : SparsePolyZZ)
    (hdivisorCanonical :
      StrictPolynomialMod.SparsePolyZZCanonical divisor)
    (hremainderCanonical :
      StrictPolynomialMod.SparsePolyZZCanonical remainder)
    (hdivisor : 0 < divisor.size) (hremainder : 0 < remainder.size)
    (hdivides : SparsePolyZZ.toPoly divisor ∣
      SparsePolyZZ.toPoly remainder) :
    divisor[0].2 ∣ remainder[0].2 := by
  have hdivisorLeading := sparsePolyZZ_leadingCoeff_eq_head divisor
    hdivisorCanonical hdivisor
  have hremainderLeading := sparsePolyZZ_leadingCoeff_eq_head remainder
    hremainderCanonical hremainder
  have hremainderNe : SparsePolyZZ.toPoly remainder ≠ 0 := by
    intro hzero
    rw [hzero] at hremainderLeading
    exact (hremainderCanonical.2 remainder[0]
      (Array.getElem_mem_toList hremainder)) (by
        simpa using hremainderLeading.symm)
  rcases hdivides with ⟨mathematicalQuotient, hfactor⟩
  have hquotientNe : mathematicalQuotient ≠ 0 := by
    intro hzero
    rw [hzero, mul_zero] at hfactor
    exact hremainderNe hfactor
  have hleading := congrArg Polynomial.leadingCoeff hfactor
  rw [Polynomial.leadingCoeff_mul] at hleading
  rw [hdivisorLeading, hremainderLeading] at hleading
  exact ⟨mathematicalQuotient.leadingCoeff, hleading⟩

/-- The checked decrease guard in one concrete exact-division iteration is
true on canonical inputs whenever the current leading coefficient is exactly
divisible.  This is the literal rank used by the generated well-founded
recursion, not a fuel counter. -/
theorem subtractScaledNormalize_divisionRank_lt
    (divisor remainder : SparsePolyZZ)
    (hdivisorCanonical :
      StrictPolynomialMod.SparsePolyZZCanonical divisor)
    (hremainderCanonical :
      StrictPolynomialMod.SparsePolyZZCanonical remainder)
    (hdivisor : 0 < divisor.size) (hremainder : 0 < remainder.size)
    (hdegree : divisor[0].1.deg ≤ remainder[0].1.deg)
    (hdivides : divisor[0].2 ∣ remainder[0].2) :
    Generated.StrictRecombine.divisionRank
        (Generated.StrictRecombine.subtractScaledNormalize remainder divisor
          (remainder[0].2 / divisor[0].2)
          (remainder[0].1.deg - divisor[0].1.deg)) <
      Generated.StrictRecombine.divisionRank remainder := by
  let degreeShift := remainder[0].1.deg - divisor[0].1.deg
  let scale := remainder[0].2 / divisor[0].2
  let next := Generated.StrictRecombine.subtractScaledNormalize remainder
    divisor scale degreeShift
  have hdivisorLeading := sparsePolyZZ_leadingCoeff_eq_head divisor
    hdivisorCanonical hdivisor
  have hremainderLeading := sparsePolyZZ_leadingCoeff_eq_head remainder
    hremainderCanonical hremainder
  have hdivisorNe : SparsePolyZZ.toPoly divisor ≠ 0 := by
    intro hzero
    rw [hzero] at hdivisorLeading
    exact (hdivisorCanonical.2 divisor[0]
      (Array.getElem_mem_toList hdivisor)) (by
        simpa using hdivisorLeading.symm)
  have hremainderNe : SparsePolyZZ.toPoly remainder ≠ 0 := by
    intro hzero
    rw [hzero] at hremainderLeading
    exact (hremainderCanonical.2 remainder[0]
      (Array.getElem_mem_toList hremainder)) (by
        simpa using hremainderLeading.symm)
  have hcancel : divisor[0].2 * scale = remainder[0].2 := by
    exact Int.mul_ediv_cancel_of_dvd hdivides
  have hscaleNe : scale ≠ 0 := by
    intro hzero
    rw [hzero, mul_zero] at hcancel
    exact (hremainderCanonical.2 remainder[0]
      (Array.getElem_mem_toList hremainder)) hcancel.symm
  let leadTerm := Polynomial.monomial degreeShift scale *
    SparsePolyZZ.toPoly divisor
  have hmonomialNe : Polynomial.monomial degreeShift scale ≠ 0 := by
    intro hzero
    have hcoeff := congrArg
      (fun poly : Polynomial Int => poly.coeff degreeShift) hzero
    simp [hscaleNe] at hcoeff
  have hleadTermNe : leadTerm ≠ 0 := by
    exact mul_ne_zero hmonomialNe hdivisorNe
  have hdivisorNatDegree := sparsePolyZZ_natDegree_eq_head divisor
    hdivisorCanonical hdivisor
  have hremainderNatDegree := sparsePolyZZ_natDegree_eq_head remainder
    hremainderCanonical hremainder
  have hleadTermNatDegree : leadTerm.natDegree =
      (SparsePolyZZ.toPoly remainder).natDegree := by
    rw [Polynomial.natDegree_mul
      hmonomialNe hdivisorNe,
      Polynomial.natDegree_monomial_eq degreeShift hscaleNe,
      hdivisorNatDegree, hremainderNatDegree]
    dsimp [degreeShift]
    omega
  have hleadTermLeading : leadTerm.leadingCoeff =
      (SparsePolyZZ.toPoly remainder).leadingCoeff := by
    rw [Polynomial.leadingCoeff_mul, Polynomial.leadingCoeff_monomial,
      hdivisorLeading, hremainderLeading]
    simpa [scale, mul_comm] using hcancel
  have hdegreeEq : (SparsePolyZZ.toPoly remainder).degree =
      leadTerm.degree := by
    rw [Polynomial.degree_eq_natDegree hremainderNe,
      Polynomial.degree_eq_natDegree hleadTermNe, hleadTermNatDegree]
  have hnextDegree : (SparsePolyZZ.toPoly next).degree <
      (SparsePolyZZ.toPoly remainder).degree := by
    rw [show SparsePolyZZ.toPoly next = SparsePolyZZ.toPoly remainder -
        leadTerm by
      exact subtractScaledNormalize_toPoly remainder divisor scale
        degreeShift]
    exact Polynomial.degree_sub_lt hdegreeEq hremainderNe
      hleadTermLeading.symm
  have hnextCanonical :
      StrictPolynomialMod.SparsePolyZZCanonical next :=
    normalization_canonical _
  unfold Generated.StrictRecombine.divisionRank
  rw [dif_pos hremainder]
  by_cases hnext : 0 < next.size
  · rw [dif_pos hnext]
    have hnextNe : SparsePolyZZ.toPoly next ≠ 0 := by
      intro hzero
      have hnextLeading := sparsePolyZZ_leadingCoeff_eq_head next
        hnextCanonical hnext
      rw [hzero] at hnextLeading
      exact (hnextCanonical.2 next[0] (Array.getElem_mem_toList hnext)) (by
        simpa using hnextLeading.symm)
    have hnatDegreeLt : (SparsePolyZZ.toPoly next).natDegree <
        (SparsePolyZZ.toPoly remainder).natDegree :=
      (Polynomial.natDegree_lt_natDegree_iff hnextNe).2 hnextDegree
    rw [sparsePolyZZ_natDegree_eq_head next hnextCanonical hnext,
      hremainderNatDegree] at hnatDegreeLt
    simpa [Nat.succ_eq_add_one] using Nat.succ_lt_succ hnatDegreeLt
  · rw [dif_neg hnext]
    exact Nat.zero_lt_succ _

/-- Completeness of the actual well-founded sparse exact-division loop on a
canonical divisible remainder.  Every branch condition is discharged from
the concrete representation and divisibility invariants, and recursion uses
the generated `divisionRank` decrease proved above. -/
theorem exactDivmodLoop_complete
    (divisor remainder quotient : SparsePolyZZ)
    (hdivisorCanonical :
      StrictPolynomialMod.SparsePolyZZCanonical divisor)
    (hremainderCanonical :
      StrictPolynomialMod.SparsePolyZZCanonical remainder)
    (hdivisor : 0 < divisor.size) :
    ∃ outputQuotient outputRemainder,
      Generated.StrictRecombine.exactDivmodLoop divisor remainder quotient =
        .ok (outputQuotient, outputRemainder) := by
  induction hmeasure : Generated.StrictRecombine.divisionRank remainder using
      Nat.strong_induction_on generalizing remainder quotient with
  | h measure ih =>
      rw [Generated.StrictRecombine.exactDivmodLoop]
      by_cases hremainder : 0 < remainder.size
      · rw [dif_pos hremainder, dif_pos hdivisor]
        dsimp only
        have hdivisorHeadNonzero : divisor[0].2 ≠ 0 :=
          hdivisorCanonical.2 divisor[0]
            (Array.getElem_mem_toList hdivisor)
        by_cases hdegree : divisor[0].1.deg ≤ remainder[0].1.deg
        · rw [dif_pos hdegree, dif_pos hdivisorHeadNonzero]
          by_cases hdivides : divisor[0].2 ∣ remainder[0].2
          · rw [dif_pos hdivides]
            let degreeShift := remainder[0].1.deg - divisor[0].1.deg
            let scale := remainder[0].2 / divisor[0].2
            let remainder' :=
              Generated.StrictRecombine.subtractScaledNormalize remainder
                divisor scale degreeShift
            let quotient' := quotient.push (⟨degreeShift⟩, scale)
            have hdecrease : Generated.StrictRecombine.divisionRank remainder' <
                Generated.StrictRecombine.divisionRank remainder := by
              exact subtractScaledNormalize_divisionRank_lt divisor remainder
                hdivisorCanonical hremainderCanonical hdivisor hremainder
                hdegree hdivides
            rw [dif_pos hdecrease]
            exact ih (Generated.StrictRecombine.divisionRank remainder')
              (by simpa [hmeasure] using hdecrease) remainder' quotient'
              (normalization_canonical _) rfl
          · rw [dif_neg hdivides]
            exact ⟨quotient, remainder, rfl⟩
        · rw [dif_neg hdegree]
          exact ⟨quotient, remainder, rfl⟩
      · rw [dif_neg hremainder]
        exact ⟨quotient, remainder, rfl⟩

/-- The public checked long-division entry never faults on a canonical
dividend and a canonical nonempty divisor.  Nondivisibility is represented by
the physical remainder, exactly as in C++, rather than by an error. -/
theorem exactDivmodRaw_complete
    (dividend divisor : SparsePolyZZ)
    (hdividendCanonical :
      StrictPolynomialMod.SparsePolyZZCanonical dividend)
    (hdivisorCanonical :
      StrictPolynomialMod.SparsePolyZZCanonical divisor)
    (hdivisor : 0 < divisor.size) :
    ∃ quotient remainder,
      Generated.StrictRecombine.exactDivmodRaw dividend divisor =
        .ok (quotient, remainder) := by
  unfold Generated.StrictRecombine.exactDivmodRaw
  exact exactDivmodLoop_complete divisor dividend #[] hdivisorCanonical
    hdividendCanonical hdivisor

theorem exactDivmodLoop_complete_of_dvd
    (divisor remainder quotient : SparsePolyZZ)
    (hdivisorCanonical :
      StrictPolynomialMod.SparsePolyZZCanonical divisor)
    (hremainderCanonical :
      StrictPolynomialMod.SparsePolyZZCanonical remainder)
    (hdivisor : 0 < divisor.size)
    (hdivides : SparsePolyZZ.toPoly divisor ∣
      SparsePolyZZ.toPoly remainder) :
    ∃ outputQuotient,
      Generated.StrictRecombine.exactDivmodLoop divisor remainder quotient =
        .ok (outputQuotient, #[]) := by
  induction hmeasure : Generated.StrictRecombine.divisionRank remainder using
      Nat.strong_induction_on generalizing remainder quotient with
  | h measure ih =>
      rw [Generated.StrictRecombine.exactDivmodLoop]
      by_cases hremainder : 0 < remainder.size
      · rw [dif_pos hremainder, dif_pos hdivisor]
        dsimp only
        have hdivisorHeadNonzero : divisor[0].2 ≠ 0 :=
          hdivisorCanonical.2 divisor[0]
            (Array.getElem_mem_toList hdivisor)
        have hremainderNe : SparsePolyZZ.toPoly remainder ≠ 0 := by
          intro hzero
          have hleading := sparsePolyZZ_leadingCoeff_eq_head remainder
            hremainderCanonical hremainder
          rw [hzero] at hleading
          exact (hremainderCanonical.2 remainder[0]
            (Array.getElem_mem_toList hremainder)) (by
              simpa using hleading.symm)
        have hdegreePoly : (SparsePolyZZ.toPoly divisor).natDegree ≤
            (SparsePolyZZ.toPoly remainder).natDegree :=
          Polynomial.natDegree_le_of_dvd hdivides hremainderNe
        have hdegree : divisor[0].1.deg ≤ remainder[0].1.deg := by
          rw [sparsePolyZZ_natDegree_eq_head divisor hdivisorCanonical
              hdivisor,
            sparsePolyZZ_natDegree_eq_head remainder hremainderCanonical
              hremainder] at hdegreePoly
          exact hdegreePoly
        rw [dif_pos hdegree, dif_pos hdivisorHeadNonzero]
        have hcoeffDivides := canonical_head_coefficient_dvd_of_poly_dvd
          divisor remainder hdivisorCanonical hremainderCanonical hdivisor
          hremainder hdivides
        rw [dif_pos hcoeffDivides]
        let degreeShift := remainder[0].1.deg - divisor[0].1.deg
        let scale := remainder[0].2 / divisor[0].2
        let remainder' :=
          Generated.StrictRecombine.subtractScaledNormalize remainder divisor
            scale degreeShift
        let quotient' := quotient.push (⟨degreeShift⟩, scale)
        have hdecrease : Generated.StrictRecombine.divisionRank remainder' <
            Generated.StrictRecombine.divisionRank remainder := by
          exact subtractScaledNormalize_divisionRank_lt divisor remainder
            hdivisorCanonical hremainderCanonical hdivisor hremainder hdegree
            hcoeffDivides
        rw [dif_pos hdecrease]
        have hremainder'Canonical :
            StrictPolynomialMod.SparsePolyZZCanonical remainder' :=
          normalization_canonical _
        have hremainder'Divides : SparsePolyZZ.toPoly divisor ∣
            SparsePolyZZ.toPoly remainder' := by
          rcases hdivides with ⟨mathematicalQuotient, hfactor⟩
          refine ⟨mathematicalQuotient -
            Polynomial.monomial degreeShift scale, ?_⟩
          rw [show SparsePolyZZ.toPoly remainder' =
              SparsePolyZZ.toPoly remainder -
                Polynomial.monomial degreeShift scale *
                  SparsePolyZZ.toPoly divisor by
            exact subtractScaledNormalize_toPoly remainder divisor scale
              degreeShift]
          rw [hfactor]
          ring
        exact ih (Generated.StrictRecombine.divisionRank remainder')
          (by simpa [hmeasure] using hdecrease) remainder' quotient'
          hremainder'Canonical hremainder'Divides rfl
      · rw [dif_neg hremainder]
        have hempty : remainder = #[] :=
          Array.size_eq_zero_iff.mp (Nat.eq_zero_of_not_pos hremainder)
        subst remainder
        exact ⟨quotient, rfl⟩

/-- Public exact-division completeness theorem for the literal generated
entry.  A canonical nonempty divisor that truly divides a canonical dividend
causes the concrete recursion to return an actual quotient and the physical
empty sparse remainder. -/
theorem exactDivmodRaw_complete_of_dvd
    (dividend divisor : SparsePolyZZ)
    (hdividendCanonical :
      StrictPolynomialMod.SparsePolyZZCanonical dividend)
    (hdivisorCanonical :
      StrictPolynomialMod.SparsePolyZZCanonical divisor)
    (hdivisor : 0 < divisor.size)
    (hdivides : SparsePolyZZ.toPoly divisor ∣
      SparsePolyZZ.toPoly dividend) :
    ∃ quotient,
      Generated.StrictRecombine.exactDivmodRaw dividend divisor =
        .ok (quotient, #[]) := by
  unfold Generated.StrictRecombine.exactDivmodRaw
  exact exactDivmodLoop_complete_of_dvd divisor dividend #[]
    hdivisorCanonical hdividendCanonical hdivisor hdivides

theorem exactDivmodLoop_toPoly (divisor remainder quotient q r : SparsePolyZZ)
    (hrun : Generated.StrictRecombine.exactDivmodLoop divisor remainder quotient =
      .ok (q, r)) :
    SparsePolyZZ.toPoly divisor * SparsePolyZZ.toPoly q +
        SparsePolyZZ.toPoly r =
      SparsePolyZZ.toPoly divisor * SparsePolyZZ.toPoly quotient +
        SparsePolyZZ.toPoly remainder := by
  induction hmeasure : Generated.StrictRecombine.divisionRank remainder using
      Nat.strong_induction_on generalizing remainder quotient q r with
  | h measure ih =>
      rw [Generated.StrictRecombine.exactDivmodLoop] at hrun
      split at hrun
      next hremainder =>
        split at hrun
        next hdivisor =>
          dsimp at hrun
          split at hrun
          next hdegree =>
            split at hrun
            next hnonzero =>
              split at hrun
              next hdivides =>
                split at hrun
                next hdecrease =>
                  let degreeShift := remainder[0].1.deg - divisor[0].1.deg
                  let scale := remainder[0].2 / divisor[0].2
                  let remainder' :=
                    Generated.StrictRecombine.subtractScaledNormalize remainder
                      divisor scale degreeShift
                  let quotient' := quotient.push (⟨degreeShift⟩, scale)
                  have hdec : Generated.StrictRecombine.divisionRank remainder' <
                      measure := by
                    rw [← hmeasure]
                    exact hdecrease
                  have htail := ih
                    (Generated.StrictRecombine.divisionRank remainder')
                    hdec remainder' quotient' q r hrun rfl
                  rw [htail, subtractScaledNormalize_toPoly]
                  simp only [quotient', SparsePolyZZ.toPoly, Array.toList_push,
                    List.map_append, List.map_singleton, List.sum_append,
                    List.sum_singleton]
                  change _ * (_ + Polynomial.monomial degreeShift scale) +
                      (_ - Polynomial.monomial degreeShift scale * _) = _
                  ring
                next hdecrease => contradiction
              next hdivides =>
                have hout := Except.ok.inj hrun
                cases hout
                rfl
            next hnonzero => contradiction
          next hdegree =>
            have hout := Except.ok.inj hrun
            cases hout
            rfl
        next hdivisor => contradiction
      next hremainder =>
        have hout := Except.ok.inj hrun
        cases hout
        rfl

/-- Recursive remainders used and returned by the checked exact-division loop
have nonzero stored coefficients.  Every recursive step obtains the property
from the concrete sparse normalization. -/
theorem exactDivmodLoop_remainder_coefficients_nonzero
    (divisor remainder quotient q r : SparsePolyZZ)
    (hremainderNonzero : ∀ term ∈ remainder.toList, term.2 ≠ 0)
    (hrun : Generated.StrictRecombine.exactDivmodLoop divisor remainder quotient =
      .ok (q, r)) :
    ∀ term ∈ r.toList, term.2 ≠ 0 := by
  induction hmeasure : Generated.StrictRecombine.divisionRank remainder using
      Nat.strong_induction_on generalizing remainder quotient q r with
  | h measure ih =>
      rw [Generated.StrictRecombine.exactDivmodLoop] at hrun
      split at hrun
      next hremainder =>
        split at hrun
        next hdivisor =>
          dsimp at hrun
          split at hrun
          next hdegree =>
            split at hrun
            next hnonzero =>
              split at hrun
              next hdivides =>
                split at hrun
                next hdecrease =>
                  let degreeShift := remainder[0].1.deg - divisor[0].1.deg
                  let scale := remainder[0].2 / divisor[0].2
                  let remainder' :=
                    Generated.StrictRecombine.subtractScaledNormalize remainder
                      divisor scale degreeShift
                  have hdec : Generated.StrictRecombine.divisionRank remainder' <
                      measure := by
                    rw [← hmeasure]
                    simpa [remainder', scale, degreeShift] using hdecrease
                  exact ih (Generated.StrictRecombine.divisionRank remainder')
                    hdec remainder'
                    (quotient.push (⟨degreeShift⟩, scale)) q r
                    (normalization_coefficients_nonzero
                      (Generated.StrictRecombine.subtractScaledTermsLoop divisor
                        scale degreeShift 0 remainder)) hrun rfl
                next hdecrease => contradiction
              next hdivides =>
                have hout := Except.ok.inj hrun
                injection hout with hq hr
                subst r
                exact hremainderNonzero
            next hnonzero => contradiction
          next hdegree =>
            have hout := Except.ok.inj hrun
            injection hout with hq hr
            subst r
            exact hremainderNonzero
        next hdivisor => contradiction
      next hremainder =>
        have hout := Except.ok.inj hrun
        injection hout with hq hr
        subst r
        exact hremainderNonzero

/-- Accumulator invariant for the source long-division quotient.  Existing
terms are canonical and lie strictly above the next degree shift whenever the
current remainder can take another division step. -/
def ExactDivmodQuotientInvariant (divisor remainder quotient : SparsePolyZZ) :
    Prop :=
  StrictPolynomialMod.SparsePolyZZCanonical quotient ∧
  ∀ (hremainder : 0 < remainder.size) (hdivisor : 0 < divisor.size),
    divisor[0].1.deg ≤ remainder[0].1.deg →
    ∀ term ∈ quotient.toList,
      remainder[0].1.deg - divisor[0].1.deg < term.1.deg

private theorem sparsePolyZZCanonical_push
    (quotient : SparsePolyZZ) (term : UMonomial × Int)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical quotient)
    (habove : ∀ existing ∈ quotient.toList,
      term.1.deg < existing.1.deg)
    (hnonzero : term.2 ≠ 0) :
    StrictPolynomialMod.SparsePolyZZCanonical (quotient.push term) := by
  unfold StrictPolynomialMod.SparsePolyZZCanonical at hcanonical ⊢
  rw [Array.toList_push]
  constructor
  · rw [List.isChain_append]
    refine ⟨hcanonical.1, by simp, ?_⟩
    intro existing hexisting candidate hcandidate
    have hcandidateEq : candidate = term := by simpa using hcandidate.symm
    subst candidate
    exact habove existing (List.mem_of_mem_getLast? hexisting)
  · intro candidate hcandidate
    rcases List.mem_append.mp hcandidate with hprefix | hlast
    · exact hcanonical.2 candidate hprefix
    · have : candidate = term := by simpa using hlast
      subst candidate
      exact hnonzero

theorem exactDivmodQuotientInvariant_empty
    (divisor remainder : SparsePolyZZ) :
    ExactDivmodQuotientInvariant divisor remainder #[] := by
  constructor
  · simp [StrictPolynomialMod.SparsePolyZZCanonical]
  · simp

/-- The quotient returned by the actual checked long-division recursion is a
canonical sparse polynomial.  Its pushed degree shifts strictly decrease by
the source `divisionRank` check, and its pushed coefficients are nonzero by
exact division of the nonzero current leading coefficient. -/
theorem exactDivmodLoop_quotient_canonical
    (divisor remainder quotient q r : SparsePolyZZ)
    (hremainderNonzero : ∀ term ∈ remainder.toList, term.2 ≠ 0)
    (hinvariant : ExactDivmodQuotientInvariant divisor remainder quotient)
    (hrun : Generated.StrictRecombine.exactDivmodLoop divisor remainder quotient =
      .ok (q, r)) :
    StrictPolynomialMod.SparsePolyZZCanonical q := by
  induction hmeasure : Generated.StrictRecombine.divisionRank remainder using
      Nat.strong_induction_on generalizing remainder quotient q r with
  | h measure ih =>
      rw [Generated.StrictRecombine.exactDivmodLoop] at hrun
      split at hrun
      next hremainder =>
        split at hrun
        next hdivisor =>
          dsimp at hrun
          split at hrun
          next hdegree =>
            split at hrun
            next hnonzero =>
              split at hrun
              next hdivides =>
                split at hrun
                next hdecrease =>
                  let degreeShift := remainder[0].1.deg - divisor[0].1.deg
                  let scale := remainder[0].2 / divisor[0].2
                  let remainder' :=
                    Generated.StrictRecombine.subtractScaledNormalize remainder
                      divisor scale degreeShift
                  let quotient' := quotient.push (⟨degreeShift⟩, scale)
                  have hscaleNonzero : scale ≠ 0 := by
                    have hleadNonzero := hremainderNonzero remainder[0]
                      (Array.getElem_mem_toList hremainder)
                    intro hscale
                    have hcancel : divisor[0].2 *
                        (remainder[0].2 / divisor[0].2) = remainder[0].2 :=
                      Int.mul_ediv_cancel_of_dvd hdivides
                    change remainder[0].2 / divisor[0].2 = 0 at hscale
                    rw [hscale] at hcancel
                    simp at hcancel
                    exact hleadNonzero hcancel.symm
                  have hquotientCanonical :
                      StrictPolynomialMod.SparsePolyZZCanonical quotient' := by
                    apply sparsePolyZZCanonical_push quotient
                      (⟨degreeShift⟩, scale) hinvariant.1
                    · intro existing hexisting
                      exact hinvariant.2 hremainder hdivisor hdegree existing
                        hexisting
                    · exact hscaleNonzero
                  have hremainder'Nonzero : ∀ term ∈ remainder'.toList,
                      term.2 ≠ 0 :=
                    normalization_coefficients_nonzero
                      (Generated.StrictRecombine.subtractScaledTermsLoop divisor
                        scale degreeShift 0 remainder)
                  have hinvariant' : ExactDivmodQuotientInvariant divisor
                      remainder' quotient' := by
                    refine ⟨hquotientCanonical, ?_⟩
                    intro hremainder' hdivisor' hdegree' term hterm
                    have hfrontDecrease : remainder'[0].1.deg <
                        remainder[0].1.deg := by
                      change Generated.StrictRecombine.divisionRank remainder' <
                        Generated.StrictRecombine.divisionRank remainder at hdecrease
                      rw [Generated.StrictRecombine.divisionRank,
                        dif_pos hremainder',
                        Generated.StrictRecombine.divisionRank,
                        dif_pos hremainder] at hdecrease
                      omega
                    have hshiftDecrease :
                        remainder'[0].1.deg - divisor[0].1.deg < degreeShift := by
                      dsimp [degreeShift]
                      omega
                    dsimp [quotient'] at hterm
                    rw [Array.toList_push] at hterm
                    rcases List.mem_append.mp hterm with hprefix | hlast
                    · exact lt_trans hshiftDecrease
                        (hinvariant.2 hremainder hdivisor hdegree term hprefix)
                    · have htermEq : term = (⟨degreeShift⟩, scale) := by
                        simpa using hlast
                      subst term
                      exact hshiftDecrease
                  have hdec : Generated.StrictRecombine.divisionRank remainder' <
                      measure := by
                    rw [← hmeasure]
                    simpa [remainder', scale, degreeShift] using hdecrease
                  exact ih (Generated.StrictRecombine.divisionRank remainder')
                    hdec remainder' quotient' q r hremainder'Nonzero hinvariant'
                    hrun rfl
                next hdecrease => contradiction
              next hdivides =>
                have hout := Except.ok.inj hrun
                injection hout with hq hr
                subst q
                exact hinvariant.1
            next hnonzero => contradiction
          next hdegree =>
            have hout := Except.ok.inj hrun
            injection hout with hq hr
            subst q
            exact hinvariant.1
        next hdivisor => contradiction
      next hremainder =>
        have hout := Except.ok.inj hrun
        injection hout with hq hr
        subst q
        exact hinvariant.1

theorem exactDivmodRaw_toPoly (dividend divisor quotient remainder : SparsePolyZZ)
    (hrun : Generated.StrictRecombine.exactDivmodRaw dividend divisor =
      .ok (quotient, remainder)) :
    SparsePolyZZ.toPoly dividend =
      SparsePolyZZ.toPoly divisor * SparsePolyZZ.toPoly quotient +
        SparsePolyZZ.toPoly remainder := by
  unfold Generated.StrictRecombine.exactDivmodRaw at hrun
  have hsemantic := exactDivmodLoop_toPoly divisor dividend #[] quotient
    remainder hrun
  simpa [SparsePolyZZ.toPoly] using hsemantic.symm

/-- Algebraic closure step for completeness of the actual exact-division
loop.  If its concrete canonical remainder is below the divisor degree, then
a genuine divisibility hypothesis forces that remainder array to be empty.
The quotient and remainder are the values returned by `exactDivmodRaw`; no
mathematical quotient is substituted for the generated execution. -/
theorem exactDivmodRaw_remainder_eq_empty_of_dvd_of_degree_lt
    (dividend divisor quotient remainder : SparsePolyZZ)
    (hdivisorNe : SparsePolyZZ.toPoly divisor ≠ 0)
    (hremainderCanonical :
      StrictPolynomialMod.SparsePolyZZCanonical remainder)
    (hdegree : SparsePolyZZ.toPoly remainder = 0 ∨
      (SparsePolyZZ.toPoly remainder).natDegree <
        (SparsePolyZZ.toPoly divisor).natDegree)
    (hdivides : SparsePolyZZ.toPoly divisor ∣
      SparsePolyZZ.toPoly dividend)
    (hrun : Generated.StrictRecombine.exactDivmodRaw dividend divisor =
      .ok (quotient, remainder)) :
    remainder = #[] := by
  have hsemantic := exactDivmodRaw_toPoly dividend divisor quotient remainder
    hrun
  rcases hdivides with ⟨mathematicalQuotient, hdividend⟩
  have hremainderDivides : SparsePolyZZ.toPoly divisor ∣
      SparsePolyZZ.toPoly remainder := by
    refine ⟨mathematicalQuotient - SparsePolyZZ.toPoly quotient, ?_⟩
    calc
      SparsePolyZZ.toPoly remainder =
          SparsePolyZZ.toPoly dividend -
            SparsePolyZZ.toPoly divisor * SparsePolyZZ.toPoly quotient := by
        rw [hsemantic]
        ring
      _ = SparsePolyZZ.toPoly divisor *
          (mathematicalQuotient - SparsePolyZZ.toPoly quotient) := by
        rw [hdividend]
        ring
  have hremainderZero : SparsePolyZZ.toPoly remainder = 0 := by
    rcases hdegree with hzero | hdegree
    · exact hzero
    · by_contra hnonzero
      have hdegreeLe := Polynomial.natDegree_le_of_dvd hremainderDivides
        hnonzero
      omega
  have hsizeZero : remainder.size = 0 := by
    by_contra hnot
    have hnonempty : 0 < remainder.size := Nat.pos_of_ne_zero hnot
    have hleading := sparsePolyZZ_leadingCoeff_eq_head remainder
      hremainderCanonical hnonempty
    rw [hremainderZero] at hleading
    have hheadNonzero := hremainderCanonical.2 remainder[0]
      (Array.getElem_mem_toList hnonempty)
    exact hheadNonzero (by simpa using hleading.symm)
  exact Array.size_eq_zero_iff.mp hsizeZero

theorem exactDivmodRaw_quotient_canonical
    (dividend divisor quotient remainder : SparsePolyZZ)
    (hdividendNonzero : ∀ term ∈ dividend.toList, term.2 ≠ 0)
    (hrun : Generated.StrictRecombine.exactDivmodRaw dividend divisor =
      .ok (quotient, remainder)) :
    StrictPolynomialMod.SparsePolyZZCanonical quotient := by
  unfold Generated.StrictRecombine.exactDivmodRaw at hrun
  exact exactDivmodLoop_quotient_canonical divisor dividend #[] quotient
    remainder hdividendNonzero
    (exactDivmodQuotientInvariant_empty divisor dividend) hrun

theorem successfulTrialExtraction_toPoly
    (fStar factor quotient quotientPrimitive : SparsePolyZZ)
    (quotientContent : Int)
    (hdivide : Generated.StrictRecombine.exactDivmodRaw fStar factor =
      .ok (quotient, #[]))
    (hprimitive : Generated.StrictRecombine.primitiveRaw quotient =
      .ok (quotientContent, quotientPrimitive)) :
    SparsePolyZZ.toPoly fStar =
      Polynomial.C quotientContent *
        (SparsePolyZZ.toPoly factor * SparsePolyZZ.toPoly quotientPrimitive) := by
  rw [exactDivmodRaw_toPoly fStar factor quotient #[] hdivide]
  rw [primitiveRaw_toPoly quotient quotientPrimitive quotientContent hprimitive]
  simp [SparsePolyZZ.toPoly]
  ring

/-- A successful concrete Zassenhaus attempt exposes the exact successful
long-division and quotient-primitive executions that produced its returned
quotient. -/
theorem zassenhausAttempt_extracted_quotient_trace
    (fStar factor quotientPrimitive : SparsePolyZZ)
    (activeLifted : Array SparsePolyZZ) (modulus : ZZ)
    (candidate : Array Nat)
    (hrun : Generated.StrictRecombine.zassenhausAttempt fStar activeLifted
      modulus candidate = .ok (.extracted factor quotientPrimitive)) :
    ∃ quotient quotientContent,
      Generated.StrictRecombine.exactDivmodRaw fStar factor =
        .ok (quotient, #[]) ∧
      Generated.StrictRecombine.primitiveRaw quotient =
        .ok (quotientContent, quotientPrimitive) := by
  unfold Generated.StrictRecombine.zassenhausAttempt at hrun
  split at hrun
  next hfstar =>
    dsimp at hrun
    cases hleading : Generated.StrictRecombine.selectedLeadingProductLoop
        candidate activeLifted 0 fStar[0].2 with
    | error fault => simp [hleading] at hrun
    | ok leadingProduct =>
      simp only [hleading] at hrun
      split at hrun
      next hpruned => simp at hrun
      next hleadingAccepted =>
        cases hconstant : Generated.StrictRecombine.selectedConstantProductLoop
            candidate activeLifted 0 fStar[0].2 with
        | error fault => simp [hconstant] at hrun
        | ok constantProduct =>
          simp only [hconstant] at hrun
          split at hrun
          next hpruned => simp at hrun
          next hconstantAccepted =>
            cases hconvert : Generated.StrictRecombine.combinationToInt32
                candidate with
            | error fault => simp [hconvert] at hrun
            | ok candidate32 =>
              simp only [hconvert] at hrun
              cases hproduct : Generated.StrictRecombine.trialProductLoop
                  ⟨()⟩ candidate32 activeLifted modulus 0
                  #[(⟨0⟩, fStar[0].2)] with
              | error fault => simp [hproduct] at hrun
              | ok product =>
                simp only [hproduct] at hrun
                cases hsymmetric : Generated.StrictRecombine.symmetricModRaw
                    product modulus with
                | error fault => simp [hsymmetric] at hrun
                | ok symmetric =>
                  simp only [hsymmetric] at hrun
                  cases hprimitive : Generated.StrictRecombine.primitiveRaw
                      symmetric with
                  | error fault => simp [hprimitive] at hrun
                  | ok primitiveResult =>
                    rcases primitiveResult with ⟨symmetricContent,
                      recoveredFactor⟩
                    simp only [hprimitive] at hrun
                    cases hdivmod : Generated.StrictRecombine.exactDivmodRaw
                        fStar recoveredFactor with
                    | error fault => simp [hdivmod] at hrun
                    | ok divResult =>
                      rcases divResult with ⟨quotient, remainder⟩
                      simp only [hdivmod] at hrun
                      by_cases hremainder : remainder.isEmpty = true
                      · simp only [hremainder, if_true] at hrun
                        have hremainderEmpty : remainder = #[] :=
                          Array.isEmpty_iff.mp hremainder
                        subst remainder
                        cases hquotientPrimitive :
                            Generated.StrictRecombine.primitiveRaw quotient with
                        | error fault => simp [hquotientPrimitive] at hrun
                        | ok quotientResult =>
                          rcases quotientResult with ⟨quotientContent,
                            recoveredQuotient⟩
                          simp only [hquotientPrimitive] at hrun
                          have hout := Except.ok.inj hrun
                          injection hout with hfactor hquotient
                          subst factor
                          subst quotientPrimitive
                          exact ⟨quotient, quotientContent, hdivmod,
                            hquotientPrimitive⟩
                      · simp only [hremainder, if_false] at hrun
                        simp at hrun
  next hfstar => contradiction

/-- A successful attempt exposes the exact candidate lowering, modular trial
product, symmetric representative, and primitive extraction that produced its
returned factor.  These witnesses are intermediate values from the generated
execution itself. -/
theorem zassenhausAttempt_extracted_candidate_trace
    (fStar factor quotientPrimitive : SparsePolyZZ)
    (activeLifted : Array SparsePolyZZ) (modulus : ZZ)
    (candidate : Array Nat)
    (hrun : Generated.StrictRecombine.zassenhausAttempt fStar activeLifted
      modulus candidate = .ok (.extracted factor quotientPrimitive)) :
    ∃ candidate32 product symmetric symmetricContent,
      Generated.StrictRecombine.combinationToInt32 candidate =
        .ok candidate32 ∧
      Generated.StrictRecombine.trialProductLoop ⟨()⟩ candidate32
        activeLifted modulus 0 #[(⟨0⟩, fStar[0]!.2)] = .ok product ∧
      Generated.StrictRecombine.symmetricModRaw product modulus =
        .ok symmetric ∧
      Generated.StrictRecombine.primitiveRaw symmetric =
        .ok (symmetricContent, factor) := by
  unfold Generated.StrictRecombine.zassenhausAttempt at hrun
  split at hrun
  next hfstar =>
    dsimp at hrun
    cases hleading : Generated.StrictRecombine.selectedLeadingProductLoop
        candidate activeLifted 0 fStar[0].2 with
    | error fault => simp [hleading] at hrun
    | ok leadingProduct =>
      simp only [hleading] at hrun
      split at hrun
      next hpruned => simp at hrun
      next hleadingAccepted =>
        cases hconstant : Generated.StrictRecombine.selectedConstantProductLoop
            candidate activeLifted 0 fStar[0].2 with
        | error fault => simp [hconstant] at hrun
        | ok constantProduct =>
          simp only [hconstant] at hrun
          split at hrun
          next hpruned => simp at hrun
          next hconstantAccepted =>
            cases hconvert : Generated.StrictRecombine.combinationToInt32
                candidate with
            | error fault => simp [hconvert] at hrun
            | ok candidate32 =>
              simp only [hconvert] at hrun
              cases hproduct : Generated.StrictRecombine.trialProductLoop
                  ⟨()⟩ candidate32 activeLifted modulus 0
                  #[(⟨0⟩, fStar[0].2)] with
              | error fault => simp [hproduct] at hrun
              | ok product =>
                simp only [hproduct] at hrun
                cases hsymmetric : Generated.StrictRecombine.symmetricModRaw
                    product modulus with
                | error fault => simp [hsymmetric] at hrun
                | ok symmetric =>
                  simp only [hsymmetric] at hrun
                  cases hprimitive : Generated.StrictRecombine.primitiveRaw
                      symmetric with
                  | error fault => simp [hprimitive] at hrun
                  | ok primitiveResult =>
                    rcases primitiveResult with
                      ⟨symmetricContent, recoveredFactor⟩
                    simp only [hprimitive] at hrun
                    cases hdivmod : Generated.StrictRecombine.exactDivmodRaw
                        fStar recoveredFactor with
                    | error fault => simp [hdivmod] at hrun
                    | ok divResult =>
                      rcases divResult with ⟨quotient, remainder⟩
                      simp only [hdivmod] at hrun
                      by_cases hremainder : remainder.isEmpty = true
                      · simp only [hremainder, if_true] at hrun
                        cases hquotientPrimitive :
                            Generated.StrictRecombine.primitiveRaw quotient with
                        | error fault => simp [hquotientPrimitive] at hrun
                        | ok quotientResult =>
                          rcases quotientResult with
                            ⟨quotientContent, recoveredQuotient⟩
                          simp only [hquotientPrimitive] at hrun
                          have hout := Except.ok.inj hrun
                          injection hout with hfactor hquotient
                          subst factor
                          refine ⟨candidate32, product, symmetric,
                            symmetricContent, rfl, ?_, hsymmetric,
                            hprimitive⟩
                          simpa [getElem!_pos fStar 0 hfstar] using hproduct
                      · simp only [hremainder, if_false] at hrun
                        simp at hrun
  next hfstar => contradiction

/-- Modulo any positive divisor of the actual Hensel modulus, the factor
returned by the generated successful attempt is related by its concrete
primitive-content scalar to the exact selected lifted-factor subproduct. -/
theorem zassenhausAttempt_extracted_factor_mod_eq_selected
    (fStar factor quotientPrimitive : SparsePolyZZ)
    (activeLifted : Array SparsePolyZZ) (modulus base : Nat)
    (candidate : Array Nat)
    (hmodulus : 0 < modulus) (hbase : 0 < base)
    (hdivides : base ∣ modulus)
    (hbound : ∀ position (hposition : position < candidate.size),
      candidate[position] < activeLifted.size)
    (hactiveFits : activeLifted.size ≤ 2 ^ 31)
    (hrun : Generated.StrictRecombine.zassenhausAttempt fStar activeLifted
      (modulus : ZZ) candidate =
        .ok (.extracted factor quotientPrimitive)) :
    ∃ content : Int,
      Polynomial.C (content : ZMod base) *
          Refinement.StrictHensel.toPolyMod base factor =
        Polynomial.C (fStar[0]!.2 : ZMod base) *
          ((selectSourceIndices activeLifted.toList candidate.toList).map
            (Refinement.StrictHensel.toPolyMod base)).prod := by
  rcases zassenhausAttempt_extracted_candidate_trace fStar factor
      quotientPrimitive activeLifted (modulus : ZZ) candidate hrun with
    ⟨candidate32, product, symmetric, content, hconvert, hproduct,
      hsymmetric, hprimitive⟩
  refine ⟨content, ?_⟩
  have htrial := trialProductLoop_source_indices_refines_of_dvd modulus base
    hmodulus hbase hdivides candidate activeLifted candidate32
    #[(⟨0⟩, fStar[0]!.2)] product hbound hactiveFits hconvert hproduct
  have hsymmetricMod := symmetricModRaw_toPolyMod_of_dvd product symmetric
    modulus base hmodulus hbase hdivides hsymmetric
  have hprimitivePoly := primitiveRaw_toPoly symmetric factor content hprimitive
  have hprimitiveMod :=
    congrArg (Polynomial.map (Int.castRingHom (ZMod base))) hprimitivePoly
  have hprimitiveMod' : Refinement.StrictHensel.toPolyMod base symmetric =
      Polynomial.C (content : ZMod base) *
        Refinement.StrictHensel.toPolyMod base factor := by
    simpa [Refinement.StrictHensel.toPolyMod] using hprimitiveMod
  have hinitial : Refinement.StrictHensel.toPolyMod base
      #[(⟨0⟩, fStar[0]!.2)] =
        Polynomial.C (fStar[0]!.2 : ZMod base) := by
    simp [Refinement.StrictHensel.toPolyMod, SparsePolyZZ.toPoly]
  rw [hsymmetricMod, htrial, hinitial] at hprimitiveMod'
  exact hprimitiveMod'.symm

/-- The same physical primitive-content scalar certifies a successful
candidate simultaneously at two divisors of the actual Hensel modulus.  Both
equations are derived from one execution trace, so their existential witnesses
cannot drift apart. -/
theorem zassenhausAttempt_extracted_factor_mod_eq_selected_pair
    (fStar factor quotientPrimitive : SparsePolyZZ)
    (activeLifted : Array SparsePolyZZ) (modulus base₁ base₂ : Nat)
    (candidate : Array Nat)
    (hmodulus : 0 < modulus) (hbase₁ : 0 < base₁) (hbase₂ : 0 < base₂)
    (hdivides₁ : base₁ ∣ modulus) (hdivides₂ : base₂ ∣ modulus)
    (hbound : ∀ position (hposition : position < candidate.size),
      candidate[position] < activeLifted.size)
    (hactiveFits : activeLifted.size ≤ 2 ^ 31)
    (hrun : Generated.StrictRecombine.zassenhausAttempt fStar activeLifted
      (modulus : ZZ) candidate =
        .ok (.extracted factor quotientPrimitive)) :
    ∃ content : Int,
      (Polynomial.C (content : ZMod base₁) *
          Refinement.StrictHensel.toPolyMod base₁ factor =
        Polynomial.C (fStar[0]!.2 : ZMod base₁) *
          ((selectSourceIndices activeLifted.toList candidate.toList).map
            (Refinement.StrictHensel.toPolyMod base₁)).prod) ∧
      (Polynomial.C (content : ZMod base₂) *
          Refinement.StrictHensel.toPolyMod base₂ factor =
        Polynomial.C (fStar[0]!.2 : ZMod base₂) *
          ((selectSourceIndices activeLifted.toList candidate.toList).map
            (Refinement.StrictHensel.toPolyMod base₂)).prod) := by
  rcases zassenhausAttempt_extracted_candidate_trace fStar factor
      quotientPrimitive activeLifted (modulus : ZZ) candidate hrun with
    ⟨candidate32, product, symmetric, content, hconvert, hproduct,
      hsymmetric, hprimitive⟩
  have equation (base : Nat) (hbase : 0 < base)
      (hdivides : base ∣ modulus) :
      Polynomial.C (content : ZMod base) *
          Refinement.StrictHensel.toPolyMod base factor =
        Polynomial.C (fStar[0]!.2 : ZMod base) *
          ((selectSourceIndices activeLifted.toList candidate.toList).map
            (Refinement.StrictHensel.toPolyMod base)).prod := by
    have htrial := trialProductLoop_source_indices_refines_of_dvd modulus base
      hmodulus hbase hdivides candidate activeLifted candidate32
      #[(⟨0⟩, fStar[0]!.2)] product hbound hactiveFits hconvert hproduct
    have hsymmetricMod := symmetricModRaw_toPolyMod_of_dvd product symmetric
      modulus base hmodulus hbase hdivides hsymmetric
    have hprimitivePoly := primitiveRaw_toPoly symmetric factor content
      hprimitive
    have hprimitiveMod :=
      congrArg (Polynomial.map (Int.castRingHom (ZMod base))) hprimitivePoly
    have hprimitiveMod' : Refinement.StrictHensel.toPolyMod base symmetric =
        Polynomial.C (content : ZMod base) *
          Refinement.StrictHensel.toPolyMod base factor := by
      simpa [Refinement.StrictHensel.toPolyMod] using hprimitiveMod
    have hinitial : Refinement.StrictHensel.toPolyMod base
        #[(⟨0⟩, fStar[0]!.2)] =
          Polynomial.C (fStar[0]!.2 : ZMod base) := by
      simp [Refinement.StrictHensel.toPolyMod, SparsePolyZZ.toPoly]
    rw [hsymmetricMod, htrial, hinitial] at hprimitiveMod'
    exact hprimitiveMod'.symm
  exact ⟨content, equation base₁ hbase₁ hdivides₁,
    equation base₂ hbase₂ hdivides₂⟩

private theorem selectedSourceProduct_ne_zero_of_irreducible
    (base : Nat) [Fact (Nat.Prime base)]
    (activeLifted : Array SparsePolyZZ) (candidate : Array Nat)
    (hbound : ∀ position (hposition : position < candidate.size),
      candidate[position] < activeLifted.size)
    (hirreducible : ∀ index (hindex : index < activeLifted.size),
      Irreducible (Refinement.StrictHensel.toPolyMod base
        activeLifted[index])) :
    ((selectSourceIndices activeLifted.toList candidate.toList).map
      (Refinement.StrictHensel.toPolyMod base)).prod ≠ 0 := by
  apply List.prod_ne_zero
  intro hzero
  rw [List.mem_map] at hzero
  rcases hzero with ⟨source, hsource, hsourceZero⟩
  unfold selectSourceIndices at hsource
  rw [List.mem_map] at hsource
  rcases hsource with ⟨index, hindex, hsourceEq⟩
  rcases List.mem_iff_getElem.mp hindex with ⟨position, hposition, hindexEq⟩
  have hpositionArray : position < candidate.size := by simpa using hposition
  have hcandidate : candidate[position] = index := by
    rw [← Array.getElem_toList hpositionArray]
    exact hindexEq
  have hactive := hbound position hpositionArray
  subst source
  rw [← hcandidate, getElem!_pos activeLifted.toList candidate[position]
    (by simpa using hactive), Array.getElem_toList hactive] at hsourceZero
  exact (hirreducible candidate[position] hactive).ne_zero hsourceZero

/-- Every nonempty legal candidate over the live Hensel array executes the
complete generated attempt without a raw fault.  A non-factor returns
`.rejected`; a factor returns `.extracted`.  Nonzeroness of the physical trial
and primitive divisor is derived modulo the selected prime from the actual
candidate product. -/
theorem zassenhausAttempt_complete
    (fStar : SparsePolyZZ) (activeLifted : Array SparsePolyZZ)
    (modulus base : Nat) [Fact (Nat.Prime base)] (candidate : Array Nat)
    (hmodulus : 0 < modulus) (hbase : 0 < base) (hdivides : base ∣ modulus)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical fStar)
    (hnonempty : 0 < fStar.size)
    (hlegal : LegalCombination activeLifted.size candidate.size candidate)
    (hcandidate : 0 < candidate.size)
    (hactiveFits : activeLifted.size ≤ 2 ^ 31)
    (hleading : (fStar[0].2 : ZMod base) ≠ 0)
    (hirreducible : ∀ index (hindex : index < activeLifted.size),
      Irreducible (Refinement.StrictHensel.toPolyMod base
        activeLifted[index])) :
    ∃ result,
      Generated.StrictRecombine.zassenhausAttempt fStar activeLifted
        (modulus : ZZ) candidate = .ok result := by
  have hbound : ∀ position (hposition : position < candidate.size),
      candidate[position] < activeLifted.size := by
    intro position hposition
    simpa [getElem!_pos candidate position hposition] using
      hlegal.2.2 position hposition
  have hactiveNonempty : ∀ position (hposition : position < candidate.size),
      activeLifted[candidate[position]]!.isEmpty = false := by
    intro position hposition
    have hactive := hbound position hposition
    have hfactorNe := (hirreducible candidate[position] hactive).ne_zero
    by_contra hempty
    have harrayEmpty : activeLifted[candidate[position]]! = #[] := by
      apply Array.isEmpty_iff.mp
      simpa using hempty
    apply hfactorNe
    have hsame : activeLifted[candidate[position]] = #[] := by
      simpa [getElem!_pos activeLifted candidate[position] hactive] using
        harrayEmpty
    rw [hsame]
    simp [Refinement.StrictHensel.toPolyMod, SparsePolyZZ.toPoly]
  have hleadingRun := selectedLeadingProductLoop_succeeds candidate
    activeLifted 0 fStar[0].2 hbound hactiveNonempty
  let leadingProduct := fStar[0].2 *
    (selectedLeadingValues candidate activeLifted 0).prod
  change Generated.StrictRecombine.selectedLeadingProductLoop candidate
      activeLifted 0 fStar[0].2 = .ok leadingProduct at hleadingRun
  by_cases hleadingPrune :
      ZZ.symmetricMod leadingProduct (modulus : ZZ) ≠ 0 ∧
        ZZ.fdiv_r 0 (fStar[0].2 * fStar[0].2)
          (ZZ.symmetricMod leadingProduct (modulus : ZZ)) ≠ 0
  · refine ⟨.rejected, ?_⟩
    unfold Generated.StrictRecombine.zassenhausAttempt
    rw [dif_pos hnonempty]
    simp only [hleadingRun]
    rw [if_pos hleadingPrune]
  · have hconstantRun := selectedConstantProductLoop_succeeds candidate
      activeLifted 0 fStar[0].2 hbound
    let constantProduct := fStar[0].2 *
      (selectedConstantValues candidate activeLifted 0).prod
    change Generated.StrictRecombine.selectedConstantProductLoop candidate
        activeLifted 0 fStar[0].2 = .ok constantProduct at hconstantRun
    by_cases hconstantPrune :
        ZZ.symmetricMod constantProduct (modulus : ZZ) ≠ 0 ∧
          ZZ.fdiv_r 0 (fStar[0].2 *
            Generated.StrictRecombine.constantTerm fStar)
            (ZZ.symmetricMod constantProduct (modulus : ZZ)) ≠ 0
    · refine ⟨.rejected, ?_⟩
      unfold Generated.StrictRecombine.zassenhausAttempt
      rw [dif_pos hnonempty]
      simp only [hleadingRun]
      rw [if_neg hleadingPrune]
      simp only [hconstantRun]
      rw [if_pos hconstantPrune]
    · have hfits : ∀ position (hposition : position < candidate.size),
          candidate[position] < 2 ^ 31 := by
        intro position hposition
        exact lt_of_lt_of_le (hbound position hposition) hactiveFits
      rcases combinationToInt32_toList candidate hfits with
        ⟨candidate32, hconvert, _⟩
      have hvalid := combinationToInt32_candidate_valid candidate
        activeLifted.size candidate32 hbound hactiveFits hconvert
      rcases trialProductLoop_complete ⟨()⟩ candidate32 activeLifted
          (modulus : ZZ) 0 #[(⟨0⟩, fStar[0].2)]
          (by exact_mod_cast hmodulus.ne') hvalid with
        ⟨product, hproduct⟩
      have hinitialCanonical : StrictPolynomialMod.SparsePolyZZCanonical
          #[(⟨0⟩, fStar[0].2)] := by
        constructor
        · simp
        · intro term hterm
          simp at hterm
          subst term
          exact hcanonical.2 fStar[0]
            (Array.getElem_mem_toList hnonempty)
      have hproductCanonical := trialProductLoop_canonical ⟨()⟩ candidate32
        activeLifted (modulus : ZZ) 0 #[(⟨0⟩, fStar[0].2)] product
        hinitialCanonical hproduct
      rcases symmetricModRaw_complete product (modulus : ZZ)
          (by exact_mod_cast hmodulus) with ⟨symmetric, hsymmetric⟩
      have hsymmetricCanonical := symmetricModRaw_canonical product symmetric
        modulus hmodulus hproductCanonical hsymmetric
      rcases primitiveRaw_complete symmetric hsymmetricCanonical with
        ⟨content, factor, hprimitive⟩
      have hselectedNe := selectedSourceProduct_ne_zero_of_irreducible base
        activeLifted candidate hbound hirreducible
      have hinitialMod : Refinement.StrictHensel.toPolyMod base
          #[(⟨0⟩, fStar[0].2)] = Polynomial.C (fStar[0].2 : ZMod base) := by
        simp [Refinement.StrictHensel.toPolyMod, SparsePolyZZ.toPoly]
      have htrialMod := trialProductLoop_source_indices_refines_of_dvd
        modulus base hmodulus hbase hdivides candidate activeLifted candidate32
        #[(⟨0⟩, fStar[0].2)] product hbound hactiveFits hconvert hproduct
      have hproductModNe : Refinement.StrictHensel.toPolyMod base product ≠ 0 := by
        rw [htrialMod, hinitialMod]
        exact mul_ne_zero (Polynomial.C_ne_zero.mpr hleading) hselectedNe
      have hsymmetricMod := symmetricModRaw_toPolyMod_of_dvd product symmetric
        modulus base hmodulus hbase hdivides hsymmetric
      have hsymmetricNe : SparsePolyZZ.toPoly symmetric ≠ 0 := by
        intro hzero
        apply hproductModNe
        rw [← hsymmetricMod]
        simp [Refinement.StrictHensel.toPolyMod, hzero]
      have hfactorNe : SparsePolyZZ.toPoly factor ≠ 0 := by
        intro hzero
        apply hsymmetricNe
        rw [primitiveRaw_toPoly symmetric factor content hprimitive, hzero]
        simp
      have hfactorCanonical := primitiveRaw_canonical symmetric factor content
        hsymmetricCanonical hprimitive
      have hfactorNonempty : 0 < factor.size := by
        by_contra hnot
        have hzero : factor.size = 0 := Nat.eq_zero_of_not_pos hnot
        have hempty : factor = #[] := Array.size_eq_zero_iff.mp hzero
        apply hfactorNe
        simp [hempty, SparsePolyZZ.toPoly]
      rcases exactDivmodRaw_complete fStar factor hcanonical hfactorCanonical
          hfactorNonempty with ⟨quotient, remainder, hdivmod⟩
      by_cases hremainder : remainder.isEmpty = true
      · have hquotientCanonical := exactDivmodRaw_quotient_canonical fStar
          factor quotient remainder hcanonical.2 hdivmod
        rcases primitiveRaw_complete quotient hquotientCanonical with
          ⟨quotientContent, quotientPrimitive, hquotientPrimitive⟩
        refine ⟨.extracted factor quotientPrimitive, ?_⟩
        unfold Generated.StrictRecombine.zassenhausAttempt
        rw [dif_pos hnonempty]
        simp only [hleadingRun]
        rw [if_neg hleadingPrune]
        simp only [hconstantRun]
        rw [if_neg hconstantPrune]
        simp only [hconvert, hproduct, hsymmetric, hprimitive, hdivmod]
        rw [if_pos hremainder]
        simp only [hquotientPrimitive]
      · refine ⟨.rejected, ?_⟩
        unfold Generated.StrictRecombine.zassenhausAttempt
        rw [dif_pos hnonempty]
        simp only [hleadingRun]
        rw [if_neg hleadingPrune]
        simp only [hconstantRun]
        rw [if_neg hconstantPrune]
        simp only [hconvert, hproduct, hsymmetric, hprimitive, hdivmod]
        rw [if_neg hremainder]

/-- The concrete fixed-size scan is total under the live Hensel invariants.
Every visited legal candidate runs the actual attempt above; recursion follows
the generated next-combination rank. -/
theorem scanZassenhausCombinations_complete
    {count : Nat} (fStar : SparsePolyZZ)
    (activeLifted : Array SparsePolyZZ) (modulus base : Nat)
    [Fact (Nat.Prime base)] (start : Array Nat)
    (hmodulus : 0 < modulus) (hbase : 0 < base) (hdivides : base ∣ modulus)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical fStar)
    (hnonempty : 0 < fStar.size) (hcount : 0 < count)
    (hactiveFits : activeLifted.size ≤ 2 ^ 31)
    (hleading : (fStar[0].2 : ZMod base) ≠ 0)
    (hirreducible : ∀ index (hindex : index < activeLifted.size),
      Irreducible (Refinement.StrictHensel.toPolyMod base
        activeLifted[index]))
    (hvalidStart : LegalCombination activeLifted.size count start) :
    ∃ result,
      Generated.StrictRecombine.scanZassenhausCombinations
        (concreteCombinationTermination activeLifted.size count)
        fStar activeLifted (modulus : ZZ) start hvalidStart = .ok result := by
  let termination := concreteCombinationTermination activeLifted.size count
  induction hmeasure : termination.rank start using Nat.strong_induction_on
      generalizing start with
  | h measure ih =>
      rw [Generated.StrictRecombine.scanZassenhausCombinations]
      have hstartPositive : 0 < start.size := by rw [hvalidStart.1]; exact hcount
      have hvalidAttempt : LegalCombination activeLifted.size start.size
          start := by
        simpa [hvalidStart.1] using hvalidStart
      rcases zassenhausAttempt_complete fStar activeLifted modulus base start
          hmodulus hbase hdivides hcanonical hnonempty hvalidAttempt
          hstartPositive hactiveFits hleading hirreducible with
        ⟨attempt, hattempt⟩
      cases attempt with
      | extracted factor quotient =>
          simp only [hattempt]
          exact ⟨.extracted factor quotient start hvalidStart.1, rfl⟩
      | rejected =>
          simp only [hattempt]
          split
          next next hnext => exact ⟨.exhausted, rfl⟩
          next next hnext =>
            have hnextValid := termination.next_valid start next
              hvalidStart hnext
            exact ih (termination.rank next)
              (by
                rw [← hmeasure]
                exact termination.next_decreases start next hvalidStart hnext)
              next hnextValid rfl

/-- If one legal fixed-size candidate is proved to extract, the literal scan
from the generated initial combination returns an actual extraction.  The
exhausted alternative is ruled out by the previously proved no-omission
theorem, not by choosing a scan result. -/
theorem scanZassenhausCombinations_extracts_of_candidate
    {count : Nat} (fStar : SparsePolyZZ)
    (activeLifted : Array SparsePolyZZ) (modulus base : Nat)
    [Fact (Nat.Prime base)] (target : Array Nat)
    (hmodulus : 0 < modulus) (hbase : 0 < base) (hdivides : base ∣ modulus)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical fStar)
    (hnonempty : 0 < fStar.size) (hcount : 0 < count)
    (hfits : count ≤ activeLifted.size)
    (hactiveFits : activeLifted.size ≤ 2 ^ 31)
    (hleading : (fStar[0].2 : ZMod base) ≠ 0)
    (hirreducible : ∀ index (hindex : index < activeLifted.size),
      Irreducible (Refinement.StrictHensel.toPolyMod base
        activeLifted[index]))
    (htarget : LegalCombination activeLifted.size count target)
    (factor quotient : SparsePolyZZ)
    (hattempt : Generated.StrictRecombine.zassenhausAttempt fStar activeLifted
      (modulus : ZZ) target = .ok (.extracted factor quotient)) :
    ∃ extractedFactor extractedQuotient candidate candidateSize,
      Generated.StrictRecombine.scanZassenhausCombinations
        (concreteCombinationTermination activeLifted.size count)
        fStar activeLifted (modulus : ZZ)
        (Generated.StrictRecombine.initialCombination count)
        (initialCombination_legal activeLifted.size count hfits) =
          .ok (.extracted extractedFactor extractedQuotient candidate
            candidateSize) := by
  let initial := Generated.StrictRecombine.initialCombination count
  let hinitial := initialCombination_legal activeLifted.size count hfits
  rcases scanZassenhausCombinations_complete fStar activeLifted modulus base
      initial hmodulus hbase hdivides hcanonical hnonempty hcount hactiveFits
      hleading hirreducible hinitial with ⟨result, hscan⟩
  cases result with
  | extracted extractedFactor extractedQuotient candidate candidateSize =>
      exact ⟨extractedFactor, extractedQuotient, candidate, candidateSize,
        hscan⟩
  | exhausted =>
      have hrejected := scanZassenhausCombinations_exhausted_rejects_all
        fStar activeLifted (modulus : ZZ) hfits hscan target htarget
      simp [hattempt] at hrejected

/-- Public form of the preceding scan theorem through the sole concrete
Zassenhaus termination bundle, without exposing the private rank definition. -/
theorem zassenhausFixedSizeScan_extracts_of_candidate
    {count : Nat} (fStar : SparsePolyZZ)
    (activeLifted : Array SparsePolyZZ) (modulus base : Nat)
    [Fact (Nat.Prime base)] (target : Array Nat)
    (hmodulus : 0 < modulus) (hbase : 0 < base) (hdivides : base ∣ modulus)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical fStar)
    (hnonempty : 0 < fStar.size) (hcount : 0 < count)
    (hfits : count ≤ activeLifted.size)
    (hactiveFits : activeLifted.size ≤ 2 ^ 31)
    (hleading : (fStar[0].2 : ZMod base) ≠ 0)
    (hirreducible : ∀ index (hindex : index < activeLifted.size),
      Irreducible (Refinement.StrictHensel.toPolyMod base
        activeLifted[index]))
    (htarget : LegalCombination activeLifted.size count target)
    (factor quotient : SparsePolyZZ)
    (hattempt : Generated.StrictRecombine.zassenhausAttempt fStar activeLifted
      (modulus : ZZ) target = .ok (.extracted factor quotient)) :
    ∃ extractedFactor extractedQuotient candidate candidateSize,
      Generated.StrictRecombine.scanZassenhausCombinations
        (concreteZassenhausTermination.combinations activeLifted.size count)
        fStar activeLifted (modulus : ZZ)
        (Generated.StrictRecombine.initialCombination count)
        (concreteZassenhausTermination.initial_valid activeLifted.size count
          hfits) = .ok (.extracted extractedFactor extractedQuotient candidate
            candidateSize) := by
  simpa [concreteZassenhausTermination] using
    scanZassenhausCombinations_extracts_of_candidate fStar activeLifted modulus
      base target hmodulus hbase hdivides hcanonical hnonempty hcount hfits
      hactiveFits hleading hirreducible htarget factor quotient hattempt

/-- Under the live selected-prime invariants, the factor returned by the
actual successful attempt is associated modulo `base` to the exact candidate
subproduct.  Both scalar units are derived from the concrete execution
equation and nonzeroness; no unit or association witness is supplied by an
oracle. -/
theorem zassenhausAttempt_extracted_factor_mod_associated_selected
    (fStar factor quotientPrimitive : SparsePolyZZ)
    (activeLifted : Array SparsePolyZZ) (modulus base : Nat)
    [Fact (Nat.Prime base)]
    (candidate : Array Nat)
    (hmodulus : 0 < modulus) (hbase : 0 < base)
    (hdivides : base ∣ modulus)
    (hbound : ∀ position (hposition : position < candidate.size),
      candidate[position] < activeLifted.size)
    (hactiveFits : activeLifted.size ≤ 2 ^ 31)
    (hleading : (fStar[0]!.2 : ZMod base) ≠ 0)
    (hirreducible : ∀ index (hindex : index < activeLifted.size),
      Irreducible (Refinement.StrictHensel.toPolyMod base
        activeLifted[index]))
    (hrun : Generated.StrictRecombine.zassenhausAttempt fStar activeLifted
      (modulus : ZZ) candidate =
        .ok (.extracted factor quotientPrimitive)) :
    Associated (Refinement.StrictHensel.toPolyMod base factor)
      ((selectSourceIndices activeLifted.toList candidate.toList).map
        (Refinement.StrictHensel.toPolyMod base)).prod := by
  rcases zassenhausAttempt_extracted_factor_mod_eq_selected fStar factor
      quotientPrimitive activeLifted modulus base candidate hmodulus hbase
      hdivides hbound hactiveFits hrun with ⟨content, heq⟩
  let selected :=
    ((selectSourceIndices activeLifted.toList candidate.toList).map
      (Refinement.StrictHensel.toPolyMod base)).prod
  have hselected : selected ≠ 0 :=
    selectedSourceProduct_ne_zero_of_irreducible base activeLifted candidate
      hbound hirreducible
  have hleadingC : Polynomial.C (fStar[0]!.2 : ZMod base) ≠ 0 := by
    exact Polynomial.C_ne_zero.mpr hleading
  have hrhs : Polynomial.C (fStar[0]!.2 : ZMod base) * selected ≠ 0 :=
    mul_ne_zero hleadingC hselected
  have hlhs : Polynomial.C (content : ZMod base) *
      Refinement.StrictHensel.toPolyMod base factor ≠ 0 := by
    rw [heq]
    exact hrhs
  have hcontent : (content : ZMod base) ≠ 0 := by
    intro hzero
    apply hlhs
    simp [hzero]
  have hcontentUnit : IsUnit
      (Polynomial.C (content : ZMod base)) :=
    Polynomial.isUnit_C.mpr (isUnit_iff_ne_zero.mpr hcontent)
  have hleadingUnit : IsUnit
      (Polynomial.C (fStar[0]!.2 : ZMod base)) :=
    Polynomial.isUnit_C.mpr (isUnit_iff_ne_zero.mpr hleading)
  have hleft : Associated
      (Polynomial.C (content : ZMod base) *
        Refinement.StrictHensel.toPolyMod base factor)
      (Refinement.StrictHensel.toPolyMod base factor) :=
    (associated_isUnit_mul_left_iff hcontentUnit).mpr (Associated.refl _)
  have hright : Associated
      (Polynomial.C (fStar[0]!.2 : ZMod base) * selected) selected :=
    (associated_isUnit_mul_left_iff hleadingUnit).mpr (Associated.refl _)
  exact hleft.symm.trans ((Associated.of_eq heq).trans hright)

/-- The factor returned by a successful generated Zassenhaus attempt is the
actual primitive result built from the canonical trial-product and
symmetric-mod traces. -/
theorem zassenhausAttempt_extracted_factor_canonical_primitive
    (fStar factor quotientPrimitive : SparsePolyZZ)
    (activeLifted : Array SparsePolyZZ) (modulus : ZZ)
    (candidate : Array Nat)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical fStar)
    (hrun : Generated.StrictRecombine.zassenhausAttempt fStar activeLifted
      modulus candidate = .ok (.extracted factor quotientPrimitive)) :
    StrictPolynomialMod.SparsePolyZZCanonical factor ∧
      (SparsePolyZZ.toPoly factor).IsPrimitive := by
  unfold Generated.StrictRecombine.zassenhausAttempt at hrun
  split at hrun
  next hfstar =>
    dsimp at hrun
    cases hleading : Generated.StrictRecombine.selectedLeadingProductLoop
        candidate activeLifted 0 fStar[0].2 with
    | error fault => simp [hleading] at hrun
    | ok leadingProduct =>
      simp only [hleading] at hrun
      split at hrun
      next hpruned => simp at hrun
      next hleadingAccepted =>
        cases hconstant : Generated.StrictRecombine.selectedConstantProductLoop
            candidate activeLifted 0 fStar[0].2 with
        | error fault => simp [hconstant] at hrun
        | ok constantProduct =>
          simp only [hconstant] at hrun
          split at hrun
          next hpruned => simp at hrun
          next hconstantAccepted =>
            cases hconvert : Generated.StrictRecombine.combinationToInt32
                candidate with
            | error fault => simp [hconvert] at hrun
            | ok candidate32 =>
              simp only [hconvert] at hrun
              cases hproduct : Generated.StrictRecombine.trialProductLoop
                  ⟨()⟩ candidate32 activeLifted modulus 0
                  #[(⟨0⟩, fStar[0].2)] with
              | error fault => simp [hproduct] at hrun
              | ok product =>
                simp only [hproduct] at hrun
                cases hsymmetric : Generated.StrictRecombine.symmetricModRaw
                    product modulus with
                | error fault => simp [hsymmetric] at hrun
                | ok symmetric =>
                  simp only [hsymmetric] at hrun
                  cases hprimitive : Generated.StrictRecombine.primitiveRaw
                      symmetric with
                  | error fault => simp [hprimitive] at hrun
                  | ok primitiveResult =>
                    rcases primitiveResult with
                      ⟨symmetricContent, recoveredFactor⟩
                    simp only [hprimitive] at hrun
                    cases hdivmod : Generated.StrictRecombine.exactDivmodRaw
                        fStar recoveredFactor with
                    | error fault => simp [hdivmod] at hrun
                    | ok divResult =>
                      rcases divResult with ⟨quotient, remainder⟩
                      simp only [hdivmod] at hrun
                      by_cases hremainder : remainder.isEmpty = true
                      · simp only [hremainder, if_true] at hrun
                        cases hquotientPrimitive :
                            Generated.StrictRecombine.primitiveRaw quotient with
                        | error fault => simp [hquotientPrimitive] at hrun
                        | ok quotientResult =>
                          rcases quotientResult with
                            ⟨quotientContent, recoveredQuotient⟩
                          simp only [hquotientPrimitive] at hrun
                          have hinitial :
                              StrictPolynomialMod.SparsePolyZZCanonical
                                #[(⟨0⟩, fStar[0].2)] := by
                            unfold StrictPolynomialMod.SparsePolyZZCanonical
                            constructor
                            · simp
                            · intro term hterm
                              simp at hterm
                              subst term
                              exact hcanonical.2 fStar[0]
                                (Array.getElem_mem_toList hfstar)
                          have hproductCanonical :=
                            trialProductLoop_canonical ⟨()⟩ candidate32
                              activeLifted modulus 0
                              #[(⟨0⟩, fStar[0].2)] product hinitial hproduct
                          have hmodulus : 0 < modulus := by
                            unfold Generated.StrictRecombine.symmetricModRaw at hsymmetric
                            split at hsymmetric
                            next hpositive => exact hpositive
                            next hpositive => contradiction
                          have hsymmetricCanonical :=
                            symmetricModRaw_canonical product symmetric
                              modulus.toNat (Int.pos_iff_toNat_pos.mp hmodulus)
                              hproductCanonical (by
                                simpa [Int.toNat_of_nonneg hmodulus.le] using
                                  hsymmetric)
                          have hfactorCanonical :=
                            primitiveRaw_canonical symmetric recoveredFactor
                              symmetricContent hsymmetricCanonical hprimitive
                          have hfactorNonempty : 0 < recoveredFactor.size := by
                            unfold Generated.StrictRecombine.exactDivmodRaw at hdivmod
                            rw [Generated.StrictRecombine.exactDivmodLoop] at hdivmod
                            rw [dif_pos hfstar] at hdivmod
                            split at hdivmod
                            next hpositive => exact hpositive
                            next hpositive => contradiction
                          have hsymmetricNonempty : 0 < symmetric.size := by
                            by_contra hnot
                            have hzero : symmetric.size = 0 :=
                              Nat.eq_zero_of_not_pos hnot
                            have hempty : symmetric.isEmpty = true := by
                              simpa [Array.isEmpty, hzero]
                            rw [Generated.StrictRecombine.primitiveRaw,
                              dif_pos hempty] at hprimitive
                            have hout := Except.ok.inj hprimitive
                            have hsame := congrArg Prod.snd hout
                            simp only [Prod.snd] at hsame
                            have : recoveredFactor.size = 0 := by
                              rw [← hsame]
                              exact hzero
                            omega
                          have hfactorPrimitive :=
                            primitiveRaw_isPrimitive symmetric recoveredFactor
                              symmetricContent hsymmetricNonempty
                              hsymmetricCanonical hprimitive
                          have hout := Except.ok.inj hrun
                          injection hout with hfactor hquotient
                          subst factor
                          exact ⟨hfactorCanonical, hfactorPrimitive⟩
                      · simp only [hremainder, if_false] at hrun
                        simp at hrun
  next hfstar => contradiction

theorem zassenhausAttempt_extracted_quotient_canonical
    (fStar factor quotientPrimitive : SparsePolyZZ)
    (activeLifted : Array SparsePolyZZ) (modulus : ZZ)
    (candidate : Array Nat)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical fStar)
    (hrun : Generated.StrictRecombine.zassenhausAttempt fStar activeLifted
      modulus candidate = .ok (.extracted factor quotientPrimitive)) :
    StrictPolynomialMod.SparsePolyZZCanonical quotientPrimitive := by
  rcases zassenhausAttempt_extracted_quotient_trace fStar factor
      quotientPrimitive activeLifted modulus candidate hrun with
    ⟨quotient, quotientContent, hdivide, hprimitive⟩
  have hquotient := exactDivmodRaw_quotient_canonical fStar factor quotient #[]
    hcanonical.2 hdivide
  exact primitiveRaw_canonical quotient quotientPrimitive quotientContent
    hquotient hprimitive

theorem zassenhausAttempt_extracted_quotient_canonical_primitive
    (fStar factor quotientPrimitive : SparsePolyZZ)
    (activeLifted : Array SparsePolyZZ) (modulus : ZZ)
    (candidate : Array Nat)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical fStar)
    (hnonempty : 0 < fStar.size)
    (hrun : Generated.StrictRecombine.zassenhausAttempt fStar activeLifted
      modulus candidate = .ok (.extracted factor quotientPrimitive)) :
    StrictPolynomialMod.SparsePolyZZCanonical quotientPrimitive ∧
      (SparsePolyZZ.toPoly quotientPrimitive).IsPrimitive := by
  rcases zassenhausAttempt_extracted_quotient_trace fStar factor
      quotientPrimitive activeLifted modulus candidate hrun with
    ⟨quotient, quotientContent, hdivide, hprimitive⟩
  have hquotientCanonical := exactDivmodRaw_quotient_canonical fStar factor
    quotient #[] hcanonical.2 hdivide
  have hfStarNe : SparsePolyZZ.toPoly fStar ≠ 0 := by
    intro hzeroPoly
    have hleading := sparsePolyZZ_leadingCoeff_eq_head fStar hcanonical
      hnonempty
    rw [hzeroPoly] at hleading
    have hheadNonzero := hcanonical.2 fStar[0]
      (Array.getElem_mem_toList hnonempty)
    exact hheadNonzero (by simpa using hleading.symm)
  have hquotientNonempty : 0 < quotient.size := by
    by_contra hnot
    have hzero : quotient.size = 0 := Nat.eq_zero_of_not_pos hnot
    have hempty : quotient = #[] := Array.size_eq_zero_iff.mp hzero
    have hsemantic := exactDivmodRaw_toPoly fStar factor quotient #[] hdivide
    rw [hempty] at hsemantic
    simp [SparsePolyZZ.toPoly] at hsemantic
    exact hfStarNe hsemantic
  exact ⟨primitiveRaw_canonical quotient quotientPrimitive quotientContent
      hquotientCanonical hprimitive,
    primitiveRaw_isPrimitive quotient quotientPrimitive quotientContent
      hquotientNonempty hquotientCanonical hprimitive⟩

theorem zassenhausAttempt_extracted_toPoly
    (fStar factor quotientPrimitive : SparsePolyZZ)
    (activeLifted : Array SparsePolyZZ) (modulus : ZZ)
    (candidate : Array Nat)
    (hrun : Generated.StrictRecombine.zassenhausAttempt fStar activeLifted
      modulus candidate = .ok
        (.extracted factor quotientPrimitive)) :
    ∃ scalar : Int,
      SparsePolyZZ.toPoly fStar = Polynomial.C scalar *
        (SparsePolyZZ.toPoly factor *
          SparsePolyZZ.toPoly quotientPrimitive) := by
  unfold Generated.StrictRecombine.zassenhausAttempt at hrun
  split at hrun
  next hfstar =>
    dsimp at hrun
    cases hleading : Generated.StrictRecombine.selectedLeadingProductLoop
        candidate activeLifted 0 fStar[0].2 with
    | error fault => simp [hleading] at hrun
    | ok leadingProduct =>
      simp only [hleading] at hrun
      split at hrun
      next hpruned => simp at hrun
      next hleadingAccepted =>
        cases hconstant :
            Generated.StrictRecombine.selectedConstantProductLoop candidate
              activeLifted 0 fStar[0].2 with
        | error fault => simp [hconstant] at hrun
        | ok constantProduct =>
          simp only [hconstant] at hrun
          split at hrun
          next hpruned => simp at hrun
          next hconstantAccepted =>
            cases hconvert : Generated.StrictRecombine.combinationToInt32
                candidate with
            | error fault => simp [hconvert] at hrun
            | ok candidate32 =>
              simp only [hconvert] at hrun
              cases hproduct : Generated.StrictRecombine.trialProductLoop
                  ⟨()⟩ candidate32 activeLifted modulus 0
                  #[(⟨0⟩, fStar[0].2)] with
              | error fault => simp [hproduct] at hrun
              | ok product =>
                simp only [hproduct] at hrun
                cases hsymmetric : Generated.StrictRecombine.symmetricModRaw
                    product modulus with
                | error fault => simp [hsymmetric] at hrun
                | ok symmetric =>
                  simp only [hsymmetric] at hrun
                  cases hprimitive : Generated.StrictRecombine.primitiveRaw
                      symmetric with
                  | error fault => simp [hprimitive] at hrun
                  | ok primitiveResult =>
                    rcases primitiveResult with ⟨symmetricContent,
                      recoveredFactor⟩
                    simp only [hprimitive] at hrun
                    cases hdivmod : Generated.StrictRecombine.exactDivmodRaw
                        fStar recoveredFactor with
                    | error fault => simp [hdivmod] at hrun
                    | ok divResult =>
                      rcases divResult with ⟨quotient, remainder⟩
                      simp only [hdivmod] at hrun
                      by_cases hremainder : remainder.isEmpty = true
                      · simp only [hremainder, if_true] at hrun
                        have hremainderEmpty : remainder = #[] :=
                          Array.isEmpty_iff.mp hremainder
                        subst remainder
                        cases hquotientPrimitive :
                            Generated.StrictRecombine.primitiveRaw quotient with
                        | error fault => simp [hquotientPrimitive] at hrun
                        | ok quotientResult =>
                          rcases quotientResult with ⟨quotientContent,
                            recoveredQuotient⟩
                          simp only [hquotientPrimitive] at hrun
                          have hout := Except.ok.inj hrun
                          injection hout with hfactor hquotient
                          subst factor
                          subst quotientPrimitive
                          exact ⟨quotientContent,
                            successfulTrialExtraction_toPoly fStar
                              recoveredFactor quotient recoveredQuotient
                              quotientContent hdivmod hquotientPrimitive⟩
                      · simp only [hremainder, if_false] at hrun
                        simp at hrun
  next hfstar => contradiction

theorem zassenhausAttempt_extracted_unit_scalar
    (fStar factor quotientPrimitive : SparsePolyZZ)
    (activeLifted : Array SparsePolyZZ) (modulus : ZZ)
    (candidate : Array Nat)
    (hprimitive : (SparsePolyZZ.toPoly fStar).IsPrimitive)
    (hrun : Generated.StrictRecombine.zassenhausAttempt fStar activeLifted
      modulus candidate = .ok
        (.extracted factor quotientPrimitive)) :
    ∃ scalar : Int, IsUnit scalar ∧
      SparsePolyZZ.toPoly fStar = Polynomial.C scalar *
        (SparsePolyZZ.toPoly factor *
          SparsePolyZZ.toPoly quotientPrimitive) := by
  rcases zassenhausAttempt_extracted_toPoly fStar factor quotientPrimitive
      activeLifted modulus candidate hrun with ⟨scalar, hproduct⟩
  have hcontent := congrArg Polynomial.content hproduct
  rw [hprimitive.content_eq_one, Polynomial.content_C_mul] at hcontent
  have hnormalizeUnit : IsUnit (normalize scalar) :=
    IsUnit.of_mul_eq_one _ hcontent.symm
  have hscalarUnit : IsUnit scalar :=
    (associated_normalize scalar).isUnit_iff.mpr hnormalizeUnit
  exact ⟨scalar, hscalarUnit, hproduct⟩

/-- Successful concrete extraction preserves the selected-prime
leading-coefficient invariant in the primitive quotient installed by the
outer loop.  The proof uses the actual exact-division/primitive trace and the
unit scalar forced by primitivity. -/
theorem zassenhausAttempt_extracted_quotient_leading_mod_ne_zero
    (fStar factor quotientPrimitive : SparsePolyZZ)
    (activeLifted : Array SparsePolyZZ) (modulus : ZZ)
    (candidate : Array Nat) (base : Nat) [Fact (Nat.Prime base)]
    (hprimitive : (SparsePolyZZ.toPoly fStar).IsPrimitive)
    (hleading : ((SparsePolyZZ.toPoly fStar).leadingCoeff : ZMod base) ≠ 0)
    (hrun : Generated.StrictRecombine.zassenhausAttempt fStar activeLifted
      modulus candidate = .ok (.extracted factor quotientPrimitive)) :
    ((SparsePolyZZ.toPoly quotientPrimitive).leadingCoeff : ZMod base) ≠ 0 := by
  rcases zassenhausAttempt_extracted_unit_scalar fStar factor
      quotientPrimitive activeLifted modulus candidate hprimitive hrun with
    ⟨scalar, hscalar, heq⟩
  have hleadingEq := congrArg Polynomial.leadingCoeff heq
  rw [Polynomial.leadingCoeff_mul, Polynomial.leadingCoeff_mul] at hleadingEq
  have hconstantLeading : (Polynomial.C scalar).leadingCoeff = scalar := by
    exact Polynomial.leadingCoeff_C scalar
  rw [hconstantLeading] at hleadingEq
  have hleadingCast := congrArg (fun coefficient : Int =>
    (coefficient : ZMod base)) hleadingEq
  intro hquotientZero
  apply hleading
  change ((SparsePolyZZ.toPoly fStar).leadingCoeff : ZMod base) =
    ((scalar * ((SparsePolyZZ.toPoly factor).leadingCoeff *
      (SparsePolyZZ.toPoly quotientPrimitive).leadingCoeff) : Int) :
        ZMod base) at hleadingCast
  rw [Int.cast_mul, Int.cast_mul, hquotientZero, mul_zero, mul_zero] at hleadingCast
  exact hleadingCast

theorem scanZassenhausCombinations_extracted_unit_scalar
    {upper count : Nat}
    (termination : Generated.StrictRecombine.CombinationTermination upper count)
    (fStar factor quotientPrimitive : SparsePolyZZ)
    (activeLifted : Array SparsePolyZZ) (modulus : ZZ)
    (start candidate : Array Nat)
    (hcandidateSize : candidate.size = count)
    (hprimitive : (SparsePolyZZ.toPoly fStar).IsPrimitive)
    (hvalidStart : termination.valid start)
    (hrun : Generated.StrictRecombine.scanZassenhausCombinations termination
      fStar activeLifted modulus start hvalidStart = .ok
        (.extracted factor quotientPrimitive candidate hcandidateSize)) :
    ∃ scalar : Int, IsUnit scalar ∧
      SparsePolyZZ.toPoly fStar = Polynomial.C scalar *
        (SparsePolyZZ.toPoly factor *
          SparsePolyZZ.toPoly quotientPrimitive) := by
  induction hmeasure : termination.rank start using Nat.strong_induction_on
      generalizing start with
  | h measure ih =>
      rw [Generated.StrictRecombine.scanZassenhausCombinations] at hrun
      cases hattempt : Generated.StrictRecombine.zassenhausAttempt fStar
          activeLifted modulus start with
      | error fault => simp [hattempt] at hrun
      | ok attempt =>
        cases attempt with
        | extracted extractedFactor extractedQuotient =>
          simp only [hattempt] at hrun
          have hout := Except.ok.inj hrun
          injection hout with hfactor hquotient hcandidate
          subst factor
          subst quotientPrimitive
          exact zassenhausAttempt_extracted_unit_scalar fStar extractedFactor
            extractedQuotient activeLifted modulus start hprimitive hattempt
        | rejected =>
          simp only [hattempt] at hrun
          split at hrun
          next next hnext => simp at hrun
          next next hnext =>
            have hdecrease := termination.next_decreases start next hvalidStart
              hnext
            have hvalidNext := termination.next_valid start next hvalidStart hnext
            rw [hmeasure] at hdecrease
            exact ih (termination.rank next)
              hdecrease next hvalidNext hrun rfl

theorem scanZassenhausCombinations_extracted_quotient_canonical
    {upper count : Nat}
    (termination : Generated.StrictRecombine.CombinationTermination upper count)
    (fStar factor quotientPrimitive : SparsePolyZZ)
    (activeLifted : Array SparsePolyZZ) (modulus : ZZ)
    (start candidate : Array Nat)
    (hcandidateSize : candidate.size = count)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical fStar)
    (hvalidStart : termination.valid start)
    (hrun : Generated.StrictRecombine.scanZassenhausCombinations termination
      fStar activeLifted modulus start hvalidStart = .ok
        (.extracted factor quotientPrimitive candidate hcandidateSize)) :
    StrictPolynomialMod.SparsePolyZZCanonical quotientPrimitive := by
  induction hmeasure : termination.rank start using Nat.strong_induction_on
      generalizing start with
  | h measure ih =>
      rw [Generated.StrictRecombine.scanZassenhausCombinations] at hrun
      cases hattempt : Generated.StrictRecombine.zassenhausAttempt fStar
          activeLifted modulus start with
      | error fault => simp [hattempt] at hrun
      | ok attempt =>
        cases attempt with
        | extracted extractedFactor extractedQuotient =>
          simp only [hattempt] at hrun
          have hout := Except.ok.inj hrun
          injection hout with hfactor hquotient hcandidate
          subst quotientPrimitive
          exact zassenhausAttempt_extracted_quotient_canonical fStar
            extractedFactor extractedQuotient activeLifted modulus start
            hcanonical hattempt
        | rejected =>
          simp only [hattempt] at hrun
          split at hrun
          next next hnext => simp at hrun
          next next hnext =>
            have hdecrease := termination.next_decreases start next hvalidStart
              hnext
            have hvalidNext := termination.next_valid start next hvalidStart hnext
            rw [hmeasure] at hdecrease
            exact ih (termination.rank next) hdecrease next hvalidNext hrun rfl

theorem scanZassenhausCombinations_extracted_canonical_primitive
    {upper count : Nat}
    (termination : Generated.StrictRecombine.CombinationTermination upper count)
    (fStar factor quotientPrimitive : SparsePolyZZ)
    (activeLifted : Array SparsePolyZZ) (modulus : ZZ)
    (start candidate : Array Nat)
    (hcandidateSize : candidate.size = count)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical fStar)
    (hnonempty : 0 < fStar.size)
    (hvalidStart : termination.valid start)
    (hrun : Generated.StrictRecombine.scanZassenhausCombinations termination
      fStar activeLifted modulus start hvalidStart = .ok
        (.extracted factor quotientPrimitive candidate hcandidateSize)) :
    (StrictPolynomialMod.SparsePolyZZCanonical factor ∧
        (SparsePolyZZ.toPoly factor).IsPrimitive) ∧
      (StrictPolynomialMod.SparsePolyZZCanonical quotientPrimitive ∧
        (SparsePolyZZ.toPoly quotientPrimitive).IsPrimitive) := by
  induction hmeasure : termination.rank start using Nat.strong_induction_on
      generalizing start with
  | h measure ih =>
      rw [Generated.StrictRecombine.scanZassenhausCombinations] at hrun
      cases hattempt : Generated.StrictRecombine.zassenhausAttempt fStar
          activeLifted modulus start with
      | error fault => simp [hattempt] at hrun
      | ok attempt =>
        cases attempt with
        | extracted extractedFactor extractedQuotient =>
          simp only [hattempt] at hrun
          have hout := Except.ok.inj hrun
          injection hout with hfactor hquotient hcandidate
          subst factor
          subst quotientPrimitive
          exact ⟨zassenhausAttempt_extracted_factor_canonical_primitive fStar
              extractedFactor extractedQuotient activeLifted modulus start
              hcanonical hattempt,
            zassenhausAttempt_extracted_quotient_canonical_primitive fStar
              extractedFactor extractedQuotient activeLifted modulus start
              hcanonical hnonempty hattempt⟩
        | rejected =>
          simp only [hattempt] at hrun
          split at hrun
          next next hnext => simp at hrun
          next next hnext =>
            have hdecrease := termination.next_decreases start next hvalidStart
              hnext
            have hvalidNext := termination.next_valid start next hvalidStart hnext
            rw [hmeasure] at hdecrease
            exact ih (termination.rank next) hdecrease next hvalidNext hrun rfl

/-- A candidate returned by the generated fixed-size scan is one of the
concrete legal combinations reached by `nextCombination`.  This is execution
provenance, not a postulated property of a successful result. -/
theorem scanZassenhausCombinations_extracted_legal
    {upper count : Nat}
    (termination : Generated.StrictRecombine.CombinationTermination upper count)
    (fStar factor quotientPrimitive : SparsePolyZZ)
    (activeLifted : Array SparsePolyZZ) (modulus : ZZ)
    (start candidate : Array Nat) (candidateSize : candidate.size = count)
    (hvalidStart : termination.valid start)
    (hrun : Generated.StrictRecombine.scanZassenhausCombinations termination
      fStar activeLifted modulus start hvalidStart = .ok
        (.extracted factor quotientPrimitive candidate candidateSize)) :
    termination.valid candidate := by
  induction hmeasure : termination.rank start using Nat.strong_induction_on
      generalizing start with
  | h measure ih =>
      rw [Generated.StrictRecombine.scanZassenhausCombinations] at hrun
      cases hattempt : Generated.StrictRecombine.zassenhausAttempt fStar
          activeLifted modulus start with
      | error fault => simp [hattempt] at hrun
      | ok attempt =>
        cases attempt with
        | extracted extractedFactor extractedQuotient =>
          simp only [hattempt] at hrun
          have hout := Except.ok.inj hrun
          injection hout with _ _ hcandidate
          subst candidate
          exact hvalidStart
        | rejected =>
          simp only [hattempt] at hrun
          split at hrun
          next next hnext => simp at hrun
          next next hnext =>
            have hdecrease := termination.next_decreases start next hvalidStart
              hnext
            have hvalidNext := termination.next_valid start next hvalidStart hnext
            rw [hmeasure] at hdecrease
            exact ih (termination.rank next) hdecrease next hvalidNext hrun rfl

/-- Concrete specialization: every candidate physically returned by the
generated scan is sorted, duplicate-free, in bounds, and has the requested
size. -/
theorem concreteScan_extracted_legal
    {count : Nat} (fStar factor quotientPrimitive : SparsePolyZZ)
    (activeLifted : Array SparsePolyZZ) (modulus : ZZ)
    (candidate : Array Nat) (candidateSize : candidate.size = count)
    (hfits : count ≤ activeLifted.size)
    (hrun : Generated.StrictRecombine.scanZassenhausCombinations
      (concreteCombinationTermination activeLifted.size count)
      fStar activeLifted modulus
      (Generated.StrictRecombine.initialCombination count)
      (initialCombination_legal activeLifted.size count hfits) = .ok
        (.extracted factor quotientPrimitive candidate candidateSize)) :
    LegalCombination activeLifted.size count candidate := by
  exact scanZassenhausCombinations_extracted_legal
    (concreteCombinationTermination activeLifted.size count) fStar factor
    quotientPrimitive activeLifted modulus
    (Generated.StrictRecombine.initialCombination count) candidate
    candidateSize (initialCombination_legal activeLifted.size count hfits) hrun

/-- A successful fixed-size scan returns the exact candidate whose literal
`zassenhausAttempt` produced the reported factor and quotient.  This follows
the generated scan recursion; the returned candidate is not reconstructed
from its semantic result. -/
theorem scanZassenhausCombinations_extracted_attempt
    {upper count : Nat}
    (termination : Generated.StrictRecombine.CombinationTermination upper count)
    (fStar factor quotientPrimitive : SparsePolyZZ)
    (activeLifted : Array SparsePolyZZ) (modulus : ZZ)
    (start candidate : Array Nat)
    (hcandidateSize : candidate.size = count)
    (hvalidStart : termination.valid start)
    (hrun : Generated.StrictRecombine.scanZassenhausCombinations termination
      fStar activeLifted modulus start hvalidStart = .ok
        (.extracted factor quotientPrimitive candidate hcandidateSize)) :
    Generated.StrictRecombine.zassenhausAttempt fStar activeLifted modulus
      candidate = .ok (.extracted factor quotientPrimitive) := by
  induction hmeasure : termination.rank start using Nat.strong_induction_on
      generalizing start with
  | h measure ih =>
      rw [Generated.StrictRecombine.scanZassenhausCombinations] at hrun
      cases hattempt : Generated.StrictRecombine.zassenhausAttempt fStar
          activeLifted modulus start with
      | error fault => simp [hattempt] at hrun
      | ok attempt =>
        cases attempt with
        | extracted extractedFactor extractedQuotient =>
          simp only [hattempt] at hrun
          have hout := Except.ok.inj hrun
          injection hout with hfactor hquotient hcandidate
          subst factor
          subst quotientPrimitive
          subst candidate
          exact hattempt
        | rejected =>
          simp only [hattempt] at hrun
          split at hrun
          next next hnext => simp at hrun
          next next hnext =>
            have hdecrease := termination.next_decreases start next
              hvalidStart hnext
            have hvalidNext := termination.next_valid start next hvalidStart
              hnext
            rw [hmeasure] at hdecrease
            exact ih (termination.rank next) hdecrease next hvalidNext hrun rfl

/-- The concrete fixed-size scan preserves the successful candidate's exact
modular association certificate and the leading-coefficient invariant of the
primitive quotient, regardless of how many preceding candidates execute and
reject. -/
theorem scanZassenhausCombinations_extracted_mod_certificate
    {count : Nat}
    (fStar factor quotientPrimitive : SparsePolyZZ)
    (activeLifted : Array SparsePolyZZ) (modulus base : Nat)
    [Fact (Nat.Prime base)]
    (start candidate : Array Nat)
    (hcandidateSize : candidate.size = count)
    (hmodulus : 0 < modulus) (hbase : 0 < base)
    (hdivides : base ∣ modulus)
    (hactiveFits : activeLifted.size ≤ 2 ^ 31)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical fStar)
    (hnonempty : 0 < fStar.size)
    (hprimitive : (SparsePolyZZ.toPoly fStar).IsPrimitive)
    (hleading : ((SparsePolyZZ.toPoly fStar).leadingCoeff : ZMod base) ≠ 0)
    (hirreducible : ∀ index (hindex : index < activeLifted.size),
      Irreducible (Refinement.StrictHensel.toPolyMod base
        activeLifted[index]))
    (hvalidStart : LegalCombination activeLifted.size count start)
    (hrun : Generated.StrictRecombine.scanZassenhausCombinations
      (concreteCombinationTermination activeLifted.size count)
      fStar activeLifted (modulus : ZZ) start hvalidStart = .ok
        (.extracted factor quotientPrimitive candidate hcandidateSize)) :
    Associated (Refinement.StrictHensel.toPolyMod base factor)
        ((selectSourceIndices activeLifted.toList candidate.toList).map
          (Refinement.StrictHensel.toPolyMod base)).prod ∧
      ((SparsePolyZZ.toPoly quotientPrimitive).leadingCoeff : ZMod base) ≠ 0 := by
  let termination := concreteCombinationTermination activeLifted.size count
  induction hmeasure : termination.rank start using Nat.strong_induction_on
      generalizing start with
  | h measure ih =>
      rw [Generated.StrictRecombine.scanZassenhausCombinations] at hrun
      cases hattempt : Generated.StrictRecombine.zassenhausAttempt fStar
          activeLifted (modulus : ZZ) start with
      | error fault => simp [hattempt] at hrun
      | ok attempt =>
        cases attempt with
        | extracted extractedFactor extractedQuotient =>
          simp only [hattempt] at hrun
          have hout := Except.ok.inj hrun
          injection hout with hfactor hquotient hcandidate
          subst factor
          subst quotientPrimitive
          subst candidate
          have hbound : ∀ position (hposition : position < start.size),
              start[position] < activeLifted.size := by
            intro position hposition
            simpa [getElem!_pos start position hposition] using
              hvalidStart.2.2 position hposition
          have hfront : (fStar[0]!.2 : ZMod base) ≠ 0 := by
            intro hzero
            apply hleading
            rw [sparsePolyZZ_leadingCoeff_eq_head fStar hcanonical hnonempty]
            simpa [getElem!_pos fStar 0 hnonempty] using hzero
          exact ⟨
            zassenhausAttempt_extracted_factor_mod_associated_selected fStar
              extractedFactor extractedQuotient activeLifted modulus base
              start hmodulus hbase hdivides hbound hactiveFits hfront
              hirreducible hattempt,
            zassenhausAttempt_extracted_quotient_leading_mod_ne_zero fStar
              extractedFactor extractedQuotient activeLifted (modulus : ZZ)
              start base hprimitive hleading hattempt⟩
        | rejected =>
          simp only [hattempt] at hrun
          split at hrun
          next next hnext => simp at hrun
          next next hnext =>
            have hdecrease := termination.next_decreases start next hvalidStart
              hnext
            have hvalidNext := termination.next_valid start next hvalidStart hnext
            rw [hmeasure] at hdecrease
            exact ih (termination.rank next) hdecrease next hvalidNext hrun rfl

/-- The modular product invariant is preserved by the exact successful scan
and removal executions used by the generated outer loop.  The proof cancels
the physically returned nonzero factor from the mapped integer extraction
equation and the occurrence-sensitive selected/complement partition. -/
theorem scanExtraction_removeCombination_preserves_mod_product
    {count : Nat}
    (fStar factor quotientPrimitive : SparsePolyZZ)
    (activeLifted active' : Array SparsePolyZZ) (modulus base : Nat)
    [Fact (Nat.Prime base)] (candidate : Array Nat)
    (candidateSize : candidate.size = count)
    (hfits : count ≤ activeLifted.size)
    (hmodulus : 0 < modulus) (hbase : 0 < base)
    (hdivides : base ∣ modulus)
    (hactiveFits : activeLifted.size ≤ 2 ^ 31)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical fStar)
    (hnonempty : 0 < fStar.size)
    (hprimitive : (SparsePolyZZ.toPoly fStar).IsPrimitive)
    (hleading : ((SparsePolyZZ.toPoly fStar).leadingCoeff : ZMod base) ≠ 0)
    (hirreducible : ∀ index (hindex : index < activeLifted.size),
      Irreducible (Refinement.StrictHensel.toPolyMod base
        activeLifted[index]))
    (hlive : Associated
      (Polynomial.map (Int.castRingHom (ZMod base))
        (SparsePolyZZ.toPoly fStar))
      ((activeLifted.toList.map
        (Refinement.StrictHensel.toPolyMod base)).prod))
    (hscan : Generated.StrictRecombine.scanZassenhausCombinations
      (concreteCombinationTermination activeLifted.size count)
      fStar activeLifted (modulus : ZZ)
      (Generated.StrictRecombine.initialCombination count)
      (initialCombination_legal activeLifted.size count hfits) = .ok
        (.extracted factor quotientPrimitive candidate candidateSize))
    (hremove : Generated.StrictRecombine.removeCombination candidate
      activeLifted = .ok active') :
    Associated
      (Polynomial.map (Int.castRingHom (ZMod base))
        (SparsePolyZZ.toPoly quotientPrimitive))
      ((active'.toList.map
        (Refinement.StrictHensel.toPolyMod base)).prod) := by
  let initial := Generated.StrictRecombine.initialCombination count
  let hinitial := initialCombination_legal activeLifted.size count hfits
  have hlegal : LegalCombination activeLifted.size count candidate :=
    concreteScan_extracted_legal fStar factor quotientPrimitive activeLifted
      (modulus : ZZ) candidate candidateSize hfits hscan
  have hlegalSize : LegalCombination activeLifted.size candidate.size
      candidate := by simpa [candidateSize] using hlegal
  have hfactorCertificate :=
    scanZassenhausCombinations_extracted_mod_certificate fStar factor
      quotientPrimitive activeLifted modulus base initial candidate
      candidateSize hmodulus hbase hdivides hactiveFits hcanonical hnonempty
      hprimitive hleading hirreducible hinitial hscan
  rcases scanZassenhausCombinations_extracted_unit_scalar
      (concreteCombinationTermination activeLifted.size count) fStar factor
      quotientPrimitive activeLifted (modulus : ZZ) initial candidate
      candidateSize hprimitive hinitial hscan with
    ⟨scalar, hscalar, hextraction⟩
  have hmappedExtraction := congrArg
    (Polynomial.map (Int.castRingHom (ZMod base))) hextraction
  simp only [Polynomial.map_mul, Polynomial.map_C] at hmappedExtraction
  have hscalarMod : IsUnit (scalar : ZMod base) :=
    hscalar.map (Int.castRingHom (ZMod base))
  have hconstantUnit : IsUnit
      (Polynomial.C (scalar : ZMod base)) :=
    Polynomial.isUnit_C.mpr hscalarMod
  have hsourceExtraction : Associated
      (Polynomial.map (Int.castRingHom (ZMod base))
        (SparsePolyZZ.toPoly fStar))
      (Refinement.StrictHensel.toPolyMod base factor *
        Polynomial.map (Int.castRingHom (ZMod base))
          (SparsePolyZZ.toPoly quotientPrimitive)) := by
    have hwithConstant := Associated.of_eq hmappedExtraction
    exact hwithConstant.trans
      ((associated_isUnit_mul_left_iff hconstantUnit).mpr (Associated.refl _))
  have hpartition := removeCombination_toPolyMod_product_partition base
    candidate activeLifted active' hlegalSize hremove
  have hcombined : Associated
      (Refinement.StrictHensel.toPolyMod base factor *
        Polynomial.map (Int.castRingHom (ZMod base))
          (SparsePolyZZ.toPoly quotientPrimitive))
      (((selectSourceIndices activeLifted.toList candidate.toList).map
          (Refinement.StrictHensel.toPolyMod base)).prod *
        (active'.toList.map
          (Refinement.StrictHensel.toPolyMod base)).prod) :=
    hsourceExtraction.symm.trans
      (hlive.trans (Associated.of_eq hpartition.symm))
  have hbound : ∀ position (hposition : position < candidate.size),
      candidate[position] < activeLifted.size := by
    intro position hposition
    simpa [getElem!_pos candidate position hposition] using
      hlegalSize.2.2 position hposition
  have hselectedNe := selectedSourceProduct_ne_zero_of_irreducible base
    activeLifted candidate hbound hirreducible
  have hfactorNe : Refinement.StrictHensel.toPolyMod base factor ≠ 0 := by
    intro hzero
    apply hselectedNe
    exact hfactorCertificate.1.eq_zero_iff.mp hzero
  exact Associated.of_mul_left hcombined hfactorCertificate.1 hfactorNe

private noncomputable def factorArrayProduct (factors : Array SparsePolyZZ) :
  Polynomial Int :=
  (factors.toList.map SparsePolyZZ.toPoly).prod

private noncomputable def vanHoeijOutputProduct
    (output : SparsePolyZZ × Array SparsePolyZZ) : Polynomial Int :=
  if output.1.isEmpty then factorArrayProduct output.2
  else SparsePolyZZ.toPoly output.1 * factorArrayProduct output.2

private theorem appendFallback_product
    (fallback result output : Array SparsePolyZZ)
    (hrun : Generated.StrictRecombine.appendFallback fallback result =
      .ok output) :
    factorArrayProduct output =
      factorArrayProduct result * factorArrayProduct fallback := by
  unfold Generated.StrictRecombine.appendFallback at hrun
  rw [appendFallbackLoop_refines] at hrun
  have hout := Except.ok.inj hrun
  subst output
  simp [factorArrayProduct]

/-- Degree sorting at either concrete recombination exit changes only factor
order and therefore preserves the represented product exactly. -/
private theorem factorArrayProduct_sortFactorsByDegree
    (factors : Array SparsePolyZZ) :
    factorArrayProduct
        (Generated.StrictRecombine.sortFactorsByDegree factors) =
      factorArrayProduct factors := by
  unfold factorArrayProduct Generated.StrictRecombine.sortFactorsByDegree
  simpa using ((List.mergeSort_perm factors.toList fun left right =>
    left[0]!.1.deg < right[0]!.1.deg).map SparsePolyZZ.toPoly).prod_eq

/-- Exact product effect of the common source finishing block, including its
conditional append of the remaining positive-degree factor. -/
private theorem factorArrayProduct_finishZassenhaus
    (fStar : SparsePolyZZ) (result : Array SparsePolyZZ) :
    factorArrayProduct
        (Generated.StrictRecombine.finishZassenhaus fStar result) =
      if hnonempty : 0 < fStar.size then
        if 0 < fStar[0].1.deg then
          factorArrayProduct result * SparsePolyZZ.toPoly fStar
        else factorArrayProduct result
      else factorArrayProduct result := by
  unfold Generated.StrictRecombine.finishZassenhaus
  rw [factorArrayProduct_sortFactorsByDegree]
  split
  next hnonempty =>
    split
    next hdegree => simp [factorArrayProduct]
    next hdegree => rfl
  next hnonempty => rfl

/-- A canonical sparse integer polynomial whose concrete front degree is zero
contains exactly that one constant term. -/
private theorem sparsePolyZZ_toPoly_eq_C_of_front_degree_zero
    (f : SparsePolyZZ)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hnonempty : 0 < f.size) (hdegree : f[0].1.deg = 0) :
    SparsePolyZZ.toPoly f = Polynomial.C f[0].2 := by
  cases hlist : f.toList with
  | nil =>
      have : f = #[] := by simpa using hlist
      subst f
      simp at hnonempty
  | cons head tail =>
      have hhead : head = f[0] := by
        have hoption := congrArg (fun list => list[0]?) hlist
        have hheadOption : f.toList[0]? = some f[0] := by simp
        exact Option.some.inj (hoption.symm.trans hheadOption)
      cases tail with
      | nil =>
          unfold SparsePolyZZ.toPoly
          rw [hlist]
          simp [hhead, hdegree]
      | cons second rest =>
          have hchain := hcanonical.1
          rw [hlist, List.isChain_cons_cons] at hchain
          rw [hhead, hdegree] at hchain
          omega

/-- The common concrete finish block represents the whole live product up to
a unit.  In the source branch that drops a degree-zero `fStar`, canonicality
identifies it as a constant and primitivity proves that constant is a unit. -/
private theorem finishZassenhaus_product_associated
    (fStar : SparsePolyZZ) (result : Array SparsePolyZZ)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical fStar)
    (hprimitive :
      (SparsePolyZZ.toPoly fStar * factorArrayProduct result).IsPrimitive) :
    Associated (SparsePolyZZ.toPoly fStar * factorArrayProduct result)
      (factorArrayProduct
        (Generated.StrictRecombine.finishZassenhaus fStar result)) := by
  rw [factorArrayProduct_finishZassenhaus]
  by_cases hnonempty : 0 < fStar.size
  · rw [dif_pos hnonempty]
    by_cases hdegree : 0 < fStar[0].1.deg
    · rw [if_pos hdegree]
      exact Associated.of_eq (mul_comm _ _)
    · rw [if_neg hdegree]
      have hdegreeZero : fStar[0].1.deg = 0 := by omega
      have hconstant := sparsePolyZZ_toPoly_eq_C_of_front_degree_zero
        fStar hcanonical hnonempty hdegreeZero
      have hcontent := Polynomial.isPrimitive_iff_content_eq_one.mp hprimitive
      rw [hconstant, Polynomial.content_mul, Polynomial.content_C] at hcontent
      have hnormalizeUnit : IsUnit (normalize fStar[0].2) :=
        IsUnit.of_mul_eq_one _ hcontent
      have hcoefficientUnit : IsUnit fStar[0].2 :=
        (associated_normalize fStar[0].2).isUnit_iff.mpr hnormalizeUnit
      have hconstantUnit : IsUnit (SparsePolyZZ.toPoly fStar) := by
        rw [hconstant]
        exact Polynomial.isUnit_C.mpr hcoefficientUnit
      exact (associated_isUnit_mul_left_iff hconstantUnit).mpr
        (Associated.refl _)
  · rw [dif_neg hnonempty]
    have hempty : fStar = #[] := by
      have hsize : fStar.size = 0 := Nat.eq_zero_of_not_pos hnonempty
      apply Array.ext (by simpa using hsize)
      intro index hindex hindexEmpty
      simp at hindexEmpty
    subst fStar
    have hcontent := Polynomial.isPrimitive_iff_content_eq_one.mp hprimitive
    simp [SparsePolyZZ.toPoly] at hcontent

private theorem isPrimitive_left_of_mul {left right : Polynomial Int}
    (hprimitive : (left * right).IsPrimitive) : left.IsPrimitive := by
  have hcontent := Polynomial.isPrimitive_iff_content_eq_one.mp hprimitive
  rw [Polynomial.content_mul] at hcontent
  have hunit : IsUnit left.content := IsUnit.of_mul_eq_one _ hcontent
  apply Polynomial.isPrimitive_iff_content_eq_one.mpr
  calc
    left.content = normalize left.content := Polynomial.normalize_content.symm
    _ = 1 := normalize_eq_one.mpr hunit

private theorem isPrimitive_of_unit_scalar_product
    {before after : Polynomial Int} {scalar : Int}
    (hscalar : IsUnit scalar) (hbefore : before.IsPrimitive)
    (heq : before = Polynomial.C scalar * after) : after.IsPrimitive := by
  have hcontent := congrArg Polynomial.content heq
  rw [hbefore.content_eq_one, Polynomial.content_C_mul,
    normalize_eq_one.mpr hscalar, one_mul] at hcontent
  exact Polynomial.isPrimitive_iff_content_eq_one.mpr hcontent.symm

/-- Under the live Hensel invariants, the complete generated Zassenhaus outer
loop cannot take a raw-fault branch.  Each fixed-size scan is executed by its
concrete well-founded recursion; a returned legal candidate is physically
removed, and the quotient/remnant invariants are obtained from that same run. -/
theorem zassenhausLoop_complete
    (modulus base : Nat) [Fact (Nat.Prime base)]
    (active : Array SparsePolyZZ) (fStar : SparsePolyZZ)
    (result : Array SparsePolyZZ) (subsetSize : Nat)
    (hmodulus : 0 < modulus) (hbase : 0 < base)
    (hdivides : base ∣ modulus)
    (hsubsetPositive : 0 < subsetSize)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical fStar)
    (hnonempty : 0 < fStar.size)
    (hprimitive : (SparsePolyZZ.toPoly fStar).IsPrimitive)
    (hleading : (fStar[0].2 : ZMod base) ≠ 0)
    (hactiveFits : active.size ≤ 2 ^ 31)
    (hirreducible : ∀ index (hindex : index < active.size),
      Irreducible (Refinement.StrictHensel.toPolyMod base active[index])) :
    ∃ output,
      Generated.StrictRecombine.zassenhausLoop concreteZassenhausTermination
        (modulus : ZZ) active fStar result subsetSize hsubsetPositive =
          .ok output := by
  rw [Generated.StrictRecombine.zassenhausLoop]
  split
  next hcontinue =>
    let count := subsetSize
    let initial := Generated.StrictRecombine.initialCombination count
    have hfits : count ≤ active.size := by omega
    let hinitial := initialCombination_legal active.size count hfits
    have hscanComplete := scanZassenhausCombinations_complete fStar active
      modulus base initial hmodulus hbase hdivides hcanonical hnonempty
      hsubsetPositive hactiveFits hleading hirreducible hinitial
    rcases hscanComplete with ⟨scanResult, hscan⟩
    cases scanResult with
    | exhausted =>
        simp only [concreteZassenhausTermination, count, initial, hinitial,
          hscan]
        exact zassenhausLoop_complete modulus base active fStar result
          (subsetSize + 1) hmodulus hbase hdivides (by omega) hcanonical
          hnonempty hprimitive hleading hactiveFits hirreducible
    | extracted factor quotient candidate candidateSize =>
        have hcandidate : LegalCombination active.size count candidate :=
          concreteScan_extracted_legal fStar factor quotient active
            (modulus : ZZ) candidate candidateSize hfits hscan
        rcases removeCombination_succeeds candidate active
            (by simpa [count, candidateSize] using hcandidate) with
          ⟨active', hremove⟩
        have hfactorQuotient :=
          scanZassenhausCombinations_extracted_canonical_primitive
            (concreteCombinationTermination active.size count) fStar factor
            quotient active (modulus : ZZ) initial candidate candidateSize
            hcanonical hnonempty hinitial hscan
        have hleadingPoly :
            ((SparsePolyZZ.toPoly fStar).leadingCoeff : ZMod base) ≠ 0 := by
          rw [sparsePolyZZ_leadingCoeff_eq_head fStar hcanonical hnonempty]
          exact hleading
        have hmodCertificate :=
          scanZassenhausCombinations_extracted_mod_certificate fStar factor
            quotient active modulus base initial candidate candidateSize
            hmodulus hbase hdivides hactiveFits hcanonical hnonempty hprimitive
            hleadingPoly hirreducible hinitial hscan
        have hquotientNe : SparsePolyZZ.toPoly quotient ≠ 0 := by
          intro hzero
          apply hmodCertificate.2
          rw [hzero]
          simp
        have hquotientNonempty : 0 < quotient.size := by
          by_contra hnot
          have hempty : quotient = #[] := Array.size_eq_zero_iff.mp
            (Nat.eq_zero_of_not_pos hnot)
          apply hquotientNe
          simp [hempty, SparsePolyZZ.toPoly]
        have hquotientLeading : (quotient[0].2 : ZMod base) ≠ 0 := by
          rw [← sparsePolyZZ_leadingCoeff_eq_head quotient
            hfactorQuotient.2.1 hquotientNonempty]
          exact hmodCertificate.2
        have hirreducible' : ∀ index (hindex : index < active'.size),
            Irreducible
              (Refinement.StrictHensel.toPolyMod base active'[index]) :=
          removeCombination_preserves_pointwise candidate active active'
            (fun poly => Irreducible
              (Refinement.StrictHensel.toPolyMod base poly))
            hirreducible hremove
        have hactiveFits' : active'.size ≤ 2 ^ 31 := by
          have hdecrease := concreteZassenhausTermination.removal_decreases
            active candidate active' (by rw [candidateSize]; exact hsubsetPositive)
            hremove
          omega
        simp only [concreteZassenhausTermination, count, initial, hinitial,
          hscan]
        rw [hremove]
        exact zassenhausLoop_complete modulus base active' quotient
          (result.push factor) 1 hmodulus hbase hdivides (by omega)
          hfactorQuotient.2.1 hquotientNonempty hfactorQuotient.2.2
          hquotientLeading hactiveFits' hirreducible'
  next hcontinue =>
    exact ⟨Generated.StrictRecombine.finishZassenhaus fStar result, rfl⟩
termination_by (active.size, active.size + 1 - subsetSize)
decreasing_by
  · exact Prod.Lex.right _ (by omega)
  · exact Prod.Lex.left _ _
      (concreteZassenhausTermination.removal_decreases active candidate active'
        (by rw [candidateSize]; exact hsubsetPositive) hremove)

/-- The complete source-shaped Zassenhaus outer loop preserves the live
product up to a unit.  Successful extraction recursively installs the exact
primitive quotient returned by the candidate scan; exhaustion and subset-size
advancement follow the generated control flow unchanged. -/
theorem zassenhausLoop_product_associated
    (termination : Generated.StrictRecombine.ZassenhausTermination)
    (modulus : ZZ) (active : Array SparsePolyZZ) (fStar : SparsePolyZZ)
    (result : Array SparsePolyZZ) (subsetSize : Nat)
    (hsubsetPositive : 0 < subsetSize)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical fStar)
    (hprimitive :
      (SparsePolyZZ.toPoly fStar * factorArrayProduct result).IsPrimitive)
    (output : Array SparsePolyZZ)
    (hrun : Generated.StrictRecombine.zassenhausLoop termination modulus active
      fStar result subsetSize hsubsetPositive = .ok output) :
    Associated (SparsePolyZZ.toPoly fStar * factorArrayProduct result)
      (factorArrayProduct output) := by
  rw [Generated.StrictRecombine.zassenhausLoop] at hrun
  split at hrun
  next hcontinue =>
    let combinationTermination :=
      termination.combinations active.size subsetSize
    let initial := Generated.StrictRecombine.initialCombination subsetSize
    let hinitial := termination.initial_valid active.size subsetSize (by omega)
    cases hscan : Generated.StrictRecombine.scanZassenhausCombinations
        combinationTermination fStar active modulus initial hinitial with
    | error fault => simp [combinationTermination, initial, hinitial, hscan] at hrun
    | ok scanResult =>
      cases scanResult with
      | exhausted =>
        simp only [combinationTermination, initial, hinitial, hscan] at hrun
        exact zassenhausLoop_product_associated termination modulus active fStar
          result (subsetSize + 1) (by omega) hcanonical hprimitive output hrun
      | extracted factor quotient candidate candidateSize =>
        simp only [combinationTermination, initial, hinitial, hscan] at hrun
        cases hremove : Generated.StrictRecombine.removeCombination candidate active with
        | error fault =>
          rw [hremove] at hrun
          contradiction
        | ok active' =>
          rw [hremove] at hrun
          have hfStarPrimitive : (SparsePolyZZ.toPoly fStar).IsPrimitive :=
            isPrimitive_left_of_mul hprimitive
          rcases scanZassenhausCombinations_extracted_unit_scalar
              combinationTermination fStar factor quotient active modulus initial
              candidate candidateSize hfStarPrimitive hinitial hscan with
            ⟨scalar, hscalar, hextraction⟩
          have hquotientCanonical :=
            scanZassenhausCombinations_extracted_quotient_canonical
              combinationTermination fStar factor quotient active modulus initial
              candidate candidateSize hcanonical hinitial hscan
          have hlive : SparsePolyZZ.toPoly fStar * factorArrayProduct result =
              Polynomial.C scalar *
                (SparsePolyZZ.toPoly quotient *
                  factorArrayProduct (result.push factor)) := by
            rw [hextraction]
            simp [factorArrayProduct]
            ring
          have hnextPrimitive :
              (SparsePolyZZ.toPoly quotient *
                factorArrayProduct (result.push factor)).IsPrimitive :=
            isPrimitive_of_unit_scalar_product hscalar hprimitive hlive
          have htail := zassenhausLoop_product_associated termination modulus
            active' quotient (result.push factor) 1 (by omega)
            hquotientCanonical hnextPrimitive output hrun
          have hconstantUnit : IsUnit (Polynomial.C scalar : Polynomial Int) :=
            Polynomial.isUnit_C.mpr hscalar
          exact (Associated.of_eq hlive).trans
            ((associated_isUnit_mul_left_iff hconstantUnit).mpr htail)
  next hcontinue =>
    have hout := Except.ok.inj hrun
    subst output
    exact finishZassenhaus_product_associated fStar result hcanonical hprimitive
termination_by (active.size, active.size + 1 - subsetSize)
decreasing_by
  · exact Prod.Lex.right _ (by omega)
  · exact Prod.Lex.left _ _
      (termination.removal_decreases active candidate active'
        (by rw [candidateSize]; exact hsubsetPositive) hremove)

theorem zassenhausRecombine_product_associated
    (termination : Generated.StrictRecombine.ZassenhausTermination)
    (f : SparsePolyZZ) (lifted output : Array SparsePolyZZ) (modulus : ZZ)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hprimitive : (SparsePolyZZ.toPoly f).IsPrimitive)
    (hrun : Generated.StrictRecombine.zassenhausRecombine termination f lifted
      modulus = .ok output) :
    Associated (SparsePolyZZ.toPoly f) (factorArrayProduct output) := by
  unfold Generated.StrictRecombine.zassenhausRecombine at hrun
  split at hrun
  next hlifted =>
    split at hrun
    next hnonempty =>
      split at hrun
      next hdegree =>
        have hout := Except.ok.inj hrun
        subst output
        simp [factorArrayProduct]
      next hdegree =>
        have hout := Except.ok.inj hrun
        subst output
        have hfinish := finishZassenhaus_product_associated f #[] hcanonical
          (by simpa [factorArrayProduct] using hprimitive)
        rw [factorArrayProduct_finishZassenhaus, dif_pos hnonempty,
          if_neg hdegree] at hfinish
        simpa [factorArrayProduct] using hfinish
    next hnonempty =>
      have hout := Except.ok.inj hrun
      subst output
      have hfinish := finishZassenhaus_product_associated f #[] hcanonical
        (by simpa [factorArrayProduct] using hprimitive)
      rw [factorArrayProduct_finishZassenhaus, dif_neg hnonempty] at hfinish
      simpa [factorArrayProduct] using hfinish
  next hlifted =>
    have hloop := zassenhausLoop_product_associated termination modulus lifted f
      #[] 1 (by omega) hcanonical
      (by simpa [factorArrayProduct] using hprimitive) output hrun
    simpa [factorArrayProduct] using hloop

/-- Public product form for consumers of the generated Zassenhaus entry. -/
theorem zassenhausRecombine_toPoly_product_associated
    (termination : Generated.StrictRecombine.ZassenhausTermination)
    (f : SparsePolyZZ) (lifted output : Array SparsePolyZZ) (modulus : ZZ)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hprimitive : (SparsePolyZZ.toPoly f).IsPrimitive)
    (hrun : Generated.StrictRecombine.zassenhausRecombine termination f lifted
      modulus = .ok output) :
    Associated (SparsePolyZZ.toPoly f)
      (output.toList.map SparsePolyZZ.toPoly).prod := by
  simpa [factorArrayProduct] using
    zassenhausRecombine_product_associated termination f lifted output modulus
      hcanonical hprimitive hrun

theorem validateCandidatesLoop_product
    (ops : Generated.StrictRecombine.CandidateValidationRawOps)
    (candidates : Array (Array Int32)) (candidateIndex : Nat)
    (activeLifted : Array SparsePolyZZ) (modulus : ZZ)
    (fStar fStar' : SparsePolyZZ) (result result' : Array SparsePolyZZ)
    (consumed consumed' : Array Bool) (remaining : Nat)
    (hrun : Generated.StrictRecombine.validateCandidatesLoop ops candidates
      candidateIndex activeLifted modulus fStar result consumed remaining =
        .ok (fStar', result', consumed')) :
    ∃ scalar : Int,
      SparsePolyZZ.toPoly fStar * factorArrayProduct result =
        Polynomial.C scalar *
          (SparsePolyZZ.toPoly fStar' * factorArrayProduct result') := by
  induction hmeasure : candidates.size - candidateIndex using Nat.strong_induction_on
      generalizing candidateIndex fStar result consumed remaining fStar' result' consumed' with
  | h measure ih =>
      rw [Generated.StrictRecombine.validateCandidatesLoop] at hrun
      split at hrun
      next hcandidates =>
        dsimp at hrun
        split at hrun
        next hempty =>
          exact ih (candidates.size - (candidateIndex + 1)) (by omega)
            (candidateIndex := candidateIndex + 1) (fStar := fStar)
            (result := result) (consumed := consumed) (remaining := remaining)
            (fStar' := fStar') (result' := result') (consumed' := consumed') hrun rfl
        next hempty =>
          split at hrun
          next htrivial =>
            exact ih (candidates.size - (candidateIndex + 1)) (by omega)
              (candidateIndex := candidateIndex + 1) (fStar := fStar)
              (result := result) (consumed := consumed) (remaining := remaining)
              (fStar' := fStar') (result' := result') (consumed' := consumed') hrun rfl
          next hnontrivial =>
            cases havailable : Generated.StrictRecombine.candidateAvailable
                candidates[candidateIndex] consumed with
            | error fault => simp [havailable] at hrun
            | ok available =>
              cases available with
              | false =>
                simp only [havailable] at hrun
                exact ih (candidates.size - (candidateIndex + 1)) (by omega)
                  (candidateIndex := candidateIndex + 1) (fStar := fStar)
                  (result := result) (consumed := consumed) (remaining := remaining)
                  (fStar' := fStar') (result' := result')
                  (consumed' := consumed') hrun rfl
              | true =>
                simp only [havailable] at hrun
                split at hrun
                next hfstar =>
                  cases hproduct : Generated.StrictRecombine.trialProductLoop
                      ops.product candidates[candidateIndex] activeLifted modulus 0
                      #[(⟨0⟩, fStar[0].2)] with
                  | error fault => simp [hproduct] at hrun
                  | ok product =>
                    simp only [hproduct] at hrun
                    cases hsymmetric : Generated.StrictRecombine.symmetricModRaw
                        product modulus with
                    | error fault => simp [hsymmetric] at hrun
                    | ok symmetric =>
                      simp only [hsymmetric] at hrun
                      cases hprimitive : Generated.StrictRecombine.primitiveRaw symmetric with
                      | error fault => simp [hprimitive] at hrun
                      | ok primitiveResult =>
                        rcases primitiveResult with ⟨content, factor⟩
                        simp only [hprimitive] at hrun
                        cases hdivmod : Generated.StrictRecombine.exactDivmodRaw
                            fStar factor with
                        | error fault => simp [hdivmod] at hrun
                        | ok divResult =>
                          rcases divResult with ⟨quotient, remainder⟩
                          simp only [hdivmod] at hrun
                          by_cases hremainder : remainder.isEmpty = true
                          · simp only [hremainder, if_true] at hrun
                            have hremainderEmpty : remainder = #[] :=
                              Array.isEmpty_iff.mp hremainder
                            subst remainder
                            cases hquotientPrimitive :
                                Generated.StrictRecombine.primitiveRaw quotient with
                            | error fault => simp [hquotientPrimitive] at hrun
                            | ok quotientResult =>
                              rcases quotientResult with ⟨quotientContent,
                                quotientPrimitive⟩
                              simp only [hquotientPrimitive] at hrun
                              cases hmark : Generated.StrictRecombine.markConsumedLoop
                                  candidates[candidateIndex] 0 consumed with
                              | error fault => simp [hmark] at hrun
                              | ok consumedNext =>
                                simp only [hmark] at hrun
                                rcases ih (candidates.size - (candidateIndex + 1))
                                    (by omega) (candidateIndex := candidateIndex + 1)
                                    (fStar := quotientPrimitive)
                                    (result := result.push factor)
                                    (consumed := consumedNext)
                                    (remaining := remaining - candidates[candidateIndex].size)
                                    (fStar' := fStar') (result' := result')
                                    (consumed' := consumed') hrun rfl with
                                  ⟨tailScalar, htail⟩
                                refine ⟨quotientContent * tailScalar, ?_⟩
                                rw [successfulTrialExtraction_toPoly fStar factor
                                  quotient quotientPrimitive quotientContent hdivmod
                                  hquotientPrimitive]
                                simp [factorArrayProduct] at htail ⊢
                                calc
                                  _ = (quotientContent : Polynomial Int) *
                                      (quotientPrimitive.toPoly *
                                        ((result.toList.map SparsePolyZZ.toPoly).prod *
                                          factor.toPoly)) := by ring
                                  _ = (quotientContent : Polynomial Int) *
                                      ((tailScalar : Polynomial Int) *
                                        (fStar'.toPoly *
                                          (result'.toList.map SparsePolyZZ.toPoly).prod)) := by
                                    rw [htail]
                                  _ = _ := by ring
                          · simp only [hremainder] at hrun
                            exact ih (candidates.size - (candidateIndex + 1))
                              (by omega) (candidateIndex := candidateIndex + 1)
                              (fStar := fStar) (result := result)
                              (consumed := consumed) (remaining := remaining)
                              (fStar' := fStar') (result' := result')
                              (consumed' := consumed') hrun rfl
                next hfstar => contradiction
      next hcandidates =>
        have hout := Except.ok.inj hrun
        cases hout
        exact ⟨1, by simp⟩

/-- Candidate validation preserves the canonical sparse representation of its
live remainder through every concrete successful extraction. -/
theorem validateCandidatesLoop_fStar_canonical
    (ops : Generated.StrictRecombine.CandidateValidationRawOps)
    (candidates : Array (Array Int32)) (candidateIndex : Nat)
    (activeLifted : Array SparsePolyZZ) (modulus : ZZ)
    (fStar fStar' : SparsePolyZZ) (result result' : Array SparsePolyZZ)
    (consumed consumed' : Array Bool) (remaining : Nat)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical fStar)
    (hrun : Generated.StrictRecombine.validateCandidatesLoop ops candidates
      candidateIndex activeLifted modulus fStar result consumed remaining =
        .ok (fStar', result', consumed')) :
    StrictPolynomialMod.SparsePolyZZCanonical fStar' := by
  induction hmeasure : candidates.size - candidateIndex using Nat.strong_induction_on
      generalizing candidateIndex fStar result consumed remaining fStar' result' consumed' with
  | h measure ih =>
      rw [Generated.StrictRecombine.validateCandidatesLoop] at hrun
      split at hrun
      next hcandidates =>
        dsimp at hrun
        split at hrun
        next hempty =>
          exact ih (candidates.size - (candidateIndex + 1)) (by omega)
            (candidateIndex := candidateIndex + 1) (fStar := fStar)
            (result := result) (consumed := consumed) (remaining := remaining)
            (fStar' := fStar') (result' := result') (consumed' := consumed')
            hcanonical hrun rfl
        next hempty =>
          split at hrun
          next htrivial =>
            exact ih (candidates.size - (candidateIndex + 1)) (by omega)
              (candidateIndex := candidateIndex + 1) (fStar := fStar)
              (result := result) (consumed := consumed) (remaining := remaining)
              (fStar' := fStar') (result' := result') (consumed' := consumed')
              hcanonical hrun rfl
          next hnontrivial =>
            cases havailable : Generated.StrictRecombine.candidateAvailable
                candidates[candidateIndex] consumed with
            | error fault => simp [havailable] at hrun
            | ok available =>
              cases available with
              | false =>
                simp only [havailable] at hrun
                exact ih (candidates.size - (candidateIndex + 1)) (by omega)
                  (candidateIndex := candidateIndex + 1) (fStar := fStar)
                  (result := result) (consumed := consumed) (remaining := remaining)
                  (fStar' := fStar') (result' := result')
                  (consumed' := consumed') hcanonical hrun rfl
              | true =>
                simp only [havailable] at hrun
                split at hrun
                next hfstar =>
                  cases hproduct : Generated.StrictRecombine.trialProductLoop
                      ops.product candidates[candidateIndex] activeLifted modulus 0
                      #[(⟨0⟩, fStar[0].2)] with
                  | error fault => simp [hproduct] at hrun
                  | ok product =>
                    simp only [hproduct] at hrun
                    cases hsymmetric : Generated.StrictRecombine.symmetricModRaw
                        product modulus with
                    | error fault => simp [hsymmetric] at hrun
                    | ok symmetric =>
                      simp only [hsymmetric] at hrun
                      cases hprimitive : Generated.StrictRecombine.primitiveRaw
                          symmetric with
                      | error fault => simp [hprimitive] at hrun
                      | ok primitiveResult =>
                        rcases primitiveResult with ⟨content, factor⟩
                        simp only [hprimitive] at hrun
                        cases hdivmod : Generated.StrictRecombine.exactDivmodRaw
                            fStar factor with
                        | error fault => simp [hdivmod] at hrun
                        | ok divResult =>
                          rcases divResult with ⟨quotient, remainder⟩
                          simp only [hdivmod] at hrun
                          by_cases hremainder : remainder.isEmpty = true
                          · simp only [hremainder, if_true] at hrun
                            have hremainderEmpty : remainder = #[] :=
                              Array.isEmpty_iff.mp hremainder
                            subst remainder
                            cases hquotientPrimitive :
                                Generated.StrictRecombine.primitiveRaw quotient with
                            | error fault => simp [hquotientPrimitive] at hrun
                            | ok quotientResult =>
                              rcases quotientResult with ⟨quotientContent,
                                quotientPrimitive⟩
                              simp only [hquotientPrimitive] at hrun
                              cases hmark : Generated.StrictRecombine.markConsumedLoop
                                  candidates[candidateIndex] 0 consumed with
                              | error fault => simp [hmark] at hrun
                              | ok consumedNext =>
                                simp only [hmark] at hrun
                                have hquotientCanonical :=
                                  exactDivmodRaw_quotient_canonical fStar factor
                                    quotient #[] hcanonical.2 hdivmod
                                have hprimitiveCanonical :=
                                  primitiveRaw_canonical quotient quotientPrimitive
                                    quotientContent hquotientCanonical
                                    hquotientPrimitive
                                exact ih
                                  (candidates.size - (candidateIndex + 1))
                                  (by omega)
                                  (candidateIndex := candidateIndex + 1)
                                  (fStar := quotientPrimitive)
                                  (result := result.push factor)
                                  (consumed := consumedNext)
                                  (remaining := remaining -
                                    candidates[candidateIndex].size)
                                  (fStar' := fStar') (result' := result')
                                  (consumed' := consumed') hprimitiveCanonical
                                  hrun rfl
                          · simp only [hremainder] at hrun
                            exact ih (candidates.size - (candidateIndex + 1))
                              (by omega) (candidateIndex := candidateIndex + 1)
                              (fStar := fStar) (result := result)
                              (consumed := consumed) (remaining := remaining)
                              (fStar' := fStar') (result' := result')
                              (consumed' := consumed') hcanonical hrun rfl
                next hfstar => contradiction
      next hcandidates =>
        have hout := Except.ok.inj hrun
        cases hout
        exact hcanonical

theorem validateCandidates_product
    (ops : Generated.StrictRecombine.CandidateValidationRawOps)
    (fStar : SparsePolyZZ) (activeLifted : Array SparsePolyZZ) (modulus : ZZ)
    (candidates : Array (Array Int32)) (result : Array SparsePolyZZ)
    (fStar' : SparsePolyZZ) (result' : Array SparsePolyZZ)
    (consumed : Array Bool)
    (hrun : Generated.StrictRecombine.validateCandidates ops fStar activeLifted
      modulus candidates result = .ok (fStar', result', consumed)) :
    ∃ scalar : Int,
      SparsePolyZZ.toPoly fStar * factorArrayProduct result =
        Polynomial.C scalar *
          (SparsePolyZZ.toPoly fStar' * factorArrayProduct result') := by
  unfold Generated.StrictRecombine.validateCandidates at hrun
  exact validateCandidatesLoop_product ops candidates 0 activeLifted modulus
    fStar fStar' result result' (Array.replicate activeLifted.size false)
    consumed activeLifted.size hrun

theorem validateCandidates_fStar_canonical
    (ops : Generated.StrictRecombine.CandidateValidationRawOps)
    (fStar : SparsePolyZZ) (activeLifted : Array SparsePolyZZ) (modulus : ZZ)
    (candidates : Array (Array Int32)) (result : Array SparsePolyZZ)
    (fStar' : SparsePolyZZ) (result' : Array SparsePolyZZ)
    (consumed : Array Bool)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical fStar)
    (hrun : Generated.StrictRecombine.validateCandidates ops fStar activeLifted
      modulus candidates result = .ok (fStar', result', consumed)) :
    StrictPolynomialMod.SparsePolyZZCanonical fStar' := by
  unfold Generated.StrictRecombine.validateCandidates at hrun
  exact validateCandidatesLoop_fStar_canonical ops candidates 0 activeLifted
    modulus fStar fStar' result result'
    (Array.replicate activeLifted.size false) consumed activeLifted.size
    hcanonical hrun

/-- A successful concrete validation run cannot hide a non-unit integer
content factor when its incoming accumulated product is primitive.  This is
the content bridge needed to turn `validateCandidates_product` into an
associated-product statement, without assuming that the generated primitive
normalization returned a semantically convenient result. -/
theorem validateCandidates_product_unit_scalar
    (ops : Generated.StrictRecombine.CandidateValidationRawOps)
    (fStar : SparsePolyZZ) (activeLifted : Array SparsePolyZZ) (modulus : ZZ)
    (candidates : Array (Array Int32)) (result : Array SparsePolyZZ)
    (fStar' : SparsePolyZZ) (result' : Array SparsePolyZZ)
    (consumed : Array Bool)
    (hprimitive :
      (SparsePolyZZ.toPoly fStar * factorArrayProduct result).IsPrimitive)
    (hrun : Generated.StrictRecombine.validateCandidates ops fStar activeLifted
      modulus candidates result = .ok (fStar', result', consumed)) :
    ∃ scalar : Int, IsUnit scalar ∧
      SparsePolyZZ.toPoly fStar * factorArrayProduct result =
        Polynomial.C scalar *
          (SparsePolyZZ.toPoly fStar' * factorArrayProduct result') := by
  rcases validateCandidates_product ops fStar activeLifted modulus candidates
      result fStar' result' consumed hrun with ⟨scalar, hproduct⟩
  have hcontent := congrArg Polynomial.content hproduct
  rw [hprimitive.content_eq_one, Polynomial.content_C_mul] at hcontent
  have hnormalizeUnit : IsUnit (normalize scalar) :=
    IsUnit.of_mul_eq_one _ hcontent.symm
  have hscalarUnit : IsUnit scalar :=
    (associated_normalize scalar).isUnit_iff.mpr hnormalizeUnit
  exact ⟨scalar, hscalarUnit, hproduct⟩

/-- The concrete extraction loop preserves primitivity of the whole live
product.  This is the induction invariant used by the source-shaped van-Hoeij
loop after it removes the consumed lifted factors. -/
theorem validateCandidates_preserves_primitive
    (ops : Generated.StrictRecombine.CandidateValidationRawOps)
    (fStar : SparsePolyZZ) (activeLifted : Array SparsePolyZZ) (modulus : ZZ)
    (candidates : Array (Array Int32)) (result : Array SparsePolyZZ)
    (fStar' : SparsePolyZZ) (result' : Array SparsePolyZZ)
    (consumed : Array Bool)
    (hprimitive :
      (SparsePolyZZ.toPoly fStar * factorArrayProduct result).IsPrimitive)
    (hrun : Generated.StrictRecombine.validateCandidates ops fStar activeLifted
      modulus candidates result = .ok (fStar', result', consumed)) :
    (SparsePolyZZ.toPoly fStar' * factorArrayProduct result').IsPrimitive := by
  rcases validateCandidates_product_unit_scalar ops fStar activeLifted modulus
      candidates result fStar' result' consumed hprimitive hrun with
    ⟨scalar, hscalarUnit, hproduct⟩
  have hcontent := congrArg Polynomial.content hproduct
  rw [hprimitive.content_eq_one, Polynomial.content_C_mul,
    normalize_eq_one.mpr hscalarUnit, one_mul] at hcontent
  exact Polynomial.isPrimitive_iff_content_eq_one.mpr hcontent.symm

private theorem natList_all_eq_one_of_pos_of_sum_eq_length
    (values : List Nat) (hpos : ∀ value ∈ values, 0 < value)
    (hsum : values.sum = values.length) :
    ∀ value ∈ values, value = 1 := by
  induction values with
  | nil => simp
  | cons head tail ih =>
      have hhead : 0 < head := hpos head (by simp)
      have htailPos : ∀ value ∈ tail, 0 < value := by
        intro value hvalue
        exact hpos value (by simp [hvalue])
      have lower : ∀ values : List Nat,
          (∀ value ∈ values, 0 < value) → values.length ≤ values.sum := by
        intro values hvalues
        induction values with
        | nil => simp
        | cons next rest restIH =>
            have hnext := hvalues next (by simp)
            have hrest : ∀ value ∈ rest, 0 < value := by
              intro value hvalue
              exact hvalues value (by simp [hvalue])
            have hrestLower := restIH hrest
            simp only [List.length_cons, List.sum_cons]
            omega
      have htailLower : tail.length ≤ tail.sum := lower tail htailPos
      have hheadOne : head = 1 := by
        simp only [List.sum_cons, List.length_cons] at hsum
        omega
      have htailSum : tail.sum = tail.length := by
        simp only [List.sum_cons, List.length_cons, hheadOne] at hsum
        omega
      intro value hvalue
      simp only [List.mem_cons] at hvalue
      rcases hvalue with rfl | hvalue
      · exact hheadOne
      · exact ih htailPos htailSum value hvalue

/-- In a UFD, a successful product decomposition with exactly as many
nonzero nonunit outputs as irreducible source atoms cannot have merged two
atoms into one output: every physical output has exactly one normalized
factor and is therefore irreducible.  This cardinality bridge is used by the
literal van-Hoeij validation/retry path; it does not replace that execution
with a semantic factorization oracle. -/
theorem irreducible_members_of_associated_products_and_equal_length
    {M : Type*} [CommMonoidWithZero M] [NormalizationMonoid M]
    [UniqueFactorizationMonoid M]
    (atoms outputs : List M)
    (hatoms : ∀ atom ∈ atoms, Irreducible atom)
    (houtputsNe : ∀ output ∈ outputs, output ≠ 0)
    (houtputsNonunit : ∀ output ∈ outputs, ¬ IsUnit output)
    (hproduct : Associated outputs.prod atoms.prod)
    (hlength : outputs.length = atoms.length) :
    ∀ output ∈ outputs, Irreducible output := by
  have normalizedFactors_list_prod : ∀ (values : List M),
      (∀ value ∈ values, value ≠ 0) →
      UniqueFactorizationMonoid.normalizedFactors values.prod =
        (values.map UniqueFactorizationMonoid.normalizedFactors).foldr
          (.+.) 0 := by
    intro values hvalues
    induction values with
    | nil => simp
    | cons head tail ih =>
        have hhead := hvalues head (by simp)
        have htail : ∀ value ∈ tail, value ≠ 0 := by
          intro value hvalue
          exact hvalues value (by simp [hvalue])
        letI := nontrivial_of_ne head 0 hhead
        have htailZero : (0 : M) ∉ tail := by
          intro hzero
          exact htail 0 hzero rfl
        rw [List.prod_cons,
          UniqueFactorizationMonoid.normalizedFactors_mul hhead
            (List.prod_ne_zero htailZero), ih htail]
        simp
  have hnormalizedProduct := hproduct.normalizedFactors_eq
  rw [normalizedFactors_list_prod outputs houtputsNe,
    normalizedFactors_list_prod atoms (fun atom hatom =>
      (hatoms atom hatom).ne_zero)] at hnormalizedProduct
  have card_normalizedFactors_fold : ∀ values : List M,
      ((values.map UniqueFactorizationMonoid.normalizedFactors).foldr
          (.+.) 0).card =
        (values.map fun value =>
          (UniqueFactorizationMonoid.normalizedFactors value).card).sum := by
    intro values
    induction values with
    | nil => simp
    | cons head tail ih => simp [ih, Multiset.card_add]
  have hcard :
      (outputs.map (fun output =>
        (UniqueFactorizationMonoid.normalizedFactors output).card)).sum =
        outputs.length := by
    have hcards := congrArg Multiset.card hnormalizedProduct
    rw [card_normalizedFactors_fold outputs,
      card_normalizedFactors_fold atoms] at hcards
    have hatomCards : atoms.map (fun atom =>
        (UniqueFactorizationMonoid.normalizedFactors atom).card) =
        atoms.map (fun _ => 1) := by
      apply List.map_congr_left
      intro atom hatom
      rw [UniqueFactorizationMonoid.normalizedFactors_irreducible
        (hatoms atom hatom)]
      simp
    simpa [hatomCards, hlength] using hcards
  have hpositive : ∀ count ∈ outputs.map (fun output =>
      (UniqueFactorizationMonoid.normalizedFactors output).card),
      0 < count := by
    intro count hcount
    rcases List.mem_map.mp hcount with ⟨output, houtput, rfl⟩
    exact Multiset.card_pos.mpr (ne_of_gt
      ((UniqueFactorizationMonoid.normalizedFactors_pos output
        (houtputsNe output houtput)).2 (houtputsNonunit output houtput)))
  have hone := natList_all_eq_one_of_pos_of_sum_eq_length
    (outputs.map fun output =>
      (UniqueFactorizationMonoid.normalizedFactors output).card)
    hpositive (by simpa using hcard)
  intro output houtput
  have houtputCard :
      (UniqueFactorizationMonoid.normalizedFactors output).card = 1 :=
    hone _ (List.mem_map.mpr ⟨output, houtput, rfl⟩)
  rcases Multiset.card_eq_one.mp houtputCard with ⟨factor, hfactor⟩
  have hfactorMem : factor ∈
      UniqueFactorizationMonoid.normalizedFactors output := by
    rw [hfactor]
    simp
  have hfactorIrreducible :=
    UniqueFactorizationMonoid.irreducible_of_normalized_factor factor
      hfactorMem
  have hassociated : Associated factor output := by
    have hnormalized := UniqueFactorizationMonoid.prod_normalizedFactors
      (houtputsNe output houtput)
    rw [hfactor] at hnormalized
    simpa using hnormalized
  exact hassociated.irreducible hfactorIrreducible

/-- The source-shaped van-Hoeij loop preserves the complete live product up
to a unit across validation extraction, precision retry, and the concrete
Zassenhaus fallback. -/
theorem vanHoeijLoop_product_associated
    (ops : Generated.StrictRecombine.VanHoeijRawOps)
    (termination : Generated.StrictRecombine.VanHoeijTermination ops)
    (lifted : Array SparsePolyZZ) (modulus : ZZ)
    (initial maximum : Nat) (hinitial : 0 < initial)
    (state : Generated.StrictRecombine.VanHoeijState)
    (hstate : Generated.StrictRecombine.VanHoeijStateValid ops state)
    (hnonempty : 0 < state.fStar.size)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical state.fStar)
    (hprimitive :
      (SparsePolyZZ.toPoly state.fStar *
        factorArrayProduct state.result).IsPrimitive)
    (output : SparsePolyZZ × Array SparsePolyZZ)
    (hrun : Generated.StrictRecombine.vanHoeijLoop ops termination lifted
      modulus initial maximum hinitial state hstate = .ok output) :
    Associated
      (SparsePolyZZ.toPoly state.fStar * factorArrayProduct state.result)
      (vanHoeijOutputProduct output) ∧
    (∀ houtputNonempty : 0 < output.1.size,
      StrictPolynomialMod.SparsePolyZZCanonical output.1 ∧
      (SparsePolyZZ.toPoly output.1 *
        factorArrayProduct output.2).IsPrimitive) := by
  rw [Generated.StrictRecombine.vanHoeijLoop] at hrun
  split at hrun
  next hdone =>
    have hout := Except.ok.inj hrun
    subst output
    have hfStarNe : state.fStar ≠ #[] := by
      intro hempty
      have hsize : state.fStar.size = 0 := by
        simpa using congrArg Array.size hempty
      omega
    constructor
    · simp [vanHoeijOutputProduct, hfStarNe]
    · intro houtputNonempty
      simpa using And.intro hcanonical hprimitive
  next hdone =>
    split at hrun
    next fault hgather => contradiction
    next activeLifted hgather =>
      let hdimension : state.matrix.size =
          activeLifted.size + state.currentColumns := by
        rw [ops.gather_size state.active lifted activeLifted hgather]
        exact hstate.dimension
      split at hrun
      next fault hprepare => contradiction
      next prepared hprepare =>
        let matrix' := prepared.1.1
        let currentColumns' := prepared.1.2.1
        let candidates := prepared.1.2.2
        dsimp only at hrun
        split at hrun
        next fault hvalidate => contradiction
        next fStar' result' consumed hvalidate =>
          by_cases hfound : ∃ index, ∃ hindex : index < consumed.size,
              consumed[index] = true
          · rw [dif_pos hfound] at hrun
            split at hrun
            next fault hremove => contradiction
            next activeNext hremove =>
              split at hrun
              next fault hreset => contradiction
              next resetMatrix resetBound hreset =>
                have hcanonical' := validateCandidates_fStar_canonical
                  ops.validation state.fStar activeLifted modulus candidates
                  state.result fStar' result' consumed hcanonical hvalidate
                have hprimitive' := validateCandidates_preserves_primitive
                  ops.validation state.fStar activeLifted modulus candidates
                  state.result fStar' result' consumed hprimitive hvalidate
                have hnonempty' : 0 < fStar'.size := by
                  by_contra hnot
                  have hempty : fStar' = #[] := by
                    have hsize : fStar'.size = 0 := Nat.eq_zero_of_not_pos hnot
                    apply Array.ext (by simpa using hsize)
                    intro index hindex hindexEmpty
                    simp at hindexEmpty
                  subst fStar'
                  have hcontent :=
                    Polynomial.isPrimitive_iff_content_eq_one.mp hprimitive'
                  simp [SparsePolyZZ.toPoly] at hcontent
                rcases validateCandidates_product_unit_scalar ops.validation
                    state.fStar activeLifted modulus candidates state.result
                    fStar' result' consumed hprimitive hvalidate with
                  ⟨scalar, hscalar, hlive⟩
                let nextState : Generated.StrictRecombine.VanHoeijState :=
                  { active := activeNext.1, fStar := fStar', result := result',
                    matrix := resetMatrix, currentColumns := 0,
                    shortBound := resetBound, target := 0 }
                have hnextValid :
                    Generated.StrictRecombine.VanHoeijStateValid ops nextState :=
                  ⟨(ops.reset_valid activeNext.1.size resetMatrix resetBound
                      hreset).1,
                    by simpa [nextState] using
                      (ops.reset_valid activeNext.1.size resetMatrix resetBound
                        hreset).2⟩
                have htail := vanHoeijLoop_product_associated ops termination
                  lifted modulus initial maximum hinitial nextState hnextValid
                  (by simpa [nextState] using hnonempty')
                  (by simpa [nextState] using hcanonical')
                  (by simpa [nextState] using hprimitive') output hrun
                have hconstantUnit :
                    IsUnit (Polynomial.C scalar : Polynomial Int) :=
                  Polynomial.isUnit_C.mpr hscalar
                constructor
                · exact (Associated.of_eq hlive).trans
                    ((associated_isUnit_mul_left_iff hconstantUnit).mpr htail.1)
                · exact htail.2
          · rw [dif_neg hfound] at hrun
            split at hrun
            next target' hprecision =>
              let nextState : Generated.StrictRecombine.VanHoeijState :=
                { active := state.active, fStar := state.fStar,
                  result := state.result, matrix := matrix',
                  currentColumns := currentColumns',
                  shortBound := state.shortBound, target := target' }
              have hnextValid :
                  Generated.StrictRecombine.VanHoeijStateValid ops nextState :=
                ⟨prepared.2.1, by
                  have hactiveSize :=
                    ops.gather_size state.active lifted activeLifted hgather
                  dsimp [nextState]
                  exact prepared.2.2.trans (by omega)⟩
              exact vanHoeijLoop_product_associated ops termination lifted modulus
                initial maximum hinitial nextState hnextValid
                (by simpa [nextState] using hnonempty)
                (by simpa [nextState] using hcanonical)
                (by simpa [nextState] using hprimitive) output hrun
            next hprecision =>
              split at hrun
              next fault hfallback => contradiction
              next fallback hfallback =>
                split at hrun
                next fault happend => contradiction
                next appended happend =>
                  have hout := Except.ok.inj hrun
                  subst output
                  have hfStarPrimitive :
                      (SparsePolyZZ.toPoly state.fStar).IsPrimitive :=
                    isPrimitive_left_of_mul hprimitive
                  have hzassenhaus := zassenhausRecombine_product_associated
                    ops.zassenhausTermination state.fStar activeLifted fallback
                    modulus hcanonical hfStarPrimitive hfallback
                  have happendProduct := appendFallback_product fallback
                    state.result appended happend
                  unfold vanHoeijOutputProduct
                  simp only [Array.isEmpty_empty, Bool.true_eq, if_true]
                  rw [happendProduct]
                  constructor
                  · exact (hzassenhaus.mul_right
                      (factorArrayProduct state.result)).trans
                      (Associated.of_eq (mul_comm _ _))
                  · intro houtputNonempty
                    simp at houtputNonempty
termination_by
  (state.active.size,
    Generated.StrictRecombine.precisionRank state.target initial maximum)
decreasing_by
  · decreasing_tactic
  · exact Prod.Lex.right _
      (Generated.StrictRecombine.nextPrecision_retry_decreases state.target
        initial maximum _ hinitial (by assumption))

/-- Product refinement of the complete generated C++
`__vanhoeij_recombine` entry.  Initialization, the lexicographically
well-founded main loop, the concrete Zassenhaus fallback, and the final
positive-degree append/sort block are all the generated execution path. -/
theorem __vanhoeij_recombine_raw_ir_product_associated
    (ops : Generated.StrictRecombine.VanHoeijRawOps)
    (termination : Generated.StrictRecombine.VanHoeijTermination ops)
    (f : SparsePolyZZ) (lifted : Array SparsePolyZZ) (modulus : ZZ)
    (hdegree : 2 ≤ (get_deg f).toNatClampNeg)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hprimitive : (SparsePolyZZ.toPoly f).IsPrimitive)
    (output : Array SparsePolyZZ)
    (hrun : Generated.StrictRecombine.__vanhoeij_recombine_raw_ir ops
      termination f lifted modulus hdegree = .ok output) :
    Associated (SparsePolyZZ.toPoly f) (factorArrayProduct output) := by
  rw [Generated.StrictRecombine.__vanhoeij_recombine_raw_ir] at hrun
  dsimp only at hrun
  split at hrun
  next fault hreset => contradiction
  next matrix bound hreset =>
    split at hrun
    next fault hloop => contradiction
    next fStar result hloop =>
      have hout := Except.ok.inj hrun
      subst output
      have hfNonempty : 0 < f.size := by
        by_contra hnot
        have hsize : f.size = 0 := Nat.eq_zero_of_not_pos hnot
        have hempty : f = #[] := by
          apply Array.ext (by simpa using hsize)
          intro index hindex hindexEmpty
          simp at hindexEmpty
        subst f
        simp [get_deg] at hdegree
      have hloopInvariant := vanHoeijLoop_product_associated ops termination
        lifted modulus
        (min (if 3 * lifted.size > (get_deg f).toNatClampNeg + 1 then 30 else 10)
          (((get_deg f).toNatClampNeg + 1) / 2))
        (((get_deg f).toNatClampNeg + 1) / 2)
        (by
          apply lt_min
          · split <;> omega
          · apply Nat.div_pos <;> omega)
        { active := (Array.range lifted.size).map
            (fun index => index.toUInt32.toInt32),
          fStar := f, result := #[], matrix := matrix,
          currentColumns := 0, shortBound := bound, target := 0 }
        ⟨(ops.reset_valid lifted.size matrix bound hreset).1, by
          rw [(ops.reset_valid lifted.size matrix bound hreset).2]
          simp⟩ hfNonempty hcanonical (by simpa [factorArrayProduct] using hprimitive)
        (fStar, result) hloop
      by_cases hfStarNonempty : 0 < fStar.size
      · rcases hloopInvariant.2 hfStarNonempty with
          ⟨hfStarCanonical, hlivePrimitive⟩
        have hfinish := finishZassenhaus_product_associated fStar result
          hfStarCanonical hlivePrimitive
        have hfStarNe : fStar ≠ #[] := by
          intro hempty
          subst fStar
          simp at hfStarNonempty
        simpa [factorArrayProduct] using hloopInvariant.1.trans (by
          simpa [vanHoeijOutputProduct, hfStarNe] using hfinish)
      · have hfStarEmpty : fStar = #[] := by
          have hsize : fStar.size = 0 := Nat.eq_zero_of_not_pos hfStarNonempty
          apply Array.ext (by simpa using hsize)
          intro index hindex hindexEmpty
          simp at hindexEmpty
        subst fStar
        rw [factorArrayProduct_finishZassenhaus]
        simp only [Array.size_empty, lt_self_iff_false, dite_false]
        simpa [factorArrayProduct, vanHoeijOutputProduct] using hloopInvariant.1

end Refinement.StrictRecombine
