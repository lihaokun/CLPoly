/-
  Strict Hensel refinement boundary.

  This module deliberately exports no candidate from the legacy generated
  corpus.  The former safe wrapper selected the L2 `hensel_lift`
  implementation whenever the candidate failed a post-hoc `HenselCorrect`
  check.  That was not an L1→L2 refinement and has been removed.
-/
import CLPoly.Algorithm.Hensel
import CLPoly.Generated.StrictHensel
import CLPoly.Math.Bigint
import CLPoly.Refinement.Basic

set_option autoImplicit false

namespace Refinement

namespace StrictHensel

noncomputable def toPolyMod (m : Nat) (f : SparsePolyZZ) :
    Polynomial (ZMod m) :=
  Polynomial.map (Int.castRingHom (ZMod m)) (SparsePolyZZ.toPoly f)

@[simp] theorem toPolyMod_empty (m : Nat) :
    toPolyMod m (#[] : SparsePolyZZ) = 0 := by
  simp [toPolyMod, SparsePolyZZ.toPoly]

@[simp] theorem toPolyMod_push (m : Nat) (f : SparsePolyZZ)
    (degree : Nat) (coefficient : Int) :
    toPolyMod m (f.push (UMonomial.mk degree, coefficient)) =
      toPolyMod m f + Polynomial.monomial degree (coefficient : ZMod m) := by
  simp [toPolyMod, SparsePolyZZ.toPoly]

/-- Exact semantic effect of the source helper which conditionally appends a
nonzero coefficient. -/
theorem pushNonzero_toPolyMod (m : Nat) (f : SparsePolyZZ)
    (degree : Nat) (coefficient : Int) :
    toPolyMod m (Generated.StrictHensel.pushNonzero f degree coefficient) =
      toPolyMod m f + Polynomial.monomial degree (coefficient : ZMod m) := by
  by_cases hzero : coefficient = 0
  · simp [Generated.StrictHensel.pushNonzero, hzero]
  · simp [Generated.StrictHensel.pushNonzero, hzero, toPolyMod_push]

/-- Polynomial represented by an arbitrary suffix of the concrete sparse
array.  Keeping this definition on lists lets the merge-loop invariant follow
the two source iterator indices exactly. -/
noncomputable def termsToPolyMod (m : Nat) (xs : List (UMonomial × Int)) :
    Polynomial (ZMod m) :=
  (xs.map fun term =>
    Polynomial.monomial term.1.deg (term.2 : ZMod m)).sum

theorem toPolyMod_eq_termsToPolyMod (m : Nat) (f : SparsePolyZZ) :
    toPolyMod m f = termsToPolyMod m f.toList := by
  unfold toPolyMod SparsePolyZZ.toPoly termsToPolyMod
  induction f.toList with
  | nil => simp
  | cons term terms ih => simp [ih]

@[simp] theorem termsToPolyMod_nil (m : Nat) :
    termsToPolyMod m [] = 0 := by
  simp [termsToPolyMod]

@[simp] theorem termsToPolyMod_cons (m : Nat) (term : UMonomial × Int)
    (terms : List (UMonomial × Int)) :
    termsToPolyMod m (term :: terms) =
      Polynomial.monomial term.1.deg (term.2 : ZMod m) +
        termsToPolyMod m terms := by
  simp [termsToPolyMod]

private theorem drop_eq_getElem_cons {α : Type} [Inhabited α] (input : Array α)
    (i : Nat) (hi : i < input.size) :
    input.toList.drop i = input[i]! :: input.toList.drop (i + 1) := by
  have hlen : i < input.toList.length := by simpa using hi
  have hdrop := List.drop_eq_getElem_cons (l := input.toList) (i := i) hlen
  have hget : input.toList[i] = input[i]! := by
    calc
      input.toList[i] = input[i] := by simp
      _ = input[i]! := by
        rw [getElem!_def, getElem?_def]
        simp [hi]
  rw [hget] at hdrop
  exact hdrop

/-- Shifting a concrete divisor term by the current quotient degree and
scaling it by the quotient coefficient is exactly multiplication by the
corresponding monomial in `ZMod m`. -/
theorem shiftedTerm_toPolyMod (m degreeShift : Nat) (coefficient : Int)
    (term : UMonomial × Int) :
    Polynomial.monomial (term.1.deg + degreeShift)
        ((coefficient * term.2 : Int) : ZMod m) =
      Polynomial.monomial degreeShift (coefficient : ZMod m) *
        Polynomial.monomial term.1.deg (term.2 : ZMod m) := by
  rw [Polynomial.monomial_mul_monomial]
  simp only [Int.cast_mul]
  rw [add_comm]

theorem shiftedMonomial (m degreeShift degree : Nat) (a b : ZMod m) :
    Polynomial.monomial (degree + degreeShift) (a * b) =
      Polynomial.monomial degreeShift a * Polynomial.monomial degree b := by
  rw [Polynomial.monomial_mul_monomial, add_comm]


private theorem intCast_fmod_natCast (m : Nat) (q : Int) :
    ((Int.fmod q (m : Int) : Int) : ZMod m) = (q : ZMod m) := by
  rw [Int.fmod_eq_emod_of_nonneg q (by omega)]
  rw [ZMod.intCast_eq_intCast_iff]
  refine (Int.modEq_iff_dvd).mpr ?_
  use q / (m : Int)
  have h := Int.mul_ediv_add_emod q (m : Int)
  omega

private theorem intCast_negative_product_remainder (m : Nat) (a b : Int) :
    ((Int.fmod ((m : Int) - (a * b % (m : Int)))
        (m : Int) : Int) : ZMod m) = -(a : ZMod m) * (b : ZMod m) := by
  rw [intCast_fmod_natCast]
  simp

private theorem intCast_difference_remainder (m : Nat) (a b c : Int) :
    ((Int.fmod (a - b * c) (m : Int) : Int) : ZMod m) =
      (a : ZMod m) - (b : ZMod m) * (c : ZMod m) := by
  rw [intCast_fmod_natCast]
  simp

set_option maxHeartbeats 0 in
/-- Exact semantic invariant of the generated ordered merge.  It follows the
well-founded recursion emitted for the C++ iterator loop and accounts for the
two unprocessed array suffixes; no polynomial-division specification is used
to define or replace the result. -/
theorem divmodMergeLoop_toPolyMod (m degreeShift : Nat)
    (r g : SparsePolyZZ) (coefficient : Int) :
    ∀ (rIndex gIndex : Nat) (result : SparsePolyZZ),
      toPolyMod m
          (Generated.StrictHensel.divmodMergeLoop r g coefficient (m : Int)
            degreeShift rIndex gIndex result) =
        toPolyMod m result + termsToPolyMod m (r.toList.drop rIndex) -
          Polynomial.monomial degreeShift (coefficient : ZMod m) *
            termsToPolyMod m (g.toList.drop gIndex) := by
  intro rIndex gIndex result
  refine Generated.StrictHensel.divmodMergeLoop.induct r g coefficient
    (m : Int) degreeShift
    (motive := fun rIndex gIndex result =>
      toPolyMod m
          (Generated.StrictHensel.divmodMergeLoop r g coefficient (m : Int)
            degreeShift rIndex gIndex result) =
        toPolyMod m result + termsToPolyMod m (r.toList.drop rIndex) -
          Polynomial.monomial degreeShift (coefficient : ZMod m) *
            termsToPolyMod m (g.toList.drop gIndex)) ?_ ?_ ?_ ?_ ?_ ?_
    rIndex gIndex result
  · intro ri gi acc hmore hgDone
    dsimp only
    intro ih
    have hri : ri < r.size := by omega
    have hrDrop := drop_eq_getElem_cons r ri hri
    have hgDrop : g.toList.drop gi = [] := by
      apply List.drop_eq_nil_of_le
      simpa using hgDone
    rw [Generated.StrictHensel.divmodMergeLoop.eq_1]
    simp [hmore, hgDone]
    simp only [ZZ.fdiv_r] at ih
    simp only [ZZ.fdiv_r]
    rw [ih, pushNonzero_toPolyMod, hrDrop, hgDrop,
      termsToPolyMod_cons, termsToPolyMod_nil,
      intCast_fmod_natCast]
    ring
  · intro ri gi acc hmore hgMore hrDone
    dsimp only
    intro ih
    have hgi : gi < g.size := by omega
    have hgDrop := drop_eq_getElem_cons g gi hgi
    have hrDrop : r.toList.drop ri = [] := by
      apply List.drop_eq_nil_of_le
      simpa using hrDone
    rw [Generated.StrictHensel.divmodMergeLoop.eq_1]
    simp [hmore, hgMore, hrDone]
    simp only [Int.emod_emod] at ih
    simp only [ZZ.fdiv_r] at ih
    simp only [ZZ.fdiv_r]
    rw [ih, pushNonzero_toPolyMod, hrDrop, hgDrop,
      termsToPolyMod_cons, termsToPolyMod_nil,
      intCast_negative_product_remainder]
    rw [show -(coefficient : ZMod m) * (g[gi]!.2 : ZMod m) =
      -((coefficient : ZMod m) * (g[gi]!.2 : ZMod m)) by ring,
      map_neg, shiftedMonomial]
    ring
  · intro ri gi acc hmore hgMore hrMore
    dsimp only
    intro hdegree
    intro ih
    have hri : ri < r.size := by omega
    have hrDrop := drop_eq_getElem_cons r ri hri
    rw [Generated.StrictHensel.divmodMergeLoop.eq_1]
    simp [hmore, hgMore, hrMore, hdegree]
    simp only [ZZ.fdiv_r] at ih
    simp only [ZZ.fdiv_r]
    rw [ih, pushNonzero_toPolyMod, hrDrop, termsToPolyMod_cons,
      intCast_fmod_natCast]
    ring
  · intro ri gi acc hmore hgMore hrMore
    dsimp only
    intro hnotGreater hless
    intro ih
    have hgi : gi < g.size := by omega
    have hgDrop := drop_eq_getElem_cons g gi hgi
    rw [Generated.StrictHensel.divmodMergeLoop.eq_1]
    simp [hmore, hgMore, hrMore, hnotGreater, hless]
    simp only [Int.emod_emod] at ih
    simp only [ZZ.fdiv_r] at ih
    simp only [ZZ.fdiv_r]
    rw [ih, pushNonzero_toPolyMod, hgDrop, termsToPolyMod_cons,
      intCast_negative_product_remainder]
    rw [show -(coefficient : ZMod m) * (g[gi]!.2 : ZMod m) =
      -((coefficient : ZMod m) * (g[gi]!.2 : ZMod m)) by ring,
      map_neg, shiftedMonomial]
    ring
  · intro ri gi acc hmore hgMore hrMore
    dsimp only
    intro hnotGreater hnotLess
    intro ih
    have hri : ri < r.size := by omega
    have hgi : gi < g.size := by omega
    have hrDrop := drop_eq_getElem_cons r ri hri
    have hgDrop := drop_eq_getElem_cons g gi hgi
    have hdegrees : r[ri]!.1.deg = g[gi]!.1.deg + degreeShift := by omega
    rw [Generated.StrictHensel.divmodMergeLoop.eq_1]
    simp [hmore, hgMore, hrMore, hnotGreater, hnotLess]
    simp only [ZZ.fdiv_r] at ih
    simp only [ZZ.fdiv_r]
    rw [ih, pushNonzero_toPolyMod, hrDrop, hgDrop,
      termsToPolyMod_cons, intCast_difference_remainder]
    rw [map_sub, hdegrees, shiftedMonomial]
    simp only [termsToPolyMod_cons]
    ring
  · intro ri gi acc hdone
    have hrDrop : r.toList.drop ri = [] := by
      apply List.drop_eq_nil_of_le
      simp only [not_or] at hdone
      simpa using hdone.1
    have hgDrop : g.toList.drop gi = [] := by
      apply List.drop_eq_nil_of_le
      simp only [not_or] at hdone
      simpa using hdone.2
    rw [Generated.StrictHensel.divmodMergeLoop.eq_1]
    simp [hdone, hrDrop, hgDrop]

/-- Direct semantic equation for the concrete source subtraction helper.  The
source deliberately starts both iterators at one because their leading terms
are cancelled by the outer long-division step. -/
theorem divmodRemainder_toPolyMod (m degreeShift : Nat)
    (r g : SparsePolyZZ) (coefficient : Int) :
    toPolyMod m
        (Generated.StrictHensel.divmodRemainder r g coefficient (m : Int)
          degreeShift) =
      termsToPolyMod m (r.toList.drop 1) -
        Polynomial.monomial degreeShift (coefficient : ZMod m) *
          termsToPolyMod m (g.toList.drop 1) := by
  rw [Generated.StrictHensel.divmodRemainder,
    divmodMergeLoop_toPolyMod]
  simp

/-- When the concrete quotient coefficient cancels the two leading terms,
the generated helper represents subtraction of that quotient monomial times
the complete divisor. -/
theorem divmodRemainder_eq_sub (m degreeShift : Nat)
    (r g : SparsePolyZZ) (coefficient : Int)
    (hr : 0 < r.size) (hg : 0 < g.size)
    (hdegree : r[0]!.1.deg = g[0]!.1.deg + degreeShift)
    (hcoefficient : (r[0]!.2 : ZMod m) =
      (coefficient : ZMod m) * (g[0]!.2 : ZMod m)) :
    toPolyMod m
        (Generated.StrictHensel.divmodRemainder r g coefficient (m : Int)
          degreeShift) =
      toPolyMod m r -
        Polynomial.monomial degreeShift (coefficient : ZMod m) *
          toPolyMod m g := by
  have hrDrop := drop_eq_getElem_cons r 0 hr
  have hgDrop := drop_eq_getElem_cons g 0 hg
  have hrFull : r.toList = r[0]! :: r.toList.drop 1 := by
    simpa using hrDrop
  have hgFull : g.toList = g[0]! :: g.toList.drop 1 := by
    simpa using hgDrop
  have hrPoly : toPolyMod m r =
      Polynomial.monomial r[0]!.1.deg (r[0]!.2 : ZMod m) +
        termsToPolyMod m (r.toList.drop 1) := by
    rw [toPolyMod_eq_termsToPolyMod]
    calc
      termsToPolyMod m r.toList =
          termsToPolyMod m (r[0]! :: r.toList.drop 1) :=
        congrArg (termsToPolyMod m) hrFull
      _ = _ := termsToPolyMod_cons m r[0]! (r.toList.drop 1)
  have hgPoly : toPolyMod m g =
      Polynomial.monomial g[0]!.1.deg (g[0]!.2 : ZMod m) +
        termsToPolyMod m (g.toList.drop 1) := by
    rw [toPolyMod_eq_termsToPolyMod]
    calc
      termsToPolyMod m g.toList =
          termsToPolyMod m (g[0]! :: g.toList.drop 1) :=
        congrArg (termsToPolyMod m) hgFull
      _ = _ := termsToPolyMod_cons m g[0]! (g.toList.drop 1)
  rw [divmodRemainder_toPolyMod, hrPoly, hgPoly, hdegree,
    hcoefficient, shiftedMonomial]
  ring

/-- The coefficient computed by the concrete `fmod(r₀ * inverse, m)` source
expression cancels the divisor's leading coefficient.  The only premise is
the success bit returned by the actual `ZZ.invert` implementation. -/
theorem divmodCoefficient_mul_leading (m : Nat) (r g : SparsePolyZZ)
    (hinvert : (ZZ.invert 0 g[0]!.2 (m : Int)).1 = true) :
    (Generated.StrictHensel.divmodCoefficient r
        (ZZ.invert 0 g[0]!.2 (m : Int)).2 (m : Int) : ZMod m) *
      (g[0]!.2 : ZMod m) = (r[0]!.2 : ZMod m) := by
  have hinv := CLPoly.Math.ZZ.invert_success_mul_eq_one g[0]!.2 m hinvert
  rw [Generated.StrictHensel.divmodCoefficient]
  simp only [ZZ.fdiv_r]
  rw [intCast_fmod_natCast]
  simp only [Int.cast_mul]
  ring_nf at hinv ⊢
  rw [hinv]
  simp

/-- Exact natural value of the signed C++ degree word for a nonempty sparse
integer polynomial whose leading exponent fits the source's signed range. -/
theorem get_deg_toNatClampNeg_eq_head (f : SparsePolyZZ)
    (hnonempty : 0 < f.size) (hdegree : f[0]!.1.deg < 2 ^ 63) :
    (get_deg f).toNatClampNeg = f[0]!.1.deg := by
  have hne : f ≠ #[] := by
    intro hempty
    subst f
    simp at hnonempty
  have hwordNat : f[0]!.1.deg.toUInt64.toNat = f[0]!.1.deg := by
    change (OfNat.ofNat f[0]!.1.deg : UInt64).toNat = f[0]!.1.deg
    rw [UInt64.toNat_ofNat, Nat.mod_eq_of_lt]
    omega
  have hsignedNat : f[0]!.1.deg.toUInt64.toInt64.toNatClampNeg =
      f[0]!.1.deg := by
    rw [UInt64_toInt64_toNatClampNeg_eq_toNat_of_lt (by
      simpa [hwordNat] using hdegree), hwordNat]
  have hget : get_deg f = f[0]!.1.deg.toUInt64.toInt64 := by
    simp [get_deg, Array.isEmpty_iff, hne, getElem!_pos f 0 hnonempty]
  rw [hget, hsignedNat]

private theorem int64_sub_toNatClampNeg (a b : Int64)
    (hb : (0 : Int64) ≤ b) (hba : b ≤ a) :
    (a - b).toNatClampNeg = a.toNatClampNeg - b.toNatClampNeg := by
  change (a - b).toInt.toNat = a.toInt.toNat - b.toInt.toNat
  rw [Int64.toInt_sub, Int.bmod_eq_of_le]
  · have hbInt : 0 ≤ b.toInt := by
      simpa [Int64.le_iff_toInt_le] using hb
    have hbaInt : b.toInt ≤ a.toInt := by
      simpa [Int64.le_iff_toInt_le] using hba
    omega
  · have hbInt : 0 ≤ b.toInt := by
      simpa [Int64.le_iff_toInt_le] using hb
    have hbaInt : b.toInt ≤ a.toInt := by
      simpa [Int64.le_iff_toInt_le] using hba
    have haUpper := Int64.toInt_lt a
    omega
  · have hbInt : 0 ≤ b.toInt := by
      simpa [Int64.le_iff_toInt_le] using hb
    have hbaInt : b.toInt ≤ a.toInt := by
      simpa [Int64.le_iff_toInt_le] using hba
    have haUpper := Int64.toInt_lt a
    omega

/-- The source's signed degree subtraction and clamp compute the exact
natural exponent shift used to cancel the current leading term. -/
theorem get_deg_sub_toNatClampNeg_eq_shift (r g : SparsePolyZZ)
    (hr : 0 < r.size) (hg : 0 < g.size)
    (hrDegree : r[0]!.1.deg < 2 ^ 63)
    (hgDegree : g[0]!.1.deg < 2 ^ 63)
    (hactive : get_deg r ≥ get_deg g) :
    r[0]!.1.deg = g[0]!.1.deg +
      (get_deg r - get_deg g).toNatClampNeg := by
  have hrNat := get_deg_toNatClampNeg_eq_head r hr hrDegree
  have hgNat := get_deg_toNatClampNeg_eq_head g hg hgDegree
  have hgNonneg : (0 : Int64) ≤ get_deg g := by
    by_cases hz : g[0]!.1.deg = 0
    · have hgne : g ≠ #[] := by
        intro hempty
        subst g
        simp at hg
      have hget : get_deg g = g[0]!.1.deg.toUInt64.toInt64 := by
        simp [get_deg, Array.isEmpty_iff, hgne, getElem!_pos g 0 hg]
      rw [hget, hz]
      decide
    · have hpos : 0 < (get_deg g).toNatClampNeg := by omega
      exact Int64.le_of_lt
        ((Int64.toNatClampNeg_pos_iff (get_deg g)).mp hpos)
  rw [int64_sub_toNatClampNeg _ _ hgNonneg hactive, hrNat, hgNat]
  have hnat : g[0]!.1.deg ≤ r[0]!.1.deg := by
    rw [← hgNat, ← hrNat]
    exact Int64.toNatClampNeg_le hactive
  omega

/-- The source's zero-quotient branch erases a leading term whose coefficient
is zero modulo `m`; therefore the decoded polynomial is unchanged. -/
theorem eraseLeading_of_divmodCoefficient_zero (m : Nat)
    (r g : SparsePolyZZ) (hr : 0 < r.size)
    (hinvert : (ZZ.invert 0 g[0]!.2 (m : Int)).1 = true)
    (hzero : Generated.StrictHensel.divmodCoefficient r
      (ZZ.invert 0 g[0]!.2 (m : Int)).2 (m : Int) = 0) :
    toPolyMod m (r.eraseIdxIfInBounds 0) = toPolyMod m r := by
  have hcancel := divmodCoefficient_mul_leading m r g hinvert
  rw [hzero] at hcancel
  simp only [Int.cast_zero, zero_mul] at hcancel
  have hrDrop := drop_eq_getElem_cons r 0 hr
  have hrFull : r.toList = r[0]! :: r.toList.drop 1 := by
    simpa using hrDrop
  rw [toPolyMod_eq_termsToPolyMod, toPolyMod_eq_termsToPolyMod]
  simp only [Array.toList_eraseIdxIfInBounds, List.eraseIdx_zero]
  rw [hrFull]
  simp [hcancel]

/-- Representation facts needed at every active state of the exact generated
division trace.  This predicate contains no result polynomial: it certifies
only that the source's signed degree arithmetic is within range. -/
def DivmodTraceDegreeSafe (g : SparsePolyZZ) (inverse : Int) (m : Nat) :
    {r q : SparsePolyZZ} →
      Generated.StrictHensel.DivmodTrace g inverse (m : Int) r q → Prop
  | _, _, .done _ _ _ => True
  | _, _, .vanished r _ _ _ next =>
      0 < r.size ∧ r[0]!.1.deg < 2 ^ 63 ∧ get_deg r ≥ get_deg g ∧
        DivmodTraceDegreeSafe g inverse m next
  | _, _, .subtract r _ _ _ next =>
      0 < r.size ∧ r[0]!.1.deg < 2 ^ 63 ∧ get_deg r ≥ get_deg g ∧
        DivmodTraceDegreeSafe g inverse m next

/-- Induction over the exact finite C++ long-division trace preserves the
decoded equation `remainder + quotient * divisor`.  Both recursive branches
use their concrete array updates. -/
theorem divmodLoop_preserves (m : Nat) (g : SparsePolyZZ)
    (hg : 0 < g.size) (hgDegree : g[0]!.1.deg < 2 ^ 63)
    (hinvert : (ZZ.invert 0 g[0]!.2 (m : Int)).1 = true) :
    ∀ {r q : SparsePolyZZ}
      (trace : Generated.StrictHensel.DivmodTrace g
        (ZZ.invert 0 g[0]!.2 (m : Int)).2 (m : Int) r q),
      DivmodTraceDegreeSafe g
          (ZZ.invert 0 g[0]!.2 (m : Int)).2 m trace →
      let output := Generated.StrictHensel.divmodLoop g
        (ZZ.invert 0 g[0]!.2 (m : Int)).2 (m : Int) trace
      toPolyMod m output.2 + toPolyMod m output.1 * toPolyMod m g =
        toPolyMod m r + toPolyMod m q * toPolyMod m g := by
  intro r q trace hsafe
  induction trace with
  | done r q inactive =>
      rfl
  | vanished r q active zero next ih =>
      simp only [DivmodTraceDegreeSafe] at hsafe
      rw [Generated.StrictHensel.divmodLoop]
      dsimp only
      rw [ih hsafe.2.2.2]
      rw [eraseLeading_of_divmodCoefficient_zero m r g hsafe.1 hinvert zero]
  | subtract r q active nonzero next ih =>
      simp only [DivmodTraceDegreeSafe] at hsafe
      rw [Generated.StrictHensel.divmodLoop]
      dsimp only
      rw [ih hsafe.2.2.2]
      have hdegree := get_deg_sub_toNatClampNeg_eq_shift r g hsafe.1 hg
        hsafe.2.1 hgDegree hsafe.2.2.1
      have hcoefficient := divmodCoefficient_mul_leading m r g hinvert
      rw [divmodRemainder_eq_sub m
        (get_deg r - get_deg g).toNatClampNeg r g
        (Generated.StrictHensel.divmodCoefficient r
          (ZZ.invert 0 g[0]!.2 (m : Int)).2 (m : Int))
        hsafe.1 hg hdegree hcoefficient.symm]
      rw [toPolyMod_push]
      ring

/-- The concrete generated `__upoly_mod_coeff` call always succeeds and
preserves the decoded polynomial modulo the requested modulus. -/
theorem __upoly_mod_coeff_raw_ir_refines (f : SparsePolyZZ) (m : Nat) :
    ∃ output,
      Generated.StrictHensel.__upoly_mod_coeff_raw_ir f (m : Int) =
        .ok output ∧
      toPolyMod m output = toPolyMod m f := by
  let output := f.filterMap fun term =>
    let coefficient := ZZ.fdiv_r term.snd term.snd (m : Int)
    if coefficient != 0 then some (term.fst, coefficient) else none
  refine ⟨output, rfl, ?_⟩
  unfold toPolyMod output
  simp only [SparsePolyZZ.toPoly, Array.toList_filterMap]
  induction f.toList with
  | nil => simp
  | cons term terms ih =>
      have ihtail := ih
      simp only [ZZ.fdiv_r, bne_iff_ne, ne_eq] at ihtail
      rw [show (fun x : UMonomial × Int =>
          if ¬x.2.fmod (m : Int) = 0 then
            some (x.1, x.2.fmod (m : Int)) else none) =
          (fun x => if x.2.fmod (m : Int) = 0 then none
            else some (x.1, x.2.fmod (m : Int))) by
            funext x
            simp [ite_not]] at ihtail
      have hcast : ((Int.fmod term.2 (m : Int) : Int) : ZMod m) =
          (term.2 : ZMod m) := intCast_fmod_natCast m term.2
      by_cases hz : Int.fmod term.2 (m : Int) = 0
      · have htermzero : (term.2 : ZMod m) = 0 := by
          rw [← hcast, hz]
          simp
        simp [ZZ.fdiv_r, hz, htermzero, ihtail]
      · simp [ZZ.fdiv_r, hz, hcast, ihtail]

/-- A valid divisor and successful concrete GMP inverse make the strict
generated modular long-division entry execute to the unique result obtained
from its exact finite source trace. -/
theorem __upoly_divmod_mod_raw_ir_terminates
    (termination : Generated.StrictHensel.DivmodTermination)
    (f g : SparsePolyZZ) (m : ZZ)
    (hg : g.isEmpty = false)
    (hinvert : (ZZ.invert 0 g[0]!.2 m).1 = true) :
    ∃ q r,
      Generated.StrictHensel.__upoly_divmod_mod_raw_ir termination f g m =
        .ok (q, r) := by
  let reduced := Generated.StrictHensel.modCoeffOutput f m
  let output := Generated.StrictHensel.divmodLoop g
    (ZZ.invert 0 g[0]!.2 m).2 m
      (termination.trace f g reduced m (by rfl))
  refine ⟨output.1, output.2, ?_⟩
  simp [Generated.StrictHensel.__upoly_divmod_mod_raw_ir, hg, hinvert,
    Generated.StrictHensel.__upoly_mod_coeff_raw_ir, reduced, output]

/-- The generated divide/reduce/compact loop represents the exact coefficient
quotient modulo `m`; removing coefficients whose residues are zero does not
change the decoded `ZMod m` polynomial. -/
theorem divideThenReduceCoeffs_toPolyMod (f : SparsePolyZZ) (m : Nat) :
    toPolyMod m (Generated.StrictHensel.divideThenReduceCoeffs f (m : Int)) =
      toPolyMod m (Generated.StrictHensel.divideCoeffs f (m : Int)) := by
  unfold toPolyMod
  simp only [Generated.StrictHensel.divideThenReduceCoeffs,
    Generated.StrictHensel.divideCoeffs, SparsePolyZZ.toPoly,
    Array.toList_filterMap, Array.toList_map, List.map_map]
  induction f.toList with
  | nil => simp
  | cons term terms ih =>
      have ihtail := ih
      simp only [ZZ.fdiv_q, ZZ.fdiv_r, bne_iff_ne, ne_eq] at ihtail
      rw [show (fun x : UMonomial × Int =>
          if ¬(x.2.fdiv (m : Int)).fmod (m : Int) = 0 then
            some (x.1, (x.2.fdiv (m : Int)).fmod (m : Int)) else none) =
          (fun x => if (x.2.fdiv (m : Int)).fmod (m : Int) = 0 then
            none else some (x.1, (x.2.fdiv (m : Int)).fmod (m : Int))) by
            funext x
            simp [ite_not]] at ihtail
      let q := Int.fdiv term.2 (m : Int)
      have hcast : ((Int.fmod q (m : Int) : Int) : ZMod m) =
          (q : ZMod m) := intCast_fmod_natCast m q
      by_cases hz : Int.fmod q (m : Int) = 0
      · have hqzero : (q : ZMod m) = 0 := by
          rw [← hcast, hz]
          simp
        simp [ZZ.fdiv_q, ZZ.fdiv_r, q, hz, hqzero, ihtail]
      · simp [ZZ.fdiv_q, ZZ.fdiv_r, q, hz, hcast, ihtail]

/-- Under the source loop's exact-divisibility precondition, multiplying the
generated coefficient quotient back by `m` reconstructs the precise raw
difference polynomial coefficient-for-coefficient. -/
theorem scaleCoeffs_divideCoeffs (f : SparsePolyZZ) (m : Nat)
    (hdivisible : ∀ term ∈ f.toList, (m : Int) ∣ term.2) :
    Generated.StrictHensel.scaleCoeffs
        (Generated.StrictHensel.divideCoeffs f (m : Int)) (m : Int) = f := by
  apply Array.toList_inj.mp
  simp only [Generated.StrictHensel.scaleCoeffs,
    Generated.StrictHensel.divideCoeffs, Array.toList_map, List.map_map]
  have listLemma : ∀ xs : List (UMonomial × Int),
      (∀ term ∈ xs, (m : Int) ∣ term.2) →
      List.map ((fun term => (term.1, term.2 * (m : Int))) ∘
        fun term => (term.1, ZZ.fdiv_q term.2 term.2 (m : Int))) xs = xs := by
    intro xs hxs
    induction xs with
    | nil => rfl
    | cons term terms ih =>
        have hhead := hxs term (by simp)
        have htail : ∀ item ∈ terms, (m : Int) ∣ item.2 := by
          intro item hitem
          exact hxs item (by simp [hitem])
        simp only [List.map_cons]
        congr 1
        · apply Prod.ext
          · rfl
          · simp only [Function.comp_apply, ZZ.fdiv_q]
            rw [Int.fdiv_eq_ediv_of_nonneg _ (by omega), mul_comm,
              Int.mul_ediv_cancel' hhead]
        · exact ih htail
  exact listLemma f.toList hdivisible

/-- Decoding the generated coefficient-scaling loop is multiplication by the
constant polynomial `m`.  This is the first raw-loop bridge used by both
factor and Bézout updates in the C++ quadratic Hensel step. -/
theorem scaleCoeffs_toPoly (f : SparsePolyZZ) (m : ZZ) :
    SparsePolyZZ.toPoly (Generated.StrictHensel.scaleCoeffs f m) =
      Polynomial.C m * SparsePolyZZ.toPoly f := by
  simp only [Generated.StrictHensel.scaleCoeffs, SparsePolyZZ.toPoly,
    Array.toList_map, List.map_map]
  induction f.toList with
  | nil => simp
  | cons term terms ih =>
      simp only [List.map_cons, List.sum_cons] at ih ⊢
      rw [ih, mul_add]
      congr 1
      change Polynomial.monomial term.1.deg (term.2 * m) =
        Polynomial.C m * Polynomial.monomial term.1.deg term.2
      rw [mul_comm term.2 m]
      exact (Polynomial.C_mul_monomial (R := ℤ) (a := m)
        (b := term.2) (n := term.1.deg)).symm

/-- Complete semantic certificate for the source error-coefficient loop: its
unreduced quotient reconstructs the input exactly, and the concrete compacted
output is that quotient modulo `m`. -/
theorem divideThenReduceCoeffs_certificate (f : SparsePolyZZ) (m : Nat)
    (hdivisible : ∀ term ∈ f.toList, (m : Int) ∣ term.2) :
    toPolyMod m (Generated.StrictHensel.divideThenReduceCoeffs f (m : Int)) =
        toPolyMod m (Generated.StrictHensel.divideCoeffs f (m : Int)) ∧
      Polynomial.C (m : Int) * SparsePolyZZ.toPoly
          (Generated.StrictHensel.divideCoeffs f (m : Int)) =
        SparsePolyZZ.toPoly f := by
  constructor
  · exact divideThenReduceCoeffs_toPolyMod f m
  · rw [← scaleCoeffs_toPoly,
      scaleCoeffs_divideCoeffs f m hdivisible]

end StrictHensel

-- No Hensel L1→L2 theorem or legacy candidate is exported until a strict
-- cpp2lean-generated entry and its direct execution proof are available.

end Refinement
