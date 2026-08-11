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
import CLPoly.Refinement.SquarefreeZp

set_option autoImplicit false

namespace Refinement

namespace StrictHensel

noncomputable def toPolyMod (m : Nat) (f : SparsePolyZZ) :
    Polynomial (ZMod m) :=
  Polynomial.map (Int.castRingHom (ZMod m)) (SparsePolyZZ.toPoly f)

/-- L2 invariant represented by a concrete C++ Hensel tree node at modulus
`m`: its two factors multiply to the target and its stored coefficients form
a Bézout certificate. Both clauses decode arrays to mathematical
polynomials; they do not call model implementations of C++ operations. -/
def HenselNodeInvariant (f : SparsePolyZZ) (m : Nat)
    (node : HenselNode) : Prop :=
  toPolyMod m node.g * toPolyMod m node.h = toPolyMod m f ∧
  toPolyMod m node.s * toPolyMod m node.g +
      toPolyMod m node.t * toPolyMod m node.h = 1

/-- Precise L2 postcondition of the original C++ `__hensel_step` entry.
Besides establishing the factor and Bézout invariants at `m²`, all four
updated polynomials reduce to their concrete input values modulo `m`. -/
def HenselStepCorrect (f : SparsePolyZZ) (m : Nat)
    (input output : HenselNode) : Prop :=
  HenselNodeInvariant f (m ^ 2) output ∧
  toPolyMod m output.g = toPolyMod m input.g ∧
  toPolyMod m output.h = toPolyMod m input.h ∧
  toPolyMod m output.s = toPolyMod m input.s ∧
  toPolyMod m output.t = toPolyMod m input.t

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

/-- Decoding commutes with the canonical projection from a stronger modulus
to a divisor modulus.  This is the raw representation bridge needed to show
that every quadratic update still reduces to its input modulo `m`. -/
theorem map_toPolyMod_of_dvd {smaller larger : Nat}
    (hdiv : smaller ∣ larger) (f : SparsePolyZZ) :
    Polynomial.map (ZMod.castHom hdiv (ZMod smaller))
        (toPolyMod larger f) =
      toPolyMod smaller f := by
  unfold toPolyMod
  rw [Polynomial.map_map]
  congr 1
  ext coefficient
  simp [ZMod.castHom_apply, ZMod.cast_intCast hdiv]

/-- In `ZMod (m²)`, multiplying a coefficient that vanishes modulo `m` by
`m` gives zero.  This is the exact algebraic fact used by quadratic Hensel
correction; it does not appeal to the L2 Hensel implementation. -/
theorem zmod_scale_eq_zero_of_cast_eq_zero (m : Nat) (hm : 0 < m)
    (a : ZMod (m ^ 2))
    (ha : ZMod.castHom (dvd_pow_self m (by omega : 2 ≠ 0))
      (ZMod m) a = 0) :
    (m : ZMod (m ^ 2)) * a = 0 := by
  haveI : NeZero (m ^ 2) := ⟨by positivity⟩
  simp only [ZMod.castHom_apply] at ha
  rw [ZMod.cast_eq_val, ZMod.natCast_eq_zero_iff] at ha
  obtain ⟨k, hk⟩ := ha
  have ha' : a = ((m * k : Nat) : ZMod (m ^ 2)) := by
    rw [← hk, ZMod.natCast_zmod_val]
  rw [ha']
  push_cast
  have hmm : (m : ZMod (m ^ 2)) * (m : ZMod (m ^ 2)) = 0 := by
    rw [← Nat.cast_mul, show m * m = m ^ 2 by ring]
    exact ZMod.natCast_self (m ^ 2)
  rw [← mul_assoc, hmm, zero_mul]

/-- Polynomial form of `zmod_scale_eq_zero_of_cast_eq_zero`. -/
theorem C_scale_mul_eq_zero_of_map_eq_zero (m : Nat) (hm : 0 < m)
    (p : Polynomial (ZMod (m ^ 2)))
    (hp : Polynomial.map
      (ZMod.castHom (dvd_pow_self m (by omega : 2 ≠ 0)) (ZMod m)) p = 0) :
    Polynomial.C (m : ZMod (m ^ 2)) * p = 0 := by
  ext degree
  rw [Polynomial.coeff_C_mul, Polynomial.coeff_zero]
  apply zmod_scale_eq_zero_of_cast_eq_zero m hm
  have hcoeff := congrArg (fun q => q.coeff degree) hp
  simpa [Polynomial.coeff_map] using hcoeff

/-- Algebraic core of the concrete factor-correction phase.  Every premise
is exactly one semantic fact supplied by a generated raw operation: the
integer error quotient, modular long division, and the `tau` construction.
There is no existentially selected correction polynomial. -/
theorem henselFactorCorrection_algebra
    (m : Nat) (hm : 0 < m)
    (F G H S T E Q R Tau : Polynomial (ZMod (m ^ 2)))
    (herror : F = G * H + Polynomial.C (m : ZMod (m ^ 2)) * E)
    (hbezout : Polynomial.map
      (ZMod.castHom (dvd_pow_self m (by omega : 2 ≠ 0)) (ZMod m))
      (S * G + T * H) = 1)
    (hdivmod : Polynomial.map
      (ZMod.castHom (dvd_pow_self m (by omega : 2 ≠ 0)) (ZMod m))
      (R + Q * H) =
      Polynomial.map
        (ZMod.castHom (dvd_pow_self m (by omega : 2 ≠ 0)) (ZMod m))
        (S * E))
    (htau : Polynomial.map
      (ZMod.castHom (dvd_pow_self m (by omega : 2 ≠ 0)) (ZMod m)) Tau =
      Polynomial.map
        (ZMod.castHom (dvd_pow_self m (by omega : 2 ≠ 0)) (ZMod m))
        (T * E + Q * G)) :
    (G + Polynomial.C (m : ZMod (m ^ 2)) * Tau) *
        (H + Polynomial.C (m : ZMod (m ^ 2)) * R) = F := by
  let π := ZMod.castHom (dvd_pow_self m (by omega : 2 ≠ 0)) (ZMod m)
  let c : Polynomial (ZMod (m ^ 2)) :=
    Polynomial.C (m : ZMod (m ^ 2))
  have hcMap : Polynomial.map π c = 0 := by
    ext degree
    simp [π, c, ZMod.castHom_apply]
  have hcc : c * c = 0 :=
    C_scale_mul_eq_zero_of_map_eq_zero m hm c hcMap
  let delta := G * R + Tau * H - E
  have hdeltaMap : Polynomial.map π delta = 0 := by
    simp only [delta, Polynomial.map_sub, Polynomial.map_add,
      Polynomial.map_mul]
    have hbezout' :
        Polynomial.map π S * Polynomial.map π G +
            Polynomial.map π T * Polynomial.map π H = 1 := by
      simpa [π, Polynomial.map_add, Polynomial.map_mul] using hbezout
    have hdivmod' :
        Polynomial.map π R + Polynomial.map π Q * Polynomial.map π H =
          Polynomial.map π S * Polynomial.map π E := by
      simpa [π, Polynomial.map_add, Polynomial.map_mul] using hdivmod
    have htau' : Polynomial.map π Tau =
        Polynomial.map π T * Polynomial.map π E +
          Polynomial.map π Q * Polynomial.map π G := by
      simpa [π, Polynomial.map_add, Polynomial.map_mul] using htau
    rw [htau']
    calc
      Polynomial.map π G * Polynomial.map π R +
            (Polynomial.map π T * Polynomial.map π E +
              Polynomial.map π Q * Polynomial.map π G) *
              Polynomial.map π H - Polynomial.map π E =
          Polynomial.map π G *
              (Polynomial.map π R +
                Polynomial.map π Q * Polynomial.map π H) +
            Polynomial.map π T * Polynomial.map π H *
              Polynomial.map π E - Polynomial.map π E := by ring
      _ = (Polynomial.map π S * Polynomial.map π G +
              Polynomial.map π T * Polynomial.map π H - 1) *
            Polynomial.map π E := by rw [hdivmod']; ring
      _ = 0 := by rw [hbezout']; simp
  have hcdelta : c * delta = 0 :=
    C_scale_mul_eq_zero_of_map_eq_zero m hm delta hdeltaMap
  change (G + c * Tau) * (H + c * R) = F
  calc
    (G + c * Tau) * (H + c * R) =
        G * H + c * delta + c * E + (c * c) * (Tau * R) := by
          simp only [delta]
          ring
    _ = G * H + c * E := by rw [hcdelta, hcc]; simp
    _ = F := herror.symm

/-- Algebraic core of the second generated phase.  It lifts the stored
Bezout identity from `m` to `m²` using exactly the correction returned by the
second concrete modular divmod. -/
theorem henselBezoutCorrection_algebra
    (m : Nat) (hm : 0 < m)
    (G H S T E Q R Tau : Polynomial (ZMod (m ^ 2)))
    (herror : (1 : Polynomial (ZMod (m ^ 2))) =
      S * G + T * H + Polynomial.C (m : ZMod (m ^ 2)) * E)
    (hbezout : Polynomial.map
      (ZMod.castHom (dvd_pow_self m (by omega : 2 ≠ 0)) (ZMod m))
      (S * G + T * H) = 1)
    (hdivmod : Polynomial.map
      (ZMod.castHom (dvd_pow_self m (by omega : 2 ≠ 0)) (ZMod m))
      (R + Q * H) =
      Polynomial.map
        (ZMod.castHom (dvd_pow_self m (by omega : 2 ≠ 0)) (ZMod m))
        (S * E))
    (htau : Polynomial.map
      (ZMod.castHom (dvd_pow_self m (by omega : 2 ≠ 0)) (ZMod m)) Tau =
      Polynomial.map
        (ZMod.castHom (dvd_pow_self m (by omega : 2 ≠ 0)) (ZMod m))
        (T * E + Q * G)) :
    (S + Polynomial.C (m : ZMod (m ^ 2)) * R) * G +
      (T + Polynomial.C (m : ZMod (m ^ 2)) * Tau) * H = 1 := by
  let π := ZMod.castHom (dvd_pow_self m (by omega : 2 ≠ 0)) (ZMod m)
  let c : Polynomial (ZMod (m ^ 2)) :=
    Polynomial.C (m : ZMod (m ^ 2))
  let delta := R * G + Tau * H - E
  have hdeltaMap : Polynomial.map π delta = 0 := by
    simp only [delta, Polynomial.map_sub, Polynomial.map_add,
      Polynomial.map_mul]
    have hbezout' :
        Polynomial.map π S * Polynomial.map π G +
            Polynomial.map π T * Polynomial.map π H = 1 := by
      simpa [π, Polynomial.map_add, Polynomial.map_mul] using hbezout
    have hdivmod' :
        Polynomial.map π R + Polynomial.map π Q * Polynomial.map π H =
          Polynomial.map π S * Polynomial.map π E := by
      simpa [π, Polynomial.map_add, Polynomial.map_mul] using hdivmod
    have htau' : Polynomial.map π Tau =
        Polynomial.map π T * Polynomial.map π E +
          Polynomial.map π Q * Polynomial.map π G := by
      simpa [π, Polynomial.map_add, Polynomial.map_mul] using htau
    rw [htau']
    calc
      Polynomial.map π R * Polynomial.map π G +
            (Polynomial.map π T * Polynomial.map π E +
              Polynomial.map π Q * Polynomial.map π G) *
              Polynomial.map π H - Polynomial.map π E =
          Polynomial.map π G *
              (Polynomial.map π R +
                Polynomial.map π Q * Polynomial.map π H) +
            Polynomial.map π T * Polynomial.map π H *
              Polynomial.map π E - Polynomial.map π E := by ring
      _ = (Polynomial.map π S * Polynomial.map π G +
              Polynomial.map π T * Polynomial.map π H - 1) *
            Polynomial.map π E := by rw [hdivmod']; ring
      _ = 0 := by rw [hbezout']; simp
  have hcdelta : c * delta = 0 :=
    C_scale_mul_eq_zero_of_map_eq_zero m hm delta hdeltaMap
  change (S + c * R) * G + (T + c * Tau) * H = 1
  calc
    (S + c * R) * G + (T + c * Tau) * H =
        S * G + T * H + c * E + c * delta := by
          simp only [delta]
          ring
    _ = S * G + T * H + c * E := by rw [hcdelta, add_zero]
    _ = 1 := herror.symm

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

set_option maxHeartbeats 0 in
/-- Semantic invariant of the generated C++ `pair_vec_add` iterator merge.
The output accumulator plus both unconsumed suffixes is preserved in every
source branch. -/
theorem pairVecAddLoop_toPolyMod (m : Nat) (a b : SparsePolyZZ) :
    ∀ (aIndex bIndex : Nat) (result : SparsePolyZZ),
      toPolyMod m
          (Generated.StrictHensel.pairVecAddLoop a b aIndex bIndex result) =
        toPolyMod m result + termsToPolyMod m (a.toList.drop aIndex) +
          termsToPolyMod m (b.toList.drop bIndex) := by
  intro aIndex bIndex result
  refine Generated.StrictHensel.pairVecAddLoop.induct a b
    (motive := fun aIndex bIndex result =>
      toPolyMod m
          (Generated.StrictHensel.pairVecAddLoop a b aIndex bIndex result) =
        toPolyMod m result + termsToPolyMod m (a.toList.drop aIndex) +
          termsToPolyMod m (b.toList.drop bIndex)) ?_ ?_ ?_ ?_ ?_ ?_
    aIndex bIndex result
  · intro ai bi acc hmore haDone ih
    have hbi : bi < b.size := by omega
    have hbDrop := drop_eq_getElem_cons b bi hbi
    have haDrop : a.toList.drop ai = [] := by
      apply List.drop_eq_nil_of_le
      simpa using haDone
    rw [Generated.StrictHensel.pairVecAddLoop.eq_1]
    simp [hmore, haDone]
    rw [ih, toPolyMod_push, haDrop, hbDrop, termsToPolyMod_cons,
      termsToPolyMod_nil]
    ring
  · intro ai bi acc hmore haMore hbDone ih
    have hai : ai < a.size := by omega
    have haDrop := drop_eq_getElem_cons a ai hai
    have hbDrop : b.toList.drop bi = [] := by
      apply List.drop_eq_nil_of_le
      simpa using hbDone
    rw [Generated.StrictHensel.pairVecAddLoop.eq_1]
    simp [hmore, haMore, hbDone]
    rw [ih, toPolyMod_push, haDrop, hbDrop, termsToPolyMod_cons,
      termsToPolyMod_nil]
    ring
  · intro ai bi acc hmore haMore hbMore hdegree ih
    have hbi : bi < b.size := by omega
    have hbDrop := drop_eq_getElem_cons b bi hbi
    rw [Generated.StrictHensel.pairVecAddLoop.eq_1]
    simp [hmore, haMore, hbMore, hdegree]
    rw [ih, toPolyMod_push, hbDrop, termsToPolyMod_cons]
    ring
  · intro ai bi acc hmore haMore hbMore hnotGreater hequal ih
    have hai : ai < a.size := by omega
    have hbi : bi < b.size := by omega
    have haDrop := drop_eq_getElem_cons a ai hai
    have hbDrop := drop_eq_getElem_cons b bi hbi
    rw [Generated.StrictHensel.pairVecAddLoop.eq_1]
    simp [hmore, haMore, hbMore, hnotGreater, hequal]
    rw [ih, pushNonzero_toPolyMod, haDrop, hbDrop,
      termsToPolyMod_cons]
    simp only [Int.cast_add, termsToPolyMod_cons]
    rw [hequal]
    rw [Polynomial.monomial_add]
    ring
  · intro ai bi acc hmore haMore hbMore hnotGreater hnotEqual ih
    have hai : ai < a.size := by omega
    have haDrop := drop_eq_getElem_cons a ai hai
    rw [Generated.StrictHensel.pairVecAddLoop.eq_1]
    simp [hmore, haMore, hbMore, hnotGreater, hnotEqual]
    rw [ih, toPolyMod_push, haDrop, termsToPolyMod_cons]
    ring
  · intro ai bi acc hdone
    have haDrop : a.toList.drop ai = [] := by
      apply List.drop_eq_nil_of_le
      simp only [not_or] at hdone
      simpa using hdone.1
    have hbDrop : b.toList.drop bi = [] := by
      apply List.drop_eq_nil_of_le
      simp only [not_or] at hdone
      simpa using hdone.2
    rw [Generated.StrictHensel.pairVecAddLoop.eq_1]
    simp [hdone, haDrop, hbDrop]

theorem __upoly_add_raw_ir_refines (m : Nat) (a b : SparsePolyZZ) :
    ∃ output,
      Generated.StrictHensel.__upoly_add_raw_ir a b = .ok output ∧
      toPolyMod m output = toPolyMod m a + toPolyMod m b := by
  refine ⟨Generated.StrictHensel.pairVecAddLoop a b 0 0 #[], rfl, ?_⟩
  rw [pairVecAddLoop_toPolyMod]
  simp [toPolyMod_eq_termsToPolyMod]

set_option maxHeartbeats 0 in
/-- Semantic invariant of the generated C++ `pair_vec_sub` iterator merge. -/
theorem pairVecSubLoop_toPolyMod (m : Nat) (a b : SparsePolyZZ) :
    ∀ (aIndex bIndex : Nat) (result : SparsePolyZZ),
      toPolyMod m
          (Generated.StrictHensel.pairVecSubLoop a b aIndex bIndex result) =
        toPolyMod m result + termsToPolyMod m (a.toList.drop aIndex) -
          termsToPolyMod m (b.toList.drop bIndex) := by
  intro aIndex bIndex result
  refine Generated.StrictHensel.pairVecSubLoop.induct a b
    (motive := fun aIndex bIndex result =>
      toPolyMod m
          (Generated.StrictHensel.pairVecSubLoop a b aIndex bIndex result) =
        toPolyMod m result + termsToPolyMod m (a.toList.drop aIndex) -
          termsToPolyMod m (b.toList.drop bIndex)) ?_ ?_ ?_ ?_ ?_ ?_
    aIndex bIndex result
  · intro ai bi acc hmore haDone ih
    have hbi : bi < b.size := by omega
    have hbDrop := drop_eq_getElem_cons b bi hbi
    have haDrop : a.toList.drop ai = [] := by
      apply List.drop_eq_nil_of_le
      simpa using haDone
    rw [Generated.StrictHensel.pairVecSubLoop.eq_1]
    simp [hmore, haDone]
    rw [ih, toPolyMod_push, haDrop, hbDrop, termsToPolyMod_cons,
      termsToPolyMod_nil]
    simp only [Int.cast_neg, Prod.fst, Prod.snd]
    rw [Polynomial.monomial_neg]
    ring
  · intro ai bi acc hmore haMore hbDone ih
    have hai : ai < a.size := by omega
    have haDrop := drop_eq_getElem_cons a ai hai
    have hbDrop : b.toList.drop bi = [] := by
      apply List.drop_eq_nil_of_le
      simpa using hbDone
    rw [Generated.StrictHensel.pairVecSubLoop.eq_1]
    simp [hmore, haMore, hbDone]
    rw [ih, toPolyMod_push, haDrop, hbDrop, termsToPolyMod_cons,
      termsToPolyMod_nil]
    ring
  · intro ai bi acc hmore haMore hbMore hdegree ih
    have hbi : bi < b.size := by omega
    have hbDrop := drop_eq_getElem_cons b bi hbi
    rw [Generated.StrictHensel.pairVecSubLoop.eq_1]
    simp [hmore, haMore, hbMore, hdegree]
    rw [ih, toPolyMod_push, hbDrop, termsToPolyMod_cons]
    simp only [Int.cast_neg, Prod.fst, Prod.snd]
    rw [Polynomial.monomial_neg]
    ring
  · intro ai bi acc hmore haMore hbMore hnotGreater hequal ih
    have hai : ai < a.size := by omega
    have hbi : bi < b.size := by omega
    have haDrop := drop_eq_getElem_cons a ai hai
    have hbDrop := drop_eq_getElem_cons b bi hbi
    rw [Generated.StrictHensel.pairVecSubLoop.eq_1]
    simp [hmore, haMore, hbMore, hnotGreater, hequal]
    rw [ih, pushNonzero_toPolyMod, haDrop, hbDrop,
      termsToPolyMod_cons]
    simp only [Int.cast_sub, termsToPolyMod_cons]
    rw [hequal, Polynomial.monomial_sub]
    ring
  · intro ai bi acc hmore haMore hbMore hnotGreater hnotEqual ih
    have hai : ai < a.size := by omega
    have haDrop := drop_eq_getElem_cons a ai hai
    rw [Generated.StrictHensel.pairVecSubLoop.eq_1]
    simp [hmore, haMore, hbMore, hnotGreater, hnotEqual]
    rw [ih, toPolyMod_push, haDrop, termsToPolyMod_cons]
    ring
  · intro ai bi acc hdone
    have haDrop : a.toList.drop ai = [] := by
      apply List.drop_eq_nil_of_le
      simp only [not_or] at hdone
      simpa using hdone.1
    have hbDrop : b.toList.drop bi = [] := by
      apply List.drop_eq_nil_of_le
      simp only [not_or] at hdone
      simpa using hdone.2
    rw [Generated.StrictHensel.pairVecSubLoop.eq_1]
    simp [hdone, haDrop, hbDrop]

theorem __upoly_sub_raw_ir_refines (m : Nat) (a b : SparsePolyZZ) :
    ∃ output,
      Generated.StrictHensel.__upoly_sub_raw_ir a b = .ok output ∧
      toPolyMod m output = toPolyMod m a - toPolyMod m b := by
  refine ⟨Generated.StrictHensel.pairVecSubLoop a b 0 0 #[], rfl, ?_⟩
  rw [pairVecSubLoop_toPolyMod]
  simp [toPolyMod_eq_termsToPolyMod]

noncomputable def mulProductsToPolyMod (m : Nat) :
    List Generated.StrictHensel.MulProduct → Polynomial (ZMod m)
  | [] => 0
  | product :: products =>
      Polynomial.monomial product.degree (product.coefficient : ZMod m) +
        mulProductsToPolyMod m products

theorem mulProductsToPolyMod_append (m : Nat)
    (left right : List Generated.StrictHensel.MulProduct) :
    mulProductsToPolyMod m (left ++ right) =
      mulProductsToPolyMod m left + mulProductsToPolyMod m right := by
  induction left with
  | nil => simp [mulProductsToPolyMod]
  | cons product products ih =>
      simp [mulProductsToPolyMod, ih, add_assoc]

/-- The `addmul` accumulator for one heap key is exactly the polynomial sum
of all frontier contributions with that degree. -/
theorem mulDegreeCoefficient_toPolyMod (m degree : Nat)
    (products : List Generated.StrictHensel.MulProduct) :
    Polynomial.monomial degree
        (Generated.StrictHensel.mulDegreeCoefficient degree products : ZMod m) =
      mulProductsToPolyMod m
        (products.filter fun product => product.degree = degree) := by
  induction products with
  | nil => simp [mulProductsToPolyMod,
      Generated.StrictHensel.mulDegreeCoefficient]
  | cons product products ih =>
      by_cases hdegree : product.degree = degree
      · simp [hdegree, mulProductsToPolyMod,
          Generated.StrictHensel.mulDegreeCoefficient,
          Polynomial.monomial_add, ih]
      · simp [hdegree, mulProductsToPolyMod,
          Generated.StrictHensel.mulDegreeCoefficient, ih]

/-- Partitioning the heap frontier at one degree preserves its represented
polynomial. -/
theorem mulProducts_partition (m degree : Nat)
    (products : List Generated.StrictHensel.MulProduct) :
    mulProductsToPolyMod m products =
      mulProductsToPolyMod m
          (products.filter fun product => product.degree = degree) +
        mulProductsToPolyMod m
          (products.filter fun product => product.degree != degree) := by
  induction products with
  | nil => simp [mulProductsToPolyMod]
  | cons product products ih =>
      by_cases hdegree : product.degree = degree
      · simp [mulProductsToPolyMod, hdegree]
        rw [ih]
        ring
      · simp [mulProductsToPolyMod, hdegree]
        rw [ih]
        ring

set_option maxHeartbeats 0 in
/-- Conservation law of the generated multiplication heap loop. -/
theorem pairVecMulHeapLoop_toPolyMod (m : Nat) :
    ∀ (products : List Generated.StrictHensel.MulProduct)
      (result : SparsePolyZZ),
      toPolyMod m
          (Generated.StrictHensel.pairVecMulHeapLoop products result) =
        toPolyMod m result + mulProductsToPolyMod m products := by
  intro products result
  refine Generated.StrictHensel.pairVecMulHeapLoop.induct
    (motive := fun products result =>
      toPolyMod m
          (Generated.StrictHensel.pairVecMulHeapLoop products result) =
        toPolyMod m result + mulProductsToPolyMod m products) ?_ ?_
    products result
  · intro result
    rw [Generated.StrictHensel.pairVecMulHeapLoop.eq_1]
    simp [mulProductsToPolyMod]
  · intro products result hempty
    dsimp only
    intro ih
    rw [Generated.StrictHensel.pairVecMulHeapLoop.eq_1]
    simp [hempty]
    have ih' :
        toPolyMod m
            (Generated.StrictHensel.pairVecMulHeapLoop
              (products.filter fun product =>
                product.degree != Generated.StrictHensel.mulMaxDegree products)
              (Generated.StrictHensel.pushNonzero result
                (Generated.StrictHensel.mulMaxDegree products)
                (Generated.StrictHensel.mulDegreeCoefficient
                  (Generated.StrictHensel.mulMaxDegree products) products))) =
          toPolyMod m
              (Generated.StrictHensel.pushNonzero result
                (Generated.StrictHensel.mulMaxDegree products)
                (Generated.StrictHensel.mulDegreeCoefficient
                  (Generated.StrictHensel.mulMaxDegree products) products)) +
            mulProductsToPolyMod m
              (products.filter fun product =>
                product.degree != Generated.StrictHensel.mulMaxDegree products) := by
      simpa using ih
    rw [ih', pushNonzero_toPolyMod,
      mulDegreeCoefficient_toPolyMod]
    rw [mulProducts_partition m
      (Generated.StrictHensel.mulMaxDegree products) products]
    ring

private theorem mulProductsForTerm_toPolyMod (m : Nat)
    (term : UMonomial × Int) (terms : List (UMonomial × Int)) :
    mulProductsToPolyMod m
        (terms.map fun other : UMonomial × Int =>
          Generated.StrictHensel.MulProduct.mk
            (term.1.deg + other.1.deg) (term.2 * other.2)) =
      Polynomial.monomial term.1.deg (term.2 : ZMod m) *
        termsToPolyMod m terms := by
  induction terms with
  | nil => simp [mulProductsToPolyMod]
  | cons other terms ih =>
      simp only [List.map_cons, mulProductsToPolyMod,
        termsToPolyMod_cons, Int.cast_mul]
      rw [ih, mul_add, Polynomial.monomial_mul_monomial]

theorem pairVecMulProducts_toPolyMod (m : Nat) (a b : SparsePolyZZ) :
    mulProductsToPolyMod m
        (Generated.StrictHensel.pairVecMulProducts a b) =
      toPolyMod m a * toPolyMod m b := by
  rw [toPolyMod_eq_termsToPolyMod, toPolyMod_eq_termsToPolyMod]
  unfold Generated.StrictHensel.pairVecMulProducts
  induction a.toList with
  | nil => simp [mulProductsToPolyMod]
  | cons term terms ih =>
      simp only [List.flatMap_cons, termsToPolyMod_cons]
      rw [show mulProductsToPolyMod m
          ((b.toList.map fun tb =>
              Generated.StrictHensel.MulProduct.mk
                (term.1.deg + tb.1.deg) (term.2 * tb.2)) ++
            terms.flatMap fun ta =>
              b.toList.map fun tb =>
                Generated.StrictHensel.MulProduct.mk
                  (ta.1.deg + tb.1.deg) (ta.2 * tb.2)) =
          mulProductsToPolyMod m
              (b.toList.map fun tb =>
                Generated.StrictHensel.MulProduct.mk
                  (term.1.deg + tb.1.deg) (term.2 * tb.2)) +
            mulProductsToPolyMod m
              (terms.flatMap fun ta =>
                b.toList.map fun tb =>
                  Generated.StrictHensel.MulProduct.mk
                    (ta.1.deg + tb.1.deg) (ta.2 * tb.2)) by
        exact mulProductsToPolyMod_append m _ _]
      rw [mulProductsForTerm_toPolyMod, ih, add_mul]

theorem __upoly_mul_raw_ir_refines (m : Nat) (a b : SparsePolyZZ) :
    ∃ output,
      Generated.StrictHensel.__upoly_mul_raw_ir a b = .ok output ∧
      toPolyMod m output = toPolyMod m a * toPolyMod m b := by
  refine ⟨Generated.StrictHensel.pairVecMulHeapLoop
      (Generated.StrictHensel.pairVecMulProducts a b) #[], rfl, ?_⟩
  rw [pairVecMulHeapLoop_toPolyMod, pairVecMulProducts_toPolyMod]
  rw [toPolyMod_empty, zero_add]

/-- The sole concrete operation bundle used by the strict Hensel entry.  All
three arithmetic fields execute generated C++ L1 control flow; none points to
the `SparsePolyZZ` model instances or an L2 polynomial operation. -/
def strictHenselRawOps
    (termination : Generated.StrictHensel.DivmodTermination) :
    Generated.StrictHensel.HenselStepRawOps where
  mul := Generated.StrictHensel.__upoly_mul_raw_ir
  add := Generated.StrictHensel.__upoly_add_raw_ir
  sub := Generated.StrictHensel.__upoly_sub_raw_ir
  divmodTermination := termination

theorem strictHenselRawOps_mul_refines
    (termination : Generated.StrictHensel.DivmodTermination)
    (m : Nat) (a b : SparsePolyZZ) :
    ∃ output,
      (strictHenselRawOps termination).mul a b = .ok output ∧
      toPolyMod m output = toPolyMod m a * toPolyMod m b :=
  __upoly_mul_raw_ir_refines m a b

theorem strictHenselRawOps_add_refines
    (termination : Generated.StrictHensel.DivmodTermination)
    (m : Nat) (a b : SparsePolyZZ) :
    ∃ output,
      (strictHenselRawOps termination).add a b = .ok output ∧
      toPolyMod m output = toPolyMod m a + toPolyMod m b :=
  __upoly_add_raw_ir_refines m a b

theorem strictHenselRawOps_sub_refines
    (termination : Generated.StrictHensel.DivmodTermination)
    (m : Nat) (a b : SparsePolyZZ) :
    ∃ output,
      (strictHenselRawOps termination).sub a b = .ok output ∧
      toPolyMod m output = toPolyMod m a - toPolyMod m b :=
  __upoly_sub_raw_ir_refines m a b

/-- Deterministic, explicit-output forms used when threading a successful raw
execution through the generated `do` sequence. -/
theorem strictHenselRawOps_mul_refines_of_run
    (termination : Generated.StrictHensel.DivmodTermination)
    (m : Nat) (a b output : SparsePolyZZ)
    (hrun : (strictHenselRawOps termination).mul a b = .ok output) :
    toPolyMod m output = toPolyMod m a * toPolyMod m b := by
  rcases strictHenselRawOps_mul_refines termination m a b with
    ⟨actual, hactual, hsemantic⟩
  rw [hrun] at hactual
  cases hactual
  exact hsemantic

theorem strictHenselRawOps_add_refines_of_run
    (termination : Generated.StrictHensel.DivmodTermination)
    (m : Nat) (a b output : SparsePolyZZ)
    (hrun : (strictHenselRawOps termination).add a b = .ok output) :
    toPolyMod m output = toPolyMod m a + toPolyMod m b := by
  rcases strictHenselRawOps_add_refines termination m a b with
    ⟨actual, hactual, hsemantic⟩
  rw [hrun] at hactual
  cases hactual
  exact hsemantic

theorem strictHenselRawOps_sub_refines_of_run
    (termination : Generated.StrictHensel.DivmodTermination)
    (m : Nat) (a b output : SparsePolyZZ)
    (hrun : (strictHenselRawOps termination).sub a b = .ok output) :
    toPolyMod m output = toPolyMod m a - toPolyMod m b := by
  rcases strictHenselRawOps_sub_refines termination m a b with
    ⟨actual, hactual, hsemantic⟩
  rw [hrun] at hactual
  cases hactual
  exact hsemantic

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

/-- Every concrete exponent fits the signed degree range used by the C++
long-division loop. -/
def DegreesBound (f : SparsePolyZZ) : Prop :=
  ∀ term ∈ f.toList, term.1.deg < 2 ^ 63

/-- The first array term is a genuine maximum exponent, as required by the
source sparse-polynomial representation. -/
def HeadDominates (f : SparsePolyZZ) : Prop :=
  ∀ term ∈ f.toList, term.1.deg ≤ f[0]!.1.deg

private theorem degreesBound_pushNonzero (f : SparsePolyZZ)
    (degree : Nat) (coefficient : Int) (hf : DegreesBound f)
    (hdegree : degree < 2 ^ 63) :
    DegreesBound (Generated.StrictHensel.pushNonzero f degree coefficient) := by
  intro term hterm
  by_cases hzero : coefficient = 0
  · simp [Generated.StrictHensel.pushNonzero, hzero] at hterm
    exact hf term (by simpa using hterm)
  · simp [Generated.StrictHensel.pushNonzero, hzero] at hterm
    rcases hterm with hterm | rfl
    · exact hf term (by simpa using hterm)
    · exact hdegree

/-- Degree-range preservation for all six branches of the concrete merge.
The assumptions describe exactly the two input suffix sources: raw remainder
terms and shifted divisor terms. -/
theorem divmodMergeLoop_degreesBound (r g : SparsePolyZZ)
    (coefficient m : Int) (degreeShift : Nat)
    (hr : DegreesBound r)
    (hgShift : ∀ term ∈ g.toList,
      term.1.deg + degreeShift < 2 ^ 63) :
    ∀ (rIndex gIndex : Nat) (result : SparsePolyZZ),
      DegreesBound result →
      DegreesBound (Generated.StrictHensel.divmodMergeLoop r g coefficient m
        degreeShift rIndex gIndex result) := by
  intro rIndex gIndex result
  refine Generated.StrictHensel.divmodMergeLoop.induct r g coefficient m
    degreeShift
    (motive := fun rIndex gIndex result => DegreesBound result →
      DegreesBound (Generated.StrictHensel.divmodMergeLoop r g coefficient m
        degreeShift rIndex gIndex result)) ?_ ?_ ?_ ?_ ?_ ?_
    rIndex gIndex result
  · intro ri gi acc hmore hgDone
    dsimp only
    intro ih hacc
    have hri : ri < r.size := by omega
    have hrmem : r[ri]! ∈ r.toList := by
      rw [getElem!_pos r ri hri]
      exact Array.getElem_mem_toList hri
    rw [Generated.StrictHensel.divmodMergeLoop.eq_1]
    simp [hmore, hgDone]
    exact ih (degreesBound_pushNonzero acc r[ri]!.1.deg _ hacc
      (hr r[ri]! hrmem))
  · intro ri gi acc hmore hgMore hrDone
    dsimp only
    intro ih hacc
    have hgi : gi < g.size := by omega
    have hgmem : g[gi]! ∈ g.toList := by
      rw [getElem!_pos g gi hgi]
      exact Array.getElem_mem_toList hgi
    rw [Generated.StrictHensel.divmodMergeLoop.eq_1]
    simp [hmore, hgMore, hrDone]
    simp only [Int.emod_emod] at ih
    exact ih (degreesBound_pushNonzero acc
      (g[gi]!.1.deg + degreeShift) _ hacc (hgShift g[gi]! hgmem))
  · intro ri gi acc hmore hgMore hrMore
    dsimp only
    intro hdegree ih hacc
    have hri : ri < r.size := by omega
    have hrmem : r[ri]! ∈ r.toList := by
      rw [getElem!_pos r ri hri]
      exact Array.getElem_mem_toList hri
    rw [Generated.StrictHensel.divmodMergeLoop.eq_1]
    simp [hmore, hgMore, hrMore, hdegree]
    exact ih (degreesBound_pushNonzero acc r[ri]!.1.deg _ hacc
      (hr r[ri]! hrmem))
  · intro ri gi acc hmore hgMore hrMore
    dsimp only
    intro hnotGreater hless ih hacc
    have hgi : gi < g.size := by omega
    have hgmem : g[gi]! ∈ g.toList := by
      rw [getElem!_pos g gi hgi]
      exact Array.getElem_mem_toList hgi
    rw [Generated.StrictHensel.divmodMergeLoop.eq_1]
    simp [hmore, hgMore, hrMore, hnotGreater, hless]
    simp only [Int.emod_emod] at ih
    exact ih (degreesBound_pushNonzero acc
      (g[gi]!.1.deg + degreeShift) _ hacc (hgShift g[gi]! hgmem))
  · intro ri gi acc hmore hgMore hrMore
    dsimp only
    intro hnotGreater hnotLess ih hacc
    have hri : ri < r.size := by omega
    have hrmem : r[ri]! ∈ r.toList := by
      rw [getElem!_pos r ri hri]
      exact Array.getElem_mem_toList hri
    rw [Generated.StrictHensel.divmodMergeLoop.eq_1]
    simp [hmore, hgMore, hrMore, hnotGreater, hnotLess]
    exact ih (degreesBound_pushNonzero acc r[ri]!.1.deg _ hacc
      (hr r[ri]! hrmem))
  · intro ri gi acc hdone hacc
    rw [Generated.StrictHensel.divmodMergeLoop.eq_1]
    simpa [hdone] using hacc

theorem divmodRemainder_degreesBound (r g : SparsePolyZZ)
    (coefficient m : Int) (hr : 0 < r.size) (hg : 0 < g.size)
    (hrBound : DegreesBound r) (hgHead : HeadDominates g)
    (hrDegree : r[0]!.1.deg < 2 ^ 63)
    (hgDegree : g[0]!.1.deg < 2 ^ 63)
    (hactive : get_deg r ≥ get_deg g) :
    DegreesBound (Generated.StrictHensel.divmodRemainder r g coefficient m
      (get_deg r - get_deg g).toNatClampNeg) := by
  have hdegree := get_deg_sub_toNatClampNeg_eq_shift r g hr hg hrDegree
    hgDegree hactive
  have hgShift : ∀ term ∈ g.toList,
      term.1.deg + (get_deg r - get_deg g).toNatClampNeg < 2 ^ 63 := by
    intro term hterm
    have hle := hgHead term hterm
    omega
  rw [Generated.StrictHensel.divmodRemainder]
  apply divmodMergeLoop_degreesBound r g coefficient m _ hrBound hgShift
  simp [DegreesBound]

theorem degreesBound_eraseIdxIfInBounds (f : SparsePolyZZ) (i : Nat)
    (hf : DegreesBound f) : DegreesBound (f.eraseIdxIfInBounds i) := by
  intro term hterm
  rw [Array.toList_eraseIdxIfInBounds] at hterm
  exact hf term (List.mem_of_mem_eraseIdx hterm)

theorem modCoeffOutput_degreesBound (f : SparsePolyZZ) (m : Int)
    (hf : DegreesBound f) :
    DegreesBound (Generated.StrictHensel.modCoeffOutput f m) := by
  intro term hterm
  simp only [Generated.StrictHensel.modCoeffOutput,
    Array.toList_filterMap, List.mem_filterMap] at hterm
  rcases hterm with ⟨source, hsource, hterm⟩
  split at hterm
  · injection hterm with hterm
    subst term
    exact hf source (by simpa using hsource)
  · simp_all

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

/-- The machine-degree certificate is derived from the concrete trace and
the ordinary sparse-representation bounds.  In particular, neither the trace
provider nor the public refinement theorem needs to supply this certificate
as an independent assumption. -/
theorem divmodTrace_degreeSafe (m : Nat) (g : SparsePolyZZ)
    (hg : 0 < g.size) (hgDegree : g[0]!.1.deg < 2 ^ 63)
    (hgHead : HeadDominates g) :
    ∀ {r q : SparsePolyZZ}
      (trace : Generated.StrictHensel.DivmodTrace g
        (ZZ.invert 0 g[0]!.2 (m : Int)).2 (m : Int) r q),
      DegreesBound r →
      DivmodTraceDegreeSafe g
        (ZZ.invert 0 g[0]!.2 (m : Int)).2 m trace := by
  intro r q trace hrBound
  induction trace with
  | done r q inactive =>
      trivial
  | vanished r q active zero next ih =>
      have hactive : r.isEmpty = false ∧ get_deg r ≥ get_deg g := by
        simpa using active
      have hr : 0 < r.size := by
        have hrne : r ≠ #[] := by
          simpa [Array.isEmpty_iff] using hactive.1
        have : r.size ≠ 0 := by
          intro hzero
          apply hrne
          exact Array.eq_empty_of_size_eq_zero hzero
        omega
      have hrmem : r[0]! ∈ r.toList := by
        rw [getElem!_pos r 0 hr]
        exact Array.getElem_mem_toList hr
      simp only [DivmodTraceDegreeSafe]
      exact ⟨hr, hrBound r[0]! hrmem, hactive.2,
        ih (degreesBound_eraseIdxIfInBounds r 0 hrBound)⟩
  | subtract r q active nonzero next ih =>
      have hactive : r.isEmpty = false ∧ get_deg r ≥ get_deg g := by
        simpa using active
      have hr : 0 < r.size := by
        have hrne : r ≠ #[] := by
          simpa [Array.isEmpty_iff] using hactive.1
        have : r.size ≠ 0 := by
          intro hzero
          apply hrne
          exact Array.eq_empty_of_size_eq_zero hzero
        omega
      have hrmem : r[0]! ∈ r.toList := by
        rw [getElem!_pos r 0 hr]
        exact Array.getElem_mem_toList hr
      have hrDegree := hrBound r[0]! hrmem
      simp only [DivmodTraceDegreeSafe]
      exact ⟨hr, hrDegree, hactive.2, ih
        (divmodRemainder_degreesBound r g
          (Generated.StrictHensel.divmodCoefficient r
            (ZZ.invert 0 g[0]!.2 (m : Int)).2 (m : Int))
          (m : Int) hr hg hrBound hgHead hrDegree hgDegree hactive.2)⟩
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

/-- Publicly usable form of trace conservation: the trace safety certificate
is reconstructed from the concrete initial remainder representation. -/
theorem divmodLoop_preserves_of_bounds (m : Nat) (g : SparsePolyZZ)
    (hg : 0 < g.size) (hgDegree : g[0]!.1.deg < 2 ^ 63)
    (hgHead : HeadDominates g)
    (hinvert : (ZZ.invert 0 g[0]!.2 (m : Int)).1 = true)
    {r q : SparsePolyZZ}
    (trace : Generated.StrictHensel.DivmodTrace g
      (ZZ.invert 0 g[0]!.2 (m : Int)).2 (m : Int) r q)
    (hrBound : DegreesBound r) :
    let output := Generated.StrictHensel.divmodLoop g
      (ZZ.invert 0 g[0]!.2 (m : Int)).2 (m : Int) trace
    toPolyMod m output.2 + toPolyMod m output.1 * toPolyMod m g =
      toPolyMod m r + toPolyMod m q * toPolyMod m g := by
  exact divmodLoop_preserves m g hg hgDegree hinvert trace
    (divmodTrace_degreeSafe m g hg hgDegree hgHead trace hrBound)

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

/-- Reducing coefficients modulo a stronger modulus also preserves the
represented polynomial at every divisor modulus.  In particular, the `m²`
stores performed by `__hensel_step` preserve all four fields modulo `m`. -/
theorem __upoly_mod_coeff_raw_ir_preserves_divisor
    {smaller larger : Nat} (hdiv : smaller ∣ larger)
    (f : SparsePolyZZ) :
    ∃ output,
      Generated.StrictHensel.__upoly_mod_coeff_raw_ir f (larger : Int) =
        .ok output ∧
      toPolyMod smaller output = toPolyMod smaller f := by
  rcases __upoly_mod_coeff_raw_ir_refines f larger with
    ⟨output, hrun, hsem⟩
  refine ⟨output, hrun, ?_⟩
  have projected := congrArg
    (Polynomial.map (ZMod.castHom hdiv (ZMod smaller))) hsem
  rw [map_toPolyMod_of_dvd hdiv, map_toPolyMod_of_dvd hdiv] at projected
  exact projected

theorem __upoly_mod_coeff_raw_ir_refines_of_run
    (f output : SparsePolyZZ) (m : Nat)
    (hrun : Generated.StrictHensel.__upoly_mod_coeff_raw_ir f (m : Int) =
      .ok output) :
    toPolyMod m output = toPolyMod m f := by
  rcases __upoly_mod_coeff_raw_ir_refines f m with
    ⟨actual, hactual, hsemantic⟩
  rw [hrun] at hactual
  cases hactual
  exact hsemantic

theorem __upoly_mod_coeff_raw_ir_preserves_divisor_of_run
    {smaller larger : Nat} (hdiv : smaller ∣ larger)
    (f output : SparsePolyZZ)
    (hrun : Generated.StrictHensel.__upoly_mod_coeff_raw_ir f
      (larger : Int) = .ok output) :
    toPolyMod smaller output = toPolyMod smaller f := by
  rcases __upoly_mod_coeff_raw_ir_preserves_divisor hdiv f with
    ⟨actual, hactual, hsemantic⟩
  rw [hrun] at hactual
  cases hactual
  exact hsemantic

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

/-- Full semantic refinement of the generated C++ modular long-division
entry.  The returned concrete arrays satisfy the quotient/remainder equation;
execution is the exact finite trace and not an L2 division call. -/
theorem __upoly_divmod_mod_raw_ir_refines
    (termination : Generated.StrictHensel.DivmodTermination)
    (f g : SparsePolyZZ) (m : Nat)
    (hg : 0 < g.size) (hgDegree : g[0]!.1.deg < 2 ^ 63)
    (hgHead : HeadDominates g) (hfBound : DegreesBound f)
    (hinvert : (ZZ.invert 0 g[0]!.2 (m : Int)).1 = true) :
    ∃ q r,
      Generated.StrictHensel.__upoly_divmod_mod_raw_ir termination f g
          (m : Int) = .ok (q, r) ∧
      toPolyMod m r + toPolyMod m q * toPolyMod m g = toPolyMod m f := by
  let reduced := Generated.StrictHensel.modCoeffOutput f (m : Int)
  let trace := termination.trace f g reduced (m : Int) (by rfl)
  let output := Generated.StrictHensel.divmodLoop g
    (ZZ.invert 0 g[0]!.2 (m : Int)).2 (m : Int) trace
  have hgFalse : g.isEmpty = false := by
    simpa [Array.isEmpty_iff] using (show g ≠ #[] by
      intro hempty
      subst g
      simp at hg)
  have hrun :
      Generated.StrictHensel.__upoly_divmod_mod_raw_ir termination f g
          (m : Int) = .ok (output.1, output.2) := by
    simp [Generated.StrictHensel.__upoly_divmod_mod_raw_ir, hgFalse,
      hinvert, Generated.StrictHensel.__upoly_mod_coeff_raw_ir,
      reduced, trace, output]
  have hpreserve := divmodLoop_preserves_of_bounds m g hg hgDegree hgHead
    hinvert trace (modCoeffOutput_degreesBound f (m : Int) hfBound)
  have hmod : toPolyMod m reduced = toPolyMod m f := by
    obtain ⟨modOutput, hmodRun, hmodPoly⟩ :=
      __upoly_mod_coeff_raw_ir_refines f m
    have : modOutput = reduced := by
      simpa [Generated.StrictHensel.__upoly_mod_coeff_raw_ir, reduced] using
        hmodRun.symm
    subst modOutput
    exact hmodPoly
  refine ⟨output.1, output.2, hrun, ?_⟩
  simpa [output, trace, hmod] using hpreserve

/-- Explicit-output form of modular long-division refinement. -/
theorem __upoly_divmod_mod_raw_ir_refines_of_run
    (termination : Generated.StrictHensel.DivmodTermination)
    (f g q r : SparsePolyZZ) (m : Nat)
    (hg : 0 < g.size) (hgDegree : g[0]!.1.deg < 2 ^ 63)
    (hgHead : HeadDominates g) (hfBound : DegreesBound f)
    (hinvert : (ZZ.invert 0 g[0]!.2 (m : Int)).1 = true)
    (hrun : Generated.StrictHensel.__upoly_divmod_mod_raw_ir termination
      f g (m : Int) = .ok (q, r)) :
    toPolyMod m r + toPolyMod m q * toPolyMod m g = toPolyMod m f := by
  rcases __upoly_divmod_mod_raw_ir_refines termination f g m hg hgDegree
      hgHead hfBound hinvert with ⟨actualQ, actualR, hactual, hsemantic⟩
  rw [hrun] at hactual
  cases hactual
  exact hsemantic

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

/-- Modular decoding of the same generated coefficient-scaling loop.  The
target modulus is intentionally independent of the integer scale: the
quadratic step uses this once at `m` and once at `m²`. -/
theorem scaleCoeffs_toPolyMod (f : SparsePolyZZ) (scale : ZZ)
    (modulus : Nat) :
    toPolyMod modulus (Generated.StrictHensel.scaleCoeffs f scale) =
      Polynomial.C (scale : ZMod modulus) * toPolyMod modulus f := by
  unfold toPolyMod
  rw [scaleCoeffs_toPoly, Polynomial.map_mul, Polynomial.map_C]
  rfl

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

/-- The exact source error quotient remains exact modulo `m²` after the C++
loop reduces and compacts the quotient coefficients modulo `m`.  This rules
out replacing the generated quotient by an independently chosen L2 witness. -/
theorem divideThenReduceCoeffs_scaled_toPolyMod_sq
    (f : SparsePolyZZ) (m : Nat) (hm : 0 < m)
    (hdivisible : ∀ term ∈ f.toList, (m : Int) ∣ term.2) :
    Polynomial.C (m : ZMod (m ^ 2)) *
        toPolyMod (m ^ 2)
          (Generated.StrictHensel.divideThenReduceCoeffs f (m : Int)) =
      toPolyMod (m ^ 2) f := by
  let quotient := Generated.StrictHensel.divideCoeffs f (m : Int)
  let reduced :=
    Generated.StrictHensel.divideThenReduceCoeffs f (m : Int)
  have hcertificate := divideThenReduceCoeffs_certificate f m hdivisible
  have hexact := congrArg
    (Polynomial.map (Int.castRingHom (ZMod (m ^ 2)))) hcertificate.2
  have hexact' :
      Polynomial.C (m : ZMod (m ^ 2)) * toPolyMod (m ^ 2) quotient =
        toPolyMod (m ^ 2) f := by
    simpa [toPolyMod, quotient, Polynomial.map_mul, Polynomial.map_C] using
      hexact
  let π := ZMod.castHom (dvd_pow_self m (by omega : 2 ≠ 0)) (ZMod m)
  have hproject : Polynomial.map π
      (toPolyMod (m ^ 2) reduced - toPolyMod (m ^ 2) quotient) = 0 := by
    rw [Polynomial.map_sub,
      map_toPolyMod_of_dvd (dvd_pow_self m (by omega : 2 ≠ 0)),
      map_toPolyMod_of_dvd (dvd_pow_self m (by omega : 2 ≠ 0))]
    exact sub_eq_zero.mpr hcertificate.1
  have hdiscarded : Polynomial.C (m : ZMod (m ^ 2)) *
      (toPolyMod (m ^ 2) reduced - toPolyMod (m ^ 2) quotient) = 0 := by
    exact C_scale_mul_eq_zero_of_map_eq_zero m hm _ hproject
  calc
    Polynomial.C (m : ZMod (m ^ 2)) * toPolyMod (m ^ 2) reduced =
        Polynomial.C (m : ZMod (m ^ 2)) * toPolyMod (m ^ 2) quotient +
          Polynomial.C (m : ZMod (m ^ 2)) *
            (toPolyMod (m ^ 2) reduced - toPolyMod (m ^ 2) quotient) := by
              ring
    _ = Polynomial.C (m : ZMod (m ^ 2)) *
        toPolyMod (m ^ 2) quotient := by rw [hdiscarded, add_zero]
    _ = toPolyMod (m ^ 2) f := hexact'

/-- Concrete error equation obtained by composing the generated heap
multiplication, generated sparse subtraction, and generated coefficient
division loops. -/
theorem factorError_from_raw_runs
    (termination : Generated.StrictHensel.DivmodTermination)
    (node : HenselNode) (f gh difference : SparsePolyZZ)
    (m : Nat) (hm : 0 < m)
    (hgh : (strictHenselRawOps termination).mul node.g node.h = .ok gh)
    (hdifference :
      (strictHenselRawOps termination).sub f gh = .ok difference)
    (hdivisible : ∀ term ∈ difference.toList, (m : Int) ∣ term.2) :
    toPolyMod (m ^ 2) f =
      toPolyMod (m ^ 2) node.g * toPolyMod (m ^ 2) node.h +
        Polynomial.C (m : ZMod (m ^ 2)) *
          toPolyMod (m ^ 2)
            (Generated.StrictHensel.divideThenReduceCoeffs
              difference (m : Int)) := by
  have hghSemantic := strictHenselRawOps_mul_refines_of_run
    termination (m ^ 2) node.g node.h gh hgh
  have hdifferenceSemantic := strictHenselRawOps_sub_refines_of_run
    termination (m ^ 2) f gh difference hdifference
  have hscaled := divideThenReduceCoeffs_scaled_toPolyMod_sq
    difference m hm hdivisible
  rw [hdifferenceSemantic, hghSemantic] at hscaled
  rw [hscaled]
  ring

/-- Concrete Bezout-error equation obtained from the generated `s*g`, `t*h`,
two sparse subtractions, and the same exact coefficient quotient loop. -/
theorem bezoutError_from_raw_runs
    (termination : Generated.StrictHensel.DivmodTermination)
    (factorNode : HenselNode)
    (sg th oneMinusSg difference : SparsePolyZZ)
    (m : Nat) (hm : 0 < m)
    (hsg : (strictHenselRawOps termination).mul factorNode.s factorNode.g =
      .ok sg)
    (hth : (strictHenselRawOps termination).mul factorNode.t factorNode.h =
      .ok th)
    (honeMinus : (strictHenselRawOps termination).sub
      (#[(UMonomial.mk 0, 1)] : SparsePolyZZ) sg = .ok oneMinusSg)
    (hdifference : (strictHenselRawOps termination).sub oneMinusSg th =
      .ok difference)
    (hdivisible : ∀ term ∈ difference.toList, (m : Int) ∣ term.2) :
    (1 : Polynomial (ZMod (m ^ 2))) =
      toPolyMod (m ^ 2) factorNode.s * toPolyMod (m ^ 2) factorNode.g +
        toPolyMod (m ^ 2) factorNode.t * toPolyMod (m ^ 2) factorNode.h +
        Polynomial.C (m : ZMod (m ^ 2)) *
          toPolyMod (m ^ 2)
            (Generated.StrictHensel.divideThenReduceCoeffs
              difference (m : Int)) := by
  have hsgSemantic := strictHenselRawOps_mul_refines_of_run termination
    (m ^ 2) factorNode.s factorNode.g sg hsg
  have hthSemantic := strictHenselRawOps_mul_refines_of_run termination
    (m ^ 2) factorNode.t factorNode.h th hth
  have honeMinusSemantic := strictHenselRawOps_sub_refines_of_run termination
    (m ^ 2) (#[(UMonomial.mk 0, 1)] : SparsePolyZZ) sg oneMinusSg
    honeMinus
  have hone : toPolyMod (m ^ 2)
      (#[(UMonomial.mk 0, 1)] : SparsePolyZZ) = 1 := by
    simp [toPolyMod, SparsePolyZZ.toPoly]
  rw [hone] at honeMinusSemantic
  have hdifferenceSemantic := strictHenselRawOps_sub_refines_of_run
    termination (m ^ 2) oneMinusSg th difference hdifference
  have hscaled := divideThenReduceCoeffs_scaled_toPolyMod_sq difference m hm
    hdivisible
  rw [hdifferenceSemantic, honeMinusSemantic, hsgSemantic, hthSemantic] at hscaled
  rw [hscaled]
  ring

/-- Full semantic composition of every concrete raw operation in the second
contiguous `__hensel_step` phase. -/
theorem henselBezoutCorrection_from_raw_runs
    (termination : Generated.StrictHensel.DivmodTermination)
    (factorNode : HenselNode) (m : Nat) (hm : 0 < m)
    (sg th oneMinusSg difference sep q r sRaw sNew tep qpg tauRaw tau
      tRaw tNew : SparsePolyZZ)
    (hsg : (strictHenselRawOps termination).mul factorNode.s factorNode.g =
      .ok sg)
    (hth : (strictHenselRawOps termination).mul factorNode.t factorNode.h =
      .ok th)
    (honeMinus : (strictHenselRawOps termination).sub
      (#[(UMonomial.mk 0, 1)] : SparsePolyZZ) sg = .ok oneMinusSg)
    (hdifference : (strictHenselRawOps termination).sub oneMinusSg th =
      .ok difference)
    (hdivisible : ∀ term ∈ difference.toList, (m : Int) ∣ term.2)
    (hsep : (strictHenselRawOps termination).mul factorNode.s
      (Generated.StrictHensel.divideThenReduceCoeffs difference (m : Int)) =
        .ok sep)
    (hh : 0 < factorNode.h.size)
    (hhDegree : factorNode.h[0]!.1.deg < 2 ^ 63)
    (hhHead : HeadDominates factorNode.h) (hsepBound : DegreesBound sep)
    (hinvert : (ZZ.invert 0 factorNode.h[0]!.2 (m : Int)).1 = true)
    (hdivmod : Generated.StrictHensel.__upoly_divmod_mod_raw_ir termination
      sep factorNode.h (m : Int) = .ok (q, r))
    (hsRaw : (strictHenselRawOps termination).add factorNode.s
      (Generated.StrictHensel.scaleCoeffs r (m : Int)) = .ok sRaw)
    (hsNew : Generated.StrictHensel.__upoly_mod_coeff_raw_ir sRaw
      (m ^ 2 : Int) = .ok sNew)
    (htep : (strictHenselRawOps termination).mul factorNode.t
      (Generated.StrictHensel.divideThenReduceCoeffs difference (m : Int)) =
        .ok tep)
    (hqpg : (strictHenselRawOps termination).mul q factorNode.g = .ok qpg)
    (htauRaw : (strictHenselRawOps termination).add tep qpg = .ok tauRaw)
    (htau : Generated.StrictHensel.__upoly_mod_coeff_raw_ir tauRaw
      (m : Int) = .ok tau)
    (htRaw : (strictHenselRawOps termination).add factorNode.t
      (Generated.StrictHensel.scaleCoeffs tau (m : Int)) = .ok tRaw)
    (htNew : Generated.StrictHensel.__upoly_mod_coeff_raw_ir tRaw
      (m ^ 2 : Int) = .ok tNew)
    (hbezoutInput :
      toPolyMod m factorNode.s * toPolyMod m factorNode.g +
        toPolyMod m factorNode.t * toPolyMod m factorNode.h = 1) :
    toPolyMod (m ^ 2) sNew * toPolyMod (m ^ 2) factorNode.g +
        toPolyMod (m ^ 2) tNew * toPolyMod (m ^ 2) factorNode.h = 1 ∧
      toPolyMod m sNew = toPolyMod m factorNode.s ∧
      toPolyMod m tNew = toPolyMod m factorNode.t := by
  let ep := Generated.StrictHensel.divideThenReduceCoeffs
    difference (m : Int)
  have herror := bezoutError_from_raw_runs termination factorNode sg th
    oneMinusSg difference m hm hsg hth honeMinus hdifference hdivisible
  have hsepSemantic := strictHenselRawOps_mul_refines_of_run termination m
    factorNode.s ep sep hsep
  have hdivmodSemantic := __upoly_divmod_mod_raw_ir_refines_of_run termination
    sep factorNode.h q r m hh hhDegree hhHead hsepBound hinvert hdivmod
  rw [hsepSemantic] at hdivmodSemantic
  have htepSemantic := strictHenselRawOps_mul_refines_of_run termination m
    factorNode.t ep tep htep
  have hqpgSemantic := strictHenselRawOps_mul_refines_of_run termination m
    q factorNode.g qpg hqpg
  have htauRawSemantic := strictHenselRawOps_add_refines_of_run termination m
    tep qpg tauRaw htauRaw
  have htauReduce := __upoly_mod_coeff_raw_ir_refines_of_run tauRaw tau m htau
  have htauSemantic : toPolyMod m tau =
      toPolyMod m factorNode.t * toPolyMod m ep +
        toPolyMod m q * toPolyMod m factorNode.g := by
    rw [htauReduce, htauRawSemantic, htepSemantic, hqpgSemantic]
  let π := ZMod.castHom (dvd_pow_self m (by omega : 2 ≠ 0)) (ZMod m)
  have hbezoutMap : Polynomial.map π
      (toPolyMod (m ^ 2) factorNode.s * toPolyMod (m ^ 2) factorNode.g +
        toPolyMod (m ^ 2) factorNode.t * toPolyMod (m ^ 2) factorNode.h) =
      1 := by
    simp only [Polynomial.map_add, Polynomial.map_mul]
    rw [map_toPolyMod_of_dvd (dvd_pow_self m (by omega : 2 ≠ 0)),
      map_toPolyMod_of_dvd (dvd_pow_self m (by omega : 2 ≠ 0)),
      map_toPolyMod_of_dvd (dvd_pow_self m (by omega : 2 ≠ 0)),
      map_toPolyMod_of_dvd (dvd_pow_self m (by omega : 2 ≠ 0))]
    exact hbezoutInput
  have hdivmodMap : Polynomial.map π
      (toPolyMod (m ^ 2) r +
        toPolyMod (m ^ 2) q * toPolyMod (m ^ 2) factorNode.h) =
      Polynomial.map π
        (toPolyMod (m ^ 2) factorNode.s * toPolyMod (m ^ 2) ep) := by
    simp only [Polynomial.map_add, Polynomial.map_mul]
    rw [map_toPolyMod_of_dvd (dvd_pow_self m (by omega : 2 ≠ 0)),
      map_toPolyMod_of_dvd (dvd_pow_self m (by omega : 2 ≠ 0)),
      map_toPolyMod_of_dvd (dvd_pow_self m (by omega : 2 ≠ 0)),
      map_toPolyMod_of_dvd (dvd_pow_self m (by omega : 2 ≠ 0)),
      map_toPolyMod_of_dvd (dvd_pow_self m (by omega : 2 ≠ 0))]
    exact hdivmodSemantic
  have htauMap : Polynomial.map π (toPolyMod (m ^ 2) tau) =
      Polynomial.map π
        (toPolyMod (m ^ 2) factorNode.t * toPolyMod (m ^ 2) ep +
          toPolyMod (m ^ 2) q * toPolyMod (m ^ 2) factorNode.g) := by
    simp only [Polynomial.map_add, Polynomial.map_mul]
    rw [map_toPolyMod_of_dvd (dvd_pow_self m (by omega : 2 ≠ 0)),
      map_toPolyMod_of_dvd (dvd_pow_self m (by omega : 2 ≠ 0)),
      map_toPolyMod_of_dvd (dvd_pow_self m (by omega : 2 ≠ 0)),
      map_toPolyMod_of_dvd (dvd_pow_self m (by omega : 2 ≠ 0)),
      map_toPolyMod_of_dvd (dvd_pow_self m (by omega : 2 ≠ 0))]
    exact htauSemantic
  have halgebra := henselBezoutCorrection_algebra m hm
    (toPolyMod (m ^ 2) factorNode.g) (toPolyMod (m ^ 2) factorNode.h)
    (toPolyMod (m ^ 2) factorNode.s) (toPolyMod (m ^ 2) factorNode.t)
    (toPolyMod (m ^ 2) ep) (toPolyMod (m ^ 2) q)
    (toPolyMod (m ^ 2) r) (toPolyMod (m ^ 2) tau) herror hbezoutMap
    hdivmodMap htauMap
  have hsRawSemantic := strictHenselRawOps_add_refines_of_run termination
    (m ^ 2) factorNode.s
      (Generated.StrictHensel.scaleCoeffs r (m : Int)) sRaw hsRaw
  rw [scaleCoeffs_toPolyMod] at hsRawSemantic
  have hsNewSemantic := __upoly_mod_coeff_raw_ir_refines_of_run sRaw sNew
    (m ^ 2) hsNew
  have htRawSemantic := strictHenselRawOps_add_refines_of_run termination
    (m ^ 2) factorNode.t
      (Generated.StrictHensel.scaleCoeffs tau (m : Int)) tRaw htRaw
  rw [scaleCoeffs_toPolyMod] at htRawSemantic
  have htNewSemantic := __upoly_mod_coeff_raw_ir_refines_of_run tRaw tNew
    (m ^ 2) htNew
  constructor
  · rw [hsNewSemantic, htNewSemantic, hsRawSemantic, htRawSemantic]
    simpa using halgebra
  constructor
  · calc
      toPolyMod m sNew = toPolyMod m sRaw :=
        __upoly_mod_coeff_raw_ir_preserves_divisor_of_run
          (dvd_pow_self m (by omega : 2 ≠ 0)) sRaw sNew hsNew
      _ = toPolyMod m factorNode.s + toPolyMod m
          (Generated.StrictHensel.scaleCoeffs r (m : Int)) :=
        strictHenselRawOps_add_refines_of_run termination m factorNode.s _
          sRaw hsRaw
      _ = toPolyMod m factorNode.s := by rw [scaleCoeffs_toPolyMod]; simp
  · calc
      toPolyMod m tNew = toPolyMod m tRaw :=
        __upoly_mod_coeff_raw_ir_preserves_divisor_of_run
          (dvd_pow_self m (by omega : 2 ≠ 0)) tRaw tNew htNew
      _ = toPolyMod m factorNode.t + toPolyMod m
          (Generated.StrictHensel.scaleCoeffs tau (m : Int)) :=
        strictHenselRawOps_add_refines_of_run termination m factorNode.t _
          tRaw htRaw
      _ = toPolyMod m factorNode.t := by rw [scaleCoeffs_toPolyMod]; simp

/-- Full semantic composition of every concrete raw operation in the first
contiguous `__hensel_step` phase.  The named intermediates are not witnesses
chosen by the proof: every one is fixed by an accompanying generated raw-run
equation. -/
theorem henselFactorCorrection_from_raw_runs
    (termination : Generated.StrictHensel.DivmodTermination)
    (node : HenselNode) (f : SparsePolyZZ) (m : Nat) (hm : 0 < m)
    (gh difference se q r te qg tauRaw tau gRaw gNew hRaw hNew :
      SparsePolyZZ)
    (hgh : (strictHenselRawOps termination).mul node.g node.h = .ok gh)
    (hdifference :
      (strictHenselRawOps termination).sub f gh = .ok difference)
    (hdivisible : ∀ term ∈ difference.toList, (m : Int) ∣ term.2)
    (hse : (strictHenselRawOps termination).mul node.s
      (Generated.StrictHensel.divideThenReduceCoeffs difference (m : Int)) =
        .ok se)
    (hh : 0 < node.h.size) (hhDegree : node.h[0]!.1.deg < 2 ^ 63)
    (hhHead : HeadDominates node.h) (hseBound : DegreesBound se)
    (hinvert : (ZZ.invert 0 node.h[0]!.2 (m : Int)).1 = true)
    (hdivmod : Generated.StrictHensel.__upoly_divmod_mod_raw_ir termination
      se node.h (m : Int) = .ok (q, r))
    (hte : (strictHenselRawOps termination).mul node.t
      (Generated.StrictHensel.divideThenReduceCoeffs difference (m : Int)) =
        .ok te)
    (hqg : (strictHenselRawOps termination).mul q node.g = .ok qg)
    (htauRaw : (strictHenselRawOps termination).add te qg = .ok tauRaw)
    (htau : Generated.StrictHensel.__upoly_mod_coeff_raw_ir tauRaw
      (m : Int) = .ok tau)
    (hgRaw : (strictHenselRawOps termination).add node.g
      (Generated.StrictHensel.scaleCoeffs tau (m : Int)) = .ok gRaw)
    (hgNew : Generated.StrictHensel.__upoly_mod_coeff_raw_ir gRaw
      (m ^ 2 : Int) = .ok gNew)
    (hhRaw : (strictHenselRawOps termination).add node.h
      (Generated.StrictHensel.scaleCoeffs r (m : Int)) = .ok hRaw)
    (hhNew : Generated.StrictHensel.__upoly_mod_coeff_raw_ir hRaw
      (m ^ 2 : Int) = .ok hNew)
    (hinvariant : HenselNodeInvariant f m node) :
    toPolyMod (m ^ 2) gNew * toPolyMod (m ^ 2) hNew =
        toPolyMod (m ^ 2) f ∧
      toPolyMod m gNew = toPolyMod m node.g ∧
      toPolyMod m hNew = toPolyMod m node.h := by
  let e := Generated.StrictHensel.divideThenReduceCoeffs
    difference (m : Int)
  have herror := factorError_from_raw_runs termination node f gh difference
    m hm hgh hdifference hdivisible
  have hseSemantic := strictHenselRawOps_mul_refines_of_run
    termination m node.s e se hse
  have hdivmodSemantic := __upoly_divmod_mod_raw_ir_refines_of_run
    termination se node.h q r m hh hhDegree hhHead hseBound hinvert hdivmod
  have hteSemantic := strictHenselRawOps_mul_refines_of_run
    termination m node.t e te hte
  have hqgSemantic := strictHenselRawOps_mul_refines_of_run
    termination m q node.g qg hqg
  have htauRawSemantic := strictHenselRawOps_add_refines_of_run
    termination m te qg tauRaw htauRaw
  have htauReduce := __upoly_mod_coeff_raw_ir_refines_of_run
    tauRaw tau m htau
  have htauSemantic : toPolyMod m tau =
      toPolyMod m node.t * toPolyMod m e +
        toPolyMod m q * toPolyMod m node.g := by
    rw [htauReduce, htauRawSemantic, hteSemantic, hqgSemantic]
  let π := ZMod.castHom (dvd_pow_self m (by omega : 2 ≠ 0)) (ZMod m)
  have hbezoutMap : Polynomial.map π
      (toPolyMod (m ^ 2) node.s * toPolyMod (m ^ 2) node.g +
        toPolyMod (m ^ 2) node.t * toPolyMod (m ^ 2) node.h) = 1 := by
    simp only [Polynomial.map_add, Polynomial.map_mul]
    rw [map_toPolyMod_of_dvd (dvd_pow_self m (by omega : 2 ≠ 0)),
      map_toPolyMod_of_dvd (dvd_pow_self m (by omega : 2 ≠ 0)),
      map_toPolyMod_of_dvd (dvd_pow_self m (by omega : 2 ≠ 0)),
      map_toPolyMod_of_dvd (dvd_pow_self m (by omega : 2 ≠ 0))]
    exact hinvariant.2
  have hdivmodMap : Polynomial.map π
      (toPolyMod (m ^ 2) r +
        toPolyMod (m ^ 2) q * toPolyMod (m ^ 2) node.h) =
      Polynomial.map π
        (toPolyMod (m ^ 2) node.s * toPolyMod (m ^ 2) e) := by
    simp only [Polynomial.map_add, Polynomial.map_mul]
    rw [map_toPolyMod_of_dvd (dvd_pow_self m (by omega : 2 ≠ 0)),
      map_toPolyMod_of_dvd (dvd_pow_self m (by omega : 2 ≠ 0)),
      map_toPolyMod_of_dvd (dvd_pow_self m (by omega : 2 ≠ 0)),
      map_toPolyMod_of_dvd (dvd_pow_self m (by omega : 2 ≠ 0)),
      map_toPolyMod_of_dvd (dvd_pow_self m (by omega : 2 ≠ 0))]
    rw [hseSemantic] at hdivmodSemantic
    exact hdivmodSemantic
  have htauMap : Polynomial.map π (toPolyMod (m ^ 2) tau) =
      Polynomial.map π
        (toPolyMod (m ^ 2) node.t * toPolyMod (m ^ 2) e +
          toPolyMod (m ^ 2) q * toPolyMod (m ^ 2) node.g) := by
    simp only [Polynomial.map_add, Polynomial.map_mul]
    rw [map_toPolyMod_of_dvd (dvd_pow_self m (by omega : 2 ≠ 0)),
      map_toPolyMod_of_dvd (dvd_pow_self m (by omega : 2 ≠ 0)),
      map_toPolyMod_of_dvd (dvd_pow_self m (by omega : 2 ≠ 0)),
      map_toPolyMod_of_dvd (dvd_pow_self m (by omega : 2 ≠ 0)),
      map_toPolyMod_of_dvd (dvd_pow_self m (by omega : 2 ≠ 0))]
    exact htauSemantic
  have halgebra := henselFactorCorrection_algebra m hm
    (toPolyMod (m ^ 2) f) (toPolyMod (m ^ 2) node.g)
    (toPolyMod (m ^ 2) node.h) (toPolyMod (m ^ 2) node.s)
    (toPolyMod (m ^ 2) node.t) (toPolyMod (m ^ 2) e)
    (toPolyMod (m ^ 2) q) (toPolyMod (m ^ 2) r)
    (toPolyMod (m ^ 2) tau) herror hbezoutMap hdivmodMap htauMap
  have hgRawSemantic := strictHenselRawOps_add_refines_of_run
    termination (m ^ 2) node.g
      (Generated.StrictHensel.scaleCoeffs tau (m : Int)) gRaw hgRaw
  rw [scaleCoeffs_toPolyMod] at hgRawSemantic
  have hgNewSemantic := __upoly_mod_coeff_raw_ir_refines_of_run
    gRaw gNew (m ^ 2) hgNew
  have hhRawSemantic := strictHenselRawOps_add_refines_of_run
    termination (m ^ 2) node.h
      (Generated.StrictHensel.scaleCoeffs r (m : Int)) hRaw hhRaw
  rw [scaleCoeffs_toPolyMod] at hhRawSemantic
  have hhNewSemantic := __upoly_mod_coeff_raw_ir_refines_of_run
    hRaw hNew (m ^ 2) hhNew
  constructor
  · rw [hgNewSemantic, hhNewSemantic, hgRawSemantic, hhRawSemantic]
    simpa using halgebra
  constructor
  · calc
      toPolyMod m gNew = toPolyMod m gRaw :=
        __upoly_mod_coeff_raw_ir_preserves_divisor_of_run
          (dvd_pow_self m (by omega : 2 ≠ 0)) gRaw gNew hgNew
      _ = toPolyMod m node.g + toPolyMod m
          (Generated.StrictHensel.scaleCoeffs tau (m : Int)) :=
        strictHenselRawOps_add_refines_of_run termination m node.g _ gRaw hgRaw
      _ = toPolyMod m node.g := by rw [scaleCoeffs_toPolyMod]; simp
  · calc
      toPolyMod m hNew = toPolyMod m hRaw :=
        __upoly_mod_coeff_raw_ir_preserves_divisor_of_run
          (dvd_pow_self m (by omega : 2 ≠ 0)) hRaw hNew hhNew
      _ = toPolyMod m node.h + toPolyMod m
          (Generated.StrictHensel.scaleCoeffs r (m : Int)) :=
        strictHenselRawOps_add_refines_of_run termination m node.h _ hRaw hhRaw
      _ = toPolyMod m node.h := by rw [scaleCoeffs_toPolyMod]; simp

set_option maxHeartbeats 0 in
/-- The complete generated factor phase refines its L2 factor invariant.
All local values below are definitionally the outputs of the concrete
generated operations, including the result of the supplied well-founded
division trace. -/
theorem __hensel_step_factor_phase_raw_ir_refines
    (termination : Generated.StrictHensel.DivmodTermination)
    (node : HenselNode) (f : SparsePolyZZ) (m : Nat) (hm : 0 < m)
    (hh : 0 < node.h.size) (hhDegree : node.h[0]!.1.deg < 2 ^ 63)
    (hhHead : HeadDominates node.h)
    (hinvert : (ZZ.invert 0 node.h[0]!.2 (m : Int)).1 = true)
    (hdivisible :
      let gh := Generated.StrictHensel.pairVecMulHeapLoop
        (Generated.StrictHensel.pairVecMulProducts node.g node.h) #[]
      let difference := Generated.StrictHensel.pairVecSubLoop f gh 0 0 #[]
      ∀ term ∈ difference.toList, (m : Int) ∣ term.2)
    (hseBound :
      let gh := Generated.StrictHensel.pairVecMulHeapLoop
        (Generated.StrictHensel.pairVecMulProducts node.g node.h) #[]
      let difference := Generated.StrictHensel.pairVecSubLoop f gh 0 0 #[]
      let e := Generated.StrictHensel.divideThenReduceCoeffs
        difference (m : Int)
      let se := Generated.StrictHensel.pairVecMulHeapLoop
        (Generated.StrictHensel.pairVecMulProducts node.s e) #[]
      DegreesBound se)
    (hinvariant : HenselNodeInvariant f m node) :
    ∃ factorNode,
      Generated.StrictHensel.__hensel_step_factor_phase_raw_ir
          (strictHenselRawOps termination) node f (m : Int) = .ok factorNode ∧
      toPolyMod (m ^ 2) factorNode.g *
          toPolyMod (m ^ 2) factorNode.h = toPolyMod (m ^ 2) f ∧
      toPolyMod m factorNode.g = toPolyMod m node.g ∧
      toPolyMod m factorNode.h = toPolyMod m node.h ∧
      factorNode.s = node.s ∧ factorNode.t = node.t := by
  let gh := Generated.StrictHensel.pairVecMulHeapLoop
    (Generated.StrictHensel.pairVecMulProducts node.g node.h) #[]
  let difference := Generated.StrictHensel.pairVecSubLoop f gh 0 0 #[]
  let e := Generated.StrictHensel.divideThenReduceCoeffs difference (m : Int)
  let se := Generated.StrictHensel.pairVecMulHeapLoop
    (Generated.StrictHensel.pairVecMulProducts node.s e) #[]
  let reduced := Generated.StrictHensel.modCoeffOutput se (m : Int)
  let trace := termination.trace se node.h reduced (m : Int) (by rfl)
  let qr := Generated.StrictHensel.divmodLoop node.h
    (ZZ.invert 0 node.h[0]!.2 (m : Int)).2 (m : Int) trace
  let te := Generated.StrictHensel.pairVecMulHeapLoop
    (Generated.StrictHensel.pairVecMulProducts node.t e) #[]
  let qg := Generated.StrictHensel.pairVecMulHeapLoop
    (Generated.StrictHensel.pairVecMulProducts qr.1 node.g) #[]
  let tauRaw := Generated.StrictHensel.pairVecAddLoop te qg 0 0 #[]
  let tau := Generated.StrictHensel.modCoeffOutput tauRaw (m : Int)
  let gRaw := Generated.StrictHensel.pairVecAddLoop node.g
    (Generated.StrictHensel.scaleCoeffs tau (m : Int)) 0 0 #[]
  let gNew := Generated.StrictHensel.modCoeffOutput gRaw (m ^ 2 : Int)
  let hRaw := Generated.StrictHensel.pairVecAddLoop node.h
    (Generated.StrictHensel.scaleCoeffs qr.2 (m : Int)) 0 0 #[]
  let hNew := Generated.StrictHensel.modCoeffOutput hRaw (m ^ 2 : Int)
  let factorNode : HenselNode := { node with g := gNew, h := hNew }
  have hhFalse : node.h.isEmpty = false := by
    simpa [Array.isEmpty_iff] using (show node.h ≠ #[] by
      intro hempty
      have hsize : node.h.size = 0 := by rw [hempty]; rfl
      exact (Nat.ne_of_gt hh) hsize)
  have hghRun : (strictHenselRawOps termination).mul node.g node.h =
      .ok gh := rfl
  have hdifferenceRun : (strictHenselRawOps termination).sub f gh =
      .ok difference := rfl
  have hseRun : (strictHenselRawOps termination).mul node.s e = .ok se := rfl
  have hdivmodRun :
      Generated.StrictHensel.__upoly_divmod_mod_raw_ir termination
        se node.h (m : Int) = .ok qr := by
    simp [Generated.StrictHensel.__upoly_divmod_mod_raw_ir, hhFalse,
      hinvert, Generated.StrictHensel.__upoly_mod_coeff_raw_ir,
      reduced, trace, qr]
  have hteRun : (strictHenselRawOps termination).mul node.t e = .ok te := rfl
  have hqgRun : (strictHenselRawOps termination).mul qr.1 node.g = .ok qg := rfl
  have htauRawRun : (strictHenselRawOps termination).add te qg =
      .ok tauRaw := rfl
  have htauRun : Generated.StrictHensel.__upoly_mod_coeff_raw_ir tauRaw
      (m : Int) = .ok tau := rfl
  have hgRawRun : (strictHenselRawOps termination).add node.g
      (Generated.StrictHensel.scaleCoeffs tau (m : Int)) = .ok gRaw := rfl
  have hgNewRun : Generated.StrictHensel.__upoly_mod_coeff_raw_ir gRaw
      (m ^ 2 : Int) = .ok gNew := rfl
  have hhRawRun : (strictHenselRawOps termination).add node.h
      (Generated.StrictHensel.scaleCoeffs qr.2 (m : Int)) = .ok hRaw := rfl
  have hhNewRun : Generated.StrictHensel.__upoly_mod_coeff_raw_ir hRaw
      (m ^ 2 : Int) = .ok hNew := rfl
  have hsemantic := henselFactorCorrection_from_raw_runs termination node f m
    hm gh difference se qr.1 qr.2 te qg tauRaw tau gRaw gNew hRaw hNew
    hghRun hdifferenceRun (by simpa [gh, difference] using hdivisible)
    hseRun hh hhDegree hhHead (by
      simpa [gh, difference, e, se] using hseBound) hinvert hdivmodRun
    hteRun hqgRun htauRawRun htauRun hgRawRun hgNewRun hhRawRun hhNewRun
    hinvariant
  refine ⟨factorNode, ?_, ?_⟩
  · have hm2 : (m : Int) * (m : Int) = (m ^ 2 : Int) := by
      norm_num [pow_two]
    rw [Generated.StrictHensel.__hensel_step_factor_phase_raw_ir]
    dsimp only [
      strictHenselRawOps, Generated.StrictHensel.__upoly_mul_raw_ir,
      Generated.StrictHensel.__upoly_add_raw_ir,
      Generated.StrictHensel.__upoly_sub_raw_ir]
    simp only [Bind.bind, Except.bind]
    rw [show Generated.StrictHensel.__upoly_divmod_mod_raw_ir termination
        (Generated.StrictHensel.pairVecMulHeapLoop
          (Generated.StrictHensel.pairVecMulProducts node.s
            (Generated.StrictHensel.divideThenReduceCoeffs
              (Generated.StrictHensel.pairVecSubLoop f
                (Generated.StrictHensel.pairVecMulHeapLoop
                  (Generated.StrictHensel.pairVecMulProducts node.g node.h)
                  #[]) 0 0 #[]) (m : Int))) #[])
        node.h (m : Int) = .ok qr by simpa [gh, difference, e, se] using
          hdivmodRun]
    simp only [Except.bind]
    rw [show Generated.StrictHensel.__upoly_mod_coeff_raw_ir
        (Generated.StrictHensel.pairVecAddLoop te qg 0 0 #[]) (m : Int) =
          .ok tau by simpa [tauRaw] using htauRun]
    simp only [Except.bind]
    rw [hm2]
    rw [show Generated.StrictHensel.__upoly_mod_coeff_raw_ir
        (Generated.StrictHensel.pairVecAddLoop node.g
          (Generated.StrictHensel.scaleCoeffs tau (m : Int)) 0 0 #[])
        (m ^ 2 : Int) = .ok gNew by simpa [gRaw] using hgNewRun]
    simp only [Except.bind]
    rw [show Generated.StrictHensel.__upoly_mod_coeff_raw_ir
        (Generated.StrictHensel.pairVecAddLoop node.h
          (Generated.StrictHensel.scaleCoeffs qr.2 (m : Int)) 0 0 #[])
        (m ^ 2 : Int) = .ok hNew by simpa [hRaw] using hhNewRun]
    rfl
  · exact ⟨hsemantic.1, hsemantic.2.1, hsemantic.2.2, rfl, rfl⟩

set_option maxHeartbeats 0 in
/-- The complete generated Bezout phase refines its L2 certificate invariant.
Every local is definitionally produced by the generated raw operations and
the supplied finite well-founded division trace. -/
theorem __hensel_step_bezout_phase_raw_ir_refines
    (termination : Generated.StrictHensel.DivmodTermination)
    (factorNode : HenselNode) (m : Nat) (hm : 0 < m)
    (hh : 0 < factorNode.h.size)
    (hhDegree : factorNode.h[0]!.1.deg < 2 ^ 63)
    (hhHead : HeadDominates factorNode.h)
    (hinvert : (ZZ.invert 0 factorNode.h[0]!.2 (m : Int)).1 = true)
    (hdivisible :
      let sg := Generated.StrictHensel.pairVecMulHeapLoop
        (Generated.StrictHensel.pairVecMulProducts factorNode.s factorNode.g)
        #[]
      let th := Generated.StrictHensel.pairVecMulHeapLoop
        (Generated.StrictHensel.pairVecMulProducts factorNode.t factorNode.h)
        #[]
      let oneMinusSg := Generated.StrictHensel.pairVecSubLoop
        (#[(UMonomial.mk 0, 1)] : SparsePolyZZ) sg 0 0 #[]
      let difference := Generated.StrictHensel.pairVecSubLoop
        oneMinusSg th 0 0 #[]
      ∀ term ∈ difference.toList, (m : Int) ∣ term.2)
    (hsepBound :
      let sg := Generated.StrictHensel.pairVecMulHeapLoop
        (Generated.StrictHensel.pairVecMulProducts factorNode.s factorNode.g)
        #[]
      let th := Generated.StrictHensel.pairVecMulHeapLoop
        (Generated.StrictHensel.pairVecMulProducts factorNode.t factorNode.h)
        #[]
      let oneMinusSg := Generated.StrictHensel.pairVecSubLoop
        (#[(UMonomial.mk 0, 1)] : SparsePolyZZ) sg 0 0 #[]
      let difference := Generated.StrictHensel.pairVecSubLoop
        oneMinusSg th 0 0 #[]
      let ep := Generated.StrictHensel.divideThenReduceCoeffs
        difference (m : Int)
      let sep := Generated.StrictHensel.pairVecMulHeapLoop
        (Generated.StrictHensel.pairVecMulProducts factorNode.s ep) #[]
      DegreesBound sep)
    (hbezoutInput :
      toPolyMod m factorNode.s * toPolyMod m factorNode.g +
        toPolyMod m factorNode.t * toPolyMod m factorNode.h = 1) :
    ∃ output,
      Generated.StrictHensel.__hensel_step_bezout_phase_raw_ir
          (strictHenselRawOps termination) factorNode (m : Int) = .ok output ∧
      toPolyMod (m ^ 2) output.s * toPolyMod (m ^ 2) output.g +
          toPolyMod (m ^ 2) output.t * toPolyMod (m ^ 2) output.h = 1 ∧
      toPolyMod m output.s = toPolyMod m factorNode.s ∧
      toPolyMod m output.t = toPolyMod m factorNode.t ∧
      output.g = factorNode.g ∧ output.h = factorNode.h := by
  let sg := Generated.StrictHensel.pairVecMulHeapLoop
    (Generated.StrictHensel.pairVecMulProducts factorNode.s factorNode.g) #[]
  let th := Generated.StrictHensel.pairVecMulHeapLoop
    (Generated.StrictHensel.pairVecMulProducts factorNode.t factorNode.h) #[]
  let one : SparsePolyZZ := #[(UMonomial.mk 0, 1)]
  let oneMinusSg := Generated.StrictHensel.pairVecSubLoop one sg 0 0 #[]
  let difference := Generated.StrictHensel.pairVecSubLoop oneMinusSg th 0 0 #[]
  let ep := Generated.StrictHensel.divideThenReduceCoeffs difference (m : Int)
  let sep := Generated.StrictHensel.pairVecMulHeapLoop
    (Generated.StrictHensel.pairVecMulProducts factorNode.s ep) #[]
  let reduced := Generated.StrictHensel.modCoeffOutput sep (m : Int)
  let trace := termination.trace sep factorNode.h reduced (m : Int) (by rfl)
  let qr := Generated.StrictHensel.divmodLoop factorNode.h
    (ZZ.invert 0 factorNode.h[0]!.2 (m : Int)).2 (m : Int) trace
  let sRaw := Generated.StrictHensel.pairVecAddLoop factorNode.s
    (Generated.StrictHensel.scaleCoeffs qr.2 (m : Int)) 0 0 #[]
  let sNew := Generated.StrictHensel.modCoeffOutput sRaw (m ^ 2 : Int)
  let tep := Generated.StrictHensel.pairVecMulHeapLoop
    (Generated.StrictHensel.pairVecMulProducts factorNode.t ep) #[]
  let qpg := Generated.StrictHensel.pairVecMulHeapLoop
    (Generated.StrictHensel.pairVecMulProducts qr.1 factorNode.g) #[]
  let tauRaw := Generated.StrictHensel.pairVecAddLoop tep qpg 0 0 #[]
  let tau := Generated.StrictHensel.modCoeffOutput tauRaw (m : Int)
  let tRaw := Generated.StrictHensel.pairVecAddLoop factorNode.t
    (Generated.StrictHensel.scaleCoeffs tau (m : Int)) 0 0 #[]
  let tNew := Generated.StrictHensel.modCoeffOutput tRaw (m ^ 2 : Int)
  let output : HenselNode := { factorNode with s := sNew, t := tNew }
  have hhFalse : factorNode.h.isEmpty = false := by
    simpa [Array.isEmpty_iff] using (show factorNode.h ≠ #[] by
      intro hempty
      have hsize : factorNode.h.size = 0 := by rw [hempty]; rfl
      exact (Nat.ne_of_gt hh) hsize)
  have hsgRun : (strictHenselRawOps termination).mul factorNode.s
      factorNode.g = .ok sg := rfl
  have hthRun : (strictHenselRawOps termination).mul factorNode.t
      factorNode.h = .ok th := rfl
  have honeMinusRun : (strictHenselRawOps termination).sub one sg =
      .ok oneMinusSg := rfl
  have hdifferenceRun : (strictHenselRawOps termination).sub oneMinusSg th =
      .ok difference := rfl
  have hsepRun : (strictHenselRawOps termination).mul factorNode.s ep =
      .ok sep := rfl
  have hdivmodRun :
      Generated.StrictHensel.__upoly_divmod_mod_raw_ir termination sep
        factorNode.h (m : Int) = .ok qr := by
    simp [Generated.StrictHensel.__upoly_divmod_mod_raw_ir, hhFalse,
      hinvert, Generated.StrictHensel.__upoly_mod_coeff_raw_ir,
      reduced, trace, qr]
  have hsRawRun : (strictHenselRawOps termination).add factorNode.s
      (Generated.StrictHensel.scaleCoeffs qr.2 (m : Int)) = .ok sRaw := rfl
  have hsNewRun : Generated.StrictHensel.__upoly_mod_coeff_raw_ir sRaw
      (m ^ 2 : Int) = .ok sNew := rfl
  have htepRun : (strictHenselRawOps termination).mul factorNode.t ep =
      .ok tep := rfl
  have hqpgRun : (strictHenselRawOps termination).mul qr.1 factorNode.g =
      .ok qpg := rfl
  have htauRawRun : (strictHenselRawOps termination).add tep qpg =
      .ok tauRaw := rfl
  have htauRun : Generated.StrictHensel.__upoly_mod_coeff_raw_ir tauRaw
      (m : Int) = .ok tau := rfl
  have htRawRun : (strictHenselRawOps termination).add factorNode.t
      (Generated.StrictHensel.scaleCoeffs tau (m : Int)) = .ok tRaw := rfl
  have htNewRun : Generated.StrictHensel.__upoly_mod_coeff_raw_ir tRaw
      (m ^ 2 : Int) = .ok tNew := rfl
  have hsemantic := henselBezoutCorrection_from_raw_runs termination factorNode
    m hm sg th oneMinusSg difference sep qr.1 qr.2 sRaw sNew tep qpg
    tauRaw tau tRaw tNew hsgRun hthRun (by simpa [one] using honeMinusRun)
    hdifferenceRun (by
      simpa [sg, th, one, oneMinusSg, difference] using hdivisible)
    hsepRun hh hhDegree hhHead (by
      simpa [sg, th, one, oneMinusSg, difference, ep, sep] using hsepBound)
    hinvert hdivmodRun hsRawRun hsNewRun htepRun hqpgRun htauRawRun htauRun
    htRawRun htNewRun hbezoutInput
  refine ⟨output, ?_, ?_⟩
  · have hm2 : (m : Int) * (m : Int) = (m ^ 2 : Int) := by
      norm_num [pow_two]
    rw [Generated.StrictHensel.__hensel_step_bezout_phase_raw_ir]
    dsimp only [strictHenselRawOps, Generated.StrictHensel.__upoly_mul_raw_ir,
      Generated.StrictHensel.__upoly_add_raw_ir,
      Generated.StrictHensel.__upoly_sub_raw_ir]
    simp only [Bind.bind, Except.bind]
    rw [show Generated.StrictHensel.__upoly_divmod_mod_raw_ir termination
        (Generated.StrictHensel.pairVecMulHeapLoop
          (Generated.StrictHensel.pairVecMulProducts factorNode.s
            (Generated.StrictHensel.divideThenReduceCoeffs
              (Generated.StrictHensel.pairVecSubLoop
                (Generated.StrictHensel.pairVecSubLoop one sg 0 0 #[])
                th 0 0 #[]) (m : Int))) #[])
        factorNode.h (m : Int) = .ok qr by
          simpa [oneMinusSg, difference, ep, sep] using hdivmodRun]
    simp only [Except.bind]
    rw [hm2]
    rw [show Generated.StrictHensel.__upoly_mod_coeff_raw_ir
        (Generated.StrictHensel.pairVecAddLoop factorNode.s
          (Generated.StrictHensel.scaleCoeffs qr.2 (m : Int)) 0 0 #[])
        (m ^ 2 : Int) = .ok sNew by simpa [sRaw] using hsNewRun]
    simp only [Except.bind]
    rw [show Generated.StrictHensel.__upoly_mod_coeff_raw_ir
        (Generated.StrictHensel.pairVecAddLoop tep qpg 0 0 #[]) (m : Int) =
          .ok tau by simpa [tauRaw] using htauRun]
    simp only [Except.bind]
    rw [show Generated.StrictHensel.__upoly_mod_coeff_raw_ir
        (Generated.StrictHensel.pairVecAddLoop factorNode.t
          (Generated.StrictHensel.scaleCoeffs tau (m : Int)) 0 0 #[])
        (m ^ 2 : Int) = .ok tNew by simpa [tRaw] using htNewRun]
    rfl
  · exact ⟨hsemantic.1, hsemantic.2.1, hsemantic.2.2, rfl, rfl⟩

/-- Safety and exact-divisibility facts required by the second source phase.
This structure contains no semantic output or correction witness. -/
structure BezoutPhasePreconditions
    (termination : Generated.StrictHensel.DivmodTermination)
    (factorNode : HenselNode) (m : Nat) : Prop where
  hNonempty : 0 < factorNode.h.size
  hDegree : factorNode.h[0]!.1.deg < 2 ^ 63
  hHead : HeadDominates factorNode.h
  hInvertible : (ZZ.invert 0 factorNode.h[0]!.2 (m : Int)).1 = true
  differenceDivisible :
    let sg := Generated.StrictHensel.pairVecMulHeapLoop
      (Generated.StrictHensel.pairVecMulProducts factorNode.s factorNode.g) #[]
    let th := Generated.StrictHensel.pairVecMulHeapLoop
      (Generated.StrictHensel.pairVecMulProducts factorNode.t factorNode.h) #[]
    let oneMinusSg := Generated.StrictHensel.pairVecSubLoop
      (#[(UMonomial.mk 0, 1)] : SparsePolyZZ) sg 0 0 #[]
    let difference := Generated.StrictHensel.pairVecSubLoop oneMinusSg th 0 0 #[]
    ∀ term ∈ difference.toList, (m : Int) ∣ term.2
  sepBound :
    let sg := Generated.StrictHensel.pairVecMulHeapLoop
      (Generated.StrictHensel.pairVecMulProducts factorNode.s factorNode.g) #[]
    let th := Generated.StrictHensel.pairVecMulHeapLoop
      (Generated.StrictHensel.pairVecMulProducts factorNode.t factorNode.h) #[]
    let oneMinusSg := Generated.StrictHensel.pairVecSubLoop
      (#[(UMonomial.mk 0, 1)] : SparsePolyZZ) sg 0 0 #[]
    let difference := Generated.StrictHensel.pairVecSubLoop oneMinusSg th 0 0 #[]
    let ep := Generated.StrictHensel.divideThenReduceCoeffs difference (m : Int)
    let sep := Generated.StrictHensel.pairVecMulHeapLoop
      (Generated.StrictHensel.pairVecMulProducts factorNode.s ep) #[]
    DegreesBound sep

/-- Complete algorithm invariant for the generated quadratic Hensel step.
It records only source safety, representation bounds, exact divisibility, and
the input L2 invariant.  The intermediate factor node is universally tied to
the unique generated first-phase execution. -/
structure HenselStepRefinementInvariant
    (termination : Generated.StrictHensel.DivmodTermination)
    (node : HenselNode) (f : SparsePolyZZ) (m : Nat) : Prop where
  positiveModulus : 0 < m
  inputInvariant : HenselNodeInvariant f m node
  hNonempty : 0 < node.h.size
  hDegree : node.h[0]!.1.deg < 2 ^ 63
  hHead : HeadDominates node.h
  hInvertible : (ZZ.invert 0 node.h[0]!.2 (m : Int)).1 = true
  factorDifferenceDivisible :
    let gh := Generated.StrictHensel.pairVecMulHeapLoop
      (Generated.StrictHensel.pairVecMulProducts node.g node.h) #[]
    let difference := Generated.StrictHensel.pairVecSubLoop f gh 0 0 #[]
    ∀ term ∈ difference.toList, (m : Int) ∣ term.2
  factorSeBound :
    let gh := Generated.StrictHensel.pairVecMulHeapLoop
      (Generated.StrictHensel.pairVecMulProducts node.g node.h) #[]
    let difference := Generated.StrictHensel.pairVecSubLoop f gh 0 0 #[]
    let e := Generated.StrictHensel.divideThenReduceCoeffs difference (m : Int)
    let se := Generated.StrictHensel.pairVecMulHeapLoop
      (Generated.StrictHensel.pairVecMulProducts node.s e) #[]
    DegreesBound se
  bezoutReady : ∀ factorNode,
    Generated.StrictHensel.__hensel_step_factor_phase_raw_ir
        (strictHenselRawOps termination) node f (m : Int) = .ok factorNode →
    BezoutPhasePreconditions termination factorNode m

set_option maxHeartbeats 0 in
/-- Final strict refinement of the generated C++ `__hensel_step` entry.
The output is obtained solely by executing its two generated raw phases. -/
theorem __hensel_step_raw_ir_refines
    (termination : Generated.StrictHensel.DivmodTermination)
    (node : HenselNode) (f : SparsePolyZZ) (m : Nat)
    (hinvariant : HenselStepRefinementInvariant termination node f m) :
    ∃ output,
      Generated.StrictHensel.__hensel_step_raw_ir
          (strictHenselRawOps termination) node f (m : Int) = .ok output ∧
      HenselStepCorrect f m node output := by
  rcases __hensel_step_factor_phase_raw_ir_refines termination node f m
      hinvariant.positiveModulus hinvariant.hNonempty hinvariant.hDegree
      hinvariant.hHead hinvariant.hInvertible
      hinvariant.factorDifferenceDivisible hinvariant.factorSeBound
      hinvariant.inputInvariant with
    ⟨factorNode, hfactorRun, hfactorProduct, hgPreserved, hhPreserved,
      hsUnchanged, htUnchanged⟩
  have hready := hinvariant.bezoutReady factorNode hfactorRun
  have hfactorBezout :
      toPolyMod m factorNode.s * toPolyMod m factorNode.g +
        toPolyMod m factorNode.t * toPolyMod m factorNode.h = 1 := by
    rw [hsUnchanged, htUnchanged, hgPreserved, hhPreserved]
    exact hinvariant.inputInvariant.2
  rcases __hensel_step_bezout_phase_raw_ir_refines termination factorNode m
      hinvariant.positiveModulus hready.hNonempty hready.hDegree hready.hHead
      hready.hInvertible hready.differenceDivisible hready.sepBound
      hfactorBezout with
    ⟨output, hbezoutRun, hbezout, hsPreserved, htPreserved,
      hgUnchanged, hhUnchanged⟩
  refine ⟨output, ?_, ?_⟩
  · change (Generated.StrictHensel.__hensel_step_factor_phase_raw_ir
        (strictHenselRawOps termination) node f (m : Int) >>= fun factorNode =>
      Generated.StrictHensel.__hensel_step_bezout_phase_raw_ir
        (strictHenselRawOps termination) factorNode (m : Int)) = .ok output
    rw [hfactorRun]
    exact hbezoutRun
  · constructor
    · constructor
      · rw [hgUnchanged, hhUnchanged]
        exact hfactorProduct
      · exact hbezout
    constructor
    · rw [hgUnchanged]
      exact hgPreserved
    constructor
    · rw [hhUnchanged]
      exact hhPreserved
    constructor
    · calc
        toPolyMod m output.s = toPolyMod m factorNode.s := hsPreserved
        _ = toPolyMod m node.s := by rw [hsUnchanged]
    · calc
        toPolyMod m output.t = toPolyMod m factorNode.t := htPreserved
        _ = toPolyMod m node.t := by rw [htUnchanged]

/-- The first contiguous source phase has a genuine raw-to-safe execution
bridge.  Its only possible source assertion is the modular division by `h`;
all generated heap arithmetic and coefficient-reduction calls are total. -/
theorem __hensel_step_factor_phase_raw_ir_terminates
    (termination : Generated.StrictHensel.DivmodTermination)
    (node : HenselNode) (f : SparsePolyZZ) (m : ZZ)
    (hh : node.h.isEmpty = false)
    (hinvert : (ZZ.invert 0 node.h[0]!.2 m).1 = true) :
    ∃ output,
      Generated.StrictHensel.__hensel_step_factor_phase_raw_ir
          (strictHenselRawOps termination) node f m = .ok output := by
  simp [Generated.StrictHensel.__hensel_step_factor_phase_raw_ir,
    strictHenselRawOps, Generated.StrictHensel.__upoly_mul_raw_ir,
    Generated.StrictHensel.__upoly_add_raw_ir,
    Generated.StrictHensel.__upoly_sub_raw_ir,
    Generated.StrictHensel.__upoly_divmod_mod_raw_ir, hh, hinvert,
    Generated.StrictHensel.__upoly_mod_coeff_raw_ir]
  exact ⟨_, rfl⟩

/-- Matching raw-to-safe bridge for the second contiguous source phase. -/
theorem __hensel_step_bezout_phase_raw_ir_terminates
    (termination : Generated.StrictHensel.DivmodTermination)
    (factorNode : HenselNode) (m : ZZ)
    (hh : factorNode.h.isEmpty = false)
    (hinvert : (ZZ.invert 0 factorNode.h[0]!.2 m).1 = true) :
    ∃ output,
      Generated.StrictHensel.__hensel_step_bezout_phase_raw_ir
          (strictHenselRawOps termination) factorNode m = .ok output := by
  simp [Generated.StrictHensel.__hensel_step_bezout_phase_raw_ir,
    strictHenselRawOps, Generated.StrictHensel.__upoly_mul_raw_ir,
    Generated.StrictHensel.__upoly_add_raw_ir,
    Generated.StrictHensel.__upoly_sub_raw_ir,
    Generated.StrictHensel.__upoly_divmod_mod_raw_ir, hh, hinvert,
    Generated.StrictHensel.__upoly_mod_coeff_raw_ir]
  exact ⟨_, rfl⟩

/-- Safety-only invariant at the proof-visible boundary between the two
contiguous statement ranges.  It contains no output specification and cannot
manufacture a result: both fields merely discharge the second concrete C++
division assertion for the uniquely produced factor-phase node. -/
structure HenselStepExecutionInvariant
    (termination : Generated.StrictHensel.DivmodTermination)
    (node : HenselNode) (f : SparsePolyZZ) (m : ZZ) : Prop where
  inputHNonempty : node.h.isEmpty = false
  inputHInvertible : (ZZ.invert 0 node.h[0]!.2 m).1 = true
  factorHReady : ∀ factorNode,
    Generated.StrictHensel.__hensel_step_factor_phase_raw_ir
        (strictHenselRawOps termination) node f m = .ok factorNode →
    factorNode.h.isEmpty = false ∧
      (ZZ.invert 0 factorNode.h[0]!.2 m).1 = true

set_option maxHeartbeats 0 in
/-- Complete raw-to-safe bridge for the generated `__hensel_step` entry.
The witness is obtained only by executing the two exact raw source phases. -/
theorem __hensel_step_raw_ir_terminates
    (termination : Generated.StrictHensel.DivmodTermination)
    (node : HenselNode) (f : SparsePolyZZ) (m : ZZ)
    (hinvariant : HenselStepExecutionInvariant termination node f m) :
    ∃ output,
      Generated.StrictHensel.__hensel_step_raw_ir
          (strictHenselRawOps termination) node f m = .ok output := by
  rcases __hensel_step_factor_phase_raw_ir_terminates termination node f m
      hinvariant.inputHNonempty hinvariant.inputHInvertible with
    ⟨factorNode, hfactor⟩
  have hready := hinvariant.factorHReady factorNode hfactor
  rcases __hensel_step_bezout_phase_raw_ir_terminates termination factorNode m
      hready.1 hready.2 with ⟨output, hbezout⟩
  refine ⟨output, ?_⟩
  change (Generated.StrictHensel.__hensel_step_factor_phase_raw_ir
      (strictHenselRawOps termination) node f m >>= fun factorNode =>
        Generated.StrictHensel.__hensel_step_bezout_phase_raw_ir
          (strictHenselRawOps termination) factorNode m) = .ok output
  rw [hfactor]
  exact hbezout

def liftChildMatches (field : Int32) :
    Option Generated.StrictHensel.HenselLiftTree → Prop
  | none => field = -1
  | some child => field.toInt = child.rootIndex

/-- Safety invariant for the exact recursive C++ tree traversal.  Its
recursive premises are universally quantified over the outputs of the
concrete raw calls, so it cannot select or manufacture an execution result.
It records only array bounds, source branch agreement, and the safety facts
needed by the two modular divisions inside each `__hensel_step`. -/
noncomputable def HenselLiftRecursiveExecutionInvariant
    (termination : Generated.StrictHensel.DivmodTermination) :
    (tree : Generated.StrictHensel.HenselLiftTree) →
      Array HenselNode → SparsePolyZZ → ZZ → Prop :=
  Generated.StrictHensel.HenselLiftTree.rec
    (motive_1 := fun _ => Array HenselNode → SparsePolyZZ → ZZ → Prop)
    (motive_2 := fun _ => Array HenselNode → SparsePolyZZ → ZZ → Prop)
    (fun index left right leftInvariant rightInvariant nodes target m =>
      ∃ node,
        nodes[index]? = some node ∧
        liftChildMatches node.left left ∧
        liftChildMatches node.right right ∧
        HenselStepExecutionInvariant termination node target m ∧
        ∀ lifted,
          Generated.StrictHensel.__hensel_step_raw_ir
              (strictHenselRawOps termination) node target m = .ok lifted →
          let stored := nodes.set! index lifted
          leftInvariant stored lifted.g m ∧
          ∀ nodesAfterLeft,
            (match left with
              | none => (.ok stored : RawExec (Array HenselNode))
              | some child =>
                  Generated.StrictHensel.__hensel_lift_recursive_raw_ir
                    (strictHenselRawOps termination) child stored lifted.g m) =
                .ok nodesAfterLeft →
            ∃ parent,
              nodesAfterLeft[index]? = some parent ∧
              liftChildMatches parent.right right ∧
              rightInvariant nodesAfterLeft parent.h m)
    (fun _ _ _ => True)
    (fun _ childInvariant => childInvariant)

/-- Source-level composition rule for one recursive node.  Splitting the four
concrete left/right branch combinations keeps the proof independent of any L2
semantics and exposes exactly the generated monadic control flow. -/
theorem henselLiftRecursive_run_of_parts
    (termination : Generated.StrictHensel.DivmodTermination)
    (index : Nat)
    (left right : Option Generated.StrictHensel.HenselLiftTree)
    (nodes nodesAfterLeft output : Array HenselNode)
    (target : SparsePolyZZ) (m : ZZ) (node lifted parent : HenselNode)
    (hnode : nodes[index]? = some node)
    (hstep : Generated.StrictHensel.__hensel_step_raw_ir
      (strictHenselRawOps termination) node target m = .ok lifted)
    (hleft :
      (match left with
        | none => (.ok (nodes.set! index lifted) : RawExec (Array HenselNode))
        | some child =>
            Generated.StrictHensel.__hensel_lift_recursive_raw_ir
              (strictHenselRawOps termination) child (nodes.set! index lifted)
                lifted.g m) = .ok nodesAfterLeft)
    (hparent : nodesAfterLeft[index]? = some parent)
    (hright :
      (match right with
        | none => (.ok nodesAfterLeft : RawExec (Array HenselNode))
        | some child =>
            Generated.StrictHensel.__hensel_lift_recursive_raw_ir
              (strictHenselRawOps termination) child nodesAfterLeft parent.h
                m) = .ok output) :
    Generated.StrictHensel.__hensel_lift_recursive_raw_ir
      (strictHenselRawOps termination) (.node index left right) nodes target m =
        .ok output := by
  cases left <;> cases right <;>
    simp_all [Generated.StrictHensel.__hensel_lift_recursive_raw_ir.eq_1,
      bind, Except.instMonad, Except.bind]

/-- Genuine raw-to-safe bridge for `__hensel_lift_recursive`.  Termination is
structural on the finite source tree certificate; there is no fuel parameter.
Every recursive call is the exact generated call appearing in the C++ order. -/
private theorem henselLiftRecursiveTerminatesAux
    (termination : Generated.StrictHensel.DivmodTermination)
    (m : ZZ) (tree : Generated.StrictHensel.HenselLiftTree)
    (nodes : Array HenselNode) (target : SparsePolyZZ) :
    HenselLiftRecursiveExecutionInvariant termination tree nodes target m →
    ∃ output,
      Generated.StrictHensel.__hensel_lift_recursive_raw_ir
          (strictHenselRawOps termination) tree nodes target m = .ok output := by
  intro hinvariant
  cases tree with
  | node index left right =>
      simp only [HenselLiftRecursiveExecutionInvariant] at hinvariant
      rcases hinvariant with
        ⟨node, hnode, hleftMatch, hrightMatch, hstepSafe, hchildren⟩
      rcases __hensel_step_raw_ir_terminates termination node target m
          hstepSafe with ⟨lifted, hstep⟩
      have hchildReady := hchildren lifted hstep
      let stored := nodes.set! index lifted
      have hleftRun : ∃ nodesAfterLeft,
          (match left with
            | none => (.ok stored : RawExec (Array HenselNode))
            | some child =>
                Generated.StrictHensel.__hensel_lift_recursive_raw_ir
                  (strictHenselRawOps termination) child stored lifted.g m) =
              .ok nodesAfterLeft := by
        cases left with
        | none => exact ⟨stored, rfl⟩
        | some child =>
            exact henselLiftRecursiveTerminatesAux termination m child stored
              lifted.g hchildReady.1
      rcases hleftRun with ⟨nodesAfterLeft, hleftRun⟩
      have hparentReady : ∃ parent,
          nodesAfterLeft[index]? = some parent ∧
          liftChildMatches parent.right right ∧
          match right with
          | none => True
          | some child =>
              HenselLiftRecursiveExecutionInvariant termination child
                nodesAfterLeft parent.h m := by
        cases left <;> cases right <;>
          simpa [HenselLiftRecursiveExecutionInvariant] using
            hchildReady.2 nodesAfterLeft hleftRun
      rcases hparentReady with
        ⟨parent, hparent, hrightStillMatches, hrightReady⟩
      have hrightRun : ∃ output,
          (match right with
            | none => (.ok nodesAfterLeft : RawExec (Array HenselNode))
            | some child =>
                Generated.StrictHensel.__hensel_lift_recursive_raw_ir
                  (strictHenselRawOps termination) child nodesAfterLeft
                    parent.h m) = .ok output := by
        cases right with
        | none => exact ⟨nodesAfterLeft, rfl⟩
        | some child =>
            exact henselLiftRecursiveTerminatesAux termination m child
              nodesAfterLeft parent.h hrightReady
      rcases hrightRun with ⟨output, hrightRun⟩
      refine ⟨output, ?_⟩
      exact henselLiftRecursive_run_of_parts termination index left right
        nodes nodesAfterLeft output target m node lifted parent hnode hstep
        hleftRun hparent hrightRun
termination_by Generated.StrictHensel.HenselLiftTree.nodeCount tree
decreasing_by
  all_goals subst tree
  all_goals simp_all [Generated.StrictHensel.HenselLiftTree.nodeCount]
  all_goals omega

theorem __hensel_lift_recursive_raw_ir_terminates
    (termination : Generated.StrictHensel.DivmodTermination)
    (tree : Generated.StrictHensel.HenselLiftTree)
    (nodes : Array HenselNode) (target : SparsePolyZZ) (m : ZZ)
    (hinvariant : HenselLiftRecursiveExecutionInvariant termination tree
      nodes target m) :
    ∃ output,
      Generated.StrictHensel.__hensel_lift_recursive_raw_ir
          (strictHenselRawOps termination) tree nodes target m = .ok output :=
  henselLiftRecursiveTerminatesAux termination m tree nodes target hinvariant

/-- Semantic trace of one complete recursive Hensel tree traversal.  Every
constructor is anchored to a concrete generated raw-step equation and carries
its L2 `HenselStepCorrect` proof.  Recursive premises describe the exact left
and right raw calls, so this relation cannot validate a substituted tree. -/
inductive HenselLiftRecursiveCorrect
    (termination : Generated.StrictHensel.DivmodTermination) (m : Nat) :
    (tree : Generated.StrictHensel.HenselLiftTree) →
      Array HenselNode → SparsePolyZZ → Array HenselNode → Prop
  | leaf
      (index : Nat) (nodes stored : Array HenselNode)
      (target : SparsePolyZZ) (inputNode lifted parent : HenselNode)
      (hnode : nodes[index]? = some inputNode)
      (hstep : Generated.StrictHensel.__hensel_step_raw_ir
        (strictHenselRawOps termination) inputNode target (m : Int) =
          .ok lifted)
      (hstepCorrect : HenselStepCorrect target m inputNode lifted)
      (hstored : stored = nodes.set! index lifted)
      (hparent : stored[index]? = some parent) :
      HenselLiftRecursiveCorrect termination m (.node index none none)
        nodes target stored
  | left
      (index : Nat) (left : Generated.StrictHensel.HenselLiftTree)
      (nodes stored nodesAfterLeft : Array HenselNode)
      (target : SparsePolyZZ) (inputNode lifted parent : HenselNode)
      (hnode : nodes[index]? = some inputNode)
      (hstep : Generated.StrictHensel.__hensel_step_raw_ir
        (strictHenselRawOps termination) inputNode target (m : Int) =
          .ok lifted)
      (hstepCorrect : HenselStepCorrect target m inputNode lifted)
      (hstored : stored = nodes.set! index lifted)
      (hleftRun :
        Generated.StrictHensel.__hensel_lift_recursive_raw_ir
          (strictHenselRawOps termination) left stored lifted.g (m : Int) =
            .ok nodesAfterLeft)
      (hleftCorrect : HenselLiftRecursiveCorrect termination m left stored
        lifted.g nodesAfterLeft)
      (hparent : nodesAfterLeft[index]? = some parent) :
      HenselLiftRecursiveCorrect termination m (.node index (some left) none)
        nodes target nodesAfterLeft
  | right
      (index : Nat) (right : Generated.StrictHensel.HenselLiftTree)
      (nodes stored output : Array HenselNode)
      (target : SparsePolyZZ) (inputNode lifted parent : HenselNode)
      (hnode : nodes[index]? = some inputNode)
      (hstep : Generated.StrictHensel.__hensel_step_raw_ir
        (strictHenselRawOps termination) inputNode target (m : Int) =
          .ok lifted)
      (hstepCorrect : HenselStepCorrect target m inputNode lifted)
      (hstored : stored = nodes.set! index lifted)
      (hparent : stored[index]? = some parent)
      (hrightRun :
        Generated.StrictHensel.__hensel_lift_recursive_raw_ir
          (strictHenselRawOps termination) right stored parent.h (m : Int) =
            .ok output)
      (hrightCorrect : HenselLiftRecursiveCorrect termination m right stored
        parent.h output) :
      HenselLiftRecursiveCorrect termination m (.node index none (some right))
        nodes target output
  | branch
      (index : Nat)
      (left right : Generated.StrictHensel.HenselLiftTree)
      (nodes stored nodesAfterLeft output : Array HenselNode)
      (target : SparsePolyZZ) (inputNode lifted parent : HenselNode)
      (hnode : nodes[index]? = some inputNode)
      (hstep : Generated.StrictHensel.__hensel_step_raw_ir
        (strictHenselRawOps termination) inputNode target (m : Int) =
          .ok lifted)
      (hstepCorrect : HenselStepCorrect target m inputNode lifted)
      (hstored : stored = nodes.set! index lifted)
      (hleftRun :
        Generated.StrictHensel.__hensel_lift_recursive_raw_ir
          (strictHenselRawOps termination) left stored lifted.g (m : Int) =
            .ok nodesAfterLeft)
      (hleftCorrect : HenselLiftRecursiveCorrect termination m left stored
        lifted.g nodesAfterLeft)
      (hparent : nodesAfterLeft[index]? = some parent)
      (hrightRun :
        Generated.StrictHensel.__hensel_lift_recursive_raw_ir
          (strictHenselRawOps termination) right nodesAfterLeft parent.h
            (m : Int) = .ok output)
      (hrightCorrect : HenselLiftRecursiveCorrect termination m right
        nodesAfterLeft parent.h output) :
      HenselLiftRecursiveCorrect termination m
        (.node index (some left) (some right)) nodes target output

/-- Full refinement invariant for a recursive tree traversal.  Like the
single-step invariant, all descendant premises are universal over uniquely
determined raw outputs.  The invariant supplies representation/safety facts,
but never an output polynomial or a replacement L2 computation. -/
noncomputable def HenselLiftRecursiveRefinementInvariant
    (termination : Generated.StrictHensel.DivmodTermination) :
    (tree : Generated.StrictHensel.HenselLiftTree) →
      Array HenselNode → SparsePolyZZ → Nat → Prop :=
  Generated.StrictHensel.HenselLiftTree.rec
    (motive_1 := fun _ => Array HenselNode → SparsePolyZZ → Nat → Prop)
    (motive_2 := fun _ => Array HenselNode → SparsePolyZZ → Nat → Prop)
    (fun index left right leftInvariant rightInvariant nodes target m =>
      ∃ node,
        nodes[index]? = some node ∧
        liftChildMatches node.left left ∧
        liftChildMatches node.right right ∧
        HenselStepRefinementInvariant termination node target m ∧
        ∀ lifted,
          Generated.StrictHensel.__hensel_step_raw_ir
              (strictHenselRawOps termination) node target (m : Int) =
                .ok lifted →
          let stored := nodes.set! index lifted
          leftInvariant stored lifted.g m ∧
          ∀ nodesAfterLeft,
            (match left with
              | none => (.ok stored : RawExec (Array HenselNode))
              | some child =>
                  Generated.StrictHensel.__hensel_lift_recursive_raw_ir
                    (strictHenselRawOps termination) child stored lifted.g
                      (m : Int)) = .ok nodesAfterLeft →
            ∃ parent,
              nodesAfterLeft[index]? = some parent ∧
              liftChildMatches parent.right right ∧
              rightInvariant nodesAfterLeft parent.h m)
    (fun _ _ _ => True)
    (fun _ childInvariant => childInvariant)

set_option maxHeartbeats 0 in
/-- Strict L1→L2 refinement of the original C++
`__hensel_lift_recursive`.  The returned array is exactly the result of the
generated tree program, and the semantic trace proves that every concrete
node update is a quadratic Hensel refinement. -/
private theorem henselLiftRecursiveRefinesAux
    (termination : Generated.StrictHensel.DivmodTermination)
    (m : Nat) (tree : Generated.StrictHensel.HenselLiftTree)
    (nodes : Array HenselNode) (target : SparsePolyZZ) :
    HenselLiftRecursiveRefinementInvariant termination tree nodes target m →
    ∃ output,
      Generated.StrictHensel.__hensel_lift_recursive_raw_ir
          (strictHenselRawOps termination) tree nodes target (m : Int) =
            .ok output ∧
      HenselLiftRecursiveCorrect termination m tree nodes target output := by
  intro hinvariant
  cases tree with
  | node index left right =>
      simp only [HenselLiftRecursiveRefinementInvariant] at hinvariant
      rcases hinvariant with
        ⟨node, hnode, hleftMatch, hrightMatch, hstepInvariant, hchildren⟩
      rcases __hensel_step_raw_ir_refines termination node target m
          hstepInvariant with ⟨lifted, hstep, hstepCorrect⟩
      have hchildReady := hchildren lifted hstep
      let stored := nodes.set! index lifted
      have hleftResult : ∃ nodesAfterLeft,
          (match left with
            | none => (.ok stored : RawExec (Array HenselNode))
            | some child =>
                Generated.StrictHensel.__hensel_lift_recursive_raw_ir
                  (strictHenselRawOps termination) child stored lifted.g
                    (m : Int)) = .ok nodesAfterLeft ∧
          match left with
          | none => nodesAfterLeft = stored
          | some child =>
              HenselLiftRecursiveCorrect termination m child stored lifted.g
                nodesAfterLeft := by
        cases left with
        | none => exact ⟨stored, rfl, rfl⟩
        | some child =>
            rcases henselLiftRecursiveRefinesAux termination m child stored
                lifted.g hchildReady.1 with
              ⟨nodesAfterLeft, hrun, hcorrect⟩
            exact ⟨nodesAfterLeft, hrun, hcorrect⟩
      rcases hleftResult with
        ⟨nodesAfterLeft, hleftRun, hleftCorrect⟩
      have hparentReady : ∃ parent,
          nodesAfterLeft[index]? = some parent ∧
          liftChildMatches parent.right right ∧
          match right with
          | none => True
          | some child =>
              HenselLiftRecursiveRefinementInvariant termination child
                nodesAfterLeft parent.h m := by
        cases left <;> cases right <;>
          simpa [HenselLiftRecursiveRefinementInvariant] using
            hchildReady.2 nodesAfterLeft hleftRun
      rcases hparentReady with
        ⟨parent, hparent, hrightStillMatches, hrightReady⟩
      have hrightResult : ∃ output,
          (match right with
            | none => (.ok nodesAfterLeft : RawExec (Array HenselNode))
            | some child =>
                Generated.StrictHensel.__hensel_lift_recursive_raw_ir
                  (strictHenselRawOps termination) child nodesAfterLeft
                    parent.h (m : Int)) = .ok output ∧
          match right with
          | none => output = nodesAfterLeft
          | some child =>
              HenselLiftRecursiveCorrect termination m child nodesAfterLeft
                parent.h output := by
        cases right with
        | none => exact ⟨nodesAfterLeft, rfl, rfl⟩
        | some child =>
            rcases henselLiftRecursiveRefinesAux termination m child
                nodesAfterLeft parent.h hrightReady with
              ⟨output, hrun, hcorrect⟩
            exact ⟨output, hrun, hcorrect⟩
      rcases hrightResult with ⟨output, hrightRun, hrightCorrect⟩
      have hfullRun :
          Generated.StrictHensel.__hensel_lift_recursive_raw_ir
            (strictHenselRawOps termination) (.node index left right) nodes
              target (m : Int) = .ok output := by
        exact henselLiftRecursive_run_of_parts termination index left right
          nodes nodesAfterLeft output target (m : Int) node lifted parent hnode
          hstep hleftRun hparent hrightRun
      refine ⟨output, hfullRun, ?_⟩
      cases left with
      | none =>
          simp only at hleftCorrect
          subst nodesAfterLeft
          cases right with
          | none =>
              simp only at hrightCorrect
              subst output
              exact .leaf index nodes stored target node lifted parent hnode
                hstep hstepCorrect rfl hparent
          | some right =>
              exact .right index right nodes stored output target node lifted
                parent hnode hstep hstepCorrect rfl hparent hrightRun
                hrightCorrect
      | some left =>
          cases right with
          | none =>
              simp only at hrightCorrect
              subst output
              exact .left index left nodes stored nodesAfterLeft target node
                lifted parent hnode hstep hstepCorrect rfl hleftRun
                hleftCorrect hparent
          | some right =>
              exact .branch index left right nodes stored nodesAfterLeft output
                target node lifted parent hnode hstep hstepCorrect rfl
                hleftRun hleftCorrect hparent hrightRun hrightCorrect
termination_by Generated.StrictHensel.HenselLiftTree.nodeCount tree
decreasing_by
  all_goals subst tree
  all_goals simp_all [Generated.StrictHensel.HenselLiftTree.nodeCount]
  all_goals omega

theorem __hensel_lift_recursive_raw_ir_refines
    (termination : Generated.StrictHensel.DivmodTermination)
    (tree : Generated.StrictHensel.HenselLiftTree)
    (nodes : Array HenselNode) (target : SparsePolyZZ) (m : Nat)
    (hinvariant : HenselLiftRecursiveRefinementInvariant termination tree
      nodes target m) :
    ∃ output,
      Generated.StrictHensel.__hensel_lift_recursive_raw_ir
          (strictHenselRawOps termination) tree nodes target (m : Int) =
            .ok output ∧
      HenselLiftRecursiveCorrect termination m tree nodes target output :=
  henselLiftRecursiveRefinesAux termination m tree nodes target hinvariant

/-- A concrete prefix of the source quadratic-precision loop.  Every edge is
anchored to an exact generated tree-traversal equation; the relation contains
no predicted output or L2 computation. -/
inductive HenselLiftLoopPrefix
    (termination : Generated.StrictHensel.DivmodTermination)
    (tree : Generated.StrictHensel.HenselLiftTree) (f : SparsePolyZZ)
    (target : Nat) (initialM : Nat) (initialNodes : Array HenselNode) :
    Nat → Array HenselNode → Prop
  | refl : HenselLiftLoopPrefix termination tree f target initialM initialNodes
      initialM initialNodes
  | step
      (m : Nat) (nodes nextNodes : Array HenselNode)
      (hprefix : HenselLiftLoopPrefix termination tree f target initialM
        initialNodes m nodes)
      (hcontinue : m ≤ target)
      (hrun : Generated.StrictHensel.__hensel_lift_recursive_raw_ir
        (strictHenselRawOps termination) tree nodes f (m : Int) =
          .ok nextNodes) :
      HenselLiftLoopPrefix termination tree f target initialM initialNodes
        (m * m) nextNodes

/-- Safety assumptions for all states actually reachable through the exact
raw loop.  Universality over `HenselLiftLoopPrefix` prevents the invariant
from choosing a convenient sequence of node arrays. -/
structure HenselLiftLoopExecutionInvariant
    (termination : Generated.StrictHensel.DivmodTermination)
    (tree : Generated.StrictHensel.HenselLiftTree) (f : SparsePolyZZ)
    (target initialM : Nat) (initialNodes : Array HenselNode) : Prop where
  initialM_ge_two : 2 ≤ initialM
  iterationReady : ∀ m nodes,
    HenselLiftLoopPrefix termination tree f target initialM initialNodes
      m nodes →
    m ≤ target →
    HenselLiftRecursiveExecutionInvariant termination tree nodes f (m : Int)

private theorem henselLiftLoopTerminatesAux
    (termination : Generated.StrictHensel.DivmodTermination)
    (tree : Generated.StrictHensel.HenselLiftTree) (f : SparsePolyZZ)
    (target initialM : Nat) (initialNodes : Array HenselNode)
    (hinvariant : HenselLiftLoopExecutionInvariant termination tree f target
      initialM initialNodes)
    (m : Nat) (hm : 2 ≤ m) (nodes : Array HenselNode)
    (hprefix : HenselLiftLoopPrefix termination tree f target initialM
      initialNodes m nodes) :
    ∃ output,
      Generated.StrictHensel.__hensel_lift_loop_raw_ir
        (strictHenselRawOps termination) tree f target m hm nodes =
          .ok output := by
  rw [Generated.StrictHensel.__hensel_lift_loop_raw_ir]
  split
  · rename_i hcontinue
    rcases __hensel_lift_recursive_raw_ir_terminates termination tree nodes f
        (m : Int) (hinvariant.iterationReady m nodes hprefix hcontinue) with
      ⟨nextNodes, hrun⟩
    simp only [hrun, bind, Except.bind]
    have hnextM : 2 ≤ m * m := by
      have hmul := Nat.mul_le_mul_left m hm
      omega
    exact henselLiftLoopTerminatesAux termination tree f target initialM
      initialNodes hinvariant (m * m) hnextM nextNodes
      (.step m nodes nextNodes hprefix hcontinue hrun)
  · exact ⟨(nodes, m), rfl⟩
termination_by target + 1 - m
decreasing_by
  have hmul := Nat.mul_le_mul_left m hm
  have hgrow : m < m * m := by omega
  omega

theorem __hensel_lift_loop_raw_ir_terminates
    (termination : Generated.StrictHensel.DivmodTermination)
    (tree : Generated.StrictHensel.HenselLiftTree) (f : SparsePolyZZ)
    (target initialM : Nat) (initialNodes : Array HenselNode)
    (hinvariant : HenselLiftLoopExecutionInvariant termination tree f target
      initialM initialNodes) :
    ∃ output,
      Generated.StrictHensel.__hensel_lift_loop_raw_ir
        (strictHenselRawOps termination) tree f target initialM
          hinvariant.initialM_ge_two initialNodes = .ok output :=
  henselLiftLoopTerminatesAux termination tree f target initialM initialNodes
    hinvariant initialM hinvariant.initialM_ge_two initialNodes .refl

/-- Semantic trace for the exact outer precision loop. -/
inductive HenselLiftLoopCorrect
    (termination : Generated.StrictHensel.DivmodTermination)
    (tree : Generated.StrictHensel.HenselLiftTree) (f : SparsePolyZZ)
    (target : Nat) :
    Nat → Array HenselNode → Array HenselNode → Nat → Prop
  | done
      (m : Nat) (nodes : Array HenselNode) (hstop : ¬ m ≤ target) :
      HenselLiftLoopCorrect termination tree f target m nodes nodes m
  | step
      (m : Nat) (nodes nextNodes outputNodes : Array HenselNode)
      (outputM : Nat) (hcontinue : m ≤ target)
      (hrun : Generated.StrictHensel.__hensel_lift_recursive_raw_ir
        (strictHenselRawOps termination) tree nodes f (m : Int) =
          .ok nextNodes)
      (hiteration : HenselLiftRecursiveCorrect termination m tree nodes f
        nextNodes)
      (htail : HenselLiftLoopCorrect termination tree f target (m * m)
        nextNodes outputNodes outputM) :
      HenselLiftLoopCorrect termination tree f target m nodes outputNodes
        outputM

structure HenselLiftLoopRefinementInvariant
    (termination : Generated.StrictHensel.DivmodTermination)
    (tree : Generated.StrictHensel.HenselLiftTree) (f : SparsePolyZZ)
    (target initialM : Nat) (initialNodes : Array HenselNode) : Prop where
  initialM_ge_two : 2 ≤ initialM
  iterationReady : ∀ m nodes,
    HenselLiftLoopPrefix termination tree f target initialM initialNodes
      m nodes →
    m ≤ target →
    HenselLiftRecursiveRefinementInvariant termination tree nodes f m

private theorem henselLiftLoopRefinesAux
    (termination : Generated.StrictHensel.DivmodTermination)
    (tree : Generated.StrictHensel.HenselLiftTree) (f : SparsePolyZZ)
    (target initialM : Nat) (initialNodes : Array HenselNode)
    (hinvariant : HenselLiftLoopRefinementInvariant termination tree f target
      initialM initialNodes)
    (m : Nat) (hm : 2 ≤ m) (nodes : Array HenselNode)
    (hprefix : HenselLiftLoopPrefix termination tree f target initialM
      initialNodes m nodes) :
    ∃ outputNodes outputM,
      Generated.StrictHensel.__hensel_lift_loop_raw_ir
        (strictHenselRawOps termination) tree f target m hm nodes =
          .ok (outputNodes, outputM) ∧
      HenselLiftLoopCorrect termination tree f target m nodes outputNodes
        outputM := by
  rw [Generated.StrictHensel.__hensel_lift_loop_raw_ir]
  split
  · rename_i hcontinue
    rcases __hensel_lift_recursive_raw_ir_refines termination tree nodes f m
        (hinvariant.iterationReady m nodes hprefix hcontinue) with
      ⟨nextNodes, hrun, hiteration⟩
    simp only [hrun, bind, Except.bind]
    have hnextM : 2 ≤ m * m := by
      have hmul := Nat.mul_le_mul_left m hm
      omega
    rcases henselLiftLoopRefinesAux termination tree f target initialM
        initialNodes hinvariant (m * m) hnextM nextNodes
        (.step m nodes nextNodes hprefix hcontinue hrun) with
      ⟨outputNodes, outputM, htailRun, htailCorrect⟩
    exact ⟨outputNodes, outputM, htailRun,
      .step m nodes nextNodes outputNodes outputM hcontinue hrun hiteration
        htailCorrect⟩
  · rename_i hstop
    exact ⟨nodes, m, rfl, .done m nodes hstop⟩
termination_by target + 1 - m
decreasing_by
  have hmul := Nat.mul_le_mul_left m hm
  have hgrow : m < m * m := by omega
  omega

theorem __hensel_lift_loop_raw_ir_refines
    (termination : Generated.StrictHensel.DivmodTermination)
    (tree : Generated.StrictHensel.HenselLiftTree) (f : SparsePolyZZ)
    (target initialM : Nat) (initialNodes : Array HenselNode)
    (hinvariant : HenselLiftLoopRefinementInvariant termination tree f target
      initialM initialNodes) :
    ∃ outputNodes outputM,
      Generated.StrictHensel.__hensel_lift_loop_raw_ir
        (strictHenselRawOps termination) tree f target initialM
          hinvariant.initialM_ge_two initialNodes =
            .ok (outputNodes, outputM) ∧
      HenselLiftLoopCorrect termination tree f target initialM initialNodes
        outputNodes outputM :=
  henselLiftLoopRefinesAux termination tree f target initialM initialNodes
    hinvariant initialM hinvariant.initialM_ge_two initialNodes .refl

/-- Bounds and branch-shape invariant for the exact factor-extraction walk. -/
noncomputable def HenselExtractInvariant :
    Generated.StrictHensel.HenselLiftTree → Array HenselNode → Prop :=
  Generated.StrictHensel.HenselLiftTree.rec
    (motive_1 := fun _ => Array HenselNode → Prop)
    (motive_2 := fun _ => Array HenselNode → Prop)
    (fun index left right leftInvariant rightInvariant nodes =>
      ∃ node,
        nodes[index]? = some node ∧
        liftChildMatches node.left left ∧
        liftChildMatches node.right right ∧
        leftInvariant nodes ∧ rightInvariant nodes)
    (fun _ => True)
    (fun _ childInvariant => childInvariant)

theorem henselExtract_run_of_parts
    (index : Nat)
    (left right : Option Generated.StrictHensel.HenselLiftTree)
    (nodes : Array HenselNode)
    (factors factorsAfterLeft output : Array SparsePolyZZ)
    (node : HenselNode) (hnode : nodes[index]? = some node)
    (hleft :
      (match left with
        | none => (.ok (factors.push node.g) : RawExec (Array SparsePolyZZ))
        | some child =>
            Generated.StrictHensel.__hensel_extract_factors_raw_ir child nodes
              factors) = .ok factorsAfterLeft)
    (hright :
      (match right with
        | none =>
            (.ok (factorsAfterLeft.push node.h) :
              RawExec (Array SparsePolyZZ))
        | some child =>
            Generated.StrictHensel.__hensel_extract_factors_raw_ir child nodes
              factorsAfterLeft) = .ok output) :
    Generated.StrictHensel.__hensel_extract_factors_raw_ir
      (.node index left right) nodes factors = .ok output := by
  cases left <;> cases right <;>
    simp_all [Generated.StrictHensel.__hensel_extract_factors_raw_ir.eq_1,
      bind, Except.bind]

private theorem henselExtractTerminatesAux
    (tree : Generated.StrictHensel.HenselLiftTree)
    (nodes : Array HenselNode) (factors : Array SparsePolyZZ) :
    HenselExtractInvariant tree nodes →
    ∃ output,
      Generated.StrictHensel.__hensel_extract_factors_raw_ir tree nodes
        factors = .ok output := by
  intro hinvariant
  cases tree with
  | node index left right =>
      simp only [HenselExtractInvariant] at hinvariant
      rcases hinvariant with
        ⟨node, hnode, hleftMatch, hrightMatch, hleftReady, hrightReady⟩
      have hleftRun : ∃ factorsAfterLeft,
          (match left with
            | none =>
                (.ok (factors.push node.g) : RawExec (Array SparsePolyZZ))
            | some child =>
                Generated.StrictHensel.__hensel_extract_factors_raw_ir child
                  nodes factors) = .ok factorsAfterLeft := by
        cases left with
        | none => exact ⟨factors.push node.g, rfl⟩
        | some child =>
            exact henselExtractTerminatesAux child nodes factors hleftReady
      rcases hleftRun with ⟨factorsAfterLeft, hleftRun⟩
      have hrightRun : ∃ output,
          (match right with
            | none =>
                (.ok (factorsAfterLeft.push node.h) :
                  RawExec (Array SparsePolyZZ))
            | some child =>
                Generated.StrictHensel.__hensel_extract_factors_raw_ir child
                  nodes factorsAfterLeft) = .ok output := by
        cases right with
        | none => exact ⟨factorsAfterLeft.push node.h, rfl⟩
        | some child =>
            exact henselExtractTerminatesAux child nodes factorsAfterLeft
              hrightReady
      rcases hrightRun with ⟨output, hrightRun⟩
      exact ⟨output, henselExtract_run_of_parts index left right nodes factors
        factorsAfterLeft output node hnode hleftRun hrightRun⟩
termination_by Generated.StrictHensel.HenselLiftTree.nodeCount tree
decreasing_by
  all_goals subst tree
  all_goals simp_all [Generated.StrictHensel.HenselLiftTree.nodeCount]
  all_goals omega

theorem __hensel_extract_factors_raw_ir_terminates
    (tree : Generated.StrictHensel.HenselLiftTree)
    (nodes : Array HenselNode) (factors : Array SparsePolyZZ)
    (hinvariant : HenselExtractInvariant tree nodes) :
    ∃ output,
      Generated.StrictHensel.__hensel_extract_factors_raw_ir tree nodes
        factors = .ok output :=
  henselExtractTerminatesAux tree nodes factors hinvariant

/-- Exact semantic trace of leaf extraction.  Constructors distinguish all
four source branch combinations and retain the concrete generated child runs. -/
inductive HenselExtractCorrect :
    Generated.StrictHensel.HenselLiftTree → Array HenselNode →
      Array SparsePolyZZ → Array SparsePolyZZ → Prop
  | leaf
      (index : Nat) (nodes : Array HenselNode)
      (input output : Array SparsePolyZZ) (node : HenselNode)
      (hnode : nodes[index]? = some node)
      (houtput : output = (input.push node.g).push node.h) :
      HenselExtractCorrect (.node index none none) nodes input output
  | left
      (index : Nat) (left : Generated.StrictHensel.HenselLiftTree)
      (nodes : Array HenselNode) (input afterLeft output : Array SparsePolyZZ)
      (node : HenselNode) (hnode : nodes[index]? = some node)
      (hleftRun : Generated.StrictHensel.__hensel_extract_factors_raw_ir left
        nodes input = .ok afterLeft)
      (hleftCorrect : HenselExtractCorrect left nodes input afterLeft)
      (houtput : output = afterLeft.push node.h) :
      HenselExtractCorrect (.node index (some left) none) nodes input output
  | right
      (index : Nat) (right : Generated.StrictHensel.HenselLiftTree)
      (nodes : Array HenselNode) (input afterLeft output : Array SparsePolyZZ)
      (node : HenselNode) (hnode : nodes[index]? = some node)
      (hafterLeft : afterLeft = input.push node.g)
      (hrightRun : Generated.StrictHensel.__hensel_extract_factors_raw_ir right
        nodes afterLeft = .ok output)
      (hrightCorrect : HenselExtractCorrect right nodes afterLeft output) :
      HenselExtractCorrect (.node index none (some right)) nodes input output
  | branch
      (index : Nat) (left right : Generated.StrictHensel.HenselLiftTree)
      (nodes : Array HenselNode) (input afterLeft output : Array SparsePolyZZ)
      (node : HenselNode) (hnode : nodes[index]? = some node)
      (hleftRun : Generated.StrictHensel.__hensel_extract_factors_raw_ir left
        nodes input = .ok afterLeft)
      (hleftCorrect : HenselExtractCorrect left nodes input afterLeft)
      (hrightRun : Generated.StrictHensel.__hensel_extract_factors_raw_ir right
        nodes afterLeft = .ok output)
      (hrightCorrect : HenselExtractCorrect right nodes afterLeft output) :
      HenselExtractCorrect (.node index (some left) (some right)) nodes input
        output

private theorem henselExtractRefinesAux
    (tree : Generated.StrictHensel.HenselLiftTree)
    (nodes : Array HenselNode) (factors : Array SparsePolyZZ) :
    HenselExtractInvariant tree nodes →
    ∃ output,
      Generated.StrictHensel.__hensel_extract_factors_raw_ir tree nodes
          factors = .ok output ∧
      HenselExtractCorrect tree nodes factors output := by
  intro hinvariant
  cases tree with
  | node index left right =>
      simp only [HenselExtractInvariant] at hinvariant
      rcases hinvariant with
        ⟨node, hnode, hleftMatch, hrightMatch, hleftReady, hrightReady⟩
      have hleftResult : ∃ afterLeft,
          (match left with
            | none =>
                (.ok (factors.push node.g) : RawExec (Array SparsePolyZZ))
            | some child =>
                Generated.StrictHensel.__hensel_extract_factors_raw_ir child
                  nodes factors) = .ok afterLeft ∧
          match left with
          | none => afterLeft = factors.push node.g
          | some child => HenselExtractCorrect child nodes factors afterLeft := by
        cases left with
        | none => exact ⟨factors.push node.g, rfl, rfl⟩
        | some child =>
            rcases henselExtractRefinesAux child nodes factors hleftReady with
              ⟨afterLeft, hrun, hcorrect⟩
            exact ⟨afterLeft, hrun, hcorrect⟩
      rcases hleftResult with ⟨afterLeft, hleftRun, hleftCorrect⟩
      have hrightResult : ∃ output,
          (match right with
            | none =>
                (.ok (afterLeft.push node.h) : RawExec (Array SparsePolyZZ))
            | some child =>
                Generated.StrictHensel.__hensel_extract_factors_raw_ir child
                  nodes afterLeft) = .ok output ∧
          match right with
          | none => output = afterLeft.push node.h
          | some child => HenselExtractCorrect child nodes afterLeft output := by
        cases right with
        | none => exact ⟨afterLeft.push node.h, rfl, rfl⟩
        | some child =>
            rcases henselExtractRefinesAux child nodes afterLeft hrightReady with
              ⟨output, hrun, hcorrect⟩
            exact ⟨output, hrun, hcorrect⟩
      rcases hrightResult with ⟨output, hrightRun, hrightCorrect⟩
      have hfullRun := henselExtract_run_of_parts index left right nodes factors
        afterLeft output node hnode hleftRun hrightRun
      refine ⟨output, hfullRun, ?_⟩
      cases left with
      | none =>
          simp only at hleftCorrect
          subst afterLeft
          cases right with
          | none =>
              simp only at hrightCorrect
              exact .leaf index nodes factors output node hnode
                hrightCorrect
          | some right =>
              exact .right index right nodes factors (factors.push node.g)
                output node hnode rfl hrightRun hrightCorrect
      | some left =>
          cases right with
          | none =>
              simp only at hrightCorrect
              exact .left index left nodes factors afterLeft output node hnode
                hleftRun hleftCorrect hrightCorrect
          | some right =>
              exact .branch index left right nodes factors afterLeft output node
                hnode hleftRun hleftCorrect hrightRun hrightCorrect
termination_by Generated.StrictHensel.HenselLiftTree.nodeCount tree
decreasing_by
  all_goals subst tree
  all_goals simp_all [Generated.StrictHensel.HenselLiftTree.nodeCount]
  all_goals omega

theorem __hensel_extract_factors_raw_ir_refines
    (tree : Generated.StrictHensel.HenselLiftTree)
    (nodes : Array HenselNode) (factors : Array SparsePolyZZ)
    (hinvariant : HenselExtractInvariant tree nodes) :
    ∃ output,
      Generated.StrictHensel.__hensel_extract_factors_raw_ir tree nodes
          factors = .ok output ∧
      HenselExtractCorrect tree nodes factors output :=
  henselExtractRefinesAux tree nodes factors hinvariant

/-- Safety premises for the final normalization block of `__hensel_lift`.
They state only that the source `front()` and `ZZ::invert` assertions are
valid on the branch where they execute. -/
structure HenselNormalizeExecutionInvariant
    (result : Array SparsePolyZZ) (m : ZZ) : Prop where
  firstNonempty : ∀ first, result[0]? = some first → ∃ leading,
    first[0]? = some leading
  inverseSucceeds : ∀ first leading,
    result[0]? = some first → first[0]? = some leading → leading.2 != 1 →
    (ZZ.invert 0 leading.2 m).1 = true

/-- Exact semantic trace of the source normalization branch. -/
inductive HenselNormalizeCorrect
    (result : Array SparsePolyZZ) (m : ZZ) :
    Array SparsePolyZZ → Prop
  | empty (hresult : result[0]? = none) :
      HenselNormalizeCorrect result m result
  | alreadyOne (first : SparsePolyZZ) (leading : UMonomial × ZZ)
      (hresult : result[0]? = some first)
      (hfirst : first[0]? = some leading)
      (hone : ¬ (leading.2 != 1) = true) :
      HenselNormalizeCorrect result m result
  | normalized (first : SparsePolyZZ) (leading : UMonomial × ZZ)
      (inverse : ZZ) (scaled normalized : SparsePolyZZ)
      (hresult : result[0]? = some first)
      (hfirst : first[0]? = some leading) (hnotOne : (leading.2 != 1) = true)
      (hinverse : ZZ.invert 0 leading.2 m = (true, inverse))
      (hscaled : scaled = Generated.StrictHensel.scaleCoeffs first inverse)
      (hnormalized : Generated.StrictHensel.__upoly_mod_coeff_raw_ir scaled m =
        .ok normalized) :
      HenselNormalizeCorrect result m (result.set! 0 normalized)

/-- Genuine raw-to-safe and semantic refinement bridge for the final source
normalization block.  The output is obtained only by executing the strict raw
program; the invariant cannot supply an output polynomial. -/
theorem __hensel_normalize_result_raw_ir_refines
    (result : Array SparsePolyZZ) (m : ZZ)
    (hinvariant : HenselNormalizeExecutionInvariant result m) :
    ∃ output,
      Generated.StrictHensel.__hensel_normalize_result_raw_ir result m =
        .ok output ∧
      HenselNormalizeCorrect result m output := by
  rw [Generated.StrictHensel.__hensel_normalize_result_raw_ir]
  generalize hresult : result[0]? = firstOption
  cases firstOption with
  | none =>
      exact ⟨result, rfl, .empty hresult⟩
  | some first =>
      rcases hinvariant.firstNonempty first hresult with ⟨leading, hfirst⟩
      simp only
      rw [hfirst]
      simp only
      split
      · rename_i hnotOne
        have hinverseTrue :=
          hinvariant.inverseSucceeds first leading hresult hfirst hnotOne
        generalize hinverse : ZZ.invert 0 leading.2 m = inverseRun at hinverseTrue ⊢
        cases inverseRun with
        | mk succeeded inverse =>
            simp only at hinverseTrue
            subst succeeded
            simp only [ite_true, Generated.StrictHensel.__upoly_mod_coeff_raw_ir,
              bind, Except.bind]
            exact ⟨_, rfl, .normalized first leading inverse
              (Generated.StrictHensel.scaleCoeffs first inverse)
              (Generated.StrictHensel.modCoeffOutput
                (Generated.StrictHensel.scaleCoeffs first inverse) m)
              hresult hfirst hnotOne hinverse rfl rfl⟩
      · rename_i hone
        exact ⟨result, rfl, .alreadyOne first leading hresult hfirst hone⟩

/-- Bounds premises for the two source reads in the leading-coefficient
baking block. -/
structure HenselAdjustFirstFactorInvariant
    (f : SparsePolyZZ) (factors : Array SparsePolyZp) : Prop where
  sourceLeading : ∃ leading, f[0]? = some leading
  firstFactor : ∃ first, factors[0]? = some first

/-- Exact trace of multiplying factor zero by `lc(f) mod p`, normalizing it,
and writing it back. -/
inductive HenselAdjustFirstFactorCorrect
    (f : SparsePolyZZ) (factors : Array SparsePolyZp) (p : UInt64) :
    Array SparsePolyZp → Prop
  | adjusted (leading : UMonomial × ZZ) (first adjusted : SparsePolyZp)
      (hsource : f[0]? = some leading)
      (hfirst : factors[0]? = some first)
      (hadjusted : adjusted = SparsePolyZp.normalization
        (Generated.StrictHensel.scaleZpCoeffs first
          (Zp.ofInt leading.2 p))) :
      HenselAdjustFirstFactorCorrect f factors p
        (factors.set! 0 adjusted)

/-- Raw-to-safe semantic refinement for the exact coefficient-baking block. -/
theorem __hensel_adjust_first_factor_raw_ir_refines
    (f : SparsePolyZZ) (factors : Array SparsePolyZp) (p : UInt64)
    (hinvariant : HenselAdjustFirstFactorInvariant f factors) :
    ∃ output,
      Generated.StrictHensel.__hensel_adjust_first_factor_raw_ir
          f factors p = .ok output ∧
      HenselAdjustFirstFactorCorrect f factors p output := by
  rcases hinvariant.sourceLeading with ⟨leading, hsource⟩
  rcases hinvariant.firstFactor with ⟨first, hfirst⟩
  rw [Generated.StrictHensel.__hensel_adjust_first_factor_raw_ir,
    hsource, hfirst]
  exact ⟨_, rfl, .adjusted leading first _ hsource hfirst rfl⟩

theorem henselExplicitTargetLoop_eq
    (p : UInt64) (remaining : Nat) (target : ZZ) :
    Generated.StrictHensel.henselExplicitTargetLoop p remaining target =
      target * (p.toNat : ZZ) ^ remaining := by
  induction remaining generalizing target with
  | zero => simp [Generated.StrictHensel.henselExplicitTargetLoop]
  | succ remaining ih =>
      rw [Generated.StrictHensel.henselExplicitTargetLoop, ih]
      rw [pow_succ]
      ring

/-- L2 meaning of the explicit-precision source branch. -/
def HenselExplicitTargetCorrect
    (p : UInt64) (aTarget : Int32) (target : ZZ) : Prop :=
  target = (p.toNat : ZZ) ^ aTarget.toNatClampNeg - 1

/-- Exact refinement of the positive-`a_target` precision computation. -/
theorem __hensel_explicit_target_raw_ir_refines
    (p : UInt64) (aTarget : Int32) (hpositive : aTarget > 0) :
    ∃ target,
      Generated.StrictHensel.__hensel_explicit_target_raw_ir p aTarget =
        .ok target ∧
      HenselExplicitTargetCorrect p aTarget target := by
  refine ⟨Generated.StrictHensel.henselExplicitTargetLoop p
      aTarget.toNatClampNeg 1 - 1,
    ?_, ?_⟩
  · simp [Generated.StrictHensel.__hensel_explicit_target_raw_ir, hpositive]
  · rw [HenselExplicitTargetCorrect, henselExplicitTargetLoop_eq]
    simp

/-- Reachability through concrete successful source EEA iterations. -/
inductive HenselEEAPrefix (ops : Generated.StrictHensel.HenselEEARawOps)
    (initial : Generated.StrictHensel.HenselEEAState) :
    Generated.StrictHensel.HenselEEAState → Prop
  | refl : HenselEEAPrefix ops initial initial
  | step (state : Generated.StrictHensel.HenselEEAState)
      (quotient remainder : SparsePolyZp)
      (hprefix : HenselEEAPrefix ops initial state)
      (hcontinue : ¬ state.r1.isEmpty = true)
      (hrun : ops.divmod state.r0 state.r1 = .ok (quotient, remainder)) :
      HenselEEAPrefix ops initial
        (Generated.StrictHensel.henselEEANextState state quotient remainder)

/-- Safety assumptions for all states reachable by actual raw EEA division
results.  Neither field supplies the final GCD or Bézout coefficients. -/
structure HenselEEAExecutionInvariant
    (ops : Generated.StrictHensel.HenselEEARawOps)
    (initial : Generated.StrictHensel.HenselEEAState) : Prop where
  divisionReady : ∀ state,
    HenselEEAPrefix ops initial state → ¬ state.r1.isEmpty = true →
    ∃ quotient remainder,
      ops.divmod state.r0 state.r1 = .ok (quotient, remainder)
  finalNonempty : ∀ state,
    HenselEEAPrefix ops initial state → state.r1.isEmpty = true →
    ∃ leading, state.r0[0]? = some leading

/-- Exact semantic execution trace of the source extended-Euclidean loop. -/
inductive HenselEEACorrect
    (ops : Generated.StrictHensel.HenselEEARawOps) :
    Generated.StrictHensel.HenselEEAState →
      (SparsePolyZp × SparsePolyZp × SparsePolyZp) → Prop
  | done (state : Generated.StrictHensel.HenselEEAState)
      (leading : UMonomial × Zp) (hdone : state.r1.isEmpty = true)
      (hleading : state.r0[0]? = some leading) :
      HenselEEACorrect ops state
        (Generated.StrictHensel.henselEEAScaleNormalize leading.2.inv state.r0,
          Generated.StrictHensel.henselEEAScaleNormalize leading.2.inv state.s0,
          Generated.StrictHensel.henselEEAScaleNormalize leading.2.inv state.t0)
  | step (state : Generated.StrictHensel.HenselEEAState)
      (quotient remainder : SparsePolyZp)
      (output : SparsePolyZp × SparsePolyZp × SparsePolyZp)
      (hcontinue : ¬ state.r1.isEmpty = true)
      (hrun : ops.divmod state.r0 state.r1 = .ok (quotient, remainder))
      (htail : HenselEEACorrect ops
        (Generated.StrictHensel.henselEEANextState state quotient remainder)
        output) :
      HenselEEACorrect ops state output

private theorem henselEEARefinesAux
    (ops : Generated.StrictHensel.HenselEEARawOps)
    (termination : Generated.StrictHensel.HenselEEATermination ops)
    (initial state : Generated.StrictHensel.HenselEEAState)
    (hinvariant : HenselEEAExecutionInvariant ops initial)
    (hprefix : HenselEEAPrefix ops initial state) :
    ∃ output,
      Generated.StrictHensel.__polynomial_GCD_eea_raw_ir ops termination state =
        .ok output ∧
      HenselEEACorrect ops state output := by
  rw [Generated.StrictHensel.__polynomial_GCD_eea_raw_ir]
  split
  · rename_i hdone
    rcases hinvariant.finalNonempty state hprefix hdone with
      ⟨leading, hleading⟩
    rw [hleading]
    exact ⟨_, rfl, .done state leading hdone hleading⟩
  · rename_i hcontinue
    rcases hinvariant.divisionReady state hprefix hcontinue with
      ⟨quotient, remainder, hrun⟩
    rw [hrun]
    rcases henselEEARefinesAux ops termination initial
        (Generated.StrictHensel.henselEEANextState state quotient remainder)
        hinvariant (.step state quotient remainder hprefix hcontinue hrun) with
      ⟨output, htailRun, htailCorrect⟩
    exact ⟨output, htailRun,
      .step state quotient remainder output hcontinue hrun htailCorrect⟩
termination_by termination.measure state
decreasing_by
  apply termination.decreases
  · assumption
  · assumption

/-- Raw-to-safe semantic bridge for the strict well-founded EEA control flow.
This theorem is intentionally generic only in the executable raw divmod
boundary; tree construction will instantiate it with the strict generated
quotient/remainder implementation before publishing an entry contract. -/
theorem __polynomial_GCD_eea_raw_ir_refines
    (ops : Generated.StrictHensel.HenselEEARawOps)
    (termination : Generated.StrictHensel.HenselEEATermination ops)
    (initial : Generated.StrictHensel.HenselEEAState)
    (hinvariant : HenselEEAExecutionInvariant ops initial) :
    ∃ output,
      Generated.StrictHensel.__polynomial_GCD_eea_raw_ir ops termination
          initial = .ok output ∧
      HenselEEACorrect ops initial output :=
  henselEEARefinesAux ops termination initial initial hinvariant .refl

/-- State of the quotient+remainder form of the current generated VHC
`pair_vec_div` loop. -/
structure HenselDivmodVHCState where
  dividendIndex : Nat
  heap : Array Nat
  nodes : Array StrictSquarefreeZp.PairVecDivVHCNode
  quotient : SparsePolyZp
  remainder : SparsePolyZp
  resetH : Nat

def HenselDivmodVHCState.quotientState (state : HenselDivmodVHCState) :
    StrictSquarefreeZp.PairVecDivVHCIterationResult where
  dividendIndex := state.dividendIndex
  heap := state.heap
  nodes := state.nodes
  quotient := state.quotient
  resetH := state.resetH

/-- One exact body of the five-argument C++ `pair_vec_div`.  It reuses the
already strict selector, equal-degree heap consumption, quotient emission,
reset activation and reverse reinsertion.  The only additional observable
action is the source remainder push when the net frontier degree is below the
divisor leading degree. -/
def henselDivmodVHCIteration (this : DenseUPolyZp)
    (state : HenselDivmodVHCState)
    (frontier : StrictSquarefreeZp.PairVecDivVHCFrontier)
    (dividend divisor : SparsePolyZp) (hdivisor : 0 < divisor.size) :
    RawExec HenselDivmodVHCState := do
  let consumed ← StrictSquarefreeZp.pairVecDivVHCConsumeEqualDegree this
    frontier.degree state.heap frontier.coefficient state.nodes #[]
    state.resetH state.quotient divisor
  let remainder :=
    if consumed.coefficient ≠ 0 ∧ frontier.degree < divisor[0].1.deg then
      state.remainder.push
        (UMonomial.mk frontier.degree, ⟨consumed.coefficient, this._p⟩)
    else state.remainder
  let emitted ← StrictSquarefreeZp.pairVecDivVHCEmit this frontier consumed
    state.quotient divisor hdivisor
  let reinserted ← StrictSquarefreeZp.pairVecDivVHCReinsertLin emitted.2.1.heap
    emitted.2.1.nodes consumed.lin
  return ⟨frontier.dividendIndex, reinserted.heap, reinserted.nodes,
    emitted.1, remainder, emitted.2.2⟩

/-- Strict well-founded general VHC loop for the quotient+remainder overload.
The recursive degree bound is the degree actually selected from the concrete
dividend/heap frontier, exactly as for the already verified quotient-only
entry. -/
def henselDivmodVHCOuterLoop (this : DenseUPolyZp)
    (degreeLimit : Nat) (state : HenselDivmodVHCState)
    (dividend divisor : SparsePolyZp) (hdivisor : 0 < divisor.size) :
    RawExec (SparsePolyZp × SparsePolyZp) :=
  if hdone : dividend.size ≤ state.dividendIndex ∧ state.heap.size = 0 then
    .ok (state.quotient, state.remainder)
  else
    match StrictSquarefreeZp.pairVecDivVHCSelectFrontier state.dividendIndex
        dividend state.heap state.nodes with
    | .error fault => .error fault
    | .ok frontier =>
        if hdecrease : frontier.degree < degreeLimit then
          match henselDivmodVHCIteration this state frontier dividend divisor
              hdivisor with
          | .error fault => .error fault
          | .ok next =>
              henselDivmodVHCOuterLoop this frontier.degree next dividend
                divisor hdivisor
        else .error .assertionFailure
termination_by degreeLimit
decreasing_by exact hdecrease

/-- Erasing the accumulated remainder from one five-argument source
iteration gives exactly the already verified quotient-only VHC iteration.
This is a control-flow projection, not a semantic replacement: both sides
execute the same selector-independent heap body and use the same concrete
quotient emission result. -/
theorem henselDivmodVHCIteration_projects_quotient
    (this : DenseUPolyZp) (state next : HenselDivmodVHCState)
    (frontier : StrictSquarefreeZp.PairVecDivVHCFrontier)
    (dividend divisor : SparsePolyZp) (hdivisor : 0 < divisor.size)
    (hselect : StrictSquarefreeZp.pairVecDivVHCSelectFrontier
      state.dividendIndex dividend state.heap state.nodes = .ok frontier)
    (hrun : henselDivmodVHCIteration this state frontier dividend divisor
      hdivisor = .ok next) :
    ∃ projected,
      StrictSquarefreeZp.pairVecDivVHCOuterIteration this
          state.dividendIndex state.heap state.nodes state.quotient dividend
          divisor state.resetH = .ok projected ∧
      projected.dividendIndex = next.dividendIndex ∧
      projected.heap = next.heap ∧ projected.nodes = next.nodes ∧
      projected.quotient = next.quotient ∧ projected.resetH = next.resetH := by
  unfold henselDivmodVHCIteration at hrun
  unfold StrictSquarefreeZp.pairVecDivVHCOuterIteration
  simp only [hdivisor, ↓reduceDIte, hselect, Bind.bind, Except.bind]
  cases hconsume : StrictSquarefreeZp.pairVecDivVHCConsumeEqualDegree this
      frontier.degree state.heap frontier.coefficient state.nodes #[]
      state.resetH state.quotient divisor with
  | error fault => simp [hconsume, bind, Except.bind] at hrun
  | ok consumed =>
      simp only [hconsume, bind, Except.bind] at hrun ⊢
      cases hemit : StrictSquarefreeZp.pairVecDivVHCEmit this frontier consumed
          state.quotient divisor hdivisor with
      | error fault => simp [hemit, bind, Except.bind] at hrun
      | ok emitted =>
          simp only [hemit, bind, Except.bind] at hrun ⊢
          cases hreinsert : StrictSquarefreeZp.pairVecDivVHCReinsertLin
              emitted.2.1.heap emitted.2.1.nodes consumed.lin with
          | error fault => simp [hreinsert, bind, Except.bind] at hrun
          | ok reinserted =>
              simp only [hreinsert, bind, Except.bind] at hrun ⊢
              have hnext :
                  HenselDivmodVHCState.mk frontier.dividendIndex
                    reinserted.heap reinserted.nodes emitted.1
                    (if consumed.coefficient ≠ 0 ∧
                        frontier.degree < divisor[0].1.deg then
                      state.remainder.push
                        (UMonomial.mk frontier.degree,
                          ⟨consumed.coefficient, this._p⟩)
                    else state.remainder)
                    emitted.2.2 = next := by
                exact Except.ok.inj hrun
              subst next
              exact ⟨_, rfl, rfl, rfl, rfl, rfl, rfl⟩

theorem henselDivmodVHCIteration_projects_quotientState
    (this : DenseUPolyZp) (state next : HenselDivmodVHCState)
    (frontier : StrictSquarefreeZp.PairVecDivVHCFrontier)
    (dividend divisor : SparsePolyZp) (hdivisor : 0 < divisor.size)
    (hselect : StrictSquarefreeZp.pairVecDivVHCSelectFrontier
      state.dividendIndex dividend state.heap state.nodes = .ok frontier)
    (hrun : henselDivmodVHCIteration this state frontier dividend divisor
      hdivisor = .ok next) :
    StrictSquarefreeZp.pairVecDivVHCOuterIteration this state.dividendIndex
        state.heap state.nodes state.quotient dividend divisor state.resetH =
      .ok next.quotientState := by
  rcases henselDivmodVHCIteration_projects_quotient this state next frontier
      dividend divisor hdivisor hselect hrun with
    ⟨projected, hprojected, hindex, hheap, hnodes, hquotient, hreset⟩
  have heq : projected = next.quotientState := by
    cases projected
    cases next
    simp_all [HenselDivmodVHCState.quotientState]
  simpa [heq] using hprojected

/-- Conversely, every successful quotient-only body lifts to the exact
five-argument body by recording the source remainder push.  No operation is
re-executed semantically: the proof unfolds the shared consume, emit, and
reinsertion calls. -/
theorem henselDivmodVHCIteration_succeeds_of_quotient
    (this : DenseUPolyZp) (state : HenselDivmodVHCState)
    (frontier : StrictSquarefreeZp.PairVecDivVHCFrontier)
    (dividend divisor : SparsePolyZp) (hdivisor : 0 < divisor.size)
    (projected : StrictSquarefreeZp.PairVecDivVHCIterationResult)
    (hselect : StrictSquarefreeZp.pairVecDivVHCSelectFrontier
      state.dividendIndex dividend state.heap state.nodes = .ok frontier)
    (hrun : StrictSquarefreeZp.pairVecDivVHCOuterIteration this
      state.dividendIndex state.heap state.nodes state.quotient dividend
      divisor state.resetH = .ok projected) :
    ∃ next,
      henselDivmodVHCIteration this state frontier dividend divisor hdivisor =
          .ok next ∧
      projected.dividendIndex = next.dividendIndex ∧
      projected.heap = next.heap ∧ projected.nodes = next.nodes ∧
      projected.quotient = next.quotient ∧ projected.resetH = next.resetH := by
  unfold StrictSquarefreeZp.pairVecDivVHCOuterIteration at hrun
  simp only [hdivisor, ↓reduceDIte, hselect, bind, Except.bind] at hrun
  cases hconsume : StrictSquarefreeZp.pairVecDivVHCConsumeEqualDegree this
          frontier.degree state.heap frontier.coefficient state.nodes #[]
          state.resetH state.quotient divisor with
  | error fault => simp [hconsume, bind, Except.bind] at hrun
  | ok consumed =>
      simp only [hconsume, bind, Except.bind] at hrun
      cases hemit : StrictSquarefreeZp.pairVecDivVHCEmit this frontier
          consumed state.quotient divisor hdivisor with
      | error fault => simp [hemit, bind, Except.bind] at hrun
      | ok emitted =>
          simp only [hemit, bind, Except.bind] at hrun
          cases hreinsert : StrictSquarefreeZp.pairVecDivVHCReinsertLin
              emitted.2.1.heap emitted.2.1.nodes consumed.lin with
          | error fault => simp [hreinsert, bind, Except.bind] at hrun
          | ok reinserted =>
              simp only [hreinsert, bind, Except.bind] at hrun
              rw [← Except.ok.inj hrun]
              refine ⟨⟨frontier.dividendIndex, reinserted.heap,
                reinserted.nodes, emitted.1,
                if consumed.coefficient ≠ 0 ∧
                    frontier.degree < divisor[0].1.deg then
                  state.remainder.push (UMonomial.mk frontier.degree,
                    ⟨consumed.coefficient, this._p⟩)
                else state.remainder,
                emitted.2.2⟩, ?_, rfl, rfl, rfl, rfl, rfl⟩
              simp [henselDivmodVHCIteration, hconsume, hemit, hreinsert,
                bind, Except.bind]
              rfl

/-- A complete successful quotient-only VHC run lifts to a successful
five-argument quotient-and-remainder run with the identical concrete quotient.
The induction follows the same selected frontier and the same well-founded
degree decrease. -/
theorem henselDivmodVHCOuterLoop_succeeds_of_quotient
    (this : DenseUPolyZp) (dividend divisor : SparsePolyZp)
    (hdivisor : 0 < divisor.size) (limit : Nat)
    (state : HenselDivmodVHCState) (quotient : SparsePolyZp)
    (hrun : StrictSquarefreeZp.pairVecDivVHCOuterLoop this limit
      state.dividendIndex state.heap state.nodes state.quotient dividend
      divisor state.resetH = .ok quotient) :
    ∃ output,
      henselDivmodVHCOuterLoop this limit state dividend divisor hdivisor =
          .ok output ∧
      output.1 = quotient := by
  induction limit using Nat.strong_induction_on generalizing state quotient with
  | h limit ih =>
      rw [StrictSquarefreeZp.pairVecDivVHCOuterLoop] at hrun
      rw [henselDivmodVHCOuterLoop]
      by_cases hdone : dividend.size ≤ state.dividendIndex ∧
          state.heap.size = 0
      · rw [dif_pos hdone] at hrun ⊢
        have hquotient : state.quotient = quotient := Except.ok.inj hrun
        exact ⟨(state.quotient, state.remainder), rfl, hquotient⟩
      · rw [dif_neg hdone] at hrun ⊢
        cases hselect : StrictSquarefreeZp.pairVecDivVHCSelectFrontier
            state.dividendIndex dividend state.heap state.nodes with
        | error fault => simp [hselect] at hrun
        | ok frontier =>
            rw [hselect] at hrun
            simp only at hrun ⊢
            by_cases hdecrease : frontier.degree < limit
            · rw [dif_pos hdecrease] at hrun ⊢
              cases hiteration :
                  StrictSquarefreeZp.pairVecDivVHCOuterIteration this
                    state.dividendIndex state.heap state.nodes state.quotient
                    dividend divisor state.resetH with
              | error fault => simp [hiteration] at hrun
              | ok projected =>
                  rw [hiteration] at hrun
                  simp only at hrun
                  rcases henselDivmodVHCIteration_succeeds_of_quotient this
                      state frontier dividend divisor hdivisor projected hselect
                      hiteration with
                    ⟨next, hnextRun, hindex, hheap, hnodes, hquotient,
                      hreset⟩
                  rw [hnextRun]
                  simp only
                  rw [hindex, hheap, hnodes, hquotient, hreset] at hrun
                  rcases ih frontier.degree hdecrease next quotient hrun with
                    ⟨output, houtput, houtputQuotient⟩
                  exact ⟨output, houtput, houtputQuotient⟩
            · rw [dif_neg hdecrease] at hrun
              contradiction

/-- Fresh-output general branch of the C++ five-argument `pair_vec_div`
overload used by Hensel EEA. -/
def henselDivmodVHCGeneralBranchIR (this : DenseUPolyZp)
    (dividend divisor : SparsePolyZp) :
    RawExec (SparsePolyZp × SparsePolyZp) :=
  if hdividend : 0 < dividend.size then
    if hdivisor : 1 < divisor.size then
      let nodes := StrictSquarefreeZp.pairVecDivVHCInit divisor
      let state : HenselDivmodVHCState :=
        ⟨0, #[], nodes, #[], #[], divisor.size - 1⟩
      henselDivmodVHCOuterLoop this (dividend[0].1.deg + 1) state dividend
        divisor (by omega)
    else
      .error .assertionFailure
  else
    .error .assertionFailure

/-- The concrete general five-argument entry is total on canonical nonempty
inputs with a genuine divisor tail.  Its quotient execution is exactly the
already verified source VHC path; the additional output is accumulated by the
same bodies. -/
theorem henselDivmodVHCGeneralBranchIR_succeeds
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (dividend divisor : SparsePolyZp)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hdividendCanonical : CLPoly.Math.SparsePolyZp.Canonical this._p.toNat
      dividend)
    (hdivisorCanonical : CLPoly.Math.SparsePolyZp.Canonical this._p.toNat
      divisor)
    (hdividend : 0 < dividend.size) (hdivisor : 1 < divisor.size) :
    ∃ output, henselDivmodVHCGeneralBranchIR this dividend divisor =
      .ok output := by
  rcases StrictSquarefreeZp.pairVecDivGeneralBranchIR_succeeds this dividend
      divisor hcfg hdividendCanonical hdivisorCanonical hdividend hdivisor with
    ⟨quotient, hquotient⟩
  let nodes := StrictSquarefreeZp.pairVecDivVHCInit divisor
  let state : HenselDivmodVHCState :=
    ⟨0, #[], nodes, #[], #[], divisor.size - 1⟩
  have hquotientRun : StrictSquarefreeZp.pairVecDivVHCOuterLoop this
      (dividend[0].1.deg + 1) 0 #[] nodes #[] dividend divisor
        (divisor.size - 1) = .ok quotient := by
    simpa [StrictSquarefreeZp.pairVecDivGeneralBranchIR, hdividend, hdivisor,
      nodes] using hquotient
  rcases henselDivmodVHCOuterLoop_succeeds_of_quotient this dividend divisor
      (by omega) (dividend[0].1.deg + 1) state quotient (by
        simpa [state] using hquotientRun) with
    ⟨output, houtput, houtputQuotient⟩
  refine ⟨output, ?_⟩
  simpa [henselDivmodVHCGeneralBranchIR, hdividend, hdivisor, state, nodes]
    using houtput

/-- Exact source-order single-term-divisor loop for the five-argument
overload.  Over `Zp`, coefficient division sets its remainder output to zero,
so only monomials not divisible by the divisor monomial are appended to the
polynomial remainder. -/
def henselDivmodVHCSingleLoopIR (this : DenseUPolyZp) (index : Nat)
    (quotient remainder dividend : SparsePolyZp)
    (divisor : UMonomial × Zp) : SparsePolyZp × SparsePolyZp :=
  if h : index < dividend.size then
    match StrictSquarefreeZp.pairVecDivSingleTermIR this divisor
        dividend[index] with
    | some term =>
        henselDivmodVHCSingleLoopIR this (index + 1)
          (quotient.push term) remainder dividend divisor
    | none =>
        henselDivmodVHCSingleLoopIR this (index + 1) quotient
          (remainder.push dividend[index]) dividend divisor
  else
    (quotient, remainder)
termination_by dividend.size - index
decreasing_by all_goals omega

/-- Erasing the separately accumulated remainder gives the already checked
single-output source loop. -/
theorem henselDivmodVHCSingleLoopIR_fst
    (this : DenseUPolyZp) (index : Nat)
    (quotient remainder dividend : SparsePolyZp)
    (divisor : UMonomial × Zp) :
    (henselDivmodVHCSingleLoopIR this index quotient remainder dividend
      divisor).1 =
      StrictSquarefreeZp.pairVecDivSingleLoopIR this index quotient dividend
        divisor := by
  rw [henselDivmodVHCSingleLoopIR,
    StrictSquarefreeZp.pairVecDivSingleLoopIR]
  split
  next hmore =>
    cases hterm : StrictSquarefreeZp.pairVecDivSingleTermIR this divisor
        dividend[index] with
    | none =>
        simp only [hterm]
        exact henselDivmodVHCSingleLoopIR_fst this (index + 1) quotient
          (remainder.push dividend[index]) dividend divisor
    | some term =>
        simp only [hterm]
        exact henselDivmodVHCSingleLoopIR_fst this (index + 1)
          (quotient.push term) remainder dividend divisor
  next hdone => rfl
termination_by dividend.size - index
decreasing_by all_goals omega

def henselDivmodVHCSingleRemainderTermIR (this : DenseUPolyZp)
    (divisor term : UMonomial × Zp) : Option (UMonomial × Zp) :=
  match StrictSquarefreeZp.pairVecDivSingleTermIR this divisor term with
  | some _ => none
  | none => some term

theorem henselDivmodVHCSingleLoopIR_snd_toList
    (this : DenseUPolyZp) (index : Nat)
    (quotient remainder dividend : SparsePolyZp)
    (divisor : UMonomial × Zp) :
    (henselDivmodVHCSingleLoopIR this index quotient remainder dividend
      divisor).2.toList =
      remainder.toList ++ (dividend.toList.drop index).filterMap
        (henselDivmodVHCSingleRemainderTermIR this divisor) := by
  rw [henselDivmodVHCSingleLoopIR]
  split
  next hmore =>
    have hdrop := List.drop_eq_getElem_cons
      (l := dividend.toList) (i := index) (by simpa using hmore)
    rw [hdrop, List.filterMap_cons, Array.getElem_toList hmore]
    cases hterm : StrictSquarefreeZp.pairVecDivSingleTermIR this divisor
        dividend[index] with
    | none =>
        simp only [hterm, henselDivmodVHCSingleRemainderTermIR]
        rw [henselDivmodVHCSingleLoopIR_snd_toList]
        simp [List.append_assoc, henselDivmodVHCSingleRemainderTermIR]
    | some term =>
        simp only [hterm, henselDivmodVHCSingleRemainderTermIR]
        exact henselDivmodVHCSingleLoopIR_snd_toList this (index + 1)
          (quotient.push term) remainder dividend divisor
  next hdone =>
    have hindex : dividend.toList.length ≤ index := by simpa using hdone
    rw [List.drop_eq_nil_iff.mpr hindex]
    simp
termination_by dividend.size - index
decreasing_by all_goals omega

theorem henselDivmodVHCSingleRemainderTermIR_some_degree_lt
    (this : DenseUPolyZp) (divisor term output : UMonomial × Zp)
    (hrun : henselDivmodVHCSingleRemainderTermIR this divisor term =
      some output) :
    output = term ∧ term.1.deg < divisor.1.deg := by
  unfold henselDivmodVHCSingleRemainderTermIR at hrun
  cases hquotient : StrictSquarefreeZp.pairVecDivSingleTermIR this divisor term
      with
  | some quotient => simp [hquotient] at hrun
  | none =>
      simp only [hquotient, Option.some.injEq] at hrun
      have hdegree : term.1.deg < divisor.1.deg := by
        unfold StrictSquarefreeZp.pairVecDivSingleTermIR at hquotient
        split at hquotient
        next hdivides => simp at hquotient
        next hnotDivides => omega
      exact ⟨hrun.symm, hdegree⟩

theorem listSum_henselDivmodVHCSingleTerms
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (divisor : UMonomial × Zp)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hdivisorReduced : CLPoly.Math.Zp.Reduced this._p.toNat divisor.2)
    (hdivisorNonzero : divisor.2.val ≠ 0) :
    ∀ terms : List (UMonomial × Zp),
      (∀ term ∈ terms, CLPoly.Math.Zp.Reduced this._p.toNat term.2) →
      (∀ term ∈ terms, term.2.val ≠ 0) →
      CLPoly.Math.listSum this._p.toNat
            (terms.filterMap
              (StrictSquarefreeZp.pairVecDivSingleTermIR this divisor)) *
          Polynomial.monomial divisor.1.deg
            (CLPoly.Math.Zp.toZMod this._p.toNat divisor.2) +
        CLPoly.Math.listSum this._p.toNat
          (terms.filterMap
            (henselDivmodVHCSingleRemainderTermIR this divisor)) =
        CLPoly.Math.listSum this._p.toNat terms := by
  intro terms hreduced hnonzero
  induction terms with
  | nil => simp [CLPoly.Math.listSum]
  | cons term rest ih =>
      have htermReduced := hreduced term List.mem_cons_self
      have htermNonzero := hnonzero term List.mem_cons_self
      have hrestReduced : ∀ item ∈ rest,
          CLPoly.Math.Zp.Reduced this._p.toNat item.2 := by
        intro item hitem
        exact hreduced item (List.mem_cons_of_mem term hitem)
      have hrestNonzero : ∀ item ∈ rest, item.2.val ≠ 0 := by
        intro item hitem
        exact hnonzero item (List.mem_cons_of_mem term hitem)
      have ih' := ih hrestReduced hrestNonzero
      unfold henselDivmodVHCSingleRemainderTermIR at ih' ⊢
      cases hrun : StrictSquarefreeZp.pairVecDivSingleTermIR this divisor term
          with
      | none =>
          simp only [List.filterMap_cons, hrun, CLPoly.Math.listSum_cons,
            zero_mul, zero_add]
          rw [CLPoly.Math.listSum_cons, CLPoly.Math.listSum_cons]
          rw [← ih']
          ac_rfl
      | some output =>
          have htermMul :=
            StrictSquarefreeZp.pairVecDivSingleTermIR_monomial_mul this divisor
              term output hcfg hdivisorReduced htermReduced hdivisorNonzero
              htermNonzero hrun
          simp only [List.filterMap_cons, hrun]
          rw [CLPoly.Math.listSum_cons, add_mul, htermMul,
            CLPoly.Math.listSum_cons]
          rw [← ih']
          ac_rfl

theorem henselDivmodVHCSingleLoopIR_equation
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (dividend : SparsePolyZp) (divisor : UMonomial × Zp)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hdivisorReduced : CLPoly.Math.Zp.Reduced this._p.toNat divisor.2)
    (hdivisorNonzero : divisor.2.val ≠ 0)
    (hdividendCanonical : CLPoly.Math.SparsePolyZp.Canonical this._p.toNat
      dividend) :
    let output := henselDivmodVHCSingleLoopIR this 0 #[] #[] dividend divisor
    CLPoly.Math.SparsePolyZp.toPoly this._p.toNat output.1 *
          Polynomial.monomial divisor.1.deg
            (CLPoly.Math.Zp.toZMod this._p.toNat divisor.2) +
        CLPoly.Math.SparsePolyZp.toPoly this._p.toNat output.2 =
      CLPoly.Math.SparsePolyZp.toPoly this._p.toNat dividend := by
  dsimp only
  unfold CLPoly.Math.SparsePolyZp.toPoly
  rw [henselDivmodVHCSingleLoopIR_fst,
    StrictSquarefreeZp.pairVecDivSingleLoopIR_toList,
    henselDivmodVHCSingleLoopIR_snd_toList]
  simp only [List.drop_zero, List.nil_append]
  exact listSum_henselDivmodVHCSingleTerms this divisor hcfg hdivisorReduced
    hdivisorNonzero dividend.toList hdividendCanonical.1
      hdividendCanonical.2.2

theorem henselDivmodVHCSingleLoopIR_remainderCanonical
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (dividend : SparsePolyZp) (divisor : UMonomial × Zp)
    (hdividendCanonical : CLPoly.Math.SparsePolyZp.Canonical this._p.toNat
      dividend) :
    CLPoly.Math.SparsePolyZp.Canonical this._p.toNat
      (henselDivmodVHCSingleLoopIR this 0 #[] #[] dividend divisor).2 := by
  rw [CLPoly.Math.SparsePolyZp.Canonical,
    CLPoly.Math.SparsePolyZp.WellFormed_arr,
    henselDivmodVHCSingleLoopIR_snd_toList]
  simp only [List.drop_zero, List.nil_append]
  refine ⟨?_, ?_, ?_⟩
  · intro output houtput
    rw [List.mem_filterMap] at houtput
    rcases houtput with ⟨term, hterm, hrun⟩
    have houtputEq :=
      (henselDivmodVHCSingleRemainderTermIR_some_degree_lt this divisor
        term output hrun).1
    rw [houtputEq]
    exact hdividendCanonical.1 term hterm
  · rw [List.isChain_iff_pairwise]
    apply List.Pairwise.filterMap
      (R := fun left right : UMonomial × Zp => left.1.deg > right.1.deg)
      (S := fun left right : UMonomial × Zp => left.1.deg > right.1.deg)
      (henselDivmodVHCSingleRemainderTermIR this divisor)
    · intro left right horder leftOut hleft rightOut hright
      rcases henselDivmodVHCSingleRemainderTermIR_some_degree_lt this divisor
          left leftOut hleft with ⟨rfl, hleftDegree⟩
      rcases henselDivmodVHCSingleRemainderTermIR_some_degree_lt this divisor
          right rightOut hright with ⟨rfl, hrightDegree⟩
      exact horder
    · exact List.isChain_iff_pairwise.mp hdividendCanonical.2.1
  · intro output houtput
    rw [List.mem_filterMap] at houtput
    rcases houtput with ⟨term, hterm, hrun⟩
    have houtputEq :=
      (henselDivmodVHCSingleRemainderTermIR_some_degree_lt this divisor
        term output hrun).1
    rw [houtputEq]
    exact hdividendCanonical.2.2 term hterm

def henselDivmodVHCSingleBranchIR (this : DenseUPolyZp)
    (dividend divisor : SparsePolyZp) :
    RawExec (SparsePolyZp × SparsePolyZp) :=
  if hone : divisor.size = 1 then
    .ok (henselDivmodVHCSingleLoopIR this 0 #[] #[] dividend divisor[0])
  else
    .error .assertionFailure

/-- Complete non-aliasing five-argument `pair_vec_div` entry in source branch
order: reject zero divisor, clear both outputs, return for zero dividend,
execute the one-term loop, otherwise enter the well-founded VHC path. -/
def henselDivmodVHCIR (this : DenseUPolyZp)
    (dividend divisor : SparsePolyZp) :
    RawExec (SparsePolyZp × SparsePolyZp) :=
  if hdivisorEmpty : divisor.size = 0 then
    .error .assertionFailure
  else if hdividendEmpty : dividend.size = 0 then
    .ok (#[], #[])
  else if hsingle : divisor.size = 1 then
    henselDivmodVHCSingleBranchIR this dividend divisor
  else
    henselDivmodVHCGeneralBranchIR this dividend divisor

/-- Totality of the complete concrete Hensel EEA division operation on a
canonical nonzero divisor. -/
theorem henselDivmodVHCIR_succeeds
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (dividend divisor : SparsePolyZp)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hdividendCanonical : CLPoly.Math.SparsePolyZp.Canonical this._p.toNat
      dividend)
    (hdivisorCanonical : CLPoly.Math.SparsePolyZp.Canonical this._p.toNat
      divisor)
    (hdivisor : 0 < divisor.size) :
    ∃ output, henselDivmodVHCIR this dividend divisor = .ok output := by
  by_cases hdividend : dividend.size = 0
  · exact ⟨(#[], #[]), by
      simp [henselDivmodVHCIR, Nat.ne_of_gt hdivisor, hdividend]⟩
  · by_cases hsingle : divisor.size = 1
    · exact ⟨henselDivmodVHCSingleLoopIR this 0 #[] #[] dividend
          divisor[0], by
        simp [henselDivmodVHCIR, Nat.ne_of_gt hdivisor, hdividend, hsingle,
          henselDivmodVHCSingleBranchIR]⟩
    · have hgeneral : 1 < divisor.size := by omega
      rcases henselDivmodVHCGeneralBranchIR_succeeds this dividend divisor hcfg
          hdividendCanonical hdivisorCanonical (Nat.pos_of_ne_zero hdividend)
          hgeneral with ⟨output, houtput⟩
      exact ⟨output, by
        simpa [henselDivmodVHCIR, Nat.ne_of_gt hdivisor, hdividend, hsingle]
          using houtput⟩

def strictHenselEEARawOps (this : DenseUPolyZp) :
    Generated.StrictHensel.HenselEEARawOps where
  divmod := henselDivmodVHCIR this

/-- The remainder component of a successful iteration is exactly the source
push decision at the selected frontier degree. -/
theorem henselDivmodVHCIteration_remainder
    (this : DenseUPolyZp) (state next : HenselDivmodVHCState)
    (frontier : StrictSquarefreeZp.PairVecDivVHCFrontier)
    (dividend divisor : SparsePolyZp) (hdivisor : 0 < divisor.size)
    (hrun : henselDivmodVHCIteration this state frontier dividend divisor
      hdivisor = .ok next) :
    ∃ consumed,
      StrictSquarefreeZp.pairVecDivVHCConsumeEqualDegree this frontier.degree
          state.heap frontier.coefficient state.nodes #[] state.resetH
          state.quotient divisor = .ok consumed ∧
      next.remainder =
        if consumed.coefficient ≠ 0 ∧
            frontier.degree < divisor[0].1.deg then
          state.remainder.push
            (UMonomial.mk frontier.degree,
              ⟨consumed.coefficient, this._p⟩)
        else state.remainder := by
  unfold henselDivmodVHCIteration at hrun
  cases hconsume : StrictSquarefreeZp.pairVecDivVHCConsumeEqualDegree this
      frontier.degree state.heap frontier.coefficient state.nodes #[]
      state.resetH state.quotient divisor with
  | error fault => simp [hconsume, bind, Except.bind] at hrun
  | ok consumed =>
      simp only [hconsume, bind, Except.bind] at hrun
      cases hemit : StrictSquarefreeZp.pairVecDivVHCEmit this frontier consumed
          state.quotient divisor hdivisor with
      | error fault => simp [hemit, bind, Except.bind] at hrun
      | ok emitted =>
          simp only [hemit, bind, Except.bind] at hrun
          cases hreinsert : StrictSquarefreeZp.pairVecDivVHCReinsertLin
              emitted.2.1.heap emitted.2.1.nodes consumed.lin with
          | error fault => simp [hreinsert, bind, Except.bind] at hrun
          | ok reinserted =>
              simp only [hreinsert, bind, Except.bind] at hrun
              have hnext :
                  HenselDivmodVHCState.mk frontier.dividendIndex
                    reinserted.heap reinserted.nodes emitted.1
                    (if consumed.coefficient ≠ 0 ∧
                        frontier.degree < divisor[0].1.deg then
                      state.remainder.push
                        (UMonomial.mk frontier.degree,
                          ⟨consumed.coefficient, this._p⟩)
                    else state.remainder)
                    emitted.2.2 = next := Except.ok.inj hrun
              subst next
              exact ⟨consumed, rfl, rfl⟩

/-- A successful quotient-and-remainder VHC run has exactly the quotient of
the already verified quotient-only source loop.  The proof follows the actual
well-founded calls and erases only the separately accumulated remainder. -/
theorem henselDivmodVHCOuterLoop_projects_quotient
    (this : DenseUPolyZp) (dividend divisor : SparsePolyZp)
    (hdivisor : 0 < divisor.size) (limit : Nat)
    (state : HenselDivmodVHCState) (output : SparsePolyZp × SparsePolyZp)
    (hrun : henselDivmodVHCOuterLoop this limit state dividend divisor
      hdivisor = .ok output) :
    StrictSquarefreeZp.pairVecDivVHCOuterLoop this limit
        state.dividendIndex state.heap state.nodes state.quotient dividend
        divisor state.resetH = .ok output.1 := by
  induction limit using Nat.strong_induction_on generalizing state output with
  | h limit ih =>
      rw [henselDivmodVHCOuterLoop] at hrun
      rw [StrictSquarefreeZp.pairVecDivVHCOuterLoop]
      by_cases hdone : dividend.size ≤ state.dividendIndex ∧
          state.heap.size = 0
      · rw [dif_pos hdone] at hrun ⊢
        exact congrArg Except.ok (congrArg Prod.fst (Except.ok.inj hrun))
      · rw [dif_neg hdone] at hrun ⊢
        cases hselect : StrictSquarefreeZp.pairVecDivVHCSelectFrontier
            state.dividendIndex dividend state.heap state.nodes with
        | error fault => simp [hselect] at hrun
        | ok frontier =>
            rw [hselect] at hrun
            simp only at hrun ⊢
            by_cases hdecrease : frontier.degree < limit
            · rw [dif_pos hdecrease] at hrun ⊢
              cases hiteration : henselDivmodVHCIteration this state frontier
                  dividend divisor hdivisor with
              | error fault => simp [hiteration] at hrun
              | ok next =>
                  rw [hiteration] at hrun
                  simp only at hrun
                  rcases henselDivmodVHCIteration_projects_quotient this state
                      next frontier dividend divisor hdivisor hselect hiteration
                      with ⟨projected, hprojected, hindex, hheap, hnodes,
                        hquotient, hreset⟩
                  rw [hprojected]
                  simp only
                  rw [hindex, hheap, hnodes, hquotient, hreset]
                  exact ih frontier.degree hdecrease next output hrun
            · rw [dif_neg hdecrease] at hrun
              contradiction

/-- The quotient produced by the five-argument loop inherits the proved
high-coefficient product semantics of the same concrete VHC execution.  No
quotient is supplied to the program: it is projected from the successful
quotient-and-remainder run. -/
theorem henselDivmodVHCOuterLoop_productAgreesAbove_lead_of_success
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (globalLimit degreeLimit : Nat) (state : HenselDivmodVHCState)
    (dividend divisor : SparsePolyZp) (output : SparsePolyZp × SparsePolyZp)
    (owners : Nat → Finset Nat) (hdivisor : 0 < divisor.size)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hcanonical : CLPoly.Math.SparsePolyZp.Canonical this._p.toNat
      state.quotient)
    (hdividendCanonical : CLPoly.Math.SparsePolyZp.Canonical this._p.toNat
      dividend)
    (hdivisorCanonical : CLPoly.Math.SparsePolyZp.Canonical this._p.toNat
      divisor)
    (hquotientReady : ∀ frontier :
        StrictSquarefreeZp.PairVecDivVHCFrontier,
      StrictSquarefreeZp.pairVecDivVHCSelectFrontier state.dividendIndex
          dividend state.heap state.nodes = .ok frontier →
      StrictSquarefreeZp.PairVecDivVHCQuotientAbove frontier.degree
        divisor[0].1.deg state.quotient)
    (hremaining : StrictSquarefreeZp.PairVecDivVHCRemainingDividendBelow
      globalLimit state.dividendIndex dividend)
    (hbelow : StrictSquarefreeZp.PairVecDivVHCAllActiveNodesBelow globalLimit
      state.nodes)
    (hdenotes : ∀ (i : Nat)
        (node : StrictSquarefreeZp.PairVecDivVHCNode),
      state.nodes[i]? = some node → node.mono ≠ none →
        StrictSquarefreeZp.PairVecDivVHCNodeDenotes state.quotient divisor node)
    (hfixed : StrictSquarefreeZp.PairVecDivVHCNodeDivisorIndicesFixed
      state.nodes)
    (hready : StrictSquarefreeZp.PairVecDivVHCResetReady state.resetH
      state.quotient.size state.nodes)
    (hownership : StrictSquarefreeZp.PairVecDivVHCHeapChainOwnership
      state.heap owners state.nodes)
    (hcovered : StrictSquarefreeZp.PairVecDivVHCStateCovered state.heap
      state.nodes #[] state.resetH)
    (hsize : state.nodes.size = divisor.size - 1)
    (hhomogeneous : StrictSquarefreeZp.PairVecDivVHCHeapChainsHomogeneous
      state.heap owners state.nodes)
    (hordered : StrictSquarefreeZp.PairVecDivVHCHeapOrdered state.heap
      state.nodes)
    (hprefix : StrictSquarefreeZp.PairVecDivVHCCursorPrefixAbove degreeLimit
      state.nodes state.quotient divisor)
    (hprocessed : StrictSquarefreeZp.PairVecDivVHCQuotientLeadAbove degreeLimit
      divisor[0].1.deg state.quotient)
    (hconsumed : StrictSquarefreeZp.PairVecDivVHCConsumedDividendAbove
      degreeLimit state.dividendIndex dividend)
    (hagrees : StrictSquarefreeZp.PairVecDivVHCProductAgreesAbove
      this._p.toNat degreeLimit state.quotient dividend divisor)
    (hrun : henselDivmodVHCOuterLoop this degreeLimit state dividend divisor
      hdivisor = .ok output) :
    StrictSquarefreeZp.PairVecDivVHCProductAgreesAbove this._p.toNat
      divisor[0].1.deg output.1 dividend divisor := by
  have hquotientRun := henselDivmodVHCOuterLoop_projects_quotient this dividend
    divisor hdivisor degreeLimit state output hrun
  exact StrictSquarefreeZp.pairVecDivVHCOuterLoop_productAgreesAbove_lead_of_success
    this globalLimit degreeLimit state.dividendIndex state.heap state.nodes
    state.quotient dividend divisor output.1 state.resetH owners hdivisor hcfg
    hcanonical hdividendCanonical hdivisorCanonical hquotientReady hremaining
    hbelow hdenotes hfixed hready hownership hcovered hsize hhomogeneous
    hordered hprefix hprocessed hconsumed hagrees hquotientRun

/-- Concrete reachable prefixes of the double-output VHC loop. -/
inductive HenselDivmodVHCPrefix (this : DenseUPolyZp)
    (dividend divisor : SparsePolyZp) (hdivisor : 0 < divisor.size)
    (initialLimit : Nat) (initial : HenselDivmodVHCState) :
    Nat → HenselDivmodVHCState → Prop
  | refl : HenselDivmodVHCPrefix this dividend divisor hdivisor initialLimit
      initial initialLimit initial
  | step (limit : Nat) (state next : HenselDivmodVHCState)
      (frontier : StrictSquarefreeZp.PairVecDivVHCFrontier)
      (hprefix : HenselDivmodVHCPrefix this dividend divisor hdivisor
        initialLimit initial limit state)
      (hnotDone : ¬ (dividend.size ≤ state.dividendIndex ∧
        state.heap.size = 0))
      (hselect : StrictSquarefreeZp.pairVecDivVHCSelectFrontier
        state.dividendIndex dividend state.heap state.nodes = .ok frontier)
      (hdecrease : frontier.degree < limit)
      (hiteration : henselDivmodVHCIteration this state frontier dividend
        divisor hdivisor = .ok next) :
      HenselDivmodVHCPrefix this dividend divisor hdivisor initialLimit initial
        frontier.degree next

/-- Raw safety for every state reached by successful concrete VHC bodies. -/
structure HenselDivmodVHCExecutionInvariant (this : DenseUPolyZp)
    (dividend divisor : SparsePolyZp) (hdivisor : 0 < divisor.size)
    (initialLimit : Nat) (initial : HenselDivmodVHCState) : Prop where
  frontierReady : ∀ limit state,
    HenselDivmodVHCPrefix this dividend divisor hdivisor initialLimit initial
      limit state →
    ¬ (dividend.size ≤ state.dividendIndex ∧ state.heap.size = 0) →
    ∃ frontier,
      StrictSquarefreeZp.pairVecDivVHCSelectFrontier state.dividendIndex
        dividend state.heap state.nodes = .ok frontier ∧
      frontier.degree < limit
  iterationReady : ∀ limit state frontier,
    HenselDivmodVHCPrefix this dividend divisor hdivisor initialLimit initial
      limit state →
    StrictSquarefreeZp.pairVecDivVHCSelectFrontier state.dividendIndex
        dividend state.heap state.nodes = .ok frontier →
    frontier.degree < limit →
    ∃ next, henselDivmodVHCIteration this state frontier dividend divisor
      hdivisor = .ok next

def HenselDivmodVHCRemainderBelow (leadDegree : Nat)
    (remainder : SparsePolyZp) : Prop :=
  ∀ term ∈ remainder.toList, term.1.deg < leadDegree

/-- Remainder terms already emitted by a descending VHC prefix are at or
above the next strict frontier bound. -/
def HenselDivmodVHCRemainderAbove (degreeLimit : Nat)
    (remainder : SparsePolyZp) : Prop :=
  ∀ term ∈ remainder.toList, degreeLimit ≤ term.1.deg

/-- Coefficient agreement for the processed low-degree interval. -/
def HenselDivmodVHCRemainderAgreesAbove (p degreeLimit leadDegree : Nat)
    (remainder quotient dividend divisor : SparsePolyZp) : Prop :=
  ∀ degree, degreeLimit ≤ degree → degree < leadDegree →
    (CLPoly.Math.SparsePolyZp.toPoly p remainder).coeff degree =
      (CLPoly.Math.SparsePolyZp.toPoly p dividend -
        CLPoly.Math.SparsePolyZp.toPoly p quotient *
          CLPoly.Math.SparsePolyZp.toPoly p divisor).coeff degree

theorem HenselDivmodVHCRemainderBelow.empty (leadDegree : Nat) :
    HenselDivmodVHCRemainderBelow leadDegree #[] := by
  intro term hterm
  simp at hterm

theorem henselDivmodVHCSingleLoopIR_remainderBelow
    (this : DenseUPolyZp) (dividend : SparsePolyZp)
    (divisor : UMonomial × Zp) :
    HenselDivmodVHCRemainderBelow divisor.1.deg
      (henselDivmodVHCSingleLoopIR this 0 #[] #[] dividend divisor).2 := by
  intro term hterm
  rw [henselDivmodVHCSingleLoopIR_snd_toList] at hterm
  simp only [List.drop_zero, List.nil_append, List.mem_filterMap] at hterm
  rcases hterm with ⟨source, hsource, hrun⟩
  rcases henselDivmodVHCSingleRemainderTermIR_some_degree_lt this divisor
      source term hrun with ⟨rfl, hdegree⟩
  exact hdegree

/-- Full semantic contract of the concrete one-term-divisor source branch. -/
theorem henselDivmodVHCSingleBranchIR_refines
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (dividend divisor : SparsePolyZp)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hdividendCanonical : CLPoly.Math.SparsePolyZp.Canonical this._p.toNat
      dividend)
    (hdivisorCanonical : CLPoly.Math.SparsePolyZp.Canonical this._p.toNat
      divisor)
    (hdivisorSize : divisor.size = 1) :
    ∃ output,
      henselDivmodVHCSingleBranchIR this dividend divisor = .ok output ∧
      CLPoly.Math.SparsePolyZp.Canonical this._p.toNat output.1 ∧
      CLPoly.Math.SparsePolyZp.Canonical this._p.toNat output.2 ∧
      HenselDivmodVHCRemainderBelow divisor[0].1.deg output.2 ∧
      CLPoly.Math.SparsePolyZp.toPoly this._p.toNat output.1 *
            CLPoly.Math.SparsePolyZp.toPoly this._p.toNat divisor +
          CLPoly.Math.SparsePolyZp.toPoly this._p.toNat output.2 =
        CLPoly.Math.SparsePolyZp.toPoly this._p.toNat dividend := by
  rcases Array.size_eq_one_iff.mp hdivisorSize with ⟨divisorTerm, rfl⟩
  have hdivisorMem : divisorTerm ∈
      (#[divisorTerm] : SparsePolyZp).toList := by simp
  have hdivisorReduced := hdivisorCanonical.1 divisorTerm hdivisorMem
  have hdivisorNonzero := hdivisorCanonical.2.2 divisorTerm hdivisorMem
  let output := henselDivmodVHCSingleLoopIR this 0 #[] #[] dividend
    divisorTerm
  refine ⟨output, by simp [henselDivmodVHCSingleBranchIR, output], ?_, ?_,
    ?_, ?_⟩
  · dsimp only [output]
    rw [henselDivmodVHCSingleLoopIR_fst]
    exact StrictSquarefreeZp.pairVecDivSingleLoopIR_zero_canonical this
      dividend divisorTerm hcfg hdividendCanonical hdivisorReduced
      hdivisorNonzero
  · exact henselDivmodVHCSingleLoopIR_remainderCanonical this dividend
      divisorTerm hdividendCanonical
  · exact henselDivmodVHCSingleLoopIR_remainderBelow this dividend
      divisorTerm
  · have hequation := henselDivmodVHCSingleLoopIR_equation this dividend
      divisorTerm hcfg hdivisorReduced hdivisorNonzero hdividendCanonical
    simpa [output, CLPoly.Math.SparsePolyZp.toPoly, CLPoly.Math.listSum]
      using hequation

theorem HenselDivmodVHCRemainderBelow.coeff_eq_zero_at_or_above
    (p leadDegree degree : Nat) (remainder : SparsePolyZp)
    (hbelow : HenselDivmodVHCRemainderBelow leadDegree remainder)
    (hdegree : leadDegree ≤ degree) :
    (CLPoly.Math.SparsePolyZp.toPoly p remainder).coeff degree = 0 := by
  unfold CLPoly.Math.SparsePolyZp.toPoly
  apply StrictSquarefreeZp.listSum_coeff_zero_of_degree_absent
  intro term hterm hequal
  have := hbelow term hterm
  omega

theorem HenselDivmodVHCRemainderAbove.empty (degreeLimit : Nat) :
    HenselDivmodVHCRemainderAbove degreeLimit #[] := by
  intro term hterm
  simp at hterm

theorem HenselDivmodVHCRemainderAbove.push
    (degreeLimit frontierDegree : Nat) (remainder : SparsePolyZp)
    (coefficient : Zp)
    (habove : HenselDivmodVHCRemainderAbove degreeLimit remainder)
    (hdecrease : frontierDegree < degreeLimit) :
    HenselDivmodVHCRemainderAbove frontierDegree
      (remainder.push (UMonomial.mk frontierDegree, coefficient)) := by
  intro term hterm
  simp only [Array.toList_push, List.mem_append, List.mem_singleton] at hterm
  rcases hterm with hold | hnew
  · exact Nat.le_trans (Nat.le_of_lt hdecrease) (habove term hold)
  · subst term
    exact Nat.le_refl _

theorem HenselDivmodVHCRemainderAbove.coeff_eq_zero_below
    (p degreeLimit degree : Nat) (remainder : SparsePolyZp)
    (habove : HenselDivmodVHCRemainderAbove degreeLimit remainder)
    (hdegree : degree < degreeLimit) :
    (CLPoly.Math.SparsePolyZp.toPoly p remainder).coeff degree = 0 := by
  unfold CLPoly.Math.SparsePolyZp.toPoly
  apply StrictSquarefreeZp.listSum_coeff_zero_of_degree_absent
  intro term hterm hequal
  have := habove term hterm
  omega

/-- A degree skipped by the concrete VHC frontier selector is also absent
from the actual division residual.  Both zero coefficients come from the
generated dividend cursor and heap/product invariants; no division
specification is consulted. -/
theorem henselDivmodVHCResidual_coeff_eq_zero_of_gap
    (p degreeLimit targetDegree dividendIndex resetH : Nat)
    (dividend quotient divisor : SparsePolyZp)
    (heap : Array Nat)
    (nodes : Array StrictSquarefreeZp.PairVecDivVHCNode)
    (owners : Nat → Finset Nat)
    (frontier : StrictSquarefreeZp.PairVecDivVHCFrontier)
    (hdivisor : 0 < divisor.size)
    (hsize : nodes.size = divisor.size - 1)
    (hfixed : StrictSquarefreeZp.PairVecDivVHCNodeDivisorIndicesFixed nodes)
    (hstate : StrictSquarefreeZp.PairVecDivVHCStateCovered heap nodes #[]
      resetH)
    (hownership : StrictSquarefreeZp.PairVecDivVHCHeapChainOwnership heap
      owners nodes)
    (hhomogeneous : StrictSquarefreeZp.PairVecDivVHCHeapChainsHomogeneous heap
      owners nodes)
    (hresetReady : StrictSquarefreeZp.PairVecDivVHCResetReady resetH
      quotient.size nodes)
    (hordered : StrictSquarefreeZp.PairVecDivVHCHeapOrdered heap nodes)
    (hdenotes : ∀ (i : Nat)
      (node : StrictSquarefreeZp.PairVecDivVHCNode),
      nodes[i]? = some node → node.mono ≠ none →
        StrictSquarefreeZp.PairVecDivVHCNodeDenotes quotient divisor node)
    (hdividendCanonical : CLPoly.Math.SparsePolyZp.Canonical p dividend)
    (hconsumed : StrictSquarefreeZp.PairVecDivVHCConsumedDividendAbove
      degreeLimit dividendIndex dividend)
    (hquotientCanonical : CLPoly.Math.SparsePolyZp.Canonical p quotient)
    (hprefix : StrictSquarefreeZp.PairVecDivVHCCursorPrefixAbove degreeLimit
      nodes quotient divisor)
    (hprocessed : StrictSquarefreeZp.PairVecDivVHCQuotientLeadAbove
      degreeLimit divisor[0].1.deg quotient)
    (hfrontier : frontier.degree < targetDegree)
    (htarget : targetDegree < degreeLimit)
    (hselect : StrictSquarefreeZp.pairVecDivVHCSelectFrontier dividendIndex
      dividend heap nodes = .ok frontier) :
    (CLPoly.Math.SparsePolyZp.toPoly p dividend -
      CLPoly.Math.SparsePolyZp.toPoly p quotient *
        CLPoly.Math.SparsePolyZp.toPoly p divisor).coeff targetDegree = 0 := by
  rw [Polynomial.coeff_sub,
    StrictSquarefreeZp.pairVecDivVHCDividend_coeff_eq_zero_of_gap p
      degreeLimit targetDegree dividendIndex dividend heap nodes frontier
      hdividendCanonical hconsumed hfrontier htarget hselect,
    StrictSquarefreeZp.pairVecDivVHCProduct_coeff_eq_zero_of_gap p
      degreeLimit targetDegree dividendIndex resetH dividend quotient divisor
      heap nodes owners frontier hdivisor hsize hfixed hstate hownership
      hhomogeneous hresetReady hordered hdenotes hquotientCanonical hprefix
      hprocessed hfrontier htarget hselect]
  simp

/-- At a concrete terminal VHC state the residual vanishes at every degree
below the current well-founded bound.  This combines exhaustion of the real
dividend cursor with exhaustion of the real product heap. -/
theorem henselDivmodVHCResidual_coeff_eq_zero_of_done
    (p degreeLimit targetDegree dividendIndex resetH : Nat)
    (dividend quotient divisor : SparsePolyZp)
    (nodes : Array StrictSquarefreeZp.PairVecDivVHCNode)
    (owners : Nat → Finset Nat)
    (hdivisor : 0 < divisor.size)
    (hsize : nodes.size = divisor.size - 1)
    (hfixed : StrictSquarefreeZp.PairVecDivVHCNodeDivisorIndicesFixed nodes)
    (hstate : StrictSquarefreeZp.PairVecDivVHCStateCovered #[] nodes #[]
      resetH)
    (hownership : StrictSquarefreeZp.PairVecDivVHCHeapChainOwnership #[]
      owners nodes)
    (hresetReady : StrictSquarefreeZp.PairVecDivVHCResetReady resetH
      quotient.size nodes)
    (hprefix : StrictSquarefreeZp.PairVecDivVHCCursorPrefixAbove degreeLimit
      nodes quotient divisor)
    (hprocessed : StrictSquarefreeZp.PairVecDivVHCQuotientLeadAbove
      degreeLimit divisor[0].1.deg quotient)
    (hdividendCanonical : CLPoly.Math.SparsePolyZp.Canonical p dividend)
    (hconsumed : StrictSquarefreeZp.PairVecDivVHCConsumedDividendAbove
      degreeLimit dividendIndex dividend)
    (hdone : dividend.size ≤ dividendIndex)
    (htarget : targetDegree < degreeLimit) :
    (CLPoly.Math.SparsePolyZp.toPoly p dividend -
      CLPoly.Math.SparsePolyZp.toPoly p quotient *
        CLPoly.Math.SparsePolyZp.toPoly p divisor).coeff targetDegree = 0 := by
  rw [Polynomial.coeff_sub,
    StrictSquarefreeZp.pairVecDivVHCDividend_coeff_eq_zero_of_done p
      degreeLimit targetDegree dividendIndex dividend hdividendCanonical
      hconsumed hdone htarget,
    StrictSquarefreeZp.pairVecDivVHCProduct_coeff_eq_zero_of_done p
      degreeLimit targetDegree resetH quotient divisor nodes owners hdivisor
      hsize hfixed hstate hownership hresetReady hprefix hprocessed htarget]
  simp

theorem henselDivmodRemainderPush_coeff
    (p degree : Nat) (remainder : SparsePolyZp) (coefficient : Zp)
    (habsent :
      (CLPoly.Math.SparsePolyZp.toPoly p remainder).coeff degree = 0) :
    (CLPoly.Math.SparsePolyZp.toPoly p
      (remainder.push (UMonomial.mk degree, coefficient))).coeff degree =
        CLPoly.Math.Zp.toZMod p coefficient := by
  unfold CLPoly.Math.SparsePolyZp.toPoly at habsent ⊢
  rw [Array.toList_push, CLPoly.Math.listSum_append,
    Polynomial.coeff_add, habsent]
  simp [CLPoly.Math.listSum, CLPoly.Math.Zp.toZMod]

theorem henselDivmodRemainderPush_coeff_ne
    (p degree frontierDegree : Nat) (remainder : SparsePolyZp)
    (coefficient : Zp) (hne : degree ≠ frontierDegree) :
    (CLPoly.Math.SparsePolyZp.toPoly p
      (remainder.push (UMonomial.mk frontierDegree, coefficient))).coeff
        degree =
      (CLPoly.Math.SparsePolyZp.toPoly p remainder).coeff degree := by
  unfold CLPoly.Math.SparsePolyZp.toPoly
  rw [Array.toList_push, CLPoly.Math.listSum_append,
    Polynomial.coeff_add]
  simp [CLPoly.Math.listSum, Polynomial.coeff_monomial, Ne.symm hne]

/-- Pure interval-extension step used by the below-leading-degree phase.  It
combines the exact source push decision with a selected residual coefficient
and the zero residual gap above the selected frontier. -/
theorem henselDivmodRemainderStep_preserves
    (p degreeLimit frontierDegree leadDegree : Nat)
    (remainder nextRemainder quotient dividend divisor : SparsePolyZp)
    (coefficient : UInt64) (prime : UInt64)
    (hdecrease : frontierDegree < degreeLimit)
    (hbelowLead : frontierDegree < leadDegree)
    (habove : HenselDivmodVHCRemainderAbove degreeLimit remainder)
    (hagrees : HenselDivmodVHCRemainderAgreesAbove p degreeLimit leadDegree
      remainder quotient dividend divisor)
    (hfrontierResidual : (coefficient.toNat : ZMod p) =
      (CLPoly.Math.SparsePolyZp.toPoly p dividend -
        CLPoly.Math.SparsePolyZp.toPoly p quotient *
          CLPoly.Math.SparsePolyZp.toPoly p divisor).coeff frontierDegree)
    (hgap : ∀ degree, frontierDegree < degree → degree < degreeLimit →
      (CLPoly.Math.SparsePolyZp.toPoly p dividend -
        CLPoly.Math.SparsePolyZp.toPoly p quotient *
          CLPoly.Math.SparsePolyZp.toPoly p divisor).coeff degree = 0)
    (hnext : nextRemainder =
      if coefficient ≠ 0 then
        remainder.push
          (UMonomial.mk frontierDegree, ⟨coefficient, prime⟩)
      else remainder) :
    HenselDivmodVHCRemainderAbove frontierDegree nextRemainder ∧
      HenselDivmodVHCRemainderAgreesAbove p frontierDegree leadDegree
        nextRemainder quotient dividend divisor := by
  by_cases hcoefficient : coefficient ≠ 0
  · rw [hnext, if_pos hcoefficient]
    refine ⟨HenselDivmodVHCRemainderAbove.push degreeLimit frontierDegree
      remainder ⟨coefficient, prime⟩ habove hdecrease, ?_⟩
    intro degree hfrontier hlead
    by_cases hequal : degree = frontierDegree
    · subst degree
      rw [henselDivmodRemainderPush_coeff p frontierDegree remainder
        ⟨coefficient, prime⟩
        (HenselDivmodVHCRemainderAbove.coeff_eq_zero_below p degreeLimit
          frontierDegree remainder habove hdecrease)]
      exact hfrontierResidual
    · rw [henselDivmodRemainderPush_coeff_ne p degree frontierDegree remainder
        ⟨coefficient, prime⟩ hequal]
      by_cases hold : degreeLimit ≤ degree
      · exact hagrees degree hold hlead
      · rw [HenselDivmodVHCRemainderAbove.coeff_eq_zero_below p degreeLimit
            degree remainder habove (by omega),
          hgap degree (by omega) (by omega)]
  · have hzero : coefficient = 0 := by
      simpa only [not_ne_iff] using hcoefficient
    subst coefficient
    rw [hnext, if_neg (by simp)]
    refine ⟨?_, ?_⟩
    · intro term hterm
      exact Nat.le_trans (Nat.le_of_lt hdecrease) (habove term hterm)
    · intro degree hfrontier hlead
      by_cases hequal : degree = frontierDegree
      · subst degree
        rw [HenselDivmodVHCRemainderAbove.coeff_eq_zero_below p degreeLimit
          frontierDegree remainder habove hdecrease]
        simpa using hfrontierResidual
      · by_cases hold : degreeLimit ≤ degree
        · exact hagrees degree hold hlead
        · rw [HenselDivmodVHCRemainderAbove.coeff_eq_zero_below p degreeLimit
              degree remainder habove (by omega),
            hgap degree (by omega) (by omega)]

/-- Exact trace of the double-output VHC source loop. -/
inductive HenselDivmodVHCCorrect (this : DenseUPolyZp)
    (dividend divisor : SparsePolyZp) (hdivisor : 0 < divisor.size) :
    Nat → HenselDivmodVHCState → SparsePolyZp × SparsePolyZp → Prop
  | done (limit : Nat) (state : HenselDivmodVHCState)
      (hdone : dividend.size ≤ state.dividendIndex ∧ state.heap.size = 0) :
      HenselDivmodVHCCorrect this dividend divisor hdivisor limit state
        (state.quotient, state.remainder)
  | step (limit : Nat) (state next : HenselDivmodVHCState)
      (frontier : StrictSquarefreeZp.PairVecDivVHCFrontier)
      (output : SparsePolyZp × SparsePolyZp)
      (hnotDone : ¬ (dividend.size ≤ state.dividendIndex ∧
        state.heap.size = 0))
      (hselect : StrictSquarefreeZp.pairVecDivVHCSelectFrontier
        state.dividendIndex dividend state.heap state.nodes = .ok frontier)
      (hdecrease : frontier.degree < limit)
      (hiteration : henselDivmodVHCIteration this state frontier dividend
        divisor hdivisor = .ok next)
      (htail : HenselDivmodVHCCorrect this dividend divisor hdivisor
        frontier.degree next output) :
      HenselDivmodVHCCorrect this dividend divisor hdivisor limit state output

theorem HenselDivmodVHCCorrect.remainderBelow
    (this : DenseUPolyZp) (dividend divisor : SparsePolyZp)
    (hdivisor : 0 < divisor.size) (limit : Nat)
    (state : HenselDivmodVHCState) (output : SparsePolyZp × SparsePolyZp)
    (hcorrect : HenselDivmodVHCCorrect this dividend divisor hdivisor limit
      state output)
    (hbelow : HenselDivmodVHCRemainderBelow divisor[0].1.deg
      state.remainder) :
    HenselDivmodVHCRemainderBelow divisor[0].1.deg output.2 := by
  induction hcorrect with
  | done limit state hdone => exact hbelow
  | step limit state next frontier output hnotDone hselect hdecrease
      hiteration htail ih =>
      rcases henselDivmodVHCIteration_remainder this state next frontier
          dividend divisor hdivisor hiteration with
        ⟨consumed, hconsume, hremainder⟩
      apply ih
      rw [hremainder]
      by_cases hpush : consumed.coefficient ≠ 0 ∧
          frontier.degree < divisor[0].1.deg
      · rw [if_pos hpush]
        intro term hterm
        simp only [Array.toList_push, List.mem_append,
          List.mem_singleton] at hterm
        rcases hterm with hold | hnew
        · exact hbelow term hold
        · subst term
          exact hpush.2
      · rw [if_neg hpush]
        exact hbelow

/-- No product with the divisor leading monomial can contribute below its
degree, independently of the current quotient contents. -/
theorem henselQuotient_mul_lead_coeff_eq_zero_below
    (p degree leadDegree : Nat) (quotient : SparsePolyZp)
    (leadCoefficient : Zp) (hbelow : degree < leadDegree) :
    (CLPoly.Math.SparsePolyZp.toPoly p quotient *
      Polynomial.monomial leadDegree
        (CLPoly.Math.Zp.toZMod p leadCoefficient)).coeff degree = 0 := by
  unfold CLPoly.Math.SparsePolyZp.toPoly
  rw [show Polynomial.monomial leadDegree
        (CLPoly.Math.Zp.toZMod p leadCoefficient) =
      CLPoly.Math.listSum p [(⟨leadDegree⟩, leadCoefficient)] by
    simp [CLPoly.Math.listSum, CLPoly.Math.Zp.toZMod]]
  rw [← StrictSquarefreeZp.pairVecDivVHCListProductCoeffValue_eq_coeff p
    degree quotient.toList [(⟨leadDegree⟩, leadCoefficient)]]
  induction quotient.toList with
  | nil => simp [StrictSquarefreeZp.pairVecDivVHCListProductCoeffValue]
  | cons term terms ih =>
      have hdegree : term.1.deg + leadDegree ≠ degree := by omega
      simp [StrictSquarefreeZp.pairVecDivVHCListProductCoeffValue, hdegree, ih]

/-- At a below-leading-degree frontier, the coefficient physically consumed
by the five-argument loop is the full polynomial residual coefficient, and
the next remainder is exactly the source's conditional push of that value. -/
theorem henselDivmodVHCIteration_below_lead_residual
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (degreeLimit : Nat) (state next : HenselDivmodVHCState)
    (frontier : StrictSquarefreeZp.PairVecDivVHCFrontier)
    (dividend divisor : SparsePolyZp)
    (owners : Nat → Finset Nat) (hdivisor : 0 < divisor.size)
    (hsize : state.nodes.size = divisor.size - 1)
    (hstate : StrictSquarefreeZp.PairVecDivVHCStateCovered state.heap
      state.nodes #[] state.resetH)
    (hownership : StrictSquarefreeZp.PairVecDivVHCHeapChainOwnership
      state.heap owners state.nodes)
    (hhomogeneous : StrictSquarefreeZp.PairVecDivVHCHeapChainsHomogeneous
      state.heap owners state.nodes)
    (hordered : StrictSquarefreeZp.PairVecDivVHCHeapOrdered state.heap
      state.nodes)
    (hdenotes : ∀ (i : Nat)
        (node : StrictSquarefreeZp.PairVecDivVHCNode),
      state.nodes[i]? = some node → node.mono ≠ none →
        StrictSquarefreeZp.PairVecDivVHCNodeDenotes state.quotient divisor node)
    (hfixed : StrictSquarefreeZp.PairVecDivVHCNodeDivisorIndicesFixed
      state.nodes)
    (hresetReady : StrictSquarefreeZp.PairVecDivVHCResetReady state.resetH
      state.quotient.size state.nodes)
    (hprefix : StrictSquarefreeZp.PairVecDivVHCCursorPrefixAbove degreeLimit
      state.nodes state.quotient divisor)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hquotientCanonical : CLPoly.Math.SparsePolyZp.Canonical this._p.toNat
      state.quotient)
    (hdividendCanonical : CLPoly.Math.SparsePolyZp.Canonical this._p.toNat
      dividend)
    (hconsumed : StrictSquarefreeZp.PairVecDivVHCConsumedDividendAbove
      degreeLimit state.dividendIndex dividend)
    (hdecrease : frontier.degree < degreeLimit)
    (hbelowLead : frontier.degree < divisor[0].1.deg)
    (hselect : StrictSquarefreeZp.pairVecDivVHCSelectFrontier
      state.dividendIndex dividend state.heap state.nodes = .ok frontier)
    (hrun : henselDivmodVHCIteration this state frontier dividend divisor
      hdivisor = .ok next) :
    ∃ consumed,
      StrictSquarefreeZp.pairVecDivVHCConsumeEqualDegree this frontier.degree
          state.heap frontier.coefficient state.nodes #[] state.resetH
          state.quotient divisor = .ok consumed ∧
      (consumed.coefficient.toNat : ZMod this._p.toNat) =
        (CLPoly.Math.SparsePolyZp.toPoly this._p.toNat dividend).coeff
            frontier.degree -
          (CLPoly.Math.SparsePolyZp.toPoly this._p.toNat state.quotient *
            CLPoly.Math.SparsePolyZp.toPoly this._p.toNat divisor).coeff
              frontier.degree ∧
      next.remainder =
        if consumed.coefficient ≠ 0 then
          state.remainder.push
            (UMonomial.mk frontier.degree,
              ⟨consumed.coefficient, this._p⟩)
        else state.remainder := by
  rcases henselDivmodVHCIteration_remainder this state next frontier dividend
      divisor hdivisor hrun with ⟨consumed, hconsume, hremainder⟩
  rcases henselDivmodVHCIteration_projects_quotient this state next frontier
      dividend divisor hdivisor hselect hrun with
    ⟨projected, hprojected, hindex, hheap, hnodes, hquotient, hreset⟩
  rcases StrictSquarefreeZp.pairVecDivVHCOuterIteration_residual_coefficient_toPoly
      this degreeLimit state.dividendIndex state.heap state.nodes
      state.quotient dividend divisor state.resetH frontier projected owners
      hdivisor hsize hstate hownership hhomogeneous hordered hdenotes hfixed
      hresetReady hprefix hcfg hquotientCanonical hdividendCanonical hconsumed
      hdecrease hselect hprojected with
    ⟨residualConsumed, hresidualConsume, hresidual⟩
  have hconsumedEq : residualConsumed = consumed := by
    rw [hconsume] at hresidualConsume
    exact (Except.ok.inj hresidualConsume).symm
  subst residualConsumed
  have hdivisorPoly :=
    StrictSquarefreeZp.sparsePolyZp_toPoly_eq_head_add_tail this._p.toNat
      divisor hdivisor
  have hleadZero := henselQuotient_mul_lead_coeff_eq_zero_below this._p.toNat
    frontier.degree divisor[0].1.deg state.quotient divisor[0].2 hbelowLead
  refine ⟨consumed, hconsume, ?_, ?_⟩
  · rw [hdivisorPoly, mul_add, Polynomial.coeff_add, hleadZero, zero_add]
    exact hresidual
  · rw [hremainder]
    by_cases hcoefficient : consumed.coefficient ≠ 0
    · simp [hcoefficient, hbelowLead]
    · simp [hcoefficient]

/-- Operational facts carried by a reachable state of the concrete
quotient-and-remainder VHC loop.  Every field describes the generated cursor,
heap, node chains, or quotient prefix; the structure contains no expected
result and no semantic division oracle. -/
structure HenselDivmodVHCOperationalInvariant (p degreeLimit : Nat)
    (state : HenselDivmodVHCState) (dividend divisor : SparsePolyZp)
    (hdivisor : 0 < divisor.size) : Prop where
  quotientCanonical : CLPoly.Math.SparsePolyZp.Canonical p state.quotient
  size : state.nodes.size = divisor.size - 1
  fixed : StrictSquarefreeZp.PairVecDivVHCNodeDivisorIndicesFixed state.nodes
  covered : StrictSquarefreeZp.PairVecDivVHCStateCovered state.heap state.nodes
    #[] state.resetH
  heapChains : ∃ owners : Nat → Finset Nat,
    StrictSquarefreeZp.PairVecDivVHCHeapChainOwnership state.heap owners
        state.nodes ∧
      StrictSquarefreeZp.PairVecDivVHCHeapChainsHomogeneous state.heap owners
        state.nodes
  ordered : StrictSquarefreeZp.PairVecDivVHCHeapOrdered state.heap state.nodes
  denotes : ∀ (i : Nat) (node : StrictSquarefreeZp.PairVecDivVHCNode),
    state.nodes[i]? = some node → node.mono ≠ none →
      StrictSquarefreeZp.PairVecDivVHCNodeDenotes state.quotient divisor node
  resetReady : StrictSquarefreeZp.PairVecDivVHCResetReady state.resetH
    state.quotient.size state.nodes
  cursorPrefix : StrictSquarefreeZp.PairVecDivVHCCursorPrefixAbove degreeLimit
    state.nodes state.quotient divisor
  quotientProcessed : StrictSquarefreeZp.PairVecDivVHCQuotientLeadAbove
    degreeLimit divisor[0].1.deg state.quotient
  dividendConsumed : StrictSquarefreeZp.PairVecDivVHCConsumedDividendAbove
    degreeLimit state.dividendIndex dividend

/-- The complete raw representation needed to transport the operational core
through another concrete body.  The extra facts are heap/frontier bounds, not
semantic predictions about the output. -/
structure HenselDivmodVHCFullOperationalInvariant
    (this : DenseUPolyZp) (globalLimit degreeLimit : Nat)
    (state : HenselDivmodVHCState) (dividend divisor : SparsePolyZp)
    (hdivisor : 0 < divisor.size) : Prop where
  core : HenselDivmodVHCOperationalInvariant this._p.toNat degreeLimit state
    dividend divisor hdivisor
  divisorCanonical : CLPoly.Math.SparsePolyZp.Canonical this._p.toNat divisor
  quotientReady : ∀ frontier : StrictSquarefreeZp.PairVecDivVHCFrontier,
    StrictSquarefreeZp.pairVecDivVHCSelectFrontier state.dividendIndex dividend
        state.heap state.nodes = .ok frontier →
      StrictSquarefreeZp.PairVecDivVHCQuotientAbove frontier.degree
        divisor[0].1.deg state.quotient
  remaining : StrictSquarefreeZp.PairVecDivVHCRemainingDividendBelow
    globalLimit state.dividendIndex dividend
  activeBelow : StrictSquarefreeZp.PairVecDivVHCAllActiveNodesBelow globalLimit
    state.nodes

theorem HenselDivmodVHCFullOperationalInvariant.step
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (globalLimit degreeLimit : Nat) (state next : HenselDivmodVHCState)
    (frontier : StrictSquarefreeZp.PairVecDivVHCFrontier)
    (dividend divisor : SparsePolyZp) (hdivisor : 0 < divisor.size)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hdividendCanonical : CLPoly.Math.SparsePolyZp.Canonical this._p.toNat
      dividend)
    (hinvariant : HenselDivmodVHCFullOperationalInvariant this globalLimit
      degreeLimit state dividend divisor hdivisor)
    (hselect : StrictSquarefreeZp.pairVecDivVHCSelectFrontier
      state.dividendIndex dividend state.heap state.nodes = .ok frontier)
    (hdecrease : frontier.degree < degreeLimit)
    (hiteration : henselDivmodVHCIteration this state frontier dividend divisor
      hdivisor = .ok next) :
    HenselDivmodVHCFullOperationalInvariant this globalLimit frontier.degree
      next dividend divisor hdivisor := by
  let projected := next.quotientState
  have hprojected : StrictSquarefreeZp.pairVecDivVHCOuterIteration this
      state.dividendIndex state.heap state.nodes state.quotient dividend divisor
      state.resetH = .ok projected := by
    exact henselDivmodVHCIteration_projects_quotientState this state next
      frontier dividend divisor hdivisor hselect hiteration
  rcases hinvariant.core.heapChains with
    ⟨owners, hownership, hhomogeneous⟩
  have hcoveredNodes := hinvariant.core.covered.covered_with state.heap
    state.nodes #[] state.resetH owners hownership
  have hnextCore :=
    StrictSquarefreeZp.pairVecDivVHCOuterIteration_preserves_heapChainsOwned
      this globalLimit state.dividendIndex state.heap state.nodes
      state.quotient dividend divisor state.resetH projected owners hdivisor
      hcfg hinvariant.core.quotientCanonical hdividendCanonical
      hinvariant.divisorCanonical hinvariant.quotientReady hinvariant.remaining
      hinvariant.activeBelow hinvariant.core.denotes hinvariant.core.fixed
      hinvariant.core.resetReady hownership hprojected
  rcases hnextCore with
    ⟨hnextCanonical, hnextOwned, hnextRemaining, hnextBelow, hnextDenotes,
      hnextFixed, hnextReady⟩
  have hnextState :=
    StrictSquarefreeZp.pairVecDivVHCOuterIteration_preserves_stateCovered this
      this._p.toNat globalLimit state.dividendIndex state.heap state.nodes
      state.quotient dividend divisor state.resetH frontier projected hdivisor
      hselect hinvariant.core.quotientCanonical hinvariant.activeBelow
      hinvariant.core.denotes hinvariant.core.fixed hinvariant.core.resetReady
      hinvariant.core.covered hprojected
  rcases
      StrictSquarefreeZp.pairVecDivVHCOuterIteration_preserves_heapChainsHomogeneous
        this this._p.toNat globalLimit state.dividendIndex state.heap state.nodes
        state.quotient dividend divisor state.resetH frontier projected owners
        hdivisor hinvariant.core.quotientCanonical hinvariant.activeBelow
        hinvariant.core.denotes hinvariant.core.fixed
        hinvariant.core.resetReady hownership hhomogeneous hselect hprojected
      with ⟨nextOwners, hnextOwnership, hnextHomogeneous⟩
  have hnextOrdered :=
    StrictSquarefreeZp.pairVecDivVHCOuterIteration_preserves_heapOrdered this
      this._p.toNat globalLimit state.dividendIndex state.heap state.nodes
      state.quotient dividend divisor state.resetH frontier projected owners
      hdivisor hinvariant.core.quotientCanonical hinvariant.activeBelow
      hinvariant.core.denotes hinvariant.core.fixed hinvariant.core.resetReady
      hownership hinvariant.core.ordered hselect hprojected
  have hnextStrict :=
    StrictSquarefreeZp.pairVecDivVHCOuterIteration_frontierBelow this
      this._p.toNat globalLimit state.dividendIndex state.heap state.nodes
      state.quotient dividend divisor state.resetH frontier projected owners
      hdivisor hselect hinvariant.core.quotientCanonical hdividendCanonical
      hinvariant.divisorCanonical hinvariant.activeBelow
      hinvariant.core.denotes hinvariant.core.fixed hinvariant.core.resetReady
      hownership hcoveredNodes hhomogeneous hinvariant.core.ordered hprojected
  have hnextQuotientReady : ∀ nextFrontier :
      StrictSquarefreeZp.PairVecDivVHCFrontier,
      StrictSquarefreeZp.pairVecDivVHCSelectFrontier next.dividendIndex dividend
          next.heap next.nodes = .ok nextFrontier →
        StrictSquarefreeZp.PairVecDivVHCQuotientAbove nextFrontier.degree
          divisor[0].1.deg next.quotient := by
    intro nextFrontier hnextSelect
    have hnextDegree :=
      StrictSquarefreeZp.pairVecDivVHCSelectFrontier_degree_lt frontier.degree
        next.dividendIndex dividend next.heap next.nodes nextFrontier (by
          simpa [projected, HenselDivmodVHCState.quotientState] using hnextStrict)
        hnextSelect
    have habove :=
      StrictSquarefreeZp.pairVecDivVHCOuterIteration_preserves_quotientAbove_of_lt
        this state.dividendIndex nextFrontier.degree state.heap state.nodes
        state.quotient dividend divisor state.resetH frontier projected
        hdivisor (hinvariant.quotientReady frontier hselect) hnextDegree hselect
        hprojected
    simpa [projected, HenselDivmodVHCState.quotientState] using habove
  have hnextPrefix :=
    StrictSquarefreeZp.pairVecDivVHCOuterIteration_preserves_cursorPrefixAbove
      this this._p.toNat globalLimit degreeLimit state.dividendIndex state.heap
      state.nodes state.quotient dividend divisor state.resetH frontier
      projected owners hdivisor hselect hdecrease
      hinvariant.core.quotientCanonical hinvariant.activeBelow
      hinvariant.core.denotes hinvariant.core.fixed hinvariant.core.resetReady
      hinvariant.core.covered hownership hhomogeneous
      hinvariant.core.cursorPrefix hprojected
  have hnextProcessed :=
    StrictSquarefreeZp.pairVecDivVHCOuterIteration_preserves_quotientLeadAbove
      this degreeLimit state.dividendIndex state.heap state.nodes
      state.quotient dividend divisor state.resetH frontier projected hdivisor
      hinvariant.core.quotientProcessed hdecrease hselect hprojected
  have hnextConsumed :=
    StrictSquarefreeZp.pairVecDivVHCOuterIteration_preserves_consumed_above this
      degreeLimit state.dividendIndex state.heap state.nodes state.quotient
      dividend divisor state.resetH frontier projected hdivisor
      hinvariant.core.dividendConsumed hdecrease hselect hprojected
  have hnextSize : next.nodes.size = divisor.size - 1 := by
    have := StrictSquarefreeZp.pairVecDivVHCOuterIteration_nodes_size this
      state.dividendIndex state.heap state.nodes state.quotient dividend divisor
      state.resetH frontier projected hdivisor hselect hprojected
    simpa [projected, HenselDivmodVHCState.quotientState,
      hinvariant.core.size] using this
  refine ⟨⟨?_, hnextSize, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩,
    hinvariant.divisorCanonical, hnextQuotientReady, ?_, ?_⟩
  · simpa [projected, HenselDivmodVHCState.quotientState] using
      hnextCanonical
  · simpa [projected, HenselDivmodVHCState.quotientState] using hnextFixed
  · simpa [projected, HenselDivmodVHCState.quotientState] using hnextState
  · exact ⟨nextOwners, by
      simpa [projected, HenselDivmodVHCState.quotientState] using
        hnextOwnership, by
      simpa [projected, HenselDivmodVHCState.quotientState] using
        hnextHomogeneous⟩
  · simpa [projected, HenselDivmodVHCState.quotientState] using hnextOrdered
  · simpa [projected, HenselDivmodVHCState.quotientState] using hnextDenotes
  · simpa [projected, HenselDivmodVHCState.quotientState] using hnextReady
  · simpa [projected, HenselDivmodVHCState.quotientState] using hnextPrefix
  · simpa [projected, HenselDivmodVHCState.quotientState] using
      hnextProcessed
  · simpa [projected, HenselDivmodVHCState.quotientState] using hnextConsumed
  · simpa [projected, HenselDivmodVHCState.quotientState] using
      hnextRemaining
  · simpa [projected, HenselDivmodVHCState.quotientState] using hnextBelow

theorem HenselDivmodVHCPrefix.fullOperationalInvariant
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (dividend divisor : SparsePolyZp) (hdivisor : 0 < divisor.size)
    (initialLimit : Nat) (initial : HenselDivmodVHCState)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hdividendCanonical : CLPoly.Math.SparsePolyZp.Canonical this._p.toNat
      dividend)
    (hinitial : HenselDivmodVHCFullOperationalInvariant this initialLimit
      initialLimit initial dividend divisor hdivisor)
    (degreeLimit : Nat) (state : HenselDivmodVHCState)
    (hprefix : HenselDivmodVHCPrefix this dividend divisor hdivisor initialLimit
      initial degreeLimit state) :
    HenselDivmodVHCFullOperationalInvariant this initialLimit degreeLimit state
      dividend divisor hdivisor := by
  induction hprefix with
  | refl => exact hinitial
  | step limit state next frontier hprefix hnotDone hselect hdecrease
      hiteration ih =>
      exact ih.step this initialLimit limit state next frontier dividend divisor
        hdivisor hcfg hdividendCanonical hselect hdecrease hiteration

theorem henselDivmodVHCFreshFullOperationalInvariant
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (dividend divisor : SparsePolyZp)
    (hdividendCanonical : CLPoly.Math.SparsePolyZp.Canonical this._p.toNat
      dividend)
    (hdivisorCanonical : CLPoly.Math.SparsePolyZp.Canonical this._p.toNat
      divisor)
    (hdividend : 0 < dividend.size) (hdivisor : 0 < divisor.size) :
    let limit := dividend[0].1.deg + 1
    let nodes := StrictSquarefreeZp.pairVecDivVHCInit divisor
    let state : HenselDivmodVHCState :=
      ⟨0, #[], nodes, #[], #[], divisor.size - 1⟩
    HenselDivmodVHCFullOperationalInvariant this limit limit state dividend
      divisor hdivisor := by
  dsimp only
  let nodes := StrictSquarefreeZp.pairVecDivVHCInit divisor
  rcases StrictSquarefreeZp.pairVecDivVHCInit_stateCovered divisor with
    ⟨owners, hownership0, hcovered0⟩
  have hownership : StrictSquarefreeZp.PairVecDivVHCHeapChainOwnership #[]
      owners nodes := by simpa [nodes] using hownership0
  have hstate : StrictSquarefreeZp.PairVecDivVHCStateCovered #[] nodes #[]
      (divisor.size - 1) := by
    exact ⟨owners, hownership, by simpa [nodes] using hcovered0⟩
  have hinactive : ∀ (i : Nat)
      (node : StrictSquarefreeZp.PairVecDivVHCNode),
      nodes[i]? = some node → node.mono = none := by
    intro i node hget
    have hi : i < divisor.size - 1 := by
      by_contra hnot
      rw [Array.getElem?_eq_none (by
        simpa [nodes, StrictSquarefreeZp.pairVecDivVHCInit_size] using hnot)]
        at hget
      contradiction
    rw [show nodes[i]? = some
        (StrictSquarefreeZp.pairVecDivVHCInitialNode i) by
      simpa [nodes] using
        StrictSquarefreeZp.pairVecDivVHCInit_get divisor i hi] at hget
    simp only [Option.some.injEq] at hget
    subst node
    rfl
  have hcanonical : CLPoly.Math.SparsePolyZp.Canonical this._p.toNat
      (#[] : SparsePolyZp) := by
    rw [CLPoly.Math.SparsePolyZp.Canonical,
      CLPoly.Math.SparsePolyZp.WellFormed_arr]
    refine ⟨?_, by simp, ?_⟩
    · intro term hterm
      simp at hterm
    · intro term hterm
      simp at hterm
  have hactiveBelow :
      StrictSquarefreeZp.PairVecDivVHCAllActiveNodesBelow
        (dividend[0].1.deg + 1) nodes := by
    intro i node mono hget hmono
    rw [hinactive i node hget] at hmono
    contradiction
  have hdenotes : ∀ (i : Nat)
      (node : StrictSquarefreeZp.PairVecDivVHCNode),
      nodes[i]? = some node → node.mono ≠ none →
        StrictSquarefreeZp.PairVecDivVHCNodeDenotes (#[] : SparsePolyZp)
          divisor node := by
    intro i node hget hmono
    exact (hmono (hinactive i node hget)).elim
  have hhomogeneous : StrictSquarefreeZp.PairVecDivVHCHeapChainsHomogeneous
      #[] owners nodes := by
    intro slot head mono hget
    simp at hget
  have hordered : StrictSquarefreeZp.PairVecDivVHCHeapOrdered #[] nodes := by
    intro child parent hchild
    simp at hchild
  refine ⟨⟨hcanonical, ?_, ?_, hstate, ⟨owners, hownership,
      hhomogeneous⟩, hordered, hdenotes, ?_, ?_, ?_, ?_⟩,
    hdivisorCanonical, ?_, ?_, hactiveBelow⟩
  · simpa [nodes] using StrictSquarefreeZp.pairVecDivVHCInit_size divisor
  · simpa [nodes] using
      StrictSquarefreeZp.pairVecDivVHCInit_divisorIndicesFixed divisor
  · simpa [nodes] using
      StrictSquarefreeZp.pairVecDivVHCInit_resetReady divisor
  · simpa [nodes] using
      (StrictSquarefreeZp.pairVecDivVHCInit_cursorPrefixAbove
        (dividend[0].1.deg + 1) (#[] : SparsePolyZp) divisor)
  · exact StrictSquarefreeZp.PairVecDivVHCQuotientLeadAbove.empty
      (dividend[0].1.deg + 1) divisor[0].1.deg
  · exact StrictSquarefreeZp.pairVecDivVHCConsumedDividendAbove_zero
      (dividend[0].1.deg + 1) dividend
  · intro frontier hselect
    exact StrictSquarefreeZp.PairVecDivVHCQuotientAbove.empty frontier.degree
      divisor[0].1.deg
  · exact StrictSquarefreeZp.pairVecDivVHCCanonicalInitialRemainingBelow
      this._p.toNat dividend hdividendCanonical hdividend

theorem henselDivmodVHCIteration_preserves_remainderCanonical
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (globalLimit degreeLimit : Nat) (state next : HenselDivmodVHCState)
    (frontier : StrictSquarefreeZp.PairVecDivVHCFrontier)
    (dividend divisor : SparsePolyZp) (hdivisor : 0 < divisor.size)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hdividendCanonical : CLPoly.Math.SparsePolyZp.Canonical this._p.toNat
      dividend)
    (hinvariant : HenselDivmodVHCFullOperationalInvariant this globalLimit
      degreeLimit state dividend divisor hdivisor)
    (hremainderCanonical : CLPoly.Math.SparsePolyZp.Canonical this._p.toNat
      state.remainder)
    (habove : HenselDivmodVHCRemainderAbove degreeLimit state.remainder)
    (hselect : StrictSquarefreeZp.pairVecDivVHCSelectFrontier
      state.dividendIndex dividend state.heap state.nodes = .ok frontier)
    (hdecrease : frontier.degree < degreeLimit)
    (hiteration : henselDivmodVHCIteration this state frontier dividend divisor
      hdivisor = .ok next) :
    CLPoly.Math.SparsePolyZp.Canonical this._p.toNat next.remainder ∧
      HenselDivmodVHCRemainderAbove frontier.degree next.remainder := by
  rcases henselDivmodVHCIteration_remainder this state next frontier dividend
      divisor hdivisor hiteration with ⟨consumed, hconsume, hremainder⟩
  have hfrontierReduced :=
    StrictSquarefreeZp.pairVecDivVHCSelectFrontier_coefficient_reduced
      this._p.toNat state.dividendIndex dividend state.heap state.nodes frontier
      (Fact.out : Nat.Prime this._p.toNat).pos hdividendCanonical hselect
  have hp : this._p ≠ 0 := by
    intro hp
    have hz : this._p.toNat = 0 := congrArg UInt64.toNat hp
    exact (Fact.out : Nat.Prime this._p.toNat).ne_zero hz
  have hcoefficientReduced :=
    StrictSquarefreeZp.pairVecDivVHCConsumeEqualDegree_coefficient_reduced this
      frontier.degree state.heap frontier.coefficient state.nodes #[]
      state.resetH state.quotient divisor consumed hp hcfg
      hinvariant.core.quotientCanonical hfrontierReduced hconsume
  rw [hremainder]
  by_cases hpush : consumed.coefficient ≠ 0 ∧
      frontier.degree < divisor[0].1.deg
  · rw [if_pos hpush]
    refine ⟨?_, HenselDivmodVHCRemainderAbove.push degreeLimit
      frontier.degree state.remainder ⟨consumed.coefficient, this._p⟩ habove
      hdecrease⟩
    exact CLPoly.Impl.StrictPolynomialGCDRefinement.canonical_push_lower
      this._p state.remainder frontier.degree consumed.coefficient
      hremainderCanonical (by
        intro term hterm
        have := habove term hterm
        omega) hcoefficientReduced hpush.1
  · rw [if_neg hpush]
    refine ⟨hremainderCanonical, ?_⟩
    intro term hterm
    exact Nat.le_trans (Nat.le_of_lt hdecrease) (habove term hterm)

theorem HenselDivmodVHCCorrect.remainderCanonical
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (dividend divisor : SparsePolyZp) (hdivisor : 0 < divisor.size)
    (initialLimit : Nat) (initial : HenselDivmodVHCState)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hdividendCanonical : CLPoly.Math.SparsePolyZp.Canonical this._p.toNat
      dividend)
    (hraw : ∀ limit state,
      HenselDivmodVHCPrefix this dividend divisor hdivisor initialLimit initial
          limit state →
        HenselDivmodVHCFullOperationalInvariant this initialLimit limit state
          dividend divisor hdivisor)
    (limit : Nat) (state : HenselDivmodVHCState)
    (output : SparsePolyZp × SparsePolyZp)
    (hcorrect : HenselDivmodVHCCorrect this dividend divisor hdivisor limit
      state output)
    (hreachable : HenselDivmodVHCPrefix this dividend divisor hdivisor
      initialLimit initial limit state)
    (hremainderCanonical : CLPoly.Math.SparsePolyZp.Canonical this._p.toNat
      state.remainder)
    (habove : HenselDivmodVHCRemainderAbove limit state.remainder) :
    CLPoly.Math.SparsePolyZp.Canonical this._p.toNat output.2 := by
  induction hcorrect with
  | done limit state hdone => exact hremainderCanonical
  | step limit state next frontier output hnotDone hselect hdecrease
      hiteration htail ih =>
      have hnextReachable : HenselDivmodVHCPrefix this dividend divisor
          hdivisor initialLimit initial frontier.degree next :=
        .step limit state next frontier hreachable hnotDone hselect hdecrease
          hiteration
      rcases henselDivmodVHCIteration_preserves_remainderCanonical this
          initialLimit limit state next frontier dividend divisor hdivisor hcfg
          hdividendCanonical (hraw limit state hreachable) hremainderCanonical
          habove hselect hdecrease hiteration with
        ⟨hnextCanonical, hnextAbove⟩
      exact ih hnextReachable hnextCanonical hnextAbove

/-- The exact generated VHC trace accumulates every low-degree residual
coefficient in its concrete remainder array.  The proof is well-founded along
the trace's selected frontier degrees.  Operational invariants are required
only for actually reachable prefixes and contain no expected quotient or
remainder. -/
theorem HenselDivmodVHCCorrect.remainderAgreesBelow
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (dividend divisor : SparsePolyZp) (hdivisor : 0 < divisor.size)
    (initialLimit : Nat) (initial : HenselDivmodVHCState)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hdividendCanonical : CLPoly.Math.SparsePolyZp.Canonical this._p.toNat
      dividend)
    (hraw : ∀ limit state,
      HenselDivmodVHCPrefix this dividend divisor hdivisor initialLimit initial
          limit state →
        HenselDivmodVHCOperationalInvariant this._p.toNat limit state dividend
          divisor hdivisor)
    (limit : Nat) (state : HenselDivmodVHCState)
    (output : SparsePolyZp × SparsePolyZp)
    (hcorrect : HenselDivmodVHCCorrect this dividend divisor hdivisor limit
      state output)
    (hreachable : HenselDivmodVHCPrefix this dividend divisor hdivisor
      initialLimit initial limit state)
    (habove : HenselDivmodVHCRemainderAbove limit state.remainder)
    (hagrees : HenselDivmodVHCRemainderAgreesAbove this._p.toNat limit
      divisor[0].1.deg state.remainder state.quotient dividend divisor) :
    ∀ degree, degree < divisor[0].1.deg →
      (CLPoly.Math.SparsePolyZp.toPoly this._p.toNat output.2).coeff degree =
        (CLPoly.Math.SparsePolyZp.toPoly this._p.toNat dividend -
          CLPoly.Math.SparsePolyZp.toPoly this._p.toNat output.1 *
            CLPoly.Math.SparsePolyZp.toPoly this._p.toNat divisor).coeff
              degree := by
  induction hcorrect with
  | done limit state hdone =>
      have hstate := hraw limit state hreachable
      rcases hstate.heapChains with ⟨owners, hownership, hhomogeneous⟩
      have hheapEmpty : state.heap = #[] :=
        Array.eq_empty_of_size_eq_zero hdone.2
      intro degree hlead
      by_cases hold : limit ≤ degree
      · exact hagrees degree hold hlead
      · rw [HenselDivmodVHCRemainderAbove.coeff_eq_zero_below
          this._p.toNat limit degree state.remainder habove (by omega),
        henselDivmodVHCResidual_coeff_eq_zero_of_done this._p.toNat limit
          degree state.dividendIndex state.resetH dividend state.quotient
          divisor state.nodes owners hdivisor hstate.size hstate.fixed
          (by simpa [hheapEmpty] using hstate.covered)
          (by simpa [hheapEmpty] using hownership) hstate.resetReady
          hstate.cursorPrefix hstate.quotientProcessed hdividendCanonical
          hstate.dividendConsumed hdone.1 (by omega)]
  | step limit state next frontier output hnotDone hselect hdecrease
      hiteration htail ih =>
      have hstate := hraw limit state hreachable
      have hnextReachable : HenselDivmodVHCPrefix this dividend divisor
          hdivisor initialLimit initial frontier.degree next :=
        .step limit state next frontier hreachable hnotDone hselect hdecrease
          hiteration
      rcases hstate.heapChains with ⟨owners, hownership, hhomogeneous⟩
      by_cases hbelowLead : frontier.degree < divisor[0].1.deg
      · rcases henselDivmodVHCIteration_below_lead_residual this limit state
            next frontier dividend divisor owners hdivisor hstate.size
            hstate.covered hownership hhomogeneous hstate.ordered hstate.denotes
            hstate.fixed hstate.resetReady hstate.cursorPrefix hcfg
            hstate.quotientCanonical hdividendCanonical
            hstate.dividendConsumed hdecrease hbelowLead hselect hiteration with
          ⟨consumed, hconsume, hresidual, hnextRemainder⟩
        have hgap : ∀ degree, frontier.degree < degree → degree < limit →
            (CLPoly.Math.SparsePolyZp.toPoly this._p.toNat dividend -
              CLPoly.Math.SparsePolyZp.toPoly this._p.toNat state.quotient *
                CLPoly.Math.SparsePolyZp.toPoly this._p.toNat divisor).coeff
                  degree = 0 := by
          intro degree hfrontier hlimit
          exact henselDivmodVHCResidual_coeff_eq_zero_of_gap
            this._p.toNat limit degree state.dividendIndex state.resetH dividend
            state.quotient divisor state.heap state.nodes owners frontier
            hdivisor hstate.size hstate.fixed hstate.covered hownership
            hhomogeneous hstate.resetReady hstate.ordered hstate.denotes
            hdividendCanonical hstate.dividendConsumed
            hstate.quotientCanonical hstate.cursorPrefix
            hstate.quotientProcessed hfrontier hlimit hselect
        rcases henselDivmodRemainderStep_preserves this._p.toNat limit
            frontier.degree divisor[0].1.deg state.remainder next.remainder
            state.quotient dividend divisor consumed.coefficient this._p
            hdecrease hbelowLead habove hagrees (by
              simpa [Polynomial.coeff_sub] using hresidual) hgap
            hnextRemainder with ⟨hnextAbove, hnextAgrees⟩
        rcases henselDivmodVHCIteration_projects_quotient this state next
            frontier dividend divisor hdivisor hselect hiteration with
          ⟨projected, hprojected, hindex, hheap, hnodes, hquotient, hreset⟩
        have hprojectedEq :=
          StrictSquarefreeZp.pairVecDivVHCOuterIteration_quotient_eq_of_frontier_lt_lead
            this state.dividendIndex state.heap state.nodes state.quotient
            dividend divisor state.resetH frontier projected hdivisor
            hbelowLead hselect hprojected
        have hnextQuotient : next.quotient = state.quotient := by
          rw [← hquotient]
          exact hprojectedEq
        apply ih hnextReachable hnextAbove
        simpa [hnextQuotient] using hnextAgrees
      · rcases henselDivmodVHCIteration_remainder this state next frontier
            dividend divisor hdivisor hiteration with
          ⟨consumed, hconsume, hnextRemainder⟩
        have hnextRemainderEq : next.remainder = state.remainder := by
          rw [hnextRemainder]
          simp [hbelowLead]
        apply ih hnextReachable
        · rw [hnextRemainderEq]
          intro term hterm
          exact Nat.le_trans (Nat.le_of_lt hdecrease) (habove term hterm)
        · intro degree hfrontier hlead
          omega

/-- High-degree quotient agreement and the exact accumulated low-degree
remainder together give the standard polynomial division equation. -/
theorem henselDivmodVHCDivisionEquation
    (p leadDegree : Nat)
    (quotient remainder dividend divisor : SparsePolyZp)
    (hproduct : StrictSquarefreeZp.PairVecDivVHCProductAgreesAbove p
      leadDegree quotient dividend divisor)
    (hremainderBelow : HenselDivmodVHCRemainderBelow leadDegree remainder)
    (hremainderAgrees : ∀ degree, degree < leadDegree →
      (CLPoly.Math.SparsePolyZp.toPoly p remainder).coeff degree =
        (CLPoly.Math.SparsePolyZp.toPoly p dividend -
          CLPoly.Math.SparsePolyZp.toPoly p quotient *
            CLPoly.Math.SparsePolyZp.toPoly p divisor).coeff degree) :
    CLPoly.Math.SparsePolyZp.toPoly p dividend =
      CLPoly.Math.SparsePolyZp.toPoly p quotient *
          CLPoly.Math.SparsePolyZp.toPoly p divisor +
        CLPoly.Math.SparsePolyZp.toPoly p remainder := by
  ext degree
  rw [Polynomial.coeff_add]
  by_cases hdegree : leadDegree ≤ degree
  · rw [HenselDivmodVHCRemainderBelow.coeff_eq_zero_at_or_above p
      leadDegree degree remainder hremainderBelow hdegree, add_zero]
    exact (hproduct degree hdegree).symm
  · have hbelow : degree < leadDegree := by omega
    rw [hremainderAgrees degree hbelow, Polynomial.coeff_sub]
    abel

/-- Every successful concrete outer-loop execution determines its exact
source trace. -/
theorem henselDivmodVHCOuterLoop_correct_of_success
    (this : DenseUPolyZp) (dividend divisor : SparsePolyZp)
    (hdivisor : 0 < divisor.size) (limit : Nat)
    (state : HenselDivmodVHCState) (output : SparsePolyZp × SparsePolyZp)
    (hrun : henselDivmodVHCOuterLoop this limit state dividend divisor
      hdivisor = .ok output) :
    HenselDivmodVHCCorrect this dividend divisor hdivisor limit state output := by
  induction limit using Nat.strong_induction_on generalizing state output with
  | h limit ih =>
      rw [henselDivmodVHCOuterLoop] at hrun
      by_cases hdone : dividend.size ≤ state.dividendIndex ∧
          state.heap.size = 0
      · rw [dif_pos hdone] at hrun
        have houtput : (state.quotient, state.remainder) = output :=
          Except.ok.inj hrun
        subst output
        exact .done limit state hdone
      · rw [dif_neg hdone] at hrun
        cases hselect : StrictSquarefreeZp.pairVecDivVHCSelectFrontier
            state.dividendIndex dividend state.heap state.nodes with
        | error fault => simp [hselect] at hrun
        | ok frontier =>
            rw [hselect] at hrun
            simp only at hrun
            by_cases hdecrease : frontier.degree < limit
            · rw [dif_pos hdecrease] at hrun
              cases hiteration : henselDivmodVHCIteration this state frontier
                  dividend divisor hdivisor with
              | error fault => simp [hiteration] at hrun
              | ok next =>
                  rw [hiteration] at hrun
                  simp only at hrun
                  exact .step limit state next frontier output hdone hselect
                    hdecrease hiteration
                    (ih frontier.degree hdecrease next output hrun)
            · rw [dif_neg hdecrease] at hrun
              contradiction

/-- Entry-level semantic theorem for the fresh general five-argument VHC
branch.  All state facts are derived from initialization and actual reachable
bodies. -/
theorem henselDivmodVHCGeneralBranchIR_refines
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (dividend divisor : SparsePolyZp) (output : SparsePolyZp × SparsePolyZp)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hdividendCanonical : CLPoly.Math.SparsePolyZp.Canonical this._p.toNat
      dividend)
    (hdivisorCanonical : CLPoly.Math.SparsePolyZp.Canonical this._p.toNat
      divisor)
    (hdividend : 0 < dividend.size) (hdivisor : 1 < divisor.size)
    (hrun : henselDivmodVHCGeneralBranchIR this dividend divisor =
      .ok output) :
    CLPoly.Math.SparsePolyZp.Canonical this._p.toNat output.1 ∧
      CLPoly.Math.SparsePolyZp.Canonical this._p.toNat output.2 ∧
      HenselDivmodVHCRemainderBelow divisor[0].1.deg output.2 ∧
      CLPoly.Math.SparsePolyZp.toPoly this._p.toNat output.1 *
            CLPoly.Math.SparsePolyZp.toPoly this._p.toNat divisor +
          CLPoly.Math.SparsePolyZp.toPoly this._p.toNat output.2 =
        CLPoly.Math.SparsePolyZp.toPoly this._p.toNat dividend := by
  let limit := dividend[0].1.deg + 1
  let nodes := StrictSquarefreeZp.pairVecDivVHCInit divisor
  let initial : HenselDivmodVHCState :=
    ⟨0, #[], nodes, #[], #[], divisor.size - 1⟩
  have houter : henselDivmodVHCOuterLoop this limit initial dividend divisor
      (by omega) = .ok output := by
    simpa [henselDivmodVHCGeneralBranchIR, hdividend, hdivisor, limit, nodes,
      initial] using hrun
  have hinitial := henselDivmodVHCFreshFullOperationalInvariant this dividend
    divisor hdividendCanonical hdivisorCanonical hdividend (by omega)
  have hinitial' : HenselDivmodVHCFullOperationalInvariant this limit limit
      initial dividend divisor (by omega) := by
    simpa [limit, nodes, initial] using hinitial
  rcases hinitial'.core.heapChains with
    ⟨owners, hownership, hhomogeneous⟩
  have hquotientRun := henselDivmodVHCOuterLoop_projects_quotient this dividend
    divisor (by omega) limit initial output houter
  have hquotientCanonical :=
    StrictSquarefreeZp.pairVecDivVHCOuterLoop_preserves_canonical_of_success
      this limit limit initial.dividendIndex initial.heap initial.nodes
      initial.quotient dividend divisor output.1 initial.resetH owners
      (by omega) hcfg hinitial'.core.quotientCanonical hdividendCanonical
      hdivisorCanonical hinitial'.quotientReady hinitial'.remaining
      hinitial'.activeBelow hinitial'.core.denotes hinitial'.core.fixed
      hinitial'.core.resetReady hownership hinitial'.core.covered hquotientRun
  have hinitialAgrees : StrictSquarefreeZp.PairVecDivVHCProductAgreesAbove
      this._p.toNat limit initial.quotient dividend divisor := by
    intro degree hdegree
    have hpolyDegree := StrictSquarefreeZp.sparsePolyZp_toPoly_degree_eq_head
      this._p.toNat dividend hdividendCanonical hdividend
    have hdegreeLt :
        (CLPoly.Math.SparsePolyZp.toPoly this._p.toNat dividend).degree <
          degree := by
      rw [hpolyDegree]
      change dividend[0].1.deg + 1 ≤ degree at hdegree
      exact_mod_cast Nat.lt_of_succ_le hdegree
    rw [Polynomial.degree_lt_iff_coeff_zero] at hdegreeLt
    have hzero := hdegreeLt degree (Nat.le_refl degree)
    simpa [initial, CLPoly.Math.SparsePolyZp.toPoly, CLPoly.Math.listSum]
      using hzero.symm
  have hproduct :=
    henselDivmodVHCOuterLoop_productAgreesAbove_lead_of_success this limit limit
      initial dividend divisor output owners (by omega) hcfg
      hinitial'.core.quotientCanonical hdividendCanonical hdivisorCanonical
      hinitial'.quotientReady hinitial'.remaining hinitial'.activeBelow
      hinitial'.core.denotes hinitial'.core.fixed hinitial'.core.resetReady
      hownership hinitial'.core.covered hinitial'.core.size hhomogeneous
      hinitial'.core.ordered hinitial'.core.cursorPrefix
      hinitial'.core.quotientProcessed hinitial'.core.dividendConsumed
      hinitialAgrees houter
  have hcorrect := henselDivmodVHCOuterLoop_correct_of_success this dividend
    divisor (by omega) limit initial output houter
  have hremainderBelow := hcorrect.remainderBelow this dividend divisor
    (by omega) limit initial output
    (HenselDivmodVHCRemainderBelow.empty divisor[0].1.deg)
  have hfullRaw : ∀ degreeLimit state,
      HenselDivmodVHCPrefix this dividend divisor (by omega) limit initial
          degreeLimit state →
        HenselDivmodVHCFullOperationalInvariant this limit degreeLimit state
          dividend divisor (by omega) := by
    intro degreeLimit state hprefix
    exact hprefix.fullOperationalInvariant this dividend divisor (by omega)
      limit initial hcfg hdividendCanonical hinitial' degreeLimit state
  have hraw : ∀ degreeLimit state,
      HenselDivmodVHCPrefix this dividend divisor (by omega) limit initial
          degreeLimit state →
        HenselDivmodVHCOperationalInvariant this._p.toNat degreeLimit state
          dividend divisor (by omega) := by
    intro degreeLimit state hprefix
    exact (hfullRaw degreeLimit state hprefix).core
  have hremainderCanonical := hcorrect.remainderCanonical this dividend divisor
    (by omega) limit initial hcfg hdividendCanonical hfullRaw limit initial
    output HenselDivmodVHCPrefix.refl (by
      simpa [initial] using hinitial'.core.quotientCanonical)
      (HenselDivmodVHCRemainderAbove.empty limit)
  have hinitialRemainderAgrees : HenselDivmodVHCRemainderAgreesAbove
      this._p.toNat limit divisor[0].1.deg initial.remainder initial.quotient
      dividend divisor := by
    intro degree hdegree hlead
    have hpolyDegree := StrictSquarefreeZp.sparsePolyZp_toPoly_degree_eq_head
      this._p.toNat dividend hdividendCanonical hdividend
    have hdegreeLt :
        (CLPoly.Math.SparsePolyZp.toPoly this._p.toNat dividend).degree <
          degree := by
      rw [hpolyDegree]
      change dividend[0].1.deg + 1 ≤ degree at hdegree
      exact_mod_cast Nat.lt_of_succ_le hdegree
    rw [Polynomial.degree_lt_iff_coeff_zero] at hdegreeLt
    have hzero := hdegreeLt degree (Nat.le_refl degree)
    simpa [initial, CLPoly.Math.SparsePolyZp.toPoly, CLPoly.Math.listSum]
      using hzero.symm
  have hremainderAgrees := hcorrect.remainderAgreesBelow this dividend divisor
    (by omega) limit initial hcfg hdividendCanonical hraw limit initial output
    HenselDivmodVHCPrefix.refl
    (HenselDivmodVHCRemainderAbove.empty limit) hinitialRemainderAgrees
  exact ⟨hquotientCanonical, hremainderCanonical, hremainderBelow,
    (henselDivmodVHCDivisionEquation this._p.toNat divisor[0].1.deg output.1
      output.2 dividend divisor hproduct hremainderBelow
      hremainderAgrees).symm⟩

structure HenselDivmodVHCResultCorrect (p : Nat)
    (dividend divisor : SparsePolyZp) (hdivisor : 0 < divisor.size)
    (output : SparsePolyZp × SparsePolyZp) : Prop where
  quotientCanonical : CLPoly.Math.SparsePolyZp.Canonical p output.1
  remainderCanonical : CLPoly.Math.SparsePolyZp.Canonical p output.2
  remainderBelow : HenselDivmodVHCRemainderBelow divisor[0].1.deg output.2
  equation : CLPoly.Math.SparsePolyZp.toPoly p output.1 *
        CLPoly.Math.SparsePolyZp.toPoly p divisor +
      CLPoly.Math.SparsePolyZp.toPoly p output.2 =
    CLPoly.Math.SparsePolyZp.toPoly p dividend

/-- Unified total semantic contract for the exact complete five-argument
source entry. -/
theorem henselDivmodVHCIR_refines
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (dividend divisor : SparsePolyZp)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hdividendCanonical : CLPoly.Math.SparsePolyZp.Canonical this._p.toNat
      dividend)
    (hdivisorCanonical : CLPoly.Math.SparsePolyZp.Canonical this._p.toNat
      divisor)
    (hdivisor : 0 < divisor.size) :
    ∃ output, henselDivmodVHCIR this dividend divisor = .ok output ∧
      HenselDivmodVHCResultCorrect this._p.toNat dividend divisor hdivisor
        output := by
  by_cases hdividend : dividend.size = 0
  · have hempty : dividend = #[] := Array.size_eq_zero_iff.mp hdividend
    subst dividend
    refine ⟨(#[], #[]), by
      simp [henselDivmodVHCIR, Nat.ne_of_gt hdivisor], ?_⟩
    exact ⟨hdividendCanonical, hdividendCanonical,
      HenselDivmodVHCRemainderBelow.empty divisor[0].1.deg, by
        simp [CLPoly.Math.SparsePolyZp.toPoly, CLPoly.Math.listSum]⟩
  · by_cases hsingle : divisor.size = 1
    · rcases henselDivmodVHCSingleBranchIR_refines this dividend divisor hcfg
          hdividendCanonical hdivisorCanonical hsingle with
        ⟨output, hrun, hquotientCanonical, hremainderCanonical,
          hremainderBelow, hequation⟩
      refine ⟨output, ?_, ⟨hquotientCanonical, hremainderCanonical,
        hremainderBelow, hequation⟩⟩
      simpa [henselDivmodVHCIR, Nat.ne_of_gt hdivisor, hdividend, hsingle]
        using hrun
    · have hgeneral : 1 < divisor.size := by omega
      rcases henselDivmodVHCGeneralBranchIR_succeeds this dividend divisor hcfg
          hdividendCanonical hdivisorCanonical (Nat.pos_of_ne_zero hdividend)
          hgeneral with ⟨output, hrun⟩
      rcases henselDivmodVHCGeneralBranchIR_refines this dividend divisor output
          hcfg hdividendCanonical hdivisorCanonical
          (Nat.pos_of_ne_zero hdividend) hgeneral hrun with
        ⟨hquotientCanonical, hremainderCanonical, hremainderBelow,
          hequation⟩
      refine ⟨output, ?_, ⟨hquotientCanonical, hremainderCanonical,
        hremainderBelow, hequation⟩⟩
      simpa [henselDivmodVHCIR, Nat.ne_of_gt hdivisor, hdividend, hsingle]
        using hrun

/-- Remainder degree is an execution property, independent of canonicality:
every successful complete source call strictly bounds stored remainder terms
below the nonempty divisor lead. -/
theorem henselDivmodVHCIR_remainderBelow_of_success
    (this : DenseUPolyZp) (dividend divisor : SparsePolyZp)
    (quotient remainder : SparsePolyZp) (hdivisor : 0 < divisor.size)
    (hrun : henselDivmodVHCIR this dividend divisor =
      .ok (quotient, remainder)) :
    HenselDivmodVHCRemainderBelow divisor[0].1.deg remainder := by
  by_cases hdividend : dividend.size = 0
  · simp [henselDivmodVHCIR, Nat.ne_of_gt hdivisor, hdividend] at hrun
    rcases hrun with ⟨rfl, rfl⟩
    exact HenselDivmodVHCRemainderBelow.empty divisor[0].1.deg
  · by_cases hsingle : divisor.size = 1
    · simp only [henselDivmodVHCIR, Nat.ne_of_gt hdivisor, hdividend,
        hsingle, henselDivmodVHCSingleBranchIR, ↓reduceDIte] at hrun
      have houtput := Except.ok.inj hrun
      have hbelow := henselDivmodVHCSingleLoopIR_remainderBelow this dividend
        divisor[0]
      rw [houtput] at hbelow
      exact hbelow
    · have hgeneral : 1 < divisor.size := by omega
      have hdividendPos : 0 < dividend.size := Nat.pos_of_ne_zero hdividend
      have houter :
          let nodes := StrictSquarefreeZp.pairVecDivVHCInit divisor
          let state : HenselDivmodVHCState :=
            ⟨0, #[], nodes, #[], #[], divisor.size - 1⟩
          henselDivmodVHCOuterLoop this (dividend[0].1.deg + 1) state dividend
            divisor (by omega) = .ok (quotient, remainder) := by
        simpa [henselDivmodVHCIR, Nat.ne_of_gt hdivisor, hdividend, hsingle,
          henselDivmodVHCGeneralBranchIR, hdividendPos, hgeneral] using hrun
      dsimp only at houter
      let nodes := StrictSquarefreeZp.pairVecDivVHCInit divisor
      let state : HenselDivmodVHCState :=
        ⟨0, #[], nodes, #[], #[], divisor.size - 1⟩
      have hcorrect := henselDivmodVHCOuterLoop_correct_of_success this dividend
        divisor (by omega) (dividend[0].1.deg + 1) state
        (quotient, remainder) (by simpa [state, nodes] using houter)
      exact hcorrect.remainderBelow this dividend divisor (by omega)
        (dividend[0].1.deg + 1) state (quotient, remainder)
        (HenselDivmodVHCRemainderBelow.empty divisor[0].1.deg)

/-- The actual source EEA decreases by the leading degree of its current
second remainder.  Empty remainders receive measure zero. -/
def henselEEAMeasure
    (state : Generated.StrictHensel.HenselEEAState) : Nat :=
  if h : 0 < state.r1.size then state.r1[0].1.deg + 1 else 0

/-- Well-founded evidence for the concrete Hensel EEA division operation.
The decrease is derived from every successful raw VHC remainder, rather than
from a fuel counter or a separately supplied specification result. -/
def strictHenselEEATermination (this : DenseUPolyZp) :
    Generated.StrictHensel.HenselEEATermination
      (strictHenselEEARawOps this) where
  measure := henselEEAMeasure
  decreases := by
    intro state quotient remainder hcontinue hrun
    have hr1Ne : state.r1 ≠ #[] := by
      simpa [Array.isEmpty_iff] using hcontinue
    have hr1 : 0 < state.r1.size := by
      have hr1Size : state.r1.size ≠ 0 := by
        intro hzero
        apply hr1Ne
        exact Array.eq_empty_of_size_eq_zero hzero
      omega
    have hbelow := henselDivmodVHCIR_remainderBelow_of_success this state.r0
      state.r1 quotient remainder hr1 (by
        simpa [strictHenselEEARawOps] using hrun)
    by_cases hremZero : remainder.size = 0
    · simp [henselEEAMeasure,
        Generated.StrictHensel.henselEEANextState, hr1, hremZero]
    · have hrem : 0 < remainder.size := Nat.pos_of_ne_zero hremZero
      have hleadLt : remainder[0].1.deg < state.r1[0].1.deg :=
        hbelow remainder[0] (Array.getElem_mem_toList hrem)
      simp only [henselEEAMeasure,
        Generated.StrictHensel.henselEEANextState]
      rw [dif_pos hrem, dif_pos hr1]
      omega

/-- Canonicality and nonemptiness facts transported by the actual EEA
quotient/remainder calls. -/
structure StrictHenselEEAStateInvariant (p : Nat)
    (state : Generated.StrictHensel.HenselEEAState) : Prop where
  r0Canonical : CLPoly.Math.SparsePolyZp.Canonical p state.r0
  r1Canonical : CLPoly.Math.SparsePolyZp.Canonical p state.r1
  r0Nonempty : 0 < state.r0.size

theorem strictHenselEEAPrefix_stateInvariant
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (initial : Generated.StrictHensel.HenselEEAState)
    (hinitial : StrictHenselEEAStateInvariant this._p.toNat initial) :
    ∀ state, HenselEEAPrefix (strictHenselEEARawOps this) initial state →
      StrictHenselEEAStateInvariant this._p.toNat state := by
  intro state hprefix
  induction hprefix with
  | refl => exact hinitial
  | step state quotient remainder hprefix hcontinue hrun ih =>
      have hr1Ne : state.r1 ≠ #[] := by
        simpa [Array.isEmpty_iff] using hcontinue
      have hr1 : 0 < state.r1.size := by
        have hr1Size : state.r1.size ≠ 0 := by
          intro hzero
          apply hr1Ne
          exact Array.eq_empty_of_size_eq_zero hzero
        omega
      rcases henselDivmodVHCIR_refines this state.r0 state.r1 hcfg
          ih.r0Canonical ih.r1Canonical hr1 with
        ⟨output, houtputRun, houtputCorrect⟩
      have hrun' : henselDivmodVHCIR this state.r0 state.r1 =
          .ok (quotient, remainder) := by
        simpa [strictHenselEEARawOps] using hrun
      have houtput : output = (quotient, remainder) :=
        Except.ok.inj (houtputRun.symm.trans hrun')
      subst output
      exact ⟨ih.r1Canonical, houtputCorrect.remainderCanonical, hr1⟩

/-- Concrete raw-to-safe readiness for the well-founded Hensel EEA.  The
invariant is derived along actual successful VHC calls. -/
def strictHenselEEAExecutionInvariant
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (initial : Generated.StrictHensel.HenselEEAState)
    (hinitial : StrictHenselEEAStateInvariant this._p.toNat initial) :
    HenselEEAExecutionInvariant (strictHenselEEARawOps this) initial where
  divisionReady := by
    intro state hprefix hcontinue
    have hstate := strictHenselEEAPrefix_stateInvariant this hcfg initial
      hinitial state hprefix
    have hr1Ne : state.r1 ≠ #[] := by
      simpa [Array.isEmpty_iff] using hcontinue
    have hr1 : 0 < state.r1.size := by
      have hr1Size : state.r1.size ≠ 0 := by
        intro hzero
        apply hr1Ne
        exact Array.eq_empty_of_size_eq_zero hzero
      omega
    rcases henselDivmodVHCIR_refines this state.r0 state.r1 hcfg
        hstate.r0Canonical hstate.r1Canonical hr1 with
      ⟨⟨quotient, remainder⟩, hrun, _⟩
    exact ⟨quotient, remainder, by
      simpa [strictHenselEEARawOps] using hrun⟩
  finalNonempty := by
    intro state hprefix _
    have hstate := strictHenselEEAPrefix_stateInvariant this hcfg initial
      hinitial state hprefix
    exact ⟨state.r0[0]'hstate.r0Nonempty,
      Array.getElem?_eq_getElem hstate.r0Nonempty⟩

/-- The generated EEA now runs with the concrete VHC divmod implementation
and its degree-based well-founded recursion. -/
theorem strictHenselEEAIR_refines
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (initial : Generated.StrictHensel.HenselEEAState)
    (hinitial : StrictHenselEEAStateInvariant this._p.toNat initial) :
    ∃ output,
      Generated.StrictHensel.__polynomial_GCD_eea_raw_ir
          (strictHenselEEARawOps this) (strictHenselEEATermination this)
          initial = .ok output ∧
      HenselEEACorrect (strictHenselEEARawOps this) initial output :=
  __polynomial_GCD_eea_raw_ir_refines (strictHenselEEARawOps this)
    (strictHenselEEATermination this) initial
    (strictHenselEEAExecutionInvariant this hcfg initial hinitial)

theorem henselEEANormalization_toPoly (p : Nat) (f : SparsePolyZp) :
    CLPoly.Math.SparsePolyZp.toPoly p (SparsePolyZp.normalization f) =
      CLPoly.Math.SparsePolyZp.toPoly p f := by
  unfold SparsePolyZp.normalization CLPoly.Math.SparsePolyZp.toPoly
  rw [Array.toList_filter]
  induction f.toList with
  | nil => simp [CLPoly.Math.listSum]
  | cons term rest ih =>
      rcases term with ⟨monomial, coefficient⟩
      by_cases hzero : coefficient.val = 0
      · simpa [hzero, CLPoly.Math.listSum, CLPoly.Math.Zp.toZMod] using ih
      · simp [hzero, CLPoly.Math.listSum, ih]

theorem henselEEANormalization_wellFormed (p : Nat) (f : SparsePolyZp)
    (hf : CLPoly.Math.SparsePolyZp.WellFormed_arr p f) :
    CLPoly.Math.SparsePolyZp.WellFormed_arr p
      (SparsePolyZp.normalization f) := by
  intro term hterm
  unfold SparsePolyZp.normalization at hterm
  rw [Array.toList_filter] at hterm
  exact hf term (List.mem_of_mem_filter hterm)

/-- The two source remainder registers are represented by the corresponding
Bézout coefficient registers throughout the actual EEA assignments. -/
structure StrictHenselEEAAlgebraicInvariant (p : Nat)
    (left right : SparsePolyZp)
    (state : Generated.StrictHensel.HenselEEAState) : Prop where
  s0WellFormed : CLPoly.Math.SparsePolyZp.WellFormed_arr p state.s0
  s1WellFormed : CLPoly.Math.SparsePolyZp.WellFormed_arr p state.s1
  t0WellFormed : CLPoly.Math.SparsePolyZp.WellFormed_arr p state.t0
  t1WellFormed : CLPoly.Math.SparsePolyZp.WellFormed_arr p state.t1
  r0Equation :
    CLPoly.Math.SparsePolyZp.toPoly p state.r0 =
      CLPoly.Math.SparsePolyZp.toPoly p state.s0 *
          CLPoly.Math.SparsePolyZp.toPoly p left +
        CLPoly.Math.SparsePolyZp.toPoly p state.t0 *
          CLPoly.Math.SparsePolyZp.toPoly p right
  r1Equation :
    CLPoly.Math.SparsePolyZp.toPoly p state.r1 =
      CLPoly.Math.SparsePolyZp.toPoly p state.s1 *
          CLPoly.Math.SparsePolyZp.toPoly p left +
        CLPoly.Math.SparsePolyZp.toPoly p state.t1 *
          CLPoly.Math.SparsePolyZp.toPoly p right

/-- One concrete EEA assignment block preserves the Bézout equations.  The
new remainder equation uses the quotient and remainder returned by the actual
VHC call. -/
theorem StrictHenselEEAAlgebraicInvariant.step
    (p : Nat) (left right : SparsePolyZp)
    (state : Generated.StrictHensel.HenselEEAState)
    (quotient remainder : SparsePolyZp)
    (h2p : 2 * p ≤ UInt64.size) (hp2 : p * p ≤ UInt64.size)
    (hr1 : 0 < state.r1.size)
    (hinvariant : StrictHenselEEAAlgebraicInvariant p left right state)
    (hdivmod : HenselDivmodVHCResultCorrect p state.r0 state.r1
      hr1 (quotient, remainder)) :
    StrictHenselEEAAlgebraicInvariant p left right
      (Generated.StrictHensel.henselEEANextState state quotient remainder) := by
  let qs := SparsePolyZp.mulImpl quotient state.s1
  let qt := SparsePolyZp.mulImpl quotient state.t1
  let sRaw := SparsePolyZp.subImpl state.s0 qs
  let tRaw := SparsePolyZp.subImpl state.t0 qt
  let s2 := SparsePolyZp.normalization sRaw
  let t2 := SparsePolyZp.normalization tRaw
  have hqWellFormed := hdivmod.quotientCanonical.1
  have hqsWellFormed : CLPoly.Math.SparsePolyZp.WellFormed_arr p qs := by
    exact CLPoly.Math.SparsePolyZp.WellFormed_arr.mul p hp2 quotient state.s1
      hqWellFormed hinvariant.s1WellFormed
  have hqtWellFormed : CLPoly.Math.SparsePolyZp.WellFormed_arr p qt := by
    exact CLPoly.Math.SparsePolyZp.WellFormed_arr.mul p hp2 quotient state.t1
      hqWellFormed hinvariant.t1WellFormed
  have hsRawWellFormed : CLPoly.Math.SparsePolyZp.WellFormed_arr p sRaw := by
    exact CLPoly.Math.SparsePolyZp.WellFormed_arr.sub p state.s0 qs
      hinvariant.s0WellFormed hqsWellFormed
  have htRawWellFormed : CLPoly.Math.SparsePolyZp.WellFormed_arr p tRaw := by
    exact CLPoly.Math.SparsePolyZp.WellFormed_arr.sub p state.t0 qt
      hinvariant.t0WellFormed hqtWellFormed
  have hs2WellFormed : CLPoly.Math.SparsePolyZp.WellFormed_arr p s2 :=
    henselEEANormalization_wellFormed p sRaw hsRawWellFormed
  have ht2WellFormed : CLPoly.Math.SparsePolyZp.WellFormed_arr p t2 :=
    henselEEANormalization_wellFormed p tRaw htRawWellFormed
  have hqsPoly : CLPoly.Math.SparsePolyZp.toPoly p qs =
      CLPoly.Math.SparsePolyZp.toPoly p quotient *
        CLPoly.Math.SparsePolyZp.toPoly p state.s1 := by
    exact CLPoly.Math.SparsePolyZp.toPoly_mul p h2p hp2 quotient state.s1
      hqWellFormed hinvariant.s1WellFormed
  have hqtPoly : CLPoly.Math.SparsePolyZp.toPoly p qt =
      CLPoly.Math.SparsePolyZp.toPoly p quotient *
        CLPoly.Math.SparsePolyZp.toPoly p state.t1 := by
    exact CLPoly.Math.SparsePolyZp.toPoly_mul p h2p hp2 quotient state.t1
      hqWellFormed hinvariant.t1WellFormed
  have hs2Poly : CLPoly.Math.SparsePolyZp.toPoly p s2 =
      CLPoly.Math.SparsePolyZp.toPoly p state.s0 -
        CLPoly.Math.SparsePolyZp.toPoly p quotient *
          CLPoly.Math.SparsePolyZp.toPoly p state.s1 := by
    rw [henselEEANormalization_toPoly p sRaw]
    change CLPoly.Math.SparsePolyZp.toPoly p (state.s0 - qs) = _
    rw [CLPoly.Math.SparsePolyZp.toPoly_sub p h2p state.s0 qs
      hinvariant.s0WellFormed hqsWellFormed, hqsPoly]
  have ht2Poly : CLPoly.Math.SparsePolyZp.toPoly p t2 =
      CLPoly.Math.SparsePolyZp.toPoly p state.t0 -
        CLPoly.Math.SparsePolyZp.toPoly p quotient *
          CLPoly.Math.SparsePolyZp.toPoly p state.t1 := by
    rw [henselEEANormalization_toPoly p tRaw]
    change CLPoly.Math.SparsePolyZp.toPoly p (state.t0 - qt) = _
    rw [CLPoly.Math.SparsePolyZp.toPoly_sub p h2p state.t0 qt
      hinvariant.t0WellFormed hqtWellFormed, hqtPoly]
  change StrictHenselEEAAlgebraicInvariant p left right
    { r0 := state.r1, r1 := remainder, s0 := state.s1, s1 := s2,
      t0 := state.t1, t1 := t2 }
  refine ⟨hinvariant.s1WellFormed, hs2WellFormed,
    hinvariant.t1WellFormed, ht2WellFormed, hinvariant.r1Equation, ?_⟩
  have hremEquation :
      CLPoly.Math.SparsePolyZp.toPoly p remainder =
        CLPoly.Math.SparsePolyZp.toPoly p state.r0 -
          CLPoly.Math.SparsePolyZp.toPoly p quotient *
            CLPoly.Math.SparsePolyZp.toPoly p state.r1 := by
    rw [← hdivmod.equation]
    ring
  rw [hremEquation, hinvariant.r0Equation, hinvariant.r1Equation,
    hs2Poly, ht2Poly]
  ring

/-- Every state reached by the concrete generated EEA maintains both Bézout
representations. -/
theorem strictHenselEEAPrefix_algebraicInvariant
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (h2p : 2 * this._p.toNat ≤ UInt64.size)
    (hp2 : this._p.toNat * this._p.toNat ≤ UInt64.size)
    (left right : SparsePolyZp)
    (initial : Generated.StrictHensel.HenselEEAState)
    (hexec : StrictHenselEEAStateInvariant this._p.toNat initial)
    (halgebra : StrictHenselEEAAlgebraicInvariant this._p.toNat left right
      initial) :
    ∀ state, HenselEEAPrefix (strictHenselEEARawOps this) initial state →
      StrictHenselEEAAlgebraicInvariant this._p.toNat left right state := by
  intro state hprefix
  induction hprefix with
  | refl => exact halgebra
  | step state quotient remainder hprefix hcontinue hrun ih =>
      have hstate := strictHenselEEAPrefix_stateInvariant this hcfg initial
        hexec state hprefix
      have hr1Ne : state.r1 ≠ #[] := by
        simpa [Array.isEmpty_iff] using hcontinue
      have hr1 : 0 < state.r1.size := by
        have hr1Size : state.r1.size ≠ 0 := by
          intro hzero
          apply hr1Ne
          exact Array.eq_empty_of_size_eq_zero hzero
        omega
      rcases henselDivmodVHCIR_refines this state.r0 state.r1 hcfg
          hstate.r0Canonical hstate.r1Canonical hr1 with
        ⟨output, houtputRun, houtputCorrect⟩
      have hrun' : henselDivmodVHCIR this state.r0 state.r1 =
          .ok (quotient, remainder) := by
        simpa [strictHenselEEARawOps] using hrun
      have houtput : output = (quotient, remainder) :=
        Except.ok.inj (houtputRun.symm.trans hrun')
      subst output
      exact StrictHenselEEAAlgebraicInvariant.step this._p.toNat left right
        state quotient remainder h2p hp2 hr1 ih houtputCorrect

private theorem henselDivmodVHCRefinesAux
    (this : DenseUPolyZp) (dividend divisor : SparsePolyZp)
    (hdivisor : 0 < divisor.size) (initialLimit : Nat)
    (initial : HenselDivmodVHCState)
    (hinvariant : HenselDivmodVHCExecutionInvariant this dividend divisor
      hdivisor initialLimit initial)
    (limit : Nat) (state : HenselDivmodVHCState)
    (hprefix : HenselDivmodVHCPrefix this dividend divisor hdivisor
      initialLimit initial limit state) :
    ∃ output,
      henselDivmodVHCOuterLoop this limit state dividend divisor hdivisor =
        .ok output ∧
      HenselDivmodVHCCorrect this dividend divisor hdivisor limit state
        output := by
  rw [henselDivmodVHCOuterLoop]
  split
  · rename_i hdone
    exact ⟨_, rfl, .done limit state hdone⟩
  · rename_i hnotDone
    rcases hinvariant.frontierReady limit state hprefix hnotDone with
      ⟨frontier, hselect, hdecrease⟩
    rw [hselect]
    simp only
    rw [dif_pos hdecrease]
    rcases hinvariant.iterationReady limit state frontier hprefix hselect
        hdecrease with ⟨next, hiteration⟩
    rw [hiteration]
    rcases henselDivmodVHCRefinesAux this dividend divisor hdivisor
        initialLimit initial hinvariant frontier.degree next
        (.step limit state next frontier hprefix hnotDone hselect hdecrease
          hiteration) with ⟨output, htailRun, htailCorrect⟩
    exact ⟨output, htailRun,
      .step limit state next frontier output hnotDone hselect hdecrease
        hiteration htailCorrect⟩
termination_by limit
decreasing_by assumption

theorem henselDivmodVHCOuterLoop_refines
    (this : DenseUPolyZp) (dividend divisor : SparsePolyZp)
    (hdivisor : 0 < divisor.size) (initialLimit : Nat)
    (initial : HenselDivmodVHCState)
    (hinvariant : HenselDivmodVHCExecutionInvariant this dividend divisor
      hdivisor initialLimit initial) :
    ∃ output,
      henselDivmodVHCOuterLoop this initialLimit initial dividend divisor
          hdivisor = .ok output ∧
      HenselDivmodVHCCorrect this dividend divisor hdivisor initialLimit
        initial output :=
  henselDivmodVHCRefinesAux this dividend divisor hdivisor initialLimit initial
    hinvariant initialLimit initial .refl

end StrictHensel

-- The discoverable public wrapper for the completed theorem above is emitted
-- by cpp2lean Pass 9 into `CLPoly.Refinement.Generated`.

end Refinement
